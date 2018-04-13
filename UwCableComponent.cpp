// Copyright 1998-2017 Epic Games, Inc. All Rights Reserved. 

#include "UwCableComponent.h"
#include "EngineGlobals.h"
#include "PrimitiveViewRelevance.h"
#include "RenderResource.h"
#include "RenderingThread.h"
#include "WorldCollision.h"
#include "PrimitiveSceneProxy.h"
#include "VertexFactory.h"
#include "MaterialShared.h"
#include "SceneManagement.h"
#include "Engine/CollisionProfile.h"
#include "Materials/Material.h"
#include "LocalVertexFactory.h"
#include "Engine/Engine.h"
#include "CableComponentStats.h"
#include "DynamicMeshBuilder.h"

DECLARE_CYCLE_STAT(TEXT("Cable Sim"), STAT_uw_Cable_SimTime, STATGROUP_UwCableComponent);
DECLARE_CYCLE_STAT(TEXT("Cable Solve"), STAT_uw_Cable_SolveTime, STATGROUP_UwCableComponent);
DECLARE_CYCLE_STAT(TEXT("Cable Collision"), STAT_uw_Cable_CollisionTime, STATGROUP_UwCableComponent);
DECLARE_CYCLE_STAT(TEXT("Cable Integrate"), STAT_uw_Cable_IntegrateTime, STATGROUP_UwCableComponent);

static FName CableEndSocketName(TEXT("CableEnd"));
static FName CableStartSocketName(TEXT("CableStart"));

//////////////////////////////////////////////////////////////////////////

/** Vertex Buffer */
class FCableVertexBuffer : public FVertexBuffer 
{
public:
        virtual void InitRHI() override
        {
                FRHIResourceCreateInfo CreateInfo;
                VertexBufferRHI = RHICreateVertexBuffer(NumVerts * sizeof(FDynamicMeshVertex), BUF_Dynamic, CreateInfo);
        }

        int32 NumVerts;
};

/** Index Buffer */
class FCableIndexBuffer : public FIndexBuffer 
{
public:
        virtual void InitRHI() override
        {
                FRHIResourceCreateInfo CreateInfo;
                IndexBufferRHI = RHICreateIndexBuffer(sizeof(int32), NumIndices * sizeof(int32), BUF_Dynamic, CreateInfo);
        }

        int32 NumIndices;
};

/** Vertex Factory */
class FCableVertexFactory : public FLocalVertexFactory
{
public:

        FCableVertexFactory()
        {}


        /** Initialization */
        void Init(const FCableVertexBuffer* VertexBuffer)
        {
                if(IsInRenderingThread())
                {
                        // Initialize the vertex factory's stream components.
                        FDataType NewData;
                        NewData.PositionComponent = STRUCTMEMBER_VERTEXSTREAMCOMPONENT(VertexBuffer,FDynamicMeshVertex,Position,VET_Float3);
                        NewData.TextureCoordinates.Add(
                                FVertexStreamComponent(VertexBuffer,STRUCT_OFFSET(FDynamicMeshVertex,TextureCoordinate),sizeof(FDynamicMeshVertex),VET_Float2)
                                );
                        NewData.TangentBasisComponents[0] = STRUCTMEMBER_VERTEXSTREAMCOMPONENT(VertexBuffer,FDynamicMeshVertex,TangentX,VET_PackedNormal);
                        NewData.TangentBasisComponents[1] = STRUCTMEMBER_VERTEXSTREAMCOMPONENT(VertexBuffer,FDynamicMeshVertex,TangentZ,VET_PackedNormal);
                        SetData(NewData);
                }
                else
                {
                        ENQUEUE_UNIQUE_RENDER_COMMAND_TWOPARAMETER(
                                InitCableVertexFactory,
                                FCableVertexFactory*,VertexFactory,this,
                                const FCableVertexBuffer*,VertexBuffer,VertexBuffer,
                        {
                                // Initialize the vertex factory's stream components.
                                FDataType NewData;
                                NewData.PositionComponent = STRUCTMEMBER_VERTEXSTREAMCOMPONENT(VertexBuffer,FDynamicMeshVertex,Position,VET_Float3);
                                NewData.TextureCoordinates.Add(
                                        FVertexStreamComponent(VertexBuffer,STRUCT_OFFSET(FDynamicMeshVertex,TextureCoordinate),sizeof(FDynamicMeshVertex),VET_Float2)
                                        );
                                NewData.TangentBasisComponents[0] = STRUCTMEMBER_VERTEXSTREAMCOMPONENT(VertexBuffer,FDynamicMeshVertex,TangentX,VET_PackedNormal);
                                NewData.TangentBasisComponents[1] = STRUCTMEMBER_VERTEXSTREAMCOMPONENT(VertexBuffer,FDynamicMeshVertex,TangentZ,VET_PackedNormal);
                                VertexFactory->SetData(NewData);
                        });
                }
        }
};

/** Dynamic data sent to render thread */
struct FCableDynamicData
{
        /** Array of points */
        TArray<FVector> CablePoints;
};

//////////////////////////////////////////////////////////////////////////
// FCableSceneProxy

class FCableSceneProxy : public FPrimitiveSceneProxy
{
public:

        FCableSceneProxy(UUwCableComponent* Component)
                : FPrimitiveSceneProxy(Component)
                , Material(NULL)
                , DynamicData(NULL)
                , MaterialRelevance(Component->GetMaterialRelevance(GetScene().GetFeatureLevel()))
                , NumSegments(Component->NumSegments)
                , CableWidth(Component->CableWidth)
                , NumSides(Component->NumSides)
                , TileMaterial(Component->TileMaterial)
        {
                VertexBuffer.NumVerts = GetRequiredVertexCount();
                IndexBuffer.NumIndices = GetRequiredIndexCount();

                // Init vertex factory
                VertexFactory.Init(&VertexBuffer);

                // Enqueue initialization of render resource
                BeginInitResource(&VertexBuffer);
                BeginInitResource(&IndexBuffer);
                BeginInitResource(&VertexFactory);

                // Grab material
                Material = Component->GetMaterial(0);
                if(Material == NULL)
                {
                        Material = UMaterial::GetDefaultMaterial(MD_Surface);
                }
        }

        virtual ~FCableSceneProxy()
        {
                VertexBuffer.ReleaseResource();
                IndexBuffer.ReleaseResource();
                VertexFactory.ReleaseResource();

                if(DynamicData != NULL)
                {
                        delete DynamicData;
                }
        }

        int32 GetRequiredVertexCount() const
        {
                return (NumSegments + 1) * (NumSides + 1);
        }

        int32 GetRequiredIndexCount() const
        {
                return (NumSegments * NumSides * 2) * 3;
        }

        int32 GetVertIndex(int32 AlongIdx, int32 AroundIdx) const
        {
                return (AlongIdx * (NumSides+1)) + AroundIdx;
        }

        void BuildCableMesh(const TArray<FVector>& InPoints, TArray<FDynamicMeshVertex>& OutVertices, TArray<int32>& OutIndices)
        {
                const FColor VertexColor(255,255,255);
                const int32 NumPoints = InPoints.Num();
                const int32 SegmentCount = NumPoints-1;

                // Build vertices

                // We double up the first and last vert of the ring, because the UVs are different
                int32 NumRingVerts = NumSides+1;

                // For each point along spline..
                for(int32 PointIdx=0; PointIdx<NumPoints; PointIdx++)
                {
                        const float AlongFrac = (float)PointIdx/(float)SegmentCount; // Distance along cable

                        // Find direction of cable at this point, by averaging previous and next points
                        const int32 PrevIndex = FMath::Max(0, PointIdx-1);
                        const int32 NextIndex = FMath::Min(PointIdx+1, NumPoints-1);
                        const FVector ForwardDir = (InPoints[NextIndex] - InPoints[PrevIndex]).GetSafeNormal();

                        // Find quat from up (Z) vector to forward
                        const FQuat DeltaQuat = FQuat::FindBetween(FVector(0, 0, -1), ForwardDir);

                        // Apply quat orth vectors
                        const FVector RightDir = DeltaQuat.RotateVector(FVector(0, 1, 0));
                        const FVector UpDir = DeltaQuat.RotateVector(FVector(1, 0, 0));

                        // Generate a ring of verts
                        for(int32 VertIdx = 0; VertIdx<NumRingVerts; VertIdx++)
                        {
                                const float AroundFrac = float(VertIdx)/float(NumSides);
                                // Find angle around the ring
                                const float RadAngle = 2.f * PI * AroundFrac;
                                // Find direction from center of cable to this vertex
                                const FVector OutDir = (FMath::Cos(RadAngle) * UpDir) + (FMath::Sin(RadAngle) * RightDir);

                                FDynamicMeshVertex Vert;
                                Vert.Position = InPoints[PointIdx] + (OutDir * 0.5f * CableWidth);
                                Vert.TextureCoordinate = FVector2D(AlongFrac * TileMaterial, AroundFrac);
                                Vert.Color = VertexColor;
                                Vert.SetTangents(ForwardDir, OutDir ^ ForwardDir, OutDir);
                                OutVertices.Add(Vert);
                        }
                }

                // Build triangles
                for(int32 SegIdx=0; SegIdx<SegmentCount; SegIdx++)
                {
                        for(int32 SideIdx=0; SideIdx<NumSides; SideIdx++)
                        {
                                int32 TL = GetVertIndex(SegIdx, SideIdx);
                                int32 BL = GetVertIndex(SegIdx, SideIdx+1);
                                int32 TR = GetVertIndex(SegIdx+1, SideIdx);
                                int32 BR = GetVertIndex(SegIdx+1, SideIdx+1);

                                OutIndices.Add(TL);
                                OutIndices.Add(BL);
                                OutIndices.Add(TR);

                                OutIndices.Add(TR);
                                OutIndices.Add(BL);
                                OutIndices.Add(BR);
                        }
                }
        }

        /** Called on render thread to assign new dynamic data */
        void SetDynamicData_RenderThread(FCableDynamicData* NewDynamicData)
        {
                check(IsInRenderingThread());

                // Free existing data if present
                if(DynamicData)
                {
                        delete DynamicData;
                        DynamicData = NULL;
                }
                DynamicData = NewDynamicData;

                // Build mesh from cable points
                TArray<FDynamicMeshVertex> Vertices;
                TArray<int32> Indices;
                BuildCableMesh(NewDynamicData->CablePoints, Vertices, Indices);

                check(Vertices.Num() == GetRequiredVertexCount());
                check(Indices.Num() == GetRequiredIndexCount());

                void* VertexBufferData = RHILockVertexBuffer(VertexBuffer.VertexBufferRHI, 0, Vertices.Num() * sizeof(FDynamicMeshVertex), RLM_WriteOnly);
                FMemory::Memcpy(VertexBufferData, &Vertices[0], Vertices.Num() * sizeof(FDynamicMeshVertex));
                RHIUnlockVertexBuffer(VertexBuffer.VertexBufferRHI);

                void* IndexBufferData = RHILockIndexBuffer(IndexBuffer.IndexBufferRHI, 0, Indices.Num() * sizeof(int32), RLM_WriteOnly);
                FMemory::Memcpy(IndexBufferData, &Indices[0], Indices.Num() * sizeof(int32));
                RHIUnlockIndexBuffer(IndexBuffer.IndexBufferRHI);
        }

        virtual void GetDynamicMeshElements(const TArray<const FSceneView*>& Views, const FSceneViewFamily& ViewFamily, uint32 VisibilityMap, FMeshElementCollector& Collector) const override
        {
                QUICK_SCOPE_CYCLE_COUNTER( STAT_uw_CableSceneProxy_GetDynamicMeshElements );

                const bool bWireframe = AllowDebugViewmodes() && ViewFamily.EngineShowFlags.Wireframe;

                auto WireframeMaterialInstance = new FColoredMaterialRenderProxy(
                        GEngine->WireframeMaterial ? GEngine->WireframeMaterial->GetRenderProxy(IsSelected()) : NULL,
                        FLinearColor(0, 0.5f, 1.f)
                        );

                Collector.RegisterOneFrameMaterialProxy(WireframeMaterialInstance);

                FMaterialRenderProxy* MaterialProxy = NULL;
                if(bWireframe)
                {
                        MaterialProxy = WireframeMaterialInstance;
                }
                else
                {
                        MaterialProxy = Material->GetRenderProxy(IsSelected());
                }

                for (int32 ViewIndex = 0; ViewIndex < Views.Num(); ViewIndex++)
                {
                        if (VisibilityMap & (1 << ViewIndex))
                        {
                                const FSceneView* View = Views[ViewIndex];
                                // Draw the mesh.
                                FMeshBatch& Mesh = Collector.AllocateMesh();
                                FMeshBatchElement& BatchElement = Mesh.Elements[0];
                                BatchElement.IndexBuffer = &IndexBuffer;
                                Mesh.bWireframe = bWireframe;
                                Mesh.VertexFactory = &VertexFactory;
                                Mesh.MaterialRenderProxy = MaterialProxy;
                                BatchElement.PrimitiveUniformBuffer = CreatePrimitiveUniformBufferImmediate(GetLocalToWorld(), GetBounds(), GetLocalBounds(), true, UseEditorDepthTest());
                                BatchElement.FirstIndex = 0;
                                BatchElement.NumPrimitives = GetRequiredIndexCount()/3;
                                BatchElement.MinVertexIndex = 0;
                                BatchElement.MaxVertexIndex = GetRequiredVertexCount();
                                Mesh.ReverseCulling = IsLocalToWorldDeterminantNegative();
                                Mesh.Type = PT_TriangleList;
                                Mesh.DepthPriorityGroup = SDPG_World;
                                Mesh.bCanApplyViewModeOverrides = false;
                                Collector.AddMesh(ViewIndex, Mesh);

#if !(UE_BUILD_SHIPPING || UE_BUILD_TEST)
                                // Render bounds
                                RenderBounds(Collector.GetPDI(ViewIndex), ViewFamily.EngineShowFlags, GetBounds(), IsSelected());
#endif
                        }
                }
        }

        virtual FPrimitiveViewRelevance GetViewRelevance(const FSceneView* View) const override
        {
                FPrimitiveViewRelevance Result;
                Result.bDrawRelevance = IsShown(View);
                Result.bShadowRelevance = IsShadowCast(View);
                Result.bDynamicRelevance = true;
                MaterialRelevance.SetPrimitiveViewRelevance(Result);
                return Result;
        }

        virtual uint32 GetMemoryFootprint( void ) const override { return( sizeof( *this ) + GetAllocatedSize() ); }

        uint32 GetAllocatedSize( void ) const { return( FPrimitiveSceneProxy::GetAllocatedSize() ); }

private:

        UMaterialInterface* Material;
        FCableVertexBuffer VertexBuffer;
        FCableIndexBuffer IndexBuffer;
        FCableVertexFactory VertexFactory;

        FCableDynamicData* DynamicData;

        FMaterialRelevance MaterialRelevance;

        int32 NumSegments;

        float CableWidth;

        int32 NumSides;

        float TileMaterial;
};



//////////////////////////////////////////////////////////////////////////

UUwCableComponent::UUwCableComponent( const FObjectInitializer& ObjectInitializer )
        : Super( ObjectInitializer )
{
        UE_LOG(LogTemp, Warning, TEXT("Created UUwCableComponent"));
        PrimaryComponentTick.bCanEverTick = true;
        bTickInEditor = true;
        bAutoActivate = true;

        bAttachStart = true;
        bAttachEnd = true;
        CableWidth = 10.f;
        NumSegments = 10;
        NumSides = 4;
        EndLocation = FVector(100.f,0,0);
        CableLength = 100.f;
        SubstepTime = 0.02f;
        SolverIterations = 1;
        TileMaterial = 1.f;
        CollisionFriction = 0.2f;
        CableGravityScale = 1.f;

        SetCollisionProfileName(UCollisionProfile::PhysicsActor_ProfileName);
}

FPrimitiveSceneProxy* UUwCableComponent::CreateSceneProxy()
{
        UE_LOG(LogTemp, Warning, TEXT("Created UUwCableComponent CreateSceneProxy"));
        return new FCableSceneProxy(this);
}

int32 UUwCableComponent::GetNumMaterials() const
{
        return 1;
}


void UUwCableComponent::InitCablePos()
{
    const int32 NumParticles = NumSegments+1;

    Particles.Reset();
    Particles.AddUninitialized(NumParticles);

    FVector CableStart, CableEnd;
    GetEndPositions(CableStart, CableEnd);

    const FVector Delta = CableEnd - CableStart;
    const float SegmentLen=CableLength/NumSegments; 

    int32 FreeSegments = Delta.Size()/SegmentLen+1;
    int32 NonFreeSegments = NumSegments - FreeSegments;
    if (FreeSegments > NumSegments) { //cable is fully streched
        FreeSegments = NumSegments;
        NonFreeSegments = 0;
    }
    //const int32 FreeSegments=NumSegments-NumOfSpooledSegments;

    UE_LOG(LogTemp, Warning, TEXT("OnRegister UWCable %d %.3f %.3f"),FreeSegments,CableLength,Delta.Size());
    UE_LOG(LogTemp, Warning, TEXT("UWCable  end start (%.3f %.3f %.3f) (%.3f %.3f %.3f) "),
            CableEnd.X,CableEnd.Y,CableEnd.Z,
            CableStart.X,CableStart.Y,CableStart.Z);

    for(int32 ParticleIdx=0; ParticleIdx<NumParticles; ParticleIdx++)
    {
        FCableParticle& Particle = Particles[ParticleIdx];
        /*
           float Alpha = (float)(ParticleIdx-NumOfSpooledSegments)/(float)FreeSegments;
           if (Alpha<0 || FreeSegments<0) Alpha=0;
           FVector InitialPosition = CableStart + (Alpha * Delta);
           if (ParticleIdx < NumOfSpooledSegments) {
           InitialPosition=CableStart;
           }

           Particle.Position = InitialPosition;
           Particle.OldPosition = InitialPosition;
           Particle.bFree = Alpha>0; // default to free, will be fixed if desired in TickComponent
           */
        if (ParticleIdx <= NonFreeSegments) {
            Particle.Position = CableStart;
            Particle.OldPosition = CableStart;
            Particle.bFree = false;
        } else {
            float Alpha = (float)(ParticleIdx-NonFreeSegments)/(float)FreeSegments;
            FVector InitialPosition = CableStart + (Alpha * Delta);
            Particle.Position = InitialPosition;
            Particle.OldPosition = InitialPosition;
            Particle.bFree = true;
        }
    }
}


void UUwCableComponent::OnRegister()
{
        Super::OnRegister();
        InitCablePos();

}

void UUwCableComponent::UpdateLength()
{
        //search for the first non free particle
        //const int32 FreeSegments=NumSegments-NumOfSpooledSegments;
        const float SegmentLen=CableLength/NumSegments;//NumOfSpooledSegments;
        for(int32 ParticleIdx=0; ParticleIdx<NumSegments; ParticleIdx++)
        {
                FCableParticle& ParticleA = Particles[ParticleIdx];
                FCableParticle& ParticleB = Particles[ParticleIdx+1];
                if(!ParticleA.bFree && ParticleB.bFree) 
                {
                        const FVector P = ParticleB.Position;
                        //UE_LOG(LogTemp, Warning, TEXT("Cable first free: %.2f,%.2f,%.2f"),P.X,P.Y,P.Z);

                        FVector Delta = ParticleB.Position - ParticleA.Position;
                        float CurrentDistance = Delta.Size();
                        // if the segment is streached more then 10% release one segment from spool
                        if (CurrentDistance>SegmentLen*1.2) {
                                ParticleA.bFree=true;
                                //give a little push outside to avoid zero devition?
                                //
                                UE_LOG(LogTemp, Warning, TEXT("Cable Release: %d/%d %f"),ParticleIdx,NumSegments,
                                                CurrentDistance/SegmentLen);
                        }
                        break;
                }

        }
}

void UUwCableComponent::VerletIntegrate(float InSubstepTime, const FVector& Gravity)
{
        SCOPE_CYCLE_COUNTER(STAT_uw_Cable_IntegrateTime);

        const int32 NumParticles = NumSegments+1;
        const float SubstepTimeSqr = InSubstepTime * InSubstepTime;
        
        int32 CountFree=0;

        for(int32 ParticleIdx=0; ParticleIdx<NumParticles; ParticleIdx++)
        {
                FCableParticle& Particle = Particles[ParticleIdx];
                if(Particle.bFree)
                {
                        // Calc overall force
                        const FVector ParticleForce = Gravity + CableForce;

                        // Find vel
                        const FVector Vel = (Particle.Position - Particle.OldPosition)*(1.0-CableDampingCoef);
                        // Update position
                        FVector NewPosition = Particle.Position + Vel + (SubstepTimeSqr * ParticleForce);
                        if (CountFree>5) { // skip the first free few segments to avoid pulling the spool 
                                NewPosition.X+=FMath::RandRange(-1.0F,1.0F)*CableRandomForceCoef;
                                NewPosition.Y+=FMath::RandRange(-1.0F,1.0F)*CableRandomForceCoef;
                                NewPosition.Z+=FMath::RandRange(-1.0F,1.0F)*CableRandomForceCoef;
                        }
                        Particle.OldPosition = Particle.Position;
                        Particle.Position = NewPosition;
                        CountFree += 1;
                }
        }
}

/** Solve a single distance constraint between a pair of particles */
static void SolveDistanceConstraint(FCableParticle& ParticleA, FCableParticle& ParticleB, float DesiredDistance)
{
        // Find current vector between particles
        FVector Delta = ParticleB.Position - ParticleA.Position;
        // 
        float CurrentDistance = Delta.Size();
        if(!ParticleA.bFree && !ParticleB.bFree) return; //both locked cant move
        if(CurrentDistance<DesiredDistance/1000) {
                UE_LOG(LogTemp, Warning, TEXT("Cable small current distance %f"),CurrentDistance);
                return;
        }
        float ErrorFactor = (CurrentDistance - DesiredDistance)/CurrentDistance;

        // Only move free particles to satisfy constraints
        if(ParticleA.bFree && ParticleB.bFree)
        {
                ParticleA.Position += ErrorFactor * 0.5f * Delta;
                ParticleB.Position -= ErrorFactor * 0.5f * Delta;
        }
        else if(ParticleA.bFree)
        {
                ParticleA.Position += ErrorFactor * Delta;
        }
        else if(ParticleB.bFree)
        {
                ParticleB.Position -= ErrorFactor * Delta;
        }
}

void UUwCableComponent::SolveConstraints()
{
        SCOPE_CYCLE_COUNTER(STAT_uw_Cable_SolveTime);

        //const float SegmentLength = CableLength/(float)NumSegments;
        //const int32 FreeSegments=NumSegments-NumOfSpooledSegments;
        const float SegmentLength=CableLength/NumSegments;

        // For each iteration..
        for (int32 IterationIdx = 0; IterationIdx < SolverIterations; IterationIdx++)
        {
                // Solve distance constraint for each segment
                for (int32 SegIdx = 0; SegIdx < NumSegments; SegIdx++)
                {
                        FCableParticle& ParticleA = Particles[SegIdx];
                        FCableParticle& ParticleB = Particles[SegIdx + 1];
                        // Solve for this pair of particles
                        SolveDistanceConstraint(ParticleA, ParticleB, SegmentLength);
                }

                // If desired, solve stiffness constraints (distance constraints between every other particle)
                if (bEnableStiffness)
                {
                        for (int32 SegIdx = 0; SegIdx < NumSegments-1; SegIdx++)
                        {
                                FCableParticle& ParticleA = Particles[SegIdx];
                                FCableParticle& ParticleB = Particles[SegIdx + 2];
                                if(ParticleA.bFree && ParticleB.bFree)
                                        SolveDistanceConstraint(ParticleA, ParticleB, 2.f*SegmentLength);
                        }
                }
        }
}

void UUwCableComponent::PerformCableCollision()
{
        SCOPE_CYCLE_COUNTER(STAT_uw_Cable_CollisionTime);

        UWorld* World = GetWorld();
        // If we have a world, and collision is not disabled
        if (World && GetCollisionEnabled() != ECollisionEnabled::NoCollision)
        {
                // Get collision settings from component
                FCollisionQueryParams Params(SCENE_QUERY_STAT(CableCollision));

                ECollisionChannel TraceChannel = GetCollisionObjectType();
                FCollisionResponseParams ResponseParams(GetCollisionResponseToChannels());

                // Iterate over each particle
                for (int32 ParticleIdx = 0; ParticleIdx < Particles.Num(); ParticleIdx++)
                {
                        FCableParticle& Particle = Particles[ParticleIdx];
                        // If particle is free
                        if (Particle.bFree)
                        {
                                // Do sphere sweep
                                FHitResult Result;
                                bool bHit = World->SweepSingleByChannel(Result, Particle.OldPosition, Particle.Position, FQuat::Identity, TraceChannel, FCollisionShape::MakeSphere(0.5f * CableWidth), Params, ResponseParams);
                                // If we got a hit, resolve it
                                if (bHit)
                                {
                                        if (Result.bStartPenetrating)
                                        {
                                                Particle.Position += (Result.Normal * Result.PenetrationDepth);
                                        }
                                        else
                                        {
                                                Particle.Position = Result.Location;
                                        }

                                        // Find new velocity, after fixing collision
                                        FVector Delta = Particle.Position - Particle.OldPosition;
                                        // Find component in normal
                                        float NormalDelta = Delta | Result.Normal;
                                        // Find component in plane
                                        FVector PlaneDelta = Delta - (NormalDelta * Result.Normal);

                                        // Zero out any positive separation velocity, basically zero restitution
                                        Particle.OldPosition += (NormalDelta * Result.Normal);

                                        // Apply friction in plane of collision if desired
                                        if (CollisionFriction > KINDA_SMALL_NUMBER)
                                        {
                                                // Scale plane delta  by 'friction'
                                                FVector ScaledPlaneDelta = PlaneDelta * CollisionFriction;

                                                // Apply delta to old position reduce implied velocity in collision plane
                                                Particle.OldPosition += ScaledPlaneDelta;
                                        }
                                }
                        }
                }
        }
}

void UUwCableComponent::PerformSubstep(float InSubstepTime, const FVector& Gravity)
{
        SCOPE_CYCLE_COUNTER(STAT_uw_Cable_SimTime);

        VerletIntegrate(InSubstepTime, Gravity);

        SolveConstraints();

        UpdateLength();

        if (bEnableCollision)
        {
                PerformCableCollision();
        }
}

void UUwCableComponent::SetAttachEndTo(AActor* Actor, FName ComponentProperty, FName SocketName)
{
        UE_LOG(LogTemp, Warning, TEXT("UUwCableComponent SetAttachEndTo"));
        AttachEndTo.OtherActor = Actor;
        AttachEndTo.ComponentProperty = ComponentProperty;
        AttachEndToSocketName = SocketName;
}

AActor* UUwCableComponent::GetAttachedActor() const
{
        return AttachEndTo.OtherActor;
}

USceneComponent* UUwCableComponent::GetAttachedComponent() const
{
        return AttachEndTo.GetComponent(GetOwner());
}

void UUwCableComponent::GetCableParticleLocations(TArray<FVector>& Locations) const
{
        Locations.Empty();
        for (const FCableParticle& Particle : Particles)
        {
                Locations.Add(Particle.Position);
        }
}


void UUwCableComponent::GetEndPositions(FVector& OutStartPosition, FVector& OutEndPosition)
{
        // Start position is just component position
        OutStartPosition = GetComponentLocation();

        // See if we want to attach the other end to some other component
        USceneComponent* EndComponent = AttachEndTo.GetComponent(GetOwner());
        if(EndComponent == NULL)
        {
                EndComponent = this;
        }

        if (AttachEndToSocketName != NAME_None)
        {
                OutEndPosition = EndComponent->GetSocketTransform(AttachEndToSocketName).TransformPosition(EndLocation);
        }
        else
        {
                OutEndPosition = EndComponent->GetComponentTransform().TransformPosition(EndLocation);
        }

}

void UUwCableComponent::TickComponent(float DeltaTime, enum ELevelTick TickType, FActorComponentTickFunction *ThisTickFunction)
{
        Super::TickComponent(DeltaTime, TickType, ThisTickFunction);

        const FVector Gravity = FVector(0, 0, GetWorld()->GetGravityZ()) * CableGravityScale;

        // Update end points
        FVector CableStart, CableEnd;
        GetEndPositions(CableStart, CableEnd);

        FCableParticle& StartParticle = Particles[0];

        if (bAttachStart)
        {
                StartParticle.Position = StartParticle.OldPosition = CableStart;
                StartParticle.bFree = false;
        }
        else
        {
                StartParticle.bFree = true;
        }

        FCableParticle& EndParticle = Particles[NumSegments];
        if (bAttachEnd)
        {
                EndParticle.Position = EndParticle.OldPosition = CableEnd;
                EndParticle.bFree = false;
        }
        else
        {
                EndParticle.bFree = true;
        }

        // Ensure a non-zero substep
        float UseSubstep = FMath::Max(SubstepTime, 0.005f);

        // Perform simulation substeps
        TimeRemainder += DeltaTime;
        while(TimeRemainder > UseSubstep)
        {
                PerformSubstep(UseSubstep, Gravity);
                TimeRemainder -= UseSubstep;
        }

        // Need to send new data to render thread
        MarkRenderDynamicDataDirty();

        // Call this because bounds have changed
        UpdateComponentToWorld();
};

void UUwCableComponent::CreateRenderState_Concurrent()
{
        Super::CreateRenderState_Concurrent();

        SendRenderDynamicData_Concurrent();
}

void UUwCableComponent::SendRenderDynamicData_Concurrent()
{
        if(SceneProxy)
        {
                // Allocate cable dynamic data
                FCableDynamicData* DynamicData = new FCableDynamicData;

                // Transform current positions from particles into component-space array
                const FTransform& ComponentTransform = GetComponentTransform();
                int32 NumPoints = NumSegments+1;
                DynamicData->CablePoints.AddUninitialized(NumPoints);
                for(int32 PointIdx=0; PointIdx<NumPoints; PointIdx++)
                {
                        DynamicData->CablePoints[PointIdx] = ComponentTransform.InverseTransformPosition(Particles[PointIdx].Position);
                }

                // Enqueue command to send to render thread
                ENQUEUE_UNIQUE_RENDER_COMMAND_TWOPARAMETER(
                        FSendCableDynamicData,
                        FCableSceneProxy*,CableSceneProxy,(FCableSceneProxy*)SceneProxy,
                        FCableDynamicData*,DynamicData,DynamicData,
                {
                        CableSceneProxy->SetDynamicData_RenderThread(DynamicData);
                });
        }
}

FBoxSphereBounds UUwCableComponent::CalcBounds(const FTransform& LocalToWorld) const
{
        // Calculate bounding box of cable points
        FBox CableBox(ForceInit);
        for(int32 ParticleIdx=0; ParticleIdx<Particles.Num(); ParticleIdx++)
        {
                const FCableParticle& Particle = Particles[ParticleIdx];
                CableBox += Particle.Position;
        }

        // Expand by cable radius (half cable width)
        return FBoxSphereBounds(CableBox.ExpandBy(0.5f * CableWidth));
}

void UUwCableComponent::QuerySupportedSockets(TArray<FComponentSocketDescription>& OutSockets) const 
{
        OutSockets.Add(FComponentSocketDescription(CableEndSocketName, EComponentSocketType::Socket));
        OutSockets.Add(FComponentSocketDescription(CableStartSocketName, EComponentSocketType::Socket));
}

FTransform UUwCableComponent::GetSocketTransform(FName InSocketName, ERelativeTransformSpace TransformSpace) const
{
        int32 NumParticles = Particles.Num();
        if ((InSocketName == CableEndSocketName || InSocketName == CableStartSocketName) && NumParticles >= 2)
        {
                FVector ForwardDir, Pos;
                if (InSocketName == CableEndSocketName)
                {
                        FVector LastPos = Particles[NumParticles - 1].Position;
                        FVector PreviousPos = Particles[NumParticles - 2].Position;

                        ForwardDir = (LastPos - PreviousPos).GetSafeNormal();
                        Pos = LastPos;
                }
                else
                {
                        FVector FirstPos = Particles[0].Position;
                        FVector NextPos = Particles[1].Position;

                        ForwardDir = (NextPos - FirstPos).GetSafeNormal();
                        Pos = FirstPos;
                }

                const FQuat RotQuat = FQuat::FindBetween(FVector(1, 0, 0), ForwardDir);
                FTransform WorldSocketTM = FTransform(RotQuat, Pos, FVector(1, 1, 1));

                switch (TransformSpace)
                {
                        case RTS_World:
                        {
                                return WorldSocketTM;
                        }
                        case RTS_Actor:
                        {
                                if (const AActor* Actor = GetOwner())
                                {
                                        return WorldSocketTM.GetRelativeTransform(GetOwner()->GetTransform());
                                }
                                break;
                        }
                        case RTS_Component:
                        {
                                return WorldSocketTM.GetRelativeTransform(GetComponentTransform());
                        }
                }
        }

        return Super::GetSocketTransform(InSocketName, TransformSpace);
}

bool UUwCableComponent::HasAnySockets() const
{
        return (Particles.Num() >= 2);
}

bool UUwCableComponent::DoesSocketExist(FName InSocketName) const
{
        return (InSocketName == CableEndSocketName) || (InSocketName == CableStartSocketName);
}
