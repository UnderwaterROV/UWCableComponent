// Fill out your copyright notice in the Description page of Project Settings.

#pragma once

#include "CoreMinimal.h"
#include "GameFramework/Actor.h"
#include "UwCableComponent.h"
#include "NewUWCableActor.generated.h"

UCLASS()
class OCEANTEST1_API ANewUWCableActor : public AActor
{
	GENERATED_BODY()
	
public:	
	// Sets default values for this actor's properties
	ANewUWCableActor();
	
	UPROPERTY(Category=Cable, VisibleAnywhere, BlueprintReadOnly)
	class UUwCableComponent* UwCableComponent;

protected:
	// Called when the game starts or when spawned
	virtual void BeginPlay() override;

public:	
	// Called every frame
	virtual void Tick(float DeltaTime) override;

	
	
};
