// Fill out your copyright notice in the Description page of Project Settings.


#include "NewUWCableActor.h"


// Sets default values
ANewUWCableActor::ANewUWCableActor()
{
 	// Set this actor to call Tick() every frame.  You can turn this off to improve performance if you don't need it.
	PrimaryActorTick.bCanEverTick = true;
	UwCableComponent = CreateDefaultSubobject<UUwCableComponent>(TEXT("UwCableComponent0"));
	//CableComponent = NewObject<UUwCableComponent>(this);
	RootComponent = UwCableComponent;

}

// Called when the game starts or when spawned
void ANewUWCableActor::BeginPlay()
{
	Super::BeginPlay();
	UE_LOG(LogTemp, Warning, TEXT("Created AUWCableActor")); 
	
}

// Called every frame
void ANewUWCableActor::Tick(float DeltaTime)
{
	Super::Tick(DeltaTime);

}

