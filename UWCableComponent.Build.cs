using UnrealBuildTool;

/// <summary>
/// Build module rules.
/// </summary>
public class Wire : ModuleRules {

	public Wire(ReadOnlyTargetRules target) : base(target)
	{
		// Public module names that this module uses.
		PublicDependencyModuleNames.AddRange(new string[] { "Core", "CoreUObject", "Engine", "InputCore" });

		// The path for the header files
		PublicIncludePaths.AddRange(new string[] {"UWCableComponent"});

		// The path for the source files
		PrivateIncludePaths.AddRange(new string[] {"UWCableComponent"});
	}
}

