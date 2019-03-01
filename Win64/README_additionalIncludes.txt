Modified files of the OpenSees structure:
- TclModelBuilderSectionCommand.cpp (added OrthotropicMembraneSection and NoTensionSection3d - requires additional includes)
- ZerolengthND.cpp (383-448): reformulated stiffness matrix composition to allow for a non-symmetric stiffness matrix.
- TclElementCommands.cpp (added Macroelement3d - requires additional includes)
- TclModelBuilderNDMaterialCommand.cpp (added BeamFrictionSupport)
- TclModelBuilderUniaxialMaterialCommand.cpp (added CompressionDamage1d and TensionDamage1d)


Additional includes to add to the material library:
- DEVELOPER\element\cpp\Macroelement3d\Macroelement3d\NoTensionSection3d.h
- DEVELOPER\element\cpp\Macroelement3d\Macroelement3d\NoTensionSection3d.cpp
- DEVELOPER\material\cpp\OrthotropicMembraneSection\OrthotropicMembraneSection.h
- DEVELOPER\material\cpp\OrthotropicMembraneSection\OrthotropicMembraneSection.cpp
- DEVELOPER\material\cpp\CompressionDamage1d\CompressionDamage1d.h
- DEVELOPER\material\cpp\CompressionDamage1d\CompressionDamage1d.cpp
- DEVELOPER\material\cpp\TensionDamage1d\TensionDamage1d.h
- DEVELOPER\material\cpp\TensionDamage1d\TensionDamage1d.cpp
- DEVELOPER\material\cpp\BeamFrictionSupport\BeamFrictionSupport.h
- DEVELOPER\material\cpp\BeamFrictionSupport\BeamFrictionSupport.cpp


Additional includes to add to the element library:
- DEVELOPER\element\cpp\Macroelement3d\Macroelement3d\NoTensionSection3d.h
- DEVELOPER\element\cpp\Macroelement3d\Macroelement3d\NoTensionSection3d.cpp
- DEVELOPER\element\cpp\Macroelement3d\Macroelement3d\CohesiveSurface.h
- DEVELOPER\element\cpp\Macroelement3d\Macroelement3d\CohesiveSurface.cpp
- DEVELOPER\element\cpp\Macroelement3d\Macroelement3d\DamageShearInterface.h
- DEVELOPER\element\cpp\Macroelement3d\Macroelement3d\DamageShearInterface.cpp
- DEVELOPER\element\cpp\Macroelement3d\Macroelement3d\GambarottaLagomarsinoModel.h
- DEVELOPER\element\cpp\Macroelement3d\Macroelement3d\GambarottaLagomarsinoModel.cpp
- DEVELOPER\element\cpp\Macroelement3d\Macroelement3d\Macroelement3d.h
- DEVELOPER\element\cpp\Macroelement3d\Macroelement3d\Macroelement3d.cpp
- DEVELOPER\element\cpp\Macroelement3d\Macroelement3d\WrappedMaterial.h
- DEVELOPER\element\cpp\Macroelement3d\Macroelement3d\WrappedMaterial.cpp