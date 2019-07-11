add_opensees_extension(Truss2D Truss2D.cpp)

add_opensees_extension(Macroelement3d
  Macroelement3d/Macroelement3d/CohesiveSurface.cpp
  Macroelement3d/Macroelement3d/DamageShearInterface.cpp
  Macroelement3d/Macroelement3d/GambarottaLagomarsinoModel.cpp
  Macroelement3d/Macroelement3d/Macroelement3d.cpp
  Macroelement3d/Macroelement3d/NoTensionSection3d.cpp
  #Macroelement3d/Macroelement3d/stdafx.cpp
  Macroelement3d/Macroelement3d/WrappedMaterial.cpp
  )

target_include_directories(Macroelement3d PRIVATE Macroelement3d/Macroelement3d)
