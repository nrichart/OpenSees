add_opensees_extension(Truss2D Truss2D.cpp)

add_opensees_extension(NoTensionSection3d
  Macroelement3d/Macroelement3d/NoTensionSection3d.cpp
  )
add_opensees_extension(CohesiveSurface
  Macroelement3d/Macroelement3d/CohesiveSurface.cpp
  )
add_opensees_extension(Macroelement3d
  Macroelement3d/Macroelement3d/Macroelement3d.cpp
  )
add_opensees_extension(DamageShearInterface
  Macroelement3d/Macroelement3d/DamageShearInterface.cpp
  )
