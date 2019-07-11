set(CMAKE_SHARED_MODULE_PREFIX "")
macro(add_opensees_extension name)
  message(STATUS "Adding extension ${name}")
  add_library(${name} MODULE ${ARGN})
  target_include_directories(${name}
    PRIVATE ${PROJECT_SOURCE_DIR}/DEVELOPER/core)
  target_link_options(${name}
    PRIVATE
    -undefined,suppress
    )
  if(APPLE)
    target_link_options(${name}
      PRIVATE
      -flat_namespace
      "-current_version 1.0"
      "-compatibility_version 1.0"
      )
  endif()
  target_compile_options(${name}
    PRIVATE
    ${OS_FLAG}
    -fvisibility=hidden
    )
  foreach(_file ${ARGN})
    if("${_file}" MATCHES ".c$")
      set_source_files_properties(${_file} PROPERTIES LANGUAGE CXX)
    endif()
  endforeach()
endmacro()

foreach(_subfolder element material recorder system integrator)
  foreach(_lang c cpp fortran)
    if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${_subfolder}/${_lang}/CMakeLists.txt)
      add_subdirectory(${_subfolder}/${_lang})
    endif()
  endforeach()

  if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${_subfolder}/CMakeLists.txt)
    add_subdirectory(${_subfolder})
  endif()
endforeach()

