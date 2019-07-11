foreach(_third_party OPS_THIRD_PARITES)
  add_subdirectory(${_third_party})
  target_link_libraries(opensees PRIVATE ${${_third_party}_LIBRARY})
endforeach()
