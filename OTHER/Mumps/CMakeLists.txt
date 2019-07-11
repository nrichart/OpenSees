find_package(Mumps)

add_library(MUMPS INTERFACE IMPORTED GLOBAL)
target_link_libraries(MUMPS INTERFACE ${MUMPS_LIBRARIES})
target_include_directories(MUMPS INTERFACE ${MUMPS_INCLUDE_DIRS})
