add_library(AMD
  amd_aat.c amd_1.c amd_2.c amd_dump.c amd_postorder.c amd_post_tree.c amd_defaults.c
  amd_order.c amd_control.c amd_info.c amd_valid.c amd_preprocess.c
  SuiteSparse_config.c
  )

target_compile_options(AMD PRIVATE -O3 -fexceptions)
target_include_directories(AMD PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
