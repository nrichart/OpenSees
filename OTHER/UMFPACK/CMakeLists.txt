set(UMFCH
  umf_assemble.c umf_blas3_update.c umf_build_tuples.c umf_create_element.c 
  umf_extend_front.c umf_garbage_collection.c umf_get_memory.c 
  umf_init_front.c umf_kernel.c umf_kernel_init.c umf_kernel_wrapup.c 
  umf_local_search.c umf_lsolve.c umf_ltsolve.c umf_mem_alloc_element.c 
  umf_mem_alloc_head_block.c umf_mem_alloc_tail_block.c 
  umf_mem_free_tail_block.c umf_mem_init_memoryspace.c 
  umf_report_vector.c umf_row_search.c umf_scale_column.c 
  umf_set_stats.c umf_solve.c umf_symbolic_usage.c umf_transpose.c 
  umf_tuple_lengths.c umf_usolve.c umf_utsolve.c umf_valid_numeric.c 
  umf_valid_symbolic.c umf_grow_front.c umf_start_front.c 
  umf_store_lu.c umf_scale.c
  )

# non-user-callable umf_*.[ch] files, int/SuiteSparse_long versions only
# (no real/complex):
set(UMFINT
  umf_analyze.c umf_apply_order.c umf_colamd.c umf_free.c umf_fsize.c
  umf_is_permutation.c umf_malloc.c umf_realloc.c umf_report_perm.c
  umf_singletons.c umf_cholmod.c
  )

# user-callable umfpack_*.[ch] files (int/SuiteSparse_long and real/complex):
set(UMFPACK
  umfpack_col_to_triplet.c umfpack_defaults.c umfpack_free_numeric.c
  umfpack_free_symbolic.c umfpack_get_numeric.c umfpack_get_lunz.c
  umfpack_get_symbolic.c umfpack_get_determinant.c umfpack_numeric.c
  umfpack_qsymbolic.c umfpack_report_control.c umfpack_report_info.c
  umfpack_report_matrix.c umfpack_report_numeric.c umfpack_report_perm.c
  umfpack_report_status.c umfpack_report_symbolic.c umfpack_report_triplet.c
  umfpack_report_vector.c umfpack_solve.c umfpack_symbolic.c
  umfpack_transpose.c umfpack_triplet_to_col.c umfpack_scale.c
  umfpack_load_numeric.c umfpack_save_numeric.c
  umfpack_load_symbolic.c umfpack_save_symbolic.c
  )

# user-callable, only one version for int/SuiteSparse_long,
# real/complex, *.[ch] files:
set(GENERIC umfpack_timer.c umfpack_tictoc.c)

set(UMFPACK_SRCs
  ${GENERIC}
  ${UMFINT}
  ${UMF}
  ${USER}
  )

include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/../AMD
  )

add_library(umf_i OBJECT ${UMFINT})
target_compile_definitions(umf_i PRIVATE -DDINT -DNCHOLMOD)

add_library(umf_di OBJECT ${UMFCH})
target_compile_definitions(umf_di PRIVATE -DDINT)

add_library(umf_di_hsolve OBJECT umf_ltsolve.c umf_utsolve.c)
target_compile_definitions(umf_di_hsolve PRIVATE -DDINT -DCONJUGATE_SOLVE)

add_library(umf_di_triplet_map_x OBJECT umf_triplet.c)
target_compile_definitions(umf_di_triplet_map_x PRIVATE -DDINT -DDO_MAP -DDO_VALUES)

add_library(umf_di_triplet_map_nox OBJECT umf_triplet.c)
target_compile_definitions(umf_di_triplet_map_nox PRIVATE -DDINT -DDO_MAP)

add_library(umf_di_triplet_nomap_x OBJECT umf_triplet.c)
target_compile_definitions(umf_di_triplet_nomap_x PRIVATE -DDINT -DDO_VALUES)

add_library(umf_di_triplet_nomap_nox OBJECT umf_triplet.c)
target_compile_definitions(umf_di_triplet_nomap_nox PRIVATE -DDINT)

add_library(umf_di_assemble_fixq OBJECT umf_assemble.c)
target_compile_definitions(umf_di_assemble_fixq PRIVATE -DDINT -DFIXQ)

add_library(umf_di_store_lu_drop OBJECT umf_store_lu.c)
target_compile_definitions(umf_di_store_lu_drop PRIVATE -DDINT -DDROP)

add_library(umfpack_di_wsolve OBJECT umfpack_solve.c)
target_compile_definitions(umfpack_di_wsolve PRIVATE -DDINT -DWSOLVE)

add_library(umfpack_di OBJECT ${UMFPACK})
target_compile_definitions(umfpack_di PRIVATE -DDINT)

add_library(umfpack_gn OBJECT ${GENERIC})

foreach(_target umf_i umf_di umf_di_hsolve umf_di_triplet_map_x umf_di_triplet_map_nox
    umf_di_triplet_nomap_x umf_di_triplet_nomap_nox umf_di_assemble_fixq umf_di_store_lu_drop
    umfpack_di_wsolve umfpack_di umfpack_gn)
  set_property(TARGET ${_target} PROPERTY POSITION_INDEPENDENT_CODE TRUE)
endforeach()

add_library(UMFPACK
  $<TARGET_OBJECTS:umf_i>
  $<TARGET_OBJECTS:umf_di>
  $<TARGET_OBJECTS:umf_di_hsolve>
  $<TARGET_OBJECTS:umf_di_triplet_map_x>
  $<TARGET_OBJECTS:umf_di_triplet_map_nox>
  $<TARGET_OBJECTS:umf_di_triplet_nomap_x>
  $<TARGET_OBJECTS:umf_di_triplet_nomap_nox>
  $<TARGET_OBJECTS:umf_di_assemble_fixq>
  $<TARGET_OBJECTS:umf_di_store_lu_drop>
  $<TARGET_OBJECTS:umfpack_di_wsolve>
  $<TARGET_OBJECTS:umfpack_di>
  $<TARGET_OBJECTS:umfpack_gn> 
  )
target_include_directories(UMFPACK PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/../AMD)
