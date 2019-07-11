set(DISTRIBUTED_ALLAUX
  sp_ienv.c etree.c sp_colorder.c get_perm_c.c
  mmd.c comm.c memory.c util.c superlu_grid.c
  pxerr_dist.c superlu_timer.c GetDiagU.c symbfact.c
  psymbfact.c psymbfact_util.c mc64ad_dist.c
  static_schedule.c xerr_dist.c smach_dist.c dmach_dist.c)

set(DISTRIBUTED_DSLUSRC
  dlangs_dist.c dgsequ_dist.c dlaqgs_dist.c dutil_dist.c
  dmemory_dist.c dmyblas2_dist.c dsp_blas2_dist.c dsp_blas3_dist.c)

set(DISTRIBUTED_DPLUSRC
  pdgssvx.c pdgssvx_ABglobal.c
  dreadhb.c dreadrb.c dreadtriple.c dreadMM.c
  pdgsequ.c pdlaqgs.c dldperm_dist.c pdlangs.c pdutil.c
  pdsymbfact_distdata.c ddistribute.c pddistribute.c
  pdgstrf.c pdgstrf2.c
  pdgstrs.c pdgstrs1.c pdgstrs_lsum.c pdgstrs_Bglobal.c
  pdgsrfs.c pdgsmv.c pdgsrfs_ABXglobal.c pdgsmv_AXglobal.c)

foreach(_file ${DISTRIBUTED_ALLAUX} ${DISTRIBUTED_DSLUSRC} ${DISTRIBUTED_DPLUSRC})
  list(APPEND DISTRIBUTED_SUPERLU_SRC ${CMAKE_CURRENT_SOURCE_DIR}/SRC/${_file})
endforeach()

find_package(MPI REQUIRED)

add_library(DistributedSuperLU ${DISTRIBUTED_SUPERLU_SRC})
target_link_libraries(DistributedSuperLU PRIVATE BLAS PRIVATE ${MPI_C_LIBRARIES})
target_include_directories(DistributedSuperLU PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/SRC ${MPI_C_INCLUDE_DIRS})
target_compile_definitions(DistributedSuperLU PRIVATE -DUSE_VENDOR_BLAS)
