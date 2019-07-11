set(METIS_SRC
  coarsen.c fm.c initpart.c match.c ccgraph.c memory.c
  pmetis.c pqueue.c refine.c util.c timing.c debug.c
  bucketsort.c graph.c stat.c kmetis.c kwayrefine.c
  kwayfm.c balance.c ometis.c srefine.c sfm.c separator.c
  mincover.c mmd.c mesh.c meshpart.c frename.c fortran.c
  myqsort.c compress.c parmetis.c estmem.c
  mpmetis.c mcoarsen.c mmatch.c minitpart.c mbalance.c
  mrefine.c mutil.c mfm.c mkmetis.c mkwayrefine.c mkwayfmh.c
  mrefine2.c minitpart2.c mbalance2.c mfm2.c
  kvmetis.c kwayvolrefine.c kwayvolfm.c subdomains.c
  )

find_package(METIS)
if(METIS_FOUND)
  add_library(METIS INTERFACE IMPORTED GLOBAL)
  target_include_directories(METIS INTERFACE ${METIS_INCLUDE_DIR})
  target_link_libraries(METIS INTERFACE ${METIS_LIBRARY})
else()
  add_library(METIS ${METIS_SRC})
  target_include_directories(METIS PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
endif()
