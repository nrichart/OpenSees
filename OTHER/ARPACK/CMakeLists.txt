set(ARPACK_SRC
  dgetv0.f dlaqrb.f dstqrb.f dsortc.f dsortr.f dstatn.f dstats.f
  dnaitr.f dnapps.f dnaup2.f dnaupd.f dnconv.f dneigh.f dngets.f
  dsaitr.f dsapps.f dsaup2.f dsaupd.f dsconv.f dseigt.f dsgets.f
  dneupd.f dseupd.f dsesrt.f dvout.f  dmout.f
  icnteq.f icopy.f iset.f iswap.f ivout.f second.f
  )

add_library(ARPACK ${ARPACK_SRC})
target_include_directories(ARPACK PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
