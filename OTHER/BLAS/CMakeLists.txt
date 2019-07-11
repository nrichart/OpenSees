set(BLAS_VENDOR "All" CACHE STRING "Version of blas to use")
mark_as_advanced(BLAS_VENDOR)
set_property(CACHE BLAS_VENDOR PROPERTY STRINGS
  All
  ACML
  ACML_GPU
  ACML_MP
  ATLAS
  Apple
  CXML
  DXML
  Generic
  Goto
  IBMESSL
  Intel
  Intel10_32
  Intel10_64lp
  Intel10_64lp_seq
  NAS
  OpenBLAS
  PhiPACK
  SCSL
  SGIMATH
  SunPerf
  )

set(ENV{BLA_VENDOR} ${BLAS_VENDOR})

find_package(BLAS)

if(BLAS_FOUND)
  add_library(BLAS INTERFACE IMPORTED GLOBAL)
  message(STATUS "Defining target BLAS to link against ${BLAS_LIBRARIES}")
  set_property(TARGET BLAS PROPERTY INTERFACE_LINK_LIBRARIES ${BLAS_LIBRARIES})
else()
  add_library(BLAS idamax.o dasum.c daxpy.c dcopy.c ddot.c  dnrm2.c drot.c  drotg.c
    dscal.c  dswap.c dgemv.c dger.c  dsymv.c dsyr2.c dtrmm.c dgbmv.c
    dgemm.c dtrsm.c dtbsv.c dsyr.c dtrsv.c dsyrk.c dtrmv.c dsymm.c
    dsyr2k.c dpotrf.c dgeev.c dtrtrs.c)
endif()
