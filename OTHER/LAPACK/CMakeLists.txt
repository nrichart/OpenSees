find_package(LAPACK)

if(LAPACK_FOUND)
  add_library(LAPACK INTERFACE IMPORTED GLOBAL)
  message(STATUS "Defining target LAPACK to link against ${LAPACK_LIBRARIES}")
  set_property(TARGET LAPACK PROPERTY INTERFACE_LINK_LIBRARIES ${LAPACK_LIBRARIES})
else()
  set(LAPACK_DOBJ
    dgeqr2.c dlabad.c dlacon.c dlacpy.c dladiv.c dlae2.c dlaev2.c
    dlaexc.c dlagtm.c dlahqr.c dlaln2.c dlamch.c dlange.c dlanhs.c
    dlanst.c dlanv2.c dlaptm.c dlapy2.c dlapy3.c dlaran.c dlarf.c
    dlarfg.c dlarfx.c dlarnd.c dlarnv.c dlartg.c dlaruv.c dlascl.c
    dlaset.c dlasr.c  dlasrt.c dlassq.c dlasy2.c dorm2r.c dsteqr.c
    dtrevc.c dtrexc.c dtrsen.c dtrsyl.c dgbtf2.c dgbtrf.c dgbtrs.c
    dgttrf.c dgttrs.c dpttrf.c dpttrs.c dlaswp.c dgbsv.c dpbtrs.c
    dgesv.c dgetrs.c dpbsv.c dgetrf.c dgetf2.c dpbtrf.c dpotf2.c
    dpbtf2.c dgels.c dgeqrf.c dgelqf.c dormqr.c dormlq.c
    dgelq2.c dlarft.c dlarfb.c dorml2.c
    dlaebz.c dlagtf.c dlagts.c dlansb.c dlar2v.c dlargv.c
    dlartv.c dsbevx.c dsbtrd.c dstebz.c dstein.c dsterf.c
    dgetri.c dtrtri.c dtrti2.c dgebak.c dgebal.c dgeev.c
    dgehd2.c dgehrd.c dhseqr.c dlahrd.c
    dorg2r.c dorghr.c dorgqr.c dpotrf.c dpotrs.c dtrtrs.c
    dggbak.c dggbal.c dggev.c dgghrd.c dhgeqz.c dlag2.c dlasv2.c dtgevc.c
    dgesvxx.c dla_gerfsx_extended.c dlascl2.c dlaqge.c dgeequb.c dgerfsx.c
    chla_transtype.c dla_gercond.c dla_wwaddw.c ilaprec.c ilatrans.c dgecon.c
    dla_geamv.c dlacn2.c dla_lin_berr.c drscl.c dlatrs.c dla_rpvgrw.c dla_gerpvgrw.c
    BLAS_dgemv2_x.c BLAS_dgemv2_x-f2c.c
    BLAS_dgemv_x.c BLAS_dgemv_x-f2c.c BLAS_error.c BLAS_error.c
    )

  set(LAPACK_IOBJ
    ilaenv.c lsame.c lsamen.c xerbla.c xlaenv.c ieeeck.c)

  add_library(LAPACK ${LAPACK_DOBJ} ${LAPACK_IOBJ})
endif()
