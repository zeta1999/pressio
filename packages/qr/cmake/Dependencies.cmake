TRIBITS_PACKAGE_DEFINE_DEPENDENCIES(
  LIB_REQUIRED_PACKAGES core
  LIB_OPTIONAL_PACKAGES
  #
  LIB_REQUIRED_TPLS  GTEST EIGEN BLAS LAPACK
  LIB_OPTIONAL_TPLS  MPI TRILINOS BLAZE ARMADILLO
  #REGRESSION_EMAIL_LIST simplecxx-regressions@someurl.none
  )
