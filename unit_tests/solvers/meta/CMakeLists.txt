

set(SRC1
  ${CMAKE_CURRENT_SOURCE_DIR}/solvers_res_jac_api.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/solvers_hess_grad_api.cc
  )
add_serial_utest(solvers_meta_serial "${SRC1}")
