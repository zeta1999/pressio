

add_subdirectory(linear)
add_subdirectory(non_linear_least_squares)
add_serial_utest(solvers_newton_raphson_eigen solvers_newton_raphson_eigen.cc)
