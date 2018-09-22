
#ifndef SOLVERS_EXPERIMENTAL_LINEAR_FACTORY_HPP
#define SOLVERS_EXPERIMENTAL_LINEAR_FACTORY_HPP

#include <iostream>
#include <memory>
#include <type_traits>
#include <Eigen/Sparse>

#include "solvers_linear_direct.hpp"
#include "solvers_linear_iterative.hpp"
#include "solvers_linear_traits.hpp"
#include "solvers_policy_linear_dense_eigen.hpp"
#include "solvers_policy_linear_iterative_eigen.hpp"
#include "solvers_policy_linear_iterative_trilinos.hpp"

#include "matrix/core_matrix_traits.hpp"


namespace solvers {


struct LinearSolvers {

  /**
   * Create a direct linear solver for dense Eigen matrices
   *
   * @return wrapped_solver direct linear solver for Eigen matrices
   */
  template <
    typename SolverT,
    typename MatrixT,
    typename std::enable_if<
      core::details::traits<MatrixT>::wrapped_package_identifier == core::details::WrappedPackageIdentifier::Eigen,
      MatrixT
    >::type* = nullptr
  >
  static auto createDirectSolver() {
    using solver_traits = linear::details::solver_traits<SolverT>;

    static_assert(solver_traits::eigen_enabled && solver_traits::direct, "Solver not available for linear systems defined by Eigen matrices");

    using wrapped_type = typename core::details::traits<MatrixT>::wrapped_t;
    using concrete_solver_type = typename solver_traits::template eigen_solver_type<wrapped_type>;
    using concrete_policy_type = SolversLinearDenseEigenPolicy<concrete_solver_type, MatrixT>;

    auto solver = std::make_shared<concrete_solver_type>();
    auto wrapped_solver = LinearDirectSolver<concrete_solver_type, MatrixT, concrete_policy_type>(solver);

    return wrapped_solver;
  }


  /**
   * Create a direct linear solver for Eigen matrices
   *
   * @param A matrix representing the linear system to be solved
   * @return solver direct linear solver for Eigen matrices
   */
  template <
    typename SolverT,
    typename MatrixT
  >
  static auto createDirectSolver(const MatrixT& A) {
    auto solver = LinearSolvers::createDirectSolver<SolverT, MatrixT>();
    solver.resetLinearSystem(A);
    return solver;
  }


  /**
   * Create a linear iterative solver for sparse Eigen matrices
   *
   * @return wrapped_solver sparse linear iterative solver for Eigen matrices
   */
  template <
    typename SolverT,
    typename MatrixT,
    typename PrecT = linear::DefaultPreconditioner,
    typename std::enable_if<
      core::details::traits<MatrixT>::wrapped_package_identifier == core::details::WrappedPackageIdentifier::Eigen,
      MatrixT
    >::type* = nullptr
  >
  static auto createIterativeSolver() {

    using solver_traits = linear::details::solver_traits<SolverT>;
    using preconditioner_traits = linear::details::preconditioner_traits<PrecT>;

    static_assert(solver_traits::eigen_enabled && !solver_traits::direct, "Solver not available for linear systems defined by Eigen matrices");
    static_assert(preconditioner_traits::eigen_enabled, "Preconditioner not available for linear systems defined by Eigen matrices");

    using wrapped_type = typename core::details::traits<MatrixT>::wrapped_t;
    using concrete_precon_type = typename preconditioner_traits::template eigen_preconditioner_type<wrapped_type>;
    using concrete_solver_type = typename solver_traits::template eigen_solver_type<wrapped_type, concrete_precon_type>;
    using concrete_policy_type = SolversLinearIterativeEigenPolicy<concrete_solver_type, MatrixT>;

    auto solver = std::make_shared<concrete_solver_type>();
    auto wrapped_solver = LinearIterativeSolver<concrete_solver_type, MatrixT, concrete_policy_type>(solver);

    return wrapped_solver;
  }


  /**
   * Create a linear iterative solver for sparse Eigen matrices
   *
   * @param  A sparse matrix representing the linear system to be solved
   * @return solver sparse linear iterative solver for Eigen matrices
    */
  template <
    typename SolverT,
    typename MatrixT,
    typename PrecT = linear::DefaultPreconditioner
  >
  static auto createIterativeSolver(const MatrixT& A) {
    auto solver = LinearSolvers::createIterativeSolver<SolverT, MatrixT, PrecT>();
    solver.resetLinearSystem(A);
    return solver;
  }


  /**
   * Create a linear iterative solver for sparse Trilinos matrices
   *
   * @return An iterative linear solver for sparse Trilinos matrices
   */
  template <
    typename SolverT,
    typename MatrixT,
    typename PrecT = linear::DefaultPreconditioner,
    typename std::enable_if<
      core::details::traits<MatrixT>::wrapped_package_identifier == core::details::WrappedPackageIdentifier::Trilinos,
      MatrixT
    >::type* = nullptr
  >
  static auto createIterativeSolver() {

    using solver_traits = linear::details::solver_traits<SolverT>;
    using preconditioner_traits = linear::details::preconditioner_traits<PrecT>;

    static_assert(solver_traits::trilinos_enabled && !solver_traits::direct, "Solver not available for linear systems defined by Trilinos matrices");
    static_assert(preconditioner_traits::trilinos_enabled, "Preconditioner not available for linear systems defined by Trilinos matrices");

    using concrete_solver_type = AztecOO;
    using concrete_policy_type = SolversLinearIterativeTrilinosPolicy<concrete_solver_type, MatrixT>;

    auto solver_flag = solver_traits::trilinos_flag;
    auto precon_flag = preconditioner_traits::trilinos_flag;

    auto solver = std::make_shared<concrete_solver_type>();
    solver->SetAztecOption(AZ_solver, solver_flag);

    switch (precon_flag) {
      case AZ_ilu: {
        solver->SetAztecOption(AZ_precond, AZ_dom_decomp);
        solver->SetAztecOption(AZ_subdomain_solve, AZ_ilu);
        break;
      }
      case AZ_ilut: {
        solver->SetAztecOption(AZ_precond, AZ_dom_decomp);
        solver->SetAztecOption(AZ_subdomain_solve, AZ_ilut);
        solver->SetAztecOption(AZ_overlap, 1);
        solver->SetAztecOption(AZ_ilut_fill, 3.0);
        break;
      }
      case AZ_icc: {
        solver->SetAztecOption(AZ_precond, AZ_dom_decomp);
        solver->SetAztecOption(AZ_subdomain_solve, AZ_icc);
        break;
      }
    }

    auto wrapped_solver = LinearIterativeSolver<concrete_solver_type, MatrixT, concrete_policy_type>(solver);

    return wrapped_solver;
  }


  /**
   * Create an iterative least square linear solver for sparse Eigen matrices
   *
   * @return iterative least square linear solver for sparse Eigen matrices
   */
  template <
    typename MatrixT,
    typename PrecT = linear::DefaultPreconditioner,
    typename std::enable_if<
      core::details::traits<MatrixT>::wrapped_package_identifier == core::details::WrappedPackageIdentifier::Eigen,
      MatrixT
    >::type* = nullptr
  >
  static auto createLeastSquareIterativeSolver() {
    return LinearSolvers::createIterativeSolver<typename linear::LSCG, MatrixT, PrecT>();
  }


  /**
   * Create an iterative least square linear solver for sparse Eigen matrices
   *
   * @param A sparse Eigen matrix defining the least square system to be solved
   * @return iterative least square linear solver for sparse Eigen matrices
   */
  template <
    typename MatrixT,
    typename PrecT = linear::DefaultPreconditioner,
    typename std::enable_if<
      core::details::traits<MatrixT>::wrapped_package_identifier == core::details::WrappedPackageIdentifier::Eigen,
      MatrixT
    >::type* = nullptr
  >
  static auto createLeastSquareIterativeSolver(const MatrixT& A) {
    return LinearSolvers::createIterativeSolver<typename linear::LSCG, MatrixT, PrecT>(A);
  }


};

} // end namespace solvers

#endif