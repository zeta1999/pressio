
#ifdef HAVE_TRILINOS
#ifndef SOLVERS_EXPERIMENTAL_POLICY_LINEAR_ITERATIVE_TRILINOS_HPP
#define SOLVERS_EXPERIMENTAL_POLICY_LINEAR_ITERATIVE_TRILINOS_HPP

#include "AztecOO.h"
#include <memory>

namespace rompp{
namespace solvers{

template <
  typename SolverT, 
  typename MatrixT
>
class SolversLinearIterativeTrilinosPolicy {

  public:

    static void resetLinearSystem(std::shared_ptr<SolverT>& solver, const MatrixT& A) {
      solver->SetUserMatrix(A.data(), true);
    }


    template <typename VectorT>
    static VectorT solve(std::shared_ptr<SolverT> solver, 
                        const VectorT& b, int maxIters, 
                        double tolerance) {
      auto x(b.data()->Map());
      solver->SetLHS(&x);
      solver->SetRHS(b.data());
      solver->Iterate(maxIters, tolerance);
      return VectorT(x);
    }

};

} // end namespace solvers

}//end namespace rompp
#endif
#endif
