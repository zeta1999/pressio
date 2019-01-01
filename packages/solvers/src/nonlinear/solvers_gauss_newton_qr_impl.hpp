
#ifndef SOLVERS_GAUSS_NEWTON_QR_IMPL_HPP
#define SOLVERS_GAUSS_NEWTON_QR_IMPL_HPP

#include "../solvers_ConfigDefs.hpp"
#include "../solvers_system_traits.hpp"
#include "../solvers_meta_static_checks.hpp"
#include "../../../CORE_OPS"
#include "../../../QR_BASIC"

namespace rompp{ namespace solvers{ namespace impl{

  template <typename system_t,
	    typename uint_t,
	    typename scalar_t,
	    typename qr_obj_t,
	    core::meta::enable_if_t<
	      ::rompp::solvers::details::system_traits<system_t>::is_system and
	      core::meta::is_core_vector_wrapper<typename system_t::state_type>::value and
	      core::meta::is_core_vector_wrapper<typename system_t::residual_type>::value
	      > * =nullptr
	    >
  void gauss_newtom_qr_solve(const system_t & sys,
			     typename system_t::state_type & x,
			     typename system_t::residual_type & resid,
			     typename system_t::jacobian_type & jacob,
			     uint_t maxNonLIt,
			     scalar_t tolerance,
			     typename system_t::state_type & QTResid,
			     typename system_t::state_type & dx,
			     qr_obj_t & qrObj,
			     scalar_t & normO,
			     scalar_t & normN)
  {

#ifdef DEBUG_PRINT
    std::cout << " starting Gauss-Newton solve "
	      << " tol = " << tolerance
	      << " maxIter = " << maxNonLIt
	      << std::endl;
#endif

    uint_t iStep = 1;
    while (iStep++ < maxNonLIt)
    {
      // QR decomposition of Jacobian
      qrObj.computeThin(jacob);

      // compute: Q^T Residual
      qrObj.project(resid, QTResid);

      // solve R dx = Q^T Residual
      qrObj.solve(QTResid, dx);

      // update solution
      x -= dx;

      normN = ::rompp::core::ops::norm2(dx);
#ifdef DEBUG_PRINT
      std::cout << " GN step=" << iStep
		<< " norm(dx)= " << normN
		<< std::endl;
#endif

      if (std::abs(normO - normN) < tolerance){
#ifdef DEBUG_PRINT
	std::cout << " GN converged! "
      		  << " final norm(dx)= " << normN
      		  << std::endl;
#endif
      	break;
      }

      normO = normN;
      sys.residual(x, resid);
      sys.jacobian(x, jacob);
    }

  }//


}}} //end namespace rompp::solvers::impl
#endif