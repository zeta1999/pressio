namespace pressio{ namespace rom{ namespace wls{ namespace impl{

//struct JTJ_JTR_policy_standard{
//public:
//  void operator()(wls_state_type wls_state, fom_t appObj) const { 
//    std::cout << "Computed JTJ and JTR through default policy" << std::endl;
//}
//};

template< typename wls_state_type, typename fom_type, typename jtj_type, typename jtr_type>
struct JTJ_JTR_policy_smart{
public:
  using jtj_t = jtj_type;
  using jtr_t = jtr_type;
  void operator()(wls_state_type  & wls_state, jtj_type & jtj, jtr_type & jtr) const { 
  std::cout << "Computed JTJ and JTR through the smart policy" << std::endl; 
}
};



template< typename wls_state_type, typename fom_type, typename residual_type>
struct residual_policy_naive{
public:
  using residual_t = residual_type;
  void operator()(wls_state_type  & wls_state, residual_t & residual) const { 
  std::cout << "Computed the (void) residual through the naive policy" << std::endl;
} 
  residual_t evaluate(wls_state_type  & wls_state) const { 
  std::cout << "Computed the residual through the naive policy" << std::endl; 
  residual_t a;// placeholder to return residual
  return a;
}
};


template< typename wls_state_type, typename fom_type, typename jacobian_type>
struct jacobian_policy_naive{
public:
  using jacobian_t = jacobian_type;
  void operator()(wls_state_type  & wls_state, jacobian_t & jacobian) const { 
  std::cout << "Computed the (void) jacobian through the naive policy" << std::endl; 
}
};




}}}}

