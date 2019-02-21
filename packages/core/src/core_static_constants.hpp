
#ifndef CORE_STATIC_CONSTANTS_HPP_
#define CORE_STATIC_CONSTANTS_HPP_


namespace rompp{ namespace core{

//template <typename T = void>
struct constants{

  // this is typically used as a template parameter
  // a positive quantity (e.g., a size of a vector)
  // is not known at compile-time, and the value
  // is defined at runtime
  static constexpr int dynamic = -1;

  template <typename scalar_t>
  static constexpr scalar_t one(){
    return static_cast<scalar_t>(1);
  }

  template <typename scalar_t>
  static constexpr scalar_t zero(){
    return static_cast<scalar_t>(0);
  }

};

}} // end of namespace rompp::core
#endif