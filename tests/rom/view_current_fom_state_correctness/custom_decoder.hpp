
#ifndef PRESSIO_TESTS_ROM_VIEW_CURRENT_FOM_STATE_CORRECTNESS_DECODER_HPP_
#define PRESSIO_TESTS_ROM_VIEW_CURRENT_FOM_STATE_CORRECTNESS_DECODER_HPP_

struct MyCustomDecoder
{
  // this is mandatory because pressio detects it
  using jacobian_type  = pressio::containers::MultiVector<Eigen::MatrixXd>;

private:
  const int romSize_ = {};
  mutable jacobian_type jac_;
  mutable int applyMappingCount_ = 0;

public:
  MyCustomDecoder(const int fomSize, const int romSize)
    : romSize_{romSize}, jac_(fomSize, romSize)
  {
    jac_.data()->setConstant(1);
  }

  const jacobian_type & getReferenceToJacobian() const{ return jac_; }

  template <typename rom_state_type>
  void updateJacobian(const rom_state_type &) const{}

  template <typename rom_state_type>
  void applyMapping(const rom_state_type & romState,
		    ::pressio::containers::Vector<Eigen::VectorXd> & result) const
  {
    ++applyMappingCount_;

    Eigen::MatrixXd A(result.extent(0), romSize_);
    A.setConstant(2);

    const auto & romStateNativeObj = *romState.data();
    auto & resultNativeObj = *result.data();
    resultNativeObj = A * romStateNativeObj;
  }
};

#endif