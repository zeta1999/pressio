
#ifndef PRESSIO_ROM_EPETRA_SKELETON_HPP_
#define PRESSIO_ROM_EPETRA_SKELETON_HPP_

#include "Epetra_MpiComm.h"
#include "Epetra_Map.h"
#include "Epetra_IntVector.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Time.h"

namespace pressio{ namespace rom{ namespace test{

class EpetraSkeleton{
protected:
  using nativeVec	= Epetra_Vector;

/* these types exposed because need to be detected */
public:
  using scalar_type	= double;
  using state_type	= Epetra_Vector;
  using velocity_type	= state_type;
  using jacobian_type	= Epetra_CrsMatrix;
  using dense_matrix_type = Epetra_MultiVector;

public:
  EpetraSkeleton() = default;
  ~EpetraSkeleton() = default;

public:
  velocity_type createVelocity() const;
  dense_matrix_type createApplyJacobianResult(const dense_matrix_type & B) const;

  void velocity(const state_type & u,
		const scalar_type /* t */,
		velocity_type & f) const;

  // computes: A = Jac B where B is a multivector
  void applyJacobian(const state_type & y,
		     const dense_matrix_type & B,
		     scalar_type t,
		     dense_matrix_type & A) const;

};//end class

}}} //namespace pressio::rom::test
#endif
