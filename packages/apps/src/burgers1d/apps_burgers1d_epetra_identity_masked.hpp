/*
//@HEADER
// ************************************************************************
//
// apps_burgers1d_epetra_identity_masked.hpp
//                     		  Pressio
//                             Copyright 2019
//    National Technology & Engineering Solutions of Sandia, LLC (NTESS)
//
// Under the terms of Contract DE-NA0003525 with NTESS, the
// U.S. Government retains certain rights in this software.
//
// Pressio is licensed under BSD-3-Clause terms of use:
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the copyright holder nor the names of its
// contributors may be used to endorse or promote products derived
// from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
// FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
// COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
// INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
// STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
// IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Francesco Rizzi (fnrizzi@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef APPS_BURGERS1D_APPS_BURGERS1D_EPETRA_IDENTITY_MASKED_HPP_
#define APPS_BURGERS1D_APPS_BURGERS1D_EPETRA_IDENTITY_MASKED_HPP_

#include "apps_burgers1d_epetra.hpp"

#include <Epetra_Import.h>

namespace pressio{ namespace apps{

class Burgers1dEpetraIdentityMask : public Burgers1dEpetra
{
public:
  using base_t = Burgers1dEpetra;
  using importer_t = Epetra_Import;
  using base_t::dense_matrix_type;
  using base_t::velocity_type;
  using base_t::state_type;
  using base_t::scalar_type;

  Burgers1dEpetraIdentityMask(std::vector<scalar_type> params,
			int Ncell, Epetra_MpiComm * comm)
    : base_t(params, Ncell, comm){
    this->setup();
  }

  ~Burgers1dEpetraIdentityMask() = default;

public:
  velocity_type createApplyMaskResult(const velocity_type & src) const
  {
    velocity_type dest(*maskMap_);
    dest.Import(src, *importer_, Insert);
    return dest;
  }

  dense_matrix_type createApplyMaskResult(const dense_matrix_type & src) const{
    dense_matrix_type dest(*maskMap_, src.NumVectors());
    dest.Import(src, *importer_, Insert);
    return dest;
  }

  template <typename T>
  void applyMask(const T & src, double time, T & dest) const{
    dest.Import(src, *importer_, Insert);
  }

private:
  void setup(){
    base_t::setup();
    // create a map to mimic the mask
    createMaskMap();
    importer_ = std::make_shared<importer_t>(*maskMap_, *dataMap_);
  };

  void createMaskMap(){
    // get # of my elements for the full map
    auto myN0 = dataMap_->NumMyElements();
    // get my global IDs
    auto myGID = dataMap_->MyGlobalElements();

    // pick all elements
    std::vector<int> myGIDnc;
    for (decltype(myN0) i=0; i<myN0; i++) {
	myGIDnc.emplace_back(myGID[i]);
    }
    maskMap_ = std::make_shared<Epetra_Map>(-1, myGIDnc.size(),
					    myGIDnc.data(), 0,
					    *comm_);
    //maskMap_->Print(std::cout);
  };

private:
  rcp<Epetra_Map> maskMap_;
  std::vector<int> myMaskGel_;
  rcp<importer_t> importer_;
};

}} //namespace pressio::apps
#endif  // APPS_BURGERS1D_APPS_BURGERS1D_EPETRA_IDENTITY_MASKED_HPP_
