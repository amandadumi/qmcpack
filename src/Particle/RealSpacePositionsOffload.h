//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


/** @file RealSpacePostionsOffload.h
 */
#ifndef QMCPLUSPLUS_REALSPACE_POSITIONS_OFFLOAD_H
#define QMCPLUSPLUS_REALSPACE_POSITIONS_OFFLOAD_H

#include "Particle/QuantumVariables.h"
#include "OhmmsSoA/Container.h"
#include "OpenMP/OMPallocator.hpp"
#include "Platforms/PinnedAllocator.h"

namespace qmcplusplus
{
/** Introduced to handle virtual moves and ratio computations, e.g. for non-local PP evaluations.
   */
class RealSpacePositionsOffload : public QuantumVariables
{
public:
  RealSpacePositionsOffload() : QuantumVariables(QuantumVariableKind::QV_POS_OFFLOAD) {}

  std::unique_ptr<QuantumVariables> makeClone() override { return std::make_unique<RealSpacePositionsOffload>(*this); }

  void resize(size_t n) override
  {
    RSoA.resize(n);
    RSoA_hostview.attachReference(RSoA.size(), RSoA.capacity(), RSoA.data());
  }

  size_t size() override { return RSoA_hostview.size(); }

  void setAllParticlePos(const ParticlePos_t& R) override
  {
    resize(R.size());
    RSoA_hostview.copyIn(R);
    RealType* data = RSoA.data();
    PRAGMA_OFFLOAD("omp target update to(data[0:RSoA.capacity()*QMCTraits::DIM])")
  }

  void setOneParticlePos(const PosType& pos, size_t iat) override
  {
    RSoA_hostview(iat) = pos;

    RealType x     = pos[0];
    RealType y     = pos[1];
    RealType z     = pos[2];
    RealType* data = RSoA.data();
    size_t offset  = RSoA.capacity();

    PRAGMA_OFFLOAD("omp target map(to : x, y, z, iat)")
    {
      data[iat]              = x;
      data[iat + offset]     = y;
      data[iat + offset * 2] = z;
    }
  }

  const PosVectorSoa& getAllParticlePos() override { return RSoA_hostview; }
  PosType getOneParticlePos(size_t iat) const override { return RSoA_hostview[iat]; }

private:
  ///particle positions in SoA layout
  VectorSoaContainer<RealType, QMCTraits::DIM, OMPallocator<RealType, PinnedAlignedAllocator<RealType>>> RSoA;

  ///host view of RSoA
  PosVectorSoa RSoA_hostview;
};
} // namespace qmcplusplus
#endif
