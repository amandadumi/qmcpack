//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_SNAP_JASTROW_BUILDER_H
#define QMCPLUSPLUS_SNAP_JASTROW_BUILDER_H
#include "QMCWaveFunctions/WaveFunctionComponentBuilder.h"

namespace qmcplusplus
{
class SNAPJastrowBuilder : public WaveFunctionComponentBuilder
{
public:
    SNAPJastrowBuilder(Communicate* comm, ParticleSet& target, ParticleSet& source);
    std::unique_ptr<WaveFunctionComponent> buildComponent(xmlNodePtr cur) override;


};
}
#endif