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
#include "SNAPJastrow.h"
namespace qmcplusplus
{
// forward declaration
class ParticleSet;
class SNAPJastrowBuilder : public WaveFunctionComponentBuilder
{
public:
    const ParticleSet& sourcePtcl;
    SNAPJastrowBuilder(Communicate* comm, ParticleSet& target, const ParticleSet& source)
    : WaveFunctionComponentBuilder(comm, target), sourcePtcl(source)
{
    ClassName = "SnapJastrowBuilder";
    NameOpt = "0";
    TypeOpt = "SNAP";
    SNAPType = "linear";
}
    std::unique_ptr<WaveFunctionComponent> buildComponent(xmlNodePtr cur) override;
    std::unique_ptr<WaveFunctionComponent> createSNAP(xmlNodePtr cur);
    bool putkids(xmlNodePtr kids, SNAPJastrow& SJ);

private:
  std::string NameOpt;
  std::string TypeOpt;
  std::string SNAPType;
  std::string RegionOpt;
  std::string SourceOpt;

};
}
#endif