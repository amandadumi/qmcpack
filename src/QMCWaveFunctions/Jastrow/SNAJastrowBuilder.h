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


#ifndef QMCPLUSPLUS_SNA_JASTROW_BUILDER_H
#define QMCPLUSPLUS_SNA_JASTROW_BUILDER_H
#include "QMCWaveFunctions/WaveFunctionComponentBuilder.h"
#include "SNAJastrow.h"
namespace qmcplusplus
{
// forward declaration
class ParticleSet;
class SNAJastrowBuilder : public WaveFunctionComponentBuilder
{
public:
    ParticleSet& jtarget;
    ParticleSet& jsource;
    SNAJastrowBuilder(Communicate* comm, ParticleSet& target, ParticleSet& source)
    : WaveFunctionComponentBuilder(comm, target), jtarget(target), jsource(source)
{
    ClassName = "SnapJastrowBuilder";
    NameOpt = "0";
    TypeOpt = "SNA";
    SNAType = "linear";
    
}
    std::unique_ptr<WaveFunctionComponent> buildComponent(xmlNodePtr cur) override;
    std::unique_ptr<WaveFunctionComponent> createSNA(xmlNodePtr cur);
    bool putkids(xmlNodePtr kids, SNAJastrow& SJ);


private:
  std::string NameOpt;
  std::string TypeOpt;
  std::string SNAType;
  std::string RegionOpt;
  std::string SourceOpt;
};
}
#endif
