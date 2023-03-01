//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Amanda Dumi, amanda.e.dumi@gmail.com, Sandia National Laboratories
//////////////////////////////////////////////////////////////////////////////////////

#include "ForceRhijn.h"
#include "Particle/DistanceTable.h"
#include "Message/Communicate.h"
#include "Utilities/ProgressReportEngine.h"
#include "Numerics/DeterminantOperators.h"
#include "Numerics/MatrixOperators.h"
#include "OhmmsData/ParameterSet.h"
#include "OhmmsData/AttributeSet.h"



namespace qmcplusplus
{
ForceRhijn::ForceRhijn(ParticleSet& ions, ParticleSet& elns)
: ForceBase(ions, elns), d_aa_ID(ions.addTable(ions)), d_ei_ID(elns.addTable(ions))
{
  // Defaults
  nstep = 10;
  ///////////////////////////////////////////////////////////////
  
}
}