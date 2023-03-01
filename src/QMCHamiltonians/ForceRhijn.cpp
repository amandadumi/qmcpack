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
    using WP = WalkerProperties::Indexes;

  // Defaults
  nstep = 10;
  ///////////////////////////////////////////////////////////////
# 

}
void ForceRhijn::calculate_gdd(){
      std::vector<RealType>::iterator Vit = values_.begin();

    int j       = 0;   // counts the steps for this walker to go back
    int FWindex = t_walker_->PHindex[p_ids_[i]] - 1;  // this is the current walkers index for a given property
    while (j < walker_lengths_[i][1])
    {
        int FWindex = t_walker_->PHindex[p_ids_[i]] - 1;
    }
void ForceRhijn::calculate_gdd(){

    }

void ForceRhijn::evaluate(ParticleSet& P){
    // for the current walker
    // find the id of the property of interest.
    // loop over nstep configurations
            // call to calc_gdd and calc_b 

    }
    
}