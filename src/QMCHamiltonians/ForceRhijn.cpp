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
        int FWi
        ndex = t_walker_->PHindex[p_ids_[i]] - 1;
    }
void ForceRhijn::calculate_gdd(){

    }

void ForceRhijn::evaluate(ParticleSet& P){
    int i = 0 //todo figure out which observable needs to be here for history point. forwardwalking.cpp line 68
    t_walker_->addPropertyHistoryPoint(p_ids_[i], P.PropertyList[h_ids_[i]]);
    // for the current walker
    // find the id of the property of interest in property history
    int j       = 0;   // counts the steps for this walker to go back
    int FWindex = t_walker_->PHindex[p_ids_[i]] - 1;  // this is the current walkers index for a given property
    //create iterator for 
    // loop over nstep configurations for this property
    // while you are less than the recorded values and less than the desired number of steps in the past
    while (j < walker_lengths_[i][1] & j != nstep){
        // take a step back in history by blockfreq to get to next recored value
        FWindex -= walker_lengths_[i][0];
        if (FWindex < 0)
            FWindex += walker_lengths_[i][2];  // is this trying to  exit the loop essentially? not sure.
        (*Vit) = t_walker_->PropertyHistory[p_ids_[i]][FWindex]
        j++;
        Vit++;
    }
        
    }

            // call to calc_gdd and calc_b 

    }
    
}