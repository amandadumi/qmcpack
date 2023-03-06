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

bool ForceRhijn::putSpecial(xmlNodePtr cur, QMCHamiltonian& h, ParticleSet& P)
    {   
        using WP = WalkerProperties::Indexes;
        first_hamiltonian_ = h.startIndex();
        my_index_ = plist.size()

        int numProps = P.PropertyList.Names.size();
        int Hindex  = WP::LOCALPOTENTIAL;
        std::string tagName = "LocPot";
        std::vector<int> pms(3);
        pms[0] = blockFreq;
        pms[1] = numT;
        pms[2] = blockSeries + 2;
        walker_lengths_.push_back(pms);
        int maxWsize = blockSeries + 2;
        int pindx    = P.addPropertyHistory(maxWsize);
        p_ids_.push_back(pindx);
    }

bool ForceRhijn::get(std::ostream& os) const
    {
    os << "ForceRhijn";
    return true;
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

void ForceRhijn::calculate_gb(){
    //pv e^{1/2(E_L(R')+ E_L(R) + E_T))}

    //E_L (R')
    float locale_rprime = 0; 

    //E_L(R)
    float locale_r = 0;
    getLocalEnergy
    //E_T
    trial_energy = Walker; 
    }

void ForceRhijn::evaluate(ParticleSet& P){
    int i = Hindex; // where observable of local energy is in hamiltonian.
    t_walker_->addPropertyHistoryPoint(pindx, P.PropertyList[i]);
    // for the current walker
    // find the id of the property of interest in property history
    int j       = 0;   // counts the steps for this walker to go back
    int FWindex = t_walker_->PHindex[pindx] - 1;  // this is the current walkers index for a given property (p_id[i])
    // loop over nstep configurations for this property
    // while you are less than the recorded values and less than the desired number of steps in the past
    while (j < walker_lengths_[i][1] & j != nstep){
        // take a step back in history by blockfreq to get to next recored value
        FWindex -= walker_lengths_[i][0];
        if (FWindex < 0)
            FWindex += walker_lengths_[i][2];  // is this trying to  exit the loop essentially? not sure.
        (*Vit) = t_walker_->PropertyHistory[pindx][FWindex]
        j++;
        Vit++;
        //accumulate forces here
        forces_ = 
    }
    copy(values_.begin(), values_.end(), t_walker_->getPropertyBase() + first_hamiltonian_ + my_index_); //todo: change first_hamiltonian and my_index as these are taken straight from forward walking
    }

            // call to calc_gdd and calc_b 


ForceRhijn::Return_t ForceRhijn::rejectedMove(ParticleSet& P)
{
  for (int i = 0; i < nobservables_; i++)
  {
    int lastindex = t_walker_->PHindex[p_ids_[i]] - 1;
    if (lastindex < 0)
      lastindex += walker_lengths_[i][2];
    t_walker_->addPropertyHistoryPoint(p_ids_[i], t_walker_->PropertyHistory[p_ids_[i]][lastindex]);
  }
  calculate(P);
  return 0.0;
}
    }
    
}

