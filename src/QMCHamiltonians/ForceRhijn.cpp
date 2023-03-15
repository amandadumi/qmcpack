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
#include "QMCHamiltonians/QMCHamiltonian.h"
#include "QMCDrivers/WalkerProperties.h"
#include <cstdlib>
#include <set>
#include <string>



namespace qmcplusplus
{

ForceRhijn::ForceRhijn()
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
        // first_hamiltonian_ = h.startIndex();

        int numProps = P.PropertyList.Names.size();
        std::string tagName = "walker_history_hash";
        walker_lengths_.push_back(nstep);
        // add property history to record the id of each walker step
        int pindx = P.addPropertyHistory(nstep);
        p_ids_.push_back(pindx);
        return true;

    }

bool ForceRhijn::get(std::ostream& os) const
    {
    os << "ForceRhijn";
    return true;
    }



ForceRhijn::Return_t ForceRhijn::evaluate(ParticleSet& P)
{
    // have walker tracker be the size of this walkers id if not already in list.
    std::vector<RealType>::iterator Vit = values_.begin();

    while (walker_tracker.size() < t_walker_->ID){
        walker_tracker.push_back(0);
    }
    // increment the element by one since we are addint to property history for it.
    walker_tracker[t_walker_->ID] += 1;
    for int j=0; j < nsteps, j++){
        int i = 0;
        t_walker_->addPropertyHistoryPoint(p_ids_[i], walker_tracker[t_walker_->ID]);
    }
    //the property will keep track of the column it needs to fill through phindex
    // app_log() << "Not a valid H element(" << Hindex << ") Valid names are:";

    return 0.0;
}


ForceRhijn::Return_t ForceRhijn::rejectedMove(ParticleSet& P){return 0.0;}


void ForceRhijn::addObservables(PropertySetType& plist)
{
  //not used
}

void ForceRhijn::addObservables(PropertySetType& plist, BufferType& collectables)
{
  my_index_ = plist.size();
  int nc    = 0;
    for (int j = 0; j < nstep; ++j)
    {
      std::stringstream sstr;
      sstr << "FR_" << "property_history" << "_" << j;
      int id = plist.add(sstr.str());
    }
  app_log() << "ForceRhijn::Observables [" << my_index_ << ", " << my_index_ + nstep << ")" << std::endl;
}


void ForceRhijn::setObservables(PropertySetType& plist)
{
  copy(values_.begin(), values_.end(), plist.begin() + my_index_);
}

void ForceRhijn::setParticlePropertyList(PropertySetType& plist, int offset)
{
  copy(values_.begin(), values_.end(), plist.begin() + my_index_ + offset);
}
}