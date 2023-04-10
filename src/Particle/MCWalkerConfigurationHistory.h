//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


/*
 * Reptile.h
 *
 * This class is a wrapper/utility class for groups of walkers within
 * MCWalkerConfiguration for Reptation Monte Carlo.
 *
 * Walkers are stored as usual in MCWalkerConfiguration.  However, this interface
 * stores the walker histories and implements functionality needed to evaluate forces with history
 */

#ifndef QMCPLUSPLUS_REPTILE_H
#define QMCPLUSPLUS_REPTILE_H

#include "QMCDrivers/DriftOperators.h"
#include "QMCDrivers/WalkerProperties.h"
#include "Configuration.h"
#include "Walker.h"

namespace qmcplusplus
{
class MCWalkerConfiguration;

class MCWalkerConfigurationHistory : public QMCTraits
{
  public:
  using WP       = WalkerProperties::Indexes;
  using Walker_t = MCWalkerConfiguration::Walker_t;
  //using Buffer_t = Walker_t::Buffer_t             ;
  //    using Walker_t = MCWalkerConfiguration::Walker_t;
  using WalkerIter_t    = MCWalkerConfiguration::iterator;
  using ReptileConfig_t = std::vector<Walker_t::ParticlePos>;

  int nhist; // the number of steps to record in history.
  std::vector<IndexType> RejectedPos; // tiny vector to store the position of a rejected move. will have to be assigned to a walker i think.
  RejectedPos.resize(3);
  IndexType acc_or_rej;
  int nhist; // the number of steps to record in history.
  

  MCWalkerConfiguration& w;
  WalkerIter_t prev_timestep_start, prev_timestep_end;
  WalkerIter_t curr_timestep_start;
  IndexType timeindex;
  Walker_t* currentfirstwalker; // the firs walker in the set at the current timestep. 
  IndexType curr_timestep; //can we get this from the size of walker configuration?

  inline History(MCWalkerConfiguration& W, WalkerIter_t start, WalkerIter_t end)
      : w(W),
        timestep_start(start),
        timestep_end(end),
        prophead(0) //, r2prop(0.0), r2accept(0.0),tau(0.0)
  {
    //TODO: find out if you can store a position as a walker property. 
    // actually this is still wrong because it will be single x value intead of vector for each particle. hm.
    RejectedPos[0] = w.addProperty("RejectedPos_x");
    RejectedPos[1] = w.addProperty("RejectedPos_y");
    RejectedPos[2] = w.addProperty("RejectedPos_z");
    acc_or_rej = w.addProperty("acc_or_rej");

  }

  ~History() {}
  // return the number of timesteps back in history that will be recorded
  inline IndexType size() { return nhist; }
  // returns a walker from the array
  inline Walker_t& operator[](IndexType i) { return getWalker(getwalkeratcurrtimeIndex(i)); }

  inline IndexType wrapIndex(IndexType repindex) { return (repindex % nbeads + nbeads) % nbeads; }

  inline Walker_t& getWalker(IndexType i)
  {
    WalkerIter_t walker = curr_timestep + i;
    return **walker;
  }
  //get a walker at a current time step?
  IndexType curr_timestep; // current time step in regards to place in nhist, not place globally
  inline void update_timestep_counter(){
    if curr_timestep == nhist{
      curr_timestep = 0;
    }
    else{curr_timestep++;}
  }

  inline IndexType getcurrtimestepIdx() { return walkers_per_step*curr_timestep); }
  inline IndexType getprevioustimestepIdx(int nsteps_back) { return walkers_per_step*(curr_timestep-nsteps_back)); }
  inline IndexType getwalkeratcurrtimeIndex(IndexType i) { return getcurrtimestepIdx()*numwalkersperstep + i); }
  inline Walker_t& getwalkeratprevtimeIndex(IndexType i, Indextype nsteps_back) { return getprevioustimestepIdx(nsteps_back)*numwalkerperstep +i); }

  // for a step copyt the current walkers as the beginning of next step. 
  inline Walker_t& set_initial_walker_state(){
    //if first create walkers
     
    // if not first add walkers?
    //todo where are iterators set?
    curr_timestep_start = curr_timestep*walkers_per_step
    W.copyWalkers(prev_timestep_start, prev_timestep_end,curr_timestep_start);
  }

  inline void printState()
  {
    app_log() << "********PRINT WALKER HISTORY INFORMATION*********\n";
    app_log() << "BeadIndex\tWrapIndex\tEnergy\tAction[0]\tAction[1]\tAction[2]\t\n";
    for (int i = 0; i < nbeads; i++)
    {
      app_log() << i << "\t" << getBeadIndex(i) << "\t" << getBead(i).Properties(WP::LOCALENERGY) << "\t"
                << getBead(i).Properties(Action[0]) << "\t" << getBead(i).Properties(Action[1]) << "\t"
                << getBead(i).Properties(Action[2]) << "\n";
    }
    app_log() << "POSITIONS===============:\n";
    for (int i = 0; i < nbeads; i++)
    {
      //  app_log()<<i<<"\t1"<<1<<"\t"<<getBead(i).R[0]<<"\n";
      //  app_log()<<i<<"\t2"<<2<<"\t"<<getBead(i).R[1]<<"\n";
      app_log() << "BEAD #" << i << " tau = " << tau * i << std::endl;
      app_log() << getBead(i).R << std::endl;
    }
    app_log() << "GVECS===============:\n";
    for (int i = 0; i < nbeads; i++)
    {
      //      app_log()<<i<<"\t1"<<1<<"\t"<<getBead(i).G[0]<<"\n";
      //      app_log()<<i<<"\t2"<<2<<"\t"<<getBead(i).G[1]<<"\n";
      app_log() << "BEAD #" << i << " tau = " << tau * i << std::endl;
      app_log() << getBead(i).G << std::endl;
    }
    app_log() << "************************************\n";
  }
  inline RealType getTau() { return tau; }
  inline void setTau(RealType t) { tau = t; }


  //This takes a value of imaginary time "t" and returns a 3N particle position vector, corresponding to a time slice extrapolated
  // from the current reptile.  If t>length of reptile, then return the last bead.  if t<0; return the first bead.
  inline Walker_t::ParticlePos linearInterp(RealType t)
  {
    IndexType nbead =
        IndexType(t / tau); //Calculate the lower bound on the timeslice.  t is between binnum*Tau and (binnum+1)Tau
    RealType beadfrac = t / tau - nbead; //the fractional coordinate between n and n+1 bead
    if (nbead <= 0)
    {
      ParticleSet::ParticlePos result = getHead().R;
      return result;
    }
    else if (nbead >= nbeads - 1)
    {
      ParticleSet::ParticlePos result = getTail().R;
      return result;
    }

    else
    {
      Walker_t::ParticlePos dR(getBead(nbead + 1).R), interpR(getBead(nbead).R);
      dR = dR - getBead(nbead).R;

      interpR = getBead(nbead).R + beadfrac * dR;
      return interpR;
    }
  }
  inline ReptileConfig_t getReptileSlicePositions(RealType tau, RealType beta)
  {
    IndexType nbeads_new = IndexType(beta / tau);
    ReptileConfig_t new_reptile_coords(0);

    for (IndexType i = 0; i < nbeads_new; i++)
      new_reptile_coords.push_back(linearInterp(tau * i));

    return new_reptile_coords;
  }

  inline void setReptileSlicePositions(ReptileConfig_t& rept)
  {
    if (rept.size() == nbeads)
    {
      for (int i = 0; i < nbeads; i++)
        getBead(i).R = rept[i];
    }
    else
      ;
  }

  inline void setReptileSlicePositions(Walker_t::ParticlePos R)
  {
    for (int i = 0; i < nbeads; i++)
      getBead(i).R = R;
  }
};


} // namespace qmcplusplus
#endif
