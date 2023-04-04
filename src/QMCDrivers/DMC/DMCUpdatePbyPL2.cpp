//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include "QMCDrivers/DMC/DMCUpdatePbyPL2.h"
#include "Particle/MCWalkerConfiguration.h"
// #include "Particle/DistanceTable.h"
#include "Particle/HDFWalkerIO.h"
#include "ParticleBase/ParticleUtility.h"
#include "ParticleBase/RandomSeqGenerator.h"
#include "QMCDrivers/DriftOperators.h"
#if !defined(REMOVE_TRACEMANAGER)
#include "Estimators/TraceManager.h"
#else
using TraceManager = int;
#endif
//#define TEST_INNERBRANCH
#include "QMCDrivers/DMC/DMCUpdatePbyP.h"

namespace qmcplusplus
{
using WP = WalkerProperties::Indexes;

/// Constructor.
DMCUpdatePbyPL2::DMCUpdatePbyPL2(MCWalkerConfiguration& w,
                                 TrialWaveFunction& psi,
                                 QMCHamiltonian& h,
                                 RandomGenerator& rg)
    : QMCUpdateBase(w, psi, h, rg)
{
  setup_timers(myTimers, DMCTimerNames, timer_level_medium);
}

/// destructor
DMCUpdatePbyPL2::~DMCUpdatePbyPL2() {}

void DMCUpdatePbyPL2::advanceWalker(Walker_t& thisWalker, bool recompute)
{
  Walker_t::WFBuffer_t& w_buffer(thisWalker.DataSet);
  {
    ScopedTimer local_timer(myTimers[DMC_buffer]);
    W.loadWalker(thisWalker, true);
    Psi.copyFromBuffer(W, w_buffer);
  }
  //create a 3N-Dimensional Gaussian with variance=1
  makeGaussRandomWithEngine(deltaR, RandomGen);
  int nAcceptTemp(0);
  int nRejectTemp(0);
  //copy the old energy and scale factor of drift
  //EstimatorRealType eold(thisWalker.Properties(LOCALENERGY));
  //EstimatorRealType enew(eold);
  FullPrecRealType eold(thisWalker.Properties(WP::LOCALENERGY));
  FullPrecRealType enew(eold);
  RealType rr_proposed = 0.0;
  RealType rr_accepted = 0.0;
  mPosType K;
  mTensorType D;
  mTensorType Dchol;
  PosType Ktmp, drtmp;
  TensorType Dtmp;
  bool L2_proj = H.has_L2();
  if (L2_proj)
  {
    Ktmp = 0.0;
    Dtmp = 0.0;
    for (int d = 0; d < DIM; d++)
      Dtmp(d, d) = 1.0;
  }
  {
    ScopedTimer local_timer(myTimers[DMC_movePbyP]);
    for (int ig = 0; ig < W.groups(); ++ig) //loop over species
    {
      RealType tauovermass = Tau * MassInvS[ig];
      RealType oneover2tau = 0.5 / (tauovermass);
      RealType sqrttau     = std::sqrt(tauovermass);
      RealType rr;
      for (int iat = W.first(ig); iat < W.last(ig); ++iat)
      {
        //W.setActive(iat);
        //get the displacement
        GradType grad_iat = Psi.evalGrad(W, iat);
        mPosType dr;
        mPosType dr_diff = deltaR[iat];
        if (!L2_proj) // normal projector
        {
          getScaledDrift(tauovermass, grad_iat, dr);
          dr += sqrttau * dr_diff;
          rr = tauovermass * dot(dr_diff, dr_diff);
          rr_proposed += rr;
          if (rr > m_r2max)
          {
            ++nRejectTemp;
            W.accept_rejectMove(iat, false);
            continue;
          }
          if (!W.makeMoveAndCheck(iat, dr))
          {
            W.accept_rejectMove(iat, false);
            continue;
          }
        }
        else // projector including L2 potential
        {
          // do a fake move (zero distance)
          // this ensures the temporary distance data is correct
          // will need to remove this later, but requires reimplementation of computeL2DK
          dr = 0.0;
          if (!W.makeMoveAndCheck(iat, dr))
          {
            W.accept_rejectMove(iat, false);
            continue;
          }

          H.computeL2DK(W, iat, Dtmp, Ktmp);
          D = Dtmp; // upcast for mixed precision
          K = Ktmp;
          getScaledDriftL2(tauovermass, grad_iat, D, K, dr);

          W.accept_rejectMove(iat, false);
          rr = tauovermass * dot(dr_diff, dr_diff);
          rr_proposed += rr;
          if (rr > m_r2max)
          {
            ++nRejectTemp;
            W.accept_rejectMove(iat, false);
            continue;
          }

          // move with just drift to update distance tables
          if (!W.makeMoveAndCheck(iat, dr))
          {
            W.accept_rejectMove(iat, false);
            continue;
          }

          // compute diffusion step
          H.computeL2D(W, iat, Dtmp);
          D       = Dtmp; // upcast for mixed precision
          Dchol   = cholesky(D);
          dr_diff = dot(Dchol, dr_diff);
          dr += sqrttau * dr_diff;

          // reverse the intermediate drift move
          W.accept_rejectMove(iat, false);
          // move with drift and diffusion together
          if (!W.makeMoveAndCheck(iat, dr))
          {
            W.accept_rejectMove(iat, false);
            continue;
          }
        }
        ValueType ratio = Psi.calcRatioGrad(W, iat, grad_iat);
        //node is crossed reject the move
        if (branchEngine->phaseChanged(Psi.getPhaseDiff()))
        {
          ++nRejectTemp;
          ++nNodeCrossing;
          W.accept_rejectMove(iat, false);
          Psi.rejectMove(iat);
        }
        else
        {
          FullPrecRealType logGf = -0.5 * dot(deltaR[iat], deltaR[iat]);
          //Use the force of the particle iat
          DriftModifier->getDrift(tauovermass, grad_iat, drtmp);
          dr                     = drtmp; // upcast for mixed precision
          dr                     = W.R[iat] - W.getActivePos() - dr;
          FullPrecRealType logGb = -oneover2tau * dot(dr, dr);
          RealType prob          = std::norm(ratio) * std::exp(logGb - logGf);
          bool is_accepted       = false;

          if (RandomGen() < prob)
          {
            is_accepted = true;

            ++nAcceptTemp;
            Psi.acceptMove(W, iat, true);
            rr_accepted += rr;
            // calculate shared observable here for new and old positions. 
            ValueType gr_rp;
            ValueType
            gr_rp = calculate_greens_function();
            thisWalker.addPropertyHistoryPoint()
          }
          else
          {
            ++nRejectTemp;
            Psi.rejectMove(iat);
          }
          W.accept_rejectMove(iat, is_accepted);
        }
      }
    }
    Psi.completeUpdates();
    W.donePbyP();
  }

  if (nAcceptTemp > 0)
  {
    //need to overwrite the walker properties
    RealType logpsi(0);
    {
      ScopedTimer local_timer(myTimers[DMC_buffer]);
      thisWalker.Age = 0;
      logpsi         = Psi.updateBuffer(W, w_buffer, recompute);
      W.saveWalker(thisWalker);
    }
    {
      ScopedTimer local_timer(myTimers[DMC_hamiltonian]);
      enew = H.evaluateWithToperator(W);
    }
    thisWalker.resetProperty(logpsi, Psi.getPhase(), enew, rr_accepted, rr_proposed, 1.0);
    thisWalker.Weight *= branchEngine->branchWeight(enew, eold);
    {
      ScopedTimer local_timer(myTimers[DMC_collectables]);
      H.auxHevaluate(W, thisWalker);
      H.saveProperty(thisWalker.getPropertyBase());
    }
  }
  else
  {
    //all moves are rejected: does not happen normally with reasonable wavefunctions
    thisWalker.Age++;
    thisWalker.Properties(WP::R2ACCEPTED) = 0.0;
    //weight is set to 0 for traces
    // consistent w/ no evaluate/auxHevaluate
    RealType wtmp     = thisWalker.Weight;
    thisWalker.Weight = 0.0;
    H.rejectedMove(W, thisWalker);
    thisWalker.Weight = wtmp;
    ++nAllRejected;
    enew = eold; //copy back old energy
    thisWalker.Weight *= branchEngine->branchWeight(enew, eold);
  }
#if !defined(REMOVE_TRACEMANAGER)
  Traces->buffer_sample(W.current_step);
#endif
  {
    ScopedTimer local_timer(myTimers[DMC_tmoves]);
    const int NonLocalMoveAcceptedTemp = H.makeNonLocalMoves(W);
    if (NonLocalMoveAcceptedTemp > 0)
    {
      RealType logpsi = Psi.updateBuffer(W, w_buffer, false);
      W.saveWalker(thisWalker);
      NonLocalMoveAccepted += NonLocalMoveAcceptedTemp;
    }
  }
  nAccept += nAcceptTemp;
  nReject += nRejectTemp;

  setMultiplicity(thisWalker);
}


//a function to do all the work for the observable without using existing code. will probably be able to clean and resuse previously calculated things
//accepted: boolean if move was accepted or not.
// parameters:
// prob: probaility that move would be accepted/rejected
//
//TODO: move this into observable class that ill be wrapped by forward walking, but called here. 
return_t calculate_greens_function(bool accepted, RealType prob, ParticleGradient& G, R_curr_,R_old){
  RealType gf //greens function
  RealType  T // transition probability
  RealType  W // branching factor.
  RealType  V // deltalnPsi in ACFOrce this was passed as a ParticleGradient G
  RealType  F // dampinf factor of volcities divergence
  RealType  S  // term for branchind depending on current energy estimate E_L  by F
  //TODO:  may need to loop over dim or understand what the pos data type is.  
  // calculate V
  V = G[iat];

  // calculate F
  F = sqrt((2*V*V*tau)-1)/(V*V*t);

  //Calculate T:
  //TODO: this can benefit from log Gb and logGf, maybe?

  T = e^[-R -R -F(R)V(R)*tau];
  
  // Calculate W and components.
  // TODO: E_L_R definitely exits within Walker Properties so can pull ou tthere and pass R' values from above calc.
  //TODO: E_est is current best esitmate of the energy, which must exists within this driver for R, not sure how to do that for R' 
  //Calculate S
  Realtype E_est_R;
  Realtype E_L_R;
  Realtype E_est_R_prime;
  Realtype E_L_R_prime;
  int N_curr;
  int N_target;

  S_R = (E_est_R - E_L_R)*F=ln(N_curr/N_target);
  S_R_prime = (E_est_R_prime -E_L_R_prime)*F

  //Calculate W
  W = exp((S_R_prime + S_R)*tau/2);
  

  //Calculate W
  if accepted{
    gf = T*prob*W;
  }
  else{
    gf = T*(1-prob)*W;
  }

  // actual expressions to store in property history:
  RealType dlambda_E_L = 0;
  RealType dlambda_G = 0;
  dlambda_E_zero = ;
  dlambda_E_L = ;

  int dlambda_E_L_prop_id;
  int dlambda_G_prop_id;
  thisWalker.addPropertyHistoryPoint(dlambda_E_L_prop_id, dlambda_E_L);
  thisWalker.addPropertyHistoryPoint(dlambda_G_prop_id,dlambda_G);




}

} // namespace qmcplusplus
