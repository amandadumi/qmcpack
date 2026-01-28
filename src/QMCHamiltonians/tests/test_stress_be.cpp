//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Yubo "Paul" Yang, yyang173@illinois.edu, University of Illinois at Urbana-Champaign
//
// File created by: Yubo "Paul" Yang, yyang173@illinois.edu, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"

#include "OhmmsData/Libxml2Doc.h"
#include "OhmmsPETE/OhmmsMatrix.h"
#include "Particle/ParticleSet.h"
#include "LongRange/EwaldHandler3D.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "QMCHamiltonians/StressPBC.h"
#include "Utilities/RuntimeOptions.h"

#include <stdio.h>
#include <string>

using std::string;

namespace qmcplusplus
{
// PBC case
TEST_CASE("Stress BCC Be Ewald3D", "[hamiltonian]")
{
  std::cout << "in Be stress test" << std::endl;

  Lattice lattice;
  lattice.BoxBConds = true; // periodic
  lattice.R(0, 0) = 4.25263466;
  lattice.R(0, 1) = 0.00000000;
  lattice.R(0, 2) = 0.00000000;
  lattice.R(1, 0) = -2.12631733;
  lattice.R(1, 1) = 3.68288965;
  lattice.R(1, 2) = 0.00000000;
  lattice.R(2, 0) = 0.00000000;
  lattice.R(2, 1) = 0.00000000;
  lattice.R(2, 2) = 6.73797666;
  lattice.LR_dim_cutoff = 40;
  lattice.reset();


  const SimulationCell simulation_cell(lattice);
  ParticleSet ions(simulation_cell);
  ParticleSet elec(simulation_cell);

  ions.setName("ion0");
  ions.create({2});
  SpeciesSet& tspecies             = ions.getSpeciesSet();
  int pIdx                         = tspecies.addSpecies("Be");
  int pchargeIdx                   = tspecies.addAttribute("charge");
  int patomicnumberIdx             = tspecies.addAttribute("atomicnumber");
  tspecies(pchargeIdx, pIdx)       = 4;
  tspecies(patomicnumberIdx, pIdx) = 4;
  ions.resetGroups();
  ions.R[0]                        = {0.0, 2.45523521, 5.05348250};
  ions.R[1]                        = {2.12628905,1.22762976,1.68449417};
  ions.createSK();
  const int ii_dist_idx = ions.addTable(ions);
  ions.update();
  double vinv = -1. / ions.getLattice().Volume;

  elec.setName("e");
  elec.create({4, 4});
  SpeciesSet& especies       = elec.getSpeciesSet();
  int upIdx                  = especies.addSpecies("u");
  int chargeIdx              = especies.addAttribute("charge");
  int massIdx                = especies.addAttribute("mass");
  especies(chargeIdx, upIdx) = -1;
  especies(massIdx, upIdx)   = 1.0;
  int dnIdx                  = especies.addSpecies("d");
  especies(chargeIdx, dnIdx) = -1;
  especies(massIdx, dnIdx)   = 1.0;
  elec.resetGroups(); // need to set Mass so
  elec.R[0] = {0.1, 0.2, 0.3};
  elec.R[1] = {1.0, 1.0, 1.0};
  elec.R[2] = {2.0, 2.0, 2.0};
  elec.R[3] = {0.8, 0.7, 0.5};
  elec.R[4] = {1.1, 1.2, 1.3};
  elec.R[5] = {1.0, 1.6, 1.4};
  elec.R[6] = {2.1, 0.2, 0.3};
  elec.R[7] = {1.3, 2.0, 1.5};
  elec.createSK();
  const int ee_dist_idx = elec.addTable(elec);
  const int ei_dist_idx = elec.addTable(ions);
  elec.update();

  
  
  RuntimeOptions runtime_options;
  TrialWaveFunction psi(runtime_options);


  LRCoulombSingleton::CoulombHandler = std::make_unique<EwaldHandler3D>(ions);
  LRCoulombSingleton::CoulombHandler->initBreakup(ions);
  LRCoulombSingleton::CoulombDerivHandler = std::make_unique<EwaldHandler3D>(ions);
  LRCoulombSingleton::CoulombDerivHandler->initBreakup(ions);
  StressPBC est(ions, elec, psi);

  std::cout << "ee_00 stresses before evaluation are "<<std::endl;
  std::cout << est.getStressEE()(0,0) <<std::endl;
  std::cout << "ii_00 stresses before evaluation are "<<std::endl;
  std::cout << est.getStressIonIon()(0,0) <<std::endl;
  
  est.evaluate(elec);

  // i-i = e-e stress is validated against Quantum Espresso's ewald method
  //  they are alseo double checked using finite-difference
  CHECK(est.getStressIonIon()(0, 0)*vinv == Approx(-0.03940591355330148));
  CHECK(est.getStressIonIon()(0, 1) == Approx(0.0));
  CHECK(est.getStressIonIon()(0, 2) == Approx(0.0));
//  CHECK(est.getStressEE()(1, 0) == Approx(0.0));
//  CHECK(est.getStressEE()(1, 1) == Approx(-0.039406049512414824));
//  CHECK(est.getStressEE()(1, 2) == Approx(0.0));
//  CHECK(est.getStressEE()(2, 0) == Approx(0.0));
//  CHECK(est.getStressEE()(2, 1) == Approx(0.0));
//  CHECK(est.getStressEE()(2, 2) == Approx(-0.03940734112399154));

  //Electron-Ion stress diagonal is internally validated using fd.
//  CHECK(est.getStressEI()(0, 0) == Approx(-0.00745376));
  //CHECK(est.getStressEI()(0, 1) == Approx(0.0));
  //CHECK(est.getStressEI()(0, 2) == Approx(0.0));
  //CHECK(est.getStressEI()(1, 0) == Approx(0.0));
//  CHECK(est.getStressEI()(1, 1) == Approx(-0.00745376));
  //CHECK(est.getStressEI()(1, 2) == Approx(0.0));
  //CHECK(est.getStressEI()(2, 0) == Approx(0.0));
  //CHECK(est.getStressEI()(2, 1) == Approx(0.0));
//  CHECK(est.getStressEI()(2, 2) == Approx(-0.00745376));

  LRCoulombSingleton::CoulombHandler.reset(nullptr);

  LRCoulombSingleton::CoulombDerivHandler.reset(nullptr);
}
} // namespace qmcplusplus
