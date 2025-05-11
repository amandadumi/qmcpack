//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by:  Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"
#include "Message/Communicate.h"
#include "Particle/ParticleSet.h"
#include "Particle/DistanceTable.h"
#include "QMCWaveFunctions/Jastrow/SNAJastrow.h"
#include "QMCWaveFunctions/Jastrow/SNAJastrowBuilder.h"
#include "OhmmsData/Libxml2Doc.h"
//lammps libraries
#include <stdio.h>
#include <string>
#include <cstring>

using std::string;
namespace qmcplusplus
{

TEST_CASE("sna_jastrow_update_rij", "[wavefunction]")
{
// short input xml check that lammps positions are correct.
  const SimulationCell simulation_cell;
  ParticleSet ions(simulation_cell), electrons(simulation_cell);

  electrons.create({1,1});
  electrons.setName("e_u");
  electrons.R[0][0] = 0.2;
  electrons.R[0][1] = 0.2;
  electrons.R[0][2] = 0.2;
  electrons.R[1][0] = 0.1;
  electrons.R[1][1] = 0.1;
  electrons.R[1][2] = 0.1;

  ions.create({1});
  ions.setName("ions");
  SpeciesSet& tspecies  = ions.getSpeciesSet();
  int ion_a  = tspecies.addSpecies("H");
  ions.R[0][0] = 0.0;
  ions.R[0][1] = 0.0;
  ions.R[0][2] = 0.0;
  ions.update();
  int ee_table = electrons.addTable(electrons);
  int ei_table = electrons.addTable(ions);
  electrons.update();

  auto jas = std::make_unique<SNAJastrow>(std::string("snap"),ions,electrons,std::string("linear"), 2, 7.0);
  std::cout<< "initialized object" <<std::endl;
  //jas->evaluateLog(electrons, electrons.G, electrons.L);
  jas->update_sna_rij(electrons,0,false);
  double true_dist = 0.1;
  double internal_dist = electrons.getDistTableAA(ee_table).getDisplRow(1)[0][0];
  REQUIRE(jas->sna_desc.rij[0][0] == Approx(true_dist));//elec1
  REQUIRE(jas->sna_desc.rij[0][0] == Approx(internal_dist));//elec1
  // what about from electron 0 to ion 0?
  true_dist = 0.2;
  internal_dist = -electrons.getDistTableAB(ei_table).getDisplRow(0)[0][0];
  REQUIRE(jas->sna_desc.rij[1][0] == Approx(true_dist));//elec1
  REQUIRE(jas->sna_desc.rij[1][0] == Approx(internal_dist));//elec1
  //update rij to reflect the ion
  jas->update_sna_rij(electrons,2,false);
  true_dist = -0.2;
  internal_dist = electrons.getDistTableAB(ei_table).getDisplRow(0)[0][0];
  REQUIRE(jas->sna_desc.rij[0][0] == Approx(true_dist));//elec1
  REQUIRE(jas->sna_desc.rij[0][0] == Approx(internal_dist));//elec1

};

TEST_CASE("snap_jastrow_init_with_coeff", "[wavefunction]")
{
  Communicate* c = OHMMS::Controller;
  std::cout<< "starting test" <<std::endl;
//   // short input xml check that lammps positions are correct.
  const SimulationCell simulation_cell;
  ParticleSet ions(simulation_cell), electrons(simulation_cell);

  electrons.create({1,1});
  electrons.setName("e");
  SpeciesSet& especies  = electrons.getSpeciesSet();
  int elec_up  = especies.addSpecies("e_u");
  int elec_down  = especies.addSpecies("e_d");
  electrons.R[0][0] = 0.2;
  electrons.R[0][1] = 0.2;
  electrons.R[0][2] = 0.2;
  electrons.R[1][0] = 0.1;
  electrons.R[1][1] = 0.1;
  electrons.R[1][2] = 0.1;

  ions.create({1});
  ions.setName("ions");
  SpeciesSet& tspecies  = ions.getSpeciesSet();
  int ion_a  = tspecies.addSpecies("He");

  ions.R[0][0] = 0.0;
  ions.R[0][1] = 0.0;
  ions.R[0][2] = 0.0;
  ions.update();
  electrons.update();



  const char * xmltext = R"XML(<tmp>
  <wavefunction name="psi0" target="e">
  <jastrow name="snap" type="snap" function="snap">
    <correlation>
     <coefficients id="e_u" type="Array"> 0.1 0.2 0.3 0.4 0.5  </coefficients>
    </correlation>
    <correlation>
     <coefficients id="e_d" type="Array"> 1.1 1.2 1.3 1.4 1.5 </coefficients>
    </correlation>
    <correlation>
     <coefficients id="He" type="Array"> 2.1 2.2 2.3 2.4 2.5 </coefficients>
    </correlation>
  </jastrow>
</wavefunction>
</tmp>)XML";

  Libxml2Document doc;
  bool okay;
  okay    = doc.parseFromString(xmltext);
  REQUIRE(okay);

  xmlNodePtr root = doc.getRoot();

  xmlNodePtr jas_node = xmlFirstElementChild(root);
  xmlNodePtr corr_node = xmlFirstElementChild(jas_node);
  SNAJastrowBuilder SJBuilder(c,electrons,ions);
  //TODO: something seems weird in the manual parsing of coeffs, however, this is not something i need so skipping this test for now
   auto sj_uptr = SJBuilder.buildComponent(corr_node);
   SNAJastrow* sj = static_cast<SNAJastrow*>(sj_uptr.get());
   REQUIRE(sj->snap_beta[0][0] == 0.1);
   REQUIRE(sj->snap_beta[0][4] == 0.5);
   REQUIRE(sj->snap_beta[1][4] == 1.5);
   REQUIRE(sj->snap_beta[2][4] == 2.5); 
};

TEST_CASE("snap_jastrow_check_bispectrum", "[wavefunction]")
{
  Communicate* c = OHMMS::Controller;
  std::cout<< "starting snap_jastrow_check_bispectrum" <<std::endl;
//   // short input xml check that lammps positions are correct.
  const SimulationCell simulation_cell;
  ParticleSet ions(simulation_cell), electrons(simulation_cell);

  electrons.setName("e");
  electrons.create({1,1});
  SpeciesSet& especies  = electrons.getSpeciesSet();
  int elec_up  = especies.addSpecies("e_u");
  int elec_down  = especies.addSpecies("e_d");
  electrons.R[0][0] = 0.2;
  electrons.R[0][1] = 0.2;
  electrons.R[0][2] = 0.2;
  electrons.R[1][0] = 0.1;
  electrons.R[1][1] = 0.1;
  electrons.R[1][2] = 0.1;

  ions.create({1});
  ions.setName("ions");
  SpeciesSet& tspecies  = ions.getSpeciesSet();
  int ion_a  = tspecies.addSpecies("H");
  ions.R[0][0] = 0.0;
  ions.R[0][1] = 0.0;
  ions.R[0][2] = 0.0;
  ions.update();

  int ee_table = electrons.addTable(electrons);
  int ei_table = electrons.addTable(ions);
  electrons.update();
  
 std::vector<double> true_bispectrum_e1 = {26.949,53.8442,80.6591,80.6323,80.5251};
 std::vector<double> true_bispectrum_e2= {26.9796,53.9054,80.7507,80.7239,80.6166};
 std::vector<double> true_bispectrum_ion = {26.949,53.8442,80.6591,80.6323,80.5251};



  auto jas = std::make_unique<SNAJastrow>(std::string("snap"), ions, electrons,std::string("linear"), 2, 7.0);
  jas->update_sna_rij(electrons,0,false);
  double true_dist = 0.1;
  double internal_dist = electrons.getDistTableAA(ee_table).getDisplRow(1)[0][0];
  REQUIRE(jas->sna_desc.rij[0][0] == Approx(internal_dist));//elec1
  REQUIRE(jas->sna_desc.rij[0][0] == Approx(true_dist));//elec1
  std::cout<< "snap_jastrow_check_bispectrum before compute_bispectrum" <<std::endl;
  jas->compute_bispectrum(0);
  REQUIRE(jas->sna[0][0] == Approx(true_bispectrum_e1[0]));//elec1
  jas->update_sna_rij(electrons,1,false);
  jas->compute_bispectrum(1);
  REQUIRE(jas->sna[1][1] == Approx(true_bispectrum_e2[1]));//elec1
  jas->update_sna_rij(electrons,2,false);
  jas->compute_bispectrum(2);
  REQUIRE(jas->sna[2][2] == Approx(true_bispectrum_ion[2]));//elec1
}

TEST_CASE("snap_jastrow_check_d_dr_bispectrum", "[wavefunction]")
{
  Communicate* c = OHMMS::Controller;
  std::cout<< "starting snap_jastrow_check_d_dr_bispectrum" <<std::endl;
  const SimulationCell simulation_cell;
  ParticleSet ions(simulation_cell), electrons(simulation_cell);

  electrons.setName("e");
  electrons.create({1,1});
  SpeciesSet& especies  = electrons.getSpeciesSet();
  int elec_up  = especies.addSpecies("e_u");
  int elec_down  = especies.addSpecies("e_d");
  electrons.R[0][0] = 0.2;
  electrons.R[0][1] = 0.2;
  electrons.R[0][2] = 0.2;
  electrons.R[1][0] = 0.1;
  electrons.R[1][1] = 0.1;
  electrons.R[1][2] = 0.1;

  ions.create({1});
  ions.setName("ions");
  SpeciesSet& tspecies  = ions.getSpeciesSet();
  int ion_a  = tspecies.addSpecies("H");
  ions.R[0][0] = 0.0;
  ions.R[0][1] = 0.0;
  ions.R[0][2] = 0.0;
  ions.update();

  int ee_table = electrons.addTable(electrons);
  int ei_table = electrons.addTable(ions);
  electrons.update();

 std::vector<double> true_d_dr_bispectrum_e1 = {0.203528, 0.585179, 1.23322, 1.32148, 1.67693, 0.203528, 0.585179, 1.23322, 1.32148, 1.67693, 0.203528, 0.585179, 1.23322, 1.32148, 1.67693, 0.0679281, 0.314384, 0.827833, 0.916299, 1.27256, 0.0679281, 0.314384, 0.827833, 0.916299, 1.27256, 0.0679281, 0.314384, 0.827833, 0.916299, 1.27256, 0.135651, 0.449426, 1.02959, 1.11785, 1.4733, 0.135651, 0.449426, 1.02959, 1.11785, 1.4733, 0.135651, 0.449426, 1.02959, 1.11785, 1.4733};
 std::vector<double> true_d_dr_bispectrum_e2= {-0.0678768, -0.135754, -0.20363, -0.20363, -0.203631, -0.0678768, -0.135754, -0.20363, -0.20363, -0.203631, -0.0678768, -0.135754, -0.20363, -0.20363, -0.203631, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0678768, 0.135754, 0.20363, 0.20363, 0.203631, 0.0678768, 0.135754, 0.20363, 0.20363, 0.203631, 0.0678768, 0.135754, 0.20363, 0.20363, 0.203631};
 std::vector<double> true_d_dr_bispectrum_ion = {-0.135651, -0.449426, -1.02959, -1.11785, -1.4733, -0.135651, -0.449426, -1.02959, -1.11785, -1.4733, -0.135651, -0.449426, -1.02959, -1.11785, -1.4733, -0.0679281, -0.314384, -0.827833, -0.916299, -1.27256, -0.0679281, -0.314384, -0.827833, -0.916299, -1.27256, -0.0679281, -0.314384, -0.827833, -0.916299, -1.27256, -0.203528, -0.585179, -1.23322, -1.32148, -1.67693, -0.203528, -0.585179, -1.23322, -1.32148, -1.67693, -0.203528, -0.585179, -1.23322, -1.32148, -1.67693};


  auto jas = std::make_unique<SNAJastrow>(std::string("snap"), ions, electrons,std::string("linear"), 2, 7.0);
  jas->update_sna_rij(electrons,0,false);
  for (int i = 0; i < 2; i ++)
    for (int j = 0; j < 3; j ++)
      jas->sna_desc.rij[i][j]*=-1.0;
  // there is a switch in the sign convention between rij for sna to snad
  double true_dist = -0.1;
  double internal_dist = -electrons.getDistTableAA(ee_table).getDisplRow(1)[0][0];
  REQUIRE(jas->sna_desc.rij[0][0] == Approx(internal_dist));//elec1
  REQUIRE(jas->sna_desc.rij[0][0] == Approx(true_dist));//elec1
  jas->compute_d_dr_bispectrum(0,jas->snad);
  jas->update_sna_rij(electrons,1,false);
  for (int i = 0; i < 2; i ++)
    for (int j = 0; j < 3; j ++)
      jas->sna_desc.rij[i][j]*=-1.0;
  jas->compute_d_dr_bispectrum(1,jas->snad);
  jas->update_sna_rij(electrons,2,false);
  for (int i = 0; i < 2; i ++)
    for (int j = 0; j < 3; j ++)
      jas->sna_desc.rij[i][j]*=-1.0;
  jas->compute_d_dr_bispectrum(2,jas->snad);
  for (int i = 0; i< 40; i++)
  REQUIRE(jas->snad[0][0] == Approx(true_d_dr_bispectrum_e1[0]));//elec1
  REQUIRE(jas->snad[0][39] == Approx(true_d_dr_bispectrum_e1[39]));//elec1
  REQUIRE(jas->snad[1][0] == Approx(true_d_dr_bispectrum_e2[0]));//elec1
  REQUIRE(jas->snad[1][39] == Approx(true_d_dr_bispectrum_e2[39]));//elec1
  internal_dist = -electrons.getDistTableAB(ei_table).getDisplRow(0)[0][0];
  // allow sign to stay the same from elec table since we are reversing it for the gradient calc
  REQUIRE(jas->sna_desc.rij[0][0] == Approx(internal_dist));//elec1
  REQUIRE(jas->snad[2][0] == Approx(true_d_dr_bispectrum_ion[0]));//elec1
  REQUIRE(jas->snad[2][39] == Approx(true_d_dr_bispectrum_ion[39]));//elec1
}


TEST_CASE("snap_jastrow_sna_rij_proposed_update", "[wavefunction]"){
  Communicate* c = OHMMS::Controller;
  std::cout<< "starting snap_jastrow_proposed_uodate" <<std::endl;
  const SimulationCell simulation_cell;
  ParticleSet ions(simulation_cell), electrons(simulation_cell);

  electrons.setName("e");
  electrons.create({1,1});
  SpeciesSet& especies  = electrons.getSpeciesSet();
  int elec_up  = especies.addSpecies("e_u");
  int elec_down  = especies.addSpecies("e_d");
  electrons.R[0][0] = 0.2;
  electrons.R[0][1] = 0.2;
  electrons.R[0][2] = 0.2;
  electrons.R[1][0] = 0.1;
  electrons.R[1][1] = 0.1;
  electrons.R[1][2] = 0.1;

  ions.create({1});
  ions.setName("ions");
  SpeciesSet& tspecies  = ions.getSpeciesSet();
  int ion_a  = tspecies.addSpecies("H");
  ions.R[0][0] = 0.0;
  ions.R[0][1] = 0.0;
  ions.R[0][2] = 0.0;
  ions.update();

  int ee_table = electrons.addTable(electrons);
  int ei_table = electrons.addTable(ions);
  electrons.update();
  std::cout <<electrons.getDistTableAA(ee_table).getDisplRow(1)[0] <<std::endl;

 std::vector<double> true_bispectrum_e1 = {26.5718, 52.7622, 78.3867, 78.2021, 77.4454};
 std::vector<double> true_bispectrum_e2 = {26.8206, 53.2578, 79.126, 78.9405, 78.1798};
 std::vector<double> true_bispectrum_ion = {26.7298, 53.0779, 78.8595, 78.6748, 77.9175};


  auto jas = std::make_unique<SNAJastrow>(std::string("snap"), ions, electrons,std::string("linear"), 2, 7.0);
  
  ParticleSet::SingleParticlePos disp(0.2, 0.3, 0.4);

  electrons.makeMove(0, disp, false);
  std::cout <<electrons.R <<std::endl;
  std::cout <<electrons.activeR(0) <<std::endl <<std::endl;
  std::cout << electrons.getDistTableAA(ee_table).getDisplRow(1)[0][0] << std::endl;
  std::cout << electrons.getDistTableAA(ee_table).getDisplRow(0)[1][0] << std::endl;
  std::cout <<electrons.getDistTableAA(ee_table).getTempDispls()[0] <<std::endl;
  std::cout <<electrons.getDistTableAA(ee_table).getTempDispls()[1] <<std::endl;
  std::cout <<electrons.getDistTableAB(ei_table).getTempDispls()[0] <<std::endl;

  jas->update_sna_rij(electrons,0,false);
  double true_dist = 0.1;
  double internal_dist = electrons.getDistTableAA(ee_table).getDisplRow(1)[0][0];
  REQUIRE(jas->sna_desc.rij[0][0] == Approx(true_dist));//elec1
  REQUIRE(jas->sna_desc.rij[0][0] == Approx(internal_dist));//elec1

  
  jas->update_sna_rij(electrons,0,true);
  true_dist = 0.3;
  internal_dist = -1.0*electrons.getDistTableAA(ee_table).getTempDispls()[1][0];
  REQUIRE(jas->sna_desc.rij[0][0] == Approx(true_dist));//elec1
  REQUIRE(jas->sna_desc.rij[0][0] == Approx(internal_dist));//elec1
  jas->compute_bispectrum(0);
  REQUIRE(jas->sna[0][0] == Approx(true_bispectrum_e1[0]));//elec1
  
  jas->update_sna_rij(electrons,1,true);
  true_dist = -0.3;
  internal_dist = electrons.getDistTableAA(ee_table).getTempDispls()[1][0];
  REQUIRE(jas->sna_desc.rij[0][0] == Approx(true_dist));//elec1
  REQUIRE(jas->sna_desc.rij[0][0] == Approx(internal_dist));//elec1
  jas->compute_bispectrum(1);
  REQUIRE(jas->sna[1][0] == Approx(true_bispectrum_e2[0]));//elec1
  REQUIRE(jas->sna[1][1] == Approx(true_bispectrum_e2[1]));//elec1
  REQUIRE(jas->sna[1][4] == Approx(true_bispectrum_e2[4]));//elec1
  
  jas->update_sna_rij(electrons,2,false);
  true_dist = -0.2;
  internal_dist = electrons.getDistTableAB(ei_table).getDisplRow(0)[0][0];
  REQUIRE(jas->sna_desc.rij[0][0] == Approx(true_dist));//elec1
  REQUIRE(jas->sna_desc.rij[0][0] == Approx(internal_dist));//elec1

  jas->update_sna_rij(electrons,2,true);
  true_dist = -0.4;
  internal_dist = electrons.getDistTableAB(ei_table).getTempDispls()[0][0];
  REQUIRE(jas->sna_desc.rij[0][0] == Approx(true_dist));//elec1
  REQUIRE(jas->sna_desc.rij[0][0] == Approx(internal_dist));//elec1
  jas->compute_bispectrum(2);
  REQUIRE(jas->sna[2][0] == Approx(true_bispectrum_ion[0]));//elec1
} 

TEST_CASE("snap_jastrow_checkinvariables ", "[wavefunction]"){
  Communicate* c = OHMMS::Controller;
  std::cout<< "starting test" <<std::endl;
  const SimulationCell simulation_cell;
  ParticleSet ions(simulation_cell), electrons(simulation_cell);

  electrons.create({1,1});
  electrons.setName("e");
  SpeciesSet& especies  = electrons.getSpeciesSet();
  int elec_up  = especies.addSpecies("e_u");
  int elec_down  = especies.addSpecies("e_d");
  electrons.R[0][0] = 0.2;
  electrons.R[0][1] = 0.2;
  electrons.R[0][2] = 0.2;
  electrons.R[1][0] = 0.1;
  electrons.R[1][1] = 0.1;
  electrons.R[1][2] = 0.1;

  ions.create({1});
  ions.setName("ions");
  SpeciesSet& tspecies  = ions.getSpeciesSet();
  int ion_a  = tspecies.addSpecies("H");
  ions.R[0][0] = 0.0;
  ions.R[0][1] = 0.0;
  ions.R[0][2] = 0.0;

  ions.update();
  electrons.update();


  const char * xmltext = R"XML(<tmp>
  <wavefunction name="psi0" target="e">
  <jastrow name="snap" type="snap" function="snap">
  </jastrow>
</wavefunction>
</tmp>)XML";

  Libxml2Document doc;
  bool okay;
  okay    = doc.parseFromString(xmltext);
  REQUIRE(okay);
  std::cout << "initialized system" << std::endl;
  xmlNodePtr root = doc.getRoot();
  xmlNodePtr jas_node = xmlFirstElementChild(root);
  xmlNodePtr corr_node = xmlFirstElementChild(jas_node);
  SNAJastrowBuilder SJBuilder(c,electrons,ions);
  auto sj_uptr = SJBuilder.buildComponent(corr_node);
  SNAJastrow* sj = static_cast<SNAJastrow*>(sj_uptr.get());

  opt_variables_type a;
  sj->checkInVariablesExclusive(a);

  a[0] = 0.5;
  a[1] = 0.3;
  a[2] = 0.02;
  a[3] = 0.256;
  a[4] = 0.311;
  a[5] = 0.034;
  a[6] = 0.37;
  a[7] = 0.41;
  a[8] = 0.21; // don't change the last one indexed  of 9 to ensure it is unchanged.
  // then call checkinparameters
  sj->resetParametersExclusive(a);
  
  REQUIRE(sj->snap_beta[0][0]==0.5 );
  REQUIRE(sj->snap_beta[0][1]==0.3 );
  REQUIRE(sj->snap_beta[0][2]==0.02 );
  REQUIRE(sj->snap_beta[0][3]==0.256 );
  REQUIRE(sj->snap_beta[0][4]==0.311 );
  REQUIRE(sj->snap_beta[1][0]==0.034 );
  REQUIRE(sj->snap_beta[1][1]==0.37 );
  REQUIRE(sj->snap_beta[1][2]==0.41 );
  REQUIRE(sj->snap_beta[1][3]==0.21 );
  REQUIRE(sj->snap_beta[1][4]==0.0 );

};

TEST_CASE("snap_jastrow_multiple_of_one_particle_type", "[wavefunction]"){
  Communicate* c = OHMMS::Controller;
  std::cout<< "starting test" <<std::endl;
  const SimulationCell simulation_cell;
  ParticleSet ions(simulation_cell), electrons(simulation_cell);
  //TODO test: in a case where you have two particles of one type, ensure youre accessing the correct type block in the lammps snap global array
  electrons.create({2,1});
  electrons.setName("e");
  SpeciesSet& especies  = electrons.getSpeciesSet();
  int elec_up  = especies.addSpecies("e_u");
  int elec_down  = especies.addSpecies("e_d");
  electrons.R[0][0] = 0.2;
  electrons.R[0][1] = 0.2;
  electrons.R[0][2] = 0.2;
  electrons.R[1][0] = 0.1;
  electrons.R[1][1] = 0.1;
  electrons.R[1][2] = 0.1;
  electrons.R[2][0] = 0.3;
  electrons.R[2][1] = 0.3;
  electrons.R[2][2] = 0.3;

  ions.create({1});
  ions.setName("ions");
  SpeciesSet& tspecies  = ions.getSpeciesSet();
  int ion_a  = tspecies.addSpecies("H");

  ions.R[0][0] = 0.0;
  ions.R[0][1] = 0.0;
  ions.R[0][2] = 0.0;
  
  ions.update();
  electrons.update();
};


TEST_CASE("snap_jastrow_linear_vs_quad_form", "[wavefunction]"){
  Communicate* c = OHMMS::Controller;
  std::cout<< "starting test to insight linear or quadratic" <<std::endl;
  const SimulationCell simulation_cell;
  ParticleSet ions(simulation_cell), electrons(simulation_cell);

  electrons.create({1,1});
  electrons.setName("e");
  SpeciesSet& especies  = electrons.getSpeciesSet();
  int elec_up  = especies.addSpecies("e_u");
  int elec_down  = especies.addSpecies("e_d");
  electrons.R[0][0] = 0.2;
  electrons.R[0][1] = 0.2;
  electrons.R[0][2] = 0.2;
  electrons.R[1][0] = 0.1;
  electrons.R[1][1] = 0.1;
  electrons.R[1][2] = 0.1;

  ions.create({1});
  ions.setName("ions");
  SpeciesSet& tspecies  = ions.getSpeciesSet();
  int ion_a  = tspecies.addSpecies("H");

  ions.R[0][0] = 0.0;
  ions.R[0][1] = 0.0;
  ions.R[0][2] = 0.0;
  ions.update();
  electrons.update();


  const char * xmltext = R"XML(<tmp>
  <wavefunction name="psi0" target="e">
  <jastrow name="snap" type="snap" function="snap" snap_type="linear">
  </jastrow>
</wavefunction>
</tmp>)XML";

  Libxml2Document doc;
  bool okay;
  okay    = doc.parseFromString(xmltext);
  REQUIRE(okay);
  xmlNodePtr root = doc.getRoot();
  xmlNodePtr jas_node = xmlFirstElementChild(root);
  xmlNodePtr corr_node = xmlFirstElementChild(jas_node);
  SNAJastrowBuilder SJBuilder(c,electrons,ions);
  auto sj_uptr = SJBuilder.buildComponent(corr_node);
  SNAJastrow* sj = static_cast<SNAJastrow*>(sj_uptr.get());
};

TEST_CASE("snap_jastrow_ion_nelec", "[wavefunction]"){
  Communicate* c = OHMMS::Controller;
  std::cout<< "starting test more than one elec of certain type" <<std::endl;
  const SimulationCell simulation_cell;
  ParticleSet ions(simulation_cell), electrons(simulation_cell);

  electrons.create({2,2});
  electrons.setName("e");
  SpeciesSet& especies  = electrons.getSpeciesSet();
  int elec_up  = especies.addSpecies("e_u");
  int elec_down  = especies.addSpecies("e_d");
  electrons.R[0][0] = 0.2;
  electrons.R[0][1] = 0.2;
  electrons.R[0][2] = 0.2;
  electrons.R[1][0] = 0.1;
  electrons.R[1][1] = 0.1;
  electrons.R[1][2] = 0.1;
  electrons.R[2][0] = 0.3;
  electrons.R[2][1] = 0.3;
  electrons.R[2][2] = 0.3;

  ions.create({1});
  ions.setName("ions");
  SpeciesSet& tspecies  = ions.getSpeciesSet();
  int ion_a  = tspecies.addSpecies("H");

  ions.R[0][0] = 0.0;
  ions.R[0][1] = 0.0;
  ions.R[0][2] = 0.0;
  ions.update();
  electrons.update();


  const char * xmltext = R"XML(<tmp>
  <wavefunction name="psi0" target="e">
  <jastrow name="snap" type="snap" function="snap" snap_type="linear" rcut="7">
  </jastrow>
</wavefunction>
</tmp>)XML";

  Libxml2Document doc;
  bool okay;
  okay    = doc.parseFromString(xmltext);
  REQUIRE(okay);
  xmlNodePtr root = doc.getRoot();
  xmlNodePtr jas_node = xmlFirstElementChild(root);
  xmlNodePtr corr_node = xmlFirstElementChild(jas_node);
  SNAJastrowBuilder SJBuilder(c,electrons,ions);
  auto sj_uptr = SJBuilder.buildComponent(corr_node);
  SNAJastrow* sj = static_cast<SNAJastrow*>(sj_uptr.get());
  REQUIRE(sj->twojmax == 2);
  REQUIRE(sj->ncoeff == 5);
  REQUIRE(sj->snap_beta[1][4]==0.0 );
  // check that size of coefficient matrix is 4, one row for each unique particle type (e_u,e_d,C,H)
  REQUIRE(sj->snap_beta.size() == 3);
  // check that the first row has right number of columns
  REQUIRE(sj->snap_beta[0].size() == 5);

};

TEST_CASE("snap_jastrow_molecule", "[wavefunction]"){
  Communicate* c = OHMMS::Controller;
  std::cout<< "starting test for molecule" <<std::endl;
  const SimulationCell simulation_cell;
  ParticleSet ions(simulation_cell), electrons(simulation_cell);

  electrons.create({1,1});
  electrons.setName("e");
  SpeciesSet& especies  = electrons.getSpeciesSet();
  int elec_up  = especies.addSpecies("e_u");
  int elec_down  = especies.addSpecies("e_d");
  electrons.R[0][0] = 0.2;
  electrons.R[0][1] = 0.2;
  electrons.R[0][2] = 0.2;
  electrons.R[1][0] = 0.1;
  electrons.R[1][1] = 0.1;
  electrons.R[1][2] = 0.1;

  ions.create({1,1});
  ions.setName("ions");
  SpeciesSet& tspecies  = ions.getSpeciesSet();
  int ion_a  = tspecies.addSpecies("H");
  int ion_b  = tspecies.addSpecies("C");
  ions.R[0][0] = 0.0;
  ions.R[0][1] = 0.0;
  ions.R[0][2] = 0.0;
  ions.R[1][0] = 0.15;
  ions.R[1][1] = 0.15;
  ions.R[1][2] = 0.15;
  ions.update();
  electrons.update();

  const char * xmltext = R"XML(<tmp>
  <wavefunction name="psi0" target="e">
  <jastrow name="snap" type="snap" function="snap" snap_type="linear" rcut="7">
  </jastrow>
</wavefunction>
</tmp>)XML";

  Libxml2Document doc;
  bool okay;
  okay    = doc.parseFromString(xmltext);
  REQUIRE(okay);
  xmlNodePtr root = doc.getRoot();
  xmlNodePtr jas_node = xmlFirstElementChild(root);
  xmlNodePtr corr_node = xmlFirstElementChild(jas_node);
  SNAJastrowBuilder SJBuilder(c,electrons,ions);
  auto sj_uptr = SJBuilder.buildComponent(corr_node);
  SNAJastrow* sj = static_cast<SNAJastrow*>(sj_uptr.get());
  REQUIRE(sj->twojmax == 2);
  REQUIRE(sj->ncoeff == 5);
  // check that size of coefficient matrix is 4, one row for each unique particle type (e_u,e_d,C,H)
  REQUIRE(sj->snap_beta.size() == 4);
  // check that the first row has right number of columns
  REQUIRE(sj->snap_beta[0].size() == 5);
};

// TEST_CASE("snap_jastrow_molecule", "[wavefunction]"){
//   Communicate* c = OHMMS::Controller;
//   std::cout<< "starting test for molecule" <<std::endl;
//   const SimulationCell simulation_cell;
//   ParticleSet ions(simulation_cell), electrons(simulation_cell);
//   electrons.create({1,1});
//   electrons.setName("e_u");
//   electrons.R[0][0] = 0.2;
//   electrons.R[0][1] = 0.2;
//   electrons.R[0][2] = 0.2;
//   electrons.R[1][0] = 0.1;
//   electrons.R[1][1] = 0.1;
//   electrons.R[1][2] = 0.1;

//   ions.create({1});
//   ions.setName("ions");
//   SpeciesSet& tspecies  = ions.getSpeciesSet();
//   int ion_a  = tspecies.addSpecies("H");

//   ions.R[0][0] = 0.0;
//   ions.R[0][1] = 0.0;
//   ions.R[0][2] = 0.0;
//   ions.update();
//   electrons.update();
  
//   const char * xmltext = R"XML(<tmp>
//   <wavefunction name="psi0" target="e">
//   <jastrow name="snap" type="snap" function="snap" snap_type="linear" rcut="7">
//   </jastrow>
// </wavefunction>
// </tmp>)XML";

//   Libxml2Document doc;
//   bool okay;
//   okay    = doc.parseFromString(xmltext);
//   REQUIRE(okay);
//   xmlNodePtr root = doc.getRoot();
//   xmlNodePtr jas_node = xmlFirstElementChild(root);
//   xmlNodePtr corr_node = xmlFirstElementChild(jas_node);
//   SNAJastrowBuilder SJBuilder(c,electrons,ions);
//   auto sj_uptr = SJBuilder.buildComponent(corr_node);
//   SNAJastrow* sj = static_cast<SNAJastrow*>(sj_uptr.get());

//   double initial_lmp_x = 0.2*0.529177;
//   double initial_lmp_y = 0.2*0.529177;
//   double initial_lmp_z = 0.2*0.529177;
// }
}
