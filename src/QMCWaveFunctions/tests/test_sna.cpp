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
#include "Particle/VirtualParticleSet.h"
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
  
 std::vector<double> true_bispectrum_e1 = {26.7969, 53.3811, 79.6481, 79.5436, 79.1201};
 std::vector<double> true_bispectrum_e2 = {26.9185, 53.6234, 80.0096, 79.9047, 79.4792};
 std::vector<double> true_bispectrum_ion = {26.7969, 53.3811, 79.6481, 79.5436, 79.1201};



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

 std::vector<double> true_d_dr_bispectrum_e1 = {0.809208, 2.32074, 4.87318, 5.21177, 6.60384, 0.809208, 2.32074, 4.87318, 5.21177, 6.60384, 0.809208, 2.32074, 4.87318, 5.21177, 6.60384, 0.271097, 1.25094, 3.28124, 3.62295, 5.02778, 0.271097, 1.25094, 3.28124, 3.62295, 5.02778, 0.271097, 1.25094, 3.28124, 3.62295, 5.02778, 0.538928, 1.78018, 4.06234, 4.40093, 5.793, 0.538928, 1.78018, 4.06234, 4.40093, 5.793, 0.538928, 1.78018, 4.06234, 4.40093, 5.793};
 std::vector<double> true_d_dr_bispectrum_e2 = {-0.27028, -0.540561, -0.810843, -0.810843, -0.810845, -0.27028, -0.540561, -0.810843, -0.810843, -0.810845, -0.27028, -0.540561, -0.810843, -0.810843, -0.810845, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.27028, 0.540561, 0.810843, 0.810843, 0.810845, 0.27028, 0.540561, 0.810843, 0.810843, 0.810845, 0.27028, 0.540561, 0.810843, 0.810843, 0.810845, };
 std::vector<double> true_d_dr_bispectrum_ion = {-0.538928, -1.78018, -4.06234, -4.40093, -5.793, -0.538928, -1.78018, -4.06234, -4.40093, -5.793, -0.538928, -1.78018, -4.06234, -4.40093, -5.793, -0.271097, -1.25094, -3.28124, -3.62295, -5.02778, -0.271097, -1.25094, -3.28124, -3.62295, -5.02778, -0.271097, -1.25094, -3.28124, -3.62295, -5.02778, -0.809208, -2.32074, -4.87318, -5.21177, -6.60384, -0.809208, -2.32074, -4.87318, -5.21177, -6.60384, -0.809208, -2.32074, -4.87318, -5.21177, -6.60384};


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
  //for (int i = 0; i< 40; i++)
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
  std::cout<< "starting snap_jastrow_proposed_update" <<std::endl;
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

 std::vector<double> true_bispectrum_e1 = {25.3282, 49.2049, 70.997, 70.3636, 67.5536};
 std::vector<double> true_bispectrum_e2= {26.2914, 51.1004, 73.7802, 73.1333, 70.2634};
 std::vector<double> true_bispectrum_ion = {25.94, 50.4234, 72.8148, 72.1791, 69.3592};


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
  for (int i=0; i< true_bispectrum_e1.size(); i++){
    REQUIRE(jas->sna[0][i] == Approx(true_bispectrum_e1[i]));//elec1
  }

  jas->update_sna_rij(electrons,1,true);
  true_dist = -0.3;
  internal_dist = electrons.getDistTableAA(ee_table).getTempDispls()[1][0];
  REQUIRE(jas->sna_desc.rij[0][0] == Approx(true_dist));//elec1
  REQUIRE(jas->sna_desc.rij[0][0] == Approx(internal_dist));//elec1
  jas->compute_bispectrum(1);
  for (int i=0; i< true_bispectrum_e2.size(); i++){
    REQUIRE(jas->sna[1][i] == Approx(true_bispectrum_e2[i]));//elec1
  }
//  REQUIRE(jas->sna[1][0] == Approx(true_bispectrum_e2[0]));//elec1
//  REQUIRE(jas->sna[1][1] == Approx(true_bispectrum_e2[1]));//elec1
//  REQUIRE(jas->sna[1][4] == Approx(true_bispectrum_e2[4]));//elec1
  
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
  for (int i=0; i< true_bispectrum_ion.size(); i++){
    REQUIRE(jas->sna[2][i] == Approx(true_bispectrum_ion[i]));//elec1
  }
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

TEST_CASE("snap_jastrow_check_bispectrum_cutoff", "[wavefunction]")
{
  Communicate* c = OHMMS::Controller;
  std::cout<< "starting snap_jastrow_check_bispectrum_cutoff" <<std::endl;
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
  electrons.R[1][0] = 7.0;
  electrons.R[1][1] = 7.0;
  electrons.R[1][2] = 7.0;

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
  
 std::vector<double> true_bispectrum_e1 = {7.92785, 15.7613, 23.4543, 23.4082, 23.2205,};
 std::vector<double> true_bispectrum_e2= {1.0, 2.0, 3.0, 3.0, 3.0};
 std::vector<double> true_bispectrum_ion = {7.92785, 15.7613, 23.4543, 23.4082, 23.2205,};



  auto jas = std::make_unique<SNAJastrow>(std::string("snap"), ions, electrons,std::string("linear"), 2, 7);
  jas->update_sna_rij(electrons,0,false);
  double internal_disp_x = -electrons.getDistTableAB(ei_table).getDisplRow(0)[0][0];
  double internal_dist_0 = electrons.getDistTableAB(ei_table).getDistRow(0)[0];
  double internal_dist_1 = electrons.getDistTableAB(ei_table).getDistRow(1)[0];
  std::cout<< "internal_disp_x" << internal_disp_x<<std::endl;
  std::cout<< "internal_dist_0" << internal_dist_0<<std::endl;
  std::cout<< "internal_dist_1" << internal_dist_1<<std::endl;
  // check that the first rij entry in sna_desc.rij is the first ion not the second electron
  REQUIRE(jas->sna_desc.rij[0][0] == Approx(internal_disp_x));//elec1
  REQUIRE(jas->sna_desc.rij[0][1] == Approx(internal_disp_x));//elec1
  REQUIRE(jas->sna_desc.rij[0][2] == Approx(internal_disp_x));//elec1
  REQUIRE(jas->num_neigh == 1);//elec1
  std::cout<< "sna_desc.rij size" << jas->sna_desc.rij.size() <<std::endl;
  std::cout<< "sna_desc.rij for non exsistent neighbor" << jas->sna_desc.rij[1][0] <<std::endl;
  jas->compute_bispectrum(0);
  REQUIRE(jas->sna[0][0] == Approx(true_bispectrum_e1[0]));//elec1
  jas->update_sna_rij(electrons,1,false);
  jas->compute_bispectrum(1);
  std::cout<< "sna 1 0 " << jas->sna[1][0] <<std::endl;
  REQUIRE(jas->sna[1][1] == Approx(true_bispectrum_e2[1]));//elec1
  jas->update_sna_rij(electrons,2,false);
  jas->compute_bispectrum(2);
  std::cout<< "sna 2 0 " << jas->sna[2][0] <<std::endl;
  REQUIRE(jas->sna[2][2] == Approx(true_bispectrum_ion[2]));//elec1
  //
  //
 std::vector<double> true_d_dr_bispectrum_e1 = {0.239283, 0.789765, 1.7996, 1.94776, 2.56271, 0.239283, 0.789765, 1.7996, 1.94776, 2.56271, 0.239283, 0.789765, 1.7996, 1.94776, 2.56271, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.239283, 0.789765, 1.7996, 1.94776, 2.56271, 0.239283, 0.789765, 1.7996, 1.94776, 2.56271, 0.239283, 0.789765, 1.7996, 1.94776, 2.56271};
 std::vector<double> true_d_dr_bispectrum_e2 = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
 std::vector<double> true_d_dr_bispectrum_h1 = {-0.239283, -0.789765, -1.7996, -1.94776, -2.56271, -0.239283, -0.789765, -1.7996, -1.94776, -2.56271, -0.239283, -0.789765, -1.7996, -1.94776, -2.56271, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.239283, -0.789765, -1.7996, -1.94776, -2.56271, -0.239283, -0.789765, -1.7996, -1.94776, -2.56271, -0.239283, -0.789765, -1.7996, -1.94776, -2.56271};

  jas->update_sna_rij(electrons,0,false);
  for (int i = 0; i < 2; i ++)
    for (int j = 0; j < 3; j ++)
      jas->sna_desc.rij[i][j]*=-1.0;
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
  
  REQUIRE(jas->snad[0][0] == Approx(true_d_dr_bispectrum_e1[0]));//elec1
  

}

TEST_CASE("snap_jastrow_ratio_check", "[wavefunction]")
{
  Communicate* c = OHMMS::Controller;
  std::cout<< "starting snap_jastrow_ratio_check" <<std::endl;
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
  electrons.R[1][1] = 0.1 ;
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
  

  //sna coefficients all equal to one.
  std::vector<std::vector<double>> snap_beta = std::vector<std::vector<double>>(3, std::vector<double>(5,.1));
  auto jas = std::make_unique<SNAJastrow>(std::string("snap"), ions, electrons,std::string("linear"), 2, 7);
  // what is initial energy with starting configurations?
  double esnap_init;
  jas->calculate_ESNA(electrons, snap_beta, esnap_init, false); 
  std::cout << "esnap is" << esnap_init << std::endl;
  // manually change one particle position.  
  electrons.R[0][0] = 0.4;
  electrons.R[0][1] = 0.5;
  electrons.R[0][2] = 0.6;
  electrons.update();
  double esnap_after; 
  jas->calculate_ESNA(electrons, snap_beta, esnap_after, false); 
  std::cout << "esnap is manually moving particle"<< esnap_after << std::endl;
  // move the particle back 
  electrons.R[0][0] = 0.2;
  electrons.R[0][1] = 0.2;
  electrons.R[0][2] = 0.2;
  electrons.update();
  //make a move.
  ParticleSet::SingleParticlePos disp(0.2, 0.3, 0.4);
  electrons.makeMove(0, disp, true);
  //ensure the energy evaluated for eold is appropriate
  double esnap_move_not_proposed; 
  jas->calculate_ESNA(electrons, snap_beta, esnap_move_not_proposed, false); 
  std::cout << "esnap is after moved by assessing not proposed positions" << esnap_move_not_proposed << std::endl;
  double esnap_move_proposed; 
  jas->calculate_ESNA(electrons, snap_beta, esnap_move_proposed, true); 
  std::cout << "esnap is after moved by assessing proposed positions" << esnap_move_proposed << std::endl;
  jas->snap_beta = snap_beta;
  SNAJastrow::LogValue esnap_eval_log = jas->evaluateLog(electrons, electrons.G, electrons.L);

  double rat_func = jas->ratio(electrons, 0);
  double rat_fd = std::exp(esnap_after - esnap_init) ;
  REQUIRE(rat_func == Approx(rat_fd));
 

  //reset elecron one position
  //move the second electron
  //compare manually computed 
  electrons.R[0][0] = 0.2;
  electrons.R[0][1] = 0.2;
  electrons.R[0][2] = 0.2;
  electrons.update();
  
  electrons.R[1][0] = 0.3;
  electrons.R[1][1] = 0.4;
  electrons.R[1][2] = 0.5;
  electrons.update();
  
  double esnap_second_particle_after; 
  jas->calculate_ESNA(electrons, snap_beta, esnap_second_particle_after, false); 
  

  electrons.R[1][0] = 0.1;
  electrons.R[1][1] = 0.1;
  electrons.R[1][2] = 0.1;
  electrons.update();
  
  electrons.makeMove(1, disp, true);
  double esnap_second_part_proposed;
  jas->calculate_ESNA(electrons, snap_beta, esnap_second_part_proposed, true); 
  double rat_func_second = jas->ratio(electrons,1);
  double rat_fd_second = std::exp(esnap_second_particle_after - esnap_init);
  REQUIRE(rat_func_second == Approx(rat_fd_second));

}

TEST_CASE("snap_jastrow_update_outside_cutoff", "[wavefunction]")
{
  Communicate* c = OHMMS::Controller;
  std::cout<< "starting snap_jastrow_ratio_check" <<std::endl;
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
  electrons.R[1][1] = 0.1 ;
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
  

  //sna coefficients all equal to one.
  std::vector<std::vector<double>> snap_beta = std::vector<std::vector<double>>(3, std::vector<double>(5,.1));
  auto jas = std::make_unique<SNAJastrow>(std::string("snap"), ions, electrons,std::string("linear"), 2, 7);
  ParticleSet::SingleParticlePos disp(0.2, 0.3, 7.0);


  double esnap_init;
  jas->calculate_ESNA(electrons, snap_beta, esnap_init, false); 
  std::cout << "esnap is" << esnap_init << std::endl;
  // manually change one particle position.  
  electrons.R[1][0] += disp[0];
  electrons.R[1][1] += disp[1];
  electrons.R[1][2] += disp[2];
  electrons.update();
  double esnap_after; 
  jas->calculate_ESNA(electrons, snap_beta, esnap_after, false); 
  std::cout << "esnap is manually moving particle"<< esnap_after << std::endl;
  // move the particle back 
  electrons.R[1][0] = 0.1;
  electrons.R[1][1] = 0.1;
  electrons.R[1][2] = 0.1;
  electrons.update();
  //make a move.
  electrons.makeMove(1, disp, true);
  //ensure the energy evaluated for eold is appropriate
  double esnap_move_not_proposed; 
  jas->calculate_ESNA(electrons, snap_beta, esnap_move_not_proposed, false); 
  std::cout << "esnap is after moved by assessing not proposed positions" << esnap_move_not_proposed << std::endl;
  double esnap_move_proposed; 
  jas->calculate_ESNA(electrons, snap_beta, esnap_move_proposed, true); 
  std::cout << "esnap is after moved by assessing proposed positions" << esnap_move_proposed << std::endl;
  jas->snap_beta = snap_beta;
  SNAJastrow::LogValue esnap_eval_log = jas->evaluateLog(electrons, electrons.G, electrons.L);

  double rat_func = jas->ratio(electrons, 0);
  double rat_fd = std::exp(esnap_after - esnap_init) ;
  REQUIRE(rat_func == Approx(rat_fd));
 

  
  
}

TEST_CASE("snap_virtual_move", "[wavefunction]")
{
  Communicate* c = OHMMS::Controller;
  std::cout<< "starting snap_jastrow_ratio_check" <<std::endl;
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
  electrons.R[1][1] = 0.1 ;
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
  


  //Pseudo code:
  //1. calculate snap at current positions.
  //2. make a virtual move wth just one point
  //4. call update_rj_vp and manually check distances
  //sna coefficients all equal to .1.
  std::vector<std::vector<double>> snap_beta = std::vector<std::vector<double>>(3, std::vector<double>(5,.1));
  auto jas = std::make_unique<SNAJastrow>(std::string("snap"), ions, electrons,std::string("linear"), 2, 7);
  jas->snap_beta = snap_beta;

  double esnap_init;
  jas->calculate_ESNA(electrons, snap_beta, esnap_init, false); 
  std::cout << "esnap is" << esnap_init << std::endl;
 
  //make a virtual particle set that has only one virtuall particle
  VirtualParticleSet vp(electrons, 1);
  //actual move the virutal partcle from the position of the first electtron
  vp.makeMoves(electrons, 0, {{0.2, 0.3, 0.4}});
  //check that the position after the move is indeed 0.3
  CHECK(Approx(vp.R[0][0]) == 0.4);
  CHECK(Approx(vp.R[0][1]) == 0.5);
  CHECK(Approx(vp.R[0][2]) == 0.6);
  const DistanceTableAB& dt_vp_ion = vp.getDistTableAB(ei_table);
  //the distance between the virtual particle and the ion should be 0.5
  REQUIRE(dt_vp_ion.getDistances().size() == 1);
  CHECK(Approx(dt_vp_ion.getDistances()[0][0]) == 0.87749);
   //check for when particle iis ref particle
  // move particle 0 according to the first ratio (which there is only one in this particlle set.
  // if done correctlly then particle 0 and the second neighbor, indexed at 1, should be the ion and as such  shoulld be 0.707106 distance
  int ratio_index = 0;
  int particle_index = 0;
  jas->update_sna_rij_vp(vp, particle_index, ratio_index);

  
  double internal_dist = -vp.getDistTableAB(ei_table).getDisplRow(ratio_index)[0][0];
  double true_dist= 0.4;
  REQUIRE(jas->sna_desc.rij[1][0] == Approx(true_dist));//elec1
  REQUIRE(jas->sna_desc.rij[1][0] == Approx(internal_dist));//elec1
 
  std::vector<double> true_bispectrum_e1 = {25.3282, 49.2049, 70.997, 70.3636, 67.5536};
  std::vector<double> true_bispectrum_e2= {26.2914, 51.1004, 73.7802, 73.1333, 70.2634};
  std::vector<double> true_bispectrum_ion = {25.94, 50.4234, 72.8148, 72.1791, 69.3592};
   //check for when particle isn't ref
  jas->compute_bispectrum(0);
  for (int i=0; i< true_bispectrum_e1.size(); i++){
    REQUIRE(jas->sna[0][i] == Approx(true_bispectrum_e1[i]));//elec1
  }

  // now test for creating rij for a particle that sn't the refPtcl in the set.
  jas->update_sna_rij_vp(vp, 1, ratio_index);
  jas->compute_bispectrum(1);
  for (int i=0; i< true_bispectrum_e2.size(); i++){
    REQUIRE(jas->sna[1][i] == Approx(true_bispectrum_e2[i]));//elec1
  }


  // test evaluateRatios
  std::vector<double> ratios(1,0.0);
  std::cout<< ratios[0] <<std::endl;
  jas->evaluateRatios(vp,ratios);
  std::cout<< ratios[0] <<std::endl;
  
  electrons.R[0][0] = 0.4;
  electrons.R[0][1] = 0.5;
  electrons.R[0][2] = 0.6;
  electrons.update();
  //ensure the energy evaluated for eold is appropriate
  double esnap_after;
  jas->calculate_ESNA(electrons, snap_beta, esnap_after, true); 
  double rat_fd = std::exp(esnap_after - esnap_init) ;
  std::cout << rat_fd<< std::endl;
  REQUIRE(ratios[0] == Approx(rat_fd));//elec1
 } 


}
