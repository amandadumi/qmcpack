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
  std::vector<double> true_bispectrum{25.94907, 51.844575,    77.659852,     77.63319,     77.52618,    25.979614,    51.905603,    77.751273, 77.72458,    77.617449,    25.949069,   0.38461147};

  const char * xmltext = R"(<tmp>
  <jastrow name="snap" type="snap" function="snap" snap_type="linear" rcut="7">
  </jastrow>
</tmp>)";
  Libxml2Document doc;
  bool okay;
  okay    = doc.parseFromString(xmltext);
  REQUIRE(okay);

  xmlNodePtr root = doc.getRoot();

  xmlNodePtr jas_node = xmlFirstElementChild(root);
  auto jas = std::make_unique<SNAJastrow>(std::string("snap"),ions,electrons,std::string("linear"), 2, 7.0);
  std::cout<< "initialiized object" <<std::endl;
  jas->put(root); 
  jas->evaluateLog(electrons, electrons.G, electrons.L);
  jas->update_sna_rij(electrons,0);
  double true_dist = 0.1;
  double internal_dist = electrons.getDistTableAA(ee_table).getDisplRow(1)[0][0];
  REQUIRE(jas->sna_desc.rij[1][0] == Approx(true_dist));//elec1
  REQUIRE(jas->sna_desc.rij[1][0] == Approx(internal_dist));//elec1
  true_dist = -0.2;
  internal_dist = electrons.getDistTableAB(ei_table).getDisplRow(0)[0][0];
  REQUIRE(jas->sna_desc.rij[2][0] == Approx(true_dist));//elec1
  REQUIRE(jas->sna_desc.rij[2][0] == Approx(internal_dist));//elec1
  //update rij to reflect the ion
  jas->update_sna_rij(electrons,2);
  true_dist = 0.2;
  internal_dist = -electrons.getDistTableAB(ei_table).getDisplRow(0)[0][0];
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

  int ee_table = electrons.addTable(electrons);
  int ei_table = electrons.addTable(ions);
  electrons.update();
  
  std::vector<double> true_bispectrum{25.94907, 51.844575,    77.659852,     77.63319,     77.52618,    25.979614,    51.905603,    77.751273, 77.72458,    77.617449,    25.949069,   0.38461147};
  const char * xmltext = R"(<tmp>
  <jastrow name="snap" type="snap" function="snap" snap_type="linear" rcut="7">
  </jastrow>
</tmp>)";



  Libxml2Document doc;
  bool okay;
  okay    = doc.parseFromString(xmltext);
  REQUIRE(okay);
  
  xmlNodePtr root = doc.getRoot();

  xmlNodePtr jas_node = xmlFirstElementChild(root);
  auto jas = std::make_unique<SNAJastrow>(std::string("snap"), ions, electrons,std::string("linear"), 2, 7.0);
  jas->put(root); 
  jas->evaluateLog(electrons, electrons.G, electrons.L);
  jas->update_sna_rij(electrons,0);
  double true_dist = 0.1;
  double internal_dist = electrons.getDistTableAA(ee_table).getDisplRow(1)[0][0];
  REQUIRE(jas->sna_desc.rij[1][0] == Approx(true_dist));//elec1
  jas->compute_bispectrum(0);
  std::cout<< "bisepctrum entry is"<< jas->sna[0][0] <<std::endl;
  std::cout<< "bisepctrum entry is"<< jas->sna[0][1] <<std::endl;
  std::cout<< "bisepctrum entry is"<< jas->sna[0][2] <<std::endl;
  REQUIRE(jas->sna[0][0] == Approx(true_bispectrum[0]));//elec1
  jas->compute_bispectrum(0);
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
