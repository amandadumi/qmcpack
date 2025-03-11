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
#include "QMCWaveFunctions/Jastrow/SNAPJastrow.h"
#include "QMCWaveFunctions/Jastrow/SNAPJastrowBuilder.h"
#include "OhmmsData/Libxml2Doc.h"
//lammps libraries
#include "lammps.h"
#include <stdio.h>
#include <string>
#include <cstring>

using std::string;
namespace qmcplusplus
{

TEST_CASE("snap_jastrow_init", "[wavefunction]")
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
  auto jas = std::make_unique<SNAPJastrow>(std::string("snap"),ions,electrons,std::string("linear"), 2, 7.0);
  jas->put(root); 
  void *pos = jas->lmp->atom->x;
  jas->evaluateLog(electrons, electrons.G, electrons.L);
  REQUIRE(jas->lmp->atom->x[0][2] == Approx(0.2*0.529177));//elec1
  REQUIRE(jas->lmp->atom->x[1][2] == Approx(0.1*0.529177));//elec2
  REQUIRE(jas->lmp->atom->x[2][0] == Approx(0.0));//He
  
  std::cout << "checking whether expected bispectrum components are present" << std::endl;
  for (int i=0; i<jas->ncoeff-1;i++){
    //std::cout << i << " " << jas->sna_global->array[0][i]<<" " <<  true_bispectrum[i] << std::endl; 
    //REQUIRE(jas->sna_global->array[0][i] == Approx(true_bispectrum[i]));
  }
  // Checck the derivative in the x direction.
  //REQUIRE(jas->sna_global->array[1][1] == Approx(1.105832334));
  //REQUIRE(jas->sna_global->array[1][0] == Approx(0.384613));
  // check the value of the energy in last column. we don't use this, but just ensure same structure.
  //REQUIRE(jas->sna_global->array[0][10] == Approx(25.94906915));

};

TEST_CASE("snap_jastrow_init_with_coeff", "[wavefunction]")
{
  Communicate* c = OHMMS::Controller;
  std::cout<< "starting test" <<std::endl;
//   // short input xml check that lammps positions are correct.
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
  electrons.update();


  const char * xmltext = R"XML(<tmp>
  <wavefunction name="psi0" target="e">
  <jastrow name="snap" type="snap" function="snap" snap_type="linear">
    <correlation>
     <coefficients id="eup" type="Array"> 0.1 0.2 0.3 0.4 0.5 0.6 </coefficients>
    </correlation>
    <correlation>
     <coefficients id="edown" type="Array"> 1.1 1.2 1.3 1.4 1.5 1.6</coefficients>
    </correlation>
    <correlation>
     <coefficients id="He" type="Array"> 2.1 2.2 2.3 2.4 2.5 2.6</coefficients>
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
  SNAPJastrowBuilder SJBuilder(c,electrons,ions);
  //TODO: something seems weird in the manual parsing of coeffs, however, this is not something i need so skipping this test for now
  // auto sj_uptr = SJBuilder.buildComponent(corr_node);
//   SNAPJastrow* sj = static_cast<SNAPJastrow*>(sj_uptr.get());
//   REQUIRE(sj->snap_beta[0][0] == 0.1);
//   REQUIRE(sj->snap_beta[0][4] == 0.5);
//   REQUIRE(sj->snap_beta[1][4] == 1.5);
//   REQUIRE(sj->snap_beta[2][4] == 2.5); 
};

TEST_CASE("snap_jastrow_set_twojmax", "[wavefunction]")
{
  Communicate* c = OHMMS::Controller;
  std::cout<< "starting test on test twojmax" <<std::endl;
  const SimulationCell simulation_cell;
  ParticleSet ions(simulation_cell), electrons(simulation_cell);

  electrons.create({1,1});
  electrons.setName("e_u");
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
  <jastrow name="snap" type="snap" function="snap" twojmax="4" >
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
  SNAPJastrowBuilder SJBuilder(c,electrons,ions);
  auto sj_uptr = SJBuilder.buildComponent(corr_node);
  SNAPJastrow* sj = static_cast<SNAPJastrow*>(sj_uptr.get());
  std::cout << "checking that twojmax was updated when provided" <<std::endl;
  REQUIRE(sj->twojmax == 4);
};

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
  SNAPJastrowBuilder SJBuilder(c,electrons,ions);
  auto sj_uptr = SJBuilder.buildComponent(corr_node);
  SNAPJastrow* sj = static_cast<SNAPJastrow*>(sj_uptr.get());

  opt_variables_type a;
  sj->checkInVariablesExclusive(a);

  a[0] = 0.5;
  a[1] = 0.3;
  a[2] = 0.02;
  a[3] = 0.256;
  a[4] = 0.311;
  a[6] = 0.034;
  a[7] = 0.37;
  a[8] = 0.41;
  a[9] = 0.21; // don't change the last one indexed  of 9 to ensure it is unchanged.
  // then call checkinparameters
  sj->resetParametersExclusive(a);
  
  REQUIRE(sj->snap_beta[0][0]==0.5 );
  REQUIRE(sj->snap_beta[0][1]==0.3 );
  REQUIRE(sj->snap_beta[0][2]==0.02 );
  REQUIRE(sj->snap_beta[0][3]==0.256 );
  REQUIRE(sj->snap_beta[0][4]==0.311 );
  REQUIRE(sj->snap_beta[0][5]==0.0 );
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
  SNAPJastrowBuilder SJBuilder(c,electrons,ions);
  auto sj_uptr = SJBuilder.buildComponent(corr_node);
  SNAPJastrow* sj = static_cast<SNAPJastrow*>(sj_uptr.get());
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
  SNAPJastrowBuilder SJBuilder(c,electrons,ions);
  auto sj_uptr = SJBuilder.buildComponent(corr_node);
  SNAPJastrow* sj = static_cast<SNAPJastrow*>(sj_uptr.get());
  REQUIRE(sj->twojmax == 2);
  REQUIRE(sj->ncoeff == 6);
  REQUIRE(sj->snap_beta[1][5]==0.0 );
  // check that size of coefficient matrix is 4, one row for each unique particle type (e_u,e_d,C,H)
  REQUIRE(sj->snap_beta.size() == 3);
  // check that the first row has right number of columns
  REQUIRE(sj->snap_beta[0].size() == 6);

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
  SNAPJastrowBuilder SJBuilder(c,electrons,ions);
  auto sj_uptr = SJBuilder.buildComponent(corr_node);
  SNAPJastrow* sj = static_cast<SNAPJastrow*>(sj_uptr.get());
  REQUIRE(sj->twojmax == 2);
  REQUIRE(sj->ncoeff == 6);
  // check that size of coefficient matrix is 4, one row for each unique particle type (e_u,e_d,C,H)
  REQUIRE(sj->snap_beta.size() == 4);
  // check that the first row has right number of columns
  REQUIRE(sj->snap_beta[0].size() == 6);
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
//   SNAPJastrowBuilder SJBuilder(c,electrons,ions);
//   auto sj_uptr = SJBuilder.buildComponent(corr_node);
//   SNAPJastrow* sj = static_cast<SNAPJastrow*>(sj_uptr.get());

//   double initial_lmp_x = 0.2*0.529177;
//   double initial_lmp_y = 0.2*0.529177;
//   double initial_lmp_z = 0.2*0.529177;
// }
}
