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
// namespace LAMMPS_NS{
{

TEST_CASE("snap_jastrow_init", "[wavefunction]")
{
  std::cout<< "starting test" <<std::endl;
  std::cout<< "starting test" <<std::endl;
// short input xml check that lammps positions are correct.
  const SimulationCell simulation_cell;
  ParticleSet ions(simulation_cell), electrons(simulation_cell);
  std::cout<< "made particle set" <<std::endl;
  std::cout<< "made particle set" <<std::endl;

  electrons.create({1,1});
  electrons.setName("e_u");
  electrons.R[0][0] = 0.2;
  electrons.R[0][1] = 0.2;
  electrons.R[0][2] = 0.2;
  electrons.R[1][0] = 0.1;
  electrons.R[1][1] = 0.1;
  electrons.R[1][2] = 0.1;

  std::cout<< "made elecs" <<std::endl;
  std::cout<< "made elecs" <<std::endl;
  ions.create({1});
  ions.setName("ions");
  SpeciesSet& tspecies  = ions.getSpeciesSet();
  int ion_a  = tspecies.addSpecies("H");
  ions.R[0][0] = 0.0;
  ions.R[0][1] = 0.0;
  ions.R[0][2] = 0.0;
  ions.update();
  electrons.update();

  std::cout << "initialized system" << std::endl;
  ions.update();
  electrons.update();

  std::cout << "initialized system" << std::endl;

   // initialize SK
  //electrons.createSK();
  //electrons.createSK();

  const char * xmltext = R"(<tmp>
  <jastrow name="snap" type="snap" function="snap" snap_type="quadratic" rcut="7">
  </jastrow>
</tmp>)";
  Libxml2Document doc;
  bool okay;
  okay    = doc.parseFromString(xmltext);
  REQUIRE(okay);
  std::cout << "initialized system" << std::endl;
  std::cout << "initialized system" << std::endl;

  xmlNodePtr root = doc.getRoot();

  xmlNodePtr jas_node = xmlFirstElementChild(root);
  auto jas = std::make_unique<SNAPJastrow>(std::string("snap"),ions,electrons,std::string("quadratic"),2, 7.0);
  jas->put(root); 
  void *pos = jas->lmp->atom->x;
  jas->evaluateLog(electrons, electrons.G, electrons.L);
  REQUIRE(jas->lmp->atom->x[0][2] == Approx(0.2*0.529177));//elec1
  REQUIRE(jas->lmp->atom->x[1][2] == Approx(0.1*0.529177));//elec1
  REQUIRE(jas->lmp->atom->x[2][0] == Approx(0.0));//elec1
  
  std::cout << "checking whether expected bispectrum components are present" << std::endl;
  jas->sna_global->compute_array();
  std::vector<double> true_bispectrum{24.934336, 49.81748, 74.623951, 74.598469, 74.4962, 24.934336, 49.81748, 74.623951, 74.4962};
  for (int i=0; i<jas->ncoeff;i++){
    std::cout<< i << " " << jas->sna_global->array[0][i]<<" " <<  true_bispectrum[i] << std::endl; 
    REQUIRE(jas->sna_global->array[0][i] == Approx(true_bispectrum[i]));
  }
  // Checck the derivative in the x direction.
  REQUIRE(jas->sna_global->array[1][0] == 0.0);
  // check the value of the energy in last collumn. we don't use this, but just ensure same structure.
  REQUIRE(jas->sna_global->array[0][10] == 26);

  std::cout << "checking energy calculation" << std::endl;
  std::vector<std::vector<double>> set_coeffs;
  set_coeffs = std::vector<std::vector<double>>(3, std::vector<double>{1.2,1.3,1.4,1.5,1.6});
  std::cout << "checking whether expected bispectrum components are present" << std::endl;
}


TEST_CASE("snap_jastrow_init_with_coeff", "[wavefunction]")
{
  Communicate* c = OHMMS::Controller;
  std::cout<< "starting test" <<std::endl;
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

   // initialize SK
  //electrons.createSK();

  const char * xmltext = R"XML(<tmp>
  <wavefunction name="psi0" target="e">
  <jastrow name="snap" type="snap" function="snap">
    <correlation>
     <coefficients id="eup" type="Array"> 0.1 0.2 0.3 0.4 0.5 </coefficients>
    </correlation>
    <correlation>
     <coefficients id="edown" type="Array"> 1.1 1.2 1.3 1.4 1.5 </coefficients>
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
  std::cout << "initialized system" << std::endl;

  xmlNodePtr root = doc.getRoot();

  xmlNodePtr jas_node = xmlFirstElementChild(root);
  xmlNodePtr corr_node = xmlFirstElementChild(jas_node);
  SNAPJastrowBuilder SJBuilder(c,electrons,ions);
  auto sj_uptr = SJBuilder.buildComponent(corr_node);
  SNAPJastrow* sj = static_cast<SNAPJastrow*>(sj_uptr.get());
  REQUIRE(sj->snap_beta[0][0] == 0.1);
  REQUIRE(sj->snap_beta[0][4] == 0.5);
  REQUIRE(sj->snap_beta[1][4] == 1.5);
  REQUIRE(sj->snap_beta[2][4] == 2.5); 
}

TEST_CASE("snap_jastrow_set_twojmax", "[wavefunction]")
{
  Communicate* c = OHMMS::Controller;
  std::cout<< "starting test" <<std::endl;
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
}

TEST_CASE("snap_jastrow_checkinvariables ", "[wavefunction]"){
  Communicate* c = OHMMS::Controller;
  std::cout<< "starting test" <<std::endl;
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
  REQUIRE(sj->snap_beta[1][4]==0.001 );

}

TEST_CASE("snap_jastrow_multiple_of_one_particle_type", "[wavefunction]"){
  //TODO test: in a case where you have two particles of one type, ensure youre accessing the correct type block in the lammps snap global array
  int x=0;
};


TEST_CASE("snap_jastrow_linear_vs_quad_form", "[wavefunction]"){
  Communicate* c = OHMMS::Controller;
  std::cout<< "starting test" <<std::endl;
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

};

TEST_CASE("snap_jastrow_ion_nelec", "[wavefunction]"){
  Communicate* c = OHMMS::Controller;
  std::cout<< "starting test" <<std::endl;
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
  <jastrow name="snap" type="snap" function="snap" snap_type="quadratic" rcut="7">
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
  REQUIRE(sj->twojmax == 2);
  REQUIRE(sj->ncoeff == 5);
  REQUIRE(sj->snap_beta[1][4]==0.001 );

};

TEST_CASE("snap_jastrow_molecule", "[wavefunction]"){
  Communicate* c = OHMMS::Controller;
  std::cout<< "starting test" <<std::endl;
  const SimulationCell simulation_cell;
  ParticleSet ions(simulation_cell), electrons(simulation_cell);

  electrons.create({1,1});
  electrons.setName("e");
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
  <jastrow name="snap" type="snap" function="snap" snap_type="quadratic" rcut="7">
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
  REQUIRE(sj->twojmax == 2);
  REQUIRE(sj->ncoeff == 5);
  REQUIRE(sj->snap_beta[3][4]==0.001 ); // see if snap beta has right number of species accessible.
};

}
