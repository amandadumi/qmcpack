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
  /*
  /*
TEST_CASE("simple_file_run", "[wavefunction]")
{
  const char *lmpargv[] {"liblammps", "-log", "none"};
  int lmpargc = sizeof(lmpargv)/sizeof(const char *);
  LAMMPS_NS::LAMMPS *lmp = new LAMMPS_NS::LAMMPS(lmpargc, (char **)lmpargv,MPI_COMM_WORLD);
  lmp->input->file("diffusing_particle.lam");
  delete lmp;
  //does it make sense to have a test just to make sure this doesn't fail? is there a REQUIRE to add?
}

TEST_CASE("lammps_one_command", "[wavefunction]")
TEST_CASE("lammps_one_command", "[wavefunction]")
{

  const char *lmpargv[] {"liblammps", "-log", "none"};
  int lmpargc = sizeof(lmpargv)/sizeof(const char *);
  LAMMPS_NS::LAMMPS *lmp = new LAMMPS_NS::LAMMPS(lmpargc, (char **)lmpargv,MPI_COMM_WORLD);
  lmp->input->one("units  metal");
  delete lmp;

}

TEST_CASE("lammps_access_pair_class", "[wavefunction]")
TEST_CASE("lammps_access_pair_class", "[wavefunction]")
{
  const char *lmpargv[] {"liblammps", "-log", "none"};
  int lmpargc = sizeof(lmpargv)/sizeof(const char *);
  LAMMPS_NS::LAMMPS *lmp = new LAMMPS_NS::LAMMPS(lmpargc, (char **)lmpargv,MPI_COMM_WORLD);
  lmp->input->file("diffusing_particle.lam");
  void *a;
  double* a_val;
  int b;
  const char *sc = "rcutmax";
  //a = static_cast<LAMMPS_NS::PairSNAP*>(lmp->force->pair)->extract(sc,b); 
  //a_val = (double*)(a);   
  //std::cout << "accessing pair_snap object, rcutmax is " << *a_val << std::endl;
  //REQUIRE(*a_val == Approx(4.61586));
  //a = static_cast<LAMMPS_NS::PairSNAP*>(lmp->force->pair)->extract(sc,b); 
  //a_val = (double*)(a);   
  //std::cout << "accessing pair_snap object, rcutmax is " << *a_val << std::endl;
  //REQUIRE(*a_val == Approx(4.61586));

  delete lmp;
}

TEST_CASE("lammps_update_pos", "[wavefunction]")
TEST_CASE("lammps_update_pos", "[wavefunction]")
{
  const char *lmpargv[] {"liblammps", "-log", "none"};
  int lmpargc = sizeof(lmpargv)/sizeof(const char *);
  LAMMPS_NS::LAMMPS *lmp = new LAMMPS_NS::LAMMPS(lmpargc, (char **)lmpargv,MPI_COMM_WORLD);
  lmp->input->file("diffusing_particle.lam");
  int b;
  int natoms = static_cast<int> (lmp->atom->natoms);
  void *pos = lmp->atom->x;
  double **x = static_cast<double **> (pos);
  std::cout <<"natom " << natoms << std::endl;
  std::cout <<"positions" << std::endl;
  // for (int i=0; i<2; i++){
  //     for (int j=0; j<3; j++){
  //       std::cout << "atom " <<i << "position "<< j << " " << (*x)[(3*i)+j] << std::endl;
  //     }
  //   }
    REQUIRE((*x)[0] == Approx(0.0));
    REQUIRE((*x)[5] == Approx(0.5));
    std::cout <<"passed requires" << std::endl;
    std::cout <<"passed requires" << std::endl;
    (*x)[(3*1)+1] = 0.6;
    std::cout <<"update" << std::endl;
    std::cout <<"update" << std::endl;
    void *pos_new = lmp->atom->x;
    double **x_new = static_cast<double **> (pos_new);
    double **x_new = static_cast<double **> (pos_new);
    
    // for (int i=0; i<2; i++){
    //     for (int j=0; j<3; j++){
    //       std::cout << "atom " <<i << "new position but reaccessed "<< j << " " << (*x_new)[(3*i)+j] << std::endl;
    //     }
    // }
  REQUIRE((*x_new)[0] == Approx(0.0));
  REQUIRE((*x_new)[4] == Approx(0.6));
  REQUIRE((*x_new)[0] == Approx(0.0));
  REQUIRE((*x_new)[4] == Approx(0.6));
  delete lmp;
}

/*
/*
TEST_CASE("pass_ions_to_lammps", "[wavefunction]")
{
  const SimulationCell simulation_cell;
  const SimulationCell simulation_cell;
  ParticleSet ions(simulation_cell);

  ions.create({2});
  ions.R[0][0] = 0.0;
  ions.R[0][1] = 0.0;
  ions.R[0][2] = 0.0;
  ions.R[1][0] = 0.0;
  ions.R[1][1] = 0.0;
  ions.R[1][2] = 0.5;
  const char *lmpargv[] {"liblammps", "-log", "none"};
  int lmpargc = sizeof(lmpargv)/sizeof(const char *);
  LAMMPS_NS::LAMMPS *lmp = new LAMMPS_NS::LAMMPS(lmpargc, (char **)lmpargv,MPI_COMM_WORLD);
  lmp->input->one("units  metal");
  lmp->input->one("atom_style  atomic");
  lmp->input->one("boundary  f f f ");
  lmp->input->one("neighbor 1.9 bin");
  lmp->input->one("neigh_modify every 1 delay 1 check yes ");
  lmp->input->one("region	mybox block -50 50 -50 50 -50 50");
  lmp->input->one("create_box 1 mybox");
  std::string var = "create_atoms 1 single " + std::to_string(ions.R[0][0]) + "  " +std::to_string(ions.R[0][1])  + " " + std::to_string(ions.R[0][2]);  
  lmp->input->one(var);
  var = "create_atoms 1 single " + std::to_string(ions.R[1][0]) + "  " +std::to_string(ions.R[1][1])  + " " + std::to_string(ions.R[1][2]);  
  lmp->input->one(var);
  lmp->input->one("mass 1 4.02");
  lmp->input->one("reset_timestep 0");
  lmp->input->one("group ions type 1");
  lmp->input->one("include Mo_Chen_PRM2017.snap");
  lmp->input->one("velocity       all create 0 4928459 rot yes mom yes dist gaussian");
  lmp->input->one("fix            ensemble all nve temp 0 298.15 100 tchain 1");
  lmp->input->one("timestep       .01");

  lmp->input->one("thermo_style   one");
  lmp->input->one("thermo         50");
  lmp->input->one("run            0");
  void *pos = lmp->atom->x;
  double **x = static_cast<double **> (pos);
  REQUIRE((*x)[0] == Approx(0.0));
  REQUIRE((*x)[5] == Approx(0.5));
  // REQUIRE(*a_val == ValueApprox(4.61586));
  double a;
  a = static_cast<double>(lmp->force->pair->eng_vdwl);
  REQUIRE(a == Approx(91.7483));
  delete lmp;
}
*/
*/

/*
/*
TEST_CASE("pass_ions_and_electrons_to_lammps", "[wavefunction]")
{
  const SimulationCell simulation_cell;
  ParticleSet ions(simulation_cell), electrons(simulation_cell);

  electrons.create({2});
  electrons.R[0][0] = 0.1;
  electrons.R[0][1] = 0.1;
  electrons.R[0][2] = 0.1;
  electrons.R[1][0] = 0.1;
  electrons.R[1][1] = 0.1;
  electrons.R[1][2] = 0.5;

    ions.create({2});
  ions.R[0][0] = 0.0;
  ions.R[0][1] = 0.0;
  ions.R[0][2] = 0.0;
  ions.R[1][0] = 0.0;
  ions.R[1][1] = 0.0;
  ions.R[1][2] = 0.5;
  const char *lmpargv[] {"liblammps", "-log", "none"};
  int lmpargc = sizeof(lmpargv)/sizeof(const char *);
  LAMMPS_NS::LAMMPS *lmp = new LAMMPS_NS::LAMMPS(lmpargc, (char **)lmpargv,MPI_COMM_WORLD);
  lmp->input->one("units  metal");
  lmp->input->one("atom_style  atomic");
  lmp->input->one("boundary  f f f ");
  lmp->input->one("neighbor 1.9 bin");
  lmp->input->one("neigh_modify every 1 delay 1 check yes ");
  lmp->input->one("region	mybox block -50 50 -50 50 -50 50");
  lmp->input->one("create_box 2 mybox");
  lmp->input->one("mass 1 4.02");
  lmp->input->one("mass 2 2.01");
  std::string var = "create_atoms 1 single " + std::to_string(ions.R[0][0]) + "  " +std::to_string(ions.R[0][1])  + " " + std::to_string(ions.R[0][2]);  
  lmp->input->one(var);
  var = "create_atoms 1 single " + std::to_string(ions.R[1][0]) + "  " +std::to_string(ions.R[1][1])  + " " + std::to_string(ions.R[1][2]);  
  lmp->input->one(var);
  void *pos = lmp->atom->x;
  double **x = static_cast<double **> (pos);
  REQUIRE((*x)[0] == Approx(0.0));
  REQUIRE((*x)[5] == Approx(0.5));

  var = "create_atoms 2 single " + std::to_string(electrons.R[0][0]) + "  " +std::to_string(electrons.R[0][1])  + " " + std::to_string(electrons.R[0][2]);  
  lmp->input->one(var);
  var = "create_atoms 2 single " + std::to_string(electrons.R[1][0]) + "  " +std::to_string(electrons.R[1][1])  + " " + std::to_string(electrons.R[1][2]);  
  lmp->input->one(var);
  lmp->input->one("reset_timestep 0");
  lmp->input->one("group ions type 1");
  lmp->input->one("group electrons type 2");
  lmp->input->one("group all type 1 2");
  lmp->input->one("include WBe_Wood_PRB2019.snap");
  lmp->input->one("velocity       all create 298.15 4928459 rot yes mom yes dist gaussian");
  lmp->input->one("fix            ensemble all nve temp 0 298.15 100 tchain 1");
  lmp->input->one("timestep       .01");

  lmp->input->one("thermo_style   one");
  lmp->input->one("thermo         50");
  lmp->input->one("run            0");
  REQUIRE((*x)[(2*3)+0] == Approx(0.1)); //access the 3rd particles x-position
  REQUIRE((*x)[(3*3)+2] == Approx(0.5)); //access the 4th particles z-position
  //should also assert group is tracked 
  // double a;
  // a = static_cast<double>(lmp->force->pair->eng_vdwl);
  }
*/
*/

/*
/*
TEST_CASE("get_dbi_drj", "[wavefunction]")
{
  const SimulationCell simulation_cell;
  ParticleSet ions(simulation_cell), electrons(simulation_cell);

  electrons.create({2});
  electrons.R[0][0] = 0.0;
  electrons.R[0][1] = 0.0;
  electrons.R[0][2] = 0.0;
  electrons.R[1][0] = 0.0;
  electrons.R[1][1] = 0.0;
  electrons.R[1][2] = 0.5;


  const char *lmpargv[] {"liblammps", "-log", "none"};
  int lmpargc = sizeof(lmpargv)/sizeof(const char *);
  LAMMPS_NS::LAMMPS *lmp = new LAMMPS_NS::LAMMPS(lmpargc, (char **)lmpargv,MPI_COMM_WORLD);
  lmp->input->file("diffusing_particle.lam");

  void *a;
  double** a_dbdrlist_vals;
  int b;
  const char *sc = "dblist";
  a = static_cast<LAMMPS_NS::PairSNAP*>(lmp->force->pair)->extract(sc,b); 
  a_dbdrlist_vals = (double**)(a);   
  
  std::cout << "made it to parsing object" << std::endl;
  for (int i=0;i<2;i++){
    for (int j = 0 ; j < 2; j++){
      std::cout << (*a_dbdrlist_vals)[i+j] << std::endl;
    }
  }
  delete lmp;
}
*/
}
*/



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
  ions.setName("ions");

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
  <jastrow name="snap" type="snap" function="snap"/>
</tmp>)";
  Libxml2Document doc;
  bool okay;
  okay    = doc.parseFromString(xmltext);
  REQUIRE(okay);
  std::cout << "initialized system" << std::endl;
  std::cout << "initialized system" << std::endl;

  xmlNodePtr root = doc.getRoot();

  xmlNodePtr jas_node = xmlFirstElementChild(root);
  std::cout << "some xml parsing stuff" << std::endl;
  auto jas = std::make_unique<SNAPJastrow>(std::string("snap"),ions,electrons,std::string("linear"),2);
  std::cout << "initial jastrow called finished" << std::endl;
  jas->put(root); 
  std::cout << "some xml parsing stuff" << std::endl;
  std::cout << "some xml parsing stuff" << std::endl;
  void *pos = jas->lmp->atom->x;
  std::cout << "get_position" << std::endl;
  jas->evaluateLog(electrons, electrons.G, electrons.L);
  std::cout << "check positions" << std::endl;
  REQUIRE(jas->lmp->atom->x[0][2] == Approx(0.2*0.529177));//elec1
  REQUIRE(jas->lmp->atom->x[1][2] == Approx(0.1*0.529177));//elec1
  REQUIRE(jas->lmp->atom->x[2][0] == Approx(0.0));//elec1
  std::cout << "check 2 position" << std::endl;
  
  std::cout << "checking whether expected bispectrum components are present" << std::endl;
  jas->sna_global->compute_array();
  std::vector<double> true_bispectrum{24.934336, 49.716745, 74.272267, 74.197307, 73.894457, 24.934336, 49.726569, 74.306526, 74.236353, 73.953025};
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

TEST_CASE("get_dbi_drj", "[particle]")
{
  const SimulationCell simulation_cell;
  ParticleSet ions(simulation_cell), electrons(simulation_cell);

  electrons.create({2});
  electrons.R[0][0] = 0.0;
  electrons.R[0][1] = 0.0;
  electrons.R[0][2] = 0.0;
  electrons.R[1][0] = 0.0;
  electrons.R[1][1] = 0.0;
  electrons.R[1][2] = 0.5;


  const char *lmpargv[] {"liblammps", "-log", "none"};
  int lmpargc = sizeof(lmpargv)/sizeof(const char *);
  LAMMPS_NS::LAMMPS *lmp = new LAMMPS_NS::LAMMPS(lmpargc, (char **)lmpargv,MPI_COMM_WORLD);
  lmp->input->file("diffusing_particle.lam");

  void *a;
  double** a_dbdrlist_vals;
  int b;
  const char *sc = "dblist";
  a = static_cast<LAMMPS_NS::PairSNAP*>(lmp->force->pair)->extract(sc,b); 
  a_dbdrlist_vals = (double**)(a);   
  
  std::cout << "made it to parsing object" << std::endl;
  for (int i=0;i<2;i++){
    for (int j = 0 ; j < 2; j++){
      std::cout << (*a_dbdrlist_vals)[i+j] << std::endl;
    }
  }
  delete lmp;

  }

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
  REQUIRE(sj->snap_beta[1][4]==0.0 );

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


}