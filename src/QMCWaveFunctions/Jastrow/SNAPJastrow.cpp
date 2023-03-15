#include "SNAPJastrow.h"


namespace qmcplusplus
{
SNAPJastorw::SNAPJastrow(ParticleSet& target) :  targetPtcl(target) {}

~SNAPJastorw::SNAPJastrow() = default;

SNAPJastorw::initialize_lammps(ParticleSet& P){
    lmp->input->one("units  metal");
    lmp->input->one("atom_style  atomic");
    // TODO: will this be set by a qmc box? probably.
    lmp->input->one("boundary  f f f ");
    lmp->input->one("neighbor 1.9 bin");
    lmp->input->one("neigh_modify every 1 delay 1 check yes ");
    //TODO: what will box region be pased on QMC object. Cell?
    lmp->input->one("region	mybox block -50 50 -50 50 -50 50");
    // create a box that will contain the number of species equal to the number of groups.
    std::string temp_command = "create_box" + P.groups() +  "mybox";
    lmp->input->one(temp_command);
    lmp->input->one("group all type 1 2");
    // add atoms to the groups
    for (int ig = 0; ig < P.groups(); ig++) { // loop over groups
        // label groups in lamps
        lmp->input->one("group electrons type 2");
    f or (int iat = P.first(ig); iat < P.last(ig); iat++) { // loop over elements in each group
    // place atom in boxes according to their group.
        temp_command = "create_atoms" + std::string(ig)+1 + "single " + std::to_string(P.R[iat][0]) + "  " +std::to_string(P.R[iat][1])  + " " + std::to_string(P.R[iat][2]);  
         lmp->input->one(temp_command);

     }
    }
            // create atoms in a certain group of particles with certain coordinates.

        }
    }
    lmp->input->one("reset_timestep 0");
    lmp->input->one("group all type 1 2");
    lmp->input->one("include WBe_Wood_PRB2019.snap");
    lmp->input->one("velocity       all create 298.15 4928459 rot yes mom yes dist gaussian");
    lmp->input->one("fix            ensemble all nve temp 0 298.15 100 tchain 1");
    lmp->input->one("timestep       .01");

    lmp->input->one("thermo_style   one");
    lmp->input->one("thermo         50");
    // loop over particle types(ions,alpha_es, beta_es?)

    // take position of first particle, create lammps container, put positions of all particles of that type in container
    
    // use particle set to initialize a lammps object with bispectrum randomly initialized.
}
SNAPJastorw::access_bispectrum(){
    pointer_to_bispectrum_coeff =  *lmp;

}
}

