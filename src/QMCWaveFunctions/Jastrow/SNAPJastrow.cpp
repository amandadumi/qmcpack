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


    
    int n_particle_types = 0; // not sure how to access this yet. // this will be through distance table!
    for (int i =0; i <n_particle_types){
        n_this_particle_type = i.getTotalNum()
    }

    // loop over particle types(ions,alpha_es, beta_es?)

    // take position of first particle, create lammps container, put positions of all particles of that type in container
    P.R;
    // use particle set to initialize a lammps object with bispectrum randomly initialized.
}
SNAPJastorw::access_bispectrum(){
    pointer_to_bispectrum_coeff =  *lmp;

}
}

