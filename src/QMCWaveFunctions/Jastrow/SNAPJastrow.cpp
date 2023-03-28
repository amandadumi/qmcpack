#include "SNAPJastrow.h"


namespace qmcplusplus
{
SNAPJastorw::SNAPJastrow(ParticleSet& target) :  targetPtcl(target) {}

~SNAPJastorw::SNAPJastrow() = default;

SNAPJastorw::initialize_lammps(){
    lmp->input->one("units  metal");
    lmp->input->one("atom_style  atomic");
    // TODO: will this be set by a qmc box? probably.
    lmp->input->one("boundary  f f f ");
    lmp->input->one("neighbor 1.9 bin");
    lmp->input->one("neigh_modify every 1 delay 1 check yes ");
    //TODO: what will box region be pased on QMC object. Cell?
    lmp->input->one("region	mybox block -50 50 -50 50 -50 50");
    // create a box that will contain the number of species equal to the number of groups.
    std::string temp_command = "create_box" + 3 +  "mybox";
    lmp->input->one(temp_command);
    // add atoms to the groups
    char ptype_snap_id; // the initial coeffs (with label from *.snapcoeff file) to use for each group of particle group
    lmp->input->one("group e_u type 1");
    lmp->input->one("group e_d type 2");
    lmp->input->one("group i type 3");
    lmp->input->one("group elecs type 1 2");
    lmp->input->one("group all type 1 2 3");
    for (int ig = 0; ig < elecs_.groups(); ig++) { // loop over groups
        // label groups in lamps
        for (int iat = elecs_.first(ig); iat < elecs_.last(ig); iat++) { // loop over elements in each group
           // place atom in boxes according to their group.
           temp_command = "create_atoms " + std::string(ig+1) + " single " + std::to_string(elecs_.R[iat][0]) + "  " +std::to_string(elecs_.R[iat][1])  + " " + std::to_string(elecs_.R[iat][2]);  
           lmp->input->one(temp_command);
         }
     }
    
    for (int iat = ions_.first(ig); iat < ions_.last(ig); iat++) { // loop over elements in each group
        temp_command = "create_atoms 3 single " + std::to_string(ions_.R[iat][0]) + "  " +std::to_string(ions_.R[iat][1])  + " " + std::to_string(ions_.R[iat][2]);  
        lmp->input->one(temp_command);
    
    }
    }
    lmp->input->one("pair_style snap")
    lmp->input->one("pair_coeff * * coeff.snapcoeff param.snapparam snap e_u e_d i" );
    lmp->input->one("compute snap all snap")
    lmp->input->one("compute snap_elec elecs snap")
    lmp->input->one("compute db elecs")
    lmp->input->one("run            0");

    lmp_pos_pntr =; 

}



SNAPJastrow::evaluateGL(ParticleSet& P,
                        ParticleSet::ParticleGradient& G,
                        ParticleSet::ParticleLaplacian& L,){
    RealType delta = 0.1; // TODO: find units
    x = new double[OHMMS_DIM*natoms];
    double G_finite_diff_forward;
    double G_finite_diff_back;
    // compute gradient
        // i.e. pull gradient out from lammps.
    auto db = lmp->modify->get_compute_by_id('db')
    auto sna_global = lmp->modify->get_compute_by_id('snap')
    
   int nelec = P.getTotalNum();
   for (int ig = 0; ig < P.groups(); ig++) {
        for (int iel = P.first(ig); iel < P.last(ig); iel++){ // loop over elements in each group
            // loop over elecs in each group
            for (int dim = 0; dim < OHMMS_DIM; dim++){
                //G[iel][dim] = lammps_object
                G[iat][dim] += db[iel][dim]
                //forward direction
                RealType r0 = P.R[iel][dim];
                RealType rp   = r0 + (delta/2);
                P.R[iel][dim] = rp;
                P.update();
                //update lammps position
                x[iel+dim] = P.R[iel][dim];
                // recalculate bispectrum stuff
                db -> compute_per_atom()
                G_finite_diff_forwared = db[iel][dim]
                //gradient
                //backward direction
                RealType r0 = P.R[iel][dim];
                RealType rm   = r0 - (delta/2);
                P.R[iel][dim] = rm;
                P.update();
                //update lammps position
                x[iel+dim] = P.R[iel][dim];
                // recalculate bispectrum stuff
                double finit_diff_lap = (G_finite_diff_front - G_finite_diff_back)/delta;
                //fill L
                L[iel][dim] -=  finite_diff_lap;
            }
        }
    }
}

SNAPJastow::evaluatelog(const ParticleSet& P,
                                  ParticleSet::ParticleGradient& G,
                                  ParticleSet::ParticleLaplacian& L,
                                  bool fromscratch){

                    }
                                  }

}

