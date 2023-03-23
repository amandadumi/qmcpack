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
    std::string temp_command = "create_box" + P.groups() +  "mybox";
    lmp->input->one(temp_command);
    lmp->input->one("group all type 1 2");
    // add atoms to the groups
    char ptype_snap_id; // the initial coeffs (with label from *.snapcoeff file) to use for each group of particle group
    for (P in {elecs_,ions_}){
         for (int ig = 0; ig < P.groups(); ig++) { // loop over groups
            // label groups in lamps
            lmp->input->one("group electrons type 2");
            for (int iat = P.first(ig); iat < P.last(ig); iat++) { // loop over elements in each group
                // place atom in boxes according to their group.
                temp_command = "create_atoms" + std::string(ig+1) + "single " + std::to_string(P.R[iat][0]) + "  " +std::to_string(P.R[iat][1])  + " " + std::to_string(P.R[iat][2]);  
                lmp->input->one(temp_command);
            }
        }
    }
    lmp->input->one("pair_style snap")
    lmp->input->one("pair_coeff * * coeff.snapcoeff param.snapparam snap e_u e_d ion" );
    lmp->input->one("run            0");

    int b;
    const char *sc = "dblist";
    pointer_to_bispectrum_deriv = static_cast<LAMMPS_NS::PairSNAP*>(lmp->force->pair)->extract(sc,b); 
}



SNAPJastrow::evaluateGL(ParticleSet& P,
                        ParticleSet::ParticleGradient& G,
                        ParticleSet::ParticleLaplacian& L,){
    RealType delta = 0.1; // TODO: find units
    x = new double[OHMMS_DIM*natoms];
    // arge: lammps object, property name x is position, type 0 int 1 double, number of values per atom, data container of correct length
    // (from library.cpp line 2079)
    // compute gradient
        // i.e. pull gradient out from lammps.
    db = lmp->modify->get_compute_by_id('db')
    
   int nelec = P.getTotalNum();
   for (int ig = 0; ig < P.groups(); ig++) {
        ids = (0,1,2,3)??? not sure how to get this....
        lammps_gather_atoms_subset(lmp,(char *) "x",1,3,P.groupsize(ig),ids,x);
        for (int iel = P.first(ig); iel < P.last(ig); iel++){ // loop over elements in each group
            // loop over elecs in each group
            for (int dim = 0; dim < OHMMS_DIM; dim++){
                //G[iel][dim] = lammps_object
                G[iati] += db[iel][dim]
                //forward direction
                RealType r0 = P.R[iel][dim];
                RealType rp   = r0 + delta;
                P.R[iel][dim] = rp;
                P.update();
                //update lammps position
                x[iel+dim] = P.R[iel][dim];
                // recalculate bispectrum stuff
                //gradient
                //backward direction
                RealType r0 = P.R[iel][dim];
                RealType rm   = r0 - delta;
                P.R[iel][dim] = rm;
                P.update();
                //update lammps position
                x[iel+dim] = P.R[iel][dim];
                // recalculate bispectrum stuff
                double finit_diff_lap;
                finite_diff_lap
                //fill L
                L[iel][dim] -=  
            }
        }
    }
    lammps_scatter_atoms(lmp,(char *) "x",1,3,x); // this probably needs to happen elsewehere.
}

SNAPJastow::evaluatelog(const ParticleSet& P,
                                  ParticleSet::ParticleGradient& G,
                                  ParticleSet::ParticleLaplacian& L,
                                  bool fromscratch){

                                  }

}

