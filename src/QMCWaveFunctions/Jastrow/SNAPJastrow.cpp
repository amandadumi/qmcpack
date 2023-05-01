#include "SNAPJastrow.h"


namespace qmcplusplus
{

SNAPJastrow::SNAPJastrow(const ParticleSet& ions, ParticleSet& els){
    Nelec = els.getTotalNum();
    Nions = ions.getTotalNum();
    NIonGroups = ions.groups() ;
    double ncoeff = 56;
    std::vector<double> default_coeff; // this should be a matrix where we resize to size of Nspecies*ncoeff
}

SNAPJastrow::~SNAPJastrow(){}

void SNAPJastrow::initialize_lammps(const ParticleSet& ions, ParticleSet& els){
    lmp->input->one("units  metal");
    lmp->input->one("atom_style  atomic");
    // TODO: will this be set by a qmc box? probably.
    lmp->input->one("boundary  f f f ");
    lmp->input->one("neighbor 1.9 bin");
    lmp->input->one("neigh_modify every 1 delay 1 check yes ");
    //TODO: what will box region be pased on QMC object. Cell?
    lmp->input->one("region	mybox block -50 50 -50 50 -50 50");
    // create a box that will contain the number of species equal to the number of groups.
    std::string temp_command = std::string("create_box") + std::to_string(3) +  "mybox";
    lmp->input->one(temp_command);
    // add atoms to the groups
    char ptype_snap_id; // the initial coeffs (with label from *.snapcoeff file) to use for each group of particle group
    lmp->input->one("group e_u type 1");
    lmp->input->one("group e_d type 2");
    lmp->input->one("group i type 3");
    lmp->input->one("group elecs type 1 2");
    lmp->input->one("group all type 1 2 3");
    for (int ig = 0; ig < els.groups(); ig++) { // loop over groups
        // label groups in lamps
        for (int iat = els.first(ig); iat < els.last(ig); iat++) { // loop over elements in each group
           // place atom in boxes according to their group.
           temp_command = std::string("create_atoms ") + std::to_string(ig+1) + " single " + std::to_string(els.R[iat][0]) + "  " + std::to_string(els.R[iat][1])  + " " + std::to_string(els.R[iat][2]);  
           lmp->input->one(temp_command);
         }
     }
    
    for (int ig = 0; ig < ions.groups(); ig++) { // loop over groups
      for (int iat = ions.first(ig); iat < ions.last(ig); iat++) { // loop over elements in each group
        temp_command = std::string("create_atoms 3 single ") + std::to_string(ions.R[iat][0]) + "  " + std::to_string(ions.R[iat][1])  + " " + std::to_string(ions.R[iat][2]);  
        lmp->input->one(temp_command);
      }
    }
    lmp->input->one("pair_style snap");
    lmp->input->one("pair_coeff * * coeff.snapcoeff param.snapparam snap e_u e_d i" );
    lmp->input->one("compute snap all snap");
    lmp->input->one("compute snap_elec elecs snap");
    lmp->input->one("compute db elecs snad/atom");
    lmp->input->one("run            0");
}

SNAPJastrow::LogValueType SNAPJastrow::evaluateGL(const ParticleSet& P,
                        ParticleSet::ParticleGradient& G,
                        ParticleSet::ParticleLaplacian& L,
                        bool fromscratch){
    double* x = new double [OHMMS_DIM*P.getTotalNum()];
    RealType delta = 0.1; // TODO: find units
    double G_finite_diff_forward;
    double G_finite_diff_back;
    // compute gradient
        // i.e. pull gradient out from lammps.
    auto db = lmp->modify->get_compute_by_id("db");
    auto sna_global = lmp->modify->get_compute_by_id("snap");
    
   int nelec = P.getTotalNum();
   for (int ig = 0; ig < P.groups(); ig++) {
        for (int iel = P.first(ig); iel < P.last(ig); iel++){ // loop over elements in each group
            // loop over elecs in each group
            for (int dim = 0; dim < OHMMS_DIM; dim++){
                //G[iel][dim] = lammps_object
                G[iel][dim] += static_cast<double>(*db->array[iel*dim]);
                //forward direction
                RealType r0 = P.R[iel][dim];
                RealType rp   = r0 + (delta/2);
                //update lammps position
                x[iel+dim] = rp;
                // recalculate bispectrum stuff
                db -> compute_peratom();
                G_finite_diff_forward = static_cast<double>(*db->array[iel*dim]);
                //gradient
                //backward direction
                r0 = P.R[iel][dim];
                RealType rm   = r0 - (delta/2);
                //update lammps position
                x[iel+dim] = rm;
                db -> compute_peratom();
                G_finite_diff_back = static_cast<double>(*db->array[iel*dim]);
                // recalculate bispectrum stuff
                double finite_diff_lap = (G_finite_diff_forward - G_finite_diff_back)/delta;
                //fill L
                L[iel] -=  finite_diff_lap;
            }
        }
    }
   log_value_ = evaluateLog(P,G,L);
   return log_value_;
}

SNAPJastrow::LogValueType SNAPJastrow::evaluateLog(const ParticleSet& P,
                                  ParticleSet::ParticleGradient& G,
                                  ParticleSet::ParticleLaplacian& L){
                         return log_value_;
                   }
}

