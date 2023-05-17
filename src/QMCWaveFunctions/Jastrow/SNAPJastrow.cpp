#include "SNAPJastrow.h"


namespace qmcplusplus
{

SNAPJastrow::SNAPJastrow(const std::string& obj_name,const ParticleSet& ions, ParticleSet& els) 
  : WaveFunctionComponent(obj_name),
    myTableID(els.addTable(ions)),
    Nelec(els.getTotalNum()),
    Nions(ions.getTotalNum()),
    NIonGroups(ions.groups()),
    Ions(ions)
{
    std::string iSpecies, eSpecies1, eSpecies2;
    ncoeff = 56;
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
    lmp->input->one("compute snap all snap dgradflag 1"); //dgradflag lets us store gradient for each atom intead of derivative!
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
                int ncols = ncoeff+1;
                G[iel][dim] += static_cast<double>(*sna_global->array[(iel*3*ncols)+(dim*ncols)+ncols]);
                //forward direction
                RealType r0 = P.R[iel][dim];
                RealType rp   = r0 + (delta/2);
                //update lammps position
                x[iel+dim] = rp;
                // recalculate bispectrum stuff
                sna_global -> compute_peratom();
                G_finite_diff_forward = static_cast<double>(*sna_global->array[(iel*3*ncols)+(dim*ncols)+ncols]);
                //gradient
                //backward direction
                RealType rm   = r0 - (delta/2);
                //update lammps position
                x[iel+dim] = rm;
                sna_global -> compute_peratom();
                G_finite_diff_back = static_cast<double>(*sna_global->array[(iel*3*ncols)+(dim*ncols)+(ncols)]);
                // recalculate bispectrum stuff
                double finite_diff_lap = (G_finite_diff_forward - G_finite_diff_back)/delta;
                //fill L
                L[iel] -=  finite_diff_lap;
                // return coordinates to original
                x[iel+dim] = r0;

              }
        }
    }
   log_value_ = evaluateLog(P,G,L);
   return log_value_;
}




SNAPJastrow::LogValueType SNAPJastrow::evaluateLog(const ParticleSet& P,
                                  ParticleSet::ParticleGradient& G,
                                  ParticleSet::ParticleLaplacian& L){
    // loop over atom types
    // for (int ig = 0; ig < P.groups(); ig++) {
    //     for (int iel = P.first(ig); iel < P.last(ig); iel++){ // loop over elements in each group
    //         double** coeffs = static_cast<double**>lmp->force->pair->coeffelem;
    //     }
        lmp->input->one("run            0");
        log_value_ =  static_cast<double>(lmp->force->pair->eng_vdwl);

    // loop over bispectrum compoenents
    // calculate linear energy 
    return log_value_;
}


void SNAPJastrow::evaluateDerivatives(ParticleSet& P, const opt_variables_type& optvars, Vector<ValueType>& dlogpsi, Vector<ValueType>& dhpsioverpsi)
{
    evaluateDerivativesWF(P, optvars, dlogpsi);
    bool recalculate(false);
    std::vector<bool> rcsingles(myVars.size(), false);
    for (int k = 0; k < myVars.size(); ++k)
    {
      int kk = myVars.where(k);
      if (kk < 0)
        continue;
      if (optvars.recompute(kk))
        recalculate = true;
      rcsingles[k] = true;
    }
    if (recalculate)
    {
      for (int k = 0; k < myVars.size(); ++k)
      {
        int kk = myVars.where(k);
        if (kk < 0)
          continue;
        if (rcsingles[k])
        {
          dhpsioverpsi[kk] = -RealType(0.5) * ValueType(Sum(lapLogPsi[k])) - ValueType(Dot(P.G, gradLogPsi[k]));
        }
      }
    }
  }



void SNAPJastrow::evaluateDerivativesWF(ParticleSet& P, const opt_variables_type& optvars, Vector<ValueType>& dlogpsi)
{
auto sna_global = lmp->modify->get_compute_by_id("snap");
bool recalculate(false);
std::vector<bool> rcsingles(myVars.size(), false);
for (int k = 0; k < myVars.size(); ++k)
    {
      int kk = myVars.where(k);
      if (kk < 0)
        continue;
      if (optvars.recompute(kk))
        recalculate = true;
      rcsingles[k] = true;
    }
    if (recalculate){
        const size_t NumVars = myVars.size();
        for (int p = 0; p < NumVars; ++p){
            gradLogPsi[p] = 0.0;
            lapLogPsi[p] = 0.0;
        }
    const auto& d_table = P.getDistTableAB(myTableID);

    constexpr RealType cone(1);
    const size_t ns = d_table.sources();
    const size_t nt = P.getTotalNum();
    
    std::vector<PosType> displ(nt);

        /*
        $\partial E_{SNAP}/\partial c$
        1. have snap energy with current coefficients at current positions. 
        2. apply a central difference to the coeffs. FD + BD/2del
        3. make sure lammps has new coeffs.
        4. recalculate lammps energy.
	      how could i do this quickly?
	      just a run 0 on lammps copy? probably just do this myself
        loop thorugh
        */
          //FD
    std::vector<RealType> fd_coeff;
    std::vector<RealType> bd_coeff;
    RealType fd_u;
    RealType bd_u;
    RealType delta = 0.1;
    for (int nv = 0; nv < NumVars; nv++){
      fd_coeff[nv] = lmp->force->pair->snap->beta[nv] + delta;
      bd_coeff[nv] = lmp->force->pair->snap->beta[nv] - delta;
      calculate_internal_ESNAP_CD(fd_coeff, fd_u);
      calculate_internal_ESNAP_CD(bd_coeff, bd_u);
      dlogpsi = fd_u - bd_u/(2*delta);
      dLogPsi[nv] -= dlogpsi;
    }
    calculate_ddc_gradlap_lammps(delta,fd_coeff,bd_coeff);
    }
}


    void SNAPJastrow::calculate_ddc_gradlap_lammps(RealType delta, std::vector<RealType> fd_coeff, std::vector<RealType> bd_coeff){
       /** brute force calculation of grad lammp where lammps object will have coefficients updated and reran.*/ 
       // 1. update coeffs to forward direction
       std::vector<SNAPJastrow::GradDerivVec> grad_forward;
       std::vector<SNAPJastrow::GradDerivVec> grad_back;
       std::vector<SNAPJastrow::GradDerivVec> grad_cd;

       std::vector<SNAPJastrow::ValueDerivVec> lap_forward;
       std::vector<SNAPJastrow::ValueDerivVec> lap_back;
       std::vector<SNAPJastrow::ValueDerivVec> lap_cd;
       std::vector<RealType> lap_back;

       coeff_pointer = lammps->coeffs;
       std::vector<std::vector<RealType,3>> true_coeffs;

       for (i=0;i<fd_coeff;i++){
       true_coeffs[i] = (RealType)*lammps->coeffs[i] // store number not pointer
       // 1. update coeffs to forward direction
       coeff_pointer[i] = fd_coeff[i]
       // 2. run lammps command to update global array
       lmp->input->one("run            0");
       // 3. calculate grad and lap
       int nelec = P.getTotalNum();
       for (int ig = 0; ig < P.groups(); ig++) {
          for (int iel = P.first(ig); iel < P.last(ig); iel++){ 
            for (int dim = 0; dim < OHMMS_DIM; dim++){
             //G[iel][dim] = lammps_object
             int ncols = ncoeff+1
             grad_forward[iel][dim] += static_cast<double>(*sna_global->array[(iel*3*ncols)+(dim*ncols)+ncols]);
             // calc lap d/dr grad
             //forward direction
             RealType r0 = P.R[iel][dim];
             RealType rp   = r0 + (delta/2);
             //update lammps position
             x[iel+dim] = rp;
             // recalculate bispectrum stuff
             sna_global -> compute_peratom();
             G_finite_diff_forward = static_cast<double>(*sna_global->array[(iel*3*ncols)+(dim*ncols)+ncols]);
             //back direction
             RealType rm   = r0 - (delta/2);
             //update lammps position
             x[iel+dim] = rp;
             // recalculate bispectrum stuff
             sna_global -> compute_peratom();
             G_finite_diff_back = static_cast<double>(*sna_global->array[(iel*3*ncols)+(dim*ncols)+ncols]);
             double finite_diff_lap = (G_finite_diff_forward - G_finite_diff_back)/delta;
             //fill L
             lap_forward[iel] -=  finite_diff_lap;
             // return coordinates to original
             x[iel+dim] = r0;
            }
          }
       }
       // 4. update coeffs to backwards direction 
       coeff_pointer[i] = bd_coeff[i]
       // 5. run lammps command to update global array
       lmp->input->one("run            0");
       // 6. calculate grad and lap
       for (int ig = 0; ig < P.groups(); ig++) {
          for (int iel = P.first(ig); iel < P.last(ig); iel++){ 
            for (int dim = 0; dim < OHMMS_DIM; dim++){
             //G[iel][dim] = lammps_object
             int ncols = ncoeff+1
             grad_backward[iel][dim] += static_cast<double>(*sna_global->array[(iel*3*ncols)+(dim*ncols)+ncols]);
             // calc lap d/dr grad
             //forward direction
             RealType r0 = P.R[iel][dim];
             RealType rp   = r0 + (delta/2);
             //update lammps position
             x[iel+dim] = rp;
             // recalculate bispectrum stuff
             sna_global -> compute_peratom();
             G_finite_diff_forward = static_cast<double>(*sna_global->array[(iel*3*ncols)+(dim*ncols)+ncols]);
             //back direction
             RealType rm   = r0 - (delta/2);
             //update lammps position
             x[iel+dim] = rp;
             // recalculate bispectrum stuff
             sna_global -> compute_peratom();
             G_finite_diff_back = static_cast<double>(*sna_global->array[(iel*3*ncols)+(dim*ncols)+ncols]);
             double finite_diff_lap = (G_finite_diff_forward - G_finite_diff_back)/delta;
             //fill L
             lap_back[iel] -=  finite_diff_lap;
             // return coordinates to original
             x[iel+dim] = r0;
             }
          lap_cd = (lap_forward-lap_back)/2*delta;
          }
       }
       // 7. use central difference to get estimate of ddc grad and lap  which populates gradLogPsi and lapLogPsi
       // centra diff grad
       //AD: can probably combine this and above loop. 
       for (int ig = 0; ig < P.groups(); ig++) {
          for (int iel = P.first(ig); iel < P.last(ig); iel++){ 
            for (int dim = 0; dim < OHMMS_DIM; dim++){
              grad_cd[iel][dim] = (grad_forward[iel][dim] - grad_background[iel][dim])/del;
            }
            lap_cd[iel] = (lap_forward[iel]-lap_back[iel])/2*delta;
          }
     
        }
      }
    gradLogPsi = grad_cd;
    lapLogPsi  = lap_cd;
  
  }

/* calculates esnap based on a set of coefficients manually in qmcpack
used to see impact of small change in coefficients on snap energy (needed to calculated d E/d beta)
without having to internally change the lammps object.
*/
void calculate_internal_ESNAP_CD(std::vector<RealType> new_coeff, RealType new_u){
  RealType ESNAP_all = 0;
  RealType ESNAP_elec;
  RealType ESNAP_ion;
  const auto& d_table = P.getDistTableAB(myTableID);

  const size_t ns = d_table.sources();
  const_size nt = P.getTotalNum();

  for (int ig = 0; ig < P.groups(); ig++) {
    for (int iel = P.first(ig); iel < P.last(ig); iel++){ // loop over elements in each group
      ESNAP_i = 0; //  maybe beta_0
      for (k=0; k < ncoeff; k++){
        // determine atom type
        RealType bispectrum_val = sna_global[ig*ncoeff+k]; //block of bispectrum + current component to add.
        ESNAP_i += new_coeff[ig*ncoeff + k] * bispectrum_val  ;
      }
    ESNAP_all += ESNAP_i;
    }
  }
  bispectrum_block = P.groups() // number of groups for electrons will be block for ions. 2*ncoeff = start of ions.
  //TODO: only written for single ion at the moment.
  for (int s = 0; s <ns; s++){
    ESNAP_i = 0;
      for (k=0; k < ncoeff; k++){
        RealType bispectrum_val = sna_global[bispectrum_block*ncoeff+k];
        ESNAP_i += beta[bispectrum_block][k] * bispectrum_val  ;
      }
    ESNAP_all += ESNAP_i;
  }
  new_u = ESNAP;
}

/////// Functions for optimization /////
void SNAPJastrow::checkOutVariables(const opt_variables_type& o ){myVars.getIndex(o);}

/////////////////////////////////// MC Related functions /////////

void SNAPJastrow::acceptMove(ParticleSet& P, int iat, bool safe_to_delay){
    int x = 0;
    // Is this where coefficients are reset or is this just for other variables?  
    // I would think optimizer handles coefficient and updating them. This would be other things that depend on the updated coefficents?

}

void SNAPJastrow::registerData(ParticleSet& P, WFBufferType& buf){
    log_value_ = evaluateLog(P, P.G, P.L);
}
SNAPJastrow::LogValueType SNAPJastrow::updateBuffer(ParticleSet& P, WFBufferType& buf, bool from_scratch){
    log_value_ = evaluateLog(P,P.G,P.L);
    return log_value_;
}

void SNAPJastrow::copyFromBuffer(ParticleSet& P, WFBufferType& buf){
    
}
SNAPJastrow::PsiValueType SNAPJastrow::ratio(ParticleSet& P, int iat){
    PsiValueType ratio;
    ratio = 0;
    return ratio;
}
}

