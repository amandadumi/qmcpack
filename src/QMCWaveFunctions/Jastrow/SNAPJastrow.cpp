#include "SNAPJastrow.h"


namespace qmcplusplus
{

SNAPJastrow::SNAPJastrow(const ParticleSet& ions, ParticleSet& els){
    Nelec = els.getTotalNum();
    Nions = ions.getTotalNum();
    NIonGroups = ions.groups() ;
    std::string iSpecies, eSpecies1, eSpecies2;
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
    // loop over atom types
    // for (int ig = 0; ig < P.groups(); ig++) {
    //     for (int iel = P.first(ig); iel < P.last(ig); iel++){ // loop over elements in each group
    //         double** coeffs = static_cast<double**>lmp->force->pair->coeffelem;
    //     }
        lmp->input->one("run            0");
        log_value_ =  static_cast<double>lmp->force->pair->eng_vdwl;

    // loop over bispectrum compoenents
    // calculate linear energy 
    return log_value_;
}


void SNAPJastrow::evaluateDerivatives(ParticleSet& P, const opt_variables_type& optvars, Vector<ValueType>& dlogpsi, Vector<ValueType>& dhpsioverpsi)
{
    evaluateDerivativesWF(P, active, dlogpsi);
    bool recalculate(false);
    std::vector<bool> rcsingles(myVars.size(), false);
    for (int k = 0; k < myVars.size(); ++k)
    {
      int kk = myVars.where(k);
      if (kk < 0)
        continue;
      if (active.recompute(kk))
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
      if (active.recompute(kk))
        recalculate = true;
      rcsingles[k] = true;
    }
    if (recalculate){
        const size_t NumVars = myVars.size();
        for (int p = 0; p < NumVars; ++p){
            gradLogPsi[p] = 0.0;
            lapLogPsi[p] = 0.0;
        }
    const auto& d_table: const auto & = P.getDistTableAB(myTableID);
    std::vector<TinyVector<RealType,3>> derivs(NumVars);

    constexpr RealType cone(1);
    const size_t ns = d_table.sources();
    const size_t nt = P.getTotalNum();
    
    std::vector<PosType> displ(nt);

    std::fill(first: derivs.begin(), last:derivs.end(),value: 0);
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
    for (nv = 0; nv < NumVars; nv++){
    fd_coeff[nv] = lmp->force->pair->snap->beta[nv] + delta;
    bd_coeff[nv] = lmp->force->pair->snap->beta[nv] - delta;
    calculate_internal_ESNAP_CD(fd_coeff, fd_u);
    calculate_internal_ESNAP_CD(bd_coeff, bd_u);
    dlogpsi = fd_u - bd_u/(2*delta);
    dLogPsi[p] -= dlogpsi;
    }

    }
     
}

/* calculates esnap based on a set of coefficients manually in qmcpack
used to see impact of small change in coefficients on snap energy (needed to calculated d E/d beta)
without having to internally change the lammps object.
*/
void calculate_internal_ESNAP_CD(std::vector<RealType> new_coeff, new_u){
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
  bispectrum_block = P.groups() // number of groups for electrons will be block for ions. 2*nbispectrum = start of ions.
  //TODO: only written for single ion at the moment.
  for (s = 0; s <ns; s++){
    ESNAP_i = 0;
      for (k=0; k < ncoeff; k++){
        RealType bispectrum_val = sna_global[bispectrum_block*n_bispectrum+k];
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

