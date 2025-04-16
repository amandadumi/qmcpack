#include "SNAJastrow.h"
#include "ResourceCollection.h"
#include "CPU/VectorOps.h"


namespace qmcplusplus
{

template<typename T>
struct SNAMultiWalkerMem : public Resource
{
  // fused buffer for fast transfer
  Vector<char, OffloadPinnedAllocator<char>> transfer_buffer;
  // multi walker result
  Vector<T, OffloadPinnedAllocator<T>> mw_vals;
  // multi walker -1
  Vector<int, OffloadPinnedAllocator<int>> mw_minus_one;

  void resize_minus_one(size_t size)
  {
    if (mw_minus_one.size() < size)
    {
      mw_minus_one.resize(size, -1);
      mw_minus_one.updateTo();
    }
  }

  SNAMultiWalkerMem() : Resource("SNAMultiWalkerMem") {}

  SNAMultiWalkerMem(const SNAMultiWalkerMem&) : SNAMultiWalkerMem() {}

  std::unique_ptr<Resource> makeClone() const override { return std::make_unique<SNAMultiWalkerMem>(*this); }
};

SNAJastrow::SNAJastrow(const std::string& obj_name, ParticleSet& ions, ParticleSet& els,const std::string input_snap_type, int input_twojmax, double input_rcut) 
  : WaveFunctionComponent(obj_name),
    OptimizableObject("snap_" + ions.getName()),
    Nions(ions.getTotalNum()),
    Nelec(els.getTotalNum()),
    NIonGroups(ions.groups()),
    ee_Table_ID_(els.addTable(els)),
    ei_Table_ID_(els.addTable(ions)),
    Ions(ions),
    timers_("SNAJastrowTimers")
    {
    
   // create object for 2snap descriptor
    ntypes = NIonGroups+els.groups(); 
    app_debug()<<"entering snap function" <<std::endl;
    sna_desc.init(Nelec+Nions,input_twojmax,input_rcut);
    ii_Table_ID_ = Ions.addTable(Ions); 

    
    twojmax = input_twojmax;
    if (twojmax%2==0){
      int m = (twojmax/2)+1;
      ncoeff = (m*(m+1)*(2*m+1))/6;
    }
    else{
      int m = (twojmax+1)/2;
      ncoeff = (m*(m+1)*(2*m))/3;
    }
    rcut = input_rcut;
    snap_type = input_snap_type;
   // create object for 2snap descriptor
    int wjkkk;
    //TODO: set this to ddefault value for everythin
    rcutij = std::vector<std::vector<double>>(Nions+Nelec,std::vector<double>(Nions+Nelec,7.0));
    radelem = std::vector<double>(ntypes,7.0);
    element = std::vector<int>(Nions+Nelec,0);
    type_map = std::vector<int>(Nions+Nelec,0);
    for (int ig = 0; ig < els.groups(); ig++) { // loop over groups
      for (int iat = els.first(ig); iat < els.last(ig); iat++) { // loop over elements in each group
       // Code that depends on group
        type_map[iat] = ig;
      }
    }
    for (int ig = 0; ig < Ions.groups(); ig++) { // loop over groups
      for (int iat = Ions.first(ig); iat < Ions.last(ig); iat++) { // loop over elements in each group
       // Code that depends on group
        type_map[Nelec+iat] = els.groups()+ig;
      }
    }
    for(int np =0;np<Nions+Nelec;np++)
      app_debug() << np << " " << type_map[np] << std::endl;

    snap_beta = std::vector<std::vector<double>>(NIonGroups+els.groups(), std::vector<double>(ncoeff,0.0));
    sna = std::vector<std::vector<double>>(Nelec+Nions, std::vector<double>(ncoeff,0.0));
    snad = std::vector<std::vector<double>>(Nions + Nelec, std::vector<double>(3*ntypes*ncoeff,0.0));
    for (int k = 0; k < ncoeff; k++){
      for (int i=0; i < NIonGroups+els.groups(); i++){
        std::stringstream name;
        name << "sna_coeff_" << i;
        name << "_"  << k ;
        myVars.insert(name.str(), snap_beta[i][k], true);
        app_debug() << "atom type "<< i<< " coeff " <<k << " "<< snap_beta[i][k] << std::endl;
      } 
    }

    resizeWFOptVectors();
    grad_u.resize(Nelec);
    lap_u.resize(Nelec);
    app_debug() << "at end of init functiion "<<std::endl;

}

SNAJastrow::~SNAJastrow(){

}

void SNAJastrow::set_coefficients(std::vector<double> id_coeffs, int id){
  if (id_coeffs.size() != ncoeff){
    app_log() << "WARNING: number of coeffs less than coeffs/particle_type" <<std::endl; 
    app_log() << "coeffs/particle_type: " << ncoeff << " , coeffs read in: " << id_coeffs.size() << std::endl;
  } 
  app_debug() << "in set coefficients" <<std::endl;
  for (int i=0; i < id_coeffs.size(); i++){
    snap_beta[id][i] = id_coeffs[i];
    app_debug()<< "snap coefficient for particle id: " << id << " is: " << id_coeffs[i]  <<std::endl;;
    myVars[(id*ncoeff) + i] = snap_beta[id][i];
  }

}



 double SNAJastrow::FD_Lap(const ParticleSet& P, int iat, int dim, int coeff, int ntype, std::vector<std::vector<double>> coeffs, bool bispectrum_only){
  double G_finite_diff_forward;
  double G_finite_diff_back;
  // if taking wave function laplacian, we just need bispectrum component but for electron laplacian we use the snap coefficients.
  double this_coeff = 1.0;
  if (not bispectrum_only){
    this_coeff = coeffs[ntype][coeff];
  }
  update_sna_rij(P, iat);
  compute_bispectrum(iat);
  int col = (ntype*(3*ncoeff))+(dim*ncoeff)+coeff;
  RealType r0 = P.R[iat][dim]/bohr_over_ang;
  //forward direction
  RealType rp = r0 + (dist_delta/bohr_over_ang);
  for (int j= 0; j < Nions+Nelec; j++){
       sna_desc.rij[j][dim] += rp;
  }
  compute_bispectrum(iat);
  compute_d_dr_bispectrum(iat);

  G_finite_diff_forward = this_coeff * snad[iat][col] * hartree_over_ev/bohr_over_ang;
  
  //backward direction
  RealType rm  = r0 - (dist_delta/bohr_over_ang);
  for (int j= 0; j < Nions+Nelec;j++){
       sna_desc.rij[j][dim] += rm;
  }
  compute_bispectrum(iat);
  compute_d_dr_bispectrum(iat);
  G_finite_diff_back = this_coeff * snad[iat][col] * hartree_over_ev/bohr_over_ang;
  //fill L
  double finite_diff_lap = (G_finite_diff_forward - G_finite_diff_back)/(2*dist_delta); 
  
  return finite_diff_lap;
 }


  SNAJastrow::LogValue SNAJastrow::evaluateGL(const ParticleSet& P,
                          ParticleSet::ParticleGradient& G,
                          ParticleSet::ParticleLaplacian& L,
                          bool fromscratch){
    app_debug() << "in evaluateGL" <<std::endl;
    return evaluateLog(P,G,L);
  }
 
  void SNAJastrow::computeGL(const ParticleSet& P, int iel){
    app_debug() << "SNAJastrow::evaluateLog entered function" <<std::endl;
    ScopedTimer local_timer(timers_.eval_gl_timer);
    double grad_val;
    int row,col;
    row = iel;             // get the derivative of a specific component with regards to the change in r of this particle. 
    // reset the internal array
    grad_u[iel] = 0;
    lap_u[iel] = 0;
    //loop over types of particles (electrons and ions)
    for (int n=0; n < ntypes; n++){
      // loop over the components
      for (int k = 0; k < ncoeff; k ++){
        //we wil need gradient in each direction.
        for (int dim = 0; dim < 3; dim++){
          // get gradient row which is 3*N rows 
          col = (n*(3*ncoeff))+(dim*ncoeff)+k;
          grad_val = snad[row][col];
          //app_debug() << "grad val is " << grad_val << std::endl;
          // in snap global, the force is stored (i.e. -dB/dr) which is the correct sign since we are taking the negative of Esnap expression
          grad_u[iel][dim] += snap_beta[n][k]*grad_val*hartree_over_ev/bohr_over_ang;
          lap_u[iel] += FD_Lap(P, iel, dim, k, n, snap_beta, false);
        } // end dim loop
      }// end k loop
    } //end n loop
    return;
   }


  SNAJastrow::LogValue SNAJastrow::evaluateLog(const ParticleSet& P,
                                    ParticleSet::ParticleGradient& G,
                                    ParticleSet::ParticleLaplacian& L){
    app_debug() << "SNAJastrow::evaluateLog entered function" <<std::endl;
    ScopedTimer local_timer(timers_.eval_log_timer);
    for (int e=0 ; e< Nelec; e++){
        app_debug() << "SNAJastrow::evaluateLog in particle " << e  <<std::endl;
        update_sna_rij(P, e);  
        compute_bispectrum(e);
        compute_d_dr_bispectrum(e);
    }
    double esnap;
    calculate_ESNA(P, snap_beta, esnap);
        app_debug() << "SNAJastrow::evaluated esnap" <<std::endl;
    log_value_ = static_cast<SNAJastrow::LogValue>(esnap);
    for (int iel = 0; iel < Nelec; iel++) {
      computeGL(P,iel);
      G[iel] += grad_u[iel];
      L[iel] += lap_u[iel];
    }
    return log_value_;
  }

    SNAJastrow::GradType SNAJastrow::evalGrad(ParticleSet& P, int iat){
    app_debug() << "in SNAJastrow::evalGrad" <<std::endl;

    for (int e=0 ; e< Nelec; e++){
      update_sna_rij(P, e);  
      compute_d_dr_bispectrum(e);
    }
    GradType grad_iat;
    for (int dim=0; dim < 3; dim++){
     int row = (3*iat) + dim;
     for (int k = 0; k < ncoeff ; k++){
       for (int n = 0; n < ntypes; n++){
         int col = (n * (3*ncoeff) ) + (dim*ncoeff) + k;
         grad_iat[dim] += snap_beta[n][k]*snad[iat][col];
       }
     }
    }
    return grad_iat;
    }


  void SNAJastrow::evaluateDerivatives(ParticleSet& P, const opt_variables_type& optvars, Vector<ValueType>& dlogpsi, Vector<ValueType>& dhpsioverpsi)
  {
    app_debug() << "in evalderivative"<<std::endl;
    ScopedTimer local_timer(timers_.eval_wf_grad_timer);
    evaluateDerivativesWF(P, optvars, dlogpsi);
    bool recalculate(false);
    std::vector<bool> rcsingles(myVars.size(), false);
    for (int k = 0; k < myVars.size(); ++k)
    {
      int k_global = myVars.where(k);
      if (k_global < 0)
       continue;
      if (optvars.recompute(k_global))
        recalculate = true;
      rcsingles[k] = true;
    }
    if (recalculate)
    {
      for (int k = 0; k < myVars.size(); ++k)
      {
        int k_global = myVars.where(k);
        if (k_global < 0)
          continue;
        if (rcsingles[k])
        {
          dhpsioverpsi[k_global] = -RealType(0.5) * RealType(Sum(lapLogPsi[k]));
          for (int i = 0; i < Nelec; i++){
            dhpsioverpsi[k_global] -= RealType(dot(P.G[i], gradLogPsi[k][i]));
          }
        }
      }
    }
  }


  void SNAJastrow::evaluateDerivativesWF(ParticleSet& P, const opt_variables_type& optvars, Vector<ValueType>& dlogpsi){
    app_debug() << "in evalderivativeWF"<<std::endl;
    for (int e=0 ; e < Nelec; e++){
        update_sna_rij(P,e);  
        compute_bispectrum(e);
        compute_d_dr_bispectrum(e);
    }
    bool recalculate(false);
    resizeWFOptVectors();
    std::vector<bool> rcsingles(myVars.size(), false);
    for (int k = 0; k < myVars.size(); ++k)
    {
      int k_global = myVars.where(k);
      if (k_global < 0)
        continue;
      if (optvars.recompute(k_global))
        recalculate = true;
      rcsingles[k] = true;
    }
    if (recalculate){
      const size_t NumVars = myVars.size();
      for (int p = 0; p < NumVars; p++){
        gradLogPsi[p] = 0.0;
        lapLogPsi[p] = 0.0;
      }
      dLogPsi = 0.0;
    
      for (int k = 0; k < myVars.size(); k++){
        int k_global = myVars.where(k);
        if (k_global < 0)
          continue;
        if (rcsingles[k]){
          evaluate_linear_derivs(P, k);
          dlogpsi[k_global] = ValueType(dLogPsi[k]);
        }
      } 
    }
  }

  void SNAJastrow::evaluate_linear_derivs(ParticleSet& P, int coeff_idx){
    app_debug() << "in linearderivs" <<std::endl;
    
    int ntype = int(coeff_idx/ncoeff);  
    int coeff = coeff_idx%ncoeff; //which coeff of this type are we on.
    if (ntype< P.groups()){
      for (int iat = P.first(ntype); iat < P.last(ntype); iat++) { // loop over elements in each group
        dLogPsi[coeff_idx] += -sna[iat][coeff];
      }
    }
    else{
     for (int iat = Ions.first(ntype-P.groups()); iat < Ions.last(ntype-P.groups()); iat++) { // loop over elements in each group
       dLogPsi[coeff_idx] += -sna[P.getTotalNum()+iat][coeff];
      }
    }
     for (int iel =0; iel < Nelec; iel++){
       for (int dim = 0; dim < OHMMS_DIM; dim++){ // loop over dim to get grad vec.
        int col = (ntype*(3*ncoeff))+(dim*ncoeff)+coeff;
        gradLogPsi[coeff_idx][iel][dim] += snad[iel][col];
        lapLogPsi[coeff_idx][iel] += FD_Lap(P, iel, dim, coeff, ntype, snap_beta, true);
       }
     }
    }


    
  void SNAJastrow::calculate_ESNA(const ParticleSet& P, const std::vector<std::vector<double>> coeff, double& new_u){
    app_debug() << "SNAJastrow::calculate_ESNA entered function" <<std::endl;
    ScopedTimer local_timer(timers_.eval_esnap_timer);
    double esnap_all=0;
    double esnap_elec=0;
    double esnap_ion=0;
    double bispectrum_val;
    // calculate electron contribution
    // the global array is summed over groups of atoms of the same type. thus we just need to sum over groups.
    for (int ig = 0; ig < P.groups(); ig++) {
      esnap_elec += coeff[ig][0]*P.groupsize(ig); // beta0 contribiution
      for (int iat = P.first(ig); iat < P.last(ig); iat++) { // loop over elements in each group
        app_debug() << "elec index is" << iat <<std::endl; 
           update_sna_rij(P, iat); 
           compute_bispectrum(iat);
        for (int k = 0; k < ncoeff; k++){
          bispectrum_val = sna[iat][k]; //block of bispectrum + current component to add.
          esnap_elec += coeff[ig][k] * bispectrum_val;
          //app_debug()<<"snap coefff for group "<<ig<< " coeff "<<k << " is " <<coeff[ig][k] <<std::endl;
        }
      }
    }
    esnap_all += esnap_elec;
    app_debug() << "SNAJastrow::calculate_ESNA calc ion center contribution" <<std::endl;

    for (int ig = 0; ig < Ions.groups(); ig++) {
      esnap_ion += coeff[P.groups()+ig][0]*Ions.groupsize(ig);
      for (int iat = Ions.first(ig); iat < Ions.last(ig); iat++) { // loop over elements in each group
           update_sna_rij(P, Nelec+iat); 
           compute_bispectrum(Nelec+iat);
        app_debug() << "ion index is" << iat <<std::endl; 
        for (int k =0; k < ncoeff; k++){
          bispectrum_val = sna[P.getTotalNum()+iat][k];
          esnap_ion += coeff[P.groups()+ig][k] * bispectrum_val;
          //app_debug()<< "snap coeff for group "<< P.groups()+ig << " coeff " << k << " is " << coeff[P.groups()+ig][k] <<std::endl;
        }
      }
    }
    esnap_all += esnap_ion;

    new_u = -esnap_all;
    return;
  }

  inline void SNAJastrow::update_sna_rij_vp(const VirtualParticleSet& VP, int iat){
 //get original distances at first and then just update info if iat is the virtual particle
  const auto& ee_displ = VP.getRefPS().getDistTableAA(ee_Table_ID_).getDisplRow(iat);
  const auto& ei_displ = VP.getRefPS().getDistTableAB(ei_Table_ID_).getDisplRow(iat);
  if (iat == VP.refPtcl){
        const auto& ee_displ = VP.getDistTableAB(ee_Table_ID_).getDisplRow(iat);
        const auto& ei_displ = VP.getDistTableAB(ei_Table_ID_).getDisplRow(iat);
  }
  const int itype = type_map[iat];
  const double radi = radelem[itype]; 
  bool elec = (itype <2);
  if (elec){ // if iat is elec
    for (int j=0; j < Nelec; j++){ //loop over all other electrons
      if (iat != j ){
         const auto& je_displ = ee_displ[j];//set the distance accordiing to the reference partcle set
         if (j == VP.refPtcl){// if were updatng thhe row for the VP, grab that distance
           const auto& je_displ = VP.getDistTableAB(ei_Table_ID_).getDisplRow(j)[iat];// sign change here?
        }
         int jtype = type_map[Nelec+j];
         int jelem = 0;
         sna_desc.rij[j][0] = je_displ[0];
         sna_desc.rij[j][1] = je_displ[1];
         sna_desc.rij[j][2] = je_displ[2];
         sna_desc.inside[j] = j;
         sna_desc.wj[j] = 1;
         sna_desc.rcutij[j] = (radi + radelem[jtype]) * rcutfac;
         sna_desc.element[j] = jelem;
      }
   }
   for (int j = 0; j< Nions; j++){
       int jtype = type_map[Nelec+j];
       int jelem = 0;
       int j_sna = Nelec+j;
       sna_desc.rij[j_sna][0] = ei_displ[0][j];
       sna_desc.rij[j_sna][1] = ei_displ[1][j];
       sna_desc.rij[j_sna][2] = ei_displ[2][j];
       sna_desc.inside[j_sna] = j;
       sna_desc.wj[j_sna] = 1;
       sna_desc.rcutij[j_sna] = (radi + radelem[jtype]) * rcutfac;
       sna_desc.element[j_sna] = jelem;
   }
  }
  else{ //iat is an ion
    for (int j = 0; j< Nelec; j++){
       const auto je_displ = ei_displ[j];
       if (j == VP.refPtcl){
          const auto je_displ = VP.getDistTableAB(ei_Table_ID_).getDisplRow(j)[iat];// sign change here?
       }
      const auto disp_ref = VP.getRefPS().getDistTableAB(ei_Table_ID_).getDisplRow(iat)[j];
      int jtype = type_map[j];
      int jelem = 0;
      sna_desc.rij[j][0] = je_displ[0];
      sna_desc.rij[j][1] = je_displ[1];
      sna_desc.rij[j][2] = je_displ[2];
      sna_desc.inside[j] = j;
      sna_desc.wj[j] = 1;
      sna_desc.rcutij[j] = (radi + radelem[jtype]) * rcutfac;
      sna_desc.element[j] = jelem;
    }
    for (int j = 0; j< Nions; j++){
      int jtype = type_map[Nelec+j];
      int jelem = 0;
      int j_sna = Nelec+j;
      const auto& ii_displ = Ions.getDistTableAA(0).getDisplRow(iat);
      sna_desc.rij[j_sna][0] = ii_displ[j][0];
      sna_desc.rij[j_sna][1] = ii_displ[j][1];
      sna_desc.rij[j_sna][2] = ii_displ[j][2];
      sna_desc.inside[j_sna] = j;
      sna_desc.wj[j_sna] = 1;
      sna_desc.rcutij[j_sna] = (radi + radelem[jtype]) * rcutfac;
      sna_desc.element[j_sna] = jelem;
   }
  }
}


  inline void SNAJastrow::update_sna_rij(const ParticleSet& P, int iat){
    const int itype = type_map[iat];
    const double radi = radelem[itype];
    bool elec = (itype <2);
    int num_neigh = 0;// TODO: temporary fix for not treating cutoff
    if (elec){
        for (int j = 0; j < Nelec; j++){
            if (iat != j){
              const auto disp_ref = iat < j  ? -1.0*P.getDistTableAA(ee_Table_ID_).getDisplRow(j)[iat] : 1.0*P.getDistTableAA(ee_Table_ID_).getDisplRow(iat)[j];
              int jtype = type_map[j];
              int jelem = 0;
              app_debug() << "SNAJastrow::update_sna_rij xmove is "<< disp_ref[0] << std::endl;
              sna_desc.rij[num_neigh][0] = -disp_ref[0];
              sna_desc.rij[num_neigh][1] = -disp_ref[1];
              sna_desc.rij[num_neigh][2] = -disp_ref[2];
              sna_desc.inside[num_neigh] = j;
              sna_desc.wj[num_neigh] = 1;
              sna_desc.rcutij[num_neigh] = (radi + radelem[jtype]) * rcutfac;
              sna_desc.element[num_neigh] = jelem;
              num_neigh +=1;
            }
        }
        for (int j = 0; j< Nions; j++){
              const auto disp_ref = P.getDistTableAB(ei_Table_ID_).getDisplRow(iat)[j];
              int j_sna = Nelec+j;
              int jtype = type_map[j_sna];
              int jelem = 0;
              sna_desc.rij[num_neigh][0] = -disp_ref[0];
              sna_desc.rij[num_neigh][1] = -disp_ref[1];
              sna_desc.rij[num_neigh][2] = -disp_ref[2];
              sna_desc.inside[num_neigh] = j_sna;
              sna_desc.wj[num_neigh] = 1;
              sna_desc.rcutij[num_neigh] = (radi + radelem[jtype]) * rcutfac;
              sna_desc.element[num_neigh] = jelem;
              num_neigh +=1;
        }
    }
    else{// ion
        for (int j = 0; j< Nelec; j++){
              const auto disp_ref = -1.0*P.getDistTableAB(ei_Table_ID_).getDisplRow(j)[iat-Nelec]; // iat is global particle index, shift for just ion
              int jtype = type_map[j];
              int jelem = 0;
              app_debug() << "SNAJastrow::update_sna_rij disp_ref"<< disp_ref[0]<< std::endl;
              sna_desc.rij[num_neigh][0] = -disp_ref[0];
              sna_desc.rij[num_neigh][1] = -disp_ref[1];
              sna_desc.rij[num_neigh][2] = -disp_ref[2];
              sna_desc.inside[num_neigh] = j;
              sna_desc.wj[num_neigh] = 1;
              sna_desc.rcutij[num_neigh] = (radi + radelem[jtype]) * rcutfac;
              sna_desc.element[num_neigh] = jelem;
              num_neigh +=1;
        }
        for (int j = 0; j< Nions; j++){
          if (iat-Nelec != j){
              int jtype = type_map[Nelec+j];
              int jelem = 0;
              int j_sna = Nelec+j;
              const auto disp_ref = iat-Nelec < j  ? -1.0*Ions.getDistTableAA(ii_Table_ID_).getDisplRow(j)[iat-Nelec] : 1.0*Ions.getDistTableAA(ii_Table_ID_).getDisplRow(iat-Nelec)[j];
              sna_desc.rij[num_neigh][0] = -disp_ref[0];
              sna_desc.rij[num_neigh][1] = -disp_ref[1];
              sna_desc.rij[num_neigh][2] = -disp_ref[2];
              sna_desc.inside[num_neigh] = j+Nelec;
              sna_desc.wj[num_neigh] = 1;
              sna_desc.rcutij[num_neigh] = (radi + radelem[jtype]) * rcutfac;
              sna_desc.element[num_neigh] = jelem;
              num_neigh +=1;
          }
        }

    }
  }
  
  void SNAJastrow::compute_bispectrum(int iat){
    //TODO
    // will be responsible for looping over particles and updatin g the bispectrum of a particular component.
    app_debug() << "SNAJastrow::compute_bispectrum "<< std::endl;
    
    int ninside = Nelec + Nions-1; // just do alll particlles for now.
    int ielem = 0;
    sna_desc.compute_ui(ninside, ielem);
    sna_desc.compute_zi();
    sna_desc.compute_bi(ielem);
    for (int icoeff = 0; icoeff < ncoeff; icoeff++){
    //std::cout<< "desc object value"<<std::endl;
    std::cout<< sna_desc.blist[icoeff] <<std::endl;
      sna[iat][icoeff] = sna_desc.blist[icoeff];
  }
}

  void SNAJastrow::compute_d_dr_bispectrum( int iat){
  app_debug() << "SNAJastrow::compute_d_dr_bispectrum "<< std::endl;
  int ntotal = Nions + Nelec;


  // clear local array
  for (int i = 0; i < ntotal; i++)
    for (int entry = 0; entry < 3*ntypes*ncoeff; entry++) {
      snad[i][entry] = 0.0;
    }
  std::cout << "SNAJastrow::compute_d_dr_bispectrum after clear array"<< std::endl;

  // invoke full neighbor list (will copy or build if necessary)

  // compute sna derivatives for each atom in group
  // use full neighbor list to count atoms less than cutoff

  //double** const x = P.R;
  int inum = iat;
  int itype = type_map[iat];
  int ielem = 0;
  const double radi = radelem[itype]; 
  //TODO make atom type mapping struct
  const int typeoffset = 3*ncoeff*(itype);

  // insure rij, inside, and typej  are of size jnum

  // rij[][3] = displacements between atom I and those neighbors
  // inside = indices of neighbors of I within cutoff
  // typej = types of neighbors of I within cutoff
  // note Rij sign convention => dU/dRij = dU/dRj = -dU/dRi

  int ninside = Nions+Nelec-1;
  for (int jj = 0; jj < ninside; jj++) { //TODO: implement cutoff, currently going over all particles.
  if (iat !=jj){
    std::cout<< "iat: " << iat << " j " << jj <<std::endl;
  
    sna_desc.compute_ui(ninside, ielem);
    sna_desc.compute_zi();
    const int j = sna_desc.inside[jj];
    int jtype = type_map[jj];
    const int jtypeoffset = 3*ncoeff*(jtype);

    sna_desc.compute_duidrj(sna_desc.rij[jj], sna_desc.wj[jj],
                                sna_desc.rcutij[jj], ninside, sna_desc.element[jj]);
    std::cout << "made it after compute_duidrj" <<std::endl;
    sna_desc.compute_dbidrj();
    std::cout << "made it after compute_dbidrj" <<std::endl;
    ninside +=1;
    // Accumulate -dBi/dRi, -dBi/dRj

    int snadi_idx = typeoffset;
    int snadj_idx = jtypeoffset;
    int yoffset = ncoeff;
    int zoffset = 2*ncoeff;
    for (int icoeff = 0; icoeff < ncoeff; icoeff++) {
      std::cout<< "icoeff is " << icoeff <<std::endl;
      snad[iat][jtypeoffset + icoeff]           += sna_desc.dblist[icoeff][0];
      snad[iat][jtypeoffset + icoeff + yoffset] += sna_desc.dblist[icoeff][1];
      snad[iat][jtypeoffset + icoeff + zoffset] += sna_desc.dblist[icoeff][2];
      snad[jj][typeoffset + icoeff]           -= sna_desc.dblist[icoeff][0];
      snad[jj][typeoffset + icoeff + yoffset] -= sna_desc.dblist[icoeff][1];
      snad[jj][typeoffset + icoeff + zoffset] -= sna_desc.dblist[icoeff][2];
    }
  }// end neighbor loop
}
}//TODO}



  void SNAJastrow::evaluateDerivRatios(const VirtualParticleSet& VP,const opt_variables_type& optvars, std::vector<ValueType>& ratios, Matrix<ValueType>& dratios){
    app_debug() << "inside evaluateDerivRatios" << std::endl;
    evaluateRatios(VP,ratios);
    Vector<RealType> dlogpsi_nlpp_ref;
    dlogpsi_nlpp_ref.resize(myVars.size());
    Vector<RealType> dlogpsi_nlpp_virt;
    dlogpsi_nlpp_virt.resize(myVars.size());

    bool recalculate(false);
    std::vector<bool> rcsingles(myVars.size(), false);
    for (int k = 0; k < myVars.size(); ++k){
      int k_global = myVars.where(k);
      if (k_global < 0)
        continue;
      if (optvars.recompute(k_global))
        recalculate = true;
      rcsingles[k] = true;
     }

    if (recalculate){
     for (int k = 0; k < myVars.size(); k++){ // local index
       int k_global = myVars.where(k); //global index
       if (k_global < 0)
         continue;
        if (rcsingles[k]){
         int ntype = int(k/ncoeff); 
         int coeff = k%ncoeff; 
         //calculate reference deriv
         // are all particles up to date?
         for (int i = 0 ; i < Nelec; i ++){
           update_sna_rij(VP.getRefPS(), i); 
           compute_bispectrum(i);
         }
         //linear derivs for ref
         if (ntype < VP.getRefPS().groups()){ // if 0 or 1 it is up or down electrons
            for (int iat = VP.getRefPS().first(ntype); iat < VP.getRefPS().last(ntype); iat++) { // loop over elements in each group
               dlogpsi_nlpp_ref[k] += -sna[iat][coeff];
              }
         }
        else{
          for (int iat = Ions.first(ntype); iat < Ions.last(ntype); iat++) { // loop over elements in each group
               dlogpsi_nlpp_ref[k] += -sna[VP.getRefPS().getTotalNum()+iat][coeff];
            }
        }


         //* perform sampling of positions. *//
         for (int r = 0; r < ratios.size(); r++){
            // update the descriptors according to this new particle moving
            for (int i = 0 ; i < Nelec; i ++){
                update_sna_rij_vp(VP,i);
                compute_bispectrum(i);
            }
           // if the coefficient belongs to an electron group
           if (ntype < VP.getRefPS().groups()){
              dlogpsi_nlpp_virt[k] = 0;
              for (int iat = VP.getRefPS().first(ntype); iat < VP.getRefPS().last(ntype); iat++) { // loop over elements in each group
                dlogpsi_nlpp_virt[k] = -sna[iat][coeff];
              }
           }
          // if the coefficient belongs to an ion group
           else{
                for (int iat = Ions.first(ntype); iat < Ions.last(ntype); iat++) { // loop over elements in each group
                   dlogpsi_nlpp_ref[k] = -sna[VP.getRefPS().getTotalNum()+iat][coeff];
                }
           }
           dratios[r][k_global] =  dlogpsi_nlpp_virt[k] - dlogpsi_nlpp_ref[k];
         } //end ratio loop
        }// end rcsingles
     } // end loop over internal coeffss
    } // end recalculate
  }

  void SNAJastrow::evaluateRatios(const VirtualParticleSet& VP, std::vector<ValueType>& ratios){
    app_debug() << "inside evaluateRatios" << std::endl;
    ScopedTimer local_timer(timers_.eval_ratio_timer);
    double Eold, Enew;
    for (int i = 0 ; i < Nelec; i ++){
      update_sna_rij(VP.getRefPS(), i); 
      compute_bispectrum(i);
    }
    calculate_ESNA(VP.getRefPS(), snap_beta, Eold);
    for (int r = 0; r < ratios.size(); r++){
        for (int e=0 ; e < Nelec; e++){
          update_sna_rij_vp(VP,e);  
          compute_bispectrum(e);
        }
      //calculate Enew
      calculate_ESNA(VP.getRefPS(), snap_beta, Enew);
      //store ratio
      ratios[r] = std::exp(static_cast<ValueType>(Enew-Eold));
    }
    return;
  }



  /////////////////////////////////// MC Related functions /////////
  void SNAJastrow::acceptMove(ParticleSet& P, int iat, bool safe_to_delay){
    app_debug() << "inside accept" << std::endl;
    for (int e=0 ; e< Nelec; e++){
      update_sna_rij(P, e);  
      compute_bispectrum(e);
      compute_d_dr_bispectrum(e);
    }
    double esnap;
    calculate_ESNA(P, snap_beta, esnap);
    current_esnap=esnap;
    log_value_ = static_cast<SNAJastrow::LogValue>(esnap);
    computeGL(P,iat);
  }

 void SNAJastrow::restore(int iat){
 }

  void SNAJastrow::registerData(ParticleSet& P, WFBufferType& buf){
  }

  SNAJastrow::LogValue SNAJastrow::updateBuffer(ParticleSet& P, WFBufferType& buf, bool from_scratch){

    app_debug() << "inside update buffer" << std::endl;

      log_value_ = evaluateLog(P, P.G, P.L);
      return log_value_;
  }

  void SNAJastrow::copyFromBuffer(ParticleSet& P, WFBufferType& buf){
    
  }

  SNAJastrow::PsiValue SNAJastrow::ratioGrad(ParticleSet& P, int iat, GradType& grad_iat){
    app_debug() << "inside ratiograd" << std::endl;
    // TODO how to get proposed move distance tables
    for (int e=0 ; e< Nelec; e++){
      update_sna_rij(P, e);  
      compute_bispectrum(e);
    }
    update_sna_rij(P, iat);  
    compute_d_dr_bispectrum(iat);// compute snad fur current particle 
    double Enew, Eold;
    int row, col;
    calculate_ESNA(P, snap_beta, Enew);
    for (int dim = 0; dim < 3; dim++){
      for (int k = 0; k < ncoeff ; k++){
        for (int n = 0; n < ntypes; n++){
          int col = (n*(3*ncoeff)) + (dim*ncoeff)+k;
          grad_iat[dim] += snap_beta[n][k]*snad[iat][col];
        }
      }
    }
    for (int e=0 ; e< Nelec; e++){
      update_sna_rij(P, e);  
      compute_bispectrum(e);
    }
    calculate_ESNA(P, snap_beta, Eold);
    SNAJastrow::PsiValue ratio = std::exp(static_cast<SNAJastrow::PsiValue>(Enew-Eold));
    return ratio;
  }

  SNAJastrow::PsiValue SNAJastrow::ratio(ParticleSet& P, int iat){
    app_debug() << "inside ratio" << std::endl;
    double Enew, Eold;
    for (int e=0 ; e< Nelec; e++){
      update_sna_rij(P, e);  
      compute_bispectrum(e);
    }
    calculate_ESNA(P, snap_beta, Eold);
    //Eold = current_esnap;
    for (int e=0 ; e< Nelec; e++){
      update_sna_rij(P, e);  
      compute_bispectrum(e);
    }
    calculate_ESNA( P, snap_beta, Enew);
    //calculate the ratio
    SNAJastrow::PsiValue ratio = std::exp(static_cast<SNAJastrow::PsiValue>(Enew-Eold));
    return ratio;
  }

void SNAJastrow::extractOptimizableObjectRefs(UniqueOptObjRefs& opt_obj_refs){opt_obj_refs.push_back(*this);}

void SNAJastrow::checkInVariablesExclusive(opt_variables_type& active){
  myVars.setIndexDefault(); // I don't actually know what this is doing?
  active.insertFrom(myVars);
}

void SNAJastrow::checkOutVariables(const opt_variables_type& active ){
    myVars.getIndex(active);
  }

void SNAJastrow::resetParametersExclusive(const opt_variables_type& active){
   app_debug() << "inside reset params" << std::endl;
  int k_global, ntype, coeff;
  for (int i=0; i < myVars.size(); i++){
    k_global = myVars.where(i);
    if (k_global >=0){
     ntype = int(i/ncoeff);
     coeff = i%ncoeff; 
     snap_beta[ntype][coeff] = myVars[i] = active[k_global];
    }
    }
  }

 
std::unique_ptr<WaveFunctionComponent> SNAJastrow::makeClone(ParticleSet& tpq) const
{
  auto snap_copy = std::make_unique<SNAJastrow>(std::string("snap"), Ions, tpq, std::string("linear"), twojmax, rcut);
  snap_copy->snap_beta = snap_beta;
  return snap_copy;
}


bool SNAJastrow::put(xmlNodePtr cur) {
  
  app_summary() << "     Number of parameters: " << myVars.size() << std::endl;
  for (int i = 0; i < myVars.size(); i++){
    app_summary() << myVars[i] <<std::endl;
  }
  return true;
}

void SNAJastrow::createResource(ResourceCollection& collection) const
{
  collection.addResource(std::make_unique<SNAMultiWalkerMem<RealType>>());
}

void SNAJastrow::acquireResource(ResourceCollection& collection,
                                       const RefVectorWithLeader<WaveFunctionComponent>& wfc_list) const
{
  auto& wfc_leader          = wfc_list.getCastedLeader<SNAJastrow>();
  wfc_leader.mw_mem_handle_ = collection.lendResource<SNAMultiWalkerMem<RealType>>();
}

void SNAJastrow::releaseResource(ResourceCollection& collection,
                                       const RefVectorWithLeader<WaveFunctionComponent>& wfc_list) const
{
  auto& wfc_leader = wfc_list.getCastedLeader<SNAJastrow>();
  collection.takebackResource(wfc_leader.mw_mem_handle_);
}


}

