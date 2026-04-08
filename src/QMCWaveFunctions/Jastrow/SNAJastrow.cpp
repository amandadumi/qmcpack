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
    sna_desc.init(Nelec+Nions, input_twojmax, input_rcut);
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
    num_neigh= Nions+Nelec-1;
    snap_type = input_snap_type;
    //TODO: set this to ddefault value for everythin
    rcutij = std::vector<std::vector<double>>(Nions+Nelec,std::vector<double>(Nions+Nelec,7.0));
    radelem = std::vector<double>(ntypes,rcut/2);
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
    snad = std::vector<std::vector<double>>(Nions+Nelec, std::vector<double>(3*ntypes*ncoeff,0.0));
    for (int i=0; i < NIonGroups+els.groups(); i++){
      for (int k = 0; k < ncoeff; k++){
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

SNAJastrow::~SNAJastrow(){}

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


 double SNAJastrow::full_FD_Lap(const ParticleSet& P, int iat, std::vector<std::vector<double>> coeffs){
  app_debug() << "SNAPJastrow::FD_Lap entered function " <<std::endl;

  double finite_diff_lap=0;
  double this_coeff; 
  int col;
  std::vector<std::vector<double>>temp_snad_forward(snad);
  std::vector<std::vector<double>>temp_snad_backward(snad);
  double G_finite_diff_forward;
  double G_finite_diff_back;
  // if taking wave function laplacian, we just need bispectrum component but for electron laplacian we use the snap coefficients.

  //forward direction
  // clear snad array
  for (int dim = 0; dim < 3; dim++){
    for (int par = 0; par < Nions + Nelec; par++){
      temp_snad_forward[par].assign(temp_snad_forward[par].size(),0.0);
      temp_snad_backward[par].assign(temp_snad_backward[par].size(),0.0);
    }

    // create the descriptor for all particles. We can optimize this, but this is lazy way for now. should reduce this to just loop over particles in ntype.
    for (int par = 0; par < Nions + Nelec; par++){
      update_sna_rij(P, par, false);// update sna_rij  
      if (par == iat) // if the rij we are creating is from par == iat, shift all of rij
        for (int j= 0; j < num_neigh; j++)
          sna_desc.rij[j][dim] += dist_delta;
      else // if rij par !=iat then we just need to update the row corresponding to iat.
        for (int j=0; j < num_neigh; j++)
            if (sna_desc.inside[j] == iat)
              sna_desc.rij[j][dim] -= dist_delta; // switch sign to go from other perspective.
      // change the sign for rij for the snad contribution.
      for (int i = 0; i < num_neigh; i ++)
        for (int j = 0; j < 3; j ++)
          sna_desc.rij[i][j]*=-1.0;
      compute_d_dr_bispectrum(par, temp_snad_forward);
    }

    for (int par = 0; par < Nions + Nelec; par++){
      update_sna_rij(P, par, false); // update sna_rij  
      if (par == iat) // if the rij we are creating is from par == iat, shift all of rij
        for (int j= 0; j < num_neigh; j++)
          sna_desc.rij[j][dim] -= dist_delta;
      else // if rij par !=iat then we just need to update the row corresponding to iat.
        for (int j= 0; j < num_neigh; j++)
            if (sna_desc.inside[j] == iat)
              sna_desc.rij[j][dim] += dist_delta; // switch sign to go from other perspective.
      // change the sign for rij for the snad contribution.
      for (int i = 0; i < num_neigh; i++)
        for (int j = 0; j < 3; j ++)
          sna_desc.rij[i][j]*=-1.0;
      compute_d_dr_bispectrum(par, temp_snad_backward);
    }
    //loop over types of particles (electrons and ions)
    for (int n=0; n < ntypes; n++){
      // loop over the components
      for (int k = 0; k < ncoeff; k ++){
        //we wil need gradient in each direction.
        col = (n*(3*ncoeff)) + (dim*ncoeff) + k;
        this_coeff = coeffs[n][k];
        G_finite_diff_forward = this_coeff * temp_snad_forward[iat][col];
        // back to normal sign for next bispectrum
        G_finite_diff_back = this_coeff * temp_snad_backward[iat][col];
        finite_diff_lap += (G_finite_diff_forward - G_finite_diff_back)/(2*dist_delta); 
      }
    }
  } // end dim


  app_debug() << "SNAPJastrow::FD_Lap at end of function " <<std::endl;
  return finite_diff_lap;
 }




 double SNAJastrow::FD_Lap(const ParticleSet& P, int iat, int dim, int coeff, int ntype){
  app_debug() << "SNAPJastrow::FD_Lap entered function " <<std::endl;
  
  std::vector<std::vector<double>>temp_snad(snad);
  double G_finite_diff_forward;
  double G_finite_diff_back;
  // if taking wave function laplacian, we just need bispectrum component but for electron laplacian we use the snap coefficients.
  int col = (ntype*(3*ncoeff))+(dim*ncoeff)+coeff;

  //forward direction
  // clear snad array
  for (int par = 0; par < Nions + Nelec; par++)
    temp_snad[par].assign(temp_snad[par].size(),0.0);
  // create the descriptor for all particles. We can optimize this, but this is lazy way for now. should reduce this to just loop over particles in ntype.
  for (int par = 0; par < Nions + Nelec; par++){
    update_sna_rij(P, par, false);// update sna_rij  
    if (par == iat) // if the rij we are creating is from par == iat, shift all of rij
      for (int j= 0; j < num_neigh; j++)
        sna_desc.rij[j][dim] += dist_delta;
    else // if rij par !=iat then we just need to update the row corresponding to iat.
      for (int j= 0; j < num_neigh; j++)
          if (sna_desc.inside[j] == iat)
            sna_desc.rij[j][dim] -= dist_delta; // switch sign to go from other perspective.
    // change the sign for rij for the snad contribution.
    for (int i = 0; i < num_neigh; i++)
      for (int j = 0; j < 3; j ++)
        sna_desc.rij[i][j]*=-1.0;
    compute_d_dr_bispectrum(par, temp_snad);
  }
  G_finite_diff_forward = temp_snad[iat][col];


  //backward direction
  for (int par = 0; par < Nelec+Nions; par++) 
    temp_snad[par].assign(temp_snad[par].size(),0.0);

  for (int par = 0; par < Nions + Nelec; par++){
    update_sna_rij(P, par,false); // update sna_rij  
    if (par == iat) // if the rij we are creating is from par == iat, shift all of rij
      for (int j= 0; j < num_neigh; j++)
        sna_desc.rij[j][dim] -= dist_delta;
    else // if rij par !=iat then we just need to update the row corresponding to iat.
      for (int j= 0; j < num_neigh; j++)
          if (sna_desc.inside[j] == iat)
            sna_desc.rij[j][dim] += dist_delta; // switch sign to go from other perspective.
    // change the sign for rij for the snad contribution.
    for (int i = 0; i < num_neigh; i++)
      for (int j = 0; j < 3; j ++)
        sna_desc.rij[i][j]*=-1.0;
    compute_d_dr_bispectrum(par,temp_snad);
  }
  G_finite_diff_back = temp_snad[iat][col];

  double finite_diff_lap = (G_finite_diff_forward - G_finite_diff_back)/(2*dist_delta); 


  app_debug() << "SNAPJastrow::FD_Lap at end of function " <<std::endl;
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
    app_debug() << "SNAJastrow::computeGL entered function" <<std::endl;
    ScopedTimer local_timer(timers_.eval_gl_timer);
    int row, col;
    row = iel; // get the derivative with respect to r of a specific particle. 
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
          col = (n*(3*ncoeff)) + (dim*ncoeff) + k;
          //app_debug() << "grad val is " << grad_val << std::endl;
          grad_u[iel][dim] += snap_beta[n][k] * snad[row][col];
        } // end dim loop
      }// end k loop
    } //end n loop
    lap_u[iel] = full_FD_Lap(P, iel, snap_beta);
    return;
   }


  SNAJastrow::LogValue SNAJastrow::evaluateLog(const ParticleSet& P,
                                    ParticleSet::ParticleGradient& G,
                                    ParticleSet::ParticleLaplacian& L){
    app_debug() << "SNAJastrow::evaluateLog entered function" <<std::endl;
    ScopedTimer local_timer(timers_.eval_log_timer);
    double esnap;
    calculate_ESNA(P, snap_beta, esnap, false);
    current_esnap = esnap;
    log_value_ = static_cast<SNAJastrow::LogValue>(esnap);
    for (int par = 0; par < Nelec+Nions; par++) 
        snad[par].assign(snad[par].size(),0.0);
    //update all of snad
    for (int par = 0; par < Nions + Nelec; par++){
      update_sna_rij(P, par, false);  
      for (int i = 0; i < num_neigh; i ++)
        for (int j = 0; j < 3; j ++)
          sna_desc.rij[i][j]*=-1.0;
      compute_d_dr_bispectrum(par,snad);
    }
    for (int e=0 ; e< Nelec; e++){
        computeGL(P,e);
        G[e] += grad_u[e];
        L[e] += lap_u[e];
    }
    app_debug() << "SNAJastrow::evaluateLog evalutated GL" <<std::endl;
    return log_value_;
  }

  SNAJastrow::GradType SNAJastrow::evalGrad(ParticleSet& P, int iat){
  
    std::vector<std::vector<double>>temp_snad(snad);
    GradType grad_iat;

    for (int par = 0; par < Nions + Nelec; par++)
      temp_snad[par].assign(temp_snad[par].size(),0.0);
      
    for (int par = 0; par < Nions + Nelec; par++){
      update_sna_rij(P, par, true);  
      for (int i = 0; i < num_neigh; i ++)
        for (int j = 0; j < 3; j ++)
          sna_desc.rij[i][j]*=-1.0;
      compute_d_dr_bispectrum(par, temp_snad);
    }

    for (int dim=0; dim < 3; dim++){
     int row = (3*iat) + dim;
     for (int k = 0; k < ncoeff ; k++){
       for (int n = 0; n < ntypes; n++){
         int col = (n * (3*ncoeff) ) + (dim*ncoeff) + k;
         grad_iat[dim] += snap_beta[n][k]*temp_snad[iat][col];
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
    for (int k = 0; k < myVars.size(); ++k)
    {
        int k_global = myVars.where(k);
        dhpsioverpsi[k_global] = -RealType(0.5) * RealType(Sum(lapLogPsi[k]));
        for (int i = 0; i < Nelec; i++)
          dhpsioverpsi[k_global] -= RealType(dot(P.G[i], gradLogPsi[k][i]));
    }
  }


  void SNAJastrow::evaluateDerivativesWF(ParticleSet& P, const opt_variables_type& optvars, Vector<ValueType>& dlogpsi){
    app_debug() << "in evalderivativeWF"<<std::endl;
    for (int par = 0; par < Nelec+Nions; par++) 
        snad[par].assign(snad[par].size(),0.0);

    for (int par = 0; par < Nions + Nelec; par++){
      update_sna_rij(P, par, false);  
      compute_bispectrum(par);
      for (int i = 0; i < num_neigh; i++)
        for (int j = 0; j < 3; j ++)
          sna_desc.rij[i][j]*=-1.0;
      compute_d_dr_bispectrum(par, snad);
    }
    resizeWFOptVectors();
    const size_t NumVars = myVars.size();
    for (int p = 0; p < NumVars; p++){
      gradLogPsi[p] = 0.0;
      lapLogPsi[p] = 0.0;
    }
    dLogPsi = 0.0;
    
    for (int k = 0; k < myVars.size(); k++){
        int k_global = myVars.where(k);
        evaluate_linear_derivs(P, k);
        dlogpsi[k_global] = ValueType(dLogPsi[k]);
    } 
  }

  void SNAJastrow::evaluate_linear_derivs(ParticleSet& P, int coeff_idx){
    app_debug() << "in linearderivs" <<std::endl;
    
    int ntype = int(coeff_idx/ncoeff);  
    int coeff = coeff_idx%ncoeff; //which coeff of this type are we on.
    if (ntype < P.groups()){ // are we considering electron derivatives?
      for (int iat = P.first(ntype); iat < P.last(ntype); iat++) // loop over elements in electron group
        dLogPsi[coeff_idx] += -sna[iat][coeff];
    }
    else{
     for (int iat = Ions.first(ntype-P.groups()); iat < Ions.last(ntype-P.groups()); iat++) // loop over elements in each group
       dLogPsi[coeff_idx] += -sna[P.getTotalNum()+iat][coeff];
    }
    for (int iel =0; iel < Nelec; iel++)
      for (int dim = 0; dim < OHMMS_DIM; dim++){ // loop over dim to get grad vec.
       int col = (ntype*(3*ncoeff))+(dim*ncoeff)+coeff;
       gradLogPsi[coeff_idx][iel][dim] += snad[iel][col];
       lapLogPsi[coeff_idx][iel] += FD_Lap(P, iel, dim, coeff, ntype);
      }
    }


  void SNAJastrow::calculate_ESNA(const ParticleSet& P, const std::vector<std::vector<double>> coeff, double& new_u, bool proposed){
    app_debug() << "SNAJastrow::calculate_ESNA entered function" <<std::endl;
    ScopedTimer local_timer(timers_.eval_esnap_timer);
    double esnap_all=0;
    double esnap_elec=0;
    double esnap_ion=0;
    double bispectrum_val;
    // calculate electron contribution
    // the global array is summed over groups of atoms of the same type. thus we just need to sum over groups.
    for (int ig = 0; ig < P.groups(); ig++) {
      for (int iat = P.first(ig); iat < P.last(ig); iat++) { // loop over elements in each group
        app_debug() << "elec index is" << iat <<std::endl; 
        update_sna_rij(P, iat, proposed); 
        app_debug() << "num neigh is" << num_neigh <<std::endl; 
        compute_bispectrum(iat);
        for (int k = 0; k < ncoeff; k++){
          bispectrum_val = sna[iat][k]; //block of bispectrum + current component to add.
          app_debug() << bispectrum_val << std::endl;
          esnap_elec += coeff[ig][k] * bispectrum_val;
          app_debug()<<"snap coefff for group "<< ig << " coeff "<<k << " is " <<coeff[ig][k] <<std::endl;
        }
      }
    }
    esnap_all += esnap_elec;
    app_debug() << "SNAJastrow::calculate_ESNA calc ion center contribution" <<std::endl;

    for (int ig = 0; ig < Ions.groups(); ig++) {
      for (int iat = Ions.first(ig); iat < Ions.last(ig); iat++) { // loop over elements in each group
        app_debug() << "ion index is" << iat <<std::endl; 
        update_sna_rij(P, Nelec+iat,proposed); 
        compute_bispectrum(Nelec+iat);
        for (int k =0; k < ncoeff; k++){
          bispectrum_val = sna[Nelec+iat][k];
          app_debug() << bispectrum_val << std::endl;
          esnap_ion += coeff[P.groups()+ig][k] * bispectrum_val;
          app_debug()<< "snap coeff for group "<< P.groups()+ig << " coeff " << k << " is " << coeff[P.groups()+ig][k] <<std::endl;
        }
      }
    }
    esnap_all += esnap_ion;

    new_u = -esnap_all;
    return;
  }

  inline void SNAJastrow::update_sna_rij_vp(const VirtualParticleSet& VP, int iat, int r){
    const int itype = type_map[iat];
    const double radi = radelem[itype];
    bool elec = (itype <2);
    PosType disp_ref;
    RealType dist_ref;
    num_neigh = 0;
    if (elec){ // if up or down elec
        for (int j = 0; j < Nelec; j++){ // loop through all other elecs
            if (iat != j){ // skip current elec
               // Displacements stored in lower triangular of ee table
               disp_ref = iat < j  ? VP.getRefPS().getDistTableAA(ee_Table_ID_).getDisplRow(j)[iat] : -1.0*VP.getRefPS().getDistTableAA(ee_Table_ID_).getDisplRow(iat)[j];
               dist_ref = iat < j  ? VP.getRefPS().getDistTableAA(ee_Table_ID_).getDistRow(j)[iat] : VP.getRefPS().getDistTableAA(ee_Table_ID_).getDistRow(iat)[j];
               if (j == VP.refPtcl){// if the neighboring electron is the virtual particle that is being moved around...
                 // get the dispalcement between the virtual particlle and the target particle
                 // we need the distance from i to r, which is this value no sign change
                 disp_ref = VP.getDistTableAB(ee_Table_ID_).getDisplRow(r)[iat];// sign change here since we are looking from perspective of i but grabbing displacement from perspective of refptcl j
                 dist_ref = VP.getDistTableAB(ee_Table_ID_).getDistRow(r)[iat];// sign change here since we are looking from perspective of i but grabbing displacement from perspective of refptcl j
               }
              else if (iat == VP.refPtcl){ // if the target is the ref particle...
                 // then it is the distance between target and neighbor j 
                 disp_ref = -VP.getDistTableAB(ee_Table_ID_).getDisplRow(r)[j];
                 dist_ref = VP.getDistTableAB(ee_Table_ID_).getDistRow(r)[j];
              }
              if (dist_ref < rcut){
                int jtype = type_map[j];
                int jelem = 0;
                sna_desc.rij[num_neigh][0] = disp_ref[0];
                sna_desc.rij[num_neigh][1] = disp_ref[1];
                sna_desc.rij[num_neigh][2] = disp_ref[2];
                sna_desc.inside[num_neigh] = j;
                sna_desc.wj[num_neigh] = 1;
                sna_desc.rcutij[num_neigh] = (radi + radelem[jtype]) * rcutfac;
                sna_desc.element[num_neigh] = jelem;
                num_neigh +=1;
              }
            }
        }
        for (int j = 0; j< Nions; j++){
              disp_ref = -VP.getRefPS().getDistTableAB(ei_Table_ID_).getDisplRow(iat)[j];
              dist_ref = VP.getRefPS().getDistTableAB(ei_Table_ID_).getDistRow(iat)[j];
              if (iat == VP.refPtcl){
                disp_ref = -VP.getDistTableAB(ei_Table_ID_).getDisplRow(r)[j];
                dist_ref = VP.getDistTableAB(ei_Table_ID_).getDistRow(r)[j];
              }
              if (dist_ref < rcut){
                int j_sna = Nelec+j;
                int jtype = type_map[j_sna];
                int jelem = 0;
                sna_desc.rij[num_neigh][0] = disp_ref[0];
                sna_desc.rij[num_neigh][1] = disp_ref[1];
                sna_desc.rij[num_neigh][2] = disp_ref[2];
                sna_desc.inside[num_neigh] = j_sna;
                sna_desc.wj[num_neigh] = 1;
                sna_desc.rcutij[num_neigh] = (radi + radelem[jtype]) * rcutfac;
                sna_desc.element[num_neigh] = jelem;
                num_neigh +=1;
              }
        }
    }
    else{// ion
        for (int j = 0; j< Nelec; j++){ // for all neighboring electrons of this ion
                disp_ref = VP.getRefPS().getDistTableAB(ei_Table_ID_).getDisplRow(j)[iat-Nelec]; // iat is global particle index, shift for just ion
                dist_ref = VP.getRefPS().getDistTableAB(ei_Table_ID_).getDistRow(j)[iat-Nelec]; // iat is global particle index, shift for just ion
              if (j == VP.refPtcl){
                 disp_ref = VP.getDistTableAB(ei_Table_ID_).getDisplRow(r)[iat-Nelec];// sign change here since we are looking from perspective of i but grabbing displacement from perspective of refptcl j
                 dist_ref = VP.getDistTableAB(ei_Table_ID_).getDistRow(r)[iat-Nelec];// sign change here since we are looking from perspective of i but grabbing displacement from perspective of refptcl j
              }
              if (dist_ref < rcut){
                int jtype = type_map[j];
                int jelem = 0;
                sna_desc.rij[num_neigh][0] = disp_ref[0];
                sna_desc.rij[num_neigh][1] = disp_ref[1];
                sna_desc.rij[num_neigh][2] = disp_ref[2];
                sna_desc.inside[num_neigh] = j;
                sna_desc.wj[num_neigh] = 1;
                sna_desc.rcutij[num_neigh] = (radi + radelem[jtype]) * rcutfac;
                sna_desc.element[num_neigh] = jelem;
                num_neigh +=1;
          }
        }
        for (int j = 0; j< Nions; j++){
          if (iat-Nelec != j){
              int jtype = type_map[Nelec+j];
              int jelem = 0;
              int j_sna = Nelec+j;
              disp_ref = iat-Nelec < j  ? -Ions.getDistTableAA(ii_Table_ID_).getDisplRow(j)[iat-Nelec] : Ions.getDistTableAA(ii_Table_ID_).getDisplRow(iat-Nelec)[j];
              dist_ref = iat-Nelec < j  ? Ions.getDistTableAA(ii_Table_ID_).getDistRow(j)[iat-Nelec] : Ions.getDistTableAA(ii_Table_ID_).getDistRow(iat-Nelec)[j];
              if (dist_ref < rcut){
                sna_desc.rij[num_neigh][0] = disp_ref[0];
                sna_desc.rij[num_neigh][1] = disp_ref[1];
                sna_desc.rij[num_neigh][2] = disp_ref[2];
                sna_desc.inside[num_neigh] = j_sna;
                sna_desc.wj[num_neigh] = 1;
                sna_desc.rcutij[num_neigh] = (radi + radelem[jtype]) * rcutfac;
                sna_desc.element[num_neigh] = jelem;
                num_neigh +=1;
              }
          }
        }
    }
  }


  inline void SNAJastrow::update_sna_rij(const ParticleSet& P, int iat, bool proposed){
    const int itype = type_map[iat];
    const double radi = radelem[itype];
    bool elec = (itype <2);
    num_neigh = 0;
    PosType disp_ref;
    RealType dist_ref;

    if (elec){
        for (int j = 0; j < Nelec; j++){
            if (iat != j){
              disp_ref = iat < j  ? P.getDistTableAA(ee_Table_ID_).getDisplRow(j)[iat] : -1.0*P.getDistTableAA(ee_Table_ID_).getDisplRow(iat)[j];
              dist_ref = iat < j  ? P.getDistTableAA(ee_Table_ID_).getDistRow(j)[iat] : P.getDistTableAA(ee_Table_ID_).getDistRow(iat)[j];
              if (proposed){
                if (iat == P.getActivePtcl()){
                  disp_ref = -1*P.getDistTableAA(ee_Table_ID_).getTempDispls()[j] ;
                  dist_ref = P.getDistTableAA(ee_Table_ID_).getTempDists()[j] ;
                }
                else if (j == P.getActivePtcl()){
                  disp_ref = P.getDistTableAA(ee_Table_ID_).getTempDispls()[iat] ;
                  dist_ref = P.getDistTableAA(ee_Table_ID_).getTempDists()[iat] ;
                }
              }        //std::cout << "x component of disp ref is" <<disp_ref[0] <<std::endl;
              if (dist_ref< rcut){
                int jtype = type_map[j];
                int jelem = 0;
                sna_desc.rij[num_neigh][0] = disp_ref[0];
                sna_desc.rij[num_neigh][1] = disp_ref[1];
                sna_desc.rij[num_neigh][2] = disp_ref[2];
                sna_desc.inside[num_neigh] = j;
                sna_desc.wj[num_neigh] = 1;
                sna_desc.rcutij[num_neigh] = (radi + radelem[jtype]) * rcutfac;
                sna_desc.element[num_neigh] = jelem;
                num_neigh +=1;
              }
            }
            }
            for (int j = 0; j< Nions; j++){
              disp_ref = -1.0*P.getDistTableAB(ei_Table_ID_).getDisplRow(iat)[j];
              dist_ref = P.getDistTableAB(ei_Table_ID_).getDistRow(iat)[j];
              if (proposed) {
                if (iat == P.getActivePtcl()){
                  disp_ref = -1.0*P.getDistTableAB(ei_Table_ID_).getTempDispls()[j] ;
                  dist_ref = P.getDistTableAB(ei_Table_ID_).getTempDists()[j] ;
                }
              }
              if (dist_ref < rcut){
                int j_sna = Nelec+j;
                int jtype = type_map[j_sna];
                int jelem = 0;
                sna_desc.rij[num_neigh][0] = disp_ref[0];
                sna_desc.rij[num_neigh][1] = disp_ref[1];
                sna_desc.rij[num_neigh][2] = disp_ref[2];
                sna_desc.inside[num_neigh] = j_sna;
                sna_desc.wj[num_neigh] = 1;
                sna_desc.rcutij[num_neigh] = (radi + radelem[jtype]) * rcutfac;
                sna_desc.element[num_neigh] = jelem;
                num_neigh +=1;
              }
        }
    }
    else{// ion
        for (int j = 0; j< Nelec; j++){
              disp_ref = P.getDistTableAB(ei_Table_ID_).getDisplRow(j)[iat-Nelec];
              dist_ref = P.getDistTableAB(ei_Table_ID_).getDistRow(j)[iat-Nelec];
              if (proposed){
                if (j == P.getActivePtcl()){
                  disp_ref = P.getDistTableAB(ei_Table_ID_).getTempDispls()[iat-Nelec] ;
                  dist_ref = P.getDistTableAB(ei_Table_ID_).getTempDists()[iat-Nelec] ;
                }
              }
              int jtype = type_map[j];
              int jelem = 0;
              if (dist_ref < rcut){
                sna_desc.rij[num_neigh][0] = disp_ref[0];
                sna_desc.rij[num_neigh][1] = disp_ref[1];
                sna_desc.rij[num_neigh][2] = disp_ref[2];
                sna_desc.inside[num_neigh] = j;
                sna_desc.wj[num_neigh] = 1;
                sna_desc.rcutij[num_neigh] = (radi + radelem[jtype]) * rcutfac;
                sna_desc.element[num_neigh] = jelem;
                num_neigh +=1;
              }
        }
        for (int j = 0; j< Nions; j++){
          if (iat-Nelec != j){
              int jtype = type_map[Nelec+j];
              int jelem = 0;
              int j_sna = Nelec+j;
              disp_ref = iat-Nelec < j  ? -1.0*Ions.getDistTableAA(ii_Table_ID_).getDisplRow(j)[iat-Nelec] : Ions.getDistTableAA(ii_Table_ID_).getDisplRow(iat-Nelec)[j];
              dist_ref = iat-Nelec < j  ? Ions.getDistTableAA(ii_Table_ID_).getDistRow(j)[iat-Nelec] : Ions.getDistTableAA(ii_Table_ID_).getDistRow(iat-Nelec)[j];
              if (dist_ref < rcut){
                sna_desc.rij[num_neigh][0] = disp_ref[0];
                sna_desc.rij[num_neigh][1] = disp_ref[1];
                sna_desc.rij[num_neigh][2] = disp_ref[2];
                sna_desc.inside[num_neigh] = j_sna;
                sna_desc.wj[num_neigh] = 1;
                sna_desc.rcutij[num_neigh] = (radi + radelem[jtype]) * rcutfac;
                sna_desc.element[num_neigh] = jelem;
                num_neigh +=1;
              }
        }
      }
    }
  }
  
  void SNAJastrow::compute_bispectrum(int iat){
    int ninside = num_neigh;
    int ielem = 0;
    sna_desc.compute_ui(ninside, ielem);
    sna_desc.compute_zi();
    sna_desc.compute_bi(ielem);
    for (int icoeff = 0; icoeff < ncoeff; icoeff++){
      sna[iat][icoeff] = sna_desc.blist[icoeff];
    }
  }

  void SNAJastrow::compute_d_dr_bispectrum( int iat, std::vector<std::vector<double>>& a_snad)
  {
    //app_debug() << "SNAJastrow::compute_d_dr_bispectrum "<< std::endl;
    int inum = iat;
    int itype = type_map[iat];
    int ielem = 0;
    const int typeoffset = 3*ncoeff*itype;

    int ninside = num_neigh;
    sna_desc.compute_ui(ninside, ielem);
    sna_desc.compute_zi();
    for (int jj = 0; jj < ninside; jj++) { 
      const int j = sna_desc.inside[jj];
      int jtype = type_map[j];
      sna_desc.compute_duidrj(sna_desc.rij[jj], sna_desc.wj[jj],
                                  sna_desc.rcutij[jj], jj, sna_desc.element[jj]);
      sna_desc.compute_dbidrj();
      // Accumulate -dBi/dRi, -dBi/dRj
      //app_debug() << "SNAJastrow::compute_d_dr_bispectrum after derivative of descriptor build"<< std::endl;

      int yoffset = ncoeff;
      int zoffset = 2*ncoeff;
      for (int icoeff = 0; icoeff < ncoeff; icoeff++) {
        //std::cout<< "icoeff is " << icoeff <<std::endl;
        a_snad[iat][typeoffset + icoeff]           += sna_desc.dblist[icoeff][0];
        a_snad[iat][typeoffset + icoeff + yoffset] += sna_desc.dblist[icoeff][1];
        a_snad[iat][typeoffset + icoeff + zoffset] += sna_desc.dblist[icoeff][2];
        a_snad[j][typeoffset + icoeff]           -= sna_desc.dblist[icoeff][0];
        a_snad[j][typeoffset + icoeff + yoffset] -= sna_desc.dblist[icoeff][1];
        a_snad[j][typeoffset + icoeff + zoffset] -= sna_desc.dblist[icoeff][2];
      }
    }// end neighbor loop
  }



  void SNAJastrow::evaluateDerivRatios(const VirtualParticleSet& VP,const opt_variables_type& optvars, std::vector<ValueType>& ratios, Matrix<ValueType>& dratios){
    evaluateRatios(VP,ratios);
    Vector<RealType> dlogpsi_nlpp_ref;
    dlogpsi_nlpp_ref.resize(myVars.size());
    Vector<RealType> dlogpsi_nlpp_virt;
    dlogpsi_nlpp_virt.resize(myVars.size());

     for (int k = 0; k < myVars.size(); k++){ // local index
         int k_global = myVars.where(k);
         int ntype = int(k/ncoeff); 
         int coeff = k%ncoeff; 
         //calculate reference deriv
         for (int i = 0; i < Nelec+Nions; i++){
           update_sna_rij(VP.getRefPS(), i, false); 
           compute_bispectrum(i);
         }
         //app_debug() << "SNAJastrow::evaluateDerivRatios after update bispectrum ref" << std::endl;
         //linear derivs for ref
         if (ntype < VP.getRefPS().groups()){ // if 0 or 1 it is up or down electrons
            for (int iat = VP.getRefPS().first(ntype); iat < VP.getRefPS().last(ntype); iat++)
               dlogpsi_nlpp_ref[k] += -sna[iat][coeff];
         }
         else{
          int internal_ntype = ntype-VP.getRefPS().groups();
          for (int iat = Ions.first(internal_ntype); iat < Ions.last(internal_ntype); iat++)
               dlogpsi_nlpp_ref[k] += -sna[Nelec+iat][coeff];
         }
         //app_debug() << "SNAJastrow::evaluateDerivRatios after dlogpsi_nlpp_ref calc" << std::endl;


         //* perform sampling of positions. *//
         for (int r = 0; r < ratios.size(); r++){
            // update the descriptors according to this new particle moving
            for (int i = 0 ; i < Nelec+Nions; i++){
                update_sna_rij_vp(VP,i,r);
                compute_bispectrum(i);
            }
           // if the coefficient belongs to an electron group
           if (ntype < VP.getRefPS().groups()){
              for (int iat = VP.getRefPS().first(ntype); iat < VP.getRefPS().last(ntype); iat++) 
                  dlogpsi_nlpp_virt[k] += -sna[iat][coeff];
           }
          // if the coefficient belongs to an ion group
           else{
              int internal_ntype = ntype-VP.getRefPS().groups();
              for (int iat = Ions.first(internal_ntype); iat < Ions.last(internal_ntype); iat++)
                 dlogpsi_nlpp_virt[k] += -sna[Nelec+iat][coeff];
             //app_debug() << "SNAJastrow::evaluateDerivRatios after dlogpsi_nlpp_virt update ion" << std::endl;
           }
           dratios[r][k_global] =  dlogpsi_nlpp_virt[k] - dlogpsi_nlpp_ref[k];
           //app_debug() << "SNAJastrow::evaluateDerivRatios after dratio update " << std::endl;
         } //end ratio loop
     } // end loop over internal coeffss
  }

  void SNAJastrow::evaluateRatios(const VirtualParticleSet& VP, std::vector<ValueType>& ratios){
    app_debug() << "inside evaluateRatios" << std::endl;
    ScopedTimer local_timer(timers_.eval_ratio_timer);
    double Eold, Enew;
    calculate_ESNA(VP.getRefPS(), snap_beta, Eold, false);
    for (int r = 0; r < ratios.size(); r++){
       for (int i = 0 ; i < Nelec+Nions; i++){
         update_sna_rij_vp(VP, i, r); 
         compute_bispectrum(i);
       }
        //calculate Enew    double esnap_all=0;
      double esnap_elec=0;
      double esnap_ion=0;
      double esnap_all=0;
      double bispectrum_val;
      // calculate electron contribution
      // the global array is summed over groups of atoms of the same type. thus we just need to sum over groups.
      for (int ig = 0; ig < VP.getRefPS().groups(); ig++) {
        for (int iat = VP.getRefPS().first(ig); iat < VP.getRefPS().last(ig); iat++) { // loop over elements in each group
          //app_debug() << "elec index is" << iat <<std::endl; 
          for (int k = 0; k < ncoeff; k++){
            bispectrum_val = sna[iat][k]; //block of bispectrum + current component to add.
            esnap_elec += snap_beta[ig][k] * bispectrum_val;
            //app_debug()<<"snap coefff for group "<<ig<< " coeff "<<k << " is " <<coeff[ig][k] <<std::endl;
          }
        }
      }
      esnap_all += esnap_elec;

      for (int ig = 0; ig < Ions.groups(); ig++) {
        for (int iat = Ions.first(ig); iat < Ions.last(ig); iat++) { // loop over elements in each group
          for (int k =0; k < ncoeff; k++){
            bispectrum_val = sna[Nelec+iat][k];
            esnap_ion += snap_beta[VP.getRefPS().groups()+ig][k] * bispectrum_val;
            //app_debug()<< "snap coeff for group "<< P.groups()+ig << " coeff " << k << " is " << coeff[P.groups()+ig][k] <<std::endl;
          }
        }
      }
      esnap_all += esnap_ion;

      Enew = -esnap_all;
      //store ratio
      ratios[r] = std::exp(static_cast<ValueType>(Enew-Eold));
    }
    return;
  }



  /////////////////////////////////// MC Related functions /////////
  void SNAJastrow::acceptMove(ParticleSet& P, int iat, bool safe_to_delay){
    app_debug() << "inside accept" << std::endl;
    for (int par = 0; par < Nelec+Nions; par++) 
        snad[par].assign(snad[par].size(),0.0);
    for (int par = 0; par < Nions + Nelec; par++){
      update_sna_rij(P, par, false);  
      for (int i = 0; i < num_neigh; i++)
        for (int j = 0; j < 3; j++)
          sna_desc.rij[i][j]*=-1.0;
      compute_d_dr_bispectrum(par, snad);
    }
    double esnap;
    calculate_ESNA(P, snap_beta, esnap, false);
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
    double Enew, Eold;
    int row, col;
    std::vector<std::vector<double>>temp_snad(snad);
    for (int par = 0; par < Nions + Nelec; par++)
      temp_snad[par].assign(temp_snad[par].size(),0.0);

    calculate_ESNA(P, snap_beta, Enew, true);
    //TODO: add update to snad here.
    for (int par = 0; par < Nions + Nelec; par++){
      update_sna_rij(P, par, true);  
      for (int i = 0; i < num_neigh; i++)
        for (int j = 0; j < 3; j++)
          sna_desc.rij[i][j]*=-1.0;
      compute_d_dr_bispectrum(par,temp_snad);
    }
    for (int dim = 0; dim < 3; dim++)
      for (int k = 0; k < ncoeff ; k++)
        for (int n = 0; n < ntypes; n++){
          int col = (n*(3*ncoeff)) + (dim*ncoeff)+k;
          grad_iat[dim] += snap_beta[n][k]*temp_snad[iat][col];
        }
    calculate_ESNA( P, snap_beta, Eold, false);
    SNAJastrow::PsiValue ratio = std::exp(static_cast<SNAJastrow::PsiValue>(Enew-Eold));
    return ratio;
  }

  SNAJastrow::PsiValue SNAJastrow::ratio(ParticleSet& P, int iat){
    app_debug() << "inside ratio" << std::endl;
    double Enew, Eold;
    calculate_ESNA( P, snap_beta, Enew, true);
    calculate_ESNA( P, snap_beta, Eold, false);
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

