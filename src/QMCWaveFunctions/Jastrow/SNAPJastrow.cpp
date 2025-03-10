#include "SNAPJastrow.h"
#include "ResourceCollection.h"


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

SNAPJastrow::SNAPJastrow(const std::string& obj_name,const ParticleSet& ions, ParticleSet& els,const std::string input_snap_type, int input_twojmax, double input_rcut) 
  : WaveFunctionComponent(obj_name),
    OptimizableObject("snap_" + ions.getName()),
    Nions(ions.getTotalNum()),
    Nelec(els.getTotalNum()),
    NIonGroups(ions.groups()),
    myTableID(els.addTable(ions)),
    Ions(ions),
    timers_("SNAPJatrowTimers")

{
    int n,me,nprocs;
    int nprocs_lammps;
    MPI_Comm_rank(MPI_COMM_WORLD,&me);
    MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
    nprocs_lammps = 1;
    if (nprocs_lammps > nprocs) {
    if (me == 0)
      printf("ERROR: LAMMPS cannot use more procs than available\n");
    MPI_Abort(MPI_COMM_WORLD,1);
  }
  // since qmc is embarassingly parallel we can split each rank into its won comm so that this lammps instance is treated right here.
  // so me is the current rank and create a specific comm for this, which willl belong to this lammps instance. 
  // 0 is the key argument since we don't really care about ranks here.
    MPI_Comm_split(MPI_COMM_WORLD,me, 0, &comm_lammps);
    MPI_Comm_rank(comm_lammps,&me);
  
    twojmax = input_twojmax;
    if (twojmax%2==0){
      int m = (twojmax/2)+1;
      ncoeff = (m*(m+1)*(2*m+1))/6;
    }
    else{
      int m = (twojmax+1)/2;
      ncoeff = (m*(m+1)*(2*m))/3;
    }
    ncoeff +=1;
    rcut = input_rcut;
    snap_type = input_snap_type;
    lmp = initialize_lammps(els, rcut);
    sna_global = static_cast<LAMMPS_NS::ComputeSnap*>(lmp->modify->get_compute_by_id("sna_global"));
    snap_beta = std::vector<std::vector<double>>(lmp->atom->ntypes, std::vector<double>(ncoeff,0.0));
    current_bispectrum = std::vector<std::vector<double>>(lmp->atom->ntypes, std::vector<double>(ncoeff-1, 0.0));
    current_bispectrum_gradient = std::vector<std::vector<double>>(3*(Nelec), std::vector<double>(lmp->atom->ntypes*(ncoeff-1),0.0)); // 3*N x ntype*ncoeff
    for (int k = 0; k < ncoeff; k++){
      for (int i=0; i < lmp->atom->ntypes; i++){
        std::stringstream name;
        name << "snap_coeff_" << i;
        name << "_"  << k ;
        myVars.insert(name.str(), snap_beta[i][k], true);
        app_debug() << "atom type "<< i<< "coeff " <<k << " "<< snap_beta[i][k] << std::endl;
      } 
    }
    //fill internal arrays with lammps information
    for (int k = 0; k < ncoeff-1; k++){
      for (int i=0; i < lmp->atom->ntypes; i++){
        current_bispectrum[i][k] = sna_global->array[0][(i*(ncoeff-1))+k];
        for(int d =0; d<3; d++){
          for(int p =0; p< Nelec; p++){
            current_bispectrum_gradient[(3*p)+d][((ncoeff-1)*i)+k] = sna_global->array[(3*p)+d+1][((ncoeff-1)*i)+k];
          }
        }
      }
    }
    resizeWFOptVectors();
    grad_u.resize(Nelec);
    lap_u.resize(Nelec);
}

SNAPJastrow::~SNAPJastrow(){
  delete lmp;

}

void SNAPJastrow::update_stored_snap(){
    app_debug() << "in update_stored_snap" <<std::endl;
  for (int k = 0; k < ncoeff-1; k++){
    for (int t=0; t < lmp->atom->ntypes; t++){
      current_bispectrum[t][k] = sna_global->array[0][(t*(ncoeff-1))+k];
      for(int d =0; d< 3; d++){
        for(int p =0; p< Nelec; p++){
          current_bispectrum_gradient[(3*p)+d][((ncoeff-1)*t)+k] = sna_global->array[(3*p)+d+1][((ncoeff-1)*t)+k];
        }
      }
    }
  }
}

void SNAPJastrow::set_coefficients(std::vector<double> id_coeffs, int id){
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

LAMMPS_NS::LAMMPS* SNAPJastrow::initialize_lammps(const ParticleSet& els, double rcut){
    app_debug() << "in initialize_lammps" <<std::endl;
    ScopedTimer local_timer(timers_.init_lammps_timer);
    const char *lmpargv[] {"liblammps","-log","lammps.out","-screen","none"};
    int lmpargc = sizeof(lmpargv)/sizeof(const char *);
    LAMMPS_NS::LAMMPS *this_lmp;
    this_lmp = new LAMMPS_NS::LAMMPS(lmpargc, (char **)lmpargv, comm_lammps);

    this_lmp->input->one("units  metal");
    this_lmp->input->one("atom_style  atomic");
    // TODO: will this be set by a qmc box? probably.
    this_lmp->input->one("boundary  f f f ");
    this_lmp->input->one("neighbor 1.9 bin");
    this_lmp->input->one("neigh_modify every 1 delay 1 check yes ");
    //TODO: what will box region be pased on QMC object. Cell?
    this_lmp->input->one("region	mybox block -50 50 -50 50 -50 50");
    // create a box that will contain the number of species equal to the number of groups.
    std::string temp_command = std::string("create_box ") + std::to_string(NIonGroups + els.groups()) +  " mybox";
    this_lmp->input->one(temp_command);
    // add electrons
    this_lmp->input->one("group e_u type 1");
    this_lmp->input->one("group e_d type 2");
    this_lmp->input->one("group elecs type 1 2");
    this_lmp->input->one("mass 1 .95");
    this_lmp->input->one("mass 2 .95");
    //Snap related variables
    temp_command = std::string("variable twojmax equal ") + std::to_string(twojmax);
    this_lmp->input->one(temp_command);
    this_lmp->input->one("variable 	rcutfac equal 1.0");
    this_lmp->input->one("variable 	rfac0 equal 0.99363");
    temp_command = std::string("variable rad_type_1 equal ") + std::to_string(rcut/bohr_over_ang);
    this_lmp->input->one(temp_command);
    temp_command = std::string("variable rad_type_2 equal ") + std::to_string(rcut/bohr_over_ang);
    this_lmp->input->one(temp_command);
    this_lmp->input->one("variable	wj1 equal 1.0");
    this_lmp->input->one("variable	wj2 equal 1.0");
    std::string snap_command =  "variable snap_options string \"${rcutfac} ${rfac0} ${twojmax}";
    std::string rad_command = " ${rad_type_1} ${rad_type_2}"; 
    std::string wj_command = " ${wj1} ${wj2}";
    std::string group_ints = "group all type 1 2";
    int group_num;
    for (int ig = 0; ig < els.groups(); ig++) { // loop over groups
      for (int iat = els.first(ig); iat < els.last(ig); iat++) { // loop over elements in each group
        temp_command = std::string("create_atoms ") + std::to_string(ig+1) + " single " + std::to_string((els.R[iat][0]+.1)/bohr_over_ang) + "  " + std::to_string((els.R[iat][1]+.01*iat)/bohr_over_ang)  + " " + std::to_string((els.R[iat][2]+.1*iat+.01)/bohr_over_ang)+ " units box";  
        this_lmp->input->one(temp_command);
      }
    }
    for (int ig = 0; ig < Ions.groups(); ig++) { // loop over groups
        group_num = els.groups() + ig+1;
        group_ints = group_ints + " " + std::to_string(group_num);
        temp_command = std::string("group ions_" + std::to_string(ig)  + " type " + std::to_string(group_num));
        this_lmp->input->one(temp_command);
        temp_command = std::string("mass "+ std::to_string(group_num)) + " 1.00";
        this_lmp->input->one(temp_command);
        temp_command = std::string("variable wj"+std::to_string(group_num) + " equal 1.0");
        this_lmp->input->one(temp_command);
        temp_command = std::string("variable rad_type_" + std::to_string(group_num)+ " equal ") + std::to_string(rcut/bohr_over_ang);
        this_lmp->input->one(temp_command);
        rad_command = rad_command + std::string(" ${rad_type_" + std::to_string(group_num) + "}");
        wj_command = wj_command + std::string(" ${wj"+ std::to_string(group_num)+ "}");
        for (int iat = Ions.first(ig); iat < Ions.last(ig); iat++) { // loop over elements in each group
          temp_command = std::string("create_atoms "  + std::to_string(group_num) + " single ") + std::to_string(Ions.R[iat][0]/bohr_over_ang) + "  " + std::to_string(Ions.R[iat][1]/bohr_over_ang)  + " " + std::to_string(Ions.R[iat][2]/bohr_over_ang) + " units box";  
          this_lmp->input->one(temp_command);
        }
    }
    this_lmp->input->one(group_ints);
    this_lmp->input->one("variable	quadratic equal 0");
    this_lmp->input->one("variable	bzero equal 0");
    this_lmp->input->one("variable	switchflag equal 0");
    snap_command = snap_command + rad_command + wj_command + " quadraticflag ${quadratic} bzeroflag ${bzero} switchflag ${switchflag}\"";
    this_lmp->input->one(snap_command);

    //snap needs some reference pair potential, but doesn't effect parts we are using. 

      temp_command = std::string("pair_style zero ") + std::to_string(rcut*2/bohr_over_ang);
      this_lmp->input->one(temp_command);
      this_lmp->input->one("pair_coeff * *");
      //TODO: generalize with loop over atom types
      this_lmp->input->one("compute sna_global all snap ${snap_options}"); 
      this_lmp->input->one("thermo 100");
      this_lmp->input->one("thermo_style   custom  c_sna_global[1][1] c_sna_global[1][3]  c_sna_global[2][5]");
      this_lmp->input->one("run            0 pre no post no");
    return this_lmp;
  }

  void SNAPJastrow::update_lmp_pos(const ParticleSet& P, LAMMPS_NS::LAMMPS* lmp_pntr, int iat, bool proposed){
      if (proposed){
        for (int dim = 0; dim < OHMMS_DIM; dim++){
          lmp_pntr->atom->x[iat][dim] = P.activeR(iat)[dim]/bohr_over_ang;
        }
      }
      else{
        for (int dim = 0; dim < OHMMS_DIM; dim++){
          lmp_pntr->atom->x[iat][dim] = P.R[iat][dim]/bohr_over_ang;
        }
      }
  }

double SNAPJastrow::FD_Lap(const ParticleSet& P, int iat, int dim, int coeff, int ntype, const std::vector<std::vector<double>> coeffs, bool bispectrum_only){
  for (int e=0 ; e< Nelec; e++)
    update_lmp_pos(P,lmp,e,false);  
  int row = (iat*3) + dim + 1;
  double G_finite_diff_forward;
  double G_finite_diff_back;
  // if taking wave function laplacian, we just need bispectrum component but for electron laplacian we use the snap coefficients.
  double this_coeff = 1.0;
  if (not bispectrum_only){
    this_coeff = coeffs[ntype][coeff];
  }
  RealType r0 = P.R[iat][dim]/bohr_over_ang;
  
  //forward direction
  RealType rp = r0 + (dist_delta/bohr_over_ang);
  lmp->atom->x[iat][dim] = rp;
  sna_global->compute_array();
  G_finite_diff_forward = this_coeff * sna_global->array[row][(ntype*(ncoeff-1))+coeff-1] * hartree_over_ev/bohr_over_ang;
  
  //backward direction
  RealType rm  = r0 - (dist_delta/bohr_over_ang);
  lmp->atom->x[iat][dim] = rm;
  sna_global->compute_array();
  G_finite_diff_back =    this_coeff * sna_global->array[row][(ntype*(ncoeff-1))+coeff-1] * hartree_over_ev/bohr_over_ang;
  //fill L
  double finite_diff_lap = (G_finite_diff_forward - G_finite_diff_back)/(2*dist_delta); 
  
  return finite_diff_lap;
}


 SNAPJastrow::LogValue SNAPJastrow::evaluateGL(const ParticleSet& P,
                          ParticleSet::ParticleGradient& G,
                          ParticleSet::ParticleLaplacian& L,
                          bool fromscratch){
                          return evaluateLog(P,G,L);
  }
 
 void SNAPJastrow::computeGL(const ParticleSet& P){
    app_debug() << "in computeGL" <<std::endl;
    ScopedTimer local_timer(timers_.eval_gl_timer);
    double grad_val;
    int row,col;
    // calculate the gradient for each electron.
    for (int iel = 0; iel < Nelec; iel++) {
        // reset the internal array
        grad_u[iel] = 0;
        lap_u[iel] = 0;
        //loop over types of particles (electrons and ions)
        for (int n=0; n < lmp->atom->ntypes; n++){
          // loop over the components
          for (int k = 1; k < ncoeff; k ++){
            //we wil need gradient in each direction.
            for (int dim = 0; dim < 3; dim++){
              // get gradient row which is 3*N rows 
              row = (iel*3) + dim; 
              // get the derivative of a specific component with regards to the change in r of this particle. 
              col = (n*(ncoeff-1)) + k-1;
              grad_val = current_bispectrum_gradient[row][col];
              //app_debug() << "grad val is " << grad_val << std::endl;
              // in snap global, the force is stored (i.e. -dB/dr) which is the correct sign since we are taking the negative of Esnap expression
              grad_u[iel][dim] += snap_beta[n][k]*grad_val*hartree_over_ev/bohr_over_ang;
              lap_u[iel] += FD_Lap(P, iel, dim, k, n, snap_beta, false);
            } // end dim loop
          }// end k loop
        } //end n loop
      }// end el loop
    return;
   }



 SNAPJastrow::LogValue SNAPJastrow::evaluateLog(const ParticleSet& P,
                                    ParticleSet::ParticleGradient& G,
                                    ParticleSet::ParticleLaplacian& L){
    ScopedTimer local_timer(timers_.eval_log_timer);
    double esnap;
    calculate_ESNAP(P, current_bispectrum, snap_beta, esnap);
    current_esnap=esnap;
    log_value_ = static_cast<SNAPJastrow::LogValue>(esnap);
    computeGL(P);
    for (int iel = 0; iel < Nelec; iel++){
      G[iel] += grad_u[iel];
      L[iel] += lap_u[iel];
    }
    return log_value_;
  }

    SNAPJastrow::GradType SNAPJastrow::evalGrad(ParticleSet& P, int iat){
    app_debug() << "in evalgrad" <<std::endl;
    GradType grad_iat;
    for (int dim=0; dim < 3; dim++){
     int row = (3*iat) + dim;
     for (int k = 1; k < ncoeff ; k++){
       for (int n = 0; n < lmp->atom->ntypes; n++){
         int col = (n*(ncoeff-1))+k-1; // lmps isn't aware of beta_0
         grad_iat[dim] += snap_beta[n][k]*current_bispectrum_gradient[row][col]*hartree_over_ev/bohr_over_ang;
       }
     }
    }
    return grad_iat;
    }


  void SNAPJastrow::evaluateDerivatives(ParticleSet& P, const opt_variables_type& optvars, Vector<ValueType>& dlogpsi, Vector<ValueType>& dhpsioverpsi)
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


  void SNAPJastrow::evaluateDerivativesWF(ParticleSet& P, const opt_variables_type& optvars, Vector<ValueType>& dlogpsi){
    app_debug() << "in evalderivativeWF"<<std::endl;

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
          if (snap_type=="linear"){
            evaluate_linear_derivs(P, k);
          } else if (snap_type=="quadratic"){
            evaluate_fd_derivs(P, k);
          }
          dlogpsi[k_global] = ValueType(dLogPsi[k]);
        }
      } 
    }
  }

  void SNAPJastrow::evaluate_linear_derivs(ParticleSet& P, int coeff_idx){
    int ntype = int(coeff_idx/ncoeff);  
    int coeff = coeff_idx%ncoeff; //which coeff of this type are we on.
    if (coeff != 0){ // anything but beta 0 is the bispectrum component
     dLogPsi[coeff_idx] = -current_bispectrum[ntype][coeff-1]*hartree_over_ev;// dlogpsi will be bispectrum component  
     for (int iel =0; iel < Nelec; iel++){
       for (int dim = 0; dim < OHMMS_DIM; dim++){ // loop over dim to get grad vec.
        gradLogPsi[coeff_idx][iel][dim] += current_bispectrum_gradient[(iel*3)+dim][(ntype*(ncoeff-1))+coeff-1]*hartree_over_ev/bohr_over_ang;
        lapLogPsi[coeff_idx][iel] += FD_Lap(P, iel, dim, coeff, ntype, snap_beta, true);
       }
     }
    }
    else{ // if we are at a beta_0 eval, then derivative is just -N_type for this beta_0
       if (ntype < P.groups()){
         dLogPsi[coeff_idx] = -P.groupsize(ntype); 
       }
        else{
         dLogPsi[coeff_idx] = -Ions.groupsize(ntype-P.groups()); //shift global list to start at ions.
        }   
    }
  }

  void SNAPJastrow::evaluate_fd_derivs(ParticleSet& P, int coeff_idx){
    /*
        ScopedTimer local_timer(timers_.eval_finite_diff_timer);

        std::vector<std::vector<double>> fd_coeff(snap_beta);
        std::vector<std::vector<double>> bd_coeff(snap_beta); 
        RealType fd_u, bd_u;
        int el = int(coeff_idx/ncoeff);
        int coeff = coeff_idx%ncoeff;
        fd_coeff[el][coeff] = snap_beta[el][coeff] + coeff_delta;
        bd_coeff[el][coeff] = snap_beta[el][coeff] - coeff_delta;
        calculate_ESNAP(P, sna_global, fd_coeff, fd_u);
        calculate_ESNAP(P, sna_global, bd_coeff, bd_u);
        dLogPsi[coeff_idx] = (fd_u - bd_u)/(2*coeff_delta); //units handled elsewhere
        if (coeff !=0){
          calculate_ddc_gradlap_lammps(P, fd_coeff, bd_coeff, coeff_idx);
        }
        */
    }

  void SNAPJastrow::calculate_ddc_gradlap_lammps(ParticleSet& P, std::vector<std::vector<double>>& fd_coeff, std::vector<std::vector<double>>& bd_coeff, int cur_val){
    /*
    SNAPJastrow::GradDerivVec ddc_grad_forward_val(Nelec);
    SNAPJastrow::GradDerivVec ddc_grad_back_val(Nelec);
    SNAPJastrow::ValueDerivVec ddc_lap_forward_val(Nelec);
    SNAPJastrow::ValueDerivVec ddc_lap_back_val(Nelec);
    for (int iel = 0; iel < Nelec; iel++){ // for each particle, step down through rows
      for (int dim = 0; dim < OHMMS_DIM; dim++){ // step down through subblock of each 
        int row = (iel*3) + dim + 1; // organized by xyz for each particle.
        for (int n = 0; n < lmp->atom->ntypes; n++){ // there will be terms from all other particles.
          for (int k = 1; k < ncoeff; k++){ // over all of the coeffs
            double grad_val = sna_global->array[row][(n*(ncoeff-1))+k-1];
            ddc_grad_forward_val[iel][dim] += fd_coeff[n][k]*grad_val*hartree_over_ev/bohr_over_ang;
            ddc_grad_back_val[iel][dim] += bd_coeff[n][k]*grad_val*hartree_over_ev/bohr_over_ang;

            ddc_lap_forward_val[iel] += FD_Lap(P, iel, dim, k, n, fd_coeff, false);
            ddc_lap_back_val[iel] += FD_Lap(P, iel, dim, k, n, bd_coeff,  false);
          } //end dim
        } //end ncoeff 
      } //end ntype
    } //end nelec
    lapLogPsi[cur_val] = (ddc_lap_forward_val - ddc_lap_back_val)/(2*coeff_delta);
    gradLogPsi[cur_val] = (ddc_grad_forward_val - ddc_grad_back_val)/(2*coeff_delta);
    */
  }
    
  void SNAPJastrow::calculate_ESNAP(const ParticleSet& P, LAMMPS_NS::ComputeSnap* snap_global, const std::vector<std::vector<double>> coeff, double& new_u){
    ScopedTimer local_timer(timers_.eval_esnap_timer);
    double esnap_all=0;
    double esnap_elec=0;
    double esnap_ion=0;
    double bispectrum_val;
    // calculate electron contribution
    // the global array is summed over groups of atoms of the same type. thus we just need to sum over groups.
    for (int ig = 0; ig < P.groups(); ig++) {
      esnap_elec += coeff[ig][0]*P.groupsize(ig); // beta0 contribiution
      for (int k = 1; k < ncoeff; k++){
        bispectrum_val = snap_global->array[0][(ig*(ncoeff-1)) + k-1]; //block of bispectrum + current component to add.
        esnap_elec += coeff[ig][k] * bispectrum_val*hartree_over_ev;
        app_debug()<<"snap coefff for group "<<ig<< " coeff "<<k << " is " <<coeff[ig][k] <<std::endl;
      }
    }
    esnap_all += esnap_elec;

    for (int ig = 0; ig < Ions.groups(); ig++) {
      esnap_ion += coeff[P.groups()+ig][0]*Ions.groupsize(ig);
      for (int k = 1; k < ncoeff; k++){
        bispectrum_val = snap_global->array[0][((P.groups()+ig)*(ncoeff-1)) + k-1];
        esnap_ion += coeff[P.groups()+ig][k] * bispectrum_val*hartree_over_ev;
        app_debug()<<"snap coeff for group "<< P.groups()+ig << " coeff "<< k << " is " <<coeff[P.groups()+ig][k] <<std::endl;
      }
    }
    esnap_all += esnap_ion;

    new_u = -esnap_all;
    return;
  }


  
  void SNAPJastrow::calculate_ESNAP(const ParticleSet& P, std::vector<std::vector<double>> current_bispectrum, const std::vector<std::vector<double>> coeff, double& new_u){
    ScopedTimer local_timer(timers_.eval_esnap_timer);
    double esnap_all=0;
    double esnap_elec=0;
    double esnap_ion=0;
    double bispectrum_val;
    // calculate electron contribution
    // the global array is summed over groups of atoms of the same type. thus we just need to sum over groups.
    for (int ig = 0; ig < P.groups(); ig++) {
      esnap_elec += coeff[ig][0]*P.groupsize(ig); // beta0 contribiution
      for (int k = 1; k < ncoeff; k++){
        bispectrum_val = current_bispectrum[ig][k-1]; //block of bispectrum + current component to add.
        esnap_elec += coeff[ig][k] * bispectrum_val * hartree_over_ev;
        app_debug()<< "snap coeff for group "<< ig << " coeff "<<k << " is " <<coeff[ig][k] <<std::endl;
      }
    }
    esnap_all += esnap_elec;

    for (int ig = 0; ig < Ions.groups(); ig++) {
      esnap_ion += coeff[P.groups()+ig][0]*Ions.groupsize(ig);
      for (int k = 1; k < ncoeff; k++){
        bispectrum_val = current_bispectrum[P.groups()+ig][k-1];
        esnap_ion += coeff[P.groups()+ig][k] * bispectrum_val*hartree_over_ev;
        app_debug()<< "snap coeff for group "<< P.groups()+ig << " coeff "<< k << " is " << coeff[P.groups()+ig][k] <<std::endl;
      }
    }
    esnap_all += esnap_ion;

    new_u = -esnap_all;
    return;
  }



  void SNAPJastrow::evaluateDerivRatios(const VirtualParticleSet& VP,const opt_variables_type& optvars, std::vector<ValueType>& ratios, Matrix<ValueType>& dratios){
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
           update_lmp_pos(VP.getRefPS(), lmp, i, false); // make sure lmps objects positions are up to date with ref set
         }
         // update descriptor
         sna_global->compute_array();
         //linear derivs for ref
         if (coeff != 0){ // lmps isn't aware of beta 0 so start at 1 for bispectrum in qmcpack, convert for lmps
             dlogpsi_nlpp_ref[k] = -current_bispectrum[ntype][coeff-1]*hartree_over_ev;// dlogpsi will be bispectrum component 
         }
         else{
           if (ntype < VP.getRefPS().groups()){ // if 0 or 1 it is up or down electrons
             dlogpsi_nlpp_ref[k] = -VP.getRefPS().groupsize(ntype);
           }
           else{
             dlogpsi_nlpp_ref[k] = -Ions.groupsize(ntype-VP.getRefPS().groups());
           }
         }  
         for (int r = 0; r < ratios.size(); r++){
           // now for virutal move, update posi6ioy
           for (int dim = 0; dim < 3; dim ++){
             // manually update position of ref particle to k position.
             lmp->atom->x[VP.refPtcl][dim] = VP.R[r][dim]/bohr_over_ang;
           }
           sna_global->compute_array(); // compute descriptor for updated positions
           if (coeff != 0){
             dlogpsi_nlpp_virt[k] = -sna_global->array[0][(ntype*(ncoeff-1))+coeff-1]*hartree_over_ev;// dlogpsi will be bispectrum component  
           }
           else{
             if (ntype < VP.getRefPS().groups()){
              dlogpsi_nlpp_virt[k] = -VP.getRefPS().groupsize(ntype);
             }
            else{
              dlogpsi_nlpp_virt[k] = -Ions.groupsize(ntype-VP.getRefPS().groups());
             }   
           }
           dratios[r][k_global] =  dlogpsi_nlpp_virt[k] - dlogpsi_nlpp_ref[k];
         } //end ratio loop
        }// end rcsingles
     } // end loop over internal coeffss
    } // end recalculate
  }

  void SNAPJastrow::evaluateRatios(const VirtualParticleSet& VP, std::vector<ValueType>& ratios){
    ScopedTimer local_timer(timers_.eval_ratio_timer);
    double Eold, Enew;
    for (int e=0 ; e < Nelec; e++)
      update_lmp_pos(VP.getRefPS(),lmp,e,false);  
    Eold = current_esnap;
    //calculate_ESNAP(VP.getRefPS(), sna_global, snap_beta, Eold);
    for (int r = 0; r < ratios.size(); r++){
     for (int dim= 0; dim < 3; dim++){
       // manually update position of ref particle to k position.
       lmp->atom->x[VP.refPtcl][dim] = VP.R[r][dim]/bohr_over_ang;
      }
      sna_global->compute_array();
      //calculate Enew
      calculate_ESNAP(VP.getRefPS(), sna_global, snap_beta, Enew);
      //store ratio
      ratios[r] = std::exp(static_cast<ValueType>(Enew-Eold));
    }
    return;
  }



  /////////////////////////////////// MC Related functions /////////
  void SNAPJastrow::acceptMove(ParticleSet& P, int iat, bool safe_to_delay){
    for (int e=0 ; e< Nelec; e++)
      update_lmp_pos(P, lmp, e, true);  
    sna_global->compute_array();
    update_stored_snap(); 
    double esnap;
    calculate_ESNAP(P, current_bispectrum, snap_beta, esnap);
    current_esnap=esnap;
    log_value_ = static_cast<SNAPJastrow::LogValue>(esnap);
    computeGL(P);
    /*grad_u[iat] = 0;
    lap_u[iat] = 0;
    int row, col;
    double grad_val;
    for (int dim = 0; dim < 3; dim++){
      row = (iat*3) + dim;
      for (int n = 0; n < lmp->atom->ntypes; n++){
        for (int k = 1; k < ncoeff; k ++){
          col = (n*(ncoeff-1))+k-1;
          grad_val = current_bispectrum_gradient[row][col];
          grad_u[iat][dim] += snap_beta[n][k]*grad_val*hartree_over_ev/bohr_over_ang;
          lap_u[iat] += FD_Lap(P, iat, dim, k, n, snap_beta, false); 
        }
      }
    }*/
  }

  void SNAPJastrow::registerData(ParticleSet& P, WFBufferType& buf){
  }

  SNAPJastrow::LogValue SNAPJastrow::updateBuffer(ParticleSet& P, WFBufferType& buf, bool from_scratch){
      log_value_ = evaluateLog(P, P.G, P.L);
      return log_value_;
  }

  void SNAPJastrow::copyFromBuffer(ParticleSet& P, WFBufferType& buf){
    
  }

  SNAPJastrow::PsiValue SNAPJastrow::ratioGrad(ParticleSet& P, int iat, GradType& grad_iat){
    for (int e=0 ; e< Nelec; e++)
      update_lmp_pos(P, lmp, e, true);  
    sna_global->compute_array();
    double Enew, Eold;
    int row, col;
    calculate_ESNAP(P, sna_global, snap_beta, Enew);
    for (int dim = 0; dim < 3; dim++){
      row = (3*iat) + dim + 1;
      for (int k = 1; k < ncoeff ; k++){
        for (int n = 0; n < lmp->atom->ntypes; n++){
          col = (n*(ncoeff-1))+ k -1;
          grad_iat[dim] += snap_beta[n][k]*sna_global->array[row][col]*hartree_over_ev/bohr_over_ang;
        }
      }
    }
    Eold=current_esnap;
    SNAPJastrow::PsiValue ratio = std::exp(static_cast<SNAPJastrow::PsiValue>(Enew-Eold));
    return ratio;
  }

  SNAPJastrow::PsiValue SNAPJastrow::ratio(ParticleSet& P, int iat){
    double Enew, Eold;
    for (int e=0 ; e< Nelec; e++)
      update_lmp_pos(P, lmp, e, true);  
    sna_global->compute_array();
    calculate_ESNAP(P, sna_global, snap_beta, Enew);
    Eold = current_esnap;
    //calculate the ratio
    SNAPJastrow::PsiValue ratio = std::exp(static_cast<SNAPJastrow::PsiValue>(Enew-Eold));
    return ratio;
  }

void SNAPJastrow::extractOptimizableObjectRefs(UniqueOptObjRefs& opt_obj_refs){opt_obj_refs.push_back(*this);}

void SNAPJastrow::checkInVariablesExclusive(opt_variables_type& active){
  myVars.setIndexDefault(); // I don't actually know what this is doing?
  active.insertFrom(myVars);
}

void SNAPJastrow::checkOutVariables(const opt_variables_type& active ){
    myVars.getIndex(active);
  }

void SNAPJastrow::resetParametersExclusive(const opt_variables_type& active){
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

 
std::unique_ptr<WaveFunctionComponent> SNAPJastrow::makeClone(ParticleSet& tpq) const
{
  auto snap_copy = std::make_unique<SNAPJastrow>(std::string("snap"), Ions, tpq, std::string("linear"), twojmax, rcut);
  snap_copy->snap_beta = snap_beta;
  return snap_copy;
}


bool SNAPJastrow::put(xmlNodePtr cur) {
  
  app_summary() << "     Number of parameters: " << myVars.size() << std::endl;
  for (int i = 0; i <myVars.size(); i++){
    app_summary() << myVars[i] <<std::endl;
  }
  return true;}

void SNAPJastrow::createResource(ResourceCollection& collection) const
{
  collection.addResource(std::make_unique<SNAMultiWalkerMem<RealType>>());
}

void SNAPJastrow::acquireResource(ResourceCollection& collection,
                                       const RefVectorWithLeader<WaveFunctionComponent>& wfc_list) const
{
  auto& wfc_leader          = wfc_list.getCastedLeader<SNAPJastrow>();
  wfc_leader.mw_mem_handle_ = collection.lendResource<SNAMultiWalkerMem<RealType>>();
}

void SNAPJastrow::releaseResource(ResourceCollection& collection,
                                       const RefVectorWithLeader<WaveFunctionComponent>& wfc_list) const
{
  auto& wfc_leader = wfc_list.getCastedLeader<SNAPJastrow>();
  collection.takebackResource(wfc_leader.mw_mem_handle_);
}


}

