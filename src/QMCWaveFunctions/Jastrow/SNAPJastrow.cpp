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
    MPI_Comm_split(MPI_COMM_WORLD,me,0,&comm_lammps);
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
    for (int i=0; i < lmp->atom->ntypes; i++){
      for (int k = 0; k < ncoeff;k++){
        std::stringstream name;
        name << "snap_coeff_" << i;
        name << "_"  << k ;
        myVars.insert(name.str(), snap_beta[i][k], true);
      }
    }
    resizeWFOptVectors();
    grad_u.resize(Nelec);
    lap_u.resize(Nelec);
}

SNAPJastrow::~SNAPJastrow(){
  delete lmp;

}

void SNAPJastrow::set_coefficients(std::vector<double> id_coeffs, int id){
  if (id_coeffs.size() != ncoeff){
    app_log() << "WARNING: number of coeffs less than coeffs/particle_type" <<std::endl; 
    app_log() << "coeffs/particle_type: " << ncoeff << " , coeffss read in: " << id_coeffs.size() << std::endl;
  } 
  app_debug() << "in set coefficients" <<std::endl;
  for (int i=0; i < id_coeffs.size(); i++){
    snap_beta[id][i] = id_coeffs[i];
    myVars[(id*ncoeff) + i] = snap_beta[id][i];
  }

}

LAMMPS_NS::LAMMPS* SNAPJastrow::initialize_lammps(const ParticleSet& els, double rcut){
    ScopedTimer local_timer(timers_.init_lammps_timer);
    const char *lmpargv[] {"liblammps","-log","lammps.out","-screen","lammps_screen.out"};
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
    // add atoms to the groups
    this_lmp->input->one("group e_u type 1");
    this_lmp->input->one("group e_d type 2");
    this_lmp->input->one("group i type 3");
    this_lmp->input->one("group elecs type 1 2");
    this_lmp->input->one("group all type 1 2 3");
    this_lmp->input->one("mass 1 .95");
    this_lmp->input->one("mass 2 .95");
    this_lmp->input->one("mass 3 1");
    for (int ig = 0; ig < els.groups(); ig++) { // loop over groups
      for (int iat = els.first(ig); iat < els.last(ig); iat++) { // loop over elements in each group
        temp_command = std::string("create_atoms ") + std::to_string(ig+1) + " single " + std::to_string((els.R[iat][0]+.1)/bohr_over_ang) + "  " + std::to_string((els.R[iat][1]+.01*iat)/bohr_over_ang)  + " " + std::to_string((els.R[iat][2]+.1*iat+.01)/bohr_over_ang)+ " units box";  
        this_lmp->input->one(temp_command);
      }
    }
      const SpeciesSet& tspecies(Ions.getSpeciesSet());
      for (int ig = 0; ig < Ions.groups(); ig++) { // loop over groups
          temp_command = std::string("group ions_"+ std::to_string(ig)  + " type " + std::to_string(els.groups()+ig+1));
          temp_command = std::string("mass "+ std::to_string(els.groups()+ig+1)) + " 1.00";
          this_lmp->input->one(temp_command);
          for (int iat = Ions.first(ig); iat < Ions.last(ig); iat++) { // loop over elements in each group
            temp_command = std::string("create_atoms "  + std::to_string(els.groups()+ig+1) + " single ") + std::to_string(Ions.R[iat][0]/bohr_over_ang) + "  " + std::to_string(Ions.R[iat][1]/bohr_over_ang)  + " " + std::to_string(Ions.R[iat][2]/bohr_over_ang) + " units box";  
            this_lmp->input->one(temp_command);
          }
      }
      temp_command = std::string("variable twojmax equal ") + std::to_string(twojmax);
      this_lmp->input->one(temp_command);
      this_lmp->input->one("variable 	rcutfac equal 1.0");
      this_lmp->input->one("variable 	rfac0 equal 0.99363");
      //setting rcut to be the same for each type, though may be interesting to try to automate this
      temp_command = std::string("variable rad_type_1 equal ") + std::to_string(rcut/bohr_over_ang);
      this_lmp->input->one(temp_command);
      temp_command = std::string("variable rad_type_2 equal ") + std::to_string(rcut/bohr_over_ang);
      this_lmp->input->one(temp_command);
      temp_command = std::string("variable rad_type_3 equal ") + std::to_string(rcut/bohr_over_ang);
      this_lmp->input->one(temp_command);
      this_lmp->input->one("variable	wj1 equal 1.0");
      this_lmp->input->one("variable	wj2 equal 1.0");
      this_lmp->input->one("variable	wj3 equal 1");
      this_lmp->input->one("variable	quadratic equal 0");
      this_lmp->input->one("variable	bzero equal 0");
      this_lmp->input->one("variable	switchflag equal 0");
      this_lmp->input->one("variable snap_options string \"${rcutfac} ${rfac0} ${twojmax} ${rad_type_1} ${rad_type_2} ${rad_type_3} ${wj1} ${wj2} ${wj3} quadraticflag ${quadratic} bzeroflag ${bzero} switchflag ${switchflag}\"");

    //snap needs some reference pair potential, but doesn't effect parts we are using. 

      temp_command = std::string("pair_style zero ") + std::to_string(rcut*2/bohr_over_ang);
      this_lmp->input->one(temp_command);
      this_lmp->input->one("pair_coeff * *");
      //TODO: generalize with loop over atom types
      this_lmp->input->one("compute sna_global all snap ${snap_options}"); 
      this_lmp->input->one("thermo 100");
      this_lmp->input->one("thermo_style   custom  c_sna_global[1][11] c_sna_global[2][1]");
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

double SNAPJastrow::FD_Lap(const ParticleSet& P,int iat, int dim, int coeff, int ntype, const std::vector<std::vector<double>> coeffs, bool bispectrum_only){
  int row = (iat*3)+dim + 1;
  double G_finite_diff_forward;
  double G_finite_diff_back;
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
  G_finite_diff_back = this_coeff * sna_global->array[row][(ntype*(ncoeff-1))+coeff-1] * hartree_over_ev/bohr_over_ang;
  //fill L
  double finite_diff_lap = (G_finite_diff_forward - G_finite_diff_back)/(2*dist_delta);
  
  // return coordinates to original
  lmp->atom->x[iat][dim] = r0;
  sna_global->compute_array();
  return finite_diff_lap;
}


 SNAPJastrow::LogValue SNAPJastrow::evaluateGL(const ParticleSet& P,
                          ParticleSet::ParticleGradient& G,
                          ParticleSet::ParticleLaplacian& L,
                          bool fromscratch){
                          // std::cout << "in evaluateGL" <<std::endl;
                          return evaluateLog(P,G,L);
  }
 
 
 void SNAPJastrow::computeGL(const ParticleSet& P){
    ScopedTimer local_timer(timers_.eval_gl_timer);
    double grad_val;
    for (int ig = 0; ig < P.groups(); ig++) {
      for (int iel = P.first(ig); iel < P.last(ig); iel++){ // loop over elements in each group
        grad_u[iel] = 0;
        lap_u[iel] = 0;
        for (int n=0; n< lmp->atom->ntypes; n++){
          for (int k =1; k < ncoeff; k ++){
            //app_debug() << "snap beta  at " << n << " " << k << " is " << snap_beta[n][k] << std::endl;
            for (int dim = 0; dim < OHMMS_DIM; dim++){
              int row = (iel*3)+dim + 1;
              int col = (n*(ncoeff-1))+k-1;
              grad_val = sna_global->array[row][col];
              //app_debug() << "grad val is " << grad_val << std::endl;
              grad_u[iel][dim] += snap_beta[n][k]*grad_val*hartree_over_ev/bohr_over_ang;
              lap_u[iel] += FD_Lap(P, iel, dim, k, n, snap_beta, false); 
            }
          }
        }
        //app_debug() << "computeGL Gradient for snap is"  << grad_u[iel] << std::endl;
        //app_debug() << "computeGL laplacian for snap is " << lap_u[iel] << std::endl;
      }
    }
    return;
   }



  SNAPJastrow::LogValue SNAPJastrow::evaluateLog(const ParticleSet& P,
                                    ParticleSet::ParticleGradient& G,
                                    ParticleSet::ParticleLaplacian& L){
    ScopedTimer local_timer(timers_.eval_log_timer);
    for (int i = 0; i < Nelec; i++){
      update_lmp_pos(P,lmp,i,false);
    }
    sna_global->compute_array();
    double esnap;
    calculate_ESNAP(P, sna_global, snap_beta, esnap);
    computeGL(P);
    for (int iel = 0; iel < Nelec; iel++){
      G[iel] += grad_u[iel];
      L[iel] += lap_u[iel];
    }
    log_value_ = static_cast<SNAPJastrow::LogValue>(esnap);
            
    return log_value_;
  }

    SNAPJastrow::GradType SNAPJastrow::evalGrad(ParticleSet& P, int iat){
    for (int i = 0; i < Nelec; i++){
      update_lmp_pos(P, lmp, i, false); //nothing in documentation suggests this is for a proposed move?
    }
    sna_global->compute_array();
    GradType grad_iat;
    for (int dim=0; dim < OHMMS_DIM; dim++){
     int row =(3*iat)+dim+1;
     for (int k = 1; k < ncoeff ; k++){
       for (int n = 0; n < lmp->atom->ntypes; n++){
         int col = (n*(ncoeff-1))+k-1; // lmps isn't aware of beta_0
         grad_iat[dim] += snap_beta[n][k]*sna_global->array[row][col]*hartree_over_ev/bohr_over_ang;
       }
     }
    }
    for (int i = 0; i < Nelec; i++){
      update_lmp_pos(P, lmp, i, false);
    }
    sna_global->compute_array();
    return grad_iat;
    }

    


  void SNAPJastrow::evaluateDerivatives(ParticleSet& P, const opt_variables_type& optvars, Vector<ValueType>& dlogpsi, Vector<ValueType>& dhpsioverpsi)
  {
    ScopedTimer local_timer(timers_.eval_wf_grad_timer);
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
          dhpsioverpsi[kk] = -RealType(0.5) * RealType(Sum(lapLogPsi[k]));
          for (int i = 0; i < Nelec; i++){
            dhpsioverpsi[kk] -= RealType(dot(P.G[i], gradLogPsi[k][i]));
          }
        }
      }
    }
  }



  void SNAPJastrow::evaluateDerivativesWF(ParticleSet& P, const opt_variables_type& optvars, Vector<ValueType>& dlogpsi)
    {
    bool recalculate(false);
    resizeWFOptVectors();
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
        for (int p = 0; p < NumVars; p++){
          gradLogPsi[p] = 0.0;
          lapLogPsi[p] = 0.0;
        }
        dLogPsi = 0.0;
      
        for (int k = 0; k < myVars.size(); k++){
          int kk = myVars.where(k);
          if (kk < 0)
            continue;
          if (rcsingles[k]){
            if (snap_type=="linear"){
              evaluate_linear_derivs(P, k);
            } else if (snap_type=="quadratic"){
              evaluate_fd_derivs(P, k);
            }
            dlogpsi[kk] = ValueType(dLogPsi[k]);
          }
        } 
      }
    }

  void SNAPJastrow::evaluate_linear_derivs(ParticleSet& P, int coeff_idx){
    int ntype = int(coeff_idx/ncoeff); // this coeff is apart of snap with el as central atom. This will also work with multple species of same type as coeeffs/ncoeff will be the type 
    int coeff = coeff_idx%ncoeff; //which coeff of this type are we on.
    if (coeff != 0){ // anything but beta 0 is just tthe bispectrum component
     dLogPsi[coeff_idx] = -sna_global->array[0][(ntype*(ncoeff-1))+coeff-1]*hartree_over_ev;// dlogpsi will be bispectrum component  
     for (int iel =0; iel <Nelec; iel++){
       for (int dim = 0; dim < OHMMS_DIM; dim++){ // loop over dim to get grad vec.
       // sign on this is that gradient is -= in qmcpack, but we have the snap derivatives stored as force, leading to +=
        gradLogPsi[coeff_idx][iel][dim] += sna_global->array[(iel*OHMMS_DIM)+dim+1][(ntype*(ncoeff-1))+coeff-1]*hartree_over_ev/bohr_over_ang;
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
    
    
   /*
     calculates esnap based on a set of coefficients manually in qmcpack
  used to see impact of small change in coefficients on snap energy (needed to calculated d E/d beta)
  without having to internally change the lammps object.
  */
  void SNAPJastrow::calculate_ESNAP(const ParticleSet& P, LAMMPS_NS::ComputeSnap* snap_global, const std::vector<std::vector<double>> coeff, double& new_u){
    ScopedTimer local_timer(timers_.eval_esnap_timer);
    double esnap_all=0;
    double esnap_elec=0;
    double esnap_ion=0;
    double bispectrum_val;
    // calculate electron contribution
    // the global array is summed over groups of atoms of the same type. thus we just need to sum over groups.
    for (int ig = 0; ig < P.groups(); ig++) {
      esnap_elec += coeff[ig][0]*P.groupsize(ig); // beta= contribiution
      for (int k = 1; k < ncoeff; k++){
        bispectrum_val = snap_global->array[0][(ig*(ncoeff-1)) + k-1]; //block of bispectrum + current component to add.
        esnap_elec += coeff[ig][k] * bispectrum_val*hartree_over_ev;
      }
    }
    esnap_all += esnap_elec;
    for (int ig = 0; ig < Ions.groups(); ig++) {
      esnap_ion += coeff[P.groups()+ig][0]*Ions.groupsize(ig);
      for (int k = 1; k < ncoeff; k++){
        bispectrum_val = snap_global->array[0][((P.groups()+ig)*(ncoeff-1)) + k-1];
        esnap_ion += coeff[P.groups()+ig][k] * bispectrum_val*hartree_over_ev;
      }
    }
    esnap_all += esnap_ion;
    new_u = -esnap_all;
    return;
  }


  void SNAPJastrow::evaluateDerivRatios(const VirtualParticleSet& VP,const opt_variables_type& optvars, std::vector<ValueType>& ratios, Matrix<ValueType>& dratios){
    evaluateRatios(VP,ratios);
    Vector<RealType> dlogpsi_nlpp_ref;
    dlogpsi_nlpp_ref.resize(myVars.size());
    Vector<RealType> dlogpsi_nlpp_virt;
    dlogpsi_nlpp_virt.resize(myVars.size());

    bool recalculate(false);
    std::vector<bool> rcsingles(myVars.size(), false);
    for (int k = 0; k < myVars.size(); ++k){
      int kk = myVars.where(k);
      if (kk < 0)
        continue;
      if (optvars.recompute(kk))
        recalculate = true;
      rcsingles[k] = true;
    }

    if (recalculate){
      for (int k = 0; k < myVars.size(); k++){ // local index
        int kk = myVars.where(k); //global index
        if (kk < 0)
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
              dlogpsi_nlpp_ref[k] = -sna_global->array[0][(ntype*(ncoeff-1))+coeff-1]*hartree_over_ev;// dlogpsi will be bispectrum component 
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
            // now for virutal move
            for (int dim= 0; dim < OHMMS_DIM; dim ++){
            // manually update position of ref particle to k position.
              lmp->atom->x[VP.refPtcl][dim] = VP.R[r][dim]/bohr_over_ang;
            }
            sna_global->compute_array();
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
            dratios[r][kk] +=  dlogpsi_nlpp_virt[k] - dlogpsi_nlpp_ref[k];
          } //end ratio loop
        }// end rcsingles
      } // end loop over internal coeffss
    } // end recalculate
    for (int dim= 0; dim < OHMMS_DIM; dim ++){
      // manually update position of ref particle to k position.
      lmp->atom->x[VP.refPtcl][dim] = VP.getRefPS().R[VP.refPtcl][dim]/bohr_over_ang;
    }
    sna_global->compute_array();
  }

  void SNAPJastrow::evaluateRatios(const VirtualParticleSet& VP, std::vector<ValueType>& ratios){
    ScopedTimer local_timer(timers_.eval_ratio_timer);
    double Eold;
    for (int i = 0 ; i < Nelec; i ++){
      update_lmp_pos(VP.getRefPS(), lmp, i, false);
    }
    sna_global->compute_array();
    calculate_ESNAP(VP.getRefPS(),sna_global, snap_beta, Eold);
    for (int r = 0; r < ratios.size(); r++){
     for (int dim= 0; dim < OHMMS_DIM; dim ++){
       // manually update position of ref particle to k position.
       lmp->atom->x[VP.refPtcl][dim] = VP.R[r][dim]/bohr_over_ang;
      }
      sna_global->compute_array();
      //calculate Enew
      double Enew;
      calculate_ESNAP(VP.getRefPS(),sna_global, snap_beta, Enew);
      //store ratio
      ratios[r] = std::exp(static_cast<ValueType>(Enew-Eold));
    }
    for (int dim= 0; dim < OHMMS_DIM; dim ++){
      lmp->atom->x[VP.refPtcl][dim] = VP.getRefPS().R[VP.refPtcl][dim]/bohr_over_ang;
    }
   sna_global->compute_array();
   return;
  }



  /////////////////////////////////// MC Related functions /////////
  void SNAPJastrow::acceptMove(ParticleSet& P, int iat, bool safe_to_delay){
    for (int i = 0 ;i <Nelec; i ++){
     update_lmp_pos(P,lmp,i,true);
    }
    sna_global->compute_array();
    double esnap;
    calculate_ESNAP(P,sna_global, snap_beta, esnap);
    grad_u[iat] = 0;
    lap_u[iat] = 0;
    double grad_val;
    for (int dim = 0; dim < OHMMS_DIM; dim++){
      int row = (iat*3)+dim + 1;
      for (int n = 0; n < lmp->atom->ntypes; n++){
        for (int k =1; k < ncoeff; k ++){
          int col = (n*(ncoeff-1))+k-1;
          grad_val = sna_global->array[row][col];
          //app_debug() << "in accpetmove snap beta is" << snap_beta[n][k] <<std::endl;
          grad_u[iat][dim] += snap_beta[n][k]*grad_val*hartree_over_ev/bohr_over_ang;
          lap_u[iat] += FD_Lap(P, iat, dim, k, n, snap_beta, false); 
        }
      }
    }
  }

  void SNAPJastrow::registerData(ParticleSet& P, WFBufferType& buf){
      log_value_ = evaluateLog(P, P.G, P.L);
  }
  SNAPJastrow::LogValue SNAPJastrow::updateBuffer(ParticleSet& P, WFBufferType& buf, bool from_scratch){
      log_value_ = evaluateLog(P,P.G,P.L);
      return log_value_;
  }

  void SNAPJastrow::copyFromBuffer(ParticleSet& P, WFBufferType& buf){
    
  }

  SNAPJastrow::PsiValue SNAPJastrow::ratioGrad(ParticleSet& P, int iat, GradType& grad_iat){
    std::cout<< "we are in snap ratiograd"<<std::endl;
    SNAPJastrow::PsiValue ratio = SNAPJastrow::ratio(P,iat);// update of proposed_lmp and proposed_sna_global happens in this function
    update_lmp_pos(P, lmp, iat, true);
    sna_global->compute_array();
     for (int dim = 0; dim < OHMMS_DIM; dim++){
      for (int k = 1; k < ncoeff ; k++){
        for (int n = 0; n < lmp->atom->ntypes; n++)
          grad_iat[dim] += snap_beta[n][k]*sna_global->array[(3*iat)+dim+1][(n*(ncoeff-1))+k-1]*hartree_over_ev/bohr_over_ang;
      }
     }
    update_lmp_pos(P, lmp, iat, false);
    sna_global->compute_array();
    return ratio;
  }

  SNAPJastrow::PsiValue SNAPJastrow::ratio(ParticleSet& P, int iat){
    for (int i = 0 ;i < Nelec; i ++){
     update_lmp_pos(P, lmp, i, false);
    }
    sna_global->compute_array();
    double Eold;
     calculate_ESNAP(P, sna_global, snap_beta, Eold);
    for (int i = 0 ;i < Nelec; i ++){
     update_lmp_pos(P, lmp, i, true);
    }
    sna_global->compute_array();
    double Enew;
    calculate_ESNAP(P, sna_global, snap_beta, Enew);
    SNAPJastrow::PsiValue ratio = std::exp(static_cast<SNAPJastrow::PsiValue>(Enew-Eold));
    for (int i = 0 ;i <Nelec; i ++){
     update_lmp_pos(P, lmp, i, false);
    }
     sna_global->compute_array();
     return ratio;
  }

void SNAPJastrow::extractOptimizableObjectRefs(UniqueOptObjRefs& opt_obj_refs){opt_obj_refs.push_back(*this);}

void SNAPJastrow::checkInVariablesExclusive(opt_variables_type& active){
  myVars.setIndexDefault(); // I don't actually know what this is doing?
  active.insertFrom(myVars);
}

void SNAPJastrow::checkOutVariables(const opt_variables_type& active ){
    myVars.getIndex(active);
   /* 
    for (int i=0; i < myVars.size(); i++){
      int loc = myVars.where(i);
      if (loc >=0){
        int ntype = int(i/ncoeff);
        int coeff = i%ncoeff;
        snap_beta[ntype][coeff] = myVars[i] = active[loc];
        app_debug() << "checkoutvariables snap beta " << ntype << " " <<coeff <<" is " << snap_beta[ntype][coeff] <<std::endl;
      }
    }
   */
  }

void SNAPJastrow::resetParametersExclusive(const opt_variables_type& active){
  for (int i=0; i < myVars.size(); i++){
    int loc = myVars.where(i);
    if (loc >=0){
     int ntype = int(i/ncoeff);
     int coeff = i%ncoeff; 
     snap_beta[ntype][coeff] = myVars[i] = active[loc];
     app_debug() << "resetParametersExclusive snap beta is " << snap_beta[ntype][coeff] <<std::endl;
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
/**
void SNAPJastrow::setCoefficients(){
const char *var = "beta";
void *snap_beta_pntr;
double** snap_beta
int b;
snap_beta_pntr = static_cast<LAMMPS_NS::PairSNAP*>(lmp->force->pair)->extract(var,b);
snap_beta = static_cast<double**>(snap_beta_pntr);
for (int i = 0; i<ncoeff; i++){
  default_coeff[i] = *snap_beta[i];
  myVars[i] = *snap_beta[i];

}
}
*/

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

