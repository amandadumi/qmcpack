#include "SNAPJastrow.h"


namespace qmcplusplus
{

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
  // snce qmc is treated as emparassingle parallel we can split eachh rank into its won comm so that this lammps instance is treated right here.
  // so me is the current rank and create a specific comm for this, which willl belong to this lammps instance. 
  // 0 is the key argument sincce we don't really care about ranks here.
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
    proposed_lmp = initialize_lammps(els, rcut);
    sna_global = static_cast<LAMMPS_NS::ComputeSnap*>(lmp->modify->get_compute_by_id("sna_global"));
    proposed_sna_global = static_cast<LAMMPS_NS::ComputeSnap*>(proposed_lmp->modify->get_compute_by_id("sna_global"));
    snap_beta = std::vector<std::vector<double>>(lmp->atom->ntypes, std::vector<double>(ncoeff,0.01));
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
  // MPI_Comm_free(&comm_lammps);
  delete lmp;
  delete proposed_lmp;
  //delete vp_lmp;
//  MPI_Comm_free(&comm_lammps);

}

void SNAPJastrow::set_coefficients(std::vector<double> id_coeffs,int id){
  if (id_coeffs.size() != ncoeff){
  }
  int kk=0;
  for (int i=0; i < id_coeffs.size(); i++){
    snap_beta[id][i] = id_coeffs[i];
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

  void SNAPJastrow::update_lmp_pos(const ParticleSet& P, LAMMPS_NS::LAMMPS* lmp_pntr, LAMMPS_NS::ComputeSnap* snap_array, int iat, bool proposed){
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
      snap_array->compute_array();
  }

double SNAPJastrow::FD_Lap(const ParticleSet& P,int iat, int dim, int coeff, int ntype, std::vector<std::vector<double>> coeffs, bool bispectrum_only){
  int row = (iat*3)+dim + 1;
  double G_finite_diff_forward;
  double G_finite_diff_back;
  double this_coeff = 1.0;
  if (not bispectrum_only){
    this_coeff = coeffs[ntype][coeff];
  }
  RealType r0 = P.R[iat][dim]/bohr_over_ang;
  
  //forward direction
  RealType rp = r0 + (dist_delta/2);
  lmp->atom->x[iat][dim] = rp;
  sna_global->compute_array();
  G_finite_diff_forward = this_coeff * sna_global->array[row][(ntype*(ncoeff-1))+coeff-1] * hartree_over_ev/bohr_over_ang;
  
  //backward direction
  RealType rm  = r0 - (dist_delta/2);
  lmp->atom->x[iat][dim] = rm;
  sna_global->compute_array();
  G_finite_diff_back = this_coeff * sna_global->array[row][(ntype*(ncoeff-1))+coeff-1] * hartree_over_ev/bohr_over_ang;
  //fill L
  double finite_diff_lap = (G_finite_diff_forward - G_finite_diff_back)/(dist_delta);
  
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
        for (int dim = 0; dim < OHMMS_DIM; dim++){
          int row = (iel*3)+dim + 1;
          for (int n = 0; n < lmp->atom->ntypes; n++){
            for (int k =1; k < ncoeff; k ++){
              int col = (n*(ncoeff-1))+k-1;
              grad_val = sna_global->array[row][col];
              grad_u[iel][dim] += snap_beta[n][k]*grad_val*hartree_over_ev/bohr_over_ang;
              lap_u[iel] += FD_Lap(P, iel, dim, k, n, snap_beta, false)/bohr_over_ang; 
            }
          }
        }
      }
    }
    return;
   }



  SNAPJastrow::LogValue SNAPJastrow::evaluateLog(const ParticleSet& P,
                                    ParticleSet::ParticleGradient& G,
                                    ParticleSet::ParticleLaplacian& L){
    ScopedTimer local_timer(timers_.eval_log_timer);
    for (int i = 0; i < Nelec; i++){
      update_lmp_pos(P,lmp,sna_global,i,false);
    }
    double esnap;
    calculate_ESNAP(P, sna_global, snap_beta, esnap);
    u_val = esnap;
    computeGL(P);
    for (int iel = 0; iel < Nelec; iel++){
      G[iel] += grad_u[iel];
      L[iel] += lap_u[iel];
    }
    log_value_ = static_cast<SNAPJastrow::LogValue>(esnap);
            
    return log_value_;
  }

    SNAPJastrow::GradType SNAPJastrow::evalGrad(ParticleSet& P, int iat){
     update_lmp_pos(P, proposed_lmp, proposed_sna_global, iat, true);
     GradType grad_iat;
     for (int dim=0; dim < OHMMS_DIM; dim++){
      int row =(3*iat)+dim+1;
      for (int k = 1; k < ncoeff ; k++){
        for (int n = 0; n < lmp->atom->ntypes; n++){
          int col = (n*(ncoeff-1))+k-1;
          grad_iat[dim] += snap_beta[n][k]*proposed_sna_global->array[row][col]*hartree_over_ev/bohr_over_ang;
        }
      }
     }
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
        for (int p = 0; p < NumVars; ++p){
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
    int coeff = coeff_idx%ncoeff; //ddd which coeff of this el are we on.
    if (coeff != 0){
     dLogPsi[coeff_idx] = -sna_global->array[0][(ntype*(ncoeff-1))+coeff-1]*hartree_over_ev;// dlogpsi will be bispectrum component  
     for (int iel =0; iel <Nelec; iel++){
       for (int dim = 0; dim < OHMMS_DIM; dim++){ // loop over dim to get grad vec.
        gradLogPsi[coeff_idx][iel][dim] += sna_global->array[(iel*OHMMS_DIM)+dim+1][(ntype*(ncoeff-1))+coeff-1]*hartree_over_ev/bohr_over_ang;
        lapLogPsi[coeff_idx][iel] += FD_Lap(P, iel, dim, coeff, ntype, snap_beta, true)/bohr_over_ang;
       }
     }
    }
    else{
      dLogPsi[coeff_idx] =-1;
    }
  }

  void SNAPJastrow::evaluate_fd_derivs(ParticleSet& P, int coeff_idx){
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
    }



  void SNAPJastrow::calculate_ddc_gradlap_lammps(ParticleSet& P, std::vector<std::vector<double>>& fd_coeff, std::vector<std::vector<double>>& bd_coeff, int cur_val){
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

            ddc_lap_forward_val[iel] += FD_Lap(P, iel, dim, k, n, fd_coeff, false)/bohr_over_ang;
            ddc_lap_back_val[iel] += FD_Lap(P, iel, dim, k, n, bd_coeff,  false)/bohr_over_ang;
          } //end dim
        } //end ncoeff 
      } //end ntype
    } //end nelec
    lapLogPsi[cur_val] = (ddc_lap_forward_val - ddc_lap_back_val)/(2*coeff_delta);
    gradLogPsi[cur_val] = (ddc_grad_forward_val - ddc_grad_back_val)/(2*coeff_delta);
  }
    
    
   /*
     calculates esnap based on a set of coefficients manually in qmcpack
  used to see impact of small change in coefficients on snap energy (needed to calculated d E/d beta)
  without having to internally change the lammps object.
  */
  void SNAPJastrow::calculate_ESNAP(const ParticleSet& P, LAMMPS_NS::ComputeSnap* snap_global, std::vector<std::vector<double>> coeff, double& new_u){
    ScopedTimer local_timer(timers_.eval_esnap_timer);
    double esnap_all=0;
    double esnap_elec=0;
    double esnap_ion=0;
    double bispectrum_val;
    // calculate electron contribution
    // the global array is summed over groups of atoms of the same type. thus here we just need to sum over groups.
    for (int ig = 0; ig < P.groups(); ig++) {
      esnap_elec += coeff[ig][0];
      for (int k = 1; k < ncoeff; k++){
        bispectrum_val = snap_global->array[0][(ig*(ncoeff-1)) + k-1]; //block of bispectrum + current component to add.
        esnap_elec += coeff[ig][k] * bispectrum_val*hartree_over_ev;
      }
    }
    esnap_all += esnap_elec;
    for (int ig = 0; ig < Ions.groups(); ig++) {
      esnap_ion += coeff[P.groups()+ig][0];
      for (int k = 1; k < ncoeff; k++){
        bispectrum_val = snap_global->array[0][((P.groups()+ig)*(ncoeff-1)) + k-1];// will need fixed for more than one ion group
        esnap_ion += coeff[P.groups()+ig][k] * bispectrum_val*hartree_over_ev;// will need fixed for more than one ion group.
      }
    }
    esnap_all += esnap_ion;
    new_u = -esnap_all;
    return;
  }


  void SNAPJastrow::evaluateRatios(const VirtualParticleSet& VP, std::vector<ValueType>& ratios){
    ScopedTimer local_timer(timers_.eval_ratio_timer);
        // create lmp object 
        // done in constructer.
        double Eold;
        calculate_ESNAP(VP.getRefPS(),sna_global, snap_beta, Eold);
        for (int i = 0 ; i < Nelec; i ++){
          update_lmp_pos(VP.getRefPS(), proposed_lmp, proposed_sna_global, i, false);
        }
        for (int r = 0; r < ratios.size(); r++){
          for (int dim= 0; dim < OHMMS_DIM; dim ++){
            // manually update posiition of ref particle to k position.
            proposed_lmp->atom->x[VP.refPtcl][dim] = VP.R[r][dim];
          }
          proposed_sna_global->compute_array();
         //calculate Enew
          double Enew;
          calculate_ESNAP(VP.getRefPS(),proposed_sna_global, snap_beta, Enew);
          //store ratio
          ratios[r] = std::exp(static_cast<ValueType>(Enew-Eold));
        }
        for (int dim= 0; dim < OHMMS_DIM; dim ++){
          proposed_lmp->atom->x[VP.refPtcl][dim] = VP.getRefPS().R[VP.refPtcl][dim];
        }
    return;
  }



  /////////////////////////////////// MC Related functions /////////
  void SNAPJastrow::acceptMove(ParticleSet& P, int iat, bool safe_to_delay){
    update_lmp_pos(P,lmp,sna_global,iat,false);
    double esnap;
    u_val=esnap;
    grad_u[iat] = 0;
    lap_u[iat] = 0;
    double grad_val;
    for (int dim = 0; dim < OHMMS_DIM; dim++){
      int row = (iat*3)+dim + 1;
      for (int n = 0; n < lmp->atom->ntypes; n++){
        for (int k =1; k < ncoeff; k ++){
          int col = (n*(ncoeff-1))+k-1;
          grad_val = sna_global->array[row][col];
          grad_u[iat][dim] += snap_beta[n][k]*grad_val*hartree_over_ev/bohr_over_ang;
          lap_u[iat] += FD_Lap(P, iat, dim, k, n, snap_beta, false)/bohr_over_ang; 
        }
      }
    }
  }

  void SNAPJastrow::registerData(ParticleSet& P, WFBufferType& buf){
      // log_value_ = evaluateLog(P, P.G, P.L);

  }
  SNAPJastrow::LogValue SNAPJastrow::updateBuffer(ParticleSet& P, WFBufferType& buf, bool from_scratch){
      // log_value_ = evaluateLog(P,P.G,P.L);
      return log_value_;
  }

  void SNAPJastrow::copyFromBuffer(ParticleSet& P, WFBufferType& buf){
    
  }

  SNAPJastrow::PsiValue SNAPJastrow::ratioGrad(ParticleSet& P, int iat, GradType& grad_iat){
     SNAPJastrow::PsiValue ratio = SNAPJastrow::ratio(P,iat);
     for (int dim = 0; dim < OHMMS_DIM; dim++){
      for (int k = 1; k < ncoeff ; k++){
        for (int n = 0; n < lmp->atom->ntypes; n++)
          grad_iat[dim] += snap_beta[n][k]*proposed_sna_global->array[(3*iat)+dim+1][(n*(ncoeff-1))+k-1]*hartree_over_ev/bohr_over_ang;
      }
     }
     return ratio;
  }

  SNAPJastrow::PsiValue SNAPJastrow::ratio(ParticleSet& P, int iat){
    update_lmp_pos(P, proposed_lmp, proposed_sna_global, iat, true);
     double Eold;
     calculate_ESNAP(P, sna_global, snap_beta, Eold);
     double Enew;
     calculate_ESNAP(P, proposed_sna_global, snap_beta, Enew);
     SNAPJastrow::PsiValue ratio = std::exp(static_cast<SNAPJastrow::PsiValue>(Enew-Eold));
     return ratio;
  }

  void SNAPJastrow::extractOptimizableObjectRefs(UniqueOptObjRefs& opt_obj_refs){opt_obj_refs.push_back(*this);}

void SNAPJastrow::checkInVariablesExclusive(opt_variables_type& active){
  active.insertFrom(myVars);
}

void SNAPJastrow::checkOutVariables(const opt_variables_type& active ){
    myVars.getIndex(active);
    /*
    for (int i=0; i < myVars.size(); i++){
      int loc = myVars.where(i);
      if (loc >=0){
       myVars[i] = o[loc];
      }
      int ntype = int(i/ncoeff);
      int coeff = i%ncoeff;
      snap_beta[ntype][coeff] = std::real(myVars[(ntype*ncoeff)+coeff]);
    }
    */
  }

void SNAPJastrow::resetParametersExclusive(const opt_variables_type& active){
  for (int i=0; i < myVars.size(); i++){
    int loc = myVars.where(i);
    if (loc >=0){
     myVars[i] = active[loc];
     int ntype = int(i/ncoeff);
     int coeff = i%ncoeff; 
     snap_beta[ntype][coeff] = myVars[i];
     }
    }
  }
  

std::unique_ptr<WaveFunctionComponent> SNAPJastrow::makeClone(ParticleSet& tpq) const
{
  auto snap_copy = std::make_unique<SNAPJastrow>(std::string("snap"),Ions, tpq, std::string("linear"), twojmax, rcut);
  return snap_copy;
}


bool SNAPJastrow::put(xmlNodePtr cur) {
  
  app_summary() << "     Number of parameters: " << myVars.size() << std::endl;
  for (int i = 0; i <myVars.size();i++){
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


}

