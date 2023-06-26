#include "SNAPJastrow.h"


namespace qmcplusplus
{

SNAPJastrow::SNAPJastrow(const std::string& obj_name,const ParticleSet& ions, ParticleSet& els) 
  : WaveFunctionComponent(obj_name),
    OptimizableObject("snap_" + ions.getName()),
    myTableID(els.addTable(ions)),
    Nelec(els.getTotalNum()),
    Nions(ions.getTotalNum()),
    NIonGroups(ions.groups()),
    Ions(ions)
{
    twojmax = 2;
    int m = (twojmax/2)+1;
    ncoeff = (m*(m+1)*(2*m+1))/6;

    lmp = initialize_lammps(els,false);
    std::cout << "lammps initialized" <<std::endl;
      proposed_lmp = initialize_lammps(els,false);
      sna_global = lmp->modify->get_compute_by_id("sna_global");
      proposed_sna_global = proposed_lmp->modify->get_compute_by_id("sna_global");
      snap_beta = std::vector<std::vector<double>>((Nions+Nelec), std::vector<double>(ncoeff,0.0));
      for (int i=0; i<(Nelec); i++){
        for (int nc = 0; nc < ncoeff;nc++){
          std::stringstream name;
          name << "snap_coeff_" << i;
          name << "_"  << nc ;
          snap_beta[i][nc] = 0.0;
        myVars.insert(name.str(), snap_beta[i][nc], true);
      }
    }
    resizeWFOptVectors(); 
}

SNAPJastrow::~SNAPJastrow(){
  delete lmp;
}

LAMMPS_NS::LAMMPS * SNAPJastrow::initialize_lammps(const ParticleSet& els, bool proposed){
    std::cout << "in initialize_lammps" <<std::endl;
    const char *lmpargv[] {"liblammps","-log","none","-screen","none"};
    int lmpargc = sizeof(lmpargv)/sizeof(const char *);
    LAMMPS_NS::LAMMPS *this_lmp;
    this_lmp = new LAMMPS_NS::LAMMPS(lmpargc, (char **)lmpargv, MPI_COMM_WORLD);
    std::cout << "created lammps instance" <<std::endl;
    this_lmp->input->one("units  metal");
    this_lmp->input->one("atom_style  atomic");
    // TODO: will this be set by a qmc box? probably.
    this_lmp->input->one("boundary  f f f ");
    this_lmp->input->one("neighbor 1.9 bin");
      this_lmp->input->one("neigh_modify every 1 delay 1 check yes ");
      //TODO: what will box region be pased on QMC object. Cell?
      this_lmp->input->one("region	mybox block -50 50 -50 50 -50 50");
      // create a box that will contain the number of species equal to the number of groups.
      std::string temp_command = std::string("create_box ") + std::to_string(3) +  " mybox";
      this_lmp->input->one(temp_command);
      // add atoms to the groups
      double bohr_over_ang = 1.88973;// to convert qmcpack storage of bohr into angstrom for lammps
      this_lmp->input->one("group e_u type 1");
      this_lmp->input->one("group e_d type 2");
      this_lmp->input->one("group i type 3");
      this_lmp->input->one("group elecs type 1 2");
      this_lmp->input->one("group all type 1 2 3");
      this_lmp->input->one("mass 1 .95");
      this_lmp->input->one("mass 2 .95");
      this_lmp->input->one("mass 3 1");
      std::cout << "before ion creation" <<std::endl;
      for (int ig = 0; ig < els.groups(); ig++) { // loop over groups
          // label groups in lamps
          for (int iat = els.first(ig); iat < els.last(ig); iat++) { // loop over elements in each group
            // place atom in boxes according to their group.
            if (proposed==true){
              temp_command = std::string("create_atoms ") + std::to_string(ig+1) + " single " + std::to_string(els.activeR(iat)[0]/bohr_over_ang) + "  " + std::to_string(els.activeR(iat)[1]/bohr_over_ang)  + " " + std::to_string(els.activeR(iat)[2]/bohr_over_ang)+ " units box";  
              this_lmp->input->one(temp_command);
            }
            else{
              temp_command = std::string("create_atoms ") + std::to_string(ig+1) + " single " + std::to_string((els.R[iat][0]+.1)/bohr_over_ang) + "  " + std::to_string((els.R[iat][1]+.01*iat)/bohr_over_ang)  + " " + std::to_string((els.R[iat][2]+.1*iat+.01)/bohr_over_ang)+ " units box";  
              this_lmp->input->one(temp_command);
            }
          }
      }
      for (int ig = 0; ig < Ions.groups(); ig++) { // loop over groups
        for (int iat = Ions.first(ig); iat < Ions.last(ig); iat++) { // loop over elements in each group
          temp_command = std::string("create_atoms 3 single ") + std::to_string(Ions.R[iat][0]/bohr_over_ang) + "  " + std::to_string(Ions.R[iat][1]/bohr_over_ang)  + " " + std::to_string(Ions.R[iat][2]/bohr_over_ang) + " units box";  
          this_lmp->input->one(temp_command);
        }
      }
      std::cout << "after ion creation" <<std::endl;
      temp_command = std::string("variable twojmax equal ") + std::to_string(twojmax);
      this_lmp->input->one(temp_command);
      this_lmp->input->one("variable 	rcutfac equal 1.0");
      this_lmp->input->one("variable 	rfac0 equal 0.99363");
      this_lmp->input->one("variable 	rmin0 equal 0");
      this_lmp->input->one("variable 	radelem1 equal 2.3");
      this_lmp->input->one("variable 	radelem2 equal 2.3");
      this_lmp->input->one("variable 	radelem3 equal 2.0");
      this_lmp->input->one("variable	wj1 equal 1.0");
      this_lmp->input->one("variable	wj2 equal 1.0");
      this_lmp->input->one("variable	wj3 equal 0.96");
      this_lmp->input->one("variable	quadratic equal 0");
      this_lmp->input->one("variable	bzero equal 0");
      this_lmp->input->one("variable	switch equal 0");
      this_lmp->input->one("variable snap_options string \"${rcutfac} ${rfac0} ${twojmax} ${radelem1} ${radelem2} ${radelem3} ${wj1} ${wj2} ${wj3} rmin0 ${rmin0} quadraticflag ${quadratic} bzeroflag ${bzero} switchflag ${switch}\"");

      this_lmp->input->one("pair_style zero ${rcutfac}");
      //TODO: generalize with loop over atom types
      this_lmp->input->one("pair_coeff * *"); 

      this_lmp->input->one("variable 	zblcutinner equal 4");
      this_lmp->input->one("variable 	zblcutouter equal 4.8");
      this_lmp->input->one("variable 	zblz equal 73");
      this_lmp->input->one("pair_style zbl ${zblcutinner} ${zblcutouter}");
      this_lmp->input->one("pair_coeff 	* * ${zblz} ${zblz}");
      this_lmp->input->one("compute sna_global all snap ${snap_options}"); 
      this_lmp->input->one("run            0");


    return this_lmp;
  }

  void SNAPJastrow::update_lmp_pos(const ParticleSet& P, LAMMPS_NS::LAMMPS* lmp_pntr, LAMMPS_NS::Compute* snap_array,int iat, bool proposed){
      void *pos = lmp_pntr->atom->x;
      double **x = static_cast<double **> (pos);
      if (proposed){
      for (int dim = 0; dim < OHMMS_DIM; dim++){
        (*x)[(iat*3)+dim] = P.activeR(iat)[dim]/bohr_over_ang;
      }
      }
      else{
      for (int dim = 0; dim < OHMMS_DIM; dim++){
        (*x)[(iat*3)+dim] = P.R[iat][dim]/bohr_over_ang;
      }
      }
      snap_array->compute_array();
  }

  SNAPJastrow::LogValueType SNAPJastrow::evaluateGL(const ParticleSet& P,
                          ParticleSet::ParticleGradient& G,
                          ParticleSet::ParticleLaplacian& L,
                          bool fromscratch){

      void *pos = lmp->atom->x;
      double **x = static_cast<double **> (pos);
      RealType delta = 0.001; // TODO: find units
      double G_finite_diff_forward;
      double G_finite_diff_back;
      // compute gradient
          // i.e. pull gradient out from lammps.
    for (int ig = 0; ig < P.groups(); ig++) {
          for (int iel = P.first(ig); iel < P.last(ig); iel++){ // loop over elements in each group
              // loop over elecs in each group
              for (int dim = 0; dim < OHMMS_DIM; dim++){
                for (int nv =0; nv <ncoeff; nv ++){
                int row = (iel*3)+dim + 1;
                int col = nv;
                //G[iel][dim] = lammps_object
                G[iel][dim] += snap_beta[iel][col]*static_cast<double>(sna_global->array[row][col])*hartree_over_ev*bohr_over_ang;
                //forward direction
                RealType r0 = P.R[iel][dim]/bohr_over_ang;
                RealType rp   = r0 + (delta/2);
                //update lammps position
                (*x)[(iel*3)+dim] = rp;
                // recalculate bispectrum stuff
                sna_global -> compute_array();
                G_finite_diff_forward = snap_beta[iel][col]*static_cast<double>(sna_global->array[row][col])*hartree_over_ev*bohr_over_ang;
                //gradient
                //backward direction
                RealType rm   = r0 - (delta/2);
                //update lammps position
                (*x)[(3*iel)+dim] = rm;
                sna_global -> compute_array();
                G_finite_diff_back = snap_beta[iel][col]*static_cast<double>(sna_global->array[row][col])*hartree_over_ev*bohr_over_ang;
                // recalculate bispectrum stuff
                double finite_diff_lap = (G_finite_diff_forward - G_finite_diff_back)/delta;
                //fill L
                L[iel] -=  finite_diff_lap;
                // return coordinates to original
                (*x)[iel+dim] = r0;
                sna_global -> compute_array();
              }
            }
          }
    }
    log_value_ = evaluateLog(P,G,L);
    return log_value_;
  }




  SNAPJastrow::LogValueType SNAPJastrow::evaluateLog(const ParticleSet& P,
                                    ParticleSet::ParticleGradient& G,
                                    ParticleSet::ParticleLaplacian& L){
      double ESNAP=0;
      for (int i = 0; i < (Nelec); i++){
        update_lmp_pos(P,lmp,sna_global,i,false);
      }
      calculate_ESNAP(P, sna_global, snap_beta, ESNAP);
      log_value_ = static_cast<SNAPJastrow::LogValueType>(ESNAP);
              
      return log_value_;
  }

    SNAPJastrow::GradType SNAPJastrow::evalGrad(ParticleSet& P, int iat){
     for (int i = 0; i < (Nelec); i++){
        update_lmp_pos(P, proposed_lmp,proposed_sna_global,i,true);
    }
     GradType grad_iat;
     for (int dim=0; dim < OHMMS_DIM; dim++){
      for (int nc = 0; nc < ncoeff ; nc++){
        grad_iat[dim] += snap_beta[iat][nc]*proposed_sna_global->array[(3*iat)+dim+1][nc]*hartree_over_ev*bohr_over_ang;
      }
     }
     return grad_iat;
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
              std::cout << "entering dlogpsi and laplogpsi[k] is " << lapLogPsi[k] <<std::endl;
              std::cout << "entering dlogpsi and gradlogpsi[k] is " << gradLogPsi[k] <<std::endl;
              std::cout << "entering dlogpsi and P.G is " << P.G <<std::endl;
              dhpsioverpsi[kk] = -RealType(0.5) * RealType(Sum(lapLogPsi[k])) - RealType(Dot(P.G, gradLogPsi[k]));
              std::cout << "dhpsioverpsi is " << dhpsioverpsi[kk] <<std::endl;
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
        evaluate_fd_derivs(P);
        for (int k = 0; k < myVars.size(); k++){
        int kk = myVars.where(k);
        if (kk < 0)
          continue;
        if (rcsingles[k]){
          dlogpsi[kk] = ValueType(dLogPsi[kk]);
        }
      }
    }
    }

  void SNAPJastrow::evaluate_fd_derivs(ParticleSet& P){
        std::vector<std::vector<double>> fd_coeff;
        std::vector<std::vector<double>> bd_coeff;
        fd_coeff = std::vector<std::vector<double>>((Nelec + Nions),std::vector<double>(ncoeff,0.0));
        bd_coeff = std::vector<std::vector<double>>((Nelec + Nions),std::vector<double>(ncoeff,0.0));
        RealType fd_u;
        RealType bd_u;
        RealType coeff_delta = 1e-6;
        RealType dist_delta = 1e-4;
        for (int nv = 0; nv < ncoeff; nv++){
          for (int i=0; i< (Nelec); i++){
            int coeff_start = i*ncoeff;
            if (nv !=0){
            // change back so dloppsi is only with respect to one parameter changing.
              fd_coeff[i][nv-1] = snap_beta[i][nv-1] ;
              bd_coeff[i][nv-1] = snap_beta[i][nv-1] ;
            }
            if (nv ==0 & i!=0){
              //if this isn't the first element then we have to change the last coeff of previous element.
              fd_coeff[i-1][ncoeff-1] = snap_beta[i-1][ncoeff-1] ;
              bd_coeff[i-1][ncoeff-1] = snap_beta[i-1][ncoeff-1] ;
            }
            fd_coeff[i][nv] = snap_beta[i][nv] + coeff_delta;
            bd_coeff[i][nv] = snap_beta[i][nv] - coeff_delta;
      
            calculate_ESNAP(P, sna_global, fd_coeff, fd_u);
            calculate_ESNAP(P, sna_global, bd_coeff, bd_u);
            dLogPsi[coeff_start+nv] += (fd_u - bd_u)/(2*coeff_delta); //units handled elsewhere
            calculate_ddc_gradlap_lammps(P, dist_delta, coeff_delta, fd_coeff, bd_coeff,nv);
        }
      }
    }



  void SNAPJastrow::calculate_ddc_gradlap_lammps(ParticleSet& P,RealType dist_delta, RealType coeff_delta, std::vector<std::vector<double>>& fd_coeff, std::vector<std::vector<double>>& bd_coeff, int cur_val){
    /** brute force calculation of grad lammp where lammps object will have coefficients updated and reran.*/ 
    // 1. update coeffs to forward direction
      // double* x = new double [OHMMS_DIM*P.getTotalNum()];
    void *pos = lmp->atom->x;
    double **x = static_cast<double **> (pos);
    double G_finite_diff_forward;
    double G_finite_diff_back;
    //std::vector<SNAPJastrow::GradDerivVec> ddc_grad_forward;
    //std::vector<SNAPJastrow::GradDerivVec> ddc_grad_back;
    //std::vector<SNAPJastrow::GradDerivVec> ddc_grad_cd;
    //ddc_grad_forward.resize(myVars.size(), GradDerivVec(2));
    //ddc_grad_back.resize(myVars.size(), GradDerivVec(2));
    SNAPJastrow::GradDerivVec ddc_grad_forward_val(2);
    SNAPJastrow::GradDerivVec ddc_grad_back_val(2);
    SNAPJastrow::ValueDerivVec ddc_lap_forward_val(2);
    SNAPJastrow::ValueDerivVec ddc_lap_back_val(2);
    //ddc_grad_forward_val.resize(GradDerivVec(2))
    //ddc_grad_back_val.resize(GradDerivVec(2))
    //ddc_lap_forward_val.resize(ValueDerivVec(2))
    //ddc_lap_back_val.resize(ValueDerivVec(2))

    //std::vector<SNAPJastrow::ValueDerivVec> ddc_lap_forward;
    //std::vector<SNAPJastrow::ValueDerivVec> ddc_lap_back;
    //std::vector<SNAPJastrow::ValueDerivVec> ddc_lap_cd;
    //ddc_lap_forward.resize(myVars.size(),SNAPJastrow::ValueDerivVec(Nelec));
    //ddc_lap_back.resize(myVars.size(),SNAPJastrow::ValueDerivVec(Nelec));
    // 3. calculate grad and lap
    for (int nv = 0; nv < ncoeff; nv++){
      int col = nv;
      for (int iel =0; iel < Nelec; iel ++){ 
        int coeff_start = (iel*ncoeff);
        for (int dim = 0; dim < OHMMS_DIM; dim++){
          int row = (iel*3)+dim + 1;
          //G[iel][dim] = lammps_object
          double grad_val = sna_global->array[row][col];
          std::cout << "grad val is " << grad_val << std::endl;
          std::cout << "for d/dc grad val is " << fd_coeff[iel][nv]*grad_val*hartree_over_ev*bohr_over_ang <<std::endl;
          std::cout << "back d/dc grad val is " << bd_coeff[iel][nv]*grad_val*hartree_over_ev*bohr_over_ang <<std::endl;
          ddc_grad_forward_val[iel][dim] += fd_coeff[iel][nv]*grad_val*hartree_over_ev*bohr_over_ang;
          ddc_grad_back_val[iel][dim] += bd_coeff[iel][nv]*grad_val*hartree_over_ev*bohr_over_ang;


          //forward direction
          RealType r0 = P.R[iel][dim];
          RealType rp   = r0 + (dist_delta/2);
          //update lammps position
          //std::cout<< "original  pos is " << (*x)[(iel*3)+dim] <<std::endl;
          (*x)[(iel*3)+dim] = rp;
          //std::cout<< "updated pos is " << (*x)[(iel*3)+dim] <<std::endl;
            // recalculate bispectrum stuff
          sna_global->compute_array();
          G_finite_diff_forward = sna_global->array[row][col];
          std::cout<< "forward grad is  " << G_finite_diff_forward <<std::endl;

          //back direction
          RealType rm = r0 - (dist_delta/2);
          //update lammps position
          (*x)[(iel*3)+dim] = rm;
          // recalculate bispectrum stuff
          sna_global->compute_array();
          G_finite_diff_back = sna_global->array[row][col];
          std::cout<< "back grad is  " << G_finite_diff_back <<std::endl;
          //finit diff d/dc lap
          double finite_diff_lap_forward = fd_coeff[iel][nv]*((G_finite_diff_forward - G_finite_diff_back)/dist_delta);
          double finite_diff_lap_back = bd_coeff[iel][nv]*((G_finite_diff_forward - G_finite_diff_back)/dist_delta);
          //fill L
          ddc_lap_forward_val[iel] +=  finite_diff_lap_forward*hartree_over_ev*bohr_over_ang;
          ddc_lap_back_val[iel] +=  finite_diff_lap_back*hartree_over_ev*bohr_over_ang;
          // return coordinates to original
          (*x)[(iel*3)+dim] = r0;
          sna_global->compute_array();
        } //end dim
      } //end iel
    } //end nv
    lapLogPsi[cur_val] = -(ddc_lap_forward_val-ddc_lap_back_val)/(2*coeff_delta);
    gradLogPsi[cur_val] = (ddc_grad_forward_val - ddc_grad_back_val)/(2*coeff_delta);
  }
    
    
   /*
     calculates esnap based on a set of coefficients manually in qmcpack
  used to see impact of small change in coefficients on snap energy (needed to calculated d E/d beta)
  without having to internally change the lammps object.
  */
  void SNAPJastrow::calculate_ESNAP(const ParticleSet& P, LAMMPS_NS::Compute* snap_global, std::vector<std::vector<double>> coeff, double& new_u ){
    double ESNAP_all = 0;
    double ESNAP_elec=0;


    for (int ig = 0; ig < P.groups(); ig++) {
      for (int iel = P.first(ig); iel < P.last(ig); iel++){ // loop over elements in each group
        ESNAP_elec = 0; //  maybe beta_0
        for (int k = 0; k < ncoeff; k++){
          // determine atom type
          double bispectrum_val = snap_global->array[0][(iel*ncoeff) + k]; //block of bispectrum + current component to add.
          ESNAP_elec += coeff[iel][k] * bispectrum_val;
          //std::cout << "esnap_elec val is" << ESNAP_elec << std::endl;
        }
        ESNAP_all += ESNAP_elec;
      }
    }
    /**
    int bispectrum_block = P.groups(); // number of groups for electrons will be block for ions. 2*ncoeff = start of ions.
    std::cout << bispectrum_block << " is the bispectrum block "<< std::endl;
    //TODO: only written for single ion at the moment.
    for (int s = 0; s <ns; s++){
        for (int k = 0; k < ncoeff; k++){
          RealType bispectrum_val = snap_global->array[0][(bispectrum_block*ncoeff) + k];
          std::cout << "bispec val for ion is " << bispectrum_val << std::endl;
          std::cout << "coeff val for ion is " << coeff[bispectrum_block][k] << std::endl;
          ESNAP_ion += coeff[bispectrum_block][k] * bispectrum_val  ;
          std::cout << "esnap_ion val is" << ESNAP_ion << std::endl;
        }
      ESNAP_all += ESNAP_ion;
     
    }
    */
    new_u = ESNAP_all*hartree_over_ev;
  }

  /////// Functions for optimization /////
  void SNAPJastrow::checkOutVariables(const opt_variables_type& o ){
    myVars.getIndex(o);
    for (int i=0; i < myVars.size(); i++){
      int loc = myVars.where(i);
      if (loc >=0){
       myVars[i] = o[loc];
     }
     for (int ipart=0; ipart<(Nelec); ipart++){
        for (int nc = 0; nc < ncoeff;nc++){
           snap_beta[ipart][nc] = std::real(myVars[(ipart*ncoeff)+nc]);
        } 
      }
    }
  }

  /////////////////////////////////// MC Related functions /////////

  void SNAPJastrow::acceptMove(ParticleSet& P, int iat, bool safe_to_delay){
      // create new lammps object with accepted positions
      update_lmp_pos(P,lmp,sna_global,iat,false);
      evaluate_fd_derivs(P);
      P.G[iat] = evalGrad(P,iat);

      // the G, dlogpsi, gradlogpsi, and laplogpsi need to be updated but only for single particle.
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

  SNAPJastrow::PsiValueType SNAPJastrow::ratioGrad(ParticleSet& P, int iat, GradType& grad_iat){
     SNAPJastrow::PsiValueType ratio = SNAPJastrow::ratio(P,iat);
     for (int dim=0; dim < OHMMS_DIM; dim++){
      for (int nc = 0; nc < ncoeff ; nc++){
        grad_iat[dim] += snap_beta[iat][nc]*proposed_sna_global->array[(3*iat)+dim+1][nc]*hartree_over_ev*bohr_over_ang;
      }
     }
     return ratio;
  }

  SNAPJastrow::PsiValueType SNAPJastrow::ratio(ParticleSet& P, int iat){
     update_lmp_pos(P,proposed_lmp, proposed_sna_global,iat, true);
     double Eold;
     double Enew;
     calculate_ESNAP(P, proposed_sna_global, snap_beta, Enew);
     calculate_ESNAP(P, sna_global, snap_beta, Eold);
     SNAPJastrow::PsiValueType ratio = std::exp(static_cast<SNAPJastrow::PsiValueType>(Eold -Enew));
     std::cout << "ratio is " << ratio <<std::endl;
     return ratio;
  }

  void SNAPJastrow::extractOptimizableObjectRefs(UniqueOptObjRefs& opt_obj_refs){opt_obj_refs.push_back(*this);}

void SNAPJastrow::checkInVariablesExclusive(opt_variables_type& active){
  active.insertFrom(myVars);
}

void SNAPJastrow::resetParametersExclusive(const opt_variables_type& active){
  for (int i=0; i < myVars.size(); i++){
    int loc = myVars.where(i);
    if (loc >=0){
     myVars[i] = active[loc];
     }
    }
  for (int ipart=0; ipart<(Nelec); ipart++){
    for (int nc = 0; nc < ncoeff;nc++){
      snap_beta[ipart][nc] = std::real(myVars[(ipart*ncoeff)+nc]);
    } 
  }
}


bool SNAPJastrow::put(xmlNodePtr cur) {return true;}
/**
void SNAPJastrow::setCoefficients(){
const char *var = "beta";
void *snap_beta_pntr;
double** snap_beta;
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

