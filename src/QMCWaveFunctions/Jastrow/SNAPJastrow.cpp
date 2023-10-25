#include "SNAPJastrow.h"


namespace qmcplusplus
{

SNAPJastrow::SNAPJastrow(const std::string& obj_name,const ParticleSet& ions, ParticleSet& els,int input_twojmax) 
  : WaveFunctionComponent(obj_name),
    OptimizableObject("snap_" + ions.getName()),
    OptimizableObject("snap_" + ions.getName()),
    myTableID(els.addTable(ions)),
    Nelec(els.getTotalNum()),
    Nions(ions.getTotalNum()),
    NIonGroups(ions.groups()),
    Ions(ions)
{

    if (NIonGroups >1){
      app_log() << "WARNING: SNAP jastrow does not currently support more than one ion group" <<std::endl;
    }
    twojmax = input_twojmax;
    int m = (twojmax/2)+1;
    ncoeff = (m*(m+1)*(2*m+1))/6;

    lmp = initialize_lammps(els);
    proposed_lmp = initialize_lammps(els);
    vp_lmp = initialize_lammps(els); //may not need, but lets be extra sure wires don't get crossed
    sna_global = static_cast<LAMMPS_NS::ComputeSnap*>(lmp->modify->get_compute_by_id("sna_global"));
    proposed_sna_global = static_cast<LAMMPS_NS::ComputeSnap*>(proposed_lmp->modify->get_compute_by_id("sna_global"));
    snap_beta = std::vector<std::vector<double>>((Nions+Nelec), std::vector<double>(ncoeff,0.1));
    for (int i=0; i<(Nelec+Nions); i++){
      for (int nc = 0; nc < ncoeff;nc++){
        std::stringstream name;
        name << "snap_coeff_" << i;
        name << "_"  << k ;
        //std::cout<< name.str() <<std::endl;
        myVars.insert(name.str(), snap_beta[i][k], true);
      }
    }
    resizeWFOptVectors(); 
    grad_u.resize(Nelec+Nions);
    lap_u.resize(Nelec+Nions);
    u_val.resize(Nelec+Nions);
}

SNAPJastrow::~SNAPJastrow(){
  delete lmp;
  delete proposed_lmp;
  delete vp_lmp;
}


void SNAPJastrow::set_coefficients(std::vector<double> id_coeffs,int id){
  //std::cout<< id_coeffs.size() << "is the size of coefficients" <<std::endl;
  //std::cout<< ncoeff << "is the number of coeffs" <<std::endl;
  if (id_coeffs.size() != ncoeff){
    app_warning() << " Warning wrong number of coefficents for snap jastrow" << std::endl;
  }
  int kk=0;
  for (int i=0; i < id_coeffs.size(); i++){
    snap_beta[id][i] = id_coeffs[i];
  }

}

LAMMPS_NS::LAMMPS * SNAPJastrow::initialize_lammps(const ParticleSet& els){
    std::cout << "in initialize_lammps" <<std::endl;
    const char *lmpargv[] {"liblammps","-log","lammps.out","-screen","lammps_screen.out"};
    const char *lmpargv[] {"liblammps","-log","lammps.out","-screen","lammps_screen.out"};
    int lmpargc = sizeof(lmpargv)/sizeof(const char *);
    LAMMPS_NS::LAMMPS *this_lmp;
    this_lmp = new LAMMPS_NS::LAMMPS(lmpargc, (char **)lmpargv, MPI_COMM_WORLD);
    std::cout << "created lammps instance" <<std::endl;
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
        temp_command = std::string("create_atoms ") + std::to_string(ig+1) + " single " + std::to_string((els.R[iat][0]+.1)/bohr_over_ang) + "  " + std::to_string((els.R[iat][1]+.01*iat)/bohr_over_ang)  + " " + std::to_string((els.R[iat][2]+.1*iat+.01)/bohr_over_ang)+ " units box";  
        std::cout << temp_command << std::endl;
        this_lmp->input->one(temp_command);
      }
    }
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
        temp_command = std::string("create_atoms ") + std::to_string(ig+1) + " single " + std::to_string((els.R[iat][0]+.1)/bohr_over_ang) + "  " + std::to_string((els.R[iat][1]+.01*iat)/bohr_over_ang)  + " " + std::to_string((els.R[iat][2]+.1*iat+.01)/bohr_over_ang)+ " units box";  
        std::cout << temp_command << std::endl;
        this_lmp->input->one(temp_command);
      }
    }
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
        temp_command = std::string("create_atoms ") + std::to_string(ig+1) + " single " + std::to_string((els.R[iat][0]+.1)/bohr_over_ang) + "  " + std::to_string((els.R[iat][1]+.01*iat)/bohr_over_ang)  + " " + std::to_string((els.R[iat][2]+.1*iat+.01)/bohr_over_ang)+ " units box";  
        //std::cout << temp_command << std::endl;
        this_lmp->input->one(temp_command);
      }
    }
      for (int ig = 0; ig < Ions.groups(); ig++) { // loop over groups
        for (int iat = Ions.first(ig); iat < Ions.last(ig); iat++) { // loop over elements in each group
          temp_command = std::string("create_atoms 3 single ") + std::to_string(Ions.R[iat][0]/bohr_over_ang) + "  " + std::to_string(Ions.R[iat][1]/bohr_over_ang)  + " " + std::to_string(Ions.R[iat][2]/bohr_over_ang) + " units box";  
          std::cout << temp_command << std::endl;
          std::cout << temp_command << std::endl;
          this_lmp->input->one(temp_command);
        }
      }
      std::cout << "after ion creation" <<std::endl;
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
      this_lmp->input->one("variable	bzero equal 1");
      this_lmp->input->one("variable	switch equal 0");
      this_lmp->input->one("variable snap_options string \"${rcutfac} ${rfac0} ${twojmax} ${radelem1} ${radelem2} ${radelem3} ${wj1} ${wj2} ${wj3} rmin0 ${rmin0} quadraticflag ${quadratic} bzeroflag ${bzero} switchflag ${switch}\"");

      this_lmp->input->one("pair_style zero ${rcutfac}");
      //TODO: generalize with loop over atom types
      this_lmp->input->one("pair_coeff * *"); 
      this_lmp->input->one("pair_style zero ${rcutfac}");
      //TODO: generalize with loop over atom types
      this_lmp->input->one("pair_coeff * *"); 

      this_lmp->input->one("variable 	zblcutinner equal 4");
      this_lmp->input->one("variable 	zblcutouter equal 4.8");
      this_lmp->input->one("variable 	zblz equal 73");
      this_lmp->input->one("pair_style zbl ${zblcutinner} ${zblcutouter}");
      this_lmp->input->one("pair_coeff 	* * ${zblz} ${zblz}");
      this_lmp->input->one("compute sna_global all snap ${snap_options}"); 
      this_lmp->input->one("thermo 100");
      this_lmp->input->one("thermo_style   custom  c_sna_global[1][11] c_sna_global[2][1]");
      this_lmp->input->one("thermo 100");
      this_lmp->input->one("thermo_style   custom  c_sna_global[1][11] c_sna_global[2][1]");
      this_lmp->input->one("run            0");
      //std::cout<< "lammps ntypes is " << this_lmp->atom->ntypes << std::endl;

    return this_lmp;
  }

  void SNAPJastrow::update_lmp_pos(const ParticleSet& P, LAMMPS_NS::LAMMPS* lmp_pntr, LAMMPS_NS::ComputeSnap* snap_array,int iat, bool proposed){
  void SNAPJastrow::update_lmp_pos(const ParticleSet& P, LAMMPS_NS::LAMMPS* lmp_pntr, LAMMPS_NS::ComputeSnap* snap_array,int iat, bool proposed){
      if (proposed){
        for (int dim = 0; dim < OHMMS_DIM; dim++){
          lmp_pntr->atom->x[iat][dim] = P.activeR(iat)[dim]/bohr_over_ang;
        }
        for (int dim = 0; dim < OHMMS_DIM; dim++){
          lmp_pntr->atom->x[iat][dim] = P.activeR(iat)[dim]/bohr_over_ang;
        }
      }
      else{
        for (int dim = 0; dim < OHMMS_DIM; dim++){
          lmp_pntr->atom->x[iat][dim] = P.R[iat][dim]/bohr_over_ang;
        }
        for (int dim = 0; dim < OHMMS_DIM; dim++){
          lmp_pntr->atom->x[iat][dim] = P.R[iat][dim]/bohr_over_ang;
        }
      }
      snap_array->compute_array();
  }

double SNAPJastrow::FD_Lap(const ParticleSet& P,int iat, int dim, int row, int coeff, int ntype, std::vector<std::vector<double>> coeffs, double dist_delta){
  //std::cout << "in FD_Lap"<<std::endl;
  double G_finite_diff_forward;
  double G_finite_diff_back;
  //forward direction
  RealType r0 = P.R[iat][dim]/bohr_over_ang;
  RealType rp   = r0 + (dist_delta/2);
  lmp->atom->x[iat][dim] = rp;
  sna_global -> compute_array();
  G_finite_diff_forward = -coeffs[P.GroupID[iat]][coeff]*sna_global->array[row][(ntype*ncoeff)+coeff]*hartree_over_ev/bohr_over_ang;
  
  //backward direction
  RealType rm  = r0 - (dist_delta/2);
  //update lammps position
  lmp->atom->x[iat][dim] = rm;
  sna_global -> compute_array();
  G_finite_diff_back = -coeffs[P.GroupID[iat]][coeff]*sna_global->array[row][(ntype*ncoeff)+coeff]*hartree_over_ev/bohr_over_ang;
  // recalculate bispectrum stuff
  double finite_diff_lap = (G_finite_diff_forward - G_finite_diff_back)/(dist_delta);
  //fill L
  // return coordinates to original
  lmp->atom->x[iat][dim] = r0;
  sna_global->compute_array();
  return finite_diff_lap;
}

 

 SNAPJastrow::LogValueType SNAPJastrow::evaluateGL(const ParticleSet& P,
                          ParticleSet::ParticleGradient& G,
                          ParticleSet::ParticleLaplacian& L,
                          bool fromscratch){

    RealType dist_delta = 0.001; // TODO: find units
    // compute gradient, i.e., pull gradient out from lammps.
    RealType dist_delta = 0.001; // TODO: find units
    std::cout << "in computeGL" <<std::endl;
    // compute gradient, i.e., pull gradient out from lammps.
    double grad_val;
    for (int ig = 0; ig < P.groups(); ig++) {
      for (int iel = P.first(ig); iel < P.last(ig); iel++){ // loop over elements in each group
        // loop over elecs in each group
        for (int dim = 0; dim < OHMMS_DIM; dim++){
          int row = (iel*3)+dim + 1;
        // loop over elecs in each group
        for (int dim = 0; dim < OHMMS_DIM; dim++){
          int row = (iel*3)+dim + 1;
          for (int nv =0; nv <ncoeff; nv ++){
            int col = nv;
            double grad_val = sna_global->array[row][col];
            double grad_val = sna_global->array[row][col];
            G[iel][dim] += snap_beta[iel][col]*grad_val*hartree_over_ev*bohr_over_ang;
            L[iel] -= FD_Lap(P, iel, dim, row, col, snap_beta, dist_delta); 
          }
        }
      }
        }
      }
    }
    log_value_ = evaluateLog(P,G,L);
    return log_value_;
   }
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
    double ESNAP=0;
    for (int i = 0; i < (Nelec); i++){
      update_lmp_pos(P,lmp,sna_global,i,false);
    }
    calculate_ESNAP(P, sna_global, snap_beta, ESNAP);
    log_value_ = static_cast<SNAPJastrow::LogValueType>(ESNAP);
            
    return log_value_;
  }

    SNAPJastrow::GradType SNAPJastrow::evalGrad(ParticleSet& P, int iat){
     update_lmp_pos(P, proposed_lmp, proposed_sna_global, iat, true);
     update_lmp_pos(P, proposed_lmp, proposed_sna_global, iat, true);
     GradType grad_iat;
     for (int dim=0; dim < OHMMS_DIM; dim++){
      int row =(3*iat)+dim+1;
      for (int k = 0; k < ncoeff ; k++){
        for (int n = 0; n < lmp->atom->ntypes; n++){
          int col = (n*ncoeff)+k;
          grad_iat[dim] += -snap_beta[P.GroupID[iat]][k]*proposed_sna_global->array[row][col]*hartree_over_ev/bohr_over_ang;
        }
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
              dhpsioverpsi[kk] = -RealType(0.5) * RealType(Sum(lapLogPsi[k])) - RealType(Dot(P.G, gradLogPsi[k]));
            }
          }
        }
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
              dhpsioverpsi[kk] = -RealType(0.5) * RealType(Sum(lapLogPsi[k])) - RealType(Dot(P.G, gradLogPsi[k]));
            }
          }
        }
      }
  void SNAPJastrow::evaluateDerivatives(ParticleSet& P, const opt_variables_type& optvars, Vector<ValueType>& dlogpsi, Vector<ValueType>& dhpsioverpsi)
  {
    //std::cout<< "in evaluateDerivatives" <<std::endl;
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
              //std::cout<< "at laplogpsi" <<std::endl;
              //std::cout << lapLogPsi[k] << std::endl;
              dhpsioverpsi[kk] = -RealType(0.5) * RealType(Sum(lapLogPsi[k]));
              for (int i =0; i < Nelec; i++){
                //std::cout<< "at P.G" <<std::endl;
                //std::cout << P.G[i] << std::endl;
                //std::cout<< "at gradlogpsi" << std::endl;
                //std::cout << gradLogPsi[k][i] << std::endl;
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
            dLogPsi[p] = 0.0;
        }
      
            dLogPsi[p] = 0.0;
        }
      
        for (int k = 0; k < myVars.size(); k++){
          int kk = myVars.where(k);
          if (kk < 0)
            continue;
          if (rcsingles[k]){
            evaluate_fd_derivs(P, kk);
            //std::cout << "we are on coeeff" << kk << " out of " << myVars.size() << std::endl;
            dlogpsi[kk] = ValueType(dLogPsi[kk]);
            //std::cout<< "finished assigning dlogpsi"<<std::endl;
          }
        } 
      }
    }

  void SNAPJastrow::evaluate_fd_derivs(ParticleSet& P, int coeff_idx){
        std::vector<std::vector<double>> fd_coeff(snap_beta);
        std::vector<std::vector<double>> bd_coeff(snap_beta);
  void SNAPJastrow::evaluate_fd_derivs(ParticleSet& P, int coeff_idx){
        std::vector<std::vector<double>> fd_coeff(snap_beta);
        std::vector<std::vector<double>> bd_coeff(snap_beta);
        RealType fd_u, bd_u;
        RealType coeff_delta = 1e-6;
        RealType dist_delta = 1e-4;
        int el = int(coeff_idx/ncoeff);
        int coeff = coeff_idx%ncoeff;
        fd_coeff[el][coeff] = snap_beta[el][coeff] + coeff_delta;
        bd_coeff[el][coeff] = snap_beta[el][coeff] - coeff_delta;
        //std::cout << "in eval fd: we are on coeeff" << coeff_idx << " out of " << ncoeff << std::endl;
        //std::cout << "and el is " << el << "and coeff is " << coeff <<std::endl;
        calculate_ESNAP(P, sna_global, fd_coeff, fd_u,false);
        calculate_ESNAP(P, sna_global, bd_coeff, bd_u,false);
        dLogPsi[coeff_idx] = (fd_u - bd_u)/(2*coeff_delta); //units handled elsewhere
        calculate_ddc_gradlap_lammps(P, dist_delta, coeff_delta, fd_coeff, bd_coeff, coeff_idx);
    }

SNAPJastow::evaluatelog(const ParticleSet& P,
                                  ParticleSet::ParticleGradient& G,
                                  ParticleSet::ParticleLaplacian& L,
                                  bool fromscratch){


  void SNAPJastrow::calculate_ddc_gradlap_lammps(ParticleSet& P,RealType dist_delta, RealType coeff_delta, std::vector<std::vector<double>>& fd_coeff, std::vector<std::vector<double>>& bd_coeff, int cur_val){
    SNAPJastrow::GradDerivVec ddc_grad_forward_val(Nions+Nelec);
    SNAPJastrow::GradDerivVec ddc_grad_back_val(Nions+Nelec);
    SNAPJastrow::ValueDerivVec ddc_lap_forward_val(Nions+Nelec);
    SNAPJastrow::ValueDerivVec ddc_lap_back_val(Nions+Nelec);
    for (int nv = 0; nv < ncoeff; nv++){
      int col = nv;
      for (int iel =0; iel < Nelec+Nions; iel ++){ 
        int coeff_start = (iel*ncoeff);
        for (int dim = 0; dim < OHMMS_DIM; dim++){
          int row = (iel*3)+dim + 1;
          double grad_val = sna_global->array[row][col];
          ddc_grad_forward_val[iel][dim] += fd_coeff[iel][nv]*grad_val*hartree_over_ev*bohr_over_ang;
          ddc_grad_back_val[iel][dim] += bd_coeff[iel][nv]*grad_val*hartree_over_ev*bohr_over_ang;

     
          double finite_diff_lap_forward = FD_Lap(P, iel, dim, row, nv,fd_coeff, dist_delta);
          double finite_diff_lap_back = FD_Lap(P, iel, dim, row, nv,bd_coeff, dist_delta);
     
          double finite_diff_lap_forward = FD_Lap(P, iel, dim, row, nv,fd_coeff, dist_delta);
          double finite_diff_lap_back = FD_Lap(P, iel, dim, row, nv,bd_coeff, dist_delta);
          ddc_lap_forward_val[iel] -=  finite_diff_lap_forward*hartree_over_ev*bohr_over_ang;
          ddc_lap_back_val[iel] -=  finite_diff_lap_back*hartree_over_ev*bohr_over_ang;
        } //end dim
      } //end iel
    } //end nv
    lapLogPsi[cur_val] = (ddc_lap_forward_val-ddc_lap_back_val)/(2*coeff_delta);
    lapLogPsi[cur_val] = (ddc_lap_forward_val-ddc_lap_back_val)/(2*coeff_delta);
    gradLogPsi[cur_val] = (ddc_grad_forward_val - ddc_grad_back_val)/(2*coeff_delta);
    //std::cout<< "finished assinging the value to dalpha array" <<std::endl;
  }
    
    
   /*
     calculates esnap based on a set of coefficients manually in qmcpack
  used to see impact of small change in coefficients on snap energy (needed to calculated d E/d beta)
  without having to internally change the lammps object.
  */
  void SNAPJastrow::calculate_ESNAP(const ParticleSet& P, LAMMPS_NS::ComputeSnap* snap_global, std::vector<std::vector<double>> coeff, double& new_u ){
  void SNAPJastrow::calculate_ESNAP(const ParticleSet& P, LAMMPS_NS::ComputeSnap* snap_global, std::vector<std::vector<double>> coeff, double& new_u ){
    double ESNAP_all = 0;
    double ESNAP_elec=0;
    double ESNAP_ion=0;


    for (int ig = 0; ig < P.groups(); ig++) {
      esnap_elec=0;
      for (int k = 0; k < ncoeff; k++){
        bispectrum_val = snap_global->array[0][(ig*ncoeff) + k]; //block of bispectrum + current component to add.
        std::cout << "bispectrum val in elec "<< bispectrum_val << std::endl;
        esnap_elec += coeff[ig][k] * bispectrum_val;
      }
      esnap_all += esnap_elec;
    }
    std::cout << "end of electron summation" << std::endl;
    esnap_ion = 0;
    for (int k = 0; k < ncoeff; k++){
      bispectrum_val = snap_global->array[0][((lmp->atom->ntypes-1)*ncoeff) + k];// will need fixed for more than one ion group
      esnap_ion += coeff[(lmp->atom->ntypes-1)][k] * bispectrum_val;// will need fixed for more than one ion group.
    }
    esnap_all += esnap_ion;
    new_u = esnap_all*hartree_over_ev;

    std::cout << "exiting calculate snap" <<std::endl; 
    return;
  }


  void SNAPJastrow::evaluateRatios(const VirtualParticleSet& VP, std::vector<ValueType>& ratios){
         //std::cout<< "we are in evaluate ratios" << std::endl;
         //std::cout<< "ratios is of size " << ratios.size() << std::endl;
         //std::cout<< "VP.refptcl is  " << VP.refPtcl << std::endl;
        // create lmp object 
        // done in constructer.
        for (int r = 0; r < ratios.size(); r++){
          for (int dim= 0; dim < OHMMS_DIM; dim ++){
            //std::cout << "lamp pos" << vp_lmp->atom->x[VP.refPtcl][dim]<<std::endl; 
            //std::cout << "vp pos " << VP.R[r][dim]<<std::endl;
            // manually update posiition of ref particle to k position.
            vp_lmp->atom->x[VP.refPtcl][dim] = VP.R[r][dim];
          }
          //std::cout<< "made it through pos reassingment" <<std::endl;
          vp_sna_global->compute_array();
          //std::cout<< "made it through updating vp array" <<std::endl;
          // calculate Eold
          double Eold;
          calculate_ESNAP(VP.getRefPS(),sna_global, snap_beta, Eold, true);
         //calculate Enew
          double Enew;
          calculate_ESNAP(VP.getRefPS(),vp_sna_global, snap_beta, Enew, false);
          //store ratio

          ratios[r] = std::exp(static_cast<ValueType>(Enew-Eold));
          //std::cout<< "made it through filling ratios" <<std::endl;
        }
        for (int dim= 0; dim < OHMMS_DIM; dim ++){
          vp_lmp->atom->x[VP.refPtcl][dim] = VP.getRefPS().R[VP.refPtcl][dim];
        }
    return;
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
  /////////////////////////////////// MC Related functions /////////

  void SNAPJastrow::acceptMove(ParticleSet& P, int iat, bool safe_to_delay){
      // create new lammps object with accepted positions
      //std::cout<< "in accept move"<<std::endl;
      update_lmp_pos(P,lmp,sna_global,iat,false);
      log_value_ = evaluateLog(P,P.G,P.L);

      // the G, dlogpsi, gradlogpsi, and laplogpsi need to be updated but only for single particle.
  }

void SNAPJastrow::registerData(ParticleSet& P, WFBufferType& buf){
    log_value_ = evaluateLog(P, P.G, P.L);
}
void SNAPJastrow::registerData(ParticleSet& P, WFBufferType& buf){
    log_value_ evaluateLog(P,P.G,P.L);
}

  void SNAPJastrow::copyFromBuffer(ParticleSet& P, WFBufferType& buf){
  void SNAPJastrow::copyFromBuffer(ParticleSet& P, WFBufferType& buf){
    
}
}

