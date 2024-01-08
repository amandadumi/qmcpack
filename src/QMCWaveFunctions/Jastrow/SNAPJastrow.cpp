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
    Ions(ions)

{
// // reserve just one process for lammps
    int n,me,nprocs;
    int nprocs_lammps;
//     // get rank and total of all procs
    MPI_Comm_rank(MPI_COMM_WORLD,&me);
    MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
//     // requesting just one process be used for lammps
    //std::cout << "nprocs from mpi_comm_world is " << nprocs << std::endl;
   //std::cout << "this rank is " << me << std::endl;
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
  
  //
  //comm_lammps = MPI_COMM_WORLD;

    if (NIonGroups >1){
      app_log() << "WARNING: SNAP jastrow does not currently support more than one ion group" <<std::endl;
    }
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
    lmp = initialize_lammps(els, rcut);
    proposed_lmp = initialize_lammps(els, rcut);
    vp_lmp = initialize_lammps(els,rcut); //may not need, but lets be extra sure wires don't get crossed
    sna_global = static_cast<LAMMPS_NS::ComputeSnap*>(lmp->modify->get_compute_by_id("sna_global"));
    proposed_sna_global = static_cast<LAMMPS_NS::ComputeSnap*>(proposed_lmp->modify->get_compute_by_id("sna_global"));
    vp_sna_global = static_cast<LAMMPS_NS::ComputeSnap*>(vp_lmp->modify->get_compute_by_id("sna_global"));
    //std::cout<< "made the lammps objects"<<std::endl;
    snap_beta = std::vector<std::vector<double>>(lmp->atom->ntypes, std::vector<double>(ncoeff,0.001));
    for (int i=0; i < lmp->atom->ntypes; i++){
      for (int k = 0; k < ncoeff;k++){
        std::stringstream name;
        name << "snap_coeff_" << i;
        name << "_"  << k ;
        //std::cout<< name.str() <<std::endl;
        myVars.insert(name.str(), snap_beta[i][k], true);
      }
    }
    resizeWFOptVectors(); 
    grad_u.resize(Nelec);
    lap_u.resize(Nelec);
    //std::cout<< "we finished initialization"<<std::endl;
}

SNAPJastrow::~SNAPJastrow(){
  // MPI_Comm_free(&comm_lammps);
  delete lmp;
  delete proposed_lmp;
  delete vp_lmp;
//  MPI_Comm_free(&comm_lammps);
  std::cout << "able to deconstruct our lammps functions" << std::endl;

}


void SNAPJastrow::set_coefficients(std::vector<double> id_coeffs,int id){
  //std::cout<< id_coeffs.size() << "is the size of coefficients" <<std::endl;
 // std::cout<< ncoeff << "is the number of coeffs" <<std::endl;
  if (id_coeffs.size() != ncoeff){
    app_warning() << " Warning wrong number of coefficents for snap jastrow" << std::endl;
  }
  int kk=0;
  for (int i=0; i < id_coeffs.size(); i++){
    snap_beta[id][i] = id_coeffs[i];
  }

}

LAMMPS_NS::LAMMPS* SNAPJastrow::initialize_lammps(const ParticleSet& els, double rcut){
    //std::cout << "in initialize_lammps" <<std::endl;
    const char *lmpargv[] {"liblammps","-log","lammps.out","-screen","lammps_screen.out"};
    int lmpargc = sizeof(lmpargv)/sizeof(const char *);
    LAMMPS_NS::LAMMPS *this_lmp;


  /* open LAMMPS input script */

    this_lmp = new LAMMPS_NS::LAMMPS(lmpargc, (char **)lmpargv, comm_lammps);

    // std::cout << "created bare lammps instance" <<std::endl;
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
   // std::cout << "before ion creation" <<std::endl;
    // std::cout << "total electron groups is" << std::to_string(els.groups()) << std::endl;
    for (int ig = 0; ig < els.groups(); ig++) { // loop over groups
      // label groups in lamps
      for (int iat = els.first(ig); iat < els.last(ig); iat++) { // loop over elements in each group
        // place atom in boxes according to their group.
        temp_command = std::string("create_atoms ") + std::to_string(ig+1) + " single " + std::to_string((els.R[iat][0]+.1)/bohr_over_ang) + "  " + std::to_string((els.R[iat][1]+.01*iat)/bohr_over_ang)  + " " + std::to_string((els.R[iat][2]+.1*iat+.01)/bohr_over_ang)+ " units box";  
        //std::cout << temp_command << std::endl;
        this_lmp->input->one(temp_command);
      }
    }
      const SpeciesSet& tspecies(Ions.getSpeciesSet());
      for (int ig = 0; ig < Ions.groups(); ig++) { // loop over groups
          //std::cout << "ion species is" << tspecies.speciesName[ Ions.GroupID[ig] ] << std::endl;
          temp_command = std::string("group ions_"+ std::to_string(ig)  + " type " + std::to_string(els.groups()+ig+1));
          //std::cout << temp_command << std::endl;
          temp_command = std::string("mass "+ std::to_string(els.groups()+ig+1)) + " 1.00";
          //std::cout << temp_command << std::endl;
          this_lmp->input->one(temp_command);
          for (int iat = Ions.first(ig); iat < Ions.last(ig); iat++) { // loop over elements in each group
            temp_command = std::string("create_atoms "  + std::to_string(els.groups()+ig+1) + " single ") + std::to_string(Ions.R[iat][0]/bohr_over_ang) + "  " + std::to_string(Ions.R[iat][1]/bohr_over_ang)  + " " + std::to_string(Ions.R[iat][2]/bohr_over_ang) + " units box";  
            //std::cout << temp_command << std::endl;
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
      //setting rcut to be the same for each type, though may be interesting to try to automate this
      temp_command = std::string("variable rad_type_1 equal ") + std::to_string(rcut/bohr_over_ang);
      this_lmp->input->one(temp_command);
      temp_command = std::string("variable rad_type_2 equal ") + std::to_string(rcut/bohr_over_ang);
      this_lmp->input->one(temp_command);
      temp_command = std::string("variable rad_type_3 equal ") + std::to_string(rcut/bohr_over_ang);
      this_lmp->input->one(temp_command);
      this_lmp->input->one("variable	wj1 equal 1.0");
      this_lmp->input->one("variable	wj2 equal 1.0");
      this_lmp->input->one("variable	wj3 equal 0.96");
      this_lmp->input->one("variable	quadratic equal 0");
      this_lmp->input->one("variable	bzero equal 1");
      this_lmp->input->one("variable	switch equal 0");
      this_lmp->input->one("variable snap_options string \"${rcutfac} ${rfac0} ${twojmax} ${rad_type_1} ${rad_type_2} ${rad_type_3} ${wj1} ${wj2} ${wj3} rmin0 ${rmin0} quadraticflag ${quadratic} bzeroflag ${bzero} switchflag ${switch}\"");

    //snap needs some reference pair potential, but doesn't effect parts we are using. 

      this_lmp->input->one(" pair_style zero 10.0");
      this_lmp->input->one("pair_coeff * *");
      //TODO: generalize with loop over atom types
      this_lmp->input->one("compute sna_global all snap ${snap_options}"); 
      this_lmp->input->one("thermo 100");
      this_lmp->input->one("thermo_style   custom  c_sna_global[1][11] c_sna_global[2][1]");
      this_lmp->input->one("run            0 pre no post no");
      //std::cout<< "lammps ntypes is " << this_lmp->atom->ntypes << std::endl;

    return this_lmp;
  }

  void SNAPJastrow::update_lmp_pos(const ParticleSet& P, LAMMPS_NS::LAMMPS* lmp_pntr, LAMMPS_NS::ComputeSnap* snap_array,int iat, bool proposed){
      if (proposed){
        //std::cout<< "in update positions: proposed" <<std::endl;
        for (int dim = 0; dim < OHMMS_DIM; dim++){
          lmp_pntr->atom->x[iat][dim] = P.activeR(iat)[dim]/bohr_over_ang;
        }
      }
      else{
        //std::cout<< "in update positions: set val" <<std::endl;
        for (int dim = 0; dim < OHMMS_DIM; dim++){
          lmp_pntr->atom->x[iat][dim] = P.R[iat][dim]/bohr_over_ang;
        }
      }
      snap_array->compute_array();
  }

double SNAPJastrow::FD_Lap(const ParticleSet& P,int iat, int dim, int coeff, int ntype, std::vector<std::vector<double>> coeffs, double dist_delta, bool bispectrum_only){
  //std::cout << "in FD_Lap"<<std::endl;
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
  G_finite_diff_forward = this_coeff * sna_global->array[row][(ntype*ncoeff)+coeff] * hartree_over_ev/bohr_over_ang;
  
  //backward direction
  RealType rm  = r0 - (dist_delta/2);
  lmp->atom->x[iat][dim] = rm;
  sna_global->compute_array();
  G_finite_diff_back = this_coeff * sna_global->array[row][(ntype*ncoeff)+coeff] * hartree_over_ev/bohr_over_ang;
  //fill L
  double finite_diff_lap = (G_finite_diff_forward - G_finite_diff_back)/(dist_delta);
  
  // return coordinates to original
  lmp->atom->x[iat][dim] = r0;
  sna_global->compute_array();
  return finite_diff_lap;
}


 SNAPJastrow::LogValueType SNAPJastrow::evaluateGL(const ParticleSet& P,
                          ParticleSet::ParticleGradient& G,
                          ParticleSet::ParticleLaplacian& L,
                          bool fromscratch){
                          // std::cout << "in evaluateGL" <<std::endl;
                          return evaluateLog(P,G,L);
                          }
 
 
 void SNAPJastrow::computeGL(const ParticleSet& P){
    RealType dist_delta = 0.001; // TODO: find units
    //std::cout << "in computeGL" <<std::endl;
    // compute gradient, i.e., pull gradient out from lammps.
    double grad_val;
    for (int ig = 0; ig < P.groups(); ig++) {
      for (int iel = P.first(ig); iel < P.last(ig); iel++){ // loop over elements in each group
        grad_u[iel] = 0;
        lap_u[iel] = 0;
        for (int dim = 0; dim < OHMMS_DIM; dim++){
          int row = (iel*3)+dim + 1;
            // std::cout << "considering electron "  << iel << std::endl;
          for (int n = 0; n < lmp->atom->ntypes; n++){
            for (int k =0; k < ncoeff; k ++){
              int col = (n*ncoeff)+k;
              // std::cout<< "row " << row << " col " << col <<std::endl;
              grad_val = sna_global->array[row][col];
              grad_u[iel][dim] += snap_beta[n][k]*grad_val*hartree_over_ev/bohr_over_ang;
              lap_u[iel] += FD_Lap(P, iel, dim, k, n, snap_beta, dist_delta, false)/bohr_over_ang; 
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
    //std::cout << "we  in evaluate log" << std::endl;
    double esnap;
    //std::cout << "we  in evaluate log" << std::endl;
    for (int i = 0; i < Nelec; i++){
      // std::cout << "iteration " << i << std::endl;
      update_lmp_pos(P,lmp,sna_global,i,false);
    }
    //std::cout << "evaluatelog updated positions" << std::endl;
    calculate_ESNAP(P, sna_global, snap_beta, esnap, true);
    //std::cout << "evaluatelog calculated esnap" << std::endl;
    computeGL(P);
    //std::cout << "evaluatelog calculated gradients" << std::endl;
    for (int iel = 0; iel < Nelec; iel++){
      G[iel] += grad_u[iel];
      L[iel] += lap_u[iel];
    }
    //std::cout << "evaluatelog reassigned gradient info" << std::endl;
    
    log_value_ = static_cast<SNAPJastrow::LogValueType>(esnap);
            
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
          grad_iat[dim] += snap_beta[n][k]*proposed_sna_global->array[row][col]*hartree_over_ev/bohr_over_ang;
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
            //std::cout<< "start assigning dlogpsi"<<std::endl;
            dlogpsi[kk] = ValueType(dLogPsi[kk]);
            //std::cout<< "finished assigning dlogpsi"<<std::endl;
          }
        } 
      }
    }

  void SNAPJastrow::evaluate_linear_derivs(ParticleSet& P, int coeff_idx){
    //TODO: this has a bug in that the el will try to be used to explore the atom info to in which case the getGroupID throws a fitting assertion error?
    //std::cout<< "coeff id is " << coeff_idx <<std::endl;
    int el = int(coeff_idx/ncoeff); // this coeff is apart of snap with el as central atom. This will also work with multple species of same type as coeeffs/ncoeff will be the type 
    int row_id;
    if (el < P.getTotalNum()){ // assuming elecs are stored first.
      row_id = P.getGroupID(el);
    }
    else{
      row_id = Ions.getGroupID(el);
    }
    int coeff = coeff_idx%ncoeff; // which coeff of this el are we on.
    computeGL(P);
    dLogPsi[coeff_idx] = sna_global->array[0][coeff_idx];// dlogpsi will be bispectrum component  
    for (int dim = 0; dim < OHMMS_DIM; dim++){ // loop over dim to get grad vec.
      gradLogPsi[coeff_idx][el][dim] += sna_global->array[(el*OHMMS_DIM)+dim+1][(row_id*ncoeff)+coeff];
      lapLogPsi[coeff_idx][el] -=0; 
      //lapLogPsi[coeff_idx] -= FD_Lap(P, el, dim, coeff, el, snap_beta, 1e-6, true)/bohr_over_ang;
    
    }
  }

  void SNAPJastrow::evaluate_fd_derivs(ParticleSet& P, int coeff_idx){
        // std::cout << "in evaluate_fd_derivs"<<std::endl;
        std::vector<std::vector<double>> fd_coeff(snap_beta);
        std::vector<std::vector<double>> bd_coeff(snap_beta);
        RealType fd_u, bd_u;
        RealType coeff_delta = 1e-6;
        RealType dist_delta = 1e-4;
        int el = int(coeff_idx/ncoeff);
        int coeff = coeff_idx%ncoeff;
        fd_coeff[el][coeff] = snap_beta[el][coeff] + coeff_delta;
        bd_coeff[el][coeff] = snap_beta[el][coeff] - coeff_delta;
        // std::cout << "in eval fd: we are on coeeff" << coeff_idx << " out of " << ncoeff << std::endl;
        // std::cout << "and el is " << el << "and coeff is " << coeff <<std::endl;
        calculate_ESNAP(P, sna_global, fd_coeff, fd_u, false);
        calculate_ESNAP(P, sna_global, bd_coeff, bd_u, false);
        dLogPsi[coeff_idx] = (fd_u - bd_u)/(2*coeff_delta); //units handled elsewhere
        calculate_ddc_gradlap_lammps(P, dist_delta, coeff_delta, fd_coeff, bd_coeff, coeff_idx);
    }

SNAPJastow::evaluatelog(const ParticleSet& P,
                                  ParticleSet::ParticleGradient& G,
                                  ParticleSet::ParticleLaplacian& L,
                                  bool fromscratch){


  void SNAPJastrow::calculate_ddc_gradlap_lammps(ParticleSet& P,RealType dist_delta, RealType coeff_delta, std::vector<std::vector<double>>& fd_coeff, std::vector<std::vector<double>>& bd_coeff, int cur_val){
    //std::cout << "in calculate_ddc_gradlap_lammps"<<std::endl;
    SNAPJastrow::GradDerivVec ddc_grad_forward_val(Nelec);
    SNAPJastrow::GradDerivVec ddc_grad_back_val(Nelec);
    SNAPJastrow::ValueDerivVec ddc_lap_forward_val(Nelec);
    SNAPJastrow::ValueDerivVec ddc_lap_back_val(Nelec);
    for (int iel = 0; iel < Nelec; iel++){ // for each particle, step down through rows
      for (int dim = 0; dim < OHMMS_DIM; dim++){ // step down through subblock of each 
        int row = (iel*3) + dim + 1; // organized by xyz for each particle.
        for (int n = 0; n < lmp->atom->ntypes; n++){ // there will be terms from all other particles.
          for (int k = 0; k < ncoeff; k++){ // over all of the coeffs
            //std::cout<< "iel " << iel << " k " << k << " n " << n << " dim "<< dim << std::endl; 
            double grad_val = sna_global->array[row][(n*ncoeff)+k];
            ddc_grad_forward_val[iel][dim] += fd_coeff[n][k]*grad_val*hartree_over_ev/bohr_over_ang;
            ddc_grad_back_val[iel][dim] += bd_coeff[n][k]*grad_val*hartree_over_ev/bohr_over_ang;

            ddc_lap_forward_val[iel] += FD_Lap(P, iel, dim, k, n, fd_coeff, dist_delta, false)/bohr_over_ang;
            ddc_lap_back_val[iel] += FD_Lap(P, iel, dim, k, n, bd_coeff, dist_delta, false)/bohr_over_ang;
          } //end dim
        } //end ncoeff 
      } //end ntype
    } //end nelec
    //std::cout<< "made it to assinging the value to dalpha array" <<std::endl;
    lapLogPsi[cur_val] = (ddc_lap_forward_val - ddc_lap_back_val)/(2*coeff_delta);
    gradLogPsi[cur_val] = (ddc_grad_forward_val - ddc_grad_back_val)/(2*coeff_delta);
    //std::cout<< "finished assinging the value to dalpha array" <<std::endl;
  }
    
    
   /*
     calculates esnap based on a set of coefficients manually in qmcpack
  used to see impact of small change in coefficients on snap energy (needed to calculated d E/d beta)
  without having to internally change the lammps object.
  */
  void SNAPJastrow::calculate_ESNAP(const ParticleSet& P, LAMMPS_NS::ComputeSnap* snap_global, std::vector<std::vector<double>> coeff, double& new_u,bool store_u=false){
    //std::cout << "in calculate snap" <<std::endl; 
    double esnap_all = 0;
    double esnap_elec;
    double esnap_ion;
    double bispectrum_val;
    // calculate electron contribution
    // the global array is summed over groups of atoms of the same type. thus here we just need to sum over groups.
    for (int ig = 0; ig < P.groups(); ig++) {
      esnap_elec=0;
      for (int k = 0; k < ncoeff; k++){
        bispectrum_val = snap_global->array[0][(ig*ncoeff) + k]; //block of bispectrum + current component to add.
        //std::cout << "bispectrum val in elec "<< bispectrum_val << std::endl;
        esnap_elec += coeff[ig][k] * bispectrum_val;
      }
      esnap_all += esnap_elec;
    }
    //std::cout << "end of electron summation" << std::endl;
    esnap_ion = 0;
    for (int ig = 0; ig < Ions.groups(); ig++) {
      for (int k = 0; k < ncoeff; k++){
        bispectrum_val = snap_global->array[0][((P.groups()+ig)*ncoeff) + k];// will need fixed for more than one ion group
        esnap_ion += coeff[P.groups()+ig][k] * bispectrum_val;// will need fixed for more than one ion group.
      }
    }
    esnap_all += esnap_ion;
    //std::cout << "end of ion summation" <<std::endl; 
    new_u = -esnap_all*hartree_over_ev;
    return;
  }


  void SNAPJastrow::evaluateRatios(const VirtualParticleSet& VP, std::vector<ValueType>& ratios){
         //std::cout<< "we are in evaluate ratios" << std::endl;
         //std::cout<< "ratios is of size " << ratios.size() << std::endl;
         //std::cout<< "VP.refptcl is  " << VP.refPtcl << std::endl;
        // create lmp object 
        // done in constructer.
        double Eold;
        calculate_ESNAP(VP.getRefPS(),sna_global, snap_beta, Eold, true);
        for (int r = 0; r < ratios.size(); r++){
          for (int dim= 0; dim < OHMMS_DIM; dim ++){
            //std::cout << "lamp pos" << vp_lmp->atom->x[VP.refPtcl][dim]<<std::endl; 
            //std::cout << "vp pos " << VP.R[r][dim]<<std::endl;
            // manually update posiition of ref particle to k position.
            vp_lmp->atom->x[VP.refPtcl][dim] = VP.R[r][dim];
          }
          //std::cout<< "made it through pos reassingment" <<std::endl;
          vp_sna_global->compute_array();
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

  SNAPJastrow::PsiValueType SNAPJastrow::ratioGrad(ParticleSet& P, int iat, GradType& grad_iat){
     SNAPJastrow::PsiValueType ratio = SNAPJastrow::ratio(P,iat);
     for (int dim = 0; dim < OHMMS_DIM; dim++){
      for (int k = 0; k < ncoeff ; k++){
        for (int n = 0; n < lmp->atom->ntypes; n++)
          grad_iat[dim] += snap_beta[n][k]*proposed_sna_global->array[(3*iat)+dim+1][(n*ncoeff)+k]*hartree_over_ev/bohr_over_ang;
      }
     }
     return ratio;
  }

  SNAPJastrow::PsiValueType SNAPJastrow::ratio(ParticleSet& P, int iat){
    //std::cout<< "in calc ratio  "  <<std::endl;
    //std::cout<< "iat is   "  <<  iat <<std::endl;
    for (int i = 0; i < (Nelec);i++){
      update_lmp_pos(P, proposed_lmp, proposed_sna_global, i, true);
    }
     double Eold;
     calculate_ESNAP(P, sna_global, snap_beta, Eold, true);
     double Enew;
     calculate_ESNAP(P, proposed_sna_global, snap_beta, Enew, true);
     //  calculate_ESNAP(P, sna_global, snap_beta, Eold, false);
     SNAPJastrow::PsiValueType ratio = std::exp(static_cast<SNAPJastrow::PsiValueType>(Enew-Eold));
      // std::cout<< "ratio in function is " << ratio <<std::endl;
      // std::cout<< "inverse  is " << std::exp(static_cast<SNAPJastrow::PsiValueType>(Eold-Enew))<<std::endl;
      // std::cout<< "particle specific  is " << std::exp(static_cast<SNAPJastrow::PsiValueType>(Eiold-Einew))<<std::endl;
      // std::cout<< "inverse particle specific is " << std::exp(static_cast<SNAPJastrow::PsiValueType>(Einew-Einew))<<std::endl;
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
  for (int n=0; n< lmp->atom->ntypes; n++){
    for (int k = 0; k < ncoeff; k++){
      snap_beta[n][k] = std::real(myVars[(n*ncoeff)+k]);
    } 
  }
}

std::unique_ptr<WaveFunctionComponent> SNAPJastrow::makeClone(ParticleSet& tpq) const
{
  auto snap_copy = std::make_unique<SNAPJastrow>(std::string("snap"),Ions, tpq, std::string("linear"), twojmax, rcut);
  snap_copy->myVars = myVars;
  snap_copy->lmp = lmp;
  snap_copy->lmp = proposed_lmp;
  snap_copy->sna_global = sna_global;
  snap_copy->proposed_sna_global = proposed_sna_global;
  snap_copy->snap_beta = snap_beta;

  return snap_copy;
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

