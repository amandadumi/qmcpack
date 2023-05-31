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
    std::string iSpecies, eSpecies1, eSpecies2;
    ncoeff = 56;
    app_log() << els.R << std::endl;
    app_log() << ions.R << std::endl;

    initialize_lammps(Ions,els);
    snap_beta = std::vector<std::vector<double>>((Nions+Nelec), std::vector<double>(ncoeff));
    std::vector<double> guess_coeffs  { -0.000000000000 , -0.001487061994 , 0.075808306870, 0.538735683870 , -0.074148039366 , 0.602629813770 , -0.147022424344,  0.117756828488 ,-0.026490439049 ,-0.035162708767 ,0.064315385091 ,-0.131936948089 ,-0.021272860272 ,-0.091171134054 ,-0.024396224398 ,-0.059813132803 ,0.069585393203 ,-0.085344044181 ,-0.155425254597 ,-0.117031758367 ,-0.040956258020 ,-0.084465000389 ,-0.020367513630 ,-0.010730484318 ,-0.054777575658 , 0.050742893747 ,-0.004686334611 ,-0.116372907121 , 0.005542497708 ,-0.126526795635 ,-0.080163926221 ,-0.082426250179 ,-0.010558777281 ,-0.001939058038 ,-0.027907949962 , 0.049483908476 , 0.005103754385 ,-0.054751505141 ,-0.055556071011 ,-0.006026917619 ,-0.060889030109 ,-0.029977673973 ,-0.014987527280 ,-0.006697686658 , 0.017369624409 , 0.047864358817 ,-0.001989812679 , 0.000153530925 ,-0.003862356345 ,-0.009754314198, 0.000777958970,-0.003031424287, 0.015612715209, 0.003210129646, -0.013088799947, 0.001465970755 };
    for (int i=0; i<(Nions+Nelec); i++){
      for (int nc = 0; nc < ncoeff;nc++){
        std::stringstream name;
        name << "snap_coeff_" << i;
        name << "_"  << nc ;
        snap_beta[i][nc] = guess_coeffs[nc];
        if (nc != 0){
          myVars.insert(name.str(), snap_beta[i][nc], true);
        }
      }
    }
    std::cout<< "finished initializing myVars" <<std::endl;
    resizeWFOptVectors(); 
}

SNAPJastrow::~SNAPJastrow(){
  delete lmp;
}

void SNAPJastrow::initialize_lammps(const ParticleSet& ions, ParticleSet& els){
    const char *lmpargv[] {"liblammps","-log","none"};
    int lmpargc = sizeof(lmpargv)/sizeof(const char *);
    lmp = new LAMMPS_NS::LAMMPS(lmpargc, (char **)lmpargv,MPI_COMM_WORLD);
    lmp->input->one("units  metal");
    lmp->input->one("atom_style  atomic");
    // TODO: will this be set by a qmc box? probably.
    lmp->input->one("boundary  f f f ");
    lmp->input->one("neighbor 1.9 bin");
    lmp->input->one("neigh_modify every 1 delay 1 check yes ");
    //TODO: what will box region be pased on QMC object. Cell?
    lmp->input->one("region	mybox block -50 50 -50 50 -50 50");
    // create a box that will contain the number of species equal to the number of groups.
    std::string temp_command = std::string("create_box ") + std::to_string(3) +  " mybox";
    lmp->input->one(temp_command);
    // add atoms to the groups
    char ptype_snap_id; // the initial coeffs (with label from *.snapcoeff file) to use for each group of particle group
    double bohr_over_ang = 1.88973;// to convert qmcpack storage of bohr into angstrom for lammps
    lmp->input->one("group e_u type 1");
    lmp->input->one("group e_d type 2");
    lmp->input->one("group i type 3");
    lmp->input->one("group elecs type 1 2");
    lmp->input->one("group all type 1 2 3");
     lmp->input->one("mass 1 .95");
  lmp->input->one("mass 2 .95");
  lmp->input->one("mass 3 1");
    for (int ig = 0; ig < els.groups(); ig++) { // loop over groups
        // label groups in lamps
        for (int iat = els.first(ig); iat < els.last(ig); iat++) { // loop over elements in each group
           // place atom in boxes according to their group.
           temp_command = std::string("create_atoms ") + std::to_string(ig+1) + " single " + std::to_string(els.R[iat][0]/bohr_over_ang) + "  " + std::to_string(els.R[iat][1]/bohr_over_ang)  + " " + std::to_string(els.R[iat][2]/bohr_over_ang)+ " units box";  
           std::cout << temp_command <<std::endl;
           std::cout << els.R[iat][0] << " " << els.R[iat][1] << " " <<els.R[iat][2] <<std::endl;
           std::cout << std::to_string(els.R[iat][0]) << " " << std::to_string(els.R[iat][1]) << " " << std::to_string(els.R[iat][2]) <<std::endl;
           lmp->input->one(temp_command);
         }
     }
    
    for (int ig = 0; ig < ions.groups(); ig++) { // loop over groups
      for (int iat = ions.first(ig); iat < ions.last(ig); iat++) { // loop over elements in each group
        temp_command = std::string("create_atoms 3 single ") + std::to_string(ions.R[iat][0]/bohr_over_ang) + "  " + std::to_string(ions.R[iat][1]/bohr_over_ang)  + " " + std::to_string(ions.R[iat][2]/bohr_over_ang) + " units box";  
        std::cout << temp_command <<std::endl;
        lmp->input->one(temp_command);
      }
    }
lmp->input->one("info all out append lammps_log.out");

lmp->input->one("variable 	twojmax equal 2");
lmp->input->one("variable 	rcutfac equal 1.0");
lmp->input->one("variable 	rfac0 equal 0.99363");
lmp->input->one("variable 	rmin0 equal 0");
lmp->input->one("variable 	radelem1 equal 2.3");
lmp->input->one("variable 	radelem2 equal 2.3");
lmp->input->one("variable 	radelem3 equal 2.0");
lmp->input->one("variable	wj1 equal 1.0");
lmp->input->one("variable	wj2 equal 1.0");
lmp->input->one("variable	wj3 equal 0.96");
lmp->input->one("variable	quadratic equal 0");
lmp->input->one("variable	bzero equal 0");
lmp->input->one("variable	switch equal 0");
lmp->input->one("variable snap_options string \"${rcutfac} ${rfac0} ${twojmax} ${radelem1} ${radelem2} ${radelem3} ${wj1} ${wj2} ${wj3} rmin0 ${rmin0} quadraticflag ${quadratic} bzeroflag ${bzero} switchflag ${switch}\"");

lmp->input->one("pair_style zero ${rcutfac}");
    //TODO: generalize with loop over atom types
lmp->input->one("pair_coeff * *"); 

lmp->input->one("variable 	zblcutinner equal 4");
lmp->input->one("variable 	zblcutouter equal 4.8");
lmp->input->one("variable 	zblz equal 73");
lmp->input->one("pair_style zbl ${zblcutinner} ${zblcutouter}");
lmp->input->one("pair_coeff 	* * ${zblz} ${zblz}");
lmp->input->one("compute sna_global all snap ${snap_options}"); 
lmp->input->one("run            0");
double val =  static_cast<double>(lmp->force->pair->eng_vdwl);
        std::cout<<"Value is " << val << std::endl;
std::cout << "we have finished intializing lammps" <<std::endl;
}


SNAPJastrow::LogValueType SNAPJastrow::evaluateGL(const ParticleSet& P,
                        ParticleSet::ParticleGradient& G,
                        ParticleSet::ParticleLaplacian& L,
                        bool fromscratch){
    std::cout << "call evaluateGL" << std::endl;
    void *pos = lmp->atom->x;
    double **x = static_cast<double **> (pos);
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
                int row = (dim * iel)+1;
                G[iel][dim] += static_cast<double>(sna_global->array[row][ncols]);
                //forward direction
                RealType r0 = P.R[iel][dim];
                RealType rp   = r0 + (delta/2);
                //update lammps position
                (*x)[(iel*3)+dim] = rp;
                // recalculate bispectrum stuff
                sna_global -> compute_peratom();
                G_finite_diff_forward = static_cast<double>(sna_global->array[(dim*iel)+1][ncols]);
                //gradient
                //backward direction
                RealType rm   = r0 - (delta/2);
                //update lammps position
                (*x)[(3*iel)+dim] = rm;
                sna_global -> compute_peratom();
                G_finite_diff_back = static_cast<double>(sna_global->array[(dim*iel)+1][ncols]);
                // recalculate bispectrum stuff
                double finite_diff_lap = (G_finite_diff_forward - G_finite_diff_back)/delta;
                //fill L
                L[iel] -=  finite_diff_lap;
                // return coordinates to original
                (*x)[iel+dim] = r0;

              }
        }
    }
   log_value_ = evaluateLog(P,G,L);
   return log_value_;
}




SNAPJastrow::LogValueType SNAPJastrow::evaluateLog(const ParticleSet& P,
                                  ParticleSet::ParticleGradient& G,
                                  ParticleSet::ParticleLaplacian& L){
    app_log() << P.R << std::endl;
    app_log() << Ions.R << std::endl;
    double ESNAP=0;
    sna_global = lmp->modify->get_compute_by_id("sna_global");
    for (int ig = 0; ig < P.groups(); ig++) {
      for (int iel = P.first(ig); iel < P.last(ig); iel++){ // loop over elements in each group
        //ESNAP += snap_beta[iel][0];
        for (int nc=1; nc < ncoeff; nc ++){
          ///lmp->input->one("sna_global->compute_array()");
          sna_global->compute_array();
          std::cout << "bispectrum val is " << sna_global->array[0][(ig*ncoeff) + nc] <<std::endl;
          ESNAP += snap_beta[iel][nc] * sna_global->array[0][(ig*ncoeff) + nc];
        }
      }
    }  
    log_value_ = ESNAP*hartree_over_ev;  
    std::cout<< "ESNAP is " << ESNAP << std::endl;
             
    return log_value_;
}


void SNAPJastrow::evaluateDerivatives(ParticleSet& P, const opt_variables_type& optvars, Vector<ValueType>& dlogpsi, Vector<ValueType>& dhpsioverpsi)
{
    std::cout << "we are in evaluateDerivatives" <<std::endl;
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



void SNAPJastrow::evaluateDerivativesWF(ParticleSet& P, const opt_variables_type& optvars, Vector<ValueType>& dlogpsi)
{
std::cout << "we are in evaluateDerivativesWF" <<std::endl;
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
    std::vector<std::vector<double>> fd_coeff;
    std::vector<std::vector<double>> bd_coeff;
    fd_coeff = std::vector<std::vector<double>>((Nelec + Nions),std::vector<double>(ncoeff));
    bd_coeff = std::vector<std::vector<double>>((Nelec + Nions),std::vector<double>(ncoeff));
    RealType fd_u;
    RealType bd_u;
    RealType delta = 0.1;
    std::cout << "sble to initialize functions." << std::endl;
    for (int i=0; i< (Nelec+Nions); i++){
      for (int nv = 1; nv < ncoeff; nv++){
        std::cout<< "particle i is " << i <<" and coeff is " << nv <<std::endl;
        fd_coeff[i][nv] = snap_beta[i][nv] + delta;
        bd_coeff[i][nv] = snap_beta[i][nv] - delta;
  
        calculate_internal_ESNAP_CD(P, fd_coeff, fd_u);
        calculate_internal_ESNAP_CD(P, bd_coeff, bd_u);
        dlogpsi[nv] = (fd_u - bd_u)/(2*delta); //units handled elsewhere

        std::cout<< "fd_u is  " << fd_u  <<std::endl;
        std::cout<< "bd_u is  " << bd_u  <<std::endl;
        std::cout<< "dlogpsi is  " << dlogpsi[nv]  <<std::endl;
        dLogPsi[nv] -= dlogpsi[nv];
      }
    }
    calculate_ddc_gradlap_lammps(P, delta, fd_coeff, bd_coeff);
    }
 }


void SNAPJastrow::calculate_ddc_gradlap_lammps(ParticleSet& P,RealType delta, std::vector<std::vector<double>> fd_coeff, std::vector<std::vector<double>> bd_coeff){
  std::cout << " in calculate_ddc_gradlap" <<std::endl;
  /** brute force calculation of grad lammp where lammps object will have coefficients updated and reran.*/ 
  // 1. update coeffs to forward direction
    // double* x = new double [OHMMS_DIM*P.getTotalNum()];
    void *pos = lmp->atom->x;
    double **x = static_cast<double **> (pos);
    double G_finite_diff_forward;
    double G_finite_diff_back;
    std::vector<SNAPJastrow::GradDerivVec> ddc_grad_forward;
    std::vector<SNAPJastrow::GradDerivVec> ddc_grad_back;
    std::vector<SNAPJastrow::GradDerivVec> ddc_grad_cd;
    std::cout<< "size of my vars is " << myVars.size() << std::endl;
    ddc_grad_forward.resize(myVars.size(), GradDerivVec(3));
    ddc_grad_back.resize(myVars.size(), GradDerivVec(3));
    ddc_grad_cd.resize(myVars.size(), GradDerivVec(3));


    std::vector<SNAPJastrow::ValueDerivVec> ddc_lap_forward;
    std::vector<SNAPJastrow::ValueDerivVec> ddc_lap_back;
    std::vector<SNAPJastrow::ValueDerivVec> ddc_lap_cd;
    ddc_lap_forward.resize(myVars.size(),SNAPJastrow::ValueDerivVec(Nelec));
    ddc_lap_back.resize(myVars.size(),SNAPJastrow::ValueDerivVec(Nelec));
    ddc_lap_cd.resize(myVars.size(),SNAPJastrow::ValueDerivVec(Nelec));
    std::cout<< "initialized arrays" << std::endl;

      // store original coeefs
      
      // 3. calculate grad and lap
    for (int nv=1; nv < ncoeff; nv++){
      for (int ig = 0; ig < P.groups(); ig++){
        for (int iel = P.first(ig); iel < P.last(ig); iel++){ 
          for (int dim = 0; dim < OHMMS_DIM; dim++){
            int row = (iel*3) + 1;
            int col = nv;
            //G[iel][dim] = lammps_object
            double grad_val = sna_global->array[row][col];
            ddc_grad_forward[iel][dim] += fd_coeff[ig][nv]*grad_val*hartree_over_ev;
            ddc_grad_back[iel][dim] += bd_coeff[ig][nv]*grad_val*hartree_over_ev;
            // calc lap d/dr grad
            //forward direction
            RealType r0 = P.R[iel][dim];
            RealType rp   = r0 + (delta/2);
            //update lammps position
            (*x)[(iel*3)+dim] = rp;
            // recalculate bispectrum stuff
            sna_global->compute_array();
            G_finite_diff_forward = sna_global->array[iel+1][dim];
            //back direction
            RealType rm   = r0 - (delta/2);
            //update lammps position
            (*x)[(iel*3)+dim] = rp;
            // recalculate bispectrum stuff
            sna_global->compute_array();
            G_finite_diff_back = sna_global->array[iel+1][dim];
            double finite_diff_lap = (G_finite_diff_forward - G_finite_diff_back)/delta;
            //fill L
            ddc_lap_forward[iel] -=  fd_coeff[ig][nv]*finite_diff_lap*hartree_over_ev;
            ddc_lap_back[iel] -=  bd_coeff[ig][nv]*finite_diff_lap*hartree_over_ev;
            // return coordinates to original
            (*x)[(iel*3)+dim] = r0;
          }
        }
      }
    }
      // 7. use central difference to get estimate of ddc grad and lap  which populates gradLogPsi and lapLogPsi
      // centra diff grad
      //AD: can probably combine this and above loop. 
      for (int ig = 0; ig < P.groups(); ig++) {
         for (int iel = P.first(ig); iel < P.last(ig); iel++){ 
           for (int dim = 0; dim < OHMMS_DIM; dim++){
             ddc_grad_cd[iel][dim] = (ddc_grad_forward[iel][dim] - ddc_grad_back[iel][dim])/(2.0*delta);
           }
           ddc_lap_cd[iel] = (ddc_lap_forward[iel]-ddc_lap_back[iel])/(2.0*delta);
         }
       }
     gradLogPsi = ddc_grad_cd;
     lapLogPsi  = ddc_lap_cd;
    }

/* calculates esnap based on a set of coefficients manually in qmcpack
used to see impact of small change in coefficients on snap energy (needed to calculated d E/d beta)
without having to internally change the lammps object.
*/
void SNAPJastrow::calculate_internal_ESNAP_CD(ParticleSet& P, std::vector<std::vector<double>> new_coeff, double& new_u){
  double ESNAP_all = 0;
  double ESNAP_elec;
  double ESNAP_ion;
  const auto& d_table = P.getDistTableAB(myTableID);

  const size_t ns = d_table.sources();
  const size_t nt = P.getTotalNum();

  for (int ig = 0; ig < P.groups(); ig++) {
    for (int iel = P.first(ig); iel < P.last(ig); iel++){ // loop over elements in each group
      ESNAP_elec = 0; //  maybe beta_0
      for (int k=1; k < ncoeff; k++){
        // determine atom type
        std::cout << "this loop is at ig " << ig << " iel " << iel << " k " << k <<std::endl;
        double bispectrum_val = sna_global->array[0][(ig*ncoeff) + k]; //block of bispectrum + current component to add.
        std::cout << "bispectrum_val is " << bispectrum_val <<std::endl;
        ESNAP_elec += new_coeff[ig][k] * bispectrum_val  ;
      }
    ESNAP_all += ESNAP_elec;
    }
  }
  std::cout<< "at ions" << std::endl;
  int bispectrum_block = P.groups(); // number of groups for electrons will be block for ions. 2*ncoeff = start of ions.
  //TODO: only written for single ion at the moment.
  for (int s = 0; s <ns; s++){
    ESNAP_ion += snap_beta[bispectrum_block][0];
      for (int k=1; k < ncoeff; k++){
        RealType bispectrum_val = sna_global->array[0][(bispectrum_block*ncoeff) + k];
        ESNAP_ion += snap_beta[bispectrum_block][k] * bispectrum_val  ;
      }
    ESNAP_all += ESNAP_ion;
  }
  new_u = ESNAP_all*hartree_over_ev;
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
void SNAPJastrow::extractOptimizableObjectRefs(UniqueOptObjRefs& opt_obj_refs){opt_obj_refs.push_back(*this);}

void SNAPJastrow::checkInVariablesExclusive(opt_variables_type& active){
  active.insertFrom(myVars);
}

void SNAPJastrow::resetParametersExclusive(const opt_variables_type& active){
  for (int i=0; i < myVars.size(); i++){
    int loc = myVars.where(i);
    if (loc >=0){
     myVars[i] = active[loc];
     for (int i=0; i<(Nions+Nelec); i++){
        for (int nc = 1; nc < ncoeff;nc++){
           snap_beta[i][nc] = std::real(myVars[(i*ncoeff)+(nc-1)]);
        }
     }
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

