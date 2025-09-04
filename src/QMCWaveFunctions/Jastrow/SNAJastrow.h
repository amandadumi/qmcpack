#ifndef QMCPLUSPLUS_SNAJASTROW
#define QMCPLUSPLUS_SNAJASTROW
#include "QMCWaveFunctions/WaveFunctionComponent.h"
#include "ParticleBase/ParticleAttrib.h"
#include "Particle/DistanceTable.h"
#include "Configuration.h"
#include <ResourceHandle.h>
#include "Message/MPIObjectBase.h"
#include "SNADesc.h"


namespace qmcplusplus
{

template<typename T>
struct SNAMultiWalkerMem;
class SNAJastrow : public WaveFunctionComponent, public OptimizableObject
{
public:

    using GradDerivVec  = ParticleAttrib<QTFull::GradType>;
    using ValueDerivVec = ParticleAttrib<QTFull::ValueType>;
    //handle d/dc info
    Vector<RealType> dLogPsi;
    std::vector<GradDerivVec> gradLogPsi;
    std::vector<ValueDerivVec> lapLogPsi;
    //handle per particle info
    RealType u_val;
    GradDerivVec grad_u;
    ValueDerivVec lap_u;
    SNAJastrow(const std::string& obj_name, ParticleSet& ions, ParticleSet& els, const std::string input_snap_type, int input_twojmax, double input_rcut);
    ~SNAJastrow();

    std::string getClassName() const override {return "SNAJastrow";}

    
    void resizeWFOptVectors(){
        dLogPsi.resize(myVars.size());
        gradLogPsi.resize(myVars.size(), GradDerivVec(Nelec));
        lapLogPsi.resize(myVars.size(), ValueDerivVec(Nelec));
    }

    /** Initialize a lammps object to get bispectrom components from current particle set.
    | * 
    | */
    void set_coefficients(std::vector<RealType>, int id);


/******MC step related functions******/

    /** Accpted move. Update Vat[iat],Grad[iat] and Lap[iat] */
    void acceptMove(ParticleSet& P, int iat, bool safe_to_delay = false) override; 

    void restore(int iat) override;

/******Optimization related functions******/
    bool isOptimizable() const override { return true; }
    /** check out optimizable variables
    */
    void checkOutVariables(const opt_variables_type& o) override;

    // void extractOptimizableObjectRefs(UniqueOptObjRefs& opt_obj_refs) override;

/****** Evaluate E_L functions ******/
    /** Calculate the ratio of proposed to current wave function element*/
    PsiValue ratio(ParticleSet& P, int iat) override;
    /** Calculate d/di U_SNap*/
    GradType evalGrad(ParticleSet& P, int iat) override;

    PsiValue ratioGrad(ParticleSet& P, int iat, GradType& grad_iat) override;
    
    LogValue evaluateGL(const ParticleSet& P,
                                  ParticleSet::ParticleGradient& G,
                                  ParticleSet::ParticleLaplacian& L,
                                  bool fromscratch) override;
    

    void computeGL(const ParticleSet& P,int iel);

    LogValue evaluateLog(const ParticleSet& P, ParticleSet::ParticleGradient& G, ParticleSet::ParticleLaplacian& L) override;
    
    
    void evaluateDerivatives(ParticleSet& P,
                                   const opt_variables_type& optvars,
                                   Vector<ValueType>& dlogpsi,
                                   Vector<ValueType>& dhpsioverpsi) override;

    void evaluateDerivativesWF(ParticleSet& P,
                                   const opt_variables_type& optvars,
                                   Vector<ValueType>& dlogpsi) override;

    /* calculates esnap based on a set of coefficients manually in qmcpack
    used to see impact of small change in coefficients on snap energy (needed to calculated d E/d beta)
    without having to internally change the lammps object.
    */
    void calculate_ESNA(const ParticleSet& P, const std::vector<std::vector<double>> new_coeff, double& new_u, bool proposed);
    void compute_bispectrum(int iat); 
    void compute_d_dr_bispectrum(int iat, std::vector<std::vector<double>>& a_snad);
    void update_sna_rij(const ParticleSet& P,int iat,bool proposed);
    void update_sna_rij_vp(const VirtualParticleSet& VP,int iat,int r);
    void evaluate_linear_derivs(ParticleSet& P, int coeff_idx);
    double FD_Lap(const ParticleSet& P,int iat, int dim, int coeff, int ntype);
    double full_FD_Lap(const ParticleSet& P,int iat, std::vector<std::vector<double>> coeffs);
    
    /****** NLPP-related functions ******/
    void evaluateRatios(const VirtualParticleSet& VP, std::vector<ValueType>& ratios) override;
    void evaluateDerivRatios(const VirtualParticleSet& VP,const opt_variables_type& optvars, std::vector<ValueType>& ratios, Matrix<ValueType>& dratios) override;
    /******Checkout-related functions ******/
    void registerData(ParticleSet& P, WFBufferType& buf) override;
    LogValue updateBuffer(ParticleSet& P, WFBufferType& buf, bool fromscratch = false) override;
    void copyFromBuffer(ParticleSet& P, WFBufferType& buf) override;
    void extractOptimizableObjectRefs(UniqueOptObjRefs& opt_obj_refs) override;
    void checkInVariablesExclusive(opt_variables_type& active) override;
    void resetParametersExclusive(const opt_variables_type& active) override;
    std::unique_ptr<WaveFunctionComponent> makeClone(ParticleSet& tpq) const override;
    bool put(xmlNodePtr cur);
    /***********Batched related options******************/
    void createResource(ResourceCollection& collection) const override;
    void acquireResource(ResourceCollection& collection,
                       const RefVectorWithLeader<WaveFunctionComponent>& wfc_list) const override;

    void releaseResource(ResourceCollection& collection,
                       const RefVectorWithLeader<WaveFunctionComponent>& wfc_list) const override;
    //variables
    const int Nions;
    const int Nelec;
    int NIonGroups;
    int ncoeff;
    int ntypes;
    int twojmax=2;
    int num_neigh; // keeps track of how many neighbors are in the current rij list
    double rcut=7;
    double rcutfac=1.0;
    double rfac0=0.99363;
    double rmin = 0;
    std::vector<std::vector<double>> rcutij;
    std::vector<int> type_map;
    std::vector<int> element;
    std::vector<double> radelem; 
    double dist_delta = 0.000000001;
    double coeff_delta = 0.0000001;
    const int ee_Table_ID_;
    const int ei_Table_ID_;
    int ii_Table_ID_;
    ParticleSet& Ions;
    std::string snap_type;
    double current_esnap=0;
    std::vector<std::vector<double>> snap_beta;
    std::vector<std::vector<double>> sna;
    std::vector<std::vector<double>> snad;
    // global arrays
    ResourceHandle<SNAMultiWalkerMem<RealType>> mw_mem_handle_;
    opt_variables_type myVars;
   // initialize the array that holds the bispectrum
   std::vector<double> bispectrum_components; //N_part x N_bispec
   // initialize the vector that holds the derivatives
   std::vector<std::vector<double>> ddr_bispectrum_components; //Nelec x N_type*N_bispec*N_dim
   SNADesc sna_desc;
  struct SNAJastrowTimers
  {
    NewTimer& eval_timer;
    NewTimer& init_lammps_timer;
    NewTimer& eval_esnap_timer;
    NewTimer& eval_log_timer;
    NewTimer& eval_ratio_timer;
    NewTimer& eval_gl_timer;
    NewTimer& eval_wf_grad_timer;
    NewTimer& eval_finite_diff_timer;
    SNAJastrowTimers(const std::string& prefix)
        : eval_timer(createGlobalTimer(prefix + "Eval", timer_level_fine)),
          init_lammps_timer(createGlobalTimer(prefix + "InitLammps", timer_level_fine)),
          eval_esnap_timer(createGlobalTimer(prefix + "evalESNA", timer_level_fine)),
          eval_log_timer(createGlobalTimer(prefix + "evalLogSNA", timer_level_fine)),
          eval_ratio_timer(createGlobalTimer(prefix + "evalRatioSNA", timer_level_fine)),
          eval_gl_timer(createGlobalTimer(prefix + "evalGLSNA", timer_level_fine)),
          eval_wf_grad_timer(createGlobalTimer(prefix + "evalWFGradSNA", timer_level_fine)),
          eval_finite_diff_timer(createGlobalTimer(prefix + "evalFiniteDiffSNA", timer_level_fine))
    {}
  }; 



  SNAJastrowTimers timers_;



 
};
}
#endif
