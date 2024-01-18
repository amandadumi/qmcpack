#ifndef QMCPLUSPLUS_SNAPJASTROW
#define QMCPLUSPLUS_SNAPJASTROW
#include "QMCWaveFunctions/WaveFunctionComponent.h"
#include "ParticleBase/ParticleAttribOps.h"
#include "Particle/DistanceTable.h"
#include "Configuration.h"
#include "Message/MPIObjectBase.h"
//lammps related headers needed.
#include "lammps.h"
#include "input.h"
#include "modify.h"
#include "compute.h"
#include "compute_snap.h"
#include "atom.h"
#include "pair.h"
#include "pair_snap.h"
#include "force.h"
#include "library.h"


namespace qmcplusplus
{
class SNAPJastrow : public WaveFunctionComponent, public OptimizableObject
{
public:

    using GradDerivVec  = ParticleAttrib<QTFull::GradType>;
    using ValueDerivVec = ParticleAttrib<QTFull::ValueType>;
    //handle d/dc info
    Vector<RealType> dLogPsi;
    std::vector<GradDerivVec> gradLogPsi;
    std::vector<ValueDerivVec> lapLogPsi;
    //handle per particle info
    Vector<RealType> u_val;
    GradDerivVec grad_u;
    ValueDerivVec lap_u;
    using PsiValueType = WaveFunctionComponent::PsiValueType;
    using LogValueType = WaveFunctionComponent::LogValueType;

    SNAPJastrow(const std::string& obj_name, const ParticleSet& ions, ParticleSet& els, const std::string input_snap_type, int input_twojmax,double input_rcut);

    ~SNAPJastrow();

    std::string getClassName() const override {return "SNAPJastrow";}

    
    void resizeWFOptVectors(){
        dLogPsi.resize(myVars.size());
        gradLogPsi.resize(myVars.size(), GradDerivVec(Nelec));
        lapLogPsi.resize(myVars.size(), ValueDerivVec(Nelec));
    }

    /** Initialize a lammps object to get bispectrom components from current particle set.
    | * 
    | */
    LAMMPS_NS::LAMMPS* initialize_lammps( const ParticleSet& els,double rcut);
    void set_coefficients(std::vector<RealType>,int id);


/******MC step related functions******/

    /** Accpted move. Update Vat[iat],Grad[iat] and Lap[iat] */
    void acceptMove(ParticleSet& P, int iat, bool safe_to_delay = false) override; 

    inline void restore(int iat) override {}
    /** From exsisting lammps object, get bispectrum components.
    | * 
    | */
    void access_bispectrum();


/******Optimization related functions******/
    bool isOptimizable() const override { return true; }
    /** check out optimizable variables
    */
    void checkOutVariables(const opt_variables_type& o) override;

    // void extractOptimizableObjectRefs(UniqueOptObjRefs& opt_obj_refs) override;

/****** Evaluate E_L functions ******/
    /** Calculate the ratio of proposed to current wave function element*/
    PsiValueType ratio(ParticleSet& P, int iat) override;
    /** Calculate d/di U_SNap*/
    GradType evalGrad(ParticleSet& P, int iat) override;

    PsiValueType ratioGrad(ParticleSet& P, int iat, GradType& grad_iat) override;
    
    LogValueType evaluateGL(const ParticleSet& P,
                                  ParticleSet::ParticleGradient& G,
                                  ParticleSet::ParticleLaplacian& L,
                                  bool fromscratch) override;
    

    void computeGL(const ParticleSet& P);

    LogValueType evaluateLog(const ParticleSet& P, ParticleSet::ParticleGradient& G, ParticleSet::ParticleLaplacian& L) override;
    
    
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
    void calculate_ESNAP(const ParticleSet& P, LAMMPS_NS::ComputeSnap* snap_global, std::vector<std::vector<double>> new_coeff, double& new_u,bool store_u);
    void calculate_ddc_gradlap_lammps(ParticleSet& P, double dist_delta, double coeff_delta,  std::vector<std::vector<double>>& fd_coeff, std::vector<std::vector<double>>& bd_coeff, int cur_val);
    void update_lmp_pos(const ParticleSet& P,LAMMPS_NS::LAMMPS* lmp_pntr, LAMMPS_NS::ComputeSnap* snap_array, int iat, bool proposed);
    void evaluate_fd_derivs(ParticleSet& P, int coeff_idx);
    void evaluate_linear_derivs(ParticleSet& P, int coeff_idx);
    double FD_Lap(const ParticleSet& P,int iat, int dim, int coeff, int ntype, std::vector<std::vector<double>> coeffs, double dist_delta, bool bispectrum_only);
    
    
    /****** NLPP-related functions ******/
    void evaluateRatios(const VirtualParticleSet& VP, std::vector<ValueType>& ratios) override;
    /******Checkout-related functions ******/
    void registerData(ParticleSet& P, WFBufferType& buf) override;
    LogValueType updateBuffer(ParticleSet& P, WFBufferType& buf, bool fromscratch = false) override;
    void copyFromBuffer(ParticleSet& P, WFBufferType& buf) override;
    void extractOptimizableObjectRefs(UniqueOptObjRefs& opt_obj_refs) override;
    void checkInVariablesExclusive(opt_variables_type& active) override;
    void resetParametersExclusive(const opt_variables_type& active) override;
    std::unique_ptr<WaveFunctionComponent> makeClone(ParticleSet& tpq) const override;
    bool put(xmlNodePtr cur);

    //variables
    const int Nions;
    const int Nelec;
    int NIonGroups;
    int ncoeff;
    int twojmax=2;
    double rcut=7;
    const int myTableID;
    const ParticleSet& Ions;
    std::string snap_type;
    std::vector<std::vector<double>> snap_beta;
    double hartree_over_ev = 1.000000589/27.211399998784;
    double bohr_over_ang = 1.88973; 
    // global arrays
    LAMMPS_NS::ComputeSnap* sna_global;
    LAMMPS_NS::ComputeSnap* proposed_sna_global;
    LAMMPS_NS::ComputeSnap* vp_sna_global;
    //lammps instance
    LAMMPS_NS::LAMMPS *lmp;
    LAMMPS_NS::LAMMPS *proposed_lmp;
    LAMMPS_NS::LAMMPS *vp_lmp;
    MPI_Comm comm_lammps;
    
    
    opt_variables_type myVars;

  struct SNAPJastrowTimers
  {
    NewTimer& eval_timer;
    NewTimer& init_lammps_timer;
    NewTimer& eval_esnap_timer;
    NewTimer& eval_log_timer;
    NewTimer& eval_ratio_timer;
    NewTimer& eval_gradient_timer;
    SNAPJastrowTimers(const std::string& prefix)
        : eval_timer(createGlobalTimer(prefix + "Eval", timer_level_fine)),
          eval_esnap_timer(createGlobalTimer(prefix + "evalESNAP", timer_level_fine)),
          init_lammps_timer(createGlobalTimer(prefix + "InitLammps", timer_level_fine)),
          eval_log_timer(createGlobalTimer(prefix + "evalLogSNAP", timer_level_fine)),
          eval_ratio_timer(createGlobalTimer(prefix + "evalRatioSNAP", timer_level_fine)),
          eval_gradient_timer(createGlobalTimer(prefix + "evalGradientSNAP", timer_level_fine))
    {}
  }; 



  SNAPJastrowTimers timers_;



 
};
}
#endif
