#ifndef QMCPLUSPLUS_SNAPJASTROW
#define QMCPLUSPLUS_SNAPJASTROW
#include "QMCWaveFunctions/WaveFunctionComponent.h"
//lammps related headers needed.
#include "lammps.h"
#include "input.h"
#include "atom.h"
#include "pair.h"
#include "pair_snap.h"
#include "force.h"
#include "library.h"


namespace qmcplusplus
{
class SNAPJastrow : public WaveFunctionComponent
{
public:

    SNAPJastrow(const ParticleSet& ions, ParticleSet& els);
        ions_(ions),
        elecs_(els),
        my_table_ee_idx_(els.addTable(els, DTModes::NEED_TEMP_DATA_ON_HOST | DTModes::NEED_VP_FULL_TABLE_ON_HOST)),
        my_table_ei_idx_(els.addTable(ions, DTModes::NEED_VP_FULL_TABLE_ON_HOST)){};


    using PtclGrpIndexes   = QMCTraits::PtclGrpIndexes;

    ~SNAPJastrow();

    std::string getClassName() const override {return "SNAPJastrow";}

    /** Initialize a lammps object to get bispectrom components from current particle set.
    | * 
    | */
    void initialize_lammps();

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

    void extractOptimizableObjectRefs(UniqueOptObjRefs& opt_obj_refs) override;


/****** Evaluate E_L functions ******/
    /** Calculate the ratio of proposed to current wave function element*/
    PsiValueType ratio(ParticleSet& P, int iat) override;
    /** Calculate d/di U_SNap*/
    GradType evalGrad(ParticleSet& P, int iat) override;

    PsiValueType ratioGrad(ParticleSet& P, int iat, GradType& grad_iat) override;
    LogValueType evaluateGL(const ParticleSet& P,
                                  ParticleSet::ParticleGradient& G,
                                  ParticleSet::ParticleLaplacian& L,
                                  bool fromscratch);

    LogValueType evaluateLog(const ParticleSet& P,
    ParticleSet::ParticleGradient& G,
    ParticleSet::ParticleLaplacian& L);

    void bispectrum_Laplacian_finite_diff();

/******Checkout related functons******/
    void registerData(ParticleSet& P, WFBufferType& buf) override;

    LogValueType updateBuffer(ParticleSet& P, WFBufferType& buf, bool fromscratch = false) override;

    void copyFromBuffer(ParticleSet& P, WFBufferType& buf) override;




void evaluateDerivatives(ParticleSet& P,
                                   const opt_variables_type& optvars,
                                   Vector<ValueType>& dlogpsi,
                                   Vector<ValueType>& dhpsioverpsi) override;

private:
 LAMMPS_NS::LAMMPS *lmp;                  // pointer to lammps object
 void* pointer_to_bispectrum_coeff;         // pointer to bispectrum coeffs.
 void** pointer_to_dbidrj;                 // pointer to location in lammps object for derivative of bispectrum components with position
 
};
}
#endif