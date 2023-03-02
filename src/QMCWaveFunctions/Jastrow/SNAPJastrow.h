#ifndef QMCPLUSPLUS_SNAPJASTROW
#define QMCPLUSPLUS_SNAPJASTROW
#include "QMCWaveFunctions/WaveFunctionComponent.h"
#include "Particle/DistanceTable.h"



namespace qmcplusplus
{
class SNAPJastrow : public WaveFunctionComponent
{

public:

    SNAPJastrow();

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
    PsiValueType ratio(ParticleSet& P, int iat) override;
    GradType evalGrad(ParticleSet& P, int iat) override;
    PsiValueType ratioGrad(ParticleSet& P, int iat, GradType& grad_iat) override;
    LogValueType evaluateGL(const ParticleSet& P,
                                  ParticleSet::ParticleGradient& G,
                                  ParticleSet::ParticleLaplacian& L,
                                  bool fromscratch);


/******Checkout related functons******/
    void registerData(ParticleSet& P, WFBufferType& buf) override;

    LogValueType updateBuffer(ParticleSet& P, WFBufferType& buf, bool fromscratch = false) override;

    void copyFromBuffer(ParticleSet& P, WFBufferType& buf) override;




void evaluateDerivatives(ParticleSet& P,
                                   const opt_variables_type& optvars,
                                   Vector<ValueType>& dlogpsi,
                                   Vector<ValueType>& dhpsioverpsi) override;

private:

};
}
#endif