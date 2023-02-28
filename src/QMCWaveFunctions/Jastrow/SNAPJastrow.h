#ifndef QMCPLUSPLUS_SNAPJASTROW
#define QMCPLUSPLUS_SNAPJASTROW
#include "Configuration.h"
#if !defined(QMC_BUILD_SANDBOX_ONLY)
#include "QMCWaveFunctions/WaveFunctionComponent.h"
#endif
#include "Particle/DistanceTable.h"
#include "CPU/SIMD/aligned_allocator.hpp"
#include "CPU/SIMD/algorithm.hpp"
#include <map>
#include <numeric>
#include <memory>



std::string std::string getClassName() const override {return "SNAPJastrow";}

/** Accpted move. Update Vat[iat],Grad[iat] and Lap[iat] */
void acceptMove(ParticleSet& P, int iat, bool safe_to_delay = false) override{
    UpdateMode = ORB_PBYP_RATIO;
    curAt      = computeU(P.getDistTableAB(myTableID).getTempDists());

};


inline void restore(int iat) override {}

PsiValueType ratio(ParticleSet& P, int iat) overrride {

}

void registerData(ParticleSet& P, WFBufferType& buf) override {}

LogValueType updateBuffer(ParticleSet& P, WFBufferType& buf, bool fromscratch = false) override {}

void copyFromBuffer(ParticleSet& P, WFBufferType& buf) override {}

/** Initialize a lammps object to get bispectrom components from current particle set.
| * 
| */
void initialize_lammps(){

}

/** From exsisting lammps object, get bispectrum components.
| * 
| */
void access_bispectrum(){


}


void evaluateDerivatives(ParticleSet& P,
                                   const opt_variables_type& optvars,
                                   Vector<ValueType>& dlogpsi,
                                   Vector<ValueType>& dhpsioverpsi) override {

                                   }