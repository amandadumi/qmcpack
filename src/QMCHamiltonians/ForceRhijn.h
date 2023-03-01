/** @file ForwardWalking.h
 * @brief Force estimator from https://pubs.acs.org/doi/10.1021/acs.jctc.1c00496
 */
#ifndef QMCPLUSPLUS_FORCE_RHIJN_HAMILTONIAN_H
#define QMCPLUSPLUS_FORCE_RHIJN_HAMILTONIAN_H
#include "QMCHamiltonians/ForceBase.h"
#include "QMCHamiltonians/OperatorBase.h"
#include "LongRange/LRCoulombSingleton.h"
#include "Numerics/OneDimGridBase.h"
#include "Numerics/OneDimGridFunctor.h"
#include "Numerics/OneDimCubicSpline.h"



namespace qmcplusplus
{
struct ForceRhijn: public OperatorBase, public ForceBase
{
public:
    int nstep;     // parameter: how many steps to consider in history
    /** constructor
    |*/
    ForceRhijn(ParticleSet& ions, ParticleSet& elns);
    ///destructor
    ~ForceRhijn();



    std::string getClassName() const override {return "ForceRhijn";}

    void resetTargetParticleSet(ParticleSet& P) override {}

    Return_t evaluate(ParticleSet& P) override;

    bool get(std::ostream& os) const override;

    std::unique_ptr<OperatorBase> makeClone(ParticleSet& qp, TrialWaveFunction& psi) final;

    bool put(xmlNodePtr cur) override;


    /**
     * implement evaluation of Pe^{-.5(E_L(R')+E_L(R)-E_T)}
    */
    void calculate_gdd(); 
    /**
     * implement evaluation of Pe^{(R'-R-2F)/2}
     * where F 2nablalnpsi v 
    */
    void calculate_gb();

private:
    const int d_aa_ID;
    const int d_ei_ID;
    std::vector<std::vector<int>> walker_lengths_;
};
}
#endif