/** @file ForwardWalking.h
 * @brief Force estimator from https://pubs.acs.org/doi/10.1021/acs.jctc.1c00496
 */
#ifndef QMCPLUSPLUS_FORCERHIJN_H
#define QMCPLUSPLUS_FORCERHIJN_H
#include "QMCHamiltonians/OperatorBase.h"

namespace qmcplusplus
{
class QMCHamiltonian;

class ForceRhijn : public OperatorBase
{
public:
        // parameter: how many steps to consider in history
    /** constructor
    |*/
    ForceRhijn();
    ///destructor
    ~ForceRhijn() override;
    int nstep; 

    std::string getClassName() const override {return "ForceRhijn";}

    void resetTargetParticleSet(ParticleSet& P) override {}

    Return_t rejectedMove(ParticleSet& P) override;

    Return_t evaluate(ParticleSet& P) override;

    bool get(std::ostream& os) const override;

    std::unique_ptr<OperatorBase> makeClone(ParticleSet& qp, TrialWaveFunction& psi) final;

    bool put(xmlNodePtr cur) override;
    bool putSpecial(xmlNodePtr cur, QMCHamiltonian& h, ParticleSet& P);

    void addObservables(PropertySetType& plist);

    void addObservables(PropertySetType& plist, BufferType& collectables) override;

    void setObservables(PropertySetType& plist) override;

    void setParticlePropertyList(PropertySetType& plist, int offset) override;

private:
    std::vector<int> h_ids_; // im guessing this is hamiltonian ids since each observable becomes a member of the hamiltonian
    std::vector<int> p_ids_; // property id stored in history. 
    std::vector<int> walker_lengths_;
    std::vector<int> walker_tracker;
    std::vector<RealType> values_;

};
}
#endif