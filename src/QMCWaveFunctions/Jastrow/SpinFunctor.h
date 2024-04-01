//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: John R. Gergely,  University of Illinois at Urbana-Champaign
//                    Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//                    Amrita Mathuriya, amrita.mathuriya@intel.com, Intel Corp.
//                    Peter W. Doak, doakpw@ornl.gov, Oak Ridge National Lab
//
// File created by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_SPIN_FUNCTOR_H
#define QMCPLUSPLUS_SPIN_FUNCTOR_H

#include <cstdio>
#include <memory>
#include "OptimizableFunctorBase.h"
#include "Utilities/ProgressReportEngine.h"
#include "OhmmsData/AttributeSet.h"
#include "OhmmsPETE/OhmmsVector.h"
#include "Numerics/LinearFit.h"
#include "OMPTarget/OffloadAlignedAllocators.hpp"
#include "Numerics/SplineBound.hpp"

namespace qmcplusplus
{

/**BsplineFunctor class for the Jastrows
 * REAL is the real type used by offload target, it is the correct type for the mw data pointers
 * and is also used to coerce/implicitly convert the Real type inherited OptimizableFunctorBase into that buffer
 * if offload is off this happens too but is just an implementation quirk.
 */
template<typename REAL>
struct BsplineFunctor : public OptimizableFunctorBase
{
  using Real = OptimizableFunctorBase::real_type;



  int NumParams;
  Real DeltaR, DeltaRInv;
  Real CuspValue;
  Real Y, dY, d2Y;
  // Stores the derivatives w.r.t. coefs
  // of the u, du/dr, and d2u/dr2
  std::vector<Real> Parameters;
  std::vector<std::string> ParameterNames;
  std::string elementType, pairType;
  std::string fileName;

  bool notOpt;
  bool periodic;

  ///constructor
  SpinFunctor(const std::string& my_name, Real cusp = 0.0)
      : OptimizableFunctorBase(my_name), NumParams(0), CuspValue(cusp), notOpt(false), periodic(true)
  {
    cutoff_radius = 0.0;
  }

  OptimizableFunctorBase* makeClone() const override { return new BsplineFunctor(*this); }

  void setCusp(Real c) override { CuspValue = c; }

  void setPeriodic(bool p) override { periodic = p; }

  /// return the max allowed beginning index to access spline coefficients
  int getMaxIndex() const { return spline_coefs_->size() - 4; }



  /** reset coefs from Parameters
   */
  void reset() override
  {
    const int numCoefs = NumParams + 4;
    const int numKnots = numCoefs - 2;
    DeltaR             = cutoff_radius / (Real)(numKnots - 1);
    DeltaRInv          = 1.0 / DeltaR;
    auto& coefs        = *spline_coefs_;
    for (int i = 0; i < coefs.size(); i++)
      coefs[i] = 0.0;
    // Ensure that cusp conditions is satisfied at the origin
    coefs[1] = Parameters[0];
    coefs[2] = Parameters[1];
    coefs[0] = Parameters[1] - 2.0 * DeltaR * CuspValue;
    for (int i = 2; i < Parameters.size(); i++)
      coefs[i + 1] = Parameters[i];
    coefs.updateTo();
  }

  /** compute value, first and second derivatives for [iStart, iEnd) pairs
   * @param iat the source particle that should be avoided (self pairs)
   * @param iStart starting particle index
   * @param iEnd ending particle index
   * @param _distArray distance arrUay
   * @param _valArray  u(r_j) for j=[iStart,iEnd)
   * @param _gradArray  du(r_j)/dr /r_j for j=[iStart,iEnd)
   * @param _lapArray  d2u(r_j)/dr2 for j=[iStart,iEnd)
   * @param distArrayCompressed temp storage to filter r_j < cutoff_radius
   * @param distIndices temp storage for the compressed index
   */
  void evaluateVGL(const int iat,
                   const int iStart,
                   const int iEnd,
                   const REAL* _distArray,
                   REAL* restrict _valArray,
                   REAL* restrict _gradArray,
                   REAL* restrict _laplArray,
                   REAL* restrict distArrayCompressed,
                   int* restrict distIndices) const;

  /** compute value, gradient and laplacian for target particles
   * This more than just a batched call of evaluateVGL
   * @param iat the source particle that should be avoided (self pairs)
   * @param num_groups the number of source particle groups
   * @param functors for the num_groups of source particles
   * @param n_src the number of source particles
   * @param grp_ids the group ids of the n_src source particles
   * @param nw batch size (number of walkers)
   * @param mw_vgl return resutls. Multi walker value, gradient and laplacian [nw][1(v)+DIM(g)+1(l)]
   * @param n_padded the padded size of source particles
   * @param mw_dist Multi walker distance table [nw][1(distance)+DIM(displacements)][n_padded]
   * @param mw_cur_allu Multi walker value, first and second derivatives of pair potentials [nw][DIM][n_padded]. if mw_cur_allu is dual space, only update device side.
   * @param transfer_buffer temporary transfer buffer.
   *
   * If mw_dist is dual space, up-to-date data is assumed on device.
   * If mw_cur_allu is dual space, data is created on the device and there is no transfer to the host
   * because it will be consumed by mw_updateVGL on the device.
   */
  static void mw_evaluateVGL(const int iat,
                             const int num_groups,
                             const BsplineFunctor* const functors[],
                             const int n_src,
                             const int* grp_ids,
                             const int nw,
                             REAL* mw_vgl, // [nw][DIM+2]
                             const int n_padded,
                             const REAL* mw_dist, // [nw][DIM+1][n_padded]
                             REAL* mw_cur_allu,   // [nw][3][n_padded]
                             Vector<char, OffloadPinnedAllocator<char>>& transfer_buffer);

  /** evaluate sum of the pair potentials for [iStart,iEnd)
   * @param iat dummy
   * @param iStart starting particle index
   * @param iEnd ending particle index
   * @param _distArray distance arrUay
   * @param distArrayCompressed temp storage to filter r_j < cutoff_radius
   * @return \f$\sum u(r_j)\f$ for r_j < cutoff_radius
   */
  REAL evaluateV(const int iat,
                 const int iStart,
                 const int iEnd,
                 const REAL* restrict _distArray,
                 REAL* restrict distArrayCompressed) const;

  /** compute value for target-source particle pair potentials
   * This more than just a batched call of evaluateV
   * @param num_groups the number of source particle groups
   * @param functors for the num_groups of source particles
   * @param n_src the number of source particles
   * @param grp_ids the group ids of the n_src source particles
   * @param nnum_pairs the number of particle pairs
   * @param ref_at the source particles that should be avoided (self pairs)
   * @param mw_vgl return resutls. Multi walker value, gradient and laplacian [nw][1(v)+DIM(g)+1(l)]
   * @param dist_stride the offset of distance pointers between to consecutive walkers
   * @param mw_dist Multi walker distance table [nw][1(distance)+DIM(displacements)][n_padded]
   * @param transfer_buffer temporary transfer buffer.
   *
   * If mw_dist is dual space, up-to-date data is assumed on device.
   */
  static void mw_evaluateV(const int num_groups,
                           const BsplineFunctor* const functors[],
                           const int n_src,
                           const int* grp_ids,
                           const int num_pairs,
                           const int* ref_at,
                           const REAL* mw_dist,
                           const int dist_stride,
                           REAL* mw_vals,
                           Vector<char, OffloadPinnedAllocator<char>>& transfer_buffer);


  inline Real evaluate(Real r) const
  {
    Real u(0);
    if (r < cutoff_radius)
      u = evaluate_impl(r, spline_coefs_->data(), DeltaRInv, getMaxIndex());
    return u;
  }

  inline bool evaluateDerivatives(Real r, std::vector<TinyVector<Real, 3>>& derivs) override
  
    return true;
  }

  inline bool evaluateDerivatives(Real r, std::vector<Real>& derivs) override
  {
    if (r >= cutoff_radius)
      return false;
    r *= DeltaRInv;

    return true;
  }

  inline Real f(Real r) override
  {
    if (r >= cutoff_radius)
      return 0.0;
    return evaluate(r);
  }

  inline Real df(Real r) override
  {
    if (r >= cutoff_radius)
      return 0.0;
    Real du, d2u;
    evaluate(r, du, d2u);
    return du;
  }


  bool put(xmlNodePtr cur) override
  {
    ReportEngine PRE("SpinFunctor", "put(xmlNodePtr)");
    //CuspValue = -1.0e10;
    NumParams = 0;
    //cutoff_radius = 0.0;
    OhmmsAttributeSet rAttrib;
    Real radius = -1.0;
    rAttrib.add(NumParams, "size");
    rAttrib.add(radius, "rcut");
    rAttrib.add(radius, "cutoff");
    rAttrib.put(cur);
    if (radius < 0.0)
      if (periodic)
      {
        app_log() << "    Jastrow cutoff unspecified.  Setting to Wigner-Seitz radius = " << cutoff_radius << std::endl;
        app_log() << std::endl;
      }
      else
      {
        APP_ABORT("  Jastrow cutoff unspecified.  Cutoff must be given when using open boundary conditions");
      }
    else if (periodic && radius > cutoff_radius)
    {
      if (radius - cutoff_radius > 1e-4)
      {
        APP_ABORT("  The Jastrow cutoff specified should not be larger than Wigner-Seitz radius.");
      }
      else
      {
        app_log() << "  The Jastrow cutoff specified is slightly larger than the Wigner-Seitz radius.";
        app_log() << "  Setting to Wigner-Seitz radius = " << cutoff_radius << ".\n";
      }
    }
    else
      cutoff_radius = radius;
    if (NumParams == 0)
    {
      PRE.error("You must specify a positive number of parameters for the Bspline jastrow function.", true);
    }
    app_summary() << "     Number of parameters: " << NumParams << std::endl;
    app_summary() << "     Cusp: " << CuspValue << std::endl;
    app_summary() << "     Cutoff radius: " << cutoff_radius << std::endl;
    resize(NumParams);
    // Now read coefficents
    xmlNodePtr xmlCoefs = cur->xmlChildrenNode;
    while (xmlCoefs != NULL)
    {
      std::string cname((const char*)xmlCoefs->name);
      if (cname == "coefficients")
      {
        std::string type("0"), id("0");
        std::string optimize("yes");
        OhmmsAttributeSet cAttrib;
        cAttrib.add(id, "id");
        cAttrib.add(type, "type");
        cAttrib.add(optimize, "optimize");
        cAttrib.put(xmlCoefs);
        if (type != "Array")
        {
          PRE.error("Unknown correlation type " + type + " in BsplineFunctor." + "Resetting to \"Array\"");
          xmlNewProp(xmlCoefs, (const xmlChar*)"type", (const xmlChar*)"Array");
        }
        std::vector<Real> params;
        putContent(params, xmlCoefs);
        if (params.size() == NumParams)
          Parameters = params;
        else
        {
          app_log() << "    Changing number of Bspline parameters from " << params.size() << " to " << NumParams
                    << ".  Performing fit:\n";
          // Fit function to new number of parameters
          const int numPoints = 500;
          SpinFunctor<REAL> tmp_func("tmp_func", CuspValue);
          tmp_func.cutoff_radius = cutoff_radius;
          tmp_func.resize(params.size());
          tmp_func.Parameters = params;
          tmp_func.reset();
          std::vector<Real> y(numPoints);
          Matrix<Real> basis(numPoints, NumParams);
          std::vector<TinyVector<Real, 3>> derivs(NumParams);
          for (int i = 0; i < numPoints; i++)
          {
            Real r = (Real)i / (Real)numPoints * cutoff_radius;
            y[i]   = tmp_func.evaluate(r);
            evaluateDerivatives(r, derivs);
            for (int j = 0; j < NumParams; j++)
              basis(i, j) = derivs[j][0];
          }
          resize(NumParams);
          LinearFit(y, basis, Parameters);
          app_log() << "New parameters are:\n";
          for (int i = 0; i < Parameters.size(); i++)
            app_log() << "   " << Parameters[i] << std::endl;
        }
        if (optimize == "yes")
        {
          notOpt = false;
        }
        else
        {
          notOpt = true;
        }
        for (int i = 0; i < NumParams; i++)
        {
          std::stringstream sstr;
          sstr << id << "_" << i;
          myVars.insert(sstr.str(), (Real)Parameters[i], !notOpt, optimize::LOGLINEAR_P);
        }
        int left_pad_space = 5;
        app_log() << std::endl;
        myVars.print(app_log(), left_pad_space, true);
      }
      xmlCoefs = xmlCoefs->next;
    }
    reset();
    Real zeros = 0;
    for (int i = 0; i < NumParams; i++)
      zeros += Parameters[i] * Parameters[i];
    return zeros > 1.0e-12; //true if Parameters are not zero
  }

  void initialize(int numPoints,
                  std::vector<Real>& x,
                  std::vector<Real>& y,
                  REAL cusp,
                  REAL rcut,
                  std::string& id,
                  std::string& optimize)
  {
    ReportEngine PRE("Spin", "initialize");
    NumParams     = numPoints;
    cutoff_radius = rcut;
    CuspValue     = cusp;
    
#if !defined(QMC_BUILD_SANDBOX_ONLY)
    if (optimize == "yes")
    {
      // Setup parameter names
      for (int i = 0; i < NumParams; i++)
      {
        std::stringstream sstr;
        sstr << id << "_" << i;
        myVars.insert(sstr.str(), (Real)Parameters[i], true, optimize::LOGLINEAR_P);
      }
      myVars.print(app_log());
    }
    else
#endif
    {
      notOpt = true;
      app_log() << "Parameters of BsplineFunctor id:" << id << " are not being optimized.\n";
    }
    reset();
  }

  void reportStatus(std::ostream& os) override
  {
    if (notOpt)
      return;
    myVars.print(os);
  }

  void checkOutVariables(const opt_variables_type& active) override
  {
    if (notOpt)
      return;
    myVars.getIndex(active);
  }

  void checkInVariablesExclusive(opt_variables_type& active) override
  {
    if (notOpt)
      return;
    myVars.setIndexDefault();
    active.insertFrom(myVars);
  }

  void resetParametersExclusive(const opt_variables_type& active) override
  {
    if (notOpt)
      return;
    for (int i = 0; i < Parameters.size(); ++i)
    {
      int loc = myVars.where(i);
      if (loc >= 0)
        Parameters[i] = std::real(myVars[i] = active[loc]);
    }
    reset();
  }

  // check if this object has active optimizable parameters
  bool isOptimizable()
  {
    if (notOpt)
      return false;
    for (int i = 0; i < Parameters.size(); ++i)
    {
      int loc = myVars.where(i);
      if (loc >= 0)
        return true;
    }
    return false;
  }
};

template<typename REAL>
inline REAL SpinFunctor<REAL>::evaluateV(const int iat,
                                            const int iStart,
                                            const int iEnd,
                                            const REAL* restrict _distArray,
                                            REAL* restrict distArrayCompressed) const
{
  int d =0
  return d;
}

template<typename REAL>
inline void SpinFunctor<REAL>::evaluateVGL(const int iat,
                                              const int iStart,
                                              const int iEnd,
                                              const REAL* _distArray,
                                              REAL* restrict _valArray,
                                              REAL* restrict _gradArray,
                                              REAL* restrict _laplArray,
                                              REAL* restrict distArrayCompressed,
                                              int* restrict distIndices) const
{
 
  }
}

extern template struct SplineFunctor<QMCTraits::RealType>;

} // namespace qmcplusplus
#endif
