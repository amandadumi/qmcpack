//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////

#include "SNAPJastrowBuilder.h"
#include "SNAPJastrow.h"
#include "OhmmsData/AttributeSet.h"
#include "Utilities/ProgressReportEngine.h"

namespace qmcplusplus
{
SNAPJastrowBuilder::SNAPJastrowBuilder(Communicate* comm, ParticleSet& target, ParticleSet& source)
    : WaveFunctionComponentBuilder(comm, target)
{
    ClassName = "SnapJastrowBuilder";
    NameOpt = "0";
    TypeOpt = "SNAP";
    SourceOpt = SourcePtcl->getName();
}

std::unique_ptr<WaveFunctionComponent> SNAPJastrowBuilder::buildComponent(xmlNodePtr cur){
  OhmmsAttributeSet oAttrib;
  oAttrib.add(TypeOpt, "type");
  oAttrib.add(NameOpt, "name");
  oAttrib.put(cur);
  return createSNAP(cur);
}

bool SNAPJastrowBuilder::putkids(xmlNodePtr kids)
{
SpeciesSet& iSet = SourcePtcl->getSpeciesSet();
SpeciesSet& eSet = targetPtcl.getSpeciesSet();

while (kids != NULL){
    std::string kidsname = (char*)kids->name;
    if (kidsname == "correlation")
    {
        RealType eecusp = 0.0;
        RealType eIcusp = 0.0;
        std::string iSpecies, eSpecies1("u"), eSpecies2("u");
        OhmmsAttributeSet rAttrib;
        rAttrib.add(iSpecies,"ispecies");
        rAttrib.add(eSpecies1,"especies1");
        rAttrib.add(eSpecies2,"especies2");
        rAttrib.add(eecusp,"ecusp");
        rAttrib.add(eIcusp,"icusp");
        rAttrib.add(kids);
        const auto coef_id = extracCoefficientsID(kids);
        auto functor = 
        std::make_unique<SNAPJastrow>(coef_id.empty() ? jname + "_" + iSpecies + eSpecies1 + eSpecies2 : coef_id, ee_cusp, eI_cusp);
        functor->iSpecies = iSpecies;
        functor->eSpecies1 = eSpecies1;
        functor->eSpecies2 = eSpecies2;
        int iNum = iSet.findSpecies(iSpecies);
        int eNum1 = eSet.findSpecies(eSpecies1);
        int eNum2 = eSet.findSpecies(eSpecies2);
        if (iNum == iSet.size()){
            APP_ABORT("ion species " + iSpecies + " requested for Jastrow " + jname + " does not exist in ParticleSet")
        }
        std::string illegal_eSpeciess;
        if (eNum1 == eSet.size()){
            illegal_eSpecies = eSpecies1;
        }
        if (eNum2 == eSet.size()){
              if (illegal_eSpecies.size())
                illegal_eSpecies += " and ";
            illegal_eSpecies += eSpecies2;
        }
        if (illegal_eSpecies.size())
            APP_ABORT("electron species " + illegal_eSpecies + " requested for Jastrow " + jname + " does not exist in ParticleSet " + targetPtcl.getName());
    }
    functor->put(kids);
    // TODO: check wigner seitz radius vs jastrow cutoff
    functor->cutoff_radius = 1e-5;
    SJ.addFunc(iNum,eNum1,eNum2,std::move(functor));
    kids = kids->next;
    SJ.check_complete();

 }
}

std::unique_ptr<WaveFunctionComponent> SNAPJastrowBuilder::createSNAP(xmlNodePtr cur)
{
ReportEngine PRE(ClassName, "createSNAP(xmlNodePtr)");
xmlNodePtr kids = cur->xmlChildrenNode;

//if ions are fet to jastrow 
if (SourcePtcl){
    std::string ftype("snap");
    OhmmsAttributeSet tAttrib;
    tAttrib.add(ftype, "function");
    tAttrib.put(cur);

    std::string input_name(getXMLAttributeValue(cur, "namae"));
    std::string jname = input_name.empty() ? "snapjastrow" : input_name;
    SpeciesSet& iSet = SourcePtcl->getSpeciesSet();
    if (ftype == "snap"){
        using SNAPJ = SNAPJastrow();
        auto SJ           = std::make_unique<SNAPJ>;
        putkids(kids, *SJ);
        return SJ;
    }
    else
    {
        std::ostringstream err_msg;
        err_msg << "Unknown function\"" << ftype << "\" in SNAPJastrowBuilder. Aborting.\n";
        APP_ABORT(err_msg.str());
    }
}
else
    APP_ABORT("You must specify the \"source\" particleset for a three-body Jastrow.\n");
return nullptr;
}


}