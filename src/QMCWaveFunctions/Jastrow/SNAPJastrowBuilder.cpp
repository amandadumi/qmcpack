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


std::unique_ptr<WaveFunctionComponent> SNAPJastrowBuilder::buildComponent(xmlNodePtr cur){
  OhmmsAttributeSet oAttrib;
  oAttrib.add(TypeOpt, "type");
  oAttrib.add(NameOpt, "name");
  oAttrib.put(cur);
  return createSNAP(cur);
}

bool SNAPJastrowBuilder::putkids(xmlNodePtr kids)
{

while (kids != NULL){
    std::string kidsname = (char*)kids->name;
    // if this section is building in correlation...
    if (kidsname == "correlation")
    {

        //create strings to denote which species are set, by default elecs are both up.
        // iSpecies isn't set yet.
        std::string iSpecies, eSpecies1("u"), eSpecies2("u");
        // I don't know what an attribute set is, im guessing its a way to create objects in the code?
        // are we able to refer to them as the value in string later?
        //rAttrib.add(kids); // this was done in eeI_jastrow builder but can't get it to work here...
        // this seems to be if there are more than one set of coeffs given in this jastrow object. We will ned this eventually. We will ned this eventually.
        //const auto coef_id = extracCoefficientsID(kids); // We can survive without this for now since we will make sure to do I + u + d only
        std::string jname("SNAP");
        //functor->put(kids);
        // TODO: check wigner seitz radius vs jastrow cutoff
        //TODO: a lot of jastrow set up will happen here, initially just hard coding values to test.
        //functor->cutoff_radius = 1e-5;
        //SJ.addFunc(iNum,eNum1,eNum2,std::move(functor));
    }
    kids = kids->next;
    //SJ.check_complete();
    return true;

 }
 return false;
}

std::unique_ptr<WaveFunctionComponent> SNAPJastrowBuilder::createSNAP(xmlNodePtr cur)
{
ReportEngine PRE(ClassName, "createSNAP(xmlNodePtr)");
xmlNodePtr kids = cur->xmlChildrenNode;

//if ions are fed to jastrow 
std::string ftype("snap");
OhmmsAttributeSet tAttrib;
tAttrib.add(ftype, "function");
tAttrib.put(cur);

std::string input_name(getXMLAttributeValue(cur, "namae"));
std::string jname = input_name.empty() ? "snapjastrow" : input_name;
if (ftype == "snap"){
    auto SJ  = std::make_unique<SNAPJastrow>(ftype,sourcePtcl, targetPtcl);
    //SJ->setCoefficients(set); // eventually, currently can happen in constructor.
    //putkids(kids, *SJ);
    return SJ;
}
else
{
    std::ostringstream err_msg;
    err_msg << "Unknown function\"" << ftype << "\" in SNAPJastrowBuilder. Aborting.\n";
    APP_ABORT(err_msg.str());
}
return nullptr;
}


}