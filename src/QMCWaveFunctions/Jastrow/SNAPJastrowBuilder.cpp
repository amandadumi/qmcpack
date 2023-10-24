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

inline bool putContent2(std::vector<double>& a, xmlNodePtr cur)
{
  std::istringstream stream(XMLNodeString{cur});
  double temp;
  a.clear();
  while (!stream.eof())
  {
    stream >> temp;
    if (stream.fail() || stream.bad())
      break;
    else
      a.push_back(temp);
  }
  return a.size() > 0;
} // adapted from kspace jastrow builder


std::unique_ptr<WaveFunctionComponent> SNAPJastrowBuilder::buildComponent(xmlNodePtr cur){
  OhmmsAttributeSet oAttrib;
  oAttrib.add(TypeOpt, "type");
  oAttrib.add(NameOpt, "name");
  oAttrib.put(cur);
  return createSNAP(cur);
}

bool SNAPJastrowBuilder::putkids(xmlNodePtr kids, SNAPJastrow& SJ)
{
    std::cout << "have entered into put kids function." <<std::endl;
    std::cout << "kisdname is " << (char*)kids->name <<std::endl;

while (kids != NULL){
    std::vector<double> snap_coeffs;
    std::vector<double>* coeffs;

    std::string kidsname = (char*)kids->name;
    std::cout << "kisdname is " << (char*)kids->name <<std::endl;

    // if this section is building in correlation...
    if (kidsname == "correlation")
    {

    std::cout << "within put kids correlation loop" <<std::endl;
    // I don't know what an attribute set is, im guessing its a way to create objects in the code?
    // are we able to refer to them as the value in string later?
    //rAttrib.add(kids); // this was done in eeI_jastrow builder but can't get it to work here...
    // this seems to be if there are more than one set of coeffs given in this jastrow object. We will ned this eventually. We will ned this eventually.
    //const auto coef_id = extracCoefficientsID(kids); // We can survive without this for now since we will make sure to do I + u + d only
    std::string jname("SNAP");
    std::string id_opt;
    if (coeffs)
      {
        xmlNodePtr xmlCoefs = kids->xmlChildrenNode;
        while (xmlCoefs != NULL)
        {
          std::string cname((const char*)xmlCoefs->name);
          if (cname == "coefficients")
          {
            std::string type("0"), id("0");
            OhmmsAttributeSet cAttrib;
            cAttrib.add(id_opt, "id");
            cAttrib.add(type, "type");
            cAttrib.put(xmlCoefs);
            coeffs = &snap_coeffs;
            if (type != "Array")
            {
              app_error() << "Unknown coefficients type "
                             ""
                          << type
                          << ""
                             " in SNAPJastrowBuilder.\n"
                          << "Resetting to "
                             "Array"
                             ".\n";
              xmlNewProp(xmlCoefs, (const xmlChar*)"type", (const xmlChar*)"Array");
            }
            //vector<T> can be read by this
            putContent2(*coeffs, xmlCoefs);
            app_log() << "  Read " << coeffs->size() << " coefficients for type " << type << std::endl;

            int this_id;
            // ## hard coded mapping for He 2 elec
            if (id_opt=="He"){
              this_id = 2;
            }
            if (id_opt=="C"){
              this_id = 2;
            }
            else if (id_opt=="eup")
            {
              this_id = 0 ;
            }
            else if (id_opt=="u")
            {
              this_id = 0 ;
            }
            else if (id_opt=="d")
            {
              this_id=1;
            }
            SJ.set_coefficients(*coeffs, this_id);
            std::cout << "done with set coeffs" <<std::endl;
          }
          xmlCoefs = xmlCoefs->next;
        }
      }
    }
    kids = kids->next;
} 
return true;
}

std::unique_ptr<WaveFunctionComponent> SNAPJastrowBuilder::createSNAP(xmlNodePtr cur)
{
ReportEngine PRE(ClassName, "createSNAP(xmlNodePtr)");
xmlNodePtr kids = cur->xmlChildrenNode;

//if ions are fed to jastrow 
int twojmax = 2;
std::string ftype("snap");
OhmmsAttributeSet tAttrib;
tAttrib.add(ftype, "function");
tAttrib.add(twojmax, "twojmax");
tAttrib.put(cur);

std::string input_name(getXMLAttributeValue(cur, "name"));
std::string jname = input_name.empty() ? "snapjastrow" : input_name;
std::cout << "creating jastrow using target: and source:" << targetPtcl.getName() << sourcePtcl.getName() <<std::endl;
if (ftype == "snap"){
    auto SJ  = std::make_unique<SNAPJastrow>(ftype, sourcePtcl, targetPtcl, twojmax);
    putkids(kids, *SJ);

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