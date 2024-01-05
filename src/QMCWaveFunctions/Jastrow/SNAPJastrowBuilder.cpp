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

  while (kids != NULL){
    std::vector<double> snap_coeffs;
    std::vector<double>* coeffs;

    std::string kidsname = (char*)kids->name;

    // if this section is building in correlation...
    if (kidsname == "correlation")
    {

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
            else if (id_opt=="C"){
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
            else if (id_opt=="ed"){
              this_id=1;
            }
            else if (id_opt=="edown"){
              this_id=1;
            }
            else{
              app_log()<< " particle type not recognized, correlation coefficient may not have been set as intended." <<std::endl;
            }
            SJ.set_coefficients(*coeffs, this_id);
            std::cout << "done with set coeffs" << std::endl;
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
double rcut=7;
std::string ftype("snap");
std::string snap_type("linear");
OhmmsAttributeSet tAttrib;
tAttrib.add(ftype, "function");
tAttrib.add(twojmax, "twojmax");
tAttrib.add(rcut, "rcut");
tAttrib.add(snap_type, "snap_type");
tAttrib.put(cur);

std::string input_name(getXMLAttributeValue(cur, "name"));
std::string jname = input_name.empty() ? "snapjastrow" : input_name;
if (ftype == "snap"){
    auto SJ  = std::make_unique<SNAPJastrow>(ftype, sourcePtcl, targetPtcl, snap_type,twojmax,rcut);

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
