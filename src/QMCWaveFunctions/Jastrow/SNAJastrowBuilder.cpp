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

#include "SNAJastrowBuilder.h"
#include "SNAJastrow.h"
#include "OhmmsData/AttributeSet.h"
#include "Utilities/ProgressReportEngine.h"

namespace qmcplusplus
{

inline bool putContent2(std::vector<double>& a, xmlNodePtr cur)
{
  std::cout << "SNAJastrowBuilder::putContent2 -> entered the function" << std::endl;
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


std::unique_ptr<WaveFunctionComponent> SNAJastrowBuilder::buildComponent(xmlNodePtr cur){
  OhmmsAttributeSet oAttrib;
  oAttrib.add(TypeOpt, "type");
  oAttrib.add(NameOpt, "name");
  oAttrib.put(cur);
  return createSNA(cur);
}

bool SNAJastrowBuilder::putkids(xmlNodePtr kids, SNAJastrow& SJ)
{

  SpeciesSet& iSet = jsource.getSpeciesSet();
  SpeciesSet& eSet = jtarget.getSpeciesSet();
  while (kids != NULL){
    std::vector<double> snap_coeffs;
    std::string kidsname = (char*)kids->name;

    // if this section is building in correlation...
    if (kidsname == "correlation")
    {
      std::string jname("SNA");
      std::string id_opt;
      xmlNodePtr xmlCoefs = kids->xmlChildrenNode;
      while (xmlCoefs != NULL)
      {
          std::string cname((const char*)xmlCoefs->name);
          if (cname == "coefficients")
          {
            std::cout<< "in coeff part of parsing"<< std::endl;
            std::string type("0"), id("0"), pset_id("0"), spec_id("0");
            OhmmsAttributeSet cAttrib;
            // which particle set do these coeffs correspond to (i.e. electrons of ions)
            cAttrib.add(pset_id, "pset_id");
            // which species coeffs correspond to since coeffs are stored by species type.
            cAttrib.add(id_opt, "id");
            cAttrib.add(type, "type");
            cAttrib.put(xmlCoefs);

            if (type != "Array")
            {
              app_error() << "Unknown coefficients type "
                             ""
                          << type
                          << ""
                             " in SNAJastrowBuilder.\n"
                          << "Resetting to "
                             "Array"
                             ".\n";
              xmlNewProp(xmlCoefs, (const xmlChar*)"type", (const xmlChar*)"Array");
            }
            //vector<T> can be read by this
            putContent2(snap_coeffs, xmlCoefs);
            app_log() << "  Read " << snap_coeffs.size() << " coefficients for type " << type << std::endl;
            std::cout << "SNAJastrowBuilder::putkids -> parsed coeffs" << std::endl;
            std::cout << "jtarget " << jtarget.getName() << " pset_id" << pset_id << std::endl;

            // ## hard coded mapping for He 2 elec

            
            int snap_beta_idx = 0;
            bool species_found;
            int num_e_species = eSet.size();
            int num_i_species = iSet.size();
            int cur_spec_id = eSet.findSpecies(id_opt);
            if (cur_spec_id < num_e_species){
              snap_beta_idx = cur_spec_id;
              species_found = true;
            }
            else if (iSet.findSpecies(id_opt)< num_i_species){
              int cur_spec_id = iSet.findSpecies(id_opt);
               std::cout << "found ion coeffs" << id_opt << std::endl;
              snap_beta_idx = cur_spec_id + eSet.size();
              species_found = true;
            }
            else{
               std::cout << "WARNING no species with label " << id_opt << std::endl;
            }
            if (species_found){
              SJ.set_coefficients(snap_coeffs, snap_beta_idx);
            }
            std::cout<< "SNAJastrowBuilder::pukids -> after set coefficients" <<std::endl;
          }
          xmlCoefs = xmlCoefs->next;
      }
    }
    kids = kids->next;
  }
  return true;
}

std::unique_ptr<WaveFunctionComponent> SNAJastrowBuilder::createSNA(xmlNodePtr cur)
{
ReportEngine PRE(ClassName, "createSNA(xmlNodePtr)");
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
    app_log() << "have created a jastrow object without any coefficients yet" <<std::endl;
    auto SJ  = std::make_unique<SNAJastrow>(ftype, jsource, jtarget, snap_type,twojmax,rcut);

    putkids(kids, *SJ);

    return SJ;
}
else
{
    std::ostringstream err_msg;
    err_msg << "Unknown function\"" << ftype << "\" in SNAJastrowBuilder. Aborting.\n";
    APP_ABORT(err_msg.str());
}
return nullptr;
}


}
