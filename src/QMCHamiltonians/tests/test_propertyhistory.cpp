#include "catch.hpp"

#include "OhmmsData/Libxml2Doc.h"
#include "OhmmsPETE/OhmmsMatrix.h"
#include "Particle/ParticleSet.h"
#include "ParticleIO/XMLParticleIO.h"
#include "QMCHamiltonians/ForceRhijn.h"
#include <stdio.h>
#include <string>
using std::string;
namespace qmcplusplus
{
using RealType = QMCTraits::RealType;
TEST_CASE("prophistory", "[hamiltonian]")
{
  Communicate* c = OHMMS::Controller;

  Libxml2Document doc;
  bool okay = doc.parse("Na2.structure.xml");
  REQUIRE(okay);
  xmlNodePtr root = doc.getRoot();

  const SimulationCell simulation_cell;

  ParticleSet ions(simulation_cell);
  XMLParticleParser parse_ions(ions);
  OhmmsXPathObject particleset_ion("//particleset[@name='ion0']", doc.getXPathContext());
  REQUIRE(particleset_ion.size() == 1);
  parse_ions.readXML(particleset_ion[0]);

  REQUIRE(ions.groups() == 1);
  REQUIRE(ions.R.size() == 2);
  ions.update();

  ParticleSet elec(simulation_cell);
  XMLParticleParser parse_elec(elec);
  OhmmsXPathObject particleset_elec("//particleset[@name='e']", doc.getXPathContext());
  REQUIRE(particleset_elec.size() == 1);
  parse_elec.readXML(particleset_elec[0]);
  REQUIRE(elec.groups() == 2);
  REQUIRE(elec.R.size() == 2);
  elec.addTable(ions);
  elec.update();

  ForceRihjn ph();
  // check that the number of times the walker has been visisted is equal to the last number recorded in property history.
 app_log() << "The first walker was progressed :  " << ph.walker_tracker[0] << std::endl;

  CHECK(ph.walker_tracker[0] == propertyhistory[0][ph.nsteps]);

}
}