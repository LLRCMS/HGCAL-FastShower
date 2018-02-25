

#ifdef STANDALONE
#include "Event.h"
#else
#include "HGCalSimulation/FastShower/interface/Event.h"
#endif


Event::
Event(uint32_t run, uint32_t event):
  run_(run),
  event_(event)
{
}

void
Event::
fillHit(uint32_t id, double energy)
{
  auto itr_insert = hits_.emplace(id, energy);
  // if id already inserted (noise and hit), add the energy
  if(!itr_insert.second) itr_insert.first->second += energy;
}

void
Event::
fillGenEn(double energy)
{
  generated_energy_.push_back(energy);
}

void
Event::
fillGenEta(double eta)
{
  generated_eta_.push_back(eta);
}

void
Event::
fillGenPhi(double phi)
{
  generated_phi_.push_back(phi);
}

void
Event::
fillPDGid(int PDGid)
{
  pdgid_.push_back(PDGid);
}

void
Event::
fillThickness(uint32_t id, int thickness)
{
  thickness_.emplace(id, thickness);
}

void
Event::
setnPart(uint32_t part)
{
  npart_ = part;
}

void
Event::
clear()
{
  hits_.clear();
}
