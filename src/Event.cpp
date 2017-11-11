

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
fillGenEn(uint32_t id, double energy)
{
  generated_energy_.emplace(id, energy);
}

void
Event::
fillGenEta(uint32_t id, double eta)
{
  generated_eta_.emplace(id, eta);
}

void
Event::
fillGenPhi(uint32_t id, double phi)
{
  generated_phi_.emplace(id, phi);
}

void
Event::
fillPDGid(uint32_t id, int PDGid)
{
  pdgid_.emplace(id, PDGid);
}

void
Event::
fillThickness(uint32_t id, int thickness)
{
  thickness_.emplace(id, thickness);
}

void
Event::
fillHitCells(uint32_t id, const Cell& hit_cell)
{
  hit_cells_.emplace(id, hit_cell);
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


int
Event::
getLayerFromId(int cell_id)
{
  int len = std::to_string(cell_id).length();

  return atoi(std::to_string(cell_id).substr(len -2, len).c_str());
}
