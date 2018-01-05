#ifndef __HGCalSimulation_FastShower_Event_h__
#define __HGCalSimulation_FastShower_Event_h__


#include <unordered_map>

class Event
{

  private:
    Event() {}; // No default constructor, use Event(run, event)

  public:
    Event(uint32_t, uint32_t);
    ~Event() {};

    void fillHit(uint32_t, double);
    void fillGenEn(uint32_t, double);
    void fillGenEta(uint32_t, double);
    void fillGenPhi(uint32_t, double);
    void fillPDGid(uint32_t, int);
    void fillThickness(uint32_t, int);
    void setnPart(uint32_t);
    void clear();

    const uint32_t run() const {return run_;}
    const uint32_t event() const {return event_;}
    const uint32_t npart() const {return npart_;}
    const std::unordered_map<uint32_t, double>& hits() const {return hits_;}
    const std::unordered_map<uint32_t, double>& generatedEnergy() const {return generated_energy_;}
    const std::unordered_map<uint32_t, double>& generatedEta() const {return generated_eta_;}
    const std::unordered_map<uint32_t, double>& generatedPhi() const {return generated_phi_;}
    const std::unordered_map<uint32_t, double>& pdg_id() const {return pdgid_;}
    const std::unordered_map<uint32_t, int>& thickness() const {return thickness_;}

  private:
    uint32_t run_;
    uint32_t event_;
    uint32_t npart_;
    std::unordered_map<uint32_t, double> hits_;
    std::unordered_map<uint32_t, double> generated_energy_;
    std::unordered_map<uint32_t, double> generated_eta_;
    std::unordered_map<uint32_t, double> generated_phi_;
    std::unordered_map<uint32_t, double> pdgid_;
    std::unordered_map<uint32_t, int> thickness_;

};

#endif
