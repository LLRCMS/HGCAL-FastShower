#ifndef __HGCalSimulation_FastShower_Event_h__
#define __HGCalSimulation_FastShower_Event_h__

#include <vector>
#include <unordered_map>

class Event
{

  private:
    Event() {}; // No default constructor, use Event(run, event)

  public:
    Event(uint32_t, uint32_t);
    ~Event() {};

    void fillHit(uint32_t, double);
    void fillThickness(uint32_t, int);
    void fillGenEn(double);
    void fillGenEta(double);
    void fillGenPhi(double);
    void fillPDGid(int);
    void setnPart(uint32_t);
    void clear();

    const uint32_t run() const {return run_;}
    const uint32_t event() const {return event_;}
    const uint32_t npart() const {return npart_;}
    const std::unordered_map<uint32_t, double>& hits() const {return hits_;}
    const std::unordered_map<uint32_t, int>& thickness() const {return thickness_;}
    const std::vector<double>& generatedEnergy() const {return generated_energy_;}
    const std::vector<double>& generatedEta() const {return generated_eta_;}
    const std::vector<double>& generatedPhi() const {return generated_phi_;}
    const std::vector<double>& pdg_id() const {return pdgid_;}

  private:
    uint32_t run_;
    uint32_t event_;
    uint32_t npart_;
    std::unordered_map<uint32_t, double> hits_;
    std::unordered_map<uint32_t, int> thickness_;
    std::vector<double> generated_energy_;
    std::vector<double> generated_eta_;
    std::vector<double> generated_phi_;
    std::vector<double> pdgid_;

};

#endif
