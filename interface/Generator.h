#ifndef __HGCalSimulation_FastShower_Generator_h__
#define __HGCalSimulation_FastShower_Generator_h__

#include <map>
#include <unordered_map>

#include "TH1F.h"
#include "TRandom3.h"
#include "TCanvas.h"

#ifdef STANDALONE
#include "Cell.h"
#include "Geometry.h"
#include "Parameters.h"
#include "ShowerParametrization.h"
#include "OutputService.h"
#else
#include "HGCalSimulation/FastShower/interface/Cell.h"
#include "HGCalSimulation/FastShower/interface/Geometry.h"
#include "HGCalSimulation/FastShower/interface/Parameters.h"
#include "HGCalSimulation/FastShower/interface/ShowerParametrization.h"
#include "HGCalSimulation/FastShower/interface/OutputService.h"
#endif


class Generator
{

  public:
    Generator(const Parameters&);
    ~Generator();

    void simulate();

    std::unique_ptr<TCanvas> display(const std::unordered_map<uint32_t,TH1F>& hCellEnergyEvtMap , std::vector<Cell>& , int ievt=0);

    std::array<std::array<double, 3>, 52> readCalibration(const std::string& filename);

  private:
    TRandom3 gun_;
    Geometry geometry_;
    OutputService output_;
    ShowerParametrization shower_;
    const Parameters& parameters_;
};


#endif
