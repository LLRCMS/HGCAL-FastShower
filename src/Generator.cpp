#include <iostream>

#include "TStopwatch.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TFile.h"
#include "TText.h"
#include "TPaveText.h"
#include <typeinfo>

#include "TMarker.h"
#include <fstream>
#include <string>
#include <algorithm>

#ifdef STANDALONE
#include "Generator.h"
#include "ShowerShapeHexagon.h"
#include "ShowerShapeTriangle.h"
#include "Event.h"
#include "Tree.h"
#else
#include "HGCalSimulation/FastShower/interface/Generator.h"
#include "HGCalSimulation/FastShower/interface/ShowerParametrization.h"
#include "HGCalSimulation/FastShower/interface/ShowerShapeHexagon.h"
#include "HGCalSimulation/FastShower/interface/ShowerShapeTriangle.h"
#include "HGCalSimulation/FastShower/interface/Event.h"
#endif

Generator::
Generator(const Parameters& params):
  geometry_(params.geometry()),
  output_(params.general().output_file),
  shower_(params.shower()),
  parameters_(params)
{
}


Generator::~Generator(){
}


std::array<std::array<double, NB_SI_THICKNESS>, NB_LAYERS> Generator::readCalibration(const std::string& filename) {

  std::array<std::array<double, NB_SI_THICKNESS>, NB_LAYERS> calib;

  std::ifstream my_calib_file(filename);
  std::stringstream ss;
  std::string line, temp;

  int col_counter = 0;
  int line_counter = 0;
  // There is one col for the layer id. MIP and noise have both same number of column as the
  // number of silicium thickness
  int dim_column = 1 + NB_SI_THICKNESS * 2;

  if (my_calib_file) {

    int layer_;
    double mip[NB_SI_THICKNESS];
    double noise[NB_SI_THICKNESS];
    double sampl, MIP;

    while(std::getline(my_calib_file, line)) {

      // Remove comment lines (beginning with #)
      if (line.find("#") != std::string::npos) {
        std::cout<<"This line begins with #. It is considered as comment -> discarded"<<std::endl;
        continue;
      }
      else {
        // Some checks for the file format : columns and lines
        if (col_counter == 0) {
          ss.clear();
          ss << line;
          while (ss >> temp) {
            col_counter++;
          }
        }
        line_counter++;
      }

      istringstream strm(line);

      strm >> layer_;
      strm >> mip[2];
      strm >> mip[1];
      strm >> mip[0];
      strm >> noise[2];
      strm >> noise[1];
      strm >> noise[0];

      for (int i = 0; i < NB_SI_THICKNESS; i++) {
        sampl = mip[i];

        if (layer_ <= parameters_.geometry().FH_limit_layer) {
          MIP = sampl/100;
        }
        else {
          MIP = sampl/10;
        }
        calib[layer_-1][i] = noise[i]*MIP/sampl;
      }
    }
  }
  else {
  std::cout << "Wrong file path or the file does not exit"<< std::endl;
  }

  if (col_counter != dim_column && line_counter != NB_LAYERS){
    std::cout << "WARNING : check the format of the calibration file"<< std::endl;
    std::cout << "You should have 7 columns (layer, 3 MIP and 3 noise for the different Si thickness)"<<std::endl;
    std::cout << "and one line per layer."<<std::endl;
  }

  my_calib_file.close();

  return calib;
}


void Generator::simulate() {
  unsigned nevents = parameters_.general().events;

  std::unordered_map<uint32_t, TH1F> hCellEnergyMap;
  std::unordered_map<uint32_t, TH1F> hCellEnergyEvtMap;

  TStopwatch t;
  t.Start();

  // some initializations
  double energygen=0.;
  double energygenincells=0.;
  double energyrec=0.;

  // Geometry
  std::cout<<"I'm building the geometry, please wait..." <<std::endl;

  // Shower parametrization
  ShowerParametrization aShowerParametrization(parameters_.shower());

  // Canvas for the geometry map in rootfile
  std::vector<std::unique_ptr<TCanvas>> canvas;

  // Build geometry
  int layer_min, layer_max;
  int display_layer;

  if (parameters_.geometry().layer==-1) {
      layer_min = 0;
      layer_max = parameters_.geometry().layers_z.size();
      display_layer = parameters_.display().layer;
  }
  else {
    layer_min = parameters_.geometry().layer;
    layer_max = parameters_.geometry().layer + 1;
    display_layer = layer_min;
  }

  for (int layer_id = layer_min; layer_id < layer_max; layer_id++) {

    std::cout << "=====> Layer " <<layer_id + 1 << " / "<< layer_max << '\r' << flush;

    if (parameters_.geometry().type!=Parameters::Geometry::Type::External) {
      geometry_.constructFromParameters(parameters_.general().debug, layer_id, display_layer);
    }
    else{
      geometry_.constructFromJson(parameters_.general().debug, layer_id);
    }

    if (display_layer == layer_id) {
      std::string hName;
      hName = "geometry_";
      hName += std::to_string(display_layer + 1);
      TH2Poly* geometry_histo = (TH2Poly*)geometry_.cellHistogram()->Clone(hName.c_str());
      geometry_histo->Write();
      delete geometry_histo;

      for (const auto& c : geometry_.getCells()) {
        if (c.second.getLayer() == display_layer) {
          int i = c.second.getIIndex();
          int j = c.second.getJIndex();

          hName="hCellEnergy_[";
          hName += std::to_string(i);
          hName += ",";
          hName += std::to_string(j);
          hName += ",";
          hName += std::to_string(display_layer);
          hName += "]";
          hCellEnergyMap.emplace(c.first, TH1F(hName.c_str(),"Energy in cell [i,j,k])",100,0.,100.));
        }
      }

      if (parameters_.display().events > 0) {
        for (const auto& c : geometry_.getCells()) {
          if (c.second.getLayer() == display_layer) {
            int i = c.second.getIIndex();
            int j = c.second.getJIndex();

            hName="hCellEnergyEvt[";
            hName += std::to_string(i);
            hName += ",";
            hName += std::to_string(j);
            hName += ",";
            hName += std::to_string(display_layer);
            hName += "]";
            hCellEnergyEvtMap.emplace(c.first,TH1F(hName.c_str(),"Event Energy in cell [i,j,k])", 100, 0., 100.));
          }
        }
      }
    }
  }

  if (parameters_.general().debug) {
    geometry_.print();
  }

  std::cout << "Done !"<<std::endl;

  // Noise calibration of each cells for all layers
  std::array<std::array<double, NB_SI_THICKNESS>, NB_LAYERS> calibratednoise;

  if (parameters_.generation().noise) {
    std::cout << "Noise calibration ..."<<std::endl;

    if (parameters_.generation().gentype == Parameters::Generation::GenType::Personnal) {

      for(int layer_id = layer_min; layer_id < layer_max; layer_id++) {

        double sigma_noise = parameters_.generation().noise_sigma;
        double mip;
        double sampl = parameters_.generation().sampling;

        if (layer_id <= parameters_.geometry().FH_limit_layer) {
          mip = sampl/100;
        }
        else {
          mip = sampl/10;
        }

        for (int i = 0; i<NB_SI_THICKNESS; i++) {
          calibratednoise[layer_id][i] = sigma_noise*mip/sampl;
        }
      }
    }
    else {
      calibratednoise = readCalibration(parameters_.generation().calib_file);
    }
  }
  std::cout << "Done !"<<std::endl;

  // start main loop on all events
  double hit_outside_geom = 0.;
  double tot_hits = 0;

  // // Histograms for transverse profile
  // TH1F *hTransverseProfile[layer_max-layer_min];
  // char* hTransProfName = new char[30];

  // for (int i = layer_min; i < layer_max; i++) {
  //   sprintf(hTransProfName, "hTransverseProfile_%d", i);
  //   hTransverseProfile[i] = new TH1F(hTransProfName, "Generated transverse profile (cm)", 600, 0., 120);
  // }

  // Build the tree map (except for external geometry) ; one key corresponds to one layer
  std::unordered_map<int, Tree*> tree_map;

  if (parameters_.geometry().type!=Parameters::Geometry::Type::External) {

    for (int layer_id = layer_min; layer_id < layer_max; layer_id++) {

      std::vector<double > cell_x;
      std::vector<double > cell_y;
      for (const auto& c : geometry_.getCells()) {
        if (c.second.getLayer() == layer_id) {
          cell_x.push_back(c.second.getX());
          cell_y.push_back(c.second.getY());
        }
      }

      double x_min = *std::minmax_element(cell_x.begin(), cell_x.end()).first;
      double x_max = *std::minmax_element(cell_x.begin(), cell_x.end()).second;
      double y_min = *std::minmax_element(cell_y.begin(), cell_y.end()).first;
      double y_max = *std::minmax_element(cell_y.begin(), cell_y.end()).second;

      // Geometry build with the cell center position, add a margin for the
      // square corners containing the boarder cells
      x_min -= parameters_.geometry().small_cell_side;
      x_max += parameters_.geometry().large_cell_side;
      y_min -= parameters_.geometry().large_cell_side * 2;
      y_max += parameters_.geometry().large_cell_side * 2;

      Rectangle plan = Rectangle(
        Point((float)x_min, (float)y_max),
        Point((float)x_max, (float)y_min)
      );

      Tree* tree = new Tree(plan, 5);
      float x_c, y_c;
      double side;
      for (const auto& c : geometry_.getCells()) {

        if (c.second.getLayer() == layer_id) {

          x_c = float(c.second.getX());
          y_c = float(c.second.getY());
          float cell_radius = sqrt((x_c*x_c)+(y_c*y_c));

          if (cell_radius <= float(parameters_.geometry().limit_first_zone)) {
            side = parameters_.geometry().small_cell_side;
          }
          else {
            side = parameters_.geometry().large_cell_side;
          }

          Point cornerA = Point(x_c - float(side*sqrt(3)/2), y_c + float(side));
          Point cornerB = Point(x_c + float(side*sqrt(3)/2), y_c + float(side));
          Point cornerC = Point(x_c + float(side*sqrt(3)/2), y_c - float(side));
          Point cornerD = Point(x_c - float(side*sqrt(3)/2), y_c - float(side));

          std::set<Tree*> alreadyAdded;

          Tree* leafA = tree->getLeaf(cornerA);
          leafA->addCell(&(c.second));
          alreadyAdded.insert(leafA);

          Tree* leafB = tree->getLeaf(cornerB);

          if(alreadyAdded.count(leafB) == 0) {
            leafB->addCell(&(c.second));
            alreadyAdded.insert(leafB);
          }

          Tree* leafC = tree->getLeaf(cornerC);

          if(alreadyAdded.count(leafC) == 0) {
            leafC->addCell(&(c.second));
            alreadyAdded.insert(leafC);
          }

          Tree* leafD = tree->getLeaf(cornerD);

          if(alreadyAdded.count(leafD) == 0) {
            leafD->addCell(&(c.second));
            alreadyAdded.insert(leafD);
          }
        }
      }

      tree_map.insert({layer_id, tree});
    }
  }

  int thickness = 0;

  for (unsigned iev=1; iev <= nevents; iev++) {
    std::cout << "================ Simulating event: " << iev << " ================" << std::endl;

    // initialize event
    Event event(0, iev); // default run number =0

    energygen = 0.;
    energygenincells = 0.;
    energyrec = 0.;

    // Particles generation
    unsigned npart = parameters_.general().part_type.size();
    event.setnPart(npart);

    for (auto& ip : parameters_.general().part_type) {

      // particle incidents parameters
      double eta_incident, phi_incident;
      double energy_incident;

      if (parameters_.generation().eta_fluctuation) {
        eta_incident = gun_.Uniform(parameters_.generation().eta_range_min, parameters_.generation().eta_range_max);
      }
      else {
        eta_incident = parameters_.generation().incident_eta;
      }

      if (parameters_.generation().phi_fluctuation) {
        phi_incident = gun_.Uniform(parameters_.generation().phi_range_min, parameters_.generation().phi_range_max);
      }
      else {
        phi_incident = parameters_.generation().incident_phi;
      }

      if (parameters_.generation().fluctuation_energy) {
        double rand_energy = gun_.Uniform(parameters_.generation().E_range_min, parameters_.generation().E_range_max);
        if (parameters_.generation().gun_type == "E") {
          energy_incident = rand_energy;
        }
        else{
          energy_incident = rand_energy*((exp(eta_incident*eta_incident)+1)/(exp(eta_incident*eta_incident)-1));
        }
      }
      else {
        if (parameters_.generation().gun_type == "E") {
          energy_incident = parameters_.generation().energy;
        }
        else {
          energy_incident = parameters_.generation().energy*((exp(eta_incident*eta_incident)+1)/(exp(eta_incident*eta_incident)-1));
        }
      }

      int position = &ip-&parameters_.general().part_type[0];
      event.fillGenEn(position, energy_incident);
      event.fillGenEta(position, eta_incident);
      event.fillGenPhi(position, phi_incident);
      event.fillPDGid(position, ip);

      // Number of hits
      // as the hits are not generated yet, the energy is computed with the smaller energy
      // resolution in order to have the bigger number of hits per events
      // energy spot
      // no fluctuations: fixed energy = 1. / nhitspergev
      // fluctuations: alpha/sqrt(E) -> Poissonian nbr hits of energy 1/alpha^2
      // where alpha is the stochastic term of the resolution

      for (int layer_id = layer_min; layer_id < layer_max; layer_id++) {

        // incident direction
        double thetainc = 2.*std::atan(std::exp(-eta_incident));
        double r = parameters_.geometry().layers_z[layer_id]*tan(thetainc); 
        double incident_x = r*cos(phi_incident);
        double incident_y = r*sin(phi_incident);

        double r0_electro = aShowerParametrization.r0_electro_(layer_id, ip);
        double r0_hadro = aShowerParametrization.r0_hadro_(layer_id, ip);

        // take longitudinal profile as mean energy per layer for fixed energy
        double layer_weight = aShowerParametrization.getLayerProfile(ip)[layer_id];

        double denrj;
        int nhits;

        if (!parameters_.generation().fluctuation){
            denrj = 1./parameters_.generation().number_of_hits_per_gev;
            nhits = int(energy_incident*layer_weight*parameters_.generation().number_of_hits_per_gev);
        }
        else {
            denrj = aShowerParametrization.spotEnergy(ip)[2];
            nhits = gun_.Poisson(energy_incident*layer_weight/denrj);
        }
        tot_hits += nhits;

        double real_energy = 0;

        // loop over hits
        for (int i = 0; i < nhits; i++) {
          // Compute the position of each hit

          double r_shower = 0;
          if (ip == 11 || ip == 22) {
            r_shower = gun_.Exp(r0_electro); // exponential exp(-r/r0)
          }
          else { // for hadronic shower templates parametrization
            double a_core = 108.;
            double b_core = 19.3;
            double a_halo = 95.;
            double b_halo = 76.;
            r_shower = a_core * gun_.Exp(r0_hadro/b_core) + a_halo * gun_.Exp(r0_hadro/b_halo);
          }

          double phi_shower = gun_.Rndm()*TMath::TwoPi();
          double x = r_shower*cos(phi_shower) + incident_x;
          double y = r_shower*sin(phi_shower) + incident_y;

          TVectorD pos(2);
          pos(0)=x;
          pos(1)=y;

          double side = 0;
          const Cell* closestCells = nullptr;

          if (parameters_.geometry().type!=Parameters::Geometry::Type::External) {

            // hit position on layer (EE, FH or BH). Energy attribution
            int hit_pos = sqrt(x*x+y*y);
            if (hit_pos <= parameters_.geometry().limit_first_zone) {
              real_energy = aShowerParametrization.spotEnergy(ip)[0];
              side = parameters_.geometry().small_cell_side;
              thickness = 100;
            }
            else if (hit_pos >= parameters_.geometry().limit_first_zone &&
                     hit_pos <= parameters_.geometry().limit_second_zone){
              real_energy = aShowerParametrization.spotEnergy(ip)[1];
              side = parameters_.geometry().large_cell_side;
              thickness = 200;
            }
            else {
              real_energy = denrj;
              side = parameters_.geometry().large_cell_side;
              thickness = 300;
            }
            energygen += real_energy;

            double rmin = numeric_limits<double>::max();

            std::vector<const Cell*>* leafCells;

            try {
              leafCells = tree_map[layer_id]->getLeaf(float(x), float(y))->getCells();
            } catch(std::string errorString) {
              std::cout << "WARNING: " << errorString << std::endl;
              std::cout << "WARNING: not processing this hit." << std::endl;
              continue;
            }


            for (const auto& leafCell : *leafCells) {
              double leafCell_x = leafCell->getX();
              double leafCell_y = leafCell->getY();
              double r = sqrt((leafCell_x - x) * (leafCell_x - x)
                          + (leafCell_y - y) * (leafCell_y - y));

              if (r < rmin) {
                rmin = r;
                closestCells = leafCell;
              }
            }
          }
          else {
            energygen += denrj;
            closestCells = geometry_.closestCell(x,y);
          }
          // for half-cell or boarder cells, check it is within the cell
          // add energy to corresponding cell
          // Note : isincell search the position of the hit in the closest cell. When it finds the
          // cell it stop : remove the limit line problem because the hit is attributed to the 
          // first cell found

          if (closestCells == nullptr) {
            std::cout << "closestCells is null; not processing this hit" << std::endl;
            continue;
          }

          const Cell& cell = *closestCells;

          bool isincell = geometry_.isInCell(pos, cell);

          double enoise = 0;
          if (!isincell) {
            if (parameters_.general().debug) {
              std::cout << "[main] point is not inside the closest cell (hit in boarder cells or in hole at the limit)  x,y " << x << " " << y <<
                      " cell position " << cell.getX() << " " << cell.getY() <<
                      " closest cell indices " << cell.getIIndex() << " " <<  cell.getJIndex()<<
                      " layer Id : "<<layer_id + 1<< std::endl;
            }
            hit_outside_geom++;
          }
          else {
            event.fillHitCells(closestCells->getId(), cell);

            energyrec += real_energy;
            energygenincells += real_energy;
            event.fillHit(closestCells->getId(), real_energy);

            // loop on cells for the noise generation from the cells calibration
            if (parameters_.generation().noise) {
              enoise = gun_.Gaus(0., calibratednoise[layer_id][thickness/100 - 1]);
              event.fillHit(closestCells->getId(), enoise);
              energyrec += enoise;
            }

            event.fillThickness(closestCells->getId(), thickness);
          }

          // // fill shower histograms
          // double deltaR = side*sqrt(3);
          // double deltaS = 2*TMath::Pi()*r_shower*deltaR;
          // if (ip == 211) {
          //   hTransverseProfile[layer_id]->Fill(abs(r_shower), abs((real_energy + enoise)/deltaS));
          // }
          // else {
          //   hTransverseProfile[layer_id]->Fill(abs(r_shower), real_energy + enoise);
          // }

          //  hPhiProfile.Fill(phi_shower,real_energy);
          //  hSpotEnergy.Fill(real_energy);
        }
      }
    }

    std::cout << "simulated energy " << energygen << std::endl;
    std::cout << "simulated energy inside cells " << energygenincells << std::endl;
    std::cout << "reconstructed energy inside cells (includes noise) " << energyrec << std::endl;

    //Not changed since all layers implementation. 
    //FIX ME: Change the following lines in order to return for each layer or one layer the 
    //bigger energy deposit.
    // std::unique_ptr<ShowerShape> aShowerShape;
    // if (parameters_.geometry().type!=Parameters::Geometry::Type::Triangles) { // hexagons
    //   aShowerShape.reset(new ShowerShapeHexagon(event.hits(), geometry_.getCells()));
    // } else { // triangles
    //   aShowerShape.reset(new ShowerShapeTriangle(event.hits(), geometry_.getCells()));
    // }  
    // std::cout << "cell max i,j " << aShowerShape->maxCell()->getIIndex() << " " << aShowerShape->maxCell()->getJIndex()
    // << " with energy " << aShowerShape->maxE1() << std::endl;
    // std::cout << "energy in first neighboors " << aShowerShape->firstNeighboors() << std::endl;

    if (!hCellEnergyMap.empty() && !hCellEnergyEvtMap.empty()) {

      for (const auto& hit : event.hits()) {

        if (event.getLayerFromId(hit.first) == display_layer) {
          hCellEnergyMap.at(hit.first).Fill(hit.second);
        }
      }

      // if requested display a few events
      if (iev<=parameters_.display().events) {
        for (const auto& hit : event.hits()) {
          if (event.getLayerFromId(hit.first) == display_layer) {
            hCellEnergyEvtMap.at(hit.first).Reset();
            hCellEnergyEvtMap.at(hit.first).Fill(hit.second);
          }
        }

        canvas.emplace_back(display(hCellEnergyEvtMap, iev));
      }
    }

    // fill global histograms
    // hEnergySum.Fill(energyrec,1.);
    // hEnergyGen.Fill(energygenincells,1.);
    output_.fillTree(event);
  }

  // // Exporting histograms to file
  // for (int i = layer_min; i < layer_max; i++) {
  //   hTransverseProfile[i]->Write();
  // }

  if (!hCellEnergyMap.empty()) {
    canvas.emplace_back(display(hCellEnergyMap));
    // Writing energy map plots
    for(const auto& canvas_ptr : canvas) {
      canvas_ptr->Write();
    }
  }

  output_.saveTree();
  // hEnergyGen.Write();
  // hPhiProfile.Write();
  // hSpotEnergy.Write();
  // hEnergySum.Write();

  std::cout<<std::endl;
  std::cout<< "---------> Simulation information : "<<std::endl;
  std::cout << nevents << " events generated " << std::endl;
  std::cout<<hit_outside_geom<<" hits are outside or at the boarder of the geometry for a total of "
      <<tot_hits<<" hits ("<<hit_outside_geom*100./tot_hits<<"\%)."<< std::endl;

  t.Stop();
  t.Print();
  std::cout << std::endl;
}


std::unique_ptr<TCanvas> Generator::display(const std::unordered_map<uint32_t,TH1F>& hCellEnergyEvtMap, int ievt) {

  // FIXME: build titles without using char[]
  std::string title1, title2, title4;
  char str[20];
  if (ievt == 0) {
    title1 = "Mean energy profile in layer ";
  }
  else {
    title1 = "Event " + std::to_string(ievt) + " energy profile in layer ";
  }

  title1 = title1 + std::to_string(parameters_.display().layer + 1);
  title2 = "E = ";
  sprintf(str,"%4.1f",parameters_.generation().energy);
  std::string string=str;
  title2 = title2 + string;
  title2 = title2 + " GeV", 
  title4 = "eta = ";
  sprintf(str,"%3.1f",parameters_.generation().incident_eta);
  string=str;
  title4 = title4 + string;
  std::string title = title1 + ", ";
  title = title + title2;
  title = title + ", ";
  title = title + title4;

  std::unique_ptr<TCanvas> c1(new TCanvas(title.c_str(),title.c_str(),700,700));
  TH2Poly* energy_map = (TH2Poly*)geometry_.cellHistogram()->Clone(std::string("test"+std::to_string(ievt)).c_str());

  for (const auto& hist : hCellEnergyEvtMap) {
    const auto& cell = geometry_.getCells().at(hist.first);
        double enrj = hist.second.GetMean();
        energy_map->Fill(cell.getX(), cell.getY(), enrj);
        // FIXME: no sprintf
        sprintf(str,"%4.1f",hist.second.GetMean());
  }

  energy_map->Draw("colz");

  TPaveText* leg1 = new TPaveText(.05,.91,.35,.97, "NDC");
  leg1->AddText(title1.c_str());
  leg1->SetFillColor(kWhite);
  leg1->SetTextSize(0.02);
  leg1->Draw();
  TPaveText* leg2 = new TPaveText(.045,.85,.18,.88, "NDC");
  leg2->AddText(title2.c_str());
  leg2->SetFillColor(kWhite);
  leg2->SetTextSize(0.02);
  leg2->SetTextColor(kBlue);
  leg2->SetBorderSize(0.0);
  leg2->Draw();
  TPaveText* leg4 = new TPaveText(.045,.79,.25,.84, "NDC");
  leg4->AddText(title4.c_str());
  leg4->SetFillColor(kWhite);
  leg4->SetTextSize(0.02);
  leg4->SetTextColor(kBlue);
  leg4->SetBorderSize(0.0);
  leg4->Draw();

  // display incident position
  double theta = 2.*std::atan(std::exp(-parameters_.generation().incident_eta));
  double r = geometry_.getZlayer()*std::tan(theta);
  double x = r*std::cos(parameters_.generation().incident_phi);
  double y = r*std::sin(parameters_.generation().incident_phi);
  TMarker* marker = new TMarker(x,y, 24);
  marker->Draw();

  return c1;

}
