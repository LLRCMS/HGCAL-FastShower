#ifdef STANDALONE
#include "OutputService.h"
#else
#include "HGCalSimulation/FastShower/interface/OutputService.h"
#endif

#include <iostream>
using namespace std;


OutputService::OutputService(const std::string& file_name):
    file_(TFile::Open(file_name.c_str(), "recreate")),
    tree_(new TTree("tree", "tree")),
    run_(0),
    event_(0),
    npart_(0),
    PDGid_(0),
    thick_(0),
    cell_n_() {
      tree_->Branch("run", &run_, "run/i");
      tree_->Branch("event", &event_, "event/i");
      tree_->Branch("npart", &npart_, "npart/i");

      tree_->Branch("PDGid", &PDGid_);
      tree_->Branch("gen_energy", &gen_energy_);
      tree_->Branch("gen_eta", &gen_eta_);
      tree_->Branch("gen_phi", &gen_phi_);

      tree_->Branch("layer", &layer_);
      tree_->Branch("thickness", &thick_);
      tree_->Branch("cell_n", &cell_n_, "cell_n/i");

      tree_->Branch("energy", &cell_energy_);
      tree_->Branch("x", &cell_x_);
      tree_->Branch("y", &cell_y_);
      tree_->Branch("z", &cell_z_);
      tree_->Branch("eta", &cell_eta_);
      tree_->Branch("phi", &cell_phi_);
}

OutputService::~OutputService() {
}


void OutputService::fillTree(const Event& event) {
    clear();
    run_ = event.run();
    event_ = event.event();
    npart_ = event.npart();

    for(auto& id_hit : event.hits()) {
        // skip zero and negative energies
        if(id_hit.second<=0.)
            continue;
        cell_energy_.emplace_back(id_hit.second);
    }
    cell_n_ = cell_energy_.size();

    int n = 0;
    for(auto& c : event.cells()) {
        double x = c.second.getX();
        double y = c.second.getY();
        double z = c.second.getZ();

        double r = std::sqrt(x*x + y*y);
        double theta = std::atan(r/z);
        double eta = -std::log(std::tan(theta/2.));
        double phi = std::copysign(std::acos(x/r),y);

        cell_x_.emplace_back(x);
        cell_y_.emplace_back(y);
        cell_z_.emplace_back(z);
        cell_eta_.emplace_back(eta);
        cell_phi_.emplace_back(phi);

        double layer = c.second.getLayer();
        cout << "hit number : "<<n++<<endl;
        cout << "cellid"<<c.first<<"  layer : "<<layer<<endl;
        layer_.emplace_back(layer);
    }
    int c_l = 0;
    for (auto& l : layer_){
        if (l == 2)
            c_l++;
        }
    cout <<"element in layer 2 : "<<c_l<<endl;


    for(const auto& id_part : event.gen_en()) {
        if(id_part.second<=0.)
            continue;
        gen_energy_.emplace_back(id_part.second);
    }
    for(const auto& id_part : event.gen_eta())
        gen_eta_.emplace_back(id_part.second);
    for(const auto& id_part : event.gen_phi())
        gen_phi_.emplace_back(id_part.second);
    for(const auto& id_part : event.pdg_id())
        PDGid_.emplace_back(id_part.second);
    for(const auto& id_part : event.thick()) {
        thick_.emplace_back(id_part.second);
    }

    tree_->Fill();
}

void OutputService::saveTree(){
    tree_->Write();
    file_->Close();
}

void OutputService::clear() {
    run_ = 0;
    event_ = 0;
    npart_ = 0.;
    PDGid_.clear();
    thick_.clear();
    gen_energy_.clear();
    gen_eta_.clear();
    gen_phi_.clear();
    layer_.clear();
    cell_n_ = 0;
    cell_energy_.clear();
    cell_x_.clear();
    cell_y_.clear();
    cell_z_.clear();
    cell_eta_.clear();
    cell_phi_.clear();
}



