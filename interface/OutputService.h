#ifndef __HGCalSimulation_FastShower_OutputService_h__
#define __HGCalSimulation_FastShower_OutputService_h__


#include <vector>

#include "TFile.h"
#include "TTree.h"

#ifdef STANDALONE
#include "Event.h"
#include "Geometry.h"
#else
#include "HGCalSimulation/FastShower/interface/Event.h"
#include "HGCalSimulation/FastShower/interface/Geometry.h"
#endif


class OutputService {
    public:
        OutputService(const std::string&);
        ~OutputService();

        void fillTree(const Event&);
        void saveTree();
        void clear();


    private:
        std::unique_ptr<TFile> file_;
        TTree* tree_; // the tree is owned by file_ and deleted when the file is closed

        // tree branches
        unsigned run_;
        unsigned event_;
        unsigned npart_;
        unsigned cell_n_;

        std::vector<float> PDGid_;
        std::vector<float> thick_;
        std::vector<float> gen_energy_;
        std::vector<float> gen_eta_;
        std::vector<float> gen_phi_;
        std::vector<float> cell_energy_;
        std::vector<float> cell_x_;
        std::vector<float> cell_y_;
        std::vector<float> cell_z_;
        std::vector<float> cell_eta_;
        std::vector<float> cell_phi_;

};

#endif
