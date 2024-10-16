// example, do:
// root -q reco.C




#include <cmath>
#include <iostream>
#include <random>
#include <vector>
using namespace std;

#include "Rtypes.h"
#include "TChain.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TSystem.h"
#include "classes/DelphesClasses.h"
#include "ExRootAnalysis/ExRootResult.h"
#include "TrackCovariance/VertexFit.h"
#include "TrackCovariance/VertexMore.h"
#include "ExRootAnalysis/ExRootTreeReader.h"

R__ADD_LIBRARY_PATH($DELPHES)
R__LOAD_LIBRARY(lib/libDelphes)



// {{{ PDG constants
TParticlePDG* muonPDG = TDatabasePDG::Instance()->GetParticle(13);
Float_t mMuonPDG = muonPDG->Mass();
// }}} 



class Features {
    // {{{ store features, export for python
public:
    Int_t iEvt;
    
    Float_t Pmup;   	// P(mu+)
    Float_t Pmum;   	// P(mu-)
    // }}}
};





void reco() {

    cout << "\n\n\n\n\n\n\n\n\n\n\n";

    // const char* inputFile;
    // string oF_st ;
    // inputFile = "../data/detector/ee2Z2Bs2PhiMuMu_1M_seed0.root";
    // oF_st = "../data/reco/ee2Z2Bs2PhiMuMu_reco.root";

    // Load lib, and read data
    gSystem->Load("libDelphes");
    TChain chain("Delphes");
    // chain.Add(inputFile);
    // ExRootTreeReader* treeReader = new ExRootTreeReader(&chain);
    // Int_t numberOfEntries = treeReader->GetEntries();

    // clone the branches
    // TClonesArray* branchParticle = treeReader->UseBranch("Particle");

    // const char* outputFile = oF_st.c_str();
    // TFile fea(outputFile, "recreate");
    // TTree tr("t", "features");
    // Features* features = new Features;
    // tr.Branch("features", &features);



    Int_t nEvt = 0; // number of true events
    


    // loop over events
    
    for (Int_t i_en = 0; i_en < numberOfEntries; i_en++) {
        treeReader->ReadEntry(i_en);  // reading the entry
        if (i_en % 1000 == 0) cout << " Event: " << i_en << "/" << numberOfEntries << "(" << float(i_en) / float(numberOfEntries) * 100 << "%)" << "\r";
        cout.flush();

        // ===================
        // ||     Truth     ||
        // ===================
        /*
        Int_t nMu = 0;
        Int_t nParticles = branchParticle->GetEntries();
        for (Int_t ip = 0; ip < nParticles; ip++) {
            GenParticle* particle = (GenParticle*)branchParticle->At(ip);
            if (abs(particle->PID) == 13) nMu += 1;
        }
        cout << nMu << endl;
        nEvt += 1;
        */



        // ===================
        // ||     Store     ||
        // ===================
        // {{{ list of features to store
        features->iEvt          =   i_en;
        /*
        features->Pmup          =   muonpV.P();
        features->Pmum          =   muonmV.P();
        */
        tr.Fill();
        // }}}


    }
    cout << endl;
    cout << nEvt << endl;
    // cout << "Simulation samples: " << numberOfEntries << endl;
    // cout << "# of truth: " << nEvt << " (" << 100*float(nEvt)/float(numberOfEntries) << "% of the simulation)" << endl;
    // cout << " Reco. eff.: " << nFS << "/" << nEvt << "(" << float(nFS)/float(nEvt)*100 << "%)" << endl;
    // cout << endl;
    // tr.Write();
    // fea.Close();

}

