// example:
// root -q 'reco.C(1)' means reconstruction signal event

//TODO: 
// check the _sameB part: why nu is not included
// change the "same direction requirement

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

TParticlePDG* KPDG = TDatabasePDG::Instance()->GetParticle(321);
Float_t mKPDG = KPDG->Mass();
// }}} 



class Features {
    // {{{ store features, export for python
public:
    Int_t iEvt;

    // Vertex info.
    Float_t Chi2;	    // Bs decay vertex fit (from the 4 tracks)
    Float_t Chi2_KK;	// Phi decay vertex fit (from the 2 Kaon tracks)
    Float_t DV_X;   	// fitted X
    Float_t DV_Y;   	// fitted Y
    Float_t DV_Z;   	// fitted Z
    Float_t DV;     	// the length from PV to fitted DV
    
    // reco info.
    Float_t mPhi;   	// m(KK)
    Float_t mDimu;  	// m(mumu)
    Float_t mBs;    	// m(KKmumu)
    Float_t EBs;    	// E(KKmumu)
    Float_t PBs;    	// P(KKmumu)
    Float_t PKp;    	// P(K+)
    Float_t PKm;    	// P(K-)
    Float_t Pmup;   	// P(mu+)
    Float_t Pmum;   	// P(mu-)
    Float_t cosTheta_dimu;  // opening angle of the pair of muons
			
    // Impact parameters
    Float_t D0Kp;   	// transverse IP of K+
    Float_t DZKp;   	// longitudinal IP of K+
    Float_t D0Km;   	// transverse IP of K-
    Float_t DZKm;   	// longitudinal IP of K-
    Float_t D0mup;   	// transverse IP of mu+
    Float_t DZmup;   	// longitudinal IP of mu+
    Float_t D0mum;   	// transverse IP of mu-
    Float_t DZmum;   	// longitudinal IP of mu-

    // truth level
    Float_t BBbar_truth;    // check if the Bs or bar{Bs}
    Float_t DV_X_truth;
    Float_t DV_Y_truth;
    Float_t DV_Z_truth;
    Float_t DV_truth;
    Float_t PV_X_truth;
    Float_t PV_Y_truth;
    Float_t PV_Z_truth;

    // }}}
};



struct iFinalStatesIndex {
    // {{{ store the indexes of the final state particles 

    // indeces of the final states
    Int_t iKp = -1;     // K+
    Int_t iKm = -1;     // K-
    Int_t iMup = -1;    // mu+
    Int_t iMum = -1;    // mu-
    // flag indicates if all final states are found
    Int_t foundAll = 0;
    // vertex info (KKmumu vertex)
    vector<Double_t> DV = { 99999, 99999, 99999 };    // fitted decay vertex
    vector<Double_t> PV = { 99999, 99999, 99999 };    // PV vertex (for truth level for now)
    Float_t Chi2 = 99999;       // fitted chi2 of the DV from 4 tracks 
    Float_t Chi2_KK = 99999;    // fitted chi2 of the vertex from only K+K- tracks

    Int_t _pass = 0;        // for truth level

    // }}}
};



Float_t calLength(Float_t X, Float_t Y, Float_t Z) { 
    // {{{ cal length of a 3D vector
    return pow(X*X + Y*Y + Z*Z, 0.5); 
    // }}}
}



Float_t calCosTheta(TLorentzVector BsV, TLorentzVector muonpV, TLorentzVector muonmV) {
    // {{{ getting the opening anlge of the pair of muons in the Bs rest frame
    TVector3 BsV3 = BsV.BoostVector();
    TLorentzVector muonpV_cp = muonpV;
    TLorentzVector muonmV_cp = muonmV;
    muonpV_cp.Boost(-BsV3);
    muonmV_cp.Boost(-BsV3);
    
    Float_t cosTheta = (muonpV_cp.Px() * muonmV_cp.Px() + muonpV_cp.Py() * muonmV_cp.Py() + muonpV_cp.Pz() * muonmV_cp.Pz())
                        / (calLength(muonpV_cp.Px(), muonpV_cp.Py(), muonpV_cp.Pz()) * calLength(muonmV_cp.Px(), muonmV_cp.Py(), muonmV_cp.Pz()));
    return cosTheta;
    // }}}
}


iFinalStatesIndex truthFindSig(TClonesArray* branchParticle, Int_t* BBbar) {
    // {{{ truth level: check and store the signal info.
    iFinalStatesIndex iFS; 
    Int_t nParticles = branchParticle->GetEntries();
    Int_t BBbar_ = -1;  // check if Bs or bar{Bs}

    // find the final state particles
    Int_t foundKaonKaon = 0;
    Int_t foundMup = 0;
    Int_t foundMum = 0;

    Int_t iBs = -99999;
    Int_t nBs = 0;  // number of Bs

    // check how many Bs
    for (Int_t ip = 0; ip < nParticles; ip++) {
        GenParticle* particle = (GenParticle*)branchParticle->At(ip);
        if (abs(particle->PID) == 531) nBs += 1;
    }
    // if (nBs != 1) return iFS;

    // first search for phi that is from the Bs
    // and check its daughters are K+ and K-
    for (Int_t ip = 0; ip < nParticles; ip++) {
        // check phi is from Bs
        GenParticle* particle = (GenParticle*)branchParticle->At(ip);
        if (abs(particle->PID) != 333) continue;
        if (particle->M1 == -1) continue; 
        GenParticle* particleM = (GenParticle*)branchParticle->At(particle->M1);
        if (abs(particleM->PID) != 531) continue;
        
        // check KK from phi
        GenParticle* particleD1 = (GenParticle*)branchParticle->At(particle->D1);
        GenParticle* particleD2 = (GenParticle*)branchParticle->At(particle->D2);
        if (not (abs(particleD1->PID) == 321 && abs(particleD2->PID) == 321)) continue;
        foundKaonKaon = 1;

        iBs = particle->M1;
        BBbar_ = Int_t(particleM->PID / abs(particleM->PID));   // sign of the Bs 
        iFS.DV = { particle->X, particle->Y, particle->Z };     // DV: decay vertex of phi
        iFS.PV = { particleM->X, particleM->Y, particleM->Z };  // PV: decay vertex of Bs
    }
    if (foundKaonKaon != 1) return iFS;
    
    // check mu+ mu- from the same Bs
    for (Int_t ip = 0; ip < nParticles; ip++) {
        GenParticle* particle = (GenParticle*)branchParticle->At(ip);
        if (particle->M1 == -1) continue;
        if (particle->PID == 13 && particle->M1 == iBs) foundMup = 1;
        if (particle->PID == -13 && particle->M1 == iBs) foundMum = 1;
    }

    // note that it may also not be always true for signal event (there are nBs=2 cases). 
    if (foundKaonKaon == 1 && foundMup == 1 && foundMum == 1) iFS._pass = 1;

    *BBbar = BBbar_;
    return iFS;
    // }}}
}



vector<Int_t> findPID(TClonesArray* branchTrack, Int_t pid, Int_t sign) {  
    // {{{ store te indeces of the PID wanted into a vector
    // loop over the tracks, and store the indexes of the PID wanted into a vector
    // 'sign' param tell if the search is sign sensitive
    // 'sign=1/-1', then only search for the specified sign
    Int_t nTracks = branchTrack->GetEntries();
    vector<Int_t> v = {};
    for (Int_t itr1 = 0; itr1 < nTracks; itr1++) {
        Track* track1 = (Track*)branchTrack->At(itr1);
        if (sign) { 
            if (track1->PID != pid) continue;
        } else {
            if (abs(track1->PID) != pid) continue;
        }
        v.push_back(itr1);
    }
    return v;
    // }}}
}



iFinalStatesIndex findFinalStatesIndex(TClonesArray* branchTrack) {
    // {{{ find the final states
    // find the final states
    // (start to reconstruct event)
    iFinalStatesIndex iFS;

    // store all kaon PID
    vector<Int_t> vKp = findPID(branchTrack, 321, 1);
    vector<Int_t> vKm = findPID(branchTrack, -321, 1);
    // store all muon PID
    vector<Int_t> vMup = findPID(branchTrack, -13, 1);
    vector<Int_t> vMum = findPID(branchTrack, 13, 1);

    // store the best combination giving min chi2
    Float_t minChi2 = 99999;
    vector<Double_t> DV_;

    // loop K+
    for (Int_t iKp : vKp) {
        Track* kaonp = (Track*)branchTrack->At(iKp);
        // loop K-
        for (Int_t iKm : vKm) {
            Track* kaonm = (Track*)branchTrack->At(iKm);
            // check the K+ and K- in the same direction
            if (kaonp->P4().Px() * kaonm->P4().Px() + \
                kaonp->P4().Py() * kaonm->P4().Py() + \
                kaonp->P4().Pz() * kaonm->P4().Pz() < 0) continue;
            
            // loop mu+
            for (Int_t iMup : vMup) {
                Track* muonp = (Track*)branchTrack->At(iMup);
                // check the K+ and mu+ in the same direction
                if (kaonp->P4().Px() * muonp->P4().Px() + \
                    kaonp->P4().Py() * muonp->P4().Py() + \
                    kaonp->P4().Pz() * muonp->P4().Pz() < 0) continue;
                // check the K- and mu+ in the same direction
                if (kaonm->P4().Px() * muonp->P4().Px() + \
                    kaonm->P4().Py() * muonp->P4().Py() + \
                    kaonm->P4().Pz() * muonp->P4().Pz() < 0) continue;
            
                // loop mu-
                for (Int_t iMum : vMum) {
                    Track* muonm = (Track*)branchTrack->At(iMum);
                    // check the K+ and mu- in the same direction
                    if (kaonp->P4().Px() * muonm->P4().Px() + \
                        kaonp->P4().Py() * muonm->P4().Py() + \
                        kaonp->P4().Pz() * muonm->P4().Pz() < 0) continue;
                    // check the K- and mu- in the same direction
                    if (kaonm->P4().Px() * muonm->P4().Px() + \
                        kaonm->P4().Py() * muonm->P4().Py() + \
                        kaonm->P4().Pz() * muonm->P4().Pz() < 0) continue;
                    // check the mu+ and mu- in the same direction
                    if (muonp->P4().Px() * muonm->P4().Px() + \
                        muonp->P4().Py() * muonm->P4().Py() + \
                        muonp->P4().Pz() * muonm->P4().Pz() < 0) continue;



                    // vertex fitting (all 4 final states tracks)
                    TVectorD* pr[4];
                    TMatrixDSym* cv[4];
                    vector<Int_t> idxes{ iKp, iKm, iMup, iMum };
                    Int_t ii = 0;
                    for (Int_t idx : idxes) {
                        Track* track = (Track*)branchTrack->At(idx);
                        pr[ii] = new TVectorD(0, 4, track->D0, track->Phi, track->C, track->DZ, track->CtgTheta, "END");
                        cv[ii] = new TMatrixDSym(track->CovarianceMatrix());
                        ii += 1;
                    }
                    VertexFit* Vtx = new VertexFit(4, pr, cv);
                    TVectorD xvtx_ = Vtx->GetVtx();
                    TMatrixDSym covX_ = Vtx->GetVtxCov();
                    Double_t Chi2_ = Vtx->GetVtxChi2();

                    // vertex fitting (2 Kaon tracks)
                    TVectorD* pr_KK[2];
                    TMatrixDSym* cv_KK[2];
                    vector<Int_t> idxes_KK{ iKp, iKm };
                    Int_t ii_KK = 0;
                    for (Int_t idx_KK : idxes_KK) {
                        Track* track = (Track*)branchTrack->At(idx_KK);
                        pr_KK[ii_KK] = new TVectorD(0, 4, track->D0, track->Phi, track->C, track->DZ, track->CtgTheta, "END");
                        cv_KK[ii_KK] = new TMatrixDSym(track->CovarianceMatrix());
                        ii_KK += 1;
                    }
                    VertexFit* Vtx_KK = new VertexFit(2, pr_KK, cv_KK);
                    Double_t Chi2_KK_ = Vtx_KK->GetVtxChi2();

                    // check if this is a better combination
                    if (iFS.Chi2 > Chi2_) {
                        iFS.Chi2 = Chi2_;
                        iFS.Chi2_KK = Chi2_KK_;
                        iFS.DV = { xvtx_[0], xvtx_[1], xvtx_[2] };
                        iFS.iKp = iKp;
                        iFS.iKm = iKm;
                        iFS.iMup = iMup;
                        iFS.iMum = iMum;
                        iFS.foundAll = 1;
                    }
                }
            }
        }
    }
    return  iFS;
    // }}}
}




void reco(Int_t type) {

    cout << "\n\n\n\n\n\n\n\n\n\n\n";

    string typeName;
    const char* inputFile;
    string oF_st ;
    if (type == 1) {
        inputFile = "../data/detector/ee2Z2Bs2PhiMuMu_1M_seed0.root";
        oF_st = "../data/reco/ee2Z2Bs2PhiMuMu_reco.root";
    } else if (type == 2) {
        inputFile = "../data/detector/ee2Z2b_comb_cutted_10M_seed0-10M_seed1-10M_seed2-10M_seed3-10M_seed4-25M_seed5-25M_seed6-25M_seed7-25M_seed8.root";
        oF_st = "../data/reco/ee2Z2b_comb_cutted_reco.root";
    }



    // Load lib, and read data
    gSystem->Load("libDelphes");
    TChain chain("Delphes");
    //TChain chain("events");
    chain.Add(inputFile);
    ExRootTreeReader* treeReader = new ExRootTreeReader(&chain);
    Int_t numberOfEntries = treeReader->GetEntries();

    // clone the branches
    TClonesArray* branchParticle = treeReader->UseBranch("Particle");
    TClonesArray* branchTrack = treeReader->UseBranch("Track");

    const char* outputFile = oF_st.c_str();
    TFile fea(outputFile, "recreate");
    TTree tr("t", "features");
    Features* features = new Features;
    tr.Branch("features", &features);



    Int_t nEvt = 0; // number of true events
    Int_t nFS = 0;  // numver of events that have tagged final states
    Int_t nFS_truth = 0;
    


    // numberOfEntries = 100;
    // loop over events
    for (Int_t i_en = 0; i_en < numberOfEntries; i_en++) {
        treeReader->ReadEntry(i_en);  // reading the entry
        if (i_en % 1000 == 0) cout << " Event: " << i_en << "/" << numberOfEntries << "(" << float(i_en) / float(numberOfEntries) * 100 << "%)" << "\r";
        cout.flush();

        // ===================
        // ||     Truth     ||
        // ===================
        Int_t BBbar;    // Bs or bar{Bs}
        iFinalStatesIndex iFS_truth ;
        iFS_truth = truthFindSig(branchParticle, &BBbar);        
        // check in truth level and see if the event is in the category that we want
        if (type != 1 && type != 2) {
            cout << "Wrong Type input" << endl;
            continue;
        } else if (type == 1 && iFS_truth._pass != 1) { // targeting signal, and we found it in truth level
            continue;   
        } else if (type == 2 && iFS_truth._pass == 1) { // targeting Z>bb bkg., and we found it is not signal event
            continue;
        }
        nEvt += 1;


      
        // ===================
        // ||     Recon     ||
        // ===================
        vector<Double_t> DV;
        iFinalStatesIndex iFS = findFinalStatesIndex(branchTrack);
        if (iFS.foundAll == 0) continue;
        nFS += 1;



        // define final state lorentz vector 
        TLorentzVector kaonpV, kaonmV, phiV, muonpV, muonmV, BsV, dimuV;
        Track* kaonp = (Track*)branchTrack->At(iFS.iKp);
        kaonpV.SetPtEtaPhiM(kaonp->PT, kaonp->Eta, kaonp->Phi, mKPDG);
        Track* kaonm = (Track*)branchTrack->At(iFS.iKm);
        kaonmV.SetPtEtaPhiM(kaonm->PT, kaonm->Eta, kaonm->Phi, mKPDG);
        Track* muonp = (Track*)branchTrack->At(iFS.iMup);
        muonpV.SetPtEtaPhiM(muonp->PT, muonp->Eta, muonp->Phi, mMuonPDG);
        Track* muonm = (Track*)branchTrack->At(iFS.iMum);
        muonmV.SetPtEtaPhiM(muonm->PT, muonm->Eta, muonm->Phi, mMuonPDG);
        // reconstruct K*, B0
        phiV = kaonpV + kaonmV;
        BsV = phiV + muonpV + muonmV;
        dimuV = muonpV + muonmV;

        Float_t cosTheta = calCosTheta(BsV, muonpV, muonmV);  


        // ===================
        // ||     Store     ||
        // ===================
        features->iEvt          =   i_en;
        features->mPhi          =   phiV.M();
        features->mBs           =   BsV.M();
        features->EBs           =   BsV.E();
        features->PBs           =   BsV.P();
        features->DV_X          =   iFS.DV[0];
        features->DV_Y          =   iFS.DV[1];
        features->DV_Z          =   iFS.DV[2];
        features->DV            =   calLength(iFS.DV[0], iFS.DV[1], iFS.DV[2]); 
        features->Chi2          =   iFS.Chi2;
        features->Chi2_KK       =   iFS.Chi2_KK;
        features->BBbar_truth   =   BBbar;
        features->mDimu         =   dimuV.M();
        features->PKp           =   kaonpV.P();
        features->PKm           =   kaonmV.P();
        features->Pmup          =   muonpV.P();
        features->Pmum          =   muonmV.P();
        features->D0Kp          =   kaonp->D0;
        features->DZKp          =   kaonp->DZ;
        features->D0Km          =   kaonm->D0;
        features->DZKm          =   kaonm->DZ;
        features->D0mup         =   muonp->D0;
        features->DZmup         =   muonp->DZ;
        features->D0mum         =   muonm->D0;
        features->DZmum         =   muonm->DZ;
        features->cosTheta_dimu =   cosTheta;
        features->DV_X_truth    =   iFS_truth.DV[0];
        features->DV_Y_truth    =   iFS_truth.DV[1];
        features->DV_Z_truth    =   iFS_truth.DV[2];
        features->DV_truth      =   calLength(iFS_truth.DV[0] - iFS_truth.PV[0], iFS_truth.DV[1] - iFS_truth.PV[1], iFS_truth.DV[2] - iFS_truth.PV[2]); 
        tr.Fill();


    }
    cout << endl;
    cout << "Simulation samples: " << numberOfEntries << endl;
    cout << "# of truth: " << nEvt << " (" << 100*float(nEvt)/float(numberOfEntries) << "% of the simulation)" << endl;
    cout << " Reco. eff.: " << nFS << "/" << nEvt << "(" << float(nFS)/float(nEvt)*100 << "%)" << endl;
    cout << endl;
    tr.Write();
    fea.Close();

}

