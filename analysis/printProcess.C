// example:
// root -q 'printProcess.C(1)' means printing out the truth level information 




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





Float_t calCosTheta(TVector3 aV, TVector3 bV) {
    // {{{ getting the opening angle of two 3 vectors
    Float_t cosTheta = (aV.Dot(bV)) / (aV.Mag() * bV.Mag());
    return cosTheta;
    // }}} 
}





Int_t checkDir(TClonesArray* branchParticle) {
    // {{{ check the direction of the 4 final states particle
    Int_t sameDir = 0;
    Int_t nParticles = branchParticle->GetEntries();

    // look for mu+
    for (Int_t ipmp = 0; ipmp < nParticles; ipmp++) {
        GenParticle* particleMp = (GenParticle*)branchParticle->At(ipmp);
        if (particleMp->PID != 13) continue;

        // look for mu-
        for (Int_t ipmm = 0; ipmm < nParticles; ipmm++) {
            GenParticle* particleMm = (GenParticle*)branchParticle->At(ipmm);
            if (particleMm->PID != -13) continue;
            // direction of mu+ & mu-
            if (particleMp->Px*particleMm->Px + particleMp->Py*particleMm->Py + particleMp->Pz*particleMm->Pz < 0) continue;

            // look for K+
            for (Int_t ipkp = 0; ipkp < nParticles; ipkp++) {
                GenParticle* particleKp = (GenParticle*)branchParticle->At(ipkp);
                if (particleKp->PID != 321) continue;
                // direction of mu+ & K+, mu- & K+
                if (particleMp->Px*particleKp->Px + particleMp->Py*particleKp->Py + particleMp->Pz*particleKp->Pz < 0 || 
                    particleMm->Px*particleKp->Px + particleMm->Py*particleKp->Py + particleMm->Pz*particleKp->Pz < 0) continue;
        
                // look for K-
                for (Int_t ipkm = 0; ipkm < nParticles; ipkm++) {
                    GenParticle* particleKm = (GenParticle*)branchParticle->At(ipkm);
                    if (particleKm->PID != -321) continue;
                    // direction of mu+ & K-, mu- & K-, K+ & K-
                    if (particleMp->Px*particleKm->Px + particleMp->Py*particleKm->Py + particleMp->Pz*particleKm->Pz < 0 || 
                        particleMm->Px*particleKm->Px + particleMm->Py*particleKm->Py + particleMm->Pz*particleKm->Pz < 0 || 
                        particleKp->Px*particleKm->Px + particleKp->Py*particleKm->Py + particleKp->Pz*particleKm->Pz < 0) continue;
                    sameDir = 1;
                }
            }
        }
    }
    return sameDir;
    // }}}
}

Int_t checkDir_det(TClonesArray* branchTrack) {
    // {{{ check the direction of the 4 final states particle
    Int_t sameDir = 0;
    Int_t nParticles = branchTrack->GetEntries();

    // look for mu+
    for (Int_t ipmp = 0; ipmp < nParticles; ipmp++) {
        Track* particleMp = (Track*)branchTrack->At(ipmp);
        if (particleMp->PID != 13) continue;

        // look for mu-
        for (Int_t ipmm = 0; ipmm < nParticles; ipmm++) {
            Track* particleMm = (Track*)branchTrack->At(ipmm);
            if (particleMm->PID != -13) continue;
            // direction of mu+ & mu-
            if (particleMp->P4().Px()*particleMm->P4().Px() + 
                particleMp->P4().Py()*particleMm->P4().Py() + 
                particleMp->P4().Pz()*particleMm->P4().Pz() < 0) continue;

            // look for K+
            for (Int_t ipkp = 0; ipkp < nParticles; ipkp++) {
                Track* particleKp = (Track*)branchTrack->At(ipkp);
                if (particleKp->PID != 321) continue;
                // direction of mu+ & K+, mu- & K+
                if (particleMp->P4().Px()*particleKp->P4().Px() + 
                    particleMp->P4().Py()*particleKp->P4().Py() + 
                    particleMp->P4().Pz()*particleKp->P4().Pz() < 0 || 
                    particleMm->P4().Px()*particleKp->P4().Px() + 
                    particleMm->P4().Py()*particleKp->P4().Py() + 
                    particleMm->P4().Pz()*particleKp->P4().Pz() < 0) continue;
        
                // look for K-
                for (Int_t ipkm = 0; ipkm < nParticles; ipkm++) {
                    Track* particleKm = (Track*)branchTrack->At(ipkm);
                    if (particleKm->PID != -321) continue;
                    // direction of mu+ & K-, mu- & K-, K+ & K-
                    if (particleMp->P4().Px()*particleKm->P4().Px() + 
                        particleMp->P4().Py()*particleKm->P4().Py() + 
                        particleMp->P4().Pz()*particleKm->P4().Pz() < 0 || 
                        particleMm->P4().Px()*particleKm->P4().Px() + 
                        particleMm->P4().Py()*particleKm->P4().Py() + 
                        particleMm->P4().Pz()*particleKm->P4().Pz() < 0 || 
                        particleKp->P4().Px()*particleKm->P4().Px() + 
                        particleKp->P4().Py()*particleKm->P4().Py() + 
                        particleKp->P4().Pz()*particleKm->P4().Pz() < 0) continue;
                    sameDir = 1;
                }
            }
        }
    }
    return sameDir;
    // }}}
}


Int_t checkPhi(TClonesArray* branchParticle) {
    // {{{ check if in the truth level, there is phi>KK 
    Int_t nParticles = branchParticle->GetEntries();
    Int_t havePhi = 0;
    for (Int_t ip = 0; ip < nParticles; ip++) {
        GenParticle* particle = (GenParticle*)branchParticle->At(ip);
        // if there is at least one phi, then can break the loop 
        if (abs(particle->PID) == 333) { 
            GenParticle* particleD1 = (GenParticle*)branchParticle->At(particle->D1);
            GenParticle* particleD2 = (GenParticle*)branchParticle->At(particle->D2);
            if ((particleD1->PID == 321 && particleD2->PID == -321) || (particleD1->PID == -321 && particleD2->PID == 321)) {
                havePhi = 1;
                break; 
            }
        }
    }
    return havePhi; 
    // }}}
}

Int_t checkMuRes(TClonesArray* branchParticle) {
    // {{{ check if the pair of muons are from resonance
    Int_t nParticles = branchParticle->GetEntries();

    Int_t muMidx = -1;
    Int_t sameM = 0;
    for (Int_t ip = 0; ip < nParticles; ip++) {
        GenParticle* particle = (GenParticle*)branchParticle->At(ip);
        // select muons
        if (not (abs(particle->PID) == 13 && particle->M1 != -1)) continue;
        // if the second muon has the same mother as the first one
        if (muMidx == particle->M1) {
            sameM = 1;
            break;
        }
        muMidx = particle->M1;
    }
    // make sure the 2 muons themselves form a resonance,
    // **NOT from radiation or what (e.g. ? > mu+ mu- gamma)
    // so the following check how many particles share the same mother tagged above
    Int_t nSame = 0;
    if (not sameM) return 0;

    for (Int_t ip = 0; ip < nParticles; ip++) {
        GenParticle* particle = (GenParticle*)branchParticle->At(ip);
        if (particle->M1 == muMidx) nSame += 1;
    }
    if (nSame == 2) return 1; 
    else return 0;
    // }}}
}

void printChain(TClonesArray* branchParticle) {
    //{{{ print out the decay chain
    Int_t nParticles = branchParticle->GetEntries();
    Int_t Midx;
    cout << endl;
    for (Int_t ip = 0; ip < nParticles; ip++) {
        GenParticle* particle = (GenParticle*)branchParticle->At(ip);
        if (particle->M1 == -1) continue;
        if ((abs(particle->PID) != 13) && abs(particle->PID) != 321) continue;
        cout << endl;
        cout << particle->PID << " (" << ip << ")";
        Midx = particle->M1;
        GenParticle* particleM;
        while (true) { 
            particleM = (GenParticle*)branchParticle->At(Midx);
            cout << " < " << particleM->PID << " (" << Midx << ")";
            Midx = particleM->M1;
            if (Midx == -1 || abs(particleM->PID) == 23) break;
        }
        if (Midx == -1 || abs(particleM->PID)/100 == 0) continue;
    }
    //}}}
}



struct flags {
    Int_t sameDir = -1;
    Int_t havePhi = -1;
    Int_t fromRes = -1;
};

flags truthLevelPhys(TClonesArray* branchParticle, Int_t print) {
    flags F;
    
    F.sameDir  = checkDir(branchParticle);
    F.havePhi  = checkPhi(branchParticle);
    F.fromRes  = checkMuRes(branchParticle);
    if (print) printChain(branchParticle);

    return F;
}






void printProcess(Int_t type, Int_t print) {

    cout << "\n\n\n\n\n\n\n\n\n\n\n";

    string typeName;
    const char* inputFile;
    string oF_st ;
    if (type == 1) {
        inputFile = "../data/detector/ee2Z2Bs2PhiMuMu_1M_seed0.root";
        oF_st = "../data/reco/ee2Z2Bs2PhiMuMu_reco.root";
    } else if (type == 2) {
        inputFile = "../data/detector/ee2Z2b_comb_cutted_500M_seed0-seed22.root";
        oF_st = "../data/reco/ee2Z2b_comb_cutted_reco.root";
    } else if (type == 3) {
        inputFile = "../data/detector/ee2Z2c_comb_cutted_100M_seed0-seed4.root";
        oF_st = "../data/reco/ee2Z2c_comb_cutted_reco.root";
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



    Int_t nEvt = 0; // number of true events
    Int_t nFS = 0;  // numver of events that have tagged final states
    Int_t nFS_truth = 0;

    Int_t nSameDir = 0;
    Int_t nSameDir_det = 0;
    Int_t nHavePhi = 0;
    Int_t nFromRes = 0;


    numberOfEntries = 10000;
    // loop over events
    for (Int_t i_en = 0; i_en < numberOfEntries; i_en++) {
        treeReader->ReadEntry(i_en);  // reading the entry
        if (i_en % 1000 == 0) cout << " Event: " << i_en << "/" << numberOfEntries << "(" << float(i_en) / float(numberOfEntries) * 100 << "%)" << "\r";
        cout.flush();


        if (print) cout << "Event: " << i_en << endl;
        // ===================
        // ||     Truth     ||
        // ===================
        flags F = truthLevelPhys(branchParticle, print);
        if (F.sameDir == 1) nSameDir += 1;
        if (checkDir_det(branchTrack) == 1) nSameDir_det += 1;
        if (F.havePhi == 1) nHavePhi += 1;
        if (F.fromRes == 1) nFromRes += 1;
    }

    cout << endl;
    cout << "\tnSameDir:      " << nSameDir << "/" << numberOfEntries << "\t(" << float(nSameDir)/float(numberOfEntries) * 100 << "%)"  << endl;
    cout << "\tnSameDir_det:  " << nSameDir_det << "/" << numberOfEntries << "\t(" << float(nSameDir_det)/float(numberOfEntries) * 100 << "%)"  << endl;
    cout << "\tnHavePhi:      " << nHavePhi << "/" << numberOfEntries << "\t(" << float(nHavePhi)/float(numberOfEntries) * 100 << "%)"  << endl;
    cout << "\tnFromRes:      " << nFromRes << "/" << numberOfEntries << "\t(" << float(nFromRes)/float(numberOfEntries) * 100 << "%)"  << endl;

}

