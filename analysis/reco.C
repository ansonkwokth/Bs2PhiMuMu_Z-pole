// example:
// root -q 'reco.C(1)' means reconstruction signal event


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



TParticlePDG* muonPDG = TDatabasePDG::Instance()->GetParticle(13);
Float_t mMuonPDG = muonPDG->Mass();

TParticlePDG* ZPDG = TDatabasePDG::Instance()->GetParticle(13);
Float_t mZPDG = ZPDG->Mass();
Float_t widthZPDG = ZPDG->Width();

TParticlePDG* KPDG = TDatabasePDG::Instance()->GetParticle(321);
Float_t mKPDG = KPDG->Mass();




// store features, export for python
class Features {
    // {{{
public:
    Int_t iEvt;

    // Vertex info.
    Float_t Chi2;	// Bs decay vertex fit (from the 4 tracks)
    Float_t Chi2_KK;	// Phi decay vertex fit (from the 2 Kaon tracks)
    Float_t DV_X;   	// fitted X
    Float_t DV_Y;   	// fitted Y
    Float_t DV_Z;   	// fitted Z
    Float_t DV;     	// the length from PV to fitted DV
    
    // Kin info.
    Float_t mPhi;   	// m(KK)
    Float_t mDimu;  	// m(mumu)
    Float_t mBs;    	// m(KKmumu)
    Float_t EBs;    	// E(KKmumu)
    Float_t PBs;    	// P(KKmumu)
    Float_t PKp;    	// P(K+)
    Float_t PKm;    	// P(K-)
    Float_t Pmup;   	// P(mu+)
    Float_t Pmum;   	// P(mu-)
			
    // Impact parameter
    Float_t D0Kp;   	// transverse IP of K+
    Float_t DZKp;   	// longitudinal IP of K+
    Float_t D0Km;   	// transverse IP of K-
    Float_t DZKm;   	// longitudinal IP of K-
    Float_t D0mup;   	// transverse IP of mu+
    Float_t DZmup;   	// longitudinal IP of mu+
    Float_t D0mum;   	// transverse IP of mu-
    Float_t DZmum;   	// longitudinal IP of mu-

    Float_t cosTheta_dimu;  // opening angle of the pair of muons

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



// store the indexes of the final state particles 
// Kaon+, Kaon-, Muon+, Muon-
// flag indicates if all final states are found
struct iFinalStatesIndex {
    // {{{
    Int_t iKp = -1;
    Int_t iKm = -1; 
    Int_t iMup = -1;
    Int_t iMum = -1;
    Int_t foundAll = 0;
    vector<Double_t> DV = { 99999, 99999, 99999 };    // fitted decay vertex
    vector<Double_t> PV = { 99999, 99999, 99999 };    // PV vertex (for truth level for now)
    Float_t Chi2 = 99999;
    Float_t Chi2_KK = 99999;
    Float_t mPhi = 99999;

    Int_t _havePhi = 0;    // for truth level (check if the bkg having real phi)
    Int_t _DimuRes = 0;    // for truth level (check if the Dimuon is from a resonace decay)
    Int_t _isCascade = 0;
    Int_t _BPID = -1;
    Int_t _sameB = 0;
    Int_t _pass = 0; // for truth level
    // }}}
};



// cal length of a 3D vector
Float_t calLength(Float_t X, Float_t Y, Float_t Z) {
    return pow(X*X + Y*Y + Z*Z, 0.5);
}



// truth level: check and store the signal info.
iFinalStatesIndex truthFindSig(TClonesArray* branchParticle, Int_t* BBbar) {
    // {{{
    iFinalStatesIndex iFS; 
    Int_t nParticles = branchParticle->GetEntries();
    Int_t BBbar_ = -1;  // check if Bs or bar{Bs}

    // find that the Bs the decays to phi
    Int_t foundKaonKaon = 0;
    Int_t foundMup = 0;
    Int_t foundMum = 0;
    Int_t iBs = -99999;
    Int_t nBs = 0;

    // check how many Bs
    for (Int_t ip = 0; ip < nParticles; ip++) {
        GenParticle* particle = (GenParticle*)branchParticle->At(ip);
        if (abs(particle->PID) == 531) nBs += 1;
    }
    if (nBs != 1) return iFS;

    // first search for phi that is from the Bs
    // and check its daughters are K+ and K-
    for (Int_t ip = 0; ip < nParticles; ip++) {
        // check phi is from Bs
        GenParticle* particle = (GenParticle*)branchParticle->At(ip);
        if (abs(particle->PID) != 333) continue;
        if (particle->M1 == -1) continue; GenParticle* particleM = (GenParticle*)branchParticle->At(particle->M1);
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

    // note that it may also not be always true for signal event. 
    // Because in the simulation at this stage didn't force to only store Bs samples.
    if (foundKaonKaon == 1 && foundMup == 1 && foundMum == 1) iFS._pass = 1;

    *BBbar = BBbar_;
    return iFS;
    // }}}
}


// identifing the resonance bkg: e.g. Bs0 > J/psi mu mu, Bs0 > psi(2S) mu mu, Bs0 > phi mu mu
// note that, the simulation is similar to the signal, which we didn't store the Bs0 events only
// so, we need to do the selection here
iFinalStatesIndex truthFindResBkg(TClonesArray* branchParticle, Int_t type) {
    iFinalStatesIndex iFS; 
    Int_t nParticles = branchParticle->GetEntries();

    // check if it is signal
    Int_t BBbar;
    iFinalStatesIndex iFS_sig = truthFindSig(branchParticle, &BBbar);
    if (iFS_sig._pass) return iFS;
    
    // check only one Bs
    Int_t nBs = 0;
    for (Int_t ip = 0; ip < nParticles; ip++) {
        GenParticle* particle = (GenParticle*)branchParticle->At(ip);
        if (abs(particle->PID) == 531) nBs += 1;
    }
    if (nBs == 0) return iFS;

    // check there are target final states
    // (and at least one kaon-kaon pair that has opposite charge), exactly 1 mu+ and 1 mu-
    Int_t foundKaonKaon = 0;
    Int_t foundMup = 0;
    Int_t foundMum = 0;
    for (Int_t ip = 0; ip < nParticles; ip++) {
        if (foundKaonKaon) break;
        GenParticle* particle = (GenParticle*)branchParticle->At(ip);
        if (particle->M1 == -1) continue;
        GenParticle* particleM = (GenParticle*)branchParticle->At(particle->M1);
        if (particleM->M1 == -1) continue;
        GenParticle* particleMM = (GenParticle*)branchParticle->At(particleM->M1);
        for (Int_t ip2 = 0; ip2 < nParticles; ip2++) {
            if (foundKaonKaon) break;
            GenParticle* particle2 = (GenParticle*)branchParticle->At(ip2);
            if (particle2->M1 == -1) continue;
            GenParticle* particle2M = (GenParticle*)branchParticle->At(particle2->M1);
            if (particle2M->M1 == -1) continue;
            if (particle->PID == 321 && abs(particleM->PID) == 333 && abs(particleMM->PID) == 531 &&    // check that the K+ is from phi and is from Bs
                particle2->PID == -321 && particle2->M1 == particle->M1 && particle2M->M1 == particleM->M1) {    // and the K- is from the same phi and same Bs
                foundKaonKaon = 1;    
            }
        }
    }
    if (foundKaonKaon == 0) return iFS;

    for (Int_t ip = 0; ip < nParticles; ip++) {
        GenParticle* particle = (GenParticle*)branchParticle->At(ip);
        if (particle->M1 == -1) continue;
        GenParticle* particleM = (GenParticle*)branchParticle->At(particle->M1);
        if (particleM->M1 == -1) continue;
        GenParticle* particleMM = (GenParticle*)branchParticle->At(particleM->M1);
        if (type == 3) { // Bs0 > Jpsi( > mu mu) phi
            if (particle->PID == 13 && abs(particleM->PID) == 443 && abs(particleMM->PID) == 531) foundMup = 1;
            if (particle->PID == -13 && abs(particleM->PID) == 443 && abs(particleMM->PID) == 531) foundMum = 1;
        } else if (type == 4) { // Bs0 > psi(2S) ( > Jpsi ( > mu mu ) X) K* 
            if (particleMM->M1 == -1) continue;
            GenParticle* particleMMM = (GenParticle*)branchParticle->At(particleMM->M1);
            // note that in the simulation: we let psi(2S) decay freely
            // we pick up the psi(2S) > J/psi X channels here
            if (particle->PID == 13 && abs(particleM->PID) == 443 && abs(particleMM->PID) == 100443 && abs(particleMMM->PID) == 531) foundMup = 1;
            if (particle->PID == -13 && abs(particleM->PID) == 443 && abs(particleMM->PID) == 100443 && abs(particleMMM->PID) == 531) foundMum = 1;
        } else if (type == 5) { // Bs0 > phi( > mu mu) K*
            if (particle->PID == 13 && abs(particleM->PID) == 333 && abs(particleMM->PID) == 531) foundMup = 1;
            if (particle->PID == -13 && abs(particleM->PID) == 333 && abs(particleMM->PID) == 531) foundMum = 1;
        } else {
            cout << " Wrong input type" << endl;
        }
    }

    if (foundKaonKaon == 1 && foundMup == 1 && foundMum == 1) iFS._pass = 1;
    return iFS ;
}
    

// checking background
// specify for Comb samples first, to understand the physics
iFinalStatesIndex truthFindCombBkg(TClonesArray* branchParticle, Int_t _print) {
    iFinalStatesIndex iFS; 
    Int_t nParticles = branchParticle->GetEntries();

    // check if it is signal
    Int_t BBbar;
    iFinalStatesIndex iFS_sig = truthFindSig(branchParticle, &BBbar);
    // if the event is signal, then stop and return defalut (not identified, and skip)
    /*
    if (iFS_sig._pass) return iFS;
    iFinalStatesIndex iFS_res3 = truthFindResBkg(branchParticle, 3);        
    if (iFS_res3._pass) return iFS;
    iFinalStatesIndex iFS_res4 = truthFindResBkg(branchParticle, 4);        
    if (iFS_res4._pass) return iFS;
    iFinalStatesIndex iFS_res5 = truthFindResBkg(branchParticle, 5);        
    if (iFS_res5._pass) return iFS;
    */
    
    // check there are target final states
    // (and at least one kaon-kaon pair that has opposite charge), exactly 1 mu+ and 1 mu-
    Int_t foundKaonKaon = 0;
    Int_t foundMup = 0;
    Int_t foundMum = 0;
    for (Int_t ip = 0; ip < nParticles; ip++) {
        GenParticle* particle = (GenParticle*)branchParticle->At(ip);
        if (abs(particle->PID) != 13 && particle->PID != 321) continue;
        // check there are two opposite charged muons
        if (particle->PID == 13) foundMup += 1;
        if (particle->PID == -13) foundMum += 1;
        for (Int_t ip2 = 0; ip2 < nParticles; ip2++) {
            GenParticle* particle2 = (GenParticle*)branchParticle->At(ip2);
            // check there are Kaon and Pion and with opposite charge
            if (particle->PID == 321 && particle2->PID == -321 && particle->Charge + particle2->Charge == 0) foundKaonKaon += 1; 
        }
    }

    // ============================
    // || To understand the phys ||
    // ============================
    // -- check if there is a real phi --
    // to estimate the fake-phi bkg (although should be very small)
    Int_t havePhi = 0;
    for (Int_t ip = 0; ip < nParticles; ip++) {
        GenParticle* particle = (GenParticle*)branchParticle->At(ip);
        // if there is at least one phi, then can break the loop 
        if (abs(particle->PID) == 333) { iFS._havePhi = 1; break; }
    }

    // -- check if the mu+ mu- are decayed from the resonance (like J/\psi) --
    // since previously has already force there is just 1 mu+, 1 mu- (the if statement is below)
    // so we dont need to count here, instead directly check their mother index
    Int_t muMidx = -1;
    Int_t sameM = 0;
    for (Int_t ip = 0; ip < nParticles; ip++) {
        GenParticle* particle = (GenParticle*)branchParticle->At(ip);
        // select muons
        if (not (abs(particle->PID) == 13 && particle->M1 != -1)) continue;

        // if the second muon has the same mother as the first one
        if (muMidx == particle->M1) {
            GenParticle* particleM = (GenParticle*)branchParticle->At(particle->M1);
            sameM = 1;
            break;
        }
        muMidx = particle->M1;
    }
    // make sure the 2 muons themselves form a resonance,
    // **NOT from radiation or what (e.g. ? > mu+ mu- gamma)
    // so the following check how many particles share the same mother tagged above
    Int_t nSame = 0;
    if (sameM) { 
        for (Int_t ip = 0; ip < nParticles; ip++) {
            GenParticle* particle = (GenParticle*)branchParticle->At(ip);
            if (particle->M1 == muMidx) nSame += 1;
        }
    }
    if (nSame == 2) iFS._DimuRes = 1; 

    // ------------------------------------------------------------------------------------------------------------------------------
    // // check if the 4 final state particles are from the same b
    // cout << "...." << endl;
    Int_t Midx;
    Int_t Midx_store;
    vector<Int_t> muMs = {};
    vector<Int_t> kaMs = {};
    for (Int_t ip = 0; ip < nParticles; ip++) {
        GenParticle* particle = (GenParticle*)branchParticle->At(ip);
        if (particle->M1 == -1) continue;
        if ((abs(particle->PID) != 13) && abs(particle->PID) != 321) continue;
        Midx = particle->M1;
        GenParticle* particleM;
        while (true) { 
            particleM = (GenParticle*)branchParticle->At(Midx);
            Midx_store = Midx;
            Midx = particleM->M1;
            if (abs(particleM->PID)/100 == 5 || Midx == -1) break;
        }
        if (Midx == -1 || abs(particleM->PID)/100 == 0) continue;
        if (abs(particle->PID) == 13) muMs.push_back(Midx_store);
        if (abs(particle->PID) == 321) kaMs.push_back(Midx_store);
        // cout << particle->PID << "; " << particleM->PID << "; " << Midx_store << endl;
    }

    Int_t Bidx;
    Int_t isCascade = 0;
    Int_t c1 = 0;
    for (Int_t midx1 : muMs) {
        Int_t c2 = 0;
        for (Int_t midx2 : muMs) {
            if (c1 == c2) continue;
            Int_t c3 = 0;
            for (Int_t midx3 : kaMs) {
                Int_t c4 = 0;
                for (Int_t midx4 : kaMs) {
                    if (c3 == c4) continue;
                    
                    if (midx1 == midx2 && midx2 == midx3 && midx3 == midx4) {
                        // cout <<"/////////////////////////////////////" << endl;
                        isCascade = 1;
                        Bidx = midx1;
                        break;
                    }
                    c4 += 1;
                }
                if (isCascade) break;
                c3 += 1;
            }
            if (isCascade) break;
            c2 += 1;
        }
        if (isCascade) break;
        c1 += 1;
    }
    iFS._isCascade = isCascade;
    if (isCascade) {
        GenParticle* particleB = (GenParticle*)branchParticle->At(Bidx);
        iFS._BPID = particleB->PID;
    }

    if (isCascade == 1) {

        Int_t sameMCount = 0;
        for (Int_t ip = 0; ip < nParticles; ip++) {
            GenParticle* particle = (GenParticle*)branchParticle->At(ip);
            if (particle->M1 == -1) continue;
            if (particle->Status != 1) continue;
            Midx = particle->M1;
            GenParticle* particleM;
            while (true) { 
                particleM = (GenParticle*)branchParticle->At(Midx);
                Midx_store = Midx;
                Midx = particleM->M1;
                if (Midx == Bidx) {
                    // cout << particle->PID << endl;
                    sameMCount += 1;
                    break;
                }
                if (abs(particleM->PID)/100 == 5 || Midx == -1) break;
            }
            
        }
        // cout << sameMCount << endl;
        // cout << endl;
        if (sameMCount == 4) iFS._sameB = 1;
    }
    


    

    // ------------------------------------------------------------------------------------------------------------------------------

    // found such final states in truth, it counts as bkg
    //if (foundKaonPion >= 1 && foundMup == 1 && foundMum == 1) iFS._pass = 1;
    if (foundKaonKaon > 0 && foundMup == 1 && foundMum == 1 && iFS._DimuRes == 0) iFS._pass = 1;


    vector<Int_t> muM = {};
    vector<Int_t> muMM = {};
    vector<Int_t> muMMM = {};
    for (Int_t ip = 0; ip < nParticles; ip++) {
        GenParticle* particle = (GenParticle*)branchParticle->At(ip);
        if (not (abs(particle->PID) == 13)) continue;
        if (particle->M1 != -1) continue;
        GenParticle* particleM = (GenParticle*)branchParticle->At(particle->M1);
        muM.push_back(abs(particleM->PID));
        if (particleM->M1 != -1) continue;
        GenParticle* particleMM = (GenParticle*)branchParticle->At(particleM->M1);
        muMM.push_back(abs(particleMM->PID));
        if (particleMM->M1 != -1) continue;
        GenParticle* particleMMM = (GenParticle*)branchParticle->At(particleMM->M1);
        muMMM.push_back(abs(particleMMM->PID));
    }
    vector<Int_t> kaM = {};
    vector<Int_t> kaMM = {};
    vector<Int_t> kaMMM = {};
    for (Int_t ip = 0; ip < nParticles; ip++) {
        GenParticle* particle = (GenParticle*)branchParticle->At(ip);
        if (not (abs(particle->PID) == 321)) continue;
        if (particle->M1 != -1) continue;
        GenParticle* particleM = (GenParticle*)branchParticle->At(particle->M1);
        muM.push_back(abs(particleM->PID));
        if (particleM->M1 != -1) continue;
        GenParticle* particleMM = (GenParticle*)branchParticle->At(particleM->M1);
        muMM.push_back(abs(particleMM->PID));
        if (particleMM->M1 != -1) continue;
        GenParticle* particleMMM = (GenParticle*)branchParticle->At(particleMM->M1);
        muMMM.push_back(abs(particleMMM->PID));
    }
    vector<Int_t> muArr[3] = { muM, muMM, muMMM };
    




    // print out where the muons decays from (to understand the physics)
    if (iFS._pass == 1 && iFS._DimuRes != 1 && _print) { 
        // check number of B mesons   
        Int_t nB = 0;
        for (Int_t ip = 0; ip < nParticles; ip++) {
            GenParticle* particle = (GenParticle*)branchParticle->At(ip);
            if (particle->M1 == -1) continue;
            GenParticle* particleM = (GenParticle*)branchParticle->At(particle->M1);
            if (abs(particleM->PID / 100) == 5) continue;  // exclude those with mother being B meson (avoid double counting)
            if (abs(particle->PID / 100) == 5 || abs(particle->PID / 1000) == 5) nB += 1;
        }
        cout << " nB: "  << nB << endl;

        Int_t Midx;
        Int_t nMu_ = 0; // number of muon, used to flag the dot product
        for (Int_t ip = 0; ip < nParticles; ip++) {
            GenParticle* particle = (GenParticle*)branchParticle->At(ip);
            GenParticle* particle_cp;
            Midx = particle->M1;

            // check where the muon from 
            if ((abs(particle->PID) == 13 || abs(particle->PID) == 321) && Midx != -1){
            // if ((particle->PID == 22 || abs(particle->PID) == 13 || abs(particle->PID) == 321) && Midx != -1){
                // momentum of the muon
                cout << "\n Particle: " << particle->PID << "; p=(" << particle->P4().Px() << ", " << particle->P4().Py() << ", "  << particle->P4().Pz() << ")"  << endl;
                /*
                if (nMu_ == 1) { 
                    cout << "     Dot: " << particle_cp->P4().Px() * particle->P4().Px() + 
                                            particle_cp->P4().Py() * particle->P4().Py() + 
                                            particle_cp->P4().Pz() * particle->P4().Pz() << endl;
                }
                */
                nMu_ += 1;
                particle_cp = particle; // store the 4 momentum of the first muon, for later dot product calculation
                // loop over the decay chain
                GenParticle* particleM;
                while (true) { 
                    particleM = (GenParticle*)branchParticle->At(Midx);
                    cout << " " << particleM->PID << " (" <<Midx << ")";
                    Midx = particleM->M1;
                    // stop when reach the B meson
                    if (abs(particleM->PID) == 5 || abs(particleM->PID)/100 == 5 || abs(particleM->PID)/1000 == 5 || Midx == -1) break;
                    cout << " <";
                }
                cout << endl;
            }
        }
        cout << "................................................."<< endl;
    }


    return iFS;
}




// loop over the tracks, and store the indexes of the PID wanted into a vector
// 'sign' param tell if the search is sign sensitive
// 'sign=1/-1', then only search for the specified sign
vector<Int_t> findPID(TClonesArray* branchTrack, Int_t pid, Int_t sign) {  
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
}



// find the final states
// (start to reconstruct event)
iFinalStatesIndex findFinalStatesIndex(TClonesArray* branchTrack) {
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

    for (Int_t iKp : vKp) {
        Track* kaonp = (Track*)branchTrack->At(iKp);
        for (Int_t iKm : vKm) {
            Track* kaonm = (Track*)branchTrack->At(iKm);
            // check the kaon and pion in the same direction
            if (kaonp->P4().Px() * kaonm->P4().Px() + \
                kaonp->P4().Py() * kaonm->P4().Py() + \
                kaonp->P4().Pz() * kaonm->P4().Pz() < 0) continue;
            
            TLorentzVector kaonpV, kaonmV, phiV;
            kaonpV.SetPtEtaPhiM(kaonp->PT, kaonp->Eta, kaonp->Phi, mKPDG);
            kaonmV.SetPtEtaPhiM(kaonm->PT, kaonm->Eta, kaonm->Phi, mKPDG);
            phiV = kaonpV + kaonmV;

            for (Int_t iMup : vMup) {
                Track* muonP = (Track*)branchTrack->At(iMup);
                for (Int_t iMum : vMum) {
                    Track* muonM = (Track*)branchTrack->At(iMum);
                    // check two muons are in the same direction
                    if (muonP->P4().Px() * muonM->P4().Px() + \
                        muonP->P4().Py() * muonM->P4().Py() + \
                        muonP->P4().Pz() * muonM->P4().Pz() < 0) continue;
                    // check each muon and the reconstructed phi are in the same direction
                    if (muonP->P4().Px() * phiV.Px() + \
                        muonP->P4().Py() * phiV.Py() + \
                        muonP->P4().Pz() * phiV.Pz() < 0) continue;
                    if (muonM->P4().Px() * phiV.Px() + \
                        muonM->P4().Py() * phiV.Py() + \
                        muonM->P4().Pz() * phiV.Pz() < 0) continue;
                   
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
                        iFS.mPhi = phiV.M();
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
        // inputFile = "../data/detector/ee2Z2b_comb_1M_seed0_1M_seed1_1M_seed2_2M_seed3.root";
        // oF_st = "../data/reco/ee2Z2b_comb_reco.root";
        inputFile = "../data/detector/ee2Z2b_comb_cutted_10M_seed0.root";
        oF_st = "../data/reco/ee2Z2b_comb_cutted_reco.root";
    } else if (type == 3) {
        inputFile = "../data/detector/ee2Z2Bs2PhiJpsi_1M_seed0.root";
        oF_st = "../data/reco/ee2Z2Bs2PhiJpsi_reco.root";
    } else if (type == 4) {
        inputFile = "../data/detector/ee2Z2Bs2PhiPsi_1M_seed0.root";
        oF_st = "../data/reco/ee2Z2Bs2PhiPsi_reco.root";
    } else if (type == 5) {
        inputFile = "../data/detector/ee2Z2Bs2PhiPhi_100k_seed0.root";
        oF_st = "../data/reco/ee2Z2Bs2PhiPhi_reco.root";
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


    // numberOfEntries = 10000;

    Int_t nEvt = 0; // number of true events
    Int_t nFS = 0;  // numver of events that have tagged final states

    Int_t nFS_truth = 0;

    Int_t nHavePhi = 0;
    Int_t nDimuRes = 0;
    Int_t nIsCascade = 0;
    Int_t nIsCascade_selected = 0;
    Int_t nB0Cascade = 0;
    Int_t nB0Cascade_selected = 0;
    Int_t nBpCascade = 0;
    Int_t nBpCascade_selected = 0;
    Int_t nBsCascade = 0;
    Int_t nBsCascade_selected = 0;

    Int_t nSameBCount = 0;
    Int_t nSameBCount_selected = 0;

    
    
    // loop over events
    for (Int_t i_en = 0; i_en < numberOfEntries; i_en++) {
        treeReader->ReadEntry(i_en);  // reading the entry
        if (i_en % 1000 == 0) cout << " Event: " << i_en << "/" << numberOfEntries << "(" << float(i_en) / float(numberOfEntries) * 100 << "%)" << "\r";


        /*
        if (i_en != 716 && i_en != 879 && i_en != 1387 && i_en != 1815 &&
            i_en != 2220 && i_en != 2501 && i_en != 2887 && i_en != 4654) continue;
        */
    
        cout.flush();


        // ===================
        // ||     Truth     ||
        // ===================
        Int_t BBbar;    // Bs or bar{Bs}
        iFinalStatesIndex iFS_truth ;
        if (type == 1) { 
            iFS_truth = truthFindSig(branchParticle, &BBbar);        
        } else if (type == 2) {
            //iFS_truth = truthFindCombBkg(branchParticle, 1);     // print out the strings   
            iFS_truth = truthFindCombBkg(branchParticle, 0);        
        } else if (type == 3 || type == 4 || type == 5) {
            iFS_truth = truthFindResBkg(branchParticle, type);        
        }
        else {
            iFS_truth._pass = 1;
        }
        
        if (iFS_truth._pass != 1) continue;
        nEvt += 1;

        if (iFS_truth._havePhi == 1) nHavePhi += iFS_truth._havePhi;
        if (iFS_truth._DimuRes == 1) nDimuRes += iFS_truth._DimuRes;
        if (iFS_truth._isCascade == 1) {
            nIsCascade += 1;
            if (abs(iFS_truth._BPID) == 511) nB0Cascade += 1; 
            if (abs(iFS_truth._BPID) == 521) nBpCascade += 1; 
            if (abs(iFS_truth._BPID) == 531) nBsCascade += 1; 
            if (iFS_truth._sameB == 1) nSameBCount += 1; 
        }

      
        // ===================
        // ||     Recon     ||
        // ===================
        vector<Double_t> DV;
        iFinalStatesIndex iFS = findFinalStatesIndex(branchTrack);
        if (iFS.foundAll == 0) continue;
        nFS += 1;
        if (iFS_truth._isCascade == 1) {
            nIsCascade_selected += 1;
            if (abs(iFS_truth._BPID) == 511) nB0Cascade_selected += 1; 
            if (abs(iFS_truth._BPID) == 521) nBpCascade_selected += 1; 
            if (abs(iFS_truth._BPID) == 531) nBsCascade_selected += 1; 
            if (iFS_truth._sameB == 1) nSameBCount_selected += 1; 
        }


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

        // getting the opening anlge of the pair of muons in the Bs rest frame
        TVector3 BsV3 = BsV.BoostVector();
        TLorentzVector muonpV_cp = muonpV;
        TLorentzVector muonmV_cp = muonmV;
        muonpV_cp.Boost(-BsV3);
        muonmV_cp.Boost(-BsV3);
        
        Float_t cosTheta = (muonpV_cp.Px() * muonmV_cp.Px() + muonpV_cp.Py() * muonmV_cp.Py() + muonpV_cp.Pz() * muonmV_cp.Pz())
                            / (calLength(muonpV_cp.Px(), muonpV_cp.Py(), muonpV_cp.Pz()) * calLength(muonmV_cp.Px(), muonmV_cp.Py(), muonmV_cp.Pz()));


        /*
        cout << kaonpV.Px() << ", "  << kaonpV.Py() << ", "  << kaonpV.Pz() << endl;
        cout << kaonmV.Px() << ", "  << kaonmV.Py() << ", "  << kaonmV.Pz() << endl;
        cout << muonpV.Px() << ", "  << muonpV.Py() << ", "  << muonpV.Pz() << endl;
        cout << muonmV.Px() << ", "  << muonmV.Py() << ", "  << muonmV.Pz() << endl;

        cout << "....................................................................................................... " << iFS.Chi2 << endl;
        */

        if (iFS_truth._sameB == 1) cout << " ///////////////////////////////////////////////////////////////////" << i_en << ";  " << iFS.Chi2 << endl;


        // ===================
        // ||     Store     ||
        // ===================
        features->iEvt          =   i_en;
        features->mPhi          =   iFS.mPhi;
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
    if (type == 2) {
        cout << "havePhi: "<< nHavePhi << " (" << 100*float(nHavePhi)/float(nEvt) << "% of the truth level)" << endl;
        cout << "nDimuRes:  "<< nDimuRes << " (" << 100*float(nDimuRes)/float(nEvt) << "% of the truth level)" << endl;
        cout << "niscascade: " << nIsCascade << " (" << 100*float(nIsCascade)/float(nEvt) << "% of the truth level)" << endl;
        cout << "niscascade_selected: " << nIsCascade_selected << " (" << 100*float(nIsCascade_selected)/float(nFS) << "%)" << endl; 
        cout << "Cascade from B0: " << nB0Cascade_selected << " (" << 100*float(nB0Cascade_selected)/float(nIsCascade_selected) << "%)" << endl; 
        cout << "Cascade from B+: " << nBpCascade_selected << " (" << 100*float(nBpCascade_selected)/float(nIsCascade_selected) << "%)" << endl; 
        cout << "Cascade from Bs: " << nBsCascade_selected << " (" << 100*float(nBsCascade_selected)/float(nIsCascade_selected) << "%)" << endl; 
        cout << "nSameBCount: " << nSameBCount << " (" << 100*float(nSameBCount)/float(nEvt) << "%)" << endl; 
        cout << "nSameBCount_selected: " << nSameBCount_selected << " (" << 100*float(nSameBCount_selected)/float(nFS) << "%)" << endl; 
    }

    cout << " Reco. eff.: " << nFS << "/" << nEvt << "(" << float(nFS)/float(nEvt)*100 << "%)" << endl;
    cout << endl;
    tr.Write();
    fea.Close();

}

