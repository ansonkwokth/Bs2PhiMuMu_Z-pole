// Simulating e+ e- > Z, Z > c \bar{c} background, at Z-pole
// Requring that there must be at least 1 mu+, 1 mu-, 1 K+, 1 K-
// To save more storing space: will apply cut: 
// - (if there are more combinations having such final state)
//   pairwise dot product between all these final state particles. Should be positive

#include "Pythia8/Pythia.h"
#include "Pythia8/Pythia.h"
#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/HepMC2.h"
#include <cmath>

using namespace Pythia8;

int main() {
    //========================= OUTPUT format ===================================================
    // Interface for conversion from Pythia8::Event to HepMC event.
    HepMC::Pythia8ToHepMC ToHepMC;
    // Specify file where HepMC events will be stored.
    HepMC::IO_GenEvent ascii_io("./hepmc_data/ee2Z2c_comb_cutted_50M_seed9.hepmc", std::ios::out);
    //========================= Settings =========================================================
    // Generator. Process selection. LHC initialization. Histogram.
    Pythia pythia;
    pythia.readString("Beams:eCM = 91.2");    //Z-pole
    pythia.readString("Beams:idA = 11");    //e-
    pythia.readString("Beams:idB = -11");   //e+
    pythia.readString("Random:setSeed =on");
    pythia.readString("Random:seed = 9");
    pythia.readString("WeakSingleBoson:ffbar2gmZ=on");
    pythia.readString("WeakSingleBoson:ffbar2gmZ=on");
    //Z^0 > c, cbar
    pythia.readString("23:onMode = off");
    pythia.readString("23:addChannel = 1 1 100 4 -4");
  
    // https://indico.cern.ch/event/1202105/contributions/5385371/attachments/2660513/4609816/BeamInstrumentation_studies_FCCee.pdf
    // page 18  
    pythia.readString("Beams:allowVertexSpread = on");
    pythia.readString("Beams:sigmaVertexZ = 0.008");
    pythia.readString("Beams:sigmaVertexX = 3.4e-5");
    pythia.readString("Beams:sigmaVertexY = 3.4e-5");

    Event& event = pythia.event; // Shorthand: pythia event record
    int nEvent = 50000000;
    //pythia.mode("Main:numberOfEvents");
    //  int Abort = pythia.mode("Main:timesAllowErrors");
    pythia.readString("PartonLevel:MPI = off");//turn of multi-particle interactions
    pythia.readString("Next:numberShowEvent = 2");
    pythia.init();

    int Ntarget = 0;
    int Flagtarget1 = 0;
    int Flagtarget2 = 0;
    int Flagtarget3 = 0;
    int Flagtarget4 = 0;
    int Flagtarget5 = 0;
    for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
        // if (iEvent%10000 == 0) cout << float(iEvent)/float(nEvent) <<endl;
        if (!pythia.next()) continue;
        Flagtarget1 = 0;
        Flagtarget2 = 0;
        Flagtarget3 = 0;
        Flagtarget4 = 0;
        Flagtarget5 = 0;
        //Count decaying Particle number
        vector<int> mu1 = {};
        vector<int> mu2 = {};
        vector<int> K1 = {};
        vector<int> K2 = {};
        for (int ipart = 0; ipart < event.size(); ipart++) {
            if (event[ipart].id() == 13) {
                Flagtarget1 = 1;
                mu1.push_back(ipart);
            }
            if (event[ipart].id() == -13) {
                Flagtarget2 = 1;
                mu2.push_back(ipart);
            }
            if (event[ipart].id() == 321) {
                Flagtarget3 = 1;
                K1.push_back(ipart);
            }
            if (event[ipart].id() == -321) {
                Flagtarget4 = 1;
                K2.push_back(ipart);
            }
        }

        for(const int& mp1 : mu1){ 
            for(const int& mp2 : mu2){ 
                if (event[mp1].px()*event[mp2].px() + event[mp1].py()*event[mp2].py() + event[mp1].pz()*event[mp2].pz() < 0) continue;
                for(const int& Kp1 : K1){ 
                    if (event[mp1].px()*event[Kp1].px() + event[mp1].py()*event[Kp1].py() + event[mp1].pz()*event[Kp1].pz() < 0 || 
                        event[mp2].px()*event[Kp1].px() + event[mp2].py()*event[Kp1].py() + event[mp2].pz()*event[Kp1].pz() < 0) continue;
                    for(const int& Kp2 : K2){ 
                        if (event[mp1].px()*event[Kp2].px() + event[mp1].py()*event[Kp2].py() + event[mp1].pz()*event[Kp2].pz() < 0 || 
                            event[mp2].px()*event[Kp2].px() + event[mp2].py()*event[Kp2].py() + event[mp2].pz()*event[Kp2].pz() < 0 ||
                            event[Kp1].px()*event[Kp2].px() + event[Kp1].py()*event[Kp2].py() + event[Kp1].pz()*event[Kp2].pz() < 0) continue;
                        Flagtarget5 = 1;
                    }
                }
            }
        }



        if (Flagtarget1 == 1 && Flagtarget2 == 1 && Flagtarget3 == 1 && Flagtarget4 == 1 && Flagtarget5 == 1) {
            Ntarget++;
            HepMC::GenEvent* hepmcevt = new HepMC::GenEvent();
            ToHepMC.fill_next_event(pythia, hepmcevt);
            ascii_io << hepmcevt;
            delete hepmcevt;
        }
        // End of event loop. Statistics. Histogram.
    }
    //  pythia.stat();
    cout << "*********************============================************************" << endl;
    cout << "Simulation finished with  " << nEvent << "  events, Ntarget=  " << Ntarget << endl;
    return 0;
}


