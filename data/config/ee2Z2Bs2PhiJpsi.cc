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
    HepMC::IO_GenEvent ascii_io("./hepmc_data/ee2Z2Bs2PhiJpsi_1M_seed0.hepmc", std::ios::out);
    //========================= Settings =========================================================
    // Generator. Process selection. LHC initialization. Histogram.
    Pythia pythia;
    pythia.readString("Beams:eCM = 91.2");    //Z-pole
    pythia.readString("Beams:idA = 11");    //e-
    pythia.readString("Beams:idB = -11");   //e+
    pythia.readString("Random:setSeed =on");
    pythia.readString("Random:seed = 0");
    pythia.readString("WeakSingleBoson:ffbar2gmZ=on");
    pythia.readString("WeakSingleBoson:ffbar2gmZ=on");
    //Z^0 > b, bbar
    pythia.readString("23:onMode = off");
    pythia.readString("23:addChannel = 1 1 100 5 -5");
    //Set relevant BR to 100%
    pythia.readString("531:onMode = off");
    pythia.readString("531:addChannel = 1 1 100 333 443");
    pythia.readString("443:onMode = off");
    pythia.readString("443:addChannel = 1 1 100 13 -13");
    pythia.readString("333:onMode = off");
    pythia.readString("333:addChannel = 1 1 100 321 -321");
  
    // https://indico.cern.ch/event/1202105/contributions/5385371/attachments/2660513/4609816/BeamInstrumentation_studies_FCCee.pdf
    // page 18  
    pythia.readString("Beams:allowVertexSpread = on");
    pythia.readString("Beams:sigmaVertexZ = 0.008");
    pythia.readString("Beams:sigmaVertexX = 3.4e-5");
    pythia.readString("Beams:sigmaVertexY = 3.4e-5");

    Event& event = pythia.event; // Shorthand: pythia event record
    int nEvent = 1000000;
    //pythia.mode("Main:numberOfEvents");
    //  int Abort = pythia.mode("Main:timesAllowErrors");
    pythia.readString("PartonLevel:MPI = off");//turn of multi-particle interactions
    pythia.readString("Next:numberShowEvent = 2");
    pythia.init();

    int Ntarget = 0;
    int Flagtarget = 0;
    for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
        // if (iEvent%10000 == 0) cout << float(iEvent)/float(nEvent) <<endl;
        if (!pythia.next()) continue;
        Flagtarget = 0;
        //Count decaying Particle number
        for (int iPart = 0; iPart < event.size(); iPart++) {
            if (abs(event[iPart].idAbs() == 531)) {
                Flagtarget += 1;
            }
        }
        if (Flagtarget != 0) {
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


