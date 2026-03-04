#include "TRandom.h"
#include "TApplication.h"
#include "TFile.h"
#include <iostream>
#include <random>
#include <TStreamerInfo.h> 
#include <TClass.h>
#include "TFile.h"
#include "TNtuple.h"
#include "TTree.h"
#include <TBranch.h>
#include <vector>
#include <TRandom3.h>
#include "TROOT.h"
#include "TH1D.h"
#include "TH1.h"
#include "TH2D.h"
#include "TH2F.h"
#include "TMath.h"
#include <TPDF.h>
#include "TF1.h"
#include "TF2.h"
#include "TCanvas.h"
#include <deque>
#include <fstream>

#include "phiTools.h"

const int k_Multiplicity = 10;
const int k_evt_max = 20000;
const int k_Num_Buffer_Size = 3; // number of previous events kept and calculate dphi with

using namespace std;

struct DCand
{
    float phi;
    float flavor; // 1 for d0 and -1 for dbar
};
struct EventPool
{
    std::vector<DCand> d;    // point is we keep only the d0 candidates for mixing in a vector and the Dcands are very light weight with only the variables above.
    std::vector<DCand> dbar; // point is we keep only the dbar candidates for mixing in a vector and the Dcands are very light weight with only the variables above.
};

void sample_hist()
{

    // --- open file ---
    TFile *f = TFile::Open("EM_output/out_40_bin_20_only.root");
    if (!f || f->IsZombie())
    {
        std::cout << "Error opening file" << std::endl;
        return;
    }

    // --- get histogram ---
    TH1 *h1 = (TH1 *)f->Get("h_single_phi");
    TH1 *h = (TH1 *)h1->Clone("h_single_phi_clone");
    // f->Close();
    cout << "h entries = " << h->GetEntries() << endl;

    // --- random generator ---
    TRandom3 rng(0); // seed to 0 for now
    TH1D *h_dphi_mixing_total = new TH1D("h_dphi_mixing_total", "h_dphi_mixing_total", 20, -1 * TMath::Pi() / 2, 3 * TMath::Pi() / 2);
    h_dphi_mixing_total->SetMinimum(0);
    TH1D *h_single_phi_out = new TH1D("h_single_phi_out", "h_single_phi_out", 20, -1 * TMath::Pi(), TMath::Pi());
    h_single_phi_out->SetMinimum(0);
    TH1D *h_single_phi_d = new TH1D("h_single_phi_d", "h_single_phi_d", 20, -1 * TMath::Pi(), TMath::Pi());
    h_single_phi_d->SetMinimum(0);
    TH1D *h_single_phi_db = new TH1D("h_single_phi_db", "h_single_phi_db", 20, -1 * TMath::Pi(), TMath::Pi());
    h_single_phi_db->SetMinimum(0);
    
    TFile *results = new TFile("toy_out_nohole.root", "recreate");
    TNtuple *t_sameEvt = new TNtuple("t_sameEvt", "t_sameEvt", "d0_phi:d0bar_phi:delta_phi");
    t_sameEvt->SetDirectory(results); 
    TNtuple *t_hole = new TNtuple("t_hole", "t_hole", "d0_phi:d0bar_phi:delta_phi");
    t_hole->SetDirectory(results); 
    // TNtuple *t_no_hole = new TNtuple("t_no_hole", "t_no_hole", "d0_phi:d0bar_phi:delta_phi");

    std::deque<EventPool> buffers[1];


    float flavor = 0;

    for (int evt = 0; evt < k_evt_max; evt++)
    {
        if (evt % 1000 == 0)
        {
            cout << "evt = " << evt << endl;
        }
        std::vector<DCand> dCands_inEvent;
        std::vector<DCand> dbarCands_inEvent;
        dCands_inEvent.reserve(k_Multiplicity);
        dbarCands_inEvent.reserve(k_Multiplicity);
        for (int i = 0; i < k_Multiplicity; i++)
        {
            // -- sample from histogram --
            double sampled_phi = h->GetRandom();

            // if (!(sampled_phi >= 0 && sampled_phi < 1))
            // if phi is not in the hole
            if (sampled_phi != 12.2)
            { 
                h_single_phi_out->Fill(sampled_phi);

                // uniform random number in (0,1)
                double r = rng.Uniform(0.0, 1.0);
                if (r < 0.5)
                    flavor = 1; // assign flavor randomly
                else
                    flavor = -1;

                DCand c{sampled_phi, flavor};
                if (flavor == 1)
                {
                    dCands_inEvent.push_back(c);
                    h_single_phi_d->Fill(sampled_phi);
                }
                else if (flavor == -1)
                {
                    dbarCands_inEvent.push_back(c);
                    h_single_phi_db->Fill(sampled_phi);
                }
            }

        } // N loop
        // f->Close();
        //--for mixing--
        for (const auto &curr_d : dCands_inEvent)
        {
            for (const auto &curr_db : dbarCands_inEvent)
            {
                float raw_phi = curr_d.phi - curr_db.phi;
                double dphi = transition_phi(raw_phi);
                // h_dphi_mixing_total->Fill(dphi);
                t_sameEvt->Fill(curr_d.phi, curr_db.phi, dphi);
            }
        }
        for (const auto &pool : buffers[0])
        {                                        // loop over each previous event in the buffer for this exact bin
            for (const auto &d : dCands_inEvent) // for all the dcands in this event
            {
                for (const auto &db_prev : pool.dbar) // for all the dbars in all the previous events
                {
                    float raw_phi = d.phi - db_prev.phi;
                    // if (d.phi > 0 && d.phi < 1) cout << "hey! NO! d.phi =" << d.phi << endl;
                    // if (db_prev.phi > 0 && db_prev.phi < 1) cout << "hey! NO! dbar_prev.phi =" << db_prev.phi << endl;
                    double dphi = transition_phi(raw_phi);
                    h_dphi_mixing_total->Fill(dphi);
                    t_hole->Fill(d.phi, db_prev.phi, dphi);
                }
            }
            for (const auto &db : dbarCands_inEvent)
            { // for all current dbars
                for (const auto &d_prev : pool.d)
                { // loop over all the previous d0s
                    float raw_phi = d_prev.phi - db.phi;
                    // if (d_prev.phi > 0 && d_prev.phi < 1) cout << "hey! NO! d_prev.phi =" << d_prev.phi << endl;
                    // if (db.phi > 0 && db.phi < 1) cout << "hey! NO! db.phi =" << db.phi << endl;
                    double dphi = transition_phi(raw_phi);
                    h_dphi_mixing_total->Fill(dphi);
                    t_hole->Fill(d_prev.phi, db.phi, dphi);
                }
            }
        }

        // --write to event pool--
        EventPool currentPool;
        currentPool.d = std::move(dCands_inEvent);
        currentPool.dbar = std::move(dbarCands_inEvent);
        if (!currentPool.dbar.empty() || !currentPool.d.empty())
        {                                                 // only events with at least 1 d0 or dbar
            buffers[0].push_back(std::move(currentPool)); // adds the current event to the end of the buffer for thiis PvtxZ Bin, we add after the mxing to make sure we do not mix events together
            if (static_cast<int>(buffers[0].size()) > k_Num_Buffer_Size)
            {                           // checks if the buffer has grown larger than the allowed number of stored past events k_num_beffer_size
                buffers[0].pop_front(); // if the buffer is now too large it removes the oldest event bc newest evetns get added at the back
            }
        }
    } // evt loop

    results->cd();    h_dphi_mixing_total->Write();
    h_single_phi_out->Write(); 
    h_single_phi_d->Write(); 
    h_single_phi_db->Write(); 
    t_hole->Write();
    t_sameEvt->Write();
    results->Close();
    f->Close();
}
int main(int argc, char *argv[])
{
    TApplication app("app", &argc, argv);
    if (argc != 1)
    {
        std::cerr << "Usage: " << argv[0];
        return 1;
    }

    sample_hist();
    return 0;
}
