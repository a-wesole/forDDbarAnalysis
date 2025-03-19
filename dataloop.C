//THIS PROGRAM GENERATES SIGNAL SAMPLE FOR BDT TRAINING
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <vector>
#include <TMath.h>
#include<fstream>
#include<string>
#include <iostream>
#include <TH1F.h>
Bool_t istest=false;
using namespace std;

void mcloop(float ptmin = 1, float ptmax =2 , int centmin = 10 , int centmax = 20, int istart=1, int iend=2 )
{

    ifstream file_stream("/home/awesole/forDDbarAnalysis/dataFiles/dataFilesTTrees.txt"); ////////////////////////////////////CHANGE
    TString outfile = TString::Format("/scratch/negishi/awesole/sidebands_complete_original_production/ROOT/sidebands_pt_%g_%g_cent_%d_%d_file_%d_%d.root",ptmin, ptmax, centmin, centmax, istart, iend);
    TFile results(outfile, "recreate");
    TTree *newtree = new TTree("T" , "Sidebands Sample");
    string filename;

    const float amass = 1.86648;
    int  nDsize=0, ev=0, ifile=0;
    int counts = 0;
    const int max_size =10000;
    int Dsize, centrality;
    float mass[max_size], phi[max_size], pT[max_size], y[max_size], x3DDecayLengthSignificance[max_size],VtxProb[max_size],x3DPointingAngle[max_size],Dtrk1Pt[max_size],Dtrk2Pt[max_size],Dtrk1PtErr[max_size],Dtrk2PtErr[max_size],Dtrk1Chi2n[max_size],Dtrk2Chi2n[max_size],xyDCASignificanceDaugther1[max_size],xyDCASignificanceDaugther2[max_size],zDCASignificanceDaugther1[max_size],zDCASignificanceDaugther2[max_size] ;

    int ncent;
    float nmass[max_size], nphi[max_size], npT[max_size], ny[max_size], n3DDecayLengthSignificance[max_size], nVtxProb[max_size],n3DPointingAngle[max_size],nDtrk1Pt[max_size],nDtrk2Pt[max_size],nDtrk1PtErr[max_size],nDtrk2PtErr[max_size],nDtrk1Chi2n[max_size],nDtrk2Chi2n[max_size], nxyDCASigD1[max_size], nxyDCASigD2[max_size], nzDCASigD1[max_size], nzDCASigD2[max_size], log3ddls[max_size];

    if(istest)
    {
        outfile = "test.root";
        ptmin = 1;
        ptmax=2;
        centmin =0;
        centmax =10;
    }

    // open the corresponding MC file to get entries and fit the mass for sigma
    TString mcfilename;
    mcfilename.Form("/scratch/negishi/awesole/mcFilesForTraining/ROOT/mcOut_pt_%g_%g_cent_%d_%d.root", ptmin, ptmax, centmin, centmax);
    TFile *mc_file = TFile::Open(mcfilename);
    if (!mc_file || mc_file->IsZombie()) {
        std::cerr << "Error: Could not open file " << mcfilename << std::endl;
        return;
    }

    TTree* mctree = (TTree*)mc_file->Get("T");
    const float max_entries = (10 * mctree->GetEntries()) / 5; //divide by 5 because we are going to break list of files in 5 parts 


    TH1F *h = new TH1F("h", "Mass Distribution", 100, 1.7, 2.1);
    mctree->Draw("nmass>>h");

    double sd = h->GetStdDev();
    cout << "sd = " << sd << endl;
    
    //CREATE A NEW FILE AND CLONE OF TREE IN IT
    newtree->Branch("nDsize",&nDsize,"nDsize/I");
    newtree->Branch("nmass", nmass, "nmass[nDsize]/F");
    newtree->Branch("nphi", nphi, "nphi[nDsize]/F");
    newtree->Branch("ncent", &ncent, "ncent/I");
    newtree->Branch("event", &ev, "ev/I");
    newtree->Branch("npT", npT, "npT[nDsize]/F");
    newtree->Branch("ny", ny, "ny[nDsize]/F");
    newtree->Branch("n3DDecayLengthSignificance", n3DDecayLengthSignificance, "n3DDecayLengthSignificance[nDsize]/F");
    newtree->Branch("log3ddls", log3ddls, "log3ddls[nDsize]/F");
    newtree->Branch("nVtxProb",nVtxProb , "nVtxProb[nDsize]/F");
    newtree->Branch("n3DPointingAngle", n3DPointingAngle,"n3DPointingAngle[nDsize]/F");
    newtree->Branch("nDtrk1Pt", nDtrk1Pt,"nDtrk1Pt[nDsize]/F");
    newtree->Branch("nDtrk2Pt", nDtrk2Pt,"nDtrk2Pt[nDsize]/F");
    newtree->Branch("nDtrk1PtErr", nDtrk1PtErr,"nDtrk1PtErr[nDsize]/F");
    newtree->Branch("nDtrk2PtErr", nDtrk2PtErr,"nDtrk2PtErr[nDsize]/F");
    newtree->Branch("nxyDCASigD1", nxyDCASigD1,"nxyDCASigD1[nDsize]/F");
    newtree->Branch("nxyDCASigD2", nxyDCASigD2,"nxyDCASigD2[nDsize]/F");
    newtree->Branch("nzDCASigD1", nzDCASigD1,"nzDCASigD1[nDsize]/F");
    newtree->Branch("nzDCASigD2", nzDCASigD2,"nzDCASigD2[nDsize]/F");


    while(true) {

        if (counts > max_entries) break;
        file_stream >> filename;
        if(file_stream.eof()) break;
        if (ifile < istart) {
            ifile++;
            continue;
        }
        if (ifile >= iend) break;


        //cout << "  " << endl;
        //cout << "-- -- new file -- --" << endl;
        //cout << "ifile = " <<  ifile <<" / 100      "<<100*ifile/100<<"%"<<endl;
        cout << "-- -- -- -- -- -- --" << endl;
        cout << "ifile = " <<  ifile << endl;
        cout << "= " <<  ifile-istart <<" / 1003      "<<100*(ifile-istart)/1003<<"%"<<endl;
        //cout << "-- -- -- -- -- -- --" << endl;
        //cout << "file path = " << filename.c_str() << endl;


        //if (istest && ifile>1) break; 

        TFile *fin = TFile::Open(filename.c_str());
        if(fin->IsZombie() || fin->Get("d0ana/VertexCompositeNtuple")==nullptr) {
            ifile ++;
            continue;
        }


        TTree* tree = (TTree*)fin->Get("d0ana/VertexCompositeNtuple");
        TTree *friendTree = (TTree*)fin->Get("eventinfoana/EventInfoNtuple");
        tree->AddFriend(friendTree);

        tree->SetBranchStatus("*", 0); // Disable all branches
        tree->SetBranchStatus("candSize", 1);
        tree->SetBranchStatus("mass", 1);
        tree->SetBranchStatus("phi", 1);
        tree->SetBranchStatus("centrality", 1); //CHECK
        tree->SetBranchStatus("pT", 1);
        tree->SetBranchStatus("y", 1);
        tree->SetBranchStatus("3DDecayLengthSignificance", 1);
        tree->SetBranchStatus("VtxProb", 1);
        tree->SetBranchStatus("3DPointingAngle", 1);
        tree->SetBranchStatus("pTD1", 1);
        tree->SetBranchStatus("pTD2", 1);
        tree->SetBranchStatus("pTerrD1", 1);
        tree->SetBranchStatus("pTerrD2", 1);
        tree->SetBranchStatus("xyDCASignificanceDaugther1", 1);
        tree->SetBranchStatus("xyDCASignificanceDaugther2", 1);
        tree->SetBranchStatus("zDCASignificanceDaugther1", 1);
        tree->SetBranchStatus("zDCASignificanceDaugther2", 1);

        //for all entries in the tree
        tree->SetBranchAddress("candSize", &Dsize);
        tree->SetBranchAddress("mass", mass);
        tree->SetBranchAddress("phi", phi);
        tree->SetBranchAddress("centrality", &centrality); //CHECK
        tree->SetBranchAddress("pT", pT);
        tree->SetBranchAddress("y", y);
        tree->SetBranchAddress("3DDecayLengthSignificance", x3DDecayLengthSignificance);
        tree->SetBranchAddress("VtxProb", VtxProb);
        tree->SetBranchAddress("3DPointingAngle", x3DPointingAngle);
        tree->SetBranchAddress("pTD1", Dtrk1Pt);
        tree->SetBranchAddress("pTD2", Dtrk2Pt);
        tree->SetBranchAddress("pTerrD1", Dtrk1PtErr);
        tree->SetBranchAddress("pTerrD2", Dtrk2PtErr);
        tree->SetBranchAddress("xyDCASignificanceDaugther1", xyDCASignificanceDaugther1);
        tree->SetBranchAddress("xyDCASignificanceDaugther2", xyDCASignificanceDaugther2);
        tree->SetBranchAddress("zDCASignificanceDaugther1", zDCASignificanceDaugther1);
        tree->SetBranchAddress("zDCASignificanceDaugther2", zDCASignificanceDaugther2);


        float nentries = tree->GetEntries();


        for (int i=0; i<nentries; i++)
            //for (int i=0; i<3000; i++)
        {
            tree->GetEntry(i);
            ev = ev+1;
            ncent=centrality;
            nDsize=0;
            /*
            if (i % 1000 == 0){
                cout <<i<<" / "<<tree->GetEntries()<< "  "<<100*i/tree->GetEntries()<<"%"<<endl;
            }
            */

            if (counts > max_entries) break;

            if (centrality >= 2*centmax || centrality < 2*centmin) continue;


            for (int j=0; j<Dsize; j++)
            {

                if (pT[j] < ptmin || pT[j] >= ptmax) continue;

                tree->GetEntry(i);

                if (TMath::IsNaN(x3DDecayLengthSignificance[j])) continue;
                if (TMath::IsNaN(nxyDCASigD2[j])) continue;
                if (TMath::IsNaN(nzDCASigD2[j])) continue;
                if (((mass[j] > amass-6*sd && mass[j] < amass-3*sd) || (mass[j] < amass+6*sd && mass[j]>amass+3*sd)) && pT[j]>=ptmin && pT[j] <ptmax && y[j]>-1 && y[j] <1){
                    counts += 1;
                    nDsize = nDsize+1;
                    nmass[nDsize-1] =  mass[j];
                    nphi[nDsize-1] =  phi[j];
                    npT[nDsize-1] =  pT[j];
                    ny[nDsize-1] =  y[j];
                    n3DDecayLengthSignificance[nDsize-1] =  x3DDecayLengthSignificance[j];
                    log3ddls[nDsize-1] = log10(x3DDecayLengthSignificance[j]);
                    nVtxProb[nDsize-1] =  VtxProb[j];
                    n3DPointingAngle[nDsize-1] =  x3DPointingAngle[j];
                    nDtrk1Pt[nDsize-1] =  Dtrk1Pt[j];
                    nDtrk2Pt[nDsize-1] =  Dtrk2Pt[j];
                    nDtrk1PtErr[nDsize-1] =  Dtrk1PtErr[j];
                    nDtrk2PtErr[nDsize-1] =  Dtrk2PtErr[j];
                    //nDtrk1Chi2n[nDsize-1] =  Dtrk1Chi2n[j];
                    //nDtrk2Chi2n[nDsize-1] =  Dtrk2Chi2n[j];
                    nxyDCASigD1[nDsize-1] =  xyDCASignificanceDaugther1[j];
                    nxyDCASigD2[nDsize-1] =  xyDCASignificanceDaugther2[j];
                    nzDCASigD1[nDsize-1] =  zDCASignificanceDaugther1[j];
                    nzDCASigD2[nDsize-1] =  zDCASignificanceDaugther2[j];

                }//if dgen==23333

            }//for dsize
            if(nDsize>0){
                newtree->Fill();
            }//if dsize>0
        }//for entries
        //cout << "counts / max_entries = " << counts << " / " << max_entries << endl;
        cout << "counts / max_entries = " << counts << " / " << max_entries << " (" << 100 * (counts / max_entries) << "%)" << endl;

        ifile++;
        fin->Close();
    } // while
    results.cd();
    newtree->Write(0,TObject::kOverwrite);
    results.Write();
    results.Close();
    cout << "DONE" << endl;
}//program

int main(int argc, char *argv[])
{
    if (argc == 7) // Expecting 5 arguments + program name
    {
        float ptmin = atof(argv[1]);
        float ptmax = atof(argv[2]);
        int centmin = atoi(argv[3]);
        int centmax = atoi(argv[4]);
        int istart = std::stoi(argv[5]);
        int iend = std::stoi(argv[6]);

        mcloop(ptmin, ptmax, centmin, centmax, istart, iend);
    }
    else
    {
        std::cout << "Usage: mergeForest <ptmin> <ptmax> <centmin> <centmax> <istart> <iend>" << std::endl;
        return 1;
    }
    return 0;
}
