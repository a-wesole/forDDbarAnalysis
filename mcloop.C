//THIS PROGRAM GENERATES SIGNAL SAMPLE FOR BDT TRAINING
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <vector>
#include <TMath.h>
#include<fstream>
#include<string>
#include <iostream>
Bool_t istest=false;
using namespace std;

void mcloop(TString outfile1="", float ptmin = 1, float ptmax =2 , int centmin = 10 , int centmax = 20 )
{

    ifstream file_stream("mcFiles.txt");
    TFile results(outfile1, "recreate");
    TTree *newtree = new TTree("T" , "MC Signal Sample");
    string filename;

    int  nDsize=0, ev=0, ifile=0;
    const int max_size =10000;
    bool matchGEN[max_size], isSwap[max_size];
    int Dsize, centrality;
    float mass[max_size], phi[max_size], pT[max_size], y[max_size], x3DDecayLengthSignificance[max_size],VtxProb[max_size],x3DPointingAngle[max_size],Dtrk1Pt[max_size],Dtrk2Pt[max_size],Dtrk1PtErr[max_size],Dtrk2PtErr[max_size],Dtrk1Chi2n[max_size],Dtrk2Chi2n[max_size],xyDCASignificanceDaugther1[max_size],xyDCASignificanceDaugther2[max_size],zDCASignificanceDaugther1[max_size],zDCASignificanceDaugther2[max_size] ;

    int ncent;
    bool nmatchGEN[max_size], nisSwap[max_size];
    float nmass[max_size], nphi[max_size], npT[max_size], ny[max_size], n3DDecayLengthSignificance[max_size], nVtxProb[max_size],n3DPointingAngle[max_size],nDtrk1Pt[max_size],nDtrk2Pt[max_size],nDtrk1PtErr[max_size],nDtrk2PtErr[max_size],nDtrk1Chi2n[max_size],nDtrk2Chi2n[max_size], nxyDCASigD1[max_size], nxyDCASigD2[max_size], nzDCASigD1[max_size], nzDCASigD2[max_size], log3ddls[max_size];

    if(istest)
    {
        outfile1 = "test.root";
        ptmin = 30;
        ptmax=50;
        centmin =0;
        centmax =100;
    }

    //CREATE A NEW FILE AND CLONE OF TREE IN IT
    newtree->Branch("nDsize",&nDsize,"nDsize/I");
    newtree->Branch("nmatchGEN",&nmatchGEN,"nmatchGEN[nDsize]/I");
    newtree->Branch("nisSwap",&nisSwap,"nisSwap[nDsize]/I");
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

        file_stream >> filename;
        if(file_stream.eof()) break;
        cout << "-- new file --" << endl;
        cout << "ifile=" << ifile << endl;
        cout << ifile <<" / 314      "<<100*ifile/314<<"%"<<endl;
        cout << "-- -- -- -- -- --" << endl;
        //cout << "file path = " << filename.c_str() << endl;


        TFile *fin = TFile::Open(filename.c_str());
        if(fin->IsZombie() || fin->Get("d0ana/VertexCompositeNtuple")==nullptr) {
            ifile ++;
            continue;
        }


    TTree* tree = (TTree*)fin->Get("d0ana/VertexCompositeNtuple");
    TTree *friendTree = (TTree*)fin->Get("eventinfoana/EventInfoNtuple");
    tree->AddFriend(friendTree);


    //for all entries in the tree
    tree->SetBranchAddress("candSize", &Dsize);
    tree->SetBranchAddress("matchGEN", matchGEN);
    tree->SetBranchAddress("isSwap", isSwap);
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
        if (i % 3000 == 0){
            cout <<i<<" / "<<tree->GetEntries()<< "  "<<100*i/tree->GetEntries()<<"%"<<endl;
        }
        if (centrality >= 2*centmax || centrality < 2*centmin) continue;
        for (int j=0; j<Dsize; j++)
        {
            if (TMath::IsNaN(x3DDecayLengthSignificance[j])) continue;
            if (TMath::IsNaN(nxyDCASigD2[j])) continue;
            if (TMath::IsNaN(nzDCASigD2[j])) continue;
            if(y[j]>-1 && y[j] <1 && matchGEN[j]==1 && isSwap[j]==0 && pT[j]>=ptmin && pT[j] <ptmax){
                nDsize = nDsize+1;
                nmatchGEN[nDsize-1]= matchGEN[j];
                nisSwap[nDsize-1]= isSwap[j];
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
    ifile++;
    fin->Close();
    } // while
    results.cd();
    newtree->Write(0,TObject::kOverwrite);
    results.Write();
    results.Close();
}//program

int main(int argc, char *argv[])
{
    if (argc == 6) // Expecting 5 arguments + program name
    {
        TString outfile1 = argv[1];
        float ptmin = atof(argv[2]);
        float ptmax = atof(argv[3]);
        int centmin = atoi(argv[4]);
        int centmax = atoi(argv[5]);

        mcloop(outfile1, ptmin, ptmax, centmin, centmax);
    }
    else
    {
        std::cout << "Usage: mergeForest <input_collection> <output_file> <ptmin> <ptmax> <centmin> <centmax>" << std::endl;
        return 1;
    }
    return 0;
}
