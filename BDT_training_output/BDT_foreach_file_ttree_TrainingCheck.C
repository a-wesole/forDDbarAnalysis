#include "header.h"
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <vector>
#include <TMath.h>
#include <iostream>
#include <fstream>
#include "TNtuple.h"
//#define MAX_XB       20000
#define MAXD0CandSize 50000
using namespace std;
//Bool_t istest = true;

int BDT_foreach_file_ttree(TString i_start="", TString i_end="")
{

  ifstream file_stream("/home/awesole/forDDbarAnalysis/dataFiles/");
  //ifstream file_stream("MC.list"); all_data.list//
  string filename;
  int ifile = 0;

  int ptMin=1;
  int ptMax=20;
  int centMax = 2*90;

  int istart = atoi(i_start);
  int iend = atoi(i_end);

  Double_t MVA[MAXD0CandSize];
  float mass[MAXD0CandSize], phi[MAXD0CandSize], pT[MAXD0CandSize], y[MAXD0CandSize], x3DDecayLengthSignificance[MAXD0CandSize],VtxProb[MAXD0CandSize],x3DPointingAngle[MAXD0CandSize],Dtrk1Pt[MAXD0CandSize],Dtrk2Pt[MAXD0CandSize],Dtrk1PtErr[MAXD0CandSize],Dtrk2PtErr[MAXD0CandSize],xyDCASignificanceDaugther1[MAXD0CandSize],xyDCASignificanceDaugther2[MAXD0CandSize],zDCASignificanceDaugther1[MAXD0CandSize],zDCASignificanceDaugther2[MAXD0CandSize] ;
  int Dsize, centrality;

  //string output_name = "/scratch/bell/chand140/ESE_output/data_output_with_BDT_ttree_data/ROOT/outfile_ttree_"+to_string(istart)+"_"+to_string(iend)+".root";
  string output_name = "test.root";
  TFile *outf = new TFile(output_name.c_str(),"recreate");
  //TTree *mvaTree = nullptr;
  TTree *mvaTree = new TTree("mvaTree","MVA");



  mvaTree->Branch("nDsize",&Dsize,"nDsize/I");
  mvaTree->Branch("nmass", mass, "nmass[nDsize]/F");
  mvaTree->Branch("nphi", phi, "nphi[nDsize]/F");
  mvaTree->Branch("ncent", &centrality, "ncent/I");
  mvaTree->Branch("npT", pT, "npT[nDsize]/F");
  mvaTree->Branch("ny", y, "ny[nDsize]/F");
  mvaTree->Branch("n3DDecayLengthSignificance", x3DDecayLengthSignificance, "n3DDecayLengthSignificance[nDsize]/F");
  mvaTree->Branch("nVtxProb",VtxProb , "nVtxProb[nDsize]/F");
  mvaTree->Branch("n3DPointingAngle", x3DPointingAngle,"n3DPointingAngle[nDsize]/F");
  mvaTree->Branch("nDtrk1Pt", Dtrk1Pt,"nDtrk1Pt[nDsize]/F");
  mvaTree->Branch("nDtrk2Pt", Dtrk2Pt,"nDtrk2Pt[nDsize]/F");
  mvaTree->Branch("nDtrk1PtErr", Dtrk1PtErr,"nDtrk1PtErr[nDsize]/F");
  mvaTree->Branch("nDtrk2PtErr", Dtrk2PtErr,"nDtrk2PtErr[nDsize]/F");
  mvaTree->Branch("nxyDCASigD1", xyDCASignificanceDaugther1,"nxyDCASigD1[nDsize]/F");
  mvaTree->Branch("nxyDCASigD2", xyDCASignificanceDaugther2,"nxyDCASigD2[nDsize]/F");
  mvaTree->Branch("nzDCASigD1", zDCASignificanceDaugther1,"nzDCASigD1[nDsize]/F");
  mvaTree->Branch("nzDCASigD2", zDCASignificanceDaugther2,"nzDCASigD2[nDsize]/F");
  mvaTree->Branch("BDT_weight", MVA, "BDT_weight[nDsize]/D");

  //   std::cout << "SET VARIABLE VECTOR" << endl;
	std::vector<std::string> theInputVars;
	string a1="log3ddls"; 
	string a2="nVtxProb"; 
	string a3="n3DPointingAngle"; 
	string a4="nDtrk1Pt"; 
	string a5="nDtrk2Pt"; 
	string a6="nxyDCASigD1";
	string a7="nxyDCASigD2";
	string a8="nzDCASigD1";
	string a9="nzDCASigD2";
  /*
	string a1="log3ddls"; 
	string a2="VtxProb"; 
	string a3="x3DPointingAngle"; 
	string a4="Dtrk1Pt"; 
	string a5="Dtrk2Pt"; 
	string a6="xyDCASignificanceDaugther1";
	string a7="xyDCASignificanceDaugther2";
	string a8="zDCASignificanceDaugther1";
	string a9="zDCASignificanceDaugther2";
  */

	theInputVars.push_back(a1);
	theInputVars.push_back(a2);
	theInputVars.push_back(a3);
	theInputVars.push_back(a4); 
	theInputVars.push_back(a5); 
	theInputVars.push_back(a6);
	theInputVars.push_back(a7);
	theInputVars.push_back(a8);
	theInputVars.push_back(a9);



  vector<double> inputValues;

  ReadBDT_pt_1_2_cent_0_10 mva_pt_1_2_cent_0_10(theInputVars);
  ReadBDT_pt_1_2_cent_10_30 mva_pt_1_2_cent_10_30(theInputVars);
  ReadBDT_pt_1_2_cent_30_50 mva_pt_1_2_cent_30_50(theInputVars);
  ReadBDT_pt_1_2_cent_50_90 mva_pt_1_2_cent_50_90(theInputVars);

  ReadBDT_pt_2_3_cent_0_10 mva_pt_2_3_cent_0_10(theInputVars);
  ReadBDT_pt_2_3_cent_10_30 mva_pt_2_3_cent_10_30(theInputVars);
  ReadBDT_pt_2_3_cent_30_50 mva_pt_2_3_cent_30_50(theInputVars);
  ReadBDT_pt_2_3_cent_50_90 mva_pt_2_3_cent_50_90(theInputVars);

  ReadBDT_pt_3_4_cent_0_10 mva_pt_3_4_cent_0_10(theInputVars);
  ReadBDT_pt_3_4_cent_10_30 mva_pt_3_4_cent_10_30(theInputVars);
  ReadBDT_pt_3_4_cent_30_50 mva_pt_3_4_cent_30_50(theInputVars);
  ReadBDT_pt_3_4_cent_50_90 mva_pt_3_4_cent_50_90(theInputVars);

  ReadBDT_pt_4_5_cent_0_10 mva_pt_4_5_cent_0_10(theInputVars);
  ReadBDT_pt_4_5_cent_10_30 mva_pt_4_5_cent_10_30(theInputVars);
  ReadBDT_pt_4_5_cent_30_50 mva_pt_4_5_cent_30_50(theInputVars);
  ReadBDT_pt_4_5_cent_50_90 mva_pt_4_5_cent_50_90(theInputVars);

  ReadBDT_pt_5_6_cent_0_10 mva_pt_5_6_cent_0_10(theInputVars);
  ReadBDT_pt_5_6_cent_10_30 mva_pt_5_6_cent_10_30(theInputVars);
  ReadBDT_pt_5_6_cent_30_50 mva_pt_5_6_cent_30_50(theInputVars);
  ReadBDT_pt_5_6_cent_50_90 mva_pt_5_6_cent_50_90(theInputVars);

  ReadBDT_pt_6_8_cent_0_10 mva_pt_6_8_cent_0_10(theInputVars);
  ReadBDT_pt_6_8_cent_10_30 mva_pt_6_8_cent_10_30(theInputVars);
  ReadBDT_pt_6_8_cent_30_50 mva_pt_6_8_cent_30_50(theInputVars);
  ReadBDT_pt_6_8_cent_50_90 mva_pt_6_8_cent_50_90(theInputVars);

  ReadBDT_pt_8_10_cent_0_10 mva_pt_8_10_cent_0_10(theInputVars);
  ReadBDT_pt_8_10_cent_10_30 mva_pt_8_10_cent_10_30(theInputVars);
  ReadBDT_pt_8_10_cent_30_50 mva_pt_8_10_cent_30_50(theInputVars);
  ReadBDT_pt_8_10_cent_50_90 mva_pt_8_10_cent_50_90(theInputVars);

  ReadBDT_pt_10_15_cent_0_10 mva_pt_10_15_cent_0_10(theInputVars);
  ReadBDT_pt_10_15_cent_10_30 mva_pt_10_15_cent_10_30(theInputVars);
  ReadBDT_pt_10_15_cent_30_50 mva_pt_10_15_cent_30_50(theInputVars);
  ReadBDT_pt_10_15_cent_50_90 mva_pt_10_15_cent_50_90(theInputVars);

  ReadBDT_pt_15_20_cent_0_10 mva_pt_15_20_cent_0_10(theInputVars);
  ReadBDT_pt_15_20_cent_10_30 mva_pt_15_20_cent_10_30(theInputVars);
  ReadBDT_pt_15_20_cent_30_50 mva_pt_15_20_cent_30_50(theInputVars);
  ReadBDT_pt_15_20_cent_50_90 mva_pt_15_20_cent_50_90(theInputVars);

  while(true)
  {
    file_stream >> filename;


    if(file_stream.eof())
      break;

    if(ifile < istart) {
      ifile++;
      continue;
    }

    if(ifile >= iend)
      break;

    
    TFile *fin = TFile::Open(filename.c_str());
    cout << "file path = " << filename.c_str() << endl;
    if (fin->IsZombie() || fin->Get("T") == nullptr)
    {
      cout << "BIG PROBLEM" << endl;
      ifile++;
      continue;
    }


    /*
    //TTree *t = (TTree *)inf.Get("d0ana/VertexCompositeNtuple");
    tree->SetBranchAddress("candSize", &Dsize);
    tree->SetBranchAddress("mass", mass);
    tree->SetBranchAddress("phi", phi);
    tree->SetBranchAddress("centrality", &centrality); // CHECK
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
    */

    TTree *t = (TTree *)fin->Get("T");
    t->SetBranchAddress("nDsize", &Dsize);
    t->SetBranchAddress("nmass", mass);
    t->SetBranchAddress("nphi", phi);
    t->SetBranchAddress("ncent", &centrality); // CHECK
    t->SetBranchAddress("npT", pT);
    t->SetBranchAddress("ny", y);
    t->SetBranchAddress("n3DDecayLengthSignificance", x3DDecayLengthSignificance);
    t->SetBranchAddress("nVtxProb", VtxProb);
    t->SetBranchAddress("n3DPointingAngle", x3DPointingAngle);
    t->SetBranchAddress("nDtrk1Pt", Dtrk1Pt);
    t->SetBranchAddress("nDtrk2Pt", Dtrk2Pt);
    t->SetBranchAddress("nDtrk1PtErr", Dtrk1PtErr);
    t->SetBranchAddress("nDtrk2PtErr", Dtrk2PtErr);
    t->SetBranchAddress("nxyDCASigD1", xyDCASignificanceDaugther1);
    t->SetBranchAddress("nxyDCASigD2", xyDCASignificanceDaugther2);
    t->SetBranchAddress("nzDCASigD2", zDCASignificanceDaugther1);
    t->SetBranchAddress("nzDCASigD2", zDCASignificanceDaugther2);



    for(int i_entry=0;i_entry<t->GetEntries();i_entry++)
    //for(int i_entry=0;i_entry<10;i_entry++)
    {

      t->GetEntry(i_entry);
      
      if(i_entry%50==0) cout <<i_entry<<" / "<<t->GetEntries()<< "  "<<100*i_entry/t->GetEntries()<<"%"<<endl;
      
      for(int j=0;j<Dsize;j++)
      {
        MVA[j]=-999;//initial value	  

        //if (pT[j]>ptMin && pT[j]<ptMax && Dtrk1Chi2n[j]<0.18 && Dtrk2Chi2n[j]<0.18 && centrality<100)
        if (pT[j] > ptMin && pT[j] < ptMax && centrality < centMax)
        {
          inputValues.clear();
          inputValues.push_back(log10(x3DDecayLengthSignificance[j])); // 00
          inputValues.push_back(VtxProb[j]);                           // 01
          inputValues.push_back(x3DPointingAngle[j]);                  // 02
          inputValues.push_back(Dtrk1Pt[j]);                           // 03
          inputValues.push_back(Dtrk2Pt[j]);                           // 04
          inputValues.push_back(xyDCASignificanceDaugther1[j]);        // 04
          inputValues.push_back(xyDCASignificanceDaugther2[j]);        // 04
          inputValues.push_back(zDCASignificanceDaugther1[j]);         // 04
          inputValues.push_back(zDCASignificanceDaugther2[j]);

          if (centrality>=0 && centrality<20){
            if (pT[j]>=1 && pT[j]< 2) MVA[j]=mva_pt_1_2_cent_0_10.GetMvaValue_pt_1_2_cent_0_10(inputValues);
            if (pT[j]>=2 && pT[j]< 3) MVA[j]=mva_pt_2_3_cent_0_10.GetMvaValue_pt_2_3_cent_0_10(inputValues);
            if (pT[j]>=3 && pT[j]< 4) MVA[j]=mva_pt_3_4_cent_0_10.GetMvaValue_pt_3_4_cent_0_10(inputValues);
            if (pT[j]>=4 && pT[j]< 5) MVA[j]=mva_pt_4_5_cent_0_10.GetMvaValue_pt_4_5_cent_0_10(inputValues);
            if (pT[j]>=5 && pT[j]< 6) MVA[j]=mva_pt_5_6_cent_0_10.GetMvaValue_pt_5_6_cent_0_10(inputValues);
            if (pT[j]>=6 && pT[j]< 8) MVA[j]=mva_pt_6_8_cent_0_10.GetMvaValue_pt_6_8_cent_0_10(inputValues);
            if (pT[j]>=8 && pT[j]< 10) MVA[j]=mva_pt_8_10_cent_0_10.GetMvaValue_pt_8_10_cent_0_10(inputValues);
            if (pT[j]>=10 && pT[j]< 15) MVA[j]=mva_pt_10_15_cent_0_10.GetMvaValue_pt_10_15_cent_0_10(inputValues);
            if (pT[j]>=15 && pT[j]< 20) MVA[j]=mva_pt_15_20_cent_0_10.GetMvaValue_pt_15_20_cent_0_10(inputValues);
          }

          if (centrality>=20 && centrality<60){
            if (pT[j]>=1 && pT[j]< 2) MVA[j]=mva_pt_1_2_cent_10_30.GetMvaValue_pt_1_2_cent_10_30(inputValues);
            if (pT[j]>=2 && pT[j]< 3) MVA[j]=mva_pt_2_3_cent_10_30.GetMvaValue_pt_2_3_cent_10_30(inputValues);
            if (pT[j]>=3 && pT[j]< 4) MVA[j]=mva_pt_3_4_cent_10_30.GetMvaValue_pt_3_4_cent_10_30(inputValues);
            if (pT[j]>=4 && pT[j]< 5) MVA[j]=mva_pt_4_5_cent_10_30.GetMvaValue_pt_4_5_cent_10_30(inputValues);
            if (pT[j]>=5 && pT[j]< 6) MVA[j]=mva_pt_5_6_cent_10_30.GetMvaValue_pt_5_6_cent_10_30(inputValues);
            if (pT[j]>=6 && pT[j]< 8) MVA[j]=mva_pt_6_8_cent_10_30.GetMvaValue_pt_6_8_cent_10_30(inputValues);
            if (pT[j]>=8 && pT[j]< 10) MVA[j]=mva_pt_8_10_cent_10_30.GetMvaValue_pt_8_10_cent_10_30(inputValues);
            if (pT[j]>=10 && pT[j]< 15) MVA[j]=mva_pt_10_15_cent_10_30.GetMvaValue_pt_10_15_cent_10_30(inputValues);
            if (pT[j]>=15 && pT[j]< 20) MVA[j]=mva_pt_15_20_cent_10_30.GetMvaValue_pt_15_20_cent_10_30(inputValues);
          }

          if (centrality>=60 && centrality<100){
            if (pT[j]>=1 && pT[j]< 2) MVA[j]=mva_pt_1_2_cent_30_50.GetMvaValue_pt_1_2_cent_30_50(inputValues);
            if (pT[j]>=2 && pT[j]< 3) MVA[j]=mva_pt_2_3_cent_30_50.GetMvaValue_pt_2_3_cent_30_50(inputValues);
            if (pT[j]>=3 && pT[j]< 4) MVA[j]=mva_pt_3_4_cent_30_50.GetMvaValue_pt_3_4_cent_30_50(inputValues);
            if (pT[j]>=4 && pT[j]< 5) MVA[j]=mva_pt_4_5_cent_30_50.GetMvaValue_pt_4_5_cent_30_50(inputValues);
            if (pT[j]>=5 && pT[j]< 6) MVA[j]=mva_pt_5_6_cent_30_50.GetMvaValue_pt_5_6_cent_30_50(inputValues);
            if (pT[j]>=6 && pT[j]< 8) MVA[j]=mva_pt_6_8_cent_30_50.GetMvaValue_pt_6_8_cent_30_50(inputValues);
            if (pT[j]>=8 && pT[j]< 10) MVA[j]=mva_pt_8_10_cent_30_50.GetMvaValue_pt_8_10_cent_30_50(inputValues);
            if (pT[j]>=10 && pT[j]< 15) MVA[j]=mva_pt_10_15_cent_30_50.GetMvaValue_pt_10_15_cent_30_50(inputValues);
            if (pT[j]>=15 && pT[j]< 20) MVA[j]=mva_pt_15_20_cent_30_50.GetMvaValue_pt_15_20_cent_30_50(inputValues);
          }
          if (centrality>=100 && centrality<180){
            if (pT[j]>=1 && pT[j]< 2) MVA[j]=mva_pt_1_2_cent_50_90.GetMvaValue_pt_1_2_cent_50_90(inputValues);
            if (pT[j]>=2 && pT[j]< 3) MVA[j]=mva_pt_2_3_cent_50_90.GetMvaValue_pt_2_3_cent_50_90(inputValues);
            if (pT[j]>=3 && pT[j]< 4) MVA[j]=mva_pt_3_4_cent_50_90.GetMvaValue_pt_3_4_cent_50_90(inputValues);
            if (pT[j]>=4 && pT[j]< 5) MVA[j]=mva_pt_4_5_cent_50_90.GetMvaValue_pt_4_5_cent_50_90(inputValues);
            if (pT[j]>=5 && pT[j]< 6) MVA[j]=mva_pt_5_6_cent_50_90.GetMvaValue_pt_5_6_cent_50_90(inputValues);
            if (pT[j]>=6 && pT[j]< 8) MVA[j]=mva_pt_6_8_cent_50_90.GetMvaValue_pt_6_8_cent_50_90(inputValues);
            if (pT[j]>=8 && pT[j]< 10) MVA[j]=mva_pt_8_10_cent_50_90.GetMvaValue_pt_8_10_cent_50_90(inputValues);
            if (pT[j]>=10 && pT[j]< 15) MVA[j]=mva_pt_10_15_cent_50_90.GetMvaValue_pt_10_15_cent_50_90(inputValues);
            if (pT[j]>=15 && pT[j]< 20) MVA[j]=mva_pt_15_20_cent_50_90.GetMvaValue_pt_15_20_cent_50_90(inputValues);
          }
        } // if
      }//for Dsize
      mvaTree->Fill();
    }//for entry
    ifile++;
    fin->Close();
  }
  outf->cd();
  mvaTree->Write(0,TObject::kOverwrite);
  outf->Close();

  return 0;

}
int main(int argc, char *argv[])
{
  if(argc==3)
  {
    BDT_foreach_file_ttree(argv[1], argv[2]);
  }
  else
  {
    std::cout << "Use <istart> <iend>" << std::endl;
    return 1;
  }
  return 0;
}
