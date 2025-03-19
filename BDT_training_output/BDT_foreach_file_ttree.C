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

  ifstream file_stream("/home/awesole/forDDbarAnalysis/dataFiles/dataFilesTTrees.txt");
  string filename;
  int ifile = 0;

  int ptMin=1;
  int ptMax=20;
  int centMax = 2*90;

  int istart = atoi(i_start);
  int iend = atoi(i_end);

  Double_t MVA[MAXD0CandSize];
  float mass[MAXD0CandSize], phi[MAXD0CandSize], pT[MAXD0CandSize], y[MAXD0CandSize], x3DDecayLengthSignificance[MAXD0CandSize],VtxProb[MAXD0CandSize],x3DPointingAngle[MAXD0CandSize],Dtrk1Pt[MAXD0CandSize],Dtrk2Pt[MAXD0CandSize],Dtrk1PtErr[MAXD0CandSize],Dtrk2PtErr[MAXD0CandSize],xyDCASignificanceDaugther1[MAXD0CandSize],xyDCASignificanceDaugther2[MAXD0CandSize],zDCASignificanceDaugther1[MAXD0CandSize],zDCASignificanceDaugther2[MAXD0CandSize] ;
  int candSize, centrality;

  string output_name = "/scratch/negishi/awesole/dataWithBdt/ROOT/outfile_ttree_"+to_string(istart)+"_"+to_string(iend)+".root";
  //string output_name = "test.root";
  TFile *outf = new TFile(output_name.c_str(),"recreate");
  //TTree *mvaTree = nullptr;
  TTree *mvaTree = new TTree("mvaTree","MVA");


  mvaTree->Branch("candSize",&candSize,"candSize/I");
  mvaTree->Branch("mass", mass, "mass[candSize]/F");
  mvaTree->Branch("phi", phi, "phi[candSize]/F");
  mvaTree->Branch("centrality", &centrality, "centrality/I");
  mvaTree->Branch("pT", pT, "pT[candSize]/F");
  mvaTree->Branch("y", y, "y[candSize]/F");
  mvaTree->Branch("3DDecayLengthSignificance", x3DDecayLengthSignificance, "3DDecayLengthSignificance[candSize]/F");
  mvaTree->Branch("VtxProb",VtxProb , "VtxProb[candSize]/F");
  mvaTree->Branch("3DPointingAngle", x3DPointingAngle,"3DPointingAngle[candSize]/F");
  mvaTree->Branch("pTD1", Dtrk1Pt,"pTD1[candSize]/F");
  mvaTree->Branch("pTD2", Dtrk2Pt,"pTD2[candSize]/F");
  mvaTree->Branch("pTerrD1", Dtrk1PtErr,"pTerrD1[candSize]/F");
  mvaTree->Branch("pTerrD2", Dtrk2PtErr,"pTerrD2[candSize]/F");
  mvaTree->Branch("xyDCASignificanceDaugther1", xyDCASignificanceDaugther1,"xyDCASigD1[candSize]/F");
  mvaTree->Branch("xyDCASignificanceDaugther2", xyDCASignificanceDaugther2,"xyDCASigD2[candSize]/F");
  mvaTree->Branch("zDCASignificanceDaugther1", zDCASignificanceDaugther1,"zDCASigD1[candSize]/F");
  mvaTree->Branch("zDCASignificanceDaugther1", zDCASignificanceDaugther2,"zDCASigD2[candSize]/F");
  mvaTree->Branch("BDT_weight", MVA, "BDT_weight[candSize]/D");

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
    if (fin->IsZombie() || fin->Get("d0ana/VertexCompositeNtuple") == nullptr)
    {
      cout << "BIG PROBLEM" << endl;
      ifile++;
      continue;
    }


    TTree *t = (TTree *)fin->Get("d0ana/VertexCompositeNtuple");
    t->SetBranchAddress("candSize", &candSize);
    t->SetBranchAddress("mass", mass);
    t->SetBranchAddress("phi", phi);
    t->SetBranchAddress("centrality", &centrality); // CHECK
    t->SetBranchAddress("pT", pT);
    t->SetBranchAddress("y", y);
    t->SetBranchAddress("3DDecayLengthSignificance", x3DDecayLengthSignificance);
    t->SetBranchAddress("VtxProb", VtxProb);
    t->SetBranchAddress("3DPointingAngle", x3DPointingAngle);
    t->SetBranchAddress("pTD1", Dtrk1Pt);
    t->SetBranchAddress("pTD2", Dtrk2Pt);
    t->SetBranchAddress("pTerrD1", Dtrk1PtErr);
    t->SetBranchAddress("pTerrD2", Dtrk2PtErr);
    t->SetBranchAddress("xyDCASignificanceDaugther1", xyDCASignificanceDaugther1);
    t->SetBranchAddress("xyDCASignificanceDaugther2", xyDCASignificanceDaugther2);
    t->SetBranchAddress("zDCASignificanceDaugther1", zDCASignificanceDaugther1);


    for(int i_entry=0;i_entry<t->GetEntries();i_entry++)
    //cout << "-------- tree entries = " << t->GetEntries() << endl;
    //for(int i_entry=0;i_entry<100;i_entry++)
    {

      t->GetEntry(i_entry);
      
      if(i_entry%5000==0) cout <<i_entry<<" / "<<t->GetEntries()<< "  "<<100*i_entry/t->GetEntries()<<"%"<<endl;
      
      for(int j=0;j<candSize;j++)
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
      }//for candSize
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
