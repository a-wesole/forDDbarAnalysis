//THIS PROGRAM GENERATES BDT VALUES FOR A SET OF DATA, OVER A LARGE RANGE OF PT VALUES
//THIS IS THE MOST UPDTED AND INCLUSIVE FILE
//NEED TO INCLUDE THE APPROPRIATE HEADER FILE THAT REFERENCES THE class.C OUTPUT FILE FROM BDT TRAINING
#include "header.h"
//#include "differentialheader.h"
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <vector>
#include <TMath.h>
#include <iostream>
#define MAX_XB       20000
#define MAXD0CandSize 1000
using namespace std;
Bool_t istest = true;
int differentialBDT(TString inputfile="", TString outfile1="")
{

	if(istest)
	{
		inputfile = "root://eos.cms.rcac.purdue.edu//rdue/store/user/wxie/D0_data_TTrees/HIPhysicsRawPrime0/D0_DATA_Run3_TTree/250306_200703/0004/analyzed_TTree_4000.root";

		outfile1 = "test.root";



	}

	TFile* inf = TFile::Open(inputfile);
	TTree* t = (TTree*)inf->Get("d0ana/VertexCompositeNtuple");

	const int max_size=100000;

	int centrality,  Dsize, Dgen[max_size];
	float mass[max_size], pT[max_size], y[max_size], x3DDecayLengthSignificance[max_size], VtxProb[max_size],x3DPointingAngle[max_size],Dtrk1Pt[max_size],Dtrk2Pt[max_size],Dtrk1PtErr[max_size],Dtrk2PtErr[max_size],Dtrk1Chi2n[max_size],Dtrk2Chi2n[max_size],xyDCASignificanceDaugther1[max_size],xyDCASignificanceDaugther2[max_size],zDCASignificanceDaugther1[max_size],zDCASignificanceDaugther2[max_size];
	double MVA[max_size];


	TFile *outf = new TFile(outfile1,"recreate");
	TTree *mvaTree = new TTree("mvaTree","MVA");


	mvaTree->Branch("candSize",&Dsize,"Dsize/I");
	mvaTree->Branch("mass", mass, "mass[Dsize]/F");
	mvaTree->Branch("phi", phi, "phi[Dsize]/F");
	mvaTree->Branch("centrality", &centrality, "centrality/I");
	mvaTree->Branch("pT", pT, "pT[Dsize]/F");
	mvaTree->Branch("y", y, "y[Dsize]/F");
	mvaTree->Branch("pTerrD1", Dtrk1PtErr,"Dtrk1PtErr[Dsize]/F");
	mvaTree->Branch("pTerrD2", Dtrk2PtErr,"Dtrk2PtErr[Dsize]/F");

	mvaTree->Branch("3DDecayLengthSignificance", x3DDecayLengthSignificance, "3DDecayLengthSignificance[Dsize]/F");
	mvaTree->Branch("VtxProb",VtxProb , "VtxProb[Dsize]/F");
	mvaTree->Branch("3DPointingAngle", x3DPointingAngle,"3DPointingAngle[Dsize]/F");
	mvaTree->Branch("pTD1", Dtrk1Pt,"Dtrk1Pt[Dsize]/F");
	mvaTree->Branch("pTD1", Dtrk2Pt,"Dtrk2Pt[Dsize]/F");
	mvaTree->Branch("xyDCASignificanceDaugther1", xyDCASignificanceDaugther1,"xyDCASignificanceDaugther1[Dsize]/F");
	mvaTree->Branch("xyDCASignificanceDaugther2", xyDCASignificanceDaugther2,"xyDCASignificanceDaugther2[Dsize]/F");
	mvaTree->Branch("zDCASignificanceDaugther1", zDCASignificanceDaugther1,"zDCASignificanceDaugther1[Dsize]/F");
	mvaTree->Branch("zDCASignificanceDaugther2", zDCASignificanceDaugther2,"zDCASignificanceDaugther2[Dsize]/F");


	mvaTree->Branch("BDT_weight",MVA,"BDT_weight[Dsize]/D");

	//   std::cout << "SET VARIABLE VECTOR" << endl;
	std::vector<std::string> theInputVars;
	string a1="log3ddls"; 
	string a2="VtxProb"; 
	string a3="x3DPointingAngle"; 
	string a4="Dtrk1Pt"; 
	string a5="Dtrk2Pt"; 
	string a6="xyDCASignificanceDaugther1";
	string a7="xyDCASignificanceDaugther2";
	string a8="zDCASignificanceDaugther1";
	string a9="zDCASignificanceDaugther2";

	theInputVars.push_back(a1);
	theInputVars.push_back(a2);
	theInputVars.push_back(a3);
	theInputVars.push_back(a4); 
	theInputVars.push_back(a5); 
	theInputVars.push_back(a6);
	theInputVars.push_back(a7);
	theInputVars.push_back(a8);
	theInputVars.push_back(a9);



	std::vector<double> inputValues;


	ReadBDTcent_0_10_pt_1_2_ mva_cent_0_10_pt_1_2_(theInputVars);
	ReadBDTcent_0_10_pt_2_3_ mva_cent_0_10_pt_2_3_(theInputVars);
	ReadBDTcent_0_10_pt_3_4_ mva_cent_0_10_pt_3_4_(theInputVars);
	ReadBDTcent_0_10_pt_4_5_ mva_cent_0_10_pt_4_5_(theInputVars);
	ReadBDTcent_0_10_pt_5_6_ mva_cent_0_10_pt_5_6_(theInputVars);
	ReadBDTcent_0_10_pt_6_8_ mva_cent_0_10_pt_6_8_(theInputVars);
	ReadBDTcent_0_10_pt_8_10_ mva_cent_0_10_pt_8_10_(theInputVars);
	ReadBDTcent_0_10_pt_10_15_ mva_cent_0_10_pt_10_15_(theInputVars);
	ReadBDTcent_0_10_pt_15_20_ mva_cent_0_10_pt_15_20_(theInputVars);
	ReadBDTcent_0_10_pt_20_30_ mva_cent_0_10_pt_20_30_(theInputVars);

	ReadBDTcent_10_30_pt_1_2_ mva_cent_10_30_pt_1_2_(theInputVars);
	ReadBDTcent_10_30_pt_2_3_ mva_cent_10_30_pt_2_3_(theInputVars);
	ReadBDTcent_10_30_pt_3_4_ mva_cent_10_30_pt_3_4_(theInputVars);
	ReadBDTcent_10_30_pt_4_5_ mva_cent_10_30_pt_4_5_(theInputVars);
	ReadBDTcent_10_30_pt_5_6_ mva_cent_10_30_pt_5_6_(theInputVars);
	ReadBDTcent_10_30_pt_6_8_ mva_cent_10_30_pt_6_8_(theInputVars);
	ReadBDTcent_10_30_pt_8_10_ mva_cent_10_30_pt_8_10_(theInputVars);
	ReadBDTcent_10_30_pt_10_15_ mva_cent_10_30_pt_10_15_(theInputVars);
	ReadBDTcent_10_30_pt_15_20_ mva_cent_10_30_pt_15_20_(theInputVars);
	ReadBDTcent_10_30_pt_20_30_ mva_cent_10_30_pt_20_30_(theInputVars);

	ReadBDTcent_30_50_pt_1_2_ mva_cent_30_50_pt_1_2_(theInputVars);
	ReadBDTcent_30_50_pt_2_3_ mva_cent_30_50_pt_2_3_(theInputVars);
	ReadBDTcent_30_50_pt_3_4_ mva_cent_30_50_pt_3_4_(theInputVars);
	ReadBDTcent_30_50_pt_4_5_ mva_cent_30_50_pt_4_5_(theInputVars);
	ReadBDTcent_30_50_pt_5_6_ mva_cent_30_50_pt_5_6_(theInputVars);
	ReadBDTcent_30_50_pt_6_8_ mva_cent_30_50_pt_6_8_(theInputVars);
	ReadBDTcent_30_50_pt_8_10_ mva_cent_30_50_pt_8_10_(theInputVars);
	ReadBDTcent_30_50_pt_10_15_ mva_cent_30_50_pt_10_15_(theInputVars);
	ReadBDTcent_30_50_pt_15_20_ mva_cent_30_50_pt_15_20_(theInputVars);
	ReadBDTcent_30_50_pt_20_30_ mva_cent_30_50_pt_20_30_(theInputVars);

	ReadBDTcent_50_90_pt_1_2_ mva_cent_50_90_pt_1_2_(theInputVars);
	ReadBDTcent_50_90_pt_2_3_ mva_cent_50_90_pt_2_3_(theInputVars);
	ReadBDTcent_50_90_pt_3_4_ mva_cent_50_90_pt_3_4_(theInputVars);
	ReadBDTcent_50_90_pt_4_5_ mva_cent_50_90_pt_4_5_(theInputVars);
	ReadBDTcent_50_90_pt_5_6_ mva_cent_50_90_pt_5_6_(theInputVars);
	ReadBDTcent_50_90_pt_6_8_ mva_cent_50_90_pt_6_8_(theInputVars);
//	ReadBDT_cent_50_90_pt_8_30_ mva_cent_50_90_pt_8_30_(theInputVars);
	t->SetBranchStatus("Dsize", 1);
	t->SetBranchStatus("mass", 1);
	t->SetBranchStatus("centrality", 1);
	t->SetBranchStatus("pT", 1);
	t->SetBranchStatus("y", 1);
	t->SetBranchStatus("Dtrk1PtErr", 1);
	t->SetBranchStatus("Dtrk2PtErr", 1);
	t->SetBranchStatus("Dtrk1Chi2n", 1);
	t->SetBranchStatus("Dtrk2Chi2n", 1);
	t->SetBranchStatus("3DDecayLengthSignificance", 1);
	t->SetBranchStatus("VtxProb", 1);
	t->SetBranchStatus("3DPointingAngle", 1);
	t->SetBranchStatus("Dtrk1Pt", 1);
	t->SetBranchStatus("Dtrk2Pt", 1);
	t->SetBranchStatus("xyDCASignificanceDaugther1", 1);
	t->SetBranchStatus("xyDCASignificanceDaugther2", 1);
	t->SetBranchStatus("zDCASignificanceDaugther1", 1);
	t->SetBranchStatus("zDCASignificanceDaugther2", 1);


	t->SetBranchAddress("Dsize", &Dsize);
	t->SetBranchAddress("mass", mass);
	t->SetBranchAddress("centrality", &centrality);
	t->SetBranchAddress("pT", pT);
	t->SetBranchAddress("y", y);
	t->SetBranchAddress("Dtrk1PtErr", Dtrk1PtErr);
	t->SetBranchAddress("Dtrk2PtErr", Dtrk2PtErr);
	t->SetBranchAddress("Dtrk1Chi2n", Dtrk1Chi2n);
	t->SetBranchAddress("Dtrk2Chi2n", Dtrk2Chi2n);
	t->SetBranchAddress("3DDecayLengthSignificance", x3DDecayLengthSignificance);
	t->SetBranchAddress("VtxProb", VtxProb);
	t->SetBranchAddress("3DPointingAngle", x3DPointingAngle);
	t->SetBranchAddress("Dtrk1Pt", Dtrk1Pt);
	t->SetBranchAddress("Dtrk2Pt", Dtrk2Pt);
	t->SetBranchAddress("xyDCASignificanceDaugther1", xyDCASignificanceDaugther1);
	t->SetBranchAddress("xyDCASignificanceDaugther2", xyDCASignificanceDaugther2);
	t->SetBranchAddress("zDCASignificanceDaugther1", zDCASignificanceDaugther1);
	t->SetBranchAddress("zDCASignificanceDaugther2", zDCASignificanceDaugther2);


	for (int i=0; i<t->GetEntries(); i++)
    //for (int i=0; i<100000; i++)
	//for (int i=0; i<100; i++)
	{
		t->GetEntry(i);
		if (i % 100000 == 0){
		//if (i % 1 == 0){
			cout << i << "/" << t->GetEntries() << endl;
		}
		//    if (i % 50 == 0){
		//oniiiily copy branches where mass is in specified range   
		for (int j=0; j<Dsize; j++)
		{
			MVA[j]=-999;

			if(TMath::IsNaN(x3DDecayLengthSignificance[j])  || TMath::IsNaN(xyDCASignificanceDaugther1[j]) || TMath::IsNaN(xyDCASignificanceDaugther1[j])) continue;
			if(pT[j]<1) continue;
			if(pT[j]>=30) continue;
			if(centrality<0) continue;
			if(centrality>=180) continue;

            if(Dtrk1Chi2n[j]>=0.18) {
                continue;
            }
            if(Dtrk2Chi2n[j]>=0.18) {
                continue;
            }


			//centrality 0-10




			//if(Dtrk1Chi2n[j]<0.18 && Dtrk2Chi2n[j]<0.18 && Dtrk1Pt[j]>=1.0 && Dtrk2Pt[j]>=1.0)
			if(Dtrk1Pt[j]>=1.0 && Dtrk2Pt[j]>=1.0)
			{
				inputValues.clear();
				inputValues.push_back(log10(x3DDecayLengthSignificance[j]));//00
				inputValues.push_back(VtxProb[j]);//01
				inputValues.push_back(x3DPointingAngle[j]);//02
				inputValues.push_back(Dtrk1Pt[j]);//03
				inputValues.push_back(Dtrk2Pt[j]);//04
				inputValues.push_back(xyDCASignificanceDaugther1[j]);//04
				inputValues.push_back(xyDCASignificanceDaugther2[j]);//04
				inputValues.push_back(zDCASignificanceDaugther1[j]);//04
				inputValues.push_back(zDCASignificanceDaugther2[j]);


				if(centrality>=0 && centrality<20)
				{
					if (pT[j]>=1 && pT[j]< 2) MVA[j]=mva_cent_0_10_pt_1_2_.GetMvaValue_cent_0_10_pt_1_2_(inputValues);
                    if (pT[j]>=2 && pT[j]< 3)MVA[j]=mva_cent_0_10_pt_2_3_.GetMvaValue_cent_0_10_pt_2_3_(inputValues);
					if (pT[j]>=3 && pT[j]< 4)MVA[j]=mva_cent_0_10_pt_3_4_.GetMvaValue_cent_0_10_pt_3_4_(inputValues);
					if (pT[j]>=4 && pT[j]< 5)MVA[j]=mva_cent_0_10_pt_4_5_.GetMvaValue_cent_0_10_pt_4_5_(inputValues);
					if (pT[j]>=5 && pT[j]< 6)MVA[j]=mva_cent_0_10_pt_5_6_.GetMvaValue_cent_0_10_pt_5_6_(inputValues);
					if (pT[j]>=6 && pT[j]< 8)MVA[j]=mva_cent_0_10_pt_6_8_.GetMvaValue_cent_0_10_pt_6_8_(inputValues);
					if (pT[j]>=8 && pT[j]< 10)MVA[j]=mva_cent_0_10_pt_8_10_.GetMvaValue_cent_0_10_pt_8_10_(inputValues);
					if (pT[j]>=10 && pT[j]< 15)MVA[j]=mva_cent_0_10_pt_10_15_.GetMvaValue_cent_0_10_pt_10_15_(inputValues);
					if (pT[j]>=15 && pT[j]< 20)MVA[j]=mva_cent_0_10_pt_15_20_.GetMvaValue_cent_0_10_pt_15_20_(inputValues);
					if (pT[j]>=20 && pT[j]< 30)MVA[j]=mva_cent_0_10_pt_20_30_.GetMvaValue_cent_0_10_pt_20_30_(inputValues);
                    }//centrality 0-10
                    if(centrality>=20 && centrality<60)
				{
					if (pT[j]>=1 && pT[j]< 2)MVA[j]=mva_cent_10_30_pt_1_2_.GetMvaValue_cent_10_30_pt_1_2_(inputValues);
					if (pT[j]>=2 && pT[j]< 3)MVA[j]=mva_cent_10_30_pt_2_3_.GetMvaValue_cent_10_30_pt_2_3_(inputValues);
					if (pT[j]>=3 && pT[j]< 4)MVA[j]=mva_cent_10_30_pt_3_4_.GetMvaValue_cent_10_30_pt_3_4_(inputValues);
					if (pT[j]>=4 && pT[j]< 5)MVA[j]=mva_cent_10_30_pt_4_5_.GetMvaValue_cent_10_30_pt_4_5_(inputValues);
					if (pT[j]>=5 && pT[j]< 6)MVA[j]=mva_cent_10_30_pt_5_6_.GetMvaValue_cent_10_30_pt_5_6_(inputValues);
					if (pT[j]>=6 && pT[j]< 8)MVA[j]=mva_cent_10_30_pt_6_8_.GetMvaValue_cent_10_30_pt_6_8_(inputValues);
					if (pT[j]>=8 && pT[j]< 10)MVA[j]=mva_cent_10_30_pt_8_10_.GetMvaValue_cent_10_30_pt_8_10_(inputValues);
					if (pT[j]>=10 && pT[j]< 15)MVA[j]=mva_cent_10_30_pt_10_15_.GetMvaValue_cent_10_30_pt_10_15_(inputValues);
					if (pT[j]>=15 && pT[j]< 20)MVA[j]=mva_cent_10_30_pt_15_20_.GetMvaValue_cent_10_30_pt_15_20_(inputValues);
					if (pT[j]>=20 && pT[j]< 30)MVA[j]=mva_cent_10_30_pt_20_30_.GetMvaValue_cent_10_30_pt_20_30_(inputValues);
				}//centrality 10-30
				if(centrality>=60 && centrality<100)
				{
					if (pT[j]>=1 && pT[j]< 2)MVA[j]=mva_cent_30_50_pt_1_2_.GetMvaValue_cent_30_50_pt_1_2_(inputValues);
					if (pT[j]>=2 && pT[j]< 3)MVA[j]=mva_cent_30_50_pt_2_3_.GetMvaValue_cent_30_50_pt_2_3_(inputValues);
					if (pT[j]>=3 && pT[j]< 4)MVA[j]=mva_cent_30_50_pt_3_4_.GetMvaValue_cent_30_50_pt_3_4_(inputValues);
					if (pT[j]>=4 && pT[j]< 5)MVA[j]=mva_cent_30_50_pt_4_5_.GetMvaValue_cent_30_50_pt_4_5_(inputValues);
					if (pT[j]>=5 && pT[j]< 6)MVA[j]=mva_cent_30_50_pt_5_6_.GetMvaValue_cent_30_50_pt_5_6_(inputValues);
					if (pT[j]>=6 && pT[j]< 8)MVA[j]=mva_cent_30_50_pt_6_8_.GetMvaValue_cent_30_50_pt_6_8_(inputValues);
					if (pT[j]>=8 && pT[j]< 10)MVA[j]=mva_cent_30_50_pt_8_10_.GetMvaValue_cent_30_50_pt_8_10_(inputValues);
					if (pT[j]>=10 && pT[j]< 15)MVA[j]=mva_cent_30_50_pt_10_15_.GetMvaValue_cent_30_50_pt_10_15_(inputValues);
					if (pT[j]>=15 && pT[j]< 20)MVA[j]=mva_cent_30_50_pt_15_20_.GetMvaValue_cent_30_50_pt_15_20_(inputValues);
					if (pT[j]>=20 && pT[j]< 30)MVA[j]=mva_cent_30_50_pt_20_30_.GetMvaValue_cent_30_50_pt_20_30_(inputValues);
				}//centrality 30-50
				if(centrality>=100 && centrality<180)
				{
					if (pT[j]>=1 && pT[j]< 2)MVA[j]=mva_cent_50_90_pt_1_2_.GetMvaValue_cent_50_90_pt_1_2_(inputValues);
					if (pT[j]>=2 && pT[j]< 3)MVA[j]=mva_cent_50_90_pt_2_3_.GetMvaValue_cent_50_90_pt_2_3_(inputValues);
					if (pT[j]>=3 && pT[j]< 4)MVA[j]=mva_cent_50_90_pt_3_4_.GetMvaValue_cent_50_90_pt_3_4_(inputValues);
					if (pT[j]>=4 && pT[j]< 5)MVA[j]=mva_cent_50_90_pt_4_5_.GetMvaValue_cent_50_90_pt_4_5_(inputValues);
					if (pT[j]>=5 && pT[j]< 6)MVA[j]=mva_cent_50_90_pt_5_6_.GetMvaValue_cent_50_90_pt_5_6_(inputValues);
					if (pT[j]>=6 && pT[j]< 8)MVA[j]=mva_cent_50_90_pt_6_8_.GetMvaValue_cent_50_90_pt_6_8_(inputValues);
//					if (pT[j]>=8 && pT[j]< 30)MVA[j]=mva_cent_50_90_pt_8_30_.GetMvaValue_cent_50_90_pt_8_30_(inputValues);
				}//centrality 50-90
                }//chi2
		}//for dsize
			mvaTree->Fill();
	}//for entries

    outf->cd();
	mvaTree->Write(0,TObject::kOverwrite);
	outf->Close();
	return 0;
	} //program

	int main(int argc, char *argv[])
	{
		if(istest) argc=3;
		if(argc==3)
		{
			differentialBDT(argv[1], argv[2]);
		}
		else
		{
			std::cout << "Usage: mergeForest <input_collection> <output_file>" << std::endl;
			return 1;
		}
		return 0;


	}


