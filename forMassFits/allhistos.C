//GENERATES HISTOGRAMS (MASS DISTRIBUTIONS) FOR ALL PT AND CENTRALITY RANGES
//all histograms will be written into one file
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <vector>
#include <TMath.h>
#include <iostream>
# include <TF1.h>
# include <TH1F.h>
# include <TH1D.h>
# include <TCanvas.h>

using namespace std;
Bool_t istest = false;
int allhistos(TString inputfile="", TString outfile1="")
{
    //    TFile *inf = TFile::Open("/scratch/brown/awesole/bdtbg.root");


    if(istest)
    {
        //inputfile="/scratch/brown/awesole/output/Dntuple_d0ana_PbPb2018_04.root_327.root";
        //inputfile="/scratch/brown/awesole/bdtsignal.root";
        ////////////inputfile="/scratch/brown/chand140/ESE_output_Sept_18/analysis_data_Dntuple_with_multiplicity_Dtrk1/ROOT/outfile_ttree_562_563.root";
        //inputfile="/scratch/bell/awesole/backup_from_brown/ptbinsoutput/bdtvalue_d0ana_PbPb2018_04.root_599.root";
        //inputfile="/scratch/bell/awesole/backup_from_brown/turtle/bdt_d0ana_PbPb2018_02.root_704.root";
        inputfile="/scratch/bell/awesole/newbackgroundsample/histo_d0ana_PbPb2018_04.root_449.root";
        outfile1="/scratch/bell/awesole/testallhisto.root";
        cout << "TEST" << endl;
    }



    TFile* inf = TFile::Open(inputfile);
    TTree *t = (TTree*)inf->Get("mvaTree");
    t->SetMakeClass(1);
    TFile *outf = new TFile(outfile1,"recreate");

    int Dsize, centrality;
    float pT[10000], mass[10000], y[10000];
    double BDT_weight[10000];


    //aalist = centrlaity 0-10 pt 1-2
    //ablist = centrality 0-10 pt 2-3 
    std::vector<TH1D*> aalist;
    aalist.resize(25);
    std::vector<TH1D*> ablist;
    ablist.resize(25);
    std::vector<TH1D*> aclist;
    aclist.resize(25);
    std::vector<TH1D*> adlist;
    adlist.resize(25);
    std::vector<TH1D*> aelist;
    aelist.resize(25);
    std::vector<TH1D*> aflist;
    aflist.resize(25);
    std::vector<TH1D*> aglist;
    aglist.resize(25);
    std::vector<TH1D*> ahlist;
    ahlist.resize(25);
    std::vector<TH1D*> ailist;
    ailist.resize(25);
    std::vector<TH1D*> ajlist;
    ajlist.resize(25);

/*
    std::vector<TH1D*> balist;
    balist.resize(25);
    std::vector<TH1D*> bblist;
    bblist.resize(25);
    std::vector<TH1D*> bclist;
    bclist.resize(25);
    std::vector<TH1D*> bdlist;
    bdlist.resize(25);
    std::vector<TH1D*> belist;
    belist.resize(25);
    std::vector<TH1D*> bflist;
    bflist.resize(25);
    std::vector<TH1D*> bglist;
    bglist.resize(25);
    std::vector<TH1D*> bhlist;
    bhlist.resize(25);
    std::vector<TH1D*> bilist;
    bilist.resize(25);
    std::vector<TH1D*> bjlist;
    bjlist.resize(25);

    std::vector<TH1D*> calist;
    calist.resize(25);
    std::vector<TH1D*> cblist;
    cblist.resize(25);
    std::vector<TH1D*> cclist;
    cclist.resize(25);
    std::vector<TH1D*> cdlist;
    cdlist.resize(25);
    std::vector<TH1D*> celist;
    celist.resize(25);
    std::vector<TH1D*> cflist;
    cflist.resize(25);
    std::vector<TH1D*> cglist;
    cglist.resize(25);
    std::vector<TH1D*> chlist;
    chlist.resize(25);
    std::vector<TH1D*> cilist;
    cilist.resize(25);
    std::vector<TH1D*> cjlist;
    cjlist.resize(25);

    std::vector<TH1D*> dalist;
    dalist.resize(25);
    std::vector<TH1D*> dblist;
    dblist.resize(25);
    std::vector<TH1D*> dclist;
    dclist.resize(25);
    std::vector<TH1D*> ddlist;
    ddlist.resize(25);
    std::vector<TH1D*> delist;
    delist.resize(25);
    std::vector<TH1D*> dflist;
    dflist.resize(25);
    std::vector<TH1D*> dglist;
    dglist.resize(25);
*/

    for (int j=0; j<25; j++)
    {

        aalist[j]=new TH1D(Form("aahisto_%d",j),Form("BDTBin_%d",j), 52, 1.74, 2);
        ablist[j]=new TH1D(Form("abhisto_%d",j),Form("BDTBin_%d",j), 52, 1.74, 2);
        aclist[j]=new TH1D(Form("achisto_%d",j),Form("BDTBin_%d",j), 52, 1.74, 2);
        adlist[j]=new TH1D(Form("adhisto_%d",j),Form("BDTBin_%d",j), 52, 1.74, 2);
        aelist[j]=new TH1D(Form("aehisto_%d",j),Form("BDTBin_%d",j), 52, 1.74, 2);
        aflist[j]=new TH1D(Form("afhisto_%d",j),Form("BDTBin_%d",j), 52, 1.74, 2);
        aglist[j]=new TH1D(Form("aghisto_%d",j),Form("BDTBin_%d",j), 52, 1.74, 2);
        ahlist[j]=new TH1D(Form("ahhisto_%d",j),Form("BDTBin_%d",j), 52, 1.74, 2);
        ailist[j]=new TH1D(Form("aihisto_%d",j),Form("BDTBin_%d",j), 52, 1.74, 2);
        ajlist[j]=new TH1D(Form("ajhisto_%d",j),Form("BDTBin_%d",j), 52, 1.74, 2);

/*
        balist[j]=new TH1D(Form("bahisto_%d",j),Form("BDTBin_%d",j), 52, 1.74, 2);
        bblist[j]=new TH1D(Form("bbhisto_%d",j),Form("BDTBin_%d",j), 52, 1.74, 2);
        bclist[j]=new TH1D(Form("bchisto_%d",j),Form("BDTBin_%d",j), 52, 1.74, 2);
        bdlist[j]=new TH1D(Form("bdhisto_%d",j),Form("BDTBin_%d",j), 52, 1.74, 2);
        belist[j]=new TH1D(Form("behisto_%d",j),Form("BDTBin_%d",j), 52, 1.74, 2);
        bflist[j]=new TH1D(Form("bfhisto_%d",j),Form("BDTBin_%d",j), 52, 1.74, 2);
        bglist[j]=new TH1D(Form("bghisto_%d",j),Form("BDTBin_%d",j), 52, 1.74, 2);
        bhlist[j]=new TH1D(Form("bhhisto_%d",j),Form("BDTBin_%d",j), 52, 1.74, 2);
        bilist[j]=new TH1D(Form("bihisto_%d",j),Form("BDTBin_%d",j), 52, 1.74, 2);
        bjlist[j]=new TH1D(Form("bjhisto_%d",j),Form("BDTBin_%d",j), 52, 1.74, 2);


        calist[j]=new TH1D(Form("cahisto_%d",j),Form("BDTBin_%d",j), 52, 1.74, 2);
        cblist[j]=new TH1D(Form("cbhisto_%d",j),Form("BDTBin_%d",j), 52, 1.74, 2);
        cclist[j]=new TH1D(Form("cchisto_%d",j),Form("BDTBin_%d",j), 52, 1.74, 2);
        cdlist[j]=new TH1D(Form("cdhisto_%d",j),Form("BDTBin_%d",j), 52, 1.74, 2);
        celist[j]=new TH1D(Form("cehisto_%d",j),Form("BDTBin_%d",j), 52, 1.74, 2);
        cflist[j]=new TH1D(Form("cfhisto_%d",j),Form("BDTBin_%d",j), 52, 1.74, 2);
        cglist[j]=new TH1D(Form("cghisto_%d",j),Form("BDTBin_%d",j), 52, 1.74, 2);
        chlist[j]=new TH1D(Form("chhisto_%d",j),Form("BDTBin_%d",j), 52, 1.74, 2);
        cilist[j]=new TH1D(Form("cihisto_%d",j),Form("BDTBin_%d",j), 52, 1.74, 2);
        cjlist[j]=new TH1D(Form("cjhisto_%d",j),Form("BDTBin_%d",j), 52, 1.74, 2);


        dalist[j]=new TH1D(Form("dahisto_%d",j),Form("BDTBin_%d",j), 52, 1.74, 2);
        dblist[j]=new TH1D(Form("dbhisto_%d",j),Form("BDTBin_%d",j), 52, 1.74, 2);
        dclist[j]=new TH1D(Form("dchisto_%d",j),Form("BDTBin_%d",j), 52, 1.74, 2);
        ddlist[j]=new TH1D(Form("ddhisto_%d",j),Form("BDTBin_%d",j), 52, 1.74, 2);
        delist[j]=new TH1D(Form("dehisto_%d",j),Form("BDTBin_%d",j), 52, 1.74, 2);
        dflist[j]=new TH1D(Form("dfhisto_%d",j),Form("BDTBin_%d",j), 52, 1.74, 2);
        dglist[j]=new TH1D(Form("dghisto_%d",j),Form("BDTBin_%d",j), 52, 1.74, 2);
*/

    }
    //cout << "blist made" << endl;
    //for(int i=0;i<10000;i++)
    for(int i=0;i<t->GetEntries();i++)
    //for(int i=0;i<1000;i++)
    {
        t->SetBranchAddress("nDsize", &Dsize); 
        t->SetBranchAddress("npT", pT);
        t->SetBranchAddress("nmass", mass);
        t->SetBranchAddress("ncent", &centrality);
        t->SetBranchAddress("ny", y);
        t->SetBranchAddress("BDT_weight", BDT_weight);
        t->GetEntry(i);
        // cout << "Dsize= " << Dsize << endl;
        if(i%100000==0) cout <<i<<" / "<<t->GetEntries()<<endl;

        //   for(int j=0; j<Dsize; j++)
        for(int j=0; j<Dsize; j++)
        { 
            if(y[j]>-1 && y[j]<1)
            {
                  //cout << "y=" << y[j] << endl;
                if (centrality>=0 && centrality<20 )
                {
                    if (pT[j]>=1 && pT[j]<2){
                        if(BDT_weight[j]>-0.6) aalist[0]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.55) aalist[1]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.5) aalist[2]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.45) aalist[3]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.4) aalist[4]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.35) aalist[5]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.3) aalist[6]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.25) aalist[7]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.2) aalist[8]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.15) aalist[9]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.1) aalist[10]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.05) aalist[11]->Fill(mass[j]);
                        if(BDT_weight[j]>0.0) aalist[12]->Fill(mass[j]);
                        if(BDT_weight[j]>0.05) aalist[13]->Fill(mass[j]);
                        if(BDT_weight[j]>0.1) aalist[14]->Fill(mass[j]);
                        if(BDT_weight[j]>0.15) aalist[15]->Fill(mass[j]);
                        if(BDT_weight[j]>0.2) aalist[16]->Fill(mass[j]);
                        if(BDT_weight[j]>0.25) aalist[17]->Fill(mass[j]);
                        if(BDT_weight[j]>0.3) aalist[18]->Fill(mass[j]);
                        if(BDT_weight[j]>0.35) aalist[19]->Fill(mass[j]);
                        if(BDT_weight[j]>0.4) aalist[20]->Fill(mass[j]);
                        if(BDT_weight[j]>0.45) aalist[21]->Fill(mass[j]);
                        if(BDT_weight[j]>0.5) aalist[22]->Fill(mass[j]);
                        if(BDT_weight[j]>0.55) aalist[23]->Fill(mass[j]);
                        if(BDT_weight[j]>0.6) aalist[24]->Fill(mass[j]);
/*
                    }//if pT

                    if (pT[j]>=2 && pT[j]<3){
                        if(BDT_weight[j]>-0.6) ablist[0]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.55) ablist[1]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.5) ablist[2]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.45) ablist[3]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.4) ablist[4]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.35) ablist[5]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.3) ablist[6]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.25) ablist[7]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.2) ablist[8]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.15) ablist[9]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.1) ablist[10]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.05) ablist[11]->Fill(mass[j]);
                        if(BDT_weight[j]>0.0) ablist[12]->Fill(mass[j]);
                        if(BDT_weight[j]>0.05) ablist[13]->Fill(mass[j]);
                        if(BDT_weight[j]>0.1) ablist[14]->Fill(mass[j]);
                        if(BDT_weight[j]>0.15) ablist[15]->Fill(mass[j]);
                        if(BDT_weight[j]>0.2) ablist[16]->Fill(mass[j]);
                        if(BDT_weight[j]>0.25) ablist[17]->Fill(mass[j]);
                        if(BDT_weight[j]>0.3) ablist[18]->Fill(mass[j]);
                        if(BDT_weight[j]>0.35) ablist[19]->Fill(mass[j]);
                        if(BDT_weight[j]>0.4) ablist[20]->Fill(mass[j]);
                        if(BDT_weight[j]>0.45) ablist[21]->Fill(mass[j]);
                        if(BDT_weight[j]>0.5) ablist[22]->Fill(mass[j]);
                        if(BDT_weight[j]>0.55) ablist[23]->Fill(mass[j]);
                        if(BDT_weight[j]>0.6) ablist[24]->Fill(mass[j]);
                    }//if pT
                    if (pT[j]>=3 && pT[j]<4){
                        if(BDT_weight[j]>-0.6) aclist[0]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.55) aclist[1]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.5) aclist[2]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.45) aclist[3]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.4) aclist[4]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.35) aclist[5]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.3) aclist[6]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.25) aclist[7]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.2) aclist[8]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.15) aclist[9]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.1) aclist[10]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.05) aclist[11]->Fill(mass[j]);
                        if(BDT_weight[j]>0.0) aclist[12]->Fill(mass[j]);
                        if(BDT_weight[j]>0.05) aclist[13]->Fill(mass[j]);
                        if(BDT_weight[j]>0.1) aclist[14]->Fill(mass[j]);
                        if(BDT_weight[j]>0.15) aclist[15]->Fill(mass[j]);
                        if(BDT_weight[j]>0.2) aclist[16]->Fill(mass[j]);
                        if(BDT_weight[j]>0.25) aclist[17]->Fill(mass[j]);
                        if(BDT_weight[j]>0.3) aclist[18]->Fill(mass[j]);
                        if(BDT_weight[j]>0.35) aclist[19]->Fill(mass[j]);
                        if(BDT_weight[j]>0.4) aclist[20]->Fill(mass[j]);
                        if(BDT_weight[j]>0.45) aclist[21]->Fill(mass[j]);
                        if(BDT_weight[j]>0.5) aclist[22]->Fill(mass[j]);
                        if(BDT_weight[j]>0.55) aclist[23]->Fill(mass[j]);
                        if(BDT_weight[j]>0.6) aclist[24]->Fill(mass[j]);
                    }//if pT
                    if (pT[j]>=4 && pT[j]<5){
                        if(BDT_weight[j]>-0.6) adlist[0]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.55) adlist[1]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.5) adlist[2]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.45) adlist[3]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.4) adlist[4]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.35) adlist[5]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.3) adlist[6]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.25) adlist[7]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.2) adlist[8]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.15) adlist[9]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.1) adlist[10]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.05) adlist[11]->Fill(mass[j]);
                        if(BDT_weight[j]>0.0) adlist[12]->Fill(mass[j]);
                        if(BDT_weight[j]>0.05) adlist[13]->Fill(mass[j]);
                        if(BDT_weight[j]>0.1) adlist[14]->Fill(mass[j]);
                        if(BDT_weight[j]>0.15) adlist[15]->Fill(mass[j]);
                        if(BDT_weight[j]>0.2) adlist[16]->Fill(mass[j]);
                        if(BDT_weight[j]>0.25) adlist[17]->Fill(mass[j]);
                        if(BDT_weight[j]>0.3) adlist[18]->Fill(mass[j]);
                        if(BDT_weight[j]>0.35) adlist[19]->Fill(mass[j]);
                        if(BDT_weight[j]>0.4) adlist[20]->Fill(mass[j]);
                        if(BDT_weight[j]>0.45) adlist[21]->Fill(mass[j]);
                        if(BDT_weight[j]>0.5) adlist[22]->Fill(mass[j]);
                        if(BDT_weight[j]>0.55) adlist[23]->Fill(mass[j]);
                        if(BDT_weight[j]>0.6) adlist[24]->Fill(mass[j]);
                    }//if pT
                    if (pT[j]>=5 && pT[j]<6){
                        if(BDT_weight[j]>-0.6) aelist[0]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.55) aelist[1]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.5) aelist[2]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.45) aelist[3]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.4) aelist[4]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.35) aelist[5]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.3) aelist[6]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.25) aelist[7]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.2) aelist[8]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.15) aelist[9]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.1) aelist[10]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.05) aelist[11]->Fill(mass[j]);
                        if(BDT_weight[j]>0.0) aelist[12]->Fill(mass[j]);
                        if(BDT_weight[j]>0.05) aelist[13]->Fill(mass[j]);
                        if(BDT_weight[j]>0.1) aelist[14]->Fill(mass[j]);
                        if(BDT_weight[j]>0.15) aelist[15]->Fill(mass[j]);
                        if(BDT_weight[j]>0.2) aelist[16]->Fill(mass[j]);
                        if(BDT_weight[j]>0.25) aelist[17]->Fill(mass[j]);
                        if(BDT_weight[j]>0.3) aelist[18]->Fill(mass[j]);
                        if(BDT_weight[j]>0.35) aelist[19]->Fill(mass[j]);
                        if(BDT_weight[j]>0.4) aelist[20]->Fill(mass[j]);
                        if(BDT_weight[j]>0.45) aelist[21]->Fill(mass[j]);
                        if(BDT_weight[j]>0.5) aelist[22]->Fill(mass[j]);
                        if(BDT_weight[j]>0.55) aelist[23]->Fill(mass[j]);
                        if(BDT_weight[j]>0.6) aelist[24]->Fill(mass[j]);
                    }//if pT
                    if (pT[j]>=6 && pT[j]<8){
                        if(BDT_weight[j]>-0.6) aflist[0]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.55) aflist[1]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.5) aflist[2]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.45) aflist[3]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.4) aflist[4]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.35) aflist[5]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.3) aflist[6]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.25) aflist[7]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.2) aflist[8]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.15) aflist[9]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.1) aflist[10]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.05) aflist[11]->Fill(mass[j]);
                        if(BDT_weight[j]>0.0) aflist[12]->Fill(mass[j]);
                        if(BDT_weight[j]>0.05) aflist[13]->Fill(mass[j]);
                        if(BDT_weight[j]>0.1) aflist[14]->Fill(mass[j]);
                        if(BDT_weight[j]>0.15) aflist[15]->Fill(mass[j]);
                        if(BDT_weight[j]>0.2) aflist[16]->Fill(mass[j]);
                        if(BDT_weight[j]>0.25) aflist[17]->Fill(mass[j]);
                        if(BDT_weight[j]>0.3) aflist[18]->Fill(mass[j]);
                        if(BDT_weight[j]>0.35) aflist[19]->Fill(mass[j]);
                        if(BDT_weight[j]>0.4) aflist[20]->Fill(mass[j]);
                        if(BDT_weight[j]>0.45) aflist[21]->Fill(mass[j]);
                        if(BDT_weight[j]>0.5) aflist[22]->Fill(mass[j]);
                        if(BDT_weight[j]>0.55) aflist[23]->Fill(mass[j]);
                        if(BDT_weight[j]>0.6) aflist[24]->Fill(mass[j]);
                    }//if pT
                    if (pT[j]>=8 && pT[j]<10){
                        if(BDT_weight[j]>-0.6) aglist[0]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.55) aglist[1]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.5) aglist[2]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.45) aglist[3]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.4) aglist[4]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.35) aglist[5]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.3) aglist[6]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.25) aglist[7]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.2) aglist[8]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.15) aglist[9]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.1) aglist[10]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.05) aglist[11]->Fill(mass[j]);
                        if(BDT_weight[j]>0.0) aglist[12]->Fill(mass[j]);
                        if(BDT_weight[j]>0.05) aglist[13]->Fill(mass[j]);
                        if(BDT_weight[j]>0.1) aglist[14]->Fill(mass[j]);
                        if(BDT_weight[j]>0.15) aglist[15]->Fill(mass[j]);
                        if(BDT_weight[j]>0.2) aglist[16]->Fill(mass[j]);
                        if(BDT_weight[j]>0.25) aglist[17]->Fill(mass[j]);
                        if(BDT_weight[j]>0.3) aglist[18]->Fill(mass[j]);
                        if(BDT_weight[j]>0.35) aglist[19]->Fill(mass[j]);
                        if(BDT_weight[j]>0.4) aglist[20]->Fill(mass[j]);
                        if(BDT_weight[j]>0.45) aglist[21]->Fill(mass[j]);
                        if(BDT_weight[j]>0.5) aglist[22]->Fill(mass[j]);
                        if(BDT_weight[j]>0.55) aglist[23]->Fill(mass[j]);
                        if(BDT_weight[j]>0.6) aglist[24]->Fill(mass[j]);
                    }//if pT
                    if (pT[j]>=10 && pT[j]<15){
                        if(BDT_weight[j]>-0.6) ahlist[0]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.55) ahlist[1]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.5) ahlist[2]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.45) ahlist[3]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.4) ahlist[4]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.35) ahlist[5]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.3) ahlist[6]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.25) ahlist[7]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.2) ahlist[8]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.15) ahlist[9]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.1) ahlist[10]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.05) ahlist[11]->Fill(mass[j]);
                        if(BDT_weight[j]>0.0) ahlist[12]->Fill(mass[j]);
                        if(BDT_weight[j]>0.05) ahlist[13]->Fill(mass[j]);
                        if(BDT_weight[j]>0.1) ahlist[14]->Fill(mass[j]);
                        if(BDT_weight[j]>0.15) ahlist[15]->Fill(mass[j]);
                        if(BDT_weight[j]>0.2) ahlist[16]->Fill(mass[j]);
                        if(BDT_weight[j]>0.25) ahlist[17]->Fill(mass[j]);
                        if(BDT_weight[j]>0.3) ahlist[18]->Fill(mass[j]);
                        if(BDT_weight[j]>0.35) ahlist[19]->Fill(mass[j]);
                        if(BDT_weight[j]>0.4) ahlist[20]->Fill(mass[j]);
                        if(BDT_weight[j]>0.45) ahlist[21]->Fill(mass[j]);
                        if(BDT_weight[j]>0.5) ahlist[22]->Fill(mass[j]);
                        if(BDT_weight[j]>0.55) ahlist[23]->Fill(mass[j]);
                        if(BDT_weight[j]>0.6) ahlist[24]->Fill(mass[j]);
                    }//if pT
                    if (pT[j]>=15 && pT[j]<20){
                        if(BDT_weight[j]>-0.6) ailist[0]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.55) ailist[1]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.5) ailist[2]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.45) ailist[3]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.4) ailist[4]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.35) ailist[5]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.3) ailist[6]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.25) ailist[7]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.2) ailist[8]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.15) ailist[9]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.1) ailist[10]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.05) ailist[11]->Fill(mass[j]);
                        if(BDT_weight[j]>0.0) ailist[12]->Fill(mass[j]);
                        if(BDT_weight[j]>0.05) ailist[13]->Fill(mass[j]);
                        if(BDT_weight[j]>0.1) ailist[14]->Fill(mass[j]);
                        if(BDT_weight[j]>0.15) ailist[15]->Fill(mass[j]);
                        if(BDT_weight[j]>0.2) ailist[16]->Fill(mass[j]);
                        if(BDT_weight[j]>0.25) ailist[17]->Fill(mass[j]);
                        if(BDT_weight[j]>0.3) ailist[18]->Fill(mass[j]);
                        if(BDT_weight[j]>0.35) ailist[19]->Fill(mass[j]);
                        if(BDT_weight[j]>0.4) ailist[20]->Fill(mass[j]);
                        if(BDT_weight[j]>0.45) ailist[21]->Fill(mass[j]);
                        if(BDT_weight[j]>0.5) ailist[22]->Fill(mass[j]);
                        if(BDT_weight[j]>0.55) ailist[23]->Fill(mass[j]);
                        if(BDT_weight[j]>0.6) ailist[24]->Fill(mass[j]);
                    }//if pT
                    if (pT[j]>=20 && pT[j]<30){
                        if(BDT_weight[j]>-0.6) ajlist[0]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.55) ajlist[1]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.5) ajlist[2]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.45) ajlist[3]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.4) ajlist[4]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.35) ajlist[5]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.3) ajlist[6]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.25) ajlist[7]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.2) ajlist[8]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.15) ajlist[9]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.1) ajlist[10]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.05) ajlist[11]->Fill(mass[j]);
                        if(BDT_weight[j]>0.0) ajlist[12]->Fill(mass[j]);
                        if(BDT_weight[j]>0.05) ajlist[13]->Fill(mass[j]);
                        if(BDT_weight[j]>0.1) ajlist[14]->Fill(mass[j]);
                        if(BDT_weight[j]>0.15) ajlist[15]->Fill(mass[j]);
                        if(BDT_weight[j]>0.2) ajlist[16]->Fill(mass[j]);
                        if(BDT_weight[j]>0.25) ajlist[17]->Fill(mass[j]);
                        if(BDT_weight[j]>0.3) ajlist[18]->Fill(mass[j]);
                        if(BDT_weight[j]>0.35) ajlist[19]->Fill(mass[j]);
                        if(BDT_weight[j]>0.4) ajlist[20]->Fill(mass[j]);
                        if(BDT_weight[j]>0.45) ajlist[21]->Fill(mass[j]);
                        if(BDT_weight[j]>0.5) ajlist[22]->Fill(mass[j]);
                        if(BDT_weight[j]>0.55) ajlist[23]->Fill(mass[j]);
                        if(BDT_weight[j]>0.6) ajlist[24]->Fill(mass[j]);
                    }//if pT

                }//if centrality
                if (centrality>=20 && centrality<60 )
                {

                    if (pT[j]>1 && pT[j]<2){
                        if(BDT_weight[j]>-0.6) balist[0]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.55) balist[1]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.5) balist[2]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.45) balist[3]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.4) balist[4]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.35) balist[5]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.3) balist[6]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.25) balist[7]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.2) balist[8]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.15) balist[9]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.1) balist[10]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.05) balist[11]->Fill(mass[j]);
                        if(BDT_weight[j]>0.0) balist[12]->Fill(mass[j]);
                        if(BDT_weight[j]>0.05) balist[13]->Fill(mass[j]);
                        if(BDT_weight[j]>0.1) balist[14]->Fill(mass[j]);
                        if(BDT_weight[j]>0.15) balist[15]->Fill(mass[j]);
                        if(BDT_weight[j]>0.2) balist[16]->Fill(mass[j]);
                        if(BDT_weight[j]>0.25) balist[17]->Fill(mass[j]);
                        if(BDT_weight[j]>0.3) balist[18]->Fill(mass[j]);
                        if(BDT_weight[j]>0.35) balist[19]->Fill(mass[j]);
                        if(BDT_weight[j]>0.4) balist[20]->Fill(mass[j]);
                        if(BDT_weight[j]>0.45) balist[21]->Fill(mass[j]);
                        if(BDT_weight[j]>0.5) balist[22]->Fill(mass[j]);
                        if(BDT_weight[j]>0.55) balist[23]->Fill(mass[j]);
                        if(BDT_weight[j]>0.6) balist[24]->Fill(mass[j]);
                    }//if pT

                    if (pT[j]>2 && pT[j]<3){
                        if(BDT_weight[j]>-0.6) bblist[0]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.55) bblist[1]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.5) bblist[2]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.45) bblist[3]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.4) bblist[4]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.35) bblist[5]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.3) bblist[6]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.25) bblist[7]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.2) bblist[8]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.15) bblist[9]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.1) bblist[10]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.05) bblist[11]->Fill(mass[j]);
                        if(BDT_weight[j]>0.0) bblist[12]->Fill(mass[j]);
                        if(BDT_weight[j]>0.05) bblist[13]->Fill(mass[j]);
                        if(BDT_weight[j]>0.1) bblist[14]->Fill(mass[j]);
                        if(BDT_weight[j]>0.15) bblist[15]->Fill(mass[j]);
                        if(BDT_weight[j]>0.2) bblist[16]->Fill(mass[j]);
                        if(BDT_weight[j]>0.25) bblist[17]->Fill(mass[j]);
                        if(BDT_weight[j]>0.3) bblist[18]->Fill(mass[j]);
                        if(BDT_weight[j]>0.35) bblist[19]->Fill(mass[j]);
                        if(BDT_weight[j]>0.4) bblist[20]->Fill(mass[j]);
                        if(BDT_weight[j]>0.45) bblist[21]->Fill(mass[j]);
                        if(BDT_weight[j]>0.5) bblist[22]->Fill(mass[j]);
                        if(BDT_weight[j]>0.55) bblist[23]->Fill(mass[j]);
                        if(BDT_weight[j]>0.6) bblist[24]->Fill(mass[j]);
                    }//if pT
                    if (pT[j]>3 && pT[j]<4){
                        if(BDT_weight[j]>-0.6) bclist[0]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.55) bclist[1]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.5) bclist[2]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.45) bclist[3]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.4) bclist[4]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.35) bclist[5]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.3) bclist[6]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.25) bclist[7]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.2) bclist[8]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.15) bclist[9]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.1) bclist[10]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.05) bclist[11]->Fill(mass[j]);
                        if(BDT_weight[j]>0.0) bclist[12]->Fill(mass[j]);
                        if(BDT_weight[j]>0.05) bclist[13]->Fill(mass[j]);
                        if(BDT_weight[j]>0.1) bclist[14]->Fill(mass[j]);
                        if(BDT_weight[j]>0.15) bclist[15]->Fill(mass[j]);
                        if(BDT_weight[j]>0.2) bclist[16]->Fill(mass[j]);
                        if(BDT_weight[j]>0.25) bclist[17]->Fill(mass[j]);
                        if(BDT_weight[j]>0.3) bclist[18]->Fill(mass[j]);
                        if(BDT_weight[j]>0.35) bclist[19]->Fill(mass[j]);
                        if(BDT_weight[j]>0.4) bclist[20]->Fill(mass[j]);
                        if(BDT_weight[j]>0.45) bclist[21]->Fill(mass[j]);
                        if(BDT_weight[j]>0.5) bclist[22]->Fill(mass[j]);
                        if(BDT_weight[j]>0.55) bclist[23]->Fill(mass[j]);
                        if(BDT_weight[j]>0.6) bclist[24]->Fill(mass[j]);
                    }//if pT
                    if (pT[j]>4 && pT[j]<5){
                        if(BDT_weight[j]>-0.6) bdlist[0]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.55) bdlist[1]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.5) bdlist[2]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.45) bdlist[3]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.4) bdlist[4]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.35) bdlist[5]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.3) bdlist[6]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.25) bdlist[7]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.2) bdlist[8]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.15) bdlist[9]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.1) bdlist[10]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.05) bdlist[11]->Fill(mass[j]);
                        if(BDT_weight[j]>0.0) bdlist[12]->Fill(mass[j]);
                        if(BDT_weight[j]>0.05) bdlist[13]->Fill(mass[j]);
                        if(BDT_weight[j]>0.1) bdlist[14]->Fill(mass[j]);
                        if(BDT_weight[j]>0.15) bdlist[15]->Fill(mass[j]);
                        if(BDT_weight[j]>0.2) bdlist[16]->Fill(mass[j]);
                        if(BDT_weight[j]>0.25) bdlist[17]->Fill(mass[j]);
                        if(BDT_weight[j]>0.3) bdlist[18]->Fill(mass[j]);
                        if(BDT_weight[j]>0.35) bdlist[19]->Fill(mass[j]);
                        if(BDT_weight[j]>0.4) bdlist[20]->Fill(mass[j]);
                        if(BDT_weight[j]>0.45) bdlist[21]->Fill(mass[j]);
                        if(BDT_weight[j]>0.5) bdlist[22]->Fill(mass[j]);
                        if(BDT_weight[j]>0.55) bdlist[23]->Fill(mass[j]);
                        if(BDT_weight[j]>0.6) bdlist[24]->Fill(mass[j]);
                    }//if pT
                    if (pT[j]>5 && pT[j]<6){
                        if(BDT_weight[j]>-0.6) belist[0]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.55) belist[1]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.5) belist[2]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.45) belist[3]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.4) belist[4]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.35) belist[5]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.3) belist[6]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.25) belist[7]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.2) belist[8]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.15) belist[9]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.1) belist[10]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.05) belist[11]->Fill(mass[j]);
                        if(BDT_weight[j]>0.0) belist[12]->Fill(mass[j]);
                        if(BDT_weight[j]>0.05) belist[13]->Fill(mass[j]);
                        if(BDT_weight[j]>0.1) belist[14]->Fill(mass[j]);
                        if(BDT_weight[j]>0.15) belist[15]->Fill(mass[j]);
                        if(BDT_weight[j]>0.2) belist[16]->Fill(mass[j]);
                        if(BDT_weight[j]>0.25) belist[17]->Fill(mass[j]);
                        if(BDT_weight[j]>0.3) belist[18]->Fill(mass[j]);
                        if(BDT_weight[j]>0.35) belist[19]->Fill(mass[j]);
                        if(BDT_weight[j]>0.4) belist[20]->Fill(mass[j]);
                        if(BDT_weight[j]>0.45) belist[21]->Fill(mass[j]);
                        if(BDT_weight[j]>0.5) belist[22]->Fill(mass[j]);
                        if(BDT_weight[j]>0.55) belist[23]->Fill(mass[j]);
                        if(BDT_weight[j]>0.6) belist[24]->Fill(mass[j]);
                    }//if pT
                    if (pT[j]>6 && pT[j]<8){
                        if(BDT_weight[j]>-0.6) bflist[0]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.55) bflist[1]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.5) bflist[2]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.45) bflist[3]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.4) bflist[4]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.35) bflist[5]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.3) bflist[6]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.25) bflist[7]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.2) bflist[8]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.15) bflist[9]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.1) bflist[10]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.05) bflist[11]->Fill(mass[j]);
                        if(BDT_weight[j]>0.0) bflist[12]->Fill(mass[j]);
                        if(BDT_weight[j]>0.05) bflist[13]->Fill(mass[j]);
                        if(BDT_weight[j]>0.1) bflist[14]->Fill(mass[j]);
                        if(BDT_weight[j]>0.15) bflist[15]->Fill(mass[j]);
                        if(BDT_weight[j]>0.2) bflist[16]->Fill(mass[j]);
                        if(BDT_weight[j]>0.25) bflist[17]->Fill(mass[j]);
                        if(BDT_weight[j]>0.3) bflist[18]->Fill(mass[j]);
                        if(BDT_weight[j]>0.35) bflist[19]->Fill(mass[j]);
                        if(BDT_weight[j]>0.4) bflist[20]->Fill(mass[j]);
                        if(BDT_weight[j]>0.45) bflist[21]->Fill(mass[j]);
                        if(BDT_weight[j]>0.5) bflist[22]->Fill(mass[j]);
                        if(BDT_weight[j]>0.55) bflist[23]->Fill(mass[j]);
                        if(BDT_weight[j]>0.6) bflist[24]->Fill(mass[j]);
                    }//if pT
                    if (pT[j]>8 && pT[j]<10){
                        if(BDT_weight[j]>-0.6) bglist[0]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.55) bglist[1]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.5) bglist[2]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.45) bglist[3]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.4) bglist[4]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.35) bglist[5]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.3) bglist[6]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.25) bglist[7]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.2) bglist[8]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.15) bglist[9]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.1) bglist[10]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.05) bglist[11]->Fill(mass[j]);
                        if(BDT_weight[j]>0.0) bglist[12]->Fill(mass[j]);
                        if(BDT_weight[j]>0.05) bglist[13]->Fill(mass[j]);
                        if(BDT_weight[j]>0.1) bglist[14]->Fill(mass[j]);
                        if(BDT_weight[j]>0.15) bglist[15]->Fill(mass[j]);
                        if(BDT_weight[j]>0.2) bglist[16]->Fill(mass[j]);
                        if(BDT_weight[j]>0.25) bglist[17]->Fill(mass[j]);
                        if(BDT_weight[j]>0.3) bglist[18]->Fill(mass[j]);
                        if(BDT_weight[j]>0.35) bglist[19]->Fill(mass[j]);
                        if(BDT_weight[j]>0.4) bglist[20]->Fill(mass[j]);
                        if(BDT_weight[j]>0.45) bglist[21]->Fill(mass[j]);
                        if(BDT_weight[j]>0.5) bglist[22]->Fill(mass[j]);
                        if(BDT_weight[j]>0.55) bglist[23]->Fill(mass[j]);
                        if(BDT_weight[j]>0.6) bglist[24]->Fill(mass[j]);
                    }//if pT
                    if (pT[j]>10 && pT[j]<15){
                        if(BDT_weight[j]>-0.6) bhlist[0]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.55) bhlist[1]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.5) bhlist[2]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.45) bhlist[3]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.4) bhlist[4]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.35) bhlist[5]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.3) bhlist[6]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.25) bhlist[7]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.2) bhlist[8]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.15) bhlist[9]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.1) bhlist[10]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.05) bhlist[11]->Fill(mass[j]);
                        if(BDT_weight[j]>0.0) bhlist[12]->Fill(mass[j]);
                        if(BDT_weight[j]>0.05) bhlist[13]->Fill(mass[j]);
                        if(BDT_weight[j]>0.1) bhlist[14]->Fill(mass[j]);
                        if(BDT_weight[j]>0.15) bhlist[15]->Fill(mass[j]);
                        if(BDT_weight[j]>0.2) bhlist[16]->Fill(mass[j]);
                        if(BDT_weight[j]>0.25) bhlist[17]->Fill(mass[j]);
                        if(BDT_weight[j]>0.3) bhlist[18]->Fill(mass[j]);
                        if(BDT_weight[j]>0.35) bhlist[19]->Fill(mass[j]);
                        if(BDT_weight[j]>0.4) bhlist[20]->Fill(mass[j]);
                        if(BDT_weight[j]>0.45) bhlist[21]->Fill(mass[j]);
                        if(BDT_weight[j]>0.5) bhlist[22]->Fill(mass[j]);
                        if(BDT_weight[j]>0.55) bhlist[23]->Fill(mass[j]);
                        if(BDT_weight[j]>0.6) bhlist[24]->Fill(mass[j]);
                    }//if pT
                    if (pT[j]>15 && pT[j]<20){
                        if(BDT_weight[j]>-0.6) bilist[0]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.55) bilist[1]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.5) bilist[2]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.45) bilist[3]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.4) bilist[4]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.35) bilist[5]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.3) bilist[6]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.25) bilist[7]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.2) bilist[8]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.15) bilist[9]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.1) bilist[10]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.05) bilist[11]->Fill(mass[j]);
                        if(BDT_weight[j]>0.0) bilist[12]->Fill(mass[j]);
                        if(BDT_weight[j]>0.05) bilist[13]->Fill(mass[j]);
                        if(BDT_weight[j]>0.1) bilist[14]->Fill(mass[j]);
                        if(BDT_weight[j]>0.15) bilist[15]->Fill(mass[j]);
                        if(BDT_weight[j]>0.2) bilist[16]->Fill(mass[j]);
                        if(BDT_weight[j]>0.25) bilist[17]->Fill(mass[j]);
                        if(BDT_weight[j]>0.3) bilist[18]->Fill(mass[j]);
                        if(BDT_weight[j]>0.35) bilist[19]->Fill(mass[j]);
                        if(BDT_weight[j]>0.4) bilist[20]->Fill(mass[j]);
                        if(BDT_weight[j]>0.45) bilist[21]->Fill(mass[j]);
                        if(BDT_weight[j]>0.5) bilist[22]->Fill(mass[j]);
                        if(BDT_weight[j]>0.55) bilist[23]->Fill(mass[j]);
                        if(BDT_weight[j]>0.6) bilist[24]->Fill(mass[j]);
                    }//if pT
                    if (pT[j]>20 && pT[j]<30){
                        if(BDT_weight[j]>-0.6) bjlist[0]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.55) bjlist[1]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.5) bjlist[2]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.45) bjlist[3]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.4) bjlist[4]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.35) bjlist[5]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.3) bjlist[6]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.25) bjlist[7]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.2) bjlist[8]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.15) bjlist[9]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.1) bjlist[10]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.05) bjlist[11]->Fill(mass[j]);
                        if(BDT_weight[j]>0.0) bjlist[12]->Fill(mass[j]);
                        if(BDT_weight[j]>0.05) bjlist[13]->Fill(mass[j]);
                        if(BDT_weight[j]>0.1) bjlist[14]->Fill(mass[j]);
                        if(BDT_weight[j]>0.15) bjlist[15]->Fill(mass[j]);
                        if(BDT_weight[j]>0.2) bjlist[16]->Fill(mass[j]);
                        if(BDT_weight[j]>0.25) bjlist[17]->Fill(mass[j]);
                        if(BDT_weight[j]>0.3) bjlist[18]->Fill(mass[j]);
                        if(BDT_weight[j]>0.35) bjlist[19]->Fill(mass[j]);
                        if(BDT_weight[j]>0.4) bjlist[20]->Fill(mass[j]);
                        if(BDT_weight[j]>0.45) bjlist[21]->Fill(mass[j]);
                        if(BDT_weight[j]>0.5) bjlist[22]->Fill(mass[j]);
                        if(BDT_weight[j]>0.55) bjlist[23]->Fill(mass[j]);
                        if(BDT_weight[j]>0.6) bjlist[24]->Fill(mass[j]);
                    }//if pT
                }//if centrality
                if (centrality>=60 && centrality<100 )
                {

                    if (pT[j]>1 && pT[j]<2){
                        if(BDT_weight[j]>-0.6) calist[0]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.55) calist[1]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.5) calist[2]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.45) calist[3]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.4) calist[4]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.35) calist[5]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.3) calist[6]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.25) calist[7]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.2) calist[8]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.15) calist[9]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.1) calist[10]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.05) calist[11]->Fill(mass[j]);
                        if(BDT_weight[j]>0.0) calist[12]->Fill(mass[j]);
                        if(BDT_weight[j]>0.05) calist[13]->Fill(mass[j]);
                        if(BDT_weight[j]>0.1) calist[14]->Fill(mass[j]);
                        if(BDT_weight[j]>0.15) calist[15]->Fill(mass[j]);
                        if(BDT_weight[j]>0.2) calist[16]->Fill(mass[j]);
                        if(BDT_weight[j]>0.25) calist[17]->Fill(mass[j]);
                        if(BDT_weight[j]>0.3) calist[18]->Fill(mass[j]);
                        if(BDT_weight[j]>0.35) calist[19]->Fill(mass[j]);
                        if(BDT_weight[j]>0.4) calist[20]->Fill(mass[j]);
                        if(BDT_weight[j]>0.45) calist[21]->Fill(mass[j]);
                        if(BDT_weight[j]>0.5) calist[22]->Fill(mass[j]);
                        if(BDT_weight[j]>0.55) calist[23]->Fill(mass[j]);
                        if(BDT_weight[j]>0.6) calist[24]->Fill(mass[j]);
                    }//if pT

                    if (pT[j]>2 && pT[j]<3){
                        if(BDT_weight[j]>-0.6) cblist[0]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.55) cblist[1]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.5) cblist[2]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.45) cblist[3]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.4) cblist[4]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.35) cblist[5]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.3) cblist[6]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.25) cblist[7]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.2) cblist[8]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.15) cblist[9]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.1) cblist[10]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.05) cblist[11]->Fill(mass[j]);
                        if(BDT_weight[j]>0.0) cblist[12]->Fill(mass[j]);
                        if(BDT_weight[j]>0.05) cblist[13]->Fill(mass[j]);
                        if(BDT_weight[j]>0.1) cblist[14]->Fill(mass[j]);
                        if(BDT_weight[j]>0.15) cblist[15]->Fill(mass[j]);
                        if(BDT_weight[j]>0.2) cblist[16]->Fill(mass[j]);
                        if(BDT_weight[j]>0.25) cblist[17]->Fill(mass[j]);
                        if(BDT_weight[j]>0.3) cblist[18]->Fill(mass[j]);
                        if(BDT_weight[j]>0.35) cblist[19]->Fill(mass[j]);
                        if(BDT_weight[j]>0.4) cblist[20]->Fill(mass[j]);
                        if(BDT_weight[j]>0.45) cblist[21]->Fill(mass[j]);
                        if(BDT_weight[j]>0.5) cblist[22]->Fill(mass[j]);
                        if(BDT_weight[j]>0.55) cblist[23]->Fill(mass[j]);
                        if(BDT_weight[j]>0.6) cblist[24]->Fill(mass[j]);
                    }//if pT
                    if (pT[j]>3 && pT[j]<4){
                        if(BDT_weight[j]>-0.6) cclist[0]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.55) cclist[1]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.5) cclist[2]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.45) cclist[3]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.4) cclist[4]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.35) cclist[5]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.3) cclist[6]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.25) cclist[7]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.2) cclist[8]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.15) cclist[9]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.1) cclist[10]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.05) cclist[11]->Fill(mass[j]);
                        if(BDT_weight[j]>0.0) cclist[12]->Fill(mass[j]);
                        if(BDT_weight[j]>0.05) cclist[13]->Fill(mass[j]);
                        if(BDT_weight[j]>0.1) cclist[14]->Fill(mass[j]);
                        if(BDT_weight[j]>0.15) cclist[15]->Fill(mass[j]);
                        if(BDT_weight[j]>0.2) cclist[16]->Fill(mass[j]);
                        if(BDT_weight[j]>0.25) cclist[17]->Fill(mass[j]);
                        if(BDT_weight[j]>0.3) cclist[18]->Fill(mass[j]);
                        if(BDT_weight[j]>0.35) cclist[19]->Fill(mass[j]);
                        if(BDT_weight[j]>0.4) cclist[20]->Fill(mass[j]);
                        if(BDT_weight[j]>0.45) cclist[21]->Fill(mass[j]);
                        if(BDT_weight[j]>0.5) cclist[22]->Fill(mass[j]);
                        if(BDT_weight[j]>0.55) cclist[23]->Fill(mass[j]);
                        if(BDT_weight[j]>0.6) cclist[24]->Fill(mass[j]);
                    }//if pT
                    if (pT[j]>4 && pT[j]<5){
                        if(BDT_weight[j]>-0.6) cdlist[0]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.55) cdlist[1]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.5) cdlist[2]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.45) cdlist[3]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.4) cdlist[4]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.35) cdlist[5]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.3) cdlist[6]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.25) cdlist[7]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.2) cdlist[8]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.15) cdlist[9]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.1) cdlist[10]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.05) cdlist[11]->Fill(mass[j]);
                        if(BDT_weight[j]>0.0) cdlist[12]->Fill(mass[j]);
                        if(BDT_weight[j]>0.05) cdlist[13]->Fill(mass[j]);
                        if(BDT_weight[j]>0.1) cdlist[14]->Fill(mass[j]);
                        if(BDT_weight[j]>0.15) cdlist[15]->Fill(mass[j]);
                        if(BDT_weight[j]>0.2) cdlist[16]->Fill(mass[j]);
                        if(BDT_weight[j]>0.25) cdlist[17]->Fill(mass[j]);
                        if(BDT_weight[j]>0.3) cdlist[18]->Fill(mass[j]);
                        if(BDT_weight[j]>0.35) cdlist[19]->Fill(mass[j]);
                        if(BDT_weight[j]>0.4) cdlist[20]->Fill(mass[j]);
                        if(BDT_weight[j]>0.45) cdlist[21]->Fill(mass[j]);
                        if(BDT_weight[j]>0.5) cdlist[22]->Fill(mass[j]);
                        if(BDT_weight[j]>0.55) cdlist[23]->Fill(mass[j]);
                        if(BDT_weight[j]>0.6) cdlist[24]->Fill(mass[j]);
                    }//if pT
                    if (pT[j]>5 && pT[j]<6){
                        if(BDT_weight[j]>-0.6) celist[0]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.55) celist[1]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.5) celist[2]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.45) celist[3]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.4) celist[4]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.35) celist[5]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.3) celist[6]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.25) celist[7]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.2) celist[8]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.15) celist[9]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.1) celist[10]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.05) celist[11]->Fill(mass[j]);
                        if(BDT_weight[j]>0.0) celist[12]->Fill(mass[j]);
                        if(BDT_weight[j]>0.05) celist[13]->Fill(mass[j]);
                        if(BDT_weight[j]>0.1) celist[14]->Fill(mass[j]);
                        if(BDT_weight[j]>0.15) celist[15]->Fill(mass[j]);
                        if(BDT_weight[j]>0.2) celist[16]->Fill(mass[j]);
                        if(BDT_weight[j]>0.25) celist[17]->Fill(mass[j]);
                        if(BDT_weight[j]>0.3) celist[18]->Fill(mass[j]);
                        if(BDT_weight[j]>0.35) celist[19]->Fill(mass[j]);
                        if(BDT_weight[j]>0.4) celist[20]->Fill(mass[j]);
                        if(BDT_weight[j]>0.45) celist[21]->Fill(mass[j]);
                        if(BDT_weight[j]>0.5) celist[22]->Fill(mass[j]);
                        if(BDT_weight[j]>0.55) celist[23]->Fill(mass[j]);
                        if(BDT_weight[j]>0.6) celist[24]->Fill(mass[j]);
                    }//if pT
                    if (pT[j]>6 && pT[j]<8){
                        if(BDT_weight[j]>-0.6) cflist[0]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.55) cflist[1]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.5) cflist[2]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.45) cflist[3]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.4) cflist[4]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.35) cflist[5]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.3) cflist[6]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.25) cflist[7]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.2) cflist[8]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.15) cflist[9]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.1) cflist[10]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.05) cflist[11]->Fill(mass[j]);
                        if(BDT_weight[j]>0.0) cflist[12]->Fill(mass[j]);
                        if(BDT_weight[j]>0.05) cflist[13]->Fill(mass[j]);
                        if(BDT_weight[j]>0.1) cflist[14]->Fill(mass[j]);
                        if(BDT_weight[j]>0.15) cflist[15]->Fill(mass[j]);
                        if(BDT_weight[j]>0.2) cflist[16]->Fill(mass[j]);
                        if(BDT_weight[j]>0.25) cflist[17]->Fill(mass[j]);
                        if(BDT_weight[j]>0.3) cflist[18]->Fill(mass[j]);
                        if(BDT_weight[j]>0.35) cflist[19]->Fill(mass[j]);
                        if(BDT_weight[j]>0.4) cflist[20]->Fill(mass[j]);
                        if(BDT_weight[j]>0.45) cflist[21]->Fill(mass[j]);
                        if(BDT_weight[j]>0.5) cflist[22]->Fill(mass[j]);
                        if(BDT_weight[j]>0.55) cflist[23]->Fill(mass[j]);
                        if(BDT_weight[j]>0.6) cflist[24]->Fill(mass[j]);
                    }//if pT
                    if (pT[j]>8 && pT[j]<10){
                        if(BDT_weight[j]>-0.6) cglist[0]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.55) cglist[1]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.5) cglist[2]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.45) cglist[3]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.4) cglist[4]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.35) cglist[5]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.3) cglist[6]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.25) cglist[7]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.2) cglist[8]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.15) cglist[9]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.1) cglist[10]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.05) cglist[11]->Fill(mass[j]);
                        if(BDT_weight[j]>0.0) cglist[12]->Fill(mass[j]);
                        if(BDT_weight[j]>0.05) cglist[13]->Fill(mass[j]);
                        if(BDT_weight[j]>0.1) cglist[14]->Fill(mass[j]);
                        if(BDT_weight[j]>0.15) cglist[15]->Fill(mass[j]);
                        if(BDT_weight[j]>0.2) cglist[16]->Fill(mass[j]);
                        if(BDT_weight[j]>0.25) cglist[17]->Fill(mass[j]);
                        if(BDT_weight[j]>0.3) cglist[18]->Fill(mass[j]);
                        if(BDT_weight[j]>0.35) cglist[19]->Fill(mass[j]);
                        if(BDT_weight[j]>0.4) cglist[20]->Fill(mass[j]);
                        if(BDT_weight[j]>0.45) cglist[21]->Fill(mass[j]);
                        if(BDT_weight[j]>0.5) cglist[22]->Fill(mass[j]);
                        if(BDT_weight[j]>0.55) cglist[23]->Fill(mass[j]);
                        if(BDT_weight[j]>0.6) cglist[24]->Fill(mass[j]);
                    }//if pT
                    if (pT[j]>10 && pT[j]<15){
                        if(BDT_weight[j]>-0.6) chlist[0]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.55) chlist[1]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.5) chlist[2]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.45) chlist[3]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.4) chlist[4]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.35) chlist[5]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.3) chlist[6]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.25) chlist[7]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.2) chlist[8]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.15) chlist[9]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.1) chlist[10]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.05) chlist[11]->Fill(mass[j]);
                        if(BDT_weight[j]>0.0) chlist[12]->Fill(mass[j]);
                        if(BDT_weight[j]>0.05) chlist[13]->Fill(mass[j]);
                        if(BDT_weight[j]>0.1) chlist[14]->Fill(mass[j]);
                        if(BDT_weight[j]>0.15) chlist[15]->Fill(mass[j]);
                        if(BDT_weight[j]>0.2) chlist[16]->Fill(mass[j]);
                        if(BDT_weight[j]>0.25) chlist[17]->Fill(mass[j]);
                        if(BDT_weight[j]>0.3) chlist[18]->Fill(mass[j]);
                        if(BDT_weight[j]>0.35) chlist[19]->Fill(mass[j]);
                        if(BDT_weight[j]>0.4) chlist[20]->Fill(mass[j]);
                        if(BDT_weight[j]>0.45) chlist[21]->Fill(mass[j]);
                        if(BDT_weight[j]>0.5) chlist[22]->Fill(mass[j]);
                        if(BDT_weight[j]>0.55) chlist[23]->Fill(mass[j]);
                        if(BDT_weight[j]>0.6) chlist[24]->Fill(mass[j]);
                    }//if pT
                    if (pT[j]>15 && pT[j]<20){
                        if(BDT_weight[j]>-0.6) cilist[0]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.55) cilist[1]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.5) cilist[2]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.45) cilist[3]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.4) cilist[4]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.35) cilist[5]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.3) cilist[6]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.25) cilist[7]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.2) cilist[8]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.15) cilist[9]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.1) cilist[10]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.05) cilist[11]->Fill(mass[j]);
                        if(BDT_weight[j]>0.0) cilist[12]->Fill(mass[j]);
                        if(BDT_weight[j]>0.05) cilist[13]->Fill(mass[j]);
                        if(BDT_weight[j]>0.1) cilist[14]->Fill(mass[j]);
                        if(BDT_weight[j]>0.15) cilist[15]->Fill(mass[j]);
                        if(BDT_weight[j]>0.2) cilist[16]->Fill(mass[j]);
                        if(BDT_weight[j]>0.25) cilist[17]->Fill(mass[j]);
                        if(BDT_weight[j]>0.3) cilist[18]->Fill(mass[j]);
                        if(BDT_weight[j]>0.35) cilist[19]->Fill(mass[j]);
                        if(BDT_weight[j]>0.4) cilist[20]->Fill(mass[j]);
                        if(BDT_weight[j]>0.45) cilist[21]->Fill(mass[j]);
                        if(BDT_weight[j]>0.5) cilist[22]->Fill(mass[j]);
                        if(BDT_weight[j]>0.55) cilist[23]->Fill(mass[j]);
                        if(BDT_weight[j]>0.6) cilist[24]->Fill(mass[j]);
                    }//if pT
                    if (pT[j]>20 && pT[j]<30){
                        if(BDT_weight[j]>-0.6) cjlist[0]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.55) cjlist[1]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.5) cjlist[2]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.45) cjlist[3]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.4) cjlist[4]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.35) cjlist[5]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.3) cjlist[6]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.25) cjlist[7]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.2) cjlist[8]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.15) cjlist[9]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.1) cjlist[10]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.05) cjlist[11]->Fill(mass[j]);
                        if(BDT_weight[j]>0.0) cjlist[12]->Fill(mass[j]);
                        if(BDT_weight[j]>0.05) cjlist[13]->Fill(mass[j]);
                        if(BDT_weight[j]>0.1) cjlist[14]->Fill(mass[j]);
                        if(BDT_weight[j]>0.15) cjlist[15]->Fill(mass[j]);
                        if(BDT_weight[j]>0.2) cjlist[16]->Fill(mass[j]);
                        if(BDT_weight[j]>0.25) cjlist[17]->Fill(mass[j]);
                        if(BDT_weight[j]>0.3) cjlist[18]->Fill(mass[j]);
                        if(BDT_weight[j]>0.35) cjlist[19]->Fill(mass[j]);
                        if(BDT_weight[j]>0.4) cjlist[20]->Fill(mass[j]);
                        if(BDT_weight[j]>0.45) cjlist[21]->Fill(mass[j]);
                        if(BDT_weight[j]>0.5) cjlist[22]->Fill(mass[j]);
                        if(BDT_weight[j]>0.55) cjlist[23]->Fill(mass[j]);
                        if(BDT_weight[j]>0.6) cjlist[24]->Fill(mass[j]);
                    }//if pT
                }//if centrality
                if (centrality>=100 && centrality<180 )
                {
                    if (pT[j]>1 && pT[j]<2){
                        if(BDT_weight[j]>-0.6) dalist[0]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.55) dalist[1]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.5) dalist[2]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.45) dalist[3]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.4) dalist[4]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.35) dalist[5]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.3) dalist[6]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.25) dalist[7]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.2) dalist[8]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.15) dalist[9]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.1) dalist[10]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.05) dalist[11]->Fill(mass[j]);
                        if(BDT_weight[j]>0.0) dalist[12]->Fill(mass[j]);
                        if(BDT_weight[j]>0.05) dalist[13]->Fill(mass[j]);
                        if(BDT_weight[j]>0.1) dalist[14]->Fill(mass[j]);
                        if(BDT_weight[j]>0.15) dalist[15]->Fill(mass[j]);
                        if(BDT_weight[j]>0.2) dalist[16]->Fill(mass[j]);
                        if(BDT_weight[j]>0.25) dalist[17]->Fill(mass[j]);
                        if(BDT_weight[j]>0.3) dalist[18]->Fill(mass[j]);
                        if(BDT_weight[j]>0.35) dalist[19]->Fill(mass[j]);
                        if(BDT_weight[j]>0.4) dalist[20]->Fill(mass[j]);
                        if(BDT_weight[j]>0.45) dalist[21]->Fill(mass[j]);
                        if(BDT_weight[j]>0.5) dalist[22]->Fill(mass[j]);
                        if(BDT_weight[j]>0.55) dalist[23]->Fill(mass[j]);
                        if(BDT_weight[j]>0.6) dalist[24]->Fill(mass[j]);
                    }//if pT

                    if (pT[j]>2 && pT[j]<3){
                        if(BDT_weight[j]>-0.6) dblist[0]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.55) dblist[1]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.5) dblist[2]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.45) dblist[3]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.4) dblist[4]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.35) dblist[5]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.3) dblist[6]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.25) dblist[7]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.2) dblist[8]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.15) dblist[9]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.1) dblist[10]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.05) dblist[11]->Fill(mass[j]);
                        if(BDT_weight[j]>0.0) dblist[12]->Fill(mass[j]);
                        if(BDT_weight[j]>0.05) dblist[13]->Fill(mass[j]);
                        if(BDT_weight[j]>0.1) dblist[14]->Fill(mass[j]);
                        if(BDT_weight[j]>0.15) dblist[15]->Fill(mass[j]);
                        if(BDT_weight[j]>0.2) dblist[16]->Fill(mass[j]);
                        if(BDT_weight[j]>0.25) dblist[17]->Fill(mass[j]);
                        if(BDT_weight[j]>0.3) dblist[18]->Fill(mass[j]);
                        if(BDT_weight[j]>0.35) dblist[19]->Fill(mass[j]);
                        if(BDT_weight[j]>0.4) dblist[20]->Fill(mass[j]);
                        if(BDT_weight[j]>0.45) dblist[21]->Fill(mass[j]);
                        if(BDT_weight[j]>0.5) dblist[22]->Fill(mass[j]);
                        if(BDT_weight[j]>0.55) dblist[23]->Fill(mass[j]);
                        if(BDT_weight[j]>0.6) dblist[24]->Fill(mass[j]);
                    }//if pT
                    if (pT[j]>3 && pT[j]<4){
                        if(BDT_weight[j]>-0.6) dclist[0]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.55) dclist[1]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.5) dclist[2]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.45) dclist[3]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.4) dclist[4]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.35) dclist[5]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.3) dclist[6]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.25) dclist[7]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.2) dclist[8]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.15) dclist[9]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.1) dclist[10]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.05) dclist[11]->Fill(mass[j]);
                        if(BDT_weight[j]>0.0) dclist[12]->Fill(mass[j]);
                        if(BDT_weight[j]>0.05) dclist[13]->Fill(mass[j]);
                        if(BDT_weight[j]>0.1) dclist[14]->Fill(mass[j]);
                        if(BDT_weight[j]>0.15) dclist[15]->Fill(mass[j]);
                        if(BDT_weight[j]>0.2) dclist[16]->Fill(mass[j]);
                        if(BDT_weight[j]>0.25) dclist[17]->Fill(mass[j]);
                        if(BDT_weight[j]>0.3) dclist[18]->Fill(mass[j]);
                        if(BDT_weight[j]>0.35) dclist[19]->Fill(mass[j]);
                        if(BDT_weight[j]>0.4) dclist[20]->Fill(mass[j]);
                        if(BDT_weight[j]>0.45) dclist[21]->Fill(mass[j]);
                        if(BDT_weight[j]>0.5) dclist[22]->Fill(mass[j]);
                        if(BDT_weight[j]>0.55) dclist[23]->Fill(mass[j]);
                        if(BDT_weight[j]>0.6) dclist[24]->Fill(mass[j]);
                    }//if pT
                    if (pT[j]>4 && pT[j]<5){
                        if(BDT_weight[j]>-0.6) ddlist[0]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.55) ddlist[1]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.5) ddlist[2]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.45) ddlist[3]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.4) ddlist[4]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.35) ddlist[5]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.3) ddlist[6]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.25) ddlist[7]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.2) ddlist[8]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.15) ddlist[9]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.1) ddlist[10]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.05) ddlist[11]->Fill(mass[j]);
                        if(BDT_weight[j]>0.0) ddlist[12]->Fill(mass[j]);
                        if(BDT_weight[j]>0.05) ddlist[13]->Fill(mass[j]);
                        if(BDT_weight[j]>0.1) ddlist[14]->Fill(mass[j]);
                        if(BDT_weight[j]>0.15) ddlist[15]->Fill(mass[j]);
                        if(BDT_weight[j]>0.2) ddlist[16]->Fill(mass[j]);
                        if(BDT_weight[j]>0.25) ddlist[17]->Fill(mass[j]);
                        if(BDT_weight[j]>0.3) ddlist[18]->Fill(mass[j]);
                        if(BDT_weight[j]>0.35) ddlist[19]->Fill(mass[j]);
                        if(BDT_weight[j]>0.4) ddlist[20]->Fill(mass[j]);
                        if(BDT_weight[j]>0.45) ddlist[21]->Fill(mass[j]);
                        if(BDT_weight[j]>0.5) ddlist[22]->Fill(mass[j]);
                        if(BDT_weight[j]>0.55) ddlist[23]->Fill(mass[j]);
                        if(BDT_weight[j]>0.6) ddlist[24]->Fill(mass[j]);
                    }//if pT
                    if (pT[j]>5 && pT[j]<6){
                        if(BDT_weight[j]>-0.6) delist[0]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.55) delist[1]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.5) delist[2]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.45) delist[3]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.4) delist[4]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.35) delist[5]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.3) delist[6]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.25) delist[7]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.2) delist[8]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.15) delist[9]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.1) delist[10]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.05) delist[11]->Fill(mass[j]);
                        if(BDT_weight[j]>0.0) delist[12]->Fill(mass[j]);
                        if(BDT_weight[j]>0.05) delist[13]->Fill(mass[j]);
                        if(BDT_weight[j]>0.1) delist[14]->Fill(mass[j]);
                        if(BDT_weight[j]>0.15) delist[15]->Fill(mass[j]);
                        if(BDT_weight[j]>0.2) delist[16]->Fill(mass[j]);
                        if(BDT_weight[j]>0.25) delist[17]->Fill(mass[j]);
                        if(BDT_weight[j]>0.3) delist[18]->Fill(mass[j]);
                        if(BDT_weight[j]>0.35) delist[19]->Fill(mass[j]);
                        if(BDT_weight[j]>0.4) delist[20]->Fill(mass[j]);
                        if(BDT_weight[j]>0.45) delist[21]->Fill(mass[j]);
                        if(BDT_weight[j]>0.5) delist[22]->Fill(mass[j]);
                        if(BDT_weight[j]>0.55) delist[23]->Fill(mass[j]);
                        if(BDT_weight[j]>0.6) delist[24]->Fill(mass[j]);
                    }//if pT
                    if (pT[j]>6 && pT[j]<8){
                        if(BDT_weight[j]>-0.6) dflist[0]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.55) dflist[1]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.5) dflist[2]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.45) dflist[3]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.4) dflist[4]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.35) dflist[5]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.3) dflist[6]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.25) dflist[7]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.2) dflist[8]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.15) dflist[9]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.1) dflist[10]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.05) dflist[11]->Fill(mass[j]);
                        if(BDT_weight[j]>0.0) dflist[12]->Fill(mass[j]);
                        if(BDT_weight[j]>0.05) dflist[13]->Fill(mass[j]);
                        if(BDT_weight[j]>0.1) dflist[14]->Fill(mass[j]);
                        if(BDT_weight[j]>0.15) dflist[15]->Fill(mass[j]);
                        if(BDT_weight[j]>0.2) dflist[16]->Fill(mass[j]);
                        if(BDT_weight[j]>0.25) dflist[17]->Fill(mass[j]);
                        if(BDT_weight[j]>0.3) dflist[18]->Fill(mass[j]);
                        if(BDT_weight[j]>0.35) dflist[19]->Fill(mass[j]);
                        if(BDT_weight[j]>0.4) dflist[20]->Fill(mass[j]);
                        if(BDT_weight[j]>0.45) dflist[21]->Fill(mass[j]);
                        if(BDT_weight[j]>0.5) dflist[22]->Fill(mass[j]);
                        if(BDT_weight[j]>0.55) dflist[23]->Fill(mass[j]);
                        if(BDT_weight[j]>0.6) dflist[24]->Fill(mass[j]);
*/ 
                   }//if pT
                }//if centrality
            }//if y
        }//for nDsize

    }//for entries 
    outf->cd();
    outf->Write();
    outf->Close();
    return 0;
}//function 


int main(int argc, char *argv[])
{
    if(istest) argc=3;
    if(argc==3)
    {
        allhistos(argv[1], argv[2]);
    }
    else
    {
        std::cout << "Usage: mergeForest <input_collection> <output_file>" << std::endl;
        return 1;
    }
    return 0;
}
