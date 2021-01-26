/* std::vectorの使い方を覚えよう
.push_back(a)    一番後ろにaを追加
.insert(i,a)    i番目にaを追加
.pop_back()  一番後ろを削除
.erase(i)   i番目を削除
.erase(i,j) i番目からj番目までを削除
.begin()    一番前を指す
.back() 一番後ろを指す
std::sort(v.begin(),v.end()) //小さい順にソートされる(vector<int>などの時)
Nbins = 500;
Xmax = 311.92;
Xmin = -88.08;
Width= 0.8;
*/

#include <TCanvas.h>
#include <TLegend.h>    //凡例
#include <TFile.h>
#include <TH1D.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <iostream> //std::
#include <sstream>
#include <fstream> //ファイル入出力
#include <TRandom.h>

void MakeNoiseDist(TH1D *noise, int noise_lambda = 0) {
    if (noise_lambda == 0) {
        for (int i=0; i<39990; ++i) {
            double x =gRandom->Gaus(0.0,5);
            noise->Fill(x);
            }
    } else {
        for (int i=0; i<39990; ++i) {
            double x = gRandom->Gaus(0.0,5);
            int npe = gRandom->Poisson(0.1*noise_lambda);
            for (int j=0; j<npe; ++j) {
                double noisecharge = gRandom->Gaus(15,5);
                (noisecharge > 0) ? x += noisecharge : x -= noisecharge;
            }
            noise->Fill(x);
        }
    }
}

void MakeOnePEDist(TH1D *OnePE, int ampgain = 10) {
    double AmpGain = 0.1 * ampgain ;
    for (int i=0; i<39990; ++i){
        double charge = gRandom->Gaus(20,10);
        if (charge < 0) charge *= -1;
        OnePE->Fill(AmpGain * charge);
    }
}

void MakeSignalDist(TH1D *signal, TH1D* noise, TH1D *OnePE,double lamb) {
    std::vector<int> EachEntries(10,0);
    for (int i=0; i<39990; ++i) {
        double charge = 0;
        int npe = gRandom->Poisson(lamb);
        for (int j = 0; j < npe ;++j) charge += OnePE->GetRandom();
        charge += noise->GetRandom();
        signal->Fill(charge);
        EachEntries[npe] += 1;
    }
    for (int i=0; i<10; ++i) {
        std::cout << std::right << std::setw(4) << i <<"PE" ;
    }
    std::cout <<std::endl;
    for (const auto& E : EachEntries) {
        std::cout << std::right << std::setw(6) << E ;
    }
    std::cout << std::endl;
    std::cout << "calc from lambda and N_entries" <<std::endl;
    double N_k = 39990 / exp(lamb);
    for (int i=1;i<10;++i) {
        std::cout << N_k << " ";
        N_k *= double(lamb / i);
    }
    std::cout << std::endl;
}

void Simulation() {
    TCanvas *can = new TCanvas("can","can", 750, 500);
    int nbins = 2500;
    double xmin = -88.08;
    double xmax = 311.92;
    double N_entries = 39990;
    double lambda = 1.0;
    int AmpGain = 15;
    int noise_lambda = 1;
    TH1D *noise= new TH1D("noise", "", nbins, xmin, xmax);  
    MakeNoiseDist(noise,noise_lambda);
    TH1D *OnePE= new TH1D("OnePE", "", nbins, xmin, xmax); 
    OnePE->SetLineColor(6);
    MakeOnePEDist(OnePE,AmpGain);
    TH1D *signal= new TH1D("signal", "", nbins, xmin, xmax); 
    signal->SetLineColor(2);
    MakeSignalDist(signal,noise,OnePE,1.0);
    gPad->SetLogy(1);
    signal->GetXaxis()->SetRangeUser(xmin,xmax);
    signal->GetYaxis()->SetRangeUser(0.8,1000);
    signal->Draw("same");
    noise->Draw("same");
    OnePE->Scale(lambda/exp(lambda));

    OnePE->Draw("same");

    std::vector<std::string> filenames = {};
    std::string nameH = "Sim/";
    std::string nameB = "normal_gain";
    std::string nameC = std::to_string(AmpGain);
    std::string nameD = "_Noise";
    std::string nameE = std::to_string(noise_lambda);
    std::string nameF = ".root";
    std::string fname = nameH + nameB + nameC + nameD + nameE + nameF;
    std::cout << fname << std::endl; 
    TFile *fout = new TFile(fname.c_str(), "recreate");
    std::cout << "Noise Mean = " << noise->GetMean() << " ± " << noise->GetMeanError() <<std::endl ;
    std::cout << "SignalMean = " << signal->GetMean() << " ± " << signal->GetMeanError() <<std::endl ;
    std::cout << "OnePE Mean = " << OnePE->GetMean() << " ± " << OnePE->GetMeanError() <<std::endl ;
    signal->Write();
    noise->Write();
    OnePE->Write();
    fout->Write();
    fout->Close();
    }