//12/10　GitHub

/*
lamb_check.Cでチェック
root [0] .L okumura1202.C+
root [1] auto anas = test()　391行目が実行
これをやると、まず ana で古田さんのデータで先ほどまでと同様に解析します。
次に、ana2 という instance には MakeMCSpectra を実行することで、
ana で使用した pedestal と抽出された 1 p.e. 分布を使って、
電荷分布を乱数で作り直し、これを再度解析しています。
TCanvas が 2 枚出てきますが、もし 1 枚目の最後の赤分布（ana2 の MC に入れた 1 p.e. 分布）と 
2 枚目の最後の赤分布（解析結果）が同じ形状になっていれば、入力の真の分布に戻せているということになります。
今のところ、まだ駄目そうです。
*/


/*
if (!hoge) {}で、hogeが存在しないというエラーの時の処理をかく

n個のbinに等分割する
TH1F h_rebin* = hist->Rebin(n);
TH1F h_rebin* = hist->Rebin(n,"h_rebin");//第二引数は新しく作るhistoの名前
好きな幅で分割する
double x[n] = {};
TH1F h_rebin* = hist->Rebin(n,"h_rebin",xbins);

std::vector::back
vector内の末尾要素への参照ポインターを渡す
std::vector<int> v = {3, 1, 4};
int& x = v.back();　//ポインターなのでint&
std::cout << x << std::endl;
//4
*/
#define READPATH "sample_data_PMT/turn_6_angle_5_CH1.root"
// #define READPATH "PMT_HVChange/turn12_CH1.root"

#include <TCanvas.h>
#include <TLegend.h>
#include <TFile.h>
#include <TH1D.h>
#include <TRandom.h>

#include <iostream>
#include <map>
#include <memory>
#include <vector>
#include <utility>
#include <fstream> // include for Savedata_Lamb()


double integrateUpToMinusOne(std::shared_ptr<TH1D> h) {
  // used for e.g. \Sigma_{-\inf} ^ {-1} n_all(i)
  int minbin = 1;
  int bin0 = h->FindBin(0);
  int maxbin = bin0 - 1;
  return h->Integral(minbin, maxbin);
}

class SinglePEAnalyzer {
 private:
  int kMaxPE = 3;  // Max num of p.e.s to be taken into account
  std::shared_ptr<TCanvas> fCanvas;
  std::shared_ptr<TH1D> fSignalH1;
  std::shared_ptr<TH1D> fNoiseH1;
  std::shared_ptr<TH1D> fSuminus0;
  double alphaH1;
  double fN0;
  double fN1;
  double fLambda;
  double fLambda_err_fromN;
  double cLambda; //Lambda calc by <Q_all>(Mean of fSignal) / <Q_1>(Mean of fNPEDist.back()[1] )
  double cLambda_err;

  std::vector<std::vector<std::shared_ptr<TH1D>>> fNPEDist;

  // std::map<std::string, double> fConfiguration;

 public:
  enum EMakeNPE {kMakeNPE_Takahashi2018, kMakeNPE_Fast};

  SinglePEAnalyzer();
  void ReadFile(const std::string& file_name, const std::string& signal_name,
                const std::string& noise_name, int rebin = 1);
  void MakeMCSpectra(std::shared_ptr<TH1D> noise, std::shared_ptr<TH1D> singlePE, double lambda, int nevents, bool ignorePedestalConvolution = false);
  double GetLambda() const {return fLambda;}
  double GetLambdaErr() const {return fLambda_err_fromN ;}
  double GetLambda_C() const {return cLambda;}
  double GetLambda_CErr() const {return cLambda_err;}
  double GetAlpha() const {return alphaH1;}
  std::shared_ptr<TH1D> GetSignalH1() const { return fSignalH1; }
  std::shared_ptr<TH1D> GetNoiseH1() const { return fNoiseH1; }
  std::shared_ptr<TH1D> GetSuminus0() const { return fSuminus0; }
  std::shared_ptr<TH1D> GetNPEDist(int n) const {
    try {
      return fNPEDist.back()[n];
    } catch (...) {
    }
    return nullptr;
  }
  // void SetConiguration(std::string name, double par) {fConfiguration[name] =
  // par;}
  void Analyze();
  void MakefSuminus0();
  void SetkMaxPE(int Iteration_MaxPE) {
    kMaxPE = Iteration_MaxPE;
  }
  void Iterate();
  void OnePEiterate(double N);
  void MakeNPEdist(EMakeNPE mode = kMakeNPE_Takahashi2018);
  void MakeNPEdistTakahashi2018();
  void MakeNPEdistFast();
  void CompareLambda_E_to_C();
  void Savedata_Lamb();
  double GetAmpGain();
};

SinglePEAnalyzer::SinglePEAnalyzer() {
  //Analyzeの先頭に移動
}

void SinglePEAnalyzer::ReadFile(const std::string& file_name,
                                const std::string& signal_name,
                                const std::string& noise_name,
                                int rebin) {
  // Get signal (LED on) and noise (LED off) spectra
  // Rebin(rebin) handles the roll of PMTbinsrough.C
  auto directory = gDirectory;
  TFile f(file_name.c_str());
  fSignalH1.reset((TH1D*)f.Get(signal_name.c_str())->Clone("h1_signal"));
  fNoiseH1.reset((TH1D*)f.Get(noise_name.c_str())->Clone("h1_noise"));
  fSignalH1->SetDirectory(directory);
  fNoiseH1->SetDirectory(directory);

  if(rebin > 1){
    fSignalH1->Rebin(rebin);
    fNoiseH1->Rebin(rebin);
  }
  f.Close();
}

void SinglePEAnalyzer::MakeMCSpectra(std::shared_ptr<TH1D> noise, std::shared_ptr<TH1D> singlePE, double lambda, int nevents, bool ignorePedestalConvolution) {
  fNPEDist.clear();
  fNoiseH1 = noise;

  fSignalH1 = std::make_shared<TH1D>("h1_signal", "", singlePE->GetNbinsX(), singlePE->GetXaxis()->GetXmin(), singlePE->GetXaxis()->GetXmax());

  for(int i = 0; i < nevents; ++i){
    int npe = gRandom->Poisson(lambda);
    double charge = 0;

    for(int j = 0; j < npe; ++j){
      charge += singlePE->GetRandom();
    }
    
    if (ignorePedestalConvolution) {
      if (npe == 0) {
        fSignalH1->Fill(fNoiseH1->GetRandom()); // use the pedestal dist
      } else {
        fSignalH1->Fill(charge);
      }
    } else {
      fSignalH1->Fill(charge + fNoiseH1->GetRandom()); // convlute pedestal
    }
  }
}

void SinglePEAnalyzer::Analyze() {
  // main method
  std::cout << "fSignalMean is " << fSignalH1->GetMean() << "±" << fSignalH1->GetMeanError() <<std::endl;

  fCanvas = std::make_shared<TCanvas>("can", "can", 2400, 800);
  fCanvas->Divide(4, 2, 1e-10, 1e-10);
  if (!fNoiseH1 || !fSignalH1) {
    std::cerr << "ERROR: Use SinglePEAnalizer::ReadFile first to fill "
                 "fSignalH1 and fNoiseH1"
              << std::endl;
    return;
  }
  auto xmin = fSignalH1->GetXaxis()->GetXmin();
  auto xmax = fSignalH1->GetXaxis()->GetXmax();
  auto binwidth = fSignalH1->GetXaxis()->GetBinWidth(1);

  // TODO: the charge offset (= binwidth / 2 here) should be automatically
  // calculated from the pedestal histogram
  double aveQ = fSignalH1->GetMean() + binwidth / 2;

  auto Nonall = fSignalH1->Integral(1, -1);
  auto Noffall = fNoiseH1->Integral(1, -1);
  auto nbins = (xmax - xmin) / binwidth;
  fSignalH1->SetLineColor(5);

  // suminus0 means distribution of N[k>0]
  MakefSuminus0();

  double Noff0 = fNoiseH1->Integral();  // by definition above Eq. (1)
  fN0 = 1.00 * Noff0 * alphaH1;  // (3)の下らへんに書いてるやつ
  double avePE0 = -log(
      fN0 / Nonall);  // <k> in Eq. (4), initial estimate of lamda
                      // NallとN0から<k>が推定可能なので、N0の推定のみで決まる
  fN1 = fN0 * avePE0;
  fLambda = avePE0;
  std::cout << "alpha: " << alphaH1 << std::endl;
  std::cout << "fN0: " << fN0 << std::endl;
  fLambda_err_fromN = sqrt((1-alphaH1)/fN0);
  std::cout << "lambda calc from -log(fN0/Nonall) : " << fLambda << " ± " << fLambda_err_fromN << std::endl;
  
  double Q1_8_err; //p12のdelta<Q1>_(8)
  for (int i = 1; i < fSignalH1->FindBin(0); ++i) {
    Q1_8_err += pow((xmin + binwidth * i) * avePE0 + aveQ, 2);
  }
  for (int i = fSignalH1->FindBin(0); i <= nbins; ++i) {
    Q1_8_err += pow((xmin + binwidth * i), 2) * fSignalH1->GetBinContent(i);
  }
  Q1_8_err /= pow(Nonall * avePE0, 2);

  fCanvas->cd(1);
  fSuminus0->SetFillColor(15);//Gray
  fSuminus0->SetTitle("fSignalH1(yellow), fNoiseH1(blue),fSuminus0(gray)");
  fSuminus0->GetYaxis()->SetRangeUser(0.8,1000);
  fSuminus0->Draw();
  fNoiseH1->Draw("same");
  fSignalH1->Draw("same");
  gPad->SetLogy(1);
  
  //iterationはここから
  for (int i = 0; i < 6; ++i) {
    Iterate();
    fCanvas->cd(i + 2);
    gPad->SetLogy(1);
    fNPEDist.back()[0]->GetYaxis()->SetRangeUser(0.8,1000);
    fNPEDist.back()[0]->Draw();//1PEを作る為の、fSuminus0から前のiterationの2PE,3PE分布を除いたものが表示
    for (int npe = 1; npe <= kMaxPE; ++npe) fNPEDist.back()[npe]->Draw("same");//NPE分布
    fSignalH1->SetLineColor(6);
    fSignalH1->Draw("same");
  }
  cLambda = fSignalH1->GetMean() / fNPEDist.back()[1]->GetMean() ;
  cLambda_err =sqrt( pow(fSignalH1->GetMeanError()/fNPEDist.back()[1]->GetMean(),2) + pow(fSignalH1->GetMean()*fNPEDist.back()[1]->GetStdDev()/(fNPEDist.back()[1]->GetMean()*fNPEDist.back()[1]->GetMean()),2) /fNPEDist.back()[1]->Integral(1,-1) );

  CompareLambda_E_to_C(); //lambdaを出す二つの方法で違いが出ないかチェック
  std::cout << "ANS lambda from entry  :" << fLambda <<" ± " << fLambda_err_fromN  << std::endl;
  std::cout << "ANS lambda from charge :" << cLambda <<" ± " << cLambda_err << std::endl;
  std::cout << "fSignalMean(平均電荷量) is " << fSignalH1->GetMean() << "±" << fSignalH1->GetMeanError() <<std::endl;

  std::cout << "1PE分布に関する情報" << std::endl; 
  std::cout << "1PE dist Entry:" << fNPEDist.back()[1]->Integral(1,-1) << std::endl;
  std::cout << "1PE dist Mean :" << fNPEDist.back()[1]->GetMean() <<" ± " << fNPEDist.back()[1]->GetMeanError() << std::endl;
  std::cout << "1PE dist StdDv:" << fNPEDist.back()[1]->GetStdDev() <<" ± " << fNPEDist.back()[1]->GetStdDevError() << std::endl;

}

void SinglePEAnalyzer::MakefSuminus0() {
  fSuminus0.reset((TH1D*)fSignalH1->Clone("suminus0"));
  auto integrate_negative_charge = [](std::shared_ptr<TH1D> h) {
    // \Sigma_{-\inf} ^ {-1} n_all(i)
    return h->Integral(1, h->FindBin(0) - 1);
  };
  double numerator = integrate_negative_charge(fSignalH1);
  double denominator = integrate_negative_charge(fNoiseH1);
  alphaH1 = numerator / denominator;  // See Eq. (3)
  fSuminus0->Add(fNoiseH1.get(),
                 -1 * alphaH1);  // n_{k > 0}(i) = n_all(i) - n_0(i) 

  // 0回りを整理する。fNoiseH1->GetStdDev()よりシグマの範囲を同じびんに放り込む
  // 1.5σの範囲にするかどうかは議論の余地がある
  double fNoiseH1_sigma = fNoiseH1->GetStdDev();
  printf("fNoiseH1 StdDev is %f\n", fNoiseH1_sigma);
  printf("this is %d bins (0bin is %d)\n", fNoiseH1->FindBin(fNoiseH1_sigma),
         fNoiseH1->FindBin(0));

  double together_bin_charge = 1;
  int together_bin = 1; // = fSuminus->FindBin(together_bin_charge);
  double together_bin_content = 0;
  double together_bin_numerator = 0;    // bunshi
  double together_bin_denominator = 0;  // bunbo

  for (int i = 1; i <= fNoiseH1->GetNbinsX(); ++i) {
    if (fNoiseH1->GetBinCenter(i) <=
        0 * fNoiseH1_sigma) {  //ここの倍率は議論の余地がある
      together_bin_content += fSuminus0->GetBinContent(i);
      together_bin_denominator += fSuminus0->GetBinContent(i);
      together_bin_numerator +=
          fSuminus0->GetBinContent(i) * fNoiseH1->GetBinCenter(i);
      fSuminus0->SetBinContent(i, 0);
    }
  }
  together_bin_charge = (together_bin_numerator / together_bin_denominator);
  together_bin = fSuminus0->FindBin(together_bin_charge);

  printf("%f / %f so, toukeigosa_matomete_ireru_bin is %d \n",
         together_bin_numerator, together_bin_denominator, together_bin);
  printf("content is %f \n", together_bin_content);
  fSuminus0->SetBinContent(together_bin, together_bin_content);
  for (int i = 1; i <= fNoiseH1->GetNbinsX() / 2; ++i) {
    if (fSuminus0->GetBinContent(i) < 0) fSuminus0->SetBinContent(i, 0);
  }  // without it, contents of some bins can become negative number
  fSuminus0->SetFillColor(1);

}

void SinglePEAnalyzer::Iterate() {
  auto n = fNPEDist.size();
  auto nbins = fSignalH1->GetNbinsX();
  auto xmin = fSignalH1->GetXaxis()->GetXmin();
  auto xmax = fSignalH1->GetXaxis()->GetXmax();
  std::cout << "Now iteration is " << n << std::endl;

  fNPEDist.push_back(std::vector<std::shared_ptr<TH1D>>());
  for (int i = 0; i <= kMaxPE; ++i) {
    fNPEDist.back().push_back(
        std::make_shared<TH1D>(Form("NPEDist_%dpe_%luth", i, n),
                               ";Charge;Entries/bin", nbins, xmin, xmax));
    fNPEDist.back().back()->SetLineColor(i + 1);
    if (i == 1) {
      fNPEDist.back().back()->SetFillColor(2);
    }
  }

  auto h0 = fNPEDist.back()[0];
  h0->Add(fSuminus0.get());

  if (n != 0) {
    for (int npe=2; npe<=kMaxPE; ++npe){
      auto prevNPE = fNPEDist[n - 1][npe].get();
      h0->Add(prevNPE, -1);
    }
    for (int i = 1; i <= h0->GetNbinsX(); ++i) {
      if (h0->GetBinContent(i) < 0) h0->SetBinContent(i, 0);
      // bin which content is negative value will be 0 content
    }
  }
  OnePEiterate(fN1);  // 4th 1PE distribution 推定
  MakeNPEdist(kMakeNPE_Fast);      // 4th 2PE,3PE distribution
}

void SinglePEAnalyzer::MakeNPEdist(EMakeNPE mode) {
  // Showing just the basic idea to switch different methods
  // TODO: Fix me for future use
  switch(mode) {
  case kMakeNPE_Takahashi2018 :
    MakeNPEdistTakahashi2018();
    break;
  case kMakeNPE_Fast :
    MakeNPEdistFast();
    break;
  default:
    MakeNPEdistTakahashi2018();
    break;
  }
}

void SinglePEAnalyzer::MakeNPEdistFast() {
  // TODO: implement me
  // 12/11　ishida implemented you
  
  std::vector<std::shared_ptr<TH1D>> h{};
  for (int i=0; i<=kMaxPE; ++i) h.push_back(fNPEDist.back()[i]);
 
  int nbins = h[1]->GetNbinsX();
  int bin0 = h[1]->FindBin(0);

  for (int npe =2; npe<=kMaxPE; ++npe) {
    for (int ibin = 1; ibin <= nbins; ++ibin) {
      int i = ibin - bin0;
      double n2i = 0;
      for (int jbin = 1; jbin <= nbins; ++jbin) { // i' = j
        int j = jbin - bin0;
        double n1j = h[npe-1]->GetBinContent(jbin);
        int i_jbin = i - j + bin0;
        double n1i_j = (1 <= i_jbin && i_jbin <= nbins) ? h[1]->GetBinContent(i_jbin) : 0;
        n2i += (n1j * n1i_j) / (npe * fN0); // Eq. (11)
      }
      h[npe]->SetBinContent(ibin, n2i);
    }
  }
  std::cout <<"N1 to N"<<kMaxPE<< " is, " <<std::endl;
  for (int i=1; i<=kMaxPE; ++i) std::cout << h[i]->Integral(1,-1) << "\t" ;
  std::cout << std::endl;
  double Q1mean_7 = 0;
  auto xmin = fSignalH1->GetXaxis()->GetXmin();
  auto binwidth = fSignalH1->GetXaxis()->GetBinWidth(1);

  for (int i = 1; i <= nbins; ++i) {
    Q1mean_7 += (xmin + i * binwidth) * (h[1]->GetBinContent(i));
  }
  Q1mean_7 /= h[1]->Integral(1, -1);
  std::cout << h[1]->GetName() << "'s <Q1>(7) is " << Q1mean_7 << std::endl;

  double Q1mean_7err;
  for (int i = 1; i < h[1]->FindBin(0); ++i) {
    Q1mean_7err += pow((xmin + binwidth * i) * fLambda + Q1mean_7, 2);
  }
  for (int i = h[1]->FindBin(0); i <= nbins; ++i) {
    Q1mean_7err += pow((xmin + binwidth * i), 2) * h[1]->GetBinContent(i);
  }
  Q1mean_7err /= pow(h[1]->Integral(1, -1) * fLambda, 2);
  Q1mean_7err = sqrt(Q1mean_7err);

  std::cout << h[1]->GetName() << "'s <Q1>(7)error is " << Q1mean_7err
            << std::endl;
  std::cout << "lambda calc from (-log(fN0/Nonall):(Mean prop) = "<< fLambda << " : "<< fSignalH1->GetMean() / h[1]->GetMean() <<std::endl;
  std::cout << "lambda from entry = " << fLambda << " ± " << fLambda_err_fromN << std::endl;
  std::cout << "lambda from chage = " << fSignalH1->GetMean() / h[1]->GetMean() << " ± " << sqrt( pow(fSignalH1->GetMeanError()/h[1]->GetMean(),2) + pow(fSignalH1->GetMean()*h[1]->GetStdDev()/(h[1]->GetMean()*h[1]->GetMean()),2) /h[1]->Integral(1,-1) ) <<std::endl;
  //std::cout << "AmpGaining is " << GetAmpGain() << std::endl; 最終的にこれを見る？
}

void SinglePEAnalyzer::MakeNPEdistTakahashi2018() {
  auto h0 = fNPEDist.back()[0];
  auto h1 = fNPEDist.back()[1];
  auto h2 = fNPEDist.back()[2];
  auto h3 = fNPEDist.back()[3];

  int nbins = h1->GetNbinsX();
  int bin0 = h1->FindBin(0);

  // 2PE
  for (int ibin = 1; ibin <= nbins; ++ibin) {
    int i = ibin - bin0;
    double n2i = 0;
    for (int jbin = 1; jbin <= nbins; ++jbin) {  // i' = j
      int j = jbin - bin0;
      double n1j = h1->GetBinContent(jbin);
      int i_jbin = i - j + bin0;
      double n1i_j =
          (1 <= i_jbin && i_jbin <= nbins) ? h1->GetBinContent(i_jbin) : 0;
      n2i += (n1j * n1i_j) / (2 * fN0);  // Eq. (11)
    }
    h2->SetBinContent(ibin, n2i);
  }

  for (int ibin = 1; ibin <= nbins; ++ibin) {
    int i = ibin - bin0;
    double n3i = 0;
    for (int jbin = 1; jbin <= nbins; ++jbin) {
      int j = jbin - bin0;
      for (int kbin = 1; kbin <= nbins; ++kbin) {
        int k = kbin - bin0;
        double n1k = h1->GetBinContent(kbin);
        double n1j = h1->GetBinContent(jbin);
        int i_j_kbin = i - j - k + bin0;
        double n1_i_j_k = h1->GetBinContent(i_j_kbin);
        n3i += n1k * n1j * n1_i_j_k / (6 * fN0 * fN0);
      }
    }
    h3->SetBinContent(ibin, n3i);
    if (ibin % 50 == 1) printf("*");
    if (ibin % 500 == 1) printf("check500\n");
  }
  std::cout << "N2 is " << h2->Integral(1,-1) << "\tN3 is " << h3->Integral(1,-1) << std::endl;
  double Q1mean_7 = 0;
  auto xmin = fSignalH1->GetXaxis()->GetXmin();
  auto binwidth = fSignalH1->GetXaxis()->GetBinWidth(1);

  for (int i = 1; i <= nbins; ++i) {
    Q1mean_7 += (xmin + i * binwidth) * (h1->GetBinContent(i));
  }
  Q1mean_7 /= h1->Integral(1, -1);
  std::cout << h1->GetName() << "'s <Q1>(7) is " << Q1mean_7 << std::endl;

  double Q1mean_7err;
  for (int i = 1; i < h1->FindBin(0); ++i) {
    Q1mean_7err += pow((xmin + binwidth * i) * fLambda + Q1mean_7, 2);
  }
  for (int i = h1->FindBin(0); i <= nbins; ++i) {
    Q1mean_7err += pow((xmin + binwidth * i), 2) * h1->GetBinContent(i);
  }
  Q1mean_7err /= pow(h1->Integral(1, -1) * fLambda, 2);
  Q1mean_7err = sqrt(Q1mean_7err);

  std::cout << h1->GetName() << "'s <Q1>(7)error is " << Q1mean_7err
            << std::endl;

  std::cout << "lambda calc from (-log(fN0/Nonall):(Mean prop) = "<< fLambda << " : "<< fSignalH1->GetMean() / h1->GetMean() <<std::endl;
  std::cout << "lambda from entry = " << fLambda << " ± " << fLambda_err_fromN << std::endl;
  std::cout << "lambda from chage = " << fSignalH1->GetMean() / h1->GetMean() << " ± " << sqrt( pow(fSignalH1->GetMeanError()/h1->GetMean(),2) + pow(fSignalH1->GetMean()*h1->GetStdDev()/(h1->GetMean()*h1->GetMean()),2) /h1->Integral(1,-1) ) <<std::endl;
  std::cout << "AmpGaining is " << GetAmpGain() << std::endl;
}

void SinglePEAnalyzer::OnePEiterate(double N) {
  double sum = 0;
  auto h0 = fNPEDist.back()[0];
  auto h1 = fNPEDist.back()[1];

  for (int i = 1; i <= h0->GetNbinsX(); ++i) {
    double c = h0->GetBinContent(i);
    sum += c;
    if (sum <= N) {
      h1->SetBinContent(i, c);
    } else {
      h1->SetBinContent(i, N - (sum - c));
      break;
    }
  }
}

void SinglePEAnalyzer::CompareLambda_E_to_C() {
  //compare_flambda : cLambda;
  std::cout << "lambdaを出す二つの方法で違いが出ないかチェック" <<std::endl;
  int number_of_N = 0;
  double Total_Entries_E = 0;
  double Total_Entries_C = 0;
  double NPE_Entry_E = fN0;
  double NPE_Entry_C = fN0;
  do {
    Total_Entries_E += NPE_Entry_E;
    Total_Entries_C += NPE_Entry_C;
    std::cout << "N"<< number_of_N <<": " << NPE_Entry_E <<"\t"<< NPE_Entry_C << std::endl;
    ++number_of_N;
    NPE_Entry_E *= double(fLambda/number_of_N);
    NPE_Entry_C *= double(cLambda/number_of_N);
  } while ( NPE_Entry_E > 1 || NPE_Entry_C > 1);
  std::cout << "from Entries:Total Entries calc PE0 to PE" << number_of_N-1 <<"  is " << Total_Entries_E << std::endl;
  std::cout << "from Charge :Total Entries calc PE0 to PE" << number_of_N-1 <<"  is " << Total_Entries_C << std::endl;
  std::cout << "Total Entries calc from fSignal->Integral() is " << fSignalH1->Integral(1,-1) << std::endl;
}

double  SinglePEAnalyzer::GetAmpGain() {
  auto h1 = fNPEDist.back()[1];
  const double EC = 1.0; //elementary charge
  const double PreampGain = 1.0;
  return h1->GetMean() / (EC * PreampGain);
}

void SinglePEAnalyzer::Savedata_Lamb() {
  std::ofstream out;
  out.open("Lambs.txt", std::ios::app);
  out << fLambda << "\t" << fLambda_err_fromN << "\t" << cLambda << "\t" << cLambda_err
      << std::endl;
  out.close();
}

void Savedata_1PECharge(std::shared_ptr<SinglePEAnalyzer> ana,std::shared_ptr<SinglePEAnalyzer> ana2) {
  //Ent,Mean,MeanErr,StdD,StdDErr,Ent...(10 properties)
  std::ofstream out;
  out.open("Charge_HV.txt", std::ios::app);
  out << ana->GetNPEDist(1)->Integral(1,-1) << "\t" << ana->GetNPEDist(1)->GetMean() << "\t" << ana->GetNPEDist(1)->GetMeanError()
  << "\t" << ana->GetNPEDist(1)->GetStdDev() << "\t" << ana->GetNPEDist(1)->GetStdDevError()
  << "\t" << ana2->GetNPEDist(1)->Integral(1,-1) << "\t" << ana2->GetNPEDist(1)->GetMean() << "\t" << ana2->GetNPEDist(1)->GetMeanError()
  << "\t" << ana2->GetNPEDist(1)->GetStdDev() << "\t" << ana2->GetNPEDist(1)->GetStdDevError() << std::endl;
  out.close();
}

void Savedata_Lamb_and_Charge(std::shared_ptr<SinglePEAnalyzer> ana,std::shared_ptr<SinglePEAnalyzer> ana2) {
  //fLamb,fLambErr,cLamb,cLambErr,[ana]Ent,Mean,MeanErr,StdD,StdDErr,[ana]Ent...(14 properties)
  std::ofstream out;
  out.open("data.txt", std::ios::app);
  out << ana->GetLambda() << "\t" << ana->GetLambdaErr() << "\t" << ana->GetLambda_C() << "\t" << ana->GetLambda_CErr() 
  << "\t" << ana->GetNPEDist(1)->Integral(1,-1) << "\t" << ana->GetNPEDist(1)->GetMean() << "\t" << ana->GetNPEDist(1)->GetMeanError()
  << "\t" << ana->GetNPEDist(1)->GetStdDev() << "\t" << ana->GetNPEDist(1)->GetStdDevError()
  << "\t" << ana2->GetNPEDist(1)->Integral(1,-1) << "\t" << ana2->GetNPEDist(1)->GetMean() << "\t" << ana2->GetNPEDist(1)->GetMeanError()
  << "\t" << ana2->GetNPEDist(1)->GetStdDev() << "\t" << ana2->GetNPEDist(1)->GetStdDevError() << std::endl;
  out.close();
}

std::pair<std::shared_ptr<SinglePEAnalyzer>, std::shared_ptr<SinglePEAnalyzer>> test(int Iteration_MaxPE = 4,const std::string& file_name = READPATH) {
  auto ana = std::make_shared<SinglePEAnalyzer>();
  ana->ReadFile(file_name, "signal", "noise", 5);
  //ReadFileでfileから取得→Rebinまでやってる
  ana->SetkMaxPE(Iteration_MaxPE);
  ana->Analyze();
  
  std::cout << "\n**border**\n" << std::endl;

  auto ana2 = std::make_shared<SinglePEAnalyzer>();
  ana2->MakeMCSpectra(ana->GetNoiseH1(), ana->GetNPEDist(1), ana->GetLambda(), ana->GetSignalH1()->GetEntries(), true);
  ana2->SetkMaxPE(Iteration_MaxPE);
  ana2->Analyze();

  std::cout << "\nanaのLambda比較 fLamb : cLamb "<<std::endl;
  std::cout << ana->GetLambda() << " ± " << ana->GetLambdaErr() << "\t" << ana->GetLambda_C() << " ± " << ana->GetLambda_CErr() << std::endl;
  std::cout << "Signalとana2のCharge比較 Signal : ana2 "<<std::endl;
  std::cout << ana->GetSignalH1()->GetMean() << " ± " << ana->GetSignalH1()->GetMeanError() << "\t" << ana2->GetSignalH1()->GetMean() << " ± " << ana2->GetSignalH1()->GetMeanError() << std::endl;
  
  std::cout << "\nanaとana2のCharge比較"<<std::endl;
  std::cout << "ana1 1PE dist Entry:" << ana->GetNPEDist(1)->Integral(1,-1) << std::endl;
  std::cout << "ana1 1PE dist Mean :" << ana->GetNPEDist(1)->GetMean() <<" ± " << ana->GetNPEDist(1)->GetMeanError() << std::endl;
  std::cout << "ana1 1PE dist StdDv:" << ana->GetNPEDist(1)->GetStdDev() <<" ± " << ana->GetNPEDist(1)->GetStdDevError() << std::endl;
  std::cout << "ana2 1PE dist Entry:" << ana2->GetNPEDist(1)->Integral(1,-1) << std::endl;
  std::cout << "ana2 1PE dist Mean :" << ana2->GetNPEDist(1)->GetMean() <<" ± " << ana2->GetNPEDist(1)->GetMeanError() << std::endl;
  std::cout << "ana2 1PE dist StdDv:" << ana2->GetNPEDist(1)->GetStdDev() <<" ± " << ana2->GetNPEDist(1)->GetStdDevError() << std::endl;
  Savedata_Lamb_and_Charge(ana,ana2);
  return std::make_pair(ana, ana2);

  //欲しい情報
  //anaのみでLambda比較: fLambda(calc from Entry = -log(fN0 / Nonall)),cLambda(calc from charge = <Q_all>/<Q_1>)
  //anaとana2でCharge比較:Ent,Mean,MeanErr,StdD,StdDErrを二つで
}

std::shared_ptr<SinglePEAnalyzer> Draw_raw_dist(const std::string& file_name = READPATH) {
  TCanvas *fcan = new TCanvas("fcan", "fcan", 1400, 600);
  gPad->SetLogy(1);
  auto ana = std::make_shared<SinglePEAnalyzer>();
  ana->ReadFile(file_name, "signal", "noise", 5);
  ana->GetSignalH1()->SetTitle("Charge Distribution");
  ana->GetSignalH1()->SetLineColor(6);
  ana->GetSignalH1()->GetYaxis()->SetRangeUser(0.8,10000);
  ana->GetSignalH1()->Draw("same");
  ana->GetNoiseH1()->Draw("same");
  ana->MakefSuminus0();
  std::cout << ana->GetAlpha() <<std::endl;
  ana->GetSuminus0()->SetFillColor(10);
  ana->GetSuminus0()->SetLineColor(12);
  ana->GetSuminus0()->Draw("same");
  TLegend *leg = new TLegend( 0.60, 0.55, 0.90, 0.70) ; 
  leg->AddEntry((TObject*)0, "Signal(pink)","");
  leg->AddEntry((TObject*)0 , "Noise(blue)","");
  leg->AddEntry((TObject*)0 , "Suminus0(gray)","");
  leg->Draw();
  return ana;
}

void test_HV() {
  //this function cannot work propery
  //TODO:fix this error
  std::vector<std::string> filenames = {};
  std::string nameH = "PMT_HVCharge/turn";
  std::string nameF = "_CH1.root";
  for (int i = 10; i<=15; ++i) {
    std::string nameB = std::to_string(i);
    std::string fname = nameH + nameB + nameF;
    std::cout << fname << std::endl; 
    // auto anas = test(5,filenames) //this is cause of the error
    filenames.push_back(fname);
  }
}