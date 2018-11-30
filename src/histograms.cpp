#ifndef HISTOGRAMS_C
#define HISTOGRAMS_C

#include "histograms.hpp"
#include "common_functions.C"
#include "TSystem.h"
#include "TCanvas.h"

#define kLOGZ  1
#define kLOGY  2
#define kLOGX  4
#define kSTATS 8

Histograms::Histograms(int nchan)
{
  nch = nchan;
  xmin.resize(nch);
  xmax.resize(nch);
  ymin.resize(nch);
  ymax.resize(nch);
  nbinsx.resize(nch);
  nbinsy.resize(nch);
};
Histograms::~Histograms()
{
  delete hChi2;
  delete hTOA;
  
  for (int i=0; i<nch; i++){
    delete hTracks[i];
    delete hC[i];
    delete hN[i];
    delete hPassed[i];
    delete hEff[i];
    delete h2Charge[i];
  }
};

void Histograms::OpenHistos()
{
  hChi2 = new TH1I("hChi2", "Track #Chi2; Track #Chi^2 / NDoF; n_{tracks}", 500, 0, 10);
  hTOA  = new TH1I("hTOA",  "Time of acknowledge signal; Time (bin); n_{events}", 500, 0, 1000);
  hCluSize = new TH1I("hCluSize", "Cluster Size; cluster size; events", 6, -0.5, 5.5);
  
//   for (int i=0; i<nch; i++){
//     nbinsx[i] = 100;
//     nbinsy[i] = 100;
//   }
  
  h2CluSize1    = new TH2F("h2CluSize1",    "Cluster size 1; x(#mum); y(#mum)", nbinsx[0], xmin[0], xmax[0], nbinsy[0], ymin[0], ymax[0]);
  h2CluSize2    = new TH2F("h2CluSize2",    "Cluster size 2; x(#mum); y(#mum)", nbinsx[0], xmin[0], xmax[0], nbinsy[0], ymin[0], ymax[0]);
  h2CluSizeGT2  = new TH2F("h2CluSizeGT2",  "Cluster size > 2; x(#mum); y(#mum)", nbinsx[0], xmin[0], xmax[0], nbinsy[0], ymin[0], ymax[0]);
  
  for (int i=0; i<nch; i++){
    hC.push_back(new TH1I(Form("hC%d",i), Form("Charge spectrum ch. %d; U (mV); count", i), 100, -5, 30)); 
    hN.push_back(new TH1I(Form("hN%d",i), Form("Noise spectrum ch. %d; U (mV); count", i), 100, -5, 5)); 
    hTracks.push_back(new TH2I(Form("hTracks%d",i), Form("Tracks ch. %d; x (#mum); y (#mum); count", i), nbinsx[i], xmin[i], xmax[i], nbinsy[i], ymin[i], ymax[i]));
    hPassed.push_back(new TH2I(Form("hPassed%d",i), Form("Passed tracks ch. %d; x (#mum); y (#mum); count", i), nbinsx[i], xmin[i], xmax[i], nbinsy[i], ymin[i], ymax[i]));
    hEff.push_back(new TH2F(Form("hEff%d",i), Form("Efficiency ch. %d; x (#mum); y (#mum); efficiency", i), nbinsx[i], xmin[i], xmax[i], nbinsy[i], ymin[i], ymax[i]));
    h2Charge.push_back(new TH2F(Form("h2Charge%d",i), Form("Charge map ch. %d; x (#mum); y (#mum); U (mV)", i), nbinsx[i], xmin[i], xmax[i], nbinsy[i], ymin[i], ymax[i]));
    h2AvgCluSize.push_back(new TH2F(Form("h2AvgCluSize%d",i),  Form("Average cluster size ch. %d; x(#mum); y(#mum)", i), nbinsx[0], xmin[0], xmax[0], nbinsy[0], ymin[0], ymax[0]));
  }
    return;
};

void Histograms::SetHistosRange(int ch, float _xmin, float _xmax, float _ymin, float _ymax, int _nbinsx, int _nbinsy)
{
  xmin[ch] = _xmin;
  xmax[ch] = _xmax;
  ymin[ch] = _ymin;
  ymax[ch] = _ymax;
  nbinsx[ch] = _nbinsx;
  nbinsy[ch] = _nbinsy;
  return;
}

void Histograms::FinalProcessing()
{
  for (int i=0; i<nch ; i++){
    hEff[i]->Reset();   // = (TH2F*)(hPassed[i]->Clone(Form("hEff%d",i)));
    hEff[i]->Add(hPassed[i]);
    hEff[i]->Divide(hTracks[i]);
    
    h2Charge[i]->Divide(hTracks[i]);
    h2AvgCluSize[i]->Divide(hTracks[i]);
  }
  
  return;
}
  
void Histograms::Plotting(TString path)
{
  vector<int> index_sig_spectrum;
  vector<int> index_noise_spectrum;
  vector<short> options;
  
  vector<TH1*> v1;
  v1.push_back(hChi2);    options.push_back(kSTATS);
  v1.push_back(hTOA);     options.push_back(kSTATS);
  v1.push_back(hCluSize); options.push_back(kSTATS);
  for(uint i=0; i<hC.size(); i++){
    index_sig_spectrum.push_back(v1.size());
    v1.push_back(hC[i]);  options.push_back(kSTATS + kLOGY);
    index_noise_spectrum.push_back(v1.size());
    v1.push_back(hN[i]);  options.push_back(kSTATS + kLOGY);
  }
  
  vector<TH2*> v2;
  v2.push_back(h2CluSize1);
  v2.push_back(h2CluSize2);
  v2.push_back(h2CluSizeGT2);
  for(uint i=0; i<hTracks.size(); i++){
    v2.push_back(hTracks[i]);
    v2.push_back(hPassed[i]);
    v2.push_back(hEff[i]);
    v2.push_back(h2Charge[i]);
    v2.push_back(h2AvgCluSize[i]);
  }
  
  gSystem->MakeDirectory(path.Data());
    
  TCanvas *c = new TCanvas("c", "c", 1);
  c->SetMargin(0.15,0.18,0.15,0.1);
  for(uint i=0; i<v1.size(); i++) {
    FormatHisto1D(v1[i], v1[i]->GetTitle());
    v1[i]->SetStats(options[i] & kSTATS);
    v1[i]->Draw();
    TString fname_pdf = path + v1[i]->GetName() + ".pdf";
    TString fname_png = path + v1[i]->GetName() + ".png";
//     cout << fname_pdf << endl;
    c->Print(fname_pdf);
    c->Print(fname_png);
  }
  
//   for(uint i=0; i<v2.size(); i++) {
//     FormatHisto2D(v2[i], v2[i]->GetTitle());
//     v2[i]->Draw("COLZ");
//     TString fname_pdf = path + v2[i]->GetName() + ".pdf";
//     TString fname_png = path + v2[i]->GetName() + ".png";
//     c->Print(fname_pdf);
//     c->Print(fname_png);
//   }
   
//   for(auto h : v1) {
//     FormatHisto1D(h, h->GetTitle());
//     h->SetStats();
//     h->Draw();
//     TString fname_pdf = path + h->GetName() + ".pdf";
//     TString fname_png = path + h->GetName() + ".png";
//     cout << fname_pdf << endl;
//     c->Print(fname_pdf);
//     c->Print(fname_png);
//   }
//   
  for(auto h : v2) {
    FormatHisto2D(h, h->GetTitle());
    h->Draw("COLZ");
    TString fname_pdf = path + h->GetName() + ".pdf";
    TString fname_png = path + h->GetName() + ".png";
    c->Print(fname_pdf);
    c->Print(fname_png);
  }
   
  // Draw spectra together with the threshold 
  for (uint i=0; i<index_sig_spectrum.size(); i++){
    if (i >= thr.size()) {
      cout << "Plotting thresholds: no entry for ch. " << i << ". Stopping" << endl; 
      break;
    }
    int is=index_sig_spectrum[i];
    int in=index_noise_spectrum[i];
    TLine l;
    TLatex tex;
    
    // signal spectrum
    v1[is]->SetStats(0);
    v1[is]->Draw();
    l.DrawLine(thr[i], 0, thr[i], 0.7*v1[is]->GetMaximum());
    tex.DrawLatex(1.15*thr[i], 0.5*v1[is]->GetMaximum(), Form("thr. = %.1lf mV", thr[i]));
    TString fname_pdf = path + v1[is]->GetName() + "_thr.pdf";
    TString fname_png = path + v1[is]->GetName() + "_thr.png";
//     cout << fname_pdf << endl;
    c->Print(fname_pdf);
    c->Print(fname_png);
    
    // noise spectrum
//     v1[in]->SetStats(0);
    v1[in]->Draw();
    l.DrawLine(thr[i], 0, thr[i], 0.7*v1[in]->GetMaximum());
    tex.DrawLatex(1.15*thr[i], 0.5*v1[in]->GetMaximum(), Form("thr. = %.1lf mV", thr[i]));
    fname_pdf = path + v1[in]->GetName() + "_thr.pdf";
    fname_png = path + v1[in]->GetName() + "_thr.png";
//     cout << fname_pdf << endl;
    c->Print(fname_pdf);
    c->Print(fname_png);
  }
   
  delete c;
  return;
}

void Histograms::SaveThresholdValue(int ch, float value)
{
  if (thr.size() > (uint)(ch)) thr[ch] = value;
  else{
    thr.resize(ch+1); thr[ch] = value;
  }
  return;
}
  
#endif