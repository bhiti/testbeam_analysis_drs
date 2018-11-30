#ifndef HISTOGRAMS_H
#define HISTOGRAMS_H

#include "TH1.h"
#include "TH2.h"
#include "TString.h"
#include "TLine.h"
#include "TLatex.h"

#include <vector>

class Histograms
{
public:
  Histograms(int nch);
  ~Histograms();
  
  int nch;
  vector<float> xmin;
  vector<float> xmax;
  vector<float> ymin;
  vector<float> ymax;
  vector<int> nbinsx;
  vector<int> nbinsy;
  
  TH1I* hChi2;
  TH1I* hTOA;
  TH1I* hCluSize;
  TH2F* h2CluSize1;
  TH2F* h2CluSize2;
  TH2F* h2CluSizeGT2;
  
  vector<TH1I*> hC;
  vector<TH1I*> hN;
  vector<TH2I*> hTracks;
  vector<TH2I*> hPassed;
  vector<TH2F*> hEff;
  vector<TH2F*> h2Charge;
  vector<TH2F*> h2AvgCluSize;
  
  vector<float> thr;

  void OpenHistos();
  void SetHistosRange(int nch, float xmin, float xmax, float ymin, float ymax, int nbinsx=100, int nbinsy=100);
  void FinalProcessing();
  void Plotting(TString path);
  void SaveThresholdValue(int ch, float value);
};

#endif