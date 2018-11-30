#ifndef RUN_H
#define RUN_H

#include "TFile.h"
#include "TTree.h"

#include "waveform.cpp"
#include "event.cpp"
#include "histograms.cpp"

#include <vector>

#ifndef MAX_CHANNELS 
#define MAX_CHANNELS 8
#endif

class run
{
public:
  run(const char* fnameDRS, int nchannels);
  run(const char* fnameDRS, const char* fnameTel, int nchannels);
  ~run();
  
  TFile* fDRS;
  TFile* fTel;
  TTree* tTracks;
  TTree* tWaveforms;
  
  Histograms *histos;
  vector<waveform*> wfavg;
  event* Event;
//   analyzed_event AnaEvent;
  
//   void FillAnalyzedEvent(int itr[2], float charge[NCH_DRS]);
  
  int nch;
  int N;
  int ts_min[MAX_CHANNELS], ts_max[MAX_CHANNELS];   //! sampling time for signal
  int tn_min[MAX_CHANNELS], tn_max[MAX_CHANNELS];   //! sampling time for noise
  float thr[MAX_CHANNELS];
  float zDUT[MAX_CHANNELS];
  int phase_correction[MAX_CHANNELS];
  
  void GetEvent(int ie); 
  void GetEventDRS(int ie); 
  void FindSamplingTime(int ch, int nsamples, int t_min, int t_max, int pulse_polarity=+1, int bin_first=0, int bin_last=DRS_N_POINTS); 
  float FindThreshold(int ch, int nsamples, float sigma, int pulse_polarity=+1, int mode=0); 
  vector<float>* GetChargeAll(int mode=0, int pulse_polarity=+1);
  vector<float>* GetNoiseAll (int mode=0, int pulse_polarity=+1);
  void FillHistograms();
  
  void GetAvgWaveformAsync(int ch, int ch_veto, int Nevents);
//   void SetSignalWindow(int ch, float min, float max);
//   void SetNoiseWindow(int ch, float min, float max);
//   float FindThreshold(int ch, int nsamples=1000, float sigma=4., int mode=0);
//   float GetEfficiency(int ch, float thr, int mode=0, int pulse_polarity=+1);
//   float GetEfficiencyClustered(float *thr, int mode=0, int pulse_polarity=+1);
//   float GetNoise(int ch, int mode=0);
// 
//   window ROI;
//   void SetROI(float xmin, float xmax, float ymin, float ymax, int nbinsx, int nbinsy=-1);
//   
  
//   DUTCorrelation dutcorr;
  
  
};

#endif