#ifndef EVENT_H
#define EVENT_H

#ifndef MAX_TRACKS 
#define MAX_TRACKS 100
#endif

#ifndef MAX_CHANNELS 
#define MAX_CHANNELS 8
#endif

#include "waveform.cpp"
// #include "event.cpp"

#include "TF1.h"
#include "TGraph.h"

class event
{
public:
  event(int nch);
  ~event();
  
  vector<waveform*> wf;
  waveform* wfAcknowledge;
  float t[DRS_N_POINTS];  // xaxis for waveforms - if needed
  
  int    ntracks;
  int    TOA;
  double trk_x[MAX_TRACKS];
  double trk_y[MAX_TRACKS];
  double trk_slopex[MAX_TRACKS];
  double trk_slopey[MAX_TRACKS];
  double trk_chi2[MAX_TRACKS];
  double trk_redchi2[MAX_TRACKS];
  int    trk_dof[MAX_TRACKS];
  
  // intercept coordinates
  double x[MAX_CHANNELS][MAX_TRACKS];   
  double y[MAX_CHANNELS][MAX_TRACKS];
  
  vector<float> charge;
  vector<float> noise;
  vector<bool>  hit;
  vector<short> clu_size;
  short clu_size_drs;
  
  void SubtractPickup(int ch_sig, int ch_sub, int phase);
  void ShiftPhase(int ch, int phase);

  float GetCharge(int ch, int t_min, int t_max, int mode=0);
  int   GetTOA(int ch, float thr);
  short GetCluSizeDRS();
  void  BaselineCorrection(int ch, int t_min, int t_max);
  
  int GetVetoPosition(int ch, float thr, float polarity=+1);
  int GetAcknowledgePosition(int ch);
  
  //   void Draw(int ch, char* option="");
//   bool VetoPresent(int ch, float thr, int mode=0); 
// //   bool VetoPresent(int ch, float thr, int first=0, int last=NCH_DRS);
};

// struct analyzed_event
// {
//   float wf[NCH_DRS][NPTS_ANAWF];
//   Float_t charge[NCH_DRS];
//   
//   Double_t xt;
//   Double_t yt;
//   Double_t chi2;
//   Double_t slopeX;
//   Double_t slopeY;
//   int nt;
//   
//   Int_t ToT;    
//   Double_t xa;    // Cluster coordinates from PixX, PixY, ToT
//   Double_t ya;
//   UShort_t clustersize;
//   int na;
// };
#endif