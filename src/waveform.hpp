#ifndef WAVEFORM_H
#define WAVEFORM_H

#ifndef DRS_N_POINTS
#define DRS_N_POINTS 1024
#endif

#include "TGraph.h"
#include "TMath.h"
#include "TAxis.h"

#include <iostream>

class waveform
{
public:
  waveform();
  ~waveform();

  Float_t wf[DRS_N_POINTS];
  Float_t GetMin(int first=-1, int last=-1);
  Float_t GetMax(int first=-1, int last=-1);
  Int_t GetMinPos(int first=-1, int last=-1);
  Int_t GetMaxPos(int first=-1, int last=-1);
  
  Int_t GetFirstAboveThr(float thr, int polarity=-1, int first=1, int last=DRS_N_POINTS);
  Float_t Integral(int first, int last);
//   Float_t MaxStep(float thr, int first=0, int last=DRS_N_POINTS-1);
//   Float_t MaxStepN(int N, int sign, float thr, int dx, int first=0, int last=DRS_N_POINTS-1);
//   Float_t SearchStep(int sign, float thr, int width, int step, int first=0, int last=DRS_N_POINTS-1);
//   Float_t GetDelta(int bin_start, int bin_end, int N=10);
  TGraph *Draw();
};

#endif