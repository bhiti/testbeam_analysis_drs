#ifndef WAVEFORM_C
#define WAVEFORM_C

#include "waveform.hpp"

waveform::waveform(){ for (int i=0; i<DRS_N_POINTS; i++) wf[i]=-1; };
waveform::~waveform(){};

Float_t waveform::GetMin(int first, int last)
{ 
  if (first < 0)      return TMath::MinElement(DRS_N_POINTS, wf);
  else if (last < 0)  return TMath::MinElement(DRS_N_POINTS-first, wf+first);
  else                return TMath::MinElement(last-first, wf+first);
};

Float_t waveform::GetMax(int first, int last)
{ 
  if (first < 0)      return TMath::MaxElement(DRS_N_POINTS, wf);
  else if (last < 0)  return TMath::MaxElement(DRS_N_POINTS-first, wf+first);
  else                return TMath::MaxElement(last-first, wf+first);
};

Int_t waveform::GetMinPos(int first, int last)
{ 
  if (first < 0)      return (TMath::LocMin(DRS_N_POINTS, wf));
  else if (last < 0)  return (TMath::LocMin(DRS_N_POINTS-first, wf+first) + first );
  else                return (TMath::LocMin(last-first, wf+first) + first );
};

Int_t waveform::GetMaxPos(int first, int last)
{ 
  if (first < 0)      return (TMath::LocMax(DRS_N_POINTS, wf) );
  else if (last < 0)  return (TMath::LocMax(DRS_N_POINTS-first, wf+first) + first);
  else                return (TMath::LocMax(last-first, wf+first) + first);
};

Float_t waveform::Integral(int bin_start, int bin_end)
{
  Float_t integral=0;
  for (int i=bin_start; i<=bin_end; i++) integral+=wf[i];
  return integral;
}
// Float_t waveform::MaxStep(float thr, int first, int last)
// {
//   // Get difference between max and min of the waveform
//   // thr:   threshold for step detection (minimal delta V)
//   // first: first sampled bin
//   // last:  last sampled bin
//   float dV=0;   // maximal delta U
//   dV = GetMax(first, last) - GetMin(first, last);
//   if (dV > thr) return dV;
//   return 0.;
// }
// Float_t waveform::MaxStepN(int N, int sign, float thr, int dx, int first, int last)
// {
//   // Get difference between max and min of the waveform
//   // Max and min are obtained by integrating over N consequtive points
//   // N:     number of integrated points
//   // sign:  step polarity (+1/-1)
//   // thr:   threshold for step detection (minimal delta V)
//   // dx:    increment step
//   // first: first sampled bin
//   // last:  last sampled bin
//   float sum=0;   // maximal delta U
//   float min=22222;
//   float max=-22222;
//   for (int i=TMath::Min(first, (int)(DRS_N_POINTS-N-10)); i<last-N; i+=dx)
//   {
//     sum=0;
//     for (int j=0; j<N; j++) sum+=wf[i+j];
//     sum/=N;
//     if (i==first) {
//       min=sum; max=sum;
//     }
//     if (sum>max) max=sum;
//     else if (sum<min) min=sum;
//   }
// //   cout << first << " " << last << " " << max << " " << min << endl; 
//   float dV=max-min;
// //   cout << first << " " << dV << endl;
//   if (dV > thr) return dV;
//   return 0.;
// }
// Float_t waveform::SearchStep(int sign, float thr, int width, int step, int first, int last)
// {
//   // search waveform for a step between two straight lines
//   // sign:  step polarity (+1/-1)
//   // thr:   threshold for step detection (minimal delta V)
//   // width: maximal width of the step (high-pass threshold) in bins
//   // step:  step (in bins) at which the waveform is sampled
//   // first: first sampled bin
//   // last:  last sampled bin
//   float dV=0;   // maximal delta U
//   dV = GetMax(first, last) - GetMin(first, last);
//   if (dV > thr) return dV;
//   return 0.;
// }
// Float_t waveform::GetDelta(int first, int last, int N)
// {
//   // Get the difference in the voltage between bins first and last
//   // N bins in the interval [first, first+N-1] and [last, last+N-1] are averaged
//   float v1=0, v2=0;
//   for (int j=0; j<N; j++) {
//     v1+=wf[first+j];
//     v2+=wf[last+j];
//   }
//   float delta=(v2-v1)/N;
// //   cout << first << " " << last << " " << delta << " " << wf[last] - wf[first] << endl;;
//   return delta;
// }

Int_t waveform::GetFirstAboveThr(float thr, int polarity, int first, int last)
{
  int p = (polarity > 0) - (polarity < 0);  // sign of polarity +1, 0, -1
  for (int i=TMath::Max(1,first); i<TMath::Min(last, DRS_N_POINTS); i++){
    if ( p*wf[i-1] < TMath::Abs(thr) && p*wf[i] > TMath::Abs(thr) ) 
      return i;
  }
  return -1111;
}

TGraph* waveform::Draw()
{
  float t[DRS_N_POINTS];  
  for(int ipt=0; ipt<DRS_N_POINTS; ipt++) 
  {
//     wf[ipt] /= r->N;
    t[ipt] = ipt;
  }  
  TGraph* g = new TGraph (DRS_N_POINTS, t, wf);  
  g->GetXaxis()->SetTitle("t (0.33 ns)");
  g->GetYaxis()->SetTitle("U (mV)");
  g->SetTitle("Average waveform");
//   DrawTH1(g->GetHistogram(), "");  // Cosmetics for plotting
  return g;
}


#endif