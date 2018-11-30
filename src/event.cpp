#ifndef EVENT_C
#define EVENT_C

#include "event.hpp"

event::event(int nch)
{
//   for (int i=0;i<DRS_N_POINTS; i++) xaxis[i]=1.0*i;
  for (int i=0;i<nch; i++) {
    wf.push_back(new waveform);
  }
  wfAcknowledge = new waveform();
  charge.resize(nch);
  noise.resize(nch);
  hit.resize(nch);
  clu_size.resize(nch);
  
  for (int i=0; i<DRS_N_POINTS; i++) t[i]=i;
};
event::~event()
{
//   for (w : wf) delete w;
  for (uint i=0; i< wf.size(); i++) delete wf[i];
  delete wfAcknowledge;
  
};

// void event::Draw(int ch, char* option)
// {
// //   Float_t x[DRS_N_POINTS];
// //   for (int i=0; i<DRS_N_POINTS; i++) x[i] = 1.0*i;
//   
//   float labelSize = 0.05;
//   if(gr[ch]) delete gr[ch];
// //   gr[ch] = new TGraph(DRS_N_POINTS, xaxis, wf[ch].wf);
//   gr[ch] = new TGraphErrors(DRS_N_POINTS, xaxis, wf[ch].wf,0,0);
//   
//   gr[ch]->SetLineColor(1+ch);
//   gr[ch]->SetLineWidth(2);
//   gr[ch]->SetMarkerColor(1+ch);
//   gr[ch]->GetYaxis()->SetRangeUser(-0.05,0.2);
//   gr[ch]->SetTitle(" ;t (ns); U (V)");
//   
//   gr[ch]->GetXaxis()->SetTitleSize(labelSize);
//   gr[ch]->GetXaxis()->SetTitleOffset(1.);
//   gr[ch]->GetXaxis()->SetLabelSize(labelSize);
//   gr[ch]->GetYaxis()->SetTitleSize(labelSize);
//   gr[ch]->GetYaxis()->SetLabelSize(labelSize);
//   gr[ch]->GetYaxis()->SetTitleOffset(1.);
//   
//   if (strlen(option)) gr[ch]->Draw(option);
//   else gr[ch]->Draw("APL");
//   return; //*/
// };
// 
  
  void event::BaselineCorrection(int ch, int t_min, int t_max)
  {
    float offset=0;
    for (int i=t_min; i<t_max; i++){
      offset += wf[ch]->wf[i];
    }
    offset /= (t_max-t_min);
    for (int i=0; i<DRS_N_POINTS; i++) wf[ch]->wf[i] -= offset;
    
//     cout << "ch. " << ch << " baseline correction: " << offset << " mV" << endl;
    return;
  }
  
  void event::ShiftPhase(int ch, int phase)
  {
    if (phase >= 0){
      for (int i=TMath::Max(0,-phase); i<DRS_N_POINTS-TMath::Max(0,phase); i++){
        wf[ch]->wf[i] = wf[ch]->wf[i+phase];
      }
    }
    else{
      for (int i=DRS_N_POINTS-TMath::Max(0,phase); i>TMath::Max(0,-phase); i--){
        wf[ch]->wf[i] = wf[ch]->wf[i+phase];
      }
    }
    
    return;
  }
  
  float event::GetCharge(int ch, int t_min, int t_max, int mode)
  {
    float _charge=-1111;
    switch(mode) {
      case 0:
        // integral
        _charge = wf[ch]->Integral(t_min, t_max) / (1+t_max-t_min);
        break;
      case 1:
        // Step size in CSA
        _charge = wf[ch]->GetMax(t_min, t_max) - wf[ch]->GetMin(t_min, t_max);
        break;
      case 2:
        // For CSA: fit two linear functions y=kx+n with the same k before and after the step. Charge is then n2-n1
        // Range of fit before the step: [tmin, tmin+dt]
        // Range of fit after  the step: [tmax, tmax+dt]
        {
          // new scope to avoid the "crosses initialization" error
          int dt=20;  // fitting range (number of points)
          TF1 f1("f1","pol1",t_min,t_min+dt);
          TF1 f2("f2","[0]+[1]*x",t_max,t_max+dt);
  
          TGraph g(DRS_N_POINTS, t, wf[ch]->wf);
          g.Fit("f1","RQ");
          f2.FixParameter(1,f1.GetParameter(1));
          g.Fit("f2","RQ");
          _charge = f1.GetParameter(0) - f2.GetParameter(0);
        }
        break;
      case 3:
        // Get the difference in the voltage between bins t_min and t_max
        // N bins in the interval [t_min, t_min+N-1] and [t_max, t_max+N-1] are averaged
        { 
          int N=10;
          float v1=0, v2=0;
          for (int j=0; j<N; j++) {
            v1+=wf[ch]->wf[t_min+j];
            v2+=wf[ch]->wf[t_max+j];
          }
          _charge=(v2-v1)/N;
        }
        break;
      default:
        // integral
        _charge = wf[ch]->Integral(t_min, t_max) / (1+t_max-t_min);
//         break;
    }
    
//     cout << _charge << endl;
    return _charge;
  }

  int event::GetTOA(int ch, float thr)
  {
//     return 111111;
    return wf[ch]->GetFirstAboveThr(thr, -1, 900, DRS_N_POINTS-10);
  }

  short event::GetCluSizeDRS()
  {
    short _clusize=0;
    for (auto h : hit){
      if (h == true) _clusize++;
    }
    clu_size_drs = _clusize;
    for (uint i=0; i<clu_size.size(); i++) clu_size[i] = _clusize;
    return _clusize;
  }
  
  int event::GetVetoPosition(int ch, float thr, float polarity)
  {
    // Searches for the rising/falling edge of veto signals in one of the channels
    // ch ... channel with the veto signal
    // thr ... threshold value 
    // polarity ... +1:rising edge, -1: falling edge
    
    return wf[ch]->GetFirstAboveThr(thr, polarity);
  }
  
  int event::GetAcknowledgePosition(int ch)
  {
    // Searches for the falling edge of the telescope acknowledge signal in one of the channels
    
    return wf[ch]->GetFirstAboveThr(0.3, -1);
  }
  
  void event::SubtractPickup(int ch_sig, int ch_sub, int phase)
  {
    for (int i=TMath::Max(0,-phase); i<DRS_N_POINTS-TMath::Max(0,phase); i++){
//       wf[2]->wf[i] = wf[ch_sig]->wf[i] - wf[ch_sub]->wf[i] - 0.05;
//       wf[ch_sig]->wf[i] = wf[ch_sig]->wf[i] - wf[ch_sub]->wf[i];
//       wf[2]->wf[i]      = wf[ch_sig]->wf[i] - wf[ch_sub]->wf[i+phase] - 0.05;
      wf[ch_sig]->wf[i] = wf[ch_sig]->wf[i] - wf[ch_sub]->wf[i+phase];
    }
    return;
  }

  
  // bool event::VetoPresent(int ch, float thr, int mode)
// {
//   // Seaches for a veto signal in one of waveforms 
//   // Veto signal has a box form so it can be discovered simply by finding the minimum and maximum of the wf
//   // mode : 
// //   0 ... get minimum
// //   1 ... get difference between min and max
//   float delta;
//   switch (mode){
//     case 1:
//       delta = (wf[ch].GetMax() - wf[ch].GetMin());
//       if ( delta > thr) return true;
//       break;
//     case 2:
//       delta = wf[ch].GetMax();
//       if ( delta > thr) return true;
//       break;
//     default:
//       delta = wf[ch].GetMin();
//       if ( delta < thr) return true;
//       break;
//   }  
//    
// //   cout << delta << endl;
//   return false;
// }
// 

#endif