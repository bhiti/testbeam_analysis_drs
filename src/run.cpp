#ifndef RUN_C
#define RUN_C
#include "run.hpp"

run::run(const char* fnameDRS, int nchannels)
{
  nch = nchannels;
  for (int i=0;i<nch; i++) wfavg.push_back(new waveform());
  Event = new event(nch);
  
  fDRS = new TFile(fnameDRS, "READ");
//   fTel = new TFile(fnameTel, "READ");
  
  tWaveforms = (TTree *)fDRS->Get("Plane0/Waveforms");
//   tTracks = (TTree *)fTel->Get("Tracks");
  
  for (int i=0; i<nch; i++){
    if(tWaveforms->SetBranchAddress(Form("waveform%d",i), Event->wf[i]) != 0){ 
      std::cout << "Branch not found, stopping script..." << std::endl;
    }
  }

  N = tWaveforms->GetEntries();
  cout << N << endl;
  histos = new Histograms(nchannels);
  histos->OpenHistos();
};

run::run(const char* fnameDRS, const char* fnameTel, int nchannels)
{
  nch = nchannels;
  for (int i=0;i<nch; i++) wfavg.push_back(new waveform());
  Event = new event(nch);
  
  fDRS = new TFile(fnameDRS, "READ");
  fTel = new TFile(fnameTel, "READ");
  
  tWaveforms = (TTree *)fDRS->Get("Plane0/Waveforms");
  tTracks = (TTree *)fTel->Get("Tracks");
  
  for (int i=0; i<nch; i++){
    if(tWaveforms->SetBranchAddress(Form("waveform%d",i), Event->wf[i]) != 0){ 
      std::cout << "Branch not found, stopping script..." << std::endl;
    }
  }
//   if(tWaveforms->SetBranchAddress("", Event->wfAcknowledge) != 0) std::cout << "Branch not found, stopping script..." << std::endl;
        
  if(tTracks->SetBranchAddress("NTracks", &Event->ntracks) != 0) std::cout << "Branch not found, stopping script..." << std::endl;
  if(tTracks->SetBranchAddress("X",       &Event->trk_x) != 0) std::cout << "Branch not found, stopping script..." << std::endl;
  if(tTracks->SetBranchAddress("Y",       &Event->trk_y) != 0) std::cout << "Branch not found, stopping script..." << std::endl;
  if(tTracks->SetBranchAddress("SlopeX",  &Event->trk_slopex) != 0) std::cout << "Branch not found, stopping script..." << std::endl;
  if(tTracks->SetBranchAddress("SlopeY",  &Event->trk_slopey) != 0) std::cout << "Branch not found, stopping script..." << std::endl;
  if(tTracks->SetBranchAddress("Chi2",    &Event->trk_chi2) != 0) std::cout << "Branch not found, stopping script..." << std::endl;
  if(tTracks->SetBranchAddress("Dof",     &Event->trk_dof) != 0) std::cout << "Branch not found, stopping script..." << std::endl;
  
  N = TMath::Min(tTracks->GetEntries(), tWaveforms->GetEntries());
  cout << "Events: " << tTracks->GetEntries() << " (tel), " << tWaveforms->GetEntries() << " (drs) --> " << N << endl; 
  histos = new Histograms(nchannels);
//   histos->Init(nch);
};

run::~run()
{
  for (uint i=0; i< wfavg.size(); i++) delete wfavg[i];
  delete Event;
  fDRS->Close();
  delete fDRS;
  fTel->Close();
  delete fTel;
//   delete histos;
};  
  
void run::GetEvent(int ie)
{
  tTracks->GetEntry(ie);
  tWaveforms->GetEntry(ie);
  for (int it=0; it<Event->ntracks; it++) {
    // Calculate reduced chi2 of the track
    Event->trk_redchi2[it] = Event->trk_chi2[it] / Event->trk_dof[it];
    // Apply slope correction to telescope track
    for (int ich=0; ich<nch; ich++){
      Event->x[ich][it] = Event->trk_x[it] + zDUT[ich] * Event->trk_slopex[it];
      Event->y[ich][it] = Event->trk_y[it] + zDUT[ich] * Event->trk_slopey[it];
    }
  }
  return;
};

void run::GetEventDRS(int ie)
{
  tWaveforms->GetEntry(ie);
  return;
};

void run::GetAvgWaveformAsync(int ch, int ch_veto, int _Nevents)
{
  int Nevents = TMath::Max(_Nevents, (int)(tWaveforms->GetEntries()));
  for(int ipt=0; ipt<DRS_N_POINTS; ipt++) wfavg[ch]->wf[ipt] = 0;     // reset average waveform;
  
  for (int ie=0; ie<Nevents; ie++){
    GetEventDRS(ie);
    int bin_veto = Event->GetVetoPosition(ch_veto, 0.1, +1);          // when the FE was reset
    
    for(int ipt=0; ipt<DRS_N_POINTS; ipt++) 
      wfavg[ch]->wf[ (ipt-bin_veto+DRS_N_POINTS)%DRS_N_POINTS ] += Event->wf[ch]->wf[ipt];
  }
  
  for(int ipt=0; ipt<DRS_N_POINTS; ipt++) wfavg[ch]->wf[ipt] /= Nevents;
  
  return;
};

void run::FindSamplingTime(int ch, int nsamples, int t_min, int t_max, int pulse_polarity, int bin_first, int bin_last)
{
  // Search for the position of signal peak in first nsamples waveforms
  // t_min/t_max = number of bins around peak in avgwf to sample for maximum in individual events
  // bin_first/bin_last = range in which to search for the signal
  // In subsequent analysis the pulse height will be sampled only in the interval [peak-t_min, peak+t_max]
  
//   int first_bin = -t_min+30;
//   int delta_noise=(-t_min+t_max) + 30;     // time before the signal sampling window in which to sample for noise
  int delta_noise = 50 + t_max;     // time after the signal sampling window in which to sample for noise
  for(int ipt=0; ipt<DRS_N_POINTS; ipt++) wfavg[ch]->wf[ipt] = 0;    // reset average waveform;

  for (int ie=0; ie<nsamples; ie++)
  {
    if (N <= ie) break;
    tWaveforms->GetEntry(ie);
//     GetEvent(ie);
//     Event->SubtractPickup(ch, 1, phase_correction[ch]);
    for(int ipt=0; ipt<DRS_N_POINTS; ipt++) wfavg[ch]->wf[ipt] += Event->wf[ch]->wf[ipt];
  }
  for(int ipt=0; ipt<DRS_N_POINTS; ipt++) wfavg[ch]->wf[ipt] /= nsamples;
  
  if(pulse_polarity > 0){
    ts_min[ch] = wfavg[ch]->GetMaxPos(bin_first-t_min, bin_last-t_max)+t_min;   // [bin_first-t_min, bin_last-t_max] to avoid index overflows due to the width of the integrating window
    ts_max[ch] = wfavg[ch]->GetMaxPos(bin_first-t_min, bin_last-t_max)+t_max;
    tn_min[ch] = ts_min[ch]+delta_noise;
    tn_max[ch] = ts_max[ch]+delta_noise;
//     tn_min[ch] = wfavg[ch]->GetMaxPos(bin_first-t_min, bin_last-t_max)+t_min-delta_noise;
//     tn_max[ch] = wfavg[ch]->GetMaxPos(bin_first-t_min, bin_last-t_max)+t_max-delta_noise;
  }
  else{
    ts_min[ch] = wfavg[ch]->GetMinPos(bin_first-t_min, bin_last-t_max)+t_min;
    ts_max[ch] = wfavg[ch]->GetMinPos(bin_first-t_min, bin_last-t_max)+t_max;
    tn_min[ch] = ts_min[ch]+delta_noise;
    tn_max[ch] = ts_max[ch]+delta_noise;
//     tn_min[ch] = wfavg[ch]->GetMinPos(bin_first-t_min, bin_last-t_max)+t_min-delta_noise;
//     tn_max[ch] = wfavg[ch]->GetMinPos(bin_first-t_min, bin_last-t_max)+t_max-delta_noise;
  }
  cout << "Sampling time ch. " << ch << ": " << ts_min[ch] << " " << ts_max[ch] << endl;
  cout << "Noise Sampling time ch. " << ch << ": " << tn_min[ch] << " " << tn_max[ch] << endl;
  return;
}

float run::FindThreshold(int ch, int nsamples, float sigma, int pulse_polarity, int mode)
{
  // Finds discrimination threshold for channel #ch
  // The function samples noise level in first [nsamples] samples and fills it in a histogram
  // The obtained distribution is fitted with a gaussian
  // The threshold is set [sigma] standard deviations above the noise level
  
  TH1I hNoise("hnoise","hnoise",400, -0.04, 0.04);  // histogram for saving noise levels
  for (int ie=0; ie<nsamples; ie++)
  {
    if (N <= ie) break; // if requested event exceeds the number of events in the file
    tWaveforms->GetEntry(ie);
    Event->BaselineCorrection(ch,150,1000);
//     GetEvent(ie);
//     Event->SubtractPickup(ch, 1, phase_correction[ch]);
    float charge = Event->GetCharge(ch, tn_min[ch], tn_max[ch], mode); 
    hNoise.Fill(charge); 
  }
//   TF1 f("fitf","gaus", hNoise.GetXaxis()->GetXmin(), hNoise.GetXaxis()->GetXmax());
//   hNoise.Fit("fitf","RQ");
//   float thr = f.GetParameter(1) + sigma * f.GetParameter(2);  // Get gaussian width and multiply by desired sigma
//   cout << ch << " " << f.GetParameter(1) << " " << f.GetParameter(2) << " " << f.GetParameter(2)*sigma << endl;
  float _thr = hNoise.GetMean() + pulse_polarity*sigma*hNoise.GetRMS();  // Get gaussian width and multiply by desired sigma
  thr[ch] = _thr;
  cout << "Threshold     ch. " << ch << ": " << _thr*1000 << " mV" << endl;
  return _thr;
};

vector<float>* run::GetChargeAll(int mode, int polarity)
{
  for (int ch=0; ch<nch; ch++) {
//     Event->charge[ch] = Event->GetCharge(ch, ts_min[ch], ts_max[ch], mode);
    float _charge = Event->GetCharge(ch, ts_min[ch], ts_max[ch], mode);
    Event->charge[ch] = _charge*polarity;
    if (_charge*polarity > thr[ch]*polarity) Event->hit[ch]=true;
    else Event->hit[ch]=false;
    
//     cout << "ch. " << ch << " (" << polarity << "): " << _charge*polarity << " " << thr[ch]*polarity << "\t" << (Event->hit[ch] ? "PASSED" : "MISS") << endl;
  }
  
  Event->GetCluSizeDRS();  // Calculate cluster size on the fly
  return &Event->charge;
}

vector<float>* run::GetNoiseAll(int mode, int polarity)
{
  for (int ch=0; ch<nch; ch++) Event->noise[ch] = polarity*Event->GetCharge(ch, tn_min[ch], tn_max[ch], mode);
  return &Event->noise;
}

void run::FillHistograms()
{
//   histos->hTracks->Fill(Event->x[c[0], Event->y[ch][0]);  // one track per good event
  for (int ich=0; ich<nch; ich++){
    histos->hC[ich]->Fill(1000*Event->charge[ich]);  
    histos->hN[ich]->Fill(1000*Event->noise[ich]);  
    histos->hTracks[ich]->Fill(Event->x[ich][0], Event->y[ich][0]);  
    if (Event->hit[ich]){
      histos->hPassed[ich]->Fill(Event->x[ich][0], Event->y[ich][0]);  
    }
    histos->h2Charge[ich]->Fill(Event->x[ich][0], Event->y[ich][0], Event->charge[ich]);  
    histos->h2AvgCluSize[ich]->Fill(Event->x[0][0], Event->y[0][0], Event->clu_size[ich]);
  }
  // fill cluster size
  short clu_size = Event->clu_size_drs;
  histos->hCluSize->Fill(clu_size);
  if      (clu_size == 1) histos->h2CluSize1->Fill(Event->x[0][0], Event->y[0][0]);
  else if (clu_size == 2) histos->h2CluSize2->Fill(Event->x[0][0], Event->y[0][0]);  
  else if (clu_size > 2)  histos->h2CluSizeGT2->Fill(Event->x[0][0], Event->y[0][0]);
  
  return;
}

// void run::FillAnalyzedEvent(int itr[2], float charge[NCH_DRS])
// {
//   IOHand->FillAnalyzedEvent(itr, charge, sigw_min);
//   return;
// };
// 

#endif