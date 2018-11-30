#include "TROOT.h"
#include "TFile.h"
#include "TPad.h"
#include "TString.h"
#include "TCanvas.h"
#include "TStyle.h"

// #include "src/waveform.cpp"
#include "src/event.cpp"
#include "src/run.cpp"
#include "src/histograms.cpp"

#include <vector>
#include <iostream>

// #ifndef DRS_N_POINTS
// #define DRS_N_POINTS 1024
// #endif


// int analyzeDRS(const char* fnameDRS="../data/conv/run087.root", const char* fnameTel="../data/conv/cosmic_004488_track-data.root")
int analyzeDRS(int first, int last=-1, TString inpath="../data/drs_conv/")
{
  gStyle->SetPalette(55);
  /// Definitions
  int nch=2;
  int charge_mode=0; // 0..integral, 1..step, 2..fit step, 3..avg step
  
  float thr_sigma = 3;
  float cut_chi2 = 3;         // reduced chi2
  float cut_time_min = 475;   // constraint on time of arrival of the acknowledgement signal
  float cut_time_max = 485;
  
  int ts_min=-20;      // integration time around the sampling point
  int ts_max=20;
  int polarity = -1;  // transient pulse polarity
    
  float xfidmin=-10000, xfidmax=10000;
  float yfidmin=-10000, yfidmax=10000;
  int nbinsx=100, nbinsy=nbinsx;
  
  float zPosDUT[] = {110000, 110000};
  
  //////////////////////////////////////////////
  if (last < 0) last = first;
  TString outFileName = inpath + (TString)(Form("res_DRS_%03d-%03d.root", first, last));
  TFile *outFile = new TFile(outFileName, "RECREATE");
  
  int nevents;
  event* Event;
  Histograms *histos = new Histograms(nch);
  for (int i=0; i<nch; i++) histos->SetHistosRange(i, xfidmin, xfidmax, yfidmin, yfidmax, nbinsx, nbinsy);
  histos->OpenHistos();

  for (int irun = first; irun <= last; irun++)
  {
    TString fnameDRS = inpath + (TString)(Form("run%03d.root", irun));
//     TString fnameTel = inpath + (TString)(Form("cosmic_%06d_track-data.root", irun));
    cout << endl << "####################" << endl << "Analyzing run " << irun << endl << "###################" << endl << endl;
//     run *r = new run(fnameDRS.Data(), fnameTel.Data(), nch);
    run *r = new run(fnameDRS.Data(), nch);
    r->histos = histos;
//     for (int ich=0; ich<nch; ich++) r->zDUT[ich]=zPosDUT[ich];
    Event = r->Event;
    outFile->cd();    
    
    cout << "Discriminating thresholds :" << endl;
    for (int i=0; i<nch; i++) {
      r->FindSamplingTime(i, 1000, ts_min, ts_max, polarity);
      r->FindThreshold(i, 1000, thr_sigma*polarity);
//       cout << r->thr[i]*1000 << " mV / ";
    }
    r->ts_min[0]=111+ts_min;
    r->ts_max[0]=111+ts_max;
    r->tn_min[0]=211+ts_min;
    r->tn_max[0]=211+ts_max;
    r->ts_min[1]=111+ts_min;
    r->ts_max[1]=111+ts_max;
    r->tn_min[1]=211+ts_min;
    r->tn_max[1]=211+ts_max;
    cout << endl;
    
    nevents = TMath::Min(r->N, 11111111);
    cout << "Analyzing " << nevents << " events" << endl;
    for (int ie=0; ie<nevents; ie++)
    {
      if (ie % 10000 == 0) cout << "\r" << ie << flush;
      r->tWaveforms->GetEntry(ie);
  
//       if (Event->ntracks != 1) continue;
      
//       histos->hChi2->Fill(Event->trk_redchi2[0]);
//       if (Event->trk_redchi2[0] > cut_chi2) continue;
      
//       int TOA = Event->GetTOA(0.2);   // time of arrival, check required due to a telescope "bug"
//       histos->hTOA->Fill(TOA);
//       if (TOA < cut_time_min || TOA > cut_time_max){
//         cout << "Event " << ie << " time of arrival out of cut window. Skipping." << endl;
//         continue;
//       }
            
      r->GetChargeAll(charge_mode, polarity);
      r->GetNoiseAll(charge_mode);
      r->FillHistograms();

      
      Event->wf[1]->Draw()->Draw("APL");
      Event->wf[0]->Draw()->Draw("PL");
      gPad->Update();
      char c;
      cin >> c;
    }
    cout << "\r" << nevents << endl;
    
    TCanvas *c2 = new TCanvas("c2", "c2", 800, 800);
    r->wfavg[0]->Draw()->Draw("APL");
    r->wfavg[1]->Draw()->Draw("PL");
    gPad->Update();
//     delete r;
  }
    
  //   hTOA->Draw();
//   histos->hChi2->Draw();
  histos->FinalProcessing();

  TCanvas *c1 = new TCanvas("c1", "c1", 800, 800);
  c1->Divide(4,1);
  int iPad=1;
//   c1->cd(iPad++);
//   histos->hTracks[0]->Draw("COLZ");
//   c1->cd(iPad++);
//   histos->hPassed[0]->Draw("COLZ");
//   c1->cd(iPad++);
//   histos->hEff[0]->Draw("COLZ");
//   c1->cd(iPad++);
//   histos->h2Charge[0]->Draw("COLZ");
//   c1->cd(iPad++);
//   histos->hChi2->Draw();
//   c1->cd(iPad++);
//   histos->hTOA->Draw();
  c1->cd(iPad++);
  histos->hC[0]->Draw();
  c1->cd(iPad++);
  histos->hN[0]->Draw();
  c1->cd(iPad++);
  histos->hC[1]->Draw();
  c1->cd(iPad++);
  histos->hN[1]->Draw();
  
  int minbin_eff = histos->hC[1]->FindBin(+0.007);
  int minbin_all = histos->hC[1]->FindBin(-0.01);
  int maxbin = histos->hC[1]->FindBin(0.99);
  cout << "Efficiency: " << histos->hC[1]->Integral(minbin_eff, maxbin) / histos->hC[1]->Integral(minbin_all, maxbin) << endl;
  cout << "S/N: " << (0.021-histos->hN[1]->GetRMS())/histos->hN[1]->GetRMS() << endl;
  
//   outFile->Write();
  outFile->Write();
  
  return 0;
}