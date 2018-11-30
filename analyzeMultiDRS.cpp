#include "TROOT.h"
#include "TFile.h"
#include "TPad.h"
#include "TString.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLine.h"

#include "src/event.cpp"
#include "src/run.cpp"
#include "src/histograms.cpp"

#include <vector>
#include <iostream>

int analyzeDRS(int first, int last=-1, TString inpath="../data/")
{
  gStyle->SetPalette(55);
  gROOT->ProcessLine( "gErrorIgnoreLevel = 1001;");  // ignore info messages
  
  /// Definitions
  int nch=2;
  int charge_mode=0; // 0..integral, 1..step, 2..fit step, 3..avg step
  
  float thr_sigma = 4;
  float cut_chi2 = 3;         // reduced chi2
  float cut_time_min = -475;   // constraint on time of arrival of the acknowledgement signal
  float cut_time_max = 1111;
  
  int ts_min=-20;      // integration time around the sampling point
  int ts_max= 20;
  int polarity = -1;  // transient pulse polarity
  int search_window_min=0;    // which bins to search for the signal
  int search_window_max=250;
  
  int nbinsx=50, nbinsy=nbinsx;
//   float xfidmin=-7500, xfidmax=-2000;
//   float yfidmin=-3500, yfidmax=2000;
  float xfidmin, xfidmax, yfidmin, yfidmax;
  if (first == 1149){
//     xfidmin=-5900, xfidmax=-4600;
//     yfidmin=-3400, yfidmax=-2000;
    xfidmin=-4500, xfidmax=-2200;
    yfidmin=  500, yfidmax= 1400;
    nbinsx=20, nbinsy=nbinsx;
  }
  else if (first == 1151){
    xfidmin=-7500, xfidmax=-2000;
    yfidmin=-3500, yfidmax= 1500;
  }
  else if (first == 1154){
    xfidmin=-5000, xfidmax=-2000;
    yfidmin=-3500, yfidmax= 1500;
  }
  else if (first == 1156){
    xfidmin=-5500, xfidmax=-4000;
    yfidmin=-1000, yfidmax= 1500;
  }
  else{
    xfidmin=-7500, xfidmax=-2000;
    yfidmin=-3500, yfidmax=2000;
  }
  
  float zPosDUT[] = {110000, 110000, 110000, 110000};
  
  //////////////////////////////////////////////
  float thr[4];
  
  if (last < 0) last = first;
  TString outFileName = inpath + (TString)(Form("results/res_DRS_%03d-%03d.root", first, last));
  TFile *outFile = new TFile(outFileName, "RECREATE");
  
  event* Event;
  Histograms *histos = new Histograms(nch);
  for (int i=0; i<nch; i++) histos->SetHistosRange(i, xfidmin, xfidmax, yfidmin, yfidmax, nbinsx, nbinsy);
  histos->OpenHistos();

  TH1I *hRatio = new TH1I("hRatio", "hRatio", 100, -10, 90);
  for (int irun = first; irun <= last; irun++)
  {
    TString fnameDRS = inpath + (TString)(Form("conv/run%03d.root", irun));
    TString fnameTel = inpath + (TString)(Form("tel/cosmic_%06d_track-data.root", irun));
    cout << endl << "####################" << endl << "Analyzing run " << irun << endl << "###################" << endl << endl;
    run *r = new run(fnameDRS.Data(), fnameTel.Data(), nch);
    r->histos = histos;
    Event = r->Event;
    for (int ich=0; ich<nch; ich++) r->zDUT[ich]=zPosDUT[ich];
    // relative phase between waveforms (in sampling bins) for fine tuning pickup subtraction
    r->phase_correction[0] = 0;
    r->phase_correction[1] = 0;
//     r->phase_correction[2] = 0;
  
    outFile->cd();    
//     Event->wfAcknowledge = Event->wf[3];
    
//     cout << "Discriminating thresholds :" << endl;
    for (int i=0; i<nch; i++) {
      r->FindSamplingTime(i, 10000, ts_min, ts_max, polarity, search_window_min, search_window_max);
      r->FindThreshold(i, 10000, thr_sigma, polarity, charge_mode);
      thr[i] = r->thr[i];
      histos->SaveThresholdValue(i, polarity*1000*thr[i]);
    }
    
    int nevents = TMath::Min(r->N, 1111111);
    cout << "Analyzing " << nevents << " events" << endl;
    for (int ie=0; ie<nevents; ie++)
    {
      if (ie % 10000 == 0) cout << "\r" << ie << flush;
      r->GetEvent(ie);
      
      if (Event->ntracks != 1) continue;
      
      histos->hChi2->Fill(Event->trk_redchi2[0]);
      if (Event->trk_redchi2[0] > cut_chi2) continue;
      
//       int TOA = Event->GetTOA(3, 0.2);   // time of arrival, check required due to a telescope "bug"
//       histos->hTOA->Fill(TOA);
//       cout << TOA << endl;
//       if (TOA < cut_time_min || TOA > cut_time_max){
//         cout << "Event " << ie << " time of arrival out of cut window. Skipping." << endl;
//         continue;
//       }
      
//       TGraph* g0 = Event->wf[0]->Draw();
//       g0->SetLineColor(4);
//       TGraph* g1 = Event->wf[1]->Draw();
//       g1->SetLineColor(2);
//       
      for (int ich=0; ich<nch; ich++) Event->BaselineCorrection(ich, 150, 1000);
      
      // waveform subtraction: due to high pickup subtract waveforms to improve S/N.
      // Below define the active area of ch. 0. 
      // If track passes through subtract w0-w1. 
      // If the track does not pass through subtract w1-w0
      float ch0_xmin, ch0_xmax, ch0_ymin, ch0_ymax;
      if (irun == 1149) {
        ch0_xmin = -111111;
        ch0_xmax = 111111;
        ch0_ymin = -111111;
        ch0_ymax = -1000;
      }
      else if (irun == 1151 || irun == 1152){
        ch0_xmin = -6800;
        ch0_xmax = -4400;
        ch0_ymin = -3200;
        ch0_ymax = -400;
      }
      else if (irun == 1154 || irun == 1155){
        ch0_xmin = -4520;
        ch0_xmax = -2000;
        ch0_ymin = -3500;
        ch0_ymax =   650;
      }
      else {
        ch0_xmin = -6800;
        ch0_xmax = -4400;
        ch0_ymin = -3200;
        ch0_ymax = -400;
      }

      float x=Event->trk_x[0], y=Event->trk_y[0];
      bool trk_in_ch0 = x > ch0_xmin && x < ch0_xmax && y > ch0_ymin && y < ch0_ymax;
//       if (y < -1000)  Event->SubtractPickup(0,1, r->phase_correction[0]);
//       else            Event->SubtractPickup(1,0, r->phase_correction[1]);      
//       if (trk_in_ch0) Event->SubtractPickup(0,1, r->phase_correction[0]);
//       else            Event->SubtractPickup(1,0, r->phase_correction[1]);
      
//       if (trk_in_ch0) hRatio->Fill(-5);
//       else hRatio->Fill(50);
      
//       TGraph* gbl0 = Event->wf[0]->Draw();
//       gbl0->SetLineColor(6);
//       gbl0->SetMarkerSize(0.5);
//       gbl0->SetMarkerStyle(20);
//       gbl0->SetMarkerColor(6);
//       TGraph* gbl1 = Event->wf[1]->Draw();
// //       gbl1->SetLineColor(4);
//       Event->ShiftPhase(1,r->phase_correction[1]);
//       TGraph* gsp1 = Event->wf[1]->Draw();
//       gsp1->SetLineColor(1);
//       gsp1->SetMarkerSize(0.5);
//       gsp1->SetMarkerStyle(20);
//       gsp1->SetMarkerColor(1);
//       
      r->GetChargeAll(charge_mode, polarity);
      r->GetNoiseAll(charge_mode);
      r->FillHistograms();
// 
//       if (Event->hit[0]){
// // //         TCanvas cc("c","c",1);
// // //         cc.Divide(2,1);
// // //         cc.cd(1);
//         TGraph* g = Event->wf[0]->Draw();
//         g->SetTitle("");
//         g->GetXaxis()->SetRangeUser(search_window_min, search_window_max);
//         gsp1->GetXaxis()->SetRangeUser(40,120);
//         g->GetYaxis()->SetRangeUser(-0.1,0.05);
// //         g->GetYaxis()->SetRangeUser(-0.02,0.02);
//         g->SetLineColor(2);
// //         gbl1->Draw("APL");
//         gsp1->Draw("APL");
// //         Event->wf[1]->Draw()->Draw("PL");
// //         g1->Draw("PL");
//         gbl0->Draw("PL");
//         g->Draw("PL");
// 
// //         if (Event->
// //         Event->SubtractPickup(0,1);
//         TGraph* gg = Event->wf[2]->Draw();
//         gg->SetLineColor(4);
// //         gg->Draw("PL");
//         
// //         Event->SubtractPickup(0,1,r->phase_correction[0]);
//         
// // //         cc.cd(2);
// // //         TH2I h("h","h",nbinsx, xfidmin, xfidmax, nbinsy, yfidmin, yfidmax);
// // //         h.Fill(Event->trk_x[0], Event->trk_y[0]);
// // //         h.Draw("COLZ");
// // //   //       Event->wf[3]->Draw()->Draw("APL");
//         gPad->Update();
//         gPad->Print(Form("../plots/%03d-%03d/run%d_ev%d_charge%.3lf.png", first, last, irun, ie, Event->charge[0]));
//         cout << ie << " " << Event->charge[0] << " " << Event->charge[1] << endl;
//         char c;
//         cin >> c;
//       }
      
//       TGraph* g = Event->wf[0]->Draw();
//       g->GetXaxis()->SetRangeUser(search_window_min, search_window_max+100);
//       g->GetYaxis()->SetRangeUser(-0.05,0.05);
//       g->Draw("APL");
// //       Event->wf[0]->Draw()->Draw("APL");
// //       Event->wf[3]->Draw()->Draw("APL");
//       gPad->Update();
//       char c;
//       cin >> c;
    }
    cout << "\r" << nevents << endl;
    
    TCanvas *c2 = new TCanvas("c2", "c2", 1);    
    TGraph* g = r->wfavg[0]->Draw();
    FormatGraph(g, (TString)("Average waveform ch. 1"));
    g->GetXaxis()->SetRangeUser(search_window_min, search_window_max);
    g->Draw("AL");
    
//     r->wfavg[0]->Draw()->Draw("APL");
//     r->wfavg[1]->Draw()->Draw("PL");
//     cout << r->wfavg[0]->GetMinPos(search_window_min,search_window_max) << endl;;
//     cout << r->wfavg[1]->GetMinPos(search_window_min,search_window_max) << endl;;
    gPad->SetMargin(0.15,0.18,0.15,0.1);
    gPad->Update();
    gPad->Print(Form("../plots/%03d-%03d/avgwf%d.png", first, last, irun));
    
//     cout << "Spectrum integral: " << histos->hC[0]->Integral(1, histos->hC[0]->GetNbinsX()) << " (" << 1 << ", " << histos->hC[0]->GetNbinsX() << "): ";
//     cout << histos->hC[0]->Integral(histos->hC[0]->FindBin(r->thr[0]*polarity*1000), histos->hC[0]->GetNbinsX()) << " (" << histos->hC[0]->FindBin(r->thr[0]*polarity*1000) << ", " << histos->hC[0]->GetNbinsX() << ")" << endl;;
  
    delete r;
  }
  
  histos->FinalProcessing();
  histos->Plotting(Form("../plots/%03d-%03d_subtr/", first, last));
  
  TCanvas *c1 = new TCanvas("c1", "c1", 1600, 1200);
  c1->Divide(4,2);
  int iPad=1;
  c1->cd(iPad++);
  histos->hTracks[1]->Draw("COLZ");
  c1->cd(iPad++);
  histos->hPassed[1]->Draw("COLZ");
  c1->cd(iPad++);
  histos->hEff[0]->Draw("COLZ");
  c1->cd(iPad++);
  histos->hEff[1]->Draw("COLZ");
  c1->cd(iPad++);
//   histos->hChi2->Draw();
  histos->hN[0]->Draw();
  c1->cd(iPad++);
//     gPad->SetLogy();
  histos->hN[1]->Draw();
  c1->cd(iPad++);
//   gPad->SetLogy();
  histos->hC[0]->Draw();
  TLine l;
  l.DrawLine(thr[0]*polarity*1000,1,thr[0]*polarity*1000, histos->hC[0]->GetMaximum());
  c1->cd(iPad++);
//   gPad->SetLogy();
  histos->hC[1]->Draw();
  l.DrawLine(thr[1]*polarity*1000,1,thr[1]*polarity*1000, histos->hC[1]->GetMaximum());
  
//   c1->cd(0);
//   gPad->SetLogy();
//   histos->hC[0]->Draw();
//   l.DrawLine(thr[0]*polarity*1000,1,thr[0]*polarity*1000, histos->hC[0]->GetMaximum());
//   gPad->Print(Form("../plots/%03d-%03d/hC0_log.png", first, last));
//   histos->hC[1]->Draw();
//   l.DrawLine(thr[1]*polarity*1000,1,thr[1]*polarity*1000, histos->hC[1]->GetMaximum());
//   gPad->Print(Form("../plots/%03d-%03d/hC1_log.png", first, last));
// //   outFile->Write();
  outFile->Write();
    
  cout << avg_in_region(histos->hEff[0], -3500, 900, -700, 1000) << endl;
  cout << avg_in_region(histos->hEff[1], -4670,  80,  450,  900) << endl;

  return 0;
}