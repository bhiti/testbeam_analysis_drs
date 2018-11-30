#include "TH2.h"
#include "TBox.h"
#include "TF1.h"
#include "TSystem.h"
#include "TBox.h"
#include "TLatex.h"

#include <vector>
#include <iostream>

using namespace std;

#ifndef COMMON_FUNCTIONS
#define COMMON_FUNCTIONS

void  FormatHisto1D(TH1* h, TString hTitle="", float ymin=-123456789, float ymax=-123456789, float labelSize=0.05);
void  FormatHisto2D(TH2* h, TString hTitle="", float zmin=-123456789, float zmax=-123456789, float labelSize=0.05);
void  draw_avg_in_region(TH2 *h, float xc, float xs, float yc, float ys, int color=1, int lineStyle=1);
float avg_in_region(TH2 *h, float xc, float xs, float yc, float ys);
void  print_avg_in_region(TH2F* h, int foldMultiplicity);
void FormatGraph(TGraph *g, TString title, double ymin=-123456789, double ymax=-123456789, int color=1, float labelSize=0.05, float MarkerSize=1, int MarkerType=20);

void FormatHisto1D(TH1* h, TString hTitle, float ymin, float ymax, float labelSize)
{
  h->SetTitle(hTitle.Data());
  h->SetStats(0);
  
  h->GetXaxis()->SetLabelSize(labelSize);
  h->GetXaxis()->SetTitleSize(labelSize);
  h->GetXaxis()->SetTitleOffset(1.1);
//   h->GetXaxis()->SetNdivisions(505);
  
  h->GetYaxis()->SetLabelSize(labelSize);
  h->GetYaxis()->SetTitleSize(labelSize);
  h->GetYaxis()->SetTitleOffset(1.45);
  if (ymin != -123456789 && ymax != -123456789)
    h->GetYaxis()->SetRangeUser(ymin, ymax);
  
  return;
}

void FormatHisto2D(TH2* h, TString hTitle, float zmin, float zmax, float labelSize)
{
  h->SetTitle(hTitle.Data());
  h->SetStats(0);
  
  h->GetXaxis()->SetLabelSize(labelSize);
  h->GetXaxis()->SetTitleSize(labelSize);
  h->GetXaxis()->SetTitleOffset(1.1);
  h->GetXaxis()->SetNdivisions(505);

  h->GetYaxis()->SetLabelSize(labelSize);
  h->GetYaxis()->SetTitleSize(labelSize);
  h->GetYaxis()->SetTitleOffset(1.45);
  
  h->GetZaxis()->SetLabelSize(labelSize);
  h->GetZaxis()->SetTitleSize(labelSize);
  h->GetZaxis()->SetTitleOffset(1);
  
  if (zmin != -123456789 && zmax != -123456789)
    h->GetZaxis()->SetRangeUser(zmin, zmax);
  
  return;
}

void draw_avg_in_region(TH2 *h, float xc, float xs, float yc, float ys, int color, int lineStyle)
{
  TBox* b = new TBox( xc-0.5*xs, yc-0.5*ys, xc+0.5*xs, yc+0.5*ys);
  b->SetFillStyle(0); 
  b->SetLineWidth(4); 
  b->SetLineStyle(lineStyle);
  b->SetLineColor(color);
  b->Draw();
  
  float eff = avg_in_region(h, xc, xs, yc, ys);
  TLatex *tex = new TLatex();
  tex->SetTextFont(42); tex->SetTextSize(0.07);  tex->SetTextColor(color);
  tex->DrawLatex(xc-6, yc-2, Form("%.0lf %%", 100*eff));
  return;
}

float avg_in_region(TH2 *h, float xc, float xs, float yc, float ys)
{
  float avgeff=0;
  int nbins=0;
  for (int ix=1; ix <= h->GetNbinsX(); ix++)
    for (int iy=1; iy <= h->GetNbinsX(); iy++){
      float x = h->GetXaxis()->GetBinCenter(ix);
      float y = h->GetYaxis()->GetBinCenter(iy);
      if ( (x > xc-0.5*xs) && (x < xc+0.5*xs) )
        if ( (y > yc-0.5*ys) && (y < yc+0.5*ys) ){
          nbins++;
          avgeff += h->GetBinContent(ix, iy);
        }
    }
  avgeff /= nbins;
//   cout << h->GetNbinsX() << " " << h->GetNbinsY() << " " << nbins << " " << avgeff << endl;
  return avgeff;
}

void print_avg_in_region(TH2F *h, int foldMultiplicity)
{
  float bcx, bcy; // box center for avg eff calculation
  float bsx, bsy; // box size for avg eff calculation
  
  float pixSizeX = h->GetXaxis()->GetBinUpEdge( h->GetXaxis()->GetLast() ) / foldMultiplicity;
  float pixSizeY = h->GetYaxis()->GetBinUpEdge( h->GetYaxis()->GetLast() ) / foldMultiplicity;
  
  cout << endl << "########\nCorners\n########" << endl;
  for (int i=0; i<4; i++){
    for (int j=0; j<4; j++){
      bcx = pixSizeX*(0.25+0.5*i);
      bcy = pixSizeY*(0.25+0.5*j);
      bsx = pixSizeX/2, bsy = pixSizeY/2;
      
      float avgeff = avg_in_region(h, bcx, bsx, bcy, bsy);
      cout << avgeff << endl;
    }
    cout << endl;
  }
  
  cout << endl << "########\nCenters\n########" << endl;
  for (int i=0; i<2; i++)
    for (int j=0; j<2; j++){
      bcx = pixSizeX*(0.5+i);
      bcy = pixSizeY*(0.5+j);
      bsx = pixSizeX/2, bsy = pixSizeY/2;
      
      float avgeff = avg_in_region(h, bcx, bsx, bcy, bsy);
      cout << avgeff << endl;
    }
  
  cout << endl << "########\nCenter average\n########" << endl;
  {
    bcx = pixSizeX;
    bcy = pixSizeY;
    bsx = pixSizeX, bsy = pixSizeY;
      
    float avgeff = avg_in_region(h, bcx, bsx, bcy, bsy);
    cout << avgeff << endl;
  }
  
  cout << endl << "########\nTotal average (" << foldMultiplicity << " pixels)\n########" << endl;
  {
    bcx = pixSizeX;
    bcy = pixSizeY;
    bsx = 2*pixSizeX, bsy = 2*pixSizeY;
        
    float avgeff = avg_in_region(h, bcx, bsx, bcy, bsy);
    cout << avgeff << endl;
  }
  
  return;
}

void FormatGraph(TGraph *g, TString title, double ymin, double ymax, int color, float labelSize, float MarkerSize, int MarkerType)
{
  if (ymin != -123456789 && ymax != -123456789)
    g->GetYaxis()->SetRangeUser(ymin, ymax);
  g->SetTitle(title.Data());
  g->SetMarkerSize(MarkerSize);
  g->SetMarkerStyle(MarkerType);
  g->SetMarkerColor(color);
  g->SetLineColor(color);
  g->SetFillColor(0);
  g->GetXaxis()->SetLabelSize(labelSize);
  g->GetXaxis()->SetTitleSize(labelSize);
  g->GetYaxis()->SetLabelSize(labelSize);
  g->GetYaxis()->SetTitleSize(labelSize);
  g->GetYaxis()->SetTitleOffset(1.3);
  g->SetLineWidth(2);
  return;
}

#endif
