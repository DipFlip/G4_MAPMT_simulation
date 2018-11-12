// Root macro to analyse the output of the Geant4 simulation SoNDe
// J.R.M. Annand  20th August 2018
// J.R.M. Annand  19th October 2018 merge sub-versions
//
#include "TStyle.h"
#include "TChain.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"

using namespace std;
// MAPMT pixel ID
Int_t idx, idy, id;
// MAPMT pixel dimensions
Double_t pix[] = {6.25, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.25};
Double_t Pixd[9];
Double_t qeMax = 0.40;             // maximum quantum efficiency of cathode
TH2D *PCqeff;                      // grid pattern mimics "dead" regions at edge of pixels
Int_t nPEpixMin = 3;               // minimum # PE to fire a pixel
void SetPixels();                  // setup pixel boundaries
void GetPixel(Double_t, Double_t); // get pixel that coord(x,y) resides in
void GenPixEff();                  // generate grid pattern
//void PlotSonnig(const Char_t* file, Double_t Ethr = 0.02, Double_t res = 0.0 )
void PlotSonnig1(const Char_t *file, Int_t imin = 0, Int_t imax = 1000000000)
{
  TChain *PD = new TChain("PD"); // in case there are several files
  TRandom3 *R3 = new TRandom3();
  PD->Add(file);
  const Int_t maxhits = 4096; // max number hits in event
  const Int_t nPD = 64;       // number pixels
  const Int_t nPix = 8;       // 8x8 grid

  Int_t Nhits; // actual # hits in event
  // arrays to hold info held in tree branches
  Int_t PD_vid[maxhits];
  Float_t PD_t[maxhits];
  Float_t PD_x[maxhits];
  Float_t PD_y[maxhits];
  Float_t PD_z[maxhits];
  Float_t PD_vx[maxhits];
  Float_t PD_vy[maxhits];
  Float_t PD_vz[maxhits];
  Float_t PD_px[maxhits];
  Float_t PD_py[maxhits];
  Float_t PD_pz[maxhits];
  Float_t PD_E[maxhits];
  Float_t PD_Ed[maxhits];

  // Set branch addresses.
  PD->SetBranchAddress("PD_Nhits", &Nhits);
  //PD->SetBranchAddress("PD_idx",PD_idx);
  //PD->SetBranchAddress("PD_idy",PD_idy);
  PD->SetBranchAddress("PD_vid", PD_vid);
  PD->SetBranchAddress("PD_t", PD_t);
  PD->SetBranchAddress("PD_E", PD_E);
  PD->SetBranchAddress("PD_Ed", PD_Ed);
  PD->SetBranchAddress("PD_x", PD_x);
  PD->SetBranchAddress("PD_y", PD_y);
  PD->SetBranchAddress("PD_z", PD_z);
  PD->SetBranchAddress("PD_vx", PD_vx);
  PD->SetBranchAddress("PD_vy", PD_vy);
  PD->SetBranchAddress("PD_vz", PD_vz);
  PD->SetBranchAddress("PD_px", PD_px);
  PD->SetBranchAddress("PD_py", PD_py);
  PD->SetBranchAddress("PD_pz", PD_pz);
  //
  // Create histograms
  TH1D *hidx = new TH1D("hidx", "Hit Pixel X", nPix, 0, nPix);
  TH1D *hidy = new TH1D("hidy", "Hit Pixel Y", nPix, 0, nPix);
  TH2D *h2id = new TH2D("h2id", "Hit Pixel X-Y", nPix, 0, nPix, nPix, 0, nPix);
  TH1D *ht1 = new TH1D("ht1", "Time1", 1000, 0, 1000);
  TH1D *ht2 = new TH1D("ht2", "Time2", 1000, 0, 1000);
  TH1D *hx = new TH1D("hx", "Hit X", 500, -30, 30);
  TH1D *hy = new TH1D("hy", "Hit Y", 500, -30, 30);
  TH1D *hz = new TH1D("hz", "Hit Z", 500, -1, 1);
  TH1D *hxpe = new TH1D("hxpe", "PE X", 500, -30, 30);
  TH1D *hype = new TH1D("hype", "PE Y", 500, -30, 30);
  TH1D *hvx = new TH1D("hvx", "Hit vX", 500, -30, 30);
  TH1D *hvy = new TH1D("hvy", "Hit vY", 500, -30, 30);
  TH1D *hvz = new TH1D("hvz", "Hit vZ", 500, -1, 1);
  TH1D *hdx = new TH1D("hdx", "Hit dX", 500, -30, 30);
  TH1D *hdy = new TH1D("hdy", "Hit dY", 500, -30, 30);
  TH1D *hdxpe = new TH1D("hdxpe", "PE dX", 500, -30, 30);
  TH1D *hdype = new TH1D("hdype", "PE dY", 500, -30, 30);
  TH1D *hed = new TH1D("hed", "Accum-E", 10000, 0, 5.0);
  TH1D *hlamda = new TH1D("hlamda", "Wavelength (nm)", 1000, 0, 1000);
  TH2D *h2xy = new TH2D("h2xy", "PC X-Y", 200, -30, 30, 100, -30, 30);
  TH2D *h2dxy = new TH2D("h2dxy", "PC dX-dY", 200, -20, 20, 200, -20, 20);
  TH2D *h2xype = new TH2D("h2xype", "PE X-Y", 200, -30, 30, 100, -30, 30);
  TH2D *h2dxype = new TH2D("h2dxype", "PE dX-dY", 200, -20, 20, 200, -20, 20);
  TH2D *h2dxvx = new TH2D("h2dxvx", "Hit dX-vX", 200, -20, 20, 200, -25, 25);
  TH2D *h2dyvy = new TH2D("h2dyvy", "Hit dY-vY", 200, -20, 20, 200, -25, 25);
  TH1D *hth = new TH1D("hth", "Hit Theta", 90, 0, 90);
  TH1D *hph = new TH1D("hph", "Hit Phi", 180, -180, 180);
  TH2D *h2ang = new TH2D("h2ang", "Hit Theta-Phi", 90, 0, 90, 180, -180, 180);
  TH1D *hNhits = new TH1D("Nhits", "Nhits", maxhits, 0, maxhits);
  TH1D *hNphot = new TH1D("Nphot", "Nphot", 1024, 0, 1024);
  TH1D *hNpe = new TH1D("Npe", "Npe", 1024, 0, 1024);
  TH1D *hNgs20 = new TH1D("Ngs20", "Ngs20", 16, 0, 16);
  TH1D *hNpix = new TH1D("Npix", "Npix", 64, 0, 64);
  TH1D *hNphPix[nPD];
  TH1D *hNpePix[nPD];
  TH2D *h2vxNpe = new TH2D("h2vxNpe", "Npe-vx", 200, -25, 25, 200, 0, 200);
  Char_t pname[64];
  for (Int_t i = 0; i < nPD; i++)
  {
    sprintf(pname, "Nph%d", i);
    hNphPix[i] = new TH1D(pname, pname, 256, 0, 256);
    sprintf(pname, "Npe%d", i);
    hNpePix[i] = new TH1D(pname, pname, 256, 0, 256);
  }
  //
  // Set up MAPMT pixel dimensions
  SetPixels();
  // # entries in ROOT tree
  Int_t nentries = PD->GetEntries();
  Int_t vid;
  Int_t ixy;
  Double_t x, y, z, vx, vy, vz, dx, dy, px, py, pz, theta, phi, tm, qeff, lamda;
  Double_t vx0, vy0, vz0;
  TVector3 p;
  const Double_t r2d = TMath::RadToDeg();
  Int_t nPCpix[nPD];
  Double_t nPEpix[nPD];
  // loop over entries in tree
  for (Int_t i = 0; i < nentries; i++)
  {
    PD->GetEntry(i); // get # hits and check for overflow
    if (i < imin)
      continue;
    if (i > imax)
      continue;
    hNhits->Fill(Nhits);
    if (Nhits > maxhits)
    {
      printf("%d hits > array dimension %d\n", Nhits, maxhits);
      Nhits = maxhits;
    }
    if (i % 1000 == 0)
      printf("Event: %d\n", i);
    // zero energy, photon and photo-electron counters
    Double_t ed = 0.0;
    Int_t nGS20 = 0;
    Int_t nPC = 0;
    Double_t nPE = 0;
    Int_t Npix = 0;
    for (Int_t k = 0; k < nPD; k++)
    {
      nPCpix[k] = 0;
      nPEpix[k] = 0;
    }
    // loop round all "hits" in one event
    for (Int_t j = 0; j < Nhits; j++)
    {
      vid = PD_vid[j]; // vid is sensitive volume ID
      x = PD_x[j];     // hit coords
      y = PD_y[j];
      z = PD_z[j];
      tm = PD_t[j]; // primary interaction time
      if (vid == 1)
      {                 // hit in GS20, ie neutron/proton/alpha/gamma
        ed += PD_Ed[j]; // accumulate energy loss
        if (nGS20 == 0)
        {                 // 1st step?
          vx0 = PD_vx[j]; // primary interaction point
          vy0 = PD_vy[j];
          vz0 = PD_vz[j];
          vx = x;
          vy = y;
          vz = z;
          hvx->Fill(vx);
          hvy->Fill(vy);
          hvz->Fill(vz);
          ht1->Fill(tm);
        }
        nGS20++; // accumulate # hadron/gamma steps
      }
      else if (vid == 2)
      {                 // hit on photocathode, ie optical photon
        GetPixel(x, y); // which pixel was hit
        // wavelength of scintillation photon
        lamda = 1239.8 / (PD_Ed[j] * 1.0E06);
        hlamda->Fill(lamda);
        // Get effective quantum efficiency
        ixy = PCqeff->FindBin(x, y);
        qeff = PCqeff->GetBinContent(ixy);
        //
        dx = x - vx; // offset of cathode hist from primary hit in GS20
        dy = y - vy;
        px = PD_px[j]; // 3 momentum of photon
        py = PD_py[j];
        pz = PD_pz[j];
        p.SetXYZ(px, py, pz);
        theta = p.Theta() * r2d; // photon angle when hits cathode
        phi = p.Phi() * r2d;
        hidx->Fill(idx); // record pixel hit
        hidy->Fill(idy);
        h2id->Fill(idx, idy);
        hdx->Fill(dx); // record position offset
        hdy->Fill(dy);
        h2dxy->Fill(dx, dy);
        hdxpe->Fill(dx, qeff); // pe histograms are photo electron
        hdype->Fill(dy, qeff);
        h2dxype->Fill(dx, dy, qeff);
        h2dxvx->Fill(dx, vx);
        h2dyvy->Fill(dy, vy);
        hx->Fill(x);
        hy->Fill(y);
        h2xy->Fill(x, y);
        hxpe->Fill(x, qeff);
        hype->Fill(y, qeff);
        h2xype->Fill(x, y, qeff);
        //
        hth->Fill(theta);
        hph->Fill(phi);
        h2ang->Fill(theta, phi);
        nPC++;
        nPCpix[id]++;
        // arrival time of photon at cathode
        ht2->Fill(tm);
        // compute photo-electron numbers
        nPE += qeff;
        nPEpix[id] += qeff;
      }
    }
    hed->Fill(ed);
    hNphot->Fill(nPC);
    hNpe->Fill(nPE);
    h2vxNpe->Fill(vx, nPE);
    for (Int_t k = 0; k < nPD; k++)
    {
      if (nPCpix[k] > 0)
        hNphPix[k]->Fill(nPCpix[k]);
      if (nPEpix[k] > 0)
      {
        hNpePix[k]->Fill(nPEpix[k]);
        if (nPEpix[k] > nPEpixMin) // over pixel PE threshold?
          Npix++;
      }
    }
    hNpix->Fill(Npix);
    hNgs20->Fill(nGS20);
  }
  Int_t jx = 1;
  printf("Pixel detected-photon totals\n");
  for (Int_t ix = 1; ix <= nPix; ix++)
  {
    printf("%6.0f %6.0f %6.0f %6.0f %6.0f %6.0f %6.0f %6.0f\n",
           h2id->GetBinContent(jx, ix), h2id->GetBinContent(jx + 1, ix),
           h2id->GetBinContent(jx + 2, ix), h2id->GetBinContent(jx + 3, ix),
           h2id->GetBinContent(jx + 4, ix), h2id->GetBinContent(jx + 5, ix),
           h2id->GetBinContent(jx + 6, ix), h2id->GetBinContent(jx + 7, ix));
  }
  return;
}

void SetPixels()
{
  // Setup MAPMT pixel dimensions
  Double_t pcd = 0.0;
  for (Int_t i = 0; i < 8; i++)
  {
    pcd = pcd + pix[i];
  }
  pcd *= 0.5;
  Pixd[0] = -pcd;
  for (Int_t i = 1; i <= 8; i++)
    Pixd[i] = Pixd[i - 1] + pix[i - 1];
  //
  // Setup a quantum-efficiency 2D map
  GenPixEff();
}

void GetPixel(Double_t x, Double_t y)
{
  // Determine the hit pixel of the MAPMT
  for (Int_t i = 0; i < 8; i++)
  {
    if ((x >= Pixd[i]) && (x <= Pixd[i + 1]))
    {
      idx = i;
      break;
    }
  }
  for (Int_t i = 0; i < 8; i++)
  {
    if ((y >= Pixd[i]) && (y <= Pixd[i + 1]))
    {
      idy = i;
      break;
    }
  }
  id = idy * 8 + idx;
}

void GenPixEff()
{
  // Generate a quantum efficiency map for the MAPMT
  //
  Int_t nx = 200;
  Int_t idx, idy;
  Double_t sig = 0.25; // width pixel "dead" boundary
  Double_t dip = 0.1;  // degree of deadness
  Double_t x, y, x0, x1, y0, y1, effx, effy, eff;
  Double_t dX = (Pixd[8] - Pixd[0]) / nx;
  PCqeff = new TH2D("QEcathode", "QEcathode",
                    nx + 1, Pixd[0], Pixd[8], nx + 1, Pixd[0], Pixd[8]);
  for (Int_t i = 0; i <= nx; i++)
  {
    x = Pixd[0] + i * dX;
    for (Int_t ii = 0; ii < 8; ii++)
    {
      if ((x >= Pixd[ii]) && (x <= Pixd[ii + 1]))
      {
        idx = ii;
        break;
      }
    }
    x0 = Pixd[idx];
    x1 = Pixd[idx + 1];
    effx = 1.0 - dip * (TMath::Gaus(x, x0, sig, true) + TMath::Gaus(x, x1, sig, true));
    for (Int_t j = 0; j <= nx; j++)
    {
      y = Pixd[0] + j * dX;
      for (Int_t ii = 0; ii < 8; ii++)
      {
        if ((y >= Pixd[ii]) && (y <= Pixd[ii + 1]))
        {
          idy = ii;
          break;
        }
      }
      y0 = Pixd[idy];
      y1 = Pixd[idy + 1];
      effy = 1.0 - dip * (TMath::Gaus(y, y0, sig, true) + TMath::Gaus(y, y1, sig, true));
      eff = qeMax * effx * effy;
      //printf("%lf %lf %lf %lf %lf\n",x,y,effx,effy,eff);
      PCqeff->SetBinContent(i + 1, j + 1, eff);
    }
  }
}
