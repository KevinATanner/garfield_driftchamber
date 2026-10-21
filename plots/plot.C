#include <TApplication.h>
#include <TCanvas.h>

#include "Garfield/MediumMagboltz.hh"
#include "Garfield/ViewMedium.hh"

using namespace Garfield;

void plot() {

  // Setup gas
  MediumMagboltz gas;
  gas.LoadGasFile("ar_90_ch4_10_11142025.gas");
  gas.LoadIonMobility("IonMobility_Ar+_Ar.txt");
  gas.PrintGas();

  ViewMedium view(&gas);

  TCanvas cV("cV", "", 600, 600);
  view.SetCanvas(&cV);
  view.PlotElectronVelocity();
  cV.SaveAs("cV.png");
  cV.Update();
  
  TCanvas cD("cD", "", 600, 600);
  view.SetCanvas(&cD);
  view.PlotElectronDiffusion();
  cD.SaveAs("cD.png");
  cD.Update();
  
  TCanvas cT("cT", "", 600, 600);
  view.SetCanvas(&cT);
  view.PlotElectronTownsend();
  cT.SaveAs("cT.png");
  cT.Update();
  
  TCanvas cA("cA", "", 600, 600);
  view.SetCanvas(&cA);
  view.PlotElectronAttachment();
  cA.SaveAs("cA.png");
  cA.Update();

  TFile* plotroot = new TFile("plotroot.root", "RECREATE");
  cV.Write();
  cD.Write();
  cT.Write();
  cA.Write();

}
