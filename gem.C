#include <iostream>
#include <vector>
#include <memory>

#include <TApplication.h>
#include <TCanvas.h>
#include <Garfield/MediumMagboltz.hh>
#include <Garfield/Sensor.hh>
#include <Garfield/ViewField.hh>
#include <Garfield/ComponentAnalyticField.hh>
#include <TH1F.h>
#include "Garfield/ViewSignal.hh"
#include "Garfield/ViewDrift.hh"
#include "Garfield/AvalancheMicroscopic.hh"
#include <TH2D.h>
#include "Garfield/ViewMedium.hh"
#include "Garfield/ViewCell.hh"
#include <TGraph.h>
#include "Garfield/AvalancheMC.hh"
#include "Garfield/AvalancheMicroscopic.hh"
#include "Garfield/TrackHeed.hh"
#include "Garfield/DriftLineRKF.hh"
#include <TH1F.h>
#include <TCanvas.h>
#include <TGraph.h>


using namespace Garfield;

int main(int argc, char *argv[]) {

  TApplication app("app", &argc, argv);

  
  double mm = 1e-3;
  double kilovolt = 1000.0;
  double volt = 1.0;
  
  // geometry paramters
  // units of cm
  
  const int nStrips = 71;
  const double pitch = 0.3; //1.0 * mm spacing;
  const double stripWidth = 0.2;  //2 mm width;
  const double topY = 1.45; //cm 
  const double midY = 0.0;
  const double botY = -1.45;
  const double leftX = -1.0; 
  const double wireAnodeRadius = 0.002; // 20 um 
  const double wireCathodeRadius = 0.01; //100 um
  // top and bottom plane cathode strip potentials

  const double V_first = -0.5 * kilovolt;
  const double V_fourth = 0.0 * volt;
  const double V_last = -4.0 * kilovolt;

  // Middle plane wires potentials:
  const double V_wire_near = -0.5 * kilovolt; // parallel to first strip
  const double V_wire_anode = 1.5 * kilovolt; // parallel to 4th strip (readout)
  const double V_wire_far = -4.0 * kilovolt;  // parallel to last strip


  MediumMagboltz* gas = new MediumMagboltz();
  gas->LoadGasFile("/u/athirak2/Gem/ar_90_ch4_10_11142025.gas");
  gas->EnableDrift();
  gas->LoadIonMobility("/u/athirak2/garfieldpp/install/share/Garfield/Data/IonMobility_Ar+_Ar.txt");  
  
  // Create analytic field component
  ComponentAnalyticField* comp = new ComponentAnalyticField();
  comp->SetMedium(gas);

  // assigning potentials along strips
    std::vector<double> potentials(nStrips, 0.);
    for (int i = 0; i <= 3; ++i) {
        double t = double(i) / 3.0;
        potentials[i] = V_first + t * (V_fourth - V_first);
    }
    for (int i = 3; i < nStrips; ++i) {
        double t = double(i - 3) / double(nStrips - 1 - 3);
        potentials[i] = V_fourth + t * (V_last - V_fourth);
    }

    std::cout << "Strip potentials:" << std::endl;
    for (int i = 0; i < nStrips; ++i) {
      std::cout << "Strip " << i << ": " << potentials[i] << " V" << std::endl;
    }


    for (int i = 0; i < nStrips; ++i) {
    double cx = leftX + i * pitch;
    comp->AddWire(cx, topY, stripWidth/2, potentials[i], "topStrip" + std::to_string(i));
    }
    for (int i = 0; i < nStrips; ++i) {
    double cx = leftX + i * pitch;
    comp->AddWire(cx, botY, stripWidth/2, potentials[i], "botStrip" + std::to_string(i));
    }

    
    // middle plane wires
    
    const double wire1_x = leftX;
    const double wire2_x = leftX + 3 * pitch;
    const double wire3_x = leftX + (nStrips - 1) * pitch;

    comp->AddWire(wire1_x, midY, wireCathodeRadius, V_wire_near, "nearCathode");
    comp->AddWire(wire2_x, midY, wireAnodeRadius, V_wire_anode, "Anode");
    comp->AddWire(wire3_x, midY, wireCathodeRadius, V_wire_far, "farCathode");

    comp->AddReadout("Anode");

    //std::cout<<" distance from anode to far cathode = "<<wire3_x - wire2_x<<std::endl;
    
    //TCanvas* c1 = new TCanvas("c1", "Drift Chamber Geometry", 800, 600);

    //comp->PlotCell(c1);

   ViewField* fieldView = new ViewField();
  fieldView->SetComponent(comp);
  fieldView->SetNumberOfContours(50);
  double xminPlot = -2.0; //  leftX - 2 * pitch;
  double xmaxPlot = 20.0 ; //leftX + (nStrips - 1) * pitch + 2 * pitch;
  double yminPlot = botY - 15.;
  double ymaxPlot = topY + 15.;

  fieldView->SetArea(xminPlot, yminPlot, xmaxPlot, ymaxPlot);
  //fieldView->PlotContour("V");
  

    Sensor* sensor = new Sensor();
    sensor->AddComponent(comp);
    sensor->AddElectrode(comp, "Anode");
    const double tstep = 0.01;
    const double tmin = -0.5 * tstep;
    const unsigned int nbins = 1000000;
    sensor->SetTimeWindow(tmin, tstep, nbins);

    // drift lines from a track
    TrackHeed track(sensor);
    track.SetParticle("muon");
    track.SetEnergy(1.e9);
    track.EnableDebugging();
  
    DriftLineRKF drift(sensor);
    drift.UseWeightingPotential(false);
    drift.EnableAvalanche();
    drift.EnableSignalCalculation();
    drift.SetGainFluctuationsPolya(0., 20000.);


    TCanvas* cD = new TCanvas("cD"," ",600,900);
    ViewDrift driftView;
    driftView.SetCanvas(cD);
    driftView.SetArea(-5, -5, 5, 5);
    drift.EnablePlotting(&driftView);
    track.EnablePlotting(&driftView);


    const double x0 = 19.0;
    const double y0 = 1.44;
    track.NewTrack(x0, y0, 0, 0, 0, -10, 0);


    std::cout << "Number of clusters: " << track.GetClusters().size() << "\n";
    for (const auto& cluster : track.GetClusters()) {
      	std::cout<<" electrons associ. w clusters = "<<cluster.electrons.size()<<"\n";
	std::cout<<" gain = "<<drift.GetGain()<<std::endl;
      for (const auto& electron : cluster.electrons) {
        drift.DriftElectron(electron.x, electron.y, electron.z, electron.t);
	double ne, ni;
	drift.GetAvalancheSize(ne, ni);
	std::cout<<" ne=   "<<ne<<" ni=   "<<ni<<std::endl;		  
      }

      for(const auto&ion:cluster.ions)
      {
          drift.DriftIon(ion.x,ion.y,ion.z,ion.t);
      }
      
    }

    //     TCanvas *c2 = new TCanvas("c2", "", 600, 600);
    cD->Clear();
    ViewCell cellView(comp);
    cellView.SetCanvas(cD);
    cellView.Plot2d();
    constexpr bool twod = true;
    constexpr bool drawaxis = false;
    driftView.Plot(twod, drawaxis);

    TCanvas *cS = new TCanvas("cS", "", 600, 600);
    sensor->PlotSignal("Anode", cS);


       // Get electric field at a specific location
    double x = 5, y = 1.0, z = 0.0;  // Position in cm
    double ex = 0., ey = 0., ez = 0.;
    int status = 0;
    Garfield::Medium* medium = nullptr;

    sensor->ElectricField(x, y, z, ex, ey, ez, medium, status);

    if (status != 0) {
      std::cerr << "Error: Could not retrieve electric field at position (" 
		<< x << ", " << y << ", " << z << ")\n";
    } else {
      std::cout << "Position: (" << x << ", " << y << ", " << z << ") cm\n";
      std::cout << "Electric field: (" << ex << ", " << ey << ", " << ez << ") V/cm\n";
    
      
      // Calculate field magnitude
      double e_mag = sqrt(ex*ex + ey*ey + ez*ez);
      std::cout << "E-field magnitude: " << e_mag << " V/cm = " 
		<< e_mag * 1.e-3 << " kV/cm\n\n";
      double bx = 0., by = 0., bz = 0.;
    
      // Get Townsend coefficient (alpha)
      double alpha = 0.;
      gas->ElectronTownsend(ex, ey, ez, bx, by, bz, alpha);
      std::cout << "Townsend coefficient (alpha): " << alpha << " / cm\n";
    
      // Get attachment coefficient (eta)
      double eta = 0.;
      gas->ElectronAttachment(ex, ey, ez, bx, by, bz, eta);
      std::cout << "Attachment coefficient (eta): " << eta << " / cm\n";
    
      // Effective Townsend coefficient (alpha - eta)
      std::cout << "Effective coefficient (alpha - eta): " 
		<< (alpha - eta) << " / cm\n";
    }
    
  
    
     
     /*
     TCanvas* c6 = new TCanvas("c6", "Cluster Gain and Signal", 1000, 500);
     c6->Divide(2,1);
     c6->cd(1);
     hGain->Draw();
     c6->cd(2);
     hSignal->Draw();
     */
     
     std::cout << "Simulation done.\n";
     
     app.Run(true);
     
     return 0;
}
     
     /*
     // Loop over the clusters along the track.
    std::cout << "Number of clusters: " << track.GetClusters().size() << "\n";
    for (const auto& cluster : track.GetClusters()) {
      // Loop over the electrons in the cluster.
      std::cout<<" electrons associ. w clusters = "<<cluster.electrons.size()<<"\n";
      for (const auto& electron : cluster.electrons) {
	drift.DriftElectron(electron.x, electron.y, electron.z, electron.t);
	aval.AvalancheElectron(electron.x, electron.y, electron.z, electron.t, electron.e, electron.dx, electron.dy, electron.dz);
	//const auto& p0 = electron.path[0];
	//drift.DriftIon(p0.x, p0.y, p0.z, p0.t);
	//drift.DriftIon(electron.x, electron.y, electron.z, electron.t);
      }
    std::cout<<" gain = "<<drift.GetGain()<<std::endl;
    std::cout<<" loss = "<<drift.GetLoss()<<std::endl;
    std::cout<<" path length = "<<drift.GetPathLength()<<std::endl;
    std::cout<<" number of driftline points = "<<drift.GetNumberOfDriftLinePoints()<<std::endl;
    std::cout<<" number of electrons = "<<drift.GetAvalancheSize().first<<std::endl;
    }
     */
    
    /*
    AvalancheMicroscopic aval(sensor);
    // position 
    double x0 = 2., y0 = 0.5, z0 = 0.0;
    double t0 = 0.;
    // initial energy eV
    double e0 = 0.1;
    // initial direction
    // in case of a null vector, the initial direction is randomized
    double dx0 = 0., dy0 = 0., dz0 = 0.;
    // calculate electron avalanche
    aval.AvalancheElectron(x0, y0, 0, t0, e0, dx0, dy0, dz0);
    
    AvalancheMC drift(sensor);
    drift.SetDistanceSteps(0.2);
    ViewDrift driftView;
    driftView.SetArea(-5, -5, 5, 5);
    //aval.EnablePlotting(&driftView);
    drift.EnablePlotting(&driftView);
  
    for (const auto& electron : aval.GetElectrons()) {
      const auto& p0 = electron.path[0];
      drift.DriftIon(p0.x, p0.y, p0.z, p0.t);
      }
    */
