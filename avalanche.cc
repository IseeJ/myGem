/**
 * avalanche.cc
 * General program flow based on example code from the Garfield++ website.
 *
 * Demonstrates electron avalanche and induced signal readout with
 * 2D finite-element visualization in Garfield++ with a LEM. 
 * LEM parameters are from:
 * C. Shalem et. al. Nucl. Instr. Meth. A, 558, 475 (2006).
 *
*/
#include <iostream>
#include <cmath>

#include <TCanvas.h>
#include <TApplication.h>
#include <TFile.h>

#include "Garfield/MediumMagboltz.hh"
#include "Garfield/ComponentElmer.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/ViewField.hh"
#include "Garfield/Plotting.hh"
#include "Garfield/ViewFEMesh.hh"
#include "Garfield/ViewSignal.hh"
#include "Garfield/GarfieldConstants.hh"
#include "Garfield/Random.hh"
#include "Garfield/AvalancheMicroscopic.hh"

using namespace Garfield;

int main(int argc, char* argv[]) {

  TApplication app("app", &argc, argv);
  plottingEngine.SetDefaultStyle();
  // Define the medium.
  // Setup the gas.
  MediumMagboltz* gas = new MediumMagboltz();
  gas->SetComposition("ar", 90., "ch4", 10.);
  gas->SetTemperature(293.15);
  gas->SetPressure(760.);
  gas->EnableDebugging();
  gas->Initialise();  
  gas->DisableDebugging();

  // Import an Elmer-created field map.
  ComponentElmer* elm = new ComponentElmer(
      "gemcell/mesh.header", "gemcell/mesh.elements", "gemcell/mesh.nodes",
      "gemcell/dielectrics.dat", "gemcell/gemcell.result", "cm");
  elm->EnablePeriodicityX();
  elm->EnableMirrorPeriodicityY();
  elm->SetMedium(0, gas);
  elm->PrintRange();
  // Import the weighting field for the readout electrode.
  // elm->SetWeightingField("gemcell/gemcell_WTlel.result", "wtlel");

  // Dimensions of the GEM [cm]
  const double pitch = 0.014;

  const bool plotField = true;
  if (plotField) {
    ViewField* fieldView = new ViewField();
    fieldView->SetComponent(elm);
    fieldView->SetPlane(0., -1., 0., 0., 0., 0.);
    fieldView->SetArea(-pitch / 2., -0.02, pitch / 2., 0.02);
    fieldView->SetVoltageRange(-160., 160.);
    TCanvas* cf = new TCanvas();
    fieldView->SetCanvas(cf);
    fieldView->PlotContour();
  }


  // Set the Penning transfer efficiency.
  const double rPenning = 0.51;
  const double lambdaPenning = 0.;
  gas->EnablePenningTransfer(rPenning, lambdaPenning, "ar");
  // Load the ion mobilities.
  const std::string path = getenv("GARFIELD_HOME");
  gas->LoadIonMobility(path + "/Data/IonMobility_Ar+_Ar.txt");
  // Associate the gas with the corresponding field map material. 
  const unsigned int nMaterials = elm->GetNumberOfMaterials();
  for (unsigned int i = 0; i < nMaterials; ++i) {
    const double eps = elm->GetPermittivity(i);
    if (eps == 1.) elm->SetMedium(i, gas);
  }
  elm->PrintMaterials();

  
  // Create the sensor.
  Sensor* sensor = new Sensor();
  sensor->AddComponent(elm);
  sensor->SetArea(-5 * pitch, -5 * pitch, -0.01,
                   5 * pitch,  5 * pitch,  0.025);

  AvalancheMicroscopic* aval = new AvalancheMicroscopic();
  aval->SetSensor(sensor);

  AvalancheMC* drift = new AvalancheMC();
  drift->SetSensor(sensor);
  drift->SetDistanceSteps(2.e-4);

  const bool plotDrift = true;
  ViewDrift* driftView = new ViewDrift();
  if (plotDrift) {
    driftView->SetArea(-2 * pitch, -2 * pitch, -0.02,
                        2 * pitch,  2 * pitch,  0.02); 
    aval->EnablePlotting(driftView);
    drift->EnablePlotting(driftView);
  }

 const int nEvents = 10;
  for (int i = nEvents; i--;) { 
    if (debug || i % 10 == 0) std::cout << i << "/" << nEvents << "\n";
    // Randomize the initial position. 
    double x0 = -pitch / 2. + RndmUniform() * pitch;
    double y0 = -pitch / 2. + RndmUniform() * pitch;
    double z0 = 0.02; 
    double t0 = 0.;
    double e0 = 0.1;
    aval->AvalancheElectron(x0, y0, z0, t0, e0, 0., 0., 0.);
    int ne = 0, ni = 0;
    aval->GetAvalancheSize(ne, ni);
    const int np = aval->GetNumberOfElectronEndpoints();
    double xe1, ye1, ze1, te1, e1;
    double xe2, ye2, ze2, te2, e2;
    double xi1, yi1, zi1, ti1;
    double xi2, yi2, zi2, ti2;
    int status;
    for (int j = np; j--;) {
      aval->GetElectronEndpoint(j, xe1, ye1, ze1, te1, e1, 
                                   xe2, ye2, ze2, te2, e2, status);
      drift->DriftIon(xe1, ye1, ze1, te1);
      drift->GetIonEndpoint(0, xi1, yi1, zi1, ti1, 
                               xi2, yi2, zi2, ti2, status);
    }
  }
  
  if (plotDrift) {
    TCanvas* cd = new TCanvas();
    driftView->SetCanvas(cd);
    driftView->Plot();
  }

  app.Run(kTRUE);

}
