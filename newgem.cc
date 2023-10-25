#include <cstdlib>
#include <iostream>
#include <fstream>

#include <TApplication.h>
#include <TCanvas.h>
#include <TH1F.h>

#include "Garfield/ComponentAnsys123.hh"
#include "Garfield/ViewField.hh"
#include "Garfield/ViewFEMesh.hh"
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/AvalancheMicroscopic.hh"
#include "Garfield/AvalancheMC.hh"
#include "Garfield/Random.hh"

using namespace Garfield;

int main(int argc, char * argv[]) {

  TApplication app("app", &argc, argv);

  // Setup the gas.
  MediumMagboltz gas("ar", 80., "co2", 20.);
  gas.SetTemperature(293.15);
  gas.SetPressure(760.);
  gas.Initialise(true);  
  // Set the Penning transfer efficiency.
  constexpr double rPenning = 0.51;
  constexpr double lambdaPenning = 0.;
  gas.EnablePenningTransfer(rPenning, lambdaPenning, "ar");
  // Load the ion mobilities.
  const std::string path = std::getenv("GARFIELD_INSTALL");
  gas.LoadIonMobility(path + "/share/Garfield/Data/IonMobility_Ar+_Ar.txt");

  // Load the field map.
  ComponentElmer* elm = new ComponentElmer("gemcell/mesh.header", "gemcell/mesh.elements", "gemcell/mesh.nodes","gemcell/dielectrics.dat", "gemcell/gemcell.result", "cm");
  elm.EnableMirrorPeriodicityX();
  elm.EnableMirrorPeriodicityY();
  elm.PrintRange();

  // Associate the gas with the corresponding field map material.
  elm.SetGas(&gas); 
  elm.PrintMaterials();
  // elm.Check();

  // Dimensions of the GEM [cm]
  constexpr double pitch = 0.014;

  ViewField fieldView;
  ViewFEMesh meshView;
  constexpr bool plotField = true;
  if (plotField) {
    fieldView.SetComponent(&elm);
    // Set the normal vector of the viewing plane (xz plane).
    fieldView.SetPlane(0, -1, 0, 0, 0, 0);
    // Set the plot limits in the current viewing plane.
    fieldView.SetArea(-0.5 * pitch, -0.02, 0.5 * pitch, 0.02);
    fieldView.SetVoltageRange(-160., 160.);
    TCanvas* cf = new TCanvas("cf", "", 600, 600);
    cf->SetLeftMargin(0.16);
    fieldView.SetCanvas(cf);
    fieldView.PlotContour();

    meshView.SetArea(-0.5 * pitch, -0.02, 0.5 * pitch, 0.02); 
    meshView.SetCanvas(cf);
    meshView.SetComponent(&elm);
    meshView.SetPlane(0, -1, 0, 0, 0, 0);
    meshView.SetFillMesh(true);
    meshView.SetColor(2, kGray);
    meshView.Plot(true);
  }

  // Create the sensor.
  Sensor sensor;
  sensor.AddComponent(&elm);
  sensor.SetArea(-5 * pitch, -5 * pitch, -0.01,
                  5 * pitch,  5 * pitch,  0.025);

  AvalancheMicroscopic aval;
  aval.SetSensor(&sensor);

  AvalancheMC drift;
  drift.SetSensor(&sensor);
  drift.SetDistanceSteps(2.e-4);

  ViewDrift driftView;
  constexpr bool plotDrift = true;
  if (plotDrift) {
    // Plot every tenth collision.
    aval.EnablePlotting(&driftView, 10);
    drift.EnablePlotting(&driftView);
  }
