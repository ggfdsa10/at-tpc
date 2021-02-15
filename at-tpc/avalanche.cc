/*
 * AT-TPC Simulation
 * (gas chamber and gem 3 layer)  
*/

#include <iostream>
#include <fstream>
#include <TCanvas.h>
#include <TApplication.h>
#include <stdio.h>
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/ComponentElmer.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/ViewField.hh"
#include "Garfield/Plotting.hh"
#include "Garfield/ViewFEMesh.hh"
#include "Garfield/GarfieldConstants.hh"
#include "Garfield/Random.hh"
#include "Garfield/AvalancheMicroscopic.hh"

using namespace Garfield;
using namespace std;


int main(int argc, char* argv[]) {

  TApplication app("app", &argc, argv);

  // Detector parameter
  double pitch = 0.014; // Gem pitch [cm]
  double gemy = 0.0121243; // Gem cell y-length (sqrt(3)*70 um) [cm]
  double gem_z = 0.41089; // Gem layer height [cm]
  double detector_x = 10; // Detector X-length [cm]
  double detector_y = 10; // Detector Y-length [cm]
  double gas_z = 14.7+0.343; // Gas chamber Z-length [cm]

  // Event parameter
  int epoch = 100; // Event epoch
  double electron_xi = 0.; // Initial electron X-position [cm]
  double electron_yi = 0.; // Initial electron Y-position [cm]
  double electron_zi = 14.; // Initial electron Z-position [cm]
  int collision_step = 1000; //Drift line plot step

  // Gas(P10) parameter
  double temperature = 293.15; // Detector temperature [K] 
  double pressure = 740.; // Detector pressure [Torr] 
  double rPenning = 0.57; // Value of transfer efficiency
  double lambdaPenning = 0.; // Penning transfer distance
  double maxenergy = 300.; // Max electron energy [eV]


  /* Variable declear
   * x0, y0, z0, t0, e0 = Initial electrons value
   * x1, y1, z1, t1, e1, dx1, dy1, dz1 = Electrons value in gem 3 layer
   * x2, y2, z2, t2, e2 = Final electron value
   */
  double x0,y0,z0,t0,e0,x1,y1,z1,t1,e1,dx1,dy1,dz1,x2,y2,z2,t2,e2;
  int status1, status2;
  

  // Gas setting
  MediumMagboltz* gas = new MediumMagboltz();
  gas -> SetTemperature(temperature);  
  gas -> SetPressure(pressure);       
  gas -> EnableDrift(true);           
  gas -> SetComposition("ar", 90., "ch4", 10.);
  gas -> EnablePenningTransfer(rPenning, lambdaPenning, "ar");
  std::string path = getenv("GARFIELD_HOME");
  gas -> LoadIonMobility(path + "/Data/IonMobility_Ar+_Ar.txt");
  gas -> EnableThermalMotion(true);
  gas -> SetMaxElectronEnergy(maxenergy);

  
  //Gas chamber 
  ComponentElmer* p10 = new ComponentElmer("gas_chamber/mesh.header", "gas_chamber/mesh.elements", "gas_chamber/mesh.nodes", "gas_chamber/dielectrics.dat", "gas_chamber/gas.result", "cm");
  p10 -> SetMedium(0, gas);

  //Gem layer region
  ComponentElmer* gem = new ComponentElmer("gem/mesh.header", "gem/mesh.elements", "gem/mesh.nodes", "gem/dielectrics.dat", "gem/gem.result", "cm");
  gem -> EnableMirrorPeriodicityX();
  gem -> EnableMirrorPeriodicityY();
  gem -> SetMedium(0, gas);

  //Gas chamber sensor
  Sensor* sensor1 = new Sensor();
  sensor1 -> AddComponent(p10);
  sensor1 -> SetArea(-detector_x/2, -detector_y/2, 0.343, detector_x/2, detector_y/2, gas_z+0.0001);

  //Gem layer sensor
  Sensor* sensor2 = new Sensor();
  sensor2 -> AddComponent(gem);
  sensor2 -> SetArea(-detector_x/2, -detector_y/2, -0.0001, detector_x/2, detector_y/2, gem_z+0.0001);

  //Gas chamber avalanche
  AvalancheMicroscopic* aval1 = new AvalancheMicroscopic();
  aval1 -> SetSensor(sensor1);
  aval1 -> SetCollisionSteps(collision_step);

  //Gem layer avalanche
  AvalancheMicroscopic* aval2 = new AvalancheMicroscopic();
  aval2 -> SetSensor(sensor2);
  aval2 -> SetCollisionSteps(collision_step);

  //Gas chamber drift
  ViewDrift* viewDrift1 = new ViewDrift();
  viewDrift1 -> SetArea(-detector_x/2, -detector_y/2, 0.39, detector_x/2, detector_y/2, gas_z+0.0001);
  aval1 -> EnablePlotting(viewDrift1);

  //Gem layer drift
  ViewDrift* viewDrift2 = new ViewDrift();
  viewDrift2 -> SetArea(-detector_x/2, -detector_y/2, -0.0001, detector_x/2, detector_y/2, gem_z+0.0001);
  aval2 -> EnablePlotting(viewDrift2);



  // simulation and Data output
  
  ofstream data;
  ofstream information;
  
  information.open("information.txt");
  
  double averge_en_total = 0;
  double averge_en_pad = 0;
  
  for (int i =0; i<epoch; i++){
    aval1 -> AvalancheElectron(electron_xi, electron_yi, electron_zi, 0., 0., 0., 0., 0.);

    int gas_en = aval1 -> GetNumberOfElectronEndpoints();
   
    cout << " epoch : " << i+1 << endl;
    cout << " The number of electrons generated in gas_chamber : " << gas_en << endl;
   
    char name[epoch];
    sprintf(name, "data_%d.txt", i+1);
    data.open(name);

    double final_en =0;
    double total_final_en =0;

    for(int j =0; j<gas_en; j++){
      aval1 -> GetElectronEndpoint(j, x0, y0, z0, t0, e0, x1, y1, z1, t1, e1, dx1, dy1, dz1, status1);
      
      if (z1 = 0.4){
	aval2 -> AvalancheElectron(x1, y1, z1, t1, e1, dx1, dy1, dz1); 
      
	int gem_en = aval2 -> GetNumberOfElectronEndpoints();
    
	for (int k = 0; k <gem_en; k++){
	  aval2 -> GetElectronEndpoint(k, x1, y1, z1, t1, e1, x2, y2, z2, t2, e2, status2);
	  if (z2 <=0.){
	    data << x0 << " " << y0 << " " << z0 << " " << x2 << " " << y2 << " " << z2 << " " << t2 << " " << e2 << endl;
	    final_en += 1;
	  }
	  total_final_en +=1;
	}
      }
    }
    data.close();

    cout << " The number of electrons generated in gem 3 layer : " << total_final_en << endl;
    cout << " The number of electrons generated on the pad : " << final_en << endl << endl;

    information << " epoch : " << i+1 << endl;
    information << " The number of electrons generated in gas_chamber : " << gas_en << endl;
    information << " The number of electrons generated in gem 3 layer : " << total_final_en << endl;
    information << " The number of electrons generated on the pad : " << final_en << endl;
    information << " Gem effieiency : " << (final_en / total_final_en)*100 << "%" << endl << endl;
    
    averge_en_total += total_final_en;
    averge_en_pad += final_en;
  }
  
  cout << "calculation finish" << endl;
  
  information << " Averge number of electrons generated in gem 3 layer : " << averge_en_total/epoch << endl;
  information << " Averge number of electrons generated on the pad : " << averge_en_pad/epoch << endl;
  information << " Averge Gem effieiency : " << (averge_en_total/averge_en_pad)*100 << "%" << endl;
 
  information.close();



  

  /*
  ViewField* vf = new ViewField();
  vf->SetSensor(sensor);
  vf->SetArea(-10*pitch, -gemy, -0.0001, 10*pitch, gemy, detector_z+0.0001);
  vf->SetNumberOfContours(100);
  vf->SetNumberOfSamples2d(30, 30);
  vf->SetPlane(0, -1, 0, 0, 0, 0);
  vf->SetVoltageRange(-2000, 0);
  vf->SetElectricFieldRange(0,80000);
  vf->PlotContour("ez");
  */
  
  // Gas chamber view
  ViewFEMesh* vFE1 = new ViewFEMesh();
  vFE1 -> SetArea(-5, -5, 0., 5, 5, gas_z+0.001);
  vFE1 -> SetComponent(p10);
  vFE1 -> SetPlane(0, -1, 0, 0, 0, 0);
  vFE1 -> SetFillMesh(true);
  vFE1 -> SetColor(1, kRed);
  vFE1 -> SetColor(2, kRed);
  vFE1 -> EnableAxes();
  vFE1 -> SetXaxisTitle("x (cm)");
  vFE1 -> SetYaxisTitle("z (cm)");
  vFE1 -> SetViewDrift(viewDrift1);
  vFE1 -> Plot();

  // Gem layer view
  ViewFEMesh* vFE2 = new ViewFEMesh();
  vFE2 -> SetArea(-4*pitch, -4*gemy, -0.001, 4*pitch, 4*gemy, gem_z+0.001);
  vFE2 -> SetComponent(gem);
  vFE2 -> SetPlane(0, -1, 0, x1, 0, 0);
  vFE2 -> SetFillMesh(true);
  vFE2 -> SetColor(1, kRed);
  vFE2 -> SetColor(2, kRed);
  vFE2 -> SetColor(3, kRed);
  vFE2 -> SetColor(4, kRed);
  vFE2 -> SetColor(5, kRed);
  vFE2 -> SetColor(6, kRed);
  vFE2 -> SetColor(7, kRed);
  vFE2 -> SetColor(8, kRed);
  vFE2 -> SetColor(9, kRed);
  vFE2 -> SetColor(10, kRed);
  vFE2 -> SetColor(11, kRed);
  vFE2 -> EnableAxes();
  vFE2 -> SetXaxisTitle("x (cm)");
  vFE2 -> SetYaxisTitle("z (cm)");
  vFE2 -> SetViewDrift(viewDrift2);
  vFE2 -> Plot();
  
  app.Run(kTRUE);

  return 0;
}
