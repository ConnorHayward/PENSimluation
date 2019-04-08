#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"
#include "SiliconPlateConstruction.hh"

#include "G4Material.hh"
#include "G4Element.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpticalSurface.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4Cons.hh"
#include "G4Torus.hh"
#include "G4Hype.hh"

#include "G4Transform3D.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "G4NistManager.hh"

#include "G4MultiFunctionalDetector.hh"
#include "G4SDManager.hh"
#include "G4PSEnergyDeposit.hh"
#include <G4VPrimitiveScorer.hh>

#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4VoxelLimits.hh"

#include "G4RunManager.hh"
#include "G4PhysicalConstants.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4Navigator.hh"
#include "G4TransportationManager.hh"

#include "G4GDMLParser.hh"

#include <G4VisAttributes.hh>
#include <iostream>
#include <fstream>

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/*
Constructs DetectorConstruction, defines default values.
*/
DetectorConstruction::DetectorConstruction()
 : G4VUserDetectorConstruction(),fPBox(nullptr), fLBox(nullptr),
  fBox(nullptr)
{
  fDetectorMessenger = new DetectorMessenger(this);
  fTargetMPT = new G4MaterialPropertiesTable();
  fExpHall_x = fExpHall_y = fExpHall_z = 2.0*m;
  fTargetName = "holder";
  fThickness = 1*mm;
  fTargetThickness = 3*mm;
  fDetectorType = 0;
  fABSL = 1;
  fRES=4.0;
  fSigAlpha = 0.5;
  fLY=10500./MeV;
  fDetectorName = "6pmt_coverage_pe";
  fVolName = "World";
  DefineMaterials();
  // SetTargetMaterial("PEN");
  SetWorldMaterial("Air");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction(){;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
/*
Sets thickness of target.
*/
void DetectorConstruction::SetSize(G4double value){
  fTargetThickness=value;
  if(fBox){
    fBox->SetZHalfLength(fTargetThickness/2);
  }
  UpdateGeometry();

  G4RunManager::GetRunManager()->PhysicsHasBeenModified();
}

void DetectorConstruction::SetLY(G4double value){
  fLY=value;
  UpdateGeometry();
  G4RunManager::GetRunManager()->PhysicsHasBeenModified();
}

void DetectorConstruction::SetRes(G4double value){
  fRES=value;
  UpdateGeometry();

}

void DetectorConstruction::SetSigAlpha(G4double value){
  fSigAlpha = value;
  UpdateGeometry();
  G4RunManager::GetRunManager()->PhysicsHasBeenModified();
}

/*
Sets which detector geometry is used.
*/
void DetectorConstruction::SetDetectorType(G4int value){
  fDetectorType=value;

  UpdateGeometry();
  G4RunManager::GetRunManager()->PhysicsHasBeenModified();
}

void DetectorConstruction::SetABS(G4double value){
  fABSL=value;

  UpdateGeometry();
  G4RunManager::GetRunManager()->PhysicsHasBeenModified();
}

void DetectorConstruction::SetDetectorName(G4String name){
  fDetectorName=name;

  UpdateGeometry();
  G4RunManager::GetRunManager()->PhysicsHasBeenModified();
}

/*
Sets material of target.
*/
void DetectorConstruction::SetTargetMaterial(G4String materialChoice)
{
  // search the material by its name
  G4Material* pttoMaterial =
     G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);

  if (pttoMaterial) {
    fTargetMaterial = pttoMaterial;
    fTargetName = fTargetMaterial->GetName();
    if ( fLBox ) { fLBox->SetMaterial(fTargetMaterial); }
  } else {
    G4cout << "\n--> warning from DetectorConstruction::SetMaterial : "
           << materialChoice << " not found" << G4endl;
  }
  G4RunManager::GetRunManager()->PhysicsHasBeenModified();
}

/*
Sets material of world volume.
*/
void DetectorConstruction::SetWorldMaterial(G4String materialChoice)
{
  // search the material by its name
  G4Material* pttoMaterial =
     G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);

  if (pttoMaterial) {
    fWorldMaterial = pttoMaterial;
    if ( fWLBox ) { fWLBox->SetMaterial(fWorldMaterial); }
  } else {
    G4cout << "\n--> warning from DetectorConstruction::SetMaterial : "
           << materialChoice << " not found" << G4endl;
  }
  G4RunManager::GetRunManager()->PhysicsHasBeenModified();
}

/*
Defines materials used in simulation. Sets material properties for PEN and other optical components.
*/
void DetectorConstruction::DefineMaterials(){// ------------- Materials -------------
  G4double a, z, density;
  G4int nelements;

  // fAir
  //
  G4Element* N = new G4Element("Nitrogen", "N", z=7 , a=14.01*g/mole);
  G4Element* O = new G4Element("Oxygen"  , "O", z=8 , a=16.00*g/mole);

  fAir = new G4Material("Air", density=1.29*mg/cm3, nelements=2);
  fAir->AddElement(N, 70.*perCent);
  fAir->AddElement(O, 30.*perCent);

  G4NistManager* man = G4NistManager::Instance();

  fGlass = man->FindOrBuildMaterial("G4_Pyrex_Glass");
  fTeflon = man->FindOrBuildMaterial("G4_TEFLON");
  fLAr = man->FindOrBuildMaterial("G4_lAr");
  fAl = man->FindOrBuildMaterial("G4_Al");
  fSi = man->FindOrBuildMaterial("G4_Si");

  // Water
  //
  G4Element* H = new G4Element("Hydrogen", "H", z=1 , a=1.01*g/mole);

  G4Material* water = new G4Material("Water", density= 1.0*g/cm3, nelements=2);
  water->AddElement(H, 2);
  water->AddElement(O, 1);

  G4Element* C = new G4Element("Carbon", "C", z=12, a=12*g/mole);
  G4Element* Pb = new G4Element("Lead", "Pb", z=87, a=207*g/mole);

  // Scintillators

  fPEN = new G4Material("PEN", density= 1.3*g/cm3, nelements=3);
  G4int number_of_atoms;
  fPEN->AddElement(O, number_of_atoms=4);
  fPEN->AddElement(H, number_of_atoms=10);
  fPEN->AddElement(C, number_of_atoms=14);

  G4double wavelength;
  char filler;
  G4double varabsorlength;
  G4double ems;
  G4double rindex;

  G4double absEnergy[102]  = {0};
  G4double abs[102]={0};
  G4double emission[102]={0};
  G4double rIndex[102]={0};
  G4double rIndex_fAir[102]={0};
  G4double ems_abs[102]={0};

  G4int absEntries = 0;
  ifstream ReadAbs;

  G4String abs_file = "../input_files/Exp4_long.csv";
  G4double emission_fibre[102]={0};
  ReadAbs.open(abs_file);
  G4double var = GetABS();
  if(ReadAbs.is_open())
  {
    while(!ReadAbs.eof())
    {
      ReadAbs>>wavelength>>filler>>varabsorlength>>filler>>ems>>filler>>rindex;
      if(ReadAbs.eof()){
        break;
      }
      absEnergy[absEntries] = (1240/wavelength)*eV;
      abs[absEntries] = 50*cm;;
      emission[absEntries] = ems;
      rIndex[absEntries] = 1.65; // 1.4906 for PMMA
      rIndex_fAir[absEntries]=1.0;
      ems_abs[absEntries]=0.02;
      emission_fibre[absEntries]=1.0;
      absEntries++;
    }
  }

  else G4cout<<"Error opening file: " <<abs_file<<G4endl;
  ReadAbs.close();
  absEntries--;

  const G4int nEntries1 = sizeof(absEnergy)/sizeof(G4double);
  assert(sizeof(rIndex) == sizeof(absEnergy));
  assert(sizeof(abs) == sizeof(absEnergy));
  assert(sizeof(emission) == sizeof(absEnergy));
  assert(sizeof(rIndex_fAir == sizeof(absEnergy)));

  G4double uv_range[2]={3.099605,5.0};
  G4double uv_abs[2]={0.00001,0.00001};

  fTargetMPT->AddProperty("WLSABSLENGTH",uv_range,abs,2)->SetSpline(true);
  fTargetMPT->AddProperty("WLSCOMPONENT",absEnergy, emission_fibre, nEntries1)->SetSpline(true);
  fTargetMPT->AddConstProperty("WLSTIMECONSTANT", 12*ns);

  fTargetMPT->AddProperty("RINDEX",       absEnergy, rIndex, nEntries1)->SetSpline(true);
  fTargetMPT->AddProperty("ABSLENGTH",    absEnergy, abs, nEntries1)->SetSpline(true); // *
  fTargetMPT->AddProperty("FASTCOMPONENT",absEnergy, emission, nEntries1)->SetSpline(true);
  fTargetMPT->AddProperty("SLOWCOMPONENT",absEnergy, emission, nEntries1)->SetSpline(true);

  fTargetMPT->AddConstProperty("SCINTILLATIONYIELD",10500./MeV); // * 2.5 * PEN = PS, 10*PEN=PS
  fTargetMPT->AddConstProperty("RESOLUTIONSCALE",4.0); // * 1, 4, 8
  fTargetMPT->AddConstProperty("FASTTIMECONSTANT", 5.198*ns);
  fTargetMPT->AddConstProperty("SLOWTIMECONSTANT",24.336*ns);
  fTargetMPT->AddConstProperty("YIELDRATIO",0.05);

  fPEN->SetMaterialPropertiesTable(fTargetMPT);

  density = universe_mean_density;    //from PhysicalConstants.h
  fVacuum = new G4Material("Galactic", z=1., a=1.008*g/mole, density,
                           kStateGas,2.73*kelvin,3.e-18*pascal);
  //
  // fAir
  G4MaterialPropertiesTable* worldMPT = new G4MaterialPropertiesTable();
  worldMPT->AddProperty("RINDEX", absEnergy, rIndex_fAir, nEntries1)->SetSpline(true);

  fAir->SetMaterialPropertiesTable(worldMPT);
  fVacuum->SetMaterialPropertiesTable(worldMPT);
}

void DetectorConstruction::SetVolName(G4ThreeVector thePoint){
  G4Navigator* theNavigator = G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();
  G4VPhysicalVolume* myVolume= theNavigator->LocateGlobalPointAndSetup(thePoint);
  fVolName =  myVolume->GetName();
}

void DetectorConstruction::UpdateGeometry(){
  G4RunManager::GetRunManager()->DefineWorldVolume(Construct());
}

/*
Clears stored geometry, then constructs all volumes that can be used in the simulation.

Builds and places volumes in world.

Defines detector sensitivities and properties.
*/
G4VPhysicalVolume* DetectorConstruction::Construct()
{
    G4GDMLParser parser;
    G4GeometryManager::GetInstance()->OpenGeometry();
    G4LogicalVolumeStore::GetInstance()->Clean();
    G4PhysicalVolumeStore::GetInstance()->Clean();
    G4SolidStore::GetInstance()->Clean();

// ------------- Volumes --------------

// The experimental Hall
  fWorldBox = new G4Box("World",fExpHall_x,fExpHall_y,fExpHall_z);

  fWLBox = new G4LogicalVolume(fWorldBox,fAir,"World",0,0,0);

  fWPBox = new G4PVPlacement(0,G4ThreeVector(),fWLBox,"World",0,false,0);

  double thickness_foil = 2.5*mm;
  double holeHeight = fThickness;
  G4Box* absorberPb_box = new G4Box("absorber",15*mm,15*mm, 1.5*mm);
  G4Tubs* hole = new G4Tubs("hole",0*mm,1.5*mm,holeHeight, 0.*deg, 360.*deg);
  G4SubtractionSolid* collimator = new G4SubtractionSolid("collimator",absorberPb_box,hole);

  G4Box* absorber_foil = new G4Box("foil",15*mm,15*mm,thickness_foil);
  G4Tubs* hole_foil = new G4Tubs("hole_foil",0*mm,1.5*mm,thickness_foil, 0.*deg, 360.*deg);
  G4SubtractionSolid* collimator_foil = new G4SubtractionSolid("collimator_foil",absorber_foil,hole_foil);

  G4Tubs* glass_plate = new G4Tubs("plate",0*mm, 5*cm, 0.75*mm, 0.*deg, 360.*deg);

  double target_width = 15*mm;
  fBox = new G4Box("target", target_width, target_width, fTargetThickness);
  fLBox = new G4LogicalVolume(glass_plate,fPEN, "target",0,0,0);
  double position = 0;

  double rod_sides = 2.5*mm;
  double rod_length = 48*mm;

  G4Box* rod = new G4Box("target",rod_sides,rod_sides,rod_length);
  G4LogicalVolume* rod_log = new G4LogicalVolume(rod,fPEN,"target",0,0,0);

  G4NistManager* man = G4NistManager::Instance();

  // --------------Detectors--------------

  char filler;
  G4double wavelength;

  G4double cath_eff;
  G4double photocath_energy[57];
  G4double perfect_EFF[57];
  G4double perfect_REFL[57];
  G4String pmt_file = "../input_files/pmtQE.csv";

  ifstream ReadEff;
  G4int effCounter = 0;
  ReadEff.open(pmt_file);

  if(ReadEff.is_open())
  {
    while(!ReadEff.eof())
    {
      ReadEff>>wavelength>>filler>>cath_eff;
      if(ReadEff.eof()){
        break;
      }
      photocath_energy[57-effCounter] = (1240/wavelength)*eV;
      perfect_EFF[57-effCounter] = 1;
      perfect_REFL[57-effCounter] = 0;
      effCounter++;
    }
  }

  else G4cout<<"Error opening file: " <<pmt_file<<G4endl;
  ReadEff.close();
  effCounter--;

  const G4int nPMT_EFF = sizeof(photocath_energy)/sizeof(G4double);

  G4OpticalSurface* perfect_optsurf = new G4OpticalSurface("perfect",glisur,polished, dielectric_metal);
  G4MaterialPropertiesTable* detector_MT = new G4MaterialPropertiesTable();
  detector_MT->AddProperty("EFFICIENCY", photocath_energy, perfect_EFF,nPMT_EFF);
  detector_MT->AddProperty("REFLECTIVITY", photocath_energy, perfect_REFL,nPMT_EFF);
  perfect_optsurf->SetMaterialPropertiesTable(detector_MT);

  G4LogicalVolume* tile_detector = new G4LogicalVolume(fBox, fSi, "tile_detector");
  new G4LogicalSkinSurface("perfect_sensor",tile_detector, perfect_optsurf);

  // Spectrometer Sensor
//  new G4LogicalSkinSurface("perfect_surf",tile_detector,perfect_optsurf);

  G4VPhysicalVolume* siPM_placement;


  G4ThreeVector* pmtVector = new G4ThreeVector(0,0,-(1+2*1+fTargetThickness));
  fSourceVector = new G4ThreeVector(0,0,(1+2*1+fTargetThickness));

  // Set Draw G4VisAttributes

  G4VisAttributes* visAttr = new G4VisAttributes();
  visAttr->SetVisibility(false);
  fWLBox->SetVisAttributes(visAttr);

  // Scintillator Materials
  G4VisAttributes* tileAttr = new G4VisAttributes(G4Colour::Blue());
  tileAttr->SetVisibility(true);
  fLBox->SetVisAttributes(tileAttr);


  // // Active Detectors
  G4VisAttributes* detectorAttr = new G4VisAttributes(G4Colour::Green());
  detectorAttr->SetVisibility(true);
  detectorAttr->SetForceSolid(false);


  // // Inactive Regions
  G4VisAttributes* supportAttr = new G4VisAttributes(G4Colour::Grey());
  supportAttr->SetVisibility(true);
  supportAttr->SetForceSolid(true);

  G4RotationMatrix* rotationMatrix = new G4RotationMatrix(0,0,0);
  /*
  0 - PMT on base of tile, collimator included.
  */
  fDetectorType = 0;

  fPBox = new G4PVPlacement(0, G4ThreeVector(0,0,0),fLBox,"rod",fWLBox,false,0,true);
  siPM_placement = new G4PVPlacement(0, G4ThreeVector(0,0,-2.5*cm),tile_detector,"spec",fWLBox,false,0,true);

  G4OpticalSurface* AirPEN = new G4OpticalSurface("AirPEN",glisur, ground, dielectric_dielectric);
  AirPEN -> SetPolish(0.8);
  AirPEN -> SetMaterialPropertiesTable(fTargetMPT);

  G4LogicalBorderSurface* surfaceAirPEN = new G4LogicalBorderSurface("AirPEN",fWPBox,fPBox,AirPEN);
  G4LogicalBorderSurface* surfacePENAir = new G4LogicalBorderSurface("AirPEN",fPBox,fWPBox,AirPEN);

  return fWPBox;
}
