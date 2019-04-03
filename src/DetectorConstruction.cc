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
  fLY=10500./MeV;
  fDetectorName = "6pmt_coverage_pe";
  fVolName = "World";
  DefineMaterials();
//  SetTargetMaterial("Scint");
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
  fPb = man->FindOrBuildMaterial("G4_Pb");
  fCu = man->FindOrBuildMaterial("G4_Cu");
  fGe = man->FindOrBuildMaterial("G4_Ge");

  G4MaterialPropertiesTable* glassMPT = new G4MaterialPropertiesTable();
  G4double photonEnergy2[] = {2.479684*eV, 2.610194*eV, 2.755204*eV, 2.883353*eV, 2.917275*eV, 3.099605*eV};
  G4double refractiveIndex4[] = {1.52, 1.52, 1.52, 1.52, 1.52, 1.52};
  G4double absorption[] = {10*m,10*m,10*m,10*m,10*m,10*m};
  assert(sizeof(refractiveIndex4) == sizeof(photonEnergy2));
  glassMPT->AddProperty("RINDEX",photonEnergy2,refractiveIndex4,6)->SetSpline(true);
  glassMPT->AddProperty("ABSLENGTH",photonEnergy2,absorption,6);
  fGlass->SetMaterialPropertiesTable(glassMPT);

  // Water
  //
  G4Element* H = new G4Element("Hydrogen", "H", z=1 , a=1.01*g/mole);

  G4Material* water = new G4Material("Water", density= 1.0*g/cm3, nelements=2);
  water->AddElement(H, 2);
  water->AddElement(O, 1);

  G4Element* C = new G4Element("Carbon", "C", z=12, a=12*g/mole);
  G4Element* Pb = new G4Element("Lead", "Pb", z=87, a=207*g/mole);

  fWood = new G4Material
    ("wood", density=0.9*g/cm3, nelements=3);
  fWood->AddElement(H , 4);
  fWood->AddElement(O , 1);
  fWood->AddElement(C , 2);

  // Plastic casings

  fPOM = new G4Material("POM",density=1.41*g/cm3,nelements=3);
  fPOM->AddElement(O,1);
  fPOM->AddElement(C,1);
  fPOM->AddElement(H,2);

  fABS = new G4Material("ABS",density=1.07*g/cm3,nelements=3);
  fABS->AddElement(C,15);
  fABS->AddElement(H,17);
  fABS->AddElement(N,1);

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
      abs[absEntries] = varabsorlength;
      emission[absEntries] = ems;
      rIndex[absEntries] = fLY; // 1.4906 for PMMA
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

  fScintilator = new G4Material("Scint", density= 1.03*g/cm3, 2);
  fScintilator->AddElement(C, 0.475);
  fScintilator->AddElement(H, 0.525);

  Pstyrene = new G4Material("Polystyrene", density= 1.03*g/cm3, 2);
  Pstyrene->AddElement(C, 8);
  Pstyrene->AddElement(H, 8);

  density = universe_mean_density;    //from PhysicalConstants.h
  fVacuum = new G4Material("Galactic", z=1., a=1.008*g/mole, density,
                           kStateGas,2.73*kelvin,3.e-18*pascal);
  //
  // fAir
  G4MaterialPropertiesTable* worldMPT = new G4MaterialPropertiesTable();
  worldMPT->AddProperty("RINDEX", absEnergy, rIndex_fAir, nEntries1)->SetSpline(true);

  fAir->SetMaterialPropertiesTable(worldMPT);
  fVacuum->SetMaterialPropertiesTable(worldMPT);

  G4MaterialPropertiesTable* lARMPT = new G4MaterialPropertiesTable();
  abs_file = "../input_files/lArScint.csv";

  ReadAbs.open(abs_file);
  absEntries = 0;
  G4double absEnergylArScint[43] = {};
  G4double emissionlArScint[43] = {};

  if(ReadAbs.is_open())
  {
    while(!ReadAbs.eof())
    {
      ReadAbs>>wavelength>>filler>>ems;
      if(ReadAbs.eof()){
        break;
      }
      absEnergylArScint[absEntries] = (1240/wavelength)*eV;
      emissionlArScint[absEntries] = ems;
      absEntries++;
    }
  }
  else G4cout<<"Error opening file: " <<abs_file<<G4endl;
  ReadAbs.close();
  absEntries--;

  G4double absEnergyLAr[]  = {1.9255*eV, 2.145*eV, 2.2704*eV, 2.4378*eV, 2.6085*eV,2.845*eV,3.0515*eV,3.397*eV};
  G4double absLAr[]={10*m,10*m,10*m,10*m,10*m,10*m,10*m,10*m};
  G4double rIndexLAr[]={1.2295, 1.2303, 1.2308,1.2316,1.2324,1.2336,1.2347,1.2367};

  lARMPT->AddProperty("RINDEX",       absEnergyLAr, rIndexLAr, 8)->SetSpline(true);
  lARMPT->AddProperty("ABSLENGTH",    absEnergyLAr, absLAr, 8)->SetSpline(true);
  lARMPT->AddProperty("FASTCOMPONENT",absEnergylArScint, emissionlArScint, 43)->SetSpline(true);
  lARMPT->AddProperty("SLOWCOMPONENT",absEnergylArScint, emissionlArScint, nEntries1)->SetSpline(true);
  lARMPT->AddConstProperty("SCINTILLATIONYIELD", 51.*MeV);
  lARMPT->AddConstProperty("FASTTIMECONSTANT", 6.2*us);
  lARMPT->AddConstProperty("SLOWTIMECONSTANT",1300*us);
  lARMPT->AddConstProperty("YIELDRATIO",0.05);
  lARMPT->AddConstProperty("RESOLUTIONSCALE",1.0);
  fLAr->SetMaterialPropertiesTable(lARMPT);

  G4double rindexEnergy[500] = {0};
  G4double scintIndex[500] = {0};

  G4int rindexEntries = 0;
  ifstream ReadRindex;

  G4String rindex_file="../input_files/rindexScint.txt";
  ReadRindex.open(rindex_file);

  if(ReadRindex.is_open())
  {
      while(!ReadRindex.eof())
      {
          ReadRindex>>wavelength>>filler>>scintIndex[76-rindexEntries];
          rindexEnergy[76 - rindexEntries] = (1240/wavelength)*eV;
          rindexEntries++;
      }
  }else G4cout<<"Error opening file: "<<rindex_file<<G4endl;
  ReadRindex.close();
  rindexEntries--;

  G4double scintEnergy[501] = {0};
  G4double scintEmit[501] = {0};
  G4double scintEmitSlow[501] = {0};

  G4int scintEntries = 0;
  ifstream ReadScint;

  G4String Scint_file="../input_files/pTP_emission.txt";
  ReadScint.open(Scint_file);

  if(ReadScint.is_open())
  {
      while(!ReadScint.eof())
      {
          ReadScint>>wavelength>>filler>>scintEmit[500-scintEntries];
          //convert wavelength to eV:
          scintEnergy[500 - scintEntries] = (1240/wavelength)*eV;
          scintEmitSlow[500 - scintEntries] = scintEmit[500 - scintEntries];
          scintEntries++;
      }
  }else G4cout<<"Error opening file: "<<Scint_file<<G4endl;
  ReadScint.close();
  scintEntries--;

  G4int absorbEntries = 0;
  G4double varabsorblength;
  G4double absorbEnergy[501] = {0};
  G4double Absorb[501] = {0};

  ifstream ReadAbsorb;
  G4String ReadAbsorbLength="../input_files/PlasticBulkAbsorb2.cfg";

  ReadAbsorb.open(ReadAbsorbLength);
  if (ReadAbsorb.is_open())
  {
      while(!ReadAbsorb.eof())
      {
          ReadAbsorb>>wavelength>>filler>>varabsorblength;
          absorbEnergy[500 - absorbEntries]=(1240/wavelength)*eV;
          Absorb[500 - absorbEntries]=varabsorblength*m;
          absorbEntries++;
      }
  }else G4cout<<"Error opening file: "<<ReadAbsorb<<G4endl;
  ReadAbsorb.close();
  absorbEntries--;

  G4double wlsEnergy[501] = {0};
  G4double wlsEmit[501] = {0};

  G4int wlsScintEntries = 0;
  ifstream ReadWLSScint;

  G4String wls_Scint_file="../input_files/full_popop_emission.cfg";
  ReadWLSScint.open(wls_Scint_file);

  if(ReadWLSScint.is_open())
  {
      while(!ReadWLSScint.eof())
      {
          ReadWLSScint>>wavelength>>filler>>wlsEmit[500-wlsScintEntries];
          //convert wavelength to eV:
          wlsEnergy[500 - wlsScintEntries] = (1240/wavelength)*eV;
          wlsScintEntries++;
      }
  }else G4cout<<"Error opening file: "<<wls_Scint_file<<G4endl;
  ReadWLSScint.close();
  wlsScintEntries--;

  G4int wlsAbsorbEntries = 0;
  G4double wlsAbsorbEnergy[501] = {0};
  G4double wlsAbsorb[501] = {0};

  ifstream ReadWLSAbsorb;
  G4String ReadWLSAbsorbLength="../input_files/scintAbsLen.txt";

  ReadWLSAbsorb.open(ReadWLSAbsorbLength);
  if (ReadWLSAbsorb.is_open())
  {
      while(!ReadWLSAbsorb.eof())
      {
          ReadWLSAbsorb>>wavelength>>filler>>varabsorblength;
          wlsAbsorbEnergy[500 - wlsAbsorbEntries]=(1240/wavelength)*eV;
          wlsAbsorb[500 - wlsAbsorbEntries]=varabsorblength*m;
          wlsAbsorbEntries++;
      }
  }else G4cout<<"Error opening file: "<<ReadWLSAbsorb<<G4endl;
  ReadWLSAbsorb.close();
  wlsAbsorbEntries--;

  G4MaterialPropertiesTable* MPT = new G4MaterialPropertiesTable();

  MPT->AddProperty("WLSABSLENGTH",wlsAbsorbEnergy,wlsAbsorb,wlsAbsorbEntries);
  MPT->AddProperty("WLSCOMPONENT",wlsEnergy,wlsEmit,wlsScintEntries);
  MPT->AddConstProperty("WLSTIMECONSTANT", 12*ns);

  MPT->AddProperty("RINDEX",        rindexEnergy,  scintIndex, rindexEntries);
  MPT->AddProperty("ABSLENGTH",     absorbEnergy, Absorb,     absorbEntries);
  MPT->AddProperty("FASTCOMPONENT", scintEnergy,  scintEmit,  scintEntries);
  MPT->AddProperty("SLOWCOMPONENT",scintEnergy, scintEmitSlow,     scintEntries);

  MPT->AddConstProperty("SCINTILLATIONYIELD",11520./MeV);
  //MPT->AddConstProperty("SCINTILLATIONYIELD",100./MeV);
  MPT->AddConstProperty("RESOLUTIONSCALE",4.0);
  MPT->AddConstProperty("FASTTIMECONSTANT", 2.1*ns);
  MPT->AddConstProperty("SLOWTIMECONSTANT",14.2*ns);
  MPT->AddConstProperty("YIELDRATIO",1.0);

  fScintilator->SetMaterialPropertiesTable(MPT);
}

void DetectorConstruction::SetVolName(G4ThreeVector thePoint){
  G4Navigator* theNavigator = G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();
  G4VPhysicalVolume* myVolume= theNavigator->LocateGlobalPointAndSetup(thePoint);
  fVolName =  myVolume->GetName();
}

void DetectorConstruction::SetPropertyTable(G4Material* material, G4MaterialPropertiesTable* table){
  material->SetMaterialPropertiesTable(table);
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

  fWLBox = new G4LogicalVolume(fWorldBox,fWorldMaterial,"World",0,0,0);

  fWPBox = new G4PVPlacement(0,G4ThreeVector(),fWLBox,"World",0,false,0);

  double thickness_foil = 2.5*mm;
  double holeHeight = fThickness;
  G4Box* absorberPb_box = new G4Box("absorber",15*mm,15*mm, 1.5*mm);
  G4Tubs* hole = new G4Tubs("hole",0*mm,1.5*mm,holeHeight, 0.*deg, 360.*deg);
  G4SubtractionSolid* collimator = new G4SubtractionSolid("collimator",absorberPb_box,hole);

  G4Box* absorber_foil = new G4Box("foil",15*mm,15*mm,thickness_foil);
  G4Tubs* hole_foil = new G4Tubs("hole_foil",0*mm,1.5*mm,thickness_foil, 0.*deg, 360.*deg);
  G4SubtractionSolid* collimator_foil = new G4SubtractionSolid("collimator_foil",absorber_foil,hole_foil);

  double target_width = 15*mm;
  fBox = new G4Box("target", target_width, target_width, fTargetThickness);
  fLBox = new G4LogicalVolume(fBox,fScintilator, "target",0,0,0);
  double position = 0;

  double rod_sides = 2.5*mm;
  double rod_length = 48*mm;

  G4Box* rod = new G4Box("target",rod_sides,rod_sides,rod_length);
  G4LogicalVolume* rod_log = new G4LogicalVolume(rod,fPEN,"target",0,0,0);

  G4NistManager* man = G4NistManager::Instance();

  // SourceHolder
  G4double outer_holder_width1 = 12.5*mm;
  G4double outer_holder_width2 = 15*mm;
  G4double outer_pom_thickness = 2*mm;
  G4double outer_al_thickness = 0.5*mm;
  G4double inner_holder_width1 = 5.5*mm;
  G4double inner_holder_width2 = 11.5*mm;
  G4double inner_pom_thickness = 1*mm;
  G4RotationMatrix* rotm = new G4RotationMatrix();

  G4Box* outer_box = new G4Box("outer_box",outer_holder_width1,outer_holder_width2,outer_pom_thickness);
  G4Box* inner_box = new G4Box("inner_box",inner_holder_width1,inner_holder_width2,inner_pom_thickness);
  G4ThreeVector box_offset = G4ThreeVector(0,0,1*mm);
  G4SubtractionSolid* pom_holder = new G4SubtractionSolid("pom_holder",outer_box,inner_box,rotm,box_offset);
  G4LogicalVolume* pom_holder_log = new G4LogicalVolume(pom_holder,fABS,"pom_source_holder");

  G4Box* al_box = new G4Box("al_box",outer_holder_width1,outer_holder_width2,outer_al_thickness);
  G4LogicalVolume* al_box_log = new G4LogicalVolume(al_box, fAl, "al_box_log");

  // Breadboard

  G4double bread_board_width1 = 50*mm;
  G4double bread_board_width2 = 100*mm;
  G4double bread_board_thickness = 2.5*mm;

  G4Box* breadboard_box = new G4Box("breadboard_box",bread_board_width1,bread_board_width2,bread_board_thickness);
  G4LogicalVolume* breadboard_log = new G4LogicalVolume(breadboard_box,fAl, "breadboard_log");

  // Table
  G4double table_width1 = 0.1*m;
  G4double table_width2 = 0.1*m;
  G4double table_thickness = 15*mm;

  G4Box* table_box = new G4Box("table_box",table_width1,table_width2,table_thickness);
  G4LogicalVolume* table_log = new G4LogicalVolume(table_box, fWood,"table_log");

  // Silicon Plates

  fSiliconPlate_h = 1.5*mm;
  fHolderWidth = 90.00*mm;

  G4double trg_b = (fHolderWidth/sqrt(2.))*mm;
  G4double trg_h = 45.*mm;
  G4double rect_x = 2.*(trg_h + 9.5*mm);
  G4double rect_y = fHolderWidth;
  G4double rect_bc_x = 28.*mm;
  G4double cut_x = 20.*mm;
  G4double cut_y = 31.*mm;

  SiliconPlateConstruction plate;
  G4VSolid* final_plate = plate.ConstructPlate();
  G4LogicalVolume* plate_log = new G4LogicalVolume(final_plate,fPEN,"plate");
  G4cout <<  G4BestUnit(plate_log->GetMass(true),"Mass") << G4endl;
  G4ThreeVector point = G4ThreeVector(0,0,5*cm);

  // Dummy plate for side detector

  G4VSolid* plate_side_det = new G4Box("side_det",rect_x/2,trg_h+1*mm,0.45*fSiliconPlate_h);
  G4VSolid* side_det = new G4SubtractionSolid("side_det",plate_side_det,final_plate);

  G4Navigator* pointNavigator = new G4Navigator();
  pointNavigator->SetWorldVolume(fWPBox);
  pointNavigator->LocateGlobalPointAndSetup(point);

  // --------------Detectors--------------

  G4VSolid* side_sphere = new G4Sphere("side_sphere",70*mm,72*mm,0*deg,360*deg,0*deg, 360*deg);
  G4LogicalVolume* sideSphereLog = new G4LogicalVolume(side_sphere, fWorldMaterial,"side_sphere");

  // PMT_Round
  char filler;
  G4double wavelength;
  G4double innerRadius_pmt = 0.*cm;
  G4double outerRadius_pmt = 53*mm;
  G4double outerRadius_cath = 46*mm;
  G4double height_pmt = 63.5*mm;
  G4double height_cath = 62.*mm;
  G4double startAngle_pmt = 0.*deg;
  G4double spanningAngle_pmt = 360.*deg;

  G4Tubs* pmt = new G4Tubs("pmt_tube",innerRadius_pmt,outerRadius_pmt, height_pmt,startAngle_pmt,spanningAngle_pmt);
  G4LogicalVolume* pmt_log = new G4LogicalVolume(pmt,fGlass, "pmt_log");
  G4Tubs* Photocath = new G4Tubs("photocath_tube",innerRadius_pmt,outerRadius_cath,
                          height_cath,startAngle_pmt,spanningAngle_pmt);
  G4LogicalVolume* photocath_log = new G4LogicalVolume(Photocath, fAl,"photocath_log");

  // PMT_Square
  G4double casing_length = 30*mm/2;
  G4double casing_height = 32.5*mm/2;
  G4double casing_hole = 29*mm / 2;
  G4double casing_hole_height = 31.5*mm/2;

  G4double photo_cath_length = 26*mm/2;
  G4double photo_cath_eff_length = 23*mm/2;
  G4double photo_cath_height = 0.8*mm;

  G4Box* pmt_case = new G4Box("case",casing_length,casing_length,casing_height);
  G4Box* pHole = new G4Box("hole",casing_hole,casing_hole,casing_hole_height);
  G4SubtractionSolid* pmt_shell = new G4SubtractionSolid("pmt_shell",pmt_case,pHole);

  G4LogicalVolume* pmt_void = new G4LogicalVolume(pHole, fVacuum,"vacuum");
  G4LogicalVolume* pmt_case_log = new G4LogicalVolume(pmt_shell,fPOM,"pmt_case_log");

  G4Box* pmt_cath = new G4Box("pmt_cathode",photo_cath_length,photo_cath_length,photo_cath_height);
  G4Box* pmt_active = new G4Box("pmt_active",photo_cath_eff_length,photo_cath_eff_length,photo_cath_height);
  G4SubtractionSolid* pmt_inactive_cath = new G4SubtractionSolid("pmt_inactive_cathode",pmt_cath,pmt_active);
  G4LogicalVolume* pmt_inactive_cath_log = new G4LogicalVolume(pmt_inactive_cath,fGlass,"pmt_inactive_cath_log");
  G4LogicalVolume* pmt_cath_log = new G4LogicalVolume(pmt_active,fGlass,"pmt_cath_log");

  G4double cath_eff;
  G4double photocath_energy[57];
  G4double photocath_EFF[57];
  G4double perfect_EFF[57];
  G4double perfect_REFL[57];
  G4double photocath_REFL[57]={0};
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
      photocath_EFF[57-effCounter] = cath_eff;
      perfect_EFF[57-effCounter] = 1;
      perfect_REFL[57-effCounter] = 0;
      effCounter++;
    }
  }

  else G4cout<<"Error opening file: " <<pmt_file<<G4endl;
  ReadEff.close();
  effCounter--;

  const G4int nPMT_EFF = sizeof(photocath_energy)/sizeof(G4double);

  G4OpticalSurface* photocath_optsurf = new G4OpticalSurface("photocath_opsurf",glisur,polished, dielectric_metal);
  G4MaterialPropertiesTable* photocath_MT = new G4MaterialPropertiesTable();
  photocath_MT->AddProperty("EFFICIENCY", photocath_energy, photocath_EFF,nPMT_EFF);
  photocath_MT->AddProperty("REFLECTIVITY", photocath_energy, photocath_REFL,nPMT_EFF);
  photocath_optsurf->SetMaterialPropertiesTable(photocath_MT);
  new G4LogicalSkinSurface("photocath_surf",photocath_log,photocath_optsurf);
  new G4LogicalSkinSurface("photocath_surf",pmt_cath_log,photocath_optsurf);

  G4OpticalSurface* perfect_optsurf = new G4OpticalSurface("perfect",glisur,polished, dielectric_metal);
  G4MaterialPropertiesTable* detector_MT = new G4MaterialPropertiesTable();
  detector_MT->AddProperty("EFFICIENCY", photocath_energy, perfect_EFF,nPMT_EFF);
  detector_MT->AddProperty("REFLECTIVITY", photocath_energy, perfect_REFL,nPMT_EFF);
  perfect_optsurf->SetMaterialPropertiesTable(detector_MT);

  G4LogicalVolume* tile_detector = new G4LogicalVolume(fBox, fSi, "tile_detector");
  new G4LogicalSkinSurface("perfect_sensor",tile_detector, perfect_optsurf);

  // Spectrometer Sensor

  G4double spectrometer_side = 2.5*mm;
  G4double spectrometer_top = 1.5*mm;
  G4double spectrometer_thickness = 0.5*mm;

  G4Box* spectrometer_sensor = new G4Box("spec_sensor",spectrometer_top,spectrometer_side,spectrometer_thickness);
  G4LogicalVolume* spec_log = new G4LogicalVolume(spectrometer_sensor,fSi,"spec_sensor_log");

  G4double spec_energy[57];
  G4double spec_EFF[57];
  G4double spec_REFL[57]={0};
  G4String spec_file = "../input_files/specQE.csv";

//  ifstream ReadEff;
  //G4int
  effCounter = 0;
  ReadEff.open(spec_file);

  if(ReadEff.is_open())
  {
    while(!ReadEff.eof())
    {
      ReadEff>>wavelength>>filler>>cath_eff;
      if(ReadEff.eof()){
        break;
      }
      spec_energy[57-effCounter] = (1240/wavelength)*eV;
      spec_EFF[57-effCounter] = cath_eff;
      effCounter++;
    }
  }

  else G4cout<<"Error opening file: " <<spec_file<<G4endl;
  ReadEff.close();
  effCounter--;

  const G4int nSPEC_EFF = sizeof(spec_energy)/sizeof(G4double);

  G4OpticalSurface* spec_optsurf = new G4OpticalSurface("photocath_opsurf",glisur,polished, dielectric_metal);
  G4MaterialPropertiesTable* spec_MT = new G4MaterialPropertiesTable();
  spec_MT->AddProperty("EFFICIENCY", spec_energy, spec_EFF,nSPEC_EFF);
  spec_MT->AddProperty("REFLECTIVITY", spec_energy, spec_REFL,nSPEC_EFF);
  spec_optsurf->SetMaterialPropertiesTable(spec_MT);
  new G4LogicalSkinSurface("spec_surf",spec_log,spec_optsurf);


  // SiPMs
  G4double SiPM_side = 1.5*mm;
  G4double SiPM_thickness = 0.225*mm;
  G4double board_size = 25*mm;

  G4double SiPM_case_side = 3.0*mm;
  G4double SiPM_case_thickness = 1*mm;

  G4Box* SiPM_window = new G4Box("SiPM_window",SiPM_side,SiPM_side,SiPM_thickness);
  G4Box* board_window = new G4Box("board_window",board_size,board_size,SiPM_case_thickness);
  G4Box* SiPM_case = new G4Box("SiPM_window",SiPM_case_side,SiPM_case_side,SiPM_case_thickness);
  G4SubtractionSolid* SiPM_case_full = new G4SubtractionSolid("SiPM_case_full",SiPM_case, SiPM_window,0,G4ThreeVector(0,0,SiPM_case_thickness-SiPM_thickness));

  G4LogicalVolume* SiPM_window_log = new G4LogicalVolume(SiPM_window,fSi,"SiPM_window_log");
  G4LogicalVolume* SiPM_case_log = new G4LogicalVolume(SiPM_case_full,fSi,"SiPM_case_log");
  G4LogicalVolume* SiPM_board_log = new G4LogicalVolume(board_window,fSi,"board_log");

  G4double fill_factor = 0.74;

  G4double SiPM_eff;
  G4double SiPM_energy[57];
  G4double SiPM_EFF[57];
  G4double SiPM_REFL[57]={0};
  G4String siPM_file = "../input_files/sipmQE.csv";

  ifstream ReadEffSiPM;
  G4int effCounterSiPM = 0;
  ReadEffSiPM.open(siPM_file);

  if(ReadEffSiPM.is_open())
  {
    while(!ReadEffSiPM.eof())
    {
      ReadEffSiPM>>wavelength>>filler>>SiPM_eff;
      if(ReadEffSiPM.eof()){
        break;
      }
      SiPM_energy[57-effCounterSiPM] = (1240/wavelength)*eV;
      SiPM_EFF[57-effCounterSiPM] = SiPM_eff/100 * fill_factor;
      effCounterSiPM++;
    }
  }

  else G4cout<<"Error opening file: " <<siPM_file<<G4endl;
  ReadEffSiPM.close();
  effCounterSiPM--;

  const G4int nSiPM_EFF = sizeof(SiPM_energy)/sizeof(G4double);

  G4OpticalSurface* SiPM_optsurf = new G4OpticalSurface("SiPM_opsurf",glisur,polished, dielectric_metal);
  G4MaterialPropertiesTable* SiPM_MT = new G4MaterialPropertiesTable();
  SiPM_MT->AddProperty("EFFICIENCY", SiPM_energy, SiPM_EFF,nSiPM_EFF);
  SiPM_MT->AddProperty("REFLECTIVITY", SiPM_energy, SiPM_REFL,nSiPM_EFF);
  SiPM_optsurf->SetMaterialPropertiesTable(SiPM_MT);
  new G4LogicalSkinSurface("sipm_surf",SiPM_window_log ,SiPM_optsurf);

  // PMT Mount
  G4double mount_base = 55*mm;
  G4double pmt_neg_length = 55.1*mm;
  G4double mount_height = 26*mm;
  G4double mount_spacer = 9*mm;

  G4Box* mount_block = new G4Box("mount_box",mount_base,mount_base,mount_height);
  G4Box* pmt_neg = new G4Box("pmt_neg",casing_length,casing_length,pmt_neg_length);
  G4Box* pmt_base_neg = new G4Box("pmt_base_neg",casing_length,casing_length,mount_spacer);

  G4SubtractionSolid* mount_box = new G4SubtractionSolid("mount_box",mount_block,pmt_base_neg,rotm,G4ThreeVector(0,0,-9*mm));
  rotm->rotateY(90.*deg);
  mount_box = new G4SubtractionSolid("mount_box",mount_box,pmt_neg,rotm,G4ThreeVector(0*mm,0,18*mm));
  rotm->rotateY(-90.*deg);
  rotm->rotateX(90.*deg);
  mount_box = new G4SubtractionSolid("mount_box",mount_box,pmt_neg,rotm,G4ThreeVector(0*mm,0,18*mm));

  G4LogicalVolume* mount_box_log = new G4LogicalVolume(mount_box,fABS, "mount_box_log");

  if (fDetectorType == 2){
      parser.Read("geCounter.gdml",false);
      fWPBox = parser.GetWorldVolume();
    }

  // Placement
  G4LogicalVolume* topDet = new G4LogicalVolume(final_plate,fWorldMaterial,"top_det");
  G4LogicalVolume* bottomDet = new G4LogicalVolume(final_plate,fWorldMaterial,"bottom_det");
  G4LogicalVolume* sideDet = new G4LogicalVolume(side_det,fWorldMaterial,"side_det");
  new G4LogicalSkinSurface("perfect_surf",topDet,perfect_optsurf);
  new G4LogicalSkinSurface("perfect_surf",bottomDet,perfect_optsurf);
  new G4LogicalSkinSurface("perfect_surf",sideDet,perfect_optsurf);
  new G4LogicalSkinSurface("perfect_surf",sideSphereLog,perfect_optsurf);
  new G4LogicalSkinSurface("perfect_surf",tile_detector,perfect_optsurf);

  G4VPhysicalVolume* topDetPlacement;
  G4VPhysicalVolume* bottomDetPlacement;
  G4VPhysicalVolume* sideDetPlacement;

  G4VPhysicalVolume* tablePlacement;
  G4VPhysicalVolume* breadboardPlacement;

  G4VPhysicalVolume* pmtPlacement;
  G4VPhysicalVolume* incathPlacement;
  G4VPhysicalVolume* cathPlacement;

  G4VPhysicalVolume* pmtPlacement1;
  G4VPhysicalVolume* incathPlacement1;
  G4VPhysicalVolume* cathPlacement1;

  G4VPhysicalVolume* pmtPlacement2;
  G4VPhysicalVolume* incathPlacement2;
  G4VPhysicalVolume* cathPlacement2;

  G4VPhysicalVolume* pmtPlacement3;
  G4VPhysicalVolume* incathPlacement3;
  G4VPhysicalVolume* cathPlacement3;

  G4VPhysicalVolume* pmtPlacement4;
  G4VPhysicalVolume* incathPlacement4;
  G4VPhysicalVolume* cathPlacement4;

  G4VPhysicalVolume* pmtPlacement5;
  G4VPhysicalVolume* incathPlacement5;
  G4VPhysicalVolume* cathPlacement5;

  G4VPhysicalVolume* collimator_phys;
  G4VPhysicalVolume* al_phys;
  G4VPhysicalVolume* collimator_phys1;
  G4VPhysicalVolume* al_phys1;
  G4VPhysicalVolume* pmt_holder_phys;

  G4VPhysicalVolume* siPM_placement;
  G4VPhysicalVolume* siPM_case_placement;
  G4VPhysicalVolume* siPM_board_placement;

  G4VPhysicalVolume* vacPlacement;

  G4double fullPMT_height = 2*casing_height+2*photo_cath_height;
  G4ThreeVector* pmtVector = new G4ThreeVector(0,0,-(casing_height+2*photo_cath_height+fTargetThickness));
  fSourceVector = new G4ThreeVector(0,0,(casing_height+2*photo_cath_height+fTargetThickness));

  // Set Draw G4VisAttributes

  G4VisAttributes* visAttr = new G4VisAttributes();
  visAttr->SetVisibility(false);
  fWLBox->SetVisAttributes(visAttr);
  pmt_void->SetVisAttributes(visAttr);
  // Scintillator Materials
  G4VisAttributes* tileAttr = new G4VisAttributes(G4Colour::Blue());
  tileAttr->SetVisibility(true);
  fLBox->SetVisAttributes(tileAttr);


  // // Active Detectors
   G4VisAttributes* detectorAttr = new G4VisAttributes(G4Colour::Green());
   detectorAttr->SetVisibility(true);
   detectorAttr->SetForceSolid(false);
   pmt_cath_log->SetVisAttributes(detectorAttr);
   SiPM_window_log->SetVisAttributes(detectorAttr);
   topDet->SetVisAttributes(detectorAttr);
   bottomDet->SetVisAttributes(detectorAttr);
   sideDet->SetVisAttributes(detectorAttr);
   sideSphereLog->SetVisAttributes(detectorAttr);

  // // Inactive Regions
  G4VisAttributes* supportAttr = new G4VisAttributes(G4Colour::Grey());
  supportAttr->SetVisibility(true);
  supportAttr->SetForceSolid(true);
  pmt_inactive_cath_log->SetVisAttributes(supportAttr);

  // Surrounding Materials

  G4VisAttributes* surroundAttr = new G4VisAttributes(G4Colour::White());
  surroundAttr->SetVisibility(true);
  table_log->SetVisAttributes(surroundAttr);
  breadboard_log->SetVisAttributes(surroundAttr);
  pmt_case_log->SetVisAttributes(surroundAttr);

  G4VisAttributes* mountAttr = new G4VisAttributes(G4Colour::Red());
  mountAttr->SetVisibility(true);
  mount_box_log->SetVisAttributes(mountAttr);

  // Active Shielding
  G4VisAttributes* shieldingAttr = new G4VisAttributes(G4Colour::Yellow());
  shieldingAttr->SetVisibility(true);
  pom_holder_log->SetVisAttributes(shieldingAttr);
  al_box_log->SetVisAttributes(shieldingAttr);

  G4RotationMatrix* rotationMatrix = new G4RotationMatrix(0,0,0);
  G4RotationMatrix* rotationMatrix1 = new G4RotationMatrix(0,0,0);
  G4RotationMatrix* rotationMatrix2 = new G4RotationMatrix(0,0,0);
  G4RotationMatrix* rotationMatrix3 = new G4RotationMatrix(0,0,0);
  G4RotationMatrix* rotationMatrix4 = new G4RotationMatrix(0,0,0);
  G4RotationMatrix* rotationMatrix5 = new G4RotationMatrix(0,0,0);
  /*
  0 - PMT on base of tile, collimator included.
  */
  fDetectorType = 0;
  switch (fDetectorType){
    case 0:
       fPBox = new G4PVPlacement(0, G4ThreeVector(0,0,0),fLBox,"rod",fWLBox,false,0,true);
       siPM_placement = new G4PVPlacement(0, G4ThreeVector(0,0,-6.001*mm),tile_detector,"spec",fWLBox,false,0,true);


       // bottomDetPlacement = new G4PVPlacement(0, G4ThreeVector(0,0,1.6*mm),bottomDet,"bottom_det",fWLBox,false,0,true);
       // topDetPlacement = new G4PVPlacement(0, G4ThreeVector(0,0,-1.6*mm),topDet,"top_det",fWLBox,false,0,true);
       //
       //
       // sideDetPlacement = new G4PVPlacement(0, G4ThreeVector(0,0,0),sideSphereLog,"side_det",fWLBox,false,0,true);
// //   Place PMT
       // pmtPlacement = new G4PVPlacement(0,*pmtVector,pmt_case_log,"pmt",fWLBox,false,0,true);
       // cathPlacement = new G4PVPlacement(0,G4ThreeVector(0,0,-(fTargetThickness+photo_cath_height)),pmt_cath_log,"active_detector",fWLBox,false,0,true);
       // incathPlacement = new G4PVPlacement(0,G4ThreeVector(0,0,-(fTargetThickness+photo_cath_height)),pmt_inactive_cath_log,"inactive_detector",fWLBox,false,0,true);
///home/iwsatlas1/hayward/Documents/LOSim/
        // vacPlacement = new G4PVPlacement(0,G4ThreeVector(),pmt_void,"vacuum",pmt_case_log,false,0);
//
//   //  Place table and breadboard
    //   breadboardPlacement = new G4PVPlacement(0,G4ThreeVector(0,0,-(fTargetThickness+bread_board_thickness+fullPMT_height+4*cm)),breadboard_log,"breadboard",fWLBox,false,0);
    //   tablePlacement = new G4PVPlacement(0,G4ThreeVector(0,0,-(fTargetThickness+2*bread_board_thickness+table_thickness+fullPMT_height+4*cm)),table_log,"table",fWLBox,false,0);
//
//   // Place source holder
    //    collimator_phys = new G4PVPlacement(0,G4ThreeVector(0,0,58*mm),pom_holder_log,"collimator",fWLBox, false, 0,true);
    //    al_phys = new G4PVPlacement(0,G4ThreeVector(0,0,56*mm-outer_al_thickness),al_box_log,"al_collimator",fWLBox,false,0,true);

    //    rotationMatrix5->rotateX(180*deg);
    //    collimator_phys1 = new G4PVPlacement(rotationMatrix5,G4ThreeVector(0,0,62*mm),pom_holder_log,"collimator1",fWLBox, false, 0,true);
    //    al_phys1 = new G4PVPlacement(rotationMatrix5,G4ThreeVector(0,0,64*mm+outer_al_thickness),al_box_log,"al_collimator1",fWLBox,false,0,true);

  //   Place PMT support structure
    //    pmt_holder_phys = new G4PVPlacement(0,G4ThreeVector(0,0,-20*mm),mount_box_log,"holder",fWLBox, false, 0,true);

          break;
    case 1:
         fPBox = new G4PVPlacement(0, G4ThreeVector(0,0,position),fLBox,"target",fWLBox,false,0,true);
    //   Place PMT
         pmtPlacement = new G4PVPlacement(0,*pmtVector,pmt_case_log,"pmt",fWLBox,false,0,true);
         cathPlacement = new G4PVPlacement(0,G4ThreeVector(0,0,-(fTargetThickness+photo_cath_height)),pmt_cath_log,"active_detector",fWLBox,false,0,true);
         incathPlacement = new G4PVPlacement(0,G4ThreeVector(0,0,-(fTargetThickness+photo_cath_height)),pmt_inactive_cath_log,"inactive_detector",fWLBox,false,0,true);

         rotationMatrix->rotateY(90.*deg);
         rotationMatrix1->rotateY(90.*deg);
         pmtPlacement1 = new G4PVPlacement(rotationMatrix,G4ThreeVector(32.9*mm,0,0*mm),pmt_case_log,"pmt1",fWLBox,false,0,true);
         cathPlacement1 = new G4PVPlacement(rotationMatrix1,G4ThreeVector(33.5*mm-casing_height,0,0*mm),pmt_cath_log,"active_detector1",fWLBox,false,0,true);
         incathPlacement1 = new G4PVPlacement(rotationMatrix1,G4ThreeVector(33.5*mm-casing_height,0,0*mm),pmt_inactive_cath_log,"inactive_detector1",fWLBox,false,0,true);

         rotationMatrix2->rotateY(-90.*deg);
         pmtPlacement2 = new G4PVPlacement(rotationMatrix,G4ThreeVector(-32.9*mm,0,0*mm),pmt_case_log,"pmt2",fWLBox,false,0,true);
         cathPlacement2 = new G4PVPlacement(rotationMatrix2,G4ThreeVector(-(33.5*mm-casing_height),0,0*mm),pmt_cath_log,"active_detector2",fWLBox,false,0,true);
         incathPlacement2 = new G4PVPlacement(rotationMatrix2,G4ThreeVector(-(33.5*mm-casing_height),0,0*mm),pmt_inactive_cath_log,"inactive_detector2",fWLBox,false,0,true);

         rotationMatrix->rotateX(90.*deg);
         rotationMatrix->rotateY(90.*deg);

         rotationMatrix3->rotateX(90.*deg);
         pmtPlacement3 = new G4PVPlacement(rotationMatrix,G4ThreeVector(0,32.9*mm,0*mm),pmt_case_log,"pmt3",fWLBox,false,0,true);
         cathPlacement3 = new G4PVPlacement(rotationMatrix3,G4ThreeVector(0,-(33.5*mm-casing_height),0*mm),pmt_cath_log,"active_detector3",fWLBox,false,0,true);
         incathPlacement3 = new G4PVPlacement(rotationMatrix3,G4ThreeVector(0,-(33.5*mm-casing_height),0*mm),pmt_inactive_cath_log,"inactive_detector3",fWLBox,false,0,true);

         rotationMatrix4->rotateX(-90.*deg);
         pmtPlacement4 = new G4PVPlacement(rotationMatrix,G4ThreeVector(0,-32.9*mm,0*mm),pmt_case_log,"pmt4",fWLBox,false,0,true);
         cathPlacement4 = new G4PVPlacement(rotationMatrix4,G4ThreeVector(0,(33.5*mm-casing_height),0*mm),pmt_cath_log,"active_detector4",fWLBox,false,0,true);
         incathPlacement4 = new G4PVPlacement(rotationMatrix4,G4ThreeVector(0,(33.5*mm-casing_height),0*mm),pmt_inactive_cath_log,"inactive_detector4",fWLBox,false,0,true);

         pmtPlacement5 = new G4PVPlacement(0,*pmtVector*-1,pmt_case_log,"pmt5",fWLBox,false,0,true);
         cathPlacement5 = new G4PVPlacement(0,G4ThreeVector(0,0,(fTargetThickness+photo_cath_height)),pmt_cath_log,"active_detector5",fWLBox,false,0,true);
         incathPlacement5 = new G4PVPlacement(0,G4ThreeVector(0,0,(fTargetThickness+photo_cath_height)),pmt_inactive_cath_log,"inactive_detector5",fWLBox,false,0,true);

         vacPlacement = new G4PVPlacement(0,G4ThreeVector(),pmt_void,"vacuum",pmt_case_log,false,0);

    //  Place table and breadboard
         breadboardPlacement = new G4PVPlacement(0,G4ThreeVector(0,0,-(fTargetThickness+bread_board_thickness+fullPMT_height+3*cm)),breadboard_log,"breadboard",fWLBox,false,0);
         tablePlacement = new G4PVPlacement(0,G4ThreeVector(0,0,-(fTargetThickness+2*bread_board_thickness+table_thickness+fullPMT_height+3*cm)),table_log,"table",fWLBox,false,0);

         // Place source holder
          collimator_phys = new G4PVPlacement(0,G4ThreeVector(0,0,58*mm),pom_holder_log,"collimator",fWLBox, false, 0,true);
          al_phys = new G4PVPlacement(0,G4ThreeVector(0,0,56*mm-outer_al_thickness),al_box_log,"al_collimator",fWLBox,false,0,true);

          rotationMatrix5->rotateX(180*deg);
          collimator_phys1 = new G4PVPlacement(rotationMatrix5,G4ThreeVector(0,0,62*mm),pom_holder_log,"collimator1",fWLBox, false, 0,true);
          al_phys1 = new G4PVPlacement(rotationMatrix5,G4ThreeVector(0,0,64*mm+outer_al_thickness),al_box_log,"al_collimator1",fWLBox,false,0,true);

    //   Place PMT support structure
        pmt_holder_phys = new G4PVPlacement(0,G4ThreeVector(0,0,-19*mm),mount_box_log,"holder",fWLBox, false, 0,true);
      break;

  }
  // G4GDMLParser* parser1 = new G4GDMLParser();
  // w = G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking()->GetWorldVolume();

  return fWPBox;
}
