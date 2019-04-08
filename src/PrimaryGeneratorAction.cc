#include "PrimaryGeneratorAction.hh"
#include "PrimaryGeneratorMessenger.hh"

#include "Randomize.hh"
#include "DetectorConstruction.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "G4IonTable.hh"
#include "G4DecayTable.hh"
#include "G4OpticalPhoton.hh"
#include "G4VPhysicalVolume.hh"

#include "G4Navigator.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction(DetectorConstruction* det)
: G4VUserPrimaryGeneratorAction(),
fParticleGun(0),
fDetector(det)
{
	G4int n_particle = 1;
	fPrimaryMessenger = new PrimaryGeneratorMessenger(this);
	fParticleGun = new G4ParticleGun(n_particle);

 	//default kinematic
 	//
	G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
	G4ParticleDefinition* particle = particleTable->FindParticle("e-");
	fSourceType = 4;
	fSourceEnergy = 1500*keV;
	fPhotonWavelength = 0;
	fParticleName = "void";
	fPoint = G4ThreeVector();
	DefineParticle();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
	delete fParticleGun;
}

void PrimaryGeneratorAction::DefineParticle(){

	// Particle Types
	// 	0 - 137Cs source
	//  1 - Perpendicular gamma with fixed position
	// 	2 - Perpendicular gamma with random position
	//	3 - Gamma generated at random point inside named volume

	G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();

  G4double ionCharge   = 0.*eplus;
  G4double excitEnergy = 0.*keV;
	G4int Z=0, A=0;
	G4double x,y,z, sum;
	G4ParticleDefinition* ion;
	fSourceType = 4;
	G4ThreeVector position = G4ThreeVector(0,0,60*mm);
	G4String name;
	float rx,ry,rz;

	fParticleGun->SetParticleDefinition(G4OpticalPhoton::Definition());
	fParticleGun->SetParticlePolarization(G4ThreeVector(2*(rand()-1),2*(rand()-1),2*(rand()-1)));
	fParticleGun->SetParticleEnergy(3.262742*eV);
	fParticleGun->SetParticleMomentumDirection(G4ThreeVector(1,0,0));

	G4double radiusX = 1.5*mm;
	G4double radiusY = 1.0*mm;

	rx = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
	ry = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
	x = (rx)*2*3.14159265;
	y= (ry)*2*3.14159265;
	G4double rradiusX = rx*radiusX;
	G4double rradiusY = ry*radiusY;

	fParticleGun->SetParticlePosition(G4ThreeVector(rradiusX*cos(y),rradiusY*sin(y),15*cm));

	G4cout<<"---------------------------------------------------------------------------------------"<<G4endl;
	//SetParticleName(Z,A,excitEnergy);
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Randomises placement and momentum vectors for required sources.

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{

	fSourceType = 4;
	G4String name;
	double x,y,z,sum;
	float rx,ry,rz;
	fParticleGun->SetParticleDefinition(G4OpticalPhoton::Definition());
	fParticleGun->SetParticlePolarization(G4ThreeVector(2*(rand()-1),2*(rand()-1),2*(rand()-1)));
	fParticleGun->SetParticleEnergy(3.262742*eV);
	fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0,0,-1));

	G4double radiusX = 1.5*mm;
	G4double radiusY = 1.0*mm;

	rx = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
	ry = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
	x = (rx)*2*3.14159265;
	y = (ry)*2*3.14159265;
	G4double rradiusX = rx*radiusX;
	G4double rradiusY = ry*radiusY;

	fParticleGun->SetParticlePosition(G4ThreeVector(rradiusX*cos(y),rradiusY*sin(y),5.5*cm));
	fParticleGun->GeneratePrimaryVertex(anEvent);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::SetParticleName(G4int Z, G4int A, G4double excitEnergy)
{
	fParticleName = G4IonTable::GetIonTable()->GetIonName(Z,A,excitEnergy);
}

void PrimaryGeneratorAction::SetPosition(G4double newPosition){
	fPosition = newPosition;
	DefineParticle();
}

void PrimaryGeneratorAction::SetSourceType(G4int newType)
{
	if (newType <= 12 && newType >= 0)
	{
		fSourceType = newType;
	}
	else
	{
		G4cerr << "The option is out of the possible values (0-5)!" << G4endl;
		G4cerr << "The default option (0) is set!" << G4endl;
		fSourceType = 0;
	}
	DefineParticle();
}

void PrimaryGeneratorAction::SetSourceEnergy(G4double newEnergy)
{
	if (newEnergy>0)
	{
		fSourceEnergy = newEnergy;
	}
	else{
		G4cerr << "New energy is < 0." << G4endl;
		G4cerr << "The default option 60 keV is set!" << G4endl;
		fSourceEnergy = 60*keV;
	}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::SetPhotonWavelength(G4double newValue)
{
	if (newValue > 200 && newValue < 700)
	{
		fPhotonWavelength = newValue;
	}
	else
	{
		G4cerr << "The new desired wavelength is out of range (200-700 nm)!" << G4endl;
		G4cerr << "The photon wavelength is set to default value (420 nm)!" << G4endl;
		fPhotonWavelength = 420;
	}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
