//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "g4matrixPrimaryGeneratorAction.hh"
// #include "g4matrixPrimaryGeneratorMessenger.hh"

#include "Randomize.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "G4ChargedGeantino.hh"
#include "G4IonTable.hh"
//#include "CreateTree.hh"
#include <iostream>
#include <fstream>


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

g4matrixPrimaryGeneratorAction::g4matrixPrimaryGeneratorAction(ConfigFile& config)
: G4VUserPrimaryGeneratorAction(),
fParticleGun(0)
//    ,fConfig(config)
{
  G4int n_particle = 1;
  fParticleGun = new G4ParticleGun(n_particle);

  plasticx = config.read<double>("plasticx");
  lysox    = config.read<double>("lysox");
  ncouples = config.read<double>("ncouples");
  crystalx = config.read<double>("crystalx");
  crystaly = config.read<double>("crystaly");
  crystalz = config.read<double>("crystalz");
  greaseBack = config.read<double>("greaseBack");
  glassBack = config.read<double>("glassBack");
  airBack = config.read<double>("airBack");
  ncrystalx = config.read<int>("ncrystalx");
  ncrystaly = config.read<int>("ncrystaly");
  esrThickness = config.read<double>("esrThickness");

  //create a messenger for this class
  //fGunMessenger = new g4matrixPrimaryGeneratorMessenger(this);

  //default kinematic
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();

  backgroudSimulation = config.read<bool>("backgroudSimulation",0); // read from config file if this is a background simulation or not - by default it is not

  if(backgroudSimulation) // if background, use a chargedgeantino to simulate the decay of lu176
  {
    G4cout << "Running a background simulation" << G4endl;
    G4ParticleDefinition* particle = particleTable->FindParticle("chargedgeantino");
    fParticleGun->SetParticleDefinition(particle);
    fParticleGun->SetParticlePosition(G4ThreeVector(0.,0.,0.));
    fParticleGun->SetParticleEnergy(1*eV);
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(1.,0.,0.));
  }
  else
  {
    //--------------------------------------------------------//
    // GAMMA SOURCE                                           //
    //--------------------------------------------------------//
    //source
    //x position
    sourcex = config.read<double>("sourcex");
    //y position of source
    sourcey = config.read<double>("sourcey");
    //x position of source
    distance = config.read<double>("distance");
    //energy of gammas
    energy = config.read<double>("energy");
    //gamma direction
    direction = config.read<int>("direction");
    userIrradiationDiagonal = config.read<double>("irradiationDiagonal",0);

    G4double fakeAir = 0.1;  //fixed to 0.1mm, it has no physical meaning in our simulation
    sourcez = -(distance + (crystalz/2.0) + greaseBack + glassBack + airBack + fakeAir);

    G4cout << "Energy of gamma source [KeV]: " << energy << G4endl;
    G4cout << "Gamma source x position [mm]: " << sourcex << G4endl;
    G4cout << "Gamma source y position [mm]: " << sourcey << G4endl;
    G4cout << "Gamma source z position [mm]: " << sourcez << G4endl;
    G4cout << "Distance of source from back ESR [mm]: " << distance << G4endl;
    // source as a sphere
    // if sphereRadius > 0 is specified in the config file, the source is no more a point but a sphere centered in sourcex,sourcey,sourcez
    sphereRadius = config.read<double>("sphereRadius",0); //radius of the source sphere, in mm
    if(sphereRadius > 0)
    {
      G4cout << "Using a sphere as a gamma source" << G4endl;
      G4cout << "Sphere radius = " << sphereRadius << " mm" << G4endl;
    }

    if(direction == 0)
    {
      G4cout << "Gammas shoot parallel to the z axis" << G4endl;
    }
    else
    {
      G4cout << "Gammas shoot randomly towards the matrix" << G4endl;
    }
    //find the x and y of matrix

    G4ParticleDefinition* particle = particleTable->FindParticle("gamma");
    fParticleGun->SetParticleDefinition(particle);
    fParticleGun->SetParticleTime(0.0*ns);

    // save the coordinates of this source in the ttree file
    fParticleGun->SetParticlePosition(G4ThreeVector(0.,0.,0.));
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
    fParticleGun->SetParticleEnergy(energy*keV);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

g4matrixPrimaryGeneratorAction::~g4matrixPrimaryGeneratorAction()
{
  delete fParticleGun;
  //   delete fGunMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void g4matrixPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  float ScintillatorsLengthX = (crystalx + esrThickness) * ncrystalx;
  float ScintillatorsLengthY = (crystaly + esrThickness) * ncrystaly;
  float ScintillatorsLengthZ = crystalz;

  if(backgroudSimulation)
  {
    // G4cout << "----------------- O -----------------" << G4endl;
    G4ParticleDefinition* particle = fParticleGun->GetParticleDefinition();
    if (particle == G4ChargedGeantino::ChargedGeantino()) {
      //lu176
      G4int Z = 71, A = 176;
      G4double ionCharge   = 0.*eplus;
      G4double excitEnergy = 0.*keV;

      G4ParticleDefinition* ion
      = G4IonTable::GetIonTable()->GetIon(Z,A,excitEnergy);
      fParticleGun->SetParticleDefinition(ion);
      fParticleGun->SetParticleCharge(ionCharge);
    }
    EventTargetx = EventTargety = EventTargetz = 0;  // target does not makes sense in this configuration
    // randomized position
    // inside the matrix volume, but only in LYSO
    // generate random position inside the matrix volume,
    //generate position
    // bool isInLYSO = false;
    // while(!isInLYSO)
    // {
      EventSourcex = (G4UniformRand()*ScintillatorsLengthX) - (ScintillatorsLengthX/2.0);
      EventSourcey = (G4UniformRand()*ScintillatorsLengthY) - (ScintillatorsLengthY/2.0);
      EventSourcez = (G4UniformRand()*ScintillatorsLengthZ) - (ScintillatorsLengthZ/2.0);
      // check if it is in the space occupied by LYSO or not
      // EventSourcez is not a problem here

    // }
    //check that it's inside the LYSO material?
    G4cout << EventSourcex << " " << EventSourcex << " " << EventSourcez << G4endl;
    fParticleGun->SetParticlePosition(G4ThreeVector(EventSourcex,EventSourcey,EventSourcez));

    //create vertex
    //
    fParticleGun->GeneratePrimaryVertex(anEvent);
  }
  else
  {
    // source as a sphere
    // if sphereRadius > 0 is specified in the config file, the source is no more a point but a sphere centered in sourcex,sourcey,sourcez
    if(sphereRadius > 0)
    {
      // G4cout << "Using a sphere as a gamma source" << G4endl;
      // G4cout << "Sphere radius = " << sphereRadius << " mm" << G4endl;
      float sphereCenterX = sourcex;
      float sphereCenterY = sourcey;
      float sphereCenterZ = sourcez;
      //generate a random point in the sphere
      float alpha = G4UniformRand()*2*CLHEP::pi;
      float ran = G4UniformRand();
      float beta = acos(2.0*ran-1.0);
      float ranradius = cbrt(G4UniformRand());
      float EventRadius =  sphereRadius*ranradius;
      //now redefine the sourcex,sourcey,sourcez on the basis of this
      EventSourcex = sphereCenterX + EventRadius*cos(beta);
      EventSourcey = sphereCenterY + EventRadius*sin(beta)*cos(alpha);
      EventSourcez = sphereCenterZ + EventRadius*sin(beta)*sin(alpha);
      //find the target point in the matrix, to derive the direction vector
      //crystals are always centered in (0,0,0), so they cover a volume defined by they physical coordinates times their number
      // float ScintillatorsLengthX = (crystalx + esrThickness) * ncrystalx;
      // float ScintillatorsLengthY = (crystaly + esrThickness) * ncrystaly;
      // float ScintillatorsLengthZ = crystalz;
      EventTargetx = (G4UniformRand()*ScintillatorsLengthX) - (ScintillatorsLengthX/2.0);
      EventTargety = (G4UniformRand()*ScintillatorsLengthY) - (ScintillatorsLengthY/2.0);
      EventTargetz = (G4UniformRand()*ScintillatorsLengthZ) - (ScintillatorsLengthZ/2.0);
    }
    else
    {
      EventSourcex = sourcex;
      EventSourcey = sourcey;
      EventSourcez = sourcez;
    }

    //save source coordinates
    // G4cout << EventSourcex << " " << EventSourcey << " " << EventSourcez << G4endl;
    // CreateTree::Instance()->SourceX = (Float_t) EventSourcex;
    // CreateTree::Instance()->SourceY = (Float_t) EventSourcey;
    // CreateTree::Instance()->SourceZ = (Float_t) EventSourcez;



    //now direction
    G4ThreeVector directionVector;
    if(direction == 0)
    {
      directionVector = G4ThreeVector(0.,0.,1.);
    }
    else
    {
      if(sphereRadius == 0)
      {
        G4double halfDiagonal;
        if(userIrradiationDiagonal == 0)
        halfDiagonal = sqrt(pow((crystalx + esrThickness) * ncrystalx,2.0) + pow( (crystaly + esrThickness) * ncrystaly ,2.0)) / 2.0;
        else
        halfDiagonal = userIrradiationDiagonal/2.0;
        G4double angleLimit = atan(halfDiagonal / distance);
        //theta = (G4UniformRand() * 2.0*angleLimit) - angleLimit; //WRONG!!!! this would make cos(theta) uniform, not sin(theta)d(theta)
        //find the limit for acos
        double acosMin = cos(angleLimit);
        //so acos will have to be generated uniformely between acosMin and +1
        double randomNum =  G4UniformRand()*(1.0 - acosMin)  + (acosMin);
        theta = acos(randomNum);
        phi = G4UniformRand() * 2.0 * CLHEP::pi;
        directionVector = G4ThreeVector(sin(theta)*sin(phi),sin(theta)*cos(phi),cos(theta) );
      }
      else
      {
        // G4ThreeVector source = G4ThreeVector(EventSourcex,EventSourcey,EventSourcez);
        // G4ThreeVector target = G4ThreeVector(EventTargetx,EventTargety,EventTargetz);

        G4ThreeVector tempVector = G4ThreeVector(EventTargetx-EventSourcex,EventTargety-EventSourcey,EventTargetz-EventSourcez) ;
        float module = sqrt(pow(tempVector.getX(),2)+pow(tempVector.getY(),2)+pow(tempVector.getX(),2));
        directionVector = tempVector /module;
      }
    }

    //TODO save direction vector
    // CreateTree::Instance()->SourceMomentumX = (Float_t) directionVector.getX();
    // CreateTree::Instance()->SourceMomentumY = (Float_t) directionVector.getY();
    // CreateTree::Instance()->SourceMomentumZ = (Float_t) directionVector.getZ();

    // CreateTree::Instance()->SourceY = EventSourcey;
    // CreateTree::Instance()->SourceZ = EventSourcez;
    fParticleGun->SetParticlePosition(G4ThreeVector(EventSourcex,EventSourcey,EventSourcez));
    fParticleGun->SetParticleMomentumDirection(directionVector);
    //G4cout << "Theta = " << theta << G4endl;
    //G4cout << "Phi = " << phi << G4endl;

    // if(direction == 0) // parallel to z axis, towards matrix
    // {
    //   fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
    // }
    // else if(sphereRadius == 0) // from pointlike source, with a given aperture
    // {
    //   fParticleGun->SetParticleMomentumDirection(G4ThreeVector(sin(theta)*sin(phi),sin(theta)*cos(phi),cos(theta) ));
    // }
    // else // from spherical source, towards the matrix
    // {
    //   G4ThreeVector source = G4ThreeVector(EventSourcex,EventSourcey,EventSourcez);
    //   G4ThreeVector target = G4ThreeVector(EventTargetx,EventTargety,EventTargetz);
    //   directionVector = target - source;
    //
    // }
    fParticleGun->GeneratePrimaryVertex(anEvent);

    // DEBUG output to check source position and target position
    // std::cout << EventSourcex << " "
    //           << EventSourcey << " "
    //           << EventSourcez << G4endl;

  }
  std::ofstream myfile;
  myfile.open ("SourcePosition.txt",std::ios::app);
  myfile << EventSourcex << " " << EventSourcey << " " << EventSourcez << " " << EventTargetx << " " << EventTargety << " " << EventTargetz << std::endl;
  myfile.close();


}
