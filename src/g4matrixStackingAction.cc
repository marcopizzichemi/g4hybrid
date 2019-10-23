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
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "g4matrixStackingAction.hh"

#include "G4VProcess.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4Track.hh"
#include "G4ios.hh"
#include "CreateTree.hh"




#include "G4SDManager.hh"
#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4HCofThisEvent.hh"
#include "G4Track.hh"
#include "G4TrackStatus.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4ios.hh"



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

g4matrixStackingAction::g4matrixStackingAction()
: G4UserStackingAction(),
fScintillationCounter(0), fCerenkovCounter(0)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

g4matrixStackingAction::~g4matrixStackingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ClassificationOfNewTrack
g4matrixStackingAction::ClassifyNewTrack(const G4Track * aTrack)
{
  //if it's a primary, check if it's a chargedgeantino, and if it is, keep it only if it has been generated in a material called LYSO
  if (aTrack->GetParentID() == 0)
  {
    G4String ParticleName = aTrack->GetDynamicParticle()->GetParticleDefinition()->GetParticleName();
    if(ParticleName == "chargedgeantino")
    {
      G4String MaterialName = aTrack->GetVolume()->GetLogicalVolume()->GetMaterial()->GetName();
      if(MaterialName == "LYSO")
      {
        return fUrgent;
      }
      else
      {
        G4cout << "Primary Lu-176 killed because generated outside of LYSO material" << G4endl;
        return fKill;
      }
    }
  }
  //kill secondary neutrino
  if (aTrack->GetDefinition() == G4NeutrinoE::NeutrinoE()) return fKill;
  // else return fUrgent;

  if(aTrack->GetDefinition() == G4OpticalPhoton::OpticalPhotonDefinition())
  { 
    // particle is secondary
    if(aTrack->GetParentID()>0)
    { 
      if(aTrack->GetCreatorProcess()->GetProcessName() == "Scintillation")
        fScintillationCounter++;
      if(aTrack->GetCreatorProcess()->GetProcessName() == "Cerenkov")
        fCerenkovCounter++;
    }

    if (CreateTree::Instance()->totalEnergyDeposited < 0.5)
    {
      return fWaiting;
    }
    
  }
  return fUrgent;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void g4matrixStackingAction::NewStage()
{
  G4cout << "Number of Scintillation photons produced in this event : "
  << fScintillationCounter << G4endl;
  G4cout << "Number of Cerenkov photons produced in this event : "
  << fCerenkovCounter << G4endl;

  //if (CreateTree::Instance()->totalEnergyDeposited < 0.5)
  //{
  //  //G4UImanager * UImanager = G4UImanager::GetUIpointer();
  //  //UImanager->ApplyCommand("/event/abort");
  //  
  //  stackManager->clear();
  //  return;
  //}
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void g4matrixStackingAction::PrepareNewEvent()
{
  fScintillationCounter = 0;
  fCerenkovCounter = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......