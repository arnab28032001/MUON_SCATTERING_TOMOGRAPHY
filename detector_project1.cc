#include "detector_project1.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4TouchableHistory.hh"
#include "G4SystemOfUnits.hh"
#include "G4RunManager.hh"
#include <iostream>
#include <fstream>
#include "event_project1.hh"
//std::ofstream outputFile("hits_output1.txt", std::ios::out);

MySensitiveDetector::MySensitiveDetector(const G4String& name)
    : G4VSensitiveDetector(name)
{
    /*if (!outputFile.is_open()) {
        G4cerr << "Error opening output file!" << G4endl;
    } */
}

MySensitiveDetector::~MySensitiveDetector()
{
    /*if (outputFile.is_open()) {
        outputFile.close();
    }*/
}

G4bool MySensitiveDetector::ProcessHits(G4Step* step, G4TouchableHistory*)
{
	auto pre = step->GetPreStepPoint();
        auto track = step->GetTrack();
	
	G4int eventID = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
	HitData hit;
	hit.position = pre->GetPosition() / cm;
	//gEventHits[eventID].clear();
	
	gEventHits[eventID].push_back(hit);
    
      
    /*G4double edep = step->GetTotalEnergyDeposit();
    if (edep == 0.) return false;

    G4StepPoint* prePoint = step->GetPreStepPoint();
    G4ThreeVector pos = prePoint->GetPosition();
    G4String particleName = step->GetTrack()->GetDefinition()->GetParticleName();
    G4String volName = prePoint->GetTouchableHandle()->GetVolume()->GetName();

    if (outputFile.is_open()) {
        outputFile << volName << "," 
                   << particleName << ","
                   << edep/MeV << ","
                   << pos.x()/cm << ","
                   << pos.y()/cm << ","
                   << pos.z()/cm << std::endl;
    }*/

    return true;
}

