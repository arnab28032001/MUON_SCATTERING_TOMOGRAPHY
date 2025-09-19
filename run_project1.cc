#include "run_project1.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
std::ofstream gHitOutputFile;
MyRunAction::MyRunAction() : G4UserRunAction() {
    // Constructor: any initial setup here
}

MyRunAction::~MyRunAction() {
    // Destructor
}

void MyRunAction::BeginOfRunAction(const G4Run* run) {
    gHitOutputFile.open("hits_output.txt", std::ios::out);
    if (!gHitOutputFile.is_open()) {
        G4cerr << "Error opening hits_output.csv" << G4endl;
    }
    G4cout << "### Run " << run->GetRunID() << " started." << G4endl;
}

void MyRunAction::EndOfRunAction(const G4Run* run) {
    if (gHitOutputFile.is_open()) {
        gHitOutputFile.close();
    }
    G4cout << "### Run " << run->GetRunID() << " ended with " 
           << run->GetNumberOfEvent() << " events." << G4endl;
}

