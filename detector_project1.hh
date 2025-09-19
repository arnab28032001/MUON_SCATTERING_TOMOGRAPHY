#ifndef DETECTOR_HH
#define DETECTOR_HH

#include "G4VSensitiveDetector.hh"
#include "G4Step.hh"
#include "G4HCofThisEvent.hh"
#include "G4TouchableHistory.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4StepPoint.hh"
#include "G4Track.hh"
#include <fstream>

class MySensitiveDetector : public G4VSensitiveDetector {
public:
    MySensitiveDetector(const G4String& name);
    virtual ~MySensitiveDetector();

    virtual G4bool ProcessHits(G4Step* step, G4TouchableHistory* hist);

private:
    static std::ofstream fOutputFile;
};

#endif

