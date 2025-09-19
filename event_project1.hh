#ifndef EVENT_HH
#define EVENT_HH

#include "G4UserEventAction.hh"
#include "G4Event.hh"
#include <map>
#include <vector>

struct HitData {
    G4String detector;
    G4ThreeVector position;
    G4ThreeVector momentum;
    G4double energy;
};

extern std::map<G4int, std::vector<HitData>> gEventHits;

class MyEventAction : public G4UserEventAction {
public:
    MyEventAction() {}
    ~MyEventAction() {}

    void EndOfEventAction(const G4Event* event) override;
};

#endif

