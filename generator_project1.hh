#ifndef GENERATOR_HH
#define GENERATOR_HH

#include "G4VUserPrimaryGeneratorAction.hh"
#include "CRYGenerator.h"
#include "CRYSetup.h"
#include "CRYParticle.h"
#include "G4ParticleGun.hh"
#include <vector>


class MyPrimaryGenerator : public G4VUserPrimaryGeneratorAction {
public:
    MyPrimaryGenerator();
    ~MyPrimaryGenerator();
    virtual void GeneratePrimaries(G4Event* anEvent);

private:
    CRYGenerator* cryGen;
    G4ParticleTable* particleTable;
    G4ParticleGun* particleGun;
};


#endif
