#include "generator_project1.hh"
#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4Event.hh"
#include "CRYSetup.h"
#include <fstream>

MyPrimaryGenerator::MyPrimaryGenerator() {
    std::ifstream inputFile("cry_setup.txt");  // your CRY setup file
    std::string setupStr((std::istreambuf_iterator<char>(inputFile)),
                         std::istreambuf_iterator<char>());

    CRYSetup* crySetup = new CRYSetup(setupStr, "/home/oem/MSC_DISSERTATION/15_08_25/Geometry_Build_4/cry_v1.7/data");
    cryGen = new CRYGenerator(crySetup);
    

    particleGun = new G4ParticleGun(1);
    particleTable = G4ParticleTable::GetParticleTable();
}

MyPrimaryGenerator::~MyPrimaryGenerator() {
   delete particleGun;
    delete cryGen;
}

void MyPrimaryGenerator::GeneratePrimaries(G4Event* anEvent) {
    std::vector<CRYParticle*>* particles = new std::vector<CRYParticle*>;
    cryGen->genEvent(particles);

    for (auto p : *particles) {
        G4int name = p->PDGid();
        G4ThreeVector pos(p->x()*cm, p->y()*cm, p->z()*cm);
        G4ThreeVector mom(p->u(), p->v(), p->w());
        mom = mom.unit();

        G4ParticleDefinition* particleDef = particleTable->FindParticle(name);
        if (!particleDef) continue;
        G4double xRange = 80.0 * cm;
	G4double yRange = 80.0 * cm;
	G4double zPos   = 66.5 * cm;
	
	// Generate random X and Y within that range
	G4double x = (G4UniformRand() - 0.5) * xRange;  // random between -75 cm and 75 cm
	G4double y = (G4UniformRand() - 0.5) * yRange;  // random between -75 cm and 75 cm
	
	// Set the randomized position
	G4ThreeVector position(x, y, zPos);
	particleGun->SetParticlePosition(position);
	
	// Momentum direction always downward (-z)
	G4ThreeVector momDir(0., 0., -1.);  // pointing down
	particleGun->SetParticleMomentumDirection(momDir);
	// Particle definition and energy as before
	particleGun->SetParticleDefinition(particleDef);
	particleGun->SetParticleEnergy(p->ke() * MeV);
	
	// Generate the primary vertex
	particleGun->GeneratePrimaryVertex(anEvent);

	
       
    }

    for (auto p : *particles) delete p;
    delete particles;
}


