// ============================================================
//  sim.cc  –  RPC muon telescope simulation with SD framework
// ============================================================

// ---------- Detector Construction ----------
#include "G4VUserDetectorConstruction.hh"
#include "G4SystemOfUnits.hh"
#include "G4LogicalVolume.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4ThreeVector.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4VisAttributes.hh"
// ---------- Sensitive Detector framework ----------
#include "G4VSensitiveDetector.hh"
#include "G4SDManager.hh"
#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4HCofThisEvent.hh"
#include "G4TouchableHistory.hh"
#include "G4Step.hh"
// ---------- Action Initialization ----------
#include "G4VUserActionInitialization.hh"
// ---------- Run Action ----------
#include "G4AccumulableManager.hh"
#include "G4UserRunAction.hh"
#include "G4Accumulable.hh"
#include "G4UnitsTable.hh"
#include "G4Run.hh"
// ---------- Event Action ----------
#include "G4UserEventAction.hh"
#include "G4Event.hh"
// ---------- Stepping Action ----------
#include "G4UserSteppingAction.hh"
#include "G4Track.hh"
// ---------- Primary Generator ----------
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleGun.hh"
//#include "CRYGenerator.h"
//#include "CRYParticle.h"
//#include "CRYSetup.h"
// ---------- Physics ----------
#include "FTFP_BERT.hh"
// ---------- Main ----------
#include "G4RunManagerFactory.hh"
#include "G4RunManager.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"
#include "G4UImanager.hh"
#include "Randomize.hh"
#include "globals.hh"
#include <fstream>
#include <vector>
#include "G4Threading.hh"


// ============================================================
//  Parameters
// ============================================================
class Parameters {
public:
    G4int    numberOfEvent  = 12648960;
    G4bool   GUI            = false;
    G4bool   checkOverlaps  = true;
    G4int    NumberOfThreads = 1;

    // All bare values in cm
    G4double RPCPosX[6]    = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    G4double RPCPosY[6]    = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    G4double RPCPosZ[6]    = {65.0, 55.0, 45.0, -45.0, -55.0, -65.0};

    G4double Det_sizeXY    = 120.0; // cm
    G4double Det_sizeZ     = 0.2;   // cm

    std::string filename   = "Arnab_Hits_120cm_2days.txt";
    G4String particleName  = "mu-";
    G4String Block_Material = "G4_Pb";

    G4double BlockSize1[3] = {120.0, 120.0, 60.0};
    G4double BlockSize2[3] = { 20.0,  20.0, 30.0};
    G4double BlockSize3[3] = { 20.0,  20.0, 30.0};
    G4double BlockSize4[3] = { 25.0,  25.0, 40.0};

    G4double BlockPos1[3]  = {  0.0,   0.0,  0.0}; // World centre
    G4double BlockPos2[3]  = {  0.0, -15.0,  0.0}; // Uranium
    G4double BlockPos3[3]  = {  0.0,  15.0,  0.0}; // Lead
    G4double BlockPos4[3]  = {-15.0,   0.0,  0.0}; // SS
    G4double BlockPos5[3]  = { 15.0,   0.0,  0.0}; // Air void
};


// ============================================================
//  RPCHit  –  one hit record stored per muon step in an RPC
// ============================================================
class RPCHit : public G4VHit {
public:
    RPCHit() : G4VHit(), fX(0.), fY(0.), fZ(0.), fCopyNo(-1), fKE(0.) {}
    virtual ~RPCHit() {}

    // Store hit data
    void SetPosition(G4double x, G4double y, G4double z) { fX=x; fY=y; fZ=z; }
    void SetCopyNo (G4int    n)  { fCopyNo = n; }
    void SetKE     (G4double ke) { fKE     = ke; }

    // Retrieve hit data
    G4double GetX()      const { return fX;      }
    G4double GetY()      const { return fY;      }
    G4double GetZ()      const { return fZ;      }
    G4int    GetCopyNo() const { return fCopyNo; }
    G4double GetKE()     const { return fKE;     }

private:
    G4double fX, fY, fZ;   // global position (Geant4 internal units = mm)
    G4int    fCopyNo;       // RPC layer index 0-5
    G4double fKE;           // muon kinetic energy
};

// Hits collection type
typedef G4THitsCollection<RPCHit> RPCHitsCollection;


// ============================================================
//  MySensitiveDetector  –  attached to logicRPC
// ============================================================
class MySensitiveDetector : public G4VSensitiveDetector {
public:
    MySensitiveDetector(G4String name)
        : G4VSensitiveDetector(name), fHitsCollectionID(-1), fHitsCollection(nullptr)
    {
        // Register the hits collection name with the SD framework
        collectionName.insert("RPCHitsCollection");
    }

    virtual ~MySensitiveDetector() {}

    // Called at the start of each event: create a fresh hits collection
    virtual void Initialize(G4HCofThisEvent* hce) override {
        fHitsCollection = new RPCHitsCollection(SensitiveDetectorName,
                                                collectionName[0]);
        // Get (or cache) the collection ID
        if (fHitsCollectionID < 0) {
            fHitsCollectionID = G4SDManager::GetSDMpointer()
                                    ->GetCollectionID(fHitsCollection);
        }
        hce->AddHitsCollection(fHitsCollectionID, fHitsCollection);
    }

    // Called for every step inside the sensitive volume
    virtual G4bool ProcessHits(G4Step* step, G4TouchableHistory*) override {
        // Only record muons (PDG 13 = mu−, −13 = mu+)
        G4int pid = step->GetTrack()->GetDefinition()->GetPDGEncoding();
        if (pid != 13 && pid != -13) return false;

        // Copy number tells us which of the 6 RPC layers was hit
        G4int copyNo = step->GetPreStepPoint()
                           ->GetTouchableHandle()->GetCopyNumber();
        if (copyNo < 0 || copyNo > 5) return false;

        // Fill a hit object and add it to the collection
        RPCHit* hit = new RPCHit();
        G4ThreeVector pos = step->GetPreStepPoint()->GetPosition();
        hit->SetPosition(pos.x(), pos.y(), pos.z());
        hit->SetCopyNo(copyNo);
        hit->SetKE(step->GetTrack()->GetDynamicParticle()->GetKineticEnergy());
        fHitsCollection->insert(hit);

        return true;
    }

    // (Optional) called at end of event – not needed here since EventAction reads the HC
    virtual void EndOfEvent(G4HCofThisEvent*) override {}

    G4int GetHitsCollectionID() const { return fHitsCollectionID; }

private:
    G4int              fHitsCollectionID;
    RPCHitsCollection* fHitsCollection;
};


// ============================================================
//  DetectorConstruction
// ============================================================
class DetectorConstruction : public G4VUserDetectorConstruction {
private:
    Parameters parameters;
    // No longer need fScoringVolume – the SD handles sensitivity
public:
    DetectorConstruction() : G4VUserDetectorConstruction() {}
    ~DetectorConstruction() {}

    // ---- geometry ----
    G4VPhysicalVolume* Construct() override {
        G4NistManager* nist   = G4NistManager::Instance();
        G4Material* world_mat = nist->FindOrBuildMaterial("G4_AIR");
        G4Material* det_mat   = nist->FindOrBuildMaterial("G4_Ar");
        G4Material* Lead      = nist->FindOrBuildMaterial("G4_Pb");
        G4Material* Uranium   = nist->FindOrBuildMaterial("G4_U");
        G4Material* Iron      = nist->FindOrBuildMaterial("G4_Fe");
        G4Material* Concrete  = nist->FindOrBuildMaterial("G4_CONCRETE");
        G4Material* SS        = nist->FindOrBuildMaterial("G4_STAINLESS-STEEL");

        G4double Det_sizeXY = parameters.Det_sizeXY * cm;
        G4double Det_sizeZ  = parameters.Det_sizeZ  * cm;

        // Concrete outer cylinder
        G4double shape1_rmina=0.*cm, shape1_rmaxa=25.*cm;
        G4double shape1_rminb=0.*cm, shape1_rmaxb=25.*cm;
        G4double shape1_hz=25.*cm;
        G4double shape1_phimin=0.*deg, shape1_phimax=360.*deg;
        // Uranium cylinder
        G4double shapeU_rmina=0.*cm, shapeU_rmaxa=5.*cm;
        G4double shapeU_rminb=0.*cm, shapeU_rmaxb=5.*cm;
        G4double shapeU_hz=15.*cm;
        G4double shapeU_phimin=0.*deg, shapeU_phimax=360.*deg;
        // Lead cylinder
        G4double shapePb_rmina=0.*cm, shapePb_rmaxa=5.*cm;
        G4double shapePb_rminb=0.*cm, shapePb_rmaxb=5.*cm;
        G4double shapePb_hz=15.*cm;
        G4double shapePb_phimin=0.*deg, shapePb_phimax=360.*deg;
        // SS cylinder
        G4double shapeSS_rmina=0.*cm, shapeSS_rmaxa=5.*cm;
        G4double shapeSS_rminb=0.*cm, shapeSS_rmaxb=5.*cm;
        G4double shapeSS_hz=15.*cm;
        G4double shapeSS_phimin=0.*deg, shapeSS_phimax=360.*deg;
        // Air void cylinder
        G4double shapeAir_rmina=0.*cm, shapeAir_rmaxa=5.*cm;
        G4double shapeAir_rminb=0.*cm, shapeAir_rmaxb=5.*cm;
        G4double shapeAir_hz=15.*cm;
        G4double shapeAir_phimin=0.*deg, shapeAir_phimax=360.*deg;

        // World
        G4double world_sizeXY = Det_sizeXY + 12.0*cm;
        G4double world_sizeZ  = ((parameters.RPCPosZ[0] - parameters.RPCPosZ[5]) + 10.0)*cm;
        G4bool   checkOverlaps = parameters.checkOverlaps;

        G4Box* solidWorld = new G4Box("World",
                                      0.5*world_sizeXY, 0.5*world_sizeXY, 0.5*world_sizeZ);
        G4LogicalVolume* logicWorld = new G4LogicalVolume(solidWorld, world_mat, "World");
        G4VPhysicalVolume* physWorld = new G4PVPlacement(
                0, G4ThreeVector(), logicWorld, "World", nullptr, false, 0, checkOverlaps);

        // RPC layers (6 copies of the same logical volume)
        G4Box* solidRPC = new G4Box("RPC", Det_sizeXY/2.0, Det_sizeXY/2.0, Det_sizeZ/2.0);
        G4LogicalVolume* logicRPC = new G4LogicalVolume(solidRPC, det_mat, "RPC");

        for (G4int i = 0; i < 6; i++) {
            G4ThreeVector pos(parameters.RPCPosX[i]*cm,
                              parameters.RPCPosY[i]*cm,
                              parameters.RPCPosZ[i]*cm);
            new G4PVPlacement(0, pos, logicRPC, "RPC", logicWorld, false, i, checkOverlaps);
        }

        // Block geometry
        G4ThreeVector bPos1(parameters.BlockPos1[0]*cm, parameters.BlockPos1[1]*cm, parameters.BlockPos1[2]*cm);
        G4ThreeVector bPos2(parameters.BlockPos2[0]*cm, parameters.BlockPos2[1]*cm, parameters.BlockPos2[2]*cm);
        G4ThreeVector bPos3(parameters.BlockPos3[0]*cm, parameters.BlockPos3[1]*cm, parameters.BlockPos3[2]*cm);
        G4ThreeVector bPos4(parameters.BlockPos4[0]*cm, parameters.BlockPos4[1]*cm, parameters.BlockPos4[2]*cm);
        G4ThreeVector bPos5(parameters.BlockPos5[0]*cm, parameters.BlockPos5[1]*cm, parameters.BlockPos5[2]*cm);

        G4Cons* solidCon  = new G4Cons("Cyl_Con",  shape1_rmina,   shape1_rmaxa,   shape1_rminb,   shape1_rmaxb,   shape1_hz,    shape1_phimin,   shape1_phimax);
        G4Cons* solidU    = new G4Cons("Cyl_U",    shapeU_rmina,   shapeU_rmaxa,   shapeU_rminb,   shapeU_rmaxb,   shapeU_hz,    shapeU_phimin,   shapeU_phimax);
        G4Cons* solidPb   = new G4Cons("Cyl_Pb",   shapePb_rmina,  shapePb_rmaxa,  shapePb_rminb,  shapePb_rmaxb,  shapePb_hz,   shapePb_phimin,  shapePb_phimax);
        G4Cons* solidSS   = new G4Cons("Cyl_SS",   shapeSS_rmina,  shapeSS_rmaxa,  shapeSS_rminb,  shapeSS_rmaxb,  shapeSS_hz,   shapeSS_phimin,  shapeSS_phimax);
        G4Cons* solidAir  = new G4Cons("Cyl_Air",  shapeAir_rmina, shapeAir_rmaxa, shapeAir_rminb, shapeAir_rmaxb, shapeAir_hz,  shapeAir_phimin, shapeAir_phimax);

        G4LogicalVolume* logicCon  = new G4LogicalVolume(solidCon,  Concrete,  "Cyl_Con");
        G4LogicalVolume* logicU    = new G4LogicalVolume(solidU,    Uranium,   "Cyl_U");
        G4LogicalVolume* logicPb   = new G4LogicalVolume(solidPb,   Lead,      "Cyl_Pb");
        G4LogicalVolume* logicSS   = new G4LogicalVolume(solidSS,   SS,        "Cyl_SS");
        G4LogicalVolume* logicAir  = new G4LogicalVolume(solidAir,  world_mat, "Cyl_Air");

        new G4PVPlacement(0, bPos1, logicCon,  "Cyl_Con", logicWorld, false, 0, checkOverlaps);
        new G4PVPlacement(0, bPos2, logicU,    "Cyl_U",   logicCon,   false, 0, checkOverlaps);
        new G4PVPlacement(0, bPos3, logicPb,   "Cyl_Pb",  logicCon,   false, 0, checkOverlaps);
        new G4PVPlacement(0, bPos4, logicSS,   "Cyl_SS",  logicCon,   false, 0, checkOverlaps);
        new G4PVPlacement(0, bPos5, logicAir,  "Cyl_Air", logicCon,   false, 0, checkOverlaps);

        // Visualisation attributes
        logicWorld->SetVisAttributes(new G4VisAttributes(G4Colour(0.0, 0.0, 1.0, 0.1)));

        auto* vCon = new G4VisAttributes(G4Colour(0.1, 0.1, 0.1, 0.4));
        vCon->SetForceSolid(true);  logicCon->SetVisAttributes(vCon);

        auto* vU = new G4VisAttributes(G4Colour(0.0, 1.0, 0.0, 1.0));
        vU->SetForceSolid(true);    logicU->SetVisAttributes(vU);

        auto* vPb = new G4VisAttributes(G4Colour(1.0, 0.0, 0.0, 1.0));
        vPb->SetForceSolid(true);   logicPb->SetVisAttributes(vPb);

        auto* vSS = new G4VisAttributes(G4Colour(0.0, 0.0, 1.0, 1.0));
        vSS->SetForceSolid(true);   logicSS->SetVisAttributes(vSS);

        auto* vAir = new G4VisAttributes(G4Colour(0.0, 0.0, 0.0, 1.0));
        vAir->SetForceSolid(true);  logicAir->SetVisAttributes(vAir);

        auto* vRPC = new G4VisAttributes(G4Colour(0.0, 0.0, 0.0, 0.1));
        vRPC->SetForceSolid(true);  logicRPC->SetVisAttributes(vRPC);

        return physWorld;
    }

    // ---- sensitive detector assignment ----
    // Geant4 calls ConstructSDandField() on every worker thread automatically.
    // This is the CORRECT place to create and attach SDs (not inside Construct()).
    void ConstructSDandField() override {
        G4SDManager* sdManager = G4SDManager::GetSDMpointer();

        MySensitiveDetector* detectorSD = new MySensitiveDetector("MyDetectorSD");
        sdManager->AddNewDetector(detectorSD);

        // Attach the SD to the RPC logical volume by name.
        // SetSensitiveDetector() looks up the logical volume in the store,
        // so it works correctly even in multi-threaded mode.
        SetSensitiveDetector("RPC", detectorSD);
    }
};


// ============================================================
//  PrimaryGeneratorAction
// ============================================================
class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
private:
    Parameters     parameters;
    G4ParticleGun* fParticleGun;
    std::ifstream  infile;
public:
    PrimaryGeneratorAction() : G4VUserPrimaryGeneratorAction(), fParticleGun(nullptr) {
        fParticleGun = new G4ParticleGun(1);
        infile.open("/home/oem/MSC_DISSERTATION_ARNAB/15_08_25/ARNAB _GEOMETRY/CRY_Data_File.txt");
        if (!infile.is_open())
            G4cerr << "ERROR: Cannot open CRY data file!" << G4endl;
    }
    ~PrimaryGeneratorAction() { delete fParticleGun; }

    void GeneratePrimaries(G4Event* anEvent) override {
        G4ParticleTable*     pt  = G4ParticleTable::GetParticleTable();
        G4ParticleDefinition* pd = pt->FindParticle(parameters.particleName);
        fParticleGun->SetParticleDefinition(pd);

        // Spread gun uniformly over (Det_sizeXY − 10 cm) in X and Y.
        // Det_sizeXY is stored in cm (bare number), so *cm is required before
        // mixing with other Geant4-unit quantities.
        G4double halfSpread = (parameters.Det_sizeXY*cm - 10.0*cm) * 0.5;
        G4double x0 = (2.0*G4UniformRand() - 1.0) * halfSpread; // already in mm
        G4double y0 = (2.0*G4UniformRand() - 1.0) * halfSpread;
        G4double z0 = (parameters.RPCPosZ[0] + parameters.Det_sizeZ + 2.0) * cm;

        // Read one CRY event
        G4int PID; G4double En, px, py, pz, xx,yy;
        infile >>PID>>En>>xx>>yy>>px>>py>>pz;
        fParticleGun->SetParticleEnergy(En*MeV);
        fParticleGun->SetParticleMomentumDirection(G4ThreeVector(px, py, pz));
        fParticleGun->SetParticlePosition(G4ThreeVector(x0, y0, z0));
        fParticleGun->GeneratePrimaryVertex(anEvent);
    }
};


// ============================================================
//  RunAction
// ============================================================
class RunAction : public G4UserRunAction {
public:
    RunAction()  : G4UserRunAction() {}
    ~RunAction() {}
    void BeginOfRunAction(const G4Run*) override {}
    void EndOfRunAction(const G4Run* run) override {
        if (run->GetNumberOfEvent() == 0) return;
        if (IsMaster())
            G4cout << G4endl << "------> End of Global Run." << G4endl;
        else
            G4cout << G4endl << "-------- End of Local Run." << G4endl;
    }
};


// ============================================================
//  EventAction  –  reads the SD hits collection and writes file
// ============================================================
class EventAction : public G4UserEventAction {
public:
    Parameters   parameters;
    std::fstream outfile;

    EventAction(RunAction* runAction) : G4UserEventAction() { (void)runAction; }
    ~EventAction() {}

    void BeginOfEventAction(const G4Event*) override {}  // nothing needed: SD resets HC automatically

    void EndOfEventAction(const G4Event* event) override {
        G4int eventID = event->GetEventID();
        if (eventID % 10000 == 0)
            G4cout << "Total Events processed: " << eventID << G4endl;

        // ---- Retrieve the hits collection from the event ----
        G4HCofThisEvent* hce = event->GetHCofThisEvent();
        if (!hce) return;

        // Look up collection ID by name (cached after first call)
        static G4int hcID = -1;
        if (hcID < 0)
            hcID = G4SDManager::GetSDMpointer()->GetCollectionID("MyDetectorSD/RPCHitsCollection");

        RPCHitsCollection* hc =
            static_cast<RPCHitsCollection*>(hce->GetHC(hcID));
        if (!hc) return;

        // ---- Consolidate: keep only the FIRST hit per RPC layer ----
        // (A single muon track may produce several steps inside one thin RPC;
        //  we want the entry point, i.e. the first step recorded per layer.)
        G4double hitX[6]    = {0,0,0,0,0,0};
        G4double hitY[6]    = {0,0,0,0,0,0};
        G4bool   hitFlag[6] = {false,false,false,false,false,false};

        G4int nHits = (G4int)hc->GetSize();
        for (G4int i = 0; i < nHits; i++) {
            RPCHit* hit = (*hc)[i];
            G4int   lay = hit->GetCopyNo();
            if (lay < 0 || lay > 5) continue;
            if (!hitFlag[lay]) {          // take only the first step per layer
                hitX[lay]    = hit->GetX();
                hitY[lay]    = hit->GetY();
                hitFlag[lay] = true;
            }
        }

        // ---- Only write events where all 6 layers fired ----
        for (G4int i = 0; i < 6; i++)
            if (!hitFlag[i]) return;

        // ---- Write hit positions to file ----
        // X,Y: divide by mm → numerical value in millimetres
        // Z:   RPCPosZ[i] is in cm; ×10 converts to mm (same unit as X,Y)
        outfile.open(parameters.filename, std::ios::out | std::ios::app);
        for (G4int i = 0; i < 6; i++) {
            outfile << hitX[i]/mm << " "
                    << hitY[i]/mm << " "
                    << parameters.RPCPosZ[i]*10.0;
            if (i < 5) outfile << " ";
        }
        outfile << G4endl;
        outfile.close();
    }
};


// ============================================================
//  SteppingAction  –  kept for structural compatibility;
//  all hit recording is now handled by MySensitiveDetector.
// ============================================================
class SteppingAction : public G4UserSteppingAction {
public:
    SteppingAction(EventAction*) : G4UserSteppingAction() {}
    ~SteppingAction() {}
    void UserSteppingAction(const G4Step*) override { /* SD does the work */ }
};


// ============================================================
//  ActionInitialization
// ============================================================
class ActionInitialization : public G4VUserActionInitialization {
public:
    ActionInitialization()  : G4VUserActionInitialization() {}
    ~ActionInitialization() {}

    void BuildForMaster() const override {
        SetUserAction(new RunAction());
    }

    void Build() const override {
        RunAction*   runAction   = new RunAction();
        EventAction* eventAction = new EventAction(runAction);
        SetUserAction(runAction);
        SetUserAction(eventAction);
        SetUserAction(new SteppingAction(eventAction));
        SetUserAction(new PrimaryGeneratorAction());
    }
};


// ============================================================
//  main
// ============================================================
int main(int argc, char** argv) {
    Parameters parameters;
    G4bool GUI = parameters.GUI;

    G4Random::setTheEngine(new CLHEP::RanecuEngine);
    auto* runManager = G4RunManagerFactory::CreateRunManager(G4RunManagerType::Default);
    runManager->SetNumberOfThreads(parameters.NumberOfThreads);

    runManager->SetUserInitialization(new DetectorConstruction());
    runManager->SetUserInitialization(new FTFP_BERT(0));
    runManager->SetUserInitialization(new ActionInitialization());

    if (!GUI) {
        runManager->Initialize();
        runManager->BeamOn(parameters.numberOfEvent);
        delete runManager;
    } else {
        G4UIExecutive* ui = nullptr;
        if (argc == 1) ui = new G4UIExecutive(argc, argv);
        G4VisManager* visManager = new G4VisExecutive;
        visManager->Initialize();
        G4UImanager* UImanager = G4UImanager::GetUIpointer();
        if (!ui) {
            UImanager->ApplyCommand("/control/execute " + G4String(argv[1]));
        } else {
            UImanager->ApplyCommand("/control/execute init_vis.mac");
            ui->SessionStart();
            delete ui;
        }
        delete visManager;
        delete runManager;
    }
}
