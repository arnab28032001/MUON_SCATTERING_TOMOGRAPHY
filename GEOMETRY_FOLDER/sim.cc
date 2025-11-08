//Detector Construction class.
#include "G4VUserDetectorConstruction.hh"
#include "G4SystemOfUnits.hh"
#include "G4LogicalVolume.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4ThreeVector.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4VisAttributes.hh"
// Action initialization class.
#include "G4VUserActionInitialization.hh"
// Run Action class.
#include "G4AccumulableManager.hh"
#include "G4UserRunAction.hh"
#include "G4Accumulable.hh"
#include "G4UnitsTable.hh"
#include "G4Run.hh"
// Event action class.
#include "G4UserEventAction.hh"
#include "G4Event.hh"
// Stepping action class.
#include "G4UserSteppingAction.hh"
#include "G4Track.hh"
#include "G4Step.hh"
// Primary generator action class.
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleGun.hh"
//#include "CRYGenerator.h"
//#include "CRYParticle.h"
//#include "CRYSetup.h"
// Physics list.
#include "FTFP_BERT.hh"
// Main.
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


class Parameters {
public:
    G4int numberOfEvent     = 7*12648960; // 126489600 for 10days data (or 10*24 hours), flux region 120*120 cm*cm;
    G4bool GUI              = false;
    G4bool checkOverlaps    = true;
    G4int NumberOfThreads   = 1;
    // All units in cm
        //Detector Position
    G4double RPCPosX[6]     = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    G4double RPCPosY[6]     = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    G4double RPCPosZ[6]     =  {65.0, 55.0, 45.0, -45.0, -55.0, -65.0};
    // Detector Size
    G4double  Det_sizeXY    = 120.0;
    G4double  Det_sizeZ     = 0.2;
    std::string filename    = "Arnab_Hits_120cm_7days.txt";
    G4String particleName   = "mu-";
    G4String Block_Material = "G4_Pb";
  
    // Block Size and pos
    G4double BlockSize1[3]  = { 120.0,  120.0, 60.0};//World size
    G4double BlockSize2[3]  = { 20.0,  20.0, 30.0};//Uranium Size
    G4double BlockSize3[3]  = { 20.0,  20.0, 30.0};// Lead Size
    G4double BlockSize4[3]  = { 25.0,  25.0, 40.0};//Air Size
  
  
    G4double BlockPos1[3]   = {0.0, 0.0, 0.0};//World
    G4double BlockPos2[3]   = {0.0, -15.0, 0.0}; //uranium
  
    G4double BlockPos3[3]   = {0.0, 15.0, 0.0}; // Lead
    G4double BlockPos4[3]   = {-15.0, 0.0, 0.0}; //SS
    G4double BlockPos5[3]   = {15.0, 0.0, 0.0}; //Air void
  
};
//Construction.cc
class DetectorConstruction : public G4VUserDetectorConstruction {
private:
    Parameters parameters;
    G4LogicalVolume*  fScoringVolume;
public:
    DetectorConstruction(): G4VUserDetectorConstruction(), fScoringVolume(0) { }
    ~DetectorConstruction(){ }
    G4VPhysicalVolume* Construct(){
        G4NistManager* nist     = G4NistManager::Instance();
        G4Material* world_mat   = nist->FindOrBuildMaterial("G4_AIR");
        G4Material* det_mat     = nist->FindOrBuildMaterial("G4_Ar");
        G4Material* Lead        = nist->FindOrBuildMaterial("G4_Pb");
        G4Material* Uranium     = nist->FindOrBuildMaterial("G4_U");
        G4Material* Iron        = nist->FindOrBuildMaterial("G4_Fe");
	G4Material* Concrete    = nist->FindOrBuildMaterial("G4_CONCRETE");
	G4Material* SS = nist->FindOrBuildMaterial("G4_STAINLESS-STEEL");
	
// concrete cylinder	

        G4double  Det_sizeXY    = parameters.Det_sizeXY *cm;
        G4double  Det_sizeZ     = parameters.Det_sizeZ *cm;
	G4double shape1_rmina =  0.*cm, shape1_rmaxa = 25.*cm;
	G4double shape1_rminb =  0.*cm, shape1_rmaxb = 25.*cm;
	G4double shape1_hz = 25.*cm;//length
	G4double shape1_phimin = 0.*deg, shape1_phimax = 360.*deg;
// Uranium cylinder
	G4double shapeU_rmina =  0.*cm, shapeU_rmaxa = 5.*cm;
	G4double shapeU_rminb =  0.*cm, shapeU_rmaxb = 5.*cm;
	G4double shapeU_hz = 15.*cm;//length
	G4double shapeU_phimin = 0.*deg, shapeU_phimax = 360.*deg;
// Lead cylinder
           G4double shapePb_rmina =  0.*cm, shapePb_rmaxa = 5.*cm;
	G4double shapePb_rminb =  0.*cm, shapePb_rmaxb = 5.*cm;
	G4double shapePb_hz = 15.*cm;//length
	G4double shapePb_phimin = 0.*deg, shapePb_phimax = 360.*deg;
	
// Stainless Steel cylinder
           
        G4double shapeSS_rmina =  0.*cm, shapeSS_rmaxa = 5.*cm;
	G4double shapeSS_rminb =  0.*cm, shapeSS_rmaxb = 5.*cm;
	G4double shapeSS_hz = 15.*cm;//length
	G4double shapeSS_phimin = 0.*deg, shapeSS_phimax = 360.*deg;
	
// Air void cylinder	
	
	G4double shapeAir_rmina =  0.*cm, shapeAir_rmaxa = 5.*cm;
	G4double shapeAir_rminb =  0.*cm, shapeAir_rmaxb = 5.*cm;
	G4double shapeAir_hz = 15.*cm;
	G4double shapeAir_phimin = 0.*deg, shapeAir_phimax = 360.*deg;
// World
        G4double world_sizeXY   = Det_sizeXY + 12.0 *cm; // must be greater than Det_sizeXY
        G4double world_sizeZ    = ((parameters.RPCPosZ[0] - parameters.RPCPosZ[5])+10.0)*cm;
        G4bool checkOverlaps    = parameters.checkOverlaps;
        G4Box* solidWorld = new G4Box("World", 0.5*world_sizeXY, 0.5*world_sizeXY, 0.5*world_sizeZ);
        
        G4LogicalVolume* logicWorld = new G4LogicalVolume(solidWorld, world_mat, "World");
        
        G4VPhysicalVolume* physWorld = new G4PVPlacement(0, G4ThreeVector(), logicWorld, "World", 0, false, 0, checkOverlaps);
// RPC
        G4ThreeVector Pos_1, Pos_2, Pos_3, Pos_4, Pos_5, Pos_6;
        Pos_1 = G4ThreeVector(parameters.RPCPosX[0]*cm, parameters.RPCPosY[0]*cm,  parameters.RPCPosZ[0]*cm);
        Pos_2 = G4ThreeVector(parameters.RPCPosX[1]*cm, parameters.RPCPosY[1]*cm,  parameters.RPCPosZ[1]*cm);
        Pos_3 = G4ThreeVector(parameters.RPCPosX[2]*cm, parameters.RPCPosY[2]*cm,  parameters.RPCPosZ[2]*cm);
        Pos_4 = G4ThreeVector(parameters.RPCPosX[3]*cm, parameters.RPCPosY[3]*cm,  parameters.RPCPosZ[3]*cm);
        Pos_5 = G4ThreeVector(parameters.RPCPosX[4]*cm, parameters.RPCPosY[4]*cm,  parameters.RPCPosZ[4]*cm);
        Pos_6 = G4ThreeVector(parameters.RPCPosX[5]*cm, parameters.RPCPosY[5]*cm,  parameters.RPCPosZ[5]*cm);
        G4Box* solidRPC = new G4Box("RPC", Det_sizeXY/2.0, Det_sizeXY/2.0, Det_sizeZ/2.0);
        G4LogicalVolume* logicRPC = new G4LogicalVolume(solidRPC, det_mat, "RPC");
        new G4PVPlacement(0, Pos_1, logicRPC, "RPC", logicWorld,  false, 0, checkOverlaps);     // RPC 1
        new G4PVPlacement(0, Pos_2, logicRPC, "RPC", logicWorld,  false, 1, checkOverlaps);     // RPC 2
        new G4PVPlacement(0, Pos_3, logicRPC, "RPC", logicWorld,  false, 2, checkOverlaps);     // RPC 3
        new G4PVPlacement(0, Pos_4, logicRPC, "RPC", logicWorld,  false, 3, checkOverlaps);     // RPC 4
        new G4PVPlacement(0, Pos_5, logicRPC, "RPC", logicWorld,  false, 4, checkOverlaps);     // RPC 5
        new G4PVPlacement(0, Pos_6, logicRPC, "RPC", logicWorld,  false, 5, checkOverlaps);     // RPC 6

// Block
        G4ThreeVector blockPos1 = G4ThreeVector(parameters.BlockPos1[0]*cm, parameters.BlockPos1[1]*cm, parameters.BlockPos1[2]*cm);
	G4ThreeVector blockPos2 = G4ThreeVector(parameters.BlockPos2[0]*cm, parameters.BlockPos2[1]*cm, parameters.BlockPos2[2]*cm);
        G4ThreeVector blockPos3 = G4ThreeVector(parameters.BlockPos3[0]*cm, parameters.BlockPos3[1]*cm, parameters.BlockPos3[2]*cm);
	G4ThreeVector blockPos4 = G4ThreeVector(parameters.BlockPos4[0]*cm, parameters.BlockPos4[1]*cm, parameters.BlockPos4[2]*cm);
	G4ThreeVector blockPos5 = G4ThreeVector(parameters.BlockPos5[0]*cm, parameters.BlockPos5[1]*cm, parameters.BlockPos5[2]*cm);

	G4Cons* solidShape1 =  new G4Cons("Cyl_Con", shape1_rmina, shape1_rmaxa, shape1_rminb, shape1_rmaxb, shape1_hz, shape1_phimin, shape1_phimax);
        G4Cons* solidBlock1 = new G4Cons("Cyl_U", shapeU_rmina, shapeU_rmaxa, shapeU_rminb, shapeU_rmaxb, shapeU_hz, shapeU_phimin, shapeU_phimax);
	G4Cons* solidBlock2 = new G4Cons("Cyl_Pb", shapePb_rmina, shapePb_rmaxa, shapePb_rminb, shapePb_rmaxb, shapePb_hz, shapePb_phimin, shapePb_phimax);
	G4Cons* solidBlock3 = new G4Cons("Cyl_SS", shapeSS_rmina, shapeSS_rmaxa, shapeSS_rminb, shapeSS_rmaxb, shapeSS_hz, shapeSS_phimin, shapeSS_phimax);
	G4Cons* solidBlock4 = new G4Cons("Cyl_Air",  shapeAir_rmina, shapeAir_rmaxa, shapeAir_rminb, shapeAir_rmaxb, shapeAir_hz, shapeAir_phimin, shapeAir_phimax);
	
	G4LogicalVolume* logicCyl1 = new G4LogicalVolume(solidShape1, Concrete, "Cyl_Con");
        G4LogicalVolume* logicBlock1 = new G4LogicalVolume(solidBlock1, Uranium, "Cyl_U");
        G4LogicalVolume* logicBlock2 = new G4LogicalVolume(solidBlock2, Lead, "Cyl_Pb");
        G4LogicalVolume* logicBlock3 = new G4LogicalVolume(solidBlock3, SS, "Cyl_SS");
        G4LogicalVolume* logicBlock4 = new G4LogicalVolume(solidBlock4, world_mat, "Cyl_Air");
	
	new G4PVPlacement(0, blockPos1, logicCyl1, "Cyl_Con", logicWorld, false, 0, checkOverlaps);
        new G4PVPlacement(0, blockPos2, logicBlock1, "Cyl_U", logicCyl1, false, 0, checkOverlaps);
        new G4PVPlacement(0, blockPos3, logicBlock2, "Cyl_Pb", logicCyl1, false, 0, checkOverlaps);
        new G4PVPlacement(0, blockPos4, logicBlock3, "Cyl_SS", logicCyl1, false, 0, checkOverlaps);
        new G4PVPlacement(0, blockPos5, logicBlock4, "Cyl_Air", logicCyl1, false, 0, checkOverlaps);
        
        auto worldVisAttributes = new G4VisAttributes(G4Colour(0.0, 0.0, 1.0, 0.2)); // Blue with transparency
        worldVisAttributes->SetVisibility(true);
        logicWorld->SetVisAttributes(worldVisAttributes);

        G4VisAttributes* visAttConcrete = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5, 0.5)); // Very transparent light gray
        visAttConcrete->SetVisibility(true);
        visAttConcrete->SetForceSolid(true);
        logicCyl1->SetVisAttributes(visAttConcrete);

        G4VisAttributes* visAtt1 = new G4VisAttributes(G4Colour(0.0, 1.0, 0.0, 1.0)); // Green (Opaque)
        G4VisAttributes* visAtt2 = new G4VisAttributes(G4Colour(1.0, 0.0, 0.0, 1.0)); // Red (Opaque)
        G4VisAttributes* visAtt3 = new G4VisAttributes(G4Colour(0.0, 0.0, 1.0, 1.0)); // Blue (Opaque)
        G4VisAttributes* visAtt4 = new G4VisAttributes(G4Colour(1.0, 1.0, 0.0, 1.0)); // Yellow (Opaque)
        G4VisAttributes* visAtt5 = new G4VisAttributes(G4Colour(0.0, 0.0, 0.0, .1)); // Black (10% Opaque)

        visAtt1->SetVisibility(true);
        visAtt1->SetForceSolid(true);
        logicBlock1->SetVisAttributes(visAtt1);

        visAtt2->SetVisibility(true);
        visAtt2->SetForceSolid(true);
        logicBlock2->SetVisAttributes(visAtt2);

        visAtt3->SetVisibility(true);
        visAtt3->SetForceSolid(true);
        logicBlock3->SetVisAttributes(visAtt3);

        visAtt4->SetVisibility(true);
        visAtt4->SetForceSolid(true);
        logicBlock4->SetVisAttributes(visAtt4);
        
        visAtt5->SetVisibility(true);
        visAtt5->SetForceSolid(true);
        logicRPC->SetVisAttributes(visAtt5);
	
        // Set ScoringVolumex
        fScoringVolume = logicRPC;
        //always return the physical World
        return physWorld;
    }
    G4LogicalVolume* GetScoringVolume() const { return fScoringVolume; }
};
//Generator.cc
class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
private:
    Parameters parameters;
    G4ParticleGun*  fParticleGun;
    std::ifstream infile;
public:
    PrimaryGeneratorAction() : G4VUserPrimaryGeneratorAction(), fParticleGun(0) {
        fParticleGun  = new G4ParticleGun(1);
        infile.open ("/home/oem/MSC_DISSERTATION/15_08_25/ARNAB _GEOMETRY/CRY_Data_File.txt");   // CRY_input_file_path
    }
    ~PrimaryGeneratorAction() { delete fParticleGun; }
    void GeneratePrimaries(G4Event* anEvent) {
        G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
        G4ParticleDefinition* particle = particleTable->FindParticle(parameters.particleName);
        fParticleGun->SetParticleDefinition(particle);
        G4double x0 = (parameters.Det_sizeXY + 4.0) * (G4UniformRand()-0.5); // x range -77.5 cm to +77.5 cm
        G4double y0 = (parameters.Det_sizeXY + 4.0) * (G4UniformRand()-0.5); // y range -77.5 cm to +77.5 cm
        G4double z0 = (parameters.RPCPosZ[0] +  parameters.Det_sizeZ + 2.0); // z position -65.0 cm
        
        //Input From Cry File
        G4int PID; G4double En, px, py, pz, xx,yy;
        infile >>PID>>En>>xx>>yy>>px>>py>>pz;
        fParticleGun->SetParticleEnergy(En*MeV);
        fParticleGun->SetParticleMomentumDirection(G4ThreeVector(px, py, pz));
        fParticleGun->SetParticlePosition(G4ThreeVector(x0*cm, y0*cm, z0*cm));
        fParticleGun->GeneratePrimaryVertex(anEvent);
    }
};

//Run.cc

class RunAction : public G4UserRunAction {
public:
    RunAction(): G4UserRunAction(){ }
    ~RunAction(){ }
    void BeginOfRunAction(const G4Run*){ }
    void EndOfRunAction(const G4Run* run){
        G4int nofEvents = run->GetNumberOfEvent();
        if (nofEvents == 0) return;
        if (IsMaster()) G4cout << G4endl << "------> " << ".------------------.End of Global Run.----------------------." << G4endl;
        else G4cout << G4endl << "--------------------End of Local Run------------------------" << G4endl;
    }
};

//Event.cc
class EventAction : public G4UserEventAction {
public:
    double RPCHits[6][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
    G4double SetKE = 0.0;
    RunAction* fRunAction;
    Parameters parameters;
    std::fstream outfile;
public:
    EventAction(RunAction* runAction) : G4UserEventAction(), fRunAction(runAction) { }
    ~EventAction() { }
    void BeginOfEventAction(const G4Event*) { }
    void EndOfEventAction(const G4Event* event) {
        G4int eventID = event->GetEventID();
        int per = 10000;
        if (eventID%per == 0) G4cout << "Total Events : " << eventID<< G4endl;
        if (RPCHits[0][2]!=0.0 && RPCHits[1][2]!=0.0 && RPCHits[2][2]!=0.0 && RPCHits[3][2]!=0.0 && RPCHits[4][2]!=0.0 && RPCHits[5][2]!=0.0){
            
            outfile.open(parameters.filename, std::ios::out | std::ios::app);
            outfile << " "
            <<RPCHits[0][0]<<" "<<RPCHits[0][1]<<" "<<parameters.RPCPosZ[0]*10.0<<" "
            <<RPCHits[1][0]<<" "<<RPCHits[1][1]<<" "<<parameters.RPCPosZ[1]*10.0<<" "
            <<RPCHits[2][0]<<" "<<RPCHits[2][1]<<" "<<parameters.RPCPosZ[2]*10.0<<" "
            <<RPCHits[3][0]<<" "<<RPCHits[3][1]<<" "<<parameters.RPCPosZ[3]*10.0<<" "
            <<RPCHits[4][0]<<" "<<RPCHits[4][1]<<" "<<parameters.RPCPosZ[4]*10.0<<" "
            <<RPCHits[5][0]<<" "<<RPCHits[5][1]<<" "<<parameters.RPCPosZ[5]*10.0<<" "
            <<G4endl;
            outfile.close();
            eventID=0;
            SetKE=0.0;
            for (G4int i=0; i<6; i++) for (G4int j=0; j<3; j++) RPCHits[i][j] = 0.0;
        }
    }
};
//Stepaction
class SteppingAction : public G4UserSteppingAction {
private:
    EventAction*  fEventAction;
    G4LogicalVolume* fScoringVolume;
public:
    SteppingAction(EventAction* eventAction): G4UserSteppingAction(), fEventAction(eventAction), fScoringVolume(0) { }
    ~SteppingAction() { }
    void UserSteppingAction(const G4Step* step) {
        if (!fScoringVolume) {
            const DetectorConstruction* detectorConstruction =
            static_cast<const DetectorConstruction*>(G4RunManager::GetRunManager()->GetUserDetectorConstruction());
            fScoringVolume = detectorConstruction->GetScoringVolume();
        }
        G4LogicalVolume* volume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume();
        
        if (volume == fScoringVolume) {
            G4ThreeVector pos = step->GetPreStepPoint()->GetPosition();
            G4double posX = pos.x(), posY = pos.y(), posZ = pos.z();
            G4int pid = step->GetTrack()->GetDefinition()->GetPDGEncoding();
            G4double KE = step->GetTrack()->GetDynamicParticle()->GetKineticEnergy();
            G4int copy_no = step->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber();
            if (pid ==13 || pid ==-13) {
                fEventAction->SetKE = KE;
                fEventAction->RPCHits[copy_no][0] = posX;
                fEventAction->RPCHits[copy_no][1] = posY;
                fEventAction->RPCHits[copy_no][2] = posZ;
            }
        }
    }
};
//action.cc
class ActionInitialization : public G4VUserActionInitialization {
public:
    ActionInitialization(): G4VUserActionInitialization() { }
    ~ActionInitialization() { }
    void BuildForMaster() const {
        RunAction* runAction = new RunAction();
        SetUserAction(runAction);
    }
    void Build() const {
        RunAction* runAction = new RunAction();
        SetUserAction(runAction);
        EventAction* eventAction = new EventAction(runAction);
        SetUserAction(eventAction);
        SteppingAction* steppingAction = new SteppingAction(eventAction);
        SetUserAction(steppingAction);
        PrimaryGeneratorAction* primaryGeneratorAction = new PrimaryGeneratorAction();
        SetUserAction(primaryGeneratorAction);
    }
};
//sim.cc

int main(int argc,char** argv){
    Parameters parameters;
    G4bool GUI = parameters.GUI;
    G4Random::setTheEngine(new CLHEP::RanecuEngine);
    auto* runManager = G4RunManagerFactory::CreateRunManager(G4RunManagerType::Default);
    runManager->SetNumberOfThreads(parameters.NumberOfThreads);
    DetectorConstruction* Detector = new DetectorConstruction();
    FTFP_BERT* Physics = new FTFP_BERT(0);
    ActionInitialization* Action = new ActionInitialization();
    runManager->SetUserInitialization(Detector);
    runManager->SetUserInitialization(Physics);
    runManager->SetUserInitialization(Action);

    if ( ! GUI ) {
        // Start a Run
        runManager->Initialize();
        int numberOfEvent = parameters.numberOfEvent;
        runManager->BeamOn(numberOfEvent);
        delete runManager;
    } else {
        // Alternative Run
        G4UIExecutive* ui = 0;
        if ( argc == 1 ) ui = new G4UIExecutive(argc, argv);
        G4VisManager* visManager = new G4VisExecutive;
        visManager->Initialize();
        G4UImanager* UImanager = G4UImanager::GetUIpointer();
        if ( ! ui ) {
            G4String command = "/control/execute ";
            G4String fileName = argv[1];
            UImanager->ApplyCommand(command+fileName);
        }
        else {
            UImanager->ApplyCommand("/control/execute init_vis.mac");
            ui->SessionStart();
            delete ui;
        }
        delete visManager;
        delete runManager;
    }
}


