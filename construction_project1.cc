#include "construction_project1.hh"
#include "G4VisAttributes.hh"
#include "G4Tubs.hh"
#include "G4Colour.hh"
#include "detector_project1.hh"
#include "G4SDManager.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4PSEnergyDeposit.hh"
#include "detector_project1.hh"

// for small detectors only
MyDetectorConstruction::MyDetectorConstruction()
{}
MyDetectorConstruction::~MyDetectorConstruction()
{}
G4VPhysicalVolume *MyDetectorConstruction :: Construct()
{
	G4NistManager *nist=G4NistManager :: Instance();
	//Materials
	G4Material* air = nist->FindOrBuildMaterial("G4_AIR");
	G4Material *detMat=nist->FindOrBuildMaterial("G4_Ar");
        G4Material* concrete = nist->FindOrBuildMaterial("G4_CONCRETE");
        G4Material* lead = nist->FindOrBuildMaterial("G4_Pb");
        G4Material* steel = nist->FindOrBuildMaterial("G4_STAINLESS-STEEL");
        G4Material* uranium = nist->FindOrBuildMaterial("G4_U");
        uranium = new G4Material("Modified_Uranium", 92., 238.03*g/mole, 19.5*g/cm3);
	lead    = new G4Material("Modified_Lead", 82., 207.2*g/mole, 11.34*g/cm3);
	steel   = new G4Material("Modified_Steel", 26., 55.845*g/mole, 7.5*g/cm3); // Approx iron-based
	air     = new G4Material("Modified_Air", 7., 14.01*g/mole, 0.0012*g/cm3);
        
        // Visual attributes
        G4VisAttributes* red = new G4VisAttributes(G4Colour::Red());
        G4VisAttributes* green = new G4VisAttributes(G4Colour::Green());
        G4VisAttributes* blue = new G4VisAttributes(G4Colour::Blue());
        G4VisAttributes* yellow=new G4VisAttributes(G4Colour::Yellow());
        G4VisAttributes* white=new G4VisAttributes(G4Colour::White());
        G4VisAttributes* gray = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5));
        red->SetForceSolid(true);
        green->SetForceSolid(true);
        blue->SetForceSolid(true);
        yellow->SetForceSolid(true);
        //gray->SetForceSolid(true);
        //white->SetForceSolid(true);
        
	/*construction of the box or the world box that will contain everything*/
	
	
	G4Box *solidWorld=new G4Box("solidWorld",1*m,1*m,1*m);/*box of length 1*2 m of each length of the cube, a cube of 200*200*200cm3 */
	G4LogicalVolume *logicWorld =new G4LogicalVolume(solidWorld,air,"logicWorld");/*This is a logical volume that actually makes the maerial*/
	G4VPhysicalVolume *physWorld=new G4PVPlacement(0, G4ThreeVector(0.,0.,0.),logicWorld, "physWorld",0,false,0,true);/*Places the logical volume and there is no Mother Volume*/
	
	/*Small detector Construction*/ 
	// Detector dimensions
	G4double detXY_small = 120.*cm;
	G4double detZ_small = 0.2*cm;

	// Material (can be plastic, silicon, etc. for real case; using Air here)
	
	G4double det_vert_gap=45.0*cm;
	G4Box* solidDetector1 = new G4Box("SmallDetector", detXY_small/2, detXY_small/2, detZ_small/2);
	G4LogicalVolume* logicDetector1_up = new G4LogicalVolume(solidDetector1,detMat, "SmallDetector_Up");
	G4LogicalVolume* logicDetector1_down = new G4LogicalVolume(solidDetector1,detMat, "SmallDetector_Down");
	logicDetector1_up->SetVisAttributes(white);
	logicDetector1_down->SetVisAttributes(white);
	//upper detector1
	new G4PVPlacement(0,G4ThreeVector(0.,0.,45.*cm),logicDetector1_up,"U1",logicWorld,false,0,true);
	//lower detector1
	new G4PVPlacement(0,G4ThreeVector(0.,0.,-45.0*cm),logicDetector1_down,"D1",logicWorld,false,0,true);
	//upper detector2
	new G4PVPlacement(0,G4ThreeVector(0.,0.,55.*cm),logicDetector1_up,"U2",logicWorld,false,0,true);
	//lower detector2
	new G4PVPlacement(0,G4ThreeVector(0.,0.,-55.0*cm),logicDetector1_down,"D2",logicWorld,false,0,true);
	//upper detector3
	new G4PVPlacement(0,G4ThreeVector(0.,0.,65.*cm),logicDetector1_up,"U3",logicWorld,false,0,true);
	//lower detector3
	new G4PVPlacement(0,G4ThreeVector(0.,0.,-65.0*cm),logicDetector1_down,"D3",logicWorld,false,0,true);
	
	// ROI - concrete cylinder
           G4Tubs* roiCyl = new G4Tubs("roi", 0, 25*cm, 25*cm, 0,360*deg);
           G4LogicalVolume* logicROI = new G4LogicalVolume(roiCyl, concrete, "ROI");
           logicROI->SetVisAttributes(gray);
           new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), logicROI, "ROI", logicWorld, false, 0);
        // ROI - uranium cylinder
            G4Tubs* uranium_box = new G4Tubs("uroi", 0, 5*cm, 15*cm, 0, 360*deg);
            G4LogicalVolume* uranium_logicROI = new G4LogicalVolume(uranium_box, uranium, "UROI");
            uranium_logicROI->SetVisAttributes(green);
            new G4PVPlacement(0, G4ThreeVector(0.*cm, -15.*cm, 0.*cm), uranium_logicROI, "UROI", logicWorld, false, 0);
        
	// ROI - lead(pb) cylinder
            G4Tubs* lead_box = new G4Tubs("pbroi", 0, 5*cm, 15*cm, 0, 360*deg);
            G4LogicalVolume* lead_logicROI = new G4LogicalVolume(lead_box, lead, "PBOI");
            lead_logicROI->SetVisAttributes(blue);
            new G4PVPlacement(0, G4ThreeVector(0.*cm, 15.*cm, 0.), lead_logicROI, "PBOI", logicWorld, false, 0);
	// ROI - stainless stell cylinder
            G4Tubs* ss_box = new G4Tubs("ssroi", 0, 5*cm, 15*cm, 0, 360*deg);
            G4LogicalVolume* ss_logicROI = new G4LogicalVolume(ss_box, steel, "SSROI");
            ss_logicROI->SetVisAttributes(red);
            new G4PVPlacement(0, G4ThreeVector(-15.*cm, 0.*cm, 0.), ss_logicROI, "SSROI", logicWorld, false, 0);
            
       // ROI - air cylinder
            G4Tubs* air_box = new G4Tubs("air_roi", 0, 5*cm, 15*cm, 0, 360*deg);
            G4LogicalVolume* air_logicROI = new G4LogicalVolume(air_box, air, "AIR_ROI");
            air_logicROI->SetVisAttributes(yellow);
            new G4PVPlacement(0, G4ThreeVector(15.*cm, 0.*cm, 0.*cm), air_logicROI, "AIR_ROI", logicWorld, false, 0);
	
	G4SDManager* sdManager = G4SDManager::GetSDMpointer();
	MySensitiveDetector* detectorSD = new MySensitiveDetector("MyDetectorSD");
	sdManager->AddNewDetector(detectorSD);
	logicDetector1_up->SetSensitiveDetector(detectorSD);
	logicDetector1_down->SetSensitiveDetector(detectorSD);
	
	

	
	
	
	return physWorld;
}

void MyDetectorConstruction::ConstructSDandField() {}



