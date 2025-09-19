//for project 1

#include<iostream>
#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4VisManager.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"
#include "construction_project1.hh"
#include "physics_project1.hh"
#include "action_project1.hh"
int main(int argc, char** argv)
{
	G4RunManager* runManager = new G4RunManager;
	runManager->SetUserInitialization(new MyDetectorConstruction());// for detector construction
	runManager->SetUserInitialization(new MyPhysicsList());// for Physics List
	runManager->SetUserInitialization(new MyActionInitialization());
	runManager->Initialize();
	//runManager->BeamOn(30*12648960);
	//runManager->BeamOn(200000);
	runManager->BeamOn(2*12648960);
	/*
	G4UIExecutive *ui =new G4UIExecutive(argc,argv);
	 
	G4VisManager *visManager = new G4VisExecutive();
	visManager-> Initialize();
	 
	G4UImanager *UImanager = G4UImanager::GetUIpointer();
	 
	 UImanager->ApplyCommand("/vis/open OGL");
	 UImanager->ApplyCommand("/vis/viewer/set/viewpointVector 1 1 1");
	 
	 //UImanager->ApplyCommand("/vis/viewer/set/targetPoint 0 0 0 cm");
	 //UImanager->ApplyCommand("/vis/viewer/set/eyePosition 0 0 600 cm");
	 //UImanager->ApplyCommand("/vis/viewer/panTo 0 0");
	 //UImanager->ApplyCommand("/vis/viewer/zoom 1");
	 //UImanager->ApplyCommand("/vis/viewer/set/style wireframe");
	 UImanager->ApplyCommand("/vis/drawVolume");
	 UImanager->ApplyCommand("vis/viewer/set/autoRefresh true");
	 UImanager->ApplyCommand("/vis/scene/add/trajectories smooth");
	 UImanager->ApplyCommand("/vis/scene/endOfEventAction accumulate"); 
	 UImanager->ApplyCommand("/vis/enable");
	 UImanager->ApplyCommand("/vis/viewer/flush");
	 UImanager->ApplyCommand("/vis/scene/add/axes 0 0 0 100 cm");
	 UImanager->ApplyCommand("/vis/viewer/update");
	 ui->SessionStart();
	 */
	 
	


	 
	return 0;
	
}
