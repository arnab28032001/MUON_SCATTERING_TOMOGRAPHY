#include "event_project1.hh"
#include "run_project1.hh"
#include <fstream>

std::map<G4int, std::vector<HitData>> gEventHits;
const G4int flushInterval = 1000; // flush every 1000 events
void MyEventAction::EndOfEventAction(const G4Event* event) {
    G4int eventID = event->GetEventID();
    if(eventID%1000==0)
    {
    G4cout<<eventID<<"is running"<<G4endl;
    }
    
    auto it = gEventHits.find(eventID);
    if (it != gEventHits.end()) {
     if (it->second.size() == 6) {
            for (const auto& hit : it->second) {
                gHitOutputFile
                    << hit.position.x() << " "
                    << hit.position.y() << " "
                    << hit.position.z() << " ";
            }
            gHitOutputFile << "\n";
        }
	gEventHits.erase(it);
        // Flush periodically to avoid data loss
        if ((eventID + 1) % 1 == 0) {
        //G4cout<<"Flusing is done"<<G4endl;
            gHitOutputFile.flush();
        }

        // Erase after writing to free memory
        //gEventHits.erase(it);
    }
    /*
    if (gEventHits.find(eventID) != gEventHits.end()&&gEventHits[eventID].size()==6) {
        for (const auto& hit : gEventHits[eventID]) {
            gHitOutputFile 
                << hit.position.x() << " "
                << hit.position.y() << " "
                << hit.position.z()<<" ";
        }
        gHitOutputFile << "\n"; // new line for next event
       
       if ((eventID + 1) % flushInterval == 0) {
            gHitOutputFile.flush();
        //gEventHits.erase(eventID); // clear after writing
    }
     gEventHits[eventID].clear();
}*/

}

