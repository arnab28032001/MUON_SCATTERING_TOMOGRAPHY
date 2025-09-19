#include "action_project1.hh"
#include "event_project1.hh"
MyActionInitialization::MyActionInitialization()
{}
MyActionInitialization::~MyActionInitialization()
{}

void MyActionInitialization::Build() const
{
	 SetUserAction(new MyPrimaryGenerator());
    	SetUserAction(new MyRunAction());
    	SetUserAction(new MyEventAction()); 
}

