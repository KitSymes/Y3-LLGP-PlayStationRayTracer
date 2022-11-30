#include "TrackerManager.h"
#include <iostream>

TrackerManager& TrackerManager::GetInstance()
{
	static TrackerManager instance;
	return instance;
}

Tracker* TrackerManager::GetDefaultTracker()
{
	if (DEBUG)
		std::cout << "Default Tracker: ";
	return &_defaultTracker;
}

Tracker* TrackerManager::GetSphereTracker()
{
	if (DEBUG)
		std::cout << "Sphere Tracker: ";
	return &_sphereTracker;
}

TrackerManager::TrackerManager()
{
	if (DEBUG)
		std::cout << "Tracker Manager Initialised" << std::endl;
}
