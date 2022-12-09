#include "TrackerManager.h"
#include <iostream>
#include <assert.h>
#include <scetypes.h>

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
	uint32_t numThreads = 16;
	uint32_t numSyncObjects = 16;
	size_t workAreaSize = sceUltWaitingQueueResourcePoolGetWorkAreaSize(numThreads, numSyncObjects);
	workArea = malloc(workAreaSize);
	sceUltWaitingQueueResourcePoolCreate(&waitingQueueResourcePool, "Thread Manager Waiting Queue", numThreads, numSyncObjects, workArea, NULL);
	_defaultTracker.Setup("Default Tracker", &waitingQueueResourcePool);
	_sphereTracker.Setup("Sphere Tracker", &waitingQueueResourcePool);
	//_defaultTracker = new Tracker("Default Tracker Mutex", &waitingQueueResourcePool);
	//_sphereTracker = new Tracker("Dummy Tracker Mutex", &waitingQueueResourcePool);
}

TrackerManager::~TrackerManager()
{
	int32_t ret = sceUltWaitingQueueResourcePoolDestroy(&waitingQueueResourcePool);
	assert(ret == SCE_OK);
	free(workArea);
	/*delete _sphereTracker;
	_sphereTracker = nullptr;
	delete _defaultTracker;
	_defaultTracker = nullptr;*/
}
