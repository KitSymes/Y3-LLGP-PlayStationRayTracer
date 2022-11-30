#pragma once
#ifndef _TRACKERMANAGER_H
#define _TRACKERMANAGER_H

#include "Tracker.h"
#include <ult.h>

class TrackerManager
{
public:
	static TrackerManager& GetInstance();

	Tracker* GetDefaultTracker();
	Tracker* GetSphereTracker();
private:
	TrackerManager();
	~TrackerManager();

	Tracker _defaultTracker;
	Tracker _sphereTracker;

	SceUltWaitingQueueResourcePool waitingQueueResourcePool;
	void* workArea;
};
#endif
