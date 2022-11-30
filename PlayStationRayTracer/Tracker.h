#pragma once
#ifndef _TRACKER_H
#define _TRACKER_H

#define DEBUG false

struct Header;
#include <ult.h>

class Tracker
{
public:
	Tracker();
	~Tracker();

	void Setup(const char* name, SceUltWaitingQueueResourcePool* waitingQueueResourcePool);
	void Add(Header* header);
	void Remove(Header* header);

	void Verify(Header* header);

	size_t GetByteCount();
private:
	size_t _bytes;
	Header* _first;
	Header* _last;
	SceUltMutex mutex;
	bool setup = false;

	void AddBytes(size_t bytes);
	void RemoveBytes(size_t bytes);
};
#endif
