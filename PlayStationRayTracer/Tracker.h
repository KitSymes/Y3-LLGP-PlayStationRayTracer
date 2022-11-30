#pragma once
#ifndef _TRACKER_H
#define _TRACKER_H

#define DEBUG false

#include <mutex>
struct Header;

class Tracker
{
public:
	Tracker();
	~Tracker();

	void Add(Header* header);
	void Remove(Header* header);

	void Verify(Header* header);

	size_t GetByteCount();
private:
	size_t _bytes;
	Header* _first;
	Header* _last;
	std::mutex m;

	void AddBytes(size_t bytes);
	void RemoveBytes(size_t bytes);
};
#endif
