#include "Tracker.h"
#include "Structures.h"
#include <iostream>

Tracker::Tracker()
{
	_bytes = 0;
	if (DEBUG)
		std::cout << "Tracker Initialised" << std::endl;
	_first = nullptr;
	_last = nullptr;
}

Tracker::~Tracker()
{
	if (!setup)
		return;

	sceUltMutexDestroy(&mutex);
}

void Tracker::Setup(const char* name, SceUltWaitingQueueResourcePool* waitingQueueResourcePool)
{
	if (setup)
		return;

	sceUltMutexCreate(&mutex, name, waitingQueueResourcePool, NULL);
	setup = true;
}

void Tracker::Add(Header* header)
{
	if (!setup)
		return;

	sceUltMutexLock(&mutex);
	header->tracker = this;
	AddBytes(header->size);

	if (header->checkvalue != 0xDEAD)
		std::cout << "Tracker#Add: Incorrect Header checkvalue: " << header->checkvalue << " not " << 0xDEAD << std::endl;
	if (header->next || header->prev)
		std::cout << "ERROR" << std::endl;

	if (!_first)
		_first = header;

	if (_last)
	{
		if (_last->checkvalue != 0xDEAD)
			std::cout << "ERROR" << std::endl;

		header->prev = _last;
		_last->next = header;
	}

	_last = header;

	if (header->next)
		if (header->next->checkvalue != 0xDEAD)
			Verify(header);
	if (header->prev)
		if (header->prev->checkvalue != 0xDEAD)
			Verify(header);

	sceUltMutexUnlock(&mutex);
}

void Tracker::Remove(Header* header)
{
	if (!setup)
		return;

	sceUltMutexLock(&mutex);
	if (header->checkvalue != 0xDEAD)
		std::cout << "Tracker#Remove: Incorrect Header checkvalue: " << header->checkvalue << " not " << 0xDEAD << std::endl;
	Footer* footer = (Footer*)(((char*)header) + sizeof(Header) + header->size);
	if (footer->checkvalue != 0xC0DE)
		std::cout << "Tracker#Remove: Incorrect Footer checkvalue: " << footer->checkvalue << " not " << 0xC0DE << std::endl;

	RemoveBytes(header->size);

	if (header->prev)
	{
		if (header->prev->checkvalue != 0xDEAD)
			Verify(header);

		if (header->next && header->next->checkvalue != 0xDEAD)
				Verify(header);

		/*if (header->next == nullptr)
		{
			header->prev->next = nullptr;
			bool b2 = header->prev->next != nullptr;
			if (b2)
				std::cout << "";
		}
		else*/
		header->prev->next = header->next;
	}

	if (header->next)
	{
		if (header->next->checkvalue != 0xDEAD)
			Verify(header);

		if (header->prev && header->prev->checkvalue != 0xDEAD)
			Verify(header);
		header->next->prev = header->prev;
	}

	if (header == _first)
		_first = header->next;

	if (header == _last)
		_last = header->prev;
	sceUltMutexUnlock(&mutex);
}

void Tracker::Verify(Header* header)
{
	Header** pointerToHeaderPointer = &_first;
	Header* temp;
	Header* prev;

	while (*pointerToHeaderPointer)
	{
		if ((*pointerToHeaderPointer)->checkvalue != 0xDEAD)
			std::cout << "Incorrect Header checkvalue: " << temp->checkvalue << " not " << 0xDEAD << std::endl;
		temp = *pointerToHeaderPointer;
		if (temp == header)
			std::cout << "Header Found" << std::endl;

		char* mem = (char*)temp;
		Footer* footer = (Footer*)(mem + sizeof(Header) + temp->size);
		if (footer->checkvalue != 0xC0DE)
			std::cout << "Incorrect Footer checkvalue: " << footer->checkvalue << " not " << 0xC0DE << std::endl;
		prev = temp;
		pointerToHeaderPointer = &(temp->next);

		if ((*pointerToHeaderPointer) && !(*pointerToHeaderPointer)->prev)
			std::cout << "ERROR Link Broke Forwards" << std::endl;
	}

	if (temp != _last)
	{
		std::cout << "Last pointer found is not the Last" << std::endl;


		Header** hPP = &_last;
		Header* h;
		Header* p;

		while (*hPP)
		{
			if ((*hPP)->checkvalue != 0xDEAD)
				std::cout << "Incorrect Header checkvalue: " << h->checkvalue << " not " << 0xDEAD << std::endl;
			h = *hPP;
			if (h == header)
				std::cout << "Header Found" << std::endl;

			char* mem = (char*)h;
			Footer* footer = (Footer*)(mem + sizeof(Header) + h->size);
			if (footer->checkvalue != 0xC0DE)
				std::cout << "Incorrect Footer checkvalue: " << footer->checkvalue << " not " << 0xC0DE << std::endl;
			p = h;
			hPP = &(h->prev);

			if (!(*hPP)->next)
				std::cout << "ERROR Link Broke Backwards" << std::endl;
		}

		if (*hPP != _first)
			std::cout << "First pointer found is not the First" << std::endl;
	}
}

size_t Tracker::GetByteCount()
{
	return _bytes;
}

void Tracker::AddBytes(size_t bytes)
{
	if (!setup)
		return;

	_bytes += bytes;
	if (DEBUG)
		std::cout << "Allocating " << bytes << std::endl;
}

void Tracker::RemoveBytes(size_t bytes)
{
	if (!setup)
		return;

	_bytes -= bytes;
	if (DEBUG)
		std::cout << "Freeing " << bytes << std::endl;
}
