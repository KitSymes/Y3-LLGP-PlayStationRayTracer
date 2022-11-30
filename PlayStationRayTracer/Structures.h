#pragma once
#ifndef _STRUCTURES_H
#define _STRUCTURES_H

#include "Tracker.h"

struct Header
{
	size_t size;
	Header* prev;
	Header* next;
	Tracker* tracker;
	int checkvalue;
};

struct Footer
{
	int reserved;
	int checkvalue;
};
#endif