#pragma once
#ifndef _STRUCTURES_H
#define _STRUCTURES_H

#include <vector>
#include "Tracker.h"
#include "Sphere.h"

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

struct TraceWrapper
{
	const std::vector<Sphere*>* spheres;
	Vec3f* image;
	unsigned int startX;
	unsigned int startY;
	unsigned int sizeX;
	unsigned int sizeY;
	unsigned imageWidth;
	unsigned imageHeight;
	float invWidth;
	float invHeight;
	float aspectRatio;
	float angle;
};
#endif