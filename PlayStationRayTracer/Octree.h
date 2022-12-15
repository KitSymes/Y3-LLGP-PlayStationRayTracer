#pragma once
#ifndef _OCTREE_H
#define _OCTREE_H

#include <vector>
#include "Vec3.h"
#include "Sphere.h"

typedef Vec3<float> Vec3f;
class Octree
{
private:
	Octree* _children[2][2][2];
	std::vector<Sphere*> _leafObjects;
	Vec3f _centre;
public:
	Octree(Vec3f centre);
	~Octree();

	void Create(std::vector<Sphere*> spheres, float sideLength);
	Sphere* Trace(const Vec3f& rayorig, const Vec3f& raydir, float tx0, float ty0, float tz0, float tx1, float ty1, float tz1, unsigned char x, float& tnear);

	float sideLength;
	float xmin, xmax;
	float ymin, ymax;
	float zmin, zmax;

	int NextNode(float txm, int yzExit, float tym, int xzExot, float tzm, int xyExit);
	Octree* ChildFromLabel(int label);
	Vec3f GetCentre() { return _centre; }
};
#endif
