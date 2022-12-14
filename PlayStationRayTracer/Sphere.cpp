#include "Sphere.h"
#include "Structures.h"
#include "TrackerManager.h"

bool Sphere::IntersectsLine(const Vec3f& rayorig, const Vec3f& raydir, float& distanceSq)
{
	// Closest point on the line to position
	Vec3f AB = raydir;
	Vec3f AP = center - rayorig;
	float lengthSqAB = AB.x * AB.x + AB.y * AB.y;
	float t = (AP.x * AB.x + AP.y * AB.y) / lengthSqAB;
	if (t < 0)
		t = 0;
	if (t > 1)
		t = 1;

	Vec3f closest = rayorig + AB * t;

	if ((closest - center).length2() > radius2)
		return false;

	distanceSq = closest.length2();

	return true;
}

bool Sphere::IntersectsPlane(Vec3f pointOnPlane, Vec3f planeNormal)
{
	return abs((center - pointOnPlane).dot(planeNormal)) < radius;
}

void* Sphere::operator new(size_t size)
{
	std::cout << "new Sphere" << std::endl;
	size_t bytes = size + sizeof(Header) + sizeof(Footer);
	char* pMem = (char*)malloc(bytes);

	Header* header = (Header*)pMem;
	header->size = size;
	header->next = nullptr;
	header->prev = nullptr;
	header->checkvalue = 0xDEAD;
	TrackerManager::GetInstance().GetSphereTracker()->Add(header);

	Footer* footer = (Footer*)(pMem + sizeof(Header) + size);
	footer->checkvalue = 0xC0DE;

	void* pStartMemBlock = pMem + sizeof(Header);
	return pStartMemBlock;
}

void Sphere::operator delete(void* pMem)
{
	std::cout << "delete Sphere" << std::endl;
	Header* header = (Header*)((char*)pMem - sizeof(Header));
	Footer* footer = (Footer*)((char*)pMem + header->size);

	header->tracker->Remove(header);
	free(header);
}
