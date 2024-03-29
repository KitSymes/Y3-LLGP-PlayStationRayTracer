#include <stdio.h>
#include <stdlib.h>
#include <scebase.h>
#include <kernel.h>
#include <gnmx.h>
#include <video_out.h>

#include <cstdio>
#include <cmath>
#include <fstream>
#include <vector>
#include <iostream>
#include <cassert>

#include <algorithm>
#include <sstream>
#include <string.h>

// Time precision
#include <chrono>

// Threading
#include <thread>

#include "allocator.h"

// Custom
#include "Structures.h"
#include "Sphere.h"
#include "TrackerManager.h"
#include <ult.h>
#include <scetypes.h>
#include <libsysmodule.h> 
#include <json.hpp>


static const size_t kOnionMemorySize = 64 * 1024 * 1024;

std::chrono::time_point<std::chrono::system_clock> start;
std::chrono::time_point<std::chrono::system_clock> end;
std::chrono::duration<double> total_elapsed_time;

static const int num_threads = 10;
#define TRACE_THREAD_PER_LINES 60
#define SUBDIVIDE_COUNT 3
#define THREADS 64 // 4 to the power of SUBDIVIDE_COUNT
#define SPLIT 8 // 2 to the power of SUBDIVIDE_COUNT
#define OCTREE true
#define COUNT 10

using namespace sce;
using namespace sce::Gnmx;

//[comment]
// This variable controls the maximum recursion depth
//[/comment]
#define MAX_RAY_DEPTH 5

float mix(const float& a, const float& b, const float& mix)
{
	return b * mix + a * (1 - mix);
}

//[comment]
// This is the main trace function. It takes a ray as argument (defined by its origin
// and direction). We test if this ray intersects any of the geometry in the scene.
// If the ray intersects an object, we compute the intersection point, the normal
// at the intersection point, and shade this point using this information.
// Shading depends on the surface property (is it transparent, reflective, diffuse).
// The function returns a color for the ray. If the ray intersects an object that
// is the color of the object at the intersection point, otherwise it returns
// the background color.
//[/comment]
Vec3f trace(
	const Vec3f& rayorig,
	const Vec3f& raydir,
	const std::vector<Sphere*>& spheres,
	const int& depth,
	Octree* oct)
{
	//if (raydir.length() != 1) std::cerr << "Error " << raydir << std::endl;
	float tnear = INFINITY;
	const Sphere* sphere = NULL;
#if OCTREE
	unsigned char x = 0;
	Vec3f dir = raydir;
	Vec3f orig = rayorig;

	if (raydir.x < 0.0)
	{
		orig.x = oct->GetCentre().x - rayorig.x;
		dir.x = -raydir.x;
		x |= 4;
	}
	if (raydir.y < 0.0)
	{
		orig.y = oct->GetCentre().y - rayorig.y;
		dir.y = -raydir.y;
		x |= 2;
	}
	if (raydir.z < 0.0)
	{
		orig.z = oct->GetCentre().z - rayorig.z;
		dir.z = -raydir.z;
		x |= 1;
	}

	float tx0 = (oct->xmin - orig.x) / dir.x;
	float tx1 = (oct->xmax - orig.x) / dir.x;
	float ty0 = (oct->ymin - orig.y) / dir.y;
	float ty1 = (oct->ymax - orig.y) / dir.y;
	float tz0 = (oct->zmin - orig.z) / dir.z;
	float tz1 = (oct->zmax - orig.z) / dir.z;

	if (std::max(std::max(tx0, ty0), tz0) < std::min(std::min(tx1, ty1), tz1))
		sphere = oct->Trace(rayorig, raydir, tx0, ty0, tz0, tx1, ty1, tz1, x, tnear);
#else

	// find intersection of this ray with the sphere in the scene
	for (unsigned i = 0; i < spheres.size(); ++i) {
		float t0 = INFINITY, t1 = INFINITY;
		if (spheres[i]->intersect(rayorig, raydir, t0, t1)) {
			if (t0 < 0) t0 = t1;
			if (t0 < tnear) {
				tnear = t0;
				sphere = spheres[i];
			}
		}
	}
#endif
	// if there's no intersection return black or background color
	if (!sphere) return Vec3f(2);
	Vec3f surfaceColor = 0; // color of the ray/surfaceof the object intersected by the ray
	Vec3f phit = rayorig + raydir * tnear; // point of intersection
	Vec3f nhit = phit - sphere->center; // normal at the intersection point
	nhit.normalize(); // normalize normal direction
	// If the normal and the view direction are not opposite to each other
	// reverse the normal direction. That also means we are inside the sphere so set
	// the inside bool to true. Finally reverse the sign of IdotN which we want
	// positive.
	float bias = 1e-4; // add some bias to the point from which we will be tracing
	bool inside = false;
	if (raydir.dot(nhit) > 0) nhit = -nhit, inside = true;
	if ((sphere->transparency > 0 || sphere->reflection > 0) && depth < MAX_RAY_DEPTH) {
		float facingratio = -raydir.dot(nhit);
		// change the mix value to tweak the effect
		float fresneleffect = mix(pow(1 - facingratio, 3), 1, 0.1);
		// compute reflection direction (not need to normalize because all vectors
		// are already normalized)
		Vec3f refldir = raydir - nhit * 2 * raydir.dot(nhit);
		refldir.normalize();
		Vec3f reflection = trace(phit + nhit * bias, refldir, spheres, depth + 1, oct);
		Vec3f refraction = 0;
		// if the sphere is also transparent compute refraction ray (transmission)
		if (sphere->transparency) {
			float ior = 1.1, eta = (inside) ? ior : 1 / ior; // are we inside or outside the surface?
			float cosi = -nhit.dot(raydir);
			float k = 1 - eta * eta * (1 - cosi * cosi);
			Vec3f refrdir = raydir * eta + nhit * (eta * cosi - sqrt(k));
			refrdir.normalize();
			refraction = trace(phit - nhit * bias, refrdir, spheres, depth + 1, oct);
		}
		// the result is a mix of reflection and refraction (if the sphere is transparent)
		surfaceColor = (
			reflection * fresneleffect +
			refraction * (1 - fresneleffect) * sphere->transparency) * sphere->surfaceColor;
	}
	else {
		// it's a diffuse object, no need to raytrace any further
		for (unsigned i = 0; i < spheres.size(); ++i) {
			if (spheres[i]->emissionColor.x > 0) {
				// this is a light
				Vec3f transmission = 1;
				Vec3f lightDirection = spheres[i]->center - phit;
				lightDirection.normalize();
				for (unsigned j = 0; j < spheres.size(); ++j) {
					if (i != j) {
						float t0, t1;
						if (spheres[j]->intersect(phit + nhit * bias, lightDirection, t0, t1)) {
							transmission = 0;
							break;
						}
					}
				}
				surfaceColor += sphere->surfaceColor * transmission *
					std::max(float(0), nhit.dot(lightDirection)) * spheres[i]->emissionColor;
			}
		}
	}

	return surfaceColor + sphere->emissionColor;
}

int32_t traceThreaded(uint64_t arg)
{
	TraceWrapper* wrapper = (TraceWrapper*)arg;
	for (int y = wrapper->startY; y < wrapper->imageHeight && y < wrapper->startY + wrapper->sizeY; y++) {
		for (int x = wrapper->startX; x < wrapper->imageWidth && x < wrapper->startX + wrapper->sizeX; x++) {
			float xx = (2 * ((x + 0.5) * wrapper->invWidth) - 1) * wrapper->angle * wrapper->aspectRatio;
			float yy = (1 - 2 * ((y + 0.5) * wrapper->invHeight)) * wrapper->angle;
			Vec3f raydir(xx, yy, -1);
			raydir.normalize();
			Vec3f pixel = trace(Vec3f(0), raydir, *wrapper->spheres, 0, wrapper->octree);
			wrapper->image[int(x + y * wrapper->imageWidth)] = pixel;
		}
	}
	return SCE_OK;
}

//[comment]
// Main rendering function. We compute a camera ray for each pixel of the image
// trace it and return a color. If the ray hits a sphere, we return the color of the
// sphere at the intersection point, else we return the background color.
//[/comment]
void render(const std::vector<Sphere*>& spheres, int iteration)
{

	//auto start = std::chrono::system_clock::now();

	// Initialize the WB_ONION memory allocator

	LinearAllocator onionAllocator;
	int ret = onionAllocator.initialize(
		kOnionMemorySize, SCE_KERNEL_WB_ONION,
		SCE_KERNEL_PROT_CPU_RW | SCE_KERNEL_PROT_GPU_ALL);

	//if (ret != SCE_OK)
	//	return ret;

	Octree* oct = new Octree(Vec3f());
	oct->Create(spheres, 50.0f);

	//unsigned width = 1920, height = 1080;
	unsigned width = 640, height = 480;
	size_t totalSize = sizeof(Vec3f) * width * height;

	void* buffer = onionAllocator.allocate(totalSize, Gnm::kAlignmentOfBufferInBytes);

	Vec3f* image = reinterpret_cast<Vec3f*>(buffer);
	Vec3f* pixel = image;

	float invWidth = 1 / float(width), invHeight = 1 / float(height);
	float fov = 30, aspectratio = width / float(height);
	float angle = tan(M_PI * 0.5 * fov / 180.);

	uint32_t maxNumUlthread = THREADS;
	uint32_t numWorkerThread = 3;
	void* runtimeWorkArea = malloc(sceUltUlthreadRuntimeGetWorkAreaSize(maxNumUlthread, numWorkerThread));
	SceUltUlthreadRuntime runtime;
	ret = sceUltUlthreadRuntimeCreate(&runtime, "Render Runtime", maxNumUlthread, numWorkerThread, runtimeWorkArea, NULL);
	assert(ret == SCE_OK);

	SceUltUlthread traceThreads[THREADS];
	TraceWrapper traceWrappers[THREADS];

	for (unsigned int i = 0; i < THREADS; i++)
	{
		traceWrappers[i].spheres = &spheres;
		traceWrappers[i].image = image;
		traceWrappers[i].sizeX = width / SPLIT;
		traceWrappers[i].sizeY = height / SPLIT;
		traceWrappers[i].startX = (i / SPLIT) * traceWrappers[i].sizeX;
		traceWrappers[i].startY = (i % SPLIT) * traceWrappers[i].sizeY;
		traceWrappers[i].imageWidth = width;
		traceWrappers[i].imageHeight = height;
		traceWrappers[i].invWidth = invWidth;
		traceWrappers[i].invHeight = invHeight;
		traceWrappers[i].aspectRatio = aspectratio;
		traceWrappers[i].angle = angle;
		traceWrappers[i].octree = oct;
		uint64_t arg = (uint64_t)&traceWrappers[i];

		sceUltUlthreadCreate(&traceThreads[i],
			"trace thread",
			traceThreaded,
			arg,
			NULL,
			NULL,
			&runtime,
			NULL);
	}

	for (unsigned int i = 0; i < THREADS; i++)
	{
		int32_t status;
		sceUltUlthreadJoin(&traceThreads[i], &status);
	}

	sceUltUlthreadRuntimeDestroy(&runtime);
	free(runtimeWorkArea);

	// Save result to a PPM image (keep these flags if you compile under Windows)
	std::stringstream ss;
	ss << "/app0/spheres" << iteration << ".ppm";
	std::string tempString = ss.str();
	char* filename = (char*)tempString.c_str();

	std::ofstream ofs(filename, std::ios::out | std::ios::binary);
	ofs << "P6\n" << width << " " << height << "\n255\n";
	for (unsigned i = 0; i < width * height; ++i) {
		ofs << (unsigned char)(std::min(float(1), image[i].x) * 255) <<
			(unsigned char)(std::min(float(1), image[i].y) * 255) <<
			(unsigned char)(std::min(float(1), image[i].z) * 255);
	}
	ofs.close();
	delete oct;
}

void LoadScene(std::vector<Sphere*>& spheres)
{
	std::ifstream f("/app0/scene2.json");
	nlohmann::json data = nlohmann::json::parse(f);

	for (nlohmann::json::iterator it = data["spheres"].begin(); it != data["spheres"].end(); ++it) {
		auto sphere = it.value();
		Vec3f centre = Vec3f(sphere["centre"][0], sphere["centre"][1], sphere["centre"][2]);
		Vec3f surfaceColor = Vec3f(sphere["surfaceColor"][0], sphere["surfaceColor"][1], sphere["surfaceColor"][2]);
		Vec3f emissionColor = Vec3f(sphere["emissionColor"][0], sphere["emissionColor"][1], sphere["emissionColor"][2]);
		spheres.push_back(new Sphere(centre,
			sphere["radius"],
			surfaceColor,
			sphere["reflection"],
			sphere["transparency"],
			emissionColor));
	}
}

void BasicRender(int iteration)
{
	std::vector<Sphere*> spheres;
	// Vector structure for Sphere (position, radius, surface color, reflectivity, transparency, emission color)

	LoadScene(spheres);
	//TrackerManager::GetInstance().GetDefaultTracker()->Verify();

	// This creates a file, titled 1.ppm in the current working directory
	render(spheres, iteration);

	for (Sphere* sphere : spheres)
		delete sphere;
}

int main(int argc, char** argv)
{
	int ret = SCE_OK;

	ret = sceSysmoduleLoadModule(SCE_SYSMODULE_ULT);
	assert(ret == SCE_OK);

	ret = sceUltInitialize();
	assert(ret == SCE_OK);

	srand(13);

	using std::chrono::high_resolution_clock;
	using std::chrono::duration_cast;
	using std::chrono::duration;
	using std::chrono::milliseconds;

	int count = COUNT;
	int total = 0;
	for (int i = 0; i < count; i++)
	{
		auto t1 = high_resolution_clock::now();
		BasicRender(i);
		auto t2 = high_resolution_clock::now();

		auto t_int = duration_cast<milliseconds>(t2 - t1);

		std::cout << t_int.count() << " ms seconds\n";
		total += t_int.count();
	}

	std::cout << (total / count) << " average ms seconds\n";

	return 0;
}

void* operator new(size_t size)
{
	size_t bytes = size + sizeof(Header) + sizeof(Footer);
	//std::cout << "new " << size << " reqeusted, allocating " << bytes << std::endl;
	char* pMem = (char*)malloc(bytes);

	Header* header = (Header*)pMem;
	header->size = size;
	header->next = nullptr;
	header->prev = nullptr;
	header->checkvalue = 0xDEAD;
	TrackerManager::GetInstance().GetDefaultTracker()->Add(header);

	Footer* footer = (Footer*)(pMem + sizeof(Header) + size);
	footer->checkvalue = 0xC0DE;

	void* pStartMemBlock = pMem + sizeof(Header);
	return pStartMemBlock;
}

void* operator new(size_t size, Tracker* tracker)
{
	size_t bytes = size + sizeof(Header) + sizeof(Footer);
	//std::cout << "new " << size << " reqeusted, allocating " << bytes << std::endl;
	char* pMem = (char*)malloc(bytes);

	Header* header = (Header*)pMem;
	header->size = size;
	header->next = nullptr;
	header->prev = nullptr;
	header->checkvalue = 0xDEAD;
	tracker->Add(header);

	Footer* footer = (Footer*)(pMem + sizeof(Header) + size);
	footer->checkvalue = 0xC0DE;

	void* pStartMemBlock = pMem + sizeof(Header);
	return pStartMemBlock;
}

void operator delete(void* pMem)
{
	Header* header = (Header*)((char*)pMem - sizeof(Header));
	Footer* footer = (Footer*)((char*)pMem + header->size);

	//std::cout << "freeing " << header->size << " (" << (header->size + sizeof(Header) + sizeof(Footer)) << ")" << std::endl;
	header->tracker->Remove(header);
	free(header);
}