//#include <tbb/tbb.h>
#include <stdio.h>
#include <math.h>
#include <gmtl/gmtl.h>
#include <imageio.h>

#include <world.h>

#include <parsers/ass_parser.h>

#include <reporter.h>

#define WIN32_LEAN_AND_MEAN
#include <Windows.h>
#include <algorithm>

#include <main.h>


#include <standard.h>

#include "lights\pointlight.h"
#include "lights\arealight.h"

#include <cmath>

#include <iostream>

#include <time.h>

int g_RenderMaxDepth = 12;

extern int g_pixel_samples;

const float AMBIENT_INTENSITY = 0.005f;

const float GO_ON_PROBABILITY = 0.1f;

World* ReadFromFile(const char* filename)
{
	World* world;
	if (!ReadAssFile(filename, world))
	{
		Error("Could not read file: %s", filename);
		return NULL;
	}
	return world;
}



float Halton(int index, int base)
{
	float result = 0.0f;
	float f = 1.0f / base;
	int i = index;
	while (i > 0)
	{
		result += f * (i % base);
		i = i / base;
		f = f / base;
	}
	return result;
}


bool calculateNextRay(float probability)
{

	float randValue = static_cast<float>(rand()) / static_cast<float>(RAND_MAX);

	//std::cout << randValue << std::endl;

	if (randValue < probability)
		return true;

	return false;
}


gmtl::Rayf generateRay(gmtl::Point3f initialPosition, gmtl::Point3f finalPosition) {
	// ray generated
	gmtl::Rayf ray;

	// set the direction of the ray

	// Direction to the light being taken into account
	gmtl::Vec3f direction(finalPosition - initialPosition);
	normalize(direction);
	ray.setOrigin(initialPosition + direction * 0.01f);
	//ray.setOrigin(initialPosition);
	ray.setDir(direction);

	return ray;
}



Spectrum calculateDiffuseComponent(World* world, PointLight* light, IntersectInfo& intersectInfo)
{
	Standard* material = (Standard *)intersectInfo.material;
	Spectrum colorDifuso = material->Kd.GetColor(intersectInfo);

	Point3f lightPosition = light->getWorldPosition();

	Vector3f lightVector = Vector3f(lightPosition[0], lightPosition[1], lightPosition[2]);
	Vector3f intersectVector = Vector3f(intersectInfo.position[0], intersectInfo.position[1], intersectInfo.position[2]);
	Vector3f distanceVector = intersectVector - lightVector;
	float squareLightDistanceToCollision = distanceVector[0] * distanceVector[0] + distanceVector[1] * distanceVector[1] + distanceVector[2] * distanceVector[2];

	IntersectInfo infoShadow;
	Ray shadowRay = generateRay(intersectInfo.position, lightPosition);
	gmtl::Vec3f direction(lightPosition - intersectInfo.position);
	world->intersect(infoShadow, shadowRay, gmtl::length(direction));
	if (infoShadow.objectID == InvalidObjectID)
	{
		Vector3f directionLight = lightPosition - intersectInfo.position;

		Vector3f normalNormalized = intersectInfo.normal;
		normalize(normalNormalized);

		normalize(directionLight);

		float productoEscalar = gmtl::dot(normalNormalized, directionLight);

		float valorCoseno = productoEscalar > 0.0f ? productoEscalar : 0.0f;

		return Spectrum(colorDifuso[0] * (light->mIntensity / squareLightDistanceToCollision) * valorCoseno,
			colorDifuso[1] * (light->mIntensity / squareLightDistanceToCollision) * valorCoseno,
			colorDifuso[2] * (light->mIntensity / squareLightDistanceToCollision) * valorCoseno);
	}

	return Spectrum(0.0f, 0.0f, 0.0f);

}


Spectrum calculateSpecularComponent(World* world, PointLight* light, IntersectInfo& intersectInfo)
{
	Standard* material = (Standard *)intersectInfo.material;
	Spectrum colorEspecular = material->Ks.GetColor(intersectInfo);

	Vector3f directionNormalized = -intersectInfo.ray.getDir();
	normalize(directionNormalized);

	Point3f lightPosition = light->getWorldPosition();
	Vector3f lightVector = Vector3f(lightPosition[0], lightPosition[1], lightPosition[2]);
	Vector3f intersectVector = Vector3f(intersectInfo.position[0], intersectInfo.position[1], intersectInfo.position[2]);
	Vector3f distanceVector = intersectVector - lightVector;
	float squareLightDistanceToCollision = distanceVector[0] * distanceVector[0] + distanceVector[1] * distanceVector[1] + distanceVector[2] * distanceVector[2];

	IntersectInfo infoShadow;
	Ray shadowRay = generateRay(intersectInfo.position, lightPosition);
	gmtl::Vec3f direction(lightPosition - intersectInfo.position);
	world->intersect(infoShadow, shadowRay, gmtl::length(direction));
	if (infoShadow.objectID == InvalidObjectID)
	{
		Vector3f directionLight = lightPosition - intersectInfo.position;

		Vector3f normalNormalized = intersectInfo.normal;
		normalize(normalNormalized);

		normalize(directionLight);

		float escalarProduct = (gmtl::dot(normalNormalized, directionLight));
		Vector3f vectorInterm1 = escalarProduct * normalNormalized;
		Vector3f vectorInterm2 = 2.0f * vectorInterm1;
		Vector3f rVector = vectorInterm2 - directionLight;
		normalize(rVector);

		float productoEscalar = gmtl::dot(rVector, directionNormalized);

		float valorCoseno = productoEscalar > 0.0f ? pow(productoEscalar, material->Kshi) : 0.0f;


		return Spectrum(colorEspecular[0] * (light->mIntensity / squareLightDistanceToCollision) * valorCoseno,
			colorEspecular[1] * (light->mIntensity / squareLightDistanceToCollision) * valorCoseno,
			colorEspecular[2] * (light->mIntensity / squareLightDistanceToCollision) * valorCoseno);
	}

	return Spectrum(0.0f, 0.0f, 0.0f);

}



gmtl::Vec3f refractedDirection(float ni, float nt, const gmtl::Vec3f& V, const gmtl::Vec3f& N)
{
	gmtl::Vec3f T;
	float eta;

	eta = ni / nt;
	float c1 = -dot(V, N);
	float c2_op = 1.0f - eta * eta*(1.0f - c1 * c1);
	if (c2_op < 0.0f)
		return gmtl::Vec3f(0.0f);

	float c2 = sqrt(c2_op);
	T = eta * V + (eta*c1 - c2)*N;

	return T;
}

gmtl::Vec3f transmittedDirection(bool& entering, const gmtl::Vec3f& N, const gmtl::Vec3f& V, const Standard* mat)
{
	gmtl::Vec3f normal_refraction = N;
	bool exitingObject = dot(N, -V) < 0.0f;
	float ni = 0.0f;
	float nt = 0.0f;

	if (exitingObject)
	{
		ni = mat->refractionIndex;
		nt = 1.0003f; // air refraction index
		normal_refraction = -normal_refraction;
	}
	else
	{
		ni = 1.0003f; // air refraction index
		nt = mat->refractionIndex;
	}

	gmtl::Vec3f T = refractedDirection(ni, nt, V, normal_refraction);

	entering = !exitingObject;
	return T;
}




// - Intersección con la escena
// - Calcular luminación emitida
// - Calcular luminación directa
// - Calcular luminación indirecta
// - iluminación total = emitida + directa + indirecta


// traceRay()
// - Intersección con la escena
// - Calcular luminación emitida
// - Calcular valor aleatorio X
// - Si X < ruleta_rusa_p
	// - Calcular una dirección para siguiente rebote
	// - Calcular iluminación para ese rayo(traceRay)
	// - Multiplicar iluminación por el BRDF
	// - Multiplicar por el coseno
	// - Dividir iluminación por el PDF
// - Iluminación total = emitida + rebote
// - Devolver iluminación total calculada


Spectrum traceRay(World* world, Ray& ray, int recursivityDepth = 0)
{
	IntersectInfo info;

	world->intersect(info, ray);

	if (info.objectID != InvalidObjectID)
	{

		Spectrum totalLight = Vector3f(0.0f);

		//Light* light = world->mLights.at(0);
		for (auto light = world->mLights.begin(); light != world->mLights.end(); ++light)
		{
			gmtl::Vec3f wi;
			float pdf;
			gmtl::Rayf visibilityRay;

			//std::cout << "calculating lights" << std::endl;

			//AreaLight* areaLight = (AreaLight*)(*light);
			//const gmtl::Point3f& position, gmtl::Vec3f& wi, float& pdf, gmtl::Rayf& visibilityRay
			Spectrum colorSample = (*light)->Sample(info.position, wi, pdf, visibilityRay);

			//gmtl::Vec3f distanceVector = info.position + wi;

			//float squareLightDistanceToCollision = gmtl::lengthSquared(distanceVector);
			// Not using emissive objects
			//totalLight += calculateLe(world, areaLight, info);

			//std::cout << "Calculated light" << std::endl;

			//colorSample /= squareLightDistanceToCollision;

			// Calculate shadows
			IntersectInfo shadowInfo;
			gmtl::Rayf newRay = visibilityRay;
			newRay.setOrigin(newRay.getOrigin() + (newRay.getDir() * 0.01f));
			//world->intersect(shadowInfo, newRay);

			//world->intersect(shadowInfo, visibilityRay);
			
			if (calculateNextRay(GO_ON_PROBABILITY))
			{
				//std::cout << "Continue" << std::endl;
				continue;
			}

			if (!world->shadow(newRay))
			{
				//info.material->Sample(wi, pdf, info);
				totalLight += info.material->BRDF(colorSample, info.position, info.position, info) * max(gmtl::dot(info.normal, wi), 0.0f) / pdf;
			}
		}

		

		return totalLight;

	}
	else
	{
		return Spectrum(0.0f, 0.0f, 0.0f);
	}
}

void render_image(World* world, unsigned int dimX, unsigned int dimY, float* image, float* alpha)
{
	Camera* camera = world->getCamera();

	unsigned int acc = 0;
#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < dimY; ++i)
	{
		for (int j = 0; j < dimX; ++j)
		{

			//std::cout << "Trace ray for pixel (" << i << ", " << j << ")" << std::endl;

			//Calcular rayo desde cámara a pixel
			Ray ray = camera->generateRay(j, i);


			// Halton para emitir varios rayos por pixel

			// Probabilidad para seguir trazando rayos por el mismo pixel

			// Calcular iluminación del rayo
			Spectrum totalColor = traceRay(world, ray, 1);

			// Guardar iluminación para el pixel en imagen
			image[(i * dimX * 3) + (j * 3)] = totalColor[0];
			image[(i * dimX * 3) + (j * 3) + 1] = totalColor[1];
			image[(i * dimX * 3) + (j * 3) + 2] = totalColor[2];

		}

		acc++;
		printf("\r%f", (float)acc / dimY);
	}


	/*for (int i = 0; i < dimX * dimY * 3; i += 3)
	{
	pixelX = i % dimX;
	pixelY = i / dimX;
	Ray ray = camera->generateRay(pixelX, pixelY);
	IntersectInfo info;
	world->intersect(info, ray);
	if (info.objectID != InvalidObjectID)
	{
	image[i] = FLT_MAX;
	image[i + 1] = FLT_MAX;
	image[i + 2] = FLT_MAX;
	}
	else
	{
	image[i] = 0.0f;
	image[i + 1] = 0.0f;
	image[i + 2] = 0.0f;
	//alpha[(i * (dimX * 3)) + j + 3] = 1;
	}
	}*/


	/*for (int i = 0; i < dimY; ++i)
	{
	for (int j = 0; j < dimX ; ++j)
	{
	Ray ray = camera->generateRay(j, i);
	IntersectInfo info;
	world->intersect(info, ray);
	if (info.objectID != InvalidObjectID)
	{
	image[(i * (dimX * 3)) + j] = FLT_MAX;
	//image[(i * (dimX * 3)) + j + 1] = FLT_MAX;
	//image[(i * (dimX * 3)) + j + 2] = FLT_MAX;
	//alpha[(i * (dimX * 3)) + j + 3] = 1;
	}
	else
	{
	image[(i * (dimX * 3)) + j] = 0.0f;
	//image[(i * (dimX * 3)) + j + 1] = 0.0f;
	//image[(i * (dimX * 3)) + j + 2] = 0.0f;
	//alpha[(i * (dimX * 3)) + j + 3] = 1;
	}
	}
	}*/
}

unsigned int g_intersectTriangleCalls;
extern "C"
{
	__declspec(dllexport) void renderSceneFromFile(float*& image, float*& alpha, World*& world, const char* filename)
	{
		google::InitGoogleLogging("rendering.dll");
		FLAGS_logtostderr = 1;

		g_intersectTriangleCalls = 0;

		// Create world from file
		world = ReadFromFile(filename);
		if (!world)
		{
			fprintf(stderr, "Error reading file %s. Press enter to exit", filename);
			getc(stdin);
			return;
		}
		INITREPORTER("report.ma", world);
		unsigned int dimX = world->getCamera()->getResolution()[0];
		unsigned int dimY = world->getCamera()->getResolution()[1];

		image = new float[dimX*dimY * 3];
		alpha = new float[dimX*dimY];


		// Compute pixel values
		clock_t tStart = clock();
		render_image(world, dimX, dimY, image, alpha);
		clock_t tEnd = clock();
		LOG(INFO) << "Time taken: " << (double)(tEnd - tStart) / CLOCKS_PER_SEC << "s";
		LOG(INFO) << "Triangles intersected: " << g_intersectTriangleCalls;

		google::ShutdownGoogleLogging();
	}

	__declspec(dllexport) void WriteImg(const char* name, float *pixels, float *alpha, int xRes,
		int yRes, int totalXRes, int totalYRes, int xOffset, int yOffset)
	{
		WriteImage(name, pixels, alpha, xRes, yRes, totalXRes, totalYRes, xOffset, yOffset);
	}
}

// dllmain.cpp : Defines the entry point for the DLL application.

BOOL APIENTRY DllMain(HMODULE hModule,
	DWORD  ul_reason_for_call,
	LPVOID lpReserved
)
{
	switch (ul_reason_for_call)
	{
	case DLL_PROCESS_ATTACH:
	case DLL_THREAD_ATTACH:
	case DLL_THREAD_DETACH:
	case DLL_PROCESS_DETACH:
		break;
	}
	return TRUE;
}