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

#include "lambert.h"

int g_RenderMaxDepth = 12;

extern int g_pixel_samples;

const float AMBIENT_INTENSITY = 0.005f;

const float GO_ON_PROBABILITY = 0.1f;

const int NUMBER_SAMPLES = 50;

const int HALTON_NUMBER_1 = 3;
const int HALTON_NUMBER_2 = 5;


Spectrum traceRay(World* world, Ray& ray, int recursivityDepth = 0);


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


float ruletaRusa()
{
	float n = static_cast<float>(rand());
	return n / static_cast<float>(RAND_MAX);
}

/*bool calculateNextRay(float probability)
{

	//std::cout << randValue << std::endl;

	if (randValue < probability)
		return true;

	return false;
}*/


gmtl::Rayf generateRay(gmtl::Point3f initialPosition, gmtl::Point3f finalPosition) {
	// ray generated
	gmtl::Rayf ray;

	// set the direction of the ray

	// Direction to the light being taken into account
	gmtl::Vec3f direction(finalPosition - initialPosition);
	normalize(direction);
	ray.setOrigin(initialPosition + direction);// *0.01f);
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


Spectrum directRadiance(World* world, Ray& ray, IntersectInfo info)
{
	Spectrum totalLight = Vector3f(0.0f);

	if (world->mLights.size() > 0)
	{
		// Elegir aleatoriamente una luz
		int lightNumber = rand() % world->mLights.size();

		// Samplear un punto de esa luz
		gmtl::Vec3f wi;
		float pdf;
		gmtl::Rayf visibilityRay;
		Spectrum lightSample = world->mLights.at(lightNumber)->Sample(info.position, wi, pdf, visibilityRay);

		//world->shadow(visibilityRay);



		// Calculate shadows
		IntersectInfo newIntersectInfo;
		gmtl::Rayf newRay = visibilityRay;
		newRay.setOrigin(newRay.getOrigin() + (newRay.getDir() * 0.01f));
		world->intersect(newIntersectInfo, newRay);

		if (!world->shadow(newRay))
		{
			info.material->Sample(wi, pdf, info);
			//totalLight += info.material->BRDF(lightSample, info.position, info.position, info);
			Spectrum intermediate = info.material->BRDF(lightSample, info.position, info.position, info);
			totalLight += intermediate / (newRay.getMaxLength() * newRay.getMaxLength());
		}


	}

	return totalLight;
}

Spectrum directRadiance2(World* world, Ray& ray, IntersectInfo info)
{
	if (world->mLights.size() > 0)
	{
		// Elegir aleatoriamente una luz
		int lightNumber = rand() % world->mLights.size();

		// Samplear un punto de esa luz
		gmtl::Vec3f wi;
		float pdf;
		gmtl::Rayf visibilityRay;
		Spectrum lightSample = world->mLights.at(lightNumber)->Sample(info.position, wi, pdf, visibilityRay);

		//world->shadow(visibilityRay);


		// Interseccón con la escena para comprobar si desde el punto, la luz se ve
		IntersectInfo newIntersectInfo;
		//visibilityRay.setOrigin(visibilityRay.getOrigin() + 0.001f * visibilityRay.getDir());
		world->intersect(newIntersectInfo, visibilityRay);

		// Si no Interseccion, la luz se ve
		if (newIntersectInfo.objectID == -1)
		{
			// Obtener emision de la luz
			//info.material->BRDF(colorSample, info.position, info.position, info) * max(gmtl::dot(info.normal, wi), 0.0f) / pdf;
			//Spectrum directLight = calculateDiffuseComponent(world, static_cast<PointLight*>(world->mLights.at(lightNumber)), newIntersectInfo);

			Spectrum directLight = world->mLights.at(lightNumber)->mColor * (world->mLights.at(lightNumber)->mIntensity / (visibilityRay.getMaxLength() * visibilityRay.getMaxLength()));

			/*gmtl::Vec3f distanceVector = info.position + wi;
			float squareLightDistanceToCollision = gmtl::lengthSquared(distanceVector);
			lightSample /= squareLightDistanceToCollision;
			info.material->Sample(wi, pdf, info);*/


			// Multiplicar por BRDF
			//directLight = directLight * newIntersectInfo.material->BRDF(directLight, info.position, info.position, info);
			directLight *= info.material->BRDF(lightSample, info.position, info.position, info);

			// Dividir por PDF
			gmtl::Vec3f wo;
			directLight = directLight / info.material->pdf(wi, wo);

			// Dividir por 1 / numero_luces
			directLight = (directLight / static_cast<float>(1 / world->mLights.size()));
			directLight = directLight / sqrt(wi[0] * wi[0] + wi[1] * wi[1] + wi[2] * wi[2]);

			return directLight;
		}

	}
	return Spectrum(0.0f);
}


Spectrum indirectRadiance(World* world, Ray& ray, IntersectInfo info, int recursivityDepth)
{
	// - Calcular valor aleatorio X
	float valorContinuacion = ruletaRusa();

	// - Si X < ruleta_rusa_p
	if (valorContinuacion < GO_ON_PROBABILITY)
	{
		Spectrum indirectLight = ((Lambert*)info.material)->Kd_color.GetColor(info);
		
		// - Calcular una dirección para siguiente rebote
		float r1 = ruletaRusa();
		float r2 = ruletaRusa();

		float anglePhi = 2 * M_PI * r1;
		//float angleTheta = acos(r2);

		float valueToSqrt = 1.0f - (r2 * r2);

		float x = cos(anglePhi) * (sqrt(valueToSqrt));
		float y = sin(anglePhi) * (sqrt(valueToSqrt));
		float z = r2;

		gmtl::Point3f finalPosition = gmtl::Point3f(x, y, z) - info.position;

		gmtl::Rayf indirectRay = generateRay(info.position, finalPosition);

		// - Calcular iluminación para ese rayo(traceRay)
		//indirectLight = indirectLight * traceRay(world, indirectRay, recursivityDepth + 1);
		indirectLight = traceRay(world, indirectRay, recursivityDepth + 1);

		// - Multiplicar iluminación por el BRDF
		indirectLight = indirectLight * info.material->BRDF(indirectLight, info.position, info.position, info);

		// - Multiplicar por el coseno
		indirectLight = indirectLight * gmtl::dot(indirectRay.getDir(), info.normal);

		// - Dividir iluminación por el PDF
		gmtl::Vec3f wi = indirectRay.getDir();
		gmtl::Vec3f wo;
		indirectLight = indirectLight / info.material->pdf(wi, wo);

		// - Dividir por ruleta_rusa_p
		indirectLight = indirectLight / GO_ON_PROBABILITY;

		return indirectLight;
	}
	return Spectrum(0.0f);
}


Spectrum traceRay(World* world, Ray& ray, int recursivityDepth)
{
	IntersectInfo info;

	world->intersect(info, ray);

	if (info.objectID != InvalidObjectID)
	{

		Spectrum totalLight = Vector3f(0.0f);

		/*//Light* light = world->mLights.at(0);
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
			world->intersect(shadowInfo, newRay);

			//world->intersect(shadowInfo, visibilityRay);

			if (!world->shadow(newRay))
			{
				info.material->Sample(wi, pdf, info);
				//totalLight += info.material->BRDF(colorSample, info.position, info.position, info);
				Spectrum intermediate = info.material->BRDF(colorSample, info.position, info.position, info);
				totalLight += intermediate / (newRay.getMaxLength() * newRay.getMaxLength());
			}
		}*/

		totalLight = directRadiance(world, ray, info) + indirectRadiance(world, ray, info, 1);



		return totalLight;

	}
	else
	{
		return Spectrum(0.0f, 0.0f, 0.0f);
	}
}



Spectrum traceRay2(World* world, Ray& ray, int recursivityDepth)
{
	// Intersección con la escena
	IntersectInfo info;

	world->intersect(info, ray);

	if (info.objectID != InvalidObjectID)
	{
		// Calcular iluminacion emitida
		// no hay materiales emisivos en la practica

		// Calcular iluminacion directa (directRadiance)
		Spectrum directLight = directRadiance(world, ray, info);
		//directLight = Spectrum(0.0f);

		// Calcular iluminacion indirecta (indirectRadiance)
		Spectrum indirectLight = indirectRadiance(world, ray, info, recursivityDepth);
		//indirectLight = Spectrum(0.0f);

		// Iluminacion total = emitida + directa + indirecta
		return directLight + indirectLight;
	}

	return Spectrum(0.0f);
}

Spectrum traceRay3(World* world, Ray& ray, int recursivityDepth = 0)
{
	IntersectInfo info;

	world->intersect(info, ray);

	if (info.objectID != InvalidObjectID)
	{

		Spectrum totalLight = Vector3f(0.0f);

		// Luz Emitida
		//if (info.material->)

		// Luz directa
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

			if (!world->shadow(newRay))
			{
				//info.material->Sample(wi, pdf, info);
				totalLight += info.material->BRDF(colorSample, info.position, info.position, info) * max(gmtl::dot(info.normal, wi), 0.0f) / pdf;
			}
		}

		// Luz indirecta
		float valorContinuacion = ruletaRusa();

		if (valorContinuacion < GO_ON_PROBABILITY)
		{
			Spectrum indirectLight = Vector3f(0.0f);

			//for (int rayNumber = 0; rayNumber < NUMBER_SAMPLES; rayNumber++)
			//{
				// Calcular dirección para el rebote
				// Trace new rays
				float r1 = ruletaRusa();
				float r2 = ruletaRusa();

				float anglePhi = 2 * M_PI * r1;
				//float angleTheta = acos(r2);

				float valueToSqrt = 1.0f - (r2 * r2);

				float x = cos(anglePhi) * (sqrt(valueToSqrt));
				float y = sin(anglePhi) * (sqrt(valueToSqrt));
				float z = r2;

				gmtl::Point3f finalPosition = gmtl::Point3f(x, y, z) - info.position;

				//std::cout << "Calcular indirecta (" << valorContinuacion << ")" << std::endl;
				gmtl::Rayf indirectRay = generateRay(info.position, finalPosition);

				float cosineValue = gmtl::dot(indirectRay.getDir(), info.normal);

				/*if (cosineValue < 0)
				{
					gmtl::Vec3f direction = indirectRay.getDir();
					direction[0] = -direction[0];
					direction[1] = -direction[1];
					direction[2] = -direction[2];
					indirectRay.setDir(direction);
				}
				cosineValue = gmtl::dot(indirectRay.getDir(), info.normal);*/

				indirectLight = traceRay(world, indirectRay, recursivityDepth + 1);
			//}

				
				// mat = info.material;
				Spectrum s = Spectrum(info.material->Ka_color.GetColor(info));
				indirectLight = indirectLight * info.material->BRDF(s, info.position, info.position, info);
				indirectLight[0] = (indirectLight[0] * cosineValue) / valorContinuacion;
				indirectLight[1] = (indirectLight[1] * cosineValue) / valorContinuacion;
				indirectLight[2] = (indirectLight[2] * cosineValue) / valorContinuacion;
			

			totalLight += indirectLight;
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

	srand(time(NULL));

#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < dimY; ++i)
	{
		for (int j = 0; j < dimX; ++j)
		{

			//std::cout << "Trace ray for pixel (" << i << ", " << j << ")" << std::endl;

			//Calcular rayo desde cámara a pixel
			Ray ray = camera->generateRay(j, i);

			Spectrum totalColor = Spectrum(0, 0, 0);

			for (int rayNumber = 0; rayNumber < NUMBER_SAMPLES; rayNumber++)
			{

				// Halton para emitir varios rayos por pixel
				float halton_1 = Halton(rayNumber, HALTON_NUMBER_1);
				float halton_2 = Halton(rayNumber, HALTON_NUMBER_2);

				//Calcular rayo desde cámara a pixel
				ray = camera->generateRay(j + halton_1, i + halton_2);

				// Calcular iluminación del rayo
				totalColor += traceRay(world, ray, 1);
			}

			
			totalColor = totalColor / static_cast<float>(NUMBER_SAMPLES);

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