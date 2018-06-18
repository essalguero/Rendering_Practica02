#ifndef __LAMBERT_H__
#define __LAMBERT_H__

#include <material.h>

class Lambert : public Material
{
public:
	Lambert(const char* name) : Material(name), Kd_color(0.0f), Kd(1.0f) {}
	Lambert(const char* name, float kd, const Color& kd_color) : Material(name), Kd(kd), Kd_color(kd_color) {}

	float Kd;		// Constante difusa
	Color Kd_color; // Color

	// color
	//Li Luz que entra
	//L direccion en la luz
	//V la dirección en la que se mira

    virtual Spectrum BRDF(const Spectrum& Li, const gmtl::Vec3f& L, const gmtl::Vec3f& V, const IntersectInfo& info) const
    {
        //return Spectrum(0.0f);
		return Li * Kd_color.GetColor(info) * (Kd / M_PI);
    }


	// Indirect light
    virtual bool Sample(gmtl::Vec3f& wi, float& pdf, const IntersectInfo& info) const
    {
		wi = info.ray.getDir();
		gmtl::Vec3f wo;

		pdf = this->pdf(wi, wo);

        return true;
    }


	// Solo sirve para multiple importance sampling
    virtual float pdf(const gmtl::Vec3f& wi, const gmtl::Vec3f& wo) const
    {
        //return 0.0f;
		//return 1.0f;

		return 1.0f / (2 * M_PI);
    }
};


#endif