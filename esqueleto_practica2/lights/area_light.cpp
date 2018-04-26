#include <lights/arealight.h>

Spectrum AreaLight::Sample(const gmtl::Point3f& position, gmtl::Vec3f& wi, float& pdf, gmtl::Rayf& visibilityRay) const
{

	//wi to be calculated. direction from a point of the area light to point
	// it is the inverse of the visibility ray

	// position of the light is the center of the light

    return Spectrum(0.0f);
}