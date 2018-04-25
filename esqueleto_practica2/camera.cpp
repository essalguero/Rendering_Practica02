#include <camera.h>

gmtl::Rayf Camera::generateRay(float pixelX, float pixelY) {
	// ray generated
	gmtl::Rayf ray;

	// Camera position in camera coordinates system
	gmtl::Point3f cameraPosition(0, 0, 0);

	// Position of the point in the raster
	gmtl::Point3f rasterPositionPoint(pixelX, pixelY, 0);

	// Matrix to transform from the raster to the world coordinates system
	Transform mvp = mCameraToWorld * mRasterToCamera;

	// Transforms the points
	gmtl::Point3f cameraInWorld = mCameraToWorld(cameraPosition);
	gmtl::Point3f pointInWorld = mvp(rasterPositionPoint);

	// set the direction of the ray
	ray.setOrigin(cameraInWorld);
	gmtl::Vec3f pixelDir(pointInWorld - cameraInWorld);
	normalize(pixelDir);
	ray.setDir(pixelDir);


	return ray;
}


gmtl::Point2ui Camera::getResolution() const
{
	return mResolution;
}

void Camera::setOutputPath(const char* filename)
{
	mOutputPath = filename;
}

const char* Camera::getOutputPath() const
{
	return mOutputPath.c_str();
}