#include "environment_light.h"

namespace CMU462 { namespace StaticScene {

EnvironmentLight::EnvironmentLight(const HDRImageBuffer* envMap)
    : envMap(envMap) {
  // TODO: initialize things here as needed
}

Spectrum EnvironmentLight::sample_L(const Vector3D& p, Vector3D* wi,
                                    float* distToLight,
                                    float* pdf) const {
  // TODO: Implement
    double Xi1 = (double)(std::rand()) / RAND_MAX;
    Xi1 = 2.0 * Xi1 - 1.0;;
    double Xi2 = (double)(std::rand()) / RAND_MAX;
    
    double theta = acos(Xi1);
    double phi = 2.0 * PI * Xi2;
    
    wi->x = sinf(theta) * cosf(phi);
    wi->y = sinf(theta) * sinf(phi);
    wi->z = cosf(theta);
    
    *pdf = 1.0 / (4.0 * PI);
    
    *distToLight = INFINITY;

  return sample_dir(*wi);
}

Spectrum EnvironmentLight::sample_dir(const Vector3D& d) const {
    
  // TODO: Implement
    double theta = acos(d.z);
    double phi = atan(d.y / d.x);
    double map_x = std::abs((phi / (2.0 * PI))* envMap->w);
    double map_y = std::abs((PI - theta) / PI * envMap->h);
    
    size_t x = floor(map_x);
    size_t y = floor(map_y);
    
    Spectrum L;
    L.r = envMap->data[x     + y       * envMap->w].r * (x + 1 - map_x) * (y + 1 - map_y) +
          envMap->data[x     + (y + 1) * envMap->w].r * (x + 1 - map_x) * (map_y - y    ) +
          envMap->data[x + 1 + y       * envMap->w].r * (map_x - x    ) * (y + 1 - map_y) +
          envMap->data[x + 1 + (y + 1) * envMap->w].r * (map_x - x    ) * (map_y - y    );
    L.g = envMap->data[x     + y       * envMap->w].g * (x + 1 - map_x) * (y + 1 - map_y) +
          envMap->data[x     + (y + 1) * envMap->w].g * (x + 1 - map_x) * (map_y - y    ) +
          envMap->data[x + 1 + y       * envMap->w].g * (map_x - x    ) * (y + 1 - map_y) +
          envMap->data[x + 1 + (y + 1) * envMap->w].g * (map_x - x    ) * (map_y - y    );
    L.b = envMap->data[x     + y       * envMap->w].b * (x + 1 - map_x) * (y + 1 - map_y) +
          envMap->data[x     + (y + 1) * envMap->w].b * (x + 1 - map_x) * (map_y - y    ) +
          envMap->data[x + 1 + y       * envMap->w].b * (map_x - x    ) * (y + 1 - map_y) +
          envMap->data[x + 1 + (y + 1) * envMap->w].b * (map_x - x    ) * (map_y - y    );

    
  return L;
}

} // namespace StaticScene
} // namespace CMU462
