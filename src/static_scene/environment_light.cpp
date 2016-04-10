#include "environment_light.h"

namespace CMU462 { namespace StaticScene {

EnvironmentLight::EnvironmentLight(const HDRImageBuffer* envMap)
    : envMap(envMap) {
  // TODO: initialize things here as needed
        double* fp = new double[envMap->h * envMap->w];
        
        // compute p(theta, phi)
        double illum_sum = 0;
        for (size_t i = 0; i < envMap->w; i++) {
            for (size_t j = 0; j < envMap->h; j++) {
                Spectrum s = envMap->data[i + j * envMap->w];
                double brightness = s.illum();
                fp[i + j * envMap->w] = brightness;
                illum_sum += brightness;
            }
        }
        double scale = 1.0 / illum_sum;
        for (size_t i = 0; i < envMap->w; i++) {
            for (size_t j = 0; j < envMap->h; j++) {
                fp[i + j * envMap->w] *= scale;
            }
        }
        
        // compute p(theta) (marginal probability) and its cumulative math function
        fmp.resize(envMap->h);
        for (size_t j = 0; j < envMap->h; j++) {
            fmp[j] = 0;
            // compute p(theta)
            double ptheta = 0;
            for (size_t i = 0; i < envMap->w; i++) {
                ptheta += fp[i + j * envMap->w];
            }
            // compute cumulative math function
            if (j != 0) {
                fmp[j] = fmp[j-1] + ptheta;
            }else {
                fmp[j] = ptheta;
            }
        }
        
        // compute p(phi|theta) (conditional probability) and its cumulative math function
        fcp.resize(envMap->w * envMap->h);
        for (size_t i = 0; i < envMap->w; i++) {
            for (size_t j = 0; j < envMap->h; j++) {
                fcp[i + j * envMap->w] = 0;
                // compute p(phi|theta)
                double pphitheta = 0;
                if (j != 0) {
                    pphitheta = fp[i + j * envMap->w] / (fmp[j] - fmp[j-1]);
                    
                    // compute cumulative math function
                    if (i == 0) {
                        fcp[i + j * envMap->w] = fcp[j * envMap->w] + pphitheta;
                    }else{
                        fcp[i + j * envMap->w] = fcp[i + j * envMap->w - 1] + pphitheta;
                    }
                    
                    
                }else {
                    pphitheta = fp[i] / fmp[0];
                    
                    // compute cumulative math function
                    if (i == 0) {
                        fcp[i + j * envMap->w] = pphitheta;
                    }else {
                        fcp[i + j * envMap->w] = fcp[i - 1] + pphitheta;
                    }
                    
                }
            }
        }
        
}
    
Spectrum EnvironmentLight::sample_L(const Vector3D& p, Vector3D* wi,
                                    float* distToLight,
                                    float* pdf) const {
  // TODO: Implement
    
    // uniform sampling
//    double Xi1 = (double)(std::rand()) / RAND_MAX;
//    Xi1 = 2.0 * Xi1 - 1.0;;
//    double Xi2 = (double)(std::rand()) / RAND_MAX;
//    
//    double theta = acos(Xi1);
//    double phi = 2.0 * PI * Xi2;
//    
//    wi->x = sinf(theta) * cosf(phi);
//    wi->y = sinf(theta) * sinf(phi);
//    wi->z = cosf(theta);
//
//    *pdf = 1.0 / (4.0 * PI);
    
    // importance sampling
    double X = (double)(std::rand()) / RAND_MAX;
    double Y = (double)(std::rand()) / RAND_MAX;
    
    // search theta
    size_t idx_y = std::lower_bound(fmp.begin(), fmp.end(), Y) - fmp.begin() + 1;
    
    // search phi
    std::vector<double> fcp_phi;
    fcp_phi.resize(envMap->w);
    for (size_t i = 0; i < envMap->w; i++) {
        fcp_phi[i] = fcp[i + idx_y * envMap->w];
    }
    size_t idx_x = std::lower_bound(fcp_phi.begin(), fcp_phi.end(), X) - fcp_phi.begin() + 1;
    
    double theta = idx_y / envMap->h * PI;
    double phi = idx_x / envMap->w * 2.0 * PI;
    
    wi->x = sinf(theta) * cosf(phi);
    wi->y = sinf(theta) * sinf(phi);
    wi->z = cosf(theta);
    
    *pdf = (fmp[idx_y] - fmp[idx_y-1]) * (fcp_phi[idx_x] - fcp_phi[idx_x-1]);
    
    *distToLight = INFINITY;

  return envMap->data[idx_x + idx_y * envMap->w];
}

Spectrum EnvironmentLight::sample_dir(const Vector3D& d) const {
    
  // TODO: Implement
    // uniform sampling
    double theta = acos(d.y);
    Vector2D xz = Vector2D(d.z, d.x);
    xz /= xz.norm();
    double phi = acos(dot(xz,Vector2D(1,0)));
    if (d.y < 0) {
        phi = 2.0 * PI - phi;
    }
    double map_x = std::abs((phi / (2.0 * PI))* envMap->w);
    double map_y = std::abs(theta / PI * envMap->h);
    
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
