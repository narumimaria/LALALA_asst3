#include "bsdf.h"

#include <iostream>
#include <algorithm>
#include <utility>

using std::min;
using std::max;
using std::swap;

namespace CMU462 {

void make_coord_space(Matrix3x3& o2w, const Vector3D& n) {

    Vector3D z = Vector3D(n.x, n.y, n.z);
    Vector3D h = z;
    if (fabs(h.x) <= fabs(h.y) && fabs(h.x) <= fabs(h.z)) h.x = 1.0;
    else if (fabs(h.y) <= fabs(h.x) && fabs(h.y) <= fabs(h.z)) h.y = 1.0;
    else h.z = 1.0;

    z.normalize();
    Vector3D y = cross(h, z);
    y.normalize();
    Vector3D x = cross(z, y);
    x.normalize();

    o2w[0] = x;
    o2w[1] = y;
    o2w[2] = z;
}

// Diffuse BSDF //

Spectrum DiffuseBSDF::f(const Vector3D& wo, const Vector3D& wi) {
    double cos_theta = std::max(0.0, wi.z);
  return albedo * (1.0 / PI) * cos_theta;
}

Spectrum DiffuseBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {
    *wi = sampler.get_sample();
    *pdf = 1.0 / (2.0 * PI);
    double cos_theta = std::max(0.0, wi->z);
  return albedo * (1.0 / PI) * cos_theta;
}

// Mirror BSDF //

Spectrum MirrorBSDF::f(const Vector3D& wo, const Vector3D& wi) {
    if (wo.z == wi.z && (wo.x / wi.x == wo.y / wi.y)) {
        return reflectance;
    }else{
        return Spectrum(0,0,0);
    }
}

Spectrum MirrorBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {

  // TODO:
  // Implement MirrorBSDF
    reflect(wo, wi);
    *pdf = 1.0f;
    
  return reflectance;
}

// Glossy BSDF //

/*
Spectrum GlossyBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  return Spectrum();
}

Spectrum GlossyBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {
  *pdf = 1.0f;
  return reflect(wo, wi, reflectance);
}
*/

// Refraction BSDF //

Spectrum RefractionBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  return Spectrum();
}

Spectrum RefractionBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {

  // TODO:
  // Implement RefractionBSDF

  return Spectrum();
}

// Glass BSDF //

Spectrum GlassBSDF::f(const Vector3D& wo, const Vector3D& wi) {
    double Rs, Rp;
    if (wo.z == wi.z && wo.z > 0 && (wo.z != 1 && wi.z != 1 && wo.x / sqrt(1 - wo.z * wo.z) == - wi.x / sqrt(1 - wi.z * wi.z))) {
        // reflect
        double cos_theta_t = sqrt(1 - (1.0 / ior) * (1.0 / ior) * (1 - wo.z * wo.z));
        Rs = (wo.z- ior * cos_theta_t) /
             (wo.z + ior * cos_theta_t);
        Rp = (ior * wo.z - cos_theta_t) /
             (ior * wo.z + cos_theta_t);
        double Fr = 0.5 * ((Rs * Rs) + (Rp * Rp));
        return reflectance * Fr;
    }else {
        if (wo.z == 1 || wo.z == -1) {
            return reflectance;
        }else if (wo.x + wi.x * ior == 0 || wo.x * ior + wi.x == 0) {
            // refrect
            if (wo.z > 0) {
                // wo is the incoming ray
                double cos_theta_t = fabs(wi.z);
                Rs = (wo.z - ior * cos_theta_t) / (wo.z + ior * cos_theta_t);
                Rp = (ior * wo.z - cos_theta_t) / (ior * wo.z + cos_theta_t);
            }else {
                // wi is the incoming ray
                double cos_theta_t = fabs(wo.z);
                Rs = (wi.z - ior * cos_theta_t) / (wi.z + ior * cos_theta_t);
                Rp = (ior * wi.z - cos_theta_t) / (ior * wi.z + cos_theta_t);
                
            }
            
            double Fr = 0.5 * ((Rs * Rs) + (Rp * Rp));
            return reflectance * (1 - Fr);
//            return Spectrum(0,0,0);
            
        }else {
            return Spectrum(0,0,0);
        }
        
    }
}

Spectrum GlassBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {

  // TODO:
  // Compute Fresnel coefficient and either reflect or refract based on it.
    if (!refract(wo, wi, ior, pdf)) {
        reflect(wo, wi);
        *pdf = 1.0f;
        return reflectance;
//        return Spectrum(0,0,0);
    }else {
        // not total internal reflection
        return *pdf * reflectance;
//        return Spectrum(0,0,0);
        
    }
}

void BSDF::reflect(const Vector3D& wo, Vector3D* wi) {

  // TODO:
  // Implement reflection of wo about normal (0,0,1) and store result in wi.
    double cos_theta = wo.z;
    Vector3D N = Vector3D(0,0,1);
    *wi = 2 * cos_theta * N - wo;

}

bool BSDF::refract(const Vector3D& wo, Vector3D* wi, float ior, float* pdf) {

  // TODO:
  // Use Snell's Law to refract wo surface and store result ray in wi.
  // Return false if refraction does not occur due to total internal reflection
  // and true otherwise. When dot(wo,n) is positive, then wo corresponds to a
  // ray entering the surface through vacuum.
    if (wo.z > 0) {
        // ray entering the surface
        double cos_theta_t = sqrt(1 - (1.0 / ior) * (1.0 / ior) * (1 - wo.z * wo.z));
        double cos_theta_i = fabs(wo.z);
        double Rs, Rp;
        Rs = (cos_theta_i - ior * cos_theta_t) / (cos_theta_i + ior * cos_theta_t);
        Rp = (ior * cos_theta_i - cos_theta_t) / (ior * cos_theta_i + cos_theta_t);
        double Fr = 0.5 * ((Rs * Rs) + (Rp * Rp));
        if ((double)(std::rand()) / RAND_MAX < Fr) {
            // reflect
            reflect(wo, wi);
            *pdf = Fr;
            
        }else {
            // refract
            if (wo.z != 1 && wo.z != -1) {
                wi->z = - sqrt(1 - (1.0 / ior) * (1.0 / ior) * (1 - wo.z * wo.z));
                wi->x = - wo.x * 1.0 / ior;
                wi->y = - wo.y * 1.0 / ior;
            }else {
                wi->z = -wo.z;
                wi->x = 0;
                wi->y = 0;
            }
            *pdf = 1 - Fr;
            
        }
    }else {
        // ray leaving the surface
            // total internal reflectance
        if (sqrt(1 - wo.z * wo.z) > 1 / ior) {
            return false;
        }else {
            double cos_theta_i = sqrt(1 - ior * ior * (1 - wo.z * wo.z));
            double cos_theta_t = fabs(wo.z);
            double Rs, Rp;
            Rs = (cos_theta_i - ior * cos_theta_t) / (cos_theta_i + ior * cos_theta_t);
            Rp = (ior * cos_theta_i - cos_theta_t) / (ior * cos_theta_i + cos_theta_t);
            double Fr = 0.5 * ((Rs * Rs) + (Rp * Rp));
            if ((double)(std::rand()) / RAND_MAX < Fr) {
                // reflect
                reflect(wo, wi);
                *pdf = Fr;
            }else {
                // refract
                if (wo.z != 1 && wo.z != -1) {
                    wi->z = sqrt(1 - ior * ior * (1 - wo.z * wo.z));
                    wi->x = - wo.x * ior;
                    wi->y = - wo.y * ior;
                }else {
                    wi->z = -wo.z;
                    wi->x = 0;
                    wi->y = 0;
                }
                *pdf = 1 - Fr;
                
            }
        }
    }

  return true;

}

// Emission BSDF //

Spectrum EmissionBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  return Spectrum();
}

Spectrum EmissionBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {
  *wi  = sampler.get_sample(pdf);
  return Spectrum();
}

} // namespace CMU462
