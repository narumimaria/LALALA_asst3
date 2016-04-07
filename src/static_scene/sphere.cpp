#include "sphere.h"

#include <cmath>

#include  "../bsdf.h"
#include "../misc/sphere_drawing.h"

namespace CMU462 { namespace StaticScene {

bool Sphere::test(const Ray& r, double& t1, double& t2) const {

  // TODO:
  // Implement ray - sphere intersection test.
  // Return true if there are intersections and writing the
  // smaller of the two intersection times in t1 and the larger in t2.
    double a, b, c;
    a = 1;
    b = 2 * dot(r.d, (r.o - o));
    c = (r.o - o).norm2() - r2;
    double delta = b * b - 4 * a * c;
    if (delta < 0) {
        return false;
    }else{
        t1 = - (b + sqrt(delta)) / (2 * a);
        t2 = - (b - sqrt(delta)) / (2 * a);
        if (t1 < 0 || t2 < 0) {
            return false;
        }
        return true;
    }
}

bool Sphere::intersect(const Ray& r) const {

  // TODO:
  // Implement ray - sphere intersection.
  // Note that you might want to use the the Sphere::test helper here.
    double t1, t2;
    if (!test(r, t1, t2)) {
        return false;
    }else{
        return true;
    }
}

bool Sphere::intersect(const Ray& r, Intersection *i) const {

  // TODO:
  // Implement ray - sphere intersection.
  // Note again that you might want to use the the Sphere::test helper here.
  // When an intersection takes place, the Intersection data should be updated
  // correspondingly.
    double t1, t2;
    if (!test(r, t1, t2)) {
        return false;
    }else{
        i->t = t1;
        i->n = normal(r.o + t1 * r.d);
        i->bsdf = get_bsdf();
        i->primitive = this;
        return true;
    }

}

void Sphere::draw(const Color& c) const {
  Misc::draw_sphere_opengl(o, r, c);
}

void Sphere::drawOutline(const Color& c) const {
    //Misc::draw_sphere_opengl(o, r, c);
}


} // namespace StaticScene
} // namespace CMU462
