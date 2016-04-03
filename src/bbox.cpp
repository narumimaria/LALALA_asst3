#include "bbox.h"

#include "GL/glew.h"

#include <algorithm>
#include <iostream>

namespace CMU462 {

bool BBox::intersect(const Ray& r, double& t0, double& t1) const {

  // TODO:
  // Implement ray - bounding box intersection test
  // If the ray intersected the bouding box within the range given by
  // t0, t1, update t0 and t1 with the new intersection times.
    
    bool itsect = false;
    double  x_tmin = (min.x - r.o.x) / r.d.x;
    if ((r.o + x_tmin * r.d).y >= min.y &&
        (r.o + x_tmin * r.d).y <= max.y &&
        (r.o + x_tmin * r.d).z >= min.z &&
        (r.o + x_tmin * r.d).z <= max.z) {
        t0 = x_tmin;
        itsect = true;
    }
    double y_tmin = (min.y - r.o.y) / r.d.y;
    if ((r.o + y_tmin * r.d).x >= min.x &&
        (r.o + y_tmin * r.d).x <= max.x &&
        (r.o + y_tmin * r.d).z >= min.z &&
        (r.o + y_tmin * r.d).z <= max.z) {
        if (itsect) {
            if (y_tmin > t0) {
                t1 = y_tmin;
            }else {
                double t = t0;
                t0 = y_tmin;
                t1 = t;
            }
            return itsect;
        }else{
            t0 = y_tmin;
            itsect = true;
        }
    }
    double z_tmin = (min.z - r.o.z) / r.d.z;
    if ((r.o + z_tmin * r.d).y >= min.y &&
        (r.o + z_tmin * r.d).y <= max.y &&
        (r.o + z_tmin * r.d).x >= min.x &&
        (r.o + z_tmin * r.d).x <= max.x) {
        if (itsect) {
            if (z_tmin > t0) {
                t1 = z_tmin;
            }else {
                double t = t0;
                t0 = z_tmin;
                t1 = t;
            }
            return itsect;
        }else {
            t0 = z_tmin;
            itsect = true;
        }
    }
    
    double  x_tmax = (max.x - r.o.x) / r.d.x;
    if ((r.o + x_tmax * r.d).y >= min.y &&
        (r.o + x_tmax * r.d).y <= max.y &&
        (r.o + x_tmax * r.d).z >= min.z &&
        (r.o + x_tmax * r.d).z <= max.z) {
        if (itsect) {
            if (x_tmax > t0) {
                t1 = x_tmax;
            }else {
                double t = t0;
                t0 = x_tmax;
                t1 = t;
            }
            return itsect;
        }else {
            t0 = x_tmax;
            itsect = true;
        }
    }
    double y_tmax = (max.y - r.o.y) / r.d.y;
    if ((r.o + y_tmax * r.d).x >= min.x &&
        (r.o + y_tmax * r.d).x <= max.x &&
        (r.o + y_tmax * r.d).z >= min.z &&
        (r.o + y_tmax * r.d).z <= max.z) {
        if (itsect) {
            if (y_tmax > t0) {
                t1 = y_tmax;
            }else {
                double t = t0;
                t0 = y_tmax;
                t1 = t;
            }
            return itsect;
        }else {
            t0 = y_tmax;
            itsect = true;
        }
    }
    double z_tmax = (max.z - r.o.z) / r.d.z;
    if ((r.o + z_tmax * r.d).y >= min.y &&
        (r.o + z_tmax * r.d).y <= max.y &&
        (r.o + z_tmax * r.d).x >= min.x &&
        (r.o + z_tmax * r.d).x <= max.x) {
        if (itsect) {
            if (z_tmax > t0) {
                t1 = z_tmax;
            }else {
                double t = t0;
                t0 = z_tmax;
                t1 = t;
            }
            return itsect;
        }else {
            t0 = z_tmax;
            itsect = true;
        }
    }
    
    return itsect;
  
}

void BBox::draw(Color c) const {

  glColor4f(c.r, c.g, c.b, c.a);

	// top
	glBegin(GL_LINE_STRIP);
	glVertex3d(max.x, max.y, max.z);
  glVertex3d(max.x, max.y, min.z);
  glVertex3d(min.x, max.y, min.z);
  glVertex3d(min.x, max.y, max.z);
  glVertex3d(max.x, max.y, max.z);
	glEnd();

	// bottom
	glBegin(GL_LINE_STRIP);
  glVertex3d(min.x, min.y, min.z);
  glVertex3d(min.x, min.y, max.z);
  glVertex3d(max.x, min.y, max.z);
  glVertex3d(max.x, min.y, min.z);
  glVertex3d(min.x, min.y, min.z);
	glEnd();

	// side
	glBegin(GL_LINES);
	glVertex3d(max.x, max.y, max.z);
  glVertex3d(max.x, min.y, max.z);
	glVertex3d(max.x, max.y, min.z);
  glVertex3d(max.x, min.y, min.z);
	glVertex3d(min.x, max.y, min.z);
  glVertex3d(min.x, min.y, min.z);
	glVertex3d(min.x, max.y, max.z);
  glVertex3d(min.x, min.y, max.z);
	glEnd();

}

std::ostream& operator<<(std::ostream& os, const BBox& b) {
  return os << "BBOX(" << b.min << ", " << b.max << ")";
}

} // namespace CMU462
