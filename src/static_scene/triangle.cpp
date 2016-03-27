#include "triangle.h"

#include "CMU462/CMU462.h"
#include "GL/glew.h"

namespace CMU462 { namespace StaticScene {

Triangle::Triangle(const Mesh* mesh, size_t v1, size_t v2, size_t v3) :
    mesh(mesh), v1(v1), v2(v2), v3(v3) { }

BBox Triangle::get_bbox() const {
  
  // TODO: 
  // compute the bounding box of the triangle
  
  return BBox();
}

bool Triangle::intersect(const Ray& r) const {
  
  // TODO: implement ray-triangle intersection
    Vector3D p0 = mesh->positions[v1];
    Vector3D p1 = mesh->positions[v2];
    Vector3D p2 = mesh->positions[v3];
    Vector3D e1 = p1 - p0;
    Vector3D e2 = p2 - p0;
    Vector3D s = r.o - p0;
    
    double u, v, t;
    double ede = dot(cross(e1, r.d), e2);
    t = - (dot(cross(s, e2), e1)) / ede;
    if (t < r.min_t || t > r.max_t) {
        return false;
    }else {
        u = - (dot(cross(s, e2), r.d)) / ede;
        if (u < 0 || u > 1) {
            return false;
        }else {
            v = (dot(cross(e1, r.d), s)) / ede;
            if (v < 0 || v > 1) {
                return false;
            }else{
                return true;
            }
        }
    }
}

bool Triangle::intersect(const Ray& r, Intersection *isect) const {
  
  // TODO: 
  // implement ray-triangle intersection. When an intersection takes
  // place, the Intersection data should be updated accordingly
    Vector3D p0 = mesh->positions[v1];
    Vector3D p1 = mesh->positions[v2];
    Vector3D p2 = mesh->positions[v3];
    Vector3D e1 = p1 - p0;
    Vector3D e2 = p2 - p0;
    Vector3D s = r.o - p0;
    
    bool tri_intersect = false;
    double u, v, t;
    double ede = 1 / dot(cross(e1, r.d), e2);
    t = - (dot(cross(s, e2), e1)) * ede;
    if (t < r.min_t || t > r.max_t) {
        return false;
    }else {
        u = - (dot(cross(s, e2), r.d)) * ede;
        if (u < 0 || u > 1) {
            return false;
        }else {
            v = (dot(cross(e1, r.d), s)) * ede;
            if (v < 0 || v > 1) {
                return false;
            }else{
                if (u + v > 1) {
                    return false;
                }else{
                    tri_intersect = true;
                }
            }
        }
    }
    isect->t = t;
    isect->bsdf = get_bsdf();
    isect->primitive = this;
    isect->n = mesh->normals[v1] + u * mesh->normals[v2] + v * mesh->normals[v3];
    r.max_t = t;
    return tri_intersect;

}

void Triangle::draw(const Color& c) const {
  glColor4f(c.r, c.g, c.b, c.a);
  glBegin(GL_TRIANGLES);
  glVertex3d(mesh->positions[v1].x,
             mesh->positions[v1].y,
             mesh->positions[v1].z);
  glVertex3d(mesh->positions[v2].x,
             mesh->positions[v2].y,
             mesh->positions[v2].z);
  glVertex3d(mesh->positions[v3].x,
             mesh->positions[v3].y,
             mesh->positions[v3].z);
  glEnd();
}

void Triangle::drawOutline(const Color& c) const {
  glColor4f(c.r, c.g, c.b, c.a);
  glBegin(GL_LINE_LOOP);
  glVertex3d(mesh->positions[v1].x,
             mesh->positions[v1].y,
             mesh->positions[v1].z);
  glVertex3d(mesh->positions[v2].x,
             mesh->positions[v2].y,
             mesh->positions[v2].z);
  glVertex3d(mesh->positions[v3].x,
             mesh->positions[v3].y,
             mesh->positions[v3].z);
  glEnd();
}



} // namespace StaticScene
} // namespace CMU462
