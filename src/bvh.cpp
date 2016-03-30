#include "bvh.h"

#include "CMU462/CMU462.h"
#include "static_scene/triangle.h"

#include <iostream>
#include <stack>

using namespace std;

namespace CMU462 { namespace StaticScene {

BVHAccel::BVHAccel(const std::vector<Primitive *> &_primitives,
                   size_t max_leaf_size) {

  this->primitives = _primitives;

  // TODO:
  // Construct a BVH from the given vector of primitives and maximum leaf
  // size configuration. The starter code build a BVH aggregate with a
  // single leaf node (which is also the root) that encloses all the
  // primitives.
    
  BBox bb;
  for (size_t i = 0; i < primitives.size(); ++i) {
    bb.expand(primitives[i]->get_bbox());
  }

  root = new BVHNode(bb, 0, primitives.size());
    
    stack<BVHNode*> bvhstack;
    bvhstack.push(root);
    const int B = 30;
    double bestX, bestY, bestZ;
    double costX, costY, costZ;
    
    BVHNode* curBVH;
    do
    {
        curBVH = bvhstack.top();
        bvhstack.pop();
        bb = curBVH->bb;
        
        double Sn = (bb.max.x - bb.min.x) * (bb.max.y - bb.min.y) +
                    (bb.max.x - bb.min.x) * (bb.max.z - bb.min.z) +
                    (bb.max.z - bb.min.z) * (bb.max.y - bb.min.y);
        double Sl, Sr;
        bool isCost;
        
        // x axis partitioning
        BBox bestXBl, bestXBr;
        isCost = false;
        double dx = bb.extent.x / (double)B;
        double curX = bb.min.x;
        bestX = curX;
        costX = 0;
        for (size_t i = 0; i < 29; i++) {
            double Nl = 0, Nr = 0;
            BBox BBl, BBr;
            curX += dx;
            Sl = (curX - bb.min.x) * (bb.max.y - bb.min.y) +
                 (curX - bb.min.x) * (bb.max.z - bb.min.z) +
                 (bb.max.z - bb.min.z) * (bb.max.y - bb.min.y);
            Sr = (bb.max.x - curX) * (bb.max.y - bb.min.y) +
                 (bb.max.x - curX) * (bb.max.z - bb.min.z) +
                 (bb.max.z - bb.min.z) * (bb.max.y - bb.min.y);
            for (size_t i = 0; i < curBVH->range; i++) {
                double centroid = primitives[curBVH->start + i]->get_bbox().centroid().x;
                if (centroid < curX) {
                    BBl.expand(primitives[curBVH->start + i]->get_bbox());
                    Nl++;
                }else {
                    BBr.expand(primitives[curBVH->start + i]->get_bbox());
                    Nr++;
                }
            }
            if (Nl !=0 && Nr != 0 && isCost == false) {
                costX = Sl * Nl / Sn + Sr * Nr / Sn;
                bestX = curX;
                bestXBl = BBl;
                bestXBr = BBr;
                isCost = true;
            }else{
                if (Sl * Nl / Sn + Sr * Nr / Sn < costX && Nl !=0 && Nr != 0) {
                    costX = Sl * Nl / Sn + Sr * Nr / Sn;
                    bestX = curX;
                    bestXBl = BBl;
                    bestXBr = BBr;
                }
            }
        }
        
        // y axis partitioning
        isCost = false;
        BBox bestYBl, bestYBr;
        double dy = bb.extent.y / (double)B;
        double curY = bb.min.y;
        bestY = curY;
        costY = 0;
        for (size_t i = 0; i < 29; i++) {
            double Nl = 0, Nr = 0;
            BBox BBl, BBr;
            curY += dy;
            Sl = (bb.max.x - bb.min.x) * (curY - bb.min.y) +
                 (bb.max.x - bb.min.x) * (bb.max.z - bb.min.z) +
                 (bb.max.z - bb.min.z) * (curY - bb.min.y);
            Sr = (bb.max.x - bb.min.x) * (bb.max.y - curY) +
                 (bb.max.x - bb.min.x) * (bb.max.z - bb.min.z) +
                 (bb.max.z - bb.min.z) * (bb.max.y - curY);
            for (size_t i = 0; i < curBVH->range; i++) {
                double centroid = primitives[curBVH->start + i]->get_bbox().centroid().y;
                if (centroid < curY) {
                    BBl.expand(primitives[curBVH->start + i]->get_bbox());
                    Nl++;
                }else {
                    BBr.expand(primitives[curBVH->start + i]->get_bbox());
                    Nr++;
                }
            }
            if (Nl !=0 && Nr != 0 && isCost == false) {
                costY = Sl * Nl / Sn + Sr * Nr / Sn;
                bestY = curY;
                bestYBl = BBl;
                bestYBr = BBr;
                isCost = true;
            }else{
                if (Sl * Nl / Sn + Sr * Nr / Sn < costY && Nl !=0 && Nr != 0) {
                    costY = Sl * Nl / Sn + Sr * Nr / Sn;
                    bestY = curY;
                    bestYBl = BBl;
                    bestYBr = BBr;
                }
            }
        }
        
        // z axis partitioning
        isCost = false;
        BBox bestZBl, bestZBr;
        double dz = bb.extent.z / (double)B;
        double curZ = bb.min.z;
        bestZ = curZ;
        costZ = 0;
        for (size_t i = 0; i < 29; i++) {
            double Nl = 0, Nr = 0;
            BBox BBl, BBr;
            curZ += dz;
            Sl = (bb.max.x - bb.min.x) * (bb.max.y - bb.min.y) +
                 (bb.max.x - bb.min.x) * (curZ - bb.min.z) +
                 (curZ - bb.min.z) * (bb.max.y - bb.min.y);
            Sr = (bb.max.x - bb.min.x) * (bb.max.y - bb.min.y) +
                 (bb.max.x - bb.min.x) * (bb.max.z - curZ) +
                 (bb.max.z - curZ) * (bb.max.y - bb.min.y);
            for (size_t i = 0; i < curBVH->range; i++) {
                double centroid = primitives[curBVH->start + i]->get_bbox().centroid().z;
                if (centroid < curZ) {
                    BBl.expand(primitives[curBVH->start + i]->get_bbox());
                    Nl++;
                }else {
                    BBr.expand(primitives[curBVH->start + i]->get_bbox());
                    Nr++;
                }
            }
            if (Nl !=0 && Nr != 0 && isCost == false) {
                costZ = Sl * Nl / Sn + Sr * Nr / Sn;
                bestZ = curZ;
                bestZBl = BBl;
                bestZBr = BBr;
                isCost = true;
            }else{
                if (Sl * Nl / Sn + Sr * Nr / Sn < costZ && Nl !=0 && Nr != 0) {
                    costZ = Sl * Nl / Sn + Sr * Nr / Sn;
                    bestZ = curZ;
                    bestZBl = BBl;
                    bestZBr = BBr;
                }
            }
        }
        
        // set current BVH's left child and right child
        int primlength = curBVH->range;
        if ((costX != 0 && costX <= costY && costX <= costZ) ||
            (costX != 0 && ((costZ == 0 && costX <= costY) || (costY == 0 && costX <= costZ))) ||
            (costX != 0 && costY == 0 && costZ == 0)) {
            
            Primitive* helpprim[primlength];
            int primcountl = 0, primcountr = 0;
            for (size_t i = 0; i < curBVH->range; i++) {
                double centroid = (primitives[curBVH->start + i]->get_bbox().centroid().x);
                if (centroid < bestX) {
                    helpprim[primcountl] = primitives[curBVH->start + i];
                    primcountl++;
                }else {
                    helpprim[primlength - primcountr - 1] = primitives[curBVH->start + i];
                    primcountr++;
                }
            }
            for (size_t i = 0; i < primlength; i++) {
                primitives[curBVH->start + i] = helpprim[i];
            }
            curBVH->l = new BVHNode(bestXBl, curBVH->start, primcountl);
            curBVH->r = new BVHNode(bestXBr, curBVH->start + primcountl, primcountr);
            
        }else if ((costY != 0 && costY <= costX && costY <= costZ) ||
                  (costY != 0 && ((costZ == 0 && costY <= costX) || (costX == 0 && costY <= costZ))) ||
                  (costY != 0 && costX == 0 && costZ == 0)) {
            
            Primitive* helpprim[primlength];
            int primcountl = 0, primcountr = 0;
            for (size_t i = 0; i < curBVH->range; i++) {
                double centroid = (primitives[curBVH->start + i]->get_bbox().centroid().y);
                if (centroid < bestY) {
                    helpprim[primcountl] = primitives[curBVH->start + i];
                    primcountl++;
                }else {
                    helpprim[primlength - primcountr - 1] = primitives[curBVH->start + i];
                    primcountr++;
                }
            }
            for (size_t i = 0; i < primlength; i++) {
                primitives[curBVH->start + i] = helpprim[i];
            }
            curBVH->l = new BVHNode(bestYBl, curBVH->start, primcountl);
            curBVH->r = new BVHNode(bestYBr, curBVH->start + primcountl, primcountr);
            
        }else if ((costZ != 0 && costZ <= costX && costZ <= costY) ||
                  (costZ != 0 && ((costX == 0 && costZ <= costY) || (costY == 0 && costZ <= costX))) ||
                  (costZ != 0 && costX == 0 && costY == 0)) {
            
            Primitive* helpprim[primlength];
            int primcountl = 0, primcountr = 0;
            for (size_t i = 0; i < curBVH->range; i++) {
                double centroid = (primitives[curBVH->start + i]->get_bbox().centroid().z);
                if (centroid < bestZ) {
                    helpprim[primcountl] = primitives[curBVH->start + i];
                    primcountl++;
                }else {
                    helpprim[primlength - primcountr - 1] = primitives[curBVH->start + i];
                    primcountr++;
                }
            }
            for (size_t i = 0; i < primlength; i++) {
                primitives[curBVH->start + i] = helpprim[i];
            }
            curBVH->l = new BVHNode(bestZBl, curBVH->start, primcountl);
            curBVH->r = new BVHNode(bestZBr, curBVH->start + primcountl, primcountr);
        }
        
        // push left child and right child to the stack
        if (curBVH->r->range > max_leaf_size) {
            bvhstack.push(curBVH->r);
        }
        if (curBVH->l->range > max_leaf_size) {
            bvhstack.push(curBVH->l);
        }
        
    }
    while (!bvhstack.empty());
    
}

BVHAccel::~BVHAccel() {

  // TODO:
  // Implement a proper destructor for your BVH accelerator aggregate

}

BBox BVHAccel::get_bbox() const {
  return root->bb;
}

bool BVHAccel::intersect(const Ray &ray) const {

  // TODO:
  // Implement ray - bvh aggregate intersection test. A ray intersects
  // with a BVH aggregate if and only if it intersects a primitive in
  // the BVH that is not an aggregate.
  bool hit = false;
    
//  for (size_t p = 0; p < primitives.size(); ++p) {
//    if(primitives[p]->intersect(ray)) hit = true;
//  }
    
    bool hitlbbox = false, hitrbbox = false;
    double t0, t1;
    if (root->bb.intersect(ray, t0, t1)) {
        stack<ClosestHit*> bvhstack_front;
        stack<ClosestHit*> bvhstack_back;
        ClosestHit* curBVH = new ClosestHit(root, t0);
        bvhstack_front.push(curBVH);
        ClosestHit* curHit = new ClosestHit(root, t1);
        while (!bvhstack_front.empty() || !bvhstack_back.empty()) {
            if (!bvhstack_front.empty()) {
                curBVH = bvhstack_front.top();
                bvhstack_front.pop();
            }else {
                curBVH = bvhstack_back.top();
                bvhstack_back.pop();
                if (curBVH->min_t > curHit->min_t) {
                    continue;
                }
            }
            
            
            if (curBVH->closestbvh->isLeaf()) {
                for (size_t j = 0; j < curBVH->closestbvh->range; j++) {
                    if (primitives[curBVH->closestbvh->start + j]->intersect(ray)) {
                        hit = true;
                    }
                }
            }else {
                
                double tl0, tl1, tr0, tr1;
                tl0 = t0; tr0 = t0; tl1 = t1; tr1 = t1;
                hitlbbox = curBVH->closestbvh->l->bb.intersect(ray, tl0, tl1);
                hitrbbox = curBVH->closestbvh->r->bb.intersect(ray, tr0, tr1);
                
                if (hitlbbox && hitrbbox) {
                    if (tl0 < tr0) {
                        ClosestHit* leftHit = new ClosestHit(curBVH->closestbvh->l, tl0);
                        ClosestHit* rightHit = new ClosestHit(curBVH->closestbvh->r, tr0);
                        bvhstack_front.push(leftHit);
                        bvhstack_back.push(rightHit);
                        t0 = tl0;
                        t1 = tl1;
                    }else {
                        ClosestHit* leftHit = new ClosestHit(curBVH->closestbvh->l, tl0);
                        ClosestHit* rightHit = new ClosestHit(curBVH->closestbvh->r, tr0);
                        bvhstack_back.push(leftHit);
                        bvhstack_front.push(rightHit);
                        t0 = tr0;
                        t1 = tr1;
                    }
                }
                if (hitrbbox && !hitlbbox) {
                    ClosestHit* rightHit = new ClosestHit(curBVH->closestbvh->r, tr0);
                    bvhstack_front.push(rightHit);
                }
                if (hitlbbox && !hitrbbox) {
                    ClosestHit* leftHit = new ClosestHit(curBVH->closestbvh->l, tl0);
                    bvhstack_front.push(leftHit);
                }
                if (!hitlbbox && !hitrbbox && bvhstack_front.empty() && bvhstack_back.empty()) {
                    return hit;
                }
            }
        }
        
        
    }else {
        return hit;
    }
    
    return hit;

}

bool BVHAccel::intersect(const Ray &ray, Intersection *i) const {

  // TODO:
  // Implement ray - bvh aggregate intersection test. A ray intersects
  // with a BVH aggregate if and only if it intersects a primitive in
  // the BVH that is not an aggregate. When an intersection does happen.
  // You should store the non-aggregate primitive in the intersection data
  // and not the BVH aggregate itself.

    bool hit = false;
//  for (size_t p = 0; p < primitives.size(); ++p) {
//    if(primitives[p]->intersect(ray, i)) hit = true;
//  }
    
    bool hitlbbox = false, hitrbbox = false;
    double t0, t1;
    if (root->bb.intersect(ray, t0, t1)) {
        stack<ClosestHit*> bvhstack_front;
        stack<ClosestHit*> bvhstack_back;
        ClosestHit* curBVH = new ClosestHit(root, t0);
        bvhstack_front.push(curBVH);
        ClosestHit* curHit = new ClosestHit(root, t1);
        while (!bvhstack_front.empty() || !bvhstack_back.empty()) {
            if (!bvhstack_front.empty()) {
                curBVH = bvhstack_front.top();
                bvhstack_front.pop();
            }else {
                curBVH = bvhstack_back.top();
                bvhstack_back.pop();
                if (curBVH->min_t > curHit->min_t) {
                    continue;
                }
            }
            
            
            if (curBVH->closestbvh->isLeaf()) {
                for (size_t j = 0; j < curBVH->closestbvh->range; j++) {
                    Intersection itsct;
                    if (primitives[curBVH->closestbvh->start + j]->intersect(ray, &itsct)) {
                        if (itsct.t < curHit->min_t) {
                            hit = true;
                            curHit->closestbvh = curBVH->closestbvh;
                            curHit->min_t = itsct.t;
                            i->t = itsct.t;
                            i->primitive = itsct.primitive;
                            i->bsdf = itsct.bsdf;
                            i->n = itsct.n;
                        }
                    }
                }
            }else {
            
                double tl0, tl1, tr0, tr1;
                tl0 = t0; tr0 = t0; tl1 = t1; tr1 = t1;
                hitlbbox = curBVH->closestbvh->l->bb.intersect(ray, tl0, tl1);
                hitrbbox = curBVH->closestbvh->r->bb.intersect(ray, tr0, tr1);
            
                if (hitlbbox && hitrbbox) {
                    if (tl0 < tr0) {
                        ClosestHit* leftHit = new ClosestHit(curBVH->closestbvh->l, tl0);
                        ClosestHit* rightHit = new ClosestHit(curBVH->closestbvh->r, tr0);
                        bvhstack_front.push(leftHit);
                        bvhstack_back.push(rightHit);
                        t0 = tl0;
                        t1 = tl1;
                    }else {
                        ClosestHit* leftHit = new ClosestHit(curBVH->closestbvh->l, tl0);
                        ClosestHit* rightHit = new ClosestHit(curBVH->closestbvh->r, tr0);
                        bvhstack_back.push(leftHit);
                        bvhstack_front.push(rightHit);
                        t0 = tr0;
                        t1 = tr1;
                    }
                }
                if (hitrbbox && !hitlbbox) {
                    ClosestHit* rightHit = new ClosestHit(curBVH->closestbvh->r, tr0);
                    bvhstack_front.push(rightHit);
                }
                if (hitlbbox && !hitrbbox) {
                    ClosestHit* leftHit = new ClosestHit(curBVH->closestbvh->l, tl0);
                    bvhstack_front.push(leftHit);
                }
                if (!hitlbbox && !hitrbbox && bvhstack_front.empty() && bvhstack_back.empty()) {
                    return hit;
                }
            }
        }
        if (curHit->closestbvh->range == 2 && curHit->closestbvh->start == 5) {
            
        }
        
        
    }else {
        return hit;
    }
    
    return hit;

}

}  // namespace StaticScene
}  // namespace CMU462
