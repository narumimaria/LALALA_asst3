#include "pathtracer.h"
#include "bsdf.h"
#include "ray.h"

#include <stack>
#include <random>
#include <algorithm>

#include "CMU462/CMU462.h"
#include "CMU462/vector3D.h"
#include "CMU462/matrix3x3.h"
#include "CMU462/lodepng.h"

#include "GL/glew.h"

#include "static_scene/sphere.h"
#include "static_scene/triangle.h"
#include "static_scene/light.h"

using namespace CMU462::StaticScene;

using std::min;
using std::max;

namespace CMU462 {

//#define ENABLE_RAY_LOGGING 1

PathTracer::PathTracer(size_t ns_aa,
                       size_t max_ray_depth, size_t ns_area_light,
                       size_t ns_diff, size_t ns_glsy, size_t ns_refr,
                       size_t num_threads, HDRImageBuffer* envmap) {
  state = INIT,
  this->ns_aa = ns_aa;
  this->max_ray_depth = max_ray_depth;
  this->ns_area_light = ns_area_light;
  this->ns_diff = ns_diff;
  this->ns_glsy = ns_diff;
  this->ns_refr = ns_refr;

  if (envmap) {
    this->envLight = new EnvironmentLight(envmap);
  } else {
    this->envLight = NULL;
  }

  bvh = NULL;
  scene = NULL;
  camera = NULL;

  gridSampler = new UniformGridSampler2D();
  hemisphereSampler = new UniformHemisphereSampler3D();

  show_rays = true;

  imageTileSize = 32;
  numWorkerThreads = num_threads;
  workerThreads.resize(numWorkerThreads);

  tm_gamma = 2.2f;
  tm_level = 1.0f;
  tm_key = 0.18;
  tm_wht = 5.0f;

}

PathTracer::~PathTracer() {

  delete bvh;
  delete gridSampler;
  delete hemisphereSampler;

}

void PathTracer::set_scene(Scene *scene) {

  if (state != INIT) {
    return;
  }

  if (this->scene != nullptr) {
    delete scene;
    delete bvh;
    selectionHistory.pop();
  }

  if (this->envLight != nullptr) {
    scene->lights.push_back(this->envLight);
  }

  this->scene = scene;
  build_accel();

  if (has_valid_configuration()) {
    state = READY;
  }
}

void PathTracer::set_camera(Camera *camera) {

  if (state != INIT) {
    return;
  }

  this->camera = camera;
  if (has_valid_configuration()) {
    state = READY;
  }

}

void PathTracer::set_frame_size(size_t width, size_t height) {
  if (state != INIT && state != READY) {
    stop();
  }
  sampleBuffer.resize(width, height);
  frameBuffer.resize(width, height);
  if (has_valid_configuration()) {
    state = READY;
  }
}

bool PathTracer::has_valid_configuration() {
  return scene && camera && gridSampler && hemisphereSampler &&
         (!sampleBuffer.is_empty());
}

void PathTracer::update_screen() {
  switch (state) {
    case INIT:
    case READY:
      break;
    case VISUALIZE:
      visualize_accel();
      break;
    case RENDERING:
      glDrawPixels(frameBuffer.w, frameBuffer.h, GL_RGBA,
                   GL_UNSIGNED_BYTE, &frameBuffer.data[0]);
      break;
    case DONE:
        //sampleBuffer.tonemap(frameBuffer, tm_gamma, tm_level, tm_key, tm_wht);
      glDrawPixels(frameBuffer.w, frameBuffer.h, GL_RGBA,
                   GL_UNSIGNED_BYTE, &frameBuffer.data[0]);
      break;
  }
}

void PathTracer::stop() {
  switch (state) {
    case INIT:
    case READY:
      break;
    case VISUALIZE:
      while (selectionHistory.size() > 1) {
        selectionHistory.pop();
      }
      state = READY;
      break;
    case RENDERING:
      continueRaytracing = false;
    case DONE:
      for (int i=0; i<numWorkerThreads; i++) {
            workerThreads[i]->join();
            delete workerThreads[i];
        }
      state = READY;
      break;
  }
}

void PathTracer::clear() {
  if (state != READY) return;
  delete bvh;
  bvh = NULL;
  scene = NULL;
  camera = NULL;
  selectionHistory.pop();
  sampleBuffer.resize(0, 0);
  frameBuffer.resize(0, 0);
  state = INIT;
}

void PathTracer::start_visualizing() {
  if (state != READY) {
    return;
  }
  state = VISUALIZE;
}

void PathTracer::start_raytracing() {
  if (state != READY) return;

  rayLog.clear();
  workQueue.clear();

  state = RENDERING;
  continueRaytracing = true;
  workerDoneCount = 0;

  sampleBuffer.clear();
  frameBuffer.clear();
  num_tiles_w = sampleBuffer.w / imageTileSize + 1;
  num_tiles_h = sampleBuffer.h / imageTileSize + 1;
  tile_samples.resize(num_tiles_w * num_tiles_h);
  memset(&tile_samples[0], 0, num_tiles_w * num_tiles_h * sizeof(int));

  // populate the tile work queue
  for (size_t y = 0; y < sampleBuffer.h; y += imageTileSize) {
      for (size_t x = 0; x < sampleBuffer.w; x += imageTileSize) {
          workQueue.put_work(WorkItem(x, y, imageTileSize, imageTileSize));
      }
  }

  // launch threads
  fprintf(stdout, "[PathTracer] Rendering... "); fflush(stdout);
  for (int i=0; i<numWorkerThreads; i++) {
      workerThreads[i] = new std::thread(&PathTracer::worker_thread, this);
  }
}


void PathTracer::build_accel() {

  // collect primitives //
  fprintf(stdout, "[PathTracer] Collecting primitives... "); fflush(stdout);
  timer.start();
  vector<Primitive *> primitives;
  for (SceneObject *obj : scene->objects) {
    const vector<Primitive *> &obj_prims = obj->get_primitives();
    primitives.reserve(primitives.size() + obj_prims.size());
    primitives.insert(primitives.end(), obj_prims.begin(), obj_prims.end());
  }
  timer.stop();
  fprintf(stdout, "Done! (%.4f sec)\n", timer.duration());

  // build BVH //
  fprintf(stdout, "[PathTracer] Building BVH... "); fflush(stdout);
  timer.start();
  bvh = new BVHAccel(primitives);
  timer.stop();
  fprintf(stdout, "Done! (%.4f sec)\n", timer.duration());

  // initial visualization //
  selectionHistory.push(bvh->get_root());
}

void PathTracer::log_ray_miss(const Ray& r) {
    rayLog.push_back(LoggedRay(r, -1.0));
}

void PathTracer::log_ray_hit(const Ray& r, double hit_t) {
    rayLog.push_back(LoggedRay(r, hit_t));
}

void PathTracer::visualize_accel() const {

  glPushAttrib(GL_ENABLE_BIT);
  glDisable(GL_LIGHTING);
  glLineWidth(1);
  glEnable(GL_DEPTH_TEST);

  // hardcoded color settings
  Color cnode = Color(.5, .5, .5, .25);
  Color cnode_hl = Color(1., .25, .0, .6);
  Color cnode_hl_child = Color(1., 1., 1., .6);

  Color cprim_hl_left = Color(.6, .6, 1., 1);
  Color cprim_hl_right = Color(.8, .8, 1., 1);
  Color cprim_hl_edges = Color(0., 0., 0., 0.5);

  BVHNode *selected = selectionHistory.top();

  // render solid geometry (with depth offset)
  glPolygonOffset(1.0, 1.0);
  glEnable(GL_POLYGON_OFFSET_FILL);

  if (selected->isLeaf()) {
    for (size_t i = 0; i < selected->range; ++i) {
       bvh->primitives[selected->start + i]->draw(cprim_hl_left);
    }
  } else {
      if (selected->l) {
          BVHNode* child = selected->l;
          for (size_t i = 0; i < child->range; ++i) {
              bvh->primitives[child->start + i]->draw(cprim_hl_left);
          }
      }
      if (selected->r) {
          BVHNode* child = selected->r;
          for (size_t i = 0; i < child->range; ++i) {
              bvh->primitives[child->start + i]->draw(cprim_hl_right);
          }
      }
  }

  glDisable(GL_POLYGON_OFFSET_FILL);

  // draw geometry outline
  for (size_t i = 0; i < selected->range; ++i) {
      bvh->primitives[selected->start + i]->drawOutline(cprim_hl_edges);
  }

  // keep depth buffer check enabled so that mesh occluded bboxes, but
  // disable depth write so that bboxes don't occlude each other.
  glDepthMask(GL_FALSE);

  // create traversal stack
  stack<BVHNode *> tstack;

  // push initial traversal data
  tstack.push(bvh->get_root());

  // draw all BVH bboxes with non-highlighted color
  while (!tstack.empty()) {

    BVHNode *current = tstack.top();
    tstack.pop();

    current->bb.draw(cnode);
    if (current->l) tstack.push(current->l);
    if (current->r) tstack.push(current->r);
  }

  // draw selected node bbox and primitives
  if (selected->l) selected->l->bb.draw(cnode_hl_child);
  if (selected->r) selected->r->bb.draw(cnode_hl_child);

  glLineWidth(3.f);
  selected->bb.draw(cnode_hl);

  // now perform visualization of the rays
  if (show_rays) {
      glLineWidth(1.f);
      glBegin(GL_LINES);

      for (size_t i=0; i<rayLog.size(); i+=500) {

          const static double VERY_LONG = 10e4;
          double ray_t = VERY_LONG;

          // color rays that are hits yellow
          // and rays this miss all geometry red
          if (rayLog[i].hit_t >= 0.0) {
              ray_t = rayLog[i].hit_t;
              glColor4f(1.f, 1.f, 0.f, 0.1f);
          } else {
              glColor4f(1.f, 0.f, 0.f, 0.1f);
          }

          Vector3D end = rayLog[i].o + ray_t * rayLog[i].d;

          glVertex3f(rayLog[i].o[0], rayLog[i].o[1], rayLog[i].o[2]);
          glVertex3f(end[0], end[1], end[2]);
      }
      glEnd();
  }

  glDepthMask(GL_TRUE);
  glPopAttrib();
}

void PathTracer::key_press(int key) {

  BVHNode *current = selectionHistory.top();
  switch (key) {
  case ']':
      ns_aa *=2;
      printf("Samples per pixel changed to %lu\n", ns_aa);
      //tm_key = clamp(tm_key + 0.02f, 0.0f, 1.0f);
      break;
  case '[':
      //tm_key = clamp(tm_key - 0.02f, 0.0f, 1.0f);
      ns_aa /=2;
      if (ns_aa < 1) ns_aa = 1;
      printf("Samples per pixel changed to %lu\n", ns_aa);
      break;
  case KEYBOARD_UP:
      if (current != bvh->get_root()) {
          selectionHistory.pop();
      }
      break;
  case KEYBOARD_LEFT:
      if (current->l) {
          selectionHistory.push(current->l);
      }
      break;
  case KEYBOARD_RIGHT:
      if (current->l) {
          selectionHistory.push(current->r);
      }
      break;
  case 'a':
  case 'A':
      show_rays = !show_rays;
  default:
      return;
  }
}

Spectrum PathTracer::trace_ray(const Ray &r) {
    
//    Vector3D test_d1 = Vector3D(-1, 0, -1);
//    Vector3D test_d2 = Vector3D(0, -0.2, 0.8);
//    test_d1.normalize();
//    test_d2.normalize();
//    Spectrum s1 = envLight->sample_dir(test_d1);
//    Spectrum s2 = envLight->sample_dir(test_d2);
//    GlassBSDF g_bsdf{ {1,1,1}, {1,1,1}, 1, 2};
//    Vector3D test_out{0.01,0.01,0.99};
//    test_out.normalize();
//    Vector3D test_in{0,0,1};
//    test_in.normalize();
//    float pdf{};
//    srand(19920810);
//    Spectrum result = g_bsdf.f(test_out, test_in);
//    result = g_bsdf.sample_f(test_out, &test_in, &pdf);
//    return envLight->sample_dir(r.d);
    
  Intersection isect;

  if (!bvh->intersect(r, &isect)) {

    // log ray miss
    #ifdef ENABLE_RAY_LOGGING
    log_ray_miss(r);
    #endif

    // TODO:
    // If you have an environment map, return the Spectrum this ray
    // samples from the environment map. If you don't return black.
//    return Spectrum(0.99,0.96,0.6);
//      return envLight->sample_L(p,&wi,&dtl,&pdf);
      if (envLight == nullptr) {
          return Spectrum(0,0,0);
      }else{
          return envLight->sample_dir(r.d);
      }
  }

  // log ray hit
  #ifdef ENABLE_RAY_LOGGING
  log_ray_hit(r, isect.t);
  #endif

    Spectrum L_e = isect.bsdf->get_emission(); // Le
    if(L_e.r > 1) {
        L_e.r = 1;
    }
    if(L_e.g > 1) L_e.g = 1;
    if(L_e.b > 1) L_e.b = 1;
    
    Spectrum L_out = L_e;
    Spectrum L_ind = Spectrum(0,0,0);
  // TODO :
  // Instead of initializing this value to a constant color, use the direct,
  // indirect lighting components calculated in the code below. The starter
  // code overwrites L_out by (.5,.5,.5) so that you can test your geometry
  // queries before you implement path tracing.
//  L_out = Spectrum(0, 0, 0);

  Vector3D hit_p = r.o + r.d * isect.t;
  

  // make a coordinate system for a hit point
  // with N aligned with the Z direction.
  Matrix3x3 o2w_1;
  make_coord_space(o2w_1, isect.n);
  Matrix3x3 w2o = o2w_1.T();

  // w_out points towards the source of the ray (e.g.,
  // toward the camera if this is a primary ray)
  Vector3D w_out = w2o * (r.o - hit_p);
  w_out.normalize();
    if (dot(isect.n, r.d) > 0) {
        isect.n *= -1;
    }
    Vector3D hit_n = isect.n;
    Matrix3x3 o2w_2;
    make_coord_space(o2w_2, isect.n);
    w2o = o2w_2.T();

  // TODO:
  // Extend the below code to compute the direct lighting for all the lights
  // in the scene, instead of just the dummy light we provided in part 1.

//  InfiniteHemisphereLight light(Spectrum(.5f, .5f, .5f));
//  DirectionalLight light(Spectrum(.5f, .5f, .5f), Vector3D(1.0, -1.0, 0.0));
    size_t num_lights = scene->lights.size();
    for (size_t j = 0; j < num_lights; j++) {
//        if(!dynamic_cast<EnvironmentLight*>(scene->lights[j]))
//            continue;
        Vector3D dir_to_light;
        float dist_to_light;
        float pdf;

        // no need to take multiple samples from a directional source
//        int num_light_samples = light.is_delta_light() ? 1 : ns_area_light;
        int num_light_samples = scene->lights[j]->is_delta_light() ? 1 : ns_area_light;

        // integrate light over the hemisphere about the normal
        double scale = 1.0 / num_light_samples;
        for (int i=0; i<num_light_samples; i++) {

            // returns a vector 'dir_to_light' that is a direction from
            // point hit_p to the point on the light source.  It also returns
            // the distance from point x to this point on the light source.
            // (pdf is the probability of randomly selecting the random
            // sample point on the light source -- more on this in part 2)
            Spectrum light_L = scene->lights[j]->sample_L(hit_p, &dir_to_light, &dist_to_light, &pdf);
//            if (light_L.r > 1) {
//                
//            }
            // convert direction into coordinate space of the surface, where
            // the surface normal is [0 0 1]
            Vector3D w_in = w2o * dir_to_light;

            // note that computing dot(n,w_in) is simple
            // in surface coordinates since the normal is [0 0 1]
//            double cos_theta = std::max(0.0, w_in[2]);

            // evaluate surface bsdf
            Spectrum f = isect.bsdf->f(w_out, w_in);

            // TODO:
            // Construct a shadow ray and compute whether the intersected surface is
            // in shadow and accumulate reflected radiance
            Intersection itsct_occlude;
            Vector3D origin = hit_p + EPS_D * dir_to_light;
            Ray r_occu = Ray(origin, dir_to_light);
            r_occu.max_t = dist_to_light;
            if (bvh->intersect(r_occu, &itsct_occlude)) {
                if (itsct_occlude.t > 0) {
              
                }else {
                    L_out += light_L * f * scale * (1 / pdf);
                }
            }else {
                L_out += light_L * f * scale * (1 / pdf);
            }
      
        }
        
    }

    // TODO:
    // Compute an indirect lighting estimate using pathtracing with Monte Carlo.
    // Note that Ray objects have a depth field now; you should use this to avoid
    // traveling down one path forever.
    if (r.depth == max_ray_depth) {
        
        
    }else {
        
        // randomly select a new ray direction (it may be
        // reflection or transmittence ray depending on
        // surface type
        Vector3D w_i;
        float pdf_i;
        Spectrum f_ind = isect.bsdf->sample_f(w_out, &w_i, &pdf_i);
        if (pdf_i < 0) {
            w_i = o2w_1 * w_i;
            pdf_i = -pdf_i;
            f_ind = -1 * f_ind;
        }else {
            w_i = o2w_2 * w_i;
        }
        
        w_i.normalize();
        
        
        double terminatep;
        float brightness = f_ind.illum();
        terminatep = 1 - brightness;
        terminatep = 0;
        if ((double)(std::rand()) / RAND_MAX < terminatep) {

        }else
        {
            Vector3D origin = hit_p + EPS_D * w_i;
            Spectrum light_ind = trace_ray(Ray(origin, w_i, INF_F, r.depth+1));
            L_ind += light_ind * f_ind * (1.0 / (pdf_i * (1 - terminatep)));
            L_out += L_ind;
        }
        
    }

  return L_out;
}

Spectrum PathTracer::raytrace_pixel(size_t x, size_t y) {

  // TODO:
  // Sample the pixel with coordinate (x,y) and return the result spectrum.
  // The sample rate is given by the number of camera rays per pixel.
//    if (x < 1300 || y < 380) {
//        return Spectrum(1,0,0);
//    }

  int num_samples = ns_aa;

  Vector2D p = Vector2D(0.5,0.5);
    
    Spectrum s;
    for (int i = 0; i < num_samples; i++) {
        Vector2D sample = gridSampler->get_sample();
        p.x = (x + sample.x) / sampleBuffer.w;
        p.y = (y + sample.y) / sampleBuffer.h;
//        if (x == 700 && y == 930) {
//            
//        }
        s += trace_ray(camera->generate_ray(p.x, p.y));
    }
    s = s * (1.0 / num_samples);
    
  return s;

}

void PathTracer::raytrace_tile(int tile_x, int tile_y,
                               int tile_w, int tile_h) {

  size_t w = sampleBuffer.w;
  size_t h = sampleBuffer.h;

  size_t tile_start_x = tile_x;
  size_t tile_start_y = tile_y;

  size_t tile_end_x = std::min(tile_start_x + tile_w, w);
  size_t tile_end_y = std::min(tile_start_y + tile_h, h);

  size_t tile_idx_x = tile_x / imageTileSize;
  size_t tile_idx_y = tile_y / imageTileSize;
  size_t num_samples_tile = tile_samples[tile_idx_x + tile_idx_y * num_tiles_w];

  for (size_t y = tile_start_y; y < tile_end_y; y++) {
    if (!continueRaytracing) return;
    for (size_t x = tile_start_x; x < tile_end_x; x++) {
        Spectrum s = raytrace_pixel(x, y);
        sampleBuffer.update_pixel(s, x, y);
    }
  }

  tile_samples[tile_idx_x + tile_idx_y * num_tiles_w] += 1;
  sampleBuffer.toColor(frameBuffer, tile_start_x, tile_start_y, tile_end_x, tile_end_y);
}

void PathTracer::worker_thread() {

  Timer timer;
  timer.start();

  WorkItem work;
  while (continueRaytracing && workQueue.try_get_work(&work)) {
    raytrace_tile(work.tile_x, work.tile_y, work.tile_w, work.tile_h);
  }

  workerDoneCount++;
  if (!continueRaytracing && workerDoneCount == numWorkerThreads) {
    timer.stop();
    fprintf(stdout, "Canceled!\n");
    state = READY;
  }

  if (continueRaytracing && workerDoneCount == numWorkerThreads) {
    timer.stop();
    fprintf(stdout, "Done! (%.4fs)\n", timer.duration());
    state = DONE;
  }
}

void PathTracer::increase_area_light_sample_count() {
  ns_area_light *= 2;
  fprintf(stdout, "[PathTracer] Area light sample count increased to %zu!\n", ns_area_light);
}

void PathTracer::decrease_area_light_sample_count() {
  if (ns_area_light > 1) ns_area_light /= 2;
  fprintf(stdout, "[PathTracer] Area light sample count decreased to %zu!\n", ns_area_light);
}

void PathTracer::save_image() {

  if (state != DONE) return;

  time_t rawtime;
  time (&rawtime);

  string filename = "Screen Shot ";
  filename += string(ctime(&rawtime));
  filename.erase(filename.end() - 1);
  filename += string(".png");

  uint32_t* frame = &frameBuffer.data[0];
  size_t w = frameBuffer.w;
  size_t h = frameBuffer.h;
  uint32_t* frame_out = new uint32_t[w * h];
  for(size_t i = 0; i < h; ++i) {
    memcpy(frame_out + i * w, frame + (h - i - 1) * w, 4 * w);
  }

  fprintf(stderr, "[PathTracer] Saving to file: %s... ", filename.c_str());
  lodepng::encode(filename, (unsigned char*) frame_out, w, h);
  fprintf(stderr, "Done!\n");
}

}  // namespace CMU462
