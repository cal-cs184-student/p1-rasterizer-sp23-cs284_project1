#include "rasterizer.h"

using namespace std;

namespace CGL {

  RasterizerImp::RasterizerImp(PixelSampleMethod psm, LevelSampleMethod lsm,
    size_t width, size_t height,
    unsigned int sample_rate) {
    this->psm = psm;
    this->lsm = lsm;
    this->width = width;
    this->height = height;
    this->sample_rate = sample_rate;
    this->dilation = sqrt(sample_rate);

    sample_buffer.resize(width * height * sample_rate, Color::White);
  }

  // Used by rasterize_point and rasterize_line
  void RasterizerImp::fill_pixel(size_t x, size_t y, Color c) {
    // TODO: Task 2: You might need this function to fix points and lines (such as the black rectangle border in test4.svg)
    // NOTE: You are not required to implement proper supersampling for points and lines
    // It is sufficient to use the same color for all supersamples of a pixel for points and lines (not triangles)

    for (int yy = 0; yy < dilation; ++yy) {
      for (int xx = 0; xx < dilation; ++xx) {
        sample_buffer[(y * dilation + yy) * width * dilation + x * dilation + xx] = c;
      } 
    }
    // sample_buffer[y * dilation * width + x] = c;
  }

  // Rasterize a point: simple example to help you start familiarizing
  // yourself with the starter code.
  //
  void RasterizerImp::rasterize_point(float x, float y, Color color) {
    // fill in the nearest pixel
    int sx = (int)floor(x);
    int sy = (int)floor(y);

    // check bounds
    if (sx < 0 || sx >= width) return;
    if (sy < 0 || sy >= height) return;

    fill_pixel(sx, sy, color);
    return;
  }

  // Rasterize a line.
  void RasterizerImp::rasterize_line(float x0, float y0,
    float x1, float y1,
    Color color) {
    if (x0 > x1) {
      swap(x0, x1); swap(y0, y1);
    }

    float pt[] = { x0,y0 };
    float m = (y1 - y0) / (x1 - x0);
    float dpt[] = { 1,m };
    int steep = abs(m) > 1;
    if (steep) {
      dpt[0] = x1 == x0 ? 0 : 1 / abs(m);
      dpt[1] = x1 == x0 ? (y1 - y0) / abs(y1 - y0) : m / abs(m);
    }

    while (floor(pt[0]) <= floor(x1) && abs(pt[1] - y0) <= abs(y1 - y0)) {
      rasterize_point(pt[0], pt[1], color);
      pt[0] += dpt[0]; pt[1] += dpt[1];
    }
  }

  // Check if sample is inside frame boundaries.
  // void RasterizerImp::isInFrame(float x0, float y0,
  //   float x1, float y1,
  //   float x2, float y2,
  //   Color color) {
    

  // }

  float RasterizerImp::line_test(float x, float y,
    float x1, float y1,
    float x2, float y2) {
    return -(x - x1) * (y2 - y1) + (y - y1) * (x2 - x1);
  }

  // Rasterize a triangle.
  void RasterizerImp::rasterize_triangle(float x0, float y0,
    float x1, float y1,
    float x2, float y2,
    Color color) {
    // // TODO: Task 1: Implement basic triangle rasterization here, no supersampling
    // int x_min = min({x0, x1, x2}) - 1, x_max = max({x0, x1, x2}) + 1;
    // int y_min = min({y0, y1, y2}) - 1, y_max = max({y0, y1, y2}) + 1;
    // for (int y = y_min; y < y_max; y++) {
    //   for (int x = x_min; x < x_max; x++) {
    //     float x_sample = x + 0.5, y_sample = y + 0.5;

    //     //TODO: check boundary conditions for this
    //     if (x_sample < 0 || x_sample >= width || y_sample < 0 || y_sample >= height) continue;
        
    //     float test1 = line_test(x_sample, y_sample, x0, y0, x1, y1);
    //     float test2 = line_test(x_sample, y_sample, x1, y1, x2, y2);
    //     float test3 = line_test(x_sample, y_sample, x2, y2, x0, y0);
    //     if ((test1 >= 0 && test2 >= 0 && test3 >= 0) || (test1 <= 0 && test2 <= 0 && test3 <= 0)) {
    //       rasterize_point(x_sample, y_sample, color);
    //     }
    //   }
    // }

    // TODO: Task 2: Update to implement super-sampled rasterization
    x0 *= dilation, x1 *= dilation, x2 *= dilation;
    y0 *= dilation, y1 *= dilation, y2 *= dilation;
    int sample_frame_width = width * dilation, sample_frame_height = height * dilation;
    int x_min = min({x0, x1, x2}) - 1, x_max = max({x0, x1, x2}) + 1;
    int y_min = min({y0, y1, y2}) - 1, y_max = max({y0, y1, y2}) + 1;
    for (int y = y_min; y < y_max; y++) {
      for (int x = x_min; x < x_max; x++) {
        float x_sample = x + 0.5, y_sample = y + 0.5;

        //TODO: check boundary conditions for this
        if (x_sample < 0 || x_sample >= sample_frame_width || y_sample < 0 || y_sample >= sample_frame_height) continue;
        
        float test1 = line_test(x_sample, y_sample, x0, y0, x1, y1);
        float test2 = line_test(x_sample, y_sample, x1, y1, x2, y2);
        float test3 = line_test(x_sample, y_sample, x2, y2, x0, y0);
        if ((test1 >= 0 && test2 >= 0 && test3 >= 0) || (test1 <= 0 && test2 <= 0 && test3 <= 0)) {
          // fill_pixel(x_sample, y_sample, color);
          sample_buffer[(int)floor(y_sample) * sample_frame_width + (int)floor(x_sample)] = color;
        }
      }
    }


  }


  void RasterizerImp::rasterize_interpolated_color_triangle(float x0, float y0, Color c0,
    float x1, float y1, Color c1,
    float x2, float y2, Color c2)
  {
    // TODO: Task 4: Rasterize the triangle, calculating barycentric coordinates and using them to interpolate vertex colors across the triangle
    // Hint: You can reuse code from rasterize_triangle
    x0 *= dilation, x1 *= dilation, x2 *= dilation;
    y0 *= dilation, y1 *= dilation, y2 *= dilation;
    int sample_frame_width = width * dilation, sample_frame_height = height * dilation;
    int x_min = min({x0, x1, x2}) - 1, x_max = max({x0, x1, x2}) + 1;
    int y_min = min({y0, y1, y2}) - 1, y_max = max({y0, y1, y2}) + 1;
    for (int y = y_min; y < y_max; y++) {
      for (int x = x_min; x < x_max; x++) {
        float x_sample = x + 0.5, y_sample = y + 0.5;

        //TODO: check boundary conditions for this
        if (x_sample < 0 || x_sample >= sample_frame_width || y_sample < 0 || y_sample >= sample_frame_height) continue;
        
        float test1 = line_test(x_sample, y_sample, x0, y0, x1, y1);
        float test2 = line_test(x_sample, y_sample, x1, y1, x2, y2);
        float test3 = line_test(x_sample, y_sample, x2, y2, x0, y0);
        if ((test1 >= 0 && test2 >= 0 && test3 >= 0) || (test1 <= 0 && test2 <= 0 && test3 <= 0)) {
          float alpha = -(x_sample - x1) * (y2 - y1) + (y_sample - y1) * (x2 - x1) / (-(x0 - x1) * (y2 - y1) + (y0 - y1) * (x2 - x1));
          float beta = -(x_sample - x2) * (y0 - y2) + (y_sample - y2) * (x0 - x2) / (-(x1 - x2) * (y0 - y2) + (y1 - y2) * (x0 - x2));
          float gamma = 1 - alpha - beta;

          sample_buffer[(int)floor(y_sample) * sample_frame_width + (int)floor(x_sample)] = alpha * c0 + beta * c1 + gamma * c2;
          // sample_buffer[(int)floor(y_sample) * sample_frame_width + (int)floor(x_sample)] = c2;
        }
      }
    }



  }


  void RasterizerImp::rasterize_textured_triangle(float x0, float y0, float u0, float v0,
    float x1, float y1, float u1, float v1,
    float x2, float y2, float u2, float v2,
    Texture& tex)
  {
    // TODO: Task 5: Fill in the SampleParams struct and pass it to the tex.sample function.
    // TODO: Task 6: Set the correct barycentric differentials in the SampleParams struct.
    // Hint: You can reuse code from rasterize_triangle/rasterize_interpolated_color_triangle




  }

  void RasterizerImp::set_sample_rate(unsigned int rate) {
    // TODO: Task 2: You may want to update this function for supersampling support

    this->sample_rate = rate;
    this->dilation = sqrt(rate);


    this->sample_buffer.resize(width * height * rate, Color::White);
  }


  void RasterizerImp::set_framebuffer_target(unsigned char* rgb_framebuffer,
    size_t width, size_t height)
  {
    // TODO: Task 2: You may want to update this function for supersampling support

    this->width = width;
    this->height = height;
    this->rgb_framebuffer_target = rgb_framebuffer;


    this->sample_buffer.resize(width * height * this->sample_rate, Color::White);
  }


  void RasterizerImp::clear_buffers() {
    std::fill(rgb_framebuffer_target, rgb_framebuffer_target + 3 * width * height, 255);
    std::fill(sample_buffer.begin(), sample_buffer.end(), Color::White);
  }


  // This function is called at the end of rasterizing all elements of the
  // SVG file.  If you use a supersample buffer to rasterize SVG elements
  // for antialising, you could use this call to fill the target framebuffer
  // pixels from the supersample buffer data.
  //
  void RasterizerImp::resolve_to_framebuffer() {
    // TODO: Task 2: You will likely want to update this function for supersampling support

    for (int x = 0; x < width; ++x) {
      for (int y = 0; y < height; ++y) {
        Color col;
        for (int xx = 0; xx < dilation; ++xx) {
          for (int yy = 0; yy < dilation; ++yy) {
            int new_cord = (y * dilation + yy) * width * dilation + x * dilation + xx;
            col += sample_buffer[new_cord];
          }
        }

        for (int k = 0; k < 3; ++k) {
          this->rgb_framebuffer_target[3 * (y * width + x) + k] = (&col.r)[k] * 255 / this->sample_rate;
        }
      }
    }

  }

  Rasterizer::~Rasterizer() { }


}// CGL
