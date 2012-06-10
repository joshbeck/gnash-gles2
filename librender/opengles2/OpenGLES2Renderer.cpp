// 
//   Copyright (C) 2005, 2006, 2007, 2008, 2009, 2010, 2011, 2012
//   Free Software Foundation, Inc.
// 
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 3 of the License, or
// (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

#ifdef HAVE_CONFIG_H
#include "gnashconfig.h"
#endif

#include <sys/time.h>
#include <cstring>
#include <cmath>
#include <iostream>
#include <boost/utility.hpp>
#include <boost/bind.hpp>

#include "log.h"
#include "RGBA.h"
#include "GnashImage.h"
#include "GnashNumeric.h"
#include "FillStyle.h"
#include "LineStyle.h"
#include "Transform.h"
#include "log.h"
#include "utility.h"
#include "Range2d.h"
#include "SWFCxForm.h"
#include "opengles2/OpenGLES2Renderer.h"
//#include "openvg/OpenVGBitmap.h"
//#include "openvg/OpenVGStyle.h"
#include "SWFMatrix.h"
#include "swf/ShapeRecord.h"
#include "CachedBitmap.h"

#include <GLES2/gl2.h>
#include <GLES2/gl2ext.h>
#define GNASH_IMAGE_QUALITY     GL_IMAGE_QUALITY_FASTER
#define GNASH_RENDER_QUALITY    GL_RENDERING_QUALITY_FASTER

#define  MAX_POINTS (4096)

/// \file Renderer_gles2.cpp
/// \brief The OpenVG renderer and related code.
///

static const int TwipsPerInch = 1440;

namespace gnash {

typedef std::vector<Path> PathVec;
typedef std::vector<geometry::Range2d<int> > ClipBounds;

namespace renderer {

namespace gles2 {

/// Transforms the current OpenVG SWFMatrix using the given SWFMatrix.
/// When it goes out of scope, the SWFMatrix will be reset to what it
/// was before the new SWFMatrix was applied.
class eglScopeMatrix : public boost::noncopyable
{
public:
    eglScopeMatrix(const SWFMatrix& m)
        {
            // GNASH_REPORT_FUNCTION;            

                log_unimpl(_("eglScopeMatrix"));
        }
  
    ~eglScopeMatrix()
        {
                log_unimpl(_("eglScopeMatrix"));
        }
private:

};

/// @note
/// A VGpath is constructed from a series of appended path
/// segments. When drawing shapes from flash, we start each path by
/// moving to a known location. Then segments are appended, and then
/// the path is closed. This is also used for fills.

#define MAX_SEG  (256)
/*
/// Start a VGPath by moving to a specified location
///
/// @param path The VGPath to start
/// @returns nothing
inline void
startpath(VGPath path, const int x, const int y)
{
    log_unimpl(_("startpath"));
}

/// Close the VGPath started by startpath()
///
/// @param path The VGPath to close
/// @returns nothing
inline void
closepath(VGPath path)
{
    log_unimpl(_("closepath"));
}

/// Add a series of edges to the existing path created by startpath()
///
/// @param path The VGPath to append segments to
/// @param edges The segments to append to the path
/// @param anchor_x The X coordinate to start from
/// @param anchor_y The Y coordinate to start from
/// @returns nothing
inline void
preparepath(VGPath path, const std::vector<Edge>& edges,
                        const float& anchor_x, const float& anchor_y)
{
    log_unimpl(_("preparepath"));
}
*/
template<typename C, typename T, typename R, typename A>
void 
for_each(C& container, R (T::*pmf)(const A&),const A& arg)
{
    std::for_each(container.begin(), container.end(),
                  boost::bind(pmf, _1, boost::ref(arg)));
}

Renderer_gles2::Renderer_gles2()
    : _display_width(0.0),
      _display_height(0.0),
      _drawing_mask(false),
/*#ifdef OPENVG_VERSION_1_1    
      _mask_layer(VG_INVALID_HANDLE),
#endif
      _fillpaint(VG_INVALID_HANDLE),
      _strokepaint(VG_INVALID_HANDLE),*/
      _aspect_ratio(0.75)       // 4:3 aspect ratio
{
    // GNASH_REPORT_FUNCTION;
}

Renderer_gles2::Renderer_gles2(renderer::GnashDevice::dtype_t /* dtype */)
    : _display_width(0.0),
      _display_height(0.0),
      _drawing_mask(false),
/*#ifdef OPENVG_VERSION_1_1    
      _mask_layer(VG_INVALID_HANDLE),
#endif
      _fillpaint(VG_INVALID_HANDLE),
      _strokepaint(VG_INVALID_HANDLE),*/
      _aspect_ratio(0.75)       // 4:3 aspect ratio
{
    // GNASH_REPORT_FUNCTION;

    set_scale(1.0f, 1.0f);

#if 0
    _fillpaint = vgCreatePaint();
    
    _strokepaint = vgCreatePaint();
    
    // this paint object is used for solid, gradient, and pattern fills.
    vgSetPaint (_fillpaint,   VG_FILL_PATH);

    // this pain object is used for paths
    vgSetPaint (_strokepaint, VG_STROKE_PATH);
#endif
}

void
Renderer_gles2::init(float x, float y)
{
    // GNASH_REPORT_FUNCTION;

    log_unimpl(_("init"));
}

Renderer_gles2::~Renderer_gles2()
{
    // GNASH_REPORT_FUNCTION;

    log_unimpl(_("~Renderer_gles2"));
}

// Given an image, returns a pointer to a CachedBitmap class
// that can later be passed to fill_styleX_bitmap(), to set a
// bitmap fill style. We only cache the GnashImage here, as a
// VGImage can't be created yet until the renderer is initialized.
CachedBitmap *
Renderer_gles2::createCachedBitmap(std::auto_ptr<image::GnashImage> im)
{
    // GNASH_REPORT_FUNCTION;

    log_unimpl(_("createCachedBitmap"));
}

// Since we store drawing operations in display lists, we take special care
// to store video frame operations in their own display list, lest they be
// anti-aliased with the rest of the drawing. Since display lists cannot be
// concatenated this means we'll add up with several display lists for normal
// drawing operations.
void
Renderer_gles2::drawVideoFrame(image::GnashImage* /* frame */, const SWFMatrix* /* m */,
                             const SWFRect* /* bounds */, bool /*smooth*/)
{
    log_unimpl(_("drawVideoFrame"));
}

void
Renderer_gles2::world_to_pixel(int& x, int& y, float world_x, float world_y) const
{
//    GNASH_REPORT_FUNCTION;

    // negative pixels seems ok here... we don't
    // clip to valid range, use world_to_pixel(rect&)
    // and Intersect() against valid range instead.
    point p(world_x, world_y);
    stage_matrix.transform(p);
    x = (int)p.x;
    y = (int)p.y;
}

geometry::Range2d<int>
Renderer_gles2::world_to_pixel(const SWFRect& wb) const
{
//    GNASH_REPORT_FUNCTION;

    using namespace gnash::geometry;
    
    if ( wb.is_null() ) return Range2d<int>(nullRange);
    if ( wb.is_world() ) return Range2d<int>(worldRange);
    
    int xmin, ymin, xmax, ymax;

    world_to_pixel(xmin, ymin, wb.get_x_min(), wb.get_y_min());
    world_to_pixel(xmax, ymax, wb.get_x_max(), wb.get_y_max());
    
    return Range2d<int>(xmin, ymin, xmax, ymax);
}

geometry::Range2d<int>
Renderer_gles2::world_to_pixel(const geometry::Range2d<float>& wb) const
{
    // GNASH_REPORT_FUNCTION;

    if (wb.isNull() || wb.isWorld()) return wb;
    
    int xmin, ymin, xmax, ymax;
    
    world_to_pixel(xmin, ymin, wb.getMinX(), wb.getMinY());
    world_to_pixel(xmax, ymax, wb.getMaxX(), wb.getMaxY());
    
    return geometry::Range2d<int>(xmin, ymin, xmax, ymax);
}

point
Renderer_gles2::pixel_to_world(int x, int y) const
{
    // GNASH_REPORT_FUNCTION;

    point p(x, y);
    SWFMatrix mat = stage_matrix;
    mat.invert().transform(p);
    return p;
}

/// Setup the renderer to display by setting the Matrix for scaling,
/// shearing, and transformations.
///
/// @param width - stage width
/// @param height - stage height
/// @param x0 - minimum frame size in X dimension in twips
/// @param x1 - maximum frame size in X dimension in twips
/// @param y0 - minimum frame size in Y dimension in twips
/// @param y1 - maximum frame size in Y dimension in twips
void
Renderer_gles2::begin_display(gnash::rgba const&, int width, int height,
                            float x0, float x1, float y0, float y1)
{
    // GNASH_REPORT_FUNCTION;

    log_unimpl(_("begin_display"));
}

void
Renderer_gles2::end_display()
{
    // GNASH_REPORT_FUNCTION;
}

/// Draw a line-strip directly, using a thin, solid line. 
//
/// Can be used to draw empty boxes and cursors.
void
Renderer_gles2::drawLine(const std::vector<point>& coords, const rgba& fill,
                       const SWFMatrix& mat)
{
    // GNASH_REPORT_FUNCTION;
    
    log_unimpl(_("drawLine"));
}

void
Renderer_gles2::draw_poly(const std::vector<point>& corners,
                        const rgba& fill, const rgba& /* outline */,
                        const SWFMatrix& mat, bool /* masked */)
{
    // GNASH_REPORT_FUNCTION;

    log_unimpl(_("draw_poly"));
}

void
Renderer_gles2::set_antialiased(bool /* enable */)
{
    log_unimpl(_("set_antialiased"));
}

void
Renderer_gles2::begin_submit_mask()
{
    // GNASH_REPORT_FUNCTION;

    PathVec mask;
    _masks.push_back(mask);
    _drawing_mask = true;
}

void
Renderer_gles2::end_submit_mask()
{
    // GNASH_REPORT_FUNCTION;

    // If masking is disabled, rhen we can't use it
    if (_drawing_mask == true) {
        _drawing_mask = false;    
        apply_mask();
    }
}

/// Apply the current mask; nesting is supported.
///
/// This method marks the stencil buffer by incrementing every stencil pixel
/// by one every time a solid from one of the current masks is drawn. When
/// all the mask solids are drawn, we change the stencil operation to permit
/// only drawing where all masks have drawn, in other words, where all masks
/// intersect, or in even other words, where the stencil pixel buffer equals
/// the number of masks.
void
Renderer_gles2::apply_mask()
{  
    // GNASH_REPORT_FUNCTION;

    log_unimpl(_("apply_mask"));
}

void
Renderer_gles2::disable_mask()
{
    // GNASH_REPORT_FUNCTION;
    
    log_unimpl(_("disable_mask"));
}

void
Renderer_gles2::add_paths(const PathVec& path_vec)
{
    // GNASH_REPORT_FUNCTION;

    SWFCxForm dummy_cx;
    
    FillStyle coloring = FillStyle(SolidFill(rgba(0, 255, 0, 255)));

    draw_submask(path_vec, SWFMatrix(), dummy_cx, coloring);
}

Path
Renderer_gles2::reverse_path(const Path& cur_path)
{
    // GNASH_REPORT_FUNCTION;

    const Edge& cur_end = cur_path.m_edges.back();    
    
    float prev_cx = cur_end.cp.x;
    float prev_cy = cur_end.cp.y;        
    
    Path newpath(cur_end.ap.x, cur_end.ap.y, cur_path.m_fill1, cur_path.m_fill0, cur_path.m_line, cur_path.m_new_shape);
    
    float prev_ax = cur_end.ap.x;
    float prev_ay = cur_end.ap.y; 
    
    for (std::vector<Edge>::const_reverse_iterator it = cur_path.m_edges.rbegin()+1, end = cur_path.m_edges.rend(); it != end; ++it) {
        const Edge& cur_edge = *it;
        
        if (prev_ax == prev_cx && prev_ay == prev_cy) {
            prev_cx = cur_edge.ap.x;
            prev_cy = cur_edge.ap.y;      
        }
        
        Edge newedge(prev_cx, prev_cy, cur_edge.ap.x, cur_edge.ap.y); 
        
        newpath.m_edges.push_back(newedge);
        
        prev_cx = cur_edge.cp.x;
        prev_cy = cur_edge.cp.y;
        prev_ax = cur_edge.ap.x;
        prev_ay = cur_edge.ap.y;
        
    }
    
    Edge newlastedge(prev_cx, prev_cy, cur_path.ap.x, cur_path.ap.y);    
    newpath.m_edges.push_back(newlastedge);
    
    return newpath;
}

const Path *
Renderer_gles2::find_connecting_path(const Path& to_connect,
                                   std::list<const Path*> path_refs)
{        
    // GNASH_REPORT_FUNCTION;

    float target_x = to_connect.m_edges.back().ap.x;
    float target_y = to_connect.m_edges.back().ap.y;
    
    if (target_x == to_connect.ap.x &&
        target_y == to_connect.ap.y) {
        return NULL;
    }
    
    for (std::list<const Path*>::const_iterator it = path_refs.begin(),
             end = path_refs.end(); it != end; ++it) {
        const Path* cur_path = *it;
        
        if (cur_path == &to_connect) {
            continue;
        }
        
        if (cur_path->ap.x == target_x && cur_path->ap.y == target_y) {
            if (cur_path->m_fill1 != to_connect.m_fill1) {
                continue;
            }
            return cur_path;
        }
    }    
    
    return NULL;  
}

PathVec
Renderer_gles2::normalize_paths(const PathVec &paths)
{
    // GNASH_REPORT_FUNCTION;

    PathVec normalized;
    
    for (PathVec::const_iterator it = paths.begin(), end = paths.end();
         it != end; ++it) {
        const Path& cur_path = *it;
        
        if (cur_path.m_edges.empty()) {
            continue;
            
        } else if (cur_path.m_fill0 && cur_path.m_fill1) {     
            
            // Two fill styles; duplicate and then reverse the left-filled one.
            normalized.push_back(cur_path);
            normalized.back().m_fill0 = 0; 
            
            Path newpath = reverse_path(cur_path);
            newpath.m_fill0 = 0;        
            
            normalized.push_back(newpath);       
            
        } else if (cur_path.m_fill0) {
            // Left fill style.
            Path newpath = reverse_path(cur_path);
            newpath.m_fill0 = 0;
            
            normalized.push_back(newpath);
        } else if (cur_path.m_fill1) {
            // Right fill style.
            normalized.push_back(cur_path);
        } else {
            // No fill styles; copy without modifying.
            normalized.push_back(cur_path);
        }
        
    }
    
    return normalized;
}


/// Analyzes a set of paths to detect real presence of fills and/or outlines
/// TODO: This should be something the character tells us and should be 
/// cached. 
void
Renderer_gles2::analyze_paths(const PathVec &paths, bool& have_shape,
                            bool& have_outline) 
{
    // GNASH_REPORT_FUNCTION;

    have_shape = false;
    have_outline = false;
    
    int pcount = paths.size();
    
    for (int pno= 0; pno<pcount; pno++) {
        
        const Path &the_path = paths[pno];
        
        // If a left or right fill is set, then this is an outline
        if ((the_path.m_fill0 > 0) || (the_path.m_fill1 > 0)) {
            have_shape = true;
            if (have_outline) return; // have both
        }
        
        // If a line is set, then it's a shape. A path can be both
        if (the_path.m_line > 0) {
            have_outline = true;
            if (have_shape) return; // have both
        }
    }    
}

bool
Renderer_gles2::apply_line_style(const LineStyle& style, const SWFCxForm& cx, 
                               const SWFMatrix& mat)
{
    // GNASH_REPORT_FUNCTION;

    log_unimpl(_("apply_line_style"));
}

typedef std::vector<const Path*> PathPtrVec;
  
void
Renderer_gles2::draw_outlines(const PathVec& path_vec, const SWFMatrix& mat,
                            const SWFCxForm& cx, const std::vector<LineStyle>& line_styles)
{
    // GNASH_REPORT_FUNCTION;

    log_unimpl(_("draw_outlines"));
}

std::list<PathPtrVec>
Renderer_gles2::get_contours(const PathPtrVec &paths)
{
    // GNASH_REPORT_FUNCTION;

    std::list<const Path*> path_refs;
    std::list<PathPtrVec> contours;
    
    for (PathPtrVec::const_iterator it = paths.begin(), end = paths.end();
         it != end; ++it) {
        const Path* cur_path = *it;
        path_refs.push_back(cur_path);
    }
    
    for (std::list<const Path*>::const_iterator it = path_refs.begin(), end = path_refs.end();
         it != end; ++it) {
        const Path* cur_path = *it;
        
        if (cur_path->m_edges.empty()) {
            continue;
        }
        
        if (!cur_path->m_fill0 && !cur_path->m_fill1) {
            continue;
        }
        
        PathPtrVec contour;
        
        contour.push_back(cur_path);
        
        const Path* connector = find_connecting_path(*cur_path, path_refs);
        
        while (connector) {       
            contour.push_back(connector);
            
            const Path* tmp = connector;
            connector = find_connecting_path(*connector, std::list<const Path*>(boost::next(it), end));
            
            // make sure we don't iterate over the connecting path in the for loop.
            path_refs.remove(tmp);
            
        } 
        
        contours.push_back(contour);   
    }
    
    return contours;
}

void
Renderer_gles2::draw_mask(const PathVec& path_vec)
{ 
    // GNASH_REPORT_FUNCTION;
   
    if (_drawing_mask == true) {
        for (PathVec::const_iterator it = path_vec.begin(), end = path_vec.end();
             it != end; ++it) {
            const Path& cur_path = *it;
            
            if (cur_path.m_fill0 || cur_path.m_fill1) {
                _masks.back().push_back(cur_path);
                _masks.back().back().m_line = 0;    
            }
        }  
    }
}

PathPtrVec
Renderer_gles2::paths_by_style(const PathVec& path_vec, unsigned int style)
{
    // GNASH_REPORT_FUNCTION;

    PathPtrVec paths;
    for (PathVec::const_iterator it = path_vec.begin(), end = path_vec.end();
         it != end; ++it) {
        const Path& cur_path = *it;
        
        if (cur_path.m_fill0 == style) {
            paths.push_back(&cur_path);
        }
        
        if (cur_path.m_fill1 == style) {
            paths.push_back(&cur_path);
        }
        
    }
    return paths;
}


std::vector<PathVec::const_iterator>
Renderer_gles2::find_subshapes(const PathVec& path_vec)
{
    // GNASH_REPORT_FUNCTION;

    std::vector<PathVec::const_iterator> subshapes;
    
    PathVec::const_iterator it = path_vec.begin(),
        end = path_vec.end();
    
    subshapes.push_back(it);
    ++it;
    
    for (;it != end; ++it) {
        const Path& cur_path = *it;
        
        if (cur_path.m_new_shape) {
            subshapes.push_back(it); 
        } 
    }
    
    if (subshapes.back() != end) {
        subshapes.push_back(end);
    }
    
    return subshapes;
}

/// Takes a path and translates it using the given SWFMatrix.
void
Renderer_gles2::apply_matrix_to_paths(std::vector<Path>& paths, const SWFMatrix& mat)
{  
    // GNASH_REPORT_FUNCTION;

    std::for_each(paths.begin(), paths.end(),
                  boost::bind(&Path::transform, _1, boost::ref(mat)));
}  

void
Renderer_gles2::draw_subshape(const PathVec& path_vec,
                            const SWFMatrix& mat,
                            const SWFCxForm& cx,
                            const std::vector<FillStyle>& fill_styles,
                            const std::vector<LineStyle>& line_styles)
{
    // GNASH_REPORT_FUNCTION;

    log_unimpl(_("draw_subshape"));
}

void
Renderer_gles2::draw_submask(const PathVec& path_vec,
                           const SWFMatrix& /* mat */,
                           const SWFCxForm& /* cx */,
                            const FillStyle& /* f_style */)
{
    log_unimpl(_("draw_submask"));
}

// Drawing procedure:
// 1. Separate paths by subshape.
// 2. Separate subshapes by fill style.
// 3. For every subshape/fill style combo:
//  a. Separate contours: find closed shapes by connecting ends.
//  b. Apply fill style.
//  c. Feed the contours in the tesselator. (Render.)
//  d. Draw outlines for every path in the subshape with a line style.
void
Renderer_gles2::drawShape(gnash::SWF::ShapeRecord const &shape, 
                        gnash::Transform const& xform)
{
    // GNASH_REPORT_FUNCTION;

    const PathVec& path_vec = shape.paths();
    
    if (!path_vec.size()) {
        // No paths. Nothing to draw...
        return;
    }
    if (_drawing_mask) {
        PathVec scaled_path_vec = path_vec;
        apply_matrix_to_paths(scaled_path_vec, xform.matrix);
        draw_mask(scaled_path_vec); 
        return;
    }    
    
    bool have_shape, have_outline;
    
    analyze_paths(path_vec, have_shape, have_outline);
    
    if (!have_shape && !have_outline) {
        return;
    }    
    
    eglScopeMatrix scope_mat(xform.matrix);
    
    std::vector<PathVec::const_iterator> subshapes = find_subshapes(path_vec);
    
    const std::vector<FillStyle>& fill_styles = shape.fillStyles();
    const std::vector<LineStyle>& line_styles = shape.lineStyles();
    
    for (size_t i = 0; i < subshapes.size()-1; ++i) {
        PathVec subshape_paths;
        
        if (subshapes[i] != subshapes[i+1]) {
            subshape_paths = PathVec(subshapes[i], subshapes[i+1]);
        } else {
            subshape_paths.push_back(*subshapes[i]);
        }
        
        draw_subshape(subshape_paths, xform.matrix, xform.colorTransform,
                      fill_styles, line_styles);
    }
}

void
Renderer_gles2::drawGlyph(const SWF::ShapeRecord& rec, const rgba& c,
                        const SWFMatrix& mat)
{
    // GNASH_REPORT_FUNCTION;

    if (_drawing_mask) {
        abort();
    }
    
    if (rec.getBounds().is_null()) {
        return;
    }
    
    SWFCxForm dummy_cx;
    std::vector<FillStyle> glyph_fs;
    
    FillStyle coloring = FillStyle(SolidFill(c));

    glyph_fs.push_back(coloring);

    std::vector<LineStyle> dummy_ls;
    
    eglScopeMatrix scope_mat(mat);
    
    draw_subshape(rec.paths(), mat, dummy_cx, glyph_fs, dummy_ls);
    
}

void
Renderer_gles2::set_scale(float xscale, float yscale)
{
    // GNASH_REPORT_FUNCTION;

    _xscale = xscale;
    _yscale = yscale;
    stage_matrix.set_identity();
    stage_matrix.set_scale(xscale/20.0f, yscale/20.0f);
}

void
Renderer_gles2::set_invalidated_regions(const InvalidatedRanges& /* ranges */)
{
    // GNASH_REPORT_FUNCTION;

    // do nothing obviously. This method is required by the base class though,
    // so something has to be here.
}
  
DSOEXPORT Renderer *
create_handler(const char */* pixelformat */)
{
    // GNASH_REPORT_FUNCTION;

    Renderer_gles2 *renderer = new Renderer_gles2;
    return renderer;
}


// These methods are only for debugging and development
void
Renderer_gles2::printGLES2Params()
{
    log_unimpl(_("printGLES2Params"));
}

/*void
Renderer_gles2::printGLES2Path(GLPath path)
{
    log_unimpl(_("printGLES2path"));
}*/

// Print an OpenGLES2 matrix, which is 3 x 3. Elements 2 and 5 are
// ignored, as they are the w0 and w1 paramaters.
// It looks like this: { sx, shy, w0, shx, sy, w1, tx, ty, w2 }
void
Renderer_gles2::printGLES2Matrix(GLfloat *mat)
{
    std::cerr << "sx, shx, tx: " << mat[0] << ", " << mat[1]<< ", "
              << std::fixed << mat[3] << std::endl;
    std::cerr << "sy, shy, ty: " << mat[4]<< ", " << mat[6] << ", "
              << std::scientific << mat[7] << std::endl;
}

void
Renderer_gles2::printGLES2Matrix(const SWFMatrix &mat)
{
    std::cerr << "a, shx, tx: " << mat.a() << ", " << mat.b() << ", " << mat.tx() << std::endl;
    std::cerr << "sy, shy, ty: " << mat.c() << ", " << mat.d() << ", " << mat.ty() << std::endl;
}

Renderer *
Renderer_gles2::startInternalRender(gnash::image::GnashImage&)
{
    // GNASH_REPORT_FUNCTION;

    return 0;
}

void
Renderer_gles2::endInternalRender()
{
    // GNASH_REPORT_FUNCTION;    
}

void
Renderer_gles2::drawVideoFrame(gnash::image::GnashImage*, gnash::Transform const&, gnash::SWFRect const*, bool)
{
    GNASH_REPORT_FUNCTION;
}

unsigned int
Renderer_gles2::getBitsPerPixel()
{
    return 0;
}

const char *
Renderer_gles2::getErrorString(GLenum error)
{
    log_unimpl(_("getErrorString"));
}

} // namespace gnash::renderer::gles2
} // namespace gnash::renderer
} // namespace gnash

// local Variables:
// mode: C++
// indent-tabs-mode: nil
// End:
