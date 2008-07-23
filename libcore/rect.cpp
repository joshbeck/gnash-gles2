// 
//   Copyright (C) 2005, 2006, 2007, 2008 Free Software Foundation, Inc.
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
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

#include "rect.h"
#include "log.h"
#include "stream.h"
#include "matrix.h"
#include "types.h" // for TWIPS_TO_PIXELS
#include "utility.h" // for flerp, clamp...

#include <sstream> // for ::print and ::toString

namespace gnash {

void    rect::read(SWFStream& in)
{
    in.align();
    in.ensureBits(5);
    int nbits = in.read_uint(5);
    in.ensureBits(nbits*4);
    
    _xMin = in.read_sint(nbits);
    _xMax = in.read_sint(nbits);
    _yMin = in.read_sint(nbits);
    _yMax = in.read_sint(nbits);

    // Check if this rect is valid.
    if (_xMax < _xMin || _yMax < _yMin)
    {
        // We set invalid rectangles to NULL, but we might instead
        // want to actually swap the values IFF the proprietary player
        // does so. TODO: check it out.
        IF_VERBOSE_MALFORMED_SWF(
        log_swferror("Invalid rectangle: "
            "xMin=%g xMax=%g yMin=%g yMax=%g", _xMin, _xMax, _yMin, _yMax);
        );
        set_null();
    } 
}

point
rect::get_point(int i) const
// Get one of the rect verts.
{
    assert( !is_null() );
    
    point p;
    switch(i)
    {
    case 0:
        p.x = _xMin; p.y = _yMin;
        break;
    case 1:
        p.x = _xMax; p.y = _yMin;
        break;
    case 2:
        p.x = _xMax; p.y = _yMax;
        break;
    case 3:
        p.x = _xMin; p.y = _yMax;
        break;
    default:
        assert(0);
        break;
    }
    return p;
}


void    rect::enclose_transformed_rect(const matrix& m, const rect& r)
// Set ourself to bound a rectangle that has been transformed by m.  
{   
    boost::int32_t  x1 = r.get_x_min();
    boost::int32_t  y1 = r.get_y_min();
    boost::int32_t  x2 = r.get_x_max();
    boost::int32_t  y2 = r.get_y_max();

    point  p0(x1, y1);
    point  p1(x2, y1);
    point  p2(x2, y2);
    point  p3(x1, y2);
    
    m.transform(p0);
    m.transform(p1);
    m.transform(p2);
    m.transform(p3);

    set_to_point(p0.x, p0.y);
    expand_to(p1.x, p1.y);
    expand_to(p2.x, p2.y);
    expand_to(p3.x, p3.y);
}

void  rect::expand_to_rect(const rect& r) 
// Expand ourself to enclose the given rect.
{    
    if( r.is_null() ) {
        return;
    }
    
    if( is_null() ) {
        *this = r;
    }else {
        _xMin = std::min(_xMin, r.get_x_min());
        _yMin = std::min(_yMin, r.get_y_min());
        _xMax = std::max(_xMax, r.get_x_max());
        _yMax = std::max(_yMax, r.get_y_max());
    }
}   

void    rect::expand_to_transformed_rect(const matrix& m, const rect& r)
// Expand ourself to a transformed rect.
{   
    if ( r.is_null() )
    {
         return;
    }

    boost::int32_t  x1 = r.get_x_min();
    boost::int32_t  y1 = r.get_y_min();
    boost::int32_t  x2 = r.get_x_max();
    boost::int32_t  y2 = r.get_y_max();

    point  p0(x1, y1);
    point  p1(x2, y1);
    point  p2(x2, y2);
    point  p3(x1, y2);
    
    m.transform(p0);
    m.transform(p1);
    m.transform(p2);
    m.transform(p3);

    if( is_null() ) {
        set_to_point(p0.x, p0.y);   
    }else {
        expand_to(p0.x, p0.y);
    }
    expand_to(p1.x, p1.y);
    expand_to(p2.x, p2.y);
    expand_to(p3.x, p3.y);
}

void    rect::set_lerp(const rect& a, const rect& b, float t)
// Set this to the lerp of a and b.
{
    assert( !a.is_null() );
    assert( !b.is_null() );
    
    using utility::flerp;   
    _xMin = (boost::int32_t)(flerp(a.get_x_min(), b.get_x_min(), t));
    _yMin = (boost::int32_t)(flerp(a.get_y_min(), b.get_y_min(), t));
    _xMax = (boost::int32_t)(flerp(a.get_x_max(), b.get_x_max(), t));
    _yMax = (boost::int32_t)(flerp(a.get_y_max(), b.get_y_max(), t));
}

void
rect::clamp(point& p) const
{
    assert( !is_null() );
    p.x = utility::clamp<boost::int32_t>(p.x, _xMin, _xMax);
    p.y = utility::clamp<boost::int32_t>(p.y, _yMin, _yMax);
}

std::string
rect::toString() const
{
    std::stringstream ss;
    ss << *this;
    return ss.str();
}

}   // end namespace gnash


// Local Variables:
// mode: C++
// c-basic-offset: 8 
// tab-width: 8
// indent-tabs-mode: t
// End: