// 
//   Copyright (C) 2007 Free Software Foundation, Inc.
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

/* $Id: BlurFilter.h,v 1.2 2007/08/27 18:13:39 cmusick Exp $ */

#ifndef GNASH_BLURFILTER_H
#define GNASH_BLURFILTER_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "BitmapFilter.h"

namespace gnash {

class BlurFilter_as; // Adapater for ActionScript

// A blur effect filter.
class BlurFilter : public BitmapFilter
{
public:
    friend class BlurFilter_as;

    // Fill from a stream. See parser/filter_factory.cpp for the implementations.
    virtual bool read(stream* in);

    virtual ~BlurFilter() { return; }

    // Clone this object and return a copy of it. (AS accessible function.)
    // Guaranteed to return an object which can be cast to BlurFilter
    Filter const clone();

    BlurFilter(as_object* obj) : BitmapFilter(obj),
        m_blurX(0.0f), m_blurY(0.0f), m_quality(0)
    { return; }

    BlurFilter() : 
        m_blurX(0.0f), m_blurY(0.0f), m_quality(0)
    { return; }

    BlurFilter(float blurX, float blurY, uint8_t quality) :
        m_blurX(blurX), m_blurY(blurY), m_quality(quality)
    { return; }

private:
    float m_blurX; // How much horizontal blur.
    float m_blurY; // How much vertical blur.
    uint8_t m_quality; // How many passes to take.
};

} // Namespace gnash

#endif // GNASH_BLURFILTER_H
