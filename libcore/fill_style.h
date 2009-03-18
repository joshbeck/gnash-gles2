// 
//   Copyright (C) 2007, 2008, 2009 Free Software Foundation, Inc.
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

// Based on work of Thatcher Ulrich <tu@tulrich.com> 2003


#ifndef GNASH_FILL_STYLE_H
#define GNASH_FILL_STYLE_H

#include "smart_ptr.h" // GNASH_USE_GC
#include "SWFMatrix.h"
#include "BitmapInfo.h"
#include "swf.h"
#include "RGBA.h" // for rgba type

#include <vector> // for composition

namespace gnash {

class SWFStream;
class movie_definition;

class gradient_record
{
public:
	gradient_record()
		:
		m_ratio(0),
		m_color()
	{}

	gradient_record(boost::uint8_t ratio, const rgba& color)
		:
		m_ratio(ratio),
		m_color(color)
	{}

	void read(SWFStream& in, SWF::TagType tag);
	
	//data:
	boost::uint8_t	m_ratio;
	rgba	m_color;
};


/// For the interior of outline shapes.
//
class DSOEXPORT fill_style 
{
public:

	/// Create a solid opaque white fill.
	fill_style();

	/// Construct a clipped bitmap fill style, for
	/// use by bitmap shape character.
	///
	/// TODO: use a subclass for this
	/// TODO: provide a setBitmap, for consisteny with other setType() methods
	///
	/// @param bitmap
	///	The bitmap character definition to use with this bitmap fill.
	///
	/// @param mat
	///	The SWFMatrix to apply to the bitmap.
	///
	fill_style(BitmapInfo* bitmap, const SWFMatrix& mat);

	void setSolid(const rgba& color);

	/// Turn this fill style into a linear gradient
	//
	/// @param gradients
	///	Gradient records.
	///
	/// @param mat
	///	Gradient SWFMatrix.
	///
	///
	void setLinearGradient(const std::vector<gradient_record>& gradients, 
			const SWFMatrix& mat);

	/// Turn this fill style into a radial gradient
	//
	/// @param gradients
	///	Gradient records.
	///
	/// @param mat
	///	Gradient SWFMatrix.
	///
	///
	void setRadialGradient(const std::vector<gradient_record>& gradients,
			const SWFMatrix& mat);

	/// Turn this fill style into a focal gradient
	//
	/// @param gradients
	///	Gradient records.
	///
	/// @param mat
	///	Gradient SWFMatrix.
	///
	/// @param fpoint
	///	Focal point.
	///
	void setRadialGradient(const std::vector<gradient_record>& gradients,
			const SWFMatrix& mat, float fpoint);

	~fill_style() {}
	
	/// Read the fill style from a stream
	//
	/// TODO: use a subclass for this (swf_fill_style?)
	///
	/// Throw a ParserException if there's no enough bytes in the
	/// currently opened tag for reading. See stream::ensureBytes()
	///
	void read(SWFStream& in, SWF::TagType t, movie_definition& m,
		fill_style *pOther = NULL);

	/// Read the fill style from a stream, morph version.
	void read_morph(SWFStream& in, SWF::TagType t, movie_definition& m,
		fill_style *pOther);

	/// \brief
	/// Make a BitmapInfo* corresponding to our gradient.
	/// We can use this to set the gradient fill style.
	BitmapInfo* create_gradient_bitmap() const;
	
	/// \brief
	/// Makes sure that _gradientBitmapInfo is not NULL. Calls 
	/// create_gradient_bitmap() if necessary and returns _gradientBitmapInfo.
	BitmapInfo* need_gradient_bitmap() const; 
	
	rgba	get_color() const { return m_color; }

	void	set_color(rgba new_color) { m_color = new_color; }

	/// Get fill type, see SWF::fill_style_type
	uint8_t	get_type() const { return m_type; }

	SWF::gradient_spread_mode get_gradient_spread_mode()
	{ return m_spread_mode; }

	SWF::gradient_interpolation_mode get_gradient_interpolation_mode()
	{ return m_interpolation; }
	
	/// Sets this style to a blend of a and b.  t = [0,1] (for shape morphing)
	void	set_lerp(const fill_style& a, const fill_style& b, float t);
	
	/// Returns the bitmap info for all styles except solid fills
	//
	/// NOTE: calling this method against a solid fill style will
	///       result in a failed assertion.
	/// 
	/// NOTE2: this function can return NULL if the character_id
	///        specified for the style in the SWF does not resolve
	///        to a character defined in the characters dictionary.
	///        (it happens..)
	///
	BitmapInfo* get_bitmap_info() const;
	
	/// Returns the bitmap transformation SWFMatrix
	const SWFMatrix& getBitmapMatrix() const; 
	
	/// Returns the gradient transformation SWFMatrix
	const SWFMatrix& getGradientMatrix() const; 
	
	/// Returns the number of color stops in the gradient
	int get_color_stop_count() const;
	
	/// Returns the color stop value at a specified index
	const gradient_record& get_color_stop(int index) const;

	/// Get and set the focal point for gradient focal fills.
	/// This should be from -1.0 to 1.0, representing the left
	/// and right edges of the rectangle.
	float get_focal_point() const { return m_focal_point; }
	void set_focal_point(float f) { m_focal_point = f; }

#ifdef GNASH_USE_GC
	/// Mark reachable resources (for the GC)
	//
	/// fill_style specific reachable resources are:
	///
	///	- gradient bitmap info (_gradientBitmapInfo)
	///	- bitmap character (m_bitmap_character)
	///
	void markReachableResources() const;
#endif // GNASH_USE_GC

private:

	/// Return the color at the specified ratio into our gradient.
	//
	/// @param ratio
	///	Ratio is in the range [0, 255].
	///
	rgba sample_gradient(boost::uint8_t ratio) const;

	friend class morph2_character_def;

	/// Fill type, see SWF::fill_style_type
	uint8_t	m_type;
	
	// For BITMAP or GRADIENT types 
	SWFMatrix	_matrix;

	// For BITMAP or GRADIENT types
	boost::intrusive_ptr<BitmapInfo> _bitmapInfo;

	// For SOLID type (and arguably GRADIENT too)
	rgba	m_color;

	// Only for GRADIENT type
	float m_focal_point; // For focal fill gradients.
	std::vector<gradient_record> m_gradients;
	SWF::gradient_spread_mode m_spread_mode;
	SWF::gradient_interpolation_mode m_interpolation;

};


} // namespace gnash


#endif // GNASH_FILL_STYLE_H


// Local Variables:
// mode: C++
// indent-tabs-mode: t
// End:
