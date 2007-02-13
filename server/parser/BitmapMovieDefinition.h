// 
//   Copyright (C) 2005, 2006 Free Software Foundation, Inc.
// 
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

// 
//


#ifndef GNASH_BITMAPMOVIEDEFINITION_H
#define GNASH_BITMAPMOVIEDEFINITION_H

#include "movie_definition.h" // for inheritance
#include "rect.h" // for composition
#include "BitmapMovieInstance.h" // for create_instance

#include <string>
#include <memory> // for auto_ptr

// Forward declarations
namespace gnash {
	class bitmap_character_def;
}

namespace gnash
{

/// A "movie" definition for a bitmap file
//
/// The create_instance function will return a BitmapMovieInstance
///
class BitmapMovieDefinition : public movie_definition
{
	int _version;
	rect _framesize;
	size_t _framecount;
	std::vector<PlayList> _playlist;
	float _framerate;
	std::string _url;
	boost::intrusive_ptr<bitmap_character_def> _bitmap;

public:


	/// Default constructor
	//
	/// Will be initialized with the following values
	///
	///  - SWF version 6
	///  - 640x480 size
	///  - Single frame (unlabeled)
	///  - 12 FPS
	///  - 0 bytes (for get_bytes_loaded()/get_bytes_total())
	///  - empty url
	///
	BitmapMovieDefinition(boost::intrusive_ptr<bitmap_character_def> bi, const std::string& url)
		:
		_version(6),
		// FIXME: extract size from the bitmap_character
		_framesize(0, 0, 640*20, 480*20),
		_framecount(1),
		_playlist(_framecount),
		_framerate(12),
		_url(url),
		_bitmap(bi)
	{
	}

	bitmap_character_def* get_bitmap_char_def()
	{
		return _bitmap.get();
	}

	virtual int	get_version() const {
		return _version;
	}

	virtual float	get_width_pixels() const {
		return _framesize.width()/20;
	}

	virtual float	get_height_pixels() const {
		return _framesize.height()/20;
	}

	virtual size_t	get_frame_count() const {
		return _framecount;
	}

	virtual float	get_frame_rate() const {
		return _framerate;
	}

	virtual const rect& get_frame_size() const {
		return _framesize;
	}

	virtual const rect& get_bound() const {
		return _framesize;
	}

	virtual size_t get_bytes_loaded() const {
		return 0;
	}

	virtual size_t get_bytes_total() const {
		return 0;
	}
	
	/// Create a playable sprite instance from a def.
	virtual sprite_instance* create_instance()
	{
		return create_movie_instance();
	}

	/// Create a playable movie_instance from this def.
	virtual movie_instance* create_movie_instance()
	{
		return new BitmapMovieInstance(this, NULL);
	}

	virtual const PlayList& get_playlist(size_t frame_number) {
		assert ( frame_number < _playlist.size() );
		return _playlist[frame_number];
	}

	virtual const std::string& get_url() const {
		return _url;
	}

	// Inheritance from movie_definition requires this
	// TODO: provide a default implementation in movie_definition instead ?
	size_t  get_loading_frame() const 
	{
		return 0;
	}
};

} // namespace gnash

#endif // GNASH_BITMAPMOVIEDEFINITION_H
