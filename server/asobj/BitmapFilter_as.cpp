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

/* $Id: BitmapFilter_as.cpp,v 1.2 2007/08/27 18:13:40 cmusick Exp $ */

#include "BitmapFilter.h"
#include "VM.h"
#include "builtin_function.h"

#define phelp_helper BitmapFilter_as
#define phelp_class BitmapFilter
#include "prophelper.h"

namespace gnash {

class BitmapFilter_as
{
    phelp_base_def;
public:
    phelp_i(bitmap_clone);
};

phelp_base_imp( , BitmapFilter);

phelp_i_attach_begin
phelp_i_attach(clone, bitmap_clone);
phelp_i_attach_end

as_value
BitmapFilter_as::ctor(const fn_call& /*fn*/)
{
    boost::intrusive_ptr<as_object> obj = new BitmapFilter(BitmapFilter_as::Interface());
    return as_value(obj.get());
}

as_value BitmapFilter_as::bitmap_clone(const fn_call& fn)
{
    boost::intrusive_ptr<BitmapFilter> filter = ensureType<BitmapFilter> (fn.this_ptr);
    boost::intrusive_ptr<as_object> retval = filter->clone();
    retval->set_prototype(filter->get_prototype());

    return as_value(retval);
}

as_object*
bitmapFilter_interface()
{
    return BitmapFilter_as::Interface();
}

} // Namespace gnash

