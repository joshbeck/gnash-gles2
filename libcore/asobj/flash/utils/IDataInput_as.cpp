// IDataInput_as.cpp:  ActionScript "IDataInput" class, for Gnash.
//
//   Copyright (C) 2009 Free Software Foundation, Inc.
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
//

#ifdef HAVE_CONFIG_H
#include "gnashconfig.h"
#endif

#include "utils/IDataInput_as.h"
#include "log.h"
#include "fn_call.h"
#include "smart_ptr.h" // for boost intrusive_ptr
#include "builtin_function.h" // need builtin_function
#include "GnashException.h" // for ActionException

namespace gnash {

// Forward declarations
namespace {
    as_value idatainput_readByte(const fn_call& fn);
    as_value idatainput_readBytes(const fn_call& fn);
    as_value idatainput_readDouble(const fn_call& fn);
    as_value idatainput_readFloat(const fn_call& fn);
    as_value idatainput_readInt(const fn_call& fn);
    as_value idatainput_readMultiByte(const fn_call& fn);
    as_value idatainput_readObject(const fn_call& fn);
    as_value idatainput_readShort(const fn_call& fn);
    as_value idatainput_readUnsignedByte(const fn_call& fn);
    as_value idatainput_readUnsignedInt(const fn_call& fn);
    as_value idatainput_readUnsignedShort(const fn_call& fn);
    as_value idatainput_readUTF(const fn_call& fn);
    as_value idatainput_readUTFBytes(const fn_call& fn);
    as_value idatainput_ctor(const fn_call& fn);
    void attachIDataInputInterface(as_object& o);
    void attachIDataInputStaticInterface(as_object& o);
    as_object* getIDataInputInterface();

}

class IDataInput_as : public as_object
{

public:

    IDataInput_as()
        :
        as_object(getIDataInputInterface())
    {}
};

// extern (used by Global.cpp)
void idatainput_class_init(as_object& global)
{
    static boost::intrusive_ptr<builtin_function> cl;

    if (!cl) {
        cl = new builtin_function(&idatainput_ctor, getIDataInputInterface());
        attachIDataInputStaticInterface(*cl);
    }

    // Register _global.IDataInput
    global.init_member("IDataInput", cl.get());
}

namespace {

void
attachIDataInputInterface(as_object& o)
{
    o.init_member("readByte", new builtin_function(idatainput_readByte));
    o.init_member("readBytes", new builtin_function(idatainput_readBytes));
    o.init_member("readDouble", new builtin_function(idatainput_readDouble));
    o.init_member("readFloat", new builtin_function(idatainput_readFloat));
    o.init_member("readInt", new builtin_function(idatainput_readInt));
    o.init_member("readMultiByte", new builtin_function(idatainput_readMultiByte));
    o.init_member("readObject", new builtin_function(idatainput_readObject));
    o.init_member("readShort", new builtin_function(idatainput_readShort));
    o.init_member("readUnsignedByte", new builtin_function(idatainput_readUnsignedByte));
    o.init_member("readUnsignedInt", new builtin_function(idatainput_readUnsignedInt));
    o.init_member("readUnsignedShort", new builtin_function(idatainput_readUnsignedShort));
    o.init_member("readUTF", new builtin_function(idatainput_readUTF));
    o.init_member("readUTFBytes", new builtin_function(idatainput_readUTFBytes));
}

void
attachIDataInputStaticInterface(as_object& o)
{

}

as_object*
getIDataInputInterface()
{
    static boost::intrusive_ptr<as_object> o;
    if ( ! o ) {
        o = new as_object();
        attachIDataInputInterface(*o);
    }
    return o.get();
}

as_value
idatainput_readByte(const fn_call& fn)
{
    boost::intrusive_ptr<IDataInput_as> ptr =
        ensureType<IDataInput_as>(fn.this_ptr);
    UNUSED(ptr);
    log_unimpl (__FUNCTION__);
    return as_value();
}

as_value
idatainput_readBytes(const fn_call& fn)
{
    boost::intrusive_ptr<IDataInput_as> ptr =
        ensureType<IDataInput_as>(fn.this_ptr);
    UNUSED(ptr);
    log_unimpl (__FUNCTION__);
    return as_value();
}

as_value
idatainput_readDouble(const fn_call& fn)
{
    boost::intrusive_ptr<IDataInput_as> ptr =
        ensureType<IDataInput_as>(fn.this_ptr);
    UNUSED(ptr);
    log_unimpl (__FUNCTION__);
    return as_value();
}

as_value
idatainput_readFloat(const fn_call& fn)
{
    boost::intrusive_ptr<IDataInput_as> ptr =
        ensureType<IDataInput_as>(fn.this_ptr);
    UNUSED(ptr);
    log_unimpl (__FUNCTION__);
    return as_value();
}

as_value
idatainput_readInt(const fn_call& fn)
{
    boost::intrusive_ptr<IDataInput_as> ptr =
        ensureType<IDataInput_as>(fn.this_ptr);
    UNUSED(ptr);
    log_unimpl (__FUNCTION__);
    return as_value();
}

as_value
idatainput_readMultiByte(const fn_call& fn)
{
    boost::intrusive_ptr<IDataInput_as> ptr =
        ensureType<IDataInput_as>(fn.this_ptr);
    UNUSED(ptr);
    log_unimpl (__FUNCTION__);
    return as_value();
}

as_value
idatainput_readObject(const fn_call& fn)
{
    boost::intrusive_ptr<IDataInput_as> ptr =
        ensureType<IDataInput_as>(fn.this_ptr);
    UNUSED(ptr);
    log_unimpl (__FUNCTION__);
    return as_value();
}

as_value
idatainput_readShort(const fn_call& fn)
{
    boost::intrusive_ptr<IDataInput_as> ptr =
        ensureType<IDataInput_as>(fn.this_ptr);
    UNUSED(ptr);
    log_unimpl (__FUNCTION__);
    return as_value();
}

as_value
idatainput_readUnsignedByte(const fn_call& fn)
{
    boost::intrusive_ptr<IDataInput_as> ptr =
        ensureType<IDataInput_as>(fn.this_ptr);
    UNUSED(ptr);
    log_unimpl (__FUNCTION__);
    return as_value();
}

as_value
idatainput_readUnsignedInt(const fn_call& fn)
{
    boost::intrusive_ptr<IDataInput_as> ptr =
        ensureType<IDataInput_as>(fn.this_ptr);
    UNUSED(ptr);
    log_unimpl (__FUNCTION__);
    return as_value();
}

as_value
idatainput_readUnsignedShort(const fn_call& fn)
{
    boost::intrusive_ptr<IDataInput_as> ptr =
        ensureType<IDataInput_as>(fn.this_ptr);
    UNUSED(ptr);
    log_unimpl (__FUNCTION__);
    return as_value();
}

as_value
idatainput_readUTF(const fn_call& fn)
{
    boost::intrusive_ptr<IDataInput_as> ptr =
        ensureType<IDataInput_as>(fn.this_ptr);
    UNUSED(ptr);
    log_unimpl (__FUNCTION__);
    return as_value();
}

as_value
idatainput_readUTFBytes(const fn_call& fn)
{
    boost::intrusive_ptr<IDataInput_as> ptr =
        ensureType<IDataInput_as>(fn.this_ptr);
    UNUSED(ptr);
    log_unimpl (__FUNCTION__);
    return as_value();
}

as_value
idatainput_ctor(const fn_call& fn)
{
    boost::intrusive_ptr<as_object> obj = new IDataInput_as;

    return as_value(obj.get()); // will keep alive
}

} // anonymous namespace 
} // gnash namespace

// local Variables:
// mode: C++
// indent-tabs-mode: t
// End:
