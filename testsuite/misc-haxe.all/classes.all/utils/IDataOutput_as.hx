// IDataOutput_as.hx:  ActionScript 3 "IDataOutput" class, for Gnash.
//
// Generated by gen-as3.sh on: 20090515 by "rob". Remove this
// after any hand editing loosing changes.
//
//   Copyright (C) 2009, 2010 Free Software Foundation, Inc.
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

// This test case must be processed by CPP before compiling to include the
//  DejaGnu.hx header file for the testing framework support.

#if flash9
import flash.utils.IDataOutput;
import flash.display.MovieClip;
#end
import flash.Lib;
import Type;

// import our testing API
import DejaGnu;

// Class must be named with the _as suffix, as that's the same name as the file.
class IDataOutput_as {
    static function main() {
#if flash9
DejaGnu.note("This class is an interface");
DejaGnu.note("Implementors: ByteArray, FileStream, Socket");
// Tests to see if all the properties exist. All these do is test for
// existance of a property, and don't test the functionality at all. This
// is primarily useful only to test completeness of the API implementation.
//	if (x1.endian == null) {
//	    DejaGnu.pass("IDataOutput.endian property exists");
//	} else {
//	    DejaGnu.fail("IDataOutput.endian property doesn't exist");
//	}
//	if (x1.objectEncoding == uint) {
//	    DejaGnu.pass("IDataOutput.objectEncoding property exists");
//	} else {
//	    DejaGnu.fail("IDataOutput.objectEncoding property doesn't exist");
//	}

// Tests to see if all the methods exist. All these do is test for
// existance of a method, and don't test the functionality at all. This
// is primarily useful only to test completeness of the API implementation.
//	if (x1.writeBoolean == null) {
//	    DejaGnu.pass("IDataOutput::writeBoolean() method exists");
//	} else {
//	    DejaGnu.fail("IDataOutput::writeBoolean() method doesn't exist");
//	}
//	if (x1.writeByte == null) {
//	    DejaGnu.pass("IDataOutput::writeByte() method exists");
//	} else {
//	    DejaGnu.fail("IDataOutput::writeByte() method doesn't exist");
//	}
//	if (x1.writeBytes == null) {
//	    DejaGnu.pass("IDataOutput::writeBytes() method exists");
//	} else {
//	    DejaGnu.fail("IDataOutput::writeBytes() method doesn't exist");
//	}
//	if (x1.writeDouble == null) {
//	    DejaGnu.pass("IDataOutput::writeDouble() method exists");
//	} else {
//	    DejaGnu.fail("IDataOutput::writeDouble() method doesn't exist");
//	}
//	if (x1.writeFloat == null) {
//	    DejaGnu.pass("IDataOutput::writeFloat() method exists");
//	} else {
//	    DejaGnu.fail("IDataOutput::writeFloat() method doesn't exist");
//	}
//	if (x1.writeInt == null) {
//	    DejaGnu.pass("IDataOutput::writeInt() method exists");
//	} else {
//	    DejaGnu.fail("IDataOutput::writeInt() method doesn't exist");
//	}
//	if (x1.writeMultiByte == null) {
//	    DejaGnu.pass("IDataOutput::writeMultiByte() method exists");
//	} else {
//	    DejaGnu.fail("IDataOutput::writeMultiByte() method doesn't exist");
//	}
//	if (x1.writeObject == null) {
//	    DejaGnu.pass("IDataOutput::writeObject() method exists");
//	} else {
//	    DejaGnu.fail("IDataOutput::writeObject() method doesn't exist");
//	}
//	if (x1.writeShort == null) {
//	    DejaGnu.pass("IDataOutput::writeShort() method exists");
//	} else {
//	    DejaGnu.fail("IDataOutput::writeShort() method doesn't exist");
//	}
//	if (x1.writeUnsignedInt == null) {
//	    DejaGnu.pass("IDataOutput::writeUnsignedInt() method exists");
//	} else {
//	    DejaGnu.fail("IDataOutput::writeUnsignedInt() method doesn't exist");
//	}
//	if (x1.writeUTF == null) {
//	    DejaGnu.pass("IDataOutput::writeUTF() method exists");
//	} else {
//	    DejaGnu.fail("IDataOutput::writeUTF() method doesn't exist");
//	}
//	if (x1.writeUTFBytes == null) {
//	    DejaGnu.pass("IDataOutput::writeUTFBytes() method exists");
//	} else {
//	    DejaGnu.fail("IDataOutput::writeUTFBytes() method doesn't exist");
//	}

        // Call this after finishing all tests. It prints out the totals.
        DejaGnu.done();
#else
	DejaGnu.note("This class (IDataOutput) is only available in flash9");
#end
    }
}

// local Variables:
// mode: C++
// indent-tabs-mode: t
// End:

