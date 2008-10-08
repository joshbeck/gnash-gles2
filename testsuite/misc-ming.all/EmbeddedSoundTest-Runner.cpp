/* 
 *   Copyright (C) 2007, 2008 Free Software Foundation, Inc.
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 *
 */ 

#define INPUT_FILENAME "EmbeddedSoundTest.swf"

#include "MovieTester.h"
#include "sprite_instance.h"
#include "character.h"
#include "DisplayList.h"
#include "log.h"
#include "GnashSleep.h"
#include "VM.h"

#include "check.h"
#include <string>
#include <cassert>

using namespace gnash;
using namespace std;

int
main(int /*argc*/, char** /*argv*/)
{
	string filename = string(TGTDIR) + string("/") + string(INPUT_FILENAME);
	MovieTester tester(filename);

	gnash::LogFile& dbglogfile = gnash::LogFile::getDefaultInstance();
	dbglogfile.setVerbosity(1);

	sprite_instance* root = tester.getRootMovie();
	assert(root);

    VM& vm = root->getVM();
    string_table& st = vm.getStringTable();

    const size_t framecount = root->get_frame_count();

    // 20 x 0.75 seconds = 15 in total. The sound should last
    // 13.74 seconds.
    for (size_t i = 0; i < 20; ++i)
    {
        while (true)
        {
            if (root->get_current_frame() + 1 == framecount) break;
            tester.advance();
	
        }
        gnashSleep(750000);
    }
    as_value eot;
    check(root->get_member(st.find("finished"), &eot));

}

