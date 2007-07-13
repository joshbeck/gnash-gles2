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
//
// along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
//
// Test case for TextField ActionScript class
// compile this test case with Ming makeswf, and then
// execute it like this gnash -1 -r 0 -v out.swf

rcsid="$Id: TextField.as,v 1.5 2007/07/13 20:23:28 strk Exp $";

#include "check.as"

#if OUTPUT_VERSION > 5

check_equals(typeof(TextField), 'function');
check_equals(typeof(TextField.prototype), 'object');
check_equals(typeof(TextField.prototype.setTextFormat), 'function');
check_equals(typeof(TextField.prototype.getTextFormat), 'function');
check_equals(typeof(TextField.prototype.setNewTextFormat), 'function');
check_equals(typeof(TextField.prototype.getNewTextFormat), 'function');
check_equals(typeof(TextField.prototype.addListener), 'function');
check_equals(typeof(TextField.prototype.removeListener), 'function');
check_equals(typeof(TextField.prototype.getDepth), 'function');
check_equals(typeof(TextField.prototype.removeTextField), 'function');
check_equals(typeof(TextField.prototype.replaceSel), 'function');

// NOTE: the following will be true after a call to createTextField ! Seek forward to see..
check( !TextField.prototype.hasOwnProperty('background'));
check( !TextField.prototype.hasOwnProperty('backgroundColor'));
check( !TextField.prototype.hasOwnProperty('autoSize') );
check( !TextField.prototype.hasOwnProperty('border') );
check( !TextField.prototype.hasOwnProperty('borderColor') );
check( !TextField.prototype.hasOwnProperty('bottomScroll') );
check( !TextField.prototype.hasOwnProperty('embedFonts') );
check( !TextField.prototype.hasOwnProperty('hscroll') );
check( !TextField.prototype.hasOwnProperty('html') );
check( !TextField.prototype.hasOwnProperty('htmlText') );
check( !TextField.prototype.hasOwnProperty('length') );
check( !TextField.prototype.hasOwnProperty('maxChars') );
check( !TextField.prototype.hasOwnProperty('maxhscroll') );
check( !TextField.prototype.hasOwnProperty('maxscroll') );
check( !TextField.prototype.hasOwnProperty('multiline') );
check( !TextField.prototype.hasOwnProperty('password') );
check( !TextField.prototype.hasOwnProperty('restrict') );
check( !TextField.prototype.hasOwnProperty('scroll') );
check( !TextField.prototype.hasOwnProperty('selectable') );
check( !TextField.prototype.hasOwnProperty('text') );
check( !TextField.prototype.hasOwnProperty('textColor') );
check( !TextField.prototype.hasOwnProperty('textHeight') );
check( !TextField.prototype.hasOwnProperty('textWidth') );
check( !TextField.prototype.hasOwnProperty('type') );
xcheck( !TextField.prototype.hasOwnProperty('variable') );
check( !TextField.prototype.hasOwnProperty('wordWrap') );

// this is a static method
check_equals(typeof(TextField.getFontList), 'function');

check_equals(typeof(TextField.prototype.getFontList), 'undefined');

#if OUTPUT_VERSION > 6
check_equals(typeof(TextField.prototype.replaceText), 'function');
#else
check_equals(typeof(TextField.prototype.replaceText), 'undefined');
#endif

tfObj = new TextField();
check_equals(typeof(tfObj), 'object');
check(tfObj instanceof TextField);

check_equals(typeof(tfObj.setTextFormat), 'function');
check_equals(typeof(tfObj.getTextFormat), 'function');
check_equals(typeof(tfObj.setNewTextFormat), 'function');
check_equals(typeof(tfObj.getNewTextFormat), 'function');
check_equals(typeof(tfObj.addListener), 'function');
check_equals(typeof(tfObj.removeListener), 'function');
check_equals(typeof(tfObj.getDepth), 'function');
check_equals(typeof(tfObj.removeTextField), 'function');
check_equals(typeof(tfObj.replaceSel), 'function');
// this is a static method, it's available as TextField.getFontList
check_equals(typeof(tfObj.getFontList), 'undefined');

//--------------------------------------------------
// Check textfield creation trough createTextField
//--------------------------------------------------

ret = createTextField("tf", 99, 10, 10, 500, 500);
#if OUTPUT_VERSION < 8
check_equals(typeof(ret), 'undefined');
#else
xcheck_equals(typeof(ret), 'object');
xcheck_equals(ret, _root.tf);
#endif

check_equals(typeof(tf), 'object');

// NOTE: the following were false before the call to createTextField ! Seek backward to see..
xcheck( TextField.prototype.hasOwnProperty('background'));
xcheck( TextField.prototype.hasOwnProperty('backgroundColor'));
xcheck( TextField.prototype.hasOwnProperty('autoSize') );
xcheck( TextField.prototype.hasOwnProperty('border') );
xcheck( TextField.prototype.hasOwnProperty('borderColor') );
xcheck( TextField.prototype.hasOwnProperty('bottomScroll') );
xcheck( TextField.prototype.hasOwnProperty('embedFonts') );
xcheck( TextField.prototype.hasOwnProperty('hscroll') );
xcheck( TextField.prototype.hasOwnProperty('html') );
xcheck( TextField.prototype.hasOwnProperty('htmlText') );
xcheck( TextField.prototype.hasOwnProperty('length') );
xcheck( TextField.prototype.hasOwnProperty('maxChars') );
xcheck( TextField.prototype.hasOwnProperty('maxhscroll') );
xcheck( TextField.prototype.hasOwnProperty('maxscroll') );
xcheck( TextField.prototype.hasOwnProperty('multiline') );
xcheck( TextField.prototype.hasOwnProperty('password') );
xcheck( TextField.prototype.hasOwnProperty('restrict') );
xcheck( TextField.prototype.hasOwnProperty('scroll') );
xcheck( TextField.prototype.hasOwnProperty('selectable') );
xcheck( TextField.prototype.hasOwnProperty('text') );
xcheck( TextField.prototype.hasOwnProperty('textColor') );
xcheck( TextField.prototype.hasOwnProperty('textHeight') );
xcheck( TextField.prototype.hasOwnProperty('textWidth') );
xcheck( TextField.prototype.hasOwnProperty('type') );
check( TextField.prototype.hasOwnProperty('variable') );
xcheck( TextField.prototype.hasOwnProperty('wordWrap') );

// Check TextField._alpha

check_equals(typeof(tf._alpha), 'number');
check( ! tf.hasOwnProperty('_alpha') ); // why ??
check( ! tf.__proto__.hasOwnProperty('_alpha') ); // why ??

// Check TextField.autoSize

xcheck_equals(typeof(tf.autoSize), 'string');
xcheck_equals(tf.autoSize, 'none'); // TODO: research which valid values we have
check(! tf.hasOwnProperty('autoSize'));

// Check TextField.background

xcheck_equals(typeof(tf.background), 'boolean');
check(!tf.hasOwnProperty('background'));

// Check TextField.backgroundColor

xcheck_equals(typeof(tf.backgroundColor), 'number');
check(!tf.hasOwnProperty('backgroundColor'));

// Check TextField.border

xcheck_equals(typeof(tf.border), 'boolean');
check(!tf.hasOwnProperty('border'));

// Check TextField.borderColor

xcheck_equals(typeof(tf.borderColor), 'number');
check(!tf.hasOwnProperty('borderColor'));

// Check TextField.bottomScroll

xcheck_equals(typeof(tf.bottomScroll), 'number');
check(!tf.hasOwnProperty('bottomScroll'));
xcheck_equals(tf.bottomScroll, 1);
tf.bottomScroll = 100; // bottomScroll is read-only
xcheck_equals(tf.bottomScroll, 1);

// Check TextField.embedFonts

xcheck_equals(typeof(tf.embedFonts), 'boolean');
check(!tf.hasOwnProperty('embedFonts'));
xcheck_equals(tf.embedFonts, false);
// TODO: do this test with really embedded fonts, in misc-ming.all/DefineEditTextTest.c

// Check TextField._highquality

xcheck_equals(typeof(tf._highquality), 'number');
check(!tf.hasOwnProperty('_highquality'));
check(!tf.__proto__.hasOwnProperty('_highquality'));
xcheck_equals(tf._highquality, 1);
tf._highquality = 0;
check_equals(tf._highquality, 0);
tf._highquality = 1;

// Check TextField._height (how is this different from textHeight?)

check_equals(typeof(tf._height), 'number');
check(!tf.hasOwnProperty('_height'));
check(!tf.__proto__.hasOwnProperty('_height'));
xcheck_equals(tf._height, 500); // as we created it, see createTextField call
tf._height = 99999;
xcheck_equals(tf._height, 99999); 
tf._height = 500;

// Check TextField.hscroll

xcheck_equals(typeof(tf.hscroll), 'number');
check(!tf.hasOwnProperty('hscroll'));
xcheck_equals(tf.hscroll, 0);
tf.hscroll = 1;
xcheck_equals(tf.hscroll, 0);
tf.hscroll = 0;

// Check TextField.html

xcheck_equals(typeof(tf.html), 'boolean');
check(!tf.hasOwnProperty('html'));
xcheck_equals(tf.html, false);
tf.html = true;
check_equals(tf.html, true);
tf.html = false;

// Check TextField.htmlText (the displayed text in explicit HTML)

xcheck_equals(typeof(tf.htmlText), 'string');
check(!tf.hasOwnProperty('htmlText'));
xcheck_equals(tf.htmlText, '');
tf.htmlText = new Array;
xcheck_equals(typeof(tf.htmlText), 'string'); // forced cast to string
xcheck_equals(tf.htmlText, ''); 
check_equals(tf.html, false);
tf.htmlText = "Hello <b>html</b> world";
check_equals(tf.html, false); // assigning to htmlText doesn't change the 'html' flag
check_equals(tf.htmlText, 'Hello <b>html</b> world');
// Changing htmlText also changes text
xcheck_equals(tf.text, 'Hello <b>html</b> world');
tf.text = "Hello world";
xcheck_equals(tf.htmlText, 'Hello world');

// Check TextField.length  (number of characters in text)

xcheck_equals(typeof(tf.length), 'number');
check(!tf.hasOwnProperty('length'));
tf.text = "";
xcheck_equals(tf.length, 0);
tf.length = 10; // you don't change lenght like this, you assign to text instead
xcheck_equals(tf.length, 0);
tf.text = "Hello world";
xcheck_equals(tf.length, 11);
tf.htmlText = "Hello <b>world</b>";
xcheck_equals(tf.length, 18); // the tags are also counted

// Check TextField.maxChars

xcheck_equals(typeof(tf.maxChars), 'null');
check(!tf.hasOwnProperty('maxChars'));
tf.maxChars = 10;
check_equals(tf.maxChars, 10);
tf.maxChars = null;

// Check TextField.maxhscroll

xcheck_equals(typeof(tf.maxhscroll), 'number');
check(!tf.hasOwnProperty('maxhscroll'));
xcheck_equals(tf.maxhscroll, 0);
tf.maxhscroll = 10;
xcheck_equals(tf.maxhscroll, 0); // read-only

// Check TextField.maxscroll

xcheck_equals(typeof(tf.maxscroll), 'number');
check(!tf.hasOwnProperty('maxscroll'));
xcheck_equals(tf.maxscroll, 1);
tf.maxscroll = 10;
xcheck_equals(tf.maxscroll, 1); // read-only

// Check TextField.multiline

xcheck_equals(typeof(tf.multiline), 'boolean');
check(!tf.hasOwnProperty('multiline'));
xcheck_equals(tf.multiline, false);
tf.multiline = true;
check_equals(tf.multiline, true);
tf.multiline = false;

// Check TextField._name

xcheck_equals(typeof(tf._name), 'string');
check(!tf.hasOwnProperty('_name'));
check(!tf.__proto__.hasOwnProperty('_name'));
xcheck_equals(tf._name, 'tf');
// TODO: see effects of changing _name

// Check TextField._parent

xcheck_equals(typeof(tf._parent), 'movieclip');
check(!tf.hasOwnProperty('_parent'));
check(!tf.__proto__.hasOwnProperty('_parent'));
xcheck_equals(tf._parent, _root);

// Check TextField.password

xcheck_equals(typeof(tf.password), 'boolean');
check(!tf.hasOwnProperty('password'));
xcheck_equals(tf.password, false);
tf.password = true;
check_equals(tf.password, true);
// TODO: check effects of setting to 'password' (should hide characters)
tf.password = false;

// Check TextField.quality

// TODO: check this, might be a string
xcheck_equals(typeof(tf._quality), 'string');
check(!tf.hasOwnProperty('quality'));
check(!tf.__proto__.hasOwnProperty('quality'));
check(!tf.__proto__.__proto__.hasOwnProperty('quality'));
check(!tf.__proto__.__proto__.__proto__.hasOwnProperty('quality'));
xcheck_equals(tf._quality, "HIGH");
tf._quality = "FAKE VALUE";
xcheck_equals(tf._quality, "HIGH");
tf._quality = "LOW";
check_equals(tf._quality, "LOW");

// Check TextField.restrict (the set of characters a user can input)

xcheck_equals(typeof(tf.restrict), 'null');
check(!tf.hasOwnProperty('restrict'));
tf.text = "Hello World";
tf.restrict = "olH";
check_equals(tf.text, "Hello World");
tf.text = "Hello World"; // override
// doesn't influence explicit setting, only GUI modification
// of the textfield (TODO: test with a MovieTester)
check_equals(tf.text, "Hello World");


// Check TextField._rotation

xcheck_equals(typeof(tf._rotation), 'number');
check(!tf.hasOwnProperty('_rotation'));
check(!tf.__proto__.hasOwnProperty('_rotation'));
xcheck_equals(tf._rotation, 0);
tf._rotation = 10;
check_equals(tf._rotation, 10);
tf._rotation = 0;

// Check TextField.scroll

// TODO: better test for this, might do nothing if there's no scrollin
xcheck_equals(typeof(tf.scroll), 'number');
check( ! tf.hasOwnProperty('scroll') ); 
xcheck_equals(tf.scroll, 1);
tf.scroll = 10;
xcheck_equals(tf.scroll, 1); // read-only

// Check TextField.selectable

xcheck_equals(typeof(tf.selectable), 'boolean');
check( ! tf.hasOwnProperty('selectable') ); 
xcheck_equals(tf.selectable, true);
tf.selectable = false;
check_equals(tf.selectable, false);
tf.selectable = true;

// Check TextField._soundbuftime

xcheck_equals(typeof(tf._soundbuftime), 'number');
check( ! tf.hasOwnProperty('_soundbuftime') ); 
check( ! tf.__proto__.hasOwnProperty('_soundbuftime') ); 
xcheck_equals(tf._soundbuftime, 5); // the default is 5, it seems

// Check TextField.tabEnabled

check_equals(typeof(tf.tabEnabled), 'undefined');
check( ! tf.hasOwnProperty('tabEnabled') ); 
check( ! tf.__proto__.hasOwnProperty('tabEnabled') ); 
tf.tabEnabled = false;
check_equals(tf.tabEnabled, false);
delete(tf.tabEnabled);

// Check TextField.tabIndex

check_equals(typeof(tf.tabIndex), 'undefined');
check( ! tf.hasOwnProperty('tabIndex') ); 
check( ! tf.__proto__hasOwnProperty('tabIndex') ); 
tf.tabIndex = 9;
check_equals(tf.tabIndex, 9);
delete(tf.tabIndex);

// Check TextField._target

xcheck_equals(typeof(tf._target), 'string');
check( ! tf.hasOwnProperty('_target') ); 
check( ! tf.__proto__.hasOwnProperty('_target') ); 
xcheck_equals(tf._target, '/tf');
// TODO: check the effect of changing _name on the _target value
tf._target = "fake_target"; // read-only
xcheck_equals(tf._target, '/tf');

// Check TextField.text

check_equals(typeof(tf.text), 'string');
check( ! tf.hasOwnProperty('text') ); 
check_equals(tf.text, 'Hello World');
tf.text = "hello world";
check_equals(tf.text, 'hello world');
xcheck_equals(tf.length, 11); // number of characters in "hello world"


// Check TextField.textColor

check_equals(typeof(tf.textColor), 'number');
check( ! tf.hasOwnProperty('textColor') ); 
xcheck_equals(tf.textColor, 0);
tf.textColor = 0xFF0000;
xcheck_equals(tf.textColor, 0xFF0000);
// TODO: check color (use misc-ming.all/DefineEditTextTest.swf and a test runner with check_pixel)

// Check TextField.textHeight (height of the bounding box)

xcheck_equals(typeof(tf.textHeight), 'number');
check( ! tf.hasOwnProperty('textHeight') ); 
currentHeight = tf.textHeight; // WARNING: this might depend on the default font height
tf.textHeight = 1000;
xcheck_equals(tf.textHeight, currentHeight); // was read-only (I think)

// Check TextField.textWidth (width of the bounding box)

check_equals(typeof(tf.textWidth), 'number');
check( ! tf.hasOwnProperty('textWidth') ); 
currentWidth = tf.textWidth; // WARNING: this might depend on the default font height
tf.textWidth = 1000;
check_equals(tf.textWidth, currentWidth); // was read-only (I think)

// Check TextField.type (input or dynamic)

xcheck_equals(typeof(tf.type), 'string');
check( ! tf.hasOwnProperty('type') ); 
xcheck_equals(tf.type, 'dynamic'); 
tf.type = "input";
check_equals(tf.type, 'input'); 
tf.type = new Array();
xcheck_equals(typeof(tf.type), 'string');  // invalid assignment
xcheck_equals(tf.type, 'input');  // keeps previous value
tf.type = "dynamic";
check_equals(tf.type, 'dynamic');  
tf.type = new Array();
xcheck_equals(tf.type, 'dynamic'); // keeps previous value 

// Check TextField._url (url of the SWF that created the textfield)

xcheck_equals(typeof(tf._url), 'string');
check( ! tf.hasOwnProperty('_url') ); 
check( ! tf.__proto__.hasOwnProperty('_url') ); 
xcheck_equals(tf._url, _root._url);
tf._url = "fake url";
xcheck_equals(tf._url, _root._url); // read-only

// Check TextField.variable (variable name associated with the textfield)

xcheck_equals(typeof(tf.variable), 'null');
check( ! tf.hasOwnProperty('variable') ); 
tf.variable = _level0.inputVar;
xcheck_equals(typeof(tf.variable), 'null'); // read-only
// TODO: test 'variable' in misc-ming.all/DefineEditTextVariableName.c,
//       as it seems the 'variable' member is not settable on a dynamic
//       text field.


// Check TextField._visible 

check_equals(typeof(tf._visible), 'boolean');
check( ! tf.hasOwnProperty('_visible') ); 
check( ! tf.__proto__.hasOwnProperty('_visible') ); 
check_equals(tf._visible, true);
tf._visible = false;
check_equals(tf._visible, false);
tf._visible = true;

// Check TextField._width (how is this different from textWidth ?)

check_equals(typeof(tf._width), 'number');
check( ! tf.hasOwnProperty('_width') ); 
check( ! tf.__proto__.hasOwnProperty('_width') ); 
xcheck_equals(tf._width, 500); // as it was set by createTextField, see above
tf._width = 99999;
xcheck_equals(tf._width, 99999); 
tf._width = 500;

// Check TextField.wordWrap (should text wrap when bbox limit is hit?)

xcheck_equals(typeof(tf.wordWrap), 'boolean');
check( ! tf.hasOwnProperty('wordWrap') ); 
xcheck_equals(tf.wordWrap, false);

// Check TextField._x 

check_equals(typeof(tf._x), 'number');
check( ! tf.hasOwnProperty('_x') );
check( ! tf.__proto__.hasOwnProperty('_x') );
check_equals(tf._x, 10); // as set by createTextField
tf._x = 20;
check_equals(tf._x, 20);

// Check TextField._xmouse

xcheck_equals(typeof(tf._xmouse), 'number');
check( ! tf.hasOwnProperty('_xmouse') );
check( ! tf.__proto__.hasOwnProperty('_xmouse') );
currXmouse = tf._xmouse; // unsafe, if user moves the mouse while running the test
tf._xmouse = "a string";
xcheck_equals(typeof(tf._xmouse), 'number');
xcheck_equals(tf._xmouse, currXmouse); // possibly unsafe, if user moves the mouse while running the test

// Check TextField._xscale

xcheck_equals(typeof(tf._xscale), 'number');
check( ! tf.hasOwnProperty('_xscale') );
check( ! tf.__proto__.hasOwnProperty('_xscale') );
xcheck_equals(tf._xscale, 100); 
// check how .textWidth and ._width change when changing _xscale
currTextWidth = tf.textWidth;
currWidth = tf._width;
tf._xscale = 200;
note("textWidth: _xscale=100: "+currTextWidth+"; _xscale=200: "+tf.textWidth);
// check_equals(tf.textWidth, currTextWidth*2); // not clear what does textWidth depend on
xcheck_equals(tf._width, currWidth*2);
tf._xscale = 100;

// Check TextField._y 

check_equals(typeof(tf._y), 'number');
check( ! tf.hasOwnProperty('_y') );
check( ! tf.__proto__.hasOwnProperty('_y') );
check_equals(tf._y, 10); // as set by createTextField
tf._y = 5;
check_equals(tf._y, 5);

// Check TextField._ymouse

xcheck_equals(typeof(tf._ymouse), 'number');
check( ! tf.hasOwnProperty('_ymouse') );
check( ! tf.__proto__.hasOwnProperty('_ymouse') );
currYmouse = tf._ymouse;
tf._ymouse = "a string";
xcheck_equals(typeof(tf._ymouse), 'number');
xcheck_equals(tf._ymouse, currYmouse); // possibly unsafe, if user moves the mouse while running the test

// Check TextField._yscale

xcheck_equals(typeof(tf._yscale), 'number');
check( ! tf.hasOwnProperty('_yscale') );
check( ! tf.__proto__.hasOwnProperty('_yscale') );
xcheck_equals(tf._yscale, 100); 
// check how .textHeight and ._height change based on _yscale
currTextHeight = tf.textHeight;
currHeight = tf._height;
tf._yscale = 200;
note("textHeight: _yscale=100: "+currTextHeight+"; _yscale=200: "+tf.textHeight);
// check_equals(tf.textHeight, currTextHeight*2); // not clear what does textHeight depend on
xcheck_equals(tf._height, currHeight*2);
tf._yscale = 100;


#endif // OUTPUT_VERSION > 5
