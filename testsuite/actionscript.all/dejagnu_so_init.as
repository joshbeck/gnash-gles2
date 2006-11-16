
// NOTE: when using ming-0.4.0-beta, a bug in 'makeswf' will
//       prevent __shared_assets clip to work (the movieclip
//       will be published with a frame-count of 0, thus
//       actions in it will *NOT* be executed)


// By default 'makeswf' makes the __shared_assets clip invisible,
// make it visible to *see* visual traces
if ( __shared_assets != undefined )
{
	if ( ! __shared_assets.dejagnu_module_initialized )
	{
		trace("No initialized dejagnu module found in __shared_assets"
			+ ": used a bogus Ming version to build Dejagnu.swf ? " 
			+ "(it's bogus up to 0.4.0beta2)");
	}
	__shared_assets._visible = true;
}
else
{
	trace("__shared_assets undefined: did you run 'makeswf' with -i<path_to_Dejagnu.swf>:dejagnu ?");
}
