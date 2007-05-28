// NetStreamFfmpeg.cpp:  Network streaming for FFMPEG video library, for Gnash.
// 
//   Copyright (C) 2005, 2006, 2007 Free Software Foundation, Inc.
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
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
//

/* $Id: NetStreamFfmpeg.cpp,v 1.59 2007/05/28 16:19:03 strk Exp $ */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef USE_FFMPEG

#include "NetStreamFfmpeg.h"
#include "log.h"
#include "fn_call.h"
#include "NetStream.h"
#include "render.h"	
#include "movie_root.h"
#include "NetConnection.h"
#include "sound_handler.h"
//#include "action.h"
#include <boost/scoped_array.hpp>

#if defined(_WIN32) || defined(WIN32)
# include <windows.h>	// for sleep()
# define usleep(x) Sleep(x/1000)
#else
# include "unistd.h" // for usleep()
#endif

// Used to free data in the AVPackets we create our self
static void avpacket_destruct(AVPacket* av) {
	delete [] av->data;
}


namespace gnash {


NetStreamFfmpeg::NetStreamFfmpeg():
	m_video_index(-1),
	m_audio_index(-1),

	m_VCodecCtx(NULL),
	m_ACodecCtx(NULL),
	m_FormatCtx(NULL),
	m_Frame(NULL),

	_decodeThread(NULL),

	m_last_video_timestamp(0),
	m_last_audio_timestamp(0),
	m_current_timestamp(0),
	m_unqueued_data(NULL),
	m_time_of_pause(0)
{

	ByteIOCxt.buffer = NULL;
}

NetStreamFfmpeg::~NetStreamFfmpeg()
{
	close();
	//delete m_parser;
}

void NetStreamFfmpeg::pause(int mode)
{
	if (mode == -1)
	{
		if (m_pause) unpauseDecoding();
		else pauseDecoding();
	}
	else
	{
		if (mode == 0) pauseDecoding();
		else unpauseDecoding();
	}
	if (!m_pause && !m_go) { 
		setStatus(playStart);
		m_go = true;
		_decodeThread = new boost::thread(boost::bind(NetStreamFfmpeg::av_streamer, this)); 
	}

}

void NetStreamFfmpeg::close()
{

	if (m_go)
	{
		// terminate thread
		m_go = false;
		decode_wait.notify_one();

		// wait till thread is complete before main continues
		_decodeThread->join();
		delete _decodeThread;

	}


	// When closing gnash before playback is finished, the soundhandler 
	// seems to be removed before netstream is destroyed.
	sound_handler* s = get_sound_handler();
	if (s != NULL)
	{
		s->detach_aux_streamer((void*) NULL);
	}

	// Make sure al decoding has stopped
	boost::mutex::scoped_lock lock(decoding_mutex);

	if (m_Frame) av_free(m_Frame);
	m_Frame = NULL;

	if (m_VCodecCtx) avcodec_close(m_VCodecCtx);
	m_VCodecCtx = NULL;

	if (m_ACodecCtx) avcodec_close(m_ACodecCtx);
	m_ACodecCtx = NULL;

	if (m_FormatCtx) {
		m_FormatCtx->iformat->flags = AVFMT_NOFILE;
		av_close_input_file(m_FormatCtx);
		m_FormatCtx = NULL;
	}

	delete m_imageframe;
	m_imageframe = NULL;
	delete m_unqueued_data;
	m_unqueued_data = NULL;

	while (m_qvideo.size() > 0)
	{
		delete m_qvideo.front();
		m_qvideo.pop();
	}

	while (m_qaudio.size() > 0)
	{
		delete m_qaudio.front();
		m_qaudio.pop();
	}

	delete [] ByteIOCxt.buffer;

}

// ffmpeg callback function
int 
NetStreamFfmpeg::readPacket(void* opaque, uint8_t* buf, int buf_size)
{

	NetStreamFfmpeg* ns = static_cast<NetStreamFfmpeg*>(opaque);
	boost::intrusive_ptr<NetConnection> nc = ns->_netCon;

	size_t ret = nc->read(static_cast<void*>(buf), buf_size);
	ns->inputPos += ret;
	return ret;

}

// ffmpeg callback function
offset_t 
NetStreamFfmpeg::seekMedia(void *opaque, offset_t offset, int whence){

	NetStreamFfmpeg* ns = static_cast<NetStreamFfmpeg*>(opaque);
	boost::intrusive_ptr<NetConnection> nc = ns->_netCon;


	// Offset is absolute new position in the file
	if (whence == SEEK_SET) {
		nc->seek(offset);
		ns->inputPos = offset;

	// New position is offset + old position
	} else if (whence == SEEK_CUR) {
		nc->seek(ns->inputPos + offset);
		ns->inputPos = ns->inputPos + offset;

	// 	// New position is offset + end of file
	} else if (whence == SEEK_END) {
		// This is (most likely) a streamed file, so we can't seek to the end!
		// Instead we seek to 50.000 bytes... seems to work fine...
		nc->seek(50000);
		ns->inputPos = 50000;
		
	}

	return ns->inputPos;
}

int
NetStreamFfmpeg::play(const std::string& c_url)
{

	// Is it already playing ?
	if (m_go)
	{
		if (m_pause) unpauseDecoding();
		return 0;
	}

	// Does it have an associated NetConnection ?
	if ( ! _netCon )
	{
		IF_VERBOSE_ASCODING_ERRORS(
		log_aserror(_("No NetConnection associated with this NetStream, won't play"));
		);
		return 0;
	}

	if (url.size() == 0) url += c_url;
	// Remove any "mp3:" prefix. Maybe should use this to mark as audio-only
	if (url.compare(0, 4, std::string("mp3:")) == 0) {
		url = url.substr(4);
	}

	m_go = true;
	pauseDecoding();

	// This starts the decoding thread
	_decodeThread = new boost::thread(boost::bind(NetStreamFfmpeg::av_streamer, this)); 

	return 0;
}

/// Finds a decoder, allocates a context and initializes it.
//
/// @param codec_id the codec ID to find
/// @return the initialized context, or NULL on failure. The caller is 
///         responsible for deallocating!
static AVCodecContext*
initContext(enum CodecID codec_id)
{

	AVCodec* codec = avcodec_find_decoder(codec_id);
	if (!codec) {
		log_error(_("libavcodec couldn't find decoder"));
		return NULL;
	}

	AVCodecContext * context = avcodec_alloc_context();
	if (!context) {
		log_error(_("libavcodec couldn't allocate context"));
		return NULL;
	}

	int rv = avcodec_open(context, codec);
	if (rv < 0) {
		avcodec_close(context);
		log_error(_("libavcodec failed to initialize codec"));
		return NULL;
	}

	return context;
}

/// Gets video info from the parser and initializes the codec.
//
/// @param parser the parser to use to get video information.
/// @return the initialized context, or NULL on failure. The caller
///         is responsible for deallocating this pointer.
static AVCodecContext* 
initFlvVideo(FLVParser* parser)
{
	// Get video info from the parser
	std::auto_ptr<FLVVideoInfo> videoInfo( parser->getVideoInfo() );
	if (!videoInfo.get()) {
		return NULL;
	}

	enum CodecID codec_id;

	// Find the decoder and init the parser
	switch(videoInfo->codec) {
		case VIDEO_CODEC_H263:
			codec_id = CODEC_ID_FLV1;
			break;
#ifdef FFMPEG_VP6
		case VIDEO_CODEC_VP6:
			codec_id = CODEC_ID_VP6F;
			break;
#endif
		case VIDEO_CODEC_SCREENVIDEO:
			codec_id = CODEC_ID_FLASHSV;
			break;
		default:
			log_error(_("Unsupported video codec %d"), (int) videoInfo->codec);
			return NULL;
	}

	return initContext(codec_id);
}


/// Like initFlvVideo, but for audio.
static AVCodecContext*
initFlvAudio(FLVParser* parser)
{
	// Get audio info from the parser
	std::auto_ptr<FLVAudioInfo> audioInfo( parser->getAudioInfo() );
	if (!audioInfo.get()) {
		return NULL;
	}

	enum CodecID codec_id;

	switch(audioInfo->codec) {
		case AUDIO_CODEC_RAW:
			codec_id = CODEC_ID_PCM_U16LE;
			break;
		case AUDIO_CODEC_ADPCM:
			codec_id = CODEC_ID_ADPCM_SWF;
			break;
		case AUDIO_CODEC_MP3:
			codec_id = CODEC_ID_MP3;
			break;
		default:
			log_error(_("Unsupported audio codec %d"), (int)audioInfo->codec);
			return NULL;
	}

	return initContext(codec_id);
}


/// Probe the stream and try to figure out what the format is.
//
/// @param ns the netstream to use for reading
/// @return a pointer to the AVInputFormat structure containing
///         information about the input format, or NULL.
static AVInputFormat*
probeStream(NetStreamFfmpeg* ns)
{
	boost::scoped_array<uint8_t> buffer(new uint8_t[2048]);

	// Probe the file to detect the format
	AVProbeData probe_data;
	probe_data.filename = "";
	probe_data.buf = buffer.get();
	probe_data.buf_size = 2048;

	if (ns->readPacket(ns, probe_data.buf, probe_data.buf_size) < 1){
 		log_error(_("Gnash could not read from movie url"));
 		return NULL;
	}

	return av_probe_input_format(&probe_data, 1);
}

bool
NetStreamFfmpeg::startPlayback()
{

	boost::intrusive_ptr<NetConnection> nc = _netCon;
	assert(nc);

	// Pass stuff from/to the NetConnection object.
	if ( !nc->openConnection(url) ) {
		log_error(_("Gnash could not open movie: %s"), url.c_str());
		setStatus(streamNotFound);
		return false;
	}

	inputPos = 0;

	// Check if the file is a FLV, in which case we use our own parser
	char head[4] = {0, 0, 0, 0};
	if (nc->read(head, 3) < 3) {
		setStatus(streamNotFound);
		return false;
	}

	nc->seek(0);
	if (std::string(head) == "FLV") {
		m_isFLV = true;
		if (!m_parser.get()) {
			m_parser.reset(new FLVParser()); // TODO: define ownership, use auto_ptr !
			if (!nc->connectParser(*m_parser)) {
				setStatus(streamNotFound);
				log_error(_("Gnash could not open FLV movie: %s"), url.c_str());
				m_parser.reset(); // release memory associated with parser
				return false;
			}
		}

		// Init the avdecoder-decoder
		avcodec_init();
		avcodec_register_all();

		m_VCodecCtx = initFlvVideo(m_parser.get());
		if (!m_VCodecCtx) {
			log_msg(_("Failed to initialize FLV video codec"));
			return false;
		}

		m_ACodecCtx = initFlvAudio(m_parser.get());
		if (!m_ACodecCtx) {
			log_msg(_("Failed to initialize FLV audio codec"));
			return false;
		}

		// We just define the indexes here, they're not really used when
		// the file format is FLV
		m_video_index = 0;
		m_audio_index = 1;

		sound_handler* s = get_sound_handler();
		if (s) s->attach_aux_streamer(audio_streamer, (void*) this);

		m_start_onbuffer = true;

		// Allocate a frame to store the decoded frame in
		m_Frame = avcodec_alloc_frame();
		return true;
	}

	// This registers all available file formats and codecs 
	// with the library so they will be used automatically when
	// a file with the corresponding format/codec is opened
	// XXX should we call avcodec_init() first?
	av_register_all();

	AVInputFormat* inputFmt = probeStream(this);
	if (!inputFmt) {
		log_error(_("Couldn't determine stream input format from URL %s"), url.c_str());
		return false;
	}

	// After the format probe, reset to the beginning of the file.
	nc->seek(0);

	// Setup the filereader/seeker mechanism. 7th argument (NULL) is the writer function,
	// which isn't needed.
	init_put_byte(&ByteIOCxt, new uint8_t[500000], 500000, 0, this, NetStreamFfmpeg::readPacket, NULL, NetStreamFfmpeg::seekMedia);
	ByteIOCxt.is_streamed = 1;

	m_FormatCtx = av_alloc_format_context();

	// Open the stream. the 4th argument is the filename, which we ignore.
	if(av_open_input_stream(&m_FormatCtx, &ByteIOCxt, "", inputFmt, NULL) < 0){
		log_error(_("Couldn't open file '%s' for decoding"), url.c_str());
		setStatus(streamNotFound);
		return false;
	}

	// Next, we need to retrieve information about the streams contained in the file
	// This fills the streams field of the AVFormatContext with valid information
	int ret = av_find_stream_info(m_FormatCtx);
	if (ret < 0)
	{
		log_error(_("Couldn't find stream information from '%s', error code: %d"), url.c_str(), ret);
		return false;
	}

//	m_FormatCtx->pb.eof_reached = 0;
//	av_read_play(m_FormatCtx);

	// Find the first video & audio stream
	m_video_index = -1;
	m_audio_index = -1;
	//assert(m_FormatCtx->nb_streams >= 0); useless assert. 
	for (unsigned int i = 0; i < m_FormatCtx->nb_streams; i++)
	{
		AVCodecContext* enc = m_FormatCtx->streams[i]->codec; 

		switch (enc->codec_type)
		{
			case CODEC_TYPE_AUDIO:
				if (m_audio_index < 0)
				{
					m_audio_index = i;
					m_audio_stream = m_FormatCtx->streams[i];
				}
				break;

			case CODEC_TYPE_VIDEO:
				if (m_video_index < 0)
				{
					m_video_index = i;
					m_video_stream = m_FormatCtx->streams[i];
				}
				break;
			default:
				break;
		}
	}

	if (m_video_index < 0)
	{
		log_error(_("Didn't find a video stream from '%s'"), url.c_str());
		return false;
	}

	// Get a pointer to the codec context for the video stream
	m_VCodecCtx = m_FormatCtx->streams[m_video_index]->codec;

	// Find the decoder for the video stream
	AVCodec* pCodec = avcodec_find_decoder(m_VCodecCtx->codec_id);
	if (pCodec == NULL)
	{
		m_VCodecCtx = NULL;
		log_error(_("Video decoder %d not found"), 
			m_VCodecCtx->codec_id);
		return false;
	}

	// Open codec
	if (avcodec_open(m_VCodecCtx, pCodec) < 0)
	{
		log_error(_("Could not open codec %d"),
			m_VCodecCtx->codec_id);
	}

	// Allocate a frame to store the decoded frame in
	m_Frame = avcodec_alloc_frame();
	
	// Determine required buffer size and allocate buffer
	if (m_videoFrameFormat == render::YUV) {
		m_imageframe = new image::yuv(m_VCodecCtx->width,	m_VCodecCtx->height);
	} else if (m_videoFrameFormat == render::RGB) {
		m_imageframe = new image::rgb(m_VCodecCtx->width,	m_VCodecCtx->height);
	}

	sound_handler* s = get_sound_handler();
	if (m_audio_index >= 0 && s != NULL)
	{
		// Get a pointer to the audio codec context for the video stream
		m_ACodecCtx = m_FormatCtx->streams[m_audio_index]->codec;

		// Find the decoder for the audio stream
		AVCodec* pACodec = avcodec_find_decoder(m_ACodecCtx->codec_id);
	    if(pACodec == NULL)
		{
			log_error(_("No available audio decoder %d to process MPEG file: '%s'"), 
				m_ACodecCtx->codec_id, url.c_str());
			return false;
		}
        
		// Open codec
		if (avcodec_open(m_ACodecCtx, pACodec) < 0)
		{
			log_error(_("Could not open audio codec %d for %s"),
				m_ACodecCtx->codec_id, url.c_str());
			return false;
		}

		s->attach_aux_streamer(audio_streamer, (void*) this);

	}

	unpauseDecoding();
	return true;
}


/// Copy RGB data from a source raw_mediadata_t to a destination image::rgb.
/// @param dst the destination image::rgb, which must already be initialized
///            with a buffer of size of at least src.m_size.
/// @param src the source raw_mediadata_t to copy data from. The m_size member
///            of this structure must be initialized.
/// @param width the width, in bytes, of a row of video data.
static void
rgbcopy(image::rgb* dst, raw_mediadata_t* src, int width)
{
	assert(src->m_size <= static_cast<uint32_t>(dst->m_width * dst->m_height * 3));

	uint8_t* dstptr = dst->m_data;

	uint8_t* srcptr = src->m_data;
	uint8_t* srcend = src->m_data + src->m_size;

	while (srcptr < srcend) {
		memcpy(dstptr, srcptr, width);
		dstptr += dst->m_pitch;
		srcptr += width;
	}
}

// decoder thread
void NetStreamFfmpeg::av_streamer(NetStreamFfmpeg* ns)
{
	GNASH_REPORT_FUNCTION;

	// This should only happen if close() is called before this thread is ready
	if (!ns->m_go)
	{
		log_debug("av_streamer: !ns->m_go, returning");
		return;
	}

	if (!ns->m_ACodecCtx && !ns->m_VCodecCtx && !ns->m_FormatCtx)
	{
		if (!ns->startPlayback())
		{
			log_debug("av_streamer: !ns->startPlayback, returning");
			return;
		}
	}
	else
	{
		// We need to restart the audio
		sound_handler* s = get_sound_handler();
		if (s) {
			s->attach_aux_streamer(audio_streamer, ns);
		}
	}

	ns->setStatus(playStart);

	ns->m_last_video_timestamp = 0;
	ns->m_last_audio_timestamp = 0;
	ns->m_current_timestamp = 0;

	ns->m_start_clock = tu_timer::ticks_to_seconds(tu_timer::get_ticks());

	ns->m_unqueued_data = NULL;

	boost::mutex::scoped_lock lock(ns->decode_wait_mutex);

	// Loop while we're playing
	while (ns->m_go)
	{
		log_debug("Decoding iteration");

		if (ns->m_isFLV) {
			// If queues are full then don't bother filling it
			if (ns->m_qvideo.size() < 20 || ns->m_qvideo.size() < 20) {

				// If we have problems with decoding - break
				if (!ns->decodeFLVFrame() && ns->m_start_onbuffer == false && ns->m_qvideo.size() == 0 && ns->m_qaudio.size() == 0)
				{
					break;
				}
			}

			if (ns->m_pause || (ns->m_qvideo.size() > 10 && ns->m_qaudio.size() > 10))
			{ 
				log_debug("Waiting on lock..");
				ns->decode_wait.wait(lock);
				log_debug("Finished waiting.");
			}
		} else {

			// If we have problems with decoding - break
			if (ns->decodeMediaFrame() == false && ns->m_start_onbuffer == false && ns->m_qvideo.size() == 0 && ns->m_qaudio.size() == 0)
			{
				break;
			}

			// If paused, wait for being unpaused, or
			// if the queue is full we wait until someone notifies us that data is needed.
			if (ns->m_pause || ((ns->m_qvideo.size() > 0 && ns->m_qaudio.size() > 0) && ns->m_unqueued_data))
			{ 
				log_debug("Waiting on lock..");
				ns->decode_wait.wait(lock);
				log_debug("Finished waiting.");
			}
		}

	}

	log_debug("Out of decoding loop");
	ns->m_go = false;

	log_debug("Setting playStop status");
	ns->setStatus(playStop);
}

// audio callback is running in sound handler thread
bool NetStreamFfmpeg::audio_streamer(void *owner, uint8_t *stream, int len)
{
	GNASH_REPORT_FUNCTION;

	NetStreamFfmpeg* ns = static_cast<NetStreamFfmpeg*>(owner);

	boost::mutex::scoped_lock  lock(ns->decoding_mutex);

	if (!ns->m_go || ns->m_pause) return false;

	while (len > 0 && ns->m_qaudio.size() > 0)
	{
		raw_mediadata_t* samples = ns->m_qaudio.front();

		// If less than 3 frames in the queue notify the decoding thread
		// so that we don't suddenly run out.
		if (ns->m_qaudio.size() < 3) {
			ns->decode_wait.notify_one();
		}

		int n = imin(samples->m_size, len);
		memcpy(stream, samples->m_ptr, n);
		stream += n;
		samples->m_ptr += n;
		samples->m_size -= n;
		len -= n;

		ns->m_current_timestamp = samples->m_pts;

		if (samples->m_size == 0)
		{
			ns->m_qaudio.pop();
			delete samples;
		}

	}
	return true;
}

bool NetStreamFfmpeg::decodeFLVFrame()
{
	AVPacket packet;

	FLVFrame* frame;
	if (m_qvideo.size() < m_qaudio.size()) {
		frame = m_parser->nextVideoFrame();
	} else {
		frame = m_parser->nextAudioFrame();
	}

	if (frame == NULL) {
		if (_netCon->loadCompleted()) {
			// Stop!
			m_go = false;
		} else {
			// We pause and load and buffer a second before continuing.
			pauseDecoding();
			m_bufferTime = static_cast<uint32_t>(m_current_timestamp) * 1000 + 1000;
			setStatus(bufferEmpty);
			m_start_onbuffer = true;
		}
		return false;
	}

	packet.destruct = avpacket_destruct;
	packet.size = frame->dataSize;
	packet.data = frame->data;
	// FIXME: is this the right value for packet.dts?
	packet.pts = packet.dts = static_cast<int64_t>(frame->timestamp);

	if (frame->tag == 9) {
		packet.stream_index = 0;
		return decodeVideo(&packet);
	} else {
		packet.stream_index = 1;
		return decodeAudio(&packet);
	}

}

bool NetStreamFfmpeg::decodeAudio(AVPacket* packet)
{
	if (!m_ACodecCtx) return false;

	int frame_size;
	unsigned int bufsize = (AVCODEC_MAX_AUDIO_FRAME_SIZE * 3) / 2;

	uint8_t* ptr = new uint8_t[bufsize];
#ifdef FFMPEG_AUDIO2
	frame_size = bufsize;
	if (avcodec_decode_audio2(m_ACodecCtx, (int16_t*) ptr, &frame_size, packet->data, packet->size) >= 0)
#else
	if (avcodec_decode_audio(m_ACodecCtx, (int16_t*) ptr, &frame_size, packet->data, packet->size) >= 0)
#endif
	{

		bool stereo = m_ACodecCtx->channels > 1 ? true : false;
		int samples = stereo ? frame_size >> 2 : frame_size >> 1;
		
		if (_resampler.init(m_ACodecCtx)){
			// Resampling is needed.
			
			uint8_t* output = new uint8_t[bufsize];
			
			samples = _resampler.resample(reinterpret_cast<int16_t*>(ptr), 
							 reinterpret_cast<int16_t*>(output), 
							 samples);
			delete [] ptr;
			ptr = reinterpret_cast<uint8_t*>(output);
		}
		
	  	raw_mediadata_t* raw = new raw_mediadata_t();
		
		raw->m_data = ptr;
		raw->m_ptr = raw->m_data;
		raw->m_size = samples * 2 * 2; // 2 for stereo and 2 for samplesize = 2 bytes
		raw->m_stream_index = m_audio_index;

		// set presentation timestamp
		if (packet->dts != static_cast<signed long>(AV_NOPTS_VALUE))
		{
			if (!m_isFLV) raw->m_pts = as_double(m_audio_stream->time_base) * packet->dts;
			else raw->m_pts = as_double(m_ACodecCtx->time_base) * packet->dts;
		}

		if (raw->m_pts != 0)
		{	
			// update audio clock with pts, if present
			m_last_audio_timestamp = raw->m_pts;
		}
		else
		{
			raw->m_pts = m_last_audio_timestamp;
		}

		// update video clock for next frame
		double frame_delay;
		if (!m_isFLV) frame_delay = as_double(m_audio_stream->codec->time_base);
		else frame_delay = static_cast<double>(m_parser->audioFrameDelay())/1000.0;

		m_last_audio_timestamp += frame_delay;

		if (m_isFLV) m_qaudio.push(raw);
		else m_unqueued_data = m_qaudio.push(raw) ? NULL : raw;
	}
	return true;
}

bool NetStreamFfmpeg::decodeVideo(AVPacket* packet)
{
	if (!m_VCodecCtx) return false;

	int got = 0;
	avcodec_decode_video(m_VCodecCtx, m_Frame, &got, packet->data, packet->size);
	if (got) {
		boost::scoped_array<uint8_t> buffer;

		if (m_imageframe == NULL) {
			if (m_videoFrameFormat == render::YUV) {
				m_imageframe = new image::yuv(m_VCodecCtx->width, m_VCodecCtx->height);
			} else if (m_videoFrameFormat == render::RGB) {
				m_imageframe = new image::rgb(m_VCodecCtx->width, m_VCodecCtx->height);
			}
		}

		if (m_videoFrameFormat == render::NONE) { // NullGui?
			return false;

		} else if (m_videoFrameFormat == render::YUV && m_VCodecCtx->pix_fmt != PIX_FMT_YUV420P) {
			assert(0);	// TODO
			//img_convert((AVPicture*) pFrameYUV, PIX_FMT_YUV420P, (AVPicture*) pFrame, pCodecCtx->pix_fmt, pCodecCtx->width, pCodecCtx->height);

		} else if (m_videoFrameFormat == render::RGB && m_VCodecCtx->pix_fmt != PIX_FMT_RGB24) {
			AVFrame* frameRGB = avcodec_alloc_frame();
			unsigned int numBytes = avpicture_get_size(PIX_FMT_RGB24, m_VCodecCtx->width, m_VCodecCtx->height);
			buffer.reset(new uint8_t[numBytes]);
			avpicture_fill((AVPicture *)frameRGB, buffer.get(), PIX_FMT_RGB24, m_VCodecCtx->width, m_VCodecCtx->height);
			img_convert((AVPicture*) frameRGB, PIX_FMT_RGB24, (AVPicture*) m_Frame, m_VCodecCtx->pix_fmt, m_VCodecCtx->width, m_VCodecCtx->height);
			av_free(m_Frame);
			m_Frame = frameRGB;
		}

		raw_mediadata_t* video = new raw_mediadata_t;
		if (m_videoFrameFormat == render::YUV) {
			video->m_data = new uint8_t[static_cast<image::yuv*>(m_imageframe)->size()];
		} else if (m_videoFrameFormat == render::RGB) {
			image::rgb* tmp = static_cast<image::rgb*>(m_imageframe);
			video->m_data = new uint8_t[tmp->m_pitch * tmp->m_height];
		}

		video->m_ptr = video->m_data;
		video->m_stream_index = m_video_index;
		video->m_pts = 0;

		// set presentation timestamp
		if (packet->dts != static_cast<signed long>(AV_NOPTS_VALUE))
		{
			if (!m_isFLV)	video->m_pts = as_double(m_video_stream->time_base) * packet->dts;
			else video->m_pts = as_double(m_VCodecCtx->time_base) * packet->dts;
		}

		if (video->m_pts != 0)
		{	
			// update video clock with pts, if present
			m_last_video_timestamp = video->m_pts;
		}
		else
		{
			video->m_pts = m_last_video_timestamp;
		}

		// update video clock for next frame
		double frame_delay;
		if (!m_isFLV) frame_delay = as_double(m_video_stream->codec->time_base);
		else frame_delay = static_cast<double>(m_parser->videoFrameDelay())/1000.0;

		// for MPEG2, the frame can be repeated, so we update the clock accordingly
		frame_delay += m_Frame->repeat_pict * (frame_delay * 0.5);

		m_last_video_timestamp += frame_delay;

		if (m_videoFrameFormat == render::YUV) {
			image::yuv* yuvframe = static_cast<image::yuv*>(m_imageframe);
			int copied = 0;
			uint8_t* ptr = video->m_data;
			for (int i = 0; i < 3 ; i++)
			{
				int shift = (i == 0 ? 0 : 1);
				uint8_t* yuv_factor = m_Frame->data[i];
				int h = m_VCodecCtx->height >> shift;
				int w = m_VCodecCtx->width >> shift;
				for (int j = 0; j < h; j++)
				{
					copied += w;
					assert(copied <= yuvframe->size());
					memcpy(ptr, yuv_factor, w);
					yuv_factor += m_Frame->linesize[i];
					ptr += w;
				}
			}
			video->m_size = copied;
		} else if (m_videoFrameFormat == render::RGB) {

			uint8_t* srcptr = m_Frame->data[0];
			uint8_t* srcend = m_Frame->data[0] + m_Frame->linesize[0] * m_VCodecCtx->height;
			uint8_t* dstptr = video->m_data;
			unsigned int srcwidth = m_VCodecCtx->width * 3;

			video->m_size = 0;

			while (srcptr < srcend) {
				memcpy(dstptr, srcptr, srcwidth);
				srcptr += m_Frame->linesize[0];
				dstptr += srcwidth;
				video->m_size += srcwidth;
			}

		}

		if (m_isFLV) m_qvideo.push(video);
		else m_unqueued_data = m_qvideo.push(video) ? NULL : video;

		return true;
	}

	return false;
}

bool NetStreamFfmpeg::decodeMediaFrame()
{
	boost::mutex::scoped_lock  lock(decoding_mutex);

	if (m_unqueued_data)
	{
		if (m_unqueued_data->m_stream_index == m_audio_index)
		{
			sound_handler* s = get_sound_handler();
			if (s)
			{
				m_unqueued_data = m_qaudio.push(m_unqueued_data) ? NULL : m_unqueued_data;
			}
		}
		else
		if (m_unqueued_data->m_stream_index == m_video_index)
		{
			m_unqueued_data = m_qvideo.push(m_unqueued_data) ? NULL : m_unqueued_data;
		}
		else
		{
			log_error(_("read_frame: not audio & video stream"));
		}
		return true;
	}

	AVPacket packet;
	int rc = av_read_frame(m_FormatCtx, &packet);

	if (rc >= 0)
	{
		if (packet.stream_index == m_audio_index && get_sound_handler())
		{
			if (!decodeAudio(&packet)) {
				log_error(_("Problems decoding audio frame"));
				return false;
			}
		}
		else
		if (packet.stream_index == m_video_index)
		{
			if (!decodeVideo(&packet)) {
				log_error(_("Problems decoding video frame"));
				return false;
			}
		}
		av_free_packet(&packet);
	}
	else
	{
		log_error(_("Problems decoding frame"));
		return false;
	}

	return true;
}

void
NetStreamFfmpeg::seek(double pos)
{
	boost::mutex::scoped_lock  lock(decoding_mutex);

	long newpos = 0;
	double timebase = 0;

	// Seek to new position
	if (m_isFLV) {
		if (m_parser.get()) {
			newpos = m_parser->seek(static_cast<uint32_t>(pos*1000));
		} else {
			newpos = 0;
		}
	} else if (m_FormatCtx) {

		timebase = static_cast<double>(m_FormatCtx->streams[m_video_index]->time_base.num) / static_cast<double>(m_FormatCtx->streams[m_video_index]->time_base.den);
		newpos = static_cast<long>(pos / timebase);
		
		if (av_seek_frame(m_FormatCtx, m_video_index, newpos, 0) < 0) {
			log_error(_("%s: seeking failed"), __FUNCTION__);
			return;
		}
	} else {
		return;
	}

	// This is kindof hackish and ugly :-(
	if (newpos == 0) {
		m_last_video_timestamp = 0;
		m_last_audio_timestamp = 0;
		m_current_timestamp = 0;

		m_start_clock = tu_timer::ticks_to_seconds(tu_timer::get_ticks());

	} else if (m_isFLV) {
		double newtime = static_cast<double>(newpos) / 1000.0;
		m_start_clock += (m_last_audio_timestamp - newtime) / 1000.0;

		m_last_audio_timestamp = newtime;
		m_last_video_timestamp = newtime;
		m_current_timestamp = newtime;
	} else {
		AVPacket Packet;
		av_init_packet(&Packet);
		double newtime = 0;
		while (newtime == 0) {
			if ( av_read_frame(m_FormatCtx, &Packet) < 0) {
				av_seek_frame(m_FormatCtx, -1, 0,AVSEEK_FLAG_BACKWARD);
				av_free_packet(&Packet);
				return;
			}

			newtime = timebase * (double)m_FormatCtx->streams[m_video_index]->cur_dts;
		}

		av_free_packet(&Packet);
		av_seek_frame(m_FormatCtx, m_video_index, newpos, 0);

		m_start_clock += (m_last_audio_timestamp - newtime) / 1000.0;

		m_last_audio_timestamp = newtime;
		m_last_video_timestamp = newtime;
		m_current_timestamp = newtime;
	}
	
	// Flush the queues
	while (m_qvideo.size() > 0)
	{
		delete m_qvideo.front();
		m_qvideo.pop();
	}

	while (m_qaudio.size() > 0)
	{
		delete m_qaudio.front();
		m_qaudio.pop();
	}

}

void
NetStreamFfmpeg::refreshVideoFrame()
{
	// If we're paused or not running, there is no need to do this
	if (!m_go || m_pause) return;

	// Loop until a good frame is found
	while(1) {
		// Get video frame from queue, will have the lowest timestamp
		raw_mediadata_t* video = m_qvideo.front();

		// If the queue is empty, we tell the decoding thread to wake up,
		// and decode some more.
		if (!video) {
			decode_wait.notify_one();
			return;
		}

		// Caclulate the current time
		double current_clock;
		if (m_ACodecCtx && get_sound_handler()) {
			current_clock = m_current_timestamp;
		} else {
			current_clock = (tu_timer::ticks_to_seconds(tu_timer::get_ticks()) - m_start_clock)*1000;
			m_current_timestamp = current_clock;
		}

		double video_clock = video->m_pts;

		// If the timestamp on the videoframe is smaller than the
		// current time, we put it in the output image.
		if (current_clock >= video_clock)
		{
			boost::mutex::scoped_lock lock(image_mutex);
			if (m_videoFrameFormat == render::YUV) {
				// XXX m_imageframe might be a byte aligned buffer, while video is not!
				static_cast<image::yuv*>(m_imageframe)->update(video->m_data);
			} else if (m_videoFrameFormat == render::RGB) {

				image::rgb* imgframe = static_cast<image::rgb*>(m_imageframe);
				rgbcopy(imgframe, video, m_VCodecCtx->width * 3);
			}

			// Delete the frame from the queue
			m_qvideo.pop();
			delete video;

			// A frame is ready for pickup
			m_newFrameReady = true;

		} else {
			// The timestamp on the first frame in the queue is greater
			// than the current time, so no need to do anything.
			return;
		}

		// If less than 3 frames in the queue notify the decoding thread
		// so that we don't suddenly run out.
		if (m_qvideo.size() < 3) {
			decode_wait.notify_one();
		}
	}
}


void
NetStreamFfmpeg::advance()
{
	// This can happen in 2 cases: 
	// 1) When playback has just started and we've been waiting for the buffer 
	//    to be filled (buffersize set by setBufferTime() and default is 100
	//    miliseconds).
	// 2) The buffer has be "starved" (not being filled as quickly as needed),
	//    and we then wait until the buffer contains some data (1 sec) again.
	if (m_go && m_pause && m_start_onbuffer && m_parser.get() && m_parser->isTimeLoaded(m_bufferTime)) {
		setStatus(bufferFull);
		unpauseDecoding();
		m_start_onbuffer = false;
	}

	// Check if there are any new status messages, and if we should
	// pass them to a event handler
	processStatusNotifications();

	// Find video frame with the most suited timestamp in the video queue,
	// and put it in the output image frame.
	refreshVideoFrame();
}

int64_t
NetStreamFfmpeg::time()
{

	if (m_FormatCtx && m_FormatCtx->nb_streams > 0) {
		double time = (double)m_FormatCtx->streams[0]->time_base.num / (double)m_FormatCtx->streams[0]->time_base.den * (double)m_FormatCtx->streams[0]->cur_dts;
		return static_cast<int64_t>(time);
	} else if (m_isFLV) {
		return static_cast<int64_t>(m_current_timestamp);
	} else {
		return 0;
	}
}

void NetStreamFfmpeg::pauseDecoding()
{
	if (m_pause) return;

	m_pause = true;

	// Save the current time so we later can tell how long the pause lasted
	m_time_of_pause = tu_timer::ticks_to_seconds(tu_timer::get_ticks());
}

void NetStreamFfmpeg::unpauseDecoding()
{
	if (!m_pause) return;

	m_pause = false;	

	if (m_current_timestamp == 0) {
		m_start_clock = tu_timer::ticks_to_seconds(tu_timer::get_ticks());
	} else {
		// Add the paused time to the start time so that the playhead doesn't
		// noticed that we have been paused
		m_start_clock += tu_timer::ticks_to_seconds(tu_timer::get_ticks()) - m_time_of_pause;
	}

	// Notify the decode thread/loop that we are running again
	decode_wait.notify_one();

	// Re-connect to the soundhandler
	sound_handler* s = get_sound_handler();
	if (s) s->attach_aux_streamer(audio_streamer, (void*) this);
}


} // gnash namespcae

#endif // USE_FFMPEG
