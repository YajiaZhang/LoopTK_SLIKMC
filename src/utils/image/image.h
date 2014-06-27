/*
    LoopTK: Protein Loop Kinematic Toolkit
    Copyright (C) 2007 Stanford University

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/

#ifndef IMAGE_IMAGE_H
#define IMAGE_IMAGE_H

#include "myfile.h"

class Image
{
public:
	Image();
	Image(const Image& img);
	virtual ~Image();
	
	const Image& operator = (const Image&);

	enum PixelFormat {
		None,
		R8G8B8,				//24 bit, 1 byte per channel
		A8R8G8B8,			//32 bit, 1 byte per channel
		R5G6B5,				//16 bit, green has 6 bits
		X1R5G5B5,			//15 bit, 1 padding
		A8,					//8 bits of alpha
		FloatRGB,			//floating point rgb
		FloatRGBA,			//floating point rgba
		FloatA,				//floating point alpha
	};
	static bool isValidFormat(PixelFormat);
	static unsigned int pixelFormatSize(PixelFormat);

	//returns -1 on failure, returns size of image on success
	virtual int initialize(int w, int h, PixelFormat format);
	virtual void unload();

	bool Read(const char*);
	bool Write(const char*) const;
	virtual bool Read(File&);
	virtual bool Write(File&) const;

	//fills in all pixels with the given color
	void clear(int dat = 0);
	//copies this image, with rect (sx,sy,w,h), onto the dest at dx, dy
	void blit(Image& dest, int sx = 0, int sy = 0, int w = -1, int h = -1, int dx = 0, int dy = 0) const;

	inline unsigned int pixelSize() const { return pixelFormatSize(format); }  //size in bytes
	inline unsigned int pixelBPP() const { return pixelSize()<<3; }
	inline unsigned int pitch() const;
	unsigned char* getData(int x, int y) const;

	PixelFormat format;			//format of image
	unsigned short w,h;			//width and height, in pixels

	//the image data
	unsigned char* data;
	unsigned int num_bytes;
};

class ImageMipmapped : public Image
{
public:
	ImageMipmapped();
	ImageMipmapped(const Image&);
	virtual ~ImageMipmapped();

	const ImageMipmapped& operator = (const Image&);
	const ImageMipmapped& operator = (const ImageMipmapped&);

	void createMipmaps();			//creates mipmaps from the first texture level
	virtual void unload();

	virtual bool Read(File&);

	//mipmap images, starting at 1 after the main image
	unsigned char** mipmap_data;
	unsigned int num_mipmap_levels;
};

typedef unsigned char color_channel_t;

struct Pixel {
	void add(const Pixel& a, const Pixel& b);
	void sub(const Pixel& a, const Pixel& b);
	void mul(const Pixel& a, const Pixel& b);
	void mix(const Pixel& a, const Pixel& b, float u);
	color_channel_t r,g,b,a;
};

class ImageOperator
{
public:
	ImageOperator();
	ImageOperator(const Image&);
	~ImageOperator();

	void initialize(int w, int h);
	void resize(int w, int h);
	void unload();
	void clear();
	void clear(const Pixel&);

	void operator = (const Image&);
	void operator = (const ImageOperator&);
	void output(Image&, Image::PixelFormat fmt) const;

	void stretchBlit(ImageOperator& dest) const;

	Pixel& getPixel(int x, int y);
	const Pixel& getPixel(int x, int y) const;
	void sampleNearest(float x, float y, Pixel& out) const;
	void sampleLinear(float x, float y, Pixel& out) const;

	Pixel* pixels;
	unsigned short w,h;
	unsigned int num_pixels;
};

#endif
