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

#include "image.h"
#include <assert.h>

bool ImportImageBMP(const char* fn, Image& image)
{
	//todo:
	//RLE compression decoding?

	FILE* f = fopen(fn, "rb");
	if(!f)
		return false;

	BITMAPFILEHEADER bfh;
	if(fread(&bfh, sizeof(bfh), 1, f) != 1)
	{
		ReportError("Couldn't load file header");
		fclose(f);
		return false;
	}

	if(memcmp(&bfh.bfType, "BM", 2) != 0) 
	{
		ReportError("This isn't a bitmap file");
		return false;
	}

	BITMAPINFOHEADER bi;
	if(fread(&bi, sizeof(bi), 1, f) != 1)
	{
		ReportError("Couldn't load info header");
		fclose(f);
		return false;
	}

	image.w = bi.biWidth;
	image.h = bi.biHeight;
	image.num_bytes = bi.biSizeImage;


	switch(bi.biCompression)
	{
	case BI_RGB:
		if(image.num_bytes == 0)
			image.num_bytes = image.w * image.h * (bi.biBitCount>>3);
		break;
	default:
		ReportError("Unsupported compression type");
		fclose(f);
		return false;
	}

	int size_to_read = image.num_bytes;

	RGBQUAD palette [256];
	int palette_size = 256;
	bool decode_8 = false;

	switch(bi.biBitCount)
	{
	case 8:
		{
		if(bi.biClrUsed != 0)
			palette_size = bi.biClrUsed;
		decode_8 = true;
		if(fread(&palette, sizeof(RGBQUAD), palette_size, f) != palette_size)
		{
			ReportError("Couldnt read palette");
			fclose(f);
			return false;
		}
		image.format = Image::R8G8B8;
		image.num_bytes = image.w * image.h * 3;
		}
		break;
	case 16:
		image.format = Image::R5G6B5;
		break;
	case 24:
		image.format = Image::R8G8B8;
		break;
	default:
		ReportError("Unsupported bit count");
		fclose(f);
		return false;
	}
	
	fseek(f, bfh.bfOffBits, SEEK_SET);

	unsigned char* bits = new unsigned char [size_to_read];
	if(fread(bits, 1, size_to_read, f) != size_to_read)
	{
		ReportError("Couldnt read bits\n");
		fclose(f);
		delete [] bits;
		return false;
	}

	if(ftell(f) != bfh.bfSize)
	{
		int cur = ftell(f);
		fseek(f, 0, SEEK_END);
		ReportError("Um, there was some stuff missed: offset is %d, size is %d, to read is %d, current is %d, end is %d\n",
			bfh.bfOffBits, bfh.bfSize, size_to_read, cur, ftell(f));
	}

	assert(image.num_bytes == image.w*image.h*image.pixelSize());
	image.initialize(image.w,image.h,image.format);
	if(decode_8)
	{
		unsigned char* bit = bits;
		unsigned char* pixel = image.data;
		for(int i=0; i<size_to_read; i++)
		{
			if(*bit >= palette_size)
			{
				ReportError("out of palette range\n");
				abort();
			}
			pixel[0] = palette[*bit].rgbBlue;
			pixel[1] = palette[*bit].rgbGreen;
			pixel[2] = palette[*bit].rgbRed;
			pixel+=3;
			bit++;
		}
	}
	else
	{
		assert(image.pixelBPP() == bi.biBitCount);
		memcpy(image.data, bits, image.num_bytes);
	}


	delete [] bits;
	fclose(f);
	return true;
}
