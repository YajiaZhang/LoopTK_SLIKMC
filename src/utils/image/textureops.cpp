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

#include "textureops.h"


bool LoadImageFromBitmap(HDC hdc, HBITMAP hbit, Image& img)
{
	BITMAP bit;

	if (!GetObject(hbit, sizeof(BITMAP), (LPSTR)&bit)) 
	{
		ReportError("couldnt get object");
		return false;
	}

	switch(bit.bmBitsPixel)
	{
	case 8:
		img.format = Image::A8;
		break;
	case 16: //why not R5G6B5?
		img.format = Image::R8G8B8;
		break;
	case 24: //why not A8R8G8B8?
	case 32:
		img.format = Image::R8G8B8;
		break;
	default:
		ReportError("Unknown bitmap color depth, %d", bit.bmBitsPixel);
		return NULL;
		break;
	}

	img.initialize(bit.bmWidth,bit.bmHeight,img.format);

	BITMAPINFO bi;
	ZeroMemory(&bi, sizeof(bi));
	bi.bmiHeader.biSize = sizeof(BITMAPINFOHEADER);
	bi.bmiHeader.biWidth = img.w;
	bi.bmiHeader.biHeight = img.h;
	bi.bmiHeader.biBitCount = img.pixelBPP();
	bi.bmiHeader.biPlanes = 1;
	bi.bmiHeader.biSizeImage = 0;
	bi.bmiHeader.biCompression = BI_RGB;

	if( 0 == GetDIBits(hdc, hbit, 0, img.h, img.data, &bi, DIB_RGB_COLORS))
	{
		ReportError("couldnt get bits!");
		return false;
	}

	return true;
}




int sumbits(unsigned int word)
{
	int count=0;
	for(int i=0; i<32; i++)
	{
		if(word & 0x1)
			count++;
		word = word >> 1;
	}
	return count;
}

int highestbit(unsigned int word)
{
	int mask = 0x80000000;
	for(int i=32; i>0; i--)
	{
		if(word & mask)
			return i;
		mask = mask >> 1;
	}
	return 0;
}

bool IsSquare(const Image& tex)
{
	if(sumbits(tex.w) != 1)
		return false;
	if(sumbits(tex.h) != 1)
		return false;
	return true;
}

void makeSquareDimensions(Image& tex)
{
	unsigned int newwidth = tex.w, newheight = tex.h;
	if(sumbits(tex.w) != 1)
		newwidth = 1 << highestbit(tex.w);
	if(sumbits(tex.h) != 1)
		newheight = 1 << highestbit(tex.h);

	tex.initialize(newwidth,newheight,tex.format);
}


void ExpandSquare(Image& tex)
{
	if(IsSquare(tex))
		return;
	Image old = tex;

	makeSquareDimensions(tex);
	tex.clear();
	old.blit(tex);
}


void StretchSquare(Image& tex)
{
	if(IsSquare(tex))
		return;
	ImageOperator oldtex = tex;
	makeSquareDimensions(tex);
	ImageOperator newtex;

	newtex.initialize(tex.w, tex.h);
	oldtex.stretchBlit(newtex);

	newtex.output(tex, tex.format);
}
