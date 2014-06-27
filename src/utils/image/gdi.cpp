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

#include "gdi.h"
#include <utils/stringutils.h>

Image::PixelFormat GdiToImagePixelFormat(Gdiplus::PixelFormat fmt, Gdiplus::PixelFormat& fmtOut)
{
	fmtOut = fmt;  //unless we need a conversion
	switch(fmt)
	{
	case PixelFormat1bppIndexed:
	case PixelFormat4bppIndexed:
	case PixelFormat8bppIndexed:
		fmtOut = PixelFormat32bppARGB;
		return Image::A8R8G8B8;
	case PixelFormat16bppARGB1555:
		return Image::X1R5G5B5;
	case PixelFormat16bppGrayScale:
		fmtOut = PixelFormat32bppARGB;
		return Image::A8R8G8B8;
	case PixelFormat16bppRGB555:
		return Image::X1R5G5B5;
	case PixelFormat16bppRGB565:
		return Image::R5G6B5;
	case PixelFormat24bppRGB:
		return Image::R8G8B8;
	case PixelFormat32bppARGB:
	case PixelFormat32bppPARGB: 
	case PixelFormat48bppRGB:
	case PixelFormat32bppRGB:
	case PixelFormat64bppARGB: 
	case PixelFormat64bppPARGB:
		fmtOut = PixelFormat32bppARGB;
		return Image::A8R8G8B8;
	default:
		fmtOut = PixelFormat32bppARGB;
		return Image::A8R8G8B8;
	}
}

Gdiplus::PixelFormat ImageToGdiPixelFormat(Image::PixelFormat fmt)
{
	switch(fmt) {
	case Image::None: return PixelFormat1bppIndexed;	
	case Image::R8G8B8: return PixelFormat24bppRGB;
	case Image::A8R8G8B8: return PixelFormat32bppARGB;
	case Image::R5G6B5: return PixelFormat16bppRGB565;
	case Image::X1R5G5B5: return PixelFormat16bppRGB555;
	case Image::A8: return PixelFormat8bppIndexed;
	default:
		ReportError("Thats an invalid GDI format");
		return PixelFormat1bppIndexed;	
	}
}

void GdiBitmapToImage(Gdiplus::Bitmap& bit, Image& img)
{
	Gdiplus::PixelFormat bmpFormat;
	Image::PixelFormat imageFormat;
	imageFormat=GdiToImagePixelFormat(bit.GetPixelFormat(),bmpFormat);

	Gdiplus::Rect rect;
	Gdiplus::RectF rectf;
	Gdiplus::Unit units;
	Gdiplus::BitmapData bitdata;
	bit.GetBounds(&rectf,&units);
	rect.X = (int)rectf.X;
	rect.Y = (int)rectf.Y;
	rect.Width = (int)rectf.Width;
	rect.Height = (int)rectf.Height;
	bit.LockBits(&rect,Gdiplus::ImageLockModeRead,bmpFormat,&bitdata);

	img.initialize(bitdata.Width,bitdata.Height,imageFormat);
	unsigned char* data = (unsigned char*)bitdata.Scan0;
	int scanline_size = bitdata.Width*img.pixelSize();
	for(UINT i=0; i<bitdata.Height;i++)
	{
		memcpy(img.getData(0,i),data,scanline_size);
		data += bitdata.Stride;
	}
	bit.UnlockBits(&bitdata);
}

Gdiplus::Bitmap* ImageToGdiBitmap(const Image& image)
{
	Gdiplus::PixelFormat fmt=ImageToGdiPixelFormat(image.format);
	return new Gdiplus::Bitmap(image.w,image.h,image.pixelSize()*image.w,fmt,image.data);
}

bool ImportImageGDIPlus(const char* fn, Image& img)
{
	//load the bitmap (requires wchars)
	int length = strlen(fn)+5;
	unsigned short* wfn = new unsigned short[length];
	ToWideChar(fn,wfn,length);
	Gdiplus::Bitmap bit(wfn);
	delete [] wfn;

	if(bit.GetWidth() == 0 || bit.GetHeight() == 0) return false;
	GdiBitmapToImage(bit,img);
	return true;
}

int GetEncoderClsid(const WCHAR* format, CLSID* pClsid)
{
   UINT  num = 0;          // number of image encoders
   UINT  size = 0;         // size of the image encoder array in bytes

   Gdiplus::ImageCodecInfo* pImageCodecInfo = NULL;

   Gdiplus::GetImageEncodersSize(&num, &size);
   if(size == 0)
      return -1;  // Failure

   pImageCodecInfo = (Gdiplus::ImageCodecInfo*)(malloc(size));
   if(pImageCodecInfo == NULL)
      return -1;  // Failure

   Gdiplus::GetImageEncoders(num, size, pImageCodecInfo);

   for(UINT j = 0; j < num; ++j)
   {
      if( wcscmp(pImageCodecInfo[j].MimeType, format) == 0 )
      {
         *pClsid = pImageCodecInfo[j].Clsid;
         free(pImageCodecInfo);
         return j;  // Success
      }    
   }

   free(pImageCodecInfo);
   return -1;  // Failure
}

bool ExportImageGDIPlus(const char* fn, Image& img)
{
	//get the proper encoder
	const char* ext=FileExtension(fn);
	if(!ext || strlen(ext) > 4) ext="bmp";
	if(0==strcmp(ext,"tif")) ext="tiff";
	char typebuf[32];
	_snprintf(typebuf,32,"image/%s", ext);
	WCHAR wtypebuf[32];
	ToWideChar(typebuf,wtypebuf,32);
	CLSID clsid;
	if(GetEncoderClsid(wtypebuf, &clsid)==-1)
		return false;

	Gdiplus::Bitmap* bmp=ImageToGdiBitmap(img);

	//save the bitmap (requires wchars)
	int length = strlen(fn)+5;
	unsigned short* wfn = new unsigned short[length];
	ToWideChar(fn,wfn,length);
	bmp->Save(wfn, &clsid, NULL);
	delete [] wfn;
	delete bmp;
	return true;
}
