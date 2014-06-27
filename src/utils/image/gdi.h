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

#ifndef IMAGE_GDI_H
#define IMAGE_GDI_H

#include <windows.h>
#include <GdiPlus.h>
#include "image.h"

//fmtOut gets the closest possible GDI+ format to be written to the image
Image::PixelFormat GdiToImagePixelFormat(Gdiplus::PixelFormat fmt, Gdiplus::PixelFormat& fmtOut);
Gdiplus::PixelFormat ImageGdiPixelFormat(Image::PixelFormat fmt);
void GdiBitmapToImage(Gdiplus::Bitmap& bit, Image& img);
Gdiplus::Bitmap* ImageToGdiBitmap(const Image& img);

bool ExportImageGDIPlus(const char* fn, Image& img);

#endif
