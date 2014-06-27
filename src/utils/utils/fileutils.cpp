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

#include "fileutils.h"
#ifdef WIN32
#define WIN32_LEAN_AND_MEAN
#include "windows.h"
#else
#include "stdlib.h"
#include "unistd.h"
#include "stdio.h"
#include "string.h"
#endif

namespace FileUtils {

bool Exists(const char* fn)
{
#ifdef WIN32
	HANDLE f = CreateFile((LPCTSTR)fn,GENERIC_READ,FILE_SHARE_READ,NULL,OPEN_ALWAYS,0,0);
	if(f == INVALID_HANDLE_VALUE) return false;
	CloseFile(f);
	return true;
#else
	FILE* f = fopen(fn,"r");
	if(!f) return false;
	fclose(f);
	return true;
#endif
}

bool Delete(const char* fn)
{
#ifdef WIN32
	return DeleteFile((LPCTSTR)fn) != FALSE;
#else
	return (unlink(fn)==0);
#endif
}

bool Rename(const char* from,const char* to)
{
#ifdef WIN32
	return MoveFile(from,to) != FALSE;
#else
	return (rename(from,to)==0);
	/*
	size_t len = strlen(from) + strlen(to) + 5;
	char* buf = new char[len];
	sprintf(buf,"mv %s %s",from,to);
	int res = system(buf);
	delete [] buf;
	return res == 0;
	*/
#endif
}

bool Copy(const char* from,const char* to,bool override)
{
#ifdef WIN32
	return CopyFile(from,to,(override?FALSE:TRUE)) != FALSE;
#else
	size_t len = strlen(from) + strlen(to) + 5;
	char* buf = new char[len];
	sprintf(buf,"cp %s %s %s",(override?"-f":""),from,to);
	int res = system(buf);
	delete [] buf;
	return res == 0;
#endif
}

bool TempName(char* out,const char* directory,const char* prefix)
{
#ifdef WIN32
	if(directory == NULL) directory = ".";
	UINT res = GetTempFileName(directory,prefix,0,out);
	return (res!=0);
#else
	char* fn=tempnam(directory,prefix);
	if(!fn) return false;
	strcpy(out,fn);
	return true;
#endif
}

} // namespace FileUtils
