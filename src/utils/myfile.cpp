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

#include "myfile.h"
#include "utils.h"
#include <assert.h>

//platform-specific defines

#ifdef WIN32

inline FILE_POINTER FileOpen(const char* fn, int openmode)
{
	if(openmode & FILEREAD)
	{
		if(openmode & FILEWRITE)
		        return CreateFile(fn, GENERIC_READ | GENERIC_WRITE, 0, NULL, OPEN_ALWAYS, FILE_ATTRIBUTE_NORMAL, NULL);
		else
		        return CreateFile(fn, GENERIC_READ, FILE_SHARE_READ, NULL, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, NULL);
	}
	else if(openmode & FILEWRITE)
	{
		return CreateFile(fn, GENERIC_WRITE, 0, NULL, CREATE_ALWAYS, FILE_ATTRIBUTE_NORMAL, NULL);
	}
	abort();
}

inline void FileClose(FILE_POINTER x)
{
  CloseHandle(x);
}

inline int FilePosition(FILE_POINTER x)
{
  return SetFilePointer(x, 0, 0, FILE_CURRENT);
}

inline bool FileSeek(FILE_POINTER x, int p, int from)
{
  return (SetFilePointer(x, p, 0, from) != 0);
}

inline int FileLength(FILE_POINTER x) {
  return GetFileSize(x, NULL);
}

inline bool FileRead(FILE_POINTER x, void* d, int size)
{
  DWORD s;
  if(!ReadFile(x, d, size, &s, NULL))
    return false;
  return (size == s);
}

inline bool FileWrite(FILE_POINTER x, const void* d, int size)
{
  DWORD s;
  if(!WriteFile(x, d, size, &s, NULL))
    return false;
  return (size == s);
}

#else

inline FILE_POINTER FileOpen(const char* fn, int openmode) {
  if(openmode & FILEREAD) {
    if(openmode & FILEWRITE)
      return fopen(fn,"r+b");
    else
      return fopen(fn,"rb");
  }
  else {
    return fopen(fn,"wb");
  }
}

inline void FileClose(FILE_POINTER x)
{
  fclose(x);
}

inline int FilePosition(FILE_POINTER x)
{
  return ftell(x);
}

inline bool FileSeek(FILE_POINTER x, int p, int from)
{
  return fseek(x,p,from);
}

inline int FileLength(FILE_POINTER x)
{
  long pos=ftell(x);
  fseek(x,0,SEEK_END);
  long len=ftell(x);
  fseek(x,pos,SEEK_SET);
  return len;
}

inline bool FileRead(FILE_POINTER x, void* d, int size)
{
  return (int)fread(d,1,size,x)==size;
}

inline bool FileWrite(FILE_POINTER x, const void* d, int size)
{
  return (int)fwrite(d,1,size,x)==size;
}

#endif




enum { MODE_NONE,
	MODE_MYFILE, //internally managed file
	MODE_EXTFILE, //external file reference
	MODE_MYDATA, //internally managed data
	MODE_EXTDATA //external data reference
};

File::File()
:mode(0),srctype(MODE_NONE),
file(INVALID_FILE_POINTER),
datafile(NULL),datapos(0),datasize(0)
{
}

File::~File()
{
	Close();
}

void File::Close()
{
        if(srctype == MODE_MYFILE && file != INVALID_FILE_POINTER) FileClose(file);
	if(srctype == MODE_MYDATA && datafile != NULL) free(datafile);

	srctype = MODE_NONE;
	mode = 0;
	file = INVALID_FILE_POINTER;
	datafile = NULL;
	datapos = 0;
	datasize = 0;
}

bool File::OpenData(void* data, int size, int openmode)
{
	Close();

	if(!data)
		return false;
	if(size < 0)
		return false;

	srctype = MODE_EXTDATA;
	if(openmode == 0)
		return false;

	datafile = (unsigned char*)data;
	datapos = 0;
	datasize = size;
	mode = openmode;
	return true;
}

bool File::OpenData(int openmode)
{
	Close();

	srctype = MODE_MYDATA;
	if(openmode == 0)
		return false;

	ResizeDataBuffer(64);
	mode = openmode;
	return true;
}

void* File::GetDataBuffer() const
{
	assert(srctype == MODE_MYDATA || srctype == MODE_EXTDATA);
	return datafile;
}

void File::ResizeDataBuffer(int size)
{
	unsigned char* olddata=datafile;
	datafile=(unsigned char*)malloc(size);
	memcpy(datafile,olddata,datasize);
	free(olddata);
	datasize = size;
}

bool File::Open(const char* fn, int openmode)
{
	Close();

	if(openmode == 0)
		return false;

	file=FileOpen(fn,openmode);

	if(file == INVALID_FILE_POINTER)
		return false;
	srctype = MODE_MYFILE;
	mode = openmode;
	return true;
}


bool File::Open(FILE_POINTER f, int openmode)
{
	Close();

	srctype = MODE_EXTFILE;
	if(openmode == 0)
		return false;

	file = f;
	mode = openmode;
	return true;
}

int File::Position() const
{
	switch(srctype)
	{
	case MODE_MYFILE:
	case MODE_EXTFILE:
	        return FilePosition(file);
	case MODE_MYDATA:
	case MODE_EXTDATA:
		return datapos;
	}
	return -1;
}

bool File::Seek(int p, int from)
{
	switch(srctype)
	{
	case MODE_MYFILE:
	case MODE_EXTFILE:
	       return FileSeek(file,p,from);
		break;
	case MODE_MYDATA:
	case MODE_EXTDATA:
		switch (from)
		{
		case SEEK_CUR:
			if(datapos + p >= datasize || datapos + p < 0)
				return false;
			datapos += p;
			break;
		case SEEK_SET:
			if(p >= datasize || p < 0)
				return false;
			datapos = p;
			break;
		case SEEK_END:
			if(datasize + p < 0 || p > 0)
				return false;
			datapos = datasize + p;
			break;
		}
	}
	return true;
}

int File::Length()
{
	switch(srctype)
	{
	case MODE_MYFILE:
	case MODE_EXTFILE:
	        return FileLength(file);
	case MODE_MYDATA:
	case MODE_EXTDATA:
		return datasize;
	}
	return -1;
}

bool File::ReadData(void* d, int size)
{
	if(mode & FILEREAD)
	{
		switch(srctype)
		{
		case MODE_MYFILE:
		case MODE_EXTFILE:
			return FileRead(file,d,size);
		case MODE_MYDATA:
		case MODE_EXTDATA:
			if(datapos + size > datasize)
				return false;
			memcpy(d, datafile+datapos, size);
			datapos += size;
			return true;
		}
	}
	return false;
}

bool File::WriteData(const void* d, int size)
{
	if(mode & FILEWRITE)
	{
		switch(srctype)
		{
		case MODE_MYFILE:
		case MODE_EXTFILE:
			return FileWrite(file,d,size);
		case MODE_MYDATA:		//resize if buffer's not big enough
			if(datapos + size > datasize) {
				int a=datapos+size,b=datasize*2;
				ResizeDataBuffer(Max(a,b));
			}
			memcpy(datafile+datapos, d, size);
			datapos += size;
			return true;
		case MODE_EXTDATA:
			if(datapos + size > datasize)
				return false;
			memcpy(datafile+datapos, d, size);
			datapos += size;
			return true;
		}
	}
	return false;
}



int ReadChar(FILE_POINTER file)
{
	char c;
	if(!FileRead(file, &c, 1)) return EOF;
	return c;
}

bool File::ReadString(char* str, int bufsize)
{
	if(mode & FILEREAD)
	{
		int i,c;
		switch(srctype)
		{
		case MODE_MYFILE:
		case MODE_EXTFILE:
			for(i=0; i<bufsize; i++)
			{
				c = ReadChar(file);
				if(c==EOF)
					return false;
				str[i]=c;
				if(c==0)
					return true;
			}
			break;
		case MODE_MYDATA:
		case MODE_EXTDATA:
			for(i=0; i<bufsize; i++)
			{
				if(datapos >= datasize)
					return false;
				str[i]=datafile[datapos];
				datapos++;
				if(str[i]==0)
					return true;
			}
			break;
		}
	}
	return false;
}

bool File::WriteString(const char* str)
{
	return WriteData(str, (int)(strlen(str)+1));
}



_DEFINE_READ_WRITE_FILE_BASIC(bool);
_DEFINE_READ_WRITE_FILE_BASIC(char);
_DEFINE_READ_WRITE_FILE_BASIC(signed char);
_DEFINE_READ_WRITE_FILE_BASIC(unsigned char);
_DEFINE_READ_WRITE_FILE_BASIC(short);
_DEFINE_READ_WRITE_FILE_BASIC(unsigned short);
_DEFINE_READ_WRITE_FILE_BASIC(int);
_DEFINE_READ_WRITE_FILE_BASIC(unsigned int);
_DEFINE_READ_WRITE_FILE_BASIC(long);
_DEFINE_READ_WRITE_FILE_BASIC(unsigned long);
_DEFINE_READ_WRITE_FILE_BASIC(long long);
_DEFINE_READ_WRITE_FILE_BASIC(unsigned long long);
_DEFINE_READ_WRITE_FILE_BASIC(float);
_DEFINE_READ_WRITE_FILE_BASIC(double);
_DEFINE_READ_WRITE_FILE_BASIC(long double);
