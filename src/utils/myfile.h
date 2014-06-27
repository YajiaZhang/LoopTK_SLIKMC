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

#ifndef BASIC_FILE_H
#define BASIC_FILE_H

#ifdef WIN32
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
typedef HANDLE FILE_POINTER;
#define INVALID_FILE_POINTER INVALID_HANDLE_VALUE

#else

#include <stdlib.h>
#include <stdio.h>
typedef FILE *FILE_POINTER;
#define INVALID_FILE_POINTER NULL
#endif // WIN32

#define FILEREAD 0x1
#define FILEWRITE 0x2
#define FILESEEKSTART 0
#define FILESEEKCURRENT 1
#define FILESEEKEND 2

/** \ingroup Standard
 * \brief A cross-platform class for reading/writing binary data.
 */
class File
{
public:
	File();
	~File();
	bool Open(const char*, int openmode = FILEREAD | FILEWRITE);
	bool Open(FILE_POINTER, int openmode = FILEREAD | FILEWRITE);
	bool OpenData(void* buf, int size, int openmode = FILEREAD | FILEWRITE);
	bool OpenData(int openmode = FILEREAD | FILEWRITE);
	void Close();

	bool Seek(int, int from = FILESEEKCURRENT);
	int Position() const;
	int Length();

	bool ReadData(void*, int size);
	bool WriteData(const void*, int size);

	bool ReadString(char*, int bufsize);
	bool WriteString(const char*);

	void* GetDataBuffer() const;
	void ResizeDataBuffer(int size);

private:
	int mode;		//file read/write mode
	int srctype;	//data source mode (file,data,etc)

	FILE_POINTER file;
	unsigned char* datafile;
	int datapos;
	int datasize;
};

/** \file myfile.h
 * \ingroup Standard
 * \brief A unified interface for reading/writing binary data to file.
 *
 * A consistent interface for reading/writing data is given by the
 * template methods ReadFile() and WriteFile().  A class can be read/written
 * to file in two ways: 1) Read/WriteFile() can be overloaded
 * to read/write the data directly, or 2) use the default Read/WriteFile(), 
 * but define the following methods in the class:
 * \code
 *   bool Read(File&)
 *   bool Write(File&) const
 * \endcode
 *
 * For structs of contiguous data, the macro 
 * _DEFINE_READ_WRITE_FILE_TEMPLATE_BASIC(T) will overload the
 * the Read/WriteFile functions to automatically read 
 * the data from file contiguously for objects of type T.
 */

template <class type>
bool ReadFile(File& f, type& t) { return t.Read(f); }
template <class type>
bool WriteFile(File& f, const type& t) { return t.Write(f); }

template <class type>
bool ReadArrayFile(File& f, type* t, int n)
{
  for(int i=0;i<n;i++)
    if(!ReadFile(f,t[i])) return false;
  return true;
}

template <class type>
bool WriteArrayFile(File& f, const type* t, int n)
{
  for(int i=0;i<n;i++)
    if(!WriteFile(f,t[i])) return false;
  return true;
}

#define _DECLARE_READ_WRITE_FILE_BASIC(type) \
  template <> bool ReadFile(File& f, type& t); \
  template <> bool WriteFile(File& f, const type& t); \
  template <> bool ReadArrayFile(File& f, type* t,int n); \
  template <> bool WriteArrayFile(File& f, const type* t,int n); \

#define _DEFINE_READ_WRITE_FILE_BASIC(type) \
  template <> bool ReadFile(File& f, type& t) \
  { return f.ReadData(&t, sizeof(t)); } \
  template <> bool WriteFile(File& f, const type& t) \
  { return f.WriteData(&t, sizeof(t)); } \
  template <> bool ReadArrayFile(File& f, type* t,int n) \
  { return f.ReadData(t, sizeof(type)*n); } \
  template <> bool WriteArrayFile(File& f, const type* t,int n) \
  { return f.WriteData(t, sizeof(type)*n); }


_DECLARE_READ_WRITE_FILE_BASIC(bool);
_DECLARE_READ_WRITE_FILE_BASIC(char);
_DECLARE_READ_WRITE_FILE_BASIC(signed char);
_DECLARE_READ_WRITE_FILE_BASIC(unsigned char);
_DECLARE_READ_WRITE_FILE_BASIC(short);
_DECLARE_READ_WRITE_FILE_BASIC(unsigned short);
_DECLARE_READ_WRITE_FILE_BASIC(int);
_DECLARE_READ_WRITE_FILE_BASIC(unsigned int);
_DECLARE_READ_WRITE_FILE_BASIC(long);
_DECLARE_READ_WRITE_FILE_BASIC(unsigned long);
_DECLARE_READ_WRITE_FILE_BASIC(long long);
_DECLARE_READ_WRITE_FILE_BASIC(unsigned long long);
_DECLARE_READ_WRITE_FILE_BASIC(float);
_DECLARE_READ_WRITE_FILE_BASIC(double);
_DECLARE_READ_WRITE_FILE_BASIC(long double);

#endif
