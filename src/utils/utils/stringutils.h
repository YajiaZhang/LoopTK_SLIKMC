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

#ifndef UTILS_STRING_H
#define UTILS_STRING_H

#include <string>

/** @file utils/stringutils.h
 * @ingroup Utils
 * @brief Utilities for string manipulation.
 */

/** @addtogroup Utils */
/*@{*/

///Returns a "close bracket" character opposite c
char CloseBracket(char c);

///Turns the string into lower/uppercase
void Lowercase(char* str);
void Uppercase(char* str);
void Lowercase(std::string& str);
void Uppercase(std::string& str);

///Replace all instances of strfind with strreplace in str
int ReplaceAll(std::string& str,const char* strfind,const char* strreplace);

bool IsValidCToken(const char* str);
bool IsValidInteger(const char* str);
bool IsValidFloat(const char* str);

///Detects a pattern in str = [prefix][digits][suffix].
///Returns the number specified by [digits], or -1 if no such pattern is found.
int DetectNumericalPattern(const char* str,char prefix[],char suffix[],int& numDigits);


///Dos-unix endline conversion
int LengthWithDOSEndlines(const char* str);
bool EndlinesToDOS(const char* str,char* out,int max);
bool EndlinesFromDOS(const char* str,char* out,int max);
void EndlinesToDOS(std::string& str);
void EndlinesFromDOS(std::string& str);

///Returns pointer to "ext" for str="filename.ext"
const char* FileExtension (const char* str);
///Replaces the file extension of str with ext, or concatenates .ext onto str
void ChangeFileExtension (char* str, const char* ext);
///Returns "file.ext" for the str="dir1/dir2/.../file.ext"
char* GetFileName(char* str);
///Extracts the path from str (formatted as above) into buf
void GetFilePath(const char* str, char* buf);
///Removes the file extension of str
void StripExtension(char* str);

#ifdef WIN32
#define WIN32_LEAN_AND_MEAN
#include "windows.h"
void ToWideChar(const char* str, WCHAR* buf, int maxBuf);
#endif //WIN32

/*@}*/

#endif
