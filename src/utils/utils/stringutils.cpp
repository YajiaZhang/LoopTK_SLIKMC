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

#include "stringutils.h"
#include <string.h>
#include <ctype.h>
#include <stdlib.h>

char CloseBracket(char c)
{
  switch(c) {
  case '[': return ']';
  case ']': return '[';
  case '(': return ')';
  case ')': return '(';
  case '<': return '>';
  case '>': return '<';
  case '`': return '\'';
  case '\'': return '`';
  case '\\': return '/';
  case '/': return '\\';
  default: return c;
  }
}

void Lowercase(char* str)
{
	for(;*str;str++)
		(*str) = tolower(*str);
}

void Uppercase(char* str)
{
	for(;*str;str++)
		(*str) = toupper(*str);
}

void Lowercase(std::string& str)
{
	for(unsigned int i=0;i<str.length();i++)
		str[i] = tolower(str[i]);
}

void Uppercase(std::string& str)
{
	for(unsigned int i=0;i<str.length();i++)
		str[i] = toupper(str[i]);
}

int ReplaceAll(std::string& str,const char* strfind,const char* strreplace)
{
	size_t lenfind=strlen(strfind);
	size_t lenrep=strlen(strreplace);
	int n=0;
	size_t j=0;
	for(;(j=str.find(strfind,j))!=str.npos;n++) {
		str.replace(j,lenfind,strreplace);
		j+=lenrep;
	}
	return n;
}

bool IsValidCToken(const char* str)
{
	//RE is (_|[alpha])(_|[alnum])*
	if(!str) return false;
	if(!*str) return false;
	if(isdigit(*str))
		return false;
	while(*str) {
		if(!(isalnum(*str) || *str=='_'))
			return false;
		str++;
	}
	return true;
}

bool IsValidInteger(const char* str)
{
	//RE is [+|-]?[digit]+
	if(!str) return false;
	if(!*str) return false;
	if(*str=='-' || *str=='+') str++;
	if(!isdigit(*str)) return false;
	while(*str) {
		if(!isdigit(*str))
			return false;
		str++;
	}
	return true;
}

bool IsValidFloat(const char* str)
{
	//RE is [+|-]?([digit]+|[digit]+.[digit]*|.[digit]+)((e|E)[+|-]?[digit]+)?
	if(!str) return false;
	if(!*str) return false;
	if(*str=='-' || *str=='+') str++;
	bool readExponent = false;
	if(isdigit(*str)) {  //branch 1 [digit]+|[digit]+.[digit]*
		str++;
		bool readDot=false;
		while(*str && !readExponent) {
			if(*str == '.') {
				if(readDot) return false;
				readDot=true;
			}
			else if(*str == 'e' || *str == 'E')
				readExponent = true;
			else if(!isdigit(*str)) return false;
			str++;
		}
	}
	else if(*str == '.') { //branch 2 .[digit]+
		str++;
		if(!isdigit(*str)) return false;
		str++;
		while(*str && !readExponent) {
			if(*str == 'e' || *str == 'E')
				readExponent = true;
			else if(!isdigit(*str)) return false;
			str++;
		}
	}
	else return false;

	if(readExponent) {
		if(!IsValidInteger(str)) return false;
	}
	return true;
}

/// Detects a pattern in str = [prefix][digits][suffix]
/// Returns the integer in [digits], or -1 if no such pattern is found.
int DetectNumericalPattern(const char* str,char prefix[],char suffix[],int& numDigits)
{
  int n=strlen(str);
  int beginDigits=n,endDigits=n;
  for(int i=0;i<n;i++) {
    if(isdigit(str[i])) {
      beginDigits=i;
      break;
    }
  }
  if(beginDigits == n) return -1;
  for(int i=beginDigits;i<n;i++) {
    if(!isdigit(str[i])) {
      endDigits=i;
      break;
    }
  }
  numDigits = endDigits-beginDigits;
  strncpy(prefix,str,beginDigits);
  prefix[beginDigits] = 0;
  strncpy(suffix,str+endDigits,n-endDigits);
  suffix[n-endDigits] = 0;
  char* buf = new char[n];
  strncpy(buf,str+beginDigits,endDigits-beginDigits);
  buf[endDigits-beginDigits] = 0;
  int val = atoi(buf);
  delete buf;
  return val;
}






int LengthWithDOSEndlines(const char* str)
{
	int i=0;
	bool return_read=false;
	while(*str) {
		switch(*str) {
		case '\r':
			return_read=true;
			break;
		case '\n':
			i+=2;
			break;
		default:
			if(return_read) {
				i+=2;
				return_read=false;
			}
			i++;
			break;
		}
		str++;
	}
	if(return_read)
		i+=2;
	return i;
}

bool EndlinesToDOS(const char* str,char* out,int max)
{
	int i=0;
	bool return_read=false;
	while(*str) {
		if(i>=max) return false;
		switch(*str) {
		case '\r':
			return_read=true;
			break;
		case '\n':
			if(i+1>=max) return false;
			out[i]='\r';
			out[i+1]='\n';
			return_read=false;
			i+=2;
			break;
		default:
			if(return_read) {
				if(i+2>=max) return false;
				out[i]='\r';
				out[i+1]='\n';
				i+=2;
				return_read=false;
			}
			out[i]=*str;
			i++;
			break;
		}
		str++;
	}
	if(return_read) {
		if(i+2>=max) return false;
		out[i]='\r';
		out[i+1]='\n';
		i+=2;
		return_read=false;
	}
	if(i>=max) return false;
	out[i]=0;
	return true;
}

bool EndlinesFromDOS(const char* str,char* out,int max)
{
	int i=0;
	while(*str) {
		if(i>=max) return false;
		if(*str != '\r') out[i]=*str;
	}
	if(i>=max) return false;
	out[i]=0;
	return true;
}

void EndlinesToDOS(std::string& str)
{
	ReplaceAll(str,"\r\n","\n");  //make sure we dont duplicate \n's
	ReplaceAll(str,"\n","\r\n");
}

void EndlinesFromDOS(std::string& str)
{
	ReplaceAll(str,"\r\n","\n");
}

const char* FileExtension (const char* str)
{
	const char* dp = strrchr(str,'.');
	if(dp == NULL) return NULL;
	dp++;

	return dp;
}

void ChangeFileExtension (char* str, const char* ext)
{
	char* dp = strrchr(str,'.');
	if(dp == NULL)
	{
		strcat(str,".");
		strcat(str,ext);
		return;
	}

	dp++;

	strcpy(dp, ext);
}

char* GetFileName(char* str)
{
	char* fp = strrchr(str,'\\');
	if(fp == NULL)
	{
		fp = strrchr(str,'/');
		if(fp == NULL)
			return str;
	}
	fp++;

	return fp;
}


void GetFilePath(const char* str, char* buf)
{
	strcpy(buf, str);
	char* fp = strrchr(buf,'\\');
	if(fp == NULL)
	{
		fp = strrchr(str,'/');
		if(fp == NULL)
			return;
	}
	fp++;
	*fp = '\0';
}

void StripExtension(char* str)
{
	char* dp = strrchr(str,'.');
	if(!dp)
		return;
	*dp = '\0';
}

#ifdef WIN32

void ToWideChar(const char* str, WCHAR* buf, int maxBuf)
{
	buf[0] = 0;
	int res = MultiByteToWideChar(CP_ACP, 0, str, -1, buf,maxBuf);
	if(res == 0)
	{
	  cerr<<"Couldnt' convert the string to wide characters"<<endl;
	  abort();
	}
}

#endif //WIN32
