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

#include "ppm.h"
#include <stdio.h>

bool WritePPM_RGB_ASCII(unsigned char image[],int m,int n,const char* file)
{
  FILE* f = fopen(file,"w");
  if(!f) {
    return false;
  }
  fprintf(f,"P3\n#%s\n",file);
  fprintf(f,"%d %d\n",m,n);
  fprintf(f,"255\n");
  int k=0;
  for(int i=0;i<m;i++) {
    for(int j=0;j<n;j++,k+=3) {
      fprintf(f,"%d %d %d  ",image[k],image[k+1],image[k+2]);
    }
    fprintf(f,"\n");
  }
  fclose(f);
  return true;
}

bool WritePPM_Grayscale_ASCII(unsigned char image[],int m,int n,const char* file)
{
  FILE* f = fopen(file,"w");
  if(!f) {
    return false;
  }
  fprintf(f,"P2\n#%s\n",file);
  fprintf(f,"%d %d\n",m,n);
  fprintf(f,"255\n");
  int k=0;
  for(int i=0;i<m;i++) {
    for(int j=0;j<n;j++,k++) {
      fprintf(f,"%d ",image[k]);
    }
    fprintf(f,"\n");
  }
  fclose(f);
  return true;
}

bool WritePPM_RGB_Binary(unsigned char image[],int m,int n,const char* file)
{
  FILE* f = fopen(file,"wb");
  if(!f) {
    return false;
  }
  fprintf(f,"P6\n#%s\n",file);
  fprintf(f,"%d %d\n",m,n);
  fprintf(f,"255\n");
  fwrite(image,m*n*3,1,f);
  fprintf(f,"\n");
  fclose(f);
  return true;
}

bool WritePPM_Grayscale_Binary(unsigned char image[],int m,int n,const char* file)
{
  FILE* f = fopen(file,"wb");
  if(!f) {
    return false;
  }
  fprintf(f,"P5\n#%s\n",file);
  fprintf(f,"%d %d\n",m,n);
  fprintf(f,"255\n");
  fwrite(image,m*n,1,f);
  fprintf(f,"\n");
  fclose(f);
  return true;
}

