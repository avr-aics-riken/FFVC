//##################################################################################
//
// FFV-C ASD module : Frontflow / violet Cartesian Active SubDomain
//
// Copyright (c) 2014 Advanced Institute for Computational Science, RIKEN.
// All rights reserved.
//
//##################################################################################

/**
 * @file   SubDomain.C
 * @brief  Active subdomain
 * @author aics
 */

#include "SubDomain.h"
#include <fstream>
#include <iostream>

namespace ASubDomain {
  
// #################################################################
// Subdomain.datの入力
ActiveSubdomain* ActiveSubdomain::LoadActiveSubdomain(const string& filepath)
{
  FILE *fp = fopen(filepath.c_str(), "rb");
  if (!fp) {
    printf("Missing Subdomain File. (%s)", filepath.c_str());
    return NULL;
  }
  
  unsigned FMT = 0;
  if ( fread(&FMT, sizeof(unsigned), 1, fp) != 1) {
    printf("File Read Error. (%s)", filepath.c_str());
    return NULL;
  }
  
  bool isNeedSwap = false;
  // Format Check
  if ( FMT != SUBDOMAIN_IDENTIFIER )
  {
    BSwap32(&FMT);
    if ( FMT != SUBDOMAIN_IDENTIFIER ) {
      printf("File Format Error. (%s)", filepath.c_str());
      return NULL;
    }
    isNeedSwap = true;
  }
  
  unsigned size[3] = {0, 0, 0};
  if ( fread(size, sizeof(unsigned), 3, fp) != 3 ) {
    printf("File Read Error. (%s)", filepath.c_str());
    return NULL;
  }
  if ( isNeedSwap ){
    BSwap32(&size[0]);
    BSwap32(&size[1]);
    BSwap32(&size[2]);
  }
  
  if ( size[0] == 0 || size[1] == 0 || size[2] == 0 ) {
    printf("Data is invalid. (%s)", filepath.c_str());
    return NULL;
  }
  
  if ( size[0] != x || size[1] != y || size[2] != z ) {
    printf("Data size is different from (%d, %d, %d)", x, y, z);
    return NULL;
  }
  
  ActiveSubdomain *sb = new ActiveSubdomain(x, y, z);
  
  if ( fread(sb->contents, sizeof(unsigned char), x*y*z, fp) != x*y*z )
  {
    printf("Data is invalid. (%s)", filepath.c_str());
    delete [] sb;
    return NULL;
  }
  
  fclose(fp);
  
  return sb;
}


// #################################################################
// Subdomain.datの出力
bool ActiveSubdomain::SaveActiveSubdomain(const string& path)
{
  FILE *fp = fopen(path.c_str(), "wb");
  if (!fp) {
    printf("File Open Error. (%s)\n", path.c_str());
    return false;
  }
  
  unsigned FMT = ('S' | ('B' << 8) | ('D' << 16) | ('M' << 24));
  if ( fwrite(&FMT, sizeof(unsigned), 1, fp) != 1)
  {
    printf("File Write Error. (%s)\n", path.c_str());
    fclose(fp);
    return false;
  }
  
  const unsigned size[3] = { x, y, z };
  if ( fwrite(size, sizeof(unsigned), 3, fp) != 3 )
  {
    printf("File Write Error. (%s)\n", path.c_str());
    fclose(fp);
    return false;
  }
  
  
  if ( fwrite(contents, sizeof(unsigned char), x*y*z, fp) !=  x*y*z )
  {
    printf("File Write Error. (%s)\n", path.c_str());
    fclose(fp);
    return false;
  }
  
  fclose(fp);
  
  return true;
}


// #################################################################
// svxフォーマットで出力する(IDのみ)
bool ActiveSubdomain::writeSVX(const string& path, const float* pch, const REAL_TYPE* org)
{
  ofstream ofs(path.c_str(), ios::out | ios::binary);
  if (!ofs) {
    cout << "\tCan't open " << path.c_str() << " file" << endl;
    return false;
  }
  
  float dx = pch[0];
  float dy = pch[1];
  float dz = pch[2];
  float ox = org[0];
  float oy = org[1];
  float oz = org[2];
  
  int *q = new int[x*y*z];
  
  for (int k=1; k<=z; k++) {
    for (int j=1; j<=y; j++) {
      for (int i=1; i<=x; i++) {
        size_t m = _F_IDX_S3D(i, j, k, x, y, z, 0);
        q[m] = (int)contents[m];
      }
    }
  }
  
  int sz;
  
  // voxel size
  sz = sizeof(int)*3;
  ofs.write( (char*)&sz, sizeof(int) );
  ofs.write( (char*)&x,  sizeof(int) );
  ofs.write( (char*)&y,  sizeof(int) );
  ofs.write( (char*)&z,  sizeof(int) );
  ofs.write( (char*)&sz, sizeof(int) );
  
  // original point of domain
  sz = sizeof(float)*3;
  ofs.write( (char*)&sz, sizeof(int) );
  ofs.write( (char*)&ox, sizeof(float) );
  ofs.write( (char*)&oy, sizeof(float) );
  ofs.write( (char*)&oz, sizeof(float) );
  ofs.write( (char*)&sz, sizeof(int) );
  
  // pitch of voxel
  ofs.write( (char*)&sz, sizeof(int) );
  ofs.write( (char*)&dx, sizeof(float) );
  ofs.write( (char*)&dy, sizeof(float) );
  ofs.write( (char*)&dz, sizeof(float) );
  ofs.write( (char*)&sz, sizeof(int) );
  
  // type of stored data
  sz = sizeof(int)*1;
  int dtype = 0;
  dtype |= ( 0x1<<2 );  // medium ID
  ofs.write( (char*)&sz,  sizeof(int) );
  ofs.write( (char*)&dtype, sizeof(int) );
  ofs.write( (char*)&sz,  sizeof(int) );
  
  // medium ID
  sz = x*y*z * sizeof(int);
  ofs.write( (char*)&sz, sizeof(int) );
  ofs.write( (char*)q,   sz );
  ofs.write( (char*)&sz, sizeof(int) );
  
  ofs.close();
  
  if (q) { delete [] q; q=NULL; }
  
  return true;
}
  
  } // namespace ASubDomain
