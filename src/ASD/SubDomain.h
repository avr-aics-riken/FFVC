#ifndef _ASD_SUBDOMAIN_H_
#define _ASD_SUBDOMAIN_H_

//##################################################################################
//
// FFV-C ASD module : Frontflow / violet Cartesian Active SubDomain
//
// Copyright (c) 2007-2011 VCAD System Research Program, RIKEN.
// All rights reserved.
//
// Copyright (c) 2011-2014 Institute of Industrial Science, The University of Tokyo.
// All rights reserved.
//
// Copyright (c) 2012-2014 Advanced Institute for Computational Science, RIKEN.
// All rights reserved.
//
//##################################################################################

/**
 * @file   SubDomain.h
 * @brief  ASD Class Header
 * @author aics
 */

#include "../FB/FB_Define.h"
#include <stdio.h>
#include <string>
#include <string.h>

using namespace std;

namespace ASubDomain {
  
#define SUBDOMAIN_IDENTIFIER ('S' | ('B' << 8) | ('D' << 16) | ('M' << 24))

  class ActiveSubdomain {
    
  private:
    unsigned x, y, z;
    unsigned char *contents;
    
  public:
    
    // default constructor
    ActiveSubdomain() {};
    
    ActiveSubdomain(unsigned _x, unsigned _y, unsigned _z) {
      x = _x;
      y = _y;
      z = _z;
      contents = new unsigned char[x * y * z];
      memset(contents, 0, sizeof(unsigned char)*x*y*z);
    }
    
    ~ActiveSubdomain() {
      delete [] contents;
    };
    
    
  public:
    ActiveSubdomain* LoadActiveSubdomain(const string& filepath);
    
    bool SaveActiveSubdomain(const string& path);
    
    bool writeSVX(const string& path, const float* pch, const REAL_TYPE* org);
    
    unsigned char* get_ptr() {
      return contents;
    }
    
    /// for 4bytes
    static inline void BSwap32(void* a) {
      unsigned* x = (unsigned*)a;
      *x = ( (((*x) & 0xff000000) >> 24) | (((*x) & 0x00ff0000) >> 8 ) |
            (((*x) & 0x0000ff00) <<  8) | (((*x) & 0x000000ff) << 24) );
    }
  };

} // namespace ASubDomain

#endif // _ASD_SUBDOMAIN_H_
