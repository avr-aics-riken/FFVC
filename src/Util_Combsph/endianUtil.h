/*
 * SolvSkl - Solver Skeleton - a Framework for Physical Simulation
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2004-2011
 *
 */

/**
 * @file include/endianUtil.h
 * エンディアンユーティリティマクロ・関数ファイル
 * @author     riken
 * @date       2010/10/01
 */

#ifndef _ENDIAN_UTIL_H_
#define _ENDIAN_UTIL_H_

#include <stdio.h>
#include <sys/types.h>
// modify(Vsphere 1.6.0) start for vsphere_win at 2008/06/20 by @hira
#ifndef _WIN32
#include <unistd.h>           // for linux
#else
#include "sph_win32_util.h"   // for win32
#endif
// modify(Vsphere 1.6.0) end for vsphere_win at 2008/06/20 by @hira
#include <iostream>
#include <string>

#ifdef SUPER_UX
  inline bool BigEndianSys() {return true;}
  inline bool LittleEndianSys() {return false;}
#else

// modify(Vsphere 1.6.0) start for vsphere_win at 2008/06/20 by @hira
//  inline bool BigEndianSys() {return (BYTE_ORDER != LITTLE_ENDIAN);}
//  inline bool LittleEndianSys() {return (BYTE_ORDER == LITTLE_ENDIAN);}
#ifdef _WIN32
  inline bool BigEndianSys() {return false;}
  inline bool LittleEndianSys() {return true;}
#else
// modify SPARC at 20110224 by @hira 
#ifdef FUJITSU_FX1            // fx1
  inline bool BigEndianSys() {return true;}
  inline bool LittleEndianSys() {return false;}
#else                         // ricc
  inline bool BigEndianSys() {return (BYTE_ORDER != LITTLE_ENDIAN);}
  inline bool LittleEndianSys() {return (BYTE_ORDER == LITTLE_ENDIAN);}
#endif

#endif
// modify(Vsphere 1.6.0) end for vsphere_win at 2008/06/20 by @hira

#endif // SUPER_UX

#ifdef SUPER_UX
  template<class X> inline void BSWAP16(X& x) {
    register unsigned char* _x_v = (unsigned char*)(&(x));
    unsigned char tmp;
    tmp = _x_v[0];
    _x_v[0] = _x_v[1];
    _x_v[1] = tmp;
  }
#else // SUPER_UX
  #ifndef BSWAP16
  #define BSWAP_X_16(x) \
    ( (((x) & 0xff00) >>  8) \
    | (((x) & 0x00ff) <<  8) )
  #define BSWAP16(x) { \
    register unsigned short& _x_v = (unsigned short&)(x); \
    _x_v = BSWAP_X_16(_x_v);}
  #endif // BSWAP16
#endif // SUPER_UX

#ifdef SUPER_UX
  template<class X> inline void BSWAP32(X& x) {
    register unsigned char* _x_v = (unsigned char*)(&(x));
    unsigned char tmp[2];
    tmp[0] = _x_v[0];
    tmp[1] = _x_v[1];
    _x_v[0] = _x_v[3];
    _x_v[1] = _x_v[2];
    _x_v[2] = tmp[1];
    _x_v[3] = tmp[0];
  }
#else // SUPER_UX
  #ifndef BSWAP32
  #define BSWAP_X_32(x) \
    ( (((x) & 0xff000000) >> 24) \
    | (((x) & 0x00ff0000) >>  8) \
    | (((x) & 0x0000ff00) <<  8) \
    | (((x) & 0x000000ff) << 24) )
  #define BSWAP32(x) \
    {register unsigned int& _x_v = (unsigned int&)(x); \
     _x_v = BSWAP_X_32(_x_v);}
  #endif // BSWAP32
#endif // SUPER_UX

#ifdef SUPER_UX
  template<class X> inline void BSWAP64(X& x) {
    register unsigned char* _x_v = (unsigned char*)(&(x));
    unsigned char tmp[4];
    tmp[0] = _x_v[0];
    tmp[1] = _x_v[1];
    tmp[2] = _x_v[2];
    tmp[3] = _x_v[3];
    _x_v[0] = _x_v[7];
    _x_v[1] = _x_v[6];
    _x_v[2] = _x_v[5];
    _x_v[3] = _x_v[4];
    _x_v[4] = tmp[3];
    _x_v[5] = tmp[2];
    _x_v[6] = tmp[1];
    _x_v[7] = tmp[0];
  }
#else // SUPER_UX
  #ifndef BSWAP64
  #define BSWAP_X_64(x) \
    ( (((x) & 0xff00000000000000ull) >> 56) \
    | (((x) & 0x00ff000000000000ull) >> 40) \
    | (((x) & 0x0000ff0000000000ull) >> 24) \
    | (((x) & 0x000000ff00000000ull) >>  8) \
    | (((x) & 0x00000000ff000000ull) <<  8) \
    | (((x) & 0x0000000000ff0000ull) << 24) \
    | (((x) & 0x000000000000ff00ull) << 40) \
    | (((x) & 0x00000000000000ffull) << 56) )
  #define BSWAP64(x) \
    {register unsigned long long& _x_v = (unsigned long long&)(x); \
     _x_v = BSWAP_X_64(_x_v);}
  #endif // BSWAP64
#endif // SUPER_UX

#ifdef SUPER_UX
  template<class X, class Y> inline void SBSWAPVEC(X* a, Y n) {
    register unsigned int nn = (unsigned int)n;
    for(register unsigned int _i=0;_i<nn;_i++){
      register unsigned short _x_v = (unsigned short)a[_i];
      BSWAP16(_x_v);
      a[_i] = _x_v;
    }
  }
#else // SUPER_UX
  #ifndef SBSWAPVEC
  #define SBSWAPVEC(a,n) do{\
    for(register unsigned int _i=0;_i<(n);_i++){BSWAP16(a[_i]);}\
  }while(0)
  #endif // SBSWAPVEC
#endif // SUPER_UX

#ifdef SUPER_UX
  template<class X, class Y> inline void BSWAPVEC(X* a, Y n) {
    register unsigned int nn = (unsigned int)n;
    for(register unsigned int _i=0;_i<nn;_i++){
      register unsigned int _x_v = (unsigned int)a[_i];
      BSWAP32(_x_v);
      a[_i] = _x_v;
    }
  }
#else // SUPER_UX
  #ifndef BSWAPVEC
   
  #define BSWAPVEC(a,n) do{\
    for(register unsigned int _i=0;_i<(n);_i++){BSWAP32(a[_i]);}\
  }while(0)
  #endif // BSWAPVEC
#endif // SUPER_UX

#ifdef SUPER_UX
  template<class X, class Y> inline void DBSWAPVEC(X* a, Y n) {
    register unsigned int nn = (unsigned int)n;
    for(register unsigned int _i=0;_i<nn;_i++){
      register unsigned long long _x_v = (unsigned long long)a[_i];
      BSWAP64(_x_v);
      a[_i] = _x_v;
    }
  }
#else // SUPER_UX
  #ifndef DBSWAPVEC
  #define DBSWAPVEC(a,n) do{\
    for(register unsigned int _i=0;_i<(n);_i++){BSWAP64(a[_i]);}\
  }while(0)
  #endif // DBSWAPVEC
#endif // SUPER_UX

typedef enum {UnKnown =0, Match, UnMatch} EMatchType;

inline EMatchType isMatchEndian(const char* fname, int val) {
  EMatchType ret = UnKnown;
  if( !fname ) return ret;

  FILE* fp = NULL;
  if( !(fp = fopen(fname, "rb")) ) {
    fprintf(stderr, "Can't open file.(%s)\n", fname);
    return ret;
  }

  union {char cBuff[4]; int iBuff;} readBuf;
  if ( fread(readBuf.cBuff, 4, 1, fp) < 1 ) {
    fclose(fp);
    return ret;
  }

  if ( readBuf.iBuff == (int)val ) ret = Match;

  BSWAP32(readBuf.iBuff);
  if ( readBuf.iBuff == (int)val ) ret = UnMatch;

  fclose(fp);
  return ret;
}

inline EMatchType isMatchEndian(const std::string& fname, int val) {
  EMatchType ret = UnKnown;
  if( fname.size() == 0 ) return ret;

  FILE* fp = NULL;
  if( !(fp = fopen(fname.c_str(), "rb")) ) {
    fprintf(stderr, "Can't open file.(%s)\n", fname.c_str());
    return ret;
  }

  union {char cBuff[4]; int iBuff;} readBuf;
  if ( fread(readBuf.cBuff, 4, 1, fp) < 1 ) {
    fclose(fp);
    return ret;
  }

  if ( readBuf.iBuff == (int)val ) ret = Match;

  BSWAP32(readBuf.iBuff);
  if ( readBuf.iBuff == (int)val ) ret = UnMatch;

  fclose(fp);
  return ret;
}

typedef enum { SKL_SEEK_SET=1, SKL_SEEK_CUR=2, SKL_SEEK_END=3} SklSeekType;
#define SKL_FSEEK_BUF_SIZE 1024

/*
 *  "_USE_SKL_FSEEK" is made active.
 *  when it uses it in the environment for which fseek cannot be used,
 */
//#define _USE_SKL_FSEEK

inline bool FSeek(FILE* fp, size_t offset, SklSeekType seek_type) {
#ifndef _USE_SKL_FSEEK
  if( !fp ) return false;
  if( seek_type == SKL_SEEK_SET ) {
    if( fseek(fp, (long)offset, SEEK_SET) != 0 ) return false;
  }
  else if( seek_type == SKL_SEEK_CUR ) {
    if( fseek(fp, (long)offset, SEEK_CUR) != 0 ) return false;
  }
  else if( seek_type == SKL_SEEK_END ) {
    if( fseek(fp, (long)offset, SEEK_END) != 0 ) return false;
  }
  else return false;
  return true;
#else
  if( !fp ) {
    fprintf(stderr, "FSeek : File pointer is NULL.\n"); fflush(stderr);
    return false;
  }
  fpos_t mypos;
  int ret;
  if( (ret = fgetpos(fp, &mypos)) != 0 ) {
    perror("ERROR : fgetpos ");
    fprintf(stderr, "%s (%d): fgetpos\n", __FILE__, __LINE__); fflush(stderr);
    return false;
  }

  register size_t n;
  size_t bufsz = SKL_FSEEK_BUF_SIZE;
  char buf[SKL_FSEEK_BUF_SIZE];
  size_t loopMax = offset / bufsz;
  size_t surplus = offset - (loopMax * bufsz);

  if( seek_type == SKL_SEEK_SET ) {
    rewind(fp);
    for(n=0; n<loopMax; n++) {
      if( fread(buf, sizeof(char), bufsz, fp) != bufsz ) {
        if( feof(fp) != 0 ) {
          fprintf(stderr, "FSeek : EOF was detected.\n"); fflush(stderr);
          return true;
        }
        if( (ret = fsetpos(fp, &mypos)) != 0 ) {
          perror("ERROR : fsetpos ");
        }
        perror("ERROR : fread ");
        return false;
      }
    }
    if( surplus > 0 ) {
      if( fread(buf, sizeof(char), surplus, fp) != surplus ) {
        if( feof(fp) != 0 ) {
          fprintf(stderr, "FSeek : EOF was detected.\n"); fflush(stderr);
          return true;
        }
        if( (ret = fsetpos(fp, &mypos)) != 0 ) {
          perror("ERROR : fsetpos ");
        }
        perror("ERROR : fread ");
        return false;
      }
    }
  }
  else if( seek_type == SKL_SEEK_CUR ) {
    for(n=0; n<loopMax; n++) {
      if( fread(buf, sizeof(char), bufsz, fp) != bufsz ) {
        if( feof(fp) != 0 ) {
          fprintf(stderr, "FSeek : EOF was detected.\n"); fflush(stderr);
          return true;
        }
        if( (ret = fsetpos(fp, &mypos)) != 0 ) {
          perror("ERROR : fsetpos ");
        }
        perror("ERROR : fread ");
        return false;
      }
    }
    if( surplus > 0 ) {
      if( fread(buf, sizeof(char), surplus, fp) != surplus ) {
        if( feof(fp) != 0 ) {
          fprintf(stderr, "FSeek : EOF was detected.\n"); fflush(stderr);
          return true;
        }
        if( (ret = fsetpos(fp, &mypos)) != 0 ) {
          perror("ERROR : fsetpos ");
        }
        perror("ERROR : fread ");
        return false;
      }
    }
  }
  else if( seek_type == SKL_SEEK_END ) {
    while( fread(buf, sizeof(char), bufsz, fp) == bufsz ) { ; }
    if( feof(fp) != 0 ) return true;
    fprintf(stderr, "FSeek : SEEK_END failed.\n"); fflush(stderr);
    return false;
  }
  else {
    fprintf(stderr, "Unknown seek type.\n"); fflush(stderr);
    return false;
  }
  return true;
#endif // _USE_SKL_FSEEK
}

/* Do not use it, caused by the specification of MPI2.
inline EMatchType isMatchEndian(FILE* fp, int val) {
  EMatchType ret = UnKnown;
  if ( ! fp ) return ret;
  long here = ftell(fp);
  if ( here < 0 ) return ret;


  union {char cBuff[4]; int iBuff;} readBuf;
  if ( fread(readBuf.cBuff, 4, 1, fp) < 1 ) return ret;

  if ( readBuf.iBuff == (int)val )
    ret = Match;

  BSWAP32(readBuf.iBuff);
  if ( readBuf.iBuff == (int)val )
    ret = UnMatch;

  (void)fseek(fp, here, SEEK_SET);
  return ret;
}

inline EMatchType isMatchEndian(int fd, int val) {
  EMatchType ret = UnKnown;
  if ( fd < 0 ) return ret;

  union {char cBuff[4]; int iBuff;} readBuf;
  if ( read(fd, readBuf.cBuff, 4) < 4 ) return ret;

  if ( readBuf.iBuff == (int)val )
    ret = Match;

  BSWAP32(readBuf.iBuff);
  if ( readBuf.iBuff == (int)val )
    ret = UnMatch;

  (void)lseek(fd, 0, SEEK_SET);
  return ret;
}

inline EMatchType isMatchEndian(std::istream& is, const int magick) {
  EMatchType ret = UnKnown;
  union {char cBuff[4]; int iBuff;} readBuf;
  std::ios::pos_type here = is.tellg();

  is.read(readBuf.cBuff, 4);
  if ( ! is.good() ) return ret;

  if ( readBuf.iBuff == magick )
    ret = Match;
  else {
    BSWAP32(readBuf.iBuff);
    if ( readBuf.iBuff == magick )
      ret = UnMatch;
  }
  is.seekg(here);

  return ret;
}
*/

#endif // _ENDIAN_UTIL_H_

