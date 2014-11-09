#ifndef _FB_UTIL_PATH_H_
#define _FB_UTIL_PATH_H_

//##################################################################################
//
// Flow Base class
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
 * @file utilPath.h
 */

#include <string>
#include <sstream>
#include <deque>
#include <string.h>

// for GetFullPathName
#if defined(_WIN32)
  #include <stdlib.h>
#else
  #include <sys/param.h>
  #include <stdlib.h>
#endif

#ifdef SUPER_UX
  extern "C" char* realpath(const char*,char*);
#endif // SUPER_UX

/**
 * path_util::ファイルパス関数ネームスペース
 */
namespace path_util {
  
  
  // #################################################################
  /**
   * パス名の絶対パスを取得する。
   * @param [in]  path               パス名
   * @param [out] resolved_path      絶対パス名
   * @param [in]  str_len            文字列長
   * @return  true=success, false=fault
   */
  inline bool GetFullPathName(const char* path,
                              char* resolved_path,
                              size_t str_len) {
    if( !path || !resolved_path || (str_len <= 0) ) return false;
    
#if defined(_WIN32)
    if( !_fullpath(resolved_path, path, str_len) ) return false;
#else
    char rpath[1024];
    if( !realpath(path, rpath) ) return false;
    if( strlen(rpath) > str_len ) return false;
    strncpy(resolved_path, rpath, str_len);
#endif
    
    return true;
  }
  
  
  // #################################################################
  /**
   * パス名からファイル名を切り出す
   * @param [in] path        パス名
   * @param [in] dc          パス区切り文字
   * @return ファイル名
   */
  inline std::string DirName(const std::string& path,
                             const char dc = '/') {
    char* name = strdup( path.c_str() );
    char* p = name;
    
    for ( ; ; ++p ) {
      if ( ! *p ) {
        if ( p > name ) {
          char rs[2] = {dc, '\0'};
          return rs;
        } else
          return(".");
      }
      if ( *p != dc ) break;
    }
    
    for ( ; *p; ++p );
    while ( *--p == dc ) continue;
    *++p = '\0';
    
    while ( --p >= name )
      if ( *p == dc ) break;
    ++p;
    if ( p == name ) return(".");
    
    while ( --p >= name )
      if ( *p != dc ) break;
    ++p;
    
    *p = '\0';
    if( p == name ) {
      char rs[2] = {dc, '\0'};
      return rs;
    } else {
      std::string s( name );
      free( name );
      return s;
    }
  }
  
  
  // #################################################################
  /**
   * パス名からディレクトリとサフィックスを取り除く。
   * @param [in] path                パス名
   * @param [in] suffix              削除末尾文字
   * @param [in] dc                  パス区切り文字
   * @return ファイル名（ディレクトリとサフィックスなし）
   */
  inline std::string BaseName(const std::string& path,
                              const std::string& suffix = std::string(""),
                              const char dc = '/') {
    char* name = strdup( path.c_str() );
    char* p = name;
    
    for ( ; ; ++p ) {
      if ( ! *p ) {
        if ( p > name ) {
          char rs[2] = {dc, '\0'};
          return rs;
        } else
          return "";
      }
      if ( *p != dc ) break;
    }
    
    for ( ; *p; ++p ) continue;
    while ( *--p == dc ) continue;
    *++p = '\0';
    
    while ( --p >= name )
      if ( *p == dc ) break;
    ++p;
    
    if ( suffix.length() > 0 ) {
      const int suffixlen = suffix.length();
      const int stringlen = strlen( p );
      if ( suffixlen < stringlen ) {
        const int off = stringlen - suffixlen;
        if ( !strcmp( p + off, suffix.c_str() ) )
          p[off] = '\0';
      }
    }
    
    std::string s( p );
    free( name );
    return s;
  }
  
  
  // #################################################################
  /**
   * パス名から相対パス文字(./, ../)を削除する。
   * @param [in] path        パス名
   * @param [in] dc          パス区切り文字
   * @return  相対パス文字(./, ../)削除パス名
   */
  inline std::string OmitDots(const std::string& path,
                              const char dc = '/') {
    using namespace std;
    if ( path.empty() ) return path;
    
    deque<string> elemLst;
    istringstream ss(path);
    std::string selStr; char c;
    
    // decomposition
    //    while ( 1 ) {
    unsigned int FLAG = 1;
    while ( FLAG ) {
      if ( ! ss.get(c) ) {
        elemLst.push_back(selStr);
        break;
      }
      if ( c == dc ) {
        elemLst.push_back(selStr);
        selStr = "";
      } else {
        selStr += c;
      }
    } // end of while(1)
    
    // remove null or '.' elem
    deque<string>::iterator it;
    for ( it = elemLst.begin(); it != elemLst.end(); ) {
      if ( it->empty() || (*it) == "." )
        it = elemLst.erase(it);
      else
        it ++;
    } // end of for(it)
    
    // remove '..' elem
    if ( ! elemLst.empty() )
      for ( it = elemLst.begin() + 1; it != elemLst.end(); ) {
        if ( (*it) == ".." ) {
          it = elemLst.erase(--it);
          it = elemLst.erase(it);
        } else
          it ++;
      } // end of for(it)
    
    // result
    string retPath;
    if ( path[0] == dc ) retPath = dc;
    for ( it = elemLst.begin(); it != elemLst.end(); it++ ) {
      retPath += (*it);
      if ( it != elemLst.end() -1 )
        retPath += dc;
    } // end of for(it)
    return retPath;
  }

} // end of namespace path_util

#endif // _FB_UTIL_PATH_H_
