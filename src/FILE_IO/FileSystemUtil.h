//##################################################################################
//
// FFV-C : Frontflow / violet Cartesian
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
//
///
/// @file FileSystemUtil.h
/// @brief ファイル操作関連ユーティリティ
///

#ifndef __FFV_FILESYSTEM_UTIL_H__
#define __FFV_FILESYSTEM_UTIL_H__

#include <string>
#include <stdlib.h>

namespace BVX_IO {

	/// ファイルパスの "\\" を "/" に変換
	/// return : "/"に変換されたファイルパス
	inline std::string ConvertPath(const std::string& path)
	{
		std::string r = path;
		size_t p = r.find("\\");
		while (p != std::string::npos)
		{
			*(r.begin() + p) = '/';
			p = r.find("\\");
		}
		return r;
	}

	/// pathからディレクトリ名を抜き出す
	/// return : ディレクトリ名 (最後の"/"は残す)
	inline std::string GetDirectory(const std::string& path)
	{
		std::string cpath = ConvertPath(path);

		std::string dir;
		size_t p = cpath.rfind("/");
		if(p != std::string::npos)
		{
			dir = cpath.substr(0, p+1);
		}else{
			dir = std::string("./");
		}
		return dir;
	}
	
	/// pathから拡張子以前のファイル名を抜き出す
	inline std::string GetFilePrefix(const std::string& path)
	{
		std::string cpath = ConvertPath(path);
		std::string filename = cpath.substr(cpath.rfind("/")+1);
		
		return filename.substr(0, filename.rfind("."));
	}

	/// dirをディレクトリ名として整形
	/// 空の場合 : "./"
	/// 文字列が入っている場合最後に "/"を追加
	inline std::string FixDirectoryPath(const std::string& dir)
	{
		if( dir == std::string("") ){
			return std::string("./");
		}else if( dir.rfind("/") != dir.length()-1 ){
			return dir + std::string("/");
		}else{
			return dir;
		}
	}

	/// pathで指定したディレクトリを作成 (mkdir -p 相当)
	/// ディレクトリ作成に失敗した場合、falseを返す
	bool CreateDirectory(const std::string& path);

  
  void CheckDir(std::string dirstr);
  
  
  //////////////////////////
  // from CDMlib
  //////////////////////////
  

  
  // A:などのドライブパスがあればtrue
  inline bool hasDrivePath(const std::string& path)
  {
    if ( path.size() < 2 ) return false;
    
    char x = path[0];
    if ( ((x >= 'A' && x <= 'Z' ) || (x >= 'a' && x <= 'z')) && path[1] == ':' )
      return true;
    return false;
  }
  
  
  //
  inline std::string emitDrivePath(std::string& path)
  {
    // returns drive (ex. 'C:')
    if ( !hasDrivePath(path) ) return std::string();
    std::string driveStr = path.substr(0, 2);
    path = path.substr(2);
    return driveStr;
  }
  
  
  // デリミタを返す
  inline char getDelimCharPath()
  {
#ifdef WIN32
    return '\\';
#else
    return '/';
#endif
  }
  
  
  // true  : Absolute Path(絶対パス)
  // false : Relative Path(相対パス)
  inline bool isAbsolutePath(const std::string& path)
  {
    std::string xpath(path);
    emitDrivePath(xpath);
    char c1, c2;
    c1 = xpath[0];
    c2 = getDelimCharPath();
    return (c1 == c2);
  }

  
  inline std::string DirNamePath(const std::string& path, const char dc = getDelimCharPath())
  {
    char* name = strdup( path.c_str() );
    char* p = name;
    
    for ( ; ; ++p ) {
      if ( ! *p ) {
        if ( p > name ) {
          char rs[2] = {dc, '\0'};
          return rs;
        } else {
          char rs[3] = {'.', dc, '\0'};
          return rs;
        }
      }
      if ( *p != dc ) break;
    }
    
    for ( ; *p; ++p );
    while ( *--p == dc ) continue;
    *++p = '\0';
    
    while ( --p >= name )
      if ( *p == dc ) break;
    ++p;
    if ( p == name )
    {
      char rs[3] = {'.', dc, '\0'};
      return rs;
    }
    
    while ( --p >= name )
      if ( *p != dc ) break;
    ++p;
    
    *p = '\0';
    if( p == name ) {
      char rs[2] = {dc, '\0'};
      return rs;
    }
    else
    {
      std::string s( name );
      free( name );
      
      if( !isAbsolutePath(s) )
      {
        const char *q = s.c_str();
        if( q[0] != '.' && q[1] != '/' )
        {
          char rs[3] = {'.', dc, '\0'};
          s = std::string(rs) + s;
        }
      }
      return s;
    }
  }


} // namespace BVX_IO


#endif // __FFV_FILESYSTEM_UTIL_H__

