/*
 * CPM - Cartesian Partition Manager
 *
 * Copyright (C) 2012 University of Tokyo.
 *
 */

/**
 * @file   cpm_Base.h
 * CPMのベースクラスのヘッダーファイル
 * @author University of Tokyo
 * @date   2012/05/31
 */

#ifndef _CPM_BASE_H_
#define _CPM_BASE_H_


#if defined(_WIN32) || defined(WIN32)
#define CPM_WINDOWS
#endif

#include "cpm_Define.h"
#include "cpm_Version.h"

#include <stdio.h>
#include <string>
#include <iostream>

#include <algorithm>
#ifndef CPM_WINDOWS
  #include <sys/time.h>    // for linux
#else
  #include "cpm_windows.h" // for win32
#endif

#ifndef _CPM_NO_INLINE_
  #define CPM_INLINE inline
#else
  #define CPM_INLINE
#endif

/** CPMのベースクラス
 */
class cpm_Base
{
////////////////////////////////////////////////////////////////////////////////
// メンバー関数
////////////////////////////////////////////////////////////////////////////////
public:
  /** NULLのランク番号を取得
   *  @return NULLのランク番号
   */
  CPM_INLINE static int getRankNull()
  {
    return MPI_PROC_NULL;
  }

  /** NULLのランクかどうかを確認
   *  @retval true  NULL
   *  @retval false NULLではない
   */
  CPM_INLINE static bool IsRankNull( int rankNo )
  {
    if( rankNo < 0 ) return true;
    return false;
  }

  /** NULLのMPI_Commを取得
   *  @return NULLのMPI_Comm
   */
  CPM_INLINE static MPI_Comm getCommNull()
  {
    return MPI_COMM_NULL;
  }

  /** NULLのMPI_Commかどうかを確認
   *  @retval true  NULL
   *  @retval false NULLではない
   */
  CPM_INLINE static bool IsCommNull( MPI_Comm comm )
  {
    if( comm == MPI_COMM_NULL ) return true;
    return false;
  }

  /** 実数型REAL_TYPEが倍精度かどうか確認
   *  @retval true  倍精度
   *  @retval false 単精度
   */
  CPM_INLINE static bool RealIsDouble()
  {
    if     ( sizeof(REAL_TYPE)==4 ) return false;
    else if( sizeof(REAL_TYPE)==8 ) return true;
    return false;
  }

	/** 時刻の取得(gettimeofday版)
   *  @retrun 時刻
   */
  CPM_INLINE
  static double GetTime()
  {
    timeval tv;
    gettimeofday(&tv, NULL);
    double t = 0.0;
    t += (double)tv.tv_sec;
    t += (double)tv.tv_usec * 1.0e-6;
    return t;
  }

  /** 経過時刻の取得(gettimeofday版)
   *  @param[in] before 計測開始時刻
   *  @return 計測開始時刻からの経過時刻
   */
  CPM_INLINE
  static double GetSpanTime(double before)
  {
    return GetTime() - before;
  }

	/** 時刻の取得(MPI_Wtime版)
   *  @retrun 時刻
   */
  CPM_INLINE
  static double GetWTime()
  {
    return MPI_Wtime();
  }

  /** 経過時刻の取得(MPI_Wtime版)
   *  @param[in] before 計測開始時刻
   *  @return 計測開始時刻からの経過時刻
   */
  CPM_INLINE
  static double GetWSpanTime(double before)
  {
    return GetWTime() - before;
  }

  /** メモリ量の文字列を返す
   *  @param[in] mem メモリ量(byte)
   *  @return メモリ量の文字列
   */
  CPM_INLINE
  static std::string GetMemString(size_t mem)
  {
    const char uu[][8] = { "[Bytes]", "[KB]", "[MB]", "[GB]", "[TB]", "[PB]" };
    double dmem = double(mem);
    int cnt = 0;
    while( dmem >= 1000.0 )
    {
      dmem /= 1024.0;
      cnt++;
    }

    char buf[10];
    sprintf( buf, "%.2f ", dmem );

    std::string str = buf;
    if( cnt <= 5 )
    {
      str += uu[cnt];
    }
    else
    {
      dmem = double(mem);
      sprintf( buf, "%.2f ", dmem );
      str = buf;
      str += uu[0];
    }

    return str;
  }

  /** バージョンを出力する
   *  @param ofs 出力ストリーム
   */
  CPM_INLINE
  static void VersionInfo()
  {
    VersionInfo( std::cout );
  }

  /** バージョンを出力する
   *  @param ofs 出力ストリーム
   */
  CPM_INLINE
  static void VersionInfo( std::ostream &ofs )
  {
    ofs << std::endl
        << " CPM - Cartesian Partition Manager " << std::endl
        << " version " << CPM_VERSION_NO << std::endl
        << std::endl;
  }

  /** 文字列の比較
   *  @param[in] str1       文字列1
   *  @param[in] str2       文字列2
   *  @param[in] ignorecase true=大文字小文字を区別しない、false=区別する
   *  @retval 0     一致する
   *  @retval 0以外 一致しない
   */
  CPM_INLINE
  int cpm_strCompare( std::string str1, std::string str2, bool ignorecase=true )
  {
    // 小文字に変換
    if( ignorecase )
    {
      std::transform(str1.begin(), str1.end(), str1.begin(), ::tolower);
      std::transform(str2.begin(), str2.end(), str2.begin(), ::tolower);
    }

    // 比較
    return str1.compare(str2);
  }

  /** 文字列の比較(文字数指定)
   *  @param[in] str1       文字列1
   *  @param[in] str2       文字列2
   *  @param[in] num        比較する文字数(先頭から)
   *  @param[in] ignorecase true=大文字小文字を区別しない、false=区別する
   *  @retval 0     一致する
   *  @retval 0以外 一致しない
   */
  CPM_INLINE
  int cpm_strCompareN( std::string str1, std::string str2, size_t num, bool ignorecase=true )
  {
    // 部分文字列のコピー
    std::string s1 = str1.substr(0,num);
    std::string s2 = str2.substr(0,num);

    // 比較
    return cpm_strCompare(s1,s2,ignorecase);
  }

protected:
  /** コンストラクタ */
  cpm_Base(){};

  /** デストラクタ */
  virtual ~cpm_Base(){};


////////////////////////////////////////////////////////////////////////////////
// メンバー変数
////////////////////////////////////////////////////////////////////////////////
public:
protected:



};

#endif /* _CPM_BASE_H_ */

