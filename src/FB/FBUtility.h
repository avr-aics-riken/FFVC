#ifndef _FB_UTY_H_
#define _FB_UTY_H_

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
 * @file   FBUtility.h
 * @brief  FlowBase FBUtility class Header
 * @author aics
 */

#include "cpm_Define.h"
#include <math.h>
#include <string>
#include <string.h>
#include <stdio.h>
#include <algorithm>
#include <errno.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <stdlib.h>

#include "FB_Define.h"
#include "Medium.h"

#include "omp.h"

using namespace std;

// #################################################################
class FBUtility {

public:
  /** コンストラクタ */
  FBUtility() {}
  
  /** デストラクタ */
  ~FBUtility() {}
 
  
public:
  
  /** 文字列を小文字にして比較
   * @param [in] str1  比較string
   * @param [in] str2  比較string
   * @return true-同じ / false-異なる
   */
  static bool compare(const string str1, const string str2)
  {
    if ( !strcasecmp(str1.c_str(), str2.c_str()) ) 
    {
      return true;
    }
    return false;
  }
  
  
  // ディレクトリがなければ作成、既存なら何もしない（単一ディレクトリ）
  static int c_mkdir(const char* path);

  
  /**
   * @brief ファイル出力時，発散値を計算する
   * @param [out]    dst  単位変換後のデータ
   * @param [in]     src  単位変換前のデータ
   * @param [in]     sz   分割数
   * @param [in]     gc   ガイドセル数
   * @param [in]     coef 係数
   */
  void cnv_Div(REAL_TYPE* dst,
               REAL_TYPE* src,
               int* sz,
               int gc,
               REAL_TYPE coef);
  
  
  /** 
   * @brief 無次元内部エネルギーから有次元/無次元温度への変換
   * @param [out]    dst      温度
   * @param [in]     size     配列長
   * @param [in]     guide    ガイドセル
   * @param [in]     src      無次元内部エネルギー
   * @param [in]     bd       BCindex B
   * @param [in]     mtbl     rho, cp, lambdaの無次元値配列 mtbl[3*(NoCompo+1)]
   * @param [in]     Base_tmp 基準温度(C)
   * @param [in]     Diff_tmp 代表温度差(C)
   * @param [in]     mode     温度の次元指定（true=dimensional, false=non-dimensional)
   * @param [in,out] flop     浮動小数演算数
   */
  void convArrayIE2Tmp(REAL_TYPE* dst,
                       const int* size,
                       const int guide,
                       const REAL_TYPE* src,
                       const int* bd,
                       const REAL_TYPE* mtbl,
                       const REAL_TYPE Base_tmp,
                       const REAL_TYPE Diff_tmp,
                       const bool mode,
                       double& flop);
  
  
  /**
   * @brief 圧力値を有次元から無次元へ変換し，scale倍
   * @param [out]    dst      圧力
   * @param [in]     size     配列長
   * @param [in]     guide    ガイドセル
   * @param [in]     Base_prs 基準圧力(Pa) 基準圧がゼロのとき，ゲージ圧
   * @param [in]     Ref_rho  代表密度(kg/m^3)
   * @param [in]     Ref_v    代表速度(m/s)
   * @param [in,out] flop     浮動小数演算数
   */
  void convArrayPrsD2ND(REAL_TYPE* dst,
                        const int* size,
                        const int guide,
                        const REAL_TYPE Base_prs,
                        const REAL_TYPE Ref_rho,
                        const REAL_TYPE Ref_v,
                        double& flop);
  
  
  /**
   * @brief 圧力値を無次元から有次元へ変換し，scale倍
   * @param [out]    dst      有次元圧力
   * @param [in]     size     配列長
   * @param [in]     guide    ガイドセル
   * @param [in]     src      無次元圧力
   * @param [in]     Base_prs 基準圧力(Pa) 基準圧がゼロのとき，ゲージ圧
   * @param [in]     Ref_rho  代表密度(kg/m^3)
   * @param [in]     Ref_v    代表速度(m/s)
   * @param [in,out] flop     浮動小数演算数
   */
  void convArrayPrsND2D(REAL_TYPE* dst,
                        const int* size,
                        const int guide,
                        const REAL_TYPE* src,
                        const REAL_TYPE Base_prs,
                        const REAL_TYPE Ref_rho,
                        const REAL_TYPE Ref_v,
                        double& flop);
  
  
  /** 
   * @brief 有次元/無次元温度から無次元内部エネルギーへの変換
   * @param [out]    dst      無次元内部エネルギー
   * @param [in]     size     配列長
   * @param [in]     guide    ガイドセル
   * @param [in]     src      温度
   * @param [in]     bd       BCindex B
   * @param [in]     mtbl     rho, cp, lambdaの無次元値配列 mtbl[3*(NoCompo+1)]
   * @param [in]     Base_tmp 基準温度(C)
   * @param [in]     Diff_tmp 代表温度差(C)
   * @param [in]     mode     温度の次元指定（true=dimensional, false=non-dimensional)
   * @param [in,out] flop     浮動小数演算数
   */
  void convArrayTmp2IE(REAL_TYPE* dst,
                       const int* size,
                       const int guide,
                       REAL_TYPE* src,
                       const int* bd,
                       const REAL_TYPE* mtbl,
                       const REAL_TYPE Base_tmp,
                       const REAL_TYPE Diff_tmp,
                       const bool mode,
                       double& flop);
  
  
  /**
   * @brief 全圧データについて，無次元から有次元単位に変換する
   * @param [out]    dst     単位変換後のデータ
   * @param [in]     src     単位変換前のデータ
   * @param [in]     sz      分割数
   * @param [in]     gc      ガイドセル数
   * @param [in]     Ref_rho 代表密度(kg/m^3)
   * @param [in]     Ref_v   代表速度(m/s)
   */
  void convArrayTpND2D(REAL_TYPE* dst,
                       REAL_TYPE* src,
                       int* sz,
                       int gc,
                       const REAL_TYPE Ref_rho,
                       const REAL_TYPE Ref_v);
  
  
  /**
   * @brief 発熱量(W/m^3)を無次元にして返す
   * @param [in] var   有次元発熱量(W/m^3)
   * @param [in] RefV  代表速度
   * @param [in] RefL  代表長さ
   * @param [in] diff  代表温度差
   * @param [in] rho   媒質密度
   * @param [in] C     媒質比熱
   */
  static REAL_TYPE convHsrcD2ND(const REAL_TYPE var,
                                const REAL_TYPE RefV,
                                const REAL_TYPE RefL,
                                const REAL_TYPE diff,
                                const REAL_TYPE rho,
                                const REAL_TYPE C)
  {
    return ( var*RefL / (RefV*diff*rho*C) );
  }
  
  
  /**
   * @brief 圧力を無次元にして返す
   * @param [in] var   有次元圧力(absolute or gauge)
   * @param [in] bp    基準圧力
   * @param [in] rho   媒質密度
   * @param [in] RefV  代表速度
   * @param [in] mode  (absolute or gauge)
   */
  static REAL_TYPE convPrsD2ND(const REAL_TYPE var,
                               const REAL_TYPE bp,
                               const REAL_TYPE rho,
                               const REAL_TYPE RefV,
                               const int mode)
  {
    const REAL_TYPE a = (mode==Unit_Absolute) ? (var-bp) : var;
    return (  a / (RefV*RefV*rho) );
  }
  
  
  /**
   * @brief 圧力を有次元(absolute or gauge)にして返す
   * @param [in] var   無次元圧力
   * @param [in] bp    基準圧力
   * @param [in] rho   媒質密度
   * @param [in] RefV  代表速度
   * @param [in] mode  (absolute or gauge)
   */
  static REAL_TYPE convPrsND2D(const REAL_TYPE var,
                               const REAL_TYPE bp,
                               const REAL_TYPE rho,
                               const REAL_TYPE RefV,
                               const int mode)
  {
    const REAL_TYPE a = var * (RefV*RefV*rho);
    return ( (mode==Unit_Absolute) ? bp+a : a );
  }
  
  
  /**
   * @brief 有次元温度var(Celsius)を無次元にして返す
   * @param [in] var  有次元温度(Celsius)
   * @param [in] base Control::BaseTemp
   * @param [in] diff Control::DiffTemp
   */
  static REAL_TYPE convTempD2ND(const REAL_TYPE var, const REAL_TYPE base, const REAL_TYPE diff) 
  {
    return ( (var - base) / (REAL_TYPE)fabs(diff) );
  }

  
  /**
   * @brief 無次元温度varを有次元(Celsius)にして返す
   * @param [in] var   無次元温度
   * @param [in] base  Control::BaseTemp
   * @param [in] diff  Control::DiffTemp
   */
  static REAL_TYPE convTempND2D(const REAL_TYPE var, const REAL_TYPE base, const REAL_TYPE diff)
  {
    return ( base + diff*var );
  }
  
  
  /**
   * @brief 有次元速度を無次元にして返す
   * @retval 無次元速度
   * @param [in] var  有次元速度
   * @param [in] refv 代表速度
   */
  static REAL_TYPE convVelD2ND(const REAL_TYPE var, const REAL_TYPE RefV) 
  {
    return ( var / RefV );
  }
  
  
  /**
   * @brief S3D配列のコピー
   * @param [out]    dst   出力
   * @param [in]     size  配列サイズ
   * @param [in]     guide ガイドセルサイズ
   * @param [in]     src   入力
   * @param [in]     scale スカラー倍数
   */
  void copyS3D (REAL_TYPE* dst, const int* size, const int guide, const REAL_TYPE* src, const REAL_TYPE scale);
  
  /**
   * @brief V3D配列のコピー
   * @param [out]    dst   出力
   * @param [in]     size  配列サイズ
   * @param [in]     guide ガイドセルサイズ
   * @param [in]     src   入力
   * @param [in]     scale スカラー倍数
   */
  void copyV3D (REAL_TYPE* dst, const int* size, const int guide, const REAL_TYPE* src, const REAL_TYPE scale);
  
  
  /**
   * @brief dirの方向ラベルを返す
   * @param [in] dir 方向コード
   * @return 方向ラベル
   */
  static string getDirection(const int dir)
  {
    string face;
    if      (dir == X_minus) face = "X-";
    else if (dir == X_plus)  face = "X+";
    else if (dir == Y_minus) face = "Y-";
    else if (dir == Y_plus)  face = "Y+";
    else if (dir == Z_minus) face = "Z-";
    else if (dir == Z_plus)  face = "Z+";
    return face;
  }
  
  
  /**
   * @brief S3D配列の初期化
   * @param [out]    dst   出力
   * @param [in]     size  配列サイズ
   * @param [in]     guide ガイドセルサイズ
   * @param [in]     init  定数
   */
  void initS3D (REAL_TYPE* dst, const int* size, const int guide, const REAL_TYPE init);
  
  
  /**
   * @brief S4DEX配列の初期化
   * @param [out]    dst   出力
   * @param [in]     size  配列サイズ
   * @param [in]     guide ガイドセルサイズ
   * @param [in]     init  定数
   */
  void initS4DEX (REAL_TYPE* dst, const int* size, const int guide, const REAL_TYPE init);
  
  
  // 階層ディレクトリの作成
  static int mkdirs(string path);
  
  
  /**
   * @brief メモリ使用量を表示する
   * @param [in] mode     処理モード
   * @param [in] Memory   必要メモリ量
   * @param [in] l_memory local
   * @param [in] fp       ファイルポインタ
   */
  static void MemoryRequirement(const char* mode, const double Memory, const double l_memory, FILE* fp);
  
  
  /** バージョン情報の表示
   * @param [in] fp   ファイルポインタ
   * @param [in] str  名称
   * @param [in] ver  バージョン番号
   */
  static void printVersion(FILE* fp, const string str, const string ver)
  {
    fprintf(fp, "\n\t%s \tVersion %s\n", str.c_str(), ver.c_str());
  }
  
};

#endif // _FB_UTY_H_
