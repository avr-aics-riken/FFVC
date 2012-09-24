#ifndef _FB_UTY_H_
#define _FB_UTY_H_

// #################################################################
//
// CAERU Library
//
// Copyright (c) 2012 All right reserved.
//
// Institute of Industrial Science, The University of Tokyo, Japan. 
//
// #################################################################

/** 
 * @file   FBUtility.h
 * @brief  FlowBase FBUtility class Header
 * @author kero
 */

#include <math.h>
#include <string>
#include <string.h>
#include <stdio.h>
#include <algorithm>

#include "FB_Define.h"
#include "cpm_Define.h"

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
  
  
  /**
   * @brief dirの方向ラベルを返す
   * @param [in] dir 方向コード
   * @return 方向ラベル
   */
  static string getDirection(const int dir)
  {
    string face;
    if      (dir == X_MINUS) face = "X-";
    else if (dir == X_PLUS)  face = "X+";
    else if (dir == Y_MINUS) face = "Y-";
    else if (dir == Y_PLUS)  face = "Y+";
    else if (dir == Z_MINUS) face = "Z-";
    else if (dir == Z_PLUS)  face = "Z+";
    return face;
  }
  
  
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
  static void printVersion(FILE* fp, const string str, const int ver)
  {
    int a, b, c;
    a = b = c = 0;
    
    a = ver / 100;
    b = (ver - a*100) / 10;
    c = ver - a*100 - b*10;
    
    fprintf(fp, "\n\t%s \tVersion %d.%d.%d\n", str.c_str(), a, b, c);
  }

  
	/**
   * @brief 無次元温度varを有次元(Kelvin)にして返す
   * @param [in] var   無次元温度
   * @param [in] base  Control::BaseTemp
   * @param [in] diff  Control::DiffTemp
   */
  static REAL_TYPE convND2Kelvin(const REAL_TYPE var, const REAL_TYPE base, const REAL_TYPE diff)
  {
    return ( base + diff*var );
  }
  
  
  /**
   * @brief 有次元温度varを無次元にして返す
   * @param [in] var  有次元温度(Kelvin or Celsius)
   * @param [in] base Control::BaseTemp
   * @param [in] diff Control::DiffTemp
   * @param [in] Unit 温度の単位
   */
  static REAL_TYPE convD2ND(const REAL_TYPE var, const REAL_TYPE base, const REAL_TYPE diff, const int Unit) 
  {
    REAL_TYPE tmp = convTemp2K(var, Unit);
    return ( (tmp - base) / (REAL_TYPE)fabs(diff) );
  }
  
  
  /**
   * @brief 有次元温度var(Kelvin)を無次元にして返す
   * @param [in] var  有次元温度(Kelvin)
   * @param [in] base Control::BaseTemp
   * @param [in] diff Control::DiffTemp
   */
  static REAL_TYPE convK2ND(const REAL_TYPE var, const REAL_TYPE base, const REAL_TYPE diff) 
  {
    return ( (var - base) / (REAL_TYPE)fabs(diff) );
  }

  
  /**
   * @brief 有次元の温度varを指定された温度単位にして返す
   * @param [in] var  有次元温度(Kelvin)
   * @param [in] Unit 温度の単位
   */
  static REAL_TYPE convK2Temp(const REAL_TYPE var, const int Unit) 
  {
    return ( (Unit==Unit_KELVIN) ? var : var-KELVIN );
  }
  
  
  /**
   * @brief 有次元の温度varを(Kelvin)にして返す
   * @param [in] var  有次元温度(Kelvin or Celsius)
   * @param [in] Unit 温度の単位
   */
  static REAL_TYPE convTemp2K(const REAL_TYPE var, const int Unit)
  {
    return ( (Unit==Unit_KELVIN) ? var : var+KELVIN );
  }
  
  
  /**
   * @brief 有次元速度を無次元にして返す
   * @retval 無次元速度
   * @param [in] var  有次元速度
   * @param [in] refv 代表速度
   */
  static REAL_TYPE convD2ND_V(const REAL_TYPE var, const REAL_TYPE RefV) 
  {
    return ( var / RefV );
  }
  
  
  /**
   * @brief 発熱量(W/m^3)を無次元にして返す
   * @param [in] var   有次元発熱量(W/m^3)
   * @param [in] RefV  代表速度
   * @param [in] RefL  代表長さ
   * @param [in] diff  代表温度差
   * @param [in] rho   媒質密度
   * @param [in] C     媒質比熱
   */
  static REAL_TYPE convD2ND_Hsrc(const REAL_TYPE var, 
                                 const REAL_TYPE RefV, 
                                 const REAL_TYPE RefL, 
                                 const REAL_TYPE diff, 
                                 const REAL_TYPE rho, 
                                 const REAL_TYPE C) 
  {
    return ( var*RefL / (RefV*diff*rho*C) );
  }
  
  
  /**
   * @brief 発熱量を有次元(W/m^3)にして返す
   * @param [in] var  無次元発熱量
   * @param [in] RefV 代表速度
   * @param [in] RefL 代表長さ
   * @param [in] diff 代表温度差
   * @param [in] rho  媒質密度
   * @param [in] C    媒質比熱
   */
  static REAL_TYPE convND2D_Hsrc(const REAL_TYPE var, 
                                 const REAL_TYPE RefV, 
                                 const REAL_TYPE RefL, 
                                 const REAL_TYPE diff, 
                                 const REAL_TYPE rho, 
                                 const REAL_TYPE C) 
  {
    return ( var* RefV*diff*rho*C / RefL );
  }
  
  /**
   * @brief 圧力を無次元にして返す
   * @param [in] var   有次元圧力(absolute or gauge)
   * @param [in] bp    基準圧力
   * @param [in] rho   媒質密度
   * @param [in] RefV  代表速度
   * @param [in] mode  (absolute or gauge)
   */
  static REAL_TYPE convD2ND_P(const REAL_TYPE var, 
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
  static REAL_TYPE convND2D_P(const REAL_TYPE var, 
                              const REAL_TYPE bp, 
                              const REAL_TYPE rho, 
                              const REAL_TYPE RefV, 
                              const int mode) 
  {
    const REAL_TYPE a = var * (RefV*RefV*rho);
    return ( (mode==Unit_Absolute) ? bp+a : a );
  }
  
  
  /**
   * @brief スカラー倍コピー
   * @param [out]    dst   出力
   * @param [in]     size  配列サイズ
   * @param [in]     guide ガイドセルサイズ
   * @param [in]     src   入力
   * @param [in]     scale スカラー倍数
   * @param [in]     mode  スカラー or ベクトル
   * @param [in,out] flop  浮動小数点演算
   */
  template<typename T>
  void xcopy (T* dst, const int* size, const int guide, const T* src, const T scale, const int mode, double& flop);
  
  /**
   * @brief 初期化
   * @param [out]    dst   出力
   * @param [in]     size  配列サイズ
   * @param [in]     guide ガイドセルサイズ
   * @param [in]     init  定数
   * @param [in]     mode  スカラー or ベクトル
   */
  template<typename T>
  void xset (T* dst, const int* size, const int guide, const T init, const int mode);
  
  
  /**
   * @brief ベクトルの初期化（内部のみ）
   * @param [out]    dst   出力
   * @param [in]     size  配列サイズ
   * @param [in]     guide ガイドセルサイズ
   * @param [in]     init  定数
   */
  template<typename T>
  void xsetv (T* dst, const int* size, const int guide, const T* init);
  
};

#endif // _FB_UTY_H_
