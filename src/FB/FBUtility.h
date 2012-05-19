#ifndef _SKL_FB_UTY_H_
#define _SKL_FB_UTY_H_

/*
 * SPHERE - Skeleton for PHysical and Engineering REsearch
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2004-2012
 *
 */

//@file FBUtility.h
//@brief FlowBase FBUtility class Header
//@author keno, FSI Team, VCAD, RIKEN

#include <math.h>
#include <string>
#include <string.h>
#include <stdio.h>
#include "FB_Define.h"

class FBUtility {

public:
  FBUtility() {}
  ~FBUtility() {}
  
protected:
  static void MemoryRequirement(const char* mode, const unsigned long Memory, const unsigned long l_memory, FILE* fp);
  
public:  
  static std::string getDirection(const unsigned dir);
  
  static void displayMemory (const char* mode, const unsigned long Memory, const unsigned long l_memory, FILE* fp, FILE* mp);
  static void printVersion  (FILE* fp, const char* str, const unsigned ver);

	/**
   @fn static REAL_TYPE convND2Kelvin(const REAL_TYPE var, const REAL_TYPE base, const REAL_TYPE diff)
   @brief 無次元温度varを有次元(Kelvin)にして返す
   @param var 無次元温度
   @param base Control::BaseTemp
   @param diff Control::DiffTemp
   */
  static REAL_TYPE convND2Kelvin(const REAL_TYPE var, const REAL_TYPE base, const REAL_TYPE diff) {
    return ( base + diff*var );
  }
  
  /**
   @fn static REAL_TYPE convD2ND(const REAL_TYPE var, const REAL_TYPE base, const REAL_TYPE diff, const unsigned Unit)
   @brief 有次元温度varを無次元にして返す
   @param var 有次元温度(Kelvin or Celsius)
   @param base Control::BaseTemp
   @param diff Control::DiffTemp
   @param Unit 温度の単位
   */
  static REAL_TYPE convD2ND(const REAL_TYPE var, const REAL_TYPE base, const REAL_TYPE diff, const unsigned Unit) {
    REAL_TYPE tmp = convTemp2K(var, Unit);
    return ( (tmp - base) / (REAL_TYPE)fabs(diff) );
  }
  
  /**
   @fn static REAL_TYPE convK2ND(const REAL_TYPE var, const REAL_TYPE base, const REAL_TYPE diff)
   @brief 有次元温度var(Kelvin)を無次元にして返す
   @param var 有次元温度(Kelvin)
   @param base Control::BaseTemp
   @param diff Control::DiffTemp
   */
  static REAL_TYPE convK2ND(const REAL_TYPE var, const REAL_TYPE base, const REAL_TYPE diff) {
    return ( (var - base) / (REAL_TYPE)fabs(diff) );
  }

  /**
   @fn static inline REAL_TYPE convK2Temp(const REAL_TYPE var, const unsigned Unit)
   @brief 有次元の温度varを指定された温度単位にして返す
   @param var 有次元温度(Kelvin)
   @param Unit 温度の単位
   */
  static REAL_TYPE convK2Temp(const REAL_TYPE var, const unsigned Unit) {
    return ( (Unit==Unit_KELVIN) ? var : var-KELVIN );
  }
  
  /**
   @fn static REAL_TYPE convTemp2K(const REAL_TYPE var, const unsigned Unit)
   @brief 有次元の温度varを(Kelvin)にして返す
   @param var 有次元温度(Kelvin or Celsius)
   @param Unit 温度の単位
   */
  static REAL_TYPE convTemp2K(const REAL_TYPE var, const unsigned Unit) {
    return ( (Unit==Unit_KELVIN) ? var : var+KELVIN );
  }
  
  /**
   @fn static REAL_TYPE convD2ND_V(const REAL_TYPE var, const REAL_TYPE refv)
   @brief 有次元速度を無次元にして返す
   @retval 無次元速度
   @param var 有次元速度
   @param refv 代表速度
   */
  static REAL_TYPE convD2ND_V(const REAL_TYPE var, const REAL_TYPE RefV) {
    return ( var / RefV );
  }
  
  /**
   @fn static REAL_TYPE convD2ND_Hsrc(REAL_TYPE var, REAL_TYPE RefV, REAL_TYPE RefL, REAL_TYPE diff, REAL_TYPE rho, REAL_TYPE C)
   @brief 発熱量(W/m^3)を無次元にして返す
   @param var 有次元発熱量(W/m^3)
   @param RefV 代表速度
   @param RefL 代表長さ
   @param diff 代表温度差
   @param rho 媒質密度
   @param C 媒質比熱
   */
  static REAL_TYPE convD2ND_Hsrc(const REAL_TYPE var, 
                                       const REAL_TYPE RefV, 
                                       const REAL_TYPE RefL, 
                                       const REAL_TYPE diff, 
                                       const REAL_TYPE rho, 
                                       const REAL_TYPE C) {
    return ( var*RefL / (RefV*diff*rho*C) );
  }
  
  /**
   @fn static REAL_TYPE convND2D_Hsrc(const REAL_TYPE var, const REAL_TYPE RefV, const REAL_TYPE RefL, const REAL_TYPE diff, const REAL_TYPE rho, const REAL_TYPE C)
   @brief 発熱量を有次元(W/m^3)にして返す
   @param var 無次元発熱量
   @param RefV 代表速度
   @param RefL 代表長さ
   @param diff 代表温度差
   @param rho 媒質密度
   @param C 媒質比熱
   */
  static REAL_TYPE convND2D_Hsrc(const REAL_TYPE var, 
                                       const REAL_TYPE RefV, 
                                       const REAL_TYPE RefL, 
                                       const REAL_TYPE diff, 
                                       const REAL_TYPE rho, 
                                       const REAL_TYPE C) {
    return ( var* RefV*diff*rho*C / RefL );
  }
  
  /**
   @fn static REAL_TYPE convD2ND_P(const REAL_TYPE var, const REAL_TYPE bp, const REAL_TYPE rho, const REAL_TYPE RefV, const unsigned mode)
   @brief 圧力を無次元にして返す
   @param var 有次元圧力(absolute or gauge)
   @param bp 基準圧力
   @param rho 媒質密度
   @param RefV 代表速度
   @param mode (absolute or gauge)
   */
  static REAL_TYPE convD2ND_P(const REAL_TYPE var, 
                                    const REAL_TYPE bp, 
                                    const REAL_TYPE rho, 
                                    const REAL_TYPE RefV, 
                                    const unsigned mode) {
    const REAL_TYPE a = (mode==Unit_Absolute) ? (var-bp) : var;
    return (  a / (RefV*RefV*rho) );
  }
  
  /**
   @fn static REAL_TYPE convND2D_P(const REAL_TYPE var, const REAL_TYPE bp, const REAL_TYPE rho, const REAL_TYPE RefV)
   @brief 圧力を有次元(absolute or gauge)にして返す
   @param var 無次元圧力
   @param bp 基準圧力
   @param rho 媒質密度
   @param RefV 代表速度
   @param mode (absolute or gauge)
   */
  static REAL_TYPE convND2D_P(const REAL_TYPE var, const REAL_TYPE bp, const REAL_TYPE rho, const REAL_TYPE RefV, const unsigned mode) {
    const REAL_TYPE a = var * (RefV*RefV*rho);
    return ( (mode==Unit_Absolute) ? bp+a : a );
  }
  
  /**
   @fn static inline unsigned getFindexS3D(const unsigned* sz, unsigned gc, int i, int j, int k)
   @brief Fortranの3次元インデックスから1次元インデックスを取得する
   @param sz    I,J,K方向サイズ（ガイドセルを含まない）
   @param gc    ガイドセル
   @param i     I方向インデックス（ガイドセルを含まない）
   @param j     J方向インデックス（ガイドセルを含まない）
   @param k     K方向インデックス（ガイドセルを含まない）
   @return      1次元インデックス
   */
  static inline unsigned getFindexS3D(const unsigned* sz, unsigned gc, int i, int j, int k) {
    //return ( (sz[0]+gc*2)*(sz[1]+gc*2)*(k+gc-1) + (sz[0]+gc*2)*(j+gc-1) + i+gc-1 );
    int t1 = gc*2;
    int t2 = gc-1;
    int t3 = sz[0]+t1;
    return ( t3*( (sz[1]+t1)*(k+t2) + j+t2 ) + i+t2 );
  }
  
  /**
   @fn static inline unsigned getFindexV3DEx(const unsigned* sz, unsigned gc, int l, int i, int j, int k)
   @brief Fortranの3次元Exベクトルインデックスから1次元インデックスを取得する
   @param sz    I,J,K方向サイズ（ガイドセルを含まない）
   @param gc    ガイドセル
   @param l     ベクトルインデックス
   @param i     I方向インデックス（ガイドセルを含まない）
   @param j     J方向インデックス（ガイドセルを含まない）
   @param k     K方向インデックス（ガイドセルを含まない）
   @return      1次元インデックス
   */
  static inline unsigned getFindexV3DEx(const unsigned* sz, unsigned gc, int l, int i, int j, int k) {
    int t1 = gc*2;
    int t2 = gc-1;
    int t3 = sz[0]+t1;
    return ( 3*(t3*(sz[1]+t1)*(k+t2) + t3*(j+t2) + i+t2) + l );
    //return ( 3*(sz[0]+gc*2)*(sz[1]+gc*2)*(k+gc-1) + 3*(sz[0]+gc*2)*(j+gc-1) + 3*(i+gc-1) + l );
  }
  
  /**
   @fn static inline unsigned getFindexS3Dcut(const unsigned* sz, unsigned gc, int l, int i, int j, int k)
   @brief Fortranの3次元cut用インデックスから1次元インデックスを取得する
   @param sz    I,J,K方向サイズ（ガイドセルを含まない）
   @param gc    ガイドセル
   @param l     方向インデックス[0,5]
   @param i     I方向インデックス（ガイドセルを含まない）
   @param j     J方向インデックス（ガイドセルを含まない）
   @param k     K方向インデックス（ガイドセルを含まない）
   @return      1次元インデックス
   */
  static inline unsigned getFindexS3Dcut(const unsigned* sz, unsigned gc, int l, int i, int j, int k) {
    int t1 = gc*2;
    int t2 = gc-1;
    int t3 = sz[0]+t1;
    return ( 6*(t3*(sz[1]+t1)*(k+t2) + t3*(j+t2) + i+t2) + l );
  }
  
  /**
   @fn static inline unsigned getFindexBID8(const unsigned* sz, unsigned gc, int l, int i, int j, int k)
   @brief Fortranの3次元cut用Bid8インデックスから1次元インデックスを取得する
   @param sz    I,J,K方向サイズ（ガイドセルを含まない）
   @param gc    ガイドセル
   @param i     I方向インデックス（ガイドセルを含まない）
   @param j     J方向インデックス（ガイドセルを含まない）
   @param k     K方向インデックス（ガイドセルを含まない）
   @return      1次元インデックス
   */
  static inline unsigned getFindexBID8(const unsigned* sz, unsigned gc, int i, int j, int k) {
    int t1 = gc*2;
    int t2 = gc-1;
    int t3 = sz[0]+t1;
    return ( 2*(t3*(sz[1]+t1)*(k+t2) + t3*(j+t2) + i+t2) );
  }
  
};

#endif // _SKL_FB_UTY_H_
