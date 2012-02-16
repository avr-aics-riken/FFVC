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

#include "Skl.h"
#include "SklSolverBase.h"
#include "FBDefine.h"
#include "Component.h"

class FBUtility{

public:
  FBUtility() {}
  ~FBUtility() {}
  
protected:
  static void MemoryRequirement(const char* mode, const unsigned long Memory, const unsigned long l_memory, FILE* fp);
  
public:  
  static string getDirection(const unsigned dir);
  
  static bool shiftVin3D  (SklVector3DEx<SKL_REAL>* dst, const SklVector3DEx<SKL_REAL>* src, SKL_REAL v00[3], unsigned stepAvr=1);
  static bool shiftVout3D (SklVector3DEx<SKL_REAL>* dst, const SklVector3DEx<SKL_REAL>* src, SKL_REAL v00[3], unsigned stepAvr=1);
  
  static void displayMemory (const char* mode, const unsigned long Memory, const unsigned long l_memory, FILE* fp, FILE* mp);
  static void printVersion  (FILE* fp, const char* str, const unsigned ver);

	/**
   @fn static SKL_REAL convND2Kelvin(const SKL_REAL var, const SKL_REAL base, const SKL_REAL diff)
   @brief 無次元温度varを有次元(Kelvin)にして返す
   @param var 無次元温度
   @param base Control::BaseTemp
   @param diff Control::DiffTemp
   */
  static SKL_REAL convND2Kelvin(const SKL_REAL var, const SKL_REAL base, const SKL_REAL diff) {
    return ( base + diff*var );
  }
  
  /**
   @fn static SKL_REAL convD2ND(const SKL_REAL var, const SKL_REAL base, const SKL_REAL diff, const unsigned Unit)
   @brief 有次元温度varを無次元にして返す
   @param var 有次元温度(Kelvin or Celsius)
   @param base Control::BaseTemp
   @param diff Control::DiffTemp
   @param Unit 温度の単位
   */
  static SKL_REAL convD2ND(const SKL_REAL var, const SKL_REAL base, const SKL_REAL diff, const unsigned Unit) {
    SKL_REAL tmp = convTemp2K(var, Unit);
    return ( (tmp - base) / (SKL_REAL)fabs(diff) );
  }
  
  /**
   @fn static SKL_REAL convK2ND(const SKL_REAL var, const SKL_REAL base, const SKL_REAL diff)
   @brief 有次元温度var(Kelvin)を無次元にして返す
   @param var 有次元温度(Kelvin)
   @param base Control::BaseTemp
   @param diff Control::DiffTemp
   */
  static SKL_REAL convK2ND(const SKL_REAL var, const SKL_REAL base, const SKL_REAL diff) {
    return ( (var - base) / (SKL_REAL)fabs(diff) );
  }

  /**
   @fn static inline SKL_REAL convK2Temp(const SKL_REAL var, const unsigned Unit)
   @brief 有次元の温度varを指定された温度単位にして返す
   @param var 有次元温度(Kelvin)
   @param Unit 温度の単位
   */
  static SKL_REAL convK2Temp(const SKL_REAL var, const unsigned Unit) {
    return ( (Unit==CompoList::Unit_KELVIN) ? var : var-KELVIN );
  }
  
  /**
   @fn static SKL_REAL convTemp2K(const SKL_REAL var, const unsigned Unit)
   @brief 有次元の温度varを(Kelvin)にして返す
   @param var 有次元温度(Kelvin or Celsius)
   @param Unit 温度の単位
   */
  static SKL_REAL convTemp2K(const SKL_REAL var, const unsigned Unit) {
    return ( (Unit==CompoList::Unit_KELVIN) ? var : var+KELVIN );
  }
  
  /**
   @fn static SKL_REAL convD2ND_V(const SKL_REAL var, const SKL_REAL refv)
   @brief 有次元速度を無次元にして返す
   @retval 無次元速度
   @param var 有次元速度
   @param refv 代表速度
   */
  static SKL_REAL convD2ND_V(const SKL_REAL var, const SKL_REAL RefV) {
    return ( var / RefV );
  }
  
  /**
   @fn static SKL_REAL convD2ND_Hsrc(SKL_REAL var, SKL_REAL RefV, SKL_REAL RefL, SKL_REAL diff, SKL_REAL rho, SKL_REAL C)
   @brief 発熱量(W/m^3)を無次元にして返す
   @param var 有次元発熱量(W/m^3)
   @param RefV 代表速度
   @param RefL 代表長さ
   @param diff 代表温度差
   @param rho 媒質密度
   @param C 媒質比熱
   */
  static SKL_REAL convD2ND_Hsrc(const SKL_REAL var, 
                                       const SKL_REAL RefV, 
                                       const SKL_REAL RefL, 
                                       const SKL_REAL diff, 
                                       const SKL_REAL rho, 
                                       const SKL_REAL C) {
    return ( var*RefL / (RefV*diff*rho*C) );
  }
  
  /**
   @fn static SKL_REAL convND2D_Hsrc(const SKL_REAL var, const SKL_REAL RefV, const SKL_REAL RefL, const SKL_REAL diff, const SKL_REAL rho, const SKL_REAL C)
   @brief 発熱量を有次元(W/m^3)にして返す
   @param var 無次元発熱量
   @param RefV 代表速度
   @param RefL 代表長さ
   @param diff 代表温度差
   @param rho 媒質密度
   @param C 媒質比熱
   */
  static SKL_REAL convND2D_Hsrc(const SKL_REAL var, 
                                       const SKL_REAL RefV, 
                                       const SKL_REAL RefL, 
                                       const SKL_REAL diff, 
                                       const SKL_REAL rho, 
                                       const SKL_REAL C) {
    return ( var* RefV*diff*rho*C / RefL );
  }
  
  /**
   @fn static SKL_REAL convD2ND_P(const SKL_REAL var, const SKL_REAL bp, const SKL_REAL rho, const SKL_REAL RefV, const unsigned mode)
   @brief 圧力を無次元にして返す
   @param var 有次元圧力(absolute or gauge)
   @param bp 基準圧力
   @param rho 媒質密度
   @param RefV 代表速度
   @param mode (absolute or gauge)
   */
  static SKL_REAL convD2ND_P(const SKL_REAL var, 
                                    const SKL_REAL bp, 
                                    const SKL_REAL rho, 
                                    const SKL_REAL RefV, 
                                    const unsigned mode) {
    const SKL_REAL a = (mode==CompoList::Absolute) ? (var-bp) : var;
    return (  a / (RefV*RefV*rho) );
  }
  
  /**
   @fn static SKL_REAL convND2D_P(const SKL_REAL var, const SKL_REAL bp, const SKL_REAL rho, const SKL_REAL RefV)
   @brief 圧力を有次元(absolute or gauge)にして返す
   @param var 無次元圧力
   @param bp 基準圧力
   @param rho 媒質密度
   @param RefV 代表速度
   @param mode (absolute or gauge)
   */
  static SKL_REAL convND2D_P(const SKL_REAL var, const SKL_REAL bp, const SKL_REAL rho, const SKL_REAL RefV, const unsigned mode) {
    const SKL_REAL a = var * (RefV*RefV*rho);
    return ( (mode==CompoList::Absolute) ? bp+a : a );
  }
  
  /**
   @fn static inline void CalcIndex(cosnt unsigned dst_ilen, const unsigned dst_jlen, const unsigned dst_klen, cosnt unsigned dst_gc, const unsigned src_gc, 
   unsigned& dst_ix, unsigned& dst_jx, unsigned& dst_kx, int& diff, unsigned& sta)
   @brief calculation of index, this method is taken from SklUtil class
   @param dst_ilen dimension size /w guide cell for dst
   @param dst_jlen 
   @param dst_klen 
   @param dst_gc size of guide cell
   @param src_gc 
   @param dst_ix return value
   @param dst_jx 
   @param dst_kx 
   @param diff 
   @param sta 
   */
  static inline void CalcIndex(const unsigned dst_ilen,
                               const unsigned dst_jlen,
                               const unsigned dst_klen,
                               const unsigned dst_gc,
                               const unsigned src_gc,
                               unsigned& dst_ix,
                               unsigned& dst_jx,
                               unsigned& dst_kx,
                               int& diff,
                               unsigned& sta) {
    diff = src_gc - dst_gc;
    if( src_gc >= dst_gc ){
      sta = 0;
      dst_ix = dst_ilen; dst_jx = dst_jlen; dst_kx = dst_klen;
    }
    else {
      sta = abs(diff);
      dst_ix = dst_ilen+diff; 
      dst_jx = dst_jlen+diff; 
      dst_kx = dst_klen+diff;
    }
  }
  
};

#endif // _SKL_FB_UTY_H_
