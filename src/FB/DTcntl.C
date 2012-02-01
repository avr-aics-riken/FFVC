/*
 * SPHERE - Skeleton for PHysical and Engineering REsearch
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2004-2010
 *
 */

//@file DTcntl.C
//@brief FlowBase DTcntl class
//@author keno, FSI Team, VCAD, RIKEN

#include "DTcntl.h"

/**
 @fn bool DTcntl::chkDtSelect(void) const
 @brief 時間積分幅とKindOfSolver種別の整合性をチェック
 */
bool DTcntl::chkDtSelect(void) const
{
  switch (KOS) {
    case FLOW_ONLY:
    case THERMAL_FLOW:
    case THERMAL_FLOW_NATURAL:
    case CONJUGATE_HEAT_TRANSFER:
      break;
      
    case SOLID_CONDUCTION:
      if ( (scheme != dt_direct) && (scheme != dt_dfn) ) {
        return false;
      }
      break;
  }
  return true;
}

/**
 @fn bool DTcntl::set_Scheme(const char* str, const double val)
 @brief Δtのスキームを設定する
 @retval 設定の成否
 @param str キーワード
 @param val 値
 @note 現時点では，Δtは一定値のみ
 */
bool DTcntl::set_Scheme(const char* str, const double val)
{
  if ( !str ) return false;

  if ( !strcasecmp(str, "Direct") ) {
    scheme = dt_direct;
  }
  else if ( !strcasecmp(str, "CFL_Reference_Velocity") ) {
    scheme = dt_cfl_ref_v;
  }
  //else if ( !strcasecmp(str, "CFL_MaxV") ) {
  //  scheme = dt_cfl_max_v;
  //}
  else if ( !strcasecmp(str, "Diffusion") ) {
    scheme = dt_dfn;
  }
  else if ( !strcasecmp(str, "CFL_Diffusion_Reference_Velocity") ) {
    scheme = dt_cfl_dfn_ref_v;
  }
  //else if ( !strcasecmp(str, "CFL_DFN_MaxV") ) {
  //  scheme = dt_cfl_dfn_max_v;
  //}
  //else if ( !strcasecmp(str, "CFL_MaxV_CP") ) {
  //  scheme = dt_cfl_max_v_cp;
  //}
  else {
    return false;
  }
  
  CFL = val;
  
  return true;
}

/**
 @fn void DTcntl::set_Vars(const unsigned m_kos, const unsigned m_mode, const double m_dh, const double re, const double pe)
 @brief 基本変数をコピー
 @param m_kos ソルバーの種類
 @param m_mode 次元モード
 @param m_dh 無次元格子幅
 @param re レイノルズ数
 @param pe ペクレ数
 */
void DTcntl::set_Vars(const unsigned m_kos, const unsigned m_mode, const double m_dh, const double re, const double pe)
{
  KOS      = m_kos;
  mode     = m_mode;
  dh       = m_dh;
  Reynolds = re;
  Peclet   = pe; 
}

/**
 @fn double DTcntl::dtCFL(const double Uref) const
 @brief CFL数で指定されるdtを計算
 @param Uref 速度の参照値（無次元）
 @retval 時間積分幅 dt
 */
double DTcntl::dtCFL(const double Uref) const
{
  double v = (Uref < 1.0) ? 1.0 : Uref; // 1.0 is non-dimensional reference velocity
  return (dh*CFL / v);
}

/**
 @fn double DTcntl::dtDFN(const double coef) const
 @brief 拡散数のdt制限
 @param coef 係数 (Reynolds number or Peclet number)
 */
double DTcntl::dtDFN(const double coef) const
{
  return coef * dh*dh/6.0;
}

/**
 @fn unsigned DTcntl::set_DT(const double vRef)
 @brief 各種モードに対応する時間積分幅を設定する
 @retval return code
 @param vRef 速度の参照値（無次元）
 @note deltaTは無次元
 */
unsigned DTcntl::set_DT(const double vRef)
{
  double dtC, dtD, a, b;

  switch ( scheme ) {
    case dt_direct:
      deltaT = CFL;
      break;
      
    case dt_cfl_ref_v:
      if ( KOS == SOLID_CONDUCTION ) return 1;
      deltaT = dtCFL( vRef );
      break;
      
    case dt_cfl_max_v:
      if ( KOS == SOLID_CONDUCTION ) return 2;
      deltaT = dtCFL(vRef);
      break;
      
    case dt_dfn:
      if ( KOS != SOLID_CONDUCTION ) return 3;
      deltaT = dtDFN( Peclet );
      break;
      
    case dt_cfl_dfn_ref_v:
      switch (KOS) {
        case FLOW_ONLY:
          dtC = dtCFL( vRef );
          dtD = dtDFN( Reynolds );
          deltaT = (dtC > dtD) ? dtD : dtC;
          break;
          
        case THERMAL_FLOW:
        case THERMAL_FLOW_NATURAL:
        case CONJUGATE_HEAT_TRANSFER:
          dtC = dtCFL( vRef );
          a = dtDFN( Reynolds );
          b = dtDFN( Peclet );
          dtD = (a > b) ? b : a;
          deltaT = (dtC > dtD) ? dtD : dtC;
          break;
          
        case SOLID_CONDUCTION:
          return 4;
          break;
      }
      break;
      
    case dt_cfl_dfn_max_v:
      return 6;
      break;
      
    case dt_cfl_max_v_cp:
      return 7;
      break;
  }
  
  return 0;
}
