//##################################################################################
//
// FFV-C : Frontflow / violet Cartesian
//
// Copyright (c) 2007-2011 VCAD System Research Program, RIKEN.
// All rights reserved.
//
// Copyright (c) 2011-2015 Institute of Industrial Science, The University of Tokyo.
// All rights reserved.
//
// Copyright (c) 2012-2015 Advanced Institute for Computational Science, RIKEN.
// All rights reserved.
//
//##################################################################################

/** 
 * @file   ffv_SetBC.C
 * @brief  FFV BC Class
 * @author aics
 */

#include "ffv_SetBC.h"




// #################################################################
// ドライバ指定のチェック
void SetBC3D::checkDriver(FILE* fp)
{
  int o_dir, c_dir;
  int o_pos, c_pos;
  int node_st_i, node_st_j, node_st_k;
  int* st = NULL;
  
  node_st_i = node_st_j = node_st_k = 0;
  
  if( numProc > 1 )
  {
    node_st_i = head[0];
    node_st_j = head[1];
    node_st_k = head[2];
  }
  
  for (int face=0; face<NOFACE; face++)
  {
    if ( (obc[face].getClass() == OBC_PERIODIC) && (obc[face].getPrdcMode() == BoundaryOuter::prdc_Driver) )
    {

      for (int n=1; n<=NoCompo; n++)
      {
        if ( cmp[n].getType() == PERIODIC )
        {
          // 方向のチェック
          o_dir = obc[face].get_DriverDir();
          c_dir = cmp[n].getPeriodicDir();
          
          if ( o_dir != c_dir )
          {
            fprintf(fp, "\tDriver direction is different between OBC[%s] and Component[%s].", 
                    FBUtility::getDirection(o_dir).c_str(),
                    FBUtility::getDirection(c_dir).c_str());
            printf("\tDriver direction is different between OBC[%s] and Component[%s].", 
                    FBUtility::getDirection(o_dir).c_str(),
                    FBUtility::getDirection(c_dir).c_str());
            Exit(0);
          }
          
          st = cmp[n].getBbox_st();
          
          // 位置のチェック  IDは，流入面に対して面直1層の指定
          switch (c_dir)
          {
            case X_minus:
            case X_plus:
              c_pos = node_st_i + st[0];
              break;
              
            case Y_minus:
            case Y_plus:
              c_pos = node_st_j + st[1];
              break;
              
            case Z_minus:
            case Z_plus:
              c_pos = node_st_k + st[2];
              break;
              
            default:
              Exit(0);
              break;
          }
          o_pos = obc[face].get_DriverIndex();
          
          if ( o_pos != c_pos )
          {
            fprintf(fp, "\tDriver Lid Position is different between OBC[%d] and Component[%d].", o_pos, c_pos);
            printf("\tDriver Lid Position is different between OBC[%d] and Component[%d].", o_pos, c_pos);
            Exit(0);
          }
        }
      }
    }
  }  
}

// #################################################################
/**
 * @brief コンポーネントの角速度成分を取り出す
 * @param [in]     n    コンポーネントのインデクス
 * @param [out]    vec  ベクトル成分
 * @param [out]    ctr  中心座標
 * @param [in]     tm   時刻
 * @param [in]     v00  格子速度
 */
void SetBC3D::extractAngularVel(const int n, REAL_TYPE* vec, REAL_TYPE* ctr, const double tm, const REAL_TYPE* v00)
{
  REAL_TYPE omg = cmp[n].ca[0] * RefL / RefV * (REAL_TYPE)tm * v00[0]; // [rad]
  vec[0] = cmp[n].nv[0] * omg; // normal vector x angular velocity
  vec[1] = cmp[n].nv[1] * omg;
  vec[2] = cmp[n].nv[2] * omg;
  ctr[0] = cmp[n].oc[0] / RefL;
  ctr[1] = cmp[n].oc[1] / RefL;
  ctr[2] = cmp[n].oc[2] / RefL;
}



// #################################################################
/**
 * @brief コンポーネントの速度境界条件の成分を取り出す
 * @param [in]     n    コンポーネントのインデクス
 * @param [out]    vec  ベクトル成分
 * @param [in]     tm   無次元時刻
 * @param [in]     v00  無次元格子速度
 * @retval 無次元速度
 */
REAL_TYPE SetBC3D::extractVelLBC(const int n, REAL_TYPE* vec, const double tm, const REAL_TYPE* v00)
{
  REAL_TYPE vel;
  const REAL_TYPE c_pai = (REAL_TYPE)(2.0*asin(1.0));
  bool policy = cmp[n].isPolicy_Massflow(); // true => Massflow
  REAL_TYPE tm2 = (REAL_TYPE)tm * RefL/RefV;
  
  if ( cmp[n].get_V_Profile() == CompoList::vel_zero )
  {
    vel = 0.0;
  }
  else if ( cmp[n].get_V_Profile() == CompoList::vel_constant )
  {
    if ( policy )
    {
      vel = cmp[n].ca[CompoList::bias] / (RefV * RefL * RefL * cmp[n].area) * v00[0];
    }
    else
    {
      vel = cmp[n].ca[CompoList::bias] / RefV * v00[0];
    }
  }
  else if ( cmp[n].get_V_Profile() == CompoList::vel_harmonic )
  {
    if ( policy )
    {
      REAL_TYPE a, b, c;
      a = cmp[n].ca[CompoList::amplitude];
      b = 2.0 * c_pai * cmp[n].ca[CompoList::frequency] * tm2  + cmp[n].ca[CompoList::initphase];
      c = cmp[n].ca[CompoList::bias];
      vel = ( a * sin(b) + c ) / (RefV * RefL * RefL * cmp[n].area) * v00[0];
    }
    else
    {
      REAL_TYPE a, b, c;
      a = cmp[n].ca[CompoList::amplitude];
      b = 2.0 * c_pai * cmp[n].ca[CompoList::frequency] * tm2 + cmp[n].ca[CompoList::initphase];
      c = cmp[n].ca[CompoList::bias];
      vel = ( a * sin(b) + c )  / RefV * v00[0];

      //a = obc[n].ca[CompoList::amplitude]/RefV; // non-dimensional velocity amplitude
      //b = 2.0*c_pai*obc[n].ca[CompoList::frequency]* RefL/RefV * (REAL_TYPE)tm + obc[n].ca[CompoList::initphase];
      //vel = ( a*sin(b) + obc[n].ca[CompoList::bias]/RefV ) * v00[0];
    }
  }
  else if ( cmp[n].get_V_Profile() == CompoList::vel_polynomial6 )
  {
    if ( policy )
    {
      vel = tm2 * ( cmp[n].ca[0]
          + tm2 * ( cmp[n].ca[1]
          + tm2 * ( cmp[n].ca[2]
          + tm2 * ( cmp[n].ca[3]
          + tm2 * ( cmp[n].ca[4]
          + tm2 * ( cmp[n].ca[5] ))))))
       / (RefV * cmp[n].area) * v00[0];
    }
    else
    {
      vel = tm2 * ( cmp[n].ca[0]
          + tm2 * ( cmp[n].ca[1]
          + tm2 * ( cmp[n].ca[2]
          + tm2 * ( cmp[n].ca[3]
          + tm2 * ( cmp[n].ca[4]
          + tm2 * ( cmp[n].ca[5] ))))))
      / RefV * v00[0];
    }
  }
  
  vec[0] = cmp[n].nv[0] * vel;
  vec[1] = cmp[n].nv[1] * vel;
  vec[2] = cmp[n].nv[2] * vel;
  
  return vel;
}


// #################################################################
/**
 * @brief 外部境界の速度境界条件の成分を取り出す
 * @param [in]     n    面のインデクス
 * @param [out]    vec  ベクトル成分
 * @param [in]     tm   時刻
 * @param [in]     v00  格子速度
 * @note OUter BCは速度規定のみ
 */
REAL_TYPE SetBC3D::extractVelOBC(const int n, REAL_TYPE* vec, const double tm, const REAL_TYPE* v00)
{
  REAL_TYPE vel;
  const REAL_TYPE c_pai = (REAL_TYPE)(2.0*asin(1.0));
  REAL_TYPE tm2 = (REAL_TYPE)tm * RefL/RefV;
  
  if ( obc[n].get_V_Profile() == CompoList::vel_zero )
  {
    vel = 0.0;
  }
  else if ( obc[n].get_V_Profile() == CompoList::vel_constant )
  {
    vel = obc[n].ca[CompoList::bias] / RefV * v00[0];
  }
  else if ( obc[n].get_V_Profile() == CompoList::vel_harmonic )
  {
    REAL_TYPE a, b, c;
    a = obc[n].ca[CompoList::amplitude];
    b = 2.0 * c_pai * obc[n].ca[CompoList::frequency] * tm2 + obc[n].ca[CompoList::initphase];
    c = obc[n].ca[CompoList::bias];
    vel = ( a * sin(b) + c )  / RefV * v00[0];
    
    //a = obc[n].ca[CompoList::amplitude]/RefV; // non-dimensional velocity amplitude
    //b = 2.0*c_pai*obc[n].ca[CompoList::frequency]* RefL/RefV * (REAL_TYPE)tm + obc[n].ca[CompoList::initphase];
    //vel = ( a*sin(b) + obc[n].ca[CompoList::bias]/RefV ) * v00[0];
  }
  else if ( obc[n].get_V_Profile() == CompoList::vel_polynomial6 )
  {
    vel = tm2 * ( obc[n].ca[0]
        + tm2 * ( obc[n].ca[1]
        + tm2 * ( obc[n].ca[2]
        + tm2 * ( obc[n].ca[3]
        + tm2 * ( obc[n].ca[4]
        + tm2 * ( obc[n].ca[5] ))))))
    / RefV * v00[0];
  }
  
  vec[0] = obc[n].nv[0] * vel;
  vec[1] = obc[n].nv[1] * vel;
  vec[2] = obc[n].nv[2] * vel;
  
  return vel;
}


// #################################################################
// 拡散部分に関する内部エネルギーの内部境界処理
void SetBC3D::InnerTBCface(REAL_TYPE* d_qbc, const int* d_cdf, const int* d_bcd, REAL_TYPE* d_ie, const REAL_TYPE* d_ie0)
{
  for (int n=1; n<=NoCompo; n++)
  {
    int H = cmp[n].getHtype();
    
    switch ( cmp[n].getType() )
    {
      case HEATFLUX:
        cmp[n].setMonitorValue( psIbcHeatflux(d_qbc, d_cdf, n) );
        break;
        
      case TRANSFER:
        if      ( H == HT_S)  cmp[n].setMonitorValue( psIbcTransferS (d_qbc, d_cdf, d_bcd, n, d_ie0) );
        else if ( H == HT_SF) cmp[n].setMonitorValue( psIbcTransferSF(d_qbc, d_cdf, d_bcd, n, d_ie0) );
        else if ( H == HT_SN) cmp[n].setMonitorValue( psIbcTransferSN(d_qbc, d_cdf, d_bcd, n, d_ie0) );
        break;
        
      case ISOTHERMAL:
        cmp[n].setMonitorValue( psIbcIsoThermal(d_qbc, d_cdf, d_bcd, n, d_ie0) );
        break;
        
      case RADIANT:
        //setRadiant(qbc, bx, n, t0);
        break;
    }
  }
}


// #################################################################
// セルに対する温度の内部境界をセットする
void SetBC3D::InnerTBCvol(REAL_TYPE* d_ie, const int* d_bcd, const REAL_TYPE dt)
{
  for (int n=1; n<=NoCompo; n++)
  {
    switch ( cmp[n].getType() )
    {
      case HEAT_SRC:
        cmp[n].setMonitorValue( psIbcHeatGen(d_ie, d_bcd, n, dt) );
        break;
        
      case CNST_TEMP:
        psIbcConstTemp(d_ie, d_bcd, n);
        break;
    }
  }
}



// #################################################################
/**
 @brief 速度ベクトルの内部周期境界条件処理
 @param d_v 速度ベクトルのデータクラス
 @param d_bd BCindex B
 */
void SetBC3D::InnerVBCperiodic(REAL_TYPE* d_v, int* d_bd)
{
  
  int st[3], ed[3];
  
  for (int n=1; n<=NoCompo; n++)
  {
    cmp[n].getBbox(st, ed);
    
    if ( cmp[n].getType() == PERIODIC )
    {
      Vibc_Prdc(d_v, st, ed, d_bd, n, cmp[n].getPeriodicDir());
    }
  }
}



// #################################################################
/**
 @brief 圧力の内部境界条件処理
 @param d_p 圧力のデータクラス
 @param d_bcd BCindex B
 */
void SetBC3D::InnerPBCperiodic(REAL_TYPE* d_p, int* d_bcd)
{
  int dir;
  int st[3], ed[3];
  REAL_TYPE pv;
  
  for (int n=1; n<=NoCompo; n++)
  {
    cmp[n].getBbox(st, ed);
    dir = cmp[n].getPeriodicDir();
    pv = FBUtility::convPrsD2ND(cmp[n].ca[0], BasePrs, rho_0, RefV, Unit_Prs);
    
    if ( cmp[n].getType() == PERIODIC )
    {
      Pibc_Prdc(d_p, st, ed, d_bcd, n, dir, pv);
    }
  }
}


// #################################################################
// 速度境界条件による速度の発散の修正ほか
void SetBC3D::modDivergence(REAL_TYPE* dv, int* d_cdf, double tm, Control* C, REAL_TYPE* v00, Gemini_R* avr, double& flop)
{
  REAL_TYPE vec[3], dummy, ctr[3];
  int st[3], ed[3];
  int typ=0;
  int gd = guide;
  double fcount = 0.0;
  
  REAL_TYPE aa[2]; // バッファ
  
  // 内部境界条件による修正
  for (int n=1; n<=NoCompo; n++)
  {
    typ = cmp[n].getType();
    
    cmp[n].getBbox(st, ed);
    
    switch (typ)
    {
      case OUTFLOW:
        div_ibc_oflow_vec_(dv, size, &gd, pitch, st, ed, d_cdf, &n, aa, &fcount);
        avr[n].p0 = aa[0]; // 積算速度
        avr[n].p1 = aa[1]; // 積算回数
        if ( aa[1] == 0.0 )
        {
          Hostonly_ printf("\tError : Number of accumulation is zero\n");
          Exit(0);
        }
        break;
        
      case SPEC_VEL:
        cmp[n].val[var_Velocity] = extractVelLBC(n, vec, tm, v00); // 指定された無次元平均流速
        div_ibc_drchlt_(dv, size, &gd, pitch, st, ed, v00, d_cdf, &n, vec, &fcount);
        break;
        
      case SOLIDREV:
        extractAngularVel(n, vec, ctr, tm, v00);
        div_ibc_sldrev_(dv, size, &gd, st, ed, pitch, v00, d_cdf, &n, vec, ctr, region, &fcount);
        break;
        
      default:
        break;
    }
  }

  
  // 外部境界条件による修正
  for (int face=0; face<NOFACE; face++)
  {
    typ = obc[face].getClass();
    
    switch (typ)
    {
      case OBC_SPEC_VEL:
      case OBC_WALL:
        dummy = extractVelOBC(face, vec, tm, v00);
        vobc_div_drchlt_(dv, size, &gd, pitch, &face, d_cdf, vec, nID);
        break;
        
      case OBC_INTRINSIC:
        if ( C->Mode.Example == id_Jet )
        {
          ((IP_Jet*)Ex)->divJetInflow(dv, face, fcount);
        }
        break;
    }

  }
  
  flop += fcount;
}


// #################################################################
/**
 @brief 圧力損失部による速度の方向修正
 @param[in,out] v 速度
 @param bd BCindex B
 @param cvf コンポーネントの体積率
 @param v00 参照速度
 @param[out] flop
 */
void SetBC3D::mod_Dir_Forcing(REAL_TYPE* d_v, int* d_bd, REAL_TYPE* d_cvf, REAL_TYPE* v00, double &flop)
{
  int st[3], ed[3];
  REAL_TYPE vec[3];
  int gd = guide;
  
  for (int n=1; n<=NoCompo; n++) {
    if ( cmp[n].isFORCING() )
    {
      cmp[n].getBbox(st, ed);
      
      vec[0] = cmp[n].nv[0];
      vec[1] = cmp[n].nv[1];
      vec[2] = cmp[n].nv[2];
      
      switch ( cmp[n].getType() )
      {
        case HEX:
          if ( cmp[n].get_sw_HexDir() ) hex_dir_ (d_v, size, &gd, st, ed, d_bd, d_cvf, &n, v00, vec, &flop);
          break;
          
        case FAN:
          Exit(0);
          break;
          
        case DARCY:
          Exit(0);
          break;
          
        default:
          break;
      }
    }
  }
}


// #################################################################
/**
 @brief 圧力損失部による疑似速度方向の修正
 @param [in,out] vc   疑似速度ベクトル
 @param [in]     v    速度ベクトル n-step
 @param [in]     bd   BCindex B
 @param [in]     cvf  コンポーネントの体積率
 @param [in]     v00  参照速度
 @param [in]     dt   時間積分幅
 @param [in,out] flop 浮動小数点演算数
 */
void SetBC3D::mod_Pvec_Forcing(REAL_TYPE* d_vc, REAL_TYPE* d_v, int* d_bd, REAL_TYPE* d_cvf, REAL_TYPE* v00, REAL_TYPE dt, double &flop)
{
  int st[3], ed[3];
  REAL_TYPE vec[3];
  int gd = guide;
  
  for (int n=1; n<=NoCompo; n++) {
    vec[0] = cmp[n].nv[0];
    vec[1] = cmp[n].nv[1];
    vec[2] = cmp[n].nv[2];
    
    cmp[n].getBbox(st, ed);
    
    switch ( cmp[n].getType() )
    {
      case HEX:
        hex_force_pvec_(d_vc, size, &gd, st, ed, d_bd, d_cvf, d_v, &n, v00, &dt, vec, &cmp[n].ca[0], &flop);
        break;
        
      case FAN:
        Exit(0);
        break;
        
      case DARCY:
        Exit(0);
        break;
        
      default:
        break;
    }

  }
}


// #################################################################
// 圧力損失部によるPoisosn式のソース項の修正とワーク用の速度を保持
void SetBC3D::mod_Psrc_Forcing(REAL_TYPE* s_1, REAL_TYPE* v, int* bd, REAL_TYPE* cvf, REAL_TYPE* v00, REAL_TYPE** c_array, double &flop)
{
  int st[3], ed[3], csz[3];
  REAL_TYPE vec[3];
  int gd = guide;
  REAL_TYPE* w_ptr=NULL;
  
  for (int n=1; n<=NoCompo; n++) {
    vec[0] = cmp[n].nv[0];
    vec[1] = cmp[n].nv[1];
    vec[2] = cmp[n].nv[2];
    
    cmp[n].getBbox(st, ed);
    cmp[n].get_cmp_sz(csz);
    w_ptr = c_array[n];
    
    // w_ptrに
    if ( cmp[n].isFORCING() ) force_keep_vec_(w_ptr, csz, st, ed, v, size, &gd);

    switch ( cmp[n].getType() )
    {
      case HEX:
        hex_psrc_(s_1, size, &gd, st, ed, bd, cvf, w_ptr, csz, &n, v00, vec, &cmp[n].ca[0], &flop);
        break;
        
      case FAN:
        Exit(0);
        break;
        
      case DARCY:
        Exit(0);
        break;
        
      default:
        break;
    }

  }
}


// #################################################################
// 圧力損失部によるセルセンタ速度の修正と速度の発散値の修正
// am[]のインデクスに注意 (Fortran <-> C)
void SetBC3D::mod_Vdiv_Forcing(REAL_TYPE* v, int* bd, REAL_TYPE* cvf, REAL_TYPE* dv, REAL_TYPE dt, REAL_TYPE* v00, Gemini_R* am, REAL_TYPE** c_array, double &flop)
{
  int st[3], ed[3], csz[3];
  REAL_TYPE vec[3];
  REAL_TYPE* w_ptr=NULL;
  int gd = guide;
  REAL_TYPE aa[2];
  
  for (int n=1; n<=NoCompo; n++) {
    if ( cmp[n].isFORCING() )
    {
      cmp[n].getBbox(st, ed);
      cmp[n].get_cmp_sz(csz);
      w_ptr = c_array[n];
      
      vec[0] = cmp[n].nv[0];
      vec[1] = cmp[n].nv[1];
      vec[2] = cmp[n].nv[2];
      
      switch ( cmp[n].getType() )
      {
        case HEX:
          hex_force_vec_(v, dv, size, &gd, st, ed, bd, cvf, w_ptr, csz, &n, v00, &dt, pitch, vec, &cmp[n].ca[0], aa, &flop);
          am[n].p0 = aa[0];
          am[n].p1 = aa[1];
          break;
          
        case FAN:
          Exit(0);
          break;
          
        case DARCY:
          Exit(0);
          break;
          
        default:
          break;
      }
    }
  }
}


// #################################################################
// 速度境界条件による流束の修正
void SetBC3D::modPvecFlux(REAL_TYPE* wv, REAL_TYPE* v, int* d_cdf, const double tm, Control* C, REAL_TYPE* v00, double& flop)
{
  REAL_TYPE vec[3], dummy, ctr[3];
  int st[3], ed[3];
  int typ;
  int gd = guide;
  
  // 内部境界（流束形式）
  for (int n=1; n<=NoCompo; n++)
  {
    
    typ = cmp[n].getType();
    cmp[n].getBbox(st, ed);
    
    if ( typ==SPEC_VEL )
    {
      dummy = extractVelLBC(n, vec, tm, v00);
      
      if ( C->CnvScheme==Control::O1_upwind || C->CnvScheme==Control::O3_muscl )
      {
        pvec_ibc_specv_fvm_(wv, size, &gd, st, ed, pitch, v00, &rei, v, d_cdf, &n, vec, &flop);
      }
      else // Central scheme
      {
        pvec_ibc_specv_fdm_(wv, size, &gd, st, ed, pitch, v00, &rei, v, d_cdf, &n, vec, &flop);
      }
    }
    else if ( typ==OUTFLOW )
    {
      vec[0] = vec[1] = vec[2] = cmp[n].val[var_Velocity]; // modDivergence()でセルフェイス流出速度がval[var_Velocity]にセット
      pvec_ibc_oflow_(wv, size, &gd, st, ed, pitch, &rei, v, d_cdf, &n, vec, &flop);
    }
    else if ( typ==SOLIDREV )
    {
      extractAngularVel(n, vec, ctr, tm, v00);
      
      if ( C->CnvScheme==Control::O1_upwind || C->CnvScheme==Control::O3_muscl )
      {
        pvec_ibc_sldrev_fvm_(wv, size, &gd, st, ed, pitch, &rei, v, d_cdf, &n, vec, ctr, origin, &flop);
      }
      else // Central scheme
      {
        //pvec_ibc_specv_fdm_(wv, size, &gd, st, ed, &dh, v00, &rei, v, d_cdf, &n, vec, &flop);
      }
    }
    
  }
  
  
  // 流束形式の外部境界条件
  for (int face=0; face<NOFACE; face++)
  {
    typ = obc[face].getClass();
    
    switch ( typ )
    {
      case OBC_SPEC_VEL:
        dummy = extractVelOBC(face, vec, tm, v00);
        
        if ( C->CnvScheme==Control::O1_upwind || C->CnvScheme==Control::O3_muscl )
        {
          vobc_pv_specv_fvm_(wv, size, &gd, &face, pitch, &rei, v, d_cdf, vec, nID, &flop);
        }
        else // Central scheme
        {
          vobc_pv_specv_fdm_(wv, size, &gd, &face, pitch, &rei, v, d_cdf, vec, nID, &flop);
        }
        
        break;
        
      case OBC_WALL:
        dummy = extractVelOBC(face, vec, tm, v00);
        vobc_pv_wall_(wv, size, &gd, &face, pitch, &rei, v, vec, nID, &flop);
        break;
        
      case OBC_INTRINSIC:
        if ( C->Mode.Example == id_Jet )
        {
          ((IP_Jet*)Ex)->vobc_pv_JetInflow(wv, face, rei, v, flop);
        }
        break;
    }
  }
  
}


// #################################################################
// 速度境界条件によるPoisosn式のソース項の修正
void SetBC3D::modPsrcVBC(REAL_TYPE* dv, int* d_cdf, const double tm, Control* C, REAL_TYPE* v00, REAL_TYPE* vf, REAL_TYPE* vc, REAL_TYPE* v0, REAL_TYPE dt, double &flop)
{
  int st[3], ed[3];
  REAL_TYPE vec[3], vel, ctr[3];
  int typ;
  int gd = guide;
  double fcount = 0.0;
  
  // 内部境界条件による修正
  for (int n=1; n<=NoCompo; n++) {
    typ = cmp[n].getType();
    cmp[n].getBbox(st, ed);
    
    switch (typ)
    {
      case SPEC_VEL:
      {
        extractVelLBC(n, vec, tm, v00);
        div_ibc_drchlt_(dv, size, &gd, pitch, st, ed, v00, d_cdf, &n, vec, &fcount);
        break;
      }
        
      case OUTFLOW:
        vel = cmp[n].val[var_Velocity]; // modDivergence()でval[var_Velocity]にセット
        div_ibc_oflow_pvec_(dv, size, &gd, st, ed, v00, &vel, &dt, pitch, d_cdf, &n, v0, vf, &fcount);
        break;
        
      case SOLIDREV:
        extractAngularVel(n, vec, ctr, tm, v00);
        div_ibc_sldrev_(dv, size, &gd, st, ed, pitch, v00, d_cdf, &n, vec, ctr, region, &fcount);
        break;
        
      default:
        break;
    }
  }

  
  
  // 外部境界条件による修正
  for (int face=0; face<NOFACE; face++)
  {
    typ = obc[face].getClass();
    
    switch ( typ )
    {
      case OBC_SPEC_VEL:
      case OBC_WALL:
      {
        extractVelOBC(face, vec, tm, v00);
        vobc_div_drchlt_(dv, size, &gd, pitch, &face, d_cdf, vec, nID);
        break;
      }
        
      case OBC_INTRINSIC:
        if ( C->Mode.Example == id_Jet )
        {
          ((IP_Jet*)Ex)->divJetInflow(dv, face, flop);
        }
        break;
        
        // 他は境界値を与え，通常スキームで計算するので不要
    }

  }
  
  flop += fcount;
}


// #################################################################
// 圧力の外部境界条件
void SetBC3D::OuterPBC(REAL_TYPE* d_p, const int* ens)
{
  int gd = guide;
  
  for (int face=0; face<NOFACE; face++)
  {
    if ( obc[face].getClass() == OBC_TRC_FREE )
    {
      REAL_TYPE pv = FBUtility::convPrsD2ND(obc[face].p, BasePrs, rho_0, RefV, Unit_Prs);
      pobc_drchlt_ (d_p, size, &gd, &face, &pv, nID);
    }
  }
  
  // 周期境界
  PobcPeriodicSimple(d_p, ens);
  PobcPeriodicDirectional(d_p, ens);
}


// #################################################################
// 拡散項計算時の温度の外部部境界処理
void SetBC3D::OuterTBCdiffusion(REAL_TYPE* d_qbc, REAL_TYPE* d_ie, const REAL_TYPE* d_ie0, const int* d_bcd, Control* C)
{
  int H;
  
  for (int face=0; face<NOFACE; face++)
  {
    // 各面の無次元流入出熱量
    REAL_TYPE va = 0.0;
    
    H = obc[face].getHTmode();
    
    if ( obc[face].getClass() == OBC_WALL )
    {
      switch ( obc[face].getHtype() ) // 熱境界条件の種類
      {
        case HEATFLUX:
          va = psObcHeatflux(d_qbc, face);
          break;
          
        case TRANSFER:
          if      ( H == HT_S)  va = psObcHeatTransferS (d_qbc, face, d_bcd, d_ie, d_ie0);
          else if ( H == HT_SF) va = psObcHeatTransferSF(d_qbc, face, d_bcd, d_ie, d_ie0);
          else if ( H == HT_SN) va = psObcHeatTransferSN(d_qbc, face, d_bcd, d_ie, d_ie0);
          break;
          
        case ISOTHERMAL:
          va = psObcIsoThermal(d_qbc, face, d_bcd, d_ie, d_ie0);
          break;
      }

    }
    
    if ( numProc > 1 )
    {
      REAL_TYPE tmp = va;
      if ( paraMngr->Allreduce(&tmp, &va, 1, MPI_SUM) != CPM_SUCCESS )
      {
        Hostonly_ printf("Allreduce Error\n");
        Exit(0);
      }
    }
    
    C->H_Dface[face] += va;
  }
}


// #################################################################
// 温度の外部周期境界条件
void SetBC3D::OuterTBCperiodic(REAL_TYPE* d_ie, const int* ens)
{
  TobcPeriodicSimple(d_ie, ens);
}



// #################################################################
// 速度の外部境界条件処理（Div反復内で値を指定する境界条件）
void SetBC3D::OuterVBC(REAL_TYPE* d_v, REAL_TYPE* d_vf, int* d_cdf, const double tm, Control* C, REAL_TYPE* v00, const int* ens)
{
  REAL_TYPE vec[3];
  int gd = guide;
  REAL_TYPE dd=0.0;

  
  // 周期境界以外
  for (int face=0; face<NOFACE; face++)
  {
    switch ( obc[face].getClass() )
    {
      case OBC_TRC_FREE:
        vobc_cf_tfree_(d_vf, size, &gd, &face, nID);
        //if ( numProc > 1 )
        //{
        //  if ( paraMngr->BndCommV3D(d_vf, size[0], size[1], size[2], gd, gd) != CPM_SUCCESS ) Exit(0);
        //}
        
        vobc_cc_tfree_(d_v, size, &gd, &face, d_vf, nID);
        if ( numProc > 1 )
        {
          if ( paraMngr->BndCommV3D(d_v, size[0], size[1], size[2], gd, gd) != CPM_SUCCESS ) Exit(0);
        }
        break;
        
      case OBC_FAR_FIELD:
        vobc_cc_neumann_(d_v, size, &gd, &face, nID);
        break;
        
      case OBC_OUTFLOW:
        // {}^c \hat{u} で計算した値　vobc_cc_neumann_(d_v, size, &gd, &face, nID);
        break;
        
      case OBC_WALL:
      case OBC_SPEC_VEL:
        extractVelOBC(face, vec, tm, v00);
        vobc_cc_drchlt_(d_v, size, &gd, &face, vec, nID);
        vobc_face_drchlt_(d_vf, size, &gd, &face, d_cdf, vec, nID);
        break;
        
      case OBC_INTRINSIC:
        if ( C->Mode.Example == id_Jet )
        {
          ((IP_Jet*)Ex)->vobcJetInflowGC(d_v, face);
        }
        break;
        
      default:
        break;
    }
  }
  
  // 周期境界
  VobcPeriodicSimple(d_v, ens); // セルフェイスの値の周期処理は不要
  
}



// #################################################################
// 速度の外部境界条件処理（セルフェイス境界のための準備）
void SetBC3D::OuterVBCfacePrep(REAL_TYPE* d_vc, REAL_TYPE* d_v, int* d_cdf, REAL_TYPE dt, Control* C, const int* ens, const unsigned m_step)
{
  int gd = guide;
  REAL_TYPE v_cnv=0.0;
  
  // 周期境界以外
  
  // 最初ｎ1ステップ目は，ガイドセル部の情報がないのでノイマンを使う
  if ( m_step == 1 )
  {
    for (int face=0; face<NOFACE; face++)
    {
      switch ( obc[face].getClass() )
      {
        case OBC_TRC_FREE:
        case OBC_FAR_FIELD:
        case OBC_OUTFLOW:
          vobc_cc_neumann_(d_vc, size, &gd, &face, nID);
          break;
          
        case OBC_INTRINSIC:
          if ( C->Mode.Example == id_Jet )
          {
            ((IP_Jet*)Ex)->vobcJetInflowGC(d_vc, face);
          }
          break;
          
        default:
          break;
      }
    }
  }
  else
  {
    for (int face=0; face<NOFACE; face++)
    {
      switch ( obc[face].getClass() )
      {
        case OBC_TRC_FREE:
          vobc_cc_neumann_(d_vc, size, &gd, &face, nID);
          break;
          
        case OBC_FAR_FIELD:
          vobc_cc_neumann_(d_vc, size, &gd, &face, nID);
          break;
          
        case OBC_OUTFLOW:
          v_cnv = C->V_Dface[face];
          vobc_cc_outflow_(d_vc, d_v, size, &gd, pitch, &dt, d_cdf, &v_cnv, &face, nID);
          break;
          
        case OBC_INTRINSIC:
          if ( C->Mode.Example == id_Jet )
          {
            ((IP_Jet*)Ex)->vobcJetInflowGC(d_vc, face);
          }
          break;
          
        default:
          break;
      }
    }
  }

  
  // 周期境界
  VobcPeriodicSimple(d_vc, ens);
  
}



// #################################################################
/**
 * @brief 圧力の外部周期境界条件（単純なコピー）
 * @param [in,out] d_p  圧力
 * @param [in]     ens  周期境界方向フラグ
 * @note 並列時には全ノードがPeriodicCommS3D()を同じ引数で評価すること
 */
void SetBC3D::PobcPeriodicSimple(REAL_TYPE* d_p, const int* ens)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
  if ( numProc > 1 ) 
  {
    // X方向
    if ( (ens[0] == ON) && (obc[X_minus].getPrdcMode() == BoundaryOuter::prdc_Simple) )
    {
      if ( paraMngr->PeriodicCommS3D(d_p, ix, jx, kx, gd, 1, X_DIR, PLUS2MINUS) != CPM_SUCCESS ) Exit(0);
      if ( paraMngr->PeriodicCommS3D(d_p, ix, jx, kx, gd, 1, X_DIR, MINUS2PLUS) != CPM_SUCCESS ) Exit(0);
    }
    
    // Y方向
    if ( (ens[1] == ON) && (obc[Y_minus].getPrdcMode() == BoundaryOuter::prdc_Simple) )
    {
      if ( paraMngr->PeriodicCommS3D(d_p, ix, jx, kx, gd, 1, Y_DIR, PLUS2MINUS) != CPM_SUCCESS ) Exit(0);
      if ( paraMngr->PeriodicCommS3D(d_p, ix, jx, kx, gd, 1, Y_DIR, MINUS2PLUS) != CPM_SUCCESS ) Exit(0);
    }
    
    // Z方向
    if ( (ens[2] == ON) && (obc[Z_minus].getPrdcMode() == BoundaryOuter::prdc_Simple) )
    {
      if ( paraMngr->PeriodicCommS3D(d_p, ix, jx, kx, gd, 1, Z_DIR, PLUS2MINUS) != CPM_SUCCESS ) Exit(0);
      if ( paraMngr->PeriodicCommS3D(d_p, ix, jx, kx, gd, 1, Z_DIR, MINUS2PLUS) != CPM_SUCCESS ) Exit(0);
    }
    
  }
  else // Serial
  {

    // X方向
    if ( (ens[0] == ON) && (obc[X_minus].getPrdcMode() == BoundaryOuter::prdc_Simple) )
    {
#pragma omp parallel for firstprivate(ix, jx, kx, gd) schedule(static)
      for (int k=1; k<=kx; k++) {
        for (int j=1; j<=jx; j++) {
          size_t m0 = _F_IDX_S3D(0,  j, k, ix, jx, kx, gd);
          size_t m1 = _F_IDX_S3D(ix, j, k, ix, jx, kx, gd);
          d_p[m0] = d_p[m1];
        }
      }
      
#pragma omp parallel for firstprivate(ix, jx, kx, gd) schedule(static)
      for (int k=1; k<=kx; k++) {
        for (int j=1; j<=jx; j++) {
          size_t m0 = _F_IDX_S3D(ix+1, j, k, ix, jx, kx, gd);
          size_t m1 = _F_IDX_S3D(1,    j, k, ix, jx, kx, gd);
          d_p[m0] = d_p[m1];
        }
      }
    }
    
    // Y方向
    if ( (ens[1] == ON) && (obc[Y_minus].getPrdcMode() == BoundaryOuter::prdc_Simple) )
    {
#pragma omp parallel for firstprivate(ix, jx, kx, gd) schedule(static)
      for (int k=1; k<=kx; k++) {
        for (int i=1; i<=ix; i++) {
          size_t m0 = _F_IDX_S3D(i, 0,  k, ix, jx, kx, gd);
          size_t m1 = _F_IDX_S3D(i, jx, k, ix, jx, kx, gd);
          d_p[m0] = d_p[m1];
        }
      }
      
#pragma omp parallel for firstprivate(ix, jx, kx, gd) schedule(static)
      for (int k=1; k<=kx; k++) {
        for (int i=1; i<=ix; i++) {
          size_t m0 = _F_IDX_S3D(i, jx+1, k, ix, jx, kx, gd);
          size_t m1 = _F_IDX_S3D(i, 1,    k, ix, jx, kx, gd);
          d_p[m0] = d_p[m1];
        }
      }
    }
    
    // Z方向
    if ( (ens[2] == ON) && (obc[Z_minus].getPrdcMode() == BoundaryOuter::prdc_Simple) )
    {
#pragma omp parallel for firstprivate(ix, jx, kx, gd) schedule(static)
      for (int j=1; j<=jx; j++) {
        for (int i=1; i<=ix; i++) {
          size_t m0 = _F_IDX_S3D(i, j, 0,  ix, jx, kx, gd);
          size_t m1 = _F_IDX_S3D(i, j, kx, ix, jx, kx, gd);
          d_p[m0] = d_p[m1];
        }
      }
      
#pragma omp parallel for firstprivate(ix, jx, kx, gd) schedule(static)
      for (int j=1; j<=jx; j++) {
        for (int i=1; i<=ix; i++) {
          size_t m0 = _F_IDX_S3D(i, j, kx+1, ix, jx, kx, gd);
          size_t m1 = _F_IDX_S3D(i, j, 1,    ix, jx, kx, gd);
          d_p[m0] = d_p[m1];
        }
      }
    }
    
  }
}


// #################################################################
/**
 * @brief 圧力の外部周期境界条件（双方向に圧力差を設定）
 * @param [in,out] d_p  圧力のデータクラス
 * @param [in]     ens  周期境界方向フラグ
 */
void SetBC3D::PobcPeriodicDirectional(REAL_TYPE* d_p, const int* ens)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
  REAL_TYPE pv; // 圧力差
  REAL_TYPE pd;
  
  if ( numProc > 1 ) 
  {
    // X方向
    if ( (ens[0] == ON) && (obc[X_minus].getPrdcMode() == BoundaryOuter::prdc_Directional) )
    {
      if ( paraMngr->PeriodicCommS3D(d_p, ix, jx, kx, gd, 1, X_DIR, PLUS2MINUS) != CPM_SUCCESS ) Exit(0);
      if ( nID[X_minus] < 0 )
      {
        // 上流面か下流面かで，圧力差の方向を逆転する
        pv = FBUtility::convPrsD2ND(obc[X_minus].p, BasePrs, rho_0, RefV, Unit_Prs);
        pd = ( obc[X_minus].get_FaceMode() == BoundaryOuter::prdc_upstream ) ? pv : -pv;
        
#pragma omp parallel for firstprivate(ix, jx, kx, gd, pd) schedule(static)
        for (int k=1; k<=kx; k++) {
          for (int j=1; j<=jx; j++) {
            size_t m0 = _F_IDX_S3D(0, j, k, ix, jx, kx, gd);
            d_p[m0] += pd;
          }
        }
      }

      
      if ( paraMngr->PeriodicCommS3D(d_p, ix, jx, kx, gd, 1, X_DIR, MINUS2PLUS) != CPM_SUCCESS ) Exit(0);
      if ( nID[X_plus] < 0 )
      {
        pv = FBUtility::convPrsD2ND(obc[X_plus].p, BasePrs, rho_0, RefV, Unit_Prs);
        pd = ( obc[X_plus].get_FaceMode() == BoundaryOuter::prdc_upstream ) ? pv : -pv;
        
#pragma omp parallel for firstprivate(ix, jx, kx, gd, pd) schedule(static)
        for (int k=1; k<=kx; k++) {
          for (int j=1; j<=jx; j++) {
            size_t m0 = _F_IDX_S3D(ix+1, j, k, ix, jx, kx, gd);
            d_p[m0] += pd;
          }
        }
      }
    }
    
    // Y方向
    if ( (ens[1] == ON) && (obc[Y_minus].getPrdcMode() == BoundaryOuter::prdc_Directional) )
    {
      if ( paraMngr->PeriodicCommS3D(d_p, ix, jx, kx, gd, 1, Y_DIR, PLUS2MINUS) != CPM_SUCCESS ) Exit(0);
      if ( nID[Y_minus] >= 0 )
      {
        pv = FBUtility::convPrsD2ND(obc[Y_minus].p, BasePrs, rho_0, RefV, Unit_Prs);
        pd = ( obc[Y_minus].get_FaceMode() == BoundaryOuter::prdc_upstream ) ? pv : -pv;
        
#pragma omp parallel for firstprivate(ix, jx, kx, gd, pd) schedule(static)
        for (int k=1; k<=kx; k++) {
          for (int i=1; i<=ix; i++) {
            size_t m0 = _F_IDX_S3D(i, 0, k, ix, jx, kx, gd);
            d_p[m0] += pd;
          }
        }
      }
      
      if ( paraMngr->PeriodicCommS3D(d_p, ix, jx, kx, gd, 1, Y_DIR, MINUS2PLUS) != CPM_SUCCESS ) Exit(0);
      if ( nID[Y_plus] < 0 )
      {
        pv = FBUtility::convPrsD2ND(obc[Y_plus].p, BasePrs, rho_0, RefV, Unit_Prs);
        pd = ( obc[Y_plus].get_FaceMode() == BoundaryOuter::prdc_upstream ) ? pv : -pv;
        
#pragma omp parallel for firstprivate(ix, jx, kx, gd, pd) schedule(static)
        for (int k=1; k<=kx; k++) {
          for (int i=1; i<=ix; i++) {
            size_t m0 = _F_IDX_S3D(i, jx+1, k, ix, jx, kx, gd);
            d_p[m0] += pd;
          }
        }
      }
      
    }
    
    // Z方向
    if ( (ens[2] == ON) && (obc[Z_minus].getPrdcMode() == BoundaryOuter::prdc_Directional) )
    {
      if ( paraMngr->PeriodicCommS3D(d_p, ix, jx, kx, gd, 1, Z_DIR, PLUS2MINUS) != CPM_SUCCESS ) Exit(0);
      if ( nID[Z_minus] < 0 )
      {
        pv = FBUtility::convPrsD2ND(obc[Z_minus].p, BasePrs, rho_0, RefV, Unit_Prs);
        pd = ( obc[Z_minus].get_FaceMode() == BoundaryOuter::prdc_upstream ) ? pv : -pv;
        
#pragma omp parallel for firstprivate(ix, jx, kx, gd, pd) schedule(static)
        for (int j=1; j<=jx; j++) {
          for (int i=1; i<=ix; i++) {
            size_t m0 = _F_IDX_S3D(i, j, 0, ix, jx, kx, gd);
            d_p[m0] += pd;
          }
        }
      }
      
      if ( paraMngr->PeriodicCommS3D(d_p, ix, jx, kx, gd, 1, Z_DIR, MINUS2PLUS) != CPM_SUCCESS ) Exit(0);
      if ( nID[Z_plus] < 0 )
      {
        pv = FBUtility::convPrsD2ND(obc[Z_plus].p, BasePrs, rho_0, RefV, Unit_Prs);
        pd = ( obc[Z_plus].get_FaceMode() == BoundaryOuter::prdc_upstream ) ? pv : -pv;
        
#pragma omp parallel for firstprivate(ix, jx, kx, gd, pd) schedule(static)
        for (int j=1; j<=jx; j++) {
          for (int i=1; i<=ix; i++) {
            size_t m0 = _F_IDX_S3D(i, j, kx+1, ix, jx, kx, gd);
            d_p[m0] += pd;
          }
        }
      }
      
    }
  }
  else // Serial
  {
    // X方向
    if ( (ens[0] == ON) && (obc[X_minus].getPrdcMode() == BoundaryOuter::prdc_Directional) )
    {
      pv = FBUtility::convPrsD2ND(obc[X_minus].p, BasePrs, rho_0, RefV, Unit_Prs);
      pd = ( obc[X_minus].get_FaceMode() == BoundaryOuter::prdc_upstream ) ? pv : -pv;
      
#pragma omp parallel for firstprivate(ix, jx, kx, gd, pd) schedule(static)
      for (int k=1; k<=kx; k++) {
        for (int j=1; j<=jx; j++) {
          size_t m0 = _F_IDX_S3D(0,  j, k, ix, jx, kx, gd);
          size_t m1 = _F_IDX_S3D(ix, j, k, ix, jx, kx, gd);
          d_p[m0] = d_p[m1] + pd;
        }
      }
      
      pv = FBUtility::convPrsD2ND(obc[X_plus].p, BasePrs, rho_0, RefV, Unit_Prs);
      pd = ( obc[X_plus].get_FaceMode() == BoundaryOuter::prdc_upstream ) ? pv : -pv;
      
#pragma omp parallel for firstprivate(ix, jx, kx, gd, pd) schedule(static)
      for (int k=1; k<=kx; k++) {
        for (int j=1; j<=jx; j++) {
          size_t m0 = _F_IDX_S3D(ix+1, j, k, ix, jx, kx, gd);
          size_t m1 = _F_IDX_S3D(1,    j, k, ix, jx, kx, gd);
          d_p[m0] = d_p[m1] + pd;
        }
      }
    }
    
    // Y方向
    if ( (ens[1] == ON) && (obc[Y_minus].getPrdcMode() == BoundaryOuter::prdc_Directional) )
    {
      pv = FBUtility::convPrsD2ND(obc[Y_minus].p, BasePrs, rho_0, RefV, Unit_Prs);
      pd = ( obc[Y_minus].get_FaceMode() == BoundaryOuter::prdc_upstream ) ? pv : -pv;
      
#pragma omp parallel for firstprivate(ix, jx, kx, gd, pd) schedule(static)
      for (int k=1; k<=kx; k++) {
        for (int i=1; i<=ix; i++) {
          size_t m0 = _F_IDX_S3D(i, 0,  k, ix, jx, kx, gd);
          size_t m1 = _F_IDX_S3D(i, jx, k, ix, jx, kx, gd);
          d_p[m0] = d_p[m1] + pd;
        }
      }
      
      pv = FBUtility::convPrsD2ND(obc[Y_plus].p, BasePrs, rho_0, RefV, Unit_Prs);
      pd = ( obc[Y_plus].get_FaceMode() == BoundaryOuter::prdc_upstream ) ? pv : -pv;
      
#pragma omp parallel for firstprivate(ix, jx, kx, gd, pd) schedule(static)
      for (int k=1; k<=kx; k++) {
        for (int i=1; i<=ix; i++) {
          size_t m0 = _F_IDX_S3D(i, jx+1, k, ix, jx, kx, gd);
          size_t m1 = _F_IDX_S3D(i, 1,    k, ix, jx, kx, gd);
          d_p[m0] = d_p[m1] + pd;
        }
      }
    }
    
    // Z方向
    if ( (ens[2] == ON) && (obc[Z_minus].getPrdcMode() == BoundaryOuter::prdc_Directional) )
    {
      pv = FBUtility::convPrsD2ND(obc[Z_minus].p, BasePrs, rho_0, RefV, Unit_Prs);
      pd = ( obc[Z_minus].get_FaceMode() == BoundaryOuter::prdc_upstream ) ? pv : -pv;
      
#pragma omp parallel for firstprivate(ix, jx, kx, gd, pd) schedule(static)
      for (int j=1; j<=jx; j++) {
        for (int i=1; i<=ix; i++) {
          size_t m0 = _F_IDX_S3D(i, j, 0,  ix, jx, kx, gd);
          size_t m1 = _F_IDX_S3D(i, j, kx, ix, jx, kx, gd);
          d_p[m0] = d_p[m1] + pd;
        }
      }
    
      pv = FBUtility::convPrsD2ND(obc[Z_plus].p, BasePrs, rho_0, RefV, Unit_Prs);
      pd = ( obc[Z_plus].get_FaceMode() == BoundaryOuter::prdc_upstream ) ? pv : -pv;
      
#pragma omp parallel for firstprivate(ix, jx, kx, gd, pd) schedule(static)
      for (int j=1; j<=jx; j++) {
        for (int i=1; i<=ix; i++) {
          size_t m0 = _F_IDX_S3D(i, j, kx+1, ix, jx, kx, gd);
          size_t m1 = _F_IDX_S3D(i, j, 1,    ix, jx, kx, gd);
          d_p[m0] = d_p[m1] + pd;
        }
      }
    }
    
  }
}


// #################################################################
/**
 @brief 圧力の内部周期境界条件（一方向の圧力差）
 @param d_p 圧力のデータクラス
 @param st コンポーネント範囲の開始インデクス
 @param ed コンポーネント範囲の終了インデクス
 @param d_bcd BCindex B
 @param odr 下流面のidが格納されているコンポーネントエントリ
 @param dir ドライバの設置方向
 @param pd 圧力差
 */
void SetBC3D::Pibc_Prdc(REAL_TYPE* d_p, int* st, int* ed, int* d_bx, int odr, int dir, REAL_TYPE pd)
{

  if ( numProc > 1 ) 
  {
    printf("Error : 'PibcPrdc' method is limited to use for serial execution\n.");
    Exit(0);
  }
  
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
  size_t m0, m1;
  
  switch (dir) 
  {
    case X_minus:
      if ( nID[dir] < 0 )
      {
        int i = st[0];
        for (int k=st[2]; k<=ed[2]; k++) {
          for (int j=st[1]; j<=ed[1]; j++) {
            m1 = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
            if ( DECODE_CMP(d_bx[m1]) == odr )
            {
              m0 = _F_IDX_S3D(0, j, k, ix, jx, kx, gd);
              d_p[m0] = d_p[m1] + pd;
            }
          }
        }
      }
      break;
      
    case X_plus:
      if ( nID[dir] < 0 )
      {
        int i = st[0];
        for (int k=st[2]; k<=ed[2]; k++) {
          for (int j=st[1]; j<=ed[1]; j++) {
            m1 = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
            if ( DECODE_CMP(d_bx[m1]) == odr )
            {
              m0 = _F_IDX_S3D(ix+1, j, k, ix, jx, kx, gd);
              d_p[m0] = d_p[m1] + pd;              
            }
          }
        }
      }
      break;
      
    case Y_minus:
      if ( nID[dir] < 0 )
      {
        int j = st[1];
        for (int k=st[2]; k<=ed[2]; k++) {
          for (int i=st[0]; i<=ed[0]; i++) {
            m1 = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
            if ( DECODE_CMP(d_bx[m1]) == odr )
            {
              m0 = _F_IDX_S3D(i, 0, k, ix, jx, kx, gd);
              d_p[m0] = d_p[m1] + pd;
            }
          }
        }
      }
      break;
      
    case Y_plus:
      if ( nID[dir] < 0 )
      {
        int j = st[1];
        for (int k=st[2]; k<=ed[2]; k++) {
          for (int i=st[0]; i<=ed[0]; i++) {
            m1 = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
            if ( DECODE_CMP(d_bx[m1]) == odr )
            {
              m0 = _F_IDX_S3D(i, jx+1, k, ix, jx, kx, gd);
              d_p[m0] = d_p[m1] + pd;
            }
          }
        }
      }
      break;
      
    case Z_minus:
      if ( nID[dir] < 0 )
      {
        int k = st[2];
        for (int j=st[1]; j<=ed[1]; j++) {
          for (int i=st[0]; i<=ed[0]; i++) {
            m1 = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
            if ( DECODE_CMP(d_bx[m1]) == odr )
            {
              m0 = _F_IDX_S3D(i, j, 0, ix, jx, kx, gd);
              d_p[m0] = d_p[m1] + pd;
            }
          }
        }
      }
      break;
      
    case Z_plus:
      if ( nID[dir] < 0 )
      {
        int k = st[2];
        for (int j=st[1]; j<=ed[1]; j++) {
          for (int i=st[0]; i<=ed[0]; i++) {
            m1 = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
            if ( DECODE_CMP(d_bx[m1]) == odr )
            {
              m0 = _F_IDX_S3D(i, j, kx+1, ix, jx, kx, gd);
              d_p[m0] = d_p[m1] + pd;
            }
          }
        }
      }
      break;
  }
}


// #################################################################
/**
 @brief 内部領域の速度のディリクレ境界をセットする
 @param v 速度ベクトル
 @param bx BCindex
 @param n 境界条件コンポーネントのエントリ番号
 @param v00 参照速度
 @retval 無次元平均流速
 @note
    - 内部境界はセル面直を仮定
 
REAL_TYPE SetBC3D::setDirectForcing(REAL_TYPE* v, int* bx, int n, REAL_TYPE v00)
{
  int i,j,k;
  int s;
  REAL_TYPE vel, vec[3];
  
  vel   = FBUtility::convVelD2ND(cmp[n].D1.Velocity, RefV)  * v00;
  vec[0]= cmp[n].nv[0]*vel;
  vec[1]= cmp[n].nv[1]*vel;
  vec[2]= cmp[n].nv[2]*vel;
  
  for (k=cmp[n].ci.st[2]-1; k<=cmp[n].ci.ed[2]; k++) {
    for (j=cmp[n].ci.st[1]-1; j<=cmp[n].ci.ed[1]; j++) {
      for (i=cmp[n].ci.st[0]-1; i<=cmp[n].ci.ed[0]; i++) {
        s = bx[ FBUtility::getFindexS3D(sz, gd, i, j, k) ];
        if ( (s & BIT_MASK_10) == n ) {
					if ( TEST_BIT(s, VFACE_I) ) {
						v[FBUtility::getFindexV3D(sz, gd, i, j, k, 0)] = vec[0];
					}
					if ( TEST_BIT(s, VFACE_J) ) {
						v[FBUtility::getFindexV3D(sz, gd, i, j, k, 1)] = vec[1];
					}
					if ( TEST_BIT(s, VFACE_K) ) {
						v[FBUtility::getFindexV3D(sz, gd, i, j, k, 2)] = vec[2];
					}
        }
      }
    }
  }

  return (vel);
}*/



// #################################################################
/**
 * @brief 温度一定の境界条件
 * @param [in,out] d_ie  内部エネルギー
 * @param [in]     d_bcd BCindex B
 * @param [in]     n     境界条件コンポーネントのエントリ番号
 */
void SetBC3D::psIbcConstTemp(REAL_TYPE* d_ie, const int* d_bcd, const int n)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  int odr = n;
  int st[3], ed[3];
  
  REAL_TYPE rho = mat[n].P[p_density] / rho_0;
  REAL_TYPE cp  = mat[n].P[p_specific_heat] / cp_0;
  REAL_TYPE tmp = FBUtility::convTempD2ND(cmp[n].getTemp(), BaseTemp, DiffTemp);
  REAL_TYPE ref = rho * cp * tmp;
  
  cmp[n].getBbox(st, ed);
  
#pragma omp parallel for firstprivate(ix, jx, kx, gd, st, ed, odr, ref) schedule(static)
  for (int k=st[2]; k<=ed[2]; k++) {
    for (int j=st[1]; j<=ed[1]; j++) {
      for (int i=st[0]; i<=ed[0]; i++) {
        size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        if ( (d_bcd[m] & MASK_5) == odr ) d_ie[m] = ref;
      }
    }
  }
}



// #################################################################
/**
 * @brief 内部領域の熱流束指定境界条件
 * @retval 無次元の単位時間あたりの熱量（有次元では[W]）
 * @param [in,out] d_qbc 境界条件熱流束
 * @param [in]     d_cdf BCindex C
 * @param [in]     n     コンポーネントリストのインデクス
 * @note
   - モニタ量の熱量vaは系に対する流入量なので，基準温度に対する熱量
   - 流入量を正にとる
 */
REAL_TYPE SetBC3D::psIbcHeatflux(REAL_TYPE* d_qbc, const int* d_cdf, const int n)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  int odr= n;
  int st[3], ed[3];
  REAL_TYPE va=0.0;
  REAL_TYPE q = cmp[n].getHeatflux() / (RefV*DiffTemp*rho_0*cp_0); // [W/m^2]を無次元化
  
  REAL_TYPE Sx = S_x;
  REAL_TYPE Sy = S_y;
  REAL_TYPE Sz = S_z;
  
  cmp[n].getBbox(st, ed);
	
#pragma omp parallel for firstprivate(ix, jx, kx, gd, st, ed, odr, q, Sx, Sy, Sz) schedule(static) reduction(+:va)
  for (int k=st[2]; k<=ed[2]; k++) {
    for (int j=st[1]; j<=ed[1]; j++) {
      for (int i=st[0]; i<=ed[0]; i++) {
        
        size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        int s = d_cdf[m];
        
        if ( GET_FACE_BC(s, BC_FACE_W) == odr )
        {
          d_qbc[_F_IDX_S4DEX(X_minus, i, j, k, NOFACE, ix, jx, kx, gd)] += q;
          va += q * Sx; // マイナス面の正の値は流入
        }
        
        if ( GET_FACE_BC(s, BC_FACE_E) == odr )
        {
          d_qbc[_F_IDX_S4DEX(X_plus, i, j, k, NOFACE, ix, jx, kx, gd)] -= q;
          va -= q * Sx; // プラス面の正の値は流出
        }
        
        if ( GET_FACE_BC(s, BC_FACE_S) == odr )
        {
          d_qbc[_F_IDX_S4DEX(Y_minus, i, j, k, NOFACE, ix, jx, kx, gd)] += q;
          va += q * Sy;
        }
        
        if ( GET_FACE_BC(s, BC_FACE_N) == odr )
        {
          d_qbc[_F_IDX_S4DEX(Y_plus, i, j, k, NOFACE, ix, jx, kx, gd)] -= q;
          va -= q * Sy;
        }
        
        if ( GET_FACE_BC(s, BC_FACE_B) == odr )
        {
          d_qbc[_F_IDX_S4DEX(Z_minus, i, j, k, NOFACE, ix, jx, kx, gd)] += q;
          va += q * Sz;
        }
        
        if ( GET_FACE_BC(s, BC_FACE_T) == odr )
        {
          d_qbc[_F_IDX_S4DEX(Z_plus, i, j, k, NOFACE, ix, jx, kx, gd)] -= q;
          va -= q * Sz;
        }
      }
    }
  }
  
  if ( numProc > 1 )
  {
		REAL_TYPE tmp = va;
    if ( paraMngr->Allreduce(&tmp, &va, 1, MPI_SUM) != CPM_SUCCESS )
    {
      Hostonly_ printf("Allreduce Error\n");
      Exit(0);
    }
	}
  
  return va;
}


// #################################################################
/**
 * @brief 発熱境界条件
 * @retval 無次元の短時間あたりの発熱量（有次元では[W]）
 * @param [in,out] d_ie  内部エネルギー
 * @param [in]     d_bcd BCindex B
 * @param [in]     n     コンポーネントのエントリ番号
 * @param [in]     dt    時間積分幅
 * @note 発熱密度はControl::setParameters()で計算 D2に発熱密度が保存されている
 */
REAL_TYPE SetBC3D::psIbcHeatGen(REAL_TYPE* d_ie, const int* d_bcd, const int n, const REAL_TYPE dt)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  int odr= n;
  int st[3], ed[3];

  REAL_TYPE c = 0.0;
  REAL_TYPE hs = FBUtility::convHsrcD2ND(cmp[n].getHeatDensity(), RefV, RefL, DiffTemp, rho_0, cp_0);
  REAL_TYPE dhs = dt * hs;

  cmp[n].getBbox(st, ed);
	
#pragma omp parallel for firstprivate(ix, jx, kx, gd, st, ed, odr, dhs, hs) schedule(static) reduction(+:c)
  for (int k=st[2]; k<=ed[2]; k++) {
    for (int j=st[1]; j<=ed[1]; j++) {
      for (int i=st[0]; i<=ed[0]; i++) {
        
        size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        if ( (d_bcd[m] & MASK_5) == odr )
        {
          d_ie[m] += dhs;
          c += hs;
        }
      }
    }
  }
  
  return c*pitch[0]*pitch[1]*pitch[2];
}


// #################################################################
/**
 * @brief 等温境界条件
 * @retval 無次元の単位時間あたりの熱量（有次元では[W]）
 * @param [in,out] d_qbc  境界条件熱流束
 * @param [in]     d_cdf  BCindex C
 * @param [in]     d_bcd  BCindex B
 * @param [in]     n      境界条件コンポーネントのエントリ番号
 * @param [in]     d_ie0  n時刻の内部エネルギー
 * @note d_qbcへ熱流束を加算（他の条件と合成）
 */
REAL_TYPE SetBC3D::psIbcIsoThermal(REAL_TYPE* d_qbc, const int* d_cdf, const int* d_bcd, const int n, const REAL_TYPE* d_ie0)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  int odr= n;
  int st[3], ed[3];
  
  REAL_TYPE va = 0.0;
  REAL_TYPE sf = FBUtility::convTempD2ND(cmp[n].getTemp(), BaseTemp, DiffTemp); // 保持されている温度はCelsius
  REAL_TYPE ppx = 2.0 / pitch[0];
  REAL_TYPE ppy = 2.0 / pitch[1];
  REAL_TYPE ppz = 2.0 / pitch[2];
  
  REAL_TYPE Sx = S_x;
  REAL_TYPE Sy = S_y;
  REAL_TYPE Sz = S_z;

  cmp[n].getBbox(st, ed);
  
#pragma omp parallel for firstprivate(ix, jx, kx, gd, st, ed, odr, ppx, ppy, ppz, Sx, Sy, Sz, sf) \
schedule(static) reduction(+:va)
  for (int k=st[2]; k<=ed[2]; k++) {
    for (int j=st[1]; j<=ed[1]; j++) {
      for (int i=st[0]; i<=ed[0]; i++) {
        
        size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        int s = d_cdf[m];
        REAL_TYPE q = 0.0;
        int l = d_bcd[m] & MASK_5;
        REAL_TYPE rho = mtbl[3*l+0];
        REAL_TYPE cp  = mtbl[3*l+1];
        REAL_TYPE lmd = mtbl[3*l+2];
        REAL_TYPE t = d_ie0[m] / (rho * cp);
        
        if ( GET_FACE_BC(s, BC_FACE_W) == odr )
        {
          q = lmd * ppx * (sf - t);
          d_qbc[_F_IDX_S4DEX(X_minus, i, j, k, NOFACE, ix, jx, kx, gd)] += q;
          va += q * Sx; // マイナス面の正の値は流入
        }
        
        if ( GET_FACE_BC(s, BC_FACE_E) == odr )
        {
          q = lmd * ppx * (t - sf);
          d_qbc[_F_IDX_S4DEX(X_plus, i, j, k, NOFACE, ix, jx, kx, gd)] += q;
          va -= q * Sx; // プラス面の正の値は流出
        }
        
        if ( GET_FACE_BC(s, BC_FACE_S) == odr )
        {
          q = lmd * ppy * (sf - t);
          d_qbc[_F_IDX_S4DEX(Y_minus, i, j, k, NOFACE, ix, jx, kx, gd)] += q;
          va += q * Sy;
        }
        
        if ( GET_FACE_BC(s, BC_FACE_N) == odr )
        {
          q = lmd * ppy * (t - sf);
          d_qbc[_F_IDX_S4DEX(Y_plus, i, j, k, NOFACE, ix, jx, kx, gd)] += q;
          va -= q * Sy;
        }
        
        if ( GET_FACE_BC(s, BC_FACE_B) == odr )
        {
          q = lmd * ppz * (sf - t);
          d_qbc[_F_IDX_S4DEX(Z_minus, i, j, k, NOFACE, ix, jx, kx, gd)] += q;
          va += q * Sz;
        }
        
        if ( GET_FACE_BC(s, BC_FACE_T) == odr )
        {
          q = lmd * ppz * (t - sf);
          d_qbc[_F_IDX_S4DEX(Z_plus, i, j, k, NOFACE, ix, jx, kx, gd)] += q;
          va -= q * Sz;
        }
      }
    }
  }
  
  
  if ( numProc > 1 )
  {
		REAL_TYPE tmp = va;
    if ( paraMngr->Allreduce(&tmp, &va, 1, MPI_SUM) != CPM_SUCCESS )
    {
      Hostonly_ printf("Allreduce Error\n");
      Exit(0);
    }
	}
  
  return va;
}



// #################################################################
/**
 * @brief 内部領域のOutflowの境界条件処理
 * @retval 無次元の単位時間あたりの熱量（有次元では[W]）
 * @param [in,out] d_ws  内部エネルギー増分
 * @param [in]     d_cdf BCindex C
 * @param [in]     n     コンポーネントリストのインデクス
 * @param [in]     d_v   速度
 * @param [in]     d_ie  内部エネルギー
 * @param [in]     v00   参照速度
 * @note 
   - モニタ量の熱量va(W)は系に対する流入量なので，基準温度に対する熱量
   - @todo 流出速度はモニター値を利用
 */
REAL_TYPE SetBC3D::psIbcOutflow (REAL_TYPE* d_ws, const int* d_cdf, const int n, const REAL_TYPE* d_v, const REAL_TYPE* d_ie, const REAL_TYPE* v00)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  int odr= n;
  int st[3], ed[3];
  
  REAL_TYPE va=0.0;
  REAL_TYPE rx = 1.0/pitch[0];
  REAL_TYPE ry = 1.0/pitch[1];
  REAL_TYPE rz = 1.0/pitch[2];
  REAL_TYPE u_ref = v00[1];
  REAL_TYPE v_ref = v00[2];
  REAL_TYPE w_ref = v00[3];
  
  REAL_TYPE Sx = S_x;
  REAL_TYPE Sy = S_y;
  REAL_TYPE Sz = S_z;
  
  cmp[n].getBbox(st, ed);
	
#pragma omp parallel for firstprivate(ix, jx, kx, gd, st, ed, odr, u_ref, v_ref, w_ref, Sx, Sy ,Sz, rx, ry, rz) \
schedule(static) reduction(+:va)
  
  for (int k=st[2]; k<=ed[2]; k++) {
    for (int j=st[1]; j<=ed[1]; j++) {
      for (int i=st[0]; i<=ed[0]; i++) {
        
        // 境界条件以外の寄与をゼロにしておく
        REAL_TYPE f_e = 0.0;
        REAL_TYPE f_w = 0.0;
        REAL_TYPE f_n = 0.0;
        REAL_TYPE f_s = 0.0;
        REAL_TYPE f_t = 0.0;
        REAL_TYPE f_b = 0.0;
        REAL_TYPE c=0.0;
        
        size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        int s = d_cdf[m];
        REAL_TYPE t_p = d_ie[m];
        
        if ( GET_FACE_BC(s, BC_FACE_W) == odr )
        {
          size_t m_0 = _F_IDX_V3D(i,   j, k, 0, ix, jx, kx, gd);
          size_t m_w = _F_IDX_V3D(i-1, j, k, 0, ix, jx, kx, gd);
          c = 0.5*(d_v[m_w]+d_v[m_0]) - u_ref;
          if ( c>0.0 ) c=0.0;
          f_w = c*t_p;
          va += f_w * Sx;
        }
        
        if ( GET_FACE_BC(s, BC_FACE_E) == odr )
        {
          size_t m_0 = _F_IDX_V3D(i,   j, k, 0, ix, jx, kx, gd);
          size_t m_e = _F_IDX_V3D(i+1, j, k, 0, ix, jx, kx, gd);
          c = 0.5*(d_v[m_e]+d_v[m_0]) - u_ref;
          if ( c<0.0 ) c=0.0;
          f_e = c*t_p;
          va += f_e * Sx;
        }
        
        if ( GET_FACE_BC(s, BC_FACE_S) == odr )
        {
          size_t m_0 = _F_IDX_V3D(i, j,   k, 1, ix, jx, kx, gd);
          size_t m_s = _F_IDX_V3D(i, j-1, k, 1, ix, jx, kx, gd);
          c = 0.5*(d_v[m_s]+d_v[m_0]) - v_ref;
          if ( c>0.0 ) c=0.0;
          f_s = c*t_p;
          va += f_s * Sy;
        }
        
        if ( GET_FACE_BC(s, BC_FACE_N) == odr )
        {
          size_t m_0 = _F_IDX_V3D(i, j,   k, 1, ix, jx, kx, gd);
          size_t m_n = _F_IDX_V3D(i, j+1, k, 1, ix, jx, kx, gd);
          c = 0.5*(d_v[m_n]+d_v[m_0]) - v_ref;
          if ( c<0.0 ) c=0.0;
          f_n = c*t_p;
          va += f_n * Sy;
        }
        
        if ( GET_FACE_BC(s, BC_FACE_B) == odr )
        {
          size_t m_0 = _F_IDX_V3D(i, j, k,   2, ix, jx, kx, gd);
          size_t m_b = _F_IDX_V3D(i, j, k-1, 2, ix, jx, kx, gd);
          c = 0.5*(d_v[m_b]+d_v[m_0]) - w_ref;
          if ( c>0.0 ) c=0.0;
          f_b = c*t_p;
          va += f_b * Sz;
        }
        
        if ( GET_FACE_BC(s, BC_FACE_T) == odr )
        {
          size_t m_0 = _F_IDX_V3D(i, j, k,   2, ix, jx, kx, gd);
          size_t m_t = _F_IDX_V3D(i, j, k+1, 2, ix, jx, kx, gd);
          c = 0.5*(d_v[m_t]+d_v[m_0]) - w_ref;
          if ( c<0.0 ) c=0.0;
          f_t = c*t_p;
          va += f_t * Sz;
        }
        
        d_ws[m] -= (f_e - f_w)*rx + (f_n - f_s)*ry + (f_t - f_b)*rz; // Outflow指定セルは，必ずFセルなのでマスク不要
      }
    }
  }
  
  
  if ( numProc > 1 )
  {
		REAL_TYPE tmp = va;
		if ( paraMngr->Allreduce(&tmp, &va, 1, MPI_SUM) != CPM_SUCCESS ) Exit(0);
	}
  
  return va;
}



// #################################################################
/**
 * @brief 内部領域の速度と温度の指定境界条件
 * @retval 無次元の単位時間あたりの熱量（有次元では[W]）
 * @param [in,out] d_ws   内部エネルギー増分
 * @param [in]     d_cdf  BCindex C
 * @param [in]     n      コンポーネントリストのインデクス
 * @param [in]     v00    参照速度
 * @param [in]     vec    指定ベクトル
 * @note モニタ量の熱量va(無次元)は系に対する流入量なので，基準温度に対する熱量
 */
REAL_TYPE SetBC3D::psIbcSpecVH (REAL_TYPE* d_ws, const int* d_cdf, const int n, const REAL_TYPE v00, const REAL_TYPE* vec)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  int odr= n;
  int st[3], ed[3];
  
  REAL_TYPE va = 0.0;
  REAL_TYPE rx = 1.0/pitch[0];
  REAL_TYPE ry = 1.0/pitch[1];
  REAL_TYPE rz = 1.0/pitch[2];
  
  REAL_TYPE rho = mat[n].P[p_density] / rho_0;
  REAL_TYPE cp  = mat[n].P[p_specific_heat] / cp_0;
  REAL_TYPE tmp = FBUtility::convTempD2ND(cmp[n].getTemp(), BaseTemp, DiffTemp);
  REAL_TYPE ref = rho * cp * tmp;
  
  REAL_TYPE hu  = vec[0]*ref;
  REAL_TYPE hv  = vec[1]*ref;
  REAL_TYPE hw  = vec[2]*ref;
  
  REAL_TYPE Sx = S_x;
  REAL_TYPE Sy = S_y;
  REAL_TYPE Sz = S_z;
  
  cmp[n].getBbox(st, ed);

#pragma omp parallel for firstprivate(ix, jx, kx, gd, st, ed, hu, hv, hw, odr, Sx, Sy ,Sz, rx, ry, rz) \
            schedule(static) reduction(+:va)
  
  for (int k=st[2]; k<=ed[2]; k++) {
    for (int j=st[1]; j<=ed[1]; j++) {
      for (int i=st[0]; i<=ed[0]; i++) {
        
        // 境界条件以外の寄与をゼロにしておく
        REAL_TYPE f_e = 0.0; 
        REAL_TYPE f_w = 0.0;
        REAL_TYPE f_n = 0.0;
        REAL_TYPE f_s = 0.0;
        REAL_TYPE f_t = 0.0;
        REAL_TYPE f_b = 0.0;
        
        size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        int s = d_cdf[m];
        
        if ( GET_FACE_BC(s, BC_FACE_W) == odr )
        {
          f_w = hu;
          va += hu * Sx;
        }
        
        if ( GET_FACE_BC(s, BC_FACE_E) == odr )
        {
          f_e = hu;
          va += hu * Sx;
        }
        
        if ( GET_FACE_BC(s, BC_FACE_S) == odr )
        {
          f_s = hv;
          va += hv * Sy;
        }
        
        if ( GET_FACE_BC(s, BC_FACE_N) == odr )
        {
          f_n = hv;
          va += hv * Sy;
        }
        
        if ( GET_FACE_BC(s, BC_FACE_B) == odr )
        {
          f_b = hw;
          va += hw * Sz;
        }
        
        if ( GET_FACE_BC(s, BC_FACE_T) == odr )
        {
          f_t = hw;
          va += hw * Sz;
        }
        
        d_ws[m] -= (f_e - f_w)*rx + (f_n - f_s)*ry + (f_t - f_b)*rz; // 流入境界は必ずFセルなのでマスク不要
      }
    }
  }
  
  if ( numProc > 1 )
  {
		REAL_TYPE tmp = va;
    if ( paraMngr->Allreduce(&tmp, &va, 1, MPI_SUM) != CPM_SUCCESS ) Exit(0);
	}
  
  return va;
}



// #################################################################
/**
 * @brief 熱伝達境界条件タイプS
 * @retval 無次元の単位時間あたりの熱量（有次元では[W]）
 * @param [in,out] d_qbc   境界条件熱流束
 * @param [in]     d_cdf   BCindex C
 * @param [in]     d_bcd   BCindex B
 * @param [in]     n       境界条件コンポーネントのエントリ番号
 * @param [in]     d_ie0   n時刻の内部エネルギー
 * @note 熱流束は加算（他の条件と合成）
 */
REAL_TYPE SetBC3D::psIbcTransferS (REAL_TYPE* d_qbc, const int* d_cdf, const int* d_bcd, const int n, const REAL_TYPE* d_ie0)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  int odr= n;
  
  REAL_TYPE va = 0.0;
  int st[3], ed[3];
  
  REAL_TYPE ht = cmp[n].getCoefHT() / (RefV*rho_0*cp_0);                        // RefMediumで指定される物性値
  REAL_TYPE bt = FBUtility::convTempD2ND(cmp[n].getTemp(), BaseTemp, DiffTemp); // 保持されている温度はCelsius

  REAL_TYPE Sx = S_x;
  REAL_TYPE Sy = S_y;
  REAL_TYPE Sz = S_z;
  
  cmp[n].getBbox(st, ed);
  
#pragma omp parallel for firstprivate(ix, jx, kx, gd, st, ed, odr, ht, bt, Sx, Sy, Sz) \
schedule(static) reduction(+:va)
  for (int k=st[2]; k<=ed[2]; k++) {
    for (int j=st[1]; j<=ed[1]; j++) {
      for (int i=st[0]; i<=ed[0]; i++) {
        
        size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        int s = d_cdf[m];
        REAL_TYPE q = 0.0;
        int l = d_bcd[m] & MASK_5;
        REAL_TYPE rho = mtbl[3*l+0];
        REAL_TYPE cp  = mtbl[3*l+1];
        REAL_TYPE t = d_ie0[m] / (rho * cp);
        
        if ( GET_FACE_BC(s, BC_FACE_W) == odr )
        {
          q = ht * (bt - t);
          d_qbc[_F_IDX_S4DEX(X_minus, i, j, k, NOFACE, ix, jx, kx, gd)] += q;
          va += q * Sx; // マイナス面の正の値は流入
        }
        
        if ( GET_FACE_BC(s, BC_FACE_E) == odr )
        {
          q = ht * (t - bt);
          d_qbc[_F_IDX_S4DEX(X_plus, i, j, k, NOFACE, ix, jx, kx, gd)] += q;
          va -= q * Sx; // プラス面の正の値は流出
        }
        
        if ( GET_FACE_BC(s, BC_FACE_S) == odr )
        {
          q = ht * (bt - t);
          d_qbc[_F_IDX_S4DEX(Y_minus, i, j, k, NOFACE, ix, jx, kx, gd)] += q;
          va += q * Sy;
        }
        
        if ( GET_FACE_BC(s, BC_FACE_N) == odr )
        {
          q = ht * (t - bt);
          d_qbc[_F_IDX_S4DEX(Y_plus, i, j, k, NOFACE, ix, jx, kx, gd)] += q;
          va -= q * Sy;
        }
        
        if ( GET_FACE_BC(s, BC_FACE_B) == odr )
        {
          q = ht * (bt - t);
          d_qbc[_F_IDX_S4DEX(Z_minus, i, j, k, NOFACE, ix, jx, kx, gd)] += q;
          va += q * Sz;
        }
        
        if ( GET_FACE_BC(s, BC_FACE_T) == odr )
        {
          q = ht * (t - bt);
          d_qbc[_F_IDX_S4DEX(Z_plus, i, j, k, NOFACE, ix, jx, kx, gd)] += q;
          va -= q * Sz;
        }
      }
    }
  }

  
  if ( numProc > 1 )
  {
		REAL_TYPE tmp = va;
    if ( paraMngr->Allreduce(&tmp, &va, 1, MPI_SUM) != CPM_SUCCESS ) Exit(0);
	}
  
  return va;
}



// #################################################################
/**
 * @brief 熱伝達境界条件タイプSF（強制対流)
 * @retval 無次元の単位時間あたりの熱量（有次元では[W]）
 * @param [in,out] d_qbc 境界条件熱流束
 * @param [in]     d_cdf BCindex C
 * @param [in]     d_bcd BCindex B
 * @param [in]     n     境界条件コンポーネントのエントリ番号
 * @param [in]     d_ie0 n時刻の内部エネルギー
 * @note 熱流束は加算（他の条件と合成）
 */
REAL_TYPE SetBC3D::psIbcTransferSF (REAL_TYPE* d_qbc, const int* d_cdf, const int* d_bcd, const int n, const REAL_TYPE* d_ie0)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  int odr= n;
  int st[3], ed[3];
  
  REAL_TYPE va = 0.0;
  REAL_TYPE a1 = cmp[n].ca[CompoList::alpha];
  REAL_TYPE b1 = cmp[n].ca[CompoList::beta];
  REAL_TYPE c1 = cmp[n].ca[CompoList::gamma];
  REAL_TYPE ht = lambda_0 / (RefV*RefL*rho_0*cp_0) * a1*pow(Reynolds, b1) * pow(Prandtl, c1); // RefMediumで指定される物性値
  REAL_TYPE sf = FBUtility::convTempD2ND(cmp[n].getTemp(), BaseTemp, DiffTemp);               // 保持されている温度はCelsius
  
  REAL_TYPE Sx = S_x;
  REAL_TYPE Sy = S_y;
  REAL_TYPE Sz = S_z;
  
  cmp[n].getBbox(st, ed);
  
#pragma omp parallel for firstprivate(ix, jx, kx, gd, st, ed, odr, ht, sf, Sx, Sy, Sz) \
schedule(static) reduction(+:va)
  for (int k=st[2]; k<=ed[2]; k++) {
    for (int j=st[1]; j<=ed[1]; j++) {
      for (int i=st[0]; i<=ed[0]; i++) {
        
        size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        int s = d_cdf[m];
        int l = d_bcd[m] & MASK_5;
        REAL_TYPE rho = mtbl[3*l+0];
        REAL_TYPE cp  = mtbl[3*l+1];
        REAL_TYPE t = d_ie0[m] / (rho * cp);
        REAL_TYPE q = 0.0;
        
        if ( GET_FACE_BC(s, BC_FACE_W) == odr )
        {
          q = ht * (sf - t);                                 // 基準温度との温度差
          d_qbc[_F_IDX_S4DEX(X_minus, i, j, k, NOFACE, ix, jx, kx, gd)] += q;
          va += q * Sx;                                                 // マイナス面の正の値は流入
        }
        
        if ( GET_FACE_BC(s, BC_FACE_E) == odr )
        {
          q = ht * (t - sf);
          d_qbc[_F_IDX_S4DEX(X_plus, i, j, k, NOFACE, ix, jx, kx, gd)] += q;
          va -= q * Sx;                                                 // プラス面の正の値は流出
        }
        
        if ( GET_FACE_BC(s, BC_FACE_S) == odr )
        {
          q = ht * (sf - t);
          d_qbc[_F_IDX_S4DEX(Y_minus, i, j, k, NOFACE, ix, jx, kx, gd)] += q;
          va += q * Sy;
        }
        
        if ( GET_FACE_BC(s, BC_FACE_N) == odr )
        {
          q = ht * (t - sf);
          d_qbc[_F_IDX_S4DEX(Y_plus, i, j, k, NOFACE, ix, jx, kx, gd)] += q;
          va -= q * Sy;
        }
        
        if ( GET_FACE_BC(s, BC_FACE_B) == odr )
        {
          q = ht * (sf - t);
          d_qbc[_F_IDX_S4DEX(Z_minus, i, j, k, NOFACE, ix, jx, kx, gd)] += q;
          va += q * Sz;
        }
        
        if ( GET_FACE_BC(s, BC_FACE_T) == odr )
        {
          q = ht * (t - sf);
          d_qbc[_F_IDX_S4DEX(Z_plus, i, j, k, NOFACE, ix, jx, kx, gd)] += q;
          va -= q * Sz;
        }
      }
    }
  }
  
  if ( numProc > 1 )
  {
    REAL_TYPE tmp = va;
    if ( paraMngr->Allreduce(&tmp, &va, 1, MPI_SUM) != CPM_SUCCESS )
    {
      Hostonly_ printf("Allreduce Error\n");
      Exit(0);
    }
  }
  
  return va;
}


// #################################################################
/**
 * @brief 熱伝達境界条件タイプSN（自然対流）
 * @retval 無次元の単位時間あたりの熱量（有次元では[W]）
 * @param [out]    d_qbc  境界条件熱流束
 * @param [in]     d_cdf  BCindex C
 * @param [in]     d_bcd  BCindex B
 * @param [in]     n      境界条件コンポーネントのエントリ番号
 * @param [in]     d_ie0  n時刻の内部エネルギー
 * @note
 *    - 熱流束は加算（他の条件と合成）
 *    - pセルは流体セル
 */
REAL_TYPE SetBC3D::psIbcTransferSN (REAL_TYPE* d_qbc, const int* d_cdf, const int* d_bcd, const int n, const REAL_TYPE* d_ie0)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  int odr= n;
  
  REAL_TYPE va = 0.0;
  int st[3], ed[3];
  
  REAL_TYPE a1 = cmp[n].ca[CompoList::vert_laminar_alpha];    // vertical, upper, laminar
  REAL_TYPE a2 = cmp[n].ca[CompoList::vert_turbulent_alpha];  // vertical, upper, turbulent
  REAL_TYPE a3 = cmp[n].cb[CompoList::lower_laminar_alpha];   // lower, laminar
  REAL_TYPE a4 = cmp[n].cb[CompoList::lower_turbulent_alpha]; // lower, turbulent
  
  REAL_TYPE b1 = cmp[n].ca[CompoList::vert_laminar_beta];     // vertical, upper, laminar
  REAL_TYPE b2 = cmp[n].ca[CompoList::vert_turbulent_beta];   // vertical, upper, turbulent
  REAL_TYPE b3 = cmp[n].cb[CompoList::lower_laminar_beta];    // lower, laminar
  REAL_TYPE b4 = cmp[n].cb[CompoList::lower_turbulent_beta];  // lower, turbulent
  
  REAL_TYPE c1 = cmp[n].ca[CompoList::vert_Ra_critial];       // vertical, upper, Ra_c
  REAL_TYPE c2 = cmp[n].cb[CompoList::lower_Ra_critial];      // lower, Ra_c
  REAL_TYPE ht = lambda_0 / (RefV*RefL*rho_0*cp_0);           // RefMediumで指定される物性値
  
  REAL_TYPE ht1 = ( Rayleigh < c1 ) ? ht*a1*pow(Rayleigh, b1) : ht*a2*pow(Rayleigh, b2); // vertical, upper
  REAL_TYPE ht3 = ( Rayleigh < c2 ) ? ht*a3*pow(Rayleigh, b3) : ht*a4*pow(Rayleigh, b4); // lower
  REAL_TYPE sf = FBUtility::convTempD2ND(cmp[n].getTemp(), BaseTemp, DiffTemp); // 保持されている温度はCelsius
  
  REAL_TYPE Sx = S_x;
  REAL_TYPE Sy = S_y;
  REAL_TYPE Sz = S_z;
  
  cmp[n].getBbox(st, ed);
  
  
#pragma omp parallel for firstprivate(ix, jx, kx, gd, st, ed, odr, ht1, ht3, sf, Sx, Sy, Sz) \
schedule(static) reduction(+:va)
  
  for (int k=st[2]; k<=ed[2]; k++) {
    for (int j=st[1]; j<=ed[1]; j++) {
      for (int i=st[0]; i<=ed[0]; i++) {
        
        size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        int s = d_cdf[m];
        int l = d_bcd[m] & MASK_5;
        REAL_TYPE rho = mtbl[3*l+0];
        REAL_TYPE cp  = mtbl[3*l+1];
        REAL_TYPE t = d_ie0[m] / (rho * cp);
        REAL_TYPE q = 0.0;
        
        if ( GET_FACE_BC(s, BC_FACE_W) == odr )
        {
          q = ht1 * (sf - t); // 基準温度との温度差
          d_qbc[_F_IDX_S4DEX(X_minus, i, j, k, NOFACE, ix, jx, kx, gd)] += q;
          va += q * Sx; // マイナス面の正の値は流入
        }
        
        if ( GET_FACE_BC(s, BC_FACE_E) == odr )
        {
          q = ht1 * (t - sf);
          d_qbc[_F_IDX_S4DEX(X_plus, i, j, k, NOFACE, ix, jx, kx, gd)] += q;
          va -= q * Sx; // プラス面の正の値は流出
        }
        
        if ( GET_FACE_BC(s, BC_FACE_S) == odr )
        {
          q = ht1 * (sf - t);
          d_qbc[_F_IDX_S4DEX(Y_minus, i, j, k, NOFACE, ix, jx, kx, gd)] += q;
          va += q * Sy;
        }
        
        if ( GET_FACE_BC(s, BC_FACE_N) == odr )
        {
          q = ht1 * (t - sf);
          d_qbc[_F_IDX_S4DEX(Y_plus, i, j, k, NOFACE, ix, jx, kx, gd)] += q;
          va -= q * Sy;
        }
        
        if ( GET_FACE_BC(s, BC_FACE_B) == odr )
        {
          q = ht3 * (sf - t);
          d_qbc[_F_IDX_S4DEX(Z_minus, i, j, k, NOFACE, ix, jx, kx, gd)] += q;
          va += q * Sz;
        }
        
        if ( GET_FACE_BC(s, BC_FACE_T) == odr )
        {
          q = ht1 * (t - sf);
          d_qbc[_F_IDX_S4DEX(Z_plus, i, j, k, NOFACE, ix, jx, kx, gd)] += q;
          va -= q * Sz;
        }
      }
    }
  }
  
  if ( numProc > 1 )
  {
    REAL_TYPE tmp = va;
    if ( paraMngr->Allreduce(&tmp, &va, 1, MPI_SUM) != CPM_SUCCESS ) Exit(0);
  }
  
  return va;
}


// #################################################################
/**
 * @brief 外部領域のOutflow, In_out, TractionFreeの境界条件処理
 * @retval 無次元の単位時間あたりの熱量（有次元では[W]）
 * @param [in,out] d_ws  内部エネルギー増分
 * @param [in]     d_bcd BCindex B
 * @param [in]     face  外部境界面番号
 * @param [in]     d_vf  セルフェイス速度
 * @param [in]     d_ie  n時刻の内部エネルギー
 * @param [in]     v00   参照速度
 * @note モニタ量の熱量 va は系に対する流入量なので，基準温度に対する熱量
 */
REAL_TYPE SetBC3D::psObcFree(REAL_TYPE* d_ws, const int* d_bcd, const int face, const REAL_TYPE* d_vf, const REAL_TYPE* d_ie, const REAL_TYPE* v00)
{
  if ( nID[face] >= 0 ) return 0.0;
  
  
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
  REAL_TYPE va = 0.0;
  REAL_TYPE rx = 1.0/pitch[0];
  REAL_TYPE ry = 1.0/pitch[1];
  REAL_TYPE rz = 1.0/pitch[2];
  REAL_TYPE u_ref = v00[1];
  REAL_TYPE v_ref = v00[2];
  REAL_TYPE w_ref = v00[3];
  
  REAL_TYPE Sx = S_x;
  REAL_TYPE Sy = S_y;
  REAL_TYPE Sz = S_z;
  
  switch (face)
  {
    case X_minus:

#pragma omp parallel for firstprivate(ix, jx, kx, gd, rx, u_ref, Sx) schedule(static) reduction(+:va)
      for (int k=1; k<=kx; k++) {
        for (int j=1; j<=jx; j++) {
          size_t m = _F_IDX_S3D(1, j, k, ix, jx, kx, gd);
          int s = d_bcd[m];
          size_t m_w = _F_IDX_V3D(0, j, k, 0, ix, jx, kx, gd);
          REAL_TYPE c = d_vf[m_w] - u_ref;
          REAL_TYPE t_p = (c < 0.0) ? d_ie[m] : 0.0;           // 流出の場合には内部の値，流入の場合には断熱
          REAL_TYPE ff = c*t_p;
          va += ff * Sx;
          d_ws[m] += ff*rx*GET_SHIFT_F(s, STATE_BIT);
        }
      }
      break;
      
    case X_plus:

#pragma omp parallel for firstprivate(ix, jx, kx, gd, rx, u_ref, Sx) schedule(static) reduction(+:va)
      for (int k=1; k<=kx; k++) {
        for (int j=1; j<=jx; j++) {
          size_t m = _F_IDX_S3D(ix, j, k, ix, jx, kx, gd);
          int s = d_bcd[m];
          size_t m_e = _F_IDX_V3D(ix, j, k, 0, ix, jx, kx, gd);
          REAL_TYPE c = d_vf[m_e] - u_ref;
          REAL_TYPE t_p = (c > 0.0) ? d_ie[m] : 0.0;
          REAL_TYPE ff = c*t_p;
          va += ff * Sx;
          d_ws[m] -= ff*rx*GET_SHIFT_F(s, STATE_BIT);
        }
      }
      break;
      
    case Y_minus:
      
#pragma omp parallel for firstprivate(ix, jx, kx, gd, ry, v_ref, Sy) schedule(static) reduction(+:va)
      for (int k=1; k<=kx; k++) {
        for (int i=1; i<=ix; i++) {
          size_t m = _F_IDX_S3D(i, 1, k, ix, jx, kx, gd);
          int s = d_bcd[m];
          size_t m_s = _F_IDX_V3D(i, 0, k, 1, ix, jx, kx, gd);
          REAL_TYPE c = d_vf[m_s] - v_ref;
          REAL_TYPE t_p = (c < 0.0) ? d_ie[m] : 0.0;
          REAL_TYPE ff = c*t_p;
          va += ff * Sy;
          d_ws[m] += ff*ry*GET_SHIFT_F(s, STATE_BIT);
        }
      }
      break;
      
    case Y_plus:
      
#pragma omp parallel for firstprivate(ix, jx, kx, gd, ry, v_ref, Sy) schedule(static) reduction(+:va)
      for (int k=1; k<=kx; k++) {
        for (int i=1; i<=ix; i++) {
          size_t m = _F_IDX_S3D(i, jx, k, ix, jx, kx, gd);
          int s = d_bcd[m];
          size_t m_n = _F_IDX_V3D(i, jx, k, 1, ix, jx, kx, gd);
          REAL_TYPE c = d_vf[m_n] - v_ref;
          REAL_TYPE t_p = (c > 0.0) ? d_ie[m] : 0.0;
          REAL_TYPE ff = c*t_p;
          va += ff * Sy;
          d_ws[m] -= ff*ry*GET_SHIFT_F(s, STATE_BIT);
        }
      }
      break;
      
    case Z_minus:
      
#pragma omp parallel for firstprivate(ix, jx, kx, gd, rz, w_ref, Sz) schedule(static) reduction(+:va)
      for (int j=1; j<=jx; j++) {
        for (int i=1; i<=ix; i++) {
          size_t m = _F_IDX_S3D(i, j, 1, ix, jx, kx, gd);
          int s = d_bcd[m];
          size_t m_b = _F_IDX_V3D(i, j, 0, 2, ix, jx, kx, gd);
          REAL_TYPE c = d_vf[m_b] - w_ref;
          REAL_TYPE t_p = (c < 0.0) ? d_ie[m] : 0.0;
          REAL_TYPE ff = c*t_p;
          va += ff * Sz;
          d_ws[m] += ff*rz*GET_SHIFT_F(s, STATE_BIT);
        }
      }
      break;
      
    case Z_plus:
      
#pragma omp parallel for firstprivate(ix, jx, kx, gd, rz, w_ref, Sz) schedule(static) reduction(+:va)
      for (int j=1; j<=jx; j++) {
        for (int i=1; i<=ix; i++) {
          size_t m = _F_IDX_S3D(i, j, kx, ix, jx, kx, gd);
          int s = d_bcd[m];
          size_t m_t = _F_IDX_V3D(i, j, kx, 2, ix, jx, kx, gd);
          REAL_TYPE c = d_vf[m_t] - w_ref;
          REAL_TYPE t_p = (c > 0.0) ? d_ie[m] : 0.0;
          REAL_TYPE ff = c*t_p;
          va += ff * Sz;
          d_ws[m] -= ff*rz*GET_SHIFT_F(s, STATE_BIT);
        }
      }
      break;
  }
  
  return va;
}


// #################################################################
/**
 * @brief 外部領域の熱流束指定の境界条件処理
 * @retval 無次元の単位時間あたりの熱量（有次元では[W]）
 * @param [in,out] d_qbc 境界条件熱流束
 * @param [in]     face  外部境界面番号
 * @note
    - モニタ量の熱量 va は系に対する流入出量
    - 流入を正にとる
 */
REAL_TYPE SetBC3D::psObcHeatflux(REAL_TYPE* d_qbc, const int face)
{
  if ( nID[face] >= 0 ) return 0.0;
  
  
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
  REAL_TYPE va = 0.0;
  REAL_TYPE q = obc[face].getHeatflux() / (RefV*DiffTemp*rho_0*cp_0); // [W/m^2]を無次元化
  
  REAL_TYPE Sx = S_x;
  REAL_TYPE Sy = S_y;
  REAL_TYPE Sz = S_z;
  
  switch (face)
  {
    case X_minus:
#pragma omp parallel for firstprivate(ix, jx, kx, gd, q, Sx) schedule(static) reduction(+:va)
      for (int k=1; k<=kx; k++) {
        for (int j=1; j<=jx; j++) {
          d_qbc[_F_IDX_S4DEX(X_minus, 1, j, k, NOFACE, ix, jx, kx, gd)] += q;
          va += q * Sx;
        }
      }
      break;
      
    case X_plus:
#pragma omp parallel for firstprivate(ix, jx, kx, gd, q, Sx) schedule(static) reduction(+:va)
      for (int k=1; k<=kx; k++) {
        for (int j=1; j<=jx; j++) {
          d_qbc[_F_IDX_S4DEX(X_plus, ix, j, k, NOFACE, ix, jx, kx, gd)] -= q; // プラス面の正の値は流出
          va -= q * Sx;
        }
      }
      break;
      
    case Y_minus:
#pragma omp parallel for firstprivate(ix, jx, kx, gd, q, Sy) schedule(static) reduction(+:va)
      for (int k=1; k<=kx; k++) {
        for (int i=1; i<=ix; i++) {
          d_qbc[_F_IDX_S4DEX(Y_minus, i, 1, k, NOFACE, ix, jx, kx, gd)] += q;
          va += q * Sy;
        }
      }
      break;
      
    case Y_plus:
#pragma omp parallel for firstprivate(ix, jx, kx, gd, q, Sy) schedule(static) reduction(+:va)
      for (int k=1; k<=kx; k++) {
        for (int i=1; i<=ix; i++) {
          d_qbc[_F_IDX_S4DEX(Y_plus, i, jx, k, NOFACE, ix, jx, kx, gd)] -= q;
          va -= q * Sy;
        }
      }
      break;
      
    case Z_minus:
#pragma omp parallel for firstprivate(ix, jx, kx, gd, q, Sz) schedule(static) reduction(+:va)
      for (int j=1; j<=jx; j++) {
        for (int i=1; i<=ix; i++) {
          d_qbc[_F_IDX_S4DEX(Z_minus, i, j, 1, NOFACE, ix, jx, kx, gd)] += q;
          va += q * Sz;
        }
      }
      break;
      
    case Z_plus:
#pragma omp parallel for firstprivate(ix, jx, kx, gd, q, Sz) schedule(static) reduction(+:va)
      for (int j=1; j<=jx; j++) {
        for (int i=1; i<=ix; i++) {
          d_qbc[_F_IDX_S4DEX(Z_plus, i, j, kx, NOFACE, ix, jx, kx, gd)] -= q;
          va -= q * Sz;
        }
      }
      break;
  }
  
  return va;
}


// #################################################################
/**
 * @brief 外部領域の熱伝達境界の境界条件処理 (TypeS)
 * @retval 無次元の単位時間あたりの熱量（有次元では[W]）
 * @param [in,out] d_qbc  境界条件熱流束
 * @param [in]     face   外部境界面番号
 * @param [in]     d_bcd  BCindex B
 * @param [out]    d_ie   n+1時刻の内部エネルギー
 * @param [in]     d_ie0  n時刻の内部エネルギー
 * @note
    - モニタ量の熱量 va は系に対する流入出量
    - 流入を正にとる
 */
REAL_TYPE SetBC3D::psObcHeatTransferS (REAL_TYPE* d_qbc, const int face, const int* d_bcd, REAL_TYPE* d_ie, const REAL_TYPE* d_ie0)
{
  if ( nID[face] >= 0 ) return 0.0;
  
  
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
  REAL_TYPE va = 0.0;
  REAL_TYPE ht = obc[face].getCoefHT() / (RefV*rho_0*cp_0);                        // 熱伝達係数
  REAL_TYPE bt = FBUtility::convTempD2ND(obc[face].getTemp(), BaseTemp, DiffTemp); // 温度．保持されている温度はCelsius
  
  int n = obc[face].getGuideMedium();
  REAL_TYPE rho_1 = mat[n].P[p_density] / rho_0;
  REAL_TYPE cp_1  = mat[n].P[p_specific_heat] / cp_0;
  REAL_TYPE ie  = rho_1 * cp_1 * bt;
  
  REAL_TYPE Sx = S_x;
  REAL_TYPE Sy = S_y;
  REAL_TYPE Sz = S_z;
  
  switch (face)
  {
    case X_minus:
#pragma omp parallel for firstprivate(ix, jx, kx, gd, ht, bt, ie, Sx) schedule(static) reduction(+:va)
      for (int k=1; k<=kx; k++) {
        for (int j=1; j<=jx; j++) {
          size_t m = _F_IDX_S3D(1, j, k, ix, jx, kx, gd);
          int l = d_bcd[m] & MASK_5;
          REAL_TYPE rho = mtbl[3*l+0];
          REAL_TYPE cp  = mtbl[3*l+1];
          REAL_TYPE t = d_ie0[m] / (rho * cp);
          REAL_TYPE q = ht * (bt - t);
          d_qbc[_F_IDX_S4DEX(X_minus, 1, j, k, NOFACE, ix, jx, kx, gd)] += q; // マイナス面の正の値は流入
          va += q * Sx;
          d_ie[_F_IDX_S3D(0, j, k, ix, jx, kx, gd)] = ie; // 隣接セルに代入
        }
      }
      break;
      
    case X_plus:
#pragma omp parallel for firstprivate(ix, jx, kx, gd, ht, bt, ie, Sx) schedule(static) reduction(+:va)
      for (int k=1; k<=kx; k++) {
        for (int j=1; j<=jx; j++) {
          size_t m = _F_IDX_S3D(ix, j, k, ix, jx, kx, gd);
          int l = d_bcd[m] & MASK_5;
          REAL_TYPE rho = mtbl[3*l+0];
          REAL_TYPE cp  = mtbl[3*l+1];
          REAL_TYPE t = d_ie0[m] / (rho * cp);
          REAL_TYPE q = ht * (t - bt);
          d_qbc[_F_IDX_S4DEX(X_plus, ix, j, k, NOFACE, ix, jx, kx, gd)] -= q; // プラス面の正の値は流出
          va -= q * Sx;
          d_ie[_F_IDX_S3D(ix+1, j, k, ix, jx, kx, gd)] = ie;
        }
      }
      break;
      
    case Y_minus:
#pragma omp parallel for firstprivate(ix, jx, kx, gd, ht, bt, ie, Sy) schedule(static) reduction(+:va)
      for (int k=1; k<=kx; k++) {
        for (int i=1; i<=ix; i++) {
          size_t m = _F_IDX_S3D(i, 1, k, ix, jx, kx, gd);
          int l = d_bcd[m] & MASK_5;
          REAL_TYPE rho = mtbl[3*l+0];
          REAL_TYPE cp  = mtbl[3*l+1];
          REAL_TYPE t = d_ie0[m] / (rho * cp);
          REAL_TYPE q = ht * (bt - t);
          d_qbc[_F_IDX_S4DEX(Y_minus, i, 1, k, NOFACE, ix, jx, kx, gd)] += q;
          va += q * Sy;
          d_ie[_F_IDX_S3D(i, 0, k, ix, jx, kx, gd)] = ie;
        }
      }
      break;
      
    case Y_plus:
#pragma omp parallel for firstprivate(ix, jx, kx, gd, ht, bt, ie, Sy) schedule(static) reduction(+:va)
      for (int k=1; k<=kx; k++) {
        for (int i=1; i<=ix; i++) {
          size_t m = _F_IDX_S3D(i, jx, k, ix, jx, kx, gd);
          int l = d_bcd[m] & MASK_5;
          REAL_TYPE rho = mtbl[3*l+0];
          REAL_TYPE cp  = mtbl[3*l+1];
          REAL_TYPE t = d_ie0[m] / (rho * cp);
          REAL_TYPE q = ht * (t - bt);
          d_qbc[_F_IDX_S4DEX(Y_plus, i, jx, k, NOFACE, ix, jx, kx, gd)] -= q;
          va -= q * Sy;
          d_ie[_F_IDX_S3D(i, jx+1, k, ix, jx, kx, gd)] = ie;
        }
      }
      break;
      
    case Z_minus:
#pragma omp parallel for firstprivate(ix, jx, kx, gd, ht, bt, ie, Sz) schedule(static) reduction(+:va)
      for (int j=1; j<=jx; j++) {
        for (int i=1; i<=ix; i++) {
          size_t m = _F_IDX_S3D(i, j, 1, ix, jx, kx, gd);
          int l = d_bcd[m] & MASK_5;
          REAL_TYPE rho = mtbl[3*l+0];
          REAL_TYPE cp  = mtbl[3*l+1];
          REAL_TYPE t = d_ie0[m] / (rho * cp);
          REAL_TYPE q = ht * (bt - t);
          d_qbc[_F_IDX_S4DEX(Z_minus, i, j, 1, NOFACE, ix, jx, kx, gd)] += q;
          va += q * Sz;
          d_ie[_F_IDX_S3D(i, j, 0, ix, jx, kx, gd)] = ie;
        }
      }
      break;
      
    case Z_plus:
#pragma omp parallel for firstprivate(ix, jx, kx, gd, ht, bt, ie, Sz) schedule(static) reduction(+:va)
      for (int j=1; j<=jx; j++) {
        for (int i=1; i<=ix; i++) {
          size_t m = _F_IDX_S3D(i, j, kx, ix, jx, kx, gd);
          int l = d_bcd[m] & MASK_5;
          REAL_TYPE rho = mtbl[3*l+0];
          REAL_TYPE cp  = mtbl[3*l+1];
          REAL_TYPE t = d_ie0[m] / (rho * cp);
          REAL_TYPE q = ht * (t - bt);
          d_qbc[_F_IDX_S4DEX(Z_plus, i, j, kx, NOFACE, ix, jx, kx, gd)] -= q;
          va -= q * Sz;
          d_ie[_F_IDX_S3D(i, j, kx+1, ix, jx, kx, gd)] = ie;
        }
      }
      break;
  }
  
  return va;
}


// #################################################################
/**
 * @brief 外部領域の熱伝達境界の境界条件処理 (Type SF)
 * @retval 無次元の単位時間あたりの熱量（有次元では[W]）
 * @param [in,out] d_qbc  境界条件熱流束
 * @param [in]     face   外部境界面番号
 * @param [in]     d_bcd  BCindex B
 * @param [out]    d_ie   n+1時刻の内部エネルギー
 * @param [in]     d_ie0  n時刻の内部エネルギー
 * @note
    - モニタ量の熱量va(-)は系に対する流入出量
    - 流入を正にとる
 */
REAL_TYPE SetBC3D::psObcHeatTransferSF(REAL_TYPE* d_qbc, const int face, const int* d_bcd, REAL_TYPE* d_ie, const REAL_TYPE* d_ie0)
{
  if ( nID[face] >= 0 ) return 0.0;
  
  
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
  REAL_TYPE va = 0.0;
  
  REAL_TYPE a1 = obc[face].ca[CompoList::alpha];
  REAL_TYPE b1 = obc[face].ca[CompoList::beta];
  REAL_TYPE c1 = obc[face].ca[CompoList::gamma];
  REAL_TYPE ht = lambda_0 / (RefV*RefL*rho_0*cp_0) * a1*pow(Reynolds, b1) * pow(Prandtl, c1);  // 熱伝達係数
  REAL_TYPE sf = FBUtility::convTempD2ND(obc[face].getTemp(), BaseTemp, DiffTemp);                 // 表面温度．保持されている温度はCelsius
  
  int n = obc[face].getGuideMedium();
  REAL_TYPE rho_1 = mat[n].P[p_density] / rho_0;
  REAL_TYPE cp_1  = mat[n].P[p_specific_heat] / cp_0;
  REAL_TYPE ie  = rho_1 * cp_1 * sf;
  
  REAL_TYPE Sx = S_x;
  REAL_TYPE Sy = S_y;
  REAL_TYPE Sz = S_z;
  
  switch (face)
  {
    case X_minus:
#pragma omp parallel for firstprivate(ix, jx, kx, gd, ht, sf, ie, Sx) schedule(static) reduction(+:va)
      for (int k=1; k<=kx; k++) {
        for (int j=1; j<=jx; j++) {
          size_t m = _F_IDX_S3D(1, j, k, ix, jx, kx, gd);
          int l = d_bcd[m] & MASK_5;
          REAL_TYPE rho = mtbl[3*l+0];
          REAL_TYPE cp  = mtbl[3*l+1];
          REAL_TYPE t = d_ie0[m] / (rho * cp);
          REAL_TYPE q = ht * (sf - t);             // 表面温度と隣接セルとの温度差
          d_qbc[_F_IDX_S4DEX(X_minus, 1, j, k, NOFACE, ix, jx, kx, gd)] += q;
          va += q * Sx;
          d_ie[_F_IDX_S3D(0, j, k, ix, jx, kx, gd)] = ie;
        }
      }
      break;
      
    case X_plus:
#pragma omp parallel for firstprivate(ix, jx, kx, gd, ht, sf, ie, Sx) schedule(static) reduction(+:va)
      for (int k=1; k<=kx; k++) {
        for (int j=1; j<=jx; j++) {
          size_t m = _F_IDX_S3D(ix, j, k, ix, jx, kx, gd);
          int l = d_bcd[m] & MASK_5;
          REAL_TYPE rho = mtbl[3*l+0];
          REAL_TYPE cp  = mtbl[3*l+1];
          REAL_TYPE t = d_ie0[m] / (rho * cp);
          REAL_TYPE q = ht * (t - sf);
          d_qbc[_F_IDX_S4DEX(X_plus, ix, j, k, NOFACE, ix, jx, kx, gd)] -= q;
          va -= q * Sx;                                                 // プラス面の正の値は流出
          d_ie[_F_IDX_S3D(ix+1, j, k, ix, jx, kx, gd)] = ie;
        }
      }
      break;
      
    case Y_minus:
#pragma omp parallel for firstprivate(ix, jx, kx, gd, ht, sf, ie, Sy) schedule(static) reduction(+:va)
      for (int k=1; k<=kx; k++) {
        for (int i=1; i<=ix; i++) {
          size_t m = _F_IDX_S3D(i, 1, k, ix, jx, kx, gd);
          int l = d_bcd[m] & MASK_5;
          REAL_TYPE rho = mtbl[3*l+0];
          REAL_TYPE cp  = mtbl[3*l+1];
          REAL_TYPE t = d_ie0[m] / (rho * cp);
          REAL_TYPE q = ht * (sf - t);
          d_qbc[_F_IDX_S4DEX(Y_minus, i, 1, k, NOFACE, ix, jx, kx, gd)] += q;
          va += q * Sy;
          d_ie[_F_IDX_S3D(i, 0, k, ix, jx, kx, gd)] = ie;
        }
      }
      break;
      
    case Y_plus:
#pragma omp parallel for firstprivate(ix, jx, kx, gd, ht, sf, ie, Sy) schedule(static) reduction(+:va)
      for (int k=1; k<=kx; k++) {
        for (int i=1; i<=ix; i++) {
          size_t m = _F_IDX_S3D(i, jx, k, ix, jx, kx, gd);
          int l = d_bcd[m] & MASK_5;
          REAL_TYPE rho = mtbl[3*l+0];
          REAL_TYPE cp  = mtbl[3*l+1];
          REAL_TYPE t = d_ie0[m] / (rho * cp);
          REAL_TYPE q = ht * (t - sf);
          d_qbc[_F_IDX_S4DEX(Y_plus, i, jx, k, NOFACE, ix, jx, kx, gd)] -= q;
          va -= q * Sy;
          d_ie[_F_IDX_S3D(i, jx+1, k, ix, jx, kx, gd)] = ie;
        }
      }
      break;
      
    case Z_minus:
#pragma omp parallel for firstprivate(ix, jx, kx, gd, ht, sf, ie, Sz) schedule(static) reduction(+:va)
      for (int j=1; j<=jx; j++) {
        for (int i=1; i<=ix; i++) {
          size_t m = _F_IDX_S3D(i, j, 1, ix, jx, kx, gd);
          int l = d_bcd[m] & MASK_5;
          REAL_TYPE rho = mtbl[3*l+0];
          REAL_TYPE cp  = mtbl[3*l+1];
          REAL_TYPE t = d_ie0[m] / (rho * cp);
          REAL_TYPE q = ht * (sf - t);
          d_qbc[_F_IDX_S4DEX(Z_minus, i, j, 1, NOFACE, ix, jx, kx, gd)] += q;
          va += q * Sz;
          d_ie[_F_IDX_S3D(i, j, 0, ix, jx, kx, gd)] = ie;
        }
      }
      break;
      
    case Z_plus:
#pragma omp parallel for firstprivate(ix, jx, kx, gd, ht, sf, ie, Sz) schedule(static) reduction(+:va)
      for (int j=1; j<=jx; j++) {
        for (int i=1; i<=ix; i++) {
          size_t m = _F_IDX_S3D(i, j, kx, ix, jx, kx, gd);
          int l = d_bcd[m] & MASK_5;
          REAL_TYPE rho = mtbl[3*l+0];
          REAL_TYPE cp  = mtbl[3*l+1];
          REAL_TYPE t = d_ie0[m] / (rho * cp);
          REAL_TYPE q = ht * (t - sf);
          d_qbc[_F_IDX_S4DEX(Z_plus, i, j, kx, NOFACE, ix, jx, kx, gd)] -= q;
          va -= q * Sz;
          d_ie[_F_IDX_S3D(i, j, kx+1, ix, jx, kx, gd)] = ie;
        }
      }
      break;
  }
  
  return va;
}


// #################################################################
/**
 * @brief 外部領域の熱伝達境界の境界条件処理 (Type SN)
 * @retval 無次元の単位時間あたりの熱量（有次元では[W]）
 * @param [in,out] d_qbc  境界条件熱流束
 * @param [in]     face   外部境界面番号
 * @param [in]     d_bcd  BCindex B
 * @param [out]    d_ie   n+1時刻の内部エネルギー
 * @param [in]     d_ie0  n時刻の内部エネルギー
 * @note
   - モニタ量の熱量 va は系に対する流入出量
   - 流入を正にとる
 */
REAL_TYPE SetBC3D::psObcHeatTransferSN(REAL_TYPE* d_qbc, const int face, const int* d_bcd, REAL_TYPE* d_ie, const REAL_TYPE* d_ie0)
{
  if ( nID[face] >= 0 ) return 0.0;
  
  
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
  REAL_TYPE va = 0.0;
  
  REAL_TYPE a1 = obc[face].ca[CompoList::vert_laminar_alpha];    // vertical, upper, laminar
  REAL_TYPE a2 = obc[face].ca[CompoList::vert_turbulent_alpha];  // vertical, upper, turbulent
  REAL_TYPE a3 = obc[face].cb[CompoList::lower_laminar_alpha];   // lower, laminar
  REAL_TYPE a4 = obc[face].cb[CompoList::lower_turbulent_alpha]; // lower, turbulent
  
  REAL_TYPE b1 = obc[face].ca[CompoList::vert_laminar_beta];     // vertical, upper, laminar
  REAL_TYPE b2 = obc[face].ca[CompoList::vert_turbulent_beta];   // vertical, upper, turbulent
  REAL_TYPE b3 = obc[face].cb[CompoList::lower_laminar_beta];    // lower, laminar
  REAL_TYPE b4 = obc[face].cb[CompoList::lower_turbulent_beta];  // lower, turbulent
  
  REAL_TYPE c1 = obc[face].ca[CompoList::vert_Ra_critial];       // vertical, upper, Ra_c
  REAL_TYPE c2 = obc[face].cb[CompoList::lower_Ra_critial];      // lower, Ra_c
  REAL_TYPE ht = lambda_0 / (RefV*RefL*rho_0*cp_0);              // RefMediumで指定される物性値
  
  REAL_TYPE ht1 = ( Rayleigh < c1 ) ? ht*a1*pow(Rayleigh, b1) : ht*a2*pow(Rayleigh, b2); // vertical, upper
  REAL_TYPE ht3 = ( Rayleigh < c2 ) ? ht*a3*pow(Rayleigh, b3) : ht*a4*pow(Rayleigh, b4); // lower
  REAL_TYPE sf = FBUtility::convTempD2ND(obc[face].getTemp(), BaseTemp, DiffTemp);           // 保持されている温度はCelsius
  
  int n = obc[face].getGuideMedium();
  REAL_TYPE rho_1 = mat[n].P[p_density] / rho_0;
  REAL_TYPE cp_1  = mat[n].P[p_specific_heat] / cp_0;
  REAL_TYPE ie  = rho_1 * cp_1 * sf;
  
  REAL_TYPE Sx = S_x;
  REAL_TYPE Sy = S_y;
  REAL_TYPE Sz = S_z;
  
  switch (face)
  {
    case X_minus:
#pragma omp parallel for firstprivate(ix, jx, kx, gd, ht1, sf, ie, Sx) schedule(static) reduction(+:va)
      for (int k=1; k<=kx; k++) {
        for (int j=1; j<=jx; j++) {
          size_t m = _F_IDX_S3D(1, j, k, ix, jx, kx, gd);
          int l = d_bcd[m] & MASK_5;
          REAL_TYPE rho = mtbl[3*l+0];
          REAL_TYPE cp  = mtbl[3*l+1];
          REAL_TYPE t = d_ie0[m] / (rho * cp);
          REAL_TYPE q = ht1 * (sf - t);                    // 基準温度との温度差
          d_qbc[_F_IDX_S4DEX(X_minus, 1, j, k, NOFACE, ix, jx, kx, gd)] += q;
          va += q * Sx;
          d_ie[_F_IDX_S3D(0, j, k, ix, jx, kx, gd)] = ie;
        }
      }
      break;
      
    case X_plus:
#pragma omp parallel for firstprivate(ix, jx, kx, gd, ht1, sf, ie, Sx) schedule(static) reduction(+:va)
      for (int k=1; k<=kx; k++) {
        for (int j=1; j<=jx; j++) {
          size_t m = _F_IDX_S3D(ix, j, k, ix, jx, kx, gd);
          int l = d_bcd[m] & MASK_5;
          REAL_TYPE rho = mtbl[3*l+0];
          REAL_TYPE cp  = mtbl[3*l+1];
          REAL_TYPE t = d_ie0[m] / (rho * cp);
          REAL_TYPE q = ht1 * (t - sf);
          d_qbc[_F_IDX_S4DEX(X_plus, ix, j, k, NOFACE, ix, jx, kx, gd)] += q;
          va -= q * Sx;                                               // プラス面の正の値は流出
          d_ie[_F_IDX_S3D(ix+1, j, k, ix, jx, kx, gd)] = ie;
        }
      }
      break;
      
    case Y_minus:
#pragma omp parallel for firstprivate(ix, jx, kx, gd, ht1, sf, ie, Sy) schedule(static) reduction(+:va)
      for (int k=1; k<=kx; k++) {
        for (int i=1; i<=ix; i++) {
          size_t m = _F_IDX_S3D(i, 1, k, ix, jx, kx, gd);
          int l = d_bcd[m] & MASK_5;
          REAL_TYPE rho = mtbl[3*l+0];
          REAL_TYPE cp  = mtbl[3*l+1];
          REAL_TYPE t = d_ie0[m] / (rho * cp);
          REAL_TYPE q = ht1 * (sf - t);
          d_qbc[_F_IDX_S4DEX(Y_minus, i, 1, k, NOFACE, ix, jx, kx, gd)] += q;
          va += q * Sy;
          d_ie[_F_IDX_S3D(i, 0, k, ix, jx, kx, gd)] = ie;
        }
      }
      break;
      
    case Y_plus:
#pragma omp parallel for firstprivate(ix, jx, kx, gd, ht1, sf, ie, Sy) schedule(static) reduction(+:va)
      for (int k=1; k<=kx; k++) {
        for (int i=1; i<=ix; i++) {
          size_t m = _F_IDX_S3D(i, jx, k, ix, jx, kx, gd);
          int l = d_bcd[m] & MASK_5;
          REAL_TYPE rho = mtbl[3*l+0];
          REAL_TYPE cp  = mtbl[3*l+1];
          REAL_TYPE t = d_ie0[m] / (rho * cp);
          REAL_TYPE q = ht1 * (t - sf);
          d_qbc[_F_IDX_S4DEX(Y_plus, i, jx, k, NOFACE, ix, jx, kx, gd)] -= q;
          va -= q * Sy;
          d_ie[_F_IDX_S3D(i, jx+1, k, ix, jx, kx, gd)] = ie;
        }
      }
      break;
      
    case Z_minus:
      // ここだけ係数が違うので注意
#pragma omp parallel for firstprivate(ix, jx, kx, gd, ht3, sf, ie, Sz) schedule(static) reduction(+:va)
      for (int j=1; j<=jx; j++) {
        for (int i=1; i<=ix; i++) {
          size_t m = _F_IDX_S3D(i, j, 1, ix, jx, kx, gd);
          int l = d_bcd[m] & MASK_5;
          REAL_TYPE rho = mtbl[3*l+0];
          REAL_TYPE cp  = mtbl[3*l+1];
          REAL_TYPE t = d_ie0[m] / (rho * cp);
          REAL_TYPE q = ht3 * (sf - t);
          d_qbc[_F_IDX_S4DEX(Z_minus, i, j, 1, NOFACE, ix, jx, kx, gd)] += q;
          va += q * Sz;
          d_ie[_F_IDX_S3D(i, j, 0, ix, jx, kx, gd)] = ie;
        }
      }
      break;
      
    case Z_plus:
#pragma omp parallel for firstprivate(ix, jx, kx, gd, ht1, sf, ie, Sz) schedule(static) reduction(+:va)
      for (int j=1; j<=jx; j++) {
        for (int i=1; i<=ix; i++) {
          size_t m = _F_IDX_S3D(i, j, kx, ix, jx, kx, gd);
          int l = d_bcd[m] & MASK_5;
          REAL_TYPE rho = mtbl[3*l+0];
          REAL_TYPE cp  = mtbl[3*l+1];
          REAL_TYPE t = d_ie0[m] / (rho * cp);
          REAL_TYPE q = ht1 * (t - sf);
          d_qbc[_F_IDX_S4DEX(Z_plus, i, j, kx, NOFACE, ix, jx, kx, gd)] -= q;
          va -= q * Sz;
          d_ie[_F_IDX_S3D(i, j, kx+1, ix, jx, kx, gd)] = ie;
        }
      }
      break;
  }
  
  return va;
}



// #################################################################
/**
 @fn REAL_TYPE SetBC3D::setHeatTransferN_SM(REAL_TYPE* qbc, REAL_TYPE* t, int* bx, int n, REAL_TYPE* t0)
 @brief 単媒質の場合の熱伝達境界条件タイプN　共役熱移動の場合
 @retval sum of heat flux (W/m^2)
 @param qbc 境界条件熱流束
 @param t n+1時刻の温度場
 @param bx BCindex
 @param n 境界条件コンポーネントのエントリ番号
 @param t0 n時刻の温度場
 @todo qsumの方向チェック
 
REAL_TYPE SetBC3D::setHeatTransferN_SM(REAL_TYPE* qbc, REAL_TYPE* t, int* bx, int n, REAL_TYPE* t0)
{
  int i,j,k;
  int s, odr, d;
  int m0, mi, mj, mk;
  REAL_TYPE ht, qsum, tmp, c;
  REAL_TYPE qi, qj, qk;

  qsum = 0.0;
  odr= cmp[n].getMatOdr();
  //ht = cmp[n].getCoefHT() / (RefV*mat[odr].P[p_density]*mat[odr].P[p_specific_heat]);
  ht = cmp[n].getCoefHT() / (RefV*rho_0*cp_0);

  for (k=cmp[n].ci.st[2]-1; k<=cmp[n].ci.ed[2]; k++) {
    for (j=cmp[n].ci.st[1]-1; j<=cmp[n].ci.ed[1]; j++) {
      for (i=cmp[n].ci.st[0]-1; i<=cmp[n].ci.ed[0]; i++) {
        m0 = FBUtility::getFindexS3D(size, guide, i  , j  , k  );
        mi = FBUtility::getFindexS3D(size, guide, i+1, j  , k  );
        mj = FBUtility::getFindexS3D(size, guide, i  , j+1, k  );
        mk = FBUtility::getFindexS3D(size, guide, i  , j  , k+1);
        s = bx[m0];
        if ( (s & BIT_MASK_10) == n ) {
          if ( TEST_BIT(s, QFACE_I) ) {
            qi = ht*(t[m0]-t[mi]);
            qbc[FBUtility::getFindexV3DEx(size, guide, 0, i, j, k)] = qi;
            d = BIT_SHIFT(s, DIR_I);
            qsum += qi * ( (0==d) ? 1.0 : -1.0 );
          }
          if ( TEST_BIT(s, QFACE_J) ) {
            qj = ht*(t[m0]-t[mj]);
            qbc[FBUtility::getFindexV3DEx(size, guide, 1, i, j, k)] = qj;
            d = BIT_SHIFT(s, DIR_J);
            qsum += qj * ( (0==d) ? 1.0 : -1.0 );
          }
          if ( TEST_BIT(s, QFACE_K) ) {
            qk = ht*(t[m0]-t[mk]);
            qbc[FBUtility::getFindexV3DEx(size, guide, 2, i, j, k)] = qk;
            d = BIT_SHIFT(s, DIR_K);
            qsum += qk * ( (0==d) ? 1.0 : -1.0 );
          }
        }
      }
    }
  }
  
  if ( numProc > 1 ) 
 {
		tmp = qsum;
		para_mng->Allreduce(&tmp, &qsum, 1, _SUM);
	}
  c = qsum * (RefV*DiffTemp*rho_0*cp_0);
  
  return (c);
}*/



// #################################################################
/**
 * @brief 外部領域の等温熱流束の境界条件処理
 * @retval 無次元の単位時間あたりの熱量（有次元では[W]）
 * @param [in,out] d_qbc  境界条件熱流束
 * @param [in]     face   外部境界面番号
 * @param [in]     d_bcd  BCindex B
 * @param [out]    d_ie   n+1時刻の内部エネルギー
 * @param [in]     d_ie0  n時刻の内部エネルギー
 * @note
   - モニタ量の熱量 va は系に対する流入出量
   - 流入を正にとる
 */
REAL_TYPE SetBC3D::psObcIsoThermal(REAL_TYPE* d_qbc, const int face, const int* d_bcd, REAL_TYPE* d_ie, const REAL_TYPE* d_ie0)
{
  if ( nID[face] >= 0 ) return 0.0;
  
  
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
  REAL_TYPE va = 0.0;
  REAL_TYPE sf = FBUtility::convTempD2ND(obc[face].getTemp(), BaseTemp, DiffTemp); // 保持されている温度はCelsius
  REAL_TYPE px = 2.0 / pitch[0];
  REAL_TYPE py = 2.0 / pitch[1];
  REAL_TYPE pz = 2.0 / pitch[2];
  
  int n = obc[face].getGuideMedium();
  REAL_TYPE rho_1 = mat[n].P[p_density] / rho_0;
  REAL_TYPE cp_1  = mat[n].P[p_specific_heat] / cp_0;
  REAL_TYPE ie  = rho_1 * cp_1 * sf;
  
  REAL_TYPE Sx = S_x;
  REAL_TYPE Sy = S_y;
  REAL_TYPE Sz = S_z;
  
  switch (face)
  {
    case X_minus:
#pragma omp parallel for firstprivate(ix, jx, kx, gd, px, sf, ie, Sx) schedule(static) reduction(+:va)
      for (int k=1; k<=kx; k++) {
        for (int j=1; j<=jx; j++) {
          size_t m = _F_IDX_S3D(1, j, k, ix, jx, kx, gd);
          int l = d_bcd[m] & MASK_5;
          REAL_TYPE rho = mtbl[3*l+0];
          REAL_TYPE cp  = mtbl[3*l+1];
          REAL_TYPE lmd = mtbl[3*l+2] * px;
          REAL_TYPE t = d_ie0[m] / (rho * cp);
          REAL_TYPE q = lmd * (sf - t);
          d_qbc[_F_IDX_S4DEX(X_minus, 1, j, k, NOFACE, ix, jx, kx, gd)] += q; // マイナス面の正の値は流入
          va += q * Sx;
          d_ie[_F_IDX_S3D(0, j, k, ix, jx, kx, gd)] = ie; // 等温セルに指定温度を代入
        }
      }
      break;
      
    case X_plus:
#pragma omp parallel for firstprivate(ix, jx, kx, gd, px, sf, ie, Sx) schedule(static) reduction(+:va)
      for (int k=1; k<=kx; k++) {
        for (int j=1; j<=jx; j++) {
          size_t m = _F_IDX_S3D(ix, j, k, ix, jx, kx, gd);
          int l = d_bcd[m] & MASK_5;
          REAL_TYPE rho = mtbl[3*l+0];
          REAL_TYPE cp  = mtbl[3*l+1];
          REAL_TYPE lmd = mtbl[3*l+2] * px;
          REAL_TYPE t = d_ie0[m] / (rho * cp);
          REAL_TYPE q = lmd * (t - sf);
          d_qbc[_F_IDX_S4DEX(X_plus, ix, j, k, NOFACE, ix, jx, kx, gd)] += q; // プラス面の正の値は流出
          va -= q * Sx;
          d_ie[_F_IDX_S3D(ix+1, j, k, ix, jx, kx, gd)] = ie;
        }
      }
      break;
      
    case Y_minus:
#pragma omp parallel for firstprivate(ix, jx, kx, gd, py, sf, ie, Sy) schedule(static) reduction(+:va)
      for (int k=1; k<=kx; k++) {
        for (int i=1; i<=ix; i++) {
          size_t m = _F_IDX_S3D(i, 1, k, ix, jx, kx, gd);
          int l = d_bcd[m] & MASK_5;
          REAL_TYPE rho = mtbl[3*l+0];
          REAL_TYPE cp  = mtbl[3*l+1];
          REAL_TYPE lmd = mtbl[3*l+2] * py;
          REAL_TYPE t = d_ie0[m] / (rho * cp);
          REAL_TYPE q = lmd * (sf - t);
          d_qbc[_F_IDX_S4DEX(Y_minus, i, 1, k, NOFACE, ix, jx, kx, gd)] += q;
          va += q * Sy;
          d_ie[_F_IDX_S3D(i, 0, k, ix, jx, kx, gd)] = ie;
        }
      }
      break;
      
    case Y_plus:
#pragma omp parallel for firstprivate(ix, jx, kx, gd, py, sf, ie, Sy) schedule(static) reduction(+:va)
      for (int k=1; k<=kx; k++) {
        for (int i=1; i<=ix; i++) {
          size_t m = _F_IDX_S3D(i, jx, k, ix, jx, kx, gd);
          int l = d_bcd[m] & MASK_5;
          REAL_TYPE rho = mtbl[3*l+0];
          REAL_TYPE cp  = mtbl[3*l+1];
          REAL_TYPE lmd = mtbl[3*l+2] * py;
          REAL_TYPE t = d_ie0[m] / (rho * cp);
          REAL_TYPE q = lmd * (t - sf);
          d_qbc[_F_IDX_S4DEX(Y_plus, i, jx, k, NOFACE, ix, jx, kx, gd)] += q;
          va -= q * Sy;
          d_ie[_F_IDX_S3D(i, jx+1, k, ix, jx, kx, gd)] = ie;
        }
      }
      break;
      
    case Z_minus:
#pragma omp parallel for firstprivate(ix, jx, kx, gd, pz, sf, ie, Sz) schedule(static) reduction(+:va)
      for (int j=1; j<=jx; j++) {
        for (int i=1; i<=ix; i++) {
          size_t m = _F_IDX_S3D(i, j, 1, ix, jx, kx, gd);
          int l = d_bcd[m] & MASK_5;
          REAL_TYPE rho = mtbl[3*l+0];
          REAL_TYPE cp  = mtbl[3*l+1];
          REAL_TYPE lmd = mtbl[3*l+2] * pz;
          REAL_TYPE t = d_ie0[m] / (rho * cp);
          REAL_TYPE q = lmd * (sf - t);
          d_qbc[_F_IDX_S4DEX(Z_minus, i, j, 1, NOFACE, ix, jx, kx, gd)] += q;
          va += q * Sz;
          d_ie[_F_IDX_S3D(i, j, 0, ix, jx, kx, gd)] = ie;
        }
      }
      break;
      
    case Z_plus:
#pragma omp parallel for firstprivate(ix, jx, kx, gd, pz, sf, ie, Sz) schedule(static) reduction(+:va)
      for (int j=1; j<=jx; j++) {
        for (int i=1; i<=ix; i++) {
          size_t m = _F_IDX_S3D(i, j, kx, ix, jx, kx, gd);
          int l = d_bcd[m] & MASK_5;
          REAL_TYPE rho = mtbl[3*l+0];
          REAL_TYPE cp  = mtbl[3*l+1];
          REAL_TYPE lmd = mtbl[3*l+2] * pz;
          REAL_TYPE t = d_ie0[m] / (rho * cp);
          REAL_TYPE q = lmd * (t - sf);
          d_qbc[_F_IDX_S4DEX(Z_plus, i, j, kx, NOFACE, ix, jx, kx, gd)] += q;
          va -= q * Sz;
          d_ie[_F_IDX_S3D(i, j, kx+1, ix, jx, kx, gd)] = ie;
        }
      }
      break;
  }
  
  return va;
}


// #################################################################
/**
 * @brief 外部領域の速度指定の境界条件処理
 * @retval 無次元の単位時間あたりの熱量（有次元では[W]）
 * @param [in,out] d_ws  内部エネルギー増分
 * @param [in]     d_cdf BCindex C
 * @param [in]     face  外部境界面番号
 * @param [in]     tm    時刻
 * @param [in]     v00   基準速度
 * @note モニタ量の熱量 va は系に対する流入量なので，基準温度に対する熱量
 */
REAL_TYPE SetBC3D::psObcSpecVH(REAL_TYPE* d_ws, const int* d_cdf, const int face, const double tm, const REAL_TYPE* v00)
{
  if ( nID[face] >= 0 ) return 0.0;
  
  
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
  REAL_TYPE va = 0.0;
  REAL_TYPE rx = 1.0/pitch[0];
  REAL_TYPE ry = 1.0/pitch[1];
  REAL_TYPE rz = 1.0/pitch[2];
  REAL_TYPE vec[3];
  REAL_TYPE u_ref = v00[1];
  REAL_TYPE v_ref = v00[2];
  REAL_TYPE w_ref = v00[3];
  REAL_TYPE dummy = extractVelOBC(face, vec, tm, v00);
  
  int n = obc[face].getGuideMedium();
  REAL_TYPE rho = mat[n].P[p_density] / rho_0;
  REAL_TYPE cp  = mat[n].P[p_specific_heat] / cp_0;
  REAL_TYPE tmp = FBUtility::convTempD2ND(obc[face].getTemp(), BaseTemp, DiffTemp); // 基準温度に対する無次元温度
  REAL_TYPE ie  = rho * cp * tmp;
  
  REAL_TYPE f_x = (vec[0] - u_ref) * ie;
  REAL_TYPE f_y = (vec[1] - v_ref) * ie;
  REAL_TYPE f_z = (vec[2] - w_ref) * ie;
  
  REAL_TYPE Sx = S_x;
  REAL_TYPE Sy = S_y;
  REAL_TYPE Sz = S_z;
  
  switch (face)
  {
    case X_minus:
      
#pragma omp parallel for firstprivate(ix, jx, kx, gd, rx, f_x, Sx) schedule(static) reduction(+:va)
      for (int k=1; k<=kx; k++) {
        for (int j=1; j<=jx; j++) {
          size_t m = _F_IDX_S3D(1, j, k, ix, jx, kx, gd);
          int s = d_cdf[m];
          REAL_TYPE ff = f_x * GET_SHIFT_F(s, STATE_BIT); // 固体マスク
          va += ff * Sx; // 系への流入
          d_ws[m] += ff * rx;
        }
      }
      break;
      
    case X_plus:
      
#pragma omp parallel for firstprivate(ix, jx, kx, gd, rx, f_x, Sx) schedule(static) reduction(+:va)
      for (int k=1; k<=kx; k++) {
        for (int j=1; j<=jx; j++) {
          size_t m = _F_IDX_S3D(ix, j, k, ix, jx, kx, gd);
          int s = d_cdf[m];
          REAL_TYPE ff = f_x * GET_SHIFT_F(s, STATE_BIT);
          va += ff * Sx;
          d_ws[m] -= ff * rx;
        }
      }
      break;
      
    case Y_minus:
      
#pragma omp parallel for firstprivate(ix, jx, kx, gd, ry, f_y, Sy) schedule(static) reduction(+:va)
      for (int k=1; k<=kx; k++) {
        for (int i=1; i<=ix; i++) {
          size_t m = _F_IDX_S3D(i, 1, k, ix, jx, kx, gd);
          int s = d_cdf[m];
          REAL_TYPE ff = f_y * GET_SHIFT_F(s, STATE_BIT);
          va += ff * Sy;
          d_ws[m] += ff * ry;
        }
      }
      break;
      
    case Y_plus:
      
#pragma omp parallel for firstprivate(ix, jx, kx, gd, ry, f_y, Sy) schedule(static) reduction(+:va)
      for (int k=1; k<=kx; k++) {
        for (int i=1; i<=ix; i++) {
          size_t m = _F_IDX_S3D(i, jx, k, ix, jx, kx, gd);
          int s = d_cdf[m];
          REAL_TYPE ff = f_y * GET_SHIFT_F(s, STATE_BIT);
          va += ff * Sy;
          d_ws[m] -= ff * ry;
        }
      }
      break;
      
    case Z_minus:
      
#pragma omp parallel for firstprivate(ix, jx, kx, gd, rz, f_z, Sz) schedule(static) reduction(+:va)
      for (int j=1; j<=jx; j++) {
        for (int i=1; i<=ix; i++) {
          size_t m = _F_IDX_S3D(i, j, 1, ix, jx, kx, gd);
          int s = d_cdf[m];
          REAL_TYPE ff = f_z * GET_SHIFT_F(s, STATE_BIT);
          va += ff * Sz;
          d_ws[m] += ff * rz;
        }
      }
      break;
      
    case Z_plus:
      
#pragma omp parallel for firstprivate(ix, jx, kx, gd, rz, f_z, Sz) schedule(static) reduction(+:va)
      for (int j=1; j<=jx; j++) {
        for (int i=1; i<=ix; i++) {
          size_t m = _F_IDX_S3D(i, j, kx, ix, jx, kx, gd);
          int s = d_cdf[m];
          REAL_TYPE ff = f_z * GET_SHIFT_F(s, STATE_BIT);
          va += ff * Sz;
          d_ws[m] -= ff * rz;
        }
      }
      break;
  }

  return va;
}



// #################################################################
// 初期温度を代入
void SetBC3D::setInitialTempCompo(const int n, const int* d_bx, REAL_TYPE* d_ie)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  int odr = n;
  
  REAL_TYPE tmp = FBUtility::convTempD2ND(cmp[n].getInitTemp(), BaseTemp, DiffTemp);
  REAL_TYPE rho = mat[n].P[p_density] / rho_0;
  REAL_TYPE cp  = mat[n].P[p_specific_heat] / cp_0;
  REAL_TYPE ref = rho * cp * tmp;
  
#pragma omp parallel for firstprivate(ix, jx, kx, gd, odr, ref) schedule(static)
  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        if ( DECODE_CMP(d_bx[m]) == odr )
        {
          d_ie[m] = ref;
        }
      }
    }
  }
}


// #################################################################
// 周期境界の場合のインデクスの同期
// @note PeriodicCommS3D()は周期境界方向に通信が発生する、外部境界面に接している全ランクについて、
//       その引数(最後の2つ)の cpm_DirFlag dir と cpm_PMFlag pm は、全てのランクで同じ値をセットしてコールすること

void SetBC3D::setBCIperiodic(int* d_bx, const int* ens)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  

  if ( numProc > 1 ) 
  {
    // X方向
    if ( ens[0] == ON )
    {
      if ( paraMngr->PeriodicCommS3D(d_bx, ix, jx, kx, gd, 1, X_DIR, PLUS2MINUS) != CPM_SUCCESS ) Exit(0);
      if ( paraMngr->PeriodicCommS3D(d_bx, ix, jx, kx, gd, 1, X_DIR, MINUS2PLUS) != CPM_SUCCESS ) Exit(0);
    }
    
    // Y方向
    if ( ens[1] == ON )
    {
      if ( paraMngr->PeriodicCommS3D(d_bx, ix, jx, kx, gd, 1, Y_DIR, PLUS2MINUS) != CPM_SUCCESS ) Exit(0);
      if ( paraMngr->PeriodicCommS3D(d_bx, ix, jx, kx, gd, 1, Y_DIR, MINUS2PLUS) != CPM_SUCCESS ) Exit(0);
    }
    
    // Z方向
    if ( ens[2] == ON )
    {
      if ( paraMngr->PeriodicCommS3D(d_bx, ix, jx, kx, gd, 1, Z_DIR, PLUS2MINUS) != CPM_SUCCESS ) Exit(0);
      if ( paraMngr->PeriodicCommS3D(d_bx, ix, jx, kx, gd, 1, Z_DIR, MINUS2PLUS) != CPM_SUCCESS ) Exit(0);
    }
  } 
  else // 逐次処理
  {
    // X方向
    if ( ens[0] == ON )
    {
#pragma omp parallel for firstprivate(ix, jx, kx, gd) schedule(static)
      for (int k=1; k<=kx; k++) {
        for (int j=1; j<=jx; j++) {
          for (int i=1-gd; i<=0; i++) {
            size_t m0 = _F_IDX_S3D(i,    j, k, ix, jx, kx, gd);
            size_t m1 = _F_IDX_S3D(i+ix, j, k, ix, jx, kx, gd);
            d_bx[m0] = d_bx[m1];
          }
        }
      }
      
#pragma omp parallel for firstprivate(ix, jx, kx, gd) schedule(static)
      for (int k=1; k<=kx; k++) {
        for (int j=1; j<=jx; j++) {
          for (int i=ix+1; i<=ix+gd; i++) {
            size_t m0 = _F_IDX_S3D(i,    j, k, ix, jx, kx, gd);
            size_t m1 = _F_IDX_S3D(i-ix, j, k, ix, jx, kx, gd);
            d_bx[m0] = d_bx[m1];
          }
        }
      }
    }
    
    // Y方向
    if ( ens[1] == ON )
    {
#pragma omp parallel for firstprivate(ix, jx, kx, gd) schedule(static)
      for (int k=1; k<=kx; k++) {
        for (int j=1-gd; j<=0; j++) {
          for (int i=1; i<=ix; i++) {
            size_t m0 = _F_IDX_S3D(i, j,    k, ix, jx, kx, gd);
            size_t m1 = _F_IDX_S3D(i, j+jx, k, ix, jx, kx, gd);
            d_bx[m0] = d_bx[m1];
          }
        }
      }
      
#pragma omp parallel for firstprivate(ix, jx, kx, gd) schedule(static)
      for (int k=1; k<=kx; k++) {
        for (int j=jx+1; j<=jx+gd; j++) {
          for (int i=1; i<=ix; i++) {
            size_t m0 = _F_IDX_S3D(i, j,    k, ix, jx, kx, gd);
            size_t m1 = _F_IDX_S3D(i, j-jx, k, ix, jx, kx, gd);
            d_bx[m0] = d_bx[m1];
          }
        }
      }
    }
    
    // Z方向
    if ( ens[2] == ON )
    {
#pragma omp parallel for firstprivate(ix, jx, kx, gd) schedule(static)
      for (int k=1-gd; k<=0; k++) {
        for (int j=1; j<=jx; j++) {
          for (int i=1; i<=ix; i++) {
            size_t m0 = _F_IDX_S3D(i, j, k,    ix, jx, kx, gd);
            size_t m1 = _F_IDX_S3D(i, j, k+kx, ix, jx, kx, gd);
            d_bx[m0] = d_bx[m1];
          }
        }
      }
      
#pragma omp parallel for firstprivate(ix, jx, kx, gd) schedule(static)
      for (int k=kx+1; k<=kx+gd; k++) {
        for (int j=1; j<=jx; j++) {
          for (int i=1; i<=ix; i++) {
            size_t m0 = _F_IDX_S3D(i, j, k,    ix, jx, kx, gd);
            size_t m1 = _F_IDX_S3D(i, j, k-kx, ix, jx, kx, gd);
            d_bx[m0] = d_bx[m1];
          }
        }
      }
    }
  }
}


// #################################################################
/**
 @brief 輻射境界条件
 @param qbc 境界条件熱流束
 @param bx BCindex
 @param n 境界条件コンポーネントのエントリ番号
 @param t n時刻の温度場
 @note
    - 未使用，形式のみ
 
void SetBC3D::setRadiant(REAL_TYPE* qbc, int* bx, int n, REAL_TYPE* t)
{
  int i,j,k;
  int s;
  int m0, mi, mj, mk;
  REAL_TYPE ht;

  ht = cmp[n].getCoefHT() / (mat[n].P[p_density]*mat[n].P[p_specific_heat]);

  for (k=cmp[n].ci.st[2]-1; k<=cmp[n].ci.ed[2]; k++) {
    for (j=cmp[n].ci.st[1]-1; j<=cmp[n].ci.ed[1]; j++) {
      for (i=cmp[n].ci.st[0]-1; i<=cmp[n].ci.ed[0]; i++) {
        m0 = FBUtility::getFindexS3D(size, guide, i  , j  , k  );
        mi = FBUtility::getFindexS3D(size, guide, i+1, j  , k  );
        mj = FBUtility::getFindexS3D(size, guide, i  , j+1, k  );
        mk = FBUtility::getFindexS3D(size, guide, i  , j  , k+1);
        s = bx[m0];
        if ( (s & BIT_MASK_10) == n ) {
          if ( TEST_BIT(s, QFACE_I) ) qbc[FBUtility::getFindexV3D(size, guide, 0, i, j, k)] = ht*(t[m0]-t[mi]);
          if ( TEST_BIT(s, QFACE_J) ) qbc[FBUtility::getFindexV3D(size, guide, 1, i, j, k)] = ht*(t[m0]-t[mj]);
          if ( TEST_BIT(s, QFACE_K) ) qbc[FBUtility::getFindexV3D(size, guide, 2, i, j, k)] = ht*(t[m0]-t[mk]);
        }
      }
    }
  }
}*/


// #################################################################
/**
 @brief CN法のときの速度境界条件による粘性項の修正
 @param[out] vc pseudo vector
 @param wv convective term and explicit viscous term of pseudo vector
 @param cf coefficient
 @param bx BCindex for V
 @param wk wroking array for jacobi
 @param tm
 @param dt
 @param omg
 @param C
 @param[out] res residual
 @param LS linear solver identifier
 @param v00
 @param[out] flop
 
void SetBC3D::mod_Vis_CN(REAL_TYPE* vc, REAL_TYPE* wv, REAL_TYPE cf, int* bx, REAL_TYPE* wk, REAL_TYPE tm, REAL_TYPE dt, REAL_TYPE omg, Control* C, REAL_TYPE* res, int LS, REAL_TYPE* v00, double& flop)
{
  int st[3], ed[3];
  REAL_TYPE vel, vec[3], a, b;
  int typ;
  const REAL_TYPE c_pai = (REAL_TYPE)(2.0*asin(1.0));
  
  for (int n=1; n<=NoCompo; n++) {
    if ( cmp[n].isVBC() ) {
      
      cmp[n].getBbox(st, ed);
      
      typ = cmp[n].getType();
      if ( typ==SPEC_VEL ) {
        if ( cmp[n].get_sw_V_profile() == CompoList::vel_constant ) { // Constant
          vel   = cmp[n].ca[CompoList::bias]/RefV * v00[0];
          vec[0]= cmp[n].nv[0]*vel;
          vec[1]= cmp[n].nv[1]*vel;
          vec[2]= cmp[n].nv[2]*vel;
        }
        else { // Harmonic oscillation
          a = cmp[n].ca[CompoList::amplitude]/RefV; // non-dimensional velocity amplitude
          b = 2.0*c_pai*cmp[n].ca[CompoList::frequency]* RefL/RefV * tm + cmp[n].ca[CompoList::initphase];
          vel = ( a*sin(b) + cmp[n].ca[CompoList::bias]/RefV ) * v00[0];
          vec[0]= cmp[n].nv[0]*vel;
          vec[1]= cmp[n].nv[1]*vel;
          vec[2]= cmp[n].nv[2]*vel;
        }
      }
      else if ( typ==OUTFLOW ) {
        vel = cmp[n].val[var_Velocity];
        vec[0]=vel;
        vec[1]=vel;
        vec[2]=vel;
      }
      
      switch (LS) {
        case JACOBI:
          vis_cn_mod_jcb_(vc, size, &gd, st, ed, &dh, &dt, v00, &rei, &omg, wv, bx, wk, &cf, res, vec, flop);
          break;
          
        case SOR:
          vis_cn_mod_sor_(vc, size, &gd, st, ed, &dh, &dt, v00, &rei, &omg, wv, bx, &cf, res, vec, flop);
          break;
          
        case SOR2SMA:
          break;
          
        default:
          printf("\tInvalid Linear Solver for Velocity_CN\n");
          Exit(0);
      }

    }
  }
}*/


// #################################################################
// 対流項計算時の流束型の境界条件処理
void SetBC3D::TBCconvection(REAL_TYPE* d_ws, const int* d_cdf, const REAL_TYPE* d_vf, const REAL_TYPE* d_ie0, const double tm, Control* C, const REAL_TYPE* v00)
{
  REAL_TYPE vec[3];
  REAL_TYPE dummy;
  
  // 外部
  for (int face=0; face<NOFACE; face++)
  {
    // 各面の無次元流入出熱量
    REAL_TYPE va = 0.0;
    
    switch ( obc[face].getClass() )
    {
      case OBC_OUTFLOW:
      case OBC_TRC_FREE:
      case OBC_FAR_FIELD:
        va = psObcFree(d_ws, d_cdf, face, d_vf, d_ie0, v00);
        break;
        
      case OBC_SPEC_VEL:
        va = psObcSpecVH(d_ws, d_cdf, face, tm, v00);
        break;
    }
    
    // vaは無次元内部エネルギー（有次元では[W]）
    if ( numProc > 1 )
    {
      REAL_TYPE tmp = va;
      if ( paraMngr->Allreduce(&tmp, &va, 1, MPI_SUM) != CPM_SUCCESS )
      {
        Hostonly_ printf("Allreduce Error\n");
        Exit(0);
      }
    }

    C->H_Dface[face] += va;
  }
  
  
  // 内部
  for (int n=1; n<=NoCompo; n++)
  {
    switch ( cmp[n].getType() )
    {
      case SPEC_VEL:
        dummy = extractVelLBC(n, vec, tm, v00);
        cmp[n].setMonitorValue( psIbcSpecVH(d_ws, d_cdf, n, v00[0], vec) );
        break;
        
      case OUTFLOW:
        cmp[n].setMonitorValue( psIbcOutflow(d_ws, d_cdf, n, d_vf, d_ie0, v00) );
        break;
        
      case SOLIDREV:
        break;
    }
    
  }
  
}



// #################################################################
/**
 * @brief 温度の外部周期境界条件（単純なコピー）
 * @param [in,out] d_ie 内部エネルギー
 * @param [in]     ens  周期境界方向フラグ
 */
void SetBC3D::TobcPeriodicSimple(REAL_TYPE* d_ie, const int* ens)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
  if ( numProc > 1 ) 
  {
    // X方向
    if ( (ens[0] == ON) &&
        ( (obc[X_minus].getPrdcMode() == BoundaryOuter::prdc_Simple) ||
          (obc[X_minus].getPrdcMode() == BoundaryOuter::prdc_Directional) ))
    {
      if ( paraMngr->PeriodicCommS3D(d_ie, ix, jx, kx, gd, gd, X_DIR, PLUS2MINUS) != CPM_SUCCESS ) Exit(0);
      if ( paraMngr->PeriodicCommS3D(d_ie, ix, jx, kx, gd, gd, X_DIR, MINUS2PLUS) != CPM_SUCCESS ) Exit(0);
    }
    
    // Y方向
    if ( (ens[1] == ON) &&
        ( (obc[Y_minus].getPrdcMode() == BoundaryOuter::prdc_Simple) ||
          (obc[Y_minus].getPrdcMode() == BoundaryOuter::prdc_Directional) ))
    {
      if ( paraMngr->PeriodicCommS3D(d_ie, ix, jx, kx, gd, gd, Y_DIR, PLUS2MINUS) != CPM_SUCCESS ) Exit(0);
      if ( paraMngr->PeriodicCommS3D(d_ie, ix, jx, kx, gd, gd, Y_DIR, MINUS2PLUS) != CPM_SUCCESS ) Exit(0);
    }
    
    // Z方向
    if ( (ens[2] == ON) &&
        ( (obc[Z_minus].getPrdcMode() == BoundaryOuter::prdc_Simple) ||
          (obc[Z_minus].getPrdcMode() == BoundaryOuter::prdc_Directional) ))
    {
      if ( paraMngr->PeriodicCommS3D(d_ie, ix, jx, kx, gd, gd, Z_DIR, PLUS2MINUS) != CPM_SUCCESS ) Exit(0);
      if ( paraMngr->PeriodicCommS3D(d_ie, ix, jx, kx, gd, gd, Z_DIR, MINUS2PLUS) != CPM_SUCCESS ) Exit(0);
    }
  }
  else  // Serial
  {
    // X方向
    if ( (ens[0] == ON) &&
        ( (obc[X_minus].getPrdcMode() == BoundaryOuter::prdc_Simple) ||
          (obc[X_minus].getPrdcMode() == BoundaryOuter::prdc_Directional) ))
    {
#pragma omp parallel for firstprivate(ix, jx, kx, gd) schedule(static)
      for (int k=1; k<=kx; k++) {
        for (int j=1; j<=jx; j++) {
          size_t m0 = _F_IDX_S3D(0,  j, k, ix, jx, kx, gd);
          size_t m1 = _F_IDX_S3D(ix, j, k, ix, jx, kx, gd);
          d_ie[m0] = d_ie[m1];
        }
      }
      
#pragma omp parallel for firstprivate(ix, jx, kx, gd) schedule(static)
      for (int k=1; k<=kx; k++) {
        for (int j=1; j<=jx; j++) {
          size_t m0 = _F_IDX_S3D(ix+1, j, k, ix, jx, kx, gd);
          size_t m1 = _F_IDX_S3D(1,    j, k, ix, jx, kx, gd);
          d_ie[m0] = d_ie[m1];
        }
      }
    }
    
    // Y方向
    if ( (ens[1] == ON) &&
        ( (obc[Y_minus].getPrdcMode() == BoundaryOuter::prdc_Simple) ||
          (obc[Y_minus].getPrdcMode() == BoundaryOuter::prdc_Directional) ))
    {
#pragma omp parallel for firstprivate(ix, jx, kx, gd) schedule(static)
      for (int k=1; k<=kx; k++) {
        for (int i=1; i<=ix; i++) {
          size_t m0 = _F_IDX_S3D(i, 0,  k, ix, jx, kx, gd);
          size_t m1 = _F_IDX_S3D(i, jx, k, ix, jx, kx, gd);
          d_ie[m0] = d_ie[m1];
        }
      }
      
#pragma omp parallel for firstprivate(ix, jx, kx, gd) schedule(static)
      for (int k=1; k<=kx; k++) {
        for (int i=1; i<=ix; i++) {
          size_t m0 = _F_IDX_S3D(i, jx+1, k, ix, jx, kx, gd);
          size_t m1 = _F_IDX_S3D(i, 1,    k, ix, jx, kx, gd);
          d_ie[m0] = d_ie[m1];
        }
      }
    }
    
    // Z方向
    if ( (ens[2] == ON) &&
        ( (obc[Z_minus].getPrdcMode() == BoundaryOuter::prdc_Simple) ||
          (obc[Z_minus].getPrdcMode() == BoundaryOuter::prdc_Directional) ))
    {
#pragma omp parallel for firstprivate(ix, jx, kx, gd) schedule(static)
      for (int j=1; j<=jx; j++) {
        for (int i=1; i<=ix; i++) {
          size_t m0 = _F_IDX_S3D(i, j, 0,  ix, jx, kx, gd);
          size_t m1 = _F_IDX_S3D(i, j, kx, ix, jx, kx, gd);
          d_ie[m0] = d_ie[m1];
        }
      }
      
#pragma omp parallel for firstprivate(ix, jx, kx, gd) schedule(static)
      for (int j=1; j<=jx; j++) {
        for (int i=1; i<=ix; i++) {
          size_t m0 = _F_IDX_S3D(i, j, kx+1, ix, jx, kx, gd);
          size_t m1 = _F_IDX_S3D(i, j, 1,    ix, jx, kx, gd);
          d_ie[m0] = d_ie[m1];
        }
      }
    }
  }
}


// #################################################################
/**
 * @brief 速度の外部周期境界条件（単純なコピー）
 * @param [in,out] d_v  速度ベクトル
 * @param [in]     ens  周期境界方向フラグ
 */
void SetBC3D::VobcPeriodicSimple(REAL_TYPE* d_v, const int* ens)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
  if ( numProc > 1 )
  {
    // X方向
    if ( (ens[0] == ON) &&
        ( (obc[X_minus].getPrdcMode() == BoundaryOuter::prdc_Simple) ||
          (obc[X_minus].getPrdcMode() == BoundaryOuter::prdc_Directional)) )
    {
      if ( paraMngr->PeriodicCommV3D(d_v, ix, jx, kx, gd, gd, X_DIR, PLUS2MINUS) != CPM_SUCCESS ) Exit(0);
      if ( paraMngr->PeriodicCommV3D(d_v, ix, jx, kx, gd, gd, X_DIR, MINUS2PLUS) != CPM_SUCCESS ) Exit(0);
    }
    
    // Y方向
    if ( (ens[1] == ON) &&
        ( (obc[Y_minus].getPrdcMode() == BoundaryOuter::prdc_Simple) ||
          (obc[Y_minus].getPrdcMode() == BoundaryOuter::prdc_Directional)) )
    {
      if ( paraMngr->PeriodicCommV3D(d_v, ix, jx, kx, gd, gd, Y_DIR, PLUS2MINUS) != CPM_SUCCESS ) Exit(0);
      if ( paraMngr->PeriodicCommV3D(d_v, ix, jx, kx, gd, gd, Y_DIR, MINUS2PLUS) != CPM_SUCCESS ) Exit(0);
    }
    
    // Z方向
    if ( (ens[2] == ON) &&
        ( (obc[Z_minus].getPrdcMode() == BoundaryOuter::prdc_Simple) ||
          (obc[Z_minus].getPrdcMode() == BoundaryOuter::prdc_Directional)) )
    {
      if ( paraMngr->PeriodicCommV3D(d_v, ix, jx, kx, gd, gd, Z_DIR, PLUS2MINUS) != CPM_SUCCESS ) Exit(0);
      if ( paraMngr->PeriodicCommV3D(d_v, ix, jx, kx, gd, gd, Z_DIR, MINUS2PLUS) != CPM_SUCCESS ) Exit(0);
    }
  }
  else // Serial
  {
    // X方向
    if ( (ens[0] == ON) &&
        ( (obc[X_minus].getPrdcMode() == BoundaryOuter::prdc_Simple) ||
          (obc[X_minus].getPrdcMode() == BoundaryOuter::prdc_Directional)) )
    {
#pragma omp parallel for firstprivate(ix, jx, kx, gd) schedule(static)
      for (int k=1; k<=kx; k++) {
        for (int j=1; j<=jx; j++) {
          for (int i=1-gd; i<=0; i++) {
            d_v[_F_IDX_V3D(i, j, k, 0, ix, jx, kx, gd)] = d_v[_F_IDX_V3D(ix+i, j, k, 0, ix, jx, kx, gd)];
            d_v[_F_IDX_V3D(i, j, k, 1, ix, jx, kx, gd)] = d_v[_F_IDX_V3D(ix+i, j, k, 1, ix, jx, kx, gd)];
            d_v[_F_IDX_V3D(i, j, k, 2, ix, jx, kx, gd)] = d_v[_F_IDX_V3D(ix+i, j, k, 2, ix, jx, kx, gd)];
          }
        }
      }
      
#pragma omp parallel for firstprivate(ix, jx, kx, gd) schedule(static)
      for (int k=1; k<=kx; k++) {
        for (int j=1; j<=jx; j++) {
          for (int i=ix+1; i<=ix+gd; i++) {
            d_v[_F_IDX_V3D(i, j, k, 0, ix, jx, kx, gd)] = d_v[_F_IDX_V3D(i-ix, j, k, 0, ix, jx, kx, gd)];
            d_v[_F_IDX_V3D(i, j, k, 1, ix, jx, kx, gd)] = d_v[_F_IDX_V3D(i-ix, j, k, 1, ix, jx, kx, gd)];
            d_v[_F_IDX_V3D(i, j, k, 2, ix, jx, kx, gd)] = d_v[_F_IDX_V3D(i-ix, j, k, 2, ix, jx, kx, gd)];
          }
        }
      }
    }
    
    // Y方向
    if ( (ens[1] == ON) &&
        ( (obc[Y_minus].getPrdcMode() == BoundaryOuter::prdc_Simple) ||
          (obc[Y_minus].getPrdcMode() == BoundaryOuter::prdc_Directional)) )
    {
#pragma omp parallel for firstprivate(ix, jx, kx, gd) schedule(static)
      for (int k=1; k<=kx; k++) {
        for (int i=1; i<=ix; i++) {
          for (int j=1-gd; j<=0; j++) {
            d_v[_F_IDX_V3D(i, j, k, 0, ix, jx, kx, gd)] = d_v[_F_IDX_V3D(i, jx+j, k, 0, ix, jx, kx, gd)];
            d_v[_F_IDX_V3D(i, j, k, 1, ix, jx, kx, gd)] = d_v[_F_IDX_V3D(i, jx+j, k, 1, ix, jx, kx, gd)];
            d_v[_F_IDX_V3D(i, j, k, 2, ix, jx, kx, gd)] = d_v[_F_IDX_V3D(i, jx+j, k, 2, ix, jx, kx, gd)];
          }
        }
      }
      
#pragma omp parallel for firstprivate(ix, jx, kx, gd) schedule(static)
      for (int k=1; k<=kx; k++) {
        for (int i=1; i<=ix; i++) {
          for (int j=jx+1; j<=jx+gd; j++) {
            d_v[_F_IDX_V3D(i, j, k, 0, ix, jx, kx, gd)] = d_v[_F_IDX_V3D(i, j-jx, k, 0, ix, jx, kx, gd)];
            d_v[_F_IDX_V3D(i, j, k, 1, ix, jx, kx, gd)] = d_v[_F_IDX_V3D(i, j-jx, k, 1, ix, jx, kx, gd)];
            d_v[_F_IDX_V3D(i, j, k, 2, ix, jx, kx, gd)] = d_v[_F_IDX_V3D(i, j-jx, k, 2, ix, jx, kx, gd)];
          }
        }
      }
    }
    
    // Z方向
    if ( (ens[2] == ON) &&
        ( (obc[Z_minus].getPrdcMode() == BoundaryOuter::prdc_Simple) ||
          (obc[Z_minus].getPrdcMode() == BoundaryOuter::prdc_Directional)) )
    {
#pragma omp parallel for firstprivate(ix, jx, kx, gd) schedule(static)
      for (int j=1; j<=jx; j++) {
        for (int i=1; i<=ix; i++) {
          for (int k=1-gd; k<=0; k++) {
            d_v[_F_IDX_V3D(i, j, k, 0, ix, jx, kx, gd)] = d_v[_F_IDX_V3D(i, j, kx+k, 0, ix, jx, kx, gd)];
            d_v[_F_IDX_V3D(i, j, k, 1, ix, jx, kx, gd)] = d_v[_F_IDX_V3D(i, j, kx+k, 1, ix, jx, kx, gd)];
            d_v[_F_IDX_V3D(i, j, k, 2, ix, jx, kx, gd)] = d_v[_F_IDX_V3D(i, j, kx+k, 2, ix, jx, kx, gd)];
          }
        }
      }
      
#pragma omp parallel for firstprivate(ix, jx, kx, gd) schedule(static)
      for (int j=1; j<=jx; j++) {
        for (int i=1; i<=ix; i++) {
          for (int k=kx+1; k<=kx+gd; k++) {
            d_v[_F_IDX_V3D(i, j, k, 0, ix, jx, kx, gd)] = d_v[_F_IDX_V3D(i, j, k-kx, 0, ix, jx, kx, gd)];
            d_v[_F_IDX_V3D(i, j, k, 1, ix, jx, kx, gd)] = d_v[_F_IDX_V3D(i, j, k-kx, 1, ix, jx, kx, gd)];
            d_v[_F_IDX_V3D(i, j, k, 2, ix, jx, kx, gd)] = d_v[_F_IDX_V3D(i, j, k-kx, 2, ix, jx, kx, gd)];
          }
        }
      }
    }
  }
  
}


// #################################################################
/**
 @brief 速度の内部周期境界条件（単純なコピー）
 @param d_v 速度ベクトル
 @param st コンポーネント範囲の開始インデクス
 @param ed コンポーネント範囲の終了インデクス
 @param d_bd BCindex B
 @param odr 下流面のidが格納されているコンポーネントエントリ
 @param dir ドライバの設置方向
 */
void SetBC3D::Vibc_Prdc(REAL_TYPE* d_v, int* st, int* ed, int* d_bx, int odr, int dir)
{
  if ( numProc > 1 )
  {
    Hostonly_ printf("Error : 'Vibc_Prdc' method is limited to use for serial execution\n.");
    Exit(0);
  }
    
  int ii, jj, kk;

  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
  
  switch (dir)
  {
    case X_minus:
      if ( nID[dir] < 0 )
      {
        int i = st[0];
        for (int k=st[2]; k<=ed[2]; k++) {
          for (int j=st[1]; j<=ed[1]; j++) {
            for (int l=0; l<3; l++) {
              size_t m1 = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
              if ( DECODE_CMP(d_bx[m1]) == odr )
              {
                for (ii=1-gd; ii<=0; ii++) {
                  size_t m0 = _F_IDX_V3D(ii,   j, k, l, ix, jx, kx, gd);
                  size_t m2 = _F_IDX_V3D(ii+i, j, k, l, ix, jx, kx, gd);
                  d_v[m0] = d_v[m2];
                }
              }
            }
          }
        }
      }
      break;
      
    case X_plus:
      if ( nID[dir] < 0 )
      {
        int i = st[0];
        for (int k=st[2]; k<=ed[2]; k++) {
          for (int j=st[1]; j<=ed[1]; j++) {
            for (int l=0; l<3; l++) {
              size_t m1 = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
              if ( DECODE_CMP(d_bx[m1]) == odr )
              {
                for (ii=ix+1; ii<=ix+gd; ii++) {
                  size_t m0 = _F_IDX_V3D(ii,        j, k, l, ix, jx, kx, gd);
                  size_t m2 = _F_IDX_V3D(i+ii-ix-1, j, k, l, ix, jx, kx, gd);
                  d_v[m0] = d_v[m2];
                }
              }              
            }
          }
        }
      }
      break;
      
    case Y_minus:
      if ( nID[dir] < 0 )
      {
        int j = st[1];
        for (int k=st[2]; k<=ed[2]; k++) {
          for (int i=st[0]; i<=ed[0]; i++) {
            for (int l=0; l<3; l++) {
              size_t m1 = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
              if ( DECODE_CMP(d_bx[m1]) == odr )
              {
                for (jj=1-gd; jj<=0; jj++) {
                  size_t m0 = _F_IDX_V3D(i, jj,   k, l, ix, jx, kx, gd);
                  size_t m2 = _F_IDX_V3D(i, jj+j, k, l, ix, jx, kx, gd);
                  d_v[m0] = d_v[m2];
                }
              }              
            }
          }
        }
      }
      break;
      
    case Y_plus:
      if ( nID[dir] < 0 )
      {
        int j = st[1];
        for (int k=st[2]; k<=ed[2]; k++) {
          for (int i=st[0]; i<=ed[0]; i++) {
            for (int l=0; l<3; l++) {
              size_t m1 = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
              if ( DECODE_CMP(d_bx[m1]) == odr ) {
                for (jj=jx+1; jj<=jx+gd; jj++)
                {
                  size_t m0 = _F_IDX_V3D(i, jj,        k, l, ix, jx, kx, gd);
                  size_t m2 = _F_IDX_V3D(i, j+jj-jx-1, k, l, ix, jx, kx, gd);
                  d_v[m0] = d_v[m2];
                }
              }
            }
          }
        }
      }
      break;
      
    case Z_minus:
      if ( nID[dir] < 0 )
      {
        int k = st[2];
        for (int j=st[1]; j<=ed[1]; j++) {
          for (int i=st[0]; i<=ed[0]; i++) {
            for (int l=0; l<3; l++) {
              size_t m1 = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
              if ( DECODE_CMP(d_bx[m1]) == odr )
              {
                for (kk=1-gd; kk<=0; kk++) {
                  size_t m0 = _F_IDX_V3D(i, j, kk,   l, ix, jx, kx, gd);
                  size_t m2 = _F_IDX_V3D(i, j, kk+k, l, ix, jx, kx, gd);
                  d_v[m0] = d_v[m2];
                }
              }         
            }
          }
        }
      }
      break;
      
    case Z_plus:
      if ( nID[dir] < 0 )
      {
        int k = st[2];
        for (int j=st[1]; j<=ed[1]; j++) {
          for (int i=st[0]; i<=ed[0]; i++) {
            for (int l=0; l<3; l++) {
              size_t m1 = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
              if ( DECODE_CMP(d_bx[m1]) == odr )
              {
                for (kk=kx+1; kk<=kx+gd; kk++) {
                  size_t m0 = _F_IDX_V3D(i, j, kk,        l, ix, jx, kx, gd);
                  size_t m2 = _F_IDX_V3D(i, j, k+kk-kx-1, l, ix, jx, kx, gd);
                  d_v[m0] = d_v[m2];
                }
              }         
            }
          }
        }
      }
      break;
  }
}
