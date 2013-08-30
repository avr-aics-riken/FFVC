//##################################################################################
//
// FFV-C : Frontflow / violet Cartesian
//
// Copyright (c) 2007-2011 VCAD System Research Program, RIKEN.
// All rights reserved.
//
// Copyright (c) 2011-2013 Institute of Industrial Science, The University of Tokyo.
// All rights reserved.
//
// Copyright (c) 2012-2013 Advanced Institute for Computational Science, RIKEN.
// All rights reserved.
//
//##################################################################################

/** 
 * @file   ffv_SetBC.C
 * @brief  FFV BC Class
 * @author kero
 */

#include "ffv_SetBC.h"


// #################################################################
// 温度指定境界条件に必要な温度をセットする
void SetBC3D::assignTemp(REAL_TYPE* d_t, int* d_bh1, const double tm, const Control* C)
{
  REAL_TYPE tc;
  
  // 内部境界条件による修正
  for (int n=1; n<=NoCompo; n++)
  {
    int st[3], ed[3];

    cmp[n].getBbox(st, ed);
    
    if ( cmp[n].getType() == SPEC_VEL_WH )
    {
      tc = FBUtility::convK2ND(cmp[n].getTemp(), BaseTemp, DiffTemp); // difference form BaseTemp
      hbc_drchlt_(d_t, size, &guide, st, ed, d_bh1, &n, &tc);
      break;
    }
  }
  /*
   // 外部境界条件
   int n = OBC_MASK;
   for (int face=0; face<NOFACE; face++) {
   typ = obc[face].getClass();
   
   // 計算領域の最外郭領域でないときに，境界処理をスキップ，次のface面を評価
   if( nID[face] >= 0 ) continue;
   
   if ( typ == OBC_SPEC_VEL ) {
   if ( !clear ) {
   extractVelOBC(face, vec, tm, flop);
   }
   else {
   vec[0] = v00[1];
   vec[1] = v00[2];
   vec[2] = v00[3];
   }
   vbc_drchlt_cc_(v, size, &gd, st, ed, v00, bv, &n, vec);
   }
   }*/
}



// #################################################################
// 速度指定境界条件に必要な参照速度をセットする
void SetBC3D::assignVelocity(REAL_TYPE* d_v, int* d_bv, const double tm, REAL_TYPE* v00, bool clear)
{
  REAL_TYPE vec[3];
  int st[3], ed[3];
  int typ;
  int gd = guide;
  REAL_TYPE dummy;
  
  // 内部境界条件による修正
  for (int n=1; n<=NoCompo; n++)
  {
    typ = cmp[n].getType();
    
    cmp[n].getBbox(st, ed);
    
    switch (typ)
    {
      case SPEC_VEL:
      case SPEC_VEL_WH:
        if ( !clear )
        {
          dummy = extractVelLBC(n, vec, tm, v00);
        }
        else
        {
          vec[0] = vec[1] = vec[2] = 0.0;
        }
        vibc_drchlt_(d_v, size, &gd, st, ed, v00, d_bv, &n, vec);
        break;
    }
  }
  
  // 外部境界条件
  for (int face=0; face<NOFACE; face++)
  {
    typ = obc[face].getClass();
    
    if ( typ == OBC_SPEC_VEL )
    {
      if ( !clear )
      {
        dummy = extractVelOBC(face, vec, tm, v00);
      }
      else
      {
        vec[0] = vec[1] = vec[2] = 0.0;
      }
      vobc_drchlt_(d_v, size, &gd, &face, d_bv, vec, nID);
    }

  }
}


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
            case X_MINUS:
            case X_PLUS:
              c_pos = node_st_i + st[0];
              break;
              
            case Y_MINUS:
            case Y_PLUS:
              c_pos = node_st_j + st[1];
              break;
              
            case Z_MINUS:
            case Z_PLUS:
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
 * @brief コンポーネントの速度境界条件の成分を取り出す
 * @param [in]     n    コンポーネントのインデクス
 * @param [out]    vec  ベクトル成分
 * @param [in]     tm   時刻
 * @param [in]     v00  格子速度
 */
REAL_TYPE SetBC3D::extractVelLBC(const int n, REAL_TYPE* vec, const double tm, const REAL_TYPE* v00)
{
  REAL_TYPE a, b, vel;
  const REAL_TYPE c_pai = (REAL_TYPE)(2.0*asin(1.0));
  
  a = cmp[n].ca[CompoList::amplitude]/RefV; // non-dimensional velocity amplitude
  b = 2.0*c_pai*cmp[n].ca[CompoList::frequency]* RefL/RefV * (REAL_TYPE)tm + cmp[n].ca[CompoList::initphase];
  vel = ( a*sin(b) + cmp[n].ca[CompoList::bias]/RefV ) * v00[0];
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
 */
REAL_TYPE SetBC3D::extractVelOBC(const int n, REAL_TYPE* vec, const double tm, const REAL_TYPE* v00)
{
  REAL_TYPE a, b, vel;
  const REAL_TYPE c_pai = (REAL_TYPE)(2.0*asin(1.0));
  
  a = obc[n].ca[CompoList::amplitude]/RefV; // non-dimensional velocity amplitude
  b = 2.0*c_pai*obc[n].ca[CompoList::frequency]* RefL/RefV * (REAL_TYPE)tm + obc[n].ca[CompoList::initphase];
  vel = ( a*sin(b) + obc[n].ca[CompoList::bias]/RefV ) * v00[0];
  vec[0] = obc[n].nv[0] * vel;
  vec[1] = obc[n].nv[1] * vel;
  vec[2] = obc[n].nv[2] * vel;
  
  return vel;
}


// #################################################################
/**
 * @brief 拡散部分に関する温度の内部境界処理
 * @param [in,out] d_qbc 境界条件熱流束
 * @param [in]     d_bh1 BCindex H1
 * @param [in,out] d_t   n+1時刻の温度場
 * @param [in]     d_t0  n時刻の温度
 */
void SetBC3D::InnerTBCface(REAL_TYPE* d_qbc, const int* d_bh1, REAL_TYPE* d_t, const REAL_TYPE* d_t0)
{
  int H;
  
  for (int n=1; n<=NoCompo; n++)
  {
    switch ( cmp[n].getType() )
    {
      case HEATFLUX:
        cmp[n].setMonHeatflux( psIbcHeatflux(d_qbc, d_bh1, n) );
        break;
        
      case TRANSFER:
        H = cmp[n].getHtype();
        //if      ( H == HT_N)  cmp[n].cmp[n].setMonCalorie( setHeatTransferN_SM(qbc, t, bx, n, t0, flop) );
        if      ( H == HT_S)  cmp[n].setMonCalorie( psIbcTransferS (d_qbc, d_bh1, n, d_t0) );
        else if ( H == HT_SN) cmp[n].setMonCalorie( psIbcTransferSN(d_qbc, d_bh1, n, d_t0) );
        else if ( H == HT_SF) cmp[n].setMonCalorie( psIbcTransferSF(d_qbc, d_bh1, n, d_t0) );
        else if ( H == HT_B)  cmp[n].setMonCalorie( psIbcTransferB (d_qbc, d_bh1, n, d_t0) );
        break;
        
      case ISOTHERMAL:
        cmp[n].setMonCalorie( psIbcIsoThermal(d_qbc, d_bh1, n, d_t0) );
        break;
        
      case RADIANT:
        //setRadiant(qbc, bx, n, t0);
        break;
    }
  }
}


// #################################################################
// セルに対する温度の内部境界をセットする
void SetBC3D::InnerTBCvol(REAL_TYPE* d_t, const int* d_bh2, const REAL_TYPE dt)
{
  for (int n=1; n<=NoCompo; n++)
  {
    switch ( cmp[n].getType() )
    {
      case HEAT_SRC:
        cmp[n].setMonCalorie( psIbcHeatGen(d_t, d_bh2, n, dt) );
        break;
        
      case CNST_TEMP:
        psIbcConstTemp(d_t, d_bh2, n);
        break;
    }
  }
}



// #################################################################
/**
 * @brief 速度ベクトルの内部境界条件処理
 * @param [in,out] v    速度ベクトル
 * @param [in]     bv   BCindex V
 * @param [in]     tm   時刻
 * @param [in]     v00  参照速度
 */
void SetBC3D::InnerVBC(REAL_TYPE* d_v, int* d_bv, const double tm, REAL_TYPE* v00)
{
  int st[3], ed[3];
  REAL_TYPE vec[3];
  int gd = guide;
  REAL_TYPE dummy;
  
  if ( isCDS ) // Cut-Distance
  {
    for (int n=1; n<=NoCompo; n++)
    {
      cmp[n].getBbox(st, ed);
      
      switch ( cmp[n].getType() )
      {
        case SPEC_VEL:
        case SPEC_VEL_WH:
          break;
          
        case OUTFLOW:
          break;
          
        default:
          break;
      }  
    }
  }
  else // Binary
  {
    for (int n=1; n<=NoCompo; n++)
    {
      cmp[n].getBbox(st, ed);
      
      switch ( cmp[n].getType() )
      {
        case SPEC_VEL:
        case SPEC_VEL_WH:
          dummy = extractVelLBC(n, vec, tm, v00);
          vibc_drchlt_(d_v, size, &gd, st, ed, v00, d_bv, &n, vec);
          break;
          
        case OUTFLOW:
          vibc_outflow_(d_v, size, &gd, st, ed, d_bv, &n);
          break;
          
        default:
          break;
      }  
    }
  }
  
}



// #################################################################
/**
 @brief 速度ベクトルの内部周期境界条件処理
 @param d_v 速度ベクトルのデータクラス
 @param d_bd BCindex ID
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
 @param d_bcd BCindex ID
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
    pv = FBUtility::convD2ND_P(cmp[n].ca[0], BasePrs, rho, RefV, Unit_Prs);
    
    if ( cmp[n].getType() == PERIODIC )
    {
      Pibc_Prdc(d_p, st, ed, d_bcd, n, dir, pv);
    }
  }
}


// #################################################################
/**
 * @brief 速度境界条件による速度の発散の修正ほか
 * @param [in,out] dv     \sum{u}
 * @param [in]     bv     BCindex V
 * @param [in]     tm     無次元時刻
 * @param [in]     v00    基準速度
 * @param [in]     avr    平均値計算のテンポラリ値
 * @param [in,out] vf     セルフェイス速度 u^{n+1}
 * @param [in,out] v      セルセンター速度 u^{n+1}
 * @param [in]     C      Controlクラス
 * @param [in]     flop   flop count
 * @note 外部境界面のdiv(u)の修正時に領域境界の流量などのモニタ値を計算し，BoundaryOuterクラスに保持 > 反復後にDomainMonitor()で集約
 *       avr[]のインデクスに注意 (Fortran <-> C)
 */

void SetBC3D::modDivergence(REAL_TYPE* dv, int* bv, double tm, REAL_TYPE* v00, Gemini_R* avr, REAL_TYPE* vf, REAL_TYPE* v, Control* C, double& flop)
{
  REAL_TYPE vec[3], dummy;
  int st[3], ed[3];
  int typ=0;
  int gd = guide;
  double fcount = 0.0;
  
  REAL_TYPE aa[2]; // バッファ
  
  // 内部境界条件による修正
  if ( isCDS ) // Cut-Distance
  {
    for (int n=1; n<=NoCompo; n++)
    {
      typ = cmp[n].getType();
      
      cmp[n].getBbox(st, ed);
      
      switch (typ)
      {
        case OUTFLOW:
          //div_ibc_oflow_vec_(dv, size, &gd, st, ed, v00, &coef, bv, &n, &avr[2*n], &fcount);
          break;
          
        case SPEC_VEL:
        case SPEC_VEL_WH:
          cmp[n].val[var_Velocity] = extractVelLBC(n, vec, tm, v00); // 指定された無次元平均流速
          div_ibc_drchlt_(dv, size, &gd, st, ed, v00, bv, &n, vec, &fcount);
          break;
          
        default:
          break;
      }
    }
  }
  else // Binary
  {
    for (int n=1; n<=NoCompo; n++)
    {
      typ = cmp[n].getType();
      
      cmp[n].getBbox(st, ed);
      
      switch (typ)
      {
        case OUTFLOW:
          div_ibc_oflow_vec_(dv, size, &gd, st, ed, bv, &n, aa, &fcount);
          avr[n].p0 = aa[0]; // 積算速度
          avr[n].p1 = aa[1]; // 積算回数
          if ( aa[1] == 0.0 )
          {
            Hostonly_ printf("\tError : Number of accumulation is zero\n");
            Exit(0);
          }
          break;
          
        case SPEC_VEL:
        case SPEC_VEL_WH:
          cmp[n].val[var_Velocity] = extractVelLBC(n, vec, tm, v00); // 指定された無次元平均流速
          div_ibc_drchlt_(dv, size, &gd, st, ed, v00, bv, &n, vec, &fcount);
          break;
          
        default:
          break;
      }
    }
  }

  
  // 外部境界条件による修正
  for (int face=0; face<NOFACE; face++)
  {
    typ = obc[face].getClass();
    
    REAL_TYPE dd;
    
    switch (typ)
    {
      case OBC_SPEC_VEL:
      case OBC_WALL:
        dummy = extractVelOBC(face, vec, tm, v00);
        vobc_div_drchlt_(dv, size, &gd, &face, bv, vec, &dd, nID, &fcount);
        //vobc_face_drchlt_(vf, size, &gd, &face, bv, vec, &dd, nID);
        obc[face].setDomainMF(dd);
        break;
        
      case OBC_INTRINSIC:
        if ( C->Mode.Example == id_Jet )
        {
          ((IP_Jet*)Ex)->divJetInflow(dv, face, bv, vf, aa, fcount);
          obc[face].setDomainV(aa);
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
 @param bd BCindex ID
 @param cvf コンポーネントの体積率
 @param v00 参照速度
 @param[out] flop
 */
void SetBC3D::mod_Dir_Forcing(REAL_TYPE* d_v, int* d_bd, float* d_cvf, REAL_TYPE* v00, double &flop)
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
 @param [in]     bd   BCindex ID
 @param [in]     cvf  コンポーネントの体積率
 @param [in]     v00  参照速度
 @param [in]     dt   時間積分幅
 @param [in,out] flop 浮動小数点演算数
 */
void SetBC3D::mod_Pvec_Forcing(REAL_TYPE* d_vc, REAL_TYPE* d_v, int* d_bd, float* d_cvf, REAL_TYPE* v00, REAL_TYPE dt, double &flop)
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
void SetBC3D::mod_Psrc_Forcing(REAL_TYPE* s_1, REAL_TYPE* v, int* bd, float* cvf, REAL_TYPE* v00, REAL_TYPE** c_array, double &flop)
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
void SetBC3D::mod_Vdiv_Forcing(REAL_TYPE* v, int* bd, float* cvf, REAL_TYPE* dv, REAL_TYPE dt, REAL_TYPE* v00, Gemini_R* am, REAL_TYPE** c_array, double &flop)
{
  int st[3], ed[3], csz[3];
  REAL_TYPE vec[3];
  REAL_TYPE* w_ptr=NULL;
  REAL_TYPE dh = deltaX;
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
          hex_force_vec_(v, dv, size, &gd, st, ed, bd, cvf, w_ptr, csz, &n, v00, &dt, &dh, vec, &cmp[n].ca[0], aa, &flop);
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
void SetBC3D::modPvecFlux(REAL_TYPE* wv, REAL_TYPE* v, int* bv, const double tm, Control* C, int v_mode, REAL_TYPE* v00, double& flop)
{
  REAL_TYPE vec[3], vel, dummy;
  int st[3], ed[3];
  int typ;
  REAL_TYPE dh = deltaX;
  int gd = guide;
  
  // 内部境界（流束形式）
  if ( isCDS ) // Cut-Distance
  {
    for (int n=1; n<=NoCompo; n++)
    {

      typ = cmp[n].getType();
      cmp[n].getBbox(st, ed);
      
      if ( (typ==SPEC_VEL) || (typ==SPEC_VEL_WH) )
      {
        dummy = extractVelLBC(n, vec, tm, v00);
        pvec_vibc_specv_(wv, size, &gd, st, ed, &dh, v00, &rei, v, bv, &n, vec, &flop);
      }
      else if ( typ==OUTFLOW )
      {
        ;
      }

    }
  }
  else // Binary
  {
    for (int n=1; n<=NoCompo; n++)
    {

      typ = cmp[n].getType();
      cmp[n].getBbox(st, ed);
      
      if ( (typ==SPEC_VEL) || (typ==SPEC_VEL_WH) )
      {
        dummy = extractVelLBC(n, vec, tm, v00);
        pvec_vibc_specv_(wv, size, &gd, st, ed, &dh, v00, &rei, v, bv, &n, vec, &flop);
      }
      else if ( typ==OUTFLOW )
      {
        vec[0] = vec[1] = vec[2] = cmp[n].val[var_Velocity]; // modDivergence()でセルフェイス流出速度がval[var_Velocity]にセット
        pvec_vibc_oflow_(wv, size, &gd, st, ed, &dh, &rei, v, bv, &n, vec, &flop);
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
        vobc_pv_specv_(wv, size, &gd, &face, &dh, &rei, v, bv, vec, nID, &flop);
        break;
        
      case OBC_WALL:
        dummy = extractVelOBC(face, vec, tm, v00);
        vobc_pv_wall_(wv, size, &gd, &face, &dh, &rei, v, vec, nID, &flop);
        break;
        
      case OBC_INTRINSIC:
        if ( C->Mode.Example == id_Jet )
        {
          ((IP_Jet*)Ex)->vobc_pv_JetInflow(wv, face, rei, v, bv, flop);
        }
        break;
    }
  }
  
}


// #################################################################
// 速度境界条件によるPoisosn式のソース項の修正
void SetBC3D::modPsrcVBC(REAL_TYPE* s_0, REAL_TYPE* vc, REAL_TYPE* v0, REAL_TYPE* vf, int* bv, const double tm, REAL_TYPE dt, Control* C, REAL_TYPE* v00, double &flop)
{
  int st[3], ed[3];
  REAL_TYPE vec[3], vel;
  REAL_TYPE dh = deltaX;
  int typ;
  int gd = guide;
  double fcount = 0.0;
  
  // 内部境界条件による修正
  if ( isCDS ) // Cut-Distance
  {
    for (int n=1; n<=NoCompo; n++) {
      typ = cmp[n].getType();
      cmp[n].getBbox(st, ed);


      switch (typ)
      {
        case SPEC_VEL:
        case SPEC_VEL_WH:
        {
          REAL_TYPE dummy = extractVelLBC(n, vec, tm, v00);
          div_ibc_drchlt_(s_0, size, &gd, st, ed, v00, bv, &n, vec, &fcount);
          break;
        }
          
        case OUTFLOW:
          break;
          
        default:
          break;
      }

    }
  }
  else // Binary
  {
    for (int n=1; n<=NoCompo; n++) {
      typ = cmp[n].getType();
      cmp[n].getBbox(st, ed);
      
      switch (typ)
      {
        case SPEC_VEL:
        case SPEC_VEL_WH:
        {
          REAL_TYPE dummy = extractVelLBC(n, vec, tm, v00);
          div_ibc_drchlt_(s_0, size, &gd, st, ed, v00, bv, &n, vec, &fcount);
          break;
        }
          
        case OUTFLOW:
          vel = cmp[n].val[var_Velocity] * dt / dh; // modDivergence()でval[var_Velocity]にセット
          div_ibc_oflow_pvec_(s_0, size, &gd, st, ed, v00, &vel, bv, &n, v0, vf, &fcount);
          break;
          
        default:
          break;
      }
    }
  }
  
  
  // 外部境界条件による修正
  for (int face=0; face<NOFACE; face++)
  {
    typ = obc[face].getClass();
    
    REAL_TYPE dd; // dummy
    
    switch ( typ )
    {
      case OBC_SPEC_VEL:
      case OBC_WALL:
      {
        REAL_TYPE dummy = extractVelOBC(face, vec, tm, v00);
        vobc_div_drchlt_(s_0, size, &gd, &face, bv, vec, &dd, nID, &fcount);
        break;
      }
        
      case OBC_INTRINSIC:
        if ( C->Mode.Example == id_Jet )
        {
          ((IP_Jet*)Ex)->divJetInflow(s_0, face, bv, vf, vec, flop);
        }
        break;
        
        // 他は境界値を与え，通常スキームで計算するので不要
    }

  }
  
  flop += fcount;
}


// #################################################################
/**
 @brief Euler陽解法のときの速度境界条件による粘性項の修正
 @param[out] vc 疑似ベクトル
 @param v0 コロケートの速度ベクトル (n step)
 @param cf 係数
 @param bx BCindex for V
 @param tm 無次元時刻
 @param dt 時間積分幅
 @param v00
 @param[out] flop
 @todo symmetricのときの修正
 */
void SetBC3D::mod_Vis_EE(REAL_TYPE* d_vc, REAL_TYPE* d_v0, REAL_TYPE cf, int* d_bx, const double tm, REAL_TYPE dt, REAL_TYPE* v00, double& flop)
{
  int st[3], ed[3];
  REAL_TYPE vec[3];
  REAL_TYPE dh = deltaX;
  int typ;
  int gd = guide;
  
  for (int n=1; n<=NoCompo; n++) {
    typ = cmp[n].getType();
    
    if ( (typ==SPEC_VEL) || (typ==SPEC_VEL_WH) )
    {
      REAL_TYPE dummy = extractVelLBC(n, vec, tm, v00);
    }
    else if ( typ==OUTFLOW )
    {
      vec[0] = vec[1] = vec[2] = cmp[n].val[var_Velocity];
    }
    
    cmp[n].getBbox(st, ed);
    vis_ee_vbc_(d_vc, size, &gd, st, ed, &dh, &dt, v00, &rei, d_v0, d_bx, &n, &cf, vec, &flop);
  }
}


// #################################################################
// 圧力の外部境界条件
void SetBC3D::OuterPBC(REAL_TYPE* d_p)
{
  int uod, F;
  REAL_TYPE pv=0.0;
  int gd = guide;
  
  for (int face=0; face<NOFACE; face++)
  {
    F = obc[face].getClass();
    
    // 周期境界条件
    if ( F == OBC_PERIODIC )
    {
      pv = FBUtility::convD2ND_P(obc[face].p, BasePrs, rho, RefV, Unit_Prs);

      switch ( obc[face].getPrdcMode() )
      {
        case BoundaryOuter::prdc_Simple:
          PobcPeriodicSimple(d_p, face);
          break;
          
        case BoundaryOuter::prdc_Directional:
          uod = obc[face].get_FaceMode();
          PobcPeriodicDirectional(d_p, face, pv, uod);
          break;
          
        case BoundaryOuter::prdc_Driver:
          // nothing
          break;
      }
    }
    else // 周期境界条件以外の処理
    {
      if ( F == OBC_TRC_FREE )
      {
        pobc_drchlt_ (d_p, size, &gd, &face, &pv, nID);
      }

    }
  }
}


// #################################################################
// 対流項計算時の流束型の境界条件処理
void SetBC3D::OuterTBCconvection(REAL_TYPE* d_ws, const int* d_bh1, const REAL_TYPE* d_vf, const REAL_TYPE* d_t0, const double tm, Control* C, const REAL_TYPE* v00)
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
        va = psObcFree(d_ws, d_bh1, face, d_vf, d_t0, v00);
        break;
        
      case OBC_SPEC_VEL:
        va = psObcSpecVH(d_ws, d_bh1, face, tm, v00);
        break;
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
    
    // 対象BCのみ
    switch ( obc[face].getClass() )
    {
      case OBC_OUTFLOW:
      case OBC_TRC_FREE:
      case OBC_FAR_FIELD:
      case OBC_SPEC_VEL:
        C->H_Dface[face] = va;
        break;
    }
    
  }
  
  
  // 内部
  for (int n=1; n<=NoCompo; n++)
  {
    switch ( cmp[n].getType() )
    {
      case SPEC_VEL_WH:
        dummy = extractVelLBC(n, vec, tm, v00);
        cmp[n].setMonCalorie( psIbcSpecVH(d_ws, d_bh1, n, v00[0], vec) );
        break;
        
      case OUTFLOW:
        cmp[n].setMonCalorie( psIbcOutflow(d_ws, d_bh1, n, d_vf, d_t0, v00) );
        break;
    }
    
  }
  
}


// #################################################################
// 拡散項計算時の温度の外部部境界処理
void SetBC3D::OuterTBCdiffusion(REAL_TYPE* d_qbc, REAL_TYPE* d_t, const REAL_TYPE* d_t0, Control* C)
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
          if      ( H == HT_S)  va = psObcHeatTransferBS(d_qbc, face, d_t, d_t0);
          else if ( H == HT_SN) va = psObcHeatTransferSN(d_qbc, face, d_t, d_t0);
          else if ( H == HT_SF) va = psObcHeatTransferSF(d_qbc, face, d_t, d_t0);
          else if ( H == HT_B)  va = psObcHeatTransferBS(d_qbc, face, d_t, d_t0);
          break;
          
        case ISOTHERMAL:
          va = psObcIsoThermal(d_qbc, face, d_t, d_t0);
          break;
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
      
      // 対象BCのみ
      if ( obc[face].getClass() == OBC_WALL )
      {
        switch ( obc[face].getHtype() ) // 熱境界条件の種類
        {
          case HEATFLUX:
          case TRANSFER:
          case ISOTHERMAL:
            C->H_Dface[face] = va;
            break;
        }
      }
    }
  }
}


// #################################################################
// 温度の外部周期境界条件
void SetBC3D::OuterTBCperiodic(REAL_TYPE* d_t)
{
  int F=0;
  
  for (int face=0; face<NOFACE; face++)
  {
    F = obc[face].getClass();
    
    // 周期境界条件
    if ( F == OBC_PERIODIC )
    {
      switch ( obc[face].getPrdcMode() )
      {
        case BoundaryOuter::prdc_Simple:
        case BoundaryOuter::prdc_Directional:
          TobcPeriodicSimple(d_t, face);
          break;
          
        case BoundaryOuter::prdc_Driver:
          // nothing
          break;
      }
    }
  }
}



// #################################################################
// 速度の外部境界条件処理（VP反復内で値を指定する境界条件）
void SetBC3D::OuterVBC(REAL_TYPE* d_v, REAL_TYPE* d_vf, int* d_bv, const double tm, Control* C, REAL_TYPE* v00, double& flop)
{
  REAL_TYPE vec[3];
  int gd = guide;
  REAL_TYPE dd=0.0;
  REAL_TYPE vsum=0.0;
  int pm;
  
  for (int face=0; face<NOFACE; face++)
  {
    switch ( obc[face].getClass() )
    {
      case OBC_TRC_FREE:
        vobc_tfree1_(d_vf, size, &guide, &face, nID);
        if ( numProc > 1 )
        {
          if ( paraMngr->BndCommV3D(d_vf, size[0], size[1], size[2], guide, guide) != CPM_SUCCESS ) Exit(0);
        }
        
        vobc_tfree2_(d_v, size, &guide, &face, d_vf, &vsum, nID, &flop);
        if ( numProc > 1 )
        {
          if ( paraMngr->BndCommV3D(d_v, size[0], size[1], size[2], guide, guide) != CPM_SUCCESS ) Exit(0);
        }
        obc[face].setDomainMF(vsum); // DomainMonitor()で集約する
        break;
        
      case OBC_FAR_FIELD:
        vobc_neumann_(d_v, size, &guide, &face, &vsum, nID);
        obc[face].setDomainMF(vsum);
        break;
        
      case OBC_OUTFLOW:
        vobc_get_massflow_(size, &gd, &face, &vsum, d_v, d_bv, nID);
        obc[face].setDomainMF(vsum);
        break;
        
      case OBC_SYMMETRIC:
        vobc_symmetric_(d_v, size, &gd, &face, nID);
        vsum = 0.0;
        obc[face].setDomainMF(vsum);
        break;
        
      case OBC_PERIODIC:
        pm = obc[face].getPrdcMode();
        
        // BoundaryOuter::prdc_Driverに対しては処理不要
        if ( (pm == BoundaryOuter::prdc_Simple) || (pm == BoundaryOuter::prdc_Directional))
        {
          VobcPeriodicSimple(d_v, face); // セルフェイスの値の周期処理は不要
          vobc_get_massflow_(size, &gd, &face, &vsum, d_v, d_bv, nID);
          obc[face].setDomainMF(vsum);
        }
        break;
        
      default:
        break;
    }
  }
  
}


// #################################################################
// 速度の外部境界処理(タイムステップに一度ガイドセルに値を設定する)
void SetBC3D::OuterVBC_GC(REAL_TYPE* d_v, int* d_bv, const double tm, const Control* C, const REAL_TYPE* v00)
{
  REAL_TYPE vec[3];
  REAL_TYPE dummy;
  int gd = guide;

  
  for (int face=0; face<NOFACE; face++)
  {
    switch ( obc[face].getClass() )
    {
      case OBC_SPEC_VEL:
        dummy = extractVelOBC(face, vec, tm, v00);
        vobc_drchlt_(d_v, size, &gd, &face, d_bv, vec, nID);
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
}


// #################################################################
/**
 * @brief 疑似速度の外部境界条件処理
 * @param [out]    d_vc   疑似速度ベクトル v^*
 * @param [in]     d_bv   BCindex V
 * @param [in]     C      Control class
 * @param [in,out] flop   浮動小数点演算数
 */
void SetBC3D::OuterVBCpseudo(REAL_TYPE* d_vc, int* d_bv, Control* C, double& flop)
{
  REAL_TYPE dd=0.0;
  int gd = guide;
  int pm;
  
  for (int face=0; face<NOFACE; face++)
  {
    switch ( obc[face].getClass() )
    {
      case OBC_OUTFLOW:
      case OBC_FAR_FIELD:
      case OBC_TRC_FREE:
        vobc_neumann_(d_vc, size, &gd, &face, &dd, nID);
        break;
        
      case OBC_SYMMETRIC:
        vobc_symmetric_(d_vc, size, &gd, &face, nID);
        break;
        
      case OBC_PERIODIC:
        pm = obc[face].getPrdcMode();
        
        // BoundaryOuter::prdc_Driverに対しては処理不要
        if ( (pm == BoundaryOuter::prdc_Simple) || (pm == BoundaryOuter::prdc_Directional))
        {
          VobcPeriodicSimple(d_vc, face); // セルフェイスの値の周期処理は不要
        }
        break;
        
    }
  }
}


// #################################################################
/**
 * @brief 圧力の外部周期境界条件（単純なコピー）
 * @param [in,out] d_p  圧力
 * @param [in]     face 面番号
 * @note 並列時には全ノードがPeriodicCommS3D()を評価すること
 */
void SetBC3D::PobcPeriodicSimple(REAL_TYPE* d_p, const int face)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
  if ( numProc > 1 ) 
  {
    switch (face)
    {
      case X_MINUS:
        if ( paraMngr->PeriodicCommS3D(d_p, ix, jx, kx, gd, 1, X_DIR, PLUS2MINUS) != CPM_SUCCESS ) Exit(0);
        break;
        
      case X_PLUS:
        if ( paraMngr->PeriodicCommS3D(d_p, ix, jx, kx, gd, 1, X_DIR, MINUS2PLUS) != CPM_SUCCESS ) Exit(0);
        break;
        
      case Y_MINUS:
        if ( paraMngr->PeriodicCommS3D(d_p, ix, jx, kx, gd, 1, Y_DIR, PLUS2MINUS) != CPM_SUCCESS ) Exit(0);
        break;
        
      case Y_PLUS:
        if ( paraMngr->PeriodicCommS3D(d_p, ix, jx, kx, gd, 1, Y_DIR, MINUS2PLUS) != CPM_SUCCESS ) Exit(0);
        break;
        
      case Z_MINUS:
        if ( paraMngr->PeriodicCommS3D(d_p, ix, jx, kx, gd, 1, Z_DIR, PLUS2MINUS) != CPM_SUCCESS ) Exit(0);
        break;
        
      case Z_PLUS:
        if ( paraMngr->PeriodicCommS3D(d_p, ix, jx, kx, gd, 1, Z_DIR, MINUS2PLUS) != CPM_SUCCESS ) Exit(0);
        break;
    }
  }
  else // Serial
  {
    if ( nID[face] >= 0 ) return;
    
    switch (face) 
    {
      case X_MINUS:
#pragma omp parallel for firstprivate(ix, jx, kx, gd) schedule(static)
        for (int k=1; k<=kx; k++) {
          for (int j=1; j<=jx; j++) {
            size_t m0 = _F_IDX_S3D(0,  j, k, ix, jx, kx, gd);
            size_t m1 = _F_IDX_S3D(ix, j, k, ix, jx, kx, gd);
            d_p[m0] = d_p[m1];
          }
        }
        break;
        
      case X_PLUS:
#pragma omp parallel for firstprivate(ix, jx, kx, gd) schedule(static)
        for (int k=1; k<=kx; k++) {
          for (int j=1; j<=jx; j++) {
            size_t m0 = _F_IDX_S3D(ix+1, j, k, ix, jx, kx, gd);
            size_t m1 = _F_IDX_S3D(1,    j, k, ix, jx, kx, gd);
            d_p[m0] = d_p[m1];
          }
        }
        break;
        
      case Y_MINUS:
#pragma omp parallel for firstprivate(ix, jx, kx, gd) schedule(static)
        for (int k=1; k<=kx; k++) {
          for (int i=1; i<=ix; i++) {
            size_t m0 = _F_IDX_S3D(i, 0,  k, ix, jx, kx, gd);
            size_t m1 = _F_IDX_S3D(i, jx, k, ix, jx, kx, gd);
            d_p[m0] = d_p[m1];
          }
        }
        break;
        
      case Y_PLUS:
#pragma omp parallel for firstprivate(ix, jx, kx, gd) schedule(static)
        for (int k=1; k<=kx; k++) {
          for (int i=1; i<=ix; i++) {
            size_t m0 = _F_IDX_S3D(i, jx+1, k, ix, jx, kx, gd);
            size_t m1 = _F_IDX_S3D(i, 1,    k, ix, jx, kx, gd);
            d_p[m0] = d_p[m1];
          }
        }
        break;
        
      case Z_MINUS:
#pragma omp parallel for firstprivate(ix, jx, kx, gd) schedule(static)
        for (int j=1; j<=jx; j++) {
          for (int i=1; i<=ix; i++) {
            size_t m0 = _F_IDX_S3D(i, j, 0,  ix, jx, kx, gd);
            size_t m1 = _F_IDX_S3D(i, j, kx, ix, jx, kx, gd);
            d_p[m0] = d_p[m1];
          }
        }
        break;
        
      case Z_PLUS:
#pragma omp parallel for firstprivate(ix, jx, kx, gd) schedule(static)
        for (int j=1; j<=jx; j++) {
          for (int i=1; i<=ix; i++) {
            size_t m0 = _F_IDX_S3D(i, j, kx+1, ix, jx, kx, gd);
            size_t m1 = _F_IDX_S3D(i, j, 1,    ix, jx, kx, gd);
            d_p[m0] = d_p[m1];
          }
        }
        break;
    }
  }
}


// #################################################################
/**
 * @brief 圧力の外部周期境界条件（双方向に圧力差を設定）
 * @param [in,out] d_p  圧力のデータクラス
 * @param [in]     face 面番号
 * @param [in]     pv   圧力差
 * @param [in]     uod  上流面 or 下流面
 */
void SetBC3D::PobcPeriodicDirectional(REAL_TYPE* d_p, const int face, REAL_TYPE pv, int uod)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
  // 上流面か下流面かで，圧力差の方向を逆転する
  REAL_TYPE pd = ( uod == BoundaryOuter::prdc_upstream ) ? pv : -pv;
  
  if ( numProc > 1 ) 
  {
    switch (face) 
    {
      case X_MINUS:
        if ( paraMngr->PeriodicCommS3D(d_p, ix, jx, kx, gd, 1, X_DIR, PLUS2MINUS) != CPM_SUCCESS ) Exit(0);
        if ( nID[face] >= 0 ) return;
        
#pragma omp parallel for firstprivate(ix, jx, kx, gd, pd) schedule(static)
        for (int k=1; k<=kx; k++) {
          for (int j=1; j<=jx; j++) {
            size_t m0 = _F_IDX_S3D(0, j, k, ix, jx, kx, gd);
            d_p[m0] += pd;
          }
        }
        break;
        
      case X_PLUS:
        if ( paraMngr->PeriodicCommS3D(d_p, ix, jx, kx, gd, 1, X_DIR, MINUS2PLUS) != CPM_SUCCESS ) Exit(0);
        if ( nID[face] >= 0 ) return;
        
#pragma omp parallel for firstprivate(ix, jx, kx, gd, pd) schedule(static)
        for (int k=1; k<=kx; k++) {
          for (int j=1; j<=jx; j++) {
            size_t m0 = _F_IDX_S3D(ix+1, j, k, ix, jx, kx, gd);
            d_p[m0] += pd;
          }
        }
        break;
        
      case Y_MINUS:
        if ( paraMngr->PeriodicCommS3D(d_p, ix, jx, kx, gd, 1, Y_DIR, PLUS2MINUS) != CPM_SUCCESS ) Exit(0);
        if ( nID[face] >= 0 ) return;
        
#pragma omp parallel for firstprivate(ix, jx, kx, gd, pd) schedule(static)
        for (int k=1; k<=kx; k++) {
          for (int i=1; i<=ix; i++) {
            size_t m0 = _F_IDX_S3D(i, 0, k, ix, jx, kx, gd);
            d_p[m0] += pd;
          }
        }
        break;
        
      case Y_PLUS:
        if ( paraMngr->PeriodicCommS3D(d_p, ix, jx, kx, gd, 1, Y_DIR, MINUS2PLUS) != CPM_SUCCESS ) Exit(0);
        if ( nID[face] >= 0 ) return;
        
#pragma omp parallel for firstprivate(ix, jx, kx, gd, pd) schedule(static)
        for (int k=1; k<=kx; k++) {
          for (int i=1; i<=ix; i++) {
            size_t m0 = _F_IDX_S3D(i, jx+1, k, ix, jx, kx, gd);
            d_p[m0] += pd;
          }
        }
        break;
        
      case Z_MINUS:
        if ( paraMngr->PeriodicCommS3D(d_p, ix, jx, kx, gd, 1, Z_DIR, PLUS2MINUS) != CPM_SUCCESS ) Exit(0);
        if ( nID[face] >= 0 ) return;
        
#pragma omp parallel for firstprivate(ix, jx, kx, gd, pd) schedule(static)
        for (int j=1; j<=jx; j++) {
          for (int i=1; i<=ix; i++) {
            size_t m0 = _F_IDX_S3D(i, j, 0, ix, jx, kx, gd);
            d_p[m0] += pd;
          }
        }
        break;
        
      case Z_PLUS:
        if ( paraMngr->PeriodicCommS3D(d_p, ix, jx, kx, gd, 1, Z_DIR, MINUS2PLUS) != CPM_SUCCESS ) Exit(0);
        if ( nID[face] >= 0 ) return;
        
#pragma omp parallel for firstprivate(ix, jx, kx, gd, pd) schedule(static)
        for (int j=1; j<=jx; j++) {
          for (int i=1; i<=ix; i++) {
            size_t m0 = _F_IDX_S3D(i, j, kx+1, ix, jx, kx, gd);
            d_p[m0] += pd;
          }
        }
        break;
    }
  }
  else // Serial
  {
    if ( nID[face] >= 0 ) return;
    
    switch (face)
    {
      case X_MINUS:
#pragma omp parallel for firstprivate(ix, jx, kx, gd, pd) schedule(static)
        for (int k=1; k<=kx; k++) {
          for (int j=1; j<=jx; j++) {
            size_t m0 = _F_IDX_S3D(0,  j, k, ix, jx, kx, gd);
            size_t m1 = _F_IDX_S3D(ix, j, k, ix, jx, kx, gd);
            d_p[m0] = d_p[m1] + pd;
          }
        }
        break;
        
      case X_PLUS:
#pragma omp parallel for firstprivate(ix, jx, kx, gd, pd) schedule(static)
        for (int k=1; k<=kx; k++) {
          for (int j=1; j<=jx; j++) {
            size_t m0 = _F_IDX_S3D(ix+1, j, k, ix, jx, kx, gd);
            size_t m1 = _F_IDX_S3D(1,    j, k, ix, jx, kx, gd);
            d_p[m0] = d_p[m1] + pd;
          }
        }
        break;
        
      case Y_MINUS:
#pragma omp parallel for firstprivate(ix, jx, kx, gd, pd) schedule(static)
        for (int k=1; k<=kx; k++) {
          for (int i=1; i<=ix; i++) {
            size_t m0 = _F_IDX_S3D(i, 0,  k, ix, jx, kx, gd);
            size_t m1 = _F_IDX_S3D(i, jx, k, ix, jx, kx, gd);
            d_p[m0] = d_p[m1] + pd;
          }
        }
        break;
        
      case Y_PLUS:
#pragma omp parallel for firstprivate(ix, jx, kx, gd, pd) schedule(static)
        for (int k=1; k<=kx; k++) {
          for (int i=1; i<=ix; i++) {
            size_t m0 = _F_IDX_S3D(i, jx+1, k, ix, jx, kx, gd);
            size_t m1 = _F_IDX_S3D(i, 1,    k, ix, jx, kx, gd);
            d_p[m0] = d_p[m1] + pd;
          }
        }
        break;
        
      case Z_MINUS:
#pragma omp parallel for firstprivate(ix, jx, kx, gd, pd) schedule(static)
        for (int j=1; j<=jx; j++) {
          for (int i=1; i<=ix; i++) {
            size_t m0 = _F_IDX_S3D(i, j, 0,  ix, jx, kx, gd);
            size_t m1 = _F_IDX_S3D(i, j, kx, ix, jx, kx, gd);
            d_p[m0] = d_p[m1] + pd;
          }
        }
        break;
        
      case Z_PLUS:
#pragma omp parallel for firstprivate(ix, jx, kx, gd, pd) schedule(static)
        for (int j=1; j<=jx; j++) {
          for (int i=1; i<=ix; i++) {
            size_t m0 = _F_IDX_S3D(i, j, kx+1, ix, jx, kx, gd);
            size_t m1 = _F_IDX_S3D(i, j, 1,    ix, jx, kx, gd);
            d_p[m0] = d_p[m1] + pd;
          }
        }
        break;
    }
  }
}


// #################################################################
/**
 @brief 圧力の内部周期境界条件（一方向の圧力差）
 @param d_p 圧力のデータクラス
 @param st コンポーネント範囲の開始インデクス
 @param ed コンポーネント範囲の終了インデクス
 @param d_bcd BCindex ID
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
    case X_MINUS:
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
      
    case X_PLUS:
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
      
    case Y_MINUS:
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
      
    case Y_PLUS:
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
      
    case Z_MINUS:
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
      
    case Z_PLUS:
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
  
  vel   = FBUtility::convD2ND_V(cmp[n].D1.Velocity, RefV)  * v00;
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
 * @param [in,out] d_t   温度場
 * @param [in]     d_bh2 BCindex H2
 * @param [in]     n     境界条件コンポーネントのエントリ番号
 */
void SetBC3D::psIbcConstTemp(REAL_TYPE* d_t, const int* d_bh2, const int n)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  int odr = n;
  int st[3], ed[3];
  
  REAL_TYPE tmp = FBUtility::convK2ND(cmp[n].getTemp(), BaseTemp, DiffTemp);
  
  cmp[n].getBbox(st, ed);
  
#pragma omp parallel for firstprivate(ix, jx, kx, gd, st, ed, odr, tmp) schedule(static)
  for (int k=st[2]; k<=ed[2]; k++) {
    for (int j=st[1]; j<=ed[1]; j++) {
      for (int i=st[0]; i<=ed[0]; i++) {
        size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        if ( (d_bh2[m] & MASK_6) == odr ) d_t[m] = tmp;
      }
    }
  }
}



// #################################################################
/**
 * @brief 内部領域の熱流束指定境界条件
 * @retval 熱量(W)
 * @param [in,out] qbc 境界条件熱流束
 * @param [in]     bh1 BCindex H1
 * @param [in]     n   コンポーネントリストのインデクス
 * @note
   - モニタ量の熱量vaは系に対する流入量なので，基準温度に対する熱量
   - 流入量を正にとる
 */
REAL_TYPE SetBC3D::psIbcHeatflux(REAL_TYPE* d_qbc, const int* d_bh1, const int n)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  int odr= n;
  
  REAL_TYPE va=0.0;
  int st[3], ed[3];
  
  REAL_TYPE q = cmp[n].getHeatflux() / (RefV*DiffTemp*rho*cp); // [W/m^2]を無次元化
  
  cmp[n].getBbox(st, ed);
	
#pragma omp parallel for firstprivate(ix, jx, kx, gd, st, ed, odr, q) schedule(static) reduction(+:va)
  for (int k=st[2]; k<=ed[2]; k++) {
    for (int j=st[1]; j<=ed[1]; j++) {
      for (int i=st[0]; i<=ed[0]; i++) {
        
        size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        int s = d_bh1[m];
        
        if ( GET_FACE_BC(s, BC_FACE_W) == odr )
        {
          d_qbc[_F_IDX_V3DEX(0, i-1, j, k, ix, jx, kx, gd)] += q;
          va += q; // マイナス面の正の値は流入
        }
        
        if ( GET_FACE_BC(s, BC_FACE_E) == odr )
        {
          d_qbc[_F_IDX_V3DEX(0, i, j, k, ix, jx, kx, gd)] -= q;
          va -= q; // プラス面の正の値は流出
        }
        
        if ( GET_FACE_BC(s, BC_FACE_S) == odr )
        {
          d_qbc[_F_IDX_V3DEX(1, i, j-1, k, ix, jx, kx, gd)] += q;
          va += q;
        }
        
        if ( GET_FACE_BC(s, BC_FACE_N) == odr )
        {
          d_qbc[_F_IDX_V3DEX(1, i, j, k, ix, jx, kx, gd)] -= q;
          va -= q;
        }
        
        if ( GET_FACE_BC(s, BC_FACE_B) == odr )
        {
          d_qbc[_F_IDX_V3DEX(2, i, j, k-1, ix, jx, kx, gd)] += q;
          va += q;
        }
        
        if ( GET_FACE_BC(s, BC_FACE_T) == odr )
        {
          d_qbc[_F_IDX_V3DEX(2, i, j, k, ix, jx, kx, gd)] -= q;
          va -= q;
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
  
  return va; // 単位面積あたりの無次元熱流束の和
}


// #################################################################
/**
 * @brief 発熱境界条件
 * @param [in,out] d_t   温度場
 * @param [in]     d_bh2 BCindex H2
 * @param [in]     n     コンポーネントのエントリ番号
 * @param [in]     dt    時間積分幅
 * @note 発熱密度はControl::setParameters()で計算 D2に発熱密度が保存されている
 */
REAL_TYPE SetBC3D::psIbcHeatGen(REAL_TYPE* d_t, const int* d_bh2, const int n, const REAL_TYPE dt)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  int odr= n;
  int st[3], ed[3];

  REAL_TYPE c = 0.0;

  REAL_TYPE hs = FBUtility::convD2ND_Hsrc(cmp[n].getHeatDensity(), RefV, RefL, DiffTemp, mat[n].P[p_density], mat[n].P[p_specific_heat]);
  //REAL_TYPE hs = FBUtility::convD2ND_Hsrc(cmp[n].getHeatDensity(), RefV, RefL, DiffTemp, rho, cp);
  REAL_TYPE dhs = dt * hs;

  cmp[n].getBbox(st, ed);
	
#pragma omp parallel for firstprivate(ix, jx, kx, gd, st, ed, odr, dhs) schedule(static) reduction(+:c)
  for (int k=st[2]; k<=ed[2]; k++) {
    for (int j=st[1]; j<=ed[1]; j++) {
      for (int i=st[0]; i<=ed[0]; i++) {
        
        size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        if ( (d_bh2[m] & MASK_6) == odr )
        {
          d_t[m] += dhs;
          c++; 
        }
      }
    }
  }
  
  return c*hs; // 無次元の総発熱量(単位体積あたり)
}


// #################################################################
/**
 * @brief 等温境界条件
 * @retval 無次元熱流束の和
 * @param [in,out] d_qbc  境界条件熱流束
 * @param [in]     d_bh1  BCindex H1
 * @param [in]     n      境界条件コンポーネントのエントリ番号
 * @param [in]     d_t0   n時刻の温度場
 * @note 熱流束は加算（他の条件と合成）
 */
REAL_TYPE SetBC3D::psIbcIsoThermal(REAL_TYPE* d_qbc, const int* d_bh1, const int n, const REAL_TYPE* d_t0)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  int odr= n;
  
  int st[3], ed[3];
  
  REAL_TYPE va = 0.0;
  REAL_TYPE dh = deltaX;
  REAL_TYPE sf = FBUtility::convK2ND(cmp[n].getTemp(), BaseTemp, DiffTemp); // 保持されている温度はKelvin
  //REAL_TYPE pp = (2.0*mat[n].P[p_thermal_conductivity]) / (dh*RefV*RefL*mat[n].P[p_density]*mat[n].P[p_specific_heat]);
  REAL_TYPE pp = (2.0*lambda) / (dh*RefV*rho*cp); // RefMediumで指定される物性値
  // 等温境界は熱流動でも固体伝熱でも使うので，

  cmp[n].getBbox(st, ed);
  
#pragma omp parallel for firstprivate(ix, jx, kx, gd, st, ed, odr, pp, sf) schedule(static) reduction(+:va)
  for (int k=st[2]; k<=ed[2]; k++) {
    for (int j=st[1]; j<=ed[1]; j++) {
      for (int i=st[0]; i<=ed[0]; i++) {
        
        size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        int s = d_bh1[m];
        REAL_TYPE q;
        
        if ( GET_FACE_BC(s, BC_FACE_W) == odr )
        {
          q = pp * (sf - d_t0[m]);
          d_qbc[_F_IDX_V3DEX(0, i-1, j, k, ix, jx, kx, gd)] += q;
          va += q; // マイナス面の正の値は流入
        }
        
        if ( GET_FACE_BC(s, BC_FACE_E) == odr )
        {
          q = pp * (d_t0[m] - sf);
          d_qbc[_F_IDX_V3DEX(0, i, j, k, ix, jx, kx, gd)] += q;
          va -= q; // プラス面の正の値は流出
        }
        
        if ( GET_FACE_BC(s, BC_FACE_S) == odr )
        {
          q = pp * (sf - d_t0[m]);
          d_qbc[_F_IDX_V3DEX(1, i, j-1, k, ix, jx, kx, gd)] += q;
          va += q;
        }
        
        if ( GET_FACE_BC(s, BC_FACE_N) == odr )
        {
          q = pp * (d_t0[m] - sf);
          d_qbc[_F_IDX_V3DEX(1, i, j, k, ix, jx, kx, gd)] += q;
          va -= q;
        }
        
        if ( GET_FACE_BC(s, BC_FACE_B) == odr )
        {
          q = pp * (sf - d_t0[m]);
          d_qbc[_F_IDX_V3DEX(2, i, j, k-1, ix, jx, kx, gd)] += q;
          va += q;
        }
        
        if ( GET_FACE_BC(s, BC_FACE_T) == odr )
        {
          q = pp * (d_t0[m] - sf);
          d_qbc[_F_IDX_V3DEX(2, i, j, k, ix, jx, kx, gd)] += q;
          va -= q;
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
  
  return va; // 単位面積あたりの無次元熱流束の和
}



// #################################################################
/**
 * @brief 内部領域のOutflowの境界条件処理
 * @retval 熱量
 * @param [in,out] d_ws  温度増分
 * @param [in]     d_bh1 BCindex H1
 * @param [in]     n     コンポーネントリストのインデクス
 * @param [in]     d_v   速度
 * @param [in]     d_t   温度
 * @param [in]     v00   参照速度
 * @note 
   - モニタ量の熱量va(W)は系に対する流入量なので，基準温度に対する熱量
   - @todo 流出速度はモニター値を利用
 */
REAL_TYPE SetBC3D::psIbcOutflow (REAL_TYPE* d_ws, const int* d_bh1, const int n, const REAL_TYPE* d_v, const REAL_TYPE* d_t, const REAL_TYPE* v00)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  int odr= n;
  
  int st[3], ed[3];
  REAL_TYPE va=0.0;
  REAL_TYPE dh1 = 1.0/deltaX;;
  REAL_TYPE u_ref = v00[1];
  REAL_TYPE v_ref = v00[2];
  REAL_TYPE w_ref = v00[3];
  
  cmp[n].getBbox(st, ed);
	
#pragma omp parallel for firstprivate(ix, jx, kx, gd, st, ed, odr, dh1, u_ref, v_ref, w_ref) \
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
        int s = d_bh1[m];
        REAL_TYPE t_p = d_t[m];
        
        if ( GET_FACE_BC(s, BC_FACE_W) == odr )
        {
          size_t m_0 = _F_IDX_V3D(i,   j, k, 0, ix, jx, kx, gd);
          size_t m_w = _F_IDX_V3D(i-1, j, k, 0, ix, jx, kx, gd);
          c = 0.5*(d_v[m_w]+d_v[m_0]) - u_ref;
          if ( c>0.0 ) c=0.0;
          f_w = c*t_p;
          va += f_w;
        }
        
        if ( GET_FACE_BC(s, BC_FACE_E) == odr )
        {
          size_t m_0 = _F_IDX_V3D(i,   j, k, 0, ix, jx, kx, gd);
          size_t m_e = _F_IDX_V3D(i+1, j, k, 0, ix, jx, kx, gd);
          c = 0.5*(d_v[m_e]+d_v[m_0]) - u_ref;
          if ( c<0.0 ) c=0.0;
          f_e = c*t_p;
          va += f_e;
        }
        
        if ( GET_FACE_BC(s, BC_FACE_S) == odr )
        {
          size_t m_0 = _F_IDX_V3D(i, j,   k, 1, ix, jx, kx, gd);
          size_t m_s = _F_IDX_V3D(i, j-1, k, 1, ix, jx, kx, gd);
          c = 0.5*(d_v[m_s]+d_v[m_0]) - v_ref;
          if ( c>0.0 ) c=0.0;
          f_s = c*t_p;
          va += f_s;
        }
        
        if ( GET_FACE_BC(s, BC_FACE_N) == odr )
        {
          size_t m_0 = _F_IDX_V3D(i, j,   k, 1, ix, jx, kx, gd);
          size_t m_n = _F_IDX_V3D(i, j+1, k, 1, ix, jx, kx, gd);
          c = 0.5*(d_v[m_n]+d_v[m_0]) - v_ref;
          if ( c<0.0 ) c=0.0;
          f_n = c*t_p;
          va += f_n;
        }
        
        if ( GET_FACE_BC(s, BC_FACE_B) == odr )
        {
          size_t m_0 = _F_IDX_V3D(i, j, k,   2, ix, jx, kx, gd);
          size_t m_b = _F_IDX_V3D(i, j, k-1, 2, ix, jx, kx, gd);
          c = 0.5*(d_v[m_b]+d_v[m_0]) - w_ref;
          if ( c>0.0 ) c=0.0;
          f_b = c*t_p;
          va += f_b;
        }
        
        if ( GET_FACE_BC(s, BC_FACE_T) == odr )
        {
          size_t m_0 = _F_IDX_V3D(i, j, k,   2, ix, jx, kx, gd);
          size_t m_t = _F_IDX_V3D(i, j, k+1, 2, ix, jx, kx, gd);
          c = 0.5*(d_v[m_t]+d_v[m_0]) - w_ref;
          if ( c<0.0 ) c=0.0;
          f_t = c*t_p;
          va += f_t;
        }
        
        REAL_TYPE ff = f_e - f_w + f_n - f_s + f_t - f_b;
        d_ws[m] -= ff*dh1; // Outflow指定セルは，必ずFセルなのでマスク不要
      }
    }
  }
  
  
  if ( numProc > 1 )
  {
		REAL_TYPE tmp = va;
		if ( paraMngr->Allreduce(&tmp, &va, 1, MPI_SUM) != CPM_SUCCESS ) Exit(0);
	}
  
  return va; // 単位面積あたりの無次元熱流束の和
}



// #################################################################
/**
 * @brief 内部領域の速度と温度の指定境界条件
 * @retval 熱量[-]
 * @param [in,out] d_ws   温度増分
 * @param [in]     d_bh1  BCindex H1
 * @param [in]     n      コンポーネントリストのインデクス
 * @param [in]     v00    参照速度
 * @param [in]     vec    指定ベクトル
 * @note モニタ量の熱量va(無次元)は系に対する流入量なので，基準温度に対する熱量
 */
REAL_TYPE SetBC3D::psIbcSpecVH (REAL_TYPE* d_ws, const int* d_bh1, const int n, const REAL_TYPE v00, const REAL_TYPE* vec)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  int odr= n;

  int st[3], ed[3];
  REAL_TYPE va = 0.0;
  
  REAL_TYPE dh = deltaX;
  REAL_TYPE dh1 = 1.0/dh;
  REAL_TYPE tc  = FBUtility::convK2ND(cmp[n].getTemp(), BaseTemp, DiffTemp); // difference form BaseTemp
  REAL_TYPE hu  = vec[0]*tc;
  REAL_TYPE hv  = vec[1]*tc;
  REAL_TYPE hw  = vec[2]*tc;
  
  cmp[n].getBbox(st, ed);

#pragma omp parallel for firstprivate(ix, jx, kx, gd, st, ed, hu, hv, hw, odr, dh1) \
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
        int s = d_bh1[m];
        
        if ( GET_FACE_BC(s, BC_FACE_W) == odr )
        {
          f_w = hu;
          va += hu;
        }
        
        if ( GET_FACE_BC(s, BC_FACE_E) == odr )
        {
          f_e = hu;
          va += hu;
        }
        
        if ( GET_FACE_BC(s, BC_FACE_S) == odr )
        {
          f_s = hv;
          va += hv;
        }
        
        if ( GET_FACE_BC(s, BC_FACE_N) == odr )
        {
          f_n = hv;
          va += hv;
        }
        
        if ( GET_FACE_BC(s, BC_FACE_B) == odr )
        {
          f_b = hw;
          va += hw;
        }
        
        if ( GET_FACE_BC(s, BC_FACE_T) == odr )
        {
          f_t = hw;
          va += hw;
        }
        
        REAL_TYPE ff = f_e - f_w + f_n - f_s + f_t - f_b;
        d_ws[m] -= ff*dh1; // 流入境界は必ずFセルなのでマスク不要
      }
    }
  }
  
  if ( numProc > 1 )
  {
		REAL_TYPE tmp = va;
    if ( paraMngr->Allreduce(&tmp, &va, 1, MPI_SUM) != CPM_SUCCESS ) Exit(0);
	}
  
  return va; // 単位面積あたりの無次元熱流束の和
}



// #################################################################
/**
 * @brief 熱伝達境界条件タイプB 固体のみを解く場合
 * @retval 無次元熱流束の和
 * @param [in,out] qbc  境界条件熱流束
 * @param [in]     bh1  BCindex H1
 * @param [in]     n    境界条件コンポーネントのエントリ番号
 * @param [in]     d_t0 n時刻の温度場
 * @note
   - 熱流束は加算（他の条件と合成）
   - pセルは固体セル
 */
REAL_TYPE SetBC3D::psIbcTransferB (REAL_TYPE* d_qbc, const int* d_bh1, const int n, const REAL_TYPE* d_t0)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  int odr= n;
  
  REAL_TYPE va = 0.0;
  int st[3], ed[3];
  REAL_TYPE q;
  
  REAL_TYPE ht = cmp[n].getCoefHT() / (RefV*rho*cp);                        // RefMediumで指定される物性値
  REAL_TYPE bt = FBUtility::convK2ND(cmp[n].getTemp(), BaseTemp, DiffTemp); // 保持されている温度はKelvin

  cmp[n].getBbox(st, ed);
  
#pragma omp parallel for firstprivate(ix, jx, kx, gd, st, ed, odr, ht, bt) private(q) schedule(static) reduction(+:va)
  for (int k=st[2]; k<=ed[2]; k++) {
    for (int j=st[1]; j<=ed[1]; j++) {
      for (int i=st[0]; i<=ed[0]; i++) {
        
        size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        int s = d_bh1[m];
        
        if ( GET_FACE_BC(s, BC_FACE_W) == odr )
        {
          q = ht * (bt - d_t0[m]);
          d_qbc[_F_IDX_V3DEX(0, i-1, j, k, ix, jx, kx, gd)] += q;
          va += q; // マイナス面の正の値は流入
        }
        
        if ( GET_FACE_BC(s, BC_FACE_E) == odr )
        {
          q = ht * (d_t0[m] - bt);
          d_qbc[_F_IDX_V3DEX(0, i, j, k, ix, jx, kx, gd)] += q;
          va -= q; // プラス面の正の値は流出
        }
        
        if ( GET_FACE_BC(s, BC_FACE_S) == odr )
        {
          q = ht * (bt - d_t0[m]);
          d_qbc[_F_IDX_V3DEX(1, i, j-1, k, ix, jx, kx, gd)] += q;
          va += q;
        }
        
        if ( GET_FACE_BC(s, BC_FACE_N) == odr )
        {
          q = ht * (d_t0[m] - bt);
          d_qbc[_F_IDX_V3DEX(1, i, j, k, ix, jx, kx, gd)] += q;
          va -= q;
        }
        
        if ( GET_FACE_BC(s, BC_FACE_B) == odr )
        {
          q = ht * (bt - d_t0[m]);
          d_qbc[_F_IDX_V3DEX(2, i, j, k-1, ix, jx, kx, gd)] += q;
          va += q;
        }
        
        if ( GET_FACE_BC(s, BC_FACE_T) == odr )
        {
          q = ht * (d_t0[m] - bt);
          d_qbc[_F_IDX_V3DEX(2, i, j, k, ix, jx, kx, gd)] += q;
          va -= q;
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
  
  return va; // 単位面積あたりの無次元熱流束の和
}


// #################################################################
/**
 * @brief 熱伝達境界条件タイプS
 * @retval 無次元熱流束の和
 * @param [in,out] d_qbc 境界条件熱流束
 * @param [in]     d_bh1 BCindex H1
 * @param [in]     n     境界条件コンポーネントのエントリ番号
 * @param [in]     d_t0  n時刻の温度場
 * @note
   - 熱流束は加算（他の条件と合成）
   - pセルは流体セル
 */
REAL_TYPE SetBC3D::psIbcTransferS (REAL_TYPE* d_qbc, const int* d_bh1, const int n, const REAL_TYPE* d_t0)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  int odr= n;
  
  REAL_TYPE va = 0.0;
  REAL_TYPE q;
  int st[3], ed[3];
  
  REAL_TYPE ht = cmp[n].getCoefHT() / (RefV*rho*cp);                        // RefMediumで指定される物性値
  REAL_TYPE sf = FBUtility::convK2ND(cmp[n].getTemp(), BaseTemp, DiffTemp); // 保持されている温度はKelvin
  
  cmp[n].getBbox(st, ed);
  
#pragma omp parallel for firstprivate(ix, jx, kx, gd, st, ed, odr, ht, sf) private(q) schedule(static) reduction(+:va)
  for (int k=st[2]; k<=ed[2]; k++) {
    for (int j=st[1]; j<=ed[1]; j++) {
      for (int i=st[0]; i<=ed[0]; i++) {
        
        size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        int s = d_bh1[m];
        
        if ( GET_FACE_BC(s, BC_FACE_W) == odr )
        {
          q = ht * (sf - d_t0[m]);
          d_qbc[_F_IDX_V3DEX(0, i-1, j, k, ix, jx, kx, gd)] += q;
          va += q; // マイナス面の正の値は流入
        }
        
        if ( GET_FACE_BC(s, BC_FACE_E) == odr )
        {
          q = ht * (d_t0[m] - sf);
          d_qbc[_F_IDX_V3DEX(0, i, j, k, ix, jx, kx, gd)] += q;
          va -= q; // プラス面の正の値は流出
        }
        
        if ( GET_FACE_BC(s, BC_FACE_S) == odr )
        {
          q = ht * (sf - d_t0[m]);
          d_qbc[_F_IDX_V3DEX(1, i, j-1, k, ix, jx, kx, gd)] += q;
          va += q;
        }
        
        if ( GET_FACE_BC(s, BC_FACE_N) == odr )
        {
          q = ht * (d_t0[m] - sf);
          d_qbc[_F_IDX_V3DEX(1, i, j, k, ix, jx, kx, gd)] += q;
          va -= q;
        }
        
        if ( GET_FACE_BC(s, BC_FACE_B) == odr )
        {
          q = ht * (sf - d_t0[m]);
          d_qbc[_F_IDX_V3DEX(2, i, j, k-1, ix, jx, kx, gd)] += q;
          va += q;
        }
        
        if ( GET_FACE_BC(s, BC_FACE_T) == odr )
        {
          q = ht * (d_t0[m] - sf);
          d_qbc[_F_IDX_V3DEX(2, i, j, k, ix, jx, kx, gd)] += q;
          va -= q;
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
  
  return va; // 単位面積あたりの無次元熱流束の和
}


// #################################################################
/**
 * @brief 熱伝達境界条件タイプSF（強制対流)
 * @retval 無次元熱流束の和
 * @param [in] d_qbc 境界条件熱流束
 * @param [in] d_bh1 BCindex H1
 * @param [in] n     境界条件コンポーネントのエントリ番号
 * @param [in] d_t0  n時刻の温度場
 * @note
   - 熱流束は加算（他の条件と合成）
   - pセルは流体セル
 */
REAL_TYPE SetBC3D::psIbcTransferSF (REAL_TYPE* d_qbc, const int* d_bh1, const int n, const REAL_TYPE* d_t0)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  int odr= n;
  
  REAL_TYPE va = 0.0;
  REAL_TYPE q;
  int st[3], ed[3];
  
  REAL_TYPE a1 = cmp[n].ca[CompoList::alpha];
  REAL_TYPE b1 = cmp[n].ca[CompoList::beta];
  REAL_TYPE c1 = cmp[n].ca[CompoList::gamma];
  REAL_TYPE ht = lambda / (RefV*RefL*rho*cp) * a1*pow(Reynolds, b1) * pow(Prandtl, c1);    // RefMediumで指定される物性値
  REAL_TYPE sf = FBUtility::convK2ND(cmp[n].getTemp(), BaseTemp, DiffTemp);                // 保持されている温度はKelvin
  
  cmp[n].getBbox(st, ed);
  
#pragma omp parallel for firstprivate(ix, jx, kx, gd, st, ed, odr, ht, sf) private(q) schedule(static) reduction(+:va)
  for (int k=st[2]; k<=ed[2]; k++) {
    for (int j=st[1]; j<=ed[1]; j++) {
      for (int i=st[0]; i<=ed[0]; i++) {
        
        size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        int s = d_bh1[m];
        
        if ( GET_FACE_BC(s, BC_FACE_W) == odr )
        {
          q = ht * (sf - d_t0[m]);                                 // 基準温度との温度差
          d_qbc[_F_IDX_V3DEX(0, i-1, j, k, ix, jx, kx, gd)] += q;
          va += q;                                                 // マイナス面の正の値は流入
        }
        
        if ( GET_FACE_BC(s, BC_FACE_E) == odr )
        {
          q = ht * (d_t0[m] - sf);
          d_qbc[_F_IDX_V3DEX(0, i, j, k, ix, jx, kx, gd)] += q;
          va -= q;                                                 // プラス面の正の値は流出
        }
        
        if ( GET_FACE_BC(s, BC_FACE_S) == odr )
        {
          q = ht * (sf - d_t0[m]);
          d_qbc[_F_IDX_V3DEX(1, i, j-1, k, ix, jx, kx, gd)] += q;
          va += q;
        }
        
        if ( GET_FACE_BC(s, BC_FACE_N) == odr )
        {
          q = ht * (d_t0[m] - sf);
          d_qbc[_F_IDX_V3DEX(1, i, j, k, ix, jx, kx, gd)] += q;
          va -= q;
        }
        
        if ( GET_FACE_BC(s, BC_FACE_B) == odr )
        {
          q = ht * (sf - d_t0[m]);
          d_qbc[_F_IDX_V3DEX(2, i, j, k-1, ix, jx, kx, gd)] += q;
          va += q;
        }
        
        if ( GET_FACE_BC(s, BC_FACE_T) == odr )
        {
          q = ht * (d_t0[m] - sf);
          d_qbc[_F_IDX_V3DEX(2, i, j, k, ix, jx, kx, gd)] += q;
          va -= q;
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
  
  return va; // 単位面積あたりの無次元熱流束の和
}



// #################################################################
/**
 * @brief 外部領域のOutflow, In_out, TractionFreeの境界条件処理
 * @retval 無次元熱量
 * @param [in,out] d_ws  温度増分
 * @param [in]     d_bh  BCindex H1
 * @param [in]     face  外部境界面番号
 * @param [in]     d_vf  セルフェイス速度
 * @param [in]     d_t   n時刻の温度
 * @param [in]     v00   参照速度
 * @note モニタ量の熱量 va は系に対する流入量なので，基準温度に対する熱量
 */
REAL_TYPE SetBC3D::psObcFree(REAL_TYPE* d_ws, const int* d_bh, const int face, const REAL_TYPE* d_vf, const REAL_TYPE* d_t, const REAL_TYPE* v00)
{
  if ( nID[face] >= 0 ) return 0.0;
  
  
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
  REAL_TYPE va = 0.0;
  REAL_TYPE dh = deltaX;
  REAL_TYPE f_e, f_w, f_n, f_s, f_t, f_b, ff;
  REAL_TYPE dh1 = 1.0/dh;
  REAL_TYPE u_ref = v00[1];
  REAL_TYPE v_ref = v00[2];
  REAL_TYPE w_ref = v00[3];
  //REAL_TYPE t_nd = FBUtility::convK2ND(obc[face].getTemp(), BaseTemp, DiffTemp);
  
  f_e = f_w = f_n = f_s = f_t = f_b = 0.0; // 境界条件以外の寄与をゼロにしておく
  
  // ff = f_e - f_w + f_n - f_s + f_t - f_b として該当流束の寄与のみを考える
  
  switch (face)
  {
    case X_MINUS:

#pragma omp parallel for firstprivate(ix, jx, kx, gd, dh1, f_w, u_ref) \
schedule(static) reduction(+:va)
      for (int k=1; k<=kx; k++) {
        for (int j=1; j<=jx; j++) {
          size_t m = _F_IDX_S3D(1, j, k, ix, jx, kx, gd);
          int s = d_bh[m];
          size_t m_w = _F_IDX_V3D(0, j, k, 0, ix, jx, kx, gd);
          REAL_TYPE c = d_vf[m_w] - u_ref;
          REAL_TYPE t_p = (c < 0.0) ? d_t[m] : 0.0;           // 流出の場合には内部の値，流入の場合には断熱
          f_w = c*t_p;
          va += f_w;
          REAL_TYPE ff = - f_w;
          d_ws[m] -= ff*dh1*GET_SHIFT_F(s, STATE_BIT);
        }
      }
      break;
      
    case X_PLUS:

#pragma omp parallel for firstprivate(ix, jx, kx, gd, dh1, f_e, u_ref) \
schedule(static) reduction(+:va)
      for (int k=1; k<=kx; k++) {
        for (int j=1; j<=jx; j++) {
          size_t m = _F_IDX_S3D(ix, j, k, ix, jx, kx, gd);
          int s = d_bh[m];
          size_t m_e = _F_IDX_V3D(ix, j, k, 0, ix, jx, kx, gd);
          REAL_TYPE c = d_vf[m_e] - u_ref;
          REAL_TYPE t_p = (c > 0.0) ? d_t[m] : 0.0;
          f_e = c*t_p;
          va += f_e;
          REAL_TYPE ff = f_e;
          d_ws[m] -= ff*dh1*GET_SHIFT_F(s, STATE_BIT);
        }
      }
      break;
      
    case Y_MINUS:
      
#pragma omp parallel for firstprivate(ix, jx, kx, gd, dh1, f_s, v_ref) \
schedule(static) reduction(+:va)
      for (int k=1; k<=kx; k++) {
        for (int i=1; i<=ix; i++) {
          size_t m = _F_IDX_S3D(i, 1, k, ix, jx, kx, gd);
          int s = d_bh[m];
          size_t m_s = _F_IDX_V3D(i, 0, k, 1, ix, jx, kx, gd);
          REAL_TYPE c = d_vf[m_s] - v_ref;
          REAL_TYPE t_p = (c < 0.0) ? d_t[m] : 0.0;
          f_s = c*t_p;
          va += f_s;
          REAL_TYPE ff = - f_s;
          d_ws[m] -= ff*dh1*GET_SHIFT_F(s, STATE_BIT);
        }
      }
      break;
      
    case Y_PLUS:
      
#pragma omp parallel for firstprivate(ix, jx, kx, gd, dh1, f_n, v_ref) \
schedule(static) reduction(+:va)
      for (int k=1; k<=kx; k++) {
        for (int i=1; i<=ix; i++) {
          size_t m = _F_IDX_S3D(i, jx, k, ix, jx, kx, gd);
          int s = d_bh[m];
          size_t m_n = _F_IDX_V3D(i, jx, k, 1, ix, jx, kx, gd);
          REAL_TYPE c = d_vf[m_n] - v_ref;
          REAL_TYPE t_p = (c > 0.0) ? d_t[m] : 0.0;
          f_n = c*t_p;
          va += f_n;
          REAL_TYPE ff = f_n;
          d_ws[m] -= ff*dh1*GET_SHIFT_F(s, STATE_BIT);
        }
      }
      break;
      
    case Z_MINUS:
      
#pragma omp parallel for firstprivate(ix, jx, kx, gd, dh1, f_b, w_ref) \
schedule(static) reduction(+:va)
      for (int j=1; j<=jx; j++) {
        for (int i=1; i<=ix; i++) {
          size_t m = _F_IDX_S3D(i, j, 1, ix, jx, kx, gd);
          int s = d_bh[m];
          size_t m_b = _F_IDX_V3D(i, j, 0, 2, ix, jx, kx, gd);
          REAL_TYPE c = d_vf[m_b] - w_ref;
          REAL_TYPE t_p = (c < 0.0) ? d_t[m] : 0.0;
          f_b = c*t_p;
          va += f_b;
          REAL_TYPE ff = - f_b;
          d_ws[m] -= ff*dh1*GET_SHIFT_F(s, STATE_BIT);
        }
      }
      break;
      
    case Z_PLUS:
      
#pragma omp parallel for firstprivate(ix, jx, kx, gd, dh1, f_t, w_ref) \
schedule(static) reduction(+:va)
      for (int j=1; j<=jx; j++) {
        for (int i=1; i<=ix; i++) {
          size_t m = _F_IDX_S3D(i, j, kx, ix, jx, kx, gd);
          int s = d_bh[m];
          size_t m_t = _F_IDX_V3D(i, j, kx, 2, ix, jx, kx, gd);
          REAL_TYPE c = d_vf[m_t] - w_ref;
          REAL_TYPE t_p = (c > 0.0) ? d_t[m] : 0.0;
          f_t = c*t_p;
          va += f_t;
          REAL_TYPE ff = f_t;
          d_ws[m] -= ff*dh1*GET_SHIFT_F(s, STATE_BIT);
        }
      }
      break;
  }
  
  return va;
}


// #################################################################
/**
 * @brief 外部領域の熱流束指定の境界条件処理
 * @retval 無次元熱量
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
  REAL_TYPE q = obc[face].getHeatflux() / (RefV*DiffTemp*rho*cp); // [W/m^2]を無次元化
  
  switch (face)
  {
    case X_MINUS:
#pragma omp parallel for firstprivate(ix, jx, kx, gd, q) schedule(static) reduction(+:va)
      for (int k=1; k<=kx; k++) {
        for (int j=1; j<=jx; j++) {
          d_qbc[_F_IDX_V3DEX(0, 0, j, k, ix, jx, kx, gd)] += q;
          va += q;
        }
      }
      break;
      
    case X_PLUS:
#pragma omp parallel for firstprivate(ix, jx, kx, gd, q) schedule(static) reduction(+:va)
      for (int k=1; k<=kx; k++) {
        for (int j=1; j<=jx; j++) {
          d_qbc[_F_IDX_V3DEX(0, ix, j, k, ix, jx, kx, gd)] -= q; // プラス面の正の値は流出
          va -= q;
        }
      }
      break;
      
    case Y_MINUS:
#pragma omp parallel for firstprivate(ix, jx, kx, gd, q) schedule(static) reduction(+:va)
      for (int k=1; k<=kx; k++) {
        for (int i=1; i<=ix; i++) {
          d_qbc[_F_IDX_V3DEX(1, i, 0, k, ix, jx, kx, gd)] += q;
          va += q;
        }
      }
      break;
      
    case Y_PLUS:
#pragma omp parallel for firstprivate(ix, jx, kx, gd, q) schedule(static) reduction(+:va)
      for (int k=1; k<=kx; k++) {
        for (int i=1; i<=ix; i++) {
          d_qbc[_F_IDX_V3DEX(1, i, jx, k, ix, jx, kx, gd)] -= q;
          va -= q;
        }
      }
      break;
      
    case Z_MINUS:
#pragma omp parallel for firstprivate(ix, jx, kx, gd, q) schedule(static) reduction(+:va)
      for (int j=1; j<=jx; j++) {
        for (int i=1; i<=ix; i++) {
          d_qbc[_F_IDX_V3DEX(2, i, j, 0, ix, jx, kx, gd)] += q;
          va += q;
        }
      }
      break;
      
    case Z_PLUS:
#pragma omp parallel for firstprivate(ix, jx, kx, gd, q) schedule(static) reduction(+:va)
      for (int j=1; j<=jx; j++) {
        for (int i=1; i<=ix; i++) {
          d_qbc[_F_IDX_V3DEX(2, i, j, kx, ix, jx, kx, gd)] -= q;
          va -= q;
        }
      }
      break;
  }
  
  return va;
}


// #################################################################
/**
 * @brief 外部領域の熱伝達境界の境界条件処理 (TypeB/S共用)
 * @retval 無次元熱量
 * @param [in,out] d_qbc 境界条件熱流束
 * @param [in]     face  外部境界面番号
 * @param [out]    d_t   n+1時刻の温度場
 * @param [in]     d_t0  n時刻の温度場
 * @note
    - モニタ量の熱量 va は系に対する流入出量
    - 流入を正にとる
 */
REAL_TYPE SetBC3D::psObcHeatTransferBS(REAL_TYPE* d_qbc, const int face, REAL_TYPE* d_t, const REAL_TYPE* d_t0)
{
  if ( nID[face] >= 0 ) return 0.0;
  
  
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
  REAL_TYPE va = 0.0;
  REAL_TYPE ht = obc[face].getCoefHT() / (RefV*rho*cp);                        // 熱伝達係数
  REAL_TYPE bt = FBUtility::convK2ND(obc[face].getTemp(), BaseTemp, DiffTemp); // 温度．保持されている温度はKelvin
  
  switch (face)
  {
    case X_MINUS:
#pragma omp parallel for firstprivate(ix, jx, kx, gd, ht, bt) schedule(static) reduction(+:va)
      for (int k=1; k<=kx; k++) {
        for (int j=1; j<=jx; j++) {
          size_t m = _F_IDX_S3D(1, j, k, ix, jx, kx, gd);
          REAL_TYPE q = ht * (bt - d_t0[m]);
          d_qbc[_F_IDX_V3DEX(0, 0, j, k, ix, jx, kx, gd)] += q; // マイナス面の正の値は流入
          va += q;
          d_t[_F_IDX_S3D(0, j, k, ix, jx, kx, gd)] = bt; // 隣接流体セルにバルク温度を代入
        }
      }
      break;
      
    case X_PLUS:
#pragma omp parallel for firstprivate(ix, jx, kx, gd, ht, bt) schedule(static) reduction(+:va)
      for (int k=1; k<=kx; k++) {
        for (int j=1; j<=jx; j++) {
          size_t m = _F_IDX_S3D(ix, j, k, ix, jx, kx, gd);
          REAL_TYPE q = ht * (d_t0[m] - bt);
          d_qbc[_F_IDX_V3DEX(0, ix, j, k, ix, jx, kx, gd)] -= q; // プラス面の正の値は流出
          va -= q;
          d_t[_F_IDX_S3D(ix+1, j, k, ix, jx, kx, gd)] = bt;
        }
      }
      break;
      
    case Y_MINUS:
#pragma omp parallel for firstprivate(ix, jx, kx, gd, ht, bt) schedule(static) reduction(+:va)
      for (int k=1; k<=kx; k++) {
        for (int i=1; i<=ix; i++) {
          size_t m = _F_IDX_S3D(i, 1, k, ix, jx, kx, gd);
          REAL_TYPE q = ht * (bt - d_t0[m]);
          d_qbc[_F_IDX_V3DEX(1, i, 0, k, ix, jx, kx, gd)] += q;
          va += q;
          d_t[_F_IDX_S3D(i, 0, k, ix, jx, kx, gd)] = bt;
        }
      }
      break;
      
    case Y_PLUS:
#pragma omp parallel for firstprivate(ix, jx, kx, gd, ht, bt) schedule(static) reduction(+:va)
      for (int k=1; k<=kx; k++) {
        for (int i=1; i<=ix; i++) {
          size_t m = _F_IDX_S3D(i, jx, k, ix, jx, kx, gd);
          REAL_TYPE q = ht * (d_t0[m] - bt);
          d_qbc[_F_IDX_V3DEX(1, i, jx, k, ix, jx, kx, gd)] -= q;
          va -= q;
          d_t[_F_IDX_S3D(i, jx+1, k, ix, jx, kx, gd)] = bt;
        }
      }
      break;
      
    case Z_MINUS:
#pragma omp parallel for firstprivate(ix, jx, kx, gd, ht, bt) schedule(static) reduction(+:va)
      for (int j=1; j<=jx; j++) {
        for (int i=1; i<=ix; i++) {
          size_t m = _F_IDX_S3D(i, j, 1, ix, jx, kx, gd);
          REAL_TYPE q = ht * (bt - d_t0[m]);
          d_qbc[_F_IDX_V3DEX(2, i, j, 0, ix, jx, kx, gd)] += q;
          va += q;
          d_t[_F_IDX_S3D(i, j, 0, ix, jx, kx, gd)] = bt;
        }
      }
      break;
      
    case Z_PLUS:
#pragma omp parallel for firstprivate(ix, jx, kx, gd, ht, bt) schedule(static) reduction(+:va)
      for (int j=1; j<=jx; j++) {
        for (int i=1; i<=ix; i++) {
          size_t m = _F_IDX_S3D(i, j, kx, ix, jx, kx, gd);
          REAL_TYPE q = ht * (d_t0[m] - bt);
          d_qbc[_F_IDX_V3DEX(2, i, j, kx, ix, jx, kx, gd)] -= q;
          va -= q;
          d_t[_F_IDX_S3D(i, j, kx+1, ix, jx, kx, gd)] = bt;
        }
      }
      break;
  }
  
  return va;
}


// #################################################################
/**
 * @brief 外部領域の熱伝達境界の境界条件処理 (Type SF)
 * @retval 無次元熱量
 * @param [in,out] d_qbc 境界条件熱流束
 * @param [in]     face  外部境界面番号
 * @param [out]    d_t   n+1時刻の温度場
 * @param [in]     d_t0  n時刻の温度場
 * @note
    - モニタ量の熱量va(-)は系に対する流入出量
    - 流入を正にとる
 */
REAL_TYPE SetBC3D::psObcHeatTransferSF(REAL_TYPE* d_qbc, const int face, REAL_TYPE* d_t, const REAL_TYPE* d_t0)
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
  REAL_TYPE ht = lambda / (RefV*RefL*rho*cp) * a1*pow(Reynolds, b1) * pow(Prandtl, c1);  // 熱伝達係数
  REAL_TYPE sf = FBUtility::convK2ND(obc[face].getTemp(), BaseTemp, DiffTemp);           // 表面温度．保持されている温度はKelvin
  
  switch (face)
  {
    case X_MINUS:
#pragma omp parallel for firstprivate(ix, jx, kx, gd, ht, sf) schedule(static) reduction(+:va)
      for (int k=1; k<=kx; k++) {
        for (int j=1; j<=jx; j++) {
          size_t m = _F_IDX_S3D(1, j, k, ix, jx, kx, gd);
          REAL_TYPE q = ht * (sf - d_t0[m]);                                 // 表面温度と隣接セルとの温度差
          d_qbc[_F_IDX_V3DEX(0, 0, j, k, ix, jx, kx, gd)] += q;
          va += q;
          d_t[_F_IDX_S3D(0, j, k, ix, jx, kx, gd)] = sf;
        }
      }
      break;
      
    case X_PLUS:
#pragma omp parallel for firstprivate(ix, jx, kx, gd, ht, sf) schedule(static) reduction(+:va)
      for (int k=1; k<=kx; k++) {
        for (int j=1; j<=jx; j++) {
          size_t m = _F_IDX_S3D(ix, j, k, ix, jx, kx, gd);
          REAL_TYPE q = ht * (d_t0[m] - sf);
          d_qbc[_F_IDX_V3DEX(0, ix, j, k, ix, jx, kx, gd)] -= q;
          va -= q;                                                 // プラス面の正の値は流出
          d_t[_F_IDX_S3D(ix+1, j, k, ix, jx, kx, gd)] = sf;
        }
      }
      break;
      
    case Y_MINUS:
#pragma omp parallel for firstprivate(ix, jx, kx, gd, ht, sf) schedule(static) reduction(+:va)
      for (int k=1; k<=kx; k++) {
        for (int i=1; i<=ix; i++) {
          size_t m = _F_IDX_S3D(i, 1, k, ix, jx, kx, gd);
          REAL_TYPE q = ht * (sf - d_t0[m]);
          d_qbc[_F_IDX_V3DEX(1, i, 0, k, ix, jx, kx, gd)] += q;
          va += q;
          d_t[_F_IDX_S3D(i, 0, k, ix, jx, kx, gd)] = sf;
        }
      }
      break;
      
    case Y_PLUS:
#pragma omp parallel for firstprivate(ix, jx, kx, gd, ht, sf) schedule(static) reduction(+:va)
      for (int k=1; k<=kx; k++) {
        for (int i=1; i<=ix; i++) {
          size_t m = _F_IDX_S3D(i, jx, k, ix, jx, kx, gd);
          REAL_TYPE q = ht * (d_t0[m] - sf);
          d_qbc[_F_IDX_V3DEX(1, i, jx, k, ix, jx, kx, gd)] -= q;
          va -= q;
          d_t[_F_IDX_S3D(i, jx+1, k, ix, jx, kx, gd)] = sf;
        }
      }
      break;
      
    case Z_MINUS:
#pragma omp parallel for firstprivate(ix, jx, kx, gd, ht, sf) schedule(static) reduction(+:va)
      for (int j=1; j<=jx; j++) {
        for (int i=1; i<=ix; i++) {
          size_t m = _F_IDX_S3D(i, j, 1, ix, jx, kx, gd);
          REAL_TYPE q = ht * (sf - d_t0[m]);
          d_qbc[_F_IDX_V3DEX(2, i, j, 0, ix, jx, kx, gd)] += q;
          va += q;
          d_t[_F_IDX_S3D(i, j, 0, ix, jx, kx, gd)] = sf;
        }
      }
      break;
      
    case Z_PLUS:
#pragma omp parallel for firstprivate(ix, jx, kx, gd, ht, sf) schedule(static) reduction(+:va)
      for (int j=1; j<=jx; j++) {
        for (int i=1; i<=ix; i++) {
          size_t m = _F_IDX_S3D(i, j, kx, ix, jx, kx, gd);
          REAL_TYPE q = ht * (d_t0[m] - sf);
          d_qbc[_F_IDX_V3DEX(2, i, j, kx, ix, jx, kx, gd)] -= q;
          va -= q;
          d_t[_F_IDX_S3D(i, j, kx+1, ix, jx, kx, gd)] = sf;
        }
      }
      break;
  }
  
  return va;
}


// #################################################################
/**
 * @brief 外部領域の熱伝達境界の境界条件処理 (Type SN)
 * @retval 無次元熱量
 * @param [in,out] d_qbc 境界条件熱流束
 * @param [in]     face  外部境界面番号
 * @param [out]    d_t   n+1時刻の温度場
 * @param [in]     d_t0  n時刻の温度場
 * @note
   - モニタ量の熱量 va は系に対する流入出量
   - 流入を正にとる
 */
REAL_TYPE SetBC3D::psObcHeatTransferSN(REAL_TYPE* d_qbc, const int face, REAL_TYPE* d_t, const REAL_TYPE* d_t0)
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
  REAL_TYPE ht = lambda / (RefV*RefL*rho*cp);                    // RefMediumで指定される物性値
  
  REAL_TYPE ht1 = ( Rayleigh < c1 ) ? ht*a1*pow(Rayleigh, b1) : ht*a2*pow(Rayleigh, b2); // vertical, upper
  REAL_TYPE ht3 = ( Rayleigh < c2 ) ? ht*a3*pow(Rayleigh, b3) : ht*a4*pow(Rayleigh, b4); // lower
  REAL_TYPE sf = FBUtility::convK2ND(obc[face].getTemp(), BaseTemp, DiffTemp);           // 保持されている温度はKelvin
  
  switch (face)
  {
    case X_MINUS:
#pragma omp parallel for firstprivate(ix, jx, kx, gd, ht1, sf) schedule(static) reduction(+:va)
      for (int k=1; k<=kx; k++) {
        for (int j=1; j<=jx; j++) {
          size_t m = _F_IDX_S3D(1, j, k, ix, jx, kx, gd);
          REAL_TYPE q = ht1 * (sf - d_t0[m]);                    // 基準温度との温度差
          d_qbc[_F_IDX_V3DEX(0, 0, j, k, ix, jx, kx, gd)] += q;
          va += q;
          d_t[_F_IDX_S3D(0, j, k, ix, jx, kx, gd)] = sf;
        }
      }
      break;
      
    case X_PLUS:
#pragma omp parallel for firstprivate(ix, jx, kx, gd, ht1, sf) schedule(static) reduction(+:va)
      for (int k=1; k<=kx; k++) {
        for (int j=1; j<=jx; j++) {
          size_t m = _F_IDX_S3D(ix, j, k, ix, jx, kx, gd);
          REAL_TYPE q = ht1 * (d_t0[m] - sf);
          d_qbc[_F_IDX_V3DEX(0, ix, j, k, ix, jx, kx, gd)] -= q;
          va -= q;                                               // プラス面の正の値は流出
          d_t[_F_IDX_S3D(ix+1, j, k, ix, jx, kx, gd)] = sf;
        }
      }
      break;
      
    case Y_MINUS:
#pragma omp parallel for firstprivate(ix, jx, kx, gd, ht1, sf) schedule(static) reduction(+:va)
      for (int k=1; k<=kx; k++) {
        for (int i=1; i<=ix; i++) {
          size_t m = _F_IDX_S3D(i, 1, k, ix, jx, kx, gd);
          REAL_TYPE q = ht1 * (sf - d_t0[m]);
          d_qbc[_F_IDX_V3DEX(1, i, 0, k, ix, jx, kx, gd)] += q;
          va += q;
          d_t[_F_IDX_S3D(i, 0, k, ix, jx, kx, gd)] = sf;
        }
      }
      break;
      
    case Y_PLUS:
#pragma omp parallel for firstprivate(ix, jx, kx, gd, ht1, sf) schedule(static) reduction(+:va)
      for (int k=1; k<=kx; k++) {
        for (int i=1; i<=ix; i++) {
          size_t m = _F_IDX_S3D(i, jx, k, ix, jx, kx, gd);
          REAL_TYPE q = ht1 * (d_t0[m] - sf);
          d_qbc[_F_IDX_V3DEX(1, i, jx, k, ix, jx, kx, gd)] -= q;
          va -= q;
          d_t[_F_IDX_S3D(i, jx+1, k, ix, jx, kx, gd)] = sf;
        }
      }
      break;
      
    case Z_MINUS:
      // ここだけ係数が違うので注意
#pragma omp parallel for firstprivate(ix, jx, kx, gd, ht3, sf) schedule(static) reduction(+:va)
      for (int j=1; j<=jx; j++) {
        for (int i=1; i<=ix; i++) {
          size_t m = _F_IDX_S3D(i, j, 1, ix, jx, kx, gd);
          REAL_TYPE q = ht3 * (sf - d_t0[m]);
          d_qbc[_F_IDX_V3DEX(2, i, j, 0, ix, jx, kx, gd)] += q;
          va += q;
          d_t[_F_IDX_S3D(i, j, 0, ix, jx, kx, gd)] = sf;
        }
      }
      break;
      
    case Z_PLUS:
#pragma omp parallel for firstprivate(ix, jx, kx, gd, ht1, sf) schedule(static) reduction(+:va)
      for (int j=1; j<=jx; j++) {
        for (int i=1; i<=ix; i++) {
          size_t m = _F_IDX_S3D(i, j, kx, ix, jx, kx, gd);
          REAL_TYPE q = ht1 * (d_t0[m] - sf);
          d_qbc[_F_IDX_V3DEX(2, i, j, kx, ix, jx, kx, gd)] -= q;
          va -= q;
          d_t[_F_IDX_S3D(i, j, kx+1, ix, jx, kx, gd)] = sf;
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
  ht = cmp[n].getCoefHT() / (RefV*rho*cp);

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
  c = qsum * (RefV*DiffTemp*rho*cp);
  
  return (c);
}*/



// #################################################################
/**
 * @brief 外部領域の等温熱流束の境界条件処理
 * @retval 無次元熱量
 * @param [in,out] d_qbc 境界条件熱流束
 * @param [in]     face  外部境界面番号
 * @param [out]    d_t   n+1時刻の温度場
 * @param [in]     d_t0  n時刻の温度場
 * @note
   - モニタ量の熱量 va は系に対する流入出量
   - 流入を正にとる
 */
REAL_TYPE SetBC3D::psObcIsoThermal(REAL_TYPE* d_qbc, const int face, REAL_TYPE* d_t, const REAL_TYPE* d_t0)
{
  if ( nID[face] >= 0 ) return 0.0;
  
  
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
  REAL_TYPE va = 0.0;
  REAL_TYPE sf = FBUtility::convK2ND(obc[face].getTemp(), BaseTemp, DiffTemp); // 保持されている温度はKelvin
  REAL_TYPE pp = (2.0*lambda) / (deltaX*RefV*rho*cp);
  
  switch (face)
  {
    case X_MINUS:
#pragma omp parallel for firstprivate(ix, jx, kx, gd, pp, sf) schedule(static) reduction(+:va)
      for (int k=1; k<=kx; k++) {
        for (int j=1; j<=jx; j++) {
          size_t m = _F_IDX_S3D(1, j, k, ix, jx, kx, gd);
          REAL_TYPE q = pp * (sf - d_t0[m]);
          d_qbc[_F_IDX_V3DEX(0, 0, j, k, ix, jx, kx, gd)] += q; // マイナス面の正の値は流入
          va += q;
          d_t[_F_IDX_S3D(0, j, k, ix, jx, kx, gd)] = sf; // 等温セルに指定温度を代入
        }
      }
      break;
      
    case X_PLUS:
#pragma omp parallel for firstprivate(ix, jx, kx, gd, pp, sf) schedule(static) reduction(+:va)
      for (int k=1; k<=kx; k++) {
        for (int j=1; j<=jx; j++) {
          size_t m = _F_IDX_S3D(ix, j, k, ix, jx, kx, gd);
          REAL_TYPE q = pp * (d_t0[m] - sf);
          d_qbc[_F_IDX_V3DEX(0, ix, j, k, ix, jx, kx, gd)] -= q; // プラス面の正の値は流出
          va -= q;
          d_t[_F_IDX_S3D(ix+1, j, k, ix, jx, kx, gd)] = sf;
        }
      }
      break;
      
    case Y_MINUS:
#pragma omp parallel for firstprivate(ix, jx, kx, gd, pp, sf) schedule(static) reduction(+:va)
      for (int k=1; k<=kx; k++) {
        for (int i=1; i<=ix; i++) {
          size_t m = _F_IDX_S3D(i, 1, k, ix, jx, kx, gd);
          REAL_TYPE q = pp * (sf - d_t0[m]);
          d_qbc[_F_IDX_V3DEX(1, i, 0, k, ix, jx, kx, gd)] += q;
          va += q;
          d_t[_F_IDX_S3D(i, 0, k, ix, jx, kx, gd)] = sf;
        }
      }
      break;
      
    case Y_PLUS:
#pragma omp parallel for firstprivate(ix, jx, kx, gd, pp, sf) schedule(static) reduction(+:va)
      for (int k=1; k<=kx; k++) {
        for (int i=1; i<=ix; i++) {
          size_t m = _F_IDX_S3D(i, jx, k, ix, jx, kx, gd);
          REAL_TYPE q = pp * (d_t0[m] - sf);
          d_qbc[_F_IDX_V3DEX(1, i, jx, k, ix, jx, kx, gd)] -= q;
          va -= q;
          d_t[_F_IDX_S3D(i, jx+1, k, ix, jx, kx, gd)] = sf;
        }
      }
      break;
      
    case Z_MINUS:
#pragma omp parallel for firstprivate(ix, jx, kx, gd, pp, sf) schedule(static) reduction(+:va)
      for (int j=1; j<=jx; j++) {
        for (int i=1; i<=ix; i++) {
          size_t m = _F_IDX_S3D(i, j, 1, ix, jx, kx, gd);
          REAL_TYPE q = pp * (sf - d_t0[m]);
          d_qbc[_F_IDX_V3DEX(2, i, j, 0, ix, jx, kx, gd)] += q;
          va += q;
          d_t[_F_IDX_S3D(i, j, 0, ix, jx, kx, gd)] = sf;
        }
      }
      break;
      
    case Z_PLUS:
#pragma omp parallel for firstprivate(ix, jx, kx, gd, pp, sf) schedule(static) reduction(+:va)
      for (int j=1; j<=jx; j++) {
        for (int i=1; i<=ix; i++) {
          size_t m = _F_IDX_S3D(i, j, kx, ix, jx, kx, gd);
          REAL_TYPE q = pp * (d_t0[m] - sf);
          d_qbc[_F_IDX_V3DEX(2, i, j, kx, ix, jx, kx, gd)] -= q;
          va -= q;
          d_t[_F_IDX_S3D(i, j, kx+1, ix, jx, kx, gd)] = sf;
        }
      }
      break;
  }
  
  return va;
}


// #################################################################
/**
 * @brief 外部領域の速度指定の境界条件処理
 * @retval 無次元熱量
 * @param [in,out] d_ws  温度増分
 * @param [in]     d_bh1 BCindex H1
 * @param [in]     face  外部境界面番号
 * @param [in]     tm    時刻
 * @param [in]     v00   基準速度
 * @note モニタ量の熱量 va は系に対する流入量なので，基準温度に対する熱量
 */
REAL_TYPE SetBC3D::psObcSpecVH(REAL_TYPE* d_ws, const int* d_bh1, const int face, const double tm, const REAL_TYPE* v00)
{
  if ( nID[face] >= 0 ) return 0.0;
  
  
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
  REAL_TYPE va = 0.0;
  REAL_TYPE f_e, f_w, f_n, f_s, f_t, f_b, c;
  REAL_TYPE dh = deltaX;
  REAL_TYPE dh1 = 1.0/dh;
  REAL_TYPE vec[3];
  REAL_TYPE u_ref = v00[1];
  REAL_TYPE v_ref = v00[2];
  REAL_TYPE w_ref = v00[3];
  REAL_TYPE t_nd = FBUtility::convK2ND(obc[face].getTemp(), BaseTemp, DiffTemp);
  REAL_TYPE dummy = extractVelOBC(face, vec, tm, v00);
  
  f_e = f_w = f_n = f_s = f_t = f_b = 0.0; // 境界条件以外の寄与をゼロにしておく
  
  switch (face)
  {
    case X_MINUS:
      c = vec[0] - u_ref;
      f_w = c * t_nd;
      
#pragma omp parallel for firstprivate(ix, jx, kx, gd, dh1, f_e, f_w, f_n, f_s, f_t, f_b) \
schedule(static) reduction(+:va)
      for (int k=1; k<=kx; k++) {
        for (int j=1; j<=jx; j++) {
          size_t m = _F_IDX_S3D(1, j, k, ix, jx, kx, gd);
          int s = d_bh1[m];
          va += f_w;
          REAL_TYPE ff = f_e - f_w + f_n - f_s + f_t - f_b;
          d_ws[m] -= ff*dh1*GET_SHIFT_F(s, STATE_BIT);
        }
      }
      break;
      
    case X_PLUS:
      c = vec[0] - u_ref;
      f_e = -c * t_nd; // ref_tが正のとき，X+方向のセルでは流出なので符号反転
      
#pragma omp parallel for firstprivate(ix, jx, kx, gd, dh1, f_e, f_w, f_n, f_s, f_t, f_b) \
schedule(static) reduction(+:va)
      for (int k=1; k<=kx; k++) {
        for (int j=1; j<=jx; j++) {
          size_t m = _F_IDX_S3D(ix, j, k, ix, jx, kx, gd);
          int s = d_bh1[m];
          va += f_e;
          REAL_TYPE ff = f_e - f_w + f_n - f_s + f_t - f_b;
          d_ws[m] -= ff*dh1*GET_SHIFT_F(s, STATE_BIT);
        }
      }
      break;
      
    case Y_MINUS:
      c = vec[1] - v_ref;
      f_s = c * t_nd;
      
#pragma omp parallel for firstprivate(ix, jx, kx, gd, dh1, f_e, f_w, f_n, f_s, f_t, f_b) \
schedule(static) reduction(+:va)
      for (int k=1; k<=kx; k++) {
        for (int i=1; i<=ix; i++) {
          size_t m = _F_IDX_S3D(i, 1, k, ix, jx, kx, gd);
          int s = d_bh1[m];
          va += f_s;
          REAL_TYPE ff = f_e - f_w + f_n - f_s + f_t - f_b;
          d_ws[m] -= ff*dh1*GET_SHIFT_F(s, STATE_BIT);
        }
      }
      break;
      
    case Y_PLUS:
      c = vec[1] - v_ref;
      f_n = -c * t_nd; // 符号反転
      
#pragma omp parallel for firstprivate(ix, jx, kx, gd, dh1, f_e, f_w, f_n, f_s, f_t, f_b) \
schedule(static) reduction(+:va)
      for (int k=1; k<=kx; k++) {
        for (int i=1; i<=ix; i++) {
          size_t m = _F_IDX_S3D(i, jx, k, ix, jx, kx, gd);
          int s = d_bh1[m];
          va += f_n;
          REAL_TYPE ff = f_e - f_w + f_n - f_s + f_t - f_b;
          d_ws[m] -= ff*dh1*GET_SHIFT_F(s, STATE_BIT);
        }
      }
      break;
      
    case Z_MINUS:
      c = vec[2] - w_ref;
      f_b = c * t_nd;
      
#pragma omp parallel for firstprivate(ix, jx, kx, gd, dh1, f_e, f_w, f_n, f_s, f_t, f_b) \
schedule(static) reduction(+:va)
      for (int j=1; j<=jx; j++) {
        for (int i=1; i<=ix; i++) {
          size_t m = _F_IDX_S3D(i, j, 1, ix, jx, kx, gd);
          int s = d_bh1[m];
          va += f_b;
          REAL_TYPE ff = f_e - f_w + f_n - f_s + f_t - f_b;
          d_ws[m] -= ff*dh1*GET_SHIFT_F(s, STATE_BIT);
        }
      }
      break;
      
    case Z_PLUS:
      c = vec[2] - w_ref;
      f_t = -c * t_nd; // 符号反転
      
#pragma omp parallel for firstprivate(ix, jx, kx, gd, dh1, f_e, f_w, f_n, f_s, f_t, f_b) \
schedule(static) reduction(+:va)
      for (int j=1; j<=jx; j++) {
        for (int i=1; i<=ix; i++) {
          size_t m = _F_IDX_S3D(i, j, kx, ix, jx, kx, gd);
          int s = d_bh1[m];
          va += f_t;
          REAL_TYPE ff = f_e - f_w + f_n - f_s + f_t - f_b;
          d_ws[m] -= ff*dh1*GET_SHIFT_F(s, STATE_BIT);
        }
      }
      break;
  }
  
  return va;
}


// #################################################################
/**
 * @brief 熱伝達境界条件タイプSN（自然対流）
 * @retval 熱流束の和 (W/m^2)
 * @param [out]    d_qbc  境界条件熱流束
 * @param [in]     d_bh1  BCindex H1
 * @param [in]     n      境界条件コンポーネントのエントリ番号
 * @param [in]     d_t0   n時刻の温度場
 * @note
 *    - 熱流束は加算（他の条件と合成）
 *    - pセルは流体セル
 */
REAL_TYPE SetBC3D::psIbcTransferSN (REAL_TYPE* d_qbc, const int* d_bh1, const int n, const REAL_TYPE* d_t0)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  int odr= n;
  
  REAL_TYPE q;
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
  REAL_TYPE ht = lambda / (RefV*RefL*rho*cp);                 // RefMediumで指定される物性値
  
  REAL_TYPE ht1 = ( Rayleigh < c1 ) ? ht*a1*pow(Rayleigh, b1) : ht*a2*pow(Rayleigh, b2); // vertical, upper
  REAL_TYPE ht3 = ( Rayleigh < c2 ) ? ht*a3*pow(Rayleigh, b3) : ht*a4*pow(Rayleigh, b4); // lower
  REAL_TYPE sf = FBUtility::convK2ND(cmp[n].getTemp(), BaseTemp, DiffTemp); // 保持されている温度はKelvin
  
  cmp[n].getBbox(st, ed);
    

#pragma omp parallel for firstprivate(ix, jx, kx, gd, st, ed, odr, ht1, ht3, sf) \
            private(q) \
            schedule(static) reduction(+:va)
    
  for (int k=st[2]; k<=ed[2]; k++) {
    for (int j=st[1]; j<=ed[1]; j++) {
      for (int i=st[0]; i<=ed[0]; i++) {
        
        size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        int s = d_bh1[m];
        
        if ( GET_FACE_BC(s, BC_FACE_W) == odr )
        {
          q = ht1 * (sf - d_t0[m]); // 基準温度との温度差
          d_qbc[_F_IDX_V3DEX(0, i-1, j, k, ix, jx, kx, gd)] += q;
          va += q; // マイナス面の正の値は流入
        }
        
        if ( GET_FACE_BC(s, BC_FACE_E) == odr )
        {
          q = ht1 * (d_t0[m] - sf);
          d_qbc[_F_IDX_V3DEX(0, i, j, k, ix, jx, kx, gd)] += q;
          va -= q; // プラス面の正の値は流出
        }
        
        if ( GET_FACE_BC(s, BC_FACE_S) == odr )
        {
          q = ht1 * (sf - d_t0[m]);
          d_qbc[_F_IDX_V3DEX(1, i, j-1, k, ix, jx, kx, gd)] += q;
          va += q;
        }
        
        if ( GET_FACE_BC(s, BC_FACE_N) == odr )
        {
          q = ht1 * (d_t0[m] - sf);
          d_qbc[_F_IDX_V3DEX(1, i, j, k, ix, jx, kx, gd)] += q;
          va -= q;
        }
        
        if ( GET_FACE_BC(s, BC_FACE_B) == odr )
        {
          q = ht3 * (sf - d_t0[m]);
          d_qbc[_F_IDX_V3DEX(2, i, j, k-1, ix, jx, kx, gd)] += q;
          va += q;
        }
        
        if ( GET_FACE_BC(s, BC_FACE_T) == odr )
        {
          q = ht1 * (d_t0[m] - sf);
          d_qbc[_F_IDX_V3DEX(2, i, j, k, ix, jx, kx, gd)] += q;
          va -= q;
        }
      }
    }
  }

  if ( numProc > 1 )
  {
    REAL_TYPE tmp = va;
    if ( paraMngr->Allreduce(&tmp, &va, 1, MPI_SUM) != CPM_SUCCESS ) Exit(0);
  }
  
  return va; // 単位面積あたりの無次元熱流束の和
}



// #################################################################
/**
 @brief 単媒質の場合の熱伝達境界条件タイプB
 @retval sum of heat flux (W/m^2)
 @param qbc 境界条件熱流束
 @param t n+1時刻の温度場
 @param bx BCindex
 @param n 境界条件コンポーネントのエントリ番号
 @param t0 n時刻の温度場
 
REAL_TYPE SetBC3D::setHeatTransferB(REAL_TYPE* qbc, REAL_TYPE* t, int* bx, int n, REAL_TYPE* t0)
{
  int i,j,k;
  int s, m0;
  int wi, wj, wk;
  int li, lj, lk;
  int ni, nj, nk;
  REAL_TYPE qi, qj, qk, pp, bt, qsum=0.0, tmp, c;
  REAL_TYPE di, dj, dk; // outer normal
  
  pp = cmp[n].getCoefHT() / RefV;
  bt = FBUtility::convK2ND(cmp[n].getTemp(), BaseTemp, DiffTemp); // 保持されている温度はKelvin
  
  for (k=cmp[n].ci.st[2]-1; k<=cmp[n].ci.ed[2]; k++) {
    for (j=cmp[n].ci.st[1]-1; j<=cmp[n].ci.ed[1]; j++) {
      for (i=cmp[n].ci.st[0]-1; i<=cmp[n].ci.ed[0]; i++) {
        m0 = FBUtility::getFindexS3D(size, guide, i  , j  , k  );
        s = bx[m0];
        if ( (s & BIT_MASK_10) == n ) {
          
          if ( TEST_BIT(s, QFACE_I) ) {
            wi = BIT_SHIFT(s, DIR_I);
            t[ FBUtility::getFindexS3D(size, guide, i+wi, j, k) ] = bt; // 流体セルにバルク温度を代入
            li = FBUtility::getFindexS3D(size, guide, i+1-wi, j, k); // 固体セルのセルインデクス
            ni = DECODE_MAT( bx[li] ); // 固体セルのマテリアルのオーダーを得る
            di = 2.0*(REAL_TYPE)wi-1.0;
            qi = pp / (mat[ni].P[p_density] * mat[ni].P[p_specific_heat]) * (t0[li]-bt) * di;
            qbc[FBUtility::getFindexV3DEx(size, guide, 0, i, j, k)] = qi;
            qsum += qi*di;
          }
          if ( TEST_BIT(s, QFACE_J) ) {
            wj = BIT_SHIFT(s, DIR_J);
            t[ FBUtility::getFindexS3D(size, guide, i, j+wj, k) ] = bt;
            lj = FBUtility::getFindexS3D(size, guide, i, j+1-wj, k);
            nj = DECODE_MAT( bx[lj] );
            dj = 2.0*(REAL_TYPE)wj-1.0;
            qj = pp / (mat[nj].P[p_density] * mat[nj].P[p_specific_heat]) * (t0[lj]-bt) * dj;
            qbc[FBUtility::getFindexV3DEx(size, guide, 1, i, j, k)] = qj;
            qsum += qj*dj;
          }
          if ( TEST_BIT(s, QFACE_K) ) {
            wk = BIT_SHIFT(s, DIR_K);
            t[ FBUtility::getFindexS3D(size, guide, i, j, k+wk) ] = bt;
            lk = FBUtility::getFindexS3D(size, guide, i, j, k+1-wk);
            nk = DECODE_MAT( bx[lk] );
            dk = 2.0*(REAL_TYPE)wk-1.0;
            qk = pp / (mat[nk].P[p_density] * mat[nk].P[p_specific_heat]) * (t0[lk]-bt) * dk;
            qbc[FBUtility::getFindexV3DEx(size, guide, 2, i, j, k)] = qk;
            qsum += qk*dk;
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
  c = qsum * (RefV*DiffTemp*rho*cp);
  
  return (c);
}*/



// #################################################################
// 初期温度を代入
void SetBC3D::setInitialTempCompo(const int n, const int* d_bx, REAL_TYPE* d_t)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
  REAL_TYPE ref = FBUtility::convK2ND(cmp[n].getInitTemp(), BaseTemp, DiffTemp);
	
  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        if ( DECODE_CMP(d_bx[m]) == n )
        {
          d_t[m] = ref; 
        }
      }
    }
  }
}


// #################################################################
/**
 * @brief 周期境界の場合のインデクスの同期
 * @param [in,out] bx BCindexのデータクラス
 */
void SetBC3D::setBCIperiodic(int* d_bx)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
  size_t m0, m1;

  if ( numProc > 1 ) 
  {
    for (int face=0; face<NOFACE; face++)
    {
      if ( obc[face].getClass() != OBC_PERIODIC ) continue; // スキップしてfaceをインクリメント

      switch (face)
      {
        case X_MINUS:
          if ( paraMngr->PeriodicCommS3D(d_bx, ix, jx, kx, gd, 1, X_DIR, PLUS2MINUS) != CPM_SUCCESS ) Exit(0);
          break;
          
        case X_PLUS:
          if ( paraMngr->PeriodicCommS3D(d_bx, ix, jx, kx, gd, 1, X_DIR, MINUS2PLUS) != CPM_SUCCESS ) Exit(0);
          break;
          
        case Y_MINUS:
          if ( paraMngr->PeriodicCommS3D(d_bx, ix, jx, kx, gd, 1, Y_DIR, PLUS2MINUS) != CPM_SUCCESS ) Exit(0);
          break;
          
        case Y_PLUS:
          if ( paraMngr->PeriodicCommS3D(d_bx, ix, jx, kx, gd, 1, Y_DIR, MINUS2PLUS) != CPM_SUCCESS ) Exit(0);
          break;
          
        case Z_MINUS:
          if ( paraMngr->PeriodicCommS3D(d_bx, ix, jx, kx, gd, 1, Z_DIR, PLUS2MINUS) != CPM_SUCCESS ) Exit(0);
          break;
          
        case Z_PLUS:
          if ( paraMngr->PeriodicCommS3D(d_bx, ix, jx, kx, gd, 1, Z_DIR, MINUS2PLUS) != CPM_SUCCESS ) Exit(0);
          break;
      }
    }
  } 
  else // 逐次処理
  {
    for (int face=0; face<NOFACE; face++)
    {
      if ( obc[face].getClass() != OBC_PERIODIC ) continue; // スキップしてfaceをインクリメント
      
      if ( nID[face] >= 0 ) continue;
      
      switch (face)
      {
        case X_MINUS:
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
          break;
          
        case X_PLUS:
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
          break;
          
        case Y_MINUS:
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
          break;
          
        case Y_PLUS:
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
          break;
          
        case Z_MINUS:
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
          break;
          
        case Z_PLUS:
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
          break;
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
      if ( (typ==SPEC_VEL) || (typ==SPEC_VEL_WH) ) {
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
/**
 * @brief 温度の外部周期境界条件（単純なコピー）
 * @param [in,out] d_t  温度のデータクラス
 * @param [in]     face 面番号
 */
void SetBC3D::TobcPeriodicSimple(REAL_TYPE* d_t, const int face)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
  if ( numProc > 1 ) 
  {
    switch (face) 
    {
      case X_MINUS:
        if ( paraMngr->PeriodicCommS3D(d_t, ix, jx, kx, gd, gd, X_DIR, PLUS2MINUS) != CPM_SUCCESS ) Exit(0);
        break;
        
      case X_PLUS:
        if ( paraMngr->PeriodicCommS3D(d_t, ix, jx, kx, gd, gd, X_DIR, MINUS2PLUS) != CPM_SUCCESS ) Exit(0);
        break;
        
      case Y_MINUS:
        if ( paraMngr->PeriodicCommS3D(d_t, ix, jx, kx, gd, gd, Y_DIR, PLUS2MINUS) != CPM_SUCCESS ) Exit(0);
        break;
        
      case Y_PLUS:
        if ( paraMngr->PeriodicCommS3D(d_t, ix, jx, kx, gd, gd, Y_DIR, MINUS2PLUS) != CPM_SUCCESS ) Exit(0);
        break;
        
      case Z_MINUS:
        if ( paraMngr->PeriodicCommS3D(d_t, ix, jx, kx, gd, gd, Z_DIR, PLUS2MINUS) != CPM_SUCCESS ) Exit(0);
        break;
        
      case Z_PLUS:
        if ( paraMngr->PeriodicCommS3D(d_t, ix, jx, kx, gd, gd, Z_DIR, MINUS2PLUS) != CPM_SUCCESS ) Exit(0);
        break;
    }
  }
  else  // Serial
  {
    if ( nID[face] >= 0 ) return;
    
    switch (face)
    {
      case X_MINUS:
#pragma omp parallel for firstprivate(ix, jx, kx, gd) schedule(static)
        for (int k=1; k<=kx; k++) {
          for (int j=1; j<=jx; j++) {
            size_t m0 = _F_IDX_S3D(0,  j, k, ix, jx, kx, gd);
            size_t m1 = _F_IDX_S3D(ix, j, k, ix, jx, kx, gd);
            d_t[m0] = d_t[m1];
          }
        }
        break;
        
      case X_PLUS:
#pragma omp parallel for firstprivate(ix, jx, kx, gd) schedule(static)
        for (int k=1; k<=kx; k++) {
          for (int j=1; j<=jx; j++) {
            size_t m0 = _F_IDX_S3D(ix+1, j, k, ix, jx, kx, gd);
            size_t m1 = _F_IDX_S3D(1,    j, k, ix, jx, kx, gd);
            d_t[m0] = d_t[m1];
          }
        }
        break;
        
      case Y_MINUS:
#pragma omp parallel for firstprivate(ix, jx, kx, gd) schedule(static)
        for (int k=1; k<=kx; k++) {
          for (int i=1; i<=ix; i++) {
            size_t m0 = _F_IDX_S3D(i, 0,  k, ix, jx, kx, gd);
            size_t m1 = _F_IDX_S3D(i, jx, k, ix, jx, kx, gd);
            d_t[m0] = d_t[m1];
          }
        }
        break;
        
      case Y_PLUS:
#pragma omp parallel for firstprivate(ix, jx, kx, gd) schedule(static)
        for (int k=1; k<=kx; k++) {
          for (int i=1; i<=ix; i++) {
            size_t m0 = _F_IDX_S3D(i, jx+1, k, ix, jx, kx, gd);
            size_t m1 = _F_IDX_S3D(i, 1,    k, ix, jx, kx, gd);
            d_t[m0] = d_t[m1];
          }
        }
        break;
        
      case Z_MINUS:
#pragma omp parallel for firstprivate(ix, jx, kx, gd) schedule(static)
        for (int j=1; j<=jx; j++) {
          for (int i=1; i<=ix; i++) {
            size_t m0 = _F_IDX_S3D(i, j, 0,  ix, jx, kx, gd);
            size_t m1 = _F_IDX_S3D(i, j, kx, ix, jx, kx, gd);
            d_t[m0] = d_t[m1];
          }
        }
        break;
        
      case Z_PLUS:
#pragma omp parallel for firstprivate(ix, jx, kx, gd) schedule(static)
        for (int j=1; j<=jx; j++) {
          for (int i=1; i<=ix; i++) {
            size_t m0 = _F_IDX_S3D(i, j, kx+1, ix, jx, kx, gd);
            size_t m1 = _F_IDX_S3D(i, j, 1,    ix, jx, kx, gd);
            d_t[m0] = d_t[m1];
          }
        }
        break;
    }
  }
}


// #################################################################
/**
 * @brief 速度の外部周期境界条件（単純なコピー）
 * @param [in,out] d_v  速度ベクトル
 * @param [in]     face 面番号
 */
void SetBC3D::VobcPeriodicSimple(REAL_TYPE* d_v, const int face)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
  if ( numProc > 1 )
  {
    switch (face)
    {
      case X_MINUS:
        if ( paraMngr->PeriodicCommV3D(d_v, ix, jx, kx, gd, gd, X_DIR, PLUS2MINUS) != CPM_SUCCESS ) Exit(0);
        break;
        
      case X_PLUS:
        if ( paraMngr->PeriodicCommV3D(d_v, ix, jx, kx, gd, gd, X_DIR, MINUS2PLUS) != CPM_SUCCESS ) Exit(0);
        break;
        
      case Y_MINUS:
        if ( paraMngr->PeriodicCommV3D(d_v, ix, jx, kx, gd, gd, Y_DIR, PLUS2MINUS) != CPM_SUCCESS ) Exit(0);
        break;
        
      case Y_PLUS:
        if ( paraMngr->PeriodicCommV3D(d_v, ix, jx, kx, gd, gd, Y_DIR, MINUS2PLUS) != CPM_SUCCESS ) Exit(0);
        break;
        
      case Z_MINUS:
        if ( paraMngr->PeriodicCommV3D(d_v, ix, jx, kx, gd, gd, Z_DIR, PLUS2MINUS) != CPM_SUCCESS ) Exit(0);
        break;
        
      case Z_PLUS:
        if ( paraMngr->PeriodicCommV3D(d_v, ix, jx, kx, gd, gd, Z_DIR, MINUS2PLUS) != CPM_SUCCESS ) Exit(0);
        break;
    }
  }
  else // Serial
  {
    if ( nID[face] >= 0 ) return;
    
    switch (face) 
    {
      case X_MINUS:
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
        break;
        
      case X_PLUS:
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
        break;
        
      case Y_MINUS:
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
        break;
        
      case Y_PLUS:
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
        break;
        
      case Z_MINUS:
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
        break;
        
      case Z_PLUS:
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
        break;
    }
  }
  
}


// #################################################################
/**
 @brief 速度の内部周期境界条件（単純なコピー）
 @param d_v 速度ベクトル
 @param st コンポーネント範囲の開始インデクス
 @param ed コンポーネント範囲の終了インデクス
 @param d_bd BCindex ID
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
    case X_MINUS:
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
      
    case X_PLUS:
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
      
    case Y_MINUS:
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
      
    case Y_PLUS:
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
      
    case Z_MINUS:
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
      
    case Z_PLUS:
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
