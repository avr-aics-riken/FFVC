/*
 * SPHERE - Skeleton for PHysical and Engineering REsearch
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2004-2012
 *
 */

//@file CBC_SetBC.C
//@brief SetBC3D class
//@author keno, FSI Team, VCAD, RIKEN

#include "CBC_SetBC.h"

/**
 @fn void SetBC3D::assign_Temp(REAL_TYPE* t, unsigned* bh1, REAL_TYPE tm, Control* C)
 @brief 温度指定境界条件に必要な温度をセットする
 @param t 温度場
 @param bh BCindex H1
 @param tm 無次元時刻
 @param C
 */
void SetBC3D::assign_Temp(REAL_TYPE* t, unsigned* bh1, REAL_TYPE tm, Control* C)
{
  REAL_TYPE flop, tc;
  int st[3], ed[3];
  unsigned typ;
  
  // 内部境界条件による修正
  for (unsigned n=1; n<=NoBC; n++) {
    typ = cmp[n].getType();
    
    cmp[n].getBbox(st, ed);
    
    switch (typ) {
      case SPEC_VEL_WH:
        tc  = FBUtility::convK2ND(cmp[n].get_Temp(), BaseTemp, DiffTemp); // difference form BaseTemp 
        cbc_hbc_drchlt_(t, dim_sz, gc, st, ed, (int*)bh1, (int*)&n, &tc);
        break;
    }
  }
  /*
   // 外部境界条件
   int n = OBC_MASK;
   for (int face=0; face<NOFACE; face++) {
   typ = obc[face].get_BCtype();
   getOuterLoopIdx(face, st, ed);
   
   // 計算領域の最外郭領域でないときに，境界処理をスキップ，次のface面を評価
   if( pn.nID[face] >= 0 ) continue;
   
   if ( typ == OBC_SPEC_VEL ) {
   if ( !clear ) {
   extractVel_OBC(face, vec, tm, flop);
   }
   else {
   vec[0] = v00[1];
   vec[1] = v00[2];
   vec[2] = v00[3];
   }
   cbc_vbc_drchlt_cc_(v, dim_sz, gc, st, ed, v00, (int*)bv, &n, vec);
   }
   }*/
}

/**
 @fn void SetBC3D::assign_Velocity(REAL_TYPE* v, unsigned* bv, REAL_TYPE tm, REAL_TYPE* v00, bool clear)
 @brief 速度指定境界条件に必要な参照速度をセットする
 @param v セルセンタ速度ベクトル (n step)
 @param bv BCindex V
 @param tm 無次元時刻
 @param v00
 @param clear trueのとき，出力時に速度を壁面速度にする（デフォルトfalse）
 */
void SetBC3D::assign_Velocity(REAL_TYPE* v, unsigned* bv, REAL_TYPE tm, REAL_TYPE* v00, bool clear)
{
  REAL_TYPE flop, vec[3];
  int st[3], ed[3];
  unsigned typ;
  
  // 内部境界条件による修正
  for (int n=1; n<=NoBC; n++) {
    typ = cmp[n].getType();
    
    cmp[n].getBbox(st, ed);
    
    switch (typ) {
      case SPEC_VEL:
      case SPEC_VEL_WH:
        if ( !clear ) {
          extractVel_IBC(n, vec, tm, v00, flop);
        }
        else {
          vec[0] = vec[1] = vec[2] = 0.0;
        }
        cbc_vibc_drchlt_(v, dim_sz, gc, st, ed, v00, (int*)bv, &n, vec);
        break;
    }
  }
  
  // 外部境界条件
  for (int face=0; face<NOFACE; face++) {
    typ = obc[face].get_BCtype();
    getOuterLoopIdx(face, st, ed);
    
    // 計算領域の最外郭領域でないときに，境界処理をスキップ，次のface面を評価
    if( pn.nID[face] >= 0 ) continue;
    
    if ( typ == OBC_SPEC_VEL ) {
      if ( !clear ) {
        extractVel_OBC(face, vec, tm, v00, flop);
      }
      else {
        vec[0] = vec[1] = vec[2] = 0.0;
      }
      cbc_vobc_drchlt_(v, dim_sz, gc, v00, (int*)bv, &face, vec);
    }
  }
}

/**
 @fn void SetBC3D::checkDriver(FILE* fp)
 @brief ドライバ指定のチェック
 @param fp
 @note コンポーネントと外部境界で指定された，方向と位置の情報が一致するかをチェック
 */
void SetBC3D::checkDriver(FILE* fp)
{
  SklParaManager* para_mng = ParaCmpo->GetParaManager();
  
  int o_dir, c_dir;
  int o_pos, c_pos;
  int node_st_i, node_st_j, node_st_k;
  int* st = NULL;
  
  node_st_i = node_st_j = node_st_k = 0;
  
  if( para_mng->IsParallel() ){
    node_st_i = para_mng->GetVoxelHeadIndex(pn.ID, 0);
    node_st_j = para_mng->GetVoxelHeadIndex(pn.ID, 1);
    node_st_k = para_mng->GetVoxelHeadIndex(pn.ID, 2);
  }
  
  for (int face=0; face<NOFACE; face++) {
    if ( (obc[face].get_BCtype() == OBC_PERIODIC) && (obc[face].get_PrdcMode() == BoundaryOuter::prdc_Driver) ) {

      for (unsigned n=1; n<=NoBC; n++) {
        if ( cmp[n].getType() == PERIODIC ) {
          
          // 方向のチェック
          o_dir = obc[face].get_DriverDir();
          c_dir = (int)cmp[n].getPeriodicDir();
          if ( o_dir != c_dir ) {
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
          switch (c_dir) {
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
          
          if ( o_pos != c_pos ) {
            fprintf(fp, "\tDriver Lid Position is different between OBC[%d] and Component[%d].", o_pos, c_pos);
            printf("\tDriver Lid Position is different between OBC[%d] and Component[%d].", o_pos, c_pos);
            Exit(0);
          }
        }
      }
    }
  }  
}

/**
 @fn REAL_TYPE SetBC3D::extractVel_IBC(int n, REAL_TYPE* vec, REAL_TYPE tm, REAL_TYPE* v00, REAL_TYPE& flop)
 @brief コンポーネントから速度境界条件の成分を取り出す
 @param n コンポーネントのインデクス
 @param[out] vec[3] ベクトル成分
 @param tm 時刻
 @param v00
 @param[in/out] flop
 */
REAL_TYPE SetBC3D::extractVel_IBC(int n, REAL_TYPE* vec, REAL_TYPE tm, REAL_TYPE* v00, REAL_TYPE& flop)
{
  REAL_TYPE a, b, vel;
  const REAL_TYPE c_pai = (REAL_TYPE)(2.0*asin(1.0));
  
  a = cmp[n].ca[CompoList::amplitude]/RefV; // non-dimensional velocity amplitude
  b = 2.0*c_pai*cmp[n].ca[CompoList::frequency]* RefL/RefV * tm + cmp[n].ca[CompoList::initphase];
  vel = ( a*sin(b) + cmp[n].ca[CompoList::bias]/RefV ) * v00[0];
  vec[0] = cmp[n].nv[0] * vel;
  vec[1] = cmp[n].nv[1] * vel;
  vec[2] = cmp[n].nv[2] * vel;
  flop += 14.0;
  
  return vel;
}

/**
 @fn REAL_TYPE SetBC3D::extractVel_OBC(int n, REAL_TYPE* vec, REAL_TYPE tm, REAL_TYPE* v00, REAL_TYPE& flop)
 @brief 外部境界条件リストから速度境界条件の成分を取り出す
 @param n コンポーネントのインデクス
 @param[out] vec[3] ベクトル成分
 @param tm 時刻
 @param v00
 @param[in/out] flop
 */
REAL_TYPE SetBC3D::extractVel_OBC(int n, REAL_TYPE* vec, REAL_TYPE tm, REAL_TYPE* v00, REAL_TYPE& flop)
{
  REAL_TYPE a, b, vel;
  const REAL_TYPE c_pai = (REAL_TYPE)(2.0*asin(1.0));
  
  a = obc[n].ca[CompoList::amplitude]/RefV; // non-dimensional velocity amplitude
  b = 2.0*c_pai*obc[n].ca[CompoList::frequency]* RefL/RefV * tm + obc[n].ca[CompoList::initphase];
  vel = ( a*sin(b) + obc[n].ca[CompoList::bias]/RefV ) * v00[0];
  vec[0] = obc[n].nv[0] * vel;
  vec[1] = obc[n].nv[1] * vel;
  vec[2] = obc[n].nv[2] * vel;
  flop += 14.0;
  
  return vel;
}

/**
 @fn void SetBC3D::flipDir_OBC(unsigned* bv, Control* C)
 @brief 外部境界のIN_OUT境界条件の時のフラグを，流入出の方向に合わせてスイッチする
 @param bv BCindex V
 @param C
 */
void SetBC3D::flipDir_OBC(unsigned* bv, Control* C)
{
  unsigned F;
  
  for (int face=0; face<NOFACE; face++) {
    F = obc[face].get_BCtype();
    
    // 計算領域の最外郭領域でないときに，境界処理をスキップ，次のface面を評価
    if( pn.nID[face] >= 0 ) continue;
    
    if ( F == OBC_IN_OUT ) {
      
      // 流入・流出の方向を決定
      switch (face) {
        case X_MINUS:
        case Y_MINUS:
        case Z_MINUS:
          obc[face].Face_inout = ( C->V_Dface[face] >= 0.0 ) ? ALT_IN : ALT_OUT;
          break;
          
        case X_PLUS:
        case Y_PLUS:
        case Z_PLUS:
          obc[face].Face_inout = ( C->V_Dface[face] <  0.0 ) ? ALT_IN : ALT_OUT;
          break;
      }
      
      // OBC_MASKをセット
      flip_ObcMask(face, bv, obc[face].Face_inout);
    }
  }
}

/**
 @fn void SetBC3D::flip_ObcMask(int face, unsigned* bv, unsigned flag)
 @brief 外部境界に接するセルにおいて，bv[]にビットフラグをセットする
 @param face 外部境界面番号
 @param bv BCindex V
 @param flag 流入か流出のフラグ
 */
void SetBC3D::flip_ObcMask(int face, unsigned* bv, unsigned flag)
{
  int i, j, k;
  unsigned m;
  unsigned register s;
  
  switch (face) {
    case X_MINUS:
      if( pn.nID[X_MINUS] < 0 ){ // 外部境界をもつノードのみ
        for (k=1; k<=kxc; k++) {
          for (j=1; j<=jxc; j++) {
            m = FBUtility::getFindexS3D(size, guide, 1  , j, k);
            s = bv[m];
            
            if ( IS_FLUID(s) ) {
              if ( flag == ALT_OUT ) {
                s |=  (OBC_MASK << BC_FACE_W); // OBC_MASK==31 外部境界条件のフラグ
              }
              else {
                s &= ~(OBC_MASK << BC_FACE_W); // Wの5ビットをクリアする
              }
              bv[m]= s;
            }
          }
        }        
      }
      break;
      
    case X_PLUS:
      if( pn.nID[X_PLUS] < 0 ){
        for (k=1; k<=kxc; k++) {
          for (j=1; j<=jxc; j++) {
            m = FBUtility::getFindexS3D(size, guide, ixc  , j, k);
            s = bv[m];
            
            if ( IS_FLUID(s) ) {
              if ( flag == ALT_OUT ) {
                s |= (OBC_MASK << BC_FACE_E);
              }
              else {
                s &= ~(OBC_MASK << BC_FACE_E);
              }
              bv[m]= s;
            }
          }
        }
      }
      break;
      
    case Y_MINUS:
      if( pn.nID[Y_MINUS] < 0 ){
        for (k=1; k<=kxc; k++) {
          for (i=1; i<=ixc; i++) {
            m = FBUtility::getFindexS3D(size, guide, i, 1  , k);
            s = bv[m];
            
            if ( IS_FLUID(s) ) {
              if ( flag == ALT_OUT ) {
                s |= (OBC_MASK << BC_FACE_S);
              }
              else {
                s &= ~(OBC_MASK << BC_FACE_S);
              }
              bv[m]= s;
            }
          }
        }
      }
      break;
      
    case Y_PLUS:
      if( pn.nID[Y_PLUS] < 0 ){
        for (k=1; k<=kxc; k++) {
          for (i=1; i<=ixc; i++) {
            m = FBUtility::getFindexS3D(size, guide, i, jxc  , k);
            s = bv[m];
            
            if ( IS_FLUID(s) ) {
              if ( flag == ALT_OUT ) {
                s |= (OBC_MASK << BC_FACE_N);
              }
              else {
                s &= ~(OBC_MASK << BC_FACE_N);
              }
              bv[m]= s;
            }
          }
        }
      }
      break;
      
    case Z_MINUS:
      if( pn.nID[Z_MINUS] < 0 ){
        for (j=1; j<=jxc; j++) {
          for (i=1; i<=ixc; i++) {
            m = FBUtility::getFindexS3D(size, guide, i, j, 1  );
            s = bv[m];
            
            if ( IS_FLUID(s) ) {
              if ( flag == ALT_OUT ) {
                s |= (OBC_MASK << BC_FACE_B);
              }
              else {
                s &= ~(OBC_MASK << BC_FACE_B);
              }
              bv[m]= s;
            }
          }
        }
      }
      break;
      
    case Z_PLUS:
      if( pn.nID[Z_PLUS] < 0 ){
        for (j=1; j<=jxc; j++) {
          for (i=1; i<=ixc; i++) {
            m = FBUtility::getFindexS3D(size, guide, i, j, kxc  );
            s = bv[m];
            
            if ( IS_FLUID(s) ) {
              if ( flag == ALT_OUT ) {
                s |= (OBC_MASK << BC_FACE_T);
              }
              else {
                s &= ~(OBC_MASK << BC_FACE_T);
              }
              bv[m]= s;
            }
          }
        }
      }
      break;
  } // end of switch
  
}

/**
 @fn void SetBC3D::InnerTBCvol(REAL_TYPE* t, unsigned* bh2, REAL_TYPE dt, REAL_TYPE& flop)
 @brief セルに対する温度の内部境界をセットする
 @param t 温度
 @param bh2 BCindex H2
 @param dt 時間積分幅
 @param flop 浮動小数演算数
 */
void SetBC3D::InnerTBCvol(REAL_TYPE* t, unsigned* bh2, REAL_TYPE dt, REAL_TYPE& flop)
{
  for (unsigned n=1; n<=NoBC; n++) {
    
    switch ( cmp[n].getType() ) {
      case HEAT_SRC:
        cmp[n].set_Mon_Calorie( ps_IBC_HeatGen_SM(t, bh2, n, dt, flop) );
        break;
        
      case CNST_TEMP:
        ps_IBC_ConstTemp(t, bh2, n);
        break;
    }
  }
}

/**
 @fn void SetBC3D::InnerVBC_Periodic(SklVector3DEx<REAL_TYPE>* d_v, SklScalar3D<unsigned>* d_bd)
 @brief 速度ベクトルの内部周期境界条件処理
 @param d_v 速度ベクトルのデータクラス
 @param d_bd BCindex ID
 */
void SetBC3D::InnerVBC_Periodic(SklVector3DEx<REAL_TYPE>* d_v, SklScalar3D<unsigned>* d_bd)
{
  REAL_TYPE *v=NULL;
  if ( !(v = d_v->GetData()) ) Exit(0);
  
  int st[3], ed[3];
  
  for (int n=1; n<=NoBC; n++) {
    cmp[n].getBbox(st, ed);
    
    if ( cmp[n].getType() == PERIODIC ) {
      mark();
      Vibc_Prdc(d_v, st, ed, d_bd, n, (int)cmp[n].getPeriodicDir());
    }    
  }
}

/**
 @fn void SetBC3D::InnerVBC(REAL_TYPE* v, unsigned* bv, SklScalar3D<unsigned>* d_bd, REAL_TYPE tm, REAL_TYPE* v00, REAL_TYPE& flop)
 @brief 速度ベクトルの内部境界条件処理(タイムステップに一度)
 @param v 速度ベクトル
 @param bv BCindex V
 @param tm
 @param v00
 @param flop
 */
void SetBC3D::InnerVBC(REAL_TYPE* v, unsigned* bv, REAL_TYPE tm, REAL_TYPE* v00, REAL_TYPE& flop)
{
  int st[3], ed[3];
  REAL_TYPE vec[3];
  
  if ( isCDS ) { // Cut-Distance
    for (int n=1; n<=NoBC; n++) {
      cmp[n].getBbox(st, ed);
      
      switch ( cmp[n].getType() ) {
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
  else { // Binary
    for (int n=1; n<=NoBC; n++) {
      cmp[n].getBbox(st, ed);
      
      switch ( cmp[n].getType() ) {
        case SPEC_VEL:
        case SPEC_VEL_WH:
          extractVel_IBC(n, vec, tm, v00, flop);
          cbc_vibc_drchlt_(v, dim_sz, gc, st, ed, v00, (int*)bv, &n, vec);
          break;
          
        case OUTFLOW:
          cbc_vibc_outflow_(v, dim_sz, gc, st, ed, (int*)bv, &n);
          break;
          
        default:
          break;
      }  
    }
  }
  
  
}

/**
 @fn void SetBC3D::InnerTBCface(REAL_TYPE* qbc, unsigned* bh1, REAL_TYPE* t, REAL_TYPE* t0, REAL_TYPE& flop)
 @brief 拡散部分に関する温度の内部境界処理
 @param qbc 境界条件熱流束
 @param bh1 BCindex H1
 @param t n+1時刻の温度場
 @param t0 温度
 @param flop 浮動小数演算数
 @note 内部境界の断熱BCは断熱マスクで処理
 */
void SetBC3D::InnerTBCface(REAL_TYPE* qbc, unsigned* bh1, REAL_TYPE* t, REAL_TYPE* t0, REAL_TYPE& flop)
{
  unsigned F, H;
  
  for (unsigned n=1; n<=NoBC; n++) {
    F = cmp[n].getType();
    
    switch ( F ) {
      case HEATFLUX:
        cmp[n].set_Mon_Heatflux( ps_IBC_Heatflux(qbc, bh1, n, flop) );
        break;
        
      case TRANSFER:
        H = cmp[n].getHtype();
        //if      ( H == HT_N)  cmp[n].cmp[n].set_Mon_Calorie( setHeatTransferN_SM(qbc, t, bx, n, t0, flop) );
        if      ( H == HT_S)  cmp[n].set_Mon_Calorie( ps_IBC_Transfer_S_SM (qbc, bh1, n, t, t0, flop) );
        else if ( H == HT_SN) cmp[n].set_Mon_Calorie( ps_IBC_Transfer_SN_SM(qbc, bh1, n, t, t0, flop) );
        else if ( H == HT_SF) cmp[n].set_Mon_Calorie( ps_IBC_Transfer_SF_SM(qbc, bh1, n, t, t0, flop) );
        else if ( H == HT_B)  cmp[n].set_Mon_Calorie( ps_IBC_Transfer_B_SM (qbc, bh1, n, t, t0, flop) );
        break;
        
      case ISOTHERMAL:
        cmp[n].set_Mon_Calorie( ps_IBC_IsoThermal_SM(qbc, bh1, n, t, t0, flop) );
        break;
        
      case RADIANT:
        //setRadiant(qbc, bx, n, t0);
        break;
    }    
  }
}

/**
 @fn void SetBC3D::InnerPBC_Periodic(SklScalar3D<REAL_TYPE>* d_p, SklScalar3D<unsigned>* d_bcd)
 @brief 圧力の内部境界条件処理
 @param d_p 圧力のデータクラス
 @param d_bcd BCindex ID
 */
void SetBC3D::InnerPBC_Periodic(SklScalar3D<REAL_TYPE>* d_p, SklScalar3D<unsigned>* d_bcd)
{
  int dir;
  int st[3], ed[3];
  REAL_TYPE pv;
  
  for (unsigned n=1; n<=NoBC; n++) {
    cmp[n].getBbox(st, ed);
    dir = (int)cmp[n].getPeriodicDir();
    pv = FBUtility::convD2ND_P(cmp[n].ca[0], BasePrs, rho, RefV, Unit_Prs);
    
    if ( cmp[n].getType() == PERIODIC ) {
      Pibc_Prdc(d_p, st, ed, d_bcd, n, dir, pv);
    }
  }
}

/**
 @fn void SetBC3D::mod_div(REAL_TYPE* div, unsigned* bv, REAL_TYPE coef, REAL_TYPE tm, REAL_TYPE* v00, REAL_TYPE* avr, REAL_TYPE& flop)
 @brief 速度境界条件による速度の発散の修正ほか
 @param div div((u)*(-h/dt)
 @param bv BCindex V
 @param coef 係数 h/dt
 @param tm 無次元時刻
 @param v00
 @param avr 平均値計算のテンポラリ値
 @param flop
 @note 外部境界面のdiv(u)の修正時に領域境界の流量などのモニタ値を計算し，BoundaryOuterクラスに保持 > 反復後にDomainMonitor()で集約
 */
void SetBC3D::mod_div(REAL_TYPE* div, unsigned* bv, REAL_TYPE coef, REAL_TYPE tm, REAL_TYPE* v00, REAL_TYPE* avr, REAL_TYPE& flop)
{
  REAL_TYPE vec[3], dummy;
  int st[3], ed[3];
  unsigned typ=0;
  
  // 内部境界条件による修正
  if ( isCDS ) { // Cut-Distance
    for (int n=1; n<=NoBC; n++) {
      typ = cmp[n].getType();
      
      cmp[n].getBbox(st, ed);
      
      switch (typ) {
        case OUTFLOW:
          //cbc_div_ibc_oflow_vec_(div, dim_sz, gc, st, ed, v00, &coef, (int*)bv, &n, avr, &flop);
          break;
          
        case SPEC_VEL:
        case SPEC_VEL_WH:
          cmp[n].val[var_Velocity] = extractVel_IBC(n, vec, tm, v00, flop); // 指定された無次元平均流速
          cbc_div_ibc_drchlt_(div, dim_sz, gc, st, ed, v00, &coef, (int*)bv, &n, vec, &flop);
          break;
          
        default:
          break;
      }
    }
  }
  else { // Binary
    for (int n=1; n<=NoBC; n++) {
      typ = cmp[n].getType();
      
      cmp[n].getBbox(st, ed);
      
      switch (typ) {
        case OUTFLOW:
          cbc_div_ibc_oflow_vec_(div, dim_sz, gc, st, ed, v00, &coef, (int*)bv, &n, avr, &flop);
          break;
          
        case SPEC_VEL:
        case SPEC_VEL_WH:
          cmp[n].val[var_Velocity] = extractVel_IBC(n, vec, tm, v00, flop); // 指定された無次元平均流速
          cbc_div_ibc_drchlt_(div, dim_sz, gc, st, ed, v00, &coef, (int*)bv, &n, vec, &flop);
          break;
          
        default:
          break;
      }
    }
  }
  
  
  // 外部境界条件による修正
  for (int face=0; face<NOFACE; face++) {
    typ = obc[face].get_BCtype();
    
    // 計算領域の最外郭領域でないときに，境界処理をスキップ，次のface面を評価
    if( pn.nID[face] >= 0 ) {
      vec[0] = 0.0;   // sum
      vec[1] = 1.0e6; // min
      vec[2] =-1.0e6; // max
      obc[face].set_DomainV(vec, face, true); // gatherする場合のダミー値を与えておく
      continue;
    }
    
    switch (typ) {
      case OBC_OUTFLOW:
      case OBC_TRC_FREE:
      case OBC_IN_OUT:
        cbc_div_obc_oflow_vec_(div, dim_sz, gc, &face, v00, &coef, (int*)bv, vec, &flop); // vecは流用
        obc[face].set_DomainV(vec, face, true); // vecは速度の和の形式で保持
        break;
        
      case OBC_SPEC_VEL:
      case OBC_WALL:
        dummy = extractVel_OBC(face, vec, tm, v00, flop);
        cbc_div_obc_drchlt_(div, dim_sz, gc, &face, v00, &coef, (int*)bv, vec, &flop);
        obc[face].set_DomainV(vec, face); // 速度の形式
        break;
        
      case OBC_SYMMETRIC:
        vec[0] = vec[1] = vec[2] = 0.0;
        // fluxはゼロなので処理不要　cbc_div_obc_drchlt_(div, dim_sz, gc, &face, v00, &coef, (int*)bv, vec, &flop);
        obc[face].set_DomainV(vec, face); // 速度の形式
        break;
    }
  }
}

/**
 @fn void SetBC3D::mod_Dir_Forcing(REAL_TYPE* v, unsigned* bd, float* cvf, REAL_TYPE* v00, REAL_TYPE &flop)
 @brief 圧力損失部による速度の方向修正
 @param[in/out] v 速度
 @param bd BCindex ID
 @param cvf コンポーネントの体積率
 @param v00 参照速度
 @param[out] flop
 */
void SetBC3D::mod_Dir_Forcing(REAL_TYPE* v, unsigned* bd, float* cvf, REAL_TYPE* v00, REAL_TYPE &flop)
{
  int st[3], ed[3];
  REAL_TYPE vec[3];
  
  for (unsigned n=1; n<=NoBC; n++) {
    if ( cmp[n].isFORCING() ) {
      
      cmp[n].getBbox(st, ed);
      
      vec[0] = cmp[n].nv[0];
      vec[1] = cmp[n].nv[1];
      vec[2] = cmp[n].nv[2];
      
      switch ( cmp[n].getType() ) {
        case HEX:
          if ( cmp[n].get_sw_HexDir() ) cbc_hex_dir_ (v, dim_sz, gc, st, ed, (int*)bd, cvf, (int*)&n, v00, vec, &flop);
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

/**
 @fn void SetBC3D::mod_Pvec_Forcing(REAL_TYPE* vc, REAL_TYPE* v, unsigned* bd, float* cvf, REAL_TYPE* v00, REAL_TYPE dt, REAL_TYPE &flop)
 @brief 圧力損失部による疑似速度方向の修正
 @param[in/out] vc 疑似速度ベクトル
 @param[in] v 速度ベクトル n-step
 @param bd BCindex ID
 @param cvf コンポーネントの体積率
 @param v00 参照速度
 @param dt 時間積分幅
 @param[in/out] flop
 */
void SetBC3D::mod_Pvec_Forcing(REAL_TYPE* vc, REAL_TYPE* v, unsigned* bd, float* cvf, REAL_TYPE* v00, REAL_TYPE dt, REAL_TYPE &flop)
{
  int st[3], ed[3];
  REAL_TYPE vec[3];
  
  for (unsigned n=1; n<=NoBC; n++) {
    vec[0] = cmp[n].nv[0];
    vec[1] = cmp[n].nv[1];
    vec[2] = cmp[n].nv[2];
    
    cmp[n].getBbox(st, ed);
    
    switch ( cmp[n].getType() ) {
      case HEX:
        cbc_hex_force_pvec_(vc, dim_sz, gc, st, ed, (int*)bd, cvf, v, (int*)&n, v00, &dt, vec, &cmp[n].ca[0], &flop);
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

/**
 @fn void SetBC3D::mod_Psrc_Forcing(REAL_TYPE* src, REAL_TYPE* v, unsigned* bd, float* cvf, REAL_TYPE dh, REAL_TYPE* v00, REAL_TYPE** c_array, REAL_TYPE &flop)
 @brief 圧力損失部によるPoisosn式のソース項 \gamma^F の修正とワーク用の速度を保持
 @param[out] src 外力項によるPoisson方程式のソース項
 @param v 速度ベクトル n+1
 @param bd BCindex ID
 @param cvf コンポーネントの体積率
 @param dh 格子幅
 @param v00 参照速度
 @param c_array コンポーネントワーク配列の管理ポインタ
 @param[out] flop
 */
void SetBC3D::mod_Psrc_Forcing(REAL_TYPE* src, REAL_TYPE* v, unsigned* bd, float* cvf, REAL_TYPE dh, REAL_TYPE* v00, REAL_TYPE** c_array, REAL_TYPE &flop)
{
  int st[3], ed[3], csz[3];
  REAL_TYPE vec[3];
  REAL_TYPE* w_ptr=NULL;
  
  for (unsigned n=1; n<=NoBC; n++) {
    vec[0] = cmp[n].nv[0];
    vec[1] = cmp[n].nv[1];
    vec[2] = cmp[n].nv[2];
    
    cmp[n].getBbox(st, ed);
    cmp[n].get_cmp_sz(csz);
    w_ptr = c_array[n];
    
    if ( cmp[n].isFORCING() ) cbc_force_keep_vec_(w_ptr, csz, st, ed, v, dim_sz, gc);

    switch ( cmp[n].getType() ) {
      case HEX:
        cbc_hex_psrc_(src, dim_sz, gc, st, ed, (int*)bd, cvf, w_ptr, csz, (int*)&n, v00, &dh, vec, &cmp[n].ca[0], &flop);
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

/**
 @fn void SetBC3D::mod_Vdiv_Forcing(REAL_TYPE* v, unsigned* bd, float* cvf, REAL_TYPE* div, REAL_TYPE dt, REAL_TYPE dh, REAL_TYPE* v00, REAL_TYPE* am, REAL_TYPE** c_array, REAL_TYPE &flop)
 @brief 圧力損失部によるセルセンタ速度の修正と速度の発散値の修正
 @param[in/out] v セルセンターの速度
 @param bd BCindex ID
 @param cvf コンポーネントの体積率
 @param div div((u)*(dh/dt)
 @param dt 時間積分幅
 @param dh 格子幅
 @param v00 参照速度
 @param am モニター用配列
 @param c_array コンポーネントワーク配列の管理ポインタ
 @param[out] flop
 */
void SetBC3D::mod_Vdiv_Forcing(REAL_TYPE* v, unsigned* bd, float* cvf, REAL_TYPE* div, REAL_TYPE dt, REAL_TYPE dh, REAL_TYPE* v00, REAL_TYPE* am, REAL_TYPE** c_array, REAL_TYPE &flop)
{
  int st[3], ed[3], csz[3];
  REAL_TYPE vec[3];
  REAL_TYPE* w_ptr=NULL;
  REAL_TYPE coef = dh/dt;
  
  for (unsigned n=1; n<=NoBC; n++) {
    if ( cmp[n].isFORCING() ) {
      
      cmp[n].getBbox(st, ed);
      cmp[n].get_cmp_sz(csz);
      w_ptr = c_array[n];
      
      vec[0] = cmp[n].nv[0];
      vec[1] = cmp[n].nv[1];
      vec[2] = cmp[n].nv[2];
      
      switch ( cmp[n].getType() ) {
        case HEX:
          cbc_hex_force_vec_(v, div, dim_sz, gc, st, ed, (int*)bd, cvf, w_ptr, csz, (int*)&n, v00, &dt, &dh, vec, &cmp[n].ca[0], am, &flop);
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

/**
 @fn void SetBC3D::mod_Pvec_Flux(REAL_TYPE* wv, SKL_RAEL* v, unsigned* bv, REAL_TYPE tm, Control* C, int v_mode, REAL_TYPE* v00, REAL_TYPE& flop)
 @brief 速度境界条件による流束の修正
 @param[in/out] wv 疑似ベクトル
 @param v 速度ベクトル u^n
 @param bv BCindex V
 @param tm 無次元時刻
 @param C
 @param v_mode 粘性項のモード (0=粘性項を計算しない, 1=粘性項を計算する, 2=壁法則)
 @param v00
 @param[out] flop
 */
void SetBC3D::mod_Pvec_Flux(REAL_TYPE* wv, REAL_TYPE* v, unsigned* bv, REAL_TYPE tm, Control* C, int v_mode, REAL_TYPE* v00, REAL_TYPE& flop)
{
  REAL_TYPE vec[3];
  int st[3], ed[3], order;
  unsigned typ;
  
  // 内部境界（流束形式）
  if ( isCDS ) { // Cut-Distance
    for (int n=1; n<=NoBC; n++) {
      
      if( cmp[n].isEns() ) {
        typ = cmp[n].getType();
        cmp[n].getBbox(st, ed);

        if ( (typ==SPEC_VEL) || (typ==SPEC_VEL_WH) ) {
          extractVel_IBC(n, vec, tm, v00, flop);
          cds_pvec_vibc_specv_(wv, dim_sz, gc, st, ed, &dh, v00, &rei, v, (int*)bv, &n, vec, &flop);
          cds_pvec_vibc_specv2_(wv, dim_sz, gc, st, ed, &dh, v00, (int*)bv, &n, vec, &flop);
        }
        else if ( typ==OUTFLOW ) {
        }      
      }
    }
  }
  else { // Binary
    for (int n=1; n<=NoBC; n++) {
      
      if( cmp[n].isEns() ) {
        typ = cmp[n].getType();
        cmp[n].getBbox(st, ed);
        
        if ( (typ==SPEC_VEL) || (typ==SPEC_VEL_WH) ) {
          extractVel_IBC(n, vec, tm, v00, flop);
          cbc_pvec_vibc_specv_(wv, dim_sz, gc, st, ed, &dh, v00, &rei, v, (int*)bv, &n, vec, &v_mode, &flop);
        }
        else if ( typ==OUTFLOW ) {
          vec[0] = vec[1] = vec[2] = cmp[n].val[var_Velocity]; // mod_div()でセルフェイス流出速度がval[var_Velocity]にセット
          cbc_pvec_vibc_oflow_(wv, dim_sz, gc, st, ed, &dh, v00, &rei, v, (int*)bv, &n, vec, &v_mode, &flop);
        }      
      }
    }
  }
  
  
  // 流束形式の外部境界条件
  for (int face=0; face<NOFACE; face++) {
    typ = obc[face].get_BCtype();
    
    // 計算領域の最外郭領域でないときに，境界処理をスキップ，次のface面を評価
    if( pn.nID[face] >= 0 ) continue;
    
    switch ( typ ) {
      case OBC_SPEC_VEL:
        extractVel_OBC(face, vec, tm, v00, flop);
        cbc_pvec_vobc_specv_(wv, dim_sz, gc, &dh, v00, &rei, v, (int*)bv, vec, &v_mode, &face, &flop);
        break;
        
      case OBC_WALL:
        extractVel_OBC(face, vec, tm, v00, flop);
        cbc_pvec_vobc_wall_(wv, dim_sz, gc, &dh, v00, &rei, v, (int*)bv, vec, &v_mode, &face, &flop);
        break;
        
      case OBC_SYMMETRIC:
        cbc_pvec_vobc_symtrc_(wv, dim_sz, gc, &dh, &rei, v, (int*)bv, &v_mode, &face, &flop);
        break;
        
      case OBC_OUTFLOW:
        vec[0] = vec[1] = vec[2] = C->V_Dface[face];
        cbc_pvec_vobc_oflow_(wv, dim_sz, gc, &dh, v00, &rei, v, (int*)bv, vec, &v_mode, &face, &flop);
        break;
        
      case OBC_IN_OUT:
        if ( obc[face].Face_inout == ALT_OUT ) { // 流出のときのみOUTFLOWと同様の処理
          vec[0] = vec[1] = vec[2] = C->V_Dface[face];
          cbc_pvec_vobc_oflow_(wv, dim_sz, gc, &dh, v00, &rei, v, (int*)bv, vec, &v_mode, &face, &flop);
        }
        break;
    }
  }
  
}

/**
 @fn void SetBC3D::mod_Psrc_VBC(REAL_TYPE* div, REAL_TYPE* vc, REAL_TYPE* v0, REAL_TYPE coef, unsigned* bv, REAL_TYPE tm, REAL_TYPE dt, Control* C, REAL_TYPE* v00, REAL_TYPE &flop)
 @brief 速度境界条件によるPoisosn式のソース項の修正
 @param[out] div divergence field
 @param vc セルセンタ疑似速度
 @param v0 セルセンタ速度 u^n
 @param coef 係数
 @param bv BCindex V
 @param tm 
 @param dt
 @param C
 @param v00
 @param[out] flop
 */
void SetBC3D::mod_Psrc_VBC(REAL_TYPE* div, REAL_TYPE* vc, REAL_TYPE* v0, REAL_TYPE coef, unsigned* bv, REAL_TYPE tm, REAL_TYPE dt, Control* C, REAL_TYPE* v00, REAL_TYPE &flop)
{
  int st[3], ed[3];
  REAL_TYPE vec[3], dummy, vel;
  unsigned typ;
  
  // 内部境界条件による修正
  if ( isCDS ) { // Cut-Distance
    for (int n=1; n<=NoBC; n++) {
      typ = cmp[n].getType();
      cmp[n].getBbox(st, ed);

      switch (typ) {
        case SPEC_VEL:
        case SPEC_VEL_WH:
          dummy = extractVel_IBC(n, vec, tm, v00, flop);
          cbc_div_ibc_drchlt_(div, dim_sz, gc, st, ed, v00, &coef, (int*)bv, &n, vec, &flop);
          break;
          
        case OUTFLOW:
          break;
          
        default:
          break;
      }
    }
  }
  else { // Binary
    for (int n=1; n<=NoBC; n++) {
      typ = cmp[n].getType();
      cmp[n].getBbox(st, ed);
      
      switch (typ) {
        case SPEC_VEL:
        case SPEC_VEL_WH:
          dummy = extractVel_IBC(n, vec, tm, v00, flop);
          cbc_div_ibc_drchlt_(div, dim_sz, gc, st, ed, v00, &coef, (int*)bv, &n, vec, &flop);
          break;
          
        case OUTFLOW:
          vel = cmp[n].val[var_Velocity] * dt / dh; // mod_div()でval[var_Velocity]にセット
          cbc_div_ibc_oflow_pvec_(div, dim_sz, gc, st, ed, v00, &vel, &coef, (int*)bv, &n, v0, &flop);
          break;
          
        default:
          break;
      }
    }
  }
  
  
  // 外部境界条件による修正
  for (int face=0; face<NOFACE; face++) {
    typ = obc[face].get_BCtype();
    
    // 計算領域の最外郭領域でないときに，境界処理をスキップ，次のface面を評価
    if( pn.nID[face] >= 0 ) continue;
    
    switch ( typ ) {
      case OBC_SPEC_VEL:
      case OBC_WALL:
        dummy = extractVel_OBC(face, vec, tm, v00, flop);
        cbc_div_obc_drchlt_(div, dim_sz, gc, &face, v00, &coef, (int*)bv, vec, &flop);
        break;
        
      case OBC_SYMMETRIC:
        // 境界面の法線速度はゼロなので，修正不要 移動格子の場合は必要
        //vec[0] = vec[1] = vec[2] = 0.0;
        //cbc_div_obc_drchlt_(div, dim_sz, gc, &face, v00, &coef, (int*)bv, vec, &flop);
        break;
        
      case OBC_OUTFLOW:
        vel = C->V_Dface[face] * dt / dh;
        cbc_div_obc_oflow_pvec_(div, dim_sz, gc, &face, v00, &vel, &coef, (int*)bv, v0, &flop);
        break;
        
      case OBC_IN_OUT:
        if ( obc[face].Face_inout == ALT_OUT ) { // 流出のときのみOUTFLOWと同様の処理
          vel = C->V_Dface[face] * dt / dh;
          cbc_div_obc_oflow_pvec_(div, dim_sz, gc, &face, v00, &vel, &coef, (int*)bv, v0, &flop);
        }
        break;
    }
  }
}

/**
 @fn void SetBC3D::mod_Vis_EE(REAL_TYPE* vc, REAL_TYPE* v0, REAL_TYPE cf, unsigned* bx, REAL_TYPE tm, REAL_TYPE dt, REAL_TYPE* v00, REAL_TYPE& flop)
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
void SetBC3D::mod_Vis_EE(REAL_TYPE* vc, REAL_TYPE* v0, REAL_TYPE cf, unsigned* bx, REAL_TYPE tm, REAL_TYPE dt, REAL_TYPE* v00, REAL_TYPE& flop)
{
  int st[3], ed[3];
  REAL_TYPE vec[3], dummy;
  unsigned typ;
  
  for (int n=1; n<=NoBC; n++) {
    typ = cmp[n].getType();
    if ( (typ==SPEC_VEL) || (typ==SPEC_VEL_WH) ) {
      dummy = extractVel_IBC(n, vec, tm, v00, flop);
    }
    else if ( typ==OUTFLOW ) {
      vec[0] = vec[1] = vec[2] = cmp[n].val[var_Velocity];
    }
    
    cmp[n].getBbox(st, ed);
    cbc_vis_ee_vbc_(vc, dim_sz, gc, st, ed, &dh, &dt, v00, &rei, v0, (int*)bx, &n, &cf, vec, &flop);
  }
}

/**
 @fn void SetBC3D::OuterPBC(SklScalar3D<REAL_TYPE>* d_p)
 @brief 圧力の外部境界条件
 @param d_p 圧力のデータクラス
 */
void SetBC3D::OuterPBC(SklScalar3D<REAL_TYPE>* d_p)
{
  unsigned uod, F;
  REAL_TYPE pv=0.0;
  REAL_TYPE *p=NULL;
  
  if ( !(p = d_p->GetData()) ) Exit(0);
  
  for (int face=0; face<NOFACE; face++) {
    F = obc[face].get_BCtype();
    
    // 周期境界条件
    if ( F == OBC_PERIODIC ) {
      pv = FBUtility::convD2ND_P(obc[face].p, BasePrs, rho, RefV, Unit_Prs);

      switch ( obc[face].get_PrdcMode() ) {
        case BoundaryOuter::prdc_Simple:
          Pobc_Prdc_Simple(d_p, face);
          break;
          
        case BoundaryOuter::prdc_Directional:
          uod = obc[face].get_FaceMode();
          Pobc_Prdc_Directional(d_p, face, pv, uod);
          break;
          
        case BoundaryOuter::prdc_Driver:
          // nothing
          break;
      }
    }
    else { // 周期境界条件以外の処理
      if( pn.nID[face] >= 0 ) continue; // @note 並列時，計算領域の最外郭領域でないときに，境界処理をスキップ，次のface面を評価
      
      if ( F == OBC_IN_OUT ) {
        if ( obc[face].Face_inout == ALT_OUT ) {
          cbc_pobc_neumann_(p, dim_sz, gc, &face);
        }
        else {
          cbc_pobc_drchlt_ (p, dim_sz, gc, &face, &pv);
        }
      }
    }
  }
}

/**
 @fn void SetBC3D::OuterVBC_Periodic(SklVector3DEx<REAL_TYPE>* d_v)
 @brief 速度の外部境界条件処理
 @param d_v 速度ベクトルのデータクラス
 */
void SetBC3D::OuterVBC_Periodic(SklVector3DEx<REAL_TYPE>* d_v)
{
  REAL_TYPE *v=NULL;
  if ( !(v = d_v->GetData()) ) Exit(0);
  
  for (int face=0; face<NOFACE; face++) {
    
    if ( obc[face].get_BCtype() == OBC_PERIODIC ) {
      mark();
      unsigned pm = obc[face].get_PrdcMode();
      if ( (pm == BoundaryOuter::prdc_Simple) || (pm == BoundaryOuter::prdc_Directional)) { // BoundaryOuter::prdc_Driverに対しては処理不要
        Vobc_Prdc(d_v, face, guide); // セルフェイスの値の周期処理は不要
      }
    }
  }  
}

/**
 @fn void SetBC3D::OuterVBC(REAL_TYPE* v, REAL_TYPE* vc, unsigned* bv, REAL_TYPE tm, REAL_TYPE dt, Control* C, REAL_TYPE* v00, REAL_TYPE& flop)
 @brief 速度の外部境界条件処理(タイムステップに一度)
 @param[out] v 速度ベクトル v^{n+1}
 @param vc 速度ベクトル v^*
 @param bv BCindex V
 @param tm
 @param dt 
 @param C
 @param v00
 @param flop
 */
void SetBC3D::OuterVBC(REAL_TYPE* v, REAL_TYPE* vc, unsigned* bv, REAL_TYPE tm, REAL_TYPE dt, Control* C, REAL_TYPE* v00, REAL_TYPE& flop)
{
  REAL_TYPE vec[3];
  REAL_TYPE v_cnv;
  
  for (int face=0; face<NOFACE; face++) {

    if( pn.nID[face] >= 0 ) continue; // @note 並列時，計算領域の最外郭領域でないときに，境界処理をスキップ，次のface面を評価
    // @note ここでスキップする場合には，tfreeの処理でMPI通信をしないこと（参加しないノードがあるためエラーとなる）
    
    switch ( obc[face].get_BCtype() ) {
      case OBC_OUTFLOW:
        cbc_vobc_update_(v, dim_sz, gc, vc, &face);
        break;
        
      case OBC_SPEC_VEL:
        extractVel_OBC(face, vec, tm, v00, flop);
        cbc_vobc_drchlt_(v, dim_sz, gc, v00, (int*)bv, &face, vec);
        break;
        
      case OBC_TRC_FREE:
        //cbc_vobc_tfree_(v, dim_sz, gc, &face, &flop);
        cbc_vobc_neumann_(v, dim_sz, gc, &face);
        break;
        
      case OBC_IN_OUT:
        if ( obc[face].Face_inout == ALT_IN ) {
          cbc_vobc_tfree_(v, dim_sz, gc, &face, &flop);
        }
        break;
        
      default:
        break;
    }
    
  }  
}

/**
 @fn void SetBC3D::OuterVBC_Pseudo(REAL_TYPE* vc, REAL_TYPE* v0, unsigned* bv, REAL_TYPE tm, REAL_TYPE dt, Control* C, REAL_TYPE* v00, REAL_TYPE& flop)
 @brief 疑似速度の外部境界条件処理
 @param[out] vc 疑似速度ベクトル v^*
 @param v0 速度ベクトル v^n
 @param bv BCindex V
 @param tm
 @param dt 
 @param C
 @param v00
 @param flop
 */
void SetBC3D::OuterVBC_Pseudo(REAL_TYPE* vc, REAL_TYPE* v0, unsigned* bv, REAL_TYPE tm, REAL_TYPE dt, Control* C, REAL_TYPE* v00, REAL_TYPE& flop)
{
  REAL_TYPE v_cnv;
  
  for (int face=0; face<NOFACE; face++) {
    
    if( pn.nID[face] >= 0 ) continue; // @note 並列時，計算領域の最外郭領域でないときに，境界処理をスキップ，次のface面を評価
    // @note ここでスキップする場合には，MPI通信をしないこと（参加しないノードがあるためエラーとなる）
    
    switch ( obc[face].get_BCtype() ) {
      case OBC_OUTFLOW:
        v_cnv = C->V_Dface[face] * dt / dh;
        cbc_vobc_outflow_(vc, dim_sz, gc, &v_cnv, (int*)bv, &face, v0, &flop);
        break;
    }
    
  }  
}

/**
 @fn void SetBC3D::OuterTBC(SklScalar3D<REAL_TYPE>* d_t)
 @brief 温度の外部境界条件処理の分岐
 @param d_t 温度のデータクラス
 @note 
    - 対流フェイズに関する境界条件はps_Convection_BC()，本メソッドは周期境界と拡散フェイズに関するもの
    - OBC_SYMMETRICは，断熱マスクで処理するため，不要
 */
void SetBC3D::OuterTBC(SklScalar3D<REAL_TYPE>* d_t)
{
  unsigned F=0;
  REAL_TYPE *t=NULL;
  
  if ( !(t = d_t->GetData()) ) Exit(0);
  
  for (int face=0; face<NOFACE; face++) {
    F = obc[face].get_BCtype();

    // 周期境界条件
    if ( F == OBC_PERIODIC ) {
      
      switch ( obc[face].get_PrdcMode() ) {
        case BoundaryOuter::prdc_Simple:
        case BoundaryOuter::prdc_Directional:
          Tobc_Prdc_Simple(d_t, face);
          break;
          
        case BoundaryOuter::prdc_Driver:
          // nothing
          break;
      }
    }
    else { // 周期境界条件以外の処理
      // 対流項寄与分はps_BC_Convection()に実装
    }    
  }
}

/**
 @fn void SetBC3D::OuterTBCface(REAL_TYPE* qbc, unsigned* bh1, REAL_TYPE* t, REAL_TYPE* t0, Control* C, REAL_TYPE& flop)
 @brief 拡散部分に関する温度の外部部境界処理（固体壁の場合）
 @param qbc 境界条件熱流束
 @param bh1 BCindex H1
 @param t n+1時刻の温度場
 @param t0 n時刻の温度場
 @param flop 浮動小数演算数
 @note 断熱BCは断熱マスクで処理
 */
void SetBC3D::OuterTBCface(REAL_TYPE* qbc, unsigned* bh1, REAL_TYPE* t, REAL_TYPE* t0, Control* C, REAL_TYPE& flop)
{
  unsigned H;
  
  for (int face=0; face<NOFACE; face++) {
    
    // 各メソッド内で通信が発生するために，計算領域の最外郭領域でないときに境界処理をスキップする処理はメソッド内で行う
    
    H = obc[face].get_HTmode();
    
    if ( obc[face].get_BCtype() == OBC_WALL ) {
      switch ( obc[face].get_hType() ) { // 熱境界条件の種類
        case HEATFLUX:
          C->H_Dface[face] = ps_OBC_Heatflux(qbc, bh1, face, flop);
          break;
          
        case TRANSFER:
          if      ( H == HT_S)  C->H_Dface[face] = ps_OBC_HeatTransfer_BS(qbc, bh1, face, t, t0, flop);
          else if ( H == HT_SN) C->H_Dface[face] = ps_OBC_HeatTransfer_SN(qbc, bh1, face, t, t0, flop);
          else if ( H == HT_SF) C->H_Dface[face] = ps_OBC_HeatTransfer_SF(qbc, bh1, face, t, t0, flop);
          else if ( H == HT_B)  C->H_Dface[face] = ps_OBC_HeatTransfer_BS(qbc, bh1, face, t, t0, flop);
          break;
          
        case ISOTHERMAL:
          C->H_Dface[face] = ps_OBC_IsoThermal(qbc, bh1, face, t, t0, flop);
          break;
      }      
    }
  }
}

/**
 @fn void SetBC3D::Pobc_Prdc_Simple(SklScalar3D<REAL_TYPE>* d_p, int face)
 @brief 圧力の外部周期境界条件（単純なコピー）
 @param d_p 圧力
 @param face 面番号
 */
void SetBC3D::Pobc_Prdc_Simple(SklScalar3D<REAL_TYPE>* d_p, int face)
{
  SklParaManager* para_mng = ParaCmpo->GetParaManager();
  
  if ( para_mng->IsParallel() ) {
    unsigned no_comm_face = 1; //通信面数
    
    switch (face) {
      case X_MINUS:
        if ( !d_p->CommPeriodicBndCell(PRDC_X_DIR, PRDC_PLUS2MINUS, no_comm_face) ) Exit(0);
        break;
        
      case X_PLUS:
        if ( !d_p->CommPeriodicBndCell(PRDC_X_DIR, PRDC_MINUS2PLUS, no_comm_face) ) Exit(0);
        break;
        
      case Y_MINUS:
        if ( !d_p->CommPeriodicBndCell(PRDC_Y_DIR, PRDC_PLUS2MINUS, no_comm_face) ) Exit(0);
        break;
        
      case Y_PLUS:
        if ( !d_p->CommPeriodicBndCell(PRDC_Y_DIR, PRDC_MINUS2PLUS, no_comm_face) ) Exit(0);
        break;
        
      case Z_MINUS:
        if ( !d_p->CommPeriodicBndCell(PRDC_Z_DIR, PRDC_PLUS2MINUS, no_comm_face) ) Exit(0);
        break;
        
      case Z_PLUS:
        if ( !d_p->CommPeriodicBndCell(PRDC_Z_DIR, PRDC_MINUS2PLUS, no_comm_face) ) Exit(0);
        break;
    }
  }
  else { // Serial
    int i,j,k;
    unsigned m0, m1;
    REAL_TYPE* p=NULL;
    if ( !(p = d_p->GetData()) ) Exit(0);
    
    switch (face) {
      case X_MINUS:
        if( pn.nID[face] < 0 ) {
          for (k=1; k<=kxc; k++) {
            for (j=1; j<=jxc; j++) {
              m0 = FBUtility::getFindexS3D(size, guide, 0  , j  , k  );
              m1 = FBUtility::getFindexS3D(size, guide, ixc, j  , k  );
              p[m0] = p[m1];
            }
          }
        }
        break;
        
      case X_PLUS:
        if( pn.nID[face] < 0 ) {
          for (k=1; k<=kxc; k++) {
            for (j=1; j<=jxc; j++) {
              m0 = FBUtility::getFindexS3D(size, guide, ixc+1, j, k);
              m1 = FBUtility::getFindexS3D(size, guide, 1,     j, k);
              p[m0] = p[m1];
            }
          }
        }
        break;
        
      case Y_MINUS:
        if( pn.nID[face] < 0 ) {
          for (k=1; k<=kxc; k++) {
            for (i=1; i<=ixc; i++) {
              m0 = FBUtility::getFindexS3D(size, guide, i, 0  , k);
              m1 = FBUtility::getFindexS3D(size, guide, i, jxc, k);
              p[m0] = p[m1];
            }
          }
        }
        break;
        
      case Y_PLUS:
        if( pn.nID[face] < 0 ) {
          for (k=1; k<=kxc; k++) {
            for (i=1; i<=ixc; i++) {
              m0 = FBUtility::getFindexS3D(size, guide, i, jxc+1, k);
              m1 = FBUtility::getFindexS3D(size, guide, i, 1    , k);
              p[m0] = p[m1];
            }
          }
        }
        break;
        
      case Z_MINUS:
        if( pn.nID[face] < 0 ) {
          for (j=1; j<=jxc; j++) {
            for (i=1; i<=ixc; i++) {
              m0 = FBUtility::getFindexS3D(size, guide, i, j, 0  );
              m1 = FBUtility::getFindexS3D(size, guide, i, j, kxc);
              p[m0] = p[m1];
            }
          }
        }
        break;
        
      case Z_PLUS:
        if( pn.nID[face] < 0 ) {
          for (j=1; j<=jxc; j++) {
            for (i=1; i<=ixc; i++) {
              m0 = FBUtility::getFindexS3D(size, guide, i, j, kxc+1);
              m1 = FBUtility::getFindexS3D(size, guide, i, j, 1    );
              p[m0] = p[m1];
            }
          }
        }
        break;
    }
  }
}

/**
 @fn void SetBC3D::Pobc_Prdc_Directional(SklScalar3D<REAL_TYPE>* d_p, int face, REAL_TYPE pv, unsigned uod)
 @brief 圧力の外部周期境界条件（双方向に圧力差を設定）
 @param d_p 圧力のデータクラス
 @param face 面番号
 @param pv 圧力差
 @param uod 上流面 or 下流面
 */
void SetBC3D::Pobc_Prdc_Directional(SklScalar3D<REAL_TYPE>* d_p, int face, REAL_TYPE pv, unsigned uod)
{
  SklParaManager* para_mng = ParaCmpo->GetParaManager();
  int i,j,k;
  unsigned m0, m1;
  unsigned no_comm_face = 1; //通信面数
  
  REAL_TYPE pd, *p=NULL;
  if ( !(p = d_p->GetData()) ) Exit(0);
  
  // 上流面か下流面かで，圧力差の方向を逆転する
  pd = ( uod == BoundaryOuter::prdc_upstream ) ? pv : -pv;
  
  if ( para_mng->IsParallel() ) {
    switch (face) {
      case X_MINUS:
        if ( !d_p->CommPeriodicBndCell(PRDC_X_DIR, PRDC_PLUS2MINUS, no_comm_face) ) Exit(0);
        if( pn.nID[face] < 0 ) {
          for (k=1; k<=kxc; k++) {
            for (j=1; j<=jxc; j++) {
              m0 = FBUtility::getFindexS3D(size, guide, 0   , j  , k  );
              p[m0] += pd;
            }
          }
        }
        break;
        
      case X_PLUS:
        if ( !d_p->CommPeriodicBndCell(PRDC_X_DIR, PRDC_MINUS2PLUS, no_comm_face) ) Exit(0);
        if( pn.nID[face] < 0 ) {
          for (k=1; k<=kxc; k++) {
            for (j=1; j<=jxc; j++) {
              m0 = FBUtility::getFindexS3D(size, guide, ixc+1, j, k);
              p[m0] += pd;
            }
          }
        }
        break;
        
      case Y_MINUS:
        if ( !d_p->CommPeriodicBndCell(PRDC_Y_DIR, PRDC_PLUS2MINUS, no_comm_face) ) Exit(0);
        if( pn.nID[face] < 0 ) {
          for (k=1; k<=kxc; k++) {
            for (i=1; i<=ixc; i++) {
              m0 = FBUtility::getFindexS3D(size, guide, i, 0   , k);
              p[m0] += pd;
            }
          }
        }
        break;
        
      case Y_PLUS:
        if ( !d_p->CommPeriodicBndCell(PRDC_Y_DIR, PRDC_MINUS2PLUS, no_comm_face) ) Exit(0);
        if( pn.nID[face] < 0 ) {
          for (k=1; k<=kxc; k++) {
            for (i=1; i<=ixc; i++) {
              m0 = FBUtility::getFindexS3D(size, guide, i, jxc+1, k);
              p[m0] += pd;
            }
          }
        }
        break;
        
      case Z_MINUS:
        if ( !d_p->CommPeriodicBndCell(PRDC_Z_DIR, PRDC_PLUS2MINUS, no_comm_face) ) Exit(0);
        if( pn.nID[face] < 0 ) {
          for (j=1; j<=jxc; j++) {
            for (i=1; i<=ixc; i++) {
              m0 = FBUtility::getFindexS3D(size, guide, i, j, 0   );
              p[m0] += pd;
            }
          }
        }
        break;
        
      case Z_PLUS:
        if ( !d_p->CommPeriodicBndCell(PRDC_Z_DIR, PRDC_MINUS2PLUS, no_comm_face) ) Exit(0);
        if( pn.nID[face] < 0 ) {
          for (j=1; j<=jxc; j++) {
            for (i=1; i<=ixc; i++) {
              m0 = FBUtility::getFindexS3D(size, guide, i, j, kxc+1);
              p[m0] += pd;
            }
          }
        }
        break;
    }
  }
  else { // Serial
    switch (face) {
      case X_MINUS:
        if( pn.nID[face] < 0 ) {
          for (k=1; k<=kxc; k++) {
            for (j=1; j<=jxc; j++) {
              m0 = FBUtility::getFindexS3D(size, guide, 0   , j  , k  );
              m1 = FBUtility::getFindexS3D(size, guide, ixc, j  , k  );
              p[m0] = p[m1] + pd;
            }
          }
        }
        break;
        
      case X_PLUS:
        if( pn.nID[face] < 0 ) {
          for (k=1; k<=kxc; k++) {
            for (j=1; j<=jxc; j++) {
              m0 = FBUtility::getFindexS3D(size, guide, ixc+1, j, k);
              m1 = FBUtility::getFindexS3D(size, guide, 1,      j, k);
              p[m0] = p[m1] + pd;
            }
          }
        }
        break;
        
      case Y_MINUS:
        if( pn.nID[face] < 0 ) {
          for (k=1; k<=kxc; k++) {
            for (i=1; i<=ixc; i++) {
              m0 = FBUtility::getFindexS3D(size, guide, i, 0   , k);
              m1 = FBUtility::getFindexS3D(size, guide, i, jxc, k);
              p[m0] = p[m1] + pd;
            }
          }
        }
        break;
        
      case Y_PLUS:
        if( pn.nID[face] < 0 ) {
          for (k=1; k<=kxc; k++) {
            for (i=1; i<=ixc; i++) {
              m0 = FBUtility::getFindexS3D(size, guide, i, jxc+1, k);
              m1 = FBUtility::getFindexS3D(size, guide, i, 1     , k);
              p[m0] = p[m1] + pd;
            }
          }
        }
        break;
        
      case Z_MINUS:
        if( pn.nID[face] < 0 ) {
          for (j=1; j<=jxc; j++) {
            for (i=1; i<=ixc; i++) {
              m0 = FBUtility::getFindexS3D(size, guide, i, j, 0   );
              m1 = FBUtility::getFindexS3D(size, guide, i, j, kxc);
              p[m0] = p[m1] + pd;
            }
          }
        }
        break;
        
      case Z_PLUS:
        if( pn.nID[face] < 0 ) {
          for (j=1; j<=jxc; j++) {
            for (i=1; i<=ixc; i++) {
              m0 = FBUtility::getFindexS3D(size, guide, i, j, kxc+1);
              m1 = FBUtility::getFindexS3D(size, guide, i, j, 1     );
              p[m0] = p[m1] + pd;
            }
          }
        }
        break;
    }
  }
}

/**
 @fn void SetBC3D::Pibc_Prdc(SklScalar3D<REAL_TYPE>* d_p, int* st, int* ed, SklScalar3D<unsigned>* d_bcd, int odr, int dir, REAL_TYPE pd)
 @brief 圧力の内部周期境界条件（一方向の圧力差）
 @param d_p 圧力のデータクラス
 @param st コンポーネント範囲の開始インデクス
 @param ed コンポーネント範囲の終了インデクス
 @param d_bcd BCindex ID
 @param odr 下流面のidが格納されているコンポーネントエントリ
 @param dir ドライバの設置方向
 @param pd 圧力差
 */
void SetBC3D::Pibc_Prdc(SklScalar3D<REAL_TYPE>* d_p, int* st, int* ed, SklScalar3D<unsigned>* d_bcd, int odr, int dir, REAL_TYPE pd)
{
  SklParaManager* para_mng = ParaCmpo->GetParaManager();
  if( para_mng->IsParallel() ){
    Hostonly_ printf("Error : 'Pibc_Prdc' method is limited to use for serial execution\n.");
    Exit(0);
  }
  
  int i,j,k;
  unsigned m0, m1, *bx=NULL;
  REAL_TYPE* p=NULL;
  
  if ( !(p = d_p->GetData()) ) Exit(0);
  if ( !(bx= d_bcd->GetData()) ) Exit(0);
  
  switch (dir) {
    case X_MINUS:
      if ( pn.nID[dir] < 0 ) {
        i = st[0];
        for (k=st[2]; k<=ed[2]; k++) {
          for (j=st[1]; j<=ed[1]; j++) {
            m1 = FBUtility::getFindexS3D(size, guide, i, j, k);
            if ( DECODE_CMP(bx[m1]) == odr ) {
              m0 = FBUtility::getFindexS3D(size, guide, 0, j, k);
              p[m0] = p[m1] + pd;
            }
          }
        }
      }
      break;
      
    case X_PLUS:
      if ( pn.nID[dir] < 0 ) {
        i = st[0];
        for (k=st[2]; k<=ed[2]; k++) {
          for (j=st[1]; j<=ed[1]; j++) {
            m1 = FBUtility::getFindexS3D(size, guide, i, j, k);
            if ( DECODE_CMP(bx[m1]) == odr ) {
              m0 = FBUtility::getFindexS3D(size, guide, ixc+1, j, k);
              p[m0] = p[m1] + pd;              
            }
          }
        }
      }
      break;
      
    case Y_MINUS:
      if ( pn.nID[dir] < 0 ) {
        j = st[1];
        for (k=st[2]; k<=ed[2]; k++) {
          for (i=st[0]; i<=ed[0]; i++) {
            m1 = FBUtility::getFindexS3D(size, guide, i, j, k);
            if ( DECODE_CMP(bx[m1]) == odr ) {
              m0 = FBUtility::getFindexS3D(size, guide, i, 0, k);
              p[m0] = p[m1] + pd;
            }
          }
        }
      }
      break;
      
    case Y_PLUS:
      if ( pn.nID[dir] < 0 ) {
        j = st[1];
        for (k=st[2]; k<=ed[2]; k++) {
          for (i=st[0]; i<=ed[0]; i++) {
            m1 = FBUtility::getFindexS3D(size, guide, i, j, k);
            if ( DECODE_CMP(bx[m1]) == odr ) {
              m0 = FBUtility::getFindexS3D(size, guide, i, jxc+1, k);
              p[m0] = p[m1] + pd;
            }
          }
        }
      }
      break;
      
    case Z_MINUS:
      if ( pn.nID[dir] < 0 ) {
        k = st[2];
        for (j=st[1]; j<=ed[1]; j++) {
          for (i=st[0]; i<=ed[0]; i++) {
            m1 = FBUtility::getFindexS3D(size, guide, i, j, k);
            if ( DECODE_CMP(bx[m1]) == odr ) {
              m0 = FBUtility::getFindexS3D(size, guide, i, j, 0);
              p[m0] = p[m1] + pd;
            }
          }
        }
      }
      break;
      
    case Z_PLUS:
      if ( pn.nID[dir] < 0 ) {
        k = st[2];
        for (j=st[1]; j<=ed[1]; j++) {
          for (i=st[0]; i<=ed[0]; i++) {
            m1 = FBUtility::getFindexS3D(size, guide, i, j, k);
            if ( DECODE_CMP(bx[m1]) == odr ) {
              m0 = FBUtility::getFindexS3D(size, guide, i, j, kxc+1);
              p[m0] = p[m1] + pd;
            }
          }
        }
      }
      break;
  }
}

/**
 @fn REAL_TYPE SetBC3D::setDirectForcing(REAL_TYPE* v, unsigned* bx, unsigned n, REAL_TYPE v00)
 @brief 内部領域の速度のディリクレ境界をセットする
 @param v 速度ベクトル
 @param bx BCindex
 @param n 境界条件コンポーネントのエントリ番号
 @param v00 参照速度
 @retval 無次元平均流速
 @note
    - 内部境界はセル面直を仮定
 
REAL_TYPE SetBC3D::setDirectForcing(REAL_TYPE* v, unsigned* bx, unsigned n, REAL_TYPE v00)
{
  int i,j,k;
  unsigned s;
  REAL_TYPE vel, vec[3];
  
  vel   = FBUtility::convD2ND_V(cmp[n].D1.Velocity, RefV)  * v00;
  vec[0]= cmp[n].nv[0]*vel;
  vec[1]= cmp[n].nv[1]*vel;
  vec[2]= cmp[n].nv[2]*vel;
  
  for (k=cmp[n].ci.st[2]-1; k<=cmp[n].ci.ed[2]; k++) {
    for (j=cmp[n].ci.st[1]-1; j<=cmp[n].ci.ed[1]; j++) {
      for (i=cmp[n].ci.st[0]-1; i<=cmp[n].ci.ed[0]; i++) {
        s = bx[ FBUtility::getFindexS3D(size, guide, i, j, k) ];
        if ( (s & BIT_MASK_10) == n ) {
					if ( BIT_IS_SHIFT(s, VFACE_I) ) {
						v[FBUtility::getFindexV3D(size, guide, i, j, k, 0)] = vec[0];
					}
					if ( BIT_IS_SHIFT(s, VFACE_J) ) {
						v[FBUtility::getFindexV3D(size, guide, i, j, k, 1)] = vec[1];
					}
					if ( BIT_IS_SHIFT(s, VFACE_K) ) {
						v[FBUtility::getFindexV3D(size, guide, i, j, k, 2)] = vec[2];
					}
        }
      }
    }
  }

  return (vel);
}*/

/**
 @fn REAL_TYPE SetBC3D::ps_IBC_HeatGen_SM(REAL_TYPE* t, unsigned* bh2, int n, REAL_TYPE dt, REAL_TYPE& flop)
 @brief 発熱境界条件
 @param t 温度場
 @param bh2 BCindex H2
 @param n 境界条件コンポーネントのエントリ番号
 @param dt 時間積分幅
 @param flop 浮動小数演算数
 @note
    - 発熱密度はControl::setParameters()で計算 D2に発熱密度が保存されている
 */
REAL_TYPE SetBC3D::ps_IBC_HeatGen_SM(REAL_TYPE* t, unsigned* bh2, int n, REAL_TYPE dt, REAL_TYPE& flop)
{
  int i,j,k;
  unsigned register s, m, c=0;
  REAL_TYPE hs;
  int st[3], ed[3];

  //odr= cmp[n].getMatOdr();
  //hs = dt * FBUtility::convD2ND_Hsrc(cmp[n].get_HeatDensity(), RefV, RefL, DiffTemp, mat[odr].P[p_density], mat[odr].P[p_specific_heat]);
  hs = FBUtility::convD2ND_Hsrc(cmp[n].get_HeatDensity(), RefV, RefL, DiffTemp, rho, cp);

  cmp[n].getBbox(st, ed);
	
  for (k=st[2]; k<=ed[2]; k++) {
    for (j=st[1]; j<=ed[1]; j++) {
      for (i=st[0]; i<=ed[0]; i++) {
        m = FBUtility::getFindexS3D(size, guide, i  , j  , k  );
        if ( (bh2[m] & MASK_6) == n ) {
          t[m] += dt*hs;
          c++; 
        }
      }
    }
  }
  flop += (REAL_TYPE)c*2.0;
  
  return (c*hs); // 無次元の総発熱量(単位体積あたり)
}

/**
 @fn void SetBC3D::ps_IBC_ConstTemp(REAL_TYPE* t, unsigned* bh2, int n)
 @brief 温度一定の境界条件
 @param t 温度場
 @param bh2 BCindex H2
 @param n 境界条件コンポーネントのエントリ番号
 */
void SetBC3D::ps_IBC_ConstTemp(REAL_TYPE* t, unsigned* bh2, int n)
{
  int i,j,k;
  unsigned m;
  REAL_TYPE tmp;
  int st[3], ed[3];

  tmp = FBUtility::convK2ND(cmp[n].get_Temp(), BaseTemp, DiffTemp);

  cmp[n].getBbox(st, ed);
	
  for (k=st[2]; k<=ed[2]; k++) {
    for (j=st[1]; j<=ed[1]; j++) {
      for (i=st[0]; i<=ed[0]; i++) {
        m = FBUtility::getFindexS3D(size, guide, i  , j  , k  );
        if ( (bh2[m] & MASK_6) == n ) t[m] = tmp;
      }
    }
  }
}

/**
 @fn void SetBC3D::ps_BC_Convection(REAL_TYPE* ws, unsigned* bh1, REAL_TYPE* v, REAL_TYPE* t, REAL_TYPE tm, Control* C, REAL_TYPE* v00, REAL_TYPE& flop)
 @brief 対流項計算時の流束型の境界条件処理
 @param ws 温度の増分
 @param bh1 BCindex H1
 @param v 速度
 @param t 温度
 @param tm 時刻
 @param C コントロールクラス
 @param v00
 @pram flop 浮動小数演算数
 @note
    - 熱量(-)はvalに保存
 */
void SetBC3D::ps_BC_Convection(REAL_TYPE* ws, unsigned* bh1, REAL_TYPE* v, REAL_TYPE* t, REAL_TYPE tm, Control* C, REAL_TYPE* v00, REAL_TYPE& flop)
{
  REAL_TYPE vec[3];
  
  // 外部
  for (int face=0; face<NOFACE; face++) {
    
    // 各メソッド内で通信が発生するために，計算領域の最外郭領域でないときに境界処理をスキップする処理はメソッド内で行う
    
    switch ( obc[face].get_BCtype() ) {
      case OBC_OUTFLOW:
        C->H_Dface[face] = ps_OBC_Free(ws, bh1, face, v, t, v00, flop);
        break;
        
      case OBC_IN_OUT:
      case OBC_TRC_FREE:
        C->H_Dface[face] = ps_OBC_Free(ws, bh1, face, v, t, v00, flop);
        break;

      case OBC_SPEC_VEL:
        C->H_Dface[face] = ps_OBC_SpecVH(ws, bh1, face, t, tm, v00, flop);
        break;
    }
  }

  // 内部
  for (int n=1; n<=NoBC; n++) {
    switch ( cmp[n].getType() ) {
      case SPEC_VEL_WH:
        extractVel_IBC(n, vec, tm, v00, flop);
        cmp[n].set_Mon_Calorie( ps_IBC_SpecVH(ws, bh1, n, v00[0], vec, flop) );
        break;
          
      case OUTFLOW:
        cmp[n].set_Mon_Calorie( ps_IBC_Outflow(ws, bh1, n, v, t, v00, flop) );
        break;
    }
  }
}

/**
 @fn REAL_TYPE SetBC3D::ps_SpecVH(REAL_TYPE* ws, unsigned* bh1, int n, REAL_TYPE v00, REAL_TYPE* vec, REAL_TYPE& flop)
 @brief 内部領域の速度と温度の指定境界条件
 @retval 熱量(-)
 @param ws 温度増分
 @param bh1 BCindex H1
 @param n コンポーネントリストのインデクス
 @param v00 参照速度
 @param vec[3] 指定ベクトル
 @param flop 浮動小数演算数
 @note
    - モニタ量の熱量va(無次元)は系に対する流入量なので，基準温度に対する熱量
 */
REAL_TYPE SetBC3D::ps_IBC_SpecVH(REAL_TYPE* ws, unsigned* bh1, int n, REAL_TYPE v00, REAL_TYPE* vec, REAL_TYPE& flop)
{
  SklParaManager* para_mng = ParaCmpo->GetParaManager();
  int i,j,k;
  unsigned register s, m;
  REAL_TYPE va=0.0;
  REAL_TYPE hu, hv, hw, tc, dhd=dh*RefL;
  int st[3], ed[3];
  REAL_TYPE f_e, f_w, f_n, f_s, f_t, f_b, ff;
  REAL_TYPE dh1 = 1.0/dh;
	
  tc  = FBUtility::convK2ND(cmp[n].get_Temp(), BaseTemp, DiffTemp); // difference form BaseTemp 
  hu  = vec[0]*tc;
  hv  = vec[1]*tc;
  hw  = vec[2]*tc;
  
  cmp[n].getBbox(st, ed);
	
  for (k=st[2]; k<=ed[2]; k++) {
    for (j=st[1]; j<=ed[1]; j++) {
      for (i=st[0]; i<=ed[0]; i++) {
        
        f_e = f_w = f_n = f_s = f_t = f_b = 0.0; // 境界条件以外の寄与をゼロにしておく
        m = FBUtility::getFindexS3D(size, guide, i, j, k);
        s = bh1[m];
        
        if ( GET_FACE_BC(s, BC_FACE_W) == n ) {
          f_w = hu;
          va += hu;
        }
        
        if ( GET_FACE_BC(s, BC_FACE_E) == n ) {
          f_e = hu;
          va += hu;
        }
        
        if ( GET_FACE_BC(s, BC_FACE_S) == n ) {
          f_s = hv;
          va += hv;
        }
        
        if ( GET_FACE_BC(s, BC_FACE_N) == n ) {
          f_n = hv;
          va += hv;
        }
        
        if ( GET_FACE_BC(s, BC_FACE_B) == n ) {
          f_b = hw;
          va += hw;
        }
        
        if ( GET_FACE_BC(s, BC_FACE_T) == n ) {
          f_t = hw;
          va += hw;
        }
        
        ff = f_e - f_w + f_n - f_s + f_t - f_b;
        ws[m] -= ff*dh1; // 流入境界は必ずFセルなのでマスク不要
      }
    }
  }
  
  flop += (REAL_TYPE)( (ed[0]-st[0]+1)*(ed[1]-st[1]+1)*(ed[2]-st[2]+1)*7 ); // 近似
  
  if( para_mng->IsParallel() ){
		REAL_TYPE tmp = va;
		para_mng->Allreduce(&tmp, &va, 1, SKL_ARRAY_DTYPE_REAL, SKL_SUM, pn.procGrp);
	}
  
  return (va); // 単位面積あたりの無次元熱流束の和
}

/**
 @fn REAL_TYPE SetBC3D::ps_IBC_Outflow(REAL_TYPE* ws, unsigned* bh1, int n, REAL_TYPE* v, REAL_TYPE* t, REAL_TYPE* v00, REAL_TYPE& flop)
 @brief 内部領域のOutflowの境界条件処理
 @retval 熱量(W)
 @param ws 温度増分
 @param bh1 BCindex H1
 @param n コンポーネントリストのインデクス
 @param v 速度
 @param t 温度
 @param v00 参照速度
 @param flop 浮動小数演算数
 @note
    - モニタ量の熱量va(W)は系に対する流入量なので，基準温度に対する熱量
 */
REAL_TYPE SetBC3D::ps_IBC_Outflow(REAL_TYPE* ws, unsigned* bh1, int n, REAL_TYPE* v, REAL_TYPE* t, REAL_TYPE* v00, REAL_TYPE& flop)
{
  SklParaManager* para_mng = ParaCmpo->GetParaManager();
  int i,j,k;
  unsigned register s, m;
  REAL_TYPE va=0.0;
  REAL_TYPE dhd=dh*RefL;
  int st[3], ed[3];
  REAL_TYPE f_e, f_w, f_n, f_s, f_t, f_b, ff, c;
  REAL_TYPE dh1 = 1.0/dh;
  REAL_TYPE u_ref, v_ref, w_ref, t_p;
  unsigned m_e, m_w, m_n, m_s, m_t, m_b, m_0;
  
  u_ref = v00[1];
  v_ref = v00[2];
  w_ref = v00[3];
  
  cmp[n].getBbox(st, ed);
	
  for (k=st[2]; k<=ed[2]; k++) {
    for (j=st[1]; j<=ed[1]; j++) {
      for (i=st[0]; i<=ed[0]; i++) {
        
        f_e = f_w = f_n = f_s = f_t = f_b = 0.0; // 境界条件以外の寄与をゼロにしておく
        m = FBUtility::getFindexS3D(size, guide, i, j, k);
        s = bh1[m];
        t_p = t[m];
        
        if ( GET_FACE_BC(s, BC_FACE_W) == n ) {
          m_0 = FBUtility::getFindexV3DEx(size, guide, 0, i  , j, k);
          m_w = FBUtility::getFindexV3DEx(size, guide, 0, i-1, j, k);
          c = 0.5*(v[m_w]+v[m_0]) - u_ref;
          if ( c>0.0 ) c=0.0;
          f_w = c*t_p;
          va += f_w;
          flop += 5.0;
        }
        
        if ( GET_FACE_BC(s, BC_FACE_E) == n ) {
          m_0 = FBUtility::getFindexV3DEx(size, guide, 0, i  , j, k);
          m_e = FBUtility::getFindexV3DEx(size, guide, 0, i+1, j, k);
          c = 0.5*(v[m_e]+v[m_0]) - u_ref;
          if ( c<0.0 ) c=0.0;
          f_e = c*t_p;
          va += f_e;
          flop += 5.0;
        }
        
        if ( GET_FACE_BC(s, BC_FACE_S) == n ) {
          m_0 = FBUtility::getFindexV3DEx(size, guide, 1, i, j  , k);
          m_s = FBUtility::getFindexV3DEx(size, guide, 1, i, j-1, k);
          c = 0.5*(v[m_s]+v[m_0]) - v_ref;
          if ( c>0.0 ) c=0.0;
          f_s = c*t_p;
          va += f_s;
          flop += 5.0;
        }
        
        if ( GET_FACE_BC(s, BC_FACE_N) == n ) {
          m_0 = FBUtility::getFindexV3DEx(size, guide, 1, i, j  , k);
          m_n = FBUtility::getFindexV3DEx(size, guide, 1, i, j+1, k);
          c = 0.5*(v[m_n]+v[m_0]) - v_ref;
          if ( c<0.0 ) c=0.0;
          f_n = c*t_p;
          va += f_n;
          flop += 5.0;
        }
        
        if ( GET_FACE_BC(s, BC_FACE_B) == n ) {
          m_0 = FBUtility::getFindexV3DEx(size, guide, 2, i, j, k  );
          m_b = FBUtility::getFindexV3DEx(size, guide, 2, i, j, k-1);
          c = 0.5*(v[m_b]+v[m_0]) - w_ref;
          if ( c>0.0 ) c=0.0;
          f_b = c*t_p;
          va += f_b;
          flop += 5.0;
        }
        
        if ( GET_FACE_BC(s, BC_FACE_T) == n ) {
          m_0 = FBUtility::getFindexV3DEx(size, guide, 2, i, j, k  );
          m_t = FBUtility::getFindexV3DEx(size, guide, 2, i, j, k+1);
          c = 0.5*(v[m_t]+v[m_0]) - w_ref;
          if ( c<0.0 ) c=0.0;
          f_t = c*t_p;
          va += f_t;
          flop += 5.0;
        }
        
        ff = f_e - f_w + f_n - f_s + f_t - f_b;
        ws[m] -= ff*dh1; // Outflow指定セルは，必ずFセルなのでマスク不要
        flop += 8.0;
      }
    }
  }
  
  if( para_mng->IsParallel() ){
		REAL_TYPE tmp = va;
		para_mng->Allreduce(&tmp, &va, 1, SKL_ARRAY_DTYPE_REAL, SKL_SUM, pn.procGrp);
	}
  
  return (va); // 単位面積あたりの無次元熱流束の和
}

/**
 @fn REAL_TYPE SetBC3D::ps_OBC_Free(REAL_TYPE* ws, unsigned* bh1, int face, REAL_TYPE* v, REAL_TYPE* t, REAL_TYPE* v00, REAL_TYPE& flop)
 @brief 外部領域のOutflow, In_out, T_Freeの境界条件処理
 @retval 熱量(-)
 @param ws 温度増分
 @param bh1 BCindex H1
 @param face 外部境界面番号
 @param v 速度
 @param t 温度
 @param v00 参照速度
 @param flop 浮動小数演算数
 @note
    - モニタ量の熱量va(-)は系に対する流入量なので，基準温度に対する熱量
 */
REAL_TYPE SetBC3D::ps_OBC_Free(REAL_TYPE* ws, unsigned* bh1, int face, REAL_TYPE* v, REAL_TYPE* t, REAL_TYPE* v00, REAL_TYPE& flop)
{
  SklParaManager* para_mng = ParaCmpo->GetParaManager();
  int i,j,k;
  unsigned register s, m;
  REAL_TYPE va=0.0, t_nd;
  REAL_TYPE dhd=dh*RefL;
  REAL_TYPE f_e, f_w, f_n, f_s, f_t, f_b, ff, c;
  REAL_TYPE dh1 = 1.0/dh;
  REAL_TYPE u_ref, v_ref, w_ref, t_p;
  unsigned m_e, m_w, m_n, m_s, m_t, m_b, m_0;
  
  u_ref = v00[1];
  v_ref = v00[2];
  w_ref = v00[3];
  
  f_e = f_w = f_n = f_s = f_t = f_b = 0.0; // 境界条件以外の寄与をゼロにしておく
  t_nd = FBUtility::convK2ND(obc[face].get_Temp(), BaseTemp, DiffTemp);
  
  switch (face) {
    case X_MINUS:
      if( pn.nID[face] < 0 ){
        i=1; // 最外層のID
        for (k=1; k<=kxc; k++) {
          for (j=1; j<=jxc; j++) {
            m   = FBUtility::getFindexS3D(size, guide, i, j, k);
            s = bh1[m];
            m_0 = FBUtility::getFindexV3DEx(size, guide, 0, i  , j, k);
            m_w = FBUtility::getFindexV3DEx(size, guide, 0, i-1, j, k);
            c = 0.5*(v[m_w]+v[m_0]) - u_ref;
            t_p = (c < 0.0) ? t[m] : t_nd; // 流出の場合には内部の値，流入の場合には参照値
            f_w = c*t_p;
            va += f_w;
            ff = f_e - f_w + f_n - f_s + f_t - f_b;
            ws[m] -= ff*dh1*GET_SHIFT_F(s, STATE_BIT);
          }
        }
      }
      flop += (REAL_TYPE)(jxc*kxc*13);
      break;
      
    case X_PLUS:
      if( pn.nID[face] < 0 ){
        i=ixc;
        for (k=1; k<=kxc; k++) {
          for (j=1; j<=jxc; j++) {
            m = FBUtility::getFindexS3D(size, guide, i, j, k);
            s = bh1[m];
            m_0 = FBUtility::getFindexV3DEx(size, guide, 0, i  , j, k);
            m_e = FBUtility::getFindexV3DEx(size, guide, 0, i+1, j, k);
            c = 0.5*(v[m_e]+v[m_0]) - u_ref;
            t_p = (c > 0.0) ? t[m] : t_nd;
            f_e = c*t_p;
            va += f_e;
            ff = f_e - f_w + f_n - f_s + f_t - f_b;
            ws[m] -= ff*dh1*GET_SHIFT_F(s, STATE_BIT);
          }
        }
      }
      flop += (REAL_TYPE)(jxc*kxc*13);
      break;
      
    case Y_MINUS:
      if( pn.nID[face] < 0 ){
        j=1;
        for (k=1; k<=kxc; k++) {
          for (i=1; i<=ixc; i++) {
            m = FBUtility::getFindexS3D(size, guide, i, j, k);
            s = bh1[m];
            m_0 = FBUtility::getFindexV3DEx(size, guide, 1, i, j  , k);
            m_s = FBUtility::getFindexV3DEx(size, guide, 1, i, j-1, k);
            c = 0.5*(v[m_s]+v[m_0]) - v_ref;
            t_p = (c < 0.0) ? t[m] : t_nd;
            f_s = c*t_p;
            va += f_s;
            ff = f_e - f_w + f_n - f_s + f_t - f_b;
            ws[m] -= ff*dh1*GET_SHIFT_F(s, STATE_BIT);
          }
        }
      }
      flop += (REAL_TYPE)(ixc*kxc*13);
      break;
      
    case Y_PLUS:
      if( pn.nID[face] < 0 ){
        j=jxc;
        for (k=1; k<=kxc; k++) {
          for (i=1; i<=ixc; i++) {
            m = FBUtility::getFindexS3D(size, guide, i, j, k);
            s = bh1[m];
            m_0 = FBUtility::getFindexV3DEx(size, guide, 1, i, j  , k);
            m_n = FBUtility::getFindexV3DEx(size, guide, 1, i, j+1, k);
            c = 0.5*(v[m_n]+v[m_0]) - v_ref;
            t_p = (c > 0.0) ? t[m] : t_nd;
            f_n = c*t_p;
            va += f_n;
            ff = f_e - f_w + f_n - f_s + f_t - f_b;
            ws[m] -= ff*dh1*GET_SHIFT_F(s, STATE_BIT);
          }
        }
      }
      flop += (REAL_TYPE)(ixc*kxc*13);
      break;
      
    case Z_MINUS:
      if( pn.nID[face] < 0 ){
        k=1;
        for (j=1; j<=jxc; j++) {
          for (i=1; i<=ixc; i++) {
            m = FBUtility::getFindexS3D(size, guide, i, j, k);
            s = bh1[m];
            m_0 = FBUtility::getFindexV3DEx(size, guide, 2, i, j, k  );
            m_b = FBUtility::getFindexV3DEx(size, guide, 2, i, j, k-1);
            c = 0.5*(v[m_b]+v[m_0]) - w_ref;
            t_p = (c < 0.0) ? t[m] : t_nd;
            f_b = c*t_p;
            va += f_b;
            ff = f_e - f_w + f_n - f_s + f_t - f_b;
            ws[m] -= ff*dh1*GET_SHIFT_F(s, STATE_BIT);
          }
        }
      }
      flop += (REAL_TYPE)(ixc*jxc*13);
      break;
      
    case Z_PLUS:
      if( pn.nID[face] < 0 ){
        k=kxc;
        for (j=1; j<=jxc; j++) {
          for (i=1; i<=ixc; i++) {
            m = FBUtility::getFindexS3D(size, guide, i, j, k);
            s = bh1[m];
            m_0 = FBUtility::getFindexV3DEx(size, guide, 2, i, j, k  );
            m_t = FBUtility::getFindexV3DEx(size, guide, 2, i, j, k+1);
            c = 0.5*(v[m_t]+v[m_0]) - w_ref;
            t_p = (c > 0.0) ? t[m] : t_nd;
            f_t = c*t_p;
            va += f_t;
            ff = f_e - f_w + f_n - f_s + f_t - f_b;
            ws[m] -= ff*dh1*GET_SHIFT_F(s, STATE_BIT);
          }
        }
      }
      flop += (REAL_TYPE)(ixc*jxc*13);
      break;
  }
  
  if( para_mng->IsParallel() ){
		REAL_TYPE tmp = va;
		para_mng->Allreduce(&tmp, &va, 1, SKL_ARRAY_DTYPE_REAL, SKL_SUM, pn.procGrp);
	}
  
  return (va);
}

/**
 @fn REAL_TYPE SetBC3D::ps_OBC_SpecVH(REAL_TYPE* ws, unsigned* bh1, int face, REAL_TYPE* t, REAL_TYPE tm, REAL_TYPE* v00, REAL_TYPE& flop)
 @brief 外部領域の速度指定の境界条件処理
 @retval 熱量(-)
 @param ws 温度増分
 @param bh1 BCindex H1
 @param face 外部境界面番号
 @param t 温度
 @param tm 時刻
 @param v00
 @param flop 浮動小数演算数
 @note
    - モニタ量の熱量va(-)は系に対する流入量なので，基準温度に対する熱量
 */
REAL_TYPE SetBC3D::ps_OBC_SpecVH(REAL_TYPE* ws, unsigned* bh1, int face, REAL_TYPE* t, REAL_TYPE tm, REAL_TYPE* v00, REAL_TYPE& flop)
{
  SklParaManager* para_mng = ParaCmpo->GetParaManager();
  int i,j,k;
  unsigned register s, m;
  REAL_TYPE va=0.0, t_nd;
  REAL_TYPE f_e, f_w, f_n, f_s, f_t, f_b, ff, c;
  REAL_TYPE dh1 = 1.0/dh;
  REAL_TYPE u_ref, v_ref, w_ref;
  REAL_TYPE vec[3];
  
  u_ref = v00[1];
  v_ref = v00[2];
  w_ref = v00[3];
  
  f_e = f_w = f_n = f_s = f_t = f_b = 0.0; // 境界条件以外の寄与をゼロにしておく
  t_nd = FBUtility::convK2ND(obc[face].get_Temp(), BaseTemp, DiffTemp);

  switch (face) {
    case X_MINUS:
      if( pn.nID[face] < 0 ){
        extractVel_OBC(face, vec, tm, v00, flop);
        c = vec[0] - u_ref;
        f_w = c * t_nd;
        
        i=1; // 最外層のID
        for (k=1; k<=kxc; k++) {
          for (j=1; j<=jxc; j++) {
            m   = FBUtility::getFindexS3D(size, guide, i, j, k);
            s = bh1[m];
            va += f_w;
            ff = f_e - f_w + f_n - f_s + f_t - f_b;
            ws[m] -= ff*dh1*GET_SHIFT_F(s, STATE_BIT);
          }
        }
      }
      flop += (REAL_TYPE)(jxc*kxc*9 + 2);
      break;
      
    case X_PLUS:
      if( pn.nID[face] < 0 ){
        extractVel_OBC(face, vec, tm, v00, flop);
        c = vec[0] - u_ref;
        f_e = -c * t_nd; // ref_tが正のとき，X+方向のセルでは流出なので符号反転
        
        i=ixc;
        for (k=1; k<=kxc; k++) {
          for (j=1; j<=jxc; j++) {
            m = FBUtility::getFindexS3D(size, guide, i, j, k);
            s = bh1[m];
            va += f_e;
            ff = f_e - f_w + f_n - f_s + f_t - f_b;
            ws[m] -= ff*dh1*GET_SHIFT_F(s, STATE_BIT);
          }
        }
      }
      flop += (REAL_TYPE)(jxc*kxc*9 + 2);
      break;
      
    case Y_MINUS:
      if( pn.nID[face] < 0 ){
        extractVel_OBC(face, vec, tm, v00, flop);
        c = vec[1] - v_ref;
        f_s = c * t_nd;
        
        j=1;
        for (k=1; k<=kxc; k++) {
          for (i=1; i<=ixc; i++) {
            m = FBUtility::getFindexS3D(size, guide, i, j, k);
            s = bh1[m];
            va += f_s;
            ff = f_e - f_w + f_n - f_s + f_t - f_b;
            ws[m] -= ff*dh1*GET_SHIFT_F(s, STATE_BIT);
          }
        }
      }
      flop += (REAL_TYPE)(ixc*kxc*9 + 2);
      break;
      
    case Y_PLUS:
      if( pn.nID[face] < 0 ){
        extractVel_OBC(face, vec, tm, v00, flop);
        c = vec[1] - v_ref;
        f_n = -c * t_nd; // 符号反転
        
        j=jxc;
        for (k=1; k<=kxc; k++) {
          for (i=1; i<=ixc; i++) {
            m = FBUtility::getFindexS3D(size, guide, i, j, k);
            s = bh1[m];
            va += f_n;
            ff = f_e - f_w + f_n - f_s + f_t - f_b;
            ws[m] -= ff*dh1*GET_SHIFT_F(s, STATE_BIT);
          }
        }
      }
      flop += (REAL_TYPE)(ixc*kxc*9 + 2);
      break;
      
    case Z_MINUS:
      if( pn.nID[face] < 0 ){
        extractVel_OBC(face, vec, tm, v00, flop);
        c = vec[2] - w_ref;
        f_b = c * t_nd;
        
        k=1;
        for (j=1; j<=jxc; j++) {
          for (i=1; i<=ixc; i++) {
            m = FBUtility::getFindexS3D(size, guide, i, j, k);
            s = bh1[m];
            va += f_b;
            ff = f_e - f_w + f_n - f_s + f_t - f_b;
            ws[m] -= ff*dh1*GET_SHIFT_F(s, STATE_BIT);
          }
        }
      }
      flop += (REAL_TYPE)(ixc*jxc*9 + 2);
      break;
      
    case Z_PLUS:
      if( pn.nID[face] < 0 ){
        extractVel_OBC(face, vec, tm, v00, flop);
        c = vec[2] - w_ref;
        f_t = -c * t_nd; // 符号反転
        
        k=kxc;
        for (j=1; j<=jxc; j++) {
          for (i=1; i<=ixc; i++) {
            m = FBUtility::getFindexS3D(size, guide, i, j, k);
            s = bh1[m];
            va += f_t;
            ff = f_e - f_w + f_n - f_s + f_t - f_b;
            ws[m] -= ff*dh1*GET_SHIFT_F(s, STATE_BIT);
          }
        }
      }
      flop += (REAL_TYPE)(ixc*jxc*9 + 2);
      break;
  }
  
  if( para_mng->IsParallel() ){
		REAL_TYPE tmp = va;
		para_mng->Allreduce(&tmp, &va, 1, SKL_ARRAY_DTYPE_REAL, SKL_SUM, pn.procGrp);
	}
  
  return (va);
}

/**
 @fn REAL_TYPE SetBC3D::ps_OBC_IsoThermal(REAL_TYPE* qbc, unsigned* bh1, int face, REAL_TYPE* t, REAL_TYPE* t0, REAL_TYPE& flop)
 @brief 外部領域の等温熱流束の境界条件処理
 @retval 熱量(-)
 @param qbc 境界条件熱流束
 @param bh1 BCindex H1
 @param face 外部境界面番号
 @param t n+1時刻の温度場
 @param t0 n時刻の温度場
 @param flop 浮動小数演算数
 @note
 - モニタ量の熱量va(-)は系に対する流入出量
 - 流入を正にとる
 */
REAL_TYPE SetBC3D::ps_OBC_IsoThermal(REAL_TYPE* qbc, unsigned* bh1, int face, REAL_TYPE* t, REAL_TYPE* t0, REAL_TYPE& flop)
{
  SklParaManager* para_mng = ParaCmpo->GetParaManager();
  int i,j,k;
  unsigned register s, m;
  REAL_TYPE q, pp, sf, va=0.0;
  
  sf = FBUtility::convK2ND(obc[face].get_Temp(), BaseTemp, DiffTemp); // 保持されている温度はKelvin
  pp = (2.0*lambda) / (dh*RefV*rho*cp);
  
  switch (face) {
    case X_MINUS:
      if( pn.nID[face] < 0 ){
        i=1; // 最外層のID
        for (k=1; k<=kxc; k++) {
          for (j=1; j<=jxc; j++) {
            m = FBUtility::getFindexS3D(size, guide, i, j, k);
            s = bh1[m];
            q = pp * (sf - t0[m]);
            qbc[FBUtility::getFindexV3DEx(size, guide, 0, i-1, j  , k  )] += q; // マイナス面の正の値は流入
            va += q; 
            t[FBUtility::getFindexS3D(size, guide, i-1, j  , k  )] = sf; // 等温セルに指定温度を代入
          }
        }
      }
      flop += (REAL_TYPE)(jxc*kxc*4);
      break;
      
    case X_PLUS:
      if( pn.nID[face] < 0 ){
        i=ixc;
        for (k=1; k<=kxc; k++) {
          for (j=1; j<=jxc; j++) {
            m = FBUtility::getFindexS3D(size, guide, i, j, k);
            s = bh1[m];
            q = pp * (t0[m] - sf);
            qbc[FBUtility::getFindexV3DEx(size, guide, 0, i  , j  , k  )] -= q; // プラス面の正の値は流出
            va -= q; 
            t[FBUtility::getFindexS3D(size, guide, i+1, j  , k  )] = sf;
          }
        }
      }
      flop += (REAL_TYPE)(jxc*kxc*4);
      break;
      
    case Y_MINUS:
      if( pn.nID[face] < 0 ){
        j=1;
        for (k=1; k<=kxc; k++) {
          for (i=1; i<=ixc; i++) {
            m = FBUtility::getFindexS3D(size, guide, i, j, k);
            s = bh1[m];
            q = pp * (sf - t0[m]);
            qbc[FBUtility::getFindexV3DEx(size, guide, 1, i  , j-1, k  )] += q;
            va += q;
            t[FBUtility::getFindexS3D(size, guide, i  , j-1, k  )] = sf;
          }
        }
      }
      flop += (REAL_TYPE)(ixc*kxc*4);
      break;
      
    case Y_PLUS:
      if( pn.nID[face] < 0 ){
        j=jxc;
        for (k=1; k<=kxc; k++) {
          for (i=1; i<=ixc; i++) {
            m = FBUtility::getFindexS3D(size, guide, i, j, k);
            s = bh1[m];
            q = pp * (t0[m] - sf);
            qbc[FBUtility::getFindexV3DEx(size, guide, 1, i  , j  , k  )] -= q;
            va -= q;
            t[FBUtility::getFindexS3D(size, guide, i  , j+1, k  )] = sf;
          }
        }
      }
      flop += (REAL_TYPE)(ixc*kxc*4);
      break;
      
    case Z_MINUS:
      if( pn.nID[face] < 0 ){
        k=1;
        for (j=1; j<=jxc; j++) {
          for (i=1; i<=ixc; i++) {
            m = FBUtility::getFindexS3D(size, guide, i, j, k);
            s = bh1[m];
            q = pp * (sf - t0[m]);
            qbc[FBUtility::getFindexV3DEx(size, guide, 2, i  , j  , k-1)] += q;
            va += q;
            t[FBUtility::getFindexS3D(size, guide, i  , j  , k-1)] = sf;
          }
        }
      }
      flop += (REAL_TYPE)(ixc*jxc*4);
      break;
      
    case Z_PLUS:
      if( pn.nID[face] < 0 ){
        k=kxc;
        for (j=1; j<=jxc; j++) {
          for (i=1; i<=ixc; i++) {
            m = FBUtility::getFindexS3D(size, guide, i, j, k);
            s = bh1[m];
            q = pp * (t0[m] - sf);
            qbc[FBUtility::getFindexV3DEx(size, guide, 2, i  , j  , k  )] -= q;
            va -= q;
            t[FBUtility::getFindexS3D(size, guide, i  , j  , k+1)] = sf;
          }
        }
      }
      flop += (REAL_TYPE)(ixc*jxc*4);
      break;
  }
  
  if( para_mng->IsParallel() ){
		REAL_TYPE tmp = va;
		para_mng->Allreduce(&tmp, &va, 1, SKL_ARRAY_DTYPE_REAL, SKL_SUM, pn.procGrp);
	}
  
  return (va);
}

/**
 @fn REAL_TYPE SetBC3D::ps_OBC_Heatflux(REAL_TYPE* qbc, unsigned* bh1, int face, REAL_TYPE& flop)
 @brief 外部領域の熱流束指定の境界条件処理
 @retval 熱量(-)
 @param qbc 境界条件熱流束
 @param bh1 BCindex H1
 @param face 外部境界面番号
 @param flop 浮動小数演算数
 @note
    - モニタ量の熱量va(-)は系に対する流入出量
    - 流入を正にとる
 */
REAL_TYPE SetBC3D::ps_OBC_Heatflux(REAL_TYPE* qbc, unsigned* bh1, int face, REAL_TYPE& flop)
{
  SklParaManager* para_mng = ParaCmpo->GetParaManager();
  int i,j,k;
  unsigned register s, m;
  REAL_TYPE q, va=0.0;
  
  q = obc[face].get_Heatflux() / (RefV*DiffTemp*rho*cp); // [W/m^2]を無次元化
  
  switch (face) {
    case X_MINUS:
      if( pn.nID[face] < 0 ){
        i=1; // 最外層のID
        for (k=1; k<=kxc; k++) {
          for (j=1; j<=jxc; j++) {
            m = FBUtility::getFindexS3D(size, guide, i, j, k);
            s = bh1[m];
            qbc[FBUtility::getFindexV3DEx(size, guide, 0, i-1, j  , k  )] += q;
            va += q;
          }
        }
      }
      flop += (REAL_TYPE)(jxc*kxc*2);
      break;
      
    case X_PLUS:
      if( pn.nID[face] < 0 ){
        i=ixc;
        for (k=1; k<=kxc; k++) {
          for (j=1; j<=jxc; j++) {
            m = FBUtility::getFindexS3D(size, guide, i, j, k);
            s = bh1[m];
            qbc[FBUtility::getFindexV3DEx(size, guide, 0, i  , j  , k  )] -= q; // プラス面の正の値は流出
            va -= q;
          }
        }
      }
      flop += (REAL_TYPE)(jxc*kxc*2);
      break;
      
    case Y_MINUS:
      if( pn.nID[face] < 0 ){
        j=1;
        for (k=1; k<=kxc; k++) {
          for (i=1; i<=ixc; i++) {
            m = FBUtility::getFindexS3D(size, guide, i, j, k);
            s = bh1[m];
            qbc[FBUtility::getFindexV3DEx(size, guide, 1, i  , j-1, k  )] += q;
            va += q;
          }
        }
      }
      flop += (REAL_TYPE)(ixc*kxc*2);
      break;
      
    case Y_PLUS:
      if( pn.nID[face] < 0 ){
        j=jxc;
        for (k=1; k<=kxc; k++) {
          for (i=1; i<=ixc; i++) {
            m = FBUtility::getFindexS3D(size, guide, i, j, k);
            s = bh1[m];
            qbc[FBUtility::getFindexV3DEx(size, guide, 1, i  , j  , k  )] -= q;
            va -= q;
          }
        }
      }
      flop += (REAL_TYPE)(ixc*kxc*2);
      break;
      
    case Z_MINUS:
      if( pn.nID[face] < 0 ){
        k=1;
        for (j=1; j<=jxc; j++) {
          for (i=1; i<=ixc; i++) {
            m = FBUtility::getFindexS3D(size, guide, i, j, k);
            s = bh1[m];
            qbc[FBUtility::getFindexV3DEx(size, guide, 2, i  , j  , k-1)] += q;
            va += q;
          }
        }
      }
      flop += (REAL_TYPE)(ixc*jxc*2);
      break;
      
    case Z_PLUS:
      if( pn.nID[face] < 0 ){
        k=kxc;
        for (j=1; j<=jxc; j++) {
          for (i=1; i<=ixc; i++) {
            m = FBUtility::getFindexS3D(size, guide, i, j, k);
            s = bh1[m];
            qbc[FBUtility::getFindexV3DEx(size, guide, 2, i  , j  , k  )] -= q;
            va -= q;
          }
        }
      }
      flop += (REAL_TYPE)(ixc*jxc*2);
      break;
  }
  
  if( para_mng->IsParallel() ){
		REAL_TYPE tmp = va;
		para_mng->Allreduce(&tmp, &va, 1, SKL_ARRAY_DTYPE_REAL, SKL_SUM, pn.procGrp);
	}
  
  return (va);
}

/**
 @fn REAL_TYPE SetBC3D::ps_OBC_HeatTransfer_BS(REAL_TYPE* qbc, unsigned* bh1, int face, REAL_TYPE* t, REAL_TYPE* t0, REAL_TYPE& flop)
 @brief 外部領域の熱伝達境界の境界条件処理
 @retval 熱量(-)
 @param qbc 境界条件熱流束
 @param bh1 BCindex H1
 @param face 外部境界面番号
 @param t n+1時刻の温度場
 @param t0 n時刻の温度場
 @param flop 浮動小数演算数
 @note
    - モニタ量の熱量va(-)は系に対する流入出量
    - 流入を正にとる
    - TypeS/B共用
 */
REAL_TYPE SetBC3D::ps_OBC_HeatTransfer_BS(REAL_TYPE* qbc, unsigned* bh1, int face, REAL_TYPE* t, REAL_TYPE* t0, REAL_TYPE& flop)
{
  SklParaManager* para_mng = ParaCmpo->GetParaManager();
  int i,j,k;
  unsigned register s, m;
  REAL_TYPE ht, bt, q, va=0.0;
  
  ht = obc[face].get_CoefHT() / (RefV*rho*cp);  // 熱伝達係数
  bt = FBUtility::convK2ND(obc[face].get_Temp(), BaseTemp, DiffTemp); // 温度．保持されている温度はKelvin
  
  switch (face) {
    case X_MINUS:
      if( pn.nID[face] < 0 ){
        i=1; // 最外層のID
        for (k=1; k<=kxc; k++) {
          for (j=1; j<=jxc; j++) {
            m = FBUtility::getFindexS3D(size, guide, i, j, k);
            s = bh1[m];
            q = ht * (bt - t0[m]);
            qbc[FBUtility::getFindexV3DEx(size, guide, 0, i-1, j  , k  )] += q; // マイナス面の正の値は流入
            va += q; 
            t[FBUtility::getFindexS3D(size, guide, i-1, j  , k  )] = bt; // 隣接流体セルにバルク温度を代入
          }
        }
      }
      flop += (REAL_TYPE)(jxc*kxc*4);
      break;
      
    case X_PLUS:
      if( pn.nID[face] < 0 ){
        i=ixc;
        for (k=1; k<=kxc; k++) {
          for (j=1; j<=jxc; j++) {
            m = FBUtility::getFindexS3D(size, guide, i, j, k);
            s = bh1[m];
            q = ht * (t0[m] - bt);
            qbc[FBUtility::getFindexV3DEx(size, guide, 0, i  , j  , k  )] -= q; // プラス面の正の値は流出
            va -= q;
            t[FBUtility::getFindexS3D(size, guide, i+1, j  , k  )] = bt;
          }
        }
      }
      flop += (REAL_TYPE)(jxc*kxc*4);
      break;
      
    case Y_MINUS:
      if( pn.nID[face] < 0 ){
        j=1;
        for (k=1; k<=kxc; k++) {
          for (i=1; i<=ixc; i++) {
            m = FBUtility::getFindexS3D(size, guide, i, j, k);
            s = bh1[m];
            q = ht * (bt - t0[m]);
            qbc[FBUtility::getFindexV3DEx(size, guide, 1, i  , j-1, k  )] += q;
            va += q;
            t[FBUtility::getFindexS3D(size, guide, i  , j-1, k  )] = bt;
          }
        }
      }
      flop += (REAL_TYPE)(ixc*kxc*4);
      break;
      
    case Y_PLUS:
      if( pn.nID[face] < 0 ){
        j=jxc;
        for (k=1; k<=kxc; k++) {
          for (i=1; i<=ixc; i++) {
            m = FBUtility::getFindexS3D(size, guide, i, j, k);
            s = bh1[m];
            q = ht * (t0[m] - bt);
            qbc[FBUtility::getFindexV3DEx(size, guide, 1, i  , j  , k  )] -= q;
            va -= q;
            t[FBUtility::getFindexS3D(size, guide, i  , j+1, k  )] = bt;
          }
        }
      }
      flop += (REAL_TYPE)(ixc*kxc*4);
      break;
      
    case Z_MINUS:
      if( pn.nID[face] < 0 ){
        k=1;
        for (j=1; j<=jxc; j++) {
          for (i=1; i<=ixc; i++) {
            m = FBUtility::getFindexS3D(size, guide, i, j, k);
            s = bh1[m];
            q = ht * (bt - t0[m]);
            qbc[FBUtility::getFindexV3DEx(size, guide, 2, i  , j  , k-1)] += q;
            va += q;
            t[FBUtility::getFindexS3D(size, guide, i  , j  , k-1)] = bt;
          }
        }
      }
      flop += (REAL_TYPE)(ixc*jxc*4);
      break;
      
    case Z_PLUS:
      if( pn.nID[face] < 0 ){
        k=kxc;
        for (j=1; j<=jxc; j++) {
          for (i=1; i<=ixc; i++) {
            m = FBUtility::getFindexS3D(size, guide, i, j, k);
            s = bh1[m];
            q = ht * (t0[m] - bt);
            qbc[FBUtility::getFindexV3DEx(size, guide, 2, i  , j  , k  )] -= q;
            va -= q;
            t[FBUtility::getFindexS3D(size, guide, i  , j  , k+1)] = bt;
          }
        }
      }
      flop += (REAL_TYPE)(ixc*jxc*4);
      break;
  }
  
  if( para_mng->IsParallel() ){
		REAL_TYPE tmp = va;
		para_mng->Allreduce(&tmp, &va, 1, SKL_ARRAY_DTYPE_REAL, SKL_SUM, pn.procGrp);
	}
  
  return (va);
}

/**
 @fn REAL_TYPE SetBC3D::ps_OBC_HeatTransfer_SF(REAL_TYPE* qbc, unsigned* bh1, int face, REAL_TYPE* t, REAL_TYPE* t0, REAL_TYPE& flop)
 @brief 外部領域の熱伝達境界の境界条件処理
 @retval 熱量(-)
 @param qbc 境界条件熱流束
 @param bh1 BCindex H1
 @param face 外部境界面番号
 @param t n+1時刻の温度場
 @param t0 n時刻の温度場
 @param flop 浮動小数演算数
 @note
    - モニタ量の熱量va(-)は系に対する流入出量
    - 流入を正にとる
 */
REAL_TYPE SetBC3D::ps_OBC_HeatTransfer_SF(REAL_TYPE* qbc, unsigned* bh1, int face, REAL_TYPE* t, REAL_TYPE* t0, REAL_TYPE& flop)
{
  SklParaManager* para_mng = ParaCmpo->GetParaManager();
  int i,j,k;
  unsigned register s, m;
  REAL_TYPE ht, sf, q, va=0.0;
  REAL_TYPE a1, b1, c1;
  
  a1 = obc[face].ca[CompoList::alpha];
  b1 = obc[face].ca[CompoList::beta];
  c1 = obc[face].ca[CompoList::gamma];
  ht = lambda / (RefV*RefL*rho*cp) * a1*pow(Reynolds, b1) * pow(Prandtl, c1);    // 熱伝達係数 
  sf = FBUtility::convK2ND(obc[face].get_Temp(), BaseTemp, DiffTemp); // 表面温度．保持されている温度はKelvin
  
  switch (face) {
    case X_MINUS:
      if( pn.nID[face] < 0 ){
        i=1; // 最外層のID
        if ( obc[face].get_HTmodeRef() == CompoList::HT_mode_bulk ) {
          for (k=1; k<=kxc; k++) {
            for (j=1; j<=jxc; j++) {
              m = FBUtility::getFindexS3D(size, guide, i, j, k);
              s = bh1[m];
              q = ht * (sf - 0.0); // 基準温度との温度差
              qbc[FBUtility::getFindexV3DEx(size, guide, 0, i-1, j  , k  )] += q; // マイナス面の正の値は流入
              va += q; // マイナス面の正の値は流入
              t[FBUtility::getFindexS3D(size, guide, i-1, j  , k  )] = sf; // 固体セルに表面温度を代入
            }
          }
        }
        else { // Local
          for (k=1; k<=kxc; k++) {
            for (j=1; j<=jxc; j++) {
              m = FBUtility::getFindexS3D(size, guide, i, j, k);
              s = bh1[m];
              q = ht * (sf - t0[m]); // 基準温度との温度差
              qbc[FBUtility::getFindexV3DEx(size, guide, 0, i-1, j  , k  )] += q;
              va += q;
              t[FBUtility::getFindexS3D(size, guide, i-1, j  , k  )] = sf;
            }
          }          
        }
      }
      flop += (REAL_TYPE)(jxc*kxc*4);
      break;
      
    case X_PLUS:
      if( pn.nID[face] < 0 ){
        i=ixc;
        if ( obc[face].get_HTmodeRef() == CompoList::HT_mode_bulk ) {
          for (k=1; k<=kxc; k++) {
            for (j=1; j<=jxc; j++) {
              m = FBUtility::getFindexS3D(size, guide, i, j, k);
              s = bh1[m];
              q = ht * (0.0 - sf);
              qbc[FBUtility::getFindexV3DEx(size, guide, 0, i  , j  , k  )] -= q; // プラス面の正の値は流出
              va -= q; // プラス面の正の値は流出
              t[FBUtility::getFindexS3D(size, guide, i+1, j  , k  )] = sf;
            }
          }          
        }
        else { // Local
          for (k=1; k<=kxc; k++) {
            for (j=1; j<=jxc; j++) {
              m = FBUtility::getFindexS3D(size, guide, i, j, k);
              s = bh1[m];
              q = ht * (t0[m] - sf);
              qbc[FBUtility::getFindexV3DEx(size, guide, 0, i  , j  , k  )] -= q;
              va -= q; // プラス面の正の値は流出
              t[FBUtility::getFindexS3D(size, guide, i+1, j  , k  )] = sf;
            }
          }          
        }
      }
      flop += (REAL_TYPE)(jxc*kxc*4);
      break;
      
    case Y_MINUS:
      if( pn.nID[face] < 0 ){
        j=1;
        if ( obc[face].get_HTmodeRef() == CompoList::HT_mode_bulk ) {
          for (k=1; k<=kxc; k++) {
            for (i=1; i<=ixc; i++) {
              m = FBUtility::getFindexS3D(size, guide, i, j, k);
              s = bh1[m];
              q = ht * (sf - 0.0);
              qbc[FBUtility::getFindexV3DEx(size, guide, 1, i  , j-1, k  )] += q;
              va += q;
              t[FBUtility::getFindexS3D(size, guide, i  , j-1, k  )] = sf;
            }
          }          
        }
        else { // Local
          for (k=1; k<=kxc; k++) {
            for (i=1; i<=ixc; i++) {
              m = FBUtility::getFindexS3D(size, guide, i, j, k);
              s = bh1[m];
              q = ht * (sf - t0[m]);
              qbc[FBUtility::getFindexV3DEx(size, guide, 1, i  , j-1, k  )] += q;
              va += q;
              t[FBUtility::getFindexS3D(size, guide, i  , j-1, k  )] = sf;
            }
          }          
        }
      }
      flop += (REAL_TYPE)(ixc*kxc*4);
      break;
      
    case Y_PLUS:
      if( pn.nID[face] < 0 ){
        j=jxc;
        if ( obc[face].get_HTmodeRef() == CompoList::HT_mode_bulk ) {
          for (k=1; k<=kxc; k++) {
            for (i=1; i<=ixc; i++) {
              m = FBUtility::getFindexS3D(size, guide, i, j, k);
              s = bh1[m];
              q = ht * (0.0 - sf);
              qbc[FBUtility::getFindexV3DEx(size, guide, 1, i  , j  , k  )] -= q;
              va -= q;
              t[FBUtility::getFindexS3D(size, guide, i  , j+1, k  )] = sf;
            }
          }          
        }
        else { // Local
          for (k=1; k<=kxc; k++) {
            for (i=1; i<=ixc; i++) {
              m = FBUtility::getFindexS3D(size, guide, i, j, k);
              s = bh1[m];
              q = ht * (t0[m] - sf);
              qbc[FBUtility::getFindexV3DEx(size, guide, 1, i  , j  , k  )] -= q;
              va -= q;
              t[FBUtility::getFindexS3D(size, guide, i  , j+1, k  )] = sf;
            }
          }          
        }
      }
      flop += (REAL_TYPE)(ixc*kxc*4);
      break;
      
    case Z_MINUS:
      if( pn.nID[face] < 0 ){
        k=1;
        if ( obc[face].get_HTmodeRef() == CompoList::HT_mode_bulk ) {
          for (j=1; j<=jxc; j++) {
            for (i=1; i<=ixc; i++) {
              m = FBUtility::getFindexS3D(size, guide, i, j, k);
              s = bh1[m];
              q = ht * (sf - 0.0);
              qbc[FBUtility::getFindexV3DEx(size, guide, 2, i  , j  , k-1)] += q;
              va += q;
              t[FBUtility::getFindexS3D(size, guide, i  , j  , k-1)] = sf;              
            }
          }
        }
        else { // Local
          for (j=1; j<=jxc; j++) {
            for (i=1; i<=ixc; i++) {
              m = FBUtility::getFindexS3D(size, guide, i, j, k);
              s = bh1[m];
              q = ht * (sf - t0[m]);
              qbc[FBUtility::getFindexV3DEx(size, guide, 2, i  , j  , k-1)] += q;
              va += q;
              t[FBUtility::getFindexS3D(size, guide, i  , j  , k-1)] = sf;
            }
          }
        }
      }
      flop += (REAL_TYPE)(ixc*jxc*4);
      break;
      
    case Z_PLUS:
      if( pn.nID[face] < 0 ){
        k=kxc;
        if ( obc[face].get_HTmodeRef() == CompoList::HT_mode_bulk ) {
          for (j=1; j<=jxc; j++) {
            for (i=1; i<=ixc; i++) {
              m = FBUtility::getFindexS3D(size, guide, i, j, k);
              s = bh1[m];
              q = ht * (0.0 - sf);
              qbc[FBUtility::getFindexV3DEx(size, guide, 2, i  , j  , k  )] -= q;
              va -= q;
              t[FBUtility::getFindexS3D(size, guide, i  , j  , k+1)] = sf;              
            }
          }
        }
        else { // Local
          for (j=1; j<=jxc; j++) {
            for (i=1; i<=ixc; i++) {
              m = FBUtility::getFindexS3D(size, guide, i, j, k);
              s = bh1[m];
              q = ht * (t0[m] - sf);
              qbc[FBUtility::getFindexV3DEx(size, guide, 2, i  , j  , k  )] -= q;
              va -= q;
              t[FBUtility::getFindexS3D(size, guide, i  , j  , k+1)] = sf;
            }
          }
        }
      }
      flop += (REAL_TYPE)(ixc*jxc*4);
      break;
  }
  
  if( para_mng->IsParallel() ){
		REAL_TYPE tmp = va;
		para_mng->Allreduce(&tmp, &va, 1, SKL_ARRAY_DTYPE_REAL, SKL_SUM, pn.procGrp);
	}
  
  return (va);
}

/**
 @fn REAL_TYPE SetBC3D::ps_OBC_HeatTransfer_SN(REAL_TYPE* qbc, unsigned* bh1, int face, REAL_TYPE* t, REAL_TYPE* t0, REAL_TYPE& flop)
 @brief 外部領域の熱伝達境界の境界条件処理
 @retval 熱量(-)
 @param qbc 境界条件熱流束
 @param bh1 BCindex H1
 @param face 外部境界面番号
 @param t n+1時刻の温度場
 @param t0 n時刻の温度場
 @param flop 浮動小数演算数
 @note
 - モニタ量の熱量va(-)は系に対する流入出量
 - 流入を正にとる
 */
REAL_TYPE SetBC3D::ps_OBC_HeatTransfer_SN(REAL_TYPE* qbc, unsigned* bh1, int face, REAL_TYPE* t, REAL_TYPE* t0, REAL_TYPE& flop)
{
  SklParaManager* para_mng = ParaCmpo->GetParaManager();
  int i,j,k;
  unsigned register s, m;
  REAL_TYPE ht, ht1, ht3, sf, q, va=0.0;
  REAL_TYPE a1, a2, a3, a4, b1, b2, b3, b4, c1, c2;
  
  a1 = obc[face].ca[CompoList::vert_laminar_alpha];    // vertical, upper, laminar
  a2 = obc[face].ca[CompoList::vert_turbulent_alpha];  // vertical, upper, turbulent
  a3 = obc[face].cb[CompoList::lower_laminar_alpha];   // lower, laminar
  a4 = obc[face].cb[CompoList::lower_turbulent_alpha]; // lower, turbulent
  
  b1 = obc[face].ca[CompoList::vert_laminar_beta];     // vertical, upper, laminar
  b2 = obc[face].ca[CompoList::vert_turbulent_beta];   // vertical, upper, turbulent
  b3 = obc[face].cb[CompoList::lower_laminar_beta];    // lower, laminar
  b4 = obc[face].cb[CompoList::lower_turbulent_beta];  // lower, turbulent
  
  c1 = obc[face].ca[CompoList::vert_Ra_critial];       // vertical, upper, Ra_c
  c2 = obc[face].cb[CompoList::lower_Ra_critial];      // lower, Ra_c
  ht = lambda / (RefV*RefL*rho*cp);      // Ref_IDで指定される物性値
  
  ht1 = ( Rayleigh < c1 ) ? ht*a1*pow(Rayleigh, b1) : ht*a2*pow(Rayleigh, b2); // vertical, upper
  ht3 = ( Rayleigh < c2 ) ? ht*a3*pow(Rayleigh, b3) : ht*a4*pow(Rayleigh, b4); // lower
  sf = FBUtility::convK2ND(obc[face].get_Temp(), BaseTemp, DiffTemp); // 保持されている温度はKelvin
  
  switch (face) {
    case X_MINUS:
      if( pn.nID[face] < 0 ){
        i=1; // 最外層のID
        if ( obc[face].get_HTmodeRef() == CompoList::HT_mode_bulk ) {
          for (k=1; k<=kxc; k++) {
            for (j=1; j<=jxc; j++) {
              m = FBUtility::getFindexS3D(size, guide, i, j, k);
              s = bh1[m];
              q = ht1 * (sf - 0.0); // 基準温度との温度差
              qbc[FBUtility::getFindexV3DEx(size, guide, 0, i-1, j  , k  )] += q; // マイナス面の正の値は流入
              va += q; // マイナス面の正の値は流入
              t[FBUtility::getFindexS3D(size, guide, i-1, j  , k  )] = sf; // 固体セルに表面温度を代入
            }
          }
        }
        else { // Local
          for (k=1; k<=kxc; k++) {
            for (j=1; j<=jxc; j++) {
              m = FBUtility::getFindexS3D(size, guide, i, j, k);
              s = bh1[m];
              q = ht1 * (sf - t0[m]); // 基準温度との温度差
              qbc[FBUtility::getFindexV3DEx(size, guide, 0, i-1, j  , k  )] += q;
              va += q;
              t[FBUtility::getFindexS3D(size, guide, i-1, j  , k  )] = sf;
            }
          }          
        }
      }
      flop += (REAL_TYPE)(jxc*kxc*4);
      break;
      
    case X_PLUS:
      if( pn.nID[face] < 0 ){
        i=ixc;
        if ( obc[face].get_HTmodeRef() == CompoList::HT_mode_bulk ) {
          for (k=1; k<=kxc; k++) {
            for (j=1; j<=jxc; j++) {
              m = FBUtility::getFindexS3D(size, guide, i, j, k);
              s = bh1[m];
              q = ht1 * (0.0 - sf);
              qbc[FBUtility::getFindexV3DEx(size, guide, 0, i  , j  , k  )] -= q; // プラス面の正の値は流出
              va -= q; // プラス面の正の値は流出
              t[FBUtility::getFindexS3D(size, guide, i+1, j  , k  )] = sf;
            }
          }          
        }
        else { // Local
          for (k=1; k<=kxc; k++) {
            for (j=1; j<=jxc; j++) {
              m = FBUtility::getFindexS3D(size, guide, i, j, k);
              s = bh1[m];
              q = ht1 * (t0[m] - sf);
              qbc[FBUtility::getFindexV3DEx(size, guide, 0, i  , j  , k  )] -= q;
              va -= q; // プラス面の正の値は流出
              t[FBUtility::getFindexS3D(size, guide, i+1, j  , k  )] = sf;
            }
          }          
        }
      }
      flop += (REAL_TYPE)(jxc*kxc*4);
      break;
      
    case Y_MINUS:
      if( pn.nID[face] < 0 ){
        j=1;
        if ( obc[face].get_HTmodeRef() == CompoList::HT_mode_bulk ) {
          for (k=1; k<=kxc; k++) {
            for (i=1; i<=ixc; i++) {
              m = FBUtility::getFindexS3D(size, guide, i, j, k);
              s = bh1[m];
              q = ht1 * (sf - 0.0);
              qbc[FBUtility::getFindexV3DEx(size, guide, 1, i  , j-1, k  )] += q;
              va += q;
              t[FBUtility::getFindexS3D(size, guide, i  , j-1, k  )] = sf;
            }
          }          
        }
        else { // Local
          for (k=1; k<=kxc; k++) {
            for (i=1; i<=ixc; i++) {
              m = FBUtility::getFindexS3D(size, guide, i, j, k);
              s = bh1[m];
              q = ht1 * (sf - t0[m]);
              qbc[FBUtility::getFindexV3DEx(size, guide, 1, i  , j-1, k  )] += q;
              va += q;
              t[FBUtility::getFindexS3D(size, guide, i  , j-1, k  )] = sf;
            }
          }          
        }
      }
      flop += (REAL_TYPE)(ixc*kxc*4);
      break;
      
    case Y_PLUS:
      if( pn.nID[face] < 0 ){
        j=jxc;
        if ( obc[face].get_HTmodeRef() == CompoList::HT_mode_bulk ) {
          for (k=1; k<=kxc; k++) {
            for (i=1; i<=ixc; i++) {
              m = FBUtility::getFindexS3D(size, guide, i, j, k);
              s = bh1[m];
              q = ht1 * (0.0 - sf);
              qbc[FBUtility::getFindexV3DEx(size, guide, 1, i  , j  , k  )] -= q;
              va -= q;
              t[FBUtility::getFindexS3D(size, guide, i  , j+1, k  )] = sf;
            }
          }          
        }
        else { // Local
          for (k=1; k<=kxc; k++) {
            for (i=1; i<=ixc; i++) {
              m = FBUtility::getFindexS3D(size, guide, i, j, k);
              s = bh1[m];
              q = ht1 * (t0[m] - sf);
              qbc[FBUtility::getFindexV3DEx(size, guide, 1, i  , j  , k  )] -= q;
              va -= q;
              t[FBUtility::getFindexS3D(size, guide, i  , j+1, k  )] = sf;
            }
          }          
        }
      }
      flop += (REAL_TYPE)(ixc*kxc*4);
      break;
      
    case Z_MINUS:
      if( pn.nID[face] < 0 ){
        k=1;
        if ( obc[face].get_HTmodeRef() == CompoList::HT_mode_bulk ) {
          for (j=1; j<=jxc; j++) {
            for (i=1; i<=ixc; i++) {
              m = FBUtility::getFindexS3D(size, guide, i, j, k);
              s = bh1[m];
              q = ht3 * (sf - 0.0);
              qbc[FBUtility::getFindexV3DEx(size, guide, 2, i  , j  , k-1)] += q;
              va += q;
              t[FBUtility::getFindexS3D(size, guide, i  , j  , k-1)] = sf;              
            }
          }
        }
        else { // Local
          for (j=1; j<=jxc; j++) {
            for (i=1; i<=ixc; i++) {
              m = FBUtility::getFindexS3D(size, guide, i, j, k);
              s = bh1[m];
              q = ht3 * (sf - t0[m]);
              qbc[FBUtility::getFindexV3DEx(size, guide, 2, i  , j  , k-1)] += q;
              va += q;
              t[FBUtility::getFindexS3D(size, guide, i  , j  , k-1)] = sf;
            }
          }
        }
      }
      flop += (REAL_TYPE)(ixc*jxc*4);
      break;
      
    case Z_PLUS:
      if( pn.nID[face] < 0 ){
        k=kxc;
        if ( obc[face].get_HTmodeRef() == CompoList::HT_mode_bulk ) {
          for (j=1; j<=jxc; j++) {
            for (i=1; i<=ixc; i++) {
              m = FBUtility::getFindexS3D(size, guide, i, j, k);
              s = bh1[m];
              q = ht1 * (0.0 - sf);
              qbc[FBUtility::getFindexV3DEx(size, guide, 2, i  , j  , k  )] -= q;
              va -= q;
              t[FBUtility::getFindexS3D(size, guide, i  , j  , k+1)] = sf;              
            }
          }
        }
        else { // Local
          for (j=1; j<=jxc; j++) {
            for (i=1; i<=ixc; i++) {
              m = FBUtility::getFindexS3D(size, guide, i, j, k);
              s = bh1[m];
              q = ht1 * (t0[m] - sf);
              qbc[FBUtility::getFindexV3DEx(size, guide, 2, i  , j  , k  )] -= q;
              va -= q;
              t[FBUtility::getFindexS3D(size, guide, i  , j  , k+1)] = sf;
            }
          }
        }
      }
      flop += (REAL_TYPE)(ixc*jxc*4);
      break;
  }
  
  if( para_mng->IsParallel() ){
		REAL_TYPE tmp = va;
		para_mng->Allreduce(&tmp, &va, 1, SKL_ARRAY_DTYPE_REAL, SKL_SUM, pn.procGrp);
	}
  
  return (va);
}

/**
 @fn REAL_TYPE SetBC3D::ps_IBC_Heatflux(REAL_TYPE* qbc, unsigned* bh1, int n, REAL_TYPE& flop)
 @brief 内部領域の熱流束指定境界条件
 @retval 熱量(W)
 @param qbc 境界条件熱流束
 @param bh1 BCindex H1
 @param n コンポーネントリストのインデクス
 @param flop 浮動小数演算数
 @note
    - モニタ量の熱量vaは系に対する流入量なので，基準温度に対する熱量
    - 流入量を正にとる
 */
REAL_TYPE SetBC3D::ps_IBC_Heatflux(REAL_TYPE* qbc, unsigned* bh1, int n, REAL_TYPE& flop)
{
  SklParaManager* para_mng = ParaCmpo->GetParaManager();
  int i,j,k;
  unsigned odr;
  unsigned register s, m;
  REAL_TYPE va=0.0, q;
  int st[3], ed[3];
  
  odr= cmp[n].getMatOdr();
  q = cmp[n].get_Heatflux() / (RefV*DiffTemp*rho*cp); // [W/m^2]を無次元化

  cmp[n].getBbox(st, ed);
	
  for (k=st[2]; k<=ed[2]; k++) {
    for (j=st[1]; j<=ed[1]; j++) {
      for (i=st[0]; i<=ed[0]; i++) {

        m = FBUtility::getFindexS3D(size, guide, i, j, k);
        s = bh1[m];
        
        if ( GET_FACE_BC(s, BC_FACE_W) == n ) {
          qbc[FBUtility::getFindexV3DEx(size, guide, 0, i-1, j  , k  )] += q;
          va += q; // マイナス面の正の値は流入
          flop += 2.0;
        }
        
        if ( GET_FACE_BC(s, BC_FACE_E) == n ) {
          qbc[FBUtility::getFindexV3DEx(size, guide, 0, i  , j  , k  )] -= q;
          va -= q; // プラス面の正の値は流出
          flop += 2.0;
        }
        
        if ( GET_FACE_BC(s, BC_FACE_S) == n ) {
          qbc[FBUtility::getFindexV3DEx(size, guide, 1, i  , j-1, k  )] += q;
          va += q;
          flop += 2.0;
        }
        
        if ( GET_FACE_BC(s, BC_FACE_N) == n ) {
          qbc[FBUtility::getFindexV3DEx(size, guide, 1, i  , j  , k  )] -= q;
          va -= q;
          flop += 2.0;
        }

        if ( GET_FACE_BC(s, BC_FACE_B) == n ) {
          qbc[FBUtility::getFindexV3DEx(size, guide, 2, i  , j  , k-1)] += q;
          va += q;
          flop += 1.0;
        }
        
        if ( GET_FACE_BC(s, BC_FACE_T) == n ) {
          qbc[FBUtility::getFindexV3DEx(size, guide, 2, i  , j  , k  )] -= q;
          va -= q;
          flop += 2.0;
        }
      }
    }
  }
  
  if( para_mng->IsParallel() ){
		REAL_TYPE tmp = va;
		para_mng->Allreduce(&tmp, &va, 1, SKL_ARRAY_DTYPE_REAL, SKL_SUM, pn.procGrp);
	}
  
  return (va); // 単位面積あたりの無次元熱流束の和
}

/**
 @fn REAL_TYPE SetBC3D::setHeatTransferN_SM(REAL_TYPE* qbc, REAL_TYPE* t, unsigned* bx, unsigned n, REAL_TYPE* t0)
 @brief 単媒質の場合の熱伝達境界条件タイプN　共役熱移動の場合
 @retval sum of heat flux (W/m^2)
 @param qbc 境界条件熱流束
 @param t n+1時刻の温度場
 @param bx BCindex
 @param n 境界条件コンポーネントのエントリ番号
 @param t0 n時刻の温度場
 @todo qsumの方向チェック
 
REAL_TYPE SetBC3D::setHeatTransferN_SM(REAL_TYPE* qbc, REAL_TYPE* t, unsigned* bx, unsigned n, REAL_TYPE* t0)
{
  SklParaManager* para_mng = ParaCmpo->GetParaManager();
  int i,j,k;
  unsigned s, odr, d;
  unsigned m0, mi, mj, mk;
  REAL_TYPE ht, qsum, tmp, c;
  REAL_TYPE qi, qj, qk;

  qsum = 0.0;
  odr= cmp[n].getMatOdr();
  //ht = cmp[n].get_CoefHT() / (RefV*mat[odr].P[p_density]*mat[odr].P[p_specific_heat]);
  ht = cmp[n].get_CoefHT() / (RefV*rho*cp);

  for (k=(int)cmp[n].ci.st[2]-1; k<=(int)cmp[n].ci.ed[2]; k++) {
    for (j=(int)cmp[n].ci.st[1]-1; j<=(int)cmp[n].ci.ed[1]; j++) {
      for (i=(int)cmp[n].ci.st[0]-1; i<=(int)cmp[n].ci.ed[0]; i++) {
        m0 = FBUtility::getFindexS3D(size, guide, i  , j  , k  );
        mi = FBUtility::getFindexS3D(size, guide, i+1, j  , k  );
        mj = FBUtility::getFindexS3D(size, guide, i  , j+1, k  );
        mk = FBUtility::getFindexS3D(size, guide, i  , j  , k+1);
        s = bx[m0];
        if ( (s & BIT_MASK_10) == n ) {
          if ( BIT_IS_SHIFT(s, QFACE_I) ) {
            qi = ht*(t[m0]-t[mi]);
            qbc[FBUtility::getFindexV3DEx(size, guide, 0, i, j, k)] = qi;
            d = BIT_SHIFT(s, DIR_I);
            qsum += qi * ( (0==d) ? 1.0 : -1.0 );
          }
          if ( BIT_IS_SHIFT(s, QFACE_J) ) {
            qj = ht*(t[m0]-t[mj]);
            qbc[FBUtility::getFindexV3DEx(size, guide, 1, i, j, k)] = qj;
            d = BIT_SHIFT(s, DIR_J);
            qsum += qj * ( (0==d) ? 1.0 : -1.0 );
          }
          if ( BIT_IS_SHIFT(s, QFACE_K) ) {
            qk = ht*(t[m0]-t[mk]);
            qbc[FBUtility::getFindexV3DEx(size, guide, 2, i, j, k)] = qk;
            d = BIT_SHIFT(s, DIR_K);
            qsum += qk * ( (0==d) ? 1.0 : -1.0 );
          }
        }
      }
    }
  }
  
  if( para_mng->IsParallel() ){
		tmp = qsum;
		para_mng->Allreduce(&tmp, &qsum, 1, SKL_ARRAY_DTYPE_REAL, SKL_SUM, pn.procGrp);
	}
  c = qsum * (RefV*DiffTemp*rho*cp);
  
  return (c);
}*/

/**
 @fn REAL_TYPE SetBC3D::ps_IBC_Transfer_S_SM(REAL_TYPE* qbc, unsigned* bh1, int n, REAL_TYPE* t, REAL_TYPE* t0, REAL_TYPE& flop)
 @brief 単媒質の場合の熱伝達温境界条件タイプS　流体の流動のみを解く場合
 @retval 熱流束の和 (W/m^2)
 @param qbc 境界条件熱流束
 @param bh1 BCindex H1
 @param n 境界条件コンポーネントのエントリ番号
 @param t n+1時刻の温度場
 @param t0 n時刻の温度場
 @param flop 浮動小数演算数
 @note 
    - 熱流束は加算（他の条件と合成）
    - pセルは流体セル
 */
REAL_TYPE SetBC3D::ps_IBC_Transfer_S_SM(REAL_TYPE* qbc, unsigned* bh1, int n, REAL_TYPE* t, REAL_TYPE* t0, REAL_TYPE& flop)
{
  SklParaManager* para_mng = ParaCmpo->GetParaManager();
  int i,j,k;
  unsigned s, odr, m;
  REAL_TYPE ht, sf, va=0.0, q;
  int st[3], ed[3];

  odr= cmp[n].getMatOdr();
  //ht = cmp[n].get_CoefHT() / (RefV*mat[odr].P[p_density]*mat[odr].P[p_specific_heat]); // 固体セルで指定される物性値
  ht = cmp[n].get_CoefHT() / (RefV*rho*cp); // Ref_IDで指定される物性値
  sf = FBUtility::convK2ND(cmp[n].get_Temp(), BaseTemp, DiffTemp); // 保持されている温度はKelvin

  cmp[n].getBbox(st, ed);
  
  for (k=st[2]; k<=ed[2]; k++) {
    for (j=st[1]; j<=ed[1]; j++) {
      for (i=st[0]; i<=ed[0]; i++) {
        
        m = FBUtility::getFindexS3D(size, guide, i  , j  , k  );
        s = bh1[m];
        
        if ( GET_FACE_BC(s, BC_FACE_W) == n ) {
          q = ht * (sf - t0[m]);
          qbc[FBUtility::getFindexV3DEx(size, guide, 0, i-1, j  , k  )] += q;
          va += q; // マイナス面の正の値は流入
          t[FBUtility::getFindexS3D(size, guide, i-1, j  , k  )] = sf; // 表面温度を代入
          flop += 4.0;
        }
        
        if ( GET_FACE_BC(s, BC_FACE_E) == n ) {
          q = ht * (t0[m] - sf);
          qbc[FBUtility::getFindexV3DEx(size, guide, 0, i  , j  , k  )] += q;
          va -= q; // プラス面の正の値は流出
          t[FBUtility::getFindexS3D(size, guide, i+1, j  , k  )] = sf;
          flop += 4.0;
        }
        
        if ( GET_FACE_BC(s, BC_FACE_S) == n ) {
          q = ht * (sf - t0[m]);
          qbc[FBUtility::getFindexV3DEx(size, guide, 1, i  , j-1, k  )] += q;
          va += q;
          t[FBUtility::getFindexS3D(size, guide, i  , j-1, k  )] = sf;
          flop += 4.0;
        }
        
        if ( GET_FACE_BC(s, BC_FACE_N) == n ) {
          q = ht * (t0[m] - sf);
          qbc[FBUtility::getFindexV3DEx(size, guide, 1, i  , j  , k  )] += q;
          va -= q;
          t[FBUtility::getFindexS3D(size, guide, i  , j+1, k  )] = sf;
          flop += 4.0;
        }
        
        if ( GET_FACE_BC(s, BC_FACE_B) == n ) {
          q = ht * (sf - t0[m]);
          qbc[FBUtility::getFindexV3DEx(size, guide, 2, i  , j  , k-1)] += q;
          va += q;
          t[FBUtility::getFindexS3D(size, guide, i  , j  , k-1)] = sf;
          flop += 4.0;
        }
        
        if ( GET_FACE_BC(s, BC_FACE_T) == n ) {
          q = ht * (t0[m] - sf);
          qbc[FBUtility::getFindexV3DEx(size, guide, 2, i  , j  , k  )] += q;
          va -= q;
          t[FBUtility::getFindexS3D(size, guide, i  , j  , k+1)] = sf;
          flop += 4.0;
        }
      }
    }
  }
  
  if( para_mng->IsParallel() ){
		REAL_TYPE tmp = va;
		para_mng->Allreduce(&tmp, &va, 1, SKL_ARRAY_DTYPE_REAL, SKL_SUM, pn.procGrp);
	}
  
  return ( va ); // 単位面積あたりの無次元熱流束の和
}

/**
 @fn REAL_TYPE SetBC3D::ps_IBC_Transfer_SN_SM(REAL_TYPE* qbc, unsigned* bh1, int n, REAL_TYPE* t, REAL_TYPE* t0, REAL_TYPE& flop)
 @brief 単媒質の場合の熱伝達温境界条件タイプSN（自然対流）　流体のみを解く場合
 @retval 熱流束の和 (W/m^2)
 @param qbc 境界条件熱流束
 @param bh1 BCindex H1
 @param n 境界条件コンポーネントのエントリ番号
 @param t n+1時刻の温度場
 @param t0 n時刻の温度場
 @param flop 浮動小数演算数
 @note 
    - 熱流束は加算（他の条件と合成）
    - pセルは流体セル
 */
REAL_TYPE SetBC3D::ps_IBC_Transfer_SN_SM(REAL_TYPE* qbc, unsigned* bh1, int n, REAL_TYPE* t, REAL_TYPE* t0, REAL_TYPE& flop)
{
  SklParaManager* para_mng = ParaCmpo->GetParaManager();
  int i,j,k;
  unsigned s, m;
  REAL_TYPE ht, ht1, ht3, sf, q, va=0.0;
  REAL_TYPE a1, a2, a3, a4, b1, b2, b3, b4, c1, c2;
  int st[3], ed[3];
  
  a1 = cmp[n].ca[CompoList::vert_laminar_alpha];    // vertical, upper, laminar
  a2 = cmp[n].ca[CompoList::vert_turbulent_alpha];  // vertical, upper, turbulent
  a3 = cmp[n].cb[CompoList::lower_laminar_alpha];   // lower, laminar
  a4 = cmp[n].cb[CompoList::lower_turbulent_alpha]; // lower, turbulent
  
  b1 = cmp[n].ca[CompoList::vert_laminar_beta];     // vertical, upper, laminar
  b2 = cmp[n].ca[CompoList::vert_turbulent_beta];   // vertical, upper, turbulent
  b3 = cmp[n].cb[CompoList::lower_laminar_beta];    // lower, laminar
  b4 = cmp[n].cb[CompoList::lower_turbulent_beta];  // lower, turbulent
  
  c1 = cmp[n].ca[CompoList::vert_Ra_critial];       // vertical, upper, Ra_c
  c2 = cmp[n].cb[CompoList::lower_Ra_critial];      // lower, Ra_c
  ht = lambda / (RefV*RefL*rho*cp);      // Ref_IDで指定される物性値
  
  ht1 = ( Rayleigh < c1 ) ? ht*a1*pow(Rayleigh, b1) : ht*a2*pow(Rayleigh, b2); // vertical, upper
  ht3 = ( Rayleigh < c2 ) ? ht*a3*pow(Rayleigh, b3) : ht*a4*pow(Rayleigh, b4); // lower
  sf = FBUtility::convK2ND(cmp[n].get_Temp(), BaseTemp, DiffTemp); // 保持されている温度はKelvin
  
  cmp[n].getBbox(st, ed);
  
  if ( cmp[n].get_sw_HTmodeRef() == CompoList::HT_mode_bulk ) {
    for (k=st[2]; k<=ed[2]; k++) {
      for (j=st[1]; j<=ed[1]; j++) {
        for (i=st[0]; i<=ed[0]; i++) {
          
          m = FBUtility::getFindexS3D(size, guide, i  , j  , k  );
          s = bh1[m];
          
          if ( GET_FACE_BC(s, BC_FACE_W) == n ) { // vertical
            q = ht1 * (sf - 0.0); // 基準温度との温度差
            qbc[FBUtility::getFindexV3DEx(size, guide, 0, i-1, j  , k  )] += q;
            va += q; // マイナス面の正の値は流入
            t[FBUtility::getFindexS3D(size, guide, i-1, j  , k  )] = sf; // 固体セルに表面温度を代入
            flop += 4.0;
          }
          
          if ( GET_FACE_BC(s, BC_FACE_E) == n ) { // vertical
            q = ht1 * (0.0 - sf);
            qbc[FBUtility::getFindexV3DEx(size, guide, 0, i  , j  , k  )] += q;
            va -= q; // プラス面の正の値は流出
            t[FBUtility::getFindexS3D(size, guide, i+1, j  , k  )] = sf;
            flop += 4.0;
          }
          
          if ( GET_FACE_BC(s, BC_FACE_S) == n ) { // vertical
            q = ht1 * (sf - 0.0);
            qbc[FBUtility::getFindexV3DEx(size, guide, 1, i  , j-1, k  )] += q;
            va += q;
            t[FBUtility::getFindexS3D(size, guide, i  , j-1, k  )] = sf;
            flop += 4.0;
          }
          
          if ( GET_FACE_BC(s, BC_FACE_N) == n ) { // vertical
            q = ht1 * (0.0 - sf);
            qbc[FBUtility::getFindexV3DEx(size, guide, 1, i  , j  , k  )] += q;
            va -= q;
            t[FBUtility::getFindexS3D(size, guide, i  , j+1, k  )] = sf;
            flop += 4.0;
          }
          
          if ( GET_FACE_BC(s, BC_FACE_B) == n ) { // horizontal
            q = ht3 * (sf - 0.0);
            qbc[FBUtility::getFindexV3DEx(size, guide, 2, i  , j  , k-1)] += q;
            va += q;
            t[FBUtility::getFindexS3D(size, guide, i  , j  , k-1)] = sf;
            flop += 4.0;
          }
          
          if ( GET_FACE_BC(s, BC_FACE_T) == n ) { // horizontal
            q = ht1 * (0.0 - sf);
            qbc[FBUtility::getFindexV3DEx(size, guide, 2, i  , j  , k  )] += q;
            va -= q;
            t[FBUtility::getFindexS3D(size, guide, i  , j  , k+1)] = sf;
            flop += 4.0;
          }
        }
      }
    }
  }
  else { // Local
    for (k=st[2]; k<=ed[2]; k++) {
      for (j=st[1]; j<=ed[1]; j++) {
        for (i=st[0]; i<=ed[0]; i++) {
          
          m = FBUtility::getFindexS3D(size, guide, i  , j  , k  );
          s = bh1[m];
          
          if ( GET_FACE_BC(s, BC_FACE_W) == n ) {
            q = ht1 * (sf - t0[m]); // 基準温度との温度差
            qbc[FBUtility::getFindexV3DEx(size, guide, 0, i-1, j  , k  )] += q;
            va += q; // マイナス面の正の値は流入
            t[FBUtility::getFindexS3D(size, guide, i-1, j  , k  )] = sf; // 固体セルに表面温度を代入
            flop += 4.0;
          }
          
          if ( GET_FACE_BC(s, BC_FACE_E) == n ) {
            q = ht1 * (t0[m] - sf);
            qbc[FBUtility::getFindexV3DEx(size, guide, 0, i  , j  , k  )] += q;
            va -= q; // プラス面の正の値は流出
            t[FBUtility::getFindexS3D(size, guide, i+1, j  , k  )] = sf;
            flop += 4.0;
          }
          
          if ( GET_FACE_BC(s, BC_FACE_S) == n ) {
            q = ht1 * (sf - t0[m]);
            qbc[FBUtility::getFindexV3DEx(size, guide, 1, i  , j-1, k  )] += q;
            va += q;
            t[FBUtility::getFindexS3D(size, guide, i  , j-1, k  )] = sf;
            flop += 4.0;
          }
          
          if ( GET_FACE_BC(s, BC_FACE_N) == n ) {
            q = ht1 * (t0[m] - sf);
            qbc[FBUtility::getFindexV3DEx(size, guide, 1, i  , j  , k  )] += q;
            va -= q;
            t[FBUtility::getFindexS3D(size, guide, i  , j+1, k  )] = sf;
            flop += 4.0;
          }
          
          if ( GET_FACE_BC(s, BC_FACE_B) == n ) {
            q = ht3 * (sf - t0[m]);
            qbc[FBUtility::getFindexV3DEx(size, guide, 2, i  , j  , k-1)] += q;
            va += q;
            t[FBUtility::getFindexS3D(size, guide, i  , j  , k-1)] = sf;
            flop += 4.0;
          }
          
          if ( GET_FACE_BC(s, BC_FACE_T) == n ) {
            q = ht1 * (t0[m] - sf);
            qbc[FBUtility::getFindexV3DEx(size, guide, 2, i  , j  , k  )] += q;
            va -= q;
            t[FBUtility::getFindexS3D(size, guide, i  , j  , k+1)] = sf;
            flop += 4.0;
          }
        }
      }
    }
  }
  
  if( para_mng->IsParallel() ){
    REAL_TYPE tmp = va;
    para_mng->Allreduce(&tmp, &va, 1, SKL_ARRAY_DTYPE_REAL, SKL_SUM, pn.procGrp);
  }
  
  return ( va ); // 単位面積あたりの無次元熱流束の和
}

/**
 @fn REAL_TYPE SetBC3D::ps_IBC_Transfer_SF_SM(REAL_TYPE* qbc, unsigned* bh1, int n, REAL_TYPE* t, REAL_TYPE* t0, REAL_TYPE& flop)
 @brief 単媒質の場合の熱伝達温境界条件タイプSF（強制対流）　流体のみを解く場合
 @retval 熱流束の和 (W/m^2)
 @param qbc 境界条件熱流束
 @param bh1 BCindex H1
 @param n 境界条件コンポーネントのエントリ番号
 @param t n+1時刻の温度場
 @param t0 n時刻の温度場
 @param flop 浮動小数演算数
 @note 
    - 熱流束は加算（他の条件と合成）
    - pセルは流体セル
 */
REAL_TYPE SetBC3D::ps_IBC_Transfer_SF_SM(REAL_TYPE* qbc, unsigned* bh1, int n, REAL_TYPE* t, REAL_TYPE* t0, REAL_TYPE& flop)
{
  SklParaManager* para_mng = ParaCmpo->GetParaManager();
  int i,j,k;
  unsigned s, m;
  REAL_TYPE ht, sf, q, va=0.0;
  REAL_TYPE a1, b1, c1;
  int st[3], ed[3];
  
  a1 = cmp[n].ca[CompoList::alpha];
  b1 = cmp[n].ca[CompoList::beta];
  c1 = cmp[n].ca[CompoList::gamma];
  ht = lambda / (RefV*RefL*rho*cp) * a1*pow(Reynolds, b1) * pow(Prandtl, c1);    // Ref_IDで指定される物性値 
  sf = FBUtility::convK2ND(cmp[n].get_Temp(), BaseTemp, DiffTemp); // 保持されている温度はKelvin
  
  cmp[n].getBbox(st, ed);
  
  if ( cmp[n].get_sw_HTmodeRef() == CompoList::HT_mode_bulk ) {
    for (k=st[2]; k<=ed[2]; k++) {
      for (j=st[1]; j<=ed[1]; j++) {
        for (i=st[0]; i<=ed[0]; i++) {
          
          m = FBUtility::getFindexS3D(size, guide, i  , j  , k  );
          s = bh1[m];
          
          if ( GET_FACE_BC(s, BC_FACE_W) == n ) {
            q = ht * (sf - 0.0); // 基準温度との温度差
            qbc[FBUtility::getFindexV3DEx(size, guide, 0, i-1, j  , k  )] += q;
            va += q; // マイナス面の正の値は流入
            t[FBUtility::getFindexS3D(size, guide, i-1, j  , k  )] = sf; // 固体セルに表面温度を代入
            flop += 4.0;
          }
          
          if ( GET_FACE_BC(s, BC_FACE_E) == n ) {
            q = ht * (0.0 - sf);
            qbc[FBUtility::getFindexV3DEx(size, guide, 0, i  , j  , k  )] += q;
            va -= q; // プラス面の正の値は流出
            t[FBUtility::getFindexS3D(size, guide, i+1, j  , k  )] = sf;
            flop += 4.0;
          }
          
          if ( GET_FACE_BC(s, BC_FACE_S) == n ) {
            q = ht * (sf - 0.0);
            qbc[FBUtility::getFindexV3DEx(size, guide, 1, i  , j-1, k  )] += q;
            va += q;
            t[FBUtility::getFindexS3D(size, guide, i  , j-1, k  )] = sf;
            flop += 4.0;
          }
          
          if ( GET_FACE_BC(s, BC_FACE_N) == n ) {
            q = ht * (0.0 - sf);
            qbc[FBUtility::getFindexV3DEx(size, guide, 1, i  , j  , k  )] += q;
            va -= q;
            t[FBUtility::getFindexS3D(size, guide, i  , j+1, k  )] = sf;
            flop += 4.0;
          }
          
          if ( GET_FACE_BC(s, BC_FACE_B) == n ) {
            q = ht * (sf - 0.0);
            qbc[FBUtility::getFindexV3DEx(size, guide, 2, i  , j  , k-1)] += q;
            va += q;
            t[FBUtility::getFindexS3D(size, guide, i  , j  , k-1)] = sf;
            flop += 4.0;
          }
          
          if ( GET_FACE_BC(s, BC_FACE_T) == n ) {
            q = ht * (0.0 - sf);
            qbc[FBUtility::getFindexV3DEx(size, guide, 2, i  , j  , k  )] += q;
            va -= q;
            t[FBUtility::getFindexS3D(size, guide, i  , j  , k+1)] = sf;
            flop += 4.0;
          }
        }
      }
    }
  }
  else { // Local
    for (k=st[2]; k<=ed[2]; k++) {
      for (j=st[1]; j<=ed[1]; j++) {
        for (i=st[0]; i<=ed[0]; i++) {
          
          m = FBUtility::getFindexS3D(size, guide, i  , j  , k  );
          s = bh1[m];
          
          if ( GET_FACE_BC(s, BC_FACE_W) == n ) {
            q = ht * (sf - t0[m]); // 基準温度との温度差
            qbc[FBUtility::getFindexV3DEx(size, guide, 0, i-1, j  , k  )] += q;
            va += q; // マイナス面の正の値は流入
            t[FBUtility::getFindexS3D(size, guide, i-1, j  , k  )] = sf; // 固体セルに表面温度を代入
            flop += 4.0;
          }
          
          if ( GET_FACE_BC(s, BC_FACE_E) == n ) {
            q = ht * (t0[m] - sf);
            qbc[FBUtility::getFindexV3DEx(size, guide, 0, i  , j  , k  )] += q;
            va -= q; // プラス面の正の値は流出
            t[FBUtility::getFindexS3D(size, guide, i+1, j  , k  )] = sf;
            flop += 4.0;
          }
          
          if ( GET_FACE_BC(s, BC_FACE_S) == n ) {
            q = ht * (sf - t0[m]);
            qbc[FBUtility::getFindexV3DEx(size, guide, 1, i  , j-1, k  )] += q;
            va += q;
            t[FBUtility::getFindexS3D(size, guide, i  , j-1, k  )] = sf;
            flop += 4.0;
          }
          
          if ( GET_FACE_BC(s, BC_FACE_N) == n ) {
            q = ht * (t0[m] - sf);
            qbc[FBUtility::getFindexV3DEx(size, guide, 1, i  , j  , k  )] += q;
            va -= q;
            t[FBUtility::getFindexS3D(size, guide, i  , j+1, k  )] = sf;
            flop += 4.0;
          }
          
          if ( GET_FACE_BC(s, BC_FACE_B) == n ) {
            q = ht * (sf - t0[m]);
            qbc[FBUtility::getFindexV3DEx(size, guide, 2, i  , j  , k-1)] += q;
            va += q;
            t[FBUtility::getFindexS3D(size, guide, i  , j  , k-1)] = sf;
            flop += 4.0;
          }
          
          if ( GET_FACE_BC(s, BC_FACE_T) == n ) {
            q = ht * (t0[m] - sf);
            qbc[FBUtility::getFindexV3DEx(size, guide, 2, i  , j  , k  )] += q;
            va -= q;
            t[FBUtility::getFindexS3D(size, guide, i  , j  , k+1)] = sf;
            flop += 4.0;
          }
        }
      }
    }
  }

  if( para_mng->IsParallel() ){
    REAL_TYPE tmp = va;
    para_mng->Allreduce(&tmp, &va, 1, SKL_ARRAY_DTYPE_REAL, SKL_SUM, pn.procGrp);
  }

  return ( va ); // 単位面積あたりの無次元熱流束の和
}

/**
 @fn REAL_TYPE SetBC3D::ps_IBC_Transfer_B_SM(REAL_TYPE* qbc, unsigned* bh1, int n, REAL_TYPE* t, REAL_TYPE* t0, REAL_TYPE& flop)
 @brief 単媒質の場合の熱伝達境界条件タイプB 固体のみを解く場合
 @retval 熱流束の和 (W/m^2)
 @param qbc 境界条件熱流束
 @param bh1 BCindex H1
 @param n 境界条件コンポーネントのエントリ番号
 @param t n+1時刻の温度場
 @param t0 n時刻の温度場
 @param flop 浮動小数演算数
 @note 
    - 熱流束は加算（他の条件と合成）
    - pセルは固体セル
 */
REAL_TYPE SetBC3D::ps_IBC_Transfer_B_SM(REAL_TYPE* qbc, unsigned* bh1, int n, REAL_TYPE* t, REAL_TYPE* t0, REAL_TYPE& flop)
{
  SklParaManager* para_mng = ParaCmpo->GetParaManager();
  int i,j,k;
  unsigned s, odr, m;
  REAL_TYPE ht, bt, q, va=0.0;
  int st[3], ed[3];

  //odr= cmp[n].getMatOdr();
  //ht = cmp[n].get_CoefHT() / (RefV*mat[odr].P[p_density]*mat[odr].P[p_specific_heat]);  // 固体セルで指定される物性値
  ht = cmp[n].get_CoefHT() / (RefV*rho*cp);  // Ref_IDで指定される物性値
  bt = FBUtility::convK2ND(cmp[n].get_Temp(), BaseTemp, DiffTemp); // 保持されている温度はKelvin

  cmp[n].getBbox(st, ed);
  
  for (k=st[2]; k<=ed[2]; k++) {
    for (j=st[1]; j<=ed[1]; j++) {
      for (i=st[0]; i<=ed[0]; i++) {
        
        m = FBUtility::getFindexS3D(size, guide, i  , j  , k  );
        s = bh1[m];
        
        if ( GET_FACE_BC(s, BC_FACE_W) == n ) {
          q = ht * (bt - t0[m]);
          qbc[FBUtility::getFindexV3DEx(size, guide, 0, i-1, j  , k  )] += q;
          va += q; // マイナス面の正の値は流入
          t[FBUtility::getFindexS3D(size, guide, i-1, j  , k  )] = bt; // 隣接流体セルにバルク温度を代入
          flop += 4.0;
        }
        
        if ( GET_FACE_BC(s, BC_FACE_E) == n ) {
          q = ht * (t0[m] - bt);
          qbc[FBUtility::getFindexV3DEx(size, guide, 0, i  , j  , k  )] += q;
          va -= q; // プラス面の正の値は流出
          t[FBUtility::getFindexS3D(size, guide, i+1, j  , k  )] = bt;
          flop += 4.0;
        }
        
        if ( GET_FACE_BC(s, BC_FACE_S) == n ) {
          q = ht * (bt - t0[m]);
          qbc[FBUtility::getFindexV3DEx(size, guide, 1, i  , j-1, k  )] += q;
          va += q;
          t[FBUtility::getFindexS3D(size, guide, i  , j-1, k  )] = bt;
          flop += 4.0;
        }
        
        if ( GET_FACE_BC(s, BC_FACE_N) == n ) {
          q = ht * (t0[m] - bt);
          qbc[FBUtility::getFindexV3DEx(size, guide, 1, i  , j  , k  )] += q;
          va -= q;
          t[FBUtility::getFindexS3D(size, guide, i  , j+1, k  )] = bt;
          flop += 4.0;
        }
        
        if ( GET_FACE_BC(s, BC_FACE_B) == n ) {
          q = ht * (bt - t0[m]);
          qbc[FBUtility::getFindexV3DEx(size, guide, 2, i  , j  , k-1)] += q;
          va += q;
          t[FBUtility::getFindexS3D(size, guide, i  , j  , k-1)] = bt;
          flop += 4.0;
        }
        
        if ( GET_FACE_BC(s, BC_FACE_T) == n ) {
          q = ht * (t0[m] - bt);
          qbc[FBUtility::getFindexV3DEx(size, guide, 2, i  , j  , k  )] += q;
          va -= q;
          t[FBUtility::getFindexS3D(size, guide, i  , j  , k+1)] = bt;
          flop += 4.0;
        }
      }
    }
  }
  
  if( para_mng->IsParallel() ){
		REAL_TYPE tmp = va;
		para_mng->Allreduce(&tmp, &va, 1, SKL_ARRAY_DTYPE_REAL, SKL_SUM, pn.procGrp);
	}
  
  return ( va ); // 単位面積あたりの無次元熱流束の和
}

/**
 @fn REAL_TYPE SetBC3D::setHeatTransferB(REAL_TYPE* qbc, REAL_TYPE* t, unsigned* bx, unsigned n, REAL_TYPE* t0)
 @brief 単媒質の場合の熱伝達境界条件タイプB
 @retval sum of heat flux (W/m^2)
 @param qbc 境界条件熱流束
 @param t n+1時刻の温度場
 @param bx BCindex
 @param n 境界条件コンポーネントのエントリ番号
 @param t0 n時刻の温度場
 
REAL_TYPE SetBC3D::setHeatTransferB(REAL_TYPE* qbc, REAL_TYPE* t, unsigned* bx, unsigned n, REAL_TYPE* t0)
{
  SklParaManager* para_mng = ParaCmpo->GetParaManager();
  int i,j,k;
  unsigned s, m0;
  int wi, wj, wk;
  unsigned li, lj, lk;
  unsigned ni, nj, nk;
  REAL_TYPE qi, qj, qk, pp, bt, qsum=0.0, tmp, c;
  REAL_TYPE di, dj, dk; // outer normal
  
  pp = cmp[n].get_CoefHT() / RefV;
  bt = FBUtility::convK2ND(cmp[n].get_Temp(), BaseTemp, DiffTemp); // 保持されている温度はKelvin
  
  for (k=(int)cmp[n].ci.st[2]-1; k<=(int)cmp[n].ci.ed[2]; k++) {
    for (j=(int)cmp[n].ci.st[1]-1; j<=(int)cmp[n].ci.ed[1]; j++) {
      for (i=(int)cmp[n].ci.st[0]-1; i<=(int)cmp[n].ci.ed[0]; i++) {
        m0 = FBUtility::getFindexS3D(size, guide, i  , j  , k  );
        s = bx[m0];
        if ( (s & BIT_MASK_10) == n ) {
          
          if ( BIT_IS_SHIFT(s, QFACE_I) ) {
            wi = BIT_SHIFT(s, DIR_I);
            t[ FBUtility::getFindexS3D(size, guide, i+wi, j, k) ] = bt; // 流体セルにバルク温度を代入
            li = FBUtility::getFindexS3D(size, guide, i+1-wi, j, k); // 固体セルのセルインデクス
            ni = DECODE_MAT( bx[li] ); // 固体セルのマテリアルのオーダーを得る
            di = 2.0*(REAL_TYPE)wi-1.0;
            qi = pp / (mat[ni].P[p_density] * mat[ni].P[p_specific_heat]) * (t0[li]-bt) * di;
            qbc[FBUtility::getFindexV3DEx(size, guide, 0, i, j, k)] = qi;
            qsum += qi*di;
          }
          if ( BIT_IS_SHIFT(s, QFACE_J) ) {
            wj = BIT_SHIFT(s, DIR_J);
            t[ FBUtility::getFindexS3D(size, guide, i, j+wj, k) ] = bt;
            lj = FBUtility::getFindexS3D(size, guide, i, j+1-wj, k);
            nj = DECODE_MAT( bx[lj] );
            dj = 2.0*(REAL_TYPE)wj-1.0;
            qj = pp / (mat[nj].P[p_density] * mat[nj].P[p_specific_heat]) * (t0[lj]-bt) * dj;
            qbc[FBUtility::getFindexV3DEx(size, guide, 1, i, j, k)] = qj;
            qsum += qj*dj;
          }
          if ( BIT_IS_SHIFT(s, QFACE_K) ) {
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
  
  if( para_mng->IsParallel() ){
		tmp = qsum;
		para_mng->Allreduce(&tmp, &qsum, 1, SKL_ARRAY_DTYPE_REAL, SKL_SUM, pn.procGrp);
	}
  c = qsum * (RefV*DiffTemp*rho*cp);
  
  return (c);
}*/

/**
 @fn REAL_TYPE SetBC3D::ps_IBC_IsoThermal_SM(REAL_TYPE* qbc, unsigned* bh1, int n, REAL_TYPE* t, REAL_TYPE* t0, REAL_TYPE& flop)
 @brief 単媒質の場合の等温境界条件
 @retval 熱流束の和 (W/m^2)
 @param qbc 境界条件熱流束
 @param bh1 BCindex H1
 @param n 境界条件コンポーネントのエントリ番号
 @param t n+1時刻の温度場
 @param t0 n時刻の温度場
 @note 
    - 熱流束は加算（他の条件と合成）
    - pセルは等温壁をもつ計算セルで，等温のセルは相手方
 */
REAL_TYPE SetBC3D::ps_IBC_IsoThermal_SM(REAL_TYPE* qbc, unsigned* bh1, int n, REAL_TYPE* t, REAL_TYPE* t0, REAL_TYPE& flop)
{
  SklParaManager* para_mng = ParaCmpo->GetParaManager();
  int i,j,k;
  unsigned s, m, odr;
  REAL_TYPE q, pp, sf, va=0.0;
  int st[3], ed[3];

  odr= cmp[n].getMatOdr();
  sf = FBUtility::convK2ND(cmp[n].get_Temp(), BaseTemp, DiffTemp); // 保持されている温度はKelvin
  //pp = (2.0*mat[odr].P[p_thermal_conductivity]) / (dh*RefV*RefL*mat[odr].P[p_density]*mat[odr].P[p_specific_heat]); // property of solid cell
  pp = (2.0*lambda) / (dh*RefV*rho*cp); // Ref_IDで指定される物性値

  cmp[n].getBbox(st, ed);
  
  for (k=st[2]; k<=ed[2]; k++) {
    for (j=st[1]; j<=ed[1]; j++) {
      for (i=st[0]; i<=ed[0]; i++) {
        
        m = FBUtility::getFindexS3D(size, guide, i  , j  , k  );
        s = bh1[m];
        
        if ( GET_FACE_BC(s, BC_FACE_W) == n ) {
          q = pp * (sf - t0[m]);
          qbc[FBUtility::getFindexV3DEx(size, guide, 0, i-1, j  , k  )] += q;
          va += q; // マイナス面の正の値は流入
          t[FBUtility::getFindexS3D(size, guide, i-1, j  , k  )] = sf; // 等温セルに指定温度を代入
          flop += 4.0;
        }
        
        if ( GET_FACE_BC(s, BC_FACE_E) == n ) {
          q = pp * (t0[m] - sf);
          qbc[FBUtility::getFindexV3DEx(size, guide, 0, i  , j  , k  )] += q;
          va -= q; // プラス面の正の値は流出
          t[FBUtility::getFindexS3D(size, guide, i+1, j  , k  )] = sf;
          flop += 4.0;
        }
        
        if ( GET_FACE_BC(s, BC_FACE_S) == n ) {
          q = pp * (sf - t0[m]);
          qbc[FBUtility::getFindexV3DEx(size, guide, 1, i  , j-1, k  )] += q;
          va += q;
          t[FBUtility::getFindexS3D(size, guide, i  , j-1, k  )] = sf;
          flop += 4.0;
        }
        
        if ( GET_FACE_BC(s, BC_FACE_N) == n ) {
          q = pp * (t0[m] - sf);
          qbc[FBUtility::getFindexV3DEx(size, guide, 1, i  , j  , k  )] += q;
          va -= q;
          t[FBUtility::getFindexS3D(size, guide, i  , j+1, k  )] = sf;
          flop += 4.0;
        }
        
        if ( GET_FACE_BC(s, BC_FACE_B) == n ) {
          q = pp * (sf - t0[m]);
          qbc[FBUtility::getFindexV3DEx(size, guide, 2, i  , j  , k-1)] += q;
          va += q;
          t[FBUtility::getFindexS3D(size, guide, i  , j  , k-1)] = sf;
          flop += 4.0;
        }
        
        if ( GET_FACE_BC(s, BC_FACE_T) == n ) {
          q = pp * (t0[m] - sf);
          qbc[FBUtility::getFindexV3DEx(size, guide, 2, i  , j  , k  )] += q;
          va -= q;
          t[FBUtility::getFindexS3D(size, guide, i  , j  , k+1)] = sf;
          flop += 4.0;
        }
      }
    }
  }
  
  if( para_mng->IsParallel() ){
		REAL_TYPE tmp = va;
		para_mng->Allreduce(&tmp, &va, 1, SKL_ARRAY_DTYPE_REAL, SKL_SUM, pn.procGrp);
	}
  
  return ( va ); // 単位面積あたりの無次元熱流束の和
}

/**
 @fn REAL_TYPE SetBC3D::setIsoThermal(REAL_TYPE* qbc, REAL_TYPE* t, unsigned* bx, unsigned n, REAL_TYPE* t0)
 @brief 単媒質の場合の等温境界条件
 @retval sum of heat flux (W/m^2)
 @param qbc 境界条件熱流束
 @param t n+1時刻の温度場
 @param bx BCindex
 @param n 境界条件コンポーネントのエントリ番号
 @param t0 n時刻の温度場
 
REAL_TYPE SetBC3D::setIsoThermal(REAL_TYPE* qbc, REAL_TYPE* t, unsigned* bx, unsigned n, REAL_TYPE* t0)
{
  SklParaManager* para_mng = ParaCmpo->GetParaManager();
  int i,j,k;
  unsigned s, m0;
  int wi, wj, wk;
  unsigned li, lj, lk;
  unsigned ni, nj, nk;
  REAL_TYPE qi, qj, qk, pp, st, qsum=0.0, tmp, c;
  REAL_TYPE di, dj, dk; // outer normal

  pp = 2.0/(dh*RefV*RefL);
  st = FBUtility::convK2ND(cmp[n].get_Temp(), BaseTemp, DiffTemp); // 保持されている温度はKelvin
  
  for (k=(int)cmp[n].ci.st[2]-1; k<=(int)cmp[n].ci.ed[2]; k++) {
    for (j=(int)cmp[n].ci.st[1]-1; j<=(int)cmp[n].ci.ed[1]; j++) {
      for (i=(int)cmp[n].ci.st[0]-1; i<=(int)cmp[n].ci.ed[0]; i++) {
        m0 = FBUtility::getFindexS3D(size, guide, i  , j  , k  );
        s = bx[m0];
        if ( (s & BIT_MASK_10) == n ) {

          if ( BIT_IS_SHIFT(s, QFACE_I) ) {
            wi = BIT_SHIFT(s, DIR_I);
            t[ FBUtility::getFindexS3D(size, guide, i+wi, j, k) ] = st; // 等温セルに指定温度を代入
            li = FBUtility::getFindexS3D(size, guide, i+1-wi, j, k); // 計算側のセルインデクス
            ni = DECODE_MAT( bx[li] ); // 計算セルのマテリアルのオーダーを得る
            di = 1.0-2.0*(REAL_TYPE)wi;
            qi = -pp * mat[ni].P[p_thermal_conductivity] / (mat[ni].P[p_density] * mat[ni].P[p_specific_heat]) * (t0[li]-st) * di;
            qbc[FBUtility::getFindexV3DEx(size, guide, 0, i, j, k)] = qi;
            qsum += qi*di;
          }
          if ( BIT_IS_SHIFT(s, QFACE_J) ) {
            wj = BIT_SHIFT(s, DIR_J);
            t[ FBUtility::getFindexS3D(size, guide, i, j+wj, k) ] = st;
            lj = FBUtility::getFindexS3D(size, guide, i, j+1-wj, k);
            nj = DECODE_MAT( bx[lj] );
            dj = 1.0-2.0*(REAL_TYPE)wj;
            qj = -pp * mat[nj].P[p_thermal_conductivity] / (mat[nj].P[p_density] * mat[nj].P[p_specific_heat]) * (t0[lj]-st) * dj;
            qbc[FBUtility::getFindexV3DEx(size, guide, 1, i, j, k)] = qj;
            qsum += qj*dj;
          }
          if ( BIT_IS_SHIFT(s, QFACE_K) ) {
            wk = BIT_SHIFT(s, DIR_K);
            t[ FBUtility::getFindexS3D(size, guide, i, j, k+wk) ] = st;
            lk = FBUtility::getFindexS3D(size, guide, i, j, k+1-wk);
            nk = DECODE_MAT( bx[lk] );
            dk = 1.0-2.0*(REAL_TYPE)wk;
            qk = -pp * mat[nk].P[p_thermal_conductivity] / (mat[nk].P[p_density] * mat[nk].P[p_specific_heat]) * (t0[lk]-st) * dk;
            qbc[FBUtility::getFindexV3DEx(size, guide, 2, i, j, k)] = qk;
            qsum += qk*dk;
          }
        }
      }
    }
  }
  
  if( para_mng->IsParallel() ){
		tmp = qsum;
		para_mng->Allreduce(&tmp, &qsum, 1, SKL_ARRAY_DTYPE_REAL, SKL_SUM, pn.procGrp);
	}
  c = qsum * (RefV*DiffTemp*rho*cp);
  
  return (c);
}*/

/**
 @fn void SetBC3D::setInitialTemp_Compo(unsigned n, unsigned* bx, SklScalar3D<REAL_TYPE>* d_t)
 @brief 初期温度を代入
 @param n エントリ
 @param bx BCindex ID
 @param d_t 
 */
void SetBC3D::setInitialTemp_Compo(unsigned n, unsigned* bx, SklScalar3D<REAL_TYPE>* d_t)
{
  int i,j,k;
  unsigned register m, id;
  REAL_TYPE ref;
  REAL_TYPE *t=NULL;
  
  if ( !(t = d_t->GetData()) ) Exit(0);
  
  ref = FBUtility::convK2ND(cmp[n].getInitTemp(), BaseTemp, DiffTemp);
  id = cmp[n].getID();
	
  for (k=1; k<=kxc; k++) {
    for (j=1; j<=jxc; j++) {
      for (i=1; i<=ixc; i++) {
        m = FBUtility::getFindexS3D(size, guide, i  , j  , k  );
        if ( DECODE_ID(bx[m]) == id ) {
          t[m] = ref; 
        }
      }
    }
  }
}

/**
 @fn void SetBC3D::setBCIperiodic(SklScalar3D<unsigned>* d_bx)
 @brief 周期境界の場合のインデクスの同期
 @param bx BCindexのデータクラス
 */
void SetBC3D::setBCIperiodic(SklScalar3D<unsigned>* d_bx)
{
  SklParaManager* para_mng = ParaCmpo->GetParaManager();
  int i,j,k, gd;
  unsigned m0, m1;
  unsigned* bx=NULL;
  
  gd = (int)guide;
  
  if ( !(bx = d_bx->GetData()) ) Exit(0);

  if( para_mng->IsParallel() ){
    for (int face=0; face<NOFACE; face++) {
      if ( obc[face].get_BCtype() != OBC_PERIODIC ) continue; // スキップしてfaceをインクリメント

      switch (face) {
        case X_MINUS:
          if ( !d_bx->CommPeriodicBndCell(PRDC_X_DIR, PRDC_PLUS2MINUS, guide) ) Exit(0);
          break;
          
        case X_PLUS:
          if ( !d_bx->CommPeriodicBndCell(PRDC_X_DIR, PRDC_MINUS2PLUS, guide) ) Exit(0);
          break;
          
        case Y_MINUS:
          if ( !d_bx->CommPeriodicBndCell(PRDC_Y_DIR, PRDC_PLUS2MINUS, guide) ) Exit(0);
          break;
          
        case Y_PLUS:
          if ( !d_bx->CommPeriodicBndCell(PRDC_Y_DIR, PRDC_MINUS2PLUS, guide) ) Exit(0);
          break;
          
        case Z_MINUS:
          if ( !d_bx->CommPeriodicBndCell(PRDC_Z_DIR, PRDC_PLUS2MINUS, guide) ) Exit(0);
          break;
          
        case Z_PLUS:
          if ( !d_bx->CommPeriodicBndCell(PRDC_Z_DIR, PRDC_MINUS2PLUS, guide) ) Exit(0);
          break;
      }
    }
  } 
  else { // 逐次処理
    for (int face=0; face<NOFACE; face++) {
      if ( obc[face].get_BCtype() != OBC_PERIODIC ) continue; // スキップしてfaceをインクリメント
      
      switch (face) {
        case X_MINUS:
          if( pn.nID[face] < 0 ) {
            for (k=1; k<=kxc; k++) {
              for (j=1; j<=jxc; j++) {
                for (i=1-gd; i<=0; i++) {
                  m0 = FBUtility::getFindexS3D(size, guide, i     , j, k);
                  m1 = FBUtility::getFindexS3D(size, guide, i+ixc, j, k);
                  bx[m0] = bx[m1];
                }
              }
            }
          }
          break;
          
        case X_PLUS:
          if( pn.nID[face] < 0 ) {
            for (k=1; k<=kxc; k++) {
              for (j=1; j<=jxc; j++) {
                for (i=ixc+1; i<=ixc+gd; i++) {
                  m0 = FBUtility::getFindexS3D(size, guide, i     , j, k);
                  m1 = FBUtility::getFindexS3D(size, guide, i-ixc, j, k);
                  bx[m0] = bx[m1];
                }
              }
            }
          }
          break;
          
        case Y_MINUS:
          if( pn.nID[face] < 0 ) {
            for (k=1; k<=kxc; k++) {
              for (j=1-gd; j<=0; j++) {
                for (i=1; i<=ixc; i++) {
                  m0 = FBUtility::getFindexS3D(size, guide, i, j     , k);
                  m1 = FBUtility::getFindexS3D(size, guide, i, j+jxc, k);
                  bx[m0] = bx[m1];
                }
              }
            }
          }
          break;
          
        case Y_PLUS:
          if( pn.nID[face] < 0 ) {
            for (k=1; k<=kxc; k++) {
              for (j=jxc+1; j<=jxc+gd; j++) {
                for (i=1; i<=ixc; i++) {
                  m0 = FBUtility::getFindexS3D(size, guide, i, j     , k);
                  m1 = FBUtility::getFindexS3D(size, guide, i, j-jxc, k);
                  bx[m0] = bx[m1];
                }
              }
            }
          }
          break;
          
        case Z_MINUS:
          if( pn.nID[face] < 0 ) {
            for (k=1-gd; k<=0; k++) {
              for (j=1; j<=jxc; j++) {
                for (i=1; i<=ixc; i++) {
                  m0 = FBUtility::getFindexS3D(size, guide, i, j, k     );
                  m1 = FBUtility::getFindexS3D(size, guide, i, j, k+kxc);
                  bx[m0] = bx[m1];
                }
              }
            }
          }
          break;
          
        case Z_PLUS:
          if( pn.nID[face] < 0 ) {
            for (k=kxc+1; k<=kxc+gd; k++) {
              for (j=1; j<=jxc; j++) {
                for (i=1; i<=ixc; i++) {
                  m0 = FBUtility::getFindexS3D(size, guide, i, j, k     );
                  m1 = FBUtility::getFindexS3D(size, guide, i, j, k-kxc);
                  bx[m0] = bx[m1];
                }
              }
            }
          }
          break;
      }    
    }
  }
}

/**
 @fn SetBC3D::setRadiant(REAL_TYPE* qbc, unsigned* bx, unsigned n, REAL_TYPE* t)
 @brief 輻射境界条件
 @param qbc 境界条件熱流束
 @param bx BCindex
 @param n 境界条件コンポーネントのエントリ番号
 @param t n時刻の温度場
 @note
    - 未使用，形式のみ
 
void SetBC3D::setRadiant(REAL_TYPE* qbc, unsigned* bx, unsigned n, REAL_TYPE* t)
{
  int i,j,k;
  unsigned s;
  unsigned m0, mi, mj, mk;
  REAL_TYPE ht;

  ht = cmp[n].get_CoefHT() / (mat[n].P[p_density]*mat[n].P[p_specific_heat]);

  for (k=(int)cmp[n].ci.st[2]-1; k<=(int)cmp[n].ci.ed[2]; k++) {
    for (j=(int)cmp[n].ci.st[1]-1; j<=(int)cmp[n].ci.ed[1]; j++) {
      for (i=(int)cmp[n].ci.st[0]-1; i<=(int)cmp[n].ci.ed[0]; i++) {
        m0 = FBUtility::getFindexS3D(size, guide, i  , j  , k  );
        mi = FBUtility::getFindexS3D(size, guide, i+1, j  , k  );
        mj = FBUtility::getFindexS3D(size, guide, i  , j+1, k  );
        mk = FBUtility::getFindexS3D(size, guide, i  , j  , k+1);
        s = bx[m0];
        if ( (s & BIT_MASK_10) == n ) {
          if ( BIT_IS_SHIFT(s, QFACE_I) ) qbc[FBUtility::getFindexV3D(size, guide, 0, i, j, k)] = ht*(t[m0]-t[mi]);
          if ( BIT_IS_SHIFT(s, QFACE_J) ) qbc[FBUtility::getFindexV3D(size, guide, 1, i, j, k)] = ht*(t[m0]-t[mj]);
          if ( BIT_IS_SHIFT(s, QFACE_K) ) qbc[FBUtility::getFindexV3D(size, guide, 2, i, j, k)] = ht*(t[m0]-t[mk]);
        }
      }
    }
  }
}*/

/**
 @fn void SetBC3D::mod_Vis_CN(REAL_TYPE* vc, REAL_TYPE* wv, REAL_TYPE cf, unsigned* bx, REAL_TYPE* wk, REAL_TYPE tm, REAL_TYPE dt, REAL_TYPE omg, Control* C, REAL_TYPE* res, REAL_TYPE* v00, REAL_TYPE& flop)
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
 
void SetBC3D::mod_Vis_CN(REAL_TYPE* vc, REAL_TYPE* wv, REAL_TYPE cf, unsigned* bx, REAL_TYPE* wk, REAL_TYPE tm, REAL_TYPE dt, REAL_TYPE omg, Control* C, REAL_TYPE* res, unsigned LS, REAL_TYPE* v00, REAL_TYPE& flop)
{
  int st[3], ed[3];
  REAL_TYPE vel, vec[3], a, b;
  unsigned typ;
  const REAL_TYPE c_pai = (REAL_TYPE)(2.0*asin(1.0));
  
  for (unsigned n=1; n<=NoBC; n++) {
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
          cbc_vis_cn_mod_jcb_(vc, dim_sz, gc, st, ed, &dh, &dt, v00, &rei, &omg, wv, (int*)bx, wk, &cf, res, vec, flop);
          break;
          
        case SOR:
          cbc_vis_cn_mod_sor_(vc, dim_sz, gc, st, ed, &dh, &dt, v00, &rei, &omg, wv, (int*)bx, &cf, res, vec, flop);
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

/**
 @fn void SetBC3D::Tobc_Prdc_Simple(SklScalar3D<REAL_TYPE>* d_pt, int face)
 @brief 温度の外部周期境界条件（単純なコピー）
 @param d_t 温度のデータクラス
 @param face 面番号
 */
void SetBC3D::Tobc_Prdc_Simple(SklScalar3D<REAL_TYPE>* d_t, int face)
{
  SklParaManager* para_mng = ParaCmpo->GetParaManager();
  
  if ( para_mng->IsParallel() ) {
    unsigned no_comm_face = guide; //通信面数
    
    switch (face) {
      case X_MINUS:
        if ( !d_t->CommPeriodicBndCell(PRDC_X_DIR, PRDC_PLUS2MINUS, no_comm_face) ) Exit(0);
        break;
        
      case X_PLUS:
        if ( !d_t->CommPeriodicBndCell(PRDC_X_DIR, PRDC_MINUS2PLUS, no_comm_face) ) Exit(0);
        break;
        
      case Y_MINUS:
        if ( !d_t->CommPeriodicBndCell(PRDC_Y_DIR, PRDC_PLUS2MINUS, no_comm_face) ) Exit(0);
        break;
        
      case Y_PLUS:
        if ( !d_t->CommPeriodicBndCell(PRDC_Y_DIR, PRDC_MINUS2PLUS, no_comm_face) ) Exit(0);
        break;
        
      case Z_MINUS:
        if ( !d_t->CommPeriodicBndCell(PRDC_Z_DIR, PRDC_PLUS2MINUS, no_comm_face) ) Exit(0);
        break;
        
      case Z_PLUS:
        if ( !d_t->CommPeriodicBndCell(PRDC_Z_DIR, PRDC_MINUS2PLUS, no_comm_face) ) Exit(0);
        break;
    }
  }
  else { // Serial
    int i,j,k;
    unsigned m0, m1;
    REAL_TYPE* t=NULL;
    if ( !(t = d_t->GetData()) ) Exit(0);
    
    switch (face) {
      case X_MINUS:
        if( pn.nID[face] < 0 ) {
          for (k=1; k<=kxc; k++) {
            for (j=1; j<=jxc; j++) {
              m0 = FBUtility::getFindexS3D(size, guide, 0   , j  , k  );
              m1 = FBUtility::getFindexS3D(size, guide, ixc, j  , k  );
              t[m0] = t[m1];
            }
          }
        }
        break;
        
      case X_PLUS:
        if( pn.nID[face] < 0 ) {
          for (k=1; k<=kxc; k++) {
            for (j=1; j<=jxc; j++) {
              m0 = FBUtility::getFindexS3D(size, guide, ixc+1, j, k);
              m1 = FBUtility::getFindexS3D(size, guide, 1,      j, k);
              t[m0] = t[m1];
            }
          }
        }
        break;
        
      case Y_MINUS:
        if( pn.nID[face] < 0 ) {
          for (k=1; k<=kxc; k++) {
            for (i=1; i<=ixc; i++) {
              m0 = FBUtility::getFindexS3D(size, guide, i, 0   , k);
              m1 = FBUtility::getFindexS3D(size, guide, i, jxc, k);
              t[m0] = t[m1];
            }
          }
        }
        break;
        
      case Y_PLUS:
        if( pn.nID[face] < 0 ) {
          for (k=1; k<=kxc; k++) {
            for (i=1; i<=ixc; i++) {
              m0 = FBUtility::getFindexS3D(size, guide, i, jxc+1, k);
              m1 = FBUtility::getFindexS3D(size, guide, i, 1     , k);
              t[m0] = t[m1];
            }
          }
        }
        break;
        
      case Z_MINUS:
        if( pn.nID[face] < 0 ) {
          for (j=1; j<=jxc; j++) {
            for (i=1; i<=ixc; i++) {
              m0 = FBUtility::getFindexS3D(size, guide, i, j, 0   );
              m1 = FBUtility::getFindexS3D(size, guide, i, j, kxc);
              t[m0] = t[m1];
            }
          }
        }
        break;
        
      case Z_PLUS:
        if( pn.nID[face] < 0 ) {
          for (j=1; j<=jxc; j++) {
            for (i=1; i<=ixc; i++) {
              m0 = FBUtility::getFindexS3D(size, guide, i, j, kxc+1);
              m1 = FBUtility::getFindexS3D(size, guide, i, j, 1     );
              t[m0] = t[m1];
            }
          }
        }
        break;
    }
  }
}

/**
 @fn void SetBC3D::Vobc_Prdc_CF(SklVector3DEx<REAL_TYPE>* d_v, int face)
 @brief 速度の外部周期境界条件（単純なコピー）
 @param d_v 速度ベクトル（セルフェイス）
 @param face 面番号
 @note cbc_update_vec_cf_()でのループ範囲が[1,ix]なので，プラス方向のみ計算している
 */
void SetBC3D::Vobc_Prdc_CF(SklVector3DEx<REAL_TYPE>* d_v, int face)
{
  SklParaManager* para_mng = ParaCmpo->GetParaManager();

  if ( para_mng->IsParallel() ) {
    unsigned no_comm_face=guide;
    
    switch (face) {
      case X_MINUS:
        if ( !d_v->CommPeriodicBndCell(PRDC_X_DIR, PRDC_PLUS2MINUS, no_comm_face) ) Exit(0);
        break;
        
      case Y_MINUS:
        if ( !d_v->CommPeriodicBndCell(PRDC_Y_DIR, PRDC_PLUS2MINUS, no_comm_face) ) Exit(0);
        break;
        
      case Z_MINUS:
        if ( !d_v->CommPeriodicBndCell(PRDC_Z_DIR, PRDC_PLUS2MINUS, no_comm_face) ) Exit(0);
        break;
    }
  }
  else { // Serial
    int i,j,k, gd;
    unsigned l, m0, m1;
    gd = (int)guide;
    REAL_TYPE* v=NULL;
    if ( !(v = d_v->GetData()) ) Exit(0);
    
    switch (face) {
      case X_MINUS:
        if ( pn.nID[face] < 0 ) {
          for (k=1; k<=kxc; k++) {
            for (j=1; j<=jxc; j++) {
              for (i=1-gd; i<=0; i++) {
                for (l=0; l<3; l++) {
                  m0 = FBUtility::getFindexV3DEx(size, guide, l, i    , j, k);
                  m1 = FBUtility::getFindexV3DEx(size, guide, l, ixc+i, j, k);
                  v[m0] = v[m1];
                }
              }
            }
          }
        }
        break;
        
      case Y_MINUS:
        if ( pn.nID[face] < 0 ) {
          for (k=1; k<=kxc; k++) {
            for (j=1-gd; j<=0; j++) {
              for (i=1; i<=ixc; i++) {
                for (l=0; l<3; l++) {
                  m0 = FBUtility::getFindexV3DEx(size, guide, l, i, j    , k);
                  m1 = FBUtility::getFindexV3DEx(size, guide, l, i, jxc+j, k);
                  v[m0] = v[m1];
                }
              }
            }
          }
        }
        break;
        
      case Z_MINUS:
        if ( pn.nID[face] < 0 ) {
          for (k=1-gd; k<=0; k++) {
            for (j=1; j<=jxc; j++) {
              for (i=1; i<=ixc; i++) {
                for (l=0; l<3; l++) {
                  m0 = FBUtility::getFindexV3DEx(size, guide, l, i, j, k    );
                  m1 = FBUtility::getFindexV3DEx(size, guide, l, i, j, kxc+k);
                  v[m0] = v[m1];
                }
              }
            }
          }
        }
        break;
    }
  }
}

/**
 @fn void SetBC3D::Vobc_Prdc(SklVector3DEx<REAL_TYPE>* d_v, int face, unsigned no_comm_face)
 @brief 速度の外部周期境界条件（単純なコピー）
 @param d_v 速度ベクトル
 @param face 面番号
 @param no_comm_face 通信面数
 */
void SetBC3D::Vobc_Prdc(SklVector3DEx<REAL_TYPE>* d_v, int face, unsigned no_comm_face)
{
  SklParaManager* para_mng = ParaCmpo->GetParaManager();
  
  if ( para_mng->IsParallel() ) {
    unsigned no_comm_face=guide;
    
    switch (face) {
      case X_MINUS:
        if ( !d_v->CommPeriodicBndCell(PRDC_X_DIR, PRDC_PLUS2MINUS, no_comm_face) ) Exit(0);
        break;
        
      case X_PLUS:
        if ( !d_v->CommPeriodicBndCell(PRDC_X_DIR, PRDC_MINUS2PLUS, no_comm_face) ) Exit(0);
        break;
        
      case Y_MINUS:
        if ( !d_v->CommPeriodicBndCell(PRDC_Y_DIR, PRDC_PLUS2MINUS, no_comm_face) ) Exit(0);
        break;
        
      case Y_PLUS:
        if ( !d_v->CommPeriodicBndCell(PRDC_Y_DIR, PRDC_MINUS2PLUS, no_comm_face) ) Exit(0);
        break;
        
      case Z_MINUS:
        if ( !d_v->CommPeriodicBndCell(PRDC_Z_DIR, PRDC_PLUS2MINUS, no_comm_face) ) Exit(0);
        break;
        
      case Z_PLUS:
        if ( !d_v->CommPeriodicBndCell(PRDC_Z_DIR, PRDC_MINUS2PLUS, no_comm_face) ) Exit(0);
        break;
    }
  }
  else { // Serial
    int i,j,k, gd;
    unsigned l, m0, m1;
    gd = (int)guide;
    REAL_TYPE* v=NULL;
    if ( !(v = d_v->GetData()) ) Exit(0);
    
    switch (face) {
      case X_MINUS:
        if ( pn.nID[face] < 0 ) {
          for (k=1; k<=kxc; k++) {
            for (j=1; j<=jxc; j++) {
              for (i=1-gd; i<=0; i++) {
                for (l=0; l<3; l++) {
                  m0 = FBUtility::getFindexV3DEx(size, guide, l, i    , j, k);
                  m1 = FBUtility::getFindexV3DEx(size, guide, l, ixc+i, j, k);
                  v[m0] = v[m1];
                }
              }
            }
          }
        }
        break;
        
      case X_PLUS:
        if ( pn.nID[face] < 0 ) {
          for (k=1; k<=kxc; k++) {
            for (j=1; j<=jxc; j++) {
              for (i=ixc+1; i<=ixc+gd; i++) {
                for (l=0; l<3; l++) {
                  m0 = FBUtility::getFindexV3DEx(size, guide, l, i    , j, k);
                  m1 = FBUtility::getFindexV3DEx(size, guide, l, i-ixc, j, k);
                  v[m0] = v[m1];
                }
              }
            }
          }
        }
        break;
        
      case Y_MINUS:
        if ( pn.nID[face] < 0 ) {
          for (k=1; k<=kxc; k++) {
            for (j=1-gd; j<=0; j++) {
              for (i=1; i<=ixc; i++) {
                for (l=0; l<3; l++) {
                  m0 = FBUtility::getFindexV3DEx(size, guide, l, i, j    , k);
                  m1 = FBUtility::getFindexV3DEx(size, guide, l, i, jxc+j, k);
                  v[m0] = v[m1];
                }
              }
            }
          }
        }
        break;
        
      case Y_PLUS:
        if ( pn.nID[face] < 0 ) {
          for (k=1; k<=kxc; k++) {
            for (j=jxc+1; j<=jxc+gd; j++) {
              for (i=1; i<=ixc; i++) {
                for (l=0; l<3; l++) {
                  m0 = FBUtility::getFindexV3DEx(size, guide, l, i, j    , k);
                  m1 = FBUtility::getFindexV3DEx(size, guide, l, i, j-jxc, k);
                  v[m0] = v[m1];
                }
              }
            }
          }
        }
        break;
        
      case Z_MINUS:
        if ( pn.nID[face] < 0 ) {
          for (k=1-gd; k<=0; k++) {
            for (j=1; j<=jxc; j++) {
              for (i=1; i<=ixc; i++) {
                for (l=0; l<3; l++) {
                  m0 = FBUtility::getFindexV3DEx(size, guide, l, i, j, k    );
                  m1 = FBUtility::getFindexV3DEx(size, guide, l, i, j, kxc+k);
                  v[m0] = v[m1];
                }
              }
            }
          }
        }
        break;
        
      case Z_PLUS:
        if ( pn.nID[face] < 0 ) {
          for (k=kxc+1; k<=kxc+gd; k++) {
            for (j=1; j<=jxc; j++) {
              for (i=1; i<=ixc; i++) {
                for (l=0; l<3; l++) {
                  m0 = FBUtility::getFindexV3DEx(size, guide, l, i, j, k    );
                  m1 = FBUtility::getFindexV3DEx(size, guide, l, i, j, k-kxc);
                  v[m0] = v[m1];
                }
              }
            }
          }
        }
        break;
    }
  }
}

/**
 @fn void SetBC3D::Vibc_Prdc(SklVector3DEx<REAL_TYPE>* d_v, int* st, int* ed, SklScalar3D<unsigned>* d_bd, int odr, int dir)
 @brief 速度の内部周期境界条件（単純なコピー）
 @param d_v 速度ベクトル
 @param st コンポーネント範囲の開始インデクス
 @param ed コンポーネント範囲の終了インデクス
 @param d_bd BCindex ID
 @param odr 下流面のidが格納されているコンポーネントエントリ
 @param dir ドライバの設置方向
 */
void SetBC3D::Vibc_Prdc(SklVector3DEx<REAL_TYPE>* d_v, int* st, int* ed, SklScalar3D<unsigned>* d_bd, int odr, int dir)
{
  SklParaManager* para_mng = ParaCmpo->GetParaManager();
  if( para_mng->IsParallel() ){
    Hostonly_ printf("Error : 'Vibc_Prdc' method is limited to use for serial execution\n.");
    Exit(0);
  }
    
  int i,j,k, gd, ii, jj, kk;
  unsigned l, m0, m1, m2, *bx=NULL;
  REAL_TYPE* v=NULL;

  gd = (int)guide;
  
  if ( !(v = d_v->GetData()) )  Exit(0);
  if ( !(bx= d_bd->GetData()) ) Exit(0);
  
  switch (dir) {
    case X_MINUS:
      if ( pn.nID[dir] < 0 ) {
        i = st[0];
        for (k=st[2]; k<=ed[2]; k++) {
          for (j=st[1]; j<=ed[1]; j++) {
            for (l=0; l<3; l++) {
              m1 = FBUtility::getFindexS3D(size, guide, i, j, k);
              if ( DECODE_CMP(bx[m1]) == odr ) {
                for (ii=1-gd; ii<=0; ii++) {
                  m0 = FBUtility::getFindexV3DEx(size, guide, l, ii,   j, k);
                  m2 = FBUtility::getFindexV3DEx(size, guide, l, i+ii, j, k);
                  v[m0] = v[m2];
                }
              }
            }
          }
        }
      }
      break;
      
    case X_PLUS:
      if ( pn.nID[dir] < 0 ) {
        i = st[0];
        for (k=st[2]; k<=ed[2]; k++) {
          for (j=st[1]; j<=ed[1]; j++) {
            for (l=0; l<3; l++) {
              m1 = FBUtility::getFindexS3D(size, guide, i, j, k);
              if ( DECODE_CMP(bx[m1]) == odr ) {
                for (ii=ixc+1; ii<=ixc+gd; ii++) {
                  m0 = FBUtility::getFindexV3DEx(size, guide, l, ii,         j, k);
                  m2 = FBUtility::getFindexV3DEx(size, guide, l, i+ii-ixc-1, j, k);
                  v[m0] = v[m2];
                }
              }              
            }
          }
        }
      }
      break;
      
    case Y_MINUS:
      if ( pn.nID[dir] < 0 ) {
        j = st[1];
        for (k=st[2]; k<=ed[2]; k++) {
          for (i=st[0]; i<=ed[0]; i++) {
            for (l=0; l<3; l++) {
              m1 = FBUtility::getFindexS3D(size, guide, i, j, k);
              if ( DECODE_CMP(bx[m1]) == odr ) {
                for (jj=1-gd; jj<=0; jj++) {
                  m0 = FBUtility::getFindexV3DEx(size, guide, l, i, jj,   k);
                  m2 = FBUtility::getFindexV3DEx(size, guide, l, i, j+jj, k);
                  v[m0] = v[m2];
                }
              }              
            }
          }
        }
      }
      break;
      
    case Y_PLUS:
      if ( pn.nID[dir] < 0 ) {
        j = st[1];
        for (k=st[2]; k<=ed[2]; k++) {
          for (i=st[0]; i<=ed[0]; i++) {
            for (l=0; l<3; l++) {
              m1 = FBUtility::getFindexS3D(size, guide, i, j, k);
              if ( DECODE_CMP(bx[m1]) == odr ) {
                for (jj=jxc+1; jj<=jxc+gd; jj++) {
                  m0 = FBUtility::getFindexV3DEx(size, guide, l, i, jj,         k);
                  m2 = FBUtility::getFindexV3DEx(size, guide, l, i, j+jj-jxc-1, k);
                  v[m0] = v[m2];
                }
              }              
            }
          }
        }
      }
      break;
      
    case Z_MINUS:
      if ( pn.nID[dir] < 0 ) {
        k = st[2];
        for (j=st[1]; j<=ed[1]; j++) {
          for (i=st[0]; i<=ed[0]; i++) {
            for (l=0; l<3; l++) {
              m1 = FBUtility::getFindexS3D(size, guide, i, j, k);
              if ( DECODE_CMP(bx[m1]) == odr ) {
                for (kk=1-gd; kk<=0; kk++) {
                  m0 = FBUtility::getFindexV3DEx(size, guide, l, i, j, kk  );
                  m2 = FBUtility::getFindexV3DEx(size, guide, l, i, j, k+kk);
                  v[m0] = v[m2];
                }
              }         
            }
          }
        }
      }
      break;
      
    case Z_PLUS:
      if ( pn.nID[dir] < 0 ) {
        k = st[2];
        for (j=st[1]; j<=ed[1]; j++) {
          for (i=st[0]; i<=ed[0]; i++) {
            for (l=0; l<3; l++) {
              m1 = FBUtility::getFindexS3D(size, guide, i, j, k);
              if ( DECODE_CMP(bx[m1]) == odr ) {
                for (kk=kxc+1; kk<=kxc+gd; kk++) {
                  m0 = FBUtility::getFindexV3DEx(size, guide, l, i, j, kk        );
                  m2 = FBUtility::getFindexV3DEx(size, guide, l, i, j, k+kk-kxc-1);
                  v[m0] = v[m2];
                }
              }         
            }
          }
        }
      }
      break;
  }
}

