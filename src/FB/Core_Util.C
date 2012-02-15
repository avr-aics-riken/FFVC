/*
 * SPHERE - Skeleton for PHysical and Engineering REsearch
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2004-2012
 *
 */

//@file Core_Util.C
//@brief Core_Utility class
//@author keno, FSI Team, VCAD, RIKEN

#include "Core_Util.h"
extern SklParaComponent* ParaCmpo;

/**
 @fn bool Core_Utility::shiftVin3D(SklVector3DEx<SKL_REAL>* dst, const SklVector3DEx<SKL_REAL>* src, SKL_REAL v00[3], unsigned stepAvr)
 @brief srcのデータにある速度成分v00[3]を加えdstにコピー
 @param dst       コピー先データ
 @param src       コピー元データ
 @param v00       減算値
 @param stepAvr   積算値
 @return    true=success, false=fault
 @note dst[index] = ( src[index] + v00 ) * stepAvr
 */
bool Core_Utility::shiftVin3D(SklVector3DEx<SKL_REAL>* dst, const SklVector3DEx<SKL_REAL>* src, SKL_REAL v00[3], unsigned stepAvr)
{
  if( !dst || !src || !v00 ) return false;
  
  const unsigned* dst_sz = dst->GetSize();  // dimension size /wo guide cell
  const unsigned* src_sz = src->GetSize();
  if(   (src_sz[0] != dst_sz[0])
     || (src_sz[1] != dst_sz[1])
     || (src_sz[2] != dst_sz[2]) ) return false;
  
  SKL_REAL* dst_data = dst->GetData();
  const SKL_REAL* src_data = src->GetData();
  if( !dst_data || !src_data ) return false;
  
  dst_sz = dst->_GetSize();  // dimension size /w guide cell
  src_sz = src->_GetSize();
  unsigned dst_gc = dst->GetVCellSize();
  unsigned src_gc = src->GetVCellSize();
  unsigned sta, ix, jx, kx;
  int diff;
  
  CalcIndex(dst_sz[0], dst_sz[1], dst_sz[2], dst_gc, src_gc, ix, jx, kx, diff, sta);
  
  SKL_REAL  ra = 1.0, va[3];
  ra=(SKL_REAL)stepAvr;     // input
  for (int i=0; i<3; i++) va[i] = v00[i]*ra;
  
  unsigned long idx, lsz[3];
  lsz[0] = src_sz[0];
  lsz[1] = src_sz[1];
  lsz[2] = src_sz[2];
  
  register unsigned i, j, k;
  
  for(k=sta; k<kx; k++){
    unsigned kk = k+diff;
    for(j=sta; j<jx; j++){
      unsigned jj = j+diff;
      for(i=sta; i<ix; i++){
        // idx : indecies for SklVector3DEx
        idx = 3*(lsz[0]*lsz[1]*kk + lsz[0]*jj + i+diff);
        dst_data[idx  ] = src_data[idx  ]*ra + va[0];
        dst_data[idx+1] = src_data[idx+1]*ra + va[1];
        dst_data[idx+2] = src_data[idx+2]*ra + va[2];
      }
    }
  }
  
  return true;
}

/**
 @fn bool Core_Utility::shiftVout3D(SklVector3DEx<SKL_REAL>* dst, const SklVector3DEx<SKL_REAL>* src, SKL_REAL v00[3], unsigned stepAvr)
 @brief shiftVin3D()の逆演算
 @param dst       コピー先データ
 @param src       コピー元データ
 @param v00       減算値
 @param stepAvr   積算値
 @return    true=success, false=fault
 @note dst[index] = ( src[index] * stepAvr ) - v00
 */
bool Core_Utility::shiftVout3D(SklVector3DEx<SKL_REAL>* dst, const SklVector3DEx<SKL_REAL>* src, SKL_REAL v00[3], unsigned stepAvr)
{
  if( !dst || !src || !v00 ) return false;
  
  const unsigned* dst_sz = dst->GetSize();  // dimension size /wo guide cell
  const unsigned* src_sz = src->GetSize();
  if(   (src_sz[0] != dst_sz[0])
     || (src_sz[1] != dst_sz[1])
     || (src_sz[2] != dst_sz[2]) ) return false;
  
  SKL_REAL* dst_data = dst->GetData();
  const SKL_REAL* src_data = src->GetData();
  if( !dst_data || !src_data ) return false;
  
  dst_sz = dst->_GetSize();  // dimension size /w guide cell
  src_sz = src->_GetSize();
  unsigned dst_gc = dst->GetVCellSize();
  unsigned src_gc = src->GetVCellSize();
  unsigned sta, ix, jx, kx;
  int diff;
  
  CalcIndex(dst_sz[0], dst_sz[1], dst_sz[2], dst_gc, src_gc, ix, jx, kx, diff, sta);
  
  SKL_REAL  ra = 1.0;
  ra=1.0/(SKL_REAL)stepAvr; // output
  
  unsigned long idx, lsz[3];
  lsz[0] = src_sz[0];
  lsz[1] = src_sz[1];
  lsz[2] = src_sz[2];
  
  register unsigned i, j, k;
  
  for(k=sta; k<kx; k++){
    unsigned kk = k+diff;
    for(j=sta; j<jx; j++){
      unsigned jj = j+diff;
      for(i=sta; i<ix; i++){
        // idx : indecies for SklVector3DEx
        idx = 3*(lsz[0]*lsz[1]*k + lsz[0]*j + i);
        dst_data[idx  ] = src_data[idx  ]*ra - v00[0];
        dst_data[idx+1] = src_data[idx+1]*ra - v00[1];
        dst_data[idx+2] = src_data[idx+2]*ra - v00[2];
      }
    }
  }
  
  return true;
}

/**
 @fn bool Core_Utility::shiftVin3D(SklVector3D<SKL_REAL>* dst, const SklVector3DEx<SKL_REAL>* src, SKL_REAL v00[3], unsigned stepAvr)
 @brief srcのデータにある速度成分v00[3]を加えdstにコピー
 @param dst       コピー先データ
 @param src       コピー元データ
 @param v00       減算値
 @param stepAvr   積算値
 @return    true=success, false=fault
 @note dst[index] = ( src[index] + v00 ) * stepAvr
 */
bool Core_Utility::shiftVin3D(SklVector3D<SKL_REAL>* dst, const SklVector3DEx<SKL_REAL>* src, SKL_REAL v00[3], unsigned stepAvr)
{
  if( !dst || !src || !v00 ) return false;
  
  const unsigned* dst_sz = dst->GetSize();  // dimension size /wo guide cell
  const unsigned* src_sz = src->GetSize();
  if(   (src_sz[0] != dst_sz[0])
     || (src_sz[1] != dst_sz[1])
     || (src_sz[2] != dst_sz[2]) ) return false;
  
  SKL_REAL* dst_data = dst->GetData();
  const SKL_REAL* src_data = src->GetData();
  if( !dst_data || !src_data ) return false;
  
  dst_sz = dst->_GetSize();  // dimension size /w guide cell
  src_sz = src->_GetSize();
  unsigned dst_gc = dst->GetVCellSize();
  unsigned src_gc = src->GetVCellSize();
  unsigned sta, ix, jx, kx;
  int diff;
  
  CalcIndex(dst_sz[0], dst_sz[1], dst_sz[2], dst_gc, src_gc, ix, jx, kx, diff, sta);
  
  SKL_REAL  ra = 1.0, va[3];
  ra=(SKL_REAL)stepAvr;     // input
  for (int i=0; i<3; i++) va[i] = v00[i]*ra;
  
  unsigned long src_idx, dst_idx, src_lsz[3], dst_lsz[3];
  src_lsz[0] = src_sz[0];
  src_lsz[1] = src_sz[1];
  src_lsz[2] = src_sz[2];
  dst_lsz[0] = dst_sz[0];
  dst_lsz[1] = dst_sz[1];
  dst_lsz[2] = dst_sz[2];
  
  unsigned long s1 = dst_sz[0]*dst_sz[1]*dst_sz[2];
  unsigned long s2 = s1*2;
  register unsigned i, j, k;
  
  for(k=sta; k<kx; k++){
    unsigned kk = k+diff;
    for(j=sta; j<jx; j++){
      unsigned jj = j+diff;
      for(i=sta; i<ix; i++){
        // src_idx : indecies for SklVector3DEx
        src_idx = 3*(src_lsz[0]*src_lsz[1]*kk + src_lsz[0]*jj + i+diff);
        // dst_idx : indecies for SklVector3D
        dst_idx = dst_lsz[0]*dst_lsz[1]*k + dst_lsz[0]*j + i;
        dst_data[dst_idx   ] = src_data[src_idx  ]*ra + va[0];
        dst_data[dst_idx+s1] = src_data[src_idx+1]*ra + va[1];
        dst_data[dst_idx+s2] = src_data[src_idx+2]*ra + va[2];
      }
    }
  }
  
  return true;
}

/**
 @fn bool Core_Utility::shiftVout3D(SklVector3DEx<SKL_REAL>* dst, const SklVector3D<SKL_REAL>* src, SKL_REAL v00[3], unsigned stepAvr)
 @brief shiftVin3D()の逆演算
 @param dst       コピー先データ
 @param src       コピー元データ
 @param v00       減算値
 @param stepAvr   積算値
 @return    true=success, false=fault
 @note dst[index] = ( src[index] * stepAvr ) - v00
 */
bool Core_Utility::shiftVout3D(SklVector3DEx<SKL_REAL>* dst, const SklVector3D<SKL_REAL>* src, SKL_REAL v00[3], unsigned stepAvr)
{
  if( !dst || !src || !v00 ) return false;
  
  const unsigned* dst_sz = dst->GetSize();  // dimension size /wo guide cell
  const unsigned* src_sz = src->GetSize();
  if(   (src_sz[0] != dst_sz[0])
     || (src_sz[1] != dst_sz[1])
     || (src_sz[2] != dst_sz[2]) ) return false;
  
  SKL_REAL* dst_data = dst->GetData();
  const SKL_REAL* src_data = src->GetData();
  if( !dst_data || !src_data ) return false;
  
  dst_sz = dst->_GetSize();  // dimension size /w guide cell
  src_sz = src->_GetSize();
  unsigned dst_gc = dst->GetVCellSize();
  unsigned src_gc = src->GetVCellSize();
  unsigned sta, ix, jx, kx;
  int diff;
  
  CalcIndex(dst_sz[0], dst_sz[1], dst_sz[2], dst_gc, src_gc, ix, jx, kx, diff, sta);
  
  SKL_REAL  ra = 1.0;
  ra=1.0/(SKL_REAL)stepAvr; // output
  
  unsigned long src_idx, dst_idx, src_lsz[3], dst_lsz[3];
  src_lsz[0] = src_sz[0];
  src_lsz[1] = src_sz[1];
  src_lsz[2] = src_sz[2];
  dst_lsz[0] = dst_sz[0];
  dst_lsz[1] = dst_sz[1];
  dst_lsz[2] = dst_sz[2];
  
  unsigned long s1 = src_sz[0]*src_sz[1]*src_sz[2];
  unsigned long s2 = s1*2;
  register unsigned i, j, k;
  
  for(k=sta; k<kx; k++){
    unsigned kk = k+diff;
    for(j=sta; j<jx; j++){
      unsigned jj = j+diff;
      for(i=sta; i<ix; i++){
        // src_idx : indecies for SklVector3D
        src_idx = src_lsz[0]*src_lsz[1]*kk + src_lsz[0]*jj + (i+diff);
        // dst_idx : indecies for SklVector3DEx
        dst_idx = 3*(dst_lsz[0]*dst_lsz[1]*k + dst_lsz[0]*j + i);
        dst_data[dst_idx  ] = src_data[src_idx   ]*ra - v00[0];
        dst_data[dst_idx+1] = src_data[src_idx+s1]*ra - v00[1];
        dst_data[dst_idx+2] = src_data[src_idx+s2]*ra - v00[2];
      }
    }
  }
  
  return true;
}

/**
 @fn SKL_REAL Core_Utility::norm_v_div_max(unsigned sz[3], unsigned guide, SKL_REAL coef, SKL_REAL* src, unsigned* bp, SKL_REAL& flop) const
 @brief 有効セルに対する発散の最大値
 @retval 発散の絶対値の最大値
 @param sz 配列サイズ
 @param guide ガイドセル
 @param coef 発散値の係数
 @param src 発散値のベース
 @param bp BCindexP
 @param flop[out] 浮動小数演算数
 */
SKL_REAL Core_Utility::norm_v_div_max(unsigned sz[3], unsigned guide, SKL_REAL coef, SKL_REAL* src, unsigned* bp, SKL_REAL& flop) const
{
  int i,j,k;
  unsigned register m;
  SKL_REAL r;  // 速度の発散値
  SKL_REAL d;  // 発散の絶対値
  SKL_REAL ds; // 最大値
  flop += (SKL_REAL)(sz[0]*sz[1]*sz[2]*4);
  
  ds = 0.0;
  
  for (k=1; k<=(int)sz[2]; k++) {
    for (j=1; j<=(int)sz[1]; j++) {
      for (i=1; i<=(int)sz[0]; i++) {
        m = SklUtil::getFindexS3D(sz, guide, i  , j  , k  );
        r = src[m] * coef * GET_SHIFT_F(bp[m], VLD_CNVG); // 有効セルの場合 1.0
        d = fabs(r);
        if ( d > ds ) ds = d;
      }
    }
  }
  return ds;
}

/**
 @fn SKL_REAL Core_Utility::norm_v_div_l2(unsigned sz[3], unsigned guide, SKL_REAL coef, SKL_REAL* src, unsigned* bp, SKL_REAL& flop) const
 @brief 有効セルに対する発散の自乗和を計算
 @retval 発散値の自乗和
 @param sz 配列サイズ
 @param guide ガイドセル
 @param coef 発散値の係数
 @param src 発散値のベース
 @param bp BCindexP
 @param flop[out] 浮動小数演算数
 */
SKL_REAL Core_Utility::norm_v_div_l2(unsigned sz[3], unsigned guide, SKL_REAL coef, SKL_REAL* src, unsigned* bp, SKL_REAL& flop) const
{
  int i,j,k;
  unsigned register m;
  SKL_REAL r;  // 速度の発散値
  SKL_REAL w;  // 発散の自乗和
  flop += (SKL_REAL)(sz[0]*sz[1]*sz[2]*4);

  w = 0.0;
  
  for (k=1; k<=(int)sz[2]; k++) {
    for (j=1; j<=(int)sz[1]; j++) {
      for (i=1; i<=(int)sz[0]; i++) {
        m = SklUtil::getFindexS3D(sz, guide, i  , j  , k  );
        r = src[m]*coef* GET_SHIFT_F(bp[m], VLD_CNVG); // 有効セルの場合 1.0
        w += r*r;
      }
    }
  }
  return w;
}

/**
 @fn void Core_Utility::norm_v_div_dbg(SKL_REAL& nrm, SKL_REAL& rm, int* index, unsigned sz[3], unsigned guide, SKL_REAL coef, SKL_REAL* src, unsigned* bp, SKL_REAL& flop) const
 @brief 有効セルに対する発散の最大値と自乗和を計算，絶対値の最大値の位置を返す
 @param nrm[out] 残差の絶対値
 @param rm [out] 残差の自乗和
 @param index[out] ノルムの最大値位置
 @param sz 配列サイズ
 @param guide ガイドセル
 @param coef 発散値の係数
 @param src 発散値のベース
 @param bp BCindexP
 @param flop[out] 浮動小数演算数
 */
void Core_Utility::norm_v_div_dbg(SKL_REAL& nrm, SKL_REAL& rm, int* index, unsigned sz[3], unsigned guide, SKL_REAL coef, SKL_REAL* src, unsigned* bp, SKL_REAL& flop) const
{
  int i,j,k, i0, j0, k0;
  unsigned register m;
  SKL_REAL r;  // 速度の発散値
  SKL_REAL d;  // 発散の絶対値
  SKL_REAL w;  // 発散の自乗和
  SKL_REAL ds; // 最大値
  flop += (SKL_REAL)(sz[0]*sz[1]*sz[2]*6);
  
  i0 = j0 = k0 = 0;
  d = r = w = ds = 0.0;
  
  for (k=1; k<=(int)sz[2]; k++) {
    for (j=1; j<=(int)sz[1]; j++) {
      for (i=1; i<=(int)sz[0]; i++) {
        m = SklUtil::getFindexS3D(sz, guide, i  , j  , k  );
        r = src[m]*coef* GET_SHIFT_F(bp[m], VLD_CNVG); // 有効セルの場合 1.0
        d = fabs(r);
        w += r*r;
        
        if ( d > ds ) {
          i0 = i;
          j0 = j;
          k0 = k;
          ds = d;
        }
      }
    }
  }
  
  nrm = ds;
  rm  = w;
  index[0] = i0;
  index[1] = j0;
  index[2] = k0;
}

/**
 @fn SKL_REAL Core_Utility::count_comm_size(unsigned sz[3], unsigned guide) const
 @brief 全ノードについて，ローカルノード1面・一層あたりの通信量の和を返す
 @retval 通信量(Byte)
 @param sz 配列サイズ
 @param guide ガイドセル
 */
SKL_REAL Core_Utility::count_comm_size(unsigned sz[3], unsigned guide) const
{
  SklParaManager* para_mng = ParaCmpo->GetParaManager();
  SKL_REAL c = 0.0;
  
  // 内部面のみをカウントする
  for (unsigned n=0; n<6; n++) {
    if ( pn.nID[n] >= 0 ) {
      
      switch (n) {
        case X_MINUS:
        case X_PLUS:
          c += (SKL_REAL)(sz[1]*sz[2]);
          break;
          
        case Y_MINUS:
        case Y_PLUS:
          c += (SKL_REAL)(sz[0]*sz[2]);
          break;
          
        case Z_MINUS:
        case Z_PLUS:
          c += (SKL_REAL)(sz[0]*sz[1]);
          break;
      }
    }
  }
  
  if( para_mng->IsParallel() ){
    SKL_REAL tmp = c;
    if ( !para_mng->Allreduce(&tmp, &c, 1, SKL_ARRAY_DTYPE_REAL, SKL_SUM, pn.procGrp) ) assert(0);
  }
  
  return c*(SKL_REAL)sizeof(SKL_REAL); // Byte
}

/**
 @fn void Core_Utility::delta_Scalar(unsigned sz[3], unsigned guide, SKL_REAL* sn, SKL_REAL* so, unsigned* bx, SKL_REAL* var, SKL_REAL& flop) const
 @brief 有効セルに対する，1タイムステップ進行時の変化量の2乗和と平均値
 @param sz 配列サイズ
 @param guide ガイドセル
 @param sn スカラー n+1ステップ
 @param so スカラー nステップ
 @param bx BCindex
 @param var[out] 平均値と変動量
 @param flop 浮動小数演算数
 */
void Core_Utility::delta_Scalar(unsigned sz[3], unsigned guide, SKL_REAL* sn, SKL_REAL* so, unsigned* bx, SKL_REAL* var, SKL_REAL& flop) const
{
  int i,j,k;
  unsigned register m;
  SKL_REAL rm, av, a, b;
  flop += (SKL_REAL)(sz[0]*sz[1]*sz[2]*6);
  
  rm = av = 0.0;
  for (k=1; k<=(int)sz[2]; k++) {
    for (j=1; j<=(int)sz[1]; j++) {
      for (i=1; i<=(int)sz[0]; i++) {
        m = SklUtil::getFindexS3D(sz, guide, i  , j  , k  );
        b = GET_SHIFT_F(bx[m], ACTIVE_BIT);
        av += sn[m]*b;
        a = (sn[m] - so[m]) * b;
        rm += a*a;
      }
    }
  }
  var[0] = rm;
  var[1] = av;
}

/**
 @fn void Core_Utility::delta_Vector(unsigned sz[3], unsigned guide, SKL_REAL* vn, SKL_REAL* vo, unsigned* bx, SKL_REAL* var, SKL_REAL& flop) const
 @brief 有効セルに対する，1タイムステップ進行時の変化量の2乗和と平均値(RootMean)
 @param sz 配列サイズ
 @param guide ガイドセル
 @param vn ベクトル値 n+1ステップ
 @param vo ベクトル値 nステップ
 @param bx BCindex
 @param var[out] 平均値と変動量
 @param flop 浮動小数演算数
 */
void Core_Utility::delta_Vector(unsigned sz[3], unsigned guide, SKL_REAL* vn, SKL_REAL* vo, unsigned* bx, SKL_REAL* var, SKL_REAL& flop) const
{
  int i,j,k;
  unsigned register m0, m1, m2, m;
  SKL_REAL rm, x, y, z, av, b, u, v, w;
  flop += (SKL_REAL)(sz[0]*sz[1]*sz[2]*17);
  
  rm = av = 0.0;
  for (k=1; k<=(int)sz[2]; k++) {
    for (j=1; j<=(int)sz[1]; j++) {
      for (i=1; i<=(int)sz[0]; i++) {
        m  = SklUtil::getFindexS3D(sz, guide, i, j, k);
        m0 = SklUtil::getFindexV3DEx(sz, guide, 0, i, j, k);
        m1 = SklUtil::getFindexV3DEx(sz, guide, 1, i, j, k);
        m2 = SklUtil::getFindexV3DEx(sz, guide, 2, i, j, k);
        b = GET_SHIFT_F(bx[m], ACTIVE_BIT);
        u = vn[m0];
        v = vn[m1];
        w = vn[m2];
        av += (u*u + v*v + w*w)*b;
        x = u - vo[m0];
        y = v - vo[m1];
        z = w - vo[m2];
        rm += (x*x + y*y + z*z)*b;
      }
    }
  }
  
  var[0] = rm;
  var[1] = av;
}

/**
 @fn void Core_Utility::CutOffRange(SKL_REAL* t, CompoList* cmp, unsigned* bx, Control* C) const
 @brief 値のフィルタ
 @param t スカラ値
 @param cmp コンポーネント
 @param bx BCIndex
 @param C Controlクラス
 @note 隠しパラメータに
 */
void Core_Utility::CutOffRange(SKL_REAL* t, CompoList* cmp, unsigned* bx, Control* C) const
{
  unsigned m_size[3], guide, st[3], ed[3];
  
  m_size[0] = C->imax;
  m_size[1] = C->jmax;
  m_size[2] = C->kmax;
  guide     = C->guide;
  
  switch ( C->Hide.Range_Limit ) {
    case Control::Range_Normal:
      break;
      
    case Control::Range_Cutoff:  // this case includes cutoff for suction
      fb_limit_scalar_(t, (int*)m_size, (int*)&guide);
      break;
  }
}

/**
 @fn void Core_Utility::TotalPressure (SKL_REAL* tp, SKL_REAL* v, SKL_REAL* p, Control* C, const char* arrangement, SKL_REAL* v00, SKL_REAL& flop) const
 @brief 全圧を計算する
 @param tp 全圧
 @param v 速度ベクトル
 @param p 圧力
 @param C Controlクラス
 @param arrangement 格子配置
 @param v00 参照速度
 @param flop 浮動小数演算数
 */
void Core_Utility::TotalPressure (SKL_REAL* tp, SKL_REAL* v, SKL_REAL* p, Control* C, const char* arrangement, SKL_REAL* v00, SKL_REAL& flop) const
{
  int flag;
  unsigned m_size[3], guide;
  
  m_size[0] = C->imax;
  m_size[1] = C->jmax;
  m_size[2] = C->kmax;
  guide     = C->guide;
  
  if ( !strcasecmp("staggered", arrangement)) {
    flag = 0;
  }
  else if ( !strcasecmp("collocated", arrangement)) {
    flag = 1;
  }
  else {
    assert(0);
  }
  
  fb_totalp_ (tp, (int*)m_size, (int*)&guide, v, p, v00, &flop);
}
