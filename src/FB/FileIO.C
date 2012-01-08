/*
 * SPHERE - Skeleton for PHysical and Engineering REsearch
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2004-2011
 *
 */

//@file FileIO.C
//@brief FlowBase FileIO class Header
//@author keno, FSI Team, VCAD, RIKEN

#include "FileIO.h"
extern SklParaComponent* ParaCmpo;

/**
 @fn void FileIO::cnv_Div(SklScalar3D<SKL_REAL>* dst, const SklScalar3D<SKL_REAL>* src, const SKL_REAL coef, SKL_REAL& flop) const
 @brief ファイル出力時，発散値を計算する
 @param dst 単位変換のデータクラス
 @param src 単位変換前のデータクラス
 @param coef 係数
 @param flop 浮動小数演算数
 @see SklUtil::cpyS3D(), Core_Utility::norm_v_div_max()
 */
void FileIO::cnv_Div(SklScalar3D<SKL_REAL>* dst, const SklScalar3D<SKL_REAL>* src, const SKL_REAL coef, SKL_REAL& flop) const
{
  if( !dst || !src ) assert(0);
  
  const unsigned* dst_sz = dst->GetSize();
  const unsigned* src_sz = src->GetSize();
  if(   (src_sz[0] != dst_sz[0])
     || (src_sz[1] != dst_sz[1])
     || (src_sz[2] != dst_sz[2]) ) assert(0);
  
  SKL_REAL* dst_data = dst->GetData();
  const SKL_REAL* src_data = src->GetData();
  if( !dst_data || !src_data ) assert(0);
  
  dst_sz = dst->_GetSize();
  src_sz = src->_GetSize();
  unsigned dst_gc = dst->GetVCellSize();
  unsigned src_gc = src->GetVCellSize();
  unsigned sta, ix, jx, kx;
  int diff;
  register unsigned i, j, k;
  
  unsigned long dst_lsz[3];
  dst_lsz[0] = dst_sz[0];
  dst_lsz[1] = dst_sz[1];
  dst_lsz[2] = dst_sz[2];
  
  // compare whole size
  if(   (src_sz[0] == dst_sz[0])
     && (src_sz[1] == dst_sz[1])
     && (src_sz[2] == dst_sz[2]) ) {
    unsigned long idx;
    
    for(k=0; k<dst_sz[2]; k++){
      for(j=0; j<dst_sz[1]; j++){
        for(i=0; i<dst_sz[0]; i++){
          idx = dst_lsz[0]*dst_lsz[1]*k  + dst_lsz[0]*j  + i;
          dst_data[idx] = src_data[idx]*coef;
        }
      }
    }
    flop += (SKL_REAL)(dst_sz[0]*dst_sz[1]*dst_sz[2]);
  }
  else {
    assert(0);
  }
}

/**
 @fn void FileIO::cnv_V_ND2D(SklVector3DEx<SKL_REAL>* dst, const SklVector3D<SKL_REAL>* src, const SKL_REAL v00[3], const SKL_REAL Ref_v, 
                             SKL_REAL& flop, unsigned stepAvr) const
 @brief 速度データについて，無次元から有次元単位に変換
 @param dst 単位変換後のデータクラス
 @param src 単位変換前のデータクラス
 @param v00 参照速度
 @param Ref_v 代表速度(m/s)
 @param flop 浮動小数演算数
 @param stepAvr 時間平均をとったステップ数
 @see Core_Util::shiftVout3D()
 */
void FileIO::cnv_V_ND2D(SklVector3DEx<SKL_REAL>* dst, const SklVector3D<SKL_REAL>* src, const SKL_REAL v00[3], const SKL_REAL Ref_v, 
                        SKL_REAL& flop, unsigned stepAvr) const
{
  if( !dst || !src || !v00 ) assert(0);
  
  const unsigned* dst_sz = dst->GetSize();  // dimension size /wo guide cell
  const unsigned* src_sz = src->GetSize();
  if(   (src_sz[0] != dst_sz[0])
     || (src_sz[1] != dst_sz[1])
     || (src_sz[2] != dst_sz[2]) ) assert(0);
  
  SKL_REAL* dst_data = dst->GetData();
  const SKL_REAL* src_data = src->GetData();
  if( !dst_data || !src_data ) assert(0);
  
  dst_sz = dst->_GetSize();  // dimension size /w guide cell
  src_sz = src->_GetSize();
  unsigned dst_gc = dst->GetVCellSize();
  unsigned src_gc = src->GetVCellSize();
  unsigned sta, ix, jx, kx;
  int diff;
  
  Core_Utility::CalcIndex(dst_sz[0], dst_sz[1], dst_sz[2], dst_gc, src_gc, ix, jx, kx, diff, sta);
  
  const SKL_REAL ra = 1.0/(SKL_REAL)stepAvr; // output
  
  unsigned long src_idx, dst_idx, src_lsz[3], dst_lsz[3];
  src_lsz[0] = src_sz[0];
  src_lsz[1] = src_sz[1];
  src_lsz[2] = src_sz[2];
  dst_lsz[0] = dst_sz[0];
  dst_lsz[1] = dst_sz[1];
  dst_lsz[2] = dst_sz[2];
  unsigned long s1 = src_lsz[0]*src_lsz[1]*src_lsz[2];
  unsigned long s2 = s1*2;
  register unsigned i, j, k;
  flop += (SKL_REAL)((kx-sta)*(jx-sta)*(ix-sta)*9);
  
  for(k=sta; k<kx; k++){
    unsigned kk = k+diff;
    for(j=sta; j<jx; j++){
      unsigned jj = j+diff;
      for(i=sta; i<ix; i++){
        // src_idx : indecies for SklVector3D
        src_idx = src_lsz[0]*src_lsz[1]*kk + src_lsz[0]*jj + (i+diff);
        // dst_idx : indecies for SklVector3DEx
        dst_idx = 3*(dst_lsz[0]*dst_lsz[1]*k + dst_lsz[0]*j + i);
        dst_data[dst_idx  ] = (src_data[src_idx   ]*ra - v00[0]) * Ref_v;
        dst_data[dst_idx+1] = (src_data[src_idx+s1]*ra - v00[1]) * Ref_v;
        dst_data[dst_idx+2] = (src_data[src_idx+s2]*ra - v00[2]) * Ref_v;
      }
    }
  }
}

/**
 @fn void FileIO::cnv_TP_ND2D(SklScalar3D<SKL_REAL>* dst, const SklScalar3D<SKL_REAL>* src, const SKL_REAL Ref_rho, const SKL_REAL Ref_v, SKL_REAL& flop) const
 @brief 全圧データについて，無次元から有次元単位に変換する
 @param dst 単位変換のデータクラス
 @param src 単位変換前のデータクラス
 @param Ref_rho 代表密度(kg/m^3)
 @param Ref_v 代表速度(m/s)
 @param flop 浮動小数演算数
 @see SklUtil::cpyS3D()
 */
void FileIO::cnv_TP_ND2D(SklScalar3D<SKL_REAL>* dst, const SklScalar3D<SKL_REAL>* src, const SKL_REAL Ref_rho, const SKL_REAL Ref_v, SKL_REAL& flop) const
{
  if( !dst || !src ) assert(0);
  
  const unsigned* dst_sz = dst->GetSize();
  const unsigned* src_sz = src->GetSize();
  if(   (src_sz[0] != dst_sz[0])
     || (src_sz[1] != dst_sz[1])
     || (src_sz[2] != dst_sz[2]) ) assert(0);
  
  SKL_REAL* dst_data = dst->GetData();
  const SKL_REAL* src_data = src->GetData();
  if( !dst_data || !src_data ) assert(0);
  
  dst_sz = dst->_GetSize();
  src_sz = src->_GetSize();
  unsigned dst_gc = dst->GetVCellSize();
  unsigned src_gc = src->GetVCellSize();
  unsigned sta, ix, jx, kx;
  int diff;
  register unsigned i, j, k;
  
  SKL_REAL c = Ref_rho*Ref_v*Ref_v;
  flop += (SKL_REAL)(dst_sz[0]*dst_sz[1]*dst_sz[2]*1);
  
  unsigned long src_idx, dst_idx, idx, src_lsz[3], dst_lsz[3];
  src_lsz[0] = src_sz[0];
  src_lsz[1] = src_sz[1];
  src_lsz[2] = src_sz[2];
  dst_lsz[0] = dst_sz[0];
  dst_lsz[1] = dst_sz[1];
  dst_lsz[2] = dst_sz[2];
  
  // compare whole size
  if(   (src_sz[0] == dst_sz[0])
     && (src_sz[1] == dst_sz[1])
     && (src_sz[2] == dst_sz[2]) ) {
    unsigned idx;
    
    for(k=0; k<dst_sz[2]; k++){
      for(j=0; j<dst_sz[1]; j++){
        for(i=0; i<dst_sz[0]; i++){
          idx = dst_lsz[0]*dst_lsz[1]*k  + dst_lsz[0]*j  + i;
          dst_data[idx] = src_data[idx]*c;
        }
      }
    }
  }
  else {
    Core_Utility::CalcIndex(dst_sz[0], dst_sz[1], dst_sz[2], dst_gc, src_gc, ix, jx, kx, diff, sta);
    
    for(k=sta; k<kx; k++){
      unsigned kk = k+diff;
      for(j=sta; j<jx; j++){
        unsigned jj = j+diff;
        for(i=sta; i<ix; i++){
          src_idx = src_lsz[0]*src_lsz[1]*kk + src_lsz[0]*jj + i+diff;
          dst_idx = dst_lsz[0]*dst_lsz[1]*k  + dst_lsz[0]*j  + i;
          dst_data[dst_idx] = src_data[src_idx]*c;
        }
      }
    }
  }
}

/**
 @fn void FileIO::cnv_V_D2ND(SklVector3D<SKL_REAL>* dst, const SKL_REAL Ref_v) const
 @brief 速度データについて，有次元単位から無次元に変換する
 @param dst 単位変換のデータクラス
 @param Ref_v 代表速度(m/s)
 */
void FileIO::cnv_V_D2ND(SklVector3D<SKL_REAL>* dst, const SKL_REAL Ref_v) const
{
  if( !dst ) assert(0);
  
  //const unsigned* dst_sz = dst->GetSize();
  SKL_REAL* dst_data = dst->GetData();
  if( !dst_data ) assert(0);
  
  const unsigned* dst_sz = dst->_GetSize(); // with guide cell -> ix+vc*2
  register unsigned i, j, k, l;
  SKL_REAL c = 1.0/Ref_v;
  
  unsigned long idx, dst_lsz[3];
  dst_lsz[0] = dst_sz[0];
  dst_lsz[1] = dst_sz[1];
  dst_lsz[2] = dst_sz[2];
    
  for (l=0; l<3; l++) {
    for(k=0; k<dst_sz[2]; k++){
      for(j=0; j<dst_sz[1]; j++){
        for(i=0; i<dst_sz[0]; i++){
          idx = dst_lsz[0]*dst_lsz[1]*dst_lsz[2]*l  + dst_lsz[0]*dst_lsz[1]*k  + dst_lsz[0]*j + i;
          dst_data[idx] *= c;
        }
      }
    }
  }
}

/**
 @fn void FileIO::cnv_T_D2ND(SklScalar3D<SKL_REAL>* dst, const SKL_REAL Base_tmp, const SKL_REAL Diff_tmp, const unsigned Unit) const
 @brief 温度データについて，Unitで指定される有次元単位から無次元に変換する
 @param dst 単位変換のデータクラス
 @param Base_tmp 基準温度(K or C)
 @param Diff_tmp 代表温度差(K or C)
 @param Unit 温度単位(K or C)
 */
void FileIO::cnv_T_D2ND(SklScalar3D<SKL_REAL>* dst, const SKL_REAL Base_tmp, const SKL_REAL Diff_tmp, const unsigned Unit) const
{
  if( !dst ) assert(0);
  
  //const unsigned* dst_sz = dst->GetSize();
  SKL_REAL* dst_data = dst->GetData();
  if( !dst_data ) assert(0);
  
  const unsigned* dst_sz = dst->_GetSize(); // with guide cell -> ix+vc*2
  register unsigned i, j, k;
  
  unsigned long idx, dst_lsz[3];
  dst_lsz[0] = dst_sz[0];
  dst_lsz[1] = dst_sz[1];
  dst_lsz[2] = dst_sz[2];
    
  for(k=0; k<dst_sz[2]; k++){
    for(j=0; j<dst_sz[1]; j++){
      for(i=0; i<dst_sz[0]; i++){
        idx = dst_lsz[0]*dst_lsz[1]*k  + dst_lsz[0]*j  + i;
        dst_data[idx] = FBUtility::convD2ND(dst_data[idx], Base_tmp, Diff_tmp, Unit);
      }
    }
  }
}

/**
 @fn void FileIO::cnv_T_ND2D(SklScalar3D<SKL_REAL>* dst, const SklScalar3D<SKL_REAL>* src, const SKL_REAL Base_tmp, constSKL_REAL Diff_tmp, 
                             const unsigned Unit, SKL_REAL& flop) const
 @brief ファイル出力時，温度データについて，無次元からUnitで指定される有次元単位に変換する
 @param dst 単位変換のデータクラス
 @param src 単位変換前のデータクラス
 @param Base_tmp 基準温度(K or C)
 @param Diff_tmp 代表温度差(K or C)
 @param Unit 温度単位(K or C)
 @param flop 浮動小数演算数
 @see SklUtil::cpyS3D()
 */
void FileIO::cnv_T_ND2D(SklScalar3D<SKL_REAL>* dst, const SklScalar3D<SKL_REAL>* src, const SKL_REAL Base_tmp, const SKL_REAL Diff_tmp, 
                        const unsigned Unit, SKL_REAL& flop) const
{
  if( !dst || !src ) assert(0);
  
  const unsigned* dst_sz = dst->GetSize();
  const unsigned* src_sz = src->GetSize();
  if(   (src_sz[0] != dst_sz[0])
     || (src_sz[1] != dst_sz[1])
     || (src_sz[2] != dst_sz[2]) ) assert(0);
  
  SKL_REAL* dst_data = dst->GetData();
  const SKL_REAL* src_data = src->GetData();
  if( !dst_data || !src_data ) assert(0);
  
  dst_sz = dst->_GetSize();
  src_sz = src->_GetSize();
  unsigned dst_gc = dst->GetVCellSize();
  unsigned src_gc = src->GetVCellSize();
  unsigned sta, ix, jx, kx;
  int diff;
  register unsigned i, j, k;

  SKL_REAL tc;
  flop += (SKL_REAL)(dst_sz[0]*dst_sz[1]*dst_sz[2]*2);
  
  unsigned long src_idx, dst_idx, idx, src_lsz[3], dst_lsz[3];
  src_lsz[0] = src_sz[0];
  src_lsz[1] = src_sz[1];
  src_lsz[2] = src_sz[2];
  dst_lsz[0] = dst_sz[0];
  dst_lsz[1] = dst_sz[1];
  dst_lsz[2] = dst_sz[2];
  
  // compare whole size
  if(   (src_sz[0] == dst_sz[0])
     && (src_sz[1] == dst_sz[1])
     && (src_sz[2] == dst_sz[2]) ) {
    unsigned idx;
    
    for(k=0; k<dst_sz[2]; k++){
      for(j=0; j<dst_sz[1]; j++){
        for(i=0; i<dst_sz[0]; i++){
          idx = dst_lsz[0]*dst_lsz[1]*k  + dst_lsz[0]*j  + i;
          tc = FBUtility::convND2Kelvin(src_data[idx], Base_tmp, Diff_tmp);
          dst_data[idx] = FBUtility::convK2Temp(tc, Unit);
        }
      }
    }
  }
  else {
    Core_Utility::CalcIndex(dst_sz[0], dst_sz[1], dst_sz[2], dst_gc, src_gc, ix, jx, kx, diff, sta);
    
    for(k=sta; k<kx; k++){
      unsigned kk = k+diff;
      for(j=sta; j<jx; j++){
        unsigned jj = j+diff;
        for(i=sta; i<ix; i++){
          src_idx = src_lsz[0]*src_lsz[1]*kk + src_lsz[0]*jj + i+diff;
          dst_idx = dst_lsz[0]*dst_lsz[1]*k  + dst_lsz[0]*j  + i;
          tc = FBUtility::convND2Kelvin(src_data[src_idx], Base_tmp, Diff_tmp);
          dst_data[dst_idx] = FBUtility::convK2Temp(tc, Unit);
        }
      }
    }
  } 
}

/**
 @fn void FileIO::cnv_P_D2ND(SklScalar3D<SKL_REAL>* dst, const SKL_REAL Base_prs, const SKL_REAL Ref_rho, const SKL_REAL Ref_v, const unsigned mode) const
 @brief 圧力データについて，有次元から無次元に変換する
 @param dst 単位変換のデータクラス
 @param Base_prs 基準圧力(Pa)
 @param Ref_rho 代表密度(kg/m^3)
 @param Ref_v 代表速度(m/s)
 @param mode unit of pressure
 */
void FileIO::cnv_P_D2ND(SklScalar3D<SKL_REAL>* dst, const SKL_REAL Base_prs, const SKL_REAL Ref_rho, const SKL_REAL Ref_v, const unsigned mode) const
{
  if( !dst ) assert(0);
  
  //const unsigned* dst_sz = dst->GetSize();
  SKL_REAL* dst_data = dst->GetData();
  if( !dst_data ) assert(0);
  
  const unsigned* dst_sz = dst->_GetSize(); // with guide cell -> ix+vc*2
  register unsigned i, j, k;
  
  unsigned long idx, dst_lsz[3];
  dst_lsz[0] = dst_sz[0];
  dst_lsz[1] = dst_sz[1];
  dst_lsz[2] = dst_sz[2];
  
  for(k=0; k<dst_sz[2]; k++){
    for(j=0; j<dst_sz[1]; j++){
      for(i=0; i<dst_sz[0]; i++){
        idx = dst_lsz[0]*dst_lsz[1]*k  + dst_lsz[0]*j  + i;
        dst_data[idx] = FBUtility::convD2ND_P(dst_data[idx], Base_prs, Ref_rho, Ref_v, mode);
      }
    }
  }
}

/**
 @fn void FileIO::cnv_P_ND2D(SklScalar3D<SKL_REAL>* dst, const SklScalar3D<SKL_REAL>* src, const SKL_REAL Base_prs, const SKL_REAL Ref_rho, 
                            const SKL_REAL Ref_v, const unsigned mode, SKL_REAL& flop) const
 @brief ファイル出力時，圧力データについて無次元から有次元単位に変換する
 @param dst 単位変換のデータクラス
 @param src 単位変換前のデータクラス
 @param Base_prs 基準圧力(Pa)
 @param Ref_rho 代表密度(kg/m^3)
 @param Ref_v 代表速度(m/s)
 @param mode unit of pressure
 @param flop 浮動小数演算数
 @see SklUtil::cpyS3D()
 */
void FileIO::cnv_P_ND2D(SklScalar3D<SKL_REAL>* dst, const SklScalar3D<SKL_REAL>* src, const SKL_REAL Base_prs, const SKL_REAL Ref_rho, 
                        const SKL_REAL Ref_v, const unsigned mode, SKL_REAL& flop) const
{
  if( !dst || !src ) assert(0);
  
  const unsigned* dst_sz = dst->GetSize();
  const unsigned* src_sz = src->GetSize();
  if(   (src_sz[0] != dst_sz[0])
     || (src_sz[1] != dst_sz[1])
     || (src_sz[2] != dst_sz[2]) ) assert(0);
  
  SKL_REAL* dst_data = dst->GetData();
  const SKL_REAL* src_data = src->GetData();
  if( !dst_data || !src_data ) assert(0);
  
  dst_sz = dst->_GetSize();
  src_sz = src->_GetSize();
  unsigned dst_gc = dst->GetVCellSize();
  unsigned src_gc = src->GetVCellSize();
  unsigned sta, ix, jx, kx;
  int diff;
  register unsigned i, j, k;
  
  unsigned long idx, dst_lsz[3];
  dst_lsz[0] = dst_sz[0];
  dst_lsz[1] = dst_sz[1];
  dst_lsz[2] = dst_sz[2];
  
  // compare whole size
  if(   (src_sz[0] == dst_sz[0])
     && (src_sz[1] == dst_sz[1])
     && (src_sz[2] == dst_sz[2]) ) {
    
    for(k=0; k<dst_sz[2]; k++){
      for(j=0; j<dst_sz[1]; j++){
        for(i=0; i<dst_sz[0]; i++){
          idx = dst_lsz[0]*dst_lsz[1]*k  + dst_lsz[0]*j  + i;
          dst_data[idx] = FBUtility::convND2D_P(src_data[idx], Base_prs, Ref_rho, Ref_v, mode);
        }
      }
    }
    flop += (SKL_REAL)(dst_sz[0]*dst_sz[1]*dst_sz[2]*4);
  }
  else {
    assert(0);
  }
}

/**
 @fn void FileIO::loadSphScalar4D(SklSolverBase* obj, FILE* fp, const char* fname, const unsigned* size, const unsigned guide,
                                  SklScalar4D<SKL_REAL>* dc_s, int& step, SKL_REAL& time, unsigned Dmode) const
 @brief sphスカラー4Dファイルをロードする
 @param solver_obj SklSolverBaseクラスのオブジェクトのポインタ (this)
 @param fp ファイルポインタ（ファイル出力）
 @param fname InFileのattrラベル名
 @param size グローバルなサイズ
 @param guide ガイドセルサイズ
 @param dc_s  結果を保持するデータクラス
 @param step ステップ
 @param time 時刻
 @param Dmode 次元（無次元-0 / 有次元-1）
 */
void FileIO::loadSphScalar4D(SklSolverBase* obj, FILE* fp, const char* fname, const unsigned* size, const unsigned guide, 
                             SklScalar4D<SKL_REAL>* dc_s, int& step, SKL_REAL& time, unsigned Dmode) const
{
  if ( !obj )   assert(0);
  // if ( !fp )    assert(0); 並列時にはfpをオープンしない．ホストのみ
  if ( !dc_s )  assert(0);
  if ( !fname ) assert(0);
  
  SklParaComponent* para_cmp = ParaCmpo;
  SklParaManager* para_mng = para_cmp->GetParaManager();
  SklVoxDataSet* sphS = NULL;

  if( !(sphS=obj->LoadFile(para_mng, fname, guide, pn.procGrp, size)) ) {
    Hostonly_ stamped_printf     ("\tError: during '%s' reading\n", sphS->GetFileName());
    Hostonly_ stamped_fprintf(fp, "\tError: during '%s' reading\n", sphS->GetFileName());
    assert(0);
  }

  if ( !SklUtil::getTimeStamp(sphS, step, time) ) {
    Hostonly_ stamped_printf     ("\tError: during getTimeStep for %s\n", sphS->GetFileName());
    Hostonly_ stamped_fprintf(fp, "\tError: during getTimeStep for %s\n", sphS->GetFileName());
    assert(0);
  }
  
  SklScalar4D<SKL_REAL>* dc_tmp;
  if( !(dc_tmp = dynamic_cast<SklScalar4D<SKL_REAL>*>(sphS->GetData(SklVoxDataSet::SPH_DATA))) ) assert(0);
  //if( !SklUtil::cpyS4D(dc_s, dc_tmp) ) assert(0);
  
  Hostonly_ printf     ("\t[%s] has read :\tstep=%d  time=%e [%s]\n", sphS->GetFileName(), step, time, (Dmode==DIMENSIONAL)?"sec.":"-");
  Hostonly_ fprintf(fp, "\t[%s] has read :\tstep=%d  time=%e [%s]\n", sphS->GetFileName(), step, time, (Dmode==DIMENSIONAL)?"sec.":"-");
  
  if (sphS) { delete sphS; sphS=NULL; }
}

/**
 @fn void FileIO::loadSphVector3DAvr(SklSolverBase* obj, FILE* fp, const char* fname, const unsigned* size, const unsigned guide, 
                                     SklVector3D<SKL_REAL>* dc_v, int& step, SKL_REAL& time, const SKL_REAL* v00, unsigned Dmode) const
 @brief 平均値のsphベクトル3Dファイルをロードする
 @param solver_obj SklSolverBaseクラスのオブジェクトのポインタ (this)
 @param fp ファイルポインタ（ファイル出力）
 @param fname InFileのattrラベル名
 @param size グローバルなサイズ
 @param guide ガイドセルサイズ
 @param dc_v  結果を保持するデータクラス
 @param step ステップ
 @param time 時刻
 @param v00 格子速度
 @param Dmode 次元（無次元-0 / 有次元-1）
 */
void FileIO::loadSphVector3DAvr(SklSolverBase* obj, FILE* fp, const char* fname, const unsigned* size, const unsigned guide, 
                                SklVector3D<SKL_REAL>* dc_v, int& step, SKL_REAL& time, SKL_REAL* v00, unsigned Dmode) const
{
  if ( !obj )   assert(0);
  // if ( !fp )    assert(0); 並列時にはfpをオープンしない．ホストのみ
  if ( !dc_v )  assert(0);
  if ( !fname ) assert(0);
  
  SklParaComponent* para_cmp = ParaCmpo;
  SklParaManager* para_mng = para_cmp->GetParaManager();
  SklVoxDataSet* sphV = NULL;
  
  if( !(sphV=obj->LoadFile(para_mng, fname, guide, pn.procGrp, size)) ) {
    Hostonly_ stamped_printf     ("\tError: during '%s' reading\n", sphV->GetFileName());
    Hostonly_ stamped_fprintf(fp, "\tError: during '%s' reading\n", sphV->GetFileName());
    assert(0);
  }

  if ( !SklUtil::getTimeStamp(sphV, step, time) ) {
    Hostonly_ stamped_printf     ("\tError: during getTimeStep for %s\n", sphV->GetFileName());
    Hostonly_ stamped_fprintf(fp, "\tError: during getTimeStep for %s\n", sphV->GetFileName());
    assert(0);
  }
  
  SklVector3DEx<SKL_REAL>* dc_tmp;
  if( !(dc_tmp = dynamic_cast<SklVector3DEx<SKL_REAL>*>(sphV->GetData(SklVoxDataSet::SPH_DATA))) ) assert(0);
  if( !Core_Utility::shiftVin3D(dc_v, dc_tmp, &v00[1], (unsigned)step) ) assert(0);
  //if( !(DataMngr.RegistData("avrv", dc_av)) ) assert(0);
  
  Hostonly_ printf     ("\t[%s] has read :\tstep=%d  time=%e [%s]\n", sphV->GetFileName(), step, time, (Dmode==DIMENSIONAL)?"sec.":"-");
  Hostonly_ fprintf(fp, "\t[%s] has read :\tstep=%d  time=%e [%s]\n", sphV->GetFileName(), step, time, (Dmode==DIMENSIONAL)?"sec.":"-");
  
  if (sphV) { delete sphV; sphV=NULL; }
}

/**
 @fn void FileIO::loadSphScalar3DAvr(SklSolverBase* obj, FILE* fp, const char* fname, const unsigned* size, const unsigned guide, 
                                     SklScalar3D<SKL_REAL>* dc_s, int& step, SKL_REAL& time, unsigned Dmode) const
 @brief 平均値のsphスカラー3Dファイルをロードする
 @param solver_obj SklSolverBaseクラスのオブジェクトのポインタ (this)
 @param fp ファイルポインタ（ファイル出力）
 @param fname InFileのattrラベル名
 @param size グローバルなサイズ
 @param guide ガイドセルサイズ
 @param dc_s  結果を保持するデータクラス
 @param step ステップ
 @param time 時刻
 @param Dmode 次元（無次元-0 / 有次元-1）
 */
void FileIO::loadSphScalar3DAvr(SklSolverBase* obj, FILE* fp, const char* fname, const unsigned* size, const unsigned guide, 
                                SklScalar3D<SKL_REAL>* dc_s, int& step, SKL_REAL& time, unsigned Dmode) const
{
  if ( !obj )   assert(0);
  // if ( !fp )    assert(0); 並列時にはfpをオープンしない．ホストのみ
  if ( !dc_s )  assert(0);
  if ( !fname ) assert(0);
  
  SklParaComponent* para_cmp = ParaCmpo;
  SklParaManager* para_mng = para_cmp->GetParaManager();
  SklVoxDataSet* sphS = NULL;

  if( !(sphS=obj->LoadFile(para_mng, fname, guide, pn.procGrp, size)) ) {
    Hostonly_ stamped_printf     ("\tError: during '%s' reading\n", sphS->GetFileName());
    Hostonly_ stamped_fprintf(fp, "\tError: during '%s' reading\n", sphS->GetFileName());
    assert(0);
  }

  if ( !SklUtil::getTimeStamp(sphS, step, time) ) {
    Hostonly_ stamped_printf     ("\tError: during getTimeStep for %s\n", sphS->GetFileName());
    Hostonly_ stamped_fprintf(fp, "\tError: during getTimeStep for %s\n", sphS->GetFileName());
    assert(0);
  }
  
  SklScalar3D<SKL_REAL>* dc_tmp;
  if( !(dc_tmp = dynamic_cast<SklScalar3D<SKL_REAL>*>(sphS->GetData(SklVoxDataSet::SPH_DATA))) ) assert(0);
  if( !SklUtil::scaleAvrS3D(dc_tmp, (unsigned)step, true) ) assert(0); // true => INPUT
  if( !SklUtil::cpyS3D(dc_s, dc_tmp) ) assert(0);
  //if( !(DataMngr.RegistData("avrp", dc_ap)) ) assert(0);
  
  Hostonly_ printf     ("\t[%s] has read :\tstep=%d  time=%e [%s]\n", sphS->GetFileName(), step, time, (Dmode==DIMENSIONAL)?"sec.":"-");
  Hostonly_ fprintf(fp, "\t[%s] has read :\tstep=%d  time=%e [%s]\n", sphS->GetFileName(), step, time, (Dmode==DIMENSIONAL)?"sec.":"-");
  
  if (sphS) { delete sphS; sphS=NULL; }
}

/**
 @fn void FileIO::loadSphVector3D(SklSolverBase* obj, FILE* fp, const char* fname, const unsigned* size, const unsigned guide, 
                                  SklVector3D<SKL_REAL>* dc_v, int& step, SKL_REAL& time, SKL_REAL *v00, unsigned Dmode) const
 @brief sphベクトル3Dファイルをロードする
 @param solver_obj SklSolverBaseクラスのオブジェクトのポインタ (this)
 @param fp ファイルポインタ（ファイル出力）
 @param fname InFileのattrラベル名
 @param size グローバルなサイズ
 @param guide ガイドセルサイズ
 @param dc_v  結果を保持するデータクラス
 @param step ステップ
 @param time 時刻
 @param v00[4]
 @param Dmode 次元（無次元-0 / 有次元-1）
 */
void FileIO::loadSphVector3D(SklSolverBase* obj, FILE* fp, const char* fname, const unsigned* size, const unsigned guide, 
                             SklVector3D<SKL_REAL>* dc_v, int& step, SKL_REAL& time, SKL_REAL *v00, unsigned Dmode) const
{
  if ( !obj )   assert(0);
  // if ( !fp )    assert(0); 並列時にはfpをオープンしない．ホストのみ
  if ( !dc_v )  assert(0);
  if ( !fname ) assert(0);
  
  SklParaComponent* para_cmp = ParaCmpo;
  SklParaManager* para_mng = para_cmp->GetParaManager();
  SklVoxDataSet* sphV = NULL;
  
  if( !(sphV=obj->LoadFile(para_mng, fname, guide, pn.procGrp, size)) ) {
    Hostonly_ stamped_printf     ("\tError: during '%s' reading\n", sphV->GetFileName());
    Hostonly_ stamped_fprintf(fp, "\tError: during '%s' reading\n", sphV->GetFileName());
    assert(0);
  }

  if ( !SklUtil::getTimeStamp(sphV, step, time) ) {
    Hostonly_ stamped_printf     ("\tError: during getTimeStep for %s\n", sphV->GetFileName());
    Hostonly_ stamped_fprintf(fp, "\tError: during getTimeStep for %s\n", sphV->GetFileName());
    assert(0);
  }
  
  SklVector3DEx<SKL_REAL>* dc_tmp;
  if( !(dc_tmp = dynamic_cast<SklVector3DEx<SKL_REAL>*>(sphV->GetData(SklVoxDataSet::SPH_DATA))) ) assert(0);
  if( !Core_Utility::shiftVin3D(dc_v, dc_tmp, &v00[1]) ) assert(0);
  
  Hostonly_ printf     ("\t[%s] has read :\tstep=%d  time=%e [%s]\n", sphV->GetFileName(), step, time, (Dmode==DIMENSIONAL)?"sec.":"-");
  Hostonly_ fprintf(fp, "\t[%s] has read :\tstep=%d  time=%e [%s]\n", sphV->GetFileName(), step, time, (Dmode==DIMENSIONAL)?"sec.":"-");
  
  if (sphV) { delete sphV; sphV=NULL; }
}

/**
 @fn void FileIO::loadSphScalar3D(SklSolverBase* obj, FILE* fp, const char* fname, const unsigned* size, const unsigned guide, 
                                  SklScalar3D<SKL_REAL>* dc_s, int& step, SKL_REAL& time, unsigned Dmode) const
 @brief sphスカラー3Dファイルをロードする
 @param solver_obj SklSolverBaseクラスのオブジェクトのポインタ (this)
 @param fp ファイルポインタ（ファイル出力）
 @param fname InFileのattrラベル名
 @param size グローバルなサイズ
 @param guide ガイドセルサイズ
 @param dc_s  結果を保持するデータクラス
 @param step[out] ステップ
 @param time[out] 時刻
 @param Dmode 次元（無次元-0 / 有次元-1）
 */
void FileIO::loadSphScalar3D(SklSolverBase* obj, FILE* fp, const char* fname, const unsigned* size, const unsigned guide, 
                             SklScalar3D<SKL_REAL>* dc_s, int& step, SKL_REAL& time, unsigned Dmode) const
{
  if ( !obj )   assert(0);
  // if ( !fp )    assert(0); 並列時にはfpをオープンしない．ホストのみ
  if ( !dc_s )  assert(0);
  if ( !fname ) assert(0);
  
  SklParaComponent* para_cmp = ParaCmpo;
  SklParaManager* para_mng = para_cmp->GetParaManager();
  SklVoxDataSet* sphS = NULL;

  if( !(sphS=obj->LoadFile(para_mng, fname, guide, pn.procGrp, size)) ) {
    Hostonly_ stamped_printf     ("\tError: during '%s' reading\n", sphS->GetFileName());
    Hostonly_ stamped_fprintf(fp, "\tError: during '%s' reading\n", sphS->GetFileName());
    assert(0);
  }
  
  if ( !SklUtil::getTimeStamp(sphS, step, time) ) {
    Hostonly_ stamped_printf     ("\tError: during getTimeStep for %s\n", sphS->GetFileName());
    Hostonly_ stamped_fprintf(fp, "\tError: during getTimeStep for %s\n", sphS->GetFileName());
    assert(0);
  }

  SklScalar3D<SKL_REAL>* dc_tmp;
  if( !(dc_tmp = dynamic_cast<SklScalar3D<SKL_REAL>*>(sphS->GetData(SklVoxDataSet::SPH_DATA))) ) assert(0);
  if( !SklUtil::cpyS3D(dc_s, dc_tmp) ) assert(0);

  Hostonly_ printf     ("\t[%s] has read :\tstep=%d  time=%e [%s]\n", sphS->GetFileName(), step, time, (Dmode==DIMENSIONAL)?"sec.":"-");
  Hostonly_ fprintf(fp, "\t[%s] has read :\tstep=%d  time=%e [%s]\n", sphS->GetFileName(), step, time, (Dmode==DIMENSIONAL)?"sec.":"-");
  
  if (sphS) { delete sphS; sphS=NULL; }
}

/**
 @fn void FileIO::readSVX(SklSolverBase* obj, FILE* fp, const char* fname, const unsigned* size, const unsigned guide,
                          SklScalar3D<int>* dc_mid, const bool vf_mode, SklScalar3D<SKL_REAL>* dc_ws) const
 @brief svxファイルをロードする
 @param solver_obj SklSolverBaseクラスのオブジェクトのポインタ (this)
 @param fp ファイルポインタ（ファイル出力）
 @param fname InFileのattrラベル名
 @param size グローバルな分割数
 @param guide ガイドセルサイズ
 @param dc_mid IDを保持するデータクラス
 @param vf_mode Volume Fractionを読み込むモード(trueのとき，デフォルト false）
 @param dc_ws  体積率を保持するデータクラス
 @todo 
    - cpyS3DCenter(dc_ws, dc_tmp_svx) が倍精度のときの挙動
 */
void FileIO::readSVX(SklSolverBase* obj, FILE* fp, const char* fname, unsigned* size, unsigned guide, 
                     SklScalar3D<int>* dc_mid, bool vf_mode, SklScalar3D<SKL_REAL>* dc_ws) const
{
  SklParaComponent* para_cmp = ParaCmpo;
  SklParaManager* para_mng = para_cmp->GetParaManager();
  
  if ( !obj )   assert(0);
  if (!dc_mid)  assert(0);
  // if ( !fp )    assert(0); 並列時にはfpをオープンしない．ホストのみ
  if ( !fname ) assert(0);
  
  SklVoxDataSet* svxSet = NULL;

  if ( !(svxSet=obj->LoadFile(para_mng, fname, guide, pn.procGrp, size)) ) { // < Global size
    Hostonly_ {
      printf     ("\tError: during SVX file reading\n");
      fprintf(fp, "\tError: during SVX file reading\n");
    }
    assert(0);
  }

  Hostonly_ {
    printf     ("\n\treading %s\n", svxSet->GetFileName());
    printf     ("\tThis file includes : ");
    fprintf(fp, "\n\treading %s\n", svxSet->GetFileName());
    fprintf(fp, "\tThis file includes : ");
  }
  
  // Medium ID
  if ( !(svxSet->GetData(4)) ) {
    Hostonly_ {
      printf     ("\tThis SVX file does not include Medium ID\n");
      fprintf(fp, "\tThis SVX file does not include Medium ID\n");
    }
    assert(0);
  }
  else {
    Hostonly_ {
      printf     ("\tMedium ID");
      fprintf(fp, "\tMedium ID");
    }
  }
  SklScalar3D<int>* dc_tmp_mid;
  if ( !(dc_tmp_mid = dynamic_cast<SklScalar3D<int>*>(svxSet->GetData(SklVoxDataSet::SVX_MID))) ) assert(0);
  if ( !SklUtil::cpyS3DCenter(dc_mid, dc_tmp_mid) ) assert(0);
  if (dc_tmp_mid) { delete dc_tmp_mid; dc_tmp_mid; }
  
  // Volume Fraction　オプション
  if (vf_mode == true) {
    if (!dc_ws)   assert(0);
    if ( svxSet->GetData(0) ) { // 含まれている場合
      if ( sizeof(SKL_REAL) == sizeof(double) ) {
        Hostonly_ printf("double is not supported at this moment\n");
        if (svxSet) { delete svxSet; svxSet=NULL; }
        return;
      }
      Hostonly_ {
        printf     ("\tVolume Fraction");
        fprintf(fp, "\tVolume Fraction");
      }
      SklScalar3D<SKL_REAL>* dc_tmp_svx;
      if ( !(dc_tmp_svx = dynamic_cast<SklScalar3D<SKL_REAL>*>(svxSet->GetData(SklVoxDataSet::SVX_VOLUME))) ) assert(0);
      if ( !SklUtil::cpyS3DCenter(dc_ws, dc_tmp_svx) ) assert(0);
      if (dc_tmp_svx) { delete dc_tmp_svx; dc_tmp_svx=NULL; }
    }
    else {
      Hostonly_ {
        printf     ("\tThis SVX file does not include Volume Fraction\n");
        fprintf(fp, "\tThis SVX file does not include Volume Fraction\n");
      }
    }
  }
  
  if (svxSet) { delete svxSet; svxSet=NULL; }
  
  Hostonly_ {
    printf     ("\n");
    fprintf(fp, "\n");
  }
}

/**
 @fn void FileIO::readSBX(SklSolverBase* obj, FILE* fp, const char* mid_str, unsigned* size, unsigned guide, 
                          SklScalar3D<int>* dc_mid, SklScalar3D<SKL_REAL>* dc_vol) const
 @brief SBXファイルの読み込み
 @param solver_obj SklSolverBaseクラスのオブジェクトのポインタ (this)
 @param fp ファイルポインタ（ファイル出力）
 @param mid_str コンフィギュレーションファイルにおけるIDの入力情報ラベル
 @param size グローバルな分割数
 @param guide ガイドセルサイズ
 @param dc_mid IDを保持するデータクラス
 @param dc_vol volume rateを保持するデータクラス
 */
void FileIO::readSBX(SklSolverBase* obj, FILE* fp, const char* mid_str, unsigned* size, unsigned guide, 
                     SklScalar3D<int>* dc_mid, SklScalar3D<SKL_REAL>* dc_vol) const
{
  SklParaComponent* para_cmp = ParaCmpo;
  SklParaManager* para_mng = para_cmp->GetParaManager();
  
  if ( !obj )   assert(0);
  if (!dc_mid)  assert(0);
  if (!mid_str) assert(0);
  // if ( !fp )    assert(0); 並列時にはfpをオープンしない．ホストのみ
  
  const unsigned int* mid_sz = dc_mid->_GetSize();
  
  unsigned long mid_lsz[3];
  mid_lsz[0] = mid_sz[0];
  mid_lsz[1] = mid_sz[1];
  mid_lsz[2] = mid_sz[2];
  
  memset(dc_mid->GetData(), 0x00, sizeof(unsigned int)*mid_lsz[0]*mid_lsz[1]*mid_lsz[2]);
  if (dc_vol != NULL) {
    const unsigned int* vol_sz = dc_vol->_GetSize();
    unsigned long vol_lsz[3];
    vol_lsz[0] = vol_sz[0];
    vol_lsz[1] = vol_sz[1];
    vol_lsz[2] = vol_sz[2];
    memset(dc_vol->GetData(), 0x00, sizeof(unsigned int)*vol_lsz[0]*vol_lsz[1]*vol_lsz[2]);
  }

  // load sbx file for MediumID
  loadSBXfile(obj, fp, mid_str, size, guide, dc_mid, dc_vol);
}

/**
 @fn void FileIO::readSBX(SklSolverBase* obj, FILE* fp, const char* file_attr, unsigned* size, unsigned guide, 
                  SklScalar3D<unsigned char>* dc_mid) const
 @brief SBXファイルの読み込み
 @param solver_obj SklSolverBaseクラスのオブジェクトのポインタ (this)
 @param fp ファイルポインタ（ファイル出力）
 @param file_str コンフィギュレーションファイルにおける体積率の入力情報ラベル
 @param size グローバルな分割数
 @param guide ガイドセルサイズ
 @param dc_mid IDを保持するデータクラス
 */
void FileIO::readSBX(SklSolverBase* obj, FILE* fp, const char* file_attr, unsigned* size, unsigned guide, 
                     SklScalar3D<unsigned char>* dc_mid) const
{
  SklParaComponent* para_cmp = ParaCmpo;
  SklParaManager* para_mng = para_cmp->GetParaManager();
  
  if ( !obj )   assert(0);
  if (!dc_mid)  assert(0);
  if (!file_attr) assert(0);
  // if ( !fp )    assert(0); 並列時にはfpをオープンしない．ホストのみ
  
  unsigned char* dataP = dc_mid->GetData();
  const unsigned int* mid_sz = dc_mid->_GetSize();
  unsigned long mid_lsz[3];
  mid_lsz[0] = mid_sz[0];
  mid_lsz[1] = mid_sz[1];
  mid_lsz[2] = mid_sz[2];
  memset(dataP, NULL, sizeof(unsigned int)*mid_lsz[0]*mid_lsz[1]*mid_lsz[2]);
  
  // load sbx file for MediumID
  loadSBXfile(obj, fp, file_attr, size, guide, dc_mid);
}

/**
 @fn void FileIO::loadSBXfile(SklSolverBase* obj, FILE* fp, const char* name, unsigned* size, unsigned guide, 
                              SklScalar3D<int>* mid_data, SklScalar3D<SKL_REAL>* vol_data) const
 @brief SBXデータのロード
 @param obj SklSolverBaseクラスのオブジェクトのポインタ (this)
 @param fp ファイルポインタ（ファイル出力）
 @param name コンフィギュレーションファイルにおける入力情報ラベル
 @param size グローバルなサイズ
 @param guide ガイドセルサイズ
 @param mid_data IDを保持するデータクラス
 @param vol_data volume rateを保持するデータクラス
 */
void FileIO::loadSBXfile(SklSolverBase* obj, FILE* fp, const char* name, unsigned* size, unsigned guide, 
                         SklScalar3D<int>* mid_data, SklScalar3D<SKL_REAL>* vol_data) const
{
  SklParaComponent* para_cmp = ParaCmpo;
  SklParaManager* para_mng = para_cmp->GetParaManager();
  
  if ( !obj )   assert(0);
  if (!name)   assert(0);
  if (!mid_data)   assert(0);
  // if ( !fp )   assert(0); 並列時にはfpをオープンしない．ホストのみ
  
  SklSbxDataSet* sbxSet = NULL;
  
  // read SBX file
  if ( !(sbxSet=dynamic_cast<SklSbxDataSet*>(obj->LoadFile(para_mng, name, guide, pn.procGrp))) ) {
    Hostonly_ printf     ("\tError: during SBX reading : %s\n", name);
    Hostonly_ fprintf(fp, "\tError: during SBX reading : %s\n", name);
    assert(0);
  }
  const char* sbx_fname = sbxSet->GetFileName();
  if (sbx_fname) {
    Hostonly_ printf     ("\n\treading %s\n", sbxSet->GetFileName());
    Hostonly_ fprintf(fp, "\n\treading %s\n", sbxSet->GetFileName());
  }
  Hostonly_ printf     ("\n\tThis file includes : ");
  Hostonly_ fprintf(fp, "\n\tThis file includes : ");
  
  // check size/org/pitch between in XML and sbx
  unsigned ix=0, jx=0, kx=0;
  if ( !(sbxSet->GetSize(ix, jx, kx)) ) {
    Hostonly_ printf     ("\tSize in SBX file is out of order : %s\n", name);
    Hostonly_ fprintf(fp, "\tSize in SBX file is out of order : %s\n", name);
    assert(0);
  }
  
  // Get SBX Data
  SklScalar3D<int>* dc_tmp_sbx = NULL;
  SklVoxDataSet::VszType vsz = sbxSet->GetVszType();
  const unsigned int sz[3] = {ix, jx, kx};
  if( !(dc_tmp_sbx = dynamic_cast<SklScalar3D<int>*>(
                                  obj->SklAllocateArray(
                                    para_mng,
                                    "sbx",
                                    SKL_CLASS_SCALAR3D,
                                    SKL_ARRAY_DTYPE_INT,
                                    sz,
                                    guide,
                                    pn.procGrp))) ) {
                                      Hostonly_ printf     ("Error : SBX AllocateArray(sbx) fault : %s\n", name);
                                      Hostonly_ fprintf(fp, "Error : SBX AllocateArray(sbx) fault : %s\n", name);
                                      assert(0);
                                    }
  
  if (vsz == SklVoxDataSet::TYPE_CHAR) {
    SklScalar3D<unsigned char>* dc_tmp_char;
    if ( !(dc_tmp_char = dynamic_cast<SklScalar3D<unsigned char>*>(sbxSet->GetData(SklVoxDataSet::SBX_DATA))) ) {
      Hostonly_ printf     ("Error : SBX Data GetData fault : %s\n", name);
      Hostonly_ fprintf(fp, "Error : SBX Data GetData fault : %s\n", name);
      assert(0);
    }
    const unsigned char* src_ds = dc_tmp_char->GetData();
    int* dest_ds = dc_tmp_sbx->GetData();
    
    // AUX data
    int lowerBits = sbxSet->GetLowerBits();
    int upperBits = sbxSet->GetUpperBits();
    int bitsFlag = sbxSet->GetValidBitsFlag();
    if (lowerBits < 0 || lowerBits > 8) {
      Hostonly_ printf     ("Error : SBX LowerBits[=%d] fault : %s\n", lowerBits, name);
      Hostonly_ fprintf(fp, "Error : SBX LowerBits[=%d] fault : %s\n", lowerBits, name);
      assert(0);
    }
    if (upperBits < 0 || upperBits > 8) {
      Hostonly_ printf     ("Error : SBX UpperBits[=%d] fault : %s\n", upperBits, name);
      Hostonly_ fprintf(fp, "Error : SBX UpperBits[=%d] fault : %s\n", upperBits, name);
      assert(0);
    }
    
    Hostonly_ printf     ("\tSBX Data type : 1 byte ");
    Hostonly_ fprintf(fp, "\tSBX Data type : 1 byte ");

    switch (bitsFlag) {
      case 2:
        Hostonly_ printf     ("-> ID only\n");
        Hostonly_ fprintf(fp, "->ID only\n");
        break;
        
      case 0:
        Hostonly_ printf     ("-> 6bits(ID) + 2bits(VF)\n");
        Hostonly_ fprintf(fp, "-> 6bits(ID) + 2bits(VF)\n");
        break;
        
      case 1:
        Hostonly_ printf     ("-> VF only\n");
        Hostonly_ fprintf(fp, "-> VF only\n");
        break;
        
      default:
        Hostonly_ printf     ("-> bit encoding error\n");
        Hostonly_ fprintf(fp, "-> bit encoding error\n");
        assert(0);
    }
    
    for (register unsigned i=0; i<dc_tmp_char->GetArrayLength(); i++) {
      if (bitsFlag == 0 || bitsFlag == 2) {
        dest_ds[i] = (int)SklUtil::shiftLowerBits(src_ds[i], lowerBits);
      }
    }
    
    if (vol_data != NULL) {
      SklScalar3D<SKL_REAL>* sbx_vol =
      dynamic_cast<SklScalar3D<SKL_REAL>*>(
                                           obj->SklAllocateArray(
                                             para_mng,
                                             "sbx_vol",
                                             SKL_CLASS_SCALAR3D,
                                             SKL_ARRAY_DTYPE_REAL,
                                             sz,
                                             guide,
                                             pn.procGrp) );
      if (!sbx_vol) {
        Hostonly_ printf     ("Error : SBX AllocateArray(sbx_vol) fault : %s\n", name);
        Hostonly_ fprintf(fp, "Error : SBX AllocateArray(sbx_vol) fault : %s\n", name);
        assert(0);
      }
      SKL_REAL* dest_vol = sbx_vol->GetData();
      for (register unsigned i=0; i<dc_tmp_char->GetArrayLength(); i++) {
        if (bitsFlag == 0 || bitsFlag == 1) {
          int bits = SklUtil::maskUpperBits(src_ds[i], upperBits);
          if (bits == 0x00) dest_vol[i] = 0.0;
          else if (bits == 0x01) dest_vol[i] = 0.5;
          else if (bits == 0x02) dest_vol[i] = 1.0;
        }
      }
      if (!SklUtil::cpyS3DOverlap(vol_data, sbx_vol)) {
        Hostonly_ printf     ("Error : SBX Data cpyS3DOverlap[SKL_REAL] fault : %s\n", name);
        Hostonly_ fprintf(fp, "Error : SBX Data cpyS3DOverlap[SKL_REAL] fault : %s\n", name);
        assert(0);
      }
    }
  }
  else if (vsz == SklVoxDataSet::TYPE_SHORT) {
    Hostonly_ stamped_printf     ("\tSBX Data type : 2 bytes\n");
    Hostonly_ stamped_fprintf(fp, "\tSBX Data type : 2 bytes\n");
    SklScalar3D<unsigned short>* dc_tmp_short;
    if ( !(dc_tmp_short = dynamic_cast<SklScalar3D<unsigned short>*>(sbxSet->GetData(SklVoxDataSet::SBX_DATA))) ) {
      Hostonly_ printf     ("Error : SBX Data GetData fault : %s\n", name);
      Hostonly_ fprintf(fp, "Error : SBX Data GetData fault : %s\n", name);
      assert(0);
    }
    const unsigned short* src_ds = dc_tmp_short->GetData();
    int* dest_ds = dc_tmp_sbx->GetData();
    for (register unsigned i=0; i<dc_tmp_short->GetArrayLength(); i++) {
      dest_ds[i] = (int)src_ds[i];
    }
  }
  else if (vsz == SklVoxDataSet::TYPE_INT) {
    Hostonly_ stamped_printf     ("\tSBX Data type : 4 bytes\n");
    Hostonly_ stamped_fprintf(fp, "\tSBX Data type : 4 bytes\n");
    SklScalar3D<unsigned int>* dc_tmp_int;
    if ( !(dc_tmp_int = dynamic_cast<SklScalar3D<unsigned int>*>(sbxSet->GetData(SklVoxDataSet::SBX_DATA))) ) {
      Hostonly_ printf     ("Error : SBX Data GetData fault : %s\n", name);
      Hostonly_ fprintf(fp, "Error : SBX Data GetData fault : %s\n", name);
      assert(0);
    }
    const unsigned int* src_ds = dc_tmp_int->GetData();
    int* dest_ds = dc_tmp_sbx->GetData();
    memcpy(dest_ds, src_ds, dc_tmp_int->GetArrayLength()*sizeof(int));
  }
  
  if (!dc_tmp_sbx) {
    Hostonly_ printf     ("Error : can't set dc_tmp_sbx : %s\n", name);
    Hostonly_ fprintf(fp, "Error : can't set dc_tmp_sbx : %s\n", name);
    assert(0);
  }
  
  // copy data
  if (!SklUtil::cpyS3DOverlap(mid_data, dc_tmp_sbx)) {
    Hostonly_ printf     ("Error : SBX Data cpyS3DOverlap[int] fault : %s\n", name);
    Hostonly_ fprintf(fp, "Error : SBX Data cpyS3DOverlap[int] fault : %s\n", name);
    assert(0);
  }
  
  Hostonly_ printf     ("\tget sbx data : %s\n", name);
  Hostonly_ fprintf(fp, "\tget sbx data : %s\n", name);

  if (dc_tmp_sbx) {
    obj->SklGetDataManager()->DeleteDataObj("sbx");
    dc_tmp_sbx=NULL;
  }
  if (sbxSet)     { delete sbxSet; sbxSet=NULL; }
  
  Hostonly_ printf     ("\t\n");
  Hostonly_ fprintf(fp, "\t\n");
}

/**
 @fn void FileIO::loadSBXfile(SklSolverBase* obj, FILE* fp, const char* name, unsigned* size, unsigned guide, SklScalar3D<unsigned char>* data) const
 @brief SBXデータのロード
 @param solver_obj SklSolverBaseクラスのオブジェクトのポインタ (this)
 @param fp ファイルポインタ（ファイル出力）
 @param name コンフィギュレーションファイルにおける入力情報ラベル
 @param size グローバルなサイズ
 @param guide ガイドセルサイズ
 @param data データクラス
 */
void FileIO::loadSBXfile( SklSolverBase* obj,
                         FILE* fp,
                         const char* name,
                         unsigned* size,
                         unsigned guide,
                         SklScalar3D<unsigned char>* data) const
{
  SklParaComponent* para_cmp = ParaCmpo;
  SklParaManager* para_mng = para_cmp->GetParaManager();
  
  if ( !obj )  assert(0);
  if (!name)   assert(0);
  if (!data)   assert(0);
  // if ( !fp )   assert(0); 並列時にはfpをオープンしない．ホストのみ
  
  SklSbxDataSet* sbxSet = NULL;
  
  // read SBX file
  if ( !(sbxSet=dynamic_cast<SklSbxDataSet*>(obj->LoadFile(para_mng, name, guide, pn.procGrp))) ) {
    Hostonly_ printf     ("\tError: during SBX reading : %s\n", name);
    Hostonly_ fprintf(fp, "\tError: during SBX reading : %s\n", name);
    assert(0);
  }
  
  const char* sbx_fname = sbxSet->GetFileName();
  if (sbx_fname) {
    Hostonly_ printf     ("\n\treading %s\n", sbxSet->GetFileName());
    Hostonly_ fprintf(fp, "\n\treading %s\n", sbxSet->GetFileName());
  }
  
  Hostonly_ printf     ("\n\tThis file includes : ");
  Hostonly_ fprintf(fp, "\n\tThis file includes : ");
  
  // Get SBX Data
  SklVoxDataSet::VszType vsz = sbxSet->GetVszType();
  
  register int i, j, k;
  if (vsz == SklVoxDataSet::TYPE_CHAR) {
    SklScalar3D<unsigned char>* dc_src_char;
    // get data
    if ( !(dc_src_char = dynamic_cast<SklScalar3D<unsigned char>*>(sbxSet->GetData(SklVoxDataSet::SBX_DATA))) ) {
      Hostonly_ printf     ("Error : SBX Data GetData fault : %s\n", name);
      Hostonly_ fprintf(fp, "Error : SBX Data GetData fault : %s\n", name);
      assert(0);
    }
    
    // copy data
    if (!SklUtil::cpyS3DOverlap(data, dc_src_char)) {
      Hostonly_ printf     ("Error : SBX Data cpyS3DOverlap fault : %s\n", name);
      Hostonly_ fprintf(fp, "Error : SBX Data cpyS3DOverlap fault : %s\n", name);
      assert(0);
    }
  }
  else if (vsz == SklVoxDataSet::TYPE_SHORT) {
    Hostonly_ printf     ("Error : SBX Data Type[short] fault : %s\n", name);
    Hostonly_ fprintf(fp, "Error : SBX Data Type[short] fault : %s\n", name);
    assert(0);
  }
  else if (vsz == SklVoxDataSet::TYPE_INT) {
    Hostonly_ printf     ("Error : SBX Data Type[int] fault : %s\n", name);
    Hostonly_ fprintf(fp, "Error : SBX Data Type[int] fault : %s\n", name);
    assert(0);
  }
  
  Hostonly_ printf     ("\tget sbx data : %s\n", name);
  Hostonly_ fprintf(fp, "\tget sbx data : %s\n", name);
  
  if (sbxSet)     { delete sbxSet; sbxSet=NULL; }
  
  Hostonly_ printf     ("\t\n");
  Hostonly_ fprintf(fp, "\t\n");
}

/**
 @fn void FileIO::writeRawSPH(const SKL_REAL *vf, const unsigned* size, const unsigned gc, const SKL_REAL* org, const SKL_REAL* ddx,
                              const unsigned m_ModePrecision) const
 @brief sphファイルの書き出し　ただし，ガイドセルは１
 @param vf スカラデータ
 @param size 配列サイズ
 @param gc ガイドセル
 @param org 基点
 @param ddx ピッチ
 @param m_ModePrecision 浮動小数点の精度
 @note 標記上，long 対応になっているが，ファイルフォーマットとの対応を確認のこと
 */
void FileIO::writeRawSPH(const SKL_REAL *vf, const unsigned* size, const unsigned gc, const SKL_REAL* org, const SKL_REAL* ddx, 
                         const unsigned m_ModePrecision) const
{
  SklParaManager* para_mng = ParaCmpo->GetParaManager();
  
  int sz, dType, stp, svType;
  int ix, jx, kx, i, j, k;
  unsigned long m, l, nx, ix_l, jx_l, kx_l;
  SKL_REAL ox, oy, oz, dx, dy, dz, tm;
  long long szl[3], stpl;
  
  char sph_fname[512];
  
  if( para_mng->IsParallel() ){
    sprintf( sph_fname, "field%010d.sph", pn.ID );
  } else {
    sprintf( sph_fname, "field.sph" );
  }
  
  ofstream ofs(sph_fname, ios::out | ios::binary);
  if (!ofs) {
    cout << "\tCan't open " << sph_fname << " file" << endl;
    assert(0);
  }
  
  ix = size[0]+2*gc;
  jx = size[1]+2*gc;
  kx = size[2]+2*gc;
  ix_l = ix;
  jx_l = jx;
  kx_l = kx;
  nx = ix_l * jx_l * kx_l;
  ox = org[0]-ddx[0]*(SKL_REAL)gc;
  oy = org[1]-ddx[1]*(SKL_REAL)gc;
  oz = org[2]-ddx[2]*(SKL_REAL)gc;
  dx = ddx[0];
  dy = ddx[1];
  dz = ddx[2];
  svType = kind_scalar;
  if ( sizeof(SKL_REAL) == sizeof(double) ) {
    for (i=0; i<3; i++)   szl[i] = (long long)size[i];
  }
  
  SKL_REAL *f = new SKL_REAL[nx];
  
  for (k=0; k<=kx; k++) {
    for (j=0; j<=jx; j++) {
      for (i=0; i<=ix; i++) {
        l = ix_l*jx_l*k + ix_l*j + i;
        m = SklUtil::getFindexS3D(size, gc, i, j, k);
        f[l] = (SKL_REAL)vf[m];
      }
    }
  }
  
  // data property
  ( m_ModePrecision == SPH_SINGLE ) ? dType=1 : dType=2;
  sz = sizeof(unsigned)*2;
  ofs.write( (char*)&sz, sizeof(int) );
  ofs.write( (char*)&svType, sizeof(int) );
  ofs.write( (char*)&dType, sizeof(int) );
  ofs.write( (char*)&sz, sizeof(int) );
  
  // voxel size
  if (dType == 1) {
    sz = sizeof(unsigned)*3;
    ofs.write( (char*)&sz, sizeof(int) );
    ofs.write( (char*)&ix, sizeof(unsigned) );
    ofs.write( (char*)&jx, sizeof(unsigned) );
    ofs.write( (char*)&kx, sizeof(unsigned) );
    ofs.write( (char*)&sz, sizeof(int) );
  }
  else {
    sz = sizeof(long long)*3;
    ofs.write( (char*)&sz, sizeof(int) );
    ofs.write( (char*)&szl[0], sizeof(long long) );
    ofs.write( (char*)&szl[1], sizeof(long long) );
    ofs.write( (char*)&szl[2], sizeof(long long) );
    ofs.write( (char*)&sz, sizeof(int) );
  }
  
  // original point of domain
  if (dType == 1) {
    sz = sizeof(float)*3;
    ofs.write( (char*)&sz, sizeof(int) );
    ofs.write( (char*)&ox, sizeof(float) );
    ofs.write( (char*)&oy, sizeof(float) );
    ofs.write( (char*)&oz, sizeof(float) );
    ofs.write( (char*)&sz, sizeof(int) );
  }
  else {
    sz = sizeof(double)*3;
    ofs.write( (char*)&sz, sizeof(int) );
    ofs.write( (char*)&ox, sizeof(double) );
    ofs.write( (char*)&oy, sizeof(double) );
    ofs.write( (char*)&oz, sizeof(double) );
    ofs.write( (char*)&sz, sizeof(int) );
  }
  
  // pitch of voxel
  if (dType == 1) {
    sz = sizeof(float)*3;
    ofs.write( (char*)&sz, sizeof(int) );
    ofs.write( (char*)&dx, sizeof(float) );
    ofs.write( (char*)&dy, sizeof(float) );
    ofs.write( (char*)&dz, sizeof(float) );
    ofs.write( (char*)&sz, sizeof(int) );
  }
  else {
    sz = sizeof(double)*3;
    ofs.write( (char*)&sz, sizeof(int) );
    ofs.write( (char*)&dx, sizeof(double) );
    ofs.write( (char*)&dy, sizeof(double) );
    ofs.write( (char*)&dz, sizeof(double) );
    ofs.write( (char*)&sz, sizeof(int) );
  }
  
  // time stamp
  if (dType == 1) {
    stp = 0;
    tm = 0.0;
    sz = sizeof(int)+sizeof(float);
    ofs.write( (char*)&sz, sizeof(int) );
    ofs.write( (char*)&stp, sizeof(int) );
    ofs.write( (char*)&tm, sizeof(float) );
    ofs.write( (char*)&sz, sizeof(int) );
  }
  else {
    stpl =0;
    tm = 0.0;
    sz = sizeof(long long)+sizeof(double);
    ofs.write( (char*)&sz, sizeof(int) );
    ofs.write( (char*)&stpl, sizeof(long long) );
    ofs.write( (char*)&tm, sizeof(double) );
    ofs.write( (char*)&sz, sizeof(int) );
  }
  
  if (svType == kind_scalar) {
    sz = (m_ModePrecision == SPH_SINGLE) ? nx * sizeof(float) : nx * sizeof(double);
    ofs.write( (char*)&sz, sizeof(int) );
    ofs.write( (char*)f,   sz );
    ofs.write( (char*)&sz, sizeof(int) );
  }
  
  ofs.close();
  
  if (f) { delete [] f; f=NULL; }
}