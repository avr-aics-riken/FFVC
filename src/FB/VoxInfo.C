// #################################################################
//
// CAERU Library
//
// Copyright (c) All right reserved. 2012
//
// Institute of Industrial Science, The University of Tokyo, Japan. 
//
// #################################################################

//@file VoxInfo.C
//@brief FlowBase VoxInfo class
//@author keno, FSI Team, VCAD, RIKEN

#include <set>
#include <algorithm>
#include <map>
#include "VoxInfo.h"

/**
 @fn void VoxInfo::adjMedium_on_GC(const int face, SklScalar3D<int>* d_mid, const int BCtype, const int c_id, const unsigned prdc_mode)
 @brief 外部境界に接するガイドセルのmid[]に媒質インデクスをエンコードする
 @param face 外部境界面番号
 @param d_mid ID配列のデータクラス
 @param BCtype 外部境界面の境界条件の種類
 @param c_id 媒質インデクス
 @param prdc_mode 周期境界条件のモード
 @note ガイドセル全てを対象
 */
void VoxInfo::adjMedium_on_GC(const int face, SklScalar3D<int>* d_mid, const int BCtype, const int c_id, const unsigned prdc_mode)
{
  unsigned m, m0, m1;
  int* mid=NULL;
  
  int ix = (int)size[0];
  int jx = (int)size[1];
  int kx = (int)size[2];
  int gd = (int)guide;
  
  if ( !(mid = d_mid->GetData()) ) Exit(0);
  
  
  // 周期境界以外
  if ( BCtype != OBC_PERIODIC ) {
    switch (face) {
      case X_MINUS:
        if( pn.nID[face] < 0 ){
          for (int k=1; k<=kx; k++) {
            for (int j=1; j<=jx; j++) {
              m = FBUtility::getFindexS3D(size, guide, 0, j, k);
              mid[m] = c_id;
            }
          }
        }
        break;
        
      case X_PLUS:
        if( pn.nID[face] < 0 ){
          for (int k=1; k<=kx; k++) {
            for (int j=1; j<=jx; j++) {
              m = FBUtility::getFindexS3D(size, guide, ix+1, j, k);
              mid[m] = c_id;
            }
          }
        }
        break;
        
      case Y_MINUS:
        if( pn.nID[face] < 0 ){
          for (int k=1; k<=kx; k++) {
            for (int i=1; i<=ix; i++) {
              m = FBUtility::getFindexS3D(size, guide, i, 0, k);
              mid[m] = c_id;
            }
          }
        }
        break;
        
      case Y_PLUS:
        if( pn.nID[face] < 0 ){
          for (int k=1; k<=kx; k++) {
            for (int i=1; i<=ix; i++) {
              m = FBUtility::getFindexS3D(size, guide, i, jx+1, k);
              mid[m] = c_id;
            }
          }
        }
        break;
        
      case Z_MINUS:
        if( pn.nID[face] < 0 ){
          for (int j=1; j<=jx; j++) {
            for (int i=1; i<=ix; i++) {
              m = FBUtility::getFindexS3D(size, guide, i, j, 0);
              mid[m] = c_id;
            }
          }
        }
        break;
        
      case Z_PLUS:
        if( pn.nID[face] < 0 ){
          for (int j=1; j<=jx; j++) {
            for (int i=1; i<=ix; i++) {
              m = FBUtility::getFindexS3D(size, guide, i, j, kx+1);
              mid[m] = c_id;
            }
          }
        }
        break;
    } // end of switch
  }
  // 周期境界のとき
  else {
    // 内部周期境界の場合には，別メソッド
    if ( prdc_mode != BoundaryOuter::prdc_Driver ) {
      // 並列時
      if ( pn.numProc > 1 ) {
        switch (face) {
          case X_MINUS:
            if ( !d_mid->CommPeriodicBndCell(PRDC_X_DIR, PRDC_PLUS2MINUS, guide) ) Exit(0);
            break;
            
          case X_PLUS:
            if ( !d_mid->CommPeriodicBndCell(PRDC_X_DIR, PRDC_MINUS2PLUS, guide) ) Exit(0);
            break;
            
          case Y_MINUS:
            if ( !d_mid->CommPeriodicBndCell(PRDC_Y_DIR, PRDC_PLUS2MINUS, guide) ) Exit(0);
            break;
            
          case Y_PLUS:
            if ( !d_mid->CommPeriodicBndCell(PRDC_Y_DIR, PRDC_MINUS2PLUS, guide) ) Exit(0);
            break;
            
          case Z_MINUS:
            if ( !d_mid->CommPeriodicBndCell(PRDC_Z_DIR, PRDC_PLUS2MINUS, guide) ) Exit(0);
            break;
            
          case Z_PLUS:
            if ( !d_mid->CommPeriodicBndCell(PRDC_Z_DIR, PRDC_MINUS2PLUS, guide) ) Exit(0);
            break;
        }      
      }
      // 非並列時
      else {
        switch (face) {
          case X_MINUS:
            if( pn.nID[face] < 0 ){
              for (int k=1; k<=kx; k++) {
                for (int j=1; j<=jx; j++) {
                  for (int i=1-gd; i<=0; i++) {
                    m0 = FBUtility::getFindexS3D(size, guide, i   , j, k);
                    m1 = FBUtility::getFindexS3D(size, guide, i+ix, j, k);
                    mid[m0] = mid[m1];
                  }
                }
              }
            }
            break;
            
          case X_PLUS:
            if( pn.nID[face] < 0 ){
              for (int k=1; k<=kx; k++) {
                for (int j=1; j<=jx; j++) {
                  for (int i=ix+1; i<=ix+gd; i++) {
                    m0 = FBUtility::getFindexS3D(size, guide, i   , j, k);
                    m1 = FBUtility::getFindexS3D(size, guide, i-ix, j, k);
                    mid[m0] = mid[m1];
                  }
                }
              }
            }
            break;
            
          case Y_MINUS:
            if( pn.nID[face] < 0 ){
              for (int k=1; k<=kx; k++) {
                for (int j=1-gd; j<=0; j++) {
                  for (int i=1; i<=ix; i++) {
                    m0 = FBUtility::getFindexS3D(size, guide, i, j   , k);
                    m1 = FBUtility::getFindexS3D(size, guide, i, j+jx, k);
                    mid[m0] = mid[m1];
                  }
                }
              }
            }
            break;
            
          case Y_PLUS:
            if( pn.nID[face] < 0 ){
              for (int k=1; k<=kx; k++) {
                for (int j=jx+1; j<=jx+gd; j++) {
                  for (int i=1; i<=ix; i++) {
                    m0 = FBUtility::getFindexS3D(size, guide, i, j   , k);
                    m1 = FBUtility::getFindexS3D(size, guide, i, j-jx, k);
                    mid[m0] = mid[m1];
                  }
                }
              }
            }
            break;
            
          case Z_MINUS:
            if( pn.nID[face] < 0 ){
              for (int k=1-gd; k<=0; k++) {
                for (int j=1; j<=jx; j++) {
                  for (int i=1; i<=ix; i++) {
                    m0 = FBUtility::getFindexS3D(size, guide, i, j, k    );
                    m1 = FBUtility::getFindexS3D(size, guide, i, j, k+kx);
                    mid[m0] = mid[m1];
                  }
                }
              }
            }
            break;
            
          case Z_PLUS:
            if( pn.nID[face] < 0 ){
              for (int k=kx+1; k<=kx+gd; k++) {
                for (int j=1; j<=jx; j++) {
                  for (int i=1; i<=ix; i++) {
                    m0 = FBUtility::getFindexS3D(size, guide, i, j, k    );
                    m1 = FBUtility::getFindexS3D(size, guide, i, j, k-kx);
                    mid[m0] = mid[m1];
                  }
                }
              }
            }
            break;
        }    
      }
    }
  }
  
}

/**
 @fn void VoxInfo::adjMediumPrdc_Inner(int face, SklScalar3D<int>* d_mid, int BCtype, unsigned c_id, unsigned prdc_mode)
 @brief 外部境界に接するガイドセルのmid[]にIDをエンコードする（内部周期境界の場合）
 @param face 外部境界面番号
 @param d_mid ID配列のデータクラス
 @param BCtype 外部境界面の境界条件の種類
 @param c_id セルID
 @param prdc_mode 周期境界条件のモード
 @note ガイドセル全てを対象
 */
void VoxInfo::adjMediumPrdc_Inner(SklScalar3D<int>* d_mid)
{
  int st[3], ed[3], dir, id;
  
  for (unsigned n=1; n<=NoBC; n++) {
    cmp[n].getBbox(st, ed);
    dir = (int)cmp[n].getPeriodicDir();
    id  = cmp[n].getMatOdr();
    
    if ( cmp[n].getType() == PERIODIC ) {

      if ( pn.numProc > 1 ) {
        Hostonly_ printf("Error : Inner Periodic condition is limited to use for serial execution on a temporary\n.");
        Exit(0);
      }
      copyID_Prdc_Inner(d_mid, st, ed, id, dir);
    }
  }
  
}


/**
 @fn int* VoxInfo::allocTable(unsigned size)
 @brief モデル情報テーブルを作成，0で初期化する
 @retval エラーの場合はNULL
 @param size アロケートするサイズ
 */
int* VoxInfo::allocTable(unsigned size)
{
  int* table=NULL;
  if( (table = (int*)malloc(size*sizeof(int))) == NULL ) return NULL;
  for (unsigned i=0; i<size; i++) {
    table[i]=0;
  }
  return table;
}

/**
 @fn void VoxInfo::alloc_voxel_nv(unsigned len)
 @brief ボクセルモデルの法線試算用の配列を確保
 @param len 配列長
 */
void VoxInfo::alloc_voxel_nv(unsigned len)
{
  vox_nv = new REAL_TYPE[len];
  if ( !vox_nv ) Exit(0);
}


/**
 @fn void VoxInfo::cal_Compo_Area_Normal(unsigned n, unsigned* bd, unsigned* bv, unsigned* bh1, REAL_TYPE dhd, int* gi)
 @brief コンポーネントの断面積と法線を求める
 @param n エントリ番号
 @param bd BCindex ID
 @param bv BCindex V
 @param bh1 BCindex H1
 @param dhd 有次元格子幅
 @param gi コンポーネントのグローバルインデクス
 */
void VoxInfo::cal_Compo_Area_Normal(unsigned n, unsigned* bd, unsigned* bv, unsigned* bh1, REAL_TYPE dhd, int* gi)
{
  REAL_TYPE a, nvx, nvy, nvz;
  int cijk[3], dir[3], area=0;
  REAL_TYPE ai, aj, ak;
  
  cijk[0] = cijk[1] = cijk[2] = 0;
  dir[0] = dir[1] = dir[2] = 0;
  
  int id  = cmp[n].getMatOdr();
  int def = cmp[n].getDef();
  
  switch ( cmp[n].getType() ) {
    case SPEC_VEL:
    case SPEC_VEL_WH:
    case OUTFLOW:
      return;
      break;
      
      
    case HEATFLUX:
    case TRANSFER:
    case ISOTHERMAL:
      countNrml_from_FaceBC(n, bh1, cijk, area);
      ai = (REAL_TYPE)cijk[0];
      aj = (REAL_TYPE)cijk[1];
      ak = (REAL_TYPE)cijk[2];
      break;
      
    case HEX:
    case FAN:
      return;
      break;
      
    case CELL_MONITOR:
    case DARCY:
      break;
  }
  
  a = sqrt( ai*ai + aj*aj + ak*ak );
  if ( a == 0.0 ) a = (REAL_TYPE)area;
  
  if ( a != 0.0 ) {
    nvx = ai / a;
    nvy = aj / a;
    nvz = ak / a;
  }
  else {
    nvx = nvy = nvz = 0.0;
  }
  
  if ( cmp[n].isFORCING() || cmp[n].isMONITOR() ) {
    nvx *= (REAL_TYPE)dir[0];
    nvy *= (REAL_TYPE)dir[1];
    nvz *= (REAL_TYPE)dir[2];
  }
  
  // following 3 lines are to avoid -0.0 for display and to take acdount suction/blowing
  vox_nv[3*n+0] = ( (cijk[0]==0) ? 0.0 : nvx );
  vox_nv[3*n+1] = ( (cijk[1]==0) ? 0.0 : nvy );
  vox_nv[3*n+2] = ( (cijk[2]==0) ? 0.0 : nvz );
  
  cmp[n].area = a * dhd * dhd;  // dhd は有次元値，近似的な面積
}

/**
 @fn void VoxInfo::checkColorTable(FILE* fp, unsigned size, int* table)
 @brief IDテーブルを表示
 @param fp 
 @param size 
 @param table ID table
 @note - デバッグ用
 */
void VoxInfo::checkColorTable(FILE* fp, unsigned size, int* table)
{
  for (unsigned i=0; i<size; i++) {
    Hostonly_ fprintf(fp, "\t\t%d : %d\n", i, table[i]);
  }
  fflush(fp);
}


/**
 @fn unsigned VoxInfo::check_fill(const int* mid)
 @brief ペイント済みかどうかをチェックする
 @param[in] mid ID配列
 @note 未ペイントセルがあれば1を返す
 */
unsigned VoxInfo::check_fill(const int* mid)
{
  for (int k=1; k<=(int)size[2]; k++) {
    for (int j=1; j<=(int)size[1]; j++) {
      for (int i=1; i<=(int)size[0]; i++) {
        
        if ( mid[FBUtility::getFindexS3D(size, guide, i  , j  , k  )] == 0 ) return 1;
        
      }
    }
  }
  
  return 0;
}


/**
 @fn bool VoxInfo::chkIDconsistency(const int m_NoMedium)
 @brief パラメータファイルとスキャンしたIDの同一性をチェック
 @param m_NoMedium Medium_Tableに記述されたIDの個数
 */
bool VoxInfo::chkIDconsistency(const int m_NoMedium)
{
  bool* chkflag = NULL;
	if ( !(chkflag = new bool[NoVoxID+1]) ) return false;
	
	for (int i=0; i<=NoVoxID; i++) chkflag[i] = false;
	
  for (int j=1; j<=m_NoMedium; j++) {
    for (int i=1; i<=NoVoxID; i++) { // スキャンしたIDのループ
			if ( colorList[i] == j ) chkflag[i] = true;
		}
	}
	
	for (int i=1; i<=NoVoxID; i++) {
		if ( !chkflag[i] ) return false;
	}
  
  if( chkflag ) delete [] chkflag;
	
	return true;
}

/**
 @fn bool VoxInfo::chkIDinside(unsigned id, int* mid, unsigned* bx)
 @brief 指定されたIDが計算領域内部にあるかどうかを判定する
 @param id サーチ対象ID
 @param mid ID配列
 @param bx BC Index
 @retval IDがあればtrue
 */
bool VoxInfo::chkIDinside(unsigned id, int* mid, unsigned* bx)
{
  int i, j, k;
  unsigned m;
  
  int ix = (int)size[0];
  int jx = (int)size[1];
  int kx = (int)size[2];
  int gd = (int)guide;
  
  for (k=1; k<=kx; k++) {
    for (j=1; j<=jx; j++) {
      for (i=1; i<=ix; i++) {
        m = FBUtility::getFindexS3D(size, guide, i, j, k);
        if ( mid[m] == (int)id ) {
          Hostonly_ stamped_printf("\tID[%d] exist inside the computational region\n", id);
          return true;
        }
      }
    }
  }
  
  return false;
}

/**
 @fn void VoxInfo::copyBCIbase(unsigned* dst, unsigned* src)
 @brief dst[]にsrc[]のstate, activeビットの情報をコピーする
 @param dst
 @param src
 */
void VoxInfo::copyBCIbase(unsigned* dst, unsigned* src)
{
  unsigned nx, m;
  
  nx = (size[0]+2*guide) * (size[1]+2*guide) * (size[2]+2*guide);
  
  // 30, 31ビットのみコピー
  for (m=0; m<nx; m++) {
    dst[m] = src[m] & 0xc0000000;
  }  
}

/**
 @fn void VoxInfo::copyID_Prdc_Inner(SklScalar3D<int>* d_mid, int* st, int* ed, int id, int dir)
 @brief 外部境界に接するガイドセルのmid[]にIDを内部周期境界からコピーする
 @param d_mid ID配列のデータクラス
 @param st コンポーネントのbbox始点
 @param ed コンポーネントのbbox終点
 @param id 対象のID
 @param dir ドライバの方向
 */
void VoxInfo::copyID_Prdc_Inner(SklScalar3D<int>* d_mid, int* st, int* ed, int id, int dir)
{
  int i, j, k, ii, jj, kk;
  unsigned m0, m1, m2;
  int* mid=NULL;
  
  int ix = (int)size[0];
  int jx = (int)size[1];
  int kx = (int)size[2];
  int gd = (int)guide;
  
  if ( !(mid = d_mid->GetData()) ) Exit(0);
  
  switch (dir) {
    case X_MINUS:
      if ( pn.nID[dir] < 0 ) {
        i = st[0];
        for (k=st[2]; k<=ed[2]; k++) {
          for (j=st[1]; j<=ed[1]; j++) {
            m1 = FBUtility::getFindexS3D(size, guide, i, j, k);
            if ( mid[m1] == id ) {
              for (ii=1-gd; ii<=0; ii++) {
                m0 = FBUtility::getFindexS3D(size, guide, ii,   j, k);
                m2 = FBUtility::getFindexS3D(size, guide, i+ii, j, k);
                mid[m0] = mid[m2];
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
            m1 = FBUtility::getFindexS3D(size, guide, i, j, k);
            if ( mid[m1] == id ) {
              for (ii=ix+1; ii<=ix+gd; ii++) {
                m0 = FBUtility::getFindexS3D(size, guide, ii,        j, k);
                m2 = FBUtility::getFindexS3D(size, guide, i+ii-ix-1, j, k);
                mid[m0] = mid[m2];
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
            m1 = FBUtility::getFindexS3D(size, guide, i, j, k);
            if ( mid[m1] == id ) {
              for (jj=1-gd; jj<=0; jj++) {
                m0 = FBUtility::getFindexS3D(size, guide, i, jj,   k);
                m2 = FBUtility::getFindexS3D(size, guide, i, j+jj, k);
                mid[m0] = mid[m2];
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
            m1 = FBUtility::getFindexS3D(size, guide, i, j, k);
            if ( mid[m1] == id ) {
              for (jj=jx+1; jj<=jx+gd; jj++) {
                m0 = FBUtility::getFindexS3D(size, guide, i, jj,        k);
                m2 = FBUtility::getFindexS3D(size, guide, i, j+jj-jx-1, k);
                mid[m0] = mid[m2];
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
            m1 = FBUtility::getFindexS3D(size, guide, i, j, k);
            if ( mid[m1] == id ) {
              for (kk=1-gd; kk<=0; kk++) {
                m0 = FBUtility::getFindexS3D(size, guide, i, j, kk  );
                m2 = FBUtility::getFindexS3D(size, guide, i, j, k+kk);
                mid[m0] = mid[m2];
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
            m1 = FBUtility::getFindexS3D(size, guide, i, j, k);
            if ( mid[m1] == id ) {
              for (kk=kx+1; kk<=kx+gd; kk++) {
                m0 = FBUtility::getFindexS3D(size, guide, i, j, kk       );
                m2 = FBUtility::getFindexS3D(size, guide, i, j, k+kk-kx-1);
                mid[m0] = mid[m2];
              }
            }
          }
        }
      }
      break;
  }  
  
}

/**
 @fn void VoxInfo::countCellState(unsigned& Lcell, unsigned& Gcell, unsigned* bx, const unsigned state)
 @brief セルの状態をカウントして，その個数をLcell, Gcellに保持する
 @param Lcell ノードローカルの値（戻り値）
 @param Gcell グローバルの値（戻り値）
 @param bx BCindex ID
 @param state カウントするセルの状態
 */
void VoxInfo::countCellState(unsigned long& Lcell, unsigned long& Gcell, unsigned* bx, const unsigned state)
{
  unsigned long cell=0;    // local
  unsigned long g_cell=0;  // global 

  unsigned m;
  
  int ix = (int)size[0];
  int jx = (int)size[1];
  int kx = (int)size[2];
  int gd = (int)guide;
  
  // described in Fortran index
  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        m = FBUtility::getFindexS3D(size, guide, i, j, k);

        if ( state == SOLID) {
          if (!IS_FLUID(bx[m]) ) cell++;  //  IS_Fluid() > 0=SOLID, 1=FLUID
        }
        else {
          if ( IS_FLUID(bx[m]) ) cell++;
        }
      }
    }
  }
  
  g_cell = cell;
  Lcell  = cell;
  
  if ( pn.numProc > 1 ) {
    unsigned long c_tmp = g_cell;
    MPI_Allreduce(&c_tmp, &g_cell, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  }
  
  Gcell = g_cell;

}


/**
 @fn void VoxInfo::countNrml_from_FaceBC(unsigned n, unsigned* bx, int* cc, int& ar)
 @brief 指定面の断面積と法線を求める
 @param n 境界条件番号
 @param bx BCindex
 @param cc[out] カウントしたセル数
 @param ar[out]
 @note
 - 法線の方向は，Fセルに向かう方向にとる
 - F-S界面のみ法線を近似．S-S界面の場合には常に正として処理し，面積のみ計算する
 - 閉じている物体は法線はゼロになる
 */
void VoxInfo::countNrml_from_FaceBC(unsigned n, unsigned* bx, int* cc, int& ar)
{
  unsigned m_p, m_e, m_w, m_n, m_s, m_t, m_b, s;
  REAL_TYPE c_p, c_e, c_w, c_n, c_s, c_t, c_b;
  REAL_TYPE dw, de, ds, dn, db, dt;
  int i,j,k, c[3];
	int st[3], ed[3];
  
  c[0] = c[1] = c[2] = 0;
  ar = 0;
  
  cmp[n].getBbox(st, ed);
  
  for (k=st[2]; k<=ed[2]; k++) {
    for (j=st[1]; j<=ed[1]; j++) {
      for (i=st[0]; i<=ed[0]; i++) {
        m_p = FBUtility::getFindexS3D(size, guide, i  , j  , k  );
        s = bx[m_p];
        
        if ( IS_INCLUDE_BC(s) ) { // セルの6面のいずれかにBCが設定されている
          m_e = FBUtility::getFindexS3D(size, guide, i+1, j  , k  );
          m_w = FBUtility::getFindexS3D(size, guide, i-1, j  , k  );
          m_n = FBUtility::getFindexS3D(size, guide, i  , j+1, k  );
          m_s = FBUtility::getFindexS3D(size, guide, i  , j-1, k  );
          m_t = FBUtility::getFindexS3D(size, guide, i  , j  , k+1);
          m_b = FBUtility::getFindexS3D(size, guide, i  , j  , k-1);
          
          // Fluid=1.0, Solod=0.0
          c_p = GET_SHIFT_F(s      , STATE_BIT);
          c_e = GET_SHIFT_F(bx[m_e], STATE_BIT);
          c_w = GET_SHIFT_F(bx[m_w], STATE_BIT);
          c_n = GET_SHIFT_F(bx[m_n], STATE_BIT);
          c_s = GET_SHIFT_F(bx[m_s], STATE_BIT);
          c_t = GET_SHIFT_F(bx[m_t], STATE_BIT);
          c_b = GET_SHIFT_F(bx[m_b], STATE_BIT);
          
          de = c_e - c_p;
          dn = c_n - c_p;
          dt = c_t - c_p;
          dw = c_p - c_w;
          ds = c_p - c_s;
          db = c_p - c_b;
          
          if ( GET_FACE_BC(s, BC_FACE_W) == n ) {
            c[0] += dw; // W面にBCが設定されている場合，dwがFセルに向かう法線方向を表す
            ar++;
          }
          
          if ( GET_FACE_BC(s, BC_FACE_E) == n ) {
            c[0] += de;
            ar++;
          }
          
          if ( GET_FACE_BC(s, BC_FACE_S) == n ) {
            c[1] += ds;
            ar++;
          }
          
          if ( GET_FACE_BC(s, BC_FACE_N) == n ) {
            c[1] += dn;
            ar++;
          }
          
          if ( GET_FACE_BC(s, BC_FACE_B) == n ) {
            c[2] += db;
            ar++;
          }
          
          if ( GET_FACE_BC(s, BC_FACE_T) == n ) {
            c[2] += dt;
            ar++;
          }
        }
        
      }
    }
  }
  
  
  if ( pn.numProc > 1 ) {
    int tmp[3];
    tmp[0] = c[0];
		tmp[1] = c[1];
		tmp[2] = c[2]; 
    MPI_Allreduce(tmp, c, 3, MPI_INT, MPI_SUM, MPI_COMM_WORLD); // このメソッドはなくなる
    
    tmp[0] = ar;
    MPI_Allreduce(tmp, &ar, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  }
  
	cc[0] = c[0];
	cc[1] = c[1];
	cc[2] = c[2];
}

/**
 @fn unsigned long VoxInfo::count_ValidCell_OBC(const int face, const unsigned* bv)
 @brief 外部境界面の有効セル数をカウントする
 @param face 外部境界面番号
 @param bv BCindex V
 @note 外部境界面の両側のセルがFluidのときのみカウント
 */
unsigned long VoxInfo::count_ValidCell_OBC(const int face, const unsigned* bv)
{
  unsigned m1, m2;
  unsigned register s1, s2;
  
  int ix = (int)size[0];
  int jx = (int)size[1];
  int kx = (int)size[2];
  int gd = (int)guide;
  
  unsigned long g=0;
  
  switch (face) {
    case X_MINUS:
      if( pn.nID[X_MINUS] < 0 ){ // 外部境界をもつノードのみ
        for (int k=1; k<=kx; k++) {
          for (int j=1; j<=jx; j++) {
            m1 = FBUtility::getFindexS3D(size, guide, 1  , j, k);
            m2 = FBUtility::getFindexS3D(size, guide, 0  , j, k);
            s1 = bv[m1];
            s2 = bv[m2];
            if ( GET_FACE_BC(s1, BC_FACE_W) == OBC_MASK ) {
              if ( IS_FLUID(s1) && IS_FLUID(s2) ) g++;
            }
          }
        }        
      }
      break;
      
    case X_PLUS:
      if( pn.nID[X_PLUS] < 0 ){
        for (int k=1; k<=kx; k++) {
          for (int j=1; j<=jx; j++) {
            m1 = FBUtility::getFindexS3D(size, guide, ix  , j, k);
            m2 = FBUtility::getFindexS3D(size, guide, ix+1, j, k);
            s1 = bv[m1];
            s2 = bv[m2];
            if ( GET_FACE_BC(s1, BC_FACE_E) == OBC_MASK ) {
              if ( IS_FLUID(s1) && IS_FLUID(s2) ) g++;
            }
          }
        }
      }
      break;
      
    case Y_MINUS:
      if( pn.nID[Y_MINUS] < 0 ){
        for (int k=1; k<=kx; k++) {
          for (int i=1; i<=ix; i++) {
            m1 = FBUtility::getFindexS3D(size, guide, i, 1  , k);
            m2 = FBUtility::getFindexS3D(size, guide, i, 0  , k);
            s1 = bv[m1];
            s2 = bv[m2];
            if ( GET_FACE_BC(s1, BC_FACE_S) == OBC_MASK ) {
              if ( IS_FLUID(s1) && IS_FLUID(s2) ) g++;
            }
          }
        }
      }
      break;
      
    case Y_PLUS:
      if( pn.nID[Y_PLUS] < 0 ){
        for (int k=1; k<=kx; k++) {
          for (int i=1; i<=ix; i++) {
            m1 = FBUtility::getFindexS3D(size, guide, i, jx  , k);
            m2 = FBUtility::getFindexS3D(size, guide, i, jx+1, k);
            s1 = bv[m1];
            s2 = bv[m2];
            if ( GET_FACE_BC(s1, BC_FACE_N) == OBC_MASK ) {
              if ( IS_FLUID(s1) && IS_FLUID(s2) ) g++;
            }
          }
        }
      }
      break;
      
    case Z_MINUS:
      if( pn.nID[Z_MINUS] < 0 ){
        for (int j=1; j<=jx; j++) {
          for (int i=1; i<=ix; i++) {
            m1 = FBUtility::getFindexS3D(size, guide, i, j, 1  );
            m2 = FBUtility::getFindexS3D(size, guide, i, j, 0  );
            s1 = bv[m1];
            s2 = bv[m2];
            if ( GET_FACE_BC(s1, BC_FACE_B) == OBC_MASK ) {
              if ( IS_FLUID(s1) && IS_FLUID(s2) ) g++;
            }
          }
        }
      }
      break;
      
    case Z_PLUS:
      if( pn.nID[Z_PLUS] < 0 ){
        for (int j=1; j<=jx; j++) {
          for (int i=1; i<=ix; i++) {
            m1 = FBUtility::getFindexS3D(size, guide, i, j, kx  );
            m2 = FBUtility::getFindexS3D(size, guide, i, j, kx+1);
            s1 = bv[m1];
            s2 = bv[m2];
            if ( GET_FACE_BC(s1, BC_FACE_T) == OBC_MASK ) {
              if ( IS_FLUID(s1) && IS_FLUID(s2) ) g++;
            }
          }
        }
      }
      break;
  } // end of switch
  
  
  if ( pn.numProc > 1 ) {
    unsigned long tmp = g;
    MPI_Allreduce(&tmp, &g, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  }
  
  return g;
}



/**
 @fn void VoxInfo::countOpenAreaOfDomain(unsigned* bx, REAL_TYPE* OpenArea)
 @brief  計算領域の外部境界で外側1層と内側の両方が流体セル数の場合にカウントする
 @param bx BCindex ID
 @param OpenArea 開口セル数
 */
void VoxInfo::countOpenAreaOfDomain(unsigned* bx, REAL_TYPE* OpenArea)
{
  unsigned m0, m1, g;
  unsigned m_area[NOFACE];
  
  int ix = (int)size[0];
  int jx = (int)size[1];
  int kx = (int)size[2];
  int gd = (int)guide;
  
  for (int i=0; i<NOFACE; i++) {
    OpenArea[i]=0.0;
    m_area[i] = 0.0;
  }
  
  // described in Fortran index
  
  // X_MINUS
  g=0;
  if( pn.nID[X_MINUS] < 0 ){
    for (int k=1; k<=kx; k++) {
      for (int j=1; j<=jx; j++) {
        m0 = FBUtility::getFindexS3D(size, guide, 0, j, k);
        m1 = FBUtility::getFindexS3D(size, guide, 1, j, k);
        if ( (  FLUID == IS_FLUID(bx[m0]))
            && (FLUID == IS_FLUID(bx[m1])) ) g++;
      }
    }
    m_area[X_MINUS] = g;
  }
  
  // X_PLUS
  g=0;
  if( pn.nID[X_PLUS] < 0 ){
    for (int k=1; k<=kx; k++) {
      for (int j=1; j<=jx; j++) {
        m0 = FBUtility::getFindexS3D(size, guide, ix+1, j, k);
        m1 = FBUtility::getFindexS3D(size, guide, ix,   j, k);
        if ( (  FLUID == IS_FLUID(bx[m0]))
            && (FLUID == IS_FLUID(bx[m1])) ) g++;
      }
    }
    m_area[X_PLUS] = g;
  }
  
  // Y_MINUS
  g=0;
  if( pn.nID[Y_MINUS] < 0 ){
    for (int k=1; k<=kx; k++) {
      for (int i=1; i<=ix; i++) {
        m0 = FBUtility::getFindexS3D(size, guide, i, 0, k);
        m1 = FBUtility::getFindexS3D(size, guide, i, 1, k);
        if ( (  FLUID == IS_FLUID(bx[m0]))
            && (FLUID == IS_FLUID(bx[m1])) ) g++;
      }
    }
    m_area[Y_MINUS] = g;
  }
  
  // Y_PLUS
  g=0;
  if( pn.nID[Y_PLUS] < 0 ){
    for (int k=1; k<=kx; k++) {
      for (int i=1; i<=ix; i++) {
        m0 = FBUtility::getFindexS3D(size, guide, i, jx+1, k);
        m1 = FBUtility::getFindexS3D(size, guide, i, jx,   k);
        if ( (  FLUID == IS_FLUID(bx[m0]))
            && (FLUID == IS_FLUID(bx[m1])) ) g++;
      }
    }
    m_area[Y_PLUS] = g;
  }
  
  // Z_MINUS
  g=0;
  if( pn.nID[Z_MINUS] < 0 ){
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        m0 = FBUtility::getFindexS3D(size, guide, i, j, 0);
        m1 = FBUtility::getFindexS3D(size, guide, i, j, 1);
        if ( (  FLUID == IS_FLUID(bx[m0]))
            && (FLUID == IS_FLUID(bx[m1])) ) g++;
      }
    }
    m_area[Z_MINUS] = g;
  }
  
  // Z_PLUS
  g=0;
  if( pn.nID[Z_PLUS] < 0 ){
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        m0 = FBUtility::getFindexS3D(size, guide, i, j, kx  );
        m1 = FBUtility::getFindexS3D(size, guide, i, j, kx+1);
        if ( (  FLUID == IS_FLUID(bx[m0]))
            && (FLUID == IS_FLUID(bx[m1])) ) g++;
      }
    }
    m_area[Z_PLUS] = g;
  }
  
  // 面素がunsignedの値域を超えることはないと仮定
  if ( pn.numProc > 1 ) {
    unsigned tmp[NOFACE];
    for (int i=0; i<NOFACE; i++) tmp[i] = m_area[i];
    MPI_Allreduce(tmp, m_area, NOFACE, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
  }
  
  for (int i=0; i<NOFACE; i++) OpenArea[i] = (REAL_TYPE)m_area[i];
}


/**
 @fn unsigned VoxInfo::countState(unsigned id, int* mid)
 @brief 媒質idの数を数え，値を返す
 @retval 計算空間内における媒質idの数
 @param id カウントするid
 @param mid ボクセルID配列
 */
unsigned long VoxInfo::countState(unsigned id, int* mid)
{
  unsigned m;
  unsigned long g=0;
  
  int ix = (int)size[0];
  int jx = (int)size[1];
  int kx = (int)size[2];
  int gd = (int)guide;
  
  // サーチ範囲はノードローカルの計算セル内
  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        m = FBUtility::getFindexS3D(size, guide, i, j, k);
        if ( mid[m] == (int)id )  g++;
      }
    }
  }
  
  
  if ( pn.numProc > 1 ) {
    unsigned long tmp = g;
    MPI_Allreduce(&tmp, &g, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  }
  
  return g;
}

/**
 @fn void VoxInfo::dbg_chkBCIndexD(unsigned* bcd, const char* fname)
 @brief BCindexを表示する（デバッグ用）
 @param bcv BCindex D
 @param fname 出力用ファイル名
 */
void VoxInfo::dbg_chkBCIndexD(unsigned* bcd, const char* fname)
{
  int i,j,k;
  unsigned m, s, d;
  FILE *fp=NULL;
  
  if ( !(fp=fopen(fname,"w")) ) {
    Hostonly_ printf("\tSorry, can't open '%s', write failed.\n", fname);
    Exit(0);
  }
  Hostonly_ printf("\n\tCheck BC Index '%s' for check\n\n", fname);
  
  int ix = (int)size[0];
  int jx = (int)size[1];
  int kx = (int)size[2];
  int gd = (int)guide;
  
  // ガイドセルを含む全領域
  for (k=1-gd; k<=kx+gd; k++) {
    for (j=1-gd; j<=jx+gd; j++) {
      for (i=1-gd; i<=ix+gd; i++) {
        m = FBUtility::getFindexS3D(size, guide, i  , j  , k  );
        s = bcd[m];
        Hostonly_ fprintf(fp, "[%4d %4d %4d], state=%1d: cmp=%3d  ID=%3d mat=%3d  vf=%3d force=%3d\n", 
                          i, j, k, IS_FLUID(s), 
                          DECODE_CMP(s),
                          DECODE_ID(s), 
                          DECODE_MAT(s), 
                          DECODE_VF(s), 
                          (s>>FORCING_BIT)&0x1);
      }
    }
  }
  fflush(fp);
  fclose(fp);
}


/**
 @fn void VoxInfo::dbg_chkBCIndexP(unsigned* bcd, unsigned* bcp, const char* fname)
 @brief BCindexを表示する（デバッグ用）
 @param bcd BCindex ID
 @param bcp BCindex P
 @param fname 出力用ファイル名
 */
void VoxInfo::dbg_chkBCIndexP(unsigned* bcd, unsigned* bcp, const char* fname)
{
  unsigned m, q, s, d;
  int id;
  int i,j,k;
  FILE *fp=NULL;
  
  if ( !(fp=fopen(fname,"w")) ) {
    Hostonly_ printf("\tSorry, can't open '%s', write failed.\n", fname);
    Exit(0);
  }
  
  Hostonly_ printf("\n\tCheck BC Index '%s' for check\n\n", fname);
  
  int ix = (int)size[0];
  int jx = (int)size[1];
  int kx = (int)size[2];
  int gd = (int)guide;
  
  // ガイドセルを含む全領域
  for (k=1-gd; k<=kx+gd; k++) {
    for (j=1-gd; j<=jx+gd; j++) {
      for (i=1-gd; i<=ix+gd; i++) {
        m = FBUtility::getFindexS3D(size, guide, i  , j  , k  );
        d = bcd[m];
        s = bcp[m];
        q = d & MASK_6;
        id = cmp[q].getMatOdr();
        Hostonly_ fprintf(fp, "[%4d %4d %4d], cmp[%3d]=%3d, mat[%3d], state=%1d: D(ewnstb) [%d %d %d %d %d %d] N [%d %d %d %d %d %d] NDAG [%d %d %d %d %d %d] DIAG=%1d :idx=%d\n", 
                          i, j, k, q, id, DECODE_MAT(d), IS_FLUID(d), 
                          (s>>BC_D_E)&0x1, (s>>BC_D_W)&0x1, (s>>BC_D_N)&0x1, (s>>BC_D_S)&0x1, (s>>BC_D_T)&0x1, (s>>BC_D_B)&0x1, 
                          (s>>BC_N_E)&0x1, (s>>BC_N_W)&0x1, (s>>BC_N_N)&0x1, (s>>BC_N_S)&0x1, (s>>BC_N_T)&0x1, (s>>BC_N_B)&0x1, 
                          (s>>BC_NDAG_E)&0x1, (s>>BC_NDAG_W)&0x1, (s>>BC_NDAG_N)&0x1, (s>>BC_NDAG_S)&0x1, (s>>BC_NDAG_T)&0x1, (s>>BC_NDAG_B)&0x1, 
                          (s>>BC_DIAG)&0x7, s);
      }
    }
  }
  fflush(fp);
  fclose(fp);
}

/**
 @fn void VoxInfo::dbg_chkBCIndexV(unsigned* bcv, const char* fname)
 @brief BCindexを表示する（デバッグ用）
 @param bcv BCindex V
 @param fname 出力用ファイル名
 */
void VoxInfo::dbg_chkBCIndexV(unsigned* bcv, const char* fname)
{
  int i,j,k;
  unsigned m, s, d;
  FILE *fp=NULL;
  
  if ( !(fp=fopen(fname,"w")) ) {
    Hostonly_ printf("\tSorry, can't open '%s', write failed.\n", fname);
    Exit(0);
  }
  Hostonly_ printf("\n\tCheck BC Index '%s' for check\n\n", fname);
  
  int ix = (int)size[0];
  int jx = (int)size[1];
  int kx = (int)size[2];
  int gd = (int)guide;
  
  // ガイドセルを含む全領域
  for (k=1-gd; k<=kx+gd; k++) {
    for (j=1-gd; j<=jx+gd; j++) {
      for (i=1-gd; i<=ix+gd; i++) {
        m = FBUtility::getFindexS3D(size, guide, i  , j  , k  );
        s = bcv[m];
        Hostonly_ fprintf(fp, "[%4d %4d %4d], state=%1d: VBC(ewnstb) [%2d %2d %2d %2d %2d %2d] : idx=%d\n", 
                          i, j, k, IS_FLUID(s), 
                          (s>>BC_FACE_E)&0x1f, // 5-bit
                          (s>>BC_FACE_W)&0x1f, 
                          (s>>BC_FACE_N)&0x1f, 
                          (s>>BC_FACE_S)&0x1f, 
                          (s>>BC_FACE_T)&0x1f, 
                          (s>>BC_FACE_B)&0x1f, s);
      }
    }
  }
  fflush(fp);
  fclose(fp);
}


/**
 @fn void VoxInfo::encActive(unsigned& Lcell, unsigned& Gcell, unsigned* bx, unsigned KOS)
 @brief BCindexにそのセルが計算に有効(active)かどうかをエンコードする
 @param Lcell[out] ノードローカルの有効セル数
 @param Gcell[out] グローバルの有効セル数
 @param bx BCindex ID
 @param KOS 解くべき方程式の種類 KIND_OF_SOLVER
 @note
 - IS_FLUID returns true if FLUID
 - 設定は内部領域のみ
 */
void VoxInfo::encActive(unsigned long& Lcell, unsigned long& Gcell, unsigned* bx, const unsigned KOS)
{
  unsigned m;
  unsigned register s;
  
  unsigned long c=0;
  
  int ix = (int)size[0];
  int jx = (int)size[1];
  int kx = (int)size[2];
  int gd = (int)guide;
  
  switch ( KOS ) {
    case FLOW_ONLY:
    case THERMAL_FLOW:
    case THERMAL_FLOW_NATURAL:
      for (int k=1; k<=kx; k++) {
        for (int j=1; j<=jx; j++) {
          for (int i=1; i<=ix; i++) {
            m = FBUtility::getFindexS3D(size, guide, i, j, k);
            s = bx[m];
            if ( IS_FLUID( s ) ) {
              s = onBit( s, ACTIVE_BIT );
              c++;
            }
            else {
              s = offBit( s, ACTIVE_BIT );
            }
            bx[m] = s;
          }
        }
      }
      break;
      
    case SOLID_CONDUCTION:
      for (int k=1; k<=kx; k++) {
        for (int j=1; j<=jx; j++) {
          for (int i=1; i<=ix; i++) {
            m = FBUtility::getFindexS3D(size, guide, i, j, k);
            s = bx[m];
            if ( !IS_FLUID( s ) ) {
              s = onBit( s, ACTIVE_BIT );
              c++;
            }
            else {
              s = offBit( s, ACTIVE_BIT );
            }
            bx[m] = s;
          }
        }
      }
      break;
      
    case CONJUGATE_HEAT_TRANSFER:
      for (int k=1; k<=kx; k++) {
        for (int j=1; j<=jx; j++) {
          for (int i=1; i<=ix; i++) {
            m = FBUtility::getFindexS3D(size, guide, i, j, k);
            bx[m] = onBit( bx[m], ACTIVE_BIT );
            c++;
          }
        }
      }
      break;
  }
  
  Lcell = c;
  
  
  if ( pn.numProc > 1 ) {
    unsigned long c_tmp = c;
    MPI_Allreduce(&c_tmp, &c, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  }
  
  Gcell = c;
}

/**
 @fn void VoxInfo::encAmask_SymtrcBC(int face, unsigned* bh2)
 @brief 外部境界に接するセルに対称境界面の断熱マスクをセットする
 @param face 外部境界面番号
 @param bh2 BCindex H2
 */
void VoxInfo::encAmask_SymtrcBC(int face, unsigned* bh2)
{
  int register i, j, k;
  unsigned register m;
  
  int ix = (int)size[0];
  int jx = (int)size[1];
  int kx = (int)size[2];
  int gd = (int)guide;
  
  switch (face) {
    case X_MINUS:
      if( pn.nID[face] < 0 ){
        i=1;
        for (k=1; k<=kx; k++) {
          for (j=1; j<=jx; j++) {
            m = FBUtility::getFindexS3D(size, guide, i, j, k); // 最外層のID
            bh2[m] = offBit(bh2[m], ADIABATIC_W); // 断熱ビットを0にする
          }
        }
      }
      break;
      
    case X_PLUS:
      if( pn.nID[face] < 0 ){
        i=ix;
        for (k=1; k<=kx; k++) {
          for (j=1; j<=jx; j++) {
            m = FBUtility::getFindexS3D(size, guide, i, j, k);
            bh2[m] = offBit(bh2[m], ADIABATIC_E);
          }
        }
      }
      break;
      
    case Y_MINUS:
      if( pn.nID[face] < 0 ){
        j=1;
        for (k=1; k<=kx; k++) {
          for (i=1; i<=ix; i++) {
            m = FBUtility::getFindexS3D(size, guide, i, j, k);
            bh2[m] = offBit(bh2[m], ADIABATIC_S);
          }
        }
      }
      break;
      
    case Y_PLUS:
      if( pn.nID[face] < 0 ){
        j=jx;
        for (k=1; k<=kx; k++) {
          for (i=1; i<=ix; i++) {
            m = FBUtility::getFindexS3D(size, guide, i, j, k);
            bh2[m] = offBit(bh2[m], ADIABATIC_N);
          }
        }
      }
      break;
      
    case Z_MINUS:
      if( pn.nID[face] < 0 ){
        k=1;
        for (j=1; j<=jx; j++) {
          for (i=1; i<=ix; i++) {
            m = FBUtility::getFindexS3D(size, guide, i, j, k);
            bh2[m] = offBit(bh2[m], ADIABATIC_B);
          }
        }
      }
      break;
      
    case Z_PLUS:
      if( pn.nID[face] < 0 ){
        k=kx;
        for (j=1; j<=jx; j++) {
          for (i=1; i<=ix; i++) {
            m = FBUtility::getFindexS3D(size, guide, i, j, k);
            bh2[m] = offBit(bh2[m], ADIABATIC_T);
          }
        }
      }
      break;
  }
}

/**
 @fn void VoxInfo::encHbit(unsigned* bh1, unsigned* bh2)
 @brief セルの各面を調べ，境界条件が設定されていれば，ビットをON
 @param bh1 BCindex H1
 @param bh2 BCindex H2
 @note 対角要素の係数をエンコードする
 */
void VoxInfo::encHbit(unsigned* bh1, unsigned* bh2)
{
  int i,j,k;
  unsigned register m, s1, s2;
  unsigned s_e, s_w, s_n, s_s, s_t, s_b, ss;
  unsigned a_e, a_w, a_n, a_s, a_t, a_b;
  
  int ix = (int)size[0];
  int jx = (int)size[1];
  int kx = (int)size[2];
  int gd = (int)guide;
  
  // 初期化
  for (k=0; k<=kx+1; k++) {
    for (j=0; j<=jx+1; j++) {
      for (i=0; i<=ix+1; i++) {
        m = FBUtility::getFindexS3D(size, guide, i  , j  , k  );
        bh2[m] |= ( 0x3f << GMA_W ); // 6bitまとめて(1)で初期化
      }
    }
  }
  
  // 境界条件ビットの設定
  for (k=1; k<=kx; k++) {
    for (j=1; j<=jx; j++) {
      for (i=1; i<=ix; i++) {
        m = FBUtility::getFindexS3D(size, guide, i  , j  , k  );
        s1 = bh1[m];
        s2 = bh2[m];
        
        if ( GET_FACE_BC(s1, BC_FACE_W) != 0 ) s2 = offBit( s2, GMA_W );
        if ( GET_FACE_BC(s1, BC_FACE_E) != 0 ) s2 = offBit( s2, GMA_E );
        if ( GET_FACE_BC(s1, BC_FACE_S) != 0 ) s2 = offBit( s2, GMA_S );
        if ( GET_FACE_BC(s1, BC_FACE_N) != 0 ) s2 = offBit( s2, GMA_N );
        if ( GET_FACE_BC(s1, BC_FACE_B) != 0 ) s2 = offBit( s2, GMA_B );
        if ( GET_FACE_BC(s1, BC_FACE_T) != 0 ) s2 = offBit( s2, GMA_T );
        
        bh2[m] = s2;
      }
    }
  }
  
  // 対角要素の係数のチェックとエンコード
  for (k=1; k<=kx; k++) {
    for (j=1; j<=jx; j++) {
      for (i=1; i<=ix; i++) {
        m = FBUtility::getFindexS3D(size, guide, i, j, k);
        s2= bh2[m];
        
        s_w = BIT_SHIFT(s2, GMA_W); // 対角要素の係数：0 or 1
        s_e = BIT_SHIFT(s2, GMA_E); 
        s_s = BIT_SHIFT(s2, GMA_S);
        s_n = BIT_SHIFT(s2, GMA_N);
        s_b = BIT_SHIFT(s2, GMA_B);
        s_t = BIT_SHIFT(s2, GMA_T);
        
        a_w = BIT_SHIFT(s2, ADIABATIC_W);
        a_e = BIT_SHIFT(s2, ADIABATIC_E);
        a_s = BIT_SHIFT(s2, ADIABATIC_S);
        a_n = BIT_SHIFT(s2, ADIABATIC_N);
        a_b = BIT_SHIFT(s2, ADIABATIC_B);
        a_t = BIT_SHIFT(s2, ADIABATIC_T);
        
        ss = s_w * a_w
        + s_e * a_e
        + s_s * a_s
        + s_n * a_n
        + s_b * a_b
        + s_t * a_t;
        bh2[m] = s2 | (ss<<H_DIAG);
        
        if ( ss == 0 ) { 
          Hostonly_ printf("\tDiagonal element is zero at (%d,%d,%d) : (Gamma:wesnbt)[%1d %1d %1d %1d %1d %1d] (A:wesnbt)[%1d %1d %1d %1d %1d %1d]\n",
                           i,j,k,
                           s_w, s_e, s_s, s_n, s_b, s_t, 
                           a_w, a_e, a_s, a_n, a_b, a_t);
        }
        
      }
    }
  }
  
  // ゼロ割防止のためのダミー係数 >> 全領域
  unsigned nx = (size[0]+2*guide) * (size[1]+2*guide) * (size[2]+2*guide);
  for (m=0; m<nx; m++) {
    s2 = bh2[m];
    if ( ((s2>>H_DIAG) & 0x7) == 0 ) { // 0x7 = 3 bit
      bh2[m] = s2 | (0x1<<H_DIAG);
    }
  }
}


/**
 @fn unsigned VoxInfo::encodeOrder(unsigned order, unsigned id, int* mid, unsigned* bx)
 @brief bx[]へCompoListのエントリをエンコードする
 @param order エンコードするエントリ
 @param id サーチ対象ID
 @param mid セルID配列
 @param bx BCindex ID/H2
 @retval エンコードした個数
 @note
 - mid[]のセルIDが指定されたidならば，bx[]に対してCompoListのエントリをエンコードする
 */
unsigned long VoxInfo::encodeOrder(const unsigned order, const unsigned id, const int* mid, unsigned* bx)
{
  int idd;
  unsigned register m;
  unsigned long g=0;
  
  int ix = (int)size[0];
  int jx = (int)size[1];
  int kx = (int)size[2];
  
  idd = (int)id;
  
  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        m = FBUtility::getFindexS3D(size, guide, i, j, k);
        if ( mid[m] == idd )  {
          bx[m] |= order; // bx[m]の下位6bitにエントリをエンコード  >> ParseBC:sertControlVars()でビット幅をチェック
          g++;
        }
      }
    }
  }
  
  
  if ( pn.numProc > 1 ) {
    unsigned long tmp = g;
    MPI_Allreduce(&tmp, &g, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  }

  return g;
}


/**
 @fn unsigned long VoxInfo::encQfaceHT_S(unsigned order, unsigned id, int* mid, unsigned* bcd, unsigned* bh1, unsigned* bh2, int deface)
 @brief TRANSFER_S/SF/SNに必要な情報をエンコードする
 @retval エンコードしたセル数
 @param order CompoListのエントリ
 @param id 対象ID
 @param mid セルID配列
 @param bcd BCindex ID
 @param bh1 BCindex H1
 @param bh2 BCindex H2
 @param deface セルフェイスを指定するための参照ID
 @note
 - この境界条件処理は，固体セルに接する流体セルの熱伝達面に対して適用される
 - 対象面の断熱ビットを非断熱(1)に戻しておく
 */
unsigned long VoxInfo::encQfaceHT_S(const unsigned order, 
                                    const unsigned id, 
                                    const int* mid, 
                                    unsigned* bcd, 
                                    unsigned* bh1, 
                                    unsigned* bh2, 
                                    const int deface)
{
  unsigned long g=0;
  unsigned m_p, m_e, m_w, m_n, m_s, m_t, m_b;
  int      c_p, c_e, c_w, c_n, c_s, c_t, c_b;
  unsigned s_e, s_w, s_n, s_s, s_t, s_b;
  unsigned register s1, s2, d;
  
  int ix = (int)size[0];
  int jx = (int)size[1];
  int kx = (int)size[2];
  int gd = (int)guide;
  
  int idd = (int)id;
  
  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        m_p = FBUtility::getFindexS3D(size, guide, i  , j  , k  );
        c_p = mid[m_p];
        
        d  = bcd[m_p];
        s1 = bh1[m_p];
        s2 = bh2[m_p];
        
        if ( (c_p == deface) && IS_FLUID(s2) ) { // 流体セルに対してのみ適用
          
          m_e = FBUtility::getFindexS3D(size, guide, i+1, j  , k  );
          m_w = FBUtility::getFindexS3D(size, guide, i-1, j  , k  );
          m_n = FBUtility::getFindexS3D(size, guide, i  , j+1, k  );
          m_s = FBUtility::getFindexS3D(size, guide, i  , j-1, k  );
          m_t = FBUtility::getFindexS3D(size, guide, i  , j  , k+1);
          m_b = FBUtility::getFindexS3D(size, guide, i  , j  , k-1);
          
          c_e = mid[m_e];
          c_w = mid[m_w];
          c_n = mid[m_n];
          c_s = mid[m_s];
          c_t = mid[m_t];
          c_b = mid[m_b];
          
          s_e = bh2[m_e];
          s_w = bh2[m_w];
          s_n = bh2[m_n];
          s_s = bh2[m_s];
          s_t = bh2[m_t];
          s_b = bh2[m_b];
          
          // X-
          if ( c_w == idd ) { // 指定IDで挟まれる面の候補
            if ( !IS_FLUID(s_w) ) {        // 隣接セルが固体であること
              d |= order; // エントリをエンコード
              s1 |= (order << BC_FACE_W);  // エントリをエンコード
              s2 = onBit(s2, ADIABATIC_W); // 断熱ビットを非断熱
              g++;
            }
          }
          
          // X+
          if ( c_e == idd ) {
            if ( !IS_FLUID(s_e) ) {
              d |= order;
              s1 |= (order << BC_FACE_E);
              s2 = onBit(s2, ADIABATIC_E);
              g++;
            }
          }
          
          // Y-
          if ( c_s == idd ) {
            if ( !IS_FLUID(s_s) ) {
              d |= order;
              s1 |= (order << BC_FACE_S);
              s2 = onBit(s2, ADIABATIC_S);
              g++;
            }
          }
          
          // Y+
          if ( c_n == idd ) {
            if ( !IS_FLUID(s_n) ) {
              d |= order;
              s1 |= (order << BC_FACE_N);
              s2 = onBit(s2, ADIABATIC_N);
              g++;
            }
          }
          
          // Z-
          if ( c_b == idd ) {
            if ( !IS_FLUID(s_b) ) {
              d |= order;
              s1 |= (order << BC_FACE_B);
              s2 = onBit(s2, ADIABATIC_B);
              g++;
            }
          }
          
          // Z+
          if ( c_t == idd ) {
            if ( !IS_FLUID(s_t) ) {
              d |= order;
              s1 |= (order << BC_FACE_T);
              s2 = onBit(s2, ADIABATIC_T);
              g++;
            }
          }
          
          bcd[m_p] = d;
          bh1[m_p] = s1;
          bh2[m_p] = s2;
        }
        
      }
    }
  }

  
  if ( pn.numProc > 1 ) {
    unsigned long tmp = g;
    MPI_Allreduce(&tmp, &g, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  }
  
  return g;
}

/**
 @fn unsigned long VoxInfo::encQfaceHT_B(unsigned order, unsigned id, int* mid, unsigned* bcd, unsigned* bh1, unsigned* bh2, int deface)
 @brief TRANSFER_Bに必要な情報をエンコードする
 @retval エンコードしたセル数
 @param order CompoListのエントリ
 @param id 対象ID
 @param mid セルID配列
 @param bcd BCindex ID
 @param bh1 BCindex H1
 @param bh2 BCindex H2
 @param deface セルフェイスを指定するための参照ID
 @note
 - この境界条件処理は，流体セルに接する固体セルの熱伝達面に対して適用される
 - 対象面の断熱ビットを非断熱(1)に戻しておく
 */
unsigned long VoxInfo::encQfaceHT_B(const unsigned order, 
                                    const unsigned id, 
                                    const int* mid, 
                                    unsigned* bcd, 
                                    unsigned* bh1, 
                                    unsigned* bh2, 
                                    const int deface)
{
  unsigned long g=0;
  unsigned m_p, m_e, m_w, m_n, m_s, m_t, m_b;
  int      c_p, c_e, c_w, c_n, c_s, c_t, c_b;
  unsigned s_e, s_w, s_n, s_s, s_t, s_b;
  unsigned register s1, s2, d;
  
  int ix = (int)size[0];
  int jx = (int)size[1];
  int kx = (int)size[2];
  int gd = (int)guide;
  
  int idd = (int)id;
  
  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        m_p = FBUtility::getFindexS3D(size, guide, i  , j  , k  );
        c_p = mid[m_p];
        
        d  = bcd[m_p];
        s1 = bh1[m_p];
        s2 = bh2[m_p];
        
        if ( (c_p == idd) && !IS_FLUID(s2) ) { // 固体セルに対してのみ適用
          m_e = FBUtility::getFindexS3D(size, guide, i+1, j  , k  );
          m_w = FBUtility::getFindexS3D(size, guide, i-1, j  , k  );
          m_n = FBUtility::getFindexS3D(size, guide, i  , j+1, k  );
          m_s = FBUtility::getFindexS3D(size, guide, i  , j-1, k  );
          m_t = FBUtility::getFindexS3D(size, guide, i  , j  , k+1);
          m_b = FBUtility::getFindexS3D(size, guide, i  , j  , k-1);
          
          c_e = mid[m_e];
          c_w = mid[m_w];
          c_n = mid[m_n];
          c_s = mid[m_s];
          c_t = mid[m_t];
          c_b = mid[m_b];
          
          s_e = bh2[m_e];
          s_w = bh2[m_w];
          s_n = bh2[m_n];
          s_s = bh2[m_s];
          s_t = bh2[m_t];
          s_b = bh2[m_b];
          
          // X-
          if ( c_w == deface ) { // 指定IDで挟まれる面の候補
            if ( IS_FLUID(s_w) ) {        // 隣接セルが流体であること
              d |= order; // エントリをエンコード
              s1 |= (order << BC_FACE_W);  // エントリをエンコード
              s2 = onBit(s2, ADIABATIC_W); // 断熱ビットを非断熱
              g++;
            }
          }
          
          // X+
          if ( c_e == deface ) {
            if ( IS_FLUID(s_e) ) {
              d |= order;
              s1 |= (order << BC_FACE_E);
              s2 = onBit(s2, ADIABATIC_E);
              g++;
            }
          }
          
          // Y-
          if ( c_s == deface ) {
            if ( IS_FLUID(s_s) ) {
              d |= order;
              s1 |= (order << BC_FACE_S);
              s2 = onBit(s2, ADIABATIC_S);
              g++;
            }
          }
          
          // Y+
          if ( c_n == deface ) {
            if ( IS_FLUID(s_n) ) {
              d |= order;
              s1 |= (order << BC_FACE_N);
              s2 = onBit(s2, ADIABATIC_N);
              g++;
            }
          }
          
          // Z-
          if ( c_b == deface ) {
            if ( IS_FLUID(s_b) ) {
              d |= order;
              s1 |= (order << BC_FACE_B);
              s2 = onBit(s2, ADIABATIC_B);
              g++;
            }
          }
          
          // Z+
          if ( c_t == deface ) {
            if ( IS_FLUID(s_t) ) {
              d |= order;
              s1 |= (order << BC_FACE_T);
              s2 = onBit(s2, ADIABATIC_T);
              g++;
            }
          }
          
          bcd[m_p] = d;
          bh1[m_p] = s1;
          bh2[m_p] = s2;
        }
        
      }
    }
  }

  
  if ( pn.numProc > 1 ) {
    unsigned long tmp = g;
    MPI_Allreduce(&tmp, &g, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  }
  
  return g;
}

/**
 @fn unsigned long VoxInfo::encQfaceISO_SF(unsigned order, unsigned id, int* mid, unsigned* bcd, unsigned* bh1, unsigned* bh2, int deface)
 @brief SF界面のISOTHERMAL処理
 @retval エンコードしたセル数
 @param order CompoListのエントリ
 @param id 対象ID
 @param mid セルID配列
 @param bcd BCindex ID
 @param bh1 BCindex H1
 @param bh2 BCindex H2
 @param deface セルフェイスを指定するための参照ID
 @note
 - この境界条件処理は，流体と固体で挟まれる面に対して適用される
 - defaceで指定されるセルは等温壁をもつ計算セルで，等温面はidで指定されるセルで挟まれる
 - 対象面の断熱ビットを非断熱(1)に戻しておく
 */
unsigned long VoxInfo::encQfaceISO_SF(const unsigned order, 
                                      const unsigned id, 
                                      const int* mid, 
                                      unsigned* bcd, 
                                      unsigned* bh1, 
                                      unsigned* bh2, 
                                      const int deface)
{
  unsigned long g=0;
  unsigned m_p, m_e, m_w, m_n, m_s, m_t, m_b;
  int      c_p, c_e, c_w, c_n, c_s, c_t, c_b;
  unsigned s_e, s_w, s_n, s_s, s_t, s_b;
  unsigned register s1, s2, d;
  
  int ix = (int)size[0];
  int jx = (int)size[1];
  int kx = (int)size[2];
  int gd = (int)guide;
  
  int idd = (int)id;
  
  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        m_p = FBUtility::getFindexS3D(size, guide, i  , j  , k  );
        c_p = mid[m_p];
        
        d  = bcd[m_p];
        s1 = bh1[m_p];
        s2 = bh2[m_p];
        
        if ( (c_p == deface) && IS_FLUID(s2) ) { // テストセルは流体でdefaceID
          m_e = FBUtility::getFindexS3D(size, guide, i+1, j  , k  );
          m_w = FBUtility::getFindexS3D(size, guide, i-1, j  , k  );
          m_n = FBUtility::getFindexS3D(size, guide, i  , j+1, k  );
          m_s = FBUtility::getFindexS3D(size, guide, i  , j-1, k  );
          m_t = FBUtility::getFindexS3D(size, guide, i  , j  , k+1);
          m_b = FBUtility::getFindexS3D(size, guide, i  , j  , k-1);
          
          s_e = bh2[m_e];
          s_w = bh2[m_w];
          s_n = bh2[m_n];
          s_s = bh2[m_s];
          s_t = bh2[m_t];
          s_b = bh2[m_b];
          
          c_e = mid[m_e];
          c_w = mid[m_w];
          c_n = mid[m_n];
          c_s = mid[m_s];
          c_t = mid[m_t];
          c_b = mid[m_b];
          
          // X-
          if ( c_w == idd ) { // 指定IDで挟まれる面
            if ( !IS_FLUID(s_w) ) {        // 隣接セルが固体であること
              d |= order; // エントリをエンコード
              s1 |= (order << BC_FACE_W);  // エントリをエンコード
              s2 = onBit(s2, ADIABATIC_W); // 断熱ビットを非断熱
              g++;
            }
          }
          
          // X+
          if ( c_e == idd ) {
            if ( !IS_FLUID(s_e) ) {
              d |= order;
              s1 |= (order << BC_FACE_E);
              s2 = onBit(s2, ADIABATIC_E);
              g++;
            }
          }
          
          // Y-
          if ( c_s == idd ) {
            if ( !IS_FLUID(s_s) ) {
              d |= order;
              s1 |= (order << BC_FACE_S);
              s2 = onBit(s2, ADIABATIC_S);
              g++;
            }
          }
          
          // Y+
          if ( c_n == idd ) {
            if ( !IS_FLUID(s_n) ) {
              d |= order;
              s1 |= (order << BC_FACE_N);
              s2 = onBit(s2, ADIABATIC_N);
              g++;
            }
          }
          
          // Z-
          if ( c_b == idd ) {
            if ( !IS_FLUID(s_b) ) {
              d |= order;
              s1 |= (order << BC_FACE_B);
              s2 = onBit(s2, ADIABATIC_B);
              g++;
            }
          }
          
          // Z+
          if ( c_t == idd ) {
            if ( !IS_FLUID(s_t) ) {
              d |= order;
              s1 |= (order << BC_FACE_T);
              s2 = onBit(s2, ADIABATIC_T);
              g++;
            }
          }
          
          bcd[m_p] = d;
          bh1[m_p] = s1;
          bh2[m_p] = s2;
        }
        
      }
    }
  }

  
  if ( pn.numProc > 1 ) {
    unsigned long tmp = g;
    MPI_Allreduce(&tmp, &g, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  }
  
  return g;
}

/**
 @fn unsigned long VoxInfo::encQfaceISO_SS(unsigned order, unsigned id, int* mid, unsigned* bcd, unsigned* bh1, unsigned* bh2, int deface)
 @brief SS界面のISOTHERMAL処理
 @retval エンコードしたセル数
 @param order CompoListのエントリ
 @param id 対象ID
 @param mid セルID配列
 @param bcd BCindex ID
 @param bh1 BCindex H1
 @param bh2 BCindex H2
 @param deface セルフェイスを指定するための参照ID
 @note
 - この境界条件処理は，固体と固体で挟まれる面に対して適用される
 -　idで指定されるセルは等温壁をもつ計算セルで，等温面はidで指定されるセルで挟まれる
 - 対象面の断熱ビットを非断熱(1)に戻しておく
 */
unsigned long VoxInfo::encQfaceISO_SS(const unsigned order, 
                                      const unsigned id, 
                                      const int* mid, 
                                      unsigned* bcd, 
                                      unsigned* bh1, 
                                      unsigned* bh2, 
                                      const int deface)
{
  unsigned long g=0;
  unsigned m_p, m_e, m_w, m_n, m_s, m_t, m_b;
  int      c_p, c_e, c_w, c_n, c_s, c_t, c_b;
  unsigned s_e, s_w, s_n, s_s, s_t, s_b;
  unsigned register s1, s2, d;
  
  int ix = (int)size[0];
  int jx = (int)size[1];
  int kx = (int)size[2];
  int gd = (int)guide;
  
  int idd = (int)id;
  
  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        m_p = FBUtility::getFindexS3D(size, guide, i  , j  , k  );
        c_p = mid[m_p];
        
        d  = bcd[m_p];
        s1 = bh1[m_p];
        s2 = bh2[m_p];
        
        if ( (c_p == idd) && !IS_FLUID(s2) ) { // テストセルは固体でキーID
          m_e = FBUtility::getFindexS3D(size, guide, i+1, j  , k  );
          m_w = FBUtility::getFindexS3D(size, guide, i-1, j  , k  );
          m_n = FBUtility::getFindexS3D(size, guide, i  , j+1, k  );
          m_s = FBUtility::getFindexS3D(size, guide, i  , j-1, k  );
          m_t = FBUtility::getFindexS3D(size, guide, i  , j  , k+1);
          m_b = FBUtility::getFindexS3D(size, guide, i  , j  , k-1);
          
          s_e = bh2[m_e];
          s_w = bh2[m_w];
          s_n = bh2[m_n];
          s_s = bh2[m_s];
          s_t = bh2[m_t];
          s_b = bh2[m_b];
          
          c_e = mid[m_e];
          c_w = mid[m_w];
          c_n = mid[m_n];
          c_s = mid[m_s];
          c_t = mid[m_t];
          c_b = mid[m_b];
          
          // X-
          if ( c_w == deface ) { // 指定IDで挟まれる面
            if ( !IS_FLUID(s_w) ) {        // 隣接セルが固体であること
              d |= order; // エントリをエンコード
              s1 |= (order << BC_FACE_W);  // エントリをエンコード
              s2 = onBit(s2, ADIABATIC_W); // 断熱ビットを非断熱
              g++;
            }
          }
          
          // X+
          if ( c_e == deface ) {
            if ( !IS_FLUID(s_e) ) {
              d |= order;
              s1 |= (order << BC_FACE_E);
              s2 = onBit(s2, ADIABATIC_E);
              g++;
            }
          }
          
          // Y-
          if ( c_s == deface ) {
            if ( !IS_FLUID(s_s) ) {
              d |= order;
              s1 |= (order << BC_FACE_S);
              s2 = onBit(s2, ADIABATIC_S);
              g++;
            }
          }
          
          // Y+
          if ( c_n == deface ) {
            if ( !IS_FLUID(s_n) ) {
              d |= order;
              s1 |= (order << BC_FACE_N);
              s2 = onBit(s2, ADIABATIC_N);
              g++;
            }
          }
          
          // Z-
          if ( c_b == deface ) {
            if ( !IS_FLUID(s_b) ) {
              d |= order;
              s1 |= (order << BC_FACE_B);
              s2 = onBit(s2, ADIABATIC_B);
              g++;
            }
          }
          
          // Z+
          if ( c_t == deface ) {
            if ( !IS_FLUID(s_t) ) {
              d |= order;
              s1 |= (order << BC_FACE_T);
              s2 = onBit(s2, ADIABATIC_T);
              g++;
            }
          }
          
          bcd[m_p] = d;
          bh1[m_p] = s1;
          bh2[m_p] = s2;
        }
        
      }
    }
  }
  
  if ( pn.numProc > 1 ) {
    unsigned long tmp = g;
    MPI_Allreduce(&tmp, &g, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  }
  
  return g;
}

/**
 @fn unsigned long VoxInfo::encQface(unsigned order, unsigned id, int* mid, unsigned* bcd, unsigned* bh1, unsigned* bh2, int deface, bool flag)
 @brief 熱境界条件のBCエントリをエンコードする
 @retval カウントしたセル数
 @param order CompoListのエントリ
 @param id 対象セルID
 @param mid セルID配列
 @param bcd BCindex ID
 @param bh1 BCindex H1
 @param bh2 BCindex H2
 @param deface 界面指定セルID
 @param flag 断熱ビットのon(true)/off(false)
 @note
 - idとdefaceで挟まれる面に対して適用
 - 対象面の断熱ビットはflagで判断
 */
unsigned long VoxInfo::encQface(const unsigned order, 
                                const unsigned id, 
                                const int* mid, 
                                unsigned* bcd, 
                                unsigned* bh1, 
                                unsigned* bh2, 
                                const int deface, 
                                const bool flag)
{
  unsigned long g=0;
  unsigned m_p, m_e, m_w, m_n, m_s, m_t, m_b;
  int      c_p, c_e, c_w, c_n, c_s, c_t, c_b;
  unsigned register s1, s2, d;
  
  int ix = (int)size[0];
  int jx = (int)size[1];
  int kx = (int)size[2];
  int gd = (int)guide;
  
  int idd = (int)id;
  
  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        m_p = FBUtility::getFindexS3D(size, guide, i  , j  , k  );
        c_p = mid[m_p];
        
        if ( c_p == deface) {
          m_e = FBUtility::getFindexS3D(size, guide, i+1, j  , k  );
          m_w = FBUtility::getFindexS3D(size, guide, i-1, j  , k  );
          m_n = FBUtility::getFindexS3D(size, guide, i  , j+1, k  );
          m_s = FBUtility::getFindexS3D(size, guide, i  , j-1, k  );
          m_t = FBUtility::getFindexS3D(size, guide, i  , j  , k+1);
          m_b = FBUtility::getFindexS3D(size, guide, i  , j  , k-1);
          
          c_e = mid[m_e];
          c_w = mid[m_w];
          c_n = mid[m_n];
          c_s = mid[m_s];
          c_t = mid[m_t];
          c_b = mid[m_b];
          
          d  = bcd[m_p];
          s1 = bh1[m_p];
          s2 = bh2[m_p];
          
          // X-
          if ( c_w == idd ) {
            d |= order; // エントリをエンコード
            s1 |= (order << BC_FACE_W);  // エントリをエンコード
            s2 = ( flag == true ) ? onBit(s2, ADIABATIC_W) : offBit( s2, ADIABATIC_W );
            g++;
          }
          
          // X+
          if ( c_e == idd ) {
            d |= order;
            s1 |= (order << BC_FACE_E);
            s2 = ( flag == true ) ? onBit(s2, ADIABATIC_E) : offBit( s2, ADIABATIC_E );
            g++;
          }
          
          // Y-
          if ( c_s == idd ) {
            d |= order;
            s1 |= (order << BC_FACE_S);
            s2 = ( flag == true ) ? onBit(s2, ADIABATIC_S) : offBit( s2, ADIABATIC_S );
            g++;
          }
          
          // Y+
          if ( c_n == idd ) {
            d |= order;
            s1 |= (order << BC_FACE_N);
            s2 = ( flag == true ) ? onBit(s2, ADIABATIC_N) : offBit( s2, ADIABATIC_N );
            g++;
          }
          
          // Z-
          if ( c_b == idd ) {
            d |= order;
            s1 |= (order << BC_FACE_B);
            s2 = ( flag == true ) ? onBit(s2, ADIABATIC_B) : offBit( s2, ADIABATIC_B );
            g++;
          }
          
          // Z+
          if ( c_t == idd ) {
            d |= order;
            s1 |= (order << BC_FACE_T);
            s2 = ( flag == true ) ? onBit(s2, ADIABATIC_T) : offBit( s2, ADIABATIC_T );
            g++;            
          }
          
          bcd[m_p] = d;
          bh1[m_p] = s1;
          bh2[m_p] = s2;
        }
        
      }
    }
  }
  
  if ( pn.numProc > 1 ) {
    unsigned long tmp = g;
    MPI_Allreduce(&tmp, &g, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  }
  
  return g;
}


/**
 @fn void VoxInfo::encQfaceSVO(unsigned order, unsigned id, int* mid, unsigned* bcd, unsigned* bh1, unsigned* bh2, int deface)
 @brief 速度指定，流出型の境界条件の計算に必要な情報をエンコードする
 @param order cmp[]のエントリ番号
 @param id  セルID
 @param mid ボクセル配列
 @param bcd BCindex ID
 @param bh1 BCindex H1
 @param bh2 BCindex H2
 @param deface 面を指定するid
 @note 
 - encVbit_IBC()でセル数はカウント済み
 - 流体かつ指定IDであり，defaceが設定されたセルと挟まれる面に対して適用
 - 対象面の断熱ビットを非断熱にする
 - encVbit_IBC()では，最外層を除外しているが，熱では除外しない
 */
void VoxInfo::encQfaceSVO(unsigned order, unsigned id, int* mid, unsigned* bcd, unsigned* bh1, unsigned* bh2, int deface)
{
  int i,j,k, idd;
  unsigned m_p, m_e, m_w, m_n, m_s, m_t, m_b;
  int      c_p, c_e, c_w, c_n, c_s, c_t, c_b;
  unsigned register s1, s2, d;
  
  int ix = (int)size[0];
  int jx = (int)size[1];
  int kx = (int)size[2];
  int gd = (int)guide;
  
  idd = (int)id;
  
  for (k=1; k<=kx; k++) {
    for (j=1; j<=jx; j++) {
      for (i=1; i<=ix; i++) {
        m_p = FBUtility::getFindexS3D(size, guide, i  , j  , k  );
        c_p = mid[m_p];
        
        d  = bcd[m_p];
        s1 = bh1[m_p];
        s2 = bh2[m_p];
        
        if ( (c_p == deface) && IS_FLUID(s1) ) { // defaceセルかつ流体の時に，以下を評価
          m_e = FBUtility::getFindexS3D(size, guide, i+1, j  , k  );
          m_w = FBUtility::getFindexS3D(size, guide, i-1, j  , k  );
          m_n = FBUtility::getFindexS3D(size, guide, i  , j+1, k  );
          m_s = FBUtility::getFindexS3D(size, guide, i  , j-1, k  );
          m_t = FBUtility::getFindexS3D(size, guide, i  , j  , k+1);
          m_b = FBUtility::getFindexS3D(size, guide, i  , j  , k-1);
          
          c_e = mid[m_e];
          c_w = mid[m_w];
          c_n = mid[m_n];
          c_s = mid[m_s];
          c_t = mid[m_t];
          c_b = mid[m_b];
          
          // X-
          if ( c_w == idd ) { // m_wセルの属性チェックはしていない
            d |= order; // エントリをエンコード
            s1 |= (order << BC_FACE_W);  // エントリをエンコード
            s2 = onBit(s2, ADIABATIC_W); // 断熱ビットを非断熱(1)
          }
          
          // X+
          if ( c_e == idd ) {
            d |= order;
            s1 |= (order << BC_FACE_E);
            s2 = onBit(s2, ADIABATIC_E);
          }
          
          // Y-
          if ( c_s == idd ) {
            d |= order;
            s1 |= (order << BC_FACE_S);
            s2 = onBit(s2, ADIABATIC_S);
          }
          
          // Y+
          if ( c_n == idd ) {
            d |= order;
            s1 |= (order << BC_FACE_N);
            s2 = onBit(s2, ADIABATIC_N);
          }
          
          // Z-
          if ( c_b == idd ) {
            d |= order;
            s1 |= (order << BC_FACE_B);
            s2 = onBit(s2, ADIABATIC_B);
          }
          
          // Z+
          if ( c_t == idd ) {
            d |= order;
            s1 |= (order << BC_FACE_T);
            s2 = onBit(s2, ADIABATIC_T);
          }
          
          bcd[m_p] = d;
          bh1[m_p] = s1;
          bh2[m_p] = s2;
        }
      } // i-loop
    }
  }
  
  // iddセルのSTATE_BITをFLUIDに変更
  for (k=1; k<=kx; k++) {
    for (j=1; j<=jx; j++) {
      for (i=1; i<=ix; i++) {
        m_p = FBUtility::getFindexS3D(size, guide, i  , j  , k  );
        c_p = mid[m_p];
        
        s1  = bh1[m_p];
        
        if ( (c_p == deface) && IS_FLUID(s1) ) { // 指定セルかつ流体の時に，以下を評価
          m_e = FBUtility::getFindexS3D(size, guide, i+1, j  , k  );
          m_w = FBUtility::getFindexS3D(size, guide, i-1, j  , k  );
          m_n = FBUtility::getFindexS3D(size, guide, i  , j+1, k  );
          m_s = FBUtility::getFindexS3D(size, guide, i  , j-1, k  );
          m_t = FBUtility::getFindexS3D(size, guide, i  , j  , k+1);
          m_b = FBUtility::getFindexS3D(size, guide, i  , j  , k-1);
          
          c_e = mid[m_e];
          c_w = mid[m_w];
          c_n = mid[m_n];
          c_s = mid[m_s];
          c_t = mid[m_t];
          c_b = mid[m_b];
          
          // X-
          if ( c_w == idd ) {
            bh1[m_w] = onBit( bh1[m_w], STATE_BIT );
          }
          
          // X+
          if ( c_e == idd ) {
            bh1[m_e] = onBit( bh1[m_e], STATE_BIT );
          }
          
          // Y-
          if ( c_s == idd ) {
            bh1[m_s] = onBit( bh1[m_s], STATE_BIT );
          }
          
          // Y+
          if ( c_n == idd ) {
            bh1[m_n] = onBit( bh1[m_n], STATE_BIT );
          }
          
          // Z-
          if ( c_b == idd ) {
            bh1[m_b] = onBit( bh1[m_b], STATE_BIT );
          }
          
          // Z+
          if ( c_t == idd ) {
            bh1[m_t] = onBit( bh1[m_t], STATE_BIT );
          }
        }
      } // i-loop
    }
  }
}

/**
 @fn void VoxInfo::encPbit(unsigned* bx)
 @brief ディリクレ条件とノイマン条件の排他性をチェックし，反復行列の非対角要素/対角要素の係数をエンコードする
 @param bx BCindex P
 @note
 - ディリクレ条件とノイマン条件の排他性のチェック
 - 非対角要素と対角要素の係数をエンコードする
 */
void VoxInfo::encPbit(unsigned* bx)
{
  int i, j, k;
  unsigned m, coef, flag, nx;
  unsigned s_e, s_w, s_n, s_s, s_t, s_b, ss;
  unsigned d_e, d_w, d_n, d_s, d_t, d_b;
  unsigned register s;
  bool exclusive;
  
  int ix = (int)size[0];
  int jx = (int)size[1];
  int kx = (int)size[2];
  int gd = (int)guide;
  
  // ディリクレ条件とノイマン条件の排他性のチェック
  exclusive = true;
  for (k=1; k<=kx; k++) {
    for (j=1; j<=jx; j++) {
      for (i=1; i<=ix; i++) {
        m = FBUtility::getFindexS3D(size, guide, i, j, k);
        s = bx[m];
        flag = 0;
        
        // ノイマン条件：値がゼロのとき，BCがセットされている
        s_e = BIT_SHIFT(s, BC_N_E); 
        s_w = BIT_SHIFT(s, BC_N_W);
        s_n = BIT_SHIFT(s, BC_N_N);
        s_s = BIT_SHIFT(s, BC_N_S);
        s_t = BIT_SHIFT(s, BC_N_T);
        s_b = BIT_SHIFT(s, BC_N_B);
        
        // ディリクレ条件：値がゼロのとき，BCがセットされている
        d_e = BIT_SHIFT(s, BC_D_E); 
        d_w = BIT_SHIFT(s, BC_D_W);
        d_n = BIT_SHIFT(s, BC_D_N);
        d_s = BIT_SHIFT(s, BC_D_S);
        d_t = BIT_SHIFT(s, BC_D_T);
        d_b = BIT_SHIFT(s, BC_D_B);
        
        // 非対角要素の係数をエンコード
        if ( (s_e * d_e) == 0 ) s = offBit( s, BC_NDAG_E );
        if ( (s_w * d_w) == 0 ) s = offBit( s, BC_NDAG_W );
        if ( (s_n * d_n) == 0 ) s = offBit( s, BC_NDAG_N );
        if ( (s_s * d_s) == 0 ) s = offBit( s, BC_NDAG_S );
        if ( (s_t * d_t) == 0 ) s = offBit( s, BC_NDAG_T );
        if ( (s_b * d_b) == 0 ) s = offBit( s, BC_NDAG_B );
        
        bx[m] = s;
        
        if ( (s_e==0) && (d_e==0) ) flag++;
        if ( (s_w==0) && (d_w==0) ) flag++;
        if ( (s_n==0) && (d_n==0) ) flag++;
        if ( (s_s==0) && (d_s==0) ) flag++;
        if ( (s_t==0) && (d_t==0) ) flag++;
        if ( (s_b==0) && (d_b==0) ) flag++;
        
        if ( flag != 0) {
          Hostonly_ printf("\tDirichlet and Neumann BC are specified on the same face in cell (%d,%d,%d)\n", i,j,k);
          exclusive = false;
        }
      }
    }
  }
  if ( !exclusive ) Exit(0);
  
  // 対角要素の係数のチェックとエンコード
  for (k=1; k<=kx; k++) {
    for (j=1; j<=jx; j++) {
      for (i=1; i<=ix; i++) {
        m = FBUtility::getFindexS3D(size, guide, i, j, k);
        s = bx[m];
        
        s_e = BIT_SHIFT(s, BC_N_E); // 0-Neumann / 1-normal
        s_w = BIT_SHIFT(s, BC_N_W);
        s_n = BIT_SHIFT(s, BC_N_N);
        s_s = BIT_SHIFT(s, BC_N_S);
        s_t = BIT_SHIFT(s, BC_N_T);
        s_b = BIT_SHIFT(s, BC_N_B);
        
        ss = s_e + s_w + s_n + s_s + s_t + s_b;
        bx[m] = s | (ss<<BC_DIAG);
        
        if ( (ss == 0) && (BIT_IS_SHIFT(s,ACTIVE_BIT)) ) { 
          Hostonly_ printf("\tError : Coefficient of diagonal element is zero at (%d,%d,%d) : (wesnbt)[%1d %1d %1d %1d %1d %1d]\n", i,j,k,
                           IS_FLUID(bx[FBUtility::getFindexS3D(size, guide, i-1, j, k)]),
                           IS_FLUID(bx[FBUtility::getFindexS3D(size, guide, i+1, j, k)]),
                           IS_FLUID(bx[FBUtility::getFindexS3D(size, guide, i, j-1, k)]),
                           IS_FLUID(bx[FBUtility::getFindexS3D(size, guide, i, j+1, k)]),
                           IS_FLUID(bx[FBUtility::getFindexS3D(size, guide, i, j, k-1)]),
                           IS_FLUID(bx[FBUtility::getFindexS3D(size, guide, i, j, k+1)]) );
        }
        
      }
    }
  }
  
  // ゼロ割防止のためのダミー係数 >> 全領域
  nx = (size[0]+2*guide) * (size[1]+2*guide) * (size[2]+2*guide);
  for (m=0; m<nx; m++) {
    s = bx[m];
    if ( ((s>>BC_DIAG) & 0x7) == 0 ) { // 0x7 = 3 bit
      bx[m] = s | (0x1<<BC_DIAG);
    }
  }
}


/**
 @fn unsigned VoxInfo::encPbit_D_IBC(unsigned order, unsigned id, int* mid, unsigned* bcd, unsigned* bcp, int deface)
 @brief 計算領域内部のコンポーネントのDirichletフラグをbcp[]にエンコードする
 @retval エンコードしたセル数
 @param order cmp[]のエントリ番号
 @param id  CellID
 @param mid ボクセル配列
 @param bcd BCindex ID
 @param bcp BCindex P
 @param deface 面を指定するid
 @note
 - 対象セルが流体セルの場合，隣接する面にDirichletフラグをエンコードする
 - 同種のBCは1セルに一つだけ
 */
unsigned long VoxInfo::encPbit_D_IBC(const unsigned order, 
                                     const unsigned id, 
                                     const int* mid, 
                                     unsigned* bcd, 
                                     unsigned* bcp, 
                                     const int deface)
{
  unsigned long g=0;
  unsigned m_e, m_w, m_n, m_s, m_t, m_b;
  unsigned s_e, s_w, s_n, s_s, s_t, s_b;
  int      c_p, c_e, c_w, c_n, c_s, c_t, c_b;
  unsigned register s, d, m;
  
  int ix = (int)size[0];
  int jx = (int)size[1];
  int kx = (int)size[2];
  int gd = (int)guide;
  
  int idd = (int)id;
  
  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        m = FBUtility::getFindexS3D(size, guide, i, j, k);
        d = bcd[m];
        s = bcp[m];
        c_p = mid[m];
        
        if ( IS_FLUID( s ) && (c_p == deface) ) { // 対象セルが流体セル，かつdefaceである
          m_e = FBUtility::getFindexS3D(size, guide, i+1, j  , k  );
          m_w = FBUtility::getFindexS3D(size, guide, i-1, j  , k  );
          m_n = FBUtility::getFindexS3D(size, guide, i  , j+1, k  );
          m_s = FBUtility::getFindexS3D(size, guide, i  , j-1, k  );
          m_t = FBUtility::getFindexS3D(size, guide, i  , j  , k+1);
          m_b = FBUtility::getFindexS3D(size, guide, i  , j  , k-1);
          
          c_e = mid[m_e];
          c_w = mid[m_w];
          c_n = mid[m_n];
          c_s = mid[m_s];
          c_t = mid[m_t];
          c_b = mid[m_b];
          
          s_e = bcp[m_e];
          s_w = bcp[m_w];
          s_n = bcp[m_n];
          s_s = bcp[m_s];
          s_t = bcp[m_t];
          s_b = bcp[m_b];
          
          // X_MINUS
          if ( c_w == idd ) { // Wセルが指定IDかどうかを先にチェック
            if ( !((pn.nID[X_MINUS] < 0) && (i == 1)) ) { // 続いて領域チェック 外部境界面は除外
              d |= order; // dにorderをエンコード
              s = offBit( s, BC_D_W );
              g++;
            }
          }
          
          // X_PLUS
          if ( c_e == idd ) {
            if ( !((pn.nID[X_PLUS] < 0) && (i == ix)) ) {
              d |= order;
              s = offBit( s, BC_D_E );
              g++;
            }
          }
          
          // Y_MINUS
          if ( c_s == idd ) {
            if ( !((pn.nID[Y_MINUS] < 0) && (j == 1)) ) {
              d |= order;
              s = offBit( s, BC_D_S );
              g++;
            }
          }
          
          // Y_PLUS
          if ( c_n == idd ) {
            if ( !((pn.nID[Y_PLUS] < 0) && (j == jx)) ) {
              d |= order;
              s = offBit( s, BC_D_N );
              g++;
            }
          }
          
          // Z_MINUS
          if ( c_b == idd ) {
            if ( !((pn.nID[Z_MINUS] < 0) && (k == 1)) ) {
              d |= order;
              s = offBit( s, BC_D_B );
              g++;
            }
          }
          
          // Z_PLUS
          if ( c_t == idd ) {
            if ( !((pn.nID[Z_PLUS] < 0) && (k == kx)) ) {
              d |= order;
              s = offBit( s, BC_D_T );
              g++;
            }
          }
          
          bcd[m] = d;
          bcp[m] = s;
        }
      }
    }
  }
  
  if ( pn.numProc > 1 ) {
    unsigned long tmp = g;
    MPI_Allreduce(&tmp, &g, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  }
  
  return g;
}


/**
 @fn unsigned long VoxInfo::encPbit_N_Binary(unsigned* bx)
 @brief bcp[]に壁面境界の圧力ノイマン条件のビットフラグと固体に隣接するFセルに方向フラグ，収束判定の有効フラグをエンコードする
 @param[in/out] bx BCindex P
 @retval 固体表面セル数
 @note 
 - 流体セルのうち，固体セルに隣接する面のノイマンフラグをゼロにする．ただし，内部領域のみ．
 - 固体セルに隣接する流体セルに方向フラグを付与する．全内部領域．
 */
unsigned long VoxInfo::encPbit_N_Binary(unsigned* bx)
{
  unsigned m_p, m_e, m_w, m_n, m_s, m_t, m_b;
  unsigned register s;
  
  int ix = (int)size[0];
  int jx = (int)size[1];
  int kx = (int)size[2];
  int gd = (int)guide;
  
  // ノイマンフラグ
  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        m_p = FBUtility::getFindexS3D(size, guide, i, j, k);
        s = bx[m_p];
        
        if ( IS_FLUID( s ) ) {
          m_e = FBUtility::getFindexS3D(size, guide, i+1, j  , k  );
          m_w = FBUtility::getFindexS3D(size, guide, i-1, j  , k  );
          m_n = FBUtility::getFindexS3D(size, guide, i  , j+1, k  );
          m_s = FBUtility::getFindexS3D(size, guide, i  , j-1, k  );
          m_t = FBUtility::getFindexS3D(size, guide, i  , j  , k+1);
          m_b = FBUtility::getFindexS3D(size, guide, i  , j  , k-1);
          
          // 外部境界面は除外 pn.nID[X_MINUS] < 0 のときが外部境界面，0-myrank, 0>other rank
          // X_MINUS
          if ( !((pn.nID[X_MINUS] < 0) && (i == 1)) ) {
            if ( !IS_FLUID(bx[m_w]) ) {
              s = offBit( s, BC_N_W );
            }
          }
          
          // X_PLUS
          if ( !((pn.nID[X_PLUS] < 0) && (i == ix)) ) {
            if ( !IS_FLUID(bx[m_e]) ) {
              s = offBit( s, BC_N_E );
            }
          }
          
          // Y_MINUS
          if ( !((pn.nID[Y_MINUS] < 0) && (j == 1)) ) {
            if ( !IS_FLUID(bx[m_s]) ) {
              s = offBit( s, BC_N_S );
            }
          }
          
          // Y_PLUS
          if ( !((pn.nID[Y_PLUS] < 0) && (j == jx)) ) {
            if ( !IS_FLUID(bx[m_n]) ) {
              s = offBit( s, BC_N_N );
            }
          }
          
          // Z_MINUS
          if ( !((pn.nID[Z_MINUS] < 0) && (k == 1)) ) {
            if ( !IS_FLUID(bx[m_b]) ) {
              s = offBit( s, BC_N_B );
            }
          }
          
          // Z_PLUS
          if ( !((pn.nID[Z_PLUS] < 0) && (k == kx)) ) {
            if ( !IS_FLUID(bx[m_t]) ) {
              s = offBit( s, BC_N_T );
            }
          }
          
          bx[m_p] = s;
        }
      }
    }
  }
  
  // wall locationフラグ
  unsigned long c = 0;
  
  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        m_p = FBUtility::getFindexS3D(size, guide, i, j, k);
        s = bx[m_p];
        
        if ( IS_FLUID( s ) ) {
          m_e = FBUtility::getFindexS3D(size, guide, i+1, j  , k  );
          m_w = FBUtility::getFindexS3D(size, guide, i-1, j  , k  );
          m_n = FBUtility::getFindexS3D(size, guide, i  , j+1, k  );
          m_s = FBUtility::getFindexS3D(size, guide, i  , j-1, k  );
          m_t = FBUtility::getFindexS3D(size, guide, i  , j  , k+1);
          m_b = FBUtility::getFindexS3D(size, guide, i  , j  , k-1);
          
          // X_MINUS
          if ( !IS_FLUID(bx[m_w]) ) { s = onBit( s, FACING_W ); c++; }
          
          // X_PLUS
          if ( !IS_FLUID(bx[m_e]) ) { s = onBit( s, FACING_E ); c++; }
          
          // Y_MINUS
          if ( !IS_FLUID(bx[m_s]) ) { s = onBit( s, FACING_S ); c++; }
          
          // Y_PLUS
          if ( !IS_FLUID(bx[m_n]) ) { s = onBit( s, FACING_N ); c++; }
          
          // Z_MINUS
          if ( !IS_FLUID(bx[m_b]) ) { s = onBit( s, FACING_B ); c++; }
          
          // Z_PLUS
          if ( !IS_FLUID(bx[m_t]) ) { s = onBit( s, FACING_T ); c++; }
          
          // 収束判定の有効フラグ，計算内部領域の流体セルのみ
          s = onBit(s, VLD_CNVG);
          
          bx[m_p] = s;
        }
      }
    }
  }
  
  if ( pn.numProc > 1 ) {
    unsigned long tmp = c;
    MPI_Allreduce(&tmp, &c, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  }
  
  return c;
}

/**
 @fn unsigned VoxInfo::encPbit_N_Cut(unsigned* bx, float* cut, const bool convergence)
 @brief bcp[]に壁面境界の圧力ノイマン条件のビットフラグと固体に隣接するFセルに方向フラグ，収束判定の有効フラグをカット情報からエンコードする
 @param[in/out] bx BCindex P
 @param[in] cut 距離情報
 @param convergence カットのあるセルは収束判定をしないオプション（trueの時）
 @retval 固体表面セル数
 @note 
 - 流体セルのうち，固体セルに隣接する面のノイマンフラグをゼロにする．ただし，内部領域のみ．
 - 固体セルに隣接する流体セルに方向フラグを付与する．全内部領域．
 */
unsigned long VoxInfo::encPbit_N_Cut(unsigned* bx, float* cut, const bool convergence)
{
  unsigned m_p, m;
  unsigned register s;
  float cp_e, cp_w, cp_n, cp_s, cp_t, cp_b;
  float* ct;
  
  int ix = (int)size[0];
  int jx = (int)size[1];
  int kx = (int)size[2];
  
  // ノイマンフラグ
  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        m_p = FBUtility::getFindexS3D(size, guide, i, j, k);
        s = bx[m_p];
        
        if ( IS_FLUID( s ) ) {
          
          m = FBUtility::getFindexS3Dcut(size, guide, 0, i, j, k);
          ct = &cut[m];

          // X_MINUS
          //if ( !((pn.nID[X_MINUS] < 0) && (i == 1)) ) {
            if (ct[0] < 1.0) {          // 交点があるなら壁面なのでノイマン条件をセット
              s = offBit( s, BC_N_W );
            }
          //}
          
          // X_PLUS
          //if ( !((pn.nID[X_PLUS] < 0) && (i == ix)) ) {
            if (ct[1] < 1.0) {
              s = offBit( s, BC_N_E );
            }
          //}
          
          // Y_MINUS
          //if ( !((pn.nID[Y_MINUS] < 0) && (j == 1)) ) {
            if (ct[2] < 1.0) {
              s = offBit( s, BC_N_S );
            }
          //}
          
          // Y_PLUS
          //if ( !((pn.nID[Y_PLUS] < 0) && (j == jx)) ) {
            if (ct[3] < 1.0) {
              s = offBit( s, BC_N_N );
            }
          //}
          
          // Z_MINUS
          //if ( !((pn.nID[Z_MINUS] < 0) && (k == 1)) ) {
            if (ct[4] < 1.0) {
              s = offBit( s, BC_N_B );
            }
          //}
          
          // Z_PLUS
          //if ( !((pn.nID[Z_PLUS] < 0) && (k == kx)) ) {
            if (ct[5] < 1.0) {
              s = offBit( s, BC_N_T );
            }
          //}
          
          bx[m_p] = s;
        }
      }
    }
  }
  
  // wall locationフラグ
  unsigned long c = 0;
  float* pos; 
  float q;
  
  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        m_p = FBUtility::getFindexS3D(size, guide, i, j, k);
        s = bx[m_p];
        
        if ( IS_FLUID( s ) ) {
          m = FBUtility::getFindexS3Dcut(size, guide, 0, i, j, k);
          pos = &cut[m];
          
          if ( pos[X_MINUS] != 1.0 ) { s = onBit( s, FACING_W ); c++; }
          if ( pos[X_PLUS]  != 1.0 ) { s = onBit( s, FACING_E ); c++; }
          if ( pos[Y_MINUS] != 1.0 ) { s = onBit( s, FACING_S ); c++; }
          if ( pos[Y_PLUS]  != 1.0 ) { s = onBit( s, FACING_N ); c++; }
          if ( pos[Z_MINUS] != 1.0 ) { s = onBit( s, FACING_B ); c++; }
          if ( pos[Z_PLUS]  != 1.0 ) { s = onBit( s, FACING_T ); c++; }
          
          bx[m_p] = s;
        }
      }
    }
  }
  
  if ( pn.numProc > 1 ) {
    unsigned long tmp = c;
    MPI_Allreduce(&tmp, &c, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  }
  
  
  // 収束判定の有効フラグ
  float q0, q1, q2, q3, q4, q5;
  unsigned long g=0;
  
  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        m_p = FBUtility::getFindexS3D(size, guide, i, j, k);
        s = bx[m_p];
        
        if ( IS_FLUID( s ) ) {
          m = FBUtility::getFindexS3Dcut(size, guide, 0, i, j, k);
          pos = &cut[m];
          
          q0 = floor(pos[0]);
          q1 = floor(pos[1]);
          q2 = floor(pos[2]);
          q3 = floor(pos[3]);
          q4 = floor(pos[4]);
          q5 = floor(pos[5]);
          
          // 収束判定の有効フラグ 全周カットがあるセルは孤立セルとして固体セルへ変更
          if ( (q0+q1+q2+q3+q4+q5) == 0.0 ) {
            s = offBit(s, VLD_CNVG);    // Out of scope
            s = offBit(s, STATE_BIT );  // Solid
            s = offBit(s, ACTIVE_BIT ); // Inactive
            g++;
          }
          else {
            s = onBit(s, VLD_CNVG);
          }
          
          bx[m_p] = s;
        }
      }
    }
  }
  
  
  if ( pn.numProc > 1 ) {
    unsigned long tmp = g;
    MPI_Allreduce(&tmp, &g, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  }
  
  Hostonly_ printf("\tThe number of cells which are changed to INACTIVE and SOLID because of all faces are cut = %ld\n\n", g);
  
  
  // カットのあるセルの収束判定をしないオプション
  if ( convergence ) {
    g = 0;
    
    for (int k=1; k<=kx; k++) {
      for (int j=1; j<=jx; j++) {
        for (int i=1; i<=ix; i++) {
          m_p = FBUtility::getFindexS3D(size, guide, i, j, k);
          s = bx[m_p];
          
          if ( IS_FLUID( s ) ) {
            m = FBUtility::getFindexS3Dcut(size, guide, 0, i, j, k);
            pos = &cut[m];
            
            q0 = floor(pos[0]);
            q1 = floor(pos[1]);
            q2 = floor(pos[2]);
            q3 = floor(pos[3]);
            q4 = floor(pos[4]);
            q5 = floor(pos[5]);
            
            // 収束判定の有効フラグ 
            if ( (q0+q1+q2+q3+q4+q5) != 6.0 ) {
              s = offBit(s, VLD_CNVG);    // Out of scope  @todo check
              g++;
            }
            
            bx[m_p] = s;
          }
        }
      }
    }
    
    if ( pn.numProc > 1 ) {
      unsigned long tmp = g;
      MPI_Allreduce(&tmp, &g, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
    }
    
    Hostonly_ printf("\tThe number of cells which are excluded to convergence judgement by cut = %ld\n\n", g);
    
  }
  
  
  return c;
}

/**
 @fn unsigned VoxInfo::encPbit_N_IBC(unsigned order, unsigned id, int* mid, unsigned* bcd, unsigned* bcp, int deface)
 @brief 計算領域内部のコンポーネントのNeumannフラグをbcp[]にエンコードする
 @retval エンコードしたセル数
 @param order cmp[]のエントリ
 @param id  CellID
 @param mid ボクセル配列
 @param bcd BCindex ID
 @param bcp BCindex P
 @param deface 面を指定するid
 @note
 - 対象セルが流体セル，かつターゲットの隣接セルが固体の場合，隣接する面にNeumannフラグをエンコードする
 - 同種のBCは1セルに一つだけ
 */
unsigned long VoxInfo::encPbit_N_IBC(const unsigned order, 
                                     const unsigned id, 
                                     const int* mid, 
                                     unsigned* bcd, 
                                     unsigned* bcp, 
                                     const int deface)
{
  unsigned long g=0; 
  unsigned m_e, m_w, m_n, m_s, m_t, m_b;
  unsigned s_e, s_w, s_n, s_s, s_t, s_b;
  int      c_p, c_e, c_w, c_n, c_s, c_t, c_b;
  unsigned register s, d, m;
  
  int ix = (int)size[0];
  int jx = (int)size[1];
  int kx = (int)size[2];
  int gd = (int)guide;
  
  int idd = (int)id;
  
  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        m = FBUtility::getFindexS3D(size, guide, i, j, k);
        d = bcd[m];
        s = bcp[m];
        c_p = mid[m];
        
        if ( IS_FLUID( s ) && (c_p == deface) ) { // 対象セルが流体セル，かつdefaceである
          m_e = FBUtility::getFindexS3D(size, guide, i+1, j  , k  );
          m_w = FBUtility::getFindexS3D(size, guide, i-1, j  , k  );
          m_n = FBUtility::getFindexS3D(size, guide, i  , j+1, k  );
          m_s = FBUtility::getFindexS3D(size, guide, i  , j-1, k  );
          m_t = FBUtility::getFindexS3D(size, guide, i  , j  , k+1);
          m_b = FBUtility::getFindexS3D(size, guide, i  , j  , k-1);
          
          c_e = mid[m_e];
          c_w = mid[m_w];
          c_n = mid[m_n];
          c_s = mid[m_s];
          c_t = mid[m_t];
          c_b = mid[m_b];
          
          s_e = bcp[m_e];
          s_w = bcp[m_w];
          s_n = bcp[m_n];
          s_s = bcp[m_s];
          s_t = bcp[m_t];
          s_b = bcp[m_b];
          
          // 外部境界面は除外
          // X_MINUS
          //if ( !((pn.nID[X_MINUS] < 0) && (i == 1)) ) {
            if ( !IS_FLUID( s_w ) && (c_w == idd) ) { // Wセルが指定ID，かつ固体である
              d |= order; // dにエントリをエンコード
              s = offBit( s, BC_N_W );
              g++;
            }
          //}
          
          // X_PLUS
          //if ( !((pn.nID[X_PLUS] < 0) && (i == ix)) ) {
            if ( !IS_FLUID( s_e ) && (c_e == idd) ) {
              d |= order;
              s = offBit( s, BC_N_E );
              g++;
            }
          //}
          
          // Y_MINUS
          //if ( !((pn.nID[Y_MINUS] < 0) && (j == 1)) ) {
            if ( !IS_FLUID( s_s ) && (c_s == idd) ) {
              d |= order;
              s = offBit( s, BC_N_S );
              g++;
            }
          //}
          
          // Y_PLUS
          //if ( !((pn.nID[Y_PLUS] < 0) && (j == jx)) ) {
            if ( !IS_FLUID( s_n ) && (c_n == idd) ) {
              d |= order;
              s = offBit( s, BC_N_N );
              g++;
            }
          //}
          
          // Z_MINUS
          //if ( !((pn.nID[Z_MINUS] < 0) && (k == 1)) ) {
            if ( !IS_FLUID( s_b ) && (c_b == idd) ) {
              d |= order;
              s = offBit( s, BC_N_B );
              g++;
            }
          //}
          
          // Z_PLUS
          //if ( !((pn.nID[Z_PLUS] < 0) && (k == kx)) ) {
            if ( !IS_FLUID( s_t ) && (c_t == idd) ) {
              d |= order;
              s = offBit( s, BC_N_T );
              g++;
            }
          //}
          
          bcd[m] = d;
          bcp[m] = s;
        }
      }
    }
  }
  
  if ( pn.numProc > 1 ) {
    unsigned long tmp = g;
    MPI_Allreduce(&tmp, &g, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  }
  
  return g;
}

/**
 @fn void VoxInfo::encPbit_OBC(int face, unsigned* bx, string key, bool dir)
 @brief 外部境界に接するセルにおいて，bx[]に圧力境界条件keyに対応するビットフラグを設定する
 @param face 外部境界面番号
 @param bx BCindex P
 @param key Dirichlet or Neumann
 @param dir 壁面の場合(true)，方向フラグをON
 @retval 固体表面セル数
 @note 
 - 流体セルに対してのみ，1-Normal, 0-BC
 - 固体セルに隣接する面のノイマンフラグをゼロにし，方向フラグを立てる
 */
void VoxInfo::encPbit_OBC(int face, unsigned* bx, string key, bool dir)
{
  int i,j,k;
  unsigned m;
  unsigned register s;
  
  int ix = (int)size[0];
  int jx = (int)size[1];
  int kx = (int)size[2];
  int gd = (int)guide;
  switch (face) {
    case X_MINUS:
      if( pn.nID[X_MINUS] < 0 ){
        i = 1;
        if ("Neumann"==key) {
          for (k=1; k<=kx; k++) {
            for (j=1; j<=jx; j++) {
              m = FBUtility::getFindexS3D(size, guide, i, j, k);
              s = bx[m];
              if ( IS_FLUID(s) ) {
                if (dir) s = onBit( s, FACING_W );
                bx[m] = offBit( s, BC_N_W );
              }
            }
          }
        }
        else {
          for (k=1; k<=kx; k++) {
            for (j=1; j<=jx; j++) {
              m = FBUtility::getFindexS3D(size, guide, i, j, k);
              s = bx[m];
              if ( IS_FLUID(s) ) {
                if (dir) s = onBit( s, FACING_W );
                bx[m] = offBit( s, BC_D_W );
              }
            }
          }
        }
      }
      break;
      
    case X_PLUS:
      if( pn.nID[X_PLUS] < 0 ){
        i = ix;
        if ("Neumann"==key) {
          for (k=1; k<=kx; k++) {
            for (j=1; j<=jx; j++) {
              m = FBUtility::getFindexS3D(size, guide, i, j, k);
              s = bx[m];
              if ( IS_FLUID(s) ) {
                if (dir) s = onBit( s, FACING_E );
                bx[m] = offBit( s, BC_N_E );
              }
            }
          }
        }
        else {
          for (k=1; k<=kx; k++) {
            for (j=1; j<=jx; j++) {
              m = FBUtility::getFindexS3D(size, guide, i, j, k);
              s = bx[m];
              if ( IS_FLUID(s) ) {
                if (dir) s = onBit( s, FACING_E );
                bx[m] = offBit( s, BC_D_E );
              }
            }
          }
        }
      }
      break;
      
    case Y_MINUS:
      if( pn.nID[Y_MINUS] < 0 ){
        j = 1;
        if ("Neumann"==key) {
          for (k=1; k<=kx; k++) {
            for (i=1; i<=ix; i++) {
              m = FBUtility::getFindexS3D(size, guide, i, j, k);
              s = bx[m];
              if ( IS_FLUID(s) ) {
                if (dir) s = onBit( s, FACING_S );
                bx[m] = offBit( s, BC_N_S );
              }
            }
          }
        }
        else {
          for (k=1; k<=kx; k++) {
            for (i=1; i<=ix; i++) {
              m = FBUtility::getFindexS3D(size, guide, i, j, k);
              s = bx[m];
              if ( IS_FLUID(s) ) {
                if (dir) s = onBit( s, FACING_S );
                bx[m] = offBit( s, BC_D_S );
              }
            }
          }
        }
      }
      break;
      
    case Y_PLUS:
      if( pn.nID[Y_PLUS] < 0 ){
        j = jx;
        if ("Neumann"==key) {
          for (k=1; k<=kx; k++) {
            for (i=1; i<=ix; i++) {
              m = FBUtility::getFindexS3D(size, guide, i, j, k);
              s = bx[m];
              if ( IS_FLUID(s) ) {
                if (dir) s = onBit( s, FACING_N );
                bx[m] = offBit( s, BC_N_N );
              }
            }
          }
        }
        else {
          for (k=1; k<=kx; k++) {
            for (i=1; i<=ix; i++) {
              m = FBUtility::getFindexS3D(size, guide, i, j, k);
              s = bx[m];
              if ( IS_FLUID(s) ) {
                if (dir) s = onBit( s, FACING_N );
                bx[m] = offBit( s, BC_D_N );
              }
            }
          }
        }
      }
      break;
      
    case Z_MINUS:
      if( pn.nID[Z_MINUS] < 0 ){
        k = 1;
        if ("Neumann"==key) {
          for (j=1; j<=jx; j++) {
            for (i=1; i<=ix; i++) {
              m = FBUtility::getFindexS3D(size, guide, i, j, k);
              s = bx[m];
              if ( IS_FLUID(s) ) {
                if (dir) s = onBit( s, FACING_B );
                bx[m] = offBit( s, BC_N_B );
              }
            }
          }
        }
        else {
          for (j=1; j<=jx; j++) {
            for (i=1; i<=ix; i++) {
              m = FBUtility::getFindexS3D(size, guide, i, j, k);
              s = bx[m];
              if ( IS_FLUID(s) ) {
                if (dir) s = onBit( s, FACING_B );
                bx[m] = offBit( s, BC_D_B );
              }
            }
          }
        }
      }
      break;
      
    case Z_PLUS:
      if( pn.nID[Z_PLUS] < 0 ){
        k = kx;
        if ("Neumann"==key) {
          for (j=1; j<=jx; j++) {
            for (i=1; i<=ix; i++) {
              m = FBUtility::getFindexS3D(size, guide, i, j, k);
              s = bx[m];
              if ( IS_FLUID(s) ) {
                if (dir) s = onBit( s, FACING_T );
                bx[m] = offBit( s, BC_N_T );
              }
            }
          }
        }
        else {
          for (j=1; j<=jx; j++) {
            for (i=1; i<=ix; i++) {
              m = FBUtility::getFindexS3D(size, guide, i, j, k);
              s = bx[m];
              if ( IS_FLUID(s) ) {
                if (dir) s = onBit( s, FACING_T );
                bx[m] = offBit( s, BC_D_T );
              }
            }
          }
        }
      }
      break;
  } // end of switch
}

/**
 @fn unsigned VoxInfo::encVbit_IBC(unsigned order, unsigned id, int* mid, unsigned* bv, int deface, unsigned* bp)
 @brief bv[]にVBCの境界条件に必要な情報をエンコードし，流入流出の場合にbp[]の隣接壁の方向フラグを除く．境界条件指定キーセルのSTATEを流体に変更する
 @retval エンコードしたセル数
 @param order cmp[]のエントリ番号
 @param id  CellID
 @param mid ボクセル配列
 @param bv BCindex V
 @param deface 面を指定するid
 @param bp BCindex P
 @note 対象セルは流体セルであること
 */
unsigned long VoxInfo::encVbit_IBC(const unsigned order, 
                                   const unsigned id, 
                                   const int* mid, 
                                   unsigned* bv, 
                                   const int deface, 
                                   unsigned* bp)
{
  unsigned long g=0;
  unsigned m_p, m_e, m_w, m_n, m_s, m_t, m_b;
  int      c_p, c_e, c_w, c_n, c_s, c_t, c_b;
  unsigned register s, q;
  
  int ix = (int)size[0];
  int jx = (int)size[1];
  int kx = (int)size[2];
  int gd = (int)guide;
  
  int idd = (int)id;
  
  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        m_p = FBUtility::getFindexS3D(size, guide, i  , j  , k  );
        c_p = mid[m_p];
        
        s = bv[m_p];
        q = bp[m_p];
        
        if ( (c_p == deface) && IS_FLUID(s) ) { // テストセルが流体，かつdeface
          m_e = FBUtility::getFindexS3D(size, guide, i+1, j  , k  );
          m_w = FBUtility::getFindexS3D(size, guide, i-1, j  , k  );
          m_n = FBUtility::getFindexS3D(size, guide, i  , j+1, k  );
          m_s = FBUtility::getFindexS3D(size, guide, i  , j-1, k  );
          m_t = FBUtility::getFindexS3D(size, guide, i  , j  , k+1);
          m_b = FBUtility::getFindexS3D(size, guide, i  , j  , k-1);
          
          c_e = mid[m_e];
          c_w = mid[m_w];
          c_n = mid[m_n];
          c_s = mid[m_s];
          c_t = mid[m_t];
          c_b = mid[m_b];
          
          // X-
          if ( c_w == idd ) {
            //if ( !((pn.nID[X_MINUS] < 0) && (i == 1)) ) {
              if ( !IS_FLUID( bv[m_w] ) ) { // 流入出セルが固体であるかをチェックする
                s |= (order << BC_FACE_W);
                q = offBit(q, FACING_W);
                g++;
              }
              else {
                Hostonly_ printf("Error : SpecVel/Outflow face should be solid at [%d,%d,%d]\n",i,j,k);
                Exit(0);
              }
            //}
          }
          
          // X+
          if ( c_e == idd ) {
            //if ( !((pn.nID[X_PLUS] < 0) && (i == ix)) ) {
              if ( !IS_FLUID( bv[m_e] ) ) { 
                s |= (order << BC_FACE_E);
                q = offBit(q, FACING_E);
                g++;
              }
              else {
                Hostonly_ printf("Error : SpecVel/Outflow face should be solid at [%d,%d,%d]\n",i,j,k);
                Exit(0);
              }
            //}
          }
          
          // Y-
          if ( c_s == idd ) {
            //if ( !((pn.nID[Y_MINUS] < 0) && (j == 1)) ) {
              if ( !IS_FLUID( bv[m_s] ) ) { 
                s |= (order << BC_FACE_S);
                q = offBit(q, FACING_S);
                g++;
              }
              else {
                Hostonly_ printf("Error : SpecVel/Outflow face should be solid at [%d,%d,%d]\n",i,j,k);
                Exit(0);
              }
            //}
          }
          
          // Y+
          if ( c_n == idd ) {
            //if ( !((pn.nID[Y_PLUS] < 0) && (j == jx)) ) {
              if ( !IS_FLUID( bv[m_n] ) ) { 
                s |= (order << BC_FACE_N);
                q = offBit(q, FACING_N);
                g++;
              }
              else {
                Hostonly_ printf("Error : SpecVel/Outflow face should be solid at [%d,%d,%d]\n",i,j,k);
                Exit(0);
              }
            //}
          }
          
          // Z-
          if ( c_b == idd ) {
            //if ( !((pn.nID[Z_MINUS] < 0) && (k == 1)) ) {
              if ( !IS_FLUID( bv[m_b] ) ) { 
                s |= (order << BC_FACE_B);
                q = offBit(q, FACING_B);
                g++;
              }
              else {
                Hostonly_ printf("Error : SpecVel/Outflow face should be solid at [%d,%d,%d]\n",i,j,k);
                Exit(0);
              }
            //}
          }
          
          // Z+
          if ( c_t == idd ) {
            //if ( !((pn.nID[Z_PLUS] < 0) && (k == kx)) ) {
              if ( !IS_FLUID( bv[m_t] ) ) { 
                s |= (order << BC_FACE_T);
                q = offBit(q, FACING_T);
                g++;
              }
              else {
                Hostonly_ printf("Error : SpecVel/Outflow face should be solid at [%d,%d,%d]\n",i,j,k);
                Exit(0);
              }
            //}
          }
          
          bv[m_p] = s;
          bp[m_p] = q;
        }
      } // i-loop
    }
  }
  
  // iddセルのSTATE_BITをFLUIDに変更 >> bvのみ違う参照をさせるため
  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        m_p = FBUtility::getFindexS3D(size, guide, i  , j  , k  );
        c_p = mid[m_p];
        
        s = bv[m_p];
        
        if ( (c_p == deface) && IS_FLUID(s) ) {
          m_e = FBUtility::getFindexS3D(size, guide, i+1, j  , k  );
          m_w = FBUtility::getFindexS3D(size, guide, i-1, j  , k  );
          m_n = FBUtility::getFindexS3D(size, guide, i  , j+1, k  );
          m_s = FBUtility::getFindexS3D(size, guide, i  , j-1, k  );
          m_t = FBUtility::getFindexS3D(size, guide, i  , j  , k+1);
          m_b = FBUtility::getFindexS3D(size, guide, i  , j  , k-1);
          
          c_e = mid[m_e];
          c_w = mid[m_w];
          c_n = mid[m_n];
          c_s = mid[m_s];
          c_t = mid[m_t];
          c_b = mid[m_b];
          
          // X-
          if ( c_w == idd ) {
            //if ( !((pn.nID[X_MINUS] < 0) && (i == 1)) ) {
              bv[m_w] = onBit( bv[m_w], STATE_BIT );
            //}
          }
          
          // X+
          if ( c_e == idd ) {
            //if ( !((pn.nID[X_PLUS] < 0) && (i == ix)) ) {
              bv[m_e] = onBit( bv[m_e], STATE_BIT );
            //}
          }
          
          // Y-
          if ( c_s == idd ) {
            //if ( !((pn.nID[Y_MINUS] < 0) && (j == 1)) ) {
              bv[m_s] = onBit( bv[m_s], STATE_BIT );
            //}
          }
          
          // Y+
          if ( c_n == idd ) {
            //if ( !((pn.nID[Y_PLUS] < 0) && (j == jx)) ) {
              bv[m_n] = onBit( bv[m_n], STATE_BIT );
            //}
          }
          
          // Z-
          if ( c_b == idd ) {
            //if ( !((pn.nID[Z_MINUS] < 0) && (k == 1)) ) {
              bv[m_b] = onBit( bv[m_b], STATE_BIT );
            //}
          }
          
          // Z+
          if ( c_t == idd ) {
            //if ( !((pn.nID[Z_PLUS] < 0) && (k == kx)) ) {
              bv[m_t] = onBit( bv[m_t], STATE_BIT );
            //}
          }
        }
      }
    }
  }
  
  if ( pn.numProc > 1 ) {
    unsigned long tmp = g;
    MPI_Allreduce(&tmp, &g, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  }
  
  return g;
}

/**
 @fn unsigned VoxInfo::encVbit_IBC_Cut(const unsigned order, const unsigned id, unsigned* bv, unsigned* bp, const float* cut, 
                                        const int* cut_id, const float* vec, const unsigned bc_dir)
 @brief bv[]にVBCの境界条件に必要な情報をエンコードし，流入流出の場合にbp[]の隣接壁の方向フラグを除く．境界条件指定キーセルのSTATEを流体に変更する
 @retval エンコードしたセル数
 @param order cmp[]のインデクス
 @param id 境界条件ID
 @param bv BCindex V
 @param bp BCindex P
 @param cut 距離情報
 @param cut_id カット点ID
 @param vec 境界面の法線
 @param bc_dir BCの位置，指定法線と同じ側(1)か反対側(2)
 @note 指定法線とセルのカット方向ベクトルの内積で判断，vspecとoutflowなのでbp[]のVBC_UWDにマスクビットを立てる
 */
unsigned long VoxInfo::encVbit_IBC_Cut(const unsigned order, 
                                       const unsigned id, 
                                       unsigned* bv, 
                                       unsigned* bp, 
                                       const float* cut, 
                                       const int* cut_id, 
                                       const float* vec, 
                                       const unsigned bc_dir)
{
  unsigned long g=0;
  
  unsigned register s, q;
  const float *pos;
  int    bid;
  unsigned m_p, m_c;
  int    id_e, id_w, id_n, id_s, id_t, id_b;
  float pp;
  
  FB::Vec3f nv(vec);
  
  // 反対方向のとき符号反転 > 内積が負のときに対象位置となる
  if ( bc_dir == CompoList::opposite_direction ) {
    nv *= -1.0;
  }
  
  FB::Vec3f e_w(-1.0,  0.0,  0.0);
  FB::Vec3f e_e(+1.0,  0.0,  0.0);
  FB::Vec3f e_s( 0.0, -1.0,  0.0);
  FB::Vec3f e_n( 0.0, +1.0,  0.0);
  FB::Vec3f e_b( 0.0,  0.0, -1.0);
  FB::Vec3f e_t( 0.0,  0.0, +1.0);
  
  int ix = (int)size[0];
  int jx = (int)size[1];
  int kx = (int)size[2];
  
  int idd = (int)id;
  
  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        
        m_c = FBUtility::getFindexS3Dcut(size, guide, 0, i, j, k);
        pos = &cut[m_c];
        pp = pos[0]+pos[1]+pos[2]+pos[3]+pos[4]+pos[5];
        
        if ( pp < 6.0 ) { // 6方向のうちいずれかにカットがある

          m_p = FBUtility::getFindexS3D(size, guide, i, j, k);
          bid = cut_id[m_p];
          s   = bv[m_p];
          q   = bp[m_p];
          
          if ( IS_FLUID(s) ) { // 流体かつdefaceであるセルがテスト対象
            // 各方向のID
            id_w = get_BID5(X_MINUS, bid);
            id_e = get_BID5(X_PLUS,  bid);
            id_s = get_BID5(Y_MINUS, bid);
            id_n = get_BID5(Y_PLUS,  bid);
            id_b = get_BID5(Z_MINUS, bid);
            id_t = get_BID5(Z_PLUS,  bid);
            
            // X-
            if ( (id_w == idd) && (dot(e_w, nv) < 0.0) ) {
              s |= (order << BC_FACE_W);
              q = offBit(q, FACING_W);
              q = onBit(q, VBC_UWD);
              g++;
            }
            
            // X+
            if ( (id_e == idd) && (dot(e_e, nv) < 0.0) ) {
              s |= (order << BC_FACE_E);
              q = offBit(q, FACING_E);
              q = onBit(q, VBC_UWD);
              g++;
            }
            
            // Y-
            if ( (id_s == idd) && (dot(e_s, nv) < 0.0) ) {
              s |= (order << BC_FACE_S);
              q = offBit(q, FACING_S);
              q = onBit(q, VBC_UWD);
              g++;
            }
            
            // Y+
            if ( (id_n == idd) && (dot(e_n, nv) < 0.0) ) {
              s |= (order << BC_FACE_N);
              q = offBit(q, FACING_N);
              q = onBit(q, VBC_UWD);
              g++;
            }
            
            // Z-
            if ( (id_b == idd) && (dot(e_b, nv) < 0.0) ) {
              s |= (order << BC_FACE_B);
              q = offBit(q, FACING_B);
              q = onBit(q, VBC_UWD);
              g++;
            }
            
            // Z+
            if ( (id_t == idd) && (dot(e_t, nv) < 0.0) ) {
              s |= (order << BC_FACE_T);
              q = offBit(q, FACING_T);
              q = onBit(q, VBC_UWD);
              g++;
            }
          } // if fluid
          
          bv[m_p] = s;
          bp[m_p] = q;
          
        } // if cut
        
      } // i-loop
    }
  }
  
  // check for debug
# if 0
  int m_flag;
  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        
        m_p = FBUtility::getFindexS3D(size, guide, i, j, k);
        s = bv[m_p];
        
        m_flag = 0;
        
        if ( GET_FACE_BC(s, BC_FACE_W) == order ) m_flag = 1;
        if ( GET_FACE_BC(s, BC_FACE_E) == order ) m_flag = 2;
        if ( GET_FACE_BC(s, BC_FACE_S) == order ) m_flag = 3;
        if ( GET_FACE_BC(s, BC_FACE_N) == order ) m_flag = 4;
        if ( GET_FACE_BC(s, BC_FACE_B) == order ) m_flag = 5;
        if ( GET_FACE_BC(s, BC_FACE_T) == order ) m_flag = 6;
        
        if ( m_flag != 0 ) printf("%3d %3d %3d >> %d\n", i, j, k, m_flag);
      }
    }
  }
#endif
  
  
  if ( pn.numProc > 1 ) {
    unsigned long tmp = g;
    MPI_Allreduce(&tmp, &g, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  }
  
  return g;
}


/**
 @fn void VoxInfo::encVbit_OBC(int face, unsigned* bv, string key, const bool enc_sw, string chk, unsigned* bp, const bool enc_uwd)
 @brief 外部境界に接するセルにおいて，各種速度境界条件に対応する媒質をチェックし，bv[]にビットフラグをセットする
 @param face 外部境界面番号
 @param bv BCindex V
 @param key fluid or solid　指定するBCが要求するガイドセルの状態 >> エラーチェックに使う
 @param enc_sw trueのとき，エンコードする．falseの場合にはガイドセルの状態チェックのみ
 @param chk ガイドセルの状態をチェックするかどうかを指定
 @param bp BCindex P
 @param enc_uwd trueのとき，1次精度のスイッチオン
 @note 
 - 外部境界条件の実装には，流束型とディリクレ型の2種類がある．
 - adjMedium_on_GC()でガイドセル上のIDを指定済み．指定BCとの適合性をチェックする
 */
void VoxInfo::encVbit_OBC(int face, unsigned* bv, string key, const bool enc_sw, string chk, unsigned* bp, const bool enc_uwd)
{
  int i, j, k;
  unsigned m, mt, sw, cw;
  unsigned register s, q;
  
  int ix = (int)size[0];
  int jx = (int)size[1];
  int kx = (int)size[2];
  
  ( "fluid" == key ) ? sw=1 : sw=0;
  ( "check" == chk ) ? cw=1 : cw=0;
  
  switch (face) {
    case X_MINUS:
      if( pn.nID[X_MINUS] < 0 ){ // 外部境界をもつノードのみ
        for (k=1; k<=kx; k++) {
          for (j=1; j<=jx; j++) {
            m = FBUtility::getFindexS3D(size, guide, 1  , j, k);
            mt= FBUtility::getFindexS3D(size, guide, 0  , j, k);
            
            s = bv[m];
            q = bp[mt];
            
            if ( IS_FLUID(s) ) {
              // エンコード処理
              if ( enc_sw ) {
                bv[m]  = s | (OBC_MASK << BC_FACE_W); // OBC_MASK==31 外部境界条件のフラグ
              }
              
              // 外部境界で安定化のため，スキームを1次精度にする
              if ( enc_uwd ) bp[mt] = onBit(q, VBC_UWD);
              
              // チェック
              if ( cw == 1 ) {
                
                if ( sw==1 ) { // ガイドセルが流体であることを要求
                  if ( !IS_FLUID(bv[mt]) ) {
                    Hostonly_ printf("Error : Fluid cell at Outer boundary X- is required\n");
                    Exit(0);
                  }
                }
                else { // ガイドセルが固体であることを要求
                  if ( IS_FLUID(bv[mt]) ) {
                    Hostonly_ printf("Error : Solid cell at Outer boundary X- is required\n");
                    Exit(0);
                  }
                }             
              }
              
            }// IS_FLUID
          }
        }        
      }
      break;
      
    case X_PLUS:
      if( pn.nID[X_PLUS] < 0 ){
        for (k=1; k<=kx; k++) {
          for (j=1; j<=jx; j++) {
            m = FBUtility::getFindexS3D(size, guide, ix  , j, k);
            mt= FBUtility::getFindexS3D(size, guide, ix+1, j, k);
            
            s = bv[m];
            q = bp[mt];
            
            if ( IS_FLUID(s) ) {
              
              // エンコード処理
              if ( enc_sw ) {
                bv[m]  = s | (OBC_MASK << BC_FACE_E);
              }
              
              // 外部境界で安定化のため，スキームを1次精度にする
              if ( enc_uwd ) bp[mt] = onBit(q, VBC_UWD);
              
              // チェック
              if ( cw == 1 ) {
                
                if ( sw==1 ) {
                  if ( !IS_FLUID(bv[mt]) ) {
                    Hostonly_ printf("Error : Fluid cell at Outer boundary X+ is required\n");
                    Exit(0);
                  }
                }
                else {
                  if ( IS_FLUID(bv[mt]) ) {
                    Hostonly_ printf("Error : Solid cell at Outer boundary X+ is required\n");
                    Exit(0);
                  }
                }
              }
              
            }// IS_FLUID
          }
        }
      }
      break;
      
    case Y_MINUS:
      if( pn.nID[Y_MINUS] < 0 ){
        for (k=1; k<=kx; k++) {
          for (i=1; i<=ix; i++) {
            m = FBUtility::getFindexS3D(size, guide, i, 1  , k);
            mt= FBUtility::getFindexS3D(size, guide, i, 0  , k);
            
            s = bv[m];
            q = bp[mt];
            
            if ( IS_FLUID(s) ) {
              
              // エンコード処理
              if ( enc_sw ) {
                bv[m]  = s | (OBC_MASK << BC_FACE_S);
              }
              
              // 外部境界で安定化のため，スキームを1次精度にする
              if ( enc_uwd ) bp[mt] = onBit(q, VBC_UWD);
              
              // チェック
              if ( cw == 1 ) {
                
                if ( sw==1 ) {
                  if ( !IS_FLUID(bv[mt]) ) {
                    Hostonly_ printf("Error : Fluid cell at Outer boundary Y- is required\n");
                    Exit(0);
                  }
                }
                else {
                  if ( IS_FLUID(bv[mt]) ) {
                    Hostonly_ printf("Error : Solid cell at Outer boundary Y- is required\n");
                    Exit(0);
                  }
                }
              }
              
            }// IS_FLUID
          }
        }
      }
      break;
      
    case Y_PLUS:
      if( pn.nID[Y_PLUS] < 0 ){
        for (k=1; k<=kx; k++) {
          for (i=1; i<=ix; i++) {
            m = FBUtility::getFindexS3D(size, guide, i, jx  , k);
            mt= FBUtility::getFindexS3D(size, guide, i, jx+1, k);
            
            s = bv[m];
            q = bp[mt];
            
            if ( IS_FLUID(s) ) {
              
              // エンコード処理
              if ( enc_sw ) {
                bv[m]  = s | (OBC_MASK << BC_FACE_N);
              }
              
              // 外部境界で安定化のため，スキームを1次精度にする
              if ( enc_uwd ) bp[mt] = onBit(q, VBC_UWD);
              
              // チェック
              if ( cw == 1 ) {
                
                if ( sw==1 ) {
                  if ( !IS_FLUID(bv[mt]) ) {
                    Hostonly_ printf("Error : Fluid cell at Outer boundary Y+ is required\n");
                    Exit(0);
                  }
                }
                else {
                  if ( IS_FLUID(bv[mt]) ) {
                    Hostonly_ printf("Error : Solid cell at Outer boundary Y+ is required\n");
                    Exit(0);
                  }
                }
              }
              
            }// IS_FLUID
          }
        }
      }
      break;
      
    case Z_MINUS:
      if( pn.nID[Z_MINUS] < 0 ){
        for (j=1; j<=jx; j++) {
          for (i=1; i<=ix; i++) {
            m = FBUtility::getFindexS3D(size, guide, i, j, 1  );
            mt= FBUtility::getFindexS3D(size, guide, i, j, 0  );
            
            s = bv[m];
            q = bp[mt];
            
            if ( IS_FLUID(s) ) {
              
              // エンコード処理
              if ( enc_sw ) {
                bv[m]  = s | (OBC_MASK << BC_FACE_B);
              }
              
              // 外部境界で安定化のため，スキームを1次精度にする
              if ( enc_uwd ) bp[mt] = onBit(q, VBC_UWD);
              
              // チェック
              if ( cw == 1 ) {
                
                if ( sw==1 ) {
                  if ( !IS_FLUID(bv[mt]) ) {
                    Hostonly_ printf("Error : Fluid cell at Outer boundary Z- is required\n");
                    Exit(0);
                  }
                }
                else {
                  if ( IS_FLUID(bv[mt]) ) {
                    Hostonly_ printf("Error : Solid cell at Outer boundary Z- is required\n");
                    Exit(0);
                  }
                }
              }
              
            }// IS_FLUID
          }
        }
      }
      break;
      
    case Z_PLUS:
      if( pn.nID[Z_PLUS] < 0 ){
        for (j=1; j<=jx; j++) {
          for (i=1; i<=ix; i++) {
            m = FBUtility::getFindexS3D(size, guide, i, j, kx  );
            mt= FBUtility::getFindexS3D(size, guide, i, j, kx+1);
            
            s = bv[m];
            q = bp[mt];
            
            if ( IS_FLUID(s) ) {
              
              // エンコード処理
              if ( enc_sw ) {
                bv[m]  = s | (OBC_MASK << BC_FACE_T);
              }
              
              // 外部境界で安定化のため，スキームを1次精度にする
              if ( enc_uwd ) bp[mt] = onBit(q, VBC_UWD);
              
              // チェック
              if ( cw == 1 ) {
                
                if ( sw==1 ) {
                  if ( !IS_FLUID(bv[mt]) ) {
                    Hostonly_ printf("Error : Fluid cell at Outer boundary Z+ is required\n");
                    Exit(0);
                  }
                }
                else {
                  if ( IS_FLUID(bv[mt]) ) {
                    Hostonly_ printf("Error : Solid cell at Outer boundary Z+ is required\n");
                    Exit(0);
                  }
                }
              }
              
            }// IS_FLUID
          }
        }
      }
      break;
  } // end of switch
}


/**
 @fn unsigned VoxInfo::fill_cell_edge(int* bid, int* mid, float* cut, const int tgt_id, const int solid_id)
 @brief フィルを実行
 @param bid カット点のID
 @param mid ID配列
 @param cut カット情報
 @param[in] tgt_id サーチするID
 @param[in] solid_id 固体ID
 */
unsigned VoxInfo::fill_cell_edge(int* bid, int* mid, float* cut, const int tgt_id, const int solid_id)
{
  int target = tgt_id;
  int sd = solid_id;
  unsigned m_p, m_e, m_w, m_n, m_s, m_t, m_b;
  int qw, qe, qs, qn, qb, qt;
  unsigned m_sz[3];
  
  int ix = (int)size[0];
  int jx = (int)size[1];
  int kx = (int)size[2];
  m_sz[0] = size[0];
  m_sz[1] = size[1];
  m_sz[2] = size[2];
  unsigned gd = guide;
  int c = 0; /// painted count
  float cpos = 0.9; // なんとなく

#pragma omp parallel for firstprivate(ix, jx, kx, m_sz, gd, target, sd, cpos) \
 private(m_p, m_e, m_w, m_n, m_s, m_t, m_b) \
 private(qw, qe, qs, qn, qb, qt) \
 schedule(static) reduction(+:c)
  
  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        
        m_p = FBUtility::getFindexS3D(m_sz, gd, i,   j,   k  );
        m_e = FBUtility::getFindexS3D(m_sz, gd, i+1, j  , k  );
        m_w = FBUtility::getFindexS3D(m_sz, gd, i-1, j  , k  );
        m_n = FBUtility::getFindexS3D(m_sz, gd, i  , j+1, k  );
        m_s = FBUtility::getFindexS3D(m_sz, gd, i  , j-1, k  );
        m_t = FBUtility::getFindexS3D(m_sz, gd, i  , j  , k+1);
        m_b = FBUtility::getFindexS3D(m_sz, gd, i  , j  , k-1);
        
        // 未ペイントの場合にテスト
        if ( mid[m_p] == 0 ) {
          
          // 隣接セルのカットID >> 0ならばカット無し
          qw = bid[m_w];
          qe = bid[m_e];
          qs = bid[m_s];
          qn = bid[m_n];
          qb = bid[m_b];
          qt = bid[m_t];
          
          // 各方向のテスト
          if ( !((i == 1) && (pn.nID[X_MINUS] < 0)) && (mid[m_w] == target) && (get_BID5(X_MINUS, bid[m_p]) == 0) ) {
            if ( (bid[m_p] == 0) && ( (qs * qn * qb * qt) == 0) ) { // テストセルはカットがなく、隣接セルの少なくともどれかはカットがない場合
              mid[m_p] = target;
              c++;
            }
            else {
              set_BID5(bid[m_w], X_PLUS, sd); // テストする方向からみて、カットIDを設定
              mid[m_p] = sd; // セルIDを固体に変更
              cut[FBUtility::getFindexS3Dcut(m_sz, gd, X_PLUS, i-1, j  , k  )] = cpos; // カット位置をセット
            }
          }
          else if ( !((j == 1) && (pn.nID[Y_MINUS] < 0)) && (mid[m_s] == target) && (get_BID5(Y_MINUS, bid[m_p]) == 0) ) {
            if ( (bid[m_p] == 0) && ( (qw * qe * qb * qt) == 0) ) {
              mid[m_p] = target;
              c++;
            }
            else {
              set_BID5(bid[m_s], Y_PLUS, sd);
              mid[m_p] = sd;
              cut[FBUtility::getFindexS3Dcut(m_sz, gd, Y_PLUS, i  , j-1, k  )] = cpos;
            }
          }
          else if ( !((k == 1) && (pn.nID[Z_MINUS] < 0)) && (mid[m_b] == target) && (get_BID5(Z_MINUS, bid[m_p]) == 0) ) {
            if ( (bid[m_p] == 0) && ( (qw * qe * qs * qn) == 0) ) {
              mid[m_p] = target;
              c++;
            }
            else {
              set_BID5(bid[m_b], Z_PLUS, sd);
              mid[m_p] = sd;
              cut[FBUtility::getFindexS3Dcut(m_sz, gd,  Z_PLUS, i  , j  , k-1)] = cpos;
            }
          }
          else if ( !((k == kx) && (pn.nID[Z_PLUS] < 0)) && (mid[m_t] == target) && (get_BID5(Z_PLUS, bid[m_p]) == 0) ) {
            if ( (bid[m_p] == 0) && ( (qw * qe * qs * qn) == 0) ) {
              mid[m_p] = target;
              c++;
            }
            else {
              set_BID5(bid[m_t], Z_MINUS, sd);
              mid[m_p] = sd;
              cut[FBUtility::getFindexS3Dcut(m_sz, gd, Z_MINUS, i  , j  , k+1)] = cpos;
            }
          }
          else if ( !((j == jx) && (pn.nID[Y_PLUS] < 0)) && (mid[m_n] == target) && (get_BID5(Y_PLUS, bid[m_p]) == 0) ) {
            if ( (bid[m_p] == 0) && ( (qw * qe * qb * qt) == 0) ) {
              mid[m_p] = target;
              c++;
            }
            else {
              set_BID5(bid[m_n], Y_MINUS, sd);
              mid[m_p] = sd;
              cut[FBUtility::getFindexS3Dcut(m_sz, gd, Y_MINUS, i  , j+1, k  )] = cpos;
            }
          }
          else if ( !((i == ix) && (pn.nID[X_PLUS] < 0)) && (mid[m_e] == target) && (get_BID5(X_PLUS, bid[m_p]) == 0) ) {
            if ( (bid[m_p] == 0) && ( (qs * qn * qb * qt) == 0) ) {
              mid[m_p] = target;
              c++;
            }
            else {
              set_BID5(bid[m_e], X_MINUS, sd);
              mid[m_p] = sd;
              cut[FBUtility::getFindexS3Dcut(m_sz, gd, X_MINUS, i+1, j  , k  )] = cpos;
            }
          }
          
        } // target        
      }
    }
  }
  
  return c;
}


/**
 @fn unsigned VoxInfo::fill_cells(int* mid, const int solid_id)
 @brief フィルを実行
 @param mid ID配列
 @param[in] solid_id 固体ID
 */
unsigned VoxInfo::fill_inside(int* mid, const int solid_id)
{
  int sd = solid_id;
  unsigned m_p;
  unsigned m_sz[3];
  
  int ix = (int)size[0];
  int jx = (int)size[1];
  int kx = (int)size[2];
  m_sz[0] = size[0];
  m_sz[1] = size[1];
  m_sz[2] = size[2];
  unsigned gd = guide;
  int c = 0; /// painted count
  
#pragma omp parallel for firstprivate(ix, jx, kx, m_sz, gd, sd) \
 private(m_p) schedule(static) reduction(+:c)
  
  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        
        m_p = FBUtility::getFindexS3D(m_sz, gd, i,   j,   k  );
        
        // 未ペイントの場合
        if ( mid[m_p] == 0 ) {
          mid[m_p] = sd;
          c++;
        }      
      }
    }
  }
  
  return c;
}


/**
 @fn unsigned VoxInfo::fill_isolated_cells(const int* bid, int* mid, const int isolated, const int solid_id)
 @brief フィルを実行
 @param[in] bid カット点のID
 @param mid ID配列
 @param[in] isolated 連結部が孤立したセルの数
 @param[in] solid_id 固体セルのID
 */
void VoxInfo::fill_isolated_cells(const int* bid, int* mid, const int isolated, const int solid_id)
{
  unsigned m_p;
  
  int ix = (int)size[0];
  int jx = (int)size[1];
  int kx = (int)size[2];
  int c = 0; /// count 
  
  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        m_p = FBUtility::getFindexS3D(size, guide, i  , j  , k  );
        
        // テストセルがマークされている場合
        if ( mid[m_p] == -1 ) {
          mid[m_p] = solid_id;
          c++;
        }       
      }
    }
  }
  
  if ( c != isolated ) {
    Hostonly_ printf("Error : The number of cells filled as solid = %d , while isolated cell = %d\n", c, isolated);
    Exit(0);
  }
  
}



/**
 @fn void VoxInfo::find_isolate_Fcell(unsigned order, int* mid, unsigned* bx)
 @brief 孤立した流体セルを探し，周囲の個体媒質で置換，BCindexを修正する
 @param order cmp[]に登録されたMediumListへのエントリ番号
 @param mid ボクセルID配列
 @param bx BCindex ID
 @attention 事前にbx[]の同期が必要 >> 隣接セルがすべて固体の場合をチェックするため
 */
void VoxInfo::find_isolated_Fcell(unsigned order, int* mid, unsigned* bx)
{
  int i, j, k;
  unsigned register s;
  unsigned m_p, m_e, m_w, m_n, m_s, m_t, m_b;
  unsigned key[6];
  int val[6];
  
  int ix = (int)size[0];
  int jx = (int)size[1];
  int kx = (int)size[2];
  
  for (k=1; k<=kx; k++) {
    for (j=1; j<=jx; j++) {
      for (i=1; i<=ix; i++) {
        m_p = FBUtility::getFindexS3D(size, guide, i, j, k);
        s = bx[m_p];
        if ( IS_FLUID(s) ) {
          m_e = FBUtility::getFindexS3D(size, guide, i+1, j  , k  );
          m_w = FBUtility::getFindexS3D(size, guide, i-1, j  , k  );
          m_n = FBUtility::getFindexS3D(size, guide, i  , j+1, k  );
          m_s = FBUtility::getFindexS3D(size, guide, i  , j-1, k  );
          m_t = FBUtility::getFindexS3D(size, guide, i  , j  , k+1);
          m_b = FBUtility::getFindexS3D(size, guide, i  , j  , k-1);
          
          // 隣接セルがすべて固体の場合
          if ( !IS_FLUID(bx[m_e]) && !IS_FLUID(bx[m_w]) &&
               !IS_FLUID(bx[m_n]) && !IS_FLUID(bx[m_s]) &&
               !IS_FLUID(bx[m_t]) && !IS_FLUID(bx[m_b]) ) {
            
            // 最頻値を求める
            for (unsigned l=0; l<6; l++) key[l]=1; // 頻度
            
            val[0] = mid[m_e]; // ID
            val[1] = mid[m_w];
            val[2] = mid[m_n];
            val[3] = mid[m_s];
            val[4] = mid[m_t];
            val[5] = mid[m_b];
            
            for (unsigned l=0; l<6; l++) {
              if ( key[l] != 0 ) {
                for (unsigned m=l+1; m<6; m++) {
                  if ( val[l] == val[m] ) {
                    key[l]++;
                    key[m] = 0;
                  }
                }
              }
            }
            unsigned max_loc=0, cnt=key[0]; // 最頻値とインデクス
            for (unsigned l=1; l<6; l++) {
              if ( key[l] > cnt ) {
                cnt = key[l];
                max_loc = l;
              }
            }
            int max_mid = val[max_loc]; // 最頻値のID
            
            Hostonly_ printf("\n\tReplace isolated fluid cell :: mid(%d,%d,%d)\t original ID=%d >> modified ID=%d (%d/6)\n", i,j,k, mid[m_p], max_mid, cnt);
            
            // 媒質オーダーの変更
            mid[m_p] = max_mid;
            s |= (order << TOP_MATERIAL); // sにMediumListのエントリorderをエンコードする
            
            // 固体セルへ状態を変更　
            bx[m_p] = offBit( s, STATE_BIT );
          }
        }
      }
    }
  }
}

/**
 @fn unsigned VoxInfo::find_mat_odr(unsigned mat_id)
 @brief cmp[]にエンコードされた媒質IDの中から対象となるIDのエントリを探す
 @param mat_id 対象とする媒質ID
 @note 媒質が格納されているオーダーの範囲は[NoBC+1,NoCompo]
 */
unsigned VoxInfo::find_mat_odr(unsigned mat_id)
{
  unsigned odr, id;
  
  id = 0;
  for (unsigned n=NoBC+1; n<=NoCompo; n++) {
    odr = cmp[n].getMatOdr();
    //id  = mat[odr].getMatID();
    if (id == mat_id) return n;
  }
  if (id == 0) {
    Hostonly_ stamped_printf("Error : Material ID [%d] is not listed in CompoList\n", mat_id);
    Exit(0);
  }
  return 0;
}

/**
 @fn void VoxInfo::findVIBCbbox(const int odr, const unsigned* bv, int* st, int* ed)
 @brief VBCのbboxを取得する
 @param odr[in] コンポーネント配列のインデクス
 @param bv[in] BCindex V
 @param st[out] コンポーネントbboxの開始セル
 @param ed[out] コンポーネントbboxの終端セル
 */
void VoxInfo::findVIBCbbox(const int odr, const unsigned* bv, int* st, int* ed)
{
  unsigned s;
  unsigned m;
  bool m_flag;
  int tmp[3];
  
  int ix = (int)size[0];
  int jx = (int)size[1];
  int kx = (int)size[2];

  st[0] = ix;
  st[1] = jx;
  st[2] = kx;
  ed[0] = 0;
  ed[1] = 0;
  ed[2] = 0;
  
  // search
  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        
        m = FBUtility::getFindexS3D(size, guide, i, j, k);
        s = bv[m];
        
        m_flag = false;
        
        // X-
        if ( GET_FACE_BC(s, BC_FACE_W) == odr ) m_flag = true;
        
        // X+
        if ( GET_FACE_BC(s, BC_FACE_E) == odr ) m_flag = true;
        
        // Y-
        if ( GET_FACE_BC(s, BC_FACE_S) == odr ) m_flag = true;
        
        // Y+
        if ( GET_FACE_BC(s, BC_FACE_N) == odr ) m_flag = true;
        
        // Z-
        if ( GET_FACE_BC(s, BC_FACE_B) == odr ) m_flag = true;
        
        // Z+
        if ( GET_FACE_BC(s, BC_FACE_T) == odr ) m_flag = true;
        
        // min, max
        if ( m_flag ) {
          tmp[0] = i;
          tmp[1] = j;
          tmp[2] = k;
          vec3_min(st, st, tmp);
          vec3_max(ed, ed, tmp);
        }
        
      }
    }
  }

}

/**
 @fn unsigned VoxInfo::flip_InActive(unsigned& L, unsigned& G, unsigned id, int* mid, unsigned* bx)
 @brief BCindexにActive/Inactiveをエンコード
 @retval 不活性化した数
 @param L[out] ノードローカルの不活性セル数
 @param G[out] グローバルの不活性セル数
 @param id セルID
 @param mid ボクセル配列
 @param bx BCindex ID
 */
unsigned long VoxInfo::flip_InActive(unsigned long& L, 
                                     unsigned long& G, 
                                     const unsigned id, 
                                     const int* mid, 
                                     unsigned* bx)
{
  int i,j,k, c_p;
  unsigned s, m;
  unsigned long c=0;
  
  int idd = (int)id;
  
  int ix = (int)size[0];
  int jx = (int)size[1];
  int kx = (int)size[2];
  int gd = (int)guide;
  
  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        m = FBUtility::getFindexS3D(size, guide, i, j, k);
        c_p = mid[m];
        s = bx[m];
        
        if ( c_p == idd ) {
          if ( BIT_IS_SHIFT(s, ACTIVE_BIT) ) { // 活性化してある場合に，不活性にし，カウント
            s = offBit( s, ACTIVE_BIT );
            bx[m] = s;
            c++;
          }
        }
        
      }
    }
  }
  
  L = c;
  
  if ( pn.numProc > 1 ) {
    unsigned long tmp = c;
    MPI_Allreduce(&tmp, &c, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  }
  
  G = c;
  
  return c;
}


/**
 @fn void VoxInfo::get_Compo_Area_Cut(unsigned n, PolylibNS::MPIPolylib* PL)
 @brief コンポーネントの断面積を求める
 @param n エントリ番号
 @param bd BCindex ID
 @param bv BCindex V
 @param bh1 BCindex H1
 @param PL MPIPolylibクラス
 @note ポリゴン情報は有次元とする
 */
void VoxInfo::get_Compo_Area_Cut(unsigned n, PolylibNS::MPIPolylib* PL)
{
  using namespace PolylibNS;
  
  REAL_TYPE a;
  
  int id  = cmp[n].getMatOdr();
  int def = cmp[n].getDef();
  
  vector<PolygonGroup*>* pg_roots = PL->get_root_groups();
  vector<PolygonGroup*>::iterator it;
  
  switch ( cmp[n].getType() ) {
    case SPEC_VEL:
    case SPEC_VEL_WH:
    case OUTFLOW:
      
      printf("\t  ID : Polygon Group : Area (m*m)\n");
      
      for (it = pg_roots->begin(); it != pg_roots->end(); it++) {
        std::string m_pg = (*it)->get_name();
        int m_id = (*it)->get_id();
        cmp[n].area = (*it)->get_group_area();
        //printf("\t %3d : %s : %e : %d\n", m_id, m_pg.c_str(), cmp[n].area, (*it)->get_group_num_tria());
      }
      
      return;
      break;
      
      
    case HEATFLUX:
    case TRANSFER:
    case ISOTHERMAL:
      break;
      
    case HEX:
    case FAN:
      return;
      break;
      
    case CELL_MONITOR:
    case DARCY:
      break;
  }
  
  delete pg_roots;
}



/**
 @fn void VoxInfo::getOffset(int* st, int* ofst)
 @brief インデクス範囲を決めるオフセットを計算する
 @param st コンポーネントの開始インデクス
 @param ofst オフセット量
 @note
 - 計算領域内部のノード間境界部でのカウント重複を防ぐため，開始点位置の判断を行う
 - 外部境界に接する面では，補正せず，常にoffset=1
 - 内部境界面では，重複しないようにoffset=0
 */
void VoxInfo::getOffset(int* st, int* ofst)
{
  if( pn.nID[X_MINUS] < 0 ){
    ofst[0] = 1;
  }
  else {
    ofst[0] = ( st[0] == 1 ) ? 0 : 1;
  }
  
  if( pn.nID[Y_MINUS] < 0 ){
    ofst[1] = 1;
  }
  else {
    ofst[1] = ( st[1] == 1 ) ? 0 : 1;
  }
  
  if( pn.nID[Z_MINUS] < 0 ){
    ofst[2] = 1;
  }
  else {
    ofst[2] = ( st[2] == 1 ) ? 0 : 1;
  }
}



/**
 @fn void VoxInfo::getQuadrant(unsigned* q, REAL_TYPE t1, REAL_TYPE t2)
 @brief 局所座標(t1, t2)がどの象限になるかを調べる
 @param q 象限番号，戻り値
 @param t1 第一座標
 @param t2 第二座標
 */
void VoxInfo::getQuadrant(unsigned* q, REAL_TYPE t1, REAL_TYPE t2)
{
  if (t1 == 0.0) { // 線上は両方の象限に加える
    if (t2 == 0.0) {
      q[1]++;
      q[2]++;
      q[3]++;
      q[4]++;
    }
    else if (t2 > 0.0) {
      q[1]++;
      q[2]++;
    }
    else { // t2<0.0
      q[3]++;
      q[4]++;
    }
  }
  else if (t1 > 0.0) {
    if (t2 == 0.0) {
      q[1]++;
      q[4]++;
    }
    else if (t2 > 0.0) {
      q[1]++;
    }
    else { // t2<0.0
      q[4]++;
    }
  }
  else { // t1<0.0
    if (t2 == 0.0) {
      q[2]++;
      q[3]++;
    }
    else if (t2 > 0.0) {
      q[2]++;
    }
    else { // t2<0.0
      q[3]++;
    }
  }
}


// CPMクラスポインタのコピー
void VoxInfo::importCPM(cpm_ParaManager* m_paraMngr)
{
  if ( !m_paraMngr ) Exit(0);
  paraMngr = m_paraMngr;
}



/**
 @fn bool VoxInfo::paint_first_seed(int* mid, const int* idx, const int target)
 @brief シード点をペイントする
 @param mid ID配列
 @param[in] idx seed点のインデクス
 @param[in] target ペイントするID
 */
bool VoxInfo::paint_first_seed(int* mid, const int* idx, const int target)
{
  unsigned m_p;
  
  m_p = FBUtility::getFindexS3D(size, guide, idx[0], idx[1], idx[2]);
  if ( mid[m_p] != 0 ) {
    return false;
  }
  else {
    mid[m_p] = target;
  }
  
  return true;
}


/**
 @fn void VoxInfo::printScanedCell(FILE* fp)
 @brief ボクセルをスキャンしたIDの数と境界条件の数，含まれるIDのリストを表示する
 @param fp ファイル出力のファイルポインタ
 */
void VoxInfo::printScanedCell(FILE* fp)
{
  for (int i=1; i<=NoVoxID; i++) {
    fprintf(fp,"\t\t%3d : ID = [%d]\n", i, colorList[i]);
    if ( colorList[i] == 0 ) {
      Hostonly_ stamped_fprintf(fp, "\t*** Voxel includes ID=0, which is invalid.  ***\n\n");
      Exit(0);
    }
  }
}


/**
 @fn int VoxInfo::scanCell(int *cell, const int* cid, const unsigned ID_replace)
 @brief cellで保持されるボクセルid配列をスキャンし，coloList[]に登録する
 @retval 含まれるセルID数
 @param cell ボクセルIDを保持する配列
 @param cid セルIDリスト 
 @param ID_replace ID[0]を置換するID
 @note 重複がないようにcolorList[]に登録する
 */ 
int VoxInfo::scanCell(int *cell, const int* cid, const unsigned ID_replace)
{
  int target;
  unsigned m;
  
  int ix = (int)size[0];
  int jx = (int)size[1];
  int kx = (int)size[2];
  int gd = (int)guide;
  
  // ID[0]を置換するオプションが指定されている場合（ID_replaceに値がセット）
  if ( ID_replace != 0 ) {
    target = (int)ID_replace;
    Hostonly_ printf("\n\tID[0] is replaced by ID[%d]\n", target);
    for (int k=1; k<=kx; k++) {
      for (int j=1; j<=jx; j++) {
        for (int i=1; i<=ix; i++) {
          m = FBUtility::getFindexS3D(size, guide, i, j, k);
          if ( cell[m] == 0 ) cell[m] = target;
        }
      }
    }
  }
  
  // 内部領域に対して，マイナスとゼロをチェック
  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        target = cell[FBUtility::getFindexS3D(size, guide, i, j, k)];
        if ( target<=0 ) {
          Hostonly_ stamped_printf("\tVoxel data includes non-positive ID [%d] at (%d, %d, %d)\n", target, i, j, k);
          Exit(0);
        }
      }
    }
  }
  
  // 内部領域の媒質IDに対して、カウント
  
  // colorSet[] ローカルなIDのカウント >> unsignedは超えないはず
  unsigned colorSet[MODEL_ID_MAX+1];
  for (int i=0; i<MODEL_ID_MAX+1; i++) colorSet[i]=0;
  
  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        m = FBUtility::getFindexS3D(size, guide, i, j, k);
        target = cell[m];
        
        if ( target <= MODEL_ID_MAX ) {
          if ( colorSet[target]++ > UINT_MAX ) {
            Hostonly_ stamped_printf("\tError : count included in Model exceeds UINT_MAX\n");
            Exit(0);
          }
        }
        else {
          Hostonly_ stamped_printf("\tError : ID included in Model is greater than %d.\n", MODEL_ID_MAX);
          Exit(0);
        }
        
      }
    }
  }
  
  
  // 外部領域の媒質IDをcolorSetに追加する
  for (int i=0; i<NOFACE; i++) {
    target = cid[i];
    colorSet[target]++;
  }
  
  // 集約時の桁あふれを回避するため、1に規格化
  for (int i=0; i<MODEL_ID_MAX+1; i++) {
    if ( colorSet[i] != 0 ) {
      colorSet[i] = 1;
    }
  }

	// colorSet[] の集約
  if ( pn.numProc > 1 ) {
		
    unsigned clist[MODEL_ID_MAX+1];
    for (int i=0; i<MODEL_ID_MAX+1; i++) clist[i]=0;
    
    for (int i=0; i<MODEL_ID_MAX+1; i++) clist[i] = colorSet[i];
    
    MPI_Allreduce(clist, colorSet, MODEL_ID_MAX+1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
    
  }
  // この時点で、存在するIDの数はpn.numProc個になっている >> unsgined の範囲内

  // colorList[]へ詰めてコピー colorList[0]は不使用
  int b=1;
  for (int i=0; i<MODEL_ID_MAX+1; i++) {
    if ( colorSet[i] != 0 ) {
      colorList[b] = i;
      b++;
    }
  }
  
  NoVoxID = b-1;
  
  return NoVoxID;
}


/**
 @fn void VoxInfo::setNoCompo_BC(unsigned m_NoBC, unsigned m_NoCompo)
 @brief コンポーネントの操作に必要な定数の設定
 @param m_NoBC 境界条件数
 @param m_NoCompo コンポーネントの総数
 */
void VoxInfo::setNoCompo_BC(unsigned m_NoBC, unsigned m_NoCompo)
{
  NoBC    = m_NoBC;
  NoCompo = m_NoCompo;
}

/**
 @fn void VoxInfo::setAmask_InActive(unsigned id, int* mid, unsigned* bh)
 @brief BCindexに不活性セルに対する断熱マスクをエンコード
 @param id セルID
 @param mid ボクセル配列
 @param bh BCindex H2
 */
void VoxInfo::setAmask_InActive(unsigned id, int* mid, unsigned* bh)
{
  int i,j,k;
  unsigned s;
  unsigned m_p, m_e, m_w, m_n, m_s, m_t, m_b;
  int      c_p, c_e, c_w, c_n, c_s, c_t, c_b;
  int idd = (int)id;
  
  int ix = (int)size[0];
  int jx = (int)size[1];
  int kx = (int)size[2];
  int gd = (int)guide;
  
  for (k=1; k<=kx; k++) {
    for (j=1; j<=jx; j++) {
      for (i=1; i<=ix; i++) {
        m_p = FBUtility::getFindexS3D(size, guide, i  , j  , k  );
        m_e = FBUtility::getFindexS3D(size, guide, i+1, j  , k  );
        m_w = FBUtility::getFindexS3D(size, guide, i-1, j  , k  );
        m_n = FBUtility::getFindexS3D(size, guide, i  , j+1, k  );
        m_s = FBUtility::getFindexS3D(size, guide, i  , j-1, k  );
        m_t = FBUtility::getFindexS3D(size, guide, i  , j  , k+1);
        m_b = FBUtility::getFindexS3D(size, guide, i  , j  , k-1);
        
        c_p = mid[m_p];
        c_e = mid[m_e];
        c_w = mid[m_w];
        c_n = mid[m_n];
        c_s = mid[m_s];
        c_t = mid[m_t];
        c_b = mid[m_b];
        
        s = bh[m_p];
        
        // 流体セル，固体セルに関わらず，隣接セルがInactive指定の場合には断熱にする
        if ( c_p != idd ) { // ただし，自セルが指定ID以外のとき
          if ( c_w == idd ) s = offBit( s, ADIABATIC_W );
          if ( c_e == idd ) s = offBit( s, ADIABATIC_E );
          if ( c_s == idd ) s = offBit( s, ADIABATIC_S );
          if ( c_n == idd ) s = offBit( s, ADIABATIC_N );
          if ( c_b == idd ) s = offBit( s, ADIABATIC_B );
          if ( c_t == idd ) s = offBit( s, ADIABATIC_T );
        }
        bh[m_p] = s;
      }
    }
  }
}



/**
 @fn void VoxInfo::setAmask_Solid(unsigned* bh)
 @brief KOSがSOLID_CONDUCTIONの場合の断熱マスクの処理
 @param bh BCindex H2
 @note 
 - S-F面のS側を断熱にする
 */
void VoxInfo::setAmask_Solid(unsigned* bh)
{
  int i,j,k;
  unsigned m_p, m_e, m_w, m_n, m_s, m_t, m_b;
  unsigned register s, s_e, s_w, s_n, s_s, s_t, s_b;
  
  int ix = (int)size[0];
  int jx = (int)size[1];
  int kx = (int)size[2];
  int gd = (int)guide;
  
  for (k=1; k<=kx; k++) {
    for (j=1; j<=jx; j++) {
      for (i=1; i<=ix; i++) {
        m_p = FBUtility::getFindexS3D(size, guide, i  , j  , k  );
        s   = bh[m_p];
        
        if ( !IS_FLUID(s) ) { // pセルが固体セルの場合が対象
          m_e = FBUtility::getFindexS3D(size, guide, i+1, j  , k  );
          m_w = FBUtility::getFindexS3D(size, guide, i-1, j  , k  );
          m_n = FBUtility::getFindexS3D(size, guide, i  , j+1, k  );
          m_s = FBUtility::getFindexS3D(size, guide, i  , j-1, k  );
          m_t = FBUtility::getFindexS3D(size, guide, i  , j  , k+1);
          m_b = FBUtility::getFindexS3D(size, guide, i  , j  , k-1);
          
          s_e = bh[m_e];
          s_w = bh[m_w];
          s_n = bh[m_n];
          s_s = bh[m_s];
          s_t = bh[m_t];
          s_b = bh[m_b];
          
          // X-
          if ( IS_FLUID(s_w) ) { // S-F界面であること
            s = offBit( s, ADIABATIC_W );
          }
          
          // X+
          if ( IS_FLUID(s_e) ) {
            s = offBit( s, ADIABATIC_E );
          }
          
          // Y-
          if ( IS_FLUID(s_s) ) {
            s = offBit( s, ADIABATIC_S );
          }
          
          // Y+
          if ( IS_FLUID(s_n) ) {
            s = offBit( s, ADIABATIC_N );
          }
          
          // Z-
          if ( IS_FLUID(s_b) ) {
            s = offBit( s, ADIABATIC_B );
          }
          
          // Z+
          if ( IS_FLUID(s_t) ) {
            s = offBit( s, ADIABATIC_T );
          }
        }
        bh[m_p] = s;
      }
    }
  }
}

/**
 @fn void VoxInfo::setAmask_Thermal(unsigned* bh)
 @brief KOSがTHERMAL_FLOW, THERMAL_FLOW_NATURALの場合の断熱マスクの処理
 @param bh BCindex H2
 @note 
 - S-F面のF側を断熱にする
 */
void VoxInfo::setAmask_Thermal(unsigned* bh)
{
  int i,j,k;
  unsigned m_p, m_e, m_w, m_n, m_s, m_t, m_b;
  unsigned register s, s_e, s_w, s_n, s_s, s_t, s_b;
  
  int ix = (int)size[0];
  int jx = (int)size[1];
  int kx = (int)size[2];
  int gd = (int)guide;
  
  for (k=1; k<=kx; k++) {
    for (j=1; j<=jx; j++) {
      for (i=1; i<=ix; i++) {
        m_p = FBUtility::getFindexS3D(size, guide, i  , j  , k  );
        s   = bh[m_p];
        
        if ( IS_FLUID(s) ) { // pセルが流体セルの場合が対象
          m_e = FBUtility::getFindexS3D(size, guide, i+1, j  , k  );
          m_w = FBUtility::getFindexS3D(size, guide, i-1, j  , k  );
          m_n = FBUtility::getFindexS3D(size, guide, i  , j+1, k  );
          m_s = FBUtility::getFindexS3D(size, guide, i  , j-1, k  );
          m_t = FBUtility::getFindexS3D(size, guide, i  , j  , k+1);
          m_b = FBUtility::getFindexS3D(size, guide, i  , j  , k-1);
          
          s_e = bh[m_e];
          s_w = bh[m_w];
          s_n = bh[m_n];
          s_s = bh[m_s];
          s_t = bh[m_t];
          s_b = bh[m_b];
          
          // X-
          if ( !IS_FLUID(s_w) ) { // S-F界面であること
            s = offBit( s, ADIABATIC_W );
          }
          
          // X+
          if ( !IS_FLUID(s_e) ) {
            s = offBit( s, ADIABATIC_E );
          }
          
          // Y-
          if ( !IS_FLUID(s_s) ) {
            s = offBit( s, ADIABATIC_S );
          }
          
          // Y+
          if ( !IS_FLUID(s_n) ) {
            s = offBit( s, ADIABATIC_N );
          }
          
          // Z-
          if ( !IS_FLUID(s_b) ) {
            s = offBit( s, ADIABATIC_B );
          }
          
          // Z+
          if ( !IS_FLUID(s_t) ) {
            s = offBit( s, ADIABATIC_T );
          }
        }
        bh[m_p] = s;
      }
    }
  }
}



/**
 @fn void VoxInfo::setBCIndex_base1(unsigned* bx, int* mid, float* cvf)
 @brief bx[]に各境界条件の共通のビット情報をエンコードする（その1）
 @param bx BCindex ID
 @param mid ID配列
 @param cvf コンポーネントの体積率
 @note 事前に，cmp[]へMediumListへのエントリ番号をエンコードしておく -> cmp[].setMatOdr()
 */
void VoxInfo::setBCIndex_base1(unsigned* bx, int* mid, float* cvf)
{
  unsigned odr;
  int id;
  unsigned register s;
  int md;
  
  size_t nx = (size[0]+2*guide) * (size[1]+2*guide) * (size[2]+2*guide); // ガイドセルを含む全領域を対象にする
  
  // セルの状態を流体で初期化
  for (unsigned m=0; m<nx; m++) {
    bx[m] = onBit( bx[m], STATE_BIT );
  }
  
  // セルIDをエンコード
  for (unsigned m=0; m<nx; m++) {
    md = mid[m];
    bx[m] |= ((unsigned)md << TOP_CELL_ID) ;
  }

  // 体積率コンポーネントのみ
  int st[3], ed[3];
  for (unsigned n=1; n<=NoBC; n++) {
    id = cmp[n].getMatOdr();
    
    // 体積率コンポーネントのオーバーラップは考慮していないので、エラーが生じる可能性あり
    unsigned m;
    switch ( cmp[n].getType() ) {
      case HEX:
      case FAN:
        cmp[n].getBbox(st, ed);

        for (int k=st[2]; k<=ed[2]; k++) {
          for (int j=st[1]; j<=ed[1]; j++) {
            for (int i=st[0]; i<=ed[0]; i++) {
              m = FBUtility::getFindexS3D(size, guide, i, j, k);
              if ( (cvf[m] > 0.0) && IS_FLUID(bx[m]) ) { // 流体でコンポーネントの体積率があるとき
                bx[m] |= (id << TOP_CELL_ID) ;
                mid[m] = id;
              }
            }
          }
        }
        break;
        
      case CELL_MONITOR:
        cmp[n].getBbox(st, ed);
        
        for (int k=st[2]; k<=ed[2]; k++) {
          for (int j=st[1]; j<=ed[1]; j++) {
            for (int i=st[0]; i<=ed[0]; i++) {
              m = FBUtility::getFindexS3D(size, guide, i, j, k);
              if ( mid[m] == (int)id ) {
                bx[m] |= (id << TOP_CELL_ID) ;
              }
            }
          }
        }
        break;
    }
  }

  // MediumListのエントリをエンコードする
  for (unsigned n=1; n<=NoCompo; n++) {
    id  = cmp[n].getMatOdr();

    for (size_t m=0; m<nx; m++) {
      if ( mid[m] == id ) bx[m] |= (id << TOP_MATERIAL);
    }
  }

  // CompoListのエントリ　セル要素のみ
  for (unsigned n=1; n<=NoBC; n++) {
    id = cmp[n].getMatOdr();
    
    switch ( cmp[n].getType() ) {
      case PERIODIC:
      case CELL_MONITOR:
      case IBM_DF:
      case DARCY:
      case HEAT_SRC: // 熱のBCもここで処理しておく
      case CNST_TEMP:
      case HEX:
      case FAN:
        cmp[n].setElement( encodeOrder(n, id, mid, bx) );
        break;
    }
  }
  
  
  // 状態のエンコード
  for (size_t m=0; m<nx; m++) {
    s = bx[m];
    odr = DECODE_MAT( s );
    if ( mat[odr].getState() == FLUID ) {
      s = onBit( s, STATE_BIT );
    }
    else {  // SOLID
      s = offBit( s, STATE_BIT );
    }
    bx[m] = s;
  }
  
}

/**
 @fn void VoxInfo::setBCIndex_base2(unsigned* bx, int* mid, SetBC* BC, unsigned& Lcell, unsigned& Gcell, unsigned KOS)
 @brief bx[]に各境界条件の共通のビット情報をエンコードする（その2）
 @param bx BCindex ID
 @param mid ID配列
 @param BC SetBCクラスのポインタ
 @param Lcell ノードローカルの有効セル数
 @param Gcell グローバルの有効セル数
 @param KOS 解くべき方程式の種類 KIND_OF_SOLVER
 @note 事前に，cmp[]へMediumListへのエントリをエンコードしておく -> cmp[].setMatOdr()
 */
void VoxInfo::setBCIndex_base2(unsigned* bx, int* mid, SetBC* BC, unsigned long & Lcell, unsigned long & Gcell, const unsigned KOS)
{
  unsigned mat_id;
  int id;
  
  // 孤立した流体セルの属性変更
  for (unsigned n=1; n<=NoCompo; n++) {
    find_isolated_Fcell( cmp[n].getMatOdr(), mid, bx);
  }
  
  // BCIndexにそのセルが計算に有効(active)かどうかをエンコードする．KindOfSolverによって異なる
  encActive(Lcell, Gcell, bx, KOS);
  
  // Inactive指定のセルを不活性にする
  unsigned long m_L, m_G;
  for (unsigned n=1; n<=NoBC; n++) {
    id = cmp[n].getMatOdr();
    
    if ( cmp[n].getType() == INACTIVE ) {
      
      m_L = m_G = 0;
      cmp[n].setElement( flip_InActive(m_L, m_G, id, mid, bx) );
      Lcell -= m_L;
      Gcell -= m_G;
    }
  }
  
  // コンポーネントに登録された媒質のセル数を数え，elementにセットする
  for (unsigned n=NoBC+1; n<=NoCompo; n++) {
    id = cmp[n].getMatOdr();
    cmp[n].setElement( countState(id, mid) );
  }
  
  
}

/**
 @fn void VoxInfo::setBCIndexH(unsigned* bcd, unsigned* bh1, unsigned* bh2, int* mid, SetBC* BC, unsigned kos)
 @brief 境界条件のビット情報をエンコードする
 @param bcd BCindex ID
 @param bh1 BCindex H1
 @param bh1 BCindex H2
 @param mid ID配列
 @param BC SetBCクラスのポインタ
 @param kos KindOfSolver
 */
void VoxInfo::setBCIndexH(unsigned* bcd, unsigned* bh1, unsigned* bh2, int* mid, SetBC* BC, unsigned kos)
{
  unsigned n;
  int id;
  int deface;
  int i, j, k;
  unsigned m;
  
  int ix = (int)size[0];
  int jx = (int)size[1];
  int kx = (int)size[2];
  int gd = (int)guide;
  
  // 断熱マスクを非断熱(1)に初期化する
  for (k=0; k<=kx+1; k++) {
    for (j=0; j<=jx+1; j++) {
      for (i=0; i<=ix+1; i++) {
        m = FBUtility::getFindexS3D(size, guide, i  , j  , k  );
        bh2[m] |= ( 0x3f << ADIABATIC_W ); // 6bitまとめて初期化
      }
    }
  }
  
  // THERMAL_FLOW, THERMAL_FLOW_NATURAL, SOLID_CONDUCTIONのときに，デフォルトとしてSolid-Fluid面を断熱にする
  switch (kos) {
    case THERMAL_FLOW:
    case THERMAL_FLOW_NATURAL:
      setAmask_Thermal(bh2);
      break;
      
    case SOLID_CONDUCTION:
      setAmask_Solid(bh2);
      break;
  }
  
  // 対称境界面に断熱マスクをセット
  for (int face=0; face<NOFACE; face++) {
    if ( BC->export_OBC(face)->get_Class() == OBC_SYMMETRIC ) {
      encAmask_SymtrcBC(face, bh2);
    }
  }
  
  // 不活性セルの場合の断熱マスク処理
  for (unsigned n=1; n<=NoBC; n++) {
    if ( cmp[n].getType() == INACTIVE ) {
      setAmask_InActive(cmp[n].getMatOdr(), mid, bh2);
    }
  }
  
  // bh2の下位5ビットにはBCのエントリのみ(1~31)エンコード
  for (n=1; n<=NoBC; n++) {
    id = cmp[n].getMatOdr();
    deface = cmp[n].getDef();
    
    switch ( cmp[n].getType() ) {
      case ADIABATIC:
        cmp[n].setElement( encQface(n, id, mid, bcd, bh1, bh2, deface, false) ); // 断熱ビット(0)
        break;
        
      case HEATFLUX:
        cmp[n].setElement( encQface(n, id, mid, bcd, bh1, bh2, deface, true) ); // 断熱ビット(1)
        break;
        
      case SPEC_VEL_WH: // 要素数については，setBCIndexV()でカウントしているので，不要
      case OUTFLOW:
        encQfaceSVO(n, id, mid, bcd, bh1, bh2, deface);
        // ?? setInactive_Compo(id, cmp[n].getDef(), mid, bh1, bh2); // 不活性セルの指定
        break;
        
      case TRANSFER:
        switch ( cmp[n].getHtype() ) {
          case HT_N:
            //
            break;
            
          case HT_S:
          case HT_SN:
          case HT_SF:
            if ( (kos == CONJUGATE_HEAT_TRANSFER) || (kos == SOLID_CONDUCTION) ) {
              Hostonly_ printf("\tHeat Transfer(_S, _SF, _SN) can be specified only in the case of 'THERMAL_FLOW' or 'THERMAL_FLOW_NATURAL'\n");
              Exit(0);
            }
            else {
              cmp[n].setElement( encQfaceHT_S(n, id, mid, bcd, bh1, bh2, deface) );
            }
            break;
            
          case HT_B:
            if ( (kos == THERMAL_FLOW) || (kos == THERMAL_FLOW_NATURAL) ) {
              Hostonly_ printf("\tHeat Transfer(_B) can be specified only in the case of 'CONJUGATE_HEAT_TRANSFER' or 'SOLID_CONDUCTION'\n");
              Exit(0);
            }
            else {
              cmp[n].setElement( encQfaceHT_B(n, id, mid, bcd, bh1, bh2, deface) );
            }
            break;
        }
        break;
        
      case ISOTHERMAL:
        if ( (kos == THERMAL_FLOW) || (kos == THERMAL_FLOW_NATURAL) ) {
          cmp[n].setElement( encQfaceISO_SF(n, id, mid, bcd, bh1, bh2, deface) );
          
        }
        else {
          cmp[n].setElement( encQfaceISO_SS(n, id, mid, bcd, bh1, bh2, deface) );
        }
        break;
        
      case RADIANT:
        break;
        
        // Q BC at Volume; idのガイドセルチェックなし
      case HEAT_SRC:
      case CNST_TEMP:
        cmp[n].setElement( encodeOrder(n, id, mid, bh2) );
        break;
    }
    
  }// end loop
  
  // set gamma coef. for Heat BC
  encHbit(bh1, bh2);
}


/**
 @fn unsigned VoxInfo::setBCIndexP(unsigned* bcd, unsigned* bcp, int* mid, SetBC* BC, bool isCDS, float* cut)
 @brief 圧力境界条件のビット情報をエンコードする
 @param bcd BCindex ID
 @param bcp BCindex P
 @param mid ID配列
 @param BC SetBCクラスのポインタ
 @param isCDS CDS->true
 @param cut 距離情報
 @retval 表面セル数
 */
unsigned VoxInfo::setBCIndexP(unsigned* bcd, unsigned* bcp, int* mid, SetBC* BC, bool isCDS, float* cut)
{

  unsigned surface = 0;

  // 初期化 @note ビットを1に初期化する．初期化範囲はガイドセルを含む全領域．セルフェイスの射影処理で必要．
  unsigned mx = (size[0]+2*guide) * (size[1]+2*guide) * (size[2]+2*guide);
  for (unsigned m=0; m<mx; m++) {
    bcp[m] |= ( 0x3ffff << BC_NDAG_W ); // BC_NDAG_W〜BC_D_Tまで18bitまとめて1に初期化
  }

  // 計算領域内の壁面のNeumannBCのマスク処理と固体に隣接するFセルに方向フラグをエンコードし，表面セル数を返す
  if ( !isCDS ) {
    surface = encPbit_N_Binary(bcp);    // Binary
  }
  else {
    surface = encPbit_N_Cut(bcp, cut, true);  // Cut-Distance
  }

  // 外部境界のビットフラグをエンコード
  BoundaryOuter* m_obc=NULL;
  unsigned F;
  
  for (int face=0; face<NOFACE; face++) {
    m_obc = BC->export_OBC(face);
    F = m_obc->get_Class();
    
    switch ( F ) {
      case OBC_WALL:
      case OBC_SYMMETRIC:
        encPbit_OBC(face, bcp, "Neumann", true);
        break;
        
      case OBC_SPEC_VEL:
        encPbit_OBC(face, bcp, "Neumann", false);
        break;
        
      case OBC_TRC_FREE:
        encPbit_OBC(face, bcp, "Dirichlet", false); // test
        break;
        
      case OBC_FAR_FIELD:
        encPbit_OBC(face, bcp, "Neumann", false);
        break;
        
      case OBC_OUTFLOW:
        if ( m_obc->get_pType() == P_DIRICHLET ) {
          encPbit_OBC(face, bcp, "Dirichlet", false);
        }
        else {
          encPbit_OBC(face, bcp, "Neumann", false);
        }
        break;
        
        // 外部境界値を指定
      case OBC_PERIODIC:
        break;
        
      default:
        Exit(0);
    }
  }

  // 内部境界のコンポーネントのエンコード
  unsigned n;
  int id;
  int deface;
  
  for (n=1; n<=NoBC; n++) {
    id = cmp[n].getMatOdr();
    deface = cmp[n].getDef();
    
    switch ( cmp[n].getType() ) {
      case SPEC_VEL:
      case SPEC_VEL_WH:
      case OUTFLOW:
        cmp[n].setElement( encPbit_N_IBC(n, id, mid, bcd, bcp, deface) );
        break;
    }
  }

  // 全周Neumannフラグのセルと排他性をチェックし，反復行列の非対角要素/対角要素をエンコードする
  encPbit(bcp);

  // debug
#if 0
  int i,j,k;
  unsigned m;
  float w, q;
  i=1;
  for (k=1; k<=kx; k++) {
    for (j=1; j<=jx; j++) {
      m = FBUtility::getFindexS3D(size, guide, i, j, k);
      w = GET_SHIFT_F(bcp[m], BC_N_W);
      q = GET_SHIFT_F(bcp[m], BC_NDAG_W);
      printf("(1, %3d, %3d) N=%f Coef_W=%f\n", j,k,w,q);
    }
  }
#endif

  return surface;
}

/**
 @fn void VoxInfo::setBCIndexV(unsigned* bv, int* mid, SetBC* BC, unsigned* bp, bool isCDS, float* cut, int* cut_id)
 @brief bv[]に境界条件のビット情報をエンコードする
 @param bv BCindex V
 @param mid ID配列
 @param BC SetBCクラスのポインタ
 @param bp BCindex P
 @param isCDS CDS->true
 @param cut 距離情報
 @param cut_id カット点ID
 */
void VoxInfo::setBCIndexV(unsigned* bv, int* mid, SetBC* BC, unsigned* bp, bool isCDS, float* cut, int* cut_id)
{
  // ガイドセルの媒質情報をチェックし，流束形式のBCの場合にビットフラグをセット
  BoundaryOuter* m_obc=NULL;
  unsigned F;
  
  // 外部境界
  for (int face=0; face<NOFACE; face++) {
    m_obc = BC->export_OBC(face);
    F = m_obc->get_Class();
    
    switch ( F ) {
      case OBC_SYMMETRIC:
      case OBC_WALL:
        encVbit_OBC(face, bv, "solid", true, "check", bp, false); // 流束形式
        break;
        
      case OBC_SPEC_VEL:
      case OBC_OUTFLOW:
        encVbit_OBC(face, bv, "fluid", true, "check", bp, false); // 流束形式
        break;
        
      case OBC_TRC_FREE:
        encVbit_OBC(face, bv, "fluid", false, "check", bp, true); // 境界値指定
        break;
        
      case OBC_FAR_FIELD:
        encVbit_OBC(face, bv, "fluid", false, "check", bp, true); // 境界値指定
        break;
        
      case OBC_PERIODIC:
        encVbit_OBC(face, bv, "fluid", false, "no_check", bp, false); // 境界値指定，内部セルの状態をコピーするので，ガイドセル状態のチェックなし
        break;
        
      default:
        Exit(0);
    }
    
    // 有効セル数をカウントし，集約
    m_obc->set_ValidCell( count_ValidCell_OBC(face, bv) );
  }
  
  // 内部境界のコンポーネントのエンコード
  unsigned n, m_dir;
  int deface, id;
  float vec[3];
  
  for (n=1; n<=NoBC; n++) {
    id = cmp[n].getMatOdr();
    deface = cmp[n].getDef();
    
    vec[0] = (float)cmp[n].nv[0];
    vec[1] = (float)cmp[n].nv[1];
    vec[2] = (float)cmp[n].nv[2];
    
    m_dir = cmp[n].getBClocation();
    
    switch ( cmp[n].getType() ) {
      case SPEC_VEL:
      case SPEC_VEL_WH:
      case OUTFLOW:
        if (isCDS) { // cut
          cmp[n].setElement( encVbit_IBC_Cut(n, id, bv, bp, cut, cut_id, vec, m_dir) );
        }
        else { // binary
          cmp[n].setElement( encVbit_IBC(n, id, mid, bv, deface, bp) ); // bdへのエンコードはsetBCindexP()で済み
        }
        break;
    }    
  }
  
}


/**
 @fn void VoxInfo::setCmpFraction(CompoList* cmp, unsigned* bx, float* vf)
 @brief bx[]のコンポーネントエントリを参照して体積率を計算し，圧力損失コンポーネントの場合にはビットを立てる
 @param cmp コンポーネントリスト
 @param bx BCindex ID
 @param vf 体積率
 */
void VoxInfo::setCmpFraction(CompoList* cmp, unsigned* bx, float* vf)
{
	unsigned m, s;
  int i,j,k;
  int st[3], ed[3];
  int f;
  
  for (unsigned n=1; n<=NoBC; n++) {
    if ( cmp[n].isVFraction() ) { // 対象のコンポーネント
      cmp[n].getBbox(st, ed);
      
      for (k=st[2]; k<=ed[2]; k++) {
        for (j=st[1]; j<=ed[1]; j++) {
          for (i=st[0]; i<=ed[0]; i++) {
            m = FBUtility::getFindexS3D(size, guide, i, j, k);
            s = bx[m];
            if ( ( s & MASK_6) == n ) {
              f = (int)floorf(vf[m]*255.0 + 0.5); // 0.0<vf<1.0 の四捨五入 > 8ビットで量子化
              if ( f > 255 ) assert(0);
              s |= (f << TOP_VF); // 8ビットで量子化
              if ( cmp[n].isFORCING() ) s = onBit(s, FORCING_BIT);
              bx[m] = s;
            }
          }
        }
      }
      
    }
  }
}

/**
 @fn void VoxInfo::setControlVars(unsigned* r_size, unsigned r_guide)
 @brief 必要な変数をコピーする
 @param r_size 分割数
 @param r_guide ガイドセル
 */
void VoxInfo::setControlVars(unsigned* r_size, unsigned r_guide)
{
  guide   = r_guide;
  size[0] = r_size[0];
  size[1] = r_size[1];
  size[2] = r_size[2];
}

/**
 @fn void VoxInfo::setInactive_Compo(unsigned id, int def, int* mid, unsigned* bh1, unsigned* bh2)
 @brief コンポーネントに接するdef_faceのIDをもつセルを不活性セルにする
 @param id 対象BCのセルID
 @param def 面指定ID
 @param mid ボクセル配列
 @param bh1 BCindex H1
 @param bh2 BCindex H2
 */
void VoxInfo::setInactive_Compo(unsigned id, int def, int* mid, unsigned* bh1, unsigned* bh2)
{
  int i,j,k;
  unsigned m_p, m_e, m_w, m_n, m_s, m_t, m_b;
  int      c_p, c_e, c_w, c_n, c_s, c_t, c_b;
  int idd = (int)id;
  
  int ix = (int)size[0];
  int jx = (int)size[1];
  int kx = (int)size[2];
  int gd = (int)guide;
  
  for (k=1; k<=kx; k++) {
    for (j=1; j<=jx; j++) {
      for (i=1; i<=ix; i++) {
        m_p = FBUtility::getFindexS3D(size, guide, i  , j  , k  );
        c_p = mid[m_p];
        
        if ( c_p == idd ) { // BCの対象セル
          m_e = FBUtility::getFindexS3D(size, guide, i+1, j  , k  );
          m_w = FBUtility::getFindexS3D(size, guide, i-1, j  , k  );
          m_n = FBUtility::getFindexS3D(size, guide, i  , j+1, k  );
          m_s = FBUtility::getFindexS3D(size, guide, i  , j-1, k  );
          m_t = FBUtility::getFindexS3D(size, guide, i  , j  , k+1);
          m_b = FBUtility::getFindexS3D(size, guide, i  , j  , k-1);
          
          c_e = mid[m_e];
          c_w = mid[m_w];
          c_n = mid[m_n];
          c_s = mid[m_s];
          c_t = mid[m_t];
          c_b = mid[m_b];
          
          // X-
          if ( c_w == def ) {
            bh1[m_w] = offBit( bh1[m_w], ACTIVE_BIT ); // defセルを不活性化
            bh2[m_w] = offBit( bh2[m_w], ACTIVE_BIT );
          }
          
          // X+
          if ( c_e == def ) {
            bh1[m_e] = offBit( bh1[m_e], ACTIVE_BIT );
            bh2[m_e] = offBit( bh2[m_e], ACTIVE_BIT );
          }
          
          // Y-
          if ( c_s == def ) {
            bh1[m_s] = offBit( bh1[m_s], ACTIVE_BIT );
            bh2[m_s] = offBit( bh2[m_s], ACTIVE_BIT );
          }
          
          // Y+
          if ( c_n == def ) {
            bh1[m_n] = offBit( bh1[m_n], ACTIVE_BIT );
            bh2[m_n] = offBit( bh2[m_n], ACTIVE_BIT );
          }
          
          // Z-
          if ( c_b == def ) {
            bh1[m_b] = offBit( bh1[m_b], ACTIVE_BIT );
            bh2[m_b] = offBit( bh2[m_b], ACTIVE_BIT );
          }
          
          // Z+
          if ( c_t == def ) {
            bh1[m_t] = offBit( bh1[m_t], ACTIVE_BIT );
            bh2[m_t] = offBit( bh2[m_t], ACTIVE_BIT );
          }
        }
        
      }
    }
  }
}

/**
 @fn void VoxInfo::setOBC_Cut(SetBC* BC, float* cut)
 @brief 外部境界のガイドセルが固体の場合に距離情報をセット
 @param BC SetBCクラスのポインタ
 @param[in/out] cut 距離情報
 */
void VoxInfo::setOBC_Cut(SetBC* BC, float* cut)
{
  BoundaryOuter* m_obc=NULL;
  unsigned F;
  
  const float pos=0.5f;
  unsigned m;
  
  const int ix = (int)size[0];
  const int jx = (int)size[1];
  const int kx = (int)size[2];
  
  for (int face=0; face<NOFACE; face++) {
    m_obc = BC->export_OBC(face);
    F = m_obc->get_Class();
    
    if ( (F == OBC_WALL) || (F == OBC_SYMMETRIC) ) {
      
      switch (face) {
        case X_MINUS:
          if( pn.nID[X_MINUS] < 0 ){ // 外部境界をもつノードのみ
            for (int k=1; k<=kx; k++) {
              for (int j=1; j<=jx; j++) {
                m = FBUtility::getFindexS3Dcut(size, guide, X_MINUS, 1, j, k);
                cut[m] = pos; 
              }
            }        
          }
          break;
          
        case X_PLUS:
          if( pn.nID[X_PLUS] < 0 ){
            for (int k=1; k<=kx; k++) {
              for (int j=1; j<=jx; j++) {
                m = FBUtility::getFindexS3Dcut(size, guide, X_PLUS, ix, j, k);
                cut[m] = pos;
              }
            }
          }
          break;
          
        case Y_MINUS:
          if( pn.nID[Y_MINUS] < 0 ){
            for (int k=1; k<=kx; k++) {
              for (int i=1; i<=ix; i++) {
                m = FBUtility::getFindexS3Dcut(size, guide, Y_MINUS, i, 1, k);
                cut[m] = pos; 
              }
            }
          }
          break;
          
        case Y_PLUS:
          if( pn.nID[Y_PLUS] < 0 ){
            for (int k=1; k<=kx; k++) {
              for (int i=1; i<=ix; i++) {
                m = FBUtility::getFindexS3Dcut(size, guide, Y_PLUS, i, jx, k);
                cut[m] = pos;
              }
            }
          }
          break;
          
        case Z_MINUS:
          if( pn.nID[Z_MINUS] < 0 ){
            for (int j=1; j<=jx; j++) {
              for (int i=1; i<=ix; i++) {
                m = FBUtility::getFindexS3Dcut(size, guide, Z_MINUS, i, j, 1);
                cut[m] = pos;
              }
            }
          }
          break;
          
        case Z_PLUS:
          if( pn.nID[Z_PLUS] < 0 ){
            for (int j=1; j<=jx; j++) {
              for (int i=1; i<=ix; i++) {
                m = FBUtility::getFindexS3Dcut(size, guide, Z_PLUS, i, j, kx);
                cut[m] = pos;
              }
            }
          }
          break;
      }
      
    }
  }
  
}



// 作業用のポインタコピー
void VoxInfo::setWorkList(CompoList* m_CMP, MediumList* m_MAT)
{
  if ( !m_CMP ) Exit(0);
  cmp = m_CMP;
  
  if ( !m_MAT ) Exit(0);
  mat = m_MAT;
}


/**
 @fn unsigned VoxInfo::test_opposite_cut(nt* bid, int* mid, const int solid_id)
 @brief Hole fillingのため対角セルのカットを調べ、挟まれていたらcut=0とする
 @param bid カット点のID
 @param mid
 @param solid_id 固体ID
 */
unsigned VoxInfo::test_opposite_cut(int* bid, int* mid, const int solid_id)
{
  unsigned m_111, m_112, m_113, m_121, m_122, m_123, m_131, m_132, m_133;
  unsigned m_211, m_212, m_213, m_221, m_222, m_223, m_231, m_232, m_233;
  unsigned m_311, m_312, m_313, m_321, m_322, m_323, m_331, m_332, m_333;
  int b_111, b_112, b_113, b_121, b_122, b_123, b_131, b_132, b_133;
  int b_211, b_212, b_213, b_221, b_222, b_223, b_231, b_232, b_233;
  int b_311, b_312, b_313, b_321, b_322, b_323, b_331, b_332, b_333;
  unsigned m_sz[3];
  m_sz[0] = size[0];
  m_sz[1] = size[1];
  m_sz[2] = size[2];
  unsigned gd = guide;
  
  int ix = (int)size[0];
  int jx = (int)size[1];
  int kx = (int)size[2];
  int c = 0; /// count
  int sd = solid_id;
  
  if ( sd >32 ) {
    Hostonly_ printf("Error : Solid ID must be under 2^5=32\n");
    Exit(0);
  }
  
  int flag = solid_id;       // shift 0
  flag |= ( solid_id << 5 ); // shift 5
  flag |= ( solid_id << 10);
  flag |= ( solid_id << 15);
  flag |= ( solid_id << 20);
  flag |= ( solid_id << 25);
  
#pragma omp parallel for firstprivate(ix, jx, kx, m_sz, gd, sd, flag) \
 private(m_111, m_112, m_113, m_121, m_122, m_123, m_131, m_132, m_133) \
 private(m_211, m_212, m_213, m_221, m_222, m_223, m_231, m_232, m_233) \
 private(m_311, m_312, m_313, m_321, m_322, m_323, m_331, m_332, m_333) \
 private(b_111, b_112, b_113, b_121, b_122, b_123, b_131, b_132, b_133) \
 private(b_211, b_212, b_213, b_221, b_222, b_223, b_231, b_232, b_233) \
 private(b_311, b_312, b_313, b_321, b_322, b_323, b_331, b_332, b_333) \
 schedule(static) reduction(+:c)
  
  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        
        m_111 = FBUtility::getFindexS3D(m_sz, gd, i-1, j-1, k-1);
        m_112 = FBUtility::getFindexS3D(m_sz, gd, i-1, j-1, k  );
        m_113 = FBUtility::getFindexS3D(m_sz, gd, i-1, j-1, k+1);
        m_121 = FBUtility::getFindexS3D(m_sz, gd, i-1, j  , k-1);
        m_122 = FBUtility::getFindexS3D(m_sz, gd, i-1, j  , k  );
        m_123 = FBUtility::getFindexS3D(m_sz, gd, i-1, j  , k+1);
        m_131 = FBUtility::getFindexS3D(m_sz, gd, i-1, j+1, k-1);
        m_132 = FBUtility::getFindexS3D(m_sz, gd, i-1, j+1, k  );
        m_133 = FBUtility::getFindexS3D(m_sz, gd, i-1, j+1, k+1);
        
        m_211 = FBUtility::getFindexS3D(m_sz, gd, i  , j-1, k-1);
        m_212 = FBUtility::getFindexS3D(m_sz, gd, i  , j-1, k  );
        m_213 = FBUtility::getFindexS3D(m_sz, gd, i  , j-1, k+1);
        m_221 = FBUtility::getFindexS3D(m_sz, gd, i  , j  , k-1);
        m_222 = FBUtility::getFindexS3D(m_sz, gd, i  , j  , k  );
        m_223 = FBUtility::getFindexS3D(m_sz, gd, i  , j  , k+1);
        m_231 = FBUtility::getFindexS3D(m_sz, gd, i  , j+1, k-1);
        m_232 = FBUtility::getFindexS3D(m_sz, gd, i  , j+1, k  );
        m_233 = FBUtility::getFindexS3D(m_sz, gd, i  , j+1, k+1);
            
        m_311 = FBUtility::getFindexS3D(m_sz, gd, i+1, j-1, k-1);
        m_312 = FBUtility::getFindexS3D(m_sz, gd, i+1, j-1, k  );
        m_313 = FBUtility::getFindexS3D(m_sz, gd, i+1, j-1, k+1);
        m_321 = FBUtility::getFindexS3D(m_sz, gd, i+1, j  , k-1);
        m_322 = FBUtility::getFindexS3D(m_sz, gd, i+1, j  , k  );
        m_323 = FBUtility::getFindexS3D(m_sz, gd, i+1, j  , k+1);
        m_331 = FBUtility::getFindexS3D(m_sz, gd, i+1, j+1, k-1);
        m_332 = FBUtility::getFindexS3D(m_sz, gd, i+1, j+1, k  );
        m_333 = FBUtility::getFindexS3D(m_sz, gd, i+1, j+1, k+1);
        
        b_111 = bid[m_111];
        b_112 = bid[m_112];
        b_113 = bid[m_113];
        b_121 = bid[m_121];
        b_122 = bid[m_122];
        b_123 = bid[m_123];
        b_131 = bid[m_131];
        b_132 = bid[m_132];
        b_133 = bid[m_133];
        
        b_211 = bid[m_211];
        b_212 = bid[m_212];
        b_213 = bid[m_213];
        b_221 = bid[m_221];
        b_222 = bid[m_222];
        b_223 = bid[m_223];
        b_231 = bid[m_231];
        b_232 = bid[m_232];
        b_233 = bid[m_233];
        
        b_311 = bid[m_311];
        b_312 = bid[m_312];
        b_313 = bid[m_313];
        b_321 = bid[m_321];
        b_322 = bid[m_322];
        b_323 = bid[m_323];
        b_331 = bid[m_331];
        b_332 = bid[m_332];
        b_333 = bid[m_333];
        
        
        if ( (b_122 != 0) && (b_322 != 0) ) {      // i
          bid[m_222] = flag;
          mid[m_222] = sd;
          c++;
        }
        else if ( (b_212 != 0) && (b_232 != 0) ) { // j
          bid[m_222] = flag;
          mid[m_222] = sd;
          c++;
        }
        else if ( (b_221 != 0) && (b_223 != 0) ) { // k
          bid[m_222] = flag;
          mid[m_222] = sd;
          c++;
        }
        else if ( (b_112 != 0) && (b_332 != 0) ) { // ij
          bid[m_222] = flag;
          mid[m_222] = sd;
          c++;
        }
        else if ( (b_132 != 0) && (b_312 != 0) ) { // ij
          bid[m_222] = flag;
          mid[m_222] = sd;
          c++;
        }
        else if ( (b_211 != 0) && (b_233 != 0) ) { // jk
          bid[m_222] = flag;
          mid[m_222] = sd;
          c++;
        }
        else if ( (b_213 != 0) && (b_231 != 0) ) { // jk
          bid[m_222] = flag;
          mid[m_222] = sd;
          c++;
        }
        else if ( (b_121 != 0) && (b_323 != 0) ) { // ki
          bid[m_222] = flag;
          mid[m_222] = sd;
          c++;
        }
        else if ( (b_123 != 0) && (b_321 != 0) ) { // ki
          bid[m_222] = flag;
          mid[m_222] = sd;
          c++;
        }
        else if ( (b_111 != 0) && (b_333 != 0) ) { // ijk
          bid[m_222] = flag;
          mid[m_222] = sd;
          c++;
        }
        else if ( (b_311 != 0) && (b_133 != 0) ) { // ijk
          bid[m_222] = flag;
          mid[m_222] = sd;
          c++;
        }
        else if ( (b_131 != 0) && (b_313 != 0) ) { // ijk
          bid[m_222] = flag;
          mid[m_222] = sd;
          c++;
        }
        else if ( (b_331 != 0) && (b_113 != 0) ) { // ijk
          bid[m_222] = flag;
          mid[m_222] = sd;
          c++;
        }
          


      }
    }
  }
  
  return c;
}
