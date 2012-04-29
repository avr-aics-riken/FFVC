/*
 * SPHERE - Skeleton for PHysical and Engineering REsearch
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2004-2012
 *
 */

//@file VoxInfo.C
//@brief FlowBase VoxInfo class
//@author keno, FSI Team, VCAD, RIKEN

#include <set>
#include <algorithm>
#include <map>
#include "VoxInfo.h"
extern SklParaComponent* ParaCmpo;

/**
 @fn void VoxInfo::adjCellID_on_GC(int face, SklScalar3D<int>* d_mid, int BCtype, int c_id, unsigned prdc_mode)
 @brief 外部境界に接するガイドセルのmid[]にIDをエンコードする
 @param face 外部境界面番号
 @param d_mid ID配列のデータクラス
 @param BCtype 外部境界面の境界条件の種類
 @param c_id セルID
 @param prdc_mode 周期境界条件のモード
 @note ガイドセル全てを対象
 */
void VoxInfo::adjCellID_on_GC(int face, SklScalar3D<int>* d_mid, int BCtype, int c_id, unsigned prdc_mode)
{
  SklParaManager* para_mng = ParaCmpo->GetParaManager();
  int i, j, k;
  register int ref_id;
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
          for (k=1; k<=kx; k++) {
            for (j=1; j<=jx; j++) {
              m0 = FBUtility::getFindexS3D(size, guide, 1, j, k); // 最外層のID
              ref_id = mid[m0];
              
              for (i=1-gd; i<=0; i++) {
                m = FBUtility::getFindexS3D(size, guide, i, j, k);
                mid[m] = ( find_ID_state(ref_id) == SOLID ) ? ref_id : c_id;
              }
            }
          }
        }
        break;
        
      case X_PLUS:
        if( pn.nID[face] < 0 ){
          for (k=1; k<=kx; k++) {
            for (j=1; j<=jx; j++) {
              m0 = FBUtility::getFindexS3D(size, guide, ix, j, k); // 最外層のID
              ref_id = mid[m0];
              
              for (i=ix+1; i<=ix+gd; i++) {
                m = FBUtility::getFindexS3D(size, guide, i, j, k);
                mid[m] = ( find_ID_state(ref_id) == SOLID ) ? ref_id : c_id;
              }
            }
          }
        }
        break;
        
      case Y_MINUS:
        if( pn.nID[face] < 0 ){
          for (k=1; k<=kx; k++) {
            for (i=1; i<=ix; i++) {
              m0 = FBUtility::getFindexS3D(size, guide, i, 1, k); // 最外層のID
              ref_id = mid[m0];
              
              for (j=1-gd; j<=0; j++) {
                m = FBUtility::getFindexS3D(size, guide, i, j, k);
                mid[m] = ( find_ID_state(ref_id) == SOLID ) ? ref_id : c_id;
              }
            }
          }
        }
        break;
        
      case Y_PLUS:
        if( pn.nID[face] < 0 ){
          for (k=1; k<=kx; k++) {
            for (i=1; i<=ix; i++) {
              m0 = FBUtility::getFindexS3D(size, guide, i, jx, k); // 最外層のID
              ref_id = mid[m0];
              
              for (j=jx+1; j<=jx+gd; j++) {
                m = FBUtility::getFindexS3D(size, guide, i, j, k);
                mid[m] = ( find_ID_state(ref_id) == SOLID ) ? ref_id : c_id;
              }
            }
          }
        }
        break;
        
      case Z_MINUS:
        if( pn.nID[face] < 0 ){
          for (j=1; j<=jx; j++) {
            for (i=1; i<=ix; i++) {
              m0 = FBUtility::getFindexS3D(size, guide, i, j, 1); // 最外層のID
              ref_id = mid[m0];
              
              for (k=1-gd; k<=0; k++) {
                m = FBUtility::getFindexS3D(size, guide, i, j, k);
                mid[m] = ( find_ID_state(ref_id) == SOLID ) ? ref_id : c_id;
              }
            }
          }
        }
        break;
        
      case Z_PLUS:
        if( pn.nID[face] < 0 ){
          for (j=1; j<=jx; j++) {
            for (i=1; i<=ix; i++) {
              m0 = FBUtility::getFindexS3D(size, guide, i, j, kx); // 最外層のID
              ref_id = mid[m0];
              
              for (k=kx+1; k<=kx+gd; k++) {
                m = FBUtility::getFindexS3D(size, guide, i, j, k);
                mid[m] = ( find_ID_state(ref_id) == SOLID ) ? ref_id : c_id;
              }
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
      if( para_mng->IsParallel() ){
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
              for (k=1; k<=kx; k++) {
                for (j=1; j<=jx; j++) {
                  for (i=1-gd; i<=0; i++) {
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
              for (k=1; k<=kx; k++) {
                for (j=1; j<=jx; j++) {
                  for (i=ix+1; i<=ix+gd; i++) {
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
              for (k=1; k<=kx; k++) {
                for (j=1-gd; j<=0; j++) {
                  for (i=1; i<=ix; i++) {
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
              for (k=1; k<=kx; k++) {
                for (j=jx+1; j<=jx+gd; j++) {
                  for (i=1; i<=ix; i++) {
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
              for (k=1-gd; k<=0; k++) {
                for (j=1; j<=jx; j++) {
                  for (i=1; i<=ix; i++) {
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
              for (k=kx+1; k<=kx+gd; k++) {
                for (j=1; j<=jx; j++) {
                  for (i=1; i<=ix; i++) {
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
 @fn void VoxInfo::adjCellID_Prdc_Inner(int face, SklScalar3D<int>* d_mid, int BCtype, unsigned c_id, unsigned prdc_mode)
 @brief 外部境界に接するガイドセルのmid[]にIDをエンコードする（内部周期境界の場合）
 @param face 外部境界面番号
 @param d_mid ID配列のデータクラス
 @param BCtype 外部境界面の境界条件の種類
 @param c_id セルID
 @param prdc_mode 周期境界条件のモード
 @note ガイドセル全てを対象
 */
void VoxInfo::adjCellID_Prdc_Inner(SklScalar3D<int>* d_mid)
{
  SklParaManager* para_mng = ParaCmpo->GetParaManager();
  int st[3], ed[3], dir, id;
  
  for (unsigned n=1; n<=NoBC; n++) {
    cmp[n].getBbox(st, ed);
    dir = (int)cmp[n].getPeriodicDir();
    id  = cmp[n].getID();
    if ( cmp[n].getType() == PERIODIC ) {
      if( para_mng->IsParallel() ){
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
  unsigned id;
  REAL_TYPE a, nvx, nvy, nvz;
  int cijk[3], dir[3], def, area=0;
  REAL_TYPE ai, aj, ak;
  
  cijk[0] = cijk[1] = cijk[2] = 0;
  dir[0] = dir[1] = dir[2] = 0;
  
  id  = cmp[n].getID();
  def = cmp[n].getDef();
  
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
      countVolumeEdge(n, bd, cijk);
      getNormalSign(n, gi, bd, dir);
      ai = 0.5*(REAL_TYPE)cijk[0]; // 両面あるので半分にする?
      aj = 0.5*(REAL_TYPE)cijk[1];
      ak = 0.5*(REAL_TYPE)cijk[2];
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
 @fn bool VoxInfo::check_fill(const int* mid)
 @brief ペイント済みかどうかをチェックする
 @param[in] mid ID配列
 */
bool VoxInfo::check_fill(const int* mid)
{
  for (int k=1; k<=(int)size[2]; k++) {
    for (int j=1; j<=(int)size[1]; j++) {
      for (int i=1; i<=(int)size[0]; i++) {
        
        if ( mid[FBUtility::getFindexS3D(size, guide, i  , j  , k  )] == 0 ) return false;
        
      }
    }
  }
  
  return true;
}


/**
 @fn bool VoxInfo::chkIDconsistency(IDtable* iTable, unsigned m_NoID)
 @brief XMLとスキャンしたボクセルIDの同一性をチェック
 @retval エラーコード
 @param iTable IDtableのリスト
 @param m_NoID XMLリストに記述されたIDの個数
 @note m_NoID >= NoVoxIDのはず
 */
bool VoxInfo::chkIDconsistency(IDtable* iTable, unsigned m_NoID)
{
  bool* chkflag = NULL;
	if( !(chkflag = new bool[NoVoxID+1]) ) return false;
	
	for (int i=0; i<=NoVoxID; i++) chkflag[i] = false;
	
	for (int i=1; i<=NoVoxID; i++) { // スキャンしたIDのループ
		for (int j=1; j<=m_NoID; j++) { // XMLのループ
			if ( colorList[i] == (int)iTable[j].getID() ) chkflag[i] = true;
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
void VoxInfo::countCellState(unsigned& Lcell, unsigned& Gcell, unsigned* bx, const unsigned state)
{
  
  SklParaManager* para_mng = ParaCmpo->GetParaManager();
  unsigned cell=0;    // local
  unsigned g_cell=0;  // global 
  int i,j,k;
  unsigned m;
  
  int ix = (int)size[0];
  int jx = (int)size[1];
  int kx = (int)size[2];
  int gd = (int)guide;
  
  // described in Fortran index
  for (k=1; k<=kx; k++) {
    for (j=1; j<=jx; j++) {
      for (i=1; i<=ix; i++) {
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
  Lcell = cell;
  
  if( para_mng->IsParallel() ) {
    unsigned c_tmp = g_cell;
    para_mng->Allreduce(&c_tmp, &g_cell, 1, SKL_ARRAY_DTYPE_UINT, SKL_SUM, pn.procGrp);
  }
  Gcell = g_cell;

}


/**
 @fn void VoxInfo::countQfaceEdge(unsigned n, unsigned* bx, int* cc)
 @brief HEATFLUX, TRANSFER, ISOTHERMAL面の断面積と法線を求める
 @param n 境界条件番号
 @param bx BCindex
 @param cc[out] カウントしたセル数
 @note
 - ボリュームのIDを見て，そのエッジを通過したときにカウントする
 - エッジの方向は，対象IDからFluidセルに向かう方向をプラス
 - 計算領域内部のノード間境界部でのカウント重複を防ぐため，開始点位置の判断を行う
 - 外部境界に接する面では，補正せず，常にoffset=1
 - 内部境界面では，重複しないようにoffset=0
 - normal is evaluated by only plus dir., because i have no good idea!
 */
void VoxInfo::countFace_S(unsigned n, unsigned* bx, int* cc)
{
	SklParaManager* para_mng = ParaCmpo->GetParaManager();
  unsigned register s, m;
  int i,j,k, c[3];
	int st[3], ed[3];
  
  c[0] = c[1] = c[2] = 0;
	
  cmp[n].getBbox(st, ed);
  
  for (k=st[2]; k<=ed[2]; k++) {
    for (j=st[1]; j<=ed[1]; j++) {
      for (i=st[0]; i<=ed[0]; i++) {
        m = FBUtility::getFindexS3D(size, guide, i, j, k);
        s = bx[m];
        
        if ( GET_FACE_BC(s, BC_FACE_W) == n ) c[0]++;
        if ( GET_FACE_BC(s, BC_FACE_E) == n ) c[0]--;
        if ( GET_FACE_BC(s, BC_FACE_S) == n ) c[1]++;
        if ( GET_FACE_BC(s, BC_FACE_N) == n ) c[1]--;
        if ( GET_FACE_BC(s, BC_FACE_B) == n ) c[2]++;
        if ( GET_FACE_BC(s, BC_FACE_T) == n ) c[2]--;
      }
    }
  }
  
	if( para_mng->IsParallel() ){
    int tmp[3];
		tmp[0] = c[0];
		tmp[1] = c[1];
		tmp[2] = c[2];
		para_mng->Allreduce(tmp, c, 3, SKL_ARRAY_DTYPE_INT, SKL_SUM, pn.procGrp);
	}
	cc[0] = c[0];
	cc[1] = c[1];
	cc[2] = c[2];
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
	SklParaManager* para_mng = ParaCmpo->GetParaManager();
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
  
	if( para_mng->IsParallel() ){
    int tmp[3];
		tmp[0] = c[0];
		tmp[1] = c[1];
		tmp[2] = c[2];
		para_mng->Allreduce(tmp, c, 3, SKL_ARRAY_DTYPE_INT, SKL_SUM, pn.procGrp);
    
    tmp[0] = ar;
    para_mng->Allreduce(tmp, &ar, 1, SKL_ARRAY_DTYPE_INT, SKL_SUM, pn.procGrp);
	}
	cc[0] = c[0];
	cc[1] = c[1];
	cc[2] = c[2];
}

/**
 @fn void VoxInfo::count_ValidCell_OBC(int face, unsigned* bv)
 @brief 外部境界面の有効セル数をカウントする
 @param face 外部境界面番号
 @param bv BCindex V
 @note 
 - 外部境界面の両側のセルがFluidのときのみカウント
 */
unsigned VoxInfo::count_ValidCell_OBC(int face, unsigned* bv)
{
  SklParaManager* para_mng = ParaCmpo->GetParaManager();
  int i, j, k;
  unsigned m1, m2, g=0;
  unsigned register s1, s2;
  
  int ix = (int)size[0];
  int jx = (int)size[1];
  int kx = (int)size[2];
  int gd = (int)guide;
  
  switch (face) {
    case X_MINUS:
      if( pn.nID[X_MINUS] < 0 ){ // 外部境界をもつノードのみ
        for (k=1; k<=kx; k++) {
          for (j=1; j<=jx; j++) {
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
        for (k=1; k<=kx; k++) {
          for (j=1; j<=jx; j++) {
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
        for (k=1; k<=kx; k++) {
          for (i=1; i<=ix; i++) {
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
        for (k=1; k<=kx; k++) {
          for (i=1; i<=ix; i++) {
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
        for (j=1; j<=jx; j++) {
          for (i=1; i<=ix; i++) {
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
        for (j=1; j<=jx; j++) {
          for (i=1; i<=ix; i++) {
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
  
  if( para_mng->IsParallel() ) {
    unsigned tmp = g;
    para_mng->Allreduce(&tmp, &g, 1, SKL_ARRAY_DTYPE_UINT, SKL_SUM, pn.procGrp);
  }
  return g;
}

/**
 @fn void VoxInfo::countVolumeEdge(unsigned n, unsigned* bx, int* cc)
 @brief モデルから法線を計算し，コンポーネントの断面積を求める
 @param n 境界条件番号
 @param bx BCindex ID
 @param cc カウントしたセル数，戻り値
 @note
 - ボリュームのIDを見て，そのエッジを通過したときにカウントする
 - ただし，F-Fセルのインターフェイスのみで，計算領域境界のエッジはカウントしない
 - 計算領域内部のノード間境界部でのカウント重複を防ぐため，開始点位置の判断を行う
 - 外部境界に接する面では，補正せず，常にoffset=1
 - 内部境界面では，重複しないようにoffset=0
 */
void VoxInfo::countVolumeEdge(unsigned n, unsigned* bx, int* cc)
{
	SklParaManager* para_mng = ParaCmpo->GetParaManager();
  unsigned m0, mi, mj, mk;
  unsigned s0, si, sj, sk;
  int i,j,k;
	int st[3], ed[3], ofst[3];
  int c[3], tmp[3];
  c[0] = c[1] = c[2] = 0;
	
  if ( cmp[n].isEns() ) {
    
    // コンポーネントのローカルインデクス
    cmp[n].getBbox(st, ed);
    getOffset(st, ofst);
    
    // Face between i and i+1
    for (k=st[2]; k<=ed[2]; k++) {
      for (j=st[1]; j<=ed[1]; j++) {
        for (i=st[0]-ofst[0]; i<=ed[0]; i++) {
          if ( (pn.nID[X_MINUS] < 0) && (i==0) ) continue;
          if ( (pn.nID[X_PLUS]  < 0) && (i==size[0]) ) continue;
          m0 = FBUtility::getFindexS3D(size, guide, i  , j  , k  );
          mi = FBUtility::getFindexS3D(size, guide, i+1, j  , k  );
          s0 = bx[m0];
          si = bx[mi];
          if ( ( ((s0 & MASK_6) == n) && ((si & MASK_6) != n) )
              ||  ( ((s0 & MASK_6) != n) && ((si & MASK_6) == n) ) ) {
            if ( IS_FLUID(si) && IS_FLUID(s0) ) c[0]++;
          }
        }
      }
    }
    
    // Face between j and j+1
    for (k=st[2]; k<=ed[2]; k++) {
      for (j=st[1]-ofst[1]; j<=ed[1]; j++) {
        for (i=st[0]; i<=ed[0]; i++) {
          if ( (pn.nID[Y_MINUS] < 0) && (j==0) ) continue;
          if ( (pn.nID[Y_PLUS]  < 0) && (j==size[1]) ) continue;
          m0 = FBUtility::getFindexS3D(size, guide, i  , j  , k  );
          mj = FBUtility::getFindexS3D(size, guide, i  , j+1, k  );
          s0 = bx[m0];
          sj = bx[mj];
          if ( ( ((s0 & MASK_6) == n) && ((sj & MASK_6) != n) )
              ||  ( ((s0 & MASK_6) != n) && ((sj & MASK_6) == n) ) ) {
            if ( IS_FLUID(sj) && IS_FLUID(s0) ) c[1]++;
          }
        }
      }
    }
    
    // Face between k and k+1
    for (k=st[2]-ofst[2]; k<=ed[2]; k++) {
      for (j=st[1]; j<=ed[1]; j++) {
        for (i=st[0]; i<=ed[0]; i++) {
          if ( (pn.nID[Z_MINUS] < 0) && (k==0) ) continue;
          if ( (pn.nID[Z_PLUS]  < 0) && (k==size[2]) ) continue;
          m0 = FBUtility::getFindexS3D(size, guide, i  , j  , k  );
          mk = FBUtility::getFindexS3D(size, guide, i  , j  , k+1);
          s0 = bx[m0];
          sk = bx[mk];
          if ( ( ((s0 & MASK_6) == n) && ((sk & MASK_6) != n) )
              ||  ( ((s0 & MASK_6) != n) && ((sk & MASK_6) == n) ) ) {
            if ( IS_FLUID(sk) && IS_FLUID(s0) ) c[2]++;
          }
        }
      }
    }
    
  } // Ens
  
	if( para_mng->IsParallel() ){
		tmp[0] = c[0];
		tmp[1] = c[1];
		tmp[2] = c[2];
		para_mng->Allreduce(tmp, c, 3, SKL_ARRAY_DTYPE_INT, SKL_SUM, pn.procGrp);
	}
	cc[0] = c[0];
	cc[1] = c[1];
	cc[2] = c[2];
}


/**
 @fn void VoxInfo::countOpenAreaOfDomain(unsigned* bx, REAL_TYPE* OpenArea)
 @brief  計算領域の外部境界で外側1層と内側の両方が流体セル数の場合にカウントする
 @param bx BCindex ID
 @param OpenArea 開口セル数
 */
void VoxInfo::countOpenAreaOfDomain(unsigned* bx, REAL_TYPE* OpenArea)
{
  SklParaManager* para_mng = ParaCmpo->GetParaManager();
  int i,j,k;
  unsigned m0, m1, g;
  unsigned m_area[NOFACE], tmp[NOFACE];
  
  int ix = (int)size[0];
  int jx = (int)size[1];
  int kx = (int)size[2];
  int gd = (int)guide;
  
  for (i=0; i<NOFACE; i++) {
    OpenArea[i]=0.0;
    m_area[i] = 0.0;
  }
  
  // described in Fortran index
  
  // X_MINUS
  g=0;
  if( pn.nID[X_MINUS] < 0 ){
    for (k=1; k<=kx; k++) {
      for (j=1; j<=jx; j++) {
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
    for (k=1; k<=kx; k++) {
      for (j=1; j<=jx; j++) {
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
    for (k=1; k<=kx; k++) {
      for (i=1; i<=ix; i++) {
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
    for (k=1; k<=kx; k++) {
      for (i=1; i<=ix; i++) {
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
    for (j=1; j<=jx; j++) {
      for (i=1; i<=ix; i++) {
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
    for (j=1; j<=jx; j++) {
      for (i=1; i<=ix; i++) {
        m0 = FBUtility::getFindexS3D(size, guide, i, j, kx  );
        m1 = FBUtility::getFindexS3D(size, guide, i, j, kx+1);
        if ( (  FLUID == IS_FLUID(bx[m0]))
            && (FLUID == IS_FLUID(bx[m1])) ) g++;
      }
    }
    m_area[Z_PLUS] = g;
  }
  
  if( para_mng->IsParallel() ){ 
    for (i=0; i<NOFACE; i++) tmp[i] = m_area[i];
    para_mng->Allreduce(tmp, m_area, NOFACE, SKL_ARRAY_DTYPE_UINT, SKL_SUM, pn.procGrp);
  }
  
  for (i=0; i<NOFACE; i++) OpenArea[i] = (REAL_TYPE)m_area[i];
}


/**
 @fn unsigned VoxInfo::countState(unsigned id, int* mid)
 @brief 媒質idの数を数え，値を返す
 @retval 計算空間内における媒質idの数
 @param id カウントするid
 @param mid ボクセルID配列
 */
unsigned VoxInfo::countState(unsigned id, int* mid)
{
  SklParaManager* para_mng = ParaCmpo->GetParaManager();
  int i,j,k;
  unsigned m;
  unsigned g=0, tmp=0;
  
  int ix = (int)size[0];
  int jx = (int)size[1];
  int kx = (int)size[2];
  int gd = (int)guide;
  
  // サーチ範囲はノードローカルの計算セル内
  for (k=1; k<=kx; k++) {
    for (j=1; j<=jx; j++) {
      for (i=1; i<=ix; i++) {
        m = FBUtility::getFindexS3D(size, guide, i, j, k);
        if ( mid[m] == (int)id )  g++;
      }
    }
  }
  
  if( para_mng->IsParallel() ) {
    tmp = g;
    para_mng->Allreduce(&tmp, &g, 1, SKL_ARRAY_DTYPE_UINT, SKL_SUM, pn.procGrp);
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
  unsigned m, id, q, s, d;
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
        id = cmp[q].getID();
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
void VoxInfo::encActive(unsigned& Lcell, unsigned& Gcell, unsigned* bx, unsigned KOS)
{
  SklParaManager* para_mng = ParaCmpo->GetParaManager();
  int i,j,k;
  unsigned m, c=0, g=0;
  unsigned register s;
  
  int ix = (int)size[0];
  int jx = (int)size[1];
  int kx = (int)size[2];
  int gd = (int)guide;
  
  switch ( KOS ) {
    case FLOW_ONLY:
    case THERMAL_FLOW:
    case THERMAL_FLOW_NATURAL:
      for (k=1; k<=kx; k++) {
        for (j=1; j<=jx; j++) {
          for (i=1; i<=ix; i++) {
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
      for (k=1; k<=kx; k++) {
        for (j=1; j<=jx; j++) {
          for (i=1; i<=ix; i++) {
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
      for (k=1; k<=kx; k++) {
        for (j=1; j<=jx; j++) {
          for (i=1; i<=ix; i++) {
            m = FBUtility::getFindexS3D(size, guide, i, j, k);
            bx[m] = onBit( bx[m], ACTIVE_BIT );
            c++;
          }
        }
      }
      break;
  }
  
  Lcell = g = c;
  if( para_mng->IsParallel() ) {
    unsigned c_tmp = g;
    para_mng->Allreduce(&c_tmp, &g, 1, SKL_ARRAY_DTYPE_UINT, SKL_SUM, pn.procGrp);
  }
  Gcell = g;
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
unsigned VoxInfo::encodeOrder(unsigned order, unsigned id, int* mid, unsigned* bx)
{
  SklParaManager* para_mng = ParaCmpo->GetParaManager();
  int idd;
  unsigned register m, g=0;
  
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
  
  if( para_mng->IsParallel() ) {
    unsigned tmp = g;
    para_mng->Allreduce(&tmp, &g, 1, SKL_ARRAY_DTYPE_UINT, SKL_SUM, pn.procGrp);
  }

  return g;
}


/**
 @fn unsigned VoxInfo::encQfaceHT_S(unsigned order, unsigned id, int* mid, unsigned* bcd, unsigned* bh1, unsigned* bh2, int deface)
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
unsigned VoxInfo::encQfaceHT_S(unsigned order, unsigned id, int* mid, unsigned* bcd, unsigned* bh1, unsigned* bh2, int deface)
{
  SklParaManager* para_mng = ParaCmpo->GetParaManager();
  int i,j,k,idd;
  unsigned g=0;
  unsigned m_p, m_e, m_w, m_n, m_s, m_t, m_b;
  int      c_p, c_e, c_w, c_n, c_s, c_t, c_b;
  unsigned s_e, s_w, s_n, s_s, s_t, s_b;
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
  
  if( para_mng->IsParallel() ) {
    unsigned tmp=g;
    para_mng->Allreduce(&tmp, &g, 1, SKL_ARRAY_DTYPE_UINT, SKL_SUM, pn.procGrp);
  }
  return g;
}

/**
 @fn unsigned VoxInfo::encQfaceHT_B(unsigned order, unsigned id, int* mid, unsigned* bcd, unsigned* bh1, unsigned* bh2, int deface)
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
unsigned VoxInfo::encQfaceHT_B(unsigned order, unsigned id, int* mid, unsigned* bcd, unsigned* bh1, unsigned* bh2, int deface)
{
  SklParaManager* para_mng = ParaCmpo->GetParaManager();
  int i,j,k,idd;
  unsigned g=0;
  unsigned m_p, m_e, m_w, m_n, m_s, m_t, m_b;
  int      c_p, c_e, c_w, c_n, c_s, c_t, c_b;
  unsigned s_e, s_w, s_n, s_s, s_t, s_b;
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
  
  if( para_mng->IsParallel() ) {
    unsigned tmp=g;
    para_mng->Allreduce(&tmp, &g, 1, SKL_ARRAY_DTYPE_UINT, SKL_SUM, pn.procGrp);
  }
  return g;
}

/**
 @fn unsigned VoxInfo::encQfaceISO_SF(unsigned order, unsigned id, int* mid, unsigned* bcd, unsigned* bh1, unsigned* bh2, int deface)
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
unsigned VoxInfo::encQfaceISO_SF(unsigned order, unsigned id, int* mid, unsigned* bcd, unsigned* bh1, unsigned* bh2, int deface)
{
  SklParaManager* para_mng = ParaCmpo->GetParaManager();
  int i,j,k,idd;
  unsigned g=0;
  unsigned m_p, m_e, m_w, m_n, m_s, m_t, m_b;
  int      c_p, c_e, c_w, c_n, c_s, c_t, c_b;
  unsigned s_e, s_w, s_n, s_s, s_t, s_b;
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
  
  if( para_mng->IsParallel() ) {
    unsigned tmp=g;
    para_mng->Allreduce(&tmp, &g, 1, SKL_ARRAY_DTYPE_UINT, SKL_SUM, pn.procGrp);
  }
  return g;
}

/**
 @fn unsigned VoxInfo::encQfaceISO_SS(unsigned order, unsigned id, int* mid, unsigned* bcd, unsigned* bh1, unsigned* bh2, int deface)
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
unsigned VoxInfo::encQfaceISO_SS(unsigned order, unsigned id, int* mid, unsigned* bcd, unsigned* bh1, unsigned* bh2, int deface)
{
  SklParaManager* para_mng = ParaCmpo->GetParaManager();
  int i,j,k,idd;
  unsigned g=0;
  unsigned m_p, m_e, m_w, m_n, m_s, m_t, m_b;
  int      c_p, c_e, c_w, c_n, c_s, c_t, c_b;
  unsigned s_e, s_w, s_n, s_s, s_t, s_b;
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
  
  if( para_mng->IsParallel() ) {
    unsigned tmp=g;
    para_mng->Allreduce(&tmp, &g, 1, SKL_ARRAY_DTYPE_UINT, SKL_SUM, pn.procGrp);
  }
  return g;
}

/**
 @fn unsigned VoxInfo::encQface(unsigned order, unsigned id, int* mid, unsigned* bcd, unsigned* bh1, unsigned* bh2, int deface, bool flag)
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
unsigned VoxInfo::encQface(unsigned order, unsigned id, int* mid, unsigned* bcd, unsigned* bh1, unsigned* bh2, int deface, bool flag)
{
  SklParaManager* para_mng = ParaCmpo->GetParaManager();
  int i,j,k,idd;
  unsigned g=0;
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
  
  if( para_mng->IsParallel() ) {
    unsigned tmp=g;
    para_mng->Allreduce(&tmp, &g, 1, SKL_ARRAY_DTYPE_UINT, SKL_SUM, pn.procGrp);
  }
  return (g);
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
unsigned VoxInfo::encPbit_D_IBC(unsigned order, unsigned id, int* mid, unsigned* bcd, unsigned* bcp, int deface)
{
  SklParaManager* para_mng = ParaCmpo->GetParaManager();
  int i,j,k, idd;
  unsigned g=0, tmp=0, m;
  unsigned m_e, m_w, m_n, m_s, m_t, m_b;
  unsigned s_e, s_w, s_n, s_s, s_t, s_b;
  int      c_p, c_e, c_w, c_n, c_s, c_t, c_b;
  unsigned register s, d;
  
  int ix = (int)size[0];
  int jx = (int)size[1];
  int kx = (int)size[2];
  int gd = (int)guide;
  
  idd = (int)id;
  
  for (k=1; k<=kx; k++) {
    for (j=1; j<=jx; j++) {
      for (i=1; i<=ix; i++) {
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
  
  if( para_mng->IsParallel() ) {
    tmp = g;
    para_mng->Allreduce(&tmp, &g, 1, SKL_ARRAY_DTYPE_UINT, SKL_SUM, pn.procGrp);
  }
  return g;
}


/**
 @fn unsigned VoxInfo::encPbit_N_Binary(unsigned* bx)
 @brief bcp[]に壁面境界の圧力ノイマン条件のビットフラグと固体に隣接するFセルに方向フラグ，収束判定の有効フラグをエンコードする
 @param[in/out] bx BCindex P
 @retval 固体表面セル数
 @note 
 - 流体セルのうち，固体セルに隣接する面のノイマンフラグをゼロにする．ただし，内部領域のみ．
 - 固体セルに隣接する流体セルに方向フラグを付与する．全内部領域．
 */
unsigned VoxInfo::encPbit_N_Binary(unsigned* bx)
{
  SklParaManager* para_mng = ParaCmpo->GetParaManager();
  int i,j,k;
  unsigned m_p, m_e, m_w, m_n, m_s, m_t, m_b;
  unsigned register s;
  
  int ix = (int)size[0];
  int jx = (int)size[1];
  int kx = (int)size[2];
  int gd = (int)guide;
  
  // ノイマンフラグ
  for (k=1; k<=kx; k++) {
    for (j=1; j<=jx; j++) {
      for (i=1; i<=ix; i++) {
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
  unsigned c = 0;
  
  for (k=1; k<=kx; k++) {
    for (j=1; j<=jx; j++) {
      for (i=1; i<=ix; i++) {
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
  
  if( para_mng->IsParallel() ) {
    unsigned tmp = c;
    para_mng->Allreduce(&tmp, &c, 1, SKL_ARRAY_DTYPE_UINT, SKL_SUM, pn.procGrp);
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
unsigned VoxInfo::encPbit_N_Cut(unsigned* bx, float* cut, const bool convergence)
{
  SklParaManager* para_mng = ParaCmpo->GetParaManager();
  int i,j,k;
  size_t m_p, m;
  unsigned register s;
  float cp_e, cp_w, cp_n, cp_s, cp_t, cp_b;
  float* ct;
  
  int ix = (int)size[0];
  int jx = (int)size[1];
  int kx = (int)size[2];
  
  // ノイマンフラグ
  for (k=1; k<=kx; k++) {
    for (j=1; j<=jx; j++) {
      for (i=1; i<=ix; i++) {
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
  unsigned c = 0;
  float* pos; 
  float q;
  
  for (k=1; k<=kx; k++) {
    for (j=1; j<=jx; j++) {
      for (i=1; i<=ix; i++) {
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
  
  if( para_mng->IsParallel() ) {
    unsigned tmp = c;
    para_mng->Allreduce(&tmp, &c, 1, SKL_ARRAY_DTYPE_UINT, SKL_SUM, pn.procGrp);
  }
  
  // 収束判定の有効フラグ
  float q0, q1, q2, q3, q4, q5;
  unsigned g=0;
  for (k=1; k<=kx; k++) {
    for (j=1; j<=jx; j++) {
      for (i=1; i<=ix; i++) {
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
  
  
  if( para_mng->IsParallel() ) {
    unsigned tmp = g;
    para_mng->Allreduce(&tmp, &g, 1, SKL_ARRAY_DTYPE_UINT, SKL_SUM, pn.procGrp);
  }
  
  Hostonly_ printf("\tThe number of cells which are changed to INACTIVE and SOLID because of all faces are cut = %d\n\n", g);
  
  
  // カットのあるセルの収束判定をしないオプション
  if ( convergence ) {
    
    for (k=1; k<=kx; k++) {
      for (j=1; j<=jx; j++) {
        for (i=1; i<=ix; i++) {
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
            if ( (q0+q1+q2+q3+q4+q5) < 6.0 ) {
              s = offBit(s, VLD_CNVG);    // Out of scope
              g++;
            }
            
            bx[m_p] = s;
          }
        }
      }
    }
    
    if( para_mng->IsParallel() ) {
      unsigned tmp = g;
      para_mng->Allreduce(&tmp, &g, 1, SKL_ARRAY_DTYPE_UINT, SKL_SUM, pn.procGrp);
    }
    
    Hostonly_ printf("\tThe number of cells which are excluded to convergence judgement by cut = %d\n\n", g);
    
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
unsigned VoxInfo::encPbit_N_IBC(unsigned order, unsigned id, int* mid, unsigned* bcd, unsigned* bcp, int deface)
{
  SklParaManager* para_mng = ParaCmpo->GetParaManager();
  int i,j,k, idd;
  unsigned g=0, tmp=0, m;
  unsigned m_e, m_w, m_n, m_s, m_t, m_b;
  unsigned s_e, s_w, s_n, s_s, s_t, s_b;
  int      c_p, c_e, c_w, c_n, c_s, c_t, c_b;
  unsigned register s, d;
  
  int ix = (int)size[0];
  int jx = (int)size[1];
  int kx = (int)size[2];
  int gd = (int)guide;
  
  idd = (int)id;
  
  for (k=1; k<=kx; k++) {
    for (j=1; j<=jx; j++) {
      for (i=1; i<=ix; i++) {
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
  
  if( para_mng->IsParallel() ) {
    tmp = g;
    para_mng->Allreduce(&tmp, &g, 1, SKL_ARRAY_DTYPE_UINT, SKL_SUM, pn.procGrp);
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
unsigned VoxInfo::encVbit_IBC(unsigned order, unsigned id, int* mid, unsigned* bv, int deface, unsigned* bp)
{
  SklParaManager* para_mng = ParaCmpo->GetParaManager();
  int i,j,k, idd;
  unsigned g=0;
  unsigned m_p, m_e, m_w, m_n, m_s, m_t, m_b;
  int      c_p, c_e, c_w, c_n, c_s, c_t, c_b;
  unsigned register s, q;
  
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
  for (k=1; k<=kx; k++) {
    for (j=1; j<=jx; j++) {
      for (i=1; i<=ix; i++) {
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
  
  // 対象面数の集約
  if( para_mng->IsParallel() ) {
    unsigned tmp = g;
    para_mng->Allreduce(&tmp, &g, 1, SKL_ARRAY_DTYPE_UINT, SKL_SUM, pn.procGrp);
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
unsigned VoxInfo::encVbit_IBC_Cut(const unsigned order, 
                                  const unsigned id, 
                                  unsigned* bv, 
                                  unsigned* bp, 
                                  const float* cut, 
                                  const int* cut_id, 
                                  const float* vec, 
                                  const unsigned bc_dir)
{
  SklParaManager* para_mng = ParaCmpo->GetParaManager();
  int    idd;
  unsigned g=0;
  
  unsigned register s, q;
  const float *pos;
  int    bid;
  size_t m_p, m_c;
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
  
  idd = (int)id;
  
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
            //id_w = (bid >> 0)  & MASK_5;
            //id_e = (bid >> 5)  & MASK_5;
            //id_s = (bid >> 10) & MASK_5;
            //id_n = (bid >> 15) & MASK_5;
            //id_b = (bid >> 20) & MASK_5;
            //id_t = (bid >> 25) & MASK_5;
            
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
  
  // 対象面数の集約
  if( para_mng->IsParallel() ) {
    unsigned tmp = g;
    para_mng->Allreduce(&tmp, &g, 1, SKL_ARRAY_DTYPE_UINT, SKL_SUM, pn.procGrp);
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
 - adjCellID_on_GC()でガイドセル上のIDを指定済み．指定BCとの適合性をチェックする
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
 @fn unsigned VoxInfo::fill_cells(const int* bid, int* mid, const int tgt_id)
 @brief フィルを実行
 @param[in] bid カット点のID
 @param[in] mid ID配列
 @param[in] tgt_id サーチするID
 @note tgt_idと接する未ペイントセルをみつけ，シードとなるセルのインデクスを返す．シード点はtgt_id
 */
unsigned VoxInfo::fill_cells(const int* bid, int* mid, const int tgt_id)
{
  int target = tgt_id;
  unsigned m_sz[3];
  unsigned m_p, m_e, m_w, m_n, m_s, m_t, m_b;
  
  int ix = (int)size[0];
  int jx = (int)size[1];
  int kx = (int)size[2];
  int c = 0;
  
  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        m_p = FBUtility::getFindexS3D(size, guide, i  , j  , k  );
        m_e = FBUtility::getFindexS3D(size, guide, i+1, j  , k  );
        m_w = FBUtility::getFindexS3D(size, guide, i-1, j  , k  );
        m_n = FBUtility::getFindexS3D(size, guide, i  , j+1, k  );
        m_s = FBUtility::getFindexS3D(size, guide, i  , j-1, k  );
        m_t = FBUtility::getFindexS3D(size, guide, i  , j  , k+1);
        m_b = FBUtility::getFindexS3D(size, guide, i  , j  , k-1);

        // テストセルがtargetの場合
        if ( mid[m_p] == target ) {

          // 隣のIDがID=0 && カットがない場合にペイント
          if ( (mid[m_e] == 0) && (get_BID5(X_PLUS, bid[m_e]) == 0) ) {
            if (i != ix) {
              mid[m_e] = target;
              c++;
            }
          }
          else if ( (mid[m_n] == 0) && (get_BID5(Y_PLUS, bid[m_n]) == 0) ) {
            if (j != jx) {
              mid[m_n] = target;
              c++;
            }
          }
          else if ( (mid[m_t] == 0) && (get_BID5(Z_PLUS, bid[m_t]) == 0) ) {
            if (k != kx) {
              mid[m_t] = target;
              c++;
            }
          }
          else if ( (mid[m_b] == 0) && (get_BID5(Z_MINUS, bid[m_b]) == 0) ) {
            if (k != 1) {
              mid[m_b] = target;
              c++;
            }
          }
          else if ( (mid[m_s] == 0) && (get_BID5(Y_MINUS, bid[m_s]) == 0) ) {
            if (j != 1) {
              mid[m_s] = target;
              c++;
            }
          }
          else if ( (mid[m_w] == 0) && (get_BID5(X_MINUS, bid[m_w]) == 0) ) {
            if (i != 1) {
              mid[m_w] = target;
              c++;
            }
          }
          
        } // target        
      }
    }
  }
  
  return c;
}


/**
 @fn void VoxInfo::find_isolate_Fcell(unsigned order, int* mid, unsigned* bx)
 @brief 孤立した流体セルを探し，周囲の個体媒質で置換，BCindexを修正する
 @param order cmp[]に登録されたMaterialListへのエントリ番号
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
  int gd = (int)guide;
  
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
            s |= (order << TOP_MATERIAL); // sにMaterialListのエントリorderをエンコードする
            
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
    id  = mat[odr].getMatID();
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
  size_t m;
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
unsigned VoxInfo::flip_InActive(unsigned& L, unsigned& G, unsigned id, int* mid, unsigned* bx)
{
  SklParaManager* para_mng = ParaCmpo->GetParaManager();
  int i,j,k, c_p;
  unsigned s, m, c=0, g=0;
  int idd = (int)id;
  
  int ix = (int)size[0];
  int jx = (int)size[1];
  int kx = (int)size[2];
  int gd = (int)guide;
  
  for (k=1; k<=kx; k++) {
    for (j=1; j<=jx; j++) {
      for (i=1; i<=ix; i++) {
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
  
  L = g = c;
  if( para_mng->IsParallel() ) {
    unsigned c_tmp = g;
    para_mng->Allreduce(&c_tmp, &g, 1, SKL_ARRAY_DTYPE_UINT, SKL_SUM, pn.procGrp);
  }
  G = g;
  
  return g;
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
  
  unsigned id;
  REAL_TYPE a;
  int def;
  
  id  = cmp[n].getID();
  def = cmp[n].getDef();
  
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
 @fn void VoxInfo::getIDrange(const CfgElem *elmL2, const char* keyword, unsigned* var)
 @brief BC Indexからマスク値を計算する
 @param elmL2 XMLツリーのポインタ
 @param keyword XMLキーワード
 @param var レンジの値
 */
void VoxInfo::getIDrange(const CfgElem *elmL2, const char* keyword, unsigned* var)
{
  const CfgParam* param=NULL;
  int value=0;
  
  if ( !(param = elmL2->GetParamFirst(keyword)) ) {
    Hostonly_ stamped_printf("\tParsing error : Missing '%s' keyword\n", keyword);
    Exit(0);
  }
  else if ( !(param->GetData(&value)) ) {
    Hostonly_ stamped_printf("\tParsing error : Invalid integer value for '%s'\n", keyword);
    Exit(0);
  }
  else if (value<0) {
    Hostonly_ stamped_printf("\tParsing error : Negative integer value for '%s' %d\n", keyword, value);
    Exit(0);
  }
  else {
    *var = (unsigned)value;
  }
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

/**
 @fn void VoxInfo::getNormalSign(unsigned n, int* gi, unsigned* bx, int* dir)
 @brief コンポーネントの法線の符号を計算する
 @param n 境界条件番号
 @param gi コンポーネントのグローバルインデクス
 @param bx BCindex ID
 @param dir 法線ベクトルの符号
 */
void VoxInfo::getNormalSign(unsigned n, int* gi, unsigned* bx, int* dir)
{
  SklParaManager* para_mng = ParaCmpo->GetParaManager();
  unsigned m, q_ij[5], q_ik[5], tmp[5];
  int st[3], ed[3];
  unsigned ni, nj, nk;
  REAL_TYPE cn[3], tx, ty, tz;
  int m_sz[3], l_dir[3];
  int i,j,k;
  
  // コンポーネントのローカルインデクス
  cmp[n].getBbox(st, ed);
  
  // 対象コンポーネント領域のBbox領域の大きさ（セル数）
  ni = (unsigned)(gi[3] - gi[0]);
  nj = (unsigned)(gi[4] - gi[1]);
  nk = (unsigned)(gi[5] - gi[2]);
  
  // グローバルインデクスでの対象コンポーネント領域の中心座標
  cn[0] = 0.5*(REAL_TYPE)ni + (REAL_TYPE)gi[0];  
  cn[1] = 0.5*(REAL_TYPE)nj + (REAL_TYPE)gi[1];
  cn[2] = 0.5*(REAL_TYPE)nk + (REAL_TYPE)gi[2];
  
  m_sz[0] = para_mng->GetVoxelHeadIndex(pn.ID, 0) + 1;
  m_sz[1] = para_mng->GetVoxelHeadIndex(pn.ID, 1) + 1;
  m_sz[2] = para_mng->GetVoxelHeadIndex(pn.ID, 2) + 1;
  
  for (int l=0; l<5; l++) q_ij[l] = q_ik[l] = 0;
  
  for (k=st[2]; k<=ed[2]; k++) {
    for (j=st[1]; j<=ed[1]; j++) {
      for (i=st[0]; i<=ed[0]; i++) {
        m = FBUtility::getFindexS3D(size, guide, i  , j  , k  );
        if ( (bx[m] & MASK_6) == n ) {
          tx = (REAL_TYPE)(i+m_sz[0]) - cn[0]; // グローバルインデクスで評価
          ty = (REAL_TYPE)(j+m_sz[1]) - cn[1];
          tz = (REAL_TYPE)(k+m_sz[2]) - cn[2];
          
          getQuadrant(q_ij, tx, ty); // Quadrant IJ
          getQuadrant(q_ik, tx, tz); // Quadrant IK
        }
      }
    }
  }
	
	if( para_mng->IsParallel() ){
		for (i=0; i<5; i++) tmp[i] = q_ij[i];
		para_mng->Allreduce(tmp, q_ij, 5, SKL_ARRAY_DTYPE_UINT, SKL_SUM, pn.procGrp);
		
		for (i=0; i<5; i++) tmp[i] = q_ik[i];
		para_mng->Allreduce(tmp, q_ik, 5, SKL_ARRAY_DTYPE_UINT, SKL_SUM, pn.procGrp);
	}
  
  // 象限の判断
  map<unsigned, unsigned, std::greater<unsigned> > QuadrantIJ, QuadrantIK;
  map<unsigned, unsigned, std::greater<unsigned> >::iterator it;
  unsigned q[4][2], p;
  l_dir[0] = 1;
  l_dir[1] = 0;
  l_dir[2] = 0;
  
  // IJ キーにカウント数を代入，ゼロは未使用なので1から代入
  for (i=1; i<5; i++) QuadrantIJ.insert( pair<unsigned, unsigned>(q_ij[i], (unsigned)i) );  
  //for (it=QuadrantIJ.begin(); it!=QuadrantIJ.end(); it++) printf("%d : %d %d\n", pn.ID, it->first, it->second); printf("\n");
  p=0;
  for (it=QuadrantIJ.begin(); it!=QuadrantIJ.end(); it++) {
    q[p][0] = it->first;
    q[p][1] = it->second;
  }
  if ( q[0][0] == q[1][0] ) { // 最初と2番目が同じ値の場合，破棄
    l_dir[1] = ( (q[2][1] == 2) || (q[2][1] == 4) ) ? -1: 1; // 3番目の候補の象限が２，４の場合はマイナス
  }
  else { // 最初の値の判断
    l_dir[1] = ( (q[0][1] == 2) || (q[0][1] == 4) ) ? -1: 1;
  }
  
  // IK キーにカウント数を代入，ゼロは未使用なので1から代入
  for (i=1; i<5; i++) QuadrantIK.insert( pair<unsigned, unsigned>(q_ik[i], (unsigned)i) );  
  //for (it=QuadrantIK.begin(); it!=QuadrantIK.end(); it++) printf("%d : %d %d\n", pn.ID, it->first, it->second); printf("\n");
  p=0;
  for (it=QuadrantIK.begin(); it!=QuadrantIK.end(); it++) {
    q[p][0] = it->first;
    q[p][1] = it->second;
  }
  if ( q[0][0] == q[1][0] ) { // 最初と2番目が同じ値の場合，破棄
    l_dir[2] = ( (q[2][1] == 2) || (q[2][1] == 4) ) ? -1: 1; // 3番目の候補の象限が２，４の場合はマイナス
  }
  else { // 最初の値の判断
    l_dir[2] = ( (q[0][1] == 2) || (q[0][1] == 4) ) ? -1: 1;
  }
  
  //printf("%d : si=%d sj=%d sk=%d\n", pn.ID, dir[0], dir[1], dir[2]);
	dir[0] = l_dir[0];
	dir[1] = l_dir[1];
	dir[2] = l_dir[2];
  
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
  for (unsigned i=1; i<=NoVoxID; i++) {
    fprintf(fp,"\t\t%3d : ID = [%d]\n", i, colorList[i]);
    if ( colorList[i] == 0 ) {
      Hostonly_ stamped_fprintf(fp, "\t*** Voxel includes ID=0, which is invalid.  ***\n\n");
      Exit(0);
    }
  }
}


/**
 @fn VoxInfo::receiveCfgPtr(SklSolverConfig* cfg)
 @brief SklSolverConfigのポインタを受け取る
 @param cfg SklSolverConfigクラスのポインタ
 @retval エラーコード
 */
bool VoxInfo::receiveCfgPtr(SklSolverConfig* cfg)
{
  if ( !cfg ) return false;
  CF = cfg;
  return true;
}

/**
 @fn unsigned VoxInfo::scanCell(int *cell, unsigned count, unsigned* cid, unsigned ID_replace)
 @brief cellで保持されるボクセルid配列をスキャンし，coloList[]に登録する
 @retval 含まれるセルID数
 @param cell ボクセルIDを保持する配列
 @param count 外部境界のセルID数
 @param cid セルIDリスト 
 @param ID_replace ID[0]を置換するID
 @note
 - 重複がないようにcolorList[]に登録する
 - IDを昇順にソート
 */ 
unsigned VoxInfo::scanCell(int *cell, unsigned count, unsigned* cid, unsigned ID_replace)
{
  int target;
  int i,j,k;
  size_t m;
  
  int ix = (int)size[0];
  int jx = (int)size[1];
  int kx = (int)size[2];
  int gd = (int)guide;
  
  // ID[0]を置換するオプションが指定されている場合（ID_replaceに値がセット）
  if ( ID_replace != 0 ) {
    target = (int)ID_replace;
    Hostonly_ printf("\n\tID[0] is replaced by ID[%d]\n", target);
    for (k=1; k<=kx; k++) {
      for (j=1; j<=jx; j++) {
        for (i=1; i<=ix; i++) {
          m = FBUtility::getFindexS3D(size, guide, i, j, k);
          if ( cell[m] == 0 ) cell[m] = target;
        }
      }
    }
  }
  
  // 内部領域に対して，マイナスとゼロをチェック
  for (k=1; k<=kx; k++) {
    for (j=1; j<=jx; j++) {
      for (i=1; i<=ix; i++) {
        target = cell[FBUtility::getFindexS3D(size, guide, i, j, k)];
        if ( target<=0 ) {
          Hostonly_ stamped_printf("\tVoxel data includes non-positive ID [%d] at (%d, %d, %d)\n", target, i, j, k);
          Exit(0);
        }
      }
    }
  }
	
  // 内部領域の媒質IDを数え，colorSetに放り込む
  set<int> colorSet;
  colorSet.clear();
  for (k=1; k<=kx; k++) {
    for (j=1; j<=jx; j++) {
      for (i=1; i<=ix; i++) {
        m = FBUtility::getFindexS3D(size, guide, i, j, k);
        target = cell[m];
        set<int>::iterator itr = colorSet.find(target);
        if( itr == colorSet.end() ) colorSet.insert(target);
      }
    }
  }
  
  // 外部領域の媒質IDをcolorSetに追加する
  for (i=0; i<count; i++) {
    target = cid[i];
    set<int>::iterator itr = colorSet.find(target);
    if( itr == colorSet.end() ) colorSet.insert(target);
  }
  
  NoVoxID = colorSet.size(); // Localの数
	
  // 並列処理時の colorList[] の取得
  SklParaManager* para_mng = ParaCmpo->GetParaManager();
  if( para_mng->IsParallel() ){
    
    int myRank  = para_mng->GetMyID();
    int nodeNum = para_mng->GetNodeNum(pn.procGrp); // 全ノード数を取得
		
    // 各ノードのID数を集約するテーブルを作成
    int* tmpArray = allocTable(nodeNum);
    if( !tmpArray ) {
      colorSet.clear();
      Exit(0);
    }
    memset(tmpArray, 0, sizeof(int)*nodeNum);
    
    int* IDNumTable = allocTable(nodeNum);
    if( !IDNumTable ) {
      colorSet.clear();
      if( tmpArray ) { delete [] tmpArray; tmpArray=NULL; }
      Exit(0);
    }
    memset(IDNumTable, 0, sizeof(int)*nodeNum);
    
    // 各ノードのID数をIDNumTable[]に保持
    tmpArray[myRank] = NoVoxID; // Localの数
    if( !para_mng->Allreduce(tmpArray, IDNumTable, nodeNum, SKL_ARRAY_DTYPE_INT, SKL_MAX, pn.procGrp) ) {
      colorSet.clear();
      if( tmpArray )   { delete [] tmpArray;   tmpArray = NULL; }
      if( IDNumTable ) { delete [] IDNumTable; IDNumTable = NULL; }
      Exit(0);
    }
    delete [] tmpArray; tmpArray = NULL; // 一旦解放
		
    // 各ノードのcolorList[]を集約するリストを作成 >> tmpArray[]
    unsigned sizeOfIdList = 0;
    for(i=0; i<nodeNum; i++) sizeOfIdList += IDNumTable[i];
    if( !(tmpArray = allocTable(sizeOfIdList)) ) { 
      colorSet.clear();
      if( IDNumTable ) { delete [] IDNumTable; IDNumTable=NULL; }
      return 0;
    }
    memset(tmpArray, 0, sizeof(int)*sizeOfIdList);
    
    // 各ノードのID番号がストアされる先頭アドレスを計算
    unsigned mytopp = 0;
    for(i=0; i<myRank; i++) mytopp += IDNumTable[i];
    
    // tmpArray[]の該当アドレスにノードのID番号を書き込む >> tmpArray[]のノード担当外はゼロ
    set<int>::iterator itr = colorSet.begin();
    set<int>::iterator itrEnd = colorSet.end();
    i = mytopp;
    for( ; itr != itrEnd; itr++) { tmpArray[i++] = *itr; }
    
    // 集約用のリスト
    int* IdList = allocTable(sizeOfIdList);
    if( !IdList ) {
      colorSet.clear();
      if( tmpArray )   { delete [] tmpArray; tmpArray=NULL; }
      if( IDNumTable ) { delete [] IDNumTable; IDNumTable=NULL; }
      return 0;
    }
    memset(IdList, 0, sizeof(int)*sizeOfIdList);
    
    if( !para_mng->Allreduce(tmpArray, IdList, sizeOfIdList, SKL_ARRAY_DTYPE_INT, SKL_MAX, pn.procGrp) ) {
      colorSet.clear();
      if( tmpArray )   { delete [] tmpArray; tmpArray=NULL; }
      if( IDNumTable ) { delete [] IDNumTable; IDNumTable=NULL; }
      if( IdList )     { delete [] IdList; IdList=NULL; }
      return 0;
    }
    delete [] tmpArray; tmpArray = NULL;
		
    // create colorSet
    for(i=0; i<sizeOfIdList; i++) {
      set<int>::iterator itr = colorSet.find(IdList[i]);
      if( itr == colorSet.end() ) colorSet.insert(IdList[i]);
    }
    NoVoxID = colorSet.size(); // Global
		
    if( IDNumTable ) { delete [] IDNumTable; IDNumTable=NULL; }
    if( IdList )     { delete [] IdList; IdList=NULL; }
  }
	
  if( NoVoxID == 0 ) {
    if( colorList ) delete [] colorList;
    colorList = NULL;
    Exit(0);
  }
	
  // allocate colorList
  if( !(colorList = allocTable(NoVoxID+1)) ) {
    Hostonly_ stamped_printf("\tAllocation error : colorList[]\n");
    Exit(0);
  }
  for (unsigned i=0; i<=NoVoxID; i++) colorList[i] = 0;
	
  // colorList[]をcolorSetから作成する
  set<int>::iterator itr = colorSet.begin();
  set<int>::iterator itrEnd = colorSet.end();
  i = 1;
  for( ; itr != itrEnd; itr++) {
    colorList[i++] = *itr;
  }
	
  // clear colorSet
  colorSet.clear();
	
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
 @note 事前に，cmp[]へMaterialListへのエントリ番号をエンコードしておく -> cmp[].setMatOdr()
 */
void VoxInfo::setBCIndex_base1(unsigned* bx, int* mid, float* cvf)
{
  unsigned odr, id;
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
    id = cmp[n].getID();
    
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
                mid[m] = (int)id;
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

  // MaterialListのエントリをエンコードする
  for (unsigned n=1; n<=NoCompo; n++) {
    odr = cmp[n].getMatOdr();
    id  = cmp[n].getID();

    for (size_t m=0; m<nx; m++) {
      if ( mid[m] == (int)id ) bx[m] |= (odr << TOP_MATERIAL);
    }
  }

  // CompoListのエントリ　セル要素のみ
  for (unsigned n=1; n<=NoBC; n++) {
    id = cmp[n].getID();
    
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
 @note 事前に，cmp[]へMaterialListへのエントリをエンコードしておく -> cmp[].setMatOdr()
 */
void VoxInfo::setBCIndex_base2(unsigned* bx, int* mid, SetBC* BC, unsigned& Lcell, unsigned& Gcell, unsigned KOS)
{
  unsigned id, mat_id;
  
  // 孤立した流体セルの属性変更
  for (unsigned n=1; n<=NoCompo; n++) {
    find_isolated_Fcell( cmp[n].getMatOdr(), mid, bx);
  }
  
  // BCIndexにそのセルが計算に有効(active)かどうかをエンコードする．KindOfSolverによって異なる
  encActive(Lcell, Gcell, bx, KOS);
  
  // Inactive指定のセルを不活性にする
  unsigned m_L, m_G;
  for (unsigned n=1; n<=NoBC; n++) {
    id = cmp[n].getID();
    
    if ( cmp[n].getType() == INACTIVE ) {
      
      m_L = m_G = 0;
      cmp[n].setElement( flip_InActive(m_L, m_G, id, mid, bx) );
      Lcell -= m_L;
      Gcell -= m_G;
    }
  }
  
  // コンポーネントに登録された媒質のセル数を数え，elementにセットする
  for (unsigned n=NoBC+1; n<=NoCompo; n++) {
    id = cmp[n].getID();
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
  unsigned n, id;
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
    if ( BC->get_OBC_Ptr(face)->get_BCtype() == OBC_SYMMETRIC ) {
      encAmask_SymtrcBC(face, bh2);
    }
  }
  
  // 不活性セルの場合の断熱マスク処理
  for (unsigned n=1; n<=NoBC; n++) {
    if ( cmp[n].getType() == INACTIVE ) {
      setAmask_InActive(cmp[n].getID(), mid, bh2);
    }
  }
  
  // bh2の下位5ビットにはBCのエントリのみ(1~31)エンコード
  for (n=1; n<=NoBC; n++) {
    id = cmp[n].getID();
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
 @fn unsigned VoxInfo::setBCIndexP(unsigned* bcd, unsigned* bcp, int* mid, SetBC* BC, float* cut, bool isCDS)
 @brief 圧力境界条件のビット情報をエンコードする
 @param bcd BCindex ID
 @param bcp BCindex P
 @param mid ID配列
 @param BC SetBCクラスのポインタ
 @param cut 距離情報
 @param isCDS CDS->true
 @retval 表面セル数
 */
unsigned VoxInfo::setBCIndexP(unsigned* bcd, unsigned* bcp, int* mid, SetBC* BC, float* cut, bool isCDS)
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
    m_obc = BC->get_OBC_Ptr(face);
    F = m_obc->get_BCtype();
    
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
  unsigned n, id;
  int deface;
  
  for (n=1; n<=NoBC; n++) {
    id = cmp[n].getID();
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
  int i,j,k,m;
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
 @fn void VoxInfo::setBCIndexV(unsigned* bv, int* mid, SetBC* BC, unsigned* bp, float* cut, int* cut_id, bool isCDS)
 @brief bv[]に境界条件のビット情報をエンコードする
 @param bv BCindex V
 @param mid ID配列
 @param BC SetBCクラスのポインタ
 @param bp BCindex P
 @param cut 距離情報
 @param cut_id カット点ID
 @param isCDS CDS->true
 */
void VoxInfo::setBCIndexV(unsigned* bv, int* mid, SetBC* BC, unsigned* bp, float* cut, int* cut_id, bool isCDS)
{
  // ガイドセルの媒質情報をチェックし，流束形式のBCの場合にビットフラグをセット
  BoundaryOuter* m_obc=NULL;
  unsigned F;
  
  // 外部境界
  for (int face=0; face<NOFACE; face++) {
    m_obc = BC->get_OBC_Ptr(face);
    F = m_obc->get_BCtype();
    
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
  unsigned n, id, m_dir;
  int deface;
  float vec[3];
  
  for (n=1; n<=NoBC; n++) {
    id = cmp[n].getID();
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
  size_t m;
  
  const int ix = (int)size[0];
  const int jx = (int)size[1];
  const int kx = (int)size[2];
  
  for (int face=0; face<NOFACE; face++) {
    m_obc = BC->get_OBC_Ptr(face);
    F = m_obc->get_BCtype();
    
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

/**
 @fn void VoxInfo::setWorklist(CompoList* m_CMP, MaterialList* m_MAT)
 @brief 作業用のポインタコピー
 @param m_CMP
 @param m_MAT 
 */
void VoxInfo::setWorkList(CompoList* m_CMP, MaterialList* m_MAT)
{
  if ( !m_CMP ) Exit(0);
  cmp = m_CMP;
  
  if ( !m_MAT ) Exit(0);
  mat = m_MAT;
}

/**
 @fn unsigned VoxInfo::Solid_from_Cut(int* mid, float* cut, const int id)
 @brief ボクセルモデルにカット情報から得られた固体情報を転写する
 @param[in/out] mid セルID
 @param[in] cut 距離情報
 @param[in] id 固体ID 
 @retval 固体セル数
 */
unsigned VoxInfo::Solid_from_Cut(int* mid, float* cut, const int id)
{
  SklParaManager* para_mng = ParaCmpo->GetParaManager();
  unsigned c=0;
  int q, m_id;
  size_t m_p, m;
  
  unsigned m_sz[3], gd;
  m_sz[0] = size[0];
  m_sz[1] = size[1];
  m_sz[2] = size[2];
  gd = guide;
  m_id= id;
  
#pragma omp parallel for firstprivate(m_sz, gd, m_id) private(q) schedule(static) reduction(+:c)
  for (int k=1; k<=(int)m_sz[2]; k++) {
    for (int j=1; j<=(int)m_sz[1]; j++) {
      for (int i=1; i<=(int)m_sz[0]; i++) {
        
        m_p = FBUtility::getFindexS3D(m_sz, gd, i, j, k);
        m = FBUtility::getFindexS3Dcut(m_sz, gd, 0, i, j, k);
        q = 0;
        
        // セル内に交点があれば，壁
        if ( cut[m+0] <= 0.5f ) q++;
        if ( cut[m+1] <= 0.5f ) q++;
        if ( cut[m+2] <= 0.5f ) q++;
        if ( cut[m+3] <= 0.5f ) q++;
        if ( cut[m+4] <= 0.5f ) q++;
        if ( cut[m+5] <= 0.5f ) q++;
        
        if ( q > 0 ) {
          mid[m_p] = m_id;
          c++;
        }
      }
    }
  }
  
  if( para_mng->IsParallel() ) {
    unsigned tmp = c;
    para_mng->Allreduce(&tmp, &c, 1, SKL_ARRAY_DTYPE_UINT, SKL_SUM, pn.procGrp);
  }
  
  return c;
}
