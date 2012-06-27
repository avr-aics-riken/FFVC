#ifndef _FB_BINVOX_H_
#define _FB_BINVOX_H_

// #################################################################
//
// CAERU Library
//
// Copyright (c) All right reserved. 2012
//
// Institute of Industrial Science, The University of Tokyo, Japan. 
//
// #################################################################

/** 
 * @file VoxInfo.h
 * @brief FlowBase VoxInfo class Header
 * @author kero
 */

#include <math.h>
#include <cstdio>

#include "cpm_Define.h"
#include "cpm_ParaManager.h"

#include "FBUtility.h"
#include "Component.h"
#include "Medium.h"
#include "SetBC.h"
#include "BndOuter.h"
#include "vec3.h"

#include "Polylib.h"
#include "MPIPolylib.h"

#include "mpi.h"
#include "limits.h"

class VoxInfo {
  
private:
  int size[3];                   ///< 計算内部領域分割数(Local)
  int guide;                     ///< ガイドセルサイズ
  int NoBC;                      ///< 境界条件数
  int NoCompo;                   ///< コンポーネントの総数
  int NoVoxID;                   ///< 含まれるIDの数(Local/Global)
  int colorList[MODEL_ID_MAX+1]; ///< ボクセルモデルに含まれるIDのリスト(Global)
  CompoList*    cmp;             ///< コンポーネントリスト
  MediumList*   mat;             ///< 媒質リスト
  cpm_ParaManager *paraMngr;     ///< Cartesian Partition Maneger

public:
  /** コンストラクタ */
  VoxInfo() {
    NoVoxID = 0;
    guide = 0;
    NoBC = NoCompo = 0;
    for (int i=0; i<3; i++) size[i]=0;
    cmp = NULL;
    mat = NULL;
    paraMngr=NULL;
    
    for (int i=0; i<MODEL_ID_MAX+1; i++) colorList[i]=0;
  }
  
  /**　デストラクタ */
  ~VoxInfo() {}
  
private:
  bool chkIDinside           (int id, int* mid, int* bx);
  bool isInTable             (int MaxSize, int* cList, int target);
  
  int* allocTable(int size);
  
  int find_mat_odr(int mat_id);
  
  unsigned long countState(int id, int* mid);
  
  unsigned long count_ValidCell_OBC (const int face, const int* bv);
  
  unsigned long encodeOrder(const int order,
                            const int id, 
                            const int* mid, 
                            int* bx);
  
  unsigned long encPbit_D_IBC(const int order, 
                              const int id, 
                              const int* mid, 
                              int* bcd, 
                              int* bcp, 
                              const int deface);
  
  unsigned long encPbit_N_Binary(int* bx);
  
  unsigned long encPbit_N_Cut(int* bx, float* cut, const bool convergence);
  
  unsigned long encPbit_N_IBC(const int order, 
                              const int id, 
                              const int* mid, 
                              int* bcd, 
                              int* bcp, 
                              const int deface);
  
  unsigned long encQface(const int order, 
                         const int id, 
                         const int* mid, 
                         int* bcd, 
                         int* bh1, 
                         int* bh2, 
                         const int deface, 
                         const bool flag);
  
  unsigned long encQfaceHT_B(const int order, 
                             const int id, 
                             const int* mid, 
                             int* bcd, 
                             int* bh1, 
                             int* bh2, 
                             const int deface);
  
  unsigned long encQfaceHT_S(const int order, 
                             const int id, 
                             const int* mid, 
                             int* bcd, 
                             int* bh1, 
                             int* bh2, 
                             const int deface);
  
  unsigned long encQfaceISO_SF(const int order, 
                               const int id, 
                               const int* mid, 
                               int* bcd, 
                               int* bh1, 
                               int* bh2, 
                               const int deface);
  
  unsigned long encQfaceISO_SS(const int order, 
                               const int id, 
                               const int* mid, 
                               int* bcd, 
                               int* bh1, 
                               int* bh2, 
                               const int deface);
  
  unsigned long encVbit_IBC(const int order, 
                            const int id, 
                            const int* mid, 
                            int* bv, 
                            const int deface, 
                            int* bp);
  
  unsigned long encVbit_IBC_Cut(const int order, 
                                const int id, 
                                int* bv, 
                                int* bp, 
                                const float* cut, 
                                const int* cut_id, 
                                const float* vec, 
                                const int bc_dir);
  
  void checkColorTable       (FILE* fp, int size, int* table);
  void copyID_Prdc_Inner     (SklScalar3D<int>* d_mid, int* st, int* ed, int id, int dir);
  void encActive             (unsigned long& Lcell, unsigned long& Gcell, int* bx, const int KOS);
  void encAmask_SymtrcBC     (int face, int* bh2);
  void encHbit               (int* bh1, int* bh2);
  void encPbit               (int* bx);
  void encPbit_OBC           (int face, int* bx, string key, bool dir);
  void encQfaceSVO           (int order, int id, int* mid, int* bcd, int* bh1, int* bh2, int deface);
  void encVbit_OBC           (int face, int* bv, string key, const bool enc_sw, string chk, int* bp, const bool enc_uwd);
  void find_isolated_Fcell   (int order, int* mid, int* bx);
  void getOffset             (int* st, int* ofst);
  void resizeBVcell          (const int* st, const int* ed, int n, int* bx, int* gcbv);
  void resizeBVface          (const int* st, const int* ed, int n, int* bx, int* gcbv);
  void setInactive_Compo     (int id, int def, int* mid, int* bh1, int* bh2);
  void setAmask_InActive     (int id, int* mid, int* bh);
  void setAmask_Solid        (int* bh);
  void setAmask_Thermal      (int* bh);
  void updateGlobalIndex     (const int* st, const int* ed, int n, int* gcbv);
  

  //@fn inline int offBit(int idx, const int shift)
  //@brief idxの第shiftビットをOFFにする
  inline int offBit(int idx, const int shift) {
    return ( idx & (~(0x1<<shift)) );
  }
  
  //@fn inline int onBit(int idx, const int shift)
  //@brief idxの第shiftビットをONにする
  inline int onBit(int idx, const int shift) {
    return ( idx | (0x1<<shift) );
  }
  
  
public:
  bool chkIDconsistency      (const int m_NoMedium);
  
  int scanCell               (int *cell, const int* cid, const int ID_replace);
  
  int check_fill        (const int* mid);
  unsigned fill_cell_edge    (int* bid, int* mid, float* cut, const int tgt_id, const int solid_id);
  unsigned fill_inside       (int* mid, const int solid_id);
  unsigned setBCIndexP       (int* bcd, int* bcp, int* mid, SetBC* BC, bool isCDS=false, float* cut=NULL);
  unsigned test_opposite_cut (int* bid, int* mid, const int solid_id);
  
  unsigned long flip_InActive(unsigned long& L, 
                              unsigned long& G, 
                              const int id, 
                              const int* mid, 
                              int* bx);
  
  void adjMedium_on_GC       (const int face, SklScalar3D<int>* d_mid, const int BCtype, const int c_id, const int prdc_mode);
  void adjMediumPrdc_Inner  (SklScalar3D<int>* d_mid);
  void copyBCIbase           (int* dst, int* src);
  void countCellState        (unsigned long& Lcell, unsigned long& Gcell, int* bx, const int state);
  void countOpenAreaOfDomain (int* bx, REAL_TYPE* OpenArea);
  void fill_isolated_cells   (const int* bid, int* mid, const int isolated, const int solid_id);
  void findVIBCbbox          (const int id, const int* bv, int* st, int* ed);
  void get_Compo_Area_Cut    (int n, PolylibNS::MPIPolylib* PL);
  bool paint_first_seed      (int* mid, const int* idx, const int target);
  void printScanedCell       (FILE* fp);
  void resizeCompoBV         (int* bd, int* bv, int* bh1, int* bh2, int kos, bool isHeat, int* gcbv);
  void setAdiabatic4SF       (int* bh);
  void setBCIndexH           (int* bd, int* bh1, int* bh2, int* mid, SetBC* BC, int kos);
  void setBCIndex_base1      (int* bd, int* mid, float* cvf);
  
  void setBCIndex_base2      (int* bx, 
                              int* mid, 
                              SetBC* BC, 
                              unsigned long & Lcell, 
                              unsigned long & Gcell, 
                              const int KOS);
  
  void setBCIndexV           (int* bv, int* mid, SetBC* BC, int* bp, bool isCDS=false, float* cut=NULL, int* cut_id=NULL);
  void setCmpFraction        (CompoList* compo, int* bx, float* vf);
  void setControlVars        (int* r_size, int r_guide);
  void setNoCompo_BC         (int m_NoBC, int m_NoCompo);
  void setOBC_Cut            (SetBC* BC, float* cut);
  void setWorkList           (CompoList* m_CMP, MediumList* m_MAT);
  
  //@fn const int* getColorList() const
  //@retval colorListのポインタ
	const int* getColorList() const { return colorList; }
	
  
  //@brief CutBid5のBoundrary IDを計算
  //@note dir = (w/X_MINUS=0, e/X_PLUS=1, s/2, n/3, b/4, t/5)
  inline int get_BID5(const int dir, const int bid) {
    return ( (bid >> dir*5) & MASK_5 );
  }
  
  //@brief CutBid5のBoundrary ID設定
  //@note dir = (w/X_MINUS=0, e/X_PLUS=1, s/2, n/3, b/4, t/5)
  inline void set_BID5(int& bid, const int dir, const int s_id) {
    bid |= (s_id << (dir*5));
  }
  

  /** CPMlibのポインタをセット 
   * @param[in] m_paraMngr  
   */
  void importCPM(cpm_ParaManager* m_paraMngr)
  
  
  // ----> debug function
  void dbg_chkBCIndexD  (int* bcd, const char* fname);
  void dbg_chkBCIndexP  (int* bcd, int* bcp, const char* fname);
  void dbg_chkBCIndexV  (int* bcv, const char* fname);
};

#endif // _FB_BINVOX_H_
