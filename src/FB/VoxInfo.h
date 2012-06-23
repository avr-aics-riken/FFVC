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
  unsigned size[3];              ///< 計算内部領域分割数(Local)
  unsigned guide;                ///< ガイドセルサイズ
  unsigned NoBC;                 ///< 境界条件数
  unsigned NoCompo;              ///< コンポーネントの総数
  int NoVoxID;                   ///< 含まれるIDの数(Local/Global)
  int colorList[MODEL_ID_MAX+1]; ///< ボクセルモデルに含まれるIDのリスト(Global)
  REAL_TYPE* vox_nv;             ///< ボクセルモデルに含まれるコンポーネントの法線を計算したもの(Global)
  CompoList*    cmp;             ///< コンポーネントリスト
  MediumList*   mat;             ///< 媒質リスト
  cpm_ParaManager *paraMngr;     ///< Cartesian Partition Maneger

public:
  /** コンストラクタ */
  VoxInfo() {
    NoVoxID = 0;
    guide = 0;
    NoBC = NoCompo = 0;
    for (unsigned i=0; i<3; i++) size[i]=0;
    vox_nv = NULL;
    cmp = NULL;
    mat = NULL;
    paraMngr=NULL;
    
    for (int i=0; i<MODEL_ID_MAX+1; i++) colorList[i]=0;
  }
  
  /**　デストラクタ */
  ~VoxInfo() {
    if (vox_nv)    delete [] vox_nv;
  }
  
private:
  bool chkIDinside           (unsigned id, int* mid, unsigned* bx);
  bool getXMLModel           (void);
  bool isInTable             (unsigned MaxSize, int* cList, int target);
  
  int* allocTable(unsigned size);
  
  unsigned find_mat_odr(unsigned mat_id);
  
  unsigned long countState(unsigned id, 
                           int* mid);
  
  unsigned long count_ValidCell_OBC (const int face, 
                                     const unsigned* bv);
  
  unsigned long encodeOrder(const unsigned order,
                            const unsigned id, 
                            const int* mid, 
                            unsigned* bx);
  
  unsigned long encPbit_D_IBC(const unsigned order, 
                              const unsigned id, 
                              const int* mid, 
                              unsigned* bcd, 
                              unsigned* bcp, 
                              const int deface);
  
  unsigned long encPbit_N_Binary(unsigned* bx);
  
  unsigned long encPbit_N_Cut(unsigned* bx, 
                              float* cut, 
                              const bool convergence);
  
  unsigned long encPbit_N_IBC(const unsigned order, 
                              const unsigned id, 
                              const int* mid, 
                              unsigned* bcd, 
                              unsigned* bcp, 
                              const int deface);
  
  unsigned long encQface(const unsigned order, 
                         const unsigned id, 
                         const int* mid, 
                         unsigned* bcd, 
                         unsigned* bh1, 
                         unsigned* bh2, 
                         const int deface, 
                         const bool flag);
  
  unsigned long encQfaceHT_B(const unsigned order, 
                             const unsigned id, 
                             const int* mid, 
                             unsigned* bcd, 
                             unsigned* bh1, 
                             unsigned* bh2, 
                             const int deface);
  
  unsigned long encQfaceHT_S(const unsigned order, 
                             const unsigned id, 
                             const int* mid, 
                             unsigned* bcd, 
                             unsigned* bh1, 
                             unsigned* bh2, 
                             const int deface);
  
  unsigned long encQfaceISO_SF(const unsigned order, 
                               const unsigned id, 
                               const int* mid, 
                               unsigned* bcd, 
                               unsigned* bh1, 
                               unsigned* bh2, 
                               const int deface);
  
  unsigned long encQfaceISO_SS(const unsigned order, 
                               const unsigned id, 
                               const int* mid, 
                               unsigned* bcd, 
                               unsigned* bh1, 
                               unsigned* bh2, 
                               const int deface);
  
  unsigned long encVbit_IBC(const unsigned order, 
                            const unsigned id, 
                            const int* mid, 
                            unsigned* bv, 
                            const int deface, 
                            unsigned* bp);
  
  unsigned long encVbit_IBC_Cut(const unsigned order, 
                                const unsigned id, 
                                unsigned* bv, 
                                unsigned* bp, 
                                const float* cut, 
                                const int* cut_id, 
                                const float* vec, 
                                const unsigned bc_dir);
  
  void checkColorTable       (FILE* fp, unsigned size, int* table);
  void copyID_Prdc_Inner     (SklScalar3D<int>* d_mid, int* st, int* ed, int id, int dir);
  void countNrml_from_FaceBC (unsigned n, unsigned* bx, int* cc, int& ar);
  void encActive             (unsigned long& Lcell, unsigned long& Gcell, unsigned* bx, const unsigned KOS);
  void encAmask_SymtrcBC     (int face, unsigned* bh2);
  void encHbit               (unsigned* bh1, unsigned* bh2);
  void encPbit               (unsigned* bx);
  void encPbit_OBC           (int face, unsigned* bx, string key, bool dir);
  void encQfaceSVO           (unsigned order, unsigned id, int* mid, unsigned* bcd, unsigned* bh1, unsigned* bh2, int deface);
  void encVbit_OBC           (int face, unsigned* bv, string key, const bool enc_sw, string chk, unsigned* bp, const bool enc_uwd);
  void find_isolated_Fcell   (unsigned order, int* mid, unsigned* bx);
  void getOffset             (int* st, int* ofst);
  void getQuadrant           (unsigned* q, REAL_TYPE t1, REAL_TYPE t2);
  void resizeBVcell          (const int* st, const int* ed, unsigned n, unsigned* bx, int* gcbv);
  void resizeBVface          (const int* st, const int* ed, unsigned n, unsigned* bx, int* gcbv);
  void setInactive_Compo     (unsigned id, int def, int* mid, unsigned* bh1, unsigned* bh2);
  void setAmask_InActive     (unsigned id, int* mid, unsigned* bh);
  void setAmask_Solid        (unsigned* bh);
  void setAmask_Thermal      (unsigned* bh);
  void updateGlobalIndex     (const int* st, const int* ed, unsigned n, int* gcbv);
  

  //@fn inline unsigned offBit(unsigned idx, const unsigned shift)
  //@brief idxの第shiftビットをOFFにする
  inline unsigned offBit(unsigned idx, const unsigned shift) {
    return ( idx & (~(0x1<<shift)) );
  }
  
  //@fn inline unsigned onBit(unsigned idx, const unsigned shift)
  //@brief idxの第shiftビットをONにする
  inline unsigned onBit(unsigned idx, const unsigned shift) {
    return ( idx | (0x1<<shift) );
  }
  
  
public:
  bool chkIDconsistency      (const int m_NoMedium);
  
  int scanCell               (int *cell, const int* cid, const unsigned ID_replace);
  
  unsigned check_fill        (const int* mid);
  unsigned fill_cell_edge    (int* bid, int* mid, float* cut, const int tgt_id, const int solid_id);
  unsigned fill_inside       (int* mid, const int solid_id);
  unsigned setBCIndexP       (unsigned* bcd, unsigned* bcp, int* mid, SetBC* BC, bool isCDS=false, float* cut=NULL);
  unsigned test_opposite_cut (int* bid, int* mid, const int solid_id);
  
  unsigned long flip_InActive(unsigned long& L, 
                              unsigned long& G, 
                              const unsigned id, 
                              const int* mid, 
                              unsigned* bx);
  
  void adjMedium_on_GC       (const int face, SklScalar3D<int>* d_mid, const int BCtype, const int c_id, const unsigned prdc_mode);
  void adjMediumPrdc_Inner  (SklScalar3D<int>* d_mid);
  void alloc_voxel_nv        (unsigned len);
  void cal_Compo_Area_Normal (unsigned n, unsigned* bd, unsigned* bv, unsigned* bh1, REAL_TYPE dhd, int* gi);
  void copyBCIbase           (unsigned* dst, unsigned* src);
  void countCellState        (unsigned long& Lcell, unsigned long& Gcell, unsigned* bx, const unsigned state);
  void countOpenAreaOfDomain (unsigned* bx, REAL_TYPE* OpenArea);
  void fill_isolated_cells   (const int* bid, int* mid, const int isolated, const int solid_id);
  void findVIBCbbox          (const int id, const unsigned* bv, int* st, int* ed);
  void get_Compo_Area_Cut    (unsigned n, PolylibNS::MPIPolylib* PL);
  bool paint_first_seed      (int* mid, const int* idx, const int target);
  void printScanedCell       (FILE* fp);
  void resizeCompoBV         (unsigned* bd, unsigned* bv, unsigned* bh1, unsigned* bh2, unsigned kos, bool isHeat, int* gcbv);
  void setAdiabatic4SF       (unsigned* bh);
  void setBCIndexH           (unsigned* bd, unsigned* bh1, unsigned* bh2, int* mid, SetBC* BC, unsigned kos);
  void setBCIndex_base1      (unsigned* bd, int* mid, float* cvf);
  
  void setBCIndex_base2      (unsigned* bx, 
                              int* mid, 
                              SetBC* BC, 
                              unsigned long & Lcell, 
                              unsigned long & Gcell, 
                              const unsigned KOS);
  
  void setBCIndexV           (unsigned* bv, int* mid, SetBC* BC, unsigned* bp, bool isCDS=false, float* cut=NULL, int* cut_id=NULL);
  void setCmpFraction        (CompoList* compo, unsigned* bx, float* vf);
  void setControlVars        (unsigned* r_size, unsigned r_guide);
  void setNoCompo_BC         (unsigned m_NoBC, unsigned m_NoCompo);
  void setOBC_Cut            (SetBC* BC, float* cut);
  void setWorkList           (CompoList* m_CMP, MediumList* m_MAT);
  
  //@fn const int* getColorList() const
  //@retval colorListのポインタ
	const int* getColorList() const { return colorList; }
	
  //@fn REAL_TYPE* get_vox_nv_ptr(void)
  //@brief vox_nvのポインタを返す
  REAL_TYPE* get_vox_nv_ptr(void) { return vox_nv; }
  
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
  void setPartitionManager(cpm_ParaManager* m_paraMngr)
  
  
  // ----> debug function
  void dbg_chkBCIndexD  (unsigned* bcd, const char* fname);
  void dbg_chkBCIndexP  (unsigned* bcd, unsigned* bcp, const char* fname);
  void dbg_chkBCIndexV  (unsigned* bcv, const char* fname);
};

#endif // _FB_BINVOX_H_
