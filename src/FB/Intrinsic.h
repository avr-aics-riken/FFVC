#ifndef _FB_INTRNSC_H_
#define _FB_INTRNSC_H_

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
 * @file Intrinsic.h
 * @brief FlowBase Intrinsic class Header
 * @author kero
 */

#include <math.h>
#include <fstream>

#include "cpm_Define.h"
#include "cpm_ParaManager.h"

#include "FB_Define.h"
#include "BndOuter.h"
#include "Control.h"
#include "Component.h"
#include "FBUtility.h"
#include "vec3.h"
#include "TPControl.h"
#include "Medium.h"


/** 組み込み例題のID */
enum Intrinsic_class {
  id_Users = 0,
  id_Duct,
  id_PPLT2D,
  id_SHC1D,
  id_PMT,
  id_Rect,
  id_Cylinder,
  id_Step,
  id_Polygon,
  id_Sphere
};


class Intrinsic {

protected:
  cpm_ParaManager *paraMngr;     ///< Cartesian Partition Maneger
  
public:
  unsigned size[3]; ///< 分割数
  unsigned imax;    ///< 分割数 i方向
  unsigned jmax;    ///< 分割数 j方向
  unsigned kmax;    ///< 分割数 k方向
  unsigned guide;   ///< ガイドセル数
  REAL_TYPE RefL;   ///< 代表長さ
  
  /** 次元のモード */
  enum dim_mode {
    dim_2d = 1,
    dim_3d
  };
  
  /** コンストラクタ */
  Intrinsic() { 
    for (unsigned i=0; i<3; i++) size[i]=0.0;
    imax = jmax = kmax = guide = 0;
    RefL = 0.0;
  }
  
  /**　デストラクタ */
  virtual ~Intrinsic() {}
    
  
public:
  virtual const char* getExampleName() { return NULL; };
  
  /**
   * @param[in] R       Controlクラスのポインタ
   * @param[in] tpCntl  TPControlクラスのポインタ
   * @return true-成功、false-失敗
   */
  virtual bool getTP(Control* R, TPControl* tpCntl) { return true; };
  
  virtual void initCond(REAL_TYPE* v, REAL_TYPE* p) {};
  virtual void PostInit(REAL_TYPE &checkTime, Control* R) {};
  
  
  /**
   @brief 例題名称の表示
   @param fp[in]  出力ファイルのファイルポインタ
   @param str[in] 表示文字列
   */
  virtual void printExample(FILE* fp, const char* str);
  
  
  virtual void printParaInfo(FILE* mp, FILE* fp, Control* R);
  virtual void printPara(FILE* fp, Control* R);
  virtual void setDomain(Control* R, unsigned sz[3], REAL_TYPE org[3], REAL_TYPE wth[3], REAL_TYPE pch[3]) {};
  virtual void setup(int* mid, Control* R, REAL_TYPE* G_org, const int Nmax, MediumList* mat) {};
  virtual void setup_cut(int* mid, Control* R, REAL_TYPE* G_org, const int Nmax, MediumList* mat, float* cut) {};
  
  
  /** 変数のコピー
   * @param[in] R  Controlクラスのポインタ
   */
  void setControlVars(Control* R);
  
  
  /** CPMlibのポインタをセット 
   * @param[in] m_paraMngr  
   */
  void setPartitionManager(cpm_ParaManager* m_paraMngr);
  
  
  void writeSVX(REAL_TYPE *vf, int *id, Control* R);
  void writeSVX(int *id, Control* R);
  
  
  /** MediumList中に登録されているkeyに対するIDを返す
   * @param[in] mat  MediumListクラス
   * @param[in] Namx リストの最大数
   * @param[in] key  探査するラベル
   * @return keyに対するIDを返す。発見できない場合はzero
   */
  int find_ID_from_Label(MediumList* mat, const int Nmax, const std::string key);
  
};


class IP_Users : public Intrinsic {
public:
  /** コンストラクタ */
  IP_Users() {}
  
  /**　デストラクタ */
  ~IP_Users() {}
  
public:
  
  /** 
   @brief ユーザー例題の名称を返す
   */
  const char* getExampleName()
  {
    return ("User's problem");
  }
};

#endif // _FB_INTRNSC_H_
