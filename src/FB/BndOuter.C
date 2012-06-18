// #################################################################
//
// CAERU Library
//
// Copyright (c) All right reserved. 2012
//
// Institute of Industrial Science, The University of Tokyo, Japan. 
//
// #################################################################

//@file BndOuter.C
//@brief FlowBase BoundaryOuter class
//@author kero

#include "BndOuter.h"

//@fn void set_DomainV(REAL_TYPE* vv, int face, bool mode)
//@brief モニタ値を保持する
//@param vv 指定速度ベクトル
//@param face 面番号
//@param mode true-outflow, false-others(default)
//@note falseの場合は，各軸方向への寄与は成分のみ
void BoundaryOuter::set_DomainV(REAL_TYPE* vv, int face, bool mode)
{
  if ( mode ) {
    dm[0] = vv[0];
    dm[1] = vv[1];
    dm[2] = vv[2];
  }
  else {
    REAL_TYPE a;
    
    switch (face) {
      case X_MINUS:
      case X_PLUS:
        a = vv[0];
        break;
        
      case Y_MINUS:
      case Y_PLUS:
        a = vv[1];
        break;
        
      case Z_MINUS:
      case Z_PLUS:
        a = vv[2];
        break;
    }
    
    dm[0] = a; // sum
    dm[1] = a; // min
    dm[2] = a; // max
  }
}

void BoundaryOuter::set_Label(std::string key)
{
  label = key;
}

//@fn void set_ValidCell(int val)
//@brief 境界面の有効セル数を保持
void BoundaryOuter::set_ValidCell(int val)
{
  valid_cell = val;
}

//@fn void set_HTmodeRef(int key)
//@brief 熱伝達境界の参照モードの保持
void BoundaryOuter::set_HTmodeRef(int key)
{
  HTref = key;
}

//@fn void set_CoefHT(REAL_TYPE val)
//@brief 熱伝達係数の保持
void BoundaryOuter::set_CoefHT(REAL_TYPE val)
{
  var1 = val;
}

//@fn void set_Temp(REAL_TYPE val)
//@brief 温度の保持
void BoundaryOuter::set_Temp(REAL_TYPE val)
{
  var1 = val;
}

//@fn void set_Heatflux(REAL_TYPE val)
//@brief 熱流束の保持
void BoundaryOuter::set_Heatflux(REAL_TYPE val)
{
  var1 = val;
}

//@fn void set_FaceMode(int key)
//@brief 周期境界のときの面の状況をセット
void BoundaryOuter::set_FaceMode(int key)
{
  Face_mode = key;
}

//@fn void set_PrdcMode(int key)
//@brief 周期境界のモードをセット
void BoundaryOuter::set_PrdcMode(int key)
{
  Prdc_mode = key;
}

//@fn void set_DriverDir(int key) 
//@brief ドライバー部分の方向をセットする
void BoundaryOuter::set_DriverDir(int key)     
{ 
  drv_dir = key;
}

//@fn void set_DriverIndex(int key) 
//@brief ドライバー部分の方向をセットする
void BoundaryOuter::set_DriverIndex(int key)     
{ 
  drv_lid = key;
}


//@fn void set_GuideMedium(int key) 
//@brief ガイドセルの媒質IDをセットする
void BoundaryOuter::set_GuideMedium(int key)     
{ 
  gc_medium = key;
}

//@fn void set_BCtype(int key) 
//@brief 境界条件の種類をセットする
void BoundaryOuter::set_BCtype(int key)     
{ 
  BCtype = key;
}

//@fn void set_HTmode(int key)
//@brief 熱伝達境界の種別をセット
void BoundaryOuter::set_HTmode(int key)
{
  HTmode = key;
}

//@fn void set_hType(int key)
//@brief 熱境界条件の種別をセット
void BoundaryOuter::set_hType(int key)
{ 
  hType  = key;
}

//@fn void set_MonRef(int key)
//@brief IN_OUT境界条件のときのBC格納番号を保持
void BoundaryOuter::set_MonRef(int key)
{
  mon_ref = key;
}

//@fn void set_ofv(int key)
//@brief 流出，流入出境界条件のときの流出対流速度の評価モードを指定
void BoundaryOuter::set_ofv(int key)
{ 
  oflowType = key;
}

//@fn void set_pType(int key)
//@brief 外部境界の圧力指定
void BoundaryOuter::set_pType(int key)
{
  pType  = key;
}

//@fn void set_vType(int key)
//@brief 速度プロファイルの指定
void BoundaryOuter::set_vType(int key)
{
  vType  = key;
}

//@fn void addVec(REAL_TYPE* vec)
//@brief ベクトルのコピー
void BoundaryOuter::addVec(REAL_TYPE* vec) 
{
  nv[0] = vec[0];
  nv[1] = vec[1];
  nv[2] = vec[2];
}

//@fn void dataCopy(BoundaryOuter* src)
//@brief メンバー変数のコピー
void BoundaryOuter::dataCopy(BoundaryOuter* src)
{
  BCtype    = src->BCtype;
  HTref     = src->HTref;
  gc_medium = src->gc_medium;
  mon_ref   = src->mon_ref;
  pType     = src->pType;
  vType     = src->vType;
  hType     = src->hType;
  oflowType = src->oflowType;
  HTmode    = src->HTmode;
  Prdc_mode = src->Prdc_mode;
  Face_mode = src->Face_mode;
  drv_dir   = src->drv_dir;
  drv_lid   = src->drv_lid;
  p         = src->p;
  var1      = src->var1;
  var2      = src->var2;
  valid_cell= src->valid_cell;
  
  label     = src->label;
  
  for (int i=0; i<3; i++) {
    nv[i] = src->nv[i];
  }
  for (int i=0; i<5; i++) {
    ca[i] = src->ca[i];
    cb[i] = src->cb[i];
  }
}
