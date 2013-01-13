// #################################################################
//
// CAERU Library
//
// Copyright (c) 2012-2013  All right reserved.
//
// Institute of Industrial Science, University of Tokyo, Japan. 
//
// #################################################################

//@file   BndOuter.C
//@brief  FlowBase BoundaryOuter class
//@author kero

#include "BndOuter.h"


/**
 * @brief モニタ値を保持する
 * @param [in] vv   指定速度ベクトル
 * @param [in] face 面番号
 * @param [in] mode BCの種類
 */
void BoundaryOuter::set_DomainV(const REAL_TYPE* vv, const int face, const int mode)
{
  REAL_TYPE a;
  
  switch (mode) {
    case OBC_OUTFLOW :
      dm[0] = vv[0]; // sum
      dm[1] = vv[1]; // min
      dm[2] = vv[2]; // max
      break;
      
    case OBC_SPEC_VEL:
      if ( (face==X_MINUS) || (face==X_PLUS) )
      {
        a = vv[0];
      }
      else if ( (face==Y_MINUS) || (face==Y_PLUS) )
      {
        a = vv[1];
      }
      else // Z_MINUS, Z_PLUS
      {
        a = vv[2];
      }
      
      dm[0] = a*(REAL_TYPE)valid_cell;
      dm[1] = a*(REAL_TYPE)valid_cell;
      dm[2] = a*(REAL_TYPE)valid_cell;
      break;
      
    default:
      a = vv[0]; // sumのみ
      dm[0] = a;
      dm[1] = a;
      dm[2] = a;
      break;
  }
}

// #################################################################
// ラベルを設定
void BoundaryOuter::set_Label(std::string key)
{
  label = key;
}


// #################################################################
// aliasラベルを設定する
void BoundaryOuter::set_Alias(std::string key)
{
  alias = key;
}


// #################################################################
// 境界面の有効セル数を保持
void BoundaryOuter::set_ValidCell(const int val)
{
  valid_cell = val;
}


// #################################################################
// 熱伝達境界の参照モードの保持
void BoundaryOuter::set_HTmodeRef(int key)
{
  HTref = key;
}


// #################################################################
// 熱伝達係数の保持
void BoundaryOuter::set_CoefHT(REAL_TYPE val)
{
  var1 = val;
}


// #################################################################
//温度の保持
void BoundaryOuter::set_Temp(REAL_TYPE val)
{
  var1 = val;
}


// #################################################################
// 熱流束の保持
void BoundaryOuter::set_Heatflux(REAL_TYPE val)
{
  var1 = val;
}


// #################################################################
// 周期境界のときの面の状況をセット
void BoundaryOuter::set_FaceMode(int key)
{
  Face_mode = key;
}


// #################################################################
// 周期境界のモードをセット
void BoundaryOuter::set_PrdcMode(int key)
{
  Prdc_mode = key;
}


// #################################################################
// ドライバー部分の方向をセットする
void BoundaryOuter::set_DriverDir(int key)     
{ 
  drv_dir = key;
}


// #################################################################
// ドライバー部分の方向をセットする
void BoundaryOuter::set_DriverIndex(int key)     
{ 
  drv_lid = key;
}


// #################################################################
// ガイドセルの媒質IDをセットする
void BoundaryOuter::set_GuideMedium(int key)     
{ 
  gc_medium = key;
}


// #################################################################
// 境界条件の種類をセットする
void BoundaryOuter::set_Class(const int key)     
{ 
  BCclass = key;
}


// #################################################################
// 壁面境界のモードをセットする
void BoundaryOuter::set_wallType(const int key)     
{ 
  wallType = key;
}


// #################################################################
// 熱伝達境界の種別をセット
void BoundaryOuter::set_HTmode(int key)
{
  HTmode = key;
}


// #################################################################
// 熱境界条件の種別をセット
void BoundaryOuter::set_hType(int key)
{ 
  hType  = key;
}


// #################################################################
// IN_OUT境界条件のときのBC格納番号を保持
void BoundaryOuter::set_MonRef(int key)
{
  mon_ref = key;
}


// #################################################################
// 流出境界条件のときの流出対流速度の評価モードを指定
void BoundaryOuter::set_ofv(int key)
{ 
  outType = key;
}


// #################################################################
// 外部境界の圧力指定
void BoundaryOuter::set_pType(int key)
{
  pType  = key;
}


// #################################################################
// 速度プロファイルの指定
void BoundaryOuter::set_V_Profile(const int key)
{
  v_profile  = key;
}


// #################################################################
// ベクトルのコピー
void BoundaryOuter::addVec(REAL_TYPE* vec) 
{
  nv[0] = vec[0];
  nv[1] = vec[1];
  nv[2] = vec[2];
}


// #################################################################
// メンバー変数のコピー
void BoundaryOuter::dataCopy(BoundaryOuter* src)
{
  BCclass   = src->BCclass;
  HTref     = src->HTref;
  gc_medium = src->gc_medium;
  mon_ref   = src->mon_ref;
  pType     = src->pType;
  v_profile = src->v_profile;
  hType     = src->hType;
  HTmode    = src->HTmode;
  Prdc_mode = src->Prdc_mode;
  Face_mode = src->Face_mode;
  drv_dir   = src->drv_dir;
  drv_lid   = src->drv_lid;
  p         = src->p;
  wallType  = src->wallType;
  var1      = src->var1;
  var2      = src->var2;
  valid_cell= src->valid_cell;
  outType   = src->outType;
  
  label     = src->label;
  alias     = src->alias;
  
  for (int i=0; i<3; i++) {
    nv[i] = src->nv[i];
  }
  for (int i=0; i<5; i++) {
    ca[i] = src->ca[i];
    cb[i] = src->cb[i];
  }
}
