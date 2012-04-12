/*
 * SPHERE - Skeleton for PHysical and Engineering REsearch
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2004-2012
 *
 */

//@file ParseBC.C
//@brief FlowBase ParseBC class
//@author keno, FSI Team, VCAD, RIKEN

#include <math.h>
#include "ParseBC.h"
extern SklParaComponent* ParaCmpo;

/**
 @fn bool ParseBC::chkBCconsistency(unsigned kos)
 @brief KOSと境界条件数の整合性をチェックする
 @param kos KindOfSolver
 */
void ParseBC::chkBCconsistency(unsigned kos)
{
  if (kos == FLOW_ONLY) {
    for (int n=1; n<NoBC; n++) {
      if ( compo[n].isHBC() ) {
        stamped_printf("\tNo consistency between 'Kind_of_Solver' and 'Inner_Boundary'\n");
        Exit(0);
      }
    }
  }
  else if (kos == SOLID_CONDUCTION) {
    for (int n=1; n<NoBC; n++) {
      if ( compo[n].isVBC() ) {
        stamped_printf("\tNo consistency between 'Kind_of_Solver' and 'Inner_Boundary'\n");
        Exit(0);
      }
    }
  }
}


/**
 @fn bool ParseBC::chkID(void)
 @brief BasicBCsのidの重複を調べる
 @retval エラーコード
 */
bool ParseBC::chkID(void)
{
  int q;
  unsigned i, j;
  bool err = true;
  
  for (i=0; i<NoBaseBC; i++) {
    q = BaseBc[i].get_BC_ID();
    for (j=i+1; j<NoBaseBC; j++) {
      if ( q == BaseBc[j].get_BC_ID() ) {
        printf("Parse Error : Reduplication of Outer Boundary ID=%d\n", q);
        err = false;
      }
    }
  }
  return err;
}

/**
 @fn void ParseBC::chkKeywordIBC(const char *keyword, unsigned m)
 @brief 内部境界条件の照合を行う
 @param keyword
 @param m エントリ番号
 @note SPEC_VEL_WHは陽には現れず，getXML_IBC_SpecVel()内で登録される
 */
void ParseBC::chkKeywordIBC(const char *keyword, unsigned m)
{
  if     ( !strcasecmp(keyword, "Adiabatic") )            compo[m].setType(ADIABATIC);
  else if( !strcasecmp(keyword, "Direct_Heat_Flux") )     compo[m].setType(HEATFLUX);
  else if( !strcasecmp(keyword, "HeatTransfer_B") ) {     compo[m].setType(TRANSFER); compo[m].setHtype(HT_B); }
  else if( !strcasecmp(keyword, "HeatTransfer_N") ) {     compo[m].setType(TRANSFER); compo[m].setHtype(HT_N); }
  else if( !strcasecmp(keyword, "HeatTransfer_S") ) {     compo[m].setType(TRANSFER); compo[m].setHtype(HT_S); }
  else if( !strcasecmp(keyword, "HeatTransfer_SF") ){     compo[m].setType(TRANSFER); compo[m].setHtype(HT_SF); }
  else if( !strcasecmp(keyword, "HeatTransfer_SN") ){     compo[m].setType(TRANSFER); compo[m].setHtype(HT_SN); }
  else if( !strcasecmp(keyword, "IsoThermal") )           compo[m].setType(ISOTHERMAL);
  else if( !strcasecmp(keyword, "Radiation") )            compo[m].setType(RADIANT);
  else if( !strcasecmp(keyword, "specified_velocity"))    compo[m].setType(SPEC_VEL);
  else if( !strcasecmp(keyword, "outflow") )              compo[m].setType(OUTFLOW);
  else if( !strcasecmp(keyword, "Forcing") )              compo[m].setType(IBM_DF);
  else if( !strcasecmp(keyword, "Heat_Source") )          compo[m].setType(HEAT_SRC);
  else if( !strcasecmp(keyword, "specified_Temperature")) compo[m].setType(CNST_TEMP);
  else if( !strcasecmp(keyword, "Pressure_Loss") )        compo[m].setType(HEX);
  else if( !strcasecmp(keyword, "Fan") )                  compo[m].setType(FAN);
  else if( !strcasecmp(keyword, "Darcy") )                compo[m].setType(DARCY);
  else if( !strcasecmp(keyword, "cell_monitor") )         compo[m].setType(CELL_MONITOR);
  else if( !strcasecmp(keyword, "inactive") )             compo[m].setType(INACTIVE);
  else if( !strcasecmp(keyword, "Periodic") )             compo[m].setType(PERIODIC);
  else {
    stamped_printf("\tInvalid keyword is described '%s'\n", keyword);
    Exit(0);
  }
}

/**
 @fn void ParseBC::chkKeywordOBC(const char *keyword, unsigned m)
 @brief 外部境界条件のキーワードを照合し，コードを登録する
 @param keyword
 @param m エントリ番号
 */
void ParseBC::chkKeywordOBC(const char *keyword, unsigned m)
{
  if     ( !strcasecmp(keyword, "Wall") )               BaseBc[m].set_BCtype(OBC_WALL);
  else if( !strcasecmp(keyword, "Outflow"))             BaseBc[m].set_BCtype(OBC_OUTFLOW);
  else if( !strcasecmp(keyword, "Specified_Velocity"))  BaseBc[m].set_BCtype(OBC_SPEC_VEL);
  else if( !strcasecmp(keyword, "Symmetric"))           BaseBc[m].set_BCtype(OBC_SYMMETRIC);
  else if( !strcasecmp(keyword, "Periodic"))            BaseBc[m].set_BCtype(OBC_PERIODIC);
  else if( !strcasecmp(keyword, "Traction_Free"))       BaseBc[m].set_BCtype(OBC_TRC_FREE);
  else if( !strcasecmp(keyword, "Far_Field"))           BaseBc[m].set_BCtype(OBC_FAR_FIELD);
  else {
    stamped_printf("\tParsing error : Invalid keyword is described '%s'\n", keyword);
    Exit(0);
  }
}

/**
 @fn unsigned ParseBC::count_Outer_Cell_ID(unsigned* cid)
 @brief XMLファイルをパースして，外部境界のセルIDとその数をカウント
 @param[out] cid セルIDのリスト
 @retval 媒質数
 */
unsigned ParseBC::count_Outer_Cell_ID(unsigned* cid)
{
  const CfgElem *elemTop, *elmL1, *elmL2;
  const CfgParam* param=NULL;
  elemTop = elmL1 = elmL2 = NULL;
  
  int md = 0;
  unsigned tmp[NOFACE];
  
  // Basic Outer BCリストの読み込み
  elemTop = CF->GetTop(OUTERBND);
  elmL1 = elemTop->GetElemFirst("Basic_BCs");
  elmL2 = elmL1->GetElemFirst();
  if ( !elmL2 ) Exit(0);
  
  // 各フェイスの媒質IDを調べる
  if ( !(elmL1 = elemTop->GetElemFirst("Face_BC")) ) {
    stamped_printf("\tParsing error : Missing 'Face_BC' description\n");
    Exit(0);
  }
  
  for (int face=0; face<NOFACE; face++) {
    // faceに対するエントリを得る
    if ( !(elmL2 = selectFace(face, elmL1)) ) Exit(0);
    
    // セルIDの取得
    if ( !(param = elmL2->GetParamFirst("Guide_Cell_ID")) ) {
      stamped_printf("\tParsing error : No entory 'Guide_Cell_ID' in 'Face_BC' : %s\n", FBUtility::getDirection(face).c_str());
      Exit(0);
    }
    else {
      if ( !param->isSetID() ) {
        stamped_printf("\tParsing error : No ID section in 'Guide_Cell_ID' : %s\n", FBUtility::getDirection(face).c_str());
        Exit(0);
      }
      if ( 1 > (md=param->GetID()) ) {
        stamped_printf("\tParsing error : Invalid Outer Guide Cell ID[%d] that shoud be > 0 : %s\n", md, FBUtility::getDirection(face).c_str());
        Exit(0);
      }
      
      if ( !isIDinTable(md) ) {
        stamped_printf("\tParsing error : ID[%d] described in '%s' is not listed on 'Model_Setting'\n", md, FBUtility::getDirection(face).c_str());
        Exit(0);
      }
      
      tmp[face] = (unsigned)md;
    }
  }
  
  // 重複チェック
  unsigned count = 1;
  unsigned tgt, flag;
  cid[0] = tmp[0];
  
  for (int i=1; i<NOFACE; i++) {
    
    tgt = tmp[i]; // テストするID
    
    // 登録したmedium[]のIDをチェック，未登録なら登録し，カウントする
    flag = 0; // flag>0 で登録済み
    for (int j=0; j<i; j++) {
      if ( tgt == cid[j] ) flag++; // 登録されていたらflagをインクリメント
    }
    if ( flag == 0 ) { // 未登録の場合
      cid[count] = tgt;
      count++;
    }
  }
  
  return count;
}

/**
 @fn int ParseBC::get_BCval_int(const CfgElem *elmL, const char* key)
 @brief 境界条件の値(int型)を取得し，返す
 @param elmL
 @param key キーワード
 */
int ParseBC::get_BCval_int(const CfgElem *elmL, const char* key)
{
  int df=0;
  
  if ( !elmL->GetValue(key, &df) ) {
		stamped_printf("\tParsing error : Invalid float value for '%s'\n", key);
		Exit(0);
	}
  return df;
}

/**
 @fn REAL_TYPE ParseBC::get_BCval_real(const CfgElem *elmL, const char* key)
 @brief 境界条件の値(REAL_TYPE型)を取得し，返す
 @param elmL
 @param key キーワード
 */
REAL_TYPE ParseBC::get_BCval_real(const CfgElem *elmL, const char* key)
{
  REAL_TYPE df=0.0f;
  
  if ( !elmL->GetValue(key, &df) ) {
		stamped_printf("\tParsing error : Invalid REAL_TYPE value for '%s'\n", key);
		Exit(0);
	}
  return df;
}


/**
 @fn void ParseBC::get_Center(const CfgElem *elmL, unsigned n, const char* str, REAL_TYPE* v)
 @brief 内部境界条件の座標値を取得し，登録する
 @param elmL
 @param n オーダー
 @param v[out] ベクトルパラメータ
 @param str エラー表示用文字列
 */
void ParseBC::get_Center(const CfgElem *elmL, unsigned n, const char* str, REAL_TYPE* v)
{
  for (unsigned i=0; i<3; i++) v[i]=0.0f;
  
  if ( !elmL->GetVctValue("Center_x", "Center_y", "Center_z", &v[0], &v[1], &v[2]) ) {
    stamped_printf("\tParsing error : fail to get vec params in '%s\n", str);
    Exit(0);
  }
}


/**
 @fn void ParseBC::getDarcy(const CfgElem *elmL, unsigned n)
 @brief Darcyのパラメータを取得する
 @param elmL Forcing_Volumeのレベル
 @param n コンポーネントリストのエントリ番号
 @note
 - 透過率[m^2]は境界条件設定時に無次元化する
 */
void ParseBC::getDarcy(const CfgElem *elmL, unsigned n)
{
  if ( !elmL ) Exit(0);
  
  REAL_TYPE v[3];
	int d;
  
  for (unsigned n=0; n<3; n++) v[n]=0.0;
  
  // check number of Elem
  if ((d = elmL->GetParamSize()) != 3) {    
    stamped_printf("\tParsing error : 3 params should be found in 'Darcy' : %d\n", d);
    Exit(0);
  }
  
  // 透過率の取得
  if ( !elmL->GetVctValue("permeability_x", "permeability_y", "permeability_z", &v[0], &v[1], &v[2]) ) {
    stamped_printf("\tParsing error : fail to get permeability params in 'Darcy'\n");
    Exit(0);
  }
  compo[n].ca[0] = v[0];
  compo[n].ca[1] = v[1];
  compo[n].ca[2] = v[2];
}


/**
 @fn void ParseBC::get_Dir(const CfgElem *elmL, unsigned n, const char* str, REAL_TYPE* v)
 @brief 内部境界条件の方向ベクトル値を取得し，登録する
 @param elmL
 @param n オーダー
 @param v[out] ベクトルパラメータ
 @param str エラー表示用文字列
 */
void ParseBC::get_Dir(const CfgElem *elmL, unsigned n, const char* str, REAL_TYPE* v)
{
  for (unsigned i=0; i<3; i++) v[i]=0.0f;
  
  if ( !elmL->GetVctValue("Dir_x", "Dir_y", "Dir_z", &v[0], &v[1], &v[2]) ) {
    stamped_printf("\tParsing error : fail to get vec params in '%s\n", str);
    Exit(0);
  }
  getUnitVec(v);
}

/**
 @fn int ParseBC::getNBC(void)
 @brief INNERBNDタグ直下のElemの個数（内部境界条件数）を返す
 @retval BCの個数
 @note 境界条件数がゼロでもエラーではない
 */
int ParseBC::getNBC(void)
{  
  int nobc=0;
  
  nobc = CF->GetSize(INNERBND);
  if ( nobc == 0 ) {
    printf("\tNo Inner Boundary Conditions\n\n");
  }
  
  return nobc;
}


/**
 @fn unsigned ParseBC::getNoMaterial(void)
 @brief Model_Settingに記述されたMaterialIDの数を数える
 @retval MaterialIDの数
 @note
 - count[]はstaticなので，変数で宣言してはいけないので，newする．deleteを忘れずに
 - NoIDは'Model_Setting'に記述されたparam文の数
 */
unsigned ParseBC::getNoMaterial(void)
{
  unsigned m=0;
  bool *count = new bool[NoID+1];
  
  if( !count ) {
    fprintf(stderr, "\tMemory allocation error.(Medium counter)\n");
    return m;
  }
  
  for (int i=1; i<=NoID; i++) count[i] = true;
  
  for (int i=1; i<=NoID; i++) {
    unsigned q = iTable[i].getMatID();
    for (int j=1; j<i; j++) {
      if ( iTable[j].getMatID() == q ) count[i] = false;
    }
  }
  
  for (int i=1; i<=NoID; i++) if ( true == count[i] ) m++;
  
  if ( count ) { delete [] count; count=NULL; }
  
  return m;
}

/**
 @fn void ParseBC::get_NV(const CfgElem *elmL, unsigned n, const char* str, REAL_TYPE* v)
 @brief 内部境界条件の法線ベクトル値を取得し，登録する
 @param elmL
 @param n オーダー
 @param v[out] ベクトルパラメータ
 @param str エラー表示用文字列
 */
void ParseBC::get_NV(const CfgElem *elmL, unsigned n, const char* str, REAL_TYPE* v)
{
  for (unsigned i=0; i<3; i++) v[i]=0.0f;
  
  if ( !elmL->GetVctValue("Normal_x", "Normal_y", "Normal_z", &v[0], &v[1], &v[2]) ) {
    stamped_printf("\tParsing error : fail to get vec params in '%s\n", str);
    Exit(0);
  }
  getUnitVec(v);
}

/**
 @fn string ParseBC::getOBCstr(unsigned id)
 @brief 外部境界条件のキーワードを照合し， BCの文字列を返す
 @param id 
 */
string ParseBC::getOBCstr(unsigned id)
{
  string bc;
  if     ( id == OBC_WALL )      bc = "Wall";
  else if( id == OBC_OUTFLOW )   bc = "Outflow";
  else if( id == OBC_SPEC_VEL )  bc = "Specified_Velocity";
  else if( id == OBC_SYMMETRIC ) bc = "Symmetric";
  else if( id == OBC_PERIODIC )  bc = "Periodic";
  else if( id == OBC_TRC_FREE )  bc = "Traction_Free";
  else if( id == OBC_FAR_FIELD ) bc = "Far_Field";
  else                           bc = "";
  return bc;
}


/**
 @fn int ParseBC::getStateinTable(unsigned id)
 @brief iTable[]にidが含まれるかどうかを調べる
 @param id サーチ対象のボクセルID
 @retval 含まれれば対応するstate，そうでなければおちる
 */
int ParseBC::getStateinTable(unsigned id)
{
  for (unsigned i=1; i<=NoID; i++) {
    if ( iTable[i].getID() == id ) { return iTable[i].getState(); }
  }
  Exit(0);
  return -1;
}

/**
 @fn void ParseBC::getUnitVec(REAL_TYPE* v)
 @brief 単位ベクトルを計算して戻す
 @param [in/out] v
 */
void ParseBC::getUnitVec(REAL_TYPE* v)
{
	REAL_TYPE a;
	
	a = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
	if ( a > 0.0 ) {
		v[0] = v[0]/a;
		v[1] = v[1]/a;
		v[2] = v[2]/a;
	}
	else {
		v[0] = 0.0;
		v[1] = 0.0;
		v[2] = 0.0;
	}
}

/**
 @fn void ParseBC::getXML_Cell_Monitor(const CfgElem *elmL, unsigned n, Control* C)
 @brief XMLファイルからMonitorの設定内容をパースし，パラメータを保持する
 @param elmL コンフィギュレーションツリーのポインタ
 @param n コンポーネントリストに登録するエントリ番号のベース
 @param C Control class
 @note Referenceは，隠しコマンドに
 */
void ParseBC::getXML_Cell_Monitor(const CfgElem *elmL, unsigned n, Control* C)
{
  int nvc=0;
  const CfgElem *elmL2=NULL;
  const CfgParam* param=NULL;
  const char *pnt=NULL;
  const char *str=NULL;
  REAL_TYPE v[3];
  
  // reference 
  if ( elmL->GetValue(CfgIdt("reference"), &pnt) ) {
    if ( !strcasecmp("yes", pnt) ) {
      compo[n].setStateCellMonitor(ON);
    }
    else if ( !strcasecmp("no", pnt) ) {
      compo[n].setStateCellMonitor(OFF);
    }
    else {
      stamped_printf("\tParsing error : Invalid string value for 'reference' : %s\n", pnt);
      Exit(0);
    }
  }
  
  // 法線ベクトル
  get_NV(elmL, n, "InnerBoundary > Cell_Monitor", v);
  copyVec(compo[n].nv, v);
  
  // Variables
  if ( !( elmL2 = elmL->GetElemFirst("Variables") ) ) {
    stamped_printf("\tParsing error : No 'Variables' keyword in 'Cell_Monitor'\n");
    Exit(0);
  }
  
  // モニタする変数と数を取得
  nvc = 0;
  
  // 速度
  if ( !elmL2->GetValue("velocity", &str) ) {
    stamped_printf("\tParsing error : fail to get 'Velocity' in 'Cell_Monitor'\n");
    Exit(0);
  }
  if     ( !strcasecmp(str, "on") )  {
    compo[n].encodeVarType(var_Velocity);
    nvc++;
  }
  else if( !strcasecmp(str, "off") ) {;} // nothing
  else {
    stamped_printf("\tInvalid keyword is described for 'Velocity'\n");
    Exit(0);
  }
  
  // 圧力
  if ( !elmL2->GetValue("pressure", &str) ) {
    stamped_printf("\tParsing error : fail to get 'Pressure' in 'Cell_Monitor'\n");
    Exit(0);
  }
  if     ( !strcasecmp(str, "on") )  {
    compo[n].encodeVarType(var_Pressure);
    nvc++;
  }
  else if( !strcasecmp(str, "off") ) {;} // nothing
  else {
    stamped_printf("\tInvalid keyword is described for 'Pressure'\n");
    Exit(0);
  }
  
  // 温度
  if ( HeatProblem ) {
    if ( !elmL2->GetValue("temperature", &str) ) {
      stamped_printf("\tParsing error : fail to get 'Temperature' in 'Cell_Monitor'\n");
      Exit(0);
    }
    if     ( !strcasecmp(str, "on") )  {
      compo[n].encodeVarType(var_Temperature);
      nvc++;
    }
    else if( !strcasecmp(str, "off") ) {;} // nothing
    else {
      stamped_printf("\tInvalid keyword is described for 'Temperature'\n");
      Exit(0);
    }
  }
  
  // 全圧
  if ( C->Mode.TP == ON ) {
    if ( !elmL2->GetValue("total_pressure", &str) ) {
      stamped_printf("\tParsing error : fail to get 'Total_Pressure' in 'Cell_Monitor'\n");
      Exit(0);
    }
    if     ( !strcasecmp(str, "on") )  {
      compo[n].encodeVarType(var_TotalP);
      nvc++;
    }
    else if( !strcasecmp(str, "off") ) {;} // nothing
    else {
      stamped_printf("\tInvalid keyword is described for 'Total_Pressure'\n");
      Exit(0);
    }
  }
  
  // モニタ面に対して指定された変数の個数（モニタの個数）を取得
  compo[n].setAttrb(nvc);
}


/**
 @fn void ParseBC::getXML_IBC_Adiabatic(const CfgElem *elmL, unsigned n)
 @brief Adiabaticのパラメータを取得する
 @param elmL
 @param n コンポーネントリストのエントリ番号
 */
void ParseBC::getXML_IBC_Adiabatic(const CfgElem *elmL, unsigned n)
{
  // check number of Param
  if (elmL->GetParamSize() != 1) {    
    stamped_printf("\tParsing error : 1 param should be found in 'InnerBoundary > Adiabatic'\n");
    Exit(0);
  }
  
  // 面指定
  set_Deface(elmL, n, "InnerBoundary > Adiabatic");
  
  // zero heat flux
  if ( Unit_Param == DIMENSIONAL ) {
    compo[n].set_Heatflux( 0.0f );
  }
  else {
    stamped_printf("\tWarning: Heat condition must be given by dimensional value\n");
    Exit(0);
  }
}


/**
 @fn void ParseBC::getXML_IBC_CnstTemp(const CfgElem *elmL, unsigned n)
 @brief Const_Temperatureのパラメータを取得する
 @param elmL Const_Temperatureのレベル
 @param n コンポーネントリストのエントリ番号
 */
void ParseBC::getXML_IBC_CnstTemp(const CfgElem *elmL, unsigned n)
{
  // check number of Param
  if (elmL->GetParamSize() != 1) {    
    stamped_printf("\tParsing error : 1 param should be found in 'Specified_Temperature'\n");
    Exit(0);
  }
  
  // 温度
  REAL_TYPE tmp = get_BCval_real(elmL, "temperature");
  compo[n].set_Temp( FBUtility::convTemp2K(tmp, Unit_Temp) );
  
  if ( Unit_Param != DIMENSIONAL ) {
    stamped_printf("\tWarning: Heat condition must be given by dimensional value\n");
    Exit(0);
  }
}

/**
 @fn void ParseBC::getXML_IBC_Fan(const CfgElem *elmL, unsigned n)
 @brief Fanのパラメータを取得する
 @param elmL Forcing_Volumeのレベル
 @param n コンポーネントリストのエントリ番号
 */
void ParseBC::getXML_IBC_Fan(const CfgElem *elmL, unsigned n)
{
  if ( !elmL ) Exit(0);
  
  REAL_TYPE v[3];
  const char *str_u=NULL;
  
  // check number of Param
  if (elmL->GetParamSize() != 10) {    
    stamped_printf("\tParsing error : 1 param should be found in Heat_Volume > Heat_Generation\n");
    Exit(0);
  }
  
  // 入力単位の指定
  if ( !elmL->GetValue("unit", &str_u) ) {
		stamped_printf("\tParsing error : Invalid float value for 'unit' in 'Pressure_Loss'\n");
		Exit(0);
	}
  if ( !strcasecmp("mmaq", str_u) ) {
    compo[n].setPrsUnit(CompoList::unit_mmAq);
  }
  else if ( !strcasecmp("mmhg", str_u) ) {
    compo[n].setPrsUnit(CompoList::unit_mmHg);
  }
  else if ( !strcasecmp("pa", str_u) ) {
    compo[n].setPrsUnit(CompoList::unit_Pa);
  }
  else if ( !strcasecmp("non_dimension", str_u) ) {
    compo[n].setPrsUnit(CompoList::unit_NonDimensional);
  }
  else {
    stamped_printf("\tDescribed unit is out of scope.\n");
    Exit(0);
  }
  
  // 法線ベクトルの取得
  get_NV(elmL, n, "InnerBoundary > Pressure_Loss : Normal Vector", v);
  copyVec(compo[n].nv, v);
  
  // 中心座標の取得
  get_Center(elmL, n, "InnerBoundary > Pressure_Loss : Center", v);
  copyVec(compo[n].oc, v);
  
  // 形状パラメータ
  compo[n].depth  = get_BCval_real(elmL, "depth");
  compo[n].shp_p1 = get_BCval_real(elmL, "fan_radius");
  compo[n].shp_p2 = get_BCval_real(elmL, "boss_radius");
  
  if ( compo[n].shp_p1 <= compo[n].shp_p2 ) {
    stamped_printf("\tError : Radius of boss is greater than fan.\n");
    Exit(0);
  }
  
}


/**
 @fn void ParseBC::getXML_IBC_HeatFlux(const CfgElem *elmL, unsigned n)
 @brief Direct_Fluxのパラメータを取得する
 @param elmL
 @param n コンポーネントリストのエントリ番号
 @note [W/m^2]
 */
void ParseBC::getXML_IBC_HeatFlux(const CfgElem *elmL, unsigned n)
{
  // check number of Param
  if (elmL->GetParamSize() != 2) {    
    stamped_printf("\tParsing error : 2 params should be found in 'InnerBoundary > Direct_Heat_Flux'\n");
    Exit(0);
  }
  
  set_Deface(elmL, n, "InnerBoundary > Direct_Heat_Flux");
  
  compo[n].set_Heatflux( get_BCval_real(elmL, "Heat_Flux") );
  
  if ( Unit_Param != DIMENSIONAL ) {
    stamped_printf("\tWarning: Heat condition must be a dimensional value\n");
    Exit(0);
  }
}


/**
 @fn void ParseBC::getXML_IBC_HeatSrc(const CfgElem *elmL, unsigned n)
 @brief Heat_Generationのパラメータを取得する
 @param elmL Pointer of Configuration tree
 @param n コンポーネントリストのエントリ番号
 @note
 - D1には発熱量を保持，D2には発熱密度を保持
 */
void ParseBC::getXML_IBC_HeatSrc(const CfgElem *elmL, unsigned n)
{
  REAL_TYPE hsrc=0.0f;
  const char *str=NULL;
  
  // check number of Param
  if (elmL->GetParamSize() != 2) {    
    stamped_printf("\tParsing error : 2 params should be found in 'Heat_Source'\n");
    Exit(0);
  }
  
  // type
  if ( !elmL->GetValue(CfgIdt("type"), &str) ) {
    stamped_printf("\tParsing error : Invalid int value for 'type' in 'InnerBoundary > Heat_Source'\n");
    Exit(0);
  }
  if ( !strcasecmp("heat_release_value", str) ) {
		compo[n].set_HSRC_policy(true);
	}
	else if ( !strcasecmp("heat_generation_density", str) ) {
		compo[n].set_HSRC_policy(false);
	}
	else {
		stamped_printf("\tParsing error : Invalid string value for 'type' : %s\n", str);
		Exit(0);
	}
  
  // 放熱量
  if ( !elmL->GetValue(CfgIdt("value"), &hsrc) ) {
    stamped_printf("\tParsing error : Invalid float value for 'Heat_Generation_Density' or 'Heat_Release_Value' in 'Heat_Generation'\n");
    Exit(0);
  }
  
  if ( compo[n].isPolicy_HeatDensity() ) {
    compo[n].set_HeatDensity( hsrc ); // 発熱密度
  }
  else {
    compo[n].set_HeatValue( hsrc ); // 発熱量
  }
  if ( Unit_Param != DIMENSIONAL ) { 
    stamped_printf("\tWarning: Heat condition must be a dimensional value\n"); 
    Exit(0);
  }
  
}


/**
 @fn void ParseBC::getXML_IBC_HT_N(const CfgElem *elmL, unsigned n)
 @brief HeatTransfer_Nのパラメータを取得する
 @param elmL レベル
 @param n コンポーネントリストのエントリ番号
 */
void ParseBC::getXML_IBC_HT_N(const CfgElem *elmL, unsigned n)
{
  if ( !elmL ) Exit(0);
  
  // check number of Param
  if (elmL->GetParamSize() != 2) {    
    stamped_printf("\tParsing error : 2 params should be found in 'InnerBoundary > HeatTransfer_N'\n");
    Exit(0);
  }
  
  // 面指定
  set_Deface(elmL, n, "InnerBoundary > HeatTransfer_N");
  
  // 熱伝達係数
	compo[n].set_CoefHT( get_BCval_real(elmL, "Coef_of_Heat_Transfer") );
  
  if ( Unit_Param != DIMENSIONAL ) {
    stamped_printf("\tWarning: Heat condition must be a dimensional value\n");
    Exit(0);
  }
}

/**
 @fn void ParseBC::getXML_IBC_HT_S(const CfgElem *elmL, unsigned n)
 @brief HeatTransfer_Sのパラメータを取得する
 @param elmL
 @param n コンポーネントリストのエントリ番号
 */
void ParseBC::getXML_IBC_HT_S(const CfgElem *elmL, unsigned n)
{
  // check number of Param
  if (elmL->GetParamSize() != 3) {    
    stamped_printf("\tParsing error : 3 params should be found in 'HeatTransfer_S'\n");
    Exit(0);
  }
  
  // 面指定
  set_Deface(elmL, n, "InnerBoundary > HeatTransfer_S");
  
  // 熱伝達係数
  compo[n].set_CoefHT( get_BCval_real(elmL, "Coef_of_Heat_Transfer") );
  
  // 表面温度
  REAL_TYPE st = get_BCval_real(elmL, "Surface_Temperature");
  compo[n].set_Temp( FBUtility::convTemp2K(st, Unit_Temp) );
  
  if ( Unit_Param != DIMENSIONAL ) {
    stamped_printf("\tWarning: Heat condition must be a dimensional value\n");
    Exit(0);
  }
}

/**
 @fn void ParseBC::getXML_IBC_HT_SN(const CfgElem *elmL, unsigned n)
 @brief HeatTransfer_SNのパラメータを取得する
 @param elmL HeatTransfer_SNのレベル
 @param n コンポーネントリストのエントリ番号
 */
void ParseBC::getXML_IBC_HT_SN(const CfgElem *elmL, unsigned n)
{
  const char *str=NULL;
  
  // check number of Param
  if (elmL->GetParamSize() != 13) {    
    stamped_printf("\tParsing error : 13 params should be found in 'InnerBoundary > HeatTransfer_SN'\n");
    Exit(0);
  }
  
  // 表面温度
  REAL_TYPE st = get_BCval_real(elmL, "Surface_Temperature");
  compo[n].set_Temp( FBUtility::convTemp2K(st, Unit_Temp) );
  
  // 面指定
  set_Deface(elmL, n, "InnerBoundary > HeatTransfer_SN");
  
  // type
  if ( !elmL->GetValue(CfgIdt("Ref_Temp_Mode"), &str) ) {
    stamped_printf("\tParsing error : Invalid int value for 'Ref_Temp_Mode' in 'InnerBoundary > HeatTransfer_SN'\n");
    Exit(0);
  }
  if ( !strcasecmp("bulk_temperature", str) ) {
		compo[n].set_sw_HTmodeRef( CompoList::HT_mode_bulk );
	}
	else if ( !strcasecmp("local_temperature", str) ) {
		compo[n].set_sw_HTmodeRef( CompoList::HT_mode_local );
	}
	else {
		stamped_printf("\tParsing error : Invalid string value for 'Ref_Temp_Mode' : %s\n", str);
		Exit(0);
	}
  
  // Vertical and upper face values
  compo[n].ca[CompoList::vert_laminar_alpha]    = get_BCval_real(elmL, "vertical_laminar_alpha");
  compo[n].ca[CompoList::vert_laminar_beta]     = get_BCval_real(elmL, "vertical_laminar_beta");
  compo[n].ca[CompoList::vert_turbulent_alpha]  = get_BCval_real(elmL, "vertical_turbulent_alpha");
  compo[n].ca[CompoList::vert_turbulent_beta]   = get_BCval_real(elmL, "vertical_turbulent_beta");
  compo[n].ca[CompoList::vert_Ra_critial]       = get_BCval_real(elmL, "vertical_Ra_critial");
  
  // Lower face values
  compo[n].cb[CompoList::lower_laminar_alpha]   = get_BCval_real(elmL, "lower_laminar_alpha");
  compo[n].cb[CompoList::lower_laminar_beta]    = get_BCval_real(elmL, "lower_laminar_beta");
  compo[n].cb[CompoList::lower_turbulent_alpha] = get_BCval_real(elmL, "lower_turbulent_alpha");
  compo[n].cb[CompoList::lower_turbulent_beta]  = get_BCval_real(elmL, "lower_turbulent_beta");
  compo[n].cb[CompoList::lower_Ra_critial]      = get_BCval_real(elmL, "lower_Ra_critial");
  
  if ( Unit_Param != DIMENSIONAL ) {
    stamped_printf("\tWarning: Heat condition must be a dimensional value\n");
    Exit(0);
  }
}

/**
 @fn void ParseBC::getXML_IBC_HT_SF(const CfgElem *elmL, unsigned n)
 @brief HeatTransfer_SFのパラメータを取得する
 @param elmL HeatTransfer_SFのレベル
 @param n コンポーネントリストのエントリ番号
 */
void ParseBC::getXML_IBC_HT_SF(const CfgElem *elmL, unsigned n)
{
  const char *str=NULL;
  
  // check number of Param
  if (elmL->GetParamSize() != 6) {    
    stamped_printf("\tParsing error : 6 params should be found in 'InnerBoundary > HeatTransfer_SF'\n");
    Exit(0);
  }
  
  // 表面温度
  REAL_TYPE st = get_BCval_real(elmL, "Surface_Temperature");
  compo[n].set_Temp( FBUtility::convTemp2K(st, Unit_Temp) );
  
  // 面指定
  set_Deface(elmL, n, "InnerBoundary > HeatTransfer_SF");
  
  // type
  if ( !elmL->GetValue(CfgIdt("Ref_Temp_Mode"), &str) ) {
    stamped_printf("\tParsing error : Invalid int value for 'Ref_Temp_Mode' in 'InnerBoundary > HeatTransfer_SF'\n");
    Exit(0);
  }
  if ( !strcasecmp("bulk_temperature", str) ) {
		compo[n].set_sw_HTmodeRef( CompoList::HT_mode_bulk );
	}
	else if ( !strcasecmp("local_temperature", str) ) {
		compo[n].set_sw_HTmodeRef( CompoList::HT_mode_local );
	}
	else {
		stamped_printf("\tParsing error : Invalid string value for 'type' : %s\n", str);
		Exit(0);
	}
  
  // coefficients
  compo[n].ca[CompoList::alpha] = get_BCval_real(elmL, "alpha");
  compo[n].ca[CompoList::beta]  = get_BCval_real(elmL, "beta");
  compo[n].ca[CompoList::gamma] = get_BCval_real(elmL, "gamma");
  
  if ( Unit_Param != DIMENSIONAL ) {
    stamped_printf("\tWarning: Heat condition must be a dimensional value\n");
    Exit(0);
  }
}

/**
 @fn void ParseBC::getXML_IBC_HT_B(const CfgElem *elmL, unsigned n)
 @brief HeatTransfer_Bのパラメータを取得する
 @param elmL HeatTransfer_Bのレベル
 @param n コンポーネントリストのエントリ番号
 */
void ParseBC::getXML_IBC_HT_B(const CfgElem *elmL, unsigned n)
{
  // check number of Param
  if (elmL->GetParamSize() != 3) {    
    stamped_printf("\tParsing error : 3 params should be found in 'InnerBoundary > HeatTransfer_B'\n");
    Exit(0);
  }
  
  // 面指定
  set_Deface(elmL, n, "InnerBoundary > HeatTransfer_B");
  
  // 熱伝達係数
  compo[n].set_CoefHT( get_BCval_real(elmL, "Coef_of_Heat_Transfer") );
  
  // バルク温度
  REAL_TYPE st = get_BCval_real(elmL, "Bulk_Temperature");
  compo[n].set_Temp( FBUtility::convTemp2K(st, Unit_Temp) );
  
  if ( Unit_Param != DIMENSIONAL ) {
    stamped_printf("\tWarning: Heat condition must be given by dimensional value\n");
    Exit(0);
  }
}


/**
 @fn void ParseBC::getXML_IBC_IBM_DF(const CfgElem *elmL, unsigned n)
 @brief Direct Forcingのパラメータを取得する
 @param elmL 
 @param n コンポーネントリストのエントリ番号
 */
void ParseBC::getXML_IBC_IBM_DF(const CfgElem *elmL, unsigned n)
{
	int d;
  REAL_TYPE v[3];
  
  // check number of Elem
  if ((d = elmL->GetParamSize()) != 4) {    
    stamped_printf("\tParsing error : 4 params should be found in 'InnerBoundary > Forcing' : %d\n", d);
    Exit(0);
  }
  
  // 法線ベクトル
  get_NV(elmL, n, "InnerBoundary > Forcing", v);
  copyVec(compo[n].nv, v);
  
  // Velocity
  REAL_TYPE ct = get_BCval_real(elmL, "Velocity");
  if ( Unit_Param == DIMENSIONAL ) {
    compo[n].set_Velocity( ct );
  }
  else {
    compo[n].set_Velocity( ct * RefVelocity );
  }
}



/**
 @fn void ParseBC::getXML_IBC_IsoTherm(const CfgElem *elmL, unsigned n)
 @brief XMLファイルから境界条件IsoThermalのパラメータを取得し保持する
 @param elmL IsoThermalのレベル
 @param n コンポーネントリストのエントリ番号
 */
void ParseBC::getXML_IBC_IsoTherm(const CfgElem *elmL, unsigned n)
{
  // check number of Param
  if (elmL->GetParamSize() != 2) {    
    stamped_printf("\tParsing error : 2 params should be found in 'InnerBoundary > IsoThermal'\n");
    Exit(0);
  }
  
  // 表面温度
  REAL_TYPE tmp = get_BCval_real(elmL, "temperature");
	compo[n].set_Temp( FBUtility::convTemp2K(tmp, Unit_Temp) );
	
  // 面指定
  set_Deface(elmL, n, "InnerBoundary > IsoThermal");
  
  if ( Unit_Param != DIMENSIONAL ) {
    stamped_printf("\tWarning: Heat condition must be given by dimensional value\n");
    Exit(0);
  }
}


/**
 @fn void ParseBC::getXML_IBC_Outflow(const CfgElem *elmL, unsigned n)
 @brief 内部の流出境界のパラメータを取得する
 @param elmL レベル
 @param n コンポーネントリストのエントリ番号
 */
void ParseBC::getXML_IBC_Outflow(const CfgElem *elmL, unsigned n)
{
  int def;
	const char *str=NULL;
  REAL_TYPE ct;
  REAL_TYPE v[3];
  
  // 圧力境界のタイプ default
  compo[n].set_sw_P_BCtype( P_GRAD_ZERO );
  
  // Hidden parameter
  if ( elmL->GetValue(CfgIdt("pressure_type"), &str) ) {
    //printf("\tParsing error : fail to get 'Pressure_Type' in 'InnerBoundary > Outflow'\n");
    //Exit(0);
    if ( !strcasecmp("dirichlet", str) ) {
      compo[n].set_sw_P_BCtype( P_DIRICHLET );
    }
    else if ( !strcasecmp("grad_zero", str) ) {
      compo[n].set_sw_P_BCtype( P_GRAD_ZERO );
    }
    else {
      printf("\tParsing error : Invalid string value for 'Pressure_Type' : %s\n", str);
      Exit(0);
    }
    if ( compo[n].get_sw_P_BCtype() == P_DIRICHLET ) {
      compo[n].set_Pressure( get_BCval_real(elmL, "pressure_value") );
    }
  }
  
  // 面指定
  set_Deface(elmL, n, "InnerBoundary > Outflow");
  
  // 流出速度のタイプ
  compo[n].setOutflowType(V_AVERAGE);
  
  // 法線ベクトルの取得
  get_NV(elmL, n, "InnerBoundary > Outflow : Normal Vector", v);
  copyVec(compo[n].nv, v);
  
  /*
   if ( !elmL->GetValue(CfgIdt("velocity_type"), &str) ) {
   printf("\tParsing error : fail to get 'Velocity_Type' in 'InnerBoundary > Outflow'\n");
   Exit(0);
   }
   if ( !strcasecmp("average", str) ) {
   compo[n].flag = V_AVERAGE;
   }
   else if ( !strcasecmp("minmax", str) ) {
   compo[n].flag = V_MINMAX;
   }
   else {
   printf("\tParsing error : Invalid string value for 'Velocity_Type' : %s\n", str);
   Exit(0);
   }*/
}

/**
 @fn void ParseBC::getXML_IBC_Periodic(const CfgElem *elmL, unsigned n)
 @brief 内部の周期境界のパラメータを取得する
 @param elmL レベル
 @param n コンポーネントリストのエントリ番号
 */
void ParseBC::getXML_IBC_Periodic(const CfgElem *elmL, unsigned n)
{
  int dir=0;
  REAL_TYPE ct=0.0;
  const char *str=NULL;
  
  // 上流側の方向
  if ( !elmL->GetValue(CfgIdt("upstream_direction"), &str) ) {
		printf("\tParsing error : fail to get 'upstream_direction' in 'InnerBoundary > Periodic'\n");
    Exit(0);
  }
  if ( !strcasecmp("x_minus", str) ) {
		dir = X_MINUS;
	}
  else if ( !strcasecmp("x_plus", str) ) {
    dir = X_PLUS;
  }
  else if ( !strcasecmp("y_minus", str) ) {
    dir = Y_MINUS;
  }
  else if ( !strcasecmp("y_plus", str) ) {
    dir = Y_PLUS;
  }
  else if ( !strcasecmp("z_minus", str) ) {
    dir = Z_MINUS;
  }
  else if ( !strcasecmp("z_plus", str) ) {
    dir = Z_PLUS;
  }
  else {
    stamped_printf("\tParsing error : Invalid direction in 'InnerBoundary > Periodic'\n");
    Exit(0);
  }
	compo[n].setPeriodicDir((unsigned)dir);
  
  // 圧力差
  if ( elmL->GetValue(CfgIdt("pressure_difference"), &ct)) {
    compo[n].ca[0] = ct;
  }
  else {
    printf("\tParsing error : Invalid value of 'Pressure difference' in 'InnerBoundary > Periodic'\n");
    Exit(0);
  }
  
  // 面指定
  set_Deface(elmL, n, "InnerBoundary > Periodic");
}

/**
 @fn void ParseBC::getXML_IBC_PrsLoss(const CfgElem *elmL, unsigned n)
 @brief HeatExchangerのパラメータを取得する
 @param elmL コンフィギュレーションツリーのポインタ
 @param n コンポーネントリストのエントリ番号
 @note この時点ではRefDensityの値が未定なので，あとでパラメータ処理
 @see Control::setParameters()
 */
void ParseBC::getXML_IBC_PrsLoss(const CfgElem *elmL, unsigned n)
{
  const char *str=NULL, *str_u=NULL;
  REAL_TYPE v[3], ct;
  
  // check number of Elem
  if ( elmL->GetParamSize() != 20) {    
    stamped_printf("\tParsing error : 11 params should be found in 'Pressure_Loss'\n");
    Exit(0);
  }
  
  // 入力単位の指定
  if ( !elmL->GetValue("unit", &str_u) ) {
		stamped_printf("\tParsing error : Invalid float value for 'unit' in 'Pressure_Loss'\n");
		Exit(0);
	}
  if ( !strcasecmp("mmaq", str_u) ) {
    compo[n].setPrsUnit(CompoList::unit_mmAq);
  }
  else if ( !strcasecmp("mmhg", str_u) ) {
    compo[n].setPrsUnit(CompoList::unit_mmHg);
  }
  else if ( !strcasecmp("pa", str_u) ) {
    compo[n].setPrsUnit(CompoList::unit_Pa);
  }
  else if ( !strcasecmp("non_dimension", str_u) ) {
    compo[n].setPrsUnit(CompoList::unit_NonDimensional);
  }
  else {
    stamped_printf("\tDescribed unit is out of scope.\n");
    Exit(0);
  }
  
  // 法線ベクトルの取得
  get_NV(elmL, n, "InnerBoundary > Pressure_Loss : Normal Vector", v);
  copyVec(compo[n].nv, v);
  
  // 方向ベクトルの取得
  get_Dir(elmL, n, "InnerBoundary > Pressure_Loss : Dir. Vector", v);
  copyVec(compo[n].dr, v);
  
  // 中心座標の取得
  get_Center(elmL, n, "InnerBoundary > Pressure_Loss : Center", v);
  copyVec(compo[n].oc, v);
  
  // 形状パラメータ
  compo[n].depth  = get_BCval_real(elmL, "depth");
  compo[n].shp_p1 = get_BCval_real(elmL, "width");
  compo[n].shp_p2 = get_BCval_real(elmL, "height");
  
  // 圧力損失パラメータ
	compo[n].ca[0] = get_BCval_real(elmL, "c1");
	compo[n].ca[1] = get_BCval_real(elmL, "c2");
	compo[n].ca[2] = get_BCval_real(elmL, "c3");
	compo[n].ca[3] = get_BCval_real(elmL, "c4");
	compo[n].ca[4] = get_BCval_real(elmL, "u_threshold");
	compo[n].ca[5] = get_BCval_real(elmL, "thickness");
  
  // 熱交換器の方向強制オプション
  if ( !elmL->GetValue("vector", &str) ) {
		stamped_printf("\tParsing error : Invalid string for 'vector' in 'Pressure_Loss'\n");
		Exit(0);
	}
  if ( !strcasecmp("directional", str) ) {
    compo[n].set_sw_HexDir( ON );
  }
  else {
    compo[n].set_sw_HexDir( OFF );
  }
}


/**
 @fn void ParseBC::getXML_IBC_SpecVel(const CfgElem *elmL, unsigned n)
 @brief 内部の流入境界のパラメータを取得する
 @param elmL レベル
 @param n コンポーネントリストのエントリ番号
 @note Control::setparameters()でcompo[].ca[]に値をセットする
 */
void ParseBC::getXML_IBC_SpecVel(const CfgElem *elmL, unsigned n)
{
  REAL_TYPE ct, vel;
  const char* str=NULL;
  REAL_TYPE v[3];
  
  // 速度指定タイプ
  compo[n].set_sw_V_profile( getXML_Vel_profile(elmL, "InnerBoundary > Specified_Velocity") );
  
  // 速度の指定モードの特定
  if ( !elmL->GetValue(CfgIdt("specified_type"), &str) ) {
    stamped_printf("\tParsing error : Invalid Specified_Type in 'InnerBoundary > Specified_Type'\n");
    Exit(0);
  }
  if ( !strcasecmp("velocity", str) ) {
		compo[n].set_VBC_policy(true);
	}
	else if ( !strcasecmp("massflow", str) ) {
		compo[n].set_VBC_policy(false);
	}
	else {
		printf("\tParsing error : Invalid string value '%s' for 'Specified_Type'\n", str);
		Exit(0);
	}
  
  // 指定値の取得
  if ( !elmL->GetValue(CfgIdt("specified_value"), &ct) ) {
    stamped_printf("\tParsing error : Invalid value in 'InnerBoundary > Specified_Value'\n");
    Exit(0);
  }
  if ( compo[n].isPolicy_Massflow() ) { // 流量指定の場合
    vel = ct; // 有次元でも無次元でも，モデルの断面積を計算して，後ほどパラメータ計算 >> Control::setParameters()
  }
  else {
    vel = ( Unit_Param == DIMENSIONAL ) ? ct : ct * RefVelocity;
  }
  
  // 面指定
  set_Deface(elmL, n, "InnerBoundary > Specified_Velocity");
  
  // 法線ベクトル
  get_NV(elmL, n, "InnerBoundary > Specified_Velocity", v);
  copyVec(compo[n].nv, v);
  
  // 速度パラメータの読み込み
  getXML_Vel_Params(elmL, compo[n].get_sw_V_profile(), compo[n].ca, vel, "InnerBoundary > Specified_Velocity");
  
  // heat problem
  if ( HeatProblem ) {
    if ( !elmL->GetValue(CfgIdt("temperature"), &ct) ) {
      stamped_printf("\tParsing error : fail to get 'Temperature' in 'InnerBoundary > SPEC_VEL'\n");
      Exit(0);
    }
    else {
      compo[n].set_Temp( FBUtility::convTemp2K(ct, Unit_Temp) );
      if ( Unit_Param != DIMENSIONAL ) {
        stamped_printf("\tWarning: Heat condition must be a dimensional value\n");
        Exit(0);
      }
      compo[n].setType(SPEC_VEL_WH); // SPEC_VELから変更
    }
  }
  
  // 境界条件位置の指定
  if ( isCDS ) {
    if ( !elmL->GetValue(CfgIdt("BC_position"), &str) ) {
      stamped_printf("\tParsing error : Invalid Specified_Type in 'InnerBoundary > BC_position'\n");
      Exit(0);
    }
    if ( !strcasecmp("same_direction", str) ) {
      compo[n].setBClocation(CompoList::same_direction);
    }
    else if ( !strcasecmp("opposite_direction", str) ) {
      compo[n].setBClocation(CompoList::opposite_direction);
    }
    else {
      printf("\tParsing error : Invalid string value '%s' for 'BC_position'\n", str);
      Exit(0);
    }
  }

}

/**
 @fn void ParseBC::getXML_IBC_Radiant(const CfgElem *elmL, unsigned n)
 @brief XMLファイルから境界条件Radiantのパラメータを取得し保持する
 @todo
 - 境界条件自体は未実装
 */
void ParseBC::getXML_IBC_Radiant(const CfgElem *elmL, unsigned n)
{
  // check number of Param
  if (elmL->GetParamSize() == 0) {    
    stamped_printf("\tParsing error : Missing Elements in 'InnerBoundary > Radiant'\n");
    Exit(0);
  }
  
  // 面指定
  set_Deface(elmL, n, "InnerBoundary > Radiant");
  
  // 係数
  compo[n].set_CoefRadEps( get_BCval_real(elmL, "epsilon") );
  
  // 射出率
  compo[n].set_CoefRadPrj( get_BCval_real(elmL, "projection") );
  
  if ( Unit_Param != DIMENSIONAL ) {
    stamped_printf("\tWarning: Heat condition must be a dimensional value\n");
    Exit(0);
  }
}


/**
 @fn void ParseBC::getXML_Medium_InitTemp(void)
 @brief 温度計算の場合の各媒質の初期値を取得する
 */
void ParseBC::getXML_Medium_InitTemp(void)
{
  const CfgElem *elemTop=NULL, *elmL1=NULL;
  const CfgParam* param=NULL;
  int id;
  const char* p=NULL;
  unsigned Cell_state;
  REAL_TYPE ct;
  
  // Check Model_Setting section
  elemTop = CF->GetTop(PARAMETER);
  if( !(elmL1 = elemTop->GetElemFirst("Init_Temp_of_Medium")) ) {
    stamped_printf("\tParsing error : Missing the section of 'Init_Temp_of_Medium'\n");
    Exit(0);
  }
  
  unsigned m_no_medium = NoCompo - NoBC;
  
  // load statement list
  param = elmL1->GetParamFirst();
  for (unsigned i=1; i<=m_no_medium; i++) {
    
    // state
    p = param->GetName();
    if ( !strcasecmp(p, "Solid") ) {
      Cell_state = SOLID;
    }
    else if ( !strcasecmp(p, "Fluid") ) {
      Cell_state = FLUID;
    }
    else {
      stamped_printf("\tParsing error : No valid keyword [SOLID/FLUID] in 'Init_Temp_of_Medium'\n");
      Exit(0);
    }
    
    // ID
    if ( !param->isSetID() ) {
      stamped_printf("\tParsing error : No ID for statement in 'Init_Temp_of_Medium'\n");
      Exit(0);
    }
    if ( -1 == (id=param->GetID()) ) {
      stamped_printf("\tParsing error : No valid ID for statement in 'Init_Temp_of_Medium'\n");
      Exit(0);
    }
    
    // マッチング
    bool m_flag = false;
    for (unsigned m=NoBC+1; m<=NoCompo; m++) {
      
      if ( compo[m].getID() == (unsigned)id ) {
        m_flag = true;
        if ( compo[m].getState() != Cell_state ) {
          stamped_printf("\tError : Inconsistent the cell state between 'Model_Setting' and 'Init_Temp_of_Medium' : Medium ID=%d\n", id);
          Exit(0);
        }
        
        if ( !param->GetData(&ct) ) {
          stamped_printf("\tParsing error : No initial temperature in 'Init_Temp_of_Medium'\n");
          Exit(0);
        }
        compo[m].setInitTemp( FBUtility::convTemp2K(ct, Unit_Temp) );
        break;
      }
    }
    //check
    if ( !m_flag ) {
      stamped_printf("\tError : could not find ID=%d in ComponentList\n", id);
      Exit(0);
    }
    
    param = elmL1->GetParamNext(param);
  }
}

/**
 @fn void ParseBC::getXML_Model(void)
 @brief XMLに記述されたモデル情報を取得し，ワークテーブルiTableを確保する
 @pre ParseBC::setControlVars()内 scanXMLmodel()でNoIDを確定
 @note
    - Model_Settingの情報をXMLから取得し，iTableに保持
    - 流体の媒質は少なくとも一つは必要
    - iTableの数の上限は，NoID
 */
void ParseBC::getXML_Model(void)
{
  const CfgElem *elemTop=NULL, *elmL1=NULL;
  const CfgParam* param=NULL;
  int id;
  const char* p=NULL;
  const char* pnt;
  
  // Check Model_Setting section
  elemTop = CF->GetTop(STEER);
  if( !(elmL1 = elemTop->GetElemFirst("Model_Setting")) ) {
    stamped_printf("\tParsing error : Missing the section of 'Model_Setting'\n");
    Exit(0);
  }

  // load statement list
  param = elmL1->GetParamFirst();
  for (int i=1; i<=NoID; i++) {

    // state
    p = param->GetName();
    if ( !strcasecmp(p, "Solid") ) {
      iTable[i].setState(SOLID);
    }
    else if ( !strcasecmp(p, "Fluid") ) {
      iTable[i].setState(FLUID);
    }
    else {
      stamped_printf("\tParsing error : No valid keyword [SOLID/FLUID] in 'Model_Setting'\n");
      Exit(0);
    }

    // ID
    if ( !param->isSetID() ) {
      stamped_printf("\tParsing error : No ID for statement in 'Model_Setting'\n");
      Exit(0);
    }
    if ( -1 == (id=param->GetID()) ) {
      stamped_printf("\tParsing error : No valid ID for statement in 'Model_Setting'\n");
      Exit(0);
    }
    iTable[i].setID((unsigned)id);

    // MatID
    if ( !param->GetData( &id ) ) {
      stamped_printf("\tParsing error : No valid Medium ID for statement in 'Model_Setting'\n");
      Exit(0);
    }
    iTable[i].setMatID((unsigned)id);
    
    // label
    pnt = NULL;
    if ( !(pnt = param->GetComment()) ) {
      stamped_printf("\tNo comment for statement in 'Model_Setting'\n");
    }
    else {
      iTable[i].setLabel(pnt);
    }

    param = elmL1->GetParamNext(param);
  }

  // sort Table
  {
    int      s, i, j;
    unsigned d, m;
    char label[LABEL];
    for (i=1; i<=NoID; i++) {
      for (j=NoID; j>i; j--) {
        if ( iTable[j].getID() < iTable[j-1].getID() ) {
          s = iTable[j].getState();
          d = iTable[j].getID();
          m = iTable[j].getMatID();
          strcpy(label, iTable[j].getLabel());
        
          iTable[j].setState( iTable[j-1].getState() );
          iTable[j].setID   ( iTable[j-1].getID()  );
          iTable[j].setMatID( iTable[j-1].getMatID() );
          iTable[j].setLabel( iTable[j-1].getLabel() );
        
          iTable[j-1].setState( s );
          iTable[j-1].setID   ( d );
          iTable[j-1].setMatID( m );
          iTable[j-1].setLabel( label );
        }
      }
    }
  };
}

/**
 @fn void ParseBC::getXML_OBC_FarField(const CfgElem *elmL, unsigned n)
 @brief 外部境界の遠方境界のパラメータを取得する
 @param elmL 
 @param n 面番号
 */
void ParseBC::getXML_OBC_FarField(const CfgElem *elmL, unsigned n)
{
  REAL_TYPE ct;
 
  BaseBc[n].set_pType(P_GRAD_ZERO);
  BaseBc[n].p = 0.0; // ダミー値
  
  // 外部雰囲気温
  if ( HeatProblem ) {
    if ( elmL->GetValue(CfgIdt("ambient_temperature"), &ct) ) {
      BaseBc[n].set_Temp( FBUtility::convTemp2K(ct, Unit_Temp) );
    }
    else {
      stamped_printf("\tParsing error : fail to get 'ambient_temperature' in 'Basic_BCs > Far_Field'\n");
      Exit(0);
    }
    if ( Unit_Param != DIMENSIONAL ) {
      stamped_printf("\tError: Heat condition must be given by dimensional value\n");
      Exit(0);
    }
  }
  
  /* 圧力の値
   if ( elmL->GetValue(CfgIdt("pressure_value"), &ct) ) {
   BaseBc[n].p = ct;
   }
   else {
   stamped_printf("\tParsing error : fail to get 'pressure_value' in 'Basic_BCs > Traction_Free'\n");
   Exit(0);
   } */
}

/*
 @fn void ParseBC::getXML_OBC_HT(const CfgElem *elmL, unsigned n, const char* kind)
 @brief 外部の壁面熱伝達境界のパラメータを取得する
 @param elmL 
 @param n 面番号
 @param kind 熱伝達境界の種類
 */ 
void ParseBC::getXML_OBC_HT(const CfgElem *elmL, unsigned n, const char* kind)
{
  const CfgParam* param=NULL;
	const char *str=NULL;
  REAL_TYPE ct;
  
  if ( !strcasecmp(kind, "HeatTransfer_B") ) {
    BaseBc[n].set_HTmode(HT_B);
    BaseBc[n].set_CoefHT( get_BCval_real(elmL, "Coef_of_Heat_Transfer") );
    ct = get_BCval_real(elmL, "Bulk_Temperature");
    BaseBc[n].set_Temp( FBUtility::convTemp2K(ct, Unit_Temp) );
  }
  else if ( !strcasecmp(kind, "HeatTransfer_N") ) {
    BaseBc[n].set_HTmode(HT_N);
    BaseBc[n].set_CoefHT( get_BCval_real(elmL, "Coef_of_Heat_Transfer") );
  }
  else if ( !strcasecmp(kind, "HeatTransfer_S") ) {
    BaseBc[n].set_HTmode(HT_S);
    BaseBc[n].set_CoefHT( get_BCval_real(elmL, "Coef_of_Heat_Transfer") );
    ct = get_BCval_real(elmL, "Surface_temperature");
    BaseBc[n].set_Temp( FBUtility::convTemp2K(ct, Unit_Temp) );
  }
  else if ( !strcasecmp(kind, "HeatTransfer_SF") ) {
    BaseBc[n].set_HTmode(HT_SF);
    BaseBc[n].set_Temp( get_BCval_real(elmL, "Surface_temperature") );
    if ( !elmL->GetValue(CfgIdt("ref_temp_mode"), &str) ) {
      stamped_printf("\tParsing error : Invalid int value for 'ref_temp_mode' in 'Basic_BCs > wall'\n");
      Exit(0);
    }
    if ( !strcasecmp("bulk_temperature", str) ) {
      BaseBc[n].set_HTmodeRef( CompoList::HT_mode_bulk );
    }
    else if ( !strcasecmp("local_temperature", str) ) {
      BaseBc[n].set_HTmodeRef( CompoList::HT_mode_local );
    }
    else {
      stamped_printf("\tParsing error : Invalid string value for 'ref_temp_mode' : %s\n", str);
      Exit(0);
    }
    // coefficients
    BaseBc[n].ca[0] = get_BCval_real(elmL, "alpha");
    BaseBc[n].ca[1] = get_BCval_real(elmL, "beta");
    BaseBc[n].ca[2] = get_BCval_real(elmL, "gamma");
  }
  else if ( !strcasecmp(kind, "HeatTransfer_SN") ) {
    BaseBc[n].set_HTmode(HT_SN);
    BaseBc[n].set_Temp( get_BCval_real(elmL, "Surface_temperature") );
    // reference mode
    if ( !elmL->GetValue(CfgIdt("ref_temp_mode"), &str) ) {
      stamped_printf("\tParsing error : Invalid int value for 'ref_temp_mode' in 'Basic_BCs > wall'\n");
      Exit(0);
    }
    if ( !strcasecmp("bulk_temperature", str) ) {
      BaseBc[n].set_HTmodeRef( CompoList::HT_mode_bulk );
    }
    else if ( !strcasecmp("local_temperature", str) ) {
      BaseBc[n].set_HTmodeRef( CompoList::HT_mode_local );
    }
    else {
      stamped_printf("\tParsing error : Invalid string value for 'ref_temp_mode' : %s\n", str);
      Exit(0);
    }
    // Vertical and upper face values
    BaseBc[n].ca[0] = get_BCval_real(elmL, "vertival_laminar_alpha");
    BaseBc[n].ca[1] = get_BCval_real(elmL, "vertival_laminar_beta");
    BaseBc[n].ca[2] = get_BCval_real(elmL, "vertival_turbulent_alpha");
    BaseBc[n].ca[3] = get_BCval_real(elmL, "vertival_turbulent_beta");
    BaseBc[n].ca[4] = get_BCval_real(elmL, "vertival_Ra_critial");
    
    // Lower face values
    BaseBc[n].cb[0] = get_BCval_real(elmL, "lower_laminar_alpha");
    BaseBc[n].cb[1] = get_BCval_real(elmL, "lower_laminar_beta");
    BaseBc[n].cb[2] = get_BCval_real(elmL, "lower_turbulent_alpha");
    BaseBc[n].cb[3] = get_BCval_real(elmL, "lower_turbulent_beta");
    BaseBc[n].cb[4] = get_BCval_real(elmL, "lower_Ra_critial");
  }
  else {
    stamped_printf("\tParsing error : fail to get HeatTransfer Type in 'Basic_BCs > wall'\n");
    Exit(0);
  }
}


/**
 @fn void ParseBC::getXML_OBC_InOut(const CfgElem *elmL, unsigned n)
 @brief 外部境界の流出条件のパラメータを取得する
 @param elmL 
 @param n 面番号
 */
void ParseBC::getXML_OBC_InOut(const CfgElem *elmL, unsigned n)
{
	const char *str=NULL;
  REAL_TYPE ct;
  
  // 流出速度のタイプ
  if ( !elmL->GetValue(CfgIdt("velocity_type"), &str) ) {
    stamped_printf("\tParsing error : fail to get 'Velocity_Type' in 'Basic_BCs > IN_OUT'\n");
    Exit(0);
  }
  if ( !strcasecmp("average", str) ) {
		BaseBc[n].set_ofv(V_AVERAGE);
	}
	else if ( !strcasecmp("minmax", str) ) {
		BaseBc[n].set_ofv(V_MINMAX);
	}
	else {
		stamped_printf("\tParsing error : Invalid string value for 'Velocity_Type' : %s\n", str);
		Exit(0);
	}
  
  // 外部雰囲気温
  if ( HeatProblem ) {
    if ( elmL->GetValue(CfgIdt("ambient_temperature"), &ct) ) {
      BaseBc[n].set_Temp( FBUtility::convTemp2K(ct, Unit_Temp) );
    }
    else {
      stamped_printf("\tParsing error : fail to get 'ambient_temperature' in 'Basic_BCs > IN_OUT'\n");
      Exit(0);
    }
    if ( Unit_Param != DIMENSIONAL ) {
      stamped_printf("\tError: Heat condition must be given by dimensional value\n");
      Exit(0);
    }
  }
}

/**
 @fn void ParseBC::getXML_OBC_Outflow(const CfgElem *elmL, unsigned n)
 @brief 外部境界の流出条件のパラメータを取得する
 @param elmL 
 @param n 面番号
 @note 圧力の値は，Control::setParameters()で無次元化する
 */
void ParseBC::getXML_OBC_Outflow(const CfgElem *elmL, unsigned n)
{
	const char *str=NULL;
  REAL_TYPE ct;
  
  // 圧力境界のタイプ  default
  BaseBc[n].set_pType(P_GRAD_ZERO);
  BaseBc[n].p = 0.0; // ダミー値
  
  /* Hidden option
   if ( !elmL->GetValue(CfgIdt("pressure_type"), &str) ) {
   ;
   }
   else {
   if ( !strcasecmp("dirichlet", str) ) {
   BaseBc[n].set_pType(P_DIRICHLET);
   }
   else if ( !strcasecmp("grad_zero", str) ) {
   BaseBc[n].set_pType(P_GRAD_ZERO);
   }
   else {
   stamped_printf("\tParsing error : Invalid string value for 'Pressure_Type' : %s\n", str);
   Exit(0);
   }
   }
   */
  
  // 流出速度のタイプ
  if ( !elmL->GetValue(CfgIdt("velocity_type"), &str) ) {
    stamped_printf("\tParsing error : fail to get 'Velocity_Type' in 'Basic_BCs > outflow'\n");
    Exit(0);
  }
  if ( !strcasecmp("average", str) ) {
		BaseBc[n].set_ofv(V_AVERAGE);
	}
	else if ( !strcasecmp("minmax", str) ) {
		BaseBc[n].set_ofv(V_MINMAX);
	}
	else {
		stamped_printf("\tParsing error : Invalid string value for 'Velocity_Type' : %s\n", str);
		Exit(0);
	}
  
  // 圧力の値
  if ( BaseBc[n].get_pType() == P_DIRICHLET ) {
    if ( !strcasecmp("dirichlet", str) ) {
      if ( elmL->GetValue(CfgIdt("pressure_value"), &ct) ) {
        BaseBc[n].p = ct;
      }
      else {
        stamped_printf("\tParsing error : fail to get 'pressure_value' in 'Basic_BCs > outflow'\n");
        Exit(0);
      }
    }
  }
}

/**
 @fn void ParseBC::getXML_OBC_Periodic(const CfgElem *elmL, unsigned n)
 @brief 外部境界の周期条件のパラメータを取得する
 @param elmL 
 @param n 面番号
 @note 圧力の値は，Control::setParameters()で無次元化する
 */
void ParseBC::getXML_OBC_Periodic(const CfgElem *elmL, unsigned n)
{
  REAL_TYPE ct;
  int def;
  const char* str=NULL;
  
  // モード
  if ( elmL->GetValue(CfgIdt("mode"), &str) ) {
    if ( !strcasecmp(str, "simple_copy") ) {
      BaseBc[n].set_PrdcMode(BoundaryOuter::prdc_Simple);
    }
    else if ( !strcasecmp(str, "directional") ) {
      BaseBc[n].set_PrdcMode(BoundaryOuter::prdc_Directional);
    }
    else if ( !strcasecmp(str, "driver") ) {
      BaseBc[n].set_PrdcMode(BoundaryOuter::prdc_Driver);
    }
  }
  else {
    printf("\tParsing error : No 'mode' section in 'Basic_BCs > periodic'\n");
    Exit(0);
  }
  
  // Directional
  if ( BaseBc[n].get_PrdcMode() == BoundaryOuter::prdc_Directional ) {
    if ( elmL->GetValue(CfgIdt("pressure_difference"), &ct)) {
      BaseBc[n].p = ct;
    }
    else {
      printf("\tParsing error : No 'Pressure_Difference' keyword in 'Basic_BCs > Periodic'\n");
      Exit(0);
    }
    if ( elmL->GetValue(CfgIdt("flow_direction"), &str)) {
      if ( !strcasecmp(str, "upstream") ) {
        BaseBc[n].set_FaceMode(BoundaryOuter::prdc_upstream);
      }
      else if ( !strcasecmp(str, "downstream") ) {
        BaseBc[n].set_FaceMode(BoundaryOuter::prdc_downstream);
      }
      else {
        printf("\tParsing error : Invalid keyword in 'Basic_BCs > Periodic > Flow_Direction'\n");
        Exit(0);
      }
    }
    else {
      printf("\tParsing error : No 'Flow_Direction' keyword in 'Basic_BCs > Periodic'\n");
      Exit(0);
    }
  }
  
  // Driver
  if ( BaseBc[n].get_PrdcMode() == BoundaryOuter::prdc_Driver ) {
    if ( elmL->GetValue(CfgIdt("driver_direction"), &str)) {
      if ( !strcasecmp("x_minus", str) ) {
        def = X_MINUS;
      }
      else if ( !strcasecmp("x_plus", str) ) {
        def = X_PLUS;
      }
      else if ( !strcasecmp("y_minus", str) ) {
        def = Y_MINUS;
      }
      else if ( !strcasecmp("y_plus", str) ) {
        def = Y_PLUS;
      }
      else if ( !strcasecmp("z_minus", str) ) {
        def = Z_MINUS;
      }
      else if ( !strcasecmp("z_plus", str) ) {
        def = Z_PLUS;
      }
      else {
        printf("\tParsing error : Invalid keyword in 'Basic_BCs > Periodic > Driver_Direction'\n");
        Exit(0);
      }
      BaseBc[n].set_DriverDir(def);
    }
    else {
      printf("\tParsing error : No 'Driver_Direction' keyword in 'Basic_BCs > Periodic'\n");
      Exit(0);
    }
    
    if ( elmL->GetValue(CfgIdt("driver_lid_index"), &def)) {
      BaseBc[n].set_DriverIndex(def);
    }
    else {
      printf("\tParsing error : Invalid 'Driver_Lid_Index' keyword in 'Basic_BCs > Periodic'\n");
      Exit(0);
    }
    
  }
}

/**
 @fn void ParseBC::getXML_OBC_SpecVH(const CfgElem *elmL, unsigned n)
 @brief 外部境界の流入条件のパラメータを取得する
 @param elmL 
 @param n 面番号
 */
void ParseBC::getXML_OBC_SpecVH(const CfgElem *elmL, unsigned n)
{
  REAL_TYPE vel, ct;
  REAL_TYPE v[3];
  const char* str=NULL;
  
  // 速度境界条件のタイプ
  BaseBc[n].set_vType( getXML_Vel_profile(elmL, "Basic_BCs > Specified_Velocity") );
  
  // 法線ベクトル
  get_NV(elmL, n, "Basic_BCs > Specified_Velocity", v);
  BaseBc[n].addVec(v);
  
  // 速度の指定モードの特定
  if ( !elmL->GetValue(CfgIdt("specified_type"), &str) ) {
    stamped_printf("\tParsing error : Invalid Specified_Type in 'Basic_BCs > Specified_Type'\n");
    Exit(0);
  }
  if ( !strcasecmp("velocity", str) ) {
		compo[n].set_VBC_policy(true);
	}
	//else if ( !strcasecmp("massflow", str) ) {
  //compo[n].set_VBC_policy(false);
	//}
	else {
		printf("\tParsing error : Invalid string value '%s' for 'Specified_Type'\n", str);
		Exit(0);
	}
	
  // 指定値の取得
  if ( !elmL->GetValue(CfgIdt("specified_value"), &ct) ) {
    stamped_printf("\tParsing error : Invalid value in 'Basic_BCs > Specified_Value'\n");
    Exit(0);
  }
  vel = ( Unit_Param == DIMENSIONAL ) ? ct : ct * RefVelocity; // 有次元値で保持
  
  // 速度のパラメータ読み込み
  getXML_Vel_Params(elmL, BaseBc[n].get_vType(), BaseBc[n].ca, vel, "Basic_BCs > Specified_Value");
  
  // heat problem
  if ( HeatProblem ) {
    if ( Unit_Param != DIMENSIONAL ) {
      stamped_printf("\tError: Heat condition must be given by dimensional value\n");
      Exit(0);
    }
    
    if ( !elmL->GetValue(CfgIdt("temperature"), &ct) ) {
      stamped_printf("\tParsing error : fail to get 'Temperature' in 'Basic_BCs > Specified_Velocity'\n");
      Exit(0);
    }
    else {
      BaseBc[n].set_Temp( FBUtility::convTemp2K(ct, Unit_Temp) );
      BaseBc[n].set_hType(CNST_TEMP);
    }
  }
}

/**
 @fn void ParseBC::getXML_OBC_Trcfree(const CfgElem *elmL, unsigned n)
 @brief 外部境界の流入条件のパラメータを取得する
 @param elmL 
 @param n 面番号
 */
void ParseBC::getXML_OBC_Trcfree(const CfgElem *elmL, unsigned n)
{
  REAL_TYPE ct;
  
  BaseBc[n].set_pType(P_DIRICHLET);
  BaseBc[n].p = 0.0; // ゲージ圧zero 固定
  
  // 外部雰囲気温
  if ( HeatProblem ) {
    if ( elmL->GetValue(CfgIdt("ambient_temperature"), &ct) ) {
      BaseBc[n].set_Temp( FBUtility::convTemp2K(ct, Unit_Temp) );
    }
    else {
      stamped_printf("\tParsing error : fail to get 'ambient_temperature' in 'Basic_BCs > IN_OUT'\n");
      Exit(0);
    }
    if ( Unit_Param != DIMENSIONAL ) {
      stamped_printf("\tError: Heat condition must be given by dimensional value\n");
      Exit(0);
    }
  }
}

/**
 @fn void ParseBC::getXML_OBC_Wall(const CfgElem *elmL, unsigned n)
 @brief 外部境界の壁面条件のパラメータを取得する
 @param elmL 
 @param n 面番号
 */
void ParseBC::getXML_OBC_Wall(const CfgElem *elmL, unsigned n)
{
  REAL_TYPE vel, ct;
  REAL_TYPE v[3];
  const char* str=NULL;
  
  // 速度境界条件のタイプ
  BaseBc[n].set_vType( getXML_Vel_profile(elmL, "Basic_BCs > Wall") );
  
  // 法線ベクトル
  get_NV(elmL, n, "Basic_BCs > Wall", v);
  BaseBc[n].addVec(v);
  
  // 速度の指定モードの特定
  if ( !elmL->GetValue(CfgIdt("specified_type"), &str) ) {
    stamped_printf("\tParsing error : Invalid Specified_Type in 'Basic_BCs > Specified_Type'\n");
    Exit(0);
  }
  if ( !strcasecmp("velocity", str) ) {
		compo[n].set_VBC_policy(true);
	}
	//else if ( !strcasecmp("massflow", str) ) {
  //compo[n].set_VBC_policy(false);
	//}
	else {
		printf("\tParsing error : Invalid string value '%s' for 'Specified_Type'\n", str);
		Exit(0);
	}
  
  // 指定値の取得
  if ( !elmL->GetValue(CfgIdt("specified_value"), &ct) ) {
    stamped_printf("\tParsing error : Invalid value in 'Basic_BCs > Specified_Value'\n");
    Exit(0);
  }
  vel = ( Unit_Param == DIMENSIONAL ) ? ct : ct * RefVelocity; // 有次元値で保持
  
  // 速度のパラメータ読み込み
  getXML_Vel_Params(elmL, BaseBc[n].get_vType(), BaseBc[n].ca, vel, "Basic_BCs > Wall");
  
  // heat problem
  if ( HeatProblem ) {
    
    char *str2=NULL;
    
    if ( !elmL->GetValue(CfgIdt("heat_type"), &str) ) {
      stamped_printf("\tParsing error : fail to get 'Heat_Type' in 'Basic_BCs > wall'\n");
      Exit(0);
    }
    if ( Unit_Param != DIMENSIONAL ) {
      stamped_printf("\tError: Heat condition must be given by dimensional value\n");
      Exit(0);
    }
    
    if     ( !strcasecmp(str, "Adiabatic") )    {
      BaseBc[n].set_hType(ADIABATIC);
      BaseBc[n].p = 0.0;
    }
    else if( !strcasecmp(strncpy(str2,str,12), "HeatTransfer") ) {
      BaseBc[n].set_hType(TRANSFER);
      getXML_OBC_HT(elmL, n, str);
    }
    else if( !strcasecmp(str, "HeatFlux") )     {
      BaseBc[n].set_hType(HEATFLUX);
      BaseBc[n].set_Heatflux( get_BCval_real(elmL, "Heat_Flux") ); // 正符号は流入
    }
    else if( !strcasecmp(str, "Isothermal") )   {
      BaseBc[n].set_hType(ISOTHERMAL);
      ct = get_BCval_real(elmL, "temperature"); // 表面温度
      BaseBc[n].set_Temp( FBUtility::convTemp2K(ct, Unit_Temp) );
    }
    /*
     else if( !strcasecmp(str, "Constant_Temperature") )   {
     BaseBc[n].set_hType(CNST_TEMP);
     ct = get_BCval_real(elmL, "temperature"); // 指定温度
     BaseBc[n].T1.temp = FBUtility::convTemp2K(ct, Unit_Temp);
     }
     */
    else {
      stamped_printf("\tParsing error : Invalid string value for 'Heat_type' : %s\n", str);
      Exit(0);
    }    
  }
}

/**
 @fn void ParseBC::getXML_Phase(void)
 @brief 2相流問題で気相か液相かを取得する
 */
void ParseBC::getXML_Phase(void)
{
	const CfgElem *elemTop=NULL, *elmL1=NULL;
  const CfgParam* param=NULL;
	const char *str=NULL;
  const char* p=NULL;
  unsigned m_phase;
  int id;
	
	elemTop = CF->GetTop(STEER);
	if( !(elmL1 = elemTop->GetElemFirst("Phase_Idetification")) ) {
    stamped_printf("\tParsing error : Missing the section of 'Phase_Idetification'\n");
    Exit(0);
  }
	
  // load statement list
  param = elmL1->GetParamFirst();
  int NoParam = elmL1->GetParamSize();
  for (int i=1; i<=NoParam; i++) {
    
    // state
    p = param->GetName();
    if ( !strcasecmp(p, "gas") ) {
      m_phase = GAS;
    }
    else if ( !strcasecmp(p, "liquid") ) {
      m_phase = LIQUID;
    }
    else {
      stamped_printf("\tParsing error : No valid keyword '%s' in 'Phase_Idetification'\n", p);
      Exit(0);
    }
    
    // ID
    if ( !param->isSetID() ) {
      stamped_printf("\tParsing error : No ID for statement in 'Phase_Idetification'\n");
      Exit(0);
    }
    if ( -1 == (id=param->GetID()) ) {
      stamped_printf("\tParsing error : No valid ID for statement in 'Phase_Idetification'\n");
      Exit(0);
    }
    
    // IDがiTableの中にリストアップされているかを調べる
    if ( !isIDinTable(id) ) {
      stamped_printf("\tParsing error : ID[%d] described in 'id' is not included in 'Model_Setting'\n", id);
      Exit(0);
    }
    
    // set phase of FLUID
    for (int n=1; n<=NoID; n++) {
      if ( compo[n].getID() == id ) {
        if ( compo[n].getState() == FLUID ) compo[n].setPhase(m_phase);
      }
    }
  }
  
  // check Phase of Fluid
  unsigned tmp;
  bool sw=true;
  for (int n=1; n<=NoID; n++) {
    if ( compo[n].getState() == FLUID ) {
      tmp = compo[n].getPhase();
      if ( (tmp!=GAS) && (tmp!=LIQUID) ) {
        stamped_printf("\tcomponent [%d] is fluid, but not identified by gas or liquid.\n", n);
        sw = false;
      }
    }
  }
  if ( sw == false ) Exit(0);
}


/**
 @fn void ParseBC::getXML_Vel_Params(const CfgElem *elmL, unsigned prof, REAL_TYPE* ca, REAL_TYPE vel, const char* err_str)
 @brief 速度のパラメータを取得する
 @param elmL 
 @param prof 速度プロファイル
 @param ca 係数パラメータの配列
 @param vel 設定する速度（有次元値）
 @param err_str エラー表示用文字列
 @note 
 - 値は，Control::setParameters()で無次元化する
 - 速度プロファイルは単振動と一定値の場合で係数の保持パターンが異なる
 - 内部境界の場合には，流量指定と速度指定があるので分岐処理（未実装）
 */
void ParseBC::getXML_Vel_Params(const CfgElem *elmL, unsigned prof, REAL_TYPE* ca, REAL_TYPE vel, const char* err_str)
{
  REAL_TYPE ct=0.0;
  
  if ( prof == CompoList::vel_harmonic) {
    ca[CompoList::amplitude] = vel;
    
    if ( !elmL->GetValue(CfgIdt("frequency"), &ct) ) {
      stamped_printf("\tParsing error : fail to get 'Frequency' in '%s'\n", err_str);
      Exit(0);
    }
    ca[CompoList::frequency] = ct;
    
    if ( !elmL->GetValue(CfgIdt("initial_phase"), &ct) ) {
      stamped_printf("\tParsing error : fail to get 'Initial_phase' in '%s'\n", err_str);
      Exit(0);
    }
    ca[CompoList::initphase] = ct;
    
    if ( !elmL->GetValue(CfgIdt("constant_bias"), &ct) ) {
      stamped_printf("\tParsing error : fail to get 'Constant_Bias' in '%s'\n", err_str);
      Exit(0);
    }
    ca[CompoList::bias] = ct;
  }
  else if ( prof == CompoList::vel_constant ) {
    ca[CompoList::amplitude] = 0.0;
    ca[CompoList::frequency] = 0.0;
    ca[CompoList::initphase] = 0.0;
    ca[CompoList::bias]      = vel;
  }
  else { // vel_zero
    ca[CompoList::amplitude] = 0.0;
    ca[CompoList::frequency] = 0.0;
    ca[CompoList::initphase] = 0.0;
    ca[CompoList::bias]      = 0.0;
  }
  
}


/**
 @fn unsigned ParseBC::getXML_Vel_profile(const CfgElem *elmL, const char* err_str)
 @brief 外部境界の速度境界条件のタイプを取得し，返す
 @param elmL 
 @param err_str 表示用ストリング
 */
unsigned ParseBC::getXML_Vel_profile(const CfgElem *elmL, const char* err_str)
{
	const char *str=NULL;
  
  if ( !elmL->GetValue(CfgIdt("Profile"), &str) ) {
    printf("\tParsing error : fail to get 'Profile' in '%s'\n", err_str);
    Exit(0);
  }
	if ( !strcasecmp("constant", str) ) {
		return CompoList::vel_constant;
	}
	else if ( !strcasecmp("harmonic", str) ) {
		return CompoList::vel_harmonic;
	}
  else if ( !strcasecmp("zero", str) ) {
		return CompoList::vel_zero;
	}
	else {
		printf("\tParsing error : Invalid string value for '%s > Profile' : %s\n", str, err_str);
		Exit(0);
	}
  return 0;
}


/**
 @fn bool ParseBC::isComponent(unsigned label)
 @brief コンポーネントが存在するかどうかを調べる
 @retval bool値
 */
bool ParseBC::isComponent(unsigned label)
{
  for (int n=1; n<=NoBC; n++) {
    if ( compo[n].getType() == label ) return true;
  }
  return false;
}

/**
 @fn bool ParseBC::isCompoTransfer(unsigned label)
 @brief HTコンポーネントが存在するかどうかを調べる
 @retval bool値
 */
bool ParseBC::isCompoTransfer(unsigned label)
{
  for (int n=1; n<=NoBC; n++) {
    if ( compo[n].getHtype() == label ) return true;
  }
  return false;
}

/**
 @fn bool ParseBC::isIDinCompo(int candidate_id, unsigned now)
 @brief candidate_idが既にコンポーネントに登録されているかを調べる
 @retval 重複していればfalseを返す
 @param candidate_id チェックの候補ID
 @param now コンポーネントリストの現在までのエントリ番号
 */
bool ParseBC::isIDinCompo(unsigned candidate_id, unsigned now)
{
  for (unsigned i=1; i<now; i++) {
    if ( candidate_id == compo[i].getID() ) {
      return false;
    }
  }
  return true;
}

/**
 @fn bool ParseBC::isIDinCompo(int candidate_id, int def, unsigned now)
 @brief candidate_idとdefの組が既にコンポーネントに登録されているかを調べる
 @retval 重複していればfalseを返す
 @param candidate_id チェックの候補ID
 @param def deface
 @param now コンポーネントリストの現在までのエントリ番号
 */
bool ParseBC::isIDinCompo(unsigned candidate_id, int def, unsigned now)
{
  for (unsigned i=1; i<now; i++) {
    if ( (candidate_id == compo[i].getID()) && (def == compo[i].getDef()) ) {
      return false;
    }
  }
  return true;
}

/**
 @fn bool ParseBC::isIDinTable(int candidate)
 @brief candidateがiTableの中にリストアップされているかを調べる
 @retval リストにあればtrue
 */
bool ParseBC::isIDinTable(int candidate)
{
  for (int n=1; n<=NoID; n++) {
    if ( iTable[n].getID() == (unsigned)candidate ) return true;
  }
  return false;
}

/**
 @fn void ParseBC::loadOuterBC(void)
 @brief XMLファイルをパースして，外部境界条件を取得，保持する
 @note
 - BasicBCsをパースする
 - idの重複をチェック
 - 各外部境界面に対して，FaceBCを設定する
 */
void ParseBC::loadOuterBC(void)
{
  const CfgElem *elemTop, *elmL1, *elmL2;
  elemTop = elmL1 = elmL2 = NULL;
  int id=0, cid=0;
  const char *Ename=NULL;
  
  // Basic Outer BCリストの読み込み
  elemTop = CF->GetTop(OUTERBND);
  elmL1 = elemTop->GetElemFirst("Basic_BCs");
  elmL2 = elmL1->GetElemFirst();
  if ( !elmL2 ) Exit(0);
  
  for (unsigned i=0; i<NoBaseBC; i++) {
    if ( !(Ename = elmL2->GetName()) ) {
      printf("\tParsing error : No Elem name in 'Basic_BCs'\n");
      Exit(0);
    }
    chkKeywordOBC(Ename, i); // BCTypeに境界条件の種別をセットする
    
    // 境界条件番号
    if ( !elmL2->isSetID() ) {
      printf("\tParsing error : No ID section in Basic_BCs\n");
      Exit(0);
    }
    if ( -1 == (id=elmL2->GetID()) ) {
      printf("\tParsing error : No valid ID for Basic_BCs\n");
      Exit(0);
    }
    
    BaseBc[i].set_BC_ID(id);
    
    switch ( BaseBc[i].get_BCtype() ) {
      case OBC_WALL:
        getXML_OBC_Wall(elmL2, i);
        break;
        
      case OBC_OUTFLOW:
        getXML_OBC_Outflow(elmL2, i);
        break;
        
      case OBC_SPEC_VEL:
        getXML_OBC_SpecVH(elmL2, i);
        break;
        
      case OBC_TRC_FREE:
        getXML_OBC_Trcfree(elmL2, i);
        break;
        
      case OBC_FAR_FIELD:
        getXML_OBC_FarField(elmL2, i);
        break;
        
      case OBC_PERIODIC:
        getXML_OBC_Periodic(elmL2, i);
        break;
        
      case OBC_SYMMETRIC:
        // nothing to do
        break;
    }
    elmL2 = elmL1->GetElemNext(elmL2);
  }
  
  // IDの重複をチェックする
  if ( !chkID() ) Exit(0);
  
  // 各フェイスに境界条件設定する
  if ( !(elmL1 = elemTop->GetElemFirst("Face_BC")) ) {
    stamped_printf("\tParsing error : Missing 'Face_BC' description\n");
    Exit(0);
  }
  
  const CfgParam* param=NULL;
  
  // 各面に与える境界条件番号を取得し，BasicListから境界情報リストに内容をコピーする．ただし，ガイドセルのセルIDと媒質番号は後で設定
  for (int face=0; face<NOFACE; face++) {
    // faceに対するエントリを得る
    if ( !(elmL2 = selectFace(face, elmL1)) ) Exit(0);
    
    // 境界条件番号を取得
    if ( !elmL2->isSetID() ) {
      printf("\tParsing error : No ID section in Basic_BCs\n");
      Exit(0);
    }
    if ( -1 == (id=elmL2->GetID()) ) {
      printf("\tParsing error : No valid ID for Basic_BCs\n");
      Exit(0);
    }
    
    // BaseBC[]からbc[]へ内容のコピー
    for (unsigned i=0; i<NoBaseBC; i++) {
      if ( BaseBc[i].get_BC_ID() == id ) bc[face].dataCopy( &BaseBc[i] );
    }
  }
  
  // ガイドセルのセルIDと媒質番号を設定する
  for (int face=0; face<NOFACE; face++) {
    
    if ( !(elmL2 = selectFace(face, elmL1)) ) Exit(0);
    
    // セルIDの取得
    if ( !(param = elmL2->GetParamFirst("Guide_Cell_ID")) ) {
      stamped_printf("\tParsing error : No entory 'Guide_Cell_ID' in 'Face_BC'\n");
      Exit(0);
    }
    else {
      if ( !param->isSetID() ) {
        stamped_printf("\tParsing error : No ID section in 'Guide_Cell_ID'\n");
        Exit(0);
      }
      if ( 1 > (cid = param->GetID()) ) {
        stamped_printf("\tParsing error : Invalid Outer Guide Cell ID[%d] that shoud be > 0\n", cid);
        Exit(0);
      }
      
      bc[face].set_GuideID(cid);
    }
    
    // セルIDから媒質番号を求める
    for (int i=1; i<=NoID; i++) {
      if ( iTable[i].getID() == (unsigned)bc[face].get_GuideID() ) {
        bc[face].set_GuideMedium( (int)iTable[i].getMatID() );
      }
    }
  }
  
  // チェック
  for (int face=0; face<NOFACE; face++) {
    if ( (id=bc[face].get_BC_ID()) == 0 ) {
      printf("\tFace BC : id=%d is not listed in 'Basic_BCs' section\n", id);
      Exit(0);
    }
  }
  
  // 周期境界条件の整合性のチェック
  
  // 部分周期境界の数
  unsigned p_flag=0;
  for (unsigned n=0; n<NOFACE; n++) {
    if (bc[n].get_PrdcMode() == BoundaryOuter::prdc_Driver) p_flag++;
  }
  
  // 部分周期条件を使わない場合，対になる外部境界のチェック
  if ( p_flag == 0 ) {
    unsigned n_pair=0;
    
    // 周期境界条件の指定をチェック
    for (unsigned n=0; n<NOFACE; n++) {
      if ( bc[n].get_BCtype() == OBC_PERIODIC ) {
        n_pair = oppositDir(n);
        if ( bc[n_pair].get_BCtype() != OBC_PERIODIC ) {
          printf("\tFace BC : No consistent Periodic Bnoudary in %s direction\n", FBUtility::getDirection(n_pair).c_str());
          Exit(0);
        }
      }
    }
    
    // 対になるモードのチェック
    for (unsigned n=0; n<NOFACE; n++) {
      if ( bc[n].get_BCtype() == OBC_PERIODIC ) {
        n_pair = oppositDir(n);
        
        switch (bc[n].get_PrdcMode()) {
          case BoundaryOuter::prdc_Simple:
            if ( bc[n_pair].get_PrdcMode() != BoundaryOuter::prdc_Simple ) { 
              printf("\tFace BC : No consistent SIMPLE Periodic Bnoudary in %s direction\n", FBUtility::getDirection(n_pair).c_str());
              Exit(0);
            }
            break;
            
          case BoundaryOuter::prdc_Directional:
            if ( bc[n_pair].get_PrdcMode() != BoundaryOuter::prdc_Directional ) {
              printf("\tFace BC : No consistent DIRECTIONAL Periodic Bnoudary in %s direction\n", FBUtility::getDirection(n_pair).c_str());
              Exit(0);
            }
            if ( bc[n].p != bc[n_pair].p ) { // 同じ値が入っていること
              printf("\tFace BC : Pressure difference value is not same in %s direction\n", FBUtility::getDirection(n_pair).c_str());
              Exit(0);
            }
            if ( (bc[n].get_FaceMode() == BoundaryOuter::prdc_upstream) && (bc[n_pair].get_FaceMode() != BoundaryOuter::prdc_downstream) ) {
              printf("\tFace BC : No consistent Upstream/Downstream relation in %s direction\n", FBUtility::getDirection(n).c_str());
              Exit(0);
            }
            if ( (bc[n].get_FaceMode() == BoundaryOuter::prdc_downstream) && (bc[n_pair].get_FaceMode() != BoundaryOuter::prdc_upstream) ) {
              printf("\tFace BC : No consistent Upstream/Downstream relation in %s direction\n", FBUtility::getDirection(n).c_str());
              Exit(0);
            }
            
            break;
        }        
      }
    }
  }
  else { // Driverが指定された場合の内部境界との整合性をチェック
    unsigned n_pair=0;
    for (unsigned n=0; n<NOFACE; n++) {
      n_pair = oppositDir(n);
      if ( bc[n].get_BCtype() == OBC_PERIODIC ) {
        if (bc[n].get_PrdcMode() == BoundaryOuter::prdc_Driver) {
          
          // 他方は周期境界以外であること
          if ( bc[n_pair].get_BCtype() == OBC_PERIODIC ) {
            printf("\tFace BC : %s direction should be non periodic BC\n", FBUtility::getDirection(n_pair).c_str());
            Exit(0);
          }
          
          unsigned cflag=0;
          for (unsigned c=1; c<=NoBC; c++) {
            if ( compo[c].getType() == PERIODIC ) {
              if ( (int)compo[c].getPeriodicDir() != bc[n].get_DriverDir() ) {
                printf("\tPeriodic Driver BC : No consistent Periodic Bnoudary in %s direction\n", FBUtility::getDirection(n_pair).c_str());
                Exit(0);
              }
              else {
                cflag++;
              }
            }
          }
          if (cflag != 1) {
            printf("\tPeriodic Driver BC can not detemine uniquely\n");
            Exit(0);
          }
        }        
      }
    }
  }
  
}


/**
 @fn unsigned ParseBC::oppositDir(unsigned dir)
 @brief 外部境界面の反対方向を返す
 @param dir 評価する方向
 */
unsigned ParseBC::oppositDir(unsigned dir)
{
  unsigned n_pair=0;
  
  if      ( dir == X_MINUS ) n_pair=X_PLUS;
  else if ( dir == X_PLUS )  n_pair=X_MINUS;
  else if ( dir == Y_MINUS ) n_pair=Y_PLUS;
  else if ( dir == Y_PLUS )  n_pair=Y_MINUS;
  else if ( dir == Z_MINUS ) n_pair=Z_PLUS;
  else if ( dir == Z_PLUS )  n_pair=Z_MINUS;
  
  return n_pair;
}


//@fn void ParseBC::printBaseOBC(FILE* fp)
//@brief 基本境界条件リストを表示する
//@param fp
void ParseBC::printBaseOBC(FILE* fp)
{
  fprintf(fp, "\n\nDEBUG : Basic Boundary Conditions\n");
  fprintf(fp, "\t\t   #   Variable     id  \n");
  for (unsigned i=0; i<NoBaseBC; i++) {
    fprintf(fp, "\t\t%4d   %3d\n", i, BaseBc[i].get_BC_ID());
  }
  fflush(fp);
}


/**
 @fn void ParseBC::printCompo(FILE* fp, REAL_TYPE* nv, int* gci, MaterialList* mat)
 @brief コンポーネントの情報を表示する
 @param nv ボクセルモデルから計算した法線
 @param gci グローバルなコンポーネントのインデクス
 @param mat MaterialList
 */
void ParseBC::printCompo(FILE* fp, REAL_TYPE* nv, int* gci, MaterialList* mat)
{
  unsigned n, m;
  bool flag;
  
  // VBC ---------------------------------------------------
  if ( isComponent(SPEC_VEL) || isComponent(SPEC_VEL_WH) ) {
    fprintf(fp, "\n\t[SPECIFIED_VELOCITY]\n");
    fprintf(fp, "\t no                    Label    ID   def    i_st    i_ed    j_st    j_ed    k_st    k_ed   Area[m*m]\n");
    
    for(n=1; n<=NoBC; n++){
      if ( (compo[n].getType() == SPEC_VEL) || (compo[n].getType() == SPEC_VEL_WH) )  {
        fprintf(fp,"\t%3d %24s %5d %5d %7d %7d %7d %7d %7d %7d %11.4e\n",
                n, compo[n].name, compo[n].getID(), compo[n].getDef(),
                getCmpGbbox_st_x(n, gci), getCmpGbbox_ed_x(n, gci), 
                getCmpGbbox_st_y(n, gci), getCmpGbbox_ed_y(n, gci), 
                getCmpGbbox_st_z(n, gci), getCmpGbbox_ed_z(n, gci), 
                compo[n].area);
      }
    }
    fprintf(fp, "\n");
    
    fprintf(fp, "\t no                    Label    ID   normal_x   normal_y   normal_z      Direction");
    
    if ( compo[n].get_sw_V_profile() == CompoList::vel_zero ) {
      fprintf(fp, "\n");
    }
    else if ( compo[n].get_sw_V_profile() == CompoList::vel_constant ) {
      fprintf(fp, "    Q[m^3/s]   Vel.[m/s]\n");
    }
    else {
      fprintf(fp, "    Q[m^3/s]   Amp.[m/s]   Freq.[Hz]  Phase[rad] Intcpt[m/s] Strauhal[-]\n");
    }
    
    for(n=1; n<=NoBC; n++){
      if ( (compo[n].getType() == SPEC_VEL) || (compo[n].getType() == SPEC_VEL_WH) )  {
        fprintf(fp,"\t%3d %24s %5d %10.3e %10.3e %10.3e %14s ",
                n, compo[n].name, compo[n].getID(), compo[n].nv[0], compo[n].nv[1], compo[n].nv[2],
                (compo[n].getBClocation()==CompoList::same_direction) ? "same dir." : "opposite dir.");
        
        if ( compo[n].get_sw_V_profile() == CompoList::vel_zero ) {
          fprintf(fp,"\n");
        }
        else if ( compo[n].get_sw_V_profile() == CompoList::vel_constant ) {
          fprintf(fp,"%11.3e %11.3e\n", compo[n].ca[CompoList::bias]*compo[n].area, compo[n].ca[CompoList::bias] );
        }
        else {
          fprintf(fp,"%11.3e %11.3e %11.3e %11.3e %11.3e %11.4e\n",
                  compo[n].ca[CompoList::amplitude]*compo[n].area,
                  compo[n].ca[CompoList::amplitude],
                  compo[n].ca[CompoList::frequency],
                  compo[n].ca[CompoList::initphase],
                  compo[n].ca[CompoList::bias],
                  compo[n].ca[CompoList::frequency]*RefLength/RefVelocity);
        }
        
      }
    }
    
    // with constant temperature
    if ( isComponent(SPEC_VEL_WH) ) {
      fprintf(fp, "\n\t[SPECIFIED_VELOCITY with constant temperature]\n");
      fprintf(fp, "\t no                    Label    ID   def    i_st    i_ed    j_st    j_ed    k_st    k_ed      Temp(%s)      Temp[-]\n", (Unit_Temp==Unit_KELVIN) ? "K" : "C");
      
      for(n=1; n<=NoBC; n++) {
        if ( compo[n].getType() == SPEC_VEL_WH )  {
          fprintf(fp, "\t%3d %24s %5d %5d %7d %7d %7d %7d %7d %7d %12.4e %12.4e\n", 
                  n, compo[n].name, compo[n].getID(), compo[n].getDef(), 
                  getCmpGbbox_st_x(n, gci), getCmpGbbox_ed_x(n, gci), 
                  getCmpGbbox_st_y(n, gci), getCmpGbbox_ed_y(n, gci), 
                  getCmpGbbox_st_z(n, gci), getCmpGbbox_ed_z(n, gci),  
									FBUtility::convK2Temp(compo[n].get_Temp(), Unit_Temp), 
									FBUtility::convK2ND(compo[n].get_Temp(), BaseTemp, DiffTemp)); // 保持されている温度はKelvin
        }
      }
      fprintf(fp, "\n");
    }
  }
  
  if ( isComponent(OUTFLOW) ) {
    fprintf(fp, "\n\t[OUTFLOW]\n");
    fprintf(fp, "\t no                    Label    ID   def    i_st    i_ed    j_st    j_ed    k_st    k_ed  outflow_vel  pressure\n");
    
    for(n=1; n<=NoBC; n++){
      if ( compo[n].getType() == OUTFLOW )  {
        fprintf(fp,"\t%3d %24s %5d %5d %7d %7d %7d %7d %7d %7d ",
                n, compo[n].name, compo[n].getID(), compo[n].getDef(),
                getCmpGbbox_st_x(n, gci), getCmpGbbox_ed_x(n, gci), 
                getCmpGbbox_st_y(n, gci), getCmpGbbox_ed_y(n, gci), 
                getCmpGbbox_st_z(n, gci), getCmpGbbox_ed_z(n, gci), 
                (compo[n].getOutflowType() == V_AVERAGE) ? "Average " : "Minmax ");
        if (compo[n].get_sw_P_BCtype() == P_DIRICHLET) {
          fprintf(fp, "%12.4e;", FBUtility::convND2D_P(compo[n].get_Pressure(), BasePrs, RefDensity, RefVelocity, Unit_Prs) );
        }
        else {
          fprintf(fp," Grad_p = 0   ---");
        }
        fprintf(fp, "\n");
      }
    }
    fprintf(fp, "\n");
  }
  
  // Forcing
  if ( isComponent(IBM_DF) ) {
    fprintf(fp, "\n\t[IBM_DF]\n");
    fprintf(fp, "\t no                    Label    ID    i_st    i_ed    j_st    j_ed    k_st    k_ed\n");
    
    for(n=1; n<=NoBC; n++){
      if ( compo[n].getType() == IBM_DF )  {
        fprintf(fp,"\t%3d %24s %5d %7d %7d %7d %7d %7d %7d\n",
                n, compo[n].name, compo[n].getID(), 
                getCmpGbbox_st_x(n, gci), getCmpGbbox_ed_x(n, gci), 
                getCmpGbbox_st_y(n, gci), getCmpGbbox_ed_y(n, gci), 
                getCmpGbbox_st_z(n, gci), getCmpGbbox_ed_z(n, gci));
      }
    }
    fprintf(fp, "\n");
    
    fprintf(fp, "\t no                    Label    ID   normal_x   normal_y   normal_z    Vel.[m/s]      Vel.[-]\n");
    
    for(n=1; n<=NoBC; n++){
      if ( compo[n].getType() == IBM_DF )  {
        fprintf(fp,"\t%3d %24s %5d %10.3e %10.3e %10.3e %12.4e %12.4e \n",
                n, compo[n].name, compo[n].getID(), compo[n].nv[0], compo[n].nv[1], compo[n].nv[2],
                compo[n].get_Velocity(), FBUtility::convD2ND_V(compo[n].get_Velocity(), RefVelocity));
      }
    }
  }
  
  // Forcing ---------------------------------------------------
	// Heat Exchanger
  if ( isComponent(HEX) ) {
    fprintf(fp, "\n\t[Heat Exchanger]\n");
    
    fprintf(fp, "\t no                    Label    ID    normal_x   normal_y   normal_z     O_x[m]     O_y[m]     O_z[m]      dir_x      dir_y      dir_z\n");
    for(n=1; n<=NoBC; n++) {
      if ( compo[n].getType() == HEX ) {
        fprintf(fp, "\t%3d %24s %5d  %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e\n", 
                n, compo[n].name, compo[n].getID(),
								compo[n].nv[0], compo[n].nv[1], compo[n].nv[2],
                compo[n].oc[0], compo[n].oc[1], compo[n].oc[2],
                compo[n].dr[0], compo[n].dr[1], compo[n].dr[2]);
      }
    }
    fprintf(fp, "\n");
    
    fprintf(fp, "\t                                     Depth[m]   Width[m]  Height[m]  Area[m*m]\n");
    for(n=1; n<=NoBC; n++) {
      if ( compo[n].getType() == HEX ) {
        fprintf(fp, "\t%3d %24s %5d %10.3e %10.3e %10.3e %10.3e\n", 
                n, compo[n].name, compo[n].getID(),
								compo[n].depth, compo[n].shp_p1, compo[n].shp_p2, compo[n].area);
      }
    }
    fprintf(fp, "\n");
    
    fprintf(fp, "\t                                      i_st    i_ed    j_st    j_ed    k_st    k_ed\n");
    for(n=1; n<=NoBC; n++) {
      if ( compo[n].getType() == HEX ) {
        fprintf(fp, "\t%3d %24s %5d %7d %7d %7d %7d %7d %7d\n", 
                n, compo[n].name, compo[n].getID(),
                getCmpGbbox_st_x(n, gci), getCmpGbbox_ed_x(n, gci), 
                getCmpGbbox_st_y(n, gci), getCmpGbbox_ed_y(n, gci), 
                getCmpGbbox_st_z(n, gci), getCmpGbbox_ed_z(n, gci));
      }
    }
    fprintf(fp, "\n");
    
    fprintf(fp, "\t no                    Label    ID         c1         c2         c3         c4  u_th[m/s]  thick[mm]     vec_forcing\n");
    for(n=1; n<=NoBC; n++) {
      if ( compo[n].getType() == HEX ) {
        fprintf(fp, "\t%3d %24s %5d %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e     %s\n", 
                n, compo[n].name, compo[n].getID(),
                compo[n].ca[0], compo[n].ca[1], compo[n].ca[2], compo[n].ca[3], compo[n].ca[4]*RefVelocity, compo[n].ca[5]*RefLength*1000.0,
                (compo[n].get_sw_HexDir()==ON) ? "Directional":"Non-directional");
      }
    }
  }
  
  // Fan
  if ( isComponent(FAN) ) {
    fprintf(fp, "\n\t[Fan]\n");
    
    fprintf(fp, "\t no                    Label    ID    normal_x   normal_y   normal_z      O_x[m]     O_y[m]     O_z[m]\n");
    for(n=1; n<=NoBC; n++) {
      if ( compo[n].getType() == FAN ) {
        fprintf(fp, "\t%3d %24s %5d  %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e\n", 
                n, compo[n].name, compo[n].getID(),
								compo[n].nv[0], compo[n].nv[1], compo[n].nv[2],
                compo[n].oc[0], compo[n].oc[1], compo[n].oc[2]);
      }
    }
    fprintf(fp, "\n");
    
    fprintf(fp, "\t                                     Depth[m]     Fan[m]    Boss[m]  Area[m*m]\n");
    for(n=1; n<=NoBC; n++) {
      if ( compo[n].getType() == FAN ) {
        fprintf(fp, "\t%3d %24s %5d %10.3e %10.3e %10.3e %10.3e\n", 
                n, compo[n].name, compo[n].getID(),
								compo[n].depth, compo[n].shp_p1, compo[n].shp_p2, compo[n].area);
      }
    }
    fprintf(fp, "\n");
    
    fprintf(fp, "\t                                      i_st    i_ed    j_st    j_ed    k_st    k_ed\n");
    for(n=1; n<=NoBC; n++) {
      if ( compo[n].getType() == FAN ) {
        fprintf(fp, "\t%3d %24s %5d %7d %7d %7d %7d %7d %7d\n", 
                n, compo[n].name, compo[n].getID(),
                getCmpGbbox_st_x(n, gci), getCmpGbbox_ed_x(n, gci), 
                getCmpGbbox_st_y(n, gci), getCmpGbbox_ed_y(n, gci), 
                getCmpGbbox_st_z(n, gci), getCmpGbbox_ed_z(n, gci));
      }
    }
    fprintf(fp, "\n");
    /*
     fprintf(fp, "\t no                    Label    ID         c1         c2         c3         c4  u_th[m/s]  thick[mm]     vec_forcing\n");
     for(n=1; n<=NoBC; n++) {
     if ( compo[n].getType() == FAN ) {
     fprintf(fp, "\t%3d %24s %5d %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e     %s\n", 
     n, compo[n].name, compo[n].getID(),
     compo[n].ca[0], compo[n].ca[1], compo[n].ca[2], compo[n].ca[3], compo[n].ca[4]*RefVelocity, compo[n].ca[5]*RefLength*1000.0,
     (compo[n].get_sw_HexDir()==ON) ? "Directional":"Non-directional");
     }
     }*/
  }
  
	// Darcy Law
  if ( isComponent(DARCY) ) {
    fprintf(fp, "\n\t[Darcy medium]\n");
    
    fprintf(fp, "\t no                    Label    ID        (Computed Normal form ID)  Area[m*m]   Prmblty_x   Prmblty_y   Prmblty_z[m^2]    Prmblty_x   Prmblty_y   Prmblty_z[-]\n");
    for(n=1; n<=NoBC; n++) {
      if ( compo[n].getType() == DARCY ) {
        fprintf(fp, "\t%3d %24s %5d %10.3e %10.3e %10.3e %10.3e %11.4e %11.4e %11.4e       %11.4e %11.4e %11.4e\n", 
                n, compo[n].name, compo[n].getID(), 
                nv[3*n+0], nv[3*n+1], nv[3*n+2], compo[n].area, 
								compo[n].ca[0], compo[n].ca[1], compo[n].ca[2], compo[n].ca[3], compo[n].ca[4], compo[n].ca[5]);
      }
    }
    
    fprintf(fp, "\t                                      i_st    i_ed    j_st    j_ed    k_st    k_ed\n");
    for(n=1; n<=NoBC; n++) {
      if ( compo[n].getType() == DARCY ) {
        fprintf(fp, "\t%3d %24s %5d %7d %7d %7d %7d %7d %7d\n", 
                n, compo[n].name, compo[n].getID(),
                getCmpGbbox_st_x(n, gci), getCmpGbbox_ed_x(n, gci), 
                getCmpGbbox_st_y(n, gci), getCmpGbbox_ed_y(n, gci), 
                getCmpGbbox_st_z(n, gci), getCmpGbbox_ed_z(n, gci));
      }
    }
  }
  
  // Heat Face ---------------------------------------------------
  // Adiabatic
  if ( HeatProblem && isComponent(ADIABATIC) ) {
    fprintf(fp, "\n\t[Adiabatic]\n");
    fprintf(fp, "\t no                    Label    ID   def\n");
    
    for(n=1; n<=NoBC; n++) {
      if ( compo[n].getType() == ADIABATIC ) {
        fprintf(fp, "\t%3d %24s %5d %5d\n", n, compo[n].name, compo[n].getID(), compo[n].getDef());
      }
    }
  }
  
  // Direct Heat Flux
  if ( HeatProblem && isComponent(HEATFLUX) ) {
    fprintf(fp, "\n\t[Direct Heat Flux]\n");
    fprintf(fp, "\t no                    Label    ID   def    i_st    i_ed    j_st    j_ed    k_st    k_ed   Area[m*m]   flux(W/m^2)        q[-]\n");
    
    for(n=1; n<=NoBC; n++) {
      if ( compo[n].getType() == HEATFLUX ) {
        fprintf(fp, "\t%3d %24s %5d %5d %7d %7d %7d %7d %7d %7d %11.4e %12.4e %12.4e\n", 
                n, compo[n].name, compo[n].getID(), compo[n].getDef(),
                getCmpGbbox_st_x(n, gci), getCmpGbbox_ed_x(n, gci), 
                getCmpGbbox_st_y(n, gci), getCmpGbbox_ed_y(n, gci), 
                getCmpGbbox_st_z(n, gci), getCmpGbbox_ed_z(n, gci), 
                compo[n].area, compo[n].get_Heatflux(), compo[n].get_Heatflux()/(RefVelocity*DiffTemp*rho*cp));
      }
    }
  }
  
  // Heat Transfer N
  if ( HeatProblem && isCompoTransfer(HT_N) ) {
    fprintf(fp, "\n\t[Heat Transfer : Type N]\n");
    fprintf(fp, "\t no                    Label    ID   def    i_st    i_ed    j_st    j_ed    k_st    k_ed   Area[m*m]   coef(W/m^2K)\n");
    
    for(n=1; n<=NoBC; n++) {
      if ( compo[n].getHtype() == HT_N ) {
        fprintf(fp, "\t%3d %24s %5d %5d %7d %7d %7d %7d %7d %7d %11.4e %12.4e\n", 
                n, compo[n].name, compo[n].getID(), compo[n].getDef(), 
                getCmpGbbox_st_x(n, gci), getCmpGbbox_ed_x(n, gci), 
                getCmpGbbox_st_y(n, gci), getCmpGbbox_ed_y(n, gci), 
                getCmpGbbox_st_z(n, gci), getCmpGbbox_ed_z(n, gci), 
                compo[n].area, compo[n].get_CoefHT());
      }
    }
  }
  
  // Heat Transfer S
  if ( HeatProblem && isCompoTransfer(HT_S) ) {
    fprintf(fp, "\n\t[Heat Transfer : Type S]\n");
    fprintf(fp, "\t no                    Label    ID   def    i_st    i_ed    j_st    j_ed    k_st    k_ed   Area[m*m]   coef(W/m^2K)   Temp(%s)   Temp[-]\n", (Unit_Temp==Unit_KELVIN) ? "K" : "C");
    
    for(n=1; n<=NoBC; n++) {
      if ( compo[n].getHtype() == HT_S ) {
        fprintf(fp, "\t%3d %24s %5d %5d %7d %7d %7d %7d %7d %7d %11.4e %12.4e %12.4e %12.4e\n", 
                n, compo[n].name, compo[n].getID(), compo[n].getDef(), 
                getCmpGbbox_st_x(n, gci), getCmpGbbox_ed_x(n, gci), 
                getCmpGbbox_st_y(n, gci), getCmpGbbox_ed_y(n, gci), 
                getCmpGbbox_st_z(n, gci), getCmpGbbox_ed_z(n, gci), 
                compo[n].area, compo[n].get_CoefHT(), FBUtility::convK2Temp(compo[n].get_Temp(), Unit_Temp),
                FBUtility::convK2ND(compo[n].get_Temp(), BaseTemp, DiffTemp) );
      }
    }
  }
  
  // Heat Transfer SN
  if ( HeatProblem && isCompoTransfer(HT_SN) ) {
    fprintf(fp, "\n\t[Heat Transfer : Type SN]\n"); 
    fprintf(fp, "\t no                    Label    ID   def    i_st    i_ed    j_st    j_ed    k_st    k_ed   Area[m*m]      Temp(%s)      Temp[-]   Type\n", (Unit_Temp==Unit_KELVIN) ? "K" : "C");
    
    for(n=1; n<=NoBC; n++) {
      if ( compo[n].getHtype() == HT_SN ) {
        fprintf(fp, "\t%3d %24s %5d %5d %7d %7d %7d %7d %7d %7d %11.4e %12.4e %12.4e   %s\n", 
                n, compo[n].name, compo[n].getID(), compo[n].getDef(), 
                getCmpGbbox_st_x(n, gci), getCmpGbbox_ed_x(n, gci), 
                getCmpGbbox_st_y(n, gci), getCmpGbbox_ed_y(n, gci), 
                getCmpGbbox_st_z(n, gci), getCmpGbbox_ed_z(n, gci), 
                compo[n].area, FBUtility::convK2Temp(compo[n].get_Temp(), Unit_Temp),
                FBUtility::convK2ND(compo[n].get_Temp(), BaseTemp, DiffTemp),
                (compo[n].get_sw_HTmodeRef()==CompoList::HT_mode_bulk) ? "Bulk" : "Local");
      }
    }
    fprintf(fp, "\t no                    Label    ID  vert_lam_a  vert_lam_b vert_turb_a vert_turb_b  vert_Ra_cr   lwr_lam_a   lwr_lam_b  lwr_turb_a  lwr_turb_b   lwr_Ra_cr\n");
    for(n=1; n<=NoBC; n++) {
      if ( compo[n].getHtype() == HT_SN ) {
        fprintf(fp, "\t%3d %24s %5d %11.4e %11.4e %11.4e %11.4e %11.4e %11.4e %11.4e %11.4e %11.4e %11.4e\n", 
                n, 
                compo[n].name, 
                compo[n].getID(), 
                compo[n].ca[CompoList::vert_laminar_alpha], 
                compo[n].ca[CompoList::vert_laminar_beta], 
                compo[n].ca[CompoList::vert_turbulent_alpha], 
                compo[n].ca[CompoList::vert_turbulent_beta], 
                compo[n].ca[CompoList::vert_Ra_critial], 
                compo[n].cb[CompoList::lower_laminar_alpha], 
                compo[n].cb[CompoList::lower_laminar_beta], 
                compo[n].cb[CompoList::lower_turbulent_alpha], 
                compo[n].cb[CompoList::lower_turbulent_beta], 
                compo[n].cb[CompoList::lower_Ra_critial] );
      }
    }
  }
  
  // Heat Transfer SF
  if ( HeatProblem && isCompoTransfer(HT_SF) ) {
    fprintf(fp, "\n\t[Heat Transfer : Type SF]\n");
    fprintf(fp, "\t no                    Label    ID   def    i_st    i_ed    j_st    j_ed    k_st    k_ed   Area[m*m]      Temp(%s)      Temp[-]   Type\n", 
            (Unit_Temp==Unit_KELVIN) ? "K" : "C");
    
    for(n=1; n<=NoBC; n++) {
      if ( compo[n].getHtype() == HT_SF ) {
        fprintf(fp, "\t%3d %24s %5d %5d %7d %7d %7d %7d %7d %7d %11.4e %12.4e %12.4e\n", 
                n, 
                compo[n].name, 
                compo[n].getID(), 
                compo[n].getDef(), 
                getCmpGbbox_st_x(n, gci), getCmpGbbox_ed_x(n, gci), 
                getCmpGbbox_st_y(n, gci), getCmpGbbox_ed_y(n, gci), 
                getCmpGbbox_st_z(n, gci), getCmpGbbox_ed_z(n, gci), 
                compo[n].area, FBUtility::convK2Temp(compo[n].get_Temp(), Unit_Temp),
                FBUtility::convK2ND(compo[n].get_Temp(), BaseTemp, DiffTemp),
                (compo[n].get_sw_HTmodeRef()==CompoList::HT_mode_bulk) ? "Bulk" : "Local");
      }
    }
    fprintf(fp, "\t no                    Label    ID       alpha        beta       gamma\n");
    for(n=1; n<=NoBC; n++) {
      if ( compo[n].getHtype() == HT_SN ) {
        fprintf(fp, "\t%3d %24s %5d %11.4e %11.4e %11.4e\n", 
                n, 
                compo[n].name, 
                compo[n].getID(), 
                compo[n].ca[CompoList::alpha], 
                compo[n].ca[CompoList::beta], 
                compo[n].ca[CompoList::gamma]);
      }
    }
  }
  
  // Heat Transfer B
  if ( HeatProblem && isCompoTransfer(HT_B)) {
    fprintf(fp, "\n\t[Heat Transfer : Type B]\n");
    fprintf(fp, "\t no                    Label    ID   def    i_st    i_ed    j_st    j_ed    k_st    k_ed   Area[m*m]  coef(W/m^2K)  BulkTemp(%s)   BulkTemp[-]\n", 
            (Unit_Temp==Unit_KELVIN) ? "K" : "C");
    
    for(n=1; n<=NoBC; n++) {
      if ( compo[n].getHtype() == HT_B ) {
        fprintf(fp, "\t%3d %24s %5d %5d %7d %7d %7d %7d %7d %7d %11.4e  %12.4e %12.4e %12.4e\n", 
                n, compo[n].name, compo[n].getID(), compo[n].getDef(), 
                getCmpGbbox_st_x(n, gci), getCmpGbbox_ed_x(n, gci), 
                getCmpGbbox_st_y(n, gci), getCmpGbbox_ed_y(n, gci), 
                getCmpGbbox_st_z(n, gci), getCmpGbbox_ed_z(n, gci), 
                compo[n].area, compo[n].get_CoefHT(), FBUtility::convK2Temp(compo[n].get_Temp(), Unit_Temp),
                FBUtility::convK2ND(compo[n].get_Temp(), BaseTemp, DiffTemp) );
      }
    }
  }
  
  // Iso Thermal
  if ( HeatProblem && isComponent(ISOTHERMAL)) {
    fprintf(fp, "\n\t[Iso-Thermal]\n");
    fprintf(fp, "\t no                    Label    ID   def    i_st    i_ed    j_st    j_ed    k_st    k_ed   Area[m*m]   Sf.Temp(%s)   Sf.Temp[-]\n", 
            (Unit_Temp==Unit_KELVIN) ? "K" : "C");
    
    for(n=1; n<=NoBC; n++) {
      if ( compo[n].getType() == ISOTHERMAL ) {
        fprintf(fp, "\t%3d %24s %5d %5d %7d %7d %7d %7d %7d %7d %11.4e %12.4e %12.4e \n", 
                n, compo[n].name, compo[n].getID(), compo[n].getDef(), 
                getCmpGbbox_st_x(n, gci), getCmpGbbox_ed_x(n, gci), 
                getCmpGbbox_st_y(n, gci), getCmpGbbox_ed_y(n, gci), 
                getCmpGbbox_st_z(n, gci), getCmpGbbox_ed_z(n, gci), 
								compo[n].area, FBUtility::convK2Temp(compo[n].get_Temp(), Unit_Temp),
                FBUtility::convK2ND(compo[n].get_Temp(), BaseTemp, DiffTemp) );
      }
    }
  }
  
  // Radiant
  if ( HeatProblem && isComponent(RADIANT)) {
    fprintf(fp, "\n\t[Radiant]\n");
    fprintf(fp, "\t no                    Label    ID    def    i_st    i_ed    j_st    j_ed    k_st    k_ed   Area[m*m]   ep[-]   pj[-]\n");
    
    for(n=1; n<=NoBC; n++) {
      if ( compo[n].getType() == RADIANT ) {
        fprintf(fp, "\t%3d %24s %5d %5d %7d %7d %7d %7d %7d %7d %11.4e %12.4e %12.4e\n", 
                n, compo[n].name, compo[n].getID(), compo[n].getDef(), 
                getCmpGbbox_st_x(n, gci), getCmpGbbox_ed_x(n, gci), 
                getCmpGbbox_st_y(n, gci), getCmpGbbox_ed_y(n, gci), 
                getCmpGbbox_st_z(n, gci), getCmpGbbox_ed_z(n, gci), 
                compo[n].area, compo[n].get_CoefRadEps(), compo[n].get_CoefRadPrj());
      }
    }
  }
  
  // Heat source ---------------------------------------------------
  // Heat generation
  if ( HeatProblem && isComponent(HEAT_SRC)) {
    fprintf(fp, "\n\t[Heat Generation]\n");
    fprintf(fp, "\t no                    Label    ID    i_st    i_ed    j_st    j_ed    k_st    k_ed     Q[W/m^3]    nrmlzd[-]\n");
    
    for(n=1; n<=NoBC; n++) {
      unsigned h_odr = compo[n].getMatOdr();
      if ( compo[n].getType() == HEAT_SRC ) {
        fprintf(fp, "\t%3d %24s %5d %7d %7d %7d %7d %7d %7d %12.4e %12.4e\n", 
                n, compo[n].name, compo[n].getID(), 
                getCmpGbbox_st_x(n, gci), getCmpGbbox_ed_x(n, gci), 
                getCmpGbbox_st_y(n, gci), getCmpGbbox_ed_y(n, gci), 
                getCmpGbbox_st_z(n, gci), getCmpGbbox_ed_z(n, gci), 
                compo[n].get_HeatValue(), 
                FBUtility::convD2ND_Hsrc(compo[n].get_HeatValue(), RefVelocity, RefLength, DiffTemp, mat[h_odr].P[p_density], mat[h_odr].P[p_specific_heat]));
      }
    }
  }
  
  // Constant Temperature
  if ( HeatProblem && isComponent(CNST_TEMP)) {
    fprintf(fp, "\n\t[Constant Temperature]\n");
    fprintf(fp, "\t no                    Label    ID    i_st    i_ed    j_st    j_ed    k_st    k_ed      Temp[%s]      Temp[-]\n", 
            (Unit_Temp==Unit_KELVIN) ? "K" : "C");
    
    for(n=1; n<=NoBC; n++) {
      if ( compo[n].getType() == CNST_TEMP ) {
        fprintf(fp, "\t%3d %24s %5d %7d %7d %7d %7d %7d %7d %12.4e %12.4e\n", 
                n, compo[n].name, compo[n].getID(), 
                getCmpGbbox_st_x(n, gci), getCmpGbbox_ed_x(n, gci), 
                getCmpGbbox_st_y(n, gci), getCmpGbbox_ed_y(n, gci), 
                getCmpGbbox_st_z(n, gci), getCmpGbbox_ed_z(n, gci), 
                FBUtility::convK2Temp(compo[n].get_Temp(), Unit_Temp), 
                FBUtility::convK2ND(compo[n].get_Temp(), BaseTemp, DiffTemp));
      }
    }
  }
  
  // Monitor ---------------------------------------------------
  if ( isComponent(CELL_MONITOR) ) {
    fprintf(fp, "\n\t[Monitor]\n");
    fprintf(fp, "\t no                    Label    ID    i_st    i_ed    j_st    j_ed    k_st    k_ed        (Computed Normal form ID)  Area[m*m]  Variables\n");
    
    for(n=1; n<=NoBC; n++) {
      if ( compo[n].getType() == CELL_MONITOR ) {
        fprintf(fp, "\t%3d %24s %5d %7d %7d %7d %7d %7d %7d %10.3e %10.3e %10.3e %10.3e  %s\n", 
                n, compo[n].name, compo[n].getID(), 
                getCmpGbbox_st_x(n, gci), getCmpGbbox_ed_x(n, gci), 
                getCmpGbbox_st_y(n, gci), getCmpGbbox_ed_y(n, gci), 
                getCmpGbbox_st_z(n, gci), getCmpGbbox_ed_z(n, gci), 
								nv[3*n+0], nv[3*n+1], nv[3*n+2], compo[n].area, compo[n].getVarStr().c_str());
      }
    }
    fprintf(fp, "\n");
    
    fprintf(fp, "\t no                    Label    ID   normal_x   normal_y   normal_z Reference\n");
    for(n=1; n<=NoBC; n++){
      if ( compo[n].getType() == CELL_MONITOR )  {
        fprintf(fp,"\t%3d %24s %5d %10.3e %10.3e %10.3e       %3s\n",
                n, compo[n].name, compo[n].getID(), compo[n].nv[0], compo[n].nv[1], compo[n].nv[2],
                (compo[n].getStateCellMonitor()==ON) ? "yes" : "no");
      }
    }
  }
  
  // Periodic ---------------------------------------------------
  if ( isComponent(PERIODIC) ) {
    fprintf(fp, "\n\t[Periodic]\n");
    fprintf(fp, "\t no                    Label    ID    i_st    i_ed    j_st    j_ed    k_st    k_ed    Pressure Difference [Pa]/[-]  Driver\n");
    
    int dir_in=0, dir_out=0, pp_in=0, pp_out=0;
    
    for(n=1; n<=NoBC; n++) {
      if ( compo[n].getType() == PERIODIC ) {
        dir_in = compo[n].getPeriodicDir();
        
        fprintf(fp, "\t%3d %24s %5d %7d %7d %7d %7d %7d %7d     ", 
                n, compo[n].name, compo[n].getID(), 
                getCmpGbbox_st_x(n, gci), getCmpGbbox_ed_x(n, gci), 
                getCmpGbbox_st_y(n, gci), getCmpGbbox_ed_y(n, gci), 
                getCmpGbbox_st_z(n, gci), getCmpGbbox_ed_z(n, gci));
        fprintf(fp,"%12.6e / %12.6e ", compo[n].ca[0], FBUtility::convD2ND_P(compo[n].ca[0], BasePrs, RefDensity, RefVelocity, Unit_Prs));
        fprintf(fp, "%7s\n", FBUtility::getDirection(dir_in).c_str());
        
        // ドライバの方向が周期境界であるかをチェック
        if ( bc[dir_in].get_BCtype() != OBC_PERIODIC ) {
          fprintf(fp, "\tError : Specified driver direction[%s] by component is not PERIODIC.", FBUtility::getDirection(dir_in).c_str());
          Exit(0);
        }
        
        // ドライバの方向が一致しているかどうかをチェック
        dir_out = bc[dir_in].get_DriverDir();
        if ( dir_in != dir_out ) {
          fprintf(fp, "\tError : The specification of driver direction is different between Outer(%d) and Inner(%d) in XML.", dir_out, dir_in);
          Exit(0);
        }
        
        // 入力のインデクス値とボクセルの位置が一致しているかどうかをチェック
        switch ( dir_in ) {
          case X_MINUS:
          case X_PLUS:
            pp_in = getCmpGbbox_st_x(n, gci);
            break;
            
          case Y_MINUS:
          case Y_PLUS:
            pp_in = getCmpGbbox_st_y(n, gci);
            break;
            
          case Z_MINUS:
          case Z_PLUS:
            pp_in = getCmpGbbox_st_z(n, gci);
            break;
        }
        
        pp_out = bc[dir_in].get_DriverIndex();
        if ( pp_in != pp_out ) {
          fprintf(fp, "\tError : The index of inner periodic cell is different in between XML(%d) and voxel model(%d).", pp_out, pp_in);
          Exit(0);
        }
      }
    }
    fprintf(fp, "\n");
  }
  
  fflush(fp);
}


/**
 @fn void ParseBC::printCompoInfo(FILE* mp, FILE* fp, REAL_TYPE* nv, int* gci, MaterialList* mat)
 @brief コンポーネントの情報を表示する
 @param mp 標準出力
 @param nv ボクセルモデルから計算した法線
 @param gci グローバルなコンポーネントのインデクス
 @param mat MaterialList
 */
void ParseBC::printCompoInfo(FILE* mp, FILE* fp, REAL_TYPE* nv, int* gci, MaterialList* mat)
{
  if( !fp || !mp ) Exit(0);
  
  printCompo(mp, nv, gci, mat);
  printCompo(fp, nv, gci, mat);
}


/**
 @fn void ParseBC::printFaceOBC(FILE* fp, REAL_TYPE* G_Lbx)
 @brief 外部境界条件の各面の情報を表示する
 @param fp
 @param G_Lbx グローバルの領域の大きさ
 */
void ParseBC::printFaceOBC(FILE* fp, REAL_TYPE* G_Lbx)
{
  for (unsigned i=0; i<NOFACE; i++) {
    fprintf(fp,"\t      Set %s up as %s : (%s)\n", Control::getDirection(i).c_str(), getOBCstr(bc[i].get_BCtype()).c_str(), OBCname[i]);
    printOBC(fp, &bc[i], G_Lbx, i);
    fprintf(fp,"\n");
  }
  fflush(fp);
}

/**
 @fn void ParseBC::printOBC(FILE* fp, BoundaryOuter* ref REAL_TYPE* G_Lbx, unsigned face)
 @brief 速度の外部境界条件処理の表示
 @param fp
 @param ref
 @param G_Lbx グローバルの領域の大きさ
 @param face 面番号
 */
void ParseBC::printOBC(FILE* fp, BoundaryOuter* ref, REAL_TYPE* G_Lbx, unsigned face)
{
  REAL_TYPE a, b, c;
  
  fprintf(fp,"\t\t\tGuide Cell ID = %d / Medium = %d\n", ref->get_GuideID(), ref->get_GuideMedium());
  
  switch ( ref->get_BCtype() ) {
    case OBC_WALL:
      if ( ref->get_vType() == CompoList::vel_harmonic ) {
        fprintf(fp,"\t\t\tVelocity  Harmonic Oscillation\n");
        fprintf(fp,"\t\t\t          nrml vec. V(%10.3e, %10.3e, %10.3e) [-]\n", ref->nv[0], ref->nv[1], ref->nv[2]);
        fprintf(fp,"\t\t\t          amp.    %12.6e [m/s] / %12.6e [-]\n", ref->ca[0], ref->ca[0]/RefVelocity);
        fprintf(fp,"\t\t\t          freq.   %12.6e [Hz]  / %12.6e [-]\n", ref->ca[1], ref->ca[1]*RefLength/RefVelocity);
        fprintf(fp,"\t\t\t          phase   %12.6e [rad] / %12.6e [-]\n", ref->ca[2], ref->ca[2]);
        fprintf(fp,"\t\t\t          intcpt. %12.6e [m/s] / %12.6e [-]\n", ref->ca[3], ref->ca[3]/RefVelocity);
      }
      else { // vel_constant, vel_zero
        c = ref->ca[CompoList::bias]; // 有次元値
        fprintf(fp,"\t\t\tVelocity V(%10.3e, %10.3e, %10.3e) [m/s] / V(%10.3e, %10.3e, %10.3e) [-]\n", 
                ref->nv[0]*c, ref->nv[1]*c, ref->nv[2]*c, ref->nv[0]*c/RefVelocity, ref->nv[1]*c/RefVelocity, ref->nv[2]*c/RefVelocity); 
      }
      fprintf(fp,"\t\t\tPressure Gradient is %s\n", ( Mode_Gradp == P_GRAD_ZERO) ? "zero" : "calculated from Navier-Stokes eqs.");
      if ( HeatProblem ) {
        unsigned htp = ref->get_hType();
        if ( htp == ADIABATIC ) {
          fprintf(fp,"\t\t\tAdiabatic\n");
        }
        else if ( htp == TRANSFER) {
          unsigned ht_mode = ref->get_HTmode();
          if ( ht_mode == HT_N ) {
            fprintf(fp, "\tHeat Transfer Type N  : H. T. Coef. = %e \n", ref->get_CoefHT());
          }
          else if ( ht_mode == HT_S ) {
            fprintf(fp, "\tHeat Transfer Type S  : H. T. Coef. = %e \n", ref->get_CoefHT());
            fprintf(fp, "\t                        Surf. temp. = %e \n", ref->get_Temp());
          }
          else if ( ht_mode == HT_SN ) {
            fprintf(fp, "\tHeat Transfer Type SN : Surf. temp. = %e \n", ref->get_Temp());
            fprintf(fp, "\t                        Ref. Temp.  = %s \n", (ref->get_HTmodeRef()==CompoList::HT_mode_bulk) ? "Bulk" : "Local");
            fprintf(fp, "\t                        vertival_laminar_alpha   = %12.6e \n", ref->ca[0]);
            fprintf(fp, "\t                        vertival_laminar_beta    = %12.6e \n", ref->ca[1]);
            fprintf(fp, "\t                        vertival_turbulent_alpha = %12.6e \n", ref->ca[2]);
            fprintf(fp, "\t                        vertival_turbulent_beta  = %12.6e \n", ref->ca[3]);
            fprintf(fp, "\t                        vertival_Ra_critial      = %12.6e \n", ref->ca[4]);
            fprintf(fp, "\t                        lower_laminar_alpha      = %12.6e \n", ref->cb[0]);
            fprintf(fp, "\t                        lower_laminar_beta       = %12.6e \n", ref->cb[1]);
            fprintf(fp, "\t                        lower_turbulent_alpha    = %12.6e \n", ref->cb[2]);
            fprintf(fp, "\t                        lower_turbulent_beta     = %12.6e \n", ref->cb[3]);
            fprintf(fp, "\t                        lower_Ra_critial         = %12.6e \n", ref->cb[4]);
          }
          else if ( ht_mode == HT_SF ) {
            fprintf(fp, "\tHeat Transfer Type SN : Surf. temp. = %e \n", ref->get_Temp());
            fprintf(fp, "\t                        Ref. Temp.  = %s \n", (ref->get_HTmodeRef()==CompoList::HT_mode_bulk) ? "Bulk" : "Local");
            fprintf(fp, "\t                        alpha       = %12.6e \n", ref->ca[0]);
            fprintf(fp, "\t                        beta        = %12.6e \n", ref->ca[1]);
            fprintf(fp, "\t                        gamma       = %12.6e \n", ref->ca[2]);
          }
          else if ( ht_mode == HT_B ) {
            fprintf(fp, "\tHeat Transfer Type B  : H. T. Coef. = %e \n", ref->get_CoefHT());
            fprintf(fp, "\t                        Bulk temp.  = %e \n", ref->get_Temp());
          }
        }
        else if ( htp == HEATFLUX ) {
          fprintf(fp,"Heat Flux %e\n", ref->get_Heatflux());
        }
        else if ( htp == ISOTHERMAL ) {
          fprintf(fp,"Isothermal\n");
        }
        //else if ( htp == CNST_TEMP ) {
        //  fprintf(fp,"Dirichlet %e [%s]\n", ref->T1.c, (Unit_Temp==CompoList::Unit_KELVIN)?"K":"C");
        //}        
      }
      break;
      
    case OBC_SPEC_VEL:
      if ( ref->get_vType() == CompoList::vel_harmonic ) {
        fprintf(fp,"\t\t\tVelocity  Harmonic Oscillation\n");
        fprintf(fp,"\t\t\t          nrml vec. V(%10.3e, %10.3e, %10.3e) [-]\n", ref->nv[0], ref->nv[1], ref->nv[2]);
        fprintf(fp,"\t\t\t          amp.    %12.6e [m/s] / %12.6e [-]\n", ref->ca[0], ref->ca[0]/RefVelocity);
        fprintf(fp,"\t\t\t          freq.   %12.6e [Hz]  / %12.6e [-]\n", ref->ca[1], ref->ca[1]*RefLength/RefVelocity);
        fprintf(fp,"\t\t\t          phase   %12.6e [rad] / %12.6e [-]\n", ref->ca[2], ref->ca[2]);
        fprintf(fp,"\t\t\t          intcpt. %12.6e [m/s] / %12.6e [-]\n", ref->ca[3], ref->ca[3]/RefVelocity);
      }
      else { // vel_zero, vel_constant
        c = ref->ca[CompoList::bias];
        fprintf(fp,"\t\t\tVelocity V(%10.3e, %10.3e, %10.3e) [m/s] / V(%10.3e, %10.3e, %10.3e) [-]\n", 
                ref->nv[0]*c, ref->nv[1]*c, ref->nv[2]*c, ref->nv[0]*c/RefVelocity, ref->nv[1]*c/RefVelocity, ref->nv[2]*c/RefVelocity);
      }
      
      fprintf(fp,"\t\t\tPressure Gradient is zero.\n");
      
      if ( HeatProblem ) {
        fprintf(fp, "\t\t\tSpecified Temperature  = %12.6e [%s] / %12.6e [-] \n", 
                FBUtility::convK2Temp(ref->get_Temp(), Unit_Temp),
                (Unit_Temp==Unit_KELVIN) ? "K" : "C", 
                FBUtility::convK2ND(ref->get_Temp(), BaseTemp, DiffTemp));
      }
      break;
      
    case OBC_OUTFLOW:
      fprintf(fp,"\t\t\tOutflow with %s convective velocity\n", (ref->get_ofv() == V_AVERAGE) ? "Average" : "Minmax" );
      if (ref->get_pType() == P_DIRICHLET) {
        fprintf(fp,"\t\t\t%12.6e [Pa]  /  %12.6e\n", ref->p, FBUtility::convD2ND_P(ref->p, BasePrs, RefDensity, RefVelocity, Unit_Prs));        
      }
      else {
        fprintf(fp,"\t\t\tPressure Gradient is zero\n");
      }
      if ( HeatProblem ) {
      }
      fprintf(fp,"\t\t\tNumber of Element = %d\n",ref->get_ValidCell());
      break;
      
    case OBC_TRC_FREE:
      //fprintf(fp,"\t\t\t%12.6e [Pa]  /  %12.6e [-]\n", ref->p, FBUtility::convD2ND_P(ref->p, BasePrs, RefDensity, RefVelocity, Unit_Prs));
      
      if ( HeatProblem ) {
        fprintf(fp, "\t\t\t    Ambient Temperature  = %12.6e \n", ref->get_Temp());
      }
      break;
      
    case OBC_FAR_FIELD:
      break;
      
    case OBC_PERIODIC:
      switch (ref->get_PrdcMode()) {
        case BoundaryOuter::prdc_Simple:
          fprintf(fp, "\t\t\tSimple periodic copy\n");
          break;
          
        case BoundaryOuter::prdc_Directional:
          fprintf(fp,"\t\t\tPressure Difference = %12.6e [Pa]    /  %12.6e [-]\n", 
                  ref->p, FBUtility::convD2ND_P(ref->p, BasePrs, RefDensity, RefVelocity, Unit_Prs));
          
          switch (face) {
            case X_MINUS:
            case X_PLUS:
              fprintf(fp,"\t\t\tPressure Gradient   = %12.6e [Pa/m]  /  %12.6e [-]\n", 
                      ref->p/(G_Lbx[0]*RefLength), FBUtility::convD2ND_P(ref->p, BasePrs, RefDensity, RefVelocity, Unit_Prs)/G_Lbx[0]);
              break;
              
            case Y_MINUS:
            case Y_PLUS:
              fprintf(fp,"\t\t\tPressure Gradient   = %12.6e [Pa/m]  /  %12.6e [-]\n", 
                      ref->p/(G_Lbx[1]*RefLength), FBUtility::convD2ND_P(ref->p, BasePrs, RefDensity, RefVelocity, Unit_Prs)/G_Lbx[1]);
              break;
              
            case Z_MINUS:
            case Z_PLUS:
              fprintf(fp,"\t\t\tPressure Gradient   = %12.6e [Pa/m]  /  %12.6e [-]\n", 
                      ref->p/(G_Lbx[2]*RefLength), FBUtility::convD2ND_P(ref->p, BasePrs, RefDensity, RefVelocity, Unit_Prs)/G_Lbx[2]);
              break;
          }          
          break;
          
        case BoundaryOuter::prdc_Driver:
          fprintf(fp,"\t\t\tDriver direction    : %s\n", FBUtility::getDirection(ref->get_DriverDir()).c_str());
          break;
          
        default:
          Exit(0);
      }
      break;
      
    case OBC_SYMMETRIC:
      break;
      
    default:
      printf("\n\tError : OuterBC\n");
      Exit(0);
      break;
  }
  
  fflush(fp);
}


/**
 @fn void ParseBC::printOBCinfo(FILE* mp, FILE* fp, REAL_TYPE* G_Lbx)
 @brief 外部境界を表示
 @param G_Lbx グローバルの領域の大きさ
 */
void ParseBC::printOBCinfo(FILE* mp, FILE* fp, REAL_TYPE* G_Lbx)
{
  
#ifdef DEBUG
  printBaseOBC(mp);
  printBaseOBC(fp);
#endif
  
  printFaceOBC(mp, G_Lbx);
  printFaceOBC(fp, G_Lbx);
}


/**
 @fn void ParseBC::printTable(FILE* fp)
 @brief iTableに登録された情報を表示
 */
void ParseBC::printTable(FILE* fp)
{
  fprintf(fp,"\t iTable order : State    ID  MatID Label\n");
  for (int i=1; i<=NoID; i++) {
    fprintf(fp,"\t        %5d : %s %5d  %5d %s\n", i, (iTable[i].getState()==SOLID) ? "Solid" : "Fluid", 
            iTable[i].getID(), iTable[i].getMatID(), iTable[i].getLabel());
  }
}

/**
 @fn void ParseBC::receiveCompoPtr(CompoList* CMP)
 @brief コンポーネントリストのポインタを受け取る
 */
void ParseBC::receiveCompoPtr(CompoList* CMP)
{
  if ( !CMP ) {
    stamped_printf("\tAn object of CompoList class is NULL\n");
    Exit(0);
  }
  compo = CMP;
}

/**
 @fn void ParseBC::receiveCfgPtr(SklSolverConfig* cfg)
 @brief コンフィギュレーションリストのポインタを受け取る
 */
void ParseBC::receiveCfgPtr(SklSolverConfig* cfg) 
{ 
  if ( !cfg ) {
    stamped_printf("\tAn object of Configuration Tree is NULL\n");
    Exit(0);
  }
  CF = cfg;
}

/**
 @fn void ParseBC::setObcPtr(BoundaryOuter* ptr)
 @brief BoundaryOuterクラスのポインタを受け取り，内部作業用のBoundaryOuterクラスをインスタンスする
 @note
  - BoundaryOuterクラスのポインタを受け取る
  - XMLファイルの外部境界条件をパースし，個数を取得する
  - 内部作業用のBoundaryOuterクラスをインスタンスする
 */
void ParseBC::setObcPtr(BoundaryOuter* ptr) 
{ 
  if ( !ptr ) Exit(0);
  bc = ptr;

  // XMLファイルの基本外部境界条件をパースし，個数を取得する
  const CfgElem *elemTop, *elmL1;
  elemTop = elmL1 = NULL;
  
  if( !(elemTop = CF->GetTop(OUTERBND)) ) {
    stamped_printf("\tParsing error : Missing OuterBoundary tree\n");
    Exit(0);
  }
  if( !(elmL1 = elemTop->GetElemFirst("Basic_BCs")) ) {
    printf("\tParsing error : No Basic_BCs was found.\n");
    Exit(0);
  }
  
  NoBaseBC = (unsigned)elmL1->GetElemSize();
  BaseBc = new BoundaryOuter[NoBaseBC];
}

/**
 @fn void ParseBC::setControlVars(Control* Cref) 
 @brief 
    - ParseBCに必要な変数をコピーする
    - Controlクラスのメンバ変数に値をコピーする
 @param Cref Controlクラスのポインタ
 */
void ParseBC::setControlVars(Control* Cref)
{
  ix          = Cref->imax;
  jx          = Cref->jmax;
  kx          = Cref->kmax;
  guide       = Cref->guide;
  RefVelocity = Cref->RefVelocity;
  BaseTemp    = Cref->BaseTemp;
  DiffTemp    = Cref->DiffTemp;
  RefDensity  = Cref->RefDensity;
  HeatProblem = Cref->isHeatProblem();
  Unit_Param  = Cref->Unit.Param;
  RefLength   = Cref->RefLength;
  KindOfSolver= Cref->KindOfSolver;
  Unit_Temp   = Cref->Unit.Temp;
  Unit_Prs    = Cref->Unit.Prs;
	BasePrs     = Cref->BasePrs;
  Mode_Gradp  = Cref->Mode.PrsNeuamnnType;
  isCDS       = Cref->isCDS();
  
  Cref->NoBC    = NoBC    = getNBC(); // XMLファイルのInnerBoundaryタグ内の各境界条件タグから境界条件の個数を取得
  Cref->NoID    = NoID    = scanXMLmodel(); // XMLファイルのModel_Settingタグ内にあるID数
  Cref->NoCompo = NoCompo = NoBC + NoID;

  unsigned m, s;
  s = MASK_6; // bit幅マスクは2^(bit幅)-1を表し，ちょうど0を除いた個数となっている
  m = log10(s+1)/log10(2);
  if ( NoCompo > s ) {
    printf("Error : No. of Component (= NoBC + NoID) must be less or equal %d(%dbit-width)\n", s, m);
    Exit(0);
  }

  s = MASK_5-1; // 0と31を除く
  m = log10(s+2)/log10(2);
  if ( NoBC > s ) {
    printf("Error : No. of BC must be less or equal %d(%dbit-width)\n", s, m);
    Exit(0);
  }

  // iTableのアロケート
  iTable = new IDtable[NoID+1];
}


/**
 @fn unsigned ParseBC::scanXMLmodel(void)
 @brief XMLに記述されたモデル情報の数を取得
 @retval モデルの数，またはエラーコード0
 */
unsigned ParseBC::scanXMLmodel(void)
{
  const CfgElem *elemTop=NULL, *elmL1=NULL;
  unsigned n;
  
  elemTop = CF->GetTop(STEER);
  if( !(elmL1 = elemTop->GetElemFirst("Model_Setting")) ) {
    stamped_printf("\tParsing error : Missing the section of 'Model_Setting'\n");
    Exit(0);
  }
  n = (unsigned)elmL1->GetParamSize();
  
  return n;
}


/**
 @fn const CfgElem* ParseBC::selectFace(int face, const CfgElem* elemTop)
 @brief 外部境界面をパースする
 @param face フェイス番号
 @param elemTop
 */
const CfgElem* ParseBC::selectFace(int face, const CfgElem* elemTop)
{
  const CfgElem *elmL1 = NULL;
  const char* cmt=NULL;
  
  switch (face) {
    case X_MINUS:
      if( !(elmL1 = elemTop->GetElemFirst("X_MINUS")) ) {
        stamped_printf("\tParsing error : Missing the section of 'X_MINUS'\n");
        Exit(0);
      }
      break;
      
    case X_PLUS:
      if( !(elmL1 = elemTop->GetElemFirst("X_PLUS")) ) {
        stamped_printf("\tParsing error : Missing the section of 'X_PLUS'\n");
        Exit(0);
      }
      break;
      
    case Y_MINUS:
      if( !(elmL1 = elemTop->GetElemFirst("Y_MINUS")) ) {
        stamped_printf("\tParsing error : Missing the section of 'Y_MINUS'\n");
        Exit(0);
      }
      break;
      
    case Y_PLUS:
      if( !(elmL1 = elemTop->GetElemFirst("Y_PLUS")) ) {
        stamped_printf("\tParsing error : Missing the section of 'Y_PLUS'\n");
        Exit(0);
      }
      break;
      
    case Z_MINUS:
      if( !(elmL1 = elemTop->GetElemFirst("Z_MINUS")) ) {
        stamped_printf("\tParsing error : Missing the section of 'Z_MINUS'\n");
        Exit(0);
      }
      break;
      
    case Z_PLUS:
      if( !(elmL1 = elemTop->GetElemFirst("Z_PLUS")) ) {
        stamped_printf("\tParsing error : Missing the section of 'Z_PLUS'\n");
        Exit(0);
      }
      break;
      
    default:
      stamped_printf("\tParsing error : the section of 'Face_BC'\n");
      Exit(0);
  }
  
  if ( !(cmt = elmL1->GetComment()) ) {
    stamped_printf("\tParsing error : No comment in 'Face_BC' section\n");
    elmL1 = NULL; // means error
    return elmL1;
  }
  strcpy(OBCname[face], cmt);
  
  return elmL1;
}
/**
 @fn void ParseBC::setCompoList(Control* C)
 @brief CompoListに内部境界条件の情報を設定する
 @param C
 @pre ParseBC::getNumberOfBC(Control* C), ParseBC::getXML_Model(Control* Cref)
 @note
    - 最初にBCの情報を登録，その後IDの情報を登録
    - XMLファイルから各内部BCのidをパースし，compoに保持する
    - 格納番号は1からスタート
 */
void ParseBC::setCompoList(Control* C)
{ 
  const CfgElem *elmL1, *elmL2;
  elmL1 = elmL2 = NULL;
  int ide;
  unsigned tp;
  const char *pnt=NULL, *ename=NULL;
  
  // 境界条件があれば，cfgツリーのポインタを辿る
  if ( NoBC > 0 ) {
    elmL1 = CF->GetTop(INNERBND);
    if ( !( elmL2  = elmL1->GetElemFirst() ) ) Exit(0);
  }
  
  // 境界条件がなければ，スキップ
  for (unsigned odr=1; odr<=NoBC; odr++) {
    
    if ( !(ename = elmL2->GetName()) ) {
      stamped_printf("\tParsing error : No Elem name in 'InnerBoundary'\n");
      Exit(0);
    }
    
    // cmp[].type, h_typeのセット
    chkKeywordIBC(ename, odr);
    
    // IDの取得
    if ( !elmL2->isSetID() ) {
      stamped_printf("\tParsing error : No ID section in '%s'\n", ename);
      Exit(0);
    }
    if ( -1 == (ide=elmL2->GetID()) ) {
      stamped_printf("\tParsing error : No valid ID in '%s'\n", ename);
      Exit(0);
    }
    
    // IDがiTableの中にリストアップされているかを調べる
    if ( !isIDinTable(ide) ) {
      stamped_printf("\tParsing error : ID[%d] described in '%s' is not included in 'Model_Setting'\n", ide, ename);
      Exit(0);
    }
    
    // Labelの取得．ラベルなしでもエラーではない
    pnt = NULL;
    if ( !(pnt = elmL2->GetComment()) ) {
      stamped_printf("\tWarning : No commnet in '%s'\n", ename);
    }
    else {
      compo[odr].setName(pnt);
    }
    
    // とりあえず登録，BCのIDの重複は続く処理で確認
    compo[odr].setID((unsigned)ide);
    
    // stateの登録
    compo[odr].setState( getStateinTable( compo[odr].getID() ) );
    
    // 各BCの処理
    tp = compo[odr].getType();
    
    if ( tp == SPEC_VEL ) {
      if ( compo[odr].getState() != SOLID ) {
        printf("\tID Error : SPEC_VEL ID must be Solid\n");
        Exit(0);
      }
      getXML_IBC_SpecVel(elmL2, odr);
      if ( !isIDinCompo(ide, compo[odr].getDef(), odr) ) {
        stamped_printf("Parse Error : Reduplication of a pair of ID[%d] and Def[%d] for BC\n", ide, compo[odr].getDef());
        Exit(0);
      }      
    }
    else if ( tp == OUTFLOW ) {
      if ( compo[odr].getState() != SOLID ) {
        printf("\tID Error : OUTFLOW ID must be Solid\n");
        Exit(0);
      }
      getXML_IBC_Outflow(elmL2, odr);
      if ( !isIDinCompo(ide, compo[odr].getDef(), odr) ) {
        stamped_printf("Parse Error : Reduplication of a pair of ID[%d] and Def[%d] for BC\n", ide, compo[odr].getDef());
        Exit(0);
      }      
    }
    else if ( tp == OUTFLOW ) {
      getXML_IBC_IBM_DF(elmL2, odr);
      if ( !isIDinCompo(ide, odr) ) {
        stamped_printf("Parse Error : Reduplication of ID[%d] for BC\n", ide);
        Exit(0);
      }
    }
    else if ( tp == HEX ) {
      getXML_IBC_PrsLoss(elmL2, odr);
      if ( !isIDinCompo(ide, odr) ) {
        stamped_printf("Parse Error : Reduplication of ID[%d] for BC\n", ide);
        Exit(0);
      }
    }
    else if ( tp == FAN ) {
      getXML_IBC_Fan(elmL2, odr);
      if ( !isIDinCompo(ide, odr) ) {
        stamped_printf("Parse Error : Reduplication of ID[%d] for BC\n", ide);
        Exit(0);
      }
    }
    else if ( tp == DARCY ) {
      getDarcy(elmL2, odr);
      if ( !isIDinCompo(ide, odr) ) {
        stamped_printf("Parse Error : Reduplication of ID[%d] for BC\n", ide);
        Exit(0);
      }
    }
    else if ( tp == CELL_MONITOR ) {
      getXML_Cell_Monitor(elmL2, odr, C);
      if ( !isIDinCompo(ide, odr) ) {
        stamped_printf("Parse Error : Reduplication of ID[%d] for BC\n", ide);
        Exit(0);
      }
    }
    else if ( tp == INACTIVE ) {
      if ( !isIDinCompo(ide, odr) ) { // IDのみの重複チェック
        stamped_printf("Parse Error : Reduplication of ID[%d] for BC\n", ide);
        Exit(0);
      }
    }
    else if ( tp == PERIODIC ) {
      getXML_IBC_Periodic(elmL2, odr);
      if ( !isIDinCompo(ide, odr) ) {
        stamped_printf("Parse Error : Reduplication of ID[%d] for BC\n", ide);
        Exit(0);
      }
    }
    else if ( HeatProblem ) { // Incase of Heat problem
      if ( C->KindOfSolver == FLOW_ONLY ) {
        stamped_printf("Parse Error : Heat BC is not allowed on FLOW_ONLY mode.\n");
        Exit(0);
      }
      
      if ( tp == ADIABATIC ) {
        getXML_IBC_Adiabatic(elmL2, odr);
        if ( !isIDinCompo(ide, compo[odr].getDef(), odr) ) {
          stamped_printf("Parse Error : Reduplication of pair of ID[%d] and Def[%d] for BC\n", ide, compo[odr].getDef());
          Exit(0);
        }        
      }
      else if ( tp == HEATFLUX ) {
        getXML_IBC_HeatFlux(elmL2, odr);
        if ( !isIDinCompo(ide, compo[odr].getDef(), odr) ) {
          stamped_printf("Parse Error : Reduplication of pair of ID[%d] and Def[%d] for BC\n", ide, compo[odr].getDef());
          Exit(0);
        }        
      }
      else if ( tp == TRANSFER ) {
        switch ( compo[odr].getHtype() ) {
          case HT_N:
            getXML_IBC_HT_N(elmL2, odr);
            if ( !isIDinCompo(ide, compo[odr].getDef(), odr) ) {
              stamped_printf("Parse Error : Reduplication of pair of ID[%d] and Def[%d] for BC\n", ide, compo[odr].getDef());
              Exit(0);
            }
            break;
            
          case HT_S:
            getXML_IBC_HT_S(elmL2, odr);
            if ( !isIDinCompo(ide, compo[odr].getDef(), odr) ) {
              stamped_printf("Parse Error : Reduplication of pair of ID[%d] and Def[%d] for BC\n", ide, compo[odr].getDef());
              Exit(0);
            }
            break;
            
          case HT_SN:
            getXML_IBC_HT_SN(elmL2, odr);
            if ( !isIDinCompo(ide, compo[odr].getDef(), odr) ) {
              stamped_printf("Parse Error : Reduplication of pair of ID[%d] and Def[%d] for BC\n", ide, compo[odr].getDef());
              Exit(0);
            }
            break;
            
          case HT_SF:
            getXML_IBC_HT_SF(elmL2, odr);
            if ( !isIDinCompo(ide, compo[odr].getDef(), odr) ) {
              stamped_printf("Parse Error : Reduplication of pair of ID[%d] and Def[%d] for BC\n", ide, compo[odr].getDef());
              Exit(0);
            }
            break;
            
          case HT_B:
            getXML_IBC_HT_B(elmL2, odr);
            if ( !isIDinCompo(ide, compo[odr].getDef(), odr) ) {
              stamped_printf("Parse Error : Reduplication of pair of ID[%d] and Def[%d] for BC\n", ide, compo[odr].getDef());
              Exit(0);
            }
            break;
        }        
      }
      else if ( tp == ISOTHERMAL ) {
        getXML_IBC_IsoTherm(elmL2, odr);
        if ( !isIDinCompo(ide, compo[odr].getDef(), odr) ) {
          stamped_printf("Parse Error : Reduplication of pair of ID[%d] and Def[%d] for BC\n", ide, compo[odr].getDef());
          Exit(0);
        }
      }
      else if ( tp == RADIANT ) {
        getXML_IBC_Radiant(elmL2, odr);
        if ( !isIDinCompo(ide, odr) ) {
          stamped_printf("Parse Error : Reduplication of ID[%d] for BC\n", ide);
          Exit(0);
        }
      }
      else if ( tp == HEAT_SRC ) {
        getXML_IBC_HeatSrc(elmL2, odr);
        if ( !isIDinCompo(ide, odr) ) {
          stamped_printf("Parse Error : Reduplication of ID[%d] for BC\n", ide);
          Exit(0);
        }
      }
      else if ( tp == CNST_TEMP ) {
        getXML_IBC_CnstTemp(elmL2, odr);
      }
      else {
        printf("\tError : Invalid Inner BC keyword [%d]\n", tp);
        Exit(0);
      }
    }
    elmL2 = elmL1->GetElemNext(elmL2);
  }
    
  // 媒質情報の登録
  for (int i=1; i<=NoID; i++) {
    compo[NoBC+i].setID   ( iTable[i].getID() );
    compo[NoBC+i].setState( iTable[i].getState() );
    compo[NoBC+i].setName ( iTable[i].getLabel() );
  }
}

/**
 @fn void ParseBC::set_Deface(const CfgElem *elmL, unsigned n, const char* str)
 @brief 内部境界条件のdef_faceを取得し，登録する
 @param elmL
 @param n コンポーネントリストのエントリ番号
 @param str エラー表示用のキーワード
 */
void ParseBC::set_Deface(const CfgElem *elmL, unsigned n, const char* str)
{
  int def=0;
  
  if ( !elmL->GetValue(CfgIdt("def_face"), &def) ) {
		stamped_printf("\tParsing error : fail to get 'def_face' in '%s'\n", str);
    Exit(0);
  }
	if ( !isIDinTable(def) ) { // IDがiTableの中にリストアップされているかを調べる
		stamped_printf("\tParsing error : ID[%d] described in 'def_face' is not included in 'Model_Setting'\n", def);
		Exit(0);
	}
	compo[n].setDef(def);
}

/**
 @fn void ParseBC::setMedium(Control* Cref)
 @brief KOSと媒質の状態の整合性をチェックし，媒質数をカウント，C.NoMediumFluid, C.NoMediumSolid, C.NoMaterialをセット
 @note
 - 流体の媒質は少なくとも一つは必要
 - iTableの数の上限は，NoID
 */
void ParseBC::setMedium(Control* Cref)
{
  // check at least one fluid
  if ( KindOfSolver != SOLID_CONDUCTION ) {
    bool check=false;
    for (int i=1; i<=NoID; i++) {
      if ( iTable[i].getState() == FLUID ) check = true;
    }
    if ( !check ) {
      stamped_printf("\tVoxel model should have at least one FLUID.\n");
      Exit(0);
    }
  }
  
  // 流体と固体の媒質数をセット
  unsigned m_fluid=0, m_solid=0;
  for (int i=1; i<=NoID; i++) {
    if ( iTable[i].getState() == SOLID ) m_solid++;
    else m_fluid++;
  }
  
  Cref->NoMediumFluid = m_fluid;
  Cref->NoMediumSolid = m_solid;
  
  // MaterialIDの数をカウント
  Cref->NoMaterial = getNoMaterial();
}


/**
 @fn void ParseBC::setRefValue(MaterialList* mat, CompoList* cmp, Control* C)
 @brief 媒質により決まる代表量をコピーする
 @param mat MaterialList class
 @param cmp CompoList class
 @param C Control class
 */
void ParseBC::setRefValue(MaterialList* mat, CompoList* cmp, Control* C)
{
  unsigned m;
  
	for (unsigned n=NoBC+1; n<=NoCompo; n++) {
    if ( cmp[n].getID() == C->RefID ) {
      m = cmp[n].getMatOdr();
      if ( mat[m].getState() == FLUID ) {
        RefDensity      = mat[m].P[p_density];
        //mu    = mat[m].P[p_viscosity];
        //nyu   = mat[m].P[p_kinematic_viscosity];
        RefSpecificHeat = mat[m].P[p_specific_heat];
        //lambda= mat[m].P[p_thermal_conductivity];
        //beta  = mat[m].P[p_vol_expansion]; // can be replaced by 1/K in the case of gas
        //snd_spd = mat[m].P[p_sound_of_speed];
      }
      else {
        RefDensity      = mat[m].P[p_density];
        RefSpecificHeat = mat[m].P[p_specific_heat];
        //lambda = mat[m].P[p_thermal_conductivity];
      }
    }
  }
}

/**
 @fn void ParseBC::setRefMedium(MaterialList* mat, Control* Cref)
 @brief 指定した媒質IDから参照物理量を求める
 @param mat MaterialList
 @param Cref Control class
 */
void ParseBC::setRefMedium(MaterialList* mat, Control* Cref)
{
  unsigned m;
  
  for (unsigned n=Cref->NoBC+1; n<=Cref->NoCompo; n++) {
    if ( compo[n].getID() == Cref->RefID ) {
      m = compo[n].getMatOdr();
      if ( mat[m].getState() == FLUID ) {
        rho    = mat[m].P[p_density];
        nyu    = mat[m].P[p_kinematic_viscosity];
        cp     = mat[m].P[p_specific_heat];
        lambda = mat[m].P[p_thermal_conductivity];
        beta   = mat[m].P[p_vol_expansion]; // can be replaced by 1/K in the case of gas
      }
      else {
        rho    = mat[m].P[p_density];
        cp     = mat[m].P[p_specific_heat];
        lambda = mat[m].P[p_thermal_conductivity];
      }
    }
  }
}
