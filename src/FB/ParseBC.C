// #################################################################
//
// CAERU Library
//
// Copyright (c) All right reserved. 2012
//
// Institute of Industrial Science, The University of Tokyo, Japan. 
//
// #################################################################

//@file ParseBC.C
//@brief FlowBase ParseBC class
//@author kero

#include <math.h>
#include "ParseBC.h"

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
        Hostonly_ stamped_printf("\tNo consistency between 'Kind_of_Solver' and 'Local_Boundary'\n");
        Exit(0);
      }
    }
  }
  else if (kos == SOLID_CONDUCTION) {
    for (int n=1; n<NoBC; n++) {
      if ( compo[n].isVBC() ) {
        Hostonly_ stamped_printf("\tNo consistency between 'Kind_of_Solver' and 'Local_Boundary'\n");
        Exit(0);
      }
    }
  }
}

//@brief ラベルの重複を調べる
bool ParseBC::chkDuplicate(const int n, const std::string m_label)
{
	for (int i=0; i<n; i++){
    if ( BaseBc[i].get_Alias() == m_label ) return false;
	}
	return true;
}


/**
 @fn void ParseBC::countMedium(Control* Cref)
 @brief KOSと媒質の状態の整合性をチェックし，媒質数をカウント，C.NoMediumFluid, C.NoMediumSolidをセット
 @note 流体の媒質は少なくとも一つは必要
 */
void ParseBC::countMedium(Control* Cref)
{
  // check at least one fluid
  if ( KindOfSolver != SOLID_CONDUCTION ) {
    bool check=false;
    for (int i=1; i<=NoMedium; i++) {
      if ( mat[i].getState() == FLUID ) check = true;
    }
    if ( !check ) {
      Hostonly_ stamped_printf("\tAnalysis model should have at least one FLUID.\n");
      Exit(0);
    }
  }
  
  // 流体と固体の媒質数をセット
  unsigned m_fluid=0, m_solid=0;
  for (int i=1; i<=NoMedium; i++) {
    if ( mat[i].getState() == SOLID ) m_solid++;
    else m_fluid++;
  }
  
  Cref->NoMediumFluid = m_fluid;
  Cref->NoMediumSolid = m_solid;
}



/**
 @fn string ParseBC::getOBCstr(const int id)
 @brief 外部境界条件のキーワードを照合し， BCの文字列を返す
 @param id 
 */
std::string ParseBC::getOBCstr(const int id)
{
  std::string bc;
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
 @fn int ParseBC::getNoLocalBC(void)
 @brief LocalBoundaryタグ直下のBCの個数（内部境界条件数）を返す
 @note 境界条件数がゼロでもエラーではない
 */
int ParseBC::getNoLocalBC(void)
{  
  int nobc=0;
  std::string str;
  std::string label;
  
  label="/BC_Table/LocalBoundary";
  
  if ( tpCntl->chkNode(label) ) { //nodeがあれば
	  nobc = tpCntl->countLabels(label);
  }
  
  return nobc;
}



/**
 @fn REAL_TYPE ParseBC::get_BCval_real(const string label)
 @brief 境界条件の値(REAL_TYPE型)を取得し，返す
 @param label
 */
REAL_TYPE ParseBC::get_BCval_real(const std::string label)
{
  REAL_TYPE df=0.0f;
  if ( !(tpCntl->GetValue(label, &df )) ) {
    stamped_printf("\tParsing error : Invalid REAL_TYPE value for '%s'\n", label.c_str());
    Exit(0);
  }
  return df;
}



/**
 @fn void ParseBC::get_Center(const std::string label_base, const int n, REAL_TYPE* v)
 @brief 内部境界条件の座標値を取得し，登録する
 @param label_base
 @param n オーダー
 @param v[out] ベクトルパラメータ
 @param str エラー表示用文字列
 */
void ParseBC::get_Center(const std::string label_base, const int n, REAL_TYPE* v)
{
  std::string label;
  for (int i=0; i<3; i++) v[i]=0.0f;
  
  label = label_base + "/Center";
  
  if( !(tpCntl->GetVector(label, v, 3)) )
  {
    stamped_printf("\tParsing error : fail to get vec params in '%s\n", label.c_str());
    Exit(0);
  }
}


/**
 @fn void ParseBC::get_Darcy(const std::string label_base, const int n)
 @brief Darcyのパラメータを取得する
 @param label_base Forcing_Volumeのレベル
 @param n コンポーネントリストのエントリ番号
 @note 透過率[m^2]は境界条件設定時に無次元化する
 */
void ParseBC::get_Darcy(const std::string label_base, const int n)
{
  REAL_TYPE v[3];
  string label;
  int nnode=0;
  
  for (unsigned n=0; n<3; n++) v[n]=0.0;
  
  // check number of Elem
  nnode = tpCntl->countLabels(label_base);
  if ( nnode != 3) {
    stamped_printf("\tParsing error : 3 params should be found in 'Darcy' : %d\n", nnode);
    Exit(0);
  }
  
  // 透過率の取得
  label = label_base + "/permeability";
  
  if( !(tpCntl->GetVector(label, v, 3)) ) {
    stamped_printf("\tParsing error : fail to get permeability params in 'Darcy'\n");
    Exit(0);
  }
  compo[n].ca[0] = v[0];
  compo[n].ca[1] = v[1];
  compo[n].ca[2] = v[2];
}



/**
 @fn void ParseBC::get_Dir(const std::string label_base, const int n, REAL_TYPE* v)
 @brief 内部境界条件の方向ベクトル値を取得し，登録する
 @param label_base
 @param n オーダー
 @param v[out] ベクトルパラメータ
 @param str エラー表示用文字列
 */
void ParseBC::get_Dir(const std::string label_base, const int n, REAL_TYPE* v)
{
  std::string label;
  for (int i=0; i<3; i++) v[i]=0.0f;
  
  label = label_base + "/Dir";
  if( !(tpCntl->GetVector(label, v, 3)) )
  {
    stamped_printf("\tParsing error : fail to get vec params in '%s\n", label.c_str());
    Exit(0);
  }
  
  //単位ベクトル化
  getUnitVec(v);
}



/**
 @fn void ParseBC::get_IBC_IBM_DF(const std::string label_base, const int n)
 @brief Direct Forcingのパラメータを取得する
 @param label_base 
 @param n コンポーネントリストのエントリ番号
 */
void ParseBC::get_IBC_IBM_DF(const std::string label_base, const int n)
{
  int d;
  int nnode=0;
  REAL_TYPE v[3];
  std::string str;
  string label;
  
  // check number of Elem
  nnode=tpCntl->countLabels(label_base);
  if ( nnode != 4) {
    stamped_printf("\tParsing error : 4 params should be found in 'LocalBoundary > Forcing' : %d\n", nnode);
    Exit(0);
  }
  
  // 法線ベクトル
  get_NV(label_base, n, v);
  copyVec(compo[n].nv, v);
  
  // Velocity
  label=label_base+"/Velocity";//
  //std::cout <<  "label : " << label << std::endl;
  REAL_TYPE ct = get_BCval_real(label);
  if ( Unit_Param == DIMENSIONAL ) {
    compo[n].set_Velocity( ct );
  }
  else {
    compo[n].set_Velocity( ct * RefVelocity );
  }
}


/**
 @fn void ParseBC::get_IBC_PrsLoss(const std::string label_base, const int n)
 @brief HeatExchangerのパラメータを取得する
 @param label_base コンフィギュレーションツリーのポインタ
 @param n コンポーネントリストのエントリ番号
 @note この時点ではRefDensityの値が未定なので，あとでパラメータ処理
 @see Control::setParameters()
 */
void ParseBC::get_IBC_PrsLoss(const std::string label_base, const int n)
{
  std::string str,str_u;
  string label;
  REAL_TYPE v[4], ct;
  int nnode=0;
  
  // check number of Elem
  nnode=tpCntl->countLabels(label_base);
  if ( nnode != 11) {
    stamped_printf("\tParsing error : 11 params should be found in 'Pressure_Loss'\n");
    Exit(0);
  }
  
  // 入力単位の指定
  label=label_base+"/unit";//
  //std::cout <<  "label : " << label << std::endl;
  if ( !(tpCntl->GetValue(label, &str_u )) ) {
		stamped_printf("\tParsing error : Invalid float value for 'unit' in 'Pressure_Loss'\n");
		Exit(0);
  }
  if ( !strcasecmp("mmaq", str_u.c_str()) ) {
    compo[n].setPrsUnit(CompoList::unit_mmAq);
  }
  else if ( !strcasecmp("mmhg", str_u.c_str()) ) {
    compo[n].setPrsUnit(CompoList::unit_mmHg);
  }
  else if ( !strcasecmp("pa", str_u.c_str()) ) {
    compo[n].setPrsUnit(CompoList::unit_Pa);
  }
  else if ( !strcasecmp("non_dimension", str_u.c_str()) ) {
    compo[n].setPrsUnit(CompoList::unit_NonDimensional);
  }
  else {
    stamped_printf("\tDescribed unit is out of scope.\n");
    Exit(0);
  }
  
  // 法線ベクトルの取得
  get_NV(label_base, n, v);
  copyVec(compo[n].nv, v);
  
  // 方向ベクトルの取得
  get_Dir(label_base, n, v);
  copyVec(compo[n].dr, v);
  
  // 中心座標の取得
  get_Center(label_base, n, v);
  copyVec(compo[n].oc, v);
  
  // 形状パラメータ
  label = label_base + "/depth";
  compo[n].depth  = get_BCval_real(label);
  
  label = label_base + "/width";
  compo[n].shp_p1 = get_BCval_real(label);
  
  label = label_base + "/height";
  compo[n].shp_p2 = get_BCval_real(label);
  
  // 圧力損失パラメータ
  //label=label_base+"/c1";
  //compo[n].ca[0] = get_BCval_real(label);
  //label=label_base+"/c2";
  //compo[n].ca[1] = get_BCval_real(label);
  //label=label_base+"/c3";
  //compo[n].ca[2] = get_BCval_real(label);
  //label=label_base+"/c4";
  //compo[n].ca[3] = get_BCval_real(label);
  label=label_base+"/c";
  if( !(tpCntl->GetVector(label, v, 4)) )
  {
    stamped_printf("\tParsing error : fail to get vec params in '%s\n", label.c_str());
    Exit(0);
  }
  compo[n].ca[0]=v[0];
  compo[n].ca[1]=v[1];
  compo[n].ca[2]=v[2];
  compo[n].ca[3]=v[3];
  
  label=label_base+"/u_threshold";
  compo[n].ca[4] = get_BCval_real(label);
  label=label_base+"/thickness";
  compo[n].ca[5] = get_BCval_real(label);
  
  // 熱交換器の方向強制オプション
  
  label=label_base+"/vector";//
  if ( !(tpCntl->GetValue(label, &str )) ) {
    stamped_printf("\tParsing error : Invalid string for 'vector' in 'Pressure_Loss'\n");
    Exit(0);
  }
  if ( !strcasecmp("directional", str.c_str()) ) {
    compo[n].set_sw_HexDir( ON );
  }
  else {
    compo[n].set_sw_HexDir( OFF );
  }
  
}



/**
 @fn void ParseBC::get_IBC_Fan(const std::string label_base, const int n)
 @brief Fanのパラメータを取得する
 @param label_base Forcing_Volumeのレベル
 @param n コンポーネントリストのエントリ番号
 */
void ParseBC::get_IBC_Fan(const std::string label_base, const int n)
{
  std::string str,str_u;
  string label;
  REAL_TYPE v[4], ct;
  int nnode=0;
  
  // check number of Elem
  nnode=tpCntl->countLabels(label_base);
  if ( nnode != 6) {
    stamped_printf("\tParsing error : 1 param should be found in Heat_Volume > Heat_Generation\n");
    Exit(0);
  }
  
  // 入力単位の指定
  label=label_base+"/unit";//
  //std::cout <<  "label : " << label << std::endl;
  if ( !(tpCntl->GetValue(label, &str_u )) ) {
		stamped_printf("\tParsing error : Invalid float value for 'unit' in 'Pressure_Loss'\n");
		Exit(0);
  }
  if ( !strcasecmp("mmaq", str_u.c_str()) ) {
    compo[n].setPrsUnit(CompoList::unit_mmAq);
  }
  else if ( !strcasecmp("mmhg", str_u.c_str()) ) {
    compo[n].setPrsUnit(CompoList::unit_mmHg);
  }
  else if ( !strcasecmp("pa", str_u.c_str()) ) {
    compo[n].setPrsUnit(CompoList::unit_Pa);
  }
  else if ( !strcasecmp("non_dimension", str_u.c_str()) ) {
    compo[n].setPrsUnit(CompoList::unit_NonDimensional);
  }
  else {
    stamped_printf("\tDescribed unit is out of scope.\n");
    Exit(0);
  }
  
  // 法線ベクトルの取得
  get_NV(label_base, n, v);
  copyVec(compo[n].nv, v);
  
  // 中心座標の取得
  get_Center(label_base, n, v);
  copyVec(compo[n].oc, v);
  
  // 形状パラメータ
  label=label_base+"/depth";
  compo[n].depth  = get_BCval_real(label);
  label=label_base+"/fan_radius";
  compo[n].shp_p1 = get_BCval_real(label);
  label=label_base+"/boss_radius";
  compo[n].shp_p2 = get_BCval_real(label);
  
  if ( compo[n].shp_p1 <= compo[n].shp_p2 ) {
    stamped_printf("\tError : Radius of boss is greater than fan.\n");
    Exit(0);
  }
  
}



/**
 @fn void ParseBC::get_IBC_Monitor(const std::string label_base, const int n, Control* C)
 @brief XMLファイルからMonitorの設定内容をパースし，パラメータを保持する
 @param label_base コンフィギュレーションツリーのポインタ
 @param n コンポーネントリストに登録するエントリ番号のベース
 @param C Control class
 @note Referenceは，隠しコマンドに
 */
void ParseBC::get_IBC_Monitor(const std::string label_base, const int n, Control* C)
{
  int nvc=0;
  REAL_TYPE v[3];
  std::string str,pnt;
  string label,label_leaf;
  
  // モードと形状
  label=label_base+"/shape";//
  //std::cout <<  "label : " << label << std::endl;
  if ( !(tpCntl->GetValue(label, &pnt )) ) {
  }
  else{
    if ( !strcasecmp("cylinder", pnt.c_str()) ) {
      compo[n].set_Shape(SHAPE_CYLINDER);
    }
    else if ( !strcasecmp("box", pnt.c_str()) ) {
      compo[n].set_Shape(SHAPE_BOX);
    }
    else if ( !strcasecmp("voxel_model", pnt.c_str()) ) {
      compo[n].set_Shape(SHAPE_VOXEL);
    }
    else {
      stamped_printf("\tParsing error : Invalid string value for 'Shape' : %s\n", pnt.c_str());
      Exit(0);
    }
  }
  
  // reference 
  label=label_base+"/reference";//
  //std::cout <<  "label : " << label << std::endl;
  if ( !(tpCntl->GetValue(label, &pnt )) ) {
  }
  else{
    if ( !strcasecmp("yes", pnt.c_str()) ) {
      compo[n].setStateCellMonitor(ON);
    }
    else if ( !strcasecmp("no", pnt.c_str()) ) {
      compo[n].setStateCellMonitor(OFF);
    }
    else {
      stamped_printf("\tParsing error : Invalid string value for 'reference' : %s\n", pnt.c_str());
      Exit(0);
    }
  }
  
  // 法線ベクトル
  get_NV(label_base, n, v);
  copyVec(compo[n].nv, v);
  
  unsigned shp = compo[n].get_Shape();
  
  if ( (shp == SHAPE_BOX) || (shp == SHAPE_CYLINDER) ) {
    // 中心座標の取得
    get_Center(label_base, n, v);
    copyVec(compo[n].oc, v);
    
  }
  
  if ( shp == SHAPE_BOX ) {
    // 方向ベクトルの取得
    get_Dir(label_base, n, v);
    copyVec(compo[n].dr, v);
    
    // 形状パラメータ
    label=label_base+"/depth";
    compo[n].depth  = get_BCval_real(label);
    label=label_base+"/width";
    compo[n].shp_p1 = get_BCval_real(label);
    label=label_base+"/height";
    compo[n].shp_p2 = get_BCval_real(label);
  }
  
  if ( shp == SHAPE_CYLINDER ) {
    // 形状パラメータ
    label=label_base+"/depth";
    compo[n].depth  = get_BCval_real(label);
    label=label_base+"/radius";
    compo[n].shp_p1 = get_BCval_real(label);
  }
  
  
  // Variables
  label_leaf=label_base+"/Variables";
  if ( !(tpCntl->chkNode(label_leaf)) ) {
    stamped_printf("\tParsing error : No 'Variables' keyword in 'Cell_Monitor'\n");
    Exit(0);
  }
  
  // サンプリングモード
  label=label_leaf+"/sampling_mode";//
  if ( !(tpCntl->GetValue(label, &pnt )) ) {
  }
  else{
    if ( !strcasecmp("all", pnt.c_str()) ) {
      compo[n].set_SamplingMode(SAMPLING_ALL);
    }
    else if ( !strcasecmp("fluid", pnt.c_str()) ) {
      compo[n].set_SamplingMode(SAMPLING_FLUID_ONLY);
    }
    else if ( !strcasecmp("solid", pnt.c_str()) ) {
      compo[n].set_SamplingMode(SAMPLING_SOLID_ONLY);
    }
    else {
      stamped_printf("\tParsing error : Invalid string value for 'Sampling_Mode' : %s\n", pnt.c_str());
      Exit(0);
    }
  }
  
  // サンプリング方法
  label=label_leaf+"/sampling_method";//
  if ( !(tpCntl->GetValue(label, &pnt )) ) {
  }
  else{
    if ( !strcasecmp("Nearest", pnt.c_str()) ) {
      compo[n].set_SamplingMethod(SAMPLING_NEAREST);
    }
    else if ( !strcasecmp("Interpolation", pnt.c_str()) ) {
      compo[n].set_SamplingMethod(SAMPLING_INTERPOLATION);
    }
    else if ( !strcasecmp("Smoothing", pnt.c_str()) ) {
      compo[n].set_SamplingMethod(SAMPLING_SMOOTHING);
    }
    else {
      stamped_printf("\tParsing error : Invalid string value for 'Sampling_Method' : %s\n", pnt.c_str());
      Exit(0);
    }
  }
  
  // モニタする変数と数を取得
  nvc = 0;
  
  // 速度
  label=label_leaf+"/velocity";//
  if ( !(tpCntl->GetValue(label, &str )) ) {
    stamped_printf("\tParsing error : fail to get 'Velocity' in 'Cell_Monitor'\n");
    Exit(0);
  }
  if     ( !strcasecmp(str.c_str(), "on") )  {
    compo[n].encodeVarType(var_Velocity);
    nvc++;
  }
  else if( !strcasecmp(str.c_str(), "off") ) {;} // nothing
  else {
    stamped_printf("\tInvalid keyword is described for 'Velocity'\n");
    Exit(0);
  }
  
  // 圧力
  label=label_leaf+"/pressure";//
  if ( !(tpCntl->GetValue(label, &str )) ) {
    stamped_printf("\tParsing error : fail to get 'Pressure' in 'Cell_Monitor'\n");
    Exit(0);
  }
  if     ( !strcasecmp(str.c_str(), "on") )  {
    compo[n].encodeVarType(var_Pressure);
    nvc++;
  }
  else if( !strcasecmp(str.c_str(), "off") ) {;} // nothing
  else {
    stamped_printf("\tInvalid keyword is described for 'Pressure'\n");
    Exit(0);
  }
  
  // 温度
  if ( HeatProblem ) {
    label=label_leaf+"/temperature";//
    if ( !(tpCntl->GetValue(label, &str )) ) {
      stamped_printf("\tParsing error : fail to get 'Temperature' in 'Cell_Monitor'\n");
      Exit(0);
    }
    if     ( !strcasecmp(str.c_str(), "on") )  {
      compo[n].encodeVarType(var_Temperature);
      nvc++;
    }
    else if( !strcasecmp(str.c_str(), "off") ) {;} // nothing
    else {
      stamped_printf("\tInvalid keyword is described for 'Temperature'\n");
      Exit(0);
    }
  }
  
  // 全圧
  label=label_leaf+"/total_pressure";//
  if ( !(tpCntl->GetValue(label, &str )) ) {
    stamped_printf("\tParsing error : fail to get 'Total_Pressure' in 'Cell_Monitor'\n");
    Exit(0);
  }
  if     ( !strcasecmp(str.c_str(), "on") )  {
    compo[n].encodeVarType(var_TotalP);
    nvc++;
  }
  else if( !strcasecmp(str.c_str(), "off") ) {;} // nothing
  else {
    stamped_printf("\tInvalid keyword is described for 'Total_Pressure'\n");
    Exit(0);
  }
  
  // モニタ面に対して指定された変数の個数（モニタの個数）を取得
  compo[n].setAttrb(nvc);
}



/**
 @fn void ParseBC::get_IBC_Periodic(const std::string label_base, const int n)
 @brief 内部の周期境界のパラメータを取得する
 @param label_base レベル
 @param n コンポーネントリストのエントリ番号
 */
void ParseBC::get_IBC_Periodic(const std::string label_base, const int n)
{
  int dir=0;
  REAL_TYPE ct=0.0;
  std::string str;
  std::string label;
  
  // 上流側の方向
  label = label_base + "/upstream_direction";

  if ( !(tpCntl->GetValue(label, &str )) ) {
	  printf("\tParsing error : fail to get 'upstream_direction' in 'LocalBoundary > Periodic'\n");
	  Exit(0);
  }
  if ( !strcasecmp("x_minus", str.c_str()) ) {
		dir = X_MINUS;
	}
  else if ( !strcasecmp("x_plus", str.c_str()) ) {
    dir = X_PLUS;
  }
  else if ( !strcasecmp("y_minus", str.c_str()) ) {
    dir = Y_MINUS;
  }
  else if ( !strcasecmp("y_plus", str.c_str()) ) {
    dir = Y_PLUS;
  }
  else if ( !strcasecmp("z_minus", str.c_str()) ) {
    dir = Z_MINUS;
  }
  else if ( !strcasecmp("z_plus", str.c_str()) ) {
    dir = Z_PLUS;
  }
  else {
    stamped_printf("\tParsing error : Invalid direction in 'LocalBoundary > Periodic'\n");
    Exit(0);
  }
	compo[n].setPeriodicDir((unsigned)dir);
  
  
  // 圧力差
  label = label_base + "/pressure_difference";
  
  if ( !(tpCntl->GetValue(label, &ct )) ) {
    printf("\tParsing error : Invalid value of 'Pressure difference' in 'LocalBoundary > Periodic'\n");
    Exit(0);
  }
  else {
    compo[n].ca[0] = ct;
  }
  
  // 面指定
  set_Deface(label_base, n);
}



/**
 @fn void ParseBC::get_IBC_Adiabatic(const std::string label_base, const int n)
 @brief Adiabaticのパラメータを取得する
 @param label_base
 @param n コンポーネントリストのエントリ番号
 */
void ParseBC::get_IBC_Adiabatic(const std::string label_base, const int n)
{
  std::string str,str_u;
  std::string label;
  REAL_TYPE v[4], ct;
  int nnode=0;
  
  // check number of Elem
  nnode = tpCntl->countLabels(label_base);
  
  if ( nnode != 1) {
    stamped_printf("\tParsing error : 1 param should be found in 'LocalBoundary > Adiabatic'\n");
    Exit(0);
  }
  
  // 面指定
  set_Deface(label_base, n);
  
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
 @fn void ParseBC::get_IBC_HeatFlux(const std::string label_base, const int n)
 @brief Direct_Fluxのパラメータを取得する
 @param label_base
 @param n コンポーネントリストのエントリ番号
 @note [W/m^2]
 */
void ParseBC::get_IBC_HeatFlux(const std::string label_base, const int n)
{
  std::string label;
  int nnode=0;
  
  // check number of Elem
  nnode = tpCntl->countLabels(label_base);
  
  if ( nnode != 2) {
    stamped_printf("\tParsing error : 2 params should be found in 'LocalBoundary > Direct_Heat_Flux'\n");
    Exit(0);
  }
  
  // 面指定
  set_Deface(label_base, n);
  
  label = label_base + "/Heat_Flux";
  
  compo[n].set_Heatflux( get_BCval_real(label) );
  
  if ( Unit_Param != DIMENSIONAL ) {
    stamped_printf("\tWarning: Heat condition must be a dimensional value\n");
    Exit(0);
  }
}


/**
 @fn void ParseBC::get_IBC_HT_N(const std::string label_base, const int n)
 @brief HeatTransfer_Nのパラメータを取得する
 @param label_base レベル
 @param n コンポーネントリストのエントリ番号
 */
void ParseBC::get_IBC_HT_N(const std::string label_base, const int n)
{
  std::string label;
  int nnode=0;
  
  // check number of Elem
  nnode = tpCntl->countLabels(label_base);
  
  if ( nnode != 2) {
    stamped_printf("\tParsing error : 2 params should be found in 'LocalBoundary > HeatTransfer_N'\n");
    Exit(0);
  }
  
  // 面指定
  set_Deface(label_base, n);
  
  // 熱伝達係数
  label = label_base + "/Coef_of_Heat_Transfer";
  
  compo[n].set_CoefHT( get_BCval_real(label) );
  
  if ( Unit_Param != DIMENSIONAL ) {
    stamped_printf("\tWarning: Heat condition must be a dimensional value\n");
    Exit(0);
  }
}


/**
 @fn void ParseBC::get_IBC_HT_S(const std::string label_base, const int n)
 @brief HeatTransfer_Sのパラメータを取得する
 @param label_base
 @param n コンポーネントリストのエントリ番号
 */
void ParseBC::get_IBC_HT_S(const std::string label_base, const int n)
{
  std::string label;
  int nnode=0;
  
  // check number of Elem
  nnode = tpCntl->countLabels(label_base);
  
  if ( nnode != 3) {
    stamped_printf("\tParsing error : 3 params should be found in 'HeatTransfer_S'\n");
    Exit(0);
  }
  
  // 面指定
  set_Deface(label_base, n);
  
  // 熱伝達係数
  label = label_base + "/Coef_of_Heat_Transfer";
  
  compo[n].set_CoefHT( get_BCval_real(label) );
  
  // 表面温度
  label = label_base + "/Surface_Temperature";
  
  REAL_TYPE st = get_BCval_real(label);
  compo[n].set_Temp( FBUtility::convTemp2K(st, Unit_Temp) );
  
  if ( Unit_Param != DIMENSIONAL ) {
    stamped_printf("\tWarning: Heat condition must be a dimensional value\n");
    Exit(0);
  }
}



/**
 @fn void ParseBC::get_IBC_HT_SN(const std::string label_base, const int n)
 @brief HeatTransfer_SNのパラメータを取得する
 @param label_base HeatTransfer_SNのレベル
 @param n コンポーネントリストのエントリ番号
 */
void ParseBC::get_IBC_HT_SN(const std::string label_base, const int n)
{
  std::string str;
  std::string label;
  int nnode=0;
  
  // check number of Elem
  nnode = tpCntl->countLabels(label_base);
  
  if ( nnode != 13) {
    stamped_printf("\tParsing error : 13 params should be found in 'LocalBoundary > HeatTransfer_SN'\n");
    Exit(0);
  }
  
  
  // 表面温度
  label = label_base + "/Surface_Temperature";
  
  REAL_TYPE st = get_BCval_real(label);
  compo[n].set_Temp( FBUtility::convTemp2K(st, Unit_Temp) );
  
  // 面指定
  set_Deface(label_base, n);
  
  // type
  label = label_base + "/Ref_Temp_Mode";
  
  if ( !(tpCntl->GetValue(label, &str )) ) {
    stamped_printf("\tParsing error : Invalid int value for 'Ref_Temp_Mode' in 'LocalBoundary > HeatTransfer_SN'\n");
    Exit(0);
  }
  if ( !strcasecmp("bulk_temperature", str.c_str()) ) {
	  compo[n].set_sw_HTmodeRef( CompoList::HT_mode_bulk );
  }
  else if ( !strcasecmp("local_temperature", str.c_str()) ) {
	  compo[n].set_sw_HTmodeRef( CompoList::HT_mode_local );
  }
  else {
	  stamped_printf("\tParsing error : Invalid string value for 'Ref_Temp_Mode' : %s\n", str.c_str());
	  Exit(0);
  }
  
  // Vertical and upper face values
  label = label_base + "/vertical_laminar_alpha";
  compo[n].ca[CompoList::vert_laminar_alpha]    = get_BCval_real(label);
  
  label = label_base + "/vertical_laminar_beta";
  compo[n].ca[CompoList::vert_laminar_beta]     = get_BCval_real(label);
  
  label = label_base + "/vertical_turbulent_alpha";
  compo[n].ca[CompoList::vert_turbulent_alpha]  = get_BCval_real(label);
  
  label = label_base + "/vertical_turbulent_beta";
  compo[n].ca[CompoList::vert_turbulent_beta]   = get_BCval_real(label);
  
  label = label_base + "/vertical_Ra_critial";
  compo[n].ca[CompoList::vert_Ra_critial]       = get_BCval_real(label);
  
  // Lower face values
  label = label_base + "/lower_laminar_alpha";
  compo[n].cb[CompoList::lower_laminar_alpha]   = get_BCval_real(label);
  
  label = label_base + "/lower_laminar_beta";
  compo[n].cb[CompoList::lower_laminar_beta]    = get_BCval_real(label);
  
  label = label_base + "/lower_turbulent_alpha";
  compo[n].cb[CompoList::lower_turbulent_alpha] = get_BCval_real(label);
  
  label = label_base + "/lower_turbulent_beta";
  compo[n].cb[CompoList::lower_turbulent_beta]  = get_BCval_real(label);
  
  label = label_base + "/lower_Ra_critial";
  compo[n].cb[CompoList::lower_Ra_critial]      = get_BCval_real(label);
  
  if ( Unit_Param != DIMENSIONAL ) {
    stamped_printf("\tWarning: Heat condition must be a dimensional value\n");
    Exit(0);
  }
}



/**
 @fn void ParseBC::get_IBC_HT_SF(const std::string label_base, const int n)
 @brief HeatTransfer_SFのパラメータを取得する
 @param label_base HeatTransfer_SFのレベル
 @param n コンポーネントリストのエントリ番号
 */
void ParseBC::get_IBC_HT_SF(const std::string label_base, const int n)
{
  std::string str;
  string label;
  int nnode=0;
  
  // check number of Elem
  nnode=tpCntl->countLabels(label_base);
  if ( nnode != 6) {
    stamped_printf("\tParsing error : 6 params should be found in 'LocalBoundary > HeatTransfer_SF'\n");
    Exit(0);
  }
  
  // 表面温度
  label=label_base+"/Surface_Temperature";
  REAL_TYPE st = get_BCval_real(label);
  compo[n].set_Temp( FBUtility::convTemp2K(st, Unit_Temp) );
  
  // 面指定
  set_Deface(label_base, n);
  
  // type
  label=label_base+"/Ref_Temp_Mode";//
  if ( !(tpCntl->GetValue(label, &str )) ) {
    stamped_printf("\tParsing error : Invalid int value for 'Ref_Temp_Mode' in 'LocalBoundary > HeatTransfer_SF'\n");
    Exit(0);
  }
  if ( !strcasecmp("bulk_temperature", str.c_str()) ) {
	  compo[n].set_sw_HTmodeRef( CompoList::HT_mode_bulk );
  }
  else if ( !strcasecmp("local_temperature", str.c_str()) ) {
	  compo[n].set_sw_HTmodeRef( CompoList::HT_mode_local );
  }
  else {
	  stamped_printf("\tParsing error : Invalid string value for 'type' : %s\n", str.c_str());
	  Exit(0);
  }
  
  // coefficients
  label=label_base+"/alpha";
  compo[n].ca[CompoList::alpha] = get_BCval_real(label);
  label=label_base+"/beta";
  compo[n].ca[CompoList::beta]  = get_BCval_real(label);
  label=label_base+"/gamma";
  compo[n].ca[CompoList::gamma] = get_BCval_real(label);
  
  if ( Unit_Param != DIMENSIONAL ) {
    stamped_printf("\tWarning: Heat condition must be a dimensional value\n");
    Exit(0);
  }
}

/**
 @fn void ParseBC::get_IBC_HT_B(const std::string label_base, const int n)
 @brief HeatTransfer_Bのパラメータを取得する
 @param label_base HeatTransfer_Bのレベル
 @param n コンポーネントリストのエントリ番号
 */
void ParseBC::get_IBC_HT_B(const std::string label_base, const int n)
{
  string label;
  int nnode=0;
  
  // check number of Elem
  nnode=tpCntl->countLabels(label_base);
  if ( nnode != 3) {
    stamped_printf("\tParsing error : 3 params should be found in 'LocalBoundary > HeatTransfer_B'\n");
    Exit(0);
  }
  
  // 面指定
  set_Deface(label_base, n);
  
  // 熱伝達係数
  label=label_base+"/Coef_of_Heat_Transfer";
  compo[n].set_CoefHT( get_BCval_real(label) );
  
  // バルク温度
  label=label_base+"/Bulk_Temperature";
  REAL_TYPE st = get_BCval_real(label);
  compo[n].set_Temp( FBUtility::convTemp2K(st, Unit_Temp) );
  
  if ( Unit_Param != DIMENSIONAL ) {
    stamped_printf("\tWarning: Heat condition must be given by dimensional value\n");
    Exit(0);
  }
}




/**
 @fn void ParseBC::get_IBC_IsoTherm(const std::string label_base, const int n)
 @brief XMLファイルから境界条件IsoThermalのパラメータを取得し保持する
 @param label_base IsoThermalのレベル
 @param n コンポーネントリストのエントリ番号
 */
void ParseBC::get_IBC_IsoTherm(const std::string label_base, const int n)
{
  string label;
  int nnode=0;
  
  // check number of Elem
  nnode=tpCntl->countLabels(label_base);
  if ( nnode != 2) {
    stamped_printf("\tParsing error : 2 params should be found in 'LocalBoundary > IsoThermal'\n");
    Exit(0);
  }
  
  // 表面温度
  label=label_base+"/temperature";
  REAL_TYPE tmp = get_BCval_real(label);
  compo[n].set_Temp( FBUtility::convTemp2K(tmp, Unit_Temp) );
	
  // 面指定
  set_Deface(label_base, n);
  
  if ( Unit_Param != DIMENSIONAL ) {
    stamped_printf("\tWarning: Heat condition must be given by dimensional value\n");
    Exit(0);
  }
}

/**
 @fn void ParseBC::get_IBC_Radiant(const std::string label_base, const int n)
 @brief XMLファイルから境界条件Radiantのパラメータを取得し保持する
 @todo
 - 境界条件自体は未実装
 */
void ParseBC::get_IBC_Radiant(const std::string label_base, const int n)
{
  string label;
  int nnode=0;
  
  // check number of Elem
  nnode=tpCntl->countLabels(label_base);
  if ( nnode == 0) {
    stamped_printf("\tParsing error : Missing Elements in 'LocalBoundary > Radiant'\n");
    Exit(0);
  }
  
  // 面指定
  set_Deface(label_base, n);
  
  // 係数
  label=label_base+"/epsilon";
  compo[n].set_CoefRadEps( get_BCval_real(label) );
  
  // 射出率
  label=label_base+"/projection";
  compo[n].set_CoefRadPrj( get_BCval_real(label) );
  
  if ( Unit_Param != DIMENSIONAL ) {
    stamped_printf("\tWarning: Heat condition must be a dimensional value\n");
    Exit(0);
  }
}



/**
 @fn void ParseBC::get_IBC_HeatSrc(const std::string label_base, const int n)
 @brief Heat_Generationのパラメータを取得する
 @param label_base Pointer of Configuration tree
 @param n コンポーネントリストのエントリ番号
 @note
 - D1には発熱量を保持，D2には発熱密度を保持
 */
void ParseBC::get_IBC_HeatSrc(const std::string label_base, const int n)
{
  REAL_TYPE hsrc=0.0f;
  std::string str;
  string label;
  int nnode=0;
  
  // check number of Elem
  nnode=tpCntl->countLabels(label_base);
  if ( nnode != 2) {
    stamped_printf("\tParsing error : 2 params should be found in 'Heat_Source'\n");
    Exit(0);
  }
  
  // type
  label=label_base+"/type";//
  if ( !(tpCntl->GetValue(label, &str )) ) {
    stamped_printf("\tParsing error : Invalid int value for 'type' in 'LocalBoundary > Heat_Source'\n");
    Exit(0);
  }
  if ( !strcasecmp("heat_release_value", str.c_str()) ) {
	  compo[n].set_HSRC_policy(true);
  }
  else if ( !strcasecmp("heat_generation_density", str.c_str()) ) {
	  compo[n].set_HSRC_policy(false);
  }
  else {
	  stamped_printf("\tParsing error : Invalid string value for 'type' : %s\n", str.c_str());
	  Exit(0);
  }
  
  // 放熱量
  label=label_base+"/value";//
  if ( !(tpCntl->GetValue(label, &hsrc )) ) {
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
 @fn void ParseBC::get_IBC_CnstTemp(const std::string label_base, const int n)
 @brief Const_Temperatureのパラメータを取得する
 @param label_base Const_Temperatureのレベル
 @param n コンポーネントリストのエントリ番号
 */
void ParseBC::get_IBC_CnstTemp(const std::string label_base, const int n)
{
  string label;
  int nnode=0;
  
  // check number of Elem
  nnode=tpCntl->countLabels(label_base);
  if ( nnode != 1) {
    stamped_printf("\tParsing error : 1 param should be found in 'Specified_Temperature'\n");
    Exit(0);
  }
  
  // 温度
  label=label_base+"/temperature";
  REAL_TYPE tmp = get_BCval_real(label);
  compo[n].set_Temp( FBUtility::convTemp2K(tmp, Unit_Temp) );
  
  if ( Unit_Param != DIMENSIONAL ) {
    stamped_printf("\tWarning: Heat condition must be given by dimensional value\n");
    Exit(0);
  }
}

/**
 @fn void ParseBC::get_IBC_Outflow(const std::string label_base, const int n)
 @brief 内部の流出境界のパラメータを取得する
 @param label_base レベル
 @param n コンポーネントリストのエントリ番号
 */
void ParseBC::get_IBC_Outflow(const std::string label_base, const int n)
{
  int def;
  REAL_TYPE ct;
  REAL_TYPE v[3];
  std::string str;
  std::string label;
  
  // 圧力境界のタイプ default
  compo[n].set_P_BCtype( P_GRAD_ZERO );
  
  // Hidden parameter
  label = label_base + "/pressure_type";

  if ( !(tpCntl->GetValue(label, &str )) ) {
    ; // なくてもOK
  }
  else {
    if ( !strcasecmp("dirichlet", str.c_str()) ) {
      compo[n].set_P_BCtype( P_DIRICHLET );
    }
    else if ( !strcasecmp("grad_zero", str.c_str()) ) {
      compo[n].set_P_BCtype( P_GRAD_ZERO );
    }
    else {
      printf("\tParsing error : Invalid string value for 'Pressure_Type' : %s\n", str.c_str());
      Exit(0);
    }
    if ( compo[n].get_P_BCtype() == P_DIRICHLET ) {
      label = label_base + "/pressure_value";
      compo[n].set_Pressure( get_BCval_real(label) );
    }
  }
  
  // 面指定
  set_Deface(label_base, n);
  
  // 流出速度のタイプ
  compo[n].setOutflowType(V_AVERAGE);
  
  // 法線ベクトルの取得
  get_NV(label_base, n, v);
  copyVec(compo[n].nv, v);
  
  /*
   label=label_base+"/velocity_type";//
   if ( !(tpCntl->GetValue(label, &str )) ) {
   printf("\tParsing error : fail to get 'Velocity_Type' in 'LocalBoundary > Outflow'\n");
   Exit(0);
   }
   if ( !strcasecmp("average", str.c_str()) ) {
   compo[n].flag = V_AVERAGE;
   }
   else if ( !strcasecmp("minmax", str.c_str()) ) {
   compo[n].flag = V_MINMAX;
   }
   else {
   printf("\tParsing error : Invalid string value for 'Velocity_Type' : %s\n", str);
   Exit(0);
   }
   */
  
}



/**
 @fn void ParseBC::get_IBC_SpecVel(const string label_base, const int n)
 @brief 内部の流入境界のパラメータを取得する
 @param label_base レベル
 @param n コンポーネントリストのエントリ番号
 @note Control::setparameters()でcompo[].ca[]に値をセットする
 */
void ParseBC::get_IBC_SpecVel(const std::string label_base, const int n)
{
  std::string str;
  std::string label;
  REAL_TYPE ct;
  REAL_TYPE v[3];
  
  // 指定タイプの特定
  label = label_base + "/type";
  
  if ( !(tpCntl->GetValue(label, &str )) ) {
    stamped_printf("\tParsing error : Invalid Specified_Type in 'LocalBoundary > Specified_Type'\n");
    Exit(0);
  }
  if ( !strcasecmp("velocity", str.c_str()) ) {
	  compo[n].set_VBC_policy(true);
  }
  else if ( !strcasecmp("massflow", str.c_str()) ) {
	  compo[n].set_VBC_policy(false);
  }
  else {
	  printf("\tParsing error : Invalid string value '%s' for 'Type'\n", str.c_str());
	  Exit(0);
  }
  
  // 速度指定タイプ
  compo[n].set_V_profile( get_Vel_profile(label_base) );
  
  
  // 速度パラメータの読み込み
  get_Vel_Params(label_base, compo[n].get_V_Profile(), compo[n].ca, "LocalBoundary", compo[n].isPolicy_Massflow());

  
  // 法線ベクトル
  get_NV(label_base, n, v);
  copyVec(compo[n].nv, v);
  
  
  // heat problem
  if ( HeatProblem ) {
    label = label_base + "/temperature";
    
    if ( !(tpCntl->GetValue(label, &ct )) ) {
      stamped_printf("\tParsing error : fail to get 'Temperature' in 'LocalBoundary > Specified_Velocity'\n");
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
    label = label_base + "/BC_direction";
    
    if ( !(tpCntl->GetValue(label, &str )) ) {
      stamped_printf("\tParsing error : Invalid Specified_Type in 'LocalBoundary > BC_direction'\n");
      Exit(0);
    }
    if ( !strcasecmp("same", str.c_str()) ) {
      compo[n].setBClocation(CompoList::same_direction);
    }
    else if ( !strcasecmp("opposite", str.c_str()) ) {
      compo[n].setBClocation(CompoList::opposite_direction);
    }
    else {
      printf("\tParsing error : Invalid string value '%s' for 'BC_direction'\n", str.c_str());
      Exit(0);
    }
  }
  
}



/**
 @fn void ParseBC::get_Medium_InitTemp(void)
 @brief 温度計算の場合の各媒質の初期値を取得する
 */
void ParseBC::get_Medium_InitTemp(void)
{
  int id;
  const char* p=NULL;
  unsigned Cell_state;
  REAL_TYPE ct;
  
  std::string label, label_base;
  std::string str;
  int counter=0;
  int nnode=0;
  
  label_base="/Parameter/Init_Temp_of_Medium";
  if ( !tpCntl->chkNode(label_base) ) {
    stamped_printf("\tParsing error : Missing the section of 'Init_Temp_of_Medium'\n");
    Exit(0);
  }
  
  unsigned m_no_medium = NoCompo - NoBC;
  
  // check number of Elem
  nnode=tpCntl->countLabels(label_base);
  if ( nnode != m_no_medium ) {
    stamped_printf("\tParsing error : at Init_Temp_of_Medium\n");
    Exit(0);
  }
  //std::cout <<  "nnode : " << nnode << std::endl;
  //std::cout <<  "m_no_medium : " << m_no_medium << std::endl;
  
  
  // load statement list
  for (unsigned i=1; i<=m_no_medium; i++) {
    
    if(!tpCntl->GetNodeStr(label_base,i,&str)){
      stamped_printf("\tGetNodeStr error\n");
      Exit(0);
    }
    
    // state
    label=label_base+"/"+str;
    if ( !(tpCntl->GetValue(label, &str )) ) {
      stamped_printf("\tParsing error : No valid keyword [SOLID/FLUID] in 'Init_Temp_of_Medium'\n");
      Exit(0);
    }
    if ( !strcasecmp("Solid", str.c_str()) ) Cell_state = SOLID;
    else if ( !strcasecmp("Fluid", str.c_str()) ) Cell_state = FLUID;
    else {
      stamped_printf("\tParsing error : No valid keyword [SOLID/FLUID] in 'Init_Temp_of_Medium'\n");
      Exit(0);
    }
    
    //// ID
    //if ( !param->isSetID() ) {
    //  stamped_printf("\tParsing error : No ID for statement in 'Init_Temp_of_Medium'\n");
    //  Exit(0);
    //}
    //if ( -1 == (id=param->GetID()) ) {
    //  stamped_printf("\tParsing error : No valid ID for statement in 'Init_Temp_of_Medium'\n");
    //  Exit(0);
    //}
    
    //// マッチング
    //bool m_flag = false;
    //for (unsigned m=NoBC+1; m<=NoCompo; m++) {
    //  
    //  if ( compo[m].getID() == (unsigned)id ) {
    //    m_flag = true;
    //    if ( compo[m].getState() != Cell_state ) {
    //      stamped_printf("\tError : Inconsistent the cell state between 'Model_Setting' and 'Init_Temp_of_Medium' : Medium ID=%d\n", id);
    //      Exit(0);
    //    }
    //    
    //    if ( !param->GetData(&ct) ) {
    //      stamped_printf("\tParsing error : No initial temperature in 'Init_Temp_of_Medium'\n");
    //      Exit(0);
    //    }
    //    compo[m].setInitTemp( FBUtility::convTemp2K(ct, Unit_Temp) );
    //    break;
    //  }
    //}
    
    ////check
    //if ( !m_flag ) {
    //  stamped_printf("\tError : could not find ID=%d in ComponentList\n", id);
    //  Exit(0);
    //}
    
  }
}



/**
 @fn void ParseBC::get_NV(const std::string label_base, const int n, REAL_TYPE* v)
 @brief 内部境界条件の法線ベクトル値を取得し，登録する
 @param label_base
 @param n オーダー
 @param v[out] ベクトルパラメータ
 @param str エラー表示用文字列
 */
void ParseBC::get_NV(const std::string label_base, const int n, REAL_TYPE* v)
{
  std::string label;
  for (unsigned i=0; i<3; i++) v[i]=0.0f;
  
  label = label_base + "/Normal";
  if( !(tpCntl->GetVector(label, v, 3)) )
  {
    stamped_printf("\tParsing error : fail to get vec params in '%s\n", label.c_str());
    Exit(0);
  }
  
  //単位ベクトル化
  getUnitVec(v);
}



/**
 @fn void ParseBC::get_OBC_Wall(const std::string label_base, const int n)
 @brief 外部境界の壁面条件のパラメータを取得する
 @param label_base 
 @param n 面番号
 */
void ParseBC::get_OBC_Wall(const std::string label_base, const int n)
{
  REAL_TYPE vel, ct;
  REAL_TYPE v[3];
  std::string str;
  std::string label, label2;
  
  // 速度のタイプの特定
  label = label_base + "/type";
  
  if ( !(tpCntl->GetValue(label, &str )) ) {
    stamped_printf("\tParsing error : Invalid Specified_Type in 'Basic_BCs > Type'\n");
    Exit(0);
  }
  
  if ( !strcasecmp("fixed", str.c_str()) ) {
	  BaseBc[n].set_Type(BoundaryOuter::fixed);
    BaseBc[n].set_V_Profile( CompoList::vel_zero );
  }
  else if ( !strcasecmp("slide", str.c_str()) ) {
	  BaseBc[n].set_Type(BoundaryOuter::slide);
    BaseBc[n].set_V_Profile( get_Vel_profile(label_base) );
  }
  else {
	  printf("\tParsing error : Invalid string value '%s' for 'Type'\n", str.c_str());
	  Exit(0);
  }

   // 法線ベクトル
  if ( BaseBc[n].get_V_Profile() != CompoList::vel_zero ) {
    get_NV(label_base, n, v);
    BaseBc[n].addVec(v); 
  }
  
  // 速度のパラメータ読み込み
  get_Vel_Params(label_base, BaseBc[n].get_V_Profile(), BaseBc[n].ca, "OuterBoundary");
  
  
  
  // heat problem
  if ( HeatProblem ) {
    
    label = label_base + "/heat_type";
    if ( !(tpCntl->GetValue(label, &str )) ) {
      stamped_printf("\tParsing error : fail to get 'Heat_Type' in 'Basic_BCs > wall'\n");
      Exit(0);
    }
    if ( Unit_Param != DIMENSIONAL ) {
      stamped_printf("\tError: Heat condition must be given by dimensional value\n");
      Exit(0);
    }
    
    if     ( !strcasecmp(str.c_str(), "Adiabatic") ){
      BaseBc[n].set_hType(ADIABATIC);
      BaseBc[n].p = 0.0;
    }
    else if( !strcasecmp(str.substr(0,12).c_str(), "HeatTransfer") ) {
      BaseBc[n].set_hType(TRANSFER);
      get_OBC_HT(label, n, str);
    }
    else if( !strcasecmp(str.c_str(), "HeatFlux") )     {
      BaseBc[n].set_hType(HEATFLUX);
      label2 = label + "/Heat_Flux";
      BaseBc[n].set_Heatflux( get_BCval_real(label2) ); // 正符号は流入
    }
    else if( !strcasecmp(str.c_str(), "Isothermal") )   {
      BaseBc[n].set_hType(ISOTHERMAL);
      label2 = label + "/temperature";
      ct = get_BCval_real(label2); // 表面温度
      BaseBc[n].set_Temp( FBUtility::convTemp2K(ct, Unit_Temp) );
    }
    else if( !strcasecmp(str.c_str(), "Constant_Temperature") )   {
      BaseBc[n].set_hType(CNST_TEMP);
      label = label_base + "/temperature";
      ct = get_BCval_real(label); // 指定温度
      BaseBc[n].set_Temp( FBUtility::convTemp2K(ct, Unit_Temp) );
    }
    else {
      stamped_printf("\tParsing error : Invalid string value for 'Heat_type' : %s\n", str.c_str());
      Exit(0);
    }    
  }
}



/**
 @fn void ParseBC::get_OBC_Outflow(const string label_base, const int n)
 @brief 外部境界の流出条件のパラメータを取得する
 @param label_base 
 @param n 面番号
 @note 圧力の値は，Control::setParameters()で無次元化する
 */
void ParseBC::get_OBC_Outflow(const std::string label_base, const int n)
{
  std::string str;
  std::string label;
  REAL_TYPE ct;
  
  
  // 流出速度のタイプ
  label = label_base + "/vel_type";
  
  if ( !(tpCntl->GetValue(label, &str )) ) {
    stamped_printf("\tParsing error : fail to get 'Vel_Type' in 'Basic_BCs > outflow'\n");
    Exit(0);
  }
  if ( !strcasecmp("average", str.c_str()) ) {
	  BaseBc[n].set_ofv(V_AVERAGE);
  }
  else if ( !strcasecmp("minmax", str.c_str()) ) {
	  BaseBc[n].set_ofv(V_MINMAX);
  }
  else {
	  stamped_printf("\tParsing error : Invalid string value for 'Vel_Type' : %s\n", str.c_str());
	  Exit(0);
  }
  
  
  // 圧力境界のタイプ  default
  BaseBc[n].set_pType(P_GRAD_ZERO);
  BaseBc[n].p = 0.0; // ダミー値
  
  // Hidden option
  label = label_base + "/prs_type";
  
  if ( !(tpCntl->GetValue(label, &str )) ) {
  }
  else {
    if ( !strcasecmp("dirichlet", str.c_str()) ) {
      BaseBc[n].set_pType(P_DIRICHLET);
    }
    else if ( !strcasecmp("grad_zero", str.c_str()) ) {
      BaseBc[n].set_pType(P_GRAD_ZERO);
    }
    else {
      stamped_printf("\tParsing error : Invalid string value for 'Pressure_Type' : %s\n", str.c_str());
      Exit(0);
    }
  }
  
  
  // 圧力の値
  if ( BaseBc[n].get_pType() == P_DIRICHLET ) {
    if ( !strcasecmp("dirichlet", str.c_str()) ) {
      label = label_base + "/prs_value";
      if ( !(tpCntl->GetValue(label, &ct )) ) {
        stamped_printf("\tParsing error : fail to get 'prs_value' in 'Basic_BCs > outflow'\n");
        Exit(0);
      }
      else {
        BaseBc[n].p = ct;
      }
    }
  }
}


/**
 @fn void ParseBC::get_OBC_SpecVH(const std::string label_base, const int n)
 @brief 外部境界の流入条件のパラメータを取得する
 @param label_base 
 @param n 面番号
 */
void ParseBC::get_OBC_SpecVH(const std::string label_base, const int n)
{
  REAL_TYPE vel, ct;
  REAL_TYPE v[3];
  std::string str;
  std::string label;
  
  
  // 速度境界条件のプロファイル
  BaseBc[n].set_V_Profile( get_Vel_profile(label_base) );
  
  // 法線ベクトル
  if ( BaseBc[n].get_V_Profile() != CompoList::vel_zero ) {
    get_NV(label_base, n, v);
    BaseBc[n].addVec(v);
  }
  
  // 速度のパラメータ読み込み
  get_Vel_Params(label_base, BaseBc[n].get_V_Profile(), BaseBc[n].ca, "OuterBoundary");

  
  // heat problem
  if ( HeatProblem ) {
    if ( Unit_Param != DIMENSIONAL ) {
      stamped_printf("\tError: Heat condition must be given by dimensional value\n");
      Exit(0);
    }
    label = label_base + "/temperature";
    if ( !(tpCntl->GetValue(label, &ct )) ) {
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
 @fn void ParseBC::get_OBC_Trcfree(const std::string label_base, const int n)
 @brief 外部境界の流入条件のパラメータを取得する
 @param label_base 
 @param n 面番号
 */
void ParseBC::get_OBC_Trcfree(const std::string label_base, const int n)
{
  REAL_TYPE ct;
  std::string str;
  std::string label;
  
  BaseBc[n].set_pType(P_DIRICHLET);
  BaseBc[n].p = 0.0; // ゲージ圧zero 固定
  
  // 外部雰囲気温
  if ( HeatProblem ) {
    label=label_base+"/ambient_temperature";
    if ( !(tpCntl->GetValue(label, &ct )) ) {
      stamped_printf("\tParsing error : fail to get 'ambient_temperature' in 'Basic_BCs > IN_OUT'\n");
      Exit(0);
    }
    else {
      BaseBc[n].set_Temp( FBUtility::convTemp2K(ct, Unit_Temp) );
    }
    if ( Unit_Param != DIMENSIONAL ) {
      stamped_printf("\tError: Heat condition must be given by dimensional value\n");
      Exit(0);
    }
  }
}


/**
 @fn void ParseBC::get_OBC_FarField(const std::string label_base, const int n)
 @brief 外部境界の遠方境界のパラメータを取得する
 @param label_base 
 @param n 面番号
 */
void ParseBC::get_OBC_FarField(const std::string label_base, const int n)
{
  REAL_TYPE ct;
  std::string str;
  std::string label;
  
  BaseBc[n].set_pType(P_GRAD_ZERO);
  BaseBc[n].p = 0.0; // ダミー値
  
  // 外部雰囲気温
  if ( HeatProblem ) {
    label=label_base+"/ambient_temperature";
    if ( !(tpCntl->GetValue(label, &ct )) ) {
      stamped_printf("\tParsing error : fail to get 'ambient_temperature' in 'Basic_BCs > IN_OUT'\n");
      Exit(0);
    }
    else {
      BaseBc[n].set_Temp( FBUtility::convTemp2K(ct, Unit_Temp) );
    }
    if ( Unit_Param != DIMENSIONAL ) {
      stamped_printf("\tError: Heat condition must be given by dimensional value\n");
      Exit(0);
    }
  }
  
}



/**
 @fn void ParseBC::get_OBC_Periodic(const std::string label_base, const int n)
 @brief 外部境界の周期条件のパラメータを取得する
 @param label_base 
 @param n 面番号
 @note 圧力の値は，Control::setParameters()で無次元化する
 */
void ParseBC::get_OBC_Periodic(const std::string label_base, const int n)
{
  REAL_TYPE ct;
  int def;
  std::string str;
  std::string label;
  
  // モード
  label = label_base + "/mode";
  if ( !(tpCntl->GetValue(label, &str )) ) {
    printf("\tParsing error : No 'mode' section in 'Basic_BCs > periodic'\n");
    Exit(0);
  }
  else {
    if ( !strcasecmp(str.c_str(), "simple_copy") ) {
      BaseBc[n].set_PrdcMode(BoundaryOuter::prdc_Simple);
    }
    else if ( !strcasecmp(str.c_str(), "directional") ) {
      BaseBc[n].set_PrdcMode(BoundaryOuter::prdc_Directional);
    }
    else if ( !strcasecmp(str.c_str(), "driver") ) {
      BaseBc[n].set_PrdcMode(BoundaryOuter::prdc_Driver);
    }
  }
  
  // Directional
  if ( BaseBc[n].get_PrdcMode() == BoundaryOuter::prdc_Directional ) {
    label = label_base + "/pressure_difference";
    if ( !(tpCntl->GetValue(label, &ct )) ) {
      printf("\tParsing error : No 'Pressure_Difference' keyword in 'Basic_BCs > Periodic'\n");
      Exit(0);
    }
    else {
      BaseBc[n].p = ct;
    }
    
    label = label_base + "/flow_direction";
    if ( !(tpCntl->GetValue(label, &str )) ) {
      printf("\tParsing error : No 'Flow_Direction' keyword in 'Basic_BCs > Periodic'\n");
      Exit(0);
    }
    else{
      if ( !strcasecmp(str.c_str(), "upstream") ) {
        BaseBc[n].set_FaceMode(BoundaryOuter::prdc_upstream);
      }
      else if ( !strcasecmp(str.c_str(), "downstream") ) {
        BaseBc[n].set_FaceMode(BoundaryOuter::prdc_downstream);
      }
      else {
        printf("\tParsing error : Invalid keyword in 'Basic_BCs > Periodic > Flow_Direction'\n");
        Exit(0);
      }
    }
  }
  
  // Driver
  if ( BaseBc[n].get_PrdcMode() == BoundaryOuter::prdc_Driver ) {
    
    label = label_base + "/driver_direction";
    if ( !(tpCntl->GetValue(label, &str )) ) {
      printf("\tParsing error : No 'Driver_Direction' keyword in 'Basic_BCs > Periodic'\n");
      Exit(0);
    }
    else {
      if ( !strcasecmp("x_minus", str.c_str()) ) {
        def = X_MINUS;
      }
      else if ( !strcasecmp("x_plus", str.c_str()) ) {
        def = X_PLUS;
      }
      else if ( !strcasecmp("y_minus", str.c_str()) ) {
        def = Y_MINUS;
      }
      else if ( !strcasecmp("y_plus", str.c_str()) ) {
        def = Y_PLUS;
      }
      else if ( !strcasecmp("z_minus", str.c_str()) ) {
        def = Z_MINUS;
      }
      else if ( !strcasecmp("z_plus", str.c_str()) ) {
        def = Z_PLUS;
      }
      else {
        printf("\tParsing error : Invalid keyword in 'Basic_BCs > Periodic > Driver_Direction'\n");
        Exit(0);
      }
      BaseBc[n].set_DriverDir(def);
    }
    
    
    label = label_base + "/driver_lid_index";
    if ( !(tpCntl->GetValue(label, &def )) ) {
      printf("\tParsing error : Invalid 'Driver_Lid_Index' keyword in 'Basic_BCs > Periodic'\n");
      Exit(0);
    }
    else {
      BaseBc[n].set_DriverIndex(def);
    }
  }
}


/*
 @fn void ParseBC::get_OBC_HT(const std::string label_base, const int n, const std::string kind)
 @brief 外部の壁面熱伝達境界のパラメータを取得する
 @param label_base 
 @param n 面番号
 @param kind 熱伝達境界の種類
 */ 
void ParseBC::get_OBC_HT(const std::string label_base, const int n, const std::string kind)
{
  std::string label;
  std::string str;
  REAL_TYPE ct;
  
  if ( !strcasecmp(kind.c_str(), "HeatTransfer_B") ) {
    BaseBc[n].set_HTmode(HT_B);
    label = label_base + "/Coef_of_Heat_Transfer";
    BaseBc[n].set_CoefHT( get_BCval_real(label) );
    
    label = label_base + "/Bulk_Temperature";
    ct = get_BCval_real(label);
    BaseBc[n].set_Temp( FBUtility::convTemp2K(ct, Unit_Temp) );
  }
  else if ( !strcasecmp(kind.c_str(), "HeatTransfer_N") ) {
    BaseBc[n].set_HTmode(HT_N);
    label = label_base + "/Coef_of_Heat_Transfer";
    BaseBc[n].set_CoefHT( get_BCval_real(label) );
  }
  else if ( !strcasecmp(kind.c_str(), "HeatTransfer_S") ) {
    BaseBc[n].set_HTmode(HT_S);
    label = label_base + "/Coef_of_Heat_Transfer";
    BaseBc[n].set_CoefHT( get_BCval_real(label) );
    
    label = label_base + "/Surface_temperature";
    ct = get_BCval_real(label);
    BaseBc[n].set_Temp( FBUtility::convTemp2K(ct, Unit_Temp) );
  }
  else if ( !strcasecmp(kind.c_str(), "HeatTransfer_SF") ) {
    BaseBc[n].set_HTmode(HT_SF);
    label = label_base + "/Surface_temperature";
    BaseBc[n].set_Temp( get_BCval_real(label) );
    
    
    label=label_base+"/ref_temp_mode";
    if ( !(tpCntl->GetValue(label, &str )) ) {
      stamped_printf("\tParsing error : Invalid int value for 'ref_temp_mode' in 'Basic_BCs > wall'\n");
      Exit(0);
    }
    if ( !strcasecmp("bulk_temperature", str.c_str()) ) {
      BaseBc[n].set_HTmodeRef( CompoList::HT_mode_bulk );
    }
    else if ( !strcasecmp("local_temperature", str.c_str()) ) {
      BaseBc[n].set_HTmodeRef( CompoList::HT_mode_local );
    }
    else {
      stamped_printf("\tParsing error : Invalid string value for 'ref_temp_mode' : %s\n", str.c_str());
      Exit(0);
    }
    // coefficients
    label = label_base + "/alpha";
    BaseBc[n].ca[0] = get_BCval_real(label);
    
    label = label_base + "/beta";
    BaseBc[n].ca[1] = get_BCval_real(label);
    
    label = label_base + "/gamma";
    BaseBc[n].ca[2] = get_BCval_real(label);
  }
  else if ( !strcasecmp(kind.c_str(), "HeatTransfer_SN") ) {
    BaseBc[n].set_HTmode(HT_SN);
    label = label_base + "/Surface_temperature";
    BaseBc[n].set_Temp( get_BCval_real(label) );
    
    // reference mode
    label = label_base + "/ref_temp_mode";
    if ( !(tpCntl->GetValue(label, &str )) ) {
      stamped_printf("\tParsing error : Invalid int value for 'ref_temp_mode' in 'Basic_BCs > wall'\n");
      Exit(0);
    }
    if ( !strcasecmp("bulk_temperature", str.c_str()) ) {
      BaseBc[n].set_HTmodeRef( CompoList::HT_mode_bulk );
    }
    else if ( !strcasecmp("local_temperature", str.c_str()) ) {
      BaseBc[n].set_HTmodeRef( CompoList::HT_mode_local );
    }
    else {
      stamped_printf("\tParsing error : Invalid string value for 'ref_temp_mode' : %s\n", str.c_str());
      Exit(0);
    }
    // Vertical and upper face values
    label = label_base + "/vertival_laminar_alpha";
    BaseBc[n].ca[0] = get_BCval_real(label);
    
    label = label_base + "/vertival_laminar_beta";
    BaseBc[n].ca[1] = get_BCval_real(label);
    
    label = label_base + "/vertival_turbulent_alpha";
    BaseBc[n].ca[2] = get_BCval_real(label);
    
    label = label_base + "/vertival_turbulent_beta";
    BaseBc[n].ca[3] = get_BCval_real(label);
    
    label = label_base + "/vertival_Ra_critial";
    BaseBc[n].ca[4] = get_BCval_real(label);
    
    // Lower face values
    label = label_base + "/lower_laminar_alpha";
    BaseBc[n].cb[0] = get_BCval_real(label);
    
    label = label_base + "/lower_laminar_beta";
    BaseBc[n].cb[1] = get_BCval_real(label);
    
    label = label_base + "/lower_turbulent_alpha";
    BaseBc[n].cb[2] = get_BCval_real(label);
    
    label = label_base + "/lower_turbulent_beta";
    BaseBc[n].cb[3] = get_BCval_real(label);
    
    label = label_base + "/lower_Ra_critial";
    BaseBc[n].cb[4] = get_BCval_real(label);
  }
  else {
    stamped_printf("\tParsing error : fail to get HeatTransfer Type in 'Basic_BCs > wall'\n");
    Exit(0);
  }
}



/**
 @fn void ParseBC::get_Phase(void)
 @brief 2相流問題で気相か液相かを取得する
 */
void ParseBC::get_Phase(void)
{
  unsigned m_phase;
  int id;
  std::string str,p;
  string label,label_base;
  int NoParam;
  
  ///////////////////////////////////////////////////////////////////////////////
  stamped_printf("\tWARNING not yet\n");
  Exit(0);
  ///////////////////////////////////////////////////////////////////////////////
  
  
  label_base="/Steer/Phase_Idetification";
  //std::cout <<  "label : " << label << std::endl;
  if ( !(tpCntl->chkNode(label_base)) ) {
    stamped_printf("\tParsing error : Missing the section of 'Phase_Idetification'\n");
    Exit(0);
  }
  
  // load statement list
  NoParam=tpCntl->countLabels(label_base);
  if ( NoParam < 0) {
    stamped_printf("\tcountLabels --- %s\n",label_base.c_str());
    Exit(0);
  }
  
  for (int i=1; i<=NoParam; i++) {
    if(!tpCntl->GetNodeStr(label_base,i,&str)){
      stamped_printf("\tParsing error : No Leaf Node \n");
      Exit(0);
    }
    label=label_base+"/"+str;//
    
    /*
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
     if ( (0 >= def) || (def > NoMedium) ) {
     stamped_printf("\tParsing error : ID[%d] described in 'id' is not included in 'Model_Setting'\n", id);
     Exit(0);
     }
     
     // set phase of FLUID
     for (int n=1; n<=NoMedium; n++) {
     if ( compo[n].getID() == id ) {
     if ( compo[n].getState() == FLUID ) compo[n].setPhase(m_phase);
     }
     }
     */
    
  }
  
  // check Phase of Fluid
  unsigned tmp;
  bool sw=true;
  for (int n=1; n<=NoMedium; n++) {
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




//@brief 外部境界の速度境界条件のタイプを取得し，返す
int ParseBC::get_Vel_profile(const std::string label_base)
{
  std::string label, str;
  
  label = label_base + "/Profile";
  
  if ( !(tpCntl->GetValue(label, &str )) ) {
    printf("\tParsing error : fail to get 'Profile' in '%s'\n", str.c_str());
    Exit(0);
  }
  if ( !strcasecmp("constant", str.c_str()) ) {
		return CompoList::vel_constant;
  }
  else if ( !strcasecmp("harmonic", str.c_str()) ) {
		return CompoList::vel_harmonic;
  }
  else {
	  printf("\tParsing error : Invalid string value for '%s' : %s\n", label.c_str(), str.c_str());
	  Exit(0);
  }
  return 0;
}



/**
 @fn void ParseBC::get_Vel_Params(const std::string label_base, const int prof, REAL_TYPE* ca, const char* str, const bool policy)
 @brief 速度のパラメータを取得する
 @param label_base 
 @param ca 係数パラメータの配列
 @param policy 流量指定(true) or 速度指定(false) > 外部境界は速度のみ
 @note 
 - 値は，Control::setParameters()で無次元化する
 - 速度プロファイルは単振動と一定値の場合で係数の保持パターンが異なる
 - 内部境界の場合には，流量指定と速度指定があるので分岐処理（未実装）
 */
void ParseBC::get_Vel_Params(const std::string label_base, const int prof, REAL_TYPE* ca, const char* str, const bool policy)
{
  std::string label;
  REAL_TYPE ct=0.0, vel;
  
  if ( prof == CompoList::vel_zero ) {
    ca[CompoList::amplitude] = 0.0;
    ca[CompoList::frequency] = 0.0;
    ca[CompoList::initphase] = 0.0;
    ca[CompoList::bias]      = 0.0;
  }
  else if ( prof == CompoList::vel_constant ) {
    label = label_base + "/velocity";
    
    if ( !(tpCntl->GetValue(label, &ct )) ) {
      stamped_printf("\tParsing error : Invalid value in '%s > Velocity'\n", str);
      Exit(0);
    }
    if ( policy ) { // 流量指定の場合 LocalBoundaryの場合のみ
      vel = ct; // 有次元でも無次元でも，モデルの断面積を計算して，後ほどパラメータ計算 >> Control::setParameters()
    }
    else {
      vel = ( Unit_Param == DIMENSIONAL ) ? ct : ct * RefVelocity;
    }
    
    ca[CompoList::amplitude] = 0.0;
    ca[CompoList::frequency] = 0.0;
    ca[CompoList::initphase] = 0.0;
    ca[CompoList::bias]      = vel;
  }
  else if ( prof == CompoList::vel_harmonic) {
    label = label_base + "/amplitude";
    if ( !(tpCntl->GetValue(label, &ct )) ) {
      stamped_printf("\tParsing error : fail to get 'Amplitude' in '%s'\n", label_base.c_str());
      Exit(0);
    }
    ca[CompoList::amplitude] = vel;
    
    label = label_base + "/frequency";
    if ( !(tpCntl->GetValue(label, &ct )) ) {
      stamped_printf("\tParsing error : fail to get 'Frequency' in '%s'\n", label_base.c_str());
      Exit(0);
    }
    ca[CompoList::frequency] = ct;
    
    label = label_base + "/initial_phase";
    if ( !(tpCntl->GetValue(label, &ct )) ) {
      stamped_printf("\tParsing error : fail to get 'Initial_phase' in '%s'\n", label_base.c_str());
      Exit(0);
    }
    ca[CompoList::initphase] = ct;
    
    
    label = label_base + "/constant_bias";
    if ( !(tpCntl->GetValue(label, &ct )) ) {
      stamped_printf("\tParsing error : fail to get 'Constant_Bias' in '%s'\n", label_base.c_str());
      Exit(0);
    }
    ca[CompoList::bias] = ct;
  }
  else {
    Exit(0);
  }
  
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
 @fn bool ParseBC::isIDinCompo(int candidate, unsigned now)
 @brief candidate_idが既にコンポーネントに登録されているかを調べる
 @retval 重複していればfalseを返す
 @param candidate チェックの候補
 @param now コンポーネントリストの現在までのエントリ番号
 */
bool ParseBC::isIDinCompo(int candidate, unsigned now)
{
  for (int i=1; i<(int)now; i++) {
    if ( candidate == compo[i].getMatOdr() ) {
      return false;
    }
  }
  return true;
}

/**
 @fn bool ParseBC::isIDinCompo(int candidate, int def, unsigned now)
 @brief candidate_idとdefの組が既にコンポーネントに登録されているかを調べる
 @retval 重複していればfalseを返す
 @param candidate_id チェックの候補ID
 @param def deface
 @param now コンポーネントリストの現在までのエントリ番号
 */
bool ParseBC::isIDinCompo(int candidate, int def, unsigned now)
{
  for (int i=1; i<(int)now; i++) {
    if ( (candidate == compo[i].getMatOdr()) && (def == compo[i].getDef()) ) {
      return false;
    }
  }
  return true;
}




/**
 @fn void ParseBC::loadBC_Local(Control* C)
 @brief CompoListに内部境界条件の情報を設定する
 @note
 - 最初にBCの情報を登録，その後IDの情報を登録
 - パラメータファイルから各内部BCのidをパースし，compoに保持する
 - 格納番号は1からスタート
 */
void ParseBC::loadBC_Local(Control* C)
{ 
  std::string str, label, ename;
  std::string label_base, label_ename, label_leaf;
  REAL_TYPE fval;
  int n=0;
  int ide;
  unsigned tp;
  
  /* Medium_Table debug print ///////////////////////////////////////////////////
   std::cout << std::endl;
   std::cout << std::endl;
   cout << "**********************" << endl;
   cout << "*****loadBC_Local*****" << endl;
   cout << "**********************" << endl;
   cout << "NoMedium : " << NoMedium << endl;
   // イテレータを生成
   std::map<string, REAL_TYPE>::iterator itr;
   for (int i=1; i<=NoMedium; i++) { //Medium_Table loop
   cout << "at glance" << i+1 << endl;
   cout << "type  : " << MTITP[i].type << endl;
   cout << "label : " << MTITP[i].label << endl;
   
   int icounter=0;
   for (itr = MTITP[i].m_fval.begin(); itr != MTITP[i].m_fval.end(); itr++)
   {
   icounter++;
   std::string a1 = itr->first;// キー取得
   REAL_TYPE a2 = itr->second;// 値取得
   cout << "i = " << i << "  icounter = " << icounter 
   << "  key:" << a1 << "  value:" << a2 << endl;
   }
   }
   std::cout << " NoBC   = " << NoBC  << std::endl;
   std::cout << std::endl;
   std::cout << std::endl;
   /////////////////////////////////////////////////////////////////////////////// */
  
  // 内部境界条件の有無を調べる
  
  label_base = "/BC_Table/LocalBoundary";
  n = tpCntl->countLabels(label_base);
  
  if ( n != NoBC) {
    stamped_printf("\tLocalBoundary error %s\n", label_base.c_str());
    Exit(0);
  }
  
  //
  //内部境界の条件設定 --- NoBC = 内部境界の数
  //
  
  // 境界条件がなければ，スキップ
  for (int odr=1; odr<=(int)NoBC; odr++) {
    
    if( !tpCntl->GetNodeStr(label_base, odr, &str)){
      stamped_printf("\tParsing error : No Leaf Node \n");
      Exit(0);
    }
    label_ename = label_base + "/" + str;
    ename = str;
    
    // cmp[].type, h_typeのセット ---> setType
    setKeywordLBC(ename, odr);
    
    // nodeの移動
    if ( !tpCntl->GetNodeStr(label_ename, 1, &str) ){
      stamped_printf("\tParsing error : No Leaf Node \n");
      Exit(0);
    }
    label_leaf = label_ename + "/" + str;
    
    // Labelの取得
    compo[odr].setLabel(str);
    
    // 各BCの処理
    tp = compo[odr].getType();
    
    if ( tp == SPEC_VEL ) {
      get_IBC_SpecVel(label_leaf, odr);      
    }
    else if ( tp == OUTFLOW ) {
      get_IBC_Outflow(label_leaf, odr);     
    }
    else if ( tp == IBM_DF ) {
      get_IBC_IBM_DF(label_leaf, odr);
      if ( !isIDinCompo(ide, odr) ) {
        stamped_printf("Parse Error : Reduplication of ID[%d] for BC\n", ide);
        Exit(0);
      }
    }
    else if ( tp == HEX ) {
      get_IBC_PrsLoss(label_leaf, odr);//debug済み
      if ( !isIDinCompo(ide, odr) ) {
        stamped_printf("Parse Error : Reduplication of ID[%d] for BC\n", ide);
        Exit(0);
      }
    }
    else if ( tp == FAN ) {
      get_IBC_Fan(label_leaf, odr);//debug済み
      if ( !isIDinCompo(ide, odr) ) {
        stamped_printf("Parse Error : Reduplication of ID[%d] for BC\n", ide);
        Exit(0);
      }
    }
    else if ( tp == DARCY ) {
      get_Darcy(label_leaf, odr);
      if ( !isIDinCompo(ide, odr) ) {
        stamped_printf("Parse Error : Reduplication of ID[%d] for BC\n", ide);
        Exit(0);
      }
    }
    else if ( tp == CELL_MONITOR ) {
      get_IBC_Monitor(label_leaf, odr, C);//debug済み
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
      get_IBC_Periodic(label_leaf, odr);
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
        get_IBC_Adiabatic(label_leaf, odr);
        if ( !isIDinCompo(ide, compo[odr].getDef(), odr) ) {
          stamped_printf("Parse Error : Reduplication of pair of ID[%d] and Def[%d] for BC\n", ide, compo[odr].getDef());
          Exit(0);
        }        
      }
      else if ( tp == HEATFLUX ) {
        get_IBC_HeatFlux(label_leaf, odr);
        if ( !isIDinCompo(ide, compo[odr].getDef(), odr) ) {
          stamped_printf("Parse Error : Reduplication of pair of ID[%d] and Def[%d] for BC\n", ide, compo[odr].getDef());
          Exit(0);
        }        
      }
      else if ( tp == TRANSFER ) {
        switch ( compo[odr].getHtype() ) {
          case HT_N:
            get_IBC_HT_N(label_leaf, odr);
            if ( !isIDinCompo(ide, compo[odr].getDef(), odr) ) {
              stamped_printf("Parse Error : Reduplication of pair of ID[%d] and Def[%d] for BC\n", ide, compo[odr].getDef());
              Exit(0);
            }
            break;
            
          case HT_S:
            get_IBC_HT_S(label_leaf, odr);
            if ( !isIDinCompo(ide, compo[odr].getDef(), odr) ) {
              stamped_printf("Parse Error : Reduplication of pair of ID[%d] and Def[%d] for BC\n", ide, compo[odr].getDef());
              Exit(0);
            }
            break;
            
          case HT_SN:
            get_IBC_HT_SN(label_leaf, odr);
            if ( !isIDinCompo(ide, compo[odr].getDef(), odr) ) {
              stamped_printf("Parse Error : Reduplication of pair of ID[%d] and Def[%d] for BC\n", ide, compo[odr].getDef());
              Exit(0);
            }
            break;
            
          case HT_SF:
            get_IBC_HT_SF(label_leaf, odr);
            if ( !isIDinCompo(ide, compo[odr].getDef(), odr) ) {
              stamped_printf("Parse Error : Reduplication of pair of ID[%d] and Def[%d] for BC\n", ide, compo[odr].getDef());
              Exit(0);
            }
            break;
            
          case HT_B:
            get_IBC_HT_B(label_leaf, odr);
            if ( !isIDinCompo(ide, compo[odr].getDef(), odr) ) {
              stamped_printf("Parse Error : Reduplication of pair of ID[%d] and Def[%d] for BC\n", ide, compo[odr].getDef());
              Exit(0);
            }
            break;
        }        
      }
      else if ( tp == ISOTHERMAL ) {
        get_IBC_IsoTherm(label_leaf, odr);
        if ( !isIDinCompo(ide, compo[odr].getDef(), odr) ) {
          stamped_printf("Parse Error : Reduplication of pair of ID[%d] and Def[%d] for BC\n", ide, compo[odr].getDef());
          Exit(0);
        }
      }
      else if ( tp == RADIANT ) {
        get_IBC_Radiant(label_leaf, odr);
        if ( !isIDinCompo(ide, odr) ) {
          stamped_printf("Parse Error : Reduplication of ID[%d] for BC\n", ide);
          Exit(0);
        }
      }
      else if ( tp == HEAT_SRC ) {
        get_IBC_HeatSrc(label_leaf, odr);
        if ( !isIDinCompo(ide, odr) ) {
          stamped_printf("Parse Error : Reduplication of ID[%d] for BC\n", ide);
          Exit(0);
        }
      }
      else if ( tp == CNST_TEMP ) {
        get_IBC_CnstTemp(label_leaf, odr);
      }
      else {
        printf("\tError : Invalid Local BC keyword [%d]\n", tp);
        Exit(0);
      }
    }
  }
  
  // 媒質情報の登録
  for (int i=1; i<=NoMedium; i++) {
    compo[NoBC+i].setState( mat[i].getState() );
    compo[NoBC+i].setLabel( mat[i].getLabel() );
    compo[NoBC+i].setMatOdr(i);
  }
}



/**
 @fn void ParseBC::loadBC_Outer(void)
 @brief パラメータファイルをパースして，外部境界条件を取得，保持する
 */
void ParseBC::loadBC_Outer(void)
{
  std::string label_base, label_leaf, label;
  std::string str;
  
  // Basic Outer BCリストの読み込み
  label_base = "/BC_Table/OuterBoundary";
  
  if ( !tpCntl->chkNode(label_base) ) {
    Hostonly_ stamped_printf("\tParsing error : Missing OuterBoundary tree\n");
    Exit(0);
  }
  
  //Basic_BCs
  for (int i=0; i<NoBaseBC; i++) {
    
    if(!tpCntl->GetNodeStr(label_base, i+1, &str)){
      Hostonly_ printf("\tParsing error : Missing 'Basic_BCs'\n");
      Exit(0);
    }
    if( strcasecmp(str.substr(0,9).c_str(), "Basic_BCs") ) continue;
    
    // alias ユニークな名称であること
    label_leaf = label_base + "/" + str;
    label = label_leaf + "/alias";
    
    if ( !(tpCntl->GetValue(label, &str )) ) {
      Hostonly_ printf("\tParsing error : No 'Alias' in 'Basic_BCs'\n");
      Exit(0);
    }
    if ( !chkDuplicate(i, str) ) {
      Hostonly_ printf("\tParsing error : 'Alias' must be unique\n");
      Exit(0);
    }
    BaseBc[i].set_Alias(str);
    
    
    // Classに境界条件の種別をセットする
    label = label_leaf + "/class";
    
    if ( !(tpCntl->GetValue(label, &str )) ) {
      Hostonly_ printf("\tParsing error : No 'Class' in 'Basic_BCs'\n");
      Exit(0);
    }
    setKeywordOBC(str, i);
    BaseBc[i].set_Label(str);
    

    
    // 各条件に応じたパラメータをロード
    switch ( BaseBc[i].get_Class() ) {
      case OBC_WALL:
        get_OBC_Wall(label_leaf, i);
        break;
        
      case OBC_OUTFLOW:
        get_OBC_Outflow(label_leaf, i);
        break;
        
      case OBC_SPEC_VEL:
        get_OBC_SpecVH(label_leaf, i);
        break;
        
      case OBC_TRC_FREE://Traction_Free
        get_OBC_Trcfree(label_leaf, i);
        break;
        
      case OBC_FAR_FIELD://Far_Field
        get_OBC_FarField(label_leaf, i);
        break;
        
      case OBC_PERIODIC://Periodic
        get_OBC_Periodic(label_leaf, i);
        break;
        
      case OBC_SYMMETRIC:
        // nothing to do
        break;
    }
  }
  
  // 各フェイスに境界条件を設定する
  label_base = "/BC_Table/OuterBoundary/Face_BC";
  
  if ( !tpCntl->chkNode(label_base) ) {
    Hostonly_ printf("\tParsing error : Missing OuterBoundary Face_BC\n");
    Exit(0);
  }
  
  // check
  int nnode = tpCntl->countLabels(label_base);
  if ( nnode != NOFACE ) {
    Hostonly_ printf("\tParsing error : OuterBoundary Face_BC count != 6\n");
    Exit(0);
  }
  
  // 各面に与える境界条件番号を取得し，BaseBcから境界情報リストに内容をコピー
  for (int face=0; face<NOFACE; face++) {
    
    // faceに対するラベルを取得
    if ( !tpCntl->GetNodeStr(label_base, face+1, &str) ){
      Hostonly_ printf("\tGetNodeStr error\n");
      Exit(0);
    }
    label_leaf = label_base + "/" + str;
    
    // 指定の境界条件を探してBaseBC[]からbc[]へ内容のコピー
    label = label_leaf + "/kind";
    
    if ( !(tpCntl->GetValue(label, &str)) ) {
      Hostonly_ printf("\tParsing error : kind cannot found : Face_BC\n");
      Exit(0);
    }
    
    // Aliasでサーチ
    for (int i=0; i<NoBaseBC; i++) {
      if ( !strcasecmp( str.c_str(), BaseBc[i].get_Alias().c_str() ) ) {
        bc[face].dataCopy( &BaseBc[i] );
        break;
      }
      else {
        if ( i == NoBaseBC-1 ) { // 最後までみつからない
          Hostonly_ printf("\tParsing error : [%d]'%s' is not listed in 'Basic_BCs'\n", i+1, str.c_str());
          Exit(0);
        }
      }
    }
  }
  
  
  // ガイドセルの媒質インデクスをセット
  label_base = "/BC_Table/OuterBoundary/Face_BC";
  
  for (int face=0; face<NOFACE; face++) {
    
    if(!tpCntl->GetNodeStr(label_base, face+1, &str)){
      Hostonly_ printf("\tGetNodeStr error\n");
      Exit(0);
    }
    label_leaf = label_base + "/" + str;
    
    // ガイドセルの媒質ラベルを取得
    label = label_leaf + "/medium_on_guide_cell";
    
    if ( !(tpCntl->GetValue(label, &str )) ) {
      Hostonly_ printf("\tParsing error : No entory 'mediun_on_guide_cell' in 'Face_BC'\n");
      Exit(0);
    }
    
    for (int i=1; i<=NoMedium; i++) {
      if( !strcasecmp( str.c_str(), MTITP[i].label.c_str() ) ){
        bc[face].set_GuideMedium(i);
        break;
      }
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
      if ( bc[n].get_Class() == OBC_PERIODIC ) {
        n_pair = oppositDir(n);
        if ( bc[n_pair].get_Class() != OBC_PERIODIC ) {
          Hostonly_ printf("\tFace BC : No consistent Periodic Bnoudary in %s direction\n", FBUtility::getDirection(n_pair).c_str());
          Exit(0);
        }
      }
    }
    
    // 対になるモードのチェック
    for (unsigned n=0; n<NOFACE; n++) {
      if ( bc[n].get_Class() == OBC_PERIODIC ) {
        n_pair = oppositDir(n);
        
        switch (bc[n].get_PrdcMode()) {
          case BoundaryOuter::prdc_Simple:
            if ( bc[n_pair].get_PrdcMode() != BoundaryOuter::prdc_Simple ) { 
              Hostonly_ printf("\tFace BC : No consistent SIMPLE Periodic Bnoudary in %s direction\n", FBUtility::getDirection(n_pair).c_str());
              Exit(0);
            }
            break;
            
          case BoundaryOuter::prdc_Directional:
            if ( bc[n_pair].get_PrdcMode() != BoundaryOuter::prdc_Directional ) {
              Hostonly_ printf("\tFace BC : No consistent DIRECTIONAL Periodic Bnoudary in %s direction\n", FBUtility::getDirection(n_pair).c_str());
              Exit(0);
            }
            if ( bc[n].p != bc[n_pair].p ) { // 同じ値が入っていること
              Hostonly_ printf("\tFace BC : Pressure difference value is not same in %s direction\n", FBUtility::getDirection(n_pair).c_str());
              Exit(0);
            }
            if ( (bc[n].get_FaceMode() == BoundaryOuter::prdc_upstream) && (bc[n_pair].get_FaceMode() != BoundaryOuter::prdc_downstream) ) {
              Hostonly_ printf("\tFace BC : No consistent Upstream/Downstream relation in %s direction\n", FBUtility::getDirection(n).c_str());
              Exit(0);
            }
            if ( (bc[n].get_FaceMode() == BoundaryOuter::prdc_downstream) && (bc[n_pair].get_FaceMode() != BoundaryOuter::prdc_upstream) ) {
              Hostonly_ printf("\tFace BC : No consistent Upstream/Downstream relation in %s direction\n", FBUtility::getDirection(n).c_str());
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
      if ( bc[n].get_Class() == OBC_PERIODIC ) {
        if (bc[n].get_PrdcMode() == BoundaryOuter::prdc_Driver) {
          
          // 他方は周期境界以外であること
          if ( bc[n_pair].get_Class() == OBC_PERIODIC ) {
            Hostonly_ printf("\tFace BC : %s direction should be non periodic BC\n", FBUtility::getDirection(n_pair).c_str());
            Exit(0);
          }
          
          unsigned cflag=0;
          for (unsigned c=1; c<=NoBC; c++) {
            if ( compo[c].getType() == PERIODIC ) {
              if ( (int)compo[c].getPeriodicDir() != bc[n].get_DriverDir() ) {
                Hostonly_ printf("\tPeriodic Driver BC : No consistent Periodic Bnoudary in %s direction\n", FBUtility::getDirection(n_pair).c_str());
                Exit(0);
              }
              else {
                cflag++;
              }
            }
          }
          if (cflag != 1) {
            Hostonly_ printf("\tPeriodic Driver BC can not detemine uniquely\n");
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



/**
 @fn void ParseBC::printCompo(FILE* fp, REAL_TYPE* nv, int* gci, MediumList* mat)
 @brief コンポーネントの情報を表示する
 @param nv ボクセルモデルから計算した法線
 @param gci グローバルなコンポーネントのインデクス
 @param mat MediumList
 */
void ParseBC::printCompo(FILE* fp, REAL_TYPE* nv, int* gci, MediumList* mat)
{
  unsigned n, m;
  bool flag;
  
  // VBC ---------------------------------------------------
  if ( isComponent(SPEC_VEL) || isComponent(SPEC_VEL_WH) ) {
    fprintf(fp, "\n\t[SPECIFIED_VELOCITY]\n");
    fprintf(fp, "\t no                    Label   Mat   def    i_st    i_ed    j_st    j_ed    k_st    k_ed   Area[m*m]\n");
    
    for(n=1; n<=NoBC; n++){
      if ( (compo[n].getType() == SPEC_VEL) || (compo[n].getType() == SPEC_VEL_WH) )  {
        fprintf(fp,"\t%3d %24s %5d %5d %7d %7d %7d %7d %7d %7d %11.4e\n",
                n, compo[n].getLabel().c_str(), compo[n].getMatOdr(), compo[n].getDef(),
                getCmpGbbox_st_x(n, gci), getCmpGbbox_ed_x(n, gci), 
                getCmpGbbox_st_y(n, gci), getCmpGbbox_ed_y(n, gci), 
                getCmpGbbox_st_z(n, gci), getCmpGbbox_ed_z(n, gci), 
                compo[n].area);
      }
    }
    fprintf(fp, "\n");
    
    fprintf(fp, "\t no                    Label   Mat   normal_x   normal_y   normal_z      Direction");
    
    if ( compo[n].get_V_Profile() == CompoList::vel_zero ) {
      fprintf(fp, "\n");
    }
    else if ( compo[n].get_V_Profile() == CompoList::vel_constant ) {
      fprintf(fp, "    Q[m^3/s]   Vel.[m/s]\n");
    }
    else {
      fprintf(fp, "    Q[m^3/s]   Amp.[m/s]   Freq.[Hz]  Phase[rad] Intcpt[m/s] Strauhal[-]\n");
    }
    
    for(n=1; n<=NoBC; n++){
      if ( (compo[n].getType() == SPEC_VEL) || (compo[n].getType() == SPEC_VEL_WH) )  {
        fprintf(fp,"\t%3d %24s %5d %10.3e %10.3e %10.3e %14s ",
                n, compo[n].getLabel().c_str(), compo[n].getMatOdr(), compo[n].nv[0], compo[n].nv[1], compo[n].nv[2],
                (compo[n].getBClocation()==CompoList::same_direction) ? "same dir." : "opposite dir.");
        
        if ( compo[n].get_V_Profile() == CompoList::vel_zero ) {
          fprintf(fp,"\n");
        }
        else if ( compo[n].get_V_Profile() == CompoList::vel_constant ) {
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
      fprintf(fp, "\t no                    Label   Mat   def    i_st    i_ed    j_st    j_ed    k_st    k_ed      Temp(%s)      Temp[-]\n", (Unit_Temp==Unit_KELVIN) ? "K" : "C");
      
      for(n=1; n<=NoBC; n++) {
        if ( compo[n].getType() == SPEC_VEL_WH )  {
          fprintf(fp, "\t%3d %24s %5d %5d %7d %7d %7d %7d %7d %7d %12.4e %12.4e\n", 
                  n, compo[n].getLabel().c_str(), compo[n].getMatOdr(), compo[n].getDef(), 
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
    fprintf(fp, "\t no                    Label   Mat   def    i_st    i_ed    j_st    j_ed    k_st    k_ed  outflow_vel  pressure\n");
    
    for(n=1; n<=NoBC; n++){
      if ( compo[n].getType() == OUTFLOW )  {
        fprintf(fp,"\t%3d %24s %5d %5d %7d %7d %7d %7d %7d %7d ",
                n, compo[n].getLabel().c_str(), compo[n].getMatOdr(), compo[n].getDef(),
                getCmpGbbox_st_x(n, gci), getCmpGbbox_ed_x(n, gci), 
                getCmpGbbox_st_y(n, gci), getCmpGbbox_ed_y(n, gci), 
                getCmpGbbox_st_z(n, gci), getCmpGbbox_ed_z(n, gci), 
                (compo[n].getOutflowType() == V_AVERAGE) ? "Average " : "Minmax ");
        if (compo[n].get_P_BCtype() == P_DIRICHLET) {
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
    fprintf(fp, "\t no                    Label   Mat    i_st    i_ed    j_st    j_ed    k_st    k_ed\n");
    
    for(n=1; n<=NoBC; n++){
      if ( compo[n].getType() == IBM_DF )  {
        fprintf(fp,"\t%3d %24s %5d %7d %7d %7d %7d %7d %7d\n",
                n, compo[n].getLabel().c_str(), compo[n].getMatOdr(), 
                getCmpGbbox_st_x(n, gci), getCmpGbbox_ed_x(n, gci), 
                getCmpGbbox_st_y(n, gci), getCmpGbbox_ed_y(n, gci), 
                getCmpGbbox_st_z(n, gci), getCmpGbbox_ed_z(n, gci));
      }
    }
    fprintf(fp, "\n");
    
    fprintf(fp, "\t no                    Label   Mat   normal_x   normal_y   normal_z    Vel.[m/s]      Vel.[-]\n");
    
    for(n=1; n<=NoBC; n++){
      if ( compo[n].getType() == IBM_DF )  {
        fprintf(fp,"\t%3d %24s %5d %10.3e %10.3e %10.3e %12.4e %12.4e \n",
                n, compo[n].getLabel().c_str(), compo[n].getMatOdr(), compo[n].nv[0], compo[n].nv[1], compo[n].nv[2],
                compo[n].get_Velocity(), FBUtility::convD2ND_V(compo[n].get_Velocity(), RefVelocity));
      }
    }
  }
  
  // Forcing ---------------------------------------------------
	// Heat Exchanger
  if ( isComponent(HEX) ) {
    fprintf(fp, "\n\t[Heat Exchanger]\n");
    
    fprintf(fp, "\t no                    Label   Mat    normal_x   normal_y   normal_z     O_x[m]     O_y[m]     O_z[m]      dir_x      dir_y      dir_z\n");
    for(n=1; n<=NoBC; n++) {
      if ( compo[n].getType() == HEX ) {
        fprintf(fp, "\t%3d %24s %5d  %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e\n", 
                n, compo[n].getLabel().c_str(), compo[n].getMatOdr(),
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
                n, compo[n].getLabel().c_str(), compo[n].getMatOdr(),
								compo[n].depth, compo[n].shp_p1, compo[n].shp_p2, compo[n].area);
      }
    }
    fprintf(fp, "\n");
    
    fprintf(fp, "\t                                      i_st    i_ed    j_st    j_ed    k_st    k_ed\n");
    for(n=1; n<=NoBC; n++) {
      if ( compo[n].getType() == HEX ) {
        fprintf(fp, "\t%3d %24s %5d %7d %7d %7d %7d %7d %7d\n", 
                n, compo[n].getLabel().c_str(), compo[n].getMatOdr(),
                getCmpGbbox_st_x(n, gci), getCmpGbbox_ed_x(n, gci), 
                getCmpGbbox_st_y(n, gci), getCmpGbbox_ed_y(n, gci), 
                getCmpGbbox_st_z(n, gci), getCmpGbbox_ed_z(n, gci));
      }
    }
    fprintf(fp, "\n");
    
    fprintf(fp, "\t no                    Label   Mat         c1         c2         c3         c4  u_th[m/s]  thick[mm]     vec_forcing\n");
    for(n=1; n<=NoBC; n++) {
      if ( compo[n].getType() == HEX ) {
        fprintf(fp, "\t%3d %24s %5d %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e     %s\n", 
                n, compo[n].getLabel().c_str(), compo[n].getMatOdr(),
                compo[n].ca[0], compo[n].ca[1], compo[n].ca[2], compo[n].ca[3], compo[n].ca[4]*RefVelocity, compo[n].ca[5]*RefLength*1000.0,
                (compo[n].get_sw_HexDir()==ON) ? "Directional":"Non-directional");
      }
    }
  }
  
  // Fan
  if ( isComponent(FAN) ) {
    fprintf(fp, "\n\t[Fan]\n");
    
    fprintf(fp, "\t no                    Label   Mat    normal_x   normal_y   normal_z      O_x[m]     O_y[m]     O_z[m]\n");
    for(n=1; n<=NoBC; n++) {
      if ( compo[n].getType() == FAN ) {
        fprintf(fp, "\t%3d %24s %5d  %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e\n", 
                n, compo[n].getLabel().c_str(), compo[n].getMatOdr(),
								compo[n].nv[0], compo[n].nv[1], compo[n].nv[2],
                compo[n].oc[0], compo[n].oc[1], compo[n].oc[2]);
      }
    }
    fprintf(fp, "\n");
    
    fprintf(fp, "\t                                     Depth[m]     Fan[m]    Boss[m]  Area[m*m]\n");
    for(n=1; n<=NoBC; n++) {
      if ( compo[n].getType() == FAN ) {
        fprintf(fp, "\t%3d %24s %5d %10.3e %10.3e %10.3e %10.3e\n", 
                n, compo[n].getLabel().c_str(), compo[n].getMatOdr(),
								compo[n].depth, compo[n].shp_p1, compo[n].shp_p2, compo[n].area);
      }
    }
    fprintf(fp, "\n");
    
    fprintf(fp, "\t                                      i_st    i_ed    j_st    j_ed    k_st    k_ed\n");
    for(n=1; n<=NoBC; n++) {
      if ( compo[n].getType() == FAN ) {
        fprintf(fp, "\t%3d %24s %5d %7d %7d %7d %7d %7d %7d\n", 
                n, compo[n].getLabel().c_str(), compo[n].getMatOdr(),
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
     n, compo[n].getLabel().c_str(), compo[n].getMatOdr(),
     compo[n].ca[0], compo[n].ca[1], compo[n].ca[2], compo[n].ca[3], compo[n].ca[4]*RefVelocity, compo[n].ca[5]*RefLength*1000.0,
     (compo[n].get_sw_HexDir()==ON) ? "Directional":"Non-directional");
     }
     }*/
  }
  
	// Darcy Law
  if ( isComponent(DARCY) ) {
    fprintf(fp, "\n\t[Darcy medium]\n");
    
    fprintf(fp, "\t no                    Label   Mat        (Computed Normal form ID)  Area[m*m]   Prmblty_x   Prmblty_y   Prmblty_z[m^2]    Prmblty_x   Prmblty_y   Prmblty_z[-]\n");
    for(n=1; n<=NoBC; n++) {
      if ( compo[n].getType() == DARCY ) {
        fprintf(fp, "\t%3d %24s %5d %10.3e %10.3e %10.3e %10.3e %11.4e %11.4e %11.4e       %11.4e %11.4e %11.4e\n", 
                n, compo[n].getLabel().c_str(), compo[n].getMatOdr(), 
                nv[3*n+0], nv[3*n+1], nv[3*n+2], compo[n].area, 
								compo[n].ca[0], compo[n].ca[1], compo[n].ca[2], compo[n].ca[3], compo[n].ca[4], compo[n].ca[5]);
      }
    }
    
    fprintf(fp, "\t                                      i_st    i_ed    j_st    j_ed    k_st    k_ed\n");
    for(n=1; n<=NoBC; n++) {
      if ( compo[n].getType() == DARCY ) {
        fprintf(fp, "\t%3d %24s %5d %7d %7d %7d %7d %7d %7d\n", 
                n, compo[n].getLabel().c_str(), compo[n].getMatOdr(),
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
    fprintf(fp, "\t no                    Label   Mat   def\n");
    
    for(n=1; n<=NoBC; n++) {
      if ( compo[n].getType() == ADIABATIC ) {
        fprintf(fp, "\t%3d %24s %5d %5d\n", n, compo[n].getLabel().c_str(), compo[n].getMatOdr(), compo[n].getDef());
      }
    }
  }
  
  // Direct Heat Flux
  if ( HeatProblem && isComponent(HEATFLUX) ) {
    fprintf(fp, "\n\t[Direct Heat Flux]\n");
    fprintf(fp, "\t no                    Label   Mat   def    i_st    i_ed    j_st    j_ed    k_st    k_ed   Area[m*m]   flux(W/m^2)        q[-]\n");
    
    for(n=1; n<=NoBC; n++) {
      if ( compo[n].getType() == HEATFLUX ) {
        fprintf(fp, "\t%3d %24s %5d %5d %7d %7d %7d %7d %7d %7d %11.4e %12.4e %12.4e\n", 
                n, compo[n].getLabel().c_str(), compo[n].getMatOdr(), compo[n].getDef(),
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
    fprintf(fp, "\t no                    Label   Mat   def    i_st    i_ed    j_st    j_ed    k_st    k_ed   Area[m*m]   coef(W/m^2K)\n");
    
    for(n=1; n<=NoBC; n++) {
      if ( compo[n].getHtype() == HT_N ) {
        fprintf(fp, "\t%3d %24s %5d %5d %7d %7d %7d %7d %7d %7d %11.4e %12.4e\n", 
                n, compo[n].getLabel().c_str(), compo[n].getMatOdr(), compo[n].getDef(), 
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
    fprintf(fp, "\t no                    Label   Mat   def    i_st    i_ed    j_st    j_ed    k_st    k_ed   Area[m*m]   coef(W/m^2K)   Temp(%s)   Temp[-]\n", (Unit_Temp==Unit_KELVIN) ? "K" : "C");
    
    for(n=1; n<=NoBC; n++) {
      if ( compo[n].getHtype() == HT_S ) {
        fprintf(fp, "\t%3d %24s %5d %5d %7d %7d %7d %7d %7d %7d %11.4e %12.4e %12.4e %12.4e\n", 
                n, compo[n].getLabel().c_str(), compo[n].getMatOdr(), compo[n].getDef(), 
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
    fprintf(fp, "\t no                    Label   Mat   def    i_st    i_ed    j_st    j_ed    k_st    k_ed   Area[m*m]      Temp(%s)      Temp[-]   Type\n", (Unit_Temp==Unit_KELVIN) ? "K" : "C");
    
    for(n=1; n<=NoBC; n++) {
      if ( compo[n].getHtype() == HT_SN ) {
        fprintf(fp, "\t%3d %24s %5d %5d %7d %7d %7d %7d %7d %7d %11.4e %12.4e %12.4e   %s\n", 
                n, compo[n].getLabel().c_str(), compo[n].getMatOdr(), compo[n].getDef(), 
                getCmpGbbox_st_x(n, gci), getCmpGbbox_ed_x(n, gci), 
                getCmpGbbox_st_y(n, gci), getCmpGbbox_ed_y(n, gci), 
                getCmpGbbox_st_z(n, gci), getCmpGbbox_ed_z(n, gci), 
                compo[n].area, FBUtility::convK2Temp(compo[n].get_Temp(), Unit_Temp),
                FBUtility::convK2ND(compo[n].get_Temp(), BaseTemp, DiffTemp),
                (compo[n].get_sw_HTmodeRef()==CompoList::HT_mode_bulk) ? "Bulk" : "Local");
      }
    }
    fprintf(fp, "\t no                    Label   Mat  vert_lam_a  vert_lam_b vert_turb_a vert_turb_b  vert_Ra_cr   lwr_lam_a   lwr_lam_b  lwr_turb_a  lwr_turb_b   lwr_Ra_cr\n");
    for(n=1; n<=NoBC; n++) {
      if ( compo[n].getHtype() == HT_SN ) {
        fprintf(fp, "\t%3d %24s %5d %11.4e %11.4e %11.4e %11.4e %11.4e %11.4e %11.4e %11.4e %11.4e %11.4e\n", 
                n, 
                compo[n].getLabel().c_str(), 
                compo[n].getMatOdr(), 
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
    fprintf(fp, "\t no                    Label   Mat   def    i_st    i_ed    j_st    j_ed    k_st    k_ed   Area[m*m]      Temp(%s)      Temp[-]   Type\n", 
            (Unit_Temp==Unit_KELVIN) ? "K" : "C");
    
    for(n=1; n<=NoBC; n++) {
      if ( compo[n].getHtype() == HT_SF ) {
        fprintf(fp, "\t%3d %24s %5d %5d %7d %7d %7d %7d %7d %7d %11.4e %12.4e %12.4e\n", 
                n, 
                compo[n].getLabel().c_str(), 
                compo[n].getMatOdr(), 
                compo[n].getDef(), 
                getCmpGbbox_st_x(n, gci), getCmpGbbox_ed_x(n, gci), 
                getCmpGbbox_st_y(n, gci), getCmpGbbox_ed_y(n, gci), 
                getCmpGbbox_st_z(n, gci), getCmpGbbox_ed_z(n, gci), 
                compo[n].area, FBUtility::convK2Temp(compo[n].get_Temp(), Unit_Temp),
                FBUtility::convK2ND(compo[n].get_Temp(), BaseTemp, DiffTemp),
                (compo[n].get_sw_HTmodeRef()==CompoList::HT_mode_bulk) ? "Bulk" : "Local");
      }
    }
    fprintf(fp, "\t no                    Label   Mat       alpha        beta       gamma\n");
    for(n=1; n<=NoBC; n++) {
      if ( compo[n].getHtype() == HT_SN ) {
        fprintf(fp, "\t%3d %24s %5d %11.4e %11.4e %11.4e\n", 
                n, 
                compo[n].getLabel().c_str(), 
                compo[n].getMatOdr(), 
                compo[n].ca[CompoList::alpha], 
                compo[n].ca[CompoList::beta], 
                compo[n].ca[CompoList::gamma]);
      }
    }
  }
  
  // Heat Transfer B
  if ( HeatProblem && isCompoTransfer(HT_B)) {
    fprintf(fp, "\n\t[Heat Transfer : Type B]\n");
    fprintf(fp, "\t no                    Label   Mat   def    i_st    i_ed    j_st    j_ed    k_st    k_ed   Area[m*m]  coef(W/m^2K)  BulkTemp(%s)   BulkTemp[-]\n", 
            (Unit_Temp==Unit_KELVIN) ? "K" : "C");
    
    for(n=1; n<=NoBC; n++) {
      if ( compo[n].getHtype() == HT_B ) {
        fprintf(fp, "\t%3d %24s %5d %5d %7d %7d %7d %7d %7d %7d %11.4e  %12.4e %12.4e %12.4e\n", 
                n, compo[n].getLabel().c_str(), compo[n].getMatOdr(), compo[n].getDef(), 
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
    fprintf(fp, "\t no                    Label   Mat   def    i_st    i_ed    j_st    j_ed    k_st    k_ed   Area[m*m]   Sf.Temp(%s)   Sf.Temp[-]\n", 
            (Unit_Temp==Unit_KELVIN) ? "K" : "C");
    
    for(n=1; n<=NoBC; n++) {
      if ( compo[n].getType() == ISOTHERMAL ) {
        fprintf(fp, "\t%3d %24s %5d %5d %7d %7d %7d %7d %7d %7d %11.4e %12.4e %12.4e \n", 
                n, compo[n].getLabel().c_str(), compo[n].getMatOdr(), compo[n].getDef(), 
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
    fprintf(fp, "\t no                    Label   Mat    def    i_st    i_ed    j_st    j_ed    k_st    k_ed   Area[m*m]   ep[-]   pj[-]\n");
    
    for(n=1; n<=NoBC; n++) {
      if ( compo[n].getType() == RADIANT ) {
        fprintf(fp, "\t%3d %24s %5d %5d %7d %7d %7d %7d %7d %7d %11.4e %12.4e %12.4e\n", 
                n, compo[n].getLabel().c_str(), compo[n].getMatOdr(), compo[n].getDef(), 
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
    fprintf(fp, "\t no                    Label   Mat    i_st    i_ed    j_st    j_ed    k_st    k_ed     Q[W/m^3]    nrmlzd[-]\n");
    
    for(n=1; n<=NoBC; n++) {
      unsigned h_odr = compo[n].getMatOdr();
      if ( compo[n].getType() == HEAT_SRC ) {
        fprintf(fp, "\t%3d %24s %5d %7d %7d %7d %7d %7d %7d %12.4e %12.4e\n", 
                n, compo[n].getLabel().c_str(), compo[n].getMatOdr(), 
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
    fprintf(fp, "\t no                    Label   Mat    i_st    i_ed    j_st    j_ed    k_st    k_ed      Temp[%s]      Temp[-]\n", 
            (Unit_Temp==Unit_KELVIN) ? "K" : "C");
    
    for(n=1; n<=NoBC; n++) {
      if ( compo[n].getType() == CNST_TEMP ) {
        fprintf(fp, "\t%3d %24s %5d %7d %7d %7d %7d %7d %7d %12.4e %12.4e\n", 
                n, compo[n].getLabel().c_str(), compo[n].getMatOdr(), 
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
    fprintf(fp, "\t no                    Label   Mat    i_st    i_ed    j_st    j_ed    k_st    k_ed   Area[m*m]  Variables\n");
    
    for(n=1; n<=NoBC; n++) {
      if ( compo[n].getType() == CELL_MONITOR ) {
        fprintf(fp, "\t%3d %24s %5d %7d %7d %7d %7d %7d %7d %11.4e  %s\n", 
                n, compo[n].getLabel().c_str(), compo[n].getMatOdr(), 
                getCmpGbbox_st_x(n, gci), getCmpGbbox_ed_x(n, gci), 
                getCmpGbbox_st_y(n, gci), getCmpGbbox_ed_y(n, gci), 
                getCmpGbbox_st_z(n, gci), getCmpGbbox_ed_z(n, gci), 
								compo[n].area, compo[n].getVarStr().c_str());
      }
    }
    fprintf(fp, "\n");
    
    fprintf(fp, "\t no                    Label   Mat   normal_x   normal_y   normal_z Reference\n");
    for(n=1; n<=NoBC; n++){
      if ( compo[n].getType() == CELL_MONITOR )  {
        fprintf(fp,"\t%3d %24s %5d %10.3e %10.3e %10.3e       %3s\n",
                n, compo[n].getLabel().c_str(), compo[n].getMatOdr(), compo[n].nv[0], compo[n].nv[1], compo[n].nv[2],
                (compo[n].getStateCellMonitor()==ON) ? "yes" : "no");
      }
    }
  }
  
  // Periodic ---------------------------------------------------
  if ( isComponent(PERIODIC) ) {
    fprintf(fp, "\n\t[Periodic]\n");
    fprintf(fp, "\t no                    Label   Mat    i_st    i_ed    j_st    j_ed    k_st    k_ed    Pressure Difference [Pa]/[-]  Driver\n");
    
    int dir_in=0, dir_out=0, pp_in=0, pp_out=0;
    
    for(n=1; n<=NoBC; n++) {
      if ( compo[n].getType() == PERIODIC ) {
        dir_in = compo[n].getPeriodicDir();
        
        fprintf(fp, "\t%3d %24s %5d %7d %7d %7d %7d %7d %7d     ", 
                n, compo[n].getLabel().c_str(), compo[n].getMatOdr(), 
                getCmpGbbox_st_x(n, gci), getCmpGbbox_ed_x(n, gci), 
                getCmpGbbox_st_y(n, gci), getCmpGbbox_ed_y(n, gci), 
                getCmpGbbox_st_z(n, gci), getCmpGbbox_ed_z(n, gci));
        fprintf(fp,"%12.6e / %12.6e ", compo[n].ca[0], FBUtility::convD2ND_P(compo[n].ca[0], BasePrs, RefDensity, RefVelocity, Unit_Prs));
        fprintf(fp, "%7s\n", FBUtility::getDirection(dir_in).c_str());
        
        // ドライバの方向が周期境界であるかをチェック
        if ( bc[dir_in].get_Class() != OBC_PERIODIC ) {
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
 @fn void ParseBC::printFaceOBC(FILE* fp, REAL_TYPE* G_Lbx)
 @brief 外部境界条件の各面の情報を表示する
 @param fp
 @param G_Lbx グローバルの領域の大きさ
 */
void ParseBC::printFaceOBC(FILE* fp, REAL_TYPE* G_Lbx)
{
  for (int i=0; i<NOFACE; i++) {
    fprintf(fp,"\t      Set %s up as %s : < %s >\n", Control::getDirection(i).c_str(), getOBCstr(bc[i].get_Class()).c_str(), bc[i].get_Alias().c_str());
    printOBC(fp, &bc[i], G_Lbx, i);
    fprintf(fp,"\n");
  }
  fflush(fp);
}

/**
 @fn void ParseBC::printOBC(FILE* fp, BoundaryOuter* ref REAL_TYPE* G_Lbx, const int face)
 @brief 速度の外部境界条件処理の表示
 @param fp
 @param ref
 @param G_Lbx グローバルの領域の大きさ
 @param face 面番号
 */
void ParseBC::printOBC(FILE* fp, BoundaryOuter* ref, REAL_TYPE* G_Lbx, const int face)
{
  REAL_TYPE a, b, c;
  
  fprintf(fp,"\t\t\tGuide Cell Medium = %s\n", mat[ref->get_GuideMedium()].getLabel().c_str());
  
  switch ( ref->get_Class() ) {
    case OBC_WALL:
      if ( ref->get_V_Profile() == CompoList::vel_harmonic ) {
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
      if ( ref->get_V_Profile() == CompoList::vel_harmonic ) {
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
 @fn void ParseBC::receiveCompoPtr(CompoList* CMP)
 @brief コンポーネントリストのポインタを受け取る
 */
void ParseBC::receiveCompoPtr(CompoList* CMP)
{
  if ( !CMP ) {
    Hostonly_ stamped_printf("\tAn object of CompoList class is NULL\n");
    Exit(0);
  }
  compo = CMP;
}


/**
 @fn void ParseBC::receive_TP_Ptr(TPControl* tp)
 @brief TPのポインタを受け取る
 */
bool ParseBC::receive_TP_Ptr(TPControl* tp) 
{ 
  if ( !tp ) return false;
  tpCntl = tp;
  return true;
}


//@brief 変数の初期化
void ParseBC::setControlVars(Control* Cref, BoundaryOuter* ptr, MediumList* m_mat)
{
  if ( !ptr ) Exit(0);
  bc = ptr;
  
  if ( !m_mat ) Exit(0);
  mat = m_mat;
  
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
  NoMedium    = Cref->NoMedium;
  NoCompo     = Cref->NoCompo;
  NoBC        = Cref->NoBC;
  
  int m;
  double s, two=2.0d;
  
  s = (double)MASK_6; // bit幅マスクは2^(bit幅)-1を表し，ちょうど0を除いた個数となっている
  m = (unsigned)(log10(s+1.0d)/log10(two) );
  
  if ( NoCompo > s ) {
    printf("Error : No. of Component (= NoBC + NoMedium) must be less or equal %d(%dbit-width)\n", (int)s, m);
    Exit(0);
  }
  
  s = (double)(MASK_5-1); // 0と31を除く
  m = (int)( log10(s+2.0d)/log10(two) );
  
  if ( NoBC > s ) {
    printf("Error : No. of BC must be less or equal %d(%dbit-width)\n", (int)s, m);
    Exit(0);
  }
  
  // 外部境界条件の種類数を取得し，内部保持配列をインスタンス
  std::string label;
  std::string str;
  
  label = "/BC_Table/OuterBoundary";
  
  if ( !tpCntl->chkNode(label) ) {
    stamped_printf("\tParsing error : Missing OuterBoundary tree\n");
    Exit(0);
  }
  
  // タグ内の数をチェック
  int nnode = tpCntl->countLabels(label);
  if ( nnode == 0 ) {
    stamped_printf("\tNo OuterBoundary defined\n");
    return;
  }
  
  int counter=0;
  for(int i=1; i<=nnode; i++){
    
    if(!tpCntl->GetNodeStr(label, i, &str)){
      stamped_printf("\tGetNodeStr error\n");
      Exit(0);
    }
    
    if( !strcasecmp(str.substr(0,9).c_str(), "Basic_BCs") ) counter++;
  }
  
  NoBaseBC = counter;
  BaseBc = new BoundaryOuter[NoBaseBC];
}



/**
 @fn void ParseBC::set_Deface(const std::string label_base, const int n)
 @brief 内部境界条件のdef_faceを取得し，登録する
 @param label_base
 @param n コンポーネントリストのエントリ番号
 */
void ParseBC::set_Deface(const std::string label_base, const int n)
{
  int def=0;
  std::string label;
  
  label=label_base+"/def_face";//
  
  if ( !(tpCntl->GetValue(label, &def )) ) {
	  stamped_printf("\tParsing error : fail to get 'def_face' in '%s'\n", label.c_str());
	  Exit(0);
  }
  
  // IDが媒質に含まれているかを調べる
  if ( (0 >= def) || (def > NoMedium) ) {
	  stamped_printf("\tParsing error : ID[%d] described in 'def_face' is not included in 'Model_Setting'\n", def);
	  Exit(0);
  }
  compo[n].setDef(def);
}



//@brief 内部境界条件の照合を行う
//@note SPEC_VEL_WHは陽には現れず，get_IBC_SpecVel()内で登録される
void ParseBC::setKeywordLBC(const std::string keyword, const int m)
{
  if     ( FBUtility::compare(keyword, "Adiabatic") )             compo[m].setType(ADIABATIC);
  else if( FBUtility::compare(keyword, "Direct_Heat_Flux") )      compo[m].setType(HEATFLUX);
  else if( FBUtility::compare(keyword, "HeatTransfer_B") ) {      compo[m].setType(TRANSFER); compo[m].setHtype(HT_B); }
  else if( FBUtility::compare(keyword, "HeatTransfer_N") ) {      compo[m].setType(TRANSFER); compo[m].setHtype(HT_N); }
  else if( FBUtility::compare(keyword, "HeatTransfer_S") ) {      compo[m].setType(TRANSFER); compo[m].setHtype(HT_S); }
  else if( FBUtility::compare(keyword, "HeatTransfer_SF") ){      compo[m].setType(TRANSFER); compo[m].setHtype(HT_SF); }
  else if( FBUtility::compare(keyword, "HeatTransfer_SN") ){      compo[m].setType(TRANSFER); compo[m].setHtype(HT_SN); }
  else if( FBUtility::compare(keyword, "IsoThermal") )            compo[m].setType(ISOTHERMAL);
  else if( FBUtility::compare(keyword, "Radiation") )             compo[m].setType(RADIANT);
  else if( FBUtility::compare(keyword, "specified_velocity") )    compo[m].setType(SPEC_VEL);
  else if( FBUtility::compare(keyword, "outflow") )               compo[m].setType(OUTFLOW);
  else if( FBUtility::compare(keyword, "Forcing") )               compo[m].setType(IBM_DF);
  else if( FBUtility::compare(keyword, "Heat_Source") )           compo[m].setType(HEAT_SRC);
  else if( FBUtility::compare(keyword, "specified_Temperature") ) compo[m].setType(CNST_TEMP);
  else if( FBUtility::compare(keyword, "Pressure_Loss") )         compo[m].setType(HEX);
  else if( FBUtility::compare(keyword, "Fan") )                   compo[m].setType(FAN);
  else if( FBUtility::compare(keyword, "Darcy") )                 compo[m].setType(DARCY);
  else if( FBUtility::compare(keyword, "cell_monitor") )          compo[m].setType(CELL_MONITOR);
  else if( FBUtility::compare(keyword, "inactive") )              compo[m].setType(INACTIVE);
  else if( FBUtility::compare(keyword, "Periodic") )              compo[m].setType(PERIODIC);
  else {
    Hostonly_ stamped_printf("\tInvalid keyword is described '%s'\n", keyword.c_str());
    Exit(0);
  }
}


//@brief 外部境界条件のキーワードを照合し，コードを登録する
void ParseBC::setKeywordOBC(const std::string keyword, const int m)
{
  if     ( FBUtility::compare(keyword, "Wall") )               BaseBc[m].set_Class(OBC_WALL);
  else if( FBUtility::compare(keyword, "Outflow") )            BaseBc[m].set_Class(OBC_OUTFLOW);
  else if( FBUtility::compare(keyword, "Specified_Velocity") ) BaseBc[m].set_Class(OBC_SPEC_VEL);
  else if( FBUtility::compare(keyword, "Symmetric") )          BaseBc[m].set_Class(OBC_SYMMETRIC);
  else if( FBUtility::compare(keyword, "Periodic") )           BaseBc[m].set_Class(OBC_PERIODIC);
  else if( FBUtility::compare(keyword, "Traction_Free") )      BaseBc[m].set_Class(OBC_TRC_FREE);
  else if( FBUtility::compare(keyword, "Far_Field") )          BaseBc[m].set_Class(OBC_FAR_FIELD);
  else {
    Hostonly_ stamped_printf("\tParsing error : Invalid keyword is described '%s'\n", keyword.c_str());
    Exit(0);
  }
}



/**
 @fn void ParseBC::setMediumPoint(MediumTableInfo *m_MTITP)
 @brief MediumTableInfoをポイント
 @param m_MTITP
 */
void ParseBC::setMediumPoint(MediumTableInfo *m_MTITP)
{
  if ( !m_MTITP ) Exit(0);
  MTITP = m_MTITP;
}


// 作業用ポインタのコピー
void ParseBC::setPartitionManager(cpm_ParaManager* m_paraMngr)
{
  if ( !m_paraMngr ) Exit(0);
  paraMngr = m_paraMngr;
}


/**
 @fn void ParseBC::setRefValue(MediumList* mat, CompoList* cmp, Control* C)
 @brief 媒質により決まる代表量をコピーする
 @param mat MediumList class
 @param cmp CompoList class
 @param C Control class
 */
void ParseBC::setRefValue(MediumList* mat, CompoList* cmp, Control* C)
{
  unsigned m;
  
	for (unsigned n=NoBC+1; n<=NoCompo; n++) {
    if ( cmp[n].getMatOdr() == C->RefMat ) {
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
 @fn void ParseBC::setRefMedium(MediumList* mat, Control* Cref)
 @brief 指定した媒質IDから参照物理量を求める
 @param mat MediumList
 @param Cref Control class
 */
void ParseBC::setRefMedium(MediumList* mat, Control* Cref)
{
  unsigned m;
  
  for (unsigned n=Cref->NoBC+1; n<=Cref->NoCompo; n++) {
    if ( compo[n].getMatOdr() == Cref->RefMat ) {
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

