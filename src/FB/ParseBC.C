// #################################################################
//
// CAERU Library
//
// Copyright (c) 2012-2013  All right reserved.
//
// Institute of Industrial Science, University of Tokyo, Japan. 
//
// #################################################################

//@file   ParseBC.C
//@brief  FlowBase ParseBC class
//@author kero

#include <math.h>
#include "ParseBC.h"


// #################################################################
/**
 * @brief KOSと境界条件数の整合性をチェックする
 * @param [in] kos KindOfSolver
 * @param [in] cmp CompListクラスのポインタ
 */
void ParseBC::chkBCconsistency(const int kos, CompoList* cmp)
{
  if (kos == FLOW_ONLY) 
  {
    for (int n=1; n<=NoBC; n++) {
      if ( cmp[n].isHBC() ) 
      {
        Hostonly_ stamped_printf("\tNo consistency between 'KindOfSolver' and 'LocalBoundary'\n");
        Exit(0);
      }
    }
  }
  else if (kos == SOLID_CONDUCTION) 
  {
    for (int n=1; n<=NoBC; n++) {
      if ( cmp[n].isVBC() ) 
      {
        Hostonly_ stamped_printf("\tNo consistency between 'KindOfSolver' and 'LocalBoundary'\n");
        Exit(0);
      }
    }
  }
}


// #################################################################
/**
 * @brief ラベルの重複を調べる
 * @param [in] n       テストするBaseBcの格納番号の最大値
 * @param [in] m_label テストラベル
 */
bool ParseBC::chkDuplicate(const int n, const string m_label)
{
	for (int i=0; i<n; i++){
    if ( BaseBc[i].get_Alias() == m_label ) return false;
	}
	return true;
}


// #################################################################
/**
 * @brief KOSと媒質の状態の整合性をチェックし，媒質数をカウント，C.NoMediumFluid, C.NoMediumSolidをセット
 * @param [in] Cref Controlクラス
 * @param [in] mat  MediumList
 */
void ParseBC::countMedium(Control* Cref, const MediumList* mat)
{
  // check at least one fluid
  if ( KindOfSolver != SOLID_CONDUCTION ) {
    bool check=false;
    for (int i=1; i<=NoMedium; i++) {
      if ( mat[i].getState() == FLUID ) check = true;
    }
    if ( !check ) {
      Hostonly_ stamped_printf("\tAnalysis model should have at least one FLUID medium in MediumTable.\n");
      Exit(0);
    }
  }
  
  // 流体と固体の媒質数をセット
  int m_fluid=0, m_solid=0;
  for (int i=1; i<=NoMedium; i++) {
    if ( mat[i].getState() == SOLID ) m_solid++;
    else m_fluid++;
  }
  
  Cref->NoMediumFluid = m_fluid;
  Cref->NoMediumSolid = m_solid;
}


// #################################################################
/**
 * @brief 外部境界条件のキーワードを照合し， BCの文字列を返す
 * @param [in] id
 */
string ParseBC::getOBCstr(const int id)
{
  string bc;
  if     ( id == OBC_WALL )      bc = "Wall";
  else if( id == OBC_OUTFLOW )   bc = "Outflow";
  else if( id == OBC_SPEC_VEL )  bc = "SpecifiedVelocity";
  else if( id == OBC_SYMMETRIC ) bc = "Symmetric";
  else if( id == OBC_PERIODIC )  bc = "Periodic";
  else if( id == OBC_TRC_FREE )  bc = "TractionFree";
  else if( id == OBC_FAR_FIELD ) bc = "FarField";
  else if( id == OBC_INTRINSIC ) bc = "Intrinsic";
  else                           bc = "";
  return bc;
}


// #################################################################
/**
 * @brief 単位ベクトルを計算して戻す
 * @param [in,out] v ベクトル値
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



// #################################################################
// LocalBoundaryタグ直下のBCの個数（内部境界条件数）を返す
// 境界条件数がゼロでもエラーではない
int ParseBC::getNoLocalBC()
{  
  int nobc=0;
  string str;
  string label;
  
  label="/BcTable/LocalBoundary";
  
  if ( tpCntl->chkNode(label) )  //nodeがあれば
  {
	  nobc = tpCntl->countLabels(label);
  }
  
  return nobc;
}


// #################################################################
/**
 * @brief 境界条件の値(REAL_TYPE型)を取得し，返す
 * @param [in] label テストラベル
 */
REAL_TYPE ParseBC::get_BCval_real(const string label)
{
  REAL_TYPE df=0.0f;
  if ( !(tpCntl->GetValue(label, &df )) ) 
  {
    stamped_printf("\tParsing error : Invalid REAL_TYPE value for '%s'\n", label.c_str());
    Exit(0);
  }
  return df;
}


// #################################################################
/**
 * @brief 内部境界条件の座標値を取得し，登録する
 * @param [in]  label_base ラベルディレクトリ
 * @param [out] v          ベクトルパラメータ
 */
void ParseBC::get_Center(const string label_base, REAL_TYPE* v)
{
  string label;
  for (int i=0; i<3; i++) v[i]=0.0f;
  
  label = label_base + "/Center";
  
  if( !(tpCntl->GetVector(label, v, 3)) )
  {
    stamped_printf("\tParsing error : fail to get vec params in '%s\n", label.c_str());
    Exit(0);
  }
}


// #################################################################
// Darcyのパラメータを取得する
void ParseBC::get_Darcy(const string label_base, const int n, CompoList* cmp)
{
  REAL_TYPE v[3];
  string label;
  
  for (int n=0; n<3; n++) v[n]=0.0;
  
  
  // 透過率の取得
  label = label_base + "/Permeability";
  
  if( !(tpCntl->GetVector(label, v, 3)) ) {
    stamped_printf("\tParsing error : fail to get permeability params in 'Darcy'\n");
    Exit(0);
  }
  cmp[n].ca[0] = v[0]; // 透過率[m^2]は境界条件設定時に無次元化する
  cmp[n].ca[1] = v[1];
  cmp[n].ca[2] = v[2];
}


// #################################################################
/**
 * @brief 内部境界条件の方向ベクトル値を取得し，登録する
 * @param [in] label_base ラベルディレクトリ
 * @param [out v          ベクトルパラメータ
 */
void ParseBC::get_Dir(const string label_base, REAL_TYPE* v)
{
  string label;
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



// #################################################################
// Const_Temperatureのパラメータを取得する
void ParseBC::get_IBC_CnstTemp(const string label_base, const int n, CompoList* cmp)
{
  string label = label_base + "/Temperature";
  
  REAL_TYPE tmp = get_BCval_real(label);
  cmp[n].set_Temp( FBUtility::convTemp2K(tmp, Unit_Temp) );
}



// #################################################################
// Fanのパラメータを取得する
void ParseBC::get_IBC_Fan(const string label_base, const int n, CompoList* cmp)
{
  string str,str_u;
  string label;
  REAL_TYPE v[4], ct;

  
  // 入力単位の指定
  label=label_base+"/Unit";//
  //cout <<  "label : " << label << endl;
  if ( !(tpCntl->GetValue(label, &str_u )) ) {
		stamped_printf("\tParsing error : Invalid float value for 'unit' in 'Pressure_Loss'\n");
		Exit(0);
  }
  if ( !strcasecmp("mmaq", str_u.c_str()) ) {
    cmp[n].setPrsUnit(CompoList::unit_mmAq);
  }
  else if ( !strcasecmp("mmhg", str_u.c_str()) ) {
    cmp[n].setPrsUnit(CompoList::unit_mmHg);
  }
  else if ( !strcasecmp("pa", str_u.c_str()) ) {
    cmp[n].setPrsUnit(CompoList::unit_Pa);
  }
  else if ( !strcasecmp("nondimension", str_u.c_str()) ) {
    cmp[n].setPrsUnit(CompoList::unit_NonDimensional);
  }
  else {
    stamped_printf("\tDescribed unit is out of scope.\n");
    Exit(0);
  }
  
  // 法線ベクトルの取得
  get_NV(label_base, v);
  copyVec(cmp[n].nv, v);
  
  // 中心座標の取得
  get_Center(label_base, v);
  copyVec(cmp[n].oc, v);
  
  // 形状パラメータ
  label=label_base+"/Depth";
  cmp[n].depth  = get_BCval_real(label);
  label=label_base+"/FanRadius";
  cmp[n].shp_p1 = get_BCval_real(label);
  label=label_base+"/BossRadius";
  cmp[n].shp_p2 = get_BCval_real(label);
  
  if ( cmp[n].shp_p1 <= cmp[n].shp_p2 ) {
    stamped_printf("\tError : Radius of boss is greater than fan.\n");
    Exit(0);
  }
  
}



// #################################################################
// Direct_Fluxのパラメータを取得する
void ParseBC::get_IBC_HeatFlux(const string label_base, const int n, CompoList* cmp)
{
  string label = label_base + "/HeatFlux";
  
  cmp[n].set_Heatflux( get_BCval_real(label) ); /// @note [W/m^2]
}


// #################################################################
// Heat_Generationのパラメータを取得する
void ParseBC::get_IBC_HeatSrc(const string label_base, const int n, CompoList* cmp)
{
  REAL_TYPE hsrc=0.0f;
  string str;
  string label;
  
  // type
  label = label_base + "/Type";
  
  if ( !(tpCntl->GetValue(label, &str )) )
  {
    stamped_printf("\tParsing error : Invalid int value for '%s\n", label.c_str());
    Exit(0);
  }
  if ( !strcasecmp("HeatReleaseValue", str.c_str()) )
  {
	  cmp[n].set_HSRC_policy(true);
  }
  else if ( !strcasecmp("HeatGenerationDensity", str.c_str()) )
  {
	  cmp[n].set_HSRC_policy(false);
  }
  else
  {
	  stamped_printf("\tParsing error : Invalid string value for 'Type' : %s\n", str.c_str());
	  Exit(0);
  }
  
  // 放熱量
  label= label_base + "/Value";
  
  if ( !(tpCntl->GetValue(label, &hsrc )) )
  {
    stamped_printf("\tParsing error : Invalid float value for '%s\n", label.c_str());
    Exit(0);
  }
  
  if ( cmp[n].isPolicy_HeatDensity() )
  {
    cmp[n].set_HeatDensity( hsrc ); // 発熱密度
  }
  else
  {
    cmp[n].set_HeatValue( hsrc ); // 発熱量
  }
  
}



// #################################################################
// HeatTransfer_Nのパラメータを取得する
void ParseBC::get_IBC_HT_N(const string label_base, const int n, CompoList* cmp)
{
  // 熱伝達係数
  string label = label_base + "/CoefOfHeatTransfer";
  
  cmp[n].set_CoefHT( get_BCval_real(label) );
}


// #################################################################
// HeatTransfer_Sのパラメータを取得す
void ParseBC::get_IBC_HT_S(const string label_base, const int n, CompoList* cmp, const MediumList* mat)
{
  string label;
  
  // 熱伝達係数
  label = label_base + "/CoefOfHeatTransfer";
  
  cmp[n].set_CoefHT( get_BCval_real(label) );
  
  // 表面温度
  label = label_base + "/SurfaceTemperature";
  
  REAL_TYPE st = get_BCval_real(label);
  cmp[n].set_Temp( FBUtility::convTemp2K(st, Unit_Temp) );
  
  // 隣接セル媒質
  get_Neighbor(label_base, n, cmp, mat);
}


// #################################################################
// HeatTransfer_SNのパラメータを取得する
void ParseBC::get_IBC_HT_SN(const string label_base, const int n, CompoList* cmp, const MediumList* mat)
{
  string str;
  string label;
  
  // 表面温度
  label = label_base + "/SurfaceTemperature";
  
  REAL_TYPE st = get_BCval_real(label);
  cmp[n].set_Temp( FBUtility::convTemp2K(st, Unit_Temp) );
  
  // type
  label = label_base + "/RefTempMode";
  
  if ( !(tpCntl->GetValue(label, &str )) )
  {
    stamped_printf("\tParsing error : Invalid int value for '%s'\n", label.c_str());
    Exit(0);
  }
  
  if ( !strcasecmp("BulkTemperature", str.c_str()) )
  {
	  cmp[n].set_sw_HTmodeRef( CompoList::HT_mode_bulk );
  }
  else if ( !strcasecmp("LocalTemperature", str.c_str()) )
  {
	  cmp[n].set_sw_HTmodeRef( CompoList::HT_mode_local );
  }
  else
  {
	  stamped_printf("\tParsing error : Invalid string value for 'RefTempMode' : %s\n", str.c_str());
	  Exit(0);
  }
  
  // 隣接セル媒質
  get_Neighbor(label_base, n, cmp, mat);
  
  
  // Vertical and upper face values
  label = label_base + "/VerticalLaminarAlpha";
  cmp[n].ca[CompoList::vert_laminar_alpha]    = get_BCval_real(label);
  
  label = label_base + "/VerticalLaminarBeta";
  cmp[n].ca[CompoList::vert_laminar_beta]     = get_BCval_real(label);
  
  label = label_base + "/VerticalTurbulentAlpha";
  cmp[n].ca[CompoList::vert_turbulent_alpha]  = get_BCval_real(label);
  
  label = label_base + "/VerticalTurbulentBeta";
  cmp[n].ca[CompoList::vert_turbulent_beta]   = get_BCval_real(label);
  
  label = label_base + "/VerticalRaCritial";
  cmp[n].ca[CompoList::vert_Ra_critial]       = get_BCval_real(label);
  
  
  // Lower face values
  label = label_base + "/LowerLaminarAlpha";
  cmp[n].cb[CompoList::lower_laminar_alpha]   = get_BCval_real(label);
  
  label = label_base + "/LowerLaminarBeta";
  cmp[n].cb[CompoList::lower_laminar_beta]    = get_BCval_real(label);
  
  label = label_base + "/LowerTurbulentAlpha";
  cmp[n].cb[CompoList::lower_turbulent_alpha] = get_BCval_real(label);
  
  label = label_base + "/LowerTurbulentBeta";
  cmp[n].cb[CompoList::lower_turbulent_beta]  = get_BCval_real(label);
  
  label = label_base + "/LowerRaCritial";
  cmp[n].cb[CompoList::lower_Ra_critial]      = get_BCval_real(label);
}


// #################################################################
// HeatTransfer_SFのパラメータを取得する
void ParseBC::get_IBC_HT_SF(const string label_base, const int n, CompoList* cmp)
{
  string str;
  string label;

  // 表面温度
  label=label_base+"/SurfaceTemperature";
  REAL_TYPE st = get_BCval_real(label);
  cmp[n].set_Temp( FBUtility::convTemp2K(st, Unit_Temp) );
  
  // type
  label=label_base+"/RefTempMode";//
  if ( !(tpCntl->GetValue(label, &str )) ) {
    stamped_printf("\tParsing error : Invalid int value for 'RefTempMode' in 'LocalBoundary > HeatTransferSF'\n");
    Exit(0);
  }
  if ( !strcasecmp("BulkTemperature", str.c_str()) ) {
	  cmp[n].set_sw_HTmodeRef( CompoList::HT_mode_bulk );
  }
  else if ( !strcasecmp("LocalTemperature", str.c_str()) ) {
	  cmp[n].set_sw_HTmodeRef( CompoList::HT_mode_local );
  }
  else {
	  stamped_printf("\tParsing error : Invalid string value for 'Type' : %s\n", str.c_str());
	  Exit(0);
  }
  
  // coefficients
  label=label_base+"/Alpha";
  cmp[n].ca[CompoList::alpha] = get_BCval_real(label);
  label=label_base+"/Beta";
  cmp[n].ca[CompoList::beta]  = get_BCval_real(label);
  label=label_base+"/Gamma";
  cmp[n].ca[CompoList::gamma] = get_BCval_real(label);
}


// #################################################################
// HeatTransfer_Bのパラメータを取得する
void ParseBC::get_IBC_HT_B(const string label_base, const int n, CompoList* cmp)
{
  string label;
  
  // 熱伝達係数
  label=label_base+"/CoefOfHeatTransfer";
  cmp[n].set_CoefHT( get_BCval_real(label) );
  
  // バルク温度
  label=label_base+"/BulkTemperature";
  REAL_TYPE st = get_BCval_real(label);
  cmp[n].set_Temp( FBUtility::convTemp2K(st, Unit_Temp) );
}


// #################################################################
// Direct Forcingのパラメータを取得する
void ParseBC::get_IBC_IBM_DF(const string label_base, const int n, CompoList* cmp)
{
  int d;
  REAL_TYPE v[3];
  string str;
  string label;
  
  
  // 法線ベクトル
  get_NV(label_base, v);
  copyVec(cmp[n].nv, v);
  
  // Velocity
  label=label_base+"/Velocity";
  
  REAL_TYPE ct = get_BCval_real(label);
  
  if ( Unit_Param == DIMENSIONAL )
  {
    cmp[n].set_Velocity( ct );
  }
  else {
    cmp[n].set_Velocity( ct * RefVelocity );
  }
}



// #################################################################
// 境界条件IsoThermalのパラメータを取得し保持する
void ParseBC::get_IBC_IsoTherm(const string label_base, const int n, CompoList* cmp)
{
  string label;
  
  // 表面温度
  label=label_base+"/Temperature";
  REAL_TYPE tmp = get_BCval_real(label);
  cmp[n].set_Temp( FBUtility::convTemp2K(tmp, Unit_Temp) );
}


// #################################################################
// Monitorの設定内容をパースし，パラメータを保持する
void ParseBC::get_IBC_Monitor(const string label_base, const int n, CompoList* cmp)
{
  int nvc = 0;
  REAL_TYPE v[3];
  string str, pnt;
  string label, label_leaf;
  
  // モードと形状
  label = label_base + "/Shape";
  
  if ( !tpCntl->GetValue(label, &pnt) )
  {
    ;
  }
  else
  {
    if ( !strcasecmp("Cylinder", pnt.c_str()) )
    {
      cmp[n].set_Shape(SHAPE_CYLINDER);
    }
    else if ( !strcasecmp("Box", pnt.c_str()) )
    {
      cmp[n].set_Shape(SHAPE_BOX);
    }
    else if ( !strcasecmp("Polygon", pnt.c_str()) )
    {
      cmp[n].set_Shape(SHAPE_VOXEL);
    }
    else
    {
      Hostonly_ stamped_printf("\tParsing error : Invalid string value for 'Shape' : %s\n", pnt.c_str());
      Exit(0);
    }
  }
  
  // reference 隠しコマンドに
  label = label_base + "/Reference";
  
  if ( !(tpCntl->GetValue(label, &pnt )) )
  {
    ;
  }
  else
  {
    if ( !strcasecmp("yes", pnt.c_str()) )
    {
      cmp[n].setStateCellMonitor(ON);
    }
    else if ( !strcasecmp("no", pnt.c_str()) )
    {
      cmp[n].setStateCellMonitor(OFF);
    }
    else
    {
      Hostonly_ stamped_printf("\tParsing error : Invalid string value for 'Reference' : %s\n", pnt.c_str());
      Exit(0);
    }
  }
  
  // 法線ベクトル
  get_NV(label_base, v);
  copyVec(cmp[n].nv, v);
  
  int shp = cmp[n].get_Shape();
  
  if ( (shp == SHAPE_BOX) || (shp == SHAPE_CYLINDER) )
  {
    // 中心座標の取得
    get_Center(label_base, v);
    copyVec(cmp[n].oc, v);
  }
  
  if ( shp == SHAPE_BOX )
  {
    // 方向ベクトルの取得
    get_Dir(label_base, v);
    copyVec(cmp[n].dr, v);
    
    // 形状パラメータ
    label=label_base+"/Depth";
    cmp[n].depth  = get_BCval_real(label);
    label=label_base+"/Width";
    cmp[n].shp_p1 = get_BCval_real(label);
    label=label_base+"/Height";
    cmp[n].shp_p2 = get_BCval_real(label);
  }
  
  if ( shp == SHAPE_CYLINDER )
  {
    // 形状パラメータ
    label=label_base+"/Depth";
    cmp[n].depth  = get_BCval_real(label);
    label=label_base+"/Radius";
    cmp[n].shp_p1 = get_BCval_real(label);
  }
  
  
  // Variables
  label_leaf = label_base + "/Variables";
  if ( !(tpCntl->chkNode(label_leaf)) )
  {
    Hostonly_ stamped_printf("\tParsing error : No 'Variables' keyword in '%s'\n", label_leaf.c_str());
    Exit(0);
  }
  
  // サンプリングモード
  label = label_leaf + "/SamplingMode";
  if ( !(tpCntl->GetValue(label, &pnt )) )
  {
    ;
  }
  else
  {
    if ( !strcasecmp("all", pnt.c_str()) )
    {
      cmp[n].set_SamplingMode(SAMPLING_ALL);
    }
    else if ( !strcasecmp("fluid", pnt.c_str()) )
    {
      cmp[n].set_SamplingMode(SAMPLING_FLUID_ONLY);
    }
    else if ( !strcasecmp("solid", pnt.c_str()) )
    {
      cmp[n].set_SamplingMode(SAMPLING_SOLID_ONLY);
    }
    else
    {
      Hostonly_ stamped_printf("\tParsing error : Invalid string value for 'SamplingMode' : %s\n", pnt.c_str());
      Exit(0);
    }
  }
  
  // サンプリング方法
  label = label_leaf + "/SamplingMethod";
  if ( !(tpCntl->GetValue(label, &pnt )) )
  {
  }
  else
  {
    if ( !strcasecmp("Nearest", pnt.c_str()) )
    {
      cmp[n].set_SamplingMethod(SAMPLING_NEAREST);
    }
    else if ( !strcasecmp("Interpolation", pnt.c_str()) )
    {
      cmp[n].set_SamplingMethod(SAMPLING_INTERPOLATION);
    }
    else if ( !strcasecmp("Smoothing", pnt.c_str()) )
    {
      cmp[n].set_SamplingMethod(SAMPLING_SMOOTHING);
    }
    else
    {
      Hostonly_ stamped_printf("\tParsing error : Invalid string value for 'SamplingMethod' : %s\n", pnt.c_str());
      Exit(0);
    }
  }
  
  // モニタする変数と数を取得
  nvc = 0;
  
  // 速度
  label = label_leaf + "/Velocity";
  if ( !(tpCntl->GetValue(label, &str )) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    Exit(0);
  }
  if ( !strcasecmp(str.c_str(), "on") )
  {
    cmp[n].encodeVarType(var_Velocity);
    nvc++;
  }
  else if( !strcasecmp(str.c_str(), "off") )
  {
    ;  // nothing
  }
  else
  {
    Hostonly_ stamped_printf("\tInvalid keyword is described for '%s'\n", label.c_str());
    Exit(0);
  }
  
  // 圧力
  label = label_leaf + "/Pressure";
  if ( !(tpCntl->GetValue(label, &str )) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    Exit(0);
  }
  if ( !strcasecmp(str.c_str(), "on") )
  {
    cmp[n].encodeVarType(var_Pressure);
    nvc++;
  }
  else if( !strcasecmp(str.c_str(), "off") )
  {
    ;  // nothing
  }
  else
  {
    Hostonly_ stamped_printf("\tInvalid keyword is described for '%s'\n", label.c_str());
    Exit(0);
  }
  
  // 温度
  if ( HeatProblem )
  {
    label = label_leaf + "/Temperature";
    if ( !(tpCntl->GetValue(label, &str )) )
    {
      Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
      Exit(0);
    }
    if ( !strcasecmp(str.c_str(), "on") )
    {
      cmp[n].encodeVarType(var_Temperature);
      nvc++;
    }
    else if( !strcasecmp(str.c_str(), "off") )
    {
      ;  // nothing
    }
    else
    {
      Hostonly_ stamped_printf("\tInvalid keyword is described for '%s'\n", label.c_str());
      Exit(0);
    }
  }
  
  // 全圧
  label = label_leaf + "/TotalPressure";
  if ( !(tpCntl->GetValue(label, &str )) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    Exit(0);
  }
  if ( !strcasecmp(str.c_str(), "on") )
  {
    cmp[n].encodeVarType(var_TotalP);
    nvc++;
  }
  else if( !strcasecmp(str.c_str(), "off") )
  {
    ;  // nothing
  }
  else
  {
    Hostonly_ stamped_printf("\tInvalid keyword is described for '%s'\n", label.c_str());
    Exit(0);
  }
  
  // モニタ面に対して指定された変数の個数（モニタの個数）を取得
  cmp[n].setAttrb(nvc);
}



// #################################################################
// 内部の流出境界のパラメータを取得する
void ParseBC::get_IBC_Outflow(const string label_base, const int n, CompoList* cmp)
{
  int def;
  REAL_TYPE ct;
  REAL_TYPE v[3];
  string str;
  string label;
  
  // 圧力境界のタイプ default
  cmp[n].set_P_BCtype( P_GRAD_ZERO );
  
  // Hidden parameter
  label = label_base + "/PressureType";
  
  if ( !(tpCntl->GetValue(label, &str )) ) {
    ; // なくてもOK
  }
  else {
    if ( !strcasecmp("Dirichlet", str.c_str()) ) {
      cmp[n].set_P_BCtype( P_DIRICHLET );
    }
    else if ( !strcasecmp("GradZero", str.c_str()) ) {
      cmp[n].set_P_BCtype( P_GRAD_ZERO );
    }
    else {
      printf("\tParsing error : Invalid string value for 'PressureType' : %s\n", str.c_str());
      Exit(0);
    }
    if ( cmp[n].get_P_BCtype() == P_DIRICHLET ) {
      label = label_base + "/PressureValue";
      cmp[n].set_Pressure( get_BCval_real(label) );
    }
  }
  
  // 流出速度のタイプ
  cmp[n].setOutflowType(V_AVERAGE);
  
  // 法線ベクトルの取得
  get_NV(label_base, v);
  copyVec(cmp[n].nv, v);
  
  /*
   label=label_base+"/velocity_type";//
   if ( !(tpCntl->GetValue(label, &str )) ) {
   printf("\tParsing error : fail to get 'Velocity_Type' in 'LocalBoundary > Outflow'\n");
   Exit(0);
   }
   if ( !strcasecmp("average", str.c_str()) ) {
   cmp[n].flag = V_AVERAGE;
   }
   else if ( !strcasecmp("minmax", str.c_str()) ) {
   cmp[n].flag = V_MINMAX;
   }
   else {
   printf("\tParsing error : Invalid string value for 'Velocity_Type' : %s\n", str);
   Exit(0);
   }
   */
  
}



// #################################################################
// 内部の周期境界のパラメータを取得する
void ParseBC::get_IBC_Periodic(const string label_base, const int n, CompoList* cmp)
{
  int dir=0;
  REAL_TYPE ct=0.0;
  string str;
  string label;
  
  // 上流側の方向
  label = label_base + "/UpstreamDirection";
  
  if ( !(tpCntl->GetValue(label, &str )) )
  {
	  printf("\tParsing error : fail to get 'upstream_direction' in 'LocalBoundary > Periodic'\n");
	  Exit(0);
  }
  if ( !strcasecmp("xminus", str.c_str()) )
  {
		dir = X_MINUS;
	}
  else if ( !strcasecmp("xplus", str.c_str()) )
  {
    dir = X_PLUS;
  }
  else if ( !strcasecmp("yminus", str.c_str()) )
  {
    dir = Y_MINUS;
  }
  else if ( !strcasecmp("yplus", str.c_str()) )
  {
    dir = Y_PLUS;
  }
  else if ( !strcasecmp("zminus", str.c_str()) )
  {
    dir = Z_MINUS;
  }
  else if ( !strcasecmp("zplus", str.c_str()) )
  {
    dir = Z_PLUS;
  }
  else
  {
    stamped_printf("\tParsing error : Invalid direction in 'LocalBoundary > Periodic'\n");
    Exit(0);
  }
	cmp[n].setPeriodicDir((int)dir);
  
  
  // 圧力差
  label = label_base + "/PressureDifference";
  
  if ( !(tpCntl->GetValue(label, &ct )) )
  {
    printf("\tParsing error : Invalid value of 'Pressure difference' in 'LocalBoundary > Periodic'\n");
    Exit(0);
  }
  else
  {
    cmp[n].ca[0] = ct;
  }
}


// #################################################################
// HeatExchangerのパラメータを取得する
///> @note この時点ではRefDensityの値が未定なので，あとでパラメータ処理
///> @see Control::setParameters()
void ParseBC::get_IBC_PrsLoss(const string label_base, const int n, CompoList* cmp)
{
  string str,str_u;
  string label;
  REAL_TYPE v[4], ct;
  
  
  // 入力単位の指定
  label = label_base + "/Unit";
  
  if ( !(tpCntl->GetValue(label, &str_u )) ) {
		stamped_printf("\tParsing error : Invalid float value for 'Unit' in 'PressureLoss'\n");
		Exit(0);
  }
  if ( !strcasecmp("mmaq", str_u.c_str()) ) {
    cmp[n].setPrsUnit(CompoList::unit_mmAq);
  }
  else if ( !strcasecmp("mmhg", str_u.c_str()) ) {
    cmp[n].setPrsUnit(CompoList::unit_mmHg);
  }
  else if ( !strcasecmp("pa", str_u.c_str()) ) {
    cmp[n].setPrsUnit(CompoList::unit_Pa);
  }
  else if ( !strcasecmp("NonDimension", str_u.c_str()) ) {
    cmp[n].setPrsUnit(CompoList::unit_NonDimensional);
  }
  else {
    stamped_printf("\tDescribed unit is out of scope.\n");
    Exit(0);
  }
  
  // 法線ベクトルの取得
  get_NV(label_base, v);
  copyVec(cmp[n].nv, v);
  
  // 方向ベクトルの取得
  get_Dir(label_base, v);
  copyVec(cmp[n].dr, v);
  
  // 中心座標の取得
  get_Center(label_base, v);
  copyVec(cmp[n].oc, v);
  
  // 形状パラメータ
  label = label_base + "/Depth";
  cmp[n].depth  = get_BCval_real(label);
  
  label = label_base + "/Width";
  cmp[n].shp_p1 = get_BCval_real(label);
  
  label = label_base + "/Height";
  cmp[n].shp_p2 = get_BCval_real(label);
  
  // 圧力損失パラメータ
  label = label_base + "/C";
  if( !(tpCntl->GetVector(label, v, 4)) )
  {
    stamped_printf("\tParsing error : fail to get vec params in '%s\n", label.c_str());
    Exit(0);
  }
  cmp[n].ca[0]=v[0];
  cmp[n].ca[1]=v[1];
  cmp[n].ca[2]=v[2];
  cmp[n].ca[3]=v[3];
  
  label = label_base + "/Uthreshold";
  cmp[n].ca[4] = get_BCval_real(label);
  
  label = label_base + "/Thickness";
  cmp[n].ca[5] = get_BCval_real(label);
  
  // 熱交換器の方向強制オプション
  
  label = label_base + "/Vector";
  if ( !(tpCntl->GetValue(label, &str )) ) {
    stamped_printf("\tParsing error : Invalid string for 'vector' in 'PressureLoss'\n");
    Exit(0);
  }
  if ( !strcasecmp("directional", str.c_str()) ) {
    cmp[n].set_sw_HexDir( ON );
  }
  else {
    cmp[n].set_sw_HexDir( OFF );
  }
  
}



// #################################################################
// 境界条件Radiantのパラメータを取得し保持する
// 境界条件自体は未実装
void ParseBC::get_IBC_Radiant(const string label_base, const int n, CompoList* cmp)
{
  string label;
  
  // 係数
  label=label_base+"/Epsilon";
  cmp[n].set_CoefRadEps( get_BCval_real(label) );
  
  // 射出率
  label=label_base+"/Projection";
  cmp[n].set_CoefRadPrj( get_BCval_real(label) );
}



// #################################################################
// 内部の流入境界のパラメータを取得する
// Control::setparameters()でcmp[].ca[]に値をセットする
void ParseBC::get_IBC_SpecVel(const string label_base, const int n, CompoList* cmp)
{
  string str;
  string label;
  REAL_TYPE ct;
  REAL_TYPE v[3];
  
  // 指定タイプの特定
  label = label_base + "/Type";
  
  if ( !(tpCntl->GetValue(label, &str )) )
  {
    Hostonly_ stamped_printf("\tParsing error : Invalid SpecifiedType in '%s'\n", label.c_str());
    Exit(0);
  }
  if ( !strcasecmp("Velocity", str.c_str()) )
  {
	  cmp[n].set_VBC_policy(true);
  }
  else if ( !strcasecmp("Massflow", str.c_str()) )
  {
	  cmp[n].set_VBC_policy(false);
  }
  else
  {
	  Hostonly_ printf("\tParsing error : Invalid string value '%s' for 'Type'\n", str.c_str());
	  Exit(0);
  }
  
  // 速度指定タイプ
  cmp[n].set_V_profile( get_Vel_profile(label_base) );
  
  
  // 速度パラメータの読み込み
  get_Vel_Params(label_base, cmp[n].get_V_Profile(), cmp[n].ca, "LocalBoundary", cmp[n].isPolicy_Massflow());

  
  // 法線ベクトル
  get_NV(label_base, v);
  copyVec(cmp[n].nv, v);
  

  
  // 境界条件の方向
  label = label_base + "/FluidDirection";
  
  if ( !(tpCntl->GetValue(label, &str )) )
  {
    stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    Exit(0);
  }
  if ( !strcasecmp("SameSideOfNormal", str.c_str()) )
  {
    cmp[n].setBClocation(CompoList::same_direction);
  }
  else if ( !strcasecmp("OppositeSideOfNormal", str.c_str()) )
  {
    cmp[n].setBClocation(CompoList::opposite_direction);
  }
  else
  {
    printf("\tParsing error : Invalid string value '%s' for 'FluidDirection'\n", str.c_str());
    Exit(0);
  }
  
  
  // heat problem
  if ( HeatProblem )
  {
    label = label_base + "/Temperature";
    
    if ( !(tpCntl->GetValue(label, &ct )) )
    {
      Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
      Exit(0);
    }

    cmp[n].set_Temp( FBUtility::convTemp2K(ct, Unit_Temp) );
    
    if ( Unit_Param != DIMENSIONAL )
    {
      Hostonly_ stamped_printf("\tWarning: Heat condition must be a dimensional value\n");
      Exit(0);
    }
    cmp[n].setType(SPEC_VEL_WH); // SPEC_VELから変更

  }
  
}


// #################################################################
/**
 * @brief 温度計算の場合の各媒質の初期値を取得する
 * @param [in,out] cmp    CompoList
 */
void ParseBC::get_Medium_InitTemp(CompoList* cmp)
{  
  string label, label_base;
  string str;
  
  label_base = "/Parameter/InitTempOfMedium";
  
  if ( !tpCntl->chkNode(label_base) )
  {
    Hostonly_ stamped_printf("\tParsing error : Missing the section of '%s'\n", label_base.c_str());
    Exit(0);
  }
  
  int m_no_medium = NoCompo - NoBC;
  
  
  for (int i=1; i<=m_no_medium; i++) {
    
    if ( !tpCntl->GetNodeStr(label_base, i, &str) )
    {
      Hostonly_ stamped_printf("\tGetNodeStr error\n");
      Exit(0);
    }
    
    // medium name
    label = label_base + "/" + str;
    REAL_TYPE ct;
    
    if ( !tpCntl->GetValue(label, &ct) )
    {
      Hostonly_ stamped_printf("\tParsing error : Invalid keyword in '%s'\n", label_base.c_str());
      Exit(0);
    }
    
    int flag = 0;
    
    for (int m=NoBC+1; m<=NoCompo; m++) {
      if ( FBUtility::compare(cmp[m].getLabel(), str) )
      {
        cmp[m].setInitTemp( FBUtility::convTemp2K(ct, Unit_Temp) );
        flag++;
      }
    }
    
    if ( flag == 0 )
    {
      Hostonly_ stamped_printf("\tParsing error : Invalid keyword in '%s'\n", label_base.c_str());
      Exit(0);
    }
    
  }
}


// #################################################################
/**
 * @brief 隣接セルラベルを取得する
 * @param [in]  label_base ラベルディレクトリ
 * @param [in]  n          コンポーネントリストの格納番号
 * @param [out] cmp        CompoList
 * @param [in]  mat        MediumList
 */
void ParseBC::get_Neighbor(const string label_base, const int n, CompoList* cmp, const MediumList* mat)
{
  string label = label_base + "/NeighborMedium";
  string str;
  
  if ( !tpCntl->GetValue(label, &str) )
  {
    Hostonly_ stamped_printf("\tParsing error : Invalid keyword in '%s'\n", label_base.c_str());
    Exit(0);
  }
  
  bool flag = false;
  for (int i=1; i<=NoMedium; i++) {
    if ( FBUtility::compare(str, mat[i].getLabel()) )
    {
      cmp[n].setDef(i);
      flag = true;
      break;
    }
  }
  
  if ( !flag )
  {
    Hostonly_ stamped_printf("\tError : Keyword '%s' is not listed on MediumTable\n", str.c_str());
    Exit(0);
  }
  
}

// #################################################################
/**
 * @brief 内部境界条件の法線ベクトル値を取得し，登録する
 * @param [in] label_base ラベルディレクトリ
 * @param [out v          ベクトルパラメータ
 */
void ParseBC::get_NV(const string label_base, REAL_TYPE* v)
{
  string label;
  for (int i=0; i<3; i++) v[i]=0.0f;
  
  label = label_base + "/Normal";
  
  if( !(tpCntl->GetVector(label, v, 3)) )
  {
    stamped_printf("\tParsing error : fail to get vec params in '%s\n", label.c_str());
    Exit(0);
  }
  
  //単位ベクトル化
  getUnitVec(v);
}


// #################################################################
/**
 * @brief Originのキーワードに対する文字列をパースし，返す
 * @param [in] label_base ラベルディレクトリ
 */
string ParseBC::get_Origin(const string label_base)
{
  string label = label_base + "/Origin";
  string str;
  
  if ( !tpCntl->GetValue(label, &str) )
  {
    Hostonly_ stamped_printf("\tParsing error : Invalid keyword of '%s'\n", label.c_str());
    Exit(0);
  }
  
  return str;
}


// #################################################################
/**
 * @brief 外部境界の遠方境界のパラメータを取得する
 * @param [in] label_base ラベルディレクトリ
 * @param [in] n          面番号
 */
void ParseBC::get_OBC_FarField(const string label_base, const int n)
{
  BaseBc[n].set_pType(P_GRAD_ZERO);
  BaseBc[n].p = 0.0; // ダミー値
  
  // 外部雰囲気温
  if ( HeatProblem )
  {
    if ( Unit_Param != DIMENSIONAL )
    {
      stamped_printf("\tError: Heat condition must be given by dimensional value\n");
      Exit(0);
    }
    
    string label = label_base + "/AmbientTemperature";
    REAL_TYPE ct;
    
    if ( !(tpCntl->GetValue(label, &ct )) )
    {
      stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
      Exit(0);
    }

    BaseBc[n].set_Temp( FBUtility::convTemp2K(ct, Unit_Temp) );
  }
  
}


// #################################################################
/*
 @fn void ParseBC::get_OBC_HT(const string label_base, const int n, const string kind)
 @brief 外部の壁面熱伝達境界のパラメータを取得する
 @param label_base
 @param n 面番号
 @param kind 熱伝達境界の種類
 */
void ParseBC::get_OBC_HT(const string label_base, const int n, const string kind)
{
  string label;
  string str;
  REAL_TYPE ct;
  
  if ( !strcasecmp(kind.c_str(), "HeatTransferB") )
  {
    BaseBc[n].set_HTmode(HT_B);
    label = label_base + "/CoefOfHeatTransfer";
    BaseBc[n].set_CoefHT( get_BCval_real(label) );
    
    label = label_base + "/BulkTemperature";
    ct = get_BCval_real(label);
    BaseBc[n].set_Temp( FBUtility::convTemp2K(ct, Unit_Temp) );
  }
  else if ( !strcasecmp(kind.c_str(), "HeatTransferN") )
  {
    BaseBc[n].set_HTmode(HT_N);
    label = label_base + "/CoefOfHeatTransfer";
    BaseBc[n].set_CoefHT( get_BCval_real(label) );
  }
  else if ( !strcasecmp(kind.c_str(), "HeatTransferS") )
  {
    BaseBc[n].set_HTmode(HT_S);
    label = label_base + "/CoefOfHeatTransfer";
    BaseBc[n].set_CoefHT( get_BCval_real(label) );
    
    label = label_base + "/SurfaceTemperature";
    ct = get_BCval_real(label);
    BaseBc[n].set_Temp( FBUtility::convTemp2K(ct, Unit_Temp) );
  }
  else if ( !strcasecmp(kind.c_str(), "HeatTransferSF") )
  {
    BaseBc[n].set_HTmode(HT_SF);
    label = label_base + "/SurfaceTemperature";
    BaseBc[n].set_Temp( get_BCval_real(label) );
    
    
    label=label_base+"/RefTempMode";
    if ( !(tpCntl->GetValue(label, &str )) )
    {
      stamped_printf("\tParsing error : Invalid int value for '%s'\n", label.c_str());
      Exit(0);
    }
    if ( !strcasecmp("BulkTemperature", str.c_str()) )
    {
      BaseBc[n].set_HTmodeRef( CompoList::HT_mode_bulk );
    }
    else if ( !strcasecmp("LocalTemperature", str.c_str()) )
    {
      BaseBc[n].set_HTmodeRef( CompoList::HT_mode_local );
    }
    else
    {
      stamped_printf("\tParsing error : Invalid string value for '%s'\n", label.c_str());
      Exit(0);
    }
    // coefficients
    label = label_base + "/Alpha";
    BaseBc[n].ca[0] = get_BCval_real(label);
    
    label = label_base + "/Beta";
    BaseBc[n].ca[1] = get_BCval_real(label);
    
    label = label_base + "/Gamma";
    BaseBc[n].ca[2] = get_BCval_real(label);
  }
  else if ( !strcasecmp(kind.c_str(), "HeatTransferSN") ) {
    BaseBc[n].set_HTmode(HT_SN);
    label = label_base + "/SurfaceTemperature";
    BaseBc[n].set_Temp( get_BCval_real(label) );
    
    // reference mode
    label = label_base + "/RefTempMode";
    if ( !(tpCntl->GetValue(label, &str )) )
    {
      stamped_printf("\tParsing error : Invalid int value for '%s'\n", label.c_str());
      Exit(0);
    }
    if ( !strcasecmp("BulkTemperature", str.c_str()) )
    {
      BaseBc[n].set_HTmodeRef( CompoList::HT_mode_bulk );
    }
    else if ( !strcasecmp("LocalTemperature", str.c_str()) )
    {
      BaseBc[n].set_HTmodeRef( CompoList::HT_mode_local );
    }
    else
    {
      stamped_printf("\tParsing error : Invalid string value for '%s'\n", label.c_str());
      Exit(0);
    }
    
    // Vertical and upper face values
    label = label_base + "/VertivalLaminarAlpha";
    BaseBc[n].ca[0] = get_BCval_real(label);
    
    label = label_base + "/VertivalLaminarBeta";
    BaseBc[n].ca[1] = get_BCval_real(label);
    
    label = label_base + "/vertivalTurbulentAlpha";
    BaseBc[n].ca[2] = get_BCval_real(label);
    
    label = label_base + "/VertivalTurbulentBeta";
    BaseBc[n].ca[3] = get_BCval_real(label);
    
    label = label_base + "/VertivalRaCritial";
    BaseBc[n].ca[4] = get_BCval_real(label);
    
    // Lower face values
    label = label_base + "/LowerLaminarAlpha";
    BaseBc[n].cb[0] = get_BCval_real(label);
    
    label = label_base + "/LowerLaminarBeta";
    BaseBc[n].cb[1] = get_BCval_real(label);
    
    label = label_base + "/LowerTurbulentAlpha";
    BaseBc[n].cb[2] = get_BCval_real(label);
    
    label = label_base + "/LowerTurbulentBeta";
    BaseBc[n].cb[3] = get_BCval_real(label);
    
    label = label_base + "/LowerRaCritial";
    BaseBc[n].cb[4] = get_BCval_real(label);
  }
  else {
    stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    Exit(0);
  }
}



// #################################################################
/**
 * @brief 外部境界の流出条件のパラメータを取得する
 * @param [in] label_base ラベルディレクトリ
 * @param [in] n          面番号
 * @note 圧力の値は，Control::setParameters()で無次元化する
 */
void ParseBC::get_OBC_Outflow(const string label_base, const int n)
{
  string str;
  string label;
  REAL_TYPE ct;
  
  
  // 流出速度のタイプ
  label = label_base + "/VelocityType";
  
  if ( !(tpCntl->GetValue(label, &str )) )
  {
    stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    Exit(0);
  }
  
  if ( !strcasecmp("Average", str.c_str()) )
  {
	  BaseBc[n].set_ofv(V_AVERAGE);
  }
  else if ( !strcasecmp("Minmax", str.c_str()) )
  {
	  BaseBc[n].set_ofv(V_MINMAX);
  }
  else
  {
	  stamped_printf("\tParsing error : Invalid string value for 'VelocityType' : %s\n", str.c_str());
	  Exit(0);
  }
  
  
  // 圧力境界のタイプ  default
  BaseBc[n].set_pType(P_GRAD_ZERO);
  BaseBc[n].p = 0.0; // ダミー値
  
  // Hidden option
  label = label_base + "/PressureType";
  
  if ( !(tpCntl->GetValue(label, &str )) )
  {
    ; // エラーではない
  }
  else
  {
    if ( !strcasecmp("dirichlet", str.c_str()) )
    {
      BaseBc[n].set_pType(P_DIRICHLET);
    }
    else if ( !strcasecmp("gradzero", str.c_str()) )
    {
      BaseBc[n].set_pType(P_GRAD_ZERO);
    }
    else
    {
      stamped_printf("\tParsing error : Invalid string value for 'PressureType' : %s\n", str.c_str());
      Exit(0);
    }
  }
  
  
  // 圧力の値
  if ( BaseBc[n].get_pType() == P_DIRICHLET )
  {
    if ( !strcasecmp("dirichlet", str.c_str()) )
    {
      label = label_base + "/PrsValue";
      
      if ( !(tpCntl->GetValue(label, &ct )) )
      {
        stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
        Exit(0);
      }
      else
      {
        BaseBc[n].p = ct;
      }
    }
  }
}


// #################################################################
/**
 * @brief 外部境界の周期条件のパラメータを取得する
 * @param [in] label_base ラベル
 * @param [in] n          面番号
 * @note 圧力の値は，Control::setParameters()で無次元化する
 */
void ParseBC::get_OBC_Periodic(const string label_base, const int n)
{
  REAL_TYPE ct;
  int def;
  string str;
  string label;
  
  // モード
  label = label_base + "/Mode";
  if ( !(tpCntl->GetValue(label, &str )) )
  {
    printf("\tParsing error : fail to get '%s'\n", label.c_str());
    Exit(0);
  }
  else
  {
    if ( !strcasecmp(str.c_str(), "SimpleCopy") )
    {
      BaseBc[n].set_PrdcMode(BoundaryOuter::prdc_Simple);
    }
    else if ( !strcasecmp(str.c_str(), "directional") )
    {
      BaseBc[n].set_PrdcMode(BoundaryOuter::prdc_Directional);
    }
    else if ( !strcasecmp(str.c_str(), "driver") )
    {
      BaseBc[n].set_PrdcMode(BoundaryOuter::prdc_Driver);
    }
  }
  
  // Directional
  if ( BaseBc[n].get_PrdcMode() == BoundaryOuter::prdc_Directional )
  {
    label = label_base + "/PressureDifference";
    if ( !(tpCntl->GetValue(label, &ct )) ) {
      printf("\tParsing error : fail to get '%s'\n", label.c_str());
      Exit(0);
    }
    else {
      BaseBc[n].p = ct;
    }
    
    label = label_base + "/FlowDirection";
    if ( !(tpCntl->GetValue(label, &str )) )
    {
      printf("\tParsing error : fail to get '%s'\n", label.c_str());
      Exit(0);
    }
    else
    {
      if ( !strcasecmp(str.c_str(), "Upstream") )
      {
        BaseBc[n].set_FaceMode(BoundaryOuter::prdc_upstream);
      }
      else if ( !strcasecmp(str.c_str(), "Downstream") ) {
        BaseBc[n].set_FaceMode(BoundaryOuter::prdc_downstream);
      }
      else
      {
        printf("\tParsing error : Invalid keyword in '%s'\n", label.c_str());
        Exit(0);
      }
    }
  }
  
  // Driver
  if ( BaseBc[n].get_PrdcMode() == BoundaryOuter::prdc_Driver )
  {
    
    label = label_base + "/DriverDirection";
    if ( !(tpCntl->GetValue(label, &str )) )
    {
      printf("\tParsing error : fail to get '%s'\n", label.c_str());
      Exit(0);
    }
    else
    {
      if ( !strcasecmp("xminus", str.c_str()) )
      {
        def = X_MINUS;
      }
      else if ( !strcasecmp("xplus", str.c_str()) )
      {
        def = X_PLUS;
      }
      else if ( !strcasecmp("yminus", str.c_str()) )
      {
        def = Y_MINUS;
      }
      else if ( !strcasecmp("yplus", str.c_str()) )
      {
        def = Y_PLUS;
      }
      else if ( !strcasecmp("zminus", str.c_str()) )
      {
        def = Z_MINUS;
      }
      else if ( !strcasecmp("zplus", str.c_str()) )
      {
        def = Z_PLUS;
      }
      else {
        printf("\tParsing error : Invalid keyword in '%s'\n", label.c_str());
        Exit(0);
      }
      BaseBc[n].set_DriverDir(def);
    }
    
    
    label = label_base + "/DriverLidIndex";
    if ( !(tpCntl->GetValue(label, &def )) ) {
      printf("\tParsing error : fail to get '%s'\n", label.c_str());
      Exit(0);
    }
    else {
      BaseBc[n].set_DriverIndex(def);
    }
  }
}



// #################################################################
/**
 * @brief 外部境界の流入条件のパラメータを取得する
 * @param [in] label_base ラベルディレクトリ
 * @param [in] n          面番号
 */
void ParseBC::get_OBC_SpecVH(const string label_base, const int n)
{
  
  // 速度境界条件のプロファイル
  BaseBc[n].set_V_Profile( get_Vel_profile(label_base) );
  
  // 法線ベクトル
  if ( BaseBc[n].get_V_Profile() != CompoList::vel_zero )
  {
    REAL_TYPE v[3];
    get_NV(label_base, v);
    BaseBc[n].addVec(v);
  }
  
  // 速度のパラメータ読み込み
  get_Vel_Params(label_base, BaseBc[n].get_V_Profile(), BaseBc[n].ca, "OuterBoundary");

  
  // heat problem
  if ( HeatProblem )
  {
    if ( Unit_Param != DIMENSIONAL )
    {
      stamped_printf("\tError: Heat condition must be given by dimensional value\n");
      Exit(0);
    }
    
    string label = label_base + "/Temperature";
    REAL_TYPE ct;
    
    if ( !(tpCntl->GetValue(label, &ct )) )
    {
      stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
      Exit(0);
    }
    else
    {
      BaseBc[n].set_Temp( FBUtility::convTemp2K(ct, Unit_Temp) );
      BaseBc[n].set_hType(CNST_TEMP);
    }
  }
}


// #################################################################
/**
 * @brief トラクションフリーの外部境界のパラメータを取得する
 * @param [in] label_base ラベルディレクトリ
 * @param [in] n          面番号
 */
void ParseBC::get_OBC_Trcfree(const string label_base, const int n)
{
  BaseBc[n].set_pType(P_DIRICHLET);
  BaseBc[n].p = 0.0; // ゲージ圧zero 固定
  
  // 外部雰囲気温
  if ( HeatProblem )
  {
    if ( Unit_Param != DIMENSIONAL )
    {
      stamped_printf("\tError: Heat condition must be given by dimensional value\n");
      Exit(0);
    }
      
    string label = label_base + "/AmbientTemperature";
    REAL_TYPE ct;
    
    if ( !(tpCntl->GetValue(label, &ct )) )
    {
      stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
      Exit(0);
    }

    BaseBc[n].set_Temp( FBUtility::convTemp2K(ct, Unit_Temp) );
  }
}


// #################################################################
/**
 * @brief 外部境界の壁面条件のパラメータを取得する
 * @param [in] label_base ラベルディレクトリ
 * @param [in] n          面番号
 */
void ParseBC::get_OBC_Wall(const string label_base, const int n)
{
  REAL_TYPE vel, ct;
  REAL_TYPE v[3];
  string str;
  string label, label2;
  
  // 速度のタイプの特定
  label = label_base + "/Type";
  
  if ( !(tpCntl->GetValue(label, &str )) )
  {
    stamped_printf("\tParsing error : Invalid Specified_Type in '%s'\n", label.c_str());
    Exit(0);
  }
  
  if ( !strcasecmp("fixed", str.c_str()) )
  {
	  BaseBc[n].set_wallType(BoundaryOuter::fixed);
    BaseBc[n].set_V_Profile( CompoList::vel_zero );
  }
  else if ( !strcasecmp("slide", str.c_str()) )
  {
	  BaseBc[n].set_wallType(BoundaryOuter::slide);
    BaseBc[n].set_V_Profile( get_Vel_profile(label_base) );
  }
  else
  {
	  printf("\tParsing error : Invalid string value '%s' for 'Type'\n", str.c_str());
	  Exit(0);
  }
  
  // 法線ベクトル
  if ( BaseBc[n].get_V_Profile() != CompoList::vel_zero )
  {
    get_NV(label_base, v);
    BaseBc[n].addVec(v);
  }
  
  // 速度のパラメータ読み込み
  get_Vel_Params(label_base, BaseBc[n].get_V_Profile(), BaseBc[n].ca, "OuterBoundary");
  
  
  
  // heat problem
  if ( HeatProblem )
  {
    label = label_base + "/ThermalOption";
    
    if ( !(tpCntl->GetValue(label, &str )) )
    {
      stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
      Exit(0);
    }
    if ( Unit_Param != DIMENSIONAL )
    {
      stamped_printf("\tError: Heat condition must be given by dimensional value\n");
      Exit(0);
    }
    
    if ( !strcasecmp(str.c_str(), "Adiabatic") )
    {
      BaseBc[n].set_hType(ADIABATIC);
      BaseBc[n].p = 0.0;
    }
    else if( !strcasecmp(str.substr(0,12).c_str(), "HeatTransfer") )
    {
      BaseBc[n].set_hType(TRANSFER);
      get_OBC_HT(label, n, str);
    }
    else if( !strcasecmp(str.c_str(), "HeatFlux") )
    {
      BaseBc[n].set_hType(HEATFLUX);
      label2 = label + "/Flux";
      BaseBc[n].set_Heatflux( get_BCval_real(label2) ); // 正符号は流入
    }
    else if( !strcasecmp(str.c_str(), "Isothermal") )
    {
      BaseBc[n].set_hType(ISOTHERMAL);
      label2 = label + "/Temperature";
      ct = get_BCval_real(label2); // 表面温度
      BaseBc[n].set_Temp( FBUtility::convTemp2K(ct, Unit_Temp) );
    }
    else if( !strcasecmp(str.c_str(), "ConstantTemperature") )
    {
      BaseBc[n].set_hType(CNST_TEMP);
      label = label_base + "/Temperature";
      ct = get_BCval_real(label); // 指定温度
      BaseBc[n].set_Temp( FBUtility::convTemp2K(ct, Unit_Temp) );
    }
    else
    {
      stamped_printf("\tParsing error : Invalid string value for 'ThermalOption' : %s\n", str.c_str());
      Exit(0);
    }
  }
}



// #################################################################
/**
 * @brief 2相流問題で気相か液相かを取得する
 * @param [out] cmp   CompoList
 */
void ParseBC::get_Phase(CompoList* cmp)
{
  int m_phase;
  int id;
  string str,p;
  string label,label_base;
  int NoParam;
  
  ///////////////////////////////////////////////////////////////////////////////
  stamped_printf("\tWARNING not yet\n");
  Exit(0);
  ///////////////////////////////////////////////////////////////////////////////
  
  
  label_base="/Steer/PhaseIdetification";
  //cout <<  "label : " << label << endl;
  if ( !(tpCntl->chkNode(label_base)) ) {
    stamped_printf("\tParsing error : Missing the section of 'PhaseIdetification'\n");
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
     if ( cmp[n].getID() == id ) {
     if ( cmp[n].getState() == FLUID ) cmp[n].setPhase(m_phase);
     }
     }
     */
    
  }
  
  // check Phase of Fluid
  int tmp;
  bool sw=true;
  for (int n=1; n<=NoMedium; n++) {
    if ( cmp[n].getState() == FLUID ) {
      tmp = cmp[n].getPhase();
      if ( (tmp!=GAS) && (tmp!=LIQUID) ) {
        stamped_printf("\tcomponent [%d] is fluid, but not identified by gas or liquid.\n", n);
        sw = false;
      }
    }
  }
  if ( sw == false ) Exit(0);
}



// #################################################################
//@brief 外部境界の速度境界条件のタイプを取得し，返す
int ParseBC::get_Vel_profile(const string label_base)
{
  string label, str;
  
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


// #################################################################
/**
 @brief 速度のパラメータを取得する
 @param label_base 
 @param ca 係数パラメータの配列
 @param policy 流量指定(true) or 速度指定(false) > 外部境界は速度のみ
 @note 
 - 値は，Control::setParameters()で無次元化する
 - 速度プロファイルは単振動と一定値の場合で係数の保持パターンが異なる
 - 内部境界の場合には，流量指定と速度指定があるので分岐処理（未実装）
 */
void ParseBC::get_Vel_Params(const string label_base, const int prof, REAL_TYPE* ca, const char* str, const bool policy)
{
  string label;
  REAL_TYPE ct=0.0, vel;
  
  if ( prof == CompoList::vel_zero )
  {
    ca[CompoList::amplitude] = 0.0;
    ca[CompoList::frequency] = 0.0;
    ca[CompoList::initphase] = 0.0;
    ca[CompoList::bias]      = 0.0;
  }
  else if ( prof == CompoList::vel_constant )
  {
    label = label_base + "/Velocity";
    
    if ( !(tpCntl->GetValue(label, &ct )) )
    {
      stamped_printf("\tParsing error : Invalid value in '%s > Velocity'\n", str);
      Exit(0);
    }
    if ( policy ) // 流量指定の場合 LocalBoundaryの場合のみ
    {
      vel = ct; // 有次元でも無次元でも，モデルの断面積を計算して，後ほどパラメータ計算 >> Control::setParameters()
    }
    else
    {
      vel = ( Unit_Param == DIMENSIONAL ) ? ct : ct * RefVelocity;
    }
    
    ca[CompoList::amplitude] = 0.0;
    ca[CompoList::frequency] = 0.0;
    ca[CompoList::initphase] = 0.0;
    ca[CompoList::bias]      = vel;
  }
  else if ( prof == CompoList::vel_harmonic)
  {
    label = label_base + "/Amplitude";
    
    if ( !(tpCntl->GetValue(label, &ct )) )
    {
      stamped_printf("\tParsing error : fail to get 'Amplitude' in '%s'\n", label_base.c_str());
      Exit(0);
    }
    ca[CompoList::amplitude] = ct;
    
    label = label_base + "/Frequency";
    
    if ( !(tpCntl->GetValue(label, &ct )) )
    {
      stamped_printf("\tParsing error : fail to get 'Frequency' in '%s'\n", label_base.c_str());
      Exit(0);
    }
    ca[CompoList::frequency] = ct;
    
    label = label_base + "/InitialPhase";
    
    if ( !(tpCntl->GetValue(label, &ct )) )
    {
      stamped_printf("\tParsing error : fail to get 'InitialPhase' in '%s'\n", label_base.c_str());
      Exit(0);
    }
    ca[CompoList::initphase] = ct;
    
    
    label = label_base + "/ConstantBias";
    
    if ( !(tpCntl->GetValue(label, &ct )) )
    {
      stamped_printf("\tParsing error : fail to get 'ConstantBias' in '%s'\n", label_base.c_str());
      Exit(0);
    }
    ca[CompoList::bias] = ct;
  }
  else
  {
    Exit(0);
  }
  
}


// #################################################################
/**
 * @brief TPのポインタを受け取る
 * @param [in] tp  TPControlクラスのポインタ
 */
void ParseBC::importTP(TPControl* tp)
{ 
  if ( !tp ) Exit(0);
  tpCntl = tp;
}


// #################################################################
/**
 * @brief コンポーネントが存在するかどうかを調べる
 * @param [in] label  テストするラベル
 * @param [in] cmp    CompoList
 * @retval bool値
 */
bool ParseBC::isComponent(const int label, const CompoList* cmp)
{
  for (int n=1; n<=NoBC; n++) {
    if ( cmp[n].getType() == label ) return true;
  }
  return false;
}


// #################################################################
/**
 * @brief HTコンポーネントが存在するかどうかを調べる
 * @param [in] label  テストするラベル
 * @param [in] cmp    CompoList
 * @retval bool値
 */
bool ParseBC::isCompoTransfer(const int label, const CompoList* cmp)
{
  for (int n=1; n<=NoBC; n++) {
    if ( cmp[n].getHtype() == label ) return true;
  }
  return false;
}


// #################################################################
/**
 * @brief 同じラベルが既にコンポーネントに登録されているかを調べる
 * @retval 重複していればfalseを返す
 * @param [in] candidate テストするラベル
 * @param [in] now       コンポーネントリストの現在までのエントリ番号
 * @param [in] cmp       CompoList
 */
bool ParseBC::isLabelinCompo(const string candidate, const int now, const CompoList* cmp)
{
  for (int i=1; i<now; i++) {
    if ( FBUtility::compare(candidate, cmp[i].getLabel()) ) return false;
  }
  return true;
}



// #################################################################
/**
 * @brief CompoListに内部境界条件の情報を設定する
 * @param [in]  C     Control
 * @param [in]  mat   MediumList
 * @param [out] cmp   CompoList
 * @param [in]  polyP ポリゴン管理構造体
 * @note 最初にBCの情報を登録，その後IDの情報を登録
 * @note パラメータファイルから各内部BCのidをパースし，cmpに保持する
 * @note 格納番号は1からスタート
 */
void ParseBC::loadBC_Local(Control* C, const MediumList* mat, CompoList* cmp, Control::Polygon_property* polyP)
{ 
  string str, label;
  string label_base, label_ename, label_leaf;
  REAL_TYPE fval;
  int nbc=0;
  int ide;
  int tp;
  
  
  // 内部境界条件の有無を調べる
  label_base = "/BcTable/LocalBoundary";
  nbc = tpCntl->countLabels(label_base);
  
  if ( nbc != NoBC)
  {
    stamped_printf("\tLocalBoundary error : '%s'\n", label_base.c_str());
    Exit(0);
  }
  
  //
  //内部境界の条件設定 --- NoBC = 内部境界の数
  //
  
  // BC[@]をサーチ
  for (int odr=1; odr<=NoBC; odr++) {
    
    if( !tpCntl->GetNodeStr(label_base, odr, &str))
    {
      stamped_printf("\tParsing error : No Leaf Node \n");
      Exit(0);
    }
    
    if( strcasecmp(str.substr(0,2).c_str(), "BC") ) continue;
    
    label_leaf = label_base + "/" + str;
    
    // class
    label = label_leaf + "/Class";
    
    if ( !(tpCntl->GetValue(label, &str )) )
    {
      Hostonly_ stamped_printf("\tParsing error : No '%s'\n", label.c_str());
      Exit(0);
    }
    
    // cmp[].type, h_typeのセット ---> setType
    setKeywordLBC(str, odr, cmp);
    
    
    // alias
    label = label_leaf + "/Alias";
    
    if ( !(tpCntl->GetValue(label, &str )) )
    {
      Hostonly_ stamped_printf("\tParsing error : No '%s'\n", label.c_str());
      Exit(0);
    }
    
    // Labelの設定
    cmp[odr].setLabel(str);
    
    
    // Aliasで設定したラベルに対する属性の取得（Polylib.tpの情報を利用する場合）
    if ( C->Mode.Example == id_Polygon )
    {
      string m_pg, m_mat;
      
      for (int i=0; i<C->num_of_polygrp; i++) {
        
        m_pg  = polyP[i].label_grp; // ポリゴンラベル
        m_mat = polyP[i].label_mat; // ポリゴンの媒質ラベル
        printf("mat=%s  grp=%s\n", m_mat.c_str(), m_pg.c_str());
        
        // ポリゴンのラベルとコンポーネントの登録ラベル(alias)が一致する場合
        if ( FBUtility::compare(m_pg, cmp[odr].getLabel()) )
        {
          // ポリゴンの媒質ラベルがMediumTableに含まれているかを調べ，媒質情報を設定
          bool flag = false;
          
          for (int i=1; i<=NoMedium; i++) {
            if ( FBUtility::compare(m_mat, mat[i].getLabel()) )
            {
              cmp[odr].setState(mat[i].getState());
              cmp[odr].setMatOdr(i);
              flag = true;
              break;
            }
          }
          
          // 含まれていない場合はエラー
          if ( !flag )
          {
            Hostonly_ stamped_printf("\tMedium label '%s' associated with Polygon group '%s' is not listed in MediumTable\n", m_mat.c_str(), m_pg.c_str());
            Exit(0);
          }
          
          break;
        }
      }
      
      // aliasで指定したラベルがPolylib.tpの中に見つけられない場合，stateは未設定
      if ( cmp[odr].getState() == -1 )
      {
        Hostonly_ stamped_printf("\tLocal boundary condition '%s' is not listed in Polygon group\n", cmp[odr].getLabel().c_str());
        Exit(0);
      }

    }
    
    
    // 各BCの処理
    tp = cmp[odr].getType();
    
    if ( tp == SPEC_VEL )
    {
      get_IBC_SpecVel(label_leaf, odr, cmp);
    }
    else if ( tp == OUTFLOW ) 
    {
      get_IBC_Outflow(label_leaf, odr, cmp);     
    }
    else if ( tp == IBM_DF ) 
    {
      get_IBC_IBM_DF(label_leaf, odr, cmp);
    }
    else if ( tp == HEX ) 
    {
      get_IBC_PrsLoss(label_leaf, odr, cmp);
    }
    else if ( tp == FAN ) 
    {
      get_IBC_Fan(label_leaf, odr, cmp);
    }
    else if ( tp == DARCY ) 
    {
      get_Darcy(label_leaf, odr, cmp);
    }
    else if ( tp == CELL_MONITOR ) 
    {
      get_IBC_Monitor(label_leaf, odr, cmp);
    }
    else if ( tp == INACTIVE ) 
    {
      ; // skip
    }
    else if ( tp == PERIODIC ) 
    {
      get_IBC_Periodic(label_leaf, odr, cmp);
    }
    else if ( HeatProblem ) // Incase of Heat problem
    {
      if ( C->KindOfSolver == FLOW_ONLY )
      {
        Hostonly_ stamped_printf("Parse Error : Heat BC is not allowed on FLOWONLY mode.\n");
        Exit(0);
      }
      
      if ( Unit_Param != DIMENSIONAL )
      {
        Hostonly_ stamped_printf("\tError : Heat condition must be given by dimensional value\n");
        Exit(0);
      }
      
      if ( tp == ADIABATIC ) 
      {
        //cmp[odr].
        cmp[odr].set_Heatflux( 0.0 );
        
      }
      else if ( tp == HEATFLUX ) 
      {
        get_IBC_HeatFlux(label_leaf, odr, cmp);
       
      }
      else if ( tp == TRANSFER ) 
      {
        switch ( cmp[odr].getHtype() )
        {
          case HT_N:
            get_IBC_HT_N(label_leaf, odr, cmp);

            break;
            
          case HT_S:
            get_IBC_HT_S(label_leaf, odr, cmp, mat);

            break;
            
          case HT_SN:
            get_IBC_HT_SN(label_leaf, odr, cmp, mat);

            break;
            
          case HT_SF:
            get_IBC_HT_SF(label_leaf, odr, cmp);

            break;
            
          case HT_B:
            get_IBC_HT_B(label_leaf, odr, cmp);

            break;
        }        
      }
      else if ( tp == ISOTHERMAL ) 
      {
        get_IBC_IsoTherm(label_leaf, odr, cmp);
      }
      else if ( tp == RADIANT )
      {
        get_IBC_Radiant(label_leaf, odr, cmp);
      }
      else if ( tp == HEAT_SRC ) 
      {
        get_IBC_HeatSrc(label_leaf, odr, cmp);
      }
      else if ( tp == CNST_TEMP ) 
      {
        get_IBC_CnstTemp(label_leaf, odr, cmp);
      }
      else 
      {
        Hostonly_ printf("\tError : Invalid Local BC keyword [%d]\n", tp);
        Exit(0);
      }
    }
  }
  
  // 媒質情報の登録
  for (int i=1; i<=NoMedium; i++) {
    cmp[NoBC+i].setState( mat[i].getState() );
    cmp[NoBC+i].setLabel( mat[i].getLabel() );
    cmp[NoBC+i].setMatOdr(i);
  }
}


// #################################################################
/**
 * @brief パラメータファイルをパースして，外部境界条件を取得，保持する
 * @param [in,out] bc     BoundaryOuter
 * @param [in]     MTITP  MediumTableInfo
 * @param [out]    cmp    CompoList
 */
void ParseBC::loadBC_Outer(BoundaryOuter* bc, const MediumTableInfo *MTITP, CompoList* cmp)
{
  string label_base, label_leaf, label;
  string str;
  
  // Basic Outer BCリストの読み込み
  label_base = "/BcTable/OuterBoundary";
  
  if ( !tpCntl->chkNode(label_base) )
  {
    Hostonly_ stamped_printf("\tParsing error : Missing '%s'\n", label.c_str());
    Exit(0);
  }
  
  
  //Basic_BCs
  for (int i=0; i<NoBaseBC; i++) {
    
    if(!tpCntl->GetNodeStr(label_base, i+1, &str))
    {
      Hostonly_ printf("\tParsing error : Missing 'BasicBCs'\n");
      Exit(0);
    }
    
    if( strcasecmp(str.substr(0,8).c_str(), "BasicBCs") ) continue;
    
    // alias ユニークな名称であること
    label_leaf = label_base + "/" + str;
    label = label_leaf + "/Alias";
    
    if ( !(tpCntl->GetValue(label, &str )) )
    {
      Hostonly_ printf("\tParsing error : No '%s'\n", label.c_str());
      Exit(0);
    }
    if ( !chkDuplicate(i, str) )
    {
      Hostonly_ printf("\tParsing error : 'Alias' must be unique\n");
      Exit(0);
    }
    BaseBc[i].set_Alias(str);
    
    
    
    // Classに境界条件の種別をセットする
    label = label_leaf + "/Class";
    
    if ( !(tpCntl->GetValue(label, &str )) )
    {
      Hostonly_ printf("\tParsing error : No 'Class' in 'BasicBCs'\n");
      Exit(0);
    }
    setKeywordOBC(str, i);
    BaseBc[i].set_Label(str);
    
    
    // 各条件に応じたパラメータをロード
    switch ( BaseBc[i].get_Class() )
    {
      case OBC_WALL:
        get_OBC_Wall(label_leaf, i);
        break;
        
      case OBC_OUTFLOW:
        get_OBC_Outflow(label_leaf, i);
        break;
        
      case OBC_SPEC_VEL:
        get_OBC_SpecVH(label_leaf, i);
        break;
        
      case OBC_TRC_FREE:
        get_OBC_Trcfree(label_leaf, i);
        break;
        
      case OBC_FAR_FIELD:
        get_OBC_FarField(label_leaf, i);
        break;
        
      case OBC_PERIODIC:
        get_OBC_Periodic(label_leaf, i);
        break;
        
      case OBC_SYMMETRIC:
      case OBC_INTRINSIC:
        // nothing to do
        break;
    }
  }
  
  // 各フェイスに境界条件を設定する
  label_base = "/BcTable/OuterBoundary/FaceBC";
  
  if ( !tpCntl->chkNode(label_base) ) 
  {
    Hostonly_ printf("\tParsing error : Missing OuterBoundary '%s'\n", label_base.c_str());
    Exit(0);
  }
  
  // check
  int nnode = tpCntl->countLabels(label_base);
  if ( nnode != NOFACE ) 
  {
    Hostonly_ printf("\tParsing error : OuterBoundary FaceBC count != 6\n");
    Exit(0);
  }
  
  // 各面に与える境界条件番号を取得し，BaseBcから境界情報リストに内容をコピー
  for (int face=0; face<NOFACE; face++) {
    
    // faceに対するラベルを取得
    if ( !tpCntl->GetNodeStr(label_base, face+1, &str) )
    {
      Hostonly_ printf("\tGetNodeStr error\n");
      Exit(0);
    }
    label_leaf = label_base + "/" + str;
    
    // 指定の境界条件を探してBaseBC[]からbc[]へ内容のコピー
    label = label_leaf + "/alias";
    
    if ( !(tpCntl->GetValue(label, &str)) ) 
    {
      Hostonly_ printf("\tParsing error : Alias cannot be found : FaceBC\n");
      Exit(0);
    }
    
    // Aliasでサーチ
    for (int i=0; i<NoBaseBC; i++) {
      if ( !strcasecmp( str.c_str(), BaseBc[i].get_Alias().c_str() ) ) 
      {
        bc[face].dataCopy( &BaseBc[i] );
        break;
      }
      else 
      {
        if ( i == NoBaseBC-1 ) // 最後までみつからない
        {
          Hostonly_ printf("\tParsing error : [%d] '%s' is not listed in 'BasicBCs'\n", i+1, str.c_str());
          Exit(0);
        }
      }
    }
  }
  
  
  // ガイドセルの媒質インデクスをセット
  label_base = "/BcTable/OuterBoundary/FaceBC";
  
  for (int face=0; face<NOFACE; face++) {
    
    if (!tpCntl->GetNodeStr(label_base, face+1, &str))
    {
      Hostonly_ printf("\tGetNodeStr error\n");
      Exit(0);
    }
    label_leaf = label_base + "/" + str;
    
    // ガイドセルの媒質ラベルを取得
    label = label_leaf + "/MediumOnGuideCell";
    
    if ( !(tpCntl->GetValue(label, &str )) ) 
    {
      Hostonly_ printf("\tParsing error : No entory 'MediumOnGuideCell' in 'FaceBC'\n");
      Exit(0);
    }
    
    // ラベル名が一致するエントリ番号をセットする
    for (int i=1; i<=NoMedium; i++) {
      if( !strcasecmp( str.c_str(), MTITP[i].label.c_str() ) )
      {
        bc[face].set_GuideMedium(i);
        break;
      }
    }
    
  }

  
  // 周期境界条件の整合性のチェック
  
  // 部分周期境界の数
  int p_flag=0;
  for (int n=0; n<NOFACE; n++) {
    if (bc[n].get_PrdcMode() == BoundaryOuter::prdc_Driver) p_flag++;
  }
  
  // 部分周期条件を使わない場合，対になる外部境界のチェック
  if ( p_flag == 0 ) {
    int n_pair=0;
    
    // 周期境界条件の指定をチェック
    for (int n=0; n<NOFACE; n++) {
      if ( bc[n].get_Class() == OBC_PERIODIC ) {
        n_pair = oppositeDir(n);
        if ( bc[n_pair].get_Class() != OBC_PERIODIC ) {
          Hostonly_ printf("\tFace BC : No consistent Periodic Bnoudary in %s direction\n", FBUtility::getDirection(n_pair).c_str());
          Exit(0);
        }
      }
    }
    
    // 対になるモードのチェック
    for (int n=0; n<NOFACE; n++) {
      if ( bc[n].get_Class() == OBC_PERIODIC ) {
        n_pair = oppositeDir(n);
        
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
    int n_pair=0;
    for (int n=0; n<NOFACE; n++) {
      n_pair = oppositeDir(n);
      if ( bc[n].get_Class() == OBC_PERIODIC ) {
        if (bc[n].get_PrdcMode() == BoundaryOuter::prdc_Driver) {
          
          // 他方は周期境界以外であること
          if ( bc[n_pair].get_Class() == OBC_PERIODIC ) {
            Hostonly_ printf("\tFace BC : %s direction should be non periodic BC\n", FBUtility::getDirection(n_pair).c_str());
            Exit(0);
          }
          
          int cflag=0;
          for (int c=1; c<=NoBC; c++) {
            if ( cmp[c].getType() == PERIODIC ) {
              if ( (int)cmp[c].getPeriodicDir() != bc[n].get_DriverDir() ) {
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


// #################################################################
/**
 * @brief 外部境界面の反対方向を返す
 * @param [in] dir 評価する方向
 * @return dirと反対方向
 */
int ParseBC::oppositeDir(const int dir)
{
  int n_pair=0;
  
  if      ( dir == X_MINUS ) n_pair=X_PLUS;
  else if ( dir == X_PLUS )  n_pair=X_MINUS;
  else if ( dir == Y_MINUS ) n_pair=Y_PLUS;
  else if ( dir == Y_PLUS )  n_pair=Y_MINUS;
  else if ( dir == Z_MINUS ) n_pair=Z_PLUS;
  else if ( dir == Z_PLUS )  n_pair=Z_MINUS;
  
  return n_pair;
}


// #################################################################
/**
 * @brief コンポーネントの情報を表示する
 * @param [in] fp  ファイルポインタ
 * @param [in] gci グローバルなコンポーネントのインデクス
 * @param [in] mat MediumList
 * @param [in] cmp CompoList
 * @param [in] bc  BoundaryOuter
 */
void ParseBC::printCompo(FILE* fp, const int* gci, const MediumList* mat, CompoList* cmp, const BoundaryOuter* bc)
{
  int n, m;
  bool flag;
  
  // VBC ---------------------------------------------------
  if ( isComponent(SPEC_VEL, cmp) || isComponent(SPEC_VEL_WH, cmp) ) {
    fprintf(fp, "\n\t[SPECIFIED_VELOCITY]\n");
    fprintf(fp, "\t no                    Label   Mat    i_st    i_ed    j_st    j_ed    k_st    k_ed   Area[m*m]\n");
    
    for(n=1; n<=NoBC; n++){
      if ( (cmp[n].getType() == SPEC_VEL) || (cmp[n].getType() == SPEC_VEL_WH) )  {
        fprintf(fp,"\t%3d %24s %5d %7d %7d %7d %7d %7d %7d %11.4e\n",
                n, cmp[n].getLabel().c_str(), cmp[n].getMatOdr(),
                getCmpGbbox_st_x(n, gci), getCmpGbbox_ed_x(n, gci), 
                getCmpGbbox_st_y(n, gci), getCmpGbbox_ed_y(n, gci), 
                getCmpGbbox_st_z(n, gci), getCmpGbbox_ed_z(n, gci), 
                cmp[n].area);
      }
    }
    fprintf(fp, "\n");
    
    fprintf(fp, "\t no                    Label   Mat   normal_x   normal_y   normal_z      Direction");
    
    if ( cmp[n].get_V_Profile() == CompoList::vel_zero ) {
      fprintf(fp, "\n");
    }
    else if ( cmp[n].get_V_Profile() == CompoList::vel_constant ) {
      fprintf(fp, "    Q[m^3/s]   Vel.[m/s]\n");
    }
    else {
      fprintf(fp, "    Q[m^3/s]   Amp.[m/s]   Freq.[Hz]  Phase[rad] Intcpt[m/s] Strauhal[-]\n");
    }
    
    for(n=1; n<=NoBC; n++){
      if ( (cmp[n].getType() == SPEC_VEL) || (cmp[n].getType() == SPEC_VEL_WH) )  {
        fprintf(fp,"\t%3d %24s %5d %10.3e %10.3e %10.3e %14s ",
                n, cmp[n].getLabel().c_str(), cmp[n].getMatOdr(), cmp[n].nv[0], cmp[n].nv[1], cmp[n].nv[2],
                (cmp[n].getBClocation()==CompoList::same_direction) ? "same dir." : "opposite dir.");
        
        if ( cmp[n].get_V_Profile() == CompoList::vel_zero ) {
          fprintf(fp,"\n");
        }
        else if ( cmp[n].get_V_Profile() == CompoList::vel_constant ) {
          fprintf(fp,"%11.3e %11.3e\n", cmp[n].ca[CompoList::bias]*cmp[n].area, cmp[n].ca[CompoList::bias] );
        }
        else {
          fprintf(fp,"%11.3e %11.3e %11.3e %11.3e %11.3e %11.4e\n",
                  cmp[n].ca[CompoList::amplitude]*cmp[n].area,
                  cmp[n].ca[CompoList::amplitude],
                  cmp[n].ca[CompoList::frequency],
                  cmp[n].ca[CompoList::initphase],
                  cmp[n].ca[CompoList::bias],
                  cmp[n].ca[CompoList::frequency]*RefLength/RefVelocity);
        }
        
      }
    }
    
    // with constant temperature
    if ( isComponent(SPEC_VEL_WH, cmp) ) {
      fprintf(fp, "\n\t[SPECIFIED_VELOCITY with constant temperature]\n");
      fprintf(fp, "\t no                    Label   Mat   def    i_st    i_ed    j_st    j_ed    k_st    k_ed      Temp(%s)      Temp[-]\n", (Unit_Temp==Unit_KELVIN) ? "K" : "C");
      
      for(n=1; n<=NoBC; n++) {
        if ( cmp[n].getType() == SPEC_VEL_WH )  {
          fprintf(fp, "\t%3d %24s %5d %5d %7d %7d %7d %7d %7d %7d %12.4e %12.4e\n", 
                  n, cmp[n].getLabel().c_str(), cmp[n].getMatOdr(), cmp[n].getDef(), 
                  getCmpGbbox_st_x(n, gci), getCmpGbbox_ed_x(n, gci), 
                  getCmpGbbox_st_y(n, gci), getCmpGbbox_ed_y(n, gci), 
                  getCmpGbbox_st_z(n, gci), getCmpGbbox_ed_z(n, gci),  
									FBUtility::convK2Temp(cmp[n].get_Temp(), Unit_Temp), 
									FBUtility::convK2ND(cmp[n].get_Temp(), BaseTemp, DiffTemp)); // 保持されている温度はKelvin
        }
      }
      fprintf(fp, "\n");
    }
  }
  
  if ( isComponent(OUTFLOW, cmp) ) {
    fprintf(fp, "\n\t[OUTFLOW]\n");
    fprintf(fp, "\t no                    Label   Mat   def    i_st    i_ed    j_st    j_ed    k_st    k_ed  outflow_vel  pressure\n");
    
    for(n=1; n<=NoBC; n++){
      if ( cmp[n].getType() == OUTFLOW )  {
        fprintf(fp,"\t%3d %24s %5d %5d %7d %7d %7d %7d %7d %7d ",
                n, cmp[n].getLabel().c_str(), cmp[n].getMatOdr(), cmp[n].getDef(),
                getCmpGbbox_st_x(n, gci), getCmpGbbox_ed_x(n, gci), 
                getCmpGbbox_st_y(n, gci), getCmpGbbox_ed_y(n, gci), 
                getCmpGbbox_st_z(n, gci), getCmpGbbox_ed_z(n, gci), 
                (cmp[n].getOutflowType() == V_AVERAGE) ? "Average " : "Minmax ");
        if (cmp[n].get_P_BCtype() == P_DIRICHLET) {
          fprintf(fp, "%12.4e;", FBUtility::convND2D_P(cmp[n].get_Pressure(), BasePrs, RefDensity, RefVelocity, Unit_Prs) );
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
  if ( isComponent(IBM_DF, cmp) ) {
    fprintf(fp, "\n\t[IBM_DF]\n");
    fprintf(fp, "\t no                    Label   Mat    i_st    i_ed    j_st    j_ed    k_st    k_ed\n");
    
    for(n=1; n<=NoBC; n++){
      if ( cmp[n].getType() == IBM_DF )  {
        fprintf(fp,"\t%3d %24s %5d %7d %7d %7d %7d %7d %7d\n",
                n, cmp[n].getLabel().c_str(), cmp[n].getMatOdr(), 
                getCmpGbbox_st_x(n, gci), getCmpGbbox_ed_x(n, gci), 
                getCmpGbbox_st_y(n, gci), getCmpGbbox_ed_y(n, gci), 
                getCmpGbbox_st_z(n, gci), getCmpGbbox_ed_z(n, gci));
      }
    }
    fprintf(fp, "\n");
    
    fprintf(fp, "\t no                    Label   Mat   normal_x   normal_y   normal_z    Vel.[m/s]      Vel.[-]\n");
    
    for(n=1; n<=NoBC; n++){
      if ( cmp[n].getType() == IBM_DF )  {
        fprintf(fp,"\t%3d %24s %5d %10.3e %10.3e %10.3e %12.4e %12.4e \n",
                n, cmp[n].getLabel().c_str(), cmp[n].getMatOdr(), cmp[n].nv[0], cmp[n].nv[1], cmp[n].nv[2],
                cmp[n].get_Velocity(), FBUtility::convD2ND_V(cmp[n].get_Velocity(), RefVelocity));
      }
    }
  }
  
  // Forcing ---------------------------------------------------
	// Heat Exchanger
  if ( isComponent(HEX, cmp) ) {
    fprintf(fp, "\n\t[Heat Exchanger]\n");
    
    fprintf(fp, "\t no                    Label   Mat    normal_x   normal_y   normal_z     O_x[m]     O_y[m]     O_z[m]      dir_x      dir_y      dir_z\n");
    for(n=1; n<=NoBC; n++) {
      if ( cmp[n].getType() == HEX ) {
        fprintf(fp, "\t%3d %24s %5d  %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e\n", 
                n, cmp[n].getLabel().c_str(), cmp[n].getMatOdr(),
								cmp[n].nv[0], cmp[n].nv[1], cmp[n].nv[2],
                cmp[n].oc[0], cmp[n].oc[1], cmp[n].oc[2],
                cmp[n].dr[0], cmp[n].dr[1], cmp[n].dr[2]);
      }
    }
    fprintf(fp, "\n");
    
    fprintf(fp, "\t                                     Depth[m]   Width[m]  Height[m]  Area[m*m]\n");
    for(n=1; n<=NoBC; n++) {
      if ( cmp[n].getType() == HEX ) {
        fprintf(fp, "\t%3d %24s %5d %10.3e %10.3e %10.3e %10.3e\n", 
                n, cmp[n].getLabel().c_str(), cmp[n].getMatOdr(),
								cmp[n].depth, cmp[n].shp_p1, cmp[n].shp_p2, cmp[n].area);
      }
    }
    fprintf(fp, "\n");
    
    fprintf(fp, "\t                                      i_st    i_ed    j_st    j_ed    k_st    k_ed\n");
    for(n=1; n<=NoBC; n++) {
      if ( cmp[n].getType() == HEX ) {
        fprintf(fp, "\t%3d %24s %5d %7d %7d %7d %7d %7d %7d\n", 
                n, cmp[n].getLabel().c_str(), cmp[n].getMatOdr(),
                getCmpGbbox_st_x(n, gci), getCmpGbbox_ed_x(n, gci), 
                getCmpGbbox_st_y(n, gci), getCmpGbbox_ed_y(n, gci), 
                getCmpGbbox_st_z(n, gci), getCmpGbbox_ed_z(n, gci));
      }
    }
    fprintf(fp, "\n");
    
    fprintf(fp, "\t no                    Label   Mat         c1         c2         c3         c4  u_th[m/s]  thick[mm]     vec_forcing\n");
    for(n=1; n<=NoBC; n++) {
      if ( cmp[n].getType() == HEX ) {
        fprintf(fp, "\t%3d %24s %5d %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e     %s\n", 
                n, cmp[n].getLabel().c_str(), cmp[n].getMatOdr(),
                cmp[n].ca[0], cmp[n].ca[1], cmp[n].ca[2], cmp[n].ca[3], cmp[n].ca[4]*RefVelocity, cmp[n].ca[5]*RefLength*1000.0,
                (cmp[n].get_sw_HexDir()==ON) ? "Directional":"Non-directional");
      }
    }
  }
  
  // Fan
  if ( isComponent(FAN, cmp) ) {
    fprintf(fp, "\n\t[Fan]\n");
    
    fprintf(fp, "\t no                    Label   Mat    normal_x   normal_y   normal_z      O_x[m]     O_y[m]     O_z[m]\n");
    for(n=1; n<=NoBC; n++) {
      if ( cmp[n].getType() == FAN ) {
        fprintf(fp, "\t%3d %24s %5d  %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e\n", 
                n, cmp[n].getLabel().c_str(), cmp[n].getMatOdr(),
								cmp[n].nv[0], cmp[n].nv[1], cmp[n].nv[2],
                cmp[n].oc[0], cmp[n].oc[1], cmp[n].oc[2]);
      }
    }
    fprintf(fp, "\n");
    
    fprintf(fp, "\t                                     Depth[m]     Fan[m]    Boss[m]  Area[m*m]\n");
    for(n=1; n<=NoBC; n++) {
      if ( cmp[n].getType() == FAN ) {
        fprintf(fp, "\t%3d %24s %5d %10.3e %10.3e %10.3e %10.3e\n", 
                n, cmp[n].getLabel().c_str(), cmp[n].getMatOdr(),
								cmp[n].depth, cmp[n].shp_p1, cmp[n].shp_p2, cmp[n].area);
      }
    }
    fprintf(fp, "\n");
    
    fprintf(fp, "\t                                      i_st    i_ed    j_st    j_ed    k_st    k_ed\n");
    for(n=1; n<=NoBC; n++) {
      if ( cmp[n].getType() == FAN ) {
        fprintf(fp, "\t%3d %24s %5d %7d %7d %7d %7d %7d %7d\n", 
                n, cmp[n].getLabel().c_str(), cmp[n].getMatOdr(),
                getCmpGbbox_st_x(n, gci), getCmpGbbox_ed_x(n, gci), 
                getCmpGbbox_st_y(n, gci), getCmpGbbox_ed_y(n, gci), 
                getCmpGbbox_st_z(n, gci), getCmpGbbox_ed_z(n, gci));
      }
    }
    fprintf(fp, "\n");
    /*
     fprintf(fp, "\t no                    Label    ID         c1         c2         c3         c4  u_th[m/s]  thick[mm]     vec_forcing\n");
     for(n=1; n<=NoBC; n++) {
     if ( cmp[n].getType() == FAN ) {
     fprintf(fp, "\t%3d %24s %5d %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e     %s\n", 
     n, cmp[n].getLabel().c_str(), cmp[n].getMatOdr(),
     cmp[n].ca[0], cmp[n].ca[1], cmp[n].ca[2], cmp[n].ca[3], cmp[n].ca[4]*RefVelocity, cmp[n].ca[5]*RefLength*1000.0,
     (cmp[n].get_sw_HexDir()==ON) ? "Directional":"Non-directional");
     }
     }*/
  }
  
	// Darcy Law
  if ( isComponent(DARCY, cmp) ) {
    fprintf(fp, "\n\t[Darcy medium]\n");
    
    fprintf(fp, "\t no                    Label   Mat        Area[m*m]   Prmblty_x   Prmblty_y   Prmblty_z[m^2]    Prmblty_x   Prmblty_y   Prmblty_z[-]\n");
    for(n=1; n<=NoBC; n++) {
      if ( cmp[n].getType() == DARCY ) {
        fprintf(fp, "\t%3d %24s %5d %10.3e %11.4e %11.4e %11.4e       %11.4e %11.4e %11.4e\n", 
                n, cmp[n].getLabel().c_str(), cmp[n].getMatOdr(), 
                cmp[n].area, 
								cmp[n].ca[0], cmp[n].ca[1], cmp[n].ca[2], cmp[n].ca[3], cmp[n].ca[4], cmp[n].ca[5]);
      }
    }
    
    fprintf(fp, "\t                                      i_st    i_ed    j_st    j_ed    k_st    k_ed\n");
    for(n=1; n<=NoBC; n++) {
      if ( cmp[n].getType() == DARCY ) {
        fprintf(fp, "\t%3d %24s %5d %7d %7d %7d %7d %7d %7d\n", 
                n, cmp[n].getLabel().c_str(), cmp[n].getMatOdr(),
                getCmpGbbox_st_x(n, gci), getCmpGbbox_ed_x(n, gci), 
                getCmpGbbox_st_y(n, gci), getCmpGbbox_ed_y(n, gci), 
                getCmpGbbox_st_z(n, gci), getCmpGbbox_ed_z(n, gci));
      }
    }
  }
  
  // Heat Face ---------------------------------------------------
  // Adiabatic
  if ( HeatProblem && isComponent(ADIABATIC, cmp) ) {
    fprintf(fp, "\n\t[Adiabatic]\n");
    fprintf(fp, "\t no                    Label   Mat   def\n");
    
    for(n=1; n<=NoBC; n++) {
      if ( cmp[n].getType() == ADIABATIC ) {
        fprintf(fp, "\t%3d %24s %5d %5d\n", n, cmp[n].getLabel().c_str(), cmp[n].getMatOdr(), cmp[n].getDef());
      }
    }
  }
  
  // Direct Heat Flux
  if ( HeatProblem && isComponent(HEATFLUX, cmp) ) {
    fprintf(fp, "\n\t[Direct Heat Flux]\n");
    fprintf(fp, "\t no                    Label   Mat   def    i_st    i_ed    j_st    j_ed    k_st    k_ed   Area[m*m]   flux(W/m^2)        q[-]\n");
    
    for(n=1; n<=NoBC; n++) {
      if ( cmp[n].getType() == HEATFLUX ) {
        fprintf(fp, "\t%3d %24s %5d %5d %7d %7d %7d %7d %7d %7d %11.4e %12.4e %12.4e\n", 
                n, cmp[n].getLabel().c_str(), cmp[n].getMatOdr(), cmp[n].getDef(),
                getCmpGbbox_st_x(n, gci), getCmpGbbox_ed_x(n, gci), 
                getCmpGbbox_st_y(n, gci), getCmpGbbox_ed_y(n, gci), 
                getCmpGbbox_st_z(n, gci), getCmpGbbox_ed_z(n, gci), 
                cmp[n].area, cmp[n].get_Heatflux(), cmp[n].get_Heatflux()/(RefVelocity*DiffTemp*rho*cp));
      }
    }
  }
  
  // Heat Transfer N
  if ( HeatProblem && isCompoTransfer(HT_N, cmp) ) {
    fprintf(fp, "\n\t[Heat Transfer : Type N]\n");
    fprintf(fp, "\t no                    Label   Mat   def    i_st    i_ed    j_st    j_ed    k_st    k_ed   Area[m*m]   coef(W/m^2K)\n");
    
    for(n=1; n<=NoBC; n++) {
      if ( cmp[n].getHtype() == HT_N ) {
        fprintf(fp, "\t%3d %24s %5d %5d %7d %7d %7d %7d %7d %7d %11.4e %12.4e\n", 
                n, cmp[n].getLabel().c_str(), cmp[n].getMatOdr(), cmp[n].getDef(), 
                getCmpGbbox_st_x(n, gci), getCmpGbbox_ed_x(n, gci), 
                getCmpGbbox_st_y(n, gci), getCmpGbbox_ed_y(n, gci), 
                getCmpGbbox_st_z(n, gci), getCmpGbbox_ed_z(n, gci), 
                cmp[n].area, cmp[n].get_CoefHT());
      }
    }
  }
  
  // Heat Transfer S
  if ( HeatProblem && isCompoTransfer(HT_S, cmp) ) {
    fprintf(fp, "\n\t[Heat Transfer : Type S]\n");
    fprintf(fp, "\t no                    Label   Mat   def    i_st    i_ed    j_st    j_ed    k_st    k_ed   Area[m*m]   coef(W/m^2K)   Temp(%s)   Temp[-]\n", (Unit_Temp==Unit_KELVIN) ? "K" : "C");
    
    for(n=1; n<=NoBC; n++) {
      if ( cmp[n].getHtype() == HT_S ) {
        fprintf(fp, "\t%3d %24s %5d %5d %7d %7d %7d %7d %7d %7d %11.4e %12.4e %12.4e %12.4e\n", 
                n, cmp[n].getLabel().c_str(), cmp[n].getMatOdr(), cmp[n].getDef(), 
                getCmpGbbox_st_x(n, gci), getCmpGbbox_ed_x(n, gci), 
                getCmpGbbox_st_y(n, gci), getCmpGbbox_ed_y(n, gci), 
                getCmpGbbox_st_z(n, gci), getCmpGbbox_ed_z(n, gci), 
                cmp[n].area, cmp[n].get_CoefHT(), FBUtility::convK2Temp(cmp[n].get_Temp(), Unit_Temp),
                FBUtility::convK2ND(cmp[n].get_Temp(), BaseTemp, DiffTemp) );
      }
    }
  }
  
  // Heat Transfer SN
  if ( HeatProblem && isCompoTransfer(HT_SN, cmp) ) {
    fprintf(fp, "\n\t[Heat Transfer : Type SN]\n"); 
    fprintf(fp, "\t no                    Label   Mat   def    i_st    i_ed    j_st    j_ed    k_st    k_ed   Area[m*m]      Temp(%s)      Temp[-]   Type\n", (Unit_Temp==Unit_KELVIN) ? "K" : "C");
    
    for(n=1; n<=NoBC; n++) {
      if ( cmp[n].getHtype() == HT_SN ) {
        fprintf(fp, "\t%3d %24s %5d %5d %7d %7d %7d %7d %7d %7d %11.4e %12.4e %12.4e   %s\n", 
                n, cmp[n].getLabel().c_str(), cmp[n].getMatOdr(), cmp[n].getDef(), 
                getCmpGbbox_st_x(n, gci), getCmpGbbox_ed_x(n, gci), 
                getCmpGbbox_st_y(n, gci), getCmpGbbox_ed_y(n, gci), 
                getCmpGbbox_st_z(n, gci), getCmpGbbox_ed_z(n, gci), 
                cmp[n].area, FBUtility::convK2Temp(cmp[n].get_Temp(), Unit_Temp),
                FBUtility::convK2ND(cmp[n].get_Temp(), BaseTemp, DiffTemp),
                (cmp[n].get_sw_HTmodeRef()==CompoList::HT_mode_bulk) ? "Bulk" : "Local");
      }
    }
    fprintf(fp, "\n");
    
    fprintf(fp, "\t no                    Label   Mat  vert_lam_a  vert_lam_b vert_turb_a vert_turb_b  vert_Ra_cr   lwr_lam_a   lwr_lam_b  lwr_turb_a  lwr_turb_b   lwr_Ra_cr\n");
    for(n=1; n<=NoBC; n++) {
      if ( cmp[n].getHtype() == HT_SN ) {
        fprintf(fp, "\t%3d %24s %5d %11.4e %11.4e %11.4e %11.4e %11.4e %11.4e %11.4e %11.4e %11.4e %11.4e\n", 
                n, 
                cmp[n].getLabel().c_str(),
                cmp[n].getMatOdr(), 
                cmp[n].ca[CompoList::vert_laminar_alpha], 
                cmp[n].ca[CompoList::vert_laminar_beta], 
                cmp[n].ca[CompoList::vert_turbulent_alpha], 
                cmp[n].ca[CompoList::vert_turbulent_beta], 
                cmp[n].ca[CompoList::vert_Ra_critial], 
                cmp[n].cb[CompoList::lower_laminar_alpha], 
                cmp[n].cb[CompoList::lower_laminar_beta], 
                cmp[n].cb[CompoList::lower_turbulent_alpha], 
                cmp[n].cb[CompoList::lower_turbulent_beta], 
                cmp[n].cb[CompoList::lower_Ra_critial] );
      }
    }
  }
  
  // Heat Transfer SF
  if ( HeatProblem && isCompoTransfer(HT_SF, cmp) ) {
    fprintf(fp, "\n\t[Heat Transfer : Type SF]\n");
    fprintf(fp, "\t no                    Label   Mat   def    i_st    i_ed    j_st    j_ed    k_st    k_ed   Area[m*m]      Temp(%s)      Temp[-]   Type\n", 
            (Unit_Temp==Unit_KELVIN) ? "K" : "C");
    
    for(n=1; n<=NoBC; n++) {
      if ( cmp[n].getHtype() == HT_SF ) {
        fprintf(fp, "\t%3d %24s %5d %5d %7d %7d %7d %7d %7d %7d %11.4e %12.4e %12.4e\n", 
                n, 
                cmp[n].getLabel().c_str(), 
                cmp[n].getMatOdr(), 
                cmp[n].getDef(), 
                getCmpGbbox_st_x(n, gci), getCmpGbbox_ed_x(n, gci), 
                getCmpGbbox_st_y(n, gci), getCmpGbbox_ed_y(n, gci), 
                getCmpGbbox_st_z(n, gci), getCmpGbbox_ed_z(n, gci), 
                cmp[n].area, FBUtility::convK2Temp(cmp[n].get_Temp(), Unit_Temp),
                FBUtility::convK2ND(cmp[n].get_Temp(), BaseTemp, DiffTemp),
                (cmp[n].get_sw_HTmodeRef()==CompoList::HT_mode_bulk) ? "Bulk" : "Local");
      }
    }
    fprintf(fp, "\t no                    Label   Mat       alpha        beta       gamma\n");
    for(n=1; n<=NoBC; n++) {
      if ( cmp[n].getHtype() == HT_SN ) {
        fprintf(fp, "\t%3d %24s %5d %11.4e %11.4e %11.4e\n", 
                n, 
                cmp[n].getLabel().c_str(), 
                cmp[n].getMatOdr(), 
                cmp[n].ca[CompoList::alpha], 
                cmp[n].ca[CompoList::beta], 
                cmp[n].ca[CompoList::gamma]);
      }
    }
  }
  
  // Heat Transfer B
  if ( HeatProblem && isCompoTransfer(HT_B, cmp)) {
    fprintf(fp, "\n\t[Heat Transfer : Type B]\n");
    fprintf(fp, "\t no                    Label   Mat   def    i_st    i_ed    j_st    j_ed    k_st    k_ed   Area[m*m]  coef(W/m^2K)  BulkTemp(%s)   BulkTemp[-]\n", 
            (Unit_Temp==Unit_KELVIN) ? "K" : "C");
    
    for(n=1; n<=NoBC; n++) {
      if ( cmp[n].getHtype() == HT_B ) {
        fprintf(fp, "\t%3d %24s %5d %5d %7d %7d %7d %7d %7d %7d %11.4e  %12.4e %12.4e %12.4e\n", 
                n, cmp[n].getLabel().c_str(), cmp[n].getMatOdr(), cmp[n].getDef(), 
                getCmpGbbox_st_x(n, gci), getCmpGbbox_ed_x(n, gci), 
                getCmpGbbox_st_y(n, gci), getCmpGbbox_ed_y(n, gci), 
                getCmpGbbox_st_z(n, gci), getCmpGbbox_ed_z(n, gci), 
                cmp[n].area, cmp[n].get_CoefHT(), FBUtility::convK2Temp(cmp[n].get_Temp(), Unit_Temp),
                FBUtility::convK2ND(cmp[n].get_Temp(), BaseTemp, DiffTemp) );
      }
    }
  }
  
  // Iso Thermal
  if ( HeatProblem && isComponent(ISOTHERMAL, cmp)) {
    fprintf(fp, "\n\t[Iso-Thermal]\n");
    fprintf(fp, "\t no                    Label   Mat   def    i_st    i_ed    j_st    j_ed    k_st    k_ed   Area[m*m]   Sf.Temp(%s)   Sf.Temp[-]\n", 
            (Unit_Temp==Unit_KELVIN) ? "K" : "C");
    
    for(n=1; n<=NoBC; n++) {
      if ( cmp[n].getType() == ISOTHERMAL ) {
        fprintf(fp, "\t%3d %24s %5d %5d %7d %7d %7d %7d %7d %7d %11.4e %12.4e %12.4e \n", 
                n, cmp[n].getLabel().c_str(), cmp[n].getMatOdr(), cmp[n].getDef(), 
                getCmpGbbox_st_x(n, gci), getCmpGbbox_ed_x(n, gci), 
                getCmpGbbox_st_y(n, gci), getCmpGbbox_ed_y(n, gci), 
                getCmpGbbox_st_z(n, gci), getCmpGbbox_ed_z(n, gci), 
								cmp[n].area, FBUtility::convK2Temp(cmp[n].get_Temp(), Unit_Temp),
                FBUtility::convK2ND(cmp[n].get_Temp(), BaseTemp, DiffTemp) );
      }
    }
  }
  
  // Radiant
  if ( HeatProblem && isComponent(RADIANT, cmp)) {
    fprintf(fp, "\n\t[Radiant]\n");
    fprintf(fp, "\t no                    Label   Mat    def    i_st    i_ed    j_st    j_ed    k_st    k_ed   Area[m*m]   ep[-]   pj[-]\n");
    
    for(n=1; n<=NoBC; n++) {
      if ( cmp[n].getType() == RADIANT ) {
        fprintf(fp, "\t%3d %24s %5d %5d %7d %7d %7d %7d %7d %7d %11.4e %12.4e %12.4e\n", 
                n, cmp[n].getLabel().c_str(), cmp[n].getMatOdr(), cmp[n].getDef(), 
                getCmpGbbox_st_x(n, gci), getCmpGbbox_ed_x(n, gci), 
                getCmpGbbox_st_y(n, gci), getCmpGbbox_ed_y(n, gci), 
                getCmpGbbox_st_z(n, gci), getCmpGbbox_ed_z(n, gci), 
                cmp[n].area, cmp[n].get_CoefRadEps(), cmp[n].get_CoefRadPrj());
      }
    }
  }
  
  // Heat source ---------------------------------------------------
  // Heat generation
  if ( HeatProblem && isComponent(HEAT_SRC, cmp)) {
    fprintf(fp, "\n\t[Heat Generation]\n");
    fprintf(fp, "\t no                    Label   Mat    i_st    i_ed    j_st    j_ed    k_st    k_ed     Q[W/m^3]    nrmlzd[-]\n");
    
    for(n=1; n<=NoBC; n++) {
      int h_odr = cmp[n].getMatOdr();
      if ( cmp[n].getType() == HEAT_SRC ) {
        fprintf(fp, "\t%3d %24s %5d %7d %7d %7d %7d %7d %7d %12.4e %12.4e\n", 
                n, cmp[n].getLabel().c_str(), cmp[n].getMatOdr(), 
                getCmpGbbox_st_x(n, gci), getCmpGbbox_ed_x(n, gci), 
                getCmpGbbox_st_y(n, gci), getCmpGbbox_ed_y(n, gci), 
                getCmpGbbox_st_z(n, gci), getCmpGbbox_ed_z(n, gci), 
                cmp[n].get_HeatValue(), 
                FBUtility::convD2ND_Hsrc(cmp[n].get_HeatValue(), RefVelocity, RefLength, DiffTemp, mat[h_odr].P[p_density], mat[h_odr].P[p_specific_heat]));
      }
    }
  }
  
  // Constant Temperature
  if ( HeatProblem && isComponent(CNST_TEMP, cmp)) {
    fprintf(fp, "\n\t[Constant Temperature]\n");
    fprintf(fp, "\t no                    Label   Mat    i_st    i_ed    j_st    j_ed    k_st    k_ed      Temp[%s]      Temp[-]\n", 
            (Unit_Temp==Unit_KELVIN) ? "K" : "C");
    
    for(n=1; n<=NoBC; n++) {
      if ( cmp[n].getType() == CNST_TEMP ) {
        fprintf(fp, "\t%3d %24s %5d %7d %7d %7d %7d %7d %7d %12.4e %12.4e\n", 
                n, cmp[n].getLabel().c_str(), cmp[n].getMatOdr(), 
                getCmpGbbox_st_x(n, gci), getCmpGbbox_ed_x(n, gci), 
                getCmpGbbox_st_y(n, gci), getCmpGbbox_ed_y(n, gci), 
                getCmpGbbox_st_z(n, gci), getCmpGbbox_ed_z(n, gci), 
                FBUtility::convK2Temp(cmp[n].get_Temp(), Unit_Temp), 
                FBUtility::convK2ND(cmp[n].get_Temp(), BaseTemp, DiffTemp));
      }
    }
  }
  
  // Monitor ---------------------------------------------------
  if ( isComponent(CELL_MONITOR, cmp) ) {
    fprintf(fp, "\n\t[Monitor]\n");
    fprintf(fp, "\t no                    Label   Mat    i_st    i_ed    j_st    j_ed    k_st    k_ed   Area[m*m]  Variables\n");
    
    for(n=1; n<=NoBC; n++) {
      if ( cmp[n].getType() == CELL_MONITOR ) 
      {
        fprintf(fp, "\t%3d %24s %5d %7d %7d %7d %7d %7d %7d %11.4e  %s\n", 
                n, cmp[n].getLabel().c_str(), cmp[n].getMatOdr(), 
                getCmpGbbox_st_x(n, gci), getCmpGbbox_ed_x(n, gci), 
                getCmpGbbox_st_y(n, gci), getCmpGbbox_ed_y(n, gci), 
                getCmpGbbox_st_z(n, gci), getCmpGbbox_ed_z(n, gci), 
								cmp[n].area, cmp[n].getVarStr() );
      }
    }
    fprintf(fp, "\n");
    
    fprintf(fp, "\t no                    Label   Mat   normal_x   normal_y   normal_z Reference\n");
    for(n=1; n<=NoBC; n++){
      if ( cmp[n].getType() == CELL_MONITOR )  
      {
        fprintf(fp,"\t%3d %24s %5d %10.3e %10.3e %10.3e       %3s\n",
                n, cmp[n].getLabel().c_str(), cmp[n].getMatOdr(), cmp[n].nv[0], cmp[n].nv[1], cmp[n].nv[2],
                (cmp[n].getStateCellMonitor()==ON) ? "yes" : "no");
      }
    }
  }
  
  // Periodic ---------------------------------------------------
  if ( isComponent(PERIODIC, cmp) ) {
    fprintf(fp, "\n\t[Periodic]\n");
    fprintf(fp, "\t no                    Label   Mat    i_st    i_ed    j_st    j_ed    k_st    k_ed    Pressure Difference [Pa]/[-]  Driver\n");
    
    int dir_in=0, dir_out=0, pp_in=0, pp_out=0;
    
    for(n=1; n<=NoBC; n++) {
      if ( cmp[n].getType() == PERIODIC ) {
        dir_in = cmp[n].getPeriodicDir();
        
        fprintf(fp, "\t%3d %24s %5d %7d %7d %7d %7d %7d %7d     ", 
                n, cmp[n].getLabel().c_str(), cmp[n].getMatOdr(), 
                getCmpGbbox_st_x(n, gci), getCmpGbbox_ed_x(n, gci), 
                getCmpGbbox_st_y(n, gci), getCmpGbbox_ed_y(n, gci), 
                getCmpGbbox_st_z(n, gci), getCmpGbbox_ed_z(n, gci));
        fprintf(fp,"%12.6e / %12.6e ", cmp[n].ca[0], FBUtility::convD2ND_P(cmp[n].ca[0], BasePrs, RefDensity, RefVelocity, Unit_Prs));
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



// #################################################################
/**
 * @brief 外部境界条件の各面の情報を表示する
 * @param [in] fp    ファイルポインタ
 * @param [in] G_reg グローバルの領域の大きさ
 * @param [in] bc    BoundaryOuter
 * @param [in] mat   MediumList
 */
void ParseBC::printFaceOBC(FILE* fp, const REAL_TYPE* G_reg, const BoundaryOuter* bc, const MediumList* mat)
{
  for (int i=0; i<NOFACE; i++) {
    fprintf(fp,"\t      Set %s up as %s : < %s >\n", 
            FBUtility::getDirection(i).c_str(), 
            getOBCstr(bc[i].get_Class()).c_str(), 
            bc[i].get_Alias().c_str());
    printOBC(fp, &bc[i], mat, G_reg, i);
    fprintf(fp,"\n");
  }
  fflush(fp);
}



// #################################################################
/**
 * @brief 速度の外部境界条件処理の表示
 * @param [in] fp    ファイルポインタ
 * @param [in] ref   BoundaryOuter
 * @param [in] mat   MediumList
 * @param [in] G_reg グローバルの領域の大きさ
 * @param [in] face  面番号
 */
void ParseBC::printOBC(FILE* fp, const BoundaryOuter* ref, const MediumList* mat, const REAL_TYPE* G_reg, const int face)
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
        int htp = ref->get_hType();
        if ( htp == ADIABATIC ) {
          fprintf(fp,"\t\t\tAdiabatic\n");
        }
        else if ( htp == TRANSFER) {
          int ht_mode = ref->get_HTmode();
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
      fprintf(fp,"\t\t\t%12.6e [Pa]  /  %12.6e [-]\n", ref->p, FBUtility::convD2ND_P(ref->p, BasePrs, RefDensity, RefVelocity, Unit_Prs));
      
      if ( HeatProblem ) {
        fprintf(fp, "\t\t\t    Ambient Temperature  = %12.6e \n", ref->get_Temp());
      }
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
                      ref->p/(G_reg[0]*RefLength), FBUtility::convD2ND_P(ref->p, BasePrs, RefDensity, RefVelocity, Unit_Prs)/G_reg[0]);
              break;
              
            case Y_MINUS:
            case Y_PLUS:
              fprintf(fp,"\t\t\tPressure Gradient   = %12.6e [Pa/m]  /  %12.6e [-]\n", 
                      ref->p/(G_reg[1]*RefLength), FBUtility::convD2ND_P(ref->p, BasePrs, RefDensity, RefVelocity, Unit_Prs)/G_reg[1]);
              break;
              
            case Z_MINUS:
            case Z_PLUS:
              fprintf(fp,"\t\t\tPressure Gradient   = %12.6e [Pa/m]  /  %12.6e [-]\n", 
                      ref->p/(G_reg[2]*RefLength), FBUtility::convD2ND_P(ref->p, BasePrs, RefDensity, RefVelocity, Unit_Prs)/G_reg[2]);
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
    case OBC_INTRINSIC:
    case OBC_FAR_FIELD:
      break;
      
      
    default:
      printf("\n\tError : OuterBC\n");
      Exit(0);
      break;
  }
  
  fflush(fp);
}



// #################################################################
//@brief 変数の初期化
void ParseBC::setControlVars(Control* Cref)
{
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
  double s, two=2.0;
  
  s = (double)MASK_6; // bit幅マスクは2^(bit幅)-1を表し，ちょうど0を除いた個数となっている
  m = (int)(log10(s+1.0)/log10(two) );
  
  if ( NoCompo > s )
  {
    printf("Error : No. of Component (= NoBC + NoMedium) must be less or equal %d(%dbit-width)\n", (int)s, m);
    Exit(0);
  }
  
  s = (double)(MASK_5-1); // 0と31を除く
  m = (int)( log10(s+2.0)/log10(two) );
  
  if ( NoBC > s )
  {
    printf("Error : No. of BC must be less or equal %d(%dbit-width)\n", (int)s, m);
    Exit(0);
  }
  
  // 外部境界条件の種類数を取得し，内部保持配列をインスタンス
  string label;
  string str;
  
  label = "/BcTable/OuterBoundary";
  
  if ( !tpCntl->chkNode(label) )
  {
    Hostonly_ stamped_printf("\tParsing error : Missing '%s\n", label.c_str());
    Exit(0);
  }
  
  // タグ内の数をチェック
  int nnode = tpCntl->countLabels(label);
  if ( nnode == 0 )
  {
    Hostonly_ stamped_printf("\tNo OuterBoundary defined\n");
    return;
  }
  
  int counter=0;
  for (int i=1; i<=nnode; i++) {
    
    if ( !tpCntl->GetNodeStr(label, i, &str) )
    {
      Hostonly_ stamped_printf("\tGetNodeStr error\n");
      Exit(0);
    }
    
    if( !strcasecmp(str.substr(0,8).c_str(), "BasicBCs") ) counter++;
  }
  
  NoBaseBC = counter;
  BaseBc = new BoundaryOuter[NoBaseBC];
}



// #################################################################
/**
 * @brief 内部境界条件の照合を行う
 * @param [in]  keyword テストキーワード
 * @param [in]  m       BaseBcの格納番号
 * @param [out] cmp     CompoList
 * @note SPEC_VEL_WHは陽には現れず，get_IBC_SpecVel()内で登録される
 */
void ParseBC::setKeywordLBC(const string keyword, const int m, CompoList* cmp)
{
  if     ( FBUtility::compare(keyword, "Adiabatic") )            cmp[m].setType(ADIABATIC);
  else if( FBUtility::compare(keyword, "DirectHeatFlux") )       cmp[m].setType(HEATFLUX);
  else if( FBUtility::compare(keyword, "HeatTransferB") ) {      cmp[m].setType(TRANSFER); cmp[m].setHtype(HT_B); }
  else if( FBUtility::compare(keyword, "HeatTransferN") ) {      cmp[m].setType(TRANSFER); cmp[m].setHtype(HT_N); }
  else if( FBUtility::compare(keyword, "HeatTransferS") ) {      cmp[m].setType(TRANSFER); cmp[m].setHtype(HT_S); }
  else if( FBUtility::compare(keyword, "HeatTransferSF") ){      cmp[m].setType(TRANSFER); cmp[m].setHtype(HT_SF); }
  else if( FBUtility::compare(keyword, "HeatTransferSN") ){      cmp[m].setType(TRANSFER); cmp[m].setHtype(HT_SN); }
  else if( FBUtility::compare(keyword, "IsoThermal") )           cmp[m].setType(ISOTHERMAL);
  else if( FBUtility::compare(keyword, "Radiation") )            cmp[m].setType(RADIANT);
  else if( FBUtility::compare(keyword, "SpecifiedVelocity") )    cmp[m].setType(SPEC_VEL);
  else if( FBUtility::compare(keyword, "Outflow") )              cmp[m].setType(OUTFLOW);
  else if( FBUtility::compare(keyword, "Forcing") )              cmp[m].setType(IBM_DF);
  else if( FBUtility::compare(keyword, "HeatSource") )           cmp[m].setType(HEAT_SRC);
  else if( FBUtility::compare(keyword, "SpecifiedTemperature") ) cmp[m].setType(CNST_TEMP);
  else if( FBUtility::compare(keyword, "PressureLoss") )         cmp[m].setType(HEX);
  else if( FBUtility::compare(keyword, "Fan") )                  cmp[m].setType(FAN);
  else if( FBUtility::compare(keyword, "Darcy") )                cmp[m].setType(DARCY);
  else if( FBUtility::compare(keyword, "CellMonitor") )          cmp[m].setType(CELL_MONITOR);
  else if( FBUtility::compare(keyword, "inactive") )             cmp[m].setType(INACTIVE);
  else if( FBUtility::compare(keyword, "Periodic") )             cmp[m].setType(PERIODIC);
  else {
    Hostonly_ stamped_printf("\tInvalid keyword is described '%s'\n", keyword.c_str());
    Exit(0);
  }
}


// #################################################################
/**
 * @brief 内部境界条件の照合を行う
 * @param [in] keyword テストキーワード
 * @param [in] m       コンポーネントの格納番号
 */
void ParseBC::setKeywordOBC(const string keyword, const int m)
{
  if     ( FBUtility::compare(keyword, "Wall") )              BaseBc[m].set_Class(OBC_WALL);
  else if( FBUtility::compare(keyword, "Outflow") )           BaseBc[m].set_Class(OBC_OUTFLOW);
  else if( FBUtility::compare(keyword, "SpecifiedVelocity") ) BaseBc[m].set_Class(OBC_SPEC_VEL);
  else if( FBUtility::compare(keyword, "Symmetric") )         BaseBc[m].set_Class(OBC_SYMMETRIC);
  else if( FBUtility::compare(keyword, "Periodic") )          BaseBc[m].set_Class(OBC_PERIODIC);
  else if( FBUtility::compare(keyword, "TractionFree") )      BaseBc[m].set_Class(OBC_TRC_FREE);
  else if( FBUtility::compare(keyword, "FarField") )          BaseBc[m].set_Class(OBC_FAR_FIELD);
  else if( FBUtility::compare(keyword, "intrinsic") )         BaseBc[m].set_Class(OBC_INTRINSIC);
  else {
    Hostonly_ stamped_printf("\tParsing error : Invalid keyword is described '%s'\n", keyword.c_str());
    Exit(0);
  }
}


// #################################################################
/**
 * @brief 指定した媒質IDから参照物理量を設定する
 * @param [in] mat MediumList
 * @param [in] cmp CompoList
 * @param [in] Ref 参照媒質番号
 */
void ParseBC::setRefMediumProperty(const MediumList* mat, const CompoList* cmp, const int Ref)
{
  int m;
  
  for (int n=NoBC+1; n<=NoCompo; n++) {
    m = cmp[n].getMatOdr();
    
    if ( m == Ref ) 
    {
      if ( mat[m].getState() == FLUID ) 
      {
        rho    = mat[m].P[p_density];
        nyu    = mat[m].P[p_kinematic_viscosity];
        cp     = mat[m].P[p_specific_heat];
        lambda = mat[m].P[p_thermal_conductivity];
        beta   = mat[m].P[p_vol_expansion]; // can be replaced by 1/K in the case of gas
        //mu    = mat[m].P[p_viscosity];
        //snd_spd = mat[m].P[p_sound_of_speed];
      }
      else 
      {
        rho    = mat[m].P[p_density];
        cp     = mat[m].P[p_specific_heat];
        lambda = mat[m].P[p_thermal_conductivity];
      }
    }
  }
  
  RefDensity      = rho;
  RefSpecificHeat = cp;
}
