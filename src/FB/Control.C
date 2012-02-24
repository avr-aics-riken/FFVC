/*
 * SPHERE - Skeleton for PHysical and Engineering REsearch
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2004-2012
 *
 */

//@file Control.C
//@brief FlowBase Control class
//@author keno, FSI Team, VCAD, RIKEN

#include "Control.h"

/**
 @fn void Control::getXML_Mon_Line(MonitorList* M, const CfgElem* elmL2, REAL_TYPE from[3], REAL_TYPE to[3], int& nDivision)
 @brief XMLに記述されたモニタ座標情報(Line)を取得
 @param M MonitorList クラスオブジェクトのポインタ
 @param elmL2 コンフィギュレーションツリーのポインタ
 @param from Line始点座標
 @param to   Line終点座標
 @param nDivision Line分割数
 @note データは無次元化して保持
 */
void Control::getXML_Mon_Line(MonitorList* M, const CfgElem* elmL2, REAL_TYPE from[3], REAL_TYPE to[3], int& nDivision)
{
  const CfgElem *elmL3=NULL;
  
  if ( !elmL2->GetValue("division", &nDivision) ) Exit(0);
  if ( nDivision == 0 ) Exit(0);
  
  // load parameter of 'from' and 'to'
  if ( !elmL2->GetValue("from", &elmL3) ) {
    Hostonly_ stamped_printf("\tParsing error : fail to get 'from' in 'line' >> %s\n", elmL3->GetName());
    Exit(0);
  }
  else {
    if ( !elmL3->GetVctValue("x", "y", "z", &from[0], &from[1], &from[2]) ) {
      Hostonly_ stamped_printf("\tParsing error : fail to get vec params in 'from'\n");
      Exit(0);
    }
    if (Sampling.unit == DIMENSIONAL) {
      normalizeCord(from);
    }
  }
  if ( !elmL2->GetValue("to", &elmL3) ) {
    Hostonly_ stamped_printf("\tParsing error : fail to get 'to' in 'line' >> %s\n", elmL3->GetName());
    Exit(0);
  }
  else {
    if ( !elmL3->GetVctValue("x", "y", "z", &to[0], &to[1], &to[2]) ) {
      Hostonly_ stamped_printf("\tParsing error : fail to get vec params in 'to'\n");
      Exit(0);
    }
    if (Sampling.unit == DIMENSIONAL) {
      normalizeCord(to);
    }
  }
}

/**
 @fn void Control::getXML_Mon_Pointset(MonitorList* M, const CfgElem *elmL2, vector<MonitorCompo::MonitorPoint>& pointSet)
 @brief XMLに記述されたモニタ座標情報を取得(PointSet)
 @param M MonitorList クラスオブジェクトのポインタ
 @param m order
 @param elmL2 コンフィギュレーションツリーのポインタ
 @param pointSet PointSet配列
 @note データは無次元化して保持
 */
void Control::getXML_Mon_Pointset(MonitorList* M, const CfgElem *elmL2, vector<MonitorCompo::MonitorPoint>& pointSet)
{
  const CfgElem *elmL3=NULL;
  REAL_TYPE v[3];
  const char* str=NULL;
  char tmpstr[20];
  
  v[0] = v[1] = v[2] = 0.0;
  
  // load parameter for a set
  elmL3 = elmL2->GetElemFirst();

  for (unsigned j=0; j<elmL2->GetElemSize(); j++) {
    
    if ( strcasecmp("set", elmL3->GetName()) ) { // not agree
      Hostonly_ stamped_printf("\tParsing error : fail to get 'set' in 'point_set' >> %s\n", elmL3->GetName());
      Exit(0);
    }
    if ( !elmL3->GetVctValue("x", "y", "z", &v[0], &v[1], &v[2]) ) {
      Hostonly_ stamped_printf("\tParsing error : fail to get vec params in 'point_set'\n");
      Exit(0);
    }
    if (Sampling.unit == DIMENSIONAL) {
      normalizeCord(v);
    }
    
    // set Labelの取得．ラベルなしでもエラーではない
    if ( !(str = elmL3->GetComment()) ) {
      Hostonly_ stamped_printf("\tParsing warning : No commnet for 'set'\n");
    }
    if ( !str ) {
      sprintf(tmpstr, "point_%d",j);
      str = tmpstr;
    }
    
    pointSet.push_back(MonitorCompo::MonitorPoint(v, str));
    
    elmL3 = elmL2->GetElemNext(elmL3); // ahead on the next pointer
  }  
}

/**
 @fn void Control::getXML_Monitor(MonitorList* M)
 @brief XMLに記述されたモニタ座標情報を取得し，リストに保持する
 @param M MonitorList クラスオブジェクトのポインタ
 */
void Control::getXML_Monitor(MonitorList* M)
{
  const CfgElem *elemTop=NULL, *elmL1=NULL, *elmL2=NULL;
  const char* p=NULL;
  const char* str=NULL;
  const char* label=NULL;
  const char* method=NULL;
  const char* mode=NULL;
  const CfgParam * param=NULL;
  REAL_TYPE f_val=0.0;
  //int nvar=0;
  MonitorCompo::Type type;
  vector<string> variables;
  
  // Monitor_List section is already confirmed
  elemTop = CF->GetTop(STEER);
  elmL1 = elemTop->GetElemFirst("Monitor_List");
  
  // ログ出力
  if ( !elmL1->GetValue(CfgIdt("log"), &str) ) {
		Hostonly_ stamped_printf("\tParsing error : Invalid string for 'Log' in 'Monitor_List'\n");
		Exit(0);
	}
  if     ( !strcasecmp(str, "on") )   Sampling.log = ON;
  else if( !strcasecmp(str, "off") )  Sampling.log = OFF;
  else {
    Hostonly_ stamped_printf("\tInvalid keyword is described for 'Monitor_List'\n");
    Exit(0);
  }
  
  if ( Sampling.log == OFF ) return;
  
  // 出力ファイル名
  if ( !elmL1->GetValue(CfgIdt("output_file"), &str) ) {
    Hostonly_ stamped_printf("\tParsing error : fail to get 'Output_File' in 'Monitor_List'\n");
    Exit(0);
  }
  strcpy(HistoryMonitorName, str);
  
  // 集約モード
  if ( !elmL1->GetValue(CfgIdt("output_mode"), &str) ) {
    Hostonly_ stamped_printf("\tParsing error : fail to get 'Output_Mode' in 'Monitor_List'\n");
    Exit(0);
  }
  if ( !strcasecmp("gather", str)) {
    Sampling.out_mode = MonitorList::GATHER;
    M->setOutputType(MonitorList::GATHER);
  }
  else if ( !strcasecmp("distribute", str)) {
    Sampling.out_mode = MonitorList::DISTRIBUTE;
    M->setOutputType(MonitorList::DISTRIBUTE);
  }
  else {
    Hostonly_ stamped_printf("\tParsing error : Invalid keyord for 'Output_Mode'\n");
    Exit(0);
  }
  
  // サンプリング間隔
  if ( !elmL1->GetValue("Sampling_Interval_Type", &str) ) {
    Hostonly_ stamped_printf("\tParsing error : fail to get 'Sampling_Interval_Type' in 'Monitor_List'\n");
    Exit(0);
  }
  else {
    if ( !strcasecmp(str, "step") ) {
      Interval[Interval_Manager::tg_sampled].setMode_Step();
    }
    else if ( !strcasecmp(str, "time") ) {
      Interval[Interval_Manager::tg_sampled].setMode_Time();
    }
    else {
      Hostonly_ stamped_printf("\tParsing error : Invalid keyword for 'Sampling_Interval_Type' in 'Monitor_List'\n");
      Exit(0);
    }
    
    if ( elmL1->GetValue("Sampling_Interval", &f_val) ) {
      Interval[Interval_Manager::tg_sampled].setInterval((double)f_val);
    }
    else {
      Hostonly_ stamped_printf("\tParsing error : fail to get 'Sampling_Interval' in 'Monitor_List'\n");
      Exit(0);
    }
  }
  
  // 単位指定
  if ( !elmL1->GetValue("Unit", &str) ) {
		Hostonly_ stamped_printf("\tParsing error : Invalid string for 'Unit' in 'Monitor_List'\n");
		Exit(0);
	}
  if     ( !strcasecmp(str, "Dimensional") ) {
    Sampling.unit = DIMENSIONAL;
    M->setSamplingUnit(DIMENSIONAL);
  }
  else if( !strcasecmp(str, "Non_Dimensional") ) {
    Sampling.unit = NONDIMENSIONAL;
    M->setSamplingUnit(NONDIMENSIONAL);
  }
  else {
    Hostonly_ stamped_printf("\tInvalid keyword is described at 'Unit' section\n");
    Exit(0);
  }
  
  // サンプリングの指定単位が有次元の場合に，無次元に変換
  if ( Sampling.unit == DIMENSIONAL ) {
    Interval[Interval_Manager::tg_sampled].normalizeInterval(Tscale);
  }
  
  // モニターリストの読み込み
  elmL2 = elmL1->GetElemFirst();

  while (elmL2) {
    
    // sampling type & param check
    p = elmL2->GetName();
    if ( !strcasecmp(p, "point_set") ) {
      type = MonitorCompo::POINT_SET;
      if ( 0 == elmL2->GetElemSize() ) {
        Hostonly_ stamped_printf("\tParsing error : At least, 1 elem of 'set' should be found in 'point_set'\n");
        Exit(0);
      }
    }
    else if ( !strcasecmp(p, "line") ) {
      type = MonitorCompo::LINE;
      if ( 2 != elmL2->GetElemSize() ) {
        Hostonly_ stamped_printf("\tParsing error : 2 elems (from/to) should be found in 'line'\n");
        Exit(0);
      }
    }
    else {
      Hostonly_ stamped_printf("\tParsing error : No valid keyword [point_set / line] in 'Monitor_List'\n");
      Exit(0);
    }
    
    // Labelの取得．ラベルなしでもエラーではない
    if ( !(label = elmL2->GetComment()) ) {
      Hostonly_ stamped_printf("\tParsing warning : No commnet in '%s'\n", p);
    }
    
    // variable
    variables.clear();
    param = elmL2->GetParamFirst("variable");
    while (param) {
      str=NULL;
      if ( !param->GetData(&str) ) {
        Hostonly_ stamped_printf("\tParsing error : fail to get 'variable' in 'Monitor_List'\n");
        Exit(0);
      }
      variables.push_back(str);
      param = elmL2->GetParamNext(param, "variable");
    }
    if (variables.size() == 0) {
      Hostonly_ stamped_printf("\tParsing error : No 'variable' in 'Monitor_List'\n");
      Exit(0);
    }
    
    // method
    if (!elmL2->GetValue("sampling_method", &method)) {
      Hostonly_ stamped_printf("\tParsing error : fail to get 'sampling_method' in 'Monitor_List'\n");
      Exit(0);
    }
    
    // mode
    if (!elmL2->GetValue("sampling_mode", &mode)) {
      Hostonly_ stamped_printf("\tParsing error : fail to get 'sampling_mode' in 'Monitor_List'\n");
      Exit(0);
    }
    
    // get coordinate
    if ( type == MonitorCompo::POINT_SET ) {
      vector<MonitorCompo::MonitorPoint> pointSet;
      getXML_Mon_Pointset(M, elmL2, pointSet);
      M->setPointSet(label, variables, method, mode, pointSet);
    }
    else {
      REAL_TYPE from[3], to[3];
      int nDivision;
      getXML_Mon_Line(M, elmL2, from, to, nDivision);
      M->setLine(label, variables, method, mode, from, to, nDivision);
    }
    
    elmL2 = elmL1->GetElemNext(elmL2); // ahead on the next pointer
  }
  
}

/**
 @fn CfgElem Control::getXML_Pointer(const char* key, string section)
 @brief section配下のkeyに対するXMLツリーのポインタを返す
 @param key 探索するキーワード
 @param section 対象セクション
 */
const CfgElem* Control::getXML_Pointer(const char* key, string section)
{
 const CfgElem *elemTop=NULL, *elmL1=NULL;
  
  if (section == "steer") {
    elemTop = CF->GetTop(STEER);
  }
  else if (section == "parameter") {
    elemTop = CF->GetTop(PARAMETER);
  }
  else {
    Hostonly_ printf("No section '%s'\n", section.c_str());
    Exit(0);
  }
  
	if( !(elmL1 = elemTop->GetElemFirst(key)) ) {
    Hostonly_ stamped_printf("\tParsing error : Missing the section of '%s'\n", key);
    Exit(0);
  }

  return elmL1;
}

/**
 @fn void Control::getXML_Time_Control(DTcntl* DT)
 @brief 時間制御に関するパラメータを取得する
 @param DT
 @note パラメータは，setParameters()で無次元して保持
 */
void Control::getXML_Time_Control(DTcntl* DT)
{
  const CfgElem *elmL1=NULL;
  REAL_TYPE ct;
  int ss=0;
  const char *str=NULL;
  
  elmL1 = getXML_Pointer("Time_Control", "steer");
  
  // 加速時間
  if ( !elmL1->GetValue("Acceleration_Type", &str) ) {
    Hostonly_ stamped_printf("\tParsing error : fail to get 'Acceleration_Type' in 'Time_Control'\n");
    Exit(0);
  }
  else {
    if ( !strcasecmp(str, "step") ) {
      Interval[Interval_Manager::tg_accelra].setMode_Step();
    }
    else if ( !strcasecmp(str, "time") ) {
      Interval[Interval_Manager::tg_accelra].setMode_Time();
    }
    else {
      Hostonly_ stamped_printf("\tParsing error : Invalid keyword for 'Acceleration_Type' in 'Time_Control'\n");
      Exit(0);
    }
    
    if ( elmL1->GetValue("Acceleration", &ct) ) {
      Interval[Interval_Manager::tg_accelra].setInterval((double)ct);
    }
    else {
      Hostonly_ stamped_printf("\tParsing error : fail to get 'Acceleration' in 'Time_Control'\n");
      Exit(0);
    }
  }
  
  // 時間積分幅を取得する
  if ( !elmL1->GetValue("Dt_Type", &str) ) {
    Hostonly_ stamped_printf("\tParsing error : fail to get 'Dt_Type' in 'Time_Control'\n");
    Exit(0);
  }
  if ( !elmL1->GetValue("Delta_t", &ct) ) {
    Hostonly_ stamped_printf("\tParsing error : fail to get 'Delta_t' in 'Time_Control'\n");
    Exit(0);
  }
  // Directで有次元の場合は，無次元化
  double ts = RefLength / RefVelocity;
  double cc;
  if ( !strcasecmp(str, "Direct") ) {
    if (Unit.Param == DIMENSIONAL) {
      cc = (double)ct / ts;
    }
  }
  else {
    cc = (double)ct;
  }
  
  if ( !DT->set_Scheme(str, cc) ) {
    Hostonly_ stamped_printf("\tParsing error : fail to set DELTA_T\n");
    Exit(0);
  }
  
  // 計算する時間を取得する
  if ( !elmL1->GetValue("Period_Type", &str) ) {
    Hostonly_ stamped_printf("\tParsing error : fail to get 'Period_Type' in 'Time_Control'\n");
    Exit(0);
  }
  else {
    if ( !strcasecmp(str, "step") ) {
      Interval[Interval_Manager::tg_compute].setMode_Step();
    }
    else if ( !strcasecmp(str, "time") ) {
      Interval[Interval_Manager::tg_compute].setMode_Time();
    }
    else {
      Hostonly_ stamped_printf("\tParsing error : Invalid keyword for 'Period_Type' in 'Time_Control'\n");
      Exit(0);
    }
    
    if ( elmL1->GetValue("Calculation_Period", &ct) ) {
      Interval[Interval_Manager::tg_compute].setInterval((double)ct);
    }
    else {
      Hostonly_ stamped_printf("\tParsing error : fail to get 'Calculation_Period' in 'Time_Control'\n");
      Exit(0);
    }
  }
}

/**
 @fn void Control::getXML_Solver_Properties(void)
 @brief ソルバーの種類を特定するパラメータを取得し，ガイドセルの値を決定する
 */
void Control::getXML_Solver_Properties(void)
{
  const CfgElem *elmL1=NULL;
  const char* str=NULL;
  int ct;
  
  // 次元の設定
  NoDimension = 3;
  
  // getXML_VarArrangement() 変数配置（stg / cc / node）を取得，下のガイドセルの値の設定に影響
  Mode.VarArrange = CELL_CENTER;
  
  if ( !(elmL1 = getXML_Pointer("Solver_Property", "steer")) ) Exit(0);
  
  // 形状近似度の取得
  if ( !elmL1->GetValue("Shape_Approximation", &str) ) {
    Hostonly_ stamped_printf("\tParsing error : fail to get 'Shape_Approximation' in 'Solver_Property'\n");
    Exit(0);
  }
  if     ( !strcasecmp(str, "Binary") )       Mode.ShapeAprx = BINARY;
  else if( !strcasecmp(str, "cut_distance") ) Mode.ShapeAprx = CUT_INFO;
  else {
    Hostonly_ stamped_printf("\tInvalid keyword is described for 'Shape_Approximation'\n");
    Exit(0);
  }
  
  // Cut-Distanceの場合
  if ( Mode.ShapeAprx == CUT_INFO ) {
    
    // ポリゴンファイル名を取得
    if ( !elmL1->GetValue("Polylib_Configuration_File", &str) ) {
      Hostonly_ stamped_printf("\tParsing error : Invalid char* value in 'Solver_Property'\n");
      Exit(0);
    }
    strcpy(PolylibConfigName, str);
    
    // 媒質指定モード
    if ( !elmL1->GetValue("Medium_File", &str) ) {
      Hostonly_ stamped_printf("\tParsing error : fail to get 'Medium_File' in 'Solver_Property'\n");
      Exit(0);
    }
    if     ( !strcasecmp(str, "Yes") ) Mode.Medium_Spec = ON;
    else if( !strcasecmp(str, "No") )  Mode.Medium_Spec = OFF;
    else {
      Hostonly_ stamped_printf("\tInvalid keyword is described for 'Medium_File'\n");
      Exit(0);
    }
    
    // 媒質ファイルを使わない場合の媒質番号の指定
    if ( Mode.Medium_Spec == OFF ) {
      if ( !elmL1->GetValue("Base_Medium", &ct) ) {
        Hostonly_ stamped_printf("\tParsing error : Invalid value for 'Base_Medium' in 'Solver_Property'\n");
        Exit(0);
      }
      if (ct<1) {
        Hostonly_ stamped_printf("\tParsing error : Base_Medium must be positive\n");
        Exit(0);
      }
      Mode.Base_Medium = (unsigned)ct;
    }
    
  }
  
  // 支配方程式の型（PDE_NS / Euler）を取得
  if ( !elmL1->GetValue("PDE_type", &str) ) {
    Hostonly_ stamped_printf("\tParsing error : fail to get 'PDE_type' in 'Solver_Property'\n");
    Exit(0);
  }
  if     ( !strcasecmp(str, "Navier_Stokes") )  Mode.PDE = PDE_NS;
  else if( !strcasecmp(str, "Euler") )          Mode.PDE = PDE_EULER;
  else {
    Hostonly_ stamped_printf("\tInvalid keyword is described for 'PDE_type'\n");
    Exit(0);
  }
  
  // 基礎方程式の種類を取得する
  if ( !elmL1->GetValue("Basic_Equation", &str) ) {
    Hostonly_ stamped_printf("\tParsing error : fail to get 'Basic_Equation' in 'Solver_Property'\n");
    Exit(0);
  }
  if     ( !strcasecmp(str, "Incompressible") )           BasicEqs = INCMP;
  else if( !strcasecmp(str, "Limited_Compressible") )     BasicEqs = LTDCMP;
  else if( !strcasecmp(str, "Compressible") )             BasicEqs = CMPRSS;
  else if( !strcasecmp(str, "Incompressible_Two_Phase") ) BasicEqs = INCMP_2PHASE;
  else {
    Hostonly_ stamped_printf("\tInvalid keyword is described for 'Basic_Equation'\n");
    Exit(0);
  }
  
  // 非定常計算，または定常計算の種別を取得する
  if ( !elmL1->GetValue("Time_Variation", &str) ) {
    Hostonly_ stamped_printf("\tParsing error : fail to get 'Time_Variation' in 'Solver_Property'\n");
    Exit(0);
  }
  if     ( !strcasecmp(str, "Steady") )    Mode.Steady = TV_Steady;
  else if( !strcasecmp(str, "Unsteady") )  Mode.Steady = TV_Unsteady;
  else {
    Hostonly_ stamped_printf("\tInvalid keyword is described for 'Time_Variation'\n");
    Exit(0);
  }
  
  // 対流項スキームの種類の取得
  getXML_Convection();
  
  // ソルバーの種類（FLOW_ONLY / THERMAL_FLOW / THERMAL_FLOW_NATURAL / CONJUGATE_HEAT_TRANSFER / SOLID_CONDUCTION）と浮力モード
  getXML_KindOfSolver(elmL1);
  
  // ガイドセルの値を決める getXML_Convection(), getXML_KindOfSolver(), getXML_VarArrangement()のあと
  if (KindOfSolver==SOLID_CONDUCTION) {
    guide = 1;
  }
  else {
    switch (CnvScheme) {
      case O1_upwind:
      case O2_central:
        guide = 1;
        break;
        
      case O3_muscl:
        guide = 2;
        break;
        
      default:
        Exit(0);
    }
  }  
}

/**
 @fn void Control::getXML_Steer_1(DTcntl* DT)
 @brief 制御，計算パラメータ群の取得
 @param DT
 @note 他のパラメータ取得に先んじて処理しておくもの
 */
void Control::getXML_Steer_1(DTcntl* DT)
{

  // ソルバーの具体的な種類を決めるパラメータ（変数配置，形状近似度，次元）を取得し，ガイドセルの値を設定する
  getXML_Solver_Properties();

  // 指定単位が有次元か無次元かを取得
  getXML_Unit();

  // Reference parameter needs to be called before setDomain();
  // パラメータの取得，代表値に関するもの．
  getXML_Para_Ref();

  // 時間制御パラメータ
  getXML_Time_Control(DT);

  // 精度
  if ( sizeof(REAL_TYPE) == sizeof(double) ) {
    Mode.Precision = SPH_DOUBLE;
  }
  else {
    Mode.Precision = SPH_SINGLE;
  }

  // ファイル入出力に関するパラメータ
  getXML_FileIO();

  // パラメータチェック
  getXML_CheckParameter();

  // スケーリングファクタの取得　***隠しパラメータ
  getXML_Scaling();

}

/**
 @fn void Control::getXML_Steer_2(ItrCtl* IC)
 @brief 制御，計算パラメータ群の取得
 */
void Control::getXML_Steer_2(ItrCtl* IC, ReferenceFrame* RF)
{
  // 流体の解法アルゴリズムを取得
  getXML_Algorithm(); 

  // パラメータを取得する．この部分は getXMLSteer()より先に記述しておくこと
  if ( Unit.Param == NONDIMENSIONAL ) {
    getXML_Para_ND();
  }
  
  if ( isHeatProblem() ) {
    getXML_Para_Temp();
  }
  //getXML_Para_Wind();
  getXML_Para_Init();

  // If a start section dose not describe, the initial start is assumed.
  StartType type = CF->GetStartType();
  switch (type) {
    case Initial:
      Start = initial_start;
      break;
      
    case Restart:
      Start = re_start;
      break;
      
    default:
      Hostonly_ stamped_printf("\tParsing error : Start section\n");
      Exit(0);
  }

  // Reference frame information : Solver defined element >> SPHERE defined
  getXML_ReferenceFrame(RF);

  // 時間平均操作
  getXML_Average_option();

  // 圧力ノイマン条件のタイプ >> getXML_Log()よりも先に
  getXML_Wall_type();
  
  // Log >> getXML_Iteration()よりも前に
  getXML_Log();
  
  // Criteria of computation
  getXML_Iteration(IC);
  
  // LES
  getXML_LES_option();
  
  // 派生変数のオプション
  getXML_Derived();
  
  // 変数範囲の処理　***隠しパラメータ
  getXML_VarRange();
  
  // Cell IDのゼロを指定IDに変更　***隠しパラメータ
  getXML_ChangeID();
  
  // 性能測定モードの処理　***隠しパラメータ
  getXML_PMtest();
}

/**
 @fn void Control::displayParams(FILE* mp, FILE* fp, ItrCtl* IC, DTcntl* DT, ReferenceFrame* RF)
 @brief 制御，計算パラメータ群の表示
 @param mp ファイルポインタ（標準出力）
 @param fp ファイルポインタ（ファイル出力）
  */
void Control::displayParams(FILE* mp, FILE* fp, ItrCtl* IC, DTcntl* DT, ReferenceFrame* RF)
{
  printSteerConditions(mp, IC, DT, RF);
  printSteerConditions(fp, IC, DT, RF);
  printParaConditions(mp);
  printParaConditions(fp);
  printInitValues(mp);
  printInitValues(fp);
}

/**
 @fn REAL_TYPE Control::OpenDomainRatio(unsigned dir, REAL_TYPE area, const unsigned Dims, unsigned* G_size)
 @brief 外部境界の各方向の開口率（流体部分の比率）
 @retval 開口率
 @param dir 方向
 @param area 計算外部領域の各面の開口率
 @param Dims 次元数
 @param G_Size 計算領域全体の分割数
 */
REAL_TYPE Control::OpenDomainRatio(unsigned dir, REAL_TYPE area, const unsigned Dims, unsigned* G_size)
{
  REAL_TYPE r = 0.0, base=0.0;
  unsigned m_imax, m_jmax, m_kmax;
  
  m_imax = G_size[0];
  m_jmax = G_size[1];
  m_kmax = G_size[2];

  switch (dir) {
    case X_MINUS:
    case X_PLUS:
      (2==Dims) ? base=(REAL_TYPE)m_jmax : base=(REAL_TYPE)(m_jmax*m_kmax);
      r = area / base*100.0;
      break;
      
    case Y_MINUS:
    case Y_PLUS:
      (2==Dims) ? base=(REAL_TYPE)m_imax : base=(REAL_TYPE)(m_imax*m_kmax);
      r = area / base*100.0;
      break;
      
    case Z_MINUS:
    case Z_PLUS:
      (2==Dims) ? base=1.0 : base=(REAL_TYPE)(m_imax*m_jmax);
      r = area / base*100.0;
      break;
  }
  return (r);
}

/**
 @fn void Control::printAreaInfo(FILE* fp, FILE* mp, unsigned G_Fcell, unsigned G_Acell, unsigned* G_size)
 @brief 有効セル数を表示する
 @param fp 出力ファイルポインタ
 @param mp 標準出力ファイルポインタ
 @param G_Fcell グローバルなFluid cell
 @param G_Acell グローバルなActive cell
 @param G_size global size
 */
void Control::printAreaInfo(FILE* fp, FILE* mp, unsigned G_Fcell, unsigned G_Acell, unsigned* G_size)
{
  printArea(mp, G_Fcell, G_Acell, G_size);
  printArea(fp, G_Fcell, G_Acell, G_size);
}

/**
 @fn void Control::printInitValues(FILE* fp)
 @brief 初期値の表示
 @note
    - Init*には無次元値が保持されている
 @see void Control::setInitialConditions(void)
 */
void Control::printInitValues(FILE* fp)
{
  if ( !fp ) {
    stamped_printf("\tFail to write into file\n");
    Exit(0);
  }

  REAL_TYPE DynamicPrs = 0.5*RefDensity * RefVelocity * RefVelocity;

  fprintf(fp,"\n---------------------------------------------------------------------------\n\n");
  fprintf(fp,"\n\t>> Initial Values for Physical Variables\n\n");

  fprintf(fp,"\tInitial  Density     [kg/m^3]/ [-]   : %12.5e / %12.5e\n", iv.Density, iv.Density/RefDensity);
  fprintf(fp,"\tInitial  Velocity.U  [m/s]   / [-]   : %12.5e / %12.5e\n", iv.VecU,    iv.VecU/RefVelocity);
  fprintf(fp,"\tInitial          .V  [m/s]   / [-]   : %12.5e / %12.5e\n", iv.VecV,    iv.VecV/RefVelocity);
  fprintf(fp,"\tInitial          .W  [m/s]   / [-]   : %12.5e / %12.5e\n", iv.VecW,    iv.VecW/RefVelocity);
  fprintf(fp,"\tDynamic  Pressure    [Pa]    / [-]   : %12.5e / %12.5e\n", DynamicPrs,             1.0);
  if (Unit.Prs == CompoList::Absolute) {
    fprintf(fp,"\tInitial  Pressure    [Pa]    / [-]   : %12.5e / %12.5e\n", iv.Pressure, (iv.Pressure-BasePrs)/DynamicPrs);
  }
  else {
    fprintf(fp,"\tInitial  Pressure    [Pa_g]  / [-]   : %12.5e / %12.5e\n", iv.Pressure, iv.Pressure/DynamicPrs);
  }
  if ( isHeatProblem() ) {
    fprintf(fp,"\tInitial  Temperature [%s]     / [-]   : %12.5e / %12.5e\n", (Unit.Temp==CompoList::Unit_KELVIN) ? "K" : "C",
						FBUtility::convK2Temp(iv.Temperature, Unit.Temp), 
						FBUtility::convK2ND(iv.Temperature, BaseTemp, DiffTemp));
  }

  fprintf(fp,"\n");
  fflush(fp);
}

/**
 @fn void Control::printNoCompo(FILE* fp)
 @brief 内部BCコンポーネントの数を表示する
 */
void Control::printNoCompo(FILE* fp)
{
  fprintf(fp,"\tNo. of Inner Boundary  : %d\n", NoBC);
  fprintf(fp,"\tNo. of Medium          : %d\n", NoMaterial);
  fprintf(fp,"\n");
  fprintf(fp,"\tNo. of Fluid ID        : %d\n", NoMediumFluid);
  fprintf(fp,"\tNo. of Solid ID        : %d\n", NoMediumSolid);
}

/**
 @fn bool Control::receiveCfgPtr(SklSolverConfig* cfg)
 @brief コンフィギュレーションのポインタを返す
 */
bool Control::receiveCfgPtr(SklSolverConfig* cfg)
{
  if ( !cfg ) return false;
  CF = cfg;
  return true;
}

/**
 @fn void Control::getXML_Algorithm(void)
 @brief 解法アルゴリズムを選択する
 */
void Control::getXML_Algorithm(void)
{
  const CfgElem *elmL1=NULL;
  const char* str=NULL;

  if ( !(elmL1 = getXML_Pointer("Algorithm", "steer")) ) Exit(0);

  // Flow
  if ( !elmL1->GetValue("Flow", &str) ) {
    Hostonly_ stamped_printf("\tParsing error : fail to get 'Flow' in 'Algorithm'\n");
    Exit(0);
  }
  if     ( !strcasecmp(str, "FS_C_EE_D_EE") )     AlgorithmF = Flow_FS_EE_EE;
  else if( !strcasecmp(str, "FS_C_RK_D_CN") )     AlgorithmF = Flow_FS_RK_CN;
  else if( !strcasecmp(str, "FS_C_AB_D_AB") )     AlgorithmF = Flow_FS_AB2;
  else if( !strcasecmp(str, "FS_C_AB_D_CN") )     AlgorithmF = Flow_FS_AB_CN;
  else {
    Hostonly_ stamped_printf("\tParsing error : Invalid keyword for 'Flow' in 'Algorithm'\n");
    Exit(0);
  }
  
  // Heat
  if ( isHeatProblem() ) {
    if ( !elmL1->GetValue("Heat", &str) ) {
      Hostonly_ stamped_printf("\tParsing error : fail to get 'Heat' in 'Algorithm'\n");
      Exit(0);
    }
    if     ( !strcasecmp(str, "C_EE_D_EE") )    AlgorithmH = Heat_EE_EE;
    else if( !strcasecmp(str, "C_EE_D_EI") )    AlgorithmH = Heat_EE_EI;
    else {
      Hostonly_ stamped_printf("\tParsing error : Invalid keyword for 'Heat' in 'Algorithm'\n");
      Exit(0);
    }
  }
}

/**
 @fn void Control::getXML_Convection(void)
 @brief 対流項スキームのパラメータを取得する
 */
void Control::getXML_Convection(void)
{
	const CfgElem *elmL1=NULL;
	const char *str=NULL;
	
  if ( !(elmL1 = getXML_Pointer("Convection_Term", "steer")) ) Exit(0);
	
	// scheme
	if ( !elmL1->GetValue("scheme", &str) ) {
		Hostonly_ stamped_printf("\tParsing error : Invalid char* value for 'Scheme' in 'Convection'\n");
		Exit(0);
	}
	
	if     ( !strcasecmp(str, "O1_Upwind") )    CnvScheme = O1_upwind;
  else if( !strcasecmp(str, "O3_muscl") )     CnvScheme = O3_muscl;
  else if( !strcasecmp(str, "O2_central") )   CnvScheme = O2_central;
  else if( !strcasecmp(str, "O4_central") ) { CnvScheme = O4_central; Exit(0); }  // not yet implemented
  else {
    Hostonly_ stamped_printf("\tInvalid keyword is described for Scheme\n");
    Exit(0);
  }
	
	// Limiter
	if ( CnvScheme == O3_muscl ) {
		if ( !elmL1->GetValue("limiter", &str) ) {
			Hostonly_ stamped_printf("\tParsing error : Invalid char* value for 'Limiter' in 'Convection'\n");
			Exit(0);
		}
		if     ( !strcasecmp(str, "No_Limiter") ) Limiter = No_Limiter;
		else if( !strcasecmp(str, "Minmod") )     Limiter = MinMod;
		else {
			Hostonly_ stamped_printf("\tInvalid keyword is described for Limiter\n");
			Exit(0);
		}
	}
}

/**
 @fn void Control::setDomainInfo(unsigned* m_sz, REAL_TYPE* m_org, REAL_TYPE* m_pch, REAL_TYPE* m_wth)
 @brief 無次元の領域情報をセットする
 @param m_sz 領域分割数（計算領域内部のみ）
 @param m_org 基点
 @param m_pch 格子幅
 @param m_wth 領域の大きさ
 @pre Control::setGiudeCell()
 @note
    - setGiudeCell()の前にコール
 */
void Control::setDomainInfo(unsigned* m_sz, REAL_TYPE* m_org, REAL_TYPE* m_pch, REAL_TYPE* m_wth)
{
  // set parameters
  imax = m_sz[0];
  jmax = m_sz[1];
  kmax = m_sz[2];

  dh    = m_pch[0];
  dx[0] = m_pch[0];
  dx[1] = m_pch[1];
  dx[2] = m_pch[2];

  org[0] = m_org[0];
  org[1] = m_org[1];
  org[2] = m_org[2];

  Lbx[0] = m_wth[0];
  Lbx[1] = m_wth[1];
  Lbx[2] = m_wth[2];
}

/**
 @fn void Control::getXML_VarArrangement(void)
 @brief 変数配置を取得
 @pre Control::setGiudeCell()
 @note
    - setGiudeCell()の前にコール
 */
void Control::getXML_VarArrangement(void)
{
  const char *keyword=NULL;
  ParseSteer Tree(CF);

  if ( !(keyword=Tree.getParam("Grid_System")) ) Exit(0);

  if     ( !strcasecmp(keyword, "staggered") )   Mode.VarArrange = STAGGERED;
  else if( !strcasecmp(keyword, "collocated") )  Mode.VarArrange = CELL_CENTER;
  else if( !strcasecmp(keyword, "node") )        Mode.VarArrange = NODE;
  else {
    Hostonly_ stamped_printf("\tInvalid keyword is described for Grid_System\n");
    Exit(0);
  }
}

/**
 @fn void Control::getXML_Para_Ref(void)
 @brief 参照パラメータを取得
 @note Ref_IDで指定される媒質を代表物性値とする
 */
void Control::getXML_Para_Ref(void)
{
  ParsePara Tree(CF);

  if ( !Tree.IsSetElem("Reference") ) Exit(0);

  if ( !Tree.getEParam("Length",    RefLength))   Exit(0);
  if ( !Tree.getEParam("Velocity",  RefVelocity)) Exit(0);
  if ( !Tree.getEParam("Gravity",   Gravity))     Exit(0);
	if ( !Tree.getEParam("Base_Pressure", BasePrs)) Exit(0);
  if ( !Tree.getEParam("Ref_ID", RefID) )         Exit(0);
}

/**
 @fn void Control::getXML_Para_ND(void)
 @brief 無次元パラメータを各種モードに応じて設定する
 @note
    - 純強制対流　有次元　（代表長さ，代表速度，動粘性係数，温度拡散係数）
    -           無次元　（Pr, Re > RefV=RefL=1）
 @see bool Control::setParameters(MaterialList* mat, CompoList* cmp)
 */
void Control::getXML_Para_ND(void)
{
  ParsePara Tree(CF);

  if ( !Tree.IsSetElem("Reference") ) Exit(0);

 if (KindOfSolver==FLOW_ONLY) {
    if ( !Tree.getEParam("Reynolds", Reynolds)) Exit(0);
    if ( !Tree.getEParam("Prandtl",  Prandtl))  Exit(0);
 }
}

/**
 @fn void Control::getXML_Para_Temp(void)
 @brief 温度の参照パラメータを取得
 */
void Control::getXML_Para_Temp(void)
{
  REAL_TYPE Base, Diff;

  const CfgElem *elmL1=NULL;
  
  if ( !(elmL1 = getXML_Pointer("Temperature", "parameter")) ) Exit(0);
  
  if( !elmL1->GetValue("Base", &Base) ) {
    Hostonly_ stamped_printf("\tParsing error : Invalid float value for 'Base'\n");
    Exit(0);
  }
  if( !elmL1->GetValue("Difference", &Diff) ) {
    Hostonly_ stamped_printf("\tParsing error : Invalid float value for 'Difference'\n");
    Exit(0);
  }

  if( Diff < 0.0f ){
    Hostonly_ stamped_printf("\tTemperature difference must be positive.\n");
    Exit(0);
  }
  DiffTemp = Diff;
  
  if ( Unit.Temp == CompoList::Unit_CELSIUS ) {
    BaseTemp = Base + KELVIN;
  }
}

/**
 @fn void Control::printDomainInfo(FILE* fp, FILE* mp, unsigned* G_size, REAL_TYPE* G_org, REAL_TYPE* G_Lbx)
 @brief 領域情報をfp, mpに出力する
 @param fp file pointer
 @param G_size global size
 @param G_org global original point
 @param G_Lbx global bbox of computational domain
 */
void Control::printDomainInfo(FILE* mp, FILE* fp, unsigned* G_size, REAL_TYPE* G_org, REAL_TYPE* G_Lbx)
{
  if ( !fp || !mp ) {
    stamped_printf("\tFail to write into file\n");
    Exit(0);
  }
  printDomain(mp, G_size, G_org, G_Lbx);
  printDomain(fp, G_size, G_org, G_Lbx);
}

/**
 @fn void Control::printDomain(FILE* fp, unsigned* G_size, REAL_TYPE* G_org, REAL_TYPE* G_Lbx)
 @brief グローバルな領域情報を表示する
 @param fp file pointer
 @param G_size global size
 @param G_org global original point
 @param G_Lbx global bbox of computational domain
 */
void Control::printDomain(FILE* fp, unsigned* G_size, REAL_TYPE* G_org, REAL_TYPE* G_Lbx)
{
  fprintf(fp,"\timax, jmax, kmax    = %13d %13d %13d     >> ", G_size[0], G_size[1], G_size[2]);
  printVoxelSize(G_size, fp);
  fprintf(fp,"\n");
  fprintf(fp,"\t(dx, dy, dz)  [m] / [-] = (%13.6e %13.6e %13.6e)  /  (%13.6e %13.6e %13.6e)\n",    dx[0]*RefLength,    dx[1]*RefLength,    dx[2]*RefLength,    dx[0],    dx[1],    dx[2]);
  fprintf(fp,"\t(ox, oy, oz)  [m] / [-] = (%13.6e %13.6e %13.6e)  /  (%13.6e %13.6e %13.6e)\n", G_org[0]*RefLength, G_org[1]*RefLength, G_org[2]*RefLength, G_org[0], G_org[1], G_org[2]);
  fprintf(fp,"\t(Lx, Ly, Lz)  [m] / [-] = (%13.6e %13.6e %13.6e)  /  (%13.6e %13.6e %13.6e)\n", G_Lbx[0]*RefLength, G_Lbx[1]*RefLength, G_Lbx[2]*RefLength, G_Lbx[0], G_Lbx[1], G_Lbx[2]);
  fprintf(fp,"\n");
  fflush(fp);
}

//@fn void Control::printDomain_debug(void)
//@brief デバッグのため，領域情報を表示する
void Control::printDomain_debug(void)
{
  printf("%d\timax, jmax, kmax    = %13d %13d %13d\n", pn.ID, imax, jmax, kmax);
  printf("%d\t(dx, dy, dz)  [m] / [-] = (%13.6e %13.6e %13.6e)  /  (%13.6e %13.6e %13.6e)\n", pn.ID,  dx[0]*RefLength,  dx[1]*RefLength,  dx[2]*RefLength,  dx[0],  dx[1],  dx[2]);
  printf("%d\t(ox, oy, oz)  [m] / [-] = (%13.6e %13.6e %13.6e)  /  (%13.6e %13.6e %13.6e)\n", pn.ID, org[0]*RefLength, org[1]*RefLength, org[2]*RefLength, org[0], org[1], org[2]);
  printf("%d\t(Lx, Ly, Lz)  [m] / [-] = (%13.6e %13.6e %13.6e)  /  (%13.6e %13.6e %13.6e)\n", pn.ID, Lbx[0]*RefLength, Lbx[1]*RefLength, Lbx[2]*RefLength, Lbx[0], Lbx[1], Lbx[2]);
  printf("\n");
  fflush(stdout);
}

/**
 @fn void Control::MemoryRequirement(char* mode, long Memory, long l_memory, FILE* fp)
 @brief メモリ使用量を表示する
 @param mode 処理モード（前処理 or ソルバー）
 @param Memory 必要メモリ量
 @param l_memory local
 */
void Control::printVoxelSize(unsigned* gs, FILE* fp)
{
  if( !fp ) Exit(0);
  
  REAL_TYPE PB=0.0, TB=0.0, GB=0.0, MB=0.0, KB=0.0, total=0.0;
  KB = 1000.0;
  MB = 1000.0*KB;
  GB = 1000.0*MB;
  TB = 1000.0*GB;
  PB = 1000.0*TB;
  
  total = (REAL_TYPE)gs[0] * (REAL_TYPE)gs[1] * (REAL_TYPE)gs[2];

  if ( total > PB ) {
    fprintf (fp,"%6.2f (P cells)\n", total / PB);
  }
  else if ( total > TB ) {
    fprintf (fp,"%6.2f (T cells)\n", total / TB);
  }
  else if ( total > GB ) {
    fprintf (fp,"%6.2f (G cells)\n", total / GB);
  }
  else if ( total > MB ) {
    fprintf (fp,"%6.2f (M cells)\n", total / MB);
  }
  else if ( total > KB ) {
    fprintf (fp,"%6.2f (K cells)\n", total / KB);
  }
  else if ( total <= KB ){
    fprintf (fp,"%6.2f (cells)\n", total);
  }
}

/**
 @fn void Control::findXMLCriteria(const CfgElem *elmL1, const char* key, unsigned order, ItrCtl* IC)
 @brief 反復の収束判定パラメータを取得
 @param elmL1 Iteration_Flow/Heat　レベル
 @param key キーワード
 @param order ItrCtl配列の格納番号
 @param IC 反復制御用クラスの配列
 */
void Control::findXMLCriteria(const CfgElem *elmL1, const char* key, unsigned order, ItrCtl* IC)
{
  const CfgElem* elmL2=NULL;
  const char *str=NULL, *slvr=NULL;
  int itr=0;
  REAL_TYPE tmp=0.0;
  unsigned LinearSolver=0;

  if ( (elmL2 = elmL1->GetElemFirst(key)) ) {
    if( !elmL2->GetValue("Iteration", &itr) ) {
      Hostonly_ stamped_printf("\tParsing error : Invalid integer value for 'Iteration' of %s in Criteria\n", key);
      Exit(0);
    }
    IC[order].set_ItrMax((unsigned)itr);
    
    if( !elmL2->GetValue("Epsilon", &tmp) ) {
      Hostonly_ stamped_printf("\tParsing error : Invalid float value for 'Epsilon' of %s in Criteria\n", key);
      Exit(0);
    }
    IC[order].set_eps((REAL_TYPE)tmp);
    
    if( !elmL2->GetValue("Omega", &tmp) ) {
      Hostonly_ stamped_printf("\tParsing error : Invalid float value for 'Omega' of %s in Criteria\n", key);
      Exit(0);
    }
    IC[order].set_omg((REAL_TYPE)tmp);
    
		if( !elmL2->GetValue("norm", &str) ) {
			Hostonly_ stamped_printf("\tParsing error : Invalid char* value for 'Norm' of %s in Criteria\n", key);
			Exit(0);
    }
    
    if( !elmL2->GetValue("Linear_Solver", &slvr) ) {
			Hostonly_ stamped_printf("\tParsing error : Invalid char* value for 'Linear_Solver' of %s in Criteria\n", key);
			Exit(0);
    }
  }
  else {
    Hostonly_ stamped_printf("\tParsing error : Invalid keyword of '%s' in Iteration_Flow/Heat\n", key);
    Exit(0);
  }
  
  // 線形ソルバーの種類
  if     ( !strcasecmp(slvr, "Jacobi") )    IC[order].set_LS(JACOBI);
  else if( !strcasecmp(slvr, "SOR") )       IC[order].set_LS(SOR);
  else if( !strcasecmp(slvr, "SOR2SMA") )   IC[order].set_LS(SOR2SMA);
  else if( !strcasecmp(slvr, "SOR2CMA") )   IC[order].set_LS(SOR2CMA);
  else {
    Hostonly_ stamped_printf("\tInvalid keyword is described for Linear_Solver\n");
    Exit(0);
  }
  
  // normのタイプ
	switch (order) {
		case ItrCtl::ic_prs_pr: // Predictor phase
		case ItrCtl::ic_prs_cr: // Corrector phase
      if (Mode.Log_Itr == ON) {
        IC[order].set_normType(ItrCtl::v_div_max_dbg);
      }
      else {
        if ( !strcasecmp(str, "v_div_max") ) {
          IC[order].set_normType(ItrCtl::v_div_max);
        }
        else if ( !strcasecmp(str, "v_div_L2") ) {
          IC[order].set_normType(ItrCtl::v_div_l2);
        }
        else if ( !strcasecmp(str, "p_res_L2_absolute") ) {
          IC[order].set_normType(ItrCtl::p_res_l2_a);
        }
        else if ( !strcasecmp(str, "p_res_L2_relative") ) {
          IC[order].set_normType(ItrCtl::p_res_l2_r);
        }
        else {
          Hostonly_ stamped_printf("\tParsing error : Invalid keyword for '%s' of Norm for Poisson iteration\n", str);
          Exit(0);
        }
      }
			break;
			
		case ItrCtl::ic_tdf_ei: // Temperature Euler Implicit
			if ( !strcasecmp(str, "t_res_L2_absolute") ) {
				IC[order].set_normType(ItrCtl::t_res_l2_a);
			}
			else if ( !strcasecmp(str, "t_res_L2_relative") ) {
				IC[order].set_normType(ItrCtl::t_res_l2_r);
			}
			else {
				Hostonly_ stamped_printf("\tParsing error : Invalid keyword for '%s' of Norm for heat iteration\n", str);
				Exit(0);
			}
			break;
      
    case ItrCtl::ic_vis_cn: // Velocity Crank-Nicolosn
			if ( !strcasecmp(str, "v_res_L2_relative") ) {
				IC[order].set_normType(ItrCtl::v_res_l2_r);
			}
			else {
				Hostonly_ stamped_printf("\tParsing error : Invalid keyword for '%s' of Norm for heat iteration\n", str);
				Exit(0);
			}
			break;
	}
}

/**
 @fn void Control::getXML_ReferenceFrame(ReferenceFrame* RF)
 @brief 参照座標系を取得する
 @todo 回転は未
 */
void Control::getXML_ReferenceFrame(ReferenceFrame* RF)
{
  const CfgElem *elemTop=NULL, *elmL1=NULL;
  const char *str=NULL;
  
  elemTop = CF->GetTop(STEER);

  if( !(elmL1 = elemTop->GetElemFirst("Reference_Frame")) ) {
    Hostonly_ stamped_printf("\tParsing error : Missing the section of Reference_Frame\n");
    Exit(0);
  }
  
  if ( !elmL1->GetValue("Reference_Frame_Type", &str) ) {
    Hostonly_ stamped_printf("\tParsing error : fail to get 'Reference_Frame_Type' in 'Reference_Frame'\n");
    Exit(0);
  }

  if ( !strcasecmp(str, "stationary") ) {
    RF->setFrame(ReferenceFrame::frm_static);
  }
  else if ( !strcasecmp(str, "translational") ) {
    RF->setFrame(ReferenceFrame::frm_translation);
    REAL_TYPE xyz[3];
    if ( !SklUtil::getVecParams(elmL1, "u", "v", "w", xyz) ) {
      Hostonly_ stamped_printf("\tParsing error : Missing param 'u', 'v', 'w' or Invalid values for Translational Reference_Frame\n");
      Exit(0);
    }
    RF->setGridVel((double*)xyz);
  }
  else {
    Hostonly_ stamped_printf("\tParsing error : Invalid keyword for 'Reference_Frame_Type' in 'Reference_Frame'\n");
    Exit(0);
  }
}

/**
 @fn void Control::getXML_Iteration(ItrCtl* IC)
 @brief 反復関連の情報を取得する
 */
void Control::getXML_Iteration(ItrCtl* IC)
{
  const CfgElem *elemTop=NULL, *elmL1=NULL, *elmL2=NULL;

  // Iteration
  elemTop = CF->GetTop(STEER);
  if( !(elmL1 = elemTop->GetElemFirst("Iteration")) ) {
    Hostonly_ stamped_printf("\tParsing error : Missing the section of 'Iteration'\n");
    Exit(0);
  }

  // Flow
  if( !(elmL2 = elmL1->GetElemFirst("Flow")) ) {
    Hostonly_ stamped_printf("\tParsing error : Missing the section of 'Flow' in 'Iteration'\n");
    Exit(0);
  }
  
  switch (AlgorithmF) {
    case Flow_FS_EE_EE:
    case Flow_FS_AB2:
      if( elmL2->GetElemSize() != 1 ) {  // check number of Elem
        Hostonly_ stamped_printf("\tOne criterion should be specified for 1st order.\n");
        Exit(0);
      }
      findXMLCriteria(elmL2, "Poisson", ItrCtl::ic_prs_pr, IC);
      break;

    case Flow_FS_AB_CN:
      if( elmL2->GetElemSize() != 2 ) {  // check number of Elem
        Hostonly_ stamped_printf("\tTwo criterions should be specified\n");
        Exit(0);
      }
      findXMLCriteria(elmL2, "Poisson", ItrCtl::ic_prs_pr, IC);
      findXMLCriteria(elmL2, "NS_CN",   ItrCtl::ic_vis_cn, IC);
      break;
      
    case Flow_FS_RK_CN:
      if( elmL2->GetElemSize() != 2 ) {  // check number of Elem
        Hostonly_ stamped_printf("\tTwo criterions should be specified for 2nd order.\n");
        Exit(0);
      }
      findXMLCriteria(elmL2, "Poisson",     ItrCtl::ic_prs_pr, IC);
      findXMLCriteria(elmL2, "Poisson_2nd", ItrCtl::ic_prs_cr, IC);
      findXMLCriteria(elmL2, "NS_CN",       ItrCtl::ic_vis_cn, IC);
      break;

    default:
      Hostonly_ stamped_printf("\tSomething wrong in 'Iteration' > 'Flow'\n");
      Exit(0);
  }

  // Heat
  if ( isHeatProblem() ) {
    
    if( !(elmL2 = elmL1->GetElemFirst("Heat")) ) {
      Hostonly_ stamped_printf("\tParsing error : Missing the section of 'Heat' in 'Iteration'\n");
      Exit(0);
    }
    
    switch (AlgorithmH) {
      case Heat_EE_EE:
        break;

      case Heat_EE_EI:
        if( elmL2->GetElemSize() != 1 ) {  // check number of Elem
          Hostonly_ stamped_printf("\tOne criteria should be specified for Euler Implicit (Temperature).\n");
          Exit(0);
        }
        findXMLCriteria(elmL2, "Euler_Implicit", ItrCtl::ic_tdf_ei, IC);
        break;

      default:
      Hostonly_ stamped_printf("\tSomething wrong in Iteration_Heat\n");
      Exit(0);
    }
  }
}

/**
 @fn void Control::select_Itr_Impl(ItrCtl* IC)
 @brief choose implementation of SOR
 */
void Control::select_Itr_Impl(ItrCtl* IC)
{
  ItrCtl* ICp1 = &IC[ItrCtl::ic_prs_pr];  /// 圧力のPoisson反復 Euler陽解法, RKの予測フェイズ
  ItrCtl* ICp2 = &IC[ItrCtl::ic_prs_cr];  /// 圧力のPoisson反復 RKの修正フェイズ
  ItrCtl* ICv  = &IC[ItrCtl::ic_vis_cn];  /// 粘性項のCrank-Nicolson反復
  
  // Flow
  switch (AlgorithmF) {
    case Flow_FS_EE_EE:
    case Flow_FS_AB2:
      if ( (ICp1->get_LS() == SOR) || (ICp1->get_LS() == SOR2SMA) || (ICp1->get_LS() == JACOBI) ) {
        ICp1->set_LoopType( (Eff_Cell_Ratio < THRESHOLD_SOR_IMPLEMENTATION) ? SKIP_LOOP : MASK_LOOP );
      }
      break;
      
    case Flow_FS_RK_CN:
      if ( (ICp1->get_LS() == SOR) || (ICp1->get_LS() == SOR2SMA) || (ICp1->get_LS() == JACOBI) ) {
        ICp1->set_LoopType( (Eff_Cell_Ratio < THRESHOLD_SOR_IMPLEMENTATION) ? SKIP_LOOP : MASK_LOOP );
      }
      if ( (ICp2->get_LS() == SOR) || (ICp2->get_LS() == SOR2SMA) || (ICp2->get_LS() == JACOBI) ) {
        ICp2->set_LoopType( (Eff_Cell_Ratio < THRESHOLD_SOR_IMPLEMENTATION) ? SKIP_LOOP : MASK_LOOP );
      }
      break;
      
    case Flow_FS_AB_CN:
      if ( (ICp1->get_LS() == SOR) || (ICp1->get_LS() == SOR2SMA) || (ICp1->get_LS() == JACOBI) ) {
        ICp1->set_LoopType( (Eff_Cell_Ratio < THRESHOLD_SOR_IMPLEMENTATION) ? SKIP_LOOP : MASK_LOOP );
      }
      if ( (ICv->get_LS() == SOR) || (ICv->get_LS() == SOR2SMA) || (ICv->get_LS() == JACOBI) ) {
        ICv->set_LoopType( (Eff_Cell_Ratio < THRESHOLD_SOR_IMPLEMENTATION) ? SKIP_LOOP : MASK_LOOP );
      }
      break;
      
    default:
      Hostonly_ stamped_printf("\tSomething wrong in Iteration_Flow\n");
      Exit(0);
  }
  
  // Heat
  if ( isHeatProblem() ) {
    ItrCtl* ICt = &IC[ItrCtl::ic_tdf_ei];  /// 温度の拡散項の反復
    switch (AlgorithmH) {
      case Heat_EE_EE:
        break;
        
      case Heat_EE_EI:
        if ( (ICt->get_LS() == SOR) || (ICt->get_LS() == SOR2SMA) || (ICt->get_LS() == JACOBI) ) {
          ICt->set_LoopType( (Eff_Cell_Ratio < THRESHOLD_SOR_IMPLEMENTATION) ? SKIP_LOOP : MASK_LOOP );
        }
        break;
        
      default:
        Hostonly_ stamped_printf("\tSomething wrong in Iteration_Heat\n");
        Exit(0);
    }
  }
}

/**
 @fn void Control::getXML_History(void)
 @brief ヒストリファイルの情報を取得する
 @pre getXML_Log()
 */
void Control::getXML_History(void)
{
  if ( Mode.Log_Base == ON ) {
    strcpy(HistoryName,      "history_base.txt");
    strcpy(HistoryCompoName, "history_compo.txt");
    strcpy(HistoryDomfxName, "history_domainflux.txt");
  }
  
  if ( Mode.Log_Wall == ON ) {
    strcpy(HistoryWallName, "history_log_wall.txt");
  }
  
  if ( Mode.Log_Itr == ON ) {
    strcpy(HistoryItrName, "history_iteration.txt");
  }
  
  /* >> HIstoryファイルの名前は固定に変更 2011.12.3
  unsigned NoHistoryList = CF->GetHistorySize();

  if ( NoHistoryList > 5 ){
    printf("Only 5 history-files can be specified.\n");
  }
	
  const SklCfgHistory* history = CF->GetHistoryFirst();
  while( history ){
    const char* history_fname = NULL;
    const char* tt = history->GetAttr();
    
    if ( Mode.Log_Base == ON ) {
      if ( !strcasecmp( tt, "log_base" ) ) {
        if ( !(history_fname = history->GetFileName()) ) {
          Hostonly_ stamped_printf("\tParsing error : Filename in HistoryFile description\n");
          Exit(0);
        }
        strcpy(HistoryName, history_fname);
      }
      else if ( !strcasecmp( tt, "log_compo" ) ) {
        if ( !(history_fname = history->GetFileName()) ) {
          Hostonly_ stamped_printf("\tParsing error : Filename in HistoryFile description\n");
          Exit(0);
        }
        strcpy(HistoryCompoName, history_fname);
      }
      else if ( !strcasecmp( tt, "log_domainflux" ) ) {
        if ( !(history_fname = history->GetFileName()) ) {
          Hostonly_ stamped_printf("\tParsing error : Filename in HistoryFile description\n");
          Exit(0);
        }
        strcpy(HistoryDomfxName, history_fname);
      }
    }
    
    if ( Mode.Log_Wall == ON ) {
      if ( !strcasecmp( tt, "log_wall_info" ) ) {
        if ( !(history_fname = history->GetFileName()) ) {
          Hostonly_ stamped_printf("\tParsing error : Filename in HistoryFile description\n");
          Exit(0);
        }
        strcpy(HistoryWallName, history_fname);
      }
    }
    
    if ( Mode.Log_Itr == ON ) {
      if ( !strcasecmp( tt, "log_iteration" ) ) {
        if ( !(history_fname = history->GetFileName()) ) {
          Hostonly_ stamped_printf("\tParsing error : Filename in HistoryFile description\n");
          Exit(0);
        }
        strcpy(HistoryItrName, history_fname);
      }
    }

    history = CF->GetHistoryNext(history);
  }
  */
}

/**
 @fn void Control::getXML_Dimension(void)
 @brief 計算する次元数を取得する
 */
void Control::getXML_Dimension(void)
{
  const char *keyword=NULL;
  ParseSteer Tree(CF);

  if ( !(keyword=Tree.getParam("Dimension")) ) Exit(0);

  if     ( !strcasecmp(keyword, "3D") )  NoDimension = 3;
  else if( !strcasecmp(keyword, "2D") )  NoDimension = 2;
  else if( !strcasecmp(keyword, "1D") )  NoDimension = 1;
  else if( !strcasecmp(keyword, "0D") )  NoDimension = 0;
  else {
    Hostonly_ stamped_printf("\tInvalid keyword is described for Dimension\n");
    Exit(0);
  }
}

/**
 @fn void Control::getXML_Average_option(void)
 @brief 平均値操作に関するパラメータを取得する
 @note パラメータは，setParameters()で無次元して保持
 */
void Control::getXML_Average_option(void)
{
  const CfgElem *elmL1=NULL;
  const char* str=NULL;
  REAL_TYPE ct;
  
  if ( !(elmL1 = getXML_Pointer("Average_option", "steer")) ) Exit(0);
  
  // 平均値操作
  if ( !elmL1->GetValue("Operation", &str) ) {
    Hostonly_ stamped_printf("\tParsing error : fail to get 'Operation' in 'Average_option'\n");
    Exit(0);
  }
  if     ( !strcasecmp(str, "on") )   Mode.Average = ON;
  else if( !strcasecmp(str, "off") )  Mode.Average = OFF;
  else {
    Hostonly_ stamped_printf("\tInvalid keyword is described for 'Operation'\n");
    Exit(0);
  }
  
  // 平均操作開始時間
  if ( Mode.Average == ON ) {
    if ( !elmL1->GetValue("Start_Type", &str) ) {
      Hostonly_ stamped_printf("\tParsing error : fail to get 'Start_Type' in 'Average_option'\n");
      Exit(0);
    }
    else {
      if ( !strcasecmp(str, "step") ) {
        Interval[Interval_Manager::tg_avstart].setMode_Step();
      }
      else if ( !strcasecmp(str, "time") ) {
        Interval[Interval_Manager::tg_avstart].setMode_Time();
      }
      else {
        Hostonly_ stamped_printf("\tParsing error : Invalid keyword for 'Start_Type' in 'Average_option'\n");
        Exit(0);
      }
      
      if ( elmL1->GetValue("Start", &ct) ) {
        Interval[Interval_Manager::tg_avstart].setInterval((double)ct);
      }
      else {
        Hostonly_ stamped_printf("\tParsing error : fail to get 'Start' in 'Average_option'\n");
        Exit(0);
      }
    }    
  }
    
}

/**
 @fn void Control::getXML_LES_option(void)
 @brief LES計算のオプションを取得する
 */
void Control::getXML_LES_option(void)
{
  const CfgElem *elmL1=NULL;
  const char* str=NULL;
  REAL_TYPE ct;
  
  if ( !(elmL1 = getXML_Pointer("LES_Option", "steer")) ) Exit(0);
  
  // 計算オプション
  if ( !elmL1->GetValue("LES_Calculation", &str) ) {
    Hostonly_ stamped_printf("\tParsing error : fail to get 'LES_Calculation' in 'LES_Option'\n");
    Exit(0);
  }
  if     ( !strcasecmp(str, "on") )    LES.Calc = ON;
  else if( !strcasecmp(str, "off") )   LES.Calc = OFF;
  else {
    Hostonly_ stamped_printf("\tParsing error : Invalid keyword for 'LES'\n");
    Exit(0);
  }
  
  if ( LES.Calc == OFF ) return;
  
  // モデル
  if( !elmL1->GetValue("model", &str) ) {
    Hostonly_ stamped_printf("\tParsing error : fail to get 'Model' in 'LES_Option'\n");
    Exit(0);
  }
  if      ( !strcasecmp(str, "smagorinsky") )  LES.Model = Smagorinsky;
  else if ( !strcasecmp(str, "Low_Reynolds") ) LES.Model = Low_Reynolds;
  else if ( !strcasecmp(str, "Dynamic") )      LES.Model = Dynamic;
  
  // Cs係数
  if( !elmL1->GetValue("Cs", &ct) ) {
    Hostonly_ stamped_printf("\tParsing error : fail to get 'Cs' in 'LES_Option'\n");
    Exit(0);
  }
  LES.Cs = ct;
  
  // Cs係数
  if( !elmL1->GetValue("Damping_Factor", &ct) ) {
    Hostonly_ stamped_printf("\tParsing error : fail to get 'Damping_Factor' in 'LES_Option'\n");
    Exit(0);
  }
  LES.damping_factor = ct;

}

/**
 @fn void Control::getXML_Derived(void)
 @brief 派生して計算する変数のオプションを取得する
 */
void Control::getXML_Derived(void)
{
  const CfgElem *elmL1=NULL;
  const char* str=NULL;
  
  if ( !(elmL1 = getXML_Pointer("Derived_Variable", "steer")) ) Exit(0);
  
  // 全圧
  if ( !elmL1->GetValue("Total_Pressure", &str) ) {
    Hostonly_ stamped_printf("\tParsing error : fail to get 'Total_Pressure' in 'Derived_Variable'\n");
    Exit(0);
  }
  if     ( !strcasecmp(str, "on") )  Mode.TP = ON;
  else if( !strcasecmp(str, "off") ) Mode.TP = OFF;
  else {
    Hostonly_ stamped_printf("\tInvalid keyword is described for 'Total_Pressure'\n");
    Exit(0);
  }
  
  // 渦度ベクトル
  if ( !elmL1->GetValue("Vorticity", &str) ) {
    Hostonly_ stamped_printf("\tParsing error : fail to get 'Vorticity' in 'Derived_Variable'\n");
    Exit(0);
  }
  if     ( !strcasecmp(str, "on") )  Mode.VRT = ON;
  else if( !strcasecmp(str, "off") ) Mode.VRT = OFF;
  else {
    Hostonly_ stamped_printf("\tInvalid keyword is described for 'Vorticity'\n");
    Exit(0);
  }
  
  // 速度勾配テンソルの第2不変量
  if ( !elmL1->GetValue("2nd_Invariant_of_VGT", &str) ) {
    Hostonly_ stamped_printf("\tParsing error : fail to get '2nd_Invariant_of_VGT' in 'Derived_Variable'\n");
    Exit(0);
  }
  if     ( !strcasecmp(str, "on") )  Mode.I2VGT = ON;
  else if( !strcasecmp(str, "off") ) Mode.I2VGT = OFF;
  else {
    Hostonly_ stamped_printf("\tInvalid keyword is described for '2nd_Invariant_of_VGT'\n");
    Exit(0);
  }
  
  // ヘリシティ
  if ( !elmL1->GetValue("Helicity", &str) ) {
    Hostonly_ stamped_printf("\tParsing error : fail to get 'Helicity' in 'Derived_Variable'\n");
    Exit(0);
  }
  if     ( !strcasecmp(str, "on") )  Mode.Helicity = ON;
  else if( !strcasecmp(str, "off") ) Mode.Helicity = OFF;
  else {
    Hostonly_ stamped_printf("\tInvalid keyword is described for 'Helicity'\n");
    Exit(0);
  }
}

/**
 @fn void Control::getXML_FileIO(void)
 @brief ファイル入出力に関するパラメータを取得し，sphフォーマットの出力の並列モードを指定する．
 @note インターバルパラメータは，setParameters()で無次元して保持
 */
void Control::getXML_FileIO(void)
{
  SklCfgInFile* infile   = (SklCfgInFile*)CF->GetInFileFirst();
  SklCfgOutFile* outfile = (SklCfgOutFile*)CF->GetOutFileFirst();

  const CfgElem *elmL1=NULL;
  const char* str=NULL;
  REAL_TYPE f_val=0.0;

  if ( !(elmL1 = getXML_Pointer("File_IO", "steer")) ) Exit(0);

  // 出力単位
  if ( !elmL1->GetValue("Unit_of_file", &str) ) {
		Hostonly_ stamped_printf("\tParsing error : Invalid string for 'Unit_of_File' in 'File_IO'\n");
		Exit(0);
	}
  if     ( !strcasecmp(str, "Dimensional") )      Unit.File = DIMENSIONAL;
  else if( !strcasecmp(str, "Non_Dimensional") )  Unit.File = NONDIMENSIONAL;
  else {
    Hostonly_ stamped_printf("\tInvalid keyword is described at 'Unit_of_File' section in 'File_IO'\n");
    Exit(0);
  }

  // 出力ガイドセルモード
  if ( !elmL1->GetValue("Guide_Out", &str) ) {
    Hostonly_ stamped_printf("\tParsing error : fail to get 'Guide_Out' in 'File_IO'\n");
    Exit(0);
  }
  if     ( !strcasecmp(str, "without") )  GuideOut = 0;
  else if( !strcasecmp(str, "with") )     GuideOut = (unsigned)guide;
  else {
    Hostonly_ stamped_printf("\tInvalid keyword is described for 'Guide_Out'\n");
    Exit(0);
  }

  // ファイル出力モード
  if ( !elmL1->GetValue("Output_Mode", &str) ) {
    Hostonly_ stamped_printf("\tParsing error : fail to get 'Output_Mode' in 'File_IO'\n");
    Exit(0);
  }
  if     ( !strcasecmp(str, "Normal") )     FIO.FileOut = IO_normal;
  else if( !strcasecmp(str, "Forced") )     FIO.FileOut = IO_forced;
  else if( !strcasecmp(str, "Every_Time") ) FIO.FileOut = IO_everytime;
  else {
    Hostonly_ stamped_printf("\tInvalid keyword is described for 'Output_Mode'\n");
    Exit(0);
  }

  // デバッグ用のdiv(u)の出力指定
  if ( !elmL1->GetValue("Debug_divergence", &str) ) {
    Hostonly_ stamped_printf("\tParsing error : fail to get 'Debug_Divergence' in 'File_IO'\n");
    Exit(0);
  }
  if     ( !strcasecmp(str, "on") )    FIO.Div_Debug = ON;
  else if( !strcasecmp(str, "off") )   FIO.Div_Debug = OFF;
  else {
    Hostonly_ stamped_printf("\tInvalid keyword is described for 'Debug_Divergence'\n");
    Exit(0);
  }

  // 入力ファイルの並列処理方式を選択
  if ( !elmL1->GetValue("Parallel_Input", &str) ) {
    Hostonly_ stamped_printf("\tParsing error : fail to get 'Parallel_Input' in 'File_IO'\n");
    Exit(0);
  }
  if     ( !strcasecmp(str, "Master") )  FIO.IO_Input = IO_GATHER;
  else if( !strcasecmp(str, "Local") )   FIO.IO_Input = IO_DISTRIBUTE;
  else {
    Hostonly_ stamped_printf("\tInvalid keyword is described for 'Parallel_Input'\n");
    Exit(0);
  }

  // 出力ファイルの並列処理方式を選択
  if ( !elmL1->GetValue("Parallel_Output", &str) ) {
    Hostonly_ stamped_printf("\tParsing error : fail to get 'Parallel_Output' in 'File_IO'\n");
    Exit(0);
  }
  if     ( !strcasecmp(str, "Master") )  FIO.IO_Output = IO_GATHER;
  else if( !strcasecmp(str, "Local") )   FIO.IO_Output = IO_DISTRIBUTE;
  else {
    Hostonly_ stamped_printf("\tInvalid keyword is described for 'Parallel_Output'\n");
    Exit(0);
  }

  // 入力ファイル情報の記述形式をチェックする
  while( infile ){          // check a valid InFile description
    const char *attr, *format, *fname;
    if ( !infile->GetData(&attr, &format, &fname) ) {
      Hostonly_ stamped_printf("\tParsing error : InFile description\n");
      Exit(0);
    }
    
    // V-Sphereフレームワークにパラメータを設定
    if ( !strcasecmp(format, "sph") ) {
      if ( FIO.IO_Input == IO_GATHER ) {
        infile->UnsetMultiInput();
      }
      else {
        infile->SetMultiInput();
      }
    }
    
    infile = (SklCfgInFile*)CF->GetInFileNext(infile);
  }

  // 出力ファイル情報の記述形式をチェックし，並列出力モードをフレームワークに通知する
  while( outfile ){          // check a valid OutFile description
    const char *attr, *format, *fname;
    unsigned interval;
    if ( !(attr = outfile->GetAttr()) ) {
      Hostonly_ stamped_printf("\tParsing error : OutFile description\n");
      Exit(0);
    }
    if ( !(format = outfile->GetFormat()) ) {
      Hostonly_ stamped_printf("\tParsing error : OutFile description\n");
      Exit(0);
    }
    if( !(fname = outfile->GetBaseName()) ) {
      Hostonly_ stamped_printf("\tParsing error : OutFile description\n");
      Exit(0);
    }
    
    // V-Sphereフレームワークにパラメータを設定
    if ( !strcasecmp(format, "sph") ) {
      // 並列出力
      if ( FIO.IO_Output == IO_GATHER ) {
        outfile->UnsetMultiOutput();
      }
      else {
        outfile->SetMultiOutput();
      }
    }
    
    outfile = (SklCfgOutFile*)CF->GetOutFileNext(outfile);
  }

  // インターバル 瞬時値
  if ( !elmL1->GetValue("Instant_Interval_Type", &str) ) {
    Hostonly_ stamped_printf("\tParsing error : fail to get 'Instant_Interval_Type' in 'File_IO'\n");
    Exit(0);
  }
  else {
    if ( !strcasecmp(str, "step") ) {
      Interval[Interval_Manager::tg_instant].setMode_Step();
    }
    else if ( !strcasecmp(str, "time") ) {
      Interval[Interval_Manager::tg_instant].setMode_Time();
    }
    else {
      Hostonly_ stamped_printf("\tParsing error : Invalid keyword for 'Instant_Interval_Type' in 'File_IO'\n");
      Exit(0);
    }
    
    if ( elmL1->GetValue("Instant_Interval", &f_val) ) {
      Interval[Interval_Manager::tg_instant].setInterval((double)f_val);
    }
    else {
      Hostonly_ stamped_printf("\tParsing error : fail to get 'Instant_Interval' in 'File_IO'\n");
      Exit(0);
    }
  }

  // インターバル　平均値
  if ( !elmL1->GetValue("Averaged_Interval_Type", &str) ) {
    Hostonly_ stamped_printf("\tParsing error : fail to get 'Averaged_Interval_Type' in 'File_IO'\n");
    Exit(0);
  }
  else {
    if ( !strcasecmp(str, "step") ) {
      Interval[Interval_Manager::tg_average].setMode_Step();
    }
    else if ( !strcasecmp(str, "time") ) {
      Interval[Interval_Manager::tg_average].setMode_Time();
    }
    else {
      Hostonly_ stamped_printf("\tParsing error : Invalid keyword for 'Averaged_Interval_Type' in 'File_IO'\n");
      Exit(0);
    }
    
    if ( elmL1->GetValue("Averaged_Interval", &f_val) ) {
      Interval[Interval_Manager::tg_average].setInterval((double)f_val);
    }
    else {
      Hostonly_ stamped_printf("\tParsing error : fail to get 'Averaged_Interval' in 'File_IO'\n");
      Exit(0);
    }
  }

}

/**
 @fn void Control::tell_Interval_2_Sphere(void)
 @brief SPHEREフレームワークにファイル出力インターバルを教える
 @note パラメータは，setParameters()で無次元して保持
 */
void Control::tell_Interval_2_Sphere(void)
{
  SklCfgOutFile* outfile = (SklCfgOutFile*)CF->GetOutFileFirst();
  
  while( outfile ){
    const char *attr;
    if ( !(attr = outfile->GetAttr()) ) {
      Hostonly_ stamped_printf("\tParsing error : OutFile description\n");
      Exit(0);
    }
    
    if      ( !strcasecmp(attr, "velocity") )       outfile->SetInterval( Interval[Interval_Manager::tg_instant].getIntervalStep() );
    else if ( !strcasecmp(attr, "pressure") )       outfile->SetInterval( Interval[Interval_Manager::tg_instant].getIntervalStep() );
    else if ( !strcasecmp(attr, "temperature") )    outfile->SetInterval( Interval[Interval_Manager::tg_instant].getIntervalStep() );
    else if ( !strcasecmp(attr, "totalpressure") )  outfile->SetInterval( Interval[Interval_Manager::tg_instant].getIntervalStep() );
    else if ( !strcasecmp(attr, "vorticity") )      outfile->SetInterval( Interval[Interval_Manager::tg_instant].getIntervalStep() );
    else if ( !strcasecmp(attr, "2ndinvrntvgt") )   outfile->SetInterval( Interval[Interval_Manager::tg_instant].getIntervalStep() );
    else if ( !strcasecmp(attr, "helicity") )       outfile->SetInterval( Interval[Interval_Manager::tg_instant].getIntervalStep() );
    else if ( !strcasecmp(attr, "vof") )            outfile->SetInterval( Interval[Interval_Manager::tg_instant].getIntervalStep() );
    else if ( !strcasecmp(attr, "divergence") )     outfile->SetInterval( Interval[Interval_Manager::tg_instant].getIntervalStep() );
    
    outfile = (SklCfgOutFile*)CF->GetOutFileNext(outfile);
  }
}

/**
 @fn void Control::tell_Avr_Interval_2_Sphere(void)
 @brief SPHEREフレームワークに平均値ファイル出力インターバルを教える
 @note パラメータは，setParameters()で無次元して保持
 */
void Control::tell_Avr_Interval_2_Sphere(void)
{
  SklCfgOutFile* outfile = (SklCfgOutFile*)CF->GetOutFileFirst();
  
  while( outfile ){
    const char *attr;
    if ( !(attr = outfile->GetAttr()) ) {
      Hostonly_ stamped_printf("\tParsing error : OutFile description\n");
      Exit(0);
    }
    
    if      ( !strcasecmp(attr, "avrvelocity") )    outfile->SetInterval( Interval[Interval_Manager::tg_average].getIntervalStep() );
    else if ( !strcasecmp(attr, "avrpressure") )    outfile->SetInterval( Interval[Interval_Manager::tg_average].getIntervalStep() );
    else if ( !strcasecmp(attr, "avrtemperature") ) outfile->SetInterval( Interval[Interval_Manager::tg_average].getIntervalStep() );
    
    outfile = (SklCfgOutFile*)CF->GetOutFileNext(outfile);
  }
}

/**
 @fn void Control::getXML_Version(void)
 @brief バージョン情報の取得
 */
void Control::getXML_Version(void)
{
  const CfgElem *elmL1=NULL;
  int ct;
  
  if ( !(elmL1 = getXML_Pointer("Version_Info", "steer")) ) Exit(0);
  
  // FlowBase
  if ( !elmL1->GetValue("Flow_Base", &ct) ) {
    Hostonly_ stamped_printf("\tParsing error : Invalid value for 'Flow_Base' in 'Version_Info'\n");
    Exit(0);
  }
  FB_version = (unsigned)ct;
  
  // FlowBase
  if ( !elmL1->GetValue("CBC", &ct) ) {
    Hostonly_ stamped_printf("\tParsing error : Invalid value for 'CBC' in 'Version_Info'\n");
    Exit(0);
  }
  version = (unsigned)ct;
}

/**
 @fn void Control::getXML_Wall_type(void)
 @brief 壁面上の扱いを指定する
 */
void Control::getXML_Wall_type(void)
{
  const CfgElem *elmL1=NULL;
  const char *str=NULL;
  
  if ( !(elmL1 = getXML_Pointer("Treatment_of_wall", "steer")) ) {
    Hostonly_ stamped_printf("\tNo 'Treatment_of_wall' tag.\n");
    Exit(0);
  }
  
  // 圧力のタイプ
  if ( !elmL1->GetValue(CfgIdt("Pressure_Gradient"), &str) ) {
		Hostonly_ stamped_printf("\tParsing error : Invalid string for 'Pressure_Gradient' in 'Treatment_of_wall'\n");
		Exit(0);
	}
  if     ( !strcasecmp(str, "grad_zero") )   Mode.PrsNeuamnnType = P_GRAD_ZERO;
  else if( !strcasecmp(str, "grad_NS") )     Mode.PrsNeuamnnType = P_GRAD_NS;
  else {
    Hostonly_ stamped_printf("\tInvalid keyword is described for 'Pressure_Gradient'\n");
    Exit(0);
  }
  
  // 壁面摩擦応力の計算モード
  if ( !elmL1->GetValue(CfgIdt("Velocity_Profile"), &str) ) {
		Hostonly_ stamped_printf("\tParsing error : Invalid string for 'Velocity_Profile' in 'Treatment_of_wall'\n");
		Exit(0);
	}
  
  if     ( !strcasecmp(str, "no_slip") )     Mode.Wall_profile = No_Slip;
  else if( !strcasecmp(str, "Slip") )        Mode.Wall_profile = Slip;
  else if( !strcasecmp(str, "Law_of_Wall") ) Mode.Wall_profile = Log_Law;
  else {
    Hostonly_ stamped_printf("\tInvalid keyword is described for 'Velocity_Profile'\n");
    Exit(0);
  }
}

/**
 @fn void Control::getXML_ChangeID(void)
 @brief Cell IDのゼロを指定IDに変更するオプションを取得する（隠しパラメータ）
 @note 'Change_ID'の文字列チェックはしないので注意して使うこと
 */
void Control::getXML_ChangeID(void)
{
  const CfgElem *elemTop=NULL;
  int ct=0;
  
  elemTop = CF->GetTop(STEER);
  if ( !elemTop->GetValue(CfgIdt("Change_ID"), &ct) ) return;
  
  if ( ct < 0 ) {
    Hostonly_ printf("Error : ID should be positive [%d]\n", ct);
    Exit(0);
  }
  else {
    Hide.Change_ID = (unsigned)ct;
  }
}

/**
 @fn void Control::getXML_PMtest(void)
 @brief 性能試験モードを取得する（隠しパラメータ）
 @note 'Performance_Test'の文字列チェックはしないので注意して使うこと
 */
void Control::getXML_PMtest(void)
{
  const CfgElem *elemTop=NULL;
  const char* str=NULL;
  
  elemTop = CF->GetTop(STEER);
  if ( !elemTop->GetValue(CfgIdt("Performance_Test"), &str) ) return;
  
  if     ( !strcasecmp(str, "on") )   Hide.PM_Test = ON;
  else if( !strcasecmp(str, "off") )  Hide.PM_Test = OFF;
  else {
    Hostonly_ stamped_printf("\tInvalid keyword is described for 'Performance_Test'\n");
    Exit(0);
  }
}

/**
 @fn void Control::getXML_Scaling(void)
 @brief スケーリングファクタを取得する（隠しパラメータ）
 @note 'Scaling_factor'の文字列チェックはしないので注意して使うこと
 */
void Control::getXML_Scaling(void)
{
  const CfgElem *elemTop=NULL;
  REAL_TYPE ct=0.0;
  Hide.Scaling_Factor = 1.0;
  
  elemTop = CF->GetTop(STEER);
  if ( !elemTop->GetValue(CfgIdt("Scaling_factor"), &ct) ) return;

  Hide.Scaling_Factor = ( ct <= 0.0 ) ? 0.0 : ct;
  if ( Hide.Scaling_Factor <= 0.0 ) {
    Hostonly_ printf("Error : Scaling factor should be positive [%f]\n", ct);
    Exit(0);
  }
}

/**
 @fn void Control::getXML_InnerItr(void)
 @brief 内部反復回数を取得する
 */
void Control::getXML_InnerItr(void)
{
  ParseSteer Tree(CF);
  if ( !Tree.getParam("Inner_Iteration", InnerItr) ) Exit(0);
}

/**
 @fn REAL_TYPE Control::getCellSize(unsigned* G_size)
 @brief 計算内部領域の全セル数を返す
 @param G_size 計算領域全体の分割数
 */
REAL_TYPE Control::getCellSize(unsigned* G_size)
{
  REAL_TYPE cell_max=0.0;
  
  switch (NoDimension) {
    case 2:
      cell_max = (REAL_TYPE)G_size[0] * (REAL_TYPE)G_size[1];
      break;
      
    case 3:
      cell_max = (REAL_TYPE)G_size[0] * (REAL_TYPE)G_size[1] * (REAL_TYPE)G_size[2];
      break;
  }
  return cell_max;
}

/**
 @fn void Control::printArea(FILE* fp, unsigned G_Fcell, unsigned G_Acell, unsigned* G_size)
 @brief 有効セル数を表示する
 @param fp 出力ファイルポインタ
 @param G_Fcell グローバルなFluid cell
 @param G_Acell グローバルなActive cell
 @param G_size global size
 */
void Control::printArea(FILE* fp, unsigned G_Fcell, unsigned G_Acell, unsigned* G_size)
{
  if( !fp ) {
    stamped_printf("\tFail to write into file\n");
    Exit(0);
  }
  
  REAL_TYPE cell_max = getCellSize(G_size);

  fprintf(fp,"\n---------------------------------------------------------------------------\n\n");
  fprintf(fp,"\n\t>> Effective cells and Open Area of Computational Domain\n\n");

  fprintf(fp,"\tFluid cell  inside whole Computational domain = %15d (%5.1f percent)\n", G_Fcell, (REAL_TYPE)G_Fcell/cell_max *100.0);
  fprintf(fp,"\tActive cell                                   = %15d (%5.1f percent)\n", G_Acell, (REAL_TYPE)G_Acell/cell_max *100.0);
  
  fprintf(fp,"\n\tFace :      Element (Open ratio)\n");
  for (unsigned i=0; i<NoDimension*2; i++) {
    fprintf(fp,"\t  %s : %12.0f (%6.2f percent)\n", getDirection(i).c_str(), OpenDomain[i], OpenDomainRatio(i, OpenDomain[i], NoDimension, G_size));
  }
  fprintf(fp,"\n");
  fflush(fp);
}

/**
 @fn string Control::getDirection(unsigned dir)
 @brief 方向を返す
 @retval 方向の文字
 */
string Control::getDirection(unsigned dir)
{
  string face;

  if      (dir == X_MINUS) face = "X-";
  else if (dir == X_PLUS)  face = "X+";
  else if (dir == Y_MINUS) face = "Y-";
  else if (dir == Y_PLUS)  face = "Y+";
  else if (dir == Z_MINUS) face = "Z-";
  else if (dir == Z_PLUS)  face = "Z+";

  return face;
}

/**
 @fn string Control::getNormString(unsigned d)
 @brief ノルムのタイプを返す
 @retval ノルムの文字
 */
string Control::getNormString(unsigned d)
{
  string nrm;
	
  if      (d == ItrCtl::v_div_max)       nrm = "V - Max. Norm of Divergence";
	else if (d == ItrCtl::v_div_max_dbg)   nrm = "V - Max. Norm of Divergence with Monitoring  ### Forced to be selected since Iteration Log is specified ###";
  else if (d == ItrCtl::v_div_l2)        nrm = "V - L2 Norm of Divergence";
  else if (d == ItrCtl::p_res_l2_a)      nrm = "P - L2 Norm of Absolute Residual";
  else if (d == ItrCtl::p_res_l2_r)      nrm = "P - L2 Norm of Relative Residual";
	else if (d == ItrCtl::v_res_l2_a)      nrm = "V - L2 Norm of Absolute Residual";
  else if (d == ItrCtl::v_res_l2_r)      nrm = "V - L2 Norm of Relative Residual";
  else if (d == ItrCtl::t_res_l2_a)      nrm = "T - L2 Norm of Absolute Residual";
  else if (d == ItrCtl::t_res_l2_r)      nrm = "T - L2 Norm of Relative Residual";
	
  return nrm;
}

/**
 @fn void Control::getXML_KindOfSolver(const CfgElem *elmL1)
 @brief ソルバーの計算対象種別と浮力モードを取得
 @param elmL1  XMLツリーのポインタ
 */
void Control::getXML_KindOfSolver(const CfgElem *elmL1)
{
	const char *str=NULL;
  
  if ( !elmL1->GetValue("Kind_of_solver", &str) ) {
		Hostonly_ stamped_printf("\tParsing error : Invalid char* value for 'Kind_of_solver' in 'Solver_Property'\n");
		Exit(0);
	}

  if     ( !strcasecmp(str, "Flow_Only") )                KindOfSolver = FLOW_ONLY;
  else if( !strcasecmp(str, "Thermal_Flow") )             KindOfSolver = THERMAL_FLOW;
  else if( !strcasecmp(str, "Thermal_Flow_Natural") )     KindOfSolver = THERMAL_FLOW_NATURAL;
  else if( !strcasecmp(str, "Conjugate_Heat_Transfer") )  KindOfSolver = CONJUGATE_HEAT_TRANSFER;
  else if( !strcasecmp(str, "Solid_Conduction") )         KindOfSolver = SOLID_CONDUCTION;
  else {
    Hostonly_ stamped_printf("\tInvalid keyword is described for Kind_of_Solver\n");
    Exit(0);
  }
  
  // Buoyancy option
  if ( (KindOfSolver==THERMAL_FLOW) || (KindOfSolver==THERMAL_FLOW_NATURAL) || (KindOfSolver==CONJUGATE_HEAT_TRANSFER) ) {
    if ( !elmL1->GetValue("Buoyancy", &str) ) {
      Hostonly_ stamped_printf("\tParsing error : Invalid char* value for 'Buoyancy' in 'Solver_Property'\n");
      Exit(0);
    }
    if     ( !strcasecmp(str, "Boussinesq") )   Mode.Buoyancy = BOUSSINESQ;
    else if( !strcasecmp(str, "Low_Mach") )     Mode.Buoyancy = LOW_MACH;
    else if( !strcasecmp(str, "No_Buoyancy") )  Mode.Buoyancy = NO_BUOYANCY;
    else {
      Hostonly_ stamped_printf("\tInvalid keyword is described for 'Buoyancy'\n");
      Exit(0);
    }
  }
}

/**
 @fn void Control::getXML_Unit(void)
 @brief 入力ファイルに記述するパラメータとファイルの有次元・無次元の指定を取得する
 */
void Control::getXML_Unit(void)
{
  const CfgElem *elmL1=NULL;
	const char *str=NULL;
  
  if ( !(elmL1 = getXML_Pointer("Unit", "steer")) ) Exit(0);
  
  if ( !elmL1->GetValue("Unit_of_input_parameter", &str) ) {
		Hostonly_ stamped_printf("\tParsing error : Invalid string for 'Unit_of_Input_Parameter' in 'Unit'\n");
		Exit(0);
	}
  if     ( !strcasecmp(str, "Dimensional") )      Unit.Param = DIMENSIONAL;
  else if( !strcasecmp(str, "Non_Dimensional") )  Unit.Param = NONDIMENSIONAL;
  else {
    Hostonly_ stamped_printf("\tInvalid keyword is described at 'Unit_of_Input_Parameter' section\n");
    Exit(0);
  }
    
  if ( !elmL1->GetValue("pressure", &str) ) {
    Hostonly_ stamped_printf("\tParsing error : Invalid string for 'Pressure' in 'Unit'\n");
    Exit(0);
  }
  if     ( !strcasecmp(str, "Gauge")  )   Unit.Prs = CompoList::Gauge;
  else if( !strcasecmp(str, "Absolute") ) Unit.Prs = CompoList::Absolute;
  else {
    Hostonly_ stamped_printf("\tInvalid keyword is described at 'Pressure' in 'Unit'\n");
    Exit(0);
  }
  
  if ( isHeatProblem() ) {
    if ( !elmL1->GetValue("temperature", &str) ) {
      Hostonly_ stamped_printf("\tParsing error : Invalid string for 'Temperature' in 'Unit'\n");
      Exit(0);
    }
    if     ( !strcasecmp(str, "Celsius") )  Unit.Temp = CompoList::Unit_CELSIUS;
    else if( !strcasecmp(str, "Kelvin") )   Unit.Temp = CompoList::Unit_KELVIN;
    else {
      Hostonly_ stamped_printf("\tInvalid keyword is described at 'Temperature' in 'Unit'\n");
      Exit(0);
    }
  }
  
}

/**
 @fn void Control::getXML_CheckParameter(void)
 @brief パラメータ入力チェックモードの取得
 */
void Control::getXML_CheckParameter(void)
{
  const char *keyword=NULL;
  ParseSteer Tree(CF);
  
  if ( !(keyword=Tree.getParam("Check_Parameter")) ) Exit(0);
  
  if     ( !strcasecmp(keyword, "On") )   CheckParam = ON;
  else if( !strcasecmp(keyword, "Off") )  CheckParam = OFF;
  else {
    Hostonly_ stamped_printf("\tInvalid keyword is described for 'Check_Parameter'\n");
    Exit(0);
  }
}

/**
 @fn void Control::getXML_Para_Wind(void)
 @brief 入射流速角度を取得
 */
void Control::getXML_Para_Wind(void)
{
  ParsePara Tree(CF);
  if ( !Tree.IsSetElem("Wind_Direction") )      Exit(0);
  if ( !Tree.getEParam("yaw_angle",     yaw))   Exit(0);
  if ( !Tree.getEParam("pitch_angle", pitch))   Exit(0);
  if ( !Tree.getEParam("roll_angle", roll))     Exit(0);
}

/**
 @fn void Control::getXML_Para_Init(void)
 @brief 初期値の値を取得する
 @note
    - このメソッド内では，初期値は無次元/有次元の判定はしない
    - 無次元指定時の値の変換は　setParameters(MaterialList* mat, CompoList* cmp)
 @see 
    - bool Control::setParameters(MaterialList* mat, CompoList* cmp)
 */
void Control::getXML_Para_Init(void)
{	
	const CfgElem *elmL1=NULL, *elmL2=NULL;

  if ( !(elmL1=getXML_Pointer("Initial_State", "parameter")) ) Exit(0);
	
  // Density
	if( !elmL1->GetValue("Density", &iv.Density) ) {
    Hostonly_ stamped_printf("\tParsing error : Invalid float value for 'Density'\n");
    Exit(0);
  }
	
  // Pressure
	if( !elmL1->GetValue("Pressure", &iv.Pressure) ) {
    Hostonly_ stamped_printf("\tParsing error : Invalid float value for 'Pressure'\n");
    Exit(0);
  }
	
  // Velocity
	if ( !(elmL2  = elmL1->GetElemFirst("Velocity")) ) {
		Hostonly_ stamped_printf("\tParsing error : No 'Velocity' section in 'Initial_State'\n");
		Exit(0);
	}
	if (elmL2->GetParamSize() != 3) {    // check number of Param
		Hostonly_ stamped_printf("\tParsing error : 3 Params should be found in Initial_State Velocity\n");
		Exit(0);
  }
	REAL_TYPE v[3];
	for (unsigned n=0; n<3; n++) v[n]=0.0;
	if ( !(elmL2->GetVctValue("u", "v", "w", &v[0], &v[1], &v[2])) ) {
    Hostonly_ stamped_printf("\tParsing error : fail to get velocity params in 'Initial_State'\n");
    Exit(0);
  }
  iv.VecU = v[0];
  iv.VecV = v[1];
  iv.VecW = v[2];

  // Temperature
  if ( isHeatProblem() ) {
    if( !elmL1->GetValue("Temperature", &iv.Temperature) ) {
			Hostonly_ stamped_printf("\tParsing error : Invalid float value for 'Temperature'\n");
			Exit(0);
		}
  }
}

/**
 @fn void Control::getXML_TimeMarching(void)
 @brief 時間積分精度
 @note
    - 未使用
 */
void Control::getXML_TimeMarching(void)
{
  const char *keyword=NULL;
  ParseSteer Tree(CF);

  if ( !(keyword=Tree.getParam("Time_Integration")) ) Exit(0);

  if     ( !strcasecmp(keyword, "1st_order") ) MarchingScheme = TM_O_1ST;
  else if( !strcasecmp(keyword, "2nd_order") ) MarchingScheme = TM_O_2ND;
  else if( !strcasecmp(keyword, "3rd_order") ) MarchingScheme = TM_O_3RD;
  else {
    Hostonly_ stamped_printf("\tInvalid keyword is described for 'Time_Integration'\n");
    Exit(0);
  }
}

/**
 @fn bool Control::chkMediumConsistency(void)
 @brief 全Voxelモデルの媒質数とKOSの整合性をチェック
 @retval エラーコード
 @note
    - NoMediumSolidなどは，ParseBC::setMedium()で取得
 */
bool Control::chkMediumConsistency(void)
{
  SklParaManager* para_mng = ParaCmpo->GetParaManager();
  unsigned long nmSolid = NoMediumSolid;
  unsigned long nmFluid = NoMediumFluid;

  if( para_mng->IsParallel() ){
    unsigned long nms = nmSolid;
    unsigned long nmf = nmFluid;
    if( !para_mng->Allreduce(&nms, &nmSolid, 1, SKL_ARRAY_DTYPE_ULONG, SKL_SUM, pn.procGrp) ) return false;
    if( !para_mng->Allreduce(&nmf, &nmFluid, 1, SKL_ARRAY_DTYPE_ULONG, SKL_SUM, pn.procGrp) ) return false;
  }

  if ( (nmFluid == 0) && (nmSolid == 0) ) {
    Hostonly_ printf("\tError : No medium\n");
    return false;
  }

  switch (KindOfSolver) {
    case FLOW_ONLY:
    case THERMAL_FLOW:
    case THERMAL_FLOW_NATURAL:

      if ( nmFluid == 0 ) {
        Hostonly_ printf("\tError : No FLUID medium\n");
        return false;
      }
      break;

    case CONJUGATE_HEAT_TRANSFER:
      if ( ( nmFluid == 0 ) || ( nmSolid == 0 ) ) {
        Hostonly_ printf("\tError : Fluid/Solid should have at least one medium.\n");
        return false;
      }
      break;

    case SOLID_CONDUCTION:
      if ( nmSolid == 0 ) {
        Hostonly_ printf("\tError : No Solid medium\n");
        return false;
      }
      break;
  };

  return true;
}

//@fn void Control::getXML_Log(void)
//@brief ログ出力モードを取得
//@note インターバルパラメータは，setParameters()で無次元して保持
void Control::getXML_Log(void)
{
  const CfgElem *elmL1=NULL;
	const char *str=NULL;
  REAL_TYPE f_val=0.0;
  
  if ( !(elmL1 = getXML_Pointer("Log", "steer")) ) Exit(0);
  
  // 出力単位
  if ( !elmL1->GetValue("Unit_of_Log", &str) ) {
		Hostonly_ stamped_printf("\tParsing error : Invalid string for 'Unit_of_Log' in 'Unit'\n");
		Exit(0);
	}
  if     ( !strcasecmp(str, "Dimensional") )      Unit.Log = DIMENSIONAL;
  else if( !strcasecmp(str, "Non_Dimensional") )  Unit.Log = NONDIMENSIONAL;
  else {
    Hostonly_ stamped_printf("\tInvalid keyword is described at 'Unit_of_Log' section\n");
    Exit(0);
  }
  
  // Log_Base
  if ( !elmL1->GetValue(CfgIdt("Log_Base"), &str) ) {
		Hostonly_ stamped_printf("\tParsing error : Invalid string for 'Log_Base' in 'Log'\n");
		Exit(0);
	}
  if     ( !strcasecmp(str, "on") )   Mode.Log_Base = ON;
  else if( !strcasecmp(str, "off") )  Mode.Log_Base = OFF;
  else {
    Hostonly_ stamped_printf("\tInvalid keyword is described for 'Log_File'\n");
    Exit(0);
  }
  
  // Log_Iteration
  if ( !elmL1->GetValue(CfgIdt("Log_Iteration"), &str) ) {
		Hostonly_ stamped_printf("\tParsing error : Invalid string for 'Log_Iteration' in 'Log'\n");
		Exit(0);
	}
  if     ( !strcasecmp(str, "on") )   Mode.Log_Itr = ON;
  else if( !strcasecmp(str, "off") )  Mode.Log_Itr = OFF;
  else {
    Hostonly_ stamped_printf("\tInvalid keyword is described for 'Log_Iteration'\n");
    Exit(0);
  }
  
  // Log_Wall_Info
  if ( Mode.Wall_profile == Log_Law ) {
    if ( !elmL1->GetValue(CfgIdt("Log_Wall_Info"), &str) ) {
      Hostonly_ stamped_printf("\tParsing error : Invalid string for 'Log_Wall_Info' in 'Log'\n");
      Exit(0);
    }
    if     ( !strcasecmp(str, "on") )   Mode.Log_Wall = ON;
    else if( !strcasecmp(str, "off") )  Mode.Log_Wall = OFF;
    else {
      Hostonly_ stamped_printf("\tInvalid keyword is described for 'Log_Wall_Info'\n");
      Exit(0);
    }
  }
  
  // Log_Profiling
  if ( !elmL1->GetValue(CfgIdt("Log_Profiling"), &str) ) {
		Hostonly_ stamped_printf("\tParsing error : Invalid string for 'Log_Profiling' in 'Log'\n");
		Exit(0);
	}
  if     ( !strcasecmp(str, "on") )     Mode.Profiling = ON;
  else if( !strcasecmp(str, "off") )    Mode.Profiling = OFF;
  else if( !strcasecmp(str, "detail") ) Mode.Profiling = DETAIL;
  else {
    Hostonly_ stamped_printf("\tInvalid keyword is described for 'Log_Profiling'\n");
    Exit(0);
  }
  
  // Interval console
  if ( !elmL1->GetValue("Console_Interval_Type", &str) ) {
    Hostonly_ stamped_printf("\tParsing error : fail to get 'Console_Interval_Type' in 'Log'\n");
    Exit(0);
  }
  else {
    if ( !strcasecmp(str, "step") ) {
      Interval[Interval_Manager::tg_console].setMode_Step();
    }
    else if ( !strcasecmp(str, "time") ) {
      Interval[Interval_Manager::tg_console].setMode_Time();
    }
    else {
      Hostonly_ stamped_printf("\tParsing error : Invalid keyword for 'Console_Interval_Type' in 'Log'\n");
      Exit(0);
    }
    
    if ( elmL1->GetValue("Console_Interval", &f_val) ) {
      Interval[Interval_Manager::tg_console].setInterval((double)f_val);
    }
    else {
      Hostonly_ stamped_printf("\tParsing error : fail to get 'Console_Interval' in 'Log'\n");
      Exit(0);
    }
  }
  
  // Interval file_history
  if ( !elmL1->GetValue("History_Interval_Type", &str) ) {
    Hostonly_ stamped_printf("\tParsing error : fail to get 'History_Interval_Type' in 'Log'\n");
    Exit(0);
  }
  else {
    if ( !strcasecmp(str, "step") ) {
      Interval[Interval_Manager::tg_history].setMode_Step();
    }
    else if ( !strcasecmp(str, "time") ) {
      Interval[Interval_Manager::tg_history].setMode_Time();
    }
    else {
      Hostonly_ stamped_printf("\tParsing error : Invalid keyword for 'History_Interval_Type' in 'Log'\n");
      Exit(0);
    }
    
    if ( elmL1->GetValue("History_Interval", &f_val) ) {
      Interval[Interval_Manager::tg_history].setInterval((double)f_val);
    }
    else {
      Hostonly_ stamped_printf("\tParsing error : fail to get 'History_Interval' in 'Log'\n");
      Exit(0);
    }
  }
  
}

/**
 @fn void Control::getXML_VarRange(void)
 @brief 変数の範囲制限モードを取得
 @note 隠しパラメータ
 */
void Control::getXML_VarRange(void)
{
  const CfgElem *elemTop=NULL;
  const char* str=NULL;

  elemTop = CF->GetTop(STEER);
  if ( !elemTop->GetValue(CfgIdt("Scaling_factor"), &str) ) return;

  if     ( !strcasecmp(str, "normal") )   Hide.Range_Limit = Range_Normal;
  else if( !strcasecmp(str, "cutoff") )   Hide.Range_Limit = Range_Cutoff;
  else {
    Hostonly_ stamped_printf("\tInvalid keyword is described for 'Variable Range'\n");
    Exit(0);
  }
}

/**
 @fn void Control::printParaConditions(FILE* fp)
 @brief 計算パラメータの表示
 @param fp
 */
void Control::printParaConditions(FILE* fp)
{
  if ( !fp ) {
    stamped_printf("\tFail to write into file\n");
    Exit(0);
  }

  fprintf(fp,"\n---------------------------------------------------------------------------\n\n");
  fprintf(fp,"\n\t>> Simulation Parameters\n\n");

  fprintf(fp,"\tReference ID              [-]         :  %d\n", RefID);
  fprintf(fp,"\n");

  // Reference values
  fprintf(fp,"\tRef. Length               [m]         : %12.5e\n", RefLength);
  fprintf(fp,"\tRef. Velocity             [m/s]       : %12.5e\n", RefVelocity);
	fprintf(fp,"\tBase Pressure             [Pa]        : %12.5e\n", BasePrs);
  fprintf(fp,"\tRef. Density              [kg/m^3]    : %12.5e\n", RefDensity);
  fprintf(fp,"\tRef. Viscosity            [Pa s]      : %12.5e\n", RefViscosity);
  fprintf(fp,"\tRef. Knmtc Viscosity      [m^2/s]     : %12.5e\n", RefKviscosity);
  fprintf(fp,"\tRef. Specific Heat        [J/(kg K)]  : %12.5e\n", RefSpecificHeat);
  fprintf(fp,"\tRef. Thermal Conductivity [W/(m K)]   : %12.5e\n", RefLambda);
  fprintf(fp,"\tRef. Sound Speed          [m/s]       : %12.5e\n", RefSoundSpeed);
  fprintf(fp,"\tGravity                   [m/s^2]     : %12.5e\n", Gravity);
  fprintf(fp,"\n");
  fprintf(fp,"\tSpacing                   [m] / [-]   : %12.5e / %12.5e\n", dh*RefLength, dh);
  fprintf(fp,"\tTime Scale                [sec]       : %12.5e\n", Tscale);
  fprintf(fp,"\n");
  if ( isHeatProblem() ) {
    fprintf(fp,"\tBase Temperature          [%s] / [-]   : %12.5e / %3.1f\n", (Unit.Temp==CompoList::Unit_KELVIN) ? "K" : "C", FBUtility::convK2Temp(BaseTemp, Unit.Temp), 0.0);
    fprintf(fp,"\tTemperature Diff.         [%s] / [-]   : %12.5e / %3.1f\n", (Unit.Temp==CompoList::Unit_KELVIN) ? "K" : "C", DiffTemp, 1.0);
  }
  fprintf(fp,"\n");

  //REAL_TYPE ay, ap;
  //ay = yaw/180.0*2.0*asin(1.0);
  //ap = pitch/180.0*2.0*asin(1.0);
  //fprintf(fp,"\tYaw   Angle               [rad]/[deg] : %12.5e / %12.5e\n", ay, yaw);
  //fprintf(fp,"\tPitch Angle               [rad]/[deg] : %12.5e / %12.5e\n", ap, pitch);

  fprintf(fp,"\n");
  fprintf(fp,"\tPrandtl  number           [-]         : %12.5e\n", Prandtl);
  if (Mode.PDE == PDE_NS) {
    fprintf(fp,"\tReynolds number           [-]         : %12.5e\n", Reynolds);
  }
  else {
    fprintf(fp,"\tMach number               [-]         : %12.5e\n", Mach);
  }
  if ( isHeatProblem() )  {
    fprintf(fp,"\tPeclet   number           [-]         : %12.5e\n", Peclet);
    fprintf(fp,"\tGrashof  number           [-]         : %12.5e\n", Grashof);
    if (KindOfSolver==THERMAL_FLOW_NATURAL)  fprintf(fp,"\tRayleigh number           [-]         : %12.5e\n", Rayleigh);
  }
  fprintf(fp,"\n");

  fflush(fp);
}

/**
 @fn void Control::printSteerConditions(FILE* fp, ItrCtl* IC, DTcntl* DT, ReferenceFrame* RF)
 @brief 制御パラメータSTEERの表示
 */
void Control::printSteerConditions(FILE* fp, ItrCtl* IC, DTcntl* DT, ReferenceFrame* RF)
{
  if( !fp ) {
    stamped_printf("\tFail to write into file\n");
    Exit(0);
  }

  REAL_TYPE dt = (REAL_TYPE)DT->get_DT();
  bool  err=true;
  double itm=0.0;

  fprintf(fp,"\n---------------------------------------------------------------------------\n\n");
  fprintf(fp,"\n\t>> Solver Control Parameters\n\n");

  fprintf(fp,"\tSolver Properties\n");
  // Basic Equation and PDE
	switch (KindOfSolver) {
    case FLOW_ONLY:
    case THERMAL_FLOW:
    case THERMAL_FLOW_NATURAL:
    case CONJUGATE_HEAT_TRANSFER:
			switch (BasicEqs) {
				case INCMP:
					fprintf(fp,"\t     Basic Equation           :   Incompressible Flow ");
					break;
				case LTDCMP:
					fprintf(fp,"\t     Basic Equation           :   Limited Compressible Flow ");
					break;
				case CMPRSS:
					fprintf(fp,"\t     Basic Equation           :   Compressible Flow ");
					break;
        case INCMP_2PHASE:
					fprintf(fp,"\t     Basic Equation           :   Incompressible Two-Phase Flow ");
					break;
				default:
					stamped_printf("\nError: Basic Equation section\n");
					err=false;
			}
			switch (Mode.PDE) {
				case PDE_EULER:
					fprintf(fp,"PDE_EULER Equations\n");
					if  ( !((KindOfSolver == FLOW_ONLY) || (KindOfSolver == THERMAL_FLOW) || (KindOfSolver == THERMAL_FLOW_NATURAL)) ) {
						fprintf(fp,"\tInvalid combination with Conjugate Heat Transfer nor Solid Conduction\n");
						err=false;
					}
					break;
				case PDE_NS:
					fprintf(fp,"Navier-Stokes Equations\n");
					break;
				default:
					stamped_printf("Error: PDE section\n");
					err=false;
			}
			break;
			
		case SOLID_CONDUCTION:
			fprintf(fp,"\t     Basic Equation           :   Heat Conduction Equation\n");
			break;
  }

  // Steady
  switch (Mode.Steady) {
    case TV_Steady:
      fprintf(fp,"\t     Time Variation           :   Steady\n");
      break;
      
    case TV_Unsteady:
      fprintf(fp,"\t     Time Variation           :   Unsteady\n");
      break;
      
    default:
      stamped_printf("Error: Time Variation[%d]\n", Mode.Steady);
      err=false;
  }
  
  /* Grid System
	if ( KindOfSolver == SOLID_CONDUCTION ) {
		fprintf(fp,"\tVariable Arrangement   :   Cell Centered\n");
	}
	else {
		switch (Mode.VarArrange) {
			case STAGGERED:
				fprintf(fp,"\tVariable Arrangement   :   Staggered\n");
				break;
			case CELL_CENTER:
				fprintf(fp,"\tVariable Arrangement   :   Collocated\n");
				break;
			case NODE:
				fprintf(fp,"\tVariable Arrangement   :   Node\n");
				break;
			default:
				stamped_printf("Error: Grid system section\n");
				err=false;
		}
  }*/
  
  // Shape approximation
  switch (Mode.ShapeAprx) {
    case BINARY:
      fprintf(fp,"\t     Shape Approximation      :   Binary\n");
      break;
      
    case CUT_INFO:
      fprintf(fp,"\t     Shape Approximation      :   Cut-Distance\n");
      break;
      
    default:
      stamped_printf("Error: Shape Approximation section\n");
      err=false;
  }

  // Precision
  if ( Mode.Precision == SPH_SINGLE )
    fprintf(fp,"\t     Precision                :   Single Precision \n");
  else
    fprintf(fp,"\t     Precision                :   Double Precision \n");
  
  // Kind Of Solver
  if (KindOfSolver==FLOW_ONLY) {
    fprintf(fp,"\t     Kind of Solver           :   Flow Only (Non Heat)\n");
  }
  else if ( (KindOfSolver==THERMAL_FLOW) && (Mode.Buoyancy==NO_BUOYANCY) ) {
    fprintf(fp,"\t     Kind of Solver           :   Forced convection without buoyancy\n");
  }
  else if ( (KindOfSolver==THERMAL_FLOW) && (Mode.Buoyancy==BOUSSINESQ) ) {
    fprintf(fp,"\t     Kind of Solver           :   Forced convection with buoyancy : Boussinesq Approximation\n");
  }
  else if ( (KindOfSolver==THERMAL_FLOW) && (Mode.Buoyancy==LOW_MACH) ) {
    fprintf(fp,"\t     Kind of Solver           :   Forced convection with buoyancy : Low Mach Approximation\n");
  }
  else if ( (KindOfSolver==THERMAL_FLOW_NATURAL) && (Mode.Buoyancy==BOUSSINESQ) ) {
    fprintf(fp,"\t     Kind of Solver           :   Natural convection : Boussinesq Approximation\n");
  }
  else if ( (KindOfSolver==THERMAL_FLOW_NATURAL) && (Mode.Buoyancy==LOW_MACH) ) {
    fprintf(fp,"\t     Kind of Solver           :   Natural convection : Low Mach Approximation\n");
  }
  else if (KindOfSolver==SOLID_CONDUCTION) {
    fprintf(fp,"\t     Kind of Solver           :   Solid Conduction\n");
  }
  else {
    fprintf(fp,"\t     Heat Solver type         :   Error\n");
    err=false;
  }
  
  // Flow Algorithm
  switch (AlgorithmF) {
    case Flow_FS_EE_EE:
      fprintf(fp,"\t     Flow Algorithm           :   Fractional Step\n");
      fprintf(fp,"\t        Time marching scheme  :   Euler Explicit O(dt1)\n");
      break;
      
    case Flow_FS_RK_CN:
      fprintf(fp,"\t     Flow Algorithm           :   Fractional Step\n");
      fprintf(fp,"\t        Time marching scheme  :   Runge-kutta and Crank-Nicholson O(dt2)\n");
      break;
      
    case Flow_FS_AB2:
      fprintf(fp,"\t     Flow Algorithm           :   Fractional Step\n");
      fprintf(fp,"\t        Time marching scheme  :   Adams-Bashforth Explicit O(dt2)\n");
      break;
      
    case Flow_FS_AB_CN:
      fprintf(fp,"\t     Flow Algorithm           :   Fractional Step\n");
      fprintf(fp,"\t        Time marching scheme  :   Adams-Bashforth Explicit O(dt2) and Crank-Nicholson O(dt2)\n");
      break;
      
    default:
      stamped_printf("No algorithm is specified for Flow\n"); // this is not error
  }
  
  // Heat Algorithm
  switch (AlgorithmH) {
    case Heat_EE_EE:
      fprintf(fp,"\t     Heat Algorithm           :   Fractional Step\n");
      fprintf(fp,"\t        Time marching scheme  :   Euler Explicit O(dt1)\n");
      break;
      
    case Heat_EE_EI:
      fprintf(fp,"\t     Heat Algorithm           :   Fractional Step\n");
      fprintf(fp,"\t        Time marching scheme  :   Euler Implicit O(dt1)\n");
      break;
      
    default:
      fprintf(fp,"\t     Heat Algorithm           :   \n");
      fprintf(fp,"\t        Time marching scheme  :   \n");
  }
  
  // Convection scheme
	if ( KindOfSolver != SOLID_CONDUCTION ) {
		switch (CnvScheme) {
			case O1_upwind:
				fprintf(fp,"\t     Convective flux scheme   :   Upwind O(dx1)\n");
        break;
			case O3_muscl:
				fprintf(fp,"\t     Convective flux scheme   :   MUSCL O(dx3)\n");
        
				switch (Limiter) {
					case No_Limiter:
						fprintf(fp,"\t         Limiter Function     :   NO\n");
						break;
            
					case MinMod:
						fprintf(fp,"\t         Limiter Function     :   Minmod\n");
						break;
            
					default:
						stamped_printf("Error: Limiter Function section\n");
						err=false;
				}
				break;
        
			default:
				stamped_printf("Error: Convection scheme section\n");
				err=false;
		}
	}
  
  // Reference Frame
  switch (RF->getFrame()) {
    case ReferenceFrame::frm_static:
      fprintf(fp,"\t     Reference Frame          :   Stationary\n");
      break;
      
    case ReferenceFrame::frm_translation:
      if (Unit.Param==DIMENSIONAL) {
        fprintf(fp,"\t     Reference Frame          :   Translational (%12.4e, %12.4e, %12.4e) [m/s]\n", GridVel[0]*RefVelocity, GridVel[1]*RefVelocity, GridVel[2]*RefVelocity);
      }
      else {
        fprintf(fp,"\t     Reference Frame          :   Translational (%12.4e, %12.4e, %12.4e) [-]\n", GridVel[0], GridVel[1], GridVel[2]);
      }
      break;
      
    case ReferenceFrame::frm_rotation:
      fprintf(fp,"\t     Reference Frame          :   Rotational\n");
      break;
      
    default:
      stamped_printf("Error: Reference frame section\n");
      err=false;
  }
  
  
  // 単位系
  fprintf(fp,"\n\tUnit\n");
  fprintf(fp,"\t     Unit of Input Parameter  :   %s\n", (Unit.Param == DIMENSIONAL) ? "Dimensional" : "Non-Dimensional");
  fprintf(fp,"\t             Pressure         :   %s\n", (Unit.Prs == CompoList::Absolute) ? "Absolute Pressure" : "Gauge Pressure");
  fprintf(fp,"\t             Temperature      :   %s\n", (Unit.Temp == CompoList::Unit_KELVIN) ? "Kelvin" : "Celsius");
  
  // 時間制御
  fprintf(fp,"\n\tTime Control\n");
  // Start
  switch (Start) {
    case initial_start:
      fprintf(fp,"\t     Start Condition          :   Impulsive start\n");
      break;
    case re_start:
      fprintf(fp,"\t     Start Condition          :   Restart from file\n");
      break;
    default:
      stamped_printf("Error: start condition section\n");
      err=false;
  }
  
  // 加速時間
  if ( !Interval[Interval_Manager::tg_accelra].isStep() ) {
    itm = Interval[Interval_Manager::tg_accelra].getIntervalTime();
    fprintf(fp,"\t     Acceleration Time        :   %12.5e [sec] / %12.5e [-]\n", itm*Tscale, itm);
  }
  else {
    fprintf(fp,"\t     Acceleration Step        :   %12d\n", Interval[Interval_Manager::tg_accelra].getIntervalStep());
  }
  
  // 時間平均
  if ( Mode.Average == ON ) {
    if ( !Interval[Interval_Manager::tg_avstart].isStep() ) {
      itm = Interval[Interval_Manager::tg_avstart].getIntervalTime();
      fprintf(fp,"\t     Averaging Operation      :   %12.5e [sec] / %12.5e [-]\n", itm*Tscale, itm);
    }
    else {
      fprintf(fp,"\t     Averaging Operation      :   %12d\n", Interval[Interval_Manager::tg_avstart].getIntervalStep());
    }
  }
  else {
    fprintf(fp,"\t     Averaging Operation      :   OFF\n");
  }
  
  // Time Increment
  REAL_TYPE d_R = dh*dh*Reynolds/6.0; // 拡散数
  REAL_TYPE d_P = dh*dh*Peclet/6.0;   // 拡散数
  REAL_TYPE cfl = (REAL_TYPE)DT->get_CFL();
  switch ( DT->get_Scheme() ) {
    case DTcntl::dt_direct:
      fprintf(fp,"\t     Time Increment dt        :   %12.5e [sec] / %12.5e [-] : Direct ", dt*Tscale, dt);
      if ( isHeatProblem() ) {
        fprintf(fp,": Diff. Num. = %7.2e\n", dt/(dh*dh*Peclet));
      }
      else {
        fprintf(fp,"\n");
      }

      break;
      
    case DTcntl::dt_cfl_ref_v:
      fprintf(fp,"\t     Time Increment dt        :   %12.5e [sec] / %12.5e [-] : CFL [%8.5f] with Reference velocity\n", dt*Tscale, dt, cfl);
      break;
      
    case DTcntl::dt_cfl_max_v:
      fprintf(fp,"\t     Time Increment dt        :   %12.5e [sec] / %12.5e [-] : CFL [%8.5f] with Maximum velocity (in case of v=1.0 for Ref.)\n", dt*Tscale, dt, cfl);
      break;
      
    case DTcntl::dt_dfn:
      fprintf(fp,"\t     Time Increment dt        :   %12.5e [sec] / %12.5e [-] : dt restricted by Diffusion number[%8.5f] (Peclet)\n", dt*Tscale, dt, d_P);
      break;
      
    case DTcntl::dt_cfl_dfn_ref_v:
    {
      fprintf(fp,"\t     Time Increment dt        :   %12.5e [sec] / %12.5e [-] : CFL & Diffusion number with Reference velocity\n",  dt*Tscale, dt);
      fprintf(fp,"\t                              :               CFL number                                     : %8.5f [-]\n", cfl);
      fprintf(fp,"\t                              :               dt restricted by Diffusion number (Reynolds)   : %8.5f [-]\n", d_R);
      if ( isHeatProblem() ) {
        fprintf(fp,"\t                              :               dt restricted by Diffusion number (Peclet)     : %8.5f [-]\n", d_P);
      }
      break;
    }
    case DTcntl::dt_cfl_dfn_max_v:
    {
      REAL_TYPE a, b, c;
      a = (REAL_TYPE)DT->dtCFL(1.0);
      b = (REAL_TYPE)DT->dtDFN((double)Reynolds);
      fprintf(fp,"\t     Time Increment dt        :   %12.5e [sec] / %12.5e [-] : CFL & Diffusion number with Maximum velocity (in case of v=1.0 for Ref.)\n", dt*Tscale, dt);
      fprintf(fp,"\t                              :             CFL number                    : %8.5f [-]\n", cfl);
      fprintf(fp,"\t                              :             dt restricted by Diffusion number (Reynolds)   : %8.5f [-]\n", d_R);
      if ( isHeatProblem() ) {
        c = (REAL_TYPE)DT->dtDFN((double)Peclet);
        fprintf(fp,"\t                              :             dt restricted by Diffusion number (Peclet)     : %8.5f [-]\n", d_P);
      }
    }
      break;
      
    case DTcntl::dt_cfl_max_v_cp:
      fprintf(fp,"\t     Time Increment dt        :   %12.5e [sec] / %12.5e [-] : CFL (/w Sound Speed) & Diffusion number with Maximum velocity\n", dt*Tscale, dt);
      break;
      
    default:
      stamped_printf("Error: Time Increment section\n");
      err=false;
  }
  
  // Calculation time/step
  if ( !Interval[Interval_Manager::tg_compute].isStep() ) {
    itm = Interval[Interval_Manager::tg_compute].getIntervalTime();
    fprintf(fp,"\t     Calculation Time         :   %12.5e [sec] / %12.5e [-]\n", itm*Tscale, itm);
  }
  else {
    fprintf(fp,"\t     Calculation Step         :   %12d\n", Interval[Interval_Manager::tg_compute].getIntervalStep());
  }
  
  
  fprintf(fp,"\n\tParallel Mode & File IO\n");

  
  // parallel mode
  SklParaManager* para_mng = ParaCmpo->GetParaManager();
  if (para_mng->GetNodeNum(pn.procGrp) == 1) {
    fprintf(fp,"\t     Parallel Mode            :   Serial\n");
  }
  else if (para_mng->IsMb()) {
    fprintf(fp,"\t     Parallel Mode            :   Multiply-Connected Partitioning\n");
  }
  else if (para_mng->IsEv()) {
    fprintf(fp,"\t     Parallel Mode            :   Equal Partitioning\n");
  }
  else {
    Exit(0);
  }
  
  fprintf(fp,"\t     Unit of File             :   %s\n", (Unit.File == DIMENSIONAL) ? "Dimensional" : "Non-Dimensional");
  
  // InputMode
  fprintf(fp,"\t     Input Mode               :   %s\n", (FIO.IO_Input==IO_GATHER) ? "Master IO" : "Local IO");
  
  // OutputMode
  switch (FIO.FileOut) {
    case IO_normal:
      fprintf(fp,"\t     Output Mode              :   %s   Normal\n", (FIO.IO_Output==IO_GATHER) ? "Master IO" : "Local IO");
      break;
      
    case IO_forced:
      fprintf(fp,"\t     Output Mode              :   %s   Forced\n", (FIO.IO_Output==IO_GATHER) ? "Master IO" : "Local IO");
      break;
      
    case IO_everytime:
      fprintf(fp,"\t     Output Mode              :   %s   Every time\n", (FIO.IO_Output==IO_GATHER) ? "Master IO" : "Local IO");
      break;
  }
  
  
  // Output guide
  fprintf(fp,"\t     Guide cell for output    :   %d\n", GuideOut);

  
  // ログ出力
  fprintf(fp,"\n\tLogs\n");
  fprintf(fp,"\t     Unit for Output          :   %s\n", (Unit.Log == DIMENSIONAL) ? "Dimensional" : "Non-Dimensional");
  fprintf(fp,"\t     Base Logs                :   %4s  %s, %s, %s\n", 
          (Mode.Log_Base == ON)?"ON >":"OFF ", (Mode.Log_Base == ON)?HistoryName:""
          , (Mode.Log_Base == ON)?HistoryCompoName:"", (Mode.Log_Base == ON)?HistoryDomfxName:"");
  fprintf(fp,"\t     Iteration Log            :   %4s  %s\n", 
          (Mode.Log_Itr == ON)?"ON >":"OFF ", (Mode.Log_Itr == ON)?HistoryItrName:"");
  fprintf(fp,"\t     Profiling report         :   %4s  %s%s\n", 
          (Mode.Profiling != OFF)?"ON >":"OFF ", 
          (Mode.Profiling == DETAIL)? "Detail mode, ":"",
          (Mode.Profiling != OFF)?"profiling.txt":"");
          
  fprintf(fp,"\t     Wall info. Log           :   %4s  %s\n", 
          (Mode.Log_Wall == ON)?"ON >":"OFF ", (Mode.Log_Wall == ON)?HistoryWallName:"");
  
  
  // Intervals
  fprintf(fp,"\n\tIntervals\n");
  
  
  // 基本履歴のコンソール出力
  if ( !Interval[Interval_Manager::tg_console].isStep() ) {
    itm = Interval[Interval_Manager::tg_console].getIntervalTime();
    fprintf(fp,"\t     Base Info.(stdout)       :   %12.6e [%s]\n", (Unit.Log==DIMENSIONAL)?itm*Tscale:itm, (Unit.Log==DIMENSIONAL)?"sec":"-");
  }
  else {
    fprintf(fp,"\t     Base Info.(stdout)       :   %12d [step]\n", Interval[Interval_Manager::tg_console].getIntervalStep());
  }
  
  // 履歴情報のファイル出力
  if ( !Interval[Interval_Manager::tg_history].isStep() ) {
    itm = Interval[Interval_Manager::tg_history].getIntervalTime();
    fprintf(fp,"\t     Other Histories          :   %12.6e [%s]\n", (Unit.Log==DIMENSIONAL)?itm*Tscale:itm, (Unit.Log==DIMENSIONAL)?"sec":"-");
  }
  else {
    fprintf(fp,"\t     Other Histories          :   %12d [step]\n", Interval[Interval_Manager::tg_history].getIntervalStep());
  }
  
  // 瞬間値のファイル出力
  if ( !Interval[Interval_Manager::tg_instant].isStep() ) {
    itm = Interval[Interval_Manager::tg_instant].getIntervalTime();
    fprintf(fp,"\t     Instant data             :   %12.6e [%s]\n", (Unit.File==DIMENSIONAL)?itm*Tscale:itm, (Unit.File==DIMENSIONAL)?"sec":"-");
  }
  else {
    fprintf(fp,"\t     Instant data             :   %12d [step]\n", Interval[Interval_Manager::tg_instant].getIntervalStep());
  }
  
  // 平均値のファイル出力
  if ( !Interval[Interval_Manager::tg_average].isStep() ) {
    itm = Interval[Interval_Manager::tg_average].getIntervalTime();
    fprintf(fp,"\t     Averaged data            :   %12.6e [%s]\n", (Unit.File==DIMENSIONAL)?itm*Tscale:itm, (Unit.File==DIMENSIONAL)?"sec":"-");
  }
  else {
    fprintf(fp,"\t     Averaged data            :   %12d [step]\n", Interval[Interval_Manager::tg_average].getIntervalStep());
  }
  
  // サンプリング情報のファイル出力
  if ( Sampling.log == ON ) {
    if ( !Interval[Interval_Manager::tg_sampled].isStep() ) {
      itm = Interval[Interval_Manager::tg_sampled].getIntervalTime();
      fprintf(fp,"\t     Sampled data             :   %12.6e [%s]\n", (Sampling.unit==DIMENSIONAL)?itm*Tscale:itm, (Sampling.unit==DIMENSIONAL)?"sec":"-");
    }
    else {
      fprintf(fp,"\t     Sampled data             :   %12d [step]\n", Interval[Interval_Manager::tg_sampled].getIntervalStep());
    }
  }

  
  // Criteria
  fprintf(fp,"\n\tParameter of Linear Equation\n");
  ItrCtl* ICp1= &IC[ItrCtl::ic_prs_pr];  /// 圧力のPoisson反復
  ItrCtl* ICp2= &IC[ItrCtl::ic_prs_cr];  /// 圧力のPoisson反復　2回目
  ItrCtl* ICv = &IC[ItrCtl::ic_vis_cn];  /// 粘性項のCrank-Nicolson反復
  
  if ( Hide.PM_Test == ON ) {
    fprintf(fp,"\t ### Performance Test Mode >> The iteration number is fixed by Iteration max.\n\n");
  }
  
	if ( KindOfSolver != SOLID_CONDUCTION ) {
		// 1st iteration
		fprintf(fp,"\t     1st Pressure Iteration \n");
		fprintf(fp,"\t       Iteration max          :   %d\n"  ,  ICp1->get_ItrMax());
		fprintf(fp,"\t       Convergence eps        :   %9.3e\n", ICp1->get_eps());
		fprintf(fp,"\t       Coef. of Relax./Accel. :   %9.3e\n", ICp1->get_omg());
		fprintf(fp,"\t       Norm type              :   %s\n", getNormString(ICp1->get_normType()).c_str());
		printLS(fp, ICp1);
    
    if ( AlgorithmF == Flow_FS_RK_CN ) {
      fprintf(fp,"\t     2nd Pressure Iteration \n");
      fprintf(fp,"\t       Iteration max          :   %d\n"  ,  ICp2->get_ItrMax());
      fprintf(fp,"\t       Convergence eps        :   %9.3e\n", ICp2->get_eps());
      fprintf(fp,"\t       Coef. of Relax./Accel. :   %9.3e\n", ICp2->get_omg());
      fprintf(fp,"\t       Norm type              :   %s\n", getNormString(ICp2->get_normType()).c_str());
      printLS(fp, ICp2);
    }
    
    // CN iteration
		if ( (AlgorithmF == Flow_FS_AB_CN) || (AlgorithmF == Flow_FS_RK_CN) ) {
      fprintf(fp,"\n");
			fprintf(fp,"\t     Velocity CN Iteration \n");
			fprintf(fp,"\t       Iteration max           :   %d\n"  ,  ICv->get_ItrMax());
			fprintf(fp,"\t       Convergence eps         :   %9.3e\n", ICv->get_eps());
			fprintf(fp,"\t       Coef. of Relax./Accel.  :   %9.3e\n", ICv->get_omg());
			fprintf(fp,"\t       Norm type               :   %s\n", getNormString(ICv->get_normType()).c_str());
			printLS(fp, ICv);
		}
	}
	
  // for Temperature
  if ( isHeatProblem() ) {
    if ( AlgorithmH == Heat_EE_EI ) {
      ItrCtl* ICt = &IC[ItrCtl::ic_tdf_ei];  /// 温度の拡散項の反復
      fprintf(fp,"\n");
      fprintf(fp,"\t     Temperature Iteration  \n");
      fprintf(fp,"\t       Iteration max          :   %d\n"  ,  ICt->get_ItrMax());
      fprintf(fp,"\t       Convergence eps        :   %9.3e\n", ICt->get_eps());
      fprintf(fp,"\t       Coef. of Relax./Accel. :   %9.3e\n", ICt->get_omg());
			fprintf(fp,"\t       Norm type              :   %s\n", getNormString(ICt->get_normType()).c_str());
      printLS(fp, ICt);
    }
  }
  
  
  // 壁面の扱い
  fprintf(fp,"\n\tCondition of Wall\n");
  fprintf(fp,"\t     No of surface            :   %d\n", NoWallSurface);
  switch (Mode.PrsNeuamnnType) {
    case P_GRAD_ZERO:
      fprintf(fp,"\t     Pressure Gradient        :   Neumann Zero\n");
      break;
      
    case P_GRAD_NS:
      fprintf(fp,"\t     Pressure Gradient        :   Navier-Stokes\n");
      break;
      
    default:
      stamped_printf("Error: Wall treatment section\n");
      err=false;
      break;
  }
  
  switch (Mode.Wall_profile) {
    case No_Slip:
      fprintf(fp,"\t     Velocity Profile         :   No Slip\n");
      break;
      
    case Slip:
      fprintf(fp,"\t     Velocity Profile         :   Slip\n");
      break;
      
    case Log_Law:
      fprintf(fp,"\t     Velocity Profile         :   Law of Wall\n");
      break;
      
    default:
      stamped_printf("Error: Wall treatment section\n");
      err=false;
      break;
  }
  
  // Inner Iteration
  //fprintf(fp,"\tInner iteration        :   %7d\n", InnerItr);
  
  // 派生変数
  fprintf(fp, "\n\tDerived variables\n");
  
  //　全圧の出力モード
  if ( Mode.TP == ON ) {
    fprintf(fp,"\t     Total Pressure           :   ON\n");
  }
  else {
    fprintf(fp,"\t     Total Pressure           :   OFF\n");
  }
  
  //　渦度の出力モード
  if ( Mode.VRT == ON ) {
    fprintf(fp,"\t     Vorticity                :   ON\n");
  }
  else {
    fprintf(fp,"\t     Vorticity                :   OFF\n");
  }
  
  //　ヘリシティの出力モード
  if ( Mode.Helicity == ON ) {
    fprintf(fp,"\t     Helicity                 :   ON\n");
  }
  else {
    fprintf(fp,"\t     Helicity                 :   OFF\n");
  }
  
  //　速度勾配テンソルの第二不変量の出力モード
  if ( Mode.I2VGT == ON ) {
    fprintf(fp,"\t     2nd Invariant of VGT     :   ON\n");
  }
  else {
    fprintf(fp,"\t     2nd Invariant of VGT     :   OFF\n");
  }
  

  // Hidden parameter
  fprintf(fp,"\n\tHidden Parameters\n");
  if ( Hide.Scaling_Factor != 1.0 ) {
    fprintf(fp,"\t     Scaling Factor           :   %f\n", Hide.Scaling_Factor);
  }
  
  if (Hide.Change_ID != 0 ) {
    fprintf(fp,"\t     Change ID from [0] to    :   %d\n", Hide.Change_ID);
  }
  fprintf(fp,"\n");
  
  switch (Hide.Range_Limit) {
    case Range_Cutoff:
      fprintf(fp,"\t     Variable Range           :   Limit value between [0,1] in normalized value\n");
      break;
      
    case Range_Normal:
      fprintf(fp,"\t     Variable Range           :   Normal \n");
      break;
      
    default:
      break;
  }
  
  fflush(fp);

  if (err==false) Exit(0);
}

/**
 @fn void Control::printLS(FILE* fp, ItrCtl* IC)
 @brief 線形ソルバー種別の表示
 @param fp
 @param IC
 */
void Control::printLS(FILE* fp, ItrCtl* IC)
{
  switch (IC->get_LS()) {
      
    case JACOBI:
      fprintf(fp,"\t       Linear Solver          :   Jacobi Relaxation\n");
      break;
      
    case SOR:
      fprintf(fp,"\t       Linear Solver          :   Point SOR method\n");
      break;
      
    case SOR2SMA:
      fprintf(fp,"\t       Linear Solver          :   2-colored SOR SMA (Stride Memory Access)\n");
      break;
      
    case SOR2CMA:
      fprintf(fp,"\t       Linear Solver          :   2-colored SOR CMA (Consecutive Memory Access)\n");
      break;
      
    default:
      stamped_printf("Error: Linear Solver section\n");
      Exit(0);
  }
}

/**
 @fn void Control::setParameters(MaterialList* mat, CompoList* cmp, unsigned NoBaseBC, BoundaryOuter* BO, ReferenceFrame* RF)
 @brief 無次元パラメータを各種モードに応じて設定する
 @param mat
 @param cmp
 @param NoBaseBC 外部境界の基本境界条件数
 @param BO 外部境界の基本リスト
 @param rf
 @note
    - 代表長さと代表速度はパラメータで必ず与えること（読み込んだ値は変更しない）
    - 純強制対流　有次元　（代表長さ，代表速度，動粘性係数，温度拡散係数）
    -           無次元　（Pr, Re > RefV=RefL=1）
    - 熱対流　　　有次元　（代表長さ，代表速度，温度差，体膨張係数，重力加速度，動粘性係数，温度拡散係数）
    - 自然対流　　有次元　（代表長さ，代表速度，温度差，体膨張係数，重力加速度，動粘性係数，温度拡散係数）
    - 固体熱伝導　有次元　（代表長さ，温度拡散係数 > Peclet=1）？
 @see 
    - bool Control::getXML_Para_ND(void)
    - void Control::getXML_Para_Init(void)
 */
void Control::setParameters(MaterialList* mat, CompoList* cmp, unsigned NoBaseBC, BoundaryOuter* BO, ReferenceFrame* RF)
{
  REAL_TYPE rho, nyu, cp, lambda, beta, mu, snd_spd=0.0;
  REAL_TYPE c1, c2, c3;
  unsigned m;

  // get reference values
  for (unsigned n=NoBC+1; n<=NoCompo; n++) {
    if ( cmp[n].getID() == RefID ) {
      m = cmp[n].getMatOdr();
      if ( mat[m].getState() == FLUID ) {
        rho   = mat[m].P[p_density];
        mu    = mat[m].P[p_viscosity];
        nyu   = mat[m].P[p_kinematic_viscosity];
        cp    = mat[m].P[p_specific_heat];
        lambda= mat[m].P[p_thermal_conductivity];
        beta  = mat[m].P[p_vol_expansion]; // can be replaced by 1/K in the case of gas
        snd_spd = mat[m].P[p_sound_of_speed];
      }
      else {
        rho   = mat[m].P[p_density];
        cp    = mat[m].P[p_specific_heat];
        lambda = mat[m].P[p_thermal_conductivity];
      }
    }
  }
  
  RefDensity      = rho;
  RefViscosity    = mu;
  RefKviscosity   = nyu;
  RefSpecificHeat = cp;
  RefLambda       = lambda;
  RefSoundSpeed   = snd_spd;

  if (KindOfSolver == SOLID_CONDUCTION) {
    if (Unit.Param == DIMENSIONAL) {
      Peclet   = 1.0;
    }
    else {
      Hostonly_ printf("Error : Solid conduction(ND)\n");
			Exit(0);
    }
  }
  else if (KindOfSolver == FLOW_ONLY) {
    if (Unit.Param == DIMENSIONAL) {
      Reynolds = RefVelocity * RefLength / nyu;
      Prandtl  = rho * cp * nyu / lambda;
    }
  }
	else if (KindOfSolver==THERMAL_FLOW) {
    switch (Mode.Buoyancy) {
      case NO_BUOYANCY:
        if (Unit.Param == DIMENSIONAL) {
          Reynolds = RefVelocity * RefLength / nyu;
          Prandtl  = rho * cp * nyu / lambda;
          Peclet   = Prandtl * Reynolds;
        }
        else {
          Hostonly_ printf("Error : Thermal flow /wo buoyancy(ND)\n");
					Exit(0);
        }
        break;

      case BOUSSINESQ:
      case LOW_MACH:
        if (Unit.Param == DIMENSIONAL) {
          Reynolds = RefVelocity * RefLength / nyu;
          Prandtl  = rho * cp * nyu / lambda;
          Peclet   = Prandtl * Reynolds;
          c1 = beta / nyu;
          c2 = Gravity *  DiffTemp / nyu;
          c3 = RefLength*RefLength*RefLength;
          Grashof  = c1 * c2 * c3;
          Rayleigh = Prandtl * Grashof;
        }
        else {
          Hostonly_ printf("Error : Thermal flow /w buoyancy(ND)\n");
					Exit(0);
        }
        break;
    }
  }
  else if (KindOfSolver==THERMAL_FLOW_NATURAL) {
    switch (Mode.Buoyancy) {
      case BOUSSINESQ:
      case LOW_MACH:
        if (Unit.Param == DIMENSIONAL) {
          Prandtl  = rho * cp * nyu / lambda;
          Reynolds = RefVelocity * RefLength / nyu;
					Peclet   = Prandtl * Reynolds;
          c1 = beta / nyu;
          c2 = Gravity *  DiffTemp / nyu;
          c3 = RefLength*RefLength*RefLength;
          Grashof  = c1 * c2 * c3;
          Rayleigh = Prandtl * Grashof;
        }
        else {
          Hostonly_ printf("Error : Natural Convection(ND)\n");
					Exit(0);
        }
        break;
      default:
        break;
		}
	}
	else { // CONJUGATE_HEAT_TRANSFER
		;
	}

  if (Mode.PDE == PDE_EULER) Reynolds=1.0e23;

  Mach = RefVelocity / RefSoundSpeed;
  
  // タイミングパラメータの無次元化
  Tscale = RefLength / RefVelocity;

  // 入力モードが有次元の場合に，無次元に変換
  if ( Unit.Param == DIMENSIONAL ) {
    Interval[Interval_Manager::tg_compute].normalizeInterval(Tscale);
    Interval[Interval_Manager::tg_console].normalizeInterval(Tscale);
    Interval[Interval_Manager::tg_history].normalizeInterval(Tscale);
    Interval[Interval_Manager::tg_instant].normalizeInterval(Tscale);
    Interval[Interval_Manager::tg_average].normalizeInterval(Tscale);
    Interval[Interval_Manager::tg_accelra].normalizeInterval(Tscale);
    Interval[Interval_Manager::tg_avstart].normalizeInterval(Tscale);
  }
  
  // Reference frame
  RF->setAccel( Interval[Interval_Manager::tg_accelra].getIntervalTime() );
  
  if ( Unit.Param == DIMENSIONAL ) {
    double g[3];
    RF->copyGridVel(g);
    g[0] /= (double)RefVelocity;
    g[1] /= (double)RefVelocity;
    g[2] /= (double)RefVelocity;
    RF->setGridVel(g);
  }
  
  // コンポーネントの指定速度
  for (int n=1; n<=NoBC; n++) {
    if ( (cmp[n].getType()==SPEC_VEL_WH) || (cmp[n].getType()==SPEC_VEL) ) {
			if ( cmp[n].isPolicy_Massflow() ) { //ポリシーが流量の場合
				if ( Unit.Param == DIMENSIONAL ) {
					cmp[n].set_Velocity( cmp[n].get_Massflow() / cmp[n].area );  // attenltion! Velocity and MassFlow use same variable
				}
				else {
					cmp[n].set_Velocity( cmp[n].get_Massflow()*RefVelocity*RefLength*RefLength / cmp[n].area );
				}
			}
      
			// 流量指定のときのみ，ca[]に有次元速度パラメータを保存  >> 速度指定の場合には，parseBC::getXML_IBC_SpecVel()で設定済み
      if ( cmp[n].isPolicy_Massflow() ) {
        if ( Unit.Param != DIMENSIONAL ) {
          Hostonly_ stamped_printf("Error: Non-dimensional condition\n");
          Exit(0);
        }
        else {
          cmp[n].ca[CompoList::amplitude] = cmp[n].get_Velocity();
          cmp[n].ca[CompoList::bias] = cmp[n].ca[CompoList::bias] / cmp[n].area; // dimensional velocity
        }
      }
		}
  } // end of NoBC
	
  // 発熱密度の計算(有次元) -- 発熱量と発熱密度
  REAL_TYPE a, vol;
  a = dh*RefLength;
  vol = a*a*a;
  
  for (int n=1; n<=NoBC; n++) {
    if ( cmp[n].getType()==HEAT_SRC ) {
      m = cmp[n].getMatOdr();
      if (cmp[n].get_sw_Heatgen() == CompoList::hsrc_watt) {
        cmp[n].set_HeatDensity( cmp[n].get_HeatValue() / ((REAL_TYPE)cmp[n].getElement()*vol) );
      }
      else { // 発熱密度
        cmp[n].set_HeatValue( cmp[n].get_HeatDensity() * ((REAL_TYPE)cmp[n].getElement()*vol) );
      }
    }
  }
  
  // Darcy係数（単媒質）
  // C[0-2]; 有次元，C[3-5]; 無次元係数
  REAL_TYPE ki;
  for (int n=1; n<=NoBC; n++) {
    if ( cmp[n].getType()==DARCY ) {
      m = cmp[n].getMatOdr();
      ki = (mat[m].P[p_viscosity]*RefLength) / (mat[m].P[p_density]*RefVelocity);
      cmp[n].ca[3] = ki / cmp[n].ca[0];
      cmp[n].ca[4] = ki / cmp[n].ca[1];
      cmp[n].ca[5] = ki / cmp[n].ca[2];
    }    
  }
  
  // Pressure Loss
  REAL_TYPE DensityOfMedium, cf[6];
  for (int n=1; n<=NoBC; n++) {
    if ( cmp[n].getType()==HEX ) {
      m = cmp[n].getMatOdr();
      for (unsigned i=0; i<6; i++) cf[i] = cmp[n].ca[i];
      
      // 流量と圧力損失量計算の有次元化の係数
      cmp[n].set_CoefMassflow( RefLength * RefLength * RefVelocity );
      cmp[n].set_CoefPrsLoss( cf[5] * RefDensity * RefVelocity * RefVelocity / RefLength );
      
      // Normalize
      if ( cmp[n].getPrsUnit() == CompoList::unit_mmAq ) {
        // Water: T=300K, p=101.325kPa > 996.62 kg/m^3
        DensityOfMedium = 996.62;
        convertHexCoef(cf, DensityOfMedium);
      }
      else if ( cmp[n].getPrsUnit() == CompoList::unit_mmHg ) {
        // Mercury: T=300K > 13538 kg/m^3
        DensityOfMedium = 13538.0;
        convertHexCoef(cf, DensityOfMedium);
      }
      else if ( cmp[n].getPrsUnit() == CompoList::unit_Pa ) {
        convertHexCoef(cf);
      }
      else if ( cmp[n].getPrsUnit() == CompoList::unit_NonDimensional ) {
        cf[4] /= RefVelocity; // 無次元の場合には単位変更のみ
        cf[5] *= (1e-3/RefLength);
      }
      
      for (unsigned i=0; i<6; i++) cmp[n].ca[i] = cf[i];
    }
    
  }

  // 外部境界面の速度の指定パラメータを有次元化
  if ( Unit.Param == NONDIMENSIONAL ) {
    for (unsigned n=0; n<NoBaseBC; n++) {
      switch ( BO[n].get_BCtype() ) {
        case OBC_WALL:
        case OBC_SPEC_VEL:
          BO[n].ca[CompoList::amplitude] *= RefVelocity;
          BO[n].ca[CompoList::frequency] *= (RefVelocity/RefLength);
          //BO[n].ca[CompoList::initphase] radは有次元化不要
          BO[n].ca[CompoList::bias]      *= RefVelocity;
          break;
          
        default:
          break;
      }
    }
  }
  
  // 外部境界面の圧力の有次元化
  if ( Unit.Param == NONDIMENSIONAL ) {
    for (unsigned n=0; n<NoBaseBC; n++) {
      switch ( BO[n].get_BCtype() ) {
        case OBC_OUTFLOW:
        case OBC_TRC_FREE:
          if ( BO[n].get_pType() == P_DIRICHLET ) {
            BO[n].p = FBUtility::convND2D_P(BO[n].p, BasePrs, RefDensity, RefVelocity, Unit.Prs); 
          }          
          break;
          
        case OBC_PERIODIC:
          if ( BO[n].get_PrdcMode() != BoundaryOuter::prdc_Simple ) { // Dirichlet or Bidirectionalを指定の場合
            BO[n].p = FBUtility::convND2D_P(BO[n].p, BasePrs, RefDensity, RefVelocity, Unit.Prs); 
          }
          break;
          
        default:
          break;
      }
    }
  }
  
  // 初期条件の値の有次元化
  if ( Unit.Param == NONDIMENSIONAL )  {
    iv.Pressure = FBUtility::convND2D_P(iv.Pressure, BasePrs, RefDensity, RefVelocity, Unit.Prs);
    iv.Density *= RefDensity;
		iv.VecU    *= RefVelocity;
		iv.VecV    *= RefVelocity;
		iv.VecW    *= RefVelocity;
    if ( isHeatProblem() ) {
			iv.Temperature = FBUtility::convND2Kelvin(iv.Temperature, BaseTemp, DiffTemp); //内部表現をKelvinに
    }
	}
	else {
		if ( isHeatProblem() ) {
			iv.Temperature = FBUtility::convTemp2K(iv.Temperature, Unit.Temp);
    }
	}
}

/**
 @fn unsigned Control::getNoInFiles(const char* key)
 @brief InFileに記述された領域数を数える
 @retval attr=keyの数
 */
unsigned Control::getNoInFiles(const char* key)
{
	const SklCfgInFile* infile = CF->GetInFileFirst();
	const char *attr;
	unsigned nd=0;
  
	// count subdomain
	while ( infile ) {
    if ( !(attr = infile->GetAttr()) ) {
      stamped_printf("\tParsing error : InFile description\n");
      return 0;
    }
    if ( !strcasecmp(attr, key) ) nd++;
    infile = CF->GetInFileNext(infile);
  }
	return nd;
}

/**
 @fn char* Control::getVoxelFileName(void)
 @brief InFileに記述されたボクセルファイル名を取得
 */
const char* Control::getVoxelFileName(void)
{
  SklCfgInFile* infile = (SklCfgInFile*)CF->GetInFileFirst();
	const char *format, *fname=NULL;
  
  // ファイル名の取得
  if ( !strcasecmp(infile->GetAttr(), "SphereSVX") ) {
    vxFormat = Sphere_SVX;
  }
  else if ( !strcasecmp(infile->GetAttr(), "SphereSBX") ) {
    vxFormat = Sphere_SBX;
  }
  else {
    stamped_printf("\tParsing error : InFile attr [%s] is invalid.\n", infile->GetAttr());
    Exit(0);
  }
  
  // フォーマットのチェック
	if ( !(format = infile->GetFormat()) ) {
		stamped_printf("\tParsing error : Invalid format at InFile description\n");
		Exit(0);
	}
  if ( !strcasecmp(format, "svx") ) {
    if ( vxFormat != Sphere_SVX ) {
      stamped_printf("\tParsing error : Specification of Voxel file format is not consistent.\n");
      Exit(0);
    }
  }
	else if ( !strcasecmp(format, "sbx") ) {
    if ( vxFormat != Sphere_SBX ) {
      stamped_printf("\tParsing error : Specification of Voxel file format is not consistent.\n");
      Exit(0);
    }
  }
  else {
		stamped_printf("\tParsing error : Invalid format for input voxel: %s\n", format);
		Exit(0);
	}
  
  //ファイル名取得
  if( !(fname = infile->GetFileName()) ) {
    stamped_printf("\tParsing error : InFile description\n");
    Exit(0);
  }
  
  // sbx/svxファイルはシリアル入力のみ対応なので，マルチモードはoff
  infile->UnsetMultiInput();

  return fname;
}

/**
 @fn void Control::convertHexCoef(REAL_TYPE* cf, REAL_TYPE Density)
 @brief 熱交換器パラメータの変換（水と水銀）
 @param cf パラメータ値
 @param Density ヘッドの単位
 */
void Control::convertHexCoef(REAL_TYPE* cf, REAL_TYPE Density)
{
  REAL_TYPE cc[6], s;
  
  s = (Density*RefLength*Gravity)/(RefDensity*cf[5]*1e-3);
  cc[0] = s*cf[0];                            // c1
  cc[1] = s*cf[1]/RefVelocity;                // c2
  cc[2] = s*cf[2]/(RefVelocity*RefVelocity);  // c3
  cc[3] = s*cf[3];                            // c4
  cc[4] = cf[4]/RefVelocity;                  // thresholdの無次元値
  cc[5] = cf[5]*1e-3/RefLength;               // thicknessの無次元値，入力はmm
  
  for (int i=0; i<6; i++) cf[i] = cc[i];
}

/**
 @fn void Control::convertHexCoef(REAL_TYPE* cf)
 @brief 熱交換器パラメータの変換（Pa）
 @param cf パラメータ値
 */
void Control::convertHexCoef(REAL_TYPE* cf)
{
  REAL_TYPE cc[6], s;
  
  s = RefLength/(RefDensity*cf[5]*1e-3);
  cc[0] = s*cf[0];
  cc[1] = s*cf[1]/RefVelocity;
  cc[2] = s*cf[2]/(RefVelocity*RefVelocity);
  cc[3] = s*cf[3];
  cc[4] = cf[4]/RefVelocity;
  cc[5] = cf[5]*1e-3/RefLength;
  
  for (int i=0; i<6; i++) cf[i] = cc[i];
}

/**
 @fn unsigned Control::countCompo(CompoList* cmp, unsigned label)
 @brief labelのコンポーネント数を返す
 @param cmp
 @param label コンポーネントID
 */
unsigned Control::countCompo(CompoList* cmp, unsigned label)
{
  unsigned cnt=0;
  for (unsigned i=1; i<=NoBC; i++) {
    if ( cmp[i].getType() == label ) cnt++;
  }
  return cnt;
}
