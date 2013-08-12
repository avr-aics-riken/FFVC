//##################################################################################
//
// Flow Base class
//
// Copyright (c) 2007-2011 VCAD System Research Program, RIKEN.
// All rights reserved.
//
// Copyright (c) 2011-2013 Institute of Industrial Science, The University of Tokyo.
// All rights reserved.
//
// Copyright (c) 2012-2013 Advanced Institute for Computational Science, RIKEN.
// All rights reserved.
//
//##################################################################################

/**
 @file   Control.C
 @brief  FlowBase Control class
 @author kero
 */

#include "Control.h"


// #################################################################
// 時間積分幅とKindOfSolver種別の整合性をチェック
bool DTcntl::chkDtSelect()
{
  switch (KOS) {
    case FLOW_ONLY:
    case THERMAL_FLOW:
    case THERMAL_FLOW_NATURAL:
    case CONJUGATE_HEAT_TRANSFER:
      break;
      
    case SOLID_CONDUCTION:
      if ( (scheme != dt_direct) && (scheme != dt_dfn) ) {
        return false;
      }
      break;
  }
  return true;
}


// #################################################################
//各種モードに対応する時間積分幅を設定する
int DTcntl::set_DT(const double vRef)
{
  double dtC, dtD, a, b;
  
  switch ( scheme ) 
  {
    case dt_direct:
      deltaT = CFL;
      break;
      
    case dt_cfl_ref_v:
      if ( KOS == SOLID_CONDUCTION ) return 1;
      deltaT = dtCFL( vRef );
      break;
      
    case dt_cfl_max_v:
      if ( KOS == SOLID_CONDUCTION ) return 2;
      deltaT = dtCFL(vRef);
      break;
      
    case dt_dfn:
      if ( KOS != SOLID_CONDUCTION ) return 3;
      deltaT = dtDFN( Peclet );
      break;
      
    case dt_cfl_dfn_ref_v:
      switch (KOS)
    {
        case FLOW_ONLY:
          dtC = dtCFL( vRef );
          dtD = dtDFN( Reynolds );
          deltaT = (dtC > dtD) ? dtD : dtC;
          break;
          
        case THERMAL_FLOW:
        case THERMAL_FLOW_NATURAL:
        case CONJUGATE_HEAT_TRANSFER:
          dtC = dtCFL( vRef );
          a = dtDFN( Reynolds );
          b = dtDFN( Peclet );
          dtD = (a > b) ? b : a;
          deltaT = (dtC > dtD) ? dtD : dtC;
          break;
          
        case SOLID_CONDUCTION:
          return 4;
          break;
      }
      break;
      
    case dt_cfl_dfn_max_v:
      return 6;
      break;
      
    case dt_cfl_max_v_cp:
      return 7;
      break;
  }
  
  return 0;
}


// #################################################################
// Δtのスキームを設定する
bool DTcntl::set_Scheme(const char* str, const double val)
{
  if ( !str ) return false;
  
  if ( !strcasecmp(str, "Direct") )
  {
    scheme = dt_direct;
  }
  else if ( !strcasecmp(str, "CflReferenceVelocity") )
  {
    scheme = dt_cfl_ref_v;
  }
  //else if ( !strcasecmp(str, "CFL_MaxV") ) {
  //  scheme = dt_cfl_max_v;
  //}
  else if ( !strcasecmp(str, "Diffusion") )
  {
    scheme = dt_dfn;
  }
  else if ( !strcasecmp(str, "CflDiffusionReferenceVelocity") )
  {
    scheme = dt_cfl_dfn_ref_v;
  }
  //else if ( !strcasecmp(str, "CFL_DFN_MaxV") ) {
  //  scheme = dt_cfl_dfn_max_v;
  //}
  //else if ( !strcasecmp(str, "CFL_MaxV_CP") ) {
  //  scheme = dt_cfl_max_v_cp;
  //}
  else
  {
    return false;
  }
  
  CFL = val;
  
  return true;
}


// #################################################################
// 基本変数をコピー
void DTcntl::set_Vars(const int m_kos, const int m_mode, const double m_dh, const double re, const double pe)
{
  KOS      = m_kos;
  mode     = m_mode;
  dh       = m_dh;
  Reynolds = re;
  Peclet   = pe; 
}







// #################################################################
// 加速時間をセットする
void ReferenceFrame::setAccel(const double m_timeAccel)
{
  TimeAccel = m_timeAccel;
}


// #################################################################
// 参照フレームの種類をセットする
void ReferenceFrame::setFrame(const int m_frame)
{
  Frame = m_frame;
}


// #################################################################
// 格子速度成分の単位方向ベクトルをセットする
void ReferenceFrame::setGridVel(const double* m_Gvel)
{
  GridVel[0] = m_Gvel[0];
  GridVel[1] = m_Gvel[1];
  GridVel[2] = m_Gvel[2];
}


// #################################################################
// 参照速度を計算する
void ReferenceFrame::setV00(const double time, const bool init) 
{
  if (init == true) 
  {
    v00[0]=1.0;
  }
  else 
  {
    if ( TimeAccel == 0.0 )
    {
      v00[0] = 1.0;
    }
    else 
    {
      const double c_pai = (double)(2.0*asin(1.0));
      v00[0] = 0.5*(1.0-cos(c_pai*time/(TimeAccel)));
      if ( time > (TimeAccel) ) v00[0] = 1.0;
    }
  }
  
  double u0 = v00[0];
  
  switch (Frame)
  {
    case frm_static:
      v00[1] = 0.0;
      v00[2] = 0.0;
      v00[3] = 0.0;
      break;
      
    case frm_translation:
      v00[1] = u0*GridVel[0];  // v0x
      v00[2] = u0*GridVel[1];  // v0y
      v00[3] = u0*GridVel[2];  // v0z
      break;
      
    case frm_rotation:
      break;
  }
  
}



// #################################################################
/**
 * @brief 熱交換器パラメータの変換（水と水銀）
 * @param [out] cf      パラメータ値
 * @param [in]  Density ヘッドの単位
 */
void Control::convertHexCoef(REAL_TYPE* cf, const REAL_TYPE Density)
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


// #################################################################
/**
 * @brief 熱交換器パラメータの変換（Pa）
 * @param [out] cf パラメータ値
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



// #################################################################
// 反復の収束判定パラメータをcopy
// @see getIteration()
void Control::copyCriteria(IterationCtl& IC, const string name)
{
  
  for (int i=0; i<NoBaseLS; i++)
  {
    if ( !strcasecmp( name.c_str(), Criteria[i].getAlias().c_str() ))
    {
      IC.copy(&Criteria[i]);
    }
  }
  
}



// #################################################################
// 制御，計算パラメータ群の表示
void Control::displayParams(FILE* mp, FILE* fp, IterationCtl* IC, DTcntl* DT, ReferenceFrame* RF, MediumList* mat)
//void Control::displayParams(FILE* mp, FILE* fp, IterationCtl* IC, DTcntl* DT, ReferenceFrame* RF, MediumList* mat, FileIO_PLOT3D_WRITE* FP3DW)
{
  /* 20130611
  printSteerConditions(mp, IC, DT, RF, FP3DW);
  printSteerConditions(fp, IC, DT, RF, FP3DW);
   */
  
  printSteerConditions(mp, IC, DT, RF);
  printSteerConditions(fp, IC, DT, RF);
  printParaConditions(mp, mat);
  printParaConditions(fp, mat);
  printInitValues(mp);
  printInitValues(fp);
}


// #################################################################
// MediumList中に登録されているkeyに対するIDを返す。発見できない場合はzero 
int Control::findIDfromLabel(const MediumList* mat, const int Nmax, const std::string key)
{
  std::string str = key;

  for (int i=1; i<=Nmax; i++) 
  {
    if ( !strcasecmp(str.c_str(), mat[i].getAlias().c_str()) ) return i;
  }

  return 0;
}



// #################################################################
/**
 * @brief 制御，計算パラメータ群の取得
 * @param [in] DT     DTctlクラス ポインタ
 //* @param [in] FP3DR  PLOT3D READクラス ポインタ
 //* @param [in] FP3DW  PLOT3D WRITEクラス ポインタ
 * @note 他のパラメータ取得に先んじて処理しておくもの
 */
//void Control::get_Steer_1(DTcntl* DT, FileIO_PLOT3D_READ* FP3DR, FileIO_PLOT3D_WRITE* FP3DW) // 20130611
void Control::get1stParameter(DTcntl* DT)
{
  
  // ソルバーの具体的な種類を決めるパラメータを取得し，ガイドセルの値を設定する
  getSolverProperties();
  
  // 指定単位を取得
  getUnit();
  
  // Reference parameter needs to be called before setDomainParameter();
  // パラメータの取得，代表値に関するもの．
  getReference();
  
  // 時間制御パラメータ
  getTimeControl(DT);
  
  
  // モニターのON/OFF 詳細パラメータはMonitorList::getMonitor()で行う
  getMonitorList();
  
  // ファイル入出力に関するパラメータ
  getFieldData();
  
  /* 20130611 PLOT3Dファイル入出力に関するパラメータ
   if (FIO.Format == plt3d_fmt) get_PLOT3D(FP3DR,FP3DW);
   */
  
}


// #################################################################
/**
 * @brief 制御，計算パラメータ群の取得
 * @param [in] RF  ReferenceFrameクラス
 */
void Control::get2ndParameter(ReferenceFrame* RF)
{
  // パラメータを取得する
  if ( Unit.Param == NONDIMENSIONAL )
  {
    if ( KindOfSolver == FLOW_ONLY ) getDimensionlessParameter();
  }
  
  // Reference frame information
  getReferenceFrame(RF);
  
  // getAverageOption(); >> getFieldData()へ
  
  // 圧力ノイマン条件のタイプ >> get_Log()よりも先に
  getWallType();
  
  // Log >> getIteration()よりも前に
  getLog();
  
  getIteration();
  
  getTurbulenceModel();
  
  
  // ラフな初期値を使い、リスタートするモード指定 >> FileIO
  getStartCondition();
  
  
  // 形状近似
  getShapeApproximation();
  
  
  getApplicationControl();

  
}




// #################################################################
/**
 * @brief アプリケーションのパラメータを取得する
 */
void Control::getApplicationControl()
{
  string str;
  string label;
  
  label = "/ApplicationControl/Operator";
  
  if ( !(tpCntl->GetValue(label, &str )) )
  {
	  Exit(0);
  }
  OperatorName = str;
  
  
  // パラメータチェックフラグ
  label = "/ApplicationControl/CheckParameter";
  
  if ( !(tpCntl->GetValue(label, &str )) )
  {
	  Exit(0);
  }
  
  if     ( !strcasecmp(str.c_str(), "On") )   CheckParam = ON;
  else if( !strcasecmp(str.c_str(), "Off") )  CheckParam = OFF;
  else
  {
    Hostonly_ printf("\tInvalid keyword is described for '%s'\n", label.c_str());
    Exit(0);
  }
  
  int ct = 0;
  
  // Performance Testの設定（隠しオプション）
  // 'PerformanceTest'の文字列チェックはしないので注意して使うこと
  Hide.PM_Test = OFF;
  label="/ApplicationControl/PerformanceTest";
  
  if ( !(tpCntl->GetValue(label, &str )) )
  {
	  ;
  }
  else
  {
    if     ( !strcasecmp(str.c_str(), "on") )   Hide.PM_Test = ON;
    else if( !strcasecmp(str.c_str(), "off") )  Hide.PM_Test = OFF;
    else
    {
      Hostonly_ stamped_printf("\tInvalid keyword is described for '%s'\n", label.c_str());
      Exit(0);
    }
  }
  
  
  // 変数の範囲制限モードを取得
  Hide.Range_Limit = Range_Normal;
  label="/ApplicationControl/VariableRange";
  
  if ( !(tpCntl->GetValue(label, &str )) )
  {
    ;
  }
  else
  {
    if     ( !strcasecmp(str.c_str(), "normal") ) Hide.Range_Limit = Range_Normal;
    else if( !strcasecmp(str.c_str(), "cutoff") ) Hide.Range_Limit = Range_Cutoff;
    else
    {
      Hostonly_ stamped_printf("\tInvalid keyword is described for '%s\n", label.c_str());
      Exit(0);
    }
  }
  
  
  // Cell IDのゼロを指定IDに変更するオプションを取得する（隠しオプション）
  label = "/ApplicationControl/ChangeID";
  if ( !(tpCntl->GetValue(label, &ct )) )
  {
	  ;
  }
  else
  {
    if ( ct < 0 )
    {
      Hostonly_ printf("Error : ID should be positive [%d]\n", ct);
      Exit(0);
    }
    else
    {
      Hide.Change_ID = ct;
    }
  }
  
}


// #################################################################
/**
 * @brief 平均値操作に関するパラメータを取得する
 * @note パラメータは，setParameters()で無次元して保持
 * @see getTimeControl()でMode.Averageをセット
 */
void Control::getAverageOption()
{
  double ct;
  string str;
  string label;
  
  
  if ( Mode.Average == ON )
  {
	  label = "/Output/Data/AveragedVariables/IntervalType";
    
	  if ( !(tpCntl->GetValue(label, &str )) )
    {
		  Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
		  Exit(0);
	  }
	  else
    {
		  if     ( !strcasecmp(str.c_str(), "step") )
      {
        Interval[Interval_Manager::tg_average].setMode_Step();
		  }
		  else if( !strcasecmp(str.c_str(), "time") )
      {
        Interval[Interval_Manager::tg_average].setMode_Time();
		  }
		  else
      {
			  Hostonly_ stamped_printf("\tParsing error : Invalid keyword for '%s'\n", label.c_str());
			  Exit(0);
		  }
    }
      
    label="/Output/Data/AveragedVariables/Interval";
    
    if ( !(tpCntl->GetValue(label, &ct )) )
    {
      Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
      Exit(0);
    }
    else
    {
      Interval[Interval_Manager::tg_average].setInterval(ct);
    }
  }
}




// #################################################################
// 計算内部領域の全セル数を返す
REAL_TYPE Control::getCellSize(const int* G_size)
{
  return (REAL_TYPE)G_size[0] * (REAL_TYPE)G_size[1] * (REAL_TYPE)G_size[2];
}



// #################################################################
// @brief 対流項スキームのパラメータを取得する
void Control::getConvection()
{
  
  int ct;
  string str;
  string label;
  
  // scheme
  label="/ConvectionTerm/scheme";
  
  if ( !(tpCntl->GetValue(label, &str )) )
  {
	  Hostonly_ stamped_printf("\tParsing error : Invalid char* value for '%s'\n", label.c_str());
	  Exit(0);
  }
  
  if     ( !strcasecmp(str.c_str(), "O1Upwind") )    CnvScheme = O1_upwind;
  else if( !strcasecmp(str.c_str(), "O3muscl") )     CnvScheme = O3_muscl;
  else if( !strcasecmp(str.c_str(), "O2central") )   CnvScheme = O2_central;
  else if( !strcasecmp(str.c_str(), "O4central") ) { CnvScheme = O4_central; Exit(0); }  // not yet implemented
  else
  {
    Hostonly_ stamped_printf("\tInvalid keyword is described for '%s'\n", label.c_str());
    Exit(0);
  }
  
  // Limiter
  if ( CnvScheme == O3_muscl )
  {
		label="/ConvectionTerm/limiter";
    
		if ( !(tpCntl->GetValue(label, &str )) )
    {
			Hostonly_ stamped_printf("\tParsing error : Invalid char* value for '%s'\n", label.c_str());
			Exit(0);
		}
    
		if     ( !strcasecmp(str.c_str(), "NoLimiter") )  Limiter = No_Limiter;
		else if( !strcasecmp(str.c_str(), "Minmod") )     Limiter = MinMod;
		else
    {
			Hostonly_ stamped_printf("\tInvalid keyword is described for '%s'\n", label.c_str());
			Exit(0);
		}
  }
}


// #################################################################
// @brief 派生して計算する変数のオプションを取得する
void Control::getDerived()
{
  string str;
  string label;
  
  // 全圧
  label="/Output/Data/DerivedVariables/TotalPressure";
  
  if ( !(tpCntl->GetValue(label, &str )) )
  {
	  Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
	  Exit(0);
  }
  
  if     ( !strcasecmp(str.c_str(), "on") )  Mode.TP = ON;
  else if( !strcasecmp(str.c_str(), "off") ) Mode.TP = OFF;
  else
  {
    Hostonly_ stamped_printf("\tInvalid keyword is described for '%s'\n", label.c_str());
    Exit(0);
  }
  
  // 渦度ベクトル
  label="/Output/Data/DerivedVariables/Vorticity";
  
  if ( !(tpCntl->GetValue(label, &str )) )
  {
	  Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
	  Exit(0);
  }
  
  if     ( !strcasecmp(str.c_str(), "on") )  Mode.VRT = ON;
  else if( !strcasecmp(str.c_str(), "off") ) Mode.VRT = OFF;
  else
  {
    Hostonly_ stamped_printf("\tInvalid keyword is described for '%s'\n", label.c_str());
    Exit(0);
  }
  
  // 速度勾配テンソルの第2不変量
  label="/Output/Data/DerivedVariables/Qcriterion";
  
  if ( !(tpCntl->GetValue(label, &str )) )
  {
	  Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
	  Exit(0);
  }
  
  if     ( !strcasecmp(str.c_str(), "on") )  Mode.I2VGT = ON;
  else if( !strcasecmp(str.c_str(), "off") ) Mode.I2VGT = OFF;
  else
  {
    Hostonly_ stamped_printf("\tInvalid keyword is described for '%s'\n", label.c_str());
    Exit(0);
  }
  
  // ヘリシティ
  label="/Output/Data/DerivedVariables/Helicity";
  
  if ( !(tpCntl->GetValue(label, &str )) )
  {
	  Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
	  Exit(0);
  }
  
  if     ( !strcasecmp(str.c_str(), "on") )  Mode.Helicity = ON;
  else if( !strcasecmp(str.c_str(), "off") ) Mode.Helicity = OFF;
  else
  {
    Hostonly_ stamped_printf("\tInvalid keyword is described for '%s'\n", label.c_str());
    Exit(0);
  }
  
  
  // FaceVelocity hidden parameter
  label="/Output/Data/DerivedVariables/FaceVelocity";
  
  if ( !(tpCntl->GetValue(label, &str )) )
  {
	  ;
  }
  else
  {
    if     ( !strcasecmp(str.c_str(), "on") )  Mode.FaceV = ON;
    else if( !strcasecmp(str.c_str(), "off") ) Mode.FaceV = OFF;
    else
    {
      Hostonly_ stamped_printf("\tInvalid keyword is described for '%s'\n", label.c_str());
      Exit(0);
    }
  }
  
}


// #################################################################
// @brief 無次元パラメータを各種モードに応じて設定する
// @note 純強制対流　有次元　（代表長さ，代表速度，動粘性係数，温度拡散係数）
//                 無次元　（Pr, Re > RefV=RefL=1）
// @see bool Control::setParameters(MediumList* mat, CompoList* cmp)
void Control::getDimensionlessParameter()
{
  REAL_TYPE ct;
  string label;
  
  label="/Reference/Reynolds";
  
  if ( !(tpCntl->GetValue(label, &ct )) )
  {
    Hostonly_ printf("\tParsing error in '%s'\n", label.c_str());
	  Exit(0);
  }
  Reynolds=(REAL_TYPE)ct;
  
  
  label="/Reference/Prandtl";
  
  if ( !(tpCntl->GetValue(label, &ct )) )
  {
    Hostonly_ printf("\tParsing error in '%s'\n", label.c_str());
	  Exit(0);
  }
  Prandtl=(REAL_TYPE)ct;
}



// #################################################################
// @brief ファイル入出力に関するパラメータを取得し，sphフォーマットの出力の並列モードを指定する．
// @note インターバルパラメータは，setParameters()で無次元して保持
void Control::getFieldData()
{
  
  REAL_TYPE f_val=0.0;
  REAL_TYPE ct;
  string str;
  string label, leaf;

  
  // Default setting
  FIO.IOmode = IO_DISTRIBUTE;
  
  // 逐次実行の場合には、強制的に IO_GATHER
  if ( (Parallelism == Serial) || (Parallelism == OpenMP) )
  {
    FIO.IOmode = IO_GATHER;
  }
  
  
  
  // 基本変数の瞬時値データ
  
  // ファイルフォーマット
  label = "/Output/Data/BasicVariables/Format";
  
  if ( !(tpCntl->GetValue(label, &str )) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
	  Exit(0);
  }
  
  if     ( !strcasecmp(str.c_str(), "sph") )     FIO.Format = sph_fmt;
  else if( !strcasecmp(str.c_str(), "bov") )     FIO.Format = bov_fmt;
  else if( !strcasecmp(str.c_str(), "plot3d") )  FIO.Format = plt3d_fmt;
  else
  {
    Hostonly_ stamped_printf("\tInvalid keyword is described for '%s'\n", label.c_str());
    Exit(0);
  }
  
  
  // インターバル
  label = "/Output/Data/BasicVariables/IntervalType";
  
  if ( !(tpCntl->GetValue(label, &str )) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    Exit(0);
  }
  else
  {
    if     ( !strcasecmp(str.c_str(), "step") )
    {
      Interval[Interval_Manager::tg_basic].setMode_Step();
    }
    else if( !strcasecmp(str.c_str(), "time") )
    {
      Interval[Interval_Manager::tg_basic].setMode_Time();
    }
    else
    {
      Hostonly_ stamped_printf("\tParsing error : Invalid keyword for '%s'\n", label.c_str());
      Exit(0);
    }
    
    label="/Output/Data/BasicVariables/Interval";
    
    if ( !(tpCntl->GetValue(label, &f_val )) )
    {
      Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
      Exit(0);
    }
    else
    {
      Interval[Interval_Manager::tg_basic].setInterval((double)f_val);
    }
  }
  
  
  
  // 派生変数
  
  /* ファイルフォーマット >> 基本変数と同じ
  label = "/Output/Data/DerivedVariables/Format";
  
  if ( !(tpCntl->GetValue(label, &str )) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
	  Exit(0);
  }
  
  if     ( !strcasecmp(str.c_str(), "sph") )
  {
    FIO.Format = sph_fmt;
    file_fmt_ext = "sph";
  }
  else if( !strcasecmp(str.c_str(), "bov") )
  {
    FIO.Format = bov_fmt;
    file_fmt_ext = "dat";
  }
  else if( !strcasecmp(str.c_str(), "plot3d") )
  {
    FIO.Format = plt3d_fmt;
    file_fmt_ext = ""; // plot3dの場合は,ファイルがいくつかあるので別に処理
  }
  else
  {
    Hostonly_ stamped_printf("\tInvalid keyword is described for '%s'\n", label.c_str());
    Exit(0);
  }*/
  
  
  // インターバル
  label = "/Output/Data/DerivedVariables/IntervalType";
  
  if ( !(tpCntl->GetValue(label, &str )) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    Exit(0);
  }
  else
  {
    if     ( !strcasecmp(str.c_str(), "step") )
    {
      Interval[Interval_Manager::tg_derived].setMode_Step();
    }
    else if( !strcasecmp(str.c_str(), "time") )
    {
      Interval[Interval_Manager::tg_derived].setMode_Time();
    }
    else
    {
      Hostonly_ stamped_printf("\tParsing error : Invalid keyword for '%s'\n", label.c_str());
      Exit(0);
    }
    
    label="/Output/Data/DerivedVariables/Interval";
    
    if ( !(tpCntl->GetValue(label, &f_val )) )
    {
      Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
      Exit(0);
    }
    else
    {
      Interval[Interval_Manager::tg_derived].setInterval((double)f_val);
    }
  }
  
  getDerived();
  
  switch ( FIO.Format )
  {
    case sph_fmt:
      getFormat_sph();
      break;
      
    case bov_fmt:
      break;
      
    case plt3d_fmt:
      getFormat_plot3d();
      break;
  }
  
  getAverageOption();
  
  
  // ボクセルファイル出力（隠しオプション）
  FIO.IO_Voxel = OFF;
  label = "/Output/VoxelOutput";
  
  if ( !(tpCntl->GetValue(label, &str )) )
  {
    ;
  }
  
  if     ( !strcasecmp(str.c_str(), "svx") )  FIO.IO_Voxel = Sphere_SVX;
  else if( !strcasecmp(str.c_str(), "off") )  FIO.IO_Voxel = OFF;
  else
  {
    Hostonly_ stamped_printf("\tInvalid keyword is described for '%s'\n", label.c_str());
    Exit(0);
  }
  
  
  // デバッグ用のdiv(u)の出力指定（隠しオプション）
  FIO.Div_Debug = OFF;
  label = "/Output/DebugDivergence";
  
  if ( !(tpCntl->GetValue(label, &str )) )
  {
    ;
  }
  
  if     ( !strcasecmp(str.c_str(), "on") )    FIO.Div_Debug = ON;
  else if( !strcasecmp(str.c_str(), "off") )   FIO.Div_Debug = OFF;
  else
  {
    Hostonly_ stamped_printf("\tInvalid keyword is described for '%s'\n", label.c_str());
    Exit(0);
  }
}


// #################################################################
// @brief ファイルフォーマットのオプションを指定する
void Control::getFormat_plot3d()
{
  string str;
  string label;
  REAL_TYPE f_val;
  
  // PLOT3Dディレクトリのチェック
  label = "/Output/FormatOption/PLOT3D";
  
  if ( !tpCntl->chkNode(label) )
  {
    Hostonly_ stamped_printf("\tParsing error : Missing '%s'\n", label.c_str());
    Exit(0);
  }
  
  
  
  // インターバル PLOT3D
  label = "/Output/FormatOption/PLOT3D/IntervalType";
  
  if ( !(tpCntl->GetValue(label, &str )) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    Exit(0);
  }
  else
  {
    if     ( !strcasecmp(str.c_str(), "step") )
    {
      Interval[Interval_Manager::tg_plot3d].setMode_Step();
    }
    else if( !strcasecmp(str.c_str(), "time") )
    {
      Interval[Interval_Manager::tg_plot3d].setMode_Time();
    }
    else
    {
      Hostonly_ stamped_printf("\tParsing error : Invalid keyword for '%s'\n", label.c_str());
      Exit(0);
    }
    
    label="/Output/FormatOption/PLOT3D/Interval";
    
    if ( !(tpCntl->GetValue(label, &f_val )) )
    {
      Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
      Exit(0);
    }
    else
    {
      Interval[Interval_Manager::tg_plot3d].setInterval((double)f_val);
    }
  }

}


// #################################################################
// @brief ファイルフォーマットのオプションを指定する．
void Control::getFormat_sph()
{
  string str;
  string label;
  
  
  // sphディレクトリのチェック
  label = "/Output/FormatOption/SPH";
  
  if ( !tpCntl->chkNode(label) )
  {
    Hostonly_ stamped_printf("\tParsing error : Missing '%s'\n", label.c_str());
    Exit(0);
  }
  
  
  // 出力ガイドセルモード
  label = "/Output/FormatOption/SPH/GuideOut";
  
  if ( !(tpCntl->GetValue(label, &str )) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    Exit(0);
	  Exit(0);
  }
  
  if     ( !strcasecmp(str.c_str(), "without") )  GuideOut = 0;
  else if( !strcasecmp(str.c_str(), "with") )     GuideOut = guide;
  else
  {
    Hostonly_ stamped_printf("\tInvalid keyword is described for '%s'\n", label.c_str());
    Exit(0);
  }
  
  // Output Directory_Path
  label = "/Output/FormatOption/SPH/DirectoryPath";
  
  if ( !(tpCntl->GetValue(label, &str)) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    Exit(0);
  }
  // 指定が無ければ，空のまま
  if ( !str.empty() )
  {
    FIO.OutDirPath = str;
  }
  
  
  // TimeSlice option
  label = "/Output/FormatOption/SPH/TimeSlice";
  
  if ( !(tpCntl->GetValue(label, &str)) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    Exit(0);
  }
  
  if ( !strcasecmp(str.c_str(), "on") )
  {
    FIO.Slice = ON;
  }
  else
  {
    FIO.Slice = OFF;
  }
  
  // 1プロセスの場合にはランク番号がないので、タイムスライス毎のディレクトリは作らない
  if ( (Parallelism == Serial) || (Parallelism == OpenMP) )
  {
    FIO.Slice = OFF;
  }
  
  
}


// #################################################################
// 計算モデルの入力ソース情報を取得
void Control::getGeometryModel()
{
  string str;
  string label;
  
  // ソース名を取得　文字列がIntrisicExampleでなければポリゴンファイル名と解釈
  label = "/GeometryModel/Source";
  
  if ( !(tpCntl->GetValue(label, &str )) )
  {
    Hostonly_ stamped_printf("\tParsing error : Invalid char* value in '%s'\n", label.c_str());
    Exit(0);
  }
  
  
  if ( !(tpCntl->GetValue(label, &str )) )
  {
    Hostonly_ stamped_printf("\tError : '%s'\n", label.c_str());
    Exit(0);
  }
  
  if     ( FBUtility::compare(str, "ParallelPlate2D") )   Mode.Example = id_PPLT2D;
  else if( FBUtility::compare(str, "Duct") )              Mode.Example = id_Duct;
  else if( FBUtility::compare(str, "SHC1D") )             Mode.Example = id_SHC1D;
  else if( FBUtility::compare(str, "PerformanceTest") )   Mode.Example = id_PMT;
  else if( FBUtility::compare(str, "Rectangular") )       Mode.Example = id_Rect;
  else if( FBUtility::compare(str, "Cylinder") )          Mode.Example = id_Cylinder;
  else if( FBUtility::compare(str, "BackStep") )          Mode.Example = id_Step;
  else if( FBUtility::compare(str, "Sphere") )            Mode.Example = id_Sphere;
  else if( FBUtility::compare(str, "Jet") )               Mode.Example = id_Jet;
  else
  {
    Mode.Example = id_Polygon;
    PolylibConfigName = str;
  }
  
  
  
  // スケーリングファクター
  REAL_TYPE ct=0.0;
  
  label = "/GeometryModel/ScalingFactor";
  
  if ( !tpCntl->GetValue(label, &ct) )
  {
	  ; // 無くても可
  }
  else
  {
    if ( ct <= 0.0 )
    {
      Hostonly_ stamped_printf("Error : Scaling factor should be positive [%f]\n", ct);
      Exit(0);
    }
    
    Scaling_Factor = ct;
  }
  
}



// #################################################################
// @brief 反復関連の情報を取得する
// @see copyCriteria()
void Control::getIteration()
{
  string str;
  string base, label, leaf;
  int i_val=0;
  double f_val=0.0;
  
  
  base = "/Iteration";
  
  // Iteration
  if( !tpCntl->chkNode(base) )
  {
    Hostonly_ stamped_printf("\tParsing error : Missing the section of '%s'\n", base.c_str());
    Exit(0);
  }
  
  // タグ内のラベル数をチェック
  int nnode = tpCntl->countLabels(base);
  if ( nnode == 0 )
  {
    Hostonly_ stamped_printf("\tNo labels inside /Iteration\n");
    return;
  }

  
  // LinearSolverの個数を取得
  int counter=0;
  for (int i=1; i<=nnode; i++)
  {
    if ( !tpCntl->GetNodeStr(base, i, &str) )
    {
      Hostonly_ stamped_printf("\tGetNodeStr error\n");
      Exit(0);
    }
    
    if( !strcasecmp(str.substr(0,12).c_str(), "LinearSolver") ) counter++;
  }

  
  // Instance of candidates
  NoBaseLS = counter;
  Criteria = new IterationCtl[NoBaseLS];
  
  
  // get criterion
  for (int i=0; i<NoBaseLS; i++)
  {
    if ( !tpCntl->GetNodeStr(base, i+1, &str) )
    {
      Hostonly_ printf("\tParsing error : Missing 'LinearSolver'\n");
      Exit(0);
    }

    if( strcasecmp(str.substr(0,12).c_str(), "LinearSolver") ) continue;
    
    
    // alias 
    leaf = base + "/" + str;
    label = leaf + "/Alias";
    
    if ( !(tpCntl->GetValue(label, &str )) )
    {
      Hostonly_ printf("\tParsing error : No '%s'\n", label.c_str());
      Exit(0);
    }

    // ユニークな名称であること
    for (int n=0; n<i; n++)
    {
      if ( !strcasecmp( Criteria[n].getAlias().c_str(), str.c_str()) )
      {
        Hostonly_ printf("\tParsing error : 'Alias' must be unique\n");
        Exit(0);
      }
    }
    Criteria[i].setAlias(str);
    
    
    // 線形ソルバーの種類
    label = leaf + "/class";
    if ( !(tpCntl->GetValue(label, &str )) )
    {
      Hostonly_ printf("\tParsing error : No '%s'\n", label.c_str());
      Exit(0);
    }
    
    if     ( !strcasecmp(str.c_str(), "SOR") )         Criteria[i].setLS(SOR);
    else if( !strcasecmp(str.c_str(), "SOR2SMA") )     Criteria[i].setLS(SOR2SMA);
    else if( !strcasecmp(str.c_str(), "SOR2CMA") )     Criteria[i].setLS(SOR2CMA);
    else if( !strcasecmp(str.c_str(), "JACOBI") )      Criteria[i].setLS(JACOBI);
    else if( !strcasecmp(str.c_str(), "GMRES") )       Criteria[i].setLS(GMRES);
    else if( !strcasecmp(str.c_str(), "RBGS") )        Criteria[i].setLS(RBGS);
    else if( !strcasecmp(str.c_str(), "PCG") )         Criteria[i].setLS(PCG);
    else if( !strcasecmp(str.c_str(), "PBiCGSTAB") )   Criteria[i].setLS(PBiCGSTAB);
    else if( !strcasecmp(str.c_str(), "VPiteration") ) Criteria[i].setLS(VP_ITERATION);
    else
    {
      Hostonly_ stamped_printf("\tInvalid keyword is described for Linear_Solver\n");
      Exit(0);
    }

    
    label = leaf + "/MaxIteration";
    if ( !(tpCntl->GetValue(label, &i_val )) )
    {
      Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
      Exit(0);
    }
    Criteria[i].setMaxIteration(i_val);


    label = leaf + "/ConvergenceCriterion";
    if ( !(tpCntl->GetValue(label, &f_val )) )
    {
      Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
      Exit(0);
    }
    Criteria[i].setCriterion(f_val);

    
    label = leaf + "/NormType";
    if ( !(tpCntl->GetValue(label, &str )) )
    {
      Hostonly_ printf("\tParsing error : No '%s'\n", label.c_str());
      Exit(0);
    }

    
    int ls = Criteria[i].getLS();
    
    // ノルム
    switch (ls)
    {
      case VP_ITERATION:
        if ( !strcasecmp(str.c_str(), "VdivMax") )
        {
          Criteria[i].setNormType(v_div_max);
        }
        else if ( !strcasecmp(str.c_str(), "VdivDbg") )
        {
          Criteria[i].setNormType(v_div_dbg);
          Mode.Log_Itr == ON;
        }
        else
        {
          Hostonly_ stamped_printf("\tParsing error : Invalid keyword for '%s' for VP iteration\n", str.c_str());
          Exit(0);
        }
        break;
        
      default:
        if ( !strcasecmp(str.c_str(), "DXbyB") )
        {
          Criteria[i].setNormType(dx_b);
        }
        else if ( !strcasecmp(str.c_str(), "RbyB") )
        {
          Criteria[i].setNormType(r_b);
        }
        else if ( !strcasecmp(str.c_str(), "RbyR0") )
        {
          Criteria[i].setNormType(r_r0);
        }
        else
        {
          Hostonly_ stamped_printf("\tParsing error : Invalid keyword for '%s' of Norm for Poisson iteration\n", str.c_str());
          Exit(0);
        }
        break;
    }

    
    // 固有パラメータ
    switch (ls)
    {
      case JACOBI:
        getParaJacobi(leaf, i);
        break;
        
      case SOR:
        getParaSOR(leaf, i);
        break;
        
      case SOR2SMA:
        getParaSOR2(leaf, i);
        break;
        
      case RBGS:
        getParaRBGS(leaf, i);
        break;
        
      case GMRES:
        getParaGmres(leaf, i);
        break;
        
      case PCG:
        getParaPCG(leaf, i);
        break;
        
      case PBiCGSTAB:
        getParaPBiCGSTAB(leaf, i);
        break;
        
      default:
        break;
    }

  }
  
}



// #################################################################
// ログ出力モードを取得
// インターバルパラメータは，setParameters()で無次元して保持
void Control::getLog()
{
  REAL_TYPE f_val=0.0;
  string str;
  string label;
  
  
  // Log_Base
  label="/Output/Log/Base";
  
  if ( !(tpCntl->GetValue(label, &str )) )
  {
	  Hostonly_ stamped_printf("\tParsing error : Invalid string for '%s'\n", label.c_str());
	  Exit(0);
  }
  
  if     ( !strcasecmp(str.c_str(), "on") )   Mode.Log_Base = ON;
  else if( !strcasecmp(str.c_str(), "off") )  Mode.Log_Base = OFF;
  else
  {
    Hostonly_ stamped_printf("\tInvalid keyword is described for '%s'\n", label.c_str());
    Exit(0);
  }
  
  
  // Log_Iteration
  label="/Output/Log/Iteration";
  
  if ( !(tpCntl->GetValue(label, &str )) )
  {
	  Hostonly_ stamped_printf("\tParsing error : Invalid string for '%s'\n", label.c_str());
	  Exit(0);
  }
  
  if     ( !strcasecmp(str.c_str(), "on") )   Mode.Log_Itr = ON;
  else if( !strcasecmp(str.c_str(), "off") )  Mode.Log_Itr = OFF;
  else
  {
    Hostonly_ stamped_printf("\tInvalid keyword is described for '%s'\n", label.c_str());
    Exit(0);
  }
  
  
  // Log_Wall_Info
  if ( Mode.Wall_profile == Log_Law )
  {
	  label="/Output/Log/WallInfo";
    
	  if ( !(tpCntl->GetValue(label, &str )) )
    {
		  Hostonly_ stamped_printf("\tParsing error : Invalid string for '%s'\n", label.c_str());
		  Exit(0);
	  }
	  
	  if     ( !strcasecmp(str.c_str(), "on") )   Mode.Log_Wall = ON;
	  else if( !strcasecmp(str.c_str(), "off") )  Mode.Log_Wall = OFF;
	  else
    {
		  Hostonly_ stamped_printf("\tInvalid keyword is described for '%s'\n", label.c_str());
		  Exit(0);
	  }
  }
  
  
  // Log_Profiling
  label="/Output/Log/Profiling";
  
  if ( !(tpCntl->GetValue(label, &str )) )
  {
	  Hostonly_ stamped_printf("\tParsing error : Invalid string for '%s'\n", label.c_str());
	  Exit(0);
  }
  
  if     ( !strcasecmp(str.c_str(), "on") )     Mode.Profiling = ON;
  else if( !strcasecmp(str.c_str(), "off") )    Mode.Profiling = OFF;
  else if( !strcasecmp(str.c_str(), "detail") ) Mode.Profiling = DETAIL;
  else
  {
    Hostonly_ stamped_printf("\tInvalid keyword is described for '%s'\n", label.c_str());
    Exit(0);
  }
  
  // Interval console
  label="/Output/Log/ConsoleIntervalType";
  
  if ( !(tpCntl->GetValue(label, &str )) )
  {
	  Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
	  Exit(0);
  }
  else
  {
	  if     ( !strcasecmp(str.c_str(), "step") )
    {
      Interval[Interval_Manager::tg_console].setMode_Step();
	  }
	  else if( !strcasecmp(str.c_str(), "time") )
    {
      Interval[Interval_Manager::tg_console].setMode_Time();
	  }
	  else
    {
		  Hostonly_ stamped_printf("\tParsing error : Invalid keyword for '%s'\n", label.c_str());
		  Exit(0);
	  }
	  
	  label="/Output/Log/ConsoleInterval";
    
	  if ( !(tpCntl->GetValue(label, &f_val )) )
    {
		  Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
		  Exit(0);
	  }
	  else
    {
		  Interval[Interval_Manager::tg_console].setInterval((double)f_val);
	  }
  }
  
  // Interval file_history
  label="/Output/Log/HistoryIntervalType";
  
  if ( !(tpCntl->GetValue(label, &str )) )
  {
	  Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
	  Exit(0);
  }
  else
  {
	  if     ( !strcasecmp(str.c_str(), "step") )
    {
		  Interval[Interval_Manager::tg_history].setMode_Step();
	  }
	  else if( !strcasecmp(str.c_str(), "time") )
    {
		  Interval[Interval_Manager::tg_history].setMode_Time();
	  }
	  else
    {
		  Hostonly_ stamped_printf("\tParsing error : Invalid keyword for '%s'\n", label.c_str());
		  Exit(0);
	  }
    
	  label="/Output/Log/HistoryInterval";
    
	  if ( !(tpCntl->GetValue(label, &f_val )) )
    {
		  Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
		  Exit(0);
	  }
	  else
    {
		  Interval[Interval_Manager::tg_history].setInterval((double)f_val);
	  }
  }
}



// #################################################################
// モニタリングのON/OFFとセルモニタの有無のみを取得　詳細パラメータは後ほど
void Control::getMonitorList()
{
  string str;
  string label;
  
  // ログ出力
  label = "/MonitorList/log";
  
  if ( !(tpCntl->GetValue(label, &str )) )
  {
	  Hostonly_ stamped_printf("\tParsing error : Invalid string for '%s'\n", label.c_str());
	  Exit(0);
  }
  
  if     ( !strcasecmp(str.c_str(), "on") )   Sampling.log = ON;
  else if( !strcasecmp(str.c_str(), "off") )  Sampling.log = OFF;
  else
  {
    Hostonly_ stamped_printf("\tInvalid keyword is described for '%s'\n", label.c_str());
    Exit(0);
  }
  
  // セルモニターの指定
  /// @see Initialize.C setExistComponent()
  if ( Sampling.log == ON )
  {
	  label = "/MonitorList/CellMonitor";
    
	  if ( !(tpCntl->GetValue(label, &str )) )
    {
		  Hostonly_ stamped_printf("\tParsing error : Invalid string for '%s'\n", label.c_str());
		  Exit(0);
	  }
    
	  if     ( !strcasecmp(str.c_str(), "on") )   EnsCompo.monitor = ON;
	  else if( !strcasecmp(str.c_str(), "off") )  EnsCompo.monitor = OFF;
	  else
    {
      Hostonly_ stamped_printf("\tInvalid keyword is described for '%s'\n", label.c_str());
		  Exit(0);
	  }
  }
  
}




// #################################################################
// コンポーネント数，媒質数，境界条件数を取得
void Control::getNoOfComponent()
{
  std::string label;
  
  // 媒質の個数
  label = "/MediumTable";
  
  vector<string> nodes_1;
  tpCntl->getLabelVector(label, nodes_1);
  
  NoMedium = nodes_1.size();
  
  if ( NoMedium <= 0)
  {
    Hostonly_ stamped_printf("\tError : Empty MediumTable\n");
    Exit(0);
  }

  
  // 境界条件数
  label = "/BcTable/LocalBoundary";
  
  if ( tpCntl->chkNode(label) )  //nodeがあれば
  {
    vector<string> nodes_2;
    tpCntl->getLabelVector(label, nodes_2);
	  NoBC = nodes_2.size();
  }

  
  NoCompo = NoMedium + NoBC;
}


// #################################################################
/**
 * @brief Gmres反復固有のパラメータを指定する
 * @param [in] base  ラベル
 * @param [in] m     反復制御用クラスの配列番号
 */
void Control::getParaGmres(const string base, const int m)
{
  ;
}


// #################################################################
/**
 * @brief Jacobi反復固有のパラメータを指定する
 * @param [in] base  ラベル
 * @param [in] m     反復制御用クラスの配列番号
 */
void Control::getParaJacobi(const string base, const int m)
{
  string str, label;
  double tmp=0.0; // 加速係数
  
  label = base + "/Omega";
  
  if ( !(tpCntl->GetValue(label, &tmp )) )
  {
    Hostonly_ stamped_printf("\tParsing error : Invalid float value for '%s'\n", label.c_str());
    Exit(0);
  }
  Criteria[m].setOmega(tmp);
  
}


// #################################################################
/**
 * @brief PCG反復固有のパラメータを指定する
 * @param [in] base  ラベル
 * @param [in] m     反復制御用クラスの配列番号
 */
void Control::getParaPCG(const string base, const int m)
{
  ;
}


// #################################################################
/**
 * @brief PBiCGSTAB反復固有のパラメータを指定する
 * @param [in] base  ラベル
 * @param [in] m     反復制御用クラスの配列番号
 */
void Control::getParaPBiCGSTAB(const string base, const int m)
{
  ;
}


// #################################################################
/**
 * @brief RBGS反復固有のパラメータを指定する
 * @param [in] base  ラベル
 * @param [in] m     反復制御用クラスの配列番号
 */
void Control::getParaRBGS(const string base, const int m)
{
  ;
}


// #################################################################
/**
 * @brief SOR反復固有のパラメータを指定する
 * @param [in] base  ラベル
 * @param [in] m     反復制御用クラスの配列番号
 */
void Control::getParaSOR(const string base, const int m)
{
  getParaJacobi(base, m);
}


// #################################################################
/**
 * @brief RB-SOR反復固有のパラメータを指定する
 * @param [in] base  ラベル
 * @param [in] m     反復制御用クラスの配列番号
 */
void Control::getParaSOR2(const string base, const int m)
{
  string str, label;
  
  getParaJacobi(base, m);
  
  label = base + "/commMode";
  if ( !(tpCntl->GetValue(label, &str )) )
  {
    Hostonly_ stamped_printf("\tParsing error : Invalid char* value for '%s'\n", label.c_str());
    Exit(0);
  }
  if ( !strcasecmp(str.c_str(), "sync") )
  {
    Criteria[m].setSyncMode(comm_sync);
  }
  else if ( !strcasecmp(str.c_str(), "async") )
  {
    Criteria[m].setSyncMode(comm_async);
  }
  else
  {
    Hostonly_ stamped_printf("\tParsing error : Invalid char* value for '%s'\n", label.c_str());
    Exit(0);
  }
  
}



// #################################################################
/**
 * @brief 参照パラメータを取得
 * @note Ref_IDで指定される媒質を代表物性値とする
 */
void Control::getReference()
{
  REAL_TYPE ct2;
  string label, str;
  
  
  label = "/Reference/Length";
  if ( !(tpCntl->GetValue(label, &ct2 )) )
  {
    Hostonly_ printf("\tParsing error in '%s'\n", label.c_str());
    Exit(0);
  }
  RefLength = ct2;
  
  
  label = "/Reference/Velocity";
  if ( !(tpCntl->GetValue(label, &ct2 )) )
  {
    Hostonly_ printf("\tParsing error in '%s'\n", label.c_str());
    Exit(0);
  }
  RefVelocity = ct2;
  
  
  label = "/Reference/MassDensity";
  if ( !(tpCntl->GetValue(label, &ct2 )) )
  {
    Hostonly_ printf("\tParsing error in '%s'\n", label.c_str());
    Exit(0);
  }
  RefDensity = ct2;
  
  
  label = "/Reference/BasePressure";
  if ( !(tpCntl->GetValue(label, &ct2 )) )
  {
    Hostonly_ printf("\tParsing error in '%s'\n", label.c_str());
    Exit(0);
  }
  BasePrs = ct2;

  
  label = "/Reference/Medium";
  if ( !tpCntl->GetValue(label, &str) )
  {
    Hostonly_ printf("\tParsing error in '%s'\n", label.c_str());
	  Exit(0);
  }
  RefMedium = str;
  
  
  if ( isHeatProblem() )
  {
    REAL_TYPE Base, Diff;
    
    label="/Reference/Temperature/Base";
    
    if ( !(tpCntl->GetValue(label, &Base )) )
    {
      Hostonly_ stamped_printf("\tParsing error : Invalid float value for '%s'\n", label.c_str());
      Exit(0);
    }
    
    
    label="/Reference/Temperature/Difference";
    
    if ( !(tpCntl->GetValue(label, &Diff )) )
    {
      Hostonly_ stamped_printf("\tParsing error : Invalid float value for '%s'\n", label.c_str());
      Exit(0);
    }
    
    
    if ( Diff < 0.0f )
    {
      Hostonly_ stamped_printf("\tTemperature difference must be positive.\n");
      Exit(0);
    }
    DiffTemp = Diff;
    
    if ( Unit.Temp == Unit_CELSIUS )
    {
      BaseTemp = Base + KELVIN;
    }
  }
  
}



// #################################################################
/**
 * @brief 参照座標系を取得する
 * @todo 回転は未
 */
void Control::getReferenceFrame(ReferenceFrame* RF)
{
  string str;
  string label;
  
  label="/ReferenceFrame/Mode";
  
  if ( !(tpCntl->GetValue(label, &str )) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
	  Exit(0);
  }
  
  if     ( !strcasecmp(str.c_str(), "stationary") )
  {
    RF->setFrame(ReferenceFrame::frm_static);
  }
  else if( !strcasecmp(str.c_str(), "translational") )
  {
    RF->setFrame(ReferenceFrame::frm_translation);
    
    REAL_TYPE xyz[3];
    for (int n=0; n<3; n++) xyz[n]=0.0;
    
    if( tpCntl->GetVector(label, xyz, 3) )
    {
      Hostonly_ stamped_printf("\tParsing error : Invalid values for '%s'\n", label.c_str());
      Exit(0);
    }
    RF->setGridVel((double*)xyz);
  }
  else
  {
    Hostonly_ stamped_printf("\tParsing error : Invalid keyword for '%s'\n", label.c_str());
    Exit(0);
  }
}



// #################################################################
// @brief ソルバーの種類を特定するパラメータを取得する
// @note ガイドセルの値に影響
void Control::getShapeApproximation()
{
  string str;
  string label;
  
  // 形状近似度の取得
  label = "/ShapeApproximation/Method";
  
  if ( !(tpCntl->GetValue(label, &str )) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    Exit(0);
  }
  
  if     ( !strcasecmp(str.c_str(), "Binary") )      Mode.ShapeAprx = BINARY;
  else if( !strcasecmp(str.c_str(), "CutDistance") ) Mode.ShapeAprx = CUT_INFO;
  else
  {
    Hostonly_ stamped_printf("\tInvalid keyword is described for '%s'\n", label.c_str());
    Exit(0);
  }
}


// #################################################################
// ソルバーの種類を特定するパラメータを取得し，ガイドセルの値を決定する
void Control::getSolverProperties()
{
  string str;
  string label;
  
  
  // 支配方程式の型（PDE_NS / Euler）を取得
  label = "/GoverningEquation/PDEType";
  
  if ( !(tpCntl->GetValue(label, &str )) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    Exit(0);
  }
  
  if     ( !strcasecmp(str.c_str(), "NavierStokes" ) ) Mode.PDE = PDE_NS;
  else if( !strcasecmp(str.c_str(), "Euler" ) )        Mode.PDE = PDE_EULER;
  else
  {
    Hostonly_ stamped_printf("\tInvalid keyword is described for '%s'\n", label.c_str());
    Exit(0);
  }
  
  
  // 基礎方程式の種類を取得する
  label = "/GoverningEquation/FlowEquation";
  
  if ( !(tpCntl->GetValue(label, &str )) ) {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    Exit(0);
  }
  
  if     ( !strcasecmp(str.c_str(), "Incompressible" ) )        BasicEqs = INCMP;
  else if( !strcasecmp(str.c_str(), "LimitedCompressible" ) )   BasicEqs = LTDCMP;
  else if( !strcasecmp(str.c_str(), "Compressible" ) )          BasicEqs = CMPRSS;
  else if( !strcasecmp(str.c_str(), "Incompressible2Phase" ) )  BasicEqs = INCMP_2PHASE;
  else
  {
    Hostonly_ stamped_printf("\tInvalid keyword is described for '%s'\n", label.c_str());
    Exit(0);
  }
  
  
  // 非定常計算，または定常計算の種別を取得する
  label = "/GoverningEquation/TimeVariation";
  
  if ( !(tpCntl->GetValue(label, &str )) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    Exit(0);
  }
  
  if     ( !strcasecmp(str.c_str(), "Steady" ) )    Mode.Steady = FB_STEADY;
  else if( !strcasecmp(str.c_str(), "Unsteady" ) )  Mode.Steady = FB_UNSTEADY;
  else
  {
    Hostonly_ stamped_printf("\tInvalid keyword is described for '%s'\n", label.c_str());
    Exit(0);
  }
  
  
  // 対流項スキームの種類の取得
  int ct;
  
  label="/ConvectionTerm/scheme";
  
  if ( !(tpCntl->GetValue(label, &str )) )
  {
	  Hostonly_ stamped_printf("\tParsing error : Invalid char* value for '%s'\n", label.c_str());
	  Exit(0);
  }
  
  if     ( !strcasecmp(str.c_str(), "O1Upwind") )    CnvScheme = O1_upwind;
  else if( !strcasecmp(str.c_str(), "O3muscl") )     CnvScheme = O3_muscl;
  else if( !strcasecmp(str.c_str(), "O2central") )   CnvScheme = O2_central;
  else if( !strcasecmp(str.c_str(), "O4central") ) { CnvScheme = O4_central; Exit(0); }  // not yet implemented
  else
  {
    Hostonly_ stamped_printf("\tInvalid keyword is described for '%s'\n", label.c_str());
    Exit(0);
  }
  
  // Limiter
  if ( CnvScheme == O3_muscl )
  {
		label="/ConvectionTerm/limiter";
    
		if ( !(tpCntl->GetValue(label, &str )) )
    {
			Hostonly_ stamped_printf("\tParsing error : Invalid char* value for '%s'\n", label.c_str());
			Exit(0);
		}
    
		if     ( !strcasecmp(str.c_str(), "NoLimiter") )  Limiter = No_Limiter;
		else if( !strcasecmp(str.c_str(), "Minmod") )     Limiter = MinMod;
		else
    {
			Hostonly_ stamped_printf("\tInvalid keyword is described for '%s'\n", label.c_str());
			Exit(0);
		}
  }
  
  
  // ソルバーの種類（FLOW_ONLY / THERMAL_FLOW / THERMAL_FLOW_NATURAL / CONJUGATE_HEAT_TRANSFER / SOLID_CONDUCTION）と浮力モード
  label="/GoverningEquation/HeatEquation";
  
  if ( !(tpCntl->GetValue(label, &str )) )
  {
	  Hostonly_ stamped_printf("\tParsing error : Invalid char* value for '%s'\n", label.c_str());
	  Exit(0);
  }
  
  if     ( !strcasecmp(str.c_str(), "FlowOnly" ) )               KindOfSolver = FLOW_ONLY;
  else if( !strcasecmp(str.c_str(), "ThermalFlow" ) )            KindOfSolver = THERMAL_FLOW;
  else if( !strcasecmp(str.c_str(), "ThermalFlowNatural" ) )     KindOfSolver = THERMAL_FLOW_NATURAL;
  else if( !strcasecmp(str.c_str(), "ConjugateHeatTransfer" ) )  KindOfSolver = CONJUGATE_HEAT_TRANSFER;
  else if( !strcasecmp(str.c_str(), "SolidConduction" ) )        KindOfSolver = SOLID_CONDUCTION;
  else
  {
    Hostonly_ stamped_printf("\tInvalid keyword is described for '%s'\n", label.c_str());
    Exit(0);
  }
  
  // Buoyancy option
  if ( (KindOfSolver==THERMAL_FLOW) || (KindOfSolver==THERMAL_FLOW_NATURAL) || (KindOfSolver==CONJUGATE_HEAT_TRANSFER) )
  {
    label="/GoverningEquation/Buoyancy";
    
    if ( !(tpCntl->GetValue(label, &str )) )
    {
      Hostonly_ stamped_printf("\tParsing error : Invalid char* value for '%s'\n", label.c_str());
      Exit(0);
    }
    
    if     ( !strcasecmp(str.c_str(), "Boussinesq" ) )  Mode.Buoyancy = BOUSSINESQ;
    else if( !strcasecmp(str.c_str(), "LowMach" ) )     Mode.Buoyancy = LOW_MACH;
    else if( !strcasecmp(str.c_str(), "NoBuoyancy" ) )  Mode.Buoyancy = NO_BUOYANCY;
    else
    {
      Hostonly_ stamped_printf("\tInvalid keyword is described for '/GovernignEquation/Buoyancy'\n");
      Exit(0);
    }
  }
  
  
  // ガイドセルの値を決める
  if (KindOfSolver==SOLID_CONDUCTION)
  {
    guide = 1;
  }
  else
  {
    switch (CnvScheme)
    {
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



// #################################################################
// @brief 初期値とリスタート条件
// @todo セルフェイスの粗格子リスタート  >> 近似なのでサボる？
// @ see getTimeControl()
void Control::getStartCondition()
{
  int ct;
  string str;
  string label;
  
  
  // Staging option
  label="/StartCondition/Restart/Staging";
  
  if ( !(tpCntl->GetValue(label, &str )) )
  {
    ;
  }
  else
  {
    if     ( !strcasecmp(str.c_str(), "on") )  Restart_staging = ON;
    else if( !strcasecmp(str.c_str(), "off") ) Restart_staging = OFF;
    else
    {
      Hostonly_ stamped_printf("\tInvalid keyword is described for '%s'\n", label.c_str());
      Exit(0);
    }
  }

  
  // リスタート時のDFIファイル名
  if ( Start != initial_start )
  {
    label="/StartCondition/Restart/DFIfiles/Pressure";
    
    if ( tpCntl->GetValue(label, &str ) )
    {
      f_dfi_in_prs = str.c_str();
    }
    if ( f_dfi_in_prs.empty() == true ) f_dfi_in_prs = "prs";
    
    
    label="/StartCondition/Restart/DFIfiles/Velocity";
    
    if ( tpCntl->GetValue(label, &str ) )
    {
      f_dfi_in_vel = str.c_str();
    }
    if ( f_dfi_in_vel.empty() == true ) f_dfi_in_vel = "vel";
    
    
    if ( Mode.FaceV == ON )
    {
      label="/StartCondition/Restart/DFIfiles/Fvelocity";
      
      if ( tpCntl->GetValue(label, &str ) )
      {
        f_dfi_in_fvel = str.c_str();
      }
      if ( f_dfi_in_fvel.empty() == true ) f_dfi_in_fvel = "fvel";
    }
    
    
    if ( isHeatProblem() )
    {
      label="/StartCondition/Restart/DFIfiles/Temperature";
      
      if ( tpCntl->GetValue(label, &str ) )
      {
        f_dfi_in_temp = str.c_str();
      }
      if ( f_dfi_in_temp.empty() == true ) f_dfi_in_temp = "tmp";
    }
    
    
    // 平均値
    if ( Mode.Average == ON )
    {
      label="/StartCondition/Restart/DFIfiles/AveragedPressure";
      
      if ( tpCntl->GetValue(label, &str ) )
      {
        f_dfi_in_prsa = str.c_str();
      }
      if ( f_dfi_in_prsa.empty() == true ) f_dfi_in_prsa = "prsa";
      
      
      label="/StartCondition/Restart/DFIfiles/AveragedVelocity";
      
      if ( tpCntl->GetValue(label, &str ) )
      {
        f_dfi_in_vela = str.c_str();
      }
      if ( f_dfi_in_vela.empty() == true ) f_dfi_in_vela = "vela";
      
      
      if ( isHeatProblem() )
      {
        label="/StartCondition/Restart/DFIfiles/AveragedTemperature";
        
        if ( tpCntl->GetValue(label, &str ) )
        {
          f_dfi_in_tempa = str.c_str();
        }
        if ( f_dfi_in_tempa.empty() == true ) f_dfi_in_tempa = "tmpa";
      }
    }
  }
  
  

  // 初期条件
  if ( Start == initial_start )
  {
    // Density
    label="/StartCondition/InitialState/MassDensity";
    
    if ( !(tpCntl->GetValue(label, &iv.Density )) )
    {
      Hostonly_ stamped_printf("\tParsing error : Invalid float value for '%s'\n", label.c_str());
      Exit(0);
    }
    
    // Pressure
    label="/StartCondition/InitialState/Pressure";
    
    if ( !(tpCntl->GetValue(label, &iv.Pressure )) )
    {
      Hostonly_ stamped_printf("\tParsing error : Invalid float value for '%s'\n", label.c_str());
      Exit(0);
    }
    
    // Velocity
    REAL_TYPE v[3];
    for (int n=0; n<3; n++) v[n]=0.0;
    label="/StartCondition/InitialState/Velocity";
    
    if( !(tpCntl->GetVector(label, v, 3)) )
    {
      Hostonly_ stamped_printf("\tParsing error : fail to get velocity in '%s'\n", label.c_str());
      Exit(0);
    }
    iv.VecU = v[0];
    iv.VecV = v[1];
    iv.VecW = v[2];
    
    // Temperature
    if ( isHeatProblem() )
    {
      label="/StartCondition/InitialState/Temperature";
      
      if ( !(tpCntl->GetValue(label, &iv.Temperature )) )
      {
        Hostonly_ stamped_printf("\tParsing error : Invalid float value for '%s'\n", label.c_str());
        Exit(0);
      }
      
    }
  }
  
}


// #################################################################
// 解法アルゴリズムを選択する
void Control::getSolvingMethod()
{
  string str;
  string label;
  
  label = "/SolvingMethod/Flow";
  
  if ( !(tpCntl->GetValue(label, &str )) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    Exit(0);
  }
  
  if     ( !strcasecmp(str.c_str(), "FS_C_EE_D_EE") )     AlgorithmF = Flow_FS_EE_EE;
  else if( !strcasecmp(str.c_str(), "FS_C_RK_D_CN") )     AlgorithmF = Flow_FS_RK_CN;
  else if( !strcasecmp(str.c_str(), "FS_C_AB_D_AB") )     AlgorithmF = Flow_FS_AB2;
  else if( !strcasecmp(str.c_str(), "FS_C_AB_D_CN") )     AlgorithmF = Flow_FS_AB_CN;
  else
  {
    Hostonly_ stamped_printf("\tParsing error : Invalid keyword for '%s'\n", label.c_str());
    Exit(0);
  }
  
  // Heat
  if ( isHeatProblem() )
  {
	  label = "/SolvingMethod/Heat";
    
	  if ( !(tpCntl->GetValue(label, &str )) )
    {
		  Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
		  Exit(0);
	  }
    
	  if     ( !strcasecmp(str.c_str(), "C_EE_D_EE") )    AlgorithmH = Heat_EE_EE;
	  else if( !strcasecmp(str.c_str(), "C_EE_D_EI") )    AlgorithmH = Heat_EE_EI;
	  else
    {
		  Hostonly_ stamped_printf("\tParsing error : Invalid keyword for '%s'\n", label.c_str());
		  Exit(0);
	  }
  }
}



// #################################################################
/**
 * @brief 時間制御に関するパラメータを取得する
 * @param [out] DT  DTcntl
 * @note パラメータは，setParameters()で無次元して保持
 */
void Control::getTimeControl(DTcntl* DT)
{
  double ct;
  int ss=0;
  
  string str;
  string label;
  
  // 加速時間
  label = "/TimeControl/AccelerationType";
  
  if ( !(tpCntl->GetValue(label, &str )) )
  {
	  Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
	  Exit(0);
  }
  else
  {
	  if     ( !strcasecmp(str.c_str(), "step") )
    {
		  Interval[Interval_Manager::tg_accelra].setMode_Step();
	  }
	  else if( !strcasecmp(str.c_str(), "time") )
    {
		  Interval[Interval_Manager::tg_accelra].setMode_Time();
	  }
	  else
    {
		  Hostonly_ stamped_printf("\tParsing error : Invalid keyword for '%s'\n", label.c_str());
		  Exit(0);
	  }
	  
	  label = "/TimeControl/Acceleration";
    
	  if ( !(tpCntl->GetValue(label, &ct )) )
    {
		  Hostonly_ stamped_printf("\tParsing error : fail to get '%s\n", label.c_str());
		  Exit(0);
	  }
	  else
    {
		  Interval[Interval_Manager::tg_accelra].setInterval(ct);
	  }
  }
  
  
  // 時間積分幅を取得する
  label = "/TimeControl/DeltaTType";
  
  if ( !(tpCntl->GetValue(label, &str )) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    Exit(0);
  }
  
  
  label = "/TimeControl/DeltaT";
  
  if ( !(tpCntl->GetValue(label, &ct )) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    Exit(0);
  }
  
  
  // Directで有次元の場合は，無次元化
  double ts = (double)RefLength / (double)RefVelocity;
  double cc;
  if ( !strcasecmp(str.c_str(), "Direct") )
  {
    if (Unit.Param == DIMENSIONAL)
    {
      cc = ct / ts;
    }
  }
  else
  {
    cc = ct;
  }
  
  if ( !DT->set_Scheme(str.c_str(), cc) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to set DeltaT\n");
    Exit(0);
  }
  
  // 計算する時間を取得する
  label = "/TimeControl/TemporalType";
  
  if ( !(tpCntl->GetValue(label, &str )) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    Exit(0);
  }
  else
  {
    if ( !strcasecmp(str.c_str(), "step") )
    {
      Interval[Interval_Manager::tg_compute].setMode_Step();
    }
    else if ( !strcasecmp(str.c_str(), "time") )
    {
      Interval[Interval_Manager::tg_compute].setMode_Time();
    }
    else
    {
      Hostonly_ stamped_printf("\tParsing error : Invalid keyword for '%s'\n", label.c_str());
      Exit(0);
    }
  }
  
  // スタート
  label = "/TimeControl/Start";
  
  if ( !(tpCntl->GetValue(label, &ct )) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    Exit(0);
  }
  double m_start = ct;
  
  // 終了
  label = "/TimeControl/End";
  
  if ( !(tpCntl->GetValue(label, &ct )) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    Exit(0);
  }
  double m_end = ct;
  
  // チェック
  if ( m_start >= m_end )
  {
    Hostonly_ stamped_printf("\tError : Start msut be less than End.\n");
    Exit(0);
  }
  
  // 丸め誤差の範囲でゼロ
  if ( fabs(m_start)< ROUND_EPS )
  {
    Start = initial_start;
  }
  else
  {
    Start = restart_sameDiv_sameRes; // リスタートタイプのデフォルト
    Restart_step = m_start;
  }
  
  Interval[Interval_Manager::tg_compute].setInterval(m_end-m_start);
  
  
  
  // 平均操作開始
  label = "/TimeControl/Average/Start";
  
  if ( !(tpCntl->GetValue(label, &ct )) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    Exit(0);
  }
  double avr_start = ct;
  
  Interval[Interval_Manager::tg_average].setStart(avr_start, avr_start);
  Restart_stepAvr = avr_start;
  
  
  // 平均操作終了
  label = "/TimeControl/Average/End";
  
  if ( !(tpCntl->GetValue(label, &ct )) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    Exit(0);
  }
  double avr_end = ct;
  
  // チェック
  if ( avr_start > avr_end )
  {
    Hostonly_ stamped_printf("\tError : Average/Start msut be less than Average/End.\n");
    Exit(0);
  }
  
  
  // 平均値操作の判断
  if ( Start == initial_start )
  {
    Mode.AverageRestart = OFF;
    
    if ( (avr_end > 0.0) && (avr_start >= 0.0) )
    {
      Mode.Average = ON;
    }
    else
    {
      Mode.Average = OFF;
    }
  }
  else
  {
    if ( (avr_start < m_start) && (avr_end < m_start) )
    {
      Mode.Average = OFF;
      Mode.AverageRestart = OFF;
    }
    else if ( (avr_start < m_start) && (avr_end > m_start) )
    {
      Mode.Average = ON;
      Mode.AverageRestart = ON;
    }
    else if ( (avr_start > m_start) && (avr_end > m_start) )
    {
      Mode.Average = ON;
      Mode.AverageRestart = OFF;
    }
  }

}



// #################################################################
// @brief 乱流計算のオプションを取得する
void Control::getTurbulenceModel()
{
  REAL_TYPE ct;
  string str;
  string label;
  
  
  LES.Calc = OFF;
  
  // モデル
  label = "/TurbulenceModeling/Model";
  
  if ( !(tpCntl->GetValue(label, &str )) )
  {
	  Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
	  Exit(0);
  }
  
  if      ( !strcasecmp(str.c_str(), "no") )
  {
    LES.Calc = OFF;
  }
  else if ( !strcasecmp(str.c_str(), "smagorinsky") )
  {
    LES.Calc = ON;
    LES.Model = Smagorinsky;
  }
  else if ( !strcasecmp(str.c_str(), "LowReynolds") )
  {
    LES.Calc = ON;
    LES.Model = Low_Reynolds;
  }
  else if ( !strcasecmp(str.c_str(), "Dynamic") )
  {
    LES.Calc = ON;
    LES.Model = Dynamic;
  }
  else
  {
    Hostonly_ stamped_printf("\tParsing error : Invalid keyword for '%s'\n", label.c_str());
    Exit(0);
  }
  
  // スマゴリンスキーモデル
  if ( LES.Model == Smagorinsky )
  {
    // Cs係数
    label = "/TurbulenceModeling/Cs";
    if ( !(tpCntl->GetValue(label, &ct )) )
    {
      Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
      Exit(0);
    }
    LES.Cs = ct;
    
    // damping factor
    label="/TurbulenceModeling/DampingFactor";
    if ( !(tpCntl->GetValue(label, &ct )) )
    {
      Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
      Exit(0);
    }
    LES.damping_factor = ct;
  }
  
}



// #################################################################
// 入力ファイルに記述するパラメータとファイルの有次元・無次元の指定を取得する
void Control::getUnit()
{
  string str;
  string label;
  
  label = "/Unit/UnitOfInputParameter";
  
  if ( !(tpCntl->GetValue(label, &str )) )
  {
		Hostonly_ stamped_printf("\tParsing error : Invalid string for '%s'\n", label.c_str());
	  Exit(0);
  }
  
  if     ( !strcasecmp(str.c_str(), "Dimensional") )     Unit.Param = DIMENSIONAL;
  else if( !strcasecmp(str.c_str(), "NonDimensional") )  Unit.Param = NONDIMENSIONAL;
  else
  {
    Hostonly_ stamped_printf("\tInvalid keyword is described at '%s'\n", label.c_str());
    Exit(0);
  }
  
  
  label = "/Unit/UnitOfOutput";
  
  if ( !(tpCntl->GetValue(label, &str )) )
  {
		Hostonly_ stamped_printf("\tParsing error : Invalid string for '%s'\n", label.c_str());
	  Exit(0);
  }
  
  if     ( !strcasecmp(str.c_str(), "Dimensional") )
  {
    Unit.Output = DIMENSIONAL;
    Unit.File   = DIMENSIONAL;
    Unit.Log    = DIMENSIONAL;
  }
  else if( !strcasecmp(str.c_str(), "NonDimensional") )
  {
    Unit.Output = NONDIMENSIONAL;
    Unit.File   = NONDIMENSIONAL;
    Unit.Log    = NONDIMENSIONAL;
  }
  else
  {
    Hostonly_ stamped_printf("\tInvalid keyword is described at '%s'\n", label.c_str());
    Exit(0);
  }
  
  
  label = "/Unit/Pressure";
  
  if ( !(tpCntl->GetValue(label, &str )) )
  {
    Hostonly_ stamped_printf("\tParsing error : Invalid string for '%s'\n", label.c_str());
	  Exit(0);
  }
  
  if     ( !strcasecmp(str.c_str(), "Gauge") )    Unit.Prs = Unit_Gauge;
  else if( !strcasecmp(str.c_str(), "Absolute") ) Unit.Prs = Unit_Absolute;
  else
  {
    Hostonly_ stamped_printf("\tInvalid keyword is described at '%s'\n", label.c_str());
    Exit(0);
  }
  
  if ( isHeatProblem() )
  {
    label = "/Unit/Temperature";
    
    if ( !(tpCntl->GetValue(label, &str )) )
    {
      Hostonly_ stamped_printf("\tParsing error : Invalid string for '%s'\n", label.c_str());
      Exit(0);
    }
    
    if     ( !strcasecmp(str.c_str(), "Celsius") )  Unit.Temp = Unit_CELSIUS;
    else if( !strcasecmp(str.c_str(), "Kelvin") )   Unit.Temp = Unit_KELVIN;
    else
    {
      Hostonly_ stamped_printf("\tInvalid keyword is described at '%s'\n", label.c_str());
      Exit(0);
    }
  }
  
}



// #################################################################
//壁面上の扱いを指定する
void Control::getWallType()
{
  string str;
  string label;
  
  // 圧力のタイプ
  label="/TreatmentOfWall/PressureGradient";
  
  if ( !(tpCntl->GetValue(label, &str )) )
  {
	  Hostonly_ stamped_printf("\tParsing error : Invalid string for '%s'\n", label.c_str());
	  Exit(0);
  }
  
  if     ( !strcasecmp(str.c_str(), "GradZero") )   Mode.PrsNeuamnnType = P_GRAD_ZERO;
  else if( !strcasecmp(str.c_str(), "GradNS") )     Mode.PrsNeuamnnType = P_GRAD_NS;
  else
  {
    Hostonly_ stamped_printf("\tInvalid keyword is described for '%s'\n", label.c_str());
    Exit(0);
  }
  
  // 壁面摩擦応力の計算モード
  label="/TreatmentOfWall/VelocityProfile";
  
  if ( !(tpCntl->GetValue(label, &str )) )
  {
	  Hostonly_ stamped_printf("\tParsing error : Invalid string for '%s'\n", label.c_str());
	  Exit(0);
  }
  
  if     ( !strcasecmp(str.c_str(), "NoSlip") )    Mode.Wall_profile = No_Slip;
  else if( !strcasecmp(str.c_str(), "Slip") )      Mode.Wall_profile = Slip;
  else if( !strcasecmp(str.c_str(), "LawOfWall") ) Mode.Wall_profile = Log_Law;
  else
  {
    Hostonly_ stamped_printf("\tInvalid keyword is described for '%s'\n", label.c_str());
    Exit(0);
  }
}



// #################################################################
// TPのポインタを受け取る
void Control::importTP(TPControl* tp) 
{ 
  if ( !tp ) Exit(0);
  tpCntl = tp;
}


// #################################################################
// 外部境界の各方向の開口率（流体部分の比率）
REAL_TYPE Control::OpenDomainRatio(const int dir, const REAL_TYPE area, const int* G_size)
{
  REAL_TYPE r = 0.0;
  
  int m_imax = G_size[0];
  int m_jmax = G_size[1];
  int m_kmax = G_size[2];
  
  switch (dir) 
  {
    case X_MINUS:
    case X_PLUS:
      r = area / (REAL_TYPE)(m_jmax*m_kmax) * 100.0;
      break;
      
    case Y_MINUS:
    case Y_PLUS:
      r = area / (REAL_TYPE)(m_imax*m_kmax) * 100.0;
      break;
      
    case Z_MINUS:
    case Z_PLUS:
      r = area / (REAL_TYPE)(m_imax*m_jmax) * 100.0;
      break;
  }
  
  return r;
}



// #################################################################
// 全計算領域の有効セル数と外部境界面の開口率を表示する
void Control::printOuterArea(FILE* fp, unsigned long G_Fcell, unsigned long G_Acell, int* G_size)
{
  if( !fp ) {
    stamped_printf("\tFail to write into file\n");
    Exit(0);
  }
  
  REAL_TYPE cell_max = getCellSize(G_size);
  
  fprintf(fp,"\n---------------------------------------------------------------------------\n\n");
  fprintf(fp,"\n\t>> Effective cells and Open Area of Computational Domain\n\n");
  
  fprintf(fp,"\tFluid  cell inside whole Computational domain = %15ld (%8.4f percent)\n", G_Fcell, (REAL_TYPE)G_Fcell/cell_max *100.0);
  fprintf(fp,"\tActive cell                                   = %15ld (%8.4f percent)\n", G_Acell, (REAL_TYPE)G_Acell/cell_max *100.0);
  
  fprintf(fp,"\n\tFace :      Element (Open ratio)\n");
  for (int i=0; i<NOFACE; i++) {
    fprintf(fp,"\t  %s : %12.0f (%6.2f percent)\n", FBUtility::getDirection(i).c_str(), OpenDomain[i], OpenDomainRatio(i, OpenDomain[i], G_size));
  }
  fprintf(fp,"\n");
  fflush(fp);
}


// #################################################################
// グローバルな領域情報を表示する
void Control::printGlobalDomain(FILE* fp, const int* G_size, const REAL_TYPE* G_org, const REAL_TYPE* G_reg, const REAL_TYPE* pch)
{
  REAL_TYPE PB=0.0, TB=0.0, GB=0.0, MB=0.0, KB=0.0, total=0.0;
  KB = 1000.0;
  MB = 1000.0*KB;
  GB = 1000.0*MB;
  TB = 1000.0*GB;
  PB = 1000.0*TB;
  
  fprintf(fp,"\timax, jmax, kmax    = %13d %13d %13d     >> ", 
          G_size[0], 
          G_size[1], 
          G_size[2]);

  total = (REAL_TYPE)G_size[0] * (REAL_TYPE)G_size[1] * (REAL_TYPE)G_size[2];
  
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
  fprintf(fp,"\n");
  
  fprintf(fp,"\t(dx, dy, dz)  [m] / [-] = (%13.6e %13.6e %13.6e)  /  (%13.6e %13.6e %13.6e)\n",    
          pch[0]*RefLength,
          pch[1]*RefLength,
          pch[2]*RefLength,
          pch[0], 
          pch[1], 
          pch[2]);
  
  fprintf(fp,"\t(ox, oy, oz)  [m] / [-] = (%13.6e %13.6e %13.6e)  /  (%13.6e %13.6e %13.6e)\n", 
          G_org[0]*RefLength, 
          G_org[1]*RefLength, 
          G_org[2]*RefLength, 
          G_org[0], 
          G_org[1], 
          G_org[2]);
  
  fprintf(fp,"\t(Lx, Ly, Lz)  [m] / [-] = (%13.6e %13.6e %13.6e)  /  (%13.6e %13.6e %13.6e)\n", 
          G_reg[0]*RefLength, 
          G_reg[1]*RefLength, 
          G_reg[2]*RefLength, 
          G_reg[0], 
          G_reg[1], 
          G_reg[2]);
  fprintf(fp,"\n");
  
  fflush(fp);
}


// #################################################################
/**
 * @brief 初期値の表示
 * @note Init*には無次元値が保持されている
 * @see void Control::setInitialConditions(void)
 */
void Control::printInitValues(FILE* fp)
{
  if ( !fp )
  {
    stamped_printf("\tFail to write into file\n");
    Exit(0);
  }
  
  // ここでは動圧を表示、方程式の無次元化は ( /rho u^2 ) 
  REAL_TYPE DynamicPrs = 0.5*RefDensity * RefVelocity * RefVelocity;
  
  fprintf(fp,"\n---------------------------------------------------------------------------\n\n");
  fprintf(fp,"\n\t>> Initial Values for Physical Variables\n\n");
  
  fprintf(fp,"\tInitial  MassDensity [kg/m^3]/ [-]   : %12.5e / %12.5e\n", iv.Density, iv.Density/RefDensity);
  fprintf(fp,"\tInitial  Velocity.U  [m/s]   / [-]   : %12.5e / %12.5e\n", iv.VecU,    iv.VecU/RefVelocity);
  fprintf(fp,"\tInitial          .V  [m/s]   / [-]   : %12.5e / %12.5e\n", iv.VecV,    iv.VecV/RefVelocity);
  fprintf(fp,"\tInitial          .W  [m/s]   / [-]   : %12.5e / %12.5e\n", iv.VecW,    iv.VecW/RefVelocity);
  fprintf(fp,"\tDynamic  Pressure    [Pa]    / [-]   : %12.5e / %12.5e\n", DynamicPrs,             1.0);
  
  if (Unit.Prs == Unit_Absolute)
  {
    fprintf(fp,"\tInitial  Pressure    [Pa]    / [-]   : %12.5e / %12.5e\n", iv.Pressure, (iv.Pressure-BasePrs)/DynamicPrs);
  }
  else
  {
    fprintf(fp,"\tInitial  Pressure    [Pa_g]  / [-]   : %12.5e / %12.5e\n", iv.Pressure, iv.Pressure/DynamicPrs);
  }
  
  if ( isHeatProblem() )
  {
    fprintf(fp,"\tInitial  Temperature [%s]     / [-]   : %12.5e / %12.5e\n", (Unit.Temp==Unit_KELVIN) ? "K" : "C",
						FBUtility::convK2Temp(iv.Temperature, Unit.Temp), 
						FBUtility::convK2ND(iv.Temperature, BaseTemp, DiffTemp));
  }
  
  fprintf(fp,"\n");
  fflush(fp);
}


// #################################################################
/**
 * @brief 線形ソルバー種別の表示
 * @param [in] fp ファイルポインタ
 * @param [in] IC IterationCtl
 */
void Control::printLS(FILE* fp, const IterationCtl* IC)
{
  switch (IC->getLS()) 
  {
    case JACOBI:
      fprintf(fp,"\t       Linear Solver          :   Jacobi method\n");
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
      
    case GMRES:
      fprintf(fp,"\t       Linear Solver          :   GMRES\n");
      break;
      
    case RBGS:
      fprintf(fp,"\t       Linear Solver          :   RBGS\n");
      break;
      
    case PCG:
      fprintf(fp,"\t       Linear Solver          :   PCG\n");
      break;
      
    case PBiCGSTAB:
      fprintf(fp,"\t       Linear Solver          :   PBiCGSTAB\n");
      break;
      
    default:
      stamped_printf("Error: Linear Solver section\n");
      Exit(0);
  }
}


// #################################################################
// 内部BCコンポーネントの数を表示する
void Control::printNoCompo(FILE* fp)
{
  fprintf(fp,"\tNo. of Local Boundary  : %d\n", NoBC);
  fprintf(fp,"\tNo. of Medium          : %d\n", NoMedium);
  fprintf(fp,"\n");
  fprintf(fp,"\tNo. of Fluid Medium    : %d\n", NoMediumFluid);
  fprintf(fp,"\tNo. of Solid Medium    : %d\n", NoMediumSolid);
}


// #################################################################
/**
 * @brief 計算パラメータの表示
 * @param [in] fp  ファイルポインタ
 * @param [in] mat MediumListクラス
 */
void Control::printParaConditions(FILE* fp, const MediumList* mat)
{
  if ( !fp )
  {
    stamped_printf("\tFail to write into file\n");
    Exit(0);
  }
  
  fprintf(fp,"\n---------------------------------------------------------------------------\n\n");
  fprintf(fp,"\n\t>> Simulation Parameters\n\n");
  
  fprintf(fp,"\tReference ID              [-]         :  %s\n", mat[RefMat].getAlias().c_str());
  fprintf(fp,"\n");
  
  // Reference values
  fprintf(fp,"\tRef. Length               [m]         : %12.5e\n", RefLength);
  fprintf(fp,"\tRef. Velocity             [m/s]       : %12.5e\n", RefVelocity);
	fprintf(fp,"\tBase Pressure             [Pa]        : %12.5e\n", BasePrs);
  fprintf(fp,"\tRef. Mass Density         [kg/m^3]    : %12.5e\n", RefDensity);
  fprintf(fp,"\tRef. Viscosity            [Pa s]      : %12.5e\n", RefViscosity);
  fprintf(fp,"\tRef. Knmtc Viscosity      [m^2/s]     : %12.5e\n", RefKviscosity);
  fprintf(fp,"\tRef. Specific Heat        [J/(kg K)]  : %12.5e\n", RefSpecificHeat);
  fprintf(fp,"\tRef. Thermal Conductivity [W/(m K)]   : %12.5e\n", RefLambda);
  fprintf(fp,"\tRef. Speed of Sound       [m/s]       : %12.5e\n", RefSoundSpeed);
  fprintf(fp,"\tGravity                   [m/s^2]     : %12.5e\n", Gravity);
  fprintf(fp,"\n");
  
  fprintf(fp,"\tSpacing                   [m] / [-]   : %12.5e / %12.5e\n", deltaX*RefLength, deltaX);
  fprintf(fp,"\tTime Scale                [sec]       : %12.5e\n", Tscale);
  fprintf(fp,"\n");
  
  if ( isHeatProblem() )
  {
    fprintf(fp,"\tBase Temperature          [%s] / [-]   : %12.5e / %3.1f\n", (Unit.Temp==Unit_KELVIN) ? "K" : "C", FBUtility::convK2Temp(BaseTemp, Unit.Temp), 0.0);
    fprintf(fp,"\tTemperature Diff.         [%s] / [-]   : %12.5e / %3.1f\n", (Unit.Temp==Unit_KELVIN) ? "K" : "C", DiffTemp, 1.0);
  }
  fprintf(fp,"\n");
  
  //REAL_TYPE ay, ap;
  //ay = yaw/180.0*2.0*asin(1.0);
  //ap = pitch/180.0*2.0*asin(1.0);
  //fprintf(fp,"\tYaw   Angle               [rad]/[deg] : %12.5e / %12.5e\n", ay, yaw);
  //fprintf(fp,"\tPitch Angle               [rad]/[deg] : %12.5e / %12.5e\n", ap, pitch);
  
  fprintf(fp,"\n");
  fprintf(fp,"\tPrandtl  number           [-]         : %12.5e\n", Prandtl);
  if (Mode.PDE == PDE_NS)
  {
    fprintf(fp,"\tReynolds number           [-]         : %12.5e\n", Reynolds);
  }
  else
  {
    fprintf(fp,"\tMach number               [-]         : %12.5e\n", Mach);
  }
  
  if ( isHeatProblem() )
  {
    fprintf(fp,"\tPeclet   number           [-]         : %12.5e\n", Peclet);
    fprintf(fp,"\tGrashof  number           [-]         : %12.5e\n", Grashof);
    if (KindOfSolver==THERMAL_FLOW_NATURAL)  fprintf(fp,"\tRayleigh number           [-]         : %12.5e\n", Rayleigh);
  }
  fprintf(fp,"\n");
  
  fflush(fp);
}


// #################################################################
// 制御パラメータSTEERの表示
// 20130611
//void Control::printSteerConditions(FILE* fp, const IterationCtl* IC, const DTcntl* DT, const ReferenceFrame* RF, FileIO_PLOT3D_WRITE* FP3DW)
void Control::printSteerConditions(FILE* fp, IterationCtl* IC, const DTcntl* DT, const ReferenceFrame* RF)
{
  if ( !fp )
  {
    stamped_printf("\tFail to write into file\n");
    Exit(0);
  }
  
  double dt = DT->get_DT();
  bool  err=true;
  double itm=0.0;
  
  fprintf(fp,"\n---------------------------------------------------------------------------\n\n");
  fprintf(fp,"\n\t>> Library Information\n\n");
  fprintf(fp,"\t     CPMlib     Version %s\n", ver_CPM.c_str());
  fprintf(fp,"\t     CIOlib     Version %s\n", ver_CIO.c_str());
  fprintf(fp,"\t     Polylib    Version %s\n", ver_Poly.c_str());
  fprintf(fp,"\t     Cutlib     Version %s\n", ver_CUT.c_str());
  fprintf(fp,"\t     PMlib      Version %s\n", ver_PM.c_str());
  fprintf(fp,"\t     TextParser Version %s\n", ver_TP.c_str());
  fprintf(fp,"\n");
  
  fprintf(fp,"\n---------------------------------------------------------------------------\n\n");
  fprintf(fp,"\n\t>> Solver Control Parameters\n\n");
  
  // ソルバープロパティ ------------------
  fprintf(fp,"\tSolver Properties\n");
  
  // Basic Equation and PDE
	switch (KindOfSolver)
  {
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
					if  ( !((KindOfSolver == FLOW_ONLY) || (KindOfSolver == THERMAL_FLOW) || (KindOfSolver == THERMAL_FLOW_NATURAL)) )
          {
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
  switch (Mode.Steady)
  {
    case FB_STEADY:
      fprintf(fp,"\t     Time Variation           :   Steady\n");
      break;
      
    case FB_UNSTEADY:
      fprintf(fp,"\t     Time Variation           :   Unsteady\n");
      break;
      
    default:
      stamped_printf("Error: Time Variation[%d]\n", Mode.Steady);
      err=false;
  }
  
  // Shape approximation
  switch (Mode.ShapeAprx)
  {
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
  if ( Mode.Precision == sizeof(float) )
  {
    fprintf(fp,"\t     Precision                :   Single Precision \n");
  }
  else
  {
    fprintf(fp,"\t     Precision                :   Double Precision \n");
  }
  
  // Kind Of Solver
  if (KindOfSolver==FLOW_ONLY)
  {
    fprintf(fp,"\t     Kind of Solver           :   Flow Only (Non Heat)\n");
  }
  else if ( (KindOfSolver==THERMAL_FLOW) && (Mode.Buoyancy==NO_BUOYANCY) )
  {
    fprintf(fp,"\t     Kind of Solver           :   Forced convection without buoyancy\n");
  }
  else if ( (KindOfSolver==THERMAL_FLOW) && (Mode.Buoyancy==BOUSSINESQ) )
  {
    fprintf(fp,"\t     Kind of Solver           :   Forced convection with buoyancy : Boussinesq Approximation\n");
  }
  else if ( (KindOfSolver==THERMAL_FLOW) && (Mode.Buoyancy==LOW_MACH) )
  {
    fprintf(fp,"\t     Kind of Solver           :   Forced convection with buoyancy : Low Mach Approximation\n");
  }
  else if ( (KindOfSolver==THERMAL_FLOW_NATURAL) && (Mode.Buoyancy==BOUSSINESQ) )
  {
    fprintf(fp,"\t     Kind of Solver           :   Natural convection : Boussinesq Approximation\n");
  }
  else if ( (KindOfSolver==THERMAL_FLOW_NATURAL) && (Mode.Buoyancy==LOW_MACH) )
  {
    fprintf(fp,"\t     Kind of Solver           :   Natural convection : Low Mach Approximation\n");
  }
  else if (KindOfSolver==SOLID_CONDUCTION)
  {
    fprintf(fp,"\t     Kind of Solver           :   Solid Conduction\n");
  }
  else
  {
    fprintf(fp,"\t     Heat Solver type         :   Error\n");
    err=false;
  }
  
  // Flow Algorithm
  switch (AlgorithmF)
  {
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
  if ( isHeatProblem() )
  {
    switch (AlgorithmH)
    {
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
  }


  // Convection scheme
	if ( KindOfSolver != SOLID_CONDUCTION )
  {
		switch (CnvScheme)
    {
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
  switch (RF->getFrame()) 
  {
    case ReferenceFrame::frm_static:
      fprintf(fp,"\t     Reference Frame          :   Stationary\n");
      break;
      
    case ReferenceFrame::frm_translation:
      if (Unit.Param==DIMENSIONAL)
      {
        fprintf(fp,"\t     Reference Frame          :   Translational (%12.4e, %12.4e, %12.4e) [m/s]\n", GridVel[0]*RefVelocity, GridVel[1]*RefVelocity, GridVel[2]*RefVelocity);
      }
      else
      {
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


  // 単位系 ------------------
  fprintf(fp,"\n\tUnit\n");
  fprintf(fp,"\t     Unit of Input Parameter  :   %s\n", (Unit.Param == DIMENSIONAL) ? "Dimensional" : "Non-Dimensional");
  fprintf(fp,"\t             Pressure         :   %s\n", (Unit.Prs == Unit_Absolute) ? "Absolute Pressure" : "Gauge Pressure");
  fprintf(fp,"\t             Temperature      :   %s\n", (Unit.Temp == Unit_KELVIN) ? "Kelvin" : "Celsius");
  
  
  // 時間制御 ------------------
  fprintf(fp,"\n\tTime Control\n");
  
  // 加速時間
  if ( !Interval[Interval_Manager::tg_accelra].isStep() )
  {
    itm = Interval[Interval_Manager::tg_accelra].getIntervalTime();
    fprintf(fp,"\t     Acceleration Time        :   %12.5e [sec] / %12.5e [-]\n", itm*Tscale, itm);
  }
  else
  {
    fprintf(fp,"\t     Acceleration Step        :   %12d\n", Interval[Interval_Manager::tg_accelra].getIntervalStep());
  }
  
  // 時間平均
  if ( Mode.Average == ON )
  {
    if ( !Interval[Interval_Manager::tg_average].isStep() )
    {
      itm = Interval[Interval_Manager::tg_average].getStartTime();
      fprintf(fp,"\t     Averaging Start          :   %12.5e [sec] / %12.5e [-]\n", itm*Tscale, itm);
    }
    else
    {
      fprintf(fp,"\t     Averaging Start          :   %12d\n", Interval[Interval_Manager::tg_average].getStartStep());
    }
  }
  else
  {
    fprintf(fp,"\t     Averaging Start          :   OFF\n");
  }
  
  // Time Increment
  REAL_TYPE d_R = deltaX*deltaX*Reynolds/6.0; // 拡散数
  REAL_TYPE d_P = deltaX*deltaX*Peclet/6.0;   // 拡散数
  REAL_TYPE cfl = (REAL_TYPE)DT->get_CFL();
  switch ( DT->get_Scheme() ) 
  {
    case DTcntl::dt_direct:
      fprintf(fp,"\t     Time Increment dt        :   %12.5e [sec] / %12.5e [-] : Direct ", dt*Tscale, dt);
      if ( isHeatProblem() )
      {
        fprintf(fp,": Diff. Num. = %7.2e\n", dt/(deltaX*deltaX*Peclet));
      }
      else
      {
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
      if ( isHeatProblem() )
      {
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
      if ( isHeatProblem() )
      {
        c = (REAL_TYPE)DT->dtDFN((double)Peclet);
        fprintf(fp,"\t                              :             dt restricted by Diffusion number (Peclet)     : %8.5f [-]\n", d_P);
      }
    }
      break;
      
    case DTcntl::dt_cfl_max_v_cp:
      fprintf(fp,"\t     Time Increment dt        :   %12.5e [sec] / %12.5e [-] : CFL (/w Speed of Sound) & Diffusion number with Maximum velocity\n", dt*Tscale, dt);
      break;
      
    default:
      stamped_printf("Error: Time Increment section\n");
      err=false;
  }
  
  // Calculation time/step
  itm = Interval[Interval_Manager::tg_compute].getIntervalTime();
  fprintf(fp,"\t     Calculation Time         :   %12.5e [sec] / %12.5e [-]\n", itm*Tscale, itm);
  fprintf(fp,"\t     Calculation Step         :   %12d\n", Interval[Interval_Manager::tg_compute].getIntervalStep());
  
  
  
  // スタートモード ------------------
  fprintf(fp,"\n\tStart Mode\n");
  
  // Start
  switch (Start)
  {
    case initial_start:
      fprintf(fp,"\t     Start Condition          :   Impulsive start\n");
      break;
      
    case restart_sameDiv_sameRes:
      fprintf(fp,"\t     Start Condition          :   Restart with same resolution and same num. of division\n");
      break;
      
    case restart_sameDiv_refinement:
      fprintf(fp,"\t     Start Condition          :   Restart with refinment and same num. of division\n");
      break;
      
    case restart_diffDiv_sameRes:
      fprintf(fp,"\t     Start Condition          :   Restart with same resolution and different division\n");
      break;
      
    case restart_diffDiv_refinement:
      fprintf(fp,"\t     Start Condition          :   Restart with refinment and different division\n");
      break;
      
    default:
      stamped_printf("Error: start condition section\n");
      err=false;
  }

  
  // parallel mode ------------------
  fprintf(fp,"\n\tParallel Mode\n");
  
  switch (Parallelism)
  {
    case Serial:
      fprintf(fp,"\t     Parallel Mode            :   Serial\n");
      break;
      
    case OpenMP:
      fprintf(fp,"\t     Parallel Mode            :   OpenMP  (%d threads)\n", num_thread);
      break;
      
    case FlatMPI:
      fprintf(fp,"\t     Parallel Mode            :   Flat MPI  (%d processes)\n", num_process);
      break;
      
    case Hybrid:
      fprintf(fp,"\t     Parallel Mode            :   Hybrid  (%d processes x %d threads)\n", num_process, num_thread);
      break;
      
    default:
      break;
  }
  
  // 空間分割
  fprintf(fp,"\t     Space Partitioning       :   Equal Partitioning\n");
  
  
  
  
  // File IO mode ------------------
  fprintf(fp,"\n\tFile IO Mode\n");
  
  fprintf(fp,"\t     Unit of File             :   %s\n", (Unit.File == DIMENSIONAL) ? "Dimensional" : "Non-Dimensional");
  
  // InputMode >> Distributed by default
  fprintf(fp,"\t     IO Mode                  :   %s\n", (FIO.IOmode==IO_GATHER) ? "Gathered" : "Distributed");
  
  
  // Output guide
  fprintf(fp,"\t     Guide cell for output    :   %d\n", GuideOut);
  
  // Voxel output
  fprintf(fp,"\t     Voxel model output       :   %s\n", (FIO.IO_Voxel==Sphere_SVX) ? "svx" : "None");
  
  
  // IO Directory
  fprintf(fp,"\t     I/O Directory Input      :   \"%s\"\n", FIO.InDirPath.c_str());
  fprintf(fp,"\t     I/O Directory Output     :   \"%s\"\n", FIO.OutDirPath.c_str());
  
  // Time Slice option
  fprintf(fp,"\t     Time Slie Directory      :   %s\n", (FIO.Slice==ON) ? "On" : "Off");
  
  
  // ログ出力 ------------------
  fprintf(fp,"\n\tLogs\n");
  fprintf(fp,"\t     Unit for Output          :   %s\n", (Unit.Log == DIMENSIONAL) ? "Dimensional" : "Non-Dimensional");
  fprintf(fp,"\t     Base Logs                :   %4s  %s, %s, %s, %s\n", 
          (Mode.Log_Base == ON) ? "ON >" : "OFF ", 
          (Mode.Log_Base == ON) ? HistoryName.c_str() : "", 
          (Mode.Log_Base == ON) ? HistoryCompoName.c_str() : "", 
          (Mode.Log_Base == ON) ? HistoryDomfxName.c_str() : "", 
          (Mode.Log_Base == ON) ? HistoryForceName.c_str() : "");
  fprintf(fp,"\t     Iteration Log            :   %4s  %s\n", 
          (Mode.Log_Itr == ON)?"ON >":"OFF ", (Mode.Log_Itr == ON) ? HistoryItrName.c_str() : "");
  fprintf(fp,"\t     Profiling report         :   %4s  %s%s\n", 
          (Mode.Profiling != OFF)?"ON >":"OFF ", 
          (Mode.Profiling == DETAIL)? "Detail mode, ":"",
          (Mode.Profiling != OFF)?"profiling.txt":"");
  
  fprintf(fp,"\t     Wall info. Log           :   %4s  %s\n", 
          (Mode.Log_Wall == ON)?"ON >":"OFF ", (Mode.Log_Wall == ON) ? HistoryWallName.c_str() : "");
  
  
  // Intervals
  fprintf(fp,"\n\tIntervals\n");
  
  
  // インターバル ------------------
  // 基本履歴のコンソール出力 
  if ( !Interval[Interval_Manager::tg_console].isStep() )
  {
    itm = Interval[Interval_Manager::tg_console].getIntervalTime();
    fprintf(fp,"\t     Base Info.(stdout)       :   %12.6e [sec] / %12.6e [-]\n", itm*Tscale, itm);
  }
  else
  {
    fprintf(fp,"\t     Base Info.(stdout)       :   %12d [step]\n", Interval[Interval_Manager::tg_console].getIntervalStep());
  }
  
  // 履歴情報のファイル出力
  if ( !Interval[Interval_Manager::tg_history].isStep() )
  {
    itm = Interval[Interval_Manager::tg_history].getIntervalTime();
    fprintf(fp,"\t     Other Histories          :   %12.6e [sec] / %12.6e [-]\n", itm*Tscale, itm);
  }
  else
  {
    fprintf(fp,"\t     Other Histories          :   %12d [step]\n", Interval[Interval_Manager::tg_history].getIntervalStep());
  }
  
  // 基本変数のファイル出力
  if ( !Interval[Interval_Manager::tg_basic].isStep() )
  {
    itm = Interval[Interval_Manager::tg_basic].getIntervalTime();
    fprintf(fp,"\t     Basic Variables          :   %12.6e [sec] / %12.6e [-]\n", itm*Tscale, itm);
  }
  else
  {
    fprintf(fp,"\t     Basic Variables          :   %12d [step]\n", Interval[Interval_Manager::tg_basic].getIntervalStep());
  }
  
  // 平均値のファイル出力
  if ( !Interval[Interval_Manager::tg_average].isStep() )
  {
    itm = Interval[Interval_Manager::tg_average].getIntervalTime();
    fprintf(fp,"\t     Averaged Variables       :   %12.6e [sec] / %12.6e [-]\n", itm*Tscale, itm);
  }
  else
  {
    fprintf(fp,"\t     Averaged Variables       :   %12d [step]\n", Interval[Interval_Manager::tg_average].getIntervalStep());
  }
  
  // 派生変数のファイル出力
  if ( !Interval[Interval_Manager::tg_derived].isStep() )
  {
    itm = Interval[Interval_Manager::tg_derived].getIntervalTime();
    fprintf(fp,"\t     Derived Variables        :   %12.6e [sec] / %12.6e [-]\n", itm*Tscale, itm);
  }
  else
  {
    fprintf(fp,"\t     Derived Variables        :   %12d [step]\n", Interval[Interval_Manager::tg_derived].getIntervalStep());
  }
  
  // サンプリング情報のファイル出力
  if ( Sampling.log == ON )
  {
    if ( !Interval[Interval_Manager::tg_sampled].isStep() )
    {
      itm = Interval[Interval_Manager::tg_sampled].getIntervalTime();
      fprintf(fp,"\t     Sampled data             :   %12.6e [sec] / %12.6e [-]\n", itm*Tscale, itm);
    }
    else
    {
      fprintf(fp,"\t     Sampled data             :   %12d [step]\n", Interval[Interval_Manager::tg_sampled].getIntervalStep());
    }
  }
  
  // PLOT3D 瞬間値のファイル出力
  if(FIO.Format == plt3d_fmt)
  {
    if ( !Interval[Interval_Manager::tg_plot3d].isStep() )
    {
      itm = Interval[Interval_Manager::tg_plot3d].getIntervalTime();
      fprintf(fp,"\t     Plot3d Instant data      :   %12.6e [sec] / %12.6e [-]\n", itm*Tscale, itm);
    }
    else
    {
      fprintf(fp,"\t     Plot3d Instant data      :   %12d [step]\n", Interval[Interval_Manager::tg_plot3d].getIntervalStep());
    }
  }


  /* 20130611 PLOT3D Options
  if(FIO.Format == plt3d_fmt)
  {
    fprintf(fp,"\n\tPLOT3D Options\n");
    fprintf(fp,"\t     grid file prefix         :   %s\n", P3Op.basename_g.c_str());
    fprintf(fp,"\t     function file prefix     :   %s\n", P3Op.basename_f.c_str());
    fprintf(fp,"\t     grid kind                :   %s\n", (FP3DW->IsGridKind()) ? "multi grid" : "single grid");
    fprintf(fp,"\t     grid mobility            :   %s\n", (FP3DW->IsMoveGrid()) ? "movable" : "immovable");
    fprintf(fp,"\t     state of time            :   %s\n", (FP3DW->IsSteady()) ? "unsteady" : "steady");
    fprintf(fp,"\t     output iblank            :   %s\n", (FP3DW->IsIBlankFlag()) ? "on" : "off");
    if (      FP3DW->GetFormat() == UNFORMATTED ) fprintf(fp,"\t     output format            :   %s\n", "Fortran Unformatted");
    else if ( FP3DW->GetFormat() == FORMATTED   ) fprintf(fp,"\t     output format            :   %s\n", "Fortran Formatted");
    else if ( FP3DW->GetFormat() == C_BINARY    ) fprintf(fp,"\t     output format            :   %s\n", "C Binary");
    fprintf(fp,"\t     output dimention         :   %iD\n", FP3DW->GetDim());
    if (      FP3DW->GetRealType() == OUTPUT_FLOAT  ) fprintf(fp,"\t     output format            :   %s\n", "float");
    else if ( FP3DW->GetRealType() == OUTPUT_DOUBLE ) fprintf(fp,"\t     output format            :   %s\n", "double");
    fprintf(fp,"\t     output xyz file          :   %s\n", (P3Op.IS_xyz) ? "on" : "off");
    fprintf(fp,"\t     output q file            :   %s\n", (P3Op.IS_q) ? "on" : "off");
    fprintf(fp,"\t     output function file     :   %s\n", (P3Op.IS_funciton) ? "on" : "off");
    fprintf(fp,"\t     output funciton name file:   %s\n", (P3Op.IS_function_name) ? "on" : "off");
    fprintf(fp,"\t     output fvbnd file        :   %s\n", (P3Op.IS_fvbnd) ? "on" : "off");
    fprintf(fp,"\t     function per item        :   %s\n", (P3Op.IS_DivideFunc) ? "on" : "off");
  }
  */

  
  // Criteria ------------------
  fprintf(fp,"\n\tParameter of Linear Equation\n");
  IterationCtl* ICp1= &IC[ic_prs1];  /// 圧力のPoisson反復
  IterationCtl* ICp2= &IC[ic_prs2];  /// 圧力のPoisson反復　2回目
  IterationCtl* ICv = &IC[ic_vel1];  /// 粘性項のCrank-Nicolson反復
  IterationCtl* ICd = &IC[ic_div];   /// V-P反復
  
  if ( Hide.PM_Test == ON )
  {
    fprintf(fp,"\t ### Performance Test Mode >> The iteration number is fixed by Iteration max.\n\n");
  }

	if ( KindOfSolver != SOLID_CONDUCTION )
  {
    // V-P iteration
		fprintf(fp,"\t     V-P Iteration \n");
		fprintf(fp,"\t       Iteration max          :   %d\n"  ,  ICd->getMaxIteration());
		fprintf(fp,"\t       Convergence eps        :   %9.3e\n", ICd->getCriterion());
		fprintf(fp,"\t       Norm type              :   %s\n",    ICd->getNormString().c_str());
    
    
		// 1st iteration
		fprintf(fp,"\t     1st Pressure Iteration \n");
		fprintf(fp,"\t       Iteration max          :   %d\n"  ,  ICp1->getMaxIteration());
		fprintf(fp,"\t       Convergence eps        :   %9.3e\n", ICp1->getCriterion());
		fprintf(fp,"\t       Coef. of Relax./Accel. :   %9.3e\n", ICp1->getOmega());
		fprintf(fp,"\t       Norm type              :   %s\n",    ICp1->getNormString().c_str());
    fprintf(fp,"\t       Communication Mode     :   %s\n",   (ICp1->getSyncMode()==comm_sync) ? "SYNC" : "ASYNC");
		printLS(fp, ICp1);
    
    if ( AlgorithmF == Flow_FS_RK_CN )
    {
      fprintf(fp,"\t     2nd Pressure Iteration \n");
      fprintf(fp,"\t       Iteration max          :   %d\n"  ,  ICp2->getMaxIteration());
      fprintf(fp,"\t       Convergence eps        :   %9.3e\n", ICp2->getCriterion());
      fprintf(fp,"\t       Coef. of Relax./Accel. :   %9.3e\n", ICp2->getOmega());
      fprintf(fp,"\t       Norm type              :   %s\n",    ICp2->getNormString().c_str());
      fprintf(fp,"\t       Communication Mode     :   %s\n",   (ICp2->getSyncMode()==comm_sync) ? "SYNC" : "ASYNC");
      printLS(fp, ICp2);
    }
    
    // CN iteration
		if ( (AlgorithmF == Flow_FS_AB_CN) || (AlgorithmF == Flow_FS_RK_CN) )
    {
      fprintf(fp,"\n");
			fprintf(fp,"\t     Velocity CN Iteration \n");
			fprintf(fp,"\t       Iteration max           :   %d\n"  ,  ICv->getMaxIteration());
			fprintf(fp,"\t       Convergence eps         :   %9.3e\n", ICv->getCriterion());
			fprintf(fp,"\t       Coef. of Relax./Accel.  :   %9.3e\n", ICv->getOmega());
			fprintf(fp,"\t       Norm type               :   %s\n",    ICv->getNormString().c_str());
      fprintf(fp,"\t       Communication Mode      :   %s\n",   (ICv->getSyncMode()==comm_sync) ? "SYNC" : "ASYNC");
			printLS(fp, ICv);
		}
	}

  // for Temperature
  if ( isHeatProblem() )
  {
    if ( AlgorithmH == Heat_EE_EI )
    {
      IterationCtl* ICt = &IC[ic_tmp1];  /// 温度の拡散項の反復
      fprintf(fp,"\n");
      fprintf(fp,"\t     Temperature Iteration  \n");
      fprintf(fp,"\t       Iteration max          :   %d\n"  ,  ICt->getMaxIteration());
      fprintf(fp,"\t       Convergence eps        :   %9.3e\n", ICt->getCriterion());
      fprintf(fp,"\t       Coef. of Relax./Accel. :   %9.3e\n", ICt->getOmega());
			fprintf(fp,"\t       Norm type              :   %s\n",    ICt->getNormString().c_str());
      fprintf(fp,"\t       Communication Mode     :   %s\n",   (ICt->getSyncMode()==comm_sync) ? "SYNC" : "ASYNC");
      printLS(fp, ICt);
    }
  }


  // 壁面の扱い ------------------
  fprintf(fp,"\n\tCondition of Wall\n");
  fprintf(fp,"\t     No of surface            :   %ld\n", NoWallSurface);
  switch (Mode.PrsNeuamnnType)
  {
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
  
  switch (Mode.Wall_profile)
  {
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
  
  
  // 派生変数 ------------------
  fprintf(fp, "\n\tDerived variables\n");
  
  //　全圧の出力モード
  if ( Mode.TP == ON )
  {
    fprintf(fp,"\t     Total Pressure           :   ON\n");
  }
  else
  {
    fprintf(fp,"\t     Total Pressure           :   OFF\n");
  }
  
  //　渦度の出力モード
  if ( Mode.VRT == ON )
  {
    fprintf(fp,"\t     Vorticity                :   ON\n");
  }
  else
  {
    fprintf(fp,"\t     Vorticity                :   OFF\n");
  }
  
  //　ヘリシティの出力モード
  if ( Mode.Helicity == ON )
  {
    fprintf(fp,"\t     Helicity                 :   ON\n");
  }
  else
  {
    fprintf(fp,"\t     Helicity                 :   OFF\n");
  }
  
  //　速度勾配テンソルの第二不変量の出力モード
  if ( Mode.I2VGT == ON ) {
    fprintf(fp,"\t     2nd Invariant of VGT     :   ON\n");
  }
  else
  {
    fprintf(fp,"\t     2nd Invariant of VGT     :   OFF\n");
  }
  
  // FaceVelocity
  if ( Mode.FaceV == ON ) {
    fprintf(fp,"\t     Face Velocity(Staggered) :   ON\n");
  }
  else
  {
    fprintf(fp,"\t     Face Velocity(Staggered) :   OFF\n");
  }
  
  
  // Hidden parameter -----------------
  
  if (Hide.Change_ID != 0 )
  {
    fprintf(fp,"\t     Change ID from [0] to    :   %d\n", Hide.Change_ID);
  }
  fprintf(fp,"\n");
  
  if (Hide.Range_Limit == Range_Cutoff)
  {
    fprintf(fp,"\t     Variable Range           :   Limit value between [0,1] in normalized value\n");
  }
  
  fflush(fp);
  
  if (err==false) Exit(0);
}



// #################################################################
// コンポーネントが存在するかを保持しておく
void Control::setExistComponent(CompoList* cmp, BoundaryOuter* OBC)
{
  int c;
  
  // Vspec
  c = 0;
  for (int n=1; n<=NoCompo; n++)
  {
    if ( cmp[n].isVBC_IO() ) c++;
  }
  if ( c>0 ) EnsCompo.vspec = ON;
  
  
  // Forcing > HEX, FAN, DARCY
  c = 0;
  for (int n=1; n<=NoCompo; n++)
  {
    if ( cmp[n].isFORCING() ) c++;
  }
  if ( c>0 ) EnsCompo.forcing = ON;
  
  
  // Heat source > HEAT_SRC, CNST_TEMP
  c = 0;
  for (int n=1; n<=NoCompo; n++)
  {
    if ( cmp[n].isHsrc() ) c++;
  }
  if ( c>0 ) EnsCompo.hsrc = ON;
  
  
  // 周期境界 > PERIODIC
  c = 0;
  for (int n=1; n<=NoCompo; n++)
  {
    if ( cmp[n].getType() == PERIODIC ) c++;
  }
  if ( c>0 ) EnsCompo.periodic = ON;
  
  
  // 流出境界 > OUTFLOW
  c = 0;
  for (int n=1; n<=NoCompo; n++)
  {
    if ( cmp[n].getType() == OUTFLOW ) c++;
  }
  if ( c>0 ) EnsCompo.outflow = ON;
  
  
  // 体積率コンポーネント
  c = 0;
  for (int n=1; n<=NoCompo; n++)
  {
    if ( cmp[n].isVFraction() ) c++;
  }
  if ( c>0 ) EnsCompo.fraction = ON;
  
  
  // トラクションフリー
  c = 0;
  for (int n=0; n<NOFACE; n++)
  {
    if ( OBC[n].getClass()== OBC_TRC_FREE ) c++;
  }
  if ( c>0 ) EnsCompo.tfree = ON;
  
  
  // モニタ
  c = 0;
  for (int n=1; n<=NoCompo; n++)
  {
    if ( cmp[n].isMONITOR() ) c++;
  }
  if ( c>0 ) EnsCompo.monitor = ON;
  
  
  // チェック　MONITOR_LISTでCELL_MONITORが指定されている場合
  if ( (existMonitor() == ON) && (c < 1) )
  {
    Hostonly_ stamped_printf("\tError : CellMonitor in MonitorList is specified, however any Monitor can not be found.\n");
    Exit(0);
  }
  
  if ( (existMonitor() == OFF) && (c > 0) )
  {
    Hostonly_ stamped_printf("\tError : CellMonitor in MonitorList is NOT specified, however Monitor section is found in LocalBoundary.\n");
    Exit(0);
  }
  
  // @see C.EnsCompo.monitor==ONは，Control::getMonitorList()
  
}


// #################################################################
/**
 @brief 無次元パラメータを各種モードに応じて設定する
 @param mat
 @param cmp
 @param rf
 @note
 - 代表長さと代表速度はパラメータで必ず与えること（読み込んだ値は変更しない）
 - 純強制対流　有次元　（代表長さ，代表速度，動粘性係数，温度拡散係数）
 -           無次元　（Pr, Re > RefV=RefL=1）
 - 熱対流　　　有次元　（代表長さ，代表速度，温度差，体膨張係数，重力加速度，動粘性係数，温度拡散係数）
 - 自然対流　　有次元　（代表長さ，代表速度，温度差，体膨張係数，重力加速度，動粘性係数，温度拡散係数）
 - 固体熱伝導　有次元　（代表長さ，温度拡散係数 > Peclet=1）？
 * @see
 *  - bool Control::getXML_Para_ND(void)
 *  - void Control::getXML_Para_Init(void)
 */
void Control::setParameters(MediumList* mat, CompoList* cmp, ReferenceFrame* RF, BoundaryOuter* BO)
{
  REAL_TYPE rho, nyu, cp, lambda, beta, mu, snd_spd=0.0;
  REAL_TYPE c1, c2, c3;
  int m;
  
  // get reference values
  for (int n=1; n<=NoMedium; n++)
  {
    if ( n == RefMat )
    {
      if ( mat[n].getState() == FLUID )
      {
        rho     = mat[n].P[p_density];
        mu      = mat[n].P[p_viscosity];
        nyu     = mat[n].P[p_kinematic_viscosity];
        cp      = mat[n].P[p_specific_heat];
        lambda  = mat[n].P[p_thermal_conductivity];
        beta    = mat[n].P[p_vol_expansion]; // can be replaced by 1/K in the case of gas
        snd_spd = mat[n].P[p_speed_of_sound];
      }
      else
      {
        rho    = mat[n].P[p_density];
        cp     = mat[n].P[p_specific_heat];
        lambda = mat[n].P[p_thermal_conductivity];
      }
    }
  }
  
  // RefDensity      = rho;
  RefViscosity    = mu;
  RefKviscosity   = nyu;
  RefSpecificHeat = cp;
  RefLambda       = lambda;
  RefSoundSpeed   = snd_spd;
  
  if (KindOfSolver == SOLID_CONDUCTION)
  {
    if (Unit.Param == DIMENSIONAL)
    {
      Peclet   = 1.0;
    }
    else
    {
      Hostonly_ printf("Error : Solid conduction(ND)\n");
			Exit(0);
    }
  }
  else if (KindOfSolver == FLOW_ONLY)
  {
    if (Unit.Param == DIMENSIONAL)
    {
      Reynolds = RefVelocity * RefLength / nyu;
      Prandtl  = rho * cp * nyu / lambda;
    }
  }
	else if (KindOfSolver==THERMAL_FLOW)
  {
    switch (Mode.Buoyancy)
    {
      case NO_BUOYANCY:
        if (Unit.Param == DIMENSIONAL)
        {
          Reynolds = RefVelocity * RefLength / nyu;
          Prandtl  = rho * cp * nyu / lambda;
          Peclet   = Prandtl * Reynolds;
        }
        else
        {
          Hostonly_ printf("Error : Thermal flow /wo buoyancy(ND)\n");
					Exit(0);
        }
        break;
        
      case BOUSSINESQ:
      case LOW_MACH:
        if (Unit.Param == DIMENSIONAL)
        {
          Reynolds = RefVelocity * RefLength / nyu;
          Prandtl  = rho * cp * nyu / lambda;
          Peclet   = Prandtl * Reynolds;
          c1 = beta / nyu;
          c2 = Gravity *  DiffTemp / nyu;
          c3 = RefLength*RefLength*RefLength;
          Grashof  = c1 * c2 * c3;
          Rayleigh = Prandtl * Grashof;
        }
        else
        {
          Hostonly_ printf("Error : Thermal flow /w buoyancy(ND)\n");
					Exit(0);
        }
        break;
    }
  }
  else if (KindOfSolver==THERMAL_FLOW_NATURAL)
  {
    switch (Mode.Buoyancy)
    {
      case BOUSSINESQ:
      case LOW_MACH:
        if (Unit.Param == DIMENSIONAL)
        {
          Prandtl  = rho * cp * nyu / lambda;
          Reynolds = RefVelocity * RefLength / nyu;
					Peclet   = Prandtl * Reynolds;
          c1 = beta / nyu;
          c2 = Gravity *  DiffTemp / nyu;
          c3 = RefLength*RefLength*RefLength;
          Grashof  = c1 * c2 * c3;
          Rayleigh = Prandtl * Grashof;
        }
        else
        {
          Hostonly_ printf("Error : Natural Convection(ND)\n");
					Exit(0);
        }
        break;
        
      default:
        break;
		}
	}
	else
  { // CONJUGATE_HEAT_TRANSFER
		;
	}
  
  if (Mode.PDE == PDE_EULER) Reynolds=1.0e23;
  
  Mach = RefVelocity / RefSoundSpeed;
  
  // タイミングパラメータの無次元化
  Tscale = (double)RefLength / (double)RefVelocity;
  
  
  if ( Unit.Param == DIMENSIONAL )
  {
    double g[3];
    RF->copyGridVel(g);
    g[0] /= (double)RefVelocity;
    g[1] /= (double)RefVelocity;
    g[2] /= (double)RefVelocity;
    RF->setGridVel(g);
  }
  
  // コンポーネントの指定速度
  for (int n=1; n<=NoCompo; n++)
  {
    if ( (cmp[n].getType()==SPEC_VEL_WH) || (cmp[n].getType()==SPEC_VEL) )
    {
			if ( cmp[n].isPolicy_Massflow() ) //ポリシーが流量の場合
      {
				if ( Unit.Param == DIMENSIONAL )
        {
					cmp[n].set_Velocity( cmp[n].get_Massflow() / cmp[n].area );  // attenltion! Velocity and MassFlow use same variable
				}
				else
        {
					cmp[n].set_Velocity( cmp[n].get_Massflow()*RefVelocity*RefLength*RefLength / cmp[n].area );
				}
			}
      
			// 流量指定のときのみ，ca[]に有次元速度パラメータを保存  >> 速度指定の場合には，parseBC::getXML_IBC_SpecVel()で設定済み
      if ( cmp[n].isPolicy_Massflow() )
      {
        if ( Unit.Param != DIMENSIONAL )
        {
          Hostonly_ stamped_printf("Error: Non-dimensional condition\n");
          Exit(0);
        }
        else
        {
          cmp[n].ca[CompoList::amplitude] = cmp[n].get_Velocity();
          cmp[n].ca[CompoList::bias] = cmp[n].ca[CompoList::bias] / cmp[n].area; // dimensional velocity
        }
      }
		}
  }
	
  // 発熱密度の計算(有次元) -- 発熱量と発熱密度
  REAL_TYPE a, vol;
  a = deltaX*RefLength;
  vol = a*a*a;
  
  for (int n=1; n<=NoCompo; n++)
  {
    if ( cmp[n].getType()==HEAT_SRC )
    {
      if (cmp[n].get_sw_Heatgen() == CompoList::hsrc_watt)
      {
        cmp[n].set_HeatDensity( cmp[n].get_HeatValue() / ((REAL_TYPE)cmp[n].getElement()*vol) );
      }
      else // 発熱密度
      {
        cmp[n].set_HeatValue( cmp[n].get_HeatDensity() * ((REAL_TYPE)cmp[n].getElement()*vol) );
      }
    }
  }
  
  // Darcy係数（単媒質）
  // C[0-2]; 有次元，C[3-5]; 無次元係数
  REAL_TYPE ki;
  for (int n=1; n<=NoCompo; n++)
  {
    if ( cmp[n].getType()==DARCY )
    {
      ki = (mat[m].P[p_viscosity]*RefLength) / (mat[n].P[p_density]*RefVelocity);
      cmp[n].ca[3] = ki / cmp[n].ca[0];
      cmp[n].ca[4] = ki / cmp[n].ca[1];
      cmp[n].ca[5] = ki / cmp[n].ca[2];
    }    
  }
  
  // Pressure Loss
  REAL_TYPE DensityOfMedium, cf[6];
  for (int n=1; n<=NoCompo; n++) {
    
    if ( cmp[n].getType()==HEX )
    {
      for (int i=0; i<6; i++) cf[i] = cmp[n].ca[i];
      
      // 流量と圧力損失量計算の有次元化の係数
      cmp[n].set_CoefMassflow( RefLength * RefLength * RefVelocity );
      cmp[n].set_CoefPrsLoss( cf[5] * RefDensity * RefVelocity * RefVelocity / RefLength );
      
      // Normalize
      if ( cmp[n].getPrsUnit() == CompoList::unit_mmAq )
      {
        // Water: T=300K, p=101.325kPa > 996.62 kg/m^3
        DensityOfMedium = 996.62;
        convertHexCoef(cf, DensityOfMedium);
      }
      else if ( cmp[n].getPrsUnit() == CompoList::unit_mmHg )
      {
        // Mercury: T=300K > 13538 kg/m^3
        DensityOfMedium = 13538.0;
        convertHexCoef(cf, DensityOfMedium);
      }
      else if ( cmp[n].getPrsUnit() == CompoList::unit_Pa )
      {
        convertHexCoef(cf);
      }
      else if ( cmp[n].getPrsUnit() == CompoList::unit_NonDimensional )
      {
        cf[4] /= RefVelocity; // 無次元の場合には単位変更のみ
        cf[5] *= (1e-3/RefLength);
      }
      
      for (int i=0; i<6; i++) cmp[n].ca[i] = cf[i];
    }
    
  }
  
  // 外部境界面の速度の指定パラメータを有次元化
  if ( Unit.Param == NONDIMENSIONAL )
  {
    for (int n=0; n<NOFACE; n++)
    {
      switch ( BO[n].getClass() )
      {
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
  if ( Unit.Param == NONDIMENSIONAL )
  {
    for (int n=0; n<NOFACE; n++) {
      
      switch ( BO[n].getClass() )
      {
        case OBC_OUTFLOW:
        case OBC_TRC_FREE:
        case OBC_FAR_FIELD:
          if ( BO[n].get_pType() == P_DIRICHLET )
          {
            BO[n].p = FBUtility::convND2D_P(BO[n].p, BasePrs, RefDensity, RefVelocity, Unit.Prs); 
          }          
          break;
          
        case OBC_PERIODIC:
          if ( BO[n].get_PrdcMode() != BoundaryOuter::prdc_Simple ) // Dirichlet or Bidirectionalを指定の場合
          {
            BO[n].p = FBUtility::convND2D_P(BO[n].p, BasePrs, RefDensity, RefVelocity, Unit.Prs);
          }
          break;
          
        default:
          break;
      }
    }
  }
  
  // 初期条件の値の有次元化
  if ( Unit.Param == NONDIMENSIONAL )
  {
    iv.Pressure = FBUtility::convND2D_P(iv.Pressure, BasePrs, RefDensity, RefVelocity, Unit.Prs);
    iv.Density *= RefDensity;
		iv.VecU    *= RefVelocity;
		iv.VecV    *= RefVelocity;
		iv.VecW    *= RefVelocity;
    
    if ( isHeatProblem() )
    {
			iv.Temperature = FBUtility::convND2Kelvin(iv.Temperature, BaseTemp, DiffTemp); //内部表現をKelvinに
    }
	}
	else
  {
		if ( isHeatProblem() )
    {
			iv.Temperature = FBUtility::convTemp2K(iv.Temperature, Unit.Temp);
    }
	}
}




/* ################################################################# 20130611
// PLOT3Dファイル入出力に関するパラメータ
void Control::get_PLOT3D(FileIO_PLOT3D_READ*  FP3DR, FileIO_PLOT3D_WRITE* FP3DW)
{
  string str;
  string label;
  
  // FileNameGrid --- option
  label = "/Steer/Plot3dOptions/FileNameGrid";
  
  if ( !(tpCntl->GetValue(label, &str)) )
  {
    P3Op.basename_g = "PLOT3DoutputGrid";
  }
  else
  {
    P3Op.basename_g = str;
  }
  if ( P3Op.basename_g.empty() )
  {
    P3Op.basename_g = "PLOT3DoutputGrid";
  }
  
  // FileNameFunc --- option
  label = "/Steer/Plot3dOptions/FileNameFunc";
  
  if ( !(tpCntl->GetValue(label, &str)) )
  {
    P3Op.basename_f = "PLOT3Doutput";
  }
  else
  {
    P3Op.basename_f = str;
  }
  if ( P3Op.basename_f.empty() )
  {
    P3Op.basename_f = "PLOT3Doutput";
  }
  
  // GridKind
  /*
   label = "/Steer/Plot3dOptions/GridKind";
   
   if ( !(tpCntl->GetValue(label, &str )) )
   {
   Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
   Exit(0);
   }
   else
   {
   if     ( !strcasecmp(str.c_str(), "SingleGrid") ) FP3DR->setSingleGrid();
   else if( !strcasecmp(str.c_str(), "MultiGrid") )  FP3DR->setMultiGrid();
   else
   {
   Hostonly_ stamped_printf("\tInvalid keyword is described for '%s'\n", label.c_str());
   Exit(0);
   }
   }
   *
  
  FP3DR->setMultiGrid(); // 常にmulti grid
  
  
  // 格子の移動
  label = "/Steer/Plot3dOptions/GridMobility";
  
  if ( !(tpCntl->GetValue(label, &str )) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    Exit(0);
  }
  else
  {
    if     ( !strcasecmp(str.c_str(), "movable") )   FP3DR->setMoveGrid(GRID_MOVE);
    else if( !strcasecmp(str.c_str(), "immovable") ) FP3DR->setMoveGrid(GRID_NOT_MOVE);
    else
    {
      Hostonly_ stamped_printf("\tInvalid keyword is described for '%s'\n", label.c_str());
      Exit(0);
    }
  }
  
  
  // 時間方向の変化
  label = "/Steer/Plot3dOptions/StateOfTime";
  
  if ( !(tpCntl->GetValue(label, &str )) )
  {
    FP3DR->setSteady(FB_UNSTEADY);
  }
  else
  {
    if     ( !strcasecmp(str.c_str(), "steady") )   FP3DR->setSteady(FB_STEADY);
    else if( !strcasecmp(str.c_str(), "unsteady") ) FP3DR->setSteady(FB_UNSTEADY);
    else
    {
      Hostonly_ stamped_printf("\tInvalid keyword is described for '%s'\n", label.c_str());
      Exit(0);
    }
  }
  
  
  // IBLANKファイル
  label = "/Steer/Plot3dOptions/SetIblankFlag";
  
  if ( !(tpCntl->GetValue(label, &str )) )
  {
    FP3DR->setIBlankFlag(SET_IBLANK);
  }
  else
  {
    if     ( !strcasecmp(str.c_str(), "on") )  FP3DR->setIBlankFlag(SET_IBLANK);
    else if( !strcasecmp(str.c_str(), "off") ) FP3DR->setIBlankFlag(NOT_SET_IBLANK);
    else
    {
      Hostonly_ stamped_printf("\tInvalid keyword is described for '%s'\n", label.c_str());
      Exit(0);
    }
  }
  
  
  // 次元数
  /*
   label = "/Steer/Plot3dOptions/Dimension";
   if ( !(tpCntl->GetValue(label, &str )) )
   {
   //Hostonly_ stamped_printf("\tParsing error : fail to get '/Steer/plot3doptions/dimension'\n");
   //Exit(0);
   FP3DR->setDimension3D();
   }
   else
   {
   if     ( !strcasecmp(str.c_str(), "2d") ) FP3DR->setDimension2D();
   else if( !strcasecmp(str.c_str(), "3d") ) FP3DR->setDimension3D();
   else
   {
   Hostonly_ stamped_printf("\tInvalid keyword is described for '%s'\n", label.c_str());
   Exit(0);
   }
   }
   *
  FP3DR->setDimension3D(); // 常に三次元
  
  // FormatType
  label = "/Steer/Plot3dOptions/FormatType";
  if ( !(tpCntl->GetValue(label, &str )) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    Exit(0);
  }
  else
  {
    if     ( !strcasecmp(str.c_str(), "unformatted") )  FP3DR->setFormat(UNFORMATTED);
    else if( !strcasecmp(str.c_str(), "formatted") )    FP3DR->setFormat(FORMATTED);
    else if( !strcasecmp(str.c_str(), "binary") )       FP3DR->setFormat(C_BINARY);
    else
    {
      Hostonly_ stamped_printf("\tInvalid keyword is described for '%s'\n", label.c_str());
      Exit(0);
    }
  }
  
  
  // 出力の単精度or倍精度指定
  label = "/Steer/Plot3dOptions/RealType";
  if ( !(tpCntl->GetValue(label, &str )) )
  {
    //Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    //Exit(0);
    int d_type = (sizeof(REAL_TYPE) == 4) ? 1 : 2;  // 1-float / 2-double
    FP3DR->setRealType(d_type);
  }
  else
  {
    if     ( !strcasecmp(str.c_str(), "float") ) FP3DR->setRealType(OUTPUT_FLOAT);
    else if( !strcasecmp(str.c_str(), "double") ) FP3DR->setRealType(OUTPUT_DOUBLE);
    else
    {
      Hostonly_ stamped_printf("\tInvalid keyword is described for '%s'\n", label.c_str());
      Exit(0);
    }
  }
  
  
  //copy options read to write
  FP3DW->setGridKind(FP3DR->IsGridKind());
  FP3DW->setMoveGrid(FP3DR->IsMoveGrid());
  FP3DW->setSteady(FP3DR->IsSteady());
  FP3DW->setIBlankFlag(FP3DR->IsIBlankFlag());
  FP3DW->setDim(FP3DR->GetDim());
  FP3DW->setFormat(FP3DR->GetFormat());
  FP3DW->setRealType(FP3DR->GetRealType());
  
  
  // OutputXyz
  label = "/Steer/Plot3dOptions/OutputXyz";
  
  if ( !(tpCntl->GetValue(label, &str )) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    Exit(0);
  }
  else
  {
    if     ( !strcasecmp(str.c_str(), "on") )  P3Op.IS_xyz = ON;
    else if( !strcasecmp(str.c_str(), "off") ) P3Op.IS_xyz = OFF;
    else
    {
      Hostonly_ stamped_printf("\tInvalid keyword is described for '%s'\n", label.c_str());
      Exit(0);
    }
  }
  
  
  // OutputQ
  
  /*
   label = "/Steer/plot3doptions/OutputQ";
   
   if ( !(tpCntl->GetValue(label, &str )) )
   {
   //Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
   //Exit(0);
   P3Op.IS_q = OFF;
   }
   else
   {
   if     ( !strcasecmp(str.c_str(), "on") )  P3Op.IS_q = ON;
   else if( !strcasecmp(str.c_str(), "off") ) P3Op.IS_q = OFF;
   else
   {
   Hostonly_ stamped_printf("\tInvalid keyword is described for '%s'\n", label.c_str());
   Exit(0);
   }
   }*
  
  P3Op.IS_q = OFF; // 常にoff
  
  
  // OutputFunction
  label = "/Steer/Plot3dOptions/OutputFunction";
  
  if ( !(tpCntl->GetValue(label, &str )) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    Exit(0);
  }
  else
  {
    if     ( !strcasecmp(str.c_str(), "on") )  P3Op.IS_funciton = ON;
    else if( !strcasecmp(str.c_str(), "off") ) P3Op.IS_funciton = OFF;
    else
    {
      Hostonly_ stamped_printf("\tInvalid keyword is described for '%s'\n", label.c_str());
      Exit(0);
    }
  }
  
  
  // OutputFuncName
  label = "/Steer/Plot3dOptions/OutputFuncName";
  
  if ( !(tpCntl->GetValue(label, &str )) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    Exit(0);
  }
  else
  {
    if     ( !strcasecmp(str.c_str(), "on") )  P3Op.IS_function_name = ON;
    else if( !strcasecmp(str.c_str(), "off") ) P3Op.IS_function_name = OFF;
    else
    {
      Hostonly_ stamped_printf("\tInvalid keyword is described for '%s'\n", label.c_str());
      Exit(0);
    }
  }
  
  
  // OutputFvbnd
  
  /*
   label = "/Steer/plot3doptions/OutputFvbnd";
   
   if ( !(tpCntl->GetValue(label, &str )) )
   {
   //Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
   //Exit(0);
   P3Op.IS_fvbnd = OFF;
   }
   else
   {
   if     ( !strcasecmp(str.c_str(), "on") )  P3Op.IS_fvbnd = ON;
   else if( !strcasecmp(str.c_str(), "off") ) P3Op.IS_fvbnd = OFF;
   else
   {
   Hostonly_ stamped_printf("\tInvalid keyword is described for '%s'\n", label.c_str());
   Exit(0);
   }
   }*
  
  P3Op.IS_fvbnd = OFF; // 常にoff
  
  // DivideFunc ---> 出力を項目別にファイル分割するオプション
  label = "/Steer/Plot3dOptions/DivideFunc";
  
  if ( !(tpCntl->GetValue(label, &str )) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    Exit(0);
  }
  else
  {
    if     ( !strcasecmp(str.c_str(), "on") )  P3Op.IS_DivideFunc = ON;
    else if( !strcasecmp(str.c_str(), "off") ) P3Op.IS_DivideFunc = OFF;
    else
    {
      Hostonly_ stamped_printf("\tInvalid keyword is described for '%s'\n", label.c_str());
      Exit(0);
    }
  }
  
}
*/
