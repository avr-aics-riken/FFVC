//##################################################################################
//
// Flow Base class
//
// Copyright (c) 2007-2011 VCAD System Research Program, RIKEN.
// All rights reserved.
//
// Copyright (c) 2011-2014 Institute of Industrial Science, The University of Tokyo.
// All rights reserved.
//
// Copyright (c) 2012-2014 Advanced Institute for Computational Science, RIKEN.
// All rights reserved.
//
//##################################################################################

/**
 @file   Control.C
 @brief  FlowBase Control class
 @author aics
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
    case CONJUGATE_HT:
    case CONJUGATE_HT_NATURAL:
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
        case CONJUGATE_HT:
        case CONJUGATE_HT_NATURAL:
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
bool DTcntl::setScheme(const char* str, const double val)
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
void Control::displayParams(FILE* mp, FILE* fp,
                            IterationCtl* IC,
                            DTcntl* DT,
                            ReferenceFrame* RF,
                            MediumList* mat,
                            CompoList* cmp,
                            const int em)
{
  printSteerConditions(mp, IC, DT, RF, em);
  printSteerConditions(fp, IC, DT, RF, em);
  printParaConditions(mp, mat);
  printParaConditions(fp, mat);
  printInitValues(mp, cmp);
  printInitValues(fp, cmp);
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
 * @note 他のパラメータ取得に先んじて処理しておくもの
 */
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
  
  
  // ファイル入出力に関するパラメータ
  getFieldData();
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

  
  // 圧力ノイマン条件のタイプ >> get_Log()よりも先に
  getWallType();
  

  getLog();
  
  
  getTurbulenceModel();
  
  
  getStartCondition();
  
  
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
  
  if ( !(tpCntl->getInspectedValue(label, str )) )
  {
	  Exit(0);
  }
  OperatorName = str;
  
  
  // パラメータチェックフラグ (Hidden)
  label = "/ApplicationControl/CheckParameter";
  
  if ( tpCntl->chkLabel(label) )
  {
    if ( tpCntl->getInspectedValue(label, str) )
    {
      if     ( !strcasecmp(str.c_str(), "On") )   CheckParam = ON;
      else if( !strcasecmp(str.c_str(), "Off") )  CheckParam = OFF;
      else
      {
        Hostonly_ printf("\tInvalid keyword is described for '%s'\n", label.c_str());
        Exit(0);
      }
    }
    else
    {
      Exit(0);
    }
    
  }
  
  
  // 変数の範囲制限モードを取得 (Hidden)
  Hide.Range_Limit = Range_Normal;
  
  label = "/ApplicationControl/VariableRange";
  
  if ( tpCntl->chkLabel(label) )
  {
    if ( tpCntl->getInspectedValue(label, str) )
    {
      if     ( !strcasecmp(str.c_str(), "on") ) Hide.Range_Limit = Range_Normal;
      else if( !strcasecmp(str.c_str(), "off") ) Hide.Range_Limit = Range_Cutoff;
      else
      {
        Hostonly_ stamped_printf("\tInvalid keyword is described for '%s'\n", label.c_str());
        Exit(0);
      }
    }
    else
    {
      Exit(0);
    }
  }
  
  
  
  // デバッグ用のdiv(u)の出力指定 (Hidden)
  FIO.Div_Debug = OFF;
  label = "/ApplicationControl/DebugDivergence";
  
  if ( tpCntl->chkLabel(label) )
  {
    if ( tpCntl->getInspectedValue(label, str) )
    {
      if     ( !strcasecmp(str.c_str(), "on") )    FIO.Div_Debug = ON;
      else if( !strcasecmp(str.c_str(), "off") )   FIO.Div_Debug = OFF;
      else
      {
        Hostonly_ stamped_printf("\tInvalid keyword is described for '%s'\n", label.c_str());
        Exit(0);
      }
    }
    else
    {
      Exit(0);
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
  string str;
  string label;
  
  // scheme
  label="/ConvectionTerm/scheme";
  
  if ( !(tpCntl->getInspectedValue(label, str)) )
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
    
		if ( !(tpCntl->getInspectedValue(label, str )) )
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
// @brief 無次元パラメータを各種モードに応じて設定する
// @note 純強制対流　有次元　（代表長さ，代表速度，動粘性係数，温度拡散係数）
//                 無次元　（Pr, Re > RefV=RefL=1）
// @see bool Control::setParameters(MediumList* mat, CompoList* cmp)
void Control::getDimensionlessParameter()
{
  REAL_TYPE ct = 0.0;
  string label;
  
  label="/Reference/Reynolds";
  
  if ( !(tpCntl->getInspectedValue(label, ct )) )
  {
    Hostonly_ printf("\tParsing error in '%s'\n", label.c_str());
	  Exit(0);
  }
  Reynolds=(REAL_TYPE)ct;
  
  
  label="/Reference/Prandtl";
  
  if ( !(tpCntl->getInspectedValue(label, ct )) )
  {
    Hostonly_ printf("\tParsing error in '%s'\n", label.c_str());
	  Exit(0);
  }
  Prandtl=(REAL_TYPE)ct;
}



// #################################################################
// 計算モデルのドライバー情報を取得
void Control::getDriver()
{
  string str;
  string label;
  REAL_TYPE ct;
  
  // ドライバの設定 値が正の値のとき，有効．ゼロの場合はドライバなし
  label = "/GeometryModel/Driver/Length";
  if ( tpCntl->getInspectedValue(label, ct ) )
  {
    drv.length = ( Unit.Param == DIMENSIONAL ) ? ct : ct * RefLength;
  }
  else
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    Exit(0);
  }
  
  if ( drv.length < 0.0 )
  {
    Hostonly_ stamped_printf("\tError : Value of '%s' must be positive.\n", label.c_str());
    Exit(0);
  }
  
  
  // Only driver is specified
  if ( drv.length > 0.0 )
  {
    label = "/GeometryModel/Driver/Medium";
    if ( !(tpCntl->getInspectedValue(label, str )) )
    {
      Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
      Exit(0);
    }
    drv.Label = str;
    
    label = "/GeometryModel/Driver/FaceMedium";
    if ( !(tpCntl->getInspectedValue(label, str )) )
    {
      Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
      Exit(0);
    }
    drv.faceLabel = str;
  }
  
}


// #################################################################
// @brief ファイル入出力に関するパラメータを取得し，sphフォーマットの出力の並列モードを指定する．
// @note インターバルパラメータは，setParameters()で無次元して保持
// @pre getTimeControl()
void Control::getFieldData()
{
  
  REAL_TYPE f_val=0.0;
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
  
  if ( !(tpCntl->getInspectedValue(label, str )) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
	  Exit(0);
  }
  
  if     ( !strcasecmp(str.c_str(), "sph") )     FIO.Format = sph_fmt;
  else if( !strcasecmp(str.c_str(), "bov") )     FIO.Format = bov_fmt;
  else
  {
    Hostonly_ stamped_printf("\tInvalid keyword is described for '%s'\n", label.c_str());
    Exit(0);
  }
  
  
  // インターバル
  label = "/Output/Data/BasicVariables/TemporalType";
  
  if ( !(tpCntl->getInspectedValue(label, str )) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    Exit(0);
  }
  else
  {
    if     ( !strcasecmp(str.c_str(), "step") )
    {
      Interval[tg_basic].setMode(IntervalManager::By_step);
    }
    else if( !strcasecmp(str.c_str(), "time") )
    {
      Interval[tg_basic].setMode(IntervalManager::By_time);
    }
    else
    {
      Hostonly_ stamped_printf("\tParsing error : Invalid keyword for '%s'\n", label.c_str());
      Exit(0);
    }
    
    label="/Output/Data/BasicVariables/Interval";
    
    if ( !(tpCntl->getInspectedValue(label, f_val )) )
    {
      Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
      Exit(0);
    }
    else
    {
      Interval[tg_basic].setInterval((double)f_val);
    }
  }
  
  
  
  // 派生変数
  
  /* ファイルフォーマット >> 基本変数と同じ
  label = "/Output/Data/DerivedVariables/Format";
  
  if ( !(tpCntl->getInspectedValue(label, str )) )
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
  else
  {
    Hostonly_ stamped_printf("\tInvalid keyword is described for '%s'\n", label.c_str());
    Exit(0);
  }*/
  
  switch ( FIO.Format )
  {
    case sph_fmt:
      getFormat_sph();
      break;
      
    case bov_fmt:
      break;
  }
  
  
  // インターバル
  label = "/Output/Data/DerivedVariables/TemporalType";
  
  if ( !(tpCntl->getInspectedValue(label, str )) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    Exit(0);
  }
  else
  {
    if     ( !strcasecmp(str.c_str(), "step") )
    {
      Interval[tg_derived].setMode(IntervalManager::By_step);
    }
    else if( !strcasecmp(str.c_str(), "time") )
    {
      Interval[tg_derived].setMode(IntervalManager::By_time);
    }
    else
    {
      Hostonly_ stamped_printf("\tParsing error : Invalid keyword for '%s'\n", label.c_str());
      Exit(0);
    }
    
    label="/Output/Data/DerivedVariables/Interval";
    
    if ( !(tpCntl->getInspectedValue(label, f_val )) )
    {
      Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
      Exit(0);
    }
    else
    {
      Interval[tg_derived].setInterval((double)f_val);
    }
  }
  
  // 全圧
  label="/Output/Data/DerivedVariables/TotalPressure";
  
  if ( !(tpCntl->getInspectedValue(label, str )) )
  {
	  Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
	  Exit(0);
  }
  
  if     ( !strcasecmp(str.c_str(), "on") )  varState[var_TotalP] = ON;
  else if( !strcasecmp(str.c_str(), "off") ) varState[var_TotalP] = OFF;
  else
  {
    Hostonly_ stamped_printf("\tInvalid keyword is described for '%s'\n", label.c_str());
    Exit(0);
  }
  
  // 渦度ベクトル
  label="/Output/Data/DerivedVariables/Vorticity";
  
  if ( !(tpCntl->getInspectedValue(label, str )) )
  {
	  Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
	  Exit(0);
  }
  
  if     ( !strcasecmp(str.c_str(), "on") )  varState[var_Vorticity] = ON;
  else if( !strcasecmp(str.c_str(), "off") ) varState[var_Vorticity] = OFF;
  else
  {
    Hostonly_ stamped_printf("\tInvalid keyword is described for '%s'\n", label.c_str());
    Exit(0);
  }
  
  // 速度勾配テンソルの第2不変量
  label="/Output/Data/DerivedVariables/Qcriterion";
  
  if ( !(tpCntl->getInspectedValue(label, str )) )
  {
	  Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
	  Exit(0);
  }
  
  if     ( !strcasecmp(str.c_str(), "on") )  varState[var_Qcr] = ON;
  else if( !strcasecmp(str.c_str(), "off") ) varState[var_Qcr] = OFF;
  else
  {
    Hostonly_ stamped_printf("\tInvalid keyword is described for '%s'\n", label.c_str());
    Exit(0);
  }
  
  // ヘリシティ
  label="/Output/Data/DerivedVariables/Helicity";
  
  if ( !(tpCntl->getInspectedValue(label, str )) )
  {
	  Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
	  Exit(0);
  }
  
  if     ( !strcasecmp(str.c_str(), "on") )  varState[var_Helicity] = ON;
  else if( !strcasecmp(str.c_str(), "off") ) varState[var_Helicity] = OFF;
  else
  {
    Hostonly_ stamped_printf("\tInvalid keyword is described for '%s'\n", label.c_str());
    Exit(0);
  }

  
  
  // 平均値操作に関するパラメータを取得
  if ( Mode.Average == ON )
  {
	  label = "/Output/Data/AveragedVariables/TemporalType";
    
	  if ( !(tpCntl->getInspectedValue(label, str )) )
    {
		  Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
		  Exit(0);
	  }
	  else
    {
		  if     ( !strcasecmp(str.c_str(), "step") )
      {
        if ( Interval[tg_average].getMode() == IntervalManager::By_time )
        {
          Hostonly_ stamped_printf("\tError : Specified temporal mode is not consistent with '/TimeControl/Average/TemporalType'\n");
          Exit(0);
        }
		  }
		  else if( !strcasecmp(str.c_str(), "time") )
      {
        if ( Interval[tg_average].getMode() == IntervalManager::By_step )
        {
          Hostonly_ stamped_printf("\tError : Specified temporal mode is not consistent with '/TimeControl/Average/TemporalType'\n");
          Exit(0);
        }
		  }
		  else
      {
			  Hostonly_ stamped_printf("\tParsing error : Invalid keyword for '%s'\n", label.c_str());
			  Exit(0);
		  }
    }
    
    double val;
    label="/Output/Data/AveragedVariables/Interval";
    
    if ( !(tpCntl->getInspectedValue(label, val )) )
    {
      Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
      Exit(0);
    }
    else
    {
      Interval[tg_average].setInterval(val);
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
  
  if ( !(tpCntl->chkNode(label)) )
  {
    Hostonly_ stamped_printf("\tParsing error : Missing '%s'\n", label.c_str());
    Exit(0);
  }
  
  
  // 出力ガイドセルモード
  label = "/Output/FormatOption/SPH/GuideOut";
  
  if ( !(tpCntl->getInspectedValue(label, str)) )
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
  
  if ( !(tpCntl->getInspectedValue(label, str)) )
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
  
  if ( !(tpCntl->getInspectedValue(label, str)) )
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
  
  if ( !(tpCntl->getInspectedValue(label, str )) )
  {
    Hostonly_ stamped_printf("\tParsing error : Invalid char* value in '%s'\n", label.c_str());
    Exit(0);
  }
  
  
  if ( !(tpCntl->getInspectedValue(label, str )) )
  {
    Hostonly_ stamped_printf("\tError : '%s'\n", label.c_str());
    Exit(0);
  }
  
  if     ( FBUtility::compare(str, "ParallelPlate2D") )   Mode.Example = id_PPLT2D;
  else if( FBUtility::compare(str, "Duct") )              Mode.Example = id_Duct;
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
  
  
  // フィルの媒質指定
  label = "/GeometryModel/FillMedium";
  if ( !(tpCntl->getInspectedValue(label, str)) )
  {
    Hostonly_ printf("\tParsing error in '%s'\n", label.c_str());
	  Exit(0);
  }
  FillMedium = str;
  
  
  // 流体セルのフィルの開始面指定
  label = "/GeometryModel/HintOfFillSeedDirection";
  
  if ( !(tpCntl->getInspectedValue(label, str )) )
  {
    Hostonly_ stamped_printf("\tParsing error : Invalid value in '%s'\n", label.c_str());
    Exit(0);
  }
  else
  {
    if     ( !strcasecmp(str.c_str(), "xminus" ) ) FillSeedDir = X_minus;
    else if( !strcasecmp(str.c_str(), "xplus" ) )  FillSeedDir = X_plus;
    else if( !strcasecmp(str.c_str(), "yminus" ) ) FillSeedDir = Y_minus;
    else if( !strcasecmp(str.c_str(), "yplus" ) )  FillSeedDir = Y_plus;
    else if( !strcasecmp(str.c_str(), "zminus" ) ) FillSeedDir = Z_minus;
    else if( !strcasecmp(str.c_str(), "zplus" ) )  FillSeedDir = Z_plus;
    else
    {
      FillSeedDir = X_minus;
      Hostonly_ printf("\tDefault 'X_minus' is set for Hint Of FillSeed direction\n");
    }
  }
  
  
  // ヒントに使うフィルの媒質指定
  label = "/GeometryModel/HintOfFillSeedMedium";
  if ( !(tpCntl->getInspectedValue(label, str)) )
  {
    Hostonly_ printf("\tParsing error in '%s'\n", label.c_str());
	  Exit(0);
  }
  SeedMedium = str;
  
  
  // フィル方向制御 (NOT mandatory)
  string dir[3];
  label = "/GeometryModel/FillDirectionControl";
  if ( tpCntl->chkLabel(label) )
  {
    if ( !(tpCntl->getInspectedVector(label, dir, 3)) )
    {
      Hostonly_ printf("\tParsing error in '%s'\n", label.c_str());
      Exit(0);
    }
    else
    {
      FillSuppress[0] = ( !strcasecmp(dir[0].c_str(), "fill" ) ) ? 1 : 0;
      FillSuppress[1] = ( !strcasecmp(dir[1].c_str(), "fill" ) ) ? 1 : 0;
      FillSuppress[2] = ( !strcasecmp(dir[2].c_str(), "fill" ) ) ? 1 : 0;
    }
  }

  
  // Geometry output (NOT mandatory)
  label = "/GeometryModel/Output";
  
  if ( tpCntl->chkLabel(label) )
  {
    if ( !(tpCntl->getInspectedValue(label, str )) )
    {
      Hostonly_ stamped_printf("\tError : '%s'\n", label.c_str());
      Exit(0);
    }
    else
    {
      if ( !strcasecmp(str.c_str(), "on") ) Hide.GeomOutput = ON;
    }
  }
  
  // Glyph output (NOT mandatory)
  label = "/GeometryModel/OutputGlyph";
  
  if ( tpCntl->chkLabel(label) )
  {
    if ( !(tpCntl->getInspectedValue(label, str )) )
    {
      Hostonly_ stamped_printf("\tError : '%s'\n", label.c_str());
      Exit(0);
    }
    else
    {
      if ( !strcasecmp(str.c_str(), "on") ) Hide.GlyphOutput = ON;
      if ( !strcasecmp(str.c_str(), "InnerOnly") ) Hide.GlyphOutput = 2; // special treatment
    }
  }
  
  
  // ボクセルファイル出力 (Hidden)
  FIO.IO_Voxel = OFF;
  label = "/GeometryModel/VoxelOutput";
  
  if ( tpCntl->chkLabel(label) )
  {
    if ( tpCntl->getInspectedValue(label, str) )
    {
      if     ( !strcasecmp(str.c_str(), "svx") )  FIO.IO_Voxel = Sphere_SVX;
      else if( !strcasecmp(str.c_str(), "off") )  FIO.IO_Voxel = OFF;
      else
      {
        Hostonly_ stamped_printf("\tInvalid keyword is described for '%s'\n", label.c_str());
        Exit(0);
      }
    }
    else
    {
      Exit(0);
    }
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
  if( !(tpCntl->chkNode(base)) )
  {
    Hostonly_ stamped_printf("\tParsing error : Missing the section of '%s'\n", base.c_str());
    Exit(0);
  }

  // タグ内のラベル数をチェック 0;
  int nnode = 0;
  
  if ( (nnode = tpCntl->countLabels(base)) == 0 )
  {
    Hostonly_ stamped_printf("\tNo labels inside /Iteration\n");
    return;
  }

  
  // LinearSolverの個数を取得
  int counter=0;
  for (int i=1; i<=nnode; i++)
  {
    if ( !(tpCntl->getNodeStr(base, i, str)) )
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
    if ( !(tpCntl->getNodeStr(base, i+1, str)) )
    {
      Hostonly_ printf("\tParsing error : Missing 'LinearSolver'\n");
      Exit(0);
    }

    if( strcasecmp(str.substr(0,12).c_str(), "LinearSolver") ) continue;
    
    
    // alias 
    leaf = base + "/" + str;
    label = leaf + "/Alias";
    
    if ( !(tpCntl->getInspectedValue(label, str )) )
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
    if ( !(tpCntl->getInspectedValue(label, str )) )
    {
      Hostonly_ printf("\tParsing error : No '%s'\n", label.c_str());
      Exit(0);
    }
    
    if ( !Criteria[i].setLS(str) )
    {
      Hostonly_ stamped_printf("\tInvalid keyword is described for Linear_Solver\n");
      Exit(0);
    }

    
    label = leaf + "/MaxIteration";
    if ( !(tpCntl->getInspectedValue(label, i_val )) )
    {
      Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
      Exit(0);
    }
    Criteria[i].setMaxIteration(i_val);


    label = leaf + "/ConvergenceCriterion";
    if ( !(tpCntl->getInspectedValue(label, f_val )) )
    {
      Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
      Exit(0);
    }
    Criteria[i].setCriterion(f_val);

    
    label = leaf + "/NormType";
    if ( !(tpCntl->getInspectedValue(label, str )) )
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
          Mode.Log_Itr = ON;
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
    if ( !Criteria[i].getInherentPara(tpCntl, leaf, ExperimentNaive) )
    {
      Hostonly_ printf("\tError : Invalid Linear Solver[%d]\n", ls);
      Exit(0);
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
  
  if ( !(tpCntl->getInspectedValue(label, str )) )
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
  
  if ( !(tpCntl->getInspectedValue(label, str )) )
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
    
	  if ( !(tpCntl->getInspectedValue(label, str )) )
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
  
  if ( !(tpCntl->getInspectedValue(label, str )) )
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
  label="/Output/Log/Console/TemporalType";
  
  if ( !(tpCntl->getInspectedValue(label, str )) )
  {
	  Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
	  Exit(0);
  }
  else
  {
	  if     ( !strcasecmp(str.c_str(), "step") )
    {
      Interval[tg_console].setMode(IntervalManager::By_step);
	  }
	  else if( !strcasecmp(str.c_str(), "time") )
    {
      Interval[tg_console].setMode(IntervalManager::By_time);
	  }
	  else
    {
		  Hostonly_ stamped_printf("\tParsing error : Invalid keyword for '%s'\n", label.c_str());
		  Exit(0);
	  }
	  
	  label="/Output/Log/Console/Interval";
    
	  if ( !(tpCntl->getInspectedValue(label, f_val )) )
    {
		  Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
		  Exit(0);
	  }
	  else
    {
		  Interval[tg_console].setInterval((double)f_val);
	  }
  }
  
  // Interval file_history
  label="/Output/Log/History/TemporalType";
  
  if ( !(tpCntl->getInspectedValue(label, str )) )
  {
	  Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
	  Exit(0);
  }
  else
  {
	  if     ( !strcasecmp(str.c_str(), "step") )
    {
		  Interval[tg_history].setMode(IntervalManager::By_step);
	  }
	  else if( !strcasecmp(str.c_str(), "time") )
    {
		  Interval[tg_history].setMode(IntervalManager::By_time);
	  }
	  else
    {
		  Hostonly_ stamped_printf("\tParsing error : Invalid keyword for '%s'\n", label.c_str());
		  Exit(0);
	  }
    
	  label="/Output/Log/History/Interval";
    
	  if ( !(tpCntl->getInspectedValue(label, f_val )) )
    {
		  Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
		  Exit(0);
	  }
	  else
    {
		  Interval[tg_history].setInterval((double)f_val);
	  }
  }
  
  // CCNV file
  label="/Output/Log/CCNVfile";
  
  if ( !(tpCntl->getInspectedValue(label, str )) )
  {
    ; // not mandatory
  }
  else{
    if     ( !strcasecmp(str.c_str(), "on") )   Mode.CCNV = ON;
    else if( !strcasecmp(str.c_str(), "off") )  Mode.CCNV = OFF;
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
 * @brief 参照パラメータを取得
 * @note Ref_IDで指定される媒質を代表物性値とする
 */
void Control::getReference()
{
  REAL_TYPE ct2;
  string label, str;
  
  
  label = "/Reference/Length";
  if ( !(tpCntl->getInspectedValue(label, ct2 )) )
  {
    Hostonly_ printf("\tParsing error in '%s'\n", label.c_str());
    Exit(0);
  }
  RefLength = ct2;
  
  
  label = "/Reference/Velocity";
  if ( !(tpCntl->getInspectedValue(label, ct2 )) )
  {
    Hostonly_ printf("\tParsing error in '%s'\n", label.c_str());
    Exit(0);
  }
  RefVelocity = ct2;
  
  /*
  label = "/Reference/MassDensity";
  if ( !(tpCntl->getInspectedValue(label, ct2 )) )
  {
    Hostonly_ printf("\tParsing error in '%s'\n", label.c_str());
    Exit(0);
  }
  RefDensity = ct2;
   */
  
  
  label = "/Reference/BasePressure";
  if ( !(tpCntl->getInspectedValue(label, ct2 )) )
  {
    Hostonly_ printf("\tParsing error in '%s'\n", label.c_str());
    Exit(0);
  }
  BasePrs = ct2;

  
  label = "/Reference/Medium";
  if ( !(tpCntl->getInspectedValue(label, str)) )
  {
    Hostonly_ printf("\tParsing error in '%s'\n", label.c_str());
	  Exit(0);
  }
  RefMedium = str;
  
  
  if ( isHeatProblem() )
  {
    REAL_TYPE Base, Diff;
    
    label="/Reference/Temperature/Base";
    
    if ( !(tpCntl->getInspectedValue(label, Base )) )
    {
      Hostonly_ stamped_printf("\tParsing error : Invalid float value for '%s'\n", label.c_str());
      Exit(0);
    }
    
    
    label="/Reference/Temperature/Difference";
    
    if ( !(tpCntl->getInspectedValue(label, Diff )) )
    {
      Hostonly_ stamped_printf("\tParsing error : Invalid float value for '%s'\n", label.c_str());
      Exit(0);
    }
    BaseTemp = Base;
    
    if ( Diff < 0.0f )
    {
      Hostonly_ stamped_printf("\tTemperature difference must be positive.\n");
      Exit(0);
    }
    DiffTemp = Diff;
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
  
  if ( !(tpCntl->getInspectedValue(label, str )) )
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
    
    if( tpCntl->getInspectedVector(label, xyz, 3) )
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
  
  if ( !(tpCntl->getInspectedValue(label, str )) )
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
  
  if ( !(tpCntl->getInspectedValue(label, str )) )
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
  
  if ( !(tpCntl->getInspectedValue(label, str )) ) {
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
  
  if ( !(tpCntl->getInspectedValue(label, str )) )
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
  label="/ConvectionTerm/scheme";
  
  if ( !(tpCntl->getInspectedValue(label, str )) )
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
    
		if ( !(tpCntl->getInspectedValue(label, str )) )
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
  
  
  // ソルバーの種類（FLOW_ONLY / THERMAL_FLOW / THERMAL_FLOW_NATURAL / CONJUGATE_HT / CONJUGATE_HT_NATURAL / SOLID_CONDUCTION）と浮力モード
  label="/GoverningEquation/HeatEquation";
  
  if ( !(tpCntl->getInspectedValue(label, str )) )
  {
	  Hostonly_ stamped_printf("\tParsing error : Invalid char* value for '%s'\n", label.c_str());
	  Exit(0);
  }
  
  if     ( !strcasecmp(str.c_str(), "FlowOnly" ) )
  {
    KindOfSolver = FLOW_ONLY;
    varState[var_Velocity]    = true;
    varState[var_Fvelocity]   = true;
    varState[var_Pressure]    = true;
  }
  else if( !strcasecmp(str.c_str(), "ThermalFlow" ) )
  {
    KindOfSolver = THERMAL_FLOW;
    varState[var_Velocity]    = true;
    varState[var_Fvelocity]   = true;
    varState[var_Pressure]    = true;
    varState[var_Temperature] = true;
  }
  else if( !strcasecmp(str.c_str(), "ThermalFlowNatural" ) )
  {
    KindOfSolver = THERMAL_FLOW_NATURAL;
    varState[var_Velocity]    = true;
    varState[var_Fvelocity]   = true;
    varState[var_Pressure]    = true;
    varState[var_Temperature] = true;
  }
  else if( !strcasecmp(str.c_str(), "ConjugateHeatTransfer" ) )
  {
    KindOfSolver = CONJUGATE_HT;
    varState[var_Velocity]    = true;
    varState[var_Fvelocity]   = true;
    varState[var_Pressure]    = true;
    varState[var_Temperature] = true;
  }
  else if( !strcasecmp(str.c_str(), "ConjugateHeatTransferNatural" ) )
  {
    KindOfSolver = CONJUGATE_HT_NATURAL;
    varState[var_Velocity]    = true;
    varState[var_Fvelocity]   = true;
    varState[var_Pressure]    = true;
    varState[var_Temperature] = true;
  }
  else if( !strcasecmp(str.c_str(), "SolidConduction" ) )
  {
    KindOfSolver = SOLID_CONDUCTION;
    varState[var_Temperature] = true;
  }
  else
  {
    Hostonly_ stamped_printf("\tInvalid keyword is described for '%s'\n", label.c_str());
    Exit(0);
  }
  
  // Buoyancy option
  if ( (KindOfSolver==THERMAL_FLOW) ||
       (KindOfSolver==THERMAL_FLOW_NATURAL) ||
       (KindOfSolver==CONJUGATE_HT) ||
       (KindOfSolver==CONJUGATE_HT_NATURAL) )
  {
    label="/GoverningEquation/Buoyancy";
    
    if ( !(tpCntl->getInspectedValue(label, str )) )
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
 
  // 熱問題の解法
  if ( isHeatProblem() )
  {
	  label = "/SolvingMethod/Heat";
    
	  if ( !(tpCntl->getInspectedValue(label, str )) )
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
// @brief 初期値とリスタート条件
// @todo セルフェイスの粗格子リスタート  >> 近似なのでサボる？
// @ see getTimeControl()
void Control::getStartCondition()
{
  string str;
  string label, leaf;
  
  
  // Staging option
  label="/StartCondition/Restart/Staging";
  
  if ( !(tpCntl->getInspectedValue(label, str )) )
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
    
    if ( tpCntl->getInspectedValue(label, str ) )
    {
      f_dfi_in_prs = str.c_str();
    }
    if ( f_dfi_in_prs.empty() == true ) f_dfi_in_prs = "prs";
    
    
    label="/StartCondition/Restart/DFIfiles/Velocity";
    
    if ( tpCntl->getInspectedValue(label, str ) )
    {
      f_dfi_in_vel = str.c_str();
    }
    if ( f_dfi_in_vel.empty() == true ) f_dfi_in_vel = "vel";
    
    
    label="/StartCondition/Restart/DFIfiles/Fvelocity";
    
    if ( tpCntl->getInspectedValue(label, str ) )
    {
      f_dfi_in_fvel = str.c_str();
    }
    if ( f_dfi_in_fvel.empty() == true ) f_dfi_in_fvel = "fvel";
    
    
    if ( isHeatProblem() )
    {
      label="/StartCondition/Restart/DFIfiles/Temperature";
      
      if ( tpCntl->getInspectedValue(label, str ) )
      {
        f_dfi_in_temp = str.c_str();
      }
      if ( f_dfi_in_temp.empty() == true ) f_dfi_in_temp = "tmp";
    }
    
    
    // 平均値
    if ( Mode.Average == ON )
    {
      label="/StartCondition/Restart/DFIfiles/AveragedPressure";
      
      if ( tpCntl->getInspectedValue(label, str ) )
      {
        f_dfi_in_prsa = str.c_str();
      }
      if ( f_dfi_in_prsa.empty() == true ) f_dfi_in_prsa = "prsa";
      
      
      label="/StartCondition/Restart/DFIfiles/AveragedVelocity";
      
      if ( tpCntl->getInspectedValue(label, str ) )
      {
        f_dfi_in_vela = str.c_str();
      }
      if ( f_dfi_in_vela.empty() == true ) f_dfi_in_vela = "vela";
      
      
      if ( isHeatProblem() )
      {
        label="/StartCondition/Restart/DFIfiles/AveragedTemperature";
        
        if ( tpCntl->getInspectedValue(label, str ) )
        {
          f_dfi_in_tempa = str.c_str();
        }
        if ( f_dfi_in_tempa.empty() == true ) f_dfi_in_tempa = "tmpa";
      }
    }
  }
  
  

  // 初期条件 温度はParseBC::getInitTempOfMedium()
  if ( Start == initial_start )
  {
    /* Density
    label="/StartCondition/InitialState/MassDensity";
    
    if ( !(tpCntl->getInspectedValue(label, iv.Density )) )
    {
      Hostonly_ stamped_printf("\tParsing error : Invalid float value for '%s'\n", label.c_str());
      Exit(0);
    }
    */
    
    // Pressure
    label="/StartCondition/InitialState/Pressure";
    
    if ( !(tpCntl->getInspectedValue(label, iv.Pressure )) )
    {
      Hostonly_ stamped_printf("\tParsing error : Invalid float value for '%s'\n", label.c_str());
      Exit(0);
    }
    
    // Velocity
    REAL_TYPE v[3];
    for (int n=0; n<3; n++) v[n]=0.0;
    label="/StartCondition/InitialState/Velocity";
    
    if( !(tpCntl->getInspectedVector(label, v, 3)) )
    {
      Hostonly_ stamped_printf("\tParsing error : fail to get velocity in '%s'\n", label.c_str());
      Exit(0);
    }
    iv.VecU = v[0];
    iv.VecV = v[1];
    iv.VecW = v[2];
  }
  
}


// #################################################################
// 流体の解法アルゴリズムを選択する
// 熱問題の解法は，getSolverProperties()でKOSを決定した後で取得
void Control::getSolvingMethod4Flow()
{
  string str;
  string label;
  
  label = "/SolvingMethod/Flow";
  
  if ( !(tpCntl->getInspectedValue(label, str )) )
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
  
}



// #################################################################
/**
 * @brief 時間制御に関するパラメータを取得する
 * @param [out] DT  DTcntl
 * @note パラメータは，setParameters()で無次元して保持
 */
void Control::getTimeControl(DTcntl* DT)
{
  double ct = 0.0;
  
  string str;
  string label;
  
  // 加速時間
  label = "/TimeControl/Acceleration/TemporalType";
  
  if ( !(tpCntl->getInspectedValue(label, str )) )
  {
	  Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
	  Exit(0);
  }
  else
  {
	  if     ( !strcasecmp(str.c_str(), "step") )
    {
		  Interval[tg_accelra].setMode(IntervalManager::By_step);
	  }
	  else if( !strcasecmp(str.c_str(), "time") )
    {
		  Interval[tg_accelra].setMode(IntervalManager::By_time);
	  }
	  else
    {
		  Hostonly_ stamped_printf("\tParsing error : Invalid keyword for '%s'\n", label.c_str());
		  Exit(0);
	  }
	  
	  label = "/TimeControl/Acceleration/AcceleratingTime";
    
	  if ( !(tpCntl->getInspectedValue(label, ct )) )
    {
		  Hostonly_ stamped_printf("\tParsing error : fail to get '%s\n", label.c_str());
		  Exit(0);
	  }
	  else
    {
		  Interval[tg_accelra].setInterval(ct);
	  }
  }
  
  
  // 時間積分幅を取得する
  label = "/TimeControl/TimeStep/Mode";
  
  if ( !(tpCntl->getInspectedValue(label, str )) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    Exit(0);
  }
  
  
  label = "/TimeControl/TimeStep/DeltaT";
  
  if ( !(tpCntl->getInspectedValue(label, ct )) )
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
  
  if ( !DT->setScheme(str.c_str(), cc) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to set DeltaT\n");
    Exit(0);
  }
  
  // 計算する時間を取得する
  label = "/TimeControl/Session/TemporalType";
  
  if ( !(tpCntl->getInspectedValue(label, str )) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    Exit(0);
  }
  else
  {
    if ( !strcasecmp(str.c_str(), "step") )
    {
      Interval[tg_compute].setMode(IntervalManager::By_step);
    }
    else if ( !strcasecmp(str.c_str(), "time") )
    {
      Interval[tg_compute].setMode(IntervalManager::By_time);
    }
    else
    {
      Hostonly_ stamped_printf("\tParsing error : Invalid keyword for '%s'\n", label.c_str());
      Exit(0);
    }
  }
  
  // スタート
  label = "/TimeControl/Session/Start";
  
  if ( !(tpCntl->getInspectedValue(label, ct )) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    Exit(0);
  }
  double m_start = ct;
  
  
  
  // 終了
  label = "/TimeControl/Session/End";
  
  if ( !(tpCntl->getInspectedValue(label, ct )) )
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
  
  // Intereval Manager への登録 >> 他はFFV::initInterval()で指定
  Interval[tg_compute].setStart(m_start);
  Interval[tg_compute].setLast(m_end);
  Interval[tg_compute].setInterval(m_end-m_start); // tg_computeのインターバルは計算するセッションの長さ
  
  
  
  // 平均値の時刻指定モード
  label = "/TimeControl/Average/TemporalType";
  
  if ( !(tpCntl->getInspectedValue(label, str)) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    Exit(0);
  }
  else
  {
    if     ( !strcasecmp(str.c_str(), "step") )
    {
      Interval[tg_average].setMode(IntervalManager::By_step);
    }
    else if( !strcasecmp(str.c_str(), "time") )
    {
      Interval[tg_average].setMode(IntervalManager::By_time);
    }
    else
    {
      Hostonly_ stamped_printf("\tParsing error : Invalid keyword for '%s'\n", label.c_str());
      Exit(0);
    }
  }
  
  
  // 平均操作開始
  label = "/TimeControl/Average/Start";
  
  if ( !(tpCntl->getInspectedValue(label, ct )) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    Exit(0);
  }
  double avr_start = ct;
  
  Restart_stepAvr = avr_start;
  
  
  // 平均操作終了
  label = "/TimeControl/Average/End";
  
  if ( !(tpCntl->getInspectedValue(label, ct )) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    Exit(0);
  }
  double avr_end = ct;
  
  
  // チェック
  if ( !(avr_start <= avr_end) )
  {
    Hostonly_ stamped_printf("\tError : Average/Start msut be less than Average/End.\n");
    Exit(0);
  }
  
  
  // Intereval Manager への登録 >> 他はFFV::initInterval()で指定
  Interval[tg_average].setStart(avr_start);
  Interval[tg_average].setLast(avr_end);
  

  /* 平均値操作の判断
   
   1    |     2      |   3     << m_start (restart step)
   ---------+------------+-------
        ^            ^
     avr_start    avr_end
   
   case 1 : 平均操作は行うが，まだ指定時刻に到達していないので，平均値ファイルは存在せず，平均値のリスタートはない
        2 : 前セッションから継続して平均操作を行うが，既に平均値ファイルが存在する（はず）ので，平均値のリスタート処理を行う
        3 : 既に平均値操作の区間は終了しているので，平均操作は行わない
   */
  
  if ( Start == initial_start )
  {
    Mode.AverageRestart = OFF; // default
    
    if ( (avr_end > 0.0) && (avr_start >= 0.0) ) // avr_end >= avr_startは既にチェック済み
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
    if ( (avr_start < m_start) && (avr_end < m_start) ) // case 3
    {
      Mode.Average = OFF;
      Mode.AverageRestart = OFF;
    }
    else if ( (avr_start < m_start) && (avr_end > m_start) ) // case 2
    {
      Mode.Average = ON;
      Mode.AverageRestart = ON;
    }
    else if ( (avr_start > m_start) && (avr_end > m_start) ) // case 1
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
  
  if ( !(tpCntl->getInspectedValue(label, str )) )
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
    if ( !(tpCntl->getInspectedValue(label, ct )) )
    {
      Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
      Exit(0);
    }
    LES.Cs = ct;
    
    // damping factor
    label="/TurbulenceModeling/DampingFactor";
    if ( !(tpCntl->getInspectedValue(label, ct )) )
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
  
  if ( !(tpCntl->getInspectedValue(label, str )) )
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
  
  if ( !(tpCntl->getInspectedValue(label, str )) )
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
  
  if ( !(tpCntl->getInspectedValue(label, str )) )
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
  
}


// #################################################################
//境界条件の値(REAL_TYPE型)を取得し，返す
REAL_TYPE Control::getValueReal(const std::string label, TextParser* tpc)
{
  REAL_TYPE df=0.0f;
  
  if ( !(tpc->getInspectedValue(label, df)) ) Exit(0);
  
  return df;
}


// #################################################################
// ベクトル値を取得し，登録する
// normalize = trueのとき，無次元化する
bool Control::getVec(const std::string label, REAL_TYPE* v, TextParser* tpc, bool normalize)
{
  for (int i=0; i<3; i++) v[i]=0.0f;
  
  if( !(tpc->getInspectedVector(label, v, 3)) ) return false;
  
  //単位ベクトル化
  if ( normalize ) UnitVec(v);
  
  return true;
}


// #################################################################
//壁面上の扱いを指定する
void Control::getWallType()
{
  string str;
  string label;
  
  // 圧力のタイプ
  label="/TreatmentOfWall/PressureGradient";
  
  if ( !(tpCntl->getInspectedValue(label, str )) )
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
  
  if ( !(tpCntl->getInspectedValue(label, str )) )
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
void Control::importTP(TextParser* tp)
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
    case X_minus:
    case X_plus:
      r = area / (REAL_TYPE)(m_jmax*m_kmax) * 100.0;
      break;
      
    case Y_minus:
    case Y_plus:
      r = area / (REAL_TYPE)(m_imax*m_kmax) * 100.0;
      break;
      
    case Z_minus:
    case Z_plus:
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
  
  fprintf(fp,"\tFluid  cell inside whole Computational domain = %15ld (%8.4f %%)\n", G_Fcell, (REAL_TYPE)G_Fcell/cell_max *100.0);
  fprintf(fp,"\tActive cell                                   = %15ld (%8.4f %%)\n", G_Acell, (REAL_TYPE)G_Acell/cell_max *100.0);
  
  fprintf(fp,"\n\tFace :      Element (Open ratio)\n");
  for (int i=0; i<NOFACE; i++) {
    fprintf(fp,"\t  %s : %12.0f (%6.2f %%)\n", FBUtility::getDirection(i).c_str(), OpenDomain[i], OpenDomainRatio(i, OpenDomain[i], G_size));
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
 * @param [in] fp  ファイルポインタ
 * @param [in] cmp コンポーネント配列
 * @see void Control::setInitialConditions(void)
 */
void Control::printInitValues(FILE* fp, CompoList* cmp)
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
  
  //fprintf(fp,"\tInitial  MassDensity [kg/m^3]/ [-]   : %12.5e / %12.5e\n", iv.Density, iv.Density/RefDensity);
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
    for (int m=1; m<=NoMedium; m++)
    {
      fprintf(fp,"\tInitial  Temperature [C]     / [-]   : %12.5e / %12.5e\n",
              cmp[m].getInitTemp(),
              FBUtility::convTempD2ND(cmp[m].getInitTemp(), BaseTemp, DiffTemp));
    }
    
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
      if (IC->getNaive()==OFF)
      {
        if ( IC->getBit3()==OFF )
        {
          fprintf(fp,"\t       Linear Solver          :   2-colored SOR SMA (Stride Memory Access, Bit compressed 1-decode)\n");
        }
        else
        {
          fprintf(fp,"\t       Linear Solver          :   2-colored SOR SMA (Stride Memory Access, Bit compressed 3-decodes)\n");
        }
      }
      else
      {
        fprintf(fp,"\t       Linear Solver          :   2-colored SOR SMA (Stride Memory Access, Naive Implementation)\n");
      }
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
  
  fprintf(fp,"\tReference Medium                      :  %s\n", mat[RefMat].getAlias().c_str());
  fprintf(fp,"\n");
  
  // Reference values
  fprintf(fp,"\tRef. Length               [m]         : %12.5e\n", RefLength);
  fprintf(fp,"\tRef. Velocity             [m/s]       : %12.5e\n", RefVelocity);
	fprintf(fp,"\tBase Pressure             [Pa]        : %12.5e\n", BasePrs);
  fprintf(fp,"\tRef. Mass Density         [kg/m^3]    : %12.5e\n", RefDensity);
  fprintf(fp,"\tRef. Specific Heat        [J/(kg K)]  : %12.5e\n", RefSpecificHeat);
  fprintf(fp,"\tRef. Speed of Sound       [m/s]       : %12.5e\n", RefSoundSpeed);
  fprintf(fp,"\tGravity                   [m/s^2]     : %12.5e\n", Gravity);
  fprintf(fp,"\n");
  
  fprintf(fp,"\tSpacing                   [m] / [-]   : %12.5e / %12.5e\n", deltaX*RefLength, deltaX);
  fprintf(fp,"\tTime Scale                [sec]       : %12.5e\n", Tscale);
  fprintf(fp,"\n");
  
  if ( isHeatProblem() )
  {
    fprintf(fp,"\tBase Temperature          [C] / [-]   : %12.5e / %3.1f\n", BaseTemp, 0.0);
    fprintf(fp,"\tTemperature Diff.         [C] / [-]   : %12.5e / %3.1f\n", DiffTemp, 1.0);
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
void Control::printSteerConditions(FILE* fp, IterationCtl* IC, const DTcntl* DT, const ReferenceFrame* RF, const int em)
{
  if ( !fp )
  {
    stamped_printf("\tFail to write into file\n");
    Exit(0);
  }
  
  double dt = DT->get_DT();
  bool  err=true;
  double itm=0.0;
  unsigned stp;
  
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
    case CONJUGATE_HT:
    case CONJUGATE_HT_NATURAL:
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
  
  else if ( (KindOfSolver==CONJUGATE_HT) && (Mode.Buoyancy==NO_BUOYANCY) )
  {
    fprintf(fp,"\t     Kind of Solver           :   Conjugate heat transfer / Forced convection without buoyancy\n");
  }
  else if ( (KindOfSolver==CONJUGATE_HT) && (Mode.Buoyancy==BOUSSINESQ) )
  {
    fprintf(fp,"\t     Kind of Solver           :   Conjugate heat transfer / Forced convection with buoyancy : Boussinesq Approximation\n");
  }
  else if ( (KindOfSolver==CONJUGATE_HT) && (Mode.Buoyancy==LOW_MACH) )
  {
    fprintf(fp,"\t     Kind of Solver           :   Conjugate heat transfer / Forced convection with buoyancy : Low Mach Approximation\n");
  }
  else if ( (KindOfSolver==CONJUGATE_HT_NATURAL) && (Mode.Buoyancy==BOUSSINESQ) )
  {
    fprintf(fp,"\t     Kind of Solver           :   Conjugate heat transfer / Natural convection : Boussinesq Approximation\n");
  }
  else if ( (KindOfSolver==CONJUGATE_HT_NATURAL) && (Mode.Buoyancy==LOW_MACH) )
  {
    fprintf(fp,"\t     Kind of Solver           :   Conjugate heat transfer / Natural convection : Low Mach Approximation\n");
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
  if ( KindOfSolver != SOLID_CONDUCTION )
  {
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
  fprintf(fp,"\t             Temperature      :   Celsius\n");
  
  
  // 時間制御 ------------------
  fprintf(fp,"\n\tTime Control\n");
  
  // 加速時間
  if ( Interval[tg_accelra].getMode() == IntervalManager::By_time )
  {
    itm = Interval[tg_accelra].getIntervalTime();
    fprintf(fp,"\t     Acceleration Time        :   %12.5e [sec] / %12.5e [-]\n", itm*Tscale, itm);
  }
  else
  {
    fprintf(fp,"\t     Acceleration Step        :   %12d\n", Interval[tg_accelra].getIntervalStep());
  }
  
  // 時間平均
  if ( Mode.Average == ON )
  {
    if ( Interval[tg_average].getMode() == IntervalManager::By_time )
    {
      itm = Interval[tg_average].getStartTime();
      stp = (unsigned)ceil(Interval[Control::tg_average].getStartTime() / dt);
      fprintf(fp,"\t     Averaging Start          :   %12.5e [sec] / %12.5e [-] : %12d [step]\n", itm*Tscale, itm, stp);
    }
    else
    {
      itm = (double)Interval[Control::tg_average].getStartStep() * dt;
      stp = Interval[tg_average].getStartStep();
      fprintf(fp,"\t     Averaging Start          :   %12.5e [sec] / %12.5e [-] : %12d [step]\n", itm*Tscale, itm, stp);
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
  
  // start & end
  if ( Interval[Control::tg_compute].getMode() == IntervalManager::By_step )
  {
    itm = (double)Interval[Control::tg_compute].getStartStep() * dt;
    stp = Interval[tg_compute].getStartStep();
  }
  else
  {
    itm = Interval[tg_compute].getStartTime();
    stp = (unsigned)ceil(Interval[Control::tg_compute].getStartTime() / dt);
  }
  fprintf(fp,"\t     Calculation Start        :   %12.5e [sec] / %12.5e [-] : %12d [step]\n", itm*Tscale, itm, stp);
  
  
  if ( Interval[Control::tg_compute].getMode() == IntervalManager::By_step )
  {
    itm = (double)Interval[Control::tg_compute].getLastStep() * dt;
    stp = Interval[tg_compute].getLastStep();
  }
  else
  {
    itm = Interval[tg_compute].getLastTime();
    stp = (unsigned)ceil(Interval[Control::tg_compute].getLastTime() / dt);
  }
  fprintf(fp,"\t     Calculation End          :   %12.5e [sec] / %12.5e [-] : %12d [step]\n", itm*Tscale, itm, stp);
  
  
  // Calculation time/step
  unsigned long m_stp;
  if ( Interval[Control::tg_compute].getMode() == IntervalManager::By_step )
  {
    itm = (double)Interval[Control::tg_compute].getIntervalStep() * dt;
    m_stp = (unsigned long)Interval[Control::tg_compute].getIntervalStep();
  }
  else
  {
    itm = Interval[tg_compute].getIntervalTime();
    m_stp = (unsigned long)(itm/dt);
  }

  fprintf(fp,"\t     Calculation Period       :   %12.5e [sec] / %12.5e [-] : %12ld [step]\n", itm*Tscale, itm, m_stp);
  
  
  
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
          (Mode.Log_Base == ON) ? "history_base.txt" : "", 
          (Mode.Log_Base == ON) ? "history_compo.txt" : "",
          (Mode.Log_Base == ON) ? "history_domainflux.txt" : "", 
          (Mode.Log_Base == ON) ? "history_force.txt" : "");
  fprintf(fp,"\t     Iteration Log            :   %4s  %s\n", 
          (Mode.Log_Itr == ON)?"ON >":"OFF ", (Mode.Log_Itr == ON) ? "history_iteration.txt" : "");
  fprintf(fp,"\t     Profiling report         :   %4s  %s%s\n", 
          (Mode.Profiling != OFF)?"ON >":"OFF ", 
          (Mode.Profiling == DETAIL)? "Detail mode, ":"",
          (Mode.Profiling != OFF)?"profiling.txt":"");
  
  fprintf(fp,"\t     Wall info. Log           :   %4s  %s\n", 
          (Mode.Log_Wall == ON)?"ON >":"OFF ", (Mode.Log_Wall == ON) ? "history_log_wall.txt" : "");
  
  
  // Intervals
  fprintf(fp,"\n\tIntervals\n");
  
  
  // インターバル ------------------
  // 基本履歴のコンソール出力 
  if ( Interval[tg_console].getMode() == IntervalManager::By_time )
  {
    itm = Interval[tg_console].getIntervalTime();
    fprintf(fp,"\t     Base Info.(stdout)       :   %12.6e [sec] / %12.6e [-]\n", itm*Tscale, itm);
  }
  else
  {
    fprintf(fp,"\t     Base Info.(stdout)       :   %12d [step]\n", Interval[tg_console].getIntervalStep());
  }
  
  // 履歴情報のファイル出力
  if ( Interval[tg_history].getMode() == IntervalManager::By_time )
  {
    itm = Interval[tg_history].getIntervalTime();
    fprintf(fp,"\t     Other Histories          :   %12.6e [sec] / %12.6e [-]\n", itm*Tscale, itm);
  }
  else
  {
    fprintf(fp,"\t     Other Histories          :   %12d [step]\n", Interval[tg_history].getIntervalStep());
  }
  
  // 基本変数のファイル出力
  if ( Interval[tg_basic].getMode() == IntervalManager::By_time )
  {
    itm = Interval[tg_basic].getIntervalTime();
    fprintf(fp,"\t     Basic Variables          :   %12.6e [sec] / %12.6e [-]\n", itm*Tscale, itm);
  }
  else
  {
    fprintf(fp,"\t     Basic Variables          :   %12d [step]\n", Interval[tg_basic].getIntervalStep());
  }
  
  // 平均値のファイル出力
  if ( Interval[tg_average].getMode() == IntervalManager::By_time )
  {
    itm = Interval[tg_average].getIntervalTime();
    fprintf(fp,"\t     Averaged Variables       :   %12.6e [sec] / %12.6e [-]\n", itm*Tscale, itm);
  }
  else
  {
    fprintf(fp,"\t     Averaged Variables       :   %12d [step]\n", Interval[tg_average].getIntervalStep());
  }
  
  // 派生変数のファイル出力
  if ( Interval[tg_derived].getMode() == IntervalManager::By_time )
  {
    itm = Interval[tg_derived].getIntervalTime();
    fprintf(fp,"\t     Derived Variables        :   %12.6e [sec] / %12.6e [-]\n", itm*Tscale, itm);
  }
  else
  {
    fprintf(fp,"\t     Derived Variables        :   %12d [step]\n", Interval[tg_derived].getIntervalStep());
  }
  
  // サンプリング情報のファイル出力
  if ( SamplingMode == ON )
  {
    if ( Interval[tg_sampled].getMode() == IntervalManager::By_time )
    {
      itm = Interval[tg_sampled].getIntervalTime();
      fprintf(fp,"\t     Sampled data             :   %12.6e [sec] / %12.6e [-]\n", itm*Tscale, itm);
    }
    else
    {
      fprintf(fp,"\t     Sampled data             :   %12d [step]\n", Interval[tg_sampled].getIntervalStep());
    }
  }
  

  // Criteria ------------------
  if ( em == ffvc_solver )
  {
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

  } // End of Criteria
  


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
  if ( varState[var_TotalP] == ON )
  {
    fprintf(fp,"\t     Total Pressure           :   ON\n");
  }
  else
  {
    fprintf(fp,"\t     Total Pressure           :   OFF\n");
  }
  
  //　渦度の出力モード
  if ( varState[var_Vorticity] == ON )
  {
    fprintf(fp,"\t     Vorticity                :   ON\n");
  }
  else
  {
    fprintf(fp,"\t     Vorticity                :   OFF\n");
  }
  
  //　ヘリシティの出力モード
  if ( varState[var_Helicity] == ON )
  {
    fprintf(fp,"\t     Helicity                 :   ON\n");
  }
  else
  {
    fprintf(fp,"\t     Helicity                 :   OFF\n");
  }
  
  //　速度勾配テンソルの第二不変量の出力モード
  if ( varState[var_Qcr] == ON ) {
    fprintf(fp,"\t     2nd Invariant of VGT     :   ON\n");
  }
  else
  {
    fprintf(fp,"\t     2nd Invariant of VGT     :   OFF\n");
  }
  
  
  // Driver ------------------
  fprintf(fp, "\n\tDriver parameter\n");
  if ( drv.length > 0.0 )
  {
    fprintf(fp,"\tDriver Length        [m] / [-]   : %12.5e / %12.5e\n", drv.length, drv.length/RefLength);
  }
  
  
  // Hidden parameter -----------------
  
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
    if ( OBC[n].getClass() == OBC_TRC_FREE ) c++;
  }
  if ( c>0 ) EnsCompo.tfree = ON;
  
  
  // コンポーネントモニター出力
  c = 0;
  for (int n=1; n<=NoCompo; n++)
  {
    if ( cmp[n].isCompoMonitor() ) c++;
  }
  if ( c>0 ) EnsCompo.monitor = ON;
  
}


// #################################################################
// コンポーネントと外部境界のパラメータを有次元に設定
void Control::setCmpParameters(MediumList* mat, CompoList* cmp, BoundaryOuter* BO)
{
  // コンポーネントの指定速度
  for (int n=1; n<=NoCompo; n++)
  {
    if ( cmp[n].getType()==SPEC_VEL )
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
        cmp[n].setHeatDensity( cmp[n].get_HeatValue() / ((REAL_TYPE)cmp[n].getElement()*vol) );
      }
      else // 発熱密度
      {
        cmp[n].setHeatValue( cmp[n].getHeatDensity() * ((REAL_TYPE)cmp[n].getElement()*vol) );
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
      ki = (mat[n].P[p_viscosity]*RefLength) / (mat[n].P[p_density]*RefVelocity);
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
    for (int n=0; n<NOFACE; n++)
    {
      switch ( BO[n].getClass() )
      {
        case OBC_OUTFLOW:
        case OBC_TRC_FREE:
        case OBC_FAR_FIELD:
          if ( BO[n].get_pType() == P_DIRICHLET )
          {
            BO[n].p = FBUtility::convPrsND2D(BO[n].p, BasePrs, RefDensity, RefVelocity, Unit.Prs); 
          }          
          break;
          
        case OBC_PERIODIC:
          if ( BO[n].getPrdcMode() != BoundaryOuter::prdc_Simple ) // Dirichlet or Bidirectionalを指定の場合
          {
            BO[n].p = FBUtility::convPrsND2D(BO[n].p, BasePrs, RefDensity, RefVelocity, Unit.Prs);
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
    iv.Pressure = FBUtility::convPrsND2D(iv.Pressure, BasePrs, RefDensity, RefVelocity, Unit.Prs);
    //iv.Density *= RefDensity;
		iv.VecU    *= RefVelocity;
		iv.VecV    *= RefVelocity;
		iv.VecW    *= RefVelocity;
    
    if ( isHeatProblem() )
    {
      for (int m=1; m<=NoMedium; m++)
      {
        cmp[m].setInitTemp( FBUtility::convTempND2D(cmp[m].getInitTemp(), BaseTemp, DiffTemp) );
      }
    }
	}
}


// #################################################################
/* 無次元パラメータを各種モードに応じて設定する
 * @note
 * - 代表長さと代表速度はパラメータで必ず与えること（読み込んだ値は変更しない）
 * - 純強制対流　有次元　（代表長さ，代表速度，動粘性係数，温度拡散係数）
 * -           無次元　（Pr, Re > RefV=RefL=1）
 * - 熱対流　　　有次元　（代表長さ，代表速度，温度差，体膨張係数，重力加速度，動粘性係数，温度拡散係数）
 * - 自然対流　　有次元　（代表長さ，代表速度，温度差，体膨張係数，重力加速度，動粘性係数，温度拡散係数）
 * - 固体熱伝導　有次元　（代表長さ，温度拡散係数)
 */
void Control::setRefParameters(MediumList* mat, ReferenceFrame* RF)
{
  REAL_TYPE rho, nyu, cp, lambda, beta, snd_spd=0.0;
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
  
  RefDensity      = rho;
  RefSpecificHeat = cp;
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
	else if ( (KindOfSolver==THERMAL_FLOW) || (KindOfSolver==CONJUGATE_HT) )
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
  else if ( (KindOfSolver==THERMAL_FLOW_NATURAL) || (KindOfSolver==CONJUGATE_HT_NATURAL) )
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

  
  if (Mode.PDE == PDE_EULER) Reynolds=1.0e23;
  
  Mach = RefVelocity / RefSoundSpeed;
  
  // タイミングパラメータの無次元化
  Tscale = (double)RefLength / (double)RefVelocity;
  
  
  // 参照速度の無次元化
  if ( Unit.Param == DIMENSIONAL )
  {
    double g[3];
    RF->copyGridVel(g);
    g[0] /= (double)RefVelocity;
    g[1] /= (double)RefVelocity;
    g[2] /= (double)RefVelocity;
    RF->setGridVel(g);
  }

}


// #################################################################
//単位ベクトルを計算して戻す
void Control::UnitVec(REAL_TYPE* v)
{
	REAL_TYPE a = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
  
	if ( a > 0.0 )
  {
		v[0] = v[0]/a;
		v[1] = v[1]/a;
		v[2] = v[2]/a;
	}
	else
  {
		v[0] = 0.0;
		v[1] = 0.0;
		v[2] = 0.0;
	}
}
