//##################################################################################
//
// Flow Base class
//
// Copyright (c) 2007-2011 VCAD System Research Program, RIKEN.
// All rights reserved.
//
// Copyright (c) 2011-2015 Institute of Industrial Science, The University of Tokyo.
// All rights reserved.
//
// Copyright (c) 2012-2015 Advanced Institute for Computational Science, RIKEN.
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
void DTcntl::set_Vars(const int m_kos, const int m_mode, const double m_min_dx, const double re, const double pe)
{
  KOS      = m_kos;
  mode     = m_mode;
  min_dx   = m_min_dx;
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
void ReferenceFrame::setGridVel(const REAL_TYPE* m_Gvel)
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
  
  REAL_TYPE u0 = v00[0];
  
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
void Control::copyCriteria(IterationCtl* IC, const string name)
{
  
  for (int i=0; i<NoBaseLS; i++)
  {
    if ( !strcasecmp( name.c_str(), Criteria[i].getAlias().c_str() ))
    {
      IC->copy(&Criteria[i]);
    }
  }
  
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
  
  
  // 乱流モデル
  getTurbulenceModel();
  
  // 初期擾乱
  getInitialPerturbation();
  
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
  

  getLog();
  
  
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
  
  
  // 安定化のフラグ (Hidden)
  label = "/ApplicationControl/StabilityControl/Control";
  
  if ( tpCntl->chkLabel(label) )
  {
    if ( tpCntl->getInspectedValue(label, str) )
    {
      if     ( !strcasecmp(str.c_str(), "On") )   Stab.control = ON;
      else if( !strcasecmp(str.c_str(), "Off") )  Stab.control = OFF;
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
  
  // 制御パラメータ >> パラメータは無次元で指定
  if (Stab.control == ON)
  {
    double ct=0.0;
    
    label = "/ApplicationControl/StabilityControl/Begin";
    
    if ( !(tpCntl->getInspectedValue(label, ct )) )
    {
      Hostonly_ stamped_printf("\tParsing error : fail to get '%s\n", label.c_str());
      Exit(0);
    }
    else
    {
      Stab.begin = ct;
    }
    
    
    label = "/ApplicationControl/StabilityControl/End";
    
    if ( !(tpCntl->getInspectedValue(label, ct )) )
    {
      Hostonly_ stamped_printf("\tParsing error : fail to get '%s\n", label.c_str());
      Exit(0);
    }
    else
    {
      Stab.end = ct;
    }
    
    // hard code
    Stab.penalty_number = 1.0e2;
  }
  
  
}




// #################################################################
// 計算内部領域の全セル数を返す
REAL_TYPE Control::getCellSize(const int* G_size)
{
  return (REAL_TYPE)G_size[0] * (REAL_TYPE)G_size[1] * (REAL_TYPE)G_size[2];
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
/**
 * @brief アプリケーションのDryRunパラメータを取得する
 */
void Control::getDryRun()
{
  string str;
  string label;
  
  // 境界条件確認のためのドライラン >> void IO_BASE::getFIOparams()でファイル出力をオフにする
  Hide.DryRun = OFF;
  
  label = "/ApplicationControl/DryrunBC";
  
  if ( tpCntl->chkLabel(label) )
  {
    if ( tpCntl->getInspectedValue(label, str) )
    {
      if     ( !strcasecmp(str.c_str(), "on") )  Hide.DryRun = ON;
      else if( !strcasecmp(str.c_str(), "off") ) Hide.DryRun = OFF;
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
// 計算モデルの入力ソース情報を取得
void Control::getGeometryModel()
{
  string str;
  string label;
  
  // ソース名を取得
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
  
  if     ( FBUtility::compare(str, "polygon") )           Mode.Example = id_Polygon;
  else if( FBUtility::compare(str, "ParallelPlate2D") )   Mode.Example = id_PPLT2D;
  else if( FBUtility::compare(str, "Duct") )              Mode.Example = id_Duct;
  else if( FBUtility::compare(str, "PerformanceTest") )   Mode.Example = id_PMT;
  else if( FBUtility::compare(str, "Rectangular") )       Mode.Example = id_Rect;
  else if( FBUtility::compare(str, "Cylinder") )          Mode.Example = id_Cylinder;
  else if( FBUtility::compare(str, "BackStep") )          Mode.Example = id_Step;
  else if( FBUtility::compare(str, "Sphere") )            Mode.Example = id_Sphere;
  else if( FBUtility::compare(str, "Jet") )               Mode.Example = id_Jet;
  else
  {
    Hostonly_ stamped_printf("\tError : Invalid source = '%s'\n", str.c_str());
    Exit(0);
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
  
}



// #################################################################
// 初期値擾乱
void Control::getInitialPerturbation()
{
  string str;
  string label, leaf;
  double ct=0.0;
  
  if ( Start == initial_start && LES.Calc == ON )
  {
    label="/StartCondition/InitialState/perturbation/DirectionOfChannelWall";
    
    if ( !(tpCntl->getInspectedValue(label, str )) )
    {
      Hostonly_ stamped_printf("\tParsing error : Invalid float value for '%s'\n", label.c_str());
      Exit(0);
    }
    
    if ( !strcasecmp("X", str.c_str()) )
    {
      LES.ChannelDir = X_minus;
    }
    else if ( !strcasecmp("Y", str.c_str()) )
    {
      LES.ChannelDir = Y_minus;
    }
    else if ( !strcasecmp("Z", str.c_str()) )
    {
      LES.ChannelDir = Z_minus;
    }
    
    
    
    label="/StartCondition/InitialState/perturbation/ChannelWidth";
    
    if ( !(tpCntl->getInspectedValue(label, ct )) )
    {
      Hostonly_ stamped_printf("\tParsing error : Invalid float value for '%s'\n", label.c_str());
      Exit(0);
    }
    LES.ChannelWidth = (REAL_TYPE)ct;
    
    
    
    label="/StartCondition/InitialState/perturbation/BulkVelocity";
    
    if ( !(tpCntl->getInspectedValue(label, ct )) )
    {
      Hostonly_ stamped_printf("\tParsing error : Invalid float value for '%s'\n", label.c_str());
      Exit(0);
    }
    LES.BulkVelocity = (REAL_TYPE)ct;
    
    
    
    label="/StartCondition/InitialState/perturbation/TubulenceReynoldsNumber";
    
    if ( !(tpCntl->getInspectedValue(label, ct )) )
    {
      Hostonly_ stamped_printf("\tParsing error : Invalid float value for '%s'\n", label.c_str());
      Exit(0);
    }
    LES.TurbulentReynoldsNum = (REAL_TYPE)ct;
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


    label = leaf + "/ResidualCriterion";
    if ( !(tpCntl->getInspectedValue(label, f_val )) )
    {
      Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
      Exit(0);
    }
    Criteria[i].setResCriterion(f_val);

    
    // 残差ノルム
    label = leaf + "/ResidualNorm";
    if ( !(tpCntl->getInspectedValue(label, str )) )
    {
      Hostonly_ printf("\tParsing error : No '%s'\n", label.c_str());
      Exit(0);
    }
    
    int ls = Criteria[i].getLS();
    
    // ノルム
    if ( !strcasecmp(str.c_str(), "RbyB") )
    {
      Criteria[i].setResType(nrm_r_b);
    }
    else if ( !strcasecmp(str.c_str(), "RbyR0") )
    {
      Criteria[i].setResType(nrm_r_r0);
    }
    else if ( !strcasecmp(str.c_str(), "RbyX") )
    {
      Criteria[i].setResType(nrm_r_x);
    }
    else
    {
      Hostonly_ stamped_printf("\tParsing error : Invalid keyword for '%s' of Norm for Poisson iteration\n", str.c_str());
      Exit(0);
    }

    
    // 誤差ノルム
    label = leaf + "/ErrorNorm";
    if ( !Criteria[i].setErrType(tpCntl, label) )
    {
      Hostonly_ printf("\tParsing error : '%s'\n", label.c_str());
      Exit(0);
    }
    
    
    // 固有パラメータ
    if ( !Criteria[i].getInherentPara(tpCntl, leaf) )
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
  
  
  // Log_Wall_Info
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
  label = "/BcTable/Boundary";
  
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
    REAL_TYPE ct;
    
    label="/Reference/Temperature/Base";
    
    if ( !(tpCntl->getInspectedValue(label, ct)) )
    {
      Hostonly_ stamped_printf("\tParsing error : Invalid float value for '%s'\n", label.c_str());
      Exit(0);
    }
    BaseTemp = ct;
    
    
    label="/Reference/Temperature/High";
    
    if ( !(tpCntl->getInspectedValue(label, ct)) )
    {
      Hostonly_ stamped_printf("\tParsing error : Invalid float value for '%s'\n", label.c_str());
      Exit(0);
    }
    HighTemp = ct;
    
    label="/Reference/Temperature/Low";
    
    if ( !(tpCntl->getInspectedValue(label, ct)) )
    {
      Hostonly_ stamped_printf("\tParsing error : Invalid float value for '%s'\n", label.c_str());
      Exit(0);
    }
    LowTemp = ct;
    
    
    if ( LowTemp > HighTemp )
    {
      Hostonly_ stamped_printf("\tTempereture difference must be High - Low > 0 \n");
      Exit(0);
    }
    DiffTemp = HighTemp - LowTemp;
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
    RF->setGridVel(xyz);
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
  
  if     ( !strcasecmp(str.c_str(), "Incompressible" ) )           BasicEqs = INCMP;
  else if( !strcasecmp(str.c_str(), "LimitedCompressibility" ) )   BasicEqs = LTDCMP;
  else if( !strcasecmp(str.c_str(), "Compressible" ) )             BasicEqs = CMPRSS;
  else if( !strcasecmp(str.c_str(), "Incompressible2Phase" ) )     BasicEqs = INCMP_2PHASE;
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
  else if( !strcasecmp(str.c_str(), "O2central") )   CnvScheme = O2_central;
  else if( !strcasecmp(str.c_str(), "O3muscl") )     CnvScheme = O3_muscl;
  else if( !strcasecmp(str.c_str(), "O4central") )   CnvScheme = O4_central;
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
    NvarsIns_plt3d = 3 + 3 + 1; // v, vf, p
    NvarsAvr_plt3d = 3 + 1;     // av, ap
  }
  else if( !strcasecmp(str.c_str(), "ThermalFlow" ) )
  {
    KindOfSolver = THERMAL_FLOW;
    varState[var_Velocity]    = true;
    varState[var_Fvelocity]   = true;
    varState[var_Pressure]    = true;
    varState[var_Temperature] = true;
    NvarsIns_plt3d = 3 + 3 + 1 + 1; // v, vf, p, t
    NvarsAvr_plt3d = 3 + 1 + 1;     // av, ap, at
  }
  else if( !strcasecmp(str.c_str(), "ThermalFlowNatural" ) )
  {
    KindOfSolver = THERMAL_FLOW_NATURAL;
    varState[var_Velocity]    = true;
    varState[var_Fvelocity]   = true;
    varState[var_Pressure]    = true;
    varState[var_Temperature] = true;
    NvarsIns_plt3d = 3 + 3 + 1 + 1; // v, vf, p, t
    NvarsAvr_plt3d = 3 + 1 + 1;     // av, ap, at
  }
  else if( !strcasecmp(str.c_str(), "ConjugateHeatTransfer" ) )
  {
    KindOfSolver = CONJUGATE_HT;
    varState[var_Velocity]    = true;
    varState[var_Fvelocity]   = true;
    varState[var_Pressure]    = true;
    varState[var_Temperature] = true;
    NvarsIns_plt3d = 3 + 3 + 1 + 1; // v, vf, p, t
    NvarsAvr_plt3d = 3 + 1 + 1;     // av, ap, at
  }
  else if( !strcasecmp(str.c_str(), "ConjugateHeatTransferNatural" ) )
  {
    KindOfSolver = CONJUGATE_HT_NATURAL;
    varState[var_Velocity]    = true;
    varState[var_Fvelocity]   = true;
    varState[var_Pressure]    = true;
    varState[var_Temperature] = true;
    NvarsIns_plt3d = 3 + 3 + 1 + 1; // v, vf, p, t
    NvarsAvr_plt3d = 3 + 1 + 1;     // av, ap, at
  }
  else if( !strcasecmp(str.c_str(), "SolidConduction" ) )
  {
    KindOfSolver = SOLID_CONDUCTION;
    varState[var_Temperature] = true;
    NvarsIns_plt3d = 1; // t
    NvarsAvr_plt3d = 1; // at
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
      case O1_upwind: // traction free条件はguide=2の必要がある
      case O2_central:
      case O3_muscl:
      case O4_central:
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
  
  // 丸め誤差の範囲でゼロの場合，イニシャルスタートとみなす
  if ( fabs(m_start)< ROUND_EPS )
  {
    Start = initial_start;
  }
  else
  {
    Start = restart_sameDiv_sameRes; // リスタートタイプのデフォルト
  }
  
  // Intereval Manager への登録 >> 他はFFV::initInterval()で指定
  Interval[tg_compute].setStart(m_start);
  Interval[tg_compute].setLast(m_end);
  Interval[tg_compute].setInterval(m_end-m_start); // tg_computeのインターバルは計算するセッションの長さ

  
  // By_time"のときのみRestartStepを指定する
  if ( Interval[tg_compute].getMode() == IntervalManager::By_time )
  {
    label="/TimeControl/Session/RestartStep";
    REAL_TYPE m_rstep;
    
    if ( !(tpCntl->getInspectedValue(label, m_rstep)) )
    {
      if ( Start != initial_start )
      {
        Hostonly_ stamped_printf("\tParsing error : fail to get '%s\n", label.c_str());
        Exit(0);
      }
    }
    else
    {
      Interval[tg_compute].restartStep = (unsigned)m_rstep;
    }
  }
  
  
  
  
  
  
  // 統計値の時刻指定モード
  label = "/TimeControl/Statistic/TemporalType";
  
  if ( !(tpCntl->getInspectedValue(label, str)) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    Exit(0);
  }
  else
  {
    if     ( !strcasecmp(str.c_str(), "step") )
    {
      Interval[tg_statistic].setMode(IntervalManager::By_step);
    }
    else if( !strcasecmp(str.c_str(), "time") )
    {
      Interval[tg_statistic].setMode(IntervalManager::By_time);
    }
    else
    {
      Hostonly_ stamped_printf("\tParsing error : Invalid keyword for '%s'\n", label.c_str());
      Exit(0);
    }
  }
  
  
  // 統計操作開始
  label = "/TimeControl/Statistic/Start";
  
  if ( !(tpCntl->getInspectedValue(label, ct )) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    Exit(0);
  }
  double stat_start = ct;
  
  
  // 統計操作終了
  label = "/TimeControl/Statistic/End";
  
  if ( !(tpCntl->getInspectedValue(label, ct )) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    Exit(0);
  }
  double stat_end = ct;
  
  
  // チェック
  if ( !(stat_start <= stat_end) )
  {
    Hostonly_ stamped_printf("\tError : Statistic/Start msut be less than Statistic/End.\n");
    Exit(0);
  }
  
  
  // Intereval Manager への登録 >> 他はFFV::initInterval()で指定
  Interval[tg_statistic].setStart(stat_start);
  Interval[tg_statistic].setLast(stat_end);
  Interval[tg_statistic].setInterval(stat_end-stat_start);
  
  
  // By_time"のときのみRestartStepを指定する
  if ( Interval[tg_statistic].getMode() == IntervalManager::By_time )
  {
    label="/TimeControl/Statistic/RestartStep";
    REAL_TYPE m_rstep;
    
    if ( !(tpCntl->getInspectedValue(label, m_rstep)) )
    {
      if ( Start != initial_start )
      {
        Hostonly_ stamped_printf("\tParsing error : fail to get '%s\n", label.c_str());
        Exit(0);
      }
    }
    else
    {
      Interval[tg_statistic].restartStep = (unsigned)m_rstep;
    }
  }
  

  /* 統計処理操作の判断
   
   1    |     2      |   3     << m_start (restart step)
   ---------+------------+-------
        ^            ^
     stat_start    stat_end
   
   case 1 : 統計操作は行うが，まだ指定時刻に到達していないので，統計値ファイルは存在せず，統計値のリスタートはない
        2 : 前セッションから継続して統計操作を行うが，既に統計値ファイルが存在する（はず）ので，統計値のリスタート処理を行う
        3 : 既に統計値操作の区間は終了しているので，統計操作は行わない
   */
  
  if ( Start == initial_start )
  {
    Mode.StatisticRestart = OFF; // default
    
    if ( (stat_end > 0.0) && (stat_start >= 0.0) ) // stat_end >= stat_startは既にチェック済み
    {
      Mode.Statistic = ON;
    }
    else
    {
      Mode.Statistic = OFF;
    }
  }
  else
  {
    if ( (stat_start < m_start) && (stat_end < m_start) ) // case 3
    {
      Mode.Statistic = OFF;
      Mode.StatisticRestart = OFF;
    }
    else if ( (stat_start < m_start) && (stat_end > m_start) ) // case 2
    {
      Mode.Statistic = ON;
      Mode.StatisticRestart = ON;
    }
    else if ( (stat_start >= m_start) && (stat_end > m_start) ) // case 1
    {
      Mode.Statistic = ON;
      Mode.StatisticRestart = OFF;
    }
  }
  
  if ( Mode.Statistic == ON )
  {
    varState[var_VelocityAvr] = true;
    varState[var_PressureAvr] = true;
    
    if ( isHeatProblem() ) varState[var_TemperatureAvr] = true;
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
    LES.Model = LES_no;
  }
  else if ( !strcasecmp(str.c_str(), "smagorinsky") )
  {
    LES.Calc = ON;
    LES.Model = LES_Smagorinsky;
  }
  else if ( !strcasecmp(str.c_str(), "csm") )
  {
    LES.Calc  = ON;
    LES.Model = LES_CSM;
  }
  else if ( !strcasecmp(str.c_str(), "wale") )
  {
    LES.Calc  = ON;
    LES.Model = LES_WALE;
  }
  else
  {
    Hostonly_ stamped_printf("\tParsing error : Invalid keyword for '%s'\n", label.c_str());
    Exit(0);
  }
  
  if ( LES.Calc == OFF ) return;
  
  
  // スマゴリンスキーモデル
  if ( LES.Model == LES_Smagorinsky )
  {
    // Cs係数
    label = "/TurbulenceModeling/Cs";
    if ( !(tpCntl->getInspectedValue(label, ct )) )
    {
      Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
      Exit(0);
    }
    LES.Cs = ct;
  }
  
  // 壁面
  label = "/TurbulenceModeling/VelocityProfile";
  
  if ( !(tpCntl->getInspectedValue(label, str)) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    Exit(0);
  }
  
  if      ( !strcasecmp(str.c_str(), "slip") )
  {
    LES.VelocityProfile = Slip;
  }
  else if ( !strcasecmp(str.c_str(), "noslip") )
  {
    LES.VelocityProfile = No_Slip;
  }
  else if ( !strcasecmp(str.c_str(), "lawofwall") )
  {
    LES.VelocityProfile = Law_of_Wall;
  }
  else if ( !strcasecmp(str.c_str(), "vandriest") )
  {
    LES.VelocityProfile = Van_Driest;
  }
  else
  {
    Hostonly_ stamped_printf("\tParsing error : Invalid keyword for '%s'\n", label.c_str());
    Exit(0);
  }
  
  LES.damping_factor = 25.0;
  
  
  // 初期擾乱
  label = "/TurbulenceModeling/InitialPerturbation";
  if ( !(tpCntl->getInspectedValue(label, str )) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    Exit(0);
  }
  if ( !strcasecmp(str.c_str(), "on") )
  {
    LES.InitialPerturbation = ON;
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
  REAL_TYPE df=0.0;
  
  if ( !(tpc->getInspectedValue(label, df)) )
  {
    stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    Exit(0);
  }
  
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
// ベクトル値(2D)を取得し，登録する
bool Control::getVec2(const std::string label, REAL_TYPE* v, TextParser* tpc)
{
  for (int i=0; i<2; i++) v[i]=0.0f;
  
  if( !(tpc->getInspectedVector(label, v, 2)) ) return false;
  
  return true;
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
  
  fprintf(fp,"\n----------\n\n");
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
  
  fprintf(fp,"\n----------\n\n");
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
  
  fprintf(fp,"\n----------\n\n");
  fprintf(fp,"\n\t>> Simulation Parameters\n\n");
  
  fprintf(fp,"\tReference Medium                      :  %s\n", mat[RefMat].alias.c_str());
  fprintf(fp,"\n");
  
  // Reference values
  fprintf(fp,"\tRef. Length               [m]         : %12.5e\n", RefLength);
  fprintf(fp,"\tRef. Velocity             [m/s]       : %12.5e\n", RefVelocity);
	fprintf(fp,"\tBase Pressure             [Pa]        : %12.5e\n", BasePrs);
  fprintf(fp,"\tRef. Mass Density         [kg/m^3]    : %12.5e\n", RefDensity);
  fprintf(fp,"\tRef. Specific Heat        [J/(kg K)]  : %12.5e\n", RefSpecificHeat);
  fprintf(fp,"\tRef. Kinematic Viscosity  [m^2/s]     : %12.5e\n", RefKviscosity);
  fprintf(fp,"\tRef. Speed of Sound       [m/s]       : %12.5e\n", RefSoundSpeed);
  fprintf(fp,"\tGravity                   [m/s^2]     : %12.5e\n", Gravity);
  fprintf(fp,"\tMach number               [-]         : %12.5e\n", Mach);
  fprintf(fp,"\n");
  
  fprintf(fp,"\tSpacing         X-dir.    [m] / [-]   : %12.5e / %12.5e\n", pitchD[0], pitch[0]);
  fprintf(fp,"\t                Y-dir.    [m] / [-]   : %12.5e / %12.5e\n", pitchD[1], pitch[1]);
  fprintf(fp,"\t                Z-dir.    [m] / [-]   : %12.5e / %12.5e\n", pitchD[2], pitch[2]);
  fprintf(fp,"\tTime Scale                [sec]       : %12.5e\n", Tscale);
  fprintf(fp,"\n");
  
  if ( isHeatProblem() )
  {
    fprintf(fp,"\tBase Temperature          [C] / [-]   : %12.5e / %3.1f\n", BaseTemp, 0.0);
    fprintf(fp,"\tHigh Temperature          [C] / [-]   : %12.5e / %3.1f\n", HighTemp, FBUtility::convTempD2ND(HighTemp, BaseTemp, DiffTemp));
    fprintf(fp,"\tLow  Temperature          [C] / [-]   : %12.5e / %3.1f\n", LowTemp,  FBUtility::convTempD2ND(LowTemp,  BaseTemp, DiffTemp));
    fprintf(fp,"\tTemperature Diff.         [C]         : %12.5e\n", DiffTemp);
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
void Control::printSteerConditions(FILE* fp,
                                   const DTcntl* DT,
                                   const ReferenceFrame* RF)
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
  
  fprintf(fp,"\n----------\n\n");
  fprintf(fp,"\n\t>> Library Information\n\n");
  fprintf(fp,"\t     CPMlib     Version %s\n", ver_CPM.c_str());
  fprintf(fp,"\t     CDMlib     Version %s\n", ver_CDM.c_str());
  fprintf(fp,"\t     Polylib    Version %s\n", ver_Poly.c_str());
  fprintf(fp,"\t     PMlib      Version %s\n", ver_PM.c_str());
  fprintf(fp,"\t     TextParser Version %s\n", ver_TP.c_str());
  fprintf(fp,"\n");
  
  fprintf(fp,"\n----------\n\n");
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
				fprintf(fp,"\t     Convective flux scheme   :   Upwind O(dx^1)\n");
        break;
        
      case O2_central:
				fprintf(fp,"\t     Convective flux scheme   :   Central O(dx^2)\n");
        break;
        
      case O4_central:
				fprintf(fp,"\t     Convective flux scheme   :   Central O(dx^4)\n");
        break;
        
			case O3_muscl:
				fprintf(fp,"\t     Convective flux scheme   :   MUSCL O(dx^3)\n");
        
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

  
  // 乱流モデル ------------------
  fprintf(fp, "\n\tTurbulence Modeling\n");
  switch (LES.Model)
  {
    case LES_no:
      fprintf(fp,"\t     Model                    :   Nothing (DNS)\n");
      break;
      
    case LES_Smagorinsky:
      fprintf(fp,"\t     Model                    :   Smagorinsky Cs=%f\n", LES.Cs);
      break;
      
    case LES_CSM:
      fprintf(fp,"\t     Model                    :   Coherent Smagorinsky Model (CSM) \n");
      break;
      
    case LES_WALE:
      fprintf(fp,"\t     Model                    :   Wall-Adapting Local Eddy-viscosity (WALE) \n");
      break;
      
    default:
      stamped_printf("Error: Turbulence Modeling section\n");
      err = false;
  }
  
  
  switch (LES.VelocityProfile)
  {
    case No_Slip:
      fprintf(fp,"\t     Velocity Profile         :   No Slip\n");
      break;
      
    case Slip:
      fprintf(fp,"\t     Velocity Profile         :   Slip\n");
      break;
      
    case Van_Driest:
      fprintf(fp,"\t     Velocity Profile         :   Van Driest function\n");
      break;
      
    case Law_of_Wall:
      fprintf(fp,"\t     Velocity Profile         :   Law of Wall\n");
      break;
      
    default:
      stamped_printf("Error: Turbulence Modeling section\n");
      err = false;
  }
  
  
  if ( LES.InitialPerturbation == ON )
  {
    fprintf(fp,"\t     Initial Perturbation     :   On \n");
    fprintf(fp,"\t     Direction of Channel wall:   %s \n", FBUtility::getDirection(LES.ChannelDir).c_str());
    fprintf(fp,"\t     Channel Width            :   %f \n", LES.ChannelWidth);
    fprintf(fp,"\t     Bulk Velocity            :   %e \n", LES.BulkVelocity);
    fprintf(fp,"\t     Turbulent Reynolds Number:   %e \n", LES.TurbulentReynoldsNum);
  }
  else
  {
    fprintf(fp,"\t     Initial Perturbation     :   Off \n");
  }

  
  // 統計 ------------------
  fprintf(fp, "\n\tStatistic Information\n");
  fprintf(fp,"\t     Pressure                 :   %s\n", ( Mode.StatPressure == ON ) ? "Yes" : "No");
  fprintf(fp,"\t     Velocity                 :   %s\n", ( Mode.StatVelocity == ON ) ? "Yes" : "No");
  fprintf(fp,"\t     Temperature              :   %s\n", ( Mode.StatTemperature == ON ) ? "Yes" : "No");
  fprintf(fp,"\t     ReynoldsStress           :   %s\n", ( Mode.ReynoldsStress == ON ) ? "Yes" : "No");
  

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
  if ( Mode.Statistic == ON )
  {
    if ( Interval[tg_statistic].getMode() == IntervalManager::By_time )
    {
      itm = Interval[tg_statistic].getStartTime();
      stp = (unsigned)ceil(Interval[Control::tg_statistic].getStartTime() / dt);
      fprintf(fp,"\t     Statistic Start          :   %12.5e [sec] / %12.5e [-] : %12d [step]\n", itm*Tscale, itm, stp);
    }
    else
    {
      itm = (double)Interval[Control::tg_statistic].getStartStep() * dt;
      stp = Interval[tg_statistic].getStartStep();
      fprintf(fp,"\t     Statistic Start          :   %12.5e [sec] / %12.5e [-] : %12d [step]\n", itm*Tscale, itm, stp);
    }
  }
  else
  {
    fprintf(fp,"\t     Statistic Start          :   OFF\n");
  }
  
  
  
  // Time Increment
  REAL_TYPE min_dx = std::min(pitch[0], std::min(pitch[1], pitch[2]));
  REAL_TYPE d_R = min_dx*min_dx*Reynolds/6.0; // 拡散数
  REAL_TYPE d_P = min_dx*min_dx*Peclet/6.0;   // 拡散数
  REAL_TYPE cfl = (REAL_TYPE)DT->get_CFL();
  switch ( DT->get_Scheme() ) 
  {
    case DTcntl::dt_direct:
      fprintf(fp,"\t     Time Increment dt        :   %12.5e [sec] / %12.5e [-] : Direct ", dt*Tscale, dt);
      if ( isHeatProblem() )
      {
        fprintf(fp,": Diff. Num. = %7.2e\n", dt/(min_dx*min_dx*Peclet));
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
  
  
  
  
  // ログ出力 ------------------
  fprintf(fp,"\n\tLogs\n");
  fprintf(fp,"\t     Unit for Output          :   %s\n", (Unit.Log == DIMENSIONAL) ? "Dimensional" : "Non-Dimensional");
  fprintf(fp,"\t     Base Logs                :   %4s  %s, %s, %s, %s\n", 
          (Mode.Log_Base == ON) ? "ON >" : "OFF ", 
          (Mode.Log_Base == ON) ? "history_base.txt" : "", 
          (Mode.Log_Base == ON) ? "history_compo.txt" : "",
          (Mode.Log_Base == ON) ? "history_domainflux.txt" : "", 
          (Mode.Log_Base == ON) ? "history_force.txt" : "");
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
    fprintf(fp,"\t     Basic/Derived Variables  :   %12.6e [sec] / %12.6e [-]\n", itm*Tscale, itm);
  }
  else
  {
    fprintf(fp,"\t     Basic/Derived Variables  :   %12d [step]\n", Interval[tg_basic].getIntervalStep());
  }
  
  // 統計値のファイル出力
  if ( Interval[tg_statistic].getMode() == IntervalManager::By_time )
  {
    itm = Interval[tg_statistic].getIntervalTime();
    fprintf(fp,"\t     Statistical Variables    :   %12.6e [sec] / %12.6e [-]\n", itm*Tscale, itm);
  }
  else
  {
    fprintf(fp,"\t     Statistical Variables    :   %12d [step]\n", Interval[tg_statistic].getIntervalStep());
  }
  
  /* 派生変数のファイル出力 >> 基本変数と同じにする 20141025
  if ( Interval[tg_derived].getMode() == IntervalManager::By_time )
  {
    itm = Interval[tg_derived].getIntervalTime();
    fprintf(fp,"\t     Derived Variables        :   %12.6e [sec] / %12.6e [-]\n", itm*Tscale, itm);
  }
  else
  {
    fprintf(fp,"\t     Derived Variables        :   %12d [step]\n", Interval[tg_derived].getIntervalStep());
  }*/
  
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
  
  
  // 派生変数 ------------------
  fprintf(fp, "\n\tDerived variables\n");
  
  //　全圧の出力モード
  fprintf(fp,"\t     Total Pressure           :   ");
  fprintf(fp, ( varState[var_TotalP] == ON ) ? "ON\n" : "OFF\n");
  
  //　渦度の出力モード
  fprintf(fp,"\t     Vorticity                :   ");
  fprintf(fp, ( varState[var_Vorticity] == ON ) ? "ON\n" : "OFF\n");
  
  //　ヘリシティの出力モード
  fprintf(fp,"\t     Helicity                 :   ");
  fprintf(fp, ( varState[var_Helicity] == ON ) ? "ON\n" : "OFF\n");
  
  //　速度勾配テンソルの第二不変量の出力モード
  fprintf(fp,"\t     2nd Invariant of VGT     :   ");
  fprintf(fp, ( varState[var_Qcr] == ON ) ? "ON\n" : "OFF\n");
  
  //　発散値
  fprintf(fp,"\t     Divergence               :   ");
  fprintf(fp, ( varState[var_Div] == ON ) ? "ON\n" : "OFF\n");

  //　RMS
  fprintf(fp,"\t     RMS of velocity          :   ");
  fprintf(fp, ( varState[var_RmsV] == ON ) ? "ON\n" : "OFF\n");
  
  //　RMSmean
  fprintf(fp,"\t     RMS Mean of velocity     :   ");
  fprintf(fp, ( varState[var_RmsMeanV] == ON ) ? "ON\n" : "OFF\n");
  
  
  // Driver ------------------
  fprintf(fp, "\n\tDriver parameter\n");
  if ( drv.length > 0.0 )
  {
    fprintf(fp,"\tDriver Length        [m] / [-]   : %12.5e / %12.5e\n", drv.length, drv.length/RefLength);
  }
  
  
  // Stability Control
  if (Stab.control == ON)
  {
    fprintf(fp,"\n\tStability Control\n");
    fprintf(fp,"\t\tBegin velocity     [-]   : %12.5e\n", Stab.begin);
    fprintf(fp,"\t\tEnd   velocity     [-]   : %12.5e\n", Stab.end);
    fprintf(fp,"\t\tPenalty value      [-]   : %12.5e\n", Stab.penalty_number);
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
int Control::setExistComponent(CompoList* cmp, BoundaryOuter* OBC, int* g_obstacle)
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
  
  
  // 物体
  c = 0;
  for (int n=1; n<=NoCompo; n++)
  {
    if ( cmp[n].getType() == OBSTACLE || cmp[n].getType() == SOLIDREV )
    {
      g_obstacle[n] = ON;
      c++;
    }
  }
  if ( c>0 ) EnsCompo.obstacle = ON;
  
  return c; // OBSTACLEの個数
}


// #################################################################
// コンポーネントと外部境界のパラメータを有次元に設定
void Control::setCmpParameters(MediumList* mat, CompoList* cmp, BoundaryOuter* BO)
{
  /* コンポーネントの指定速度
  for (int n=1; n<=NoCompo; n++)
  {
    if ( cmp[n].getType()==SPEC_VEL )
    {
			if ( cmp[n].isPolicy_Massflow() ) //ポリシーが流量の場合
      {
				if ( Unit.Param == DIMENSIONAL )
        {
					cmp[n].set_Velocity( cmp[n].get_Massflow() / cmp[n].area );
				}
				else
        {
					cmp[n].set_Velocity( cmp[n].get_Massflow()*RefVelocity*RefLength*RefLength / cmp[n].area );
				}
			}
      
			// 流量指定のときのみ，ca[]に有次元速度パラメータを保存  >> 速度指定の場合には，parseBC::getIbcSpecVel()で設定済み
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
   */
  
	
  // 発熱密度の計算(有次元) -- 発熱量と発熱密度
  REAL_TYPE vol = pitchD[0]*pitchD[1]*pitchD[2];
  
  for (int n=1; n<=NoCompo; n++)
  {
    if ( cmp[n].getType()==HEAT_SRC ) cmp[n].setHeatSrcValue(vol);
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
          if ( BO[n].get_V_Profile() == CompoList::vel_harmonic )
          {
            BO[n].ca[CompoList::amplitude] *= RefVelocity;
            BO[n].ca[CompoList::frequency] *= (RefVelocity/RefLength);
            //BO[n].ca[CompoList::initphase] radは有次元化不要
            BO[n].ca[CompoList::bias]      *= RefVelocity;
          }
          else if ( BO[n].get_V_Profile() == CompoList::vel_polynomial6 )
          {
            for (int i=0; i<6; i++) BO[n].ca[i] *= RefVelocity;
          }
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
  RefKviscosity   = nyu;
  
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
    REAL_TYPE g[3];
    RF->copyGridVel(g);
    g[0] /= (double)RefVelocity;
    g[1] /= (double)RefVelocity;
    g[2] /= (double)RefVelocity;
    RF->setGridVel(g);
  }
  
  
  // Mach
  Mach = RefVelocity / 340.0;
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
