// #################################################################
//
// CAERU Library
//
// Copyright (c) 2012 All right reserved.
//
// Institute of Industrial Science, The University of Tokyo, Japan. 
//
// #################################################################

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
  else if ( !strcasecmp(str, "CFL_Reference_Velocity") )
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
  else if ( !strcasecmp(str, "CFL_Diffusion_Reference_Velocity") )
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
// 熱交換器パラメータの変換（水と水銀）
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
// 熱交換器パラメータの変換（Pa）
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
// labelのコンポーネント数を返す
int Control::countCompo(CompoList* cmp, const int label)
{
  int cnt=0;
  
  for (int i=1; i<=NoBC; i++) {
    if ( cmp[i].getType() == label ) cnt++;
  }
  return cnt;
}


// #################################################################
// 制御，計算パラメータ群の表示
void Control::displayParams(FILE* mp, FILE* fp, ItrCtl* IC, DTcntl* DT, ReferenceFrame* RF, MediumList* mat, FileIO_PLOT3D_WRITE* FP3DW)
{
  printSteerConditions(mp, IC, DT, RF, FP3DW);
  printSteerConditions(fp, IC, DT, RF, FP3DW);
  printParaConditions(mp, mat);
  printParaConditions(fp, mat);
  printInitValues(mp);
  printInitValues(fp);
}


// #################################################################
// MediumList中に登録されているkeyに対するIDを返す。発見できない場合はzero 
int Control::find_ID_from_Label(MediumList* mat, const int Nmax, const std::string key)
{
  std::string str = key;
  
  for (int i=1; i<=Nmax; i++) 
  {
    if ( !strcasecmp(str.c_str(), mat[i].getLabel().c_str()) ) return i;
  }
  
  return 0;
}


// #################################################################
// 反復の収束判定パラメータを取得
void Control::findCriteria(const string label0, const int order, ItrCtl* IC)
{
  int itr=0;
  REAL_TYPE tmp=0.0;
  string str, label;
  
  if ( tpCntl->chkNode(label0) )
  {
    label = label0 + "/Iteration";
    if ( !(tpCntl->GetValue(label, &itr )) ) 
    {
      Hostonly_ stamped_printf("\tParsing error : Invalid integer value for '%s'\n", label.c_str());
      Exit(0);
    }
    IC->set_ItrMax(itr);
    
    if (order == ItrCtl::ic_div)
    {
      if (itr == 1)
      {
        IC->set_ItrMax(2); // minimum 2
      }
    }
    
    
    label = label0 + "/Epsilon";
    if ( !(tpCntl->GetValue(label, &tmp )) ) 
    {
      Hostonly_ stamped_printf("\tParsing error : Invalid float value for '%s'\n", label.c_str());
      Exit(0);
    }
    IC->set_eps((double)tmp);
    
    
    label = label0 + "/norm";
    if ( !(tpCntl->GetValue(label, &str )) )
    {
      Hostonly_ stamped_printf("\tParsing error : Invalid char* value for '%s'\n", label.c_str());
      Exit(0);
    }
    
    
    
    // normのタイプ
    switch (order)
    {
      case ItrCtl::ic_prs_pr: // Predictor phase
      case ItrCtl::ic_prs_cr: // Corrector phase
      case ItrCtl::ic_vis_cn: // Velocity Crank-Nicolosn
      case ItrCtl::ic_tdf_ei: // Temperature Euler Implicit
        
        if ( !strcasecmp(str.c_str(), "dx_b") )
        {
          IC->set_normType(ItrCtl::dx_b);
        }
        else if ( !strcasecmp(str.c_str(), "r_b") )
        {
          IC->set_normType(ItrCtl::r_b);
        }
        else if ( !strcasecmp(str.c_str(), "r_r0") )
        {
          IC->set_normType(ItrCtl::r_r0);
        }
        else
        {
          Hostonly_ stamped_printf("\tParsing error : Invalid keyword for '%s' of Norm for Poisson iteration\n", str.c_str());
          Exit(0);
        }
        break;
        
      case ItrCtl::ic_div: // VP iteration
        if ( !strcasecmp(str.c_str(), "v_div_max") )
        {
          IC->set_normType(ItrCtl::v_div_max);
        }
        else if ( !strcasecmp(str.c_str(), "v_div_dbg") )
        {
          IC->set_normType(ItrCtl::v_div_dbg);
          Mode.Log_Itr == ON;
        }
        else
        {
          Hostonly_ stamped_printf("\tParsing error : Invalid keyword for '%s' for heat iteration\n", str.c_str());
          Exit(0);
        }
        break;
        
        default:
          Hostonly_ stamped_printf("\tParsing error : Invalid keyword for '%s' for heat iteration\n", str.c_str());
          Exit(0);
    }
    
    
    if ( order != ItrCtl::ic_div )
    {
      label = label0 + "/Omega";
      if ( !(tpCntl->GetValue(label, &tmp )) )
      {
        Hostonly_ stamped_printf("\tParsing error : Invalid float value for '%s'\n", label.c_str());
        Exit(0);
      }
      IC->set_omg(tmp);
      
      
      label = label0 + "/comm_mode";
      if ( !(tpCntl->GetValue(label, &str )) )
      {
        Hostonly_ stamped_printf("\tParsing error : Invalid char* value for '%s'\n", label.c_str());
        Exit(0);
      }
      if ( !strcasecmp(str.c_str(), "sync") )
      {
        IC->set_SyncMode(comm_sync);
      }
      else if ( !strcasecmp(str.c_str(), "async") )
      {
        IC->set_SyncMode(comm_async);
      }
      else
      {
        Hostonly_ stamped_printf("\tParsing error : Invalid char* value for '%s'\n", label.c_str());
        Exit(0);
      }
      

      label = label0 + "/Linear_Solver";
      if ( !(tpCntl->GetValue(label, &str )) )
      {
        Hostonly_ stamped_printf("\tParsing error : Invalid char* value for '%s'\n", label.c_str());
        Exit(0);
      }
      
      // 線形ソルバーの種類
      if     ( !strcasecmp(str.c_str(), "SOR") )         IC->set_LS(SOR);
      else if( !strcasecmp(str.c_str(), "SOR2SMA") )     IC->set_LS(SOR2SMA);
      else if( !strcasecmp(str.c_str(), "SOR2CMA") )     IC->set_LS(SOR2CMA);
      else if( !strcasecmp(str.c_str(), "JACOBI") )      IC->set_LS(JACOBI);
      else if( !strcasecmp(str.c_str(), "GMRES") )       IC->set_LS(GMRES);
      else
      {
        Hostonly_ stamped_printf("\tInvalid keyword is described for Linear_Solver\n");
        Exit(0);
      }
    }
    
  }
  else 
  {
    Hostonly_ stamped_printf("\tParsing error : Invalid keyword of '%s'\n", label.c_str());
    Exit(0);
  }

}



// #################################################################
// 解法アルゴリズムを選択する
void Control::get_Algorithm()
{
  string str;
  string label;
  
  label = "/Steer/Algorithm/Flow";
  
  if ( !(tpCntl->GetValue(label, &str )) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '/Steer/Algorithm/Flow'\n");
    Exit(0);
  }
  
  if     ( !strcasecmp(str.c_str(), "FS_C_EE_D_EE") )     AlgorithmF = Flow_FS_EE_EE;
  else if( !strcasecmp(str.c_str(), "FS_C_RK_D_CN") )     AlgorithmF = Flow_FS_RK_CN;
  else if( !strcasecmp(str.c_str(), "FS_C_AB_D_AB") )     AlgorithmF = Flow_FS_AB2;
  else if( !strcasecmp(str.c_str(), "FS_C_AB_D_CN") )     AlgorithmF = Flow_FS_AB_CN;
  else
  {
    Hostonly_ stamped_printf("\tParsing error : Invalid keyword for '/Steer/Algorithm/Flow'\n");
    Exit(0);
  }
  
  // Heat
  if ( isHeatProblem() )
  {
	  label = "/Steer/Algorithm/Heat";
    
	  if ( !(tpCntl->GetValue(label, &str )) )
    {
		  Hostonly_ stamped_printf("\tParsing error : fail to get '/Steer/Algorithm/Heat'\n");
		  Exit(0);
	  }
    
	  if     ( !strcasecmp(str.c_str(), "C_EE_D_EE") )    AlgorithmH = Heat_EE_EE;
	  else if( !strcasecmp(str.c_str(), "C_EE_D_EI") )    AlgorithmH = Heat_EE_EI;
	  else
    {
		  Hostonly_ stamped_printf("\tParsing error : Invalid keyword for '/Steer/Algorithm/Heat'\n");
		  Exit(0);
	  }
  }
}


// #################################################################
// 平均値操作に関するパラメータを取得する
// パラメータは，setParameters()で無次元して保持
void Control::get_Average_option()
{
  REAL_TYPE ct;
  string str;
  string label;
  
  label="/Steer/Average_option/Operation";
  
  if ( !(tpCntl->GetValue(label, &str )) ) {
    Hostonly_ stamped_printf("\tParsing error : fail to get '/Steer/Average_option/Operation'\n");
	  Exit(0);
  }
  
  if     ( !strcasecmp(str.c_str(), "on") )   Mode.Average = ON;
  else if( !strcasecmp(str.c_str(), "off") )  Mode.Average = OFF;
  else
  {
    Hostonly_ stamped_printf("\tInvalid keyword is described for '/Steer/Average_option/Operation'\n");
    Exit(0);
  }
  
  // 平均操作開始時間
  if ( Mode.Average == ON )
  {
	  label = "/Steer/Average_option/Start_Type";
    
	  if ( !(tpCntl->GetValue(label, &str )) )
    {
		  Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
		  Exit(0);
	  }
	  else
    {
		  if     ( !strcasecmp(str.c_str(), "step") )
      {
			  Interval[Interval_Manager::tg_avstart].setMode_Step();
		  }
		  else if( !strcasecmp(str.c_str(), "time") )
      {
			  Interval[Interval_Manager::tg_avstart].setMode_Time();
		  }
		  else
      {
			  Hostonly_ stamped_printf("\tParsing error : Invalid keyword for '%s'\n", label.c_str());
			  Exit(0);
		  }
		  
		  label = "/Steer/Average_option/Start";
		  if ( !(tpCntl->GetValue(label, &ct )) )
      {
			  Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
			  Exit(0);
		  }
		  else
      {
			  Interval[Interval_Manager::tg_avstart].setInterval((double)ct);
		  }
	  }
  }
}



// #################################################################
// Cell IDのゼロを指定IDに変更するオプションを取得する（隠しパラメータ）
// 'Change_ID'の文字列チェックはしないので注意して使うこと
void Control::get_ChangeID()
{
  int ct=0;
  string label;
  
  label="/Steer/Change_ID";
  if ( !(tpCntl->GetValue(label, &ct )) )
  {
	  return;
  }
  
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



// #################################################################
// パラメータ入力チェックモードの取得
void Control::get_CheckParameter()
{
  string str;
  string label;
  
  label = "/Steer/Check_Parameter";
  
  if ( !(tpCntl->GetValue(label, &str )) )
  {
	  Exit(0);
  }
  
  if     ( !strcasecmp(str.c_str(), "On") )   CheckParam = ON;
  else if( !strcasecmp(str.c_str(), "Off") )  CheckParam = OFF;
  else
  {
    Hostonly_ printf("\tInvalid keyword is described for '/Steer/Check_Parameter'\n");
    Exit(0);
  }
}


// #################################################################
// 計算内部領域の全セル数を返す
REAL_TYPE Control::getCellSize(const int* G_size)
{
  return (REAL_TYPE)G_size[0] * (REAL_TYPE)G_size[1] * (REAL_TYPE)G_size[2];
}



// #################################################################
// 対流項スキームのパラメータを取得する
void Control::get_Convection()
{
  
  int ct;
  string str;
  string label;
  
  // scheme
  label="/Steer/Convection_Term/scheme";
  
  if ( !(tpCntl->GetValue(label, &str )) ) {
	  Hostonly_ stamped_printf("\tParsing error : Invalid char* value for '/Steer/Convection_Term/scheme'\n");
	  Exit(0);
  }
  
  if     ( !strcasecmp(str.c_str(), "O1_Upwind") )    CnvScheme = O1_upwind;
  else if( !strcasecmp(str.c_str(), "O3_muscl") )     CnvScheme = O3_muscl;
  else if( !strcasecmp(str.c_str(), "O2_central") )   CnvScheme = O2_central;
  else if( !strcasecmp(str.c_str(), "O4_central") ) { CnvScheme = O4_central; Exit(0); }  // not yet implemented
  else
  {
    Hostonly_ stamped_printf("\tInvalid keyword is described for '/Steer/Convection_Term/scheme'\n");
    Exit(0);
  }
  
  // Limiter
  if ( CnvScheme == O3_muscl )
  {
		label="/Steer/Convection_Term/limiter";
    
		if ( !(tpCntl->GetValue(label, &str )) )
    {
			Hostonly_ stamped_printf("\tParsing error : Invalid char* value for '/Steer/Convection_Term/limiter'\n");
			Exit(0);
		}
    
		if     ( !strcasecmp(str.c_str(), "No_Limiter") )  Limiter = No_Limiter;
		else if( !strcasecmp(str.c_str(), "Minmod") )      Limiter = MinMod;
		else
    {
			Hostonly_ stamped_printf("\tInvalid keyword is described for '/Steer/Convection_Term/limiter'\n");
			Exit(0);
		}
  }
}


// #################################################################
// 派生して計算する変数のオプションを取得する
void Control::get_Derived()
{
  string str;
  string label;
  
  // 全圧
  label="/Steer/Derived_Variable/Total_Pressure";
  
  if ( !(tpCntl->GetValue(label, &str )) )
  {
	  Hostonly_ stamped_printf("\tParsing error : fail to get '/Steer/Derived_Variable/Total_Pressure'\n");
	  Exit(0);
  }
  
  if     ( !strcasecmp(str.c_str(), "on") )  Mode.TP = ON;
  else if( !strcasecmp(str.c_str(), "off") ) Mode.TP = OFF;
  else
  {
    Hostonly_ stamped_printf("\tInvalid keyword is described for '/Steer/Derived_Variable/Total_Pressure\n");
    Exit(0);
  }
  
  // 渦度ベクトル
  label="/Steer/Derived_Variable/Vorticity";
  
  if ( !(tpCntl->GetValue(label, &str )) )
  {
	  Hostonly_ stamped_printf("\tParsing error : fail to get '/Steer/Derived_Variable/Vorticity'\n");
	  Exit(0);
  }
  
  if     ( !strcasecmp(str.c_str(), "on") )  Mode.VRT = ON;
  else if( !strcasecmp(str.c_str(), "off") ) Mode.VRT = OFF;
  else
  {
    Hostonly_ stamped_printf("\tInvalid keyword is described for '/Steer/Derived_Variable/Vorticity'\n");
    Exit(0);
  }
  
  // 速度勾配テンソルの第2不変量
  label="/Steer/Derived_Variable/2nd_Invariant_of_VGT";
  
  if ( !(tpCntl->GetValue(label, &str )) )
  {
	  Hostonly_ stamped_printf("\tParsing error : fail to get '/Steer/Derived_Variable/2nd_Invariant_of_VGT'\n");
	  Exit(0);
  }
  
  if     ( !strcasecmp(str.c_str(), "on") )  Mode.I2VGT = ON;
  else if( !strcasecmp(str.c_str(), "off") ) Mode.I2VGT = OFF;
  else
  {
    Hostonly_ stamped_printf("\tInvalid keyword is described for '/Steer/Derived_Variable/2nd_Invariant_of_VGT'\n");
    Exit(0);
  }
  
  // ヘリシティ
  label="/Steer/Derived_Variable/Helicity";
  
  if ( !(tpCntl->GetValue(label, &str )) )
  {
	  Hostonly_ stamped_printf("\tParsing error : fail to get '/Steer/Derived_Variable/Helicity'\n");
	  Exit(0);
  }
  
  if     ( !strcasecmp(str.c_str(), "on") )  Mode.Helicity = ON;
  else if( !strcasecmp(str.c_str(), "off") ) Mode.Helicity = OFF;
  else
  {
    Hostonly_ stamped_printf("\tInvalid keyword is described for '/Steer/Derived_Variable/Helicity\n");
    Exit(0);
  }
  
}


// #################################################################
// ファイル入出力に関するパラメータを取得し，sphフォーマットの出力の並列モードを指定する．
// インターバルパラメータは，setParameters()で無次元して保持
void Control::get_FileIO()
{
  
  REAL_TYPE f_val=0.0;
  REAL_TYPE ct;
  string str;
  string label, leaf;
  
  // 出力単位
  label = "/Steer/File_IO/Unit_of_file";
  
  if ( !(tpCntl->GetValue(label, &str )) )
  {
	  Hostonly_ stamped_printf("\tParsing error : Invalid string for '%s'\n", label.c_str());
	  Exit(0);
  }
  
  if     ( !strcasecmp(str.c_str(), "Dimensional") )      Unit.File = DIMENSIONAL;
  else if( !strcasecmp(str.c_str(), "Non_Dimensional") )  Unit.File = NONDIMENSIONAL;
  else
  {
    Hostonly_ stamped_printf("\tInvalid keyword is described for '%s'\n", label.c_str());
    Exit(0);
  }
  
  
  // 出力ガイドセルモード
  label = "/Steer/File_IO/Guide_Out";
  
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
  
  
  // デバッグ用のdiv(u)の出力指定
  label = "/Steer/File_IO/Debug_divergence";
  
  if ( !(tpCntl->GetValue(label, &str )) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
	  Exit(0);
  }
  
  if     ( !strcasecmp(str.c_str(), "on") )    FIO.Div_Debug = ON;
  else if( !strcasecmp(str.c_str(), "off") )   FIO.Div_Debug = OFF;
  else
  {
    Hostonly_ stamped_printf("\tInvalid keyword is described for '%s'\n", label.c_str());
    Exit(0);
  }
  
  
  // ボクセルファイル出力
  label = "/Steer/File_IO/Voxel_Output";
  
  if ( !(tpCntl->GetValue(label, &str )) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
	  Exit(0);
  }
  
  if     ( !strcasecmp(str.c_str(), "svx") )  FIO.IO_Voxel = Sphere_SVX;
  else if( !strcasecmp(str.c_str(), "off") )  FIO.IO_Voxel = OFF;
  else
  {
    Hostonly_ stamped_printf("\tInvalid keyword is described for '%s'\n", label.c_str());
    Exit(0);
  }
  
  
  // インターバル 瞬時値
  label = "/Steer/File_IO/Instant_Interval_Type";
  
  if ( !(tpCntl->GetValue(label, &str )) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    Exit(0);
  }
  else
  {
    if     ( !strcasecmp(str.c_str(), "step") )
    {
      Interval[Interval_Manager::tg_instant].setMode_Step();
    }
    else if( !strcasecmp(str.c_str(), "time") )
    {
      Interval[Interval_Manager::tg_instant].setMode_Time();
    }
    else
    {
      Hostonly_ stamped_printf("\tParsing error : Invalid keyword for '%s'\n", label.c_str());
      Exit(0);
    }
    
    label="/Steer/File_IO/Instant_Interval";
    
    if ( !(tpCntl->GetValue(label, &f_val )) )
    {
      Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
      Exit(0);
    }
    else
    {
      Interval[Interval_Manager::tg_instant].setInterval((double)f_val);
    }
  }
  
  
  // インターバル　平均値
  label = "/Steer/File_IO/Averaged_Interval_Type";
  
  if ( !(tpCntl->GetValue(label, &str )) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    Exit(0);
  }
  else
  {
    if ( !strcasecmp(str.c_str(), "step") )
    {
      Interval[Interval_Manager::tg_average].setMode_Step();
    }
    else if ( !strcasecmp(str.c_str(), "time") )
    {
      Interval[Interval_Manager::tg_average].setMode_Time();
    }
    else
    {
      Hostonly_ stamped_printf("\tParsing error : Invalid keyword for '%s'\n", label.c_str());
      Exit(0);
    }
    
    label="/Steer/File_IO/Averaged_Interval";
    
    if ( !(tpCntl->GetValue(label, &f_val )) )
    {
      Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
      Exit(0);
    }
    else
    {
      Interval[Interval_Manager::tg_average].setInterval((double)f_val);
    }
  }
  
  
  // ファイル出力モード
  label = "/Steer/File_IO/Output_Mode";
  
  if( !tpCntl->chkNode(label) )
  {
    Hostonly_ stamped_printf("\tParsing error : Missing the section of '%s'\n", label.c_str());
    Exit(0);
  }
  
  leaf = label + "/Mode";
  
  if ( !(tpCntl->GetValue(leaf, &str )) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", leaf.c_str());
    Exit(0);
  }
  else
  {
    if ( !strcasecmp(str.c_str(), "current") )
    {
      FIO.IO_Mode = io_current;
    }
    else if ( !strcasecmp(str.c_str(), "specified") )
    {
      FIO.IO_Mode = io_specified;
    }
    else if ( !strcasecmp(str.c_str(), "time_slice") )
    {
      FIO.IO_Mode = io_time_slice;
    }
    else
    {
      Hostonly_ stamped_printf("\tParsing error : Invalid keyword for '%s'\n", leaf.c_str());
      Exit(0);
    }
    
    // Output Directory_Path
    if ( FIO.IO_Mode == io_specified )
    {
      label = "/Steer/File_IO/Output_Mode/Directory_Path";
      
      if ( !(tpCntl->GetValue(label, &str)) )
      {
        Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
        Exit(0);
      }
      
      // 指定が無ければ，空のまま
      if ( !str.empty() )
      {
        FIO.IO_DirPath = str;
      }
    }

  }
  
  
  // PLOT3D Option
  label = "/Steer/File_IO/PLOT3D";
  
  if( !tpCntl->chkNode(label) )
  {
    Hostonly_ stamped_printf("\tParsing error : Missing the section of '%s'\n", label.c_str());
    Exit(0);
  }
  
  leaf = label + "/Option";
  
  if ( !tpCntl->GetValue(leaf, &str) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", leaf.c_str());
    Exit(0);
  }
  if     ( !strcasecmp(str.c_str(), "off") ) FIO.PLOT3D_OUT = OFF;
  else if( !strcasecmp(str.c_str(), "on") )  FIO.PLOT3D_OUT = ON;
  else
  {
    Hostonly_ stamped_printf("\tInvalid keyword is described for '%s'\n", leaf.c_str());
    Exit(0);
  }
  
  // インターバル PLOT3D
  if(FIO.PLOT3D_OUT == ON)
  {
    label = "/Steer/File_IO/PLOT3D/Interval_Type";
    
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
      
      label="/Steer/File_IO/PLOT3D/Interval";
      
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
  
}


// #################################################################
// モデル形状情報
void Control::get_Geometry(const MediumTableInfo *MTITP)
{
  string str;
  string label;
  
  // ファイル名を取得
  label = "/Steer/Geometry_Model/Polylib_File";
  
  if ( !(tpCntl->GetValue(label, &str )) )
  {
    Hostonly_ stamped_printf("\tParsing error : Invalid char* value in '%s'\n", label.c_str());
    Exit(0);
  }
  PolylibConfigName = str;
  
  
  // フィルの流体媒質番号の指定
  label = "/Steer/Geometry_Model/Fluid_Medium_for_Fill";
  
  if ( !(tpCntl->GetValue(label, &str )) )
  {
    Hostonly_ stamped_printf("\tParsing error : Invalid value for '%s'\n", label.c_str());
    Exit(0);
  }
  
  // ラベル名が媒質リストにあるか否かを確認
  for (int i=1; i<=NoMedium; i++) {
    
    if ( !strcasecmp( str.c_str(), MTITP[i].label.c_str() ) )
    {
      Fill_Fluid = i;
      break;
    }
  }
  
  // チェック
  if ( Fill_Fluid == 0 )
  {
    Hostonly_ stamped_printf("\tError : Medium '%s' is not listed in '%s'\n", str.c_str(), label.c_str());
    Exit(0);
  }
  
  
  // フィルの固体媒質番号の指定
  label = "/Steer/Geometry_Model/Solid_Medium_for_Fill";
  
  if ( !(tpCntl->GetValue(label, &str )) )
  {
    Hostonly_ stamped_printf("\tParsing error : Invalid value for '%s'\n", label.c_str());
    Exit(0);
  }
  
  // ラベル名が媒質リストにあるか否かを確認
  for (int i=1; i<=NoMedium; i++) {
    
    if ( !strcasecmp( str.c_str(), MTITP[i].label.c_str() ) )
    {
      Fill_Solid = i;
      break;
    }
  }
  
  // チェック
  if ( Fill_Solid == 0 )
  {
    Hostonly_ stamped_printf("\tError : Medium '%s' in not listed in '%s'\n", str.c_str(), label.c_str());
    Exit(0);
  }
  
  
  // 流体セルのフィルのヒント
  label = "/Steer/Geometry_Model/Hint_of_Filling_Fluid";
  
  if ( !(tpCntl->GetValue(label, &str )) )
  {
    Hostonly_ stamped_printf("\tParsing error : Invalid value in '%s'\n", label.c_str());
    Exit(0);
  }
  else
  {
    if     ( !strcasecmp(str.c_str(), "no" ) )      Fill_Hint = -1;
    else if( !strcasecmp(str.c_str(), "x_minus" ) ) Fill_Hint = X_MINUS;
    else if( !strcasecmp(str.c_str(), "x_plus" ) )  Fill_Hint = X_PLUS;
    else if( !strcasecmp(str.c_str(), "y_minus" ) ) Fill_Hint = Y_MINUS;
    else if( !strcasecmp(str.c_str(), "y_plus" ) )  Fill_Hint = Y_PLUS;
    else if( !strcasecmp(str.c_str(), "z_minus" ) ) Fill_Hint = Z_MINUS;
    else if( !strcasecmp(str.c_str(), "z_plus" ) )  Fill_Hint = Z_PLUS;
    else
    {
      Hostonly_ stamped_printf("\tInvalid keyword is described for '%s'\n", label.c_str());
      Exit(0);
    }
  }
  

  // スケーリングファクター
  REAL_TYPE ct=0.0;
  
  label = "/Steer/Geometry_Model/Scaling_Factor";
  
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
// 反復関連の情報を取得する
void Control::get_Iteration(ItrCtl* IC)
{
  string str;
  string label;
  
  label = "/Steer/Iteration";
  
  // Iteration
  if( !tpCntl->chkNode(label) )
  {
    Hostonly_ stamped_printf("\tParsing error : Missing the section of '%s'\n", label.c_str());
    Exit(0);
  }
  
  label = "/Steer/Iteration/Flow";
  
  // Flow
  if( !tpCntl->chkNode(label) )
  {
    Hostonly_ stamped_printf("\tParsing error : Missing the section of '%s'\n", label.c_str());
    Exit(0);
  }
  
  
  // Flow
  switch (AlgorithmF)
  {
    case Flow_FS_EE_EE:
    case Flow_FS_AB2:
      findCriteria("/Steer/Iteration/Flow/Poisson", ItrCtl::ic_prs_pr, &IC[ItrCtl::ic_prs_pr]);
      findCriteria("/Steer/Iteration/Flow/VP",      ItrCtl::ic_div,    &IC[ItrCtl::ic_div]);
      break;
      
    case Flow_FS_AB_CN:
      findCriteria("/Steer/Iteration/Flow/Poisson", ItrCtl::ic_prs_pr, &IC[ItrCtl::ic_prs_pr]);
      findCriteria("/Steer/Iteration/Flow/NS_CN",   ItrCtl::ic_vis_cn, &IC[ItrCtl::ic_vis_cn]);
      findCriteria("/Steer/Iteration/Flow/VP",      ItrCtl::ic_div,    &IC[ItrCtl::ic_div]);
      break;
      
    case Flow_FS_RK_CN:
      findCriteria("/Steer/Iteration/Flow/Poisson",     ItrCtl::ic_prs_pr, &IC[ItrCtl::ic_prs_pr]);
      findCriteria("/Steer/Iteration/Flow/Poisson_2nd", ItrCtl::ic_prs_cr, &IC[ItrCtl::ic_prs_cr]);
      findCriteria("/Steer/Iteration/Flow/NS_CN",       ItrCtl::ic_vis_cn, &IC[ItrCtl::ic_vis_cn]);
      findCriteria("/Steer/Iteration/Flow/VP",          ItrCtl::ic_div,    &IC[ItrCtl::ic_div]);
      break;
      
    default:
      Hostonly_ stamped_printf("\tSomething wrong in '%s'\n", label.c_str());
      Exit(0);
  }
  
  // Heat
  label = "/Steer/Iteration/Heat";
  
  if ( isHeatProblem() )
  {
    if( !tpCntl->chkNode(label) )
    {
      Hostonly_ stamped_printf("\tParsing error : Missing the section of '%s'\n", label.c_str());
      Exit(0);
    }
    
    switch (AlgorithmH)
    {
      case Heat_EE_EE:
        break;
        
      case Heat_EE_EI:
        findCriteria("/Steer/Iteration/Heat/Euler_Implicit", ItrCtl::ic_tdf_ei, &IC[ItrCtl::ic_tdf_ei]);
        break;
        
      default:
        Hostonly_ stamped_printf("\tSomething wrong in '%s'\n", label.c_str());
        Exit(0);
    }
  }
}



// #################################################################
// ソルバーの計算対象種別と浮力モードを取得
void Control::get_KindOfSolver()
{
  string str;
  string label;
  
  label="/Steer/Solver_Property/Kind_of_solver";
  
  if ( !(tpCntl->GetValue(label, &str )) )
  {
	  Hostonly_ stamped_printf("\tParsing error : Invalid char* value for '/Steer/Solver_Property/Kind_of_solver\n");
	  Exit(0);
  }
  
  if     ( !strcasecmp(str.c_str(), "Flow_Only" ) )                KindOfSolver = FLOW_ONLY;
  else if( !strcasecmp(str.c_str(), "Thermal_Flow" ) )             KindOfSolver = THERMAL_FLOW;
  else if( !strcasecmp(str.c_str(), "Thermal_Flow_Natural" ) )     KindOfSolver = THERMAL_FLOW_NATURAL;
  else if( !strcasecmp(str.c_str(), "Conjugate_Heat_Transfer" ) )  KindOfSolver = CONJUGATE_HEAT_TRANSFER;
  else if( !strcasecmp(str.c_str(), "Solid_Conduction" ) )         KindOfSolver = SOLID_CONDUCTION;
  else
  {
    Hostonly_ stamped_printf("\tInvalid keyword is described for '/Steer/Solver_Property/Kind_of_solver'\n");
    Exit(0);
  }
  
  // Buoyancy option
  if ( (KindOfSolver==THERMAL_FLOW) || (KindOfSolver==THERMAL_FLOW_NATURAL) || (KindOfSolver==CONJUGATE_HEAT_TRANSFER) )
  {
    label="/Steer/Solver_Property/Buoyancy";
    
    if ( !(tpCntl->GetValue(label, &str )) )
    {
      Hostonly_ stamped_printf("\tParsing error : Invalid char* value for '/Steer/Solver_Property/Buoyancy'\n");
      Exit(0);
    }
    
    if     ( !strcasecmp(str.c_str(), "Boussinesq" ) )   Mode.Buoyancy = BOUSSINESQ;
    else if( !strcasecmp(str.c_str(), "Low_Mach" ) )     Mode.Buoyancy = LOW_MACH;
    else if( !strcasecmp(str.c_str(), "No_Buoyancy" ) )  Mode.Buoyancy = NO_BUOYANCY;
    else
    {
      Hostonly_ stamped_printf("\tInvalid keyword is described for '/Steer/Solver_Property/Buoyancy'\n");
      Exit(0);
    }
  }
}



// #################################################################
// LES計算のオプションを取得する
void Control::get_LES_option()
{
  REAL_TYPE ct;
  string str;
  string label;
  
  // 計算オプション
  label="/Steer/LES_Option/LES_Calculation";
  
  if ( !(tpCntl->GetValue(label, &str )) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '/Steer/LES_Option/LES_Calculation'\n");
	  Exit(0);
  }
  
  if     ( !strcasecmp(str.c_str(), "on") )    LES.Calc = ON;
  else if( !strcasecmp(str.c_str(), "off") )   LES.Calc = OFF;
  else
  {
    Hostonly_ stamped_printf("\tParsing error : Invalid keyword for '/Steer/LES_Option/LES_Calculation'\n");
    Exit(0);
  }
  
  if ( LES.Calc == OFF ) return;
  
  
  // モデル
  label="/Steer/LES_Option/Model";
  
  if ( !(tpCntl->GetValue(label, &str )) )
  {
	  Hostonly_ stamped_printf("\tParsing error : fail to get '/Steer/LES_Option/Model'\n");
	  Exit(0);
  }
  
  if      ( !strcasecmp(str.c_str(), "smagorinsky") )  LES.Model = Smagorinsky;
  else if ( !strcasecmp(str.c_str(), "Low_Reynolds") ) LES.Model = Low_Reynolds;
  else if ( !strcasecmp(str.c_str(), "Dynamic") )      LES.Model = Dynamic;
  else
  {
    Hostonly_ stamped_printf("\tParsing error : Invalid keyword for '/Steer/LES_Option/Model'\n");
    Exit(0);
  }
  
  
  // Cs係数
  label="/Steer/LES_Option/Cs";
  if ( !(tpCntl->GetValue(label, &ct )) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '/Steer/LES_Option/Cs'\n");
    Exit(0);
  }
  LES.Cs = ct;
  
  // Cs係数
  label="/Steer/LES_Option/Damping_Factor";
  if ( !(tpCntl->GetValue(label, &ct )) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '/Steer/LES_Option/Damping_Factor'\n");
    Exit(0);
  }
  LES.damping_factor = ct;
  
}



// #################################################################
// ログ出力モードを取得
// インターバルパラメータは，setParameters()で無次元して保持
void Control::get_Log()
{
  REAL_TYPE f_val=0.0;
  string str;
  string label;
  
  // 出力単位
  label="/Steer/Log/Unit_of_Log";
  
  if ( !(tpCntl->GetValue(label, &str )) )
  {
	  Hostonly_ stamped_printf("\tParsing error : Invalid string for '/Steer/Log/Unit_of_Log'\n");
	  Exit(0);
  }
  
  if     ( !strcasecmp(str.c_str(), "Dimensional") )      Unit.Log = DIMENSIONAL;
  else if( !strcasecmp(str.c_str(), "Non_Dimensional") )  Unit.Log = NONDIMENSIONAL;
  else
  {
    Hostonly_ stamped_printf("\tInvalid keyword is described at 'Unit_of_Log' section\n");
    Exit(0);
  }
  
  
  // Log_Base
  label="/Steer/Log/Log_Base";
  
  if ( !(tpCntl->GetValue(label, &str )) )
  {
	  Hostonly_ stamped_printf("\tParsing error : Invalid string for '/Steer/Log/Log_Base'\n");
	  Exit(0);
  }
  
  if     ( !strcasecmp(str.c_str(), "on") )   Mode.Log_Base = ON;
  else if( !strcasecmp(str.c_str(), "off") )  Mode.Log_Base = OFF;
  else
  {
    Hostonly_ stamped_printf("\tInvalid keyword is described for '/Steer/Log/Log_Base'\n");
    Exit(0);
  }
  
  
  // Log_Iteration
  label="/Steer/Log/Log_Iteration";
  
  if ( !(tpCntl->GetValue(label, &str )) )
  {
	  Hostonly_ stamped_printf("\tParsing error : Invalid string for '/Steer/Log/Log_Iteration'\n");
	  Exit(0);
  }
  
  if     ( !strcasecmp(str.c_str(), "on") )   Mode.Log_Itr = ON;
  else if( !strcasecmp(str.c_str(), "off") )  Mode.Log_Itr = OFF;
  else
  {
    Hostonly_ stamped_printf("\tInvalid keyword is described for '/Steer/Log/Log_Iteration'\n");
    Exit(0);
  }
  
  
  // Log_Wall_Info
  if ( Mode.Wall_profile == Log_Law )
  {
	  label="/Steer/Log/Log_Wall_Info";
    
	  if ( !(tpCntl->GetValue(label, &str )) )
    {
		  Hostonly_ stamped_printf("\tParsing error : Invalid string for '/Steer/Log/Log_Wall_Info'\n");
		  Exit(0);
	  }
	  
	  if     ( !strcasecmp(str.c_str(), "on") )   Mode.Log_Wall = ON;
	  else if( !strcasecmp(str.c_str(), "off") )  Mode.Log_Wall = OFF;
	  else
    {
		  Hostonly_ stamped_printf("\tInvalid keyword is described for '/Steer/Log/Log_Wall_Info'\n");
		  Exit(0);
	  }
  }
  
  
  // Log_Profiling
  label="/Steer/Log/Log_Profiling";
  
  if ( !(tpCntl->GetValue(label, &str )) )
  {
	  Hostonly_ stamped_printf("\tParsing error : Invalid string for '/Steer/Log/Log_Profiling'\n");
	  Exit(0);
  }
  
  if     ( !strcasecmp(str.c_str(), "on") )     Mode.Profiling = ON;
  else if( !strcasecmp(str.c_str(), "off") )    Mode.Profiling = OFF;
  else if( !strcasecmp(str.c_str(), "detail") ) Mode.Profiling = DETAIL;
  else
  {
    Hostonly_ stamped_printf("\tInvalid keyword is described for '/Steer/Log/Log_Profiling'\n");
    Exit(0);
  }
  
  // Interval console
  label="/Steer/Log/Console_Interval_Type";
  
  if ( !(tpCntl->GetValue(label, &str )) )
  {
	  Hostonly_ stamped_printf("\tParsing error : fail to get '/Steer/Log/Console_Interval_Type'\n");
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
		  Hostonly_ stamped_printf("\tParsing error : Invalid keyword for '/Steer/Log/Console_Interval_Type'\n");
		  Exit(0);
	  }
	  
	  label="/Steer/Log/Console_Interval";
    
	  if ( !(tpCntl->GetValue(label, &f_val )) )
    {
		  Hostonly_ stamped_printf("\tParsing error : fail to get '/Steer/Log/Console_Interval'\n");
		  Exit(0);
	  }
	  else
    {
		  Interval[Interval_Manager::tg_console].setInterval((double)f_val);
	  }
  }
  
  // Interval file_history
  label="/Steer/Log/History_Interval_Type";
  
  if ( !(tpCntl->GetValue(label, &str )) )
  {
	  Hostonly_ stamped_printf("\tParsing error : fail to get '/Steer/Log/History_Interval_Type'\n");
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
		  Hostonly_ stamped_printf("\tParsing error : Invalid keyword for '/Steer/Log/History_Interval_Type'\n");
		  Exit(0);
	  }
    
	  label="/Steer/Log/History_Interval";
    
	  if ( !(tpCntl->GetValue(label, &f_val )) )
    {
		  Hostonly_ stamped_printf("\tParsing error : fail to get '/Steer/Log/History_Interval'\n");
		  Exit(0);
	  }
	  else
    {
		  Interval[Interval_Manager::tg_history].setInterval((double)f_val);
	  }
  }
}



// #################################################################
// ノルムのラベルを返す
string Control::getNormString(const int d)
{
  string nrm;
	
	if      (d == ItrCtl::v_div_dbg) nrm = "Max. Norm : Divergence of velocity with Monitoring  ### Forced to be selected since Iteration Log is specified ###";
  else if (d == ItrCtl::v_div_max) nrm = "Max. Norm : Divergence of velocity";
  else if (d == ItrCtl::dx_b)      nrm = "dx_b : Increment of vector x divided by RHS vector b";
  else if (d == ItrCtl::r_b)       nrm = "r_b  : Residual vector divided by RHS vector b";
	else if (d == ItrCtl::r_r0)      nrm = "r_r0 : Residual vector divided by initial residual vector";
	
  return nrm;
}


// #################################################################
// 作業者情報の取得
void Control::get_Operator()
{
  string str;
  string label;
  
  label = "/Steer/Operator";
  
  if ( !(tpCntl->GetValue(label, &str )) )
  {
	  Exit(0);
  }
  
  OperatorName = str;
  
}



// #################################################################
// 無次元パラメータを各種モードに応じて設定する
// 純強制対流　有次元　（代表長さ，代表速度，動粘性係数，温度拡散係数）
//           無次元　（Pr, Re > RefV=RefL=1）
// @see bool Control::setParameters(MediumList* mat, CompoList* cmp)
void Control::get_Para_ND()
{
  REAL_TYPE ct;
  string label;
  
  label="/Parameter/Reference/Reynolds";
  
  if ( !(tpCntl->GetValue(label, &ct )) )
  {
    Hostonly_ printf("\tParsing error in '%s'\n", label.c_str());
	  Exit(0);
  }
  Reynolds=(REAL_TYPE)ct;
  
  
  label="/Parameter/Reference/Prandtl";
  
  if ( !(tpCntl->GetValue(label, &ct )) )
  {
    Hostonly_ printf("\tParsing error in '%s'\n", label.c_str());
	  Exit(0);
  }
  Prandtl=(REAL_TYPE)ct;
}



// #################################################################
// 参照パラメータを取得
void Control::get_Para_Ref()
{
  REAL_TYPE ct2;
  string label, str;
  
  label = "/Parameter/Reference/Length";
  if ( !(tpCntl->GetValue(label, &ct2 )) )
  {
    Hostonly_ printf("\tParsing error in '%s'\n", label.c_str());
    Exit(0);
  }
  RefLength = ct2;
  
  label = "/Parameter/Reference/Velocity";
  if ( !(tpCntl->GetValue(label, &ct2 )) )
  {
    Hostonly_ printf("\tParsing error in '%s'\n", label.c_str());
    Exit(0);
  }
  RefVelocity = ct2;
  
  label = "/Parameter/Reference/Gravity";
  if ( !(tpCntl->GetValue(label, &ct2 )) )
  {
    Hostonly_ printf("\tParsing error in '%s'\n", label.c_str());
    Exit(0);
  }
  Gravity = ct2;
  
  label = "/Parameter/Reference/Base_Pressure";
  if ( !(tpCntl->GetValue(label, &ct2 )) )
  {
    Hostonly_ printf("\tParsing error in '%s'\n", label.c_str());
    Exit(0);
  }
  BasePrs = ct2;

  label = "/Parameter/Reference/Medium";
  if ( !tpCntl->GetValue(label, &str) )
  {
    Hostonly_ printf("\tParsing error in '%s'\n", label.c_str());
	  Exit(0);
  }
  Ref_Medium = str;
}



// #################################################################
// 温度の参照パラメータを取得
void Control::get_Para_Temp()
{
  REAL_TYPE Base, Diff;
  string label;
  
  label="/Parameter/Temperature/Base";
  
  if ( !(tpCntl->GetValue(label, &Base )) )
  {
    Hostonly_ stamped_printf("\tParsing error : Invalid float value for '%s'\n", label.c_str());
	  Exit(0);
  }
  
  
  label="/Parameter/Temperature/Difference";
  
  if ( !(tpCntl->GetValue(label, &Diff )) )
  {
    Hostonly_ stamped_printf("\tParsing error : Invalid float value for '%s'\n", label.c_str());
	  Exit(0);
  }
  
  
  if( Diff < 0.0f )
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



// #################################################################
// PLOT3Dファイル入出力に関するパラメータ
void Control::get_PLOT3D(FileIO_PLOT3D_READ*  FP3DR, FileIO_PLOT3D_WRITE* FP3DW)
{
  string str;
  string label;
  
  // File_plot3d_filename --- option
  label = "/Steer/plot3d_options/Filename";
  
  if ( !(tpCntl->GetValue(label, &str)) )
  {
    P3Op.basename = "PLOT3Doutput_";
  }
  else
  {
    P3Op.basename = str;
  }
  if ( P3Op.basename.empty() )
  {
    P3Op.basename = "PLOT3Doutput_";
  }
  
  // File_plot3d_grid_kind
  label = "/Steer/plot3d_options/Grid_kind";
  
  if ( !(tpCntl->GetValue(label, &str )) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    Exit(0);
  }
  else
  {
    if     ( !strcasecmp(str.c_str(), "single_grid") ) FP3DR->setSingleGrid();
    else if( !strcasecmp(str.c_str(), "multi_grid") )  FP3DR->setMultiGrid();
    else
    {
      Hostonly_ stamped_printf("\tInvalid keyword is described for '%s'\n", label.c_str());
      Exit(0);
    }
  }
  
  
  // 格子の移動
  label = "/Steer/plot3d_options/Grid_Mobility";
  
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
  label = "/Steer/plot3d_options/state_of_time";
  
  if ( !(tpCntl->GetValue(label, &str )) )
  {
    //Hostonly_ stamped_printf("\tParsing error : fail to get '/Steer/plot3d_options/state_of_time'\n");
    //Exit(0);
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
  label = "/Steer/plot3d_options/set_iblank_flag";
  
  if ( !(tpCntl->GetValue(label, &str )) )
  {
    //Hostonly_ stamped_printf("\tParsing error : fail to get '/Steer/plot3d_options/set_iblank_flag'\n");
    //Exit(0);
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
  label = "/Steer/plot3d_options/dimension";
  if ( !(tpCntl->GetValue(label, &str )) )
  {
    //Hostonly_ stamped_printf("\tParsing error : fail to get '/Steer/plot3d_options/dimension'\n");
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
  
  
  //File_plot3d_format
  label = "/Steer/plot3d_options/Format_Type";
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
  label = "/Steer/plot3d_options/real_type";
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
  
  
  // Output_xyz
  label = "/Steer/plot3d_options/output_xyz";
  
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
  
  
  // Output_q
  
  /*
  label = "/Steer/plot3d_options/Output_q";
   
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
  }*/
  
  P3Op.IS_q = OFF; // 常にoff
  
  
  // Output_function
  label = "/Steer/plot3d_options/Output_function";
  
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
  
  
  // Output_function_name
  label = "/Steer/plot3d_options/Output_func_name";
  
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
  
  
  // Output_fvbnd
  
  /*
  label = "/Steer/plot3d_options/Output_fvbnd";
   
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
  }*/
  
  P3Op.IS_fvbnd = OFF; // 常にoff
  
  // divide_func ---> 出力を項目別にファイル分割するオプション
  label = "/Steer/plot3d_options/divide_func";
  
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



// #################################################################
// 性能試験モードを取得する（隠しパラメータ）
// 'Performance_Test'の文字列チェックはしないので注意して使うこと
void Control::get_PMtest()
{
  string str;
  string label;
  
  label="/Steer/Performance_Test";
  
  if ( !(tpCntl->GetValue(label, &str )) )
  {
	  return;
  }
  
  if     ( !strcasecmp(str.c_str(), "on") )   Hide.PM_Test = ON;
  else if( !strcasecmp(str.c_str(), "off") )  Hide.PM_Test = OFF;
  else
  {
    Hostonly_ stamped_printf("\tInvalid keyword is described for '/Steer/Performance_Test'\n");
    Exit(0);
  }
}


// #################################################################
// 参照座標系を取得する
void Control::get_ReferenceFrame(ReferenceFrame* RF)
{
  string str;
  string label;
  
  label="/Steer/Reference_Frame/Reference_Frame_Type";
  
  if ( !(tpCntl->GetValue(label, &str )) ) {
    Hostonly_ stamped_printf("\tParsing error : fail to get '/Steer/Reference_Frame/Reference_Frame_Type'\n");
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
      Hostonly_ stamped_printf("\tParsing error : Invalid values for Translational Reference_Frame\n");
      Exit(0);
    }
    RF->setGridVel((double*)xyz);
  }
  else
  {
    Hostonly_ stamped_printf("\tParsing error : Invalid keyword for '/Steer/Reference_Frame/Reference_Frame_Type'\n");
    Exit(0);
  }
}



// #################################################################
// モニタリングのON/OFFとセルモニタの有無のみを取得　詳細パラメータは後ほど
void Control::get_Sampling()
{
  string str;
  string label;
  
  // ログ出力
  label = "/Steer/Monitor_List/log";
  
  if ( !(tpCntl->GetValue(label, &str )) )
  {
	  Hostonly_ stamped_printf("\tParsing error : Invalid string for '/Steer/Monitor_List/Log'\n");
	  Exit(0);
  }
  
  if     ( !strcasecmp(str.c_str(), "on") )   Sampling.log = ON;
  else if( !strcasecmp(str.c_str(), "off") )  Sampling.log = OFF;
  else
  {
    Hostonly_ stamped_printf("\tInvalid keyword is described for '/Steer/Monitor_List/Log'\n");
    Exit(0);
  }
  
  // セルモニターの指定
  /// @see Initialize.C setEnsComponent()
  if ( Sampling.log == ON )
  {
	  label = "/Steer/Monitor_List/Cell_Monitor";
    
	  if ( !(tpCntl->GetValue(label, &str )) )
    {
		  Hostonly_ stamped_printf("\tParsing error : Invalid string for '/Steer/Monitor_List/Cell_Monitor'\n");
		  Exit(0);
	  }
    
	  if     ( !strcasecmp(str.c_str(), "on") )   EnsCompo.monitor = ON;
	  else if( !strcasecmp(str.c_str(), "off") )  EnsCompo.monitor = OFF;
	  else
    {
      Hostonly_ stamped_printf("\tInvalid keyword is described for '/Steer/Monitor_List/Cell_Monitor'\n");
		  Exit(0);
	  }
  }
  
}



// #################################################################
// ソルバーの種類を特定するパラメータを取得し，ガイドセルの値を決定する
void Control::get_Solver_Properties()
{
  string str;
  string label;
  
  // 形状近似度の取得
  label = "/Steer/Solver_Property/Shape_Approximation";
  
  if ( !(tpCntl->GetValue(label, &str )) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '/Steer/Solver_Property/Shape_Approximation'\n");
    Exit(0);
  }
  
  if     ( !strcasecmp(str.c_str(), "Binary") )       Mode.ShapeAprx = BINARY;
  else if( !strcasecmp(str.c_str(), "cut_distance") ) Mode.ShapeAprx = CUT_INFO;
  else
  {
    Hostonly_ stamped_printf("\tInvalid keyword is described for '/Steer/Solver_Property/Shape_Approximation'\n");
    Exit(0);
  }
  
  
  // 支配方程式の型（PDE_NS / Euler）を取得
  label = "/Steer/Solver_Property/PDE_type";
  
  if ( !(tpCntl->GetValue(label, &str )) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '/Steer/Solver_Property/PDE_type'\n");
    Exit(0);
  }
  
  if     ( !strcasecmp(str.c_str(), "Navier_Stokes" ) ) Mode.PDE = PDE_NS;
  else if( !strcasecmp(str.c_str(), "Euler" ) )         Mode.PDE = PDE_EULER;
  else
  {
    Hostonly_ stamped_printf("\tInvalid keyword is described for '/Steer/Solver_Property/PDE_type'\n");
    Exit(0);
  }
  
  
  // 基礎方程式の種類を取得する
  label = "/Steer/Solver_Property/Basic_Equation";
  
  if ( !(tpCntl->GetValue(label, &str )) ) {
    Hostonly_ stamped_printf("\tParsing error : fail to get '/Steer/Solver_Property/Basic_Equation'\n");
    Exit(0);
  }
  
  if     ( !strcasecmp(str.c_str(), "Incompressible" ) )            BasicEqs = INCMP;
  else if( !strcasecmp(str.c_str(), "Limited_Compressible" ) )      BasicEqs = LTDCMP;
  else if( !strcasecmp(str.c_str(), "Compressible" ) )              BasicEqs = CMPRSS;
  else if( !strcasecmp(str.c_str(), "Incompressible_Two_Phase" ) )  BasicEqs = INCMP_2PHASE;
  else
  {
    Hostonly_ stamped_printf("\tInvalid keyword is described for '/Steer/Solver_Property/Basic_Equation'\n");
    Exit(0);
  }
  
  
  // 非定常計算，または定常計算の種別を取得する
  label = "/Steer/Solver_Property/Time_Variation";
  
  if ( !(tpCntl->GetValue(label, &str )) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '/Steer/Solver_Property/Time_Variation'\n");
    Exit(0);
  }
  
  if     ( !strcasecmp(str.c_str(), "Steady" ) )    Mode.Steady = FB_STEADY;
  else if( !strcasecmp(str.c_str(), "Unsteady" ) )  Mode.Steady = FB_UNSTEADY;
  else
  {
    Hostonly_ stamped_printf("\tInvalid keyword is described for '/Steer/Solver_Property/Time_Variation'\n");
    Exit(0);
  }
  
  
  // 対流項スキームの種類の取得
  get_Convection();
  
  
  // ソルバーの種類（FLOW_ONLY / THERMAL_FLOW / THERMAL_FLOW_NATURAL / CONJUGATE_HEAT_TRANSFER / SOLID_CONDUCTION）と浮力モード
  get_KindOfSolver();
  
  
  // ガイドセルの値を決める get_Convection(), get_KindOfSolver()のあと
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
  
  
  // 平均値の引き戻しオプション
  label = "/Steer/Solver_Property/Pressure_Shift";
  
  if ( !(tpCntl->GetValue(label, &str )) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
	  Exit(0);
  }
  
  if     ( !strcasecmp(str.c_str(), "off" ) )     Mode.Pshift = -1;
  else if( !strcasecmp(str.c_str(), "x_minus" ) ) Mode.Pshift = X_MINUS;
  else if( !strcasecmp(str.c_str(), "x_plus" ) )  Mode.Pshift = X_PLUS;
  else if( !strcasecmp(str.c_str(), "y_minus" ) ) Mode.Pshift = Y_MINUS;
  else if( !strcasecmp(str.c_str(), "y_plus" ) )  Mode.Pshift = Y_PLUS;
  else if( !strcasecmp(str.c_str(), "z_minus" ) ) Mode.Pshift = Z_MINUS;
  else if( !strcasecmp(str.c_str(), "z_plus" ) )  Mode.Pshift = Z_PLUS;
  else
  {
    stamped_printf("\tInvalid keyword is described for '%s'\n", label.c_str());
    Exit(0);
  }
}




// #################################################################
//初期値とリスタート条件
void Control::get_start_condition()
{
  int ct;
  string str;
  string label;
  
  label="/Steer/Start_condition/start_type";
  
  if ( !(tpCntl->GetValue(label, &str )) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    Exit(0);
  }
  
  if      ( !strcasecmp(str.c_str(), "initial") )                  Start = initial_start;
  else if ( !strcasecmp(str.c_str(), "restart") )                  Start = restart;
  else if ( !strcasecmp(str.c_str(), "restart_from_coarse_data") ) Start = coarse_restart;
  else if ( !strcasecmp(str.c_str(), "restart_different_nproc") )  Start = restart_different_nproc;
  else
  {
    Hostonly_ stamped_printf("\tInvalid keyword is described for '%s'\n", label.c_str());
    Exit(0);
  }
  
  
  // リスタート時のタイムスタンプ
  if ( (Start == restart) || (Start == coarse_restart) || (Start == restart_different_nproc) )
  {
    label="/Steer/Start_condition/restart_step";
    
    if ( !(tpCntl->GetValue(label, &ct )) )
    {
      Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
      Exit(0);
    }
    Restart_step = ct;
    
  }
  
  // リスタート時のタイムスタンプ Coarse_restart
  if ( Start == coarse_restart )
  {
    label="/Steer/Start_condition/Prefix_of_Coarse_Pressure";
    
    if ( !(tpCntl->GetValue(label, &str )) )
    {
      Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
      Exit(0);
    }
    
    f_Coarse_pressure = str.c_str();
    
    
    label="/Steer/Start_condition/Prefix_of_Coarse_Velocity";
    
    if ( !(tpCntl->GetValue(label, &str )) )
    {
      Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
      Exit(0);
    }
    
    f_Coarse_velocity = str.c_str();
    
    
    if ( isHeatProblem() )
    {
      label="/Steer/Start_condition/Prefix_of_Coarse_Temperature";
      
      if ( !(tpCntl->GetValue(label, &str )) )
      {
        Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
        Exit(0);
      }
      
      f_Coarse_temperature = str.c_str();
    }
    
  }
  
  
  // 粗い格子を使う場合の初期値データのプリフィクス
  if ( Start == coarse_restart )
  {
    // プロセス並列時にローカルでのファイル入力を指定した場合
    if ( FIO.IO_Input == IO_DISTRIBUTE )
    {
      label="/Steer/Start_condition/dfi_file_pressure";
      
      if ( !(tpCntl->GetValue(label, &str )) )
      {
        Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
        Exit(0);
      }
      
      f_Coarse_dfi_prs = str.c_str();
      
      
      label="/Steer/Start_condition/dfi_file_velocity";
      
      if ( !(tpCntl->GetValue(label, &str )) )
      {
        Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
        Exit(0);
      }
      
      f_Coarse_dfi_vel = str.c_str();
      
      
      if ( isHeatProblem() )
      {
        label="/Steer/Start_condition/dfi_file_temperature";
        
        if ( !(tpCntl->GetValue(label, &str )) )
        {
          Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
          Exit(0);
        }
        
        f_Coarse_dfi_temp = str.c_str();
        
      }
    }
    
  }
  
  
  // 異なる並列数からリスタートする場合の初期値データのプリフィクス
  if ( Start == restart_different_nproc )
  {
    
#ifdef _STAGING_
    
    label="/Steer/Start_condition/Prefix_of_different_nproc_dir";
    
    if ( !(tpCntl->GetValue(label, &str )) )
    {
      Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
      Exit(0);
    }
    
    f_different_nproc_dir_prefix = str.c_str();
    
#endif
    
    label="/Steer/Start_condition/Prefix_of_different_nproc_Pressure";
    
    if ( !(tpCntl->GetValue(label, &str )) )
    {
      Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
      Exit(0);
    }
    
    f_different_nproc_pressure = str.c_str();
    
    
    label="/Steer/Start_condition/Prefix_of_different_nproc_Velocity";
    
    if ( !(tpCntl->GetValue(label, &str )) )
    {
      Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
      Exit(0);
    }
    
    f_different_nproc_velocity = str.c_str();
    
    
    if ( isHeatProblem() )
    {
      label="/Steer/Start_condition/Prefix_of_different_nproc_Temperature";
      
      if ( !(tpCntl->GetValue(label, &str )) )
      {
        Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
        Exit(0);
      }
      
      f_different_nproc_temperature = str.c_str();
    }
    
    // プロセス並列時にローカルでのファイル入力を指定した場合 --- always IO_DISTRIBUTE
    if ( FIO.IO_Input == IO_DISTRIBUTE )
    {
      label="/Steer/Start_condition/dfi_file_pressure";
      
      if ( !(tpCntl->GetValue(label, &str )) )
      {
        Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
        Exit(0);
      }
      
      f_different_nproc_dfi_prs = str.c_str();
      
      
      label="/Steer/Start_condition/dfi_file_velocity";
      
      if ( !(tpCntl->GetValue(label, &str )) )
      {
        Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
        Exit(0);
      }
      
      f_different_nproc_dfi_vel = str.c_str();
      
      
      if ( isHeatProblem() )
      {
        label="/Steer/Start_condition/dfi_file_temperature";
        
        if ( !(tpCntl->GetValue(label, &str )) )
        {
          Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
          Exit(0);
        }
        
        f_different_nproc_dfi_temp = str.c_str();
        
      }
    }
    
  }
  
  
  // 初期条件
  if ( Start == initial_start )
  {
    // Density
    label="/Steer/Start_condition/Initial_State/Density";
    
    if ( !(tpCntl->GetValue(label, &iv.Density )) )
    {
      Hostonly_ stamped_printf("\tParsing error : Invalid float value for '%s'\n", label.c_str());
      Exit(0);
    }
    
    // Pressure
    label="/Steer/Start_condition/Initial_State/Pressure";
    
    if ( !(tpCntl->GetValue(label, &iv.Pressure )) )
    {
      Hostonly_ stamped_printf("\tParsing error : Invalid float value for '%s'\n", label.c_str());
      Exit(0);
    }
    
    // Velocity
    REAL_TYPE v[3];
    for (int n=0; n<3; n++) v[n]=0.0;
    label="/Steer/Start_condition/Initial_State/Velocity";
    
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
      label="/Steer/Start_condition/Initial_State/Temperature";
      
      if ( !(tpCntl->GetValue(label, &iv.Temperature )) )
      {
        Hostonly_ stamped_printf("\tParsing error : Invalid float value for '%s'\n", label.c_str());
        Exit(0);
      }
      
    }
  }
  
}



// #################################################################
// 制御，計算パラメータ群の取得
void Control::get_Steer_1(DTcntl* DT, FileIO_PLOT3D_READ* FP3DR, FileIO_PLOT3D_WRITE* FP3DW)
{
  
  // ソルバーの具体的な種類を決めるパラメータを取得し，ガイドセルの値を設定する
  get_Solver_Properties();
  
  // 指定単位が有次元か無次元かを取得
  get_Unit();
  
  // Reference parameter needs to be called before setDomain();
  // パラメータの取得，代表値に関するもの．
  get_Para_Ref();
  
  // 時間制御パラメータ
  get_Time_Control(DT);
  
  // パラメータチェック
  get_CheckParameter();
  
  // モニターのON/OFF 詳細パラメータはget_Monitor()で行う
  get_Sampling();
  
  // ファイル入出力に関するパラメータ
  get_FileIO();
  
  // PLOT3Dファイル入出力に関するパラメータ
  if (FIO.PLOT3D_OUT == ON) get_PLOT3D(FP3DR,FP3DW);
  
}


// #################################################################
// 制御，計算パラメータ群の取得
void Control::get_Steer_2(ItrCtl* IC, ReferenceFrame* RF)
{
  // 流体の解法アルゴリズムを取得
  get_Algorithm();
  
  // パラメータを取得する
  if ( Unit.Param == NONDIMENSIONAL )
  {
    if ( KindOfSolver == FLOW_ONLY ) get_Para_ND();
  }
  
  if ( isHeatProblem() )
  {
    get_Para_Temp();
  }
  
  // Reference frame information
  get_ReferenceFrame(RF);
  
  // 時間平均操作
  get_Average_option();
  
  // 圧力ノイマン条件のタイプ >> get_Log()よりも先に
  get_Wall_type();
  
  // Log >> get_Iteration()よりも前に
  get_Log();
  
  // Criteria of computation
  get_Iteration(IC);
  
  // LES
  get_LES_option();
  
  
  // 派生変数のオプション
  get_Derived();
  
  // 変数範囲の処理　***隠しパラメータ
  get_VarRange();
  
  // Cell IDのゼロを指定IDに変更　***隠しパラメータ
  get_ChangeID();
  
  // 性能測定モードの処理　***隠しパラメータ
  get_PMtest();
  
  // ラフな初期値を使い、リスタートするモード指定 >> FileIO
  get_start_condition();
  
  // 作業者情報
  get_Operator();
  
}



// #################################################################
// 時間制御に関するパラメータを取得する
// パラメータは，setParameters()で無次元して保持
void Control::get_Time_Control(DTcntl* DT)
{
  REAL_TYPE ct;
  int ss=0;
  
  string str;
  string label;
  
  // 加速時間
  label = "/Steer/Time_Control/Acceleration_Type";
  
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
	  
	  label = "/Steer/Time_Control/Acceleration";
    
	  if ( !(tpCntl->GetValue(label, &ct )) )
    {
		  Hostonly_ stamped_printf("\tParsing error : fail to get '%s\n", label.c_str());
		  Exit(0);
	  }
	  else
    {
		  Interval[Interval_Manager::tg_accelra].setInterval((double)ct);
	  }
  }
  
  
  // 時間積分幅を取得する
  label = "/Steer/Time_Control/Dt_Type";
  
  if ( !(tpCntl->GetValue(label, &str )) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    Exit(0);
  }
  
  
  label = "/Steer/Time_Control/Delta_t";
  
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
      cc = (double)ct / ts;
    }
  }
  else
  {
    cc = (double)ct;
  }
  
  if ( !DT->set_Scheme(str.c_str(), cc) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to set DELTA_T\n");
    Exit(0);
  }
  
  // 計算する時間を取得する
  label = "/Steer/Time_Control/Period_Type";
  
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
    
    
    label = "/Steer/Time_Control/Calculation_Period";
    
    if ( !(tpCntl->GetValue(label, &ct )) )
    {
      Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
      Exit(0);
    }
    else
    {
      Interval[Interval_Manager::tg_compute].setInterval((double)ct);
    }
    
  }
}



// #################################################################
// 入力ファイルに記述するパラメータとファイルの有次元・無次元の指定を取得する
void Control::get_Unit()
{
  string str;
  string label;
  
  label = "/Steer/Unit/Unit_of_Input_Parameter";
  
  if ( !(tpCntl->GetValue(label, &str )) )
  {
		Hostonly_ stamped_printf("\tParsing error : Invalid string for '%s'\n", label.c_str());
	  Exit(0);
  }
  
  if     ( !strcasecmp(str.c_str(), "Dimensional") )      Unit.Param = DIMENSIONAL;
  else if( !strcasecmp(str.c_str(), "Non_Dimensional") )  Unit.Param = NONDIMENSIONAL;
  else
  {
    Hostonly_ stamped_printf("\tInvalid keyword is described at '%s'\n", label.c_str());
    Exit(0);
  }
  
  
  label = "/Steer/Unit/pressure";
  
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
    label = "/Steer/Unit/temperature";
    
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
// 変数の範囲制限モードを取得
void Control::get_VarRange()
{
  string str;
  string label;

  label="/Steer/Variable_Range";
  
  if ( !(tpCntl->GetValue(label, &str )) )
  {
	  return;
  }
  
  if     ( !strcasecmp(str.c_str(), "normal") )   Hide.Range_Limit = Range_Normal;
  else if( !strcasecmp(str.c_str(), "cutoff") )   Hide.Range_Limit = Range_Cutoff;
  else
  {
    Hostonly_ stamped_printf("\tInvalid keyword is described for '%s\n", label.c_str());
    Exit(0);
  }
  
}



// #################################################################
// バージョン情報の取得
void Control::get_Version()
{
  int ct;
  string label;
  
  // FlowBase
  label = "/Steer/Version_Info/Flow_Base";
  
  if ( !(tpCntl->GetValue(label, &ct )) )
  {
    Hostonly_ stamped_printf("\tParsing error : Invalid value for '%s'\n", label.c_str());
    Exit(0);
  }
  FB_version = ct;
  
  // FFV
  label = "/Steer/Version_Info/FFV";
  
  if ( !(tpCntl->GetValue(label, &ct)) )
  {
    Hostonly_ stamped_printf("\tParsing error : Invalid value for '%s'\n", label.c_str());
    Exit(0);
  }
  version = ct;
}



// #################################################################
//壁面上の扱いを指定する
void Control::get_Wall_type()
{
  string str;
  string label;
  
  // 圧力のタイプ
  label="/Steer/Treatment_of_wall/Pressure_Gradient";
  
  if ( !(tpCntl->GetValue(label, &str )) )
  {
	  Hostonly_ stamped_printf("\tParsing error : Invalid string for '%s'\n", label.c_str());
	  Exit(0);
  }
  
  if     ( !strcasecmp(str.c_str(), "grad_zero") )   Mode.PrsNeuamnnType = P_GRAD_ZERO;
  else if( !strcasecmp(str.c_str(), "grad_NS") )     Mode.PrsNeuamnnType = P_GRAD_NS;
  else
  {
    Hostonly_ stamped_printf("\tInvalid keyword is described for '%s'\n", label.c_str());
    Exit(0);
  }
  
  // 壁面摩擦応力の計算モード
  label="/Steer/Treatment_of_wall/Velocity_Profile";
  
  if ( !(tpCntl->GetValue(label, &str )) )
  {
	  Hostonly_ stamped_printf("\tParsing error : Invalid string for '%s'\n", label.c_str());
	  Exit(0);
  }
  
  if     ( !strcasecmp(str.c_str(), "no_slip") )     Mode.Wall_profile = No_Slip;
  else if( !strcasecmp(str.c_str(), "Slip") )        Mode.Wall_profile = Slip;
  else if( !strcasecmp(str.c_str(), "Law_of_Wall") ) Mode.Wall_profile = Log_Law;
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
 @fn void Control::printInitValues(FILE* fp)
 @brief 初期値の表示
 @note Init*には無次元値が保持されている
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
  if (Unit.Prs == Unit_Absolute) {
    fprintf(fp,"\tInitial  Pressure    [Pa]    / [-]   : %12.5e / %12.5e\n", iv.Pressure, (iv.Pressure-BasePrs)/DynamicPrs);
  }
  else {
    fprintf(fp,"\tInitial  Pressure    [Pa_g]  / [-]   : %12.5e / %12.5e\n", iv.Pressure, iv.Pressure/DynamicPrs);
  }
  if ( isHeatProblem() ) {
    fprintf(fp,"\tInitial  Temperature [%s]     / [-]   : %12.5e / %12.5e\n", (Unit.Temp==Unit_KELVIN) ? "K" : "C",
						FBUtility::convK2Temp(iv.Temperature, Unit.Temp), 
						FBUtility::convK2ND(iv.Temperature, BaseTemp, DiffTemp));
  }
  
  fprintf(fp,"\n");
  fflush(fp);
}


// #################################################################
// 線形ソルバー種別の表示
void Control::printLS(FILE* fp, const ItrCtl* IC)
{
  switch (IC->get_LS()) 
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
// 計算パラメータの表示
void Control::printParaConditions(FILE* fp, const MediumList* mat)
{
  if ( !fp )
  {
    stamped_printf("\tFail to write into file\n");
    Exit(0);
  }
  
  fprintf(fp,"\n---------------------------------------------------------------------------\n\n");
  fprintf(fp,"\n\t>> Simulation Parameters\n\n");
  
  fprintf(fp,"\tReference ID              [-]         :  %s\n", mat[RefMat].getLabel().c_str());
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
void Control::printSteerConditions(FILE* fp, const ItrCtl* IC, const DTcntl* DT, const ReferenceFrame* RF, FileIO_PLOT3D_WRITE* FP3DW)
{
  if ( !fp )
  {
    stamped_printf("\tFail to write into file\n");
    Exit(0);
  }
  
  REAL_TYPE dt = (REAL_TYPE)DT->get_DT();
  bool  err=true;
  double itm=0.0;
  
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
  
  // Pressure shift
  if (Mode.Pshift == -1)
  {
    fprintf(fp,"\t     Pressure Shift           :   Off\n");
  }
  else
  {
    fprintf(fp,"\t     Pressure Shift           :   %s\n", FBUtility::getDirection(Mode.Pshift).c_str());
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
    if ( !Interval[Interval_Manager::tg_avstart].isStep() )
    {
      itm = Interval[Interval_Manager::tg_avstart].getIntervalTime();
      fprintf(fp,"\t     Averaging Operation      :   %12.5e [sec] / %12.5e [-]\n", itm*Tscale, itm);
    }
    else
    {
      fprintf(fp,"\t     Averaging Operation      :   %12d\n", Interval[Interval_Manager::tg_avstart].getIntervalStep());
    }
  }
  else
  {
    fprintf(fp,"\t     Averaging Operation      :   OFF\n");
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
      fprintf(fp,"\t     Time Increment dt        :   %12.5e [sec] / %12.5e [-] : CFL (/w Sound Speed) & Diffusion number with Maximum velocity\n", dt*Tscale, dt);
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
      
    case restart:
      fprintf(fp,"\t     Start Condition          :   Restart from previous session\n");
      break;
      
    case coarse_restart:
      fprintf(fp,"\t     Start Condition          :   Restart from coarse grid data\n");
      break;
      
    case restart_different_nproc:
      fprintf(fp,"\t     Start Condition          :   Restart from previous session that nproc differ from\n");
      break;
      
    default:
      stamped_printf("Error: start condition section\n");
      err=false;
  }
  
  // 粗い格子の計算結果を使ったリスタート
  if ( Start == coarse_restart )
  {
    if ( FIO.IO_Input == IO_GATHER )
    {
      fprintf(fp,"\t     with Coarse Initial data files\n");
      fprintf(fp,"\t          Pressure            :   %s\n", f_Coarse_pressure.c_str());
      fprintf(fp,"\t          Velocity            :   %s\n", f_Coarse_velocity.c_str());
      if ( isHeatProblem() )
      {
        fprintf(fp,"\t          Temperature         :   %s\n", f_Coarse_temperature.c_str());
      }
    }
    else
    {
      fprintf(fp,"\t     with Coarse Initial data files\n");
      fprintf(fp,"\t          DFI file Pressure   :   %s\n", f_Coarse_dfi_prs.c_str());
      fprintf(fp,"\t          DFI file Velocity   :   %s\n", f_Coarse_dfi_vel.c_str());
      if ( isHeatProblem() )
      {
        fprintf(fp,"\t          DFI file Temperature:   %s\n", f_Coarse_dfi_temp.c_str());
      }
      fprintf(fp,"\t          Prefix of Pressure  :   %s\n", f_Coarse_pressure.c_str());
      fprintf(fp,"\t          Prefix of Velocity  :   %s\n", f_Coarse_velocity.c_str());
      if ( isHeatProblem() )
      {
        fprintf(fp,"\t          Prefix of Temp.     :   %s\n", f_Coarse_temperature.c_str());
      }
    }
    
  }
  
  // 異なる並列数からリスタート
  if ( Start == restart_different_nproc )
  {
    if ( FIO.IO_Input == IO_GATHER )
    {
      fprintf(fp,"\t     with different_nproc Initial data files\n");
      fprintf(fp,"\t          Pressure            :   %s\n", f_different_nproc_pressure.c_str());
      fprintf(fp,"\t          Velocity            :   %s\n", f_different_nproc_velocity.c_str());
      if ( isHeatProblem() )
      {
        fprintf(fp,"\t          Temperature         :   %s\n", f_different_nproc_temperature.c_str());
      }
    }
    else
    {
      fprintf(fp,"\t     with different_nproc Initial data files\n");
      fprintf(fp,"\t          DFI file Pressure   :   %s\n", f_different_nproc_dfi_prs.c_str());
      fprintf(fp,"\t          DFI file Velocity   :   %s\n", f_different_nproc_dfi_vel.c_str());
      if ( isHeatProblem() )
      {
        fprintf(fp,"\t          DFI file Temperature:   %s\n", f_different_nproc_dfi_temp.c_str());
      }
      fprintf(fp,"\t          Prefix of Pressure  :   %s\n", f_different_nproc_pressure.c_str());
      fprintf(fp,"\t          Prefix of Velocity  :   %s\n", f_different_nproc_velocity.c_str());
      if ( isHeatProblem() )
      {
        fprintf(fp,"\t          Prefix of Temp.     :   %s\n", f_different_nproc_temperature.c_str());
      }
    }
    
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
  
  // InputMode >> デフォルトでローカル
  //fprintf(fp,"\t     Input Mode               :   %s\n", (FIO.IO_Input==IO_GATHER) ? "Master IO" : "Local IO");
  
  
  // Output guide
  fprintf(fp,"\t     Guide cell for output    :   %d\n", GuideOut);
  
  // Voxel output
  fprintf(fp,"\t     Voxel model output       :   %s\n", (FIO.IO_Voxel==Sphere_SVX) ? "svx" : "None");
  
  
  // Output mode
  switch (FIO.IO_Mode)
  {
    case io_current:
      fprintf(fp,"\t     Output Mode              :   Current Directory\n");
      break;
      
    case io_specified:
      fprintf(fp,"\t     Output Mode              :   Specified Directory \"%s\"\n", FIO.IO_DirPath.c_str());
      break;
      
    case io_time_slice:
      fprintf(fp,"\t     Output Mode              :   TIme Slice\n");
      break;
      
    default:
      break;
  }
  
  
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
  
  // 瞬間値のファイル出力
  if ( !Interval[Interval_Manager::tg_instant].isStep() )
  {
    itm = Interval[Interval_Manager::tg_instant].getIntervalTime();
    fprintf(fp,"\t     Instant data             :   %12.6e [sec] / %12.6e [-]\n", itm*Tscale, itm);
  }
  else
  {
    fprintf(fp,"\t     Instant data             :   %12d [step]\n", Interval[Interval_Manager::tg_instant].getIntervalStep());
  }
  
  // 平均値のファイル出力
  if ( !Interval[Interval_Manager::tg_average].isStep() )
  {
    itm = Interval[Interval_Manager::tg_average].getIntervalTime();
    fprintf(fp,"\t     Averaged data            :   %12.6e [sec] / %12.6e [-]\n", itm*Tscale, itm);
  }
  else
  {
    fprintf(fp,"\t     Averaged data            :   %12d [step]\n", Interval[Interval_Manager::tg_average].getIntervalStep());
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
  if(FIO.PLOT3D_OUT == ON)
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

  
  // PLOT3D Options
  if(FIO.PLOT3D_OUT == ON)
  {
    fprintf(fp,"\n\tPLOT3D Options\n");
    fprintf(fp,"\t     file prefix              :   %s\n", P3Op.basename.c_str());
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
  

  
  // Criteria ------------------
  fprintf(fp,"\n\tParameter of Linear Equation\n");
  const ItrCtl* ICp1= &IC[ItrCtl::ic_prs_pr];  /// 圧力のPoisson反復
  const ItrCtl* ICp2= &IC[ItrCtl::ic_prs_cr];  /// 圧力のPoisson反復　2回目
  const ItrCtl* ICv = &IC[ItrCtl::ic_vis_cn];  /// 粘性項のCrank-Nicolson反復
  const ItrCtl* ICd = &IC[ItrCtl::ic_div];     /// V-P反復
  
  if ( Hide.PM_Test == ON )
  {
    fprintf(fp,"\t ### Performance Test Mode >> The iteration number is fixed by Iteration max.\n\n");
  }
  
	if ( KindOfSolver != SOLID_CONDUCTION )
  {
    // V-P iteration
		fprintf(fp,"\t     V-P Iteration \n");
		fprintf(fp,"\t       Iteration max          :   %d\n"  ,  ICd->get_ItrMax());
		fprintf(fp,"\t       Convergence eps        :   %9.3e\n", ICd->get_eps());
		fprintf(fp,"\t       Norm type              :   %s\n", getNormString(ICd->get_normType()).c_str() );
    
    
		// 1st iteration
		fprintf(fp,"\t     1st Pressure Iteration \n");
		fprintf(fp,"\t       Iteration max          :   %d\n"  ,  ICp1->get_ItrMax());
		fprintf(fp,"\t       Convergence eps        :   %9.3e\n", ICp1->get_eps());
		fprintf(fp,"\t       Coef. of Relax./Accel. :   %9.3e\n", ICp1->get_omg());
		fprintf(fp,"\t       Norm type              :   %s\n", getNormString(ICp1->get_normType()).c_str() );
    fprintf(fp,"\t       Communication Mode     :   %s\n", (ICp1->get_SyncMode()==comm_sync) ? "SYNC" : "ASYNC");
		printLS(fp, ICp1);
    
    if ( AlgorithmF == Flow_FS_RK_CN )
    {
      fprintf(fp,"\t     2nd Pressure Iteration \n");
      fprintf(fp,"\t       Iteration max          :   %d\n"  ,  ICp2->get_ItrMax());
      fprintf(fp,"\t       Convergence eps        :   %9.3e\n", ICp2->get_eps());
      fprintf(fp,"\t       Coef. of Relax./Accel. :   %9.3e\n", ICp2->get_omg());
      fprintf(fp,"\t       Norm type              :   %s\n", getNormString(ICp2->get_normType()).c_str() );
      fprintf(fp,"\t       Communication Mode     :   %s\n", (ICp1->get_SyncMode()==comm_sync) ? "SYNC" : "ASYNC");
      printLS(fp, ICp2);
    }
    
    // CN iteration
		if ( (AlgorithmF == Flow_FS_AB_CN) || (AlgorithmF == Flow_FS_RK_CN) )
    {
      fprintf(fp,"\n");
			fprintf(fp,"\t     Velocity CN Iteration \n");
			fprintf(fp,"\t       Iteration max           :   %d\n"  ,  ICv->get_ItrMax());
			fprintf(fp,"\t       Convergence eps         :   %9.3e\n", ICv->get_eps());
			fprintf(fp,"\t       Coef. of Relax./Accel.  :   %9.3e\n", ICv->get_omg());
			fprintf(fp,"\t       Norm type               :   %s\n", getNormString(ICv->get_normType()).c_str() );
      fprintf(fp,"\t       Communication Mode      :   %s\n", (ICp1->get_SyncMode()==comm_sync) ? "SYNC" : "ASYNC");
			printLS(fp, ICv);
		}
	}
	
  // for Temperature
  if ( isHeatProblem() )
  {
    if ( AlgorithmH == Heat_EE_EI )
    {
      const ItrCtl* ICt = &IC[ItrCtl::ic_tdf_ei];  /// 温度の拡散項の反復
      fprintf(fp,"\n");
      fprintf(fp,"\t     Temperature Iteration  \n");
      fprintf(fp,"\t       Iteration max          :   %d\n"  ,  ICt->get_ItrMax());
      fprintf(fp,"\t       Convergence eps        :   %9.3e\n", ICt->get_eps());
      fprintf(fp,"\t       Coef. of Relax./Accel. :   %9.3e\n", ICt->get_omg());
			fprintf(fp,"\t       Norm type              :   %s\n", getNormString(ICt->get_normType()).c_str() );
      fprintf(fp,"\t       Communication Mode     :   %s\n", (ICp1->get_SyncMode()==comm_sync) ? "SYNC" : "ASYNC");
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
  
  
  // Hidden parameter -----------------
  
  if (Hide.Change_ID != 0 )
  {
    fprintf(fp,"\t     Change ID from [0] to    :   %d\n", Hide.Change_ID);
  }
  fprintf(fp,"\n");
  
  switch (Hide.Range_Limit)
  {
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
 */
void Control::setParameters(MediumList* mat, CompoList* cmp, ReferenceFrame* RF, BoundaryOuter* BO)
{
  REAL_TYPE rho, nyu, cp, lambda, beta, mu, snd_spd=0.0;
  REAL_TYPE c1, c2, c3;
  int m;
  
  // get reference values
  for (int n=NoBC+1; n<=NoCompo; n++) {
    
    if ( cmp[n].getMatOdr() == RefMat )
    {
      m = cmp[n].getMatOdr();
      if ( mat[m].getState() == FLUID )
      {
        rho   = mat[m].P[p_density];
        mu    = mat[m].P[p_viscosity];
        nyu   = mat[m].P[p_kinematic_viscosity];
        cp    = mat[m].P[p_specific_heat];
        lambda= mat[m].P[p_thermal_conductivity];
        beta  = mat[m].P[p_vol_expansion]; // can be replaced by 1/K in the case of gas
        snd_spd = mat[m].P[p_sound_of_speed];
      }
      else
      {
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
	else if (KindOfSolver==THERMAL_FLOW) {
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
  Tscale = RefLength / RefVelocity;
  
  
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
  for (int n=1; n<=NoBC; n++) {
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
  } // end of NoBC
	
  // 発熱密度の計算(有次元) -- 発熱量と発熱密度
  REAL_TYPE a, vol;
  a = deltaX*RefLength;
  vol = a*a*a;
  
  for (int n=1; n<=NoBC; n++) {
    
    if ( cmp[n].getType()==HEAT_SRC )
    {
      m = cmp[n].getMatOdr();
      
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
  for (int n=1; n<=NoBC; n++) {
    
    if ( cmp[n].getType()==DARCY )
    {
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
    
    if ( cmp[n].getType()==HEX )
    {
      m = cmp[n].getMatOdr();
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
  if ( Unit.Param == NONDIMENSIONAL ) {
    for (int n=0; n<NOFACE; n++) {
      
      switch ( BO[n].get_Class() )
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
      
      switch ( BO[n].get_Class() )
      {
        case OBC_OUTFLOW:
        case OBC_TRC_FREE:
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
