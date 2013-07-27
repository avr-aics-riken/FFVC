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
 * @file   History.C
 * @brief  FlowBase History class
 * @author kero
 */

#include "History.h"



// #################################################################
/**
 * @brief タイムスタンプの更新
 * @param [in] m_stp ステップ数
 * @param [in] m_tm  時刻
 * @param [in] vMax  速度最大値成分
 */
void History::updateTimeStamp(const int m_stp, const REAL_TYPE m_tm, const REAL_TYPE vMax)
{
  step  = m_stp;
  time  = m_tm;
  v_max = vMax;
}


// #################################################################
/**
 * @brief 反復過程の状況モニタのヘッダー出力
 * @param [in] fp 出力ファイルポインタ
 */
void History::printHistoryItrTitle(FILE* fp)
{
  fprintf(fp, "step=%16d  time=%13.6e    Itr_VP         Div_V    Itr_P          Norm (           i,            j,            k)\n", step, printTime());
}


// #################################################################
/**
 * @brief コンポーネントモニタの履歴出力
 * @param [in] fp  出力ファイルポインタ
 * @param [in] IC  反復管理クラス
 * @param [in] idx divの最大値の発生セルインデクス
 */
void History::printHistoryItr(FILE* fp, const ItrCtl* IC, const FB::Vec3i idx)
{
  const ItrCtl* ICd = &IC[ItrCtl::ic_div];     ///< 圧力-速度反復
  const ItrCtl* ICp = &IC[ItrCtl::ic_prs_pr];  ///< 圧力のPoisson反復
	fprintf(fp, "                                           %8d %13.6e %8d %13.6e (%12d, %12d, %12d)\n",
          ICd->LoopCount, ICd->get_normValue(),
          ICp->LoopCount, ICp->get_normValue(),
          idx.x, idx.y, idx.z);
}


// #################################################################
/**
 * @brief 反復過程の状況モニタのヘッダー出力
 * @param [in] fp 出力ファイルポインタ
 * @param [in] IC 反復管理クラス
 * @param [in] C  制御クラス
 * @param [in] disp  計算時間表示の有無
 */
void History::printHistoryTitle(FILE* fp, const ItrCtl* IC, const Control* C, const bool disp)
{
  const ItrCtl* ICp1 = &IC[ItrCtl::ic_prs_pr];  /// 圧力のPoisson反復
  const ItrCtl* ICv  = &IC[ItrCtl::ic_vis_cn];  /// 粘性項のCrank-Nicolson反復
  const ItrCtl* ICt  = &IC[ItrCtl::ic_tdf_ei];  /// 温度の拡散項の反復
  
  if ( Unit_Log == DIMENSIONAL )
  {
    fprintf(fp, "    step      time[sec]  v_max[m/s] ItrVP v_div_max[-]");
  }
  else
  {
    fprintf(fp, "    step        time[-]    v_max[-] ItrVP v_div_max[-]");
  }
  
  if ( C->KindOfSolver != SOLID_CONDUCTION )
  {
    switch (C->AlgorithmF)
    {
      case Control::Flow_FS_EE_EE:
      case Control::Flow_FS_AB2:
      case Control::Flow_FS_AB_CN:
        fprintf(fp, "  ItrP");
        if      (ICp1->get_normType() == ItrCtl::dx_b)       fprintf(fp, "        dx_b");
        else if (ICp1->get_normType() == ItrCtl::r_b)        fprintf(fp, "         r_b");
        else if (ICp1->get_normType() == ItrCtl::r_r0)       fprintf(fp, "        r_r0");
        break;
    }
    
    if (C->AlgorithmF == Control::Flow_FS_AB_CN)
    {
      fprintf(fp, "  ItrV");
      if      (ICv->get_normType() == ItrCtl::dx_b)       fprintf(fp, "        dx_b");
      else if (ICv->get_normType() == ItrCtl::r_b)        fprintf(fp, "         r_b");
      else if (ICv->get_normType() == ItrCtl::r_r0)       fprintf(fp, "        r_r0");
      else
      {
        printf("\n\tError : Norm selection type=%d\n", ICv->get_normType());
        Exit(0);
      }
    }
    
    if ( C->isHeatProblem() )
    {
      switch (C->AlgorithmH)
      {
        case Control::Heat_EE_EI:
          fprintf(fp, "  ItrT");
          if      (ICt->get_normType() == ItrCtl::dx_b)       fprintf(fp, "        dx_b");
          else if (ICt->get_normType() == ItrCtl::r_b)        fprintf(fp, "         r_b");
          else if (ICt->get_normType() == ItrCtl::r_r0)       fprintf(fp, "        r_r0");
          else
          {
            printf("\n\tError : Norm selection type=%d\n", ICt->get_normType());
            Exit(0);
          }
          break;
      }
    }
    
    fprintf(fp, "     deltaP       avrP     deltaV       avrV");
    if ( C->isHeatProblem() ) fprintf(fp, "     deltaT       avrT");
  }
  else
  {
    switch (C->AlgorithmH)
    {
      case Control::Heat_EE_EI:
        fprintf(fp, "  ItrT");
        if      (ICt->get_normType() == ItrCtl::dx_b)       fprintf(fp, "        dx_b");
        else if (ICt->get_normType() == ItrCtl::r_b)        fprintf(fp, "         r_b");
        else if (ICt->get_normType() == ItrCtl::r_r0)       fprintf(fp, "        r_r0");
        else
        {
          printf("\n\tError : Norm selection type=%d\n", ICt->get_normType());
          Exit(0);
        }
        break;
    }
    
    fprintf(fp, "     deltaT       avrT");
  }
  
  if ( disp )
  {
    fprintf(fp, "     time[sec]");
  }
  
  fprintf(fp, "\n");
}


// #################################################################
/**
 * @brief 標準履歴の出力
 * @param [in] fp    出力ファイルポインタ
 * @param [in] avr   1タイムステップの平均値　（0-pressure, 1-velocity, 2-temperature)
 * @param [in] rms   1タイムステップの変化量　（0-pressure, 1-velocity, 2-temperature)
 * @param [in] IC    ItrCtlクラスのポインタ
 * @param [in] C     Controlクラスへのポインタ
 * @param [in] stptm 1タイムステップの計算時間
 * @param [in] disp  計算時間表示の有無
 */
void History::printHistory(FILE* fp, const double* avr, const double* rms, const ItrCtl* IC, const Control* C, const double stptm, const bool disp)
{
  const ItrCtl* ICp1 = &IC[ItrCtl::ic_prs_pr];  ///< 圧力のPoisson反復
  const ItrCtl* ICv  = &IC[ItrCtl::ic_vis_cn];  ///< 粘性項のCrank-Nicolson反復
  const ItrCtl* ICt  = &IC[ItrCtl::ic_tdf_ei];  ///< 温度の拡散項の反復
  const ItrCtl* ICd  = &IC[ItrCtl::ic_div];     ///< 圧力-速度反復
  
  fprintf(fp, "%8d %14.6e %11.4e %5d  %11.4e",
          step, printTime(), printVmax(), ICd->LoopCount, ICd->get_normValue() );

  if ( C->KindOfSolver != SOLID_CONDUCTION ) 
  {
    switch (C->AlgorithmF)
    {
      case Control::Flow_FS_EE_EE:
      case Control::Flow_FS_AB2:
      case Control::Flow_FS_AB_CN:
        fprintf(fp, " %5d %11.4e", ICp1->LoopCount, ICp1->get_normValue());
        break;
    }
    
    if (C->AlgorithmF == Control::Flow_FS_AB_CN) 
    {
      fprintf(fp, " %5d %11.4e", ICv->LoopCount, ICv->get_normValue());
    }
    
    if ( C->isHeatProblem() ) 
    {
      switch (C->AlgorithmH) 
      {				
        case Control::Heat_EE_EI:
          fprintf(fp, " %5d %11.4e", ICt->LoopCount, ICt->get_normValue());
          break;
      }
    }
    
    fprintf(fp, " %10.3e %10.3e %10.3e %10.3e", rms[var_Pressure], avr[var_Pressure], rms[var_Velocity], avr[var_Velocity]);
    if ( C->isHeatProblem() ) fprintf(fp, " %10.3e %10.3e", rms[var_Temperature], avr[var_Temperature]);
  }
  else 
  {
    switch (C->AlgorithmH) 
    {				
      case Control::Heat_EE_EI:
        fprintf(fp, " %5d %11.4e", ICt->LoopCount, ICt->get_normValue());
        break;
    }
    fprintf(fp, " %10.3e %10.3e", rms[var_Temperature], avr[var_Temperature]);
  }
  
  if ( disp )
  {
    fprintf(fp, "%14.6e", stptm);
  }
  fprintf(fp, "\n");
  fflush(fp);
}



// #################################################################
/**
 * @brief コンポーネントモニタのヘッダー出力
 * @param [in] fp  出力ファイルポインタ
 * @param [in] cmp CompoListクラスのポインタ
 * @param [in] C   Controlクラスへのポインタ
 */
void History::printHistoryCompoTitle(FILE* fp, const CompoList* cmp, const Control* C)
{
  int cid;
  int def;

  if ( Unit_Log == DIMENSIONAL ) {
    fprintf(fp, "    step      time[sec]");
  }
  else {
    fprintf(fp, "    step        time[-]");
  }

  for (int i=1; i<=C->NoCompo; i++) {
    def = cmp[i].getDef();
    
    switch ( cmp[i].getType() ) {
      case SPEC_VEL:
        fprintf(fp, "  V[%03d:%03d]", i, def);
        break;
        
      case SPEC_VEL_WH:
        fprintf(fp, "  V[%03d:%03d]  Q[%03d:%03d] qa[%03d:%03d]", i, def, cid, def, cid, def);
        break;
        
      case OUTFLOW:
        fprintf(fp, "  V[%03d:%03d]", i, def);
        if ( C->isHeatProblem() ) {
          fprintf(fp, "  Q[%03d:%03d] qa[%03d:%03d]", i, def, i, def);
        }
        break;
        
      case HEX:
        fprintf(fp, "     Va[%03d]    DPa[%03d]", i, i);
        break;
        
      case DARCY:
        fprintf(fp, "      U[%03d]      V[%03d]      W[%03d]", i, i, i);
        break;
        
      case HEATFLUX:
      case TRANSFER:
      case ISOTHERMAL:
      case RADIANT:
        fprintf(fp, "  Q[%03d:%03d] qa[%03d:%03d]", i, def, i, def);
        break;
      
      case CELL_MONITOR:
        if ( cmp[i].isVarEncoded(var_Velocity) )     fprintf(fp, "      V[%03d]", i);
        if ( cmp[i].isVarEncoded(var_Pressure) )     fprintf(fp, "      P[%03d]", i);
        if ( cmp[i].isVarEncoded(var_Temperature) )  fprintf(fp, "      T[%03d]", i);
        if ( cmp[i].isVarEncoded(var_TotalP) )       fprintf(fp, "     TP[%03d]", i);
        break;
    }
  }

  fprintf(fp, "\n");
}


// #################################################################
/**
 * @brief コンポーネントモニタの履歴出力(dimensional value)
 * @param [in] fp  出力ファイルポインタ
 * @param [in] cmp CompoListクラスのポインタ
 * @param [in] C   Controlクラスへのポインタ
 */
void History::printHistoryCompo(FILE* fp, const CompoList* cmp, const Control* C)
{
  REAL_TYPE c, pp, dr, dp;
  const REAL_TYPE p0 = RefDensity * RefVelocity * RefVelocity;
  
  fprintf(fp, "%8d %14.6e", step, printTime());
  
  for (int i=1; i<=C->NoCompo; i++) {
    switch ( cmp[i].getType() ) {
      case SPEC_VEL:
        fprintf(fp, " %11.4e", printVel(cmp[i].val[var_Velocity]) );
        break;
        
      case SPEC_VEL_WH:
        pp = printQF(cmp[i].get_Mon_Calorie());
        fprintf(fp, " %11.4e %11.4e %11.4e", printVel(cmp[i].val[var_Velocity]), pp, pp/(REAL_TYPE)cmp[i].getElement() );
        break;
      
      case OUTFLOW:
        fprintf(fp, " %11.4e %11.4e", printVel(cmp[i].val[var_Velocity]) );
        if ( C->isHeatProblem() ) {
          pp = printQF(cmp[i].get_Mon_Calorie());
          fprintf(fp, " %11.4e %11.4e", pp, pp/(REAL_TYPE)cmp[i].getElement() ); // [W], [W/m^2]
        }
        break;
        
      case HEX:
        dr = cmp[i].ca[4]; // 熱交換器の無次元厚さ
        dp = cmp[i].val[var_Pressure] * p0 * dr/dh * RefLength;
        fprintf(fp, " %11.4e %11.4e", printVel(cmp[i].val[var_Velocity]), dp);
        break;
        
      case DARCY:
        fprintf(fp, " %11.4e %11.4e %11.4e", cmp[i].val[0], cmp[i].val[1], cmp[i].val[2]); // ? unit
        break;
        
      case HEATFLUX:
      case TRANSFER:
      case ISOTHERMAL:
      case RADIANT:
        pp = printQF(cmp[i].get_Mon_Calorie());
        fprintf(fp, " %11.4e %11.4e", pp, pp/cmp[i].area);
        break;
      
      case HEAT_SRC:
        fprintf(fp, " %11.4e", printQV(cmp[i].get_Mon_Calorie()));
        break;
        
      case CELL_MONITOR:
        if ( cmp[i].isVarEncoded(var_Velocity) )     fprintf(fp, " %11.4e", printVel(cmp[i].val[var_Velocity]) );
        if ( cmp[i].isVarEncoded(var_Pressure) )     fprintf(fp, " %11.4e", printPrs(cmp[i].val[var_Pressure]) );
        if ( cmp[i].isVarEncoded(var_Temperature) )  fprintf(fp, " %11.4e", printTmp(cmp[i].val[var_Temperature]) );
        if ( cmp[i].isVarEncoded(var_TotalP) )       fprintf(fp, " %11.4e", printTP(cmp[i].val[var_TotalP]) );
        break;
    }
  }

  fprintf(fp, "\n");
  fflush(fp);
}


// #################################################################
/**
 * @brief 計算領域の流束履歴のヘッダー出力
 * @param [in] fp 出力ファイルポインタ
 * @param [in] C  コントロールクラス
 */
void History::printHistoryDomfxTitle(FILE* fp, const Control* C)
{
  if ( Unit_Log == DIMENSIONAL ) {
    fprintf(fp, "    step      time[sec]");
  }
  else {
    fprintf(fp, "    step        time[-]");
  }
  
  for (int i=0; i<NOFACE; i++) fprintf(fp, "         Q:%s", FBUtility::getDirection(i).c_str());
  
  fprintf(fp, " >>      Balance : ");
  
  for (int i=0; i<NOFACE; i++) fprintf(fp, "         V:%s", FBUtility::getDirection(i).c_str());
  
  if (C->isHeatProblem()) {
    for (int i=0; i<NOFACE; i++) fprintf(fp, "         H:%s", FBUtility::getDirection(i).c_str());
  }
  fprintf(fp, "\n");
}


// #################################################################
/**
 * @brief 計算領域の流束履歴の出力
 * @param [in] fp 出力ファイルポインタ
 * @param [in] C  Controlクラスへのポインタ
 */
void History::printHistoryDomfx(FILE* fp, const Control* C)
{
  const REAL_TYPE sgn=-1.0;
  REAL_TYPE balance=0.0, s=1.0;
  
  for (int i=0; i<NOFACE; i++) {
    s *= sgn;
    balance += C->Q_Dface[i]*s;
  }
  
  fprintf(fp, "%8d %14.6e", step, printTime());
  
  for (int i=0; i<NOFACE; i++) fprintf(fp, " %12.4e", printMF(C->Q_Dface[i]) );
  fprintf(fp, " >> %12.4e : ", printMF(balance) );
  
  for (int i=0; i<NOFACE; i++) fprintf(fp, " %12.4e", printVel(C->V_Dface[i]) );
  
  if (C->isHeatProblem()) {
    for (int i=0; i<NOFACE; i++) fprintf(fp, " %12.4e", printQF(C->H_Dface[i]) ); // Watt
  }

  fprintf(fp, "\n");
  fflush(fp);
}


// #################################################################
/**
 * @brief 物体に働く力の履歴のヘッダー出力
 * @param [in] fp 出力ファイルポインタ
 */
void History::printHistoryForceTitle(FILE* fp)
{
  if ( Unit_Log == DIMENSIONAL ) {
    fprintf(fp, "    step      time[sec]");
  }
  else {
    fprintf(fp, "    step        time[-]");
  }
  
  if ( Unit_Log == DIMENSIONAL ) {
    fprintf(fp, "        Fx[N]         Fy[N]         Fz[N]");
  }
  else {
    fprintf(fp, "        Fx[-]         Fy[-]         Fz[-]");
  }

  fprintf(fp, "\n");
}


// #################################################################
/**
 * @brief 物体に働く力の履歴の出力
 * @param [in] fp 出力ファイルポインタ
 */
void History::printHistoryForce(FILE* fp, const REAL_TYPE* force)
{
  fprintf(fp, "%8d %14.6e", step, printTime());
  fprintf(fp, " %12.4e  %12.4e  %12.4e", printForce(force[0]), printForce(force[1]), printForce(force[2]) );
  
  fprintf(fp, "\n");
  fflush(fp);
}


// #################################################################
/**
 * @brief 反復過程の状況モニタのヘッダー出力
 * @param [in] fp 出力ファイルポインタ
 */
void History::printHistoryWallTitle(FILE* fp)
{
  if ( Unit_Log == DIMENSIONAL ) {
    fprintf(fp, "    step      time[sec]    Yp_Min[m]     Yp_Max[m]  Ut_Min[m/s]   Ut_Max[m/s]");
  }
  else {
    fprintf(fp, "    step        time[-]    Yp_Min[-]     Yp_Max[-]    Ut_Min[-]     Ut_Max[-]");
  }
  
  fprintf(fp, "\n");
}


// #################################################################
/**
 * @brief 壁面履歴の出力
 * @param [in] fp        出力ファイルポインタ
 * @param [in] range_Yp  壁座標の最小最大値
 * @param [in] range_Ut  摩擦速度の最小最大値
 */
void History::printHistoryWall(FILE* fp, const REAL_TYPE* range_Yp, const REAL_TYPE* range_Ut)
{
  fprintf(fp, "%8d %14.6e ", step, printTime() );
          
  fprintf(fp, "%12.6e %12.6e %12.6e %12.6e",
          printLen(range_Yp[0]), 
          printLen(range_Yp[1]),
          printVel(range_Yp[0]), 
          printVel(range_Yp[1]) );
  
  fprintf(fp, "\n");
  fflush(fp);
}
