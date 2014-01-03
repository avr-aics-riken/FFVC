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
 * @file   History.C
 * @brief  FlowBase History class
 * @author kero
 */

#include "History.h"


// #################################################################
// 履歴のCCNVファイルへの出力
void History::printCCNV(const double* avr, const double* rms, const IterationCtl* IC, const Control* C, const double stptm)
{
  const IterationCtl* ICp1 = &IC[ic_prs1];  ///< 圧力のPoisson反復
  const IterationCtl* ICv  = &IC[ic_vel1];  ///< 粘性項のCrank-Nicolson反復
  const IterationCtl* ICt  = &IC[ic_tmp1];  ///< 温度の拡散項の反復
  const IterationCtl* ICd  = &IC[ic_div];   ///< 圧力-速度反復
  
  
  // データ出力用
	double data[CCNV_MAX];
    
  // データ出力
  int c=0;
  data[c++] = (double)step;
  data[c++] = (double)printTime();
  
  
  if ( (C->KindOfSolver==FLOW_ONLY) ||
       (C->KindOfSolver==THERMAL_FLOW) ||
       (C->KindOfSolver==THERMAL_FLOW_NATURAL) ||
       (C->KindOfSolver==CONJUGATE_HEAT_TRANSFER) )
  {
    // Max Velocity
    data[c++] = (double)printVmax();
    
    // Iteration VP
    data[c++] = (double)ICd->getLoopCount()+1.0;
    
    // max divergence
    data[c++] = (double)ICd->getNormValue();
    
    switch (C->AlgorithmF)
    {
      case Flow_FS_EE_EE:
      case Flow_FS_AB2:
      case Flow_FS_AB_CN:
        // Iteration Pressure
        data[c++] = (double)ICp1->getLoopCount()+1.0;
        
        // Norm Pressure
        data[c++] = (double)ICp1->getNormValue();
        break;
    }
    
    if (C->AlgorithmF == Flow_FS_AB_CN)
    {
      // Iteration Velocity
      data[c++] = (double)ICv->getLoopCount()+1.0;
      
      // Norm Velocity
      data[c++] = (double)ICv->getNormValue();
    }
    
    // Delta Pressure[-];
    data[c++] = rms[var_Pressure];
    
    // Average Pressure[-]
    data[c++] = avr[var_Pressure];
    
    // Delta Velocity[-]
    data[c++] = rms[var_Velocity];
    
    // Average Velocity[-]
    data[c++] = avr[var_Velocity];
    
    
    if ( C->isHeatProblem() )
    {
      switch (C->AlgorithmH)
      {
        case Heat_EE_EE:
          break;
          
        case Heat_EE_EI:
          // Iteration Energy
          data[c++] = (double)ICt->getLoopCount()+1.0;
          
          // Norm Energy
          data[c++] = (double)ICt->getNormValue();
          break;
      }
      
      // Delta Energy[-]
      data[c++] = rms[var_Temperature];
      
      // Average Energy[-]
      data[c++] = avr[var_Temperature];
    }
  }
  else if (C->KindOfSolver==SOLID_CONDUCTION)
  {
    switch (C->AlgorithmH)
    {
      case Heat_EE_EE:
        break;
        
      case Heat_EE_EI:
        // Iteration Energy
        data[c++] = (double)ICt->getLoopCount()+1.0;
        
        // Norm Energy
        data[c++] = (double)ICt->getNormValue();
        break;
    }
    
    // Delta Energy[-]
    data[c++] = rms[var_Temperature];
    
    // Average Energy[-]
    data[c++] = avr[var_Temperature];
  }
  
  // time cost for 1-step
  data[c++] = stptm;
  
  if (c > CCNV_MAX) Exit(0);
  
  
  // データタイプ
	//  4: 単精度
	//  8: 倍精度
	long ccnv_pr = 8; // 出力はdoubleで固定 (long)C->Mode.Precision;
  
  
  // データを追加で書込みます
  FILE* fp = fopen(ccnvfile, "ab");
  fwrite(data, ccnv_pr, c, fp);
  fclose(fp);
}


// #################################################################
// CCNV履歴のヘッダ出力
void History::printCCNVtitle(const IterationCtl* IC, const Control* C)
{
  const IterationCtl* ICp1 = &IC[ic_prs1];  ///< 圧力のPoisson反復
  const IterationCtl* ICv  = &IC[ic_vel1];  ///< 粘性項のCrank-Nicolson反復
  const IterationCtl* ICt  = &IC[ic_tmp1];  ///< 温度の拡散項の反復
  
  // X軸タイトル数 最大2
	long num_x = 2;
	char x_title[2][128];
  for (int i=0; i<2; i++) memset(x_title[i], 0, sizeof(char)*128);
  
  std::string str("Step");
  sprintf(x_title[0], "%s", str.c_str());
  
  str = (Unit_Log == DIMENSIONAL) ? "Time[sec]" : "Time[-]";
  sprintf(x_title[1], "%s", str.c_str());
  
  
  // Y軸タイトル数 最大15
	char y_title[CCNV_MAX - 2][128];
  for (int i=0; i<CCNV_MAX-2; i++) memset(y_title[i], 0, sizeof(char)*128);
  int c=0;
  
  if ( (C->KindOfSolver==FLOW_ONLY) ||
       (C->KindOfSolver==THERMAL_FLOW) ||
       (C->KindOfSolver==THERMAL_FLOW_NATURAL) ||
       (C->KindOfSolver==CONJUGATE_HEAT_TRANSFER) )
  {
    str = (Unit_Log == DIMENSIONAL) ? "Maximum Velocity[m/s]" : "Maximum Velocity[-]";
    sprintf(y_title[c++], "%s", str.c_str()); // c=0
    
    str = "Iteration VP";
    sprintf(y_title[c++], "%s", str.c_str()); // c=1
    
    str = "Maximum Divergence[-]";
    sprintf(y_title[c++], "%s", str.c_str()); // c=2
    
    switch (C->AlgorithmF)
    {
      case Flow_FS_EE_EE:
      case Flow_FS_AB2:
      case Flow_FS_AB_CN:
        str = "Iteration Pressure";
        sprintf(y_title[c++], "%s", str.c_str());
        
        if      (ICp1->getNormType() == dx_b)  str = "dx_b";
        else if (ICp1->getNormType() == r_b)   str = "r_b";
        else if (ICp1->getNormType() == r_r0)  str = "r_r0";
        sprintf(y_title[c++], "%s", str.c_str());
        break;
    }
    
    if (C->AlgorithmF == Flow_FS_AB_CN)
    {
      str = "Iteration Velocity";
      sprintf(y_title[c++], "%s", str.c_str());
      
      if      (ICv->getNormType() == dx_b)  str = "dx_b";
      else if (ICv->getNormType() == r_b)   str = "r_b";
      else if (ICv->getNormType() == r_r0)  str = "r_r0";
      sprintf(y_title[c++], "%s", str.c_str());
    }
    
    str = "Delta Pressure[-]";
    sprintf(y_title[c++], "%s", str.c_str());
    
    str = "Average Pressure[-]";
    sprintf(y_title[c++], "%s", str.c_str());
    
    str = "Delta Velocity[-]";
    sprintf(y_title[c++], "%s", str.c_str());
    
    str = "Average Velocity[-]";
    sprintf(y_title[c++], "%s", str.c_str());
    
    if ( C->isHeatProblem() )
    {
      switch (C->AlgorithmH)
      {
        case Heat_EE_EE:
          break;
          
        case Heat_EE_EI:
          str = "Iteration Energy";
          sprintf(y_title[c++], "%s", str.c_str());
          
          if      (ICt->getNormType() == dx_b)   str = "dx_b";
          else if (ICt->getNormType() == r_b)    str = "r_b";
          else if (ICt->getNormType() == r_r0)   str = "r_r0";
          sprintf(y_title[c++], "%s", str.c_str());
          break;
      }
      
      str = "Delta Energy[-]";
      sprintf(y_title[c++], "%s", str.c_str());
      
      str = "Average Energy[-]";
      sprintf(y_title[c++], "%s", str.c_str());
    }
  }
  else if (C->KindOfSolver==SOLID_CONDUCTION)
  {
    switch (C->AlgorithmH)
    {
      case Heat_EE_EE:
        break;
        
      case Heat_EE_EI:
        str = "Iteration Energy";
        sprintf(y_title[c++], "%s", str.c_str());
        
        if      (ICt->getNormType() == dx_b)   str = "dx_b";
        else if (ICt->getNormType() == r_b)    str = "r_b";
        else if (ICt->getNormType() == r_r0)   str = "r_r0";
        sprintf(y_title[c++], "%s", str.c_str());
        break;
    }
    
    str = "Delta Energy[-]";
    sprintf(y_title[c++], "%s", str.c_str());
    
    str = "Average Energy[-]";
    sprintf(y_title[c++], "%s", str.c_str());
  }
  
  sprintf(y_title[c++], "time[s]/Step");
  
  if (c > CCNV_MAX) Exit(0);
  
  
  
  // CCNV ID (固定値)
	long ccnv_id = 64408111;
  
	// CCNV VERSION (固定値)
	char ccnv_ver[16] = {"1.0.0"};
  
	// Balnk
	char ccnv_actid[128] = {""};
	char ccnv_jobid[128] = {""};
	char ccnv_pass[128] = {""};
  
	// SOLVER ID for FFV-C provided by VINAS
	char ccnv_solid[128] = {"OS-0081-0003-000000-0000-0000-0000"};
  
	// データタイプ
	//  4: 単精度
	//  8: 倍精度
	long ccnv_pr = 8; // (long)C->Mode.Precision;
  
	// ファイル名（固定）
	char ccnvfile[16] = {"ccnv.log"};
  
  
  // ヘッダーの出力
	FILE* fp = fopen(ccnvfile, "wb");
  
	fwrite(&ccnv_id, 4, 1, fp);
	fwrite(ccnv_ver, 16, 1, fp);
	fwrite(ccnv_actid, 128, 1, fp);
	fwrite(ccnv_jobid, 128, 1, fp);
	fwrite(ccnv_pass, 128, 1, fp);
	fwrite(ccnv_solid, 128, 1, fp);
  
	fwrite(&num_x, 4, 1, fp);
	fwrite(x_title[0], 128, 1, fp);
	fwrite(x_title[1], 128, 1, fp);
  
	fwrite(&c, 4, 1, fp);
  
  // Y軸個数分
  for (int i=0; i<c; i++)
  {
    fwrite(y_title[i], 128, 1, fp);
  }
  
	fwrite(&ccnv_pr, 4, 1, fp);
	
	fclose(fp);
}


// #################################################################
// 標準履歴の出力
void History::printHistory(FILE* fp, const double* avr, const double* rms, const IterationCtl* IC, const Control* C, const double stptm, const bool disp)
{
  const IterationCtl* ICp1 = &IC[ic_prs1];  ///< 圧力のPoisson反復
  const IterationCtl* ICv  = &IC[ic_vel1];  ///< 粘性項のCrank-Nicolson反復
  const IterationCtl* ICt  = &IC[ic_tmp1];  ///< 温度の拡散項の反復
  const IterationCtl* ICd  = &IC[ic_div];   ///< 圧力-速度反復
  
  fprintf(fp, "%8d %14.6e", step, printTime());
  
  if ( (C->KindOfSolver==FLOW_ONLY) ||
      ( C->KindOfSolver==THERMAL_FLOW) ||
      ( C->KindOfSolver==THERMAL_FLOW_NATURAL) ||
      ( C->KindOfSolver==CONJUGATE_HEAT_TRANSFER) )
  {
    fprintf(fp, " %11.4e %5d  %11.4e",
            printVmax(), ICd->getLoopCount(), ICd->getNormValue() );
    
    switch (C->AlgorithmF)
    {
      case Flow_FS_EE_EE:
      case Flow_FS_AB2:
      case Flow_FS_AB_CN:
        fprintf(fp, " %5d %11.4e", ICp1->getLoopCount()+1, ICp1->getNormValue());
        break;
    }
    
    if (C->AlgorithmF == Flow_FS_AB_CN)
    {
      fprintf(fp, " %5d %11.4e", ICv->getLoopCount()+1, ICv->getNormValue());
    }
    
    fprintf(fp, " %10.3e %10.3e %10.3e %10.3e",
            rms[var_Pressure],
            avr[var_Pressure],
            rms[var_Velocity],
            avr[var_Velocity]);
    
    if ( C->isHeatProblem() )
    {
      switch (C->AlgorithmH)
      {
        case Heat_EE_EE:
          break;
          
        case Heat_EE_EI:
          fprintf(fp, " %5d %11.4e", ICt->getLoopCount()+1, ICt->getNormValue());
          break;
      }
      
      fprintf(fp, " %10.3e %10.3e", rms[var_Temperature], avr[var_Temperature]);
    }
  }
  else if (C->KindOfSolver==SOLID_CONDUCTION)
  {
    switch (C->AlgorithmH)
    {
      case Heat_EE_EE:
        break;
        
      case Heat_EE_EI:
        fprintf(fp, " %5d %11.4e", ICt->getLoopCount()+1, ICt->getNormValue());
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
// 標準履歴モニタのヘッダー出力
void History::printHistoryTitle(FILE* fp, const IterationCtl* IC, const Control* C, const bool disp)
{
  const IterationCtl* ICp1 = &IC[ic_prs1];  /// 圧力のPoisson反復
  const IterationCtl* ICv  = &IC[ic_vel1];  /// 粘性項のCrank-Nicolson反復
  const IterationCtl* ICt  = &IC[ic_tmp1];  /// 温度の拡散項の反復
  
  if ( Unit_Log == DIMENSIONAL )
  {
    fprintf(fp, "    step      time[sec]");
  }
  else
  {
    fprintf(fp, "    step        time[-]");
  }
  
  if ( (C->KindOfSolver == FLOW_ONLY) ||
      ( C->KindOfSolver == THERMAL_FLOW) ||
      ( C->KindOfSolver == THERMAL_FLOW_NATURAL) )
  {
    if ( Unit_Log == DIMENSIONAL )
    {
      fprintf(fp, "  v_max[m/s] ItrVP v_div_max[-]");
    }
    else
    {
      fprintf(fp, "    v_max[-] ItrVP v_div_max[-]");
    }
    
    switch (C->AlgorithmF)
    {
      case Flow_FS_EE_EE:
      case Flow_FS_AB2:
      case Flow_FS_AB_CN:
        fprintf(fp, "  ItrP");
        if      (ICp1->getNormType() == dx_b)       fprintf(fp, "        dx_b");
        else if (ICp1->getNormType() == r_b)        fprintf(fp, "         r_b");
        else if (ICp1->getNormType() == r_r0)       fprintf(fp, "        r_r0");
        break;
    }
    
    if (C->AlgorithmF == Flow_FS_AB_CN)
    {
      fprintf(fp, "  ItrV");
      if      (ICv->getNormType() == dx_b)       fprintf(fp, "        dx_b");
      else if (ICv->getNormType() == r_b)        fprintf(fp, "         r_b");
      else if (ICv->getNormType() == r_r0)       fprintf(fp, "        r_r0");
      else
      {
        printf("\n\tError : Norm selection type=%d\n", ICv->getNormType());
        Exit(0);
      }
    }
    
    fprintf(fp, "     deltaP       avrP     deltaV       avrV");
    
    if ( C->isHeatProblem() )
    {
      switch (C->AlgorithmH)
      {
        case Heat_EE_EE:
          break;
          
        case Heat_EE_EI:
          fprintf(fp, "  ItrT");
          if      (ICt->getNormType() == dx_b)       fprintf(fp, "        dx_b");
          else if (ICt->getNormType() == r_b)        fprintf(fp, "         r_b");
          else if (ICt->getNormType() == r_r0)       fprintf(fp, "        r_r0");
          else
          {
            printf("\n\tError : Norm selection type=%d\n", ICt->getNormType());
            Exit(0);
          }
          break;
      }
      
      fprintf(fp, "     deltaT       avrT");
    }
  }
  else if ( C->KindOfSolver == SOLID_CONDUCTION )
  {
    switch (C->AlgorithmH)
    {
      case Heat_EE_EE:
        break;
        
      case Heat_EE_EI:
        fprintf(fp, "  ItrT");
        if      (ICt->getNormType() == dx_b)       fprintf(fp, "        dx_b");
        else if (ICt->getNormType() == r_b)        fprintf(fp, "         r_b");
        else if (ICt->getNormType() == r_r0)       fprintf(fp, "        r_r0");
        else
        {
          printf("\n\tError : Norm selection type=%d\n", ICt->getNormType());
          Exit(0);
        }
        break;
    }
    
    fprintf(fp, "     deltaT       avrT");
  }
  else // CONJUGATE_HEAT_TRANSFER
  {
    
  }
  
  if ( disp )
  {
    fprintf(fp, "     time[sec]");
  }
  
  fprintf(fp, "\n");
}


// #################################################################
// コンポーネントモニタの履歴出力
void History::printHistoryCompo(FILE* fp, const CompoList* cmp, const Control* C, const REAL_TYPE dt)
{
  REAL_TYPE dr, dp;
  const REAL_TYPE p0 = RefDensity * RefVelocity * RefVelocity;
  
  fprintf(fp, "%8d %14.6e", step, printTime());
  
  for (int i=1; i<=C->NoCompo; i++)
  {
    switch ( cmp[i].getType() )
    {
      case SPEC_VEL:
        if ( !cmp[i].isHeatMode() )
        {
          fprintf(fp, " %11.4e", printVel(cmp[i].val[var_Velocity]) );
        }
        else
        {
          fprintf(fp, " %11.4e %11.4e", printVel(cmp[i].val[var_Velocity]), printQF(cmp[i].getMonCalorie()) );
        }
        break;
      
      case OUTFLOW:
        fprintf(fp, " %11.4e", printVel(cmp[i].val[var_Velocity]) );
        if ( C->isHeatProblem() )
        {
          fprintf(fp, " %11.4e", printQF(cmp[i].getMonCalorie()) ); // [W]
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
        fprintf(fp, " %11.4e", printQF(cmp[i].getMonCalorie()) );
        break;
      
      case HEAT_SRC:
        fprintf(fp, " %11.4e", printQV(cmp[i].getMonCalorie()));
        break;
    }
  }

  fprintf(fp, "\n");
  fflush(fp);
}


// #################################################################
// コンポーネントモニタのヘッダー出力
void History::printHistoryCompoTitle(FILE* fp, const CompoList* cmp, const Control* C)
{
  if ( Unit_Log == DIMENSIONAL )
  {
    fprintf(fp, "    step      time[sec]");
  }
  else
  {
    fprintf(fp, "    step        time[-]");
  }
  
  for (int i=1; i<=C->NoCompo; i++)
  {
    
    switch ( cmp[i].getType() )
    {
      case SPEC_VEL:
        if ( !cmp[i].isHeatMode() ) fprintf(fp, "       V[%02d]", i);
        else                        fprintf(fp, "       V[%02d]       Q[%02d]", i, i);
        break;
        
      case OUTFLOW:
        fprintf(fp, "       V[%02d]", i);
        if ( C->isHeatProblem() )
        {
          fprintf(fp, "       Q[%02d]", i);
        }
        break;
        
      case HEX:
        fprintf(fp, "      Va[%02d]     DPa[%02d]", i, i);
        break;
        
      case DARCY:
        fprintf(fp, "      U[%02d]      V[%02d]      W[%02d]", i, i, i);
        break;
        
      case HEATFLUX:
      case TRANSFER:
      case ISOTHERMAL:
      case RADIANT:
        fprintf(fp, "       Q[%02d]", i);
        break;
    }
  }
  
  fprintf(fp, "\n");
}


// #################################################################
// 計算領域の流束履歴の出力
void History::printHistoryDomfx(FILE* fp, const Control* C, const REAL_TYPE dt)
{
  const REAL_TYPE sgn=-1.0;
  REAL_TYPE balance=0.0, s=1.0;
  
  for (int i=0; i<NOFACE; i++)
  {
    s *= sgn;
    balance += C->Q_Dface[i]*s;
  }
  
  fprintf(fp, "%8d %14.6e", step, printTime());
  
  for (int i=0; i<NOFACE; i++) fprintf(fp, " %12.4e", printMF(C->Q_Dface[i]) );
  fprintf(fp, " >> %12.4e : ", printMF(balance) );
  
  for (int i=0; i<NOFACE; i++) fprintf(fp, " %12.4e", printVel(C->V_Dface[i]) );
  
  if (C->isHeatProblem())
  {
    for (int i=0; i<NOFACE; i++) fprintf(fp, " %12.4e", printQF(C->H_Dface[i]) ); // [W]
  }

  fprintf(fp, "\n");
  fflush(fp);
}


// #################################################################
// 計算領域の流束履歴のヘッダー出力
void History::printHistoryDomfxTitle(FILE* fp, const Control* C)
{
  if ( Unit_Log == DIMENSIONAL )
  {
    fprintf(fp, "    step      time[sec]");
  }
  else
  {
    fprintf(fp, "    step        time[-]");
  }
  
  for (int i=0; i<NOFACE; i++) fprintf(fp, "         Q:%s", FBUtility::getDirection(i).c_str());
  
  fprintf(fp, " >>      Balance : ");
  
  for (int i=0; i<NOFACE; i++) fprintf(fp, "         V:%s", FBUtility::getDirection(i).c_str());
  
  if (C->isHeatProblem())
  {
    for (int i=0; i<NOFACE; i++) fprintf(fp, "         H:%s", FBUtility::getDirection(i).c_str());
  }
  fprintf(fp, "\n");
}


// #################################################################
// 物体に働く力の履歴の出力
void History::printHistoryForce(FILE* fp, const REAL_TYPE* force)
{
  fprintf(fp, "%8d %14.6e", step, printTime());
  fprintf(fp, " %12.4e  %12.4e  %12.4e", printForce(force[0]), printForce(force[1]), printForce(force[2]) );
  
  fprintf(fp, "\n");
  fflush(fp);
}


// #################################################################
// 物体に働く力の履歴のヘッダー出力
void History::printHistoryForceTitle(FILE* fp)
{
  if ( Unit_Log == DIMENSIONAL )
  {
    fprintf(fp, "    step      time[sec]");
  }
  else
  {
    fprintf(fp, "    step        time[-]");
  }
  
  if ( Unit_Log == DIMENSIONAL )
  {
    fprintf(fp, "        Fx[N]         Fy[N]         Fz[N]");
  }
  else
  {
    fprintf(fp, "        Fx[-]         Fy[-]         Fz[-]");
  }

  fprintf(fp, "\n");
}


// #################################################################
// 反復履歴出力
void History::printHistoryItr(FILE* fp, const IterationCtl* IC, const FB::Vec3i idx)
{
  const IterationCtl* ICd = &IC[ic_div];   ///< 圧力-速度反復
  const IterationCtl* ICp = &IC[ic_prs1];  ///< 圧力のPoisson反復
	fprintf(fp, "                                           %8d %13.6e %8d %13.6e (%12d, %12d, %12d)\n",
          ICd->getLoopCount(), ICd->getNormValue(),
          ICp->getLoopCount()+1, ICp->getNormValue(),
          idx.x, idx.y, idx.z);
}


// #################################################################
// 反復過程の状況モニタのヘッダー出力
void History::printHistoryItrTitle(FILE* fp)
{
  fprintf(fp, "step=%16d  time=%13.6e    Itr_VP         Div_V    Itr_P          Norm (           i,            j,            k)\n", step, printTime());
}


// #################################################################
// 壁面履歴の出力
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


// #################################################################
// 壁面履歴のヘッダー出力
void History::printHistoryWallTitle(FILE* fp)
{
  if ( Unit_Log == DIMENSIONAL )
  {
    fprintf(fp, "    step      time[sec]    Yp_Min[m]     Yp_Max[m]  Ut_Min[m/s]   Ut_Max[m/s]");
  }
  else
  {
    fprintf(fp, "    step        time[-]    Yp_Min[-]     Yp_Max[-]    Ut_Min[-]     Ut_Max[-]");
  }
  
  fprintf(fp, "\n");
}


// #################################################################
// タイムスタンプの更新
void History::updateTimeStamp(const int m_stp, const REAL_TYPE m_tm, const REAL_TYPE vMax)
{
  step  = m_stp;
  time  = m_tm;
  v_max = vMax;
}
