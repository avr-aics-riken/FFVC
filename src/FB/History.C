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
 * @file   History.C
 * @brief  FlowBase History class
 * @author aics
 */

#include "History.h"


// #################################################################
// 履歴のCCNVファイルへの出力
void History::printCCNV(const double* rms, const double* avr, const double* container, const Control* C, const double divergence, const double stptm)
{
  // データ出力用
	double data[CCNV_MAX];
    
  // データ出力
  int c=0;
  data[c++] = (double)step;
  data[c++] = (double)printTime();
  
  
  if ( (C->KindOfSolver==FLOW_ONLY) ||
       (C->KindOfSolver==THERMAL_FLOW) ||
       (C->KindOfSolver==THERMAL_FLOW_NATURAL) ||
       (C->KindOfSolver==CONJUGATE_HT) ||
       (C->KindOfSolver==CONJUGATE_HT_NATURAL) )
  {
    // Max Velocity
    data[c++] = (double)printVmax();

    
    switch (C->AlgorithmF)
    {
      case Flow_FS_EE_EE:
      case Flow_FS_AB2:
      case Flow_FS_AB_CN:
        // Iteration Pressure
        data[c++] = container[3*ic_prs1+0]; //getLoopCount();
        
        // Residual Pressure
        data[c++] = container[3*ic_prs1+1]; //getResidual();
        
        // Error Pressure
        data[c++] = container[3*ic_prs1+2]; //getError();
        break;
    }
    
    if (C->AlgorithmF == Flow_FS_AB_CN)
    {
      // Iteration Velocity
      data[c++] = container[3*ic_vel1+0]; //getLoopCount();
      
      // Residual Velocity
      data[c++] = container[3*ic_vel1+1]; //getResidual();
      
      // Error Velocity
      data[c++] = container[3*ic_vel1+2]; //getError();
    }
    
    // max divergence
    data[c++] = divergence;
    
    // Delta Pressure[-];
    data[c++] = rms[var_Pressure];
    
    // Average Pressure[-]
    data[c++] = avr[var_Pressure];
    
    // Delta Velocity[-]
    data[c++] = rms[var_Velocity];
    
    
    if ( C->isHeatProblem() )
    {
      switch (C->AlgorithmH)
      {
        case Heat_EE_EE:
          break;
          
        case Heat_EE_EI:
          // Iteration Energy
          data[c++] = container[3*ic_tmp1+0]; //getLoopCount();
          
          // Residual Energy
          data[c++] = container[3*ic_tmp1+1]; //getResidual();
          
          // Error Energy
          data[c++] = container[3*ic_tmp1+2]; //getError();
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
        data[c++] = container[3*ic_tmp1+0]; //getLoopCount();
        
        // Residual Energy
        data[c++] = container[3*ic_tmp1+1]; //getResidual();
        
        // Error Energy
        data[c++] = container[3*ic_tmp1+2]; //getError();
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
void History::printCCNVtitle(const int* container, const Control* C)
{
  
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
       (C->KindOfSolver==CONJUGATE_HT) ||
       (C->KindOfSolver==CONJUGATE_HT_NATURAL) )
  {
    str = (Unit_Log == DIMENSIONAL) ? "Maximum Velocity[m/s]" : "Maximum Velocity[-]";
    sprintf(y_title[c++], "%s", str.c_str()); // c=0

    
    switch (C->AlgorithmF)
    {
      case Flow_FS_EE_EE:
      case Flow_FS_AB2:
      case Flow_FS_AB_CN:
        str = "Iteration Pressure";
        sprintf(y_title[c++], "%s", str.c_str());
        
        if      (container[2*ic_prs1+0] == nrm_r_b)     str = "r_b";
        else if (container[2*ic_prs1+0] == nrm_r_x)     str = "r_x";
        else if (container[2*ic_prs1+0] == nrm_r_r0)    str = "r_r0";
        sprintf(y_title[c++], "%s", str.c_str());
        
        if      (container[2*ic_prs1+1] == nrm_dx)      str = "deltaP";
        else if (container[2*ic_prs1+1] == nrm_dx_x)    str = "deltaP_P";
        sprintf(y_title[c++], "%s", str.c_str());
        break;
    }
    
    if (C->AlgorithmF == Flow_FS_AB_CN)
    {
      str = "Iteration Velocity";
      sprintf(y_title[c++], "%s", str.c_str());
      
      if      (container[2*ic_vel1+0] == nrm_r_b)     str = "r_b";
      else if (container[2*ic_vel1+0] == nrm_r_x)     str = "r_x";
      else if (container[2*ic_vel1+0] == nrm_r_r0)    str = "r_r0";
      sprintf(y_title[c++], "%s", str.c_str());
      
      if      (container[2*ic_vel1+1] == nrm_dx)      str = "deltaV";
      else if (container[2*ic_vel1+1] == nrm_dx_x)    str = "deltaV_V";
      sprintf(y_title[c++], "%s", str.c_str());
    }
    
    str = "Maximum Divergence[-]";
    sprintf(y_title[c++], "%s", str.c_str()); // c=2
    
    str = "RMS Pressure[-]";
    sprintf(y_title[c++], "%s", str.c_str());
    
    str = "Average Pressure[-]";
    sprintf(y_title[c++], "%s", str.c_str());
    
    str = "Delta Velocity[-]";
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
          
          if      (container[2*ic_tmp1+0] == nrm_r_b)     str = "r_b";
          else if (container[2*ic_tmp1+0] == nrm_r_x)     str = "r_x";
          else if (container[2*ic_tmp1+0] == nrm_r_r0)    str = "r_r0";
          sprintf(y_title[c++], "%s", str.c_str());
          
          if      (container[2*ic_tmp1+1] == nrm_dx)      str = "deltaE";
          else if (container[2*ic_tmp1+1] == nrm_dx_x)    str = "deltaE_E";
          sprintf(y_title[c++], "%s", str.c_str());
          break;
      }
      
      str = "RMS Energy[-]";
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
        
        if      (container[2*ic_tmp1+0] == nrm_r_b)     str = "r_b";
        else if (container[2*ic_tmp1+0] == nrm_r_x)     str = "r_x";
        else if (container[2*ic_tmp1+0] == nrm_r_r0)    str = "r_r0";
        sprintf(y_title[c++], "%s", str.c_str());
        
        if      (container[2*ic_tmp1+1] == nrm_dx)      str = "deltaE";
        else if (container[2*ic_tmp1+1] == nrm_dx_x)    str = "deltaE_E";
        sprintf(y_title[c++], "%s", str.c_str());
        break;
    }
    
    str = "RMS Energy[-]";
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
// 物体に働く力の平均値出力（コンポーネント毎）
bool History::printForceAvr(const CompoList* cmp, const REAL_TYPE* frc)
{
  for (int n=1; n<=NoCompo; n++)
  {
    if ( cmp[n].getType()==OBSTACLE || cmp[n].getType()==SOLIDREV )
    {
      char fname[128];
      sprintf( fname, "force_avr_%s.txt", cmp[n].alias.c_str() );
      
      FILE* fp;
      
      if ( !(fp=fopen(fname, "w")) )
      {
        stamped_printf("\tSorry, can't open '%s' file.\n", fname);
        return false;
      }
      
      if ( Unit_Log == DIMENSIONAL )
      {
        fprintf(fp, " AvrStep   AvrTime[sec]          F_x           F_y           F_z\n");
      }
      else
      {
        fprintf(fp, " AvrStep     AvrTime[-]          F_x           F_y           F_z\n");
      }
      
      fprintf(fp, "%8d %14.6e %12.4e  %12.4e  %12.4e\n",
              step,
              printTime(),
              printForce(frc[3*n+0]),
              printForce(frc[3*n+1]),
              printForce(frc[3*n+2])
              );
      
      fclose(fp);
    }
  }
  return true;
}


// #################################################################
// 標準履歴の出力
void History::printHistory(FILE* fp,
                           const double* rms,
                           const double* avr,
                           const double* container,
                           const Control* C,
                           const DivConvergence* DC,
                           const double stptm,
                           const bool disp)
{
  fprintf(fp, "%8d %14.6e", step, printTime());
  
  if ( (C->KindOfSolver==FLOW_ONLY) ||
      ( C->KindOfSolver==THERMAL_FLOW) ||
      ( C->KindOfSolver==THERMAL_FLOW_NATURAL) ||
      ( C->KindOfSolver==CONJUGATE_HT) ||
      ( C->KindOfSolver==CONJUGATE_HT_NATURAL) )
  {
    fprintf(fp, " %11.4e %5d    %12.5e", printVmax(), DC->Iteration, DC->divergence);
    
    switch (C->AlgorithmF)
    {
      case Flow_FS_EE_EE:
      case Flow_FS_AB2:
      case Flow_FS_AB_CN:
        fprintf(fp, " %5d %11.4e %11.4e", (int)container[3*ic_prs1+0], container[3*ic_prs1+1], container[3*ic_prs1+2]);
        break;
    }
    
    if (C->AlgorithmF == Flow_FS_AB_CN)
    {
      fprintf(fp, " %5d %11.4e %11.4e", (int)container[3*ic_vel1+0], container[3*ic_vel1+1], container[3*ic_vel1+2]);
    }
    
    fprintf(fp, "  %10.3e %10.3e %10.3e",
            rms[var_Pressure],
            avr[var_Pressure],
            rms[var_Velocity]);
    
    if ( C->isHeatProblem() )
    {
      switch (C->AlgorithmH)
      {
        case Heat_EE_EE:
          break;
          
        case Heat_EE_EI:
          fprintf(fp, " %5d %11.4e %11.4e", (int)container[3*ic_tmp1+0], container[3*ic_tmp1+1], container[3*ic_tmp1+2]);
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
        fprintf(fp, " %5d %11.4e %11.4e", (int)container[3*ic_tmp1+0], container[3*ic_tmp1+1], container[3*ic_tmp1+2]);
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
void History::printHistoryTitle(FILE* fp, const int* container, const Control* C, const DivConvergence* DC, const bool disp)
{
  fprintf(fp, "Column_Data_00\n");
  
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
      ( C->KindOfSolver == THERMAL_FLOW_NATURAL) ||
      ( C->KindOfSolver == CONJUGATE_HT) ||
      ( C->KindOfSolver == CONJUGATE_HT_NATURAL) )
  {
    if ( Unit_Log == DIMENSIONAL )
    {
      fprintf(fp, "  v_max[m/s]");
    }
    else
    {
      fprintf(fp, "    v_max[-]");
    }
    
    fprintf(fp, "  ItrD");
    if ( DC->divType == nrm_div_max )
    {
      fprintf(fp, " Divergence[max]");
    }
    else
    {
      fprintf(fp, "  Divergence[L2]");
    }
    
    switch (C->AlgorithmF)
    {
      case Flow_FS_EE_EE:
      case Flow_FS_AB2:
      case Flow_FS_AB_CN:
        fprintf(fp, "  ItrP");
        if      (container[2*ic_prs1+0] == nrm_r_b)     fprintf(fp, "         r_b");
        else if (container[2*ic_prs1+0] == nrm_r_x)     fprintf(fp, "         r_x");
        else if (container[2*ic_prs1+0] == nrm_r_r0)    fprintf(fp, "        r_r0");
        
        if      (container[2*ic_prs1+1] == nrm_dx)      fprintf(fp, "      deltaP");
        else if (container[2*ic_prs1+1] == nrm_dx_x)    fprintf(fp, "    deltaP_P");
        break;
    }
    
    if (C->AlgorithmF == Flow_FS_AB_CN)
    {
      fprintf(fp, "  ItrV");
      if      (container[2*ic_vel1+0] == nrm_r_b)     fprintf(fp, "         r_b");
      else if (container[2*ic_vel1+0] == nrm_r_x)     fprintf(fp, "         r_x");
      else if (container[2*ic_vel1+0] == nrm_r_r0)    fprintf(fp, "        r_r0");

      if      (container[2*ic_vel1+1] == nrm_dx)      fprintf(fp, "      deltaV");
      else if (container[2*ic_vel1+1] == nrm_dx_x)    fprintf(fp, "    deltaV_V");
    }
    
    fprintf(fp, "        rmsP       avrP       rmsV");
    
    if ( C->isHeatProblem() )
    {
      switch (C->AlgorithmH)
      {
        case Heat_EE_EE:
          break;
          
        case Heat_EE_EI:
          fprintf(fp, "  ItrE");
          if      (container[2*ic_tmp1+0] == nrm_r_b)     fprintf(fp, "         r_b");
          else if (container[2*ic_tmp1+0] == nrm_r_x)     fprintf(fp, "         r_x");
          else if (container[2*ic_tmp1+0] == nrm_r_r0)    fprintf(fp, "        r_r0");
          
          if      (container[2*ic_tmp1+1] == nrm_dx)      fprintf(fp, "      deltaT");
          else if (container[2*ic_tmp1+1] == nrm_dx_x)    fprintf(fp, "    deltaT_T");
          break;
      }
      
      fprintf(fp, "       rmsE       avrE");
    }
  }
  else if ( C->KindOfSolver == SOLID_CONDUCTION )
  {
    switch (C->AlgorithmH)
    {
      case Heat_EE_EE:
        break;
        
      case Heat_EE_EI:
        fprintf(fp, "  ItrE");
        if      (container[2*ic_tmp1+0] == nrm_r_b)     fprintf(fp, "         r_b");
        else if (container[2*ic_tmp1+0] == nrm_r_x)     fprintf(fp, "         r_x");
        else if (container[2*ic_tmp1+0] == nrm_r_r0)    fprintf(fp, "        r_r0");

        if      (container[2*ic_tmp1+1] == nrm_dx)      fprintf(fp, "      deltaT");
        else if (container[2*ic_tmp1+1] == nrm_dx_x)    fprintf(fp, "    deltaT_T");
        break;
    }
    
    fprintf(fp, "       rmsE       avrE");
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
        dp = cmp[i].val[var_Pressure] * p0 * dr * RefLength;
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
  fprintf(fp, "Column_Data_00\n");
  
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
      case HEAT_SRC:
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
  REAL_TYPE balance=0.0, s=-1.0;
  
  fprintf(fp, "%8d %14.6e", step, printTime());
  
  for (int i=0; i<NOFACE; i++) {
    s *= sgn;
    fprintf(fp, " %12.4e", printMF(C->Q_Dface[i]*s) );
    balance += C->Q_Dface[i]*s;
  }
  fprintf(fp, " >   %12.4e : ", printMF(balance) );
  
  
  //for (int i=0; i<NOFACE; i++) fprintf(fp, " %12.4e", printVel(C->V_Dface[i]) );
  
  
  if (C->isHeatProblem())
  {
    for (int i=0; i<NOFACE; i++) fprintf(fp, " %12.4e", printQF(C->H_Dface[i]) ); // [W]
    
    balance = 0.0;
    for (int i=0; i<NOFACE; i++) {
      balance += C->H_Dface[i];
    }
    fprintf(fp, " >   %12.4e", printQF(balance) );
  }

  fprintf(fp, "\n");
  fflush(fp);
}


// #################################################################
// 計算領域の流束履歴のヘッダー出力
void History::printHistoryDomfxTitle(FILE* fp, const Control* C)
{
  fprintf(fp, "Column_Data_00\n");
  
  if ( Unit_Log == DIMENSIONAL )
  {
    fprintf(fp, "    step      time[sec]");
  }
  else
  {
    fprintf(fp, "    step        time[-]");
  }
  
  for (int i=0; i<NOFACE; i++) fprintf(fp, "         Q:%s", FBUtility::getDirection(i).c_str());
  
  ( Unit_Log == DIMENSIONAL ) ? fprintf(fp, " > Balance[m^3/s] : ") : fprintf(fp, " > Balance[-]     : ");
  
  //for (int i=0; i<NOFACE; i++) fprintf(fp, "         V:%s", FBUtility::getDirection(i).c_str());
  
  if (C->isHeatProblem())
  {
    for (int i=0; i<NOFACE; i++) fprintf(fp, "         H:%s", FBUtility::getDirection(i).c_str());
    
    fprintf(fp, " >     Balance[W]");
  }
  
  fprintf(fp, "\n");
}


// #################################################################
// 物体に働く力の履歴の出力（コンポーネント毎）
bool History::printHistoryForce(const CompoList* cmp, const REAL_TYPE* frc)
{
  for (int n=1; n<=NoCompo; n++)
  {
    if ( cmp[n].getType()==OBSTACLE || cmp[n].getType()==SOLIDREV )
    {
      char fname[128];
      sprintf( fname, "history_force_%s.txt", cmp[n].alias.c_str() );
      
      FILE* fp;
      
      if ( !(fp=fopen(fname, "a")) )
      {
        stamped_printf("\tSorry, can't open '%s' file.\n", fname);
        return false;
      }
      
      fprintf(fp, "%8d %14.6e %12.4e  %12.4e  %12.4e\n",
              step,
              printTime(),
              printForce(frc[3*n+0]),
              printForce(frc[3*n+1]),
              printForce(frc[3*n+2])
              );
      
      fclose(fp);
    }
  }
  return true;
}


// #################################################################
// 物体に働く力の履歴のヘッダー出力（コンポーネント毎）
bool History::printHistoryForceTitle(const CompoList* cmp)
{
  for (int n=1; n<=NoCompo; n++)
  {
    if ( cmp[n].getType()==OBSTACLE || cmp[n].getType()==SOLIDREV )
    {
      char fname[128];
      sprintf( fname, "history_force_%s.txt", cmp[n].alias.c_str() );
      
      FILE* fp;
      
      if ( !(fp=fopen(fname, "w")) )
      {
        stamped_printf("\tSorry, can't open '%s' file.\n", fname);
        return false;
      }
      
      fprintf(fp, "Column_Data_00\n");
      
      if ( Unit_Log == DIMENSIONAL )
      {
        fprintf(fp, "    step      time[sec]          F_x           F_y           F_z\n");
      }
      else
      {
        fprintf(fp, "    step        time[-]          F_x           F_y           F_z\n");
      }

      fclose(fp);
    }
  }
  return true;
}


// #################################################################
// 物体に働く力の統計量出力
bool History::printCompoStatistics(const CompoList* cmp, const REAL_TYPE* frc)
{
  FILE* fp = NULL;
  
  if ( !(fp=fopen("component_statistic.txt", "w")) )
  {
    stamped_printf("\tSorry, can't open 'component_statistic.txt' file.\n");
    return false;
  }
  
  fprintf(fp, "   No           Compo. Label          F_x           F_y           F_z\n");
  
  for (int i=1; i<=NoCompo; i++)
  {
    if ( cmp[i].getType() )
    {
      fprintf(fp, "%5d %22s %12.4e  %12.4e  %12.4e\n",
              i,
              cmp[i].alias.c_str(),
              printForce(frc[3*i+0]),
              printForce(frc[3*i+1]),
              printForce(frc[3*i+2]) );
    }
  }
  fclose(fp);
  return true;
}


// #################################################################
// 物体に働く力の統計量のリスタート
void History::loadCompoStatistics(FILE* fp, const CompoList* cmp, const Control* C, REAL_TYPE* frc)
{
  char buf[1024];
  int idummy;
  
  // 1行目を読み込み
  fgets(buf, 1, fp);
  
  // 2行目以降
  for (int i=1; i<=C->NoCompo; i++)
  {
    if ( cmp[i].getType() )
    {
      fscanf(fp, "%d %s %e %e %e",
              &idummy,
              buf,
              &frc[3*i+0],
              &frc[3*i+1],
              &frc[3*i+2] );
    }
  }
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
