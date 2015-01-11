//##################################################################################
//
// FFV-C : Frontflow / violet Cartesian
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
 * @file   ffv_io_base.C
 * @brief  File IO base Class
 * @author aics
 */

#include "ffv_io_base.h"

#include "FileCommon.h"
#include "BitVoxel.h"
#include "RLE.h"
#include "FileSystemUtil.h"
#include "type.h"
#include "BlockSaver.h"

#include "ffv_LSfunc.h"
#include "ffv_Ffunc.h"
                             
                             
// #################################################################
// ファイル入出力に関するパラメータを取得し，出力の並列モードを指定, PLOT3Dバッファサイズ
// @pre Control::getTimeControl()
void IO_BASE::getFIOparams()
{
  
  REAL_TYPE f_val=0.0;
  string str;
  string label, leaf;
  
  
  // Default setting
  IOmode = IO_DISTRIBUTE;
  
  // 逐次実行の場合には、強制的に IO_GATHER
  if ( (C->Parallelism == Control::Serial) || (C->Parallelism == Control::OpenMP) )
  {
    IOmode = IO_GATHER;
  }
  
  
  
  // 基本変数の瞬時値データ
  
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
      C->Interval[Control::tg_basic].setMode(IntervalManager::By_step);
      C->Interval[Control::tg_derived].setMode(IntervalManager::By_step);
    }
    else if( !strcasecmp(str.c_str(), "time") )
    {
      C->Interval[Control::tg_basic].setMode(IntervalManager::By_time);
      C->Interval[Control::tg_derived].setMode(IntervalManager::By_time);
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
      C->Interval[Control::tg_basic].setInterval((double)f_val);
      C->Interval[Control::tg_derived].setInterval((double)f_val);
    }
  }
  
  
  switch ( Format )
  {
    case sph_fmt:
      getFormatOption("sph");
      break;
      
    case bov_fmt:
      getFormatOption("bov");
      break;
      
    case plt3d_fun_fmt:
      getFormatOption("plot3d");
      break;
  }
  
  
  // 固有のオプション
  getInherentOption();
  
  
  
  // 全圧
  label="/Output/Data/BasicVariables/TotalPressure";
  if ( tpCntl->chkLabel(label) )
  {
    if ( !(tpCntl->getInspectedValue(label, str )) )
    {
      Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
      Exit(0);
    }
    
    if     ( !strcasecmp(str.c_str(), "on") )
    {
      C->varState[var_TotalP] = ON;
      C->NvarsIns_plt3d += 1;
    }
    else if( !strcasecmp(str.c_str(), "off") ) C->varState[var_TotalP] = OFF;
    else
    {
      Hostonly_ stamped_printf("\tInvalid keyword is described for '%s'\n", label.c_str());
      Exit(0);
    }
  }

  
  // 渦度ベクトル
  label="/Output/Data/BasicVariables/Vorticity";
  if ( tpCntl->chkLabel(label) )
  {
    if ( !(tpCntl->getInspectedValue(label, str )) )
    {
      Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
      Exit(0);
    }
    
    if     ( !strcasecmp(str.c_str(), "on") )
    {
      C->varState[var_Vorticity] = ON;
      C->NvarsIns_plt3d += 3;
    }
    else if( !strcasecmp(str.c_str(), "off") ) C->varState[var_Vorticity] = OFF;
    else
    {
      Hostonly_ stamped_printf("\tInvalid keyword is described for '%s'\n", label.c_str());
      Exit(0);
    }
  }

  
  // 速度勾配テンソルの第2不変量
  label="/Output/Data/BasicVariables/Qcriterion";
  if ( tpCntl->chkLabel(label) )
  {
    if ( !(tpCntl->getInspectedValue(label, str )) )
    {
      Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
      Exit(0);
    }
    
    if     ( !strcasecmp(str.c_str(), "on") )
    {
      C->varState[var_Qcr] = ON;
      C->NvarsIns_plt3d += 1;
    }
    else if( !strcasecmp(str.c_str(), "off") ) C->varState[var_Qcr] = OFF;
    else
    {
      Hostonly_ stamped_printf("\tInvalid keyword is described for '%s'\n", label.c_str());
      Exit(0);
    }
  }

  
  // ヘリシティ
  label="/Output/Data/BasicVariables/Helicity";
  if ( tpCntl->chkLabel(label) )
  {
    if ( !(tpCntl->getInspectedValue(label, str )) )
    {
      Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
      Exit(0);
    }
    
    if     ( !strcasecmp(str.c_str(), "on") )
    {
      C->varState[var_Helicity] = ON;
      C->NvarsIns_plt3d += 1;
    }
    else if( !strcasecmp(str.c_str(), "off") ) C->varState[var_Helicity] = OFF;
    else
    {
      Hostonly_ stamped_printf("\tInvalid keyword is described for '%s'\n", label.c_str());
      Exit(0);
    }
  }

  
  // 発散値
  label="/Output/Data/BasicVariables/Divergence";
  if ( tpCntl->chkLabel(label) )
  {
    if ( !(tpCntl->getInspectedValue(label, str )) )
    {
      Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
      Exit(0);
    }
    
    if     ( !strcasecmp(str.c_str(), "on") )
    {
      C->varState[var_Div] = ON;
      C->NvarsIns_plt3d += 1;
    }
    else if( !strcasecmp(str.c_str(), "off") ) C->varState[var_Div] = OFF;
    else
    {
      Hostonly_ stamped_printf("\tInvalid keyword is described for '%s'\n", label.c_str());
      Exit(0);
    }
  }
  
  

  
  
  // 統計値操作に関するパラメータを取得
  if ( C->Mode.Statistic == ON )
  {
    label = "/Output/Data/StatisticalVariables/TemporalType";
    
    if ( !(tpCntl->getInspectedValue(label, str )) )
    {
      Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
      Exit(0);
    }
    else
    {
      if     ( !strcasecmp(str.c_str(), "step") )
      {
        if ( C->Interval[Control::tg_statistic].getMode() == IntervalManager::By_time )
        {
          Hostonly_ stamped_printf("\tError : Specified temporal mode is not consistent with '/TimeControl/Statistic/TemporalType'\n");
          Exit(0);
        }
      }
      else if( !strcasecmp(str.c_str(), "time") )
      {
        if ( C->Interval[Control::tg_statistic].getMode() == IntervalManager::By_step )
        {
          Hostonly_ stamped_printf("\tError : Specified temporal mode is not consistent with '/TimeControl/Statistic/TemporalType'\n");
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
    label="/Output/Data/StatisticalVariables/Interval";
    
    if ( !(tpCntl->getInspectedValue(label, val )) )
    {
      Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
      Exit(0);
    }
    else
    {
      C->Interval[Control::tg_statistic].setInterval(val);
    }
    
    // Statistic for velocity
    label="/Output/Data/StatisticalVariables/VelocityStat";
    if ( tpCntl->chkLabel(label) )
    {
      if ( !(tpCntl->getInspectedValue(label, str )) )
      {
        Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
        Exit(0);
      }
      
      if ( !strcasecmp(str.c_str(), "on") )
      {
        C->Mode.StatVelocity = ON;
        
        C->varState[var_RmsV] = ON;
        C->NvarsAvr_plt3d += 3;
        
        C->varState[var_RmsMeanV] = ON;
        C->NvarsAvr_plt3d += 3;
      }
      else if( !strcasecmp(str.c_str(), "off") )
      {
        C->Mode.StatVelocity = OFF;
        C->varState[var_RmsV] = OFF;
        C->varState[var_RmsMeanV] = OFF;
      }
      else
      {
        Hostonly_ stamped_printf("\tInvalid keyword is described for '%s'\n", label.c_str());
        Exit(0);
      }
    }
    
    // Statistic for pressure
    label="/Output/Data/StatisticalVariables/PressureStat";
    if ( tpCntl->chkLabel(label) )
    {
      if ( !(tpCntl->getInspectedValue(label, str )) )
      {
        Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
        Exit(0);
      }
      
      if ( !strcasecmp(str.c_str(), "on") )
      {
        C->Mode.StatPressure = ON;
        
        C->varState[var_RmsP] = ON;
        C->NvarsAvr_plt3d += 1;
        
        C->varState[var_RmsMeanP] = ON;
        C->NvarsAvr_plt3d += 1;
      }
      else if( !strcasecmp(str.c_str(), "off") )
      {
        C->Mode.StatPressure = OFF;
        C->varState[var_RmsP] = OFF;
        C->varState[var_RmsMeanP] = OFF;
      }
      else
      {
        Hostonly_ stamped_printf("\tInvalid keyword is described for '%s'\n", label.c_str());
        Exit(0);
      }
    }
    
    // Statistic for temperature
    if ( C->isHeatProblem() )
    {
      label="/Output/Data/StatisticalVariables/TemperatureStat";
      if ( tpCntl->chkLabel(label) )
      {
        if ( !(tpCntl->getInspectedValue(label, str )) )
        {
          Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
          Exit(0);
        }
        
        if ( !strcasecmp(str.c_str(), "on") )
        {
          C->Mode.StatTemperature = ON;
          
          C->varState[var_RmsT] = ON;
          C->NvarsAvr_plt3d += 1;
          
          C->varState[var_RmsMeanT] = ON;
          C->NvarsAvr_plt3d += 1;
        }
        else if( !strcasecmp(str.c_str(), "off") )
        {
          C->Mode.StatTemperature = OFF;
          C->varState[var_RmsT] = OFF;
          C->varState[var_RmsMeanT] = OFF;
        }
        else
        {
          Hostonly_ stamped_printf("\tInvalid keyword is described for '%s'\n", label.c_str());
          Exit(0);
        }
      }
    }
    
  }
  

  
  
  // ボクセルファイル出力 (Hidden)
  IO_Voxel = OFF;
  label = "/GeometryModel/VoxelOutput";
  
  if ( tpCntl->chkLabel(label) )
  {
    if ( tpCntl->getInspectedValue(label, str) )
    {
      if     ( !strcasecmp(str.c_str(), "svx") )  IO_Voxel = voxel_SVX;
      else if( !strcasecmp(str.c_str(), "bvx") )  IO_Voxel = voxel_BVX;
      else if( !strcasecmp(str.c_str(), "off") )  IO_Voxel = OFF;
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
  
  
  // BCflag出力 (Hidden)
  IO_BCflag = OFF;
  label = "/GeometryModel/BCflagOutput";
  
  if ( tpCntl->chkLabel(label) )
  {
    if ( tpCntl->getInspectedValue(label, str) )
    {
      if     ( !strcasecmp(str.c_str(), "on") )  IO_BCflag = ON;
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
  
  
  
  // デバッグ用出力 (Hidden)

  label="/Output/Data/VTKoption";
  
  if ( tpCntl->chkLabel(label) )
  {
    if ( tpCntl->getInspectedValue(label, str) )
    {
      if     ( !strcasecmp(str.c_str(), "on") )  output_vtk = ON;
      else if( !strcasecmp(str.c_str(), "off") ) output_vtk = OFF;
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
  
  
  label="/Output/Data/OutputDebug";
  
  if ( tpCntl->chkLabel(label) )
  {
    if ( tpCntl->getInspectedValue(label, str) )
    {
      if     ( !strcasecmp(str.c_str(), "on") )  output_debug = ON;
      else if( !strcasecmp(str.c_str(), "off") ) output_debug = OFF;
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
// ファイルフォーマットのオプションを指定する
void IO_BASE::getFormatOption(const string form)
{
  string str;
  string label;
  string dir;
  int ct;
  
  
  // ディレクトリのチェック
  dir = "/Output/FormatOption/" + form;
  
  if ( !(tpCntl->chkNode(dir)) )
  {
    Hostonly_ stamped_printf("\tParsing error : Missing '%s'\n", dir.c_str());
    Exit(0);
  }
  
  
  // 出力ガイドセルモード
  label = dir + "/GuideOut";
  
  if ( !(tpCntl->getInspectedValue(label, ct)) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    Exit(0);
    Exit(0);
  }
  if ( ct < 0 || ct > guide )
  {
    Hostonly_ stamped_printf("\tInvalid range of guide out /Output/FormatOption/*/GuideOut : 0 <= guideout <= %d\n", guide);
    Exit(0);
  }
  C->GuideOut = GuideOut = ct;

  
  // Output Directory_Path
  label = dir + "/DirectoryPath";
  
  if ( !(tpCntl->getInspectedValue(label, str)) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    Exit(0);
  }
  // 指定が無ければ，空のまま
  if ( !str.empty() )
  {
    OutDirPath = str;
  }
  
  
  // TimeSlice option
  label = dir + "/TimeSlice";
  
  if ( !(tpCntl->getInspectedValue(label, str)) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    Exit(0);
  }
  
  if ( !strcasecmp(str.c_str(), "on") )
  {
    Slice = ON;
  }
  else
  {
    Slice = OFF;
  }
  
  // 1プロセスの場合にはランク番号がないので、タイムスライス毎のディレクトリは作らない
  if ( (C->Parallelism == Control::Serial) || (C->Parallelism == Control::OpenMP) )
  {
    Slice = OFF;
  }
  
}



// #################################################################
// ステージングオプション
void IO_BASE::getStagingOption()
{
  string str;
  string label;
  
  
  // Staging option
  label="/StartCondition/Restart/Staging";
  
  if ( !(tpCntl->getInspectedValue(label, str )) )
  {
    ;
  }
  else
  {
    if     ( !strcasecmp(str.c_str(), "on") )  C->Restart_staging = ON;
    else if( !strcasecmp(str.c_str(), "off") ) C->Restart_staging = OFF;
    else
    {
      Hostonly_ stamped_printf("\tInvalid keyword is described for '%s'\n", label.c_str());
      Exit(0);
    }
  }
  
}


// #################################################################
// 初期値
// @todo セルフェイスの粗格子リスタート  >> 近似なのでサボる？
// @ see getTimeControl()
void IO_BASE::getStartCondition()
{
  string str;
  string label, leaf;
  
  
  // 初期条件 温度はParseBC::getInitTempOfMedium()
  if ( C->Start == initial_start )
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
    
    if ( !(tpCntl->getInspectedValue(label, C->iv.Pressure )) )
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
    C->iv.VecU = v[0];
    C->iv.VecV = v[1];
    C->iv.VecW = v[2];
  }
  
}




// #################################################################
// 制御パラメータSTEERの表示
void IO_BASE::printSteerConditions(FILE* fp)
{
  // File IO mode ------------------
  fprintf(fp,"\n\tFile IO Mode\n");
  
  fprintf(fp,"\t     Unit of File             :   %s\n", (C->Unit.File == DIMENSIONAL) ? "Dimensional" : "Non-Dimensional");
  
  // InputMode >> Distributed by default
  fprintf(fp,"\t     IO Mode                  :   %s\n", (IOmode==IO_GATHER) ? "Gathered" : "Distributed");
  
  // Output guide
  fprintf(fp,"\t     Guide cell for output    :   %d\n", C->GuideOut);
  
  // IO Directory
  fprintf(fp,"\t     I/O Directory Input      :   \"%s\"\n", InDirPath.c_str());
  fprintf(fp,"\t     I/O Directory Output     :   \"%s\"\n", OutDirPath.c_str());
  
  // Time Slice option
  fprintf(fp,"\t     Time Slice Directory     :   %s\n", (Slice==ON) ? "On" : "Off");
  
  
  // 固有パラメータ
  printSteerConditionsInherent(fp);
  
  
  // Hidden
  fprintf(fp,"\t     Voxel model output       :   ");
  switch (IO_Voxel) {
    case voxel_SVX:
      fprintf(fp,"svx\n");
      break;
      
    case voxel_BVX:
      fprintf(fp,"bvx\n");
      break;
      
    default:
      fprintf(fp,"none\n");
      break;
  }
  fprintf(fp,"\t     BC flag output           :   %s\n", (IO_BCflag==ON) ? "On" : "Off");
  fprintf(fp,"\t     VTK output               :   %s\n", (output_vtk==ON) ? "On" : "Off");
  fprintf(fp,"\t     Debug output             :   %s\n", (output_debug==ON) ? "On" : "Off");
}




// #################################################################
/**
 * @brief リスタートの最大値と最小値の表示
 * @param [in]  fp   ファイルポインタ
 * @param [out] flop 浮動小数点演算数
 */
void IO_BASE::RestartDisplayMinmax(FILE* fp, double& flop)
{
  Hostonly_ fprintf(stdout, "\n\tNon-dimensional value\n");
  Hostonly_ fprintf(fp, "\n\tNon-dimensional value\n");
  REAL_TYPE f_min, f_max, vec_min[3], vec_max[3];
  
  // Velocity
  fb_minmax_v_ (vec_min, vec_max, size, &guide, RF->getV00(), d_v, &flop); // allreduceすること
  
  if ( numProc > 1 )
  {
    REAL_TYPE vmin_tmp[3] = {vec_min[0], vec_min[1], vec_min[2]};
    if( paraMngr->Allreduce(vmin_tmp, vec_min, 3, MPI_MIN) != CPM_SUCCESS ) Exit(0);
    
    REAL_TYPE vmax_tmp[3] = {vec_max[0], vec_max[1], vec_max[2]};
    if( paraMngr->Allreduce(vmax_tmp, vec_max, 3, MPI_MAX) != CPM_SUCCESS ) Exit(0);
  }
  
  Hostonly_ fprintf(stdout, "\t\tVelocity U      : min=%13.6e / max=%13.6e\n", vec_min[0], vec_max[0]);
  Hostonly_ fprintf(stdout, "\t\t         V      : min=%13.6e / max=%13.6e\n", vec_min[1], vec_max[1]);
  Hostonly_ fprintf(stdout, "\t\t         W      : min=%13.6e / max=%13.6e\n", vec_min[2], vec_max[2]);
  Hostonly_ fprintf(fp,     "\t\tVelocity U      : min=%13.6e / max=%13.6e\n", vec_min[0], vec_max[0]);
  Hostonly_ fprintf(fp,     "\t\t         V      : min=%13.6e / max=%13.6e\n", vec_min[1], vec_max[1]);
  Hostonly_ fprintf(fp,     "\t\t         W      : min=%13.6e / max=%13.6e\n", vec_min[2], vec_max[2]);
  
  
  // FVelocity
  fb_minmax_v_ (vec_min, vec_max, size, &guide, RF->getV00(), d_vf, &flop); // allreduceすること
  
  if ( numProc > 1 )
  {
    REAL_TYPE vmin_tmp[3] = {vec_min[0], vec_min[1], vec_min[2]};
    if( paraMngr->Allreduce(vmin_tmp, vec_min, 3, MPI_MIN) != CPM_SUCCESS ) Exit(0);
    
    REAL_TYPE vmax_tmp[3] = {vec_max[0], vec_max[1], vec_max[2]};
    if( paraMngr->Allreduce(vmax_tmp, vec_max, 3, MPI_MAX) != CPM_SUCCESS ) Exit(0);
  }
  
  Hostonly_ fprintf(stdout, "\t\tFace Velocity U : min=%13.6e / max=%13.6e\n", vec_min[0], vec_max[0]);
  Hostonly_ fprintf(stdout, "\t\t              V : min=%13.6e / max=%13.6e\n", vec_min[1], vec_max[1]);
  Hostonly_ fprintf(stdout, "\t\t              W : min=%13.6e / max=%13.6e\n", vec_min[2], vec_max[2]);
  Hostonly_ fprintf(fp,     "\t\tFace Velocity U : min=%13.6e / max=%13.6e\n", vec_min[0], vec_max[0]);
  Hostonly_ fprintf(fp,     "\t\t              V : min=%13.6e / max=%13.6e\n", vec_min[1], vec_max[1]);
  Hostonly_ fprintf(fp,     "\t\t              W : min=%13.6e / max=%13.6e\n", vec_min[2], vec_max[2]);
  
  
  // Pressure
  fb_minmax_s_ (&f_min, &f_max, size, &guide, d_p, &flop);
  
  if ( numProc > 1 )
  {
    REAL_TYPE min_tmp = f_min;
    if( paraMngr->Allreduce(&min_tmp, &f_min, 1, MPI_MIN) != CPM_SUCCESS ) Exit(0);
    
    REAL_TYPE max_tmp = f_max;
    if( paraMngr->Allreduce(&max_tmp, &f_max, 1, MPI_MAX) != CPM_SUCCESS ) Exit(0);
  }
  
  Hostonly_ fprintf(stdout, "\t\tPressure        : min=%13.6e / max=%13.6e\n", f_min, f_max);
  Hostonly_ fprintf(fp,     "\t\tPressure        : min=%13.6e / max=%13.6e\n", f_min, f_max);
  
  
  // temperature
  if ( C->isHeatProblem() )
  {
    fb_minmax_s_ (&f_min, &f_max, size, &guide, d_ie, &flop);
    
    if ( numProc > 1 )
    {
      REAL_TYPE min_tmp = f_min;
      if( paraMngr->Allreduce(&min_tmp, &f_min, 1, MPI_MIN) != CPM_SUCCESS ) Exit(0);
      
      REAL_TYPE max_tmp = f_max;
      if( paraMngr->Allreduce(&max_tmp, &f_max, 1, MPI_MAX) != CPM_SUCCESS ) Exit(0);
    }
    
    Hostonly_ fprintf(stdout, "\t\tTemperature    : min=%13.6e / max=%13.6e\n", f_min, f_max);
    Hostonly_ fprintf(fp,     "\t\tTemperature    : min=%13.6e / max=%13.6e\n", f_min, f_max);
  }
  
}



// #################################################################
// 配列ポインタのセット
void IO_BASE::setVarPointers(REAL_TYPE* m_d_p,
                             REAL_TYPE* m_d_v,
                             REAL_TYPE* m_d_vf,
                             REAL_TYPE* m_d_ie,
                             REAL_TYPE* m_d_ws,
                             REAL_TYPE* m_d_wv,
                             REAL_TYPE* m_d_ap,
                             REAL_TYPE* m_d_av,
                             REAL_TYPE* m_d_ae,
                             REAL_TYPE* m_d_dv,
                             REAL_TYPE* m_d_rms_v,
                             REAL_TYPE* m_d_rms_mean_v,
                             REAL_TYPE* m_d_rms_p,
                             REAL_TYPE* m_d_rms_mean_p,
                             REAL_TYPE* m_d_rms_t,
                             REAL_TYPE* m_d_rms_mean_t,
                             int* m_d_bcd,
                             int* m_d_cdf,
                             double* m_mat_tbl,
                             int* m_d_mid,
                             REAL_TYPE* m_d_iob)
{
  d_p          = m_d_p;
  d_v          = m_d_v;
  d_vf         = m_d_vf;
  d_ie         = m_d_ie;
  d_ws         = m_d_ws;
  d_wv         = m_d_wv;
  d_ap         = m_d_ap;
  d_av         = m_d_av;
  d_ae         = m_d_ae;
  d_dv         = m_d_dv;
  d_rms_v      = m_d_rms_v;
  d_rms_mean_v = m_d_rms_mean_v;
  d_rms_p      = m_d_rms_p;
  d_rms_mean_p = m_d_rms_mean_p;
  d_rms_t      = m_d_rms_t;
  d_rms_mean_t = m_d_rms_mean_t;
  d_bcd        = m_d_bcd;
  d_cdf        = m_d_cdf;
  mat_tbl      = m_mat_tbl;
  d_mid        = m_d_mid;
  d_iobuf      = m_d_iob;
}


// #################################################################
// BCflagをbvxで出力する
int IO_BASE::writeBCflag(const int out_gc)
{
  if (IO_BCflag != ON) return true;
  
  unsigned bitWidth = 5;
  int rank = paraMngr->GetMyRankID();
  
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gc = out_gc;
  int gd = guide;
  
  size_t nx = (ix+2*gc) * (jx+2*gc) * (kx+2*gc);
  
  // unsignd int
  unsigned* buf = new unsigned[nx];
  
  // start index Active, State bitはマスクする >> 30bitのみ
  unsigned val = (unsigned)(d_cdf[ _F_IDX_S3D(1-gc, 1-gc, 1-gc, ix, jx, kx, gd) ] & 0x3fffffff);
  int c=0;
  
#pragma omp parallel for firstprivate(ix, jx, kx, gc, gd, val) schedule(static) reduction(+:c)
  for (int k=1-gc; k<=kx+gc; k++) {
    for (int j=1-gc; j<=jx+gc; j++) {
      for (int i=1-gc; i<=ix+gc; i++) {
        size_t m0 = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        size_t m1 = _F_IDX_S3D(i, j, k, ix, jx, kx, gc);
        unsigned tmp = (unsigned)(d_cdf[m0] & 0x3fffffff);
        buf[m1] = tmp;
        if ( tmp != val ) c++;
      }
    }
  }
  
  bool ret = false;
  int ret_val;
  
  // サブドメイン内が同じ値の時には、BCflag配列を書き出さない
  if ( c>0 )
  {
    ret = BVX_IO::Save_Block_BCflag(size, gc, bitWidth, rank, OutDirPath, buf, BVXcomp);
    if ( !ret )
    {
      stamped_printf("\tError : when saving BCflag\n");
      Exit(0);
    }
    ret_val = -1;
  }
  else
  {
    ret_val = (int)val;
  }
  
  
  if (buf)
  {
    delete [] buf;
    buf = NULL;
  }
  
  return ret_val;
}


// #################################################################
// Cell IDをbvxで出力する
int IO_BASE::writeCellID(const int out_gc)
{
  if (IO_Voxel != voxel_BVX) return true;
  
  unsigned bitWidth = 5;
  int rank = paraMngr->GetMyRankID();
  
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gc = out_gc;
  int gd = guide;
  
  size_t nx = (ix+2*gc) * (jx+2*gc) * (kx+2*gc);
  
  // unsignd char
  u8 *buf = new u8[nx];
  
  // start index 下位5bitの値のみ
  u8 val = DECODE_CMP(d_bcd[ _F_IDX_S3D(1-gc, 1-gc, 1-gc, ix, jx, kx, gd) ]);
  int c=0;
  
#pragma omp parallel for firstprivate(ix, jx, kx, gc, gd, val) schedule(static) reduction(+:c)
  for (int k=1-gc; k<=kx+gc; k++) {
    for (int j=1-gc; j<=jx+gc; j++) {
      for (int i=1-gc; i<=ix+gc; i++) {
        size_t m0 = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        size_t m1 = _F_IDX_S3D(i, j, k, ix, jx, kx, gc);
        u8 tmp = DECODE_CMP(d_bcd[m0]);
        buf[m1] = tmp;
        if ( tmp != val ) c++;
      }
    }
  }
  
  bool ret = false;
  int ret_val;
  
  // サブドメイン内が同じ値の時には、CellID配列を書き出さない
  if ( c>0 )
  {
    ret = BVX_IO::Save_Block_CellID(size, gc, bitWidth, rank, OutDirPath, buf, BVXcomp);
    if ( !ret )
    {
      stamped_printf("\tError : when saving CellID\n");
      Exit(0);
    }
    ret_val = -1;
  }
  else
  {
    ret_val = (int)val;
  }
  

  if (buf)
  {
    delete [] buf;
    buf = NULL;
  }

  return ret_val;
}


// #################################################################
// sphファイルの書き出し（内部領域のみ）
void IO_BASE::writeRawSPH(const REAL_TYPE *vf,
                          const int* sz,
                          const int gc,
                          const int gc_out,
                          const REAL_TYPE* org,
                          const REAL_TYPE* ddx,
                          const int m_ModePrecision)
{
  int pad, dType, stp, svType;
  REAL_TYPE ox, oy, oz, dx, dy, dz, tm;
  long long szl[3], stpl;
  
  
  char sph_fname[512];
  
  if ( paraMngr->IsParallel() )
  {
    sprintf( sph_fname, "field%010d.sph", paraMngr->GetMyRankID() );
  }
  else
  {
    sprintf( sph_fname, "field.sph" );
  }
  
  ofstream ofs(sph_fname, ios::out | ios::binary);
  if (!ofs)
  {
    printf("\tCan't open %s file\n", sph_fname);
    Exit(0);
  }
  
  int m_sz[3];
  m_sz[0] = sz[0]+2*gc_out;
  m_sz[1] = sz[1]+2*gc_out;
  m_sz[2] = sz[2]+2*gc_out;
  int gd = gc;
  
  size_t nx = m_sz[0] * m_sz[1] * m_sz[2];
  
  ox = org[0]-ddx[0]*(REAL_TYPE)gc_out;
  oy = org[1]-ddx[1]*(REAL_TYPE)gc_out;
  oz = org[2]-ddx[2]*(REAL_TYPE)gc_out;
  dx = ddx[0];
  dy = ddx[1];
  dz = ddx[2];
  //printf("org: %f %f %f\n", ox, oy, oz);
  //printf("dx : %f %f %f\n", dx, dy, dz);
  
  svType = kind_scalar;
  if ( sizeof(REAL_TYPE) == sizeof(double) )
  {
    for (int i=0; i<3; i++)   szl[i] = (long long)m_sz[i];
  }
  
  REAL_TYPE *f = new REAL_TYPE[nx];
  
  size_t m, l;
  
  for (int k=1; k<=m_sz[2]; k++) {
    for (int j=1; j<=m_sz[1]; j++) {
      for (int i=1; i<=m_sz[0]; i++) {
        l = _F_IDX_S3D(i, j, k, m_sz[0], m_sz[1], m_sz[2], gc_out);
        m = _F_IDX_S3D(i, j, k, m_sz[0], m_sz[1], m_sz[2], gd);
        f[l] = (REAL_TYPE)vf[m];
      }
    }
  }
  
  // data property
  ( m_ModePrecision == sizeof(float) ) ? dType=1 : dType=2;
  pad = sizeof(int)*2;
  ofs.write( (char*)&pad, sizeof(int) );
  ofs.write( (char*)&svType, sizeof(int) );
  ofs.write( (char*)&dType, sizeof(int) );
  ofs.write( (char*)&pad, sizeof(int) );
  
  // voxel size
  if (dType == 1) {
    pad = sizeof(int)*3;
    ofs.write( (char*)&pad, sizeof(int) );
    ofs.write( (char*)&m_sz[0], sizeof(int) );
    ofs.write( (char*)&m_sz[1], sizeof(int) );
    ofs.write( (char*)&m_sz[2], sizeof(int) );
    ofs.write( (char*)&pad, sizeof(int) );
  }
  else {
    pad = sizeof(long long)*3;
    ofs.write( (char*)&pad, sizeof(int) );
    ofs.write( (char*)&szl[0], sizeof(long long) );
    ofs.write( (char*)&szl[1], sizeof(long long) );
    ofs.write( (char*)&szl[2], sizeof(long long) );
    ofs.write( (char*)&pad, sizeof(int) );
  }
  
  // original point of domain
  if (dType == 1) {
    pad = sizeof(float)*3;
    ofs.write( (char*)&pad, sizeof(int) );
    ofs.write( (char*)&ox, sizeof(float) );
    ofs.write( (char*)&oy, sizeof(float) );
    ofs.write( (char*)&oz, sizeof(float) );
    ofs.write( (char*)&pad, sizeof(int) );
  }
  else {
    pad = sizeof(double)*3;
    ofs.write( (char*)&pad, sizeof(int) );
    ofs.write( (char*)&ox, sizeof(double) );
    ofs.write( (char*)&oy, sizeof(double) );
    ofs.write( (char*)&oz, sizeof(double) );
    ofs.write( (char*)&pad, sizeof(int) );
  }
  
  // pitch of voxel
  if (dType == 1) {
    pad = sizeof(float)*3;
    ofs.write( (char*)&pad, sizeof(int) );
    ofs.write( (char*)&dx, sizeof(float) );
    ofs.write( (char*)&dy, sizeof(float) );
    ofs.write( (char*)&dz, sizeof(float) );
    ofs.write( (char*)&pad, sizeof(int) );
  }
  else {
    pad = sizeof(double)*3;
    ofs.write( (char*)&pad, sizeof(int) );
    ofs.write( (char*)&dx, sizeof(double) );
    ofs.write( (char*)&dy, sizeof(double) );
    ofs.write( (char*)&dz, sizeof(double) );
    ofs.write( (char*)&pad, sizeof(int) );
  }
  
  // time stamp
  if (dType == 1) {
    stp = 0;
    tm = 0.0;
    pad = sizeof(int)+sizeof(float);
    ofs.write( (char*)&pad, sizeof(int) );
    ofs.write( (char*)&stp, sizeof(int) );
    ofs.write( (char*)&tm, sizeof(float) );
    ofs.write( (char*)&pad, sizeof(int) );
  }
  else {
    stpl =0;
    tm = 0.0;
    pad = sizeof(long long)+sizeof(double);
    ofs.write( (char*)&pad, sizeof(int) );
    ofs.write( (char*)&stpl, sizeof(long long) );
    ofs.write( (char*)&tm, sizeof(double) );
    ofs.write( (char*)&pad, sizeof(int) );
  }
  
  if (svType == kind_scalar) {
    int pad = (m_ModePrecision == sizeof(float)) ? nx * sizeof(float) : nx * sizeof(double);
    ofs.write( (char*)&pad, sizeof(int) );
    ofs.write( (char*)f,   pad );
    ofs.write( (char*)&pad, sizeof(int) );
  }
  
  ofs.close();
  
  if (f) { delete [] f; f=NULL; }
}

// #################################################################
// 例題のモデルをsvxフォーマットで出力する(体積率とID)
bool IO_BASE::writeSVX(REAL_TYPE *vf, int *id)
{
  if ( IO_Voxel != voxel_SVX ) return false;
  
  int    sz, ix, jx, kx;
  size_t m, l;
  float  ox, oy, oz, dx, dy, dz;
  
  char svx_fname[512];
  
  if ( paraMngr->IsParallel() ) {
    sprintf( svx_fname, "example_%06d.svx", paraMngr->GetMyRankID() );
  } else {
    sprintf( svx_fname, "example.svx" );
  }
  ofstream ofs(svx_fname, ios::out | ios::binary);
  if (!ofs) {
    cout << "\tCan't open " << svx_fname << " file" << endl;
    Exit(0);
  }
  
  int imax = size[0];
  int jmax = size[1];
  int kmax = size[2];
  int gd = guide;
  
  ix = imax+2;  // +2 means guide cell for IP model
  jx = jmax+2;
  kx = kmax+2;
  
  size_t nx = ix*jx*kx;
  
  dx = (float)pitch[0]*C->RefLength;
  dy = (float)pitch[1]*C->RefLength;
  dz = (float)pitch[2]*C->RefLength;
  ox = (float)origin[0]*C->RefLength - dx; // 片側1層分をシフト
  oy = (float)origin[1]*C->RefLength - dy;
  oz = (float)origin[2]*C->RefLength - dz;
  
  //stamped_printf("example out org(%e %e %e) dimensional\n", ox, oy, oz);
  
  float *f = new float[nx];
  int   *q = new int[nx];
  
  for (int k=0; k<=(kmax+1); k++) {
    for (int j=0; j<=(jmax+1); j++) {
      for (int i=0; i<=(imax+1); i++) {
        l = (size_t)(ix*jx*k + ix*j + i);
        m = _F_IDX_S3D(i, j, k, imax, jmax, kmax, gd);
        q[l] = id[m];
        f[l] = (float)vf[m];
      }
    }
  }
  
  // voxel size
  sz = sizeof(int)*3;
  ofs.write( (char*)&sz, sizeof(int) );
  ofs.write( (char*)&ix, sizeof(int) );
  ofs.write( (char*)&jx, sizeof(int) );
  ofs.write( (char*)&kx, sizeof(int) );
  ofs.write( (char*)&sz, sizeof(int) );
  
  // original point of domain
  sz = sizeof(float)*3;
  ofs.write( (char*)&sz, sizeof(int) );
  ofs.write( (char*)&ox, sizeof(float) );
  ofs.write( (char*)&oy, sizeof(float) );
  ofs.write( (char*)&oz, sizeof(float) );
  ofs.write( (char*)&sz, sizeof(int) );
  
  // pitch of voxel
  ofs.write( (char*)&sz, sizeof(int) );
  ofs.write( (char*)&dx, sizeof(float) );
  ofs.write( (char*)&dy, sizeof(float) );
  ofs.write( (char*)&dz, sizeof(float) );
  ofs.write( (char*)&sz, sizeof(int) );
  
  // type of stored data
  sz = sizeof(int)*1;
  int dtype = 0;
  dtype |= ( 0x1<<0 );  // volume fraction
  dtype |= ( 0x1<<2 );  // medium ID
  ofs.write( (char*)&sz,  sizeof(int) );
  ofs.write( (char*)&dtype, sizeof(int) );
  ofs.write( (char*)&sz,  sizeof(int) );
  
  // volume fraction
  sz = nx * sizeof(float);
  ofs.write( (char*)&sz, sizeof(int) );
  ofs.write( (char*)f,   sz );
  ofs.write( (char*)&sz, sizeof(int) );
  
  // medium ID
  sz = nx * sizeof(int);
  ofs.write( (char*)&sz, sizeof(int) );
  ofs.write( (char*)q,   sz );
  ofs.write( (char*)&sz, sizeof(int) );
  
  ofs.close();
  
  if (f) { delete [] f; f=NULL; }
  if (q) { delete [] q; q=NULL; }
  
  return true;
}


// #################################################################
// 例題のモデルをsvxフォーマットで出力する(ID)
bool IO_BASE::writeSVX(const int* bcd)
{
  if ( IO_Voxel != voxel_SVX ) return false;
  
  int   sz, ix, jx, kx;
  size_t m, l;
  float ox, oy, oz, dx, dy, dz;
  
  char svx_fname[512];
  
  if ( paraMngr->IsParallel() ) {
    sprintf( svx_fname, "example_%06d.svx", paraMngr->GetMyRankID() );
  } else {
    sprintf( svx_fname, "example.svx" );
  }
  ofstream ofs(svx_fname, ios::out | ios::binary);
  if (!ofs) {
    cout << "\tCan't open " << svx_fname << " file" << endl;
    Exit(0);
  }
  
  int imax = size[0];
  int jmax = size[1];
  int kmax = size[2];
  int gd = guide;
  
  ix = imax+2;  // +2 means guide cell for IP model
  jx = jmax+2;
  kx = kmax+2;
  
  size_t nx = (size_t)(ix*jx*kx);
  
  dx = (float)pitch[0]*C->RefLength;
  dy = (float)pitch[1]*C->RefLength;
  dz = (float)pitch[2]*C->RefLength;
  ox = (float)origin[0]*C->RefLength - dx; // 片側1層分をシフト
  oy = (float)origin[1]*C->RefLength - dy;
  oz = (float)origin[2]*C->RefLength - dz;
  
  //stamped_printf("example out org(%e %e %e) dimensional\n", ox, oy, oz);
  
  int   *q = new int[nx];
  
  for (int k=0; k<=(kmax+1); k++) {
    for (int j=0; j<=(jmax+1); j++) {
      for (int i=0; i<=(imax+1); i++) {
        l = (size_t)(ix*jx*k + ix*j + i);
        m = _F_IDX_S3D(i, j, k, imax, jmax, kmax, gd);
        q[l] = DECODE_CMP(bcd[m]);
      }
    }
  }
  
  // voxel size
  sz = sizeof(int)*3;
  ofs.write( (char*)&sz, sizeof(int) );
  ofs.write( (char*)&ix, sizeof(int) );
  ofs.write( (char*)&jx, sizeof(int) );
  ofs.write( (char*)&kx, sizeof(int) );
  ofs.write( (char*)&sz, sizeof(int) );
  
  // original point of domain
  sz = sizeof(float)*3;
  ofs.write( (char*)&sz, sizeof(int) );
  ofs.write( (char*)&ox, sizeof(float) );
  ofs.write( (char*)&oy, sizeof(float) );
  ofs.write( (char*)&oz, sizeof(float) );
  ofs.write( (char*)&sz, sizeof(int) );
  
  // pitch of voxel
  ofs.write( (char*)&sz, sizeof(int) );
  ofs.write( (char*)&dx, sizeof(float) );
  ofs.write( (char*)&dy, sizeof(float) );
  ofs.write( (char*)&dz, sizeof(float) );
  ofs.write( (char*)&sz, sizeof(int) );
  
  // type of stored data
  sz = sizeof(int)*1;
  int dtype = 0;
  dtype |= ( 0x1<<2 );  // medium ID
  ofs.write( (char*)&sz,  sizeof(int) );
  ofs.write( (char*)&dtype, sizeof(int) );
  ofs.write( (char*)&sz,  sizeof(int) );
  
  // medium ID
  sz = nx * sizeof(int);
  ofs.write( (char*)&sz, sizeof(int) );
  ofs.write( (char*)q,   sz );
  ofs.write( (char*)&sz, sizeof(int) );
  
  ofs.close();
  
  if (q) { delete [] q; q=NULL; }
  
  return true;
}

