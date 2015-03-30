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
 * @file   ffv_sph.C
 * @brief  File IO of SPH Class
 * @author aics
 */

#include "ffv_sph.h"

#include "ffv_LSfunc.h"
#include "ffv_Ffunc.h"


// #################################################################
// リスタートのDFIファイル
// @todo セルフェイスの粗格子リスタート  >> 近似なのでサボる？
// @see getTimeControl()
void SPH::getRestartDFI()
{
  string str;
  string label, leaf;
  
  
  // リスタート時のDFIファイル名
  if ( C->Start != initial_start )
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
    
    
    if ( C->isHeatProblem() )
    {
      label="/StartCondition/Restart/DFIfiles/Temperature";
      
      if ( tpCntl->getInspectedValue(label, str ) )
      {
        f_dfi_in_temp = str.c_str();
      }
      if ( f_dfi_in_temp.empty() == true ) f_dfi_in_temp = "tmp";
    }
    
    
    // 統計値
    if ( C->Mode.Statistic == ON )
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
      
      
      if ( C->isHeatProblem() )
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
  
}


// #################################################################
/* @brief ファイル出力の初期化
 */
void SPH::initFileOut(const int id_cell, const int id_bcf)
{
  // Format
  CDM::E_CDM_FORMAT cdm_format = CDM::E_CDM_FMT_SPH;

  
  // Datatype
  CDM::E_CDM_DTYPE datatype;
  
  if ( sizeof(REAL_TYPE) == 4 )
  {
    datatype = CDM::E_CDM_FLOAT32;
  }
  else if ( sizeof(REAL_TYPE) == 8 )
  {
    datatype = CDM::E_CDM_FLOAT64;
  }
  else
  {
    Exit(0);
  }
  
  
  // 出力ファイルヘッダ
  int cdm_tail[3], cdm_div[3];
  for (int i=0; i<3; i++) cdm_tail[i]=size[i];
  for (int i=0; i<3; i++) cdm_div[i]=1;
  
  if ( numProc > 1)
  {
    const int* p_tail = paraMngr->GetVoxelTailIndex();
    for (int i=0; i<3; i++ ) cdm_tail[i]=p_tail[i]+1;
    
    const int* p_div = paraMngr->GetDivNum();
    for (int i=0; i<3; i++ ) cdm_div[i] = p_div[i];
  }
  
  
  int gc_out = C->GuideOut;
  REAL_TYPE cdm_org[3], cdm_pit[3];
  
  for (int i=0; i<3; i++)
  {
    cdm_org[i] = origin[i]; // CDM用オリジナルポイントのセット
    cdm_pit[i] = pitch[i];
  }
  
  
  // 出力ファイルの指定が有次元の場合
  if ( C->Unit.File == DIMENSIONAL )
  {
    for (int i=0; i<3; i++)
    {
      cdm_org[i] *= C->RefLength;
      cdm_pit[i] *= C->RefLength;
    }
  }
  
  // make output directory
  std::string path = OutDirPath;
  
  
  // タイムスライス出力オプション
  CDM::E_CDM_ONOFF TimeSliceDir;
  if ( Slice == ON )
  {
    TimeSliceDir = CDM::E_CDM_ON;
  }
  else
  {
    TimeSliceDir = CDM::E_CDM_OFF;
  }
  
  
  std::string process="./proc.dfi";
  int comp = 1;
  
  
  // 単位の設定
  
  std::string UnitL;
  switch (C->Unit.Length)
  {
    case LTH_ND:
      UnitL = "NonDimensional";
      break;
      
    case LTH_m:
      UnitL = "M";
      break;
  }
  
  
  std::string UnitV;
  if ( C->Unit.File == DIMENSIONAL )
  {
    UnitV = "m/s";
  }
  else
  {
    UnitV = "NonDimensional";
  }
  
  
  std::string UnitP;
  if ( C->Unit.File == DIMENSIONAL )
  {
    UnitP = "Pa";
  }
  else
  {
    UnitP = "NonDimensional";
  }
  
  double DiffPrs = (double)C->RefDensity * (double)C->RefVelocity * (double)C->RefVelocity;
  
  
  std::string UnitT = "Celsius";
  
  
  // Pressure
  comp = 1;
  DFI_OUT_PRS = cdm_DFI::WriteInit(MPI_COMM_WORLD,
                                   cdm_DFI::Generate_DFI_Name(f_Pressure),
                                   path,
                                   f_Pressure,
                                   cdm_format,
                                   gc_out,
                                   datatype,
                                   comp,
                                   process,
                                   G_size,
                                   cdm_pit,
                                   cdm_org,
                                   cdm_div,
                                   head,
                                   cdm_tail,
                                   HostName,
                                   TimeSliceDir);
  
  if ( DFI_OUT_PRS == NULL )
  {
    Hostonly_ stamped_printf("\tFails to instance Pressure dfi.\n");
    Exit(0);
  }
  
  if( Slice == ON )
  {
    DFI_OUT_PRS->SetTimeSliceFlag(CDM::E_CDM_ON);
  }
  else
  {
    DFI_OUT_PRS->SetTimeSliceFlag(CDM::E_CDM_OFF);
  }
  
  
  DFI_OUT_PRS->AddUnit("Length",   UnitL, (double)C->RefLength);
  DFI_OUT_PRS->AddUnit("Velocity", UnitV, (double)C->RefVelocity);
  DFI_OUT_PRS->AddUnit("Pressure", UnitP, (double)C->BasePrs, DiffPrs, true);
  
  DFI_OUT_PRS->WriteProcDfiFile(MPI_COMM_WORLD, true, id_cell, id_bcf);
  
  
  
  // Velocity
  comp = 3;
  DFI_OUT_VEL = cdm_DFI::WriteInit(MPI_COMM_WORLD,
                                   cdm_DFI::Generate_DFI_Name(f_Velocity),
                                   path,
                                   f_Velocity,
                                   cdm_format,
                                   gc_out,
                                   datatype,
                                   comp,
                                   process,
                                   G_size,
                                   cdm_pit,
                                   cdm_org,
                                   cdm_div,
                                   head,
                                   cdm_tail,
                                   HostName,
                                   TimeSliceDir);
  
  if ( DFI_OUT_VEL == NULL )
  {
    Hostonly_ stamped_printf("\tFails to instance Velocity dfi.\n");
    Exit(0);
  }
  
  if ( Slice == ON )
  {
    DFI_OUT_VEL->SetTimeSliceFlag(CDM::E_CDM_ON);
  }
  else
  {
    DFI_OUT_VEL->SetTimeSliceFlag(CDM::E_CDM_OFF);
  }
  
  DFI_OUT_VEL->AddUnit("Length"  , UnitL, (double)C->RefLength);
  DFI_OUT_VEL->AddUnit("Velocity", UnitV, (double)C->RefVelocity);
  DFI_OUT_VEL->AddUnit("Pressure", UnitP, (double)C->BasePrs, DiffPrs, true);
  
  
  
  // Fvelocity
  comp = 3;
  DFI_OUT_FVEL = cdm_DFI::WriteInit(MPI_COMM_WORLD,
                                    cdm_DFI::Generate_DFI_Name(f_Fvelocity),
                                    path,
                                    f_Fvelocity,
                                    cdm_format,
                                    gc_out,
                                    datatype,
                                    comp,
                                    process,
                                    G_size,
                                    cdm_pit,
                                    cdm_org,
                                    cdm_div,
                                    head,
                                    cdm_tail,
                                    HostName,
                                    TimeSliceDir);
  
  if ( DFI_OUT_FVEL == NULL )
  {
    Hostonly_ stamped_printf("\tFails to instance Fvelocity dfi.\n");
    Exit(0);
  }
  
  if ( Slice == ON )
  {
    DFI_OUT_FVEL->SetTimeSliceFlag(CDM::E_CDM_ON);
  }
  else
  {
    DFI_OUT_FVEL->SetTimeSliceFlag(CDM::E_CDM_OFF);
  }
  
  DFI_OUT_FVEL->AddUnit("Length"  , UnitL, (double)C->RefLength);
  DFI_OUT_FVEL->AddUnit("Velocity", UnitV, (double)C->RefVelocity);
  DFI_OUT_FVEL->AddUnit("Pressure", UnitP, (double)C->BasePrs, DiffPrs, true);
  
  
  
  // Temperature
  if ( C->isHeatProblem() )
  {
    comp = 1;
    DFI_OUT_TEMP = cdm_DFI::WriteInit(MPI_COMM_WORLD,
                                      cdm_DFI::Generate_DFI_Name(f_Temperature),
                                      path,
                                      f_Temperature,
                                      cdm_format,
                                      gc_out,
                                      datatype,
                                      comp,
                                      process,
                                      G_size,
                                      cdm_pit,
                                      cdm_org,
                                      cdm_div,
                                      head,
                                      cdm_tail,
                                      HostName,
                                      TimeSliceDir);
    
    if ( DFI_OUT_TEMP == NULL )
    {
      Hostonly_ stamped_printf("\tFails to instance Temperature dfi.\n");
      Exit(0);
    }
    
    if ( Slice == ON )
    {
      DFI_OUT_TEMP->SetTimeSliceFlag(CDM::E_CDM_ON);
    }
    else
    {
      DFI_OUT_TEMP->SetTimeSliceFlag(CDM::E_CDM_OFF);
    }
    
    DFI_OUT_TEMP->AddUnit("Length"  , UnitL, (double)C->RefLength);
    DFI_OUT_TEMP->AddUnit("Velocity", UnitV, (double)C->RefVelocity);
    DFI_OUT_TEMP->AddUnit("Pressure", UnitP, (double)C->BasePrs, DiffPrs, true);
    DFI_OUT_TEMP->AddUnit("Temperature", UnitT, (double)C->BaseTemp, (double)C->DiffTemp, true);
  }
  
  
  
  // 平均値
  if ( C->Mode.Statistic == ON )
  {
    
    // Pressure
    comp = 1;
    DFI_OUT_PRSA = cdm_DFI::WriteInit(MPI_COMM_WORLD,
                                      cdm_DFI::Generate_DFI_Name(f_AvrPressure),
                                      path,
                                      f_AvrPressure,
                                      cdm_format,
                                      gc_out,
                                      datatype,
                                      comp,
                                      process,
                                      G_size,
                                      cdm_pit,
                                      cdm_org,
                                      cdm_div,
                                      head,
                                      cdm_tail,
                                      HostName,
                                      TimeSliceDir);
    
    if ( DFI_OUT_PRSA == NULL )
    {
      Hostonly_ stamped_printf("\tFails to instance AvrPressure dfi.\n");
      Exit(0);
    }
    
    if ( Slice == ON )
    {
      DFI_OUT_PRSA->SetTimeSliceFlag(CDM::E_CDM_ON);
    }
    else
    {
      DFI_OUT_PRSA->SetTimeSliceFlag(CDM::E_CDM_OFF);
    }
    
    DFI_OUT_PRSA->AddUnit("Length"  , UnitL, (double)C->RefLength);
    DFI_OUT_PRSA->AddUnit("Velocity", UnitV, (double)C->RefVelocity);
    DFI_OUT_PRSA->AddUnit("Pressure", UnitP, (double)C->BasePrs, DiffPrs, true);
    
    
    
    // Velocity
    comp = 3;
    DFI_OUT_VELA = cdm_DFI::WriteInit(MPI_COMM_WORLD,
                                      cdm_DFI::Generate_DFI_Name(f_AvrVelocity),
                                      path,
                                      f_AvrVelocity,
                                      cdm_format,
                                      gc_out,
                                      datatype,
                                      comp,
                                      process,
                                      G_size,
                                      cdm_pit,
                                      cdm_org,
                                      cdm_div,
                                      head,
                                      cdm_tail,
                                      HostName,
                                      TimeSliceDir);
    
    if ( DFI_OUT_VELA == NULL )
    {
      Hostonly_ stamped_printf("\tcan not instance AvrVelocity dfi\n");
      Exit(0);
    }
    
    if ( Slice == ON )
    {
      DFI_OUT_VELA->SetTimeSliceFlag(CDM::E_CDM_ON);
    }
    else
    {
      DFI_OUT_VELA->SetTimeSliceFlag(CDM::E_CDM_OFF);
    }
    
    DFI_OUT_VELA->AddUnit("Length"  , UnitL, (double)C->RefLength);
    DFI_OUT_VELA->AddUnit("Velocity", UnitV, (double)C->RefVelocity);
    DFI_OUT_VELA->AddUnit("Pressure", UnitP, (double)C->BasePrs, DiffPrs, true);
    
    
    
    // Temperature
    if ( C->isHeatProblem() )
    {
      comp = 1;
      DFI_OUT_TEMPA = cdm_DFI::WriteInit(MPI_COMM_WORLD,
                                         cdm_DFI::Generate_DFI_Name(f_AvrTemperature),
                                         path,
                                         f_AvrTemperature,
                                         cdm_format,
                                         gc_out,
                                         datatype,
                                         comp,
                                         process,
                                         G_size,
                                         cdm_pit,
                                         cdm_org,
                                         cdm_div,
                                         head,
                                         cdm_tail,
                                         HostName,
                                         TimeSliceDir);
      
      if ( DFI_OUT_TEMPA == NULL )
      {
        Hostonly_ stamped_printf("\tFails to instance AvrTemperature dfi.\n");
        Exit(0);
      }
      
      if ( Slice == ON )
      {
        DFI_OUT_TEMPA->SetTimeSliceFlag(CDM::E_CDM_ON);
      }
      else
      {
        DFI_OUT_TEMPA->SetTimeSliceFlag(CDM::E_CDM_OFF);
      }
      
      DFI_OUT_TEMPA->AddUnit("Length"  , UnitL, (double)C->RefLength);
      DFI_OUT_TEMPA->AddUnit("Velocity", UnitV, (double)C->RefVelocity);
      DFI_OUT_TEMPA->AddUnit("Pressure", UnitP, (double)C->BasePrs, DiffPrs, true);
      DFI_OUT_TEMP->AddUnit("Temperature", UnitT, (double)C->BaseTemp, (double)C->DiffTemp, true);
    }
  } // statistic
  
  
  
  // Derived Variables
  
  // Total Pressure
  if (C->varState[var_TotalP] == ON )
  {
    comp = 1;
    DFI_OUT_TP = cdm_DFI::WriteInit(MPI_COMM_WORLD,
                                    cdm_DFI::Generate_DFI_Name(f_TotalP),
                                    path,
                                    f_TotalP,
                                    cdm_format,
                                    gc_out,
                                    datatype,
                                    comp,
                                    process,
                                    G_size,
                                    cdm_pit,
                                    cdm_org,
                                    cdm_div,
                                    head,
                                    cdm_tail,
                                    HostName,
                                    TimeSliceDir);
    
    if ( DFI_OUT_TP == NULL )
    {
      Hostonly_ stamped_printf("\tFails to instance TotalPressure dfi.\n");
      Exit(0);
    }
    
    if ( Slice == ON )
    {
      DFI_OUT_TP->SetTimeSliceFlag(CDM::E_CDM_ON);
    }
    else
    {
      DFI_OUT_TP->SetTimeSliceFlag(CDM::E_CDM_OFF);
    }
    
    DFI_OUT_TP->AddUnit("Length"  , UnitL, (double)C->RefLength);
    DFI_OUT_TP->AddUnit("Velocity", UnitV, (double)C->RefVelocity);
    DFI_OUT_TP->AddUnit("Pressure", UnitP, (double)C->BasePrs, DiffPrs, true);
  }
  
  
  
  // Vorticity
  if (C->varState[var_Vorticity] == ON )
  {
    comp = 3;
    DFI_OUT_VRT = cdm_DFI::WriteInit(MPI_COMM_WORLD,
                                     cdm_DFI::Generate_DFI_Name(f_Vorticity),
                                     path,
                                     f_Vorticity,
                                     cdm_format,
                                     gc_out,
                                     datatype,
                                     comp,
                                     process,
                                     G_size,
                                     cdm_pit,
                                     cdm_org,
                                     cdm_div,
                                     head,
                                     cdm_tail,
                                     HostName,
                                     TimeSliceDir);
    if( DFI_OUT_VRT == NULL )
    {
      Hostonly_ stamped_printf("\tFails to instance Vorticity dfi.\n");
      Exit(0);
    }
    if( Slice == ON )
    {
      DFI_OUT_VRT->SetTimeSliceFlag(CDM::E_CDM_ON);
    }
    else
    {
      DFI_OUT_VRT->SetTimeSliceFlag(CDM::E_CDM_OFF);
    }
    
    DFI_OUT_VRT->AddUnit("Length"  , UnitL, (double)C->RefLength);
    DFI_OUT_VRT->AddUnit("Velocity", UnitV, (double)C->RefVelocity);
    DFI_OUT_VRT->AddUnit("Pressure", UnitP, (double)C->BasePrs, DiffPrs, true);
  }
  
  
  
  // 2nd Invariant of VGT
  if ( C->varState[var_Qcr] == ON )
  {
    comp = 1;
    DFI_OUT_I2VGT = cdm_DFI::WriteInit(MPI_COMM_WORLD,
                                       cdm_DFI::Generate_DFI_Name(f_I2VGT),
                                       path,
                                       f_I2VGT,
                                       cdm_format,
                                       gc_out,
                                       datatype,
                                       comp,
                                       process,
                                       G_size,
                                       cdm_pit,
                                       cdm_org,
                                       cdm_div,
                                       head,
                                       cdm_tail,
                                       HostName,
                                       TimeSliceDir);
    
    if ( DFI_OUT_I2VGT == NULL )
    {
      Hostonly_ stamped_printf("\tFails to instance 2nd Invariant of VGT dfi.\n");
      Exit(0);
    }
    
    if ( Slice == ON )
    {
      DFI_OUT_I2VGT->SetTimeSliceFlag(CDM::E_CDM_ON);
    }
    else
    {
      DFI_OUT_I2VGT->SetTimeSliceFlag(CDM::E_CDM_OFF);
    }
    
    DFI_OUT_I2VGT->AddUnit("Length"  , UnitL, (double)C->RefLength);
    DFI_OUT_I2VGT->AddUnit("Velocity", UnitV, (double)C->RefVelocity);
    DFI_OUT_I2VGT->AddUnit("Pressure", UnitP, (double)C->BasePrs, DiffPrs, true);
  }
  
  
  
  // Helicity
  if ( C->varState[var_Helicity] == ON )
  {
    comp = 1;
    DFI_OUT_HLT = cdm_DFI::WriteInit(MPI_COMM_WORLD,
                                     cdm_DFI::Generate_DFI_Name(f_Helicity),
                                     path,
                                     f_Helicity,
                                     cdm_format,
                                     gc_out,
                                     datatype,
                                     comp,
                                     process,
                                     G_size,
                                     cdm_pit,
                                     cdm_org,
                                     cdm_div,
                                     head,
                                     cdm_tail,
                                     HostName,
                                     TimeSliceDir);
    
    if ( DFI_OUT_HLT == NULL )
    {
      Hostonly_ stamped_printf("\tFails to instance Helicity dfi.\n");
      Exit(0);
    }
    
    if ( Slice == ON )
    {
      DFI_OUT_HLT->SetTimeSliceFlag(CDM::E_CDM_ON);
    }
    else
    {
      DFI_OUT_HLT->SetTimeSliceFlag(CDM::E_CDM_OFF);
    }
    
    DFI_OUT_HLT->AddUnit("Length",   UnitL, (double)C->RefLength);
    DFI_OUT_HLT->AddUnit("Velocity", UnitV, (double)C->RefVelocity);
    DFI_OUT_HLT->AddUnit("Pressure", UnitP, (double)C->BasePrs, DiffPrs, true);
  }
  
  
  // Divergence
  if ( C->varState[var_Div] == ON )
  {
    comp = 1;
    DFI_OUT_DIV = cdm_DFI::WriteInit(MPI_COMM_WORLD,
                                     cdm_DFI::Generate_DFI_Name(f_DivDebug),
                                     path,
                                     f_DivDebug,
                                     cdm_format,
                                     gc_out,
                                     datatype,
                                     comp,
                                     process,
                                     G_size,
                                     cdm_pit,
                                     cdm_org,
                                     cdm_div,
                                     head,
                                     cdm_tail,
                                     HostName,
                                     TimeSliceDir);
    
    if ( DFI_OUT_DIV == NULL )
    {
      Hostonly_ stamped_printf("\tFails to instance Divergence dfi.\n");
      Exit(0);
    }
    
    if( Slice == ON )
    {
      DFI_OUT_DIV->SetTimeSliceFlag(CDM::E_CDM_ON);
    }
    else
    {
      DFI_OUT_DIV->SetTimeSliceFlag(CDM::E_CDM_OFF);
    }
    
    DFI_OUT_DIV->AddUnit("Length"  , UnitL, (double)C->RefLength);
    DFI_OUT_DIV->AddUnit("Velocity", UnitV, (double)C->RefVelocity);
    DFI_OUT_DIV->AddUnit("Pressure", UnitP, (double)C->BasePrs, DiffPrs, true);
  }

}




// #################################################################
// 時間平均値のファイル出力
void SPH::OutputStatisticalVarables(const unsigned m_CurrentStep,
                                 const double m_CurrentTime,
                                 const unsigned m_CurrentStepStat,
                                 const double m_CurrentTimeStat,
                                 double& flop)
{
  REAL_TYPE f_min, f_max, min_tmp, max_tmp, vec_min[4], vec_max[4];
  REAL_TYPE minmax[2];
  REAL_TYPE cdm_minmax[8];
  
  // エラーコード
  CDM::E_CDM_ERRORCODE ret;
  
  
  // 出力ファイルの指定が有次元の場合
  double timeAvr;
  
  if (C->Unit.File == DIMENSIONAL)
  {
    timeAvr = m_CurrentTimeStat * C->Tscale;
  }
  else
  {
    timeAvr = m_CurrentTimeStat;
  }
  
  // 平均操作の母数
  unsigned stepAvr = m_CurrentStepStat;
  REAL_TYPE scale = 1.0;
  
  // ガイドセル出力
  int gc_out = C->GuideOut;
  
  // ファイル出力のタイムスタンプに使うステップ数
  unsigned m_step = m_CurrentStep;
  
  
  // ファイル出力のタイムスタンプの次元変換
  REAL_TYPE m_time;
  
  if (C->Unit.File == DIMENSIONAL)
  {
    m_time = (REAL_TYPE)(m_CurrentTime * C->Tscale);
  }
  else
  {
    m_time = (REAL_TYPE)m_CurrentTime;
  }
  
  
  if ( C->KindOfSolver != SOLID_CONDUCTION )
  {
    // Pressure
    if (C->Unit.File == DIMENSIONAL)
    {
      REAL_TYPE bp = ( C->Unit.Prs == Unit_Absolute ) ? C->BasePrs : 0.0;
      U.convArrayPrsND2D(d_ws, size, guide, d_ap, bp, C->RefDensity, C->RefVelocity, flop);
    }
    else
    {
      U.copyS3D(d_ws, size, guide, d_ap, scale);
    }
    
    // 最大値と最小値
    fb_minmax_s_ (&f_min, &f_max, size, &guide, d_ws, &flop);
    
    if ( numProc > 1 )
    {
      min_tmp = f_min;
      if ( paraMngr->Allreduce(&min_tmp, &f_min, 1, MPI_MIN) != CPM_SUCCESS ) Exit(0);
      
      max_tmp = f_max;
      if ( paraMngr->Allreduce(&max_tmp, &f_max, 1, MPI_MAX) != CPM_SUCCESS ) Exit(0);
    }
    
    minmax[0] = f_min;
    minmax[1] = f_max;
    
    
    if ( !DFI_OUT_PRSA )
    {
      printf("[%d] DFI_OUT_PRSA Pointer Error\n", paraMngr->GetMyRankID());
      Exit(-1);
    }
    
    DFI_OUT_PRSA->setVariableName(0, "AvrPressure");
  
    
    ret = DFI_OUT_PRSA->WriteData(m_step,   // 出力step番号
                                  m_time,   // 出力時刻
                                  size,     // d_wsの実ボクセル数
                                  1,        // d_wsの成分数
                                  guide,    // d_wsの仮想セル数
                                  d_ws,     // フィールドデータポインタ
                                  minmax,   // 最小値と最大値
                                  false,    // 平均出力指示 false:出力あり
                                  stepAvr,  // 平均をとったステップ数
                                  timeAvr); // 平均をとった時刻
    
    if ( ret != CDM::E_CDM_SUCCESS ) Exit(0);
    
    
    
    // Velocity
    REAL_TYPE unit_velocity = (C->Unit.File == DIMENSIONAL) ? C->RefVelocity : 1.0;
    
    fb_vout_nijk_(d_wv, d_av, size, &guide, RF->getV00(), &unit_velocity, &flop); // 配列並びを変換
    fb_minmax_vex_ (vec_min, vec_max, size, &guide, RF->getV00(), d_wv, &flop);

    
    
    if ( numProc > 1 )
    {
      REAL_TYPE vmin_tmp[4] = {vec_min[0], vec_min[1], vec_min[2], vec_min[3]};
      if ( paraMngr->Allreduce(vmin_tmp, vec_min, 4, MPI_MIN) != CPM_SUCCESS ) Exit(0);
      
      REAL_TYPE vmax_tmp[4] = {vec_max[0], vec_max[1], vec_max[2], vec_max[3]};
      if ( paraMngr->Allreduce(vmax_tmp, vec_max, 4, MPI_MAX) != CPM_SUCCESS ) Exit(0);
    }
    
    
    if ( !DFI_OUT_VELA )
    {
      printf("[%d] DFI_OUT_VELA Pointer Error\n", paraMngr->GetMyRankID());
      Exit(-1);
    }
    
    cdm_minmax[0] = vec_min[1]; ///<<< vec_u min
    cdm_minmax[1] = vec_max[1]; ///<<< vec_u max
    cdm_minmax[2] = vec_min[2]; ///<<< vec_v min
    cdm_minmax[3] = vec_max[2]; ///<<< vec_v max
    cdm_minmax[4] = vec_min[3]; ///<<< vec_w min
    cdm_minmax[5] = vec_max[3]; ///<<< vec_w max
    cdm_minmax[6] = vec_min[0]; ///<<< u,v,wの合成値のmin
    cdm_minmax[7] = vec_max[0]; ///<<< u,v,wの合成値のmax
    
    DFI_OUT_VELA->setVariableName(0, "Avr_U");
    DFI_OUT_VELA->setVariableName(1, "Avr_V");
    DFI_OUT_VELA->setVariableName(2, "Avr_W");
    
    ret = DFI_OUT_VELA->WriteData(m_step,
                                  m_time,
                                  size,
                                  3,
                                  guide,
                                  d_wv,
                                  cdm_minmax,
                                  false,
                                  stepAvr,
                                  timeAvr);
    
    if ( ret != CDM::E_CDM_SUCCESS ) Exit(0);
  }
  
  
  // Temperature
  if ( C->isHeatProblem() )
  {
    U.convArrayIE2Tmp(d_ws, size, guide, d_ae, d_bcd, mat_tbl, C->BaseTemp, C->DiffTemp, C->Unit.File, flop);
    
    fb_minmax_s_ (&f_min, &f_max, size, &guide, d_ws, &flop);
    
    if ( numProc > 1 )
    {
      min_tmp = f_min;
      if ( paraMngr->Allreduce(&min_tmp, &f_min, 1, MPI_MIN) != CPM_SUCCESS ) Exit(0);
      
      max_tmp = f_max;
      if ( paraMngr->Allreduce(&max_tmp, &f_max, 1, MPI_MAX) != CPM_SUCCESS ) Exit(0);
    }
    
    minmax[0] = f_min;
    minmax[1] = f_max;
    
    if ( !DFI_OUT_TEMPA )
    {
      printf("[%d] DFI_OUT_TEMPA Pointer Error\n", paraMngr->GetMyRankID());
      Exit(-1);
    }
    
    DFI_OUT_TEMPA->setVariableName(0, "Avr_Temp");
    
    ret = DFI_OUT_TEMPA->WriteData(m_step,
                                   m_time,
                                   size,
                                   1,
                                   guide,
                                   d_ws,
                                   minmax,
                                   false,
                                   stepAvr,
                                   timeAvr);
    
    if ( ret != CDM::E_CDM_SUCCESS ) Exit(0);
  }
}


// #################################################################
// 基本変数のファイル出力
void SPH::OutputBasicVariables(const unsigned m_CurrentStep,
                               const double m_CurrentTime,
                               double& flop)
{
  REAL_TYPE scale = 1.0;
  
  // ステップ数
  unsigned m_step = m_CurrentStep;
  
  // 時間の次元変換
  REAL_TYPE m_time;
  if (C->Unit.File == DIMENSIONAL)
  {
    m_time = (REAL_TYPE)(m_CurrentTime * C->Tscale);
  }
  else
  {
    m_time = (REAL_TYPE)m_CurrentTime;
  }
  
  // ガイドセル出力
  int gc_out = C->GuideOut;
  
  
  // 最大値と最小値
  REAL_TYPE f_min, f_max, min_tmp, max_tmp, vec_min[4], vec_max[4];
  REAL_TYPE minmax[2];
  REAL_TYPE cdm_minmax[8];
  
  
  // エラーコード
  CDM::E_CDM_ERRORCODE ret;
  
  // Velocity
  REAL_TYPE unit_velocity = (C->Unit.File == DIMENSIONAL) ? C->RefVelocity : 1.0;
  
  
  if ( C->KindOfSolver != SOLID_CONDUCTION )
  {
    // Pressure
    if (C->Unit.File == DIMENSIONAL)
    {
      REAL_TYPE bp = ( C->Unit.Prs == Unit_Absolute ) ? C->BasePrs : 0.0;
      U.convArrayPrsND2D(d_ws, size, guide, d_p, bp, C->RefDensity, C->RefVelocity, flop);
    }
    else
    {
      U.copyS3D(d_ws, size, guide, d_p, scale);
    }
    
    fb_minmax_s_ (&f_min, &f_max, size, &guide, d_ws, &flop);
    
    if ( numProc > 1 )
    {
      min_tmp = f_min;
      if( paraMngr->Allreduce(&min_tmp, &f_min, 1, MPI_MIN) != CPM_SUCCESS ) Exit(0);
      
      max_tmp = f_max;
      if( paraMngr->Allreduce(&max_tmp, &f_max, 1, MPI_MAX) != CPM_SUCCESS ) Exit(0);
    }
    minmax[0] = f_min;
    minmax[1] = f_max;
    
    if( !DFI_OUT_PRS )
    {
      printf("[%d] DFI_OUT_PRS Pointer Error\n", paraMngr->GetMyRankID());
      Exit(-1);
    }
    
    DFI_OUT_PRS->setVariableName(0, "Pressure");
    
    ret = DFI_OUT_PRS->WriteData(m_step,
                                 m_time,
                                 size,
                                 1,
                                 guide,
                                 d_ws,
                                 minmax,
                                 true,
                                 0,
                                 0.0);
    
    if( ret != CDM::E_CDM_SUCCESS )
    {
      Hostonly_ printf("CDMlib error code = %d\n", ret);
      Exit(0);
    }
    
    
    
    fb_vout_nijk_(d_wv, d_v, size, &guide, RF->getV00(), &unit_velocity, &flop);
    
    fb_minmax_vex_ (vec_min, vec_max, size, &guide, RF->getV00(), d_wv, &flop);
    
    
    if ( numProc > 1 )
    {
      REAL_TYPE vmin_tmp[4] = {vec_min[0], vec_min[1], vec_min[2], vec_min[3]};
      if( paraMngr->Allreduce(vmin_tmp, vec_min, 4, MPI_MIN) != CPM_SUCCESS ) Exit(0);
      
      REAL_TYPE vmax_tmp[4] = {vec_max[0], vec_max[1], vec_max[2], vec_max[3]};
      if( paraMngr->Allreduce(vmax_tmp, vec_max, 4, MPI_MAX) != CPM_SUCCESS ) Exit(0);
    }
    
    if ( !DFI_OUT_VEL )
    {
      printf("[%d] DFI_OUT_VEL Pointer Error\n", paraMngr->GetMyRankID());
      Exit(-1);
    }
    
    cdm_minmax[0] = vec_min[1]; ///<<< vec_u min
    cdm_minmax[1] = vec_max[1]; ///<<< vec_u max
    cdm_minmax[2] = vec_min[2]; ///<<< vec_v min
    cdm_minmax[3] = vec_max[2]; ///<<< vec_v max
    cdm_minmax[4] = vec_min[3]; ///<<< vec_w min
    cdm_minmax[5] = vec_max[3]; ///<<< vec_w max
    cdm_minmax[6] = vec_min[0]; ///<<< u,v,wの合成値のmin
    cdm_minmax[7] = vec_max[0]; ///<<< u,v,wの合成値のmax
    
    DFI_OUT_VEL->setVariableName(0, "u");
    DFI_OUT_VEL->setVariableName(1, "v");
    DFI_OUT_VEL->setVariableName(2, "w");
    
    ret = DFI_OUT_VEL->WriteData(m_step,
                                 m_time,
                                 size,
                                 3,
                                 guide,
                                 d_wv,
                                 cdm_minmax,
                                 true,
                                 0,
                                 0.0);
    
    if( ret != CDM::E_CDM_SUCCESS ) Exit(0);
    
    
    // Face Velocity

    fb_vout_nijk_(d_wv, d_vf, size, &guide, RF->getV00(), &unit_velocity, &flop);
    fb_minmax_vex_ (vec_min, vec_max, size, &guide, RF->getV00(), d_wv, &flop);

    
    if ( numProc > 1 )
    {
      REAL_TYPE vmin_tmp[4] = {vec_min[0], vec_min[1], vec_min[2], vec_min[3]};
      if( paraMngr->Allreduce(vmin_tmp, vec_min, 4, MPI_MIN) != CPM_SUCCESS ) Exit(0);
      
      REAL_TYPE vmax_tmp[4] = {vec_max[0], vec_max[1], vec_max[2], vec_max[3]};
      if( paraMngr->Allreduce(vmax_tmp, vec_max, 4, MPI_MAX) != CPM_SUCCESS ) Exit(0);
    }
    
    cdm_minmax[0] = vec_min[1]; ///<<< vec_u min
    cdm_minmax[1] = vec_max[1]; ///<<< vec_u max
    cdm_minmax[2] = vec_min[2]; ///<<< vec_v min
    cdm_minmax[3] = vec_max[2]; ///<<< vec_v max
    cdm_minmax[4] = vec_min[3]; ///<<< vec_w min
    cdm_minmax[5] = vec_max[3]; ///<<< vec_w max
    cdm_minmax[6] = vec_min[0]; ///<<< u,v,wの合成値のmin
    cdm_minmax[7] = vec_max[0]; ///<<< u,v,wの合成値のmax
    
    DFI_OUT_FVEL->setVariableName(0, "fu");
    DFI_OUT_FVEL->setVariableName(1, "fv");
    DFI_OUT_FVEL->setVariableName(2, "fw");
    
    ret = DFI_OUT_FVEL->WriteData(m_step,
                                  m_time,
                                  size,
                                  3,
                                  guide,
                                  d_wv,
                                  cdm_minmax,
                                  true,
                                  0,
                                  0.0);
    
    if ( ret != CDM::E_CDM_SUCCESS ) Exit(0);
  }
  
  
  // Tempearture
  if ( C->isHeatProblem() )
  {
    
    U.convArrayIE2Tmp(d_ws, size, guide, d_ie, d_bcd, mat_tbl, C->BaseTemp, C->DiffTemp, C->Unit.File, flop);
    
    fb_minmax_s_ (&f_min, &f_max, size, &guide, d_ws, &flop);
    
    if ( numProc > 1 )
    {
      min_tmp = f_min;
      if( paraMngr->Allreduce(&min_tmp, &f_min, 1, MPI_MIN) != CPM_SUCCESS ) Exit(0);
      
      max_tmp = f_max;
      if( paraMngr->Allreduce(&max_tmp, &f_max, 1, MPI_MAX) != CPM_SUCCESS ) Exit(0);
    }
    minmax[0] = f_min;
    minmax[1] = f_max;
    
    if ( !DFI_OUT_TEMP )
    {
      printf("[%d] DFI_OUT_TEMP Pointer Error\n", paraMngr->GetMyRankID());
      Exit(-1);
    }
    
    DFI_OUT_TEMP->setVariableName(0, "Temperature");
    
    ret = DFI_OUT_TEMP->WriteData(m_step,
                                  m_time,
                                  size,
                                  1,
                                  guide,
                                  d_ws,
                                  minmax,
                                  true,
                                  0,
                                  0.0);
    
    if ( ret != CDM::E_CDM_SUCCESS ) Exit(0);
  }
  

  
  
  
  // Total Pressure
  if (C->varState[var_TotalP] == ON )
  {
    fb_totalp_ (d_ws, size, &guide, d_v, d_p, RF->getV00(), &flop);
    
    // convert non-dimensional to dimensional, iff file is dimensional
    if (C->Unit.File == DIMENSIONAL)
    {
      U.convArrayTpND2D(d_ws, size, guide, C->RefDensity, C->RefVelocity);
    }
    
    
    fb_minmax_s_ (&f_min, &f_max, size, &guide, d_ws, &flop);
    
    if ( numProc > 1 )
    {
      min_tmp = f_min;
      if( paraMngr->Allreduce(&min_tmp, &f_min, 1, MPI_MIN) != CPM_SUCCESS ) Exit(0);
      
      max_tmp = f_max;
      if( paraMngr->Allreduce(&max_tmp, &f_max, 1, MPI_MAX) != CPM_SUCCESS ) Exit(0);
    }
    minmax[0] = f_min;
    minmax[1] = f_max;
    
    if ( !DFI_OUT_TP )
    {
      printf("[%d] DFI_OUT_TP Pointer Error\n",paraMngr->GetMyRankID());
      Exit(-1);
    }
    
    DFI_OUT_TP->setVariableName(0, "TotalPressure");
    
    ret = DFI_OUT_TP->WriteData(m_step,
                                m_time,
                                size,
                                1,
                                guide,
                                d_ws,
                                minmax,
                                true,
                                0,
                                0.0);
    
    if ( ret != CDM::E_CDM_SUCCESS ) Exit(0);
  }
  
  
  // Vorticity
  if (C->varState[var_Vorticity] == ON )
  {
    rot_v_(d_wv, size, &guide, pitch, d_v, d_cdf, RF->getV00(), &flop);
    
    REAL_TYPE  vz[3];
    vz[0] = vz[1] = vz[2] = 0.0;
    unit_velocity = (C->Unit.File == DIMENSIONAL) ? C->RefVelocity/C->RefLength : 1.0;
    
    fb_vout_nijk_(d_iobuf, d_wv, size, &guide, vz, &unit_velocity, &flop);
    fb_minmax_vex_ (vec_min, vec_max, size, &guide, RF->getV00(), d_iobuf, &flop);
    
    
    if ( numProc > 1 )
    {
      REAL_TYPE vmin_tmp[4] = {vec_min[0], vec_min[1], vec_min[2], vec_min[3]};
      if( paraMngr->Allreduce(vmin_tmp, vec_min, 4, MPI_MIN) != CPM_SUCCESS ) Exit(0);
      
      REAL_TYPE vmax_tmp[4] = {vec_max[0], vec_max[1], vec_max[2], vec_max[3]};
      if( paraMngr->Allreduce(vmax_tmp, vec_max, 4, MPI_MAX) != CPM_SUCCESS ) Exit(0);
    }
    
    if ( !DFI_OUT_VRT )
    {
      printf("[%d] DFI_OUT_VRT Pointer Error\n", paraMngr->GetMyRankID());
      Exit(-1);
    }
    
    cdm_minmax[0] = vec_min[1]; ///<<< vec_u min
    cdm_minmax[1] = vec_max[1]; ///<<< vec_u max
    cdm_minmax[2] = vec_min[2]; ///<<< vec_v min
    cdm_minmax[3] = vec_max[2]; ///<<< vec_v max
    cdm_minmax[4] = vec_min[3]; ///<<< vec_w min
    cdm_minmax[5] = vec_max[3]; ///<<< vec_w max
    cdm_minmax[6] = vec_min[0]; ///<<< u,v,wの合成値のmin
    cdm_minmax[7] = vec_max[0]; ///<<< u,v,wの合成値のmax
    
    DFI_OUT_VRT->setVariableName(0, "vrt_u");
    DFI_OUT_VRT->setVariableName(1, "vrt_v");
    DFI_OUT_VRT->setVariableName(2, "vrt_w");
    
    ret = DFI_OUT_VRT->WriteData(m_step,
                                 m_time,
                                 size,
                                 3,
                                 guide,
                                 d_iobuf,
                                 cdm_minmax,
                                 true,
                                 0,
                                 0.0);
    
    if ( ret != CDM::E_CDM_SUCCESS ) Exit(0);
  }
  
  
  // 2nd Invariant of Velocity Gradient Tensor
  if (C->varState[var_Qcr] == ON )
  {
    i2vgt_ (d_iobuf, size, &guide, pitch, d_v, d_cdf, RF->getV00(), &flop);
    
    // 無次元で出力
    U.copyS3D(d_ws, size, guide, d_iobuf, scale);
    
    fb_minmax_s_ (&f_min, &f_max, size, &guide, d_ws, &flop);
    
    if ( numProc > 1 )
    {
      min_tmp = f_min;
      if( paraMngr->Allreduce(&min_tmp, &f_min, 1, MPI_MIN) != CPM_SUCCESS ) Exit(0);
      
      max_tmp = f_max;
      if( paraMngr->Allreduce(&max_tmp, &f_max, 1, MPI_MAX) != CPM_SUCCESS ) Exit(0);
    }
    minmax[0] = f_min;
    minmax[1] = f_max;
    
    if ( !DFI_OUT_I2VGT )
    {
      printf("[%d] DFI_OUT_I2VGT Pointer Error\n", paraMngr->GetMyRankID());
      Exit(-1);
    }
    
    DFI_OUT_I2VGT->setVariableName(0, "Qcriterion");
    
    ret = DFI_OUT_I2VGT->WriteData(m_step,
                                   m_time,
                                   size,
                                   1,
                                   guide,
                                   d_ws,
                                   minmax,
                                   true,
                                   0,
                                   0.0);
    
    if( ret != CDM::E_CDM_SUCCESS ) Exit(0);
  }
  
  
  // Helicity
  if (C->varState[var_Helicity] == ON )
  {
    helicity_(d_iobuf, size, &guide, pitch, d_v, d_cdf, RF->getV00(), &flop);
    
    // 無次元で出力
    U.copyS3D(d_ws, size, guide, d_iobuf, scale);
    
    fb_minmax_s_ (&f_min, &f_max, size, &guide, d_ws, &flop);
    
    if ( numProc > 1 )
    {
      min_tmp = f_min;
      if( paraMngr->Allreduce(&min_tmp, &f_min, 1, MPI_MIN) != CPM_SUCCESS ) Exit(0);
      
      max_tmp = f_max;
      if( paraMngr->Allreduce(&max_tmp, &f_max, 1, MPI_MAX) != CPM_SUCCESS ) Exit(0);
    }
    minmax[0] = f_min;
    minmax[1] = f_max;
    
    if ( !DFI_OUT_HLT )
    {
      printf("[%d] DFI_OUT_HLT Pointer Error\n", paraMngr->GetMyRankID());
      Exit(-1);
    }
    
    DFI_OUT_HLT->setVariableName(0, "Helicity");
    
    ret = DFI_OUT_HLT->WriteData(m_step,
                                 m_time,
                                 size,
                                 1,
                                 guide,
                                 d_ws,
                                 minmax,
                                 true,
                                 0,
                                 0.0);
    
    if( ret != CDM::E_CDM_SUCCESS ) Exit(0);
  }
  
  
  // Divergence for Debug
  if (C->varState[var_Div] == ON )
  {
    U.cnv_Div(d_ws, d_dv, size, guide);
    
    fb_minmax_s_ (&f_min, &f_max, size, &guide, d_ws, &flop);
    
    if ( numProc > 1 )
    {
      min_tmp = f_min;
      if( paraMngr->Allreduce(&min_tmp, &f_min, 1, MPI_MIN) != CPM_SUCCESS ) Exit(0);
      
      max_tmp = f_max;
      if( paraMngr->Allreduce(&max_tmp, &f_max, 1, MPI_MAX) != CPM_SUCCESS ) Exit(0);
    }
    minmax[0] = f_min;
    minmax[1] = f_max;
    
    if( !DFI_OUT_DIV )
    {
      printf("[%d] DFI_OUT_DIV Pointer Error\n", paraMngr->GetMyRankID());
      Exit(-1);
    }
    
    DFI_OUT_DIV->setVariableName(0, "Divergence");
    
    ret = DFI_OUT_DIV->WriteData(m_step,
                                 m_time,
                                 size,
                                 1,
                                 guide,
                                 d_ws,
                                 minmax,
                                 true,
                                 0,
                                 0.0);
    
    if( ret != CDM::E_CDM_SUCCESS ) Exit(0);
  }

}



// #################################################################
// リスタート時の平均値ファイル読み込み
void SPH::RestartStatistic(FILE* fp,
                           const unsigned m_CurrentStep,
                           const double m_CurrentTime,
                           unsigned& m_CurrentStepStat,
                           double& m_CurrentTimeStat,
                           double& flop)
{
  std::string fname;
  std::string fmt(file_fmt_ext);
  
  unsigned m_Session_step = C->Interval[Control::tg_compute].getStartStep(); ///< セッションの開始ステップ
  double   m_Session_time = C->Interval[Control::tg_compute].getStartTime(); ///< セッションの開始時刻
  
  // リスタートステップ
  unsigned m_RestartStep;
  if ( C->Interval[Control::tg_compute].getMode() == IntervalManager::By_step )
  {
    m_RestartStep = C->Interval[Control::tg_compute].getStartStep();
  }
  else // By_time
  {
    m_RestartStep = C->Interval[Control::tg_compute].restartStep;
  }
  
  
  // ガイド出力
  int gs = C->GuideOut;
  
  
  // まだ平均値開始時刻になっていなければ，何もしない
  if ( C->Interval[Control::tg_statistic].getMode() == IntervalManager::By_step )
  {
    if ( m_Session_step >= C->Interval[Control::tg_statistic].getStartStep() )
    {
      Hostonly_ printf     ("\tRestart from Previous Calculation Results of Statistical field\n");
      Hostonly_ fprintf(fp, "\tRestart from Previous Calculation Results of Statistical field\n");
      Hostonly_ printf     ("\tStep : base=%u current=%u\n", m_Session_step, m_CurrentStep);
      Hostonly_ fprintf(fp, "\tStep : base=%u current=%u\n", m_Session_step, m_CurrentStep);
    }
    else
    {
      return;
    }
  }
  else if ( C->Interval[Control::tg_statistic].getMode() == IntervalManager::By_time )
  {
    if ( m_Session_time >= C->Interval[Control::tg_statistic].getStartTime() )
    {
      Hostonly_ printf     ("\tRestart from Previous Calculation Results of Statistical field\n");
      Hostonly_ fprintf(fp, "\tRestart from Previous Calculation Results of Statistical field\n");
      Hostonly_ printf     ("\tTime : base=%e[sec.]/%e[-] current=%e[-]\n", m_Session_time*C->Tscale, m_Session_time, m_CurrentTime);
      Hostonly_ fprintf(fp, "\tTime : base=%e[sec.]/%e[-] current=%e[-]\n", m_Session_time*C->Tscale, m_Session_time, m_CurrentTime);
    }
    else
    {
      return;
    }
  }
  else
  {
    Exit(0);
  }
  
  
  // 現在のセッションの領域分割数の取得
  int gdiv[3] = {1, 1, 1};
  
  if ( numProc > 1)
  {
    const int* m_div = paraMngr->GetDivNum();
    for (int i=0; i<3; i++ ) gdiv[i]=m_div[i];
  }
  
  
  // エラーコード
  CDM::E_CDM_ERRORCODE cdm_error;
  
  
  // Statistical dataの初期化
  if ( C->Mode.Statistic == ON && C->Interval[Control::tg_statistic].isStarted(m_CurrentStep, m_CurrentTime) )
  {
    DFI_IN_PRSA = cdm_DFI::ReadInit(MPI_COMM_WORLD, f_dfi_in_prsa, G_size, gdiv, cdm_error);
    if ( cdm_error != CDM::E_CDM_SUCCESS ) Exit(0);
    
    DFI_IN_VELA = cdm_DFI::ReadInit(MPI_COMM_WORLD, f_dfi_in_vela, G_size, gdiv, cdm_error);
    if ( cdm_error != CDM::E_CDM_SUCCESS ) Exit(0);
    
    if ( DFI_IN_PRSA == NULL || DFI_IN_VELA == NULL ) Exit(0);
    
    
    if ( C->isHeatProblem() )
    {
      DFI_IN_TEMPA = cdm_DFI::ReadInit(MPI_COMM_WORLD, f_dfi_in_tempa, G_size, gdiv, cdm_error);
      if ( cdm_error != CDM::E_CDM_SUCCESS ) Exit(0);
      if ( DFI_IN_TEMPA == NULL ) Exit(0);
    }
  }
  
  
  
  
  unsigned step_stat = 0;
  double time_stat = 0.0;
  
  
  // Pressure
  REAL_TYPE bp = ( C->Unit.Prs == Unit_Absolute ) ? C->BasePrs : 0.0;
  
  
  //自身の領域終点インデックス
  int tail[3];
  for (int i=0; i<3; i++) tail[i] = head[i]+size[i]-1;
  
  
  double r_time;
  if ( DFI_IN_PRSA->ReadData(d_ap,
                             m_RestartStep,
                             guide,
                             G_size,
                             gdiv,
                             head,
                             tail,
                             r_time,
                             false,
                             step_stat,
                             time_stat) != CDM::E_CDM_SUCCESS ) Exit(0);
  
  if( d_ap == NULL ) Exit(0);
  
  
  m_CurrentStepStat = step_stat;
  m_CurrentTimeStat = time_stat;
  
  if ( DFI_IN_VELA->ReadData(d_wv,
                             m_RestartStep,
                             guide,
                             G_size,
                             gdiv,
                             head,
                             tail,
                             r_time,
                             false,
                             step_stat,
                             time_stat) != CDM::E_CDM_SUCCESS ) Exit(0);
  
  if( d_wv == NULL ) Exit(0);
  
  REAL_TYPE refv = (C->Unit.File == DIMENSIONAL) ? C->RefVelocity : 1.0;
  REAL_TYPE scale = (REAL_TYPE)step_stat;
  REAL_TYPE u0[4];
  
  RF->copyV00(u0);
  
  
  fb_vin_nijk_(d_av, size, &guide, d_wv, u0, &refv, &flop);
  
  if ( (step_stat != m_CurrentStepStat) || (time_stat != m_CurrentTimeStat) ) // 圧力とちがう場合
  {
    Hostonly_ printf     ("\n\tTime stamp is different between files\n");
    Hostonly_ fprintf(fp, "\n\tTime stamp is different between files\n");
    Exit(0);
  }
  
  
  // Temperature
  if ( C->isHeatProblem() )
  {
    if ( DFI_IN_TEMPA->ReadData(d_ae,
                                m_RestartStep,
                                guide,
                                G_size,
                                gdiv,
                                head,
                                tail,
                                r_time,
                                false,
                                step_stat,
                                time_stat) != CDM::E_CDM_SUCCESS ) Exit(0);
    
    if ( d_ae == NULL ) Exit(0);
    
    if ( (step_stat != m_CurrentStepStat) || (time_stat != m_CurrentTimeStat) )
    {
      Hostonly_ printf     ("\n\tTime stamp is different between files\n");
      Hostonly_ fprintf(fp, "\n\tTime stamp is different between files\n");
      Exit(0);
    }
  }
}



// #################################################################
// リスタート時の瞬時値ファイル読み込み
void SPH::RestartInstantaneous(FILE* fp,
                               unsigned& m_CurrentStep,
                               double& m_CurrentTime,
                               double& flop)
{
  double time, r_time;
  std::string fname;
  std::string fmt(file_fmt_ext);
  
  REAL_TYPE bp = ( C->Unit.Prs == Unit_Absolute ) ? C->BasePrs : 0.0;
  REAL_TYPE refD = C->RefDensity;
  REAL_TYPE refV = C->RefVelocity;
  
  
  // リスタートステップ
  unsigned m_RestartStep;
  if ( C->Interval[Control::tg_compute].getMode() == IntervalManager::By_step )
  {
    m_RestartStep = C->Interval[Control::tg_compute].getStartStep();
  }
  else // By_time
  {
    m_RestartStep = C->Interval[Control::tg_compute].restartStep;
  }
  
  
  // ガイド出力
  int gs = C->GuideOut;
  
  // dummy
  unsigned i_dummy=0;
  double f_dummy=0.0;
  
  
  const int* m_div = paraMngr->GetDivNum();
  
  // 自身の領域終点インデクス
  int tail[3];
  for (int i=0;i<3;i++) tail[i]=head[i]+size[i]-1;
  
  // Pressure
  if ( DFI_IN_PRS->ReadData(d_p,
                            m_RestartStep,
                            guide,
                            G_size,
                            (int *)m_div,
                            head,
                            tail,
                            r_time,
                            true,
                            i_dummy,
                            f_dummy) != CDM::E_CDM_SUCCESS ) Exit(0);
  
  if ( d_p == NULL ) Exit(0);
  time = r_time;
  
  // 有次元の場合，無次元に変換する
  if ( C->Unit.File == DIMENSIONAL )
  {
    REAL_TYPE bp = ( C->Unit.Prs == Unit_Absolute ) ? C->BasePrs : 0.0;
    U.convArrayPrsD2ND(d_p, size, guide, bp, C->RefDensity, C->RefVelocity, flop);
  }
  
  Hostonly_ printf     ("\tPressure has read :\tstep=%d  time=%e [%s]\n",
                        m_RestartStep, time, (C->Unit.File == DIMENSIONAL)?"sec.":"-");
  Hostonly_ fprintf(fp, "\tPressure has read :\tstep=%d  time=%e [%s]\n",
                    m_RestartStep, time, (C->Unit.File == DIMENSIONAL)?"sec.":"-");
  
  
  // ここでタイムスタンプを得る
  if (C->Unit.File == DIMENSIONAL) time /= C->Tscale;
  m_CurrentStep = m_RestartStep;
  m_CurrentTime = time;
  
  // v00[]に値をセット
  RF->setV00(time);
  
  
  if ( DFI_IN_VEL->ReadData(d_wv,
                            m_RestartStep,
                            guide,
                            G_size,
                            (int *)m_div,
                            head,
                            tail,
                            r_time,
                            true,
                            i_dummy,
                            f_dummy) != CDM::E_CDM_SUCCESS ) Exit(0);
  
  if( d_wv == NULL ) Exit(0);
  
  REAL_TYPE refv = (C->Unit.File == DIMENSIONAL) ? refV : 1.0;
  REAL_TYPE u0[4];
  RF->copyV00(u0);
  
  
  Hostonly_ printf     ("\tVelocity has read :\tstep=%d  time=%e [%s]\n",
                        m_RestartStep, time, (C->Unit.File == DIMENSIONAL)?"sec.":"-");
  Hostonly_ fprintf(fp, "\tVelocity has read :\tstep=%d  time=%e [%s]\n",
                    m_RestartStep, time, (C->Unit.File == DIMENSIONAL)?"sec.":"-");
  
  time = r_time;
  
  if (C->Unit.File == DIMENSIONAL) time /= C->Tscale;
  
  if ( time != m_CurrentTime )
  {
    Hostonly_ printf     ("\n\tTime stamp is different between files\n");
    Hostonly_ fprintf(fp, "\n\tTime stamp is different between files\n");
    Exit(0);
  }
  
  // indexの変換と無次元化
  fb_vin_nijk_(d_v, size, &guide, d_wv, u0, &refv, &flop);
  
  
  
  if ( !C->isHeatProblem() ) return;
  
  
  
  // Instantaneous Temperature fields
  if ( DFI_IN_TEMP->ReadData(d_ws,
                             m_RestartStep,
                             guide,
                             G_size,
                             (int *)m_div,
                             head,
                             tail,
                             r_time,
                             true,
                             i_dummy,
                             f_dummy) != CDM::E_CDM_SUCCESS ) Exit(0);
  
  if( d_ws == NULL ) Exit(0);
  
  time = r_time;
  
  if (C->Unit.File == DIMENSIONAL) time /= C->Tscale;
  
  Hostonly_ printf     ("\tTemperature has read :\tstep=%d  time=%e [%s]\n",
                        m_RestartStep, time, (C->Unit.File == DIMENSIONAL)?"sec.":"-");
  Hostonly_ fprintf(fp, "\tTemperature has read :\tstep=%d  time=%e [%s]\n",
                    m_RestartStep, time, (C->Unit.File == DIMENSIONAL)?"sec.":"-");
  
  if ( time != m_CurrentTime )
  {
    Hostonly_ printf     ("\n\tTime stamp is different between files\n");
    Hostonly_ fprintf(fp, "\n\tTime stamp is different between files\n");
    Exit(0);
  }
  
  U.convArrayTmp2IE(d_ie, size, guide, d_ws, d_bcd, mat_tbl, C->BaseTemp, C->DiffTemp, C->Unit.File, flop);
}




// #################################################################
// リスタートモードを判定し，瞬時値のデータをロードする
void SPH::Restart(FILE* fp, unsigned& m_CurrentStep, double& m_CurrentTime)
{
  if ( C->Start != initial_start)
  {
    // エラーコード
    CDM::E_CDM_ERRORCODE cdm_error;
    
    
    // 現在のセッションの領域分割数の取得
    int gdiv[3] = {1, 1, 1};
    
    if ( numProc > 1)
    {
      const int* p_div = paraMngr->GetDivNum();
      for (int i=0; i<3; i++ ) gdiv[i]=p_div[i];
    }
    
    
    // Instantaneous dataの初期化
    
    // Pressure
    DFI_IN_PRS = cdm_DFI::ReadInit(MPI_COMM_WORLD, f_dfi_in_prs, G_size, gdiv, cdm_error);
    if ( cdm_error != CDM::E_CDM_SUCCESS ) Exit(0);
    
    
    // Velocity
    DFI_IN_VEL = cdm_DFI::ReadInit(MPI_COMM_WORLD, f_dfi_in_vel, G_size, gdiv, cdm_error);
    if ( cdm_error != CDM::E_CDM_SUCCESS ) Exit(0);
    
    if ( DFI_IN_PRS == NULL || DFI_IN_VEL == NULL ) Exit(0);
    
    
    // Fvelocity
    DFI_IN_FVEL = cdm_DFI::ReadInit(MPI_COMM_WORLD, f_dfi_in_fvel, G_size, gdiv, cdm_error);
    if ( cdm_error != CDM::E_CDM_SUCCESS ) Exit(0);
    if ( DFI_IN_FVEL == NULL ) Exit(0);
    
    // Temperature
    if ( C->isHeatProblem() )
    {
      DFI_IN_TEMP = cdm_DFI::ReadInit(MPI_COMM_WORLD, f_dfi_in_temp, G_size, gdiv, cdm_error);
      if ( cdm_error != CDM::E_CDM_SUCCESS ) Exit(0);
      if ( DFI_IN_TEMP == NULL ) Exit(0);
    }
    
    
    ///CDM FileInfoの成分名の取得、成分名が登録されていないときは、空白が戻される
    /*
     std::string VCompVariable[3];
     VCompVariable[0]=DFI_IN_VEL->getComponentVariable(0);
     VCompVariable[1]=DFI_IN_VEL->getComponentVariable(1);
     VCompVariable[2]=DFI_IN_VEL->getComponentVariable(2);
     
     ///CDM TimeSliceのVectorMinMaxの取得、取得出来たときはCDM::E_CDM_SUCCESS
     double vec_minmax[2];
     cdm_error = DFI_IN_VEL->getVectorMinMax(C->Restart_step,vec_minmax[0],vec_minmax[1]);
     
     ///CDM TimeSlice minmaxの取得、取得出来たときはCDM::E_CDM_SUCCESS
     double minmax[6];
     cdm_error =  DFI_IN_VEL->getMinMax(C->Restart_step,0,minmax[0],minmax[1]);
     cdm_error =  DFI_IN_VEL->getMinMax(C->Restart_step,1,minmax[2],minmax[3]);
     cdm_error =  DFI_IN_VEL->getMinMax(C->Restart_step,2,minmax[4],minmax[5]);
     */
    
    
    
    /* Statistical dataの初期化
     if ( C->Mode.Statistic == ON && C->Interval[Control::tg_statistic].isStarted(CurrentStep, CurrentTime) )
     {
     DFI_IN_PRSA = cdm_DFI::ReadInit(MPI_COMM_WORLD, C->f_dfi_in_prsa, G_size, gdiv, cdm_error);
     if ( cdm_error != CDM::E_CDM_SUCCESS ) Exit(0);
     
     DFI_IN_VELA = cdm_DFI::ReadInit(MPI_COMM_WORLD, C->f_dfi_in_vela, G_size, gdiv, cdm_error);
     if ( cdm_error != CDM::E_CDM_SUCCESS ) Exit(0);
     
     if ( DFI_IN_PRSA == NULL || DFI_IN_VELA == NULL ) Exit(0);
     
     
     if ( C->isHeatProblem() )
     {
     DFI_IN_TEMPA = cdm_DFI::ReadInit(MPI_COMM_WORLD, C->f_dfi_in_tempa, G_size, gdiv, cdm_error);
     if ( cdm_error != CDM::E_CDM_SUCCESS ) Exit(0);
     if ( DFI_IN_TEMPA == NULL ) Exit(0);
     }
     }
     */
    
    
    bool isSameDiv = true; // 同一分割数
    bool isSameRes = true; // 同一解像度
    
    
    // 前のセッションの領域分割数の取得
    const int* DFI_div=NULL;
    
    if ( C->KindOfSolver != SOLID_CONDUCTION )
    {
      DFI_div = DFI_IN_PRS->GetDFIGlobalDivision();
    }
    else
    {
      DFI_div = DFI_IN_TEMP->GetDFIGlobalDivision();
    }
    
    
    
    // 前セッションと領域分割数が異なる場合
    for (int i=0; i<3; i++ )
    {
      if ( gdiv[i] != DFI_div[i] )
      {
        isSameDiv = false;
      }
    }
    
    // 前のセッションの全要素数の取得
    const int* DFI_G_size = DFI_IN_PRS->GetDFIGlobalVoxel();
    
    
    // 前セッションと全要素数が異なる場合
    for (int i=0; i<3; i++ )
    {
      if ( G_size[i] != DFI_G_size[i] )
      {
        isSameRes = false;
      }
    }
    
    //  ボクセル数が2倍のチェック
    if ( !isSameRes )
    {
      for(int i=0; i<3; i++)
      {
        if ( G_size[i] != DFI_G_size[i]*2 )
        {
          printf("\tDimension size error (%d %d %d)\n", G_size[0], G_size[1], G_size[2]);
          Exit(0);
        }
      }
    }
    
    
    // モード判定と登録
    if ( isSameDiv )
    {
      if ( isSameRes ) // 同一解像度、同一分割数
      {
        C->Start = restart_sameDiv_sameRes;
      }
      else // Refinement、同一分割数
      {
        C->Start = restart_sameDiv_refinement;
      }
    }
    else
    {
      if ( isSameRes ) // 同一解像度、異なる分割数
      {
        C->Start = restart_diffDiv_sameRes;
      }
      else // Refinement、異なる分割数
      {
        C->Start = restart_diffDiv_refinement;
      }
    }

  }
  
  //
  //
  // 初期スタートのステップ，時間を設定する
  if ( C->Start == initial_start || C->Hide.PM_Test == ON  )
  {
    m_CurrentStep = 0;
    m_CurrentTime = 0.0;
    
    // V00の値のセット．モードがONの場合はV00[0]=1.0に設定，そうでなければtmに応じた値
    if ( C->CheckParam == ON ) RF->setV00(m_CurrentTime, true);
    else                       RF->setV00(m_CurrentTime);
    
    return;
  }
  
  
  
  switch (C->Start)
  {
      
      // 同一解像度・同一分割数のリスタート
    case restart_sameDiv_sameRes:
      Hostonly_ fprintf(stdout, "\t>> Restart with same resolution and same num. of division\n\n");
      Hostonly_ fprintf(fp, "\t>> Restart with same resolution and same num. of division\n\n");
      break;
      
    case restart_sameDiv_refinement: // 同一分割数・リファインメント
      Hostonly_ fprintf(stdout, "\t>> Restart with refinemnt and same num. of division\n\n");
      Hostonly_ fprintf(fp, "\t>> Restart with refinemnt and same num. of division\n\n");
      break;
      
    case restart_diffDiv_sameRes:    // 異なる分割数・同一解像度
      Hostonly_ fprintf(stdout, "\t>> Restart with same resolution and different division\n\n");
      Hostonly_ fprintf(fp, "\t>> Restart with same resolution and different division\n\n");
      break;
      
    case restart_diffDiv_refinement: // 異なる分割数・リファインメント
      Hostonly_ fprintf(stdout, "\t>> Restart with refinement and different division\n\n");
      Hostonly_ fprintf(fp, "\t>> Restart with refinement and different division\n\n");
      break;
      
    default:
      Exit(0);
      break;
  }
  
  double flop_task = 0.0;
  RestartInstantaneous(fp, m_CurrentStep, m_CurrentTime, flop_task);
}





// #################################################################
// チャネル乱流統計量の出力
void SPH::OutputMean(REAL_TYPE*         d_av,
                     REAL_TYPE*         d_rms_mean_v,
                     REAL_TYPE*         d_aR,
                     REAL_TYPE*         d_aP,
                     REAL_TYPE*         d_aE,
                     REAL_TYPE*         d_aT,
                     REAL_TYPE*         d_aPI,
                     int                myRank,
                     int*               sz,
                     unsigned long int  CurrentStepStat,
                     REAL_TYPE*         dh,
                     int*               g,
                     double&            flop)
{
   int i;
   FILE *fp1, *fp2;
   char str1[512], str2[512];

   // 主流方向・スパン方向で平均化
   averaging_xz_plane_(d_av_mean, d_arms_mean, d_aR_mean, d_aP_mean, d_aE_mean, d_aT_mean, d_aPI_mean, sz, g, d_av, d_rms_mean_v, d_aR, d_aP, d_aE, d_aT, d_aPI);

   sprintf(str1, "channel_base_p%d_latest.txt", myRank);
   fp1 = fopen(str1, "w");
   fprintf(fp1, "# step, y, umean, vmean, wmean, urms, vrms, wrms, uvmean \n");
   for(i = 0; i < sz[1]; i++)
   {
       fprintf(
                  fp1, "%lu %f %e %e %e %e %e %e %e \n", 
                  // 積算ステップ
                  CurrentStepStat,
                  // 壁面鉛直座標
                  (i+0.5)*dh[1],
                  // 平均速度
                  d_av_mean[i*3],   d_av_mean[i*3+1],   d_av_mean[i*3+2], 
                  // 乱流強度
                  d_arms_mean[i*3], d_arms_mean[i*3+1], d_arms_mean[i*3+2],
                  // レイノルズせん断応力
                  d_aR_mean[i*6+1]
             );
   }
   fclose(fp1);

   sprintf(str1, "channel_reynolds_p%d_latest.txt", myRank);
   fp1 = fopen(str1, "w");
   fprintf(fp1, "# step, y, uu, uv, uw, vv, vw, ww \n");
   for(i = 0; i < sz[1]; i++)
   {
       fprintf(
                  fp1, "%lu %f %e %e %e %e %e %e \n", 
                  // 積算ステップ
                  CurrentStepStat,
                  // 壁面鉛直座標
                  (i+0.5)*dh[1],
                  // レイノルズ応力
                  d_aR_mean[i*6],   d_aR_mean[i*6+1], d_aR_mean[i*6+2], 
                  d_aR_mean[i*6+3], d_aR_mean[i*6+4], d_aR_mean[i*6+5]
             );
   }
   fclose(fp1);

   sprintf(str1, "channel_budget_uu_p%d_latest.txt", myRank);
   fp1 = fopen(str1, "w");
   fprintf(fp1, "# step, y, Prod, Diss, Turb, PressGrad \n");
   for(i = 0; i < sz[1]; i++)
   {
       fprintf(
                  fp1, "%lu %f %e %e %e %e \n", 
                  // 積算ステップ
                  CurrentStepStat,
                  // 壁面鉛直座標
                  (i+0.5)*dh[1],
                  // レイノルズ応力収支
                  d_aP_mean[i*6], d_aE_mean[i*6], d_aT_mean[i*6], d_aPI_mean[i*6]
             );
   }
   fclose(fp1);


   sprintf(str2, "channel_base_p%d_step%lu.txt", myRank, CurrentStepStat);
   fp2 = fopen(str2, "w");
   fprintf(fp2, "# step, y, umean, vmean, wmean, urms, vrms, wrms, uvmean \n");
   for(i = 0; i < sz[1]; i++)
   {
       fprintf(
                  fp2, "%lu %f %e %e %e %e %e %e %e \n", 
                  // 積算ステップ
                  CurrentStepStat,
                  // 壁面鉛直座標
                  (i+0.5)*dh[1],
                  // 平均速度
                  d_av_mean[i*3],   d_av_mean[i*3+1],   d_av_mean[i*3+2], 
                  // 乱流強度
                  d_arms_mean[i*3], d_arms_mean[i*3+1], d_arms_mean[i*3+2],
                  // レイノルズせん断応力
                  d_aR_mean[i*6+1]
             );
   }
   fclose(fp2);

   sprintf(str2, "channel_reynolds_p%d_step%lu.txt", myRank, CurrentStepStat);
   fp2 = fopen(str2, "w");
   fprintf(fp2, "# step, y, uu, uv, uw, vv, vw, ww \n");
   for(i = 0; i < sz[1]; i++)
   {
       fprintf(
                  fp2, "%lu %f %e %e %e %e %e %e \n", 
                  // 積算ステップ
                  CurrentStepStat,
                  // 壁面鉛直座標
                  (i+0.5)*dh[1],
                  // レイノルズ応力
                  d_aR_mean[i*6],   d_aR_mean[i*6+1], d_aR_mean[i*6+2], 
                  d_aR_mean[i*6+3], d_aR_mean[i*6+4], d_aR_mean[i*6+5]
             );
   }
   fclose(fp2);

   sprintf(str2, "channel_budget_uu_p%d_step%lu.txt", myRank, CurrentStepStat);
   fp2 = fopen(str2, "w");
   fprintf(fp2, "# step, y, Prod, Diss, Turb, PressGrad \n");
   for(i = 0; i < sz[1]; i++)
   {
       fprintf(
                  fp2, "%lu %f %e %e %e %e \n", 
                  // 積算ステップ
                  CurrentStepStat,
                  // 壁面鉛直座標
                  (i+0.5)*dh[1],
                  // レイノルズ応力収支 (uu)
                  d_aP_mean[i*6], d_aE_mean[i*6], d_aT_mean[i*6], d_aPI_mean[i*6]
             );
   }
   fclose(fp2);

}

