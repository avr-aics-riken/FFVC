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
 * @file   ffv_plot3d.C
 * @brief  File IO of PLOT3D Class
 * @author aics
 */

#include "ffv_plot3d.h"

#include "ffv_LSfunc.h"
#include "ffv_Ffunc.h"

// #################################################################
// リスタートのDFIファイル
void PLT3D::getRestartDFI()
{
  string str;
  string label;
  
  
  // リスタート時のDFIファイル名
  if ( C->Start != initial_start )
  {
    label="/StartCondition/Restart/DFIfiles/Instantaneous";
    
    if ( tpCntl->getInspectedValue(label, str ) )
    {
      f_dfi_in_ins = str.c_str();
    }
    if ( f_dfi_in_ins.empty() == true ) f_dfi_in_ins = "field.dfi";
    
    
    // リスタート時のガイドセル数取得 >> field.dfiのパラメータを取得
    TextParser* tp_index = new TextParser;
    
    int ierror=0;
    
    if ( (ierror = tp_index->read(f_dfi_in_ins)) != TP_NO_ERROR )
    {
      Hostonly_ stamped_printf("\tError at reading '%s' file : %d\n", f_dfi_in_ins.c_str(), ierror);
      Exit(0);
    }
    
    int ct;
    
    label = "/FileInfo/GuideCell";
    if ( !tp_index->getInspectedValue(label, ct) )
    {
      Hostonly_ stamped_printf("\tInvalid keyword is described for '%s'\n", label.c_str());
      Exit(0);
    }
    GuideIn = ct;
    
    delete tp_index;
    
    
    
    // 統計値
    if ( C->Mode.Statistic == ON )
    {
      label="/StartCondition/Restart/DFIfiles/Statistic";
      
      if ( tpCntl->getInspectedValue(label, str ) )
      {
        f_dfi_in_stat = str.c_str();
      }
      if ( f_dfi_in_stat.empty() == true ) f_dfi_in_stat = "field_stat.dfi";
      
    }
  }
  
}


// #################################################################
// 固有オプションをロード
void PLT3D::getInherentOption()
{
  string label, str;

  // XYZファイル出力オプション
  label = "/Output/FormatOption/PLOT3D/XYZfile";
  
  if ( !tpCntl->getInspectedValue(label, str) )
  {
    Hostonly_ stamped_printf("\tInvalid keyword is described for '%s'\n", label.c_str());
    Exit(0);
  }
  if     ( !strcasecmp(str.c_str(), "on") )   XYZfile = ON;
  else if( !strcasecmp(str.c_str(), "off") )  XYZfile = OFF;
  else
  {
    Hostonly_ stamped_printf("\tInvalid keyword is described for '%s'\n", label.c_str());
    Exit(0);
  }
  
  
  // IBLANKファイルを使うオプション
  label = "/Output/FormatOption/PLOT3D/IblankFile";
  
  if ( !tpCntl->getInspectedValue(label, str) )
  {
    Hostonly_ stamped_printf("\tInvalid keyword is described for '%s'\n", label.c_str());
    Exit(0);
  }
  if     ( !strcasecmp(str.c_str(), "on") )   Iblank = ON;
  else if( !strcasecmp(str.c_str(), "off") )  Iblank = OFF;
  else
  {
    Hostonly_ stamped_printf("\tInvalid keyword is described for '%s'\n", label.c_str());
    Exit(0);
  }

}


// #################################################################
/* @brief ファイル出力の初期化
 */
void PLT3D::initFileOut(const int id_cell, const int id_bcf)
{
  // バッファサイズチェック  バッファサイズdnumはFALLOC::allocArray_Main()の数をハードコード
  // 書き出しサイズNvarsIns_plt3d, NvarsAvr_plt3dは、以下のメソッドに依存
  // Control::getSolverProperty()
  // IO_BASE::getFIOparams()
  
  
  // ソルバーでアロケートしたスカラー配列サイズ
  size_t dims[3];
  dims[0] = (size_t)(size[0] + 2*guide);
  dims[1] = (size_t)(size[1] + 2*guide);
  dims[2] = (size_t)(size[2] + 2*guide);
  size_allocated = dims[0] * dims[1] * dims[2];
  
  
  // 出力するスカラのバッファサイズ
  dims[0] = (size_t)(size[0] + 2*GuideOut);
  dims[1] = (size_t)(size[1] + 2*GuideOut);
  dims[2] = (size_t)(size[2] + 2*GuideOut);
  size_OutBuffer = dims[0] * dims[1] * dims[2];
  
  
  
  //printf("alloc=%zu out=%zu\n", size_allocated, size_OutBuffer);
  
  
  // IO バッファサイズ
  int dnum;
  if ( C->isHeatProblem() )
  {
    dnum = IO_BLOCK_SIZE_HEAT;
  }
  else
  {
    dnum = IO_BLOCK_SIZE_FLOW;
  }
  
  int NumVars     = C->NvarsIns_plt3d;
  int NumVarsStat = C->NvarsAvr_plt3d;
  
  if ( (NumVars > dnum) || (NumVarsStat > dnum) )
  {
    Hostonly_ printf("Error : Buffer size of PLOT3D exseeds predetermined value. Control the number of derived arrays for output.\n");
    Exit(0);
  }
  

  
  // Format
  CDM::E_CDM_FORMAT cdm_format = CDM::E_CDM_FMT_PLOT3D;

  
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
  
  
  std::string procfile="./proc.dfi";
  
  
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
  
  
  
  // 瞬時値と派生データ
  DFI_OUT_INS = cdm_DFI::WriteInit(MPI_COMM_WORLD, ///<MPI コミュニケータ
                                   cdm_DFI::Generate_DFI_Name(f_dfi_out_ins), ///<dfi ファイル名
                                   path,          ///<出力ディレクトリ
                                   f_dfi_out_ins, ///<ベースファイル名
                                   cdm_format,    ///<出力フォーマット
                                   GuideOut,      ///<出力仮想セル数
                                   datatype,      ///<データ型
                                   NumVars,       ///<データの変数の個数
                                   procfile,      ///<proc ファイル名
                                   G_size,        ///<計算空間全体のボクセルサイズ
                                   cdm_pit,       ///<ピッチ
                                   cdm_org,       ///<原点座標値
                                   cdm_div,       ///<領域分割数
                                   head,          ///<計算領域の開始位置
                                   cdm_tail,      ///<計算領域の終了位置
                                   HostName,      ///<ホスト名
                                   TimeSliceDir); ///<タイムスライス出力オプション
  
  if ( DFI_OUT_INS == NULL )
  {
    Hostonly_ stamped_printf("\tFails to instance Index.dfi.\n");
    Exit(0);
  }
  
  if( Slice == ON )
  {
    DFI_OUT_INS->SetTimeSliceFlag(CDM::E_CDM_ON);
  }
  else
  {
    DFI_OUT_INS->SetTimeSliceFlag(CDM::E_CDM_OFF);
  }
  
  
  DFI_OUT_INS->AddUnit("Length",   UnitL, (double)C->RefLength);
  DFI_OUT_INS->AddUnit("Velocity", UnitV, (double)C->RefVelocity);
  DFI_OUT_INS->AddUnit("Pressure", UnitP, (double)C->BasePrs, DiffPrs, true);
  
  
  // Proc file
  DFI_OUT_INS->WriteProcDfiFile(MPI_COMM_WORLD, true, id_cell, id_bcf);
  

  
  // XYZ file
  if ( XYZfile == ON )
  {
    
    // IBLANKファイルのモード
    if ( Iblank == OFF )
    {
      DFI_OUT_INS->WriteGridFile();
      Hostonly_ printf("\tXYZ file was successfully generated.\n");
    }
    else
    {
      if ( !d_mid ) Exit(0);
      generate_iblank_(d_mid, size, &guide, d_bcd, &GuideOut);
      
      DFI_OUT_INS->WriteGridFile(d_mid);
      Hostonly_ printf("\tXYZ file with IBLANK was successfully generated from BCindex.\n");
    }
  }
  
  
  
  // 統計値
  if ( C->Mode.Statistic == ON )
  {
    DFI_OUT_STAT = cdm_DFI::WriteInit(MPI_COMM_WORLD, ///< MPI コミュニケータ
                                      cdm_DFI::Generate_DFI_Name(f_dfi_out_stat), ///<dfi ファイル名
                                      path,           ///< 出力ディレクトリ
                                      f_dfi_out_stat, ///< ベースファイル名
                                      cdm_format,     ///< 出力フォーマット
                                      GuideOut,       ///< 出力仮想セル数
                                      datatype,       ///< データ型
                                      NumVarsStat,    ///< データの変数の個数
                                      procfile,       ///< proc ファイル名
                                      G_size,         ///< 計算空間全体のボクセルサイズ
                                      cdm_pit,        ///< ピッチ
                                      cdm_org,        ///< 原点座標値
                                      cdm_div,        ///< 領域分割数
                                      head,           ///< 計算領域の開始位置
                                      cdm_tail,       ///< 計算領域の終了位置
                                      HostName,       ///< ホスト名
                                      TimeSliceDir);  ///< タイムスライス出力オプション

    
    if ( DFI_OUT_STAT == NULL )
    {
      Hostonly_ stamped_printf("\tFails to instance AvrPressure dfi.\n");
      Exit(0);
    }
    
    if ( Slice == ON )
    {
      DFI_OUT_STAT->SetTimeSliceFlag(CDM::E_CDM_ON);
    }
    else
    {
      DFI_OUT_STAT->SetTimeSliceFlag(CDM::E_CDM_OFF);
    }
    
    DFI_OUT_STAT->AddUnit("Length"  , UnitL, (double)C->RefLength);
    DFI_OUT_STAT->AddUnit("Velocity", UnitV, (double)C->RefVelocity);
    DFI_OUT_STAT->AddUnit("Pressure", UnitP, (double)C->BasePrs, DiffPrs, true);
  }
  
}




// #################################################################
// 時間統計値のファイル出力
// @note d_ws,d_wvをワークに使用
void PLT3D::OutputStatisticalVarables(const unsigned m_CurrentStep,
                                      const double m_CurrentTime,
                                      const unsigned m_CurrentStepStat,
                                      const double m_CurrentTimeStat,
                                      double& flop)
{
  
  // packing  >> ap, av, atの順
  
  REAL_TYPE f_min, f_max, vec_min[3], vec_max[3];
  REAL_TYPE minmax[15*2];
  
  
  // エラーコード
  CDM::E_CDM_ERRORCODE ret;
  
  
  // 出力ファイルの指定が有次元の場合
  double timeStat;
  
  if (C->Unit.File == DIMENSIONAL)
  {
    timeStat = m_CurrentTimeStat * C->Tscale;
  }
  else
  {
    timeStat = m_CurrentTimeStat;
  }
  
  // 統計操作の母数
  unsigned stepStat = m_CurrentStepStat;
  
  
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
  
  REAL_TYPE unit_velocity = (C->Unit.File == DIMENSIONAL) ? C->RefVelocity : 1.0;
  
  int varN = 0; // 変数登録のインデクス
  int var  = 0; // minmaxのインデクス
  
  
  if ( !DFI_OUT_STAT )
  {
    printf("[%d] DFI_OUT_TEMPA Pointer Error\n", paraMngr->GetMyRankID());
    Exit(-1);
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
      blas_copy_(d_ws, d_ap, size, &guide);
    }
    
    // 最大値と最小値
    fb_minmax_s_ (&f_min, &f_max, size, &guide, d_ws, &flop);
    
    if ( numProc > 1 )
    {
      REAL_TYPE min_tmp = f_min;
      if ( paraMngr->Allreduce(&min_tmp, &f_min, 1, MPI_MIN) != CPM_SUCCESS ) Exit(0);
      
      REAL_TYPE max_tmp = f_max;
      if ( paraMngr->Allreduce(&max_tmp, &f_max, 1, MPI_MAX) != CPM_SUCCESS ) Exit(0);
    }
    
    minmax[var++] = f_min;  ///<<< p min
    minmax[var++] = f_max;  ///<<< p max
    
    
    pack_scalar_(&d_iobuf[size_OutBuffer*varN], size, &guide, d_ws, &GuideOut);
    DFI_OUT_STAT->setVariableName(varN++, l_avr_pressure);
    

    
    // Velocity
    fb_vout_ijkn_(d_wv, d_av, size, &guide, RF->getV00(), &unit_velocity, &flop);
    
    fb_minmax_v_ (vec_min, vec_max, size, &guide, RF->getV00(), d_wv, &flop);
    
    if ( numProc > 1 )
    {
      REAL_TYPE vmin_tmp[3] = {vec_min[0], vec_min[1], vec_min[2]};
      if ( paraMngr->Allreduce(vmin_tmp, vec_min, 3, MPI_MIN) != CPM_SUCCESS ) Exit(0);
      
      REAL_TYPE vmax_tmp[3] = {vec_max[0], vec_max[1], vec_max[2]};
      if ( paraMngr->Allreduce(vmax_tmp, vec_max, 3, MPI_MAX) != CPM_SUCCESS ) Exit(0);
    }
    
    minmax[var++] = vec_min[0]; ///<<< vec_u min
    minmax[var++] = vec_max[0]; ///<<< vec_u max
    minmax[var++] = vec_min[1]; ///<<< vec_v min
    minmax[var++] = vec_max[1]; ///<<< vec_v max
    minmax[var++] = vec_min[2]; ///<<< vec_w min
    minmax[var++] = vec_max[2]; ///<<< vec_w max
    

    pack_vector_(&d_iobuf[size_OutBuffer*varN], size, &guide, d_wv, &GuideOut);
    DFI_OUT_STAT->setVariableName(varN++, l_avr_velocity_x);
    DFI_OUT_STAT->setVariableName(varN++, l_avr_velocity_y);
    DFI_OUT_STAT->setVariableName(varN++, l_avr_velocity_z);
  }

  
  
  
  // Temperature
  if ( C->isHeatProblem() )
  {
    U.convArrayIE2Tmp(d_ws, size, guide, d_ae, d_bcd, mat_tbl, C->BaseTemp, C->DiffTemp, C->Unit.File, flop);
    
    fb_minmax_s_ (&f_min, &f_max, size, &guide, d_ws, &flop);
    
    if ( numProc > 1 )
    {
      REAL_TYPE min_tmp = f_min;
      if ( paraMngr->Allreduce(&min_tmp, &f_min, 1, MPI_MIN) != CPM_SUCCESS ) Exit(0);
      
      REAL_TYPE max_tmp = f_max;
      if ( paraMngr->Allreduce(&max_tmp, &f_max, 1, MPI_MAX) != CPM_SUCCESS ) Exit(0);
    }
    
    minmax[var++] = f_min; ///<<< t min
    minmax[var++] = f_max; ///<<< t max
    
    // ae
    pack_scalar_(&d_iobuf[size_OutBuffer*varN], size, &guide, d_ws, &GuideOut);
    DFI_OUT_STAT->setVariableName(varN++, l_avr_temperature);
  }
  
  
  // 統計
  if ( C->Mode.StatVelocity == ON )
  {
    // rms
    if (C->varState[var_RmsV] == ON )
    {
      fb_vout_ijkn_(d_wv, d_rms_v, size, &guide, RF->getV00(), &unit_velocity, &flop);
      
      fb_minmax_v_ (vec_min, vec_max, size, &guide, RF->getV00(), d_wv, &flop);
      
      if ( numProc > 1 )
      {
        REAL_TYPE vmin_tmp[3] = {vec_min[0], vec_min[1], vec_min[2]};
        if( paraMngr->Allreduce(vmin_tmp, vec_min, 3, MPI_MIN) != CPM_SUCCESS ) Exit(0);
        
        REAL_TYPE vmax_tmp[3] = {vec_max[0], vec_max[1], vec_max[2]};
        if( paraMngr->Allreduce(vmax_tmp, vec_max, 3, MPI_MAX) != CPM_SUCCESS ) Exit(0);
      }
      
      minmax[var++] = vec_min[0]; ///<<< vec_u min
      minmax[var++] = vec_max[0]; ///<<< vec_u max
      minmax[var++] = vec_min[1]; ///<<< vec_v min
      minmax[var++] = vec_max[1]; ///<<< vec_v max
      minmax[var++] = vec_min[2]; ///<<< vec_w min
      minmax[var++] = vec_max[2]; ///<<< vec_w max
      
      pack_vector_(&d_iobuf[size_OutBuffer*varN], size, &guide, d_wv, &GuideOut);
      DFI_OUT_STAT->setVariableName(varN++, l_rmsV_x);
      DFI_OUT_STAT->setVariableName(varN++, l_rmsV_y);
      DFI_OUT_STAT->setVariableName(varN++, l_rmsV_z);
    }
    
    // rms mean
    if (C->varState[var_RmsMeanV] == ON )
    {
      fb_vout_ijkn_(d_wv, d_rms_mean_v, size, &guide, RF->getV00(), &unit_velocity, &flop);
      
      fb_minmax_v_ (vec_min, vec_max, size, &guide, RF->getV00(), d_wv, &flop);
      
      if ( numProc > 1 )
      {
        REAL_TYPE vmin_tmp[3] = {vec_min[0], vec_min[1], vec_min[2]};
        if( paraMngr->Allreduce(vmin_tmp, vec_min, 3, MPI_MIN) != CPM_SUCCESS ) Exit(0);
        
        REAL_TYPE vmax_tmp[3] = {vec_max[0], vec_max[1], vec_max[2]};
        if( paraMngr->Allreduce(vmax_tmp, vec_max, 3, MPI_MAX) != CPM_SUCCESS ) Exit(0);
      }
      
      minmax[var++] = vec_min[0]; ///<<< vec_u min
      minmax[var++] = vec_max[0]; ///<<< vec_u max
      minmax[var++] = vec_min[1]; ///<<< vec_v min
      minmax[var++] = vec_max[1]; ///<<< vec_v max
      minmax[var++] = vec_min[2]; ///<<< vec_w min
      minmax[var++] = vec_max[2]; ///<<< vec_w max
      
      pack_vector_(&d_iobuf[size_OutBuffer*varN], size, &guide, d_wv, &GuideOut);
      DFI_OUT_STAT->setVariableName(varN++, l_rmsmeanV_x);
      DFI_OUT_STAT->setVariableName(varN++, l_rmsmeanV_y);
      DFI_OUT_STAT->setVariableName(varN++, l_rmsmeanV_z);
    }
  }
  
  
  if ( C->Mode.StatPressure == ON )
  {
    // rms
    if (C->varState[var_RmsP] == ON )
    {
      if (C->Unit.File == DIMENSIONAL)
      {
        REAL_TYPE bp = ( C->Unit.Prs == Unit_Absolute ) ? C->BasePrs : 0.0;
        U.convArrayPrsND2D(d_ws, size, guide, d_rms_p, bp, C->RefDensity, C->RefVelocity, flop);
      }
      else
      {
        blas_copy_(d_ws, d_rms_p, size, &guide);
      }
      
      // 最大値と最小値
      fb_minmax_s_ (&f_min, &f_max, size, &guide, d_ws, &flop);
      
      if ( numProc > 1 )
      {
        REAL_TYPE min_tmp = f_min;
        if ( paraMngr->Allreduce(&min_tmp, &f_min, 1, MPI_MIN) != CPM_SUCCESS ) Exit(0);
        
        REAL_TYPE max_tmp = f_max;
        if ( paraMngr->Allreduce(&max_tmp, &f_max, 1, MPI_MAX) != CPM_SUCCESS ) Exit(0);
      }
      
      minmax[var++] = f_min;  ///<<< p min
      minmax[var++] = f_max;  ///<<< p max
      
      pack_scalar_(&d_iobuf[size_OutBuffer*varN], size, &guide, d_ws, &GuideOut);
      DFI_OUT_STAT->setVariableName(varN++, l_rmsP);
    }
    
    // rms mean
    if (C->varState[var_RmsMeanP] == ON )
    {
      if (C->Unit.File == DIMENSIONAL)
      {
        REAL_TYPE bp = ( C->Unit.Prs == Unit_Absolute ) ? C->BasePrs : 0.0;
        U.convArrayPrsND2D(d_ws, size, guide, d_rms_mean_p, bp, C->RefDensity, C->RefVelocity, flop);
      }
      else
      {
        blas_copy_(d_ws, d_rms_mean_p, size, &guide);
      }
      
      // 最大値と最小値
      fb_minmax_s_ (&f_min, &f_max, size, &guide, d_ws, &flop);
      
      if ( numProc > 1 )
      {
        REAL_TYPE min_tmp = f_min;
        if ( paraMngr->Allreduce(&min_tmp, &f_min, 1, MPI_MIN) != CPM_SUCCESS ) Exit(0);
        
        REAL_TYPE max_tmp = f_max;
        if ( paraMngr->Allreduce(&max_tmp, &f_max, 1, MPI_MAX) != CPM_SUCCESS ) Exit(0);
      }
      
      minmax[var++] = f_min;  ///<<< p min
      minmax[var++] = f_max;  ///<<< p max
      
      pack_scalar_(&d_iobuf[size_OutBuffer*varN], size, &guide, d_ws, &GuideOut);
      DFI_OUT_STAT->setVariableName(varN++, l_rmsmeanP);
    }
  }

  
  if ( C->Mode.StatTemperature == ON )
  {
    // rms
    if (C->varState[var_RmsT] == ON )
    {
      U.convArrayIE2Tmp(d_ws, size, guide, d_rms_t, d_bcd, mat_tbl, C->BaseTemp, C->DiffTemp, C->Unit.File, flop);
      
      fb_minmax_s_ (&f_min, &f_max, size, &guide, d_ws, &flop);
      
      if ( numProc > 1 )
      {
        REAL_TYPE min_tmp = f_min;
        if ( paraMngr->Allreduce(&min_tmp, &f_min, 1, MPI_MIN) != CPM_SUCCESS ) Exit(0);
        
        REAL_TYPE max_tmp = f_max;
        if ( paraMngr->Allreduce(&max_tmp, &f_max, 1, MPI_MAX) != CPM_SUCCESS ) Exit(0);
      }
      
      minmax[var++] = f_min; ///<<< t min
      minmax[var++] = f_max; ///<<< t max
      
      pack_scalar_(&d_iobuf[size_OutBuffer*varN], size, &guide, d_ws, &GuideOut);
      DFI_OUT_STAT->setVariableName(varN++, l_rmsT);
    }
    
    // rms mean
    if (C->varState[var_RmsMeanT] == ON )
    {
      U.convArrayIE2Tmp(d_ws, size, guide, d_rms_mean_t, d_bcd, mat_tbl, C->BaseTemp, C->DiffTemp, C->Unit.File, flop);
      
      fb_minmax_s_ (&f_min, &f_max, size, &guide, d_ws, &flop);
      
      if ( numProc > 1 )
      {
        REAL_TYPE min_tmp = f_min;
        if ( paraMngr->Allreduce(&min_tmp, &f_min, 1, MPI_MIN) != CPM_SUCCESS ) Exit(0);
        
        REAL_TYPE max_tmp = f_max;
        if ( paraMngr->Allreduce(&max_tmp, &f_max, 1, MPI_MAX) != CPM_SUCCESS ) Exit(0);
      }
      
      minmax[var++] = f_min; ///<<< t min
      minmax[var++] = f_max; ///<<< t max
      
      pack_scalar_(&d_iobuf[size_OutBuffer*varN], size, &guide, d_ws, &GuideOut);
      DFI_OUT_STAT->setVariableName(varN++, l_rmsmeanT);
    }
  }
  
  
  if ( C->NvarsAvr_plt3d != varN )
  {
    Hostonly_ printf("\tThe number of variables to be written is inconsistent. NvarsAvr=%d, varN=%d\n",
                     C->NvarsAvr_plt3d, varN);
    Exit(0);
  }

  
  ret = DFI_OUT_STAT->WriteData(m_step,     // 出力step番号
                                m_time,     // 出力時刻
                                size,       // 配列サイズ
                                varN,       // 成分数
                                GuideOut,   // 出力仮想セル数
                                d_iobuf,    // フィールドデータポインタ
                                minmax,     // 最小値と最大値
                                false,      // 統計出力指示 false:出力あり
                                stepStat,   // 統計をとったステップ数
                                timeStat);  // 統計をとった時刻
  
  if ( ret != CDM::E_CDM_SUCCESS ) Exit(0);
}


// #################################################################
// 基本変数のファイル出力
// @note d_ws,d_wvをワークに使用
void PLT3D::OutputBasicVariables(const unsigned m_CurrentStep,
                                 const double m_CurrentTime,
                                 double& flop)
{
  // packing >> p, v, vf, t, tp, vor, q, hlt, divの順（全15 scalarあるが，バッファ以上はinitOutFile()でチェック）
  
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
  
  
  // 最大値と最小値
  REAL_TYPE f_min, f_max, vec_min[3], vec_max[3];
  REAL_TYPE minmax[15*2];
  
  
  // エラーコード
  CDM::E_CDM_ERRORCODE ret;
  
  REAL_TYPE unit_velocity = (C->Unit.File == DIMENSIONAL) ? C->RefVelocity : 1.0;
  
  
  if( !DFI_OUT_INS )
  {
    printf("[%d] DFI_OUT_INS Pointer Error\n", paraMngr->GetMyRankID());
    Exit(-1);
  }
  
  
  int varN = 0; // 変数登録のインデクス
  int var  = 0; // minmaxのインデクス
  
  
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
      blas_copy_(d_ws, d_p, size, &guide);
    }
    
    fb_minmax_s_ (&f_min, &f_max, size, &guide, d_ws, &flop);
    
    if ( numProc > 1 )
    {
      REAL_TYPE min_tmp = f_min;
      if( paraMngr->Allreduce(&min_tmp, &f_min, 1, MPI_MIN) != CPM_SUCCESS ) Exit(0);
      
      REAL_TYPE max_tmp = f_max;
      if( paraMngr->Allreduce(&max_tmp, &f_max, 1, MPI_MAX) != CPM_SUCCESS ) Exit(0);
    }
    minmax[var++] = f_min; ///<<< p min
    minmax[var++] = f_max; ///<<< p max

    
    // p
    pack_scalar_(&d_iobuf[size_OutBuffer*varN], size, &guide, d_ws, &GuideOut);
    DFI_OUT_INS->setVariableName(varN++, l_pressure);
    
    
    
    // Velocity
    fb_vout_ijkn_(d_wv, d_v, size, &guide, RF->getV00(), &unit_velocity, &flop);
    
    fb_minmax_v_ (vec_min, vec_max, size, &guide, RF->getV00(), d_wv, &flop);
    
    if ( numProc > 1 )
    {
      REAL_TYPE vmin_tmp[3] = {vec_min[0], vec_min[1], vec_min[2]};
      if( paraMngr->Allreduce(vmin_tmp, vec_min, 3, MPI_MIN) != CPM_SUCCESS ) Exit(0);
      
      REAL_TYPE vmax_tmp[3] = {vec_max[0], vec_max[1], vec_max[2]};
      if( paraMngr->Allreduce(vmax_tmp, vec_max, 3, MPI_MAX) != CPM_SUCCESS ) Exit(0);
    }
    
    minmax[var++] = vec_min[0]; ///<<< vec_u min
    minmax[var++] = vec_max[0]; ///<<< vec_u max
    minmax[var++] = vec_min[1]; ///<<< vec_v min
    minmax[var++] = vec_max[1]; ///<<< vec_v max
    minmax[var++] = vec_min[2]; ///<<< vec_w min
    minmax[var++] = vec_max[2]; ///<<< vec_w max


    pack_vector_(&d_iobuf[size_OutBuffer*varN], size, &guide, d_wv, &GuideOut);
    DFI_OUT_INS->setVariableName(varN++, l_velocity_x);
    DFI_OUT_INS->setVariableName(varN++, l_velocity_y);
    DFI_OUT_INS->setVariableName(varN++, l_velocity_z);
    
    
    
    // Face Velocity
    fb_vout_ijkn_(d_wv, d_vf, size, &guide, RF->getV00(), &unit_velocity, &flop);
    
    fb_minmax_v_ (vec_min, vec_max, size, &guide, RF->getV00(), d_wv, &flop);
    
    if ( numProc > 1 )
    {
      REAL_TYPE vmin_tmp[3] = {vec_min[0], vec_min[1], vec_min[2]};
      if( paraMngr->Allreduce(vmin_tmp, vec_min, 3, MPI_MIN) != CPM_SUCCESS ) Exit(0);
      
      REAL_TYPE vmax_tmp[3] = {vec_max[0], vec_max[1], vec_max[2]};
      if( paraMngr->Allreduce(vmax_tmp, vec_max, 3, MPI_MAX) != CPM_SUCCESS ) Exit(0);
    }
    
    minmax[var++] = vec_min[0]; ///<<< vec_u min
    minmax[var++] = vec_max[0]; ///<<< vec_u max
    minmax[var++] = vec_min[1]; ///<<< vec_v min
    minmax[var++] = vec_max[1]; ///<<< vec_v max
    minmax[var++] = vec_min[2]; ///<<< vec_w min
    minmax[var++] = vec_max[2]; ///<<< vec_w max

    
    pack_vector_(&d_iobuf[size_OutBuffer*varN], size, &guide, d_wv, &GuideOut);
    DFI_OUT_INS->setVariableName(varN++, l_fvelocity_x);
    DFI_OUT_INS->setVariableName(varN++, l_fvelocity_y);
    DFI_OUT_INS->setVariableName(varN++, l_fvelocity_z);
  }
  
  
  if ( C->isHeatProblem() )
  {
    // Tempearture
    U.convArrayIE2Tmp(d_ws, size, guide, d_ie, d_bcd, mat_tbl, C->BaseTemp, C->DiffTemp, C->Unit.File, flop);
    
    fb_minmax_s_ (&f_min, &f_max, size, &guide, d_ws, &flop);
    
    if ( numProc > 1 )
    {
      REAL_TYPE min_tmp = f_min;
      if( paraMngr->Allreduce(&min_tmp, &f_min, 1, MPI_MIN) != CPM_SUCCESS ) Exit(0);
      
      REAL_TYPE max_tmp = f_max;
      if( paraMngr->Allreduce(&max_tmp, &f_max, 1, MPI_MAX) != CPM_SUCCESS ) Exit(0);
    }
    minmax[var++] = f_min;
    minmax[var++] = f_max;
    
    pack_scalar_(&d_iobuf[size_OutBuffer*varN], size, &guide, d_ws, &GuideOut);
    DFI_OUT_INS->setVariableName(varN++, l_temperature);
  }
  
  
  
  // 派生変数
  
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
      REAL_TYPE min_tmp = f_min;
      if( paraMngr->Allreduce(&min_tmp, &f_min, 1, MPI_MIN) != CPM_SUCCESS ) Exit(0);
      
      REAL_TYPE max_tmp = f_max;
      if( paraMngr->Allreduce(&max_tmp, &f_max, 1, MPI_MAX) != CPM_SUCCESS ) Exit(0);
    }
    minmax[var++] = f_min;
    minmax[var++] = f_max;
    
    pack_scalar_(&d_iobuf[size_OutBuffer*varN], size, &guide, d_ws, &GuideOut);
    DFI_OUT_INS->setVariableName(varN++, l_totalp);
  }
  
  
  // Vorticity
  if (C->varState[var_Vorticity] == ON )
  {
    rot_v_ (d_wv, size, &guide, &deltaX, d_v, d_cdf, RF->getV00(), &flop);
    
    REAL_TYPE vz[3]; // dummy
    vz[0] = vz[1] = vz[2] = 0.0;
    REAL_TYPE unit_vort = (C->Unit.File == DIMENSIONAL) ? C->RefVelocity/C->RefLength : 1.0;
    
    fb_vout_ijkn_(d_wv, d_wv, size, &guide, vz, &unit_vort, &flop);
    fb_minmax_v_ (vec_min, vec_max, size, &guide, RF->getV00(), d_wv, &flop);
    
    if ( numProc > 1 )
    {
      REAL_TYPE vmin_tmp[3] = {vec_min[0], vec_min[1], vec_min[2]};
      if( paraMngr->Allreduce(vmin_tmp, vec_min, 3, MPI_MIN) != CPM_SUCCESS ) Exit(0);
      
      REAL_TYPE vmax_tmp[3] = {vec_max[0], vec_max[1], vec_max[2]};
      if( paraMngr->Allreduce(vmax_tmp, vec_max, 3, MPI_MAX) != CPM_SUCCESS ) Exit(0);
    }
    
    minmax[var++] = vec_min[0]; ///<<< vec_u min
    minmax[var++] = vec_max[0]; ///<<< vec_u max
    minmax[var++] = vec_min[1]; ///<<< vec_v min
    minmax[var++] = vec_max[1]; ///<<< vec_v max
    minmax[var++] = vec_min[2]; ///<<< vec_w min
    minmax[var++] = vec_max[2]; ///<<< vec_w max
    
    pack_vector_(&d_iobuf[size_OutBuffer*varN], size, &guide, d_wv, &GuideOut);
    DFI_OUT_INS->setVariableName(varN++, l_vorticity_x);
    DFI_OUT_INS->setVariableName(varN++, l_vorticity_y);
    DFI_OUT_INS->setVariableName(varN++, l_vorticity_z);
  }

  
  // 2nd Invariant of Velocity Gradient Tensor 無次元で出力
  if (C->varState[var_Qcr] == ON )
  {
    i2vgt_ (d_ws, size, &guide, &deltaX, d_v, d_cdf, RF->getV00(), &flop);
    
    fb_minmax_s_ (&f_min, &f_max, size, &guide, d_ws, &flop);
    
    if ( numProc > 1 )
    {
      REAL_TYPE min_tmp = f_min;
      if( paraMngr->Allreduce(&min_tmp, &f_min, 1, MPI_MIN) != CPM_SUCCESS ) Exit(0);
      
      REAL_TYPE max_tmp = f_max;
      if( paraMngr->Allreduce(&max_tmp, &f_max, 1, MPI_MAX) != CPM_SUCCESS ) Exit(0);
    }
    minmax[var++] = f_min;
    minmax[var++] = f_max;
    
    pack_scalar_(&d_iobuf[size_OutBuffer*varN], size, &guide, d_ws, &GuideOut);
    DFI_OUT_INS->setVariableName(varN++, l_invariantQ);
  }
  
  
  // Helicity 無次元で出力
  if (C->varState[var_Helicity] == ON )
  {
    helicity_ (d_ws, size, &guide, &deltaX, d_v, d_cdf, RF->getV00(), &flop);
    
    fb_minmax_s_ (&f_min, &f_max, size, &guide, d_ws, &flop);
    
    if ( numProc > 1 )
    {
      REAL_TYPE min_tmp = f_min;
      if( paraMngr->Allreduce(&min_tmp, &f_min, 1, MPI_MIN) != CPM_SUCCESS ) Exit(0);
      
      REAL_TYPE max_tmp = f_max;
      if( paraMngr->Allreduce(&max_tmp, &f_max, 1, MPI_MAX) != CPM_SUCCESS ) Exit(0);
    }

    minmax[var++] = f_min;
    minmax[var++] = f_max;
    
    pack_scalar_(&d_iobuf[size_OutBuffer*varN], size, &guide, d_ws, &GuideOut);
    DFI_OUT_INS->setVariableName(varN++, l_helicity);
  }
  
  
  // Divergence for Debug
  if (C->varState[var_Div] == ON )
  {
    REAL_TYPE coef = 1.0/deltaX; /// 発散値を計算するための係数　1/h
    U.cnv_Div(d_ws, d_dv, size, guide, coef);
    
    fb_minmax_s_ (&f_min, &f_max, size, &guide, d_ws, &flop);
    
    if ( numProc > 1 )
    {
      REAL_TYPE min_tmp = f_min;
      if( paraMngr->Allreduce(&min_tmp, &f_min, 1, MPI_MIN) != CPM_SUCCESS ) Exit(0);
      
      REAL_TYPE max_tmp = f_max;
      if( paraMngr->Allreduce(&max_tmp, &f_max, 1, MPI_MAX) != CPM_SUCCESS ) Exit(0);
    }

    minmax[var++] = f_min;
    minmax[var++] = f_max;
    
    pack_scalar_(&d_iobuf[size_OutBuffer*varN], size, &guide, d_ws, &GuideOut); // << check
    DFI_OUT_INS->setVariableName(varN++, l_divergence);
  }
  
  if ( C->NvarsIns_plt3d != varN )
  {
    Hostonly_ printf("\tThe number of variables to be written is inconsistent. NvarsIns=%d, varN=%d\n",
                     C->NvarsIns_plt3d, varN);
    Exit(0);
  }
  

  ret = DFI_OUT_INS->WriteData(m_step,
                               m_time,
                               size,
                               varN,
                               GuideOut, // 出力するガイドセル数 << check
                               d_iobuf,
                               minmax,
                               true,
                               0,
                               0.0);
  
  if ( ret != CDM::E_CDM_SUCCESS ) Exit(0);

}


// #################################################################
// 固有の制御パラメータSTEERの表示
void PLT3D::printSteerConditionsInherent(FILE* fp)
{
  fprintf(fp,"\t     XYZ file                 :   %s\n", (XYZfile==ON) ? "On" : "Off");
  fprintf(fp,"\t     Iblank option            :   %s\n", (Iblank==ON) ? "On" : "Off");
}


// #################################################################
// リスタート時の統計値ファイル読み込み
void PLT3D::RestartStatistic(FILE* fp,
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
  
  
  
  // まだ統計値開始時刻になっていなければ，何もしない
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
    // 読込み用インスタンスのポインタ取得
    DFI_IN_STAT = cdm_DFI::ReadInit(MPI_COMM_WORLD,
                                   f_dfi_in_stat,
                                   G_size,
                                   gdiv,
                                   cdm_error);
    
    if ( cdm_error != CDM::E_CDM_SUCCESS || DFI_IN_STAT == NULL  ) Exit(0);
  }
  
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
  
  // 読込みフィールドデータ型のチェック
  if( DFI_IN_STAT->GetDataType() != datatype )
  {
    Hostonly_ printf("Error Datatype unmatch.\n");
    Exit(0);
  }
  

  
  int dnum;
  if ( C->isHeatProblem() )
  {
    dnum = IO_BLOCK_SIZE_HEAT;
  }
  else
  {
    dnum = IO_BLOCK_SIZE_FLOW;
  }
  
  
  // 読込みフィールドデータの成分数を取得
  int NumVarsStat = DFI_IN_STAT->GetNumVariables();
  
  
  if ( NumVarsStat > dnum )
  {
    Hostonly_ printf("Error : The number of stored variables exceeds current limitation.\n");
    Exit(0);
  }
  
  
  // 領域クリア
  memset(d_iobuf, 0, sizeof(REAL_TYPE)*size_allocated*dnum);
  
  
  
  unsigned step_stat = 0;
  double time_stat = 0.0;
  
  
  //自身の領域終点インデックス
  int tail[3];
  for (int i=0; i<3; i++) tail[i] = head[i]+size[i]-1;
  
  
  double r_time;
  if ( DFI_IN_STAT->ReadData(d_iobuf,        ///< 読込み先配列のポインタ
                             m_RestartStep,  ///< 読込みフィールドデータのステップ番号
                             GuideIn,        ///< 計算空間の仮想セル数
                             G_size,         ///< 計算空間全体のボクセルサイズ
                             gdiv,           ///< 領域分割数
                             head,           ///< 自サブドメインの開始インデクス
                             tail,           ///< 自サブドメインの終端インデクス
                             r_time,         ///< dfi から読込んだ時間
                             false,          ///< 統計値指定
                             step_stat,
                             time_stat) != CDM::E_CDM_SUCCESS ) Exit(0);
  
  if( d_iobuf == NULL ) Exit(0);
  
  
  m_CurrentStepStat = step_stat;
  m_CurrentTimeStat = time_stat;
  
  
  REAL_TYPE refv = (C->Unit.File == DIMENSIONAL) ? C->RefVelocity : 1.0;
  REAL_TYPE u0[4];
  
  RF->copyV00(u0);
  
  REAL_TYPE bp = ( C->Unit.Prs == Unit_Absolute ) ? C->BasePrs : 0.0;
  
  // 変数名をキーにしてデータロード
  int block = 0;
  while ( block < NumVarsStat )
  {
    string variable = DFI_IN_STAT->getVariableName(block);
    
    if ( variable.empty() == true || variable == "" ) Exit(0);
    
    // pressure
    if ( !strcasecmp(variable.c_str(), l_avr_pressure.c_str()) )
    {
      unpack_scalar_(d_ap, size, &guide, &d_iobuf[size_InBuffer*block], &GuideIn);
      block++;
      
      if ( C->Unit.File == DIMENSIONAL ) // 有次元の場合，無次元に変換する
      {
        U.convArrayPrsD2ND(d_ap, size, guide, bp, C->RefDensity, C->RefVelocity, flop);
      }
    }
    
    // velocity
    else if ( !strcasecmp(variable.c_str(), l_avr_velocity_x.c_str()) ) // 先頭はvec_Xが書かれていると仮定
    {
      unpack_vector_(d_av, size, &guide, &d_iobuf[size_InBuffer*block], &GuideIn);
      block += 3;
      
      fb_vin_ijkn_(d_av, size, &guide, u0, &refv, &flop);
    }
    
    // temperature
    else if ( !strcasecmp(variable.c_str(), l_avr_temperature.c_str()) )
    {
      if ( C->isHeatProblem() )
      {
        unpack_scalar_(d_ws, size, &guide, &d_iobuf[size_InBuffer*block], &GuideIn);
        block++;
        
        U.convArrayTmp2IE(d_ae, size, guide, d_ws, d_bcd, mat_tbl, C->BaseTemp, C->DiffTemp, C->Unit.File, flop);
      }
    }
    
    // rms mean V
    else if ( !strcasecmp(variable.c_str(), l_rmsmeanV_x.c_str()) && C->Mode.StatVelocity==ON )
    {
      unpack_vector_(d_rms_mean_v, size, &guide, &d_iobuf[size_InBuffer*block], &GuideIn);
      block += 3;

      fb_vin_ijkn_(d_rms_mean_v, size, &guide, u0, &refv, &flop);
    }
    
    // rms mean P
    else if ( !strcasecmp(variable.c_str(), l_rmsmeanP.c_str()) && C->Mode.StatPressure==ON )
    {
      unpack_scalar_(d_rms_mean_p, size, &guide, &d_iobuf[size_InBuffer*block], &GuideIn);
      block ++;
      
      if ( C->Unit.File == DIMENSIONAL ) // 有次元の場合，無次元に変換する
      {
        U.convArrayPrsD2ND(d_rms_mean_p, size, guide, bp, C->RefDensity, C->RefVelocity, flop);
      }
    }
    
    // rms mean T
    else if ( !strcasecmp(variable.c_str(), l_rmsmeanT.c_str()) && C->Mode.StatTemperature==ON )
    {
      unpack_scalar_(d_ws, size, &guide, &d_iobuf[size_InBuffer*block], &GuideIn);
      block ++;
      
      U.convArrayTmp2IE(d_rms_mean_t, size, guide, d_ws, d_bcd, mat_tbl, C->BaseTemp, C->DiffTemp, C->Unit.File, flop);
    }
    
  }
  
  Hostonly_ printf     ("\tStatistical fields have read :\tRestart step=%d time=%e  : Statistical step=%d time=%e [%s]\n",
                        m_RestartStep,
                        r_time,
                        step_stat,
                        time_stat, (C->Unit.File == DIMENSIONAL)?"sec.":"-");
  Hostonly_ fprintf(fp, "\tStatistical fields have read :\tRestart step=%d time=%e  : Statistical step=%d time=%e [%s]\n",
                    m_RestartStep,
                    r_time,
                    step_stat,
                    time_stat, (C->Unit.File == DIMENSIONAL)?"sec.":"-");
  
}



// #################################################################
// リスタート時の瞬時値ファイル読み込み
void PLT3D::RestartInstantaneous(FILE* fp,
                                 unsigned& m_CurrentStep,
                                 double& m_CurrentTime,
                                 double& flop)
{
  std::string fmt(file_fmt_ext);
  
  REAL_TYPE bp = ( C->Unit.Prs == Unit_Absolute ) ? C->BasePrs : 0.0;
  
  
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
  
  
  
  // 現在のセッションの領域分割数の取得
  int gdiv[3] = {1, 1, 1};
  if ( numProc > 1)
  {
    const int* m_div = paraMngr->GetDivNum();
    for (int i=0; i<3; i++ ) gdiv[i]=m_div[i];
  }
  
  
  // 自身の領域終点インデクス
  int tail[3];
  for (int i=0;i<3;i++) tail[i]=head[i]+size[i]-1;
  
  
  
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
  
  // 読込みフィールドデータ型のチェック
  if( DFI_IN_INS->GetDataType() != datatype )
  {
    Hostonly_ printf("Error Datatype unmatch.\n");
    Exit(0);
  }
  
  
  int dnum;
  if ( C->isHeatProblem() )
  {
    dnum = IO_BLOCK_SIZE_HEAT;
  }
  else
  {
    dnum = IO_BLOCK_SIZE_FLOW;
  }
  
  
  // 読込みフィールドデータの成分数を取得
  int NumVars = DFI_IN_INS->GetNumVariables();
  
  
  if ( NumVars > dnum )
  {
    Hostonly_ printf("Error : The number of stored variables exceeds current limitation.\n");
    Exit(0);
  }
  if ( NumVars != C->NvarsIns_plt3d )
  {
    Hostonly_ printf("Error : The number of stored variables is not consistent with specified value %d.\n", C->NvarsIns_plt3d);
    Exit(0);
  }
  
    
    
  // バッファクリア
  memset(d_iobuf, 0, sizeof(REAL_TYPE)*size_allocated*dnum);

  
  
  // dummy
  unsigned i_dummy=0;
  double f_dummy=0.0;
  
  double r_time;

  
  if ( DFI_IN_INS->ReadData(d_iobuf,        ///< 読込み先配列のポインタ
                            m_RestartStep,  ///< 読込みフィールドデータのステップ番号
                            GuideIn,        ///< 計算空間の仮想セル数
                            G_size,         ///< 計算空間全体のボクセルサイズ
                            gdiv,           ///< 領域分割数
                            head,           ///< 自サブドメインの開始インデクス
                            tail,           ///< 自サブドメインの終端インデクス
                            r_time,         ///< dfi から読込んだ時間
                            true,           ///< 瞬時値値指定
                            i_dummy,
                            f_dummy) != CDM::E_CDM_SUCCESS ) Exit(0);
  
  if( d_iobuf == NULL ) Exit(0);

  
  // ここでタイムスタンプを得る
  if (C->Unit.File == DIMENSIONAL) r_time /= C->Tscale;
  m_CurrentStep = m_RestartStep;
  m_CurrentTime = r_time;
  
  // v00[]に値をセット
  RF->setV00(r_time);
  
  
  REAL_TYPE refv = (C->Unit.File == DIMENSIONAL) ? C->RefVelocity  : 1.0;
  REAL_TYPE u0[4];
  RF->copyV00(u0);
  
  
  
  // 変数名をキーにしてデータロード @todo 書きこまれた変数を確実に読むアルゴに変更
  int block = 0;
  while ( block < NumVars )
  {
    string variable = DFI_IN_INS->getVariableName(block);
    
    if ( variable.empty() == true || variable == "" ) Exit(0);
    

    // pressure
    if ( !strcasecmp(variable.c_str(), l_pressure.c_str()) )
    {
      unpack_scalar_(d_p, size, &guide, &d_iobuf[size_InBuffer*block], &GuideIn);
      block++;
      
      if ( C->Unit.File == DIMENSIONAL ) // 有次元の場合，無次元に変換する
      {
        U.convArrayPrsD2ND(d_p, size, guide, bp, C->RefDensity, C->RefVelocity, flop);
      }
    }
    // Velocity
    else if ( !strcasecmp(variable.c_str(), l_velocity_x.c_str()) ) // 先頭はvec_Xが書かれていると仮定
    {
      unpack_vector_(d_v, size, &guide, &d_iobuf[size_InBuffer*block], &GuideIn);
      block += 3;
      
      fb_vin_ijkn_(d_v, size, &guide, u0, &refv, &flop);
    }
    
    // Fvelocity
    else if ( !strcasecmp(variable.c_str(), l_fvelocity_x.c_str()) ) // 先頭はvec_Xが書かれていると仮定
    {
      unpack_vector_(d_vf, size, &guide, &d_iobuf[size_InBuffer*block], &GuideIn);
      block += 3;

      fb_vin_ijkn_(d_vf, size, &guide, u0, &refv, &flop);
    }
    
    // Temperature
    else if ( !strcasecmp(variable.c_str(), l_temperature.c_str()) )
    {
      if ( C->isHeatProblem() )
      {
        unpack_scalar_(d_ws, size, &guide, &d_iobuf[size_InBuffer*block], &GuideIn);
        block++;
        U.convArrayTmp2IE(d_ie, size, guide, d_ws, d_bcd, mat_tbl, C->BaseTemp, C->DiffTemp, C->Unit.File, flop);
      }
    }
    
  }
  
  Hostonly_ printf     ("\n\tField data have read :\tstep=%d  time=%e [%s]\n",
                        m_RestartStep, r_time, (C->Unit.File == DIMENSIONAL)?"sec.":"-");
  Hostonly_ fprintf(fp, "\n\tField data have read :\tstep=%d  time=%e [%s]\n",
                    m_RestartStep, r_time, (C->Unit.File == DIMENSIONAL)?"sec.":"-");
  
}



// #################################################################
// リスタートモードを判定し，瞬時値のデータをロードする
void PLT3D::Restart(FILE* fp, unsigned& m_CurrentStep, double& m_CurrentTime)
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
    
    DFI_IN_INS = cdm_DFI::ReadInit(MPI_COMM_WORLD,
                                   f_dfi_in_ins,
                                   G_size,
                                   gdiv,
                                   cdm_error);
    
    if ( cdm_error != CDM::E_CDM_SUCCESS || DFI_IN_INS == NULL ) Exit(0);
    
    
    
    bool isSameDiv = true; // 同一分割数
    bool isSameRes = true; // 同一解像度
    
    
    // 前のセッションの領域分割数の取得
    const int* DFI_div=NULL;
    
    DFI_div = DFI_IN_INS->GetDFIGlobalDivision();
    
    
    
    // 前セッションと領域分割数が異なる場合
    for (int i=0; i<3; i++ )
    {
      if ( gdiv[i] != DFI_div[i] )
      {
        isSameDiv = false;
      }
    }
    
    // 前のセッションの全要素数の取得
    const int* DFI_G_size = DFI_IN_INS->GetDFIGlobalVoxel();
    
    
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
  
  size_t dims[3];
  dims[0] = (size_t)(size[0] + 2*GuideIn);
  dims[1] = (size_t)(size[1] + 2*GuideIn);
  dims[2] = (size_t)(size[2] + 2*GuideIn);
  size_InBuffer = dims[0] * dims[1] * dims[2];
  
  
  double flop_task = 0.0;
  RestartInstantaneous(fp, m_CurrentStep, m_CurrentTime, flop_task);

}

