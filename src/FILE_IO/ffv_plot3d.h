#ifndef _FFV_PLT3D_H_
#define _FFV_PLT3D_H_

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
// Copyright (c) 2012-2016 Advanced Institute for Computational Science, RIKEN.
// All rights reserved.
//
// Copyright (c) 2016-2020 Research Institute for Information Research, Kyushu University.
// All rights reserved.
//
//##################################################################################
//

/*
 * @file   ffv_plot3d.h
 * @brief  File IO of PLOT3D Class Header
 * @author aics
 */


#include "ffv_io_base.h"


class PLT3D : public IO_BASE {
  
private:
  
  
  size_t size_allocated; ///< ソルバーでアロケートされたスカラ配列のサイズ >> InitFileOut()で計算
  size_t size_OutBuffer; ///< 出力時のバッファサイズ（スカラ配列）
  size_t size_InBuffer;  ///< リスタート入力時のバッファサイズ
  int XYZfile;           ///< XYZ fike出力オプション
  int Iblank;            ///< PLOT3DのIBLANKオプション
  
  
  // 入力dfiファイルのプレフィックス
  string f_dfi_in_ins;
  string f_dfi_in_stat;

  // 出力ファイルのプレフィックス
  string f_dfi_out_ins;
  string f_dfi_out_stat;
  
  // InFile
  cdm_DFI *DFI_IN_INS;       ///< Instantaneus field
  cdm_DFI *DFI_IN_STAT;      ///< Statistical field
  
  // OutFile
  cdm_DFI *DFI_OUT_INS;      ///< Instantaneus field
  cdm_DFI *DFI_OUT_STAT;     ///< Statistical field
  
  // ラベル名
  string l_divergence;
  string l_helicity;
  string l_vorticity_x;
  string l_vorticity_y;
  string l_vorticity_z;
  string l_invariantQ;
  string l_totalp;
  string l_velocity_x;
  string l_velocity_y;
  string l_velocity_z;
  string l_fvelocity_x;
  string l_fvelocity_y;
  string l_fvelocity_z;
  string l_pressure;
  string l_temperature;
  string l_avr_pressure;
  string l_avr_temperature;
  string l_avr_velocity_x;
  string l_avr_velocity_y;
  string l_avr_velocity_z;
  string l_rmsV_x;
  string l_rmsV_y;
  string l_rmsV_z;
  string l_rmsmeanV_x;
  string l_rmsmeanV_y;
  string l_rmsmeanV_z;
  string l_rmsP;
  string l_rmsmeanP;
  string l_rmsT;
  string l_rmsmeanT;
  
  
public:
  
  PLT3D() {
    size_allocated = 0;
    size_OutBuffer = 0;
    size_InBuffer  = 0;
    Iblank   = 0;
    XYZfile  = 0;
    
    // ファイル入出力
    DFI_IN_INS   = NULL;
    DFI_IN_STAT   = NULL;
    DFI_OUT_INS  = NULL;
    DFI_OUT_STAT  = NULL;
    
    // ファイル名
    f_dfi_out_ins = "field";
    f_dfi_out_stat= "field_stat";
    
    l_divergence  = "Divergence_V";
    l_helicity    = "Helicity";
    l_vorticity_x = "Vorticity_X";
    l_vorticity_y = "Vorticity_Y";
    l_vorticity_z = "Vorticity_Z";
    l_invariantQ  = "Invariant_Q";
    l_totalp      = "TotalPressure";
    l_velocity_x  = "Velocity_X";
    l_velocity_y  = "Velocity_Y";
    l_velocity_z  = "Velocity_Z";
    l_fvelocity_x = "Face_Velocity_X";
    l_fvelocity_y = "Face_Velocity_Y";
    l_fvelocity_z = "Face_Velocity_Z";
    l_pressure    = "Pressure";
    l_temperature = "Temperature";
    l_avr_pressure    = "Avr_Pressure";
    l_avr_temperature = "Avr_Temperature";
    l_avr_velocity_x  = "Avr_Velocity_X";
    l_avr_velocity_y  = "Avr_Velocity_Y";
    l_avr_velocity_z  = "Avr_Velocity_Z";
    l_rmsV_x     = "RMS_Velocity_X";
    l_rmsV_y     = "RMS_Velocity_Y";
    l_rmsV_z     = "RMS_Velocity_Z";
    l_rmsmeanV_x = "RMS_mean_Velocity_X";
    l_rmsmeanV_y = "RMS_mean_Velocity_Y";
    l_rmsmeanV_z = "RMS_mean_Velocity_Z";
    l_rmsP       = "RMS_Pressure";
    l_rmsmeanP   = "RMS_mean_Pressure";
    l_rmsT       = "RMS_Temperature";
    l_rmsmeanT   = "RMS_mean_Temperature";
  }
  
  ~PLT3D() {
    if( DFI_IN_INS   != NULL ) delete DFI_IN_INS;
    if( DFI_IN_STAT  != NULL ) delete DFI_IN_STAT;
    if( DFI_OUT_INS  != NULL ) delete DFI_OUT_INS;
    if( DFI_OUT_STAT != NULL ) delete DFI_OUT_STAT;
  }
  
 
private:
  
  
  // 固有のオプション
  virtual void getInherentOption();
  
  
  /*
   * @brief リスタート時の瞬時値ファイル読み込み
   * @param [in]  fp             ファイルポインタ
   * @param [out] m_CurrentStep  CurrentStep
   * @param [out] m_CurrentTime  CurrentTime
   * @param [out] flop           浮動小数点演算数
   */
  virtual void RestartInstantaneous(FILE* fp,
                                    unsigned& m_CurrentStep,
                                    double& m_CurrentTime,
                                    double& flop);
  
  
  // 固有パラメータの表示
  virtual void printSteerConditionsInherent(FILE* fp);
  
  
public:
  
  // Iblankの利用の有無を返す
  int getIblank() const
  {
    return Iblank;
  }
  
  
  // リスタートに必要なDFIファイルを取得
  virtual void getRestartDFI();
  
  
  /**
   * @brief ファイル出力の初期化
   * @param [in] id_cell   CellID
   * @param [in] id_bcf    BCflagID
   */
  virtual void initFileOut(const int id_cell, const int id_bcf);
  
  
  /**
   * @brief 時間統計値のファイル出力
   * @param [in]     m_CurrentStep     CurrentStep
   * @param [in]     m_CurrentTime     CurrentTime
   * @param [in]     m_CurrentStepStat CurrentStepStat
   * @param [in]     m_CurrentTimeStat CurrentTimeStat
   * @param [in,out] flop              浮動小数点演算数
   */
  virtual void OutputStatisticalVarables(const unsigned m_CurrentStep,
                                         const double m_CurrentTime,
                                         const unsigned m_CurrentStepStat,
                                         const double m_CurrentTimeStat,
                                         double& flop);
  
  
  /**
   * @brief 基本変数のファイル出力
   * @param [in]     m_CurrentStep     CurrentStep
   * @param [in]     m_CurrentTime     CurrentTime
   * @param [in,out] flop              浮動小数点演算数
   */
  virtual void OutputBasicVariables(const unsigned m_CurrentStep,
                                    const double m_CurrentTime,
                                    double& flop);
  
  
  /**
   * @brief リスタートプロセス
   * @param [in]     fp                ファイルポインタ
   * @param [out]    m_CurrentStep     CurrentStep
   * @param [out]    m_CurrentTime     CurrentTime
   */
  virtual void Restart(FILE* fp,
                       unsigned& m_CurrentStep,
                       double& m_CurrentTime);
  
  
  /**
   * @brief リスタート時の統計値ファイル読み込み
   * @param [in]  fp                ファイルポインタ
   * @param [in]  m_CurrentStep     CurrentStep
   * @param [in]  m_CurrentTime     CurrentTime
   * @param [out] m_CurrentStepStat CurrentStepStat
   * @param [out] m_CurrentTimeStat CurrentTimeStat
   * @param [out] flop              浮動小数点演算数
   */
  virtual void RestartStatistic(FILE* fp,
                                const unsigned m_CurrentStep,
                                const double m_CurrentTime,
                                unsigned& m_CurrentStepStat,
                                double& m_CurrentTimeStat,
                                double& flop);

};

#endif // _FFV_PLT3D_H_

