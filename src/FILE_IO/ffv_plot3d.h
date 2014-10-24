#ifndef _FFV_PLT3D_H_
#define _FFV_PLT3D_H_

//##################################################################################
//
// FFV-C : Frontflow / violet Cartesian
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
//

/**
 * @file   ffv_plot3d.h
 * @brief  File IO of PLOT3D Class Header
 * @author aics
 */


#include "ffv_io_base.h"


class PLT3D : public IO_BASE {
  
private:
  // 入力dfiファイルのプレフィックス
  string f_dfi_in_ins;
  string f_dfi_in_avr;

  // 出力ファイルのプレフィックス
  string f_dfi_out_ins;
  string f_dfi_out_avr;
  
  // InFile
  cdm_DFI *DFI_IN_INS;      ///< Instantaneus field
  cdm_DFI *DFI_IN_AVR;      ///< Averaged field
  
  // OutFile
  cdm_DFI *DFI_OUT_INS;     ///< Instantaneus field
  cdm_DFI *DFI_OUT_AVR;     ///< Averaged field
  
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
  
  
public:
  
  PLT3D() {
    // ファイル入出力
    DFI_IN_INS   = NULL;
    DFI_IN_AVR   = NULL;
    DFI_OUT_INS  = NULL;
    DFI_OUT_AVR  = NULL;
    
    // ファイル名
    f_dfi_out_ins = "field";
    f_dfi_out_avr = "field_avr";
    
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
  }
  
  ~PLT3D() {
    if( DFI_IN_INS    != NULL ) delete DFI_IN_INS;
    if( DFI_IN_AVR    != NULL ) delete DFI_IN_AVR;
    if( DFI_OUT_INS   != NULL ) delete DFI_OUT_INS;
    if( DFI_OUT_AVR   != NULL ) delete DFI_OUT_AVR;
  }
  
 
private:
  
  /**
   * @brief スカラデータをバッファから抜き出す
   * @param [in]  blk  先頭ブロック数
   * @apram [out] vec  スカラ配列
   * @retval 次の先頭ブロック数
   */
  int extract_scalar(int blk, REAL_TYPE* scl);
  
  
  /**
   * @brief ベクトルデータをバッファから抜き出す
   * @param [in]  blk  先頭ブロック数
   * @apram [out] vec  ベクトル配列
   * @retval 次の先頭ブロック数
   */
  int extract_vector(int blk, REAL_TYPE* vec);
  
  
  // @brief IBLANK ファイルを作成
  void generateIBLANK();
  
  
  /**
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
  
  
public:
  
  // リスタートに必要なDFIファイルを取得
  virtual void getRestartDFI();
  
  
  // @brief ファイル出力の初期化
  virtual void initFileOut();
  
  
  /**
   * @brief 時間平均値のファイル出力
   * @param [in]     m_CurrentStep     CurrentStep
   * @param [in]     m_CurrentTime     CurrentTime
   * @param [in]     m_CurrentStepAvr  CurrentStepAvr
   * @param [in]     m_CurrentTimeAvr  CurrentTimeAvr
   * @param [in,out] flop              浮動小数点演算数
   */
  virtual void OutputAveragedVarables(const unsigned m_CurrentStep,
                                      const double m_CurrentTime,
                                      const unsigned m_CurrentStepAvr,
                                      const double m_CurrentTimeAvr,
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
   * @brief リスタート時の平均値ファイル読み込み
   * @param [in]  fp                ファイルポインタ
   * @param [in]  m_CurrentStep     CurrentStep
   * @param [in]  m_CurrentTime     CurrentTime
   * @param [out] m_CurrentStepAvr  CurrentStepAvr
   * @param [out] m_CurrentTimeAvr  CurrentTimeAvr
   * @param [out] flop              浮動小数点演算数
   */
  virtual void RestartAvrerage(FILE* fp,
                               const unsigned m_CurrentStep,
                               const double m_CurrentTime,
                               unsigned& m_CurrentStepAvr,
                               double& m_CurrentTimeAvr,
                               double& flop);

};

#endif // _FFV_PLT3D_H_