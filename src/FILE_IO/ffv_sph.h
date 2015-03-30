#ifndef _FFV_SPH_H_
#define _FFV_SPH_H_

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
//

/**
 * @file   ffv_sph.h
 * @brief  File IO of SPH Class Header
 * @author aics
 */


#include "ffv_io_base.h"


class SPH : public IO_BASE {
  
private:
  // 入力dfiファイルのプレフィックス
  string f_dfi_in_prs;
  string f_dfi_in_vel;
  string f_dfi_in_fvel;
  string f_dfi_in_temp;
  string f_dfi_in_prsa;
  string f_dfi_in_vela;
  string f_dfi_in_tempa;
  
  // 出力ファイルのプレフィックス
  string f_Velocity;
  string f_Pressure;
  string f_Temperature;
  string f_AvrPressure;
  string f_AvrVelocity;
  string f_AvrTemperature;
  string f_DivDebug;
  string f_Helicity;
  string f_TotalP;
  string f_I2VGT;
  string f_Vorticity;
  string f_Fvelocity;
  
  // pointers to CDM class
  
  // InFile
  cdm_DFI *DFI_IN_PRS;      ///< Pressure
  cdm_DFI *DFI_IN_VEL;      ///< Velocity
  cdm_DFI *DFI_IN_FVEL;     ///< Face velocity
  cdm_DFI *DFI_IN_TEMP;     ///< Temperature
  cdm_DFI *DFI_IN_PRSA;     ///< Averaged pressure
  cdm_DFI *DFI_IN_VELA;     ///< Averaged velocity
  cdm_DFI *DFI_IN_TEMPA;    ///< Averaged temperature
  
  // OutFile
  cdm_DFI *DFI_OUT_PRS;     ///< Pressure
  cdm_DFI *DFI_OUT_VEL;     ///< Velocity
  cdm_DFI *DFI_OUT_FVEL;    ///< Face velocity
  cdm_DFI *DFI_OUT_TEMP;    ///< Temperature
  cdm_DFI *DFI_OUT_PRSA;    ///< Averaged Pressure
  cdm_DFI *DFI_OUT_VELA;    ///< Averaged velocity
  cdm_DFI *DFI_OUT_TEMPA;   ///< Averaged temperature
  cdm_DFI *DFI_OUT_TP;      ///< Total Pressure
  cdm_DFI *DFI_OUT_VRT;     ///< Vorticity
  cdm_DFI *DFI_OUT_I2VGT;   ///< 2nd Invariant of Velocity Gradient Tensor
  cdm_DFI *DFI_OUT_HLT;     ///< Helicity
  cdm_DFI *DFI_OUT_DIV;     ///< Divergence for debug
  
  
public:
  
  SPH() {
    // ファイル入出力
    DFI_IN_PRS   = NULL;
    DFI_IN_VEL   = NULL;
    DFI_IN_FVEL  = NULL;
    DFI_IN_TEMP  = NULL;
    DFI_IN_PRSA  = NULL;
    DFI_IN_VELA  = NULL;
    DFI_IN_TEMPA = NULL;
    DFI_OUT_PRS  = NULL;
    DFI_OUT_VEL  = NULL;
    DFI_OUT_FVEL = NULL;
    DFI_OUT_TEMP = NULL;
    DFI_OUT_PRSA = NULL;
    DFI_OUT_VELA = NULL;
    DFI_OUT_TEMPA= NULL;
    DFI_OUT_TP   = NULL;
    DFI_OUT_VRT  = NULL;
    DFI_OUT_I2VGT= NULL;
    DFI_OUT_HLT  = NULL;
    DFI_OUT_DIV  = NULL;
    
    // ファイル名
    f_Pressure       = "prs";
    f_Velocity       = "vel";
    f_Fvelocity      = "fvel";
    f_Temperature    = "tmp";
    f_AvrPressure    = "prsa";
    f_AvrVelocity    = "vela";
    f_AvrTemperature = "tmpa";
    f_DivDebug       = "div";
    f_Helicity       = "hlt";
    f_TotalP         = "tp";
    f_I2VGT          = "qcr";
    f_Vorticity      = "vrt";
  }
  
  ~SPH() {
    if( DFI_IN_PRS    != NULL ) delete DFI_IN_PRS;
    if( DFI_IN_VEL    != NULL ) delete DFI_IN_VEL;
    if( DFI_IN_FVEL   != NULL ) delete DFI_IN_FVEL;
    if( DFI_IN_TEMP   != NULL ) delete DFI_IN_TEMP;
    if( DFI_IN_PRSA   != NULL ) delete DFI_IN_PRSA;
    if( DFI_IN_VELA   != NULL ) delete DFI_IN_VELA;
    if( DFI_IN_TEMPA  != NULL ) delete DFI_IN_TEMPA;
    if( DFI_OUT_PRS   != NULL ) delete DFI_OUT_PRS;
    if( DFI_OUT_VEL   != NULL ) delete DFI_OUT_VEL;
    if( DFI_OUT_FVEL  != NULL ) delete DFI_OUT_FVEL;
    if( DFI_OUT_TEMP  != NULL ) delete DFI_OUT_TEMP;
    if( DFI_OUT_PRSA  != NULL ) delete DFI_OUT_PRSA;
    if( DFI_OUT_VELA  != NULL ) delete DFI_OUT_VELA;
    if( DFI_OUT_TEMPA != NULL ) delete DFI_OUT_TEMPA;
    if( DFI_OUT_TP    != NULL ) delete DFI_OUT_TP;
    if( DFI_OUT_VRT   != NULL ) delete DFI_OUT_VRT;
    if( DFI_OUT_I2VGT != NULL ) delete DFI_OUT_I2VGT;
    if( DFI_OUT_HLT   != NULL ) delete DFI_OUT_HLT;
    if( DFI_OUT_DIV   != NULL ) delete DFI_OUT_DIV;
  }
  
  
private:
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
  
  
  /**
   * @brief ファイル出力の初期化
   * @param [in] id_cell   CellID
   * @param [in] id_bcf    BCflagID
   */
  virtual void initFileOut(const int id_cell, const int id_bcf);
  

  /**
   * @brief 時間平均値のファイル出力
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
   * @brief リスタート時の平均値ファイル読み込み
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
  

  /**
   * @brief // チャネル乱流統計量の出力
   * @param [in]     d_av              速度 (時間平均値)            
   * @param [in]     d_rms_mean_v      乱流強度 (時間平均値)            
   * @param [in]     d_aR              レイノルズ応力 (時間平均値)            
   * @param [in]     d_aP              生成項 (時間平均値)            
   * @param [in]     d_aE              散逸項 (時間平均値)            
   * @param [in]     d_aT              乱流拡散項 (時間平均値)            
   * @param [in]     d_aPI             速度圧力勾配相関項 (時間平均値)            
   * @param [in]     myRank            自ノードのランク番号 
   * @param [in]     sz                領域サイズ (時間平均値)            
   * @param [in]     CurrentStepStat   統計操作の積算ステップ数 (時間平均値)            
   * @param [in]     dh                格子幅 (時間平均値)            
   * @param [in]     g                 ガイドセル長 (時間平均値)            
   * @param [in,out] flop              浮動小数点演算数
   */
  virtual void OutputMean(REAL_TYPE*         d_av,
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
                          double&            flop);
};

#endif // _FFV_SPH_H_
