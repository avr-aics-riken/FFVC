#ifndef _PLT3D_PLOT3D_H_
#define _PLT3D_PLOT3D_H_

// #################################################################
//
// FFV : Frontflow / violet
//
// Copyright (c) 2012 All right reserved.
//
// Institute of Industrial Science, University of Tokyo, Japan. 
//
// #################################################################

/** 
 * @file   PLOT3D.h
 * @brief  FlowBase FileIO_PLOT3D class Header
 * @author kero
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <time.h>
#include <fcntl.h>
#include <string>

#include "cpm_Define.h"
#include "cpm_ParaManager.h"

#include "DomainInfo.h"
#include "FB_Define.h"
#include "FB_Ffunc.h"
#include "mydebug.h"


#ifdef _WIN32

// PLT3D.f90
#define open_plot3d_file_       OPEN_PLOT3D_FILE
#define open_plot3d_outputfile_ OPEN_PLOT3D_OUTPUTFILE
#define open_plot3d_inputfile_  OPEN_PLOT3D_INPUTFILE
#define close_plot3d_file_      CLOSE_PLOT3D_FILE
#define write_line_             WRITE_LINE
#define read_line_              READ_LINE
#define write_fvbnd_boundary_   WRITE_FVBND_BOUNDARY
#define read_fvbnd_boundary_    READ_FVBND_BOUNDARY

#endif // _WIN32


extern "C" {

// PLT3D.f90
  void open_plot3d_file_ (int* iflag,
                          char* fname,
                          int* ifl,
                          int* fnsize,
                          int* ierror);
  void open_plot3d_outputfile_ (int* iflag,
                                char* fname,
                                int* ifl,
                                int* fnsize,
                                int* ierror);
  void open_plot3d_inputfile_ (int* iflag,
                               char* fname,
                               int* ifl,
                               int* fnsize,
                               int* ierror);
  void close_plot3d_file_ (int* ifl);
  void write_line_ (char* buff,
                    int* ifl,
                    int* fnsize);
  void read_line_ (char* buff,
                   int* ifl,
                   int* fnsize,
                   int* iend);
  void write_fvbnd_boundary_ (int* tp,
                              int* gn,
                              int* imin,
                              int* imsx,
                              int* jmin,
                              int* jmax,
                              int* kmin,
                              int* kmax,
                              char* flag,
                              int* dir,
                              int* ifl);
  void read_fvbnd_boundary_ (int* tp,
                             int* gn,
                             int* imin,
                             int* imsx,
                             int* jmin,
                             int* jmax,
                             int* kmin,
                             int* kmax,
                             char* flag,
                             int* dir,
                             int* ifl);
}


using namespace std;


class FileIO_PLOT3D : public DomainInfo {

public:
  /** コンストラクタ */
  FileIO_PLOT3D() 
  {

    P3Op.GridKind   = SINGLE_GRID;
    P3Op.MoveGrid   = GRID_NOT_MOVE;
    P3Op.Steady     = FB_UNSTEADY;
    P3Op.IBlankFlag = NOT_SET_IBLANK;
    P3Op.DimIs      = DIMENSION_3D;
    P3Op.Format     = UNFORMATTED; 
    P3Op.realtype   = OUTPUT_FLOAT; 

    fp = NULL;
	  ifl = 31;
  }

  /**　デストラクタ */
  ~FileIO_PLOT3D() {}
  
//private:
  string fname;

  FILE* fp;
  int ifl;

  int id;
  int jd;
  int kd;
  int ngrid;
  int* iblank;
  int nvar;
  string funcname;

  float* x;
  float* y;
  float* z;
  float fsmach;
  float alpha;
  float re;
  float time;
  float* q;
  float* d;

  double* dx;
  double* dy;
  double* dz;
  double dfsmach;
  double dalpha;
  double dre;
  double dtime;
  double* dq;
  double* dd;

protected:
  
  // オプション
  typedef struct {
    int GridKind;    ///=0:single grid, =1:multi grid
    int MoveGrid;    ///=0:not move, =1:move
    int Steady;      ///=0:Steady, =1:Unsteady
    int IBlankFlag;  ///=0:not set IBlank, =1 :set IBlank
    int DimIs;       ///=2:2D, =3:3D
    int Format;      ///
    int realtype;    ///=1:float,2=double output flag
  } Plot3D_Option;

  Plot3D_Option P3Op;

public:

  //@fn int IsGridKind(void) const
  //@brief return 0 : single grid, return 1 : multi grid
  int IsGridKind(void) const {
    return P3Op.GridKind;
  };

  //@fn int IsMoveGrid(void) const
  //@brief return 0 : not move, return 1 : move
  int IsMoveGrid(void) const {
    return P3Op.MoveGrid;
  };

  //@fn int IsSteady(void) const
  //@brief return 0 : Steady, return 1 : Unsteady
  int IsSteady(void) const {
    return P3Op.Steady;
  };

  //@fn int IsIBlankFlag(void) const
  //@brief return 0 : not set IBlank, return 1 : set IBlank
  int IsIBlankFlag(void) const {
    return P3Op.IBlankFlag;
  };

  //@fn int GetDim(void) const
  //@brief 
  int GetDim(void) const { 
    return P3Op.DimIs; ///=2:2D, =3:3D
  };

  //@fn int GetFormat(void) const
  //@brief 
  int GetFormat(void) const { 
    return P3Op.Format; 
  };

  //@fn int GetRealType(void) const
  //@brief 
  int GetRealType(void) const { 
    return P3Op.realtype; 
  };

  /**
   * @brief MoveGridを設定する
   * @retval 設定の成否
   * @param is フラグ
   */
  void setMoveGrid(const int is);

  /**
   * @brief Steadyを設定する
   * @retval 設定の成否
   * @param is フラグ
   */
  void setSteady(const int is);

  /**
   * @brief IBlankFlagを設定する
   * @retval 設定の成否
   * @param is フラグ
   */
  void setIBlankFlag(const int is);

  /**
   * @brief 次元を設定する
   */
  void setDim(const int is){ P3Op.DimIs=is; };

  /**
   * @brief 次元を2Dに設定する
   */
  void setDimension2D(){ P3Op.DimIs=DIMENSION_2D; };

  /**
   * @brief 次元を3Dに設定する
   */
  void setDimension3D(){ P3Op.DimIs=DIMENSION_3D; };

  /**
   * @brief GridKindを設定する
   */
  void setGridKind(const int is){ P3Op.GridKind=is; };

  /**
   * @brief SINGLE_GRIDに設定する
   */
  void setSingleGrid(){ P3Op.GridKind=SINGLE_GRID; };

  /**
   * @brief MULTI_GRIDに設定する
   */
  void setMultiGrid(){ P3Op.GridKind=MULTI_GRID; };

  /**
   * @brief PLOT3Dのファイルフォーマットを設定する
   * @param is フラグ
   */
  void setFormat(const int is);

  /**
   * @brief fortarn出力時のファイル装置番号のセット
   * @param is ファイル装置番号
   */
  void setFilePortNumber(const int is){ ifl=is; };

  /**
   * @brief 出力時の単精度or倍精度フラグ
   * @param is フラグ
   */
  void setRealType(const int is){ P3Op.realtype=is; };

  /**
   * @brief Gridデータのセット
   */
  void setGridData(
  int m_id,
  int m_jd,
  int m_kd,
  int m_ngrid);

  /**
   * @brief Gridデータのセット（2D）
   */
  void setGridData2D(
  int m_id,
  int m_jd,
  int m_ngrid);

  /**
   * @brief 形状データのセット
   */
  void setXYZData(
  float* m_x,
  float* m_y,
  float* m_z,
  int* m_iblank);

  /**
   * @brief 形状データのセット
   */
  void setXYZData(
  double* m_x,
  double* m_y,
  double* m_z,
  int* m_iblank);

  /**
   * @brief 形状データのセット（2D）
   */
  void setXYZData2D(
  float* m_x,
  float* m_y,
  int* m_iblank);

  /**
   * @brief 形状データのセット（2D）
   */
  void setXYZData2D(
  double* m_x,
  double* m_y,
  int* m_iblank);

  /**
   * @brief 結果Qデータのセット
   */
  void setQData(
  float m_fsmach,
  float m_alpha,
  float m_re,
  float m_time,
  float* m_q);

  /**
   * @brief 結果Qデータのセット
   */
  void setQData(
  double m_fsmach,
  double m_alpha,
  double m_re,
  double m_time,
  double* m_q);

  /**
   * @brief 結果Funcデータ項目数のセット
   */
  void setFuncDataNum(
  int m_nvar);

  /**
   * @brief 結果Funcデータのセット
   */
  void setFuncData(
  float* m_d);

  /**
   * @brief 結果Funcデータのセット
   */
  void setFuncData(
  double* m_d);


  /**
   * @brief PLOT3Dのファイル名保持
   */
  void setFileName(const char* tmp){ fname=tmp; };

  /**
   * @brief PLOT3Dのファイルオープン
   */
  virtual bool OpenFile();

  /**
   * @brief PLOT3Dのファイルクローズ
   */
  void CloseFile();

};

#endif // _PLT3D_PLOT3D_H_
