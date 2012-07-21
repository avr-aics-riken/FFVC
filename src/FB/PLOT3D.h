#ifndef _FB_FILE_IO_PLOT3D_H_
#define _FB_FILE_IO_PLOT3D_H_

// #################################################################
//
// FFV : Frontflow / violet
//
// Copyright (c) All right reserved. 2012
//
// Institute of Industrial Science, The University of Tokyo, Japan. 
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

//#ifndef _WIN32
//#include <unistd.h>
//#include <strings.h>
//#else
//#include "sph_win32_util.h"
//#endif
//#include <sys/types.h>
//
//#if defined(IA32_LINUX) || defined(IA64_LINUX) || defined(SGI_ALTIX)
//#include <sys/stat.h>
//#endif

//#ifdef MacOSX
//#include <sys/uio.h>
//#endif
//
//#ifdef REAL_TYPE_DOUBLE
//#define REAL_TYPE float
//#else
//#define REAL_TYPE double
//#endif

// general
//#define FB_FILE_PATH_LENGTH 64
#define FB_BUFF_LENGTH      80


#ifdef _WIN32

//FileIO_PLOT3D_fortran.f90
#define open_plot3d_file_       OPEN_PLOT3D_FILE
#define open_plot3d_outputfile_ OPEN_PLOT3D_OUTPUTFILE
#define open_plot3d_inputfile_  OPEN_PLOT3D_INPUTFILE
#define close_plot3d_file_      CLOSE_PLOT3D_FILE
#define write_line_             WRITE_FUNCTION_NAME
#define read_line_              READ_FUNCTION_NAME
#define write_fvbnd_boundary_   WRITE_FVBND_BOUNDARY
#define read_fvbnd_boundary_    READ_FVBND_BOUNDARY

#endif // _WIN32


extern "C" {

//FileIO_PLOT3D_fortran.f90
  void open_plot3d_file_ (int* iflag, char* fname, int* ifl, int* fnsize, int* ierror);
  void open_plot3d_outputfile_ (int* iflag, char* fname, int* ifl, int* fnsize, int* ierror);
  void open_plot3d_inputfile_ (int* iflag, char* fname, int* ifl, int* fnsize, int* ierror);
  void close_plot3d_file_ (int* ifl);
  void write_line_ (char* buff, int* ifl, int* fnsize);
  void read_line_ (char* buff, int* ifl, int* fnsize, int* iend);
  void write_fvbnd_boundary_ (
	  int* tp, int* gn, int* imin, int* imsx, int* jmin, int* jmax, int* kmin, int* kmax, char* flag, int* dir, int* ifl);
  void read_fvbnd_boundary_ (
	  int* tp, int* gn, int* imin, int* imsx, int* jmin, int* jmax, int* kmin, int* kmax, char* flag, int* dir, int* ifl);
}


using namespace std;

//class FileIO_PLOT3D {
class FileIO_PLOT3D : public DomainInfo {

public:
  /** コンストラクタ */
  FileIO_PLOT3D() 
  {
    //paraMngr = NULL;

    P3Op.GridKind   = SINGLE_GRID;
    P3Op.MoveGrid   = GRID_NOT_MOVE;
    P3Op.Steady     = FB_UNSTEADY;
    P3Op.IBlankFlag = NOT_SET_IBLANK;
    P3Op.DimIs      = DIMENSION_3D;
    P3Op.Format     = UNFORMATTED; 

    fp_xyz = NULL;
    fp_q = NULL;
	ifl = 31;
  }

  /**　デストラクタ */
  ~FileIO_PLOT3D() {}
  
//private:
  string fname;

  FILE* fp_xyz;
  FILE* fp_q;
  int ifl;

  ////データ変数ポインタ
  int id;
  int jd;
  int kd;
  int ngrid;
  REAL_TYPE* x;
  REAL_TYPE* y;
  REAL_TYPE* z;
  int* iblank;
  REAL_TYPE fsmach;
  REAL_TYPE alpha;
  REAL_TYPE re;
  REAL_TYPE time;
  REAL_TYPE* q;
  int nvar;
  REAL_TYPE* d;
  string funcname;

protected:
//  cpm_ParaManager *paraMngr; ///< Cartesian Partition Maneger

  // オプション
  typedef struct {
    int GridKind;    ///=0:single grid, =1:multi grid
    int MoveGrid;    ///=0:not move, =1:move
    int Steady;      ///=0:Steady, =1:Unsteady
    int IBlankFlag;  ///=0:not set IBlank, =1 :set IBlank
    int DimIs;       ///=2:2D, =3:3D
    int Format;      ///
  } Plot3D_Option;

  //// GRIDの種類
  //enum PLOT3D_GridKind {
  //  SINGLE_GRID,
  //  MULTI_GRID
  //};

  //// *.xyzファイルの出力を時刻暦にするかどうか
  //enum PLOT3D_MoveGrid {
  //  GRID_NOT_MOVE,
  //  GRID_MOVE
  //};

  //// 定常or非定常解析 *.qファイルの出力
  //enum PLOT3D_Steady {
  //  STEADY,
  //  UNSTEADY
  //};

  //// IBLANKフラグのセットの有無
  //enum PLOT3D_IBlankFlag {
  //  NOT_SET_IBLANK,
  //  SET_IBLANK
  //};

  //// DIMENSITON
  //enum PLOT3D_DimIs {
  //  DIMENSION2D=2,
  //  DIMENSION3D
  //};

  //// PLOT3D File Format
  //enum PLOT3D_Fiel_Format {
  //  UNFORMATTED=1,
  //  FORMATTED,
  //  UNFORMATTED_SPECIAL
  //};

  Plot3D_Option P3Op;

public:
  
  /** CPMlibのポインタをセット 
   * @param [in] m_paraMngr  CPMクラスのポインタ
   */
  //void importCPM(cpm_ParaManager* m_paraMngr);

  //@fn int IsGridKind(void) const
  //@brief return 0 : single grid, return 1 : multi grid
  int IsGridKind(void) const {
    return P3Op.GridKind;
  }

  //@fn int IsMoveGrid(void) const
  //@brief return 0 : not move, return 1 : move
  int IsMoveGrid(void) const {
    return P3Op.MoveGrid;
  }

  //@fn int IsSteady(void) const
  //@brief return 0 : Steady, return 1 : Unsteady
  int IsSteady(void) const {
    return P3Op.Steady;
  }

  //@fn int IsIBlankFlag(void) const
  //@brief return 0 : not set IBlank, return 1 : set IBlank
  int IsIBlankFlag(void) const {
    return P3Op.IBlankFlag;
  }

  //@fn int GetDim(void) const
  //@brief 
  int GetDim(void) const { 
    return P3Op.DimIs; ///=2:2D, =3:3D
  }

  //@fn int GetFormat(void) const
  //@brief 
  int GetFormat(void) const { 
    return P3Op.Format; 
  }

  /**
   * @brief MoveGridを設定する
   * @retval 設定の成否
   * @param is フラグ
   */
  bool setMoveGrid(const int is);

  /**
   * @brief Steadyを設定する
   * @retval 設定の成否
   * @param is フラグ
   */
  bool setSteady(const int is);

  /**
   * @brief IBlankFlagを設定する
   * @retval 設定の成否
   * @param is フラグ
   */
  bool setIBlankFlag(const int is);

  /**
   * @brief 次元を2Dに設定する
   */
  void setDimension2D(){ P3Op.DimIs=DIMENSION_2D; };

  /**
   * @brief 次元を3Dに設定する
   */
  void setDimension3D(){ P3Op.DimIs=DIMENSION_3D; };

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
   * @retval 設定の成否
   * @param is フラグ
   */
  bool setFormat(const int is);

  /**
   * @brief IBlankFlagを設定する
   * @retval 設定の成否
   * @param is フラグ
   */
  void setFilePortNumver(const int is){ ifl=is; };

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
  REAL_TYPE* m_x,
  REAL_TYPE* m_y,
  REAL_TYPE* m_z,
  int* m_iblank);

  /**
   * @brief 形状データのセット（2D）
   */
  void setXYZData2D(
  REAL_TYPE* m_x,
  REAL_TYPE* m_y,
  int* m_iblank);

  /**
   * @brief 結果Qデータのセット
   */
  void setQData(
  REAL_TYPE m_fsmach,
  REAL_TYPE m_alpha,
  REAL_TYPE m_re,
  REAL_TYPE m_time,
  REAL_TYPE* m_q);

  /**
   * @brief 結果Funcデータのセット
   */
  void setFuncData(
  int m_nvar,
  REAL_TYPE* m_d);

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

#endif // _FB_FILE_IO_PLOT3D_H_
