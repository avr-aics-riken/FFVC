#ifndef _COMB_H_
#define _COMB_H_

// #################################################################
//
// Combine sph files and output 
//
// Copyright (c) All right reserved. 2012
//
// Institute of Industrial Science, The University of Tokyo, Japan. 
//
// #################################################################

/** 
 * @file   comb.h
 * @brief  COMB Class Header
 * @author kero
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <iostream>
#include <fstream>
#include <errno.h>

#include "mpi.h"

#ifndef _WIN32
#include <dirent.h>
#else
#include "sph_win32_util.h"   // for win32
#endif

#include "PerfMonitor.h"
#include "cpm_ParaManager.h"

//#include "DomainInfo.h"
#include "FB_Define.h"
#include "FBUtility.h"
#include "mydebug.h"
#include "TPControl.h"
#include "dfi.h"
#include "PLOT3D_read.h"
#include "PLOT3D_write.h"
#include "FileIO_sph.h"
//#include "omp.h"

#include "dfiinfo.h"

#include "limits.h" // for UBUNTU

//// FX10 profiler
//#if defined __K_FPCOLL
//#include "fjcoll.h"
//#elif defined __FX_FAPP
//#include "fj_tool/fjcoll.h"
//#endif

#include "COMB_Define.h"

using namespace std;


class COMB {

public:
  cpm_ParaManager* paraMngr; ///< Cartesian Partition Manager
  
public:
  int procGrp;         ///< プロセスグループ番号
  int myRank;          ///< 自ノードのランク番号
  int numProc;         ///< 全ランク数
  std::string HostName;  ///< ホスト名

  vector<int> index;

private:

  // dfi ファイル管理用 -> Kind_of_vars in FB_Define.h
  // 同じ解像度のリスタート時には、既にdfiファイルが存在する場合には、その内容を継続する
  // ラフリスタートの場合には、新規dfiファイルを生成する >> dfi.C
  int dfi_mng[var_END];

public:
  
  string filename;
  string dirname;
  int pflag;
  int pflagv;
  int lflag;
  int lflagv;
  bool thin_out;
  int thin_count;

  /** PLOT3D オプション */
  typedef struct 
  {
    string basename;
    int IS_xyz;
    int IS_q;
    int IS_funciton;
    int IS_function_name;
    int IS_fvbnd;
    int IS_DivideFunc;
    //Initializeでセット
    int ngrid; //出力ブロック数
    int nvar;  //出力項目数
  } Plot3D_Option;
  Plot3D_Option  P3Op;

  int output_real_type;
  int out_format;//combine sph or output plot3d
  int ndfi;//number of dfi file list
  vector<string> dfi_name;

  DFI DFI;                   ///< 分散ファイルインデクス管理クラス
  FileIO_PLOT3D_READ  FP3DR; ///< PLOT3D READクラス
  FileIO_PLOT3D_WRITE FP3DW; ///< PLOT3D WRITEクラス
  DfiInfo *DI;
  //FileIO_SPH FSPH;

public:
  /** コンストラクタ */
  COMB();

  /**　デストラクタ */
  ~COMB();
  
protected:

  FILE* fplog;

public:

  /**
   * @brief CPMのポインタをコピーし、ランク情報を設定
   * @param [in] m_paraMngr  cpm_ParaManagerクラス
   * @return  エラーコード
   */
  bool importCPM(cpm_ParaManager* m_paraMngr)
  {
    if ( !m_paraMngr ) return false;
    paraMngr = m_paraMngr;
    setRankInfo();
    return true;
  }

  /**
   * @brief ランク情報をセットする
   * @param [in] m_paraMngr  CPMlibポインタ
   * @param [in] m_proGrp    プロセスグループ番号
   */
  void setRankInfo()
  {
    myRank  = paraMngr->GetMyRankID();
    numProc = paraMngr->GetNumRank();
    HostName= paraMngr->GetHostName();
  }


  /**  
   * @brief 入力ファイルの読み込み
   * @param [in] fname  入力ファイル名
   */
  void ReadInit(string input_file);

  /**  
   * @brief 連結用入力ファイルの読み込み
   * @param [in] tpCntl  tpCntlクラスポインタ
   */
  void ReadInputFile(TPControl* tpCntl);
  
  /**
   * @brief PLOT3Dファイル入出力に関するパラメータ
   * @param [in] tpCntl  tpCntlクラスポインタ
   */
  void get_PLOT3D(TPControl* tpCntl);

  /**
   * @brief dfiファイルの読み込みとDfiInfoクラスデータの作成
   */
  void ReadDfiFiles();

  /**
   * @brief 出力指定ディレクトリのチェック
   */
  void CheckDir(string dirstr);

  /**
   * @brief ログファイルのオープン
   */
  void OpenLogFile();

  /**
   * @brief ログファイルのクローズ
   */
  void CloseLogFile();

  /**
   * @brief 所要時間の記述
   */
  void WriteTime(double* tt);

  /**
   * @brief sphファイルの読み込みとcombine sph or plot3d output
   */
  void CombineFiles();

  /**
   * @brief sphファイルのデータタイプの読み込み
   * @param[out] m_sv_type    データ種別
   * @param[out] m_d_type     データ型タイプ
   * @param[in] fp_in         FORTRAN入力用ファイル装置番号
   * @param[in] fname         ファイル名
   */
  bool ReadSphDataType(int* m_sv_type, int* m_d_type, int fp_in, string fname);

  /**
   * @brief sphファイルのheaderの読み込み（単精度）
   * @param[out] m_step       ステップ数
   * @param[out] m_sv_type    データ種別
   * @param[out] m_d_type     データ型タイプ
   * @param[out] m_imax       x方向ボクセルサイズ
   * @param[out] m_jmax       y方向ボクセルサイズ
   * @param[out] m_kmax       z方向ボクセルサイズ
   * @param[out] m_time       時間
   * @param[out] m_org        原点座標
   * @param[out] m_pit        ピッチ
   * @param[in] fp_in         FORTRAN入力用ファイル装置番号
   * @param[in] fname         ファイル名
   */
  bool ReadSphHeader(int* m_step,
                     int* m_sv_type,
                     int* m_d_type,
                     int* m_imax,
                     int* m_jmax,
                     int* m_kmax,
                     float* m_time,
                     float* m_org,
                     float* m_pit,
                     int fp_in,
                     string fname);

  /**
   * @brief sphファイルのheaderの読み込み（倍精度）
   * @param[out] m_step       ステップ数
   * @param[out] m_sv_type    データ種別
   * @param[out] m_d_type     データ型タイプ
   * @param[out] m_imax       x方向ボクセルサイズ
   * @param[out] m_jmax       y方向ボクセルサイズ
   * @param[out] m_kmax       z方向ボクセルサイズ
   * @param[out] m_time       時間
   * @param[out] m_org        原点座標
   * @param[out] m_pit        ピッチ
   * @param[in] fp_in         FORTRAN入力用ファイル装置番号
   * @param[in] fname         ファイル名
   */
  bool ReadSphHeader(long long* m_step,
                     int* m_sv_type,
                     int* m_d_type,
                     long long* m_imax,
                     long long* m_jmax,
                     long long* m_kmax,
                     double* m_time,
                     double* m_org,
                     double* m_pit,
                     int fp_in,
                     string fname);

  /**
   * @brief sphファイルのdataの読み込み（単精度）
   * @param[out] wk           データポインタ
   * @param[in] wksize        データサイズ
   * @param[in] dim           次元
   * @param[in] fp_in         FORTRAN入力用ファイル装置番号
   * @param[in] fname         ファイル名
   */
  bool ReadSphData(float* wk,
                   int wksize,
                   int* size,
                   int dim,
                   int fp_in,
                   string fname);

  /**
   * @brief sphファイルのdataの読み込み（倍精度）
   * @param[out] wk           データポインタ
   * @param[in] wksize        データサイズ
   * @param[in] dim           次元
   * @param[in] fp_in         FORTRAN入力用ファイル装置番号
   * @param[in] fname         ファイル名
   */
  bool ReadSphData(double* wk,
                   int wksize,
                   int* size,
                   int dim,
                   int fp_in,
                   string fname);


  /**
   * @brief 出力DFIファイル名を作成する
   * @param [in] prefix ファイル接頭文字
   */
  std::string Generate_DFI_Name(const std::string prefix);

  /**
   * @brief ファイル名を作成する
   * @param [in] prefix ファイル接頭文字
   * @param [in] m_step ステップ数
   * @param [in] m_id   ランク番号
   * @param [in] mio    出力時の分割指定　 true = local / false = gather(default)
   */
  std::string Generate_FileName(const std::string prefix, const unsigned m_step, const int m_id, const bool mio=false);
  

  /**
   * @brief ファイル名を作成する。（拡張子自由）
   * @param [in] prefix ファイル接頭文字
   * @param [in] xxx    拡張子
   * @param [in] m_step ステップ数
   * @param [in] m_id   ランク番号
   * @param [in] mio    出力時の分割指定　 true = local / false = gather(default)
   */
  std::string Generate_FileName_Free(const std::string prefix, const std::string xxx, const unsigned m_step, const int m_id, const bool mio=false);


  /**
   * @brief DFIファイルをコピーする
   * @param [in] base_name  コピー元ファイル
   * @param [out] new_name  コピー先ファイル
   */
  bool Copy_DFIFile(const std::string base_name, const std::string new_name, const std::string prefix, int& dfi_mng);

  
  /**
   * @brief Tab(space２つ)を出力する
   * @param [in] fp      ファイルポインタ
   * @param [in] tab     インデント数
   */
  void Write_Tab(FILE* fp, const unsigned tab);


  /**
   * @brief メモリ使用量を表示する
   * @param [in] Memory メモリ量
   * @param [in] fp     ファイルポインタ
   */
  void MemoryRequirement(const double Memory, FILE* fp);


  /**
   * @brief メモリ使用量を表示する
   * @param [in] TotalMemory トータルメモリ使用量最大値
   * @param [in] sphMemory sphファイル読み込みのためのwkメモリ使用量最大値
   * @param [in] plot3dMemory plot3dファイル書き込みのためのメモリ使用量最大値
   * @param [in] fp     ファイルポインタ
   */
  void MemoryRequirement(const double TotalMemory, const double sphMemory, const double plot3dMemory, FILE* fp);


  /**
   * @brief 現在のエンディアンをチェックする
   */
  int tt_check_machine_endian();


///////////////////////////////////////////////////////////////////////////////
// comb_sph.C

  /**
   * @brief sphファイルの連結
   */
  void output_sph();

#if 0

  /**
   * @brief sphファイルのheaderの書き込み（REAL_TYPE）
   * @param[in] step       ステップ数
   * @param[in] sv_type    データ種別
   * @param[in] d_type     データ型タイプ
   * @param[in] imax       x方向ボクセルサイズ
   * @param[in] jmax       y方向ボクセルサイズ
   * @param[in] kmax       z方向ボクセルサイズ
   * @param[in] time       時間
   * @param[in] org        原点座標
   * @param[in] pit        ピッチ
   * @param[in] fp         ファイルポインタ
   */
  bool WriteSphHeader(
    int step, int sv_type, int d_type, int imax, int jmax, int kmax,
    REAL_TYPE time, REAL_TYPE* org, REAL_TYPE* pit, FILE *fp);

#endif

  /**
   * @brief sphファイルのheaderの書き込み（単精度）
   * @param[in] step       ステップ数
   * @param[in] sv_type    データ種別
   * @param[in] d_type     データ型タイプ
   * @param[in] imax       x方向ボクセルサイズ
   * @param[in] jmax       y方向ボクセルサイズ
   * @param[in] kmax       z方向ボクセルサイズ
   * @param[in] time       時間
   * @param[in] org        原点座標
   * @param[in] pit        ピッチ
   * @param[in] fp         ファイルポインタ
   */
  bool WriteSphHeader(
    int step, int sv_type, int d_type, int imax, int jmax, int kmax,
    float time, float* org, float* pit, FILE *fp);

  /**
   * @brief sphファイルのheaderの書き込み（倍精度）
   * @param[in] step       ステップ数
   * @param[in] sv_type    データ種別
   * @param[in] d_type     データ型タイプ
   * @param[in] imax       x方向ボクセルサイズ
   * @param[in] jmax       y方向ボクセルサイズ
   * @param[in] kmax       z方向ボクセルサイズ
   * @param[in] time       時間
   * @param[in] org        原点座標
   * @param[in] pit        ピッチ
   * @param[in] fp         ファイルポインタ
   */
  bool WriteSphHeader(
    int step, int sv_type, int d_type, int imax, int jmax, int kmax,
    double time, double* org, double* pit, FILE *fp);

  /**
   * @brief マーカーの書き込み
   * @param[in] dmy         マーカー
   * @param[in] fp          ファイルポインタ
   */
  bool WriteCombineDataMarker(int dmy, FILE* fp);

  /**
   * @brief sphファイルのデータの書き込み（単精度）
   * @param[in] data       データ
   * @param[in] dLen       データサイズ
   * @param[in] fp         ファイルポインタ
   */
  bool WriteCombineData(float* data, size_t dLen, FILE* fp);

  /**
   * @brief sphファイルのデータの書き込み（倍精度）
   * @param[in] data       データ
   * @param[in] dLen       データサイズ
   * @param[in] fp         ファイルポインタ
   */
  bool WriteCombineData(double* data, size_t dLen, FILE* fp);


  /**
   * @brief sphファイルのデータの連結（float ---> float）
   * @param[in] d       連結した一層分データ
   * @param[in] size_d  dのデータサイズ
   * @param[in] s       連結する領域データ
   * @param[in] size_s  sのデータサイズ
   * @param[in] xsize   x方向サイズ
   * @param[in] ysize   y方向サイズ
   * @param[in] zsize   z方向サイズ
   * @param[in] z       z方向位置（層の高さ）
   * @param[in] sz      領域サイズ
   * @param[in] dim     =1:scalar　=3:vector
   * @param[in] xs      x方向位置
   * @param[in] ys      y方向位置
   */
  void CombineLayerData(
    float* d, int size_d, float* s, int size_s,
    int xsize, int ysize, int zsize,
    int z, int* sz, int dim, int xs, int ys); //z=zzz ; z!=zh


  /**
   * @brief sphファイルのデータの連結（double ---> double）
   * @param[in] d       連結した一層分データ
   * @param[in] size_d  dのデータサイズ
   * @param[in] s       連結する領域データ
   * @param[in] size_s  sのデータサイズ
   * @param[in] xsize   x方向サイズ
   * @param[in] ysize   y方向サイズ
   * @param[in] zsize   z方向サイズ
   * @param[in] z       z方向位置（層の高さ）
   * @param[in] sz      領域サイズ
   * @param[in] dim     =1:scalar　=3:vector
   * @param[in] xs      x方向位置
   * @param[in] ys      y方向位置
   */
  void CombineLayerData(
    double* d, int size_d, double* s, int size_s,
    int xsize, int ysize, int zsize,
    int z, int* sz, int dim, int xs, int ys); //z=zzz ; z!=zh


  /**
   * @brief sphファイルのデータの連結（float ---> double）
   * @param[in] d       連結した一層分データ
   * @param[in] size_d  dのデータサイズ
   * @param[in] s       連結する領域データ
   * @param[in] size_s  sのデータサイズ
   * @param[in] xsize   x方向サイズ
   * @param[in] ysize   y方向サイズ
   * @param[in] zsize   z方向サイズ
   * @param[in] z       z方向位置（層の高さ）
   * @param[in] sz      領域サイズ
   * @param[in] dim     =1:scalar　=3:vector
   * @param[in] xs      x方向位置
   * @param[in] ys      y方向位置
   */
  void CombineLayerData(
    double* d, int size_d, float* s, int size_s,
    int xsize, int ysize, int zsize,
    int z, int* sz, int dim, int xs, int ys); //z=zzz ; z!=zh


  /**
   * @brief sphファイルのデータの連結（double ---> float）
   * @param[in] d       連結した一層分データ
   * @param[in] size_d  dのデータサイズ
   * @param[in] s       連結する領域データ
   * @param[in] size_s  sのデータサイズ
   * @param[in] xsize   x方向サイズ
   * @param[in] ysize   y方向サイズ
   * @param[in] zsize   z方向サイズ
   * @param[in] z       z方向位置（層の高さ）
   * @param[in] sz      領域サイズ
   * @param[in] dim     =1:scalar　=3:vector
   * @param[in] xs      x方向位置
   * @param[in] ys      y方向位置
   */
  void CombineLayerData(
    float* d, int size_d, double* s, int size_s,
    int xsize, int ysize, int zsize,
    int z, int* sz, int dim, int xs, int ys); //z=zzz ; z!=zh


///////////////////////////////////////////////////////////////////////////////
// comb_plot3d.C

  /**
   * @brief PLOT3Dファイルの出力
   */
  void output_plot3d();

  /**
   * @brief xyzファイルの出力（sph:float,plot3d:float）
   * @param [in] m_step     ステップ 
   * @param [in] m_rank     ランク
   * @param [in] guide      ガイドセル数
   * @param [in] origin     基点座標
   * @param [in] pitch      ピッチ
   * @param [in] size       セルサイズ
   * @param [in] x          x方向座標ワーク
   * @param [in] y          y方向座標ワーク
   * @param [in] z          z方向座標ワーク
   */
  void OutputPlot3D_xyz(int m_step, int m_rank, int guide, float* origin, float* pitch, int* size, float* x, float* y, float* z);

  /**
   * @brief xyzファイルの出力（sph:double,plot3d:double）
   * @param [in] m_step     ステップ 
   * @param [in] m_rank     ランク
   * @param [in] guide      ガイドセル数
   * @param [in] origin     基点座標
   * @param [in] pitch      ピッチ
   * @param [in] size       セルサイズ
   * @param [in] x          x方向座標ワーク
   * @param [in] y          y方向座標ワーク
   * @param [in] z          z方向座標ワーク
   */
  void OutputPlot3D_xyz(int m_step, int m_rank, int guide, double* origin, double* pitch, int* size, double* x, double* y, double* z);

  /**
   * @brief xyzファイルの出力（sph:float,plot3d:double）
   * @param [in] m_step     ステップ 
   * @param [in] m_rank     ランク
   * @param [in] guide      ガイドセル数
   * @param [in] origin     基点座標
   * @param [in] pitch      ピッチ
   * @param [in] size       セルサイズ
   * @param [in] x          x方向座標ワーク
   * @param [in] y          y方向座標ワーク
   * @param [in] z          z方向座標ワーク
   */
  void OutputPlot3D_xyz(int m_step, int m_rank, int guide, float* origin, float* pitch, int* size, double* x, double* y, double* z);

  /**
   * @brief xyzファイルの出力（sph:double,plot3d:float）
   * @param [in] m_step     ステップ 
   * @param [in] m_rank     ランク
   * @param [in] guide      ガイドセル数
   * @param [in] origin     基点座標
   * @param [in] pitch      ピッチ
   * @param [in] size       セルサイズ
   * @param [in] x          x方向座標ワーク
   * @param [in] y          y方向座標ワーク
   * @param [in] z          z方向座標ワーク
   */
  void OutputPlot3D_xyz(int m_step, int m_rank, int guide, double* origin, double* pitch, int* size, float* x, float* y, float* z);


  /**
   * @brief Scalarの格子点での値をセット（ガイドセルの値は計算対象にいれていない）（float ---> float）
   * @param [out]    d        格子点data
   * @param [in]     data     セル中心data
   * @param [in]     id       セル中心data x方向サイズ
   * @param [in]     jd       セル中心data y方向サイズ
   * @param [in]     kd       セル中心data z方向サイズ
   */
  void setScalarGridData(int* size, int guide, float* d, float* data, int id, int jd, int kd);


  /**
   * @brief Vectorの格子点での値をセット（ガイドセルの値は計算対象にいれていない）（float ---> float）
   * @param [out]    d        格子点data
   * @param [in]     data     セル中心data
   * @param [in]     id       セル中心data x方向サイズ
   * @param [in]     jd       セル中心data y方向サイズ
   * @param [in]     kd       セル中心data z方向サイズ
   */
  void setVectorGridData(int* size, int guide, float* d, float* data, int id, int jd, int kd);


  /**
   * @brief 成分別Vectorの格子点での値をセット（ガイドセルの値は計算対象にいれていない）（float ---> float）
   * @param [out]    d        格子点data
   * @param [in]     data     セル中心data
   * @param [in]     id       セル中心data x方向サイズ
   * @param [in]     jd       セル中心data y方向サイズ
   * @param [in]     kd       セル中心data z方向サイズ
   * @param [in]     ivar     ベクトル成分 =0:x =1:y =2:z
   */
  void setVectorComponentGridData(int* size, int guide, float* d, float* data, int id, int jd, int kd, int ivar);


  /**
   * @brief Scalarの格子点での値をセット（ガイドセルの値は計算対象にいれていない）（double ---> double）
   * @param [out]    d        格子点data
   * @param [in]     data     セル中心data
   * @param [in]     id       セル中心data x方向サイズ
   * @param [in]     jd       セル中心data y方向サイズ
   * @param [in]     kd       セル中心data z方向サイズ
   */
  void setScalarGridData(int* size, int guide, double* d, double* data, int id, int jd, int kd);


  /**
   * @brief Vectorの格子点での値をセット（ガイドセルの値は計算対象にいれていない）（double ---> double）
   * @param [out]    d        格子点data
   * @param [in]     data     セル中心data
   * @param [in]     id       セル中心data x方向サイズ
   * @param [in]     jd       セル中心data y方向サイズ
   * @param [in]     kd       セル中心data z方向サイズ
   */
  void setVectorGridData(int* size, int guide, double* d, double* data, int id, int jd, int kd);


  /**
   * @brief 成分別Vectorの格子点での値をセット（ガイドセルの値は計算対象にいれていない）（double ---> double）
   * @param [out]    d        格子点data
   * @param [in]     data     セル中心data
   * @param [in]     id       セル中心data x方向サイズ
   * @param [in]     jd       セル中心data y方向サイズ
   * @param [in]     kd       セル中心data z方向サイズ
   * @param [in]     ivar     ベクトル成分 =0:x =1:y =2:z
   */
  void setVectorComponentGridData(int* size, int guide, double* d, double* data, int id, int jd, int kd, int ivar);


  /**
   * @brief Scalarの格子点での値をセット（ガイドセルの値は計算対象にいれていない）（double ---> float）
   * @param [out]    d        格子点data
   * @param [in]     data     セル中心data
   * @param [in]     id       セル中心data x方向サイズ
   * @param [in]     jd       セル中心data y方向サイズ
   * @param [in]     kd       セル中心data z方向サイズ
   */
  void setScalarGridData(int* size, int guide, float* d, double* data, int id, int jd, int kd);


  /**
   * @brief Vectorの格子点での値をセット（ガイドセルの値は計算対象にいれていない）（double ---> float）
   * @param [out]    d        格子点data
   * @param [in]     data     セル中心data
   * @param [in]     id       セル中心data x方向サイズ
   * @param [in]     jd       セル中心data y方向サイズ
   * @param [in]     kd       セル中心data z方向サイズ
   */
  void setVectorGridData(int* size, int guide, float* d, double* data, int id, int jd, int kd);


  /**
   * @brief 成分別Vectorの格子点での値をセット（ガイドセルの値は計算対象にいれていない）（double ---> float）
   * @param [out]    d        格子点data
   * @param [in]     data     セル中心data
   * @param [in]     id       セル中心data x方向サイズ
   * @param [in]     jd       セル中心data y方向サイズ
   * @param [in]     kd       セル中心data z方向サイズ
   * @param [in]     ivar     ベクトル成分 =0:x =1:y =2:z
   */
  void setVectorComponentGridData(int* size, int guide, float* d, double* data, int id, int jd, int kd, int ivar);


  /**
   * @brief Scalarの格子点での値をセット（ガイドセルの値は計算対象にいれていない）（float ---> double）
   * @param [out]    d        格子点data
   * @param [in]     data     セル中心data
   * @param [in]     id       セル中心data x方向サイズ
   * @param [in]     jd       セル中心data y方向サイズ
   * @param [in]     kd       セル中心data z方向サイズ
   */
  void setScalarGridData(int* size, int guide, double* d, float* data, int id, int jd, int kd);


  /**
   * @brief Vectorの格子点での値をセット（ガイドセルの値は計算対象にいれていない）（float ---> double）
   * @param [out]    d        格子点data
   * @param [in]     data     セル中心data
   * @param [in]     id       セル中心data x方向サイズ
   * @param [in]     jd       セル中心data y方向サイズ
   * @param [in]     kd       セル中心data z方向サイズ
   */
  void setVectorGridData(int* size, int guide, double* d, float* data, int id, int jd, int kd);


  /**
   * @brief 成分別Vectorの格子点での値をセット（ガイドセルの値は計算対象にいれていない）（float ---> double）
   * @param [out]    d        格子点data
   * @param [in]     data     セル中心data
   * @param [in]     id       セル中心data x方向サイズ
   * @param [in]     jd       セル中心data y方向サイズ
   * @param [in]     kd       セル中心data z方向サイズ
   * @param [in]     ivar     ベクトル成分 =0:x =1:y =2:z
   */
  void setVectorComponentGridData(int* size, int guide, double* d, float* data, int id, int jd, int kd, int ivar);


  /**
   * @brief 内部の格子点のデータを8で割る（float）
   * @param [out]    d        格子点data
   * @param [in]     id       セル中心d x方向サイズ
   * @param [in]     jd       セル中心d y方向サイズ
   * @param [in]     kd       セル中心d z方向サイズ
   */
  void VolumeDataDivideBy8(float* d, int id, int jd, int kd);


  /**
   * @brief 面上の格子点のデータを4で割る（float）
   * @param [out]    d        格子点data
   * @param [in]     id       セル中心d x方向サイズ
   * @param [in]     jd       セル中心d y方向サイズ
   * @param [in]     kd       セル中心d z方向サイズ
   */
  void FaceDataDivideBy4(float* d, int id, int jd, int kd);


  /**
   * @brief 辺上の格子点のデータを2で割る（float）
   * @param [out]    d        格子点data
   * @param [in]     id       セル中心d x方向サイズ
   * @param [in]     jd       セル中心d y方向サイズ
   * @param [in]     kd       セル中心d z方向サイズ
   */
  void LineDataDivideBy2(float* d, int id, int jd, int kd);


  /**
   * @brief 内部の格子点のデータを8で割る（double）
   * @param [out]    d        格子点data
   * @param [in]     id       セル中心d x方向サイズ
   * @param [in]     jd       セル中心d y方向サイズ
   * @param [in]     kd       セル中心d z方向サイズ
   */
  void VolumeDataDivideBy8(double* d, int id, int jd, int kd);


  /**
   * @brief 面上の格子点のデータを4で割る（double）
   * @param [out]    d        格子点data
   * @param [in]     id       セル中心d x方向サイズ
   * @param [in]     jd       セル中心d y方向サイズ
   * @param [in]     kd       セル中心d z方向サイズ
   */
  void FaceDataDivideBy4(double* d, int id, int jd, int kd);


  /**
   * @brief 辺上の格子点のデータを2で割る（double）
   * @param [out]    d        格子点data
   * @param [in]     id       セル中心d x方向サイズ
   * @param [in]     jd       セル中心d y方向サイズ
   * @param [in]     kd       セル中心d z方向サイズ
   */
  void LineDataDivideBy2(double* d, int id, int jd, int kd);


  /**
   * @brief PLOT3Dファイルの出力（セル中心での出力）
   */
  void output_plot3d_on_cell();

  /**
   * @brief xyzファイルの出力（sph:float,plot3d:float）
   * @param [in] m_step     ステップ 
   * @param [in] m_rank     ランク
   * @param [in] guide      ガイドセル数
   * @param [in] origin     基点座標
   * @param [in] pitch      ピッチ
   * @param [in] size       セルサイズ
   * @param [in] x          x方向座標ワーク
   * @param [in] y          y方向座標ワーク
   * @param [in] z          z方向座標ワーク
   */
  void OutputPlot3D_xyz_on_cell(int m_step, int m_rank, int guide, float* origin, float* pitch, int* size, float* x, float* y, float* z);

};

#endif // _COMB_H_
