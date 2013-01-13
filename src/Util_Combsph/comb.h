#ifndef _COMB_H_
#define _COMB_H_

// #################################################################
//
// Combine sph files and output
//
// Copyright (c) 2012-2013  All right reserved.
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
  
  //並列実行時のSTAGINGのON/OFF
  unsigned staging;
  
public:
  
  string filename;
  string out_dirname;
  string in_dirname;
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
  
  ::DFI DFI;                 ///< 分散ファイルインデクス管理クラス
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
    procGrp = 0;
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
   * @param [in] 間引きオプションのためのメモリ使用量最大値
   * @param [in] fp     ファイルポインタ
   */
  void MemoryRequirement(const double TotalMemory, const double sphMemory, const double plot3dMemory, const double thinMemory, FILE* fp);
  
  
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
   * @brief sphファイルのデータの連結
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
  template<class T1, class T2>
  void CombineLayerData(
                        T1* d, int size_d, T2* s, int size_s,
                        int xsize, int ysize, int zsize,
                        int z, int* sz, int dim, int xs, int ys) //z=zzz ; z!=zh
  {
    int ix = sz[0];
    int jy = sz[1];
    int kz = sz[2];
    
    for(int j=0;j<jy;j++){
      int yp=ys+j;
      int yrest=yp%thin_count;
      if(yrest != 0) continue;
      
      int ipss=z*dim*ix*jy+j*dim*ix;
      int ipds=(yp/thin_count)*dim*xsize+(xs/thin_count)*dim;
      
      int ic=0;
      for(int i=0;i<ix;i++){
        int xp=xs+i;
        int xrest=xp%thin_count;
        if(xrest != 0) continue;
        
        for(int idim=0;idim<dim;idim++){
          d[ipds+ic]=(T1)s[ipss+i*dim+idim];
          ic++;
        }
      }
    }
    
    return;
  };
  
  
  ///////////////////////////////////////////////////////////////////////////////
  // comb_plot3d.C
  
  /**
   * @brief PLOT3Dファイルの出力
   */
  void output_plot3d();
  
  
  /**
   * @brief 辺上の格子点のデータを2で割る（float）
   * @param [out]    dt       間引き後格子点data
   * @param [in]     d        格子点data
   * @param [in]     idt      間引き後セル中心dt x方向サイズ
   * @param [in]     jdt      間引き後セル中心dt y方向サイズ
   * @param [in]     kdt      間引き後セル中心dt z方向サイズ
   * @param [in]     id       セル中心d x方向サイズ
   * @param [in]     jd       セル中心d y方向サイズ
   * @param [in]     kd       セル中心d z方向サイズ
   */
  void thinout_plot3d(float* dt, float* d, int idt, int jdt, int kdt, int id, int jd, int kd);
  
  
  /**
   * @brief 辺上の格子点のデータを2で割る（double）
   * @param [out]    dt       間引き後格子点data
   * @param [in]     d        格子点data
   * @param [in]     idt      間引き後セル中心dt x方向サイズ
   * @param [in]     jdt      間引き後セル中心dt y方向サイズ
   * @param [in]     kdt      間引き後セル中心dt z方向サイズ
   * @param [in]     id       セル中心d x方向サイズ
   * @param [in]     jd       セル中心d y方向サイズ
   * @param [in]     kd       セル中心d z方向サイズ
   */
  void thinout_plot3d(double* dt, double* d, int idt, int jdt, int kdt, int id, int jd, int kd);
  
  
  /**
   * @brief xyzファイルの出力
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
  template<class T1, class T2>
  void OutputPlot3D_xyz(int m_step, int m_rank, int guide, T1* origin, T1* pitch, int* size, T2* x, T2* y, T2* z)
  {
    //value
    int ngrid=1;
    int ix = size[0];
    int jx = size[1];
    int kx = size[2];
    int gd = guide;
    int gc_out = 0;//plot3dは常にガイドセルを出力しない
    
    int *iblank=NULL;//dummy
    int id,jd,kd;//出力サイズ
    id=size[0]+1;//+2*gc_out
    jd=size[1]+1;//+2*gc_out
    kd=size[2]+1;//+2*gc_out
    
    //間引きのための処理
    int irest=(id-1)%thin_count;
    int jrest=(jd-1)%thin_count;
    int krest=(kd-1)%thin_count;
    id=(id-1)/thin_count;
    jd=(jd-1)/thin_count;
    kd=(kd-1)/thin_count;
    id=id+1;
    jd=jd+1;
    kd=kd+1;
    if(irest!=0) id=id+1;
    if(jrest!=0) jd=jd+1;
    if(krest!=0) kd=kd+1;
    
    // ガイドセル出力があった場合オリジナルポイントを調整しておく
    T2 m_org[3], m_pit[3];
    for (int i=0; i<3; i++)
    {
      m_org[i] = (T2)origin[i] + (T2)pitch[i]*(T2)gd;
      m_pit[i] = (T2)pitch[i];
    }
    
    // 出力ファイル名
    std::string tmp;
    tmp = Generate_FileName_Free(P3Op.basename, "xyz", m_step, m_rank, true);
    tmp = out_dirname + tmp;
    
    //open file
    FP3DW.setFileName(tmp.c_str());
    if(!FP3DW.OpenFile()){
      printf("Error : error OpenFile\n");
      Exit(0);
    }
    
    //write block data
    FP3DW.WriteNgrid(ngrid);//if multi grid
    FP3DW.WriteBlockData(id,jd,kd);
    
    for(int k=0;k<kd;k++){
      for(int j=0;j<jd;j++){
        for(int i=0;i<id;i++){
          size_t ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
          x[ip]=m_org[0]+(T2)thin_count*m_pit[0]*(T2)i;//-pitch[0]*(float)gc_out;
          y[ip]=m_org[1]+(T2)thin_count*m_pit[1]*(T2)j;//-pitch[1]*(float)gc_out;
          z[ip]=m_org[2]+(T2)thin_count*m_pit[2]*(T2)k;//-pitch[2]*(float)gc_out;
        }
      }
    }
    
    //x direction modify
    if(irest!=0 && (id-2)>=0 ){
      for(int k=0;k<kd;k++){
        for(int j=0;j<jd;j++){
          //for(int i=id-1;i<id;i++){
          size_t ip = _F_IDX_S3D(id, j+1, k+1, id, jd, kd, 0);
          //size_t ip=k*id*jd+j*id+id-1;
          x[ip]=m_org[0]+(T2)thin_count*m_pit[0]*(T2)(id-2)+(T2)irest*m_pit[0];//-pitch[0]*(float)gc_out;
          //}
        }
      }
    }
    
    //y direction modify
    if(jrest!=0 && (jd-2)>=0 ){
      for(int k=0;k<kd;k++){
        //for(int j=jd-1;j<jd;j++){
        for(int i=0;i<id;i++){
          size_t ip = _F_IDX_S3D(i+1, jd, k+1, id, jd, kd, 0);
          //size_t ip=k*id*jd+(jd-1)*id+i;
          y[ip]=m_org[1]+(T2)thin_count*m_pit[1]*(T2)(jd-2)+(T2)jrest*m_pit[1];//-pitch[1]*(float)gc_out;
        }
        //}
      }
    }
    
    //z direction modify
    if(krest!=0 && (kd-2)>=0 ){
      //for(int k=kd-1;k<kd;k++){
      for(int j=0;j<jd;j++){
        for(int i=0;i<id;i++){
          size_t ip = _F_IDX_S3D(i+1, j+1, kd, id, jd, kd, 0);
          //size_t ip=(kd-1)*id*jd+j*id+i;
          z[ip]=m_org[2]+(T2)thin_count*m_pit[2]*(T2)(kd-2)+(T2)krest*m_pit[2];//-pitch[2]*(float)gc_out;
        }
      }
      //}
    }
    
    //z direction modify
    if(krest!=0){
      for(int k=kd-1;k<kd;k++){
        for(int j=0;j<jd;j++){
          for(int i=0;i<id;i++){
            size_t ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
            z[ip]=m_org[2]+(T2)krest*m_pit[2]*(T2)k;//-pitch[2]*(float)gc_out;
          }
        }
      }
    }
    
    //write
    FP3DW.setGridData(id,jd,kd,ngrid);
    FP3DW.setXYZData(x,y,z,iblank);
    if(!FP3DW.WriteXYZData()) printf("\terror WriteXYZData\n");
    
    //close file
    FP3DW.CloseFile();
    
  };
  
  /**
   * @brief 内部の格子点のデータを8で割る
   * @param [out]    d        格子点data
   * @param [in]     id       セル中心d x方向サイズ
   * @param [in]     jd       セル中心d y方向サイズ
   * @param [in]     kd       セル中心d z方向サイズ
   */
  template<class T>
  void VolumeDataDivideBy8(T* d, int id, int jd, int kd)
  {
    int i,j,k;
    size_t ip;
    
    for (k=1; k<kd-1; k++){
      for (j=1; j<jd-1; j++){
        for (i=1; i<id-1; i++){
          ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
          d[ip]=d[ip]*0.125;
        }
      }
    }
  };
  
  /**
   * @brief 面上の格子点のデータを4で割る
   * @param [out]    d        格子点data
   * @param [in]     id       セル中心d x方向サイズ
   * @param [in]     jd       セル中心d y方向サイズ
   * @param [in]     kd       セル中心d z方向サイズ
   */
  template<class T>
  void FaceDataDivideBy4(T* d, int id, int jd, int kd)
  {
    int i,j,k;
    size_t ip;
    
    i=0;
    for (k=1; k<kd-1; k++){
      for (j=1; j<jd-1; j++){
        ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
        d[ip]=d[ip]*0.25;
      }
    }
    
    i=id-1;
    for (k=1; k<kd-1; k++){
      for (j=1; j<jd-1; j++){
        ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
        d[ip]=d[ip]*0.25;
      }
    }
    
    j=0;
    for (k=1; k<kd-1; k++){
      for (i=1; i<id-1; i++){
        ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
        d[ip]=d[ip]*0.25;
      }
    }
    
    j=jd-1;
    for (k=1; k<kd-1; k++){
      for (i=1; i<id-1; i++){
        ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
        d[ip]=d[ip]*0.25;
      }
    }
    
    k=0;
    for (j=1; j<jd-1; j++){
      for (i=1; i<id-1; i++){
        ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
        d[ip]=d[ip]*0.25;
      }
    }
    
    k=kd-1;
    for (j=1; j<jd-1; j++){
      for (i=1; i<id-1; i++){
        ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
        d[ip]=d[ip]*0.25;
      }
    }
  };
  
  
  /**
   * @brief 辺上の格子点のデータを2で割る
   * @param [out]    d        格子点data
   * @param [in]     id       セル中心d x方向サイズ
   * @param [in]     jd       セル中心d y方向サイズ
   * @param [in]     kd       セル中心d z方向サイズ
   */
  template<class T>
  void LineDataDivideBy2(T* d, int id, int jd, int kd)
  {
    int i,j,k;
    size_t ip;
    
    i=0; j=0;
    for (k=1; k<kd-1; k++){
      ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
      d[ip]=d[ip]*0.5;
    }
    
    i=0; j=jd-1;
    for (k=1; k<kd-1; k++){
      ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
      d[ip]=d[ip]*0.5;
    }
    
    i=0; k=0;
    for (j=1; j<jd-1; j++){
      ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
      d[ip]=d[ip]*0.5;
    }
    
    i=0; k=kd-1;
    for (j=1; j<jd-1; j++){
      ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
      d[ip]=d[ip]*0.5;
    }
    
    j=0; k=0;
    for (i=1; i<id-1; i++){
      ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
      d[ip]=d[ip]*0.5;
    }
    
    j=0; k=kd-1;
    for (i=1; i<id-1; i++){
      ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
      d[ip]=d[ip]*0.5;
    }
    
    j=jd-1; k=0;
    for (i=1; i<id-1; i++){
      ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
      d[ip]=d[ip]*0.5;
    }
    
    j=jd-1; k=kd-1;
    for (i=1; i<id-1; i++){
      ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
      d[ip]=d[ip]*0.5;
    }
    
    i=id-1; j=0;
    for (k=1; k<kd-1; k++){
      ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
      d[ip]=d[ip]*0.5;
    }
    
    i=id-1; j=jd-1;
    for (k=1; k<kd-1; k++){
      ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
      d[ip]=d[ip]*0.5;
    }
    
    i=id-1; k=0;
    for (j=1; j<jd-1; j++){
      ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
      d[ip]=d[ip]*0.5;
    }
    
    i=id-1; k=kd-1;
    for (j=1; j<jd-1; j++){
      ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
      d[ip]=d[ip]*0.5;
    }
    
  };
  
  
  /**
   * @brief Scalarの格子点での値をセット
   * @param [out]    d        格子点data
   * @param [in]     data     セル中心data
   * @param [in]     id       セル中心data x方向サイズ
   * @param [in]     jd       セル中心data y方向サイズ
   * @param [in]     kd       セル中心data z方向サイズ
   */
  template<class T1, class T2>
  void setScalarGridData(int* size, int guide, T1* d, T2* data, int id, int jd, int kd, int gc_out)
  {
    int ix = size[0];
    int jx = size[1];
    int kx = size[2];
    int gd = guide;
    
    size_t mip;
    float ddd;
    int i,j,k;
    
    size_t dsize = (size_t)(id*jd*kd);
    
    for (size_t l=0; l<dsize; l++) d[l]=0.0;
    
    for (int km=1-gc_out; km<=kx+gc_out; km++) {
      for (int jm=1-gc_out; jm<=jx+gc_out; jm++) {
        for (int im=1-gc_out; im<=ix+gc_out; im++) {
          mip = _F_IDX_S3D(im, jm, km, ix, jx, kx, gd);
          ddd=(T1)data[mip];
          i=im-1+gc_out;
          j=jm-1+gc_out;
          k=km-1+gc_out;
          size_t ip1 = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
          size_t ip2 = _F_IDX_S3D(i+2, j+1, k+1, id, jd, kd, 0);
          size_t ip3 = _F_IDX_S3D(i+2, j+2, k+1, id, jd, kd, 0);
          size_t ip4 = _F_IDX_S3D(i+1, j+2, k+1, id, jd, kd, 0);
          size_t ip5 = _F_IDX_S3D(i+1, j+1, k+2, id, jd, kd, 0);
          size_t ip6 = _F_IDX_S3D(i+2, j+1, k+2, id, jd, kd, 0);
          size_t ip7 = _F_IDX_S3D(i+2, j+2, k+2, id, jd, kd, 0);
          size_t ip8 = _F_IDX_S3D(i+1, j+2, k+2, id, jd, kd, 0);
          d[ip1]=d[ip1]+ddd;
          d[ip2]=d[ip2]+ddd;
          d[ip3]=d[ip3]+ddd;
          d[ip4]=d[ip4]+ddd;
          d[ip5]=d[ip5]+ddd;
          d[ip6]=d[ip6]+ddd;
          d[ip7]=d[ip7]+ddd;
          d[ip8]=d[ip8]+ddd;
        }
      }
    }
    
    //内部の格子点のデータを8で割る
    VolumeDataDivideBy8(d, id, jd, kd);
    
    //面上の格子点のデータを4で割る
    FaceDataDivideBy4(d, id, jd, kd);
    
    //辺上の格子点のデータを2で割る
    LineDataDivideBy2(d, id, jd, kd);
    
  };
  
  
  /**
   * @brief Vectorの格子点での値をセット
   * @param [out]    d        格子点data
   * @param [in]     data     セル中心data
   * @param [in]     id       セル中心data x方向サイズ
   * @param [in]     jd       セル中心data y方向サイズ
   * @param [in]     kd       セル中心data z方向サイズ
   */
  template<class T1, class T2>
  void setVectorGridData(int* size, int guide, T1* d, T2* data, int id, int jd, int kd, int gc_out)
  {
    int ix = size[0];
    int jx = size[1];
    int kx = size[2];
    int gd = guide;
    
    size_t mip;
    T1 ddd;
    int i,j,k;
    
    size_t dsize = (size_t)(id*jd*kd);
    size_t dsize3 = (size_t)(id*jd*kd*3);
    
    for (size_t l=0; l<dsize3; i++) d[l]=0.0;
    
    for (size_t ivar=0;ivar<3;ivar++){
      
      for (int km=1-gc_out; km<=kx+gc_out; km++) {
        for (int jm=1-gc_out; jm<=jx+gc_out; jm++) {
          for (int im=1-gc_out; im<=ix+gc_out; im++) {
            mip = _F_IDX_V3D(im, jm, km, ivar, ix, jx, kx, gd);
            ddd=(T1)data[mip];
            i=im-1+gc_out;
            j=jm-1+gc_out;
            k=km-1+gc_out;
            size_t ip1 = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
            size_t ip2 = _F_IDX_S3D(i+2, j+1, k+1, id, jd, kd, 0);
            size_t ip3 = _F_IDX_S3D(i+2, j+2, k+1, id, jd, kd, 0);
            size_t ip4 = _F_IDX_S3D(i+1, j+2, k+1, id, jd, kd, 0);
            size_t ip5 = _F_IDX_S3D(i+1, j+1, k+2, id, jd, kd, 0);
            size_t ip6 = _F_IDX_S3D(i+2, j+1, k+2, id, jd, kd, 0);
            size_t ip7 = _F_IDX_S3D(i+2, j+2, k+2, id, jd, kd, 0);
            size_t ip8 = _F_IDX_S3D(i+1, j+2, k+2, id, jd, kd, 0);
            d[ip1+dsize*ivar]=d[ip1+dsize*ivar]+ddd;
            d[ip2+dsize*ivar]=d[ip2+dsize*ivar]+ddd;
            d[ip3+dsize*ivar]=d[ip3+dsize*ivar]+ddd;
            d[ip4+dsize*ivar]=d[ip4+dsize*ivar]+ddd;
            d[ip5+dsize*ivar]=d[ip5+dsize*ivar]+ddd;
            d[ip6+dsize*ivar]=d[ip6+dsize*ivar]+ddd;
            d[ip7+dsize*ivar]=d[ip7+dsize*ivar]+ddd;
            d[ip8+dsize*ivar]=d[ip8+dsize*ivar]+ddd;
          }
        }
      }
      
      //内部の格子点のデータを8で割る
      VolumeDataDivideBy8(&d[dsize*ivar], id, jd, kd);
      
      //面上の格子点のデータを4で割る
      FaceDataDivideBy4(&d[dsize*ivar], id, jd, kd);
      
      //辺上の格子点のデータを2で割る
      LineDataDivideBy2(&d[dsize*ivar], id, jd, kd);
      
      //境界条件処理
      
      
    }//loop ivar
    
  };
  
  
  /**
   * @brief 成分別Vectorの格子点での値をセット
   * @param [out]    d        格子点data
   * @param [in]     data     セル中心data
   * @param [in]     id       セル中心data x方向サイズ
   * @param [in]     jd       セル中心data y方向サイズ
   * @param [in]     kd       セル中心data z方向サイズ
   * @param [in]     ivar     ベクトル成分 =0:x =1:y =2:z
   */
  template<class T1, class T2>
  void setVectorComponentGridData(int* size, int guide, T1* d, T2* data, int id, int jd, int kd, int ivar, int gc_out)
  {
    int ix = size[0];
    int jx = size[1];
    int kx = size[2];
    int gd = guide;
    
    size_t mip;
    T1 ddd;
    int i,j,k;
    
    size_t dsize = (size_t)(id*jd*kd);
    
    for (size_t l=0; l<dsize; l++) d[l]=0.0;
    
    for (int km=1-gc_out; km<=kx+gc_out; km++) {
      for (int jm=1-gc_out; jm<=jx+gc_out; jm++) {
        for (int im=1-gc_out; im<=ix+gc_out; im++) {
          mip = _F_IDX_V3D(im, jm, km, ivar, ix, jx, kx, gd);
          ddd=(T1)data[mip];
          i=im-1+gc_out;
          j=jm-1+gc_out;
          k=km-1+gc_out;
          size_t ip1 = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
          size_t ip2 = _F_IDX_S3D(i+2, j+1, k+1, id, jd, kd, 0);
          size_t ip3 = _F_IDX_S3D(i+2, j+2, k+1, id, jd, kd, 0);
          size_t ip4 = _F_IDX_S3D(i+1, j+2, k+1, id, jd, kd, 0);
          size_t ip5 = _F_IDX_S3D(i+1, j+1, k+2, id, jd, kd, 0);
          size_t ip6 = _F_IDX_S3D(i+2, j+1, k+2, id, jd, kd, 0);
          size_t ip7 = _F_IDX_S3D(i+2, j+2, k+2, id, jd, kd, 0);
          size_t ip8 = _F_IDX_S3D(i+1, j+2, k+2, id, jd, kd, 0);
          d[ip1]=d[ip1]+ddd;
          d[ip2]=d[ip2]+ddd;
          d[ip3]=d[ip3]+ddd;
          d[ip4]=d[ip4]+ddd;
          d[ip5]=d[ip5]+ddd;
          d[ip6]=d[ip6]+ddd;
          d[ip7]=d[ip7]+ddd;
          d[ip8]=d[ip8]+ddd;
        }
      }
    }
    
    //内部の格子点のデータを8で割る
    VolumeDataDivideBy8(d, id, jd, kd);
    
    //面上の格子点のデータを4で割る
    FaceDataDivideBy4(d, id, jd, kd);
    
    //辺上の格子点のデータを2で割る
    LineDataDivideBy2(d, id, jd, kd);
    
  };
  
};

#endif // _COMB_H_
