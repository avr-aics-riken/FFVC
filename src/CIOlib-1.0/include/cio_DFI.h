#ifndef _CIO_DFI_H_
#define _CIO_DFI_H_

/* #################################################################
 *
 * CIOlib - Cartesian Input / Output library
 *
 * Copyright (c) AICS, RIKEN. All right reserved. 2013
 *
 * #################################################################
 */

/** 
 * @file   cio_DFI.h
 * @brief  cio_DFI Class Header
 * @author kero    
 */
#include <stdlib.h>
#include <errno.h>
#include <sys/stat.h>

#include "cio_TextParser.h"
#include "cio_Define.h"
#include "cio_DFI_Ffunc.h"

using namespace std;

/** CIO main class */
class cio_DFI {


public:

  MPI_Comm m_comm;
  string   m_indexDfiName;
  int m_dfi_mng;

  /** index.dfi ファイルの FileInfo */
  struct cio_FileInfo
  {
    string DirectoryPath;                 ///<ディレクトリパス
    string Prefix;                        ///<ファイル接頭文字
    string FileFormat;                    ///<ファイルフォーマット "bov","sph",,,
    int    GuideCell;                     ///<仮想セルの数
    string DataType;                      ///<配列のデータタイプ "S3D","V3D"
    string Endian;                        ///<エンディアンタイプ "big","little"
    string ArrayShape;                    ///<配列形状
    int    Component;                     ///<成分数

    cio_FileInfo() 
    {
      DirectoryPath="";
      Prefix       ="";
      FileFormat   ="";
      GuideCell    =0;
      DataType     ="";
      Endian       ="";
      ArrayShape   ="";
      Component    =0;
    }

    cio_FileInfo(string _DirectoryPath, string _Prefix, string _FileFormat,
             int _GuideCell, string _DataType, string _Endian, 
             string _ArrayShape, int _Component)
    {
      DirectoryPath=_DirectoryPath;
      Prefix       =_Prefix;
      FileFormat   =_FileFormat;
      GuideCell    =_GuideCell;
      DataType     =_DataType;
      Endian       =_Endian;
      ArrayShape   =_ArrayShape;
      Component    =_Component;
    }   

  };


  /** index.dfi ファイルの FilePath */
  struct cio_FilePath
  {
    string Process;                       ///<proc.dfi ファイル名

    cio_FilePath()
    {
      Process="";
    }

    cio_FilePath(string _Process)
    {
      Process=_Process;
    }

  };


  /** index.dfi ファイルの Unit */
  struct cio_Unit  
  {
    bool out_Length;                      ///<Length,L0 出力フラグ
    string Length;                        ///<(NonDimensional, m, cm, mm)
    REAL_TYPE L0;                         ///<規格化に用いた長さスケール

    bool out_Velocity;                    ///<Velocity,V0 出力フラグ
    string Velocity;                      ///<(NonDimensional, m/s)     
    REAL_TYPE V0;                         ///<代表速度(m/s)

    bool out_Pressure;                    ///<Presuure,P0,DiffPrs出力フラグ
    string Pressure;                      ///<(NonDimensional, Pa)      
    REAL_TYPE P0;                         ///<基準圧力(Pa)
    REAL_TYPE DiffPrs;                    ///<圧力差(Pa)

    bool out_Temperature;                 ///<Temperatur,BaseTemp,DiffTemp出力フラグ
    string Temperatur;                    ///<(NonDimensional, C, K)    
    REAL_TYPE BaseTemp;                   ///<指定単力　　　　
    REAL_TYPE DiffTemp;                   ///<指定単位　　　

    cio_Unit()
    {
      out_Length=false;
      Length    ="";
      L0        =0.0;

      out_Velocity=false;
      Velocity    ="";
      V0          =0.0;

      out_Pressure=false;
      Pressure    ="";
      P0          =0.0;
      DiffPrs     =0.0;

      out_Temperature=false;
      Temperatur     ="";
      BaseTemp       =0.0;
      DiffTemp       =0.0;
    }

    cio_Unit(bool _out_Length, string _Length, REAL_TYPE _L0,
         bool _out_Velocity, string _Velocity, REAL_TYPE _V0,
         bool _out_Pressure, string _Pressure, REAL_TYPE _P0, REAL_TYPE _DiffPrs,
         bool _out_Temperature, string _Temperatur, REAL_TYPE _BaseTemp, REAL_TYPE _DiffTemp) 
    {

      out_Length=_out_Length;
      Length    =_Length;
      L0        =_L0;

      out_Velocity=_out_Velocity;
      Velocity    =_Velocity;
      V0          =_V0;

      out_Pressure=_out_Pressure;
      Pressure    =_Pressure;
      P0          =_P0;
      DiffPrs     =_DiffPrs;

      out_Temperature=_out_Temperature;
      Temperatur     =_Temperatur;
      BaseTemp       =_BaseTemp;
      DiffTemp       =_DiffTemp;
    }

  };

  /** index.dfi ファイルの Slice */
  struct cio_Slice
  {
    int    step;                          ///<ステップ番号
    REAL_TYPE time;                       ///<時刻
    vector<REAL_TYPE> Min;                ///<最小値
    vector<REAL_TYPE> Max;                ///<最大値
  };
  
  /** proc.dfi ファイルの Domain */
  struct cio_Domain
  {
    REAL_TYPE GlobalOrigin[3];             ///<計算空間の起点座標
    REAL_TYPE GlobalRegion[3];             ///<計算空間の各軸方向の長さ
    int GlobalVoxel[3];                    ///<計算領域全体のボクセル数
    int GlobalDivision[3];                 ///<計算領域の分割数
    string ActiveSubdomain;                ///<ActiveSubdomainファイル名

    cio_Domain()
    {
      for(int i=0; i<3; i++) GlobalOrigin[i]=0.0;
      for(int i=0; i<3; i++) GlobalRegion[i]=0.0;
      for(int i=0; i<3; i++) GlobalVoxel[i]=0;
      for(int i=0; i<3; i++) GlobalDivision[i]=0;
      ActiveSubdomain="";
    }

    cio_Domain(REAL_TYPE* _GlobalOrigin, REAL_TYPE* _GlobalRegion, int* _GlobalVoxel, 
           int* _GlobalDivision)
    {
      GlobalOrigin[0]=_GlobalOrigin[0];
      GlobalOrigin[1]=_GlobalOrigin[1];
      GlobalOrigin[2]=_GlobalOrigin[2];

      GlobalRegion[0]=_GlobalRegion[0];
      GlobalRegion[1]=_GlobalRegion[1];
      GlobalRegion[2]=_GlobalRegion[2];

      GlobalVoxel[0]=_GlobalVoxel[0];
      GlobalVoxel[1]=_GlobalVoxel[1];
      GlobalVoxel[2]=_GlobalVoxel[2];

      GlobalDivision[0]=_GlobalDivision[0];
      GlobalDivision[1]=_GlobalDivision[1];
      GlobalDivision[2]=_GlobalDivision[2];
    }

  };

  /** proc.dfi ファイルの MPI */
  struct cio_MPI
  {
    int NumberOfRank;                      ///<プロセス数
    int NumberOfGroup;                     ///<グループ数

    cio_MPI()
    {
       NumberOfRank=0;
       NumberOfGroup=1;
    }

    cio_MPI(int _NumberOfRank)
    {
       NumberOfRank=_NumberOfRank;
    }

  };

  /** proc.dfi ファイルの Process */
  struct cio_Rank
  {
    int RankID;                           ///<ランク番号
    string HostName;                      ///<ホスト名 
    int VoxelSize[3];                     ///<ボクセルサイズ
    int HeadIndex[3];                     ///<始点インデックス
    int TailIndex[3];                     ///<終点インデックス
  };

  cio_FileInfo DFI_Finfo;
  cio_FilePath DFI_Fpath;
  cio_Unit     DFI_Unit;
  cio_Domain   DFI_Domain;
  cio_MPI      DFI_MPI;
  vector<cio_Slice> TimeSlice;
  vector<cio_Rank> RankInfo;

public:
  /** コンストラクタ */
  cio_DFI();
  
  /**　デストラクタ */
  ~cio_DFI();

  /**
   * @brief read インスタンス
   * @param [in] comm    MPIコミュニケータ
   * @param [in] dfifile DFIファイル名
   * @return インスタンスされたクラスのポインタ
   */
  static cio_DFI* ReadInit(MPI_Comm comm, string dfifile); 

  /**
   * @brief write インスタンス
   * @param [in] comm          MPIコミュニケータ
   * @param [in] DfiName       DFIファイル名
   * @param [in] DirectoryPath フィールドデータのディレクトリ
   * @param [in] Prefix        ベースファイル名
   * @param [in] FileFormat    ファイルフォーマット
   * @param [in] GuideCell     出力仮想セル数　　　
   * @param [in] DataType      データタイプ　　　　
   * @param [in] ArrayShape    配列形状　　　　　　
   * @param [in] Component     成分数　　　　　　　
   * @param [in] Process       proc.dfiファイル名
   * @param [in] G_size[3]     グローバルボクセルサイズ　
   * @param [in] pitch[3]      ピッチ　　　　　　　　　　
   * @param [in] G_origin[3]   原点座標値　　　　　　　　
   * @param [in] division[3]   領域分割数　　　　　　　　
   * @param [in] head[3]       計算領域の開始位置　　　　
   * @param [in] tail[3]       計算領域の終了位置　　　　
   * @param [in] hostname      ホスト名　　　　　　　　　
   * @return インスタンスされたクラスのポインタ
   */
  static cio_DFI* WriteInit(MPI_Comm comm, 
                            string DfiName,  
                            string DirectoryPath,                
                            string Prefix,                
                            string FileFormat,           
                            int    GuideCell,           
                            string DataType,           
                            string ArrayShape,      
                            int    Component,
                            string Process,
                            int    G_size[3],
                            REAL_TYPE pitch[3],
                            REAL_TYPE G_origin[3],
                            int    division[3],
                            int    head[3],
                            int    tail[3],
                            string hostname);   

  /**
   * @brief field data  ファイル名の作成
   * @param [in] RankID ランク番号
   * @param [in] step   読込みステップ番号
   * @param [in] mio    並列判定フラグ（逐次or並列の判定用）
   * @return 生成されたファイル名　　　　　　　
   */
  std::string Generate_FileName(int RankID,int step, const bool mio);

  /**
   * @brief ディレクトリパスの作成
   * @param [in] path パス
   * @return error code　　　　　　　
   */ 
  int MakeDirectory(string path);

  /**
   * @brief initialise dfi
   */ 
  void InitDFI();

  /**
   * @brief read FileInfo(inde.dfi)
   * @param [in]   dfifile index.dfiファイル名
   * @param [in]   tpCntl  cio_TextParserクラス 
   * @param [out]  finfo   読込んだcio_FileInfo 
   * @return error code
   */
  static int readFileInfo(string dfifile, cio_TextParser tpCntl, cio_FileInfo &finfo);

  /**
   * @brief read FilePath(inde.dfi)
   * @param [in]   dfifile index.dfiファイル名
   * @param [in]   tpCntl  cio_TextParserクラス 
   * @param [out]  fpath   読込んだFilePath 
   * @return error code
   */
  static int readFilePath(string dfifile, cio_TextParser tpCntl, cio_FilePath &fpath);

  /**
   * @brief read Unit(inde.dfi)
   * @param [in]   dfifile index.dfiファイル名
   * @param [in]   tpCntl  cio_TextParserクラス 
   * @param [out]  unit    読込んだUnit 
   * @return error code
   */
  static int readUnit(string dfifile, cio_TextParser tpCntl, cio_Unit &unit);

  /**
   * @brief read Slice(inde.dfi)
   * @param [in]      dfifile   index.dfiファイル名
   * @param [in]      tpCntl    cio_TextParserクラス 
   * @param [out]     TimeSlice 読込んだSliceを格納した領域 
   * @param [out,out] slice     TimeSlice読込み用領域 
   * @return error code
   */
  static int readSlice(string dfifile, cio_TextParser tpCntl, vector<cio_Slice> &TimeSlice, cio_Slice  slice);

  /**
   * @brief read Domain(proc.dfi)
   * @param [in]   dfifile proc.dfiファイル名
   * @param [in]   tpCntl  cio_TextParserクラス 
   * @param [out]  domain  読込んだDomain
   * @return error code
   */
  static int readDomain(string dfifile, cio_TextParser tpCntl, cio_Domain &domain);

  /**
   * @brief read MPI(proc.dfi)
   * @param [in]   dfifile proc.dfiファイル名
   * @param [in]   tpCntl  cio_TextParserクラス 
   * @param [out]  mpi     読込んだMPI
   * @return error code
   */
  static int readMPI(string dfifile, cio_TextParser tpCntl, cio_MPI &mpi);

  /**
   * @brief read Rank(proc.dfi)
   * @param [in]   dfifile proc.dfiファイル名
   * @param [in]   tpCntl  cio_TextParserクラス 
   * @param [out]  rank    読込んだProcess
   * @return error code
   */
  static int readRank(string dfifile, cio_TextParser tpCntl, vector<cio_Rank> &RankInfo, cio_Rank rank);

  /**
   * @brief get ArrayShape （配列形状の取り出し関数）
   * @return 配列形状
   */
  std::string getArrayShape();

  /**
   * @brief get DataType （データタイプの取り出し関数）
   * @return データタイプ
   */
  std::string getDataType();

  /**
   * @brief get Component （成分数の取り出し関数）
   * @return 成分数
   */
  int getComponent();

  /**
   * @brief Create Domain & Process 
   * @param [in] comm          MPIコミュニケータ
   * @param [in] G_voxel[3]    グローバルボクセルサイズ　
   * @param [in] G_division[3] 領域分割数　　　　　　　　
   * @param [in] head[3]       計算領域の開始位置　　　　
   * @param [in] tail[3]       計算領域の終了位置　　　　
   * @param [out]G_domain      Domain情報(構造体)　　　　
   * @param [out]G_RankInfo    Process情報(vector)　　　
   * @param [in] G_Rank        Process情報(構造体)　　　
   */
  static void cio_Create_Domain(MPI_Comm comm,
                                int G_voxel[3], int G_division[3], int head[3], int tail[3],
                                cio_Domain &G_domain, vector<cio_Rank> &G_RankInfo, 
                                cio_Rank G_Rank);

  /**
   * @brief read field data record
   * @param [in] step 入力ステップ番号
   * @param [in] gc   仮想セル数　　　
   * @param [in] Gvoxel[3]    グローバルボクセルサイズ　
   * @param [in] Gdivision[3] 領域分割数　　　　　　　　
   * @param [in] head[3]       計算領域の開始位置　　　　
   * @param [in] tail[3]       計算領域の終了位置　　　　
   * @param [out]val           読み込んだデータポインタ　
   *
   */
  virtual void ReadData(int step, int gc, 
                        int Gvoxel[3], int Gdivision[3], int head[3], int tail[3],
                        REAL_TYPE *val, REAL_TYPE &time ){};

  /**
   * @brief 粗密データ判定
   * @param [in] Gvoxel     計算空間全体のボクセルサイズ（自）
   * @param [in] DFI_Gvoxel 計算空間全体のボクセルサイズ（DFI）
   * @return CIO_E_GV_SAME:密 CIO_E_GVX2_SAME:粗 CIO_E_OTHER:その他
   */
  cio_EGlobalVoxel CheckGlobalVoxel(int Gvoxel[3], int DFI_Gvoxel[3]); 

  /**
   * @brief 粗密ファイルのMxN判定
   * @param [in] rankList 読込みに必要なランク番号リスト
   * @param [in] head[3]  計算領域の開始位置　　　　
   * @param [in] tail[3]  計算領域の終了位置　　　　
   * @param [in] gc       仮想セル数(自)　　　
   * @param [in] dfi_gc   仮想セル数(DFI)　　　
   * @return  true:密 false:粗
   */
  bool CheckMxN(vector<int> &rankList, int head[3], int tail[3], int gc, int dfi_gc);
 
  /**
   * @brief 読込みランクファイルリストの作成
   * @param [in] head[3]  計算領域の開始位置　　　　
   * @param [in] tail[3]  計算領域の終了位置　　　　
   * @param [in] gc       仮想セル数(自)　　　
   * @param [in] readflag 粗密データ判定フラグ
   * @param [out]rankList 読込みに必要なランク番号リスト 
   */
  void CreateRankList(int head[3], int tail[3], int gc, cio_EGlobalVoxel readflag,
                      vector<int> &rankList);

  /**
   * @brief 読込み範囲を求める
   * @param [in] head[3]     計算領域の開始位置(自)　
   * @param [in] tail[3]     計算領域の終了位置(自)　
   * @param [in] gc          仮想セル数(自)　
   * @param [in] DEF_head[3] 計算領域の開始位置(DFI)　　
   * @param [in] DEF_tail[3] 計算領域の終了位置(DFI)　　
   * @param [in] DFI_gc      仮想セル数(DFI)　
   * @param [out]sta[3]      読込み開始位置
   * @param [out]end[3]      読込み終了位置　　
   * @return true:１対１
   */
  bool CheckReadArea(int head[3], int tail[3], int gc, int DEF_head[3], int DFI_tail[3],
                     int DFI_gc, cio_EGlobalVoxel readflag, int sta[3], int end[3]);


  /**
   * @brief write field data record
   * @param [in] step 出力ステップ番号
   * @param [in] gc   仮想セル数　　　
   * @param [in] tiem 出力時刻　　　　
   * @param [in] val  出力データポインタ
   */ 
  virtual void WriteData(int step, int gc, REAL_TYPE time, 
                        REAL_TYPE *val, REAL_TYPE *minmax){};

  /**
   * @brief index DFIファイル出力コントロール
   * @param [in] prefix  ファイル接頭文字
   * @param [in] step    ステップ
   * @param [in] time    時間　　
   * @param [in] mio     出力時の分割指定　 true = local / false = gather
   * @return true:出力成功 false:出力失敗
   */
  bool WriteIndexDfiFile(string DfiName, int RabkID,const std::string prefix, const unsigned step, REAL_TYPE time, REAL_TYPE *minmax, const bool mio); 

  /**
   * @brief proc DFIファイル出力コントロール
   * @param [in] RankID ランク番号　　　　　　　　　
   * @return true:出力成功 false:出力失敗
   */
  bool WriteProcDfiFile(int RankID); 

  /**
   * @brief proc DFIファイル出力コントロール
   * @param [in] comm          MPIコミュニケータ
   * @param [in] procFileName  出力proc.dfiファイル名
   * @param [in] G_size[3]     グローバルボクセルサイズ　
   * @param [in] division[3]   領域分割数　　　　　　　　
   * @param [in] head[3]       計算領域の開始位置　　　　
   * @param [in] tail[3]       計算領域の終了位置　　　　
   * @param [in] hostname      ホスト名　　　　　　　　　
   * @param [in] out_host      ホスト名出力フラグ　　　　
   * @return true:出力成功 false:出力失敗
   */
  static bool WriteProcDfiFile(MPI_Comm comm, string procFileName, int G_size[3],
                               int division[3], int head[3], int tail[3], 
                               REAL_TYPE org[3], REAL_TYPE pch[3], string hostname,
                               bool out_host); 

  /**
   * @brief 出力DFIファイル名を作成する
   * @param [in] prefix ファイル接頭文字
   * @return DFIファイル名
   */ 
  std::string Generate_DFI_Name(const std::string prefix);

  /**
   * @brief DFIファイルを出力する
   * @param [in] dfi_name  DFIファイル名
   * @param [in] prefix    ファイル接頭文字
   * @param [in] step      ステップ数
   * @param [in] time      時間　　
   * @param [in] dfi_mng   出力管理カウンタ
   * @param [in] mio       出力時の分割指定　 true = local / false = gather
   * @return true:出力成功 false:出力失敗
   */
  bool Write_Index_File(const std::string dfi_name, const std::string prefix, const unsigned step, REAL_TYPE time, int& dfi_mng, REAL_TYPE *minmax, const bool mio); 

  /**
   * @brief DFIファイルを出力する
   * @param [in] dfi_name  DFIファイル名
   * @return true:出力成功 false:出力失敗
   */
  bool Write_Proc_File(const std::string dfi_name); 

  /**
   * @brief DFIファイルを出力する
   * @param [in] dfi_name   DFIファイル名
   * @param [in] out_domain 出力Domain
   * @param [in] out_mpi    出力MPI
   * @param [in] out_RankInfo 出力Process
   * @return true:出力成功 false:出力失敗
   */
  static bool Write_Proc_File(const std::string dfi_name, cio_Domain out_domain, 
                              cio_MPI out_mpi, vector<cio_Rank> out_RankInfo); 

  /**
   * @brief DFIファイル:FileInfo要素を出力する
   * @param [in] fp      ファイルポインタ
   * @param [in] tab     インデント
   * @param [in] prefix  ファイル接頭文字
   * @return true:出力成功 false:出力失敗
   */
   bool Write_FileInfo(FILE* fp, const unsigned tab, const std::string prefix);

  /**
   * @brief DFIファイル:Unit要素を出力する
   * @param [in] fp      ファイルポインタ
   * @param [in] tab     インデント
   * @param [in] prefix  ファイル接頭文字
   * @return true:出力成功 false:出力失敗
   */
   bool Write_Unit(FILE* fp, const unsigned tab, const std::string prefix);

  /**
   * @brief DFIファイル:TimeSlice要素を出力する
   * @param [in] fp      ファイルポインタ
   * @param [in] tab     インデント
   * @param [in] step    ステップ数
   * @param [in] time    時間　　
   * @param [in] mio     出力時の分割指定　 true = local / false = gather
   * @return true:出力成功 false:出力失敗
   */
   bool Write_TimeSlice(FILE* fp, const unsigned tab, const unsigned step, REAL_TYPE time,
                        REAL_TYPE* minmax);

  /**
   * @brief Tab(space２つ)を出力する
   * @param [in] fp      ファイルポインタ
   * @param [in] tab     インデント数
   */
  void Write_Tab(FILE* fp, const unsigned tab);


  /**
   * @brief DFIファイル:出力ファイル情報要素を出力する
   * @param [in] fp      ファイルポインタ
   * @param [in] prefix  ファイル接頭文字
   * @param [in] tab     インデント
   * @param [in] step    ステップ数
   * @param [in] mio     出力時の分割指定　 true = local / false = gather
   * @return true:出力成功 false:出力失敗
   */
  bool Write_OutFileInfo(FILE* fp, const unsigned tab, const std::string prefix, const unsigned step, const bool mio);


  /**
   * @brief DFIファイル:BaseName要素を出力する
   * @param [in] fp      ファイルポインタ
   * @param [in] tab     インデント
   * @param [in] dirpath ディレクトリ　　
   */
 void Write_DirectoryPath(FILE* fp, const unsigned tab, const std::string dirpath); 

  /**
   * @brief DFIファイル:BaseName要素を出力する
   * @param [in] fp      ファイルポインタ
   * @param [in] tab     インデント
   * @param [in] prefix  ファイル接頭文字
   */
 void Write_BaseName(FILE* fp, const unsigned tab, const std::string prefix); 

  /**
   * @brief DFIファイル:ファイルフォーマット要素を出力する
   * @param [in] fp      ファイルポインタ
   * @param [in] tab     インデント
   */
  void Write_FileFormat(FILE* fp, const unsigned tab); 

  /**
   * @brief DFIファイル:ガイドセル要素を出力する
   * @param [in] fp     ファイルポインタ
   * @param [in] tab    インデント
   */
  void Write_GuideCell(FILE* fp, const unsigned tab);

  /**
   * @brief DFIファイル:データタイプを出力する
   * @param [in] fp     ファイルポインタ
   * @param [in] tab    インデント
   */
  void Write_DataType(FILE* fp, const unsigned tab);

  /**
   * @brief DFIファイル:Endianを出力する
   * @param [in] fp     ファイルポインタ
   * @param [in] tab    インデント
   */
  void Write_Endian(FILE* fp, const unsigned tab);

  /**
   * @brief DFIファイル:ArrayShapeを出力する
   * @param [in] fp     ファイルポインタ
   * @param [in] tab    インデント
   */
  void Write_ArrayShape(FILE* fp, const unsigned tab);

  /**
   * @brief DFIファイル:Componentを出力する
   * @param [in] fp     ファイルポインタ
   * @param [in] tab    インデント
   */
  void Write_Component(FILE* fp, const unsigned tab);

  /**
   * @brief DFIファイル:Processを出力する
   * @param [in] fp     ファイルポインタ
   * @param [in] tab    インデント
   */
  void Write_Process(FILE* fp, const unsigned tab);

  /**
   * @brief DFIファイル:Lengthを出力する
   * @param [in] fp     ファイルポインタ
   * @param [in] tab    インデント
   */
  void Write_Length(FILE* fp, const unsigned tab);

  /**
   * @brief DFIファイル:L0を出力する
   * @param [in] fp     ファイルポインタ
   * @param [in] tab    インデント
   */
  void Write_L0(FILE* fp, const unsigned tab);

  /**
   * @brief DFIファイル:Velocityを出力する
   * @param [in] fp     ファイルポインタ
   * @param [in] tab    インデント
   */
  void Write_Velocity(FILE* fp, const unsigned tab);

  /**
   * @brief DFIファイル:V0を出力する
   * @param [in] fp     ファイルポインタ
   * @param [in] tab    インデント
   */
  void Write_V0(FILE* fp, const unsigned tab);

  /**
   * @brief DFIファイル:Pressureを出力する
   * @param [in] fp     ファイルポインタ
   * @param [in] tab    インデント
   */
  void Write_Pressure(FILE* fp, const unsigned tab);

  /**
   * @brief DFIファイル:P0を出力する
   * @param [in] fp     ファイルポインタ
   * @param [in] tab    インデント
   */
  void Write_P0(FILE* fp, const unsigned tab);

  /**
   * @brief DFIファイル:P0を出力する
   * @param [in] fp     ファイルポインタ
   * @param [in] tab    インデント
   */
  void Write_DiffPrs(FILE* fp, const unsigned tab);

  /**
   * @brief DFIファイル:Temperaturを出力する
   * @param [in] fp     ファイルポインタ
   * @param [in] tab    インデント
   */
  void Write_Temperatur(FILE* fp, const unsigned tab);

  /**
   * @brief DFIファイル:BaseTempを出力する
   * @param [in] fp     ファイルポインタ
   * @param [in] tab    インデント
   */
  void Write_BaseTemp(FILE* fp, const unsigned tab);

  /**
   * @brief DFIファイル:DiffTempを出力する
   * @param [in] fp     ファイルポインタ
   * @param [in] tab    インデント
   */
  void Write_DiffTemp(FILE* fp, const unsigned tab);

  /**
   * @brief DFIファイル:Stepを出力する
   * @param [in] fp     ファイルポインタ
   * @param [in] tab    インデント
   * @param [in] step   step番号
   */
  void Write_Step(FILE* fp, const unsigned tab, int step);

  /**
   * @brief DFIファイル:Timeを出力する
   * @param [in] fp     ファイルポインタ
   * @param [in] tab    インデント
   * @param [in] time   time 
   */
  void Write_Time(FILE* fp, const unsigned tab, REAL_TYPE time);

  /**
   * @brief DFIファイル:成分を出力する
   * @param [in] fp       ファイルポインタ
   * @param [in] tab      インデント
   * @param [in] compname 成分名
   * @param [in] comp     出力成分 
   */
  void Write_Comp(FILE* fp, const unsigned tab, const std::string compname,
                  REAL_TYPE comp);

  /**
   * @brief DFIファイル:Domainを出力する
   * @param [in] fp     ファイルポインタ
   * @param [in] tab    インデント
   * @return true:出力成功 false:出力失敗
   */
  bool Write_Domain(FILE* fp, const unsigned tab);

  /**
   * @brief DFIファイル:Domainを出力する
   * @param [in] fp         ファイルポインタ
   * @param [in] tab        インデント
   * @param [in] out_domain 出力Domain
   * @return true:出力成功 false:出力失敗
   */
  static bool Write_Domain(FILE* fp, const unsigned tab, cio_Domain out_domain);

  /**
   * @brief DFIファイル:MPIを出力する
   * @param [in] fp     ファイルポインタ
   * @param [in] tab    インデント
   * @return true:出力成功 false:出力失敗
   */
  bool Write_MPI(FILE* fp, const unsigned tab);

  /**
   * @brief DFIファイル:MPIを出力する
   * @param [in] fp      ファイルポインタ
   * @param [in] tab     インデント
   * @param [in] out_mpi 出力MPI
   * @return true:出力成功 false:出力失敗
   */
  static bool Write_MPI(FILE* fp, const unsigned tab, cio_MPI out_mpi);

  /**
   * @brief DFIファイル:Processを出力する
   * @param [in] fp     ファイルポインタ
   * @param [in] tab    インデント
   * @return true:出力成功 false:出力失敗
   */
  bool Write_Process_Rank(FILE* fp, const unsigned tab);

  /**
   * @brief DFIファイル:Processを出力する
   * @param [in] fp           ファイルポインタ
   * @param [in] tab          インデント
   * @param [in] out_RankInfo 出力Process
   * @return true:出力成功 false:出力失敗
   */
  static bool Write_Process_Rank(FILE* fp, const unsigned tabi, vector<cio_Rank> out_RankInfo);

  /**
   * @brief DFIファイル:Processを出力する
   * @param [in] fp     ファイルポインタ
   * @param [in] tab    インデント
   * @return true:出力成功 false:出力失敗
   */
  bool Write_Rank(FILE* fp, const unsigned tab, const int n);

  /**
   * @brief DFIファイル:Processを出力する
   * @param [in] fp     ファイルポインタ
   * @param [in] tab    インデント
   * @param [in] rank   出力Process
   * @return true:出力成功 false:出力失敗
   */
  static bool Write_Rank(FILE* fp, const unsigned tab, cio_Rank rank);

  /**
   * @brief DFIファイル:Originを出力する
   * @param [in] fp     ファイルポインタ
   * @param [in] tab    インデント
   */
  void Write_Origin(FILE* fp, const unsigned tab);

  /**
   * @brief DFIファイル:Originを出力する
   * @param [in] fp     ファイルポインタ
   * @param [in] tab    インデント
   * @param [in] org[3] 出力Origin
   */
  static void Write_Origin(FILE* fp, const unsigned tab, REAL_TYPE org[3]);

  /**
   * @brief DFIファイル:Regionを出力する
   * @param [in] fp     ファイルポインタ
   * @param [in] tab    インデント
   */
  void Write_Region(FILE* fp, const unsigned tab);

  /**
   * @brief DFIファイル:Regionを出力する
   * @param [in] fp        ファイルポインタ
   * @param [in] tab       インデント
   * @param [in] Region[3] 出力Region
   */
  static void Write_Region(FILE* fp, const unsigned tab, REAL_TYPE Region[3]);

  /**
   * @brief DFIファイル:ノード番号要素を出力する
   * @param [in] fp      ファイルポインタ
   * @param [in] tab     インデント
  */
  void Write_MyID(FILE* fp, const unsigned tab, const int n);

  /**
   * @brief DFIファイル:ノード数要素を出力する
   * @param [in] fp  ファイルポインタ
   * @param [in] tab インデント
   */
  void Write_NodeNum(FILE* fp, const unsigned tab); 

  /**
   * @brief DFIファイル:ノード数要素を出力する
   * @param [in] fp           ファイルポインタ
   * @param [in] tab          インデント
   * @param [in] NumberOfRank ノード数
   */
  static void Write_NodeNum(FILE* fp, const unsigned tab, int NumberOfRank); 

  /**
   * @brief DFIファイル:全体ボクセルサイズ要素を出力する
   * @param [in] fp      ファイルポインタ
   * @param [in] tab     インデント
   */
  void Write_WholeSize(FILE* fp, const unsigned tab); 

  /**
   * @brief DFIファイル:全体ボクセルサイズ要素を出力する
   * @param [in] fp     　　    ファイルポインタ
   * @param [in] tab    　　    インデント
   * @param [in] GlobalVoxel[3] 全体ボクセルサイズ
   */
  static void Write_WholeSize(FILE* fp, const unsigned tab, int GlobalVoxel[3]); 

  /**
   * @brief DFIファイル:I,J,K分割数要素を出力する
   * @param [in] fp      ファイルポインタ
   * @param [in] tab     インデント
   */
  void Write_NumDivDomain(FILE* fp, const unsigned tab); 

  /**
   * @brief DFIファイル:I,J,K分割数要素を出力する
   * @param [in] fp                ファイルポインタ
   * @param [in] tab               インデント
   * @param [in] GlobalDivision[3] 分割数
   */
  static void Write_NumDivDomain(FILE* fp, const unsigned tab, int GlobalDivision[3]); 

  /**
   * @brief DFIファイル:IDを出力する
   * @param [in] fp  ファイルポインタ
   * @param [in] tab インデント
   * @param [in] n   ランク番号     
  */
  void Write_ID(FILE* fp, const unsigned tab, const int n);

  /**
   * @brief DFIファイル:IDを出力する
   * @param [in] fp  ファイルポインタ
   * @param [in] tab インデント
   * @param [in] n   ランク番号     
  */
  static void Write_RankID(FILE* fp, const unsigned tab, const int n);

  /**
   * @brief DFIファイル:Hostnameを出力する
   * @param [in] fp  ファイルポインタ
   * @param [in] tab インデント
   * @param [in] n   ランク番号
  */
  void Write_Hostname(FILE* fp, const unsigned tab, const int n);

  /**
   * @brief DFIファイル:Hostnameを出力する
   * @param [in] fp       ファイルポインタ
   * @param [in] tab      インデント
   * @param [in] hostname ホスト名
  */
  static void Write_Hostname(FILE* fp, const unsigned tab, string hostname);

  /**
   * @brief DFIファイル:VoxelSizeを出力する
   * @param [in] fp  ファイルポインタ
   * @param [in] tab インデント
   * @param [in] n   ランク番号
  */
  void Write_L_VoxelSize(FILE* fp, const unsigned tab, const int n);

  /**
   * @brief DFIファイル:VoxelSizeを出力する
   * @param [in] fp           ファイルポインタ
   * @param [in] tab          インデント
   * @param [in] VoxelSize[3] VoxelSize
  */
  static void Write_L_VoxelSize(FILE* fp, const unsigned tab, int VoxelSize[3]);

  /**
   * @brief DFIファイル:HeadIndexを出力する
   * @param [in] fp  ファイルポインタ
   * @param [in] tab インデント
   * @param [in] n   ランク番号
  */
  void Write_HeadIndex(FILE* fp, const unsigned tab, const int n);

  /**
   * @brief DFIファイル:HeadIndexを出力する
   * @param [in] fp　　　　   ファイルポインタ
   * @param [in] tab          インデント
   * @param [in] HeadIndex[3] HeadIndex
  */
  static void Write_HeadIndex(FILE* fp, const unsigned tab, int HeadIndex[3]);

  /**
   * @brief DFIファイル:HeadIndexを出力する
   * @param [in] fp  ファイルポインタ
   * @param [in] tab インデント
   * @param [in] n   ランク番号
  */
  void Write_TailIndex(FILE* fp, const unsigned tab, const int n);

  /**
   * @brief DFIファイル:HeadIndexを出力する
   * @param [in] fp           ファイルポインタ
   * @param [in] tab          インデント
   * @param [in] TailIndex[3] TailIndex
  */
  static void Write_TailIndex(FILE* fp, const unsigned tab, int TailIndex[3]);

  /**
   * 
   *
   */
  bool dbwrite(int RankID);

};

#endif // _cio_DFI_H_
