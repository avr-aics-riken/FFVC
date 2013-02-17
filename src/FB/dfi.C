// #################################################################
//
// CAERU Library
//
// Copyright (c) 2012-2013  All right reserved.
//
// Institute of Industrial Science, University of Tokyo, Japan. 
//
// #################################################################

//@file   dfi.C
//@brief  DFIファイル生成
//@author kero


#include "dfi.h"
#include "util_Path.h"


// #################################################################
/**
 * @brief マシンのエンディアンを調べる
 * @ret true-Big, false-Little
 */
bool DFI::chekcEndian()
{
	int a = 1;
  bool ret;
  
	if (*((char *)&a)) {
		ret = false;
	} else if (*((char *)&a + (sizeof(int) - 1))) {
		ret = true;
	} else {
		printf("Error : Unknown Endian\n");
    Exit(-1);
	}
  
  return ret;
}


// #################################################################
/**
 * @brief ディレクトリ名を作成する
 * @param [in] path      ディレクトリ名
 * @param [in] m_step    ステップ数
 * @param [in] slice     時系列出力モード (ON / OFF)
 */
std::string DFI::GenerateDirName(const std::string path, const unsigned m_step, const int slice)
{
  char digit[11];
  memset(digit, 0, sizeof(char)*11);
  std::string tmp;
  
  if (slice == OFF)
  {
    tmp = path + "/";
  }
  else
  {
    sprintf(digit, "%010u", m_step);
    tmp = path + "/" + digit + "/";
  }
  
  std::string fname(tmp);
  
  return fname;
}


// #################################################################
/**
 * @brief 出力DFIファイル名を作成する
 * @param [in] prefix ファイル接頭文字
 */
std::string DFI::GenerateDFIname(const std::string prefix)
{
  if ( prefix.empty() ) return NULL;
  
  int len = prefix.size() + 5; // postfix(4) + 1(\0)
  char* tmp = new char[len];
  memset(tmp, 0, sizeof(char)*len);
  
  sprintf(tmp, "%s.%s", prefix.c_str(), "dfi");
  
  std::string fname(tmp);
  if ( tmp ) delete [] tmp;
  
  return fname;
}



// #################################################################
/**
 * @brief ファイル名を作成する
 * @param [in] prefix ファイル接頭文字
 * @param [in] fmt    ファイルフォーマット拡張子, 0:"sph", 2:"dat", 1:others(plot3dとして認識)
 * @param [in] m_step ステップ数
 * @param [in] m_id   ランク番号
 * @param [in] divide 出力時の分割指定　 single:false / divide:true
 */
std::string DFI::GenerateFileName(const std::string prefix, const std::string fmt, const unsigned m_step, const int m_id, bool divide)
{
  if ( prefix.empty() ) return NULL;
  if ( fmt.empty() ) return NULL;
  
  int len = prefix.size() + fmt.size() + 22; // 1(_) + 10(step) + 3(_id) + 6(rank) + 1(.) + 1(\0)
  char* tmp = new char[len];
  memset(tmp, 0, sizeof(char)*len);
  
  
  if ( divide )
  {
    if ( !strcasecmp(fmt.c_str(), "sph") || !strcasecmp(fmt.c_str(), "dat") )
    {
      sprintf(tmp, "%s_%010u_id%06d.%s", prefix.c_str(), m_step, m_id, fmt.c_str());
    }
    else // PLOT3D:FieldView がランク番号+ステップ数の記述のため
    {
      sprintf(tmp, "%s_%06d_%010u.%s", prefix.c_str(), m_id, m_step, fmt.c_str());
    }
  }
  else
  {
    if ( !strcasecmp(fmt.c_str(), "sph") || !strcasecmp(fmt.c_str(), "dat") )
    {
      sprintf(tmp, "%s_%010u.%s", prefix.c_str(), m_step, fmt.c_str());
    }
    else // PLOT3D:FieldView がランク番号+ステップ数の記述のため
    {
      sprintf(tmp, "%s.%s", prefix.c_str(), fmt.c_str());
    }
  }

  std::string fname(tmp);
  if ( tmp ) delete [] tmp;
  
  return fname;
}


// #################################################################
/**
 * @brief 初期化
 * @param [in] g_size  グローバルサイズ
 * @param [in] m_div   ノード分割数
 * @param [in] gc      ガイドセル
 * @param [in] stype   スタートタイプ
 * @param [in] m_refL  代表長さ
 * @param [in] m_refV  代表速度
 * @param [in] m_BaseP 基準圧力
 * @param [in] m_DiffP 圧力差
 * @param [in] Unit_L  長さの単位
 * @param [in] Unit_V  速度の単位
 * @param [in] Unit_P  圧力の単位
 * @param [in] hidx    開始インデクス
 * @param [in] tidx    終端インデクス
 * @param [in] m_host  ホスト名
 */
bool DFI::init(const int* g_size,
               const int* m_div,
               const int gc,
               const int stype,
               const REAL_TYPE m_refL,
               const REAL_TYPE m_refV,
               const REAL_TYPE m_BaseP,
               const REAL_TYPE m_DiffP,
               const std::string m_UnitL,
               const std::string m_UnitV,
               const std::string m_UnitP,
               const int* hidx,
               const int* tidx,
               const std::string m_host)
{
  MPI_Comm_size(MPI_COMM_WORLD, &Num_Node);
  
  MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
  
  if ( my_id < 0 ) return false;
  
  // global size
  Gsize[0] = g_size[0];
  Gsize[1] = g_size[1];
  Gsize[2] = g_size[2];
  
  // ノード分割数
  div_domain[0] = m_div[0];
  div_domain[1] = m_div[1];
  div_domain[2] = m_div[2];
  
  guide = gc;
  
  start_type = stype;
  
  RefLength   = m_refL;
  RefVelocity = m_refV;
  BasePrs     = m_BaseP;
  DiffPrs     = m_DiffP;
  
  Unit_L = m_UnitL;
  Unit_V = m_UnitV;
  Unit_P = m_UnitP;
  
  head = new int[3*Num_Node];
  tail = new int[3*Num_Node];
  
  hostname = m_host;
  
  for (int i=0; i<Num_Node*3; i++) {
    head[i] = hidx[i];
    tail[i] = tidx[i];
  }
  
  return true;
}



// #################################################################
/**
 * @brief Rankの情報を出力する
 * @param [in] fp      ファイルポインタ
 * @param [in] tab     インデント
 * @param [in] n       対象ノードID
 */
bool DFI::WriteRank(FILE* fp, const unsigned tab, const int n)
{
  WriteTab(fp, tab); 
  fprintf(fp, "Rank[@] {\n");
  
  // ID
  WriteTab(fp, tab+1);
  fprintf(fp, "ID        = %d\n", n);
  
  // Hostname
  WriteTab(fp, tab+1);
  fprintf(fp, "HostName  = \"%s\"\n", hostname.c_str());
  
  // Voxel Size
  WriteTab(fp, tab+1); 
  fprintf(fp, "VoxelSize = (%d, %d, %d)\n",
          tail[3*n+0] - head[3*n+0] + 1,
          tail[3*n+1] - head[3*n+1] + 1,
          tail[3*n+2] - head[3*n+2] + 1);
  
  // Head Index
  WriteTab(fp, tab+1); 
  fprintf(fp, "HeadIndex = (%d, %d, %d)\n",
          head[3*n+0],
          head[3*n+1],
          head[3*n+2]);
  
  // Tail Index
  WriteTab(fp, tab+1); 
  fprintf(fp, "TailIndex = (%d, %d, %d)\n",
          tail[3*n+0],
          tail[3*n+1],
          tail[3*n+2]);
  
  WriteTab(fp, tab);
  fprintf(fp, "}\n");
  
  return true;
}



// #################################################################
/**
 * @brief Tab(space２つ)を出力する
 * @param [in] fp      ファイルポインタ
 * @param [in] tab     インデント数
 */
void DFI::WriteTab(FILE* fp, const unsigned tab)
{
  for(int n=0; n<tab; n++) fprintf(fp, "  ");
}


// #################################################################
/**
 * @brief DFI procファイルを生成する
 * @param [in]  g_org  計算領域の基点
 * @param [in]  g_reg  計算領域の大きさ
 */
bool DFI::WriteDFIproc(const REAL_TYPE* g_org, const REAL_TYPE* g_reg)
{

  // check master node only
  int mm;
  MPI_Comm_rank(MPI_COMM_WORLD, &mm);
  if ( mm != 0 ) return false;

  
  std::string dfi_name = GenerateDFIname("proc");
  
  if ( dfi_name.empty() ) return false;

  
  FILE* fp = NULL;
  
  if ( !(fp = fopen(dfi_name.c_str(), "w")) )
  {
    fprintf(stderr, "Can't open file.(%s)\n", dfi_name.c_str());
    return false;
  }
  
  // Domain {} -------------------------------
  fprintf(fp, "Domain {\n");
  
  
  // GlobalOrigin
  WriteTab(fp, 1);
  fprintf(fp, "GlobalOrigin   = (%e, %e, %e)\n", g_org[0], g_org[1], g_org[2]);
  
  // GlobalRegion
  WriteTab(fp, 1);
  fprintf(fp, "GlobalRegion   = (%e, %e, %e)\n", g_reg[0], g_reg[1], g_reg[2]);
  
  // GlobalVoxel
  WriteTab(fp, 1);
  fprintf(fp, "GlobalVoxel    = (%d, %d, %d)\n", Gsize[0], Gsize[1], Gsize[2]);
  
  // GlobalDivision
  WriteTab(fp, 1);
  fprintf(fp, "GlobalDivision = (%d, %d, %d)\n", div_domain[0], div_domain[1], div_domain[2]);
  
  // end of Domain {}
  fprintf(fp, "}\n");
  
  
  
  // MPI {} -------------------------------
  fprintf(fp, "MPI {\n");
  
  // NumberOfRank
  WriteTab(fp, 1);
  fprintf(fp, "NumberOfRank   = %d\n", Num_Node);
  
  // NumberOfGroup
  WriteTab(fp, 1);
  fprintf(fp, "NumberOfGroup  = %d\n", 1);
  
  // end of MPI {}
  fprintf(fp, "}\n");
  
  
  
  // Process {} -------------------------------
  fprintf(fp, "Process {\n");
  
  
  for (int n=0; n<Num_Node; n++) {
    if ( !WriteRank(fp, 1, n) ) return false;
  }
  
  
  // end of Process {}
  fprintf(fp, "}\n");
  
  if (fp) fclose(fp);
  
  return true;
}


// #################################################################
/**
 * @brief DFI indexファイルを生成する
 * @param [in]     prefix   ファイル接頭文字
 * @param [in]     dir      ディレクトリパス
 * @param [in]     fmt      ファイル拡張子
 * @param [in]     step     ステップ
 * @param [in]     time     時間
 * @param [in,out] dfi_mng  出力管理カウンタ
 * @param [in]     shape    配列の形式 ("nijk" / "ijkn")
 * @param [in]     compo    データの成分数(n)
 * @param [in]     minmax   最小値、最大値
 * @param [in]     mio      出力時の分割指定　 single:false / divide:true
 * @param [in]     avr_mode 平均値出力の場合、false　デフォルトtrue
 * @param [in]     a_step   平均ステップ数
 * @param [in]     a_time   平均時間
 */
bool DFI::WriteDFIindex(const std::string prefix,
                        const std::string dir,
                        const std::string fmt,
                        const unsigned step,
                        const double time,
                        int& dfi_mng,
                        const std::string shape,
                        const int compo,
                        const REAL_TYPE* minmax,
                        const bool mio,
                        bool avr_mode,
                        unsigned a_step,
                        double a_time)
{
  if ( prefix.empty() ) return NULL;
  if ( dir.empty() ) return NULL;
  
  // check master node only
  int mm;
  MPI_Comm_rank(MPI_COMM_WORLD, &mm);
  if ( mm != 0 ) return false;
  
  std::string dfi_name;
  
  dfi_name = GenerateDFIname(prefix);
  
  if ( dfi_name.empty() ) return false;
  
  if ( !WriteIndex(dfi_name, prefix, dir, fmt, step, time, dfi_mng, shape, compo, minmax, avr_mode, a_step, a_time) ) return false;
  
  return true;
}


// #################################################################
/**
 * @brief Indexファイルの内容を書き出す
 * @param [in]     dfi_name  DFIファイル名
 * @param [in]     prefix    ファイル接頭文字
 * @param [in]     dir       ディレクトリパス
 * @param [in]     fmt       ファイル拡張子
 * @param [in]     step      ステップ数
 * @param [in]     time      時間
 * @param [in,out] dfi_mng   出力管理カウンタ
 * @param [in]     shape     配列の形式 ("nijk" / "ijkn")
 * @param [in]     compo     データの成分数(n)
 * @param [in]     minmax    最小値、最大値 
 * @param [in]     avr_mode  平均値出力の場合、false
 * @param [in]     a_step    平均ステップ数
 * @param [in]     a_time    平均時間
 */
bool DFI::WriteIndex(const std::string dfi_name,
                     const std::string prefix,
                     const std::string dir,
                     const std::string fmt,
                     const unsigned step,
                     const double time,
                     int& dfi_mng,
                     const std::string shape,
                     const int compo,
                     const REAL_TYPE* minmax,
                     const bool avr_mode,
                     const unsigned a_step,
                     const double a_time)
{
  if ( dfi_name.empty() ) return false;
  if ( prefix.empty() ) return false;
  
  FILE* fp = NULL;
  
  // File exist ?
  bool flag = false;
  if ( fp = fopen(dfi_name.c_str(), "r") )
  {
    flag = true;
    fclose(fp);
  }


  // カウントゼロ=セッションの開始、または既存ファイルが存在しない、または粗格子リスタート
  if ( (dfi_mng == 0) || !flag || (start_type == restart_refinement) ) 
  {
    if ( !(fp = fopen(dfi_name.c_str(), "w")) )
    {
      fprintf(stderr, "Can't open file.(%s)\n", dfi_name.c_str());
      return false;
    }
    
    // FileInfo {} -------------------------------
    fprintf(fp, "FileInfo {\n");
    
    
    // DirectoryPath
    WriteTab(fp, 1);
    fprintf(fp, "DirectoryPath = \"%s\"\n", dir.c_str());
    
    
    // Prefix
    WriteTab(fp, 1);
    fprintf(fp, "Prefix        = \"%s\"\n", prefix.c_str());
    
    
    // FileFormat
    WriteTab(fp, 1);
    fprintf(fp, "FileFormat    = \"%s\"\n", fmt.c_str());
    
    
    // GuideCell
    WriteTab(fp, 1);
    fprintf(fp, "GuideCell     = %d\n", guide);
    
    
    // DataType
    WriteTab(fp, 1);
    if ( sizeof(REAL_TYPE) == 4 )
    {
      fprintf(fp, "DataType      = \"Float32\"\n");
    }
    else
    {
      fprintf(fp, "DataType      = \"Float64\"\n");
    }
    
    
    // Endian
    WriteTab(fp, 1);
    if ( chekcEndian() )
    {
      fprintf(fp, "Endian        = \"Big\"\n");
    }
    else
    {
      fprintf(fp, "Endian        = \"Little\"\n");
    }
    
    
    // ArrayShape
    WriteTab(fp, 1);
    if ( !strcasecmp(shape.c_str(), "nijk") )
    {
      fprintf(fp, "ArrayShape    = \"nijk\"\n");
    }
    else
    {
      fprintf(fp, "ArrayShape    = \"ijkn\"\n");
    }
    
    
    // Component
    WriteTab(fp, 1);
    fprintf(fp, "Component     = %d\n", compo);
    
    
    // end of FileInfo {}
    fprintf(fp, "}\n\n");
    
    
    
    // Unit {} -------------------------------
    fprintf(fp, "Unit {\n");
    
    
    // Length
    WriteTab(fp, 1);
    fprintf(fp, "Length        = \"%s\"\n", Unit_L.c_str());
    
    // L0
    WriteTab(fp, 1);
    fprintf(fp, "L0            = %e\n", RefLength);
    
    // Velocity
    WriteTab(fp, 1);
    fprintf(fp, "Velocity      = \"%s\"\n", Unit_V.c_str());
    
    // V0
    WriteTab(fp, 1);
    fprintf(fp, "V0            = %e\n", RefVelocity);
    
    // Pressure
    WriteTab(fp, 1);
    fprintf(fp, "Pressure      = \"%s\"\n", Unit_P.c_str());
    
    // P0
    WriteTab(fp, 1);
    fprintf(fp, "P0            = %e\n", BasePrs);
    
    // DiffPrs
    WriteTab(fp, 1);
    fprintf(fp, "DiffPrs       = %e\n", DiffPrs);
    
    // end of Unit {}
    fprintf(fp, "}\n\n");
    
    
    
    // FilePath {} -------------------------------
    fprintf(fp, "FilePath {\n");
    
    
    // proc file
    WriteTab(fp, 1);
    fprintf(fp, "Process       = \"proc.dfi\"\n");
    
    
    // end of FilePath {}
    fprintf(fp, "}\n\n");
    
    
    
    // TimeSlice {} -------------------------------
    fprintf(fp, "TimeSlice {\n");
    
    WriteTimeSlice(fp, 1, prefix, step, time, minmax, avr_mode, a_step, a_time);
    
    // end of TimeSlice {}
    fprintf(fp, "}\n\n");
    
    if (fp) fclose(fp);
 
  }
  else // 既存ファイルが存在する、あるいはセッションが始まり既に書き込み済み >> 追記
  {
    
    // ファイルの内容をバッファ
    fp = fopen(dfi_name.c_str(), "r");
    
    std::string str;
    while ( !feof(fp) ) {
      int c = fgetc(fp);
      if ( !feof(fp) ) str += c;
    }
    fclose(fp);
    
    register int i = str.size() - 1;
    while( --i > 0 ) {
      if ( str[i] == '\n' ) { str[i+1] = '\0'; break; }
    }
    
    while( --i > 0 ) {
      if ( str[i] == '\n' ) { str[i+1] = '\0'; break; }
    }
    
    // 新規ファイルを生成し、バッファを書きだしたあとにファイル情報を追記
    if( !(fp = fopen(dfi_name.c_str(), "w")) )
    {
      fprintf(stderr, "Can't open file.(%s)\n", dfi_name.c_str());
      return false;
    }
    
    if ( fp && fwrite(str.c_str(), sizeof(char), strlen(str.c_str()), fp) != strlen(str.c_str()) )
    {
      if (fp) fclose(fp);
      return false;
    }
    
    WriteTimeSlice(fp, 1, prefix, step, time, minmax, avr_mode, a_step, a_time);
    
    // end of TimeSlice {}
    fprintf(fp, "}\n\n");
    
    if (fp) fclose(fp);
    
  }
  
  dfi_mng++;
  
  return true;
}


// #################################################################
/**
 * @brief 時系列情報を出力
 * @param [in] fp       ファイルポインタ
 * @param [in] tab      インデント
 * @param [in] prefix   ファイル接頭文字
 * @param [in] step     ステップ数
 * @param [in] time     時間
 * @param [in] minmax   最小値、最大値
 * @param [in] avr_mode 平均値出力の場合、false
 * @param [in] a_step   平均ステップ数
 * @param [in] a_time   平均時間
 */
void DFI::WriteTimeSlice(FILE* fp,
                         const unsigned tab,
                         const std::string prefix,
                         const unsigned step,
                         const double time,
                         const REAL_TYPE* minmax,
                         const bool avr_mode,
                         const unsigned a_step,
                         const double a_time)
{
  WriteTab(fp, tab);
  fprintf(fp, "Slice[@] {\n");
  
  WriteTab(fp, tab+1);
  fprintf(fp, "Step         = %u\n", step);
  
  WriteTab(fp, tab+1);
  fprintf(fp, "Time         = %e\n", time);
  
  
  if ( !avr_mode )
  {
    WriteTab(fp, tab+1);
    fprintf(fp, "AveragedStep = %u\n", a_step);
    
    WriteTab(fp, tab+1);
    fprintf(fp, "AveragedTime = %e\n", a_time);
    
  }
  
  WriteTab(fp, tab+1);
  fprintf(fp, "MinMax[@] {\n");
  
  WriteTab(fp, tab+2);
  fprintf(fp, "Min  = %e\n", minmax[0]);
  
  WriteTab(fp, tab+2);
  fprintf(fp, "Max  = %e\n", minmax[1]);
  
  WriteTab(fp, tab+1);
  fprintf(fp, "}\n");
  
  
  // you can write any annotation here
  
  
  WriteTab(fp, tab);
  fprintf(fp, "}\n");
}
