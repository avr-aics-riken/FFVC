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
// 出力ディレクトリ名を作成する
std::string DFI::Generate_DirName(const std::string prefix, const unsigned m_step)
{
  if ( prefix.empty() ) return NULL;
  
  int len = prefix.size() + 11; // postfix(10) + 1(\0)
  char* tmp = new char[len];
  memset(tmp, 0, sizeof(char)*len);
  
  sprintf(tmp, "%s_%010d", prefix.c_str(), m_step);
  
  std::string fname(tmp);
  if ( tmp ) delete [] tmp;
  
  return fname;
}


// #################################################################
// 出力DFIファイル名を作成する
std::string DFI::Generate_DFI_Name(const std::string prefix)
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
 * @param [in] m_step ステップ数
 * @param [in] m_id   ランク番号
 * @param [in] mio    出力時の分割指定　 local / gather
 */
std::string DFI::Generate_FileName(const std::string prefix, const unsigned m_step, const int m_id, const char* mio)
{
  if ( prefix.empty() ) return NULL;
  
  int len = prefix.size() + 24; // step(10) + id(9) + postfix(4) + 1(\0)
  char* tmp = new char[len];
  memset(tmp, 0, sizeof(char)*len);
  
  // local出力が指定された場合、分割出力
  if ( !strcasecmp(mio, "local") )
  {
    sprintf(tmp, "%s_%010d_id%06d.%s", prefix.c_str(), m_step, m_id, "sph");
  }
  else
  {
    sprintf(tmp, "%s_%010d.%s", prefix.c_str(), m_step, "sph");
  }
  
  std::string fname(tmp);
  if ( tmp ) delete [] tmp;
  
  return fname;
}


// #################################################################
// ファイル名を作成する（拡張子自由）
std::string DFI::Generate_FileName_Free(const std::string prefix, const std::string xxx, const unsigned m_step, const int m_id, const bool mio)
{
  if ( prefix.empty() ) return NULL;
  if ( xxx.empty() ) return NULL;
  
  int len = prefix.size() + xxx.size() + 18; // id(7) + step(10) + 1(.拡張子) + 1(\0)
  
  char* tmp = new char[len];
  memset(tmp, 0, sizeof(char)*len);
  
  // local出力が指定された場合、分割出力
  if ( mio )
  {
    // FieldView がランク番号+ステップ数の記述のため
    sprintf(tmp, "%s_%06d_%010d.%s", prefix.c_str(), m_id, m_step, xxx.c_str());
  }
  else
  {
    sprintf(tmp, "%s_%010d.%s", prefix.c_str(), m_step, xxx.c_str());
  }
  
  std::string fname(tmp);
  if ( tmp ) delete [] tmp;
  
  return fname;
}


// #################################################################
// 初期化
bool DFI::init(const int* g_size, const int* m_div, const int gc, const int stype, const int* hidx, const int* tidx, const std::string m_host)
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
// DFIファイル:BaseName要素を出力する
void DFI::Write_BaseName(FILE* fp, const unsigned tab, const std::string prefix)
{
  Write_Tab(fp, tab);
  fprintf(fp, "Prefix = \"%s\"\n", prefix.c_str());
}



// #################################################################
// データをファイルに書き込む
bool DFI::Write_DFI_File(const std::string prefix, const unsigned step, int& dfi_mng, const bool mio)
{
  if ( prefix.empty() ) return NULL;
  
  // master node only
  int mm;
  MPI_Comm_rank(MPI_COMM_WORLD, &mm);
  
  if ( mm != 0 ) return false;
  
  std::string dfi_name;
  
  if ( mio )
  {
    dfi_name = Generate_DFI_Name(prefix);
    
    if( dfi_name.empty() )
    {
      return false;
    }
    
    if( !Write_File(dfi_name, prefix, step, dfi_mng, mio) )
    {
      return false;
    }
  }
  
  return true;
}



// #################################################################
// DFIファイルを出力する
bool DFI::Write_File(const std::string dfi_name, const std::string prefix, const unsigned step, int& dfi_mng, const char* mio)
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

  
  if ( (dfi_mng == 0) || !flag || (start_type == coarse_restart) ) // カウントゼロ=セッションの開始、または既存ファイルが存在しない、または粗格子リスタート
  {
    if( !(fp = fopen(dfi_name.c_str(), "w")) )
    {
      fprintf(stderr, "Can't open file.(%s)\n", dfi_name.c_str());
      return false;
    }
    
    if (fp) fprintf(fp, "DistributedFileInfo {\n");
    if (fp) fprintf(fp, "\n");
    
    if( !Write_Header(fp, 0, prefix) )
    {
      if (fp) fclose(fp);
      return false;
    }
    
    if (fp) fprintf(fp, "\n");
    
    if (fp) Write_Tab(fp, 1);
    if (fp) fprintf(fp, "FileInfo {\n");

    if ( !Write_OutFileInfo(fp, 1, prefix, step) )
    {
      if (fp) fclose(fp);
      return false;
    }
    
    if (fp) Write_Tab(fp, 1);
    if (fp) fprintf(fp, "}\n");
    if (fp) fprintf(fp, "}\n");
    if (fp) fclose(fp);
    
  }
  else // 既存ファイルが存在する、あるいはセッションが始まり既に書き込み済み >> 追記
  {
    
    // ファイルの内容をバッファ
    fp = fopen(dfi_name.c_str(), "r");
    
    std::string str;
    while( !feof(fp) ){
      int c = fgetc(fp);
      if( !feof(fp) ) str += c;
    }
    fclose(fp);
    
    register int i = str.size() - 1;
    while( --i > 0 ) {
      if( str[i] == '\n' ) { str[i+1] = '\0'; break; }
    }
    
    while( --i > 0 ) {
      if( str[i] == '\n' ) { str[i+1] = '\0'; break; }
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
    
    if ( !Write_OutFileInfo(fp, 1, prefix, step) )
    {
      if (fp) fclose(fp); 
      return false;
    }
    
    if (fp) Write_Tab(fp, 1);
    if (fp) fprintf(fp, "}\n");
    if (fp) fprintf(fp, "}\n");
    if (fp) fclose(fp);
    
  }
  
  dfi_mng++;
  
  return true;
}



// #################################################################
// DFIファイル:ファイルフォーマット要素を出力する
void DFI::Write_FileFormat(FILE* fp, const unsigned tab)
{
  Write_Tab(fp, tab);
  fprintf(fp, "FileFormat = \"sph\"\n");
}



// #################################################################
// DFIファイル:ガイドセル要素を出力する
void DFI::Write_GuideCell(FILE* fp, const unsigned tab)
{
  Write_Tab(fp, tab);
  fprintf(fp, "GuideCell = %d\n", guide);
}



// #################################################################
// DFIファイル:ヘッダー要素を出力する
bool DFI::Write_Header(FILE* fp, const unsigned tab, const std::string prefix)
{
  Write_BaseName(fp, tab+1, prefix);
  if (fp) fprintf(fp, "\n");
  
  Write_MyID(fp, tab+1);
  if (fp) fprintf(fp, "\n");
  
  Write_NodeNum(fp, tab+1);
  if (fp) fprintf(fp, "\n");
  
  Write_WholeSize(fp, tab+1);
  if (fp) fprintf(fp, "\n");
  
  Write_NumDivDomain(fp, tab+1);
  if (fp) fprintf(fp, "\n");
  
  Write_FileFormat(fp, tab+1);
  if (fp) fprintf(fp, "\n");
  
  Write_GuideCell(fp, tab+1);
  if (fp) fprintf(fp, "\n");
  
  if( !Write_NodeInfo(fp, tab+1, prefix) ) return false;
  
  return true;
}



// #################################################################
// DFIファイル:ノード番号要素を出力する
void DFI::Write_MyID(FILE* fp, const unsigned tab)
{
  Write_Tab(fp, tab);
  fprintf(fp, "RankIDinMPIworld = %d\n", my_id);
  
  Write_Tab(fp, tab);
  fprintf(fp, "GroupIDinMPIworld = %d\n", my_id);
}



// #################################################################
// DFIファイル:ボクセル情報要素を出力する
bool DFI::Write_Node(FILE* fp, const unsigned tab, const int n, const std::string prefix)
{
  Write_Tab(fp, tab); 
  fprintf(fp, "Node[@] {\n");
  
  // ID
  Write_Tab(fp, tab+1);
  fprintf(fp, "RankID = %d\n", n);
  
  // Hostname
  Write_Tab(fp, tab+1);
  fprintf(fp, "HostName = \"%s\"\n", hostname.c_str());
  
  // Voxel Size
  Write_Tab(fp, tab+1); 
  fprintf(fp, "VoxelSize = (%d, %d, %d)\n",
          tail[3*n+0] - head[3*n+0] + 1,
          tail[3*n+1] - head[3*n+1] + 1,
          tail[3*n+2] - head[3*n+2] + 1);
  
  // Head Index
  Write_Tab(fp, tab+1); 
  fprintf(fp, "HeadIndex = (%d, %d, %d)\n",
          head[3*n+0],
          head[3*n+1],
          head[3*n+2]);
  
  // Tail Index
  Write_Tab(fp, tab+1); 
  fprintf(fp, "TailIndex = (%d, %d, %d)\n",
          tail[3*n+0],
          tail[3*n+1],
          tail[3*n+2]);
  
  Write_Tab(fp, tab);
  fprintf(fp, "}\n");
  
  return true;
}


// #################################################################
// DFIファイル:ノード情報要素を出力する
bool DFI::Write_NodeInfo(FILE* fp, const unsigned tab, const std::string prefix)
{
  if (fp)
  {
    Write_Tab(fp, tab);
    fprintf(fp, "NodeInfo {\n");
  }
  
  for (int n=0; n<Num_Node; n++) {
    if ( !Write_Node(fp, tab+1, n, prefix) ) return false;
  }
  
  if (fp)
  {
    Write_Tab(fp, tab);
    fprintf(fp, "}\n");
  }
  
  return true;
}


// #################################################################
// DFIファイル:ノード数要素を出力する
void DFI::Write_NodeNum(FILE* fp, const unsigned tab)
{
  Write_Tab(fp, tab);
  fprintf(fp, "NumberOfRankInMPIworld = %d\n", Num_Node);
  
  Write_Tab(fp, tab);
  fprintf(fp, "NumberOfGroupInMPIworld = %d\n", 1);
}



// #################################################################
// DFIファイル:I,J,K分割数要素を出力する
void DFI::Write_NumDivDomain(FILE* fp, const unsigned tab)
{
  Write_Tab(fp, tab);
  fprintf(fp, "GlobalDivision = (%d, %d, %d)\n", div_domain[0], div_domain[1], div_domain[2]);
}



// #################################################################
/**
 * @brief DFIファイル:出力ファイル情報要素を出力する
 * @param [in] fp      ファイルポインタ
 * @param [in] tab     インデント
 * @param [in] prefix  ファイル接頭文字
 * @param [in] step    ステップ数
 */
bool DFI::Write_OutFileInfo(FILE* fp, const unsigned tab, const std::string prefix, const unsigned step)
{
  if (fp)
  {
    Write_Tab(fp, tab+1);
    fprintf(fp, "Step[@] = %d\n", step);
  }
  
  return true;
}


// #################################################################
/**
 * @brief DFIファイル:ファイル名要素を出力する
 * @param [in] fp     ファイルポインタ
 * @param [in] tab    インデント
 * @param [in] prefix ファイル接頭文字
 * @param [in] step   ステップ数
 * @param [in] id     対象ノードID
 * @param [in] mio    出力時の分割指定　 local / gather
 */
bool DFI::Write_OutFileName(FILE* fp, const unsigned tab, const std::string prefix, const unsigned step, const int id, const char* mio)
{
  char fname[FB_FILE_PATH_LENGTH];
  memset(fname, 0, sizeof(char)*FB_FILE_PATH_LENGTH);
  std::string tmp = Generate_FileName(prefix, step, id, mio);

  //if( !path_util::GetFullPathName(tmp.c_str(), fname, FB_FILE_PATH_LENGTH) ) {
  //  return false;
  //}
  strcpy(fname, tmp.c_str());

  if (fp) Write_Tab(fp, tab+1);
  if (fp) fprintf(fp, "FileName[@] = \"%s\"\n", fname);

  return true;
}


// #################################################################
// DFIファイル:ステップ数を出力する
void DFI::Write_Step(FILE* fp, const unsigned tab, const unsigned step)
{
  Write_Tab(fp, tab+1);
  fprintf(fp, "Step = %d\n", step);
}


// #################################################################
// Tab(space２つ)を出力する
void DFI::Write_Tab(FILE* fp, const unsigned tab)
{
  for(int n=0; n<tab; n++) fprintf(fp, "  ");
}


// #################################################################
// DFIファイル:全体ボクセルサイズ要素を出力する
void DFI::Write_WholeSize(FILE* fp, const unsigned tab)
{
  Write_Tab(fp, tab);
  fprintf(fp, "GlobalVoxel = (%d, %d, %d)\n", Gsize[0], Gsize[1], Gsize[2]);
}



// #################################################################
// データをファイルに書き込む
bool DFI::Write_DFI_File(const std::string prefix, const unsigned step, const double time, int& dfi_mng, const bool mio)
{
  if ( prefix.empty() ) return NULL;
  
  // master node only
  int mm;
  MPI_Comm_rank(MPI_COMM_WORLD, &mm);
  
  if ( mm != 0 ) return false;
  
  std::string dfi_name;
  
  if ( mio )
  {
    dfi_name = Generate_DFI_Name(prefix);
    
    if( dfi_name.empty() )
    {
      return false;
    }
    
    if( !Write_File(dfi_name, prefix, step, time, dfi_mng, mio) )
    {
      return false;
    }
  }
  
  return true;
}


// #################################################################
// DFIファイルを出力する
bool DFI::Write_File(const std::string dfi_name, const std::string prefix, const unsigned step, const double time, int& dfi_mng, const bool mio)
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
  
  
  if ( (dfi_mng == 0) || !flag )// || (start_type == coarse_restart) ) // カウントゼロ=セッションの開始、または既存ファイルが存在しない、または粗格子リスタート
  {
    if( !(fp = fopen(dfi_name.c_str(), "w")) )
    {
      fprintf(stderr, "Can't open file.(%s)\n", dfi_name.c_str());
      return false;
    }
    
    if (fp) fprintf(fp, "DistributedFileInfo {\n");
    if (fp) fprintf(fp, "\n");
    
    if( !Write_Header(fp, 0, prefix) )
    {
      if (fp) fclose(fp);
      return false;
    }
    
    if (fp) fprintf(fp, "\n");
    
    if (fp) Write_Tab(fp, 1);
    if (fp) fprintf(fp, "FileInfo {\n");
    
    if ( !Write_OutFileInfo(fp, 1, prefix, step, time, mio) )
    {
      if (fp) fclose(fp);
      return false;
    }
    
    if (fp) Write_Tab(fp, 1);
    if (fp) fprintf(fp, "}\n");
    if (fp) fprintf(fp, "}\n");
    if (fp) fclose(fp);
    
  }
  else // 既存ファイルが存在する、あるいはセッションが始まり既に書き込み済み >> 追記
  {
    
    // ファイルの内容をバッファ
    fp = fopen(dfi_name.c_str(), "r");
    
    std::string str;
    while( !feof(fp) ){
      int c = fgetc(fp);
      if( !feof(fp) ) str += c;
    }
    fclose(fp);
    
    register int i = str.size() - 1;
    while( --i > 0 ) {
      if( str[i] == '\n' ) { str[i+1] = '\0'; break; }
    }
    
    while( --i > 0 ) {
      if( str[i] == '\n' ) { str[i+1] = '\0'; break; }
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
    
    if ( !Write_OutFileInfo(fp, 1, prefix, step, time, mio) )
    {
      if (fp) fclose(fp);
      return false;
    }
    
    if (fp) Write_Tab(fp, 1);
    if (fp) fprintf(fp, "}\n");
    if (fp) fprintf(fp, "}\n");
    if (fp) fclose(fp);
    
  }
  
  dfi_mng++;
  
  return true;
}


// #################################################################
// DFIファイル:出力ファイル情報要素を出力する
bool DFI::Write_OutFileInfo(FILE* fp, const unsigned tab, const std::string prefix, const unsigned step, const double time, const bool mio)
{
  if (fp)
  {
    Write_Tab(fp, tab+1);
    fprintf(fp, "Step[@] = %d\n", step);
  }
  
  if (fp)
  {
    Write_Tab(fp, tab+1);
    fprintf(fp, "Time[@] = %f\n", time);
  }
  
  return true;
}
