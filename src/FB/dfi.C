/*
 * SPHERE - Skeleton for PHysical and Engineering REsearch
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2004-2012
 *
 */

//@file DFI.C
//@brief DFIファイル生成
//@author keno, Advanced Vis Team, AICS, RIKEN

#include "dfi.h"
#include "util_Path.h"

//@fn void DFI::init()
bool DFI::init(const int* g_size, const int* m_div, const int gc, const int* hidx, const int* tidx)
{
  MPI_Comm_size(MPI_COMM_WORLD, &Num_Node);
  if ( Num_Node < 2 ) {
    return false;
  }

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
  
  head = new int[3*Num_Node];
  tail = new int[3*Num_Node];
  
  hostname = new char*[Num_Node];
  for (int i=0; i<Num_Node; i++) {
    hostname[i] = new char[LABEL];
  }
  
  for (int i=0; i<Num_Node*3; i++) {
    head[i] = hidx[i];
    tail[i] = tidx[i];
  }
  
  return true;
}
           

/**
 * データをファイルに書き込む。
 * @param prefix ファイル接頭文字
 * @param step   ステップ
 * @param mio    出力時の分割指定　 true = local / false = gather
 */
bool DFI::Write_DFI_File(const std::string prefix, const int step, const bool mio)
{
  if ( prefix.empty() ) return NULL;

  // master node only
  int mm;
  MPI_Comm_rank(MPI_COMM_WORLD, &mm);

  if ( mm != 0 ) return false;

  std::string dfi_name;

  if( mio ) {

    dfi_name = Generate_DFI_Name(prefix, my_id);

    if( dfi_name.empty() ) {
      return false;
    }

    if( !Write_File(dfi_name, prefix, step, mio) ) {
      return false;
    }

  }

  return true;
}


/**
 * ファイル名を作成する。
 * @param prefix ファイル接頭文字
 * @param m_step
 * @param m_id 
 * @param mio    出力時の分割指定　 true = local / false = gather(default)
 */
std::string DFI::Generate_FileName(const std::string prefix, const int m_step, const int m_id, const bool mio)
{
  if ( prefix.empty() ) return NULL;
  
  int len = prefix.size() + 24; // step(10) + id(9) + postfix(4) + 1
  char* tmp = new char[len];
  memset(tmp, 0, sizeof(char)*len);
  
  // local出力が指定された場合、分割出力
  if( mio ){
    sprintf(tmp, "%s%010d_id%06d.%s", prefix.c_str(), m_step, m_id, "sph");
  }
  else {
    sprintf(tmp, "%s%010d.%s", prefix.c_str(), m_step, "sph");
  }
  
  std::string fname(tmp);
  if ( tmp ) delete [] tmp;
  
  return fname;
}


/**
 * 出力DFIファイル名を作成する。
 *
 * @param prefix ファイル接頭文字
 * @param m_id   ランク番号
 */
std::string DFI::Generate_DFI_Name(const std::string prefix, const int m_id)
{
  if ( prefix.empty() ) return NULL;
  
  int len = prefix.size() + 14; // id(9) + postfix(4) + 1
  char* tmp = new char[len];
  memset(tmp, 0, sizeof(char)*len);
  
  sprintf(tmp, "%sid%06d.%s", prefix.c_str(), m_id, "dfi");
  
  std::string fname(tmp);
  if ( tmp ) delete [] tmp;
  
  return fname;
}

/**
 * DFIファイルを出力する。
 *
 * DFIファイル=ボクセル分割情報ファイル
 * @param dfi_name  DFIファイル名
 * @param prefix    ファイル接頭文字
 * @param step      ステップ数
 * @param mio    出力時の分割指定　 true = local / false = gather
 */
bool DFI::Write_File(const std::string dfi_name, const std::string prefix, const int step, const bool mio)
{
  if ( dfi_name.empty() ) return false;
  if ( prefix.empty() ) return false;

  FILE* fp = NULL;
  
  // File exist ?
  bool flag = false;
  if ( fp = fopen(dfi_name.c_str(), "r") ) {
    flag = true;
    fclose(fp);
  }
      
  if ( WriteCount == 0 ) { // カウントゼロのとき（リスタート時も） << bug

    if( !(fp = fopen(dfi_name.c_str(), "w")) ) {
      fprintf(stderr, "Can't open file.(%s)\n", dfi_name.c_str());
      return false;
    }
    
    if (fp) fprintf(fp, "<SphereDispersedFileInfo>\n");
    if (fp) fprintf(fp, "\n");
    
    if( !Write_Header(fp, 0, prefix) ) {
      if (fp) fclose(fp);
      return false;
    }
    
    if (fp) fprintf(fp, "\n");
    
    if (fp) Write_Tab(fp, 1);
    if (fp) fprintf(fp, "<Elem name=\"FileInfo\">\n");

    if( !Write_OutFileInfo(fp, 1, prefix, step, mio) ){
      if (fp) fclose(fp);
      return false;
    }
    
    if (fp) Write_Tab(fp, 1);
    if (fp) fprintf(fp, "</Elem>\n");
    if (fp) fprintf(fp, "</SphereDispersedFileInfo>\n");
    if (fp) fclose(fp);
    
  }
  else { // file exist
    
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
    
    if( !(fp = fopen(dfi_name.c_str(), "w")) ) {
      fprintf(stderr, "Can't open file.(%s)\n", dfi_name.c_str());
      return false;
    }
    
    if( fp && fwrite(str.c_str(), sizeof(char), strlen(str.c_str()), fp) != strlen(str.c_str()) ){
      if (fp) fclose(fp);
      return false;
    }
    
    if( !Write_OutFileInfo(fp, 1, prefix, step, mio) ){
      if (fp) fclose(fp); 
      return false;
    }
    
    if (fp) Write_Tab(fp, 1);
    if (fp) fprintf(fp, "</Elem>\n");
    if (fp) fprintf(fp, "</SphereDispersedFileInfo>\n");
    if (fp) fclose(fp);
    
  }
  
  WriteCount++; // bug
  
  return true;
}

/**
 * Tab(space２つ)を出力する。
 * @param fp      ファイルポインタ
 * @param tab     インデント数
 */
void DFI::Write_Tab(FILE* fp, const unsigned tab)
{
  for(int n=0; n<tab; n++) fprintf(fp, "  ");
}

/**
 * DFIファイル:ヘッダー要素を出力する。
 *
 * @param fp      ファイルポインタ
 * @param tab     インデント
 * @param prefix  ファイル接頭文字
 */
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
  
  if( !Write_NodeInfo(fp, tab+1, prefix) ) return false;
  
  return true;
}

/**
 * DFIファイル:BaseName要素を出力する。
 *
 * @param fp      ファイルポインタ
 * @param tab     インデント
 * @param prefix  ファイル接頭文字
 */
void DFI::Write_BaseName(FILE* fp, const unsigned tab, const std::string prefix)
{
  Write_Tab(fp, tab);
  fprintf(fp, "<Param name=\"BaseName\" dtype=\"STRING\" value=\"%s\" />\n", prefix.c_str());
}

/**
 * DFIファイル:ノード番号要素を出力する。
 *
 * @param fp      ファイルポインタ
 * @param tab     インデント
 */
void DFI::Write_MyID(FILE* fp, const unsigned tab)
{
  Write_Tab(fp, tab);
  fprintf(fp, "<Param name=\"WorldID\" dtype=\"INT\" value=\"%d\" />\n", my_id);
  
  Write_Tab(fp, tab);
  fprintf(fp, "<Param name=\"GroupID\" dtype=\"INT\" value=\"%d\" />\n", my_id);
}

/**
 * DFIファイル:ノード数要素を出力する。
 *
 * @param fp  ファイルポインタ
 * @param tab インデント
 */
void DFI::Write_NodeNum(FILE* fp, const unsigned tab)
{
  Write_Tab(fp, tab);
  fprintf(fp, "<Param name=\"WorldNodeNum\" dtype=\"INT\" value=\"%d\" />\n", Num_Node);
  
  Write_Tab(fp, tab);
  fprintf(fp, "<Param name=\"GroupNodeNum\" dtype=\"INT\" value=\"%d\" />\n", Num_Node);
}

/**
 * DFIファイル:全体ボクセルサイズ要素を出力する。
 *
 * @param fp      ファイルポインタ
 * @param tab     インデント
 */
void DFI::Write_WholeSize(FILE* fp, const unsigned tab)
{
  Write_Tab(fp, tab);
  fprintf(fp, "<Elem name=\"WholeVoxelSize\">\n");
  
  Write_Tab(fp, tab+1);
  fprintf(fp, "<Param name=\"I\" dtype=\"INT\" value=\"%d\" />\n", Gsize[0]);
  
  Write_Tab(fp, tab+1);
  fprintf(fp, "<Param name=\"J\" dtype=\"INT\" value=\"%d\" />\n", Gsize[1]);
  
  Write_Tab(fp, tab+1);
  fprintf(fp, "<Param name=\"K\" dtype=\"INT\" value=\"%d\" />\n", Gsize[2]);
  
  Write_Tab(fp, tab);
  fprintf(fp, "</Elem>\n");
}

/**
 * DFIファイル:I,J,K分割数要素を出力する。
 *
 * @param fp      ファイルポインタ
 * @param tab     インデント
 */
void DFI::Write_NumDivDomain(FILE* fp, const unsigned tab)
{
  Write_Tab(fp, tab);
  fprintf(fp, "<Elem name=\"VoxelDivMethod\">\n");
  
  Write_Tab(fp, tab+1);
  fprintf(fp, "<Param name=\"I\" dtype=\"INT\" value=\"%d\" />\n", div_domain[0]);
  
  Write_Tab(fp, tab+1);
  fprintf(fp, "<Param name=\"J\" dtype=\"INT\" value=\"%d\" />\n", div_domain[1]);
  
  Write_Tab(fp, tab+1);
  fprintf(fp, "<Param name=\"K\" dtype=\"INT\" value=\"%d\" />\n", div_domain[2]);
  
  Write_Tab(fp, tab);
  fprintf(fp, "</Elem>\n");
}

/**
 * DFIファイル:ファイルフォーマット要素を出力する。
 *
 * @param fp      ファイルポインタ
 * @param tab     インデント
 */
void DFI::Write_FileFormat(FILE* fp, const unsigned tab)
{
  Write_Tab(fp, tab);
  fprintf(fp, "<Param name=\"Format\" dtype=\"STRING\" value=\"sph");
  fprintf(fp, "\" />\n");
}

/**
 * DFIファイル:ノード情報要素を出力する。
 * @param fp      ファイルポインタ
 * @param tab     インデント
 * @param prefix  ファイル接頭文字
 */
bool DFI::Write_NodeInfo(FILE* fp, const unsigned tab, const std::string prefix)
{
  if (fp) {
    Write_Tab(fp, tab); 
    fprintf(fp, "<Elem name=\"NodeInfo\">\n");
  }
  
  for (int n=0; n<Num_Node; n++){
    if ( !Write_Node(fp, tab+1, n, prefix) ) return false;
  }
  
  if (fp) {
    Write_Tab(fp, tab); 
    fprintf(fp, "</Elem>\n");
  }
  
  return true;
}

/**
 * DFIファイル:ボクセル情報要素を出力する。
 *
 * @param fp      ファイルポインタ
 * @param tab     インデント
 * @param n       対象ノードID
 * @param prefix  ファイル接頭文字
 */
bool DFI::Write_Node(FILE* fp, const unsigned tab, const int n, const std::string prefix)
{
  Write_Tab(fp, tab); 
  fprintf(fp, "<Elem name=\"Node\">\n");
  
  // ID
  Write_Tab(fp, tab+1);
  fprintf(fp, "<Param name=\"GroupID\" dtype=\"INT\" value=\"%d\" />\n", n);
  
  // Hostname
  Write_Tab(fp, tab+1);
  fprintf(fp, "<Param name=\"HostName\" dtype=\"STRING\" value=\"%s\" />\n", hostname[n]);
  
  // VoxelSize
  Write_Tab(fp, tab+1); 
  fprintf(fp, "<Elem name=\"VoxelSize\">\n");
  
  Write_Tab(fp, tab+2);
  fprintf(fp, "<Param name=\"I\" dtype=\"INT\" value=\"%d\" />\n", tail[3*n+0] - head[3*n+0] + 1);
  
  Write_Tab(fp, tab+2);
  fprintf(fp, "<Param name=\"J\" dtype=\"INT\" value=\"%d\" />\n", tail[3*n+1] - head[3*n+1] + 1);
  
  Write_Tab(fp, tab+2);
  fprintf(fp, "<Param name=\"K\" dtype=\"INT\" value=\"%d\" />\n", tail[3*n+2] - head[3*n+2] + 1);
  
  Write_Tab(fp, tab+1); 
  fprintf(fp, "</Elem>\n");
  
  // Head Index
  Write_Tab(fp, tab+1); 
  fprintf(fp, "<Elem name=\"HeadIndex\">\n");
  
  Write_Tab(fp, tab+2);
  fprintf(fp, "<Param name=\"I\" dtype=\"INT\" value=\"%d\" />\n", head[3*n+0]);
  
  Write_Tab(fp, tab+2);
  fprintf(fp, "<Param name=\"J\" dtype=\"INT\" value=\"%d\" />\n", head[3*n+1]);
  
  Write_Tab(fp, tab+2);
  fprintf(fp, "<Param name=\"K\" dtype=\"INT\" value=\"%d\" />\n", head[3*n+2]);
  
  Write_Tab(fp, tab+1); 
  fprintf(fp, "</Elem>\n");
  
  // Tail Index
  Write_Tab(fp, tab+1); 
  fprintf(fp, "<Elem name=\"TailIndex\">\n");
  
  Write_Tab(fp, tab+2);
  fprintf(fp, "<Param name=\"I\" dtype=\"INT\" value=\"%d\" />\n", tail[3*n+0]);
  
  Write_Tab(fp, tab+2);
  fprintf(fp, "<Param name=\"J\" dtype=\"INT\" value=\"%d\" />\n", tail[3*n+1]);
  
  Write_Tab(fp, tab+2);
  fprintf(fp, "<Param name=\"K\" dtype=\"INT\" value=\"%d\" />\n", tail[3*n+2]);
  
  Write_Tab(fp, tab+1); 
  fprintf(fp, "</Elem>\n");
  
  Write_Tab(fp, tab); fprintf(fp, "</Elem>\n");
  
  return true;
}

/**
 * DFIファイル:出力ファイル情報要素を出力する。
 *
 * @param fp      ファイルポインタ
 * @param prefix  ファイル接頭文字
 * @param tab     インデント
 * @param step    ステップ数
 * @param mio    出力時の分割指定　 true = local / false = gather
 */
bool DFI::Write_OutFileInfo(FILE* fp, const unsigned tab, const std::string prefix, const int step, const bool mio)
{
  if (fp) {
    Write_Tab(fp, tab+1); 
    fprintf(fp, "<Elem name=\"File\" id=\"%d\" >\n", step);
  }
  
  Write_GuideCell(fp, tab+1);
  
  for(int n=0; n<Num_Node; n++) {
    if ( !Write_OutFileName(fp, tab+1, prefix, step, n, mio) ) return false;
  }

  if (fp) {
    Write_Tab(fp, tab+1); 
    fprintf(fp, "</Elem>\n");
  }
  
  return true;
}


/**
 * DFIファイル:ガイドセル要素を出力する。
 *
 * @param fp      ファイルポインタ
 * @param tab     インデント
 * @param gc     ガイドセル
 */
void DFI::Write_GuideCell(FILE* fp, const unsigned tab)
{
  Write_Tab(fp, tab+1);
  fprintf(fp, "<Param name=\"GuideCell\" dtype=\"INT\" value=\"%d\" />\n", guide);
}


/**
 * DFIファイル:ファイル名要素を出力する。
 *
 * @param fp     ファイルポインタ
 * @param tab    インデント
 * @param prefix ファイル接頭文字
 * @param step   ステップ数
 * @param id     対象ノードID
 * @param mio    出力時の分割指定　 true = local / false = gather
 */
bool DFI::Write_OutFileName(FILE* fp, const unsigned tab, const std::string prefix, const int step, const int id, const bool mio)
{
  char fname[512];
  memset(fname, 0, sizeof(char)*512);
  std::string tmp = Generate_FileName(prefix, step, id, mio);

  if( !path_util::GetFullPathName(tmp.c_str(), fname, 512) ) {
    return false;
  }

  if (fp) Write_Tab(fp, tab+1);
  if (fp) fprintf(fp, "<Param name=\"FileName\" dtype=\"STRING\" value=\"%s\" id=\"%d\" />\n", fname, id);

  return true;
}

