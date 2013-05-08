/* #################################################################
 *
 * CIOlib - Cartesian Input / Output library
 *
 * Copyright (c) 2013 Advanced Institute for Computational Science, RIKEN.
 * All right reserved.
 *
 * #################################################################
 */

/** 
 * @file   cio_DFI_Write.C
 * @brief  cio_DFI Class
 * @author kero     
 */

#include "cio_DFI.h"

// #################################################################
// Index DFIファイルの出力のコントロール
bool cio_DFI::WriteIndexDfiFile(string DfiName, int RankID, const std::string prefix, const unsigned step, REAL_TYPE time, REAL_TYPE* minmax, const bool mio)
{

  if ( prefix.empty() || RankID != 0 ) return NULL;

  if ( mio )
  {
    //dfi_name = Generate_DFI_Name(prefix);

    if( DfiName.empty() )
    {
      return false;
    }

    if( !Write_Index_File(DfiName, prefix, step, time, m_dfi_mng, minmax, mio) )
    {
      return false;
    }
  }

 // printf("WriteIndexDfiFile step : %d  m_dfi_mng : %d\n",step ,m_dfi_mng);

  return true;

}

// #################################################################
// Proc DFIファイルの出力のコントロール
bool cio_DFI::WriteProcDfiFile(int RankID)
{

  if( DFI_Fpath.Process.empty() || RankID != 0 ) return NULL;
  //printf("Process : %s\n",DFI_Fpath.Process.c_str());

  std::string dfi_name = DFI_Fpath.Process;

  if( !Write_Proc_File(dfi_name) )
  {
    return false;
  }

  return true;

}


// #################################################################
// Proc DFIファイルの出力のコントロール(static)
bool cio_DFI::WriteProcDfiFile(MPI_Comm comm, string procFileName, int G_size[3],
                               int division[3], int head[3], int tail[3], REAL_TYPE org[3],
                               REAL_TYPE pch[3], string hostname, bool out_host)
{

  if( procFileName.empty() ) return false;

  int RankID;
  MPI_Comm_rank( comm, &RankID );

  cio_MPI out_mpi;
  int nrank;
  MPI_Comm_size( comm, &nrank );
  out_mpi.NumberOfRank = nrank;
  out_mpi.NumberOfGroup = 1;

  cio_Domain out_domain;
  vector<cio_Rank> out_RankInfo;
  cio_Rank out_Rank;

  cio_Create_Domain(comm, G_size, division, head, tail, out_domain, out_RankInfo, out_Rank);
  for(int i=0; i<3; i++) {
    out_domain.GlobalOrigin[i] = org[i];
    out_domain.GlobalRegion[i] = pch[i]*G_size[i];
  }

  if( out_host ) {
  const int LEN=256;
  char *recbuf = new char[out_RankInfo.size()*LEN];
  char  sedbuf[LEN];
  sprintf(sedbuf,"%s",hostname.c_str());
  MPI_Gather(sedbuf,LEN,MPI_CHAR,recbuf,LEN,MPI_CHAR,0,MPI_COMM_WORLD);

    for( int i=0; i<out_RankInfo.size(); i++ ) {
     char* hn =&(recbuf[i*LEN]);
     out_RankInfo[i].HostName=(string(hn));
    }
  }

  if(RankID != 0) return NULL;

  if( !Write_Proc_File(procFileName,out_domain,out_mpi,out_RankInfo) )
  {
    return false;
  }

  return true;

}

// #################################################################
// 出力DFIファイル名を作成する
std::string cio_DFI::Generate_DFI_Name(const std::string prefix)
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
// DFIファイルを出力する
bool cio_DFI::Write_Index_File(const std::string dfi_name, const std::string prefix, const unsigned step, REAL_TYPE time, int& dfi_mng, REAL_TYPE *minmax, const bool mio)
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

  if ( (dfi_mng == 0) || !flag  ) // カウントゼロ=セッションの開始、または既存ファイルが存在しない
  {
    if( !(fp = fopen(dfi_name.c_str(), "w")) )
    {
      fprintf(stderr, "Can't open file.(%s)\n", dfi_name.c_str());
      return false;
    }

    if (fp) fprintf(fp, "FileInfo {\n");
    if (fp) fprintf(fp, "\n");

    if( !Write_FileInfo(fp, 0, prefix) )
    {
      if (fp) fclose(fp);
      return false;
    }

    if (fp) fprintf(fp, "\n");
    if (fp) fprintf(fp, "}\n");
    if (fp) fprintf(fp, "\n");

    if (fp) fprintf(fp, "FilePath {\n");
    if (fp) fprintf(fp, "\n");
    if (fp) Write_Process(fp, 1);
    if (fp) fprintf(fp, "\n");
    if (fp) fprintf(fp, "}\n");
    if (fp) fprintf(fp, "\n");

    if (fp) fprintf(fp, "Unit {\n");
    if (fp) fprintf(fp, "\n");

    if( !Write_Unit(fp, 0, prefix) )
    {
      if (fp) fclose(fp);
      return false;
    }
    if (fp) fprintf(fp, "\n");
    if (fp) fprintf(fp, "}\n");
    if (fp) fprintf(fp, "\n");

    if (fp) fprintf(fp, "TimeSlice {\n");
    if (fp) fprintf(fp, "\n");
    //if ( !Write_OutFileInfo(fp, 1, prefix, step, mio) )
    if ( !Write_TimeSlice(fp, 1, step, time, minmax) )
    {
      if (fp) fclose(fp);
      return false;
    }

    if (fp) fprintf(fp, "}\n\n");
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

    //if ( !Write_OutFileInfo(fp, 1, prefix, step, mio) )
    if ( !Write_TimeSlice(fp, 1, step, time, minmax) )
    {
      if (fp) fclose(fp);
      return false;
    }
    //if (fp) Write_Tab(fp, 1);
    if (fp) fprintf(fp, "}\n\n");
    if (fp) fclose(fp);

  }

  dfi_mng++;    

  return true;

}
// #################################################################
// DFIファイル:FileInfo要素を出力する
bool cio_DFI::Write_FileInfo(FILE* fp, const unsigned tab, const std::string prefix)
{
  Write_DirectoryPath(fp, tab+1, DFI_Finfo.DirectoryPath);
  Write_BaseName(fp, tab+1, prefix);
  Write_FileFormat(fp, tab+1);
  Write_GuideCell(fp, tab+1);
  Write_DataType(fp, tab+1);
  Write_Endian(fp, tab+1);
  Write_ArrayShape(fp, tab+1);
  Write_Component(fp, tab+1);
  return true;
}

// #################################################################
// DFIファイル:Unit要素を出力する
bool cio_DFI::Write_Unit(FILE* fp, const unsigned tab, const std::string prefix)
{
  Write_Length(fp, tab+1);
  Write_L0(fp, tab+1);
  Write_Velocity(fp, tab+1);
  Write_V0(fp, tab+1);
  Write_Pressure(fp, tab+1);
  Write_P0(fp, tab+1);
  Write_DiffPrs(fp, tab+1);
  Write_Temperatur(fp, tab+1);
  Write_BaseTemp(fp, tab+1);
  Write_DiffTemp(fp, tab+1);
  return true;
}

// #################################################################
// DFIファイル:TimeSlice要素を出力する
bool cio_DFI::Write_TimeSlice(FILE* fp, const unsigned tab, const unsigned step, REAL_TYPE time,
                              REAL_TYPE* minmax)
{

  REAL_TYPE comp1,comp2;
  string compname;

  comp1=0.0;
  comp2=100.0;

  Write_Tab(fp, tab);
  fprintf(fp, "Slice[@] {\n");

 // printf("Slice[@] {\n");
 // printf("step %d\n",step);

  Write_Step(fp,tab+1,step);

  Write_Time(fp,tab+1,time);

  if( DFI_Finfo.Component ) {
    Write_Tab(fp, tab+1);
    fprintf(fp, "MinMax[@] {\n");
    //MinMax
    //for(int i=0; i<DFI_Finfo.Component; i++){
    for(int i=0; i<1; i++){
      compname="Min";
      Write_Comp(fp,tab+2,compname,minmax[i*2]);
      compname="Max";
      Write_Comp(fp,tab+2,compname,minmax[i*2+1]);
    }
    Write_Tab(fp, tab+1);
    fprintf(fp, "}\n");
  }

  Write_Tab(fp, tab);
  fprintf(fp, "}\n");

  return true;
}
// #################################################################
// DFIファイル:BaseName要素を出力する
void cio_DFI::Write_DirectoryPath(FILE* fp, const unsigned tab, const std::string dirpath)
{
  Write_Tab(fp, tab);
  fprintf(fp, "DirectoryPath = \"%s\"\n", dirpath.c_str());
}

// #################################################################
// DFIファイル:BaseName要素を出力する
void cio_DFI::Write_BaseName(FILE* fp, const unsigned tab, const std::string prefix)
{
  Write_Tab(fp, tab);
  fprintf(fp, "Prefix = \"%s\"\n", prefix.c_str());
}

// #################################################################
// DFIファイル:ファイルフォーマット要素を出力する
void cio_DFI::Write_FileFormat(FILE* fp, const unsigned tab)
{
  Write_Tab(fp, tab);
  fprintf(fp, "FileFormat =\"%s\"\n",DFI_Finfo.FileFormat.c_str());
}

// #################################################################
// DFIファイル:ガイドセル要素を出力する
void cio_DFI::Write_GuideCell(FILE* fp, const unsigned tab)
{
  Write_Tab(fp, tab);
  fprintf(fp, "GuideCell = %d\n", DFI_Finfo.GuideCell);
}

// #################################################################
// DFIファイル:データタイプを出力する
void cio_DFI::Write_DataType(FILE* fp, const unsigned tab)
{
  Write_Tab(fp, tab);
  fprintf(fp, "DataType =\"%s\"\n",DFI_Finfo.DataType.c_str());
}

// #################################################################
// DFIファイル:Endianを出力する
void cio_DFI::Write_Endian(FILE* fp, const unsigned tab)
{
  Write_Tab(fp, tab);
  fprintf(fp, "Endian =\"%s\"\n",DFI_Finfo.Endian.c_str());
}

// #################################################################
// DFIファイル:ArrayShapeを出力する
void cio_DFI::Write_ArrayShape(FILE* fp, const unsigned tab)
{
  Write_Tab(fp, tab);
  fprintf(fp, "ArrayShape =\"%s\"\n",DFI_Finfo.ArrayShape.c_str());
}

// #################################################################
// DFIファイル:Componentを出力する
void cio_DFI::Write_Component(FILE* fp, const unsigned tab)
{
  Write_Tab(fp, tab);
  fprintf(fp, "Component =%d\n",DFI_Finfo.Component);
}

// #################################################################
// DFIファイル:processを出力する
void cio_DFI::Write_Process(FILE* fp, const unsigned tab)
{
  Write_Tab(fp, tab);
  fprintf(fp, "Process =\"%s\"\n",DFI_Fpath.Process.c_str());
}

// #################################################################
// DFIファイル:Lengthを出力する
void cio_DFI::Write_Length(FILE* fp, const unsigned tab)
{
  Write_Tab(fp, tab);
  fprintf(fp, "Length =\"%s\"\n",DFI_Unit.Length.c_str());
}

// #################################################################
// DFIファイル:L0を出力する
void cio_DFI::Write_L0(FILE* fp, const unsigned tab)
{
  Write_Tab(fp, tab);
  fprintf(fp, "L0 =%e\n",DFI_Unit.L0);
}

// #################################################################
// DFIファイル:Velocityhを出力する
void cio_DFI::Write_Velocity(FILE* fp, const unsigned tab)
{
  Write_Tab(fp, tab);
  fprintf(fp, "Velocity =\"%s\"\n",DFI_Unit.Velocity.c_str());
}

// #################################################################
// DFIファイル:V0を出力する
void cio_DFI::Write_V0(FILE* fp, const unsigned tab)
{
  Write_Tab(fp, tab);
  fprintf(fp, "V0 =%e\n",DFI_Unit.V0);
}

// #################################################################
// DFIファイル:Pressureを出力する
void cio_DFI::Write_Pressure(FILE* fp, const unsigned tab)
{
  Write_Tab(fp, tab);
  fprintf(fp, "Pressure =\"%s\"\n",DFI_Unit.Pressure.c_str());
}

// #################################################################
// DFIファイル:P0を出力する
void cio_DFI::Write_P0(FILE* fp, const unsigned tab)
{
  Write_Tab(fp, tab);
  fprintf(fp, "P0 =%e\n",DFI_Unit.P0);
}

// #################################################################
// DFIファイル:DiffPrsを出力する
void cio_DFI::Write_DiffPrs(FILE* fp, const unsigned tab)
{
  Write_Tab(fp, tab);
  fprintf(fp, "DiffPrs =%e\n",DFI_Unit.DiffPrs);
}

// #################################################################
// DFIファイル:Temperaturを出力する
void cio_DFI::Write_Temperatur(FILE* fp, const unsigned tab)
{
  Write_Tab(fp, tab);
  fprintf(fp, "Temperatur =\"%s\"\n",DFI_Unit.Temperatur.c_str());
}

// #################################################################
// DFIファイル:BaseTempを出力する
void cio_DFI::Write_BaseTemp(FILE* fp, const unsigned tab)
{
  Write_Tab(fp, tab);
  fprintf(fp, "BaseTemp =%e\n",DFI_Unit.BaseTemp);
}

// #################################################################
// DFIファイル:DiffTempを出力する
void cio_DFI::Write_DiffTemp(FILE* fp, const unsigned tab)
{
  Write_Tab(fp, tab);
  fprintf(fp, "DiffTemp =%e\n",DFI_Unit.DiffTemp);
}

// #################################################################
// DFIファイル:Stepを出力する
void cio_DFI::Write_Step(FILE* fp, const unsigned tab, int step)
{
  Write_Tab(fp, tab);
  fprintf(fp, "Step =%d\n",step);
}

// #################################################################
// DFIファイル:Timeを出力する
void cio_DFI::Write_Time(FILE* fp, const unsigned tab, REAL_TYPE time)
{
  Write_Tab(fp, tab);
  fprintf(fp, "Time =%e\n",time);
}

// #################################################################
// DFIファイル:MinMax等を出力する
void cio_DFI::Write_Comp(FILE* fp, const unsigned tab, const std::string compname, REAL_TYPE comp)
{
  Write_Tab(fp, tab);
  fprintf(fp, "%s =%e\n",compname.c_str(),comp);
}

// #################################################################
// DFIファイルを出力する
bool cio_DFI::Write_Proc_File(const std::string dfi_name)
{

  FILE* fp = NULL;
  if( !(fp = fopen(dfi_name.c_str(), "w")) )
  {
    fprintf(stderr, "Can't open file.(%s)\n", dfi_name.c_str());
    return false;
  }

  if (fp) fprintf(fp, "Domain {\n");
  if (fp) fprintf(fp, "\n");
  if( !Write_Domain(fp, 0) )
  {
    if (fp) fclose(fp);
    return false;
  }
  if (fp) fprintf(fp, "\n");
  if (fp) fprintf(fp, "}\n");
  if (fp) fprintf(fp, "\n");

  if (fp) fprintf(fp, "MPI {\n");
  if (fp) fprintf(fp, "\n");
  if( !Write_MPI(fp, 0) )
  {
    if (fp) fclose(fp);
    return false;
  }
  if (fp) fprintf(fp, "\n");
  if (fp) fprintf(fp, "}\n");
  if (fp) fprintf(fp, "\n");

  if (fp) fprintf(fp, "Process {\n");
  if (fp) fprintf(fp, "\n");
  if( !Write_Process_Rank(fp, 0) )
  {
    if (fp) fclose(fp);
    return false;
  }
  if (fp) fprintf(fp, "\n");
  if (fp) fprintf(fp, "}\n");

  if (fp) fclose(fp);
  return true;
}


// #################################################################
// DFIファイルを出力する(static)
bool cio_DFI::Write_Proc_File(const std::string dfi_name, cio_Domain out_domain,
                              cio_MPI out_mpi, vector<cio_Rank> out_RankInfo)
{

  FILE* fp = NULL;
  if( !(fp = fopen(dfi_name.c_str(), "w")) )
  {
    fprintf(stderr, "Can't open file.(%s)\n", dfi_name.c_str());
    return false;
  }
  if (fp) fprintf(fp, "Domain {\n");
  if (fp) fprintf(fp, "\n");
  if( !Write_Domain(fp, 0, out_domain) )
  {
    if (fp) fclose(fp);
    return false;
  }
  if (fp) fprintf(fp, "\n");
  if (fp) fprintf(fp, "}\n");
  if (fp) fprintf(fp, "\n");

  if (fp) fprintf(fp, "MPI {\n");
  if (fp) fprintf(fp, "\n");
  if( !Write_MPI(fp, 0, out_mpi) )
  {
    if (fp) fclose(fp);
    return false;
  }
  if (fp) fprintf(fp, "\n");
  if (fp) fprintf(fp, "}\n");
  if (fp) fprintf(fp, "\n");

  if (fp) fprintf(fp, "Process {\n");
  if (fp) fprintf(fp, "\n");
  if( !Write_Process_Rank(fp, 0, out_RankInfo) )
  {
    if (fp) fclose(fp);
    return false;
  }
  if (fp) fprintf(fp, "\n");
  if (fp) fprintf(fp, "}\n");

  if (fp) fclose(fp);
  return true;
}

// #################################################################
// DFIファイル:Domain要素を出力する
bool cio_DFI::Write_Domain(FILE* fp, const unsigned tab)
{
  Write_Origin(fp, tab+1);
  Write_Region(fp, tab+1);
  Write_WholeSize(fp, tab+1);
  Write_NumDivDomain(fp, tab+1);
  return true;
}

// #################################################################
// DFIファイル:Domain要素を出力する (static)
bool cio_DFI::Write_Domain(FILE* fp, const unsigned tab, cio_Domain out_domain)
{
  Write_Origin(fp, tab+1, out_domain.GlobalOrigin);
  Write_Region(fp, tab+1, out_domain.GlobalRegion);
  Write_WholeSize(fp, tab+1, out_domain.GlobalVoxel);
  Write_NumDivDomain(fp, tab+1, out_domain.GlobalDivision);
  return true;
}

// #################################################################
// DFIファイル:MPI要素を出力する
bool cio_DFI::Write_MPI(FILE* fp, const unsigned tab)
{
  Write_NodeNum(fp, tab+1);
  return true;
}

// #################################################################
// DFIファイル:MPI要素を出力する(static)
bool cio_DFI::Write_MPI(FILE* fp, const unsigned tab, cio_MPI out_mpi)
{
  Write_NodeNum(fp, tab+1, out_mpi.NumberOfRank);
  return true;
}


// #################################################################
// DFIファイル:Process要素を出力する
bool cio_DFI::Write_Process_Rank(FILE* fp, const unsigned tab)
{

  for(int i=0; i<RankInfo.size(); i++) {
    if( !Write_Rank(fp, tab+1, i) ) return false;
  }

  return true;
}

// #################################################################
// DFIファイル:Process要素を出力する(static)
bool cio_DFI::Write_Process_Rank(FILE* fp, const unsigned tab, vector<cio_Rank> out_RankInfo)
{

  for(int i=0; i<out_RankInfo.size(); i++) {
    if( !Write_Rank(fp, tab+1,out_RankInfo[i] ) ) return false;
  }

  return true;
}



// #################################################################
// DFIファイル:Origin要素を出力する
void cio_DFI::Write_Origin(FILE* fp, const unsigned tab)
{
  Write_Tab(fp, tab);
  fprintf(fp, "GlobalOrigin   = (%e, %e, %e)\n", DFI_Domain.GlobalOrigin[0], 
                                                 DFI_Domain.GlobalOrigin[1],
                                                 DFI_Domain.GlobalOrigin[2]);
}

// #################################################################
// DFIファイル:Origin要素を出力する (static)
void cio_DFI::Write_Origin(FILE* fp, const unsigned tab, REAL_TYPE org[3])
{
  //Write_Tab(fp, tab);
  for(int n=0; n<tab; n++) fprintf(fp, "  ");
  fprintf(fp, "GlobalOrigin   = (%e, %e, %e)\n", org[0], 
                                                 org[1],
                                                 org[2]);
}

// #################################################################
// DFIファイル:Resion要素を出力する
void cio_DFI::Write_Region(FILE* fp, const unsigned tab)
{
  Write_Tab(fp, tab);
  fprintf(fp, "GlobalRegion   = (%e, %e, %e)\n",DFI_Domain.GlobalRegion[0],
                                                DFI_Domain.GlobalRegion[1],
                                                DFI_Domain.GlobalRegion[2]);
}

// #################################################################
// DFIファイル:Resion要素を出力する (static)
void cio_DFI::Write_Region(FILE* fp, const unsigned tab, REAL_TYPE Region[3])
{
  //Write_Tab(fp, tab);
  for(int n=0; n<tab; n++) fprintf(fp, "  ");
  fprintf(fp, "GlobalRegion   = (%e, %e, %e)\n",Region[0],
                                                Region[1],
                                                Region[2]);
}

// #################################################################
// DFIファイル:全体ボクセルサイズ要素を出力する
void cio_DFI::Write_WholeSize(FILE* fp, const unsigned tab)
{
  Write_Tab(fp, tab);
  fprintf(fp, "GlobalVoxel    = (%d, %d, %d)\n", DFI_Domain.GlobalVoxel[0],
                                                 DFI_Domain.GlobalVoxel[1],
                                                 DFI_Domain.GlobalVoxel[2]);
}

// #################################################################
// DFIファイル:全体ボクセルサイズ要素を出力する(static)
void cio_DFI::Write_WholeSize(FILE* fp, const unsigned tab, int GlobalVoxel[3])
{
  //Write_Tab(fp, tab);
  for(int n=0; n<tab; n++) fprintf(fp, "  ");
  fprintf(fp, "GlobalVoxel    = (%d, %d, %d)\n", GlobalVoxel[0],
                                                 GlobalVoxel[1],
                                                 GlobalVoxel[2]);
}

// #################################################################
// DFIファイル:I,J,K分割数要素を出力する
void cio_DFI::Write_NumDivDomain(FILE* fp, const unsigned tab)
{
  Write_Tab(fp, tab);
  fprintf(fp, "GlobalDivision = (%d, %d, %d)\n", DFI_Domain.GlobalDivision[0],
                                                 DFI_Domain.GlobalDivision[1],
                                                 DFI_Domain.GlobalDivision[2]);
}

// #################################################################
// DFIファイル:I,J,K分割数要素を出力する(static)
void cio_DFI::Write_NumDivDomain(FILE* fp, const unsigned tab, int GlobalDivision[3])
{
  //Write_Tab(fp, tab);
  for(int n=0; n<tab; n++) fprintf(fp, "  ");
  fprintf(fp, "GlobalDivision = (%d, %d, %d)\n", GlobalDivision[0],
                                                 GlobalDivision[1],
                                                 GlobalDivision[2]);
}

// #################################################################
// DFIファイル:ノード数要素を出力する
void cio_DFI::Write_NodeNum(FILE* fp, const unsigned tab)
{
  Write_Tab(fp, tab);
  fprintf(fp, "NumberOfRank   = %d\n", DFI_MPI.NumberOfRank);

  Write_Tab(fp, tab);
  fprintf(fp, "NumberOfGroup  = %d\n", 1);
}

// #################################################################
// DFIファイル:ノード数要素を出力する
void cio_DFI::Write_NodeNum(FILE* fp, const unsigned tab, int NumberOfRank)
{
  //Write_Tab(fp, tab);
  for(int n=0; n<tab; n++) fprintf(fp, "  ");
  fprintf(fp, "NumberOfRank   = %d\n", NumberOfRank);

  //Write_Tab(fp, tab);
  for(int n=0; n<tab; n++) fprintf(fp, "  ");
  fprintf(fp, "NumberOfGroup  = %d\n", 1);
}

// #################################################################
// DFIファイル:Rank情報要素を出力する
bool cio_DFI::Write_Rank(FILE* fp, const unsigned tab, const int n)
{
  Write_Tab(fp, tab);
  fprintf(fp, "Rank[@] {\n");
  fprintf(fp, "\n");

  Write_ID(fp,tab+1,n);
  Write_Hostname(fp,tab+1,n);
  Write_L_VoxelSize(fp,tab+1,n);
  Write_HeadIndex(fp,tab+1,n);
  Write_TailIndex(fp,tab+1,n);

  fprintf(fp, "\n");
  Write_Tab(fp, tab);
  fprintf(fp, "}\n");

  return true;
}

// #################################################################
// DFIファイル:Rank情報要素を出力する (static)
bool cio_DFI::Write_Rank(FILE* fp, const unsigned tab, cio_Rank rank)
{
  //Write_Tab(fp, tab);
  for(int n=0; n<tab; n++) fprintf(fp, "  ");
  fprintf(fp, "Rank[@] {\n");
  fprintf(fp, "\n");

  Write_RankID(fp,tab+1,rank.RankID);
  Write_Hostname(fp,tab+1,rank.HostName);
  Write_L_VoxelSize(fp,tab+1,rank.VoxelSize);
  Write_HeadIndex(fp,tab+1,rank.HeadIndex);
  Write_TailIndex(fp,tab+1,rank.TailIndex);

  fprintf(fp, "\n");
  //Write_Tab(fp, tab);
  for(int n=0; n<tab; n++) fprintf(fp, "  ");
  fprintf(fp, "}\n");

  return true;
}


// #################################################################
// DFIファイル:ノード番号要素を出力する
void cio_DFI::Write_ID(FILE* fp, const unsigned tab, const int n)
{
  Write_Tab(fp, tab);
  fprintf(fp, "ID        = %d\n", n);
}

// #################################################################
// DFIファイル:ノード番号要素を出力する
void cio_DFI::Write_RankID(FILE* fp, const unsigned tab, const int n)
{
  //Write_Tab(fp, tab);
  for(int n=0; n<tab; n++) fprintf(fp, "  ");
  fprintf(fp, "ID        = %d\n", n);
}

// #################################################################
// DFIファイル:ノード番号要素を出力する
void cio_DFI::Write_Hostname(FILE* fp, const unsigned tab, const int n)
{
  Write_Tab(fp, tab);
  fprintf(fp, "HostName  = \"%s\"\n", RankInfo[n].HostName.c_str());
}

// #################################################################
// DFIファイル:ノード番号要素を出力する(static)
void cio_DFI::Write_Hostname(FILE* fp, const unsigned tab, string HostName)
{
  //Write_Tab(fp, tab);
  for(int n=0; n<tab; n++) fprintf(fp, "  ");
  fprintf(fp, "HostName  = \"%s\"\n", HostName.c_str());
}

// #################################################################
// DFIファイル:ノードのボクセルサイズ要素を出力する
void cio_DFI::Write_L_VoxelSize(FILE* fp, const unsigned tab, const int n)
{
  Write_Tab(fp, tab);
  fprintf(fp, "VoxelSize = (%d, %d, %d)\n", RankInfo[n].VoxelSize[0],
                                            RankInfo[n].VoxelSize[1],
                                            RankInfo[n].VoxelSize[2]);
}

// #################################################################
// DFIファイル:ノードのボクセルサイズ要素を出力する(static)
void cio_DFI::Write_L_VoxelSize(FILE* fp, const unsigned tab, int VoxelSize[3])
{
  //Write_Tab(fp, tab);
  for(int n=0; n<tab; n++) fprintf(fp, "  ");
  fprintf(fp, "VoxelSize = (%d, %d, %d)\n", VoxelSize[0],
                                            VoxelSize[1],
                                            VoxelSize[2]);
}


// #################################################################
// DFIファイル:ノードのheadIndex要素を出力する
void cio_DFI::Write_HeadIndex(FILE* fp, const unsigned tab, const int n)
{
  Write_Tab(fp, tab);
  fprintf(fp, "HeadIndex = (%d, %d, %d)\n", RankInfo[n].HeadIndex[0],
                                            RankInfo[n].HeadIndex[1],
                                            RankInfo[n].HeadIndex[2]);
}

// #################################################################
// DFIファイル:ノードのheadIndex要素を出力する
void cio_DFI::Write_HeadIndex(FILE* fp, const unsigned tab, int HeadIndex[3])
{
  //Write_Tab(fp, tab);
  for(int n=0; n<tab; n++) fprintf(fp, "  ");
  fprintf(fp, "HeadIndex = (%d, %d, %d)\n", HeadIndex[0],
                                            HeadIndex[1],
                                            HeadIndex[2]);
}

// #################################################################
// DFIファイル:ノードのTailIndex要素を出力する
void cio_DFI::Write_TailIndex(FILE* fp, const unsigned tab, const int n)
{
  Write_Tab(fp, tab);
  fprintf(fp, "TailIndex = (%d, %d, %d)\n", RankInfo[n].TailIndex[0],
                                            RankInfo[n].TailIndex[1],
                                            RankInfo[n].TailIndex[2]);
}

// #################################################################
// DFIファイル:ノードのTailIndex要素を出力する
void cio_DFI::Write_TailIndex(FILE* fp, const unsigned tab, int TailIndex[3])
{
  //Write_Tab(fp, tab);
  for(int n=0; n<tab; n++) fprintf(fp, "  ");
  fprintf(fp, "TailIndex = (%d, %d, %d)\n", TailIndex[0],
                                            TailIndex[1],
                                            TailIndex[2]);
}

// #################################################################
// DFIファイル:ノード番号要素を出力する
void cio_DFI::Write_MyID(FILE* fp, const unsigned tab, const int n)
{
  Write_Tab(fp, tab);
  fprintf(fp, "ID = %d\n", RankInfo[n].RankID);
}

// #################################################################
// Tab(space２つ)を出力する
void cio_DFI::Write_Tab(FILE* fp, const unsigned tab)
{
  for(int n=0; n<tab; n++) fprintf(fp, "  ");
}

// #################################################################
// DFIファイル:出力ファイル情報要素を出力する
bool cio_DFI::Write_OutFileInfo(FILE* fp, const unsigned tab, const std::string prefix, const unsigned step, const bool mio)
{
  if (fp)
  {
    Write_Tab(fp, tab+1);
    fprintf(fp, "Step[@] = %d\n", step);
  }

  return true;
}
