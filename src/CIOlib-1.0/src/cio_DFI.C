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
 * @file   cio_DFI.C
 * @brief  cio_DFI Class
 * @author kero    
 */

#include <unistd.h> // for gethostname() of FX10/K

#include "cio_DFI.h"
#include "cio_DFI_SPH.h"
#include "cio_DFI_BOV.h"

// #################################################################
// コンストラクタ
cio_DFI::cio_DFI()
{

 m_dfi_mng = 0;

}


// #################################################################
// デストラクタ
cio_DFI::~cio_DFI()
{

}

// #################################################################
//
cio_DFI* cio_DFI::ReadInit(MPI_Comm comm, string DfiName)
{

  int RankID=0;
  cio_TextParser tpCntl;

  //index.dfi read
  //TPインスタンス
  tpCntl.getTPinstance();

  FILE*fp = NULL;
  if( !(fp=fopen(DfiName.c_str(),"rb")) ) {
    printf("Can't open file. (%s)\n",DfiName.c_str());
    return false;
  }
  fclose(fp);

  //入力ファイル index.dfi をセット
  int ierror = 0;
  ierror = tpCntl.readTPfile(DfiName);
  if ( ierror )
  {
    printf("\tinput file not found '%s'\n",DfiName.c_str());
    return NULL;
  }

  cio_FileInfo F_info;
  if( readFileInfo(DfiName,tpCntl,F_info) != CIO_SUCCESS ) 
  {
    printf("\tFileInfo Data Read error %s\n",DfiName.c_str());
    return NULL;
  }

  cio_FilePath F_path;
  if( readFilePath(DfiName,tpCntl,F_path) != CIO_SUCCESS )
  {
    printf("\tFilePath Data Read error %s\n",DfiName.c_str());
    return NULL;
  }

  cio_Unit unit;
  if( readUnit(DfiName,tpCntl,unit) != CIO_SUCCESS )
  {
    printf("\tUnit Data Read error %s\n",DfiName.c_str());
    return NULL;
  }

  vector<cio_Slice> TimeSlice;
  cio_Slice slice;
  if( readSlice(DfiName,tpCntl,TimeSlice,slice) != CIO_SUCCESS )
  {
    printf("\tTimeSlice Data Read error %s\n",DfiName.c_str());
    return NULL;
  }

  //TextParserの破棄
  tpCntl.remove();

  //proc file name set
  string procfile = F_path.Process;

  //proc.dfi read
  //TPインスタンス
  tpCntl.getTPinstance();

  fp = NULL;
  if( !(fp=fopen(procfile.c_str(),"rb")) ) {
    printf("Can't open file. (%s)\n",procfile.c_str());
    return false;
  }
  fclose(fp);

  //入力ファイル proc.dfi をセット
  ierror = tpCntl.readTPfile(procfile);
  if ( ierror )
  {
    printf("\tinput file not found '%s'\n",procfile.c_str());
    return NULL;
  }

  cio_Domain domain;
  if( readDomain(procfile,tpCntl,domain) != CIO_SUCCESS ) 
  {
    printf("\tDomain Data Read error %s\n",procfile.c_str());
    return NULL;
  }
    

  cio_MPI mpi;
  if( readMPI(procfile,tpCntl,mpi) != CIO_SUCCESS )
  {
    printf("\tMPI Data Read error %s\n",procfile.c_str());
    return NULL;
  }

  vector<cio_Rank> RankInfo;
  cio_Rank rank;
  if( readRank(procfile,tpCntl,RankInfo,rank) != CIO_SUCCESS )
  {
    printf("\tProcess Data Read error %s\n",procfile.c_str());
    return NULL;
  }

  //TextParserの破棄
  tpCntl.remove();

  std::string fmt = F_info.FileFormat;
  //printf("FileFormat : %s\n",fmt.c_str());

  cio_DFI *dfi = NULL;
  //if( fmt == "sph" ) {
  if( !strcasecmp(fmt.c_str() , "sph" ) ) {
    dfi = new cio_DFI_SPH(F_info, F_path, unit, domain, mpi, TimeSlice, RankInfo);
    dfi->m_comm = comm;
    dfi->m_indexDfiName = DfiName;
  } else if( !strcasecmp(fmt.c_str() , "bov" ) ) {
    dfi = new cio_DFI_BOV();
  }

  return dfi;

}

// #################################################################
// 
cio_DFI* cio_DFI::WriteInit(MPI_Comm comm, string DfiName, string Path, string prefix, 
                            string format, int GCell, string DataType, string ArrayShape,
                            int Comp, string process, int G_size[3], REAL_TYPE pitch[3],
                            REAL_TYPE G_origin[3], int division[3], int head[3], int tail[3],
                            string hostname)
{


  cio_DFI *dfi = NULL;

  int RankID;
  MPI_Comm_rank( comm, &RankID );

  cio_FileInfo out_F_info;
  out_F_info.DirectoryPath = Path;
  out_F_info.Prefix        = prefix;
  out_F_info.FileFormat    = format;
  out_F_info.GuideCell     = GCell;
  out_F_info.DataType      = DataType;
  out_F_info.ArrayShape    = ArrayShape;
  out_F_info.Component     = Comp;
 
  int idumy = 1;
  char* cdumy = (char*)(&idumy);
  if( cdumy[0] == 0x01 ) out_F_info.Endian = "little";
  if( cdumy[0] == 0x00 ) out_F_info.Endian = "big";

  cio_FilePath out_F_path;
  out_F_path.Process = process;

  cio_Unit out_unit;

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
    out_domain.GlobalOrigin[i] = G_origin[i];
    out_domain.GlobalRegion[i] = pitch[i]*G_size[i];
  }

  vector<cio_Slice> out_TSlice;

  //if( format == "sph" ) { 
  if( !strcasecmp(format.c_str() , "sph" ) ) {
    //printf("out data format sph\n");
    dfi = new cio_DFI_SPH(out_F_info, out_F_path, out_unit, out_domain, out_mpi, 
                          out_TSlice, out_RankInfo);
  } else if( !strcasecmp(format.c_str(), "bov" ) ) {
    dfi = new cio_DFI_BOV();
  }


  dfi->m_comm = comm;
  dfi->m_indexDfiName = DfiName;
  //REAL_TYPE minmax[2];
  //minmax[0]=0.0:
  //minmax[1]=0.0;
  //dfi->WriteIndexDfiFile(DfiName,RankID,prefix,0,0.0,minmax,true);
  //dfi->WriteProcDfiFile(RankID);


  char tmpname[512];
  memset(tmpname,0x00,sizeof(char)*512);
  if( gethostname(tmpname, 512) != 0 ) printf("*** error gethostname() \n");
  //printf("tmpname : %s\n",tmpname);

  

  return dfi;

}

// #################################################################
// 初期化
void cio_DFI::InitDFI()
{
}
//

// #################################################################
//
int cio_DFI::readFileInfo(string dfifile, cio_TextParser tpCntl, cio_FileInfo &finfo)
{

  string str;
  string label,label_base,label_leaf;
  int ct;
  REAL_TYPE dt;

  //Directorypath
  label = "/FileInfo/DirectoryPath";
  if ( !(tpCntl.GetValue(label, &str )) )
  {
    printf("\tParsing error : fail to get '%s'\n",label.c_str());
    return CIO_ERROR;
  }
  finfo.DirectoryPath=str;

  //Prefix
  label = "/FileInfo/Prefix";
  if ( !(tpCntl.GetValue(label, &str )) )
  {
    printf("\tParsing error : fail to get '%s'\n",label.c_str());
    return CIO_ERROR;
  }
  finfo.Prefix=str;

  //FileFormat
  label = "/FileInfo/FileFormat";
  if ( !(tpCntl.GetValue(label, &str )) )
  {
    printf("\tParsing error : fail to get '%s'\n",label.c_str());
    return CIO_ERROR;
  }
  finfo.FileFormat=str;

  //GuidCell
  label = "/FileInfo/GuideCell";
  if ( !(tpCntl.GetValue(label, &ct )) )
  {
    printf("\tParsing error : fail to get '%s'\n",label.c_str());
    return CIO_ERROR;
  }
  finfo.GuideCell=ct;

  //DataType
  label = "/FileInfo/DataType";
  if ( !(tpCntl.GetValue(label, &str )) )
  {
    printf("\tParsing error : fail to get '%s'\n",label.c_str());
    return CIO_ERROR;
  }
  finfo.DataType=str;

  //Endian
  label = "/FileInfo/Endian";
  if ( !(tpCntl.GetValue(label, &str )) )
  {
    printf("\tParsing error : fail to get '%s'\n",label.c_str());
    return CIO_ERROR;
  }
  finfo.Endian=str;

  //ArrayShape  
  label = "/FileInfo/ArrayShape";
  if ( !(tpCntl.GetValue(label, &str )) )
  {
    printf("\tParsing error : fail to get '%s'\n",label.c_str());
    return CIO_ERROR;
  }
  finfo.ArrayShape=str;

  //Componet  
  label = "/FileInfo/Component";
  if ( !(tpCntl.GetValue(label, &ct )) )
  {
    printf("\tParsing error : fail to get '%s'\n",label.c_str());
    return CIO_ERROR;
  }
  finfo.Component=ct;

  return CIO_SUCCESS;
}

// #################################################################
//
int cio_DFI::readFilePath(string dfifile, cio_TextParser tpCntl, cio_FilePath &fpath)
{

  string str;
  string label;

  //Process
  label = "/FilePath/Process";
  if ( !(tpCntl.GetValue(label, &str )) )
  {
    printf("\tParsing error : fail to get '%s'\n",label.c_str());
    return CIO_ERROR;
  }
  fpath.Process=str;

  return CIO_SUCCESS;

}
// #################################################################
//
int cio_DFI::readUnit(string dfifile, cio_TextParser tpCntl, cio_Unit &unit)
{

  string str;
  string label;
  int ct;
  REAL_TYPE dt;

  //Length
  label = "/Unit/Length";
  if ( !(tpCntl.GetValue(label, &str )) )
  {
    printf("\tParsing error : fail to get '%s'\n",label.c_str());
    return CIO_ERROR; 
  }
  unit.Length=str;

  //L0
  label = "/Unit/L0";
  if ( !(tpCntl.GetValue(label, &dt )) )
  {
    printf("\tParsing error : fail to get '%s'\n",label.c_str());
    return CIO_ERROR; 
  }
  unit.L0=dt;

  //Velocity
  label = "/Unit/Velocity";
  if ( !(tpCntl.GetValue(label, &str )) )
  {
    printf("\tParsing error : fail to get '%s'\n",label.c_str());
    return CIO_ERROR; 
  }
  unit.Velocity=str;

  //L0
  label = "/Unit/V0";
  if ( !(tpCntl.GetValue(label, &dt )) )
  {
    printf("\tParsing error : fail to get '%s'\n",label.c_str());
    return CIO_ERROR; 
  }
  unit.V0=dt;

  //Pressure
  label = "/Unit/Pressure";
  if ( !(tpCntl.GetValue(label, &str )) )
  {
    printf("\tParsing error : fail to get '%s'\n",label.c_str());
    return CIO_ERROR; 
  }
  unit.Pressure=str;

  //P0
  label = "/Unit/P0";
  if ( !(tpCntl.GetValue(label, &dt )) )
  {
    printf("\tParsing error : fail to get '%s'\n",label.c_str());
    return CIO_ERROR; 
  }
  unit.P0=dt;

  //DiffPrs
  label = "/Unit/DiffPrs";
  if ( !(tpCntl.GetValue(label, &dt )) )
  {
    printf("\tParsing error : fail to get '%s'\n",label.c_str());
    return CIO_ERROR; 
  }
  unit.DiffPrs=dt;

  //Temperatur
  label = "/Unit/Temperatur";
  if ( !(tpCntl.GetValue(label, &str )) )
  {
    //printf("\tParsing error : fail to get '%s'\n",label.c_str());
    //return CIO_ERROR; 
    unit.Temperatur="";
  } else {
    unit.Temperatur=str;
  }

  //BaseTemp
  label = "/Unit/BaseTemp";
  if ( !(tpCntl.GetValue(label, &dt )) )
  {
    //printf("\tParsing error : fail to get '%s'\n",label.c_str());
    //return CIO_ERROR; 
    unit.BaseTemp=0.0;
  } else {
    unit.BaseTemp=dt;
  }

  //DiffTemp
  label = "/Unit/DiffTemp";
  if ( !(tpCntl.GetValue(label, &dt )) )
  {
    //printf("\tParsing error : fail to get '%s'\n",label.c_str());
    //return CIO_ERROR; 
    unit.DiffTemp = 0.0;
  } else {
    unit.DiffTemp=dt;
  }

  return CIO_SUCCESS; 

}

// #################################################################
//
int cio_DFI::readSlice(string dfifile, cio_TextParser tpCntl, vector<cio_Slice> &TimeSlice, cio_Slice  slice)
{

  string str;
  string label,label_base,label_leaf,label_leaf_leaf;
  int ct;
  REAL_TYPE dt;
  int nnode=0;

  //TimeSlice
  nnode=0;
  label_base = "/TimeSlice";
  if ( tpCntl.chkNode(label_base) )  //があれば
  {
    nnode = tpCntl.countLabels(label_base);
  }


  for (int i=0; i<nnode; i++) {

    if(!tpCntl.GetNodeStr(label_base,i+1,&str))
    {
      printf("\tParsing error : No Elem name\n");
      return CIO_ERROR;
    }
    if( strcasecmp(str.substr(0,5).c_str(), "Slice") ) continue;
    label_leaf=label_base+"/"+str;

    //Step
    label = label_leaf + "/Step";
    if ( !(tpCntl.GetValue(label, &ct )) ) {
      printf("\tParsing error : fail to get '%s'\n",label.c_str());
      return CIO_ERROR;
    }
    else {
      slice.step=ct;
    }

    //Time
    label = label_leaf + "/Time";
    if ( !(tpCntl.GetValue(label, &dt )) ) {
      printf("\tParsing error : fail to get '%s'\n",label.c_str());
      return CIO_ERROR;
    }
    else {
      slice.time= dt;
    }

    //MinMax
    int ncomp=0;
    label_leaf_leaf = label_leaf + "/MinMax";
    if ( tpCntl.chkNode(label_leaf_leaf) )  //があれば
    {
      ncomp = tpCntl.countLabels(label_leaf_leaf);
    }
   
 
    slice.Min.clear();
    slice.Max.clear();

    for ( int j=0; j<ncomp; j++ ) {

      if(!tpCntl.GetNodeStr(label_leaf,j+3,&str))
      {
        printf("\tParsing error : No Elem name\n");
        return CIO_ERROR;
      } 
      if( strcasecmp(str.substr(0,6).c_str(), "minmax") ) continue;
      label_leaf_leaf = label_leaf+"/"+str;

      label = label_leaf_leaf + "/Min";
      if ( !(tpCntl.GetValue(label, &dt )) ) {
        printf("\tParsing error : fail to get '%s'\n",label.c_str());
        return CIO_ERROR;
      }
      else {
        slice.Min.push_back(dt);
      }

      label = label_leaf_leaf + "/Max";
      if ( !(tpCntl.GetValue(label, &dt )) ) {
        printf("\tParsing error : fail to get '%s'\n",label.c_str());
        return CIO_ERROR;
      }
      else {
        slice.Max.push_back(dt);
      }

    }

   TimeSlice.push_back(slice); 

  }

  return CIO_SUCCESS;

}

// #################################################################
//
int cio_DFI::readDomain(string dfifile, cio_TextParser tpCntl, cio_Domain &domain)
{

  string str;
  string label;
  int ct;
  REAL_TYPE dt;
  REAL_TYPE v[3];

  //GlobalOrign
  label = "/Domain/GlobalOrigin";
  for (int n=0; n<3; n++) v[n]=0.0;
  if ( !(tpCntl.GetVector(label, v, 3 )) ) 
  {
    printf("\tParsing error : fail to get '%s'\n",label.c_str());
    return CIO_ERROR;
  }
  domain.GlobalOrigin[0]=v[0];
  domain.GlobalOrigin[1]=v[1];
  domain.GlobalOrigin[2]=v[2];

  //GlobalRegion
  label = "/Domain/GlobalRegion";
  for (int n=0; n<3; n++) v[n]=0.0;
  if ( !(tpCntl.GetVector(label, v, 3 )) ) 
  {
    printf("\tParsing error : fail to get '%s'\n",label.c_str());
    return CIO_ERROR;
  }
  domain.GlobalRegion[0]=v[0];
  domain.GlobalRegion[1]=v[1];
  domain.GlobalRegion[2]=v[2];

  //Global_Voxel
  label = "/Domain/GlobalVoxel";
  for (int n=0; n<3; n++) v[n]=0.0;
  if ( !(tpCntl.GetVector(label, v, 3 )) ) 
  {
    printf("\tParsing error : fail to get '%s'\n",label.c_str());
    return CIO_ERROR;
  }
  domain.GlobalVoxel[0]=v[0];
  domain.GlobalVoxel[1]=v[1];
  domain.GlobalVoxel[2]=v[2];

  //Global_Division
  label = "/Domain/GlobalDivision";
  for (int n=0; n<3; n++) v[n]=0.0;
  if ( !(tpCntl.GetVector(label, v, 3 )) ) 
  {
    printf("\tParsing error : fail to get '%s'\n",label.c_str());
    return CIO_ERROR;
  }
  domain.GlobalDivision[0]=v[0];
  domain.GlobalDivision[1]=v[1];
  domain.GlobalDivision[2]=v[2];

  return CIO_SUCCESS;

}

// #################################################################
//
int cio_DFI::readMPI(string dfifile, cio_TextParser tpCntl, cio_MPI &mpi)
{

  string str;
  string label;
  int ct;

  //NumberOfRank
  label = "/MPI/NumberOfRank";
  if ( !(tpCntl.GetValue(label, &ct )) ) {
    printf("\tParsing error : fail to get '%s'\n",label.c_str());
    return CIO_ERROR;
  }
  else {
    mpi.NumberOfRank = ct;
  }

  //NumberOfGroup
  label = "/MPI/NumberOfGroup";
  if ( !(tpCntl.GetValue(label, &ct )) ) {
    printf("\tParsing error : fail to get '%s'\n",label.c_str());
    return CIO_ERROR;
  }
  else {
    mpi.NumberOfGroup = ct;
  }

  return CIO_SUCCESS;

}

// #################################################################
//
int cio_DFI::readRank(string dfifile, cio_TextParser tpCntl, vector<cio_Rank> &RankInfo, cio_Rank rank)
{

  string str;
  string label,label_base,label_leaf;
  int ct;
  REAL_TYPE dt;
  REAL_TYPE v[3];
  int nnode=0;

  //Process 
  nnode=0;
  label_base = "/Process";
  if ( tpCntl.chkNode(label_base) )  //nodeがあれば
  {
    nnode = tpCntl.countLabels(label_base);
  }

  for (int i=0; i<nnode; i++) {

    if(!tpCntl.GetNodeStr(label_base,i+1,&str))
    {
      printf("\tParsing error : No Elem name\n");
      return CIO_ERROR;
    }
    if( strcasecmp(str.substr(0,4).c_str(), "Rank") ) continue;
    label_leaf=label_base+"/"+str;

    //ID
    label = label_leaf + "/ID";
    if ( !(tpCntl.GetValue(label, &ct )) ) {
      printf("\tParsing error : fail to get '%s'\n",label.c_str());
      return CIO_ERROR;
    }
    else {
      rank.RankID= ct;
    }

    //HostName
    label = label_leaf + "/HostName";
    if ( !(tpCntl.GetValue(label, &str )) ) {
      printf("\tParsing error : fail to get '%s'\n", label.c_str());
      return CIO_ERROR;
    }
    rank.HostName= str;

    //VoxelSize
    label = label_leaf + "/VoxelSize";
    for (int n=0; n<3; n++) v[n]=0.0;
    if ( !(tpCntl.GetVector(label, v, 3 )) ) 
    {
      printf("\tParsing error : fail to get '%s'\n",label.c_str());
      return CIO_ERROR;
    }
    rank.VoxelSize[0]=v[0];
    rank.VoxelSize[1]=v[1];
    rank.VoxelSize[2]=v[2];

    //HeadIndex
    label = label_leaf + "/HeadIndex";
    for (int n=0; n<3; n++) v[n]=0.0;
    if ( !(tpCntl.GetVector(label, v, 3 )) ) 
    {
      printf("\tParsing error : fail to get '%s'\n",label.c_str());
      return CIO_ERROR;
    }
    rank.HeadIndex[0]=v[0];
    rank.HeadIndex[1]=v[1];
    rank.HeadIndex[2]=v[2];

    //TailIndex
    label = label_leaf + "/TailIndex";
    for (int n=0; n<3; n++) v[n]=0.0;
    if ( !(tpCntl.GetVector(label, v, 3 )) ) 
    {
      printf("\tParsing error : fail to get '%s'\n",label.c_str());
      return CIO_ERROR;
    }
    rank.TailIndex[0]=v[0];
    rank.TailIndex[1]=v[1];
    rank.TailIndex[2]=v[2];

    RankInfo.push_back(rank); 

  }

  return CIO_SUCCESS;

}

// #################################################################
// 配列形状の取り出し
std::string cio_DFI::getArrayShape()
{
  return DFI_Finfo.ArrayShape;
}

// #################################################################
// データタイプの取り出し
std::string cio_DFI::getDataType()
{
  return DFI_Finfo.DataType;
}

// #################################################################
// 成分数の取り出し
int cio_DFI::getComponent()
{
  return DFI_Finfo.Component;
}

// #################################################################
// Create Domain & Process  
void cio_DFI::cio_Create_Domain(MPI_Comm comm,
                                int G_voxel[3], int G_division[3], int head[3], int tail[3],
                                cio_Domain &G_domain, vector<cio_Rank> &G_RankInfo,
                                cio_Rank G_Rank)
{

  for(int i=0; i<3; i++ ) {
    G_domain.GlobalVoxel[i]   =G_voxel[i];
    G_domain.GlobalDivision[i]=G_division[i];
  }

  int RankID;
  MPI_Comm_rank( comm, &RankID );

  int nrank;
  MPI_Comm_size( comm, &nrank );

  int *head_All  = new int[3*nrank];
  int *tail_All  = new int[3*nrank];
  for(int i=0; i<3; i++) {
     head_All[RankID*3+i] =head[i];
     tail_All[RankID*3+i] =tail[i];
  }

  if( nrank > 1 ) {
     MPI_Allgather(head,3,MPI_INT,head_All,3,MPI_INT,comm);
     MPI_Allgather(tail,3,MPI_INT,tail_All,3,MPI_INT,comm);
  }

  for(int i=0; i<nrank; i++) {
     G_Rank.RankID=i;
     for(int j=0; j<3; j++) G_Rank.HeadIndex[j]=head_All[i*3+j];
     for(int j=0; j<3; j++) G_Rank.TailIndex[j]=tail_All[i*3+j];
     for(int j=0; j<3; j++) G_Rank.VoxelSize[j]=G_Rank.TailIndex[j]-G_Rank.HeadIndex[j]+1;
     G_RankInfo.push_back(G_Rank);
  }
}

// #################################################################
// 粗密判定
cio_EGlobalVoxel cio_DFI::CheckGlobalVoxel(int Gvoxel[3], int DFI_Gvoxel[3])
{

  if( Gvoxel[0] == DFI_Gvoxel[0]   &&
      Gvoxel[1] == DFI_Gvoxel[1]   &&
      Gvoxel[2] == DFI_Gvoxel[2]   ) return CIO_E_GV_SAME;

  if( Gvoxel[0] == DFI_Gvoxel[0]*2 &&
      Gvoxel[1] == DFI_Gvoxel[1]*2 &&
      Gvoxel[2] == DFI_Gvoxel[2]*2 ) return CIO_E_GVX2_SAME;
  
  return CIO_E_OTHER;
}

// #################################################################
// 粗密ファイルのMxN判定
bool cio_DFI::CheckMxN(vector<int> &rankList, int head[3], int tail[3], int gc, int dfi_gc)
{

  int dfi_head[3],dfi_tail[3];

  for(int i=0; i<rankList.size(); i++) {
    int n=rankList[i];
    for(int j=0; j<3; j++) {
      dfi_head[j]=RankInfo[n].HeadIndex[j]*2-1;
      dfi_tail[j]=RankInfo[n].TailIndex[j]*2;
    }

  if( head[0] != dfi_head[0] || head[1] != dfi_head[1] || head[2] != dfi_head[2] ||
      tail[0] != dfi_tail[0] || tail[1] != dfi_tail[1] || tail[2] != dfi_tail[2] ||
      gc != dfi_gc ) return false;
  }
  return true;
}

// #################################################################
// 読込み範囲を求める
bool cio_DFI::CheckReadArea(int head[3], int tail[3], int gc, int DFI_head[3],
           int DFI_tail[3], int DFI_gc, cio_EGlobalVoxel readflag, int sta[3], int end[3])
{

  int cal_gc=0;
  if( DFI_gc > 0 || gc > 0 ) {
    if( gc <= DFI_gc ) {
      cal_gc = gc;
    } else {
      cal_gc = DFI_gc;
    }
  } 

  //printf("gc : %d DFI_gc : %d cal_gc : %d\n",gc,DFI_gc,cal_gc);

  for(int i=0; i<3; i++) {
   sta[i]=0;
   end[i]=0;
  }

  for( int i=0; i<3; i++ ) {
    sta[i] = max(head[i],DFI_head[i]) -1;
    end[i] = min(tail[i],DFI_tail[i]) -1;

    //printf("i : %d sta[i] : %d end[i] : %d\n",i,sta[i],end[i]);

    if( sta[i] == 0 ) sta[i] -= cal_gc;
    else if( head[i]>DFI_head[i] ) sta[i] -= gc;

    //printf("i : %d head[i] %d DFI_head[i] %d sta[i] %d gc %d\n",
    //        i,head[i],DFI_head[i],sta[i],gc);

    if( readflag == CIO_E_GV_SAME ) {
      if( (end[i]+1) == DFI_Domain.GlobalVoxel[i] ) end[i] += cal_gc;
      else if( tail[i]<DFI_tail[i] ) end[i] += gc;
    } else {
      if( (end[i]+1) == DFI_Domain.GlobalVoxel[i]*2 ) end[i] += cal_gc;
      else if( tail[i]<DFI_tail[i] ) end[i] += gc;
    }

    //printf("i : %d tail[i] %d DFI_tail[i] %d end[i] %d gc %d\n",
    //        i,tail[i],DFI_tail[i],end[i],gc);

  }

  if( head[0] == DFI_head[0] && head[1] == DFI_head[1] && head[2] == DFI_head[2] &&
      tail[0] == DFI_tail[0] && tail[1] == DFI_tail[1] && tail[2] == DFI_tail[2] &&
      gc == DFI_gc ) return true;

  return false;

}

// #################################################################
// 読込みランクファイルリストの作成 
void cio_DFI::CreateRankList(int head[3], int tail[3], int gc, cio_EGlobalVoxel readflag,
                             vector<int> &rankList) 
{
  rankList.clear();

  int dfi_head[3],dfi_tail[3];
  int sta_x,end_x,sta_y,end_y,sta_z,end_z;

  for(int i=0; i<RankInfo.size(); i++) {

    if( readflag == CIO_E_GV_SAME ) {
      for(int n=0; n<3; n++) {
        dfi_head[n]=RankInfo[i].HeadIndex[n];
        dfi_tail[n]=RankInfo[i].TailIndex[n];
      }
    }else if( readflag == CIO_E_GVX2_SAME ) {
      for(int n=0; n<3; n++) {
        dfi_head[n]=RankInfo[i].HeadIndex[n]*2-1;
        dfi_tail[n]=RankInfo[i].TailIndex[n]*2;
      }
    }

    //printf("rank %d dfi_head : %d %d %d\n",i,dfi_head[0],dfi_head[1],dfi_head[2]);
    //printf("rank %d dfi_tail : %d %d %d\n",i,dfi_tail[0],dfi_tail[1],dfi_tail[2]);

    //x 方向のスタートエンド
    sta_x=max(head[0],dfi_head[0]);
    end_x=min(tail[0],dfi_tail[0]);

    //y 方向のスタートエンド
    sta_y=max(head[1],dfi_head[1]);
    end_y=min(tail[1],dfi_tail[1]);

    //z 方向のスタートエンド
    sta_z=max(head[2],dfi_head[2]);
    end_z=min(tail[2],dfi_tail[2]);

    if( sta_x <= end_x && sta_y <= end_y && sta_z <= end_z ) rankList.push_back(i);
  }
}

// #################################################################
// ファイル名を作成
std::string cio_DFI::Generate_FileName(int RankID, int step, const bool mio)
{

  if( DFI_Finfo.DirectoryPath.empty() ) return NULL;
  if( DFI_Finfo.Prefix.empty() ) return NULL;

  int len = DFI_Finfo.DirectoryPath.size() + DFI_Finfo.Prefix.size() +DFI_Finfo.FileFormat.size() + 25; 
  // id(6) + step(10) + 1(\0) + "_"(2) + "."(1)+"id"(2)

  char* tmp = new char[len];
  memset(tmp, 0, sizeof(char)*len);

  if( mio ) {
    sprintf(tmp, "%s/%s_%010d_id%06d.%s",DFI_Finfo.DirectoryPath.c_str(),DFI_Finfo.Prefix.c_str(), 
            step,RankID,DFI_Finfo.FileFormat.c_str());
  } else {
    sprintf(tmp, "%s/%s_%010d.%s",DFI_Finfo.DirectoryPath.c_str(),DFI_Finfo.Prefix.c_str(), 
            step,DFI_Finfo.FileFormat.c_str());
  }  
  std::string fname(tmp);
  if( tmp ) delete [] tmp;

  return fname;
}


// #################################################################
// ディレクトリがなければ作成、既存なら何もしない
int cio_DFI::MakeDirectory(string path)
{
  // 標準ライブラリ呼び出し
  // パーミッションはumaskで指定
  umask(022);

  int ret = mkdir(path.c_str(), 0777); // rwx

  if ( 0 != ret )
  {
    // 既存以外のエラー
    if ( EEXIST != errno )
    {
      printf( "\tError(errno)=[%s]\n", strerror(errno) );
      return 0;
    }
  }

  return 1;
}



// #################################################################
// debug write index.dfi & proc.dfi
bool cio_DFI::dbwrite(int RankID)
{

  if( RankID != 0 ) return true;

  printf("**** input index.dfi file infomation ****\n"); 
  printf("FileInfo {\n");
  printf("  DirectoryPath               :%s\n",DFI_Finfo.DirectoryPath.c_str()); 
  printf("  Prefix                      :%s\n",DFI_Finfo.Prefix.c_str()); 
  printf("  FileFormat                  :%s\n",DFI_Finfo.FileFormat.c_str()); 
  printf("  GuideCell                   :%d\n",DFI_Finfo.GuideCell); 
  printf("  DataType                    :%s\n",DFI_Finfo.DataType.c_str()); 
  printf("  Endian                      :%s\n",DFI_Finfo.Endian.c_str()); 
  printf("  ArrayShape                  :%s\n",DFI_Finfo.ArrayShape.c_str()); 
  printf("  Component                   :%d\n",DFI_Finfo.Component); 
  printf("}\n");
  printf("\n");

  printf("FilePath {\n");
  printf("  Process                     :%s\n",DFI_Fpath.Process.c_str());
  printf("}\n");
  printf("\n");

  printf("Unit {\n");
  printf("  Length                      :%s\n",DFI_Unit.Length.c_str());
  printf("  L0                          :%e\n",DFI_Unit.L0);
  printf("  Velocity                    :%s\n",DFI_Unit.Velocity.c_str());
  printf("  V0                          :%e\n",DFI_Unit.V0);
  printf("  Pressure                    :%s\n",DFI_Unit.Pressure.c_str());
  printf("  P0                          :%e\n",DFI_Unit.P0);
  printf("  DiffPrs                     :%e\n",DFI_Unit.DiffPrs);
  printf("  Temperatur                  :%s\n",DFI_Unit.Temperatur.c_str());
  printf("  BaseTemp                    :%e\n",DFI_Unit.BaseTemp);
  printf("  DiffTemp                    :%e\n",DFI_Unit.DiffTemp);
  printf("}\n");
  printf("\n");

  printf("TimeSlice {\n");
  for(int i=0; i<TimeSlice.size(); i++) {
    printf("  Slice[%d] {\n",i);
    printf("    Step                        :%d\n",TimeSlice[i].step);
    printf("    Time                        :%e\n",TimeSlice[i].time);
    for(int j=0; j<TimeSlice[i].Min.size(); j++) {
      printf("    MinMax[%d] {\n",j+1);
      printf("      Min                         :%e\n",TimeSlice[i].Min[j]);
      printf("      Max                         :%e\n",TimeSlice[i].Max[j]);
      printf("    }\n");
    }
    printf("  }\n");
  }
  printf("}\n");
  printf("\n");

  printf("**** input proc.dfi file infomation ****\n"); 
  printf("Domian {\n");
  printf("  GlovalOrigin                :%e %e %e\n",DFI_Domain.GlobalOrigin[0],
                                                     DFI_Domain.GlobalOrigin[1],
                                                     DFI_Domain.GlobalOrigin[2]); 
  printf("  GlovalRegion                :%e %e %e\n",DFI_Domain.GlobalRegion[0],
                                                     DFI_Domain.GlobalRegion[1],
                                                     DFI_Domain.GlobalRegion[2]); 
  printf("  GlobalVoxel                 :%d %d %d\n",DFI_Domain.GlobalVoxel[0],
                                                     DFI_Domain.GlobalVoxel[1],
                                                     DFI_Domain.GlobalVoxel[2]); 
  printf("  GlobalDivision              :%d %d %d\n",DFI_Domain.GlobalDivision[0],
                                                     DFI_Domain.GlobalDivision[1],
                                                     DFI_Domain.GlobalDivision[2]); 
  printf("}\n");
  printf("\n");

  printf("MPI {\n");
  printf("  NumberofRank                :%d\n",DFI_MPI.NumberOfRank); 
  printf("  NumberofGroup               :%d\n",DFI_MPI.NumberOfGroup); 
  printf("}\n");
  printf("\n");

  printf("Process {\n");
  for (int i=0; i<RankInfo.size(); i++) {
    printf("  Rank[%d] {\n",i);
    printf("    ID                          :%d\n",RankInfo[i].RankID);
    printf("    HostName                    :%s\n",RankInfo[i].HostName.c_str()); 
    printf("    VoxelSize                   :%d %d %d\n",RankInfo[i].VoxelSize[0],RankInfo[i].VoxelSize[1],RankInfo[i].VoxelSize[2]); 
    printf("    HeadIndex                   :%d %d %d\n",RankInfo[i].HeadIndex[0],RankInfo[i].HeadIndex[1],RankInfo[i].HeadIndex[2]); 
    printf("    TailIndex                   :%d %d %d\n",RankInfo[i].TailIndex[0],RankInfo[i].TailIndex[1],RankInfo[i].TailIndex[2]); 
    printf("  }\n");
  }
  printf("}\n");

  return true;

}

