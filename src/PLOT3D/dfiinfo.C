//##################################################################################
//
// Flow Base class
//
// Copyright (c) 2007-2011 VCAD System Research Program, RIKEN.
// All rights reserved.
//
// Copyright (c) 2011-2014 Institute of Industrial Science, The University of Tokyo.
// All rights reserved.
//
// Copyright (c) 2012-2014 Advanced Institute for Computational Science, RIKEN.
// All rights reserved.
//
//##################################################################################

/** 
 * @file   dfiinfo.C
 * @brief  DfiInfo Class
 * @author kero
 */

#include "dfiinfo.h"

// #################################################################
// コンストラクタ
DfiInfo::DfiInfo()
{
  GuideCell=0;
  Component=0;
  
  L0=0.0;
  V0=0.0;
  P0=0.0;
  DiffPrs=0.0;
  
  Node=NULL;
  //RankID_in_MPIworld=0;
  //GroupID_in_MPIworld=0;
  NumberOfRank=0;
  NumberOfGroup=0;
  for(int i=0;i<3;i++) Global_Voxel[i]=0;
  for(int i=0;i<3;i++) Global_Division[i]=0;
  
  Sc.clear();
  
  dim=0;
  NodeInfoSize=0;
  index_y.clear();
  index_z.clear();
}


// #################################################################
// デストラクタ
DfiInfo::~DfiInfo()
{
  if (Node) delete[] Node;
}


// #################################################################
// 
void DfiInfo::ReadDfiFile(string fname)
{
  TextParser tpCntl;
  
  string str;
  string label,label_base,label_leaf;
  int ct;
  REAL_TYPE ct2;
  int nnode=0;
  
  
  //入力ファイルをセット
  int ierror = tpCntl.read(fname);
  if ( ierror )
  {
    Hostonly_ stamped_printf("\tinput file not found '%s'\n",fname.c_str());
    Exit(0);
  }
  
  
  //FileInfo
  
  //Prefix
  //label = "/DistributedFileInfo/Prefix";
  label = "/FileInfo/Prefix";
  if ( !(tpCntl.getInspectedValue(label, str )) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n",label.c_str());
    Exit(0);
  }
  Prefix=str;
  
  //FileFormat
  //label = "/DistributedFileInfo/FileFormat";
  label = "/FileInfo/FileFormat";
  if ( !(tpCntl.getInspectedValue(label, str )) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n",label.c_str());
    Exit(0);
  }
  FileFormat=str;
  
  //GuideCell
  //label = "/DistributedFileInfo/GuideCell";
  label = "/FileInfo/GuideCell";
  if ( !(tpCntl.getInspectedValue(label, ct )) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n",label.c_str());
    Exit(0);
  }
  else {
    GuideCell = ct;
  }
  
  //DataType
  label = "/FileInfo/DataType";
  if ( !(tpCntl.getInspectedValue(label, str )) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n",label.c_str());
    Exit(0);
  }
  DataType=str;
  
  //Endian
  label = "/FileInfo/Endian";
  if ( !(tpCntl.getInspectedValue(label, str )) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n",label.c_str());
    Exit(0);
  }
  Endian=str;
  
  //ArrayShape
  label = "/FileInfo/ArrayShape";
  if ( !(tpCntl.getInspectedValue(label, str )) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n",label.c_str());
    Exit(0);
  }
  ArrayShape=str;
  
  //Component
  label = "/FileInfo/Component";
  if ( !(tpCntl.getInspectedValue(label, ct )) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n",label.c_str());
    Exit(0);
  }
  else {
    Component = ct;
  }
  
  //Unit
  
  //Length
  label = "/Unit/Length";
  if ( !(tpCntl.getInspectedValue(label, str )) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n",label.c_str());
    Exit(0);
  }
  Length=str;
  
  //L0
  label = "/Unit/L0";
  if ( !(tpCntl.getInspectedValue(label, ct2 )) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n",label.c_str());
    Exit(0);
  }
  else {
    L0 = ct2;
  }
  
  //Velocity
  label = "/Unit/Velocity";
  if ( !(tpCntl.getInspectedValue(label, str )) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n",label.c_str());
    Exit(0);
  }
  Velocity=str;
  
  //V0
  label = "/Unit/V0";
  if ( !(tpCntl.getInspectedValue(label, ct2 )) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n",label.c_str());
    Exit(0);
  }
  else {
    V0 = ct2;
  }
  
  //Length
  label = "/Unit/Pressure";
  if ( !(tpCntl.getInspectedValue(label, str )) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n",label.c_str());
    Exit(0);
  }
  Pressure=str;
  
  //P0
  label = "/Unit/P0";
  if ( !(tpCntl.getInspectedValue(label, ct2 )) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n",label.c_str());
    Exit(0);
  }
  else {
    P0 = ct2;
  }
  
  //DiffPrs
  label = "/Unit/DiffPrs";
  if ( !(tpCntl.getInspectedValue(label, ct2 )) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n",label.c_str());
    Exit(0);
  }
  else {
    DiffPrs = ct2;
  }
  
  
  //FilePath
  label = "/FilePath/Process";
  if ( !(tpCntl.getInspectedValue(label, str )) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n",label.c_str());
    Exit(0);
  }
  dfi_proc=str;
  
  //read process
  ReadDfiProc(dfi_proc);
  
  
  //TimeSlice
  nnode=0;
  label_base = "/TimeSlice/Slice";
  if ( tpCntl.chkNode(label_base) )  //nodeがあれば
  {
    nnode = tpCntl.countLabels(label_base);
  }
  
  label_base = "/TimeSlice";
  for (int i=0; i<nnode; i++)
  {
    if ( !(tpCntl.getNodeStr(label_base, i+1, str)) )
    {
      printf("\tParsing error : No Elem name\n");
      Exit(0);
    }
    
    if( strcasecmp(str.substr(0,5).c_str(), "Slice") ) continue;
    
    //step
    label=label_base+"/"+str+"/Step";
    if ( !(tpCntl.getInspectedValue(label, ct )) )
    {
      Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n",label.c_str());
      Exit(0);
    }
    
    //time
    label=label_base+"/"+str+"/Time";
    if ( !(tpCntl.getInspectedValue(label, ct2 )) )
    {
      Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n",label.c_str());
      Exit(0);
    }
    
    SetSlice(ct,ct2);
  }
  
  //TextParserの破棄
  tpCntl.remove();
  
  //内部変数の計算
  SetValue();
  
  
#if 0
  cout << endl;
  cout << endl;
  //for(int i=0;i<ndfi;i++){
  cout << "" << endl;
  cout << "Prefix = " << this->Prefix << endl;
  //cout << "RankIDinMPIworld = " << this->RankID_in_MPIworld << endl;
  //cout << "GroupIDinMPIworld = " << this->GroupID_in_MPIworld << endl;
  cout << "NumberOfRankIn = " << this->NumberOfRank << endl;
  cout << "NumberOfGroup  = " << this->NumberOfGroup << endl;
  cout << "GlobalOrigin[0] = " << this->Global_Origin[0] << endl;
  cout << "GlobalOrigin[1] = " << this->Global_Origin[1] << endl;
  cout << "GlobalOrigin[2] = " << this->Global_Origin[2] << endl;
  cout << "GlobalRegion[0] = " << this->Global_Region[0] << endl;
  cout << "GlobalRegion[1] = " << this->Global_Region[1] << endl;
  cout << "GlobalRegion[2] = " << this->Global_Region[2] << endl;
  cout << "GlobalVoxel[0] = " << this->Global_Voxel[0] << endl;
  cout << "GlobalVoxel[1] = " << this->Global_Voxel[1] << endl;
  cout << "GlobalVoxel[2] = " << this->Global_Voxel[2] << endl;
  cout << "GlobalDivision[0] = " << this->Global_Division[0] << endl;
  cout << "GlobalDivision[1] = " << this->Global_Division[1] << endl;
  cout << "GlobalDivision[2] = " << this->Global_Division[2] << endl;
  cout << "FileFormat = " << this->FileFormat << endl;
  cout << "GuideCell = " << this->GuideCell << endl;
  cout << "" << endl;
  cout << "NodeInfoSize = " << this->NodeInfoSize << endl;
  for(int j=0; j< this->NodeInfoSize; j++ ) {
    cout << "" << endl;
    cout << "Node[" << j << "].RankID = " << this->Node[j].RankID << endl;
    cout << "Node[" << j << "].HostName = " << this->Node[j].HostName << endl;
    cout << "Node[" << j << "].VoxelSize[0] = " << this->Node[j].VoxelSize[0] << endl;
    cout << "Node[" << j << "].VoxelSize[1] = " << this->Node[j].VoxelSize[1] << endl;
    cout << "Node[" << j << "].VoxelSize[2] = " << this->Node[j].VoxelSize[2] << endl;
    cout << "Node[" << j << "].HeadIndex[0] = " << this->Node[j].HeadIndex[0] << endl;
    cout << "Node[" << j << "].HeadIndex[1] = " << this->Node[j].HeadIndex[1] << endl;
    cout << "Node[" << j << "].HeadIndex[2] = " << this->Node[j].HeadIndex[2] << endl;
    cout << "Node[" << j << "].TailIndex[0] = " << this->Node[j].TailIndex[0] << endl;
    cout << "Node[" << j << "].TailIndex[1] = " << this->Node[j].TailIndex[1] << endl;
    cout << "Node[" << j << "].TailIndex[2] = " << this->Node[j].TailIndex[2] << endl;
    cout << "Node[" << j << "].IJK = " << this->Node[j].IJK << endl;
    cout << "Node[" << j << "].IJK_JK = " << this->Node[j].IJK_JK << endl;
    cout << "Node[" << j << "].IJK_K = " << this->Node[j].IJK_K << endl;
  }
  cout << "" << endl;
  //cout << "step.size() = " << this->step.size() << endl;
  //for(int j=0; j< this->step.size(); j++ ) {
  //  cout << "step[" << j << "] = " << this->step[j] << endl;
  //}
  cout << "Sc.size() = " << this->Sc.size() << endl;
  for(int j=0; j< this->Sc.size(); j++ ) {
    cout << "step[" << j << "] = " << this->Sc[j]->step << endl;
    cout << "time[" << j << "] = " << this->Sc[j]->time << endl;
  }
  cout << "" << endl;
  cout << "index_y.size() = " << this->index_y.size() << endl;
  for(int j=0; j< this->index_y.size(); j++ ) {
    cout << "index_y[" << j << "] = " << this->index_y[j] << endl;
  }
  cout << "" << endl;
  cout << "index_z.size() = " << this->index_z.size() << endl;
  for(int j=0; j< this->index_z.size(); j++ ) {
    cout << "index_z[" << j << "] = " << this->index_z[j] << endl;
  }
  //}
  cout << endl;
  cout << endl;
  
  //Exit(0);
  
#endif
  
}


// #################################################################
//
void DfiInfo::SetSlice(int m_step, REAL_TYPE m_time)
{
  Slice* s = new DfiInfo::Slice;
  s->step=m_step;
  s->time=m_time;
  DfiInfo::Sc.push_back(s);
}


// #################################################################
//
void DfiInfo::ReadDfiProc(string fname)
{
  TextParser tpCntl;
  string str;
  string label,label_base,label_leaf;
  int ct;
  int nnode=0;
  int iv[3];
  REAL_TYPE v[3];
  
  
  //入力ファイルをセット
  int ierror = tpCntl.read(fname);
  if ( ierror )
  {
    Hostonly_ stamped_printf("\tinput file not found '%s'\n",fname.c_str());
    Exit(0);
  }
  
  
  //Domain
  
  //Global_Origin
  label = "/Domain/GlobalOrigin";
  for (int n=0; n<3; n++) v[n]=0.0;
  if ( !(tpCntl.getInspectedVector(label, v, 3 )) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n",label.c_str());
    Exit(0);
  }
  Global_Origin[0]=v[0];
  Global_Origin[1]=v[1];
  Global_Origin[2]=v[2];
  
  //Global_Region
  label = "/Domain/GlobalRegion";
  for (int n=0; n<3; n++) v[n]=0.0;
  if ( !(tpCntl.getInspectedVector(label, v, 3 )) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n",label.c_str());
    Exit(0);
  }
  Global_Region[0]=v[0];
  Global_Region[1]=v[1];
  Global_Region[2]=v[2];
  
  //Global_Voxel
  label = "/Domain/GlobalVoxel";
  for (int n=0; n<3; n++) iv[n]=0;
  if ( !(tpCntl.getInspectedVector(label, iv, 3 )) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n",label.c_str());
    Exit(0);
  }
  Global_Voxel[0]=iv[0];
  Global_Voxel[1]=iv[1];
  Global_Voxel[2]=iv[2];
  
  //Global_Division
  label = "/Domain/GlobalDivision";
  for (int n=0; n<3; n++) iv[n]=0.0;
  if ( !(tpCntl.getInspectedVector(label, iv, 3 )) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n",label.c_str());
    Exit(0);
  }
  Global_Division[0]=iv[0];
  Global_Division[1]=iv[1];
  Global_Division[2]=iv[2];
  
  
  //MPI
  
  //NumberOfRank
  label = "/MPI/NumberOfRank";
  if ( !(tpCntl.getInspectedValue(label, ct )) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n",label.c_str());
    Exit(0);
  }
  else {
    NumberOfRank = ct;
  }
  
  //NumberOfGroup
  label = "/MPI/NumberOfGroup";
  if ( !(tpCntl.getInspectedValue(label, ct )) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n",label.c_str());
    Exit(0);
  }
  else {
    NumberOfGroup = ct;
  }
  
  
  ////RankID_in_MPIworld
  //label = "/DistributedFileInfo/RankIDinMPIworld";
  //if ( !(tpCntl.getInspectedValue(label, &ct )) ) {
  //  Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n",label.c_str());
  //  Exit(0);
  //}
  //else {
  //  RankID_in_MPIworld = ct;
  //}
  
  ////GroupID_in_MPIworld
  //label = "/DistributedFileInfo/GroupIDinMPIworld";
  //if ( !(tpCntl.getInspectedValue(label, &ct )) ) {
  //  Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n",label.c_str());
  //  Exit(0);
  //}
  //else {
  //  GroupID_in_MPIworld = ct;
  //}
  
  
  //Process <--- NodeInfo
  
  nnode=0;
  //label_base = "/DistributedFileInfo/NodeInfo";
  label_base = "/Process";
  if ( tpCntl.chkNode(label_base) )  //nodeがあれば
  {
    nnode = tpCntl.countLabels(label_base);
  }
  NodeInfoSize=nnode;
  Node = new DfiInfo::NodeInfo[nnode];
  
  for (int i=0; i<NodeInfoSize; i++)
  {
    if ( !(tpCntl.getNodeStr(label_base, i+1, str)) )
    {
      printf("\tParsing error : No Elem name\n");
      Exit(0);
    }
    if( strcasecmp(str.substr(0,4).c_str(), "Rank") ) continue;
    label_leaf=label_base+"/"+str;
    
    //RankID
    label = label_leaf + "/ID";
    if ( !(tpCntl.getInspectedValue(label, ct )) )
    {
      Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n",label.c_str());
      Exit(0);
    }
    else {
      Node[i].RankID= ct;
    }
    
    //HostName
    label = label_leaf + "/HostName";
    if ( !(tpCntl.getInspectedValue(label, str )) )
    {
      Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
      Exit(0);
    }
    Node[i].HostName= str;
    
    //VoxelSize
    label = label_leaf + "/VoxelSize";
    for (int n=0; n<3; n++) v[n]=0.0;
    if ( !(tpCntl.getInspectedVector(label, v, 3 )) )
    {
      Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n",label.c_str());
      Exit(0);
    }
    Node[i].VoxelSize[0]=v[0];
    Node[i].VoxelSize[1]=v[1];
    Node[i].VoxelSize[2]=v[2];
    
    //HeadIndex
    label = label_leaf + "/HeadIndex";
    for (int n=0; n<3; n++) v[n]=0.0;
    if ( !(tpCntl.getInspectedVector(label, v, 3 )) )
    {
      Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n",label.c_str());
      Exit(0);
    }
    Node[i].HeadIndex[0]=v[0];
    Node[i].HeadIndex[1]=v[1];
    Node[i].HeadIndex[2]=v[2];
    
    //TailIndex
    label = label_leaf + "/TailIndex";
    for (int n=0; n<3; n++) v[n]=0.0;
    if ( !(tpCntl.getInspectedVector(label, v, 3 )) )
    {
      Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n",label.c_str());
      Exit(0);
    }
    Node[i].TailIndex[0]=v[0];
    Node[i].TailIndex[1]=v[1];
    Node[i].TailIndex[2]=v[2];
    
  }
  
  
  //TextParserの破棄
  tpCntl.remove();
  
}



// #################################################################
// 
void DfiInfo::SetValue()
{
  //一次元アドレス変換の保持
  long NI=(long)this->Global_Voxel[0];
  long NJ=(long)this->Global_Voxel[1];
  long NK=(long)this->Global_Voxel[2];
  for(int j=0; j< this->NodeInfoSize; j++ ) {
    long hi=(long)this->Node[j].HeadIndex[0];
    long hj=(long)this->Node[j].HeadIndex[1];
    long hk=(long)this->Node[j].HeadIndex[2];

    this->Node[j].IJK    = (hk-1)*NI*NJ+(hj-1)*NI+(hi-1)+1;
    this->Node[j].IJK_JK = (hk-1)*NI*NJ+(hj-1)*NI;
    this->Node[j].IJK_K  = (hk-1)*NI*NJ;
  }

  //同一IJK_JKを持つNodeのindexの作成
  this->index_y.clear();
  long keep_IJK_JK=0;
  this->index_y.push_back(0);
  for(int j=0; j< this->NodeInfoSize; j++ ) {
    if( this->Node[j].IJK_JK == keep_IJK_JK ){
    }else{
      keep_IJK_JK=this->Node[j].IJK_JK;
      this->index_y.push_back(j);
    }
  }
  this->index_y.push_back(this->NodeInfoSize);

  //同一IJK_Kを持つNodeのindexの作成
  this->index_z.clear();
  long keep_IJK_K=0;
  this->index_z.push_back(0);
  for(int j=0; j< this->NodeInfoSize; j++ ) {
    if( this->Node[j].IJK_K == keep_IJK_K ){
    }else{
      keep_IJK_K=this->Node[j].IJK_K;
      this->index_z.push_back(j);
    }
  }
  this->index_z.push_back(this->NodeInfoSize);

}

