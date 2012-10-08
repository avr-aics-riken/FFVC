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
 * @file   dfiinfo.C
 * @brief  DfiInfo Class
 * @author kero
 */

#include "dfiinfo.h"

// #################################################################
// コンストラクタ
DfiInfo::DfiInfo()
{
  RankID_in_MPIworld=0;
  GroupID_in_MPIworld=0;
  Number_of_Rank_in_MPIworld=0;
  Number_of_Group_in_MPIworld=0;
  for(int i=0;i<3;i++) Global_Voxel[i]=0;
  for(int i=0;i<3;i++) Global_Division[i]=0;
  GuideCell=0;
  NodeInfoSize=0;

  Node=NULL;

  step.clear();
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
  TPControl tpCntl;
  string str;
  string label,label_base,label_leaf;
  int ct;
  REAL_TYPE v[3];
  int nnode=0;

  //TPインスタンス生成
  tpCntl.getTPinstance();

  //入力ファイルをセット
  int ierror = tpCntl.readTPfile(fname);
  if ( ierror )
  {
    Hostonly_ stamped_printf("\tinput file not found '%s'\n",fname.c_str());
    Exit(0);
  }

  //Prefix
  label = "/Distributed_File_Info/Prefix";
  if ( !(tpCntl.GetValue(label, &str )) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n",label.c_str());
    Exit(0);
  }
  Prefix=str;

  //RankID_in_MPIworld
  label = "/Distributed_File_Info/RankID_in_MPIworld";
  if ( !(tpCntl.GetValue(label, &ct )) ) {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n",label.c_str());
    Exit(0);
  }
  else {
    RankID_in_MPIworld = ct;
  }

  //GroupID_in_MPIworld
  label = "/Distributed_File_Info/GroupID_in_MPIworld";
  if ( !(tpCntl.GetValue(label, &ct )) ) {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n",label.c_str());
    Exit(0);
  }
  else {
    GroupID_in_MPIworld = ct;
  }

  //Number_of_Rank_in_MPIworld
  label = "/Distributed_File_Info/Number_of_Rank_in_MPIworld";
  if ( !(tpCntl.GetValue(label, &ct )) ) {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n",label.c_str());
    Exit(0);
  }
  else {
    Number_of_Rank_in_MPIworld = ct;
  }

  //Number_of_Group_in_MPIworld
  label = "/Distributed_File_Info/Number_of_Group_in_MPIworld";
  if ( !(tpCntl.GetValue(label, &ct )) ) {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n",label.c_str());
    Exit(0);
  }
  else {
    Number_of_Group_in_MPIworld = ct;
  }

  //Global_Voxel
  label = "/Distributed_File_Info/Global_Voxel";
  for (int n=0; n<3; n++) v[n]=0.0;
  if ( !(tpCntl.GetVector(label, v, 3 )) ) 
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n",label.c_str());
    Exit(0);
  }
  Global_Voxel[0]=v[0];
  Global_Voxel[1]=v[1];
  Global_Voxel[2]=v[2];

  //Global_Division
  label = "/Distributed_File_Info/Global_Division";
  for (int n=0; n<3; n++) v[n]=0.0;
  if ( !(tpCntl.GetVector(label, v, 3 )) ) 
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n",label.c_str());
    Exit(0);
  }
  Global_Division[0]=v[0];
  Global_Division[1]=v[1];
  Global_Division[2]=v[2];

  //FileFormat
  label = "/Distributed_File_Info/FileFormat";
  if ( !(tpCntl.GetValue(label, &str )) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n",label.c_str());
    Exit(0);
  }
  FileFormat=str;

  //GuideCell
  label = "/Distributed_File_Info/GuideCell";
  if ( !(tpCntl.GetValue(label, &ct )) ) {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n",label.c_str());
    Exit(0);
  }
  else {
    GuideCell = ct;
  }

  //NodeInfo
  nnode=0;
  label_base = "/Distributed_File_Info/NodeInfo";
  if ( tpCntl.chkNode(label_base) )  //nodeがあれば
  {
    nnode = tpCntl.countLabels(label_base);
  }
  NodeInfoSize=nnode;
  Node = new DfiInfo::NodeInfo[nnode];

  for (int i=0; i<NodeInfoSize; i++) {

    if(!tpCntl.GetNodeStr(label_base,i+1,&str))
    {
      printf("\tParsing error : No Elem name\n");
      Exit(0);
    }
    if( strcasecmp(str.substr(0,4).c_str(), "node") ) continue;
    label_leaf=label_base+"/"+str;

	//RankID
    label = label_leaf + "/RankID";
    if ( !(tpCntl.GetValue(label, &ct )) ) {
      Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n",label.c_str());
      Exit(0);
    }
    else {
      Node[i].RankID= ct;
    }

	//HostName
    label = label_leaf + "/HostName";
    if ( !(tpCntl.GetValue(label, &str )) ) {
      Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
      Exit(0);
    }
    Node[i].HostName= str;

    //VoxelSize
    label = label_leaf + "/VoxelSize";
    for (int n=0; n<3; n++) v[n]=0.0;
    if ( !(tpCntl.GetVector(label, v, 3 )) ) 
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
    if ( !(tpCntl.GetVector(label, v, 3 )) ) 
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
    if ( !(tpCntl.GetVector(label, v, 3 )) ) 
    {
      Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n",label.c_str());
      Exit(0);
    }
    Node[i].TailIndex[0]=v[0];
    Node[i].TailIndex[1]=v[1];
    Node[i].TailIndex[2]=v[2];

  }

  //FileInfo
  nnode=0;
  label_base = "/Distributed_File_Info/FileInfo";
  if ( tpCntl.chkNode(label_base) )  //nodeがあれば
  {
    nnode = tpCntl.countLabels(label_base);
  }

  step.clear();
  for (int i=0; i<nnode; i++) {

    if(!tpCntl.GetNodeStr(label_base,i+1,&str))
    {
      printf("\tParsing error : No Elem name\n");
      Exit(0);
    }
    if( strcasecmp(str.substr(0,4).c_str(), "step") ) continue;

    //step
    label=label_base+"/"+str;
    if ( !(tpCntl.GetValue(label, &ct )) ) {
      Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n",label.c_str());
      Exit(0);
    }
    else {
      step.push_back(ct);
    }
  }

  //TextParserの破棄
  tpCntl.remove();

  //内部変数の計算
  SetValue();

#if 0
  cout << endl;
  cout << endl;
  for(int i=0;i<ndfi;i++){
    cout << "" << endl;
    cout << "Prefix = " << this->Prefix << endl;
    cout << "RankID_in_MPIworld = " << this->RankID_in_MPIworld << endl;
    cout << "GroupID_in_MPIworld = " << this->GroupID_in_MPIworld << endl;
    cout << "Number_of_Rank_in_MPIworld = " << this->Number_of_Rank_in_MPIworld << endl;
    cout << "Number_of_Group_in_MPIworld = " << this->Number_of_Group_in_MPIworld << endl;
    cout << "Global_Voxel[0] = " << this->Global_Voxel[0] << endl;
    cout << "Global_Voxel[1] = " << this->Global_Voxel[1] << endl;
    cout << "Global_Voxel[2] = " << this->Global_Voxel[2] << endl;
    cout << "Global_Division[0] = " << this->Global_Division[0] << endl;
    cout << "Global_Division[1] = " << this->Global_Division[1] << endl;
    cout << "Global_Division[2] = " << this->Global_Division[2] << endl;
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
    cout << "step.size() = " << this->step.size() << endl;
    for(int j=0; j< this->step.size(); j++ ) {
      cout << "step[" << j << "] = " << this->step[j] << endl;
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
  }
  cout << endl;
  cout << endl;
#endif

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

