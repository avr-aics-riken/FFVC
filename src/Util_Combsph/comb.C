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
 * @file   comb.C
 * @brief  COMB Class
 * @author kero
 */

#include "comb.h"

// #################################################################
// コンストラクタ
COMB::COMB()
{
  procGrp = 0;
  myRank  = -1;
  numProc = 0;
  
  filename="";
  out_dirname="./";
  in_dirname="./";
  pflag=0;
  pflagv=0;
  lflag=0;
  lflagv=0;
  thin_out = false;
  
  P3Op.IS_xyz = ON;
  P3Op.IS_q = OFF;
  P3Op.IS_funciton = ON;
  P3Op.IS_function_name = ON;
  P3Op.IS_fvbnd = OFF;
  P3Op.IS_DivideFunc=OFF;
  P3Op.ngrid=0; //出力ブロック数
  P3Op.nvar=0;  //出力項目数
  
  output_real_type=0;
  out_format=0;
  ndfi=0;
  DI=NULL;
  
  // dfi管理
  for (int i=0; i<var_END; i++) dfi_mng[i]=0;
  
  staging=0;
  
}


// #################################################################
// デストラクタ
COMB::~COMB()
{
  if( DI ) delete [] DI;
}

// #################################################################
//
void COMB::ReadInit(string input_file)
{
  
  // ------------------------------------
  FILE* fp = NULL;
  
  // TPインスタンス生成
  TPControl tpCntl;
  tpCntl.getTPinstance();
  
  //入力ファイルをセット
  int ierror = tpCntl.readTPfile(input_file);
  
  //入力ファイルの読み込み--->パラメータのセット
  ReadInputFile(&tpCntl);
  
  //TextParserの破棄
  tpCntl.remove();
  
  return;
}

// #################################################################
//
void COMB::ReadInputFile(TPControl* tpCntl)
{
  string str;
  string label,label_base,label_leaf;
  
  // node数の取得
  int nnode=0;
  label_base = "/CombData";
  if ( tpCntl->chkNode(label_base) )  //nodeがあれば
  {
    nnode = tpCntl->countLabels(label_base);
  }
  
  // dfi_nameの取得
  dfi_name.clear();
  label_base = "/CombData";
  for (int i=0; i<nnode; i++) {
    
    if(!tpCntl->GetNodeStr(label_base,i+1,&str))
    {
      printf("\tParsing error : No Elem name\n");
      Exit(0);
    }
    if( strcasecmp(str.substr(0,4).c_str(), "list") ) continue;
    label=label_base+"/"+str;
    
    if ( !(tpCntl->GetValue(label, &str )) ) {
      printf("\tParsing error : fail to get '%s'\n", label.c_str());
      Exit(0);
    }
    //FList[ilist].name = str;
    dfi_name.push_back(str.c_str());
    
  }
  
#if 0
  cout << "dfi_name.size() = " << dfi_name.size() << endl;
  vector<string>::const_iterator it;
  for (it = dfi_name.begin(); it != dfi_name.end(); it++) {
    cout << "name = " << (*it).c_str() << endl;
  }
#endif
  
  
  //出力ディレクトリの指定 ---> 実行オプションよりこちらが優先される
  label = "/CombData/OutputDir";
  if ( (tpCntl->GetValue(label, &str )) )
  {
    out_dirname=str;
    LOG_OUT_ fprintf(fplog,"\tReset Output Directory '%s'\n", out_dirname.c_str());
    STD_OUT_ printf("\tReset Output Directory '%s'\n", out_dirname.c_str());
    CheckDir(out_dirname);
    if( out_dirname.size() != 0 ) out_dirname=out_dirname+"/";
  }
  
  
  //入力ディレクトリの指定
  label = "/CombData/InputDir";
  if ( (tpCntl->GetValue(label, &str )) )
  {
    in_dirname=str;
    LOG_OUT_ fprintf(fplog,"\tReset Input Directory '%s'\n", in_dirname.c_str());
    STD_OUT_ printf("\tReset Input Directory '%s'\n", in_dirname.c_str());
    CheckDir(in_dirname);
    if( in_dirname.size() != 0 ) in_dirname=in_dirname+"/";
  }
  
  
  //並列実行時のSTAGINGのON/OFF
  label = "/CombData/Staging";
  if ( !(tpCntl->GetValue(label, &str )) )
  {
    //Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    //Exit(0);
    staging = OFF;
  }
  else
  {
    if     ( !strcasecmp(str.c_str(), "on") )  staging = ON;
    else if( !strcasecmp(str.c_str(), "off") ) staging = OFF;
    else
    {
      Hostonly_ stamped_printf("\tInvalid keyword is described for '%s'\n", label.c_str());
      Exit(0);
    }
  }
  
  
  //出力の単精度or倍精度指定 ---> PLOT3Dの場合は、optionに記述があればそちらを優先
  label = "/CombData/OutputRealType";
  if ( !(tpCntl->GetValue(label, &str )) )
  {
    //printf("\tParsing error : fail to get '%s'\n", label.c_str());
    //Exit(0);
    output_real_type = OUTPUT_REAL_UNKNOWN;
  }
  if     ( !strcasecmp(str.c_str(), "float" ) )  output_real_type = OUTPUT_FLOAT;
  else if( !strcasecmp(str.c_str(), "double" ) ) output_real_type = OUTPUT_DOUBLE;
  else
  {
    printf("\tInvalid keyword is described for '%s'\n", label.c_str());
    Exit(0);
  }
  
  // 連結ファイルの出力フォーマット
  label = "/CombData/OutFormat";
  if ( !(tpCntl->GetValue(label, &str )) )
  {
    printf("\tParsing error : fail to get '%s'\n", label.c_str());
    Exit(0);
  }
  if     ( !strcasecmp(str.c_str(), "sph" ) )    out_format = OUTFORMAT_IS_SPH;
  else if( !strcasecmp(str.c_str(), "plot3d" ) ) out_format = OUTFORMAT_IS_PLOT3D;
  else
  {
    printf("\tInvalid keyword is described for  '%s'\n", label.c_str());
    Exit(0);
  }
  
  // PLOT3Dオプションの読み込み
  if( out_format == OUTFORMAT_IS_PLOT3D ) get_PLOT3D(tpCntl);
  
}


// #################################################################
// PLOT3Dファイル入出力に関するパラメータ
void COMB::get_PLOT3D(TPControl* tpCntl)
{
  string str;
  string label;
  
  // Filename
  //label = "/Plot3dOptions/Filename";
  //
  //if ( !(tpCntl->GetValue(label, &str )) )
  //{
  //  P3Op.basename = "PLOT3Doutput";
  //}
  //else
  //{
  //  P3Op.basename = str;
  //}
  
  // FileNameGrid --- option
  label = "/Steer/Plot3dOptions/FileNameGrid";
  
  if ( !(tpCntl->GetValue(label, &str)) )
  {
    P3Op.basename_g = "PLOT3DoutputGrid";
  }
  else
  {
    P3Op.basename_g = str;
  }
  if ( P3Op.basename_g.empty() )
  {
    P3Op.basename_g = "PLOT3DoutputGrid";
  }
  
  // FileNameFunc --- option
  label = "/Steer/Plot3dOptions/FileNameFunc";
  
  if ( !(tpCntl->GetValue(label, &str)) )
  {
    P3Op.basename_f = "PLOT3Doutput";
  }
  else
  {
    P3Op.basename_f = str;
  }
  if ( P3Op.basename_f.empty() )
  {
    P3Op.basename_f = "PLOT3Doutput";
  }
  
  
  /* GridKind
   label = "/Plot3dOptions/GridKind";
   
   if ( !(tpCntl->GetValue(label, &str )) )
   {
   printf("\tParsing error : fail to get '%s'\n", label.c_str());
   Exit(0);
   }
   else
   {
   if     ( !strcasecmp(str.c_str(), "single_grid") ) FP3DR.setSingleGrid();
   else if( !strcasecmp(str.c_str(), "multi_grid") )  FP3DR.setMultiGrid();
   else
   {
   printf("\tInvalid keyword is described for '%s'\n", label.c_str());
   Exit(0);
   }
   }
   */
  FP3DR.setMultiGrid(); // 常にmulti grid
  
  // 格子の移動
  label = "/Plot3dOptions/GridMobility";
  
  if ( !(tpCntl->GetValue(label, &str )) )
  {
    printf("\tParsing error : fail to get '%s'\n", label.c_str());
    Exit(0);
  }
  else
  {
    if     ( !strcasecmp(str.c_str(), "movable") )   FP3DR.setMoveGrid(GRID_MOVE);
    else if( !strcasecmp(str.c_str(), "immovable") ) FP3DR.setMoveGrid(GRID_NOT_MOVE);
    else
    {
      printf("\tInvalid keyword is described for '%s'\n", label.c_str());
      Exit(0);
    }
  }
  
  
  // 時間方向の変化
  label = "/Plot3dOptions/StateOfTime";
  
  if ( !(tpCntl->GetValue(label, &str )) )
  {
    FP3DR.setSteady(FB_UNSTEADY);
  }
  else
  {
    if     ( !strcasecmp(str.c_str(), "steady") )   FP3DR.setSteady(FB_STEADY);
    else if( !strcasecmp(str.c_str(), "unsteady") ) FP3DR.setSteady(FB_UNSTEADY);
    else
    {
      printf("\tInvalid keyword is described for '%s'\n", label.c_str());
      Exit(0);
    }
  }
  
  
  /* IBLANKファイル
   label = "/Plot3dOptions/SetIblankFlag";
   
   if ( !(tpCntl->GetValue(label, &str )) )
   {
   FP3DR.setIBlankFlag(NOT_SET_IBLANK);
   }
   else
   {
   if     ( !strcasecmp(str.c_str(), "on") )  FP3DR.setIBlankFlag(SET_IBLANK);
   else if( !strcasecmp(str.c_str(), "off") ) FP3DR.setIBlankFlag(NOT_SET_IBLANK);
   else
   {
   printf("\tInvalid keyword is described for '%s'\n", label.c_str());
   Exit(0);
   }
   }*/
  FP3DR.setIBlankFlag(NOT_SET_IBLANK); //sphファイルの情報からIblankは作れないので常にoff
  
  /* 次元数
   label = "/Plot3dOptions/Dimension";
   if ( !(tpCntl->GetValue(label, &str )) )
   {
   FP3DR.setDimension3D();
   }
   else
   {
   if     ( !strcasecmp(str.c_str(), "2d") ) FP3DR.setDimension2D();
   else if( !strcasecmp(str.c_str(), "3d") ) FP3DR.setDimension3D();
   else
   {
   printf("\tInvalid keyword is described for '%s'\n", label.c_str());
   Exit(0);
   }
   }*/
  FP3DR.setDimension3D(); // 常に三次元
  
  
  //FormatType
  label = "/Plot3dOptions/FormatType";
  if ( !(tpCntl->GetValue(label, &str )) )
  {
    printf("\tParsing error : fail to get '%s'\n", label.c_str());
    Exit(0);
  }
  else
  {
    if     ( !strcasecmp(str.c_str(), "unformatted") )  FP3DR.setFormat(UNFORMATTED);
    else if( !strcasecmp(str.c_str(), "formatted") )    FP3DR.setFormat(FORMATTED);
    else if( !strcasecmp(str.c_str(), "binary") )       FP3DR.setFormat(C_BINARY);
    else
    {
      printf("\tInvalid keyword is described for '%s'\n", label.c_str());
      Exit(0);
    }
  }
  
  
  // 出力の単精度or倍精度指定
  label = "/Plot3dOptions/RealType";
  if ( !(tpCntl->GetValue(label, &str )) )
  {
    //int d_type = (sizeof(REAL_TYPE) == 4) ? 1 : 2;  // 1-float / 2-double
    //FP3DR.setRealType(d_type);
    FP3DR.setRealType(output_real_type);
  }
  else
  {
    if     ( !strcasecmp(str.c_str(), "float") ) FP3DR.setRealType(1);
    else if( !strcasecmp(str.c_str(), "double") ) FP3DR.setRealType(2);
    else
    {
      printf("\tInvalid keyword is described for '%s'\n", label.c_str());
      Exit(0);
    }
  }
  
  //copy options read to write
  FP3DW.setGridKind(FP3DR.IsGridKind());
  FP3DW.setMoveGrid(FP3DR.IsMoveGrid());
  FP3DW.setSteady(FP3DR.IsSteady());
  FP3DW.setIBlankFlag(FP3DR.IsIBlankFlag());
  FP3DW.setDim(FP3DR.GetDim());
  FP3DW.setFormat(FP3DR.GetFormat());
  FP3DW.setRealType(FP3DR.GetRealType());
  
  // OutputXyz
  label = "/Plot3dOptions/OutputXyz";
  
  if ( !(tpCntl->GetValue(label, &str )) )
  {
    P3Op.IS_xyz = ON;
  }
  else
  {
    if     ( !strcasecmp(str.c_str(), "on") )  P3Op.IS_xyz = ON;
    else if( !strcasecmp(str.c_str(), "off") ) P3Op.IS_xyz = OFF;
    else
    {
      printf("\tInvalid keyword is described for '%s'\n", label.c_str());
      Exit(0);
    }
  }
  
  
  /*
   // OutputQ
   label = "/Plot3dOptions/OutputQ";
   
   if ( !(tpCntl->GetValue(label, &str )) )
   {
   P3Op.IS_q = OFF;
   }
   else
   {
   if     ( !strcasecmp(str.c_str(), "on") )  P3Op.IS_q = ON;
   else if( !strcasecmp(str.c_str(), "off") ) P3Op.IS_q = OFF;
   else
   {
   printf("\tInvalid keyword is described for '%s'\n", label.c_str());
   Exit(0);
   }
   }
   */
  P3Op.IS_q = OFF;
  
  // OutputFunction
  label = "/Plot3dOptions/OutputFunction";
  
  if ( !(tpCntl->GetValue(label, &str )) )
  {
    //printf("\tParsing error : fail to get '%s'\n", label.c_str());
    //Exit(0);
    P3Op.IS_funciton = ON;
  }
  else
  {
    if     ( !strcasecmp(str.c_str(), "on") )  P3Op.IS_funciton = ON;
    else if( !strcasecmp(str.c_str(), "off") ) P3Op.IS_funciton = OFF;
    else
    {
      printf("\tInvalid keyword is described for '%s'\n", label.c_str());
      Exit(0);
    }
  }
  
  
  // OutputFuncName
  label = "/Plot3dOptions/OutputFuncName";
  
  if ( !(tpCntl->GetValue(label, &str )) )
  {
    P3Op.IS_function_name = ON;
  }
  else
  {
    if     ( !strcasecmp(str.c_str(), "on") )  P3Op.IS_function_name = ON;
    else if( !strcasecmp(str.c_str(), "off") ) P3Op.IS_function_name = OFF;
    else
    {
      printf("\tInvalid keyword is described for '%s'\n", label.c_str());
      Exit(0);
    }
  }
  
  /*
   // Output_fvbnd
   label = "/Plot3dOptions/Output_fvbnd";
   
   if ( !(tpCntl->GetValue(label, &str )) )
   {
   P3Op.IS_fvbnd = OFF;
   }
   else
   {
   if     ( !strcasecmp(str.c_str(), "on") )  P3Op.IS_fvbnd = ON;
   else if( !strcasecmp(str.c_str(), "off") ) P3Op.IS_fvbnd = OFF;
   else
   {
   printf("\tInvalid keyword is described for '%s'\n", label.c_str());
   Exit(0);
   }
   }
   */
  P3Op.IS_fvbnd = OFF;
  
  
  // DivideFunc ---> 出力を項目別にファイル分割するオプション
  label = "/Plot3dOptions/DivideFunc";
  
  if ( !(tpCntl->GetValue(label, &str )) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    Exit(0);
  }
  else
  {
    if     ( !strcasecmp(str.c_str(), "on") )  P3Op.IS_DivideFunc = ON;
    else if( !strcasecmp(str.c_str(), "off") ) P3Op.IS_DivideFunc = OFF;
    else
    {
      Hostonly_ stamped_printf("\tInvalid keyword is described for '%s'\n", label.c_str());
      Exit(0);
    }
  }
  
}


// #################################################################
//
void COMB::ReadDfiFiles()
{
  int ic=0;
  vector<string>::const_iterator it;
  
  // allocate dfi info class
  ndfi=dfi_name.size();
  DI = new DfiInfo[ndfi];
  
  // ランク情報をセット
  ic=0;
  for (it = dfi_name.begin(); it != dfi_name.end(); it++) {
    DI[ic].setRankInfo(paraMngr, procGrp);
    ic++;
  }
  
  // set dfi info class
  ic=0;
  for (it = dfi_name.begin(); it != dfi_name.end(); it++) {
    string fname=(*it).c_str();
    DI[ic].ReadDfiFile(fname);
    ic++;
  }
  
  LOG_OUTV_ {
    fprintf(fplog,"\n");
    fprintf(fplog,"*** dfi file info ***\n");
    fprintf(fplog,"\n");
    for(int i=0;i<ndfi;i++){
      fprintf(fplog,"\tDI[%2d].Prefix                      = %s\n",i,DI[i].Prefix.c_str());
      //fprintf(fplog,"\tDI[%2d].RankID_in_MPIworld          = %d\n",i,DI[i].RankID_in_MPIworld);
      //fprintf(fplog,"\tDI[%2d].GroupID_in_MPIworld         = %d\n",i,DI[i].GroupID_in_MPIworld);
      fprintf(fplog,"\tDI[%2d].NumberOfRank                = %d\n",i,DI[i].NumberOfRank);
      fprintf(fplog,"\tDI[%2d].NumberOfGroup               = %d\n",i,DI[i].NumberOfGroup);
      fprintf(fplog,"\tDI[%2d].Global_Voxel[0]             = %d\n",i,DI[i].Global_Voxel[0]);
      fprintf(fplog,"\tDI[%2d].Global_Voxel[1]             = %d\n",i,DI[i].Global_Voxel[1]);
      fprintf(fplog,"\tDI[%2d].Global_Voxel[2]             = %d\n",i,DI[i].Global_Voxel[2]);
      fprintf(fplog,"\tDI[%2d].Global_Division[0]          = %d\n",i,DI[i].Global_Division[0]);
      fprintf(fplog,"\tDI[%2d].Global_Division[1]          = %d\n",i,DI[i].Global_Division[1]);
      fprintf(fplog,"\tDI[%2d].Global_Division[2]          = %d\n",i,DI[i].Global_Division[2]);
      fprintf(fplog,"\tDI[%2d].FileFormat                  = %s\n",i,DI[i].FileFormat.c_str());
      fprintf(fplog,"\tDI[%2d].GuideCell                   = %d\n",i,DI[i].GuideCell);
      fprintf(fplog,"\n");
      fprintf(fplog,"\tDI[%2d].NodeInfoSizet = %d\n",i,DI[i].NodeInfoSize);
      for(int j=0; j< DI[i].NodeInfoSize; j++ ) {
        fprintf(fplog,"\t  DI[%2d].Node[%4d].RankID = %d\n",i,j,DI[i].Node[j].RankID);
        fprintf(fplog,"\t  DI[%2d].Node[%4d].HostName = %s\n",i,j,DI[i].Node[j].HostName.c_str());
        fprintf(fplog,"\t  DI[%2d].Node[%4d].VoxelSize[0] = %d\n",i,j,DI[i].Node[j].VoxelSize[0]);
        fprintf(fplog,"\t  DI[%2d].Node[%4d].VoxelSize[1] = %d\n",i,j,DI[i].Node[j].VoxelSize[1]);
        fprintf(fplog,"\t  DI[%2d].Node[%4d].VoxelSize[2] = %d\n",i,j,DI[i].Node[j].VoxelSize[2]);
        fprintf(fplog,"\t  DI[%2d].Node[%4d].HeadIndex[0] = %d\n",i,j,DI[i].Node[j].HeadIndex[0]);
        fprintf(fplog,"\t  DI[%2d].Node[%4d].HeadIndex[1] = %d\n",i,j,DI[i].Node[j].HeadIndex[1]);
        fprintf(fplog,"\t  DI[%2d].Node[%4d].HeadIndex[2] = %d\n",i,j,DI[i].Node[j].HeadIndex[2]);
        fprintf(fplog,"\t  DI[%2d].Node[%4d].TailIndex[0] = %d\n",i,j,DI[i].Node[j].TailIndex[0]);
        fprintf(fplog,"\t  DI[%2d].Node[%4d].TailIndex[1] = %d\n",i,j,DI[i].Node[j].TailIndex[1]);
        fprintf(fplog,"\t  DI[%2d].Node[%4d].TailIndex[2] = %d\n",i,j,DI[i].Node[j].TailIndex[2]);
        fprintf(fplog,"\t  DI[%2d].Node[%4d].IJK          = %ld\n",i,j,DI[i].Node[j].IJK);
        fprintf(fplog,"\t  DI[%2d].Node[%4d].IJK_JK       = %ld\n",i,j,DI[i].Node[j].IJK_JK);
        fprintf(fplog,"\t  DI[%2d].Node[%4d].IJK_K        = %ld\n",i,j,DI[i].Node[j].IJK_K);
      }
      
      fprintf(fplog,"\n");
      //for(int j=0; j< DI[i].step.size(); j++ ) {
      //  fprintf(fplog,"\t  DI[%2d].step[%4d] = %d\n",i,j,DI[i].step[j]);
      //}
      for(int j=0; j< DI[i].Sc.size(); j++ ) {
        fprintf(fplog,"\t  DI[%2d].step[%4d] = %d\n",i,j,DI[i].Sc[j]->step);
        fprintf(fplog,"\t  DI[%2d].time[%4d] = %f\n",i,j,DI[i].Sc[j]->time);
      }
      
      fprintf(fplog,"\n");
      for(int j=0; j< DI[i].index_y.size(); j++ ) {
        fprintf(fplog,"\t  DI[%2d].index_y[%4d] = %d\n",i,j,DI[i].index_y[j]);
      }
      
      fprintf(fplog,"\n");
      for(int j=0; j< DI[i].index_y.size(); j++ ) {
        fprintf(fplog,"\t  DI[%2d].index_z[%4d] = %d\n",i,j,DI[i].index_z[j]);
      }
    }
    fprintf(fplog,"\n");
  }
  
  STD_OUTV_ {
    printf("\n");
    printf("*** dfi file info ***\n");
    printf("\n");
    for(int i=0;i<ndfi;i++){
      printf("\tDI[%2d].Prefix                      = %s\n",i,DI[i].Prefix.c_str());
      //printf("\tDI[%2d].RankID_in_MPIworld          = %d\n",i,DI[i].RankID_in_MPIworld);
      //printf("\tDI[%2d].GroupID_in_MPIworld         = %d\n",i,DI[i].GroupID_in_MPIworld);
      printf("\tDI[%2d].NumberOfRank                = %d\n",i,DI[i].NumberOfRank);
      printf("\tDI[%2d].NumberOfGroup               = %d\n",i,DI[i].NumberOfGroup);
      printf("\tDI[%2d].Global_Voxel[0]             = %d\n",i,DI[i].Global_Voxel[0]);
      printf("\tDI[%2d].Global_Voxel[1]             = %d\n",i,DI[i].Global_Voxel[1]);
      printf("\tDI[%2d].Global_Voxel[2]             = %d\n",i,DI[i].Global_Voxel[2]);
      printf("\tDI[%2d].Global_Division[0]          = %d\n",i,DI[i].Global_Division[0]);
      printf("\tDI[%2d].Global_Division[1]          = %d\n",i,DI[i].Global_Division[1]);
      printf("\tDI[%2d].Global_Division[2]          = %d\n",i,DI[i].Global_Division[2]);
      printf("\tDI[%2d].FileFormat                  = %s\n",i,DI[i].FileFormat.c_str());
      printf("\tDI[%2d].GuideCell                   = %d\n",i,DI[i].GuideCell);
      printf("\n");
      printf("\tDI[%2d].NodeInfoSizet = %d\n",i,DI[i].NodeInfoSize);
      for(int j=0; j< DI[i].NodeInfoSize; j++ ) {
        printf("\t  DI[%2d].Node[%4d].RankID = %d\n",i,j,DI[i].Node[j].RankID);
        printf("\t  DI[%2d].Node[%4d].HostName = %s\n",i,j,DI[i].Node[j].HostName.c_str());
        printf("\t  DI[%2d].Node[%4d].VoxelSize[0] = %d\n",i,j,DI[i].Node[j].VoxelSize[0]);
        printf("\t  DI[%2d].Node[%4d].VoxelSize[1] = %d\n",i,j,DI[i].Node[j].VoxelSize[1]);
        printf("\t  DI[%2d].Node[%4d].VoxelSize[2] = %d\n",i,j,DI[i].Node[j].VoxelSize[2]);
        printf("\t  DI[%2d].Node[%4d].HeadIndex[0] = %d\n",i,j,DI[i].Node[j].HeadIndex[0]);
        printf("\t  DI[%2d].Node[%4d].HeadIndex[1] = %d\n",i,j,DI[i].Node[j].HeadIndex[1]);
        printf("\t  DI[%2d].Node[%4d].HeadIndex[2] = %d\n",i,j,DI[i].Node[j].HeadIndex[2]);
        printf("\t  DI[%2d].Node[%4d].TailIndex[0] = %d\n",i,j,DI[i].Node[j].TailIndex[0]);
        printf("\t  DI[%2d].Node[%4d].TailIndex[1] = %d\n",i,j,DI[i].Node[j].TailIndex[1]);
        printf("\t  DI[%2d].Node[%4d].TailIndex[2] = %d\n",i,j,DI[i].Node[j].TailIndex[2]);
        printf("\t  DI[%2d].Node[%4d].IJK          = %ld\n",i,j,DI[i].Node[j].IJK);
        printf("\t  DI[%2d].Node[%4d].IJK_JK       = %ld\n",i,j,DI[i].Node[j].IJK_JK);
        printf("\t  DI[%2d].Node[%4d].IJK_K        = %ld\n",i,j,DI[i].Node[j].IJK_K);
      }
      
      printf("\n");
      //for(int j=0; j< DI[i].step.size(); j++ ) {
      //  printf("\t  DI[%2d].step[%4d] = %d\n",i,j,DI[i].step[j]);
      //}
      for(int j=0; j< DI[i].Sc.size(); j++ ) {
        printf("\t  DI[%2d].step[%4d] = %d\n",i,j,DI[i].Sc[j]->step);
        printf("\t  DI[%2d].time[%4d] = %f\n",i,j,DI[i].Sc[j]->time);
      }
      
      printf("\n");
      for(int j=0; j< DI[i].index_y.size(); j++ ) {
        printf("\t  DI[%2d].index_y[%4d] = %d\n",i,j,DI[i].index_y[j]);
      }
      
      printf("\n");
      for(int j=0; j< DI[i].index_y.size(); j++ ) {
        printf("\t  DI[%2d].index_z[%4d] = %d\n",i,j,DI[i].index_z[j]);
      }
    }
    printf("\n");
  }
  
}

// #################################################################
//
void COMB::CheckDir(string dirstr)
{
  Hostonly_
  {
    
#ifndef _WIN32
    
    if( dirstr.size() == 0 ) {
      //printf("\toutput current directory\n");
      return;
    }
    
    DIR* dir;
    if( !(dir = opendir(dirstr.c_str())) ) {
      if( errno == ENOENT ) {
        mode_t mode = S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH;
        if ( !FBUtility::mkdirs(dirstr.c_str()) )
        {
          printf("\tCan't generate directory(%s).\n", dirstr.c_str());
          Exit(0);
        }
      }
      else {
        printf("Directory open error.(%s)", dirstr.c_str());
        Exit(0);
      }
    }
    else {
      if( closedir(dir) == -1 ) {
        printf("Directory close error.(%s)", dirstr.c_str());
        Exit(0);
      }
    }
    
#else // for windows
    
    if( dirstr.size() == 0 ) {
      printf("\toutput current directory\n");
      return;
    }
    
    // check to exist directory
    if (IsDirExsist(dirstr)) {
      // exist directory
      return;
    }
    
    // make directory
    if(!CreateDirectory(dirstr.c_str(), NULL)){
      printf("\tCan't generate directory(%s).\n", dirstr.c_str());
      Exit(0);
    }
    
#endif  // _WIN32
    
  }
  
  return;
}


// #################################################################
//
void COMB::OpenLogFile()
{
  //log file open
  string prefix="log_comb_id";
  int len = prefix.size()+10;//+6+4
  char* tmp = new char[len];
  memset(tmp, 0, sizeof(char)*len);
  sprintf(tmp, "%s_%06d.%s", prefix.c_str(), myRank, "txt");
  
  std::string logname(tmp);
  if ( tmp ) delete [] tmp;
  
  if ( !(fplog = fopen(logname.c_str(), "w")) ){
    printf("\tFile Open Error : '%s'\n",logname.c_str());
    Exit(0);
  }
  fprintf(fplog,"####################\n",logname.c_str());
  fprintf(fplog,"### log_comb.txt ###\n",logname.c_str());
  fprintf(fplog,"####################\n",logname.c_str());
  fprintf(fplog,"\n");
  
  fprintf(fplog,"procGrp  = %d\n", procGrp);
  fprintf(fplog,"myRank   = %d\n", myRank);
  fprintf(fplog,"numProc  = %d\n", numProc);
  fprintf(fplog,"HostName = %s\n", HostName.c_str());
  fprintf(fplog,"\n");
  
}

// #################################################################
//
void COMB::CloseLogFile()
{
  fclose(fplog);
}

// #################################################################
//
void COMB::WriteTime(double* tt)
{
  fprintf(fplog,"\n\n");
  fprintf(fplog,"TIME : ReadInit      %10.3f sec.\n", tt[0]);
  fprintf(fplog,"TIME : ReadDfiFiles  %10.3f sec.\n", tt[1]);
  fprintf(fplog,"TIME : CombineFiles  %10.3f sec.\n", tt[2]);
  fprintf(fplog,"TIME : Total Time    %10.3f sec.\n", tt[3]);
}

// #################################################################
//
void COMB::CombineFiles()
{
  if( out_format == OUTFORMAT_IS_SPH )
  {
    output_sph();
  }
  else if( out_format == OUTFORMAT_IS_PLOT3D )
  {
    output_plot3d();
  }
}

// #################################################################
// 出力DFIファイル名を作成する
std::string COMB::Generate_DFI_Name(const std::string prefix)
{
  if ( prefix.empty() ) return NULL;
  
  int len = prefix.size() + 5; // postfix(4) + 1
  char* tmp = new char[len];
  memset(tmp, 0, sizeof(char)*len);
  
  sprintf(tmp, "%s.%s", prefix.c_str(), "dfi");
  
  std::string fname(tmp);
  if ( tmp ) delete [] tmp;
  
  return fname;
}


// #################################################################
// ファイル名を作成する
std::string COMB::Generate_FileName(const std::string prefix, const unsigned m_step, const int m_id, const bool mio)
{
  if ( prefix.empty() ) return NULL;
  
  int len = prefix.size() + 24; // step(10) + id(9) + postfix(4) + 1
  char* tmp = new char[len];
  memset(tmp, 0, sizeof(char)*len);
  
  // local出力が指定された場合、分割出力
  if ( mio )
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
std::string COMB::Generate_FileName_Free(const std::string prefix, const std::string xxx, const unsigned m_step, const int m_id, const bool mio)
{
  if ( prefix.empty() ) return NULL;
  
  int len = prefix.size() + xxx.size() + 21; // step(10) + id(9) + 1(.拡張子) + 1(\0)
  char* tmp = new char[len];
  memset(tmp, 0, sizeof(char)*len);
  
  // local出力が指定された場合、分割出力
  if ( mio )
  {
    sprintf(tmp, "%s_%06d_%010d.%s", prefix.c_str(), m_id, m_step, xxx.c_str());
    //sprintf(tmp, "%s_%010d_id%06d.%s", prefix.c_str(), m_step, m_id, xxx.c_str());
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
// DFIファイルをコピーする
bool COMB::Copy_DFIFile(const std::string base_name, const std::string new_name, const std::string prefix)//, int& dfi_mng)
{
  if ( base_name.empty() ) return false;
  if ( new_name.empty() ) return false;
  if ( prefix.empty() ) return false;
  
  FILE* fpin;
  FILE* fpout;
  char buff[DFI_LINE_LENGTH];
  string buffs;
  
  if ( !(fpin = fopen(base_name.c_str(), "r")) ) return false;
  if ( !(fpout = fopen(new_name.c_str(), "w")) ) return false;
  
  while( fgets(buff,DFI_LINE_LENGTH,fpin) != NULL ) {
    buffs=buff;
    //if( !strcasecmp(buffs.c_str(), "  FileInfo {\n" ) ) break;
    //if( !strcasecmp(buffs.substr(0,12).c_str(), "  FileInfo {") ) break;
    if( !strcasecmp(buffs.substr(0,8).c_str(), "  Prefix") ){
      Write_Tab(fpout, 1);
      fprintf(fpout, "Prefix        = \"%s\"\n", prefix.c_str());
    }
    else if( !strcasecmp(buffs.substr(0,9).c_str(), "      Min") ){
      Write_Tab(fpout, 3);
      fprintf(fpout, "Min  = %e\n", 0.0);
    }
    else if( !strcasecmp(buffs.substr(0,9).c_str(), "      Max") ){
      Write_Tab(fpout, 3);
      fprintf(fpout, "Max  = %e\n", 0.0);
    }
    else{
      fputs(buff,fpout);
    }
  }
  
  //if (fpout) Write_Tab(fpout, 1);
  //if (fpout) fprintf(fpout, "FileInfo {\n");
  //if (fpout) Write_Tab(fpout, 1);
  //if (fpout) fprintf(fpout, "}\n");
  //if (fpout) fprintf(fpout, "}\n");
  
  fclose(fpin);
  fclose(fpout);
  
  //dfi_mng++;
  
  return true;
}


// #################################################################
// Tab(space２つ)を出力する
void COMB::Write_Tab(FILE* fp, const unsigned tab)
{
  for(int n=0; n<tab; n++) fprintf(fp, "  ");
}

// #################################################################
// メモリ使用量を表示する
void COMB::MemoryRequirement(const double Memory, FILE* fp)
{
  const double mem = Memory;
  const double KB = 1024.0;
  const double MB = 1024.0*KB;
  const double GB = 1024.0*MB;
  const double TB = 1024.0*GB;
  const double PB = 1024.0*TB;
  const double factor = 1.05; // estimate 5% for addtional
  
  // Global memory
  fprintf (fp," MemorySize = ");
  if ( mem > PB ) {
    fprintf (fp,"%6.2f (PB)\n", mem / PB *factor);
  }
  else if ( mem > TB ) {
    fprintf (fp,"%6.2f (TB)\n", mem / TB *factor);
  }
  else if ( mem > GB ) {
    fprintf (fp,"%6.2f (GB)\n", mem / GB *factor);
  }
  else if ( mem > MB ) {
    fprintf (fp,"%6.2f (MB)\n", mem / MB *factor);
  }
  else if ( mem > KB ) {
    fprintf (fp,"%6.2f (KB)\n", mem / KB *factor);
  }
  else if ( mem <= KB ){
    fprintf (fp,"%6.2f (B)\n", mem *factor);
  }
  else {
    fprintf (fp,"Caution! Memory required : %d (Byte)", (int)(mem *factor) );
  }
  
  fflush(fp);
}

// #################################################################
// メモリ使用量を表示する
void COMB::MemoryRequirement(const double TotalMemory, const double sphMemory, const double plot3dMemory, const double thinMemory, FILE* fp)
{
  double mem;
  const double KB = 1024.0;
  const double MB = 1024.0*KB;
  const double GB = 1024.0*MB;
  const double TB = 1024.0*GB;
  const double PB = 1024.0*TB;
  const double factor = 1.05; // estimate 5% for addtional
  
  fprintf (fp,"*** Required MemorySize ***");
  fprintf (fp,"\n");
  
  mem = sphMemory;
  fprintf (fp,"  read SPH MemorySize = ");
  if ( mem > PB ) {
    fprintf (fp,"%6.2f (PB)", mem / PB *factor);
  }
  else if ( mem > TB ) {
    fprintf (fp,"%6.2f (TB)", mem / TB *factor);
  }
  else if ( mem > GB ) {
    fprintf (fp,"%6.2f (GB)", mem / GB *factor);
  }
  else if ( mem > MB ) {
    fprintf (fp,"%6.2f (MB)", mem / MB *factor);
  }
  else if ( mem > KB ) {
    fprintf (fp,"%6.2f (KB)", mem / KB *factor);
  }
  else if ( mem <= KB ){
    fprintf (fp,"%6.2f (B)", mem *factor);
  }
  else {
    fprintf (fp,"Caution! Memory required : %d (Byte)", (int)(mem *factor) );
  }
  
  mem = plot3dMemory;
  fprintf (fp,"  write PLOT3D MemorySize = ");
  if ( mem > PB ) {
    fprintf (fp,"%6.2f (PB)", mem / PB *factor);
  }
  else if ( mem > TB ) {
    fprintf (fp,"%6.2f (TB)", mem / TB *factor);
  }
  else if ( mem > GB ) {
    fprintf (fp,"%6.2f (GB)", mem / GB *factor);
  }
  else if ( mem > MB ) {
    fprintf (fp,"%6.2f (MB)", mem / MB *factor);
  }
  else if ( mem > KB ) {
    fprintf (fp,"%6.2f (KB)", mem / KB *factor);
  }
  else if ( mem <= KB ){
    fprintf (fp,"%6.2f (B)", mem *factor);
  }
  else {
    fprintf (fp,"Caution! Memory required : %d (Byte)", (int)(mem *factor) );
  }
  
  mem = thinMemory;
  fprintf (fp,"  write thin out MemorySize = ");
  if ( mem > PB ) {
    fprintf (fp,"%6.2f (PB)", mem / PB *factor);
  }
  else if ( mem > TB ) {
    fprintf (fp,"%6.2f (TB)", mem / TB *factor);
  }
  else if ( mem > GB ) {
    fprintf (fp,"%6.2f (GB)", mem / GB *factor);
  }
  else if ( mem > MB ) {
    fprintf (fp,"%6.2f (MB)", mem / MB *factor);
  }
  else if ( mem > KB ) {
    fprintf (fp,"%6.2f (KB)", mem / KB *factor);
  }
  else if ( mem <= KB ){
    fprintf (fp,"%6.2f (B)", mem *factor);
  }
  else {
    fprintf (fp,"Caution! Memory required : %d (Byte)", (int)(mem *factor) );
  }
  
  mem = TotalMemory;
  fprintf (fp,"  TotalMemorySize = ");
  if ( mem > PB ) {
    fprintf (fp,"%6.2f (PB)", mem / PB *factor);
  }
  else if ( mem > TB ) {
    fprintf (fp,"%6.2f (TB)", mem / TB *factor);
  }
  else if ( mem > GB ) {
    fprintf (fp,"%6.2f (GB)", mem / GB *factor);
  }
  else if ( mem > MB ) {
    fprintf (fp,"%6.2f (MB)", mem / MB *factor);
  }
  else if ( mem > KB ) {
    fprintf (fp,"%6.2f (KB)", mem / KB *factor);
  }
  else if ( mem <= KB ){
    fprintf (fp,"%6.2f (B)", mem *factor);
  }
  else {
    fprintf (fp,"Caution! Memory required : %d (Byte)", (int)(mem *factor) );
  }
  
  fprintf (fp,"\n");
  fprintf (fp,"\n");
  
  fflush(fp);
}


// #################################################################
//
int COMB::tt_check_machine_endian()
{
	int v = 1;
	char* p = (char*)&v;
  
	if (p[0])					return TT_LITTLE_ENDIAN;
	else if (p[sizeof(int)-1])	return TT_BIG_ENDIAN;
	else						return TT_OTHER_ENDIAN;
}


// #################################################################
//
bool COMB::ReadSphDataType(int* m_sv_type, int* m_d_type, int fp_in, string fname)
{
  if( !(FileIO_SPH::GetDataType(fname.c_str(), m_d_type, m_sv_type)) ){
    printf("\terror : FileIO_SPH::GetDataType\n");
    Exit(0);
  }
  return true;
}

// #################################################################
//
bool COMB::ReadSphHeader(int* m_step,
                         int* m_sv_type,
                         int* m_d_type,
                         int* m_imax,
                         int* m_jmax,
                         int* m_kmax,
                         float* m_time,
                         float* m_org,
                         float* m_pit,
                         int fp_in,
                         string fname)
{
  unsigned int voxsize[3];
  if( !(FileIO_SPH::GetHeader(fname.c_str(), m_sv_type, voxsize, m_org, m_pit, m_step, m_time)) ){
    printf("\terror : FileIO_SPH::GetHeader\n");
    Exit(0);
  }
  *m_imax=(int)voxsize[0];
  *m_jmax=(int)voxsize[1];
  *m_kmax=(int)voxsize[2];
  return true;
}

// #################################################################
//
bool COMB::ReadSphHeader(long long* m_step,
                         int* m_sv_type,
                         int* m_d_type,
                         long long* m_imax,
                         long long* m_jmax,
                         long long* m_kmax,
                         double* m_time,
                         double* m_org,
                         double* m_pit,
                         int fp_in,
                         string fname)
{
  unsigned long long voxsize[3];
  if( !(FileIO_SPH::GetHeader(fname.c_str(), m_sv_type, voxsize, m_org, m_pit, m_step, m_time)) ){
    printf("\terror : FileIO_SPH::GetHeader\n");
    Exit(0);
  }
  *m_imax=(long long)voxsize[0];
  *m_jmax=(long long)voxsize[1];
  *m_kmax=(long long)voxsize[2];
  return true;
}


// #################################################################
//
bool COMB::ReadSphData(float* wk,
                       int wksize,
                       int* size,
                       int dim,
                       int fp_in,
                       string fname)
{
  int m_sv_type;
  int m_d_type;
  int m_step;
  float m_time;
  float m_org[3];
  float m_pit[3];
  unsigned int voxsize[3];
  unsigned int dsize=(unsigned int)wksize;
  if( !(FileIO_SPH::GetData(fname.c_str(), dsize, &m_sv_type, voxsize, m_org, m_pit, &m_step, &m_time, wk)) ){
    printf("\terror : FileIO_SPH::GetData\n");
    Exit(0);
  }
  return true;
}

// #################################################################
//
bool COMB::ReadSphData(double* wk,
                       int wksize,
                       int* size,
                       int dim,
                       int fp_in,
                       string fname)
{
  int m_sv_type;
  int m_d_type;
  long long m_step;
  double m_time;
  double m_org[3];
  double m_pit[3];
  unsigned long long voxsize[3];
  unsigned long long dsize=(unsigned long long)wksize;
  if( !(FileIO_SPH::GetData(fname.c_str(), dsize, &m_sv_type, voxsize, m_org, m_pit, &m_step, &m_time, wk)) ){
    printf("\terror : FileIO_SPH::GetData\n");
    Exit(0);
  }
  return true;
}
