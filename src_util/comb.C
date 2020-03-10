//##################################################################################
//
// FFV-C : Frontflow / violet Cartesian
//
// Copyright (c) 2007-2011 VCAD System Research Program, RIKEN.
// All rights reserved.
//
// Copyright (c) 2011-2015 Institute of Industrial Science, The University of Tokyo.
// All rights reserved.
//
// Copyright (c) 2012-2016 Advanced Institute for Computational Science, RIKEN.
// All rights reserved.
//
//##################################################################################

/**
 * @file   comb.C
 * @brief  COMB Class
 * @author aics
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
  
  /* PLOT3Dfunctoins_20131005
  P3Op.IS_xyz = ON;
  P3Op.IS_q = OFF;
  P3Op.IS_funciton = ON;
  P3Op.IS_function_name = ON;
  P3Op.IS_fvbnd = OFF;
  P3Op.IS_DivideFunc=OFF;
  P3Op.ngrid=0; //出力ブロック数
  P3Op.nvar=0;  //出力項目数
  */
  
  out_format=0;
  
  output_real_type=0;

  ndfi=0;
  dfi.clear();
  
  // dfi管理
  for (int i=0; i<var_END; i++) dfi_mng[i]=0;
  
  staging=0;
  
  output_format_type = 0;
  output_divfunc     = 0;
  
}


// #################################################################
// デストラクタ
COMB::~COMB()
{
  for(int i=0; i<dfi.size(); i++ ) if( !dfi[i] ) delete dfi[i];
  
}

// #################################################################
//
void COMB::ReadInit(string input_file)
{
  
  // ------------------------------------
  FILE* fp = NULL;
  
  // TPインスタンス生成
  TextParser tpCntl;
  
  // 入力ファイルをセット
  int ierror = tpCntl.read(input_file);
  
  // 入力ファイルの読み込み--->パラメータのセット
  ReadInputFile(&tpCntl);
  
  // TextParserの破棄
  tpCntl.remove();
  
  return;
}

// #################################################################
//
void COMB::ReadInputFile(TextParser* tpCntl)
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
    
    if ( !(tpCntl->getNodeStr(label_base, i+1, str)) )
    {
      printf("\tParsing error : No Elem name\n");
      Exit(0);
    }
    if( strcasecmp(str.substr(0,4).c_str(), "list") ) continue;
    label=label_base+"/"+str;
    
    if ( !(tpCntl->getInspectedValue(label, str )) ) {
      printf("\tParsing error : fail to get '%s'\n", label.c_str());
      Exit(0);
    }

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
  if ( (tpCntl->getInspectedValue(label, str )) )
  {
    out_dirname=str;
    LOG_OUT_ fprintf(fplog,"\tReset Output Directory '%s'\n", out_dirname.c_str());
    STD_OUT_ printf("\tReset Output Directory '%s'\n", out_dirname.c_str());
    CheckDir(out_dirname);
    if( out_dirname.size() != 0 ) out_dirname=out_dirname+"/";
  }
  
  
  //入力ディレクトリの指定
  label = "/CombData/InputDir";
  if ( (tpCntl->getInspectedValue(label, str )) )
  {
    in_dirname=str;
    LOG_OUT_ fprintf(fplog,"\tReset Input Directory '%s'\n", in_dirname.c_str());
    STD_OUT_ printf("\tReset Input Directory '%s'\n", in_dirname.c_str());
    CheckDir(in_dirname);
    if( in_dirname.size() != 0 ) in_dirname=in_dirname+"/";
  }
  
  
  //並列実行時のSTAGINGのON/OFF
  label = "/CombData/Staging";
  if ( !(tpCntl->getInspectedValue(label, str )) )
  {
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
  if ( !(tpCntl->getInspectedValue(label, str )) )
  {
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
  if ( !(tpCntl->getInspectedValue(label, str )) )
  {
    printf("\tParsing error : fail to get '%s'\n", label.c_str());
    Exit(0);
  }
  if     ( !strcasecmp(str.c_str(), "sph" ) )    out_format = OUTFORMAT_IS_SPH;
  else if( !strcasecmp(str.c_str(), "plot3d" ) ) out_format = OUTFORMAT_IS_PLOT3D;
  else if( !strcasecmp(str.c_str(), "avs" ) ) out_format = OUTFORMAT_IS_AVS;
  else
  {
    printf("\tInvalid keyword is described for  '%s'\n", label.c_str());
    Exit(0);
  }
  
  // PLOT3Dオプションの読み込み
  // PLOT3Dfunctions_20131005 if( out_format == OUTFORMAT_IS_PLOT3D ) get_PLOT3D(tpCntl);
  
  // AVSオプションの読み込み
  if( out_format == OUTFORMAT_IS_AVS ) get_AVSoptions(tpCntl);
  
}

/* PLOT3Dfunctions_20131005
// #################################################################
// PLOT3Dファイル入出力に関するパラメータ
void COMB::get_PLOT3D(TextParser* tpCntl)
{
  string str;
  string label;
  
  // Filename
  //label = "/Plot3dOptions/Filename";
  //
  //if ( !(tpCntl->getInspectedValue(label, &str )) )
  //{
  //  P3Op.basename = "PLOT3Doutput";
  //}
  //else
  //{
  //  P3Op.basename = str;
  //}
  
  // FileNameGrid --- option
  label = "/Steer/Plot3dOptions/FileNameGrid";
  
  if ( !(tpCntl->getInspectedValue(label, str)) )
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
  
  if ( !(tpCntl->getInspectedValue(label, str)) )
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
   
   if ( !(tpCntl->getInspectedValue(label, &str )) )
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
   
  FP3DR.setMultiGrid(); // 常にmulti grid
  
  // 格子の移動
  label = "/Plot3dOptions/GridMobility";
  
  if ( !(tpCntl->getInspectedValue(label, str )) )
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
  
  if ( !(tpCntl->getInspectedValue(label, str )) )
  {
    FP3DR.setSteady(UNSTEADY);
  }
  else
  {
    if     ( !strcasecmp(str.c_str(), "steady") )   FP3DR.setSteady(STEADY);
    else if( !strcasecmp(str.c_str(), "unsteady") ) FP3DR.setSteady(UNSTEADY);
    else
    {
      printf("\tInvalid keyword is described for '%s'\n", label.c_str());
      Exit(0);
    }
  }
  
  
  /* IBLANKファイル
   label = "/Plot3dOptions/SetIblankFlag";
   
   if ( !(tpCntl->getInspectedValue(label, &str )) )
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
   }
  FP3DR.setIBlankFlag(NOT_SET_IBLANK); //sphファイルの情報からIblankは作れないので常にoff
  
  /* 次元数
   label = "/Plot3dOptions/Dimension";
   if ( !(tpCntl->getInspectedValue(label, &str )) )
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
   }
  FP3DR.setDimension3D(); // 常に三次元
  
  
  //FormatType
  label = "/Plot3dOptions/FormatType";
  if ( !(tpCntl->getInspectedValue(label, str )) )
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
  if ( !(tpCntl->getInspectedValue(label, str )) )
  {
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
  
  if ( !(tpCntl->getInspectedValue(label, str )) )
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
   
   if ( !(tpCntl->getInspectedValue(label, &str )) )
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
   
  P3Op.IS_q = OFF;
  
  // OutputFunction
  label = "/Plot3dOptions/OutputFunction";
  
  if ( !(tpCntl->getInspectedValue(label, str )) )
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
  
  if ( !(tpCntl->getInspectedValue(label, str )) )
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
   
   if ( !(tpCntl->getInspectedValue(label, &str )) )
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
   
  P3Op.IS_fvbnd = OFF;
  
  
  // DivideFunc ---> 出力を項目別にファイル分割するオプション
  label = "/Plot3dOptions/DivideFunc";
  
  if ( !(tpCntl->getInspectedValue(label, str )) )
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
 */

// #################################################################
// AVSファイル入出力に関するパラメータ
void COMB::get_AVSoptions(TextParser* tpCntl)
{
  string str;
  string label;
  
  //FormatType
  label = "/AVSoptions/FormatType";
  if ( !(tpCntl->getInspectedValue(label, str )) )
  {
    output_format_type = OUTPUT_FORMAT_TYPE_BINARY;
  }
  else
  {
    if     ( !strcasecmp(str.c_str(), "binary") ) output_format_type = OUTPUT_FORMAT_TYPE_BINARY;
    else if( !strcasecmp(str.c_str(), "ascii") )  output_format_type = OUTPUT_FORMAT_TYPE_ASCII;
    else
    {
      Hostonly_ stamped_printf("\tInvalid keyword is described for '%s'\n", label.c_str());
      Exit(0);
    }
  }
  
  //DivideFunc
  label = "/AVSoptions/DivideFunc";
  if(  !(tpCntl->getInspectedValue(label, str )) )
  {
    output_divfunc = OUTPUT_DIV_FUNC_OFF;
  }
  else
  {
    if     ( !strcasecmp(str.c_str(), "on") )  output_divfunc = OUTPUT_DIV_FUNC_ON;
    else if( !strcasecmp(str.c_str(), "uff") ) output_divfunc = OUTPUT_DIV_FUNC_OFF;
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
  
  
  // dfiファイルの読込み
  ic=0;
  int tempg[3];
  int tempd[3];
  CDM::E_CDM_ERRORCODE ret = CDM::E_CDM_SUCCESS;
  for (it = dfi_name.begin(); it != dfi_name.end(); it++) {
    string fname=(*it).c_str();
    cdm_DFI* dfi_in = cdm_DFI::ReadInit(paraMngr->GetMPI_Comm(procGrp),
                                        fname,
                                        tempg,
                                        tempd,
                                        ret);
    if( dfi_in == NULL ) exit(0);
    if( ret != CDM::E_CDM_SUCCESS && ret != CDM::E_CDM_ERROR_INVALID_DIVNUM ) exit(0);
    dfi.push_back(dfi_in);
    ic++;
  }
  
  LOG_OUTV_ {
    fprintf(fplog,"\n");
    fprintf(fplog,"*** dfi file info ***\n");
    fprintf(fplog,"\n");
    //for(int i=0;i<ndfi;i++){
    for(int i=0;i<dfi.size();i++){
      const cdm_FileInfo *DFI_Info = dfi[i]->GetcdmFileInfo();
      fprintf(fplog,"\tDFI_Info->DirectoryPath            = %s\n",DFI_Info->DirectoryPath.c_str());
      fprintf(fplog,"\tDFI_Info->TimeSliceDirFlag         = %d\n",DFI_Info->TimeSliceDirFlag);
      fprintf(fplog,"\tDFI_Info->Prefix                   = %s\n",DFI_Info->Prefix.c_str());
      fprintf(fplog,"\tDFI_Info->FileFormat               = %d\n",DFI_Info->FileFormat);
      fprintf(fplog,"\tDFI_Info->GuideCell                = %d\n",DFI_Info->GuideCell);
      fprintf(fplog,"\tDFI_Info->DataType                 = %d\n",DFI_Info->DataType);
      fprintf(fplog,"\tDFI_Info->Endian                   = %d\n",DFI_Info->Endian);
      fprintf(fplog,"\tDFI_Info->ArrayShape               = %d\n",DFI_Info->ArrayShape);
      fprintf(fplog,"\tDFI_Info->NumVariables             = %d\n",DFI_Info->NumVariables);
      
      const cdm_MPI *DFI_MPI = dfi[i]->GetcdmMPI();
      fprintf(fplog,"\tDFI_MPI->NumberOfRank              = %d\n",DFI_MPI->NumberOfRank);
      fprintf(fplog,"\tDFI_MPI->NumberOfGroup             = %d\n",DFI_MPI->NumberOfGroup);
      
      const cdm_Domain *DFI_Domain = dfi[i]->GetcdmDomain();
      fprintf(fplog,"\tDFI_Domain->GlobalVoxel[0]         = %d\n",DFI_Domain->GlobalVoxel[0]);
      fprintf(fplog,"\tDFI_Domain->GlobalVoxel[1]         = %d\n",DFI_Domain->GlobalVoxel[1]);
      fprintf(fplog,"\tDFI_Domain->GlobalVoxel[2]         = %d\n",DFI_Domain->GlobalVoxel[2]);
      fprintf(fplog,"\tDFI_Domain->GlobalDivision[0]      = %d\n",DFI_Domain->GlobalDivision[0]);
      fprintf(fplog,"\tDFI_Domain->GlobalDivision[1]      = %d\n",DFI_Domain->GlobalDivision[1]);
      fprintf(fplog,"\tDFI_Domain->GlobalDivision[2]      = %d\n",DFI_Domain->GlobalDivision[2]);
      
      const cdm_Process *DFI_Process = dfi[i]->GetcdmProcess();
      fprintf(fplog,"\n");
      fprintf(fplog,"\tDFI_Process->RankList.size()       = %d\n",DFI_Process->RankList.size());
      for(int j=0; j< DFI_Process->RankList.size(); j++ ) {
        fprintf(fplog,"\t  DFI_Process->RankList[%d].RankID       = %d\n",j,DFI_Process->RankList[j].RankID);
        fprintf(fplog,"\t  DFI_Process->RankList[%d].HostName     = %s\n",j,DFI_Process->RankList[j].HostName.c_str());
        fprintf(fplog,"\t  DFI_Process->RankList[%d].VoxelSize[0] = %d\n",j,DFI_Process->RankList[j].VoxelSize[0]);
        fprintf(fplog,"\t  DFI_Process->RankList[%d].VoxelSize[1] = %d\n",j,DFI_Process->RankList[j].VoxelSize[1]);
        fprintf(fplog,"\t  DFI_Process->RankList[%d].VoxelSize[2] = %d\n",j,DFI_Process->RankList[j].VoxelSize[2]);
        fprintf(fplog,"\t  DFI_Process->RankList[%d].HeadIndex[0] = %d\n",j,DFI_Process->RankList[j].HeadIndex[0]);
        fprintf(fplog,"\t  DFI_Process->RankList[%d].HeadIndex[1] = %d\n",j,DFI_Process->RankList[j].HeadIndex[1]);
        fprintf(fplog,"\t  DFI_Process->RankList[%d].HeadIndex[2] = %d\n",j,DFI_Process->RankList[j].HeadIndex[2]);
        fprintf(fplog,"\t  DFI_Process->RankList[%d].TailIndex[0] = %d\n",j,DFI_Process->RankList[j].TailIndex[0]);
        fprintf(fplog,"\t  DFI_Process->RankList[%d].TailIndex[1] = %d\n",j,DFI_Process->RankList[j].TailIndex[1]);
        fprintf(fplog,"\t  DFI_Process->RankList[%d].TailIndex[2] = %d\n",j,DFI_Process->RankList[j].TailIndex[2]);
      }
      
      fprintf(fplog,"\n");
      const cdm_TimeSlice* DFI_TSlice = dfi[i]->GetcdmTimeSlice();
      for(int j=0; j< DFI_TSlice->SliceList.size(); j++ ) {
        fprintf(fplog,"\t  DFI_TSlice->SliceList[%d].step         = %d\n",j,DFI_TSlice->SliceList[j].step);
        fprintf(fplog,"\t  DFI_TSlice->SliceList[%d}.time         = %f\n",j,DFI_TSlice->SliceList[j].time);
      }
    }
    fprintf(fplog,"\n");
  }
  
  STD_OUTV_ {
    printf("\n");
    printf("*** dfi file info ***\n");
    printf("\n");
    for(int i=0;i<dfi.size();i++){
      const cdm_FileInfo *DFI_Info = dfi[i]->GetcdmFileInfo();
      printf("\tDFI_Info->DirectoryPath            = %s\n",DFI_Info->DirectoryPath.c_str());
      printf("\tDFI_Info->TimeSliceDirFlag         = %d\n",DFI_Info->TimeSliceDirFlag);
      printf("\tDFI_Info->Prefix                   = %s\n",DFI_Info->Prefix.c_str());
      printf("\tDFI_Info->FileFormat               = %d\n",DFI_Info->FileFormat);
      printf("\tDFI_Info->GuideCell                = %d\n",DFI_Info->GuideCell);
      printf("\tDFI_Info->DataType                 = %d\n",DFI_Info->DataType);
      printf("\tDFI_Info->Endian                   = %d\n",DFI_Info->Endian);
      printf("\tDFI_Info->ArrayShape               = %d\n",DFI_Info->ArrayShape);
      printf("\tDFI_Info->NumVariables             = %d\n",DFI_Info->NumVariables);
      
      const cdm_MPI *DFI_MPI = dfi[i]->GetcdmMPI();
      printf("\tDFI_MPI->NumberOfRank              = %d\n",DFI_MPI->NumberOfRank);
      printf("\tDFI_MPI->NumberOfGroup             = %d\n",DFI_MPI->NumberOfGroup);
      
      const cdm_Domain *DFI_Domain = dfi[i]->GetcdmDomain();
      printf("\tDFI_Domain->GlobalVoxel[0]         = %d\n",DFI_Domain->GlobalVoxel[0]);
      printf("\tDFI_Domain->GlobalVoxel[1]         = %d\n",DFI_Domain->GlobalVoxel[1]);
      printf("\tDFI_Domain->GlobalVoxel[2]         = %d\n",DFI_Domain->GlobalVoxel[2]);
      printf("\tDFI_Domain->GlobalDivision[0]      = %d\n",DFI_Domain->GlobalDivision[0]);
      printf("\tDFI_Domain->GlobalDivision[1]      = %d\n",DFI_Domain->GlobalDivision[1]);
      printf("\tDFI_Domain->GlobalDivision[2]      = %d\n",DFI_Domain->GlobalDivision[2]);
      
      const cdm_Process *DFI_Process = dfi[i]->GetcdmProcess();
      printf("\n");
      printf("\tDFI_Process->RankList.size()       = %d\n",(int)DFI_Process->RankList.size());
      for(int j=0; j< DFI_Process->RankList.size(); j++ ) {
        printf("\t  DFI_Process->RankList[%d].RankID       = %d\n",j,DFI_Process->RankList[j].RankID);
        printf("\t  DFI_Process->RankList[%d].HostName     = %s\n",j,DFI_Process->RankList[j].HostName.c_str());
        printf("\t  DFI_Process->RankList[%d].VoxelSize[0] = %d\n",j,DFI_Process->RankList[j].VoxelSize[0]);
        printf("\t  DFI_Process->RankList[%d].VoxelSize[1] = %d\n",j,DFI_Process->RankList[j].VoxelSize[1]);
        printf("\t  DFI_Process->RankList[%d].VoxelSize[2] = %d\n",j,DFI_Process->RankList[j].VoxelSize[2]);
        printf("\t  DFI_Process->RankList[%d].HeadIndex[0] = %d\n",j,DFI_Process->RankList[j].HeadIndex[0]);
        printf("\t  DFI_Process->RankList[%d].HeadIndex[1] = %d\n",j,DFI_Process->RankList[j].HeadIndex[1]);
        printf("\t  DFI_Process->RankList[%d].HeadIndex[2] = %d\n",j,DFI_Process->RankList[j].HeadIndex[2]);
        printf("\t  DFI_Process->RankList[%d].TailIndex[0] = %d\n",j,DFI_Process->RankList[j].TailIndex[0]);
        printf("\t  DFI_Process->RankList[%d].TailIndex[1] = %d\n",j,DFI_Process->RankList[j].TailIndex[1]);
        printf("\t  DFI_Process->RankList[%d].TailIndex[2] = %d\n",j,DFI_Process->RankList[j].TailIndex[2]);
      }
      
      printf("\n");
      const cdm_TimeSlice* DFI_TSlice = dfi[i]->GetcdmTimeSlice();
      for(int j=0; j< DFI_TSlice->SliceList.size(); j++ ) {
        printf("\t  DFI_TSlice->SliceList[%d].step         = %d\n",j,DFI_TSlice->SliceList[j].step);
        printf("\t  DFI_TSlice->SliceList[%d}.time         = %f\n",j,DFI_TSlice->SliceList[j].time);
      }
      
    }
  }
  
}

// #################################################################
//
void COMB::CheckDir(string dirstr)
{
  Hostonly_
  {
    
    if( dirstr.size() == 0 ) {
      //printf("\toutput current directory\n");
      return;
    }
    
    DIR* dir;
    if( !(dir = opendir(dirstr.c_str())) ) {
      if( errno == ENOENT ) {
        mode_t mode = S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH;
        //if ( !FBUtility::mkdirs(dirstr.c_str()) )
//CDM.20131008.s
        //if ( !cdm_DFI::MakeDirectorySub(dirstr) )
        if ( cdm_DFI::MakeDirectorySub(dirstr) != 0 )
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
  fprintf(fplog,"####################\n");
  fprintf(fplog,"### log_comb.txt ###\n");
  fprintf(fplog,"####################\n");
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
    // PLTT3Dfunctions_20131005 output_plot3d();
  }
  else if( out_format == OUTFORMAT_IS_AVS )
  {
//CDM.20131008.s
    //output_avs();
    output_sph();
    output_avs_header();
//CDM.20131008.e
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
//CDM.20131008.s
   if( !strcasecmp(xxx.c_str(), "sph") || !strcasecmp(xxx.c_str(), "dat") )
   { 
    sprintf(tmp, "%s_%010d.%s", prefix.c_str(), m_step, xxx.c_str());
   }
    else
    {
      sprintf(tmp, "%s.%s", prefix.c_str(), xxx.c_str());
    }
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
  
  fclose(fpin);
  fclose(fpout);
  
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
bool COMB::ReadSphHeader(
                         double* m_org,
                         double* m_pit,
                         string fname)
{
  if( !(FileIO_SPH::GetHeader(fname.c_str(), m_org, m_pit)) ){
    printf("\terror : FileIO_SPH::GetHeader\n");
    Exit(0);
  }
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

// #################################################################
//
bool COMB::read_HeaderRecord(FILE* fp,
                             EMatchType eType,
                             const int m_d_type)
{
  
  unsigned int dmy, type_dmy;
  int data_type,real_type;
  
  if( fread(&dmy, sizeof(int), 1, fp) != 1 ) { fclose(fp); return false; }
  if( eType == UnMatch ) BSWAP32(dmy);
  if( dmy != 8 ) { fclose(fp); return false; }
  if( fread(&data_type, sizeof(int), 1, fp) != 1 ) { fclose(fp); return false; }
  if( eType == UnMatch ) BSWAP32(data_type);
  if( fread(&real_type, sizeof(int), 1, fp) != 1 ) { fclose(fp); return false; }
  if( eType == UnMatch ) BSWAP32(real_type);
  if( fread(&dmy, sizeof(int), 1, fp) != 1 ) { fclose(fp); return false; }
  if( eType == UnMatch ) BSWAP32(dmy);
  if( dmy != 8 ) { fclose(fp); return false; }
  
  if( real_type == SPH_FLOAT ) {
    type_dmy = 12;
  } else if( real_type == SPH_DOUBLE ) {
    type_dmy = 24;
  } else {
    return false;
  }
  
  //voxcel
  if( fread(&dmy, sizeof(int), 1, fp) != 1 ) { fclose(fp); return false; }
  if( eType == UnMatch ) BSWAP32(dmy);
  if( dmy != type_dmy ) { fclose(fp); return false; }
  if( real_type == SPH_FLOAT ) {
    unsigned int tmp[3];
    if( fread(tmp, sizeof(int), 3, fp) != 3 ){fclose(fp);return false;}
    if( eType == UnMatch ) {
      BSWAP32(tmp[0]);
      BSWAP32(tmp[1]);
      BSWAP32(tmp[2]);
    }
  } else {
    unsigned long long tmp[3];
    if( fread(tmp, sizeof(long long), 3, fp) != 3 ){fclose(fp);return false;}
    if( eType == UnMatch ) {
      BSWAP64(tmp[0]);
      BSWAP64(tmp[1]);
      BSWAP64(tmp[2]);
    }
  }
  if( fread(&dmy, sizeof(int), 1, fp) != 1 ) { fclose(fp); return false; }
  if( eType == UnMatch ) BSWAP32(dmy);
  if( dmy != type_dmy ) { fclose(fp); return false; }
  
  //org
  if( fread(&dmy, sizeof(int), 1, fp) != 1 ) { fclose(fp); return false; }
  if( eType == UnMatch ) BSWAP32(dmy);
  if( dmy != type_dmy ) { fclose(fp); return false; }
  if( real_type == SPH_FLOAT ) {
    float ftmp[3];
    if( fread(ftmp, sizeof(float), 3, fp) != 3 ){fclose(fp);return false;}
    if( eType == UnMatch ) {
      BSWAP32(ftmp[0]);
      BSWAP32(ftmp[1]);
      BSWAP32(ftmp[2]);
    }
  } else {
    double dtmp[3];
    if( fread(dtmp, sizeof(double), 3, fp) != 3 ){fclose(fp);return false;}
    if( eType == UnMatch ) {
      BSWAP64(dtmp[0]);
      BSWAP64(dtmp[1]);
      BSWAP64(dtmp[2]);
    }
  }
  if( fread(&dmy, sizeof(int), 1, fp) != 1 ) { fclose(fp); return false; }
  if( eType == UnMatch ) BSWAP32(dmy);
  if( dmy != type_dmy ) { fclose(fp); return false; }
  
  //pit
  if( fread(&dmy, sizeof(int), 1, fp) != 1 ) { fclose(fp); return false; }
  if( eType == UnMatch ) BSWAP32(dmy);
  if( dmy != type_dmy ) { fclose(fp); return false; }
  if( real_type == SPH_FLOAT ) {
    float ftmp[3];
    if( fread(ftmp, sizeof(float), 3, fp) != 3 ){fclose(fp);return false;}
    if( eType == UnMatch ) {
      BSWAP32(ftmp[0]);
      BSWAP32(ftmp[1]);
      BSWAP32(ftmp[2]);
    }
  } else {
    double dtmp[3];
    if( fread(dtmp, sizeof(double), 3, fp) != 3 ){fclose(fp);return false;}
    if( eType == UnMatch ) {
      BSWAP64(dtmp[0]);
      BSWAP64(dtmp[1]);
      BSWAP64(dtmp[2]);
    }
  }
  if( fread(&dmy, sizeof(int), 1, fp) != 1 ) { fclose(fp); return false; }
  if( eType == UnMatch ) BSWAP32(dmy);
  if( dmy != type_dmy ) { fclose(fp); return false; }
  
  //step & time
  if( real_type == SPH_FLOAT ) {
    type_dmy = 8;
  } else {
    type_dmy = 16;
  }
  if( fread(&dmy, sizeof(int), 1, fp) != 1 ) { fclose(fp); return false; }
  if( eType == UnMatch ) BSWAP32(dmy);
  if( dmy != type_dmy ) { fclose(fp); return false; }
  if( real_type == SPH_FLOAT ) {
    int t_step;
    if( fread(&t_step, sizeof(int), 1, fp) != 1 ) {fclose(fp); return false;}
    if( eType == UnMatch ) BSWAP32(t_step);
    float t_time;
    if( fread(&t_time, sizeof(float), 1, fp) != 1 ) { fclose(fp); return false; }
    if( eType == UnMatch ) BSWAP32(t_time);
  } else {
    long long t_step;
    if( fread(&t_step, sizeof(long long), 1, fp) != 1 ) {fclose(fp); return false;}
    if( eType == UnMatch ) BSWAP64(t_step);
    double t_time;
    if( fread(&t_time, sizeof(double), 1, fp) != 1 ) { fclose(fp); return false; }
    if( eType == UnMatch ) BSWAP64(t_time);
  }
  if( fread(&dmy, sizeof(int), 1, fp) != 1 ) { fclose(fp); return false; }
  if( eType == UnMatch ) BSWAP32(dmy);
  if( dmy != type_dmy ) { fclose(fp); return false; }
  
  return true;
}


// #################################################################
//
bool COMB::combineXY(bool matchEndian,
                     cdm_Array* buf,
                     cdm_Array* &src,
                     int headS[3],
                     int tailS[3])
{
  
  //copy
  int gcB          = buf->getGcInt();
  const int *headB = buf->getHeadIndex();
  const int *tailB = buf->getTailIndex();
  int       gcS    = src->getGcInt();
  int sta[3],end[3];
  for( int i=0;i<3;i++ )
  {
    sta[i] = (headB[i]-gcB>=headS[i]-gcS) ? headB[i]-gcB : headS[i]-gcS;
    end[i] = (tailB[i]+gcB<=tailS[i]+gcS) ? tailB[i]+gcB : tailS[i]+gcS;
  }
  
  
  //同じデータ型のコピー
  if( buf->getDataType() == src->getDataType() ) {
    // float to float
    if( buf->getDataType() == CDM::E_CDM_FLOAT32 ) {
      cdm_TypeArray<float> *B = dynamic_cast<cdm_TypeArray<float>*>(buf);
      cdm_TypeArray<float> *S = dynamic_cast<cdm_TypeArray<float>*>(src);
      copyArray(B, S, sta, end);
      //copyArray(buf, src, sta, end);
      // double to double
    } else {
      cdm_TypeArray<double> *B = dynamic_cast<cdm_TypeArray<double>*>(buf);
      cdm_TypeArray<double> *S = dynamic_cast<cdm_TypeArray<double>*>(src);
      copyArray(B, S, sta, end);
    }
    //違う型のコピー
  } else {
    // float to double
    if( buf->getDataType() == CDM::E_CDM_FLOAT32 &&
       src->getDataType() == CDM::E_CDM_FLOAT64 ) {
      cdm_TypeArray<float>  *B = dynamic_cast<cdm_TypeArray<float>*>(buf);
      cdm_TypeArray<double> *S = dynamic_cast<cdm_TypeArray<double>*>(src);
      copyArray(B, S, sta, end);
      //doubel to float
    } else {
      cdm_TypeArray<double> *B = dynamic_cast<cdm_TypeArray<double>*>(buf);
      cdm_TypeArray<float>  *S = dynamic_cast<cdm_TypeArray<float>*>(src);
      copyArray(B, S, sta, end);
    }
  }
  
  return true;
}
