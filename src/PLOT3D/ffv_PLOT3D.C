//##################################################################################
//
// FFV-C : Frontflow / violet Cartesian
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
 * @file   ffv_PLOT3D.C
 * @brief  Plot3D Class
 * @author aics
 */

#include "ffv_PLOT3D.h"


// #################################################################
// PLOT3Dファイル入出力に関するパラメータ取得
void Plot3D::getParameter(TextParser* tpCntl)
{
  string str;
  string label;
  
  // Output Directory_Path
  label = "/Output/FormatOption/PLOT3D/DirectoryPath";
  
  if ( !(tpCntl->getInspectedValue(label, str)) )
  {
    if (myRank==0) stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    Exit(0);
  }
  // 指定が無ければ，空のまま
  if ( !str.empty() )
  {
    C->FIO.OutDirPath = str;
  }
  
  // FileNameGrid --- option
  label = "/Output/FormatOption/PLOT3D/FileNameGrid";
  
  if ( !(tpCntl->getInspectedValue(label, str)) )
  {
    C->P3Op.basename_g = "PLOT3DoutputGrid";
  }
  else
  {
    C->P3Op.basename_g = str;
  }
  
  if ( C->P3Op.basename_g.empty() )
  {
    C->P3Op.basename_g = "PLOT3DoutputGrid";
  }
  
  // FileNameFunc --- option
  label = "/Output/FormatOption/PLOT3D/FileNameFunc";
  
  if ( !(tpCntl->getInspectedValue(label, str)) )
  {
    C->P3Op.basename_f = "PLOT3Doutput";
  }
  else
  {
    C->P3Op.basename_f = str;
  }
  
  if ( C->P3Op.basename_f.empty() )
  {
    C->P3Op.basename_f = "PLOT3Doutput";
  }
  
  // GridKind
  /*
   label = "/Steer/Plot3dOptions/GridKind";
   
   if ( !(tpCntl->getInspectedValue(label, str )) )
   {
   Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
   Exit(0);
   }
   else
   {
   if     ( !strcasecmp(str.c_str(), "SingleGrid") ) FP3DR->setSingleGrid();
   else if( !strcasecmp(str.c_str(), "MultiGrid") )  FP3DR->setMultiGrid();
   else
   {
   Hostonly_ stamped_printf("\tInvalid keyword is described for '%s'\n", label.c_str());
   Exit(0);
   }
   }
   */
  
  FP3DR->setMultiGrid(); // 常にmulti grid
  
  
  // 格子の移動
  label = "/Output/FormatOption/PLOT3D/GridMobility";
  
  if ( !(tpCntl->getInspectedValue(label, str )) )
  {
    if (myRank==0) stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    Exit(0);
  }
  else
  {
    if     ( !strcasecmp(str.c_str(), "movable") )   FP3DR->setMoveGrid(GRID_MOVE);
    else if( !strcasecmp(str.c_str(), "immovable") ) FP3DR->setMoveGrid(GRID_NOT_MOVE);
    else
    {
      if (myRank==0) stamped_printf("\tInvalid keyword is described for '%s'\n", label.c_str());
      Exit(0);
    }
  }
  
  
  // 時間方向の変化
  label = "/Output/FormatOption/PLOT3D/StateOfTime";
  
  if ( !(tpCntl->getInspectedValue(label, str )) )
  {
    FP3DR->setSteady(FB_UNSTEADY);
  }
  else
  {
    if     ( !strcasecmp(str.c_str(), "steady") )   FP3DR->setSteady(FB_STEADY);
    else if( !strcasecmp(str.c_str(), "unsteady") ) FP3DR->setSteady(FB_UNSTEADY);
    else
    {
      if (myRank==0) stamped_printf("\tInvalid keyword is described for '%s'\n", label.c_str());
      Exit(0);
    }
  }
  
  
  // IBLANKファイル
  label = "/Output/FormatOption/PLOT3D/SetIblankFlag";
  
  if ( !(tpCntl->getInspectedValue(label, str )) )
  {
    FP3DR->setIBlankFlag(SET_IBLANK);
  }
  else
  {
    if     ( !strcasecmp(str.c_str(), "on") )  FP3DR->setIBlankFlag(SET_IBLANK);
    else if( !strcasecmp(str.c_str(), "off") ) FP3DR->setIBlankFlag(NOT_SET_IBLANK);
    else
    {
      if (myRank==0) stamped_printf("\tInvalid keyword is described for '%s'\n", label.c_str());
      Exit(0);
    }
  }
  
  
  // 次元数
  /*
   label = "/Steer/Plot3dOptions/Dimension";
   if ( !(tpCntl->getInspectedValue(label, str )) )
   {
   //Hostonly_ stamped_printf("\tParsing error : fail to get '/Steer/plot3doptions/dimension'\n");
   //Exit(0);
   FP3DR->setDimension3D();
   }
   else
   {
   if     ( !strcasecmp(str.c_str(), "2d") ) FP3DR->setDimension2D();
   else if( !strcasecmp(str.c_str(), "3d") ) FP3DR->setDimension3D();
   else
   {
   Hostonly_ stamped_printf("\tInvalid keyword is described for '%s'\n", label.c_str());
   Exit(0);
   }
   }
   */
  FP3DR->setDimension3D(); // 常に三次元
  
  // FormatType
  label = "/Output/FormatOption/PLOT3D/FormatType";
  
  if ( !(tpCntl->getInspectedValue(label, str )) )
  {
    if (myRank==0) stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    Exit(0);
  }
  else
  {
    if     ( !strcasecmp(str.c_str(), "unformatted") )  FP3DR->setFormat(UNFORMATTED);
    else if( !strcasecmp(str.c_str(), "formatted") )    FP3DR->setFormat(FORMATTED);
    else if( !strcasecmp(str.c_str(), "binary") )       FP3DR->setFormat(C_BINARY);
    else
    {
      if (myRank==0) stamped_printf("\tInvalid keyword is described for '%s'\n", label.c_str());
      Exit(0);
    }
  }
  
  
  // 出力の単精度or倍精度指定
  label = "/Output/FormatOption/PLOT3D/RealType";
  
  if ( !(tpCntl->getInspectedValue(label, str )) )
  {
    //Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    //Exit(0);
    int d_type = (sizeof(REAL_TYPE) == 4) ? 1 : 2;  // 1-float / 2-double
    FP3DR->setRealType(d_type);
  }
  else
  {
    if     ( !strcasecmp(str.c_str(), "float") )  FP3DR->setRealType(OUTPUT_FLOAT);
    else if( !strcasecmp(str.c_str(), "double") ) FP3DR->setRealType(OUTPUT_DOUBLE);
    else
    {
      if (myRank==0) stamped_printf("\tInvalid keyword is described for '%s'\n", label.c_str());
      Exit(0);
    }
  }
  
  
  //copy options read to write
  FP3DW->setGridKind(FP3DR->IsGridKind());
  FP3DW->setMoveGrid(FP3DR->IsMoveGrid());
  FP3DW->setSteady(FP3DR->IsSteady());
  FP3DW->setIBlankFlag(FP3DR->IsIBlankFlag());
  FP3DW->setDim(FP3DR->GetDim());
  FP3DW->setFormat(FP3DR->GetFormat());
  FP3DW->setRealType(FP3DR->GetRealType());
  
  
  // OutputXyz
  label = "/Output/FormatOption/PLOT3D/OutputXyz";
  
  if ( !(tpCntl->getInspectedValue(label, str )) )
  {
    if (myRank==0) stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    Exit(0);
  }
  else
  {
    if     ( !strcasecmp(str.c_str(), "on") )  C->P3Op.IS_xyz = ON;
    else if( !strcasecmp(str.c_str(), "off") ) C->P3Op.IS_xyz = OFF;
    else
    {
      if (myRank==0) stamped_printf("\tInvalid keyword is described for '%s'\n", label.c_str());
      Exit(0);
    }
  }
  
  
  // OutputQ
  
  /*
   label = "/Steer/plot3doptions/OutputQ";
   
   if ( !(tpCntl->getInspectedValue(label, str )) )
   {
   //Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
   //Exit(0);
   P3Op.IS_q = OFF;
   }
   else
   {
   if     ( !strcasecmp(str.c_str(), "on") )  P3Op.IS_q = ON;
   else if( !strcasecmp(str.c_str(), "off") ) P3Op.IS_q = OFF;
   else
   {
   Hostonly_ stamped_printf("\tInvalid keyword is described for '%s'\n", label.c_str());
   Exit(0);
   }
   }*/
  
  C->P3Op.IS_q = OFF; // 常にoff
  
  
  // OutputFunction
  label = "/Output/FormatOption/PLOT3D/OutputFunction";
  
  if ( !(tpCntl->getInspectedValue(label, str )) )
  {
    if (myRank==0) stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    Exit(0);
  }
  else
  {
    if     ( !strcasecmp(str.c_str(), "on") )  C->P3Op.IS_funciton = ON;
    else if( !strcasecmp(str.c_str(), "off") ) C->P3Op.IS_funciton = OFF;
    else
    {
      if (myRank==0) stamped_printf("\tInvalid keyword is described for '%s'\n", label.c_str());
      Exit(0);
    }
  }
  
  
  // OutputFuncName
  label = "/Output/FormatOption/PLOT3D/OutputFuncName";
  
  if ( !(tpCntl->getInspectedValue(label, str )) )
  {
    if (myRank==0) stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    Exit(0);
  }
  else
  {
    if     ( !strcasecmp(str.c_str(), "on") )  C->P3Op.IS_function_name = ON;
    else if( !strcasecmp(str.c_str(), "off") ) C->P3Op.IS_function_name = OFF;
    else
    {
      if (myRank==0) stamped_printf("\tInvalid keyword is described for '%s'\n", label.c_str());
      Exit(0);
    }
  }
  
  
  // OutputFvbnd
  
  /*
   label = "/Steer/plot3doptions/OutputFvbnd";
   
   if ( !(tpCntl->getInspectedValue(label, str )) )
   {
   //Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
   //Exit(0);
   P3Op.IS_fvbnd = OFF;
   }
   else
   {
   if     ( !strcasecmp(str.c_str(), "on") )  P3Op.IS_fvbnd = ON;
   else if( !strcasecmp(str.c_str(), "off") ) P3Op.IS_fvbnd = OFF;
   else
   {
   Hostonly_ stamped_printf("\tInvalid keyword is described for '%s'\n", label.c_str());
   Exit(0);
   }
   }*/
  
  C->P3Op.IS_fvbnd = OFF; // 常にoff
  
  // DivideFunc ---> 出力を項目別にファイル分割するオプション
  label = "/Output/FormatOption/PLOT3D/DivideFunc";
  
  if ( !(tpCntl->getInspectedValue(label, str )) )
  {
    if (myRank==0) stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    Exit(0);
  }
  else
  {
    if     ( !strcasecmp(str.c_str(), "on") )  C->P3Op.IS_DivideFunc = ON;
    else if( !strcasecmp(str.c_str(), "off") ) C->P3Op.IS_DivideFunc = OFF;
    else
    {
      if (myRank==0) stamped_printf("\tInvalid keyword is described for '%s'\n", label.c_str());
      Exit(0);
    }
  }
  
}



// #################################################################
void Plot3D::Initialize(const int* m_size,
                        const int m_guide,
                        const REAL_TYPE m_deltaX,
                        Control* m_C,
                        FileIO_PLOT3D_READ* m_FP3DR,
                        FileIO_PLOT3D_WRITE* m_FP3DW,
                        DFI* m_dfi,
                        REAL_TYPE* m_d_ws,
                        REAL_TYPE* m_d_p,
                        REAL_TYPE* m_d_wo,
                        REAL_TYPE* m_d_v,
                        REAL_TYPE* m_d_ie,
                        REAL_TYPE* m_d_p0,
                        REAL_TYPE* m_d_wv,
                        int*       m_d_cdf,
                        int*       m_d_bcd)
{
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
  
  size[0] = m_size[0];
  size[1] = m_size[1];
  size[2] = m_size[2];
  guide   = m_guide;
  deltaX  = m_deltaX;
  
  C    = m_C;
  FP3DR= m_FP3DR;
  FP3DW= m_FP3DW;
  dfi  = m_dfi;
  
  d_ws = m_d_ws;
  d_p  = m_d_p;
  d_wo = m_d_wo;
  d_v  = m_d_v;
  d_ie = m_d_ie;
  d_p0 = m_d_p0;
  d_wv = m_d_wv;
  d_cdf= m_d_cdf;
  d_bcd= m_d_bcd;
}


// #################################################################
void Plot3D::function(const unsigned CurrentStep,
                      const double CurrentTime,
                      REAL_TYPE* v00,
                      double& flop)
{
  //value
  int igrid;
  int ngrid;
  int *id,*jd,*kd;
  int *nvar;
  int ix,jx,kx;
  float *d;
  double *dd;
  
  //
  REAL_TYPE scale = 1.0;
  
  // ガイドセル出力
  //int gc_out = C->GuideOut;
  int gc_out = 0;//plot3dは常にガイドセルを出力しない
  
  // ステップ数
  int m_step = (int)CurrentStep;
  
  // 時間の次元変換
  REAL_TYPE m_time;
  if (C->Unit.File == DIMENSIONAL)
  {
    m_time = (REAL_TYPE)(CurrentTime * C->Tscale);
  }
  else
  {
    m_time = (REAL_TYPE)CurrentTime;
  }
  
  //allocate
  ngrid = C->P3Op.ngrid;
  id = new int[ngrid];
  jd = new int[ngrid];
  kd = new int[ngrid];
  nvar = new int[ngrid];
  
  //set grid data and nvar and work area size
  int maxid=0;
  int maxjd=0;
  int maxkd=0;
  int maxnvar=0;
  
  //for(igrid=0;igrid<ngrid;igrid++){ //--->BCMでループが必要になる？
  
  igrid=0;
  id[igrid]=size[0]+1+2*gc_out;
  jd[igrid]=size[1]+1+2*gc_out;
  kd[igrid]=size[2]+1+2*gc_out;
  nvar[igrid]=C->P3Op.nvar;
  
  if(maxid<id[igrid]) maxid=id[igrid];
  if(maxjd<jd[igrid]) maxjd=jd[igrid];
  if(maxkd<kd[igrid]) maxkd=kd[igrid];
  if(maxnvar<nvar[igrid]) maxnvar=nvar[igrid];
  
  
  //allocate workarea
  if(FP3DW->GetFormat() == C_BINARY){//C_BINARYでの出力は項目ごとに書き出し
    if( FP3DW->GetRealType() == OUTPUT_FLOAT ){
      if (!(d = new float[ maxid*maxjd*maxkd ])){
        if (myRank==0) printf(    "\t>> cannot allocate work area : function()\n\n");
        Exit(0);
      }
    }else if( FP3DW->GetRealType() == OUTPUT_DOUBLE ){
      if (!(dd = new double[ maxid*maxjd*maxkd ])){
        if (myRank==0) printf(    "\t>> cannot allocate work area : function()\n\n");
        Exit(0);
      }
    }
  }
  else//Fortranによる出力では出力項目すべてを一度に書き出し
  {
    if( FP3DW->GetRealType() == OUTPUT_FLOAT ){
      if (!(d = new float[ maxid*maxjd*maxkd*maxnvar ])){
        if (myRank==0) printf(    "\t>> cannot allocate work area : function()\n\n");
        Exit(0);
      }
    }else if( FP3DW->GetRealType() == OUTPUT_DOUBLE ){
      if (!(dd = new double[ maxid*maxjd*maxkd*maxnvar ])){
        if (myRank==0) printf(    "\t>> cannot allocate work area : function()\n\n");
        Exit(0);
      }
    }
  }

  
  // 出力ファイル名
  std::string dtmp = dfi->GenerateDirName(C->FIO.OutDirPath, m_step, C->FIO.Slice);
  
  // 出力ディレクトリの作成
  if ( !FBUtility::mkdirs(dtmp) ) {
    if (myRank==0) printf("Error : create directory \"%s\"\n", dtmp.c_str());
    Exit(-1);
  }
  
  std::string tmp = dfi->GenerateFileName(C->P3Op.basename_f, "func", m_step, myRank, true);
  
  //open file
  FP3DW->setFileName((dtmp+tmp).c_str());
  if(!FP3DW->OpenFile()){
    std::cout << "error OpenFile" << std::endl;
    Exit(0);
  }
  
  //write block data
  FP3DW->WriteNgrid(ngrid);
  FP3DW->WriteFuncBlockData(id,jd,kd,nvar,ngrid);
  
  //write function data
  
  //for(igrid=0;igrid<ngrid;igrid++){ //--->BCMでループが必要になる？
  igrid=0;
  
  size_t out_length = (size_t)(id[igrid]*jd[igrid]*kd[igrid]);
  
  //set grid data
  FP3DW->setGridData(id[igrid],jd[igrid],kd[igrid],ngrid);
  FP3DW->setFuncDataNum(nvar[igrid]);
  
  //output start
  int ivar=0;
  
  // Pressure
  if (C->Unit.File == DIMENSIONAL)
  {
    REAL_TYPE bp = ( C->Unit.Prs == Unit_Absolute ) ? C->BasePrs : 0.0;
    U.convArrayPrsND2D(d_ws, size, guide, d_p, bp, C->RefDensity, C->RefVelocity, flop);
  }
  else
  {
    U.copyS3D(d_ws, size, guide, d_p, scale);
  }
  
  if(FP3DW->GetFormat() == C_BINARY){//C_BINARYでの出力は項目ごとに書き出し
    if( FP3DW->GetRealType() == OUTPUT_FLOAT ){
      setScalarGridData(d,d_ws,id[igrid],jd[igrid],kd[igrid],gc_out);
      FP3DW->setFuncData(d);
    }else if( FP3DW->GetRealType() == OUTPUT_DOUBLE ){
      setScalarGridData(dd,d_ws,id[igrid],jd[igrid],kd[igrid],gc_out);
      FP3DW->setFuncData(dd);
    }
    if(!FP3DW->WriteFuncData()) std::cout << "error WriteFuncData" << std::endl;
  }
  else//Fortranによる出力では出力項目すべてを一度に書き出し
  {
    if( FP3DW->GetRealType() == OUTPUT_FLOAT ){
      setScalarGridData(&d[out_length*ivar],d_ws,id[igrid],jd[igrid],kd[igrid],gc_out);
    }else if( FP3DW->GetRealType() == OUTPUT_DOUBLE ){
      setScalarGridData(&dd[out_length*ivar],d_ws,id[igrid],jd[igrid],kd[igrid],gc_out);
    }
  }
  ivar++;
  
  // Velocity
  REAL_TYPE unit_velocity = (C->Unit.File == DIMENSIONAL) ? C->RefVelocity : 1.0;
  fb_vout_nijk_(d_wo, d_v, size, &guide, v00, &unit_velocity, &flop);
  
  if(FP3DW->GetFormat() == C_BINARY){//C_BINARYでの出力は項目ごとに書き出し
    if( FP3DW->GetRealType() == OUTPUT_FLOAT ){
      setVectorComponentGridData(d,d_wo,id[igrid],jd[igrid],kd[igrid],0,gc_out);//0:x
      FP3DW->setFuncData(d);
    }else if( FP3DW->GetRealType() == OUTPUT_DOUBLE ){
      setVectorComponentGridData(dd,d_wo,id[igrid],jd[igrid],kd[igrid],0,gc_out);//0:x
      FP3DW->setFuncData(dd);
    }
    if(!FP3DW->WriteFuncData()) std::cout << "error WriteFuncData" << std::endl;
    if( FP3DW->GetRealType() == OUTPUT_FLOAT ){
      setVectorComponentGridData(d,d_wo,id[igrid],jd[igrid],kd[igrid],1,gc_out);//1:y
      FP3DW->setFuncData(d);
    }else if( FP3DW->GetRealType() == OUTPUT_DOUBLE ){
      setVectorComponentGridData(dd,d_wo,id[igrid],jd[igrid],kd[igrid],1,gc_out);//1:y
      FP3DW->setFuncData(dd);
    }
    if(!FP3DW->WriteFuncData()) std::cout << "error WriteFuncData" << std::endl;
    if( FP3DW->GetRealType() == OUTPUT_FLOAT ){
      setVectorComponentGridData(d,d_wo,id[igrid],jd[igrid],kd[igrid],2,gc_out);//2:z
      FP3DW->setFuncData(d);
    }else if( FP3DW->GetRealType() == OUTPUT_DOUBLE ){
      setVectorComponentGridData(dd,d_wo,id[igrid],jd[igrid],kd[igrid],2,gc_out);//2:z
      FP3DW->setFuncData(dd);
    }
    if(!FP3DW->WriteFuncData()) std::cout << "error WriteFuncData" << std::endl;
  }
  else//Fortranによる出力では出力項目すべてを一度に書き出し
  {
    if( FP3DW->GetRealType() == OUTPUT_FLOAT ){
      setVectorGridData(&d[out_length*ivar],d_wo,id[igrid],jd[igrid],kd[igrid],gc_out);
    }else if( FP3DW->GetRealType() == OUTPUT_DOUBLE ){
      setVectorGridData(&dd[out_length*ivar],d_wo,id[igrid],jd[igrid],kd[igrid],gc_out);
    }
  }
  ivar=ivar+3;
  
  // Tempearture
  if( C->isHeatProblem() )
  {
    if (C->Unit.File == DIMENSIONAL)
    {
      //U.convArrayIE2Tmp(d_ws, size, guide, d_ie, d_bcd, mat_tbl, C->BaseTemp, C->DiffTemp, true, flop);
    }
    else
    {
      //U.convArrayIE2Tmp(d_ws, size, guide, d_ie, d_bcd, mat_tbl, C->BaseTemp, C->DiffTemp, false, flop);
    }
    
    if(FP3DW->GetFormat() == C_BINARY){//C_BINARYでの出力は項目ごとに書き出し
      if( FP3DW->GetRealType() == OUTPUT_FLOAT ){
        setScalarGridData(d,d_ws,id[igrid],jd[igrid],kd[igrid],gc_out);
        FP3DW->setFuncData(d);
      }else if( FP3DW->GetRealType() == OUTPUT_DOUBLE ){
        setScalarGridData(dd,d_ws,id[igrid],jd[igrid],kd[igrid],gc_out);
        FP3DW->setFuncData(dd);
      }
      if(!FP3DW->WriteFuncData()) std::cout << "error WriteFuncData" << std::endl;
    }
    else//Fortranによる出力では出力項目すべてを一度に書き出し
    {
      if( FP3DW->GetRealType() == OUTPUT_FLOAT ){
        setScalarGridData(&d[out_length*ivar],d_ws,id[igrid],jd[igrid],kd[igrid],gc_out);
      }else if( FP3DW->GetRealType() == OUTPUT_DOUBLE ){
        setScalarGridData(&dd[out_length*ivar],d_ws,id[igrid],jd[igrid],kd[igrid],gc_out);
      }
    }
    ivar++;
    
  }
  
  // Total Pressure
  if (C->varState[var_TotalP] == ON ) {
    fb_totalp_ (d_p0, size, &guide, d_v, d_p, v00, &flop);
    
    // convert non-dimensional to dimensional, iff file is dimensional
    if (C->Unit.File == DIMENSIONAL)
    {
      U.convArrayTpND2D(d_ws, d_p0, size, guide, C->RefDensity, C->RefVelocity);
    }
    else
    {
      REAL_TYPE* tp;
      tp = d_ws; d_ws = d_p0; d_p0 = tp;
    }
    
    if(FP3DW->GetFormat() == C_BINARY){//C_BINARYでの出力は項目ごとに書き出し
      if( FP3DW->GetRealType() == OUTPUT_FLOAT ){
        setScalarGridData(d,d_ws,id[igrid],jd[igrid],kd[igrid],gc_out);
        FP3DW->setFuncData(d);
      }else if( FP3DW->GetRealType() == OUTPUT_DOUBLE ){
        setScalarGridData(dd,d_ws,id[igrid],jd[igrid],kd[igrid],gc_out);
        FP3DW->setFuncData(dd);
      }
      if(!FP3DW->WriteFuncData()) std::cout << "error WriteFuncData" << std::endl;
    }
    else//Fortranによる出力では出力項目すべてを一度に書き出し
    {
      if( FP3DW->GetRealType() == OUTPUT_FLOAT ){
        setScalarGridData(&d[out_length*ivar],d_ws,id[igrid],jd[igrid],kd[igrid],gc_out);
      }else if( FP3DW->GetRealType() == OUTPUT_DOUBLE ){
        setScalarGridData(&dd[out_length*ivar],d_ws,id[igrid],jd[igrid],kd[igrid],gc_out);
      }
    }
    ivar++;
  }
  
  // Vorticity
  if (C->varState[var_Vorticity] == ON )
  {
    rot_v_(d_wv, size, &guide, &deltaX, d_v, d_cdf, v00, &flop);
    
    REAL_TYPE  vz[3];
    vz[0] = vz[1] = vz[2] = 0.0;
    unit_velocity = (C->Unit.File == DIMENSIONAL) ? C->RefVelocity/C->RefLength : 1.0;
    fb_vout_nijk_(d_wo, d_wv, size, &guide, vz, &unit_velocity, &flop);
    
    if(FP3DW->GetFormat() == C_BINARY){//C_BINARYでの出力は項目ごとに書き出し
      if( FP3DW->GetRealType() == OUTPUT_FLOAT ){
        setVectorComponentGridData(d,d_wo,id[igrid],jd[igrid],kd[igrid],0,gc_out);//0:x
        FP3DW->setFuncData(d);
      }else if( FP3DW->GetRealType() == OUTPUT_DOUBLE ){
        setVectorComponentGridData(dd,d_wo,id[igrid],jd[igrid],kd[igrid],0,gc_out);//0:x
        FP3DW->setFuncData(dd);
      }
      if(!FP3DW->WriteFuncData()) std::cout << "error WriteFuncData" << std::endl;
      if( FP3DW->GetRealType() == OUTPUT_FLOAT ){
        setVectorComponentGridData(d,d_wo,id[igrid],jd[igrid],kd[igrid],1,gc_out);//1:y
        FP3DW->setFuncData(d);
      }else if( FP3DW->GetRealType() == OUTPUT_DOUBLE ){
        setVectorComponentGridData(dd,d_wo,id[igrid],jd[igrid],kd[igrid],1,gc_out);//1:y
        FP3DW->setFuncData(dd);
      }
      if(!FP3DW->WriteFuncData()) std::cout << "error WriteFuncData" << std::endl;
      if( FP3DW->GetRealType() == OUTPUT_FLOAT ){
        setVectorComponentGridData(d,d_wo,id[igrid],jd[igrid],kd[igrid],2,gc_out);//2:z
        FP3DW->setFuncData(d);
      }else if( FP3DW->GetRealType() == OUTPUT_DOUBLE ){
        setVectorComponentGridData(dd,d_wo,id[igrid],jd[igrid],kd[igrid],2,gc_out);//2:z
        FP3DW->setFuncData(dd);
      }
      if(!FP3DW->WriteFuncData()) std::cout << "error WriteFuncData" << std::endl;
    }
    else//Fortranによる出力では出力項目すべてを一度に書き出し
    {
      if( FP3DW->GetRealType() == OUTPUT_FLOAT ){
        setVectorGridData(&d[out_length*ivar],d_wo,id[igrid],jd[igrid],kd[igrid],gc_out);
      }else if( FP3DW->GetRealType() == OUTPUT_DOUBLE ){
        setVectorGridData(&dd[out_length*ivar],d_wo,id[igrid],jd[igrid],kd[igrid],gc_out);
      }
    }
    ivar=ivar+3;
  }
  
  // 2nd Invariant of Velocity Gradient Tensor
  if (C->varState[var_Qcr] == ON ) {
    i2vgt_ (d_p0, size, &guide, &deltaX, d_v, d_cdf, v00, &flop);
    
    U.copyS3D(d_ws, size, guide, d_p0, scale);
    
    if(FP3DW->GetFormat() == C_BINARY){//C_BINARYでの出力は項目ごとに書き出し
      if( FP3DW->GetRealType() == OUTPUT_FLOAT ){
        setScalarGridData(d,d_ws,id[igrid],jd[igrid],kd[igrid],gc_out);
        FP3DW->setFuncData(d);
      }else if( FP3DW->GetRealType() == OUTPUT_DOUBLE ){
        setScalarGridData(dd,d_ws,id[igrid],jd[igrid],kd[igrid],gc_out);
        FP3DW->setFuncData(dd);
      }
      if(!FP3DW->WriteFuncData()) std::cout << "error WriteFuncData" << std::endl;
    }
    else//Fortranによる出力では出力項目すべてを一度に書き出し
    {
      if( FP3DW->GetRealType() == OUTPUT_FLOAT ){
        setScalarGridData(&d[out_length*ivar],d_ws,id[igrid],jd[igrid],kd[igrid],gc_out);
      }else if( FP3DW->GetRealType() == OUTPUT_DOUBLE ){
        setScalarGridData(&dd[out_length*ivar],d_ws,id[igrid],jd[igrid],kd[igrid],gc_out);
      }
    }
    ivar++;
  }
  
  // Helicity
  if (C->varState[var_Helicity] == ON )
  {
    helicity_(d_p0, size, &guide, &deltaX, d_v, d_cdf, v00, &flop);
    
    U.copyS3D(d_ws, size, guide, d_p0, scale);
    
    if(FP3DW->GetFormat() == C_BINARY){//C_BINARYでの出力は項目ごとに書き出し
      if( FP3DW->GetRealType() == OUTPUT_FLOAT ){
        setScalarGridData(d,d_ws,id[igrid],jd[igrid],kd[igrid],gc_out);
        FP3DW->setFuncData(d);
      }else if( FP3DW->GetRealType() == OUTPUT_DOUBLE ){
        setScalarGridData(dd,d_ws,id[igrid],jd[igrid],kd[igrid],gc_out);
        FP3DW->setFuncData(dd);
      }
      if(!FP3DW->WriteFuncData()) std::cout << "error WriteFuncData" << std::endl;
    }
    else//Fortranによる出力では出力項目すべてを一度に書き出し
    {
      if( FP3DW->GetRealType() == OUTPUT_FLOAT ){
        setScalarGridData(&d[out_length*ivar],d_ws,id[igrid],jd[igrid],kd[igrid],gc_out);
      }else if( FP3DW->GetRealType() == OUTPUT_DOUBLE ){
        setScalarGridData(&dd[out_length*ivar],d_ws,id[igrid],jd[igrid],kd[igrid],gc_out);
      }
    }
    ivar++;
  }
  
  //write all
  if(FP3DW->GetFormat() != C_BINARY){//C_BINARY以外は出力項目すべてを一度に書き出し
    if( FP3DW->GetRealType() == OUTPUT_FLOAT ){
      FP3DW->setFuncData(d);
    }else if( FP3DW->GetRealType() == OUTPUT_DOUBLE ){
      FP3DW->setFuncData(dd);
    }
    if(!FP3DW->WriteFuncData()) std::cout << "error WriteFuncData" << std::endl;
  }
  
  //}//igrid loop
  
  //close file
  FP3DW->CloseFile();
  
  //deallocate
  if( FP3DW->GetRealType() == OUTPUT_FLOAT ){
    delete [] d;
  }else if( FP3DW->GetRealType() == OUTPUT_DOUBLE ){
    delete [] dd;
  }
  delete [] id;
  delete [] jd;
  delete [] kd;
}


// #################################################################
void Plot3D::function_divide(const unsigned CurrentStep,
                             const double CurrentTime,
                             REAL_TYPE* v00,
                             double& flop)
{
  //value
  int igrid;
  int ngrid;
  int *id,*jd,*kd;
  int *nvar;
  int ix,jx,kx;
  float *d;
  double *dd;
  int *nvar_scalar;
  int *nvar_vector;
  
  //
  REAL_TYPE scale = 1.0;
  
  // ガイドセル出力
  //int gc_out = C->GuideOut;
  int gc_out = 0;//plot3dは常にガイドセルを出力しない
  
  // ステップ数
  int m_step = (int)CurrentStep;
  
  // 時間の次元変換
  REAL_TYPE m_time;
  if (C->Unit.File == DIMENSIONAL)
  {
    m_time = (REAL_TYPE)(CurrentTime * C->Tscale);
  }
  else
  {
    m_time = (REAL_TYPE)CurrentTime;
  }
  
  //allocate
  ngrid=C->P3Op.ngrid;
  id = new int[ngrid];
  jd = new int[ngrid];
  kd = new int[ngrid];
  nvar = new int[ngrid];
  nvar_scalar = new int[ngrid];
  nvar_vector = new int[ngrid];
  
  //set grid data and nvar and work area size
  int maxid=0;
  int maxjd=0;
  int maxkd=0;
  
  //for(igrid=0;igrid<ngrid;igrid++){ //--->BCMでループが必要になる？
  
  igrid=0;
  id[igrid]=size[0]+1+2*gc_out;
  jd[igrid]=size[1]+1+2*gc_out;
  kd[igrid]=size[2]+1+2*gc_out;
  nvar[igrid]=C->P3Op.nvar;
  nvar_scalar[igrid]=1;
  nvar_vector[igrid]=3;
  
  if(maxid<id[igrid]) maxid=id[igrid];
  if(maxjd<jd[igrid]) maxjd=jd[igrid];
  if(maxkd<kd[igrid]) maxkd=kd[igrid];
  
  //}//igrid loop
  
  //allocate workarea
  if( FP3DW->GetRealType() == OUTPUT_FLOAT ){
    if (!(d = new float[ maxid*maxjd*maxkd ])){
      if (myRank==0)  printf(    "\t>> cannot allocate work area : function()\n\n");
      Exit(0);
    }
  }else if( FP3DW->GetRealType() == OUTPUT_DOUBLE ){
    if (!(dd = new double[ maxid*maxjd*maxkd ])){
      if (myRank==0)  printf(    "\t>> cannot allocate work area : function()\n\n");
      Exit(0);
    }
  }
  
  
  // 出力ファイル名
  std::string fname,dfi_name;
  std::string dtmp = dfi->GenerateDirName(C->FIO.OutDirPath, m_step, C->FIO.Slice);
  
  // 出力ディレクトリの作成
  if ( !FBUtility::mkdirs(dtmp) ) {
    if (myRank==0) printf("Error : create directory \"%s\"\n", dtmp.c_str());
    Exit(-1);
  }
  
  std::string tmp  = dfi->GenerateFileName(C->P3Op.basename_f, "func", m_step, myRank, true);
  
  
  // Pressure
  if (C->Unit.File == DIMENSIONAL)
  {
    REAL_TYPE bp = ( C->Unit.Prs == Unit_Absolute ) ? C->BasePrs : 0.0;
    U.convArrayPrsND2D(d_ws, size, guide, d_p, bp, C->RefDensity, C->RefVelocity, flop);
  }
  else
  {
    U.copyS3D(d_ws, size, guide, d_p, scale);
  }
  
  fname = "prs_" + tmp;
  FP3DW->setFileName((dtmp+fname).c_str());
  if(!FP3DW->OpenFile()){
    std::cout << "error OpenFile" << std::endl;
    Exit(0);
  }
  
  FP3DW->WriteNgrid(ngrid);
  FP3DW->WriteFuncBlockData(id,jd,kd,nvar_scalar,ngrid);
  
  //for(igrid=0;igrid<ngrid;igrid++){ //--->BCMでループが必要になる？
  igrid=0;
  FP3DW->setGridData(id[igrid],jd[igrid],kd[igrid],ngrid);
  FP3DW->setFuncDataNum(nvar_scalar[igrid]);
  if( FP3DW->GetRealType() == OUTPUT_FLOAT ){
    setScalarGridData(d,d_ws,id[igrid],jd[igrid],kd[igrid],gc_out);
    FP3DW->setFuncData(d);
  }else if( FP3DW->GetRealType() == OUTPUT_DOUBLE ){
    setScalarGridData(dd,d_ws,id[igrid],jd[igrid],kd[igrid],gc_out);
    FP3DW->setFuncData(dd);
  }
  if(!FP3DW->WriteFuncData()) std::cout << "error WriteFuncData" << std::endl;
  //}//igrid loop
  FP3DW->CloseFile();
  dfi_name = "prs_" + C->P3Op.basename_f;

  
  // Velocity
  REAL_TYPE unit_velocity = (C->Unit.File == DIMENSIONAL) ? C->RefVelocity : 1.0;
  fb_vout_nijk_(d_wo, d_v, size, &guide, v00, &unit_velocity, &flop);
  
  fname = "vel_" + tmp;
  FP3DW->setFileName((dtmp+fname).c_str());
  if(!FP3DW->OpenFile()){
    std::cout << "error OpenFile" << std::endl;
    Exit(0);
  }
  
  FP3DW->WriteNgrid(ngrid);
  FP3DW->WriteFuncBlockData(id,jd,kd,nvar_vector,ngrid);
  
  //for(igrid=0;igrid<ngrid;igrid++){ //--->BCMでループが必要になる？
  igrid=0;
  FP3DW->setGridData(id[igrid],jd[igrid],kd[igrid],ngrid);
  FP3DW->setFuncDataNum(nvar_vector[igrid]);
  if( FP3DW->GetRealType() == OUTPUT_FLOAT ){
    setVectorComponentGridData(d,d_wo,id[igrid],jd[igrid],kd[igrid],0,gc_out);//0:x
    FP3DW->setFuncData(d);
  }else if( FP3DW->GetRealType() == OUTPUT_DOUBLE ){
    setVectorComponentGridData(dd,d_wo,id[igrid],jd[igrid],kd[igrid],0,gc_out);//0:x
    FP3DW->setFuncData(dd);
  }
  if(!FP3DW->WriteFuncData()) std::cout << "error WriteFuncData" << std::endl;
  if( FP3DW->GetRealType() == OUTPUT_FLOAT ){
    setVectorComponentGridData(d,d_wo,id[igrid],jd[igrid],kd[igrid],1,gc_out);//1:y
    FP3DW->setFuncData(d);
  }else if( FP3DW->GetRealType() == OUTPUT_DOUBLE ){
    setVectorComponentGridData(dd,d_wo,id[igrid],jd[igrid],kd[igrid],1,gc_out);//1:y
    FP3DW->setFuncData(dd);
  }
  if(!FP3DW->WriteFuncData()) std::cout << "error WriteFuncData" << std::endl;
  if( FP3DW->GetRealType() == OUTPUT_FLOAT ){
    setVectorComponentGridData(d,d_wo,id[igrid],jd[igrid],kd[igrid],2,gc_out);//2:z
    FP3DW->setFuncData(d);
  }else if( FP3DW->GetRealType() == OUTPUT_DOUBLE ){
    setVectorComponentGridData(dd,d_wo,id[igrid],jd[igrid],kd[igrid],2,gc_out);//2:z
    FP3DW->setFuncData(dd);
  }
  if(!FP3DW->WriteFuncData()) std::cout << "error WriteFuncData" << std::endl;
  //}//igrid loop
  FP3DW->CloseFile();
  dfi_name = "vel_" + C->P3Op.basename_f;

  
  // Tempearture
  if( C->isHeatProblem() ){
    if (C->Unit.File == DIMENSIONAL)
    {
      //U.convArrayTmp2IE(d_ws, size, guide, d_ie, d_bcd, mat_tbl, C->BaseTemp, C->DiffTemp, true, flop);
    }
    else
    {
      //U.convArrayTmp2IE(d_ws, size, guide, d_ie, d_bcd, mat_tbl, C->BaseTemp, C->DiffTemp, false, flop);
    }
    
    fname = "tmp_" + tmp;
    FP3DW->setFileName((dtmp+fname).c_str());
    if(!FP3DW->OpenFile()){
      std::cout << "error OpenFile" << std::endl;
      Exit(0);
    }
    
    FP3DW->WriteNgrid(ngrid);
    FP3DW->WriteFuncBlockData(id,jd,kd,nvar_scalar,ngrid);
    //for(igrid=0;igrid<ngrid;igrid++){ //--->BCMでループが必要になる？
    igrid=0;
    FP3DW->setGridData(id[igrid],jd[igrid],kd[igrid],ngrid);
    FP3DW->setFuncDataNum(nvar_scalar[igrid]);
    if( FP3DW->GetRealType() == OUTPUT_FLOAT ){
      setScalarGridData(d,d_ws,id[igrid],jd[igrid],kd[igrid],gc_out);
      FP3DW->setFuncData(d);
    }else if( FP3DW->GetRealType() == OUTPUT_DOUBLE ){
      setScalarGridData(dd,d_ws,id[igrid],jd[igrid],kd[igrid],gc_out);
      FP3DW->setFuncData(dd);
    }
    if(!FP3DW->WriteFuncData()) std::cout << "error WriteFuncData" << std::endl;
    //}//igrid loop
    FP3DW->CloseFile();
    dfi_name = "tmp_" + C->P3Op.basename_f;

  }
  
  // Total Pressure
  if (C->varState[var_TotalP] == ON ){
    fb_totalp_ (d_p0, size, &guide, d_v, d_p, v00, &flop);
    
    // convert non-dimensional to dimensional, iff file is dimensional
    if (C->Unit.File == DIMENSIONAL)
    {
      U.convArrayTpND2D(d_ws, d_p0, size, guide, C->RefDensity, C->RefVelocity);
    }
    else
    {
      REAL_TYPE* tp;
      tp = d_ws; d_ws = d_p0; d_p0 = tp;
    }
    
    fname = "tp_" + tmp;
    FP3DW->setFileName((dtmp+fname).c_str());
    if(!FP3DW->OpenFile()){
      std::cout << "error OpenFile" << std::endl;
      Exit(0);
    }
    
    FP3DW->WriteNgrid(ngrid);
    FP3DW->WriteFuncBlockData(id,jd,kd,nvar_scalar,ngrid);
    //for(igrid=0;igrid<ngrid;igrid++){ //--->BCMでループが必要になる？
    igrid=0;
    FP3DW->setGridData(id[igrid],jd[igrid],kd[igrid],ngrid);
    FP3DW->setFuncDataNum(nvar_scalar[igrid]);
    if( FP3DW->GetRealType() == OUTPUT_FLOAT ){
      setScalarGridData(d,d_ws,id[igrid],jd[igrid],kd[igrid],gc_out);
      FP3DW->setFuncData(d);
    }else if( FP3DW->GetRealType() == OUTPUT_DOUBLE ){
      setScalarGridData(dd,d_ws,id[igrid],jd[igrid],kd[igrid],gc_out);
      FP3DW->setFuncData(dd);
    }
    if(!FP3DW->WriteFuncData()) std::cout << "error WriteFuncData" << std::endl;
    //}//igrid loop
    FP3DW->CloseFile();
    dfi_name = "tp_" + C->P3Op.basename_f;
  }
  
  // Vorticity
  if (C->varState[var_Vorticity] == ON ){
    rot_v_(d_wv, size, &guide, &deltaX, d_v, d_cdf, v00, &flop);
    
    REAL_TYPE  vz[3];
    vz[0] = vz[1] = vz[2] = 0.0;
    unit_velocity = (C->Unit.File == DIMENSIONAL) ? C->RefVelocity/C->RefLength : 1.0;
    fb_vout_nijk_(d_wo, d_wv, size, &guide, vz, &unit_velocity, &flop);
    
    fname = "vrt_" + tmp;
    FP3DW->setFileName((dtmp+fname).c_str());
    if(!FP3DW->OpenFile()){
      std::cout << "error OpenFile" << std::endl;
      Exit(0);
    }
    
    FP3DW->WriteNgrid(ngrid);
    FP3DW->WriteFuncBlockData(id,jd,kd,nvar_vector,ngrid);
    //for(igrid=0;igrid<ngrid;igrid++){ //--->BCMでループが必要になる？
    igrid=0;
    FP3DW->setGridData(id[igrid],jd[igrid],kd[igrid],ngrid);
    FP3DW->setFuncDataNum(nvar_vector[igrid]);
    if( FP3DW->GetRealType() == OUTPUT_FLOAT ){
      setVectorComponentGridData(d,d_wo,id[igrid],jd[igrid],kd[igrid],0,gc_out);//0:x
      FP3DW->setFuncData(d);
    }else if( FP3DW->GetRealType() == OUTPUT_DOUBLE ){
      setVectorComponentGridData(dd,d_wo,id[igrid],jd[igrid],kd[igrid],0,gc_out);//0:x
      FP3DW->setFuncData(dd);
    }
    if(!FP3DW->WriteFuncData()) std::cout << "error WriteFuncData" << std::endl;
    if( FP3DW->GetRealType() == OUTPUT_FLOAT ){
      setVectorComponentGridData(d,d_wo,id[igrid],jd[igrid],kd[igrid],1,gc_out);//1:y
      FP3DW->setFuncData(d);
    }else if( FP3DW->GetRealType() == OUTPUT_DOUBLE ){
      setVectorComponentGridData(dd,d_wo,id[igrid],jd[igrid],kd[igrid],1,gc_out);//1:y
      FP3DW->setFuncData(dd);
    }
    if(!FP3DW->WriteFuncData()) std::cout << "error WriteFuncData" << std::endl;
    if( FP3DW->GetRealType() == OUTPUT_FLOAT ){
      setVectorComponentGridData(d,d_wo,id[igrid],jd[igrid],kd[igrid],2,gc_out);//2:z
      FP3DW->setFuncData(d);
    }else if( FP3DW->GetRealType() == OUTPUT_DOUBLE ){
      setVectorComponentGridData(dd,d_wo,id[igrid],jd[igrid],kd[igrid],2,gc_out);//2:z
      FP3DW->setFuncData(dd);
    }
    
    //}//igrid loop
    FP3DW->CloseFile();
    dfi_name = "vrt_" + C->P3Op.basename_f;
  }
  
  // 2nd Invariant of Velocity Gradient Tensor
  if (C->varState[var_Qcr] == ON ) {
    
    i2vgt_ (d_p0, size, &guide, &deltaX, d_v, d_cdf, v00, &flop);
    
    U.copyS3D(d_ws, size, guide, d_p0, scale);
    
    fname = "iv2gt_" + tmp;
    FP3DW->setFileName((dtmp+fname).c_str());
    if(!FP3DW->OpenFile()){
      std::cout << "error OpenFile" << std::endl;
      Exit(0);
    }
    
    FP3DW->WriteNgrid(ngrid);
    FP3DW->WriteFuncBlockData(id,jd,kd,nvar_scalar,ngrid);
    //for(igrid=0;igrid<ngrid;igrid++){ //--->BCMでループが必要になる？
    igrid=0;
    FP3DW->setGridData(id[igrid],jd[igrid],kd[igrid],ngrid);
    FP3DW->setFuncDataNum(nvar_scalar[igrid]);
    if( FP3DW->GetRealType() == OUTPUT_FLOAT ){
      setScalarGridData(d,d_ws,id[igrid],jd[igrid],kd[igrid],gc_out);
      FP3DW->setFuncData(d);
    }else if( FP3DW->GetRealType() == OUTPUT_DOUBLE ){
      setScalarGridData(dd,d_ws,id[igrid],jd[igrid],kd[igrid],gc_out);
      FP3DW->setFuncData(dd);
    }
    if(!FP3DW->WriteFuncData()) std::cout << "error WriteFuncData" << std::endl;
    //}//igrid loop
    FP3DW->CloseFile();
    dfi_name = "iv2gt_" + C->P3Op.basename_f;
  }
  
  
  // Helicity
  if (C->varState[var_Helicity] == ON ){
    helicity_(d_p0, size, &guide, &deltaX, d_v, d_cdf, v00, &flop);
    
    U.copyS3D(d_ws, size, guide, d_p0, scale);
    
    fname = "hlt_" + tmp;
    FP3DW->setFileName((dtmp+fname).c_str());
    if(!FP3DW->OpenFile()){
      std::cout << "error OpenFile" << std::endl;
      Exit(0);
    }
    
    FP3DW->WriteNgrid(ngrid);
    FP3DW->WriteFuncBlockData(id,jd,kd,nvar_scalar,ngrid);
    //for(igrid=0;igrid<ngrid;igrid++){ //--->BCMでループが必要になる？
    igrid=0;
    FP3DW->setGridData(id[igrid],jd[igrid],kd[igrid],ngrid);
    FP3DW->setFuncDataNum(nvar_scalar[igrid]);
    if( FP3DW->GetRealType() == OUTPUT_FLOAT ){
      setScalarGridData(d,d_ws,id[igrid],jd[igrid],kd[igrid],gc_out);
      FP3DW->setFuncData(d);
    }else if( FP3DW->GetRealType() == OUTPUT_DOUBLE ){
      setScalarGridData(dd,d_ws,id[igrid],jd[igrid],kd[igrid],gc_out);
      FP3DW->setFuncData(dd);
    }
    if(!FP3DW->WriteFuncData()) std::cout << "error WriteFuncData" << std::endl;
    //}//igrid loop
    FP3DW->CloseFile();
    dfi_name = "hlt_" + C->P3Op.basename_f;
  }
  
  //deallocate
  delete [] d;
  delete [] id;
  delete [] jd;
  delete [] kd;
}


// #################################################################
void Plot3D::function_name()
{  
  //function_nameファイルはかならずformatted形式
  int keep_format=FP3DW->GetFormat();
  FP3DW->setFormat(FORMATTED);
  
  //set filename
  
  // 出力ファイル名
  std::string dtmp = dfi->GenerateDirName(C->FIO.OutDirPath, 0, C->FIO.Slice);
  
  // 出力ディレクトリの作成
  if ( !FBUtility::mkdirs(dtmp) ) {
    if (myRank==0) printf("Error : create directory \"%s\"\n", dtmp.c_str());
    Exit(-1);
  }
  
  std::string tmp = dfi->GenerateFileName(C->P3Op.basename_f, "nam", 0, myRank, false);
  
  //open file
  FP3DW->setFileName((dtmp+tmp).c_str());
  if(!FP3DW->OpenFile()){
    if (myRank==0) printf("Error : error OpenFile\n");
    Exit(0);
  }
  
  // Pressure
  FP3DW->WriteFunctionName("Pressure");
  
  // Velocity
  FP3DW->WriteFunctionName("U-Velocity ; Velocity");
  FP3DW->WriteFunctionName("V-Velocity");
  FP3DW->WriteFunctionName("W-Velocity");
  
  // Tempearture
  if( C->isHeatProblem() ) FP3DW->WriteFunctionName("Tempearture");
  
  // Total Pressure
  if (C->varState[var_TotalP] == ON ) FP3DW->WriteFunctionName("Total_Pressure");
  
  // Vorticity
  if (C->varState[var_Vorticity] == ON ){
    FP3DW->WriteFunctionName("U-Vorticity ; Vorticity");
    FP3DW->WriteFunctionName("V-Vorticity");
    FP3DW->WriteFunctionName("W-Vorticity");
  }
  
  // 2nd Invariant of Velocity Gradient Tensor
  if (C->varState[var_Qcr] == ON ) FP3DW->WriteFunctionName("2nd Invariant of Velocity Gradient Tensor");
  
  // Helicity
  if (C->varState[var_Helicity] == ON ) FP3DW->WriteFunctionName("Helicity");
  
  //close file
  FP3DW->CloseFile();
  
  //reset option
  FP3DW->setFormat(keep_format);
  
}


// #################################################################
void Plot3D::function_name_divide()
{  
  //function_nameファイルはかならずformatted形式
  int keep_format=FP3DW->GetFormat();
  FP3DW->setFormat(FORMATTED);
  
  
  // 出力ファイル名
  std::string fname;
  std::string dtmp = dfi->GenerateDirName(C->FIO.OutDirPath, 0, C->FIO.Slice);
  
  // 出力ディレクトリの作成
  if ( !FBUtility::mkdirs(dtmp) ) {
    if (myRank==0) printf("Error : create directory \"%s\"\n", dtmp.c_str());
    Exit(-1);
  }
  
  std::string tmp = dfi->GenerateFileName(C->P3Op.basename_f, "nam", 0, myRank, false);
  
  // Pressure
  fname = "prs_" + tmp;
  FP3DW->setFileName((dtmp+fname).c_str());
  if(!FP3DW->OpenFile()){
    if (myRank==0) printf("Error : error OpenFile\n");
    Exit(0);
  }
  FP3DW->WriteFunctionName("Pressure");
  FP3DW->CloseFile();
  
  // Velocity
  fname = "vel_" + tmp;
  FP3DW->setFileName((dtmp+fname).c_str());
  if(!FP3DW->OpenFile()){
    if (myRank==0) printf("Error : error OpenFile\n");
    Exit(0);
  }
  FP3DW->WriteFunctionName("U-Velocity ; Velocity");
  FP3DW->WriteFunctionName("V-Velocity");
  FP3DW->WriteFunctionName("W-Velocity");
  FP3DW->CloseFile();
  
  // Tempearture
  if( C->isHeatProblem() ){
    fname = "tmp_" + tmp;
    FP3DW->setFileName((dtmp+fname).c_str());
    if(!FP3DW->OpenFile()){
      if (myRank==0) printf("Error : error OpenFile\n");
      Exit(0);
    }
    FP3DW->WriteFunctionName("Tempearture");
    FP3DW->CloseFile();
  }
  
  // Total Pressure
  if (C->varState[var_TotalP] == ON ){
    fname = "tp_" + tmp;
    FP3DW->setFileName((dtmp+fname).c_str());
    if(!FP3DW->OpenFile()){
      if (myRank==0) printf("Error : error OpenFile\n");
      Exit(0);
    }
    FP3DW->WriteFunctionName("Total_Pressure");
    FP3DW->CloseFile();
  }
  
  // Vorticity
  if (C->varState[var_Vorticity] == ON ){
    fname = "vrt_" + tmp;
    FP3DW->setFileName((dtmp+fname).c_str());
    if(!FP3DW->OpenFile()){
      if (myRank==0) printf("Error : error OpenFile\n");
      Exit(0);
    }
    FP3DW->WriteFunctionName("U-Vorticity ; Vorticity");
    FP3DW->WriteFunctionName("V-Vorticity");
    FP3DW->WriteFunctionName("W-Vorticity");
    FP3DW->CloseFile();
  }
  
  // 2nd Invariant of Velocity Gradient Tensor
  if (C->varState[var_Qcr] == ON ){
    fname = "qcr_" + tmp;
    FP3DW->setFileName((dtmp+fname).c_str());
    if(!FP3DW->OpenFile()){
      if (myRank==0) printf("Error : error OpenFile\n");
      Exit(0);
    }
    FP3DW->WriteFunctionName("2nd Invariant of Velocity Gradient Tensor");
    FP3DW->CloseFile();
  }
  
  
  // Helicity
  if (C->varState[var_Helicity] == ON ){
    fname = "hlt_" + tmp;
    FP3DW->setFileName((dtmp+fname).c_str());
    if(!FP3DW->OpenFile()){
      if (myRank==0) printf("Error : error OpenFile\n");
      Exit(0);
    }
    FP3DW->WriteFunctionName("Helicity");
    FP3DW->CloseFile();
  }
  
  //reset option
  FP3DW->setFormat(keep_format);
  
}


// #################################################################
void Plot3D::fvbnd()
{

}


// #################################################################
void Plot3D::post(const unsigned CurrentStep,
                  const double CurrentTime,
                  REAL_TYPE* v00,
                  const REAL_TYPE* origin,
                  const REAL_TYPE* pitch,
                  int& dfi_mng,
                  double& flop)
{
  if ( C->P3Op.IS_q == ON ) Q(flop);
  
  if ( C->P3Op.IS_DivideFunc == ON )
  {
    if ( C->P3Op.IS_funciton == ON ) function_divide(CurrentStep, CurrentTime, v00, flop);
  }
  else
  {
    if ( C->P3Op.IS_funciton == ON ) function(CurrentStep, CurrentTime, v00, flop);
  }
  
  if ( C->P3Op.IS_xyz == ON )
  {
    if ( FP3DW->IsMoveGrid() )
    {
      xyz(CurrentStep, origin, pitch);
    }
  }
  
  //output dfi file for plot3d
  
  // ステップ数
  int m_step = (int)CurrentStep;
  
  // 時間の次元変換
  REAL_TYPE m_time;
  if (C->Unit.File == DIMENSIONAL)
  {
    m_time = (REAL_TYPE)(CurrentTime * C->Tscale);
  }
  else
  {
    m_time = (REAL_TYPE)CurrentTime;
  }
  
  // 出力モード
  bool mio = (bool)C->FIO.IOmode;
  
  // 最大値と最小値
  REAL_TYPE minmax[2];
  
  minmax[0] = 0.0;
  minmax[1] = 0.0;
  if (myRank==0) if ( !dfi->WriteDFIindex(C->P3Op.basename_f,
                                          C->FIO.OutDirPath,
                                          C->file_fmt_ext,
                                          m_step,
                                          m_time,
                                          dfi_mng,
                                          "ijkn",
                                          1,
                                          minmax,
                                          mio) ) Exit(0);
  
}


// #################################################################
void Plot3D::Q(double& flop)
{
  
}


// #################################################################
void Plot3D::xyz(const unsigned CurrentStep, const REAL_TYPE* origin, const REAL_TYPE* pitch)
{
  //value
  int igrid;
  int ngrid;
  int *id,*jd,*kd;//出力サイズ
  float *x,*y,*z;
  double *dx,*dy,*dz;
  int *iblank,*iblankr;
  
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
  // ガイドセル出力
  //int gc_out = C->GuideOut;
  int gc_out = 0;//plot3dは常にガイドセルを出力しない
  
  //allocate
  ngrid=C->P3Op.ngrid;
  id = new int[ngrid];
  jd = new int[ngrid];
  kd = new int[ngrid];
  
#if 0 //iblank作るときガイドセルを見る
  
  int *idr,*jdr,*kdr;//実際のサイズ
  idr = new int[ngrid];
  jdr = new int[ngrid];
  kdr = new int[ngrid];
  
#endif
  
  //set grid data
  //for(igrid=0;igrid<ngrid;igrid++){ //--->BCMでループが必要になる？
  
  igrid=0;
  id[igrid]=size[0]+1;//+2*gc_out
  jd[igrid]=size[1]+1;//+2*gc_out
  kd[igrid]=size[2]+1;//+2*gc_out
  
#if 0 //iblank作るときガイドセルを見る
  
  idr[igrid]=size[0]+2*gd +1;
  jdr[igrid]=size[1]+2*gd +1;
  kdr[igrid]=size[2]+2*gd +1;
  
#endif
  
  // ステップ数
  int m_step = (int)CurrentStep;
  
  //}//igrid loop
  
  
  // 出力ファイル名
  std::string dtmp = dfi->GenerateDirName(C->FIO.OutDirPath, m_step, C->FIO.Slice);
  
  // 出力ディレクトリの作成
  if ( !FBUtility::mkdirs(dtmp) ) {
    if (myRank==0) printf("Error : create directory \"%s\"\n", dtmp.c_str());
    Exit(-1);
  }
  
  std::string tmp = dfi->GenerateFileName(C->P3Op.basename_g, "xyz", m_step, myRank, true);
  
  //open file
  FP3DW->setFileName((dtmp+tmp).c_str());
  if(!FP3DW->OpenFile()){
    if (myRank==0) printf("Error : error OpenFile\n");
    Exit(0);
  }
  
  //write block data
  FP3DW->WriteNgrid(ngrid);//if multi grid
  for(igrid=0;igrid<ngrid;igrid++){
    FP3DW->WriteBlockData(id[igrid],jd[igrid],kd[igrid]);
  }
  
  //write xyz and iblank data
  
  //for(igrid=0;igrid<ngrid;igrid++){ //--->BCMでループが必要になる？
  
  igrid=0;//igrid=0
  
  //set xyz
  if( FP3DW->GetRealType() == OUTPUT_FLOAT ){
    
    if (!(x = new float[ id[igrid]*jd[igrid]*kd[igrid] ])){
      if (myRank==0)  printf(    "\t>> cannot allocate work area : xyz()\n\n");
      Exit(0);
    }
    if (!(y = new float[ id[igrid]*jd[igrid]*kd[igrid] ])){
      if (myRank==0)  printf(    "\t>> cannot allocate work area : xyz()\n\n");
      Exit(0);
    }
    if (!(z = new float[ id[igrid]*jd[igrid]*kd[igrid] ])){
      if (myRank==0)  printf(    "\t>> cannot allocate work area : xyz()\n\n");
      Exit(0);
    }
    for(int k=0;k<kd[igrid];k++){
      for(int j=0;j<jd[igrid];j++){
        for(int i=0;i<id[igrid];i++){
          size_t ip = _F_IDX_S3D(i+1, j+1, k+1, id[igrid], jd[igrid], kd[igrid], 0);
          x[ip]=(float)origin[0]+(float)pitch[0]*(float)i;//-pitch[0]*(float)gc_out;
          y[ip]=(float)origin[1]+(float)pitch[1]*(float)j;//-pitch[1]*(float)gc_out;
          z[ip]=(float)origin[2]+(float)pitch[2]*(float)k;//-pitch[2]*(float)gc_out;
        }
      }
    }
  }else if( FP3DW->GetRealType() == OUTPUT_DOUBLE ){
    
    if (!(dx = new double[ id[igrid]*jd[igrid]*kd[igrid] ])){
      if (myRank==0)  printf(    "\t>> cannot allocate work area : xyz()\n\n");
      Exit(0);
    }
    if (!(dy = new double[ id[igrid]*jd[igrid]*kd[igrid] ])){
      if (myRank==0)  printf(    "\t>> cannot allocate work area : xyz()\n\n");
      Exit(0);
    }
    if (!(dz = new double[ id[igrid]*jd[igrid]*kd[igrid] ])){
      if (myRank==0)  printf(    "\t>> cannot allocate work area : xyz()\n\n");
      Exit(0);
    }
    for(int k=0;k<kd[igrid];k++){
      for(int j=0;j<jd[igrid];j++){
        for(int i=0;i<id[igrid];i++){
          size_t ip = _F_IDX_S3D(i+1, j+1, k+1, id[igrid], jd[igrid], kd[igrid], 0);
          dx[ip]=(double)origin[0]+(double)pitch[0]*(double)i;//-pitch[0]*(double)gc_out;
          dy[ip]=(double)origin[1]+(double)pitch[1]*(double)j;//-pitch[1]*(double)gc_out;
          dz[ip]=(double)origin[2]+(double)pitch[2]*(double)k;//-pitch[2]*(double)gc_out;
        }
      }
    }
  }
  
  //set iblank
  if(FP3DW->IsIBlankFlag()){
    
#if 0 //iblank作るときガイドセルを見る
    
    //iblank = new int[ id[igrid]*jd[igrid]*kd[igrid] ];
    if (!(iblankr = new int[ idr[igrid]*jdr[igrid]*kdr[igrid] ])){
      if (myRank==0)  printf(    "\t>> cannot allocate work area : xyz()\n\n");
      Exit(0);
    }
    
    setIblankGuide(iblankr,idr[igrid],jdr[igrid],kdr[igrid]);
    
    if(gc_out != 0){// ガイドセル出力する場合
      iblank = iblankr;
    }
    else// ガイドセル出力しない場合
    {
      if (!(iblank = new int[ id[igrid]*jd[igrid]*kd[igrid] ])){
        if (myRank==0)  printf(    "\t>> cannot allocate work area : xyz()\n\n");
        Exit(0);
      }
      for(int kr=gd;kr<kdr[igrid]-gd;kr++){
        for(int jr=gd;jr<jdr[igrid]-gd;jr++){
          for(int ir=gd;ir<idr[igrid]-gd;ir++){
            int i=ir-gd;
            int j=jr-gd;
            int k=kr-gd;
            size_t ipr = _F_IDX_S3D(ir+1, jr+1, kr+1, idr[igrid], jdr[igrid], kdr[igrid], 0);
            size_t ip  = _F_IDX_S3D(i+1, j+1, k+1, id[igrid], jd[igrid], kd[igrid], 0);
            iblank[ip]=iblankr[ipr];
          }
        }
      }
    }
    
#else
    
    if (!(iblank = new int[ id[igrid]*jd[igrid]*kd[igrid] ]) ){
      if (myRank==0)  printf(    "\t>> cannot allocate work area : xyz()\n\n");
      Exit(0);
    }
    setIblank(iblank,id[igrid],jd[igrid],kd[igrid]);
    
#endif
    
  }
  
  //write
  FP3DW->setGridData(id[igrid],jd[igrid],kd[igrid],ngrid);
  if( FP3DW->GetRealType() == OUTPUT_FLOAT ) FP3DW->setXYZData(x,y,z,iblank);
  else if( FP3DW->GetRealType() == OUTPUT_DOUBLE ) FP3DW->setXYZData(dx,dy,dz,iblank);
  if(!FP3DW->WriteXYZData()) std::cout << "error WriteXYZData" << std::endl;
  
  //deallocate
  if( FP3DW->GetRealType() == OUTPUT_FLOAT ){
    delete [] x;
    delete [] y;
    delete [] z;
  }else if( FP3DW->GetRealType() == OUTPUT_DOUBLE ){
    delete [] dx;
    delete [] dy;
    delete [] dz;
  }
  if(FP3DW->IsIBlankFlag()){
#if 0 //iblank作るときガイドセルを見る
    delete [] iblankr;
    if(gc_out == 0) delete [] iblank;
#else
    delete [] iblank;
#endif
  }
  
  //}//igrid loop
  
  //close file
  FP3DW->CloseFile();
  
}


// #################################################################
void Plot3D::printParameters(FILE* fp)
{
  fprintf(fp,"\n\tPLOT3D Options\n");
  fprintf(fp,"\t     grid file prefix         :   %s\n", C->P3Op.basename_g.c_str());
  fprintf(fp,"\t     function file prefix     :   %s\n", C->P3Op.basename_f.c_str());
  fprintf(fp,"\t     grid kind                :   %s\n", (FP3DW->IsGridKind()) ? "multi grid" : "single grid");
  fprintf(fp,"\t     grid mobility            :   %s\n", (FP3DW->IsMoveGrid()) ? "movable" : "immovable");
  fprintf(fp,"\t     state of time            :   %s\n", (FP3DW->IsSteady()) ? "unsteady" : "steady");
  fprintf(fp,"\t     output iblank            :   %s\n", (FP3DW->IsIBlankFlag()) ? "on" : "off");
  if (      FP3DW->GetFormat() == UNFORMATTED ) fprintf(fp,"\t     output format            :   %s\n", "Fortran Unformatted");
  else if ( FP3DW->GetFormat() == FORMATTED   ) fprintf(fp,"\t     output format            :   %s\n", "Fortran Formatted");
  else if ( FP3DW->GetFormat() == C_BINARY    ) fprintf(fp,"\t     output format            :   %s\n", "C Binary");
  fprintf(fp,"\t     output dimention         :   %iD\n", FP3DW->GetDim());
  if (      FP3DW->GetRealType() == OUTPUT_FLOAT  ) fprintf(fp,"\t     output format            :   %s\n", "float");
  else if ( FP3DW->GetRealType() == OUTPUT_DOUBLE ) fprintf(fp,"\t     output format            :   %s\n", "double");
  fprintf(fp,"\t     output xyz file          :   %s\n", (C->P3Op.IS_xyz) ? "on" : "off");
  fprintf(fp,"\t     output q file            :   %s\n", (C->P3Op.IS_q) ? "on" : "off");
  fprintf(fp,"\t     output function file     :   %s\n", (C->P3Op.IS_funciton) ? "on" : "off");
  fprintf(fp,"\t     output funciton name file:   %s\n", (C->P3Op.IS_function_name) ? "on" : "off");
  fprintf(fp,"\t     output fvbnd file        :   %s\n", (C->P3Op.IS_fvbnd) ? "on" : "off");
  fprintf(fp,"\t     function per item        :   %s\n", (C->P3Op.IS_DivideFunc) ? "on" : "off");
}


// #################################################################
void Plot3D::setIblank(int* iblank, int id, int jd, int kd)
{
  //iblank = 1 : 計算グリッド
  //       = 0 : 非計算グリッド
  //       = 2 : 壁面グリッド
  
  int i,j,k;
  size_t ip,mip,ipxm,ipxp,ipym,ipyp,ipzm,ipzp,ipcr;
  size_t m;
  int s, odr;
  
#if 1 //壁面iblankなし
  
  //すべて非計算グリッド（iblank=0）で初期化
  for(k=0; k<kd; k++){
    for(j=0; j<jd; j++){
      for(i=0; i<id; i++){
        ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
        iblank[ip]=0;
      }
    }
  }
  
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
  for (int km=1; km<=kx; km++) {
    for (int jm=1; jm<=jx; jm++) {
      for (int im=1; im<=ix; im++) {
        mip = _F_IDX_S3D(im, jm, km, ix, jx, kx, gd);
        s = d_bcd[mip];
        if(IS_FLUID(s)){//流体であれば
          i=im-1;
          j=jm-1;
          k=km-1;
          size_t ip1 = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
          size_t ip2 = _F_IDX_S3D(i+2, j+1, k+1, id, jd, kd, 0);
          size_t ip3 = _F_IDX_S3D(i+2, j+2, k+1, id, jd, kd, 0);
          size_t ip4 = _F_IDX_S3D(i+1, j+2, k+1, id, jd, kd, 0);
          size_t ip5 = _F_IDX_S3D(i+1, j+1, k+2, id, jd, kd, 0);
          size_t ip6 = _F_IDX_S3D(i+2, j+1, k+2, id, jd, kd, 0);
          size_t ip7 = _F_IDX_S3D(i+2, j+2, k+2, id, jd, kd, 0);
          size_t ip8 = _F_IDX_S3D(i+1, j+2, k+2, id, jd, kd, 0);
          iblank[ip1]=1;
          iblank[ip2]=1;
          iblank[ip3]=1;
          iblank[ip4]=1;
          iblank[ip5]=1;
          iblank[ip6]=1;
          iblank[ip7]=1;
          iblank[ip8]=1;
        }
      }
    }
  }
  
#else //壁面ibrankあり
  
  //すべて計算グリッド（iblank=1）で初期化
  for(k=0;k<kd;k++){
    for(j=0;j<jd;j++){
      for(i=0;i<id;i++){
        ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
        iblank[ip]=1;
      }
    }
  }
  
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
  for (int km=1; km<=kx; km++) {
    for (int jm=1; jm<=jx; jm++) {
      for (int im=1; im<=ix; im++) {
        mip = _F_IDX_S3D(im, jm, km, ix, jx, kx, gd);
        s = d_bcd[mip];
        if(!IS_FLUID(s)){//流体でなければ
          i=im-1;
          j=jm-1;
          k=km-1;
          size_t ip1 = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
          size_t ip2 = _F_IDX_S3D(i+2, j+1, k+1, id, jd, kd, 0);
          size_t ip3 = _F_IDX_S3D(i+2, j+2, k+1, id, jd, kd, 0);
          size_t ip4 = _F_IDX_S3D(i+1, j+2, k+1, id, jd, kd, 0);
          size_t ip5 = _F_IDX_S3D(i+1, j+1, k+2, id, jd, kd, 0);
          size_t ip6 = _F_IDX_S3D(i+2, j+1, k+2, id, jd, kd, 0);
          size_t ip7 = _F_IDX_S3D(i+2, j+2, k+2, id, jd, kd, 0);
          size_t ip8 = _F_IDX_S3D(i+1, j+2, k+2, id, jd, kd, 0);
          iblank[ip1]=0;
          iblank[ip2]=0;
          iblank[ip3]=0;
          iblank[ip4]=0;
          iblank[ip5]=0;
          iblank[ip6]=0;
          iblank[ip7]=0;
          iblank[ip8]=0;
        }
      }
    }
  }
  
  //内部
  for(k=1;k<kd-1;k++){
    for(j=1;j<jd-1;j++){
      for(i=1;i<id-1;i++){
        ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
        if(iblank[ip]!=0) continue;//非計算グリッドでない場合は飛ばす
        
        //非計算グリッドの場合、計算グリッドと非計算グリッドにはさまれているか調べる
        for(int n=-1;n<=1;n++){
          for(int m=-1;m<=1;m++){
            for(int l=-1;l<=1;l++){
              if((l+m+n)==0) continue;
              size_t ipwk = _F_IDX_S3D(i+l+1, j+m+1, k+n+1, id, jd, kd, 0);
              if( iblank[ipwk]==1){
                iblank[ip]=2;
                break;
              }
            }
            if( iblank[ip]==2) break;
          }
          if( iblank[ip]==2) break;
        }
      }
    }
  }
  
  //外面（6面）
  i=0;
  for(k=1;k<kd-1;k++){
    for(j=1;j<jd-1;j++){
      ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
      if(iblank[ip]!=0) continue;//非計算グリッドでない場合は飛ばす
      
      //非計算グリッドの場合、計算グリッドと非計算グリッドにはさまれているか調べる
      for(int n=-1;n<=1;n++){
        for(int m=-1;m<=1;m++){
          int l=0;
          if((l+m+n)==0) continue;
          size_t ipwk = _F_IDX_S3D(i+l+1, j+m+1, k+n+1, id, jd, kd, 0);
          if( iblank[ipwk]==1){
            iblank[ip]=2;
            break;
          }
        }
        if( iblank[ip]==2) break;
      }
    }
  }
  
  i=id-1;
  for(k=1;k<kd-1;k++){
    for(j=1;j<jd-1;j++){
      ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
      if(iblank[ip]!=0) continue;//非計算グリッドでない場合は飛ばす
      
      //非計算グリッドの場合、計算グリッドと非計算グリッドにはさまれているか調べる
      for(int n=-1;n<=1;n++){
        for(int m=-1;m<=1;m++){
          int l=0;
          if((l+m+n)==0) continue;
          size_t ipwk = _F_IDX_S3D(i+l+1, j+m+1, k+n+1, id, jd, kd, 0);
          if( iblank[ipwk]==1){
            iblank[ip]=2;
            break;
          }
        }
        if( iblank[ip]==2) break;
      }
    }
  }
  
  j=0;
  for(k=1;k<kd-1;k++){
    for(i=1;i<id-1;i++){
      ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
      if(iblank[ip]!=0) continue;//非計算グリッドでない場合は飛ばす
      
      //非計算グリッドの場合、計算グリッドと非計算グリッドにはさまれているか調べる
      for(int n=-1;n<=1;n++){
        int m=0;
        for(int l=-1;l<=1;l++){
          if((l+m+n)==0) continue;
          size_t ipwk = _F_IDX_S3D(i+l+1, j+m+1, k+n+1, id, jd, kd, 0);
          if( iblank[ipwk]==1){
            iblank[ip]=2;
            break;
          }
        }
        if( iblank[ip]==2) break;
      }
    }
  }
  
  j=jd-1;
  for(k=1;k<kd-1;k++){
    for(i=1;i<id-1;i++){
      ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
      if(iblank[ip]!=0) continue;//非計算グリッドでない場合は飛ばす
      
      //非計算グリッドの場合、計算グリッドと非計算グリッドにはさまれているか調べる
      for(int n=-1;n<=1;n++){
        int m=0;
        for(int l=-1;l<=1;l++){
          if((l+m+n)==0) continue;
          size_t ipwk = _F_IDX_S3D(i+l+1, j+m+1, k+n+1, id, jd, kd, 0);
          if( iblank[ipwk]==1){
            iblank[ip]=2;
            break;
          }
        }
        if( iblank[ip]==2) break;
      }
    }
  }
  
  k=0;
  for(j=1;j<jd-1;j++){
    for(i=1;i<id-1;i++){
      ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
      if(iblank[ip]!=0) continue;//非計算グリッドでない場合は飛ばす
      
      //非計算グリッドの場合、計算グリッドと非計算グリッドにはさまれているか調べる
      int n=0;
      for(int m=-1;m<=1;m++){
        for(int l=-1;l<=1;l++){
          if((l+m+n)==0) continue;
          size_t ipwk = _F_IDX_S3D(i+l+1, j+m+1, k+n+1, id, jd, kd, 0);
          if( iblank[ipwk]==1){
            iblank[ip]=2;
            break;
          }
        }
        if( iblank[ip]==2) break;
      }
    }
  }
  
  k=id-1;
  for(j=1;j<jd-1;j++){
    for(i=1;i<id-1;i++){
      ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
      if(iblank[ip]!=0) continue;//非計算グリッドでない場合は飛ばす
      
      //非計算グリッドの場合、計算グリッドと非計算グリッドにはさまれているか調べる
      int n=0;
      for(int m=-1;m<=1;m++){
        for(int l=-1;l<=1;l++){
          if((l+m+n)==0) continue;
          size_t ipwk = _F_IDX_S3D(i+l+1, j+m+1, k+n+1, id, jd, kd, 0);
          if( iblank[ipwk]==1){
            iblank[ip]=2;
            break;
          }
        }
        if( iblank[ip]==2) break;
      }
    }
  }
  
  //辺（12辺）
  i=0; j=0;
  for(k=1;k<kd-1;k++){
    ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
    if(iblank[ip]!=0) continue;//非計算グリッドでない場合は飛ばす
    
    //非計算グリッドの場合、計算グリッドと非計算グリッドにはさまれているか調べる
    int l=0; int m=0;
    for(int n=-1;n<=1;n++){
      if((l+m+n)==0) continue;
      size_t ipwk = _F_IDX_S3D(i+l+1, j+m+1, k+n+1, id, jd, kd, 0);
      if( iblank[ipwk]==1){
        iblank[ip]=2;
        break;
      }
    }
  }
  
  i=0; j=jd;
  for(k=1;k<kd-1;k++){
    ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
    if(iblank[ip]!=0) continue;//非計算グリッドでない場合は飛ばす
    
    //非計算グリッドの場合、計算グリッドと非計算グリッドにはさまれているか調べる
    int l=0; int m=0;
    for(int n=-1;n<=1;n++){
      if((l+m+n)==0) continue;
      size_t ipwk = _F_IDX_S3D(i+l+1, j+m+1, k+n+1, id, jd, kd, 0);
      if( iblank[ipwk]==1){
        iblank[ip]=2;
        break;
      }
    }
  }
  
  i=0; k=0;
  for(j=1;j<jd-1;j++){
    ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
    if(iblank[ip]!=0) continue;//非計算グリッドでない場合は飛ばす
    
    //非計算グリッドの場合、計算グリッドと非計算グリッドにはさまれているか調べる
    int l=0; int n=0;
    for(int m=-1;m<=1;m++){
      if((l+m+n)==0) continue;
      size_t ipwk = _F_IDX_S3D(i+l+1, j+m+1, k+n+1, id, jd, kd, 0);
      if( iblank[ipwk]==1){
        iblank[ip]=2;
        break;
      }
    }
  }
  
  i=0; k=kd-1;
  for(j=1;j<jd-1;j++){
    ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
    if(iblank[ip]!=0) continue;//非計算グリッドでない場合は飛ばす
    
    //非計算グリッドの場合、計算グリッドと非計算グリッドにはさまれているか調べる
    int l=0; int n=0;
    for(int m=-1;m<=1;m++){
      if((l+m+n)==0) continue;
      size_t ipwk = _F_IDX_S3D(i+l+1, j+m+1, k+n+1, id, jd, kd, 0);
      if( iblank[ipwk]==1){
        iblank[ip]=2;
        break;
      }
    }
  }
  
  j=0; k=0;
  for(i=1;i<id-1;i++){
    ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
    if(iblank[ip]!=0) continue;//非計算グリッドでない場合は飛ばす
    
    //非計算グリッドの場合、計算グリッドと非計算グリッドにはさまれているか調べる
    int m=0; int n=0;
    for(int l=-1;l<=1;l++){
      if((l+m+n)==0) continue;
      size_t ipwk = _F_IDX_S3D(i+l+1, j+m+1, k+n+1, id, jd, kd, 0);
      if( iblank[ipwk]==1){
        iblank[ip]=2;
        break;
      }
    }
  }
  
  j=0; k=kd-1;
  for(i=1;i<id-1;i++){
    ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
    if(iblank[ip]!=0) continue;//非計算グリッドでない場合は飛ばす
    
    //非計算グリッドの場合、計算グリッドと非計算グリッドにはさまれているか調べる
    int m=0; int n=0;
    for(int l=-1;l<=1;l++){
      if((l+m+n)==0) continue;
      size_t ipwk = _F_IDX_S3D(i+l+1, j+m+1, k+n+1, id, jd, kd, 0);
      if( iblank[ipwk]==1){
        iblank[ip]=2;
        break;
      }
    }
  }
  
  j=jd; k=0;
  for(i=1;i<id-1;i++){
    ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
    if(iblank[ip]!=0) continue;//非計算グリッドでない場合は飛ばす
    
    //非計算グリッドの場合、計算グリッドと非計算グリッドにはさまれているか調べる
    int m=0; int n=0;
    for(int l=-1;l<=1;l++){
      if((l+m+n)==0) continue;
      size_t ipwk = _F_IDX_S3D(i+l+1, j+m+1, k+n+1, id, jd, kd, 0);
      if( iblank[ipwk]==1){
        iblank[ip]=2;
        break;
      }
    }
  }
  
  j=jd-1; k=kd-1;
  for(i=1;i<id-1;i++){
    ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
    if(iblank[ip]!=0) continue;//非計算グリッドでない場合は飛ばす
    
    //非計算グリッドの場合、計算グリッドと非計算グリッドにはさまれているか調べる
    int m=0; int n=0;
    for(int l=-1;l<=1;l++){
      if((l+m+n)==0) continue;
      size_t ipwk = _F_IDX_S3D(i+l+1, j+m+1, k+n+1, id, jd, kd, 0);
      if( iblank[ipwk]==1){
        iblank[ip]=2;
        break;
      }
    }
  }
  
  i=id-1; j=0;
  for(k=1;k<kd-1;k++){
    ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
    if(iblank[ip]!=0) continue;//非計算グリッドでない場合は飛ばす
    
    //非計算グリッドの場合、計算グリッドと非計算グリッドにはさまれているか調べる
    int l=0; int m=0;
    for(int n=-1;n<=1;n++){
      if((l+m+n)==0) continue;
      size_t ipwk = _F_IDX_S3D(i+l+1, j+m+1, k+n+1, id, jd, kd, 0);
      if( iblank[ipwk]==1){
        iblank[ip]=2;
        break;
      }
    }
  }
  
  i=id-1; j=jd-1;
  for(k=1;k<kd-1;k++){
    ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
    if(iblank[ip]!=0) continue;//非計算グリッドでない場合は飛ばす
    
    //非計算グリッドの場合、計算グリッドと非計算グリッドにはさまれているか調べる
    int l=0; int m=0;
    for(int n=-1;n<=1;n++){
      if((l+m+n)==0) continue;
      size_t ipwk = _F_IDX_S3D(i+l+1, j+m+1, k+n+1, id, jd, kd, 0);
      if( iblank[ipwk]==1){
        iblank[ip]=2;
        break;
      }
    }
  }
  
  i=id-1; k=0;
  for(j=1;j<jd-1;j++){
    ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
    if(iblank[ip]!=0) continue;//非計算グリッドでない場合は飛ばす
    
    //非計算グリッドの場合、計算グリッドと非計算グリッドにはさまれているか調べる
    int l=0; int n=0;
    for(int m=-1;m<=1;m++){
      if((l+m+n)==0) continue;
      size_t ipwk = _F_IDX_S3D(i+l+1, j+m+1, k+n+1, id, jd, kd, 0);
      if( iblank[ipwk]==1){
        iblank[ip]=2;
        break;
      }
    }
  }
  
  i=id-1; k=kd-1;
  for(j=1;j<jd-1;j++){
    ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
    if(iblank[ip]!=0) continue;//非計算グリッドでない場合は飛ばす
    
    //非計算グリッドの場合、計算グリッドと非計算グリッドにはさまれているか調べる
    int l=0; int n=0;
    for(int m=-1;m<=1;m++){
      if((l+m+n)==0) continue;
      size_t ipwk = _F_IDX_S3D(i+l+1, j+m+1, k+n+1, id, jd, kd, 0);
      if( iblank[ipwk]==1){
        iblank[ip]=2;
        break;
      }
    }
  }
  
  //頂点（8点）
  i=0; j=0; k=0;
  ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
  if(iblank[ip]==0){ //非計算グリッドの場合
    ipxp = _F_IDX_S3D(i+2, j+1, k+1, id, jd, kd, 0);
    ipyp = _F_IDX_S3D(i+1, j+2, k+1, id, jd, kd, 0);
    ipzp = _F_IDX_S3D(i+1, j+1, k+2, id, jd, kd, 0);
    ipcr = _F_IDX_S3D(i+2, j+2, k+2, id, jd, kd, 0);
    if( iblank[ipxp]==1 || iblank[ipyp]==1 || iblank[ipzp]==1 || iblank[ipcr]==1){
      iblank[ip]=2;
    }
  }
  
  i=id-1;j=0; k=0;
  ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
  if(iblank[ip]==0){ //非計算グリッドの場合
    ipxm = _F_IDX_S3D(i  , j+1, k+1, id, jd, kd, 0);
    ipyp = _F_IDX_S3D(i+1, j+2, k+1, id, jd, kd, 0);
    ipzp = _F_IDX_S3D(i+1, j+1, k+2, id, jd, kd, 0);
    ipcr = _F_IDX_S3D(i  , j+2, k+2, id, jd, kd, 0);
    if( iblank[ipxm]==1 || iblank[ipyp]==1 || iblank[ipzp]==1 || iblank[ipcr]==1){
      iblank[ip]=2;
    }
  }
  
  i=id-1;j=jd-1;k=0;
  ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
  if(iblank[ip]==0){ //非計算グリッドの場合
    ipxm = _F_IDX_S3D(i  , j+1, k+1, id, jd, kd, 0);
    ipym = _F_IDX_S3D(i+1, j  , k+1, id, jd, kd, 0);
    ipzp = _F_IDX_S3D(i+1, j+1, k+2, id, jd, kd, 0);
    ipcr = _F_IDX_S3D(i  , j  , k+2, id, jd, kd, 0);
    if( iblank[ipxm]==1 || iblank[ipym]==1 || iblank[ipzp]==1 || iblank[ipcr]==1){
      iblank[ip]=2;
    }
  }
  
  i=0; j=jd-1;k=0;
  ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
  if(iblank[ip]==0){ //非計算グリッドの場合
    ipxp = _F_IDX_S3D(i+2, j+1, k+1, id, jd, kd, 0);
    ipym = _F_IDX_S3D(i+1, j  , k+1, id, jd, kd, 0);
    ipzp = _F_IDX_S3D(i+1, j+1, k+2, id, jd, kd, 0);
    ipcr = _F_IDX_S3D(i+2, j  , k+2, id, jd, kd, 0);
    if( iblank[ipxp]==1 || iblank[ipym]==1 || iblank[ipzp]==1 || iblank[ipcr]==1){
      iblank[ip]=2;
    }
  }
  
  i=0; j=0; k=kd-1;
  ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
  if(iblank[ip]==0){ //非計算グリッドの場合
    ipxp = _F_IDX_S3D(i+2, j+1, k+1, id, jd, kd, 0);
    ipyp = _F_IDX_S3D(i+1, j+2, k+1, id, jd, kd, 0);
    ipzm = _F_IDX_S3D(i+1, j+1, k  , id, jd, kd, 0);
    ipcr = _F_IDX_S3D(i+2, j+2, k  , id, jd, kd, 0);;
    if( iblank[ipxp]==1 || iblank[ipyp]==1 || iblank[ipzm]==1 || iblank[ipcr]==1){
      iblank[ip]=2;
    }
  }
  
  i=id-1;j=0; k=kd-1;
  ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
  if(iblank[ip]==0){ //非計算グリッドの場合
    ipxm = _F_IDX_S3D(i  , j+1, k+1, id, jd, kd, 0);
    ipyp = _F_IDX_S3D(i+1, j+2, k+1, id, jd, kd, 0);
    ipzm = _F_IDX_S3D(i+1, j+1, k  , id, jd, kd, 0);
    ipcr = _F_IDX_S3D(i  , j+2, k  , id, jd, kd, 0);
    if( iblank[ipxm]==1 || iblank[ipyp]==1 || iblank[ipzm]==1 || iblank[ipcr]==1){
      iblank[ip]=2;
    }
  }
  
  i=id-1;j=jd-1;k=kd-1;
  ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
  if(iblank[ip]==0){ //非計算グリッドの場合
    ipxm = _F_IDX_S3D(i  , j+1, k+1, id, jd, kd, 0);
    ipym = _F_IDX_S3D(i+1, j  , k+1, id, jd, kd, 0);
    ipzm = _F_IDX_S3D(i+1, j+1, k  , id, jd, kd, 0);
    ipcr = _F_IDX_S3D(i  , j  , k  , id, jd, kd, 0);
    if( iblank[ipxm]==1 || iblank[ipym]==1 || iblank[ipzm]==1 || iblank[ipcr]==1){
      iblank[ip]=2;
    }
  }
  
  i=0; j=jd-1;k=kd-1;
  ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
  if(iblank[ip]==0){ //非計算グリッドの場合
    ipxp = _F_IDX_S3D(i+2, j+1, k+1, id, jd, kd, 0);
    ipym = _F_IDX_S3D(i+1, j  , k+1, id, jd, kd, 0);
    ipzm = _F_IDX_S3D(i+1, j+1, k  , id, jd, kd, 0);
    ipcr = _F_IDX_S3D(i  , j  , k  , id, jd, kd, 0);
    if( iblank[ipxp]==1 || iblank[ipym]==1 || iblank[ipzm]==1 || iblank[ipcr]==1){
      iblank[ip]=2;
    }
  }
  
#endif
  
}


// #################################################################
void Plot3D::setIblankGuide(int* iblank, int id, int jd, int kd)
{
  //iblank = 1 : 計算グリッド
  //       = 0 : 非計算グリッド
  //       = 2 : 壁面グリッド
  
  int i,j,k;
  size_t ip,mip,ipxm,ipxp,ipym,ipyp,ipzm,ipzp,ipcr;
  size_t m;
  int s, odr;
  
#if 1 //壁面iblankなし
  
  //すべて非計算グリッド（iblank=0）で初期化
  for(k=0;k<kd;k++){
    for(j=0;j<jd;j++){
      for(i=0;i<id;i++){
        ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
        iblank[ip]=0;
      }
    }
  }
  
  //int id_of_solid; // Geometry Direct Interfaceでテスト的に固定ID=2を与える defined ffv.h
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
  for (int km=1-gd; km<=kx+gd; km++) {
    for (int jm=1-gd; jm<=jx+gd; jm++) {
      for (int im=1-gd; im<=ix+gd; im++) {
        mip = _F_IDX_S3D(im, jm, km, ix, jx, kx, gd);
        s = d_bcd[mip];
        if(IS_FLUID(s)){//流体であれば
          i=im-1+gd;
          j=jm-1+gd;
          k=km-1+gd;
          size_t ip1 = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
          size_t ip2 = _F_IDX_S3D(i+2, j+1, k+1, id, jd, kd, 0);
          size_t ip3 = _F_IDX_S3D(i+2, j+2, k+1, id, jd, kd, 0);
          size_t ip4 = _F_IDX_S3D(i+1, j+2, k+1, id, jd, kd, 0);
          size_t ip5 = _F_IDX_S3D(i+1, j+1, k+2, id, jd, kd, 0);
          size_t ip6 = _F_IDX_S3D(i+2, j+1, k+2, id, jd, kd, 0);
          size_t ip7 = _F_IDX_S3D(i+2, j+2, k+2, id, jd, kd, 0);
          size_t ip8 = _F_IDX_S3D(i+1, j+2, k+2, id, jd, kd, 0);
          iblank[ip1]=1;
          iblank[ip2]=1;
          iblank[ip3]=1;
          iblank[ip4]=1;
          iblank[ip5]=1;
          iblank[ip6]=1;
          iblank[ip7]=1;
          iblank[ip8]=1;
        }
      }
    }
  }
  
#else //壁面ibrankあり
  
  //すべて計算グリッド（iblank=1）で初期化
  for(k=0;k<kd;k++){
    for(j=0;j<jd;j++){
      for(i=0;i<id;i++){
        ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
        iblank[ip]=1;
      }
    }
  }
  
  //int id_of_solid; // Geometry Direct Interfaceでテスト的に固定ID=2を与える defined ffv.h
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
  for (int km=1-gd; km<=kx+gd; km++) {
    for (int jm=1-gd; jm<=jx+gd; jm++) {
      for (int im=1-gd; im<=ix+gd; im++) {
        mip = _F_IDX_S3D(im, jm, km, ix, jx, kx, gd);
        s = d_bcd[mip];
        if(!IS_FLUID(s)){//流体でなければ
          i=im-1+gd;
          j=jm-1+gd;
          k=km-1+gd;
          size_t ip1 = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
          size_t ip2 = _F_IDX_S3D(i+2, j+1, k+1, id, jd, kd, 0);
          size_t ip3 = _F_IDX_S3D(i+2, j+2, k+1, id, jd, kd, 0);
          size_t ip4 = _F_IDX_S3D(i+1, j+2, k+1, id, jd, kd, 0);
          size_t ip5 = _F_IDX_S3D(i+1, j+1, k+2, id, jd, kd, 0);
          size_t ip6 = _F_IDX_S3D(i+2, j+1, k+2, id, jd, kd, 0);
          size_t ip7 = _F_IDX_S3D(i+2, j+2, k+2, id, jd, kd, 0);
          size_t ip8 = _F_IDX_S3D(i+1, j+2, k+2, id, jd, kd, 0);
          iblank[ip1]=0;
          iblank[ip2]=0;
          iblank[ip3]=0;
          iblank[ip4]=0;
          iblank[ip5]=0;
          iblank[ip6]=0;
          iblank[ip7]=0;
          iblank[ip8]=0;
        }
      }
    }
  }
  
  //内部
  for(k=1;k<kd-1;k++){
    for(j=1;j<jd-1;j++){
      for(i=1;i<id-1;i++){
        ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
        if(iblank[ip]!=0) continue;//非計算グリッドでない場合は飛ばす
        
        //非計算グリッドの場合、計算グリッドと非計算グリッドにはさまれているか調べる
        for(int n=-1;n<=1;n++){
          for(int m=-1;m<=1;m++){
            for(int l=-1;l<=1;l++){
              if((l+m+n)==0) continue;
              size_t ipwk = _F_IDX_S3D(i+l+1, j+m+1, k+n+1, id, jd, kd, 0);
              if( iblank[ipwk]==1){
                iblank[ip]=2;
                break;
              }
            }
            if( iblank[ip]==2) break;
          }
          if( iblank[ip]==2) break;
        }
      }
    }
  }
  
  //外面（6面）
  i=0;
  for(k=1;k<kd-1;k++){
    for(j=1;j<jd-1;j++){
      ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
      if(iblank[ip]!=0) continue;//非計算グリッドでない場合は飛ばす
      
      //非計算グリッドの場合、計算グリッドと非計算グリッドにはさまれているか調べる
      for(int n=-1;n<=1;n++){
        for(int m=-1;m<=1;m++){
          int l=0;
          if((l+m+n)==0) continue;
          size_t ipwk = _F_IDX_S3D(i+l+1, j+m+1, k+n+1, id, jd, kd, 0);
          if( iblank[ipwk]==1){
            iblank[ip]=2;
            break;
          }
        }
        if( iblank[ip]==2) break;
      }
    }
  }
  
  i=id-1;
  for(k=1;k<kd-1;k++){
    for(j=1;j<jd-1;j++){
      ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
      if(iblank[ip]!=0) continue;//非計算グリッドでない場合は飛ばす
      
      //非計算グリッドの場合、計算グリッドと非計算グリッドにはさまれているか調べる
      for(int n=-1;n<=1;n++){
        for(int m=-1;m<=1;m++){
          int l=0;
          if((l+m+n)==0) continue;
          size_t ipwk = _F_IDX_S3D(i+l+1, j+m+1, k+n+1, id, jd, kd, 0);
          if( iblank[ipwk]==1){
            iblank[ip]=2;
            break;
          }
        }
        if( iblank[ip]==2) break;
      }
    }
  }
  
  j=0;
  for(k=1;k<kd-1;k++){
    for(i=1;i<id-1;i++){
      ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
      if(iblank[ip]!=0) continue;//非計算グリッドでない場合は飛ばす
      
      //非計算グリッドの場合、計算グリッドと非計算グリッドにはさまれているか調べる
      for(int n=-1;n<=1;n++){
        int m=0;
        for(int l=-1;l<=1;l++){
          if((l+m+n)==0) continue;
          size_t ipwk = _F_IDX_S3D(i+l+1, j+m+1, k+n+1, id, jd, kd, 0);
          if( iblank[ipwk]==1){
            iblank[ip]=2;
            break;
          }
        }
        if( iblank[ip]==2) break;
      }
    }
  }
  
  j=jd-1;
  for(k=1;k<kd-1;k++){
    for(i=1;i<id-1;i++){
      ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
      if(iblank[ip]!=0) continue;//非計算グリッドでない場合は飛ばす
      
      //非計算グリッドの場合、計算グリッドと非計算グリッドにはさまれているか調べる
      for(int n=-1;n<=1;n++){
        int m=0;
        for(int l=-1;l<=1;l++){
          if((l+m+n)==0) continue;
          size_t ipwk = _F_IDX_S3D(i+l+1, j+m+1, k+n+1, id, jd, kd, 0);
          if( iblank[ipwk]==1){
            iblank[ip]=2;
            break;
          }
        }
        if( iblank[ip]==2) break;
      }
    }
  }
  
  k=0;
  for(j=1;j<jd-1;j++){
    for(i=1;i<id-1;i++){
      ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
      if(iblank[ip]!=0) continue;//非計算グリッドでない場合は飛ばす
      
      //非計算グリッドの場合、計算グリッドと非計算グリッドにはさまれているか調べる
      int n=0;
      for(int m=-1;m<=1;m++){
        for(int l=-1;l<=1;l++){
          if((l+m+n)==0) continue;
          size_t ipwk = _F_IDX_S3D(i+l+1, j+m+1, k+n+1, id, jd, kd, 0);
          if( iblank[ipwk]==1){
            iblank[ip]=2;
            break;
          }
        }
        if( iblank[ip]==2) break;
      }
    }
  }
  
  k=id-1;
  for(j=1;j<jd-1;j++){
    for(i=1;i<id-1;i++){
      ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
      if(iblank[ip]!=0) continue;//非計算グリッドでない場合は飛ばす
      
      //非計算グリッドの場合、計算グリッドと非計算グリッドにはさまれているか調べる
      int n=0;
      for(int m=-1;m<=1;m++){
        for(int l=-1;l<=1;l++){
          if((l+m+n)==0) continue;
          size_t ipwk = _F_IDX_S3D(i+l+1, j+m+1, k+n+1, id, jd, kd, 0);
          if( iblank[ipwk]==1){
            iblank[ip]=2;
            break;
          }
        }
        if( iblank[ip]==2) break;
      }
    }
  }
  
  //辺（12辺）
  i=0; j=0;
  for(k=1;k<kd-1;k++){
    ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
    if(iblank[ip]!=0) continue;//非計算グリッドでない場合は飛ばす
    
    //非計算グリッドの場合、計算グリッドと非計算グリッドにはさまれているか調べる
    int l=0; int m=0;
    for(int n=-1;n<=1;n++){
      if((l+m+n)==0) continue;
      size_t ipwk = _F_IDX_S3D(i+l+1, j+m+1, k+n+1, id, jd, kd, 0);
      if( iblank[ipwk]==1){
        iblank[ip]=2;
        break;
      }
    }
  }
  
  i=0; j=jd;
  for(k=1;k<kd-1;k++){
    ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
    if(iblank[ip]!=0) continue;//非計算グリッドでない場合は飛ばす
    
    //非計算グリッドの場合、計算グリッドと非計算グリッドにはさまれているか調べる
    int l=0; int m=0;
    for(int n=-1;n<=1;n++){
      if((l+m+n)==0) continue;
      size_t ipwk = _F_IDX_S3D(i+l+1, j+m+1, k+n+1, id, jd, kd, 0);
      if( iblank[ipwk]==1){
        iblank[ip]=2;
        break;
      }
    }
  }
  
  i=0; k=0;
  for(j=1;j<jd-1;j++){
    ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
    if(iblank[ip]!=0) continue;//非計算グリッドでない場合は飛ばす
    
    //非計算グリッドの場合、計算グリッドと非計算グリッドにはさまれているか調べる
    int l=0; int n=0;
    for(int m=-1;m<=1;m++){
      if((l+m+n)==0) continue;
      size_t ipwk = _F_IDX_S3D(i+l+1, j+m+1, k+n+1, id, jd, kd, 0);
      if( iblank[ipwk]==1){
        iblank[ip]=2;
        break;
      }
    }
  }
  
  i=0; k=kd-1;
  for(j=1;j<jd-1;j++){
    ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
    if(iblank[ip]!=0) continue;//非計算グリッドでない場合は飛ばす
    
    //非計算グリッドの場合、計算グリッドと非計算グリッドにはさまれているか調べる
    int l=0; int n=0;
    for(int m=-1;m<=1;m++){
      if((l+m+n)==0) continue;
      size_t ipwk = _F_IDX_S3D(i+l+1, j+m+1, k+n+1, id, jd, kd, 0);
      if( iblank[ipwk]==1){
        iblank[ip]=2;
        break;
      }
    }
  }
  
  j=0; k=0;
  for(i=1;i<id-1;i++){
    ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
    if(iblank[ip]!=0) continue;//非計算グリッドでない場合は飛ばす
    
    //非計算グリッドの場合、計算グリッドと非計算グリッドにはさまれているか調べる
    int m=0; int n=0;
    for(int l=-1;l<=1;l++){
      if((l+m+n)==0) continue;
      size_t ipwk = _F_IDX_S3D(i+l+1, j+m+1, k+n+1, id, jd, kd, 0);
      if( iblank[ipwk]==1){
        iblank[ip]=2;
        break;
      }
    }
  }
  
  j=0; k=kd-1;
  for(i=1;i<id-1;i++){
    ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
    if(iblank[ip]!=0) continue;//非計算グリッドでない場合は飛ばす
    
    //非計算グリッドの場合、計算グリッドと非計算グリッドにはさまれているか調べる
    int m=0; int n=0;
    for(int l=-1;l<=1;l++){
      if((l+m+n)==0) continue;
      size_t ipwk = _F_IDX_S3D(i+l+1, j+m+1, k+n+1, id, jd, kd, 0);
      if( iblank[ipwk]==1){
        iblank[ip]=2;
        break;
      }
    }
  }
  
  j=jd; k=0;
  for(i=1;i<id-1;i++){
    ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
    if(iblank[ip]!=0) continue;//非計算グリッドでない場合は飛ばす
    
    //非計算グリッドの場合、計算グリッドと非計算グリッドにはさまれているか調べる
    int m=0; int n=0;
    for(int l=-1;l<=1;l++){
      if((l+m+n)==0) continue;
      size_t ipwk = _F_IDX_S3D(i+l+1, j+m+1, k+n+1, id, jd, kd, 0);
      if( iblank[ipwk]==1){
        iblank[ip]=2;
        break;
      }
    }
  }
  
  j=jd-1; k=kd-1;
  for(i=1;i<id-1;i++){
    ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
    if(iblank[ip]!=0) continue;//非計算グリッドでない場合は飛ばす
    
    //非計算グリッドの場合、計算グリッドと非計算グリッドにはさまれているか調べる
    int m=0; int n=0;
    for(int l=-1;l<=1;l++){
      if((l+m+n)==0) continue;
      size_t ipwk = _F_IDX_S3D(i+l+1, j+m+1, k+n+1, id, jd, kd, 0);
      if( iblank[ipwk]==1){
        iblank[ip]=2;
        break;
      }
    }
  }
  
  i=id-1; j=0;
  for(k=1;k<kd-1;k++){
    ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
    if(iblank[ip]!=0) continue;//非計算グリッドでない場合は飛ばす
    
    //非計算グリッドの場合、計算グリッドと非計算グリッドにはさまれているか調べる
    int l=0; int m=0;
    for(int n=-1;n<=1;n++){
      if((l+m+n)==0) continue;
      size_t ipwk = _F_IDX_S3D(i+l+1, j+m+1, k+n+1, id, jd, kd, 0);
      if( iblank[ipwk]==1){
        iblank[ip]=2;
        break;
      }
    }
  }
  
  i=id-1; j=jd-1;
  for(k=1;k<kd-1;k++){
    ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
    if(iblank[ip]!=0) continue;//非計算グリッドでない場合は飛ばす
    
    //非計算グリッドの場合、計算グリッドと非計算グリッドにはさまれているか調べる
    int l=0; int m=0;
    for(int n=-1;n<=1;n++){
      if((l+m+n)==0) continue;
      size_t ipwk = _F_IDX_S3D(i+l+1, j+m+1, k+n+1, id, jd, kd, 0);
      if( iblank[ipwk]==1){
        iblank[ip]=2;
        break;
      }
    }
  }
  
  i=id-1; k=0;
  for(j=1;j<jd-1;j++){
    ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
    if(iblank[ip]!=0) continue;//非計算グリッドでない場合は飛ばす
    
    //非計算グリッドの場合、計算グリッドと非計算グリッドにはさまれているか調べる
    int l=0; int n=0;
    for(int m=-1;m<=1;m++){
      if((l+m+n)==0) continue;
      size_t ipwk = _F_IDX_S3D(i+l+1, j+m+1, k+n+1, id, jd, kd, 0);
      if( iblank[ipwk]==1){
        iblank[ip]=2;
        break;
      }
    }
  }
  
  i=id-1; k=kd-1;
  for(j=1;j<jd-1;j++){
    ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
    if(iblank[ip]!=0) continue;//非計算グリッドでない場合は飛ばす
    
    //非計算グリッドの場合、計算グリッドと非計算グリッドにはさまれているか調べる
    int l=0; int n=0;
    for(int m=-1;m<=1;m++){
      if((l+m+n)==0) continue;
      size_t ipwk = _F_IDX_S3D(i+l+1, j+m+1, k+n+1, id, jd, kd, 0);
      if( iblank[ipwk]==1){
        iblank[ip]=2;
        break;
      }
    }
  }
  
  //頂点（8点）
  i=0; j=0; k=0;
  ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
  if(iblank[ip]==0){ //非計算グリッドの場合
    ipxp = _F_IDX_S3D(i+2, j+1, k+1, id, jd, kd, 0);
    ipyp = _F_IDX_S3D(i+1, j+2, k+1, id, jd, kd, 0);
    ipzp = _F_IDX_S3D(i+1, j+1, k+2, id, jd, kd, 0);
    ipcr = _F_IDX_S3D(i+2, j+2, k+2, id, jd, kd, 0);
    if( iblank[ipxp]==1 || iblank[ipyp]==1 || iblank[ipzp]==1 || iblank[ipcr]==1){
      iblank[ip]=2;
    }
  }
  
  i=id-1;j=0; k=0;
  ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
  if(iblank[ip]==0){ //非計算グリッドの場合
    ipxm = _F_IDX_S3D(i  , j+1, k+1, id, jd, kd, 0);
    ipyp = _F_IDX_S3D(i+1, j+2, k+1, id, jd, kd, 0);
    ipzp = _F_IDX_S3D(i+1, j+1, k+2, id, jd, kd, 0);
    ipcr = _F_IDX_S3D(i  , j+2, k+2, id, jd, kd, 0);
    if( iblank[ipxm]==1 || iblank[ipyp]==1 || iblank[ipzp]==1 || iblank[ipcr]==1){
      iblank[ip]=2;
    }
  }
  
  i=id-1;j=jd-1;k=0;
  ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
  if(iblank[ip]==0){ //非計算グリッドの場合
    ipxm = _F_IDX_S3D(i  , j+1, k+1, id, jd, kd, 0);
    ipym = _F_IDX_S3D(i+1, j  , k+1, id, jd, kd, 0);
    ipzp = _F_IDX_S3D(i+1, j+1, k+2, id, jd, kd, 0);
    ipcr = _F_IDX_S3D(i  , j  , k+2, id, jd, kd, 0);
    if( iblank[ipxm]==1 || iblank[ipym]==1 || iblank[ipzp]==1 || iblank[ipcr]==1){
      iblank[ip]=2;
    }
  }
  
  i=0; j=jd-1;k=0;
  ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
  if(iblank[ip]==0){ //非計算グリッドの場合
    ipxp = _F_IDX_S3D(i+2, j+1, k+1, id, jd, kd, 0);
    ipym = _F_IDX_S3D(i+1, j  , k+1, id, jd, kd, 0);
    ipzp = _F_IDX_S3D(i+1, j+1, k+2, id, jd, kd, 0);
    ipcr = _F_IDX_S3D(i+2, j  , k+2, id, jd, kd, 0);
    if( iblank[ipxp]==1 || iblank[ipym]==1 || iblank[ipzp]==1 || iblank[ipcr]==1){
      iblank[ip]=2;
    }
  }
  
  i=0; j=0; k=kd-1;
  ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
  if(iblank[ip]==0){ //非計算グリッドの場合
    ipxp = _F_IDX_S3D(i+2, j+1, k+1, id, jd, kd, 0);
    ipyp = _F_IDX_S3D(i+1, j+2, k+1, id, jd, kd, 0);
    ipzm = _F_IDX_S3D(i+1, j+1, k  , id, jd, kd, 0);
    ipcr = _F_IDX_S3D(i+2, j+2, k  , id, jd, kd, 0);
    if( iblank[ipxp]==1 || iblank[ipyp]==1 || iblank[ipzm]==1 || iblank[ipcr]==1){
      iblank[ip]=2;
    }
  }
  
  i=id-1;j=0; k=kd-1;
  ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
  if(iblank[ip]==0){ //非計算グリッドの場合
    ipxm = _F_IDX_S3D(i  , j+1, k+1, id, jd, kd, 0);
    ipyp = _F_IDX_S3D(i+1, j+2, k+1, id, jd, kd, 0);
    ipzm = _F_IDX_S3D(i+1, j+1, k  , id, jd, kd, 0);
    ipcr = _F_IDX_S3D(i  , j+2, k  , id, jd, kd, 0);
    if( iblank[ipxm]==1 || iblank[ipyp]==1 || iblank[ipzm]==1 || iblank[ipcr]==1){
      iblank[ip]=2;
    }
  }
  
  i=id-1;j=jd-1;k=kd-1;
  ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
  if(iblank[ip]==0){ //非計算グリッドの場合
    ipxm = _F_IDX_S3D(i  , j+1, k+1, id, jd, kd, 0);
    ipym = _F_IDX_S3D(i+1, j  , k+1, id, jd, kd, 0);
    ipzm = _F_IDX_S3D(i+1, j+1, k  , id, jd, kd, 0);
    ipcr = _F_IDX_S3D(i  , j  , k  , id, jd, kd, 0);
    if( iblank[ipxm]==1 || iblank[ipym]==1 || iblank[ipzm]==1 || iblank[ipcr]==1){
      iblank[ip]=2;
    }
  }
  
  i=0; j=jd-1;k=kd-1;
  ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
  if(iblank[ip]==0){ //非計算グリッドの場合
    ipxp = _F_IDX_S3D(i+2, j+1, k+1, id, jd, kd, 0);
    ipym = _F_IDX_S3D(i+1, j  , k+1, id, jd, kd, 0);
    ipzm = _F_IDX_S3D(i+1, j+1, k  , id, jd, kd, 0);
    ipcr = _F_IDX_S3D(i  , j  , k  , id, jd, kd, 0);
    if( iblank[ipxp]==1 || iblank[ipym]==1 || iblank[ipzm]==1 || iblank[ipcr]==1){
      iblank[ip]=2;
    }
  }
  
#endif
  
}


// #################################################################
void Plot3D::setValuePlot3D()
{
  //set ngrid
  C->P3Op.ngrid=1;
  
  //set nvar
  int nvar=4;//pressure + velocity(3)
  
  if ( C->isHeatProblem() )     nvar++;
  if ( C->varState[var_TotalP] == ON )       nvar++;
  if ( C->varState[var_Vorticity] == ON )      nvar+3;
  if ( C->varState[var_Qcr] == ON )    nvar++;
  if ( C->varState[var_Helicity] == ON ) nvar++;
  
  C->P3Op.nvar=nvar;
}
