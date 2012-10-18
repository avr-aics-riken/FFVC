// #################################################################
//
// FFV : Frontflow / violet Cartesian
//
// Copyright (c) 2012 All right reserved.
//
// Institute of Industrial Science, The University of Tokyo, Japan. 
//
// #################################################################

/** 
 * @file   ffv_plot3d.C
 * @brief  FFV Class
 * @author kero
 */

#include "ffv.h"


//
void FFV::setValuePlot3D()
{
  //set ngrid
  C.P3Op.ngrid=1;

  //set nvar
  int nvar=4;//pressure + velocity(3)
  if( C.isHeatProblem() ) nvar++;
  if (C.Mode.TP == ON ) nvar++;
  if (C.Mode.VRT == ON ) nvar+3;
  if (C.Mode.I2VGT == ON ) nvar++;
  if (C.Mode.Helicity == ON ) nvar++;
  C.P3Op.nvar=nvar;
}

//
void FFV::OutputPlot3D_post(double& flop)
{
  if(C.P3Op.IS_q == ON) OutputPlot3D_q(flop);

  if(C.P3Op.IS_DivideFunc == ON) {
    if(C.P3Op.IS_funciton == ON) OutputPlot3D_function_divide(flop);
  }else{
    if(C.P3Op.IS_funciton == ON) OutputPlot3D_function(flop);
  }

  if(C.P3Op.IS_xyz == ON){
    if(FP3DW.IsMoveGrid()){
      OutputPlot3D_xyz();
    }
  }

}


//
void FFV::OutputPlot3D_xyz()
{
  //value
  int igrid;
  int ngrid;
  int *id,*jd,*kd;//出力サイズ
  //REAL_TYPE *x,*y,*z;
  float *x,*y,*z;
  double *dx,*dy,*dz;
  int *iblank,*iblankr;

  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;

#if 0
  size_t m;
  int s, odr;
  for (int k=1-gd; k<=kx+gd; k++) {  // ガイドセルを含む全領域
    for (int j=1-gd; j<=jx+gd; j++) {
      for (int i=1-gd; i<=ix+gd; i++) {
        m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        s = d_bcd[m];
        Hostonly_ fprintf(stdout, "[%4d %4d %4d], state=%1d: cmp=%3d  ID=%3d mat=%3d  vf=%3d force=%3d\n", 
                          i, j, k, IS_FLUID(s), 
                          DECODE_CMP(s),
                          DECODE_ID(s), 
                          DECODE_MAT(s), 
                          DECODE_VF(s), 
                          (s>>FORCING_BIT)&0x1);
      }
    }
  }
#endif

  // ガイドセル出力
  //int gc_out = C.GuideOut;
  int gc_out = 0;//plot3dは常にガイドセルを出力しない

  //allocate
  ngrid=C.P3Op.ngrid;
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

  // 並列出力モード
  bool pout = ( C.FIO.IO_Output == IO_GATHER ) ? false : true;

  // 出力ファイル名
  std::string tmp;
  std::string dtmp;
  tmp = DFI.Generate_FileName_Free(C.P3Op.basename, "xyz", m_step, myRank, pout);
  //tmp = DFI.Generate_FileName_Free(C.P3Op.basename, "xyz", 0, myRank, pout);
  dtmp = (C.FIO.IO_Mode == Control::io_time_slice) ? DFI.Generate_DirName(C.f_DivDebug, m_step) : C.FIO.IO_DirPath;
  tmp = directory_prefix(dtmp, tmp, C.FIO.IO_Mode, C.Parallelism);

  //open file
  FP3DW.setFileName(tmp.c_str());
  if(!FP3DW.OpenFile()){
    Hostonly_ printf("Error : error OpenFile\n");
    Exit(0);
  }

  //write block data
  FP3DW.WriteNgrid(ngrid);//if multi grid
  for(igrid=0;igrid<ngrid;igrid++){
    FP3DW.WriteBlockData(id[igrid],jd[igrid],kd[igrid]);
  }

  //write xyz and iblank data
  
  //for(igrid=0;igrid<ngrid;igrid++){ //--->BCMでループが必要になる？
  
  igrid=0;//igrid=0
  
  //set xyz
  if( FP3DW.GetRealType() == OUTPUT_FLOAT ){

    if (!(x = new float[ id[igrid]*jd[igrid]*kd[igrid] ])){
      Hostonly_  printf(    "\t>> cannot allocate work area : OutputPlot3D_xyz\n\n");
      Exit(0);
    }
    if (!(y = new float[ id[igrid]*jd[igrid]*kd[igrid] ])){
      Hostonly_  printf(    "\t>> cannot allocate work area : OutputPlot3D_xyz\n\n");
      Exit(0);
    }
    if (!(z = new float[ id[igrid]*jd[igrid]*kd[igrid] ])){
      Hostonly_  printf(    "\t>> cannot allocate work area : OutputPlot3D_xyz\n\n");
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
  }else if( FP3DW.GetRealType() == OUTPUT_DOUBLE ){

    if (!(dx = new double[ id[igrid]*jd[igrid]*kd[igrid] ])){
      Hostonly_  printf(    "\t>> cannot allocate work area : OutputPlot3D_xyz\n\n");
      Exit(0);
    }
    if (!(dy = new double[ id[igrid]*jd[igrid]*kd[igrid] ])){
      Hostonly_  printf(    "\t>> cannot allocate work area : OutputPlot3D_xyz\n\n");
      Exit(0);
    }
    if (!(dz = new double[ id[igrid]*jd[igrid]*kd[igrid] ])){
      Hostonly_  printf(    "\t>> cannot allocate work area : OutputPlot3D_xyz\n\n");
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
  if(FP3DW.IsIBlankFlag()){

#if 0 //iblank作るときガイドセルを見る

    //iblank = new int[ id[igrid]*jd[igrid]*kd[igrid] ];
    if (!(iblankr = new int[ idr[igrid]*jdr[igrid]*kdr[igrid] ])){
      Hostonly_  printf(    "\t>> cannot allocate work area : OutputPlot3D_xyz\n\n");
      Exit(0);
    }

    setIblankGuide(iblankr,idr[igrid],jdr[igrid],kdr[igrid]);

    if(gc_out != 0){// ガイドセル出力する場合
      iblank = iblankr;
    }
    else// ガイドセル出力しない場合
    {
      if (!(iblank = new int[ id[igrid]*jd[igrid]*kd[igrid] ])){
        Hostonly_  printf(    "\t>> cannot allocate work area : OutputPlot3D_xyz\n\n");
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
      Hostonly_  printf(    "\t>> cannot allocate work area : OutputPlot3D_xyz\n\n");
      Exit(0);
    }
    setIblank(iblank,id[igrid],jd[igrid],kd[igrid]);

#endif

  }

  //write
  FP3DW.setGridData(id[igrid],jd[igrid],kd[igrid],ngrid);
  if( FP3DW.GetRealType() == OUTPUT_FLOAT ) FP3DW.setXYZData(x,y,z,iblank);
  else if( FP3DW.GetRealType() == OUTPUT_DOUBLE ) FP3DW.setXYZData(dx,dy,dz,iblank);
  if(!FP3DW.WriteXYZData()) std::cout << "error WriteXYZData" << std::endl;

  //deallocate
  if( FP3DW.GetRealType() == OUTPUT_FLOAT ){
    delete [] x;
    delete [] y;
    delete [] z;
  }else if( FP3DW.GetRealType() == OUTPUT_DOUBLE ){
    delete [] dx;
    delete [] dy;
    delete [] dz;
  }
  if(FP3DW.IsIBlankFlag()){
#if 0 //iblank作るときガイドセルを見る
    delete [] iblankr;
    if(gc_out == 0) delete [] iblank;
#else
    delete [] iblank;
#endif
  }

  //}//igrid loop
  
  //close file
  FP3DW.CloseFile();

}


// 圧縮性流体のための計算結果ファイル（*.q）出力
void FFV::OutputPlot3D_q(double& flop)
{
  
}


// 計算結果ファイル（*.func）出力
void FFV::OutputPlot3D_function(double& flop)
{
  //value
  int igrid;
  int ngrid;
  int *id,*jd,*kd;
  int *nvar;
  int ix,jx,kx;
  //REAL_TYPE *d;
  float *d;
  double *dd;

  //
  REAL_TYPE scale = 1.0;
  //int d_length = (size[0]+2*guide) * (size[1]+2*guide) * (size[2]+2*guide);

  // ガイドセル出力
  //int gc_out = C.GuideOut;
  int gc_out = 0;//plot3dは常にガイドセルを出力しない

  // ステップ数
  int m_step = (int)CurrentStep;
  
  // 時間の次元変換
  REAL_TYPE m_time;
  if (C.Unit.File == DIMENSIONAL)
  {
    m_time = (REAL_TYPE)CurrentTime * C.Tscale;
  }
  else
  {
    m_time = (REAL_TYPE)CurrentTime;
  }
  
  //allocate
  ngrid=C.P3Op.ngrid;
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
  nvar[igrid]=C.P3Op.nvar;
  if(maxid<id[igrid]) maxid=id[igrid];
  if(maxjd<jd[igrid]) maxjd=jd[igrid];
  if(maxkd<kd[igrid]) maxkd=kd[igrid];
  if(maxnvar<nvar[igrid]) maxnvar=nvar[igrid];

  //}//igrid loop

  //allocate workarea
  if(FP3DW.GetFormat() == C_BINARY){//C_BINARYでの出力は項目ごとに書き出し
    if( FP3DW.GetRealType() == OUTPUT_FLOAT ){
      if (!(d = new float[ maxid*maxjd*maxkd ])){
        Hostonly_  printf(    "\t>> cannot allocate work area : OutputPlot3D_function\n\n");
        Exit(0);
      }
    }else if( FP3DW.GetRealType() == OUTPUT_DOUBLE ){
      if (!(dd = new double[ maxid*maxjd*maxkd ])){
        Hostonly_  printf(    "\t>> cannot allocate work area : OutputPlot3D_function\n\n");
        Exit(0);
      }
    }
  }
  else//Fortranによる出力では出力項目すべてを一度に書き出し
  {
    if( FP3DW.GetRealType() == OUTPUT_FLOAT ){
      if (!(d = new float[ maxid*maxjd*maxkd*maxnvar ])){
        Hostonly_  printf(    "\t>> cannot allocate work area : OutputPlot3D_function\n\n");
        Exit(0);
      }
    }else if( FP3DW.GetRealType() == OUTPUT_DOUBLE ){
      if (!(dd = new double[ maxid*maxjd*maxkd*maxnvar ])){
        Hostonly_  printf(    "\t>> cannot allocate work area : OutputPlot3D_function\n\n");
        Exit(0);
      }
    }
  }

  // 並列出力モード
  bool pout = ( C.FIO.IO_Output == IO_GATHER ) ? false : true;

  // 出力ファイル名
  std::string tmp;
  std::string dtmp;
  tmp = DFI.Generate_FileName_Free(C.P3Op.basename, "func", m_step, myRank, pout);
  dtmp = (C.FIO.IO_Mode == Control::io_time_slice) ? DFI.Generate_DirName(C.f_DivDebug, m_step) : C.FIO.IO_DirPath;
  tmp = directory_prefix(dtmp, tmp, C.FIO.IO_Mode, C.Parallelism);

  //open file
  FP3DW.setFileName(tmp.c_str());
  if(!FP3DW.OpenFile()){
    std::cout << "error OpenFile" << std::endl;
    Exit(0);
  }

  //write block data
  FP3DW.WriteNgrid(ngrid);
  FP3DW.WriteFuncBlockData(id,jd,kd,nvar,ngrid);
  
  //write function data

  //for(igrid=0;igrid<ngrid;igrid++){ //--->BCMでループが必要になる？
  igrid=0;
  
  size_t out_length = (size_t)(id[igrid]*jd[igrid]*kd[igrid]);

  //set grid data
  FP3DW.setGridData(id[igrid],jd[igrid],kd[igrid],ngrid);
  FP3DW.setFuncDataNum(nvar[igrid]);

  //output start
  int ivar=0;

  // Pressure
  //if (C.Unit.File == DIMENSIONAL)
  //{
  //  REAL_TYPE bp = ( C.Unit.Prs == Unit_Absolute ) ? C.BasePrs : 0.0;
  //  fb_prs_nd2d_(d_ws, d_p, &d_length, &bp, &C.RefDensity, &C.RefVelocity, &scale, &flop);
  //}
  //else
  //{
  //  fb_xcopy_(d_ws, d_p, &d_length, &scale, &flop);
  //}
  if (C.Unit.File == DIMENSIONAL) 
  {
    REAL_TYPE bp = ( C.Unit.Prs == Unit_Absolute ) ? C.BasePrs : 0.0;
    U.prs_array_ND2D(d_ws, size, guide, d_p, bp, C.RefDensity, C.RefVelocity, scale, flop);
  }
  else 
  {
    U.xcopy(d_ws, size, guide, d_p, scale, kind_scalar, flop);
  }

  if(FP3DW.GetFormat() == C_BINARY){//C_BINARYでの出力は項目ごとに書き出し
    if( FP3DW.GetRealType() == OUTPUT_FLOAT ){
      setScalarGridData(d,d_ws,id[igrid],jd[igrid],kd[igrid]);
      FP3DW.setFuncData(d);
    }else if( FP3DW.GetRealType() == OUTPUT_DOUBLE ){
      setScalarGridData(dd,d_ws,id[igrid],jd[igrid],kd[igrid]);
      FP3DW.setFuncData(dd);
    }
    if(!FP3DW.WriteFuncData()) std::cout << "error WriteFuncData" << std::endl;
  }
  else//Fortranによる出力では出力項目すべてを一度に書き出し
  {
    if( FP3DW.GetRealType() == OUTPUT_FLOAT ){
      setScalarGridData(&d[out_length*ivar],d_ws,id[igrid],jd[igrid],kd[igrid]);
    }else if( FP3DW.GetRealType() == OUTPUT_DOUBLE ){
      setScalarGridData(&dd[out_length*ivar],d_ws,id[igrid],jd[igrid],kd[igrid]);
    }
  }
  ivar++;
  
  // Velocity
  REAL_TYPE unit_velocity = (C.Unit.File == DIMENSIONAL) ? C.RefVelocity : 1.0;
  fb_shift_refv_out_(d_wo, d_v, size, &guide, v00, &scale, &unit_velocity, &flop);

  if(FP3DW.GetFormat() == C_BINARY){//C_BINARYでの出力は項目ごとに書き出し
    if( FP3DW.GetRealType() == OUTPUT_FLOAT ){
      setVectorComponentGridData(d,d_wo,id[igrid],jd[igrid],kd[igrid],0);//0:x
      FP3DW.setFuncData(d);
    }else if( FP3DW.GetRealType() == OUTPUT_DOUBLE ){
      setVectorComponentGridData(dd,d_wo,id[igrid],jd[igrid],kd[igrid],0);//0:x
      FP3DW.setFuncData(dd);
    }
    if(!FP3DW.WriteFuncData()) std::cout << "error WriteFuncData" << std::endl;
    if( FP3DW.GetRealType() == OUTPUT_FLOAT ){
      setVectorComponentGridData(d,d_wo,id[igrid],jd[igrid],kd[igrid],1);//1:y
      FP3DW.setFuncData(d);
    }else if( FP3DW.GetRealType() == OUTPUT_DOUBLE ){
      setVectorComponentGridData(dd,d_wo,id[igrid],jd[igrid],kd[igrid],1);//1:y
      FP3DW.setFuncData(dd);
    }
    if(!FP3DW.WriteFuncData()) std::cout << "error WriteFuncData" << std::endl;
    if( FP3DW.GetRealType() == OUTPUT_FLOAT ){
      setVectorComponentGridData(d,d_wo,id[igrid],jd[igrid],kd[igrid],2);//2:z
      FP3DW.setFuncData(d);
    }else if( FP3DW.GetRealType() == OUTPUT_DOUBLE ){
      setVectorComponentGridData(dd,d_wo,id[igrid],jd[igrid],kd[igrid],2);//2:z
      FP3DW.setFuncData(dd);
    }
    if(!FP3DW.WriteFuncData()) std::cout << "error WriteFuncData" << std::endl;
  }
  else//Fortranによる出力では出力項目すべてを一度に書き出し
  {
    if( FP3DW.GetRealType() == OUTPUT_FLOAT ){
      setVectorGridData(&d[out_length*ivar],d_wo,id[igrid],jd[igrid],kd[igrid]);
    }else if( FP3DW.GetRealType() == OUTPUT_DOUBLE ){
      setVectorGridData(&dd[out_length*ivar],d_wo,id[igrid],jd[igrid],kd[igrid]);
    }
  }
  ivar=ivar+3;

  // Tempearture
  if( C.isHeatProblem() )
  {
    //if (C.Unit.File == DIMENSIONAL)
    //{
    //  REAL_TYPE klv = ( C.Unit.Temp == Unit_KELVIN ) ? 0.0 : KELVIN;
    //  fb_tmp_nd2d_(d_ws, d_t, &d_length, &C.BaseTemp, &C.DiffTemp, &klv, &scale, &flop);
    //}
    //else
    //{
    //  fb_xcopy_(d_ws, d_t, &d_length, &scale, &flop);
    //}
    if (C.Unit.File == DIMENSIONAL) 
    {
      REAL_TYPE klv = ( C.Unit.Temp == Unit_KELVIN ) ? 0.0 : KELVIN;
      U.tmp_array_ND2D(d_ws, size, guide, d_t, C.BaseTemp, C.DiffTemp, klv, scale, flop);
    }
    else 
    {
      U.xcopy(d_ws, size, guide, d_t, scale, kind_scalar, flop);
    }

    if(FP3DW.GetFormat() == C_BINARY){//C_BINARYでの出力は項目ごとに書き出し
      if( FP3DW.GetRealType() == OUTPUT_FLOAT ){
        setScalarGridData(d,d_ws,id[igrid],jd[igrid],kd[igrid]);
        FP3DW.setFuncData(d);
      }else if( FP3DW.GetRealType() == OUTPUT_DOUBLE ){
        setScalarGridData(dd,d_ws,id[igrid],jd[igrid],kd[igrid]);
        FP3DW.setFuncData(dd);
      }
      if(!FP3DW.WriteFuncData()) std::cout << "error WriteFuncData" << std::endl;
    }
    else//Fortranによる出力では出力項目すべてを一度に書き出し
    {
      if( FP3DW.GetRealType() == OUTPUT_FLOAT ){
        setScalarGridData(&d[out_length*ivar],d_ws,id[igrid],jd[igrid],kd[igrid]);
      }else if( FP3DW.GetRealType() == OUTPUT_DOUBLE ){
        setScalarGridData(&dd[out_length*ivar],d_ws,id[igrid],jd[igrid],kd[igrid]);
      }
    }
    ivar++;

  }

  // Total Pressure
  if (C.Mode.TP == ON ) {
    fb_totalp_ (d_p0, size, &guide, d_v, d_p, v00, &flop);
    
    // convert non-dimensional to dimensional, iff file is dimensional
    if (C.Unit.File == DIMENSIONAL)
    {
      F.cnv_TP_ND2D(d_ws, d_p0, size, guide, C.RefDensity, C.RefVelocity, flop);
    }
    else
    {
      REAL_TYPE* tp;
      tp = d_ws; d_ws = d_p0; d_p0 = tp;
    }
    
    if(FP3DW.GetFormat() == C_BINARY){//C_BINARYでの出力は項目ごとに書き出し
      if( FP3DW.GetRealType() == OUTPUT_FLOAT ){
        setScalarGridData(d,d_ws,id[igrid],jd[igrid],kd[igrid]);
        FP3DW.setFuncData(d);
      }else if( FP3DW.GetRealType() == OUTPUT_DOUBLE ){
        setScalarGridData(dd,d_ws,id[igrid],jd[igrid],kd[igrid]);
        FP3DW.setFuncData(dd);
      }
      if(!FP3DW.WriteFuncData()) std::cout << "error WriteFuncData" << std::endl;
    }
    else//Fortranによる出力では出力項目すべてを一度に書き出し
    {
      if( FP3DW.GetRealType() == OUTPUT_FLOAT ){
        setScalarGridData(&d[out_length*ivar],d_ws,id[igrid],jd[igrid],kd[igrid]);
      }else if( FP3DW.GetRealType() == OUTPUT_DOUBLE ){
        setScalarGridData(&dd[out_length*ivar],d_ws,id[igrid],jd[igrid],kd[igrid]);
      }
    }
    ivar++;
  }

  // Vorticity
  if (C.Mode.VRT == ON )
  {
    rot_v_(d_wv, size, &guide, &deltaX, d_v, d_bcv, v00, &flop);
    
    REAL_TYPE  vz[3];
    vz[0] = vz[1] = vz[2] = 0.0;
    unit_velocity = (C.Unit.File == DIMENSIONAL) ? C.RefVelocity/C.RefLength : 1.0;
    fb_shift_refv_out_(d_wo, d_wv, size, &guide, vz, &scale, &unit_velocity, &flop);

    if(FP3DW.GetFormat() == C_BINARY){//C_BINARYでの出力は項目ごとに書き出し
      if( FP3DW.GetRealType() == OUTPUT_FLOAT ){
        setVectorComponentGridData(d,d_wo,id[igrid],jd[igrid],kd[igrid],0);//0:x
        FP3DW.setFuncData(d);
      }else if( FP3DW.GetRealType() == OUTPUT_DOUBLE ){
        setVectorComponentGridData(dd,d_wo,id[igrid],jd[igrid],kd[igrid],0);//0:x
        FP3DW.setFuncData(dd);
      }
      if(!FP3DW.WriteFuncData()) std::cout << "error WriteFuncData" << std::endl;
      if( FP3DW.GetRealType() == OUTPUT_FLOAT ){
        setVectorComponentGridData(d,d_wo,id[igrid],jd[igrid],kd[igrid],1);//1:y
        FP3DW.setFuncData(d);
      }else if( FP3DW.GetRealType() == OUTPUT_DOUBLE ){
        setVectorComponentGridData(dd,d_wo,id[igrid],jd[igrid],kd[igrid],1);//1:y
        FP3DW.setFuncData(dd);
      }
      if(!FP3DW.WriteFuncData()) std::cout << "error WriteFuncData" << std::endl;
      if( FP3DW.GetRealType() == OUTPUT_FLOAT ){
        setVectorComponentGridData(d,d_wo,id[igrid],jd[igrid],kd[igrid],2);//2:z
        FP3DW.setFuncData(d);
      }else if( FP3DW.GetRealType() == OUTPUT_DOUBLE ){
        setVectorComponentGridData(dd,d_wo,id[igrid],jd[igrid],kd[igrid],2);//2:z
        FP3DW.setFuncData(dd);
      }
      if(!FP3DW.WriteFuncData()) std::cout << "error WriteFuncData" << std::endl;
    }
    else//Fortranによる出力では出力項目すべてを一度に書き出し
    {
      if( FP3DW.GetRealType() == OUTPUT_FLOAT ){
        setVectorGridData(&d[out_length*ivar],d_wo,id[igrid],jd[igrid],kd[igrid]);
      }else if( FP3DW.GetRealType() == OUTPUT_DOUBLE ){
        setVectorGridData(&dd[out_length*ivar],d_wo,id[igrid],jd[igrid],kd[igrid]);
      }
    }
    ivar=ivar+3;
  }

  // 2nd Invariant of Velocity Gradient Tensor
  if (C.Mode.I2VGT == ON ) {
    i2vgt_ (d_p0, size, &guide, &deltaX, d_v, d_bcv, v00, &flop);
    
    // 無次元で出力
    //d_length = (size[0]+2*guide) * (size[1]+2*guide) * (size[2]+2*guide);
    //fb_xcopy_(d_ws, d_p0, &d_length, &scale, &flop);
    U.xcopy(d_ws, size, guide, d_p0, scale, kind_scalar, flop);

    if(FP3DW.GetFormat() == C_BINARY){//C_BINARYでの出力は項目ごとに書き出し
      if( FP3DW.GetRealType() == OUTPUT_FLOAT ){
        setScalarGridData(d,d_ws,id[igrid],jd[igrid],kd[igrid]);
        FP3DW.setFuncData(d);
      }else if( FP3DW.GetRealType() == OUTPUT_DOUBLE ){
        setScalarGridData(dd,d_ws,id[igrid],jd[igrid],kd[igrid]);
        FP3DW.setFuncData(dd);
      }
      if(!FP3DW.WriteFuncData()) std::cout << "error WriteFuncData" << std::endl;
    }
    else//Fortranによる出力では出力項目すべてを一度に書き出し
    {
      if( FP3DW.GetRealType() == OUTPUT_FLOAT ){
        setScalarGridData(&d[out_length*ivar],d_ws,id[igrid],jd[igrid],kd[igrid]);
      }else if( FP3DW.GetRealType() == OUTPUT_DOUBLE ){
        setScalarGridData(&dd[out_length*ivar],d_ws,id[igrid],jd[igrid],kd[igrid]);
      }
    }
    ivar++;
  }
  
  // Helicity
  if (C.Mode.Helicity == ON )
  {
    helicity_(d_p0, size, &guide, &deltaX, d_v, d_bcv, v00, &flop);
    
    // 無次元で出力
    //d_length = (size[0]+2*guide) * (size[1]+2*guide) * (size[2]+2*guide);
    //fb_xcopy_(d_ws, d_p0, &d_length, &scale, &flop);
    U.xcopy(d_ws, size, guide, d_p0, scale, kind_scalar, flop);

    if(FP3DW.GetFormat() == C_BINARY){//C_BINARYでの出力は項目ごとに書き出し
      if( FP3DW.GetRealType() == OUTPUT_FLOAT ){
        setScalarGridData(d,d_ws,id[igrid],jd[igrid],kd[igrid]);
        FP3DW.setFuncData(d);
      }else if( FP3DW.GetRealType() == OUTPUT_DOUBLE ){
        setScalarGridData(dd,d_ws,id[igrid],jd[igrid],kd[igrid]);
        FP3DW.setFuncData(dd);
      }
      if(!FP3DW.WriteFuncData()) std::cout << "error WriteFuncData" << std::endl;
    }
    else//Fortranによる出力では出力項目すべてを一度に書き出し
    {
      if( FP3DW.GetRealType() == OUTPUT_FLOAT ){
        setScalarGridData(&d[out_length*ivar],d_ws,id[igrid],jd[igrid],kd[igrid]);
      }else if( FP3DW.GetRealType() == OUTPUT_DOUBLE ){
        setScalarGridData(&dd[out_length*ivar],d_ws,id[igrid],jd[igrid],kd[igrid]);
      }
    }
    ivar++;
  }
  
  //write all
  if(FP3DW.GetFormat() != C_BINARY){//C_BINARY以外は出力項目すべてを一度に書き出し
    if( FP3DW.GetRealType() == OUTPUT_FLOAT ){
      FP3DW.setFuncData(d);
    }else if( FP3DW.GetRealType() == OUTPUT_DOUBLE ){
      FP3DW.setFuncData(dd);
    }
    if(!FP3DW.WriteFuncData()) std::cout << "error WriteFuncData" << std::endl;
  }

  //}//igrid loop
  
  //close file
  FP3DW.CloseFile();
  
  //deallocate
  if( FP3DW.GetRealType() == OUTPUT_FLOAT ){
    delete [] d;
  }else if( FP3DW.GetRealType() == OUTPUT_DOUBLE ){
    delete [] dd;
  }
  delete [] id;
  delete [] jd;
  delete [] kd;
  
  //dfiファイルの出力
  Hostonly_ if ( !DFI.Write_DFI_File(C.P3Op.basename, m_step, (double)m_time, dfi_mng[var_Plot3D], pout) ) Exit(0);

}


// 項目別計算結果ファイル（*.func）出力
void FFV::OutputPlot3D_function_divide(double& flop)
{
  //value
  int igrid;
  int ngrid;
  int *id,*jd,*kd;
  int *nvar;
  int ix,jx,kx;
  //REAL_TYPE *d;
  float *d;
  double *dd;
  int *nvar_scalar;
  int *nvar_vector;

  //
  REAL_TYPE scale = 1.0;
  //int d_length = (size[0]+2*guide) * (size[1]+2*guide) * (size[2]+2*guide);

  // ガイドセル出力
  //int gc_out = C.GuideOut;
  int gc_out = 0;//plot3dは常にガイドセルを出力しない

  // ステップ数
  int m_step = (int)CurrentStep;
  
  // 時間の次元変換
  REAL_TYPE m_time;
  if (C.Unit.File == DIMENSIONAL)
  {
    m_time = (REAL_TYPE)CurrentTime * C.Tscale;
  }
  else
  {
    m_time = (REAL_TYPE)CurrentTime;
  }
  
  //allocate
  ngrid=C.P3Op.ngrid;
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
  nvar[igrid]=C.P3Op.nvar;
  nvar_scalar[igrid]=1;
  nvar_vector[igrid]=3;

  if(maxid<id[igrid]) maxid=id[igrid];
  if(maxjd<jd[igrid]) maxjd=jd[igrid];
  if(maxkd<kd[igrid]) maxkd=kd[igrid];

  //}//igrid loop

  //allocate workarea
  if( FP3DW.GetRealType() == OUTPUT_FLOAT ){
    if (!(d = new float[ maxid*maxjd*maxkd ])){
      Hostonly_  printf(    "\t>> cannot allocate work area : OutputPlot3D_function\n\n");
      Exit(0);
    }
  }else if( FP3DW.GetRealType() == OUTPUT_DOUBLE ){
    if (!(dd = new double[ maxid*maxjd*maxkd ])){
      Hostonly_  printf(    "\t>> cannot allocate work area : OutputPlot3D_function\n\n");
      Exit(0);
    }
  }

  // 並列出力モード
  bool pout = ( C.FIO.IO_Output == IO_GATHER ) ? false : true;

  // 出力ファイル名
  std::string tmp,fname,dfi_name;
  std::string dtmp;
  tmp = DFI.Generate_FileName_Free(C.P3Op.basename, "func", m_step, myRank, pout);


  // Pressure
  //if (C.Unit.File == DIMENSIONAL)
  //{
  //  REAL_TYPE bp = ( C.Unit.Prs == Unit_Absolute ) ? C.BasePrs : 0.0;
  //  fb_prs_nd2d_(d_ws, d_p, &d_length, &bp, &C.RefDensity, &C.RefVelocity, &scale, &flop);
  //}
  //else
  //{
  //  fb_xcopy_(d_ws, d_p, &d_length, &scale, &flop);
  //}
  if (C.Unit.File == DIMENSIONAL) 
  {
    REAL_TYPE bp = ( C.Unit.Prs == Unit_Absolute ) ? C.BasePrs : 0.0;
    U.prs_array_ND2D(d_ws, size, guide, d_p, bp, C.RefDensity, C.RefVelocity, scale, flop);
  }
  else 
  {
    U.xcopy(d_ws, size, guide, d_p, scale, kind_scalar, flop);
  }

  fname = "prs_" + tmp;
  dtmp = (C.FIO.IO_Mode == Control::io_time_slice) ? DFI.Generate_DirName(C.f_DivDebug, m_step) : C.FIO.IO_DirPath;
  fname = directory_prefix(dtmp, fname, C.FIO.IO_Mode, C.Parallelism);
  FP3DW.setFileName(fname.c_str());
  if(!FP3DW.OpenFile()){
    std::cout << "error OpenFile" << std::endl;
    Exit(0);
  }

  FP3DW.WriteNgrid(ngrid);
  FP3DW.WriteFuncBlockData(id,jd,kd,nvar_scalar,ngrid);

  //for(igrid=0;igrid<ngrid;igrid++){ //--->BCMでループが必要になる？
  igrid=0;
  FP3DW.setGridData(id[igrid],jd[igrid],kd[igrid],ngrid);
  FP3DW.setFuncDataNum(nvar_scalar[igrid]);
  if( FP3DW.GetRealType() == OUTPUT_FLOAT ){
    setScalarGridData(d,d_ws,id[igrid],jd[igrid],kd[igrid]);
    FP3DW.setFuncData(d);
  }else if( FP3DW.GetRealType() == OUTPUT_DOUBLE ){
    setScalarGridData(dd,d_ws,id[igrid],jd[igrid],kd[igrid]);
    FP3DW.setFuncData(dd);
  }
  if(!FP3DW.WriteFuncData()) std::cout << "error WriteFuncData" << std::endl;
  //}//igrid loop
  FP3DW.CloseFile();
  dfi_name = "prs_" + C.P3Op.basename;
  Hostonly_ if ( !DFI.Write_DFI_File(dfi_name, m_step, (double)m_time, dfi_mng[var_Plot3D], pout) ) Exit(0);

  // Velocity
  REAL_TYPE unit_velocity = (C.Unit.File == DIMENSIONAL) ? C.RefVelocity : 1.0;
  fb_shift_refv_out_(d_wo, d_v, size, &guide, v00, &scale, &unit_velocity, &flop);

  fname = "vel_" + tmp;
  dtmp = (C.FIO.IO_Mode == Control::io_time_slice) ? DFI.Generate_DirName(C.f_DivDebug, m_step) : C.FIO.IO_DirPath;
  fname = directory_prefix(dtmp, fname, C.FIO.IO_Mode, C.Parallelism);
  FP3DW.setFileName(fname.c_str());
  if(!FP3DW.OpenFile()){
    std::cout << "error OpenFile" << std::endl;
    Exit(0);
  }

  FP3DW.WriteNgrid(ngrid);
  FP3DW.WriteFuncBlockData(id,jd,kd,nvar_vector,ngrid);

  //for(igrid=0;igrid<ngrid;igrid++){ //--->BCMでループが必要になる？
  igrid=0;
  FP3DW.setGridData(id[igrid],jd[igrid],kd[igrid],ngrid);
  FP3DW.setFuncDataNum(nvar_vector[igrid]);
  if( FP3DW.GetRealType() == OUTPUT_FLOAT ){
    setVectorComponentGridData(d,d_wo,id[igrid],jd[igrid],kd[igrid],0);//0:x
    FP3DW.setFuncData(d);
  }else if( FP3DW.GetRealType() == OUTPUT_DOUBLE ){
    setVectorComponentGridData(dd,d_wo,id[igrid],jd[igrid],kd[igrid],0);//0:x
    FP3DW.setFuncData(dd);
  }
  if(!FP3DW.WriteFuncData()) std::cout << "error WriteFuncData" << std::endl;
  if( FP3DW.GetRealType() == OUTPUT_FLOAT ){
    setVectorComponentGridData(d,d_wo,id[igrid],jd[igrid],kd[igrid],1);//1:y
    FP3DW.setFuncData(d);
  }else if( FP3DW.GetRealType() == OUTPUT_DOUBLE ){
    setVectorComponentGridData(dd,d_wo,id[igrid],jd[igrid],kd[igrid],1);//1:y
    FP3DW.setFuncData(dd);
  }
  if(!FP3DW.WriteFuncData()) std::cout << "error WriteFuncData" << std::endl;
  if( FP3DW.GetRealType() == OUTPUT_FLOAT ){
    setVectorComponentGridData(d,d_wo,id[igrid],jd[igrid],kd[igrid],2);//2:z
    FP3DW.setFuncData(d);
  }else if( FP3DW.GetRealType() == OUTPUT_DOUBLE ){
    setVectorComponentGridData(dd,d_wo,id[igrid],jd[igrid],kd[igrid],2);//2:z
    FP3DW.setFuncData(dd);
  }
  if(!FP3DW.WriteFuncData()) std::cout << "error WriteFuncData" << std::endl;
  //}//igrid loop
  FP3DW.CloseFile();
  dfi_name = "vel_" + C.P3Op.basename;
  Hostonly_ if ( !DFI.Write_DFI_File(dfi_name, m_step, (double)m_time, dfi_mng[var_Plot3D], pout) ) Exit(0);

  // Tempearture
  if( C.isHeatProblem() ){
    //if (C.Unit.File == DIMENSIONAL)
    //{
    //  REAL_TYPE klv = ( C.Unit.Temp == Unit_KELVIN ) ? 0.0 : KELVIN;
    //  fb_tmp_nd2d_(d_ws, d_t, &d_length, &C.BaseTemp, &C.DiffTemp, &klv, &scale, &flop);
    //}
    //else
    //{
    //  fb_xcopy_(d_ws, d_t, &d_length, &scale, &flop);
    //}
    if (C.Unit.File == DIMENSIONAL) 
    {
      REAL_TYPE klv = ( C.Unit.Temp == Unit_KELVIN ) ? 0.0 : KELVIN;
      U.tmp_array_ND2D(d_ws, size, guide, d_t, C.BaseTemp, C.DiffTemp, klv, scale, flop);
    }
    else 
    {
      U.xcopy(d_ws, size, guide, d_t, scale, kind_scalar, flop);
    }

    fname = "tmp_" + tmp;
    dtmp = (C.FIO.IO_Mode == Control::io_time_slice) ? DFI.Generate_DirName(C.f_DivDebug, m_step) : C.FIO.IO_DirPath;
    fname = directory_prefix(dtmp, fname, C.FIO.IO_Mode, C.Parallelism);
    FP3DW.setFileName(fname.c_str());
    if(!FP3DW.OpenFile()){
      std::cout << "error OpenFile" << std::endl;
      Exit(0);
    }

    FP3DW.WriteNgrid(ngrid);
    FP3DW.WriteFuncBlockData(id,jd,kd,nvar_scalar,ngrid);
    //for(igrid=0;igrid<ngrid;igrid++){ //--->BCMでループが必要になる？
    igrid=0;
    FP3DW.setGridData(id[igrid],jd[igrid],kd[igrid],ngrid);
    FP3DW.setFuncDataNum(nvar_scalar[igrid]);
    if( FP3DW.GetRealType() == OUTPUT_FLOAT ){
      setScalarGridData(d,d_ws,id[igrid],jd[igrid],kd[igrid]);
      FP3DW.setFuncData(d);
    }else if( FP3DW.GetRealType() == OUTPUT_DOUBLE ){
      setScalarGridData(dd,d_ws,id[igrid],jd[igrid],kd[igrid]);
      FP3DW.setFuncData(dd);
    }
    if(!FP3DW.WriteFuncData()) std::cout << "error WriteFuncData" << std::endl;
    //}//igrid loop
    FP3DW.CloseFile();
    dfi_name = "tmp_" + C.P3Op.basename;
    Hostonly_ if ( !DFI.Write_DFI_File(dfi_name, m_step, (double)m_time, dfi_mng[var_Plot3D], pout) ) Exit(0);
  }

  // Total Pressure
  if (C.Mode.TP == ON ){
    fb_totalp_ (d_p0, size, &guide, d_v, d_p, v00, &flop);
    
    // convert non-dimensional to dimensional, iff file is dimensional
    if (C.Unit.File == DIMENSIONAL)
    {
      F.cnv_TP_ND2D(d_ws, d_p0, size, guide, C.RefDensity, C.RefVelocity, flop);
    }
    else
    {
      REAL_TYPE* tp;
      tp = d_ws; d_ws = d_p0; d_p0 = tp;
    }

    fname = "tp_" + tmp;
    dtmp = (C.FIO.IO_Mode == Control::io_time_slice) ? DFI.Generate_DirName(C.f_DivDebug, m_step) : C.FIO.IO_DirPath;
    fname = directory_prefix(dtmp, fname, C.FIO.IO_Mode, C.Parallelism);
    FP3DW.setFileName(fname.c_str());
    if(!FP3DW.OpenFile()){
      std::cout << "error OpenFile" << std::endl;
      Exit(0);
    }

    FP3DW.WriteNgrid(ngrid);
    FP3DW.WriteFuncBlockData(id,jd,kd,nvar_scalar,ngrid);
    //for(igrid=0;igrid<ngrid;igrid++){ //--->BCMでループが必要になる？
    igrid=0;
    FP3DW.setGridData(id[igrid],jd[igrid],kd[igrid],ngrid);
    FP3DW.setFuncDataNum(nvar_scalar[igrid]);
    if( FP3DW.GetRealType() == OUTPUT_FLOAT ){
      setScalarGridData(d,d_ws,id[igrid],jd[igrid],kd[igrid]);
      FP3DW.setFuncData(d);
    }else if( FP3DW.GetRealType() == OUTPUT_DOUBLE ){
      setScalarGridData(dd,d_ws,id[igrid],jd[igrid],kd[igrid]);
      FP3DW.setFuncData(dd);
    }
    if(!FP3DW.WriteFuncData()) std::cout << "error WriteFuncData" << std::endl;
    //}//igrid loop
    FP3DW.CloseFile();
    dfi_name = "tp_" + C.P3Op.basename;
    Hostonly_ if ( !DFI.Write_DFI_File(dfi_name, m_step, (double)m_time, dfi_mng[var_Plot3D], pout) ) Exit(0);
  }

  // Vorticity
  if (C.Mode.VRT == ON ){
    rot_v_(d_wv, size, &guide, &deltaX, d_v, d_bcv, v00, &flop);
    
    REAL_TYPE  vz[3];
    vz[0] = vz[1] = vz[2] = 0.0;
    unit_velocity = (C.Unit.File == DIMENSIONAL) ? C.RefVelocity/C.RefLength : 1.0;
    fb_shift_refv_out_(d_wo, d_wv, size, &guide, vz, &scale, &unit_velocity, &flop);

    fname = "vrt_" + tmp;
    dtmp = (C.FIO.IO_Mode == Control::io_time_slice) ? DFI.Generate_DirName(C.f_DivDebug, m_step) : C.FIO.IO_DirPath;
    fname = directory_prefix(dtmp, fname, C.FIO.IO_Mode, C.Parallelism);
    FP3DW.setFileName(fname.c_str());
    if(!FP3DW.OpenFile()){
      std::cout << "error OpenFile" << std::endl;
      Exit(0);
    }

    FP3DW.WriteNgrid(ngrid);
    FP3DW.WriteFuncBlockData(id,jd,kd,nvar_vector,ngrid);
    //for(igrid=0;igrid<ngrid;igrid++){ //--->BCMでループが必要になる？
    igrid=0;
    FP3DW.setGridData(id[igrid],jd[igrid],kd[igrid],ngrid);
    FP3DW.setFuncDataNum(nvar_vector[igrid]);
    if( FP3DW.GetRealType() == OUTPUT_FLOAT ){
      setVectorComponentGridData(d,d_wo,id[igrid],jd[igrid],kd[igrid],0);//0:x
      FP3DW.setFuncData(d);
    }else if( FP3DW.GetRealType() == OUTPUT_DOUBLE ){
      setVectorComponentGridData(dd,d_wo,id[igrid],jd[igrid],kd[igrid],0);//0:x
      FP3DW.setFuncData(dd);
    }
    if(!FP3DW.WriteFuncData()) std::cout << "error WriteFuncData" << std::endl;
    if( FP3DW.GetRealType() == OUTPUT_FLOAT ){
      setVectorComponentGridData(d,d_wo,id[igrid],jd[igrid],kd[igrid],1);//1:y
      FP3DW.setFuncData(d);
    }else if( FP3DW.GetRealType() == OUTPUT_DOUBLE ){
      setVectorComponentGridData(dd,d_wo,id[igrid],jd[igrid],kd[igrid],1);//1:y
      FP3DW.setFuncData(dd);
    }
    if(!FP3DW.WriteFuncData()) std::cout << "error WriteFuncData" << std::endl;
    if( FP3DW.GetRealType() == OUTPUT_FLOAT ){
      setVectorComponentGridData(d,d_wo,id[igrid],jd[igrid],kd[igrid],2);//2:z
      FP3DW.setFuncData(d);
    }else if( FP3DW.GetRealType() == OUTPUT_DOUBLE ){
      setVectorComponentGridData(dd,d_wo,id[igrid],jd[igrid],kd[igrid],2);//2:z
      FP3DW.setFuncData(dd);
    }

    //}//igrid loop
    FP3DW.CloseFile();
    dfi_name = "vrt_" + C.P3Op.basename;
    Hostonly_ if ( !DFI.Write_DFI_File(dfi_name, m_step, (double)m_time, dfi_mng[var_Plot3D], pout) ) Exit(0);
  }

  // 2nd Invariant of Velocity Gradient Tensor
  if (C.Mode.I2VGT == ON ) {

    i2vgt_ (d_p0, size, &guide, &deltaX, d_v, d_bcv, v00, &flop);
    
    // 無次元で出力
    //d_length = (size[0]+2*guide) * (size[1]+2*guide) * (size[2]+2*guide);
    //fb_xcopy_(d_ws, d_p0, &d_length, &scale, &flop);
    U.xcopy(d_ws, size, guide, d_p0, scale, kind_scalar, flop);

    fname = "iv2gt_" + tmp;
    dtmp = (C.FIO.IO_Mode == Control::io_time_slice) ? DFI.Generate_DirName(C.f_DivDebug, m_step) : C.FIO.IO_DirPath;
    fname = directory_prefix(dtmp, fname, C.FIO.IO_Mode, C.Parallelism);
    FP3DW.setFileName(fname.c_str());
    if(!FP3DW.OpenFile()){
      std::cout << "error OpenFile" << std::endl;
      Exit(0);
    }

    FP3DW.WriteNgrid(ngrid);
    FP3DW.WriteFuncBlockData(id,jd,kd,nvar_scalar,ngrid);
    //for(igrid=0;igrid<ngrid;igrid++){ //--->BCMでループが必要になる？
    igrid=0;
    FP3DW.setGridData(id[igrid],jd[igrid],kd[igrid],ngrid);
    FP3DW.setFuncDataNum(nvar_scalar[igrid]);
    if( FP3DW.GetRealType() == OUTPUT_FLOAT ){
      setScalarGridData(d,d_ws,id[igrid],jd[igrid],kd[igrid]);
      FP3DW.setFuncData(d);
    }else if( FP3DW.GetRealType() == OUTPUT_DOUBLE ){
      setScalarGridData(dd,d_ws,id[igrid],jd[igrid],kd[igrid]);
      FP3DW.setFuncData(dd);
    }
    if(!FP3DW.WriteFuncData()) std::cout << "error WriteFuncData" << std::endl;
    //}//igrid loop
    FP3DW.CloseFile();
    dfi_name = "iv2gt_" + C.P3Op.basename;
    Hostonly_ if ( !DFI.Write_DFI_File(dfi_name, m_step, (double)m_time, dfi_mng[var_Plot3D], pout) ) Exit(0);

  }


  // Helicity
  if (C.Mode.Helicity == ON ){
    helicity_(d_p0, size, &guide, &deltaX, d_v, d_bcv, v00, &flop);

    // 無次元で出力
    //d_length = (size[0]+2*guide) * (size[1]+2*guide) * (size[2]+2*guide);
    //fb_xcopy_(d_ws, d_p0, &d_length, &scale, &flop);
    U.xcopy(d_ws, size, guide, d_p0, scale, kind_scalar, flop);

    fname = "hlt_" + tmp;
    dtmp = (C.FIO.IO_Mode == Control::io_time_slice) ? DFI.Generate_DirName(C.f_DivDebug, m_step) : C.FIO.IO_DirPath;
    fname = directory_prefix(dtmp, fname, C.FIO.IO_Mode, C.Parallelism);
    FP3DW.setFileName(fname.c_str());
    if(!FP3DW.OpenFile()){
      std::cout << "error OpenFile" << std::endl;
      Exit(0);
    }

    FP3DW.WriteNgrid(ngrid);
    FP3DW.WriteFuncBlockData(id,jd,kd,nvar_scalar,ngrid);
    //for(igrid=0;igrid<ngrid;igrid++){ //--->BCMでループが必要になる？
    igrid=0;
    FP3DW.setGridData(id[igrid],jd[igrid],kd[igrid],ngrid);
    FP3DW.setFuncDataNum(nvar_scalar[igrid]);
    if( FP3DW.GetRealType() == OUTPUT_FLOAT ){
      setScalarGridData(d,d_ws,id[igrid],jd[igrid],kd[igrid]);
      FP3DW.setFuncData(d);
    }else if( FP3DW.GetRealType() == OUTPUT_DOUBLE ){
      setScalarGridData(dd,d_ws,id[igrid],jd[igrid],kd[igrid]);
      FP3DW.setFuncData(dd);
    }
    if(!FP3DW.WriteFuncData()) std::cout << "error WriteFuncData" << std::endl;
    //}//igrid loop
    FP3DW.CloseFile();
    dfi_name = "hlt_" + C.P3Op.basename;
    Hostonly_ if ( !DFI.Write_DFI_File(dfi_name, m_step, (double)m_time, dfi_mng[var_Plot3D], pout) ) Exit(0);
  }

  //deallocate
  delete [] d;
  delete [] id;
  delete [] jd;
  delete [] kd;
  
  //dfiファイルの出力
  Hostonly_ if ( !DFI.Write_DFI_File(C.P3Op.basename, m_step, (double)m_time, dfi_mng[var_Plot3D], pout) ) Exit(0);

}


// 計算結果ファイルの項目（*.nam）出力
void FFV::OutputPlot3D_function_name()
{
  //HostRankのみ出力
  if(myRank != 0 ) return;
  
  //function_nameファイルはかならずformatted形式
  int keep_format=FP3DW.GetFormat();
  FP3DW.setFormat(FORMATTED);
  
  //set filename

  // 出力ファイル名
  std::string tmp;
  std::string dtmp;

  // 並列出力モード
  bool pout = ( C.FIO.IO_Output == IO_GATHER ) ? false : true;
  
  //tmp = DFI.Generate_FileName_Free(C.P3Op.basename, "nam", 0, myRank, pout);
  tmp = DFI.Generate_FileName_Free(C.P3Op.basename, "nam", 0, myRank, false);
  dtmp = (C.FIO.IO_Mode == Control::io_time_slice) ? DFI.Generate_DirName(C.f_DivDebug, 0) : C.FIO.IO_DirPath;
  tmp = directory_prefix(dtmp, tmp, C.FIO.IO_Mode, C.Parallelism);

  //open file
  FP3DW.setFileName(tmp.c_str());
  if(!FP3DW.OpenFile()){
    Hostonly_ printf("Error : error OpenFile\n");
    Exit(0);
  }
  
  // Pressure
  FP3DW.WriteFunctionName("Pressure");
  
  // Velocity
  FP3DW.WriteFunctionName("U-Velocity ; Velocity");
  FP3DW.WriteFunctionName("V-Velocity");
  FP3DW.WriteFunctionName("W-Velocity");
  
  // Tempearture
  if( C.isHeatProblem() ) FP3DW.WriteFunctionName("Tempearture");
  
  // Total Pressure
  if (C.Mode.TP == ON ) FP3DW.WriteFunctionName("Total_Pressure");
  
  // Vorticity
  if (C.Mode.VRT == ON ){
    FP3DW.WriteFunctionName("U-Vorticity ; Vorticity");
    FP3DW.WriteFunctionName("V-Vorticity");
    FP3DW.WriteFunctionName("W-Vorticity");
  }

  // 2nd Invariant of Velocity Gradient Tensor
  if (C.Mode.I2VGT == ON ) FP3DW.WriteFunctionName("2nd Invariant of Velocity Gradient Tensor");
  
  // Helicity
  if (C.Mode.Helicity == ON ) FP3DW.WriteFunctionName("Helicity");
  
  //close file
  FP3DW.CloseFile();
  
  //reset option
  FP3DW.setFormat(keep_format);
  
}


// 項目別計算結果ファイルの項目（*.nam）出力
void FFV::OutputPlot3D_function_name_divide()
{
  //HostRankのみ出力
  if(myRank != 0 ) return;
  
  //function_nameファイルはかならずformatted形式
  int keep_format=FP3DW.GetFormat();
  FP3DW.setFormat(FORMATTED);
  
  // 並列出力モード
  bool pout = ( C.FIO.IO_Output == IO_GATHER ) ? false : true;
  
  // 出力ファイル名
  std::string tmp,fname;
  std::string dtmp;
  //tmp = DFI.Generate_FileName_Free(C.P3Op.basename, "nam", 0, myRank, pout);
  tmp = DFI.Generate_FileName_Free(C.P3Op.basename, "nam", 0, myRank, false);
  
  // Pressure
  fname = "prs_" + tmp;
  dtmp = (C.FIO.IO_Mode == Control::io_time_slice) ? DFI.Generate_DirName(C.f_DivDebug, 0) : C.FIO.IO_DirPath;
  fname = directory_prefix(dtmp, fname, C.FIO.IO_Mode, C.Parallelism);
  FP3DW.setFileName(fname.c_str());
  if(!FP3DW.OpenFile()){
    Hostonly_ printf("Error : error OpenFile\n");
    Exit(0);
  }
  FP3DW.WriteFunctionName("Pressure");
  FP3DW.CloseFile();

  // Velocity
  fname = "vel_" + tmp;
  dtmp = (C.FIO.IO_Mode == Control::io_time_slice) ? DFI.Generate_DirName(C.f_DivDebug, 0) : C.FIO.IO_DirPath;
  fname = directory_prefix(dtmp, fname, C.FIO.IO_Mode, C.Parallelism);
  FP3DW.setFileName(fname.c_str());
  if(!FP3DW.OpenFile()){
    Hostonly_ printf("Error : error OpenFile\n");
    Exit(0);
  }
  FP3DW.WriteFunctionName("U-Velocity ; Velocity");
  FP3DW.WriteFunctionName("V-Velocity");
  FP3DW.WriteFunctionName("W-Velocity");
  FP3DW.CloseFile();

  // Tempearture
  if( C.isHeatProblem() ){
    fname = "tmp_" + tmp;
    dtmp = (C.FIO.IO_Mode == Control::io_time_slice) ? DFI.Generate_DirName(C.f_DivDebug, 0) : C.FIO.IO_DirPath;
    fname = directory_prefix(dtmp, fname, C.FIO.IO_Mode, C.Parallelism);
    FP3DW.setFileName(fname.c_str());
    if(!FP3DW.OpenFile()){
      Hostonly_ printf("Error : error OpenFile\n");
      Exit(0);
    }
    FP3DW.WriteFunctionName("Tempearture");
    FP3DW.CloseFile();
  }
  
  // Total Pressure
  if (C.Mode.TP == ON ){
    fname = "tp_" + tmp;
    dtmp = (C.FIO.IO_Mode == Control::io_time_slice) ? DFI.Generate_DirName(C.f_DivDebug, 0) : C.FIO.IO_DirPath;
    fname = directory_prefix(dtmp, fname, C.FIO.IO_Mode, C.Parallelism);
    FP3DW.setFileName(fname.c_str());
    if(!FP3DW.OpenFile()){
      Hostonly_ printf("Error : error OpenFile\n");
      Exit(0);
    }
    FP3DW.WriteFunctionName("Total_Pressure");
    FP3DW.CloseFile();
  }

  // Vorticity
  if (C.Mode.VRT == ON ){
    fname = "vrt_" + tmp;
    dtmp = (C.FIO.IO_Mode == Control::io_time_slice) ? DFI.Generate_DirName(C.f_DivDebug, 0) : C.FIO.IO_DirPath;
    fname = directory_prefix(dtmp, fname, C.FIO.IO_Mode, C.Parallelism);
    FP3DW.setFileName(fname.c_str());
    if(!FP3DW.OpenFile()){
      Hostonly_ printf("Error : error OpenFile\n");
      Exit(0);
    }
    FP3DW.WriteFunctionName("U-Vorticity ; Vorticity");
    FP3DW.WriteFunctionName("V-Vorticity");
    FP3DW.WriteFunctionName("W-Vorticity");
    FP3DW.CloseFile();
  }

  // 2nd Invariant of Velocity Gradient Tensor
  if (C.Mode.I2VGT == ON ){
    fname = "i2vgt_" + tmp;
    dtmp = (C.FIO.IO_Mode == Control::io_time_slice) ? DFI.Generate_DirName(C.f_DivDebug, 0) : C.FIO.IO_DirPath;
    fname = directory_prefix(dtmp, fname, C.FIO.IO_Mode, C.Parallelism);
    FP3DW.setFileName(fname.c_str());
    if(!FP3DW.OpenFile()){
      Hostonly_ printf("Error : error OpenFile\n");
      Exit(0);
    }
    FP3DW.WriteFunctionName("2nd Invariant of Velocity Gradient Tensor");
    FP3DW.CloseFile();
  }


  // Helicity
  if (C.Mode.Helicity == ON ){
    fname = "hlt_" + tmp;
    dtmp = (C.FIO.IO_Mode == Control::io_time_slice) ? DFI.Generate_DirName(C.f_DivDebug, 0) : C.FIO.IO_DirPath;
    fname = directory_prefix(dtmp, fname, C.FIO.IO_Mode, C.Parallelism);
    FP3DW.setFileName(fname.c_str());
    if(!FP3DW.OpenFile()){
      Hostonly_ printf("Error : error OpenFile\n");
      Exit(0);
    }
    FP3DW.WriteFunctionName("Helicity");
    FP3DW.CloseFile();
  }

  //reset option
  FP3DW.setFormat(keep_format);
  
}


// 境界面定義ファイル（*.fvbnd）出力
void FFV::OutputPlot3D_fvbnd()
{
  //領域を超えてfvbndを出力するのは結構難しい＆あまりメリットはない--->当面保留
  return;

  //HostRankのみ出力
  if(myRank != 0 ) return;

//////
//////  //fvbndファイルはかならずformatted形式
//////  int keep_format=FP3DW.GetFormat();
//////  FP3DW.setFormat(FORMATTED);
//////
//////  // 出力ファイル名
//////  std::string tmp;
//////  std::string dtmp;
//////
//////  // 並列出力モード
//////  bool pout = ( C.FIO.IO_Output == IO_GATHER ) ? false : true;
//////  
//////  //tmp = DFI.Generate_FileName_Free(C.P3Op.basename, "xyz", 0, myRank, pout);
//////  tmp = DFI.Generate_FileName_Free(C.P3Op.basename, "xyz", 0, myRank, false);
//////  tmp = tmp + ".fvbnd";
//////  dtmp = (C.FIO.IO_Mode == Control::io_time_slice) ? DFI.Generate_DirName(C.f_DivDebug, 0) : C.FIO.IO_DirPath;
//////  tmp = directory_prefix(dtmp, tmp, C.FIO.IO_Mode, C.Parallelism);
//////
//////  //open file
//////  FP3DW.setFileName(tmp.c_str());
//////  if(!FP3DW.OpenFile()){
//////    Hostonly_ printf("Error : error OpenFile\n");
//////    Exit(0);
//////  }
//////
//////  //境界名の取得--->将来的に独立ルーチンにして名前だけ保持し続ける？*.fvbndが複数の場合に対処
//////  vector<string> bcname;
//////  bcname.clear();
//////  B.GetBoundaryNameforPLOT3D(bcname,cmp);
//////
//////#if 0
//////  vector<string>::const_iterator it;
//////  for (it = bcname.begin(); it != bcname.end(); it++) {
//////    cout << "name = " << (*it).c_str() << endl;
//////  }
//////#endif
//////
//////  //write boundary
//////  BC.WriteBoundaryPLOT3D(&FP3DW,bcname);
//////
//////  //close file
//////  FP3DW.CloseFile();
//////  
//////  //reset option
//////  FP3DW.setFormat(keep_format);
//////
//////  return;
}


// Iblankのセット（ガイドセルの値は計算対象にいれていない）
void FFV::setIblank(int* iblank, int id, int jd, int kd)
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

  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;

  for (int km=1; km<=kx; km++) {
    for (int jm=1; jm<=jx; jm++) {
      for (int im=1; im<=ix; im++) {
        mip = _F_IDX_S3D(im, jm, km, ix, jx, kx, gd);
        s = d_bcd[mip];
        //if(d_mid[mip]==1){//流体であれば
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
          //size_t ip1= k   *id*jd+ j   *id+i;
          //size_t ip2= k   *id*jd+ j   *id+i+1;
          //size_t ip3= k   *id*jd+(j+1)*id+i+1;
          //size_t ip4= k   *id*jd+(j+1)*id+i;
          //size_t ip5=(k+1)*id*jd+ j   *id+i;
          //size_t ip6=(k+1)*id*jd+ j   *id+i+1;
          //size_t ip7=(k+1)*id*jd+(j+1)*id+i+1;
          //size_t ip8=(k+1)*id*jd+(j+1)*id+i;
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

#if 0

  cout << endl;
  cout << "iblank debug" << endl;
  for(k=0;k<kd;k++){
    for(j=0;j<jd;j++){
      //cout << "j = " << j << " k = " << k << endl;
      for(i=0;i<id;i++){
        ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
        //cout << "iblank[" << ip << "] = " << iblank[ip];
        cout << " " << iblank[ip];
      }
        cout << endl;
    }
      cout << endl;
  }

#endif

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
        //if(d_mid[mip]==id_of_solid){//固体であれば
        //if(d_mid[mip]!=target_id){//流体でなければ
        //if(d_mid[mip]!=1){//流体でなければ
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
          //size_t ip1= k   *id*jd+ j   *id+i;
          //size_t ip2= k   *id*jd+ j   *id+i+1;
          //size_t ip3= k   *id*jd+(j+1)*id+i+1;
          //size_t ip4= k   *id*jd+(j+1)*id+i;
          //size_t ip5=(k+1)*id*jd+ j   *id+i;
          //size_t ip6=(k+1)*id*jd+ j   *id+i+1;
          //size_t ip7=(k+1)*id*jd+(j+1)*id+i+1;
          //size_t ip8=(k+1)*id*jd+(j+1)*id+i;
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
              //size_t ipwk = _F_IDX_S3D(i+l+1, j+m+1, k+n+1, id, jd, kd, 0);
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
          //size_t ipwk = _F_IDX_S3D(i+l+1, j+m+1, k+n+1, id, jd, kd, 0);
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
          //size_t ipwk = _F_IDX_S3D(i+l+1, j+m+1, k+n+1, id, jd, kd, 0);
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
          //size_t ipwk = _F_IDX_S3D(i+l+1, j+m+1, k+n+1, id, jd, kd, 0);
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
          //size_t ipwk = _F_IDX_S3D(i+l+1, j+m+1, k+n+1, id, jd, kd, 0);
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
          //size_t ipwk = _F_IDX_S3D(i+l+1, j+m+1, k+n+1, id, jd, kd, 0);
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
          //size_t ipwk = _F_IDX_S3D(i+l+1, j+m+1, k+n+1, id, jd, kd, 0);
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
      //size_t ipwk = _F_IDX_S3D(i+l+1, j+m+1, k+n+1, id, jd, kd, 0);
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
      //size_t ipwk = _F_IDX_S3D(i+l+1, j+m+1, k+n+1, id, jd, kd, 0);
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
      //size_t ipwk = _F_IDX_S3D(i+l+1, j+m+1, k+n+1, id, jd, kd, 0);
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
      //size_t ipwk = _F_IDX_S3D(i+l+1, j+m+1, k+n+1, id, jd, kd, 0);
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
      //size_t ipwk = _F_IDX_S3D(i+l+1, j+m+1, k+n+1, id, jd, kd, 0);
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
      //size_t ipwk = _F_IDX_S3D(i+l+1, j+m+1, k+n+1, id, jd, kd, 0);
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
      //size_t ipwk = _F_IDX_S3D(i+l+1, j+m+1, k+n+1, id, jd, kd, 0);
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
      //size_t ipwk = _F_IDX_S3D(i+l+1, j+m+1, k+n+1, id, jd, kd, 0);
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
      //size_t ipwk = _F_IDX_S3D(i+l+1, j+m+1, k+n+1, id, jd, kd, 0);
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
      //size_t ipwk = _F_IDX_S3D(i+l+1, j+m+1, k+n+1, id, jd, kd, 0);
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
      //size_t ipwk = _F_IDX_S3D(i+l+1, j+m+1, k+n+1, id, jd, kd, 0);
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
    //ipxp= k   *id*jd+ j   *id+i+1;
    //ipyp= k   *id*jd+(j+1)*id+i;
    //ipzp=(k+1)*id*jd+ j   *id+i;
    //ipcr=(k+1)*id*jd+(j+1)*id+i+1;
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
    //ipxm= k   *id*jd+ j   *id+i-1;
    //ipyp= k   *id*jd+(j+1)*id+i;
    //ipzp=(k+1)*id*jd+ j   *id+i;
    //ipcr=(k+1)*id*jd+(j+1)*id+i-1;
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
    //ipxm= k   *id*jd+ j   *id+i-1;
    //ipym= k   *id*jd+(j-1)*id+i;
    //ipzp=(k+1)*id*jd+ j   *id+i;
    //ipcr=(k+1)*id*jd+(j-1)*id+i-1;
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
    //ipxp= k   *id*jd+ j   *id+i+1;
    //ipym= k   *id*jd+(j-1)*id+i;
    //ipzp=(k+1)*id*jd+ j   *id+i;
    //ipcr=(k+1)*id*jd+(j-1)*id+i+1;
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
    //ipxp= k   *id*jd+ j   *id+i+1;
    //ipyp= k   *id*jd+(j+1)*id+i;
    //ipzm=(k-1)*id*jd+ j   *id+i;
    //ipcr=(k-1)*id*jd+(j+1)*id+i+1;
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
    //ipxm= k   *id*jd+ j   *id+i-1;
    //ipyp= k   *id*jd+(j+1)*id+i;
    //ipzm=(k-1)*id*jd+ j   *id+i;
    //ipcr=(k-1)*id*jd+(j+1)*id+i-1;
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
    //ipxm= k   *id*jd+ j   *id+i-1;
    //ipym= k   *id*jd+(j-1)*id+i;
    //ipzm=(k-1)*id*jd+ j   *id+i;
    //ipcr=(k-1)*id*jd+(j-1)*id+i-1;
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
    //ipxp= k   *id*jd+ j   *id+i+1;
    //ipym= k   *id*jd+(j-1)*id+i;
    //ipzm=(k-1)*id*jd+ j   *id+i;
    //ipcr=(k-1)*id*jd+(j-1)*id+i-1;
    if( iblank[ipxp]==1 || iblank[ipym]==1 || iblank[ipzm]==1 || iblank[ipcr]==1){
      iblank[ip]=2;
    }
  }

#endif

#if 0

  cout << endl;
  cout << "iblank debug" << endl;
  for(k=0;k<kd;k++){
    for(j=0;j<jd;j++){
      //cout << "j = " << j << " k = " << k << endl;
      for(i=0;i<id;i++){
        ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
        //cout << "iblank[" << ip << "] = " << iblank[ip];
        cout << " " << iblank[ip];
      }
        cout << endl;
    }
      cout << endl;
  }

#endif

}


//*****************************************************************************
//*****************************************************************************
//*****************************************************************************

// Scalarの格子点での値をセット（ガイドセルの値は計算対象にいれていない）
void FFV::setScalarGridData(float* d, float* data, int id, int jd, int kd)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;

  size_t dsize = (size_t)(id*jd*kd);

  for (size_t i=0; i<dsize; i++) d[i]=0.0;

  for (int km=1; km<=kx; km++) {
    for (int jm=1; jm<=jx; jm++) {
      for (int im=1; im<=ix; im++) {
        size_t mip = _F_IDX_S3D(im, jm, km, ix, jx, kx, gd);
        float ddd=data[mip];
        int i=im-1;
        int j=jm-1;
        int k=km-1;
        size_t ip1 = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
        size_t ip2 = _F_IDX_S3D(i+2, j+1, k+1, id, jd, kd, 0);
        size_t ip3 = _F_IDX_S3D(i+2, j+2, k+1, id, jd, kd, 0);
        size_t ip4 = _F_IDX_S3D(i+1, j+2, k+1, id, jd, kd, 0);
        size_t ip5 = _F_IDX_S3D(i+1, j+1, k+2, id, jd, kd, 0);
        size_t ip6 = _F_IDX_S3D(i+2, j+1, k+2, id, jd, kd, 0);
        size_t ip7 = _F_IDX_S3D(i+2, j+2, k+2, id, jd, kd, 0);
        size_t ip8 = _F_IDX_S3D(i+1, j+2, k+2, id, jd, kd, 0);
        //size_t ip1= k   *id*jd+ j   *id+i;
        //size_t ip2= k   *id*jd+ j   *id+i+1;
        //size_t ip3= k   *id*jd+(j+1)*id+i+1;
        //size_t ip4= k   *id*jd+(j+1)*id+i;
        //size_t ip5=(k+1)*id*jd+ j   *id+i;
        //size_t ip6=(k+1)*id*jd+ j   *id+i+1;
        //size_t ip7=(k+1)*id*jd+(j+1)*id+i+1;
        //size_t ip8=(k+1)*id*jd+(j+1)*id+i;
        d[ip1]=d[ip1]+ddd;
        d[ip2]=d[ip2]+ddd;
        d[ip3]=d[ip3]+ddd;
        d[ip4]=d[ip4]+ddd;
        d[ip5]=d[ip5]+ddd;
        d[ip6]=d[ip6]+ddd;
        d[ip7]=d[ip7]+ddd;
        d[ip8]=d[ip8]+ddd;
      }
    }
  }

  //内部の格子点のデータを8で割る
  VolumeDataDivideBy8(d, id, jd, kd);

  //面上の格子点のデータを4で割る
  FaceDataDivideBy4(d, id, jd, kd);

  //辺上の格子点のデータを2で割る
  LineDataDivideBy2(d, id, jd, kd);

  //境界条件処理（wkdにガイドセルを含む格子点上のデータ）

}

// Vectorの格子点での値をセット（ガイドセルの値は計算対象にいれていない）
void FFV::setVectorGridData(float* d, float* data, int id, int jd, int kd)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;

  size_t dsize = (size_t)(id*jd*kd);
  size_t dsize3 = (size_t)(id*jd*kd*3);

  for (size_t i=0; i<dsize3; i++) d[i]=0.0;

  for (size_t ivar=0;ivar<3;ivar++){

    for (int km=1; km<=kx; km++) {
      for (int jm=1; jm<=jx; jm++) {
        for (int im=1; im<=ix; im++) {
          size_t mip = _F_IDX_V3DEX(ivar, im, jm, km, ix, jx, kx, gd);//ivar=0:x方向、ivar=1:y方向、ivar=2:z方向
          float ddd=data[mip];
          int i=im-1;
          int j=jm-1;
          int k=km-1;
          size_t ip1 = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
          size_t ip2 = _F_IDX_S3D(i+2, j+1, k+1, id, jd, kd, 0);
          size_t ip3 = _F_IDX_S3D(i+2, j+2, k+1, id, jd, kd, 0);
          size_t ip4 = _F_IDX_S3D(i+1, j+2, k+1, id, jd, kd, 0);
          size_t ip5 = _F_IDX_S3D(i+1, j+1, k+2, id, jd, kd, 0);
          size_t ip6 = _F_IDX_S3D(i+2, j+1, k+2, id, jd, kd, 0);
          size_t ip7 = _F_IDX_S3D(i+2, j+2, k+2, id, jd, kd, 0);
          size_t ip8 = _F_IDX_S3D(i+1, j+2, k+2, id, jd, kd, 0);
          //size_t ip1= k   *id*jd+ j   *id+i;
          //size_t ip2= k   *id*jd+ j   *id+i+1;
          //size_t ip3= k   *id*jd+(j+1)*id+i+1;
          //size_t ip4= k   *id*jd+(j+1)*id+i;
          //size_t ip5=(k+1)*id*jd+ j   *id+i;
          //size_t ip6=(k+1)*id*jd+ j   *id+i+1;
          //size_t ip7=(k+1)*id*jd+(j+1)*id+i+1;
          //size_t ip8=(k+1)*id*jd+(j+1)*id+i;
          d[ip1+dsize*ivar]=d[ip1+dsize*ivar]+ddd;
          d[ip2+dsize*ivar]=d[ip2+dsize*ivar]+ddd;
          d[ip3+dsize*ivar]=d[ip3+dsize*ivar]+ddd;
          d[ip4+dsize*ivar]=d[ip4+dsize*ivar]+ddd;
          d[ip5+dsize*ivar]=d[ip5+dsize*ivar]+ddd;
          d[ip6+dsize*ivar]=d[ip6+dsize*ivar]+ddd;
          d[ip7+dsize*ivar]=d[ip7+dsize*ivar]+ddd;
          d[ip8+dsize*ivar]=d[ip8+dsize*ivar]+ddd;
        }
      }
    }

    //内部の格子点のデータを8で割る
    VolumeDataDivideBy8(&d[dsize*ivar], id, jd, kd);

    //面上の格子点のデータを4で割る
    FaceDataDivideBy4(&d[dsize*ivar], id, jd, kd);

    //辺上の格子点のデータを2で割る
    LineDataDivideBy2(&d[dsize*ivar], id, jd, kd);

    //境界条件処理（wkdにガイドセルを含む格子点上のデータ）
 

  }//loop ivar

}

// Vectorの格子点での値をセット（ガイドセルの値は計算対象にいれていない）
void FFV::setVectorComponentGridData(float* d, float* data, int id, int jd, int kd, int ivar)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;

  size_t dsize = (size_t)(id*jd*kd);

  for (size_t i=0; i<dsize; i++) d[i]=0.0;
 
  for (int km=1; km<=kx; km++) {
    for (int jm=1; jm<=jx; jm++) {
      for (int im=1; im<=ix; im++) {
        size_t mip = _F_IDX_V3DEX(ivar, im, jm, km, ix, jx, kx, gd);//ivar=0:x方向、ivar=1:y方向、ivar=2:z方向
        float ddd=data[mip];
        int i=im-1;
        int j=jm-1;
        int k=km-1;
        size_t ip1 = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
        size_t ip2 = _F_IDX_S3D(i+2, j+1, k+1, id, jd, kd, 0);
        size_t ip3 = _F_IDX_S3D(i+2, j+2, k+1, id, jd, kd, 0);
        size_t ip4 = _F_IDX_S3D(i+1, j+2, k+1, id, jd, kd, 0);
        size_t ip5 = _F_IDX_S3D(i+1, j+1, k+2, id, jd, kd, 0);
        size_t ip6 = _F_IDX_S3D(i+2, j+1, k+2, id, jd, kd, 0);
        size_t ip7 = _F_IDX_S3D(i+2, j+2, k+2, id, jd, kd, 0);
        size_t ip8 = _F_IDX_S3D(i+1, j+2, k+2, id, jd, kd, 0);
        //size_t ip1= k   *id*jd+ j   *id+i;
        //size_t ip2= k   *id*jd+ j   *id+i+1;
        //size_t ip3= k   *id*jd+(j+1)*id+i+1;
        //size_t ip4= k   *id*jd+(j+1)*id+i;
        //size_t ip5=(k+1)*id*jd+ j   *id+i;
        //size_t ip6=(k+1)*id*jd+ j   *id+i+1;
        //size_t ip7=(k+1)*id*jd+(j+1)*id+i+1;
        //size_t ip8=(k+1)*id*jd+(j+1)*id+i;
        d[ip1]=d[ip1]+ddd;
        d[ip2]=d[ip2]+ddd;
        d[ip3]=d[ip3]+ddd;
        d[ip4]=d[ip4]+ddd;
        d[ip5]=d[ip5]+ddd;
        d[ip6]=d[ip6]+ddd;
        d[ip7]=d[ip7]+ddd;
        d[ip8]=d[ip8]+ddd;
      }
    }
  }

  //内部の格子点のデータを8で割る
  VolumeDataDivideBy8(d, id, jd, kd);

  //面上の格子点のデータを4で割る
  FaceDataDivideBy4(d, id, jd, kd);

  //辺上の格子点のデータを2で割る
  LineDataDivideBy2(d, id, jd, kd);

  //境界条件処理（wkdにガイドセルを含む格子点上のデータ）

}

//*****************************************************************************
//*****************************************************************************
//*****************************************************************************

// Scalarの格子点での値をセット（ガイドセルの値は計算対象にいれていない）
void FFV::setScalarGridData(double* d, double* data, int id, int jd, int kd)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;

  size_t dsize = (size_t)(id*jd*kd);

  for (size_t i=0; i<dsize; i++) d[i]=0.0;

  for (int km=1; km<=kx; km++) {
    for (int jm=1; jm<=jx; jm++) {
      for (int im=1; im<=ix; im++) {
        size_t mip = _F_IDX_S3D(im, jm, km, ix, jx, kx, gd);
        double ddd=data[mip];
        int i=im-1;
        int j=jm-1;
        int k=km-1;
        size_t ip1 = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
        size_t ip2 = _F_IDX_S3D(i+2, j+1, k+1, id, jd, kd, 0);
        size_t ip3 = _F_IDX_S3D(i+2, j+2, k+1, id, jd, kd, 0);
        size_t ip4 = _F_IDX_S3D(i+1, j+2, k+1, id, jd, kd, 0);
        size_t ip5 = _F_IDX_S3D(i+1, j+1, k+2, id, jd, kd, 0);
        size_t ip6 = _F_IDX_S3D(i+2, j+1, k+2, id, jd, kd, 0);
        size_t ip7 = _F_IDX_S3D(i+2, j+2, k+2, id, jd, kd, 0);
        size_t ip8 = _F_IDX_S3D(i+1, j+2, k+2, id, jd, kd, 0);
        //size_t ip1= k   *id*jd+ j   *id+i;
        //size_t ip2= k   *id*jd+ j   *id+i+1;
        //size_t ip3= k   *id*jd+(j+1)*id+i+1;
        //size_t ip4= k   *id*jd+(j+1)*id+i;
        //size_t ip5=(k+1)*id*jd+ j   *id+i;
        //size_t ip6=(k+1)*id*jd+ j   *id+i+1;
        //size_t ip7=(k+1)*id*jd+(j+1)*id+i+1;
        //size_t ip8=(k+1)*id*jd+(j+1)*id+i;
        d[ip1]=d[ip1]+ddd;
        d[ip2]=d[ip2]+ddd;
        d[ip3]=d[ip3]+ddd;
        d[ip4]=d[ip4]+ddd;
        d[ip5]=d[ip5]+ddd;
        d[ip6]=d[ip6]+ddd;
        d[ip7]=d[ip7]+ddd;
        d[ip8]=d[ip8]+ddd;
      }
    }
  }

  //内部の格子点のデータを8で割る
  VolumeDataDivideBy8(d, id, jd, kd);

  //面上の格子点のデータを4で割る
  FaceDataDivideBy4(d, id, jd, kd);

  //辺上の格子点のデータを2で割る
  LineDataDivideBy2(d, id, jd, kd);

  //境界条件処理（wkdにガイドセルを含む格子点上のデータ）

}

// Vectorの格子点での値をセット（ガイドセルの値は計算対象にいれていない）
void FFV::setVectorGridData(double* d, double* data, int id, int jd, int kd)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;

  size_t dsize = (size_t)(id*jd*kd);
  size_t dsize3 = (size_t)(id*jd*kd*3);

  for (size_t i=0; i<dsize3; i++) d[i]=0.0;

  for (size_t ivar=0;ivar<3;ivar++){

    for (int km=1; km<=kx; km++) {
      for (int jm=1; jm<=jx; jm++) {
        for (int im=1; im<=ix; im++) {
          size_t mip = _F_IDX_V3DEX(ivar, im, jm, km, ix, jx, kx, gd);//ivar=0:x方向、ivar=1:y方向、ivar=2:z方向
          double ddd=data[mip];
          int i=im-1;
          int j=jm-1;
          int k=km-1;
          size_t ip1 = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
          size_t ip2 = _F_IDX_S3D(i+2, j+1, k+1, id, jd, kd, 0);
          size_t ip3 = _F_IDX_S3D(i+2, j+2, k+1, id, jd, kd, 0);
          size_t ip4 = _F_IDX_S3D(i+1, j+2, k+1, id, jd, kd, 0);
          size_t ip5 = _F_IDX_S3D(i+1, j+1, k+2, id, jd, kd, 0);
          size_t ip6 = _F_IDX_S3D(i+2, j+1, k+2, id, jd, kd, 0);
          size_t ip7 = _F_IDX_S3D(i+2, j+2, k+2, id, jd, kd, 0);
          size_t ip8 = _F_IDX_S3D(i+1, j+2, k+2, id, jd, kd, 0);
          //size_t ip1= k   *id*jd+ j   *id+i;
          //size_t ip2= k   *id*jd+ j   *id+i+1;
          //size_t ip3= k   *id*jd+(j+1)*id+i+1;
          //size_t ip4= k   *id*jd+(j+1)*id+i;
          //size_t ip5=(k+1)*id*jd+ j   *id+i;
          //size_t ip6=(k+1)*id*jd+ j   *id+i+1;
          //size_t ip7=(k+1)*id*jd+(j+1)*id+i+1;
          //size_t ip8=(k+1)*id*jd+(j+1)*id+i;
          d[ip1+dsize*ivar]=d[ip1+dsize*ivar]+ddd;
          d[ip2+dsize*ivar]=d[ip2+dsize*ivar]+ddd;
          d[ip3+dsize*ivar]=d[ip3+dsize*ivar]+ddd;
          d[ip4+dsize*ivar]=d[ip4+dsize*ivar]+ddd;
          d[ip5+dsize*ivar]=d[ip5+dsize*ivar]+ddd;
          d[ip6+dsize*ivar]=d[ip6+dsize*ivar]+ddd;
          d[ip7+dsize*ivar]=d[ip7+dsize*ivar]+ddd;
          d[ip8+dsize*ivar]=d[ip8+dsize*ivar]+ddd;
        }
      }
    }

    //内部の格子点のデータを8で割る
    VolumeDataDivideBy8(&d[dsize*ivar], id, jd, kd);

    //面上の格子点のデータを4で割る
    FaceDataDivideBy4(&d[dsize*ivar], id, jd, kd);

    //辺上の格子点のデータを2で割る
    LineDataDivideBy2(&d[dsize*ivar], id, jd, kd);

    //境界条件処理（wkdにガイドセルを含む格子点上のデータ）
 
  }//loop ivar

}

// Vectorの格子点での値をセット（ガイドセルの値は計算対象にいれていない）
void FFV::setVectorComponentGridData(double* d, double* data, int id, int jd, int kd, int ivar)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;

  size_t dsize = (size_t)(id*jd*kd);

  for (size_t i=0; i<dsize; i++) d[i]=0.0;
 
  for (int km=1; km<=kx; km++) {
    for (int jm=1; jm<=jx; jm++) {
      for (int im=1; im<=ix; im++) {
        size_t mip = _F_IDX_V3DEX(ivar, im, jm, km, ix, jx, kx, gd);//ivar=0:x方向、ivar=1:y方向、ivar=2:z方向
        double ddd=data[mip];
        int i=im-1;
        int j=jm-1;
        int k=km-1;
        size_t ip1 = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
        size_t ip2 = _F_IDX_S3D(i+2, j+1, k+1, id, jd, kd, 0);
        size_t ip3 = _F_IDX_S3D(i+2, j+2, k+1, id, jd, kd, 0);
        size_t ip4 = _F_IDX_S3D(i+1, j+2, k+1, id, jd, kd, 0);
        size_t ip5 = _F_IDX_S3D(i+1, j+1, k+2, id, jd, kd, 0);
        size_t ip6 = _F_IDX_S3D(i+2, j+1, k+2, id, jd, kd, 0);
        size_t ip7 = _F_IDX_S3D(i+2, j+2, k+2, id, jd, kd, 0);
        size_t ip8 = _F_IDX_S3D(i+1, j+2, k+2, id, jd, kd, 0);
        //size_t ip1= k   *id*jd+ j   *id+i;
        //size_t ip2= k   *id*jd+ j   *id+i+1;
        //size_t ip3= k   *id*jd+(j+1)*id+i+1;
        //size_t ip4= k   *id*jd+(j+1)*id+i;
        //size_t ip5=(k+1)*id*jd+ j   *id+i;
        //size_t ip6=(k+1)*id*jd+ j   *id+i+1;
        //size_t ip7=(k+1)*id*jd+(j+1)*id+i+1;
        //size_t ip8=(k+1)*id*jd+(j+1)*id+i;
        d[ip1]=d[ip1]+ddd;
        d[ip2]=d[ip2]+ddd;
        d[ip3]=d[ip3]+ddd;
        d[ip4]=d[ip4]+ddd;
        d[ip5]=d[ip5]+ddd;
        d[ip6]=d[ip6]+ddd;
        d[ip7]=d[ip7]+ddd;
        d[ip8]=d[ip8]+ddd;
      }
    }
  }

  //内部の格子点のデータを8で割る
  VolumeDataDivideBy8(d, id, jd, kd);

  //面上の格子点のデータを4で割る
  FaceDataDivideBy4(d, id, jd, kd);

  //辺上の格子点のデータを2で割る
  LineDataDivideBy2(d, id, jd, kd);

  //境界条件処理（wkdにガイドセルを含む格子点上のデータ）


}

//*****************************************************************************
//*****************************************************************************
//*****************************************************************************

// Scalarの格子点での値をセット（ガイドセルの値は計算対象にいれていない）
void FFV::setScalarGridData(float* d, double* data, int id, int jd, int kd)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;

  size_t dsize = (size_t)(id*jd*kd);

  for (size_t i=0; i<dsize; i++) d[i]=0.0;

  for (int km=1; km<=kx; km++) {
    for (int jm=1; jm<=jx; jm++) {
      for (int im=1; im<=ix; im++) {
        size_t mip = _F_IDX_S3D(im, jm, km, ix, jx, kx, gd);
        float ddd=(float)data[mip];
        int i=im-1;
        int j=jm-1;
        int k=km-1;
        size_t ip1 = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
        size_t ip2 = _F_IDX_S3D(i+2, j+1, k+1, id, jd, kd, 0);
        size_t ip3 = _F_IDX_S3D(i+2, j+2, k+1, id, jd, kd, 0);
        size_t ip4 = _F_IDX_S3D(i+1, j+2, k+1, id, jd, kd, 0);
        size_t ip5 = _F_IDX_S3D(i+1, j+1, k+2, id, jd, kd, 0);
        size_t ip6 = _F_IDX_S3D(i+2, j+1, k+2, id, jd, kd, 0);
        size_t ip7 = _F_IDX_S3D(i+2, j+2, k+2, id, jd, kd, 0);
        size_t ip8 = _F_IDX_S3D(i+1, j+2, k+2, id, jd, kd, 0);
        //size_t ip1= k   *id*jd+ j   *id+i;
        //size_t ip2= k   *id*jd+ j   *id+i+1;
        //size_t ip3= k   *id*jd+(j+1)*id+i+1;
        //size_t ip4= k   *id*jd+(j+1)*id+i;
        //size_t ip5=(k+1)*id*jd+ j   *id+i;
        //size_t ip6=(k+1)*id*jd+ j   *id+i+1;
        //size_t ip7=(k+1)*id*jd+(j+1)*id+i+1;
        //size_t ip8=(k+1)*id*jd+(j+1)*id+i;
        d[ip1]=d[ip1]+ddd;
        d[ip2]=d[ip2]+ddd;
        d[ip3]=d[ip3]+ddd;
        d[ip4]=d[ip4]+ddd;
        d[ip5]=d[ip5]+ddd;
        d[ip6]=d[ip6]+ddd;
        d[ip7]=d[ip7]+ddd;
        d[ip8]=d[ip8]+ddd;
      }
    }
  }

  //内部の格子点のデータを8で割る
  VolumeDataDivideBy8(d, id, jd, kd);

  //面上の格子点のデータを4で割る
  FaceDataDivideBy4(d, id, jd, kd);

  //辺上の格子点のデータを2で割る
  LineDataDivideBy2(d, id, jd, kd);

  //境界条件処理（wkdにガイドセルを含む格子点上のデータ）

}

// Vectorの格子点での値をセット（ガイドセルの値は計算対象にいれていない）
void FFV::setVectorGridData(float* d, double* data, int id, int jd, int kd)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;

  size_t dsize = (size_t)(id*jd*kd);
  size_t dsize3 = (size_t)(id*jd*kd*3);

  for (size_t i=0; i<dsize3; i++) d[i]=0.0;

  for (size_t ivar=0;ivar<3;ivar++){

    for (int km=1; km<=kx; km++) {
      for (int jm=1; jm<=jx; jm++) {
        for (int im=1; im<=ix; im++) {
          size_t mip = _F_IDX_V3DEX(ivar, im, jm, km, ix, jx, kx, gd);//ivar=0:x方向、ivar=1:y方向、ivar=2:z方向
          float ddd=(float)data[mip];
          int i=im-1;
          int j=jm-1;
          int k=km-1;
          size_t ip1 = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
          size_t ip2 = _F_IDX_S3D(i+2, j+1, k+1, id, jd, kd, 0);
          size_t ip3 = _F_IDX_S3D(i+2, j+2, k+1, id, jd, kd, 0);
          size_t ip4 = _F_IDX_S3D(i+1, j+2, k+1, id, jd, kd, 0);
          size_t ip5 = _F_IDX_S3D(i+1, j+1, k+2, id, jd, kd, 0);
          size_t ip6 = _F_IDX_S3D(i+2, j+1, k+2, id, jd, kd, 0);
          size_t ip7 = _F_IDX_S3D(i+2, j+2, k+2, id, jd, kd, 0);
          size_t ip8 = _F_IDX_S3D(i+1, j+2, k+2, id, jd, kd, 0);
          //size_t ip1= k   *id*jd+ j   *id+i;
          //size_t ip2= k   *id*jd+ j   *id+i+1;
          //size_t ip3= k   *id*jd+(j+1)*id+i+1;
          //size_t ip4= k   *id*jd+(j+1)*id+i;
          //size_t ip5=(k+1)*id*jd+ j   *id+i;
          //size_t ip6=(k+1)*id*jd+ j   *id+i+1;
          //size_t ip7=(k+1)*id*jd+(j+1)*id+i+1;
          //size_t ip8=(k+1)*id*jd+(j+1)*id+i;
          d[ip1+dsize*ivar]=d[ip1+dsize*ivar]+ddd;
          d[ip2+dsize*ivar]=d[ip2+dsize*ivar]+ddd;
          d[ip3+dsize*ivar]=d[ip3+dsize*ivar]+ddd;
          d[ip4+dsize*ivar]=d[ip4+dsize*ivar]+ddd;
          d[ip5+dsize*ivar]=d[ip5+dsize*ivar]+ddd;
          d[ip6+dsize*ivar]=d[ip6+dsize*ivar]+ddd;
          d[ip7+dsize*ivar]=d[ip7+dsize*ivar]+ddd;
          d[ip8+dsize*ivar]=d[ip8+dsize*ivar]+ddd;
        }
      }
    }

    //内部の格子点のデータを8で割る
    VolumeDataDivideBy8(&d[dsize*ivar], id, jd, kd);

    //面上の格子点のデータを4で割る
    FaceDataDivideBy4(&d[dsize*ivar], id, jd, kd);

    //辺上の格子点のデータを2で割る
    LineDataDivideBy2(&d[dsize*ivar], id, jd, kd);

    //境界条件処理（wkdにガイドセルを含む格子点上のデータ）

  }//loop ivar

}

// Vectorの格子点での値をセット（ガイドセルの値は計算対象にいれていない）
void FFV::setVectorComponentGridData(float* d, double* data, int id, int jd, int kd, int ivar)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;

  size_t dsize = (size_t)(id*jd*kd);

  for (size_t i=0; i<dsize; i++) d[i]=0.0;
 
  for (int km=1; km<=kx; km++) {
    for (int jm=1; jm<=jx; jm++) {
      for (int im=1; im<=ix; im++) {
        size_t mip = _F_IDX_V3DEX(ivar, im, jm, km, ix, jx, kx, gd);//ivar=0:x方向、ivar=1:y方向、ivar=2:z方向
        float ddd=(float)data[mip];
        int i=im-1;
        int j=jm-1;
        int k=km-1;
        size_t ip1 = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
        size_t ip2 = _F_IDX_S3D(i+2, j+1, k+1, id, jd, kd, 0);
        size_t ip3 = _F_IDX_S3D(i+2, j+2, k+1, id, jd, kd, 0);
        size_t ip4 = _F_IDX_S3D(i+1, j+2, k+1, id, jd, kd, 0);
        size_t ip5 = _F_IDX_S3D(i+1, j+1, k+2, id, jd, kd, 0);
        size_t ip6 = _F_IDX_S3D(i+2, j+1, k+2, id, jd, kd, 0);
        size_t ip7 = _F_IDX_S3D(i+2, j+2, k+2, id, jd, kd, 0);
        size_t ip8 = _F_IDX_S3D(i+1, j+2, k+2, id, jd, kd, 0);
        //size_t ip1= k   *id*jd+ j   *id+i;
        //size_t ip2= k   *id*jd+ j   *id+i+1;
        //size_t ip3= k   *id*jd+(j+1)*id+i+1;
        //size_t ip4= k   *id*jd+(j+1)*id+i;
        //size_t ip5=(k+1)*id*jd+ j   *id+i;
        //size_t ip6=(k+1)*id*jd+ j   *id+i+1;
        //size_t ip7=(k+1)*id*jd+(j+1)*id+i+1;
        //size_t ip8=(k+1)*id*jd+(j+1)*id+i;
        d[ip1]=d[ip1]+ddd;
        d[ip2]=d[ip2]+ddd;
        d[ip3]=d[ip3]+ddd;
        d[ip4]=d[ip4]+ddd;
        d[ip5]=d[ip5]+ddd;
        d[ip6]=d[ip6]+ddd;
        d[ip7]=d[ip7]+ddd;
        d[ip8]=d[ip8]+ddd;
      }
    }
  }

  //内部の格子点のデータを8で割る
  VolumeDataDivideBy8(d, id, jd, kd);

  //面上の格子点のデータを4で割る
  FaceDataDivideBy4(d, id, jd, kd);

  //辺上の格子点のデータを2で割る
  LineDataDivideBy2(d, id, jd, kd);

  //境界条件処理（wkdにガイドセルを含む格子点上のデータ）

}


//*****************************************************************************
//*****************************************************************************
//*****************************************************************************

// Scalarの格子点での値をセット（ガイドセルの値は計算対象にいれていない）
void FFV::setScalarGridData(double* d, float* data, int id, int jd, int kd)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;

  size_t dsize = (size_t)(id*jd*kd);

  for (size_t i=0; i<dsize; i++) d[i]=0.0;

  for (int km=1; km<=kx; km++) {
    for (int jm=1; jm<=jx; jm++) {
      for (int im=1; im<=ix; im++) {
        size_t mip = _F_IDX_S3D(im, jm, km, ix, jx, kx, gd);
        double ddd=(double)data[mip];
        int i=im-1;
        int j=jm-1;
        int k=km-1;
        size_t ip1 = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
        size_t ip2 = _F_IDX_S3D(i+2, j+1, k+1, id, jd, kd, 0);
        size_t ip3 = _F_IDX_S3D(i+2, j+2, k+1, id, jd, kd, 0);
        size_t ip4 = _F_IDX_S3D(i+1, j+2, k+1, id, jd, kd, 0);
        size_t ip5 = _F_IDX_S3D(i+1, j+1, k+2, id, jd, kd, 0);
        size_t ip6 = _F_IDX_S3D(i+2, j+1, k+2, id, jd, kd, 0);
        size_t ip7 = _F_IDX_S3D(i+2, j+2, k+2, id, jd, kd, 0);
        size_t ip8 = _F_IDX_S3D(i+1, j+2, k+2, id, jd, kd, 0);
        //size_t ip1= k   *id*jd+ j   *id+i;
        //size_t ip2= k   *id*jd+ j   *id+i+1;
        //size_t ip3= k   *id*jd+(j+1)*id+i+1;
        //size_t ip4= k   *id*jd+(j+1)*id+i;
        //size_t ip5=(k+1)*id*jd+ j   *id+i;
        //size_t ip6=(k+1)*id*jd+ j   *id+i+1;
        //size_t ip7=(k+1)*id*jd+(j+1)*id+i+1;
        //size_t ip8=(k+1)*id*jd+(j+1)*id+i;
        d[ip1]=d[ip1]+ddd;
        d[ip2]=d[ip2]+ddd;
        d[ip3]=d[ip3]+ddd;
        d[ip4]=d[ip4]+ddd;
        d[ip5]=d[ip5]+ddd;
        d[ip6]=d[ip6]+ddd;
        d[ip7]=d[ip7]+ddd;
        d[ip8]=d[ip8]+ddd;
      }
    }
  }

  //内部の格子点のデータを8で割る
  VolumeDataDivideBy8(d, id, jd, kd);

  //面上の格子点のデータを4で割る
  FaceDataDivideBy4(d, id, jd, kd);

  //辺上の格子点のデータを2で割る
  LineDataDivideBy2(d, id, jd, kd);

  //境界条件処理（wkdにガイドセルを含む格子点上のデータ）

}

// Vectorの格子点での値をセット（ガイドセルの値は計算対象にいれていない）
void FFV::setVectorGridData(double* d, float* data, int id, int jd, int kd)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;

  size_t dsize = (size_t)(id*jd*kd);
  size_t dsize3 = (size_t)(id*jd*kd*3);

  for (size_t i=0; i<dsize3; i++) d[i]=0.0;

  for (size_t ivar=0;ivar<3;ivar++){

    for (int km=1; km<=kx; km++) {
      for (int jm=1; jm<=jx; jm++) {
        for (int im=1; im<=ix; im++) {
          size_t mip = _F_IDX_V3DEX(ivar, im, jm, km, ix, jx, kx, gd);//ivar=0:x方向、ivar=1:y方向、ivar=2:z方向
          double ddd=(double)data[mip];
          int i=im-1;
          int j=jm-1;
          int k=km-1;
          size_t ip1 = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
          size_t ip2 = _F_IDX_S3D(i+2, j+1, k+1, id, jd, kd, 0);
          size_t ip3 = _F_IDX_S3D(i+2, j+2, k+1, id, jd, kd, 0);
          size_t ip4 = _F_IDX_S3D(i+1, j+2, k+1, id, jd, kd, 0);
          size_t ip5 = _F_IDX_S3D(i+1, j+1, k+2, id, jd, kd, 0);
          size_t ip6 = _F_IDX_S3D(i+2, j+1, k+2, id, jd, kd, 0);
          size_t ip7 = _F_IDX_S3D(i+2, j+2, k+2, id, jd, kd, 0);
          size_t ip8 = _F_IDX_S3D(i+1, j+2, k+2, id, jd, kd, 0);
          //size_t ip1= k   *id*jd+ j   *id+i;
          //size_t ip2= k   *id*jd+ j   *id+i+1;
          //size_t ip3= k   *id*jd+(j+1)*id+i+1;
          //size_t ip4= k   *id*jd+(j+1)*id+i;
          //size_t ip5=(k+1)*id*jd+ j   *id+i;
          //size_t ip6=(k+1)*id*jd+ j   *id+i+1;
          //size_t ip7=(k+1)*id*jd+(j+1)*id+i+1;
          //size_t ip8=(k+1)*id*jd+(j+1)*id+i;
          d[ip1+dsize*ivar]=d[ip1+dsize*ivar]+ddd;
          d[ip2+dsize*ivar]=d[ip2+dsize*ivar]+ddd;
          d[ip3+dsize*ivar]=d[ip3+dsize*ivar]+ddd;
          d[ip4+dsize*ivar]=d[ip4+dsize*ivar]+ddd;
          d[ip5+dsize*ivar]=d[ip5+dsize*ivar]+ddd;
          d[ip6+dsize*ivar]=d[ip6+dsize*ivar]+ddd;
          d[ip7+dsize*ivar]=d[ip7+dsize*ivar]+ddd;
          d[ip8+dsize*ivar]=d[ip8+dsize*ivar]+ddd;
        }
      }
    }

    //内部の格子点のデータを8で割る
    VolumeDataDivideBy8(&d[dsize*ivar], id, jd, kd);

    //面上の格子点のデータを4で割る
    FaceDataDivideBy4(&d[dsize*ivar], id, jd, kd);

    //辺上の格子点のデータを2で割る
    LineDataDivideBy2(&d[dsize*ivar], id, jd, kd);

    //境界条件処理（wkdにガイドセルを含む格子点上のデータ）

  }//loop ivar

}

// Vectorの格子点での値をセット（ガイドセルの値は計算対象にいれていない）
void FFV::setVectorComponentGridData(double* d, float* data, int id, int jd, int kd, int ivar)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;

  size_t dsize = (size_t)(id*jd*kd);

  for (size_t i=0; i<dsize; i++) d[i]=0.0;

  for (int km=1; km<=kx; km++) {
    for (int jm=1; jm<=jx; jm++) {
      for (int im=1; im<=ix; im++) {
        size_t mip = _F_IDX_V3DEX(ivar, im, jm, km, ix, jx, kx, gd);//ivar=0:x方向、ivar=1:y方向、ivar=2:z方向
        double ddd=(double)data[mip];
        int i=im-1;
        int j=jm-1;
        int k=km-1;
        size_t ip1 = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
        size_t ip2 = _F_IDX_S3D(i+2, j+1, k+1, id, jd, kd, 0);
        size_t ip3 = _F_IDX_S3D(i+2, j+2, k+1, id, jd, kd, 0);
        size_t ip4 = _F_IDX_S3D(i+1, j+2, k+1, id, jd, kd, 0);
        size_t ip5 = _F_IDX_S3D(i+1, j+1, k+2, id, jd, kd, 0);
        size_t ip6 = _F_IDX_S3D(i+2, j+1, k+2, id, jd, kd, 0);
        size_t ip7 = _F_IDX_S3D(i+2, j+2, k+2, id, jd, kd, 0);
        size_t ip8 = _F_IDX_S3D(i+1, j+2, k+2, id, jd, kd, 0);
        //size_t ip1= k   *id*jd+ j   *id+i;
        //size_t ip2= k   *id*jd+ j   *id+i+1;
        //size_t ip3= k   *id*jd+(j+1)*id+i+1;
        //size_t ip4= k   *id*jd+(j+1)*id+i;
        //size_t ip5=(k+1)*id*jd+ j   *id+i;
        //size_t ip6=(k+1)*id*jd+ j   *id+i+1;
        //size_t ip7=(k+1)*id*jd+(j+1)*id+i+1;
        //size_t ip8=(k+1)*id*jd+(j+1)*id+i;
        d[ip1]=d[ip1]+ddd;
        d[ip2]=d[ip2]+ddd;
        d[ip3]=d[ip3]+ddd;
        d[ip4]=d[ip4]+ddd;
        d[ip5]=d[ip5]+ddd;
        d[ip6]=d[ip6]+ddd;
        d[ip7]=d[ip7]+ddd;
        d[ip8]=d[ip8]+ddd;
      }
    }
  }

  //内部の格子点のデータを8で割る
  VolumeDataDivideBy8(d, id, jd, kd);

  //面上の格子点のデータを4で割る
  FaceDataDivideBy4(d, id, jd, kd);

  //辺上の格子点のデータを2で割る
  LineDataDivideBy2(d, id, jd, kd);

  //境界条件処理（wkdにガイドセルを含む格子点上のデータ）

}

//*****************************************************************************
//*****************************************************************************
//*****************************************************************************

//内部の格子点のデータを8で割る
void FFV::VolumeDataDivideBy8(float* d, int id, int jd, int kd)
{
  int i,j,k;
  size_t ip;

  for (k=1; k<kd-1; k++){
    for (j=1; j<jd-1; j++){
      for (i=1; i<id-1; i++){
        ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
        d[ip]=d[ip]*0.125;
      }
    }
  }
}

//面上の格子点のデータを4で割る
void FFV::FaceDataDivideBy4(float* d, int id, int jd, int kd)
{
  int i,j,k;
  size_t ip;

  i=0;
  for (k=1; k<kd-1; k++){
    for (j=1; j<jd-1; j++){
      ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
      d[ip]=d[ip]*0.25;
    }
  }

  i=id-1;
  for (k=1; k<kd-1; k++){
    for (j=1; j<jd-1; j++){
      ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
      d[ip]=d[ip]*0.25;
    }
  }

  j=0;
  for (k=1; k<kd-1; k++){
    for (i=1; i<id-1; i++){
      ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
      d[ip]=d[ip]*0.25;
    }
  }

  j=jd-1;
  for (k=1; k<kd-1; k++){
    for (i=1; i<id-1; i++){
      ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
      d[ip]=d[ip]*0.25;
    }
  }

  k=0;
  for (j=1; j<jd-1; j++){
    for (i=1; i<id-1; i++){
      ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
      d[ip]=d[ip]*0.25;
    }
  }

  k=kd-1;
  for (j=1; j<jd-1; j++){
    for (i=1; i<id-1; i++){
      ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
      d[ip]=d[ip]*0.25;
    }
  }
}

//辺上の格子点のデータを2で割る
void FFV::LineDataDivideBy2(float* d, int id, int jd, int kd)
{
  int i,j,k;
  size_t ip;

  i=0; j=0;
  for (k=1; k<kd-1; k++){
    ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
    d[ip]=d[ip]*0.5;
  }

  i=0; j=jd-1;
  for (k=1; k<kd-1; k++){
    ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
    d[ip]=d[ip]*0.5;
  }

  i=0; k=0;
  for (j=1; j<jd-1; j++){
    ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
    d[ip]=d[ip]*0.5;
  }

  i=0; k=kd-1;
  for (j=1; j<jd-1; j++){
    ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
    d[ip]=d[ip]*0.5;
  }

  j=0; k=0;
  for (i=1; i<id-1; i++){
    ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
    d[ip]=d[ip]*0.5;
  }

  j=0; k=kd-1;
  for (i=1; i<id-1; i++){
    ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
    d[ip]=d[ip]*0.5;
  }

  j=jd-1; k=0;
  for (i=1; i<id-1; i++){
    ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
    d[ip]=d[ip]*0.5;
  }

  j=jd-1; k=kd-1;
  for (i=1; i<id-1; i++){
    ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
    d[ip]=d[ip]*0.5;
  }

  i=id-1; j=0;
  for (k=1; k<kd-1; k++){
    ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
    d[ip]=d[ip]*0.5;
  }

  i=id-1; j=jd-1;
  for (k=1; k<kd-1; k++){
    ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
    d[ip]=d[ip]*0.5;
  }

  i=id-1; k=0;
  for (j=1; j<jd-1; j++){
    ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
    d[ip]=d[ip]*0.5;
  }

  i=id-1; k=kd-1;
  for (j=1; j<jd-1; j++){
    ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
    d[ip]=d[ip]*0.5;
  }

}

//*****************************************************************************
//*****************************************************************************
//*****************************************************************************

//内部の格子点のデータを8で割る
void FFV::VolumeDataDivideBy8(double* d, int id, int jd, int kd)
{
  int i,j,k;
  size_t ip;

  for (k=1; k<kd-1; k++){
    for (j=1; j<jd-1; j++){
      for (i=1; i<id-1; i++){
        ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
        d[ip]=d[ip]*0.125;
      }
    }
  }
}

//面上の格子点のデータを4で割る
void FFV::FaceDataDivideBy4(double* d, int id, int jd, int kd)
{
  int i,j,k;
  size_t ip;

  i=0;
  for (k=1; k<kd-1; k++){
    for (j=1; j<jd-1; j++){
      ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
      d[ip]=d[ip]*0.25;
    }
  }

  i=id-1;
  for (k=1; k<kd-1; k++){
    for (j=1; j<jd-1; j++){
      ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
      d[ip]=d[ip]*0.25;
    }
  }

  j=0;
  for (k=1; k<kd-1; k++){
    for (i=1; i<id-1; i++){
      ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
      d[ip]=d[ip]*0.25;
    }
  }

  j=jd-1;
  for (k=1; k<kd-1; k++){
    for (i=1; i<id-1; i++){
      ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
      d[ip]=d[ip]*0.25;
    }
  }

  k=0;
  for (j=1; j<jd-1; j++){
    for (i=1; i<id-1; i++){
      ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
      d[ip]=d[ip]*0.25;
    }
  }

  k=kd-1;
  for (j=1; j<jd-1; j++){
    for (i=1; i<id-1; i++){
      ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
      d[ip]=d[ip]*0.25;
    }
  }
}

//辺上の格子点のデータを2で割る
void FFV::LineDataDivideBy2(double* d, int id, int jd, int kd)
{
  int i,j,k;
  size_t ip;

  i=0; j=0;
  for (k=1; k<kd-1; k++){
    ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
    d[ip]=d[ip]*0.5;
  }

  i=0; j=jd-1;
  for (k=1; k<kd-1; k++){
    ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
    d[ip]=d[ip]*0.5;
  }

  i=0; k=0;
  for (j=1; j<jd-1; j++){
    ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
    d[ip]=d[ip]*0.5;
  }

  i=0; k=kd-1;
  for (j=1; j<jd-1; j++){
    ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
    d[ip]=d[ip]*0.5;
  }

  j=0; k=0;
  for (i=1; i<id-1; i++){
    ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
    d[ip]=d[ip]*0.5;
  }

  j=0; k=kd-1;
  for (i=1; i<id-1; i++){
    ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
    d[ip]=d[ip]*0.5;
  }

  j=jd-1; k=0;
  for (i=1; i<id-1; i++){
    ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
    d[ip]=d[ip]*0.5;
  }

  j=jd-1; k=kd-1;
  for (i=1; i<id-1; i++){
    ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
    d[ip]=d[ip]*0.5;
  }

  i=id-1; j=0;
  for (k=1; k<kd-1; k++){
    ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
    d[ip]=d[ip]*0.5;
  }

  i=id-1; j=jd-1;
  for (k=1; k<kd-1; k++){
    ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
    d[ip]=d[ip]*0.5;
  }

  i=id-1; k=0;
  for (j=1; j<jd-1; j++){
    ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
    d[ip]=d[ip]*0.5;
  }

  i=id-1; k=kd-1;
  for (j=1; j<jd-1; j++){
    ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
    d[ip]=d[ip]*0.5;
  }

}

//*****************************************************************************
//*****************************************************************************
//*****************************************************************************
// float ---> float

// Scalarの格子点での値をセット（ガイドセルに値があることを想定しているバージョン）
void FFV::setScalarGridDataGuide(float* d, float* data, int id, int jd, int kd, int gc_out)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;

  float *wkd;
  int iw=size[0]+2*guide+1;
  int jw=size[1]+2*guide+1;
  int kw=size[2]+2*guide+1;

  if (!(wkd = new float[ iw*jw*kw ])){
    Hostonly_  printf(    "\t>> cannot allocate work area : setScalarGridData\n\n");
    Exit(0);
  }

  for (size_t i=0; i<iw*jw*kw; i++) wkd[i]=0.0;

  for (int km=1-gd; km<=kx+gd; km++) {
    for (int jm=1-gd; jm<=jx+gd; jm++) {
      for (int im=1-gd; im<=ix+gd; im++) {
        size_t mip = _F_IDX_S3D(im, jm, km, ix, jx, kx, gd);
        float ddd=data[mip];
        int i=im-1+gd;
        int j=jm-1+gd;
        int k=km-1+gd;
        size_t ip1 = _F_IDX_S3D(i+1, j+1, k+1, iw, jw, kw, 0);
        size_t ip2 = _F_IDX_S3D(i+2, j+1, k+1, iw, jw, kw, 0);
        size_t ip3 = _F_IDX_S3D(i+2, j+2, k+1, iw, jw, kw, 0);
        size_t ip4 = _F_IDX_S3D(i+1, j+2, k+1, iw, jw, kw, 0);
        size_t ip5 = _F_IDX_S3D(i+1, j+1, k+2, iw, jw, kw, 0);
        size_t ip6 = _F_IDX_S3D(i+2, j+1, k+2, iw, jw, kw, 0);
        size_t ip7 = _F_IDX_S3D(i+2, j+2, k+2, iw, jw, kw, 0);
        size_t ip8 = _F_IDX_S3D(i+1, j+2, k+2, iw, jw, kw, 0);
        //size_t ip1= k   *iw*jw+ j   *iw+i;
        //size_t ip2= k   *iw*jw+ j   *iw+i+1;
        //size_t ip3= k   *iw*jw+(j+1)*iw+i+1;
        //size_t ip4= k   *iw*jw+(j+1)*iw+i;
        //size_t ip5=(k+1)*iw*jw+ j   *iw+i;
        //size_t ip6=(k+1)*iw*jw+ j   *iw+i+1;
        //size_t ip7=(k+1)*iw*jw+(j+1)*iw+i+1;
        //size_t ip8=(k+1)*iw*jw+(j+1)*iw+i;
        wkd[ip1]=wkd[ip1]+ddd;
        wkd[ip2]=wkd[ip2]+ddd;
        wkd[ip3]=wkd[ip3]+ddd;
        wkd[ip4]=wkd[ip4]+ddd;
        wkd[ip5]=wkd[ip5]+ddd;
        wkd[ip6]=wkd[ip6]+ddd;
        wkd[ip7]=wkd[ip7]+ddd;
        wkd[ip8]=wkd[ip8]+ddd;
      }
    }
  }

  //内部の格子点のデータを8で割る
  VolumeDataDivideBy8(wkd, iw, jw, kw);

  //面上の格子点のデータを4で割る
  FaceDataDivideBy4(wkd, iw, jw, kw);

  //辺上の格子点のデータを2で割る
  LineDataDivideBy2(wkd, iw, jw, kw);

  //境界条件処理（wkdにガイドセルを含む格子点上のデータ）

  //set d <--- wkd
  if(gc_out==0){
    for(int k=0;k<kd;k++){
      for(int j=0;j<jd;j++){
        for(int i=0;i<id;i++){
          size_t ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
          size_t ipw = _F_IDX_S3D(i+gd+1, j+gd+1, k+gd+1, iw, jw, kw, 0);
          d[ip]=wkd[ipw];
        }
      }
    }
  }
  else{
    for (size_t i=0; i<iw*jw*kw; i++) d[i]=wkd[i];
  }

  //delete [] wkd;
  if (wkd) delete [] wkd;

}

// Vectorの格子点での値をセット（ガイドセルに値があることを想定しているバージョン）
void FFV::setVectorGridDataGuide(float* d, float* data, int id, int jd, int kd, int gc_out)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;

  float *wkd;
  int iw=size[0]+2*guide+1;
  int jw=size[1]+2*guide+1;
  int kw=size[2]+2*guide+1;
  if (!(wkd = new float[ iw*jw*kw ])){
    Hostonly_  printf(    "\t>> cannot allocate work area : setVectorGridData\n\n");
    Exit(0);
  }

  for (size_t ivar=0;ivar<3;ivar++){

    for (size_t i=0; i<iw*jw*kw; i++) wkd[i]=0.0;
 
    for (int km=1-gd; km<=kx+gd; km++) {
      for (int jm=1-gd; jm<=jx+gd; jm++) {
        for (int im=1-gd; im<=ix+gd; im++) {
          size_t mip = _F_IDX_V3DEX(ivar, im, jm, km, ix, jx, kx, gd);//ivar=0:x方向、ivar=1:y方向、ivar=2:z方向
          float ddd=data[mip];
          int i=im-1+gd;
          int j=jm-1+gd;
          int k=km-1+gd;
          size_t ip1 = _F_IDX_S3D(i+1, j+1, k+1, iw, jw, kw, 0);
          size_t ip2 = _F_IDX_S3D(i+2, j+1, k+1, iw, jw, kw, 0);
          size_t ip3 = _F_IDX_S3D(i+2, j+2, k+1, iw, jw, kw, 0);
          size_t ip4 = _F_IDX_S3D(i+1, j+2, k+1, iw, jw, kw, 0);
          size_t ip5 = _F_IDX_S3D(i+1, j+1, k+2, iw, jw, kw, 0);
          size_t ip6 = _F_IDX_S3D(i+2, j+1, k+2, iw, jw, kw, 0);
          size_t ip7 = _F_IDX_S3D(i+2, j+2, k+2, iw, jw, kw, 0);
          size_t ip8 = _F_IDX_S3D(i+1, j+2, k+2, iw, jw, kw, 0);
          //size_t ip1= k   *iw*jw+ j   *iw+i;
          //size_t ip2= k   *iw*jw+ j   *iw+i+1;
          //size_t ip3= k   *iw*jw+(j+1)*iw+i+1;
          //size_t ip4= k   *iw*jw+(j+1)*iw+i;
          //size_t ip5=(k+1)*iw*jw+ j   *iw+i;
          //size_t ip6=(k+1)*iw*jw+ j   *iw+i+1;
          //size_t ip7=(k+1)*iw*jw+(j+1)*iw+i+1;
          //size_t ip8=(k+1)*iw*jw+(j+1)*iw+i;
          wkd[ip1]=wkd[ip1]+ddd;
          wkd[ip2]=wkd[ip2]+ddd;
          wkd[ip3]=wkd[ip3]+ddd;
          wkd[ip4]=wkd[ip4]+ddd;
          wkd[ip5]=wkd[ip5]+ddd;
          wkd[ip6]=wkd[ip6]+ddd;
          wkd[ip7]=wkd[ip7]+ddd;
          wkd[ip8]=wkd[ip8]+ddd;
        }
      }
    }
 
    //内部の格子点のデータを8で割る
    VolumeDataDivideBy8(wkd, iw, jw, kw);
 
    //面上の格子点のデータを4で割る
    FaceDataDivideBy4(wkd, iw, jw, kw);
 
    //辺上の格子点のデータを2で割る
    LineDataDivideBy2(wkd, iw, jw, kw);
 
    //境界条件処理（wkdにガイドセルを含む格子点上のデータ）
 
 
    //set d <--- wkd
    size_t dsize = (size_t)(id*jd*kd);
    if(gc_out==0){
      for(int k=0;k<kd;k++){
        for(int j=0;j<jd;j++){
          for(int i=0;i<id;i++){
            size_t ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
            size_t ipw = _F_IDX_S3D(i+gd+1, j+gd+1, k+gd+1, iw, jw, kw, 0);
            d[ip+dsize*ivar]=wkd[ipw];
          }
        }
      }
    }
    else{
      for (size_t i=0; i<iw*jw*kw; i++) d[i+dsize*ivar]=wkd[i];
    }

  }//loop ivar

  //delete [] wkd;
  if (wkd) delete [] wkd;
}


// Vectorの格子点での値をセット（ガイドセルに値があることを想定しているバージョン）
void FFV::setVectorComponentGridDataGuide(float* d, float* data, int id, int jd, int kd, int gc_out, int ivar)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;

  float *wkd;
  int iw=size[0]+2*guide+1;
  int jw=size[1]+2*guide+1;
  int kw=size[2]+2*guide+1;
  if (!(wkd = new float[ iw*jw*kw ])){
    Hostonly_  printf(    "\t>> cannot allocate work area : setVectorGridData\n\n");
    Exit(0);
  }

  for (size_t i=0; i<iw*jw*kw; i++) wkd[i]=0.0;
 
  for (int km=1-gd; km<=kx+gd; km++) {
    for (int jm=1-gd; jm<=jx+gd; jm++) {
      for (int im=1-gd; im<=ix+gd; im++) {
        size_t mip = _F_IDX_V3DEX(ivar, im, jm, km, ix, jx, kx, gd);//ivar=0:x方向、ivar=1:y方向、ivar=2:z方向
        float ddd=data[mip];
        int i=im-1+gd;
        int j=jm-1+gd;
        int k=km-1+gd;
        size_t ip1 = _F_IDX_S3D(i+1, j+1, k+1, iw, jw, kw, 0);
        size_t ip2 = _F_IDX_S3D(i+2, j+1, k+1, iw, jw, kw, 0);
        size_t ip3 = _F_IDX_S3D(i+2, j+2, k+1, iw, jw, kw, 0);
        size_t ip4 = _F_IDX_S3D(i+1, j+2, k+1, iw, jw, kw, 0);
        size_t ip5 = _F_IDX_S3D(i+1, j+1, k+2, iw, jw, kw, 0);
        size_t ip6 = _F_IDX_S3D(i+2, j+1, k+2, iw, jw, kw, 0);
        size_t ip7 = _F_IDX_S3D(i+2, j+2, k+2, iw, jw, kw, 0);
        size_t ip8 = _F_IDX_S3D(i+1, j+2, k+2, iw, jw, kw, 0);
        //size_t ip1= k   *iw*jw+ j   *iw+i;
        //size_t ip2= k   *iw*jw+ j   *iw+i+1;
        //size_t ip3= k   *iw*jw+(j+1)*iw+i+1;
        //size_t ip4= k   *iw*jw+(j+1)*iw+i;
        //size_t ip5=(k+1)*iw*jw+ j   *iw+i;
        //size_t ip6=(k+1)*iw*jw+ j   *iw+i+1;
        //size_t ip7=(k+1)*iw*jw+(j+1)*iw+i+1;
        //size_t ip8=(k+1)*iw*jw+(j+1)*iw+i;
        wkd[ip1]=wkd[ip1]+ddd;
        wkd[ip2]=wkd[ip2]+ddd;
        wkd[ip3]=wkd[ip3]+ddd;
        wkd[ip4]=wkd[ip4]+ddd;
        wkd[ip5]=wkd[ip5]+ddd;
        wkd[ip6]=wkd[ip6]+ddd;
        wkd[ip7]=wkd[ip7]+ddd;
        wkd[ip8]=wkd[ip8]+ddd;
      }
    }
  }
 
  //内部の格子点のデータを8で割る
  VolumeDataDivideBy8(wkd, iw, jw, kw);
 
  //面上の格子点のデータを4で割る
  FaceDataDivideBy4(wkd, iw, jw, kw);
 
  //辺上の格子点のデータを2で割る
  LineDataDivideBy2(wkd, iw, jw, kw);
 
  //境界条件処理（wkdにガイドセルを含む格子点上のデータ）
 

  //set d <--- wkd
  size_t dsize = (size_t)(id*jd*kd);
  if(gc_out==0){
    for(int k=0;k<kd;k++){
      for(int j=0;j<jd;j++){
        for(int i=0;i<id;i++){
          size_t ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
          size_t ipw = _F_IDX_S3D(i+gd+1, j+gd+1, k+gd+1, iw, jw, kw, 0);
          d[ip]=wkd[ipw];
        }
      }
    }
  }
  else{
    for (size_t i=0; i<iw*jw*kw; i++) d[i]=wkd[i];
  }

  //delete [] wkd;
  if (wkd) delete [] wkd;
}


//*****************************************************************************
//*****************************************************************************
//*****************************************************************************
// double ---> double

// Scalarの格子点での値をセット（ガイドセルに値があることを想定しているバージョン）
void FFV::setScalarGridDataGuide(double* d, double* data, int id, int jd, int kd, int gc_out)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;

  double *wkd;
  int iw=size[0]+2*guide+1;
  int jw=size[1]+2*guide+1;
  int kw=size[2]+2*guide+1;

  if (!(wkd = new double[ iw*jw*kw ])){
    Hostonly_  printf(    "\t>> cannot allocate work area : setScalarGridData\n\n");
    Exit(0);
  }

  for (size_t i=0; i<iw*jw*kw; i++) wkd[i]=0.0;

  for (int km=1-gd; km<=kx+gd; km++) {
    for (int jm=1-gd; jm<=jx+gd; jm++) {
      for (int im=1-gd; im<=ix+gd; im++) {
        size_t mip = _F_IDX_S3D(im, jm, km, ix, jx, kx, gd);
        double ddd=data[mip];
        int i=im-1+gd;
        int j=jm-1+gd;
        int k=km-1+gd;
        size_t ip1 = _F_IDX_S3D(i+1, j+1, k+1, iw, jw, kw, 0);
        size_t ip2 = _F_IDX_S3D(i+2, j+1, k+1, iw, jw, kw, 0);
        size_t ip3 = _F_IDX_S3D(i+2, j+2, k+1, iw, jw, kw, 0);
        size_t ip4 = _F_IDX_S3D(i+1, j+2, k+1, iw, jw, kw, 0);
        size_t ip5 = _F_IDX_S3D(i+1, j+1, k+2, iw, jw, kw, 0);
        size_t ip6 = _F_IDX_S3D(i+2, j+1, k+2, iw, jw, kw, 0);
        size_t ip7 = _F_IDX_S3D(i+2, j+2, k+2, iw, jw, kw, 0);
        size_t ip8 = _F_IDX_S3D(i+1, j+2, k+2, iw, jw, kw, 0);
        //size_t ip1= k   *iw*jw+ j   *iw+i;
        //size_t ip2= k   *iw*jw+ j   *iw+i+1;
        //size_t ip3= k   *iw*jw+(j+1)*iw+i+1;
        //size_t ip4= k   *iw*jw+(j+1)*iw+i;
        //size_t ip5=(k+1)*iw*jw+ j   *iw+i;
        //size_t ip6=(k+1)*iw*jw+ j   *iw+i+1;
        //size_t ip7=(k+1)*iw*jw+(j+1)*iw+i+1;
        //size_t ip8=(k+1)*iw*jw+(j+1)*iw+i;
        wkd[ip1]=wkd[ip1]+ddd;
        wkd[ip2]=wkd[ip2]+ddd;
        wkd[ip3]=wkd[ip3]+ddd;
        wkd[ip4]=wkd[ip4]+ddd;
        wkd[ip5]=wkd[ip5]+ddd;
        wkd[ip6]=wkd[ip6]+ddd;
        wkd[ip7]=wkd[ip7]+ddd;
        wkd[ip8]=wkd[ip8]+ddd;
      }
    }
  }

  //内部の格子点のデータを8で割る
  VolumeDataDivideBy8(wkd, iw, jw, kw);

  //面上の格子点のデータを4で割る
  FaceDataDivideBy4(wkd, iw, jw, kw);

  //辺上の格子点のデータを2で割る
  LineDataDivideBy2(wkd, iw, jw, kw);

  //境界条件処理（wkdにガイドセルを含む格子点上のデータ）

  //set d <--- wkd
  if(gc_out==0){
    for(int k=0;k<kd;k++){
      for(int j=0;j<jd;j++){
        for(int i=0;i<id;i++){
          size_t ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
          size_t ipw = _F_IDX_S3D(i+gd+1, j+gd+1, k+gd+1, iw, jw, kw, 0);
          d[ip]=wkd[ipw];
        }
      }
    }
  }
  else{
    for (size_t i=0; i<iw*jw*kw; i++) d[i]=wkd[i];
  }

  //delete [] wkd;
  if (wkd) delete [] wkd;

}

// Vectorの格子点での値をセット（ガイドセルに値があることを想定しているバージョン）
void FFV::setVectorGridDataGuide(double* d, double* data, int id, int jd, int kd, int gc_out)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;

  double *wkd;
  int iw=size[0]+2*guide+1;
  int jw=size[1]+2*guide+1;
  int kw=size[2]+2*guide+1;
  if (!(wkd = new double[ iw*jw*kw ])){
    Hostonly_  printf(    "\t>> cannot allocate work area : setVectorGridData\n\n");
    Exit(0);
  }

  for (size_t ivar=0;ivar<3;ivar++){

    for (size_t i=0; i<iw*jw*kw; i++) wkd[i]=0.0;
 
    for (int km=1-gd; km<=kx+gd; km++) {
      for (int jm=1-gd; jm<=jx+gd; jm++) {
        for (int im=1-gd; im<=ix+gd; im++) {
          size_t mip = _F_IDX_V3DEX(ivar, im, jm, km, ix, jx, kx, gd);//ivar=0:x方向、ivar=1:y方向、ivar=2:z方向
          double ddd=data[mip];
          int i=im-1+gd;
          int j=jm-1+gd;
          int k=km-1+gd;
          size_t ip1 = _F_IDX_S3D(i+1, j+1, k+1, iw, jw, kw, 0);
          size_t ip2 = _F_IDX_S3D(i+2, j+1, k+1, iw, jw, kw, 0);
          size_t ip3 = _F_IDX_S3D(i+2, j+2, k+1, iw, jw, kw, 0);
          size_t ip4 = _F_IDX_S3D(i+1, j+2, k+1, iw, jw, kw, 0);
          size_t ip5 = _F_IDX_S3D(i+1, j+1, k+2, iw, jw, kw, 0);
          size_t ip6 = _F_IDX_S3D(i+2, j+1, k+2, iw, jw, kw, 0);
          size_t ip7 = _F_IDX_S3D(i+2, j+2, k+2, iw, jw, kw, 0);
          size_t ip8 = _F_IDX_S3D(i+1, j+2, k+2, iw, jw, kw, 0);
          //size_t ip1= k   *iw*jw+ j   *iw+i;
          //size_t ip2= k   *iw*jw+ j   *iw+i+1;
          //size_t ip3= k   *iw*jw+(j+1)*iw+i+1;
          //size_t ip4= k   *iw*jw+(j+1)*iw+i;
          //size_t ip5=(k+1)*iw*jw+ j   *iw+i;
          //size_t ip6=(k+1)*iw*jw+ j   *iw+i+1;
          //size_t ip7=(k+1)*iw*jw+(j+1)*iw+i+1;
          //size_t ip8=(k+1)*iw*jw+(j+1)*iw+i;
          wkd[ip1]=wkd[ip1]+ddd;
          wkd[ip2]=wkd[ip2]+ddd;
          wkd[ip3]=wkd[ip3]+ddd;
          wkd[ip4]=wkd[ip4]+ddd;
          wkd[ip5]=wkd[ip5]+ddd;
          wkd[ip6]=wkd[ip6]+ddd;
          wkd[ip7]=wkd[ip7]+ddd;
          wkd[ip8]=wkd[ip8]+ddd;
        }
      }
    }
 
    //内部の格子点のデータを8で割る
    VolumeDataDivideBy8(wkd, iw, jw, kw);
 
    //面上の格子点のデータを4で割る
    FaceDataDivideBy4(wkd, iw, jw, kw);
 
    //辺上の格子点のデータを2で割る
    LineDataDivideBy2(wkd, iw, jw, kw);
 
    //境界条件処理（wkdにガイドセルを含む格子点上のデータ）
 
 
    //set d <--- wkd
    size_t dsize = (size_t)(id*jd*kd);
    if(gc_out==0){
      for(int k=0;k<kd;k++){
        for(int j=0;j<jd;j++){
          for(int i=0;i<id;i++){
            size_t ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
            size_t ipw = _F_IDX_S3D(i+gd+1, j+gd+1, k+gd+1, iw, jw, kw, 0);
            d[ip+dsize*ivar]=wkd[ipw];
          }
        }
      }
    }
    else{
      for (size_t i=0; i<iw*jw*kw; i++) d[i+dsize*ivar]=wkd[i];
    }

  }//loop ivar

  //delete [] wkd;
  if (wkd) delete [] wkd;
}


// Vectorの格子点での値をセット（ガイドセルに値があることを想定しているバージョン）
void FFV::setVectorComponentGridDataGuide(double* d, double* data, int id, int jd, int kd, int gc_out, int ivar)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;

  double *wkd;
  int iw=size[0]+2*guide+1;
  int jw=size[1]+2*guide+1;
  int kw=size[2]+2*guide+1;
  if (!(wkd = new double[ iw*jw*kw ])){
    Hostonly_  printf(    "\t>> cannot allocate work area : setVectorGridData\n\n");
    Exit(0);
  }

  for (size_t i=0; i<iw*jw*kw; i++) wkd[i]=0.0;
 
  for (int km=1-gd; km<=kx+gd; km++) {
    for (int jm=1-gd; jm<=jx+gd; jm++) {
      for (int im=1-gd; im<=ix+gd; im++) {
        size_t mip = _F_IDX_V3DEX(ivar, im, jm, km, ix, jx, kx, gd);//ivar=0:x方向、ivar=1:y方向、ivar=2:z方向
        double ddd=data[mip];
        int i=im-1+gd;
        int j=jm-1+gd;
        int k=km-1+gd;
        size_t ip1 = _F_IDX_S3D(i+1, j+1, k+1, iw, jw, kw, 0);
        size_t ip2 = _F_IDX_S3D(i+2, j+1, k+1, iw, jw, kw, 0);
        size_t ip3 = _F_IDX_S3D(i+2, j+2, k+1, iw, jw, kw, 0);
        size_t ip4 = _F_IDX_S3D(i+1, j+2, k+1, iw, jw, kw, 0);
        size_t ip5 = _F_IDX_S3D(i+1, j+1, k+2, iw, jw, kw, 0);
        size_t ip6 = _F_IDX_S3D(i+2, j+1, k+2, iw, jw, kw, 0);
        size_t ip7 = _F_IDX_S3D(i+2, j+2, k+2, iw, jw, kw, 0);
        size_t ip8 = _F_IDX_S3D(i+1, j+2, k+2, iw, jw, kw, 0);
        //size_t ip1= k   *iw*jw+ j   *iw+i;
        //size_t ip2= k   *iw*jw+ j   *iw+i+1;
        //size_t ip3= k   *iw*jw+(j+1)*iw+i+1;
        //size_t ip4= k   *iw*jw+(j+1)*iw+i;
        //size_t ip5=(k+1)*iw*jw+ j   *iw+i;
        //size_t ip6=(k+1)*iw*jw+ j   *iw+i+1;
        //size_t ip7=(k+1)*iw*jw+(j+1)*iw+i+1;
        //size_t ip8=(k+1)*iw*jw+(j+1)*iw+i;
        wkd[ip1]=wkd[ip1]+ddd;
        wkd[ip2]=wkd[ip2]+ddd;
        wkd[ip3]=wkd[ip3]+ddd;
        wkd[ip4]=wkd[ip4]+ddd;
        wkd[ip5]=wkd[ip5]+ddd;
        wkd[ip6]=wkd[ip6]+ddd;
        wkd[ip7]=wkd[ip7]+ddd;
        wkd[ip8]=wkd[ip8]+ddd;
      }
    }
  }
 
  //内部の格子点のデータを8で割る
  VolumeDataDivideBy8(wkd, iw, jw, kw);
 
  //面上の格子点のデータを4で割る
  FaceDataDivideBy4(wkd, iw, jw, kw);
 
  //辺上の格子点のデータを2で割る
  LineDataDivideBy2(wkd, iw, jw, kw);
 
  //境界条件処理（wkdにガイドセルを含む格子点上のデータ）
 



  //set d <--- wkd
  size_t dsize = (size_t)(id*jd*kd);
  if(gc_out==0){
    for(int k=0;k<kd;k++){
      for(int j=0;j<jd;j++){
        for(int i=0;i<id;i++){
          size_t ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
          size_t ipw = _F_IDX_S3D(i+gd+1, j+gd+1, k+gd+1, iw, jw, kw, 0);
          d[ip]=wkd[ipw];
        }
      }
    }
  }
  else{
    for (size_t i=0; i<iw*jw*kw; i++) d[i]=wkd[i];
  }

  //delete [] wkd;
  if (wkd) delete [] wkd;
}


//*****************************************************************************
//*****************************************************************************
//*****************************************************************************
// float ---> double

// Scalarの格子点での値をセット（ガイドセルに値があることを想定しているバージョン）
void FFV::setScalarGridDataGuide(double* d, float* data, int id, int jd, int kd, int gc_out)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;

  float *wkd;
  int iw=size[0]+2*guide+1;
  int jw=size[1]+2*guide+1;
  int kw=size[2]+2*guide+1;

  if (!(wkd = new float[ iw*jw*kw ])){
    Hostonly_  printf(    "\t>> cannot allocate work area : setScalarGridData\n\n");
    Exit(0);
  }

  for (size_t i=0; i<iw*jw*kw; i++) wkd[i]=0.0;

  for (int km=1-gd; km<=kx+gd; km++) {
    for (int jm=1-gd; jm<=jx+gd; jm++) {
      for (int im=1-gd; im<=ix+gd; im++) {
        size_t mip = _F_IDX_S3D(im, jm, km, ix, jx, kx, gd);
        float ddd=data[mip];
        int i=im-1+gd;
        int j=jm-1+gd;
        int k=km-1+gd;
        size_t ip1 = _F_IDX_S3D(i+1, j+1, k+1, iw, jw, kw, 0);
        size_t ip2 = _F_IDX_S3D(i+2, j+1, k+1, iw, jw, kw, 0);
        size_t ip3 = _F_IDX_S3D(i+2, j+2, k+1, iw, jw, kw, 0);
        size_t ip4 = _F_IDX_S3D(i+1, j+2, k+1, iw, jw, kw, 0);
        size_t ip5 = _F_IDX_S3D(i+1, j+1, k+2, iw, jw, kw, 0);
        size_t ip6 = _F_IDX_S3D(i+2, j+1, k+2, iw, jw, kw, 0);
        size_t ip7 = _F_IDX_S3D(i+2, j+2, k+2, iw, jw, kw, 0);
        size_t ip8 = _F_IDX_S3D(i+1, j+2, k+2, iw, jw, kw, 0);
        //size_t ip1= k   *iw*jw+ j   *iw+i;
        //size_t ip2= k   *iw*jw+ j   *iw+i+1;
        //size_t ip3= k   *iw*jw+(j+1)*iw+i+1;
        //size_t ip4= k   *iw*jw+(j+1)*iw+i;
        //size_t ip5=(k+1)*iw*jw+ j   *iw+i;
        //size_t ip6=(k+1)*iw*jw+ j   *iw+i+1;
        //size_t ip7=(k+1)*iw*jw+(j+1)*iw+i+1;
        //size_t ip8=(k+1)*iw*jw+(j+1)*iw+i;
        wkd[ip1]=wkd[ip1]+ddd;
        wkd[ip2]=wkd[ip2]+ddd;
        wkd[ip3]=wkd[ip3]+ddd;
        wkd[ip4]=wkd[ip4]+ddd;
        wkd[ip5]=wkd[ip5]+ddd;
        wkd[ip6]=wkd[ip6]+ddd;
        wkd[ip7]=wkd[ip7]+ddd;
        wkd[ip8]=wkd[ip8]+ddd;
      }
    }
  }

  //内部の格子点のデータを8で割る
  VolumeDataDivideBy8(wkd, iw, jw, kw);

  //面上の格子点のデータを4で割る
  FaceDataDivideBy4(wkd, iw, jw, kw);

  //辺上の格子点のデータを2で割る
  LineDataDivideBy2(wkd, iw, jw, kw);

  //境界条件処理（wkdにガイドセルを含む格子点上のデータ）

  //set d <--- wkd
  if(gc_out==0){
    for(int k=0;k<kd;k++){
      for(int j=0;j<jd;j++){
        for(int i=0;i<id;i++){
          size_t ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
          size_t ipw = _F_IDX_S3D(i+gd+1, j+gd+1, k+gd+1, iw, jw, kw, 0);
          d[ip]=(double)wkd[ipw];
        }
      }
    }
  }
  else{
    for (size_t i=0; i<iw*jw*kw; i++) d[i]=(double)wkd[i];
  }

  //delete [] wkd;
  if (wkd) delete [] wkd;

}

// Vectorの格子点での値をセット（ガイドセルに値があることを想定しているバージョン）
void FFV::setVectorGridDataGuide(double* d, float* data, int id, int jd, int kd, int gc_out)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;

  float *wkd;
  int iw=size[0]+2*guide+1;
  int jw=size[1]+2*guide+1;
  int kw=size[2]+2*guide+1;
  if (!(wkd = new float[ iw*jw*kw ])){
    Hostonly_  printf(    "\t>> cannot allocate work area : setVectorGridData\n\n");
    Exit(0);
  }

  for (size_t ivar=0;ivar<3;ivar++){

    for (size_t i=0; i<iw*jw*kw; i++) wkd[i]=0.0;
 
    for (int km=1-gd; km<=kx+gd; km++) {
      for (int jm=1-gd; jm<=jx+gd; jm++) {
        for (int im=1-gd; im<=ix+gd; im++) {
          size_t mip = _F_IDX_V3DEX(ivar, im, jm, km, ix, jx, kx, gd);//ivar=0:x方向、ivar=1:y方向、ivar=2:z方向
          float ddd=data[mip];
          int i=im-1+gd;
          int j=jm-1+gd;
          int k=km-1+gd;
          size_t ip1 = _F_IDX_S3D(i+1, j+1, k+1, iw, jw, kw, 0);
          size_t ip2 = _F_IDX_S3D(i+2, j+1, k+1, iw, jw, kw, 0);
          size_t ip3 = _F_IDX_S3D(i+2, j+2, k+1, iw, jw, kw, 0);
          size_t ip4 = _F_IDX_S3D(i+1, j+2, k+1, iw, jw, kw, 0);
          size_t ip5 = _F_IDX_S3D(i+1, j+1, k+2, iw, jw, kw, 0);
          size_t ip6 = _F_IDX_S3D(i+2, j+1, k+2, iw, jw, kw, 0);
          size_t ip7 = _F_IDX_S3D(i+2, j+2, k+2, iw, jw, kw, 0);
          size_t ip8 = _F_IDX_S3D(i+1, j+2, k+2, iw, jw, kw, 0);
          //size_t ip1= k   *iw*jw+ j   *iw+i;
          //size_t ip2= k   *iw*jw+ j   *iw+i+1;
          //size_t ip3= k   *iw*jw+(j+1)*iw+i+1;
          //size_t ip4= k   *iw*jw+(j+1)*iw+i;
          //size_t ip5=(k+1)*iw*jw+ j   *iw+i;
          //size_t ip6=(k+1)*iw*jw+ j   *iw+i+1;
          //size_t ip7=(k+1)*iw*jw+(j+1)*iw+i+1;
          //size_t ip8=(k+1)*iw*jw+(j+1)*iw+i;
          wkd[ip1]=wkd[ip1]+ddd;
          wkd[ip2]=wkd[ip2]+ddd;
          wkd[ip3]=wkd[ip3]+ddd;
          wkd[ip4]=wkd[ip4]+ddd;
          wkd[ip5]=wkd[ip5]+ddd;
          wkd[ip6]=wkd[ip6]+ddd;
          wkd[ip7]=wkd[ip7]+ddd;
          wkd[ip8]=wkd[ip8]+ddd;
        }
      }
    }
 
    //内部の格子点のデータを8で割る
    VolumeDataDivideBy8(wkd, iw, jw, kw);
 
    //面上の格子点のデータを4で割る
    FaceDataDivideBy4(wkd, iw, jw, kw);
 
    //辺上の格子点のデータを2で割る
    LineDataDivideBy2(wkd, iw, jw, kw);
 
    //境界条件処理（wkdにガイドセルを含む格子点上のデータ）
 
 
    //set d <--- wkd
    size_t dsize = (size_t)(id*jd*kd);
    if(gc_out==0){
      for(int k=0;k<kd;k++){
        for(int j=0;j<jd;j++){
          for(int i=0;i<id;i++){
            size_t ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
            size_t ipw = _F_IDX_S3D(i+gd+1, j+gd+1, k+gd+1, iw, jw, kw, 0);
            d[ip+dsize*ivar]=(double)wkd[ipw];
          }
        }
      }
    }
    else{
      for (size_t i=0; i<iw*jw*kw; i++) d[i+dsize*ivar]=(double)wkd[i];
    }

  }//loop ivar

  //delete [] wkd;
  if (wkd) delete [] wkd;
}


// Vectorの格子点での値をセット（ガイドセルに値があることを想定しているバージョン）
void FFV::setVectorComponentGridDataGuide(double* d, float* data, int id, int jd, int kd, int gc_out, int ivar)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;

  float *wkd;
  int iw=size[0]+2*guide+1;
  int jw=size[1]+2*guide+1;
  int kw=size[2]+2*guide+1;
  if (!(wkd = new float[ iw*jw*kw ])){
    Hostonly_  printf(    "\t>> cannot allocate work area : setVectorGridData\n\n");
    Exit(0);
  }

  for (size_t i=0; i<iw*jw*kw; i++) wkd[i]=0.0;
 
  for (int km=1-gd; km<=kx+gd; km++) {
    for (int jm=1-gd; jm<=jx+gd; jm++) {
      for (int im=1-gd; im<=ix+gd; im++) {
        size_t mip = _F_IDX_V3DEX(ivar, im, jm, km, ix, jx, kx, gd);//ivar=0:x方向、ivar=1:y方向、ivar=2:z方向
        float ddd=data[mip];
        int i=im-1+gd;
        int j=jm-1+gd;
        int k=km-1+gd;
        size_t ip1 = _F_IDX_S3D(i+1, j+1, k+1, iw, jw, kw, 0);
        size_t ip2 = _F_IDX_S3D(i+2, j+1, k+1, iw, jw, kw, 0);
        size_t ip3 = _F_IDX_S3D(i+2, j+2, k+1, iw, jw, kw, 0);
        size_t ip4 = _F_IDX_S3D(i+1, j+2, k+1, iw, jw, kw, 0);
        size_t ip5 = _F_IDX_S3D(i+1, j+1, k+2, iw, jw, kw, 0);
        size_t ip6 = _F_IDX_S3D(i+2, j+1, k+2, iw, jw, kw, 0);
        size_t ip7 = _F_IDX_S3D(i+2, j+2, k+2, iw, jw, kw, 0);
        size_t ip8 = _F_IDX_S3D(i+1, j+2, k+2, iw, jw, kw, 0);
        //size_t ip1= k   *iw*jw+ j   *iw+i;
        //size_t ip2= k   *iw*jw+ j   *iw+i+1;
        //size_t ip3= k   *iw*jw+(j+1)*iw+i+1;
        //size_t ip4= k   *iw*jw+(j+1)*iw+i;
        //size_t ip5=(k+1)*iw*jw+ j   *iw+i;
        //size_t ip6=(k+1)*iw*jw+ j   *iw+i+1;
        //size_t ip7=(k+1)*iw*jw+(j+1)*iw+i+1;
        //size_t ip8=(k+1)*iw*jw+(j+1)*iw+i;
        wkd[ip1]=wkd[ip1]+ddd;
        wkd[ip2]=wkd[ip2]+ddd;
        wkd[ip3]=wkd[ip3]+ddd;
        wkd[ip4]=wkd[ip4]+ddd;
        wkd[ip5]=wkd[ip5]+ddd;
        wkd[ip6]=wkd[ip6]+ddd;
        wkd[ip7]=wkd[ip7]+ddd;
        wkd[ip8]=wkd[ip8]+ddd;
      }
    }
  }
 
  //内部の格子点のデータを8で割る
  VolumeDataDivideBy8(wkd, iw, jw, kw);
 
  //面上の格子点のデータを4で割る
  FaceDataDivideBy4(wkd, iw, jw, kw);
 
  //辺上の格子点のデータを2で割る
  LineDataDivideBy2(wkd, iw, jw, kw);
 
  //境界条件処理（wkdにガイドセルを含む格子点上のデータ）
 



  //set d <--- wkd
  size_t dsize = (size_t)(id*jd*kd);
  if(gc_out==0){
    for(int k=0;k<kd;k++){
      for(int j=0;j<jd;j++){
        for(int i=0;i<id;i++){
          size_t ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
          size_t ipw = _F_IDX_S3D(i+gd+1, j+gd+1, k+gd+1, iw, jw, kw, 0);
          d[ip]=(double)wkd[ipw];
        }
      }
    }
  }
  else{
    for (size_t i=0; i<iw*jw*kw; i++) d[i]=(double)wkd[i];
  }

  //delete [] wkd;
  if (wkd) delete [] wkd;
}


//*****************************************************************************
//*****************************************************************************
//*****************************************************************************
// double ---> float

// Scalarの格子点での値をセット（ガイドセルに値があることを想定しているバージョン）
void FFV::setScalarGridDataGuide(float* d, double* data, int id, int jd, int kd, int gc_out)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;

  double *wkd;
  int iw=size[0]+2*guide+1;
  int jw=size[1]+2*guide+1;
  int kw=size[2]+2*guide+1;

  if (!(wkd = new double[ iw*jw*kw ])){
    Hostonly_  printf(    "\t>> cannot allocate work area : setScalarGridData\n\n");
    Exit(0);
  }

  for (size_t i=0; i<iw*jw*kw; i++) wkd[i]=0.0;

  for (int km=1-gd; km<=kx+gd; km++) {
    for (int jm=1-gd; jm<=jx+gd; jm++) {
      for (int im=1-gd; im<=ix+gd; im++) {
        size_t mip = _F_IDX_S3D(im, jm, km, ix, jx, kx, gd);
        double ddd=data[mip];
        int i=im-1+gd;
        int j=jm-1+gd;
        int k=km-1+gd;
        size_t ip1 = _F_IDX_S3D(i+1, j+1, k+1, iw, jw, kw, 0);
        size_t ip2 = _F_IDX_S3D(i+2, j+1, k+1, iw, jw, kw, 0);
        size_t ip3 = _F_IDX_S3D(i+2, j+2, k+1, iw, jw, kw, 0);
        size_t ip4 = _F_IDX_S3D(i+1, j+2, k+1, iw, jw, kw, 0);
        size_t ip5 = _F_IDX_S3D(i+1, j+1, k+2, iw, jw, kw, 0);
        size_t ip6 = _F_IDX_S3D(i+2, j+1, k+2, iw, jw, kw, 0);
        size_t ip7 = _F_IDX_S3D(i+2, j+2, k+2, iw, jw, kw, 0);
        size_t ip8 = _F_IDX_S3D(i+1, j+2, k+2, iw, jw, kw, 0);
        //size_t ip1= k   *iw*jw+ j   *iw+i;
        //size_t ip2= k   *iw*jw+ j   *iw+i+1;
        //size_t ip3= k   *iw*jw+(j+1)*iw+i+1;
        //size_t ip4= k   *iw*jw+(j+1)*iw+i;
        //size_t ip5=(k+1)*iw*jw+ j   *iw+i;
        //size_t ip6=(k+1)*iw*jw+ j   *iw+i+1;
        //size_t ip7=(k+1)*iw*jw+(j+1)*iw+i+1;
        //size_t ip8=(k+1)*iw*jw+(j+1)*iw+i;
        wkd[ip1]=wkd[ip1]+ddd;
        wkd[ip2]=wkd[ip2]+ddd;
        wkd[ip3]=wkd[ip3]+ddd;
        wkd[ip4]=wkd[ip4]+ddd;
        wkd[ip5]=wkd[ip5]+ddd;
        wkd[ip6]=wkd[ip6]+ddd;
        wkd[ip7]=wkd[ip7]+ddd;
        wkd[ip8]=wkd[ip8]+ddd;
      }
    }
  }

  //内部の格子点のデータを8で割る
  VolumeDataDivideBy8(wkd, iw, jw, kw);

  //面上の格子点のデータを4で割る
  FaceDataDivideBy4(wkd, iw, jw, kw);

  //辺上の格子点のデータを2で割る
  LineDataDivideBy2(wkd, iw, jw, kw);

  //境界条件処理（wkdにガイドセルを含む格子点上のデータ）

  //set d <--- wkd
  if(gc_out==0){
    for(int k=0;k<kd;k++){
      for(int j=0;j<jd;j++){
        for(int i=0;i<id;i++){
          size_t ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
          size_t ipw = _F_IDX_S3D(i+gd+1, j+gd+1, k+gd+1, iw, jw, kw, 0);
          d[ip]=(float)wkd[ipw];
        }
      }
    }
  }
  else{
    for (size_t i=0; i<iw*jw*kw; i++) d[i]=(float)wkd[i];
  }

  //delete [] wkd;
  if (wkd) delete [] wkd;

}

// Vectorの格子点での値をセット（ガイドセルに値があることを想定しているバージョン）
void FFV::setVectorGridDataGuide(float* d, double* data, int id, int jd, int kd, int gc_out)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;

  double *wkd;
  int iw=size[0]+2*guide+1;
  int jw=size[1]+2*guide+1;
  int kw=size[2]+2*guide+1;
  if (!(wkd = new double[ iw*jw*kw ])){
    Hostonly_  printf(    "\t>> cannot allocate work area : setVectorGridData\n\n");
    Exit(0);
  }

  for (size_t ivar=0;ivar<3;ivar++){

    for (size_t i=0; i<iw*jw*kw; i++) wkd[i]=0.0;
 
    for (int km=1-gd; km<=kx+gd; km++) {
      for (int jm=1-gd; jm<=jx+gd; jm++) {
        for (int im=1-gd; im<=ix+gd; im++) {
          size_t mip = _F_IDX_V3DEX(ivar, im, jm, km, ix, jx, kx, gd);//ivar=0:x方向、ivar=1:y方向、ivar=2:z方向
          double ddd=data[mip];
          int i=im-1+gd;
          int j=jm-1+gd;
          int k=km-1+gd;
          size_t ip1= k   *iw*jw+ j   *iw+i;
          size_t ip2= k   *iw*jw+ j   *iw+i+1;
          size_t ip3= k   *iw*jw+(j+1)*iw+i+1;
          size_t ip4= k   *iw*jw+(j+1)*iw+i;
          size_t ip5=(k+1)*iw*jw+ j   *iw+i;
          size_t ip6=(k+1)*iw*jw+ j   *iw+i+1;
          size_t ip7=(k+1)*iw*jw+(j+1)*iw+i+1;
          size_t ip8=(k+1)*iw*jw+(j+1)*iw+i;
          wkd[ip1]=wkd[ip1]+ddd;
          wkd[ip2]=wkd[ip2]+ddd;
          wkd[ip3]=wkd[ip3]+ddd;
          wkd[ip4]=wkd[ip4]+ddd;
          wkd[ip5]=wkd[ip5]+ddd;
          wkd[ip6]=wkd[ip6]+ddd;
          wkd[ip7]=wkd[ip7]+ddd;
          wkd[ip8]=wkd[ip8]+ddd;
        }
      }
    }
 
    //内部の格子点のデータを8で割る
    VolumeDataDivideBy8(wkd, iw, jw, kw);
 
    //面上の格子点のデータを4で割る
    FaceDataDivideBy4(wkd, iw, jw, kw);
 
    //辺上の格子点のデータを2で割る
    LineDataDivideBy2(wkd, iw, jw, kw);
 
    //境界条件処理（wkdにガイドセルを含む格子点上のデータ）
 
 
    //set d <--- wkd
    size_t dsize = (size_t)(id*jd*kd);
    if(gc_out==0){
      for(int k=0;k<kd;k++){
        for(int j=0;j<jd;j++){
          for(int i=0;i<id;i++){
            size_t ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
            size_t ipw = _F_IDX_S3D(i+gd+1, j+gd+1, k+gd+1, iw, jw, kw, 0);
            d[ip+dsize*ivar]=(float)wkd[ipw];
          }
        }
      }
    }
    else{
      for (size_t i=0; i<iw*jw*kw; i++) d[i+dsize*ivar]=(float)wkd[i];
    }

  }//loop ivar

  //delete [] wkd;
  if (wkd) delete [] wkd;
}


// Vectorの格子点での値をセット（ガイドセルに値があることを想定しているバージョン）
void FFV::setVectorComponentGridDataGuide(float* d, double* data, int id, int jd, int kd, int gc_out, int ivar)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;

  double *wkd;
  int iw=size[0]+2*guide+1;
  int jw=size[1]+2*guide+1;
  int kw=size[2]+2*guide+1;
  if (!(wkd = new double[ iw*jw*kw ])){
    Hostonly_  printf(    "\t>> cannot allocate work area : setVectorGridData\n\n");
    Exit(0);
  }

  for (size_t i=0; i<iw*jw*kw; i++) wkd[i]=0.0;
 
  for (int km=1-gd; km<=kx+gd; km++) {
    for (int jm=1-gd; jm<=jx+gd; jm++) {
      for (int im=1-gd; im<=ix+gd; im++) {
        size_t mip = _F_IDX_V3DEX(ivar, im, jm, km, ix, jx, kx, gd);//ivar=0:x方向、ivar=1:y方向、ivar=2:z方向
        double ddd=data[mip];
        int i=im-1+gd;
        int j=jm-1+gd;
        int k=km-1+gd;
        size_t ip1 = _F_IDX_S3D(i+1, j+1, k+1, iw, jw, kw, 0);
        size_t ip2 = _F_IDX_S3D(i+2, j+1, k+1, iw, jw, kw, 0);
        size_t ip3 = _F_IDX_S3D(i+2, j+2, k+1, iw, jw, kw, 0);
        size_t ip4 = _F_IDX_S3D(i+1, j+2, k+1, iw, jw, kw, 0);
        size_t ip5 = _F_IDX_S3D(i+1, j+1, k+2, iw, jw, kw, 0);
        size_t ip6 = _F_IDX_S3D(i+2, j+1, k+2, iw, jw, kw, 0);
        size_t ip7 = _F_IDX_S3D(i+2, j+2, k+2, iw, jw, kw, 0);
        size_t ip8 = _F_IDX_S3D(i+1, j+2, k+2, iw, jw, kw, 0);
        //size_t ip1= k   *iw*jw+ j   *iw+i;
        //size_t ip2= k   *iw*jw+ j   *iw+i+1;
        //size_t ip3= k   *iw*jw+(j+1)*iw+i+1;
        //size_t ip4= k   *iw*jw+(j+1)*iw+i;
        //size_t ip5=(k+1)*iw*jw+ j   *iw+i;
        //size_t ip6=(k+1)*iw*jw+ j   *iw+i+1;
        //size_t ip7=(k+1)*iw*jw+(j+1)*iw+i+1;
        //size_t ip8=(k+1)*iw*jw+(j+1)*iw+i;
        wkd[ip1]=wkd[ip1]+ddd;
        wkd[ip2]=wkd[ip2]+ddd;
        wkd[ip3]=wkd[ip3]+ddd;
        wkd[ip4]=wkd[ip4]+ddd;
        wkd[ip5]=wkd[ip5]+ddd;
        wkd[ip6]=wkd[ip6]+ddd;
        wkd[ip7]=wkd[ip7]+ddd;
        wkd[ip8]=wkd[ip8]+ddd;
      }
    }
  }
 
  //内部の格子点のデータを8で割る
  VolumeDataDivideBy8(wkd, iw, jw, kw);
 
  //面上の格子点のデータを4で割る
  FaceDataDivideBy4(wkd, iw, jw, kw);
 
  //辺上の格子点のデータを2で割る
  LineDataDivideBy2(wkd, iw, jw, kw);
 
  //境界条件処理（wkdにガイドセルを含む格子点上のデータ）
 



  //set d <--- wkd
  size_t dsize = (size_t)(id*jd*kd);
  if(gc_out==0){
    for(int k=0;k<kd;k++){
      for(int j=0;j<jd;j++){
        for(int i=0;i<id;i++){
          size_t ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
          size_t ipw = _F_IDX_S3D(i+gd+1, j+gd+1, k+gd+1, iw, jw, kw, 0);
          d[ip]=(float)wkd[ipw];
        }
      }
    }
  }
  else{
    for (size_t i=0; i<iw*jw*kw; i++) d[i]=(float)wkd[i];
  }

  //delete [] wkd;
  if (wkd) delete [] wkd;
}


// Iblankのセット（ガイドセルに値があることを想定しているバージョン）
void FFV::setIblankGuide(int* iblank, int id, int jd, int kd)
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

  //障害物セルを構成する格子点すべてiblank=0
  //int target_id = C.Fill_Medium;
  //cout << "target_id = " << target_id << endl;

  for (int km=1-gd; km<=kx+gd; km++) {
    for (int jm=1-gd; jm<=jx+gd; jm++) {
      for (int im=1-gd; im<=ix+gd; im++) {
        mip = _F_IDX_S3D(im, jm, km, ix, jx, kx, gd);
        s = d_bcd[mip];
        //if(d_mid[mip]==1){//流体であれば
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
          //size_t ip1= k   *id*jd+ j   *id+i;
          //size_t ip2= k   *id*jd+ j   *id+i+1;
          //size_t ip3= k   *id*jd+(j+1)*id+i+1;
          //size_t ip4= k   *id*jd+(j+1)*id+i;
          //size_t ip5=(k+1)*id*jd+ j   *id+i;
          //size_t ip6=(k+1)*id*jd+ j   *id+i+1;
          //size_t ip7=(k+1)*id*jd+(j+1)*id+i+1;
          //size_t ip8=(k+1)*id*jd+(j+1)*id+i;
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

#if 0

  cout << endl;
  cout << "iblank debug" << endl;
  for(k=0;k<kd;k++){
    for(j=0;j<jd;j++){
      //cout << "j = " << j << " k = " << k << endl;
      for(i=0;i<id;i++){
        ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
        //cout << "iblank[" << ip << "] = " << iblank[ip];
        cout << " " << iblank[ip];
      }
        cout << endl;
    }
      cout << endl;
  }

#endif

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

  //障害物セルを構成する格子点すべてiblank=0
  //int target_id = C.Fill_Medium;
  //cout << "target_id = " << target_id << endl;

  for (int km=1-gd; km<=kx+gd; km++) {
    for (int jm=1-gd; jm<=jx+gd; jm++) {
      for (int im=1-gd; im<=ix+gd; im++) {
        mip = _F_IDX_S3D(im, jm, km, ix, jx, kx, gd);
        s = d_bcd[mip];
        //if(d_mid[mip]==id_of_solid){//固体であれば
        //if(d_mid[mip]!=target_id){//流体でなければ
        //if(d_mid[mip]!=1){//流体でなければ
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
          //size_t ip1= k   *id*jd+ j   *id+i;
          //size_t ip2= k   *id*jd+ j   *id+i+1;
          //size_t ip3= k   *id*jd+(j+1)*id+i+1;
          //size_t ip4= k   *id*jd+(j+1)*id+i;
          //size_t ip5=(k+1)*id*jd+ j   *id+i;
          //size_t ip6=(k+1)*id*jd+ j   *id+i+1;
          //size_t ip7=(k+1)*id*jd+(j+1)*id+i+1;
          //size_t ip8=(k+1)*id*jd+(j+1)*id+i;
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
    //ipxp= k   *id*jd+ j   *id+i+1;
    //ipyp= k   *id*jd+(j+1)*id+i;
    //ipzp=(k+1)*id*jd+ j   *id+i;
    //ipcr=(k+1)*id*jd+(j+1)*id+i+1;
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
    //ipxm= k   *id*jd+ j   *id+i-1;
    //ipyp= k   *id*jd+(j+1)*id+i;
    //ipzp=(k+1)*id*jd+ j   *id+i;
    //ipcr=(k+1)*id*jd+(j+1)*id+i-1;
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
    //ipxm= k   *id*jd+ j   *id+i-1;
    //ipym= k   *id*jd+(j-1)*id+i;
    //ipzp=(k+1)*id*jd+ j   *id+i;
    //ipcr=(k+1)*id*jd+(j-1)*id+i-1;
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
    //ipxp= k   *id*jd+ j   *id+i+1;
    //ipym= k   *id*jd+(j-1)*id+i;
    //ipzp=(k+1)*id*jd+ j   *id+i;
    //ipcr=(k+1)*id*jd+(j-1)*id+i+1;
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
    //ipxp= k   *id*jd+ j   *id+i+1;
    //ipyp= k   *id*jd+(j+1)*id+i;
    //ipzm=(k-1)*id*jd+ j   *id+i;
    //ipcr=(k-1)*id*jd+(j+1)*id+i+1;
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
    //ipxm= k   *id*jd+ j   *id+i-1;
    //ipyp= k   *id*jd+(j+1)*id+i;
    //ipzm=(k-1)*id*jd+ j   *id+i;
    //ipcr=(k-1)*id*jd+(j+1)*id+i-1;
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
    //ipxm= k   *id*jd+ j   *id+i-1;
    //ipym= k   *id*jd+(j-1)*id+i;
    //ipzm=(k-1)*id*jd+ j   *id+i;
    //ipcr=(k-1)*id*jd+(j-1)*id+i-1;
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
    //ipxp= k   *id*jd+ j   *id+i+1;
    //ipym= k   *id*jd+(j-1)*id+i;
    //ipzm=(k-1)*id*jd+ j   *id+i;
    //ipcr=(k-1)*id*jd+(j-1)*id+i-1;
    if( iblank[ipxp]==1 || iblank[ipym]==1 || iblank[ipzm]==1 || iblank[ipcr]==1){
      iblank[ip]=2;
    }
  }

#endif

#if 0

  cout << endl;
  cout << "iblank debug" << endl;
  for(k=0;k<kd;k++){
    for(j=0;j<jd;j++){
      //cout << "j = " << j << " k = " << k << endl;
      for(i=0;i<id;i++){
        ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
        //cout << "iblank[" << ip << "] = " << iblank[ip];
        cout << " " << iblank[ip];
      }
        cout << endl;
    }
      cout << endl;
  }

#endif

}
