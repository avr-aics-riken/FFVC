// #################################################################
//
// FFV : Frontflow / violet Cartesian
//
// Copyright (c) All right reserved. 2012
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
  if (C.Mode.VRT == ON ) nvar++;
  //if (C.Mode.I2VGT == ON ) nvar++;
  if (C.Mode.Helicity == ON ) nvar++;
  C.P3Op.nvar=nvar;
  
}

//
void FFV::OutputPlot3D_xyz(const bool restart)
{
  //value
  int igrid;
  int ngrid;
  int *id,*jd,*kd;
  REAL_TYPE *x,*y,*z;
  int *iblank;
  
  //allocate
  ngrid=C.P3Op.ngrid;
  id = new int[ngrid];
  jd = new int[ngrid];
  kd = new int[ngrid];
  
  //set grid data
  //for(igrid=0;igrid<ngrid;igrid++;){ //--->BCMでループが必要になる？
  
  igrid=0;
  id[igrid]=size[0]+2*guide;//将来的にGRIDはBCMのブロックになる？
  jd[igrid]=size[1]+2*guide;
  kd[igrid]=size[2]+2*guide;
  
  //}//igrid loop
  
  
  //set filename
  
  // 出力ファイル名
  std::string tmp;
  
  // 並列出力モード
  bool pout = ( C.FIO.IO_Output == IO_GATHER ) ? false : true;
  
  tmp = DFI.Generate_FileName_Free(C.P3Op.basename, "xyz", 0, myRank, pout);
  if ( restart ) tmp = "restart_" + tmp; // リスタート用
  
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
  
  //for(igrid=0;igrid<ngrid;igrid++;){ //--->BCMでループが必要になる？
  
  igrid=0;//igrid=0
  
  //set xyz
  x = new REAL_TYPE[ id[igrid]*jd[igrid]*kd[igrid] ];
  y = new REAL_TYPE[ id[igrid]*jd[igrid]*kd[igrid] ];
  z = new REAL_TYPE[ id[igrid]*jd[igrid]*kd[igrid] ];
  for(int k=0;k<kd[igrid];k++){
    for(int j=0;j<jd[igrid];j++){
      for(int i=0;i<id[igrid];i++){
        int ip=k*id[igrid]*jd[igrid]+j*id[igrid]+i;
        x[ip]=origin[0]+pitch[0]*(REAL_TYPE)i-pitch[0]*(REAL_TYPE)guide;
        y[ip]=origin[1]+pitch[1]*(REAL_TYPE)j-pitch[1]*(REAL_TYPE)guide;
        z[ip]=origin[2]+pitch[2]*(REAL_TYPE)k-pitch[2]*(REAL_TYPE)guide;
      }
    }
  }
  
  //set iblank
  if(FP3DW.IsIBlankFlag()){
    iblank = new int[ id[igrid]*jd[igrid]*kd[igrid] ];
    for(int k=0;k<kd[igrid];k++){
      for(int j=0;j<jd[igrid];j++){
        for(int i=0;i<id[igrid];i++){
          int ip=k*id[igrid]*jd[igrid]+j*id[igrid]+i;
          iblank[ip]=1;
        }
      }
    }
  }
  
  //write
  FP3DW.setGridData(id[igrid],jd[igrid],kd[igrid],ngrid);
  FP3DW.setXYZData(x,y,z,iblank);
  if(!FP3DW.WriteXYZData()) std::cout << "error WriteXYZData" << std::endl;
  delete [] x;
  delete [] y;
  delete [] z;
  if(FP3DW.IsIBlankFlag()){
    delete [] iblank;
  }
  
  //}//igrid loop
  
  //close file
  FP3DW.CloseFile();
  
}

//
void FFV::OutputPlot3D_post(double& flop, const bool restart)
{
  if(C.P3Op.IS_q == ON) OutputPlot3D_q(flop, restart);
  if(C.P3Op.IS_funciton == ON) OutputPlot3D_function(flop, restart);
}

//
void FFV::OutputPlot3D_q(double& flop, const bool restart)
{
  
}


//
void FFV::OutputPlot3D_function(double& flop, const bool restart)
{
  //value
  int igrid;
  int ngrid;
  int *id,*jd,*kd;
  REAL_TYPE *x,*y,*z;
  int *nvar;
  REAL_TYPE *d;
  
  //
  REAL_TYPE scale = 1.0;
  int d_length = (size[0]+2*guide) * (size[1]+2*guide) * (size[2]+2*guide);
  
  // 出力用のヘッダ
  REAL_TYPE m_org[3], m_pit[3];
  
  //  ガイドセルがある場合(GuideOut != 0)にオリジナルポイントを調整
  for (int i=0; i<3; i++)
  {
    m_org[i] = origin[i] - pitch[i]*(REAL_TYPE)C.GuideOut;
    m_pit[i] = pitch[i];
  }
  
  // 出力ファイルの指定が有次元の場合
  if ( C.Unit.File == DIMENSIONAL )
  {
    for (int i=0; i<3; i++)
    {
      m_org[i] *= C.RefLength;
      m_pit[i] *= C.RefLength;
    }
  }
  
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
  
  // ガイドセル出力
  int gc_out = C.GuideOut;
  
  //allocate
  ngrid=C.P3Op.ngrid;
  id = new int[ngrid];
  jd = new int[ngrid];
  kd = new int[ngrid];
  nvar = new int[ngrid];
  
  //set grid data and nvar
  
  //for(igrid=0;igrid<ngrid;igrid++;){ //--->BCMでループが必要になる？
  
  igrid=0;
  id[igrid]=size[0]+2*guide;//将来的にGRIDはBCMのブロックになる？
  jd[igrid]=size[1]+2*guide;
  kd[igrid]=size[2]+2*guide;
  nvar[igrid]=C.P3Op.nvar;
  
  //}//igrid loop
  
  //set filename
  
  // 出力ファイル名
  std::string tmp;
  
  // 並列出力モード
  bool pout = ( C.FIO.IO_Output == IO_GATHER ) ? false : true;
  
  tmp = DFI.Generate_FileName_Free(C.P3Op.basename, "func", m_step, myRank, pout);
  if ( restart ) tmp = "restart_" + tmp; // リスタート用
  
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
  
  //for(igrid=0;igrid<ngrid;igrid++;){ //--->BCMでループが必要になる？
  
  igrid=0;
  
  d = new REAL_TYPE[ id[igrid]*jd[igrid]*kd[igrid]*nvar[igrid] ];
  
  //set grid data
  FP3DW.setGridData(id[igrid],jd[igrid],kd[igrid],ngrid);
  
  int ivar=0;
  
  // Pressure
  if (C.Unit.File == DIMENSIONAL)
  {
    REAL_TYPE bp = ( C.Unit.Prs == Unit_Absolute ) ? C.BasePrs : 0.0;
    fb_prs_nd2d_(d_ws, d_p, &d_length, &bp, &C.RefDensity, &C.RefVelocity, &scale, &flop);
  }
  else
  {
    fb_xcopy_(d_ws, d_p, &d_length, &scale, &flop);
  }
  
  //F.writeScalar(tmp, size, guide, d_ws, m_step, m_time, m_org, m_pit, gc_out);
  for(int k=0;k<kd[igrid];k++){
    for(int j=0;j<jd[igrid];j++){
      for(int i=0;i<id[igrid];i++){
        int ip=k*id[igrid]*jd[igrid]+j*id[igrid]+i;
        d[ip+d_length*ivar]=d_ws[ip];
      }
    }
  }
  ivar++;
  
  // Velocity
  REAL_TYPE unit_velocity = (C.Unit.File == DIMENSIONAL) ? C.RefVelocity : 1.0;
  fb_shift_refv_out_(d_wo, d_v, size, &guide, v00, &scale, &unit_velocity, &flop);
  
  //F.writeVector(tmp, size, guide, d_wo, m_step, m_time, m_org, m_pit, gc_out);
  for(int k=0;k<kd[igrid];k++){
    for(int j=0;j<jd[igrid];j++){
      for(int i=0;i<id[igrid];i++){
        int ip=k*id[igrid]*jd[igrid]+j*id[igrid]+i;
        d[ip+d_length*(ivar+0)]=d_wo[ip*3+0];
        d[ip+d_length*(ivar+1)]=d_wo[ip*3+1];
        d[ip+d_length*(ivar+2)]=d_wo[ip*3+2];
      }
    }
  }
  ivar=ivar+3;
  
  // Tempearture
  if( C.isHeatProblem() )
  {
    if (C.Unit.File == DIMENSIONAL)
    {
      REAL_TYPE klv = ( C.Unit.Temp == Unit_KELVIN ) ? 0.0 : KELVIN;
      fb_tmp_nd2d_(d_ws, d_t, &d_length, &C.BaseTemp, &C.DiffTemp, &klv, &scale, &flop);
    }
    else
    {
      fb_xcopy_(d_ws, d_t, &d_length, &scale, &flop);
    }
    
    for(int k=0;k<kd[igrid];k++){
      for(int j=0;j<jd[igrid];j++){
        for(int i=0;i<id[igrid];i++){
          int ip=k*id[igrid]*jd[igrid]+j*id[igrid]+i;
          d[ip+d_length*ivar]=d_ws[ip];
        }
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
    
    for(int k=0;k<kd[igrid];k++){
      for(int j=0;j<jd[igrid];j++){
        for(int i=0;i<id[igrid];i++){
          int ip=k*id[igrid]*jd[igrid]+j*id[igrid]+i;
          d[ip+d_length*ivar]=d_ws[ip];
        }
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
    
    for(int k=0;k<kd[igrid];k++){
      for(int j=0;j<jd[igrid];j++){
        for(int i=0;i<id[igrid];i++){
          int ip=k*id[igrid]*jd[igrid]+j*id[igrid]+i;
          d[ip+d_length*(ivar+0)]=d_wo[ip*3+0];
          d[ip+d_length*(ivar+1)]=d_wo[ip*3+1];
          d[ip+d_length*(ivar+2)]=d_wo[ip*3+2];
        }
      }
    }
    ivar=ivar+3;
  }
  
  //// 2nd Invariant of Velocity Gradient Tensor
  //if (C.Mode.I2VGT == ON ) {
  //}
  
  // Helicity
  if (C.Mode.Helicity == ON )
  {
    helicity_(d_p0, size, &guide, &deltaX, d_v, d_bcv, v00, &flop);
    
    // 無次元で出力
    d_length = (size[0]+2*guide) * (size[1]+2*guide) * (size[2]+2*guide);
    fb_xcopy_(d_ws, d_p0, &d_length, &scale, &flop);
    
    for(int k=0;k<kd[igrid];k++){
      for(int j=0;j<jd[igrid];j++){
        for(int i=0;i<id[igrid];i++){
          int ip=k*id[igrid]*jd[igrid]+j*id[igrid]+i;
          d[ip+d_length*ivar]=d_ws[ip];
        }
      }
    }
    ivar++;
  }
  
  //write all
  FP3DW.setFuncData(ivar,d);
  if(!FP3DW.WriteFuncData()) std::cout << "error WriteFuncData" << std::endl;
  
  delete [] d;
  
  //}//igrid loop
  
  //close file
  FP3DW.CloseFile();
  
  //deallocate
  delete [] id;
  delete [] jd;
  delete [] kd;
  
}

//
void FFV::OutputPlot3D_function_name()
{
  //HostRankでデータ集約する？
  //if(myRank != 0 ) return;
  
  //function_nameファイルはかならずformatted形式
  int keep_format=FP3DW.GetFormat();
  if(!FP3DW.setFormat(2)) std::cout << "error set" << std::endl;
  
  //set filename
  
  // 出力ファイル名
  std::string tmp;
  
  // 並列出力モード
  bool pout = ( C.FIO.IO_Output == IO_GATHER ) ? false : true;
  
  tmp = DFI.Generate_FileName_Free(C.P3Op.basename, "nam", 0, myRank, pout);
  
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
  if (C.Mode.VRT == ON ) FP3DW.WriteFunctionName("Vorticity");
  
  
  // 2nd Invariant of Velocity Gradient Tensor
  //if (C.Mode.I2VGT == ON ) FP3DW.WriteFunctionName("");
  
  // Helicity
  if (C.Mode.Helicity == ON ) FP3DW.WriteFunctionName("Helicity");
  
  //close file
  FP3DW.CloseFile();
  
  //reset option
  if(!FP3DW.setFormat(keep_format)) std::cout << "error set" << std::endl;
  
}


//
void FFV::OutputPlot3D_fvbnd()
{
  return;
  
  //HostRankでデータ集約する？
  //if(myRank != 0 ) return;
  
  //fvbndファイルはかならずformatted形式
  int keep_format=FP3DW.GetFormat();
  if(!FP3DW.setFormat(2)) std::cout << "error set" << std::endl;
  
  //set filename
  
  // 出力ファイル名
  std::string tmp;
  
  // 並列出力モード
  bool pout = ( C.FIO.IO_Output == IO_GATHER ) ? false : true;
  
  tmp = DFI.Generate_FileName_Free(C.P3Op.basename, "fvbnd", 0, myRank, pout);
  
  //open file
  FP3DW.setFileName(tmp.c_str());
  if(!FP3DW.OpenFile()){
    Hostonly_ printf("Error : error OpenFile\n");
    Exit(0);
  }
  
  //set data
  
  int nbname=5;
  int nb=3;
  string boundary_name[5];
  int type[3];
  int gridnum[3];
  int Imin[3];
  int Imax[3];
  int Jmin[3];
  int Jmax[3];
  int Kmin[3];
  int Kmax[3];
  string ResultFlag[3];
  int dir[3];
  
  boundary_name[0]="wall";
  boundary_name[1]="slide wall";
  boundary_name[2]="outflow";
  boundary_name[3]="inflow";
  boundary_name[4]="synmmetry";
  
  type[0]=1;
  gridnum[0]=1;
  Imin[0]=1;
  Imax[0]=1;
  Jmin[0]=1;
  Jmax[0]=11;
  Kmin[0]=12;
  Kmax[0]=441;
  ResultFlag[0]="F";
  dir[0]=1;
  
  type[1]=4;
  gridnum[1]=2;
  Imin[1]=1;
  Imax[1]=1;
  Jmin[1]=1;
  Jmax[1]=3;
  Kmin[1]=6;
  Kmax[1]=66;
  ResultFlag[1]="T";
  dir[1]=0;
  
  type[2]=5;
  gridnum[2]=3;
  Imin[2]=1;
  Imax[2]=1;
  Jmin[2]=1;
  Jmax[2]=44;
  Kmin[2]=53;
  Kmax[2]=3;
  ResultFlag[2]="T";
  dir[2]=0;
  
  //write fvbnd
  FP3DW.WriteFVBND(
                   nbname, nb, boundary_name, type, gridnum,
                   Imin, Imax, Jmin, Jmax, Kmin, Kmax, ResultFlag, dir);
  
  //close file
  FP3DW.CloseFile();
  
  //reset option
  if(!FP3DW.setFormat(keep_format)) std::cout << "error set" << std::endl;
  
}
