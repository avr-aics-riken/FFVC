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

  // ガイドセル出力
  int gc_out = C.GuideOut;

  //allocate
  ngrid=C.P3Op.ngrid;
  id = new int[ngrid];
  jd = new int[ngrid];
  kd = new int[ngrid];
  
  //set grid data
  //for(igrid=0;igrid<ngrid;igrid++;){ //--->BCMでループが必要になる？
  
  igrid=0;
  id[igrid]=size[0]+2*gc_out;
  jd[igrid]=size[1]+2*gc_out;
  kd[igrid]=size[2]+2*gc_out;
  
  int out_length = id[igrid]*jd[igrid]*kd[igrid];

  //}//igrid loop
  
  
  //set filename
  
  // 出力ファイル名
  std::string tmp,rtmp;

  // 並列出力モード
  bool pout = ( C.FIO.IO_Output == IO_GATHER ) ? false : true;

  rtmp = C.P3Op.basename;
  if ( restart ) rtmp = "restart_" + rtmp;
  tmp = DFI.Generate_FileName_Free(rtmp, "xyz", 0, myRank, pout);
  //tmp = DFI.Generate_FileName_Free(C.P3Op.basename, "xyz", 0, myRank, pout);
  //if ( restart ) tmp = "restart_" + tmp; // リスタート用
  
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
        x[ip]=origin[0]+pitch[0]*(REAL_TYPE)i-pitch[0]*(REAL_TYPE)gc_out;
        y[ip]=origin[1]+pitch[1]*(REAL_TYPE)j-pitch[1]*(REAL_TYPE)gc_out;
        z[ip]=origin[2]+pitch[2]*(REAL_TYPE)k-pitch[2]*(REAL_TYPE)gc_out;
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
  id[igrid]=size[0]+2*gc_out;
  jd[igrid]=size[1]+2*gc_out;
  kd[igrid]=size[2]+2*gc_out;
  nvar[igrid]=C.P3Op.nvar;
  
  int out_length = id[igrid]*jd[igrid]*kd[igrid];

  //}//igrid loop

  //set filename
  
  // 出力ファイル名
  std::string tmp,rtmp;
  
  // 並列出力モード
  bool pout = ( C.FIO.IO_Output == IO_GATHER ) ? false : true;
  
  rtmp = C.P3Op.basename;
  if ( restart ) rtmp = "restart_" + rtmp;
  tmp = DFI.Generate_FileName_Free(rtmp, "func", m_step, myRank, pout);
  //tmp = DFI.Generate_FileName_Free(C.P3Op.basename, "func", m_step, myRank, pout);
  //if ( restart ) tmp = "restart_" + tmp; // リスタート用
  
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
  
  if(gc_out==0){

    for(int k=0+guide;k<kd[igrid]+guide;k++){
      for(int j=0+guide;j<jd[igrid]+guide;j++){
        for(int i=0+guide;i<id[igrid]+guide;i++){
          int ip1=k*(id[igrid]+2*guide)*(jd[igrid]+2*guide)+j*(id[igrid]+2*guide)+i;
          int ip2=(k-guide)*id[igrid]*jd[igrid]+(j-guide)*id[igrid]+(i-guide);
          d[ip2+out_length*ivar]=d_ws[ip1];
        }
      }
    }
  }
  else{
    for(int k=0;k<kd[igrid];k++){
      for(int j=0;j<jd[igrid];j++){
        for(int i=0;i<id[igrid];i++){
          int ip=k*id[igrid]*jd[igrid]+j*id[igrid]+i;
          d[ip+out_length*ivar]=d_ws[ip];
        }
      }
    }
  }
  ivar++;
  
  // Velocity
  REAL_TYPE unit_velocity = (C.Unit.File == DIMENSIONAL) ? C.RefVelocity : 1.0;
  fb_shift_refv_out_(d_wo, d_v, size, &guide, v00, &scale, &unit_velocity, &flop);

  //F.writeVector(tmp, size, guide, d_wo, m_step, m_time, m_org, m_pit, gc_out);
  if(gc_out==0){
    for(int k=0+guide;k<kd[igrid]+guide;k++){
      for(int j=0+guide;j<jd[igrid]+guide;j++){
        for(int i=0+guide;i<id[igrid]+guide;i++){
          int ip1=k*(id[igrid]+2*guide)*(jd[igrid]+2*guide)+j*(id[igrid]+2*guide)+i;
          int ip2=(k-guide)*id[igrid]*jd[igrid]+(j-guide)*id[igrid]+(i-guide);
          d[ip2+out_length*(ivar+0)]=d_wo[ip1*3+0];
          d[ip2+out_length*(ivar+1)]=d_wo[ip1*3+1];
          d[ip2+out_length*(ivar+2)]=d_wo[ip1*3+2];
        }
      }
    }
  }
  else{
    for(int k=0;k<kd[igrid];k++){
      for(int j=0;j<jd[igrid];j++){
        for(int i=0;i<id[igrid];i++){
          int ip=k*id[igrid]*jd[igrid]+j*id[igrid]+i;
          d[ip+out_length*(ivar+0)]=d_wo[ip*3+0];
          d[ip+out_length*(ivar+1)]=d_wo[ip*3+1];
          d[ip+out_length*(ivar+2)]=d_wo[ip*3+2];
        }
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
    
    if(gc_out==0){
      for(int k=0+guide;k<kd[igrid]+guide;k++){
        for(int j=0+guide;j<jd[igrid]+guide;j++){
          for(int i=0+guide;i<id[igrid]+guide;i++){
            int ip1=k*(id[igrid]+2*guide)*(jd[igrid]+2*guide)+j*(id[igrid]+2*guide)+i;
            int ip2=(k-guide)*id[igrid]*jd[igrid]+(j-guide)*id[igrid]+(i-guide);
            d[ip2+out_length*ivar]=d_ws[ip1];
          }
        }
      }
    }
    else{
      for(int k=0;k<kd[igrid];k++){
        for(int j=0;j<jd[igrid];j++){
          for(int i=0;i<id[igrid];i++){
            int ip=k*id[igrid]*jd[igrid]+j*id[igrid]+i;
            d[ip+out_length*ivar]=d_ws[ip];
          }
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
    
    if(gc_out==0){
      for(int k=0+guide;k<kd[igrid]+guide;k++){
        for(int j=0+guide;j<jd[igrid]+guide;j++){
          for(int i=0+guide;i<id[igrid]+guide;i++){
            int ip1=k*(id[igrid]+2*guide)*(jd[igrid]+2*guide)+j*(id[igrid]+2*guide)+i;
            int ip2=(k-guide)*id[igrid]*jd[igrid]+(j-guide)*id[igrid]+(i-guide);
            d[ip2+out_length*ivar]=d_ws[ip1];
          }
        }
      }
    }
    else{
      for(int k=0;k<kd[igrid];k++){
        for(int j=0;j<jd[igrid];j++){
          for(int i=0;i<id[igrid];i++){
            int ip=k*id[igrid]*jd[igrid]+j*id[igrid]+i;
            d[ip+out_length*ivar]=d_ws[ip];
          }
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
    
    if(gc_out==0){
      for(int k=0+guide;k<kd[igrid]+guide;k++){
        for(int j=0+guide;j<jd[igrid]+guide;j++){
          for(int i=0+guide;i<id[igrid]+guide;i++){
            int ip1=k*(id[igrid]+2*guide)*(jd[igrid]+2*guide)+j*(id[igrid]+2*guide)+i;
            int ip2=(k-guide)*id[igrid]*jd[igrid]+(j-guide)*id[igrid]+(i-guide);
            d[ip2+out_length*(ivar+0)]=d_wo[ip1*3+0];
            d[ip2+out_length*(ivar+1)]=d_wo[ip1*3+1];
            d[ip2+out_length*(ivar+2)]=d_wo[ip1*3+2];
          }
        }
      }
    }
    else{
      for(int k=0;k<kd[igrid];k++){
        for(int j=0;j<jd[igrid];j++){
          for(int i=0;i<id[igrid];i++){
            int ip=k*id[igrid]*jd[igrid]+j*id[igrid]+i;
            d[ip+out_length*(ivar+0)]=d_wo[ip*3+0];
            d[ip+out_length*(ivar+1)]=d_wo[ip*3+1];
            d[ip+out_length*(ivar+2)]=d_wo[ip*3+2];
          }
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
    
    if(gc_out==0){
      for(int k=0+guide;k<kd[igrid]+guide;k++){
        for(int j=0+guide;j<jd[igrid]+guide;j++){
          for(int i=0+guide;i<id[igrid]+guide;i++){
            int ip1=k*(id[igrid]+2*guide)*(jd[igrid]+2*guide)+j*(id[igrid]+2*guide)+i;
            int ip2=(k-guide)*id[igrid]*jd[igrid]+(j-guide)*id[igrid]+(i-guide);
            d[ip2+out_length*ivar]=d_ws[ip1];
          }
        }
      }
    }
    else{
      for(int k=0;k<kd[igrid];k++){
        for(int j=0;j<jd[igrid];j++){
          for(int i=0;i<id[igrid];i++){
            int ip=k*id[igrid]*jd[igrid]+j*id[igrid]+i;
            d[ip+out_length*ivar]=d_ws[ip];
          }
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
  //HostRankでデータ集約する？
  //if(myRank != 0 ) return;

  //fvbndファイルはかならずformatted形式
  int keep_format=FP3DW.GetFormat();
  if(!FP3DW.setFormat(2)) std::cout << "error set" << std::endl;
  
  // 出力ファイル名
  std::string tmp;
  
  // 並列出力モード
  bool pout = ( C.FIO.IO_Output == IO_GATHER ) ? false : true;
  
  tmp = DFI.Generate_FileName_Free(C.P3Op.basename, "xyz", 0, myRank, pout);
  tmp = tmp + ".fvbnd";
  
  //open file
  FP3DW.setFileName(tmp.c_str());
  if(!FP3DW.OpenFile()){
    Hostonly_ printf("Error : error OpenFile\n");
    Exit(0);
  }

  //境界名の取得--->将来的に独立ルーチンにして名前だけ保持し続ける？*.fvbndが複数の場合に対処
  vector<string> bcname;
  bcname.clear();
  B.GetBoundaryNameforPLOT3D(bcname,cmp);

#if 0
  vector<string>::const_iterator it;
  for (it = bcname.begin(); it != bcname.end(); it++) {
    cout << "name = " << (*it).c_str() << endl;
  }
#endif

  //write boundary
  BC.WriteBoundaryPLOT3D(&FP3DW,bcname);

  //close file
  FP3DW.CloseFile();
  
  //reset option
  if(!FP3DW.setFormat(keep_format)) std::cout << "error set" << std::endl;

  return;
}
