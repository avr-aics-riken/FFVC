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
 * @file   comb_plot3d.C
 * @brief  COMB Class
 * @author kero
 */


#include "comb.h"


// #################################################################
// 
void COMB::output_plot3d()
{
  if( !DI ) return;

  char tmp[FB_FILE_PATH_LENGTH];
  string prefix,outfile,infile;
  int d_type;

  int m_sv_type, m_d_type;
  int m_rank, m_step, m_imax, m_jmax, m_kmax;
  long long m_dstep, m_dimax, m_djmax, m_dkmax;
  int sz[3];
  long long dsz[3];
  float m_org[3], m_pit[3];
  float m_time;
  double m_dorg[3], m_dpit[3];
  double m_dtime;
  int dim;

  int fp_in =8;
  FILE *fp;

  // work area for read sph
  int vsize,wksize;
  int maxsize=0;
  int maxdim=0;
  float *wk;
  double *dwk;

  // work area for write plot3d
  int outsize;
  int maxid=0;
  int maxjd=0;
  int maxkd=0;
  int maxnvar=3;
  int nvar=ndfi;
  float *d;
  double *dd;
  int id,jd,kd;
  int ngrid=1;

// set loop count
  int nstep=DI[0].step.size();
  int nnode=DI[0].NodeInfoSize;


#if defined (_STAGING_)
#elif defined (_NO_STAGING_)
//並列処理のためのインデックス作成
  int iproc=0;
  int ips=0;
  int ip;
  index.clear();
  for(int istep=0; istep< nstep; istep++ ) { // step loop
    for(int inode=0; inode< nnode; inode++ ) { // node loop
      ip=ips+inode;
      if(myRank==iproc) index.push_back(ip);
      iproc++;
      if(iproc==numProc) iproc=0;
    }
    ips=ips+nstep;
  }
  LOG_OUTV_ {
    fprintf(fplog,"\n");
    for(int j=0; j< index.size(); j++ ) {
      fprintf(fplog,"\tindex[%4d] = %d\n",j,index[j]);
    }
  }
  STD_OUTV_ {
    fprintf(fplog,"\n");
    for(int j=0; j< index.size(); j++ ) {
      printf("\tindex[%4d] = %d\n",j,index[j]);
    }
  }
#else
#endif


// copy dfi file
  string dfi_in = Generate_DFI_Name(DI[0].Prefix);
  string dfi_out = Generate_DFI_Name(P3Op.basename);
  dfi_out = dirname + dfi_out;
  Copy_DFIFile(dfi_in, dfi_out, P3Op.basename, dfi_mng[var_Plot3D]);
  if(P3Op.IS_DivideFunc == ON){ // 分割出力のとき
    for(int i=0;i<ndfi;i++){
      dfi_mng[var_Plot3D]=0;//本来はカウンターだが初期判定にのみ利用
      std::string dfipre;
      dfipre = DI[i].Prefix + P3Op.basename;
      cout << "dfipre = " << dfipre << endl;
      dfi_in = Generate_DFI_Name(DI[i].Prefix);
      dfi_out = Generate_DFI_Name(dfipre);
      dfi_out = dirname + dfi_out;
      Copy_DFIFile(dfi_in, dfi_out, dfipre, dfi_mng[var_Plot3D]);
    }
  }

// check data type
  m_step = DI[0].step[0];
  m_rank = DI[0].Node[0].RankID;
  prefix = DI[0].Prefix;
  infile = Generate_FileName(prefix, m_step, m_rank, true);
  ReadSphDataType (&m_sv_type, &m_d_type, fp_in, infile);
  d_type=m_d_type;
  if( d_type == SPH_FLOAT ){// sphファイルがfloatの時、PLOT3D出力は必ずfloatにする
    FP3DR.setRealType(d_type);
    FP3DW.setRealType(d_type);
  }

// set work area size
#if defined (_STAGING_)
#elif defined (_NO_STAGING_)
  ips=0;
#else
#endif
  for(int istep=0; istep< nstep; istep++ ) { // step loop
    m_step=DI[0].step[istep];

    for(int inode=0; inode< nnode; inode++ ) { // node loop
#if defined (_STAGING_)
#elif defined (_NO_STAGING_)
      //並列処理
      int iskip=0;
      ip=ips+inode;
      for(int ic=0; ic< index.size(); ic++ ) {
        if(ip==index[ic]) iskip=1;
      }
      if(!iskip) continue;
#else
#endif

      int ivar=0;
      for(int i=0;i<ndfi;i++){
        prefix=DI[i].Prefix;
        if     ( !strcasecmp(prefix.c_str(), "prs_" ) ) dim=1;
        else if( !strcasecmp(prefix.c_str(), "vel_" ) ) dim=3;
        else if( !strcasecmp(prefix.c_str(), "tmp_" ) ) dim=1;
        else if( !strcasecmp(prefix.c_str(), "tp_"  ) ) dim=1;
        else if( !strcasecmp(prefix.c_str(), "vrt_" ) ) dim=3;
        else if( !strcasecmp(prefix.c_str(), "i2vgt_" ) ) dim=1;
        else if( !strcasecmp(prefix.c_str(), "hlt_" ) ) dim=1;
	    else{
          m_step = DI[i].step[istep];
          m_rank = DI[i].Node[inode].RankID;
          infile = Generate_FileName(prefix, m_step, m_rank, true);
          ReadSphDataType (&m_sv_type, &m_d_type, fp_in, infile);
          dim=m_sv_type;
        }
        if(maxdim<dim) maxdim=dim;
        ivar=ivar+dim;
      }
      nvar=ivar;

      if(FP3DW.GetFormat() == C_BINARY) if(maxnvar<dim) maxnvar=dim;  // C_BINARYでの出力は項目ごとに書き出し
      else                              if(maxnvar<ivar) maxnvar=ivar;// Fortranによる出力では出力項目すべてを一度に書き出し
      if( P3Op.IS_DivideFunc == ON ) maxnvar=maxdim; //項目別出力onの時

      // read sph header
      prefix = DI[0].Prefix;
      m_rank = DI[0].Node[inode].RankID;
      infile = Generate_FileName(prefix, m_step, m_rank, true);
      if( d_type == SPH_FLOAT ){
        ReadSphHeader (&m_step, &m_sv_type, &m_d_type, &m_imax, &m_jmax, &m_kmax, &m_time, m_org, m_pit, fp_in, infile);
      }else{
        ReadSphHeader (&m_dstep, &m_sv_type, &m_d_type, &m_dimax, &m_djmax, &m_dkmax, &m_dtime, m_dorg, m_dpit, fp_in, infile);
        m_step=(int)m_dstep;
        m_imax=(int)m_dimax;
        m_jmax=(int)m_djmax;
        m_kmax=(int)m_dkmax;
      }

      // set size
      int guide = DI[0].GuideCell;
      sz[0]=m_imax;
      sz[1]=m_jmax;
      sz[2]=m_kmax;
      int id,jd,kd;
      id=sz[0]+1-2*guide;//PLOT3D出力はguideセルを書かないのでguideセル分減じておく
      jd=sz[1]+1-2*guide;
      kd=sz[2]+1-2*guide;
      wksize=sz[0]*sz[1]*sz[2];
      if(maxsize<wksize) maxsize=wksize;
      if(maxid<id) maxid=id;
      if(maxjd<jd) maxjd=jd;
      if(maxkd<kd) maxkd=kd;

    }
#if defined (_STAGING_)
#elif defined (_NO_STAGING_)
    ips=ips+nnode;
#else
#endif
  }

// write memory size
  double mc; 

  // sph読み込みwk用メモリ計算
  mc = (double)maxsize*(double)maxdim;
  if(mc>(double)INT_MAX){// 整数値あふれ出しチェック --- 参考 894*894*894*3=2143550952 INT_MAX 2147483647
    printf("\tsize error : wkmaxsize>INT_MAX\n");
    Exit(0);
  }
  double sphMemory;
  if( d_type == SPH_FLOAT ) sphMemory = mc * (double)sizeof(float);
  else                      sphMemory = mc * (double)sizeof(double);

  // PLOT3D出力用wksizeのmax予測
  mc = (double)maxid*(double)maxjd*(double)maxkd*(double)maxnvar;
  if(mc>(double)INT_MAX){// 整数値あふれ出しチェック --- 参考 894*894*894*3=2143550952 INT_MAX 2147483647
    printf("\tsize error : psize>INT_MAX\n");
    Exit(0);
  }
  double plot3dMemory;
  if( FP3DW.GetRealType() == OUTPUT_FLOAT ) plot3dMemory = mc * (double)sizeof(float);
  else                                      plot3dMemory = mc * (double)sizeof(double);

  // メモリチェック
  double TotalMemory=sphMemory+plot3dMemory;
  LOG_OUT_ MemoryRequirement(TotalMemory,sphMemory,plot3dMemory,fplog);
  STD_OUT_ MemoryRequirement(TotalMemory,sphMemory,plot3dMemory,stdout);

// allocate work area
  unsigned long wkmaxsize = maxsize*maxdim;
  unsigned long psize = maxid*maxjd*maxkd*maxnvar;
  if( d_type == SPH_FLOAT ){
    if ( !(wk = new float[wkmaxsize]) ){
      printf("\tallocate error : wk\n");
      Exit(0);
    }
  }else{
    if ( !(dwk = new double[wkmaxsize]) ){
      printf("\tallocate error : dwk\n");
      Exit(0);
    }
  }
  if( FP3DW.GetRealType() == OUTPUT_FLOAT ){
    if ( !(d = new float[psize]) ){
      printf("\tallocate error : wk\n");
      Exit(0);
    }
  }else{
    if ( !(dd = new double[psize]) ){
      printf("\tallocate error : wk\n");
      Exit(0);
    }
  }

// output xyz file
  LOG_OUT_ fprintf(fplog,"\t*** output xyz file ***\n");
  STD_OUT_ printf("\t*** output xyz file ***\n");

  // step loop
#if defined (_STAGING_)
#elif defined (_NO_STAGING_)
  ips=0;
#else
#endif
  for(int istep=0; istep< nstep; istep++ ) {
    m_step=DI[0].step[istep];
    LOG_OUT_ fprintf(fplog,"\tstep = %d\n", m_step);
    STD_OUT_ printf("\tstep = %d\n", m_step);

    // node loop
    for(int inode=0; inode< nnode; inode++ ) {
#if defined (_STAGING_)
#elif defined (_NO_STAGING_)
      //並列処理
      int iskip=0;
      ip=ips+inode;
      for(int ic=0; ic< index.size(); ic++ ) {
        if(ip==index[ic]) iskip=1;
      }
      if(!iskip) continue;
#else
#endif

      // read sph header
      prefix=DI[0].Prefix;
      m_rank = DI[0].Node[inode].RankID;
      infile = Generate_FileName(prefix, m_step, m_rank, true);
      if( d_type == SPH_FLOAT ){
        ReadSphHeader (&m_step, &m_sv_type, &m_d_type, &m_imax, &m_jmax, &m_kmax, &m_time, m_org, m_pit, fp_in, infile);
      }else{
        ReadSphHeader (&m_dstep, &m_sv_type, &m_d_type, &m_dimax, &m_djmax, &m_dkmax, &m_dtime, m_dorg, m_dpit, fp_in, infile);
        m_step=(int)m_dstep;
        m_imax=(int)m_dimax;
        m_jmax=(int)m_djmax;
        m_kmax=(int)m_dkmax;
      }

      // set size
      int guide = DI[0].GuideCell;
      sz[0]=m_imax;
      sz[1]=m_jmax;
      sz[2]=m_kmax;
      id=sz[0]+1-2*guide;//+2*gc_out
      jd=sz[1]+1-2*guide;//+2*gc_out
      kd=sz[2]+1-2*guide;//+2*gc_out
      outsize=id*jd*kd;
      if(outsize>psize){
        printf("\tsize error : outsize>psize\n");
        Exit(0);
      }

      // write plot 3d xyz file
      if(  FP3DW.IsMoveGrid()){
        if( d_type == SPH_FLOAT ){
          if( FP3DW.GetRealType() == OUTPUT_FLOAT ){ // float ---> float
            OutputPlot3D_xyz(m_step, m_rank, guide, m_org, m_pit, sz, &d[0], &d[outsize], &d[outsize*2]);
          }else{ // float ---> double
            OutputPlot3D_xyz(m_step, m_rank, guide, m_org, m_pit, sz, &dd[0], &dd[outsize], &dd[outsize*2]);
          }
        }else{
          if( FP3DW.GetRealType() == OUTPUT_FLOAT ){ // double ---> float
            OutputPlot3D_xyz(m_step, m_rank, guide, m_dorg, m_dpit, sz, &d[0], &d[outsize], &d[outsize*2]);
          }else{ // double ---> double
            OutputPlot3D_xyz(m_step, m_rank, guide, m_dorg, m_dpit, sz, &dd[0], &dd[outsize], &dd[outsize*2]);
          }
        }
      }else{
        if(istep==0){
          if( d_type == SPH_FLOAT ){
            if( FP3DW.GetRealType() == OUTPUT_FLOAT ){ // float ---> float
              OutputPlot3D_xyz(m_step, m_rank, guide, m_org, m_pit, sz, &d[0], &d[outsize], &d[outsize*2]);
            }else{ // float ---> double
              OutputPlot3D_xyz(m_step, m_rank, guide, m_org, m_pit, sz, &dd[0], &dd[outsize], &dd[outsize*2]);
            }
          }else{
            if( FP3DW.GetRealType() == OUTPUT_FLOAT ){ // double ---> float
              OutputPlot3D_xyz(m_step, m_rank, guide, m_dorg, m_dpit, sz, &d[0], &d[outsize], &d[outsize*2]);
            }else{ // double ---> double
              OutputPlot3D_xyz(m_step, m_rank, guide, m_dorg, m_dpit, sz, &dd[0], &dd[outsize], &dd[outsize*2]);
            }
          }
        }
      }

    } // node
#if defined (_STAGING_)
#elif defined (_NO_STAGING_)
    ips=ips+nnode;
#else
#endif
  } // step

// output function file
  LOG_OUT_ fprintf(fplog,"\t*** output func file ***\n");
  STD_OUT_ printf("\t*** output func file ***\n");

  // step loop
#if defined (_STAGING_)
#elif defined (_NO_STAGING_)
  ips=0;
#else
#endif
  for(int istep=0; istep< nstep; istep++ ) {
    m_step=DI[0].step[istep];
    LOG_OUT_ fprintf(fplog,"\tstep = %d\n", m_step);
    STD_OUT_ printf("\tstep = %d\n", m_step);

    // node loop
    for(int inode=0; inode< nnode; inode++ ) {
      m_rank = DI[0].Node[inode].RankID;
      LOG_OUT_ fprintf(fplog,"\t  rank = %d\n", m_rank);
      STD_OUT_ printf("\t  rank = %d\n", m_rank);

#if defined (_STAGING_)
#elif defined (_NO_STAGING_)
      //並列処理
      int iskip=0;
      ip=ips+inode;
      for(int ic=0; ic< index.size(); ic++ ) {
        if(ip==index[ic]) iskip=1;
      }
      if(!iskip) continue;
#else
#endif

      if(P3Op.IS_DivideFunc == OFF){ // 一括出力のとき、先にヘッダーを呼んでgridなどを書き出す

        // read sph header
        prefix=DI[0].Prefix;
        infile = Generate_FileName(prefix, m_step, m_rank, true);
        if( d_type == SPH_FLOAT ){
          ReadSphHeader (&m_step, &m_sv_type, &m_d_type, &m_imax, &m_jmax, &m_kmax, &m_time, m_org, m_pit, fp_in, infile);
        }else{
          ReadSphHeader (&m_dstep, &m_sv_type, &m_d_type, &m_dimax, &m_djmax, &m_dkmax, &m_dtime, m_dorg, m_dpit, fp_in, infile);
          m_step=(int)m_dstep;
          m_imax=(int)m_dimax;
          m_jmax=(int)m_djmax;
          m_kmax=(int)m_dkmax;
        }

        // set size
        int guide = DI[0].GuideCell;
        sz[0]=m_imax;
        sz[1]=m_jmax;
        sz[2]=m_kmax;
        id=sz[0]+1-2*guide;//+2*gc_out
        jd=sz[1]+1-2*guide;//+2*gc_out
        kd=sz[2]+1-2*guide;//+2*gc_out
        outsize=id*jd*kd;
        if(outsize>psize){
          printf("\tsize error : outsize>psize\n");
          Exit(0);
        }

        // plot3d func file open
        std::string ptmp;
        ptmp = Generate_FileName_Free(P3Op.basename, "func", m_step, m_rank, true);
        ptmp = dirname + ptmp;
        FP3DW.setFileName(ptmp.c_str());
        if(!FP3DW.OpenFile()){
          printf("\terror OpenFile : %s\n", ptmp.c_str());
          Exit(0);
        }

        // write block data
        FP3DW.WriteNgrid(ngrid);
        FP3DW.WriteFuncBlockData(&id,&jd,&kd,&nvar,ngrid);
      }

      // component loop --- dfi file ---> prs_, vel, ,,,
      int ivar=0;
      for(int i=0;i<ndfi;i++){
        prefix=DI[i].Prefix;
        LOG_OUTV_ fprintf(fplog,"\t    COMPONENT : %s\n", prefix.c_str());
        STD_OUTV_ printf("\t    COMPONENT : %s\n", prefix.c_str());

        // Scalar or Vector
        if     ( !strcasecmp(prefix.c_str(), "prs_" ) ) dim=1;
        else if( !strcasecmp(prefix.c_str(), "vel_" ) ) dim=3;
        else if( !strcasecmp(prefix.c_str(), "tmp_" ) ) dim=1;
        else if( !strcasecmp(prefix.c_str(), "tp_"  ) ) dim=1;
        else if( !strcasecmp(prefix.c_str(), "vrt_" ) ) dim=3;
        else if( !strcasecmp(prefix.c_str(), "i2vgt_" ) ) dim=1;
        else if( !strcasecmp(prefix.c_str(), "hlt_" ) ) dim=1;
	    else{
          m_step = DI[i].step[istep];
          m_rank = DI[i].Node[inode].RankID;
          infile = Generate_FileName(prefix, m_step, m_rank, true);
          ReadSphDataType (&m_sv_type, &m_d_type, fp_in, infile);
          dim=m_sv_type;
        }

        // read sph header
        infile = Generate_FileName(prefix, m_step, m_rank, true);
        if( d_type == SPH_FLOAT ){
          ReadSphHeader (&m_step, &m_sv_type, &m_d_type, &m_imax, &m_jmax, &m_kmax, &m_time, m_org, m_pit, fp_in, infile);
        }else{
          ReadSphHeader (&m_dstep, &m_sv_type, &m_d_type, &m_dimax, &m_djmax, &m_dkmax, &m_dtime, m_dorg, m_dpit, fp_in, infile);
          m_step=(int)m_dstep;
          m_time=(float)m_dtime;
          m_imax=(int)m_dimax;
          m_jmax=(int)m_djmax;
          m_kmax=(int)m_dkmax;
        }

        // set size
        int guide = DI[0].GuideCell;
        sz[0]=m_imax;
        sz[1]=m_jmax;
        sz[2]=m_kmax;
        id=sz[0]+1-2*guide;//+2*gc_out
        jd=sz[1]+1-2*guide;//+2*gc_out
        kd=sz[2]+1-2*guide;//+2*gc_out
        vsize=sz[0]*sz[1]*sz[2];
        wksize = vsize*dim;
        outsize=id*jd*kd;
        if(wksize>maxsize*maxdim){
          printf("\tsize error : wksize>maxsize*maxdim\n");
          Exit(0);
        }
        if(outsize>psize){
          printf("\tsize error : outsize>psize\n");
          Exit(0);
        }

        // read sph data
        if( d_type == SPH_FLOAT ) ReadSphData ( wk, wksize, sz, dim, fp_in, infile );
        else                      ReadSphData ( dwk, wksize, sz, dim, fp_in, infile );

        if(P3Op.IS_DivideFunc == ON){ // 分割出力のとき
       
          // plot3d func file open
          std::string ptmp;
          ptmp = Generate_FileName_Free(P3Op.basename, "func", m_step, m_rank, true);
          ptmp = DI[i].Prefix + ptmp;
          ptmp = dirname + ptmp;
          FP3DW.setFileName(ptmp.c_str());
          if(!FP3DW.OpenFile()){
            printf("\terror OpenFile : %s\n", ptmp.c_str());
            Exit(0);
          }
       
          // write block data
          ivar=0;
          nvar=dim;
          FP3DW.WriteNgrid(ngrid);
          FP3DW.WriteFuncBlockData(&id,&jd,&kd,&nvar,ngrid);
        }

        // write function file
        FP3DW.setGridData(id,jd,kd,ngrid);
        FP3DW.setFuncDataNum(nvar);
        if( dim==1 ){
          if(FP3DW.GetFormat() == C_BINARY){ // C_BINARYでの出力は項目ごとに書き出し
            if( d_type == SPH_FLOAT ){
              if( FP3DW.GetRealType() == OUTPUT_FLOAT ){ // float ---> float
                setScalarGridData(sz,guide,d,wk,id,jd,kd);
                FP3DW.setFuncData(d);
              }else{ // float ---> double
                setScalarGridData(sz,guide,dd,wk,id,jd,kd);
                FP3DW.setFuncData(dd);
              }
              if(!FP3DW.WriteFuncData()) printf("\terror WriteFuncData\n");
            }else{
              if( FP3DW.GetRealType() == OUTPUT_FLOAT ){ // double ---> float
                setScalarGridData(sz,guide,d,dwk,id,jd,kd);
                FP3DW.setFuncData(d);
              }else{ // double ---> double
                setScalarGridData(sz,guide,dd,dwk,id,jd,kd);
                FP3DW.setFuncData(dd);
              }
              if(!FP3DW.WriteFuncData()) printf("\terror WriteFuncData\n");
            }
          }
          else // Fortranによる出力では出力項目すべてを一度に書き出し
          {
            if( d_type == SPH_FLOAT ){
              if( FP3DW.GetRealType() == OUTPUT_FLOAT ){ // float ---> float
                setScalarGridData(sz,guide,&d[outsize*ivar],wk,id,jd,kd);
              }else{ // float ---> double
                setScalarGridData(sz,guide,&dd[outsize*ivar],wk,id,jd,kd);
              }
            }else{
              if( FP3DW.GetRealType() == OUTPUT_FLOAT ){ // double ---> float
                setScalarGridData(sz,guide,&d[outsize*ivar],dwk,id,jd,kd);
              }else{ // double ---> double
                setScalarGridData(sz,guide,&dd[outsize*ivar],dwk,id,jd,kd);
              }
            }
          }
        }else if( dim==3 ){
          if(FP3DW.GetFormat() == C_BINARY){ // C_BINARYでの出力は項目ごとに書き出し
            if( d_type == SPH_FLOAT ){
              if( FP3DW.GetRealType() == OUTPUT_FLOAT ){ // float ---> float
                for(int ixyz=0;ixyz<3;ixyz++){
                  setVectorComponentGridData(sz,guide,d,wk,id,jd,kd,ixyz);
                  FP3DW.setFuncData(d);
                  if(!FP3DW.WriteFuncData()) printf("\terror WriteFuncData\n");
                }
              }else{ // float ---> double
                for(int ixyz=0;ixyz<3;ixyz++){
                  setVectorComponentGridData(sz,guide,dd,wk,id,jd,kd,ixyz);
                  FP3DW.setFuncData(dd);
                  if(!FP3DW.WriteFuncData()) printf("\terror WriteFuncData\n");
                }
              }
            }else{ // double ---> float
              if( FP3DW.GetRealType() == OUTPUT_FLOAT ){
                for(int ixyz=0;ixyz<3;ixyz++){
                  setVectorComponentGridData(sz,guide,d,dwk,id,jd,kd,ixyz);
                  FP3DW.setFuncData(d);
                  if(!FP3DW.WriteFuncData()) printf("\terror WriteFuncData\n");
                }
              }else{ // double ---> double
                for(int ixyz=0;ixyz<3;ixyz++){
                  setVectorComponentGridData(sz,guide,dd,dwk,id,jd,kd,ixyz);
                  FP3DW.setFuncData(dd);
                  if(!FP3DW.WriteFuncData()) printf("\terror WriteFuncData\n");
                }
              }
            }
          }
          else // Fortranによる出力では出力項目すべてを一度に書き出し
          {
            if( d_type == SPH_FLOAT ){
              if( FP3DW.GetRealType() == OUTPUT_FLOAT ){ // float ---> float
                setVectorGridData(sz,guide,&d[outsize*ivar],wk,id,jd,kd);
              }else{ // float ---> double
                setVectorGridData(sz,guide,&dd[outsize*ivar],wk,id,jd,kd);
              }
            }else{
              if( FP3DW.GetRealType() == OUTPUT_FLOAT ){ // double ---> float
                setVectorGridData(sz,guide,&d[outsize*ivar],dwk,id,jd,kd);
              }else{ // double ---> double
                setVectorGridData(sz,guide,&dd[outsize*ivar],dwk,id,jd,kd);
              }
            }
          }
        }
        
        if(P3Op.IS_DivideFunc == ON){ // 分割出力のとき
          // write all
          if(FP3DW.GetFormat() != C_BINARY){ // C_BINARY以外は出力項目すべてを一度に書き出し
            if( FP3DW.GetRealType() == OUTPUT_FLOAT ) FP3DW.setFuncData(d);
            else                                      FP3DW.setFuncData(dd);
            if(!FP3DW.WriteFuncData()) printf("\terror WriteFuncData\n");
          }
          // plot3d func file close
          FP3DW.CloseFile();

          //dfiファイルの出力
          std::string dfipre;
          dfipre = DI[i].Prefix + P3Op.basename;
          dfipre = dirname + dfipre;
          if(inode == 0) if ( !DFI.Write_DFI_File(dfipre, m_step, (double)m_time, dfi_mng[var_Plot3D], true) ) Exit(0);
        }

        ivar=ivar+dim;
      }

      if(P3Op.IS_DivideFunc == OFF){ // 一括出力のとき
        // write all
        if(FP3DW.GetFormat() != C_BINARY){//C_BINARY以外は出力項目すべてを一度に書き出し
          if( FP3DW.GetRealType() == OUTPUT_FLOAT ) FP3DW.setFuncData(d);
          else                                      FP3DW.setFuncData(dd);
          if(!FP3DW.WriteFuncData()) printf("\terror WriteFuncData\n");
        }
        // plot3d func file close
        FP3DW.CloseFile();
      }

      //dfiファイルの出力
      std::string dfipre;
      dfipre = dirname + P3Op.basename;
      if(inode == 0) if ( !DFI.Write_DFI_File(dfipre, m_step, (double)m_time, dfi_mng[var_Plot3D], true) ) Exit(0);

    }//node loop
#if defined (_STAGING_)
#elif defined (_NO_STAGING_)
    ips=ips+nnode;
#else
#endif
  }//step loop


// deallocate work area
  if( d_type == SPH_FLOAT ) delete[] wk;
  else                      delete[] dwk;
  if( FP3DW.GetRealType() == OUTPUT_FLOAT ) delete[] d;
  else                                      delete[] dd;

// output function name file
  LOG_OUT_ fprintf(fplog,"\t*** output func name file ***\n");
  STD_OUT_ printf("\t*** output func name file ***\n");

  // function_nameファイルはかならずformatted形式
  int keep_format=FP3DW.GetFormat();
  FP3DW.setFormat(FORMATTED);
  
  // open file
  if(P3Op.IS_DivideFunc == OFF){ // 一括出力のとき
    std::string ptmp;
    ptmp = DFI.Generate_FileName_Free(P3Op.basename, "nam", 0, 0, false);
    ptmp = dirname + ptmp;
    FP3DW.setFileName(ptmp.c_str());
    if(!FP3DW.OpenFile()){
      printf("Error : error OpenFile\n");
      Exit(0);
    }
  }

  int unknowncomp=0;
  for(int i=0;i<ndfi;i++){
    string comp = DI[i].Prefix;
    LOG_OUT_ fprintf(fplog,"\t    COMPONENT : %s\n", comp.c_str());
    STD_OUT_ printf("\t    COMPONENT : %s\n", comp.c_str());

    // open file
    if(P3Op.IS_DivideFunc == ON){ // 分割出力のとき
      std::string ptmp;
      ptmp = Generate_FileName_Free(P3Op.basename, "nam", 0, 0, true);
      ptmp = DI[i].Prefix + ptmp;
      ptmp = dirname + ptmp;
      FP3DW.setFileName(ptmp.c_str());
      if(!FP3DW.OpenFile()){
        printf("\terror OpenFile : %s\n", ptmp.c_str());
        Exit(0);
      }
    }

    // write comp name
    if     ( !strcasecmp(comp.c_str(), "prs_" ) ){
      FP3DW.WriteFunctionName("Pressure");
    }else if( !strcasecmp(comp.c_str(), "vel_" ) ){
      FP3DW.WriteFunctionName("U-Velocity ; Velocity");
      FP3DW.WriteFunctionName("V-Velocity");
      FP3DW.WriteFunctionName("W-Velocity");
    }else if( !strcasecmp(comp.c_str(), "tmp_" ) ){
      FP3DW.WriteFunctionName("Tempearture");
    }else if( !strcasecmp(comp.c_str(), "tp_"  ) ){
      FP3DW.WriteFunctionName("Total_Pressure");
    }else if( !strcasecmp(comp.c_str(), "vrt_" ) ){
      FP3DW.WriteFunctionName("U-Vorticity ; Vorticity");
      FP3DW.WriteFunctionName("V-Vorticity");
      FP3DW.WriteFunctionName("W-Vorticity");
    }else if( !strcasecmp(comp.c_str(), "i2vgt_"  ) ){
      FP3DW.WriteFunctionName("2nd Invariant of Velocity Gradient Tensor");
    }else if( !strcasecmp(comp.c_str(), "hlt_" ) ){
      FP3DW.WriteFunctionName("Helicity");
    }else{
#if 0
      unknowncomp++;
      char* numbuff = new char[DFI_LINE_LENGTH];
      memset(numbuff, 0, sizeof(char)*DFI_LINE_LENGTH);
      sprintf(numbuff, "%02d", unknowncomp);
      string compbuff=(numbuff);
      compbuff ="F-"+compbuff;
      FP3DW.WriteFunctionName(compbuff.c_str());
      if ( numbuff ) delete [] numbuff;
#else
      FP3DW.WriteFunctionName(comp.c_str());
#endif
    }

    // close file
    if(P3Op.IS_DivideFunc == ON) FP3DW.CloseFile();
  }

  // close file
  if(P3Op.IS_DivideFunc == OFF) FP3DW.CloseFile();
  
  // reset option
  FP3DW.setFormat(keep_format);
  
}


//
void COMB::OutputPlot3D_xyz(int m_step, int m_rank, int guide, float* origin, float* pitch, int* size, float* x, float* y, float* z)
{
  //value
  int ngrid=1;
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  int gc_out = 0;//plot3dは常にガイドセルを出力しない

  int *iblank=NULL;//dummy
  int id,jd,kd;//出力サイズ
  id=size[0]+1;//+2*gc_out
  jd=size[1]+1;//+2*gc_out
  kd=size[2]+1;//+2*gc_out

  // ガイドセル出力があった場合オリジナルポイントを調整しておく
  float m_org[3], m_pit[3];
  for (int i=0; i<3; i++) 
  {
    m_org[i] = (float)origin[i] + (float)pitch[i]*(float)gd;
    m_pit[i] = (float)pitch[i];
  }

  // 出力ファイル名
  std::string tmp;
  tmp = Generate_FileName_Free(P3Op.basename, "xyz", m_step, m_rank, true);
  tmp = dirname + tmp;

  //open file
  FP3DW.setFileName(tmp.c_str());
  if(!FP3DW.OpenFile()){
    printf("Error : error OpenFile\n");
    Exit(0);
  }

  //write block data
  FP3DW.WriteNgrid(ngrid);//if multi grid
  FP3DW.WriteBlockData(id,jd,kd);

  for(int k=0;k<kd;k++){
    for(int j=0;j<jd;j++){
      for(int i=0;i<id;i++){
        int ip=k*id*jd+j*id+i;
        //int ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
        x[ip]=m_org[0]+m_pit[0]*(float)i;//-pitch[0]*(float)gc_out;
        y[ip]=m_org[1]+m_pit[1]*(float)j;//-pitch[1]*(float)gc_out;
        z[ip]=m_org[2]+m_pit[2]*(float)k;//-pitch[2]*(float)gc_out;
      }
    }
  }

  //write
  FP3DW.setGridData(id,jd,kd,ngrid);
  FP3DW.setXYZData(x,y,z,iblank);
  if(!FP3DW.WriteXYZData()) printf("\terror WriteXYZData\n");
  
  //close file
  FP3DW.CloseFile();
}


//
void COMB::OutputPlot3D_xyz(int m_step, int m_rank, int guide, double* origin, double* pitch, int* size, double* x, double* y, double* z)
{
  //value
  int ngrid=1;
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  int gc_out = 0;//plot3dは常にガイドセルを出力しない

  //int *id,*jd,*kd;//出力サイズ
  //double *x,*y,*z;
  int *iblank=NULL;//dummy
  int id,jd,kd;//出力サイズ
  id=size[0]+1;//+2*gc_out
  jd=size[1]+1;//+2*gc_out
  kd=size[2]+1;//+2*gc_out

  // ガイドセル出力があった場合オリジナルポイントを調整しておく
  double m_org[3], m_pit[3];
  for (int i=0; i<3; i++) 
  {
    m_org[i] = (double)origin[i] + (double)pitch[i]*(double)gd;
    m_pit[i] = (double)pitch[i];
  }

  // 出力ファイル名
  std::string tmp;
  tmp = Generate_FileName_Free(P3Op.basename, "xyz", m_step, m_rank, true);
  tmp = dirname + tmp;

  //open file
  FP3DW.setFileName(tmp.c_str());
  if(!FP3DW.OpenFile()){
    printf("Error : error OpenFile\n");
    Exit(0);
  }

  //write block data
  FP3DW.WriteNgrid(ngrid);//if multi grid
  FP3DW.WriteBlockData(id,jd,kd);

  for(int k=0;k<kd;k++){
    for(int j=0;j<jd;j++){
      for(int i=0;i<id;i++){
        int ip=k*id*jd+j*id+i;
        //int ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
        x[ip]=m_org[0]+m_pit[0]*(double)i;//-pitch[0]*(double)gc_out;
        y[ip]=m_org[1]+m_pit[1]*(double)j;//-pitch[1]*(double)gc_out;
        z[ip]=m_org[2]+m_pit[2]*(double)k;//-pitch[2]*(double)gc_out;
      }
    }
  }

  //write
  FP3DW.setGridData(id,jd,kd,ngrid);
  FP3DW.setXYZData(x,y,z,iblank);
  if(!FP3DW.WriteXYZData()) printf("\terror WriteXYZData\n");
  
  //close file
  FP3DW.CloseFile();
}




//
void COMB::OutputPlot3D_xyz(int m_step, int m_rank, int guide, float* origin, float* pitch, int* size, double* x, double* y, double* z)
{
  //value
  int ngrid=1;
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  int gc_out = 0;//plot3dは常にガイドセルを出力しない

  int *iblank=NULL;//dummy
  int id,jd,kd;//出力サイズ
  id=size[0]+1;//+2*gc_out
  jd=size[1]+1;//+2*gc_out
  kd=size[2]+1;//+2*gc_out

  // ガイドセル出力があった場合オリジナルポイントを調整しておく
  double m_org[3], m_pit[3];
  for (int i=0; i<3; i++) 
  {
    m_org[i] = (double)origin[i] + (double)pitch[i]*(double)gd;
    m_pit[i] = (double)pitch[i];
  }

  // 出力ファイル名
  std::string tmp;
  tmp = Generate_FileName_Free(P3Op.basename, "xyz", m_step, m_rank, true);
  tmp = dirname + tmp;

  //open file
  FP3DW.setFileName(tmp.c_str());
  if(!FP3DW.OpenFile()){
    printf("Error : error OpenFile\n");
    Exit(0);
  }

  //write block data
  FP3DW.WriteNgrid(ngrid);//if multi grid
  FP3DW.WriteBlockData(id,jd,kd);

  for(int k=0;k<kd;k++){
    for(int j=0;j<jd;j++){
      for(int i=0;i<id;i++){
        int ip=k*id*jd+j*id+i;
        //int ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
        x[ip]=m_org[0]+m_pit[0]*(double)i;//-pitch[0]*(double)gc_out;
        y[ip]=m_org[1]+m_pit[1]*(double)j;//-pitch[1]*(double)gc_out;
        z[ip]=m_org[2]+m_pit[2]*(double)k;//-pitch[2]*(double)gc_out;
      }
    }
  }

  //write
  FP3DW.setGridData(id,jd,kd,ngrid);
  FP3DW.setXYZData(x,y,z,iblank);
  if(!FP3DW.WriteXYZData()) printf("\terror WriteXYZData\n");
  
  //close file
  FP3DW.CloseFile();
}



//
void COMB::OutputPlot3D_xyz(int m_step, int m_rank, int guide, double* origin, double* pitch, int* size, float* x, float* y, float* z)
{
  //value
  int ngrid=1;
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  int gc_out = 0;//plot3dは常にガイドセルを出力しない

  int *iblank=NULL;//dummy
  int id,jd,kd;//出力サイズ
  id=size[0]+1;//+2*gc_out
  jd=size[1]+1;//+2*gc_out
  kd=size[2]+1;//+2*gc_out

  // ガイドセル出力があった場合オリジナルポイントを調整しておく
  float m_org[3], m_pit[3];
  for (int i=0; i<3; i++) 
  {
    m_org[i] = (float)origin[i] + (float)pitch[i]*(float)gd;
    m_pit[i] = (float)pitch[i];
  }

  // 出力ファイル名
  std::string tmp;
  tmp = Generate_FileName_Free(P3Op.basename, "xyz", m_step, m_rank, true);
  tmp = dirname + tmp;

  //open file
  FP3DW.setFileName(tmp.c_str());
  if(!FP3DW.OpenFile()){
    printf("Error : error OpenFile\n");
    Exit(0);
  }

  //write block data
  FP3DW.WriteNgrid(ngrid);//if multi grid
  FP3DW.WriteBlockData(id,jd,kd);

  for(int k=0;k<kd;k++){
    for(int j=0;j<jd;j++){
      for(int i=0;i<id;i++){
        int ip=k*id*jd+j*id+i;
        //int ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
        x[ip]=m_org[0]+m_pit[0]*(float)i;//-pitch[0]*(float)gc_out;
        y[ip]=m_org[1]+m_pit[1]*(float)j;//-pitch[1]*(float)gc_out;
        z[ip]=m_org[2]+m_pit[2]*(float)k;//-pitch[2]*(float)gc_out;
      }
    }
  }

  //write
  FP3DW.setGridData(id,jd,kd,ngrid);
  FP3DW.setXYZData(x,y,z,iblank);
  if(!FP3DW.WriteXYZData()) printf("\terror WriteXYZData\n");
  
  //close file
  FP3DW.CloseFile();
}

//*****************************************************************************
//*****************************************************************************
//*****************************************************************************

// Scalarの格子点での値をセット（ガイドセルの値は計算対象にいれていない）
void COMB::setScalarGridData(int* size, int guide, float* d, float* data, int id, int jd, int kd)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;

  for (int i=0; i<id*jd*kd; i++) d[i]=0.0;

  for (int km=1; km<=kx; km++) {
    for (int jm=1; jm<=jx; jm++) {
      for (int im=1; im<=ix; im++) {
        int mip = _F_IDX_S3D(im, jm, km, ix, jx, kx, gd);
        float ddd=data[mip];
        int i=im-1;
        int j=jm-1;
        int k=km-1;
        int ip1= k   *id*jd+ j   *id+i;
        int ip2= k   *id*jd+ j   *id+i+1;
        int ip3= k   *id*jd+(j+1)*id+i+1;
        int ip4= k   *id*jd+(j+1)*id+i;
        int ip5=(k+1)*id*jd+ j   *id+i;
        int ip6=(k+1)*id*jd+ j   *id+i+1;
        int ip7=(k+1)*id*jd+(j+1)*id+i+1;
        int ip8=(k+1)*id*jd+(j+1)*id+i;
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
void COMB::setVectorGridData(int* size, int guide, float* d, float* data, int id, int jd, int kd)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;

  int dsize=id*jd*kd;

  for (int i=0; i<id*jd*kd*3; i++) d[i]=0.0;

  for (int ivar=0;ivar<3;ivar++){

    for (int km=1; km<=kx; km++) {
      for (int jm=1; jm<=jx; jm++) {
        for (int im=1; im<=ix; im++) {
          int mip = _F_IDX_V3DEX(ivar, im, jm, km, ix, jx, kx, gd);//ivar=0:x方向、ivar=1:y方向、ivar=2:z方向
          float ddd=data[mip];
          int i=im-1;
          int j=jm-1;
          int k=km-1;
          int ip1= k   *id*jd+ j   *id+i;
          int ip2= k   *id*jd+ j   *id+i+1;
          int ip3= k   *id*jd+(j+1)*id+i+1;
          int ip4= k   *id*jd+(j+1)*id+i;
          int ip5=(k+1)*id*jd+ j   *id+i;
          int ip6=(k+1)*id*jd+ j   *id+i+1;
          int ip7=(k+1)*id*jd+(j+1)*id+i+1;
          int ip8=(k+1)*id*jd+(j+1)*id+i;
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
void COMB::setVectorComponentGridData(int* size, int guide, float* d, float* data, int id, int jd, int kd, int ivar)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;

  for (int i=0; i<id*jd*kd; i++) d[i]=0.0;
 
  for (int km=1; km<=kx; km++) {
    for (int jm=1; jm<=jx; jm++) {
      for (int im=1; im<=ix; im++) {
        int mip = _F_IDX_V3DEX(ivar, im, jm, km, ix, jx, kx, gd);//ivar=0:x方向、ivar=1:y方向、ivar=2:z方向
        float ddd=data[mip];
        int i=im-1;
        int j=jm-1;
        int k=km-1;
        int ip1= k   *id*jd+ j   *id+i;
        int ip2= k   *id*jd+ j   *id+i+1;
        int ip3= k   *id*jd+(j+1)*id+i+1;
        int ip4= k   *id*jd+(j+1)*id+i;
        int ip5=(k+1)*id*jd+ j   *id+i;
        int ip6=(k+1)*id*jd+ j   *id+i+1;
        int ip7=(k+1)*id*jd+(j+1)*id+i+1;
        int ip8=(k+1)*id*jd+(j+1)*id+i;
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
void COMB::setScalarGridData(int* size, int guide, double* d, double* data, int id, int jd, int kd)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;

  for (int i=0; i<id*jd*kd; i++) d[i]=0.0;

  for (int km=1; km<=kx; km++) {
    for (int jm=1; jm<=jx; jm++) {
      for (int im=1; im<=ix; im++) {
        int mip = _F_IDX_S3D(im, jm, km, ix, jx, kx, gd);
        double ddd=data[mip];
        int i=im-1;
        int j=jm-1;
        int k=km-1;
        int ip1= k   *id*jd+ j   *id+i;
        int ip2= k   *id*jd+ j   *id+i+1;
        int ip3= k   *id*jd+(j+1)*id+i+1;
        int ip4= k   *id*jd+(j+1)*id+i;
        int ip5=(k+1)*id*jd+ j   *id+i;
        int ip6=(k+1)*id*jd+ j   *id+i+1;
        int ip7=(k+1)*id*jd+(j+1)*id+i+1;
        int ip8=(k+1)*id*jd+(j+1)*id+i;
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
void COMB::setVectorGridData(int* size, int guide, double* d, double* data, int id, int jd, int kd)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;

  int dsize=id*jd*kd;

  for (int i=0; i<id*jd*kd*3; i++) d[i]=0.0;

  for (int ivar=0;ivar<3;ivar++){

    for (int km=1; km<=kx; km++) {
      for (int jm=1; jm<=jx; jm++) {
        for (int im=1; im<=ix; im++) {
          int mip = _F_IDX_V3DEX(ivar, im, jm, km, ix, jx, kx, gd);//ivar=0:x方向、ivar=1:y方向、ivar=2:z方向
          double ddd=data[mip];
          int i=im-1;
          int j=jm-1;
          int k=km-1;
          int ip1= k   *id*jd+ j   *id+i;
          int ip2= k   *id*jd+ j   *id+i+1;
          int ip3= k   *id*jd+(j+1)*id+i+1;
          int ip4= k   *id*jd+(j+1)*id+i;
          int ip5=(k+1)*id*jd+ j   *id+i;
          int ip6=(k+1)*id*jd+ j   *id+i+1;
          int ip7=(k+1)*id*jd+(j+1)*id+i+1;
          int ip8=(k+1)*id*jd+(j+1)*id+i;
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
void COMB::setVectorComponentGridData(int* size, int guide, double* d, double* data, int id, int jd, int kd, int ivar)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;

  for (int i=0; i<id*jd*kd; i++) d[i]=0.0;
 
  for (int km=1; km<=kx; km++) {
    for (int jm=1; jm<=jx; jm++) {
      for (int im=1; im<=ix; im++) {
        int mip = _F_IDX_V3DEX(ivar, im, jm, km, ix, jx, kx, gd);//ivar=0:x方向、ivar=1:y方向、ivar=2:z方向
        double ddd=data[mip];
        int i=im-1;
        int j=jm-1;
        int k=km-1;
        int ip1= k   *id*jd+ j   *id+i;
        int ip2= k   *id*jd+ j   *id+i+1;
        int ip3= k   *id*jd+(j+1)*id+i+1;
        int ip4= k   *id*jd+(j+1)*id+i;
        int ip5=(k+1)*id*jd+ j   *id+i;
        int ip6=(k+1)*id*jd+ j   *id+i+1;
        int ip7=(k+1)*id*jd+(j+1)*id+i+1;
        int ip8=(k+1)*id*jd+(j+1)*id+i;
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
void COMB::setScalarGridData(int* size, int guide, float* d, double* data, int id, int jd, int kd)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;

  for (int i=0; i<id*jd*kd; i++) d[i]=0.0;

  for (int km=1; km<=kx; km++) {
    for (int jm=1; jm<=jx; jm++) {
      for (int im=1; im<=ix; im++) {
        int mip = _F_IDX_S3D(im, jm, km, ix, jx, kx, gd);
        float ddd=(float)data[mip];
        int i=im-1;
        int j=jm-1;
        int k=km-1;
        int ip1= k   *id*jd+ j   *id+i;
        int ip2= k   *id*jd+ j   *id+i+1;
        int ip3= k   *id*jd+(j+1)*id+i+1;
        int ip4= k   *id*jd+(j+1)*id+i;
        int ip5=(k+1)*id*jd+ j   *id+i;
        int ip6=(k+1)*id*jd+ j   *id+i+1;
        int ip7=(k+1)*id*jd+(j+1)*id+i+1;
        int ip8=(k+1)*id*jd+(j+1)*id+i;
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
void COMB::setVectorGridData(int* size, int guide, float* d, double* data, int id, int jd, int kd)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;

  int dsize=id*jd*kd;

  for (int i=0; i<id*jd*kd*3; i++) d[i]=0.0;

  for (int ivar=0;ivar<3;ivar++){

    for (int km=1; km<=kx; km++) {
      for (int jm=1; jm<=jx; jm++) {
        for (int im=1; im<=ix; im++) {
          int mip = _F_IDX_V3DEX(ivar, im, jm, km, ix, jx, kx, gd);//ivar=0:x方向、ivar=1:y方向、ivar=2:z方向
          float ddd=(float)data[mip];
          int i=im-1;
          int j=jm-1;
          int k=km-1;
          int ip1= k   *id*jd+ j   *id+i;
          int ip2= k   *id*jd+ j   *id+i+1;
          int ip3= k   *id*jd+(j+1)*id+i+1;
          int ip4= k   *id*jd+(j+1)*id+i;
          int ip5=(k+1)*id*jd+ j   *id+i;
          int ip6=(k+1)*id*jd+ j   *id+i+1;
          int ip7=(k+1)*id*jd+(j+1)*id+i+1;
          int ip8=(k+1)*id*jd+(j+1)*id+i;
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
void COMB::setVectorComponentGridData(int* size, int guide, float* d, double* data, int id, int jd, int kd, int ivar)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;

  for (int i=0; i<id*jd*kd; i++) d[i]=0.0;
 
  for (int km=1; km<=kx; km++) {
    for (int jm=1; jm<=jx; jm++) {
      for (int im=1; im<=ix; im++) {
        int mip = _F_IDX_V3DEX(ivar, im, jm, km, ix, jx, kx, gd);//ivar=0:x方向、ivar=1:y方向、ivar=2:z方向
        float ddd=(float)data[mip];
        int i=im-1;
        int j=jm-1;
        int k=km-1;
        int ip1= k   *id*jd+ j   *id+i;
        int ip2= k   *id*jd+ j   *id+i+1;
        int ip3= k   *id*jd+(j+1)*id+i+1;
        int ip4= k   *id*jd+(j+1)*id+i;
        int ip5=(k+1)*id*jd+ j   *id+i;
        int ip6=(k+1)*id*jd+ j   *id+i+1;
        int ip7=(k+1)*id*jd+(j+1)*id+i+1;
        int ip8=(k+1)*id*jd+(j+1)*id+i;
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
void COMB::setScalarGridData(int* size, int guide, double* d, float* data, int id, int jd, int kd)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;

  for (int i=0; i<id*jd*kd; i++) d[i]=0.0;

  for (int km=1; km<=kx; km++) {
    for (int jm=1; jm<=jx; jm++) {
      for (int im=1; im<=ix; im++) {
        int mip = _F_IDX_S3D(im, jm, km, ix, jx, kx, gd);
        double ddd=(double)data[mip];
        int i=im-1;
        int j=jm-1;
        int k=km-1;
        int ip1= k   *id*jd+ j   *id+i;
        int ip2= k   *id*jd+ j   *id+i+1;
        int ip3= k   *id*jd+(j+1)*id+i+1;
        int ip4= k   *id*jd+(j+1)*id+i;
        int ip5=(k+1)*id*jd+ j   *id+i;
        int ip6=(k+1)*id*jd+ j   *id+i+1;
        int ip7=(k+1)*id*jd+(j+1)*id+i+1;
        int ip8=(k+1)*id*jd+(j+1)*id+i;
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
void COMB::setVectorGridData(int* size, int guide, double* d, float* data, int id, int jd, int kd)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;

  int dsize=id*jd*kd;

  for (int i=0; i<id*jd*kd*3; i++) d[i]=0.0;

  for (int ivar=0;ivar<3;ivar++){

    for (int km=1; km<=kx; km++) {
      for (int jm=1; jm<=jx; jm++) {
        for (int im=1; im<=ix; im++) {
          int mip = _F_IDX_V3DEX(ivar, im, jm, km, ix, jx, kx, gd);//ivar=0:x方向、ivar=1:y方向、ivar=2:z方向
          double ddd=(double)data[mip];
          int i=im-1;
          int j=jm-1;
          int k=km-1;
          int ip1= k   *id*jd+ j   *id+i;
          int ip2= k   *id*jd+ j   *id+i+1;
          int ip3= k   *id*jd+(j+1)*id+i+1;
          int ip4= k   *id*jd+(j+1)*id+i;
          int ip5=(k+1)*id*jd+ j   *id+i;
          int ip6=(k+1)*id*jd+ j   *id+i+1;
          int ip7=(k+1)*id*jd+(j+1)*id+i+1;
          int ip8=(k+1)*id*jd+(j+1)*id+i;
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
void COMB::setVectorComponentGridData(int* size, int guide, double* d, float* data, int id, int jd, int kd, int ivar)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;

  for (int i=0; i<id*jd*kd; i++) d[i]=0.0;
 
  for (int km=1; km<=kx; km++) {
    for (int jm=1; jm<=jx; jm++) {
      for (int im=1; im<=ix; im++) {
        int mip = _F_IDX_V3DEX(ivar, im, jm, km, ix, jx, kx, gd);//ivar=0:x方向、ivar=1:y方向、ivar=2:z方向
        double ddd=(double)data[mip];
        int i=im-1;
        int j=jm-1;
        int k=km-1;
        int ip1= k   *id*jd+ j   *id+i;
        int ip2= k   *id*jd+ j   *id+i+1;
        int ip3= k   *id*jd+(j+1)*id+i+1;
        int ip4= k   *id*jd+(j+1)*id+i;
        int ip5=(k+1)*id*jd+ j   *id+i;
        int ip6=(k+1)*id*jd+ j   *id+i+1;
        int ip7=(k+1)*id*jd+(j+1)*id+i+1;
        int ip8=(k+1)*id*jd+(j+1)*id+i;
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
void COMB::VolumeDataDivideBy8(float* d, int id, int jd, int kd)
{
  int i,j,k,ip;
  for (k=1; k<kd-1; k++){
    for (j=1; j<jd-1; j++){
      for (i=1; i<id-1; i++){
        ip=k*id*jd+j*id+i;
        d[ip]=d[ip]*0.125;
      }
    }
  }
}

//面上の格子点のデータを4で割る
void COMB::FaceDataDivideBy4(float* d, int id, int jd, int kd)
{
  int i,j,k,ip;

  i=0;
  for (k=1; k<kd-1; k++){
    for (j=1; j<jd-1; j++){
      ip=k*id*jd+j*id+i;
      d[ip]=d[ip]*0.25;
    }
  }

  i=id-1;
  for (k=1; k<kd-1; k++){
    for (j=1; j<jd-1; j++){
      ip=k*id*jd+j*id+i;
      d[ip]=d[ip]*0.25;
    }
  }

  j=0;
  for (k=1; k<kd-1; k++){
    for (i=1; i<id-1; i++){
      ip=k*id*jd+j*id+i;
      d[ip]=d[ip]*0.25;
    }
  }

  j=jd-1;
  for (k=1; k<kd-1; k++){
    for (i=1; i<id-1; i++){
      ip=k*id*jd+j*id+i;
      d[ip]=d[ip]*0.25;
    }
  }

  k=0;
  for (j=1; j<jd-1; j++){
    for (i=1; i<id-1; i++){
      ip=k*id*jd+j*id+i;
      d[ip]=d[ip]*0.25;
    }
  }

  k=kd-1;
  for (j=1; j<jd-1; j++){
    for (i=1; i<id-1; i++){
      ip=k*id*jd+j*id+i;
      d[ip]=d[ip]*0.25;
    }
  }
}

//辺上の格子点のデータを2で割る
void COMB::LineDataDivideBy2(float* d, int id, int jd, int kd)
{
  int i,j,k,ip;

  i=0; j=0;
  for (k=1; k<kd-1; k++){
    ip=k*id*jd+j*id+i;
    d[ip]=d[ip]*0.5;
  }

  i=0; j=jd-1;
  for (k=1; k<kd-1; k++){
    ip=k*id*jd+j*id+i;
    d[ip]=d[ip]*0.5;
  }

  i=0; k=0;
  for (j=1; j<jd-1; j++){
    ip=k*id*jd+j*id+i;
    d[ip]=d[ip]*0.5;
  }

  i=0; k=kd-1;
  for (j=1; j<jd-1; j++){
    ip=k*id*jd+j*id+i;
    d[ip]=d[ip]*0.5;
  }

  j=0; k=0;
  for (i=1; i<id-1; i++){
    ip=k*id*jd+j*id+i;
    d[ip]=d[ip]*0.5;
  }

  j=0; k=kd-1;
  for (i=1; i<id-1; i++){
    ip=k*id*jd+j*id+i;
    d[ip]=d[ip]*0.5;
  }

  j=jd-1; k=0;
  for (i=1; i<id-1; i++){
    ip=k*id*jd+j*id+i;
    d[ip]=d[ip]*0.5;
  }

  j=jd-1; k=kd-1;
  for (i=1; i<id-1; i++){
    ip=k*id*jd+j*id+i;
    d[ip]=d[ip]*0.5;
  }

  i=id-1; j=0;
  for (k=1; k<kd-1; k++){
    ip=k*id*jd+j*id+i;
    d[ip]=d[ip]*0.5;
  }

  i=id-1; j=jd-1;
  for (k=1; k<kd-1; k++){
    ip=k*id*jd+j*id+i;
    d[ip]=d[ip]*0.5;
  }

  i=id-1; k=0;
  for (j=1; j<jd-1; j++){
    ip=k*id*jd+j*id+i;
    d[ip]=d[ip]*0.5;
  }

  i=id-1; k=kd-1;
  for (j=1; j<jd-1; j++){
    ip=k*id*jd+j*id+i;
    d[ip]=d[ip]*0.5;
  }

}

//*****************************************************************************
//*****************************************************************************
//*****************************************************************************

//内部の格子点のデータを8で割る
void COMB::VolumeDataDivideBy8(double* d, int id, int jd, int kd)
{
  int i,j,k,ip;
  for (k=1; k<kd-1; k++){
    for (j=1; j<jd-1; j++){
      for (i=1; i<id-1; i++){
        ip=k*id*jd+j*id+i;
        d[ip]=d[ip]*0.125;
      }
    }
  }
}

//面上の格子点のデータを4で割る
void COMB::FaceDataDivideBy4(double* d, int id, int jd, int kd)
{
  int i,j,k,ip;

  i=0;
  for (k=1; k<kd-1; k++){
    for (j=1; j<jd-1; j++){
      ip=k*id*jd+j*id+i;
      d[ip]=d[ip]*0.25;
    }
  }

  i=id-1;
  for (k=1; k<kd-1; k++){
    for (j=1; j<jd-1; j++){
      ip=k*id*jd+j*id+i;
      d[ip]=d[ip]*0.25;
    }
  }

  j=0;
  for (k=1; k<kd-1; k++){
    for (i=1; i<id-1; i++){
      ip=k*id*jd+j*id+i;
      d[ip]=d[ip]*0.25;
    }
  }

  j=jd-1;
  for (k=1; k<kd-1; k++){
    for (i=1; i<id-1; i++){
      ip=k*id*jd+j*id+i;
      d[ip]=d[ip]*0.25;
    }
  }

  k=0;
  for (j=1; j<jd-1; j++){
    for (i=1; i<id-1; i++){
      ip=k*id*jd+j*id+i;
      d[ip]=d[ip]*0.25;
    }
  }

  k=kd-1;
  for (j=1; j<jd-1; j++){
    for (i=1; i<id-1; i++){
      ip=k*id*jd+j*id+i;
      d[ip]=d[ip]*0.25;
    }
  }
}

//辺上の格子点のデータを2で割る
void COMB::LineDataDivideBy2(double* d, int id, int jd, int kd)
{
  int i,j,k,ip;

  i=0; j=0;
  for (k=1; k<kd-1; k++){
    ip=k*id*jd+j*id+i;
    d[ip]=d[ip]*0.5;
  }

  i=0; j=jd-1;
  for (k=1; k<kd-1; k++){
    ip=k*id*jd+j*id+i;
    d[ip]=d[ip]*0.5;
  }

  i=0; k=0;
  for (j=1; j<jd-1; j++){
    ip=k*id*jd+j*id+i;
    d[ip]=d[ip]*0.5;
  }

  i=0; k=kd-1;
  for (j=1; j<jd-1; j++){
    ip=k*id*jd+j*id+i;
    d[ip]=d[ip]*0.5;
  }

  j=0; k=0;
  for (i=1; i<id-1; i++){
    ip=k*id*jd+j*id+i;
    d[ip]=d[ip]*0.5;
  }

  j=0; k=kd-1;
  for (i=1; i<id-1; i++){
    ip=k*id*jd+j*id+i;
    d[ip]=d[ip]*0.5;
  }

  j=jd-1; k=0;
  for (i=1; i<id-1; i++){
    ip=k*id*jd+j*id+i;
    d[ip]=d[ip]*0.5;
  }

  j=jd-1; k=kd-1;
  for (i=1; i<id-1; i++){
    ip=k*id*jd+j*id+i;
    d[ip]=d[ip]*0.5;
  }

  i=id-1; j=0;
  for (k=1; k<kd-1; k++){
    ip=k*id*jd+j*id+i;
    d[ip]=d[ip]*0.5;
  }

  i=id-1; j=jd-1;
  for (k=1; k<kd-1; k++){
    ip=k*id*jd+j*id+i;
    d[ip]=d[ip]*0.5;
  }

  i=id-1; k=0;
  for (j=1; j<jd-1; j++){
    ip=k*id*jd+j*id+i;
    d[ip]=d[ip]*0.5;
  }

  i=id-1; k=kd-1;
  for (j=1; j<jd-1; j++){
    ip=k*id*jd+j*id+i;
    d[ip]=d[ip]*0.5;
  }

}


// #################################################################
// 
void COMB::output_plot3d_on_cell()
{
  if( !DI ) return;

  char tmp[FB_FILE_PATH_LENGTH];
  string prefix,outfile,infile;
  int d_type;

  int m_sv_type, m_d_type;
  int m_rank, m_step, m_imax, m_jmax, m_kmax;
  long long m_dstep, m_dimax, m_djmax, m_dkmax;
  int sz[3];
  long long dsz[3];
  float m_org[3], m_pit[3];
  float m_time;
  double m_dorg[3], m_dpit[3];
  double m_dtime;
  int dim;

  int fp_in =8;
  FILE *fp;

  // work area for read sph
  int vsize,wksize;
  int maxsize=0;
  int maxdim=0;
  float *wk;
  double *dwk;

  // work area for write plot3d
  int outsize;
  int maxid=0;
  int maxjd=0;
  int maxkd=0;
  int maxnvar=3;
  int nvar=ndfi;
  float *d;
  double *dd;
  int id,jd,kd;
  int ngrid=1;

// set loop count
  int nstep=DI[0].step.size();
  int nnode=DI[0].NodeInfoSize;

// copy dfi file
  string dfi_in = Generate_DFI_Name(DI[0].Prefix);
  string dfi_out = Generate_DFI_Name(P3Op.basename);
  dfi_out = dirname + dfi_out;
  Copy_DFIFile(dfi_in, dfi_out, P3Op.basename, dfi_mng[var_Plot3D]);
  if(P3Op.IS_DivideFunc == ON){ // 分割出力のとき
    for(int i=0;i<ndfi;i++){
      dfi_mng[var_Plot3D]=0;//本来はカウンターだが初期判定にのみ利用
      std::string dfipre;
      dfipre = DI[i].Prefix + P3Op.basename;
      dfi_in = Generate_DFI_Name(DI[i].Prefix);
      dfi_out = Generate_DFI_Name(dfipre);
      dfi_out = dirname + dfi_out;
      Copy_DFIFile(dfi_in, dfi_out, dfipre, dfi_mng[var_Plot3D]);
    }
  }

// check data type
  m_step = DI[0].step[0];
  m_rank = DI[0].Node[0].RankID;
  prefix = DI[0].Prefix;
  infile = Generate_FileName(prefix, m_step, m_rank, true);
  ReadSphDataType (&m_sv_type, &m_d_type, fp_in, infile);
  d_type=m_d_type;
  if( d_type == SPH_FLOAT ){// sphファイルがfloatの時、PLOT3D出力は必ずfloatにする
    FP3DR.setRealType(d_type);
    FP3DW.setRealType(d_type);
  }

// set work area size
  for(int istep=0; istep< nstep; istep++ ) { // step loop
    m_step=DI[0].step[istep];

    for(int inode=0; inode< nnode; inode++ ) { // node loop

      int ivar=0;
      for(int i=0;i<ndfi;i++){
        prefix=DI[i].Prefix;
        if     ( !strcasecmp(prefix.c_str(), "prs_" ) ) dim=1;
        else if( !strcasecmp(prefix.c_str(), "vel_" ) ) dim=3;
        else if( !strcasecmp(prefix.c_str(), "tmp_" ) ) dim=1;
        else if( !strcasecmp(prefix.c_str(), "tp_"  ) ) dim=1;
        else if( !strcasecmp(prefix.c_str(), "vrt_" ) ) dim=3;
        else if( !strcasecmp(prefix.c_str(), "i2vgt_" ) ) dim=1;
        else if( !strcasecmp(prefix.c_str(), "hlt_" ) ) dim=1;
	    else{
          m_step = DI[i].step[istep];
          m_rank = DI[i].Node[inode].RankID;
          infile = Generate_FileName(prefix, m_step, m_rank, true);
          ReadSphDataType (&m_sv_type, &m_d_type, fp_in, infile);
          dim=m_sv_type;
        }
        if(maxdim<dim) maxdim=dim;
        ivar=ivar+dim;
      }
      nvar=ivar;

      if(FP3DW.GetFormat() == C_BINARY) if(maxnvar<dim) maxnvar=dim;// C_BINARYでの出力は項目ごとに書き出し
      else                              if(maxnvar<ivar) maxnvar=ivar;// Fortranによる出力では出力項目すべてを一度に書き出し
      if( P3Op.IS_DivideFunc == ON ) maxnvar=maxdim; //項目別出力onの時

      // read sph header
      prefix = DI[0].Prefix;
      m_rank = DI[0].Node[inode].RankID;
      infile = Generate_FileName(prefix, m_step, m_rank, true);
      if( d_type == SPH_FLOAT ){
        ReadSphHeader (&m_step, &m_sv_type, &m_d_type, &m_imax, &m_jmax, &m_kmax, &m_time, m_org, m_pit, fp_in, infile);
      }else{
        ReadSphHeader (&m_dstep, &m_sv_type, &m_d_type, &m_dimax, &m_djmax, &m_dkmax, &m_dtime, m_dorg, m_dpit, fp_in, infile);
        m_step=(int)m_dstep;
        m_imax=(int)m_dimax;
        m_jmax=(int)m_djmax;
        m_kmax=(int)m_dkmax;
      }

      // set size
      int guide = DI[0].GuideCell;
      sz[0]=m_imax;
      sz[1]=m_jmax;
      sz[2]=m_kmax;
      int id,jd,kd;
      id=sz[0];//+1-2*guide;//PLOT3D出力はguideセルを書かないのでguideセル分減じておく
      jd=sz[1];//+1-2*guide;
      kd=sz[2];//+1-2*guide;
      wksize=sz[0]*sz[1]*sz[2];
      if(maxsize<wksize) maxsize=wksize;
      if(maxid<id) maxid=id;
      if(maxjd<jd) maxjd=jd;
      if(maxkd<kd) maxkd=kd;
    }
  }

// write memory size
  double mc; 

  // sph読み込みwk用メモリ計算
  mc = (double)maxsize*(double)maxdim;
  if(mc>(double)INT_MAX){// 整数値あふれ出しチェック //参考 894*894*894*3=2143550952 INT_MAX 2147483647
    printf("\tsize error : wkmaxsize>INT_MAX\n");
    Exit(0);
  }
  double sphMemory;
  if( d_type == SPH_FLOAT ) sphMemory = mc * (double)sizeof(float);
  else                      sphMemory = mc * (double)sizeof(double);

  // PLOT3D出力用wksizeのmax予測
  mc = (double)maxid*(double)maxjd*(double)maxkd*(double)maxnvar;
  if(mc>(double)INT_MAX){// 整数値あふれ出しチェック //参考 894*894*894*3=2143550952 INT_MAX 2147483647
    printf("\tsize error : psize>INT_MAX\n");
    Exit(0);
  }
  double plot3dMemory;
  if( FP3DW.GetRealType() == 1 ) plot3dMemory = mc * (double)sizeof(float);
  else                           plot3dMemory = mc * (double)sizeof(double);

  // メモリチェック
  double TotalMemory=sphMemory+plot3dMemory;
  LOG_OUT_ MemoryRequirement(TotalMemory,sphMemory,plot3dMemory,fplog);
  STD_OUT_ MemoryRequirement(TotalMemory,sphMemory,plot3dMemory,stdout);

// allocate work area
  unsigned long wkmaxsize = maxsize*maxdim;
  unsigned long psize = maxid*maxjd*maxkd*maxnvar;
  if( d_type == SPH_FLOAT ){
    if ( !(wk = new float[wkmaxsize]) ){
      printf("\tallocate error : wk\n");
      Exit(0);
    }
  }else{
    if ( !(dwk = new double[wkmaxsize]) ){
      printf("\tallocate error : dwk\n");
      Exit(0);
    }
  }
  if( FP3DW.GetRealType() == 1 ){
    if ( !(d = new float[psize]) ){
      printf("\tallocate error : wk\n");
      Exit(0);
    }
  }else{
    if ( !(dd = new double[psize]) ){
      printf("\tallocate error : wk\n");
      Exit(0);
    }
  }

// output xyz file
  LOG_OUT_ fprintf(fplog,"\t*** output xyz file ***\n");
  STD_OUT_ printf("\t*** output xyz file ***\n");

  // step loop
  for(int istep=0; istep< nstep; istep++ ) {
    m_step=DI[0].step[istep];
    LOG_OUT_ fprintf(fplog,"\tstep = %d\n", m_step);
    STD_OUT_ printf("\tstep = %d\n", m_step);

    // node loop
    for(int inode=0; inode< nnode; inode++ ) {

      // read sph header
      prefix=DI[0].Prefix;
      m_rank = DI[0].Node[inode].RankID;
      infile = Generate_FileName(prefix, m_step, m_rank, true);
      if( d_type == SPH_FLOAT ){
        ReadSphHeader (&m_step, &m_sv_type, &m_d_type, &m_imax, &m_jmax, &m_kmax, &m_time, m_org, m_pit, fp_in, infile);
      }else{
        ReadSphHeader (&m_dstep, &m_sv_type, &m_d_type, &m_dimax, &m_djmax, &m_dkmax, &m_dtime, m_dorg, m_dpit, fp_in, infile);
        m_step=(int)m_dstep;
        m_imax=(int)m_dimax;
        m_jmax=(int)m_djmax;
        m_kmax=(int)m_dkmax;
      }

      // set size
      int guide = DI[0].GuideCell;
      sz[0]=m_imax;
      sz[1]=m_jmax;
      sz[2]=m_kmax;
      id=sz[0];//1-2*guide;//+2*gc_out
      jd=sz[1];//+1-2*guide;//+2*gc_out
      kd=sz[2];//+1-2*guide;//+2*gc_out
      outsize=id*jd*kd;
      if(outsize>psize){
        printf("\tsize error : outsize>psize\n");
        Exit(0);
      }

      // write plot 3d xyz file
      if( d_type == SPH_FLOAT ){
        if( FP3DW.GetRealType() == 1 ){ // float ---> float
          OutputPlot3D_xyz_on_cell(m_step, m_rank, guide, m_org, m_pit, sz, &d[0], &d[outsize], &d[outsize*2]);
        }
      }

    } // node
  } // step

// output function file
  LOG_OUT_ fprintf(fplog,"\t*** output func file ***\n");
  STD_OUT_ printf("\t*** output func file ***\n");

  // step loop
  for(int istep=0; istep< nstep; istep++ ) {
    m_step=DI[0].step[istep];
    LOG_OUT_ fprintf(fplog,"\tstep = %d\n", m_step);
    STD_OUT_ printf("\tstep = %d\n", m_step);

    // node loop
    for(int inode=0; inode< nnode; inode++ ) {
      m_rank = DI[0].Node[inode].RankID;
      LOG_OUT_ fprintf(fplog,"\t  rank = %d\n", m_rank);
      STD_OUT_ printf("\t  rank = %d\n", m_rank);

      if(P3Op.IS_DivideFunc == OFF){ // 一括出力のとき、先にヘッダーを呼んでgridなどを書き出す

        // read sph header
        prefix=DI[0].Prefix;
        infile = Generate_FileName(prefix, m_step, m_rank, true);
        if( d_type == SPH_FLOAT ){
          ReadSphHeader (&m_step, &m_sv_type, &m_d_type, &m_imax, &m_jmax, &m_kmax, &m_time, m_org, m_pit, fp_in, infile);
        }else{
          ReadSphHeader (&m_dstep, &m_sv_type, &m_d_type, &m_dimax, &m_djmax, &m_dkmax, &m_dtime, m_dorg, m_dpit, fp_in, infile);
          m_step=(int)m_dstep;
          m_imax=(int)m_dimax;
          m_jmax=(int)m_djmax;
          m_kmax=(int)m_dkmax;
        }

        // set size
        int guide = DI[0].GuideCell;
        sz[0]=m_imax;
        sz[1]=m_jmax;
        sz[2]=m_kmax;
        id=sz[0];//+1-2*guide;//+2*gc_out
        jd=sz[1];//+1-2*guide;//+2*gc_out
        kd=sz[2];//+1-2*guide;//+2*gc_out
        outsize=id*jd*kd;
        if(outsize>psize){
          printf("\tsize error : outsize>psize\n");
          Exit(0);
        }

        // plot3d func file open
        std::string ptmp;
        ptmp = Generate_FileName_Free(P3Op.basename, "func", m_step, m_rank, true);
        ptmp = dirname + ptmp;
        FP3DW.setFileName(ptmp.c_str());
        if(!FP3DW.OpenFile()){
          printf("\terror OpenFile : %s\n", ptmp.c_str());
          Exit(0);
        }

        // write block data
        FP3DW.WriteNgrid(ngrid);
        FP3DW.WriteFuncBlockData(&id,&jd,&kd,&nvar,ngrid);
      }

      // component loop --- dfi file ---> prs_, vel, ,,,
      int ivar=0;
      for(int i=0;i<ndfi;i++){
        prefix=DI[i].Prefix;
        LOG_OUTV_ fprintf(fplog,"\t    COMPONENT : %s\n", prefix.c_str());
        STD_OUTV_ printf("\t    COMPONENT : %s\n", prefix.c_str());
  
        // Scalar or Vector
        if     ( !strcasecmp(prefix.c_str(), "prs_" ) ) dim=1;
        else if( !strcasecmp(prefix.c_str(), "vel_" ) ) dim=3;
        else if( !strcasecmp(prefix.c_str(), "tmp_" ) ) dim=1;
        else if( !strcasecmp(prefix.c_str(), "tp_"  ) ) dim=1;
        else if( !strcasecmp(prefix.c_str(), "vrt_" ) ) dim=3;
        else if( !strcasecmp(prefix.c_str(), "i2vgt_" ) ) dim=1;
        else if( !strcasecmp(prefix.c_str(), "hlt_" ) ) dim=1;
	    else{
          m_step = DI[i].step[istep];
          m_rank = DI[i].Node[inode].RankID;
          infile = Generate_FileName(prefix, m_step, m_rank, true);
          ReadSphDataType (&m_sv_type, &m_d_type, fp_in, infile);
          dim=m_sv_type;
        }

        // read sph header
        infile = Generate_FileName(prefix, m_step, m_rank, true);
        if( d_type == SPH_FLOAT ){
          ReadSphHeader (&m_step, &m_sv_type, &m_d_type, &m_imax, &m_jmax, &m_kmax, &m_time, m_org, m_pit, fp_in, infile);
        }else{
          ReadSphHeader (&m_dstep, &m_sv_type, &m_d_type, &m_dimax, &m_djmax, &m_dkmax, &m_dtime, m_dorg, m_dpit, fp_in, infile);
          m_step=(int)m_dstep;
          m_time=(float)m_dtime;
          m_imax=(int)m_dimax;
          m_jmax=(int)m_djmax;
          m_kmax=(int)m_dkmax;
        }

        // set size
        int guide = DI[0].GuideCell;
        sz[0]=m_imax;
        sz[1]=m_jmax;
        sz[2]=m_kmax;
        id=sz[0];//+1-2*guide;//+2*gc_out
        jd=sz[1];//+1-2*guide;//+2*gc_out
        kd=sz[2];//+1-2*guide;//+2*gc_out
        vsize=sz[0]*sz[1]*sz[2];
        wksize = vsize*dim;
        outsize=id*jd*kd;
        if(wksize>maxsize*maxdim){
          printf("\tsize error : wksize>maxsize*maxdim\n");
          Exit(0);
        }
        if(outsize>psize){
          printf("\tsize error : outsize>psize\n");
          Exit(0);
        }

        // read sph data
        if( d_type == SPH_FLOAT ) ReadSphData ( wk, wksize, sz, dim, fp_in, infile );
        else                      ReadSphData ( dwk, wksize, sz, dim, fp_in, infile );

        if(P3Op.IS_DivideFunc == ON){ // 分割出力のとき
       
          // plot3d func file open
          std::string ptmp;
          ptmp = Generate_FileName_Free(P3Op.basename, "func", m_step, m_rank, true);
          ptmp = DI[i].Prefix + ptmp;
          ptmp = dirname + ptmp;
          FP3DW.setFileName(ptmp.c_str());
          if(!FP3DW.OpenFile()){
            printf("\terror OpenFile : %s\n", ptmp.c_str());
            Exit(0);
          }
       
          // write block data
          ivar=0;
          nvar=dim;
          FP3DW.WriteNgrid(ngrid);
          FP3DW.WriteFuncBlockData(&id,&jd,&kd,&nvar,ngrid);
        }

        // write function file
        FP3DW.setGridData(id,jd,kd,ngrid);
        FP3DW.setFuncDataNum(nvar);
        if( dim==1 ){
          if(FP3DW.GetFormat() == C_BINARY){ // C_BINARYでの出力は項目ごとに書き出し
            if( d_type == SPH_FLOAT ){
              if( FP3DW.GetRealType() == 1 ){ // float ---> float
                for(id=0;id<outsize;id++) { d[id]=wk[id]; }
                FP3DW.setFuncData(d);
              }
              if(!FP3DW.WriteFuncData()) printf("\terror WriteFuncData\n");
            }
          }
          else // Fortranによる出力では出力項目すべてを一度に書き出し
          {
            if( d_type == SPH_FLOAT ){
              if( FP3DW.GetRealType() == 1 ){ // float ---> float
                for(id=0;id<outsize;id++) { d[outsize*ivar+id]=wk[id]; }
              }
            }
          }
        }else if( dim==3 ){
          if(FP3DW.GetFormat() == C_BINARY){ // C_BINARYでの出力は項目ごとに書き出し
            if( d_type == SPH_FLOAT ){
              if( FP3DW.GetRealType() == 1 ){ // float ---> float
                for(int ixyz=0;ixyz<3;ixyz++){
                  for(id=0;id<outsize;id++) { 
                    d[id]=wk[id*3+ixyz]; 
                  }
                  FP3DW.setFuncData(d);
                  if(!FP3DW.WriteFuncData()) printf("\terror WriteFuncData\n");
                }
              }
            }
          }
          else // Fortranによる出力では出力項目すべてを一度に書き出し
          {
            if( d_type == SPH_FLOAT ){
              if( FP3DW.GetRealType() == 1 ){ // float ---> float
                for(int ixyz=0;ixyz<3;ixyz++){
                  for(id=0;id<outsize;id++) { 
                    d[outsize*(ivar+ixyz)+id]=wk[id*3+ixyz]; 
                  }
                }
              }
            }
          }
        }
        
        if(P3Op.IS_DivideFunc == ON){ // 分割出力のとき
          // write all
          if(FP3DW.GetFormat() != C_BINARY){ // C_BINARY以外は出力項目すべてを一度に書き出し
            if( FP3DW.GetRealType() == 1 ) FP3DW.setFuncData(d);
            else                           FP3DW.setFuncData(dd);
            if(!FP3DW.WriteFuncData()) printf("\terror WriteFuncData\n");
          }
          // plot3d func file close
          FP3DW.CloseFile();

          //dfiファイルの出力
          std::string dfipre;
          dfipre = DI[i].Prefix + P3Op.basename;
          dfipre = dirname + dfipre;
          if(inode == 0) if ( !DFI.Write_DFI_File(dfipre, m_step, (double)m_time, dfi_mng[var_Plot3D], true) ) Exit(0);
        }

        ivar=ivar+dim;
      }

      if(P3Op.IS_DivideFunc == OFF){ // 一括出力のとき
        // write all
        if(FP3DW.GetFormat() != C_BINARY){//C_BINARY以外は出力項目すべてを一度に書き出し
          if( FP3DW.GetRealType() == 1 ) FP3DW.setFuncData(d);
          else FP3DW.setFuncData(dd);
          if(!FP3DW.WriteFuncData()) printf("\terror WriteFuncData\n");
        }
        // plot3d func file close
        FP3DW.CloseFile();
      }

      //dfiファイルの出力
      std::string dfipre;
      dfipre = dirname + P3Op.basename;
      if(inode == 0) if ( !DFI.Write_DFI_File(dfipre, m_step, (double)m_time, dfi_mng[var_Plot3D], true) ) Exit(0);

    }//node loop
  }//step loop

// deallocate work area
  if( d_type == SPH_FLOAT ) delete[] wk;
  else                      delete[] dwk;
  if( FP3DW.GetRealType() == 1 ) delete[] d;
  else                           delete[] dd;

// output function name file
  LOG_OUT_ fprintf(fplog,"\t*** output func name file ***\n");
  STD_OUT_ printf("\t*** output func name file ***\n");

  // function_nameファイルはかならずformatted形式
  int keep_format=FP3DW.GetFormat();
  FP3DW.setFormat(FORMATTED);
  
  // open file
  if(P3Op.IS_DivideFunc == OFF){ // 一括出力のとき
    std::string ptmp;
    ptmp = DFI.Generate_FileName_Free(P3Op.basename, "nam", 0, 0, false);
    ptmp = dirname + ptmp;
    FP3DW.setFileName(ptmp.c_str());
    if(!FP3DW.OpenFile()){
      printf("Error : error OpenFile\n");
      Exit(0);
    }
  }

  int unknowncomp=0;
  for(int i=0;i<ndfi;i++){
    string comp = DI[i].Prefix;
    LOG_OUT_ fprintf(fplog,"\t    COMPONENT : %s\n", comp.c_str());
    STD_OUT_ printf("\t    COMPONENT : %s\n", comp.c_str());

    // open file
    if(P3Op.IS_DivideFunc == ON){ // 分割出力のとき
      std::string ptmp;
      ptmp = Generate_FileName_Free(P3Op.basename, "nam", 0, 0, true);
      ptmp = DI[i].Prefix + ptmp;
      ptmp = dirname + ptmp;
      FP3DW.setFileName(ptmp.c_str());
      if(!FP3DW.OpenFile()){
        printf("\terror OpenFile : %s\n", ptmp.c_str());
        Exit(0);
      }
    }

    // write comp name
    if     ( !strcasecmp(comp.c_str(), "prs_" ) ){
      FP3DW.WriteFunctionName("Pressure");
    }else if( !strcasecmp(comp.c_str(), "vel_" ) ){
      FP3DW.WriteFunctionName("U-Velocity ; Velocity");
      FP3DW.WriteFunctionName("V-Velocity");
      FP3DW.WriteFunctionName("W-Velocity");
    }else if( !strcasecmp(comp.c_str(), "tmp_" ) ){
      FP3DW.WriteFunctionName("Tempearture");
    }else if( !strcasecmp(comp.c_str(), "tp_"  ) ){
      FP3DW.WriteFunctionName("Total_Pressure");
    }else if( !strcasecmp(comp.c_str(), "vrt_" ) ){
      FP3DW.WriteFunctionName("U-Vorticity ; Vorticity");
      FP3DW.WriteFunctionName("V-Vorticity");
      FP3DW.WriteFunctionName("W-Vorticity");
    }else if( !strcasecmp(comp.c_str(), "i2vgt_"  ) ){
      FP3DW.WriteFunctionName("2nd Invariant of Velocity Gradient Tensor");
    }else if( !strcasecmp(comp.c_str(), "hlt_" ) ){
      FP3DW.WriteFunctionName("Helicity");
    }else{
#if 0
      unknowncomp++;
      char* numbuff = new char[DFI_LINE_LENGTH];
      memset(numbuff, 0, sizeof(char)*DFI_LINE_LENGTH);
      sprintf(numbuff, "%02d", unknowncomp);
      string compbuff=(numbuff);
      compbuff ="F-"+compbuff;
      FP3DW.WriteFunctionName(compbuff.c_str());
      if ( numbuff ) delete [] numbuff;
#else
      FP3DW.WriteFunctionName(comp.c_str());
#endif
    }

    // close file
    if(P3Op.IS_DivideFunc == ON) FP3DW.CloseFile();

  }

  // close file
  if(P3Op.IS_DivideFunc == OFF) FP3DW.CloseFile();
  
  // reset option
  FP3DW.setFormat(keep_format);
  
}


//
void COMB::OutputPlot3D_xyz_on_cell(int m_step, int m_rank, int guide, float* origin, float* pitch, int* size, float* x, float* y, float* z)
{
  //value
  int ngrid=1;
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  int gc_out = 0;//plot3dは常にガイドセルを出力しない

  int *iblank=NULL;//dummy
  int id,jd,kd;//出力サイズ
  id=size[0];//+1;//+2*gc_out
  jd=size[1];//+1;//+2*gc_out
  kd=size[2];//+1;//+2*gc_out

  // ガイドセル出力があった場合オリジナルポイントを調整しておく
  float m_org[3], m_pit[3];
  for (int i=0; i<3; i++) 
  {
    m_org[i] = (float)origin[i] + (float)pitch[i]*(float)gd;
    m_pit[i] = (float)pitch[i];
  }

  // 出力ファイル名
  std::string tmp;
  tmp = Generate_FileName_Free(P3Op.basename, "xyz", m_step, m_rank, true);
  tmp = dirname + tmp;

  //open file
  FP3DW.setFileName(tmp.c_str());
  if(!FP3DW.OpenFile()){
    printf("Error : error OpenFile\n");
    Exit(0);
  }

  //write block data
  FP3DW.WriteNgrid(ngrid);//if multi grid
  FP3DW.WriteBlockData(id,jd,kd);

  for(int k=0;k<kd;k++){
    for(int j=0;j<jd;j++){
      for(int i=0;i<id;i++){
        int ip=k*id*jd+j*id+i;
        //int ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
        x[ip]=m_org[0]+m_pit[0]*(float)i;//-pitch[0]*(float)gc_out;
        y[ip]=m_org[1]+m_pit[1]*(float)j;//-pitch[1]*(float)gc_out;
        z[ip]=m_org[2]+m_pit[2]*(float)k;//-pitch[2]*(float)gc_out;
      }
    }
  }

  //write
  FP3DW.setGridData(id,jd,kd,ngrid);
  FP3DW.setXYZData(x,y,z,iblank);
  if(!FP3DW.WriteXYZData()) printf("\terror WriteXYZData\n");
  
  //close file
  FP3DW.CloseFile();
}

