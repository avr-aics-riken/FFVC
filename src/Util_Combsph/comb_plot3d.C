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
  
  int guide_out=0;
  
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
  size_t outsize;
  int maxid=0;
  int maxjd=0;
  int maxkd=0;
  int maxnvar=3;
  int nvar=ndfi;
  float *d;
  double *dd;
  int id,jd,kd;
  int ngrid=1;
  
  // 間引きのための変数
  size_t outsize_thin;
  int id_thin,jd_thin,kd_thin;
  int irest,jrest,krest;
  float *d_thin;
  double *dd_thin;
  
  
  // set loop count
  int nstep=DI[0].step.size();
  int nnode=DI[0].NodeInfoSize;
  
  
  //並列処理のためのインデックス作成
  int iproc=0;
  int ips=0;
  int ip;
  if(numProc > 1) {
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
  }
  
  // copy dfi file
  string dfi_in = Generate_DFI_Name(DI[0].Prefix);
  string dfi_out = Generate_DFI_Name(P3Op.basename);
  //dfi_out = out_dirname + dfi_out; // dfiファイルは他のdfiファイルと同じディレクトリにコピー
  Copy_DFIFile(dfi_in, dfi_out, P3Op.basename, dfi_mng[var_Plot3D]);
  if(P3Op.IS_DivideFunc == ON){ // 分割出力のとき
    for(int i=0;i<ndfi;i++){
      dfi_mng[var_Plot3D]=0;//本来はカウンターだが初期判定にのみ利用
      std::string dfipre;
      dfipre = DI[i].Prefix + P3Op.basename;
      cout << "dfipre = " << dfipre << endl;
      dfi_in = Generate_DFI_Name(DI[i].Prefix);
      dfi_out = Generate_DFI_Name(dfipre);
      //dfi_out = out_dirname + dfi_out; // dfiファイルは他のdfiファイルと同じディレクトリにコピー
      Copy_DFIFile(dfi_in, dfi_out, dfipre, dfi_mng[var_Plot3D]);
    }
  }
  
  // check data type
  m_step = DI[0].step[0];
  m_rank = DI[0].Node[0].RankID;
  prefix = DI[0].Prefix;
  infile = Generate_FileName(prefix, m_step, m_rank, true);
  infile = in_dirname + infile;
  ReadSphDataType (&m_sv_type, &m_d_type, fp_in, infile);
  d_type=m_d_type;
  if( d_type == SPH_FLOAT ){// sphファイルがfloatの時、PLOT3D出力は必ずfloatにする
    FP3DR.setRealType(d_type);
    FP3DW.setRealType(d_type);
  }
  
  // set work area size
  if(numProc > 1) ips=0;
  for(int istep=0; istep< nstep; istep++ ) { // step loop
    m_step=DI[0].step[istep];
    
    for(int inode=0; inode< nnode; inode++ ) { // node loop
      
      //並列処理
      int iskip=1;
      if(numProc > 1) {
        iskip=0;
        ip=ips+inode;
        for(int ic=0; ic< index.size(); ic++ ) {
          if(ip==index[ic]) iskip=1;
        }
      }
      if(!iskip) continue;
      
      int ivar=0;
      for(int i=0;i<ndfi;i++){
        prefix=DI[i].Prefix;
        if     ( !strcasecmp(prefix.c_str(), "prs" ) ) dim=1;
        else if( !strcasecmp(prefix.c_str(), "vel" ) ) dim=3;
        else if( !strcasecmp(prefix.c_str(), "fvel" )) dim=3;
        else if( !strcasecmp(prefix.c_str(), "tmp" ) ) dim=1;
        else if( !strcasecmp(prefix.c_str(), "tp"  ) ) dim=1;
        else if( !strcasecmp(prefix.c_str(), "vrt" ) ) dim=3;
        else if( !strcasecmp(prefix.c_str(), "i2vgt")) dim=1;
        else if( !strcasecmp(prefix.c_str(), "hlt" ) ) dim=1;
        else{
          m_step = DI[i].step[istep];
          m_rank = DI[i].Node[inode].RankID;
          infile = Generate_FileName(prefix, m_step, m_rank, true);
          infile = in_dirname + infile;
          ReadSphDataType (&m_sv_type, &m_d_type, fp_in, infile);
          dim=1;
          if( m_sv_type == SPH_VECTOR ) dim = 3;
        }
        if(maxdim<dim) maxdim=dim;
        ivar=ivar+dim;
        DI[i].dim=dim;
      }
      nvar=ivar;
      
      if(FP3DW.GetFormat() == C_BINARY) if(maxnvar<dim) maxnvar=dim;  // C_BINARYでの出力は項目ごとに書き出し
      else                              if(maxnvar<ivar) maxnvar=ivar;// Fortranによる出力では出力項目すべてを一度に書き出し
      if( P3Op.IS_DivideFunc == ON ) maxnvar=maxdim; //項目別出力onの時
      
      // read sph header
      prefix = DI[0].Prefix;
      m_rank = DI[0].Node[inode].RankID;
      infile = Generate_FileName(prefix, m_step, m_rank, true);
      infile = in_dirname + infile;
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
    if(numProc > 1) ips=ips+nnode;
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
    printf("\tsize error : mc>INT_MAX\n");
    Exit(0);
  }
  double plot3dMemory;
  if( FP3DW.GetRealType() == OUTPUT_FLOAT ) plot3dMemory = mc * (double)sizeof(float);
  else                                      plot3dMemory = mc * (double)sizeof(double);
  
  // 間引きのためのワークエリア
  int maxid_thin=(maxid-1)/thin_count+1;
  int maxjd_thin=(maxjd-1)/thin_count+1;
  int maxkd_thin=(maxkd-1)/thin_count+1;
  mc = (double)maxjd_thin*(double)maxjd_thin*(double)maxkd_thin*(double)maxnvar;
  if(mc>(double)INT_MAX){// 整数値あふれ出しチェック --- 参考 894*894*894*3=2143550952 INT_MAX 2147483647
    printf("\tsize error : mc>INT_MAX\n");
    Exit(0);
  }
  double thinMemory=0.0;
  if(thin_out){ // 間引きオプションが指定されたときのみアロケート
    if( FP3DW.GetRealType() == OUTPUT_FLOAT ) thinMemory = mc * (double)sizeof(float);
    else                                      thinMemory = mc * (double)sizeof(double);
  }
  
  // メモリチェック
  double TotalMemory=sphMemory+plot3dMemory+thinMemory;
  LOG_OUT_ MemoryRequirement(TotalMemory,sphMemory,plot3dMemory,thinMemory,fplog);
  STD_OUT_ MemoryRequirement(TotalMemory,sphMemory,plot3dMemory,thinMemory,stdout);
  
  // allocate work area
  size_t wkmaxsize = (size_t)(maxsize*maxdim);
  size_t psize = (size_t)(maxid*maxjd*maxkd*maxnvar);
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
  
  size_t psize_thin = (size_t)(maxid_thin*maxjd_thin*maxkd_thin*maxnvar);
  if(thin_out){
    if( FP3DW.GetRealType() == OUTPUT_FLOAT ){
      if ( !(d_thin = new float[psize_thin]) ){
        printf("\tallocate error : wk\n");
        Exit(0);
      }
    }else{
      if ( !(dd_thin = new double[psize_thin]) ){
        printf("\tallocate error : wk\n");
        Exit(0);
      }
    }
  }
  
  // output xyz file
  LOG_OUT_ fprintf(fplog,"\t*** output xyz file ***\n");
  STD_OUT_ printf("\t*** output xyz file ***\n");
  
  // step loop
  if(numProc > 1) ips=0;
  for(int istep=0; istep< nstep; istep++ ) {
    m_step=DI[0].step[istep];
    LOG_OUT_ fprintf(fplog,"\tstep = %d\n", m_step);
    STD_OUT_ printf("\tstep = %d\n", m_step);
    
    // node loop
    for(int inode=0; inode< nnode; inode++ ) {
      
      //並列処理
      int iskip=1;
      if(numProc > 1) {
        iskip=0;
        ip=ips+inode;
        for(int ic=0; ic< index.size(); ic++ ) {
          if(ip==index[ic]) iskip=1;
        }
      }
      if(!iskip) continue;
      
      // read sph header
      prefix=DI[0].Prefix;
      m_rank = DI[0].Node[inode].RankID;
      infile = Generate_FileName(prefix, m_step, m_rank, true);
      infile = in_dirname + infile;
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
      outsize=(size_t)id*(size_t)jd*(size_t)kd;
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
    if(numProc > 1) ips=ips+nnode;
  } // step
  
  // output function file
  LOG_OUT_ fprintf(fplog,"\t*** output func file ***\n");
  STD_OUT_ printf("\t*** output func file ***\n");
  
  // step loop
  if(numProc > 1) ips=0;
  for(int istep=0; istep< nstep; istep++ ) {
    m_step=DI[0].step[istep];
    LOG_OUT_ fprintf(fplog,"\tstep = %d\n", m_step);
    STD_OUT_ printf("\tstep = %d\n", m_step);
    
    // node loop
    for(int inode=0; inode< nnode; inode++ ) {
      m_rank = DI[0].Node[inode].RankID;
      LOG_OUT_ fprintf(fplog,"\t  rank = %d\n", m_rank);
      STD_OUT_ printf("\t  rank = %d\n", m_rank);
      
      //並列処理
      int iskip=1;
      if(numProc > 1) {
        iskip=0;
        ip=ips+inode;
        for(int ic=0; ic< index.size(); ic++ ) {
          if(ip==index[ic]) iskip=1;
        }
      }
      if(!iskip) continue;
      
      if(P3Op.IS_DivideFunc == OFF){ // 一括出力のとき、先にヘッダーを呼んでgridなどを書き出す
        
        // read sph header
        prefix=DI[0].Prefix;
        infile = Generate_FileName(prefix, m_step, m_rank, true);
        infile = in_dirname + infile;
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
        outsize=(size_t)id*(size_t)jd*(size_t)kd;
        if(outsize>psize){
          printf("\tsize error : outsize>psize\n");
          Exit(0);
        }
        
        //間引きのための処理
        irest=(id-1)%thin_count;
        jrest=(jd-1)%thin_count;
        krest=(kd-1)%thin_count;
        id_thin=(id-1)/thin_count;
        jd_thin=(jd-1)/thin_count;
        kd_thin=(kd-1)/thin_count;
        id_thin=id_thin+1;
        jd_thin=jd_thin+1;
        kd_thin=kd_thin+1;
        if(irest!=0) id_thin=id_thin+1;
        if(jrest!=0) jd_thin=jd_thin+1;
        if(krest!=0) kd_thin=kd_thin+1;
        outsize_thin=(size_t)id_thin*(size_t)jd_thin*(size_t)kd_thin;
        if(outsize_thin>psize_thin){
          printf("\tsize error : outsize_thin>psize_thin\n");
          Exit(0);
        }
        
        // plot3d func file open
        std::string ptmp;
        ptmp = Generate_FileName_Free(P3Op.basename, "func", m_step, m_rank, true);
        ptmp = out_dirname + ptmp;
        FP3DW.setFileName(ptmp.c_str());
        if(!FP3DW.OpenFile()){
          printf("\terror OpenFile : %s\n", ptmp.c_str());
          Exit(0);
        }
        
        // write block data
        FP3DW.WriteNgrid(ngrid);
        //FP3DW.WriteFuncBlockData(&id,&jd,&kd,&nvar,ngrid);
        FP3DW.WriteFuncBlockData(&id_thin,&jd_thin,&kd_thin,&nvar,ngrid);
      }
      
      // component loop --- dfi file ---> prs_, vel, ,,,
      int ivar=0;
      for(int i=0;i<ndfi;i++){
        prefix=DI[i].Prefix;
        LOG_OUTV_ fprintf(fplog,"\t    COMPONENT : %s\n", prefix.c_str());
        STD_OUTV_ printf("\t    COMPONENT : %s\n", prefix.c_str());
        
        // Scalar or Vector
        if     ( !strcasecmp(prefix.c_str(), "prs" ) ) dim=1;
        else if( !strcasecmp(prefix.c_str(), "vel" ) ) dim=3;
        else if( !strcasecmp(prefix.c_str(), "fvel" )) dim=3;
        else if( !strcasecmp(prefix.c_str(), "tmp" ) ) dim=1;
        else if( !strcasecmp(prefix.c_str(), "tp"  ) ) dim=1;
        else if( !strcasecmp(prefix.c_str(), "vrt" ) ) dim=3;
        else if( !strcasecmp(prefix.c_str(), "i2vgt")) dim=1;
        else if( !strcasecmp(prefix.c_str(), "hlt" ) ) dim=1;
        else{
          m_step = DI[i].step[istep];
          m_rank = DI[i].Node[inode].RankID;
          infile = Generate_FileName(prefix, m_step, m_rank, true);
          infile = in_dirname + infile;
          ReadSphDataType (&m_sv_type, &m_d_type, fp_in, infile);
          dim=1;
          if( m_sv_type == SPH_VECTOR ) dim = 3;
        }
        
        // read sph header
        infile = Generate_FileName(prefix, m_step, m_rank, true);
        infile = in_dirname + infile;
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
        outsize=(size_t)id*(size_t)jd*(size_t)kd;
        if(wksize>maxsize*maxdim){
          printf("\tsize error : wksize>maxsize*maxdim\n");
          Exit(0);
        }
        if(outsize>psize){
          printf("\tsize error : outsize>psize\n");
          Exit(0);
        }
        
        //間引きのための処理
        irest=(id-1)%thin_count;
        jrest=(jd-1)%thin_count;
        krest=(kd-1)%thin_count;
        id_thin=(id-1)/thin_count;
        jd_thin=(jd-1)/thin_count;
        kd_thin=(kd-1)/thin_count;
        id_thin=id_thin+1;
        jd_thin=jd_thin+1;
        kd_thin=kd_thin+1;
        if(irest!=0) id_thin=id_thin+1;
        if(jrest!=0) jd_thin=jd_thin+1;
        if(krest!=0) kd_thin=kd_thin+1;
        outsize_thin=(size_t)id_thin*(size_t)jd_thin*(size_t)kd_thin;
        if(outsize_thin>psize_thin){
          printf("\tsize error : outsize_thin>psize_thin\n");
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
          ptmp = out_dirname + ptmp;
          FP3DW.setFileName(ptmp.c_str());
          if(!FP3DW.OpenFile()){
            printf("\terror OpenFile : %s\n", ptmp.c_str());
            Exit(0);
          }
          
          // write block data
          ivar=0;
          nvar=dim;
          FP3DW.WriteNgrid(ngrid);
          //FP3DW.WriteFuncBlockData(&id,&jd,&kd,&nvar,ngrid);
          FP3DW.WriteFuncBlockData(&id_thin,&jd_thin,&kd_thin,&nvar,ngrid);
        }
        
        // write function file
        //FP3DW.setGridData(id,jd,kd,ngrid);
        FP3DW.setGridData(id_thin,jd_thin,kd_thin,ngrid);
        FP3DW.setFuncDataNum(nvar);
        if( dim==1 ){
          if(FP3DW.GetFormat() == C_BINARY){ // C_BINARYでの出力は項目ごとに書き出し
            if( d_type == SPH_FLOAT ){
              if( FP3DW.GetRealType() == OUTPUT_FLOAT ){ // float ---> float
                setScalarGridData(sz,guide,d,wk,id,jd,kd,guide_out);
                if(thin_out){
                  thinout_plot3d(d_thin,d,id_thin,jd_thin,kd_thin,id,jd,kd);
                  FP3DW.setFuncData(d_thin);
                }else{
                  FP3DW.setFuncData(d);
                }
              }else{ // float ---> double
                setScalarGridData(sz,guide,dd,wk,id,jd,kd,guide_out);
                if(thin_out){
                  thinout_plot3d(dd_thin,dd,id_thin,jd_thin,kd_thin,id,jd,kd);
                  FP3DW.setFuncData(dd_thin);
                }else{
                  FP3DW.setFuncData(dd);
                }
              }
              if(!FP3DW.WriteFuncData()) printf("\terror WriteFuncData\n");
            }else{
              if( FP3DW.GetRealType() == OUTPUT_FLOAT ){ // double ---> float
                setScalarGridData(sz,guide,d,dwk,id,jd,kd,guide_out);
                if(thin_out){
                  thinout_plot3d(d_thin,d,id_thin,jd_thin,kd_thin,id,jd,kd);
                  FP3DW.setFuncData(d_thin);
                }else{
                  FP3DW.setFuncData(d);
                }
              }else{ // double ---> double
                setScalarGridData(sz,guide,dd,dwk,id,jd,kd,guide_out);
                if(thin_out){
                  thinout_plot3d(dd_thin,dd,id_thin,jd_thin,kd_thin,id,jd,kd);
                  FP3DW.setFuncData(dd_thin);
                }else{
                  FP3DW.setFuncData(dd);
                }
              }
              if(!FP3DW.WriteFuncData()) printf("\terror WriteFuncData\n");
            }
          }
          else // Fortranによる出力では出力項目すべてを一度に書き出し
          {
            if( d_type == SPH_FLOAT ){
              if( FP3DW.GetRealType() == OUTPUT_FLOAT ){ // float ---> float
                setScalarGridData(sz,guide,&d[outsize*ivar],wk,id,jd,kd,guide_out);
                if(thin_out){
                  thinout_plot3d(&d_thin[outsize_thin*ivar],&d[outsize*ivar],id_thin,jd_thin,kd_thin,id,jd,kd);
                }
              }else{ // float ---> double
                setScalarGridData(sz,guide,&dd[outsize*ivar],wk,id,jd,kd,guide_out);
                if(thin_out){
                  thinout_plot3d(&dd_thin[outsize_thin*ivar],&dd[outsize*ivar],id_thin,jd_thin,kd_thin,id,jd,kd);
                }
              }
            }else{
              if( FP3DW.GetRealType() == OUTPUT_FLOAT ){ // double ---> float
                setScalarGridData(sz,guide,&d[outsize*ivar],dwk,id,jd,kd,guide_out);
                if(thin_out){
                  thinout_plot3d(&d_thin[outsize_thin*ivar],&d[outsize*ivar],id_thin,jd_thin,kd_thin,id,jd,kd);
                }
              }else{ // double ---> double
                setScalarGridData(sz,guide,&dd[outsize*ivar],dwk,id,jd,kd,guide_out);
                if(thin_out){
                  thinout_plot3d(&dd_thin[outsize_thin*ivar],&dd[outsize*ivar],id_thin,jd_thin,kd_thin,id,jd,kd);
                }
              }
            }
          }
        }else if( dim==3 ){
          if(FP3DW.GetFormat() == C_BINARY){ // C_BINARYでの出力は項目ごとに書き出し
            if( d_type == SPH_FLOAT ){
              if( FP3DW.GetRealType() == OUTPUT_FLOAT ){ // float ---> float
                for(int ixyz=0;ixyz<3;ixyz++){
                  setVectorComponentGridData(sz,guide,d,wk,id,jd,kd,ixyz,guide_out);
                  if(thin_out){
                    thinout_plot3d(d_thin,d,id_thin,jd_thin,kd_thin,id,jd,kd);
                    FP3DW.setFuncData(d_thin);
                  }else{
                    FP3DW.setFuncData(d);
                  }
                  if(!FP3DW.WriteFuncData()) printf("\terror WriteFuncData\n");
                }
              }else{ // float ---> double
                for(int ixyz=0;ixyz<3;ixyz++){
                  setVectorComponentGridData(sz,guide,dd,wk,id,jd,kd,ixyz,guide_out);
                  if(thin_out){
                    thinout_plot3d(dd_thin,dd,id_thin,jd_thin,kd_thin,id,jd,kd);
                    FP3DW.setFuncData(dd_thin);
                  }else{
                    FP3DW.setFuncData(dd);
                  }
                  if(!FP3DW.WriteFuncData()) printf("\terror WriteFuncData\n");
                }
              }
            }else{ // double ---> float
              if( FP3DW.GetRealType() == OUTPUT_FLOAT ){
                for(int ixyz=0;ixyz<3;ixyz++){
                  setVectorComponentGridData(sz,guide,d,dwk,id,jd,kd,ixyz,guide_out);
                  if(thin_out){
                    thinout_plot3d(d_thin,d,id_thin,jd_thin,kd_thin,id,jd,kd);
                    FP3DW.setFuncData(d_thin);
                  }else{
                    FP3DW.setFuncData(d);
                  }
                  if(!FP3DW.WriteFuncData()) printf("\terror WriteFuncData\n");
                }
              }else{ // double ---> double
                for(int ixyz=0;ixyz<3;ixyz++){
                  setVectorComponentGridData(sz,guide,dd,dwk,id,jd,kd,ixyz,guide_out);
                  if(thin_out){
                    thinout_plot3d(dd_thin,dd,id_thin,jd_thin,kd_thin,id,jd,kd);
                    FP3DW.setFuncData(dd_thin);
                  }else{
                    FP3DW.setFuncData(dd);
                  }
                  if(!FP3DW.WriteFuncData()) printf("\terror WriteFuncData\n");
                }
              }
            }
          }
          else // Fortranによる出力では出力項目すべてを一度に書き出し
          {
            if( d_type == SPH_FLOAT ){
              if( FP3DW.GetRealType() == OUTPUT_FLOAT ){ // float ---> float
                setVectorGridData(sz,guide,&d[outsize*ivar],wk,id,jd,kd,guide_out);
                if(thin_out){
                  thinout_plot3d(&d_thin[outsize_thin*ivar],&d[outsize*ivar],id_thin,jd_thin,kd_thin,id,jd,kd);
                }
              }else{ // float ---> double
                setVectorGridData(sz,guide,&dd[outsize*ivar],wk,id,jd,kd,guide_out);
                if(thin_out){
                  thinout_plot3d(&dd_thin[outsize_thin*ivar],&dd[outsize*ivar],id_thin,jd_thin,kd_thin,id,jd,kd);
                }
              }
            }else{
              if( FP3DW.GetRealType() == OUTPUT_FLOAT ){ // double ---> float
                setVectorGridData(sz,guide,&d[outsize*ivar],dwk,id,jd,kd,guide_out);
                if(thin_out){
                  thinout_plot3d(&d_thin[outsize_thin*ivar],&d[outsize*ivar],id_thin,jd_thin,kd_thin,id,jd,kd);
                }
              }else{ // double ---> double
                setVectorGridData(sz,guide,&dd[outsize*ivar],dwk,id,jd,kd,guide_out);
                if(thin_out){
                  thinout_plot3d(&dd_thin[outsize_thin*ivar],&dd[outsize*ivar],id_thin,jd_thin,kd_thin,id,jd,kd);
                }
              }
            }
          }
        }
        
        if(P3Op.IS_DivideFunc == ON){ // 分割出力のとき
          // write all
          if(FP3DW.GetFormat() != C_BINARY){ // C_BINARY以外は出力項目すべてを一度に書き出し
            if(thin_out){
              if( FP3DW.GetRealType() == OUTPUT_FLOAT ) FP3DW.setFuncData(d_thin);
              else                                      FP3DW.setFuncData(dd_thin);
            }
            else
            {
              if( FP3DW.GetRealType() == OUTPUT_FLOAT ) FP3DW.setFuncData(d);
              else                                      FP3DW.setFuncData(dd);
            }
            if(!FP3DW.WriteFuncData()) printf("\terror WriteFuncData\n");
            
          }
          // plot3d func file close
          FP3DW.CloseFile();
          
          //dfiファイルの出力
          std::string dfipre;
          dfipre = DI[i].Prefix + P3Op.basename;
          //dfipre = out_dirname + dfipre; // dfiファイルは他のdfiファイルと同じディレクトリにコピー
          if(inode == 0) if ( !DFI.Write_DFI_File(dfipre, (unsigned)m_step, (double)m_time, dfi_mng[var_Plot3D], true) ) Exit(0);
        }
        
        ivar=ivar+dim;
      }
      
      if(P3Op.IS_DivideFunc == OFF){ // 一括出力のとき
        // write all
        if(FP3DW.GetFormat() != C_BINARY){//C_BINARY以外は出力項目すべてを一度に書き出し
          if(thin_out){
            if( FP3DW.GetRealType() == OUTPUT_FLOAT ) FP3DW.setFuncData(d_thin);
            else                                      FP3DW.setFuncData(dd_thin);
          }
          else
          {
            if( FP3DW.GetRealType() == OUTPUT_FLOAT ) FP3DW.setFuncData(d);
            else                                      FP3DW.setFuncData(dd);
          }
          if(!FP3DW.WriteFuncData()) printf("\terror WriteFuncData\n");
        }
        // plot3d func file close
        FP3DW.CloseFile();
      }
      
      //dfiファイルの出力
      std::string dfipre;
      dfipre = P3Op.basename;
      //dfipre = out_dirname + P3Op.basename; // dfiファイルは他のdfiファイルと同じディレクトリにコピー
      if(inode == 0) if ( !DFI.Write_DFI_File(dfipre, (unsigned)m_step, (double)m_time, dfi_mng[var_Plot3D], true) ) Exit(0);
      
    }//node loop
    if(numProc > 1) ips=ips+nnode;
  }//step loop
  
  
  // deallocate work area
  if( d_type == SPH_FLOAT ) delete[] wk;
  else                      delete[] dwk;
  if( FP3DW.GetRealType() == OUTPUT_FLOAT ) delete[] d;
  else                                      delete[] dd;
  if(thin_out){
    if( FP3DW.GetRealType() == OUTPUT_FLOAT ) delete[] d_thin;
    else                                      delete[] dd_thin;
  }
  
  
  // output function name file
  if( myRank!=0 ) return; // myrank==0のときのみ処理
  LOG_OUT_ fprintf(fplog,"\t*** output func name file ***\n");
  STD_OUT_ printf("\t*** output func name file ***\n");
  
  // function_nameファイルはかならずformatted形式
  int keep_format=FP3DW.GetFormat();
  FP3DW.setFormat(FORMATTED);
  
  // open file
  if(P3Op.IS_DivideFunc == OFF){ // 一括出力のとき
    std::string ptmp;
    ptmp = DFI.Generate_FileName_Free(P3Op.basename, "nam", 0, 0, false);
    ptmp = out_dirname + ptmp;
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
      ptmp = out_dirname + ptmp;
      FP3DW.setFileName(ptmp.c_str());
      if(!FP3DW.OpenFile()){
        printf("\terror OpenFile : %s\n", ptmp.c_str());
        Exit(0);
      }
    }
    
    // write comp name
    if     ( !strcasecmp(comp.c_str(), "prs" ) )
    {
      FP3DW.WriteFunctionName("Pressure");
    }
    else if( !strcasecmp(comp.c_str(), "vel" ) )
    {
      FP3DW.WriteFunctionName("U-Velocity ; Velocity");
      FP3DW.WriteFunctionName("V-Velocity");
      FP3DW.WriteFunctionName("W-Velocity");
    }
    else if( !strcasecmp(comp.c_str(), "tmp" ) )
    {
      FP3DW.WriteFunctionName("Tempearture");
    }
    else if( !strcasecmp(comp.c_str(), "tp"  ) )
    {
      FP3DW.WriteFunctionName("Total_Pressure");
    }
    else if( !strcasecmp(comp.c_str(), "vrt" ) )
    {
      FP3DW.WriteFunctionName("U-Vorticity ; Vorticity");
      FP3DW.WriteFunctionName("V-Vorticity");
      FP3DW.WriteFunctionName("W-Vorticity");
    }
    else if( !strcasecmp(comp.c_str(), "i2vgt"  ) )
    {
      FP3DW.WriteFunctionName("2nd Invariant of Velocity Gradient Tensor");
    }
    else if( !strcasecmp(comp.c_str(), "hlt" ) )
    {
      FP3DW.WriteFunctionName("Helicity");
    }
    else
    {
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
      if( DI[i].dim == 1){
        FP3DW.WriteFunctionName(comp.c_str());
      }else if( DI[i].dim == 3){
        string compbuff;
        compbuff = "X-" + comp + " ; " + comp;
        FP3DW.WriteFunctionName(compbuff.c_str());
        compbuff = "Y-" + comp;
        FP3DW.WriteFunctionName(compbuff.c_str());
        compbuff = "Z-" + comp;
        FP3DW.WriteFunctionName(compbuff.c_str());
      }
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


void COMB::thinout_plot3d(float* dt, float* d, int idt, int jdt, int kdt, int id, int jd, int kd)
{
  int irest=(id-1)%thin_count;
  int jrest=(jd-1)%thin_count;
  int krest=(kd-1)%thin_count;
  int tc=thin_count;
  if(irest==0) irest=tc;
  if(jrest==0) jrest=tc;
  if(krest==0) krest=tc;
  
  int i,j,k,ix,jy,kz;
  size_t ipt,ip;
  
  for (k=0; k<kdt; k++){
    for (j=0; j<jdt; j++){
      for (i=0; i<idt; i++){
        ipt=k*idt*jdt+j*idt+i;
        ix=tc*i;
        jy=tc*j;
        kz=tc*k;
        if(i==idt-1) ix=tc*(i-1)+irest;
        if(j==jdt-1) jy=tc*(j-1)+jrest;
        if(k==kdt-1) kz=tc*(k-1)+krest;
        ip = _F_IDX_S3D(ix+1, jy+1, kz+1, id, jd, kd, 0);
        dt[ipt]=d[ip];
      }
    }
  }
}


void COMB::thinout_plot3d(double* dt, double* d, int idt, int jdt, int kdt, int id, int jd, int kd)
{
  int irest=(id-1)%thin_count;
  int jrest=(jd-1)%thin_count;
  int krest=(kd-1)%thin_count;
  int tc=thin_count;
  if(irest==0) irest=tc;
  if(jrest==0) jrest=tc;
  if(krest==0) krest=tc;
  
  int i,j,k,ix,jy,kz;
  size_t ipt,ip;
  
  for (k=0; k<kdt; k++){
    for (j=0; j<jdt; j++){
      for (i=0; i<idt; i++){
        ipt=k*idt*jdt+j*idt+i;
        ix=tc*i;
        jy=tc*j;
        kz=tc*k;
        if(i==idt-1) ix=tc*(i-1)+irest;
        if(j==jdt-1) jy=tc*(j-1)+jrest;
        if(k==kdt-1) kz=tc*(k-1)+krest;
        ip = _F_IDX_S3D(ix+1, jy+1, kz+1, id, jd, kd, 0);
        dt[ipt]=d[ip];
      }
    }
  }
}
