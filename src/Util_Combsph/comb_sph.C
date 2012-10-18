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
 * @file   comb_sph.C
 * @brief  COMB Class
 * @author kero
 */

#include "comb.h"

// #################################################################
// 
void COMB::output_sph()
{
  char tmp[FB_FILE_PATH_LENGTH];
  string prefix,outfile,infile;
  int fp_in =8;
  FILE *fp;
  FILE *fpin;

  int m_rank;
  int m_sv_type, m_d_type, d_type;
  int m_step, m_imax, m_jmax, m_kmax;
  long long m_dstep, m_dimax, m_djmax, m_dkmax;
  int sz[3];
  long long dsz[3];
  float m_org[3], m_pit[3];
  float m_time;
  double m_dorg[3], m_dpit[3];
  double m_dtime;
  int xsize,ysize,zsize,asize,vsize;
  long long dxsize,dysize,dzsize,dasize,dvsize;
  int dim;
  int wk1size,wk2size;
  float *wk1, *wk2;
  double *dwk1, *dwk2;
  int jsNode,jeNode,ksNode,keNode;

  //for thin_out option
  int m_imax_th,m_jmax_th,m_kmax_th;
  long long m_dimax_th,m_djmax_th,m_dkmax_th;


  //並列処理のためのインデックス作成
  int iproc=0;
  int ips=0;
  int ip;
  if(numProc > 1) {
    index.clear();
    for(int i=0;i<ndfi;i++){
      for(int j=0; j< DI[i].step.size(); j++ ) {
        ip=ips+j;
        if(myRank==iproc) index.push_back(ip);
        iproc++;
        if(iproc==numProc) iproc=0;
      }
      ips=ips+DI[i].step.size();
    }
    LOG_OUTV_ {
      fprintf(fplog,"\n");
      for(int j=0; j< index.size(); j++ ) {
        fprintf(fplog,"\tindex[%4d] = %d\n",j,index[j]);
      }
    }
    STD_OUTV_ {
      printf("\n");
      for(int j=0; j< index.size(); j++ ) {
        printf("\tindex[%4d] = %d\n",j,index[j]);
      }
    }
  }

  //dfi file loop ---> prs_, vel, ,,,
  if(numProc > 1) ips=0;
  for(int i=0;i<ndfi;i++){
    prefix=DI[i].Prefix;
    LOG_OUTV_ fprintf(fplog,"  COMBINE SPH START : %s\n", prefix.c_str());
    STD_OUTV_ printf("  COMBINE SPH START : %s\n", prefix.c_str());

    //Scalar or Vector
    if     ( !strcasecmp(prefix.c_str(), "prs_" ) ) dim=1;
    else if( !strcasecmp(prefix.c_str(), "vel_" ) ) dim=3;
    else if( !strcasecmp(prefix.c_str(), "tmp_" ) ) dim=1;
    else if( !strcasecmp(prefix.c_str(), "tp_"  ) ) dim=1;
    else if( !strcasecmp(prefix.c_str(), "vrt_" ) ) dim=3;
    else if( !strcasecmp(prefix.c_str(), "i2vgt_" ) ) dim=1;
    else if( !strcasecmp(prefix.c_str(), "hlt_" ) ) dim=1;
	else{
      m_step = DI[i].step[0];
      m_rank = DI[i].Node[0].RankID;
      infile = Generate_FileName(prefix, m_step, m_rank, true);
      infile = in_dirname + infile;
      ReadSphDataType (&m_sv_type, &m_d_type, fp_in, infile);
      dim=m_sv_type;
    }

    //step loop
    for(int j=0; j< DI[i].step.size(); j++ ) {

      //並列処理
      int iskip=1;
      if(numProc > 1) {
        iskip=0;
        ip=ips+j;
        for(int ic=0; ic< index.size(); ic++ ) {
          if(ip==index[ic]) iskip=1;
        }
      }
      if(!iskip) continue;

      m_step=DI[i].step[j];
      LOG_OUTV_ fprintf(fplog,"\tstep = %d\n", m_step);
      STD_OUTV_ printf("\tstep = %d\n", m_step);

      //連結出力ファイルオープン
      outfile = Generate_FileName(prefix, m_step, 0, false);
      outfile = out_dirname + outfile;
      memset(tmp, 0, sizeof(char)*FB_FILE_PATH_LENGTH);
      strcpy(tmp, outfile.c_str());
      if( (fp = fopen(outfile.c_str(), "wb")) == NULL ) {
        printf("\tCan't open file.(%s)\n",outfile.c_str());
        Exit(0);
      }

      //rank0のheader読み込み ---> 最初のファイルを開いてheaderだけ先に記述
      infile = Generate_FileName(prefix, m_step, 0, true);
      infile = in_dirname + infile;
      ReadSphDataType (&m_sv_type, &m_d_type, fp_in, infile);
      if( output_real_type == OUTPUT_REAL_UNKNOWN )
      {
        if( m_d_type == SPH_FLOAT ) output_real_type = OUTPUT_FLOAT;
        else                        output_real_type = OUTPUT_DOUBLE;
      }
      d_type=output_real_type;

      if( m_d_type == SPH_FLOAT ){
        ReadSphHeader (&m_step, &m_sv_type, &m_d_type, &m_imax, &m_jmax, &m_kmax, &m_time, m_org, m_pit, fp_in, infile);
      }else{
        ReadSphHeader (&m_dstep, &m_sv_type, &m_d_type, &m_dimax, &m_djmax, &m_dkmax, &m_dtime, m_dorg, m_dpit, fp_in, infile);
        m_step=(int)m_dstep;
        m_time=(float)m_dtime;
      }

      //全体サイズのキープ
      m_imax= DI[i].Global_Voxel[0] + 2*DI[i].GuideCell;
      m_jmax= DI[i].Global_Voxel[1] + 2*DI[i].GuideCell;
      m_kmax= DI[i].Global_Voxel[2] + 2*DI[i].GuideCell;

      //間引きを考慮
      m_imax_th=m_imax/thin_count;//間引き後のxサイズ
      m_jmax_th=m_jmax/thin_count;//間引き後のyサイズ
      m_kmax_th=m_kmax/thin_count;//間引き後のzサイズ
      if(m_imax%thin_count != 0) m_imax_th++;
      if(m_jmax%thin_count != 0) m_jmax_th++;
      if(m_kmax%thin_count != 0) m_kmax_th++;

      //連結ファイルのheader書き込み
      if( m_d_type == SPH_FLOAT && output_real_type == OUTPUT_FLOAT)
      { // float ---> float
        float out_pit[3];
        for(int ic=0;ic<3;ic++) out_pit[ic]=m_pit[ic]*float(thin_count);
        if( !(WriteSphHeader(m_step, m_sv_type, d_type, m_imax_th, m_jmax_th, m_kmax_th, m_time, m_org, out_pit, fp)) ) {
          printf("\twrite header error\n");
          Exit(0);
        }
      }
      else if( m_d_type == SPH_DOUBLE && output_real_type == OUTPUT_DOUBLE)
      { // double ---> double
        m_dimax_th=(long long)m_imax_th;
        m_djmax_th=(long long)m_jmax_th;
        m_dkmax_th=(long long)m_kmax_th;
        double out_dpit[3];
        for(int ic=0;ic<3;ic++) out_dpit[ic]=m_dpit[ic]*double(thin_count);
        if( !(WriteSphHeader(m_dstep, m_sv_type, d_type, m_dimax_th, m_djmax_th, m_dkmax_th, m_dtime, m_dorg, out_dpit, fp)) ) {
          printf("\twrite header error\n");
          Exit(0);
        }
      }
      if( m_d_type == SPH_FLOAT && output_real_type == OUTPUT_DOUBLE)
      { // float ---> double
        m_dimax_th=(long long)m_imax_th;
        m_djmax_th=(long long)m_jmax_th;
        m_dkmax_th=(long long)m_kmax_th;
        m_dstep=(long long)m_step;
        m_dtime=(double)m_time;
        m_dorg[0]=(double)m_org[0];
        m_dorg[1]=(double)m_org[1];
        m_dorg[2]=(double)m_org[2];
        m_dpit[0]=(double)m_pit[0];
        m_dpit[1]=(double)m_pit[1];
        m_dpit[2]=(double)m_pit[2];
        double out_dpit[3];
        for(int ic=0;ic<3;ic++) out_dpit[ic]=m_dpit[ic]*double(thin_count);
        if( !(WriteSphHeader(m_dstep, m_sv_type, d_type, m_dimax_th, m_djmax_th, m_dkmax_th, m_dtime, m_dorg, out_dpit, fp)) ) {
          printf("\twrite header error\n");
          Exit(0);
        }
      }
      if( m_d_type == SPH_DOUBLE && output_real_type == OUTPUT_FLOAT)
      { // double ---> float
        m_org[0]=(float)m_dorg[0];
        m_org[1]=(float)m_dorg[1];
        m_org[2]=(float)m_dorg[2];
        m_pit[0]=(float)m_dpit[0];
        m_pit[1]=(float)m_dpit[1];
        m_pit[2]=(float)m_dpit[2];
        float out_pit[3];
        for(int ic=0;ic<3;ic++) out_pit[ic]=m_pit[ic]*float(thin_count);
        if( !(WriteSphHeader(m_step, m_sv_type, d_type, m_imax_th, m_jmax_th, m_kmax_th, m_time, m_org, out_pit, fp)) ) {
          printf("\twrite header error\n");
          Exit(0);
        }
      }

      //全体の大きさの計算とデータのヘッダ書き込み
      int dummy;
      size_t dLen;
      //dLen = (size_t)(m_imax * m_jmax * m_kmax);
      dLen = (size_t)(m_imax_th * m_jmax_th * m_kmax_th);
      if( m_sv_type == SPH_VECTOR ) dLen *= 3;
      //if( m_d_type == SPH_FLOAT ) dummy = dLen * sizeof(float);
      if( output_real_type == OUTPUT_FLOAT ) dummy = dLen * sizeof(float);
      else                                   dummy = dLen * sizeof(double);
      if( !(WriteCombineDataMarker(dummy, fp)) ) {
        printf("\twrite data header error\n");
        Exit(0);
      }

      //同一Z面ループ
      int zstart=0;
      for(int kp=0; kp< DI[i].index_z.size()-1; kp++ ) {
        ksNode = DI[i].index_z[kp];
        keNode = DI[i].index_z[kp+1];

        //書き込みworkareaのサイズ決め
        m_imax_th=m_imax/thin_count;//間引き後のxサイズ
        m_jmax_th=m_jmax/thin_count;//間引き後のyサイズ
        if(m_imax%thin_count != 0) m_imax_th++;
        if(m_jmax%thin_count != 0) m_jmax_th++;
        //xsize=m_imax;
        //ysize=m_jmax;
        xsize=m_imax_th;
        ysize=m_jmax_th;
        zsize=DI[i].Node[ksNode].VoxelSize[2]; 
        asize=xsize*ysize;
        vsize=0;
        for(int n=ksNode; n< keNode; n++ ) {
          int szx,szy,szz;
          szx=DI[i].Node[n].VoxelSize[0];
          szy=DI[i].Node[n].VoxelSize[1];
          szz=DI[i].Node[n].VoxelSize[2];
          int vdum=szx*szy*szz;
          if(vsize < vdum) vsize=vdum;
        }

        // メモリチェック
        LOG_OUTV_ fprintf(fplog,"\tNode %4d - Node %4d\n", ksNode, keNode-1);
        STD_OUTV_ printf("\tNode %4d - Node %4d\n", ksNode, keNode-1);
        double mc1 = (double)asize*(double)dim;
        double mc2 = (double)vsize*(double)dim;
        if(mc1>(double)INT_MAX){// 整数値あふれ出しチェック //参考 894*894*894*3=2143550952 INT_MAX 2147483647
          printf("\tsize error : mc1>INT_MAX\n");
          Exit(0);
        }
        if(mc2>(double)INT_MAX){// 整数値あふれ出しチェック //参考 894*894*894*3=2143550952 INT_MAX 2147483647
          printf("\tsize error : mc2>INT_MAX\n");
          Exit(0);
        }
        double TotalMemory=0.0; // = mc * (double)sizeof(REAL_TYPE);
        if( output_real_type == OUTPUT_FLOAT ) TotalMemory = TotalMemory + mc1 * (double)sizeof(float);
        else                                   TotalMemory = TotalMemory + mc1 * (double)sizeof(double);
        if( m_d_type == SPH_FLOAT ) TotalMemory = TotalMemory + mc2 * (double)sizeof(float);
        else                        TotalMemory = TotalMemory + mc2 * (double)sizeof(double);
        LOG_OUT_ MemoryRequirement(TotalMemory,fplog);
        STD_OUT_ MemoryRequirement(TotalMemory,stdout);

        //allocate work area
        wk1size=asize*dim;
        wk2size=vsize*dim;

        // work for read
        if( m_d_type == SPH_FLOAT ){
          if ( !(wk2 = new float[wk2size]) ){
            printf("\tallocate error : wk2\n");
            Exit(0);
          }
        }else{
          if ( !(dwk2 = new double[wk2size]) ){
            printf("\tallocate error : dwk2\n");
            Exit(0);
          }
        }

        // work for write
        if( output_real_type == OUTPUT_FLOAT ){
          if ( !(wk1 = new float[wk1size]) ){
            printf("\tallocate error : wk1\n");
            Exit(0);
          }
        }else if( output_real_type == OUTPUT_DOUBLE ){
          if ( !(dwk1 = new double[wk1size]) ){
            printf("\tallocate error : dwk1\n");
            Exit(0);
          }
        }

        //同一Y面のindexのスタート、エンドを決める
        int index_y_s=0;
        int index_y_e=0;
        for(int jp=0; jp< DI[i].index_y.size()-1; jp++ ) {
          jsNode = DI[i].index_y[jp];
          jeNode = DI[i].index_y[jp+1];
          if(jsNode == ksNode) index_y_s=jp;
          if(jeNode == keNode) index_y_e=jp+1;
        }

        //同一Y面ループをz方向高さ分だけ繰り返しファイルをcombine
        int zzz=0;
        while(zzz != zsize){
          int zh = zstart + zzz;
          int zrest=zh%thin_count;

          if(zrest == 0){

            int ystart=0;
            for(int jp=index_y_s; jp< index_y_e; jp++ ) {
              jsNode = DI[i].index_y[jp];
              jeNode = DI[i].index_y[jp+1];
              LOG_OUTV_ fprintf(fplog,"\t  Node %4d - Node %4d\n", jsNode, jeNode-1);
              STD_OUTV_ printf("\t  Node %4d - Node %4d\n", jsNode, jeNode-1);
           
              int xstart=0;
              for(int n=jsNode; n< jeNode; n++ ) {
           
                //連結対象ファイルオープン
                m_rank=DI[i].Node[n].RankID;
                infile = Generate_FileName(prefix, m_step, m_rank, true);
                infile = in_dirname + infile;

                //combine
                sz[0]=DI[i].Node[n].VoxelSize[0];
                sz[1]=DI[i].Node[n].VoxelSize[1];
                sz[2]=DI[i].Node[n].VoxelSize[2];
                if( m_d_type == SPH_FLOAT ){
                  ReadSphData ( wk2, wk2size, sz, dim, fp_in, infile );
                  if( output_real_type == OUTPUT_FLOAT ){
                    CombineLayerData(
                      wk1, wk1size, wk2, wk2size, xsize, ysize, zsize,
                      zzz, sz, dim, xstart, ystart );
                  }else if( output_real_type == OUTPUT_DOUBLE ){
                    CombineLayerData(
                      dwk1, wk1size, wk2, wk2size, xsize, ysize, zsize,
                      zzz, sz, dim, xstart, ystart );
                  }
                }else{
                  ReadSphData ( dwk2, wk2size, sz, dim, fp_in, infile );
                  if( output_real_type == OUTPUT_FLOAT ){
                    CombineLayerData(
                      wk1, wk1size, dwk2, wk2size, xsize, ysize, zsize,
                      zzz, sz, dim, xstart, ystart );
                  }else if( output_real_type == OUTPUT_DOUBLE ){
                    CombineLayerData(
                      dwk1, wk1size, dwk2, wk2size, xsize, ysize, zsize,
                      zzz, sz, dim, xstart, ystart );
                  }
                }
           
                //書き込みデータスタート位置の更新
                xstart=xstart+DI[i].Node[n].VoxelSize[0];
           
              }
           
              //書き込みデータスタート位置の更新
              ystart=ystart+DI[i].Node[jsNode].VoxelSize[1];
            }
           
            //write
            if( output_real_type == OUTPUT_FLOAT ){
              if( !(WriteCombineData(wk1, (size_t)wk1size, fp)) ) {
                printf("\twrite data footer error\n");
                Exit(0);
              }
            }else{
              if( !(WriteCombineData(dwk1, (size_t)wk1size, fp)) ) {
                printf("\twrite data footer error\n");
                Exit(0);
              }
            }

          }

          zzz++;
        }

        //deallocate workarea
        if( m_d_type == SPH_FLOAT ) delete[] wk2;
        else                        delete[] dwk2;
        if( output_real_type == OUTPUT_FLOAT ) delete[] wk1;
        else                                   delete[] dwk1;

        //書き込みデータスタート位置の更新
        zstart=zstart+DI[i].Node[jsNode].VoxelSize[2];

	  }//z --- kp

      //平均値識別の記述 ---> なし？

      //データのフッタ書き込み
      if( !(WriteCombineDataMarker(dummy, fp)) ) {
        printf("\twrite data error\n");
        Exit(0);
      }

      //出力ファイルクローズ
      fclose(fp);
    }
    if(numProc > 1) ips=ips+DI[i].step.size();
  }

}


#if 0

// #################################################################
// 
bool COMB::WriteSphHeader(
    int step, int sv_type, int d_type, int imax, int jmax, int kmax,
    REAL_TYPE time, REAL_TYPE* org, REAL_TYPE* pit, FILE *fp)
{
  if( !fp ) return false;

  unsigned int dmy;
  dmy = 8;
  if( fwrite(&dmy, sizeof(int), 1, fp) != 1 ) return false;
  if( fwrite(&sv_type, sizeof(int), 1, fp) != 1 ) return false;
  if( fwrite(&d_type, sizeof(int), 1, fp) != 1 ) return false;
  if( fwrite(&dmy, sizeof(int), 1, fp) != 1 ) return false;

  if( d_type == OUTPUT_FLOAT ){
    dmy = 12;
    if( fwrite(&dmy, sizeof(int), 1, fp) != 1 ) return false;
    if( fwrite(&imax, sizeof(int), 1, fp) != 1 ) return false;
    if( fwrite(&jmax, sizeof(int), 1, fp) != 1 ) return false;
    if( fwrite(&kmax, sizeof(int), 1, fp) != 1 ) return false;
    if( fwrite(&dmy, sizeof(int), 1, fp) != 1 ) return false;

    if( fwrite(&dmy, sizeof(int), 1, fp) != 1 ) return false;
    if( fwrite(org, sizeof(float), 3, fp) != 3 ) return false;
    if( fwrite(&dmy, sizeof(int), 1, fp) != 1 ) return false;

    if( fwrite(&dmy, sizeof(int), 1, fp) != 1 ) return false;
    if( fwrite(pit, sizeof(float), 3, fp) != 3 ) return false;
    if( fwrite(&dmy, sizeof(int), 1, fp) != 1 ) return false;
  } else { // OUTPUT_DOUBLE
    dmy = 24;
    if( fwrite(&dmy, sizeof(int), 1, fp) != 1 ) return false;
    if( fwrite(&imax, sizeof(long long), 1, fp) != 1 ) return false;
    if( fwrite(&jmax, sizeof(long long), 1, fp) != 1 ) return false;
    if( fwrite(&kmax, sizeof(long long), 1, fp) != 1 ) return false;
    if( fwrite(&dmy, sizeof(int), 1, fp) != 1 ) return false;

    if( fwrite(&dmy, sizeof(int), 1, fp) != 1 ) return false;
    if( fwrite(org, sizeof(double), 3, fp) != 3 ) return false;
    if( fwrite(&dmy, sizeof(int), 1, fp) != 1 ) return false;

    if( fwrite(&dmy, sizeof(int), 1, fp) != 1 ) return false;
    if( fwrite(pit, sizeof(double), 3, fp) != 3 ) return false;
    if( fwrite(&dmy, sizeof(int), 1, fp) != 1 ) return false;
  }

  if( d_type == OUTPUT_FLOAT ) dmy = 8;
  else dmy = 16;
  if( fwrite(&dmy, sizeof(int), 1, fp) != 1 ) return false;
  if( d_type == OUTPUT_FLOAT ){
    if( fwrite(&step, sizeof(int), 1, fp) != 1 ) return false;
    if( fwrite(&time, sizeof(float), 1, fp) != 1 ) return false;
  }else{
    if( fwrite(&step, sizeof(long long), 1, fp) != 1 ) return false;
    if( fwrite(&time, sizeof(double), 1, fp) != 1 ) return false;
  }
  if( fwrite(&dmy, sizeof(int), 1, fp) != 1 ) return false;

  return true;

}

#endif


// #################################################################
// 
bool COMB::WriteSphHeader(
    int step, int sv_type, int d_type, int imax, int jmax, int kmax,
    float time, float* org, float* pit, FILE *fp)
{
  if( !fp ) return false;

  unsigned int dmy;

  dmy = 2 * sizeof(int);
  if( fwrite(&dmy, sizeof(int), 1, fp) != 1 ) return false;
  if( fwrite(&sv_type, sizeof(int), 1, fp) != 1 ) return false;
  if( fwrite(&d_type, sizeof(int), 1, fp) != 1 ) return false;
  if( fwrite(&dmy, sizeof(int), 1, fp) != 1 ) return false;

  dmy = 3 * sizeof(int);
  if( fwrite(&dmy, sizeof(int), 1, fp) != 1 ) return false;
  if( fwrite(&imax, sizeof(int), 1, fp) != 1 ) return false;
  if( fwrite(&jmax, sizeof(int), 1, fp) != 1 ) return false;
  if( fwrite(&kmax, sizeof(int), 1, fp) != 1 ) return false;
  if( fwrite(&dmy, sizeof(int), 1, fp) != 1 ) return false;

  dmy = 3 * sizeof(float);
  if( fwrite(&dmy, sizeof(int), 1, fp) != 1 ) return false;
  if( fwrite(org, sizeof(float), 3, fp) != 3 ) return false;
  if( fwrite(&dmy, sizeof(int), 1, fp) != 1 ) return false;

  if( fwrite(&dmy, sizeof(int), 1, fp) != 1 ) return false;
  if( fwrite(pit, sizeof(float), 3, fp) != 3 ) return false;
  if( fwrite(&dmy, sizeof(int), 1, fp) != 1 ) return false;

  dmy  = sizeof(int) + sizeof(float);
  if( fwrite(&dmy, sizeof(int), 1, fp) != 1 ) return false;
  if( fwrite(&step, sizeof(int), 1, fp) != 1 ) return false;
  if( fwrite(&time, sizeof(float), 1, fp) != 1 ) return false;
  if( fwrite(&dmy, sizeof(int), 1, fp) != 1 ) return false;

  return true;

}

// #################################################################
// 
bool COMB::WriteSphHeader(
    int step, int sv_type, int d_type, int imax, int jmax, int kmax,
    double time, double* org, double* pit, FILE *fp)
{
  if( !fp ) return false;

  unsigned int dmy;
  dmy = 2 * sizeof(int);
  if( fwrite(&dmy, sizeof(int), 1, fp) != 1 ) return false;
  if( fwrite(&sv_type, sizeof(int), 1, fp) != 1 ) return false;
  if( fwrite(&d_type, sizeof(int), 1, fp) != 1 ) return false;
  if( fwrite(&dmy, sizeof(int), 1, fp) != 1 ) return false;

  dmy = 3 * sizeof(long long);
  if( fwrite(&dmy, sizeof(int), 1, fp) != 1 ) return false;
  if( fwrite(&imax, sizeof(long long), 1, fp) != 1 ) return false;
  if( fwrite(&jmax, sizeof(long long), 1, fp) != 1 ) return false;
  if( fwrite(&kmax, sizeof(long long), 1, fp) != 1 ) return false;
  if( fwrite(&dmy, sizeof(int), 1, fp) != 1 ) return false;

  dmy = 3 * sizeof(double);
  if( fwrite(&dmy, sizeof(int), 1, fp) != 1 ) return false;
  if( fwrite(org, sizeof(double), 3, fp) != 3 ) return false;
  if( fwrite(&dmy, sizeof(int), 1, fp) != 1 ) return false;

  if( fwrite(&dmy, sizeof(int), 1, fp) != 1 ) return false;
  if( fwrite(pit, sizeof(double), 3, fp) != 3 ) return false;
  if( fwrite(&dmy, sizeof(int), 1, fp) != 1 ) return false;

  dmy = sizeof(long long) + sizeof(double);
  if( fwrite(&dmy, sizeof(int), 1, fp) != 1 ) return false;
  if( fwrite(&step, sizeof(long long), 1, fp) != 1 ) return false;
  if( fwrite(&time, sizeof(double), 1, fp) != 1 ) return false;
  if( fwrite(&dmy, sizeof(int), 1, fp) != 1 ) return false;

  return true;
}


// #################################################################
// 
bool COMB::WriteCombineDataMarker(int dmy, FILE* fp)
{
  if( fwrite(&dmy, sizeof(int), 1, fp) != 1 ) return false;
  return true;
}

#if 0

// #################################################################
// 
bool COMB::WriteCombineData(REAL_TYPE* data, size_t dLen, int d_type, FILE* fp)
{
  if( d_type == OUTPUT_FLOAT ){
    if( fwrite(data, sizeof(float), dLen, fp) != dLen ) return false;
  }else{
    if( fwrite(data, sizeof(double), dLen, fp) != dLen ) return false;
  }
  return true;
}

#endif

// #################################################################
// 
bool COMB::WriteCombineData(float* data, size_t dLen, FILE* fp)
{
  if( fwrite(data, sizeof(float), dLen, fp) != dLen ) return false;
  return true;
}

// #################################################################
// 
bool COMB::WriteCombineData(double* data, size_t dLen, FILE* fp)
{
  if( fwrite(data, sizeof(double), dLen, fp) != dLen ) return false;
  return true;
}

