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
 * @file   comb_sph.C
 * @brief  COMB Class
 * @author aics
 */

#include "comb.h"

// #################################################################
//
void COMB::output_sph()
{
  char tmp[FILE_PATH_LENGTH];
//CDM.20131008.s
  string prefix,outfile,infile,inPath;
//CDM.20131008.e
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
  
  //for thin_out option
  int m_imax_th,m_jmax_th,m_kmax_th;
  long long m_dimax_th,m_djmax_th,m_dkmax_th;
  
  // 出力モード
  bool mio = false;
  const cdm_MPI* DFI_MPI = dfi[0]->GetcdmMPI();
  if( DFI_MPI->NumberOfRank > 1) mio=true;
  
  //並列処理のためのインデックス作成
  int iproc=0;
  int ips=0;
  int ip;
  if(numProc > 1) {
    index.clear();
    for(int i=0;i<ndfi;i++){
      const cdm_TimeSlice* TSlice = dfi[i]->GetcdmTimeSlice();
      for(int j=0; j< TSlice->SliceList.size(); j++ ) {
        ip=ips+j;
        if(myRank==iproc) index.push_back(ip);
        iproc++;
        if(iproc==numProc) iproc=0;
      }
      ips=ips+TSlice->SliceList.size();
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
  if (numProc > 1) ips=0;
  
  for (int i=0;i<ndfi;i++) {

//CDM.20131008.s
    inPath = CDM::cdmPath_DirName(dfi_name[i]);
//CDM.20131008.e

    const cdm_FileInfo* DFI_FInfo = dfi[i]->GetcdmFileInfo();
    prefix=DFI_FInfo->Prefix;
    LOG_OUTV_ fprintf(fplog,"  COMBINE SPH START : %s\n", prefix.c_str());
    STD_OUTV_ printf("  COMBINE SPH START : %s\n", prefix.c_str());
    
    //Scalar or Vector
    dim=DFI_FInfo->Component;
    
    const cdm_Domain* DFI_Domian = dfi[i]->GetcdmDomain();
    
    cdm_Process* DFI_Process = (cdm_Process *)dfi[i]->GetcdmProcess();
    
    int div[3];
    div[0]=DFI_Domian->GlobalDivision[0];
    div[1]=DFI_Domian->GlobalDivision[1];
    div[2]=DFI_Domian->GlobalDivision[2];
    
    int *rankMap;
    typedef std::map<int,int> headT;
    headT mapHeadX;
    headT mapHeadY;
    headT mapHeadZ;
    rankMap = DFI_Process->CreateRankMap(div,
                                         mapHeadX,
                                         mapHeadY,
                                         mapHeadZ);
    
    const cdm_TimeSlice* TSlice = dfi[i]->GetcdmTimeSlice();
    for(int j=0; j< TSlice->SliceList.size(); j++ ) {
      
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
      
      m_step=TSlice->SliceList[j].step;
      m_time=(float)TSlice->SliceList[j].time;
      
      LOG_OUTV_ fprintf(fplog,"\tstep = %d\n", m_step);
      STD_OUTV_ printf("\tstep = %d\n", m_step);
      
      //連結出力ファイルオープン
      outfile = Generate_FileName(prefix, m_step, 0, false);
      outfile = out_dirname + outfile;
      memset(tmp, 0, sizeof(char)*FILE_PATH_LENGTH);
      strcpy(tmp, outfile.c_str());
      if( (fp = fopen(outfile.c_str(), "wb")) == NULL ) {
        printf("\tCan't open file.(%s)\n",outfile.c_str());
        Exit(0);
      }
      
      //rank0のheader読み込み ---> 最初のファイルを開いてheaderだけ先に記述
//CDM.20131008.s
      //infile = dfi[i]->Generate_FieldFileName(0,m_step,mio);
      infile = CDM::cdmPath_ConnectPath(inPath,dfi[i]->Generate_FieldFileName(0,m_step,mio));
//CDM.20131008.e
      
      EMatchType eType;
      eType = isMatchEndian(infile, 8);
      bool matchEndian = true;
      if( eType == UnMatch ) matchEndian = false;
      
      //m_sv_typeのセット (スカラー or ベクター)
      if( dfi[i]->GetNumComponent() == 1 ) {
        m_sv_type = SPH_SCALAR;
      } else if( dfi[i]->GetNumComponent() >= 3 ) {
        m_sv_type = SPH_VECTOR;
      } else Exit(0);
      
      //m_d_typeのセット (float or double)
      if( dfi[i]->GetDataType() == CDM::E_CDM_FLOAT32 ) {
        m_d_type = SPH_FLOAT;
      } else if( dfi[i]->GetDataType() == CDM::E_CDM_FLOAT64 ) {
        m_d_type = SPH_DOUBLE;
      } else Exit(0);
      
      if( output_real_type == OUTPUT_REAL_UNKNOWN )
      {
        if( m_d_type == SPH_FLOAT ) output_real_type = OUTPUT_FLOAT;
        else                        output_real_type = OUTPUT_DOUBLE;
      }
      d_type=output_real_type;
      
      //sphのヘッダーレコードからオリジンとピッチを読込む
      ReadSphHeader(m_dorg,  m_dpit, infile);

//CDM.20131008.s
      //ガイドセルがあるときオリジンを修正
      if( DFI_FInfo->GuideCell>0 ) {
        m_dorg[0]=DFI_Domian->GlobalOrigin[0]+0.5*m_dpit[0];
        m_dorg[1]=DFI_Domian->GlobalOrigin[1]+0.5*m_dpit[1];
        m_dorg[2]=DFI_Domian->GlobalOrigin[2]+0.5*m_dpit[2];
      }
      
      //全体サイズのキープ
      //m_imax= DFI_Domian->GlobalVoxel[0] + 2*DFI_FInfo->GuideCell;
      //m_jmax= DFI_Domian->GlobalVoxel[1] + 2*DFI_FInfo->GuideCell;
      //m_kmax= DFI_Domian->GlobalVoxel[2] + 2*DFI_FInfo->GuideCell;
      m_imax= DFI_Domian->GlobalVoxel[0];
      m_jmax= DFI_Domian->GlobalVoxel[1];
      m_kmax= DFI_Domian->GlobalVoxel[2];
//CDM.20131008.e
      
      //間引きを考慮
      m_imax_th=m_imax/thin_count;//間引き後のxサイズ
      m_jmax_th=m_jmax/thin_count;//間引き後のyサイズ
      m_kmax_th=m_kmax/thin_count;//間引き後のzサイズ
      if(m_imax%thin_count != 0) m_imax_th++;
      if(m_jmax%thin_count != 0) m_jmax_th++;
      if(m_kmax%thin_count != 0) m_kmax_th++;
      
      //出力sphのヘッダーレコードを出力
      double out_dpit[3];
      for(int ic=0;ic<3;ic++) out_dpit[ic]=m_dpit[ic]*double(thin_count);
      if( !(WriteSphHeader(m_step, m_sv_type, d_type, m_imax_th, m_jmax_th, m_kmax_th,
                           m_time, m_dorg, out_dpit, fp)) ) {
        printf("\twrite header error\n");
        Exit(0);
      }
      
      //全体の大きさの計算とデータのヘッダ書き込み
      int dummy;
      size_t dLen;

      dLen = size_t(m_imax_th) * size_t(m_jmax_th) * size_t(m_kmax_th);
      if( m_sv_type == SPH_VECTOR ) dLen *= 3;
      if( output_real_type == OUTPUT_FLOAT ) dummy = dLen * sizeof(float);
      else                                   dummy = dLen * sizeof(double);
      if( !(WriteCombineDataMarker(dummy, fp)) ) {
        printf("\twrite data header error\n");
        Exit(0);
      }
      
      //書き込みworkareaのサイズ決め
      m_imax_th=m_imax/thin_count;//間引き後のxサイズ
      m_jmax_th=m_jmax/thin_count;//間引き後のyサイズ
      if(m_imax%thin_count != 0) m_imax_th++;
      if(m_jmax%thin_count != 0) m_jmax_th++;
      xsize=m_imax_th;
      ysize=m_jmax_th;
      asize=xsize*ysize;
      vsize=0;
      for(int n=0; n< DFI_Process->RankList.size(); n++ ) {
        int szx,szy,szz;
        szx=DFI_Process->RankList[n].VoxelSize[0];
        szy=DFI_Process->RankList[n].VoxelSize[1];
        szz=DFI_Process->RankList[n].VoxelSize[2];
        int vdum=szx*szy*szz;
        if(vsize < vdum) vsize=vdum;
      }
      // メモリチェック
//CDM.20131008.s
      /*
      LOG_OUTV_ fprintf(fplog,"\tNode %4d - Node %4d\n", 0,
                        DFI_Domian->GlobalVoxel[2] + 2*DFI_FInfo->GuideCell-1);
      STD_OUTV_ printf("\tNode %4d - Node %4d\n", 0,
                       DFI_Domian->GlobalVoxel[2] + 2*DFI_FInfo->GuideCell-1);
      */
      LOG_OUTV_ fprintf(fplog,"\tNode %4d - Node %4d\n", 0,
                        DFI_Domian->GlobalVoxel[2] -1);
      STD_OUTV_ printf("\tNode %4d - Node %4d\n", 0,
                       DFI_Domian->GlobalVoxel[2] -1);
//CDM.20131008.e
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
      
      //出力バッファのインスタンス
      int szS[3];
      int headS[3],tailS[3];
      szS[0]=m_imax_th;
      szS[1]=m_jmax_th;
      szS[2]=1;
//CDM.20131008.s
      //headS[0]=0-DFI_FInfo->GuideCell;
      //headS[1]=0-DFI_FInfo->GuideCell;
      headS[0]=0;
      headS[1]=0;
//CDM.20131008.e
      tailS[0]=headS[0]+m_imax-1;
      tailS[1]=headS[1]+m_jmax-1;
      
      CDM::E_CDM_DTYPE cdm_d_type;
      if( output_real_type == OUTPUT_FLOAT ) {
        cdm_d_type = CDM::E_CDM_FLOAT32;
      } else if( output_real_type == OUTPUT_DOUBLE ) {
        cdm_d_type = CDM::E_CDM_FLOAT64;
      } else Exit(0);
      cdm_Array* src = cdm_Array::instanceArray
      ( cdm_d_type
       , DFI_FInfo->ArrayShape
       , szS
       , 0
       , DFI_FInfo->Component );
      
      int kdiv,jdiv,idiv;
      //z方向の分割数回のループ
      for( headT::iterator itz=mapHeadZ.begin(); itz!= mapHeadZ.end(); itz++ ) {
        kdiv = itz->second;
        int kp_sta,kp_end;
        kp_sta = itz->first;

        int nrank = _CDM_IDX_IJK(0,0,kdiv,div[0],div[1],div[2],0);
        kp_end = kp_sta + DFI_Process->RankList[nrank].VoxelSize[2];
       
//CDM.20131008.s 
        //if( kdiv == 0 ) kp_sta -= DFI_FInfo->GuideCell;

        //if( kp_end == DFI_Domian->GlobalVoxel[2]+1 ) kp_end += DFI_FInfo->GuideCell;
//CDM.20131008.e 
        
        //同一Z面のループ
        for(int kp=kp_sta; kp< kp_end; kp++) {
//CDM.20131008.s
          //int kk = kp+DFI_FInfo->GuideCell-1;
          int kk = kp-1;
//CDM.20131008.e
          //間引きの層のときスキップ
          if( kk%thin_count != 0 ) continue;
          
          //y方向の分割数のループ
          for( headT::iterator ity=mapHeadY.begin(); ity!= mapHeadY.end(); ity++ ) {
            jdiv = ity->second;
            int jp_sta,jp_end;
            jp_sta = ity->first;
            int nrank = _CDM_IDX_IJK(0,jdiv,kdiv,div[0],div[1],div[2],0);
            jp_end = jp_sta + DFI_Process->RankList[nrank].VoxelSize[1];

//CDM.20131008.s            
            //if( jdiv == 0 ) jp_sta -= DFI_FInfo->GuideCell;

            //if( jp_end == DFI_Domian->GlobalVoxel[1]+1 ) jp_end += DFI_FInfo->GuideCell;
//CDM.20131008.e            
            
            //x方向の分割数のループ
            for( headT::iterator itx=mapHeadX.begin(); itx!= mapHeadX.end(); itx++ ) {
              idiv = itx->second;
              int ip_sta,ip_end;
              ip_sta = itx->first;
              int nrank = _CDM_IDX_IJK(idiv,jdiv,kdiv,div[0],div[1],div[2],0);
              ip_end = ip_sta + DFI_Process->RankList[nrank].VoxelSize[0];
             
//CDM.20131008.s 
              //if( idiv == 0 ) ip_sta -= DFI_FInfo->GuideCell;

              //if( ip_end == DFI_Domian->GlobalVoxel[0]+1 ) ip_end += DFI_FInfo->GuideCell;
//CDM.20131008.e 
              
              int RankID = _CDM_IDX_IJK(idiv,jdiv,kdiv,div[0],div[1],div[2],0);
              
              int read_sta[3],read_end[3];
              read_sta[0]=ip_sta;
              read_sta[1]=jp_sta;
              read_sta[2]=kp;
              read_end[0]=ip_end-1;
              read_end[1]=jp_end-1;
              read_end[2]=kp;
              
              m_rank=DFI_Process->RankList[RankID].RankID;
//CDM.20131008.s
              //infile = dfi[i]->Generate_FieldFileName(m_rank,m_step,mio);
              infile = CDM::cdmPath_ConnectPath(inPath,dfi[i]->Generate_FieldFileName(m_rank,m_step,mio));
//CDM.20131008.e
              unsigned int avr_step;
              double avr_time;
              CDM::E_CDM_ERRORCODE ret;
              //連結対象ファイルの読込み
              cdm_Array* buf = dfi[i]->ReadFieldData(infile, m_step, m_dtime,
                                                     read_sta, read_end,
                                                     DFI_Process->RankList[RankID].HeadIndex,
                                                     DFI_Process->RankList[RankID].TailIndex,
                                                     true, avr_step, avr_time, ret);
              
              //headIndexを０スタートにしてセット
              int headB[3];
//CDM.20131008.s
              //headB[0]=read_sta[0]+DFI_FInfo->GuideCell-1;
              //headB[1]=read_sta[1]+DFI_FInfo->GuideCell-1;
              //headB[2]=read_sta[2]+DFI_FInfo->GuideCell-1;
              headB[0]=read_sta[0]-1;
              headB[1]=read_sta[1]-1;
              headB[2]=read_sta[2]-1;
//CDM.20131008.e
              buf->setHeadIndex( headB );
              
              int headS0[3];
//CDM.20131008.s
              //headS0[0]=headS[0]+DFI_FInfo->GuideCell;
              //headS0[1]=headS[1]+DFI_FInfo->GuideCell;
              headS0[0]=headS[0];
              headS0[1]=headS[1];
//CDM.20131008.s
              headS0[2]=kk;
              src->setHeadIndex( headS0 );
              
              tailS[2]=headS0[2];
              
              
              //出力配列へのコンバイン
              combineXY(matchEndian,buf,src,headS0,tailS);
              
              delete buf;
              
            } //x --- itx
          } //y --- ity
          
          //一層分出力
          size_t dLen = szS[0]*szS[1]*szS[2]*DFI_FInfo->Component;
          if( src->writeBinary(fp) != dLen ) Exit(0);
          
        } //z --- kp
      } //z --- itz
      
      delete src;
      
      //平均値識別の記述 ---> なし？
      
      //データのフッタ書き込み
      if( !(WriteCombineDataMarker(dummy, fp)) ) {
        printf("\twrite data error\n");
        Exit(0);
      }
      
      //出力ファイルクローズ
      fclose(fp);
    }
    if(numProc > 1) ips=ips+TSlice->SliceList.size();
  }
  
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
  
  if( d_type == SPH_FLOAT ) {
    dmy = 3 * sizeof(int);
  } else if( d_type == SPH_DOUBLE ) {
    dmy = 3 * sizeof(long long);
  }
  
  //dmy = 3 * sizeof(long long);
  if( fwrite(&dmy, sizeof(int), 1, fp) != 1 ) return false;
  if( d_type == SPH_FLOAT ) {
    if( fwrite(&imax, sizeof(int), 1, fp) != 1 ) return false;
    if( fwrite(&jmax, sizeof(int), 1, fp) != 1 ) return false;
    if( fwrite(&kmax, sizeof(int), 1, fp) != 1 ) return false;
  } else {
    if( fwrite(&imax, sizeof(long long), 1, fp) != 1 ) return false;
    if( fwrite(&jmax, sizeof(long long), 1, fp) != 1 ) return false;
    if( fwrite(&kmax, sizeof(long long), 1, fp) != 1 ) return false;
  }
  if( fwrite(&dmy, sizeof(int), 1, fp) != 1 ) return false;
  
  if( d_type == SPH_FLOAT ) {
    dmy = 3 * sizeof(float);
  } else {
    dmy = 3 * sizeof(double);
  }
  //dmy = 3 * sizeof(double);
  if( fwrite(&dmy, sizeof(int), 1, fp) != 1 ) return false;
  if( d_type == SPH_FLOAT ) {
    float tmp[3];
    tmp[0] = (float)org[0];
    tmp[1] = (float)org[1];
    tmp[2] = (float)org[2];
    if( fwrite(tmp, sizeof(float), 3, fp) != 3 ) return false;
  } else {
    if( fwrite(org, sizeof(double), 3, fp) != 3 ) return false;
  }
  if( fwrite(&dmy, sizeof(int), 1, fp) != 1 ) return false;
  
  if( fwrite(&dmy, sizeof(int), 1, fp) != 1 ) return false;
  if( d_type == SPH_FLOAT ) {
    float tmp[3];
    tmp[0] = (float)pit[0];
    tmp[1] = (float)pit[1];
    tmp[2] = (float)pit[2];
    if( fwrite(tmp, sizeof(float), 3, fp) != 3 ) return false;
  } else {
    if( fwrite(pit, sizeof(double), 3, fp) != 3 ) return false;
  }
  if( fwrite(&dmy, sizeof(int), 1, fp) != 1 ) return false;
  
  if( d_type == SPH_FLOAT ) {
    dmy = sizeof(int) + sizeof(float);
  } else {
    dmy = sizeof(long long) + sizeof(double);
  }
  //dmy = sizeof(long long) + sizeof(double);
  if( fwrite(&dmy, sizeof(int), 1, fp) != 1 ) return false;
  if( d_type == SPH_FLOAT ) {
    float ftmp = (float)time;
    if( fwrite(&step, sizeof(int), 1, fp) != 1 ) return false;
    if( fwrite(&ftmp, sizeof(float), 1, fp) != 1 ) return false;
  } else {
    if( fwrite(&step, sizeof(long long), 1, fp) != 1 ) return false;
    if( fwrite(&time, sizeof(double), 1, fp) != 1 ) return false;
  }
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

