/* #################################################################
 *
 * CIOlib - Cartesian Input / Output library
 *
 * Copyright (c) AICS, RIKEN. All right reserved. 2013
 *
 * #################################################################
 */

/** 
 * @file   cio_DFI_SPH.C
 * @brief  cio_DFI_SPH Class
 * @author kero    
 */

#include "cio_DFI.h"
#include "cio_DFI_SPH.h"

// #################################################################
// コンストラクタ
cio_DFI_SPH::cio_DFI_SPH()
{

}


// #################################################################
// デストラクタ
cio_DFI_SPH::~cio_DFI_SPH()
{

}


// #################################################################
// 
void cio_DFI_SPH::ReadData(int step, int gc, 
                           int Gvoxel[3], int Gdivision[3], int head[3], int tail[3],
                           REAL_TYPE *val, REAL_TYPE &time)
{

  int RankID;
  MPI_Comm_rank( m_comm, &RankID );

  int nrank;
  MPI_Comm_size( m_comm, &nrank );

  bool mio = false;
//  if( nrank > 1 ) mio = true;
  //printf("nrank : %d RankInfo.size() : %d\n",nrank,RankInfo.size());
  //printf("nrank : %d Gvoxel                 : %d %d %d\n",nrank,Gvoxel[0],Gvoxel[1],Gvoxel[2]);
  //printf("nrank : %d DFI_Domain.GlobalVoxel : %d %d %d\n",nrank,
  //              DFI_Domain.GlobalVoxel[0],DFI_Domain.GlobalVoxel[1],DFI_Domain.GlobalVoxel[2]);
  
  if( RankInfo.size() > 1 ) mio = true;

  cio_Domain G_Domain;
  vector<cio_Rank> G_RankInfo;
  cio_Rank G_Rank;
  cio_Create_Domain(m_comm,Gvoxel, Gdivision, head, tail,G_Domain,G_RankInfo,G_Rank);

  cio_EGlobalVoxel readflag;  //CIO_E_GV_SAME:密　CIO_E_GVX2_SAME:粗 CIO_E_OTHER:その他
  readflag = CheckGlobalVoxel(Gvoxel, DFI_Domain.GlobalVoxel);
  //if( readflag == CIO_E_GV_SAME )   printf("same file\n");
  //if( readflag == CIO_E_GVX2_SAME ) printf("x2 same file\n");
  //if( readflag == CIO_E_OTHER )     printf("other file\n");

  vector<int> RList;
  CreateRankList(head, tail, gc, readflag, RList);
  //printf("RList.size() : %d\n",RList.size());

//密データの読込み処理
  if( readflag == CIO_E_GV_SAME ) { 
    for(int i=0; i<RList.size(); i++) {
      std::string fname = Generate_FileName(RList[i],step,mio);

      //printf("fname : %s\n",fname.c_str());
      int n = RList[i];

      int sta[3],end[3];
      if( CheckReadArea(head, tail, gc, RankInfo[n].HeadIndex, RankInfo[n].TailIndex,
                        DFI_Finfo.GuideCell, readflag, sta, end) )
      { 

        //printf("*** 1x1 read ***\n");

        //1対1の読込み
        if( !read_S3D(fname.c_str(), step, DFI_Finfo.GuideCell, val, time ) ) {
          val = NULL;
          printf("**** error read s3d ****");
          return;
        }

      } else {

        //printf("*** mxn read ***\n");
        //printf("ID : %d dfi_ID %d sta : %d %d %d end : %d %d %d\n",RankID,n,
        //        sta[0],sta[1],sta[2],end[0],end[1],end[2]);

        //1対多の読込み
        if( !read_MxN(fname.c_str(), step, gc, head, tail, RankInfo[n].HeadIndex,
                    RankInfo[n].TailIndex, sta, end, val, time, n) ) {
          val = NULL;
          printf("**** error read mxn ID %d dfi_ID : %d fname : %s\n",RankID,n,
                  fname.c_str());
          return;
        }
      }
    }
    return;
//粗い粒子の処理
  }else if( readflag == CIO_E_GVX2_SAME ) { 

    for(int i=0; i<RList.size(); i++) {
      std::string fname = Generate_FileName(RList[i],step,mio);
      //printf("fname : %s\n",fname.c_str());
      int sta[3],end[3],dfi_head[3],dfi_tail[3];

      int n = RList[i];

      for(int j=0; j<3; j++){
        dfi_head[j]=RankInfo[n].HeadIndex[j]*2-1;
        dfi_tail[j]=RankInfo[n].TailIndex[j]*2;
      }
      //1対1の読込み&補間
      if( CheckReadArea(head, tail, gc, dfi_head, dfi_tail,DFI_Finfo.GuideCell, readflag, 
                        sta, end)) 
      {

        //printf("*** 1x1 read ***\n");
        //printf("ID : %d dfi_ID %d sta : %d %d %d end : %d %d %d\n",RankID,n,
        //        sta[0],sta[1],sta[2],end[0],end[1],end[2]);

        if( !read_Coarse(fname.c_str(), step, gc, head, tail, Gvoxel, RankInfo[n].HeadIndex,
                          RankInfo[n].TailIndex, RankInfo[n].VoxelSize, DFI_Finfo.GuideCell,
                          DFI_Finfo.Component, sta, end, val, time) ) {
           val = NULL;
           return ;
        }
      } else {
      //1対多の読込み&補間
      //
        //printf("*** mxn ***\n");
        //printf("ID : %d dfi_ID %d sta : %d %d %d end : %d %d %d\n",RankID,n,
        //        sta[0],sta[1],sta[2],end[0],end[1],end[2]);

        if( !read_Coarse_MxN(fname.c_str(), step, gc, head, tail, Gvoxel, RankInfo[n].HeadIndex,
                          RankInfo[n].TailIndex, RankInfo[n].VoxelSize, DFI_Finfo.GuideCell,
                          DFI_Finfo.Component, sta, end, val, time, n) ) {
           val = NULL;
           return ;
        }
      }
    }
  } else {
    printf("Error field data format\n");
    return;
  }

}

// #################################################################
// 
bool cio_DFI_SPH::read_Head(FILE* fp, int Etype, int step, INT_TYPE voxsize[3], REAL_TYPE &time)
{

  unsigned int dmy;
  if( fread(&dmy, sizeof(int), 1, fp) != 1 ) { fclose(fp); return false; }
  //if( Etype = CIO_UnMatch )  
  if( dmy != 8 ) { fclose(fp); return false; }

  DataDims data_dims;
  if( fread(&data_dims, sizeof(int), 1, fp) != 1 ) { fclose(fp); return false; }
  if( data_dims == _SCALAR && DFI_Finfo.Component != 1 ) { fclose(fp); return false; } 
  if( data_dims == _VECTOR && DFI_Finfo.Component <= 1 ) { fclose(fp); return false; } 

  RealType real_type;
  if( fread(&real_type, sizeof(int), 1, fp) != 1 ) { fclose(fp); return false; }

  if( fread(&dmy, sizeof(int), 1, fp) != 1 ) { fclose(fp); return false; }
  if( dmy != 8 ) { fclose(fp); return false; }

//ボクセルサイズ
  if( fread(&dmy, sizeof(int), 1, fp) != 1 ) { fclose(fp); return false; }
  if( fread(voxsize, sizeof(INT_TYPE), 3, fp) != 3 ){fclose(fp);return false;}
  if( fread(&dmy, sizeof(int), 1, fp) != 1 ) { fclose(fp); return false; }

//原点座標
  if( fread(&dmy, sizeof(int), 1, fp) != 1 ) { fclose(fp); return false; }
  REAL_TYPE voxorg[3];
  if( fread(voxorg, sizeof(REAL_TYPE), 3, fp) != 3 ){fclose(fp);return false;}
  if( fread(&dmy, sizeof(int), 1, fp) != 1 ) { fclose(fp); return false; }
  
//pit
  if( fread(&dmy, sizeof(int), 1, fp) != 1 ) { fclose(fp); return false; }
  REAL_TYPE voxpit[3];
  if( fread(voxpit, sizeof(REAL_TYPE), 3, fp) != 3 ){fclose(fp);return false;}
  if( fread(&dmy, sizeof(int), 1, fp) != 1 ) { fclose(fp); return false; }

//step,time
  INT_TYPE r_step;
  REAL_TYPE r_time;
  if( fread(&dmy, sizeof(int), 1, fp) != 1 ) { fclose(fp); return false; }
  if( fread(&r_step, sizeof(INT_TYPE), 1, fp) != 1 ) { fclose(fp); return false; }
  if( fread(&r_time, sizeof(REAL_TYPE), 1, fp) != 1 ) { fclose(fp); return false; }
  if( fread(&dmy, sizeof(int), 1, fp) != 1 ) { fclose(fp); return false; }
  if( r_step != step ) { fclose(fp); return false; }

  time = r_time;

return true;

}

// #################################################################
// 
bool cio_DFI_SPH::read_S3D(const char* fname, int step, int gc, REAL_TYPE*val, REAL_TYPE &time)
{

  if( !fname || !DFI_Finfo.Component ) return false;

  //printf("fname      : %s\n",fname);
  //printf("DataType   : %s\n",DFI_Finfo.DataType.c_str());
  //printf("ArrayShape : %s\n",DFI_Finfo.ArrayShape.c_str());
  //printf("Component  : %d\n",DFI_Finfo.Component);

  //Endian check....
  int idumy = 1;
  char* cdumy = (char*)(&idumy);
  string Endian="";
  if( cdumy[0] == 0x01 ) Endian = "little";
  if( cdumy[0] == 0x00 ) Endian = "big";

  int Etype = CIO_UnKnown;
  if( DFI_Finfo.Endian == Endian ) {
    Etype = CIO_Match;
  }else {
    Etype = CIO_UnMatch;
  }

  if( Etype == CIO_UnKnown ) return false;

  FILE* fp;
  if( !(fp=fopen(fname,"rb")) ) {
    printf("Can't open file. (%s)\n",fname);
    return false;
  }

  INT_TYPE voxsize[3];
  if( !read_Head(fp, Etype, step, voxsize, time) ) { fclose(fp); return false; }
  //printf("voxsize : %d %d %d\n",voxsize[0],voxsize[1],voxsize[2]);

  unsigned int dmy;
//data 
  int arraylen = voxsize[0]*voxsize[1]*voxsize[2]*DFI_Finfo.Component;
  if( fread(&dmy, sizeof(int), 1, fp) != 1 ) { fclose(fp); return false; }
  if( dmy != arraylen*sizeof(REAL_TYPE ) ) { fclose(fp); return false; }
  if( fread(val, sizeof(REAL_TYPE), arraylen, fp) != arraylen ) { fclose(fp); return false; }
  if( fread(&dmy, sizeof(int), 1, fp) != 1 ) { fclose(fp); return false; }
 
  //printf("*******************\n");
  //printf("*** read sph ok ***\n");
  //printf("*******************\n");

  fclose(fp);
  return true;

}

// #################################################################
// 
bool cio_DFI_SPH::read_MxN(const char* fname, int step, int gc, int head[3], int tail[3], 
                           int R_head[3], int R_tail[3], int sta[3], int end[3], REAL_TYPE*val, 
                           REAL_TYPE &time, int n_dfi)
{

  if( !fname || !DFI_Finfo.Component ) return false;


  int RankID;
  MPI_Comm_rank( m_comm, &RankID );

  //printf("RankID : %d n_dfi : %d head   : %d %d %d\n",RankID,n_dfi,head[0],head[1],head[2]);
 // printf("RankID : %d n_dfi : %d tail   : %d %d %d\n",RankID,n_dfi,tail[0],tail[1],tail[2]);
  //printf("RankID : %d n_dfi : %d R_head : %d %d %d\n",RankID,n_dfi,R_head[0],R_head[1],R_head[2]);
  //printf("RankID : %d n_dfi : %d R_tail : %d %d %d\n",RankID,n_dfi,R_tail[0],R_tail[1],R_tail[2]);
  //printf("RankID : %d n_dfi : %d sta    : %d %d %d\n",RankID,n_dfi,sta[0],sta[1],sta[2]);
  //printf("RankID : %d n_dfi : %d end    : %d %d %d\n",RankID,n_dfi,end[0],end[1],end[2]);
  //printf("RankID : %d n_dfi : %d step   : %d \n",RankID,n_dfi,step);
  //printf("RankID : %d n_dfi : %d gc     : %d \n",RankID,n_dfi,gc);
  //printf("RankID : %d n_dfi : %d dfi_gc : %d \n",RankID,n_dfi,DFI_Finfo.GuideCell);

/*
  char x_fname[128];
  sprintf( x_fname, "dfi_%d_%d.txt",RankID,n_dfi);
  printf("x_fname : %s\n",x_fname);  
  FILE* xfp;
  if( !(xfp=fopen(x_fname,"w")) ) {
    printf("Can't open file. (%s)\n",x_fname);
    return false;
  }
*/

  //Endian check....
  int idumy = 1;
  char* cdumy = (char*)(&idumy);
  string Endian="";
  if( cdumy[0] == 0x01 ) Endian = "little";
  if( cdumy[0] == 0x00 ) Endian = "big";

  int Etype = CIO_UnKnown;
  if( DFI_Finfo.Endian == Endian ) {
    Etype = CIO_Match;
  }else {
    Etype = CIO_UnMatch;
  }

  if( Etype == CIO_UnKnown ) return false;

  FILE* fp;
  if( !(fp=fopen(fname,"rb")) ) {
    printf("Can't open file. (%s)\n",fname);
    return false;
  }

  INT_TYPE voxsize[3];
  if( !read_Head(fp, Etype, step, voxsize, time) ) { fclose(fp); return false; }

  //printf("voxsize : %d %d %d\n",voxsize[0],voxsize[1],voxsize[2]);

  int off_set;
  long off_setleng;
//ijkn
  //if( DFI_Finfo.ArrayShape == "ijkn" ) { 
  if( !strcasecmp( DFI_Finfo.ArrayShape.c_str(), "ijkn" ) || DFI_Finfo.Component == 1 ) { 
    off_set = (head[2]-gc)-(R_head[2]-DFI_Finfo.GuideCell);
    //off_setleng = (voxsize[0]+2*gc)*(voxsize[1]+2*gc)*(off_set+2*gc);
    off_setleng = (voxsize[0])*(voxsize[1])*(off_set);
  } else if( !strcasecmp(DFI_Finfo.ArrayShape.c_str(), "nijk" ) ) {
//nijk
    off_set = (head[2]-gc)-(R_head[2]-DFI_Finfo.GuideCell);
    //off_setleng = (voxsize[0]+2*gc)*(voxsize[1]+2*gc)*(off_set+2*gc)*DFI_Finfo.Component;
    off_setleng = (voxsize[0])*(voxsize[1])*(off_set)*DFI_Finfo.Component;
  }

  int nkrec = end[2]-sta[2]+1;
  //printf("nkrec : %d\n",nkrec);

  unsigned int dmy;
  if( fread(&dmy, sizeof(int), 1, fp) != 1 ) { fclose(fp); return false; }

  long long idx,idx2;
  int imax,jmax,kmax;

  imax = tail[0]-head[0]+1;
  jmax = tail[1]-head[1]+1;
  kmax = tail[2]-head[2]+1;

  //printf("imax : %d jmax : %d kmax : %d\n",imax,jmax,kmax);

  // ijkn
  if( !strcasecmp(DFI_Finfo.ArrayShape.c_str(), "ijkn") || DFI_Finfo.Component == 1 ) {

    //printf("ArrayShape : ijkn\n");

    //int arraylen = (voxsize[0]+2*gc)*(voxsize[1]+2*gc);
    int arraylen = (voxsize[0])*(voxsize[1]);
    REAL_TYPE *rbuff = new REAL_TYPE[arraylen];

    for(int ncomp=0; ncomp<DFI_Finfo.Component; ncomp++) {
      if( off_set > 0 ) {
        //printf("ID : %d off_set : %d off_setleng : %d\n",RankID,off_set,off_setleng);
        if( fseek(fp,off_setleng*sizeof(REAL_TYPE),SEEK_CUR) != 0 ) 
        {
          printf("*** Error fseek length offset length : %d\n", off_set*sizeof(REAL_TYPE));
          fclose(fp); return false; 
        }
      }

      //for(int n=0; n<nkrec+2*gc; n++) {
      int iimax,jjmax;
      iimax=voxsize[0]-2*DFI_Finfo.GuideCell;
      jjmax=voxsize[1]-2*DFI_Finfo.GuideCell;
      if( sta[2]>0 && R_head[2]>1) {
        for(int n=0; n<DFI_Finfo.GuideCell; n++) {
           if( fread(rbuff, sizeof(REAL_TYPE), arraylen, fp) != arraylen ) { fclose(fp); return false; }
        }
        nkrec-=DFI_Finfo.GuideCell;
      }

      //printf("read nkrec : %d\n",nkrec);

      for(int n=0; n<nkrec; n++) {
        if( fread(rbuff, sizeof(REAL_TYPE), arraylen, fp) != arraylen ) { fclose(fp); return false; }
        int k=sta[2]+n;
        for( int j=sta[1]; j<=end[1]; j++ ) {
        for( int i=sta[0]; i<=end[0]; i++ ) {
           idx = _IDX_IJKN(i-(head[0]-1),j-(head[1]-1),k-(head[2]-1),ncomp,imax,jmax,kmax,gc);
           //idx2= _IDX_IJ(i-(R_head[0]-1),j-(R_head[1]-1),voxsize[0],voxsize[1],DFI_Finfo.GuideCell);
           idx2= _IDX_IJ(i-(R_head[0]-1),j-(R_head[1]-1),iimax,jjmax,DFI_Finfo.GuideCell);
           val[idx]=rbuff[idx2];

           //fprintf(xfp,"idx : %d idx2 : %d\n",idx,idx2);

        }}
      }
    }

    delete [] rbuff;

  // nijk
  } else if( !strcasecmp(DFI_Finfo.ArrayShape.c_str(), "nijk" ) ) {

    //printf("ArrayShape : nijk\n");

    //int arraylen = (voxsize[0]+2*gc)*(voxsize[1]+2*gc)*DFI_Finfo.Component;
    int arraylen = (voxsize[0])*(voxsize[1])*DFI_Finfo.Component;
    REAL_TYPE *rbuff = new REAL_TYPE[arraylen];

    //printf("nkrec : %d voxsize[0] : %d voxsize[1] : %d DFI_Finfo.Component : %d arraylen : %d\n",
    //       nkrec,voxsize[0],voxsize[1],DFI_Finfo.Component,arraylen);

    if( off_set > 0 ) {
      //printf("off_setleng : %d\n",off_setleng);
      if( fseek(fp,off_setleng*sizeof(REAL_TYPE),SEEK_CUR) != 0 ) 
      {
        printf("*** Error fseek length offset length : %d\n", off_set*sizeof(REAL_TYPE));
        fclose(fp); return false; 
      }
    }

    //for(int n=0; n<nkrec+2*gc; n++) {
    int iimax,jjmax;
    iimax=voxsize[0]-2*DFI_Finfo.GuideCell;
    jjmax=voxsize[1]-2*DFI_Finfo.GuideCell;
    if( sta[2]>0 && R_head[2]>1 ) {
      for(int n=0; n<DFI_Finfo.GuideCell; n++) {
         if( fread(rbuff, sizeof(REAL_TYPE), arraylen, fp) != arraylen ) { fclose(fp); return false; }
      }
      nkrec-=DFI_Finfo.GuideCell;
    }

    //printf("read nkrec : %d\n",nkrec);

    for(int n=0; n<nkrec; n++) {
      if( fread(rbuff, sizeof(REAL_TYPE), arraylen, fp) != arraylen ) { fclose(fp); return false; }
      int k=sta[2]+n;
      for( int j=sta[1]; j<=end[1]; j++ ) {
      for( int i=sta[0]; i<=end[0]; i++ ) {
      for(int ncomp=0; ncomp<DFI_Finfo.Component; ncomp++) {
         idx = _IDX_NIJK(ncomp,i-(head[0]-1),j-(head[1]-1),k-(head[2]-1),DFI_Finfo.Component,imax,jmax,kmax,gc);
         //idx2= _IDX_NIJ(ncomp,i-(R_head[0]-1),j-(R_head[1]-1),voxsize[0],voxsize[1],DFI_Finfo.Component,
         //               DFI_Finfo.GuideCell);
         idx2= _IDX_NIJ(ncomp,i-(R_head[0]-1),j-(R_head[1]-1),iimax,jjmax,DFI_Finfo.Component,
                          DFI_Finfo.GuideCell);
         val[idx]=rbuff[idx2];

         //fprintf(xfp,"idx : %d idx2 : %d\n",idx,idx2);

      }}}
    }

    delete [] rbuff;

  }


  //printf("*******************\n");
  //printf("*** read sph ok ***\n");
  //printf("*******************\n");

  fclose(fp);
  //fclose(xfp);

  return true;

}

// #################################################################
// 
bool cio_DFI_SPH::read_Coarse(const char* fname, int step, int gc, int head[3], int tail[3], 
                           int G_Voxelsize[3], int R_head[3], int R_tail[3], int Voxelsize[3], 
                           int dfi_gc, int ncomp, int sta[3], int end[3], REAL_TYPE*val, 
                           REAL_TYPE &time)
{

  REAL_TYPE *src = NULL;
  int size_val = (Voxelsize[0]+2*dfi_gc)*(Voxelsize[1]+2*dfi_gc)*(Voxelsize[2]+2*dfi_gc)*ncomp;
  src = new REAL_TYPE[size_val];

  if( !read_S3D(fname, step, dfi_gc, src, time) ) {
    delete [] src;
    val = NULL;
    printf("**** error read s3d ****");
    return false;
  }

  int RankID;
  MPI_Comm_rank( m_comm, &RankID );  

  if( !strcasecmp(DFI_Finfo.ArrayShape.c_str(), "ijkn") || ncomp == 1 ) {
    cio_dfi_coarse_ijkn_(val, head, tail, &gc, src, R_head, R_tail, &dfi_gc, &ncomp, sta, end, &RankID);
  } else {
    cio_dfi_coarse_nijk_(val, head, tail, &gc, src, R_head, R_tail, &dfi_gc, &ncomp, sta, end, &RankID);
  }
  delete [] src;

  return true;

}

// #################################################################
// 
bool cio_DFI_SPH::read_Coarse_MxN(const char* fname, int step, int gc, int head[3], int tail[3], 
                           int G_Voxelsize[3], int R_head[3], int R_tail[3], int Voxelsize[3], 
                           int dfi_gc, int ncomp, int sta[3], int end[3], REAL_TYPE *val, 
                           REAL_TYPE &time, int n_dfi)
{

  int read_s[3],read_e[3];
  int temp_head[3],temp_tail[3],t_voxsize[3];
  for(int i=0; i<3; i++ ) {
    temp_head[i]=(head[i]-1)/2+1;
    temp_tail[i]=(tail[i]-1)/2+1;
    if( temp_head[i] > R_head[i] ) temp_head[i] = R_head[i];
    if( temp_tail[i] < R_tail[i] ) temp_tail[i] = R_tail[i];
    t_voxsize[i]=(temp_tail[i]-temp_head[i])+1;
  }

  //cio_EGlobalVoxel readflag = CIO_E_GVX2_SAME;
  cio_EGlobalVoxel readflag = CIO_E_GV_SAME;
  //CheckReadArea(temp_head,temp_tail,gc,R_head,R_tail,dfi_gc,readflag,read_s,read_e);
  CheckReadArea(temp_head,temp_tail,dfi_gc,R_head,R_tail,dfi_gc,readflag,read_s,read_e);

  REAL_TYPE *src = NULL;
  int size_val = (t_voxsize[0]+2*dfi_gc)*(t_voxsize[1]+2*dfi_gc)*(t_voxsize[2]+2*dfi_gc)*ncomp;
  //int size_val = (t_voxsize[0]+2*gc)*(t_voxsize[1]+2*gc)*(t_voxsize[2]+2*gc)*ncomp;
  src = new REAL_TYPE[size_val];

  //if( !read_MxN(fname,step,gc,temp_head,temp_tail,R_head,R_tail,read_s,read_e,src,n_dfi) ) { 
  if( !read_MxN(fname,step,dfi_gc,temp_head,temp_tail,R_head,R_tail,read_s,read_e,src,time,n_dfi) ) { 
    delete [] src;
    val = NULL;
    printf("**** error read Mxn ****");
    return false;
  }

  int RankID;
  MPI_Comm_rank( m_comm, &RankID );  

  if( !strcasecmp(DFI_Finfo.ArrayShape.c_str(), "ijkn") || ncomp == 1 ) {
    cio_dfi_coarse_ijkn_(val, head, tail, &gc, src, temp_head, temp_tail, &dfi_gc, &ncomp, sta, end, &RankID);
  } else {
    cio_dfi_coarse_nijk_(val, head, tail, &gc, src, temp_head, temp_tail, &dfi_gc, &ncomp, sta, end, &RankID);
  }

  delete [] src;

  return true;

}

// #################################################################
// 
void cio_DFI_SPH::WriteData(int step, int gc, REAL_TYPE time,
                           REAL_TYPE *val, REAL_TYPE *minmax)
{

  int RankID;
  MPI_Comm_rank( m_comm, &RankID );

  int nrank;
  MPI_Comm_size( m_comm, &nrank );
  bool mio=false;
  if( nrank>1 ) mio=true;

  std::string outFile = Generate_FileName(RankID,step,mio);
  if( MakeDirectory(DFI_Finfo.DirectoryPath) != 1 ) return;

  FILE* fp;
  if( (fp = fopen(outFile.c_str(),"wb")) == NULL ) {
    fprintf(stderr,"Can't open file.(%s)\n",outFile.c_str());
    return;
  }

  //if( !write_head(fp, step, time, gc, RankID) ) { fclose(fp); return; }
  if( !write_head(fp, step, time, RankID) ) { fclose(fp); return; }
  //printf(" sph hedaer record write ok!!!\n");

  if( gc == DFI_Finfo.GuideCell ) { 
    if( !write_data(fp, val, gc, RankID) ) { fclose(fp); return; }
    //printf(" sph data record write ok!!!\n");
  } else {
    REAL_TYPE *src = NULL;
    int sz[3];
    for(int i=0; i<3; i++ ) sz[i] = RankInfo[RankID].VoxelSize[i]+2*DFI_Finfo.GuideCell;
    int size_val = sz[0]*sz[1]*sz[2]*DFI_Finfo.Component;
    //printf("sz : %d %d %d size_val : %d\n",sz[0],sz[1],sz[2],size_val);
    src = new REAL_TYPE[size_val];
    if( !strcasecmp(DFI_Finfo.ArrayShape.c_str(), "ijkn") || DFI_Finfo.Component == 1 ) {
      cio_dfi_copy_ijkn_(src,sz,&DFI_Finfo.GuideCell,val,&gc,&DFI_Finfo.Component); 
    } else {
      cio_dfi_copy_nijk_(src,sz,&DFI_Finfo.GuideCell,val,&gc,&DFI_Finfo.Component); 
    }
    if( !write_data(fp, src, DFI_Finfo.GuideCell, RankID) ) { 
      fclose(fp); 
      delete [] src;
      return;
    }
    //printf(" sph data record write ok!!!\n");

    delete [] src;
  }

  fclose(fp);

  //if( step > 0 ) {
    WriteIndexDfiFile(m_indexDfiName,RankID,DFI_Finfo.Prefix,step,time,minmax,true);
  //}
}

// #################################################################
// 
//bool cio_DFI_SPH::write_head(FILE* fp, int step, REAL_TYPE time, int gc, int n)
bool cio_DFI_SPH::write_head(FILE* fp, int step, REAL_TYPE time, int n)
{

  int svType = 0;
  if( DFI_Finfo.Component == 1 ) svType = 1;
  if( DFI_Finfo.Component > 1  ) svType = 2;
  if( svType == 0 ) return false;

  int dType = 0;
  if( sizeof(INT_TYPE) == 4 ) dType = 1;
  if( sizeof(INT_TYPE) == 8 ) dType = 2;
  if( dType == 0 ) return false;

  unsigned int dmy;
  dmy = 8;
  if( fwrite(&dmy, sizeof(int), 1, fp) != 1 ) return false;
  if( fwrite(&svType, sizeof(int), 1, fp) != 1 ) return false;
  if( fwrite(&dType, sizeof(int), 1, fp) != 1 ) return false;
  if( fwrite(&dmy, sizeof(int), 1, fp) != 1 ) return false;


  if( dType == 1 ) dmy = 12; //float
  else             dmy = 24; //double

  //voxel size
  INT_TYPE size[3];
  //for(int i=0; i<3; i++ ) size[i] = (INT_TYPE)RankInfo[n].VoxelSize[i]+(INT_TYPE)(2*gc);
  for(int i=0; i<3; i++ ) size[i] = (INT_TYPE)RankInfo[n].VoxelSize[i]+(INT_TYPE)(2*DFI_Finfo.GuideCell);
  if( fwrite(&dmy, sizeof(int), 1, fp) != 1 ) return false;
  if( fwrite(size, sizeof(INT_TYPE), 3, fp) !=3 ) return false;
  if( fwrite(&dmy, sizeof(int), 1, fp) != 1 ) return false;

  REAL_TYPE pch[3],org[3];
  for(int i=0; i<3; i++ ) pch[i]=DFI_Domain.GlobalRegion[i]/DFI_Domain.GlobalVoxel[i];
  //for(int i=0; i<3; i++ ) org[i]=(DFI_Domain.GlobalOrigin[i]+pch[i]*(RankInfo[n].HeadIndex[i]-1));
  for(int i=0; i<3; i++ ) org[i]=DFI_Domain.GlobalOrigin[i];
  //origin
  if( fwrite(&dmy, sizeof(int), 1, fp) != 1 ) return false;
  if( fwrite(org, sizeof(REAL_TYPE), 3, fp) !=3 ) return false;
  if( fwrite(&dmy, sizeof(int), 1, fp) != 1 ) return false;

  //pitch
  if( fwrite(&dmy, sizeof(int), 1, fp) != 1 ) return false;
  if( fwrite(pch, sizeof(REAL_TYPE), 3, fp) !=3 ) return false;
  if( fwrite(&dmy, sizeof(int), 1, fp) != 1 ) return false;

  //step&time
  INT_TYPE dstep = (INT_TYPE)step;
  if ( dType == 1 ) dmy = 8;
  else              dmy = 16;
  if( fwrite(&dmy, sizeof(int), 1, fp) != 1 ) return false;
  if( fwrite(&dstep, sizeof(INT_TYPE), 1, fp) != 1 ) return false;
  if( fwrite(&time, sizeof(REAL_TYPE), 1, fp) != 1 ) return false;
  if( fwrite(&dmy, sizeof(int), 1, fp) != 1 ) return false;
  
  return true;
}

// #################################################################
// 
bool cio_DFI_SPH::write_data(FILE* fp, REAL_TYPE* val, int gc, int n)
{

  INT_TYPE size[3];
  for(int i=0; i<3; i++ ) size[i] = (INT_TYPE)RankInfo[n].VoxelSize[i]+(INT_TYPE)(2*gc);

  size_t dLen = (size_t)(size[0] * size[1] * size[2]);
  if( DFI_Finfo.Component > 1 ) dLen *= 3;  

  unsigned int dmy = dLen * sizeof(REAL_TYPE);

  if( fwrite(&dmy, sizeof(int), 1, fp) != 1 ) return false;
  if( fwrite((const REAL_TYPE*)val, sizeof(REAL_TYPE), dLen, fp) != dLen ) return false; 
  if( fwrite(&dmy, sizeof(int), 1, fp) != 1 ) return false;

  return true;
}
