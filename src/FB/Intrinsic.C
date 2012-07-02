// #################################################################
//
// CAERU Library
//
// Copyright (c) All right reserved. 2012
//
// Institute of Industrial Science, The University of Tokyo, Japan. 
//
// #################################################################

/** 
 * @file   Intrinsic.C
 * @brief  FlowBase Intrinsic class
 * @author kero
 */

#include "Intrinsic.h"



// 例題名称の表示
void Intrinsic::printExample(FILE* fp, const char* str)
{
  if( !fp ) {
    stamped_printf("\tFail to write into file\n");
    Exit(0);
  }
  fprintf(fp,"\n\tExample : %s\n\n", str); fflush(fp);

  printf("\n\tExample : %s\n\n", str);
}



// パラメータの表示
void Intrinsic::printPara(FILE* fp, const Control* R)
{
  if ( !fp ) {
    stamped_printf("\tFail to write into file\n");
    Exit(0);
  }
  // describe method
}



// 例題のモデルをsvxフォーマットで出力する(体積率とID)
void Intrinsic::writeSVX(REAL_TYPE *vf, int *id, Control* R)
{

  int    sz, ix, jx, kx;
  size_t m, l;
  float  ox, oy, oz, dx, dy, dz;

  char svx_fname[512];

  if ( paraMngr->IsParallel() ) {
    sprintf( svx_fname, "example_%06d.svx", paraMngr->GetMyRankID() );
  } else {
    sprintf( svx_fname, "example.svx" );
  }
  ofstream ofs(svx_fname, ios::out | ios::binary);
  if (!ofs) {
    cout << "\tCan't open " << svx_fname << " file" << endl;
    Exit(0);
  }

  int imax = size[0];
  int jmax = size[1];
  int kmax = size[2];
  int gd = guide;
  
  ix = imax+2;  // +2 means guide cell for IP model
  jx = jmax+2;
  kx = kmax+2;
  
  size_t nx = ix*jx*kx;
  
  dx = (float)R->dx[0]*RefL;
  dy = (float)R->dx[1]*RefL;
  dz = (float)R->dx[2]*RefL;
  ox = (float)R->org[0]*RefL - dx; // 片側1層分をシフト
  oy = (float)R->org[1]*RefL - dy;
  oz = (float)R->org[2]*RefL - dz;
  
  //stamped_printf("example out org(%e %e %e) dimensional\n", ox, oy, oz);

  float *f = new float[nx];
  int   *q = new int[nx];

  for (int k=0; k<=(kmax+1); k++) {
    for (int j=0; j<=(jmax+1); j++) {
      for (int i=0; i<=(imax+1); i++) {
        l = (size_t)(ix*jx*k + ix*j + i);
        m = _F_IDX_S3D(i, j, k, imax, jmax, kmax, gd); //FBUtility::getFindexS3D(size, guide, i, j, k);
        q[l] = id[m];
        f[l] = (float)vf[m];
      }
    }
  }

  // voxel size
  sz = sizeof(int)*3;
  ofs.write( (char*)&sz, sizeof(int) );
  ofs.write( (char*)&ix, sizeof(int) );
  ofs.write( (char*)&jx, sizeof(int) );
  ofs.write( (char*)&kx, sizeof(int) );
  ofs.write( (char*)&sz, sizeof(int) );

  // original point of domain
  sz = sizeof(float)*3;
  ofs.write( (char*)&sz, sizeof(int) );
  ofs.write( (char*)&ox, sizeof(float) );
  ofs.write( (char*)&oy, sizeof(float) );
  ofs.write( (char*)&oz, sizeof(float) );
  ofs.write( (char*)&sz, sizeof(int) );

  // pitch of voxel
  ofs.write( (char*)&sz, sizeof(int) );
  ofs.write( (char*)&dx, sizeof(float) );
  ofs.write( (char*)&dy, sizeof(float) );
  ofs.write( (char*)&dz, sizeof(float) );
  ofs.write( (char*)&sz, sizeof(int) );

  // type of stored data
  sz = sizeof(int)*1;
  int dtype = 0;
  dtype |= ( 0x1<<0 );  // volume fraction
  dtype |= ( 0x1<<2 );  // medium ID
  ofs.write( (char*)&sz,  sizeof(int) );
  ofs.write( (char*)&dtype, sizeof(int) );
  ofs.write( (char*)&sz,  sizeof(int) );

  // volume fraction
  sz = nx * sizeof(float);
  ofs.write( (char*)&sz, sizeof(int) );
  ofs.write( (char*)f,   sz );
  ofs.write( (char*)&sz, sizeof(int) );

  // medium ID
  sz = nx * sizeof(int);
  ofs.write( (char*)&sz, sizeof(int) );
  ofs.write( (char*)q,   sz );
  ofs.write( (char*)&sz, sizeof(int) );

  ofs.close();

  if (f) { delete [] f; f=NULL; }
  if (q) { delete [] q; q=NULL; }
}


// 例題のモデルをsvxフォーマットで出力する(ID)
void Intrinsic::writeSVX(int *id, Control* R)
{
  
  int   sz, ix, jx, kx;
  size_t m, l;
  float ox, oy, oz, dx, dy, dz;
  
  char svx_fname[512];
  
  if ( paraMngr->IsParallel() ) {
    sprintf( svx_fname, "example_%06d.svx", paraMngr->GetMyRankID() );
  } else {
    sprintf( svx_fname, "example.svx" );
  }
  ofstream ofs(svx_fname, ios::out | ios::binary);
  if (!ofs) {
    cout << "\tCan't open " << svx_fname << " file" << endl;
    Exit(0);
  }
  
  int imax = size[0];
  int jmax = size[1];
  int kmax = size[2];
  int gd = guide;
  
  ix = imax+2;  // +2 means guide cell for IP model
  jx = jmax+2;
  kx = kmax+2;
  
  size_t nx = (size_t)(ix*jx*kx);
  
  dx = (float)R->dx[0]*RefL;
  dy = (float)R->dx[1]*RefL;
  dz = (float)R->dx[2]*RefL;
  ox = (float)R->org[0]*RefL - dx; // 片側1層分をシフト
  oy = (float)R->org[1]*RefL - dy;
  oz = (float)R->org[2]*RefL - dz;
  
  //stamped_printf("example out org(%e %e %e) dimensional\n", ox, oy, oz);
  
  int   *q = new int[nx];
  
  for (int k=0; k<=(kmax+1); k++) {
    for (int j=0; j<=(jmax+1); j++) {
      for (int i=0; i<=(imax+1); i++) {
        l = (size_t)(ix*jx*k + ix*j + i);
        m = _F_IDX_S3D(i, j, k, imax, jmax, kmax, gd); //FBUtility::getFindexS3D(size, guide, i, j, k);
        q[l] = id[m];
      }
    }
  }
  
  // voxel size
  sz = sizeof(int)*3;
  ofs.write( (char*)&sz, sizeof(int) );
  ofs.write( (char*)&ix, sizeof(int) );
  ofs.write( (char*)&jx, sizeof(int) );
  ofs.write( (char*)&kx, sizeof(int) );
  ofs.write( (char*)&sz, sizeof(int) );
  
  // original point of domain
  sz = sizeof(float)*3;
  ofs.write( (char*)&sz, sizeof(int) );
  ofs.write( (char*)&ox, sizeof(float) );
  ofs.write( (char*)&oy, sizeof(float) );
  ofs.write( (char*)&oz, sizeof(float) );
  ofs.write( (char*)&sz, sizeof(int) );
  
  // pitch of voxel
  ofs.write( (char*)&sz, sizeof(int) );
  ofs.write( (char*)&dx, sizeof(float) );
  ofs.write( (char*)&dy, sizeof(float) );
  ofs.write( (char*)&dz, sizeof(float) );
  ofs.write( (char*)&sz, sizeof(int) );
  
  // type of stored data
  sz = sizeof(int)*1;
  int dtype = 0;
  dtype |= ( 0x1<<2 );  // medium ID
  ofs.write( (char*)&sz,  sizeof(int) );
  ofs.write( (char*)&dtype, sizeof(int) );
  ofs.write( (char*)&sz,  sizeof(int) );
  
  // medium ID
  sz = nx * sizeof(int);
  ofs.write( (char*)&sz, sizeof(int) );
  ofs.write( (char*)q,   sz );
  ofs.write( (char*)&sz, sizeof(int) );
  
  ofs.close();
  
  if (q) { delete [] q; q=NULL; }
}
