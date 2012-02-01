/*
 * SPHERE - Skeleton for PHysical and Engineering REsearch
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2004-2010
 *
 */

//@file Intrinsic.C
//@brief FlowBase Intrinsic class
//@author keno, FSI Team, VCAD, RIKEN

#include "Intrinsic.h"

/**
 @fn void Intrinsic::printExample(FILE* fp, const char* str)
 @brief 例題名称の表示
 @param fp 出力ファイルのファイルポインタ
 @param str 表示文字列
 */
void Intrinsic::printExample(FILE* fp, const char* str)
{
  if( !fp ) {
    stamped_printf("\tFail to write into file\n");
    assert(0);
  }
  fprintf(fp,"\n\tExample : %s\n\n", str); fflush(fp);

  printf("\n\tExample : %s\n\n", str);
}

/**
 @fn void Intrinsic::setControlVars(Control* R)
 @brief 基本情報のコピー
 @param R コントロールクラスのポインタ
 */
void Intrinsic::setControlVars(Control* R)
{
  guide          = R->guide;
  imax = size[0] = R->imax;
  jmax = size[1] = R->jmax;
  kmax = size[2] = R->kmax;
  RefL = R->RefLength;
}

/**
 @fn void Intrinsic::printParaInfo(FILE* mp, FILE* fp, Control* R)
 @brief パラメータの表示
 @param mp 標準出力のファイルポインタ
 @param fp ファイルポインタ
 @param R コントロールクラスのポインタ
 */
void Intrinsic::printParaInfo(FILE* mp, FILE* fp, Control* R)
{
  printPara(mp, R);
  printPara(fp, R);
}

/**
 @fn void Intrinsic::printPara(FILE* fp, Control* R)
 @brief パラメータの表示
 @param fp ファイルポインタ
 @param R コントロールクラスのポインタ
 */
void Intrinsic::printPara(FILE* fp, Control* R)
{
  if ( !fp ) {
    stamped_printf("\tFail to write into file\n");
    assert(0);
  }
  // describe method
}

/**
 @fn void Intrinsic::genVFfromBcx(SKL_REAL* VF, unsigned* bx)
 @brief BCindexから体積率を計算する．Fluid=1.0, Solid=0.0
 @param VF 体積占有率（マスク）
 @param bx BCindex
 */
void Intrinsic::genVFfromBcx(SKL_REAL* VF, unsigned* bx)
{
  int i,j,k;
  unsigned m;

  for (k=0; k<=(int)(kmax+1); k++) {
    for (j=0; j<=(int)(jmax+1); j++) {
      for (i=0; i<=(int)(imax+1); i++) {
        m = SklUtil::getFindexS3D(size, guide, i, j, k);
        VF[m] = GET_SHIFT_F( bx[m], STATE_BIT );
      }
    }
  }
}

/**
 @fn void Intrinsic::writeSVX(SKL_REAL *vf, int *id, Control* R)
 @brief 例題のモデルをsvxフォーマットで出力する(体積率とID)
 @param vf 体積占有率
 @param id ID情報
 @param R コントロールクラスのポインタ
 */
void Intrinsic::writeSVX(SKL_REAL *vf, int *id, Control* R)
{
  SklParaManager* para_mng = ParaCmpo->GetParaManager();

  int      nx, sz, i,j,k;
  unsigned ix, jx, kx, m, l;
  float    ox, oy, oz, dx, dy, dz;

  char svx_fname[512];

  if( para_mng->IsParallel() ){
    sprintf( svx_fname, "example_%06d.svx", pn.ID );
  } else {
    sprintf( svx_fname, "example.svx" );
  }
  ofstream ofs(svx_fname, ios::out | ios::binary);
  if (!ofs) {
    cout << "\tCan't open " << svx_fname << " file" << endl;
    assert(0);
  }

  ix = imax+2;  // +2 means guide cell for IP model
  jx = jmax+2;
  kx = kmax+2;
  nx = (int)ix*jx*kx;
  dx = (float)R->dx[0]*RefL;
  dy = (float)R->dx[1]*RefL;
  dz = (float)R->dx[2]*RefL;
  ox = (float)R->org[0]*RefL - dx; // 片側1層分をシフト
  oy = (float)R->org[1]*RefL - dy;
  oz = (float)R->org[2]*RefL - dz;
  
  //stamped_printf("example out org(%e %e %e) dimensional\n", ox, oy, oz);

  float *f = new float[nx];
  int   *q = new int[nx];

  for (k=0; k<=(int)(kmax+1); k++) {
    for (j=0; j<=(int)(jmax+1); j++) {
      for (i=0; i<=(int)(imax+1); i++) {
        l = (unsigned)(ix*jx*k + ix*j + i);
        m = SklUtil::getFindexS3D(size, guide, i, j, k);
        q[l] = id[m];
        f[l] = (float)vf[m];
      }
    }
  }

  // voxel size
  sz = sizeof(unsigned)*3;
  ofs.write( (char*)&sz, sizeof(int) );
  ofs.write( (char*)&ix, sizeof(unsigned) );
  ofs.write( (char*)&jx, sizeof(unsigned) );
  ofs.write( (char*)&kx, sizeof(unsigned) );
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

/**
 @fn void Intrinsic::writeSVX(int *id, Control* R)
 @brief 例題のモデルをsvxフォーマットで出力する(ID)
 @param id ID情報
 @param R コントロールクラスのポインタ
 */
void Intrinsic::writeSVX(int *id, Control* R)
{
  SklParaManager* para_mng = ParaCmpo->GetParaManager();
  
  int      nx, sz, i,j,k;
  unsigned ix, jx, kx, m, l;
  float    ox, oy, oz, dx, dy, dz;
  
  char svx_fname[512];
  
  if( para_mng->IsParallel() ){
    sprintf( svx_fname, "example_%06d.svx", pn.ID );
  } else {
    sprintf( svx_fname, "example.svx" );
  }
  ofstream ofs(svx_fname, ios::out | ios::binary);
  if (!ofs) {
    cout << "\tCan't open " << svx_fname << " file" << endl;
    assert(0);
  }
  
  ix = imax+2;  // +2 means guide cell for IP model
  jx = jmax+2;
  kx = kmax+2;
  nx = (int)ix*jx*kx;
  dx = (float)R->dx[0]*RefL;
  dy = (float)R->dx[1]*RefL;
  dz = (float)R->dx[2]*RefL;
  ox = (float)R->org[0]*RefL - dx; // 片側1層分をシフト
  oy = (float)R->org[1]*RefL - dy;
  oz = (float)R->org[2]*RefL - dz;
  
  //stamped_printf("example out org(%e %e %e) dimensional\n", ox, oy, oz);
  
  int   *q = new int[nx];
  
  for (k=0; k<=(int)(kmax+1); k++) {
    for (j=0; j<=(int)(jmax+1); j++) {
      for (i=0; i<=(int)(imax+1); i++) {
        l = (unsigned)(ix*jx*k + ix*j + i);
        m = SklUtil::getFindexS3D(size, guide, i, j, k);
        q[l] = id[m];
      }
    }
  }
  
  // voxel size
  sz = sizeof(unsigned)*3;
  ofs.write( (char*)&sz, sizeof(int) );
  ofs.write( (char*)&ix, sizeof(unsigned) );
  ofs.write( (char*)&jx, sizeof(unsigned) );
  ofs.write( (char*)&kx, sizeof(unsigned) );
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
