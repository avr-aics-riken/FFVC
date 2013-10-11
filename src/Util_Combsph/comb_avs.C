//##################################################################################
//
// FFV-C : Frontflow / violet Cartesian
//
// Copyright (c) 2007-2011 VCAD System Research Program, RIKEN.
// All rights reserved.
//
// Copyright (c) 2011-2013 Institute of Industrial Science, The University of Tokyo.
// All rights reserved.
//
// Copyright (c) 2012-2013 Advanced Institute for Computational Science, RIKEN.
// All rights reserved.
//
//##################################################################################

/**
 * @file   comb_avs.C
 * @brief  COMB Class
 * @author kero
 */

#include "comb.h"

// #################################################################
// output avs file
void COMB::output_avs_header()
{

  if( myRank != 0 ) return; //myRank==0のときのみヘッダーレコードを出力

  std::string prefix;
  std::string fld_fname;
  std::string out_fname;
  FILE* fp;

  int ndim=3;
  int nspace=3;
  int dims[3];
  std::string dType;
  double min_ext[3],max_ext[3];
  double min_val[3],max_val[3];

  for(int i=0; i<ndfi; i++) {

    //ヘッダー出力に必用なcioクラスを取得しておく
    const cio_Domain* DFI_Domain  = dfi[i]->GetcioDomain();
    const cio_FileInfo* DFI_FInfo = dfi[i]->GetcioFileInfo();
    const cio_TimeSlice* TSlice   = dfi[i]->GetcioTimeSlice();

    //間引きを考慮しての計算空間サイズをセット
    dims[0]=(DFI_Domain->GlobalVoxel[0]+2*DFI_FInfo->GuideCell)/thin_count;
    dims[1]=(DFI_Domain->GlobalVoxel[1]+2*DFI_FInfo->GuideCell)/thin_count;
    dims[2]=(DFI_Domain->GlobalVoxel[2]+2*DFI_FInfo->GuideCell)/thin_count;
    if((DFI_Domain->GlobalVoxel[0]+2*DFI_FInfo->GuideCell)%thin_count != 0 ) dims[0]++;
    if((DFI_Domain->GlobalVoxel[1]+2*DFI_FInfo->GuideCell)%thin_count != 0 ) dims[1]++;
    if((DFI_Domain->GlobalVoxel[2]+2*DFI_FInfo->GuideCell)%thin_count != 0 ) dims[2]++;

    //データタイプのセット
    if( dfi[i]->GetDataType() == CIO::E_CIO_FLOAT32 ) 
    {
      dType = "float";
    } else if( dfi[i]->GetDataType() == CIO::E_CIO_FLOAT64 )
    {
      dType = "double";
    } else {
      dType = dfi[i]->GetDataTypeString();
      printf("\tillergal data type.(%s)\n",dType.c_str());
      Exit(0);
    }

    //座標値の最小値、最大値をセット
    for(int j=0; j<3; j++) {
      double pit=(DFI_Domain->GlobalRegion[j])/(double)(DFI_Domain->GlobalVoxel[j]);
      min_ext[j]=DFI_Domain->GlobalOrigin[j]+0.5*pit;
      max_ext[j]=(DFI_Domain->GlobalOrigin[j]+DFI_Domain->GlobalRegion[j])-0.5*pit;
    }

    prefix = DFI_FInfo->Prefix;

    fld_fname = Generate_FileName_Free(prefix,"fld",0,0,false);
    fld_fname = out_dirname + fld_fname;
    //出力ヘッダーファイルオープン
    if( (fp = fopen(fld_fname.c_str(),"w")) == NULL ) {
      printf("\tCan't open file.(%s)\n",fld_fname.c_str());
      Exit(0);
    }

    //先頭レコードの出力
    fprintf(fp,"# AVS field file\n");

    //計算空間の次元数を出力
    fprintf(fp,"ndim=%d\n",ndim);

    //計算空間サイズを出力
    fprintf(fp,"dim1=%d\n",dims[0]);
    fprintf(fp,"dim2=%d\n",dims[1]);
    fprintf(fp,"dim3=%d\n",dims[2]);

    //物理空間の次元数を出力
    fprintf(fp,"nspace=%d\n",nspace);

    //成分数の出力
    fprintf(fp,"veclen=%d\n",DFI_FInfo->Component);

    //データのタイプ出力
    fprintf(fp,"data=%s\n",dType.c_str());

    //座標定義情報の出力
    fprintf(fp,"field=uniform\n");

    //座標値の最小値、最大値出力
    fprintf(fp,"min_ext=%.6f %.6f %.6f\n",min_ext[0],min_ext[1],min_ext[2]);
    fprintf(fp,"max_ext=%.6f %.6f %.6f\n",max_ext[0],max_ext[1],max_ext[2]);

    //labelの出力
    for(int j=0; j<DFI_FInfo->Component; j++) {
     std::string label=dfi[i]->getComponentVariable(j);
     if( label == "" ) continue;
     fprintf(fp,"label=%s\n",label.c_str());
    }

    //step毎の出力
    if( TSlice->SliceList.size()>1 ) {
      fprintf(fp,"nstep=%d\n",TSlice->SliceList.size());
    }
    for(int j=0; j<TSlice->SliceList.size(); j++ ) {
      fprintf(fp,"time value=%.6f\n",TSlice->SliceList[j].time);
      for(int n=1; n<=DFI_FInfo->Component; n++) {
        int skip;
        if( dType == "float" ) {
          skip=96+(n-1)*4;
        } else {
          skip=140+(n-1)*8;
        }
        out_fname=Generate_FileName(prefix, TSlice->SliceList[j].step, 0, false);
        fprintf(fp,"variable %d file=%s filetype=binary skip=%d stride=%d\n",
                n,out_fname.c_str(),skip,DFI_FInfo->Component);
      }
      fprintf(fp,"EOT\n");
    }

    //出力ヘッダーファイルクローズ
    fclose(fp);

  }

  return;

}
