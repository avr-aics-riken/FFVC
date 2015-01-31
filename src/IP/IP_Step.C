//##################################################################################
//
// FFV-C : Frontflow / violet Cartesian
//
// Copyright (c) 2007-2011 VCAD System Research Program, RIKEN.
// All rights reserved.
//
// Copyright (c) 2011-2015 Institute of Industrial Science, The University of Tokyo.
// All rights reserved.
//
// Copyright (c) 2012-2015 Advanced Institute for Computational Science, RIKEN.
// All rights reserved.
//
//##################################################################################

/** 
 * @file   IP_Step.C
 * @brief  IP_Step class
 * @author aics
 */

#include "IP_Step.h"


// #################################################################
// パラメータをロード
bool IP_Step::getTP(Control* R, TextParser* tpCntl)
{
  std::string str;
  std::string label;
  REAL_TYPE ct;
  
  // 2D or 3D mode
  label = "/IntrinsicExample/Dimension";
  
  if ( !(tpCntl->getInspectedValue(label, str )) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    return false;
  }
  
  if     ( !strcasecmp(str.c_str(), "2d") )
  {
    mode = dim_2d;
  }
  else if( !strcasecmp(str.c_str(), "3d") )
  {
    mode = dim_3d;
  }
  else
  {
    Hostonly_ stamped_printf("\tParsing error : Invalid '%s'\n", label.c_str());
    return false;
  }
  
  
  // 媒質指定
  label = "/IntrinsicExample/FluidMedium";
  if ( !(tpCntl->getInspectedValue(label, str )) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    return false;
  }
  m_fluid = str;
  
  label = "/IntrinsicExample/SolidMedium";
  if ( !(tpCntl->getInspectedValue(label, str )) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    return false;
  }
  m_solid = str;
  
  
  
  // x-dir step
  label = "/IntrinsicExample/StepLength";
  if ( !(tpCntl->getInspectedValue(label, ct )) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    return false;
  }
  else
  {
	  width = ( R->Unit.Param == DIMENSIONAL ) ? ct : ct * RefL;
  }
  
  // z-dir step
  label = "/IntrinsicExample/StepHeight";
  if ( !(tpCntl->getInspectedValue(label, ct )) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    return false;
  }
  else
  {
	  height = ( R->Unit.Param == DIMENSIONAL ) ? ct : ct * RefL;
  }
  
  return true;
}


// #################################################################
// パラメータの表示
void IP_Step::printPara(FILE* fp, const Control* R)
{
  if ( !fp )
  {
    stamped_printf("\tFail to write into file\n");
    Exit(0);
  }
  
  fprintf(fp,"\n----------\n\n");
  fprintf(fp,"\n\t>> Intrinsic Backstep Class Parameters\n\n");
  
  fprintf(fp,"\tDimension Mode                     :  %s\n", (mode==dim_2d)?"2 Dimensional":"3 Dimensional");
  fprintf(fp,"\tStep Width (x-dir.)    [m] / [-]   : %12.5e / %12.5e\n", width, width/RefL);
  fprintf(fp,"\tStep Height(z-dir.)    [m] / [-]   : %12.5e / %12.5e\n", height, height/RefL);
  
}


// #################################################################
// 計算領域のセルIDとカット情報を設定する
void IP_Step::setup(int* bcd, Control* R, const int NoMedium, const MediumList* mat, long long* cut, int* bid)
{
  int mid_fluid;        /// 流体
  int mid_solid;        /// 固体

  // 流体
  if ( (mid_fluid = FBUtility::findIDfromLabel(mat, NoMedium, m_fluid)) == 0 )
  {
    Hostonly_ printf("\tLabel '%s' is not listed in MediumList\n", m_fluid.c_str());
    Exit(0);
  }
  
  // 固体
  if ( (mid_solid = FBUtility::findIDfromLabel(mat, NoMedium, m_solid)) == 0 )
  {
    Hostonly_ printf("\tLabel '%s' is not listed in MediumList\n", m_solid.c_str());
    Exit(0);
  }
  


  // ローカル
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
  REAL_TYPE dx = pitch[0];
  //REAL_TYPE dy = pitch[1];
  REAL_TYPE dz = pitch[2];
  
  // ローカルな無次元位置
  REAL_TYPE ox, oz;
  ox = origin[0];
  oz = origin[2];
  

  // step length 有次元値
  REAL_TYPE len = G_origin[0] + width/R->RefLength; // グローバルな無次元位置
  
  // step height 有次元値
  REAL_TYPE ht  = G_origin[2] + height/R->RefLength; // グローバルな無次元位置

  
#pragma omp parallel for firstprivate(ix, jx, kx, gd, mid_solid, ox, oz, dx, dz, len, ht) schedule(static)
  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        REAL_TYPE x = ox + 0.5*dx + dx*(i-1); // position of cell center
        REAL_TYPE z = oz + 0.5*dz + dz*(k-1); // position of cell center
        REAL_TYPE s = (len - x)/dx;
        
        if ( z <= ht )
        {
          if ( (x <= len) && (len < x+dx) )
          {
            setBit5(bid[m], mid_solid, X_plus);
            int r = quantize9(s);
            setCut9(cut[m], r, X_plus);
            
            size_t m1 = _F_IDX_S3D(i+1, j, k, ix, jx, kx, gd);
            setBit5(bid[m1], mid_solid, X_minus);
            int rr = quantize9(1.0-s);
            setCut9(cut[m1], rr, X_minus);
          }
          else if ( (x-dx < len) && (len < x) )
          {
            setBit5(bid[m], mid_solid, X_minus);
            int r = quantize9(-s);
            setCut9(cut[m], r, X_minus);
            
            size_t m1 = _F_IDX_S3D(i-1, j, k, ix, jx, kx, gd);
            setBit5(bid[m1], mid_solid, X_plus);
            int rr = quantize9(1.0+s);
            setCut9(cut[m1], rr, X_plus);
          }
        }

      }
    }
  }
  
#pragma omp parallel for firstprivate(ix, jx, kx, gd, mid_solid, ox, oz, dx, dz, len, ht) schedule(static)
  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        REAL_TYPE x = ox + 0.5*dx + dx*(i-1); // position of cell center
        REAL_TYPE z = oz + 0.5*dz + dz*(k-1); // position of cell center
        REAL_TYPE c = (ht - z)/dz;
        
        if ( x <= len )
        {
          if ( (z <= ht) && (ht < z+dz) )
          {
            setBit5(bid[m], mid_solid, Z_plus);
            int r = quantize9(c);
            setCut9(cut[m], r, Z_plus);

            size_t m1 = _F_IDX_S3D(i, j, k+1, ix, jx, kx, gd);
            setBit5(bid[m1], mid_solid, Z_minus);
            int rr = quantize9(1.0-c);
            setCut9(cut[m1], rr, Z_minus);
          }
          else if ( (z-dz < ht) && (ht < z) )
          {
            setBit5(bid[m], mid_solid, Z_minus);
            int r = quantize9(-c);
            setCut9(cut[m], r, Z_minus);

            size_t m1 = _F_IDX_S3D(i, j, k-1, ix, jx, kx, gd);
            setBit5(bid[m1], mid_solid, Z_plus);
            int rr = quantize9(1.0+c);
            setCut9(cut[m1], rr, Z_plus);
          }
        }
        
      }
    }
  }
  
  // ステップ部のガイドセルの設定
  
  // 隣接ランクのIDを取得 nID[6]
  const int* nID = paraMngr->GetNeighborRankID();
  
  if ( nID[X_minus] < 0 )
  {
    // デフォルトでガイドセルをSolidにする
#pragma omp parallel for firstprivate(ix, jx, kx, gd, mid_solid) schedule(static)
    for (int k=1; k<=kx; k++) {
      for (int j=1; j<=jx; j++) {
        
        // 媒質エントリ
        size_t m = _F_IDX_S3D(0, j, k, ix, jx, kx, gd);
        setBitID(bcd[m], mid_solid);
        
        // 交点
        size_t l = _F_IDX_S3D(1  , j  , k  , ix, jx, kx, gd);
        int r = quantize9(0.5);
        setCut9(cut[l], r, X_minus);
        
        // 境界ID
        setBit5(bid[l], mid_solid, X_minus);
      }
    }
    
    
    // Channel
#pragma omp parallel for firstprivate(ix, jx, kx, gd, mid_fluid, ht, oz, dz) schedule(static)
    for (int k=1; k<=kx; k++) {
      for (int j=1; j<=jx; j++) {
        
        REAL_TYPE z = oz + ( (REAL_TYPE)k-0.5 ) * dz;
        
        if ( z > ht )
        {
          size_t m = _F_IDX_S3D(0, j, k, ix, jx, kx, gd);
          setBitID(bcd[m], mid_fluid);
          
          size_t l = _F_IDX_S3D(1  , j  , k  , ix, jx, kx, gd);
          int r = quantize9(1.0);
          setCut9(cut[l], r, X_minus);
          setBit5(bid[l], 0, X_minus);
        }
      }
    }
    
  }
  
}
