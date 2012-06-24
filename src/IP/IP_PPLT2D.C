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
 * @file IP_PPLT2D.C
 * @brief IP_PPLT2D class
 * @author kero
 */

#include "IP_PPLT2D.h"

//@brief パラメータを取得する
bool IP_PPLT2D::getTP(Control* R, TPControl* tpCntl)
{
  std::string str;
  std::string label;
  
  // 媒質指定
  label="/Parameter/Intrinsic_Example/Fluid_medium";
  if ( !(tpCntl->GetValue(label, &str )) ) {
    Hostonly_ stamped_printf("\tParsing error : fail to get 'Fluid_medium' in 'Intrinsic_Example'\n");
    return false;
  }
  m_fluid = str;
  
  label="/Parameter/Intrinsic_Example/Solid_medium";
  if ( !(tpCntl->GetValue(label, &str )) ) {
    Hostonly_ stamped_printf("\tParsing error : fail to get 'Solid_medium' in 'Intrinsic_Example'\n");
    return false;
  }
  m_solid = str;
  
  return true;
}



// PPLT2Dの領域情報を設定する
void IP_PPLT2D::setDomain(Control* R, unsigned sz[3], REAL_TYPE org[3], REAL_TYPE wth[3], REAL_TYPE pch[3])
{
  // forced
  if (R->Unit.Param != NONDIMENSIONAL) {
    Hostonly_ printf("\tError : PPLT2D class is designed for only non-dimensional parameter\n");
    Exit(0);
  }
  
  // Z方向のチェック
  if ( sz[2] != 3 ) {
    Hostonly_ printf("\tError : The size of Z-direction must be 3.\n");
    Exit(0);
  }
  
  // 領域アスペクト比のチェック
  if ( sz[0] != sz[1]*2 ) {
    Hostonly_ printf("\tError : The number of division must be 2:1 (=X:Y)\n");
    Exit(0);
  }

  // 二次元の領域設定
  wth[0] =  2.0;
  wth[1] =  1.0;
  org[0] = -1.0;
  org[1] =  0.0;
  
  pch[0] = wth[0] / (REAL_TYPE)sz[0];
  pch[1] = wth[1] / (REAL_TYPE)sz[1];
  
  if ( pch[0] != pch[1] ) {
    Hostonly_ printf("\tVoxel width must be same between X and Y direction.\n");
    Exit(0);
  }
  
  // Z方向は3層に合わせて調整
  pch[2] =  pch[0];
  wth[2] =  pch[0]*3.0;
  org[2] = -pch[2]*1.5;
}


// PPLT2Dの計算領域のセルIDを設定する
void IP_PPLT2D::setup(int* mid, Control* R, REAL_TYPE* G_org, const int Nmax, MediumList* mat)
{
  int i,j,k;
  unsigned m;

  // Inner
  for (k=1; k<=(int)kmax; k++) {
    for (j=1; j<=(int)jmax; j++) {
      for (i=1; i<=(int)imax; i++) {
        m = FBUtility::getFindexS3D(size, guide, i, j, k);
        mid[m] = 1;
      }
    }
  }
}
