//##################################################################################
//
// Flow Base class
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

//@file   CompoFraction.C
//@brief  CompoFraction class
//@author aics

#include "CompoFraction.h"


// #################################################################
// コンポーネントの属するセルインデクスを求める
void CompoFraction::bbox_index(int* st, int* ed)
{
  find_index(st, box_min);
  find_index(ed, box_max);
}


// #################################################################
// 直方体領域のbboxを計算し、投影面積を計算
REAL_TYPE CompoFraction::bbox_rect_cylinder(Vec3r& mn, Vec3r& mx)
{
  Vec3r p[8], o, u, v, w;
  
  o = center;
  u = dir;
  w = nv;
  
  v = cross(w, u).normalize();
  
  // 直方体の8頂点の生成
  p[0] = o - 0.5f*width*u - 0.5f*height*v;
  p[1] = p[0] + depth *w;
  p[2] = p[1] + width *u;
  p[3] = p[0] + width *u;
  p[4] = p[0] + height*v;
  p[5] = p[4] + depth *w;
  p[6] = p[5] + width *u;
  p[7] = p[4] + width *u;
  
  for (int i=0; i<8; i++)
  {
    get_min(mn, p[i]);
    get_max(mx, p[i]);
  }
  
  // 投影面積
  REAL_TYPE a = (p[0]-p[3]).length();
  REAL_TYPE b = (p[0]-p[4]).length();

  return a*b;
}


// #################################################################
// 円筒領域のbboxを計算
// 標準位置で円周上の点をサンプリングし，逆変換後，min/max
REAL_TYPE CompoFraction::bbox_circ_cylinder(Vec3r& mn, Vec3r& mx)
{
  Vec3r r, q;
  int div_r = 1000; // 周上の分割数
  REAL_TYPE x, y;
  REAL_TYPE d;
  REAL_TYPE pi = 2.0*asin(1.0);
  REAL_TYPE dth = 2.0*pi/(REAL_TYPE)div_r; // radian
  
  for (int i=0; i<div_r; i++)
  {
    d = dth * (REAL_TYPE)i;
    x = r_fan * cos(d);
    y = r_fan * sin(d);
    
    // 表面
    q = rotate_inv(angle, r.assign(x, y, 0.0)) + center;
    get_min(mn, q);
    get_max(mx, q);
    
    // 裏面
    q = rotate_inv(angle, r.assign(x, y, depth)) + center;
    get_min(mn, q);
    get_max(mx, q);
  }
  
  // 投影面積
  REAL_TYPE a = (REAL_TYPE)(pi*(r_fan*r_fan-r_boss*r_boss));
  
  return a;
}


// #################################################################
// 点pの属するセルインデクスを求める
// Fortran index
void CompoFraction::find_index(int* w, const Vec3r p)
{
  Vec3r q = (p-org)/pch;
  
  w[0] = (int)ceil(q.x);
  w[1] = (int)ceil(q.y);
  w[2] = (int)ceil(q.z);
}


// #################################################################
// 指定法線nvがz軸の方向ベクトルに向かう回転角を計算する
// 回転の符号はz軸に向かう回転が右ねじ方向の場合を正にとる
void CompoFraction::get_angle()
{
  REAL_TYPE alpha, beta, c, d, c_alp, c_bta;
  REAL_TYPE eps = 1.0e-5, f_yz, f_xz;
  Vec3r p, q;
  Vec3r z(0.0, 0.0, 1.0);
  
  // 単位ベクトルnvがz軸の単位ベクトルと作る角度を返す
  // yz面への射影
  p.x = 0.0;
  p.y = nv.y;
  p.z = nv.z;
  c = p.length();
  
  if ( c != 0.0 )
  {
    c_alp = dot(z, p)/c;
    d = acos( c_alp );
    f_yz = c_alp+1.0;
    alpha = (nv.y >= 0.0) ? d : -d;
  }
  else
  {
    alpha = 0.0; // yz面への射影ベクトルがゼロの場合には回転しない
  }
  //printf("yz=%f : (%f %f %f)\n", alpha, p.x, p.y, p.z);
  
  // 参照ベクトルをalphaだけ回転して評価ベクトルを生成 > xz面への射影
  q.assign(alpha, 0.0, 0.0);
  p = rotate(q, nv);
  c = p.length();
  
  if ( c != 0.0 )
  {
    c_bta = dot(z, p)/c;
    d = acos( c_bta );
    f_xz = c_bta+1.0;
    beta = (nv.x >= 0.0) ? -d : d;
  }
  else
  {
    beta = 0.0;
  }
  
  // x軸とy軸の両方から見てz軸と反対の場合にだけ，y軸回りのみ回転する
  if ( (f_yz<eps) && (f_xz<eps) )
  {
    alpha = 0.0;
    beta  = acos(-1.0);
  }
  
  angle.assign(alpha, beta, 0.0);
  
  
  // 矩形の場合，単位ベクトルdirが回転した後，x軸の単位ベクトルへ回転する角度を計算
  if ( smode == mon_BOX )
  {
    REAL_TYPE c_gma, f_xy;
    Vec3r x(1.0, 0.0, 0.0);
    
    q = rotate(angle, dir); // 回転によりxy平面上に射影される > q.z=0.0
    c = q.length();
    
    if ( c != 0.0 )
    {
      c_gma = dot(x, q)/c;
      d = acos( c_gma );
      f_xy = c_gma+1.0;
      
      if ( f_xy<eps )
      {
        angle.z = 2.0*asin(1.0); // 反対方向なのでπ
      }
      else
      {
        angle.z = (q.y >= 0.0) ? -d : d;
      }
    }
    else
    {
      stamped_printf("\tInvalid Parameter of Heat exchanger : lateral vector is zero\n");
      Exit(0);
    }
  }
  
// ##########
#if 0
  stamped_printf("angle = (%f %f %f)\n", angle.x, angle.y, angle.z);
#endif
// ##########
  
}



// #################################################################
// 形状のbboxと投影面積を求める
REAL_TYPE CompoFraction::get_BboxArea()
{
  box_min.assign(1.0e6, 1.0e6, 1.0e6);
  box_max.assign(-1.0e6, -1.0e6, -1.0e6);
  REAL_TYPE a;
  
  if ( smode == mon_BOX )
  {
    a = bbox_rect_cylinder(box_min, box_max);
  }
  else
  {
    a = bbox_circ_cylinder(box_min, box_max);
  }

// ##########
#if 0
  stamped_printf("bbox min : %f %f %f\n", box_min.x, box_min.y, box_min.z);
  stamped_printf("bbox max : %f %f %f\n", box_max.x, box_max.y, box_max.z);
#endif
// ##########
  
  return a;
}


// #################################################################
// 矩形の形状パラメータをセットする
void CompoFraction::setShapeParam (const REAL_TYPE m_nv[3],
                                   const REAL_TYPE m_ctr[3],
                                   const REAL_TYPE m_dir[3],
                                   const REAL_TYPE m_depth,
                                   const REAL_TYPE m_width,
                                   const REAL_TYPE m_height)
{
  smode  = mon_BOX;
  nv     = m_nv;
  center = m_ctr;
  dir    = m_dir;
  depth  = m_depth;
  width  = m_width;
  height = m_height;
  
  nv.normalize();
  dir.normalize();
  
  if ( (nv.length() == 0.0) || (dir.length() == 0.0) )
  {
    stamped_printf("\tError : Invalid parameter of Heat Exchanger : zero vector\n");
    Exit(0);
  }
  
  if ( dot(nv, dir) != 0.0 )
  {
    stamped_printf("\tError : Invalid parameter of Heat Exchanger : non-orthogonal vectors\n");
    Exit(0);
  }
}


// #################################################################
// 円筒の形状パラメータをセットする
void CompoFraction::setShapeParam (const REAL_TYPE m_nv[3],
                                   const REAL_TYPE m_ctr[3],
                                   const REAL_TYPE m_depth,
                                   const REAL_TYPE m_r_fan,
                                   const REAL_TYPE m_r_boss)
{
  smode  = mon_CYLINDER;
  nv     = m_nv;
  center = m_ctr;
  depth  = m_depth;
  r_fan  = m_r_fan;
  r_boss = m_r_boss;
  
  nv.normalize();
  
  if ( nv.length() == 0.0 )
  {
    stamped_printf("\tError : Invalid parameter of Fan : zero vector\n");
    Exit(0);
  }
}


// #################################################################
// 体積率が(0,1)の間のセルに対してサブディビジョンを実施
void CompoFraction::subdivision(const int st[], const int ed[], REAL_TYPE* vf, double& flop)
{
  Vec3r base, b;
  Vec3r p, o;
  size_t m;
  REAL_TYPE c, h, ff, ph;
  int ix, jx, kx, gd, dv;
  REAL_TYPE r;
  
  // for optimization > variables defined outside
  ix = size[0];
  jx = size[1];
  kx = size[2];
  gd = guide;
  dv = division;
  h  = pch.x/(REAL_TYPE)dv;
  ff = 1.0/(REAL_TYPE)(dv*dv*dv);
  o  = org;
  ph = pch.x;
  
  if ( smode == mon_BOX )
  {
    for (int k=st[2]; k<=ed[2]; k++) {
      for (int j=st[1]; j<=ed[1]; j++) {
        for (int i=st[0]; i<=ed[0]; i++) {
          m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
          r = vf[m];
          
          if ( (r>0.0) && (r<1.0) ) {
            base.assign((REAL_TYPE)i-1.0, (REAL_TYPE)j-1.0, (REAL_TYPE)k-1.0);
            b = o + base * ph;
            
            c = 0.0;
            for (int k1=0; k1<dv; k1++) {
              p.z = b.z + ((REAL_TYPE)k1+0.5)*h;
              flop += 9.0;
              
              for (int j1=0; j1<dv; j1++) {
                p.y = b.y + ((REAL_TYPE)j1+0.5)*h;
                flop += 9.0;
                
                for (int i1=0; i1<dv; i1++) {
                  p.x = b.x + ((REAL_TYPE)i1+0.5)*h;
                  
                  c += judge_rect(p);
                  flop += 10.0 + 183.0;
                }
              }
            } // k1
            vf[m] = (REAL_TYPE)(c*ff);
            flop += 10.0;
          }
          
        }
      }
    }
  }
  else  // Cylinder
  {
    for (int k=st[2]; k<=ed[2]; k++) {
      for (int j=st[1]; j<=ed[1]; j++) {
        for (int i=st[0]; i<=ed[0]; i++) {
          m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
          r = vf[m];
          
          if ( (r>0.0) && (r<1.0) ) {
            base.assign((REAL_TYPE)i-1.0, (REAL_TYPE)j-1.0, (REAL_TYPE)k-1.0);
            base = o + base * ph;
            
            c = 0.0;
            for (int k1=0; k1<dv; k1++) {
              p.z = base.z + ((REAL_TYPE)k1+0.5)*h;
              flop += 9.0;
              
              for (int j1=0; j1<dv; j1++) {
                p.y = base.y + ((REAL_TYPE)j1+0.5)*h;
                flop += 9.0;
                
                for (int i1=0; i1<dv; i1++) {
                  p.x = base.x + ((REAL_TYPE)i1+0.5)*h;
                  
                  c += judge_cylinder(p);
                  flop += 10.0 + 186.0;
                }
              }
            } // k1
            vf[m] = (REAL_TYPE)(c*ff);
            flop += 10.0;
          }
          
        }
      }
    } 
    
  } // endif
  
} 



// #################################################################
// セルの8頂点の内外判定を行い，0, 1, otherに分類
// テスト候補のループ範囲（st[], ed[]）内で，テストセルの8頂点座標を生成し，形状の範囲内かどうかを判定する
// vfは加算するので、初期化しておく
void CompoFraction::vertex8(const int st[], const int ed[], REAL_TYPE* vf, double& flop)
{
  Vec3r base, o, b;
  Vec3r p[8];
  REAL_TYPE c, ph;
  
  // for optimization > variables defined outside
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  o  = org;
  ph = pch.x;
  
  if ( smode == mon_BOX )
  {
    // Rect cylinder
    for (int k=st[2]; k<=ed[2]; k++) {
      for (int j=st[1]; j<=ed[1]; j++) {
        for (int i=st[0]; i<=ed[0]; i++) {
          base.assign((REAL_TYPE)i-1.0, (REAL_TYPE)j-1.0, (REAL_TYPE)k-1.0);
          b    = o + base * ph;
          
          p[0] =          b;      // (0,0,0) 
          p[1] = shift_f1(b, ph); // (1,0,0)
          p[2] = shift_f2(b, ph); // (0,1,0)
          p[3] = shift_f3(b, ph); // (1,1,0)
          p[4] = shift_f4(b, ph); // (0,0,1)
          p[5] = shift_f5(b, ph); // (1,0,1)
          p[6] = shift_f6(b, ph); // (0,1,1)
          p[7] = shift_f7(b, ph); // (1,1,1)
          
          c = 0.0;
          
          for (int l=0; l<8; l++) {
            c += judge_rect(p[l]);
          }
          size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
          vf[m] += (REAL_TYPE)(c*0.125); // 1/8
        }
      }
    }
    flop += (double)( (ed[0]-st[0]+1)*(ed[1]-st[1]+1)*(ed[2]-st[2]+1) )*
    (3.0 + 2.0*3.0 + 6.0 + 8.0*183.0 + 2.0);
  }
  else
  {
    // Circular cylinder
    for (int k=st[2]; k<=ed[2]; k++) {
      for (int j=st[1]; j<=ed[1]; j++) {
        for (int i=st[0]; i<=ed[0]; i++) {
          base.assign((REAL_TYPE)i-1.0, (REAL_TYPE)j-1.0, (REAL_TYPE)k-1.0);
          b    = o + base * ph;
          
          p[0] =          b;      // (0,0,0) 
          p[1] = shift_f1(b, ph); // (1,0,0)
          p[2] = shift_f2(b, ph); // (0,1,0)
          p[3] = shift_f3(b, ph); // (1,1,0)
          p[4] = shift_f4(b, ph); // (0,0,1)
          p[5] = shift_f5(b, ph); // (1,0,1)
          p[6] = shift_f6(b, ph); // (0,1,1)
          p[7] = shift_f7(b, ph); // (1,1,1)
          
          c = 0.0;
          
          for (int l=0; l<8; l++) {
            c += judge_cylinder(p[l]);
          }
          size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
          vf[m] += (REAL_TYPE)(c*0.125); // 1/8
        }
      }
    }
    flop += (double)( (ed[0]-st[0]+1)*(ed[1]-st[1]+1)*(ed[2]-st[2]+1) )*
    (3.0 + 2.0*3.0 + 6.0 + 8.0*186.0 + 2.0);
  }
}


// #################################################################
// セルの8頂点の内外判定より50%以上のセルにIDを設定する
void ShapeMonitor::setID(const int st[], const int ed[], int* bcd, const int id)
{
  Vec3r base, o, b;
  Vec3r p[8];
  REAL_TYPE c, ph;
  
  // for optimization > variables defined outside
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  int odr = id;
  o  = org;
  ph = pch.x;
  
  if ( smode == mon_BOX )
  {
    // Rect cylinder
    for (int k=st[2]; k<=ed[2]; k++) {
      for (int j=st[1]; j<=ed[1]; j++) {
        for (int i=st[0]; i<=ed[0]; i++) {
          base.assign((REAL_TYPE)i-1.0, (REAL_TYPE)j-1.0, (REAL_TYPE)k-1.0);
          b    = o + base * ph;
          
          p[0] =          b;      // (0,0,0) 
          p[1] = shift_f1(b, ph); // (1,0,0)
          p[2] = shift_f2(b, ph); // (0,1,0)
          p[3] = shift_f3(b, ph); // (1,1,0)
          p[4] = shift_f4(b, ph); // (0,0,1)
          p[5] = shift_f5(b, ph); // (1,0,1)
          p[6] = shift_f6(b, ph); // (0,1,1)
          p[7] = shift_f7(b, ph); // (1,1,1)
          
          c = 0.0;
          
          for (int l=0; l<8; l++)
          {
            c += judge_rect(p[l]);
          }

          if ( c>=4.0)
          {
            size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
            setBitID(bcd[m], odr);
          }
        }
      }
    }

  }
  else
  {
    // Circular cylinder
    for (int k=st[2]; k<=ed[2]; k++) {
      for (int j=st[1]; j<=ed[1]; j++) {
        for (int i=st[0]; i<=ed[0]; i++) {
          base.assign((REAL_TYPE)i-1.0, (REAL_TYPE)j-1.0, (REAL_TYPE)k-1.0);
          b    = o + base * ph;
          
          p[0] =          b;      // (0,0,0) 
          p[1] = shift_f1(b, ph); // (1,0,0)
          p[2] = shift_f2(b, ph); // (0,1,0)
          p[3] = shift_f3(b, ph); // (1,1,0)
          p[4] = shift_f4(b, ph); // (0,0,1)
          p[5] = shift_f5(b, ph); // (1,0,1)
          p[6] = shift_f6(b, ph); // (0,1,1)
          p[7] = shift_f7(b, ph); // (1,1,1)
          
          c = 0.0;
          
          for (int l=0; l<8; l++)
          {
            c += judge_cylinder(p[l]);
          }
          
          if ( c>=4.0)
          {
            size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
            setBitID(bcd[m], odr);
          }
        }
      }
    }

  }
  /*
  for (int k=st[2]; k<=ed[2]; k++) {
    for (int j=st[1]; j<=ed[1]; j++) {
      for (int i=st[0]; i<=ed[0]; i++) {
        size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        int b = DECODE_CMP(bcd[m]);
        printf("%d ", b);
      }
    }
  }
  */
  
}
