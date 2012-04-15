/*
 * SPHERE - Skeleton for PHysical and Engineering REsearch
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2004-2012
 *
 */

//@file CompoFraction.C
//@brief CompoFraction class
//@author keno, AICS, RIKEN

#include "CompoFraction.h"

//@fn void CompoFraction::bbox_index(int* st, int* ed)
//@brief コンポーネントの属するセルインデクスを求める
void CompoFraction::bbox_index(int* st, int* ed)
{
  find_index(st, box_min);
  find_index(ed, box_max);
}

//@fn void CompoFraction::find_index(int* w, const FB::Vec3f p)
//@brief 点pの属するセルインデクスを求める
//@note Fortran index
void CompoFraction::find_index(int* w, const FB::Vec3f p)
{
  FB::Vec3f q = (p-org)/pch;
  
  w[0] = (int)ceil(q.x);
  w[1] = (int)ceil(q.y);
  w[2] = (int)ceil(q.z);
}

//@fn void CompoFraction::bbox_rect_cylinder(FB::Vec3f& mn, FB::Vec3f& mx)
//@brief 直方体領域のbboxを計算し、投影面積を計算
float CompoFraction::bbox_rect_cylinder(FB::Vec3f& mn, FB::Vec3f& mx)
{
  FB::Vec3f p[8], o, u, v, w;
  
  o = center;
  u = dir;
  w = nv;
  
  v = cross(w, u).normalize();
  
  // 直方体の8頂点の生成
  p[0] = o - 0.5*width*u - 0.5*height*v;
  p[1] = p[0] + depth *w;
  p[2] = p[1] + width *u;
  p[3] = p[0] + width *u;
  p[4] = p[0] + height*v;
  p[5] = p[4] + depth *w;
  p[6] = p[5] + width *u;
  p[7] = p[4] + width *u;
  
  for (int i=0; i<8; i++) {
    get_min(mn, p[i]);
    get_max(mx, p[i]);
  }
  
  // 投影面積
  float a = (p[0]-p[3]).length();
  float b = (p[0]-p[4]).length();

  return a*b;
}

//@fn void CompoFraction::bbox_circ_cylinder(FB::Vec3f& mn, FB::Vec3f& mx)
//@brief 円筒領域のbboxを計算
//@note 標準位置で円周上の点をサンプリングし，逆変換後，min/max
float CompoFraction::bbox_circ_cylinder(FB::Vec3f& mn, FB::Vec3f& mx)
{
  FB::Vec3f r, q;
  int div_r = 1000; // 周上の分割数
  float x, y;
  double d;
  double pi = 2.0*asin(1.0);
  double dth = 2.0*pi/(double)div_r; // radian
  
  for (int i=0; i<div_r; i++) {
    d = dth * (double)i;
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
  float a = (float)(pi*(r_fan*r_fan-r_boss*r_boss));
  
  return a;
}

//@fn void CompoFraction::get_BboxArea(void)
//@brief 形状のbboxと投影面積を求める
float CompoFraction::get_BboxArea(void)
{
  box_min.assign(1.0e6, 1.0e6, 1.0e6);
  box_max.assign(-1.0e6, -1.0e6, -1.0e6);
  float a;
  
  if ( smode == SHAPE_BOX ) {
    a = bbox_rect_cylinder(box_min, box_max);
  }
  else {
    a = bbox_circ_cylinder(box_min, box_max);
  }
  
#if 0
  stamped_printf("bbox min : %f %f %f\n", box_min.x, box_min.y, box_min.z);
  stamped_printf("bbox max : %f %f %f\n", box_max.x, box_max.y, box_max.z);
#endif
  
  return a;
}


//@fn void CompoFraction::setShapeParam (const float m_nv[3], const float m_ctr[3], const float m_dir[3], const float m_depth, const float m_width, const float m_height)
//@brief 矩形の形状パラメータをセットする
void CompoFraction::setShapeParam (const float m_nv[3], const float m_ctr[3], const float m_dir[3], const float m_depth, const float m_width, const float m_height)
{
  smode  = SHAPE_BOX;
  nv     = m_nv;
  center = m_ctr;
  dir    = m_dir;
  depth  = m_depth;
  width  = m_width;
  height = m_height;
  
  nv.normalize();
  dir.normalize();
  if ( (nv.length() == 0.0) || (dir.length() == 0.0) ) {
    stamped_printf("\tError : Invalid parameter of Heat Exchanger : zero vector\n");
    Exit(0);
  }
  
  if ( dot(nv, dir) != 0.0 ) {
    stamped_printf("\tError : Invalid parameter of Heat Exchanger : non-orthogonal vectors\n");
    Exit(0);
  }
}

//@fn void CompoFraction::setShapeParam (const float m_nv[3], const float m_ctr[3], const float m_depth, const float m_r_fan, const float m_r_boss)
//@brief 円筒の形状パラメータをセットする
void CompoFraction::setShapeParam (const float m_nv[3], const float m_ctr[3], const float m_depth, const float m_r_fan, const float m_r_boss)
{
  smode  = SHAPE_CYLINDER;
  nv     = m_nv;
  center = m_ctr;
  depth  = m_depth;
  r_fan  = m_r_fan;
  r_boss = m_r_boss;
  
  nv.normalize();
  if ( nv.length() == 0.0 ) {
    stamped_printf("\tError : Invalid parameter of Fan : zero vector\n");
    Exit(0);
  }
}

//@fn void CompoFraction::subdivision(const int st[], const int ed[], float* vf, double& flop)
//@brief 体積率が(0,1)の間のセルに対してサブディビジョンを実施
//@param st 開始インデクス
//@param ed 終了インデクス
//@param vf フラクション
void CompoFraction::subdivision(const int st[], const int ed[], float* vf, double& flop)
{
  FB::Vec3f base, b;
  FB::Vec3f p, o;
  unsigned m;
  float c, r, h, ff, ph;
  int ix, jx, kx, gc, dv;
  
  // for optimization > variables defined outside
  ix = size[0];
  jx = size[1];
  kx = size[2];
  gc = guide;
  dv = division;
  h  = pch.x/(float)dv;
  ff = 1.0/(float)(dv*dv*dv);
  o  = org;
  ph = pch.x;
  
  if ( smode == SHAPE_BOX ) {
    
    for (int k=st[2]; k<=ed[2]; k++) {
      for (int j=st[1]; j<=ed[1]; j++) {
        for (int i=st[0]; i<=ed[0]; i++) {
          m = F_INDEX_S3D(ix, jx, kx, gc, i, j, k);
          r = vf[m];
          
          if ( (r>0.0) && (r<1.0) ) {
            base.assign((float)i-1.0, (float)j-1.0, (float)k-1.0);
            b = o + base * ph;
            
            c = 0.0;
            for (int k1=0; k1<dv; k1++) {
              p.z = b.z + ((float)k1+0.5)*h;
              flop += 9.0;
              
              for (int j1=0; j1<dv; j1++) {
                p.y = b.y + ((float)j1+0.5)*h;
                flop += 9.0;
                
                for (int i1=0; i1<dv; i1++) {
                  p.x = b.x + ((float)i1+0.5)*h;
                  
                  c += judge_rect(p);
                  flop += 10.0 + 183.0;
                }
              }
            } // k1
            vf[m] = c*ff;
            flop += 10.0;
          }
          
        }
      }
    }
  }
  else { // Cylinder
    
    for (int k=st[2]; k<=ed[2]; k++) {
      for (int j=st[1]; j<=ed[1]; j++) {
        for (int i=st[0]; i<=ed[0]; i++) {
          m = F_INDEX_S3D(ix, jx, kx, gc, i, j, k);
          r = vf[m];
          
          if ( (r>0.0) && (r<1.0) ) {
            base.assign((float)i-1.0, (float)j-1.0, (float)k-1.0);
            base = o + base * ph;
            
            c = 0.0;
            for (int k1=0; k1<dv; k1++) {
              p.z = base.z + ((float)k1+0.5)*h;
              flop += 9.0;
              
              for (int j1=0; j1<dv; j1++) {
                p.y = base.y + ((float)j1+0.5)*h;
                flop += 9.0;
                
                for (int i1=0; i1<dv; i1++) {
                  p.x = base.x + ((float)i1+0.5)*h;
                  
                  c += judge_cylinder(p);
                  flop += 10.0 + 186.0;
                }
              }
            } // k1
            vf[m] = c*ff;
            flop += 10.0;
          }
          
        }
      }
    } 
    
  } // endif
  
} 

//@fn void CompoFraction::vertex8(const int st[], const int ed[], float* vf, double& flop)
//@brief セルの8頂点の内外判定を行い，0, 1, otherに分類
//@param st 開始インデクス
//@param en 終了インデクス
//@param vf フラクション
//@note テスト候補のループ範囲（st[], ed[]）内で，テストセルの8頂点座標を生成し，形状の範囲内かどうかを判定する
//@note vfは加算するので、初期化しておく
void CompoFraction::vertex8(const int st[], const int ed[], float* vf, double& flop)
{
  FB::Vec3f base, o, b;
  FB::Vec3f p[8];
  unsigned m;
  float c, ph;
  int ix, jx, kx, gc;
  
  // for optimization > variables defined outside
  ix = size[0];
  jx = size[1];
  kx = size[2];
  gc = guide;
  o  = org;
  ph = pch.x;
  
  if ( smode == SHAPE_BOX ) {
    // Rect cylinder
    for (int k=st[2]; k<=ed[2]; k++) {
      for (int j=st[1]; j<=ed[1]; j++) {
        for (int i=st[0]; i<=ed[0]; i++) {
          base.assign((float)i-1.0, (float)j-1.0, (float)k-1.0);
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
          m = F_INDEX_S3D(ix, jx, kx, gc, i, j, k);
          vf[m] += c*0.125; // 1/8
        }
      }
    }
    flop += (double)( (ed[0]-st[0]+1)*(ed[1]-st[1]+1)*(ed[2]-st[2]+1) )*
            (3.0 + 2.0*3.0 + 6.0 + 8.0*183.0 + 2.0);
  }
  else {
    // Circular cylinder
    for (int k=st[2]; k<=ed[2]; k++) {
      for (int j=st[1]; j<=ed[1]; j++) {
        for (int i=st[0]; i<=ed[0]; i++) {
          base.assign((float)i-1.0, (float)j-1.0, (float)k-1.0);
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
          m = F_INDEX_S3D(ix, jx, kx, gc, i, j, k);
          vf[m] += c*0.125; // 1/8
        }
      }
    }
    flop += (double)( (ed[0]-st[0]+1)*(ed[1]-st[1]+1)*(ed[2]-st[2]+1) )*
            (3.0 + 2.0*3.0 + 6.0 + 8.0*186.0 + 2.0);
  }
}

//@fn void CompoFraction::get_angle(void)
//@brief 指定法線nvがz軸の方向ベクトルに向かう回転角を計算する
//@note 回転の符号はz軸に向かう回転が右ねじ方向の場合を正にとる
void CompoFraction::get_angle(void)
{
  float alpha, beta, c, d, c_alp, c_bta;
  float eps = 1.0e-5, f_yz, f_xz;
  FB::Vec3f p, q;
  FB::Vec3f z(0.0, 0.0, 1.0);
  
  // 単位ベクトルnvがz軸の単位ベクトルと作る角度を返す
  // yz面への射影
  p.x = 0.0;
  p.y = nv.y;
  p.z = nv.z;
  c = p.length();
  
  if ( c != 0.0 ) {
    c_alp = dot(z, p)/c;
    d = acos( c_alp );
    f_yz = c_alp+1.0;
    alpha = (nv.y >= 0.0) ? d : -d;
  }
  else {
    alpha = 0.0; // yz面への射影ベクトルがゼロの場合には回転しない
  }
  
  // 参照ベクトルをalphaだけ回転して評価ベクトルを生成 > xz面への射影
  q.assign(alpha, 0.0, 0.0);
  p = rotate(q, nv);
  c = p.length();
  
  if ( c != 0.0 ) {
    c_bta = dot(z, p)/c;
    d = acos( c_bta );
    f_xz = c_bta+1.0;
    beta = (nv.x >= 0.0) ? -d : d;
  }
  else {
    beta = 0.0;
  }
  
  // x軸とy軸の両方から見てz軸と反対の場合にだけ，y軸回りのみ回転する
  if ( (f_yz<eps) && (f_xz<eps) ) {
    alpha = 0.0;
    beta  = acos(-1.0);
  }
  
  angle.assign(alpha, beta, 0.0);
  
  
  // 矩形の場合，単位ベクトルdirが回転した後，x軸の単位ベクトルへ回転する角度を計算
  if ( smode == SHAPE_BOX ) {
    float c_gma, f_xy;
    FB::Vec3f x(1.0, 0.0, 0.0);
    
    q = rotate(angle, dir); // 回転によりxy平面上に射影される > q.z=0.0
    c = q.length();
    
    if ( c != 0.0 ) {
      c_gma = dot(x, q)/c;
      d = acos( c_gma );
      f_xy = c_gma+1.0;
      if ( f_xy<eps ) {
        angle.z = 2.0*asin(1.0); // 反対方向なのでπ
      }
      else {
        angle.z = (q.y >= 0.0) ? -d : d;
      }
    }
    else {
      stamped_printf("\tInvalid Parameter of Heat exchanger : lateral vector is zero\n");
      Exit(0);
    }
  }
#if 0
  stamped_printf("angle = (%f %f %f)\n", angle.x, angle.y, angle.z);
#endif
  
}

//@fn void ShapeMonitor::setID(const int st[], const int ed[], int* mid, int id)
//@brief セルの8頂点の内外判定より50%以上のセルにIDを設定する
//@param st 開始インデクス
//@param en 終了インデクス
//@param mid
//@pram id
void ShapeMonitor::setID(const int st[], const int ed[], int* mid, int id)
{
  FB::Vec3f base, o, b;
  FB::Vec3f p[8];
  unsigned m;
  float c, ph;
  int ix, jx, kx, gc;
  
  // for optimization > variables defined outside
  ix = size[0];
  jx = size[1];
  kx = size[2];
  gc = guide;
  o  = org;
  ph = pch.x;
  
  if ( smode == SHAPE_BOX ) {
    // Rect cylinder
    for (int k=st[2]; k<=ed[2]; k++) {
      for (int j=st[1]; j<=ed[1]; j++) {
        for (int i=st[0]; i<=ed[0]; i++) {
          base.assign((float)i-1.0, (float)j-1.0, (float)k-1.0);
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
          if ( c>=4.0) {
            m = F_INDEX_S3D(ix, jx, kx, gc, i, j, k);
            mid[m] = id;
          }
        }
      }
    }

  }
  else {
    // Circular cylinder
    for (int k=st[2]; k<=ed[2]; k++) {
      for (int j=st[1]; j<=ed[1]; j++) {
        for (int i=st[0]; i<=ed[0]; i++) {
          base.assign((float)i-1.0, (float)j-1.0, (float)k-1.0);
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
          if ( c>=4.0) {
            m = F_INDEX_S3D(ix, jx, kx, gc, i, j, k);
            mid[m] = id;
          }
        }
      }
    }

  }
}
