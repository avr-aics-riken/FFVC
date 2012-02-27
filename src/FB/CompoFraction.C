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

//@fn void CompoFraction::setShapeParam (FB::Vec3f m_nv, FB::Vec3f m_ctr, FB::Vec3f m_dir, float m_depth, float m_width, float m_height)
//@brief 矩形の形状パラメータをセットする
void CompoFraction::setShapeParam (FB::Vec3f m_nv, FB::Vec3f m_ctr, FB::Vec3f m_dir, float m_depth, float m_width, float m_height)
{
  smode  = RECT_CYL;
  nv     = m_nv;
  center = m_ctr;
  dir    = m_dir;
  depth  = m_depth;
  width  = m_width;
  height = m_height;
}

//@fn void CompoFraction::setShapeParam (FB::Vec3f m_nv, FB::Vec3f m_ctr, float m_depth, float m_radius)
//@brief 円筒の形状パラメータをセットする
void CompoFraction::setShapeParam (FB::Vec3f m_nv, FB::Vec3f m_ctr, float m_depth, float m_radius)
{
  smode  = CIRC_CYL;
  nv     = m_nv;
  center = m_ctr;
  depth  = m_depth;
  radius = m_radius;
}

//@fn void CompoFraction::get_fraction(int st[], int ed[], float* vf)
//@brief コンポーネントの体積率を計算
//@param st 開始インデクス
//@param ed 終了インデクス
//@param vf フラクション
void CompoFraction::get_fraction(int st[], int ed[], float* vf)
{
  get_angle(); printf("angle = (%f %f %f)\n", angle.x, angle.y, angle.z);
  
  vertex8(st, ed, vf);
  
  subdivision(st, ed, vf);
}

//@fn void CompoFraction::subdivision(int st[], int ed[], float* vf)
//@brief 体積率が(0,1)の間のセルに対してサブディビジョンを実施
//@param st 開始インデクス
//@param ed 終了インデクス
//@param vf フラクション
void CompoFraction::subdivision(int st[], int ed[], float* vf)
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
  
  if ( smode == RECT_CYL ) {
    
    for (int k=st[2]; k<ed[2]; k++) {
      for (int j=st[1]; j<ed[1]; j++) {
        for (int i=st[0]; i<ed[0]; i++) {
          m = F_INDEX_S3D(ix, jx, kx, gc, i, j, k);
          r = vf[m];
          
          if ( (r>0.0) && (r<1.0) ) {
            base.assign((float)i-1.0, (float)j-1.0, (float)k-1.0);
            b = o + base * ph;
            
            c = 0.0;
            for (int k1=0; k1<dv; k1++) {
              p.z = b.z + ((float)k1+0.5)*h;
              
              for (int j1=0; j1<dv; j1++) {
                p.y = b.y + ((float)j1+0.5)*h;
                
                for (int i1=0; i1<dv; i1++) {
                  p.x = b.x + ((float)i1+0.5)*h;
                  
                  c += judge_rect(p);
                }
              }
            } // k1
            vf[m] = c*ff; // 1/(dv*dv)
          }
          
        }
      }
    }
    
  }
  else { // Cylinder
    
    for (int k=st[2]; k<ed[2]; k++) {
      for (int j=st[1]; j<ed[1]; j++) {
        for (int i=st[0]; i<ed[0]; i++) {
          m = F_INDEX_S3D(ix, jx, kx, gc, i, j, k);
          r = vf[m];
          
          if ( (r>0.0) && (r<1.0) ) {
            base.assign((float)i-1.0, (float)j-1.0, (float)k-1.0);
            base = o + base * ph;
            
            c = 0.0;
            for (int k1=0; k1<dv; k1++) {
              p.z = base.z + ((float)k1+0.5)*h;
              
              for (int j1=0; j1<dv; j1++) {
                p.y = base.y + ((float)j1+0.5)*h;
                
                for (int i1=0; i1<dv; i1++) {
                  p.x = base.x + ((float)i1+0.5)*h;
                  
                  c += judge_cylinder(p);
                }
              }
            } // k1
            vf[m] = c*ff; // 1/(dv*dv)
          }
          
        }
      }
    } 
    
  } // endif
  
} 

//@fn void CompoFraction::vertex8(int st[], int ed[], float* vf)
//@brief セルの8頂点の内外判定を行い，0, 1, otherに分類
//@param st 開始インデクス
//@param en 終了インデクス
//@param vf フラクション
//@note テスト候補のループ範囲（st[], ed[]）内で，テストセルの8頂点座標を生成し，形状の範囲内かどうかを判定する
void CompoFraction::vertex8(int st[], int ed[], float* vf)
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
  
  if ( smode == RECT_CYL ) {
    // Rect cylinder
    for (int k=st[2]; k<ed[2]; k++) {
      for (int j=st[1]; j<ed[1]; j++) {
        for (int i=st[0]; i<ed[0]; i++) {
          base.assign((float)i-1.0, (float)j-1.0, (float)k-1.0);
          b    = o    + base * ph;
          
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
          vf[m] = c*0.125; // 1/8
        }
      }
    }
    
  }
  else {
    // Circular cylinder
    for (int k=st[2]; k<ed[2]; k++) {
      for (int j=st[1]; j<ed[1]; j++) {
        for (int i=st[0]; i<ed[0]; i++) {
          base.assign((float)i-1.0, (float)j-1.0, (float)k-1.0);
          b    = o +    base * ph;

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
          vf[m] = c*0.125; // 1/8
        }
      }
    }
    
  }
}

//@fn void CompoFraction::get_angle(void)
//@brief 指定法線nvがz軸の方向ベクトルに向かう回転角を計算する
//@note 回転の符号はz軸に向かう回転が右ねじ方向の場合を正にとる
void CompoFraction::get_angle(void)
{
  float alpha, beta, c, d, d_yz, d_xz;
  float eps = 1.0e-5, f_yz, f_xz;
  FB::Vec3f p;
  FB::Vec3f x(1.0, 0.0, 0.0);
  FB::Vec3f z(0.0, 0.0, 1.0);
  
  // 単位ベクトルnvがz軸の単位ベクトルと作る角度を返す
  // yz面への射影
  p.x = 0.0;
  p.y = nv.y;
  p.z = nv.z;
  c = p.length();
  d_yz = (c==0.0) ? 0.0 : dot(z, p)/c;
  d = acos( d_yz );
  f_yz = fabs(d_yz+1.0);
  if ( f_yz<eps ) { // x軸から見てnvとzが反対方向の場合には，xz面の回転を採用する．判定は誤差を考慮
    alpha = 0.0;
  }
  else {
    alpha = (nv.y >= 0.0) ? d : -d;
  }
  
  // xz面への射影
  p.x = nv.x;
  p.y = 0.0;
  p.z = nv.z;
  c = p.length();
  d_xz = (c==0.0) ? 0.0 : dot(z, p)/c;
  d = acos( d_xz );
  f_xz = fabs(d_xz+1.0);
  if ( f_xz<eps ) { // y軸から見てnvとzが反対方向の場合には，d_yzで判断
    beta = (f_yz<eps) ? d : 0.0; // 既にx軸から見て反対方向である場合にはπとする
  }
  else {
    beta = (nv.x >= 0.0) ? -d : d;
  }
  
  
  angle.assign(alpha, beta, 0.0);
  printf("A: %f %f %f\n", angle.x, angle.y, angle.z);
  
  // 矩形の場合，　単位ベクトルdirが回転した後，x軸の単位ベクトルと作る角度を返す
  if ( smode == RECT_CYL ) {
    float d_xy;
    FB::Vec3f q = transform(angle, dir); // 回転によりxy平面上に射影される
    c = q.length();
    d_xy = (c==0.0) ? 0.0 : dot(x, q)/c;
    d = acos( d_xy );
    angle.z = (q.y >= 0.0) ? -d : d;
  }
}

//@fn FB::Vec3f CompoFraction::transform(const Vec3f p, const Vec3f u)
//@brief 回転ベクトルp(alpha, beta, gamma)でベクトルuを回転する
//@param p 回転角度
//@param u 方向ベクトル
//@ret 角度
FB::Vec3f CompoFraction::transform(const FB::Vec3f p, const FB::Vec3f u)
{
  FB::Vec3f a, b, c;
  
  // line vector expression
  a.x =  cos(p.y)*cos(p.z);
  a.y =  sin(p.x)*sin(p.y)*cos(p.z) - cos(p.x)*sin(p.z);
  a.z =  cos(p.x)*sin(p.y)*cos(p.z) + sin(p.x)*sin(p.z);
  
  b.x =  cos(p.y)*sin(p.z);
  b.y =  sin(p.x)*sin(p.y)*sin(p.z) + cos(p.x)*cos(p.z);
  b.z =  cos(p.x)*sin(p.y)*sin(p.z) - sin(p.x)*cos(p.z);
  
  c.x = -sin(p.y);
  c.y =  sin(p.x)*cos(p.y);
  c.z =  cos(p.x)*cos(p.y);
  
  return FB::Vec3f( dot(a, u), dot(b, u), dot(c, u) );
}
