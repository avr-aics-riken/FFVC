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

//@fn void CompoFraction::get_fraction(int st[], int ed[], float* vf)
//@brief コンポーネントの体積率を計算
//@param st 開始インデクス
//@param ed 終了インデクス
//@param vf フラクション
void CompoFraction::get_fraction(int st[], int ed[], float* vf)
{
  get_angle();
  
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
  Vec3f base;
  Vec3f p;
  unsigned m;
  float c, r, h, ff;
  int ix, jx, kx, gc, dv;
  
  // for optimization > variables defined outside
  ix = (int)size[0];
  jx = (int)size[1];
  kx = (int)size[2];
  gc = (int)guide;
  dv = division;
  h  = pch.x/(float)dv;
  ff = 1.0/(float)(dv*dv);
  
  for (int k=st[2]; k<ed[2]; k++) {
    for (int j=st[1]; j<ed[1]; j++) {
      for (int i=st[0]; i<ed[0]; i++) {
        m = F_INDEX_S3D(ix, jx, kx, gc, i, j, k);
        r = vf[m];
        
        if ( (r>0.0) && (r<1.0) ) {
          base.assign((float)i-1.0, (float)j-1.0, (float)k-1.0);
          base = org +  base * pch;
          
          if ( smode == RECT ) {
            c = 0.0;
            for (int k1=0; k1<dv; k1++) {
              p.z = base.z + ((float)k1+0.5)*h;
              
              for (int j1=0; j1<dv; j1++) {
                p.y = base.y + ((float)j1+0.5)*h;
                
                for (int i1=0; i1<dv; i1++) {
                  p.x = base.x + ((float)i1+0.5)*h;
                  
                  c += judge_rect(p);
                }
              }
            } // k1
            vf[m] = c*ff; // 1/(dv*dv)
          }
          else { // Cylinder
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
        } // if sub-divide
        
      }
    }
  } // k
}

//@fn float CompoFraction::judge_cylinder(Vec3f p)
//@brief 円筒形状の内外判定
//@param p テスト点座標
//@ret 内部のときに1.0を返す
float CompoFraction::judge_cylinder(Vec3f p)
{
  Vec3f q = transform(angle, p);
  
  if (q.z > depth) return 0.0;
  float r = sqrtf(q.x*q.x + q.y*q.y);
  
  return (r<=radius) ? 1.0 : 0.0;
}

//@fn float CompoFraction::judge_rect(Vec3f p)
//@brief 矩形形状の内外判定
//@param p テスト点座標
//@ret 内部のときに1.0を返す
float CompoFraction::judge_rect(Vec3f p)
{
  Vec3f q = transform(angle, p);
  
  if (q.z > depth) return 0.0;
  if (q.x > fabs(0.5*width)) return 0.0;
  if (q.y > fabs(0.5*height)) return 0.0;
  
  return 1.0;
}

//@fn void CompoFraction::vertex8(int st[], int ed[], float* vf)
//@brief セルの8頂点の内外判定を行い，0, 1, otherに分類
//@param st 開始インデクス
//@param en 終了インデクス
//@param vf フラクション
void CompoFraction::vertex8(int st[], int ed[], float* vf)
{
  Vec3f base, o;
  Vec3f p[8];
  unsigned m;
  float c, h;
  int ix, jx, kx, gc;
  
  // for optimization > variables defined outside
  ix = (int)size[0];
  jx = (int)size[1];
  kx = (int)size[2];
  gc = (int)guide;
  o = org;
  h = pch.x;
  
  if ( smode == RECT ) {
    // Rect cylinder
    for (int k=st[2]; k<ed[2]; k++) {
      for (int j=st[1]; j<ed[1]; j++) {
        for (int i=st[0]; i<ed[0]; i++) {
          base.assign((float)i-1.0, (float)j-1.0, (float)k-1.0);
          base = o +    base * h;
          p[0] =        base;     // (0,0,0) 
          p[1] = shift_f1(base, h); // (1,0,0)
          p[2] = shift_f2(base, h); // (0,1,0)
          p[3] = shift_f3(base, h); // (1,1,0)
          p[4] = shift_f4(base, h); // (0,0,1)
          p[5] = shift_f5(base, h); // (1,0,1)
          p[6] = shift_f6(base, h); // (0,1,1)
          p[7] = shift_f7(base, h); // (1,1,1)
          
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
          base = o +    base * h;
          p[0] =        base;     // (0,0,0) 
          p[1] = shift_f1(base, h); // (1,0,0)
          p[2] = shift_f2(base, h); // (0,1,0)
          p[3] = shift_f3(base, h); // (1,1,0)
          p[4] = shift_f4(base, h); // (0,0,1)
          p[5] = shift_f5(base, h); // (1,0,1)
          p[6] = shift_f6(base, h); // (0,1,1)
          p[7] = shift_f7(base, h); // (1,1,1)
          
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
//@brief 回転角度を計算する
void CompoFraction::get_angle(void)
{
  float alpha, beta, c;
  Vec3f p;
  Vec3f x(1.0, 0.0, 0.0);
  Vec3f z(0.0, 0.0, 1.0);
  
  // 単位ベクトルnvがz軸の単位ベクトルと作る角度を返す
  // yz面への射影
  p.x = 0.0;
  p.y = nv.y;
  p.z = nv.z;
  
  c = p.length();
  beta  = acos( dot(nv, p)/c );
  alpha = acos( dot(z,  p)/c );
  
  angle.assign(alpha, beta, 0.0);
  
  // 矩形の場合，　単位ベクトルdirが回転した後，x軸の単位ベクトルと作る角度を返す
  if ( smode == RECT ) {

    Vec3f p = transform(angle, dir); // dirの回転
    angle.z = acos( dot(x, p)/p.length() );
  }
}

//@fn Vec3f CompoFraction::transform(const Vec3f p, const Vec3f u)
//@brief 角度p(alpha, beta, 0.0)でベクトルuを回転する
//@param p 回転角度
//@param u 方向ベクトル
//@ret 角度
Vec3f CompoFraction::transform(const Vec3f p, const Vec3f u)
{
  Vec3f a, b, c;
  
  // line vector expression
  a.x =  cos(p.y)*cos(p.z);
  a.y =  sin(p.x)*sin(p.y)*cos(p.z) - cos(p.x)*sin(p.z);
  a.z =  sin(p.x)*sin(p.z) + cos(p.x)*sin(p.y)*cos(p.z);
  
  b.x =  cos(p.y)*sin(p.z);
  b.y =  sin(p.x)*sin(p.y)*sin(p.z) - cos(p.x)*cos(p.z);
  b.z = -sin(p.x)*cos(p.z) + cos(p.x)*sin(p.y)*sin(p.z);
  
  c.x = -sin(p.y);
  c.y =  sin(p.x)*cos(p.y);
  c.z =  cos(p.x)*cos(p.y);
  
  return Vec3f( dot(a, u), dot(b, u), dot(c, u) );
}
