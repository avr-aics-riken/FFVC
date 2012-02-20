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


//@fn Vec3f CompoFraction::get_alpha_beta(const Vec3f w)
//@brief 単位ベクトルwがz軸の単位ベクトルと作る角度を返す
//@param w z軸にあわせるベクトル
//@ret 回転角度
Vec3f CompoFraction::get_alpha_beta(const Vec3f& w)
{
  float alpha, beta, c;
  Vec3f p;
  Vec3f z(0.0, 0.0, 1.0);
  
  // yz面への射影
  p.x = 0.0;
  p.y = w.y;
  p.z = w.z;
  
  c = p.length();
  beta  = acos( vec3f_dot(&w, &p)/c );
  alpha = acos( vec3f_dot(&z, &p)/c );
  
  return Vec3f(alpha, beta, 0.0);
}

//@fn float CompoFraction::get_gamma(const Vec3f u, const Vec3f angle)
//@brief 単位ベクトルuが回転した後，x軸の単位ベクトルと作る角度を返す
//@param u 方向ベクトル
//@param angle 回転角度
//@ret 回転角度
float CompoFraction::get_gamma(const Vec3f& u, const Vec3f& angle)
{
  float gamma;
  Vec3f p;
  Vec3f x(1.0, 0.0, 0.0);
  
  // uの射影
  p = transform(angle, u);
  
  gamma = acos( vec3f_dot(&x, &p)/p.length() );
  
  return gamma;
}

//@fn Vec3f CompoFraction::transform(const Vec3f angle, const Vec3f u)
//@brief 角度angle(alpha, beta, 0.0)でベクトルuを回転する
//@param angle 回転角度
//@param u 方向ベクトル
//@ret 角度
Vec3f CompoFraction::transform(const Vec3f angle, const Vec3f u)
{
  Vec3f a, b, c, *p;
  
  p = &angle;
  
  a.x =  cos(p.y)*cos(p.z);
  a.y =  cos(p.y)*sin(p.z);
  a.z = -sin(p.y);
  
  b.x = sin(p.x)*sin(p.y)*cos(p.z) - cos(p.x)*sin(p.z);
  b.y = sin(p.x)*sin(p.y)*sin(p.z) - cos(p.x)*cos(p.z);
  b.z = sin(p.x)*cos(p.y);
  
  c.x =  sin(p.x)*sin(p.z) + cos(p.x)*sin(p.y)*cos(p.z);
  c.y = -sin(p.x)*cos(p.z) + cos(p.x)*sin(p.y)*sin(p.z);
  c.z =  cos(p.x)*cos(p.y);
  
  return Vec3f(
  a.x * u.x + b.x * u.y + c.x * u.z, 
  a.y * u.x + b.y * u.y + c.y * u.z, 
  a.z * u.x + b.z * u.y + c.z * u.z);
}
