//##################################################################################
//
// Flow Base class
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

//@file   CompoFraction.C
//@brief  CompoFraction class
//@author aics

#include "CompoFraction.h"

// #################################################################
// 円筒領域のbboxを計算
// 標準位置で円周上の点をサンプリングし，逆変換後，min/max
REAL_TYPE CompoFraction::calcBboxCircle(Vec3r& mn, Vec3r& mx)
{
  Vec3r r, q;
  int div_r = 1000; // 周上の分割数

  REAL_TYPE pi = 2.0*asin(1.0);
  REAL_TYPE dth = 2.0*pi/(REAL_TYPE)div_r; // radian
  
  
  for (int i=0; i<div_r; i++)
  {
    REAL_TYPE d = dth * (REAL_TYPE)i;
    REAL_TYPE x = R1 * cos(d);
    REAL_TYPE y = R1 * sin(d);
    
    // 下面
    if (shape == shape_cylinder)
    {
      q = rotate_inv(angle, r.assign(x, y, 0.0)) + center;
      Geometry::get_min(mn, q);
      Geometry::get_max(mx, q);
    }

    
    // 上面
    q = rotate_inv(angle, r.assign(x, y, depth)) + center;
    Geometry::get_min(mn, q);
    Geometry::get_max(mx, q);
  }
  
  // 投影面積
  REAL_TYPE a = (REAL_TYPE)(pi*(R1*R1-R2*R2));
  
  return a;
}


// #################################################################
// 直方体領域のbboxを計算し、投影面積を計算
REAL_TYPE CompoFraction::calcBboxRect(Vec3r& mn, Vec3r& mx)
{
  Vec3r p[8];
  Vec3r o = center;
  Vec3r u = dir;
  Vec3r w = nv;
  
  Vec3r v = cross(w, u).normalize();
  
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
    Geometry::get_min(mn, p[i]);
    Geometry::get_max(mx, p[i]);
  }
  
  // 投影面積
  REAL_TYPE a = (p[0]-p[3]).length();
  REAL_TYPE b = (p[0]-p[4]).length();
  
  return a*b;
}


// #################################################################
/**
 * @brief 円柱底面と線分の交点を求める
 * @param [in] st         bbox開始インデクス
 * @param [in] ed         bbox終了インデクス
 * @param [in,out] bid    境界ID（5ビット幅x6方向）
 * @param [in,out] cut    カット情報
 * @param [in]     pl     テストする平面方程式の係数
 * @param [in]     s_id   エンコードするID
 * @param [in]     Dsize  サイズ
 */
int CompoFraction::CylinderPlane(const int st[],
                                 const int ed[],
                                 int* bid,
                                 long long* cut,
                                 const REAL_TYPE pl[4],
                                 const int s_id,
                                 const int* Dsize)
{
  int ix, jx, kx, gd;

  ix = Dsize[0];
  jx = Dsize[1];
  kx = Dsize[2];
  gd = Dsize[3];
  
  Vec3r base, o, b;
  int count = 0;
  
  o  = org;
  
  /*
   * 上下面テスト z=z1, z2の面を下面、上面とする (z1 < z2)
   *    テスト基準点b(i,j,k)の内外判定により、外部点であるとき、単位セルの6方向に対してテスト
   *    b_z<=z1 かつ p>=z1の場合、intersectLineByPlane()でテストし、交点があれば円内かどうかを判定、円内であれば交点とする
   *    b_z>=z2 かつ p<=z2の場合も同様にテストする
   *    得られた交点までの距離は無次元
   */
  
  
  for (int k=st[2]; k<=ed[2]; k++) {
    for (int j=st[1]; j<=ed[1]; j++) {
      for (int i=st[0]; i<=ed[0]; i++) {
        
        // position of cell center
        base.assign((REAL_TYPE)i-0.5, (REAL_TYPE)j-0.5, (REAL_TYPE)k-0.5);
        b = o + base * pch;
        
        
        // 基点bが円柱のz軸方向の外側にあることが条件
        if ( b.z < 0.0 || b.z > depth )
        {
          // 隣接点の座標値
          Vec3r p[6];
          p[0] = shift_W(b, pch); // (-1, 0, 0)
          p[1] = shift_E(b, pch); // ( 1, 0, 0)
          p[2] = shift_S(b, pch); // ( 0,-1, 0)
          p[3] = shift_N(b, pch); // ( 0, 1, 0)
          p[4] = shift_B(b, pch); // ( 0, 0,-1)
          p[5] = shift_T(b, pch); // ( 0, 0, 1)
          
          // テンポラリに保持
          size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
          int bb = bid[m];
          long long cc = cut[m];
          
          // 各方向の交点を評価、短い距離を記録する。新規記録の場合のみカウント
          count += updateCutPlane(b, p[0], pl, cc, bb, X_minus, s_id);
          count += updateCutPlane(b, p[1], pl, cc, bb, X_plus,  s_id);
          count += updateCutPlane(b, p[2], pl, cc, bb, Y_minus, s_id);
          count += updateCutPlane(b, p[3], pl, cc, bb, Y_plus,  s_id);
          count += updateCutPlane(b, p[4], pl, cc, bb, Z_minus, s_id);
          count += updateCutPlane(b, p[5], pl, cc, bb, Z_plus,  s_id);
          
          bid[m] = bb;
          cut[m] = cc;
        }
      }
    }
  }
  
  return count;
}


// #################################################################
/**
 * @brief 円柱側面と線分の交点を求める
 * @param [in] st         bbox開始インデクス
 * @param [in] ed         bbox終了インデクス
 * @param [in,out] bid    境界ID（5ビット幅x6方向）
 * @param [in,out] cut    カット情報
 * @param [in]     s_id   エンコードするID
 * @param [in]     Dsize  サイズ
 */
int CompoFraction::CylinderSide(const int st[],
                                const int ed[],
                                int* bid,
                                long long* cut,
                                const int s_id,
                                const int* Dsize)
{
  int ix, jx, kx, gd;
  
  ix = Dsize[0];
  jx = Dsize[1];
  kx = Dsize[2];
  gd = Dsize[3];
  
  Vec3r base, o, b;
  int count = 0;
  
  o  = org;
  
  REAL_TYPE r_2 = R1*R1;
  
  
  /* 円筒側面テスト
   *    テスト基準点b(i,j,k)の内外判定により、外部点のとき、単位セルの6方向に対してテスト
   *    itersectLineByCylinder()でテスト
   *    交点Xが円筒内かどうかをチェック
   *    無次元化
   *    必ず、円柱の外側の点を基準とする。パラメータtで線分を表すと、解が2点あるときには近い方の交点を選択する
   */
  
  for (int k=st[2]; k<=ed[2]; k++) {
    for (int j=st[1]; j<=ed[1]; j++) {
      for (int i=st[0]; i<=ed[0]; i++) {
        
        // position of cell center
        base.assign((REAL_TYPE)i-0.5, (REAL_TYPE)j-0.5, (REAL_TYPE)k-0.5);
        b = o + base * pch;
        
        
        // 基点bが円の外側にあることが条件
        if ( b.x*b.x+b.y*b.y > r_2 )
        {
          // 隣接点の座標値
          Vec3r p[6];
          p[0] = shift_W(b, pch); // (-1, 0, 0)
          p[1] = shift_E(b, pch); // ( 1, 0, 0)
          p[2] = shift_S(b, pch); // ( 0,-1, 0)
          p[3] = shift_N(b, pch); // ( 0, 1, 0)
          p[4] = shift_B(b, pch); // ( 0, 0,-1)
          p[5] = shift_T(b, pch); // ( 0, 0, 1)
          
          
          // テンポラリに保持
          size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
          int bb = bid[m];
          long long cc = cut[m];
          
          
          // 各方向の交点を評価、短い距離を記録する。新規記録の場合のみカウント
          count += updateCutCylinder(b, p[0], cc, bb, X_minus, s_id);
          count += updateCutCylinder(b, p[1], cc, bb, X_plus,  s_id);
          count += updateCutCylinder(b, p[2], cc, bb, Y_minus, s_id);
          count += updateCutCylinder(b, p[3], cc, bb, Y_plus,  s_id);
          count += updateCutCylinder(b, p[4], cc, bb, Z_minus, s_id);
          count += updateCutCylinder(b, p[5], cc, bb, Z_plus,  s_id);

          bid[m] = bb;
          cut[m] = cc;
        }
        
      }
    }
  }
  
  return count;
}


// #################################################################
// 指定法線nvがz軸の方向ベクトルに向かう回転角を計算する
// 回転の符号はz軸に向かう回転が右ねじ方向の場合を正にとる
void CompoFraction::getAngle()
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
    alpha = (nv.y >= 0.0) ? d : -d; // alpahはx軸周りの回転角
  }
  else
  {
    alpha = 0.0; // yz面への射影ベクトルがゼロの場合には回転しない
  }

  
  // 参照ベクトルをalphaだけ回転して評価ベクトルを生成 > xz面への射影
  q.assign(alpha, 0.0, 0.0);
  p = rotate(q, nv);
  c = p.length();
  
  if ( c != 0.0 )
  {
    c_bta = dot(z, p)/c;
    d = acos( c_bta );
    f_xz = c_bta+1.0;
    beta = (nv.x >= 0.0) ? -d : d; // betaはy軸まわりの回転角
  }
  else
  {
    beta = 0.0;
  }
  
  // pがz-方向の場合にだけ，y軸回りのみ回転させてz+方向にする
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
  stamped_printf("rotation angle = (%f %f %f)\n", angle.x, angle.y, angle.z);
#endif
// ##########
  
}



// #################################################################
// 形状のbboxと投影面積を求める
REAL_TYPE CompoFraction::getBboxArea(int* st, int* ed)
{
  box_min.assign(1.0e6, 1.0e6, 1.0e6);
  box_max.assign(-1.0e6, -1.0e6, -1.0e6);
  REAL_TYPE a;
  
  if ( smode == mon_BOX )
  {
    a = calcBboxRect(box_min, box_max);
  }
  else
  {
    a = calcBboxCircle(box_min, box_max);
  }

// ##########
#if 0
  stamped_printf("bbox min : %f %f %f\n", box_min.x, box_min.y, box_min.z);
  stamped_printf("bbox max : %f %f %f\n", box_max.x, box_max.y, box_max.z);
#endif
// ##########
  
  
  // コンポーネントの属するセルインデクスを求める
  findIndex(st, box_min);
  findIndex(ed, box_max);
  
  // >> 1層外側まで拡大
  st[0]--;
  st[1]--;
  st[2]--;
  ed[0]++;
  ed[1]++;
  ed[2]++;
  
  if ( st[0] < 1 ) st[0] = 1;
  if ( st[1] < 1 ) st[1] = 1;
  if ( st[2] < 1 ) st[2] = 1;
  
  if ( ed[0] > size[0] ) ed[0] = size[0];
  if ( ed[1] > size[1] ) ed[1] = size[1];
  if ( ed[2] > size[2] ) ed[2] = size[2];
  
  // ここで得られたst[],ed[]の値は、VoxInfo.CのencVIBCrev()で書き換えられる
  
  return a;
}


// #################################################################
/** 
 * @brief 円柱と線分の交点を求める
 * @param [in] st         bbox開始インデクス
 * @param [in] ed         bbox終了インデクス
 * @param [in,out] bid    境界ID（5ビット幅x6方向）
 * @param [in,out] cut    カット情報
 * @param [in]     tgt_id 固体媒質ID
 * @param [in]     Dsize  サイズ
 */
bool CompoFraction::intersectCylinder(const int st[],
                                      const int ed[],
                                      int* bid,
                                      long long* cut,
                                      const int tgt_id,
                                      const int* Dsize)
{
  int ix, jx, kx, gd;
  
  if ( !Dsize )
  {
    ix = size[0];
    jx = size[1];
    kx = size[2];
    gd = guide;
  }
  else // ASD module用
  {
    ix = Dsize[0];
    jx = Dsize[1];
    kx = Dsize[2];
    gd = 1;
  }
  
  int m_sz[4] = {ix, jx, kx, gd};
  
  
  int count = 0;
  
  // 平面の方程式の係数 ax+by+cz+d=0
  // 円筒の法線をz軸方向に回転させ、原点に平行移動した状態なので、xy平面
  REAL_TYPE pl[4];
  pl[0] = 0.0;
  pl[1] = 0.0;
  pl[2] = 1.0;
  
  
  // 上面
  pl[3] = depth;
  count += CylinderPlane(st, ed, bid, cut, pl, tgt_id, m_sz);
  
  
  if (shape == shape_cylinder)
  {
    // 下面
    pl[3] = 0.0;
    count += CylinderPlane(st, ed, bid, cut, pl, tgt_id, m_sz);
    
    // 側面
    count += CylinderSide(st, ed, bid, cut, tgt_id, m_sz);
  }
  
  return true;
}



// #################################################################
/** 
 * @brief 円と線分の交点を求める
 * @param [out]  X      交点座標
 * @param [in]   A      始点座標ベクトル
 * @param [in]   B      終点座標ベクトル
 * @retval  交点距離の小さい方を返す。負値の場合は交点が無い
 */
REAL_TYPE CompoFraction::intersectLineByCylinder(Vec3r& X,
                                                 const Vec3r A,
                                                 const Vec3r B)
{
  Vec3r v(B.x-A.x, B.y-A.y, B.z-A.z);
  REAL_TYPE a = v.x * v.x + v.y * v.y;
  REAL_TYPE b = v.x * A.x + v.y * A.y;
  REAL_TYPE c = A.x * A.x + A.y * A.y - R1 * R1;
  REAL_TYPE D = b * b - a * c;
  REAL_TYPE t = -1.0;
  
  if ( D < 0.0 )
  {
    return -1.0;
  }
  else if ( D == 0.0 )
  {
    t = - b / a;
  }
  else
  {
    // 2つの解を調べ、小さい（Aから近い方）を選択
    REAL_TYPE t1 = (-b - sqrt(D)) / a;
    REAL_TYPE t2 = (-b + sqrt(D)) / a;
    t = (t1<=t2) ? t1 : t2;
    
    // 0.0 <= t <= 1.0
    if ( t<0.0 || t>1.0 ) return -1.0;
  }
  
  // 交点を求める
  X = A + v * t;
  
  return t;
}


// #################################################################
/**
 * @brief 平面と線分の交点を求める
 * @param [out]  X   平面P上の交点Xの座標
 * @param [in]   A   線分ABの端点
 * @param [in]   B   線分ABの端点
 * @param [in]   PL  平面PLの係数 ax+by+cz+d=0
 * @retval 平面P上の交点XのAからの距離(0<=r<=1), 負値の場合は交点が無い
 * @see http://www.sousakuba.com/Programming/gs_plane_line_intersect.html
 */

REAL_TYPE CompoFraction::intersectLineByPlane(Vec3r& X, const Vec3r A, const Vec3r B, const REAL_TYPE PL[4])
{
  // 平面PLの法線
  Vec3r n(PL[0], PL[1], PL[2]);
  
  // 平面PL上の点P
  Vec3r P = n * PL[3];
  
  //printf("(%f %f %f %f) (%f %f %f) (%f %f %f)\n", PL[0], PL[1], PL[2], PL[3], A.x, A.y, A.z, B.x, B.y, B.z);
  
  // vectors
  Vec3r a(A.x-P.x, A.y-P.y, A.z-P.z);
  Vec3r b(B.x-P.x, B.y-P.y, B.z-P.z);
  
  // 平面法線と内積をとり交点があるかを符号から判別
  REAL_TYPE da = dot(a, n);
  REAL_TYPE db = dot(b, n);
  REAL_TYPE r;
  
  if ( da * db <= 0.0 )
  {
    // 比率を求める
    REAL_TYPE abs_da = fabs(da);
    REAL_TYPE abs_db = fabs(db);
    r = abs_da / ( abs_da + abs_db );
    
    // 交点を求める
    Vec3r ab(B.x-A.x, B.y-A.y, B.z-A.z);
    X = A + ab * r;
    //printf("(%f %f %f) da=%f : %f db=%f : %f r=%f\n", n.x, n.y, n.z, da, abs_da, db, abs_db, r);
  }
  else
  {
    r = -1.0;
  }
  
  return r;
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
                                   const REAL_TYPE m_R1,
                                   const REAL_TYPE m_R2,
                                   const int m_shape)
{
  smode  = mon_CYLINDER;
  nv     = m_nv;
  center = m_ctr;
  depth  = m_depth;
  R1     = m_R1;
  R2     = m_R2;
  shape  = m_shape;
  
  nv.normalize();
  
  if ( nv.length() == 0.0 )
  {
    stamped_printf("\tError : Invalid parameter of Fan : zero vector\n");
    Exit(0);
  }
  
  printf("\tSolid Revolution : %s\n\n", (shape == shape_plate)? "PLATE":"CYLINDER");
}


// #################################################################
// 体積率が(0,1)の間のセルに対してサブディビジョンを実施
void CompoFraction::subdivision(const int st[], const int ed[], REAL_TYPE* vf, double& flop)
{
  Vec3r base, b;
  Vec3r p, o, h;
  
  // for optimization > variables defined outside
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  int dv = division;
  h  = pch/(REAL_TYPE)dv;
  REAL_TYPE ff = 1.0/(REAL_TYPE)(dv*dv*dv);
  o  = org;
  
  if ( smode == mon_BOX )
  {
    for (int k=st[2]; k<=ed[2]; k++) {
      for (int j=st[1]; j<=ed[1]; j++) {
        for (int i=st[0]; i<=ed[0]; i++) {
          size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
          REAL_TYPE r = vf[m];
          
          if ( (r>0.0) && (r<1.0) ) {
            base.assign((REAL_TYPE)i-1.0, (REAL_TYPE)j-1.0, (REAL_TYPE)k-1.0);
            b = o + base * pch;
            
            int c = 0;
            for (int k1=0; k1<dv; k1++) {
              p.z = b.z + ((REAL_TYPE)k1+0.5)*h.z;
              flop += 9.0;
              
              for (int j1=0; j1<dv; j1++) {
                p.y = b.y + ((REAL_TYPE)j1+0.5)*h.y;
                flop += 9.0;
                
                for (int i1=0; i1<dv; i1++) {
                  p.x = b.x + ((REAL_TYPE)i1+0.5)*h.x;
                  
                  c += judgeRect(p);
                  flop += 10.0 + 183.0;
                }
              }
            } // k1
            vf[m] = (REAL_TYPE)c*ff;
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
          size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
          REAL_TYPE r = vf[m];
          
          if ( (r>0.0) && (r<1.0) ) {
            base.assign((REAL_TYPE)i-1.0, (REAL_TYPE)j-1.0, (REAL_TYPE)k-1.0);
            base = o + base * pch;
            
            int c = 0;
            for (int k1=0; k1<dv; k1++) {
              p.z = base.z + ((REAL_TYPE)k1+0.5)*h.z;
              flop += 9.0;
              
              for (int j1=0; j1<dv; j1++) {
                p.y = base.y + ((REAL_TYPE)j1+0.5)*h.y;
                flop += 9.0;
                
                for (int i1=0; i1<dv; i1++) {
                  p.x = base.x + ((REAL_TYPE)i1+0.5)*h.x;
                  
                  c += judgeCylider(p);
                  flop += 10.0 + 186.0;
                }
              }
            } // k1
            vf[m] = (REAL_TYPE)c*ff;
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
  
  // for optimization > variables defined outside
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  o  = org;
  
  if ( smode == mon_BOX )
  {
    // Rect cylinder
    for (int k=st[2]; k<=ed[2]; k++) {
      for (int j=st[1]; j<=ed[1]; j++) {
        for (int i=st[0]; i<=ed[0]; i++) {
          base.assign((REAL_TYPE)i-1.0, (REAL_TYPE)j-1.0, (REAL_TYPE)k-1.0);
          b    = o + base * pch;
          
          p[0] =          b;       // (0,0,0)
          p[1] = shift_f1(b, pch); // (1,0,0)
          p[2] = shift_f2(b, pch); // (0,1,0)
          p[3] = shift_f3(b, pch); // (1,1,0)
          p[4] = shift_f4(b, pch); // (0,0,1)
          p[5] = shift_f5(b, pch); // (1,0,1)
          p[6] = shift_f6(b, pch); // (0,1,1)
          p[7] = shift_f7(b, pch); // (1,1,1)
          
          int c = 0;
          
          for (int l=0; l<8; l++) {
            c += judgeRect(p[l]);
          }
          size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
          vf[m] += (REAL_TYPE)c*0.125; // 1/8
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
          b    = o + base * pch;
          
          p[0] =          b;       // (0,0,0)
          p[1] = shift_f1(b, pch); // (1,0,0)
          p[2] = shift_f2(b, pch); // (0,1,0)
          p[3] = shift_f3(b, pch); // (1,1,0)
          p[4] = shift_f4(b, pch); // (0,0,1)
          p[5] = shift_f5(b, pch); // (1,0,1)
          p[6] = shift_f6(b, pch); // (0,1,1)
          p[7] = shift_f7(b, pch); // (1,1,1)
          
          int c = 0;
          
          for (int l=0; l<8; l++) {
            c += judgeCylider(p[l]);
          }
          size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
          vf[m] += (REAL_TYPE)c*0.125; // 1/8
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
  
  // for optimization > variables defined outside
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  int odr = id;
  o  = org;
  
  if ( smode == mon_BOX )
  {
    // Rect cylinder
    for (int k=st[2]; k<=ed[2]; k++) {
      for (int j=st[1]; j<=ed[1]; j++) {
        for (int i=st[0]; i<=ed[0]; i++) {
          base.assign((REAL_TYPE)i-1.0, (REAL_TYPE)j-1.0, (REAL_TYPE)k-1.0);
          b    = o + base * pch;
          
          p[0] =          b;       // (0,0,0)
          p[1] = shift_f1(b, pch); // (1,0,0)
          p[2] = shift_f2(b, pch); // (0,1,0)
          p[3] = shift_f3(b, pch); // (1,1,0)
          p[4] = shift_f4(b, pch); // (0,0,1)
          p[5] = shift_f5(b, pch); // (1,0,1)
          p[6] = shift_f6(b, pch); // (0,1,1)
          p[7] = shift_f7(b, pch); // (1,1,1)
          
          int c = 0;
          
          for (int l=0; l<8; l++)
          {
            c += judgeRect(p[l]);
          }

          if ( c>=4 )
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
          b    = o + base * pch;
          
          p[0] =          b;       // (0,0,0)
          p[1] = shift_f1(b, pch); // (1,0,0)
          p[2] = shift_f2(b, pch); // (0,1,0)
          p[3] = shift_f3(b, pch); // (1,1,0)
          p[4] = shift_f4(b, pch); // (0,0,1)
          p[5] = shift_f5(b, pch); // (1,0,1)
          p[6] = shift_f6(b, pch); // (0,1,1)
          p[7] = shift_f7(b, pch); // (1,1,1)
          
          int c = 0;
          
          for (int l=0; l<8; l++)
          {
            c += judgeCylider(p[l]);
          }
          
          if ( c>=4 )
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
