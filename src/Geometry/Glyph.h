#ifndef _FB_GLYPH_H_
#define _FB_GLYPH_H_

//##################################################################################
//
// Glyph class
//
// Copyright (c) 2007-2011 VCAD System Research Program, RIKEN.
// All rights reserved.
//
// Copyright (c) 2011-2015 Institute of Industrial Science, The University of Tokyo.
// All rights reserved.
//
// Copyright (c) 2012-2016 Advanced Institute for Computational Science, RIKEN.
// All rights reserved.
//
//##################################################################################

//@file   Glyph.h
//@brief  Glyph class Header
//@author aics

#ifndef DISABLE_MPI
  #include "mpi.h" // add header explicitly to avoid compile error for Intel MPI
#else
  #include "cpm_mpistub.h"
#endif

#include <iostream>
#include <fstream>
#include <string>
#include "string.h"
#include <stdlib.h>
#include <iomanip>
#include "FB_Define.h"
#include "common/Vec3.h" // defined in Polylib

using namespace std;
using namespace Vec3class;

class Glyph {

private:
  Vec3f pch;         ///< セル幅
  Vec3f org;         ///< 計算領域の基点
  unsigned element;  ///< ポリゴン数
  unsigned poly;     ///< ポリゴンの番号
  int myRank;        ///< 自ノードのランク番号

  Vec3f* nvc; // array for normal vector of element
  Vec3f* xyz; // array for vertex coordiante
  int* b_id;


public:

  // default constructor
	Glyph() {
    element = 0;
    myRank = 0;
    poly = 0;
    nvc = xyz = NULL;
    b_id = NULL;
  }

  /** constructor
   * @param [in] m_pch  セル幅
   * @param [in] m_org  ローカル領域基点座標
   * @param [in] m_elm  交点の数
   * @param [in] m_rank ランク番号
   */
  Glyph(const float* m_pch, const float* m_org, const unsigned m_elem, const int m_rank)
  {
    this->pch.x    = m_pch[0];
    this->pch.y    = m_pch[1];
    this->pch.z    = m_pch[2];
    this->org.x    = m_org[0];
    this->org.y    = m_org[1];
    this->org.z    = m_org[2];
    this->element  = m_elem * 12; // 6-face x 2Polygons = 12
    this->myRank   = m_rank;

    poly = 0;

    if (0 != element)
    {
      if ( !(nvc=new Vec3f[element]) )
      {
        cout << "\tMemory allocation error!" << endl;
        Exit(-1);
      }

      if ( !(xyz=new Vec3f[element*3]) )
      {
        cout << "\tMemory allocation error!" << endl;
        Exit(-1);
      }

      if ( !(b_id=new int[element]) )
      {
        cout << "\tMemory allocation error!" << endl;
        Exit(-1);
      }
    }

  }

  // destructor
	~Glyph() {
    if (nvc) delete [] nvc;
    if (xyz) delete [] xyz;
    if (b_id) delete [] b_id;
  }


  /**
   * @brief グリフ作成のための頂点を生成
   * @param [in] idx    セルインデクス
   * @param [in] pos    カット距離
   * @param [in] dir    方向
   * @param [in] m_bid  (i,j,k)の境界ID
   */
  void generateVertex(const Vec3i idx, const long long pos, const int dir, const int m_bid);


  /**
   * @brief アスキー出力
   * @param [in] ofs 出力ストリーム
   */
  void writeAscii (ofstream &ofs);


  /**
   * @brief バイナリー出力
   * @param [in] outFile 出力ファイル名
   */
  void writeBinary(const string outFile);




private:

  /**
   * @brief グリフポリゴンの登録
   * @param [in] p      頂点座標
   * @param [in] m_bid  境界ID
   */
  void registerPolygon(const Vec3f p[8], const int m_bid);


  /**
   * @brief ファイル書き出し
   * @param [in] os    ostream
   * @param [in] _data データポインタ
   * @param [in] size  データ型の大きさ
   * @param [in] n     データ長
   */
  inline void tt_write(ostream& os, const void* _data, int size, int n)
  {
    const char* data = (const char*)_data;

    os.write(data, size * n);
  }

};


#endif // _FB_GLYPH_H_
