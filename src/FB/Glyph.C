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

//@file   Glyph.C
//@brief  FlowBase Glyph class Header
//@author kero

#include "Glyph.h"

#define STL_HEAD 80	// header size for STL binary

// #################################################################
// グリフ作成のための頂点を生成
void Glyph::generateVertex(const FB::Vec3i idx, const float* pos, const char* str, const int m_bid)
{
  FB::Vec3f b;     // セルセンターのシフトインデクス
  FB::Vec3f c;     // セルセンター座標
  FB::Vec3f p[8];  // グリフの8頂点
  
  float w = 0.25;   // Glyph幅の係数
  float d = 0.05;   // Glyph厚さの係数
  
  b.assign((float)idx.x-0.5, (float)idx.y-0.5, (float)idx.z-0.5);
  c = org + b * pch;
  
  if ( !strcasecmp("x-", str))
  {
    b.assign( -pos[0], -w, -w);  p[0] = c + b * pch;
    b.assign(d-pos[0], -w, -w);  p[1] = c + b * pch;
    b.assign(d-pos[0],  w, -w);  p[2] = c + b * pch;
    b.assign( -pos[0],  w, -w);  p[3] = c + b * pch;
    b.assign( -pos[0], -w,  w);  p[4] = c + b * pch;
    b.assign(d-pos[0], -w,  w);  p[5] = c + b * pch;
    b.assign(d-pos[0],  w,  w);  p[6] = c + b * pch;
    b.assign( -pos[0],  w,  w);  p[7] = c + b * pch;
    registerPolygon(p, m_bid);
  }
  else if ( !strcasecmp("x+", str))
  {
    b.assign(pos[1]-d, -w, -w);  p[0] = c + b * pch;
    b.assign(pos[1]  , -w, -w);  p[1] = c + b * pch;
    b.assign(pos[1]  ,  w, -w);  p[2] = c + b * pch;
    b.assign(pos[1]-d,  w, -w);  p[3] = c + b * pch;
    b.assign(pos[1]-d, -w,  w);  p[4] = c + b * pch;
    b.assign(pos[1]  , -w,  w);  p[5] = c + b * pch;
    b.assign(pos[1]  ,  w,  w);  p[6] = c + b * pch;
    b.assign(pos[1]-d,  w,  w);  p[7] = c + b * pch;
    registerPolygon(p, m_bid);
  }
  else if ( !strcasecmp("y-", str))
  {
    b.assign(-w,  -pos[2], -w);  p[0] = c + b * pch;
    b.assign( w,  -pos[2], -w);  p[1] = c + b * pch;
    b.assign( w, d-pos[2], -w);  p[2] = c + b * pch;
    b.assign(-w, d-pos[2], -w);  p[3] = c + b * pch;
    b.assign(-w,  -pos[2],  w);  p[4] = c + b * pch;
    b.assign( w,  -pos[2],  w);  p[5] = c + b * pch;
    b.assign( w, d-pos[2],  w);  p[6] = c + b * pch;
    b.assign(-w, d-pos[2],  w);  p[7] = c + b * pch;
    registerPolygon(p, m_bid);
  }
  else if ( !strcasecmp("y+", str))
  {
    b.assign(-w, pos[3]-d, -w);  p[0] = c + b * pch;
    b.assign( w, pos[3]-d, -w);  p[1] = c + b * pch;
    b.assign( w, pos[3]  , -w);  p[2] = c + b * pch;
    b.assign(-w, pos[3]  , -w);  p[3] = c + b * pch;
    b.assign(-w, pos[3]-d,  w);  p[4] = c + b * pch;
    b.assign( w, pos[3]-d,  w);  p[5] = c + b * pch;
    b.assign( w, pos[3]  ,  w);  p[6] = c + b * pch;
    b.assign(-w, pos[3]  ,  w);  p[7] = c + b * pch;
    registerPolygon(p, m_bid);
  }
  else if ( !strcasecmp("z-", str))
  {
    b.assign(-w, -w,  -pos[4]);  p[0] = c + b * pch;
    b.assign( w, -w,  -pos[4]);  p[1] = c + b * pch;
    b.assign( w,  w,  -pos[4]);  p[2] = c + b * pch;
    b.assign(-w,  w,  -pos[4]);  p[3] = c + b * pch;
    b.assign(-w, -w, d-pos[4]);  p[4] = c + b * pch;
    b.assign( w, -w, d-pos[4]);  p[5] = c + b * pch;
    b.assign( w,  w, d-pos[4]);  p[6] = c + b * pch;
    b.assign(-w,  w, d-pos[4]);  p[7] = c + b * pch;
    registerPolygon(p, m_bid);
  }
  else if ( !strcasecmp("z+", str))
  {
    b.assign(-w, -w, pos[5]-d);  p[0] = c + b * pch;
    b.assign( w, -w, pos[5]-d);  p[1] = c + b * pch;
    b.assign( w,  w, pos[5]-d);  p[2] = c + b * pch;
    b.assign(-w,  w, pos[5]-d);  p[3] = c + b * pch;
    b.assign(-w, -w, pos[5]  );  p[4] = c + b * pch;
    b.assign( w, -w, pos[5]  );  p[5] = c + b * pch;
    b.assign( w,  w, pos[5]  );  p[6] = c + b * pch;
    b.assign(-w,  w, pos[5]  );  p[7] = c + b * pch;
    registerPolygon(p, m_bid);
  }
  else
  {
    Exit(0);
  }
  
}


// #################################################################
// ポリゴンを登録する．各方向２ポリゴン
void Glyph::registerPolygon(const FB::Vec3f p[8], const int m_bid)
{
  FB::Vec3f b;
  
  size_t m1;
  size_t m2;

  // X-
  b.assign(-1.0, 0.0, 0.0);
  
  m1 = poly+0;
  m2 = m1*3;
  b_id[m1]  = m_bid;
  nvc[m1]   = b;
  xyz[m2+0] = p[0];
  xyz[m2+1] = p[7];
  xyz[m2+2] = p[3];

  m1 = poly+1;
  m2 = m1*3;
  b_id[m1]  = m_bid;
  nvc[m1]   = b;
  xyz[m2+0] = p[0];
  xyz[m2+1] = p[4];
  xyz[m2+2] = p[7];
  
  
  // X+
  b.assign(1.0, 0.0, 0.0);
  
  m1 = poly+2;
  m2 = m1*3;
  b_id[m1]  = m_bid;
  nvc[m1]   = b;
  xyz[m2+0] = p[1];
  xyz[m2+1] = p[2];
  xyz[m2+2] = p[6];
  
  m1 = poly+3;
  m2 = m1*3;
  b_id[m1]  = m_bid;
  nvc[m1]   = b;
  xyz[m2+0] = p[1];
  xyz[m2+1] = p[6];
  xyz[m2+2] = p[5];
  
  
  // Y-
  b.assign(0.0, -1.0, 0.0);
  
  m1 = poly+4;
  m2 = m1*3;
  b_id[m1]  = m_bid;
  nvc[m1]   = b;
  xyz[m2+0] = p[0];
  xyz[m2+1] = p[1];
  xyz[m2+2] = p[5];
  
  m1 = poly+5;
  m2 = m1*3;
  b_id[m1]  = m_bid;
  nvc[m1]   = b;
  xyz[m2+0] = p[0];
  xyz[m2+1] = p[5];
  xyz[m2+2] = p[4];
  
  
  // Y+
  b.assign(0.0, 1.0, 0.0);
  
  m1 = poly+6;
  m2 = m1*3;
  b_id[m1]  = m_bid;
  nvc[m1]   = b;
  xyz[m2+0] = p[3];
  xyz[m2+1] = p[7];
  xyz[m2+2] = p[6];
  
  m1 = poly+7;
  m2 = m1*3;
  b_id[m1]  = m_bid;
  nvc[m1]   = b;
  xyz[m2+0] = p[3];
  xyz[m2+1] = p[6];
  xyz[m2+2] = p[2];
  
  
  // Z-
  b.assign(0.0, 0.0, -1.0);
  
  m1 = poly+8;
  m2 = m1*3;
  b_id[m1]  = m_bid;
  nvc[m1]   = b;
  xyz[m2+0] = p[0];
  xyz[m2+1] = p[3];
  xyz[m2+2] = p[1];
  
  m1 = poly+9;
  m2 = m1*3;
  b_id[m1]  = m_bid;
  nvc[m1]   = b;
  xyz[m2+0] = p[3];
  xyz[m2+1] = p[2];
  xyz[m2+2] = p[1];
  
  
  // Z+
  b.assign(0.0, 0.0, 1.0);
  
  m1 = poly+10;
  m2 = m1*3;
  b_id[m1]  = m_bid;
  nvc[m1]   = b;
  xyz[m2+0] = p[4];
  xyz[m2+1] = p[5];
  xyz[m2+2] = p[7];
  
  m1 = poly+11;
  m2 = m1*3;
  b_id[m1]  = m_bid;
  nvc[m1]   = b;
  xyz[m2+0] = p[5];
  xyz[m2+1] = p[6];
  xyz[m2+2] = p[7];
  
  poly += 12;
}



// #################################################################
// アスキー出力
void Glyph::writeAscii(ofstream &ofs)
{
	
	for (unsigned m=0; m<element; m++)
  {
    for (int k=0; k<12; k++)
    {
      ofs << "facet normal";
      
      for (int i=0; i<3; i++)
      {
        ofs << " " << setprecision(6) << scientific << *nvc;
        nvc++;
      }
      ofs << endl;
      
      ofs << "outer loop" << endl;
      
      for (int j=0; j<3; j++)
      {
        ofs << "vertex";
        
        for (int i=0; i<3; i++)
        {
          ofs << " " << setprecision(6) << scientific << *xyz;
          xyz++;
        }
        ofs << fixed << endl;
      }
      ofs << "endloop" << endl;
      ofs << "endfacet" << endl;
    }
	}
  
}



// #################################################################
// バイナリー出力
void Glyph::writeBinary(const string outFile)
{
  if (0 == element) return;
  if ( outFile.empty() ) return;
  
  unsigned long len = outFile.size() + 14; // id(9) + postfix(4) + 1
  char* tmp = new char[len];
  memset(tmp, 0, sizeof(char)*len);
  
  sprintf(tmp, "%s_id%06d.%s", outFile.c_str(), myRank, "stl");

  
  ofstream ofs(tmp, ios::out | ios::binary);
  
	if (!ofs)
  {
		cout << "\tFile '" << outFile.c_str() << "' could not open." << endl;
		Exit(0);
  }
  
  
  unsigned short c=0;
  char head[STL_HEAD];
  char *pch=NULL;
  
  pch = &head[0];
  for (int i=0; i<STL_HEAD; i++) *pch++ = 0x00;
  
  strcpy(head, "CutGlyph");

  tt_write(ofs, head, sizeof(char), STL_HEAD);
	tt_write(ofs, &element, sizeof(unsigned), 1);
  
  float v[3];
  
	for (unsigned m=0; m<element; m++)
  {
    int q = b_id[m];
    c = 0;
    c |= (0x1 << 15); // extend color format
    c |= (q << 10);   // R
    c |= (q << 5);    // G
    c |= (q << 0);    // B
    
    v[0] = nvc[m].x;
    v[1] = nvc[m].y;
    v[2] = nvc[m].z;
    tt_write(ofs, v, sizeof(float), 3);
    
    v[0] = xyz[3*m  ].x;
    v[1] = xyz[3*m  ].y;
    v[2] = xyz[3*m  ].z;
    tt_write(ofs, v, sizeof(float), 3);
    
    v[0] = xyz[3*m+1].x;
    v[1] = xyz[3*m+1].y;
    v[2] = xyz[3*m+1].z;
    tt_write(ofs, v, sizeof(float), 3);
    
    v[0] = xyz[3*m+2].x;
    v[1] = xyz[3*m+2].y;
    v[2] = xyz[3*m+2].z;
    tt_write(ofs, v, sizeof(float), 3);
    
    tt_write(ofs, &c, sizeof(unsigned short), 1);
	}
  ofs.close();
}
