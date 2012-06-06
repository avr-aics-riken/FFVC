#ifndef _FB_TPCONTROL_H_
#define _FB_TPCONTROL_H_

/*
 * SPHERE - Skeleton for PHysical and Engineering REsearch
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2012
 *
 */

/**
 @file TPControl.h
 @brief TextParser Control class Header
 @author keno, FSI Team, VCAD, RIKEN
 */

#include <math.h>

#include "FB_Define.h"

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include "string.h"
#include "TextParser.h"

using namespace std;


class SubDomain {
public:
  SubDomain()
  {
    m_pos[0] = m_pos[1] = m_pos[2] = 0;
    m_bcid[0] = m_bcid[1] = m_bcid[2] = m_bcid[3] = m_bcid[4] = m_bcid[5] = 0;
  }
  ~SubDomain()
  {
  }

  int m_pos[3];
  int m_bcid[6];
};

class DomainInfo {
public:
  DomainInfo()
  {
    m_globalOrigin[0] = m_globalOrigin[1] = m_globalOrigin[2] = REAL_TYPE(0.0);
    m_globalRegion[0] = m_globalRegion[1] = m_globalRegion[2] = REAL_TYPE(0.0);
    m_globalPitch[0]  = m_globalPitch[1]  = m_globalPitch[2]  = REAL_TYPE(0.0);
    m_domainDiv[0]    = m_domainDiv[1]    = m_domainDiv[2]    = 0;
  }
  ~DomainInfo()
  {
  }

  REAL_TYPE m_VoxelOrigin[3];
  REAL_TYPE m_VoxelSize[3];
  REAL_TYPE m_VoxelPitch[3];

  REAL_TYPE m_globalOrigin[3];
  REAL_TYPE m_globalRegion[3];
  REAL_TYPE m_globalPitch[3];
  int       m_domainDiv[3];
  vector<SubDomain> m_subDomain;
};


class MediumTableInfo {
private:
  //IDtable*          iTable;
  
public:
  int type;
  string label;
  map<string, REAL_TYPE> m_fval;

public:
  MediumTableInfo() 
  {
  }
  ~MediumTableInfo() 
  {
  }
};

class TPControl {

protected:

private:

	TextParser* tp;

public:

	/// コンストラクタ。
	TPControl(){};

	/// デストラクタ。
	~TPControl(){};

	//
	int getTPinstance();                      //TextParserをインスタンスする
	int readTPfile( const string filename );  //入力ファイルをTextParserへ読み込ませる

	//変数取得関数
	bool GetVector(string label, int *vec, const int nvec);
	bool GetVector(string label, REAL_TYPE *vec, const int nvec);
	bool GetVector(string label, string *vec, const int nvec);
	bool GetValue(const string label, int *ct);
	bool GetValue(const string label, REAL_TYPE *ct);
	bool GetValue(const string label, string *ct);

	//Label、Nodeのチェック関数
	bool chkLabel(const string label);
	bool chkNode(const string label);
	int countLabels(const string label);
	bool GetNodeStr(const string label, const int nnode, string *ct);

};

#endif // _FB_TPCONTROL_H_
