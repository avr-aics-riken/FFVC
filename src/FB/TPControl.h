#ifndef _FB_TPCONTROL_H_
#define _FB_TPCONTROL_H_

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
 @file TPControl.h
 @brief TextParser Control class Header
 @author kero
 */

#include <math.h>

#include "FB_Define.h"

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include "string.h"
#include "TextParser.h"


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
  std::vector<SubDomain> m_subDomain;
};



class MediumTableInfo {
public:
  int type;
  std::string label;
  std::map<std::string, REAL_TYPE> m_fval;

public:
  MediumTableInfo() 
  {
    type = -1;
  }
  ~MediumTableInfo() 
  {
  }
  
};

class TPControl {

private:
	TextParser* tp;

public:
	TPControl(){};

	~TPControl(){};
	

	//変数取得関数
	bool GetVector(std::string label, int *vec, const int nvec);
	bool GetVector(std::string label, REAL_TYPE *vec, const int nvec);
	bool GetVector(std::string label, std::string *vec, const int nvec);
	bool GetValue(const std::string label, int *ct);
	bool GetValue(const std::string label, REAL_TYPE *ct);
	bool GetValue(const std::string label, std::string *ct);

	//Label、Nodeのチェック関数
	bool chkLabel(const std::string label);
	bool chkNode(const std::string label);
	bool GetNodeStr(const std::string label, const int nnode, std::string *ct);
  
  int countLabels(const std::string label);
  int getTPinstance(void);                     //TextParserをインスタンスする
	int readTPfile(const std::string filename);  //入力ファイルをTextParserへ読み込ませる

};

#endif // _FB_TPCONTROL_H_
