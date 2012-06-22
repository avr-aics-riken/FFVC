#ifndef _FFV_H_
#define _FFV_H_

// #################################################################
//
// FFV : Frontflow / violet
//
// Copyright (c) All right reserved. 2012
//
// Institute of Industrial Science, The University of Tokyo, Japan. 
//
// #################################################################

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <iostream>
#include <fstream>

#include "cpm_ParaManager.h"
#include "cpm_TextParserDomain.h"

#include "../FB/FB_Define.h"
#include "ffv_Define.h"

using namespace std;

class FFV {
protected:
  FILE *mp;      /// 標準出力
  
  int session_maxStep;     /// セッションのステップ数
  int session_currentStep; /// セッションの現在のステップ
  
public:
  FFV();
  virtual ~FFV();
  
  
protected:
  virtual bool stepPost(void);
  virtual bool Post(cpm_ParaManager *paraMngr);
  
  virtual int Initialize(void);
  virtual int Loop(int m_step);
  virtual int MainLoop(void);
  
  virtual void Usage(void);
  
  
  //@brief マスターノードのみ trueを返す
  bool IsMaster(cpm_ParaManager *paraMngr) {
    return ( !paraMngr->IsParallel() || (paraMngr->GetMyRankID() == 0) ) ? true : false;
  }
  
};

#endif // _FFV_H_