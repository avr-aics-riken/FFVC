#ifndef _SKL_IP_USERS_H_
#define _SKL_IP_USERS_H_

/*
 * SPHERE - Skeleton for PHysical and Engineering REsearch
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2004-2012
 *
 */

//@file IP_Users.h
//@brief FlowBase IP_Users class Header
//@author keno, FSI Team, VCAD, RIKEN

#include "Intrinsic.h"

class IP_Users : public Intrinsic {
public:
  IP_Users(){}
  ~IP_Users() {}
  
protected:
  
public:
  // @fn const char* getExampleName(void)
  // @brief ユーザー例題の名称を返す
  const char* getExampleName(void) {
    return ("User's problem");
  }
};

#endif // _SKL_IP_USERS_H_
