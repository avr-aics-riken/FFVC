#ifndef _SKL_FB_PARA_NODE_H_
#define _SKL_FB_PARA_NODE_H_

/*
 * SPHERE - Skeleton for PHysical and Engineering REsearch
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2004-2010
 *
 */

//@file Parallel_node.h
//@brief FlowBase Parallel_Node class Header
//@author keno, FSI Team, VCAD, RIKEN

#include "Skl.h"
#include "SklSolverBase.h"

typedef struct {
  int procGrp;    // プロセスグループ番号
  int ID;         // 自ノードのランク番号
  int nID[6];     // 隣接ブロックのランク番号
  int st_idx[3];  // 開始インデクス（グローバルインデクス）
  int solverID;   // ソルバーのID
} Parallel_Info;

class Parallel_Node {
public:
  Parallel_Info pn;
  
public:
  Parallel_Node() {
    pn.procGrp = 0;
    pn.ID = 0; 
    for (int i=0; i<6; i++) pn.nID[i] = -1;
    pn.st_idx[0] = pn.st_idx[1] = pn.st_idx[2] = 0;
  }
  virtual ~Parallel_Node() {}
  
public:
  void setParallelInfo(Parallel_Info& ref);
  
};

#endif // _SKL_FB_PARA_NODE_H_