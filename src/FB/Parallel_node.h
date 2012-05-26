#ifndef _FB_PARA_NODE_H_
#define _FB_PARA_NODE_H_

/*
 * SPHERE - Skeleton for PHysical and Engineering REsearch
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2004-2012
 *
 */

//@file Parallel_node.h
//@brief FlowBase Parallel_Node class Header
//@author keno, FSI Team, VCAD, RIKEN

typedef struct {
  int procGrp;    // プロセスグループ番号
  int myrank;     // 自ノードのランク番号
  int numProc;    // 全ランク数
  int nID[6];     // 隣接ブロックのランク番号
  int st_idx[3];  // 開始インデクス（グローバルインデクス）
} Parallel_Info;

class Parallel_Node {
public:
  Parallel_Info pn;
  
public:
  Parallel_Node() {
    pn.procGrp = 0;
    pn.myrank  = -1;
    pn.numProc = 0;
    for (int i=0; i<6; i++) pn.nID[i] = -1;
    pn.st_idx[0] = pn.st_idx[1] = pn.st_idx[2] = 0;
  }
  virtual ~Parallel_Node() {}
  
public:

  /// @brif ランクとグループをセットする
  void set_Rank(Parallel_Info& ref)
  {
    pn.procGrp   = ref.procGrp;
    pn.myrank    = ref.myrank;
    pn.numProc   = ref.numProc;
  }
  
  /// @brif ランクとグループをセットする
  void set_Neighbor(Parallel_Info& ref)
  {
    pn.st_idx[0] = ref.st_idx[0];
    pn.st_idx[1] = ref.st_idx[1];
    pn.st_idx[2] = ref.st_idx[2];
    
    for(int i=0; i<6; i++) pn.nID[i] = ref.nID[i];
  }
  
};

#endif // _FB_PARA_NODE_H_
