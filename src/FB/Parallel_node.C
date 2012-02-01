/*
 * SPHERE - Skeleton for PHysical and Engineering REsearch
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2004-2010
 *
 */

//@file Parallel_node.C
//@brief FlowBase Parallel_node class
//@author keno, FSI Team, VCAD, RIKEN

#include "Parallel_node.h"

/** 
 @fn void Parallel_Node::setParallelInfo(Parallel_Info& ref)
 @brif 並列情報をセットする
 @param ref Parallel_Infoクラスの参照
 */
void Parallel_Node::setParallelInfo(Parallel_Info& ref)
{
  pn.procGrp   = ref.procGrp;
  pn.ID        = ref.ID;
  pn.st_idx[0] = ref.st_idx[0];
  pn.st_idx[1] = ref.st_idx[1];
  pn.st_idx[2] = ref.st_idx[2];
  
  for(int i=0; i<6; i++) pn.nID[i] = ref.nID[i];
}
