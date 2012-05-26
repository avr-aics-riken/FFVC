/*
 * SPHERE - Skeleton for PHysical and Engineering REsearch
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2004-2012
 *
 */

//@file Monitor.C
//@brief FlowBase MonitorList class
//@author keno, FSI Team, VCAD, RIKEN

#include "Monitor.h"

/// Line指定の端点座標をグローバル計算領域内にクリッピング.
///
///   @param[in,out] from Line始点
///   @param[in,out] to   Line終点
///
void MonitorList::clipLine(REAL_TYPE from[3], REAL_TYPE to[3])
{
  const char* OUT_OF_REGION = "out of region";
  Vec3r st(from);
  Vec3r ed(to);
  Vec3r r0 = g_org;
  Vec3r r1 = r0 + g_box;
  Vec3r d = ed - st;
  
  
  // 境界上に端点が存在する場合は、少し計算領域内にずらす
  REAL_TYPE eps = 1.0e-5;
  r0 += eps * pch;
  r1 -= eps * pch;
  
  REAL_TYPE t_st = 0.0;
  REAL_TYPE t_ed = 1.0;
  
  try {
    if (d.x == 0) {
      if (st.x < r0.x || r1.x < st.x) throw OUT_OF_REGION;
    }
    else if (d.x > 0) {
      if (ed.x < r0.x || r1.x < st.x) throw OUT_OF_REGION;
      t_st = max(t_st, (r0.x-st.x)/d.x);
      t_ed = min(t_ed, (r1.x-st.x)/d.x);
    }
    else if (d.x < 0) {
      if (st.x < r0.x || r1.x < ed.x) throw OUT_OF_REGION;
      t_st = max(t_st, (r1.x-st.x)/d.x);
      t_ed = min(t_ed, (r0.x-st.x)/d.x);
    }
    
    if (d.y == 0) {
      if (st.y < r0.y || r1.y < st.y) throw OUT_OF_REGION;
    }
    else if (d.y > 0) {
      if (ed.y < r0.y || r1.y < st.y) throw OUT_OF_REGION;
      t_st = max(t_st, (r0.y-st.y)/d.y);
      t_ed = min(t_ed, (r1.y-st.y)/d.y);
    }
    else if (d.y < 0) {
      if (st.y < r0.y || r1.y < ed.y) throw OUT_OF_REGION;
      t_st = max(t_st, (r1.y-st.y)/d.y);
      t_ed = min(t_ed, (r0.y-st.y)/d.y);
    }
    
    if (d.z == 0) {
      if (st.z < r0.z || r1.z < st.z) throw OUT_OF_REGION;
    }
    else if (d.z > 0) {
      if (ed.z < r0.z || r1.z < st.z) throw OUT_OF_REGION;
      t_st = max(t_st, (r0.z-st.z)/d.z);
      t_ed = min(t_ed, (r1.z-st.z)/d.z);
    }
    else if (d.z < 0) {
      if (st.z < r0.z || r1.z < ed.z) throw OUT_OF_REGION;
      t_st = max(t_st, (r1.z-st.z)/d.z);
      t_ed = min(t_ed, (r0.z-st.z)/d.z);
    }
    
    if (t_st >= t_ed) throw OUT_OF_REGION;
  }
  catch (const char *str) {
    if (pn.myrank == 0) {
      printf("Line [%14.6e %14.6e %14.6e]-[%14.6e %14.6e %14.6e] is %s\n", // %12.4 >> %14.6
             from[0], from[1], from[2], to[0], to[1], to[2], str);
    }
    Exit(0);
  }
  
  ed = st + t_ed * d;
  st = st + t_st * d;
  
  from[0] = st.x;
  from[1] = st.y;
  from[2] = st.z;
  to[0]   = ed.x;
  to[1]   = ed.y;
  to[2]   = ed.z;
}


/// 出力ファイルクローズ.
///
///    @note 出力タイプによらず全プロセスから呼んでも問題ない
///
void MonitorList::closeFile()
{
  if ((outputType == GATHER && pn.myrank == 0) ||
      outputType == DISTRIBUTE) {
    for (int i = 0; i < nGroup; i++) {
      if (monGroup[i]->getType() != MonitorCompo::INNER_BOUNDARY) monGroup[i]->closeFile();
    }
  }
}


/// 出力タイプ文字列の取得.
string MonitorList::getOutputTypeStr()
{
  string str;
  
  switch (outputType) {
    case GATHER:
      str = "Gather";
      break;
      
    case DISTRIBUTE:
      str = "Distribute";
      break;
      
    case NONE:
      str = "Non";
      break;
      
    default:
      Exit(0);
  }
  
  return str;
}


/// 出力ファイルオープン.
///
///    @param str ファイル名テンプレート
///
///    @note 出力タイプによらず全プロセスから呼んでも問題ない
///
void MonitorList::openFile(const char* str)
{
  if ((outputType == GATHER && pn.myrank == 0) ||
      outputType == DISTRIBUTE) {
    for (int i = 0; i < nGroup; i++) {
      if (monGroup[i]->getType() != MonitorCompo::INNER_BOUNDARY) {
        string fileName(str);
        string label("_");
        label = label + monGroup[i]->getLabel();
        string::size_type pos = fileName.rfind(".");
        fileName.insert(pos, label);
        if (outputType == GATHER) {
          monGroup[i]->openFile(fileName.c_str(), true);
        }
        else {
          monGroup[i]->openFile(fileName.c_str(), false);
        }
      }
    }
  }
}

/// モニタ結果出力(Line, PointSet指定).
///
///   @param[in] step サンプリング時の計算ステップ
///   @param[in] tm   サンプリング時の計算時刻
///
///   @note gatherの場合も全プロセスから呼ぶこと
///
void MonitorList::print(unsigned step, REAL_TYPE tm)
{
  for (int i = 0; i < nGroup; i++) {
    if (monGroup[i]->getType() != MonitorCompo::INNER_BOUNDARY) {
      switch (outputType) {
        case GATHER:
          monGroup[i]->print_gather(step, tm);
          break;
          
        case DISTRIBUTE:
          monGroup[i]->print_distribute(step, tm);
          break;
          
        default:
          Exit(0);
      }
    }
  }
  
}


/// モニタ情報を出力.
void MonitorList::printMonitorInfo(FILE* fp, const char* str, const bool verbose)
{
  fprintf(fp,"\n---------------------------------------------------------------------------\n\n");
  fprintf(fp, "\n\t>> Monitor Information\n\n");
  
  fprintf(fp, "\t  Output Type              : %s\n", getOutputTypeStr().c_str());
  fprintf(fp, "\t  Base Name of Output File : %s\n", str);
  fprintf(fp, "\t  Specified Unit           : %s\n\n", (refVar.modeUnit == DIMENSIONAL) ? "Dimensional" : "Non Dimensional");
  
  if ( verbose ) {
    for (int i = 0; i < nGroup; i++) {
      monGroup[i]->printInfo(fp, i);
    }
  }
  
  fprintf(fp,"\n");
  fflush(fp);
}




/// 必要なパラメータのコピー.
///
///   @param[in] bcd BCindex ID
///   @param[in] g_org,g_lbx  グローバル領域基点座標，領域サイズ
///   @param[in] org,dx,lbx   ローカル領域基点座標，セル幅，領域サイズ
///   @param[in] rs, gc       ローカルセルサイズ，ガイドセル数
///   @param[in] refVelocity    代表速度
///   @param[in] baseTemp       基準温度
///   @param[in] diffTemp       代表温度差
///   @param[in] refDensity     基準密度
///   @param[in] refLength      代表長さ
///   @param[in] basePrs        基準圧力
///   @param[in] unitTemp       温度単位指定フラグ (Kelvin / Celsius)
///   @param[in] modePrecision  出力精度指定フラグ (単精度，倍精度)
///   @param[in] unitPrs        圧力単位指定フラグ (絶対値，ゲージ圧)
///
void MonitorList::setControlVars(unsigned* bcd,
                                 REAL_TYPE g_org[3], REAL_TYPE g_lbx[3],
                                 REAL_TYPE org[3], REAL_TYPE dx[3], REAL_TYPE lbx[3],
                                 unsigned rs[3], unsigned gc,
                                 REAL_TYPE refVelocity, REAL_TYPE baseTemp, REAL_TYPE diffTemp, 
                                 REAL_TYPE refDensity, REAL_TYPE refLength, REAL_TYPE basePrs, 
                                 unsigned unitTemp, unsigned modePrecision, unsigned unitPrs) 
{
  this->bcd = bcd;
  this->g_org = g_org;
  this->g_box = g_lbx;
  this->org = org;
  this->pch = dx;
  this->box = lbx;
  
  size[0] = rs[0];
  size[1] = rs[1];
  size[2] = rs[2];
  guide = gc;
  
  refVar.refVelocity = refVelocity;
  refVar.baseTemp    = baseTemp;
  refVar.diffTemp    = diffTemp;
  refVar.refDensity  = refDensity;
  refVar.refLength   = refLength;
  refVar.basePrs     = basePrs;
  refVar.unitTemp    = unitTemp;
  refVar.unitPrs     = unitPrs;
  refVar.modePrecision = modePrecision;
}

/// 参照速度のコピー
///   @param[in] v00 参照（座標系移動）速度
///
void MonitorList::set_V00(REAL_TYPE v00[4]) 
{
  refVar.v00.x = v00[1];
  refVar.v00.y = v00[2];
  refVar.v00.z = v00[3];
}



/// PointSet登録.
///
///   @param[in] str ラベル文字列
///   @param[in] variables モニタ変数vector
///   @param[in] method method文字列
///   @param[in] mode   mode文字列
///   @param[in] pointSet  PointSet
///
void MonitorList::setPointSet(const char* str, vector<string>& variables,
                              const char* method, const char* mode,
                              vector<MonitorCompo::MonitorPoint>& pointSet)
{
  MonitorCompo* m = new MonitorCompo(pn, org, pch, box, g_org, g_box,
                                     size, guide, refVar, bcd);
  m->setPointSet(str, variables, method, mode, pointSet);
  
  monGroup.push_back(m);
  nGroup = monGroup.size();
}


/// Line点登録.
///
///   @param[in] str ラベル文字列
///   @param[in] variables モニタ変数vector
///   @param[in] method method文字列
///   @param[in] mode   mode文字列
///   @param[in] from Line始点
///   @param[in] to   Line終点
///   @param[in] nDivision 分割数(モニタ点数-1)
///
void MonitorList::setLine(const char* str, vector<string>& variables,
                          const char* method, const char* mode,
                          REAL_TYPE from[3], REAL_TYPE to[3], int nDivision)
{
  clipLine(from, to);
  
  MonitorCompo* m = new MonitorCompo(pn, org, pch, box, g_org, g_box,
                                     size, guide, refVar, bcd);
  m->setLine(str, variables, method, mode, from, to, nDivision);
  
  monGroup.push_back(m);
  nGroup = monGroup.size();
}


/// 内部境界条件としてモニタ点を登録.
///
///   @param[in] cmp コンポーネント配列
///   @param[in] nBC コンポーネント数
///
void MonitorList::setInnerBoundary(CompoList *cmp, int nBC)
{
  for (int i = 1; i <= nBC; i++) {
    if (cmp[i].isMONITOR()) {
      MonitorCompo* m = new MonitorCompo(pn, org, pch, box, g_org, g_box, size, guide, refVar, bcd);
      m->setInnerBoundary(i, cmp[i]);
      
      monGroup.push_back(m);
    }
  }
  nGroup = monGroup.size();
}



/// XMLにより指定されるモニタ点位置にID=255を書き込む.
///
///   @param[in] id セルID配列
///
void MonitorList::write_ID(int* id)
{
  for (int i = 0; i < nGroup; i++) {
    for (int m = 0; m < monGroup[i]->getSize(); m++) {
      if (monGroup[i]->getType() != MonitorCompo::INNER_BOUNDARY) {
        //<<<
        if (!monGroup[i]->check_region(m, org, box)) continue;
        //>>>
        Vec3i index = monGroup[i]->getSamplingCellIndex(m);
        id[FBUtility::getFindexS3D(size, guide, index.x, index.y, index.z)] = 255;
      }
    }
  }
}
