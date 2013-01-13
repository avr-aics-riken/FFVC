// #################################################################
//
// CAERU Library
//
// Copyright (c) 2012-2013  All right reserved.
//
// Institute of Industrial Science, University of Tokyo, Japan. 
//
// #################################################################

/**
 * @file   Monitor.C
 * @brief  FlowBase MonitorList class
 * @author kero
 */

#include "Monitor.h"

/// Line指定の端点座標をグローバル計算領域内にクリッピング.
///
///   @param[in,out] from Line始点
///   @param[in,out] to   Line終点
///
void MonitorList::clipLine(REAL_TYPE from[3], REAL_TYPE to[3])
{
  const char* OUT_OF_REGION = "out of region";
  FB::Vec3r st(from);
  FB::Vec3r ed(to);
  FB::Vec3r r0 = g_org;
  FB::Vec3r r1 = r0 + g_box;
  FB::Vec3r d = ed - st;
  
  
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
    if (myRank == 0) {
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
  if ((outputType == GATHER && myRank == 0) ||
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
  if ((outputType == GATHER && myRank == 0) ||
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
///   @param [in] bcd BCindex ID
///   @param [in] refVelocity    代表速度
///   @param [in] baseTemp       基準温度
///   @param [in] diffTemp       代表温度差
///   @param [in] refDensity     基準密度
///   @param [in] refLength      代表長さ
///   @param [in] basePrs        基準圧力
///   @param [in] unitTemp       温度単位指定フラグ (Kelvin / Celsius)
///   @param [in] modePrecision  出力精度指定フラグ (単精度，倍精度)
///   @param [in] unitPrs        圧力単位指定フラグ (絶対値，ゲージ圧)
///
void MonitorList::setControlVars(int* bcd,
                                 REAL_TYPE refVelocity, REAL_TYPE baseTemp, REAL_TYPE diffTemp, 
                                 REAL_TYPE refDensity, REAL_TYPE refLength, REAL_TYPE basePrs, 
                                 int unitTemp, int modePrecision, int unitPrs, int num_process) 
{
  this->bcd = bcd;
  this->g_org = G_origin;
  this->g_box = G_region;
  this->org = origin;
  this->pch = deltaX;
  this->box = region;
  this->num_process = num_process;
  
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
  MonitorCompo* m = new MonitorCompo(org, pch, box, g_org, g_box, refVar, bcd, num_process);
  
  m->setRankInfo(paraMngr, procGrp);
  m->setNeighborInfo(guide);
  
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
  
  MonitorCompo* m = new MonitorCompo(org, pch, box, g_org, g_box, refVar, bcd, num_process);
  
  m->setRankInfo(paraMngr, procGrp);
  m->setNeighborInfo(guide);
  
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
    if (cmp[i].isMONITOR()) 
    {
      MonitorCompo* m = new MonitorCompo(org, pch, box, g_org, g_box, refVar, bcd, num_process);
      m->setRankInfo(paraMngr, procGrp);
      m->setNeighborInfo(guide);
      
      m->setInnerBoundary(i, cmp[i]);
      
      monGroup.push_back(m);
    }
  }
  nGroup = monGroup.size();
}



/// 入力ファイルにより指定されるモニタ点位置にID=255を書き込む.
///
///   @param [in] id セルID配列
///
void MonitorList::write_ID(int* id)
{
  // for optimization > variables defined outside
  size_t q;
  int ix, jx, kx, gd;
  ix = size[0];
  jx = size[1];
  kx = size[2];
  gd = guide;
  
  for (int i = 0; i < nGroup; i++) {
    for (int m = 0; m < monGroup[i]->getSize(); m++) {
      if (monGroup[i]->getType() != MonitorCompo::INNER_BOUNDARY) 
      {
        if (!monGroup[i]->check_region(m, org, box)) continue;
        
        FB::Vec3i index = monGroup[i]->getSamplingCellIndex(m);
        q = _F_IDX_S3D(index.x, index.y, index.z, ix, jx, kx, gd);
        id[q] = 255;
      }
    }
  }
}

// TPのポインタを受け取る
void MonitorList::importTP(TPControl* tp) 
{ 
  if ( !tp ) Exit(0);
  tpCntl = tp;
}



// TPに記述されたモニタ座標情報を取得し，リストに保持する
void MonitorList::get_Monitor(Control* C)
{
  REAL_TYPE f_val=0.0;
  MonitorCompo::Type type;
  vector<string> variables;
  
  std::string str,label;
  string label_base,label_leaf;
  string name;
  string method;
  string mode;
  
  REAL_TYPE fval;
  
  // ログ出力のON/OFFはControl::getTP_Sampling()で取得済み
  
  // 集約モード
  label = "/Steer/MonitorList/OutputMode";
  
  if ( !(tpCntl->GetValue(label, &str )) ) 
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    Exit(0);
  }
  
  if ( !strcasecmp(str.c_str(), "gather") ) 
  {
    C->Sampling.out_mode = MonitorList::GATHER;
    setOutputType(MonitorList::GATHER);
  }
  else if( !strcasecmp(str.c_str(), "distribute") ) 
  {
    C->Sampling.out_mode = MonitorList::DISTRIBUTE;
    setOutputType(MonitorList::DISTRIBUTE);
  }
  else 
  {
    Hostonly_ stamped_printf("\tParsing error : Invalid keyord for 'Output_Mode' in 'MonitorList'\n");
    Exit(0);
  }
  
  
  // サンプリング間隔
  label="/Steer/MonitorList/SamplingIntervalType";
  
  if ( !(tpCntl->GetValue(label, &str )) ) 
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    Exit(0);
  }
  else 
  {
    if ( !strcasecmp(str.c_str(), "step") ) 
    {
      C->Interval[Interval_Manager::tg_sampled].setMode_Step();
    }
    else if ( !strcasecmp(str.c_str(), "time") ) 
    {
      C->Interval[Interval_Manager::tg_sampled].setMode_Time();
    }
    else 
    {
      Hostonly_ stamped_printf("\tParsing error : Invalid keyword for '%s'\n", label.c_str());
      Exit(0);
    }
    
    label="/Steer/MonitorList/SamplingInterval";
    
    if ( !(tpCntl->GetValue(label, &f_val )) ) 
    {
      Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
      Exit(0);
    }
    else 
    {
      C->Interval[Interval_Manager::tg_sampled].setInterval((double)f_val);
    }
  }
  
  // 単位指定
  label="/Steer/MonitorList/Unit";
  
  if ( !(tpCntl->GetValue(label, &str )) ) 
  {
    Hostonly_ stamped_printf("\tParsing error : Invalid string for '%s'\n", label.c_str());
    Exit(0);
  }
  
  if ( !strcasecmp(str.c_str(), "Dimensional") ) 
  {
    C->Sampling.unit = DIMENSIONAL;
    setSamplingUnit(DIMENSIONAL);
  }
  else if( !strcasecmp(str.c_str(), "NonDimensional") ) 
  {
    C->Sampling.unit = NONDIMENSIONAL;
    setSamplingUnit(NONDIMENSIONAL);
  }
  else 
  {
    Hostonly_ stamped_printf("\tInvalid keyword is described at '%s'\n", label.c_str());
    Exit(0);
  }
  
  // サンプリングの指定単位が有次元の場合に，無次元に変換
  if ( C->Sampling.unit == DIMENSIONAL ) 
  {
    C->Interval[Interval_Manager::tg_sampled].normalizeInterval(C->Tscale);
  }
  
  
  // 指定モニタ個数のチェック
  int nnode=0;
  int nlist=0;
  label_base = "/Steer/MonitorList";
  
  nnode = tpCntl->countLabels(label_base);
  if ( nnode == 0 ) 
  {
    stamped_printf("\tcountLabels --- %s\n",label_base.c_str());
    Exit(0);
  }
  
  for (int i=0; i<nnode; i++) {
    
    if(!tpCntl->GetNodeStr(label_base, i+1, &str))
    {
      printf("\tParsing error : No No List[@]\n");
      Exit(0);
    }
    
    if( strcasecmp(str.substr(0,4).c_str(), "List") ) continue;
    
    nlist++;
  }
  
  if (nlist==0 && C->isMonitor() == OFF) 
  {
    Hostonly_ stamped_printf("\tError : No monitoring points. Please confirm 'MonitorList' and 'InnerBoundary' in Input parameter file. \n");
    Exit(0);
  }
  
  // モニターリストの読み込み
  label_base = "/Steer/MonitorList";
  
  for (int i=0; i<nnode; i++) 
  {
    if(!tpCntl->GetNodeStr(label_base, i+1, &str))
    {
      printf("\tParsing error : No List[@]\n");
      Exit(0);
    }
    if( strcasecmp(str.substr(0,4).c_str(), "List") ) continue;
    
    
    label_leaf = label_base + "/" + str;
    
    // sampling type & param check
    label = label_leaf + "/Type";
    
    if ( !(tpCntl->GetValue(label, &str )) ) 
    {
      stamped_printf("\tParsing error : No entory '%s'\n", label.c_str());
      Exit(0);
    }
    
    if ( !strcasecmp(str.c_str(), "PointSet") ) 
    {
      type = MonitorCompo::POINT_SET;
    }
    else if ( !strcasecmp(str.c_str(), "Line") ) 
    {
      type = MonitorCompo::LINE;
    }
    else 
    {
      Hostonly_ stamped_printf("\tParsing error : No valid keyword [PointSet / Line] in '%s'\n", label.c_str());
      Exit(0);
    }
    
    // Labelの取得
    label = label_leaf + "/Label";
    
    if ( !(tpCntl->GetValue(label, &name )) ) 
    {
      Hostonly_ stamped_printf("\tParsing warning : No Label in '%s'\n", label.c_str());
      Exit(0);
    }
    
    // variable ---> 複数setできるような仕様にする？
    variables.clear();
    
    label = label_leaf + "/Variable";
    
    if ( !(tpCntl->GetValue(label, &str )) ) 
    {
      Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
      Exit(0);
    }
    variables.push_back(str.c_str());

    
    if (variables.size() == 0) 
    {
      Hostonly_ stamped_printf("\tParsing error : No 'Variable' in '%s'\n", label.c_str());
      Exit(0);
    }
    
    // method
    label = label_leaf + "/SamplingMethod";
    
    if ( !(tpCntl->GetValue(label, &method )) ) 
    {
      Hostonly_ stamped_printf("\tParsing error : fail to get 'SamplingMethod' in '%s'\n", label.c_str());
      Exit(0);
    }
    
    // mode
    label = label_leaf + "/SamplingMode";
    
    if ( !(tpCntl->GetValue(label, &mode )) ) 
    {
      Hostonly_ stamped_printf("\tParsing error : fail to get 'SamplingMode' in '%s'\n", label.c_str());
      Exit(0);
    }
    
    // get coordinate
    if ( type == MonitorCompo::POINT_SET )
    {
      vector<MonitorCompo::MonitorPoint> pointSet;
      get_Mon_Pointset(C, label_leaf, pointSet);
      setPointSet(name.c_str(), variables, method.c_str(), mode.c_str(), pointSet);
    }
    else 
    {
      REAL_TYPE from[3], to[3];
      int nDivision;
      get_Mon_Line(C, label_leaf, from, to, nDivision);
      setLine(name.c_str(), variables, method.c_str(), mode.c_str(), from, to, nDivision);
    }
  }
  
  
}


// TPに記述されたモニタ座標情報(Line)を取得
void MonitorList::get_Mon_Line(Control* C,
                               const string label_base,
                               REAL_TYPE from[3],
                               REAL_TYPE to[3],
                               int& nDivision)
{
  std::string str,label;
  REAL_TYPE v[3];
  
  label=label_base+"/Division";

  if ( !(tpCntl->GetValue(label, &nDivision )) ) 
  {
	  Hostonly_ stamped_printf("\tParsing error : No Division\n");
	  Exit(0);
  }
  if ( nDivision == 0 ) Exit(0);
  
  // load parameter of 'from' and 'to'
  label=label_base+"/From";

  for (int n=0; n<3; n++) v[n]=0.0;
  if ( !(tpCntl->GetVector(label, v, 3 )) ) 
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get 'From' in 'Line'\n");
    Exit(0);
  }
  from[0]=v[0];
  from[1]=v[1];
  from[2]=v[2];
  if (C->Sampling.unit == DIMENSIONAL) 
  {
    normalizeCord(C->RefLength,from);
  }
  
  label=label_base+"/To";

  for (int n=0; n<3; n++) v[n]=0.0;
  if ( !(tpCntl->GetVector(label, v, 3 )) ) 
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get 'To' in 'Line'\n");
    Exit(0);
  }
  to[0]=v[0];
  to[1]=v[1];
  to[2]=v[2];
  if (C->Sampling.unit == DIMENSIONAL) 
  {
    normalizeCord(C->RefLength,to);
  }
}

// TPに記述されたモニタ座標情報を取得(PointSet)
void MonitorList::get_Mon_Pointset(Control* C,
                                   const string label_base,
                                   vector<MonitorCompo::MonitorPoint>& pointSet)
{
  REAL_TYPE v[3];
  char tmpstr[20];
  std::string str,label;
  string label_leaf;
  
  // PointSet個数のチェック
  int nnode=0;
  int nlist=0;
  
  nnode=tpCntl->countLabels(label_base);
  if ( nnode == 0 )
  {
    stamped_printf("\tcountLabels --- %s\n",label_base.c_str());
    Exit(0);
  }
  
  for (int i=0; i<nnode; i++) {
    if(!tpCntl->GetNodeStr(label_base,i+1,&str)){
      printf("\tParsing error : No Elem name\n");
      Exit(0);
    }
    if( strcasecmp(str.substr(0,3).c_str(), "set") ) continue;
    nlist++;
  }
  
  // PointSet取得
  int pc=0;
  for (int i=0; i<nnode; i++)
  {
    if(!tpCntl->GetNodeStr(label_base,i+1,&str))
    {
      printf("\tParsing error : No Elem name\n");
      Exit(0);
    }
    if( strcasecmp(str.substr(0,3).c_str(), "set") ) continue;
    pc++;
    
    label_leaf = label_base + "/" + str;
    
    // set coordinate
    label=label_leaf+"/Coordinate";
    for (int n=0; n<3; n++) v[n]=0.0;
    if ( !(tpCntl->GetVector(label, v, 3 )) )
    {
      Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
      Exit(0);
    }
    if (C->Sampling.unit == DIMENSIONAL) {
      normalizeCord(C->RefLength,v);
    }
    
    // set Labelの取得．ラベルなしでもエラーではない
    label=label_leaf+"/Label";
    if ( !(tpCntl->GetValue(label, &str )) ) {
      Hostonly_ stamped_printf("\tParsing warning : No commnet for '%s'\n", label.c_str());
    }
    if ( !strcasecmp(str.c_str(), "") ) {
      sprintf(tmpstr, "point_%d",pc);
      str = tmpstr;
    }
    
    pointSet.push_back(MonitorCompo::MonitorPoint(v, str.c_str()));
  }
  
}
