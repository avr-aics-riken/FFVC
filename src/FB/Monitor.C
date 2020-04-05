//##################################################################################
//
// Flow Base class
//
// Copyright (c) 2007-2011 VCAD System Research Program, RIKEN.
// All rights reserved.
//
// Copyright (c) 2011-2015 Institute of Industrial Science, The University of Tokyo.
// All rights reserved.
//
// Copyright (c) 2012-2016 Advanced Institute for Computational Science, RIKEN.
// All rights reserved.
//
// Copyright (c) 2016-2020 Research Institute for Information Research, Kyushu University.
// All rights reserved.
//
//##################################################################################

/*
 * @file   Monitor.C
 * @brief  FlowBase MonitorList class
 * @author aics
 */

#include "Monitor.h"

// #################################################################
/// サンプリングの変数指定と計算モードの整合性をチェックする
bool MonitorList::checkConsistency(const int KOS)
{
  int count[var_END];
  
  for (int k=var_Velocity; k<var_END; k++)
  {
    count[k]=0;
  }
  
  // モニタで使用する変数を取得
  for (int k=var_Velocity; k<var_END; k++)
  {
    for (int i = 0; i < nGroup; i++)
    {
      if (monGroup[i]->getStateVariable(k) == true) count[k]++;
    }
  }
  
  
  // 利用しない変数のサンプリングは中止
  switch (KOS)
  {
    case COLD_FLOW:
      if ( count[var_Temperature]>0 )
      {
        Hostonly_{
          printf("\tError : Temperature is not able to sample in 'COLD_FLOW' mode\n");
        }
        return false;
      }
      break;
      
    case SOLID_CONDUCTION:
      if ( (count[var_Velocity]>0) ||
           (count[var_Pressure]>0) ||
           (count[var_TotalP]>0)   ||
           (count[var_Helicity]>0) ||
           (count[var_Vorticity]>0) )
      {
        Hostonly_{
          printf("\tError : In 'SOLID_CONDUCTION' mode, sampling variable is temperature only\n");
        }
        return false;
      }
      break;
  };

  return true;
}



// #################################################################
/// Polygonモニターの交点の状況を確認
void MonitorList::checkStatus()
{
  for (int i = 0; i < nGroup; i++)
  {
    monGroup[i]->checkMonitorPoints();
  }
}


// #################################################################
/// Line指定の端点座標をグローバル計算領域内にクリッピング
///
///   @param [in,out] from Line始点
///   @param [in,out] to   Line終点
///
void MonitorList::clipLine(REAL_TYPE from[3], REAL_TYPE to[3])
{
  //>> Graph Ploter
  int dbg_pos = -1;
  const char* OUT_OF_REGION1 = "out of region 1";
  const char* OUT_OF_REGION2 = "out of region 2";
  const char* OUT_OF_REGION3 = "out of region 3";
  const char* OUT_OF_REGION4 = "out of region 4";
  const char* OUT_OF_REGION5 = "out of region 5";
  const char* OUT_OF_REGION6 = "out of region 6";
  const char* OUT_OF_REGION7 = "out of region 7";
  const char* OUT_OF_REGION8 = "out of region 8";
  const char* OUT_OF_REGION9 = "out of region 9";
  const char* OUT_OF_REGION0 = "out of region 0";
  //<< Graph Ploter

  Vec3r st(from);
  Vec3r ed(to);
  
  st.x=from[0]; st.y=from[1]; st.z=from[2];
  ed.x=to[0];   ed.y=to[1];   ed.z=to[2];
  
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
      if (st.x < r0.x || r1.x < st.x) throw OUT_OF_REGION1;
    }
    else if (d.x > 0) {
      if (ed.x < r0.x || r1.x < st.x) throw OUT_OF_REGION2;
      t_st = (std::max)(t_st, (r0.x-st.x)/d.x);
      t_ed = (std::min)(t_ed, (r1.x-st.x)/d.x);
      dbg_pos = 1;
    }
    else if (d.x < 0) {
      if (st.x < r0.x || r1.x < ed.x) throw OUT_OF_REGION3;
      t_st = (std::max)(t_st, (r1.x-st.x)/d.x);
      t_ed = (std::min)(t_ed, (r0.x-st.x)/d.x);
      dbg_pos = 2;
    }
    
    if (d.y == 0) {
      if (st.y < r0.y || r1.y < st.y) throw OUT_OF_REGION4;
    }
    else if (d.y > 0) {
      if (ed.y < r0.y || r1.y < st.y) throw OUT_OF_REGION5;
      t_st = (std::max)(t_st, (r0.y-st.y)/d.y);
      t_ed = (std::min)(t_ed, (r1.y-st.y)/d.y);
      dbg_pos = 3;
    }
    else if (d.y < 0) {
      if (st.y < r0.y || r1.y < ed.y) throw OUT_OF_REGION6;
      t_st = (std::max)(t_st, (r1.y-st.y)/d.y);
      t_ed = (std::min)(t_ed, (r0.y-st.y)/d.y);
      dbg_pos = 4;
    }
    
    if (d.z == 0) {
      if (st.z < r0.z || r1.z < st.z) throw OUT_OF_REGION7;
    }
    else if (d.z > 0) {
      if (ed.z < r0.z || r1.z < st.z) throw OUT_OF_REGION8;
      t_st = (std::max)(t_st, (r0.z-st.z)/d.z);
      t_ed = (std::min)(t_ed, (r1.z-st.z)/d.z);
      dbg_pos = 5;
    }
    else if (d.z < 0) {
      if (st.z < r0.z || r1.z < ed.z) throw OUT_OF_REGION9;
      t_st = (std::max)(t_st, (r1.z-st.z)/d.z);
      t_ed = (std::min)(t_ed, (r0.z-st.z)/d.z);
      dbg_pos = 6;
    }
    
    if (t_st >= t_ed) throw OUT_OF_REGION0;
  }
  catch (const char *str)
  {
    if (myRank == 0)
    {
      printf("Line [%e %e %e]-[%e %e %e] is %s\n",
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



// #################################################################
/// Polygonモニターの交点と境界IDの除去
void MonitorList::clearCut()
{
  for (int i = 0; i < nGroup; i++)
  {
    if ( monGroup[i]->getType() == mon_POLYGON )
    {
      unsigned long c = monGroup[i]->clearMonitorCut();
      
      Hostonly_
      {
        printf("\tPolygon Monitor : cleared cell = %ld\n", c);
      }
    }
  }
}


// #################################################################
/// 出力ファイルクローズ
///
///    @note 出力タイプによらず全プロセスから呼んでも問題ない
///
void MonitorList::closeFile()
{
  if ((outputType == GATHER && myRank == 0) || outputType == DISTRIBUTE)
  {
    for (int i = 0; i < nGroup; i++)
    {
      monGroup[i]->closeFile();
    }
  }
}


// #################################################################
/**
 * @brief プリミティブ形状を作成 >> bcd[]にIDをペイント
 * @param [in]  odr       ID番号
 * @param [in]  montyp    モニタのタイプ
 * @param [in]  nv        法線ベクトル
 * @param [in]  ctr       中心座標
 * @param [in]  depth     厚さ
 * @param [in]  tmp       幅 / 半径
 * @param [in]  height    高さ
 * @param [in]  dr        方向ベクトル
 */
void MonitorList::generatePrimitiveShape(const int odr,
                                         const Monitor_Type montyp,
                                         REAL_TYPE nv[3],
                                         REAL_TYPE ctr[3],
                                         REAL_TYPE depth,
                                         REAL_TYPE tmp,
                                         REAL_TYPE height,
                                         REAL_TYPE dr[3])
{
  ctr[0] /= refVar.refLength;
  ctr[1] /= refVar.refLength;
  ctr[2] /= refVar.refLength;
  
  depth  /= refVar.refLength;
  tmp    /= refVar.refLength;
  height /= refVar.refLength;
  
  // ShapeMonitorのインスタンス >> float引数
  ShapeMonitor SM(size, guide, pitch, origin);
  
  
  switch ( montyp )
  {
    case mon_BOX:
      SM.setShapeParam(nv, ctr, dr, depth, tmp, height);
      break;
      
    case mon_CYLINDER:
      SM.setShapeParam(nv, ctr, depth, tmp);
      break;
      
    default:
      Exit(0);
      break;
  }
  
  // 回転角度の計算
  SM.getAngle();
  
  
  // bboxと投影面積の計算 [m^2]とインデクスの計算
  int f_st[3], f_ed[3];
  area = SM.getBboxArea(f_st, f_ed) * refVar.refLength * refVar.refLength;
  //printf("area = %f\n", area);
  //printf("(%d %d %d) - (%d %d %d)\n", f_st[0], f_st[1], f_st[2], f_ed[0], f_ed[1], f_ed[2]);
  
  SM.setID(f_st, f_ed, bcd, odr);
}


// #################################################################
/**
 * @brief モニタ座標情報(Line)を取得
 * @param [in]  C         Control クラスオブジェクトのポインタ
 * @param [out] from      Line始点座標
 * @param [out] to        Line終点座標
 * @param [out] nDivision Line分割数
 * @note データは無次元化して保持
 */
void MonitorList::getLine(const Control* C,
                          const string label_base,
                          REAL_TYPE from[3],
                          REAL_TYPE to[3],
                          int& nDivision)
{
  string str,label;
  REAL_TYPE v[3];
  
  label=label_base+"/Division";
  
  if ( !(tpCntl->getInspectedValue(label, nDivision )) )
  {
	  Hostonly_ stamped_printf("\tParsing error : No Division\n");
	  Exit(0);
  }
  if ( nDivision == 0 ) Exit(0);
  
  // load parameter of 'from' and 'to'
  label=label_base+"/From";
  
  for (int n=0; n<3; n++) v[n]=0.0;
  if ( !(tpCntl->getInspectedVector(label, v, 3 )) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get 'From' in 'Line'\n");
    Exit(0);
  }
  from[0]=v[0];
  from[1]=v[1];
  from[2]=v[2];
  
  // 入力パラメータの次元が有次元のとき，無次元化する
  if (C->Unit.Output == DIMENSIONAL) normalizeCord(from);

  
  label=label_base+"/To";
  
  for (int n=0; n<3; n++) v[n]=0.0;
  if ( !(tpCntl->getInspectedVector(label, v, 3 )) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get 'To' in 'Line'\n");
    Exit(0);
  }
  to[0]=v[0];
  to[1]=v[1];
  to[2]=v[2];
  
  // 入力パラメータの次元が有次元のとき，無次元化する
  if (C->Unit.Output == DIMENSIONAL) normalizeCord(to);

}



// #################################################################
// モニタ座標情報を取得し，リストに保持する
bool MonitorList::getMonitor(Control* C, CompoList* cmp)
{
  Monitor_Type mon_type;
  vector<string> variables;
  
  string str, label;
  string label_base,label_leaf;
  string name;
  string method;
  string mode;
  
  // ログ出力
  label = "/MonitorList/log";
  
  if ( !(tpCntl->getInspectedValue(label, str )) )
  {
	  Hostonly_ stamped_printf("\tParsing error : Invalid string for '%s'\n", label.c_str());
	  Exit(0);
  }
  
  if     ( !strcasecmp(str.c_str(), "on") )   C->SamplingMode = ON;
  else if( !strcasecmp(str.c_str(), "off") )  C->SamplingMode = OFF;
  else
  {
    Hostonly_ stamped_printf("\tInvalid keyword is described for '%s'\n", label.c_str());
    Exit(0);
  }
  
  if (C->SamplingMode == OFF)
  {
    return false;
  }
  
  
  
  // 集約モード
  label = "/MonitorList/OutputMode";
  
  if ( !(tpCntl->getInspectedValue(label, str )) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    Exit(0);
  }
  
  if ( !strcasecmp(str.c_str(), "gather") )
  {
    setOutputType(MonitorList::GATHER);
  }
  else if( !strcasecmp(str.c_str(), "distribute") )
  {
    setOutputType(MonitorList::DISTRIBUTE);
  }
  else
  {
    Hostonly_ stamped_printf("\tParsing error : Invalid keyord for 'OutputMode' in 'MonitorList'\n");
    Exit(0);
  }
  
  
  // サンプリング間隔
  label="/MonitorList/Sampling/TemporalType";
  
  if ( !(tpCntl->getInspectedValue(label, str )) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    Exit(0);
  }
  else
  {
    if ( !strcasecmp(str.c_str(), "step") )
    {
      C->Interval[Control::tg_sampled].setMode(IntervalManager::By_step);
    }
    else if ( !strcasecmp(str.c_str(), "time") )
    {
      C->Interval[Control::tg_sampled].setMode(IntervalManager::By_time);
    }
    else
    {
      Hostonly_ stamped_printf("\tParsing error : Invalid keyword for '%s'\n", label.c_str());
      Exit(0);
    }
    
    label="/MonitorList/Sampling/Interval";
    double f_val=0.0;
    
    if ( !(tpCntl->getInspectedValue(label, f_val )) )
    {
      Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
      Exit(0);
    }
    else
    {
      C->Interval[Control::tg_sampled].setInterval(f_val);
    }
  }

  
  // 出力と入力の単位モード 
  setSamplingUnit(C->Unit.Output);
  
  
  
  // 指定モニタ個数のチェック
  int nnode=0;
  int nlist=0;
  label_base = "/MonitorList";
  
  nnode = tpCntl->countLabels(label_base);
  if ( nnode == 0 )
  {
    stamped_printf("\tcountLabels --- %s\n", label_base.c_str());
    Exit(0);
  }
  
  for (int i=0; i<nnode; i++)
  {
    if ( !(tpCntl->getNodeStr(label_base, i+1, str)) )
    {
      printf("\tParsing error : No No List[@]\n");
      Exit(0);
    }
    
    if( strcasecmp(str.substr(0,4).c_str(), "List") ) continue;
    
    nlist++;
  }
  
  if ( nlist==0 )
  {
    Hostonly_ stamped_printf("\tError : No monitoring points. Please confirm 'MonitorList' in Input parameter file. \n");
  }
  
  
  // モニターリストの読み込み
  label_base = "/MonitorList";
  
  for (int i=0; i<nnode; i++)
  {
    if ( !(tpCntl->getNodeStr(label_base, i+1, str)) )
    {
      printf("\tParsing error : No List[@]\n");
      Exit(0);
    }
    if( strcasecmp(str.substr(0,4).c_str(), "List") ) continue;
    
    
    label_leaf = label_base + "/" + str;
    
    // sampling type & param check
    label = label_leaf + "/Type";
    
    if ( !(tpCntl->getInspectedValue(label, str )) )
    {
      stamped_printf("\tParsing error : No entory '%s'\n", label.c_str());
      Exit(0);
    }
    
    if ( !strcasecmp(str.c_str(), "PointSet") )
    {
      mon_type = mon_POINT_SET;
    }
    else if ( !strcasecmp(str.c_str(), "Line") )
    {
      mon_type = mon_LINE;
    }
    else if ( !strcasecmp(str.c_str(), "Cylinder") )
    {
      mon_type = mon_CYLINDER;
    }
    else if ( !strcasecmp(str.c_str(), "Box") )
    {
      mon_type = mon_BOX;
    }
    else if ( !strcasecmp(str.c_str(), "Polygon") )
    {
      mon_type = mon_POLYGON;
    }
    else if ( !strcasecmp(str.c_str(), "Plane") )
    {
      mon_type = mon_PLANE;
    }
    else
    {
      Hostonly_ stamped_printf("\tParsing error : Invalid keyword in label='%s', str=%s\n", label.c_str(), str.c_str());
      Exit(0);
    }
    
    
    // Labelの取得
    label = label_leaf + "/Label";
    
    if ( !(tpCntl->getInspectedValue(label, name )) )
    {
      Hostonly_ stamped_printf("\tParsing warning : No Label in '%s'\n", label.c_str());
      Exit(0);
    }
    
    
    // Variables
    label = label_leaf + "/Variables";
    
    if ( !(tpCntl->chkNode(label)) )
    {
      Hostonly_ stamped_printf("\tParsing error : No 'Variables' keyword in '%s'\n", label.c_str());
      Exit(0);
    }
    
    // 登録変数の数
    int n_vars = 0;
    n_vars = tpCntl->countLabels(label);
    
    // モニタする変数名と数を取得
    variables.clear();
    
    for (int nc=1; nc<=n_vars; nc++)
    {
      registVars(label, variables, "Velocity");
      registVars(label, variables, "Pressure");
      

      if ( C->isHeatProblem() )
      {
        registVars(label, variables, "Temperature");
      }
      
      registVars(label, variables, "TotalPressure");
      registVars(label, variables, "Helicity");
      registVars(label, variables, "Vorticity");
    }
    
    
    // variable
    if (variables.size() == 0)
    {
      Hostonly_ stamped_printf("\tParsing error : No 'Variable' in '%s'\n", label.c_str());
      Exit(0);
    }
    
    
    // method
    label = label_leaf + "/SamplingMethod";
    
    if ( !(tpCntl->getInspectedValue(label, method )) )
    {
      Hostonly_ stamped_printf("\tParsing error : fail to get 'SamplingMethod' in '%s'\n", label.c_str());
      Exit(0);
    }
    
    // mode
    label = label_leaf + "/SamplingMode";
    
    if ( !(tpCntl->getInspectedValue(label, mode )) )
    {
      Hostonly_ stamped_printf("\tParsing error : fail to get 'SamplingMode' in '%s'\n", label.c_str());
      Exit(0);
    }

    
    // set coordinate values
    if ( mon_type == mon_POINT_SET )
    {
      vector<MonitorCompo::MonitorPoint> pointSet;
      getPointset(C, label_leaf, pointSet);
      MonitorCompo *m = setPointSet(name.c_str(), variables, method.c_str(), mode.c_str(), pointSet, mon_type);
      if( m != NULL )
      {
        m->setObjType(mon_POINT_SET);
      }
    }
    
    else if ( mon_type == mon_LINE )
    {
      REAL_TYPE from[3], to[3];
      int nDivision;
      getLine(C, label_leaf, from, to, nDivision);
      MonitorCompo *m = setLine(name.c_str(), variables, method.c_str(), mode.c_str(), from, to, nDivision, mon_type);
      if( m != NULL )
      {
        m->setObjType(mon_LINE);
      }

    }
    
    else if( mon_type == mon_PLANE )
    {
      REAL_TYPE c[3], z[3], x[3], dim[2], div[2];
      getPlane( mon_type, label_leaf, c, z, x, dim, div );
      
      vector<MonitorCompo::MonitorPoint> pointSet;
      vector<int> grid_flags;
      
      Vec3r orig_o( c[0], c[1], c[2] );
      Vec3r axis_z( z[0], z[1], z[2] );
      Vec3r axis_x( x[0], x[1], x[2] );
      getPointsOnPlane(orig_o, axis_z, axis_x, dim, div, &pointSet, &grid_flags);
      
      MonitorCompo *m = setPointSet(name.c_str(), variables, method.c_str(), mode.c_str(), pointSet, mon_POINT_SET);
      
      if( m != NULL )
      {
        m->setPlaneData(c, z, x, dim, div);
        m->setPlaneFlags(&grid_flags);
        m->setObjType(mon_PLANE);
      }
    }
    
    
    else if ( mon_type == mon_POLYGON )
    {
      // nameに対するエントリ番号の取得
      int odr = -1;
      
      for (int k=1; k<=C->NoCompo; k++)
      {
        if ( !strcasecmp( cmp[k].alias.c_str(), name.c_str()) )
        {
          odr = k;
          break;
        }
      }
      
      if ( (odr < 1) || (odr > C->NoCompo) )
      {
        Hostonly_ stamped_printf("\tSomthing wrong %d\n", odr);
        Exit(0);
      }

#if 0
      printf("%s %d\n", name.c_str(), odr);
#endif
      
      // 法線ベクトル
      REAL_TYPE nv[3];
      label = label_leaf + "/Normal";
      if ( !Control::getVec(label, nv, tpCntl, true) ) Exit(0);
      
      setPolygon(name.c_str(), variables, method.c_str(), mode.c_str(), odr, nv, mon_type);
      setOutputType(MonitorList::GATHER); // 強制
    }
    
    else if ( (mon_type == mon_CYLINDER) || (mon_type == mon_BOX) )
    {
      REAL_TYPE nv[3], ctr[3], dir[3];
      REAL_TYPE depth, tmp, height;
      getPrimitive(mon_type, label_leaf, nv, ctr, depth, tmp, height, dir);
      generatePrimitiveShape(i+1, mon_type, nv, ctr, depth, tmp, height, dir);
      setPrimitive(name.c_str(), variables, method.c_str(), mode.c_str(), i+1, nv, mon_type);
      setOutputType(MonitorList::GATHER); // 強制
    }
    
    else
    {
      Exit(0);
    }
  }
  
  nGroup = (int)monGroup.size();
  
  return true;
}


// #################################################################
/// 出力タイプ文字列の取得
string MonitorList::getOutputTypeStr()
{
  string str;
  
  switch (outputType)
  {
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



// #################################################################
/**
 * @brief モニタ座標情報を取得(PointSet)
 * @param [in]     C          Control クラスオブジェクトのポインタ
 * @param [in]     label_base テキストパーサーのノード
 * @param [in,out] pointSet   PointSet配列
 * @note データは無次元化して保持
 */
void MonitorList::getPointset(const Control* C,
                              const string label_base,
                              vector<MonitorCompo::MonitorPoint>& pointSet)
{
  char tmpstr[20];
  string str,label;
  string label_leaf;
  
  // PointSet個数のチェック
  int nnode=0;
  int nlist=0;
  
  nnode=tpCntl->countLabels(label_base);
  
  if ( nnode == 0 )
  {
    stamped_printf("\tcountLabels --- %s\n", label_base.c_str());
    Exit(0);
  }
  
  for (int i=0; i<nnode; i++)
  {
    if ( !(tpCntl->getNodeStr(label_base, i+1, str)) )
    {
      printf("\tParsing error : No node name\n");
      Exit(0);
    }
    
    if ( strcasecmp(str.substr(0,3).c_str(), "set") ) continue;
    nlist++;
  }
  
  // PointSet取得
  int pc=0;
  for (int i=0; i<nnode; i++)
  {
    if ( !(tpCntl->getNodeStr(label_base, i+1, str)) )
    {
      printf("\tParsing error : No Elem name\n");
      Exit(0);
    }
    
    if ( strcasecmp(str.substr(0,3).c_str(), "set") ) continue;
    pc++;
    
    label_leaf = label_base + "/" + str;
    
    // set coordinate
    label = label_leaf + "/Coordinate";
    REAL_TYPE v[3];
    for (int n=0; n<3; n++) v[n]=0.0;
    
    if ( !(tpCntl->getInspectedVector(label, v, 3 )) )
    {
      Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
      Exit(0);
    }
    
    // 入力パラメータの次元が有次元のとき，無次元化する
    if (C->Unit.Param == DIMENSIONAL) normalizeCord(v);

    
    // set tagの取得．ラベルなしでもエラーではない
    label = label_leaf + "/tag";
    
    if ( !(tpCntl->getInspectedValue(label, str)) )
    {
      Hostonly_ stamped_printf("\tParsing warning : No tag for '%s'\n", label.c_str());
    }
    if ( !strcasecmp(str.c_str(), "") )
    {
      sprintf(tmpstr, "point_%d", pc);
      str = tmpstr;
    }
    
    pointSet.push_back(MonitorCompo::MonitorPoint(v, str.c_str()));
  }
  
}


// #################################################################
/**
 * @brief プリミティブのモニタ形状情報を取得
 * @param [in]  mon_type   モニタのタイプ
 * @param [in]  label_base テキストパーサーのノード
 * @param [out] nv         法線ベクトル
 * @param [out] center     中心座標
 * @param [out] depth      厚さ
 * @param [out] tmp        幅 / 半径
 * @param [out] height     高さ
 * @param [out] dir        方向ベクトル
 */
void MonitorList::getPrimitive(Monitor_Type mon_type,
                               const string label_base,
                               REAL_TYPE nv[3],
                               REAL_TYPE center[3],
                               REAL_TYPE& depth,
                               REAL_TYPE& tmp,
                               REAL_TYPE& height,
                               REAL_TYPE dir[3])
{
  string label;
  
  // 法線ベクトル
  label = label_base + "/Normal";
  if ( !Control::getVec(label, nv, tpCntl, true) ) Exit(0);
  
  
  // 中心座標の取得
  label = label_base + "/Center";
  if ( !Control::getVec(label, center, tpCntl, false) ) Exit(0);
  
  label = label_base + "/Depth";
  depth = Control::getValueReal(label, tpCntl);
  
  
  switch (mon_type)
  {
    case mon_BOX:
      label = label_base + "/OrientationVector";
      if ( !Control::getVec(label, dir, tpCntl, true) ) Exit(0);
      
      label = label_base + "/Width";
      tmp = Control::getValueReal(label, tpCntl);
      
      label = label_base + "/Height";
      height= Control::getValueReal(label, tpCntl);
      break;
      
    case mon_CYLINDER:
      label = label_base + "/Radius";
      tmp = Control::getValueReal(label, tpCntl);
      break;
      
    default:
      break;
  }
  
}

// #################################################################
// Graph Ploter
void MonitorList::getBox(Monitor_Type mon_type,
                         const string label_base,
                         REAL_TYPE center[3],
                         REAL_TYPE main_vec[3],
                         REAL_TYPE ref_vec[3],
                         REAL_TYPE dim_vec[3])
{
  string label;
  
  // 中心座標の取得
  label = label_base + "/Center";
  if ( !Control::getVec(label, center, tpCntl, false) ) Exit(0);
  
  // 主座標軸ベクトル　ローカル　Ｚ
  label = label_base + "/MainDirection";
  if ( !Control::getVec(label, main_vec, tpCntl, true) ) Exit(0);
  
  // 参考座標軸ベクトル　ローカル　Ｘ
  label = label_base + "/RefDirection";
  if ( !Control::getVec(label, ref_vec, tpCntl, true) ) Exit(0);
  
  // 箱のサイズ　ローカル座標系　ｄｘ、ｄｙ、ｄｚ
  label = label_base + "/Dimension";
  if ( !Control::getVec(label, dim_vec, tpCntl, false) ) Exit(0);
}

// #################################################################
// Graph Ploter
void MonitorList::getPointsInBox( Vec3r orig_o, Vec3r axis_z, Vec3r axis_x, REAL_TYPE dim[3], vector<MonitorCompo::MonitorPoint> *pointSet )
{
  char   tmpstr[256]="";
  
  REAL_TYPE x0 = origin[0],  y0 = origin[1],     z0 = origin[2];
  REAL_TYPE dx = pitch[0],   dy = pitch[1],      dz = pitch[2];
  REAL_TYPE hdx= dx*0.5,     hdy= dy*0.5,        hdz= dz*0.5;
  
  Vec3r axis_y = cross(axis_z, axis_x);
  axis_y.normalize();
  
  int ii=0;
  for( int _k=0; _k<size[2]; _k++ )
  {
    REAL_TYPE _z = z0 + _k*dz + hdz;
    
    for( int _j=0; _j<size[1]; _j++ )
    {
      REAL_TYPE _y = y0 + _j*dy + hdy;
      
      for( int _i=0; _i<size[0]; _i++ )
      {
        REAL_TYPE _x = x0 + _i*dx + hdx;
        
        Vec3r global_pt(_x, _y, _z);
        
        Vec3r local_pt = localPt(orig_o, axis_z, axis_x, axis_y, global_pt);
        
        if( isInBox(local_pt, dim) == true )
        {
          REAL_TYPE v[3] = {_x, _y, _z};
          sprintf(tmpstr, "point_in_box_%d", ii++);
          pointSet->push_back( MonitorCompo::MonitorPoint(v, (const char*)tmpstr) );
        }
      }
    }
  }
  
  return;
}

// #################################################################
// Graph Ploter
void MonitorList::getPointsInCyl( Vec3r orig_o, Vec3r axis_z, Vec3r axis_x, REAL_TYPE dim[3], vector<MonitorCompo::MonitorPoint> *pointSet )
{
  char   tmpstr[256]="";
  
  REAL_TYPE x0 = origin[0],  y0 = origin[1],     z0 = origin[2];
  REAL_TYPE dx = pitch[0],   dy = pitch[1],      dz = pitch[2];
  REAL_TYPE hdx= dx*0.5,     hdy= dy*0.5,        hdz= dz*0.5;
  
  Vec3r axis_y = cross(axis_z, axis_x);
  axis_y.normalize();
  
  int ii=0;
  for( int _k=0; _k<size[2]; _k++ )
  {
    REAL_TYPE _z = z0 + _k*dz + hdz;
    
    for( int _j=0; _j<size[1]; _j++ )
    {
      REAL_TYPE _y = y0 + _j*dy + hdy;
      
      for( int _i=0; _i<size[0]; _i++ )
      {
        REAL_TYPE _x = x0 + _i*dx + hdx;
        
        Vec3r global_pt(_x, _y, _z);
        
        Vec3r local_pt = localPt(orig_o, axis_z, axis_x, axis_y, global_pt);
        
        if( isInCyl(local_pt, dim) == true )
        {
          REAL_TYPE v[3] = {_x, _y, _z};
          sprintf(tmpstr, "point_in_cyl_%d", ii++);
          pointSet->push_back( MonitorCompo::MonitorPoint(v, (const char*)tmpstr) );
        }
      }
    }
  }
  
  return;
}

// #################################################################
// Graph Ploter
void MonitorList::getPointsOnPlane( Vec3r orig_o, Vec3r axis_z, Vec3r axis_x, REAL_TYPE dim[2], REAL_TYPE div[2],
                                   vector<MonitorCompo::MonitorPoint> *pointSet, std::vector<int> *grid_flags )
{
  char   tmpstr[256]="";
  
  int m = (int)div[0];//分割数
  int n = (int)div[1];//分割数
  
  REAL_TYPE x0 = - dim[0] * 0.5;
  REAL_TYPE y0 = - dim[1] * 0.5;
  REAL_TYPE dx = dim[0] / m;
  REAL_TYPE dy = dim[1] / n;
  
  Vec3r axis_y = cross(axis_z, axis_x);
  axis_y.normalize();
  
  int jj=0;
  for( int j=0; j<n+1; j++ )
  {
    REAL_TYPE y = y0 + j * dy;
    
    for( int i=0; i<m+1; i++ )
    {
      REAL_TYPE x = x0 + i * dx;
      
      Vec3r local_pt(x, y, 0.0);
      Vec3r global_pt = globalPt(orig_o, axis_z, axis_x, axis_y, local_pt);
      
      REAL_TYPE lp[3] = {local_pt.x, local_pt.y, local_pt.z};
      REAL_TYPE gp[3] = {global_pt.x, global_pt.y, global_pt.z};
      
      int  index = -1;
      
      if ((global_pt.x < origin[0]) || (global_pt.x > (origin[0]+region[0]))  ||
          (global_pt.y < origin[1]) || (global_pt.y > (origin[1]+region[1]))  ||
          (global_pt.z < origin[2]) || (global_pt.z > (origin[2]+region[2]))  )
      {
        index = -1;
      }
      else
      {
        index = jj;
        
        sprintf(tmpstr, "point_in_pln_%d", index);
        pointSet->push_back( MonitorCompo::MonitorPoint(gp, (const char*)tmpstr) );
        
        jj++;
      }
      
      grid_flags->push_back(index);
    }
  }
}

// #################################################################
// Graph Ploter
void MonitorList::getCylinder(Monitor_Type mon_type,
                              const string label_base,
                              REAL_TYPE center[3],
                              REAL_TYPE main_vec[3],
                              REAL_TYPE ref_vec[3],
                              REAL_TYPE dim_vec[3])
{
  string label;
  
  // 中心座標の取得
  label = label_base + "/Center";
  if ( !Control::getVec(label, center, tpCntl, false) ) Exit(0);
  
  // 主座標軸ベクトル　ローカル　Ｚ
  label = label_base + "/MainDirection";
  if ( !Control::getVec(label, main_vec, tpCntl, true) ) Exit(0);
  
  // 参考座標軸ベクトル　ローカル　Ｘ
  label = label_base + "/RefDirection";
  if ( !Control::getVec(label, ref_vec, tpCntl, true) ) Exit(0);
  
  // 円柱（円錐）の寸法　dim_vec[0]---底面半径、dim_vec[1]---上面半径、dim_vec[3]---高さ
  label = label_base + "/Dimension";
  if ( !Control::getVec(label, dim_vec, tpCntl, false) ) Exit(0);
}

void MonitorList::getPlane(Monitor_Type mon_type,
                           const string label_base,
                           REAL_TYPE center[3],
                           REAL_TYPE main_vec[3],
                           REAL_TYPE ref_vec[3],
                           REAL_TYPE dim_vec[2],
                           REAL_TYPE div_vec[2])
{
  string label;
  
  // 中心座標の取得
  label = label_base + "/Center";
  if ( !Control::getVec(label, center, tpCntl, false) ) Exit(0);
  
  // 主座標軸ベクトル　ローカル　Ｚ
  label = label_base + "/MainDirection";
  if ( !Control::getVec(label, main_vec, tpCntl, true) ) Exit(0);
  
  // 参考座標軸ベクトル　ローカル　Ｘ
  label = label_base + "/RefDirection";
  if ( !Control::getVec(label, ref_vec, tpCntl, true) ) Exit(0);
  
  // ローカル座標系　ｄｘ、ｄｙ
  label = label_base + "/Dimension";
  if ( !Control::getVec2(label, dim_vec, tpCntl) ) Exit(0);
  
  // ローカル座標系　グリッド分割数　nx、ny
  label = label_base + "/Grid";
  if ( !Control::getVec2(label, div_vec, tpCntl) ) Exit(0);
  
  //matplotlib は、ｍ＝ｎが要求していますので、ここで、最大分割数を使用します。
  int m = (int)div_vec[0];
  int n = (int)div_vec[1];
  int n_div = (std::max)(m, n);
  div_vec[0] = (REAL_TYPE) n_div;
  div_vec[1] = (REAL_TYPE) n_div;
  
  return;
}

// #################################################################
// Graph Ploter
void MonitorList::getPolygon(Monitor_Type mon_type,
                             const string label_base,
                             string & stl_filename )
{
  std::string label=label_base+"/DefFile";
  std::string val_str = "";
  
  if ( !(tpCntl->getInspectedValue(label, val_str )) )
  {
    Hostonly_ stamped_printf("\tParsing error : /DefFile\n");
    Exit(0);
  }
  stl_filename = val_str;
}



// #################################################################
/// VorticityとHelicityのサンプリングが指定されている場合にtrueを返す
bool MonitorList::getStateVorticity()
{
  int count[var_END];
  
  for (int k=var_Velocity; k<var_END; k++)
  {
    count[k]=0;
  }
  
  // モニタで使用する変数を取得
  for (int k=var_Velocity; k<var_END; k++)
  {
    for (int i = 0; i < nGroup; i++)
    {
      if (monGroup[i]->getStateVariable(k) == true) count[k]++;
    }
  }
  
  // SamplingでHelicityあるいはVorticityが指定されている場合には，VorticityをON
  if ( count[var_Helicity]>0 || count[var_Vorticity]>0 ) return true;
  
  return false;
}


// #################################################################
// TPのポインタを受け取る
void MonitorList::importTP(TextParser* tp)
{
  if ( !tp ) Exit(0);
  tpCntl = tp;
}



// #################################################################
/// 出力ファイルオープン
///
///    @note 出力タイプによらず全プロセスから呼んでも問題ない
///
void MonitorList::openFile()
{
  if( nGroup != monGroup.size() )
  {
    printf("\t\t Graph Ploter --- MonitorList::openFile(), ***CAUTION: monGroup.size()=%d, nGroup=%d, we do nGroup = monGroup.size();\n", monGroup.size(), nGroup );
    nGroup = monGroup.size();
  }
  
  if ((outputType == GATHER && myRank == 0) || outputType == DISTRIBUTE)
  {
    for (int i = 0; i < nGroup; i++)
    {
      string fileName(fname_sampling);
      
      string label("_");
      label = label + monGroup[i]->getLabel();
      string::size_type pos = fileName.rfind(".");
      fileName.insert(pos, label);
      
      if (outputType == GATHER)
      {
        monGroup[i]->openFile(fileName.c_str(), true);
      }
      else
      {
        monGroup[i]->openFile(fileName.c_str(), false);
      }
    }
  }
}


// #################################################################
/// モニタ結果出力
///
///   @param [in] step サンプリング時の計算ステップ
///   @param [in] tm   サンプリング時の計算時刻
///
///   @note gatherの場合も全プロセスから呼ぶこと
///
void MonitorList::print(const unsigned step, const double tm)
{
  std::string type_str = getOutputTypeStr().c_str();
  
  if( nGroup != monGroup.size() )
  {
    printf("\t\t Graph ploter --- MonitorList::print(), ***CAUTION: this=%ld, monGroup.size()=%d, nGroup=%d, we do nGroup = monGroup.size();\n",(long)this, monGroup.size(), nGroup );
    nGroup = monGroup.size();
    
  }
  
  for (int i = 0; i < nGroup; i++)
  {
    switch (outputType)
    {
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



// #################################################################
/// モニタ情報を出力
void MonitorList::printMonitorInfo(FILE* fp, const char* str, const bool verbose)
{
  fprintf(fp,"\n----------\n\n");
  fprintf(fp, "\n\t>> Sampling Information\n\n");
  
  fprintf(fp, "\t  Output Type              : %s\n", getOutputTypeStr().c_str());
  fprintf(fp, "\t  Base Name of Output File : %s\n", str);
  fprintf(fp, "\t  Unit                     : %s\n\n", (refVar.modeUnitOutput == DIMENSIONAL) ? "Dimensional" : "Non Dimensional");
  
  if ( verbose )
  {
    for (int i = 0; i < nGroup; i++)
    {
      monGroup[i]->printInfo(fp, i);
    }
  }
  
  fprintf(fp,"\n");
  fflush(fp);
}


// #################################################################
// Graph Ploter
void MonitorList::printTextParser(std::string given_label)
{
  int ierror=0;
  TextParser *tp = tpCntl;
  if( tp != NULL )
  {
    // 全てのパラメータのラベルを取得
    std::vector<std::string> labels;
    ierror = tp->getAllLabels(labels);
    
    for (int i = 0; i < labels.size(); i++)
    {
      std::string label_str = labels[i];
      
      if( given_label.empty() == false )
      {
        if( label_str.find(given_label.c_str()) == string::npos )
        {
          continue;
        }
      }
      
      std::string value_str = "";
      ierror = tp->getValue(labels[i], value_str);
      
      TextParserValueType type = tp->getType(labels[i], &ierror);
      /*
       typedef enum {
       TP_UNDEFFINED_VALUE = 0,  //!< 不定
       TP_NUMERIC_VALUE = 1,     //!< 数値
       TP_STRING_VALUE = 2,      //!< 文字列
       TP_DEPENDENCE_VALUE = 3,  //!< 依存関係付き値
       TP_VECTOR_UNDEFFINED = 4, //!< ベクトル型不定
       TP_VECTOR_NUMERIC = 5,    //!< ベクトル型数値
       TP_VECTOR_STRING = 6,     //!< ベクトル型文字列
       TP_RANGE_NUMERIC = 7,     //!< 値領域指定型 @range
       TP_LIST_NUMERIC = 8,      //!< 値領域指定型 @list
       } TextParserValueType;
       */
      
    }
  }
}




// #################################################################
/**
 * @brief サンプリング対象変数の登録
 * @param [in] label     ノード
 * @param [in] variables モニタ変数vector
 * @param [in] key       変数名
 */
void MonitorList::registVars(const string label, vector<string>& variables, const string key)
{
  string str;
  string leaf = label + "/" + key;
  
  if ( tpCntl->getInspectedValue(leaf, str) )
  {
    // onで登録，なければ登録なし
    if ( !strcasecmp(str.c_str(), "on") )
    {
      variables.push_back(key.c_str());
    }
  }
}


// #################################################################
/// 必要なパラメータのコピー
void MonitorList::setControlVars(int* bid,
                                 int* cut_l,
                                 int* cut_u,
                                 int* bcd,
                                 const REAL_TYPE refVelocity,
                                 const REAL_TYPE baseTemp,
                                 const REAL_TYPE diffTemp,
                                 const REAL_TYPE refDensity,
                                 const REAL_TYPE refLength,
                                 const REAL_TYPE basePrs,
                                 const REAL_TYPE m_RefL,
                                 const int modePrecision,
                                 const int unitPrs,
                                 const int num_process,
                                 const int m_NoCompo,
                                 double* m_mtbl)
{
  this->bid      = bid;
  this->bcd      = bcd;
  this->cutL     = cut_l;
  this->cutU     = cut_u;
  this->g_org.x  = G_origin[0];
  this->g_org.y  = G_origin[1];
  this->g_org.z  = G_origin[2];
  this->g_box.x  = G_region[0];
  this->g_box.y  = G_region[1];
  this->g_box.z  = G_region[2];
  this->org.x    = origin[0];
  this->org.y    = origin[1];
  this->org.z    = origin[2];
  this->pch.x    = pitch[0];
  this->pch.y    = pitch[1];
  this->pch.z    = pitch[2];
  this->box.x    = region[0];
  this->box.y    = region[1];
  this->box.z    = region[2];
  this->num_process = num_process;
  this->NoCompo     = m_NoCompo;
  this->mtbl        = m_mtbl;
  this->RefL     = m_RefL;
  
  refVar.refVelocity   = refVelocity;
  refVar.baseTemp      = baseTemp;
  refVar.diffTemp      = diffTemp;
  refVar.refDensity    = refDensity;
  refVar.refLength     = refLength;
  refVar.basePrs       = basePrs;
  refVar.unitPrs       = unitPrs;
  refVar.modePrecision = modePrecision;
}


// #################################################################
/// Line点登録
///
///   @param [in] str       ラベル文字列
///   @param [in] variables モニタ変数vector
///   @param [in] method    method文字列
///   @param [in] mode      mode文字列
///   @param [in] from      Line始点
///   @param [in] to        Line終点
///   @param [in] nDivision 分割数(モニタ点数-1)
///   @param [in] mon_type  モニタタイプ
///
MonitorCompo* MonitorList::setLine(const char* str,
                          vector<string>& variables,
                          const char* method,
                          const char* mode,
                          REAL_TYPE from[3],
                          REAL_TYPE to[3],
                          int nDivision,
                          Monitor_Type mon_type)
{
  clipLine(from, to);
  
  MonitorCompo* m = new MonitorCompo(org, pch, box, g_org, g_box, refVar, bid, bcd, cutL, cutU, num_process, NoCompo, mtbl);
  
  m->setRankInfo(paraMngr, procGrp);
  m->setDomainInfo(guide, RefL);
  m->setLine(str, variables, method, mode, from, to, nDivision, mon_type);
  
  monGroup.push_back(m);
  
  if( nGroup != nGroup )
  {
    printf("\t\t\tGraph Ploter --- MonitorList::setLine(), ***CAUTION: monGroup.size()=%d, nGroup=%d, we do nGroup = monGroup.size();\n", monGroup.size(), nGroup );
    nGroup = monGroup.size();
  }
  return m;
}


// #################################################################
/// Polygonモニターの交点の数をcmp[]にセット
void MonitorList::setMonitorNpoint(CompoList* cmp, const int NoCompo)
{
  for ( int m=1; m<=NoCompo; m++)
  {
    if ( cmp[m].getType() == MONITOR )
    {
      for (int i = 0; i < nGroup; i++)
      {
        if (monGroup[i]->getPolyID() == m)  cmp[m].setElement(monGroup[i]->getSize());
      }
    }
  }
  
}



// #################################################################
/// PointSet登録
///
///   @param [in] str       ラベル文字列
///   @param [in] variables モニタ変数vector
///   @param [in] method    method文字列
///   @param [in] mode      mode文字列
///   @param [in] pointSet  PointSet
///   @param [in] mon_type  モニタタイプ
///
MonitorCompo* MonitorList::setPointSet(const char* str,
                              vector<string>& variables,
                              const char* method,
                              const char* mode,
                              vector<MonitorCompo::MonitorPoint>& pointSet,
                              Monitor_Type mon_type)
{
  int n_pts = pointSet.size();
  
  if( n_pts <= 0 ) return NULL;
  
  MonitorCompo* m = new MonitorCompo(org, pch, box, g_org, g_box, refVar, bid, bcd, cutL, cutU, num_process, NoCompo, mtbl);
  
  m->setRankInfo(paraMngr, procGrp);
  m->setDomainInfo(guide, RefL);
  m->setPointSet(str, variables, method, mode, pointSet, mon_type);
  
  monGroup.push_back(m);
  return m;
}


// #################################################################
/**
 * @brief Polygon登録
 * @param [in] str       ラベル文字列
 * @param [in] variables モニタ変数vector
 * @param [in] method    method文字列
 * @param [in] mode      mode文字列
 * @param [in] order     エントリ番号
 * @param [in] nv        法線ベクトル
 * @param [in] mon_type  モニタタイプ
 * @note orderはLocalBCのエントリ
 */
void MonitorList::setPolygon(const char* str,
                             vector<string>& variables,
                             const char* method,
                             const char* mode,
                             const int order,
                             const REAL_TYPE nv[3],
                             Monitor_Type mon_type)
{
  MonitorCompo* m = new MonitorCompo(org, pch, box, g_org, g_box, refVar, bid, bcd, cutL, cutU, num_process, NoCompo, mtbl);
  
  m->setRankInfo(paraMngr, procGrp);
  m->setDomainInfo(guide, RefL);
  m->setPolygon(str, variables, method, mode, order, nv, mon_type);
  
  monGroup.push_back(m);
}


// #################################################################
/**
 * @brief Primitiveの登録
 * @param [in]  str       ラベル文字列
 * @param [in]  variables モニタ変数vector
 * @param [in]  method    method文字列
 * @param [in]  mode      mode文字列
 * @param [in]  mon_type  モニタタイプ
 * @param [in]  nv        法線ベクトル
 * @param [in]  center    中心座標
 * @param [in]  depth     厚さ
 * @param [in]  tmp       幅 / 半径
 * @param [in]  height    高さ
 * @param [in]  dir       方向ベクトル
 * @param [in]  order     bcdにペイントする識別ID（パラメータの出現番号を利用）
 */
void MonitorList::setPrimitive(const char* str,
                               vector<string>& variables,
                               const char* method,
                               const char* mode,
                               const int order,
                               const REAL_TYPE nv[3],
                               Monitor_Type mon_type)
{
  MonitorCompo* m = new MonitorCompo(org, pch, box, g_org, g_box, refVar, bid, bcd, cutL, cutU, num_process, NoCompo, mtbl);

  m->setRankInfo(paraMngr, procGrp);
  m->setDomainInfo(guide, RefL);
  m->setPrimitive(str, variables, method, mode, order, nv, mon_type);
  
  monGroup.push_back(m);
}
