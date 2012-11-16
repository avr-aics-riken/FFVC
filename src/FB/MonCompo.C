// #################################################################
//
// FFV : Frontflow / violet
//
// Copyright (c) 2012 All right reserved.
//
// Institute of Industrial Science, University of Tokyo, Japan. 
//
// #################################################################

/** 
 * @file   MonCompo.h
 * @brief  FlowBase MonitorCompo class
 * @author kero
 */

#include "MonCompo.h"
#include <sstream>


/// モニタリング管理用配列の確保.
void MonitorCompo::allocArray()
{
  if (nPoint == 0) Exit(0); // サンプリング点数
  
  if (!(crd  = new FB::Vec3r[nPoint]))    Exit(0);
  if (!(rank = new int[nPoint]))      Exit(0);
  if (!(comment = new string[nPoint])) Exit(0);
  if (!(pointStatus = new int[nPoint])) Exit(0);
  
  if (!(mon = new Sampling*[nPoint])) Exit(0);
  for (int i = 0; i < nPoint; i++) mon[i] = NULL;
}


/// サンプリング値を格納する配列の確保.
void MonitorCompo::allocSamplingArray()
{
  if ( nPoint == 0 ) Exit(0); // サンプリング点数
  
  const REAL_TYPE DUMMY = 1.0e10;
  //const REAL_TYPE DUMMY = 0.0;
  
  if (variable[VELOCITY]) {
    if (!(vel = new FB::Vec3r[nPoint])) Exit(0);
    for (int i = 0; i < nPoint; i++) vel[i] = FB::Vec3r(DUMMY, DUMMY, DUMMY);
  }
  if (variable[PRESSURE]) {
    if (!(prs = new REAL_TYPE[nPoint])) Exit(0);
    for (int i = 0; i < nPoint; i++) prs[i] = DUMMY;
  }
  if (variable[TEMPERATURE]) {
    if (!(tmp = new REAL_TYPE[nPoint])) Exit(0);
    for (int i = 0; i < nPoint; i++) tmp[i] = DUMMY;
  }
  if (variable[TOTAL_PRESSURE]) {
    if (!(tp = new REAL_TYPE[nPoint])) Exit(0);
    for (int i = 0; i < nPoint; i++) tp[i] = DUMMY;
  }
  //if (variable[VORTICITY]) {
  //  if (!(vor = new FB::Vec3r[nPoint])) Exit(0);
  //  for (int i = 0; i < nPoint; i++) vor[i] = FB::Vec3r(DUMMY, DUMMY, DUMMY);
  //}
}


/// Allreduceによる総和(実数配列上書き，work配列指定).
bool MonitorCompo::allReduceSum(REAL_TYPE* array, int n, REAL_TYPE* sendBuf)
{
  if ( numProc <= 1 ) return true;
  
  for (int i = 0; i < n; i++) sendBuf[i] = array[i];
  

  if( sizeof(REAL_TYPE) == 8 )
  {
    if( MPI_Allreduce(sendBuf, array, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD) != MPI_SUCCESS ) return false;
  }
  else
  {
    if( MPI_Allreduce(sendBuf, array, n, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD) != MPI_SUCCESS ) return false;
  }
  
  return true;
}


/// Allreduceによる総和(実数配列上書き).
bool MonitorCompo::allReduceSum(REAL_TYPE* array, int n)
{
  if ( numProc <= 1 ) return true;
  
  REAL_TYPE* sendBuf = new REAL_TYPE[n];
  bool ret = allReduceSum(array, n, sendBuf);
  delete[] sendBuf;
  
  return ret;
}


/// Allreduceによる総和(整数配列上書き，work配列指定).
bool MonitorCompo::allReduceSum(int* array, int n, unsigned long* sendBuf)
{
  if ( numProc <= 1 ) return true;
  
  //for (int i = 0; i < n; i++) sendBuf[i] = array[i];
  //if ( MPI_Allreduce(sendBuf, array, n, MPI_INT, MPI_SUM, MPI_COMM_WORLD) != MPI_SUCCESS ) return false;
  
  // 一度ulongへキャストし、チェック
  unsigned long* recvBuf = new unsigned long[n];
  for (int i = 0; i < n; i++) sendBuf[i] = (unsigned long)array[i];
  if ( MPI_Allreduce(sendBuf, recvBuf, n, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD) != MPI_SUCCESS ) return false;
  for (int i = 0; i < n; i++)
  {
    if (recvBuf[i] > INT_MAX)
    {
      Hostonly_ stamped_printf("\tError : allReduceSum OverFlow\n");
      Exit(0);
    }
    array[i] = (int)recvBuf[i];
  }
  
  return true;
}


/// Allreduceによる総和(整数配列上書き).
bool MonitorCompo::allReduceSum(int* array, int n)
{
  if ( numProc <= 1 ) return true;
  
  unsigned long* sendBuf = new unsigned long[n];
  bool ret = allReduceSum(array, n, sendBuf);
  delete[] sendBuf;
  
  return ret;
}


/// 内部境界条件として指定されたモニタ領域内でスカラー変数を平均.
///
///   @param[in] s スカラー変数配列
///   @return モニタ領域内平均値
///
REAL_TYPE MonitorCompo::averageScalar(REAL_TYPE* s)
{
  REAL_TYPE sum= 0.0;
  for (int i = 0; i < nPoint; i++) {
    if (rank[i] == myRank) sum+= s[i];
  }
  allReduceSum(&sum, 1);
  
  return sum / nPoint;
}


/// 内部境界条件として指定されたモニタ領域内でベクトル変数を平均.
///
///   @param[in] v ベクトル変数配列
///   @return モニタ領域内平均値
///
FB::Vec3r MonitorCompo::averageVector(FB::Vec3r* v)
{
  REAL_TYPE sum[3] = { 0.0, 0.0, 0.0 };
  for (int i = 0; i < nPoint; i++) {
    if (rank[i] == myRank) {
      sum[0] += v[i].x;
      sum[1] += v[i].y;
      sum[2] += v[i].z;
    }
  }
  allReduceSum(sum, 3);
  
  return FB::Vec3r(sum) / nPoint;
}



/// モニタ点の状態を調べ，不正モニタ点フラグ配列pointStatusを設定.
void MonitorCompo::checkMonitorPoints()
{
  for (int m = 0; m < nPoint; m++) {
    pointStatus[m] = 0;
    if (rank[m] == myRank) {
      pointStatus[m] = mon[m]->checkMonitorPoint();
    }
  }
  
  if (!allReduceSum(pointStatus, nPoint)) Exit(0);
  
}



/// モニタ点が指定領域内にあるかを判定.
///
///   @param[in] m モニタ点番号
///   @param[in] org 調査領域の基点
///   @param[in] box 調査領域のサイズ
///   @param[in] flag メッセージ出力フラグ(trueの時出力)
///   @return true=領域内/false=領域外
///
bool MonitorCompo::check_region(int m, FB::Vec3r org, FB::Vec3r box, bool flag) 
{
  if ( (crd[m].x< org.x)        ||
      (crd[m].x>(org.x+box.x))  ||
      (crd[m].y< org.y)         ||
      (crd[m].y>(org.y+box.y))  ||
      (crd[m].z< org.z)         ||
      (crd[m].z>(org.z+box.z)) ) { 
        if (flag) { stamped_printf("\trank=%d : [%14.6e %14.6e %14.6e] is out of region\n", myRank, crd[m].x, crd[m].y, crd[m].z); }
    return false;
  }
  return true;
}


/// 出力ファイルクローズ.
void MonitorCompo::closeFile()
{ 
  assert(fp);
  fclose(fp);
}




/// サンプリングした変数をノード0に集約.
void MonitorCompo::gatherSampled()
{
  REAL_TYPE* vSendBuf = NULL;
  REAL_TYPE* vRecvBuf = NULL;
  REAL_TYPE* sRecvBuf = NULL;
  
  if (variable[VELOCITY])       gatherSampledVector(vel, vSendBuf, vRecvBuf);
  if (variable[PRESSURE])       gatherSampledScalar(prs, sRecvBuf);
  if (variable[TEMPERATURE])    gatherSampledScalar(tmp, sRecvBuf);
  if (variable[TOTAL_PRESSURE]) gatherSampledScalar(tp, sRecvBuf);
  //if (variable[VORTICITY])      gatherSampledVector(vor, vSendBuf, vRecvBuf);
  
  if (vSendBuf) delete[] vSendBuf;
  if (vRecvBuf) delete[] vRecvBuf;
  if (sRecvBuf) delete[] sRecvBuf;
}


/// サンプリングしたスカラー変数をノード0に集約.
///
///   @param [in,out] s         スカラー変数配列
///   @param          sRecvBuf  通信用work領域
///
void MonitorCompo::gatherSampledScalar(REAL_TYPE* s, REAL_TYPE* sRecvBuf)
{
  int np = num_process;
  if ( numProc <= 1 ) return;
  
  if (myRank == 0 && !sRecvBuf)
  {
    if (!(sRecvBuf = new REAL_TYPE[nPoint*np])) Exit(0);
  }
  
  if ( numProc > 1 ) 
  {
    if ( sizeof(REAL_TYPE) == 8 )
    {
      if( MPI_Gather(s, nPoint, MPI_DOUBLE, sRecvBuf, nPoint, MPI_DOUBLE, 0, MPI_COMM_WORLD) != MPI_SUCCESS ) Exit(0);
    }
    else
    {
      if( MPI_Gather(s, nPoint, MPI_FLOAT, sRecvBuf, nPoint, MPI_FLOAT, 0, MPI_COMM_WORLD) != MPI_SUCCESS ) Exit(0);
    }
  }
  
  if (myRank == 0) 
  {
    for (int i = 0; i < np; i++) {
      for (int m = 0; m < nPoint; m++) {
        if (rank[m] == i) 
        {
          s[m] = sRecvBuf[nPoint*i+m];
        }
      }
    }
  }
}


/// サンプリングしたベクトル変数をノード0に集約.
///
///   @param[in,out] v ベクトル変数配列
///   @param  vSendBuf,vRecvBuf  通信用work領域
///
void MonitorCompo::gatherSampledVector(FB::Vec3r* v, REAL_TYPE* vSendBuf, REAL_TYPE* vRecvBuf)
{
  int np = num_process;
  if ( numProc <= 1 ) return;
  
  if (!vSendBuf) 
  {
    if (!(vSendBuf = new REAL_TYPE[nPoint*3*np])) Exit(0);
  }
  
  if (myRank == 0 && !vRecvBuf) 
  {
    if (!(vRecvBuf = new REAL_TYPE[nPoint*3*np])) Exit(0);
  }
  
  for (int m = 0; m < nPoint; m++) {
    vSendBuf[3*m  ] = v[m].x;
    vSendBuf[3*m+1] = v[m].y;
    vSendBuf[3*m+2] = v[m].z;
  }
  
  if ( numProc > 1 )
  {
    if ( sizeof(REAL_TYPE) == 8 )
    {
      if( MPI_Gather(vSendBuf, nPoint*3, MPI_DOUBLE, vRecvBuf, nPoint*3, MPI_DOUBLE, 0, MPI_COMM_WORLD) != MPI_SUCCESS ) Exit(0);
    }
    else
    {
      if( MPI_Gather(vSendBuf, nPoint*3, MPI_FLOAT, vRecvBuf, nPoint*3, MPI_FLOAT, 0, MPI_COMM_WORLD) != MPI_SUCCESS ) Exit(0);
    }
  }
  
  if (myRank == 0) 
  {
    for (int i = 0; i < np; i++) {
      for (int m = 0; m < nPoint; m++) {
        if (rank[m] == i) {
          v[m].x = vRecvBuf[nPoint*3*i+m*3  ];
          v[m].y = vRecvBuf[nPoint*3*i+m*3+1];
          v[m].z = vRecvBuf[nPoint*3*i+m*3+2];
        }
      }
    }
  }
}



/// モニタ点指定方法文字列の取得.
string MonitorCompo::getTypeStr()
{
  string str;
  
  switch (type) {
    case POINT_SET:
      str = "PointSet"; break;
    case LINE:
      str = "Line"; break;
    case INNER_BOUNDARY:
      str = "LocalBoundary"; break;
    default:
      Exit(0);
  }
  
  return str;
}


/// サンプリング方法文字列の取得.
string MonitorCompo::getMethodStr()
{
  string str;
  
  switch (method) {
    case SAMPLING_NEAREST:
      str = "Nearest"; break;
    case SAMPLING_INTERPOLATION:
      str = "Interpolation"; break;
    case SAMPLING_SMOOTHING:
      str = "Smoothing"; break;
    default:
      Exit(0);
  }
  
  return str;
}


/// サンプリングモード文字列の取得.
string MonitorCompo::getModeStr()
{
  string str;
  
  switch (mode) {
    case SAMPLING_ALL:
      str = "All"; break;
    case SAMPLING_FLUID_ONLY:
      str = "Fluid"; break;
    case SAMPLING_SOLID_ONLY:
      str = "Solid"; break;
    default:
      Exit(0);
  }
  
  return str;
}


/// モニタ変数を結合した文字列の取得.
string MonitorCompo::getVarStr() 
{
  string var;
  
  if (variable[VELOCITY])        var += "Velocity ";
  if (variable[PRESSURE])        var += "Pressure ";
  if (variable[TEMPERATURE])     var += "Temperature ";
  if (variable[TOTAL_PRESSURE])  var += "TotalPressure ";
  //if (variable[VORTICITY])       var += "Vorticity ";
  
  return var;
}



/// 出力ファイルオープン.
///
///    @param str ファイル名テンプレート
///    @param gathered true=gather出力/false=disutribute出力
///
void MonitorCompo::openFile(const char* str, bool gathered)
{
  if (gathered) 
  {
    if (!(fp = fopen(str, "w"))) 
    {
      perror(str); Exit(0);
    }
    writeHeader(true);
  }
  else 
  {
    string fileName(str);
    ostringstream rankStr;
    rankStr << "_" << myRank;
    string::size_type pos = fileName.rfind(".");
    fileName.insert(pos, rankStr.str());
    if (!(fp = fopen(fileName.c_str(), "w"))) 
    {
      perror(fileName.c_str()); Exit(0);
    }
    writeHeader(false);
  }
}


/// モニタ結果出力.
///
///   @param[in] step サンプリング時の計算ステップ
///   @param[in] tm   サンプリング時の計算時刻
///   @param[in] gathered 出力モードフラグ(true=gather出力/false=disutribute出力)
///
void MonitorCompo::print(unsigned step, REAL_TYPE tm, bool gathered)
{
  assert(fp);
  
  char* sFmtSingle = "%15.7e ";
  char* vFmtSingle = "%15.7e %15.7e %15.7e ";
  char* sFmtDouble = "%24.16e ";
  char* vFmtDouble = "%24.16e %24.16e %24.16e ";
  char* sFmt;
  char* vFmt;
  
  if (refVar.modePrecision == sizeof(float))
  {
    sFmt = sFmtSingle;
    vFmt = vFmtSingle;
  }
  else 
  {
    sFmt = sFmtDouble;
    vFmt = vFmtDouble;
  }
  
  fprintf(fp, "\n");
  if (refVar.modePrecision == sizeof(float)) 
  {
    fprintf(fp, "%d %14.6e\n", step, convTime(tm)); // %12.4 >> %14.6
  }
  else 
  {
    fprintf(fp, "%d %24.16e\n", step, convTime(tm));
  }
  
  for (int i = 0; i < nPoint; i++) {
    if (!gathered && rank[i] != myRank) continue;
    if (pointStatus[i] != Sampling::POINT_STATUS_OK) 
    {
      fprintf(fp, "%s\n", "  *NA*");
      continue;
    }
    if (variable[VELOCITY]) 
    {
      fprintf(fp, vFmt, convVel(vel[i].x), convVel(vel[i].y), convVel(vel[i].z));
    }
    if (variable[PRESSURE]) 
    {
      fprintf(fp, sFmt, convPrs(prs[i]));
    }
    if (variable[TEMPERATURE]) 
    {
      fprintf(fp, sFmt, convTmp(tmp[i]));
    }
    if (variable[TOTAL_PRESSURE]) 
    {
      fprintf(fp, sFmt, convTP(tp[i]));
    }
    //if (variable[VORTICITY]) {
    //  fprintf(fp, vFmt, convVor(vor[i].x), convVor(vor[i].y), convVor(vor[i].z));
    //}
    fprintf(fp, "\n");
  }
}


/// モニタ情報を出力.
///
///    @param[in] fp 出力ファイルポインタ
///    @param[in] no モニタグループ通し番号
//
void MonitorCompo::printInfo(FILE* fp, int no)
{
  fprintf(fp,"\t%3d : %s\t division=%d  [%s]\n", no+1, getTypeStr().c_str(), nPoint, label.c_str());
  fprintf(fp,"\t\tVariables : %s\n", getVarStr().c_str());
  fprintf(fp,"\t\t   Method : %s\n", getMethodStr().c_str());
  fprintf(fp,"\t\t     Mode : %s\n", getModeStr().c_str());
  
  fprintf(fp, "\t\t    order :            X              Y              Z    :   rank : comment\n");
  for (int j = 0; j < nPoint; j++) {
    fprintf(fp,"\t\t%9d : %14.6e %14.6e %14.6e  : %6d : %s", // %12.4 >> %14.6
            j+1, convCrd(crd[j].x), convCrd(crd[j].y), convCrd(crd[j].z), rank[j], comment[j].c_str());
    if (pointStatus[j] == Sampling::UNEXPECTED_SOLID) 
    {
      fprintf(fp, "  *skip(unexpected solid)*\n");
    }
    else if (pointStatus[j] == Sampling::UNEXPECTED_FLUID) 
    {
      fprintf(fp, "  *skip(unexpected fluid)*\n");
    }
    else 
    {
      fprintf(fp, "\n");
    }
  }
  fprintf(fp,"\n");
}


/// サンプリング(Line, PointSet).
void MonitorCompo::sampling()
{
  for (int i = 0; i < nPoint; i++) {
    //  if (!(mon[i] && pointStatus[i] == Sampling::POINT_STATUS_OK)) continue;
    if (!mon[i]) continue;
    if (variable[VELOCITY])       vel[i] = mon[i]->samplingVelocity(vSource);
    if (variable[PRESSURE])       prs[i] = mon[i]->samplingPressure(pSource);
    if (variable[TEMPERATURE])    tmp[i] = mon[i]->samplingTemperature(tSource);
    if (variable[TOTAL_PRESSURE]) tp[i]  = mon[i]->samplingTotalPressure(vSource, pSource);
    //if (variable[VORTICITY])      vor[i] = mon[i]->samplingVorticity(vSource);
  }
}



/// 内部境界条件モニタ点でのサンプリング.
///
///   サンプリング結果を集計，コンポーネント領域での平均値を計算.
///   速度は法線ベクトルとの内積をとる．
///   結果はコンポーネントcmpに格納.
///
void MonitorCompo::samplingInnerBoundary()
{
  assert(type == INNER_BOUNDARY);
  sampling();
  
  if (variable[VELOCITY]) 
  {
    FB::Vec3r velAve = averageVector(vel);
    cmp->val[var_Velocity]
    = velAve.x * cmp->nv[0] + velAve.y * cmp->nv[1] + velAve.z * cmp->nv[2];
  }
  if (variable[PRESSURE])       cmp->val[var_Pressure]    = averageScalar(prs);
  if (variable[TEMPERATURE])    cmp->val[var_Temperature] = averageScalar(tmp);
  if (variable[TOTAL_PRESSURE]) cmp->val[var_TotalP]      = averageScalar(tp);
}



/// 内部境界条件として指定されたモニタ領域のセル中心座標をcrd[]に設定.
///
///   @param[in] n コンポーネントエントリ
///   @param[in] cmp コンポーネント
///
void MonitorCompo::setIBPoints(int n, CompoList& cmp)
{
  int np = num_process;
  int st[3], ed[3];
  int s;
  
  int* nPointList;
  if (!(nPointList = new int[np])) Exit(0);
  for (int i = 0; i < np; i++) nPointList[i] = 0;
  
  cmp.getBbox(st, ed);

  // for optimization > variables defined outside
  size_t mm;
  int ix, jx, kx, gd;
  ix = size[0];
  jx = size[1];
  kx = size[2];
  gd = guide;
  
  if (cmp.isEns()) 
  {
    for (int k = st[2]; k <= ed[2]; k++) {
      for (int j = st[1]; j <= ed[1]; j++) {
        for (int i = st[0]; i <= ed[0]; i++) {
          mm = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
          s = bcd[mm];
          
          if ((s & MASK_6) == n) 
          {
            nPointList[myRank]++;
          }
        }
      }
    }
  }
  
  if (!allReduceSum(nPointList, np)) Exit(0);
  
  
  //check
  int sum = 0;
  for (int i = 0; i < np; i++) sum += nPointList[i];
  assert(sum == nPoint);
  
  REAL_TYPE* buf;
  if (!(buf = new REAL_TYPE[nPoint*3])) Exit(0);
  for (int i = 0; i < nPoint*3; i++) buf[i] = 0.0;
  
  int m0 = 0;
  for (int i = 0; i < myRank; i++) m0 += nPointList[i];
  
  int m = 0;
  if (cmp.isEns()) 
  {
    int i0 = head[0] - 1; //fortran index なので1からカウントアップしている ---> Cのindexは0からカウントアップ
    int j0 = head[1] - 1;
    int k0 = head[2] - 1;

    for (int k = st[2]; k <= ed[2]; k++) {
      for (int j = st[1]; j <= ed[1]; j++) {
        for (int i = st[0]; i <= ed[0]; i++) {
          mm=_F_IDX_S3D(i, j, k, ix, jx, kx, gd);
          if ((bcd[mm] & MASK_6) == n) 
          {
            buf[(m0+m)*3+0] = g_org.x + (i + i0 - 0.5) * pch.x;
            buf[(m0+m)*3+1] = g_org.y + (j + j0 - 0.5) * pch.y;
            buf[(m0+m)*3+2] = g_org.z + (k + k0 - 0.5) * pch.z;
            m++;
          }
        }
      }
    }
  }
  
  if (!allReduceSum(buf, nPoint*3)) Exit(0);
  
  for (m = 0; m < nPoint; m++) {
    crd[m].x = buf[m*3+0];
    crd[m].y = buf[m*3+1];
    crd[m].z = buf[m*3+2];
  }
  
  delete[] buf;
  delete[] nPointList;
}


/// 内部境界条件としてモニタ点を登録.
///
///   @param[in] n コンポーネントエントリ
///   @param[in] cmp コンポーネント
///
void MonitorCompo::setInnerBoundary(int n, CompoList& cmp)
{
  ostringstream oss;
  oss << "LocalBoundary" << n;
  type = INNER_BOUNDARY;
  label = oss.str();
  method = SAMPLING_NEAREST;
  mode = SAMPLING_ALL;
  
  this->cmp = &cmp;
  
  if (cmp.isVarEncoded(var_Velocity))    variable[VELOCITY]       = true;
  if (cmp.isVarEncoded(var_Pressure))    variable[PRESSURE]       = true;
  if (cmp.isVarEncoded(var_Temperature)) variable[TEMPERATURE]    = true;
  if (cmp.isVarEncoded(var_TotalP))      variable[TOTAL_PRESSURE] = true;
  
  
  nPoint = (int)cmp.getElement();
  allocArray();
  allocSamplingArray();
  
  setIBPoints(n, cmp);
  
  setRankArray();
  
  for (int m = 0; m < nPoint; m++) {
    ostringstream oss;
    oss << "point_" << m;
    comment[m] = oss.str();
    if (!check_region(m, g_org, g_box, true)) Exit(0);
    if (rank[m] == myRank) 
    {
      mon[m] = new Nearest(mode, size, guide, crd[m], org, pch, refVar.v00, bcd);
    }
  }
  
  checkMonitorPoints();
}



/// Line登録.
///
///   @param[in] labelStr ラベル文字列
///   @param[in] variables モニタ変数vector
///   @param[in] methodStr method文字列
///   @param[in] modeStr   mode文字列
///   @param[in] from Line始点
///   @param[in] to   Line終点
///   @param[in] nDivision 分割数(モニタ点数-1)
///
void MonitorCompo::setLine(const char* labelStr, vector<string>& variables,
                           const char* methodStr, const char* modeStr,
                           REAL_TYPE from[3], REAL_TYPE to[3], int nDivision)
{
  type = LINE;
  label = labelStr;
  setSamplingMethod(methodStr);
  setSamplingMode(modeStr);
  
  vector<string>::const_iterator it;
  for (it = variables.begin(); it != variables.end(); it++) {
    setMonitorVar((*it).c_str());
  }
  
  nPoint = nDivision + 1;
  if (nPoint < 2) Exit(0);
  
  allocArray();
  allocSamplingArray();
  
  FB::Vec3r st(from), ed(to);
  FB::Vec3r dd = ed - st;
  dd /= nPoint - 1;
  
  for (int m = 0; m < nPoint; m++) {
    ostringstream oss;
    oss << "point_" << m;
    crd[m] = st + dd * (REAL_TYPE)m;
    if (!check_region(m, g_org, g_box, true)) Exit(0);
    comment[m] = oss.str();
  }
  
  setRankArray();
  
  for (int m = 0; m < nPoint; m++) {
    if (rank[m] == myRank) 
    {
      switch (method) {
        case SAMPLING_NEAREST:
          mon[m] = new Nearest(mode, size, guide, crd[m], org, pch, refVar.v00, bcd);
          break;
          
        case SAMPLING_INTERPOLATION:
          mon[m] = new Interpolation(mode, size, guide, crd[m], org, pch, refVar.v00, bcd);
          break;
          
        case SAMPLING_SMOOTHING:
          mon[m] = new Smoothing(mode, size, guide, crd[m], org, pch, refVar.v00, bcd);
          break;
          
        default:
          Exit(0);
      }
    }
  }
  
  checkMonitorPoints();
}



/// モニタ対象物理量の設定.
///
///    @param[in] str モニタ対象物理量文字列
///
void MonitorCompo::setMonitorVar(const char* str) 
{
  if      (!strcasecmp("velocity", str))       variable[VELOCITY] = true;
  else if (!strcasecmp("pressure", str))       variable[PRESSURE] = true;
  else if (!strcasecmp("temperature", str))    variable[TEMPERATURE] = true;
  else if (!strcasecmp("totalpressure", str))  variable[TOTAL_PRESSURE] = true;
  //else if (!strcasecmp("vorticity", str))       variable[VORTICITY] = true;
  else {
    Hostonly_ stamped_printf("\tError : Invalid variable keyword [%s]\n", str);
    Exit(0);
  }
}


/// PointSet登録.
///
///   @param[in] labelStr ラベル文字列
///   @param[in] variables モニタ変数vector
///   @param[in] methodStr method文字列
///   @param[in] modeStr   mode文字列
///   @param[in] pointSet  PointSet
///
void MonitorCompo::setPointSet(const char* labelStr, vector<string>& variables,
                               const char* methodStr, const char* modeStr,
                               vector<MonitorPoint>& pointSet)
{
  type = POINT_SET;
  label = labelStr;
  setSamplingMethod(methodStr);
  setSamplingMode(modeStr);
  
  vector<string>::const_iterator it;
  for (it = variables.begin(); it != variables.end(); it++) {
    setMonitorVar((*it).c_str());
  }
  
  nPoint = pointSet.size();
  if (nPoint == 0) Exit(0); // サンプリング点数
  
  allocArray();
  allocSamplingArray();
  
  for (int m = 0; m < nPoint; m++) {
    crd[m] = pointSet[m].crd;
    if (!check_region(m, g_org, g_box, true)) Exit(0);
    comment[m] = pointSet[m].label;
  }
  
  setRankArray();
  
  for (int m = 0; m < nPoint; m++) {
    if (rank[m] == myRank) 
    {
      switch (method) {
        case SAMPLING_NEAREST:
          mon[m] = new Nearest(mode, size, guide, crd[m], org, pch, refVar.v00, bcd);
          break;
          
        case SAMPLING_INTERPOLATION:
          mon[m] = new Interpolation(mode, size, guide, crd[m], org, pch, refVar.v00, bcd);
          break;
          
        case SAMPLING_SMOOTHING:
          mon[m] = new Smoothing(mode, size, guide, crd[m], org, pch, refVar.v00, bcd);
          break;
          
        default:
          Exit(0);
      }
    }
  }
  
  checkMonitorPoints();
}



/// 各モニタ点を担当するランク番号を配列rank[]にセット.
///
///   @note 領域境界上のモニタ点は，ランク番号の大きい方の領域が担当
///
void MonitorCompo::setRankArray()
{
  int* sendBuf =  new int[nPoint];
  if (!sendBuf) Exit(0);
  
  for (int m = 0; m < nPoint; m++) {
    sendBuf[m] = -1;
    if (check_region(m, org, box)) sendBuf[m] = myRank;
  }
  
  // gather max rank number
  if ( numProc > 1 ) 
  {
    if ( MPI_Allreduce(sendBuf, rank, nPoint, MPI_INT, MPI_MAX, MPI_COMM_WORLD) != MPI_SUCCESS ) Exit(0);
  }
  else 
  {
    for (int m = 0; m < nPoint; m++) rank[m] = sendBuf[m];
  }
  
  delete[] sendBuf;
}




/// サンプリング方法の設定.
///
///    @param[in] str サンプリング方法文字列
///
void MonitorCompo::setSamplingMethod(const char* str) 
{
  if      (!strcasecmp("nearest", str))       method = SAMPLING_NEAREST;
  else if (!strcasecmp("interpolation", str)) method = SAMPLING_INTERPOLATION;
  else if (!strcasecmp("smoothing", str))     method = SAMPLING_SMOOTHING;
  else {
    Hostonly_ stamped_printf("\tError : Invalid samping_method keyword [%s]\n", str);
    Exit(0);
  }
}


/// サンプリングモードの設定.
///
///    @param[in] str サンプリングモード文字列
///
void MonitorCompo::setSamplingMode(const char* str) 
{
  if      (!strcasecmp("all",   str)) mode = SAMPLING_ALL;
  else if (!strcasecmp("fluid", str)) mode = SAMPLING_FLUID_ONLY;
  else if (!strcasecmp("solid", str)) mode = SAMPLING_SOLID_ONLY;
  else {
    Hostonly_ stamped_printf("\tError : Invalid samping_mode keyword [%s]\n", str);
    Exit(0);
  }
}


/// モニタ結果出力ファイルにヘッダ部を出力.
///
///   @param[in] gathered 出力モードフラグ(true=gather出力/false=disutribute出力)
///
void MonitorCompo::writeHeader(bool gathered)
{
  assert(fp);
  
  int n;
  if (gathered) {
    n = nPoint;
  }
  else {
    n = 0;
    for (int i = 0; i < nPoint; i++) if (rank[i] == myRank) n++;
  }
  
  fprintf(fp, "%d %s\n", n, getVarStr().c_str());
  
  for (int i = 0; i < nPoint; i++) {
    if (gathered || rank[i] == myRank) 
    {
      fprintf(fp, "%14.6e %14.6e %14.6e  #%s", // %12.4 >> %14.6
              convCrd(crd[i].x), convCrd(crd[i].y), convCrd(crd[i].z), comment[i].c_str());
      
      if (pointStatus[i] == Sampling::UNEXPECTED_SOLID) 
      {
        fprintf(fp, "  *skip(unexpected solid)*\n");
      }
      else if (pointStatus[i] == Sampling::UNEXPECTED_FLUID) 
      {
        fprintf(fp, "  *skip(unexpected fluid)*\n");
      }
      else 
      {
        fprintf(fp, "\n");
      }
    }
  }
  
  fflush(fp);
}
