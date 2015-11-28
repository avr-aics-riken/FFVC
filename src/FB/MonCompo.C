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
// Copyright (c) 2012-2015 Advanced Institute for Computational Science, RIKEN.
// All rights reserved.
//
//##################################################################################

/** 
 * @file   MonCompo.h
 * @brief  FlowBase MonitorCompo class
 * @author aics
 */

#include "MonCompo.h"
#include <sstream>


// #################################################################
/// モニタリング管理用配列の確保
void MonitorCompo::allocArray()
{
  if (nPoint == 0) Exit(0); // サンプリング点数
  
  if (!(crd  = new Vec3r[nPoint]))    Exit(0);
  if (!(rank = new int[nPoint]))      Exit(0);
  if (!(comment = new string[nPoint])) Exit(0);
  if (!(pointStatus = new int[nPoint])) Exit(0);
  
  if (!(mon = new Sampling*[nPoint])) Exit(0);
  for (int i = 0; i < nPoint; i++) mon[i] = NULL;
}


// #################################################################
/// サンプリング値を格納する配列の確保
void MonitorCompo::allocSamplingArray()
{
  if ( nPoint == 0 ) Exit(0); // サンプリング点数
  
  const REAL_TYPE DUMMY = 1.0e10;

  
  if (variable[var_Velocity])
  {
    if (!(vel = new Vec3r[nPoint])) Exit(0);
    for (int i = 0; i < nPoint; i++) vel[i] = Vec3r(DUMMY, DUMMY, DUMMY);
  }
  
  if (variable[var_Pressure])
  {
    if (!(prs = new REAL_TYPE[nPoint])) Exit(0);
    for (int i = 0; i < nPoint; i++) prs[i] = DUMMY;
  }
  
  if (variable[var_Temperature])
  {
    if (!(tmp = new REAL_TYPE[nPoint])) Exit(0);
    for (int i = 0; i < nPoint; i++) tmp[i] = DUMMY;
  }

  
  
  if (variable[var_TotalP])
  {
    if (!(tp = new REAL_TYPE[nPoint])) Exit(0);
    for (int i = 0; i < nPoint; i++) tp[i] = DUMMY;
  }
  
  if (variable[var_Helicity])
  {
    if (!(hlt = new REAL_TYPE[nPoint])) Exit(0);
    for (int i = 0; i < nPoint; i++) hlt[i] = DUMMY;
  }
  
  if (variable[var_Vorticity])
  {
    if (!(vor = new Vec3r[nPoint])) Exit(0);
    for (int i = 0; i < nPoint; i++) vor[i] = Vec3r(DUMMY, DUMMY, DUMMY);
  }
}


// #################################################################
/// Allreduceによる総和(実数配列上書き，work配列指定)
bool MonitorCompo::allReduceSum(REAL_TYPE* array, int n, REAL_TYPE* sendBuf)
{
  if ( numProc <= 1 ) return true;
  
  for (int i = 0; i < n; i++) sendBuf[i] = array[i];
  
  if( sizeof(REAL_TYPE) == 8 )
  {
    if( MPI_Allreduce(sendBuf, array, n, MPI_DOUBLE, MPI_SUM, paraMngr->GetMPI_Comm(procGrp)) != MPI_SUCCESS ) return false;
  }
  else
  {
    if( MPI_Allreduce(sendBuf, array, n, MPI_FLOAT, MPI_SUM, paraMngr->GetMPI_Comm(procGrp)) != MPI_SUCCESS ) return false;
  }
  
  return true;
}


// #################################################################
/// Allreduceによる総和(実数配列上書き)
bool MonitorCompo::allReduceSum(REAL_TYPE* array, int n)
{
  if ( numProc <= 1 ) return true;
  
  REAL_TYPE* sBuf=NULL;
  if ( !(sBuf = new REAL_TYPE[n])) Exit(0);
  
  bool ret = allReduceSum(array, n, sBuf);

  if ( sBuf ) delete[] sBuf; sBuf=NULL;
  
  return ret;
}


// #################################################################
/// Allreduceによる総和(整数配列上書き，work配列指定)
bool MonitorCompo::allReduceSum(int* array, int n, unsigned long* sendBuf)
{
  if ( numProc <= 1 ) return true;
  
  //for (int i = 0; i < n; i++) sendBuf[i] = array[i];
  //if ( MPI_Allreduce(sendBuf, array, n, MPI_INT, MPI_SUM, GetMPI_Comm(procGrp)) != MPI_SUCCESS ) return false;
  
  // 一度ulongへキャストし、チェック
  unsigned long* recvBuf = new unsigned long[n];
  for (int i = 0; i < n; i++) sendBuf[i] = (unsigned long)array[i];
  if ( MPI_Allreduce(sendBuf, recvBuf, n, MPI_UNSIGNED_LONG, MPI_SUM, paraMngr->GetMPI_Comm(procGrp)) != MPI_SUCCESS ) return false;
  
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


// #################################################################
/// Allreduceによる総和(整数配列上書き)
bool MonitorCompo::allReduceSum(int* array, int n)
{
  if ( numProc <= 1 ) return true;
  
  unsigned long* sendBuf = new unsigned long[n];
  bool ret = allReduceSum(array, n, sendBuf);
  delete[] sendBuf;
  
  return ret;
}


// #################################################################
/// 指定されたモニタ領域内でスカラー変数を平均
///
///   @param [in] s スカラー変数配列
///   @return モニタ領域内平均値
///
REAL_TYPE MonitorCompo::averageScalar(REAL_TYPE* s)
{
  REAL_TYPE sum= 0.0;
  
  for (int i = 0; i < nPoint; i++)
  {
    if (rank[i] == myRank) sum+= s[i];
  }
  allReduceSum(&sum, 1);
  
  return sum / nPoint;
}


// #################################################################
/// 指定されたモニタ領域内でベクトル変数を平均
///
///   @param [in] v ベクトル変数配列
///   @return モニタ領域内平均値
///
Vec3r MonitorCompo::averageVector(Vec3r* v)
{
  REAL_TYPE sum[3] = { 0.0, 0.0, 0.0 };
  
  for (int i = 0; i < nPoint; i++)
  {
    if (rank[i] == myRank)
    {
      sum[0] += v[i].x;
      sum[1] += v[i].y;
      sum[2] += v[i].z;
    }
  }
  allReduceSum(sum, 3);
  
  return Vec3r(sum) / nPoint;
}


// #################################################################
/// モニタ点の状態を調べ，不正モニタ点フラグ配列pointStatusを設定
void MonitorCompo::checkMonitorPoints()
{
  for (int m = 0; m < nPoint; m++)
  {
    pointStatus[m] = 0;
    if (rank[m] == myRank)
    {
      pointStatus[m] = mon[m]->checkMonitorPoint();
    }
  }
  
  if ( !allReduceSum(pointStatus, nPoint) ) Exit(0);
  
}



// #################################################################
/**
 * @brief セルモニタの場合の交点情報をクリアする
 */
unsigned long MonitorCompo::clearMonitorCut()
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  int odr = polyID;
  
  long long* ct = cut;
  int* bd = bid;
  
  unsigned long c = 0;
  
#pragma omp parallel for firstprivate(ix, jx, kx, gd, odr) schedule(static) reduction(+:c)
  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        
        size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        int qq = bd[m];
        
        if ( IS_CUT(qq) ) // 6面のいずれかにIDがある
        {
          long long pos = ct[m];
          
          size_t m_w = _F_IDX_S3D(i-1, j,   k,   ix, jx, kx, gd);
          size_t m_e = _F_IDX_S3D(i+1, j,   k,   ix, jx, kx, gd);
          size_t m_s = _F_IDX_S3D(i,   j-1, k,   ix, jx, kx, gd);
          size_t m_n = _F_IDX_S3D(i,   j+1, k,   ix, jx, kx, gd);
          size_t m_b = _F_IDX_S3D(i,   j,   k-1, ix, jx, kx, gd);
          size_t m_t = _F_IDX_S3D(i,   j,   k+1, ix, jx, kx, gd);
          
          
          // 隣接セルの方向に対する境界ID
          int qw = getBit5(qq, X_minus);
          int qe = getBit5(qq, X_plus);
          int qs = getBit5(qq, Y_minus);
          int qn = getBit5(qq, Y_plus);
          int qb = getBit5(qq, Z_minus);
          int qt = getBit5(qq, Z_plus);
          
          int rw = getBit5(bd[m_w], X_plus);
          int re = getBit5(bd[m_e], X_minus);
          int rs = getBit5(bd[m_s], Y_plus);
          int rn = getBit5(bd[m_n], Y_minus);
          int rb = getBit5(bd[m_b], Z_plus);
          int rt = getBit5(bd[m_t], Z_minus);

          
          int flag = 0;
          
          REAL_TYPE p_xm = getCut9(pos, X_minus);
          REAL_TYPE p_xp = getCut9(pos, X_plus);
          REAL_TYPE p_ym = getCut9(pos, Y_minus);
          REAL_TYPE p_yp = getCut9(pos, Y_plus);
          REAL_TYPE p_zm = getCut9(pos, Z_minus);
          REAL_TYPE p_zp = getCut9(pos, Z_plus);
          
          int r1 = quantize9(1.0);
          
          // X-
          if ( (p_xm <= 0.5) && (qw == odr) )
          {
            REAL_TYPE dd = 1.0 - getCut9(ct[m_w], X_plus);
            
            if ( (fabs(dd-p_xm) < ROUND_EPS) && ( qw == rw ) ) // 流体にする
            {
              // iセルのX-方向
              setCut9(ct[m], r1, X_minus);
              setBit5(qq, 0, X_minus);
              
              // i-1セルのX+方向
              setCut9(ct[m_w], r1, X_plus);
              setBit5(bd[m_w], 0, X_plus);
            }
            else // 反対側の固体で置き換え
            {
              setCut9(ct[m], quantize9(dd), X_minus);
              setBit5(qq, rw, X_minus);
            }
            flag++;
          }
          
          // X+
          if ( (p_xp <= 0.5) && (qe == odr) )
          {
            REAL_TYPE dd = 1.0 - getCut9(ct[m_e], X_minus);
            
            if ( (fabs(dd-p_xp) < ROUND_EPS) && ( qe == re ) )
            {
              setCut9(ct[m], r1, X_plus);
              setBit5(qq, 0, X_plus);
              
              setCut9(ct[m_e], r1, X_minus);
              setBit5(bd[m_e], 0, X_minus);
            }
            else
            {
              setCut9(ct[m], quantize9(dd), X_plus);
              setBit5(qq, re, X_plus);
            }
            flag++;
          }
          
          // Y-
          if ( (p_ym <= 0.5) && (qs == odr) )
          {
            REAL_TYPE dd = 1.0 - getCut9(ct[m_s], Y_plus);
            
            if ( (fabs(dd-p_ym) < ROUND_EPS) && ( qs == rs ) )
            {
              setCut9(ct[m], r1, Y_minus);
              setBit5(qq, 0, Y_minus);
              
              setCut9(ct[m_s], r1, Y_plus);
              setBit5(bd[m_s], 0, Y_plus);
            }
            else
            {
              setCut9(ct[m], quantize9(dd), Y_minus);
              setBit5(qq, rs, Y_minus);
            }
            flag++;
          }
          
          // Y+
          if ( (p_yp <= 0.5) && (qn == odr) )
          {
            REAL_TYPE dd = 1.0 - getCut9(ct[m_n], Y_minus);
            
            if ( (fabs(dd-p_yp) < ROUND_EPS) && ( qn == rn ) )
            {
              setCut9(ct[m], r1, Y_plus);
              setBit5(qq, 0, Y_plus);
              
              setCut9(ct[m_n], r1, Y_minus);
              setBit5(bd[m_n], 0, Y_minus);
            }
            else
            {
              setCut9(ct[m], quantize9(dd), Y_plus);
              setBit5(qq, rn, Y_plus);
            }
            flag++;
          }
          
          // Z-
          if ( (p_zm <= 0.5) && (qb == odr) )
          {
            REAL_TYPE dd = 1.0 - getCut9(ct[m_b], Z_plus);
            
            if ( (fabs(dd-p_zm) < ROUND_EPS) && ( qb == rb ) )
            {
              setCut9(ct[m], r1, Z_minus);
              setBit5(qq, 0, Z_minus);
              
              setCut9(ct[m_b], r1, Z_plus);
              setBit5(bd[m_b], 0, Z_plus);
            }
            else
            {
              setCut9(ct[m], quantize9(dd), Z_minus);
              setBit5(qq, rb, Z_minus);
            }
            flag++;
          }
          
          // Z+
          if ( (p_zp <= 0.5) && (qt == odr) )
          {
            REAL_TYPE dd = 1.0 - getCut9(ct[m_t], Z_minus);
            
            if ( (fabs(dd-p_zp) < ROUND_EPS) && ( qt == rt ) )
            {
              setCut9(ct[m], r1, Z_plus);
              setBit5(qq, 0, Z_plus);
              
              setCut9(ct[m_t], r1, Z_minus);
              setBit5(bd[m_t], 0, Z_minus);
            }
            else
            {
              setCut9(ct[m], quantize9(dd), Z_plus);
              setBit5(qq, rt, Z_plus);
            }
            flag++;
          }
          
          bd[m] = qq;
          
          if ( flag > 0 ) c++;
        }
        
      }
    }
  }
  
  if ( numProc > 1 )
  {
    unsigned long c_tmp = c;
    if ( paraMngr->Allreduce(&c_tmp, &c, 1, MPI_SUM, procGrp) != CPM_SUCCESS ) Exit(0);
  }
  
  return c;
}



// #################################################################
/// 出力ファイルクローズ
void MonitorCompo::closeFile()
{ 
  assert(fp);
  fclose(fp);
}



// #################################################################
/// サンプリングした変数をノード0に集約
void MonitorCompo::gatherSampled()
{
  REAL_TYPE* vSendBuf = NULL;
  REAL_TYPE* vRecvBuf = NULL;
  REAL_TYPE* sRecvBuf = NULL;
  
  if (variable[var_Velocity])      gatherSampledVector(vel, vSendBuf, vRecvBuf);
  if (variable[var_Pressure])      gatherSampledScalar(prs, sRecvBuf);
  if (variable[var_Temperature])   gatherSampledScalar(tmp, sRecvBuf);
  
  if (variable[var_TotalP])        gatherSampledScalar(tp, sRecvBuf);
  if (variable[var_Helicity])      gatherSampledScalar(hlt, sRecvBuf);
  if (variable[var_Vorticity])     gatherSampledVector(vor, vSendBuf, vRecvBuf);
  
  if (vSendBuf) delete[] vSendBuf;
  if (vRecvBuf) delete[] vRecvBuf;
  if (sRecvBuf) delete[] sRecvBuf;
}


// #################################################################
/// サンプリングしたスカラー変数をノード0に集約
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
      if( MPI_Gather(s, nPoint, MPI_DOUBLE, sRecvBuf, nPoint, MPI_DOUBLE, 0, paraMngr->GetMPI_Comm(procGrp)) != MPI_SUCCESS ) Exit(0);
    }
    else
    {
      if( MPI_Gather(s, nPoint, MPI_FLOAT, sRecvBuf, nPoint, MPI_FLOAT, 0, paraMngr->GetMPI_Comm(procGrp)) != MPI_SUCCESS ) Exit(0);
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


// #################################################################
/// サンプリングしたベクトル変数をノード0に集約
///
///   @param [in,out] v                  ベクトル変数配列
///   @param [in,out] vSendBuf,vRecvBuf  通信用work領域
///
void MonitorCompo::gatherSampledVector(Vec3r* v, REAL_TYPE* vSendBuf, REAL_TYPE* vRecvBuf)
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
  
  for (int m = 0; m < nPoint; m++)
  {
    vSendBuf[3*m  ] = v[m].x;
    vSendBuf[3*m+1] = v[m].y;
    vSendBuf[3*m+2] = v[m].z;
  }
  
  if ( numProc > 1 )
  {
    if ( sizeof(REAL_TYPE) == 8 )
    {
      if( MPI_Gather(vSendBuf, nPoint*3, MPI_DOUBLE, vRecvBuf, nPoint*3, MPI_DOUBLE, 0, paraMngr->GetMPI_Comm(procGrp)) != MPI_SUCCESS ) Exit(0);
    }
    else
    {
      if( MPI_Gather(vSendBuf, nPoint*3, MPI_FLOAT, vRecvBuf, nPoint*3, MPI_FLOAT, 0, paraMngr->GetMPI_Comm(procGrp)) != MPI_SUCCESS ) Exit(0);
    }
  }
  
  if (myRank == 0) 
  {
    for (int i = 0; i < np; i++) {
      for (int m = 0; m < nPoint; m++) {
        if (rank[m] == i)
        {
          v[m].x = vRecvBuf[nPoint*3*i+m*3  ];
          v[m].y = vRecvBuf[nPoint*3*i+m*3+1];
          v[m].z = vRecvBuf[nPoint*3*i+m*3+2];
        }
      }
    }
  }
}



// #################################################################
/// モニタ点指定方法文字列の取得
string MonitorCompo::getTypeStr()
{
  string str;
  
  switch (monitor_type)
  {
    case mon_POINT_SET:
      str = "PointSet";
      break;
      
    case mon_LINE:
      str = "Line";
      break;
      
    case mon_CYLINDER:
      str = "Cylinder";
      break;
      
    case mon_BOX:
      str = "Box";
      break;
      
    case mon_POLYGON:
      str = "Polygon";
      break;
      
    default:
      printf("Monitor Type = %d\n", monitor_type);
      Exit(0);
  }
  
  return str;
}


// #################################################################
/// サンプリング方法文字列の取得
string MonitorCompo::getMethodStr()
{
  string str;
  
  switch (method)
  {
    case SAMPLING_NEAREST:
      str = "Nearest";
      break;
      
    case SAMPLING_INTERPOLATION:
      str = "Interpolation";
      break;
      
    case SAMPLING_SMOOTHING:
      str = "Smoothing";
      break;
      
    default:
      Exit(0);
  }
  
  return str;
}


// #################################################################
/// サンプリングモード文字列の取得
string MonitorCompo::getModeStr()
{
  string str;
  
  switch (mode)
  {
    case SAMPLING_ALL:
      str = "All";
      break;
      
    case SAMPLING_FLUID_ONLY:
      str = "Fluid";
      break;
      
    case SAMPLING_SOLID_ONLY:
      str = "Solid";
      break;
      
    default:
      Exit(0);
  }
  
  return str;
}


// #################################################################
/// モニタ変数を結合した文字列の取得
string MonitorCompo::getVarStr() 
{
  string var;
  
  if (variable[var_Velocity])       var += "Velocity ";
  if (variable[var_Pressure])       var += "Pressure ";
  if (variable[var_Temperature])    var += "Temperature ";
  
  if (variable[var_TotalP])         var += "TotalPressure ";
  if (variable[var_Helicity])       var += "Helicity ";
  if (variable[var_Vorticity])      var += "Vorticity ";
  
  return var;
}



// #################################################################
/// 出力ファイルオープン
///
///    @param str ファイル名テンプレート
///    @param gathered true=gather出力/false=disutribute出力
///
void MonitorCompo::openFile(const char* str, const bool gathered)
{
  if ( (monitor_type == mon_LINE) || (monitor_type == mon_POINT_SET) )
  {
    if (gathered)
    {
      if (!(fp = fopen(str, "w"))) { perror(str); Exit(0); }
      
      writeHeader(true);
      fflush(fp);
    }
    else
    {
      string fileName(str);
      ostringstream rankStr;
      rankStr << "_" << myRank;
      string::size_type pos = fileName.rfind(".");
      fileName.insert(pos, rankStr.str());
      
      if (!(fp = fopen(fileName.c_str(), "w"))) { perror(fileName.c_str()); Exit(0); }
      
      writeHeader(false);
      fflush(fp);
    }
  }
  else // mon_POLYGON, mon_CYLINDER, mon_BOX, mon_PLANE>> gatherd only
  {
    if (!(fp = fopen(str, "w"))) { perror(str); Exit(0); }
    
    writeHeaderCompo();
    fflush(fp);
  }
}


// #################################################################
/// モニタ結果出力
///
///   @param [in] step     サンプリング時の計算ステップ
///   @param [in] tm       サンプリング時の計算時刻
///   @param [in] gathered 出力モードフラグ(true=gather出力/false=disutribute出力)
///
void MonitorCompo::print(unsigned step, double tm, bool gathered)
{
  assert(fp);
  
  if ( (monitor_type == mon_LINE) ||
       (monitor_type == mon_POINT_SET && m_ObjType == mon_POINT_SET)
      )
  {
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
    
    // GraphPloter.py does not allow empty line.
    //fprintf(fp, "\n");
    
    if (refVar.modePrecision == sizeof(float))
    {
      fprintf(fp, "%d %14.6e\n", step, convTime(tm));
    }
    else
    {
      fprintf(fp, "%d %24.16e\n", step, convTime(tm));
    }
    
    for (int i = 0; i < nPoint; i++)
    {
      if (!gathered && rank[i] != myRank) continue;
      
      if (pointStatus[i] != Sampling::POINT_STATUS_OK)
      {
        fprintf(fp, "%s\n", "  *NA*");
        continue;
      }
      
      if (variable[var_Velocity])    fprintf(fp, vFmt, convVel(vel[i].x), convVel(vel[i].y), convVel(vel[i].z));
      if (variable[var_Pressure])    fprintf(fp, sFmt, convPrs(prs[i]));
      if (variable[var_Temperature]) fprintf(fp, sFmt, convTmp(tmp[i]));
      
      if (variable[var_TotalP])      fprintf(fp, sFmt, convTP(tp[i]));
      if (variable[var_Helicity])    fprintf(fp, sFmt, convHlt(hlt[i]));
      if (variable[var_Vorticity])   fprintf(fp, vFmt, convVor(vor[i].x), convVor(vor[i].y), convVor(vor[i].z));
      
      fprintf(fp, "\n");
    }
  }
  else if( m_ObjType == mon_PLANE || monitor_type == mon_PLANE )
  {
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
    
    //GraphPloter.py does not allow empty line.
    //fprintf(fp, "\n");
    
    //int dummy_i = -99999;
    
    if (refVar.modePrecision == sizeof(float))
    {
      fprintf(fp, "%d %14.6e\n", step, convTime(tm));
    }
    else
    {
      fprintf(fp, "%d %24.16e\n", step, convTime(tm));
    }
    
    int m = (int)m_Div[0];//分割数
    int n = (int)m_Div[1];//分割数
    
    int i = -1;
    for( int jj=0; jj<n+1; jj++ )
    {
      for( int ii=0; ii<m+1; ii++ )
      {
        i++;
        
        int pi = plane_grid_flags.at(i);
        
        if( pi < 0 )
        {
          fprintf(fp, "%s\n", "  *NA*");
          continue;
        }
        
        if (variable[var_Velocity])    fprintf(fp, vFmt, convVel(vel[pi].x), convVel(vel[pi].y), convVel(vel[pi].z));
        if (variable[var_Pressure])    fprintf(fp, sFmt, convPrs(prs[pi]));
        if (variable[var_Temperature]) fprintf(fp, sFmt, convTmp(tmp[pi]));
        
        if (variable[var_TotalP])      fprintf(fp, sFmt, convTP(tp[pi]));
        if (variable[var_Helicity])    fprintf(fp, sFmt, convHlt(hlt[pi]));
        if (variable[var_Vorticity])   fprintf(fp, vFmt, convVor(vor[pi].x), convVor(vor[pi].y), convVor(vor[pi].z));
        
        fprintf(fp, "\n");
      }
    }
  }
  else if( m_ObjType == mon_CYLINDER || m_ObjType == mon_BOX || monitor_type == mon_CYLINDER || monitor_type == mon_BOX ) //gatherd only
  {
    if ( !gathered ) Exit(0);
    
    char* sFmtSingle = "%15.7e ";
    char* sFmtDouble = "%24.16e ";
    char* vFmtSingle = "%15.7e %15.7e %15.7e ";
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
    
    if (refVar.modePrecision == sizeof(float))
    {
      fprintf(fp, "%d %14.6e\n", step, convTime(tm));
    }
    else
    {
      fprintf(fp, "%d %24.16e\n", step, convTime(tm));
    }

    Vec3r velo, vort;
    REAL_TYPE pres=0.0, temp=0.0, tolp=0.0, helt=0.0, velx=0.0, vely=0.0, velz=0.0, vorx=0.0, vory=0.0, vorz=0.0;
    int count=0;
    for (int i = 0; i < nPoint; i++)
    {
      if (!gathered && rank[i] != myRank) continue;
      
      if (pointStatus[i] != Sampling::POINT_STATUS_OK)
      {
        fprintf(fp, "%s\n", "  *NA*");
        continue;
      }
      
      if (variable[var_Velocity])     { velx += vel[i].x; vely += vel[i].y; velz += vel[i].z; }
      if (variable[var_Pressure])     pres = pres + convPrs(prs[i]);
      if (variable[var_Temperature])  temp = temp + convTmp(tmp[i]);
      
      if (variable[var_TotalP])       tolp = tolp + convTP(tp[i]);
      if (variable[var_Helicity])     helt = helt + convHlt(hlt[i]);
      if (variable[var_Vorticity])   { vorx += vor[i].x; vory += vor[i].y; vorz += vor[i].z; }
      
      count++;
    }
    
    if( count < 0 ){fflush(fp); return;}
    
    velo = velo / (REAL_TYPE) count;
    vort = vort / (REAL_TYPE) count;
    pres /= (REAL_TYPE) count;
    temp /= (REAL_TYPE) count;
    tolp /= (REAL_TYPE) count;
    helt /= (REAL_TYPE) count;
    
    if (variable[var_Velocity])    fprintf(fp, vFmt, convVel(velo.x), convVel(velo.y), convVel(velo.z));
    if (variable[var_Pressure])    fprintf(fp, sFmt, pres);
    if (variable[var_Temperature]) fprintf(fp, sFmt, temp);
    
    if (variable[var_TotalP])      fprintf(fp, sFmt, tolp);
    if (variable[var_Helicity])    fprintf(fp, sFmt, helt);
    if (variable[var_Vorticity])   fprintf(fp, vFmt, convVor(vort.x), convVor(vort.y), convVor(vort.z));
    
    fprintf(fp, "\n");
  }
  else if( monitor_type == mon_POLYGON )
  {
    int aaa=0;
  }

  fflush(fp);
}


// #################################################################
/// モニタ情報を出力
///
///    @param [in] fp 出力ファイルポインタ
///    @param [in] no モニタグループ通し番号
//
void MonitorCompo::printInfo(FILE* fp, int no)
{
  fprintf(fp,"\t%3d : %s\t division=%d  [%s]\n", no+1, getTypeStr().c_str(), nPoint, label.c_str());
  fprintf(fp,"\t\tVariables : %s\n", getVarStr().c_str());
  fprintf(fp,"\t\t   Method : %s\n", getMethodStr().c_str());
  fprintf(fp,"\t\t     Mode : %s\n\n", getModeStr().c_str());
  
  fprintf(fp,"\t\t    order :            X              Y              Z    :   rank : comment\n");
  for (int j = 0; j < nPoint; j++)
  {
    fprintf(fp,"\t\t%9d : %14.6e %14.6e %14.6e  : %6d : %s", // %12.4 >> %14.6
            j+1, convCrd(crd[j].x), convCrd(crd[j].y), convCrd(crd[j].z), rank[j], comment[j].c_str());
    
    if (pointStatus[j] == Sampling::UNEXPECTED_SOLID)      // 流体セルを指定したが固体だった
    {
      fprintf(fp, "  *skip(unexpected solid)*\n");
    }
    else if (pointStatus[j] == Sampling::UNEXPECTED_FLUID) // 固体セルを指定したが流体だった
    {
      fprintf(fp, "  *skip(unexpected fluid)*\n");
    }
    else
    {
      fprintf(fp, "\n");
    }
  }
  fprintf(fp,"\n\n");
}


// #################################################################
/// サンプリング(Line, PointSet)
void MonitorCompo::sampling()
{
  for (int i = 0; i < nPoint; i++)
  {
    //  if (!(mon[i] && pointStatus[i] == Sampling::POINT_STATUS_OK)) continue;
    if (!mon[i]) continue;

    if (variable[var_Velocity])     vel[i] = mon[i]->samplingVelocity(vSource);
    if (variable[var_Pressure])     prs[i] = mon[i]->samplingPressure(pSource);
    if (variable[var_Temperature])  tmp[i] = mon[i]->samplingTemperature(tSource);

    if (variable[var_TotalP])       tp[i]  = mon[i]->samplingTotalPressure(vSource, pSource);
    if (variable[var_Vorticity])    vor[i] = mon[i]->samplingVorticity(vSource);
    if (variable[var_Helicity])     hlt[i] = mon[i]->samplingHelicity(vSource);
  }
}


// #################################################################
/// モニタ点での平均値サンプリング
///
///   サンプリング結果を集計，領域での平均値を計算
///   速度は法線ベクトルとの内積をとる
///
void MonitorCompo::samplingAverage()
{
  sampling();
  
  if (variable[var_Velocity]) 
  {
    Vec3r velAve = averageVector(vel);
    val[var_Velocity] = velAve.x * nv[0] + velAve.y * nv[1] + velAve.z * nv[2];
  }
  if (variable[var_Pressure])     val[var_Pressure]    = averageScalar(prs);
  if (variable[var_Temperature])  val[var_Temperature] = averageScalar(tmp);
  
  
  if (variable[var_Vorticity])
  {
    Vec3r vrtAve = averageVector(vor);
    val[var_Vorticity] = vrtAve.x * nv[0] + vrtAve.y * nv[1] + vrtAve.z * nv[2];
  }
  if (variable[var_TotalP])       val[var_TotalP]      = averageScalar(tp);
  if (variable[var_Helicity])     val[var_Helicity]    = averageScalar(hlt);
}


// #################################################################
/// Line登録
void MonitorCompo::setLine(const char* labelStr,
                           vector<string>& variables,
                           const char* methodStr,
                           const char* modeStr,
                           const REAL_TYPE from[3],
                           const REAL_TYPE to[3],
                           const int nDivision,
                           Monitor_Type m_type)
{
  label = labelStr;
  setSamplingMethod(methodStr);
  setSamplingMode(modeStr);
  monitor_type = m_type;
  
  vector<string>::const_iterator it;
  for (it = variables.begin(); it != variables.end(); it++)
  {
    setMonitorVar((*it).c_str());
  }
  
  nPoint = nDivision + 1;
  if (nPoint < 2) Exit(0);
  
  // アロケート
  allocArray();
  allocSamplingArray();
  
  Vec3r st(from);
  Vec3r ed(to);
  Vec3r dd = ed - st;
  dd /= (REAL_TYPE)nPoint - 1.0;
  
  
  for (int m = 0; m < nPoint; m++)
  {
    ostringstream oss;
    oss << "point_" << m;
    
    crd[m] = st + dd * (REAL_TYPE)m;
    
    // 計算領域全体でのチェック
    if (!checkRegion(m, g_org, g_box, true)) Exit(0);
    
    comment[m] = oss.str();
  }
  
  // サンプリングポイントの担当ランクを決める
  setRankArray();
  
  // モニタ点の登録
  for (int m = 0; m < nPoint; m++)
  {
    if (rank[m] == myRank) 
    {
      switch (method)
      {
        case SAMPLING_NEAREST:
          mon[m] = new Nearest(mode, size, guide, crd[m], org, pch, refVar.v00, bcd, NoCompo, mtbl);
          break;
          
        case SAMPLING_INTERPOLATION:
          mon[m] = new Interpolation(mode, size, guide, crd[m], org, pch, refVar.v00, bcd, NoCompo, mtbl);
          break;
          
        case SAMPLING_SMOOTHING:
          mon[m] = new Smoothing(mode, size, guide, crd[m], org, pch, refVar.v00, bcd, NoCompo, mtbl);
          break;
          
        default:
          Exit(0);
      }
    }
  }
}


// #################################################################
/// モニタ対象物理量の設定
///
///    @param[in] str モニタ対象物理量文字列
///
void MonitorCompo::setMonitorVar(const char* str) 
{
  if      (!strcasecmp("velocity", str))             variable[var_Velocity]       = true;
  else if (!strcasecmp("pressure", str))             variable[var_Pressure]       = true;
  else if (!strcasecmp("temperature", str))          variable[var_Temperature]    = true;

  else if (!strcasecmp("totalpressure", str))        variable[var_TotalP]         = true;
  else if (!strcasecmp("helicity", str))             variable[var_Helicity]       = true;
  else if (!strcasecmp("vorticity", str))            variable[var_Vorticity]      = true;
  else
  {
    Hostonly_ stamped_printf("\tError : Invalid variable keyword [%s]\n", str);
    Exit(0);
  }
}


// #################################################################
/// PointSet登録
void MonitorCompo::setPointSet(const char* labelStr,
                               vector<string>& variables,
                               const char* methodStr,
                               const char* modeStr,
                               vector<MonitorPoint>& pointSet,
                               Monitor_Type m_type)
{
  label = labelStr;
  setSamplingMethod(methodStr);
  setSamplingMode(modeStr);
  monitor_type = m_type;
  
  vector<string>::const_iterator it;
  for (it = variables.begin(); it != variables.end(); it++)
  {
    setMonitorVar((*it).c_str());
  }
  
  nPoint = pointSet.size();
  
  if (nPoint == 0)
  {
    printf("\tError : No sampling point\n");
    Exit(0);
  }
  
  // アロケート
  allocArray();
  allocSamplingArray();
  
  for (int m = 0; m < nPoint; m++)
  {
    crd[m] = pointSet[m].crd;
    
    // 計算領域全体でのチェック
    if (!checkRegion(m, g_org, g_box, true))
    {
      Exit(0);
    }
    
    comment[m] = pointSet[m].label;
  }
  
  // サンプリングポイントの担当ランクを決める
  setRankArray();
  
  // モニタ点の登録
  for (int m = 0; m < nPoint; m++)
  {
    if (rank[m] == myRank) 
    {
      switch (method)
      {
        case SAMPLING_NEAREST:
          mon[m] = new Nearest(mode, size, guide, crd[m], org, pch, refVar.v00, bcd, NoCompo, mtbl);
          break;
          
        case SAMPLING_INTERPOLATION:
          mon[m] = new Interpolation(mode, size, guide, crd[m], org, pch, refVar.v00, bcd, NoCompo, mtbl);
          break;
          
        case SAMPLING_SMOOTHING:
          mon[m] = new Smoothing(mode, size, guide, crd[m], org, pch, refVar.v00, bcd, NoCompo, mtbl);
          break;
          
        default:
          Exit(0);
      }
    }
  }
  
}


// #################################################################
/// Polygon登録
void MonitorCompo::setPolygon(const char* labelStr,
                              vector<string>& variables,
                              const char* methodStr,
                              const char* modeStr,
                              const int order,
                              const REAL_TYPE m_nv[3],
                              Monitor_Type m_type)
{
  label = labelStr;
  setSamplingMethod(methodStr);
  setSamplingMode(modeStr);
  monitor_type = m_type;
  nv[0] = m_nv[0];
  nv[1] = m_nv[1];
  nv[2] = m_nv[2];
  polyID = order;
  int odr = order;

  
  vector<string>::const_iterator it;
  for (it = variables.begin(); it != variables.end(); it++)
  {
    setMonitorVar((*it).c_str());
  }
  

  // ローカルのサンプリング点数 nPointList[myRank] を求める
  int np = num_process;
  
  int* nPointList;
  if (!(nPointList = new int[np])) Exit(0);
  for (int i = 0; i < np; i++) nPointList[i] = 0;
  
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
  // Serial
  for (int k = 1; k <= kx; k++) {
    for (int j = 1; j <= jx; j++) {
      for (int i = 1; i <= ix; i++) {
        
        size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        int bd = bid[m];
        
        if ( IS_CUT(bd) ) // 6面のいずれかにIDがある
        {
          const long long pos = cut[m];
          
          for (int i=0; i<6; i++)
          {
            int d = (bd >> i*5) & MASK_5;
            if ( (getCut9(pos, i) <= 0.5) && (d == odr) ) nPointList[myRank]++; // セル内部にカットが存在し，境界IDがエントリ番号
          }
        }
      }
    }
  }
  
  
  // 全サンプリング数を求める
  if ( !allReduceSum(nPointList, np) )
  {
    if (nPointList) delete [] nPointList;
    nPointList = NULL;
    Exit(0);
  }
  
  int sum = 0;
  for (int i = 0; i < np; i++) sum += nPointList[i];
  
  // Number of sampling points
  nPoint = sum;
  
  
#if 0
  printf("Polygon Monitor[rank=%d] : sum=%d LocalPoint=%d\n", myRank, sum, nPointList[myRank]);
#endif
  
  if (nPoint == 0)
  {
    printf("\tError : No sampling point\n");
    Exit(0);
  }

  
  // アロケート
  allocArray();
  allocSamplingArray();
  
  
  // 座標値を保持するための配列，バッファ共用
  REAL_TYPE* buf;
  if (!(buf = new REAL_TYPE[nPoint*3])) Exit(0);
  
  // 初期値にゼロをいれておき，MPI_SUMで集約
  for (int m = 0; m < nPoint*3; m++) buf[m] = 0.0;
  
  // 各ランクの開始点offset
  int offset = 0;
  for (int i = 0; i < myRank; i++) offset += nPointList[i];
  
  
  // 座標値を計算
  int lp = 0;
  
  // serial
  for (int k = 1; k <= kx; k++) {
    for (int j = 1; j <= jx; j++) {
      for (int i = 1; i <= ix; i++) {
        
        size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        int bd = bid[m];
        int flag = 0;
        
        float cx = 0.0;
        float cy = 0.0;
        float cz = 0.0;
        
        if ( IS_CUT(bd) ) // 6面のいずれかにIDがある
        {
          const long long pos = cut[m];
          
          for (int i=0; i<6; i++)
          {
            int d = (bd >> i*5) & MASK_5;
            
            REAL_TYPE pp = getCut9(pos, i);
            
            if ( (pp <= 0.5) && (d == odr) ) // セル内部に存在する
            {
              if      (i == 0) cx = -pp;
              else if (i == 1) cx =  pp;
              else if (i == 2) cy = -pp;
              else if (i == 3) cy =  pp;
              else if (i == 4) cz = -pp;
              else if (i == 5) cz =  pp;
              flag++;
            }
          }
        }
        
        if ( flag > 0 )
        {
          buf[(offset+lp)*3+0] = org.x + (i - 0.5 + cx) * pch.x;
          buf[(offset+lp)*3+1] = org.y + (j - 0.5 + cy) * pch.y;
          buf[(offset+lp)*3+2] = org.z + (k - 0.5 + cz) * pch.z;
          lp++;
        }
        
      }
    }
  }
  
  if ( lp != nPointList[myRank] )
  {
    printf("Rank[%d] : number of sampling points is inconsistent (%d, %d)\n", myRank, lp, nPointList[myRank]);
    Exit(0);
  }
  

  
  // 座標値を集める
  if (!allReduceSum(buf, nPoint*3)) Exit(0);
  
  
  // serial
  for (int m = 0; m < nPoint; m++)
  {
    ostringstream oss;
    oss << "point_" << m;
    
    crd[m].x = buf[m*3+0];
    crd[m].y = buf[m*3+1];
    crd[m].z = buf[m*3+2];
    
    // 計算領域全体でのチェック
    if (!checkRegion(m, g_org, g_box, true)) Exit(0);
    
    comment[m] = oss.str();
  }

  
  // サンプリングポイントの担当ランクを決める
  setRankArray();
  
  
  // モニタ点の登録
  for (int m = 0; m < nPoint; m++)
  {
    if (rank[m] == myRank)
    {
      switch (method)
      {
        case SAMPLING_NEAREST:
          mon[m] = new Nearest(mode, size, guide, crd[m], org, pch, refVar.v00, bcd, NoCompo, mtbl);
          break;
          
        case SAMPLING_INTERPOLATION:
          mon[m] = new Interpolation(mode, size, guide, crd[m], org, pch, refVar.v00, bcd, NoCompo, mtbl);
          break;
          
        case SAMPLING_SMOOTHING:
          mon[m] = new Smoothing(mode, size, guide, crd[m], org, pch, refVar.v00, bcd, NoCompo, mtbl);
          break;
          
        default:
          Exit(0);
      }
    }
  }

  // テンポラリバッファの削除
  if (buf) delete [] buf; buf = NULL;
  if (nPointList) delete [] nPointList; nPointList = NULL;
  
}


// #################################################################
/// Primitiveを登録し，bcd[]をゼロクリア
void MonitorCompo::setPrimitive(const char* labelStr,
                                vector<string>& variables,
                                const char* methodStr,
                                const char* modeStr,
                                const int order,
                                const REAL_TYPE m_nv[3],
                                Monitor_Type m_type)
{
  label = labelStr;
  setSamplingMethod(methodStr);
  setSamplingMode(modeStr);
  monitor_type = m_type;
  nv[0] = m_nv[0];
  nv[1] = m_nv[1];
  nv[2] = m_nv[2];
  int odr = order;
  
  vector<string>::const_iterator it;
  for (it = variables.begin(); it != variables.end(); it++)
  {
    setMonitorVar((*it).c_str());
  }
  
  
  // ローカルのサンプリング点数 nPointList[myRank] を求める
  int np = num_process;
  
  int* nPointList;
  if (!(nPointList = new int[np])) Exit(0);
  for (int i = 0; i < np; i++) nPointList[i] = 0;
  
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
  // serial
  for (int k = 1; k <= kx; k++) {
    for (int j = 1; j <= jx; j++) {
      for (int i = 1; i <= ix; i++) {

        if ( DECODE_CMP( bcd[_F_IDX_S3D(i, j, k, ix, jx, kx, gd)] ) == odr )
        {
          nPointList[myRank]++;
        }
        
      }
    }
  }
  
  
  // 全サンプリング数を求める
  if ( !allReduceSum(nPointList, np) )
  {
    if (nPointList) delete [] nPointList;
    nPointList = NULL;
    Exit(0);
  }
  
  int sum = 0;
  for (int i = 0; i < np; i++) sum += nPointList[i];
  
  // Number of sampling points
  nPoint = sum;
  
  
#if 0
  printf("Polygon Monitor[rank=%d] : sum=%d LocalPoint=%d\n", myRank, sum, nPointList[myRank]);
#endif
  
  
  if (nPoint == 0)
  {
    printf("\tError : No sampling point\n");
    Exit(0);
  }
  
  
  // アロケート
  allocArray();
  allocSamplingArray();
  
  
  // 座標値を保持
  REAL_TYPE* buf;
  if (!(buf = new REAL_TYPE[nPoint*3])) Exit(0);

  // 初期値にゼロをいれておき，MPI_SUMで集約
  for (int m = 0; m < nPoint*3; m++) buf[m] = 0.0;
  
  // 各ランクの開始点offset
  int offset = 0;
  for (int i = 0; i < myRank; i++) offset += nPointList[i];
  
  
  // 座標値を計算
  int lp = 0;
  
  // serial
  for (int k = 1; k <= kx; k++) {
    for (int j = 1; j <= jx; j++) {
      for (int i = 1; i <= ix; i++) {
        
        size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        
        if ( DECODE_CMP( bcd[m] ) == odr )
        {
          buf[(offset+lp)*3+0] = org.x + (i - 0.5) * pch.x;
          buf[(offset+lp)*3+1] = org.y + (j - 0.5) * pch.y;
          buf[(offset+lp)*3+2] = org.z + (k - 0.5) * pch.z;
          lp++;
          
          setMediumID(bcd[m], 0); // クリア
        }
        
      }
    }
  }
  
  if ( lp != nPointList[myRank] )
  {
    printf("Rank[%d] : number of sampling points is inconsistent (%d, %d)\n", myRank, lp, nPointList[myRank]);
    Exit(0);
  }
  
  
  
  // 座標値を集める
  if (!allReduceSum(buf, nPoint*3)) Exit(0);
  
  
  // serial
  for (int m = 0; m < nPoint; m++)
  {
    ostringstream oss;
    oss << "point_" << m;
    
    crd[m].x = buf[m*3+0];
    crd[m].y = buf[m*3+1];
    crd[m].z = buf[m*3+2];
    
    // 計算領域全体でのチェック
    if (!checkRegion(m, g_org, g_box, true)) Exit(0);
    
    comment[m] = oss.str();
  }

  
  // サンプリングポイントの担当ランクを決める
  setRankArray();
  
  
  // モニタ点の登録
  for (int m = 0; m < nPoint; m++)
  {
    if (rank[m] == myRank)
    {
      switch (method)
      {
        case SAMPLING_NEAREST:
          mon[m] = new Nearest(mode, size, guide, crd[m], org, pch, refVar.v00, bcd, NoCompo, mtbl);
          break;
          
        case SAMPLING_INTERPOLATION:
          mon[m] = new Interpolation(mode, size, guide, crd[m], org, pch, refVar.v00, bcd, NoCompo, mtbl);
          break;
          
        case SAMPLING_SMOOTHING:
          mon[m] = new Smoothing(mode, size, guide, crd[m], org, pch, refVar.v00, bcd, NoCompo, mtbl);
          break;
          
        default:
          Exit(0);
      }
    }
  }

  // テンポラリバッファの削除
  if (buf) delete [] buf; buf = NULL;
  if (nPointList) delete [] nPointList; nPointList = NULL;
  
}


// #################################################################
/// 各モニタ点を担当するランク番号を配列rank[]にセット
///
///   @note 領域境界上のモニタ点は，ランク番号の大きい方の領域が担当
///
void MonitorCompo::setRankArray()
{
  int* sendBuf =  new int[nPoint];
  if (!sendBuf) Exit(0);
  
  for (int m = 0; m < nPoint; m++)
  {
    sendBuf[m] = -1;
    
    // ローカルノードの担当領域をチェック
    if (checkRegion(m, org, box)) sendBuf[m] = myRank;
  }
  
  // gather max rank number
  if ( numProc > 1 ) 
  {
    if ( MPI_Allreduce(sendBuf, rank, nPoint, MPI_INT, MPI_MAX, paraMngr->GetMPI_Comm(procGrp)) != MPI_SUCCESS ) Exit(0);
  }
  else 
  {
    for (int m = 0; m < nPoint; m++) rank[m] = sendBuf[m];
  }
  
  if (sendBuf) delete[] sendBuf;
}



// #################################################################
/// サンプリング方法の設定
///
///    @param[in] str サンプリング方法文字列
///
void MonitorCompo::setSamplingMethod(const char* str) 
{
  if      (!strcasecmp("nearest", str))       method = SAMPLING_NEAREST;
  else if (!strcasecmp("interpolation", str)) method = SAMPLING_INTERPOLATION;
  else if (!strcasecmp("smoothing", str))     method = SAMPLING_SMOOTHING;
  else
  {
    Hostonly_ stamped_printf("\tError : Invalid samping_method keyword [%s]\n", str);
    Exit(0);
  }
}


// #################################################################
/// サンプリングモードの設定
///
///    @param[in] str サンプリングモード文字列
///
void MonitorCompo::setSamplingMode(const char* str) 
{
  if      (!strcasecmp("all",   str)) mode = SAMPLING_ALL;
  else if (!strcasecmp("fluid", str)) mode = SAMPLING_FLUID_ONLY;
  else if (!strcasecmp("solid", str)) mode = SAMPLING_SOLID_ONLY;
  else
  {
    Hostonly_ stamped_printf("\tError : Invalid samping_mode keyword [%s]\n", str);
    Exit(0);
  }
}

// #################################################################
/// Graph Ploter
void MonitorCompo::setPlaneData(REAL_TYPE c[3], REAL_TYPE z[3], REAL_TYPE x[3], REAL_TYPE dim[2], REAL_TYPE div[2] )
{
  VEC3_EQUATE( m_Center, c );
  VEC3_EQUATE( m_MainDir, z );
  VEC3_EQUATE( m_RefDir, x );
  VEC2_EQUATE( m_Dim2, dim );
  VEC2_EQUATE( m_Div, div );
  setObjType(mon_PLANE);
  
}

// #################################################################
/// Graph Ploter
void MonitorCompo::setCylinderData(REAL_TYPE c[3], REAL_TYPE z[3], REAL_TYPE x[3], REAL_TYPE dim[3])
{
  VEC3_EQUATE( m_Center, c );
  VEC3_EQUATE( m_MainDir, z );
  VEC3_EQUATE( m_RefDir, x );
  VEC3_EQUATE( m_Dim3, dim );
  setObjType(mon_CYLINDER);
}

// #################################################################
/// Graph Ploter
void MonitorCompo::setBoxData(REAL_TYPE c[3], REAL_TYPE z[3], REAL_TYPE x[3], REAL_TYPE dim[3])
{
  VEC3_EQUATE( m_Center, c );
  VEC3_EQUATE( m_MainDir, z );
  VEC3_EQUATE( m_RefDir, x );
  VEC3_EQUATE( m_Dim3, dim );
  setObjType(mon_BOX);
}

// #################################################################
/// モニタ結果出力ファイルにヘッダ部を出力
///
///   @param[in] gathered 出力モードフラグ(true=gather出力/false=disutribute出力)
///
void MonitorCompo::writeHeader(bool gathered)
{
  assert(fp);
  
  int n;
  if (gathered)
  {
    n = nPoint;
  }
  else
  {
    n = 0;
    for (int i = 0; i < nPoint; i++) if (rank[i] == myRank) n++;
  }
  
  std::string type  = getTypeStr();
  std::string label = getLabel();
  
  std::string s_val = "", s_comp = "";
  int n_var = 0;
  if (variable[var_Velocity])     { n_var++; s_val +=" Velocity";     s_comp += " 3";}
  if (variable[var_Pressure])     { n_var++; s_val +=" TotalPresure"; s_comp += " 1";}
  if (variable[var_Temperature])  { n_var++; s_val +=" Temperature";  s_comp += " 1";}
  if (variable[var_TotalP])       { n_var++; s_val +=" TotalIP";      s_comp += " 1";}
  if (variable[var_Helicity])     { n_var++; s_val +=" Helicity";     s_comp += " 1";}
  if (variable[var_Vorticity])    { n_var++; s_val +=" Vorticity";    s_comp += " 3";}
  
  if( monitor_type == mon_LINE )
  {
    REAL_TYPE p1[3] = { convCrd(crd[0].x), convCrd(crd[0].y), convCrd(crd[0].z) };
    REAL_TYPE p2[3] = { convCrd(crd[nPoint-1].x), convCrd(crd[nPoint-1].y), convCrd(crd[nPoint-1].z) };
    
    fprintf(fp, "Line    %s\n", label.c_str());
    fprintf(fp, "%d     %14.6e %14.6e %14.6e %14.6e %14.6e %14.6e\n", n, p1[0],p1[1],p1[2],  p2[0],p2[1],p2[2]);
    fprintf(fp, "%d %s %s\n", n_var, s_val.c_str(), s_comp.c_str() );
    
    for (int i = 0; i < nPoint; i++)
    {
      if (gathered || rank[i] == myRank)
      {
        //例： 0.000000e+00   0.000000e+00   0.000000e+00
        fprintf(fp, "%14.6e %14.6e %14.6e %s", convCrd(crd[i].x), convCrd(crd[i].y), convCrd(crd[i].z), comment[i].c_str());
        
        if (pointStatus[i] == Sampling::UNEXPECTED_SOLID)      // 流体セルを指定したが固体だった
        {
          fprintf(fp, "  *skip(unexpected solid)*\n");
        }
        else if (pointStatus[i] == Sampling::UNEXPECTED_FLUID) // 固体セルを指定したが流体だった
        {
          fprintf(fp, "  *skip(unexpected fluid)*\n");
        }
        else
        {
          fprintf(fp, "\n");
        }
      }
    }
  }
  else if( monitor_type == mon_POINT_SET &&  m_ObjType == mon_POINT_SET )
  {
    fprintf(fp, "Pointset    %s\n", label.c_str());
    fprintf(fp, "%d\n", n );
    fprintf(fp, "%d %s %s\n", n_var, s_val.c_str(), s_comp.c_str() );
    
    for (int i = 0; i < nPoint; i++)
    {
      if (gathered || rank[i] == myRank)
      {
        //例： 0.000000e+00   0.000000e+00   0.000000e+00
        fprintf(fp, "%14.6e %14.6e %14.6e %s", convCrd(crd[i].x), convCrd(crd[i].y), convCrd(crd[i].z), comment[i].c_str());
        
        if (pointStatus[i] == Sampling::UNEXPECTED_SOLID)      // 流体セルを指定したが固体だった
        {
          fprintf(fp, "  *skip(unexpected solid)*\n");
        }
        else if (pointStatus[i] == Sampling::UNEXPECTED_FLUID) // 固体セルを指定したが流体だった
        {
          fprintf(fp, "  *skip(unexpected fluid)*\n");
        }
        else
        {
          fprintf(fp, "\n");
        }
      }
    }
  }
  else if( m_ObjType ==mon_PLANE ||  monitor_type == mon_PLANE )
  {
    fprintf(fp, "Plane    %s\n", label.c_str());
    fprintf(fp, "%14.6e %14.6e %14.6e\n", m_Center[0], m_Center[1], m_Center[2]);     //Center
    fprintf(fp, "%14.6e %14.6e %14.6e\n", m_MainDir[0], m_MainDir[1], m_MainDir[2]);  //MainDirection
    fprintf(fp, "%14.6e %14.6e %14.6e\n", m_RefDir[0], m_RefDir[1], m_RefDir[2]);     //RefDirection
    fprintf(fp, "%d %d\n", (int)m_Div[0]+1, (int)m_Div[1]+1);                         //uGrid, vGrid
    fprintf(fp, "%14.6e %14.6e\n", m_Dim2[0]/m_Div[0], m_Dim2[1]/m_Div[1]);           //uSize, vSize
    fprintf(fp, "%d %s %s\n", n_var, s_val.c_str(), s_comp.c_str() );
    
    int m = (int)m_Div[0];//分割数
    int n = (int)m_Div[1];//分割数
    REAL_TYPE x0 = - m_Dim2[0] * 0.5;
    REAL_TYPE y0 = - m_Dim2[1] * 0.5;
    REAL_TYPE dx = m_Dim2[0] / m;
    REAL_TYPE dy = m_Dim2[1] / n;
    
    Vec3r orig_o( m_Center[0], m_Center[1], m_Center[2] );
    Vec3r axis_z( m_MainDir[0], m_MainDir[1], m_MainDir[2] );
    Vec3r axis_x( m_RefDir[0], m_RefDir[1], m_RefDir[2] );
    Vec3r axis_y = cross(axis_z, axis_x);
    axis_y.normalize();
    REAL_TYPE axis_yy[3] = {axis_y.x, axis_y.y, axis_y.z};
    
    for( int j=0; j<n+1; j++ )
    {
      REAL_TYPE y = y0 + j * dy;
      for( int i=0; i<m+1; i++ )
      {
        REAL_TYPE x = x0 + i * dx;
        Vec3r local_pt(x, y, 0.0);
        Vec3r gp = globalPt(orig_o, axis_z, axis_x, axis_y, local_pt);
        
        fprintf(fp, "%14.6e %14.6e %14.6e\n", convCrd(gp.x), convCrd(gp.y), convCrd(gp.z));
      }
    }
  }
  else if( m_ObjType ==mon_BOX || m_ObjType ==mon_CYLINDER || monitor_type == mon_BOX  || monitor_type == mon_CYLINDER )
  {
    fprintf(fp, "Pointet    %s\n", label.c_str());
    fprintf(fp, "%d\n", 1 );
    fprintf(fp, "%d %s %s\n", n_var, s_val.c_str(), s_comp.c_str() );
    
    int count=0;
    Vec3r sum_vec, vec;
    for (int i = 0; i < nPoint; i++)
    {
      if (gathered || rank[i] == myRank)
      {
        Vec3r vec(crd[i].x, crd[i].y, crd[i].z);
        sum_vec = sum_vec + vec;
        count++;
      }
    }
    vec = sum_vec / (REAL_TYPE)count;
    fprintf(fp, "%14.6e %14.6e %14.6e   %s\n", convCrd(vec.x), convCrd(vec.y), convCrd(vec.z), "averaged_pt");
  }
  
  fflush(fp);
}


// #################################################################
/// 平均値モニタ結果出力ファイルにヘッダ部を出力
void MonitorCompo::writeHeaderCompo()
{
  if ( refVar.modeUnitOutput == DIMENSIONAL )
  {
    fprintf(fp, "      step      time[sec]");
  }
  else
  {
    fprintf(fp, "      step        time[-]");
  }
  
  if (refVar.modePrecision == sizeof(float))
  {
    if ( variable[var_Velocity] )     fprintf(fp, "     Velocity [m/s]");
    if ( variable[var_Pressure] )     fprintf(fp, "      Pressure [pa]");
    if ( variable[var_Temperature] )  fprintf(fp, "    Temperature [C]");
    if ( variable[var_TotalP] )       fprintf(fp, "  TotalPressure[pa]");
    if ( variable[var_Helicity] )     fprintf(fp, "     Helicity [m/s]");
    if ( variable[var_Vorticity] )    fprintf(fp, "    Vorticity [m/s]");
  }
  else
  {
    if ( variable[var_Velocity] )     fprintf(fp, "              Velocity [m/s]");
    if ( variable[var_Pressure] )     fprintf(fp, "               Pressure [pa]");
    if ( variable[var_Temperature] )  fprintf(fp, "             Temperature [C]");
    if ( variable[var_TotalP] )       fprintf(fp, "           TotalPressure[pa]");
    if ( variable[var_Helicity] )     fprintf(fp, "              Helicity [m/s]");
    if ( variable[var_Vorticity] )    fprintf(fp, "             Vorticity [m/s]");
  }
  
  fprintf(fp, "\n");
}
