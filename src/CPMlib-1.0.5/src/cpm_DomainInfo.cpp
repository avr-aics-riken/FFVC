/*
 * CPM - Cartesian Partition Manager
 *
 * Copyright (C) 2012 University of Tokyo.
 *
 */

/**
 * @file   cpm_DomainInfo.cpp
 * DomainInfoクラスのソースファイル
 * @author University of Tokyo
 * @date   2012/05/31
 */
#include "cpm_DomainInfo.h"

////////////////////////////////////////////////////////////////////////////////
// コンストラクタ
cpm_DomainInfo::cpm_DomainInfo()
  : cpm_Base()
{
  clear();
}

////////////////////////////////////////////////////////////////////////////////
// デストラクタ
cpm_DomainInfo::~cpm_DomainInfo()
{
}

////////////////////////////////////////////////////////////////////////////////
// 情報のクリア
void
cpm_DomainInfo::clear()
{
  REAL_TYPE RZERO = REAL_TYPE(0);
  for( int i=0;i<3;i++ )
  {
    m_origin[i] = RZERO;
    m_region[i] = RZERO;
    m_pitch[i]  = RZERO;
    m_voxNum[i] = 0;
  }
}

////////////////////////////////////////////////////////////////////////////////
// 原点のセット
void
cpm_DomainInfo::SetOrigin( REAL_TYPE org[3] )
{
  m_origin[0] = org[0];
  m_origin[1] = org[1];
  m_origin[2] = org[2];
}

////////////////////////////////////////////////////////////////////////////////
// 原点の取得
const REAL_TYPE*
cpm_DomainInfo::GetOrigin() const
{
  return m_origin;
}

////////////////////////////////////////////////////////////////////////////////
// ピッチのセット
void
cpm_DomainInfo::SetPitch( REAL_TYPE pch[3] )
{
  m_pitch[0] = pch[0];
  m_pitch[1] = pch[1];
  m_pitch[2] = pch[2];
}

////////////////////////////////////////////////////////////////////////////////
// ピッチの取得
const REAL_TYPE*
cpm_DomainInfo::GetPitch() const
{
  return m_pitch;
}

////////////////////////////////////////////////////////////////////////////////
// 空間サイズのセット
void
cpm_DomainInfo::SetRegion( REAL_TYPE rgn[3] )
{
  m_region[0] = rgn[0];
  m_region[1] = rgn[1];
  m_region[2] = rgn[2];
}

////////////////////////////////////////////////////////////////////////////////
// 空間サイズの取得
const REAL_TYPE*
cpm_DomainInfo::GetRegion() const
{
  return m_region;
}

////////////////////////////////////////////////////////////////////////////////
// VOXEL数のセット
void
cpm_DomainInfo::SetVoxNum( int vox[3] )
{
  m_voxNum[0] = vox[0];
  m_voxNum[1] = vox[1];
  m_voxNum[2] = vox[2];
}

////////////////////////////////////////////////////////////////////////////////
// VOXEL数の取得
const int*
cpm_DomainInfo::GetVoxNum() const
{
  return m_voxNum;
}

////////////////////////////////////////////////////////////////////////////////
// 領域情報のチェック
cpm_ErrorCode cpm_DomainInfo::CheckData()
{
  REAL_TYPE RZERO = REAL_TYPE(0);

  if( m_region[0] <= RZERO || m_region[1] <= RZERO || m_region[2] <= RZERO )
  {
    return CPM_ERROR_INVALID_REGION;
  }
  if( m_voxNum[0] <= 0 || m_voxNum[1] <=0 || m_voxNum[2] <= 0 )
  {
    return CPM_ERROR_INVALID_VOXELSIZE;
  }

  return CPM_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
// デフォルトコンストラクタ
cpm_ActiveSubdomainInfo::cpm_ActiveSubdomainInfo()
  : cpm_Base()
{
  clear();
}

////////////////////////////////////////////////////////////////////////////////
// コンストラクタ
cpm_ActiveSubdomainInfo::cpm_ActiveSubdomainInfo( int pos[3] )
  : cpm_Base()
{
  SetPos(pos);
}

////////////////////////////////////////////////////////////////////////////////
// デストラクタ
cpm_ActiveSubdomainInfo::~cpm_ActiveSubdomainInfo()
{
}

////////////////////////////////////////////////////////////////////////////////
// 情報のクリア
void
cpm_ActiveSubdomainInfo::clear()
{
  m_pos[0] = 0;
  m_pos[1] = 0;
  m_pos[2] = 0;
}

////////////////////////////////////////////////////////////////////////////////
// 位置情報のセット
void
cpm_ActiveSubdomainInfo::SetPos( int pos[3] )
{
  m_pos[0] = pos[0];
  m_pos[1] = pos[1];
  m_pos[2] = pos[2];
}

////////////////////////////////////////////////////////////////////////////////
// 位置情報の取得
const int*
cpm_ActiveSubdomainInfo::GetPos() const
{
  return m_pos;
}

////////////////////////////////////////////////////////////////////////////////
// 比較演算子
bool
cpm_ActiveSubdomainInfo::operator==( cpm_ActiveSubdomainInfo dom )
{
  if( m_pos[0] != dom.m_pos[0] ) return false;
  if( m_pos[1] != dom.m_pos[1] ) return false;
  if( m_pos[2] != dom.m_pos[2] ) return false;
  return true;
}

////////////////////////////////////////////////////////////////////////////////
// 比較演算子
bool
cpm_ActiveSubdomainInfo::operator!=( cpm_ActiveSubdomainInfo dom )
{
  if( m_pos[0] == dom.m_pos[0] ) return false;
  if( m_pos[1] == dom.m_pos[1] ) return false;
  if( m_pos[2] == dom.m_pos[2] ) return false;
  return true;
}

////////////////////////////////////////////////////////////////////////////////
// コンストラクタ
cpm_GlobalDomainInfo::cpm_GlobalDomainInfo()
  : cpm_DomainInfo()
{
  clear();
}

////////////////////////////////////////////////////////////////////////////////
// デストラクタ
cpm_GlobalDomainInfo::~cpm_GlobalDomainInfo()
{
}

////////////////////////////////////////////////////////////////////////////////
// 情報のクリア
void
cpm_GlobalDomainInfo::clear()
{
  cpm_DomainInfo::clear();
  for( int i=0;i<3;i++ )
  {
    m_divNum[i] = 0;
  }
  m_subDomainInfo.clear();
}

////////////////////////////////////////////////////////////////////////////////
// 領域分割数のセット
void
cpm_GlobalDomainInfo::SetDivNum( int div[3] )
{
  m_divNum[0] = div[0];
  m_divNum[1] = div[1];
  m_divNum[2] = div[2];
}

////////////////////////////////////////////////////////////////////////////////
// 領域分割数の取得
const int*
cpm_GlobalDomainInfo::GetDivNum() const
{
  return m_divNum;
}

////////////////////////////////////////////////////////////////////////////////
// 活性サブドメイン情報の存在チェック
bool
cpm_GlobalDomainInfo::IsExistSubdomain( cpm_ActiveSubdomainInfo subDomain )
{
  for( size_t i=0;i<m_subDomainInfo.size();i++ )
  {
    cpm_ActiveSubdomainInfo dom = m_subDomainInfo[i];
    if( dom == subDomain ) return true;
  }
  return false;
}

////////////////////////////////////////////////////////////////////////////////
// 活性サブドメイン情報の追加
bool
cpm_GlobalDomainInfo::AddSubdomain( cpm_ActiveSubdomainInfo subDomain )
{
  //既存チェック
  if( IsExistSubdomain(subDomain) ) return false;

  //追加
  m_subDomainInfo.push_back(subDomain);
  return true;
}

////////////////////////////////////////////////////////////////////////////////
// 活性サブドメインの数を取得
int
cpm_GlobalDomainInfo::GetSubdomainNum() const
{
  if( m_subDomainInfo.size() > 0 )
  {
    return (int)m_subDomainInfo.size();
  }
  if( m_divNum[0] <= 0 || m_divNum[1] <= 0 || m_divNum[2] <= 0 )
  {
    return 0;
  }
  return m_divNum[0] * m_divNum[1] * m_divNum[2];
}

////////////////////////////////////////////////////////////////////////////////
// 活性サブドメインの数を取得(情報数)
int
cpm_GlobalDomainInfo::GetSubdomainArraySize() const
{
  return (int)m_subDomainInfo.size();
}

////////////////////////////////////////////////////////////////////////////////
// 活性サブドメイン情報を取得
const cpm_ActiveSubdomainInfo*
cpm_GlobalDomainInfo::GetSubdomainInfo( size_t idx ) const
{
  if( int(idx) >= GetSubdomainNum() ) return NULL;
  return &(m_subDomainInfo[idx]);
}

////////////////////////////////////////////////////////////////////////////////
// 領域情報のチェック
cpm_ErrorCode cpm_GlobalDomainInfo::CheckData( int nRank )
{
  cpm_ErrorCode ret;

  // 親クラス
  if( (ret = cpm_DomainInfo::CheckData()) != CPM_SUCCESS )
  {
    return ret;
  }

  // 領域分割数
  if( m_divNum[0] <= 0 || m_divNum[1] <= 0 || m_divNum[2] <= 0 )
  {
    return CPM_ERROR_INVALID_DIVNUM;
  }

  // 活性サブドメイン情報
  int ndom = m_subDomainInfo.size();
  if( ndom == 0 )
  {
    //活性サブドメイン情報が空のとき、全領域を活性サブドメインとする
    if( nRank != m_divNum[0]*m_divNum[1]*m_divNum[2] )
    {
      return CPM_ERROR_MISMATCH_NP_SUBDOMAIN;
    }
    for( int k=0;k<m_divNum[2];k++ ){
    for( int j=0;j<m_divNum[1];j++ ){
    for( int i=0;i<m_divNum[0];i++ ){
      int pos[3] = {i,j,k};
      cpm_ActiveSubdomainInfo dom( pos );
      AddSubdomain( dom );
    }}}
  }
  else
  {
    if( nRank != ndom )
    {
      return CPM_ERROR_MISMATCH_NP_SUBDOMAIN;
    }
  }

  return CPM_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
// コンストラクタ
cpm_LocalDomainInfo::cpm_LocalDomainInfo()
  : cpm_DomainInfo(), cpm_ActiveSubdomainInfo()
{
}

////////////////////////////////////////////////////////////////////////////////
// デストラクタ
cpm_LocalDomainInfo::~cpm_LocalDomainInfo()
{
}

////////////////////////////////////////////////////////////////////////////////
// 情報のクリア
void
cpm_LocalDomainInfo::clear()
{
  cpm_DomainInfo::clear();
  cpm_ActiveSubdomainInfo::clear();
}

