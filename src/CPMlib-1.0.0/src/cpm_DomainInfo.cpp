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
  for( int i=0;i<3;i++ )
  {
    m_origin[i] = 0.0;
    m_pitch[i]  = 0.0;
    m_region[i] = 0.0;
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
// デフォルトコンストラクタ
cpm_ActiveSubDomainInfo::cpm_ActiveSubDomainInfo()
  : cpm_Base()
{
  clear();
}

////////////////////////////////////////////////////////////////////////////////
// コンストラクタ
cpm_ActiveSubDomainInfo::cpm_ActiveSubDomainInfo(int pos[3], int bcid[6])
  : cpm_Base()
{
  Set(pos, bcid);
}

////////////////////////////////////////////////////////////////////////////////
// デストラクタ
cpm_ActiveSubDomainInfo::~cpm_ActiveSubDomainInfo()
{
}

////////////////////////////////////////////////////////////////////////////////
// 情報のクリア
void
cpm_ActiveSubDomainInfo::clear()
{
  m_pos[0] = 0;
  m_pos[1] = 0;
  m_pos[2] = 0;
  m_bcid[X_MINUS] = 0;
  m_bcid[Y_MINUS] = 0;
  m_bcid[Z_MINUS] = 0;
  m_bcid[X_PLUS]  = 0;
  m_bcid[Y_PLUS]  = 0;
  m_bcid[Z_PLUS]  = 0;
}

////////////////////////////////////////////////////////////////////////////////
// 値のセット
void
cpm_ActiveSubDomainInfo::Set( int pos[3], int bcid[6] )
{
  SetPos(pos);
  SetBCID(bcid);
}

////////////////////////////////////////////////////////////////////////////////
// 位置情報のセット
void
cpm_ActiveSubDomainInfo::SetPos( int pos[3] )
{
  m_pos[0] = pos[0];
  m_pos[1] = pos[1];
  m_pos[2] = pos[2];
}

////////////////////////////////////////////////////////////////////////////////
// 位置情報の取得
const int*
cpm_ActiveSubDomainInfo::GetPos() const
{
  return m_pos;
}

////////////////////////////////////////////////////////////////////////////////
// 境界条件IDのセット
void
cpm_ActiveSubDomainInfo::SetBCID( int bcid[6] )
{
  m_bcid[X_MINUS] = bcid[X_MINUS];
  m_bcid[Y_MINUS] = bcid[Y_MINUS];
  m_bcid[Z_MINUS] = bcid[Z_MINUS];
  m_bcid[X_PLUS]  = bcid[X_PLUS];
  m_bcid[Y_PLUS]  = bcid[Y_PLUS];
  m_bcid[Z_PLUS]  = bcid[Z_PLUS];
}

////////////////////////////////////////////////////////////////////////////////
// 境界条件IDの取得
const int*
cpm_ActiveSubDomainInfo::GetBCID() const
{
  return m_bcid;
}

////////////////////////////////////////////////////////////////////////////////
// 比較演算子
bool
cpm_ActiveSubDomainInfo::operator==( cpm_ActiveSubDomainInfo dom )
{
  if( m_pos[0] != dom.m_pos[0] ) return false;
  if( m_pos[1] != dom.m_pos[1] ) return false;
  if( m_pos[2] != dom.m_pos[2] ) return false;
  return true;
}

////////////////////////////////////////////////////////////////////////////////
// 比較演算子
bool
cpm_ActiveSubDomainInfo::operator!=( cpm_ActiveSubDomainInfo dom )
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
cpm_GlobalDomainInfo::IsExistSubDomain( cpm_ActiveSubDomainInfo subDomain )
{
  for( size_t i=0;i<m_subDomainInfo.size();i++ )
  {
    cpm_ActiveSubDomainInfo dom = m_subDomainInfo[i];
    if( dom == subDomain ) return true;
  }
  return false;
}

////////////////////////////////////////////////////////////////////////////////
// 活性サブドメイン情報の追加
bool
cpm_GlobalDomainInfo::AddSubDomain( cpm_ActiveSubDomainInfo subDomain )
{
  //既存チェック
  if( IsExistSubDomain(subDomain) ) return false;

  //追加
  m_subDomainInfo.push_back(subDomain);
  return true;
}

////////////////////////////////////////////////////////////////////////////////
// 活性サブドメインの数を取得
int
cpm_GlobalDomainInfo::GetSubDomainNum() const
{
  return (int)m_subDomainInfo.size();
}

////////////////////////////////////////////////////////////////////////////////
// 活性サブドメイン情報を取得
const cpm_ActiveSubDomainInfo*
cpm_GlobalDomainInfo::GetSubDomainInfo( size_t idx ) const
{
  if( int(idx) >= GetSubDomainNum() ) return NULL;
  return &(m_subDomainInfo[idx]);
}

////////////////////////////////////////////////////////////////////////////////
// コンストラクタ
cpm_LocalDomainInfo::cpm_LocalDomainInfo()
  : cpm_DomainInfo(), cpm_ActiveSubDomainInfo()
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
  cpm_ActiveSubDomainInfo::clear();
}

