/*
 * CPM - Cartesian Partition Manager
 *
 * Copyright (C) 2012 University of Tokyo.
 *
 */

/**
 * @file   cpm_TextParserDomain.cpp
 * CPM領域情報のTextParserクラスのソースファイル
 * @author University of Tokyo
 * @date   2012/05/31
 */
#include "cpm_TextParserDomain.h"

////////////////////////////////////////////////////////////////////////////////
// コンストラクタ
cpm_TextParserDomain::cpm_TextParserDomain()
  : cpm_TextParser()
{
}

////////////////////////////////////////////////////////////////////////////////
// デストラクタ
cpm_TextParserDomain::~cpm_TextParserDomain()
{
}

////////////////////////////////////////////////////////////////////////////////
// 読み込み(静的関数)
cpm_GlobalDomainInfo*
cpm_TextParserDomain::Read( std::string filename, int &errorcode )
{
  // インスタンス
  cpm_TextParserDomain tp;

  // 読み込みメイン
  return tp.ReadMain( filename, errorcode );
}

////////////////////////////////////////////////////////////////////////////////
// 読み込み処理のメイン
cpm_GlobalDomainInfo*
cpm_TextParserDomain::ReadMain( std::string filename, int &errorcode )
{
  errorcode = TP_NO_ERROR;

  // デフォルトの読み込み処理
  if( (errorcode = cpm_TextParser::Read( filename )) != TP_NO_ERROR )
  {
    return NULL;
  }

  // 領域情報のインスタンス
  cpm_GlobalDomainInfo *dInfo = new cpm_GlobalDomainInfo();
  if( !dInfo )
  {
    errorcode = CPM_ERROR_INVALID_PTR;
    return NULL;
  }

  // DomainInfoの読み込み
  if( (errorcode = ReadDomainInfo( dInfo )) != TP_NO_ERROR )
  {
    delete dInfo;
    return NULL;
  }

  // SubdomainInfoの読み込み
  if( (errorcode = ReadSubdomainInfo( dInfo )) != TP_NO_ERROR )
  {
    delete dInfo;
    return NULL;
  }

  return dInfo;
}

////////////////////////////////////////////////////////////////////////////////
// DomainInfoの読み込み
int
cpm_TextParserDomain::ReadDomainInfo( cpm_GlobalDomainInfo* dInfo )
{
  int ret;

  if( !dInfo )
  {
    return CPM_ERROR_INVALID_PTR;
  }

  // カレントノードを変更
  std::string oldpos;
  if( (ret = m_tp->currentNode( oldpos )) != TP_NO_ERROR )
  {
    return ret;
  }
  if( (ret = m_tp->changeNode( "/DomainInfo" )) != TP_NO_ERROR )
  {
    return ret;
  }

  // リーフのラベルリストを取得
  std::vector<std::string> labels;
  if( (ret = m_tp->getLabels( labels )) != TP_NO_ERROR )
  {
    return ret;
  }

  REAL_TYPE RZERO  = REAL_TYPE(0);
  REAL_TYPE org[3] = {RZERO,RZERO,RZERO}; bool borg = false;
  int       vox[3] = {0,0,0};             bool bvox = false;
  REAL_TYPE pch[3] = {RZERO,RZERO,RZERO}; bool bpch = false;
  REAL_TYPE rgn[3] = {RZERO,RZERO,RZERO}; bool brgn = false;
  int       div[3] = {0,0,0};             bool bdiv = false;

  for( size_t i=0;i<labels.size();i++ )
  {
    std::string label = labels[i];

    // G_origin
    if( !borg && cpm_strCompare( label, "G_origin" ) == 0 )
    {
      if( readVector( label, org, 3 ) != TP_NO_ERROR )
      {
        return CPM_ERROR_TP_INVALID_G_ORG;
      }
      borg = true;
      continue;
    }

    // G_voxel
    if( !bvox && cpm_strCompare( label, "G_voxel" ) == 0 )
    {
      if( readVector( label, vox, 3 ) != TP_NO_ERROR )
      {
        return CPM_ERROR_TP_INVALID_G_VOXEL;
      }
      bvox=true;
      continue;
    }

    // G_pitch
    if( !bpch && cpm_strCompare( label, "G_pitch" ) == 0 )
    {
      if( readVector( label, pch, 3 ) != TP_NO_ERROR )
      {
        return CPM_ERROR_TP_INVALID_G_PITCH;
      }
      bpch = true;
      continue;
    }

    // G_region
    if( !brgn && cpm_strCompare( label, "G_region" ) == 0 )
    {
      if( readVector( label, rgn, 3 ) != TP_NO_ERROR )
      {
        return CPM_ERROR_TP_INVALID_G_PITCH;
      }
      brgn = true;
      continue;
    }

    // G_div
    if( !bdiv && cpm_strCompare( label, "G_div" ) == 0 )
    {
      if( readVector( label, div, 3 ) != TP_NO_ERROR )
      {
        return CPM_ERROR_TP_INVALID_G_DIV;
      }
      bdiv = true;
      continue;
    }
  }

  // G_orgをセット(必須)
  if( borg )
    dInfo->SetOrigin( org );
  else
    return CPM_ERROR_TP_INVALID_G_ORG;

  // G_voxelをセット(オプション)
  if( bvox )
  {
    if( vox[0] <= 0 || vox[1] <= 0 || vox[2] <= 0 )
      return CPM_ERROR_TP_INVALID_G_VOXEL;
  }
  dInfo->SetVoxNum( vox );

  // G_region、G_pitchをセット(G_region優先)
  // regionが必ず決定するように指定されている必要がある
  // 1.G_regionが指定されているとき
  //   G_voxが指定されていればpitchを自動計算
  // 2.G_regionが指定されておらずG_pitchが指定されているとき
  //   G_voxが指定されていればregionを自動計算
  //   G_voxが指定されていないときはエラー
  if( brgn )
  {
    if( rgn[0] <= RZERO || rgn[0] <= RZERO || rgn[0] <= RZERO )
    {
      return CPM_ERROR_TP_INVALID_G_RGN;
    }
    dInfo->SetRegion( rgn );
    if( bvox )
    {
      for( int i=0;i<3;i++ ) pch[i] = rgn[i] / REAL_TYPE(vox[i]);
      dInfo->SetPitch( pch );
    }
  }
  else if( bpch && bvox)
  {
    if( pch[0] <= RZERO || pch[0] <= RZERO || pch[0] <= RZERO )
    {
      return CPM_ERROR_TP_INVALID_G_PITCH;
    }
    for( int i=0;i<3;i++ ) rgn[i] = pch[i] * REAL_TYPE(vox[i]);
    dInfo->SetPitch( pch );
    dInfo->SetRegion( rgn );
  }
  else
  {
    return CPM_ERROR_TP_INVALID_G_RGN;
  }

  // G_divをセット(オプション)
  if( bdiv )
  {
    if( div[0] <= 0 || div[1] <= 0 || div[2] <= 0 )
    return CPM_ERROR_TP_INVALID_G_DIV;
  }
  dInfo->SetDivNum( div );

  // 位置を元に戻す
  if( (ret = m_tp->changeNode( oldpos )) != TP_NO_ERROR )
  {
    return ret;
  }

  return CPM_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
// SubdomainInfoの読み込み
int
cpm_TextParserDomain::ReadSubdomainInfo( cpm_GlobalDomainInfo* dInfo )
{
  int ret;

  if( !dInfo )
  {
    return CPM_ERROR_INVALID_PTR;
  }

  // 領域分割数の取得
  const int *div = dInfo->GetDivNum();
  if( !div )
  {
    return CPM_ERROR_TP_INVALID_G_DIV;
  }

  // カレントノードを変更
  std::string oldpos;
  if( (ret = m_tp->currentNode( oldpos )) != TP_NO_ERROR )
  {
    return ret;
  }
  if( (ret = m_tp->changeNode( "SubdomainInfo" )) != TP_NO_ERROR )
  {
    return ret;
  }

  // リーフのラベルリストを取得
  std::vector<std::string> labels;
  if( (ret = m_tp->getLabels( labels )) != TP_NO_ERROR )
  {
    return ret;
  }

  // 各ラベルをパース
  std::string filename = ""; bool bfname = false;
  for( size_t i=0;i<labels.size();i++ )
  {
    // ラベル
    std::string label = labels[i];

    // filename
    if( !bfname && cpm_strCompare( label, "filename" ) == 0 )
    {
      if( (ret = m_tp->getValue( label, filename )) != TP_NO_ERROR )
      {
        return ret;
      }
      bfname = true;
#ifdef _DEBUG
      std::cout << "ActiveSubdomainInfo filename = " << filename << std::endl;
#endif
    }
  }

  // 位置を元に戻す
  if( (ret = m_tp->changeNode( oldpos )) != TP_NO_ERROR )
  {
    return ret;
  }

  return CPM_SUCCESS;
}

