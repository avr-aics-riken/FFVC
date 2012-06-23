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

  // ActiveSubDomainsの読み込み
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

  REAL_TYPE org[3]; bool borg = false;
  int       vox[3]; bool bvox = false;
  REAL_TYPE pch[3]; bool bpch = false;
  REAL_TYPE rgn[3]; bool brgn = false;
  int       div[3]; bool bdiv = false;

  for( size_t i=0;i<labels.size();i++ )
  {
    std::string label = labels[i];

    // G_org
    if( !borg && cpm_strCompare( label, "G_org" ) == 0 )
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

  // G_orgをセット
  if( borg )
    dInfo->SetOrigin( org );
  else
    return CPM_ERROR_TP_INVALID_G_ORG;

  // G_voxelをセット
  if( bvox ) 
    dInfo->SetVoxNum( vox );
  else
    return CPM_ERROR_TP_INVALID_G_VOXEL;

  // G_pitch, G_regionをセット
  // (G_pitch優先)
  if( bpch )
  {
    for( int i=0;i<3;i++ ) rgn[i] = pch[i] * REAL_TYPE(vox[i]);
    dInfo->SetPitch( pch );
    dInfo->SetRegion( rgn );
  }
  else if( brgn )
  {
    for( int i=0;i<3;i++ ) pch[i] = rgn[i] / REAL_TYPE(vox[i]);
    dInfo->SetPitch( pch );
    dInfo->SetRegion( rgn );
  }
  else
  {
    return CPM_ERROR_TP_INVALID_G_PITCH;
  }

  // G_divをセット
  if( bdiv )
    dInfo->SetDivNum( div );
  else
    return CPM_ERROR_TP_INVALID_G_DIV;

  // 位置を元に戻す
  if( (ret = m_tp->changeNode( oldpos )) != TP_NO_ERROR )
  {
    return ret;
  }

  return CPM_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
// ActiveSubDomainsの読み込み
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
  if( (ret = m_tp->changeNode( "ActiveSubDomains" )) != TP_NO_ERROR )
  {
    return ret;
  }

  // 子ノードのラベルリストを取得
  std::vector<std::string> subNodes;
  if( (ret = m_tp->getNodes( subNodes )) != TP_NO_ERROR )
  {
    return ret;
  }

  // 各ノードをパース
  for( size_t i=0;i<subNodes.size();i++ )
  {
    // ラベル
    std::string node = subNodes[i];

    // domain[xx]をチェック
    if( cpm_strCompareN( node, "domain[", 7 ) != 0 ) continue;

    // 子ノードに移動
    if( (ret = m_tp->changeNode( node )) != TP_NO_ERROR )
    {
      return ret;
    }

    // リーフのラベルリストを取得
    std::vector<std::string> labels;
    if( (ret = m_tp->getLabels( labels )) != TP_NO_ERROR )
    {
      return ret;
    }

    // リーフラベルを検索
    int pos[3];  bool bpos = false;
    int bcid[6]; bool bbcid = false;
    for( size_t j=0;j<labels.size();j++ )
    {
      std::string label = labels[j];

      // pos
      if( !bpos && cpm_strCompare( label, "pos" ) == 0 )
      {
        if( readVector( label, pos, 3 ) != TP_NO_ERROR )
        {
          return CPM_ERROR_TP_INVALID_POS;
        }
        if( pos[0] < 0 || pos[0] >=div[0] ||
            pos[1] < 0 || pos[1] >=div[1] ||
            pos[2] < 0 || pos[2] >=div[2] )
        {
          return CPM_ERROR_TP_INVALID_POS;
        }
        bpos = true;
        continue;
      }

      // bcid
      if( !bbcid && cpm_strCompare( label, "bcid" ) == 0 )
      {
        if( readVector( label, bcid, 6 ) != TP_NO_ERROR )
        {
          return CPM_ERROR_TP_INVALID_BCID;
        }
        bbcid = true;
        continue;
      }
    }

    // サブドメイン情報にセット
    cpm_ActiveSubDomainInfo dom;
    if( bpos )
      dom.SetPos(pos);
    else
      return CPM_ERROR_TP_INVALID_POS;
    if( bbcid )
      dom.SetBCID(bcid);
//    else
//      return CPM_ERROR_TP_INVALID_BCID;
      
    // ドメイン情報に追加
    dInfo->AddSubDomain(dom);

    // 親ノードに戻る
    if( (ret = m_tp->changeNode( ".." )) != TP_NO_ERROR )
    {
      return ret;
    }
  }

  // 位置を元に戻す
  if( (ret = m_tp->changeNode( oldpos )) != TP_NO_ERROR )
  {
    return ret;
  }

  return CPM_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
// 関数名

