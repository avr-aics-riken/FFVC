/*
 * CPM - Cartesian Partition Manager
 *
 * Copyright (C) 2012 University of Tokyo.
 *
 */

/**
 * @file   cpm_DomainInfo.h
 * 領域情報クラスのヘッダーファイル
 * @author University of Tokyo
 * @date   2012/05/31
 */

#ifndef _CPM_DomainInfo_H_
#define _CPM_DomainInfo_H_

#include <vector>
#include "cpm_Base.h"

/** CPMの領域情報クラス */
class cpm_DomainInfo : public cpm_Base
{
////////////////////////////////////////////////////////////////////////////////
// メンバー関数
////////////////////////////////////////////////////////////////////////////////
public:
  /** コンストラクタ */
  cpm_DomainInfo();

  /** デストラクタ */
  virtual ~cpm_DomainInfo();

  /** 情報のクリア */
  virtual void clear();

  /** 原点のセット
   *  @param[in] org 原点情報
   */
  void SetOrigin( REAL_TYPE org[3] ); 

  /** 原点の取得
   *  @return 原点情報実数配列のポインタ
   */
  const REAL_TYPE* GetOrigin() const;

  /** ピッチのセット
   *  @param[in] pch ピッチ情報
   */
  void SetPitch( REAL_TYPE pch[3] ); 

  /** ピッチの取得
   *  @return ピッチ情報実数配列のポインタ
   */
  const REAL_TYPE* GetPitch() const;

  /** 空間サイズのセット
   *  @param[in] rgn 空間サイズ情報
   */
  void SetRegion( REAL_TYPE rgn[3] ); 

  /** 空間サイズの取得
   *  @return 空間サイズ情報実数配列のポインタ
   */
  const REAL_TYPE* GetRegion() const;

  /** VOXEL数のセット
   *  @param[in] vox VOXEL数情報
   */
  void SetVoxNum( int vox[3] ); 

  /** VOXEL数の取得
   *  @return VOXEL数情報実数配列のポインタ
   */
  const int* GetVoxNum() const;


private:


////////////////////////////////////////////////////////////////////////////////
// メンバー変数
////////////////////////////////////////////////////////////////////////////////
public:


private:
  REAL_TYPE m_origin[3]; ///< 原点
  REAL_TYPE m_pitch[3];  ///< ピッチ
  REAL_TYPE m_region[3]; ///< 空間サイズ
  int       m_voxNum[3]; ///< VOXEL数


};

/** CPMのサブ領域情報クラス */
class cpm_ActiveSubDomainInfo : public cpm_Base
{
////////////////////////////////////////////////////////////////////////////////
// メンバー関数
////////////////////////////////////////////////////////////////////////////////
public:
  /** デフォルトコンストラクタ */
  cpm_ActiveSubDomainInfo();

  /** コンストラクタ
   *  @param[in] pos  領域分割内での位置
   *  @param[in] bcid 面の境界条件ID
   */
  cpm_ActiveSubDomainInfo( int pos[3], int bcid[6] );

  /** デストラクタ */
  virtual ~cpm_ActiveSubDomainInfo();

  /** 情報のクリア */
  virtual void clear();

  /** 値のセット
   *  @param[in] pos  領域分割内での位置
   *  @param[in] bcid 面の境界条件ID
   */
  void Set( int pos[3], int bcid[6] );

  /** 位置のセット
   *  @param[in] pos  領域分割内での位置
   */
  void SetPos( int pos[3] );

  /** 位置の取得
   *  @return 位置情報整数配列のポインタ
   */
  const int* GetPos() const;

  /** 境界条件IDのセット
   *  @param[in] bcid 面の境界条件ID
   */
  void SetBCID( int bcid[6] );

  /** BCIDの取得
   *  @return 境界条件ID情報整数配列のポインタ
   */
  const int* GetBCID() const;

  /** 比較演算子
   *  @param[in] dom   比較対象の活性サブドメイン情報
   *  @retval    true  同じ位置情報を持つ
   *  @retval    false 違う位置情報を持つ
   */
  bool operator==( cpm_ActiveSubDomainInfo dom );

  /** 比較演算子
   *  @param[in] dom   比較対象の活性サブドメイン情報
   *  @retval    true  違う位置情報を持つ
   *  @retval    false 同じ位置情報を持つ
   */
  bool operator!=( cpm_ActiveSubDomainInfo dom );


private:


////////////////////////////////////////////////////////////////////////////////
// メンバー変数
////////////////////////////////////////////////////////////////////////////////
public:


private:
  int m_pos[3];  ///< 領域分割内での位置
  int m_bcid[6]; ///< 6面の境界条件ID


};

/** CPMの全体領域情報クラス */
class cpm_GlobalDomainInfo : public cpm_DomainInfo
{
////////////////////////////////////////////////////////////////////////////////
// メンバー関数
////////////////////////////////////////////////////////////////////////////////
public:
  /** コンストラクタ */
  cpm_GlobalDomainInfo();

  /** デストラクタ */
  virtual ~cpm_GlobalDomainInfo();

  /** 情報のクリア */
  virtual void clear();

  /** 領域分割数のセット
   *  @param[in] div 領域分割数
   */
  void SetDivNum( int div[3] );

  /** 領域分割数の取得
   *  @return 領域分割数整数配列のポインタ
   */
  const int* GetDivNum() const;

  /** 活性サブドメイン情報の存在チェック
   *  @param[in] subDomain チェックする活性サブドメイン情報 
   *  @retval    true      存在する
   *  @retval    false     存在しない
   */
  bool IsExistSubDomain( cpm_ActiveSubDomainInfo subDomain );

  /** 活性サブドメイン情報の追加
   *  @param[in] subDomain 追加する活性サブドメイン情報 
   *  @retval    true      追加した
   *  @retval    false     追加に失敗(同じ領域分割位置で追加済み)
   */
  bool AddSubDomain( cpm_ActiveSubDomainInfo subDomain );

  /** 活性サブドメインの数を取得
   *  @return 活性サブドメインの数
   */
  int GetSubDomainNum() const;

  /** 活性サブドメイン情報を取得
   *  @param[in]  idx 登録順番号
   *  @return 活性サブドメイン情報ポインタ
   */
  const cpm_ActiveSubDomainInfo* GetSubDomainInfo( size_t idx ) const;


private:


////////////////////////////////////////////////////////////////////////////////
// メンバー変数
////////////////////////////////////////////////////////////////////////////////
public:


private:
  int m_divNum[3]; ///< 領域分割数
  std::vector<cpm_ActiveSubDomainInfo> m_subDomainInfo; ///< 活性サブドメイン情報


};

/** CPMのローカル領域情報クラス */
class cpm_LocalDomainInfo : public cpm_DomainInfo, public cpm_ActiveSubDomainInfo
{
////////////////////////////////////////////////////////////////////////////////
// メンバー関数
////////////////////////////////////////////////////////////////////////////////
public:
  /** コンストラクタ */
  cpm_LocalDomainInfo();

  /** デストラクタ */
  virtual ~cpm_LocalDomainInfo();

  /** 情報のクリア */
  virtual void clear();


private:


////////////////////////////////////////////////////////////////////////////////
// メンバー変数
////////////////////////////////////////////////////////////////////////////////
public:


private:


};

#endif /* _CPM_DomainInfo_H_ */

