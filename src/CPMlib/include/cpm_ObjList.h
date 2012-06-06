/*
 * CPM - Cartesian Partition Manager
 *
 * Copyright (C) 2012 University of Tokyo.
 *
 */

/**
 * @file   cpm_ObjList.h
 * 汎用オブジェクトの管理クラスのヘッダーファイル
 * @author University of Tokyo
 * @date   2012/05/31
 */

#ifndef _CPM_OBJLIST_H_
#define _CPM_OBJLIST_H_

#include <map>
#include <list>
#include "cpm_Base.h"


/** プロセスグループ毎のランク番号マップ */
typedef std::map<int, int*> RankNoMap;

/** CPMの汎用オブジェクト管理クラス
 */
template<class T>
class cpm_ObjList : public cpm_Base
{
////////////////////////////////////////////////////////////////////////////////
// メンバー変数
////////////////////////////////////////////////////////////////////////////////
public:


private:
  /** オブジェクトのマップ */
  typedef std::map<int, void*> ObjectMap;
  ObjectMap m_ObjectMap;

  /** 削除済み登録番号のリスト */
  typedef std::list<int> DelKeyList;
  DelKeyList m_DelKeyList;

  /** 使用可能な登録番号 */
  int m_newKey;

  
////////////////////////////////////////////////////////////////////////////////
// メンバー関数
////////////////////////////////////////////////////////////////////////////////
public:

  /** コンストラクタ */
  cpm_ObjList()
  {
    m_newKey = 1;
    m_ObjectMap.clear();
    m_DelKeyList.clear();
  };

  /** デストラクタ */
  ~cpm_ObjList()
  {
    ObjectMap::iterator it = m_ObjectMap.begin();
    for( ; it != m_ObjectMap.end(); it++ )
    {
      if( it->second ) delete (T*)it->second;
    }
    m_ObjectMap.clear();
    m_DelKeyList.clear();
  };

  /** オブジェクトの生成
   *  デフォルトコンストラクタが必要
   *  @return 生成したオブジェクトのポインタ
   */
  T* Create()
  {
    T* obj = new T;
    return obj;
  };

  /** オブジェクトの追加
   *  @param[in] obj 追加するオブジェクト
   *  @return 登録番号(負のとき登録失敗)
   */
  int Add( T *obj )
  {
    int key = -1;
    if( m_DelKeyList.size() == 0 )
    {
      if( !(m_ObjectMap.insert(std::make_pair(m_newKey, obj)).second) )
      {
        return -1;
      }
      key = m_newKey;
      m_newKey++;
    }
    else
    {
      DelKeyList::iterator it = m_DelKeyList.begin();
      if( !(m_ObjectMap.insert(std::make_pair(*it, obj)).second) )
      {
        return -1;
      }
      key = *it;
      m_DelKeyList.erase(m_DelKeyList.begin());
    }

    return key;
  };

  /** オブジェクトの削除
   *  @param[in] key Addの戻り値である登録番号
   *  @return CPM終了コード(0,CPM_SUCCESS=正常終了)
   */
  cpm_ErrorCode Delete( int key )
  {
    T *obj = Get(key);
    if( !obj )
    {
      return CPM_ERROR_INVALID_OBJKEY;
    }

    delete obj;
    m_ObjectMap.erase(key);
    m_DelKeyList.push_back(key);

    return CPM_SUCCESS;
  };

  /** オブジェクトの取得
   *  @param[in] key Addの戻り値である登録番号
   *  @return オブジェクトのポインタ
   */
  T* Get( int key )
  {
    ObjectMap::iterator it = m_ObjectMap.find(key);
    if( it == m_ObjectMap.end() )
    {
      return NULL;
    }

    return (T*)it->second;
  };


private:




};

#endif /* _CPM_OBJLIST_H_ */
