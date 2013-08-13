#ifndef _FB_PARSE_M_H_
#define _FB_PARSE_M_H_

//##################################################################################
//
// Flow Base class
//
// Copyright (c) 2007-2011 VCAD System Research Program, RIKEN.
// All rights reserved.
//
// Copyright (c) 2011-2013 Institute of Industrial Science, The University of Tokyo.
// All rights reserved.
//
// Copyright (c) 2012-2013 Advanced Institute for Computational Science, RIKEN.
// All rights reserved.
//
//##################################################################################

/**
 * @file   ParseMat.h
 * @brief  FlowBase ParseMat class Header
 * @author kero
 */

#include "FB_Define.h"
#include "Medium.h"
#include "Component.h"
#include "FBUtility.h"
#include "TextParser.h"

class  ParseMat {
private:
  bool ChkList[property_END];  // # of parameters in MediumList must be less than # of property_END
  
  int NoMedium;    ///< コンポーネントの数
  
  TextParser* tpCntl;
  
  
public:
  
  /** コンストラクタ */
  ParseMat() {
    NoMedium   = 0;
    tpCntl     = NULL;
    
    for (int i=0; i<property_END; i++) ChkList[i]=false;
  }
  
  /**　デストラクタ */
  ~ParseMat() {}
  
private:
  
  // mat[]に値を格納する
  void addVaues(MediumList* mat, const int m, const property_list key, const REAL_TYPE fval);
  
  // ラベルの重複を調べる
  bool chkDuplicateLabel(const MediumList* mat, const int n, const std::string m_label);
  
  
  // 媒質情報の内容物をチェックする
  bool chkList4Solver(const MediumList* mat, const int m);
  
  
  // 警告メッセージの表示
  int missingMessage(const MediumList* mat, const int m, const int key);
  
  
public:

  /**
   * @brief MediumListのチェック
   * @param [in] mat      MediumList
   */
  bool check(const MediumList* mat);
  
  
  /**
   * @brief TPのポインタを受け取る
   * @param [in] tp  TextParser
   */
  bool importTP(TextParser* tp);
  
  
  /**
   * @brief MediumTableをロードしてmat[]に保持
   * @param [in]     m_NoMedium  媒質数
   * @param [in,out] mat         MediumList
   */
  void getMediumTable(const int m_NoMedium, MediumList* mat);
  
  
  /**
   * @brief 取得したCompoList[]の内容を表示する
   * @param [in] fp      ファイルポインタ
   * @param [in] compo   CompoList
   * @param [in] basicEq 基礎方程式の種類
   * @note Hostonly
   */
  void chkList(FILE* fp, const CompoList* compo, const int basicEq);
  
 
  /**
   * @brief 媒質情報の表示
   * @param [in] fp       ファイルポインタ
   * @param [in] mat      MediumList
   * @note Hostonly
   */
  void printMatList(FILE* fp, const MediumList* mat);

};

#endif // _FB_PARSE_M_H_
