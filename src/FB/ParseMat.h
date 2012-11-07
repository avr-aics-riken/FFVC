#ifndef _FB_PARSE_M_H_
#define _FB_PARSE_M_H_

// #################################################################
//
// CAERU Library
//
// Copyright (c) 2012 All right reserved.
//
// Institute of Industrial Science, University of Tokyo, Japan. 
//
// #################################################################

/**
 * @file   ParseMat.h
 * @brief  FlowBase ParseMat class Header
 * @author kero
 */

#include "FB_Define.h"
#include "Medium.h"
#include "Component.h"
#include "TPControl.h"
#include "FBUtility.h"

class  ParseMat {
private:
  bool ChkList[property_END];  // # of parameters in MediumList must be less than # of property_END
  
  int NoCompo;   ///< コンポーネントの数
  int NoBC;      ///< 境界条件の数
  int Unit_Temp; ///< 温度単位
  int KOS;       ///< 支配方程式の種類
  
  MediumTableInfo *MTITP;
  TPControl* tpCntl;
  
public:
  
  /** コンストラクタ */
  ParseMat() {
    NoCompo    = 0;
    NoBC       = 0;
    Unit_Temp  = 0;
    KOS        = 0;
    MTITP      = NULL;
    tpCntl     = NULL;
    
    for (int i=0; i<property_END; i++) ChkList[i]=false;
  }
  
  /**　デストラクタ */
  ~ParseMat() 
  {
    if ( MTITP ) delete [] MTITP;
  }
  
private:

  /**
   * @brief ラベルの重複を調べる
   * @param [in] mat     MediumList
   * @param [in] n       現在までに登録した媒質リストの格納番号
   * @param [in] m_label 検査するラベル名
   */
  bool chkDuplicateLabel(MediumList* mat, const int n, const std::string m_label);
  
  /**
   * @brief 媒質情報の内容物をチェックする
   * @param [in] mat MediumList
   * @param [in] m   媒質リストの格納番号
   */
  bool chkList4Solver(MediumList* mat, const int m);
  
  
  /**
   * @brief 警告メッセージの表示
   * @param [in] mat MediumList
   * @param [in] m   媒質リストの格納番号
   * @param [in] key キーワードの登録番号
   */
  int missingMessage(MediumList* mat, const int m, const int key);
  
  
  /**
   * @brief matの変数値を格納する
   * @param [in] mat MediumList
   * @param [in] n   媒質リストの格納番号
   */
  void copyProperty(MediumList* mat, const int n);
    
public:
  /**
   * @brief MTITPのポインタをエキスポート
   */
  MediumTableInfo* export_MTI() 
  {
    return MTITP;
  }
  

  /**
   * @brief MediumListを作成する
   * @param [in] mat      MediumList
   * @param [in] NoMedium 媒質リストの数
   */
  bool makeMediumList(MediumList* mat, const int NoMedium);
  
  
  /**
   * @brief TPのポインタを受け取る
   * @param [in] tp  TPControl
   */
  bool importTP(TPControl* tp);
  
  
  /**
   * @brief Medium_Tableを読んでMediumTableInfoクラスオブジェクトに格納
   */
  int get_MediumTable();
  
  
  /**
   * @brief 取得したCompoList[]の内容を表示する
   * @param [in] fp      ファイルポインタ
   * @param [in] compo   CompoList
   * @param [in] basicEq 基礎方程式の種類
   * @note Hostonly
   */
  void chkList(FILE* fp, CompoList* compo, const int basicEq);
  
 
  /**
   * @brief 媒質情報の表示
   * @param [in] fp       ファイルポインタ
   * @param [in] mat      MediumList
   * @param [in] NoMedium 媒質リストの数
   * @note Hostonly
   */
  void printMatList(FILE* fp, MediumList* mat, const int NoMedium);
  
  
  /**
   * @brief CompoList[]とMediumList[]のチェック
   * @param [in] fp    ファイルポインタ
   * @param [in] compo CompoList
   * @param [in] mat   MediumList
   * @note debug function
   */
  void printRelation(FILE* fp, CompoList* compo, MediumList* mat);
  
  /**
   * @brief CompoList[]とMediumList[]のチェック
   * @param [in] m_NoCompo         コンポーネントの数
   * @param [in] m_NoBC CompoList  境界条件の数
   * @param [in] m_Unit_Temp       温度単位
   * @param [in] m_KOS             支配方程式の種類
   */
  void setControlVars(const int m_NoCompo, const int m_NoBC, const int m_Unit_Temp, const int m_KOS);

};

#endif // _FB_PARSE_M_H_
