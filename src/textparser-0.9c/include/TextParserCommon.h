/****************************************************************************
**
** Copyright (C) 2012 Tokyo University.
**
****************************************************************************/
/** @file TextParserCommon.h
 * ここにはシステムで共通のパラメータが定義されています。
 *
 */

#ifndef __TEXTPARSER_COMMON_H__
#define __TEXTPARSER_COMMON_H__

#define	TP_DEBUG	1
#define TP_INDENT_STRING    "    "

/** エレメントのタイプ
 *
 */
typedef enum {
    TP_UNKNOWN_ELEMENT = 0,   //!< 不明
    TP_NODE_ELEMENT = 1, //!< ノード
    TP_LEAF_ELEMENT = 2,      //!< リーフ
    TP_VALUE_ELEMENT = 3,     //!< 値
} TextParserElementType;

/** 値のタイプ
 *
 */
typedef enum {
    TP_UNDEFFINED_VALUE = 0,  //!< 不定
    TP_NUMERIC_VALUE = 1,     //!< 数値
    TP_STRING_VALUE = 2,      //!< 文字列
    TP_DEPENDENCE_VALUE = 3,  //!< 依存関係付き値
    TP_VECTOR_UNDEFFINED = 4, //!< ベクトル型不定
    TP_VECTOR_NUMERIC = 5,    //!< ベクトル型数値
    TP_VECTOR_STRING = 6,     //!< ベクトル型文字列
} TextParserValueType;

/** 論理値
 *
 */
typedef enum {
    TP_UNDEFINED_BOOL = 0,    //!< 不定
    TP_TRUE = 1,              //!< True
    TP_FALSE = 2,             //!< False
} TextParserBool;

/** パーサーの解析モード
 *
 */
typedef enum {
    TP_NO_PARSE = 0,          //!< 解析中ではない
    TP_NODE_PARSE = 1,   //!< ノード解析
    TP_LEAF_PARSE = 2,        //!< リーフ解析
} TextParserParseMode;

/** エラーコード
 *
 */
typedef enum {
    TP_NO_ERROR = 0,                          //!< エラーなし
    TP_ERROR = 100,                           //!< エラー
    TP_DATABASE_NOT_READY_ERROR = 101,        //!< データベースにアクセス出来ない
    TP_DATABASE_ALREADY_SET_ERROR = 102,      //!< データベースが既に読み込まれている
    TP_FILEOPEN_ERROR = 103,                  //!< ファイルオープンエラー
    TP_FILEINPUT_ERROR = 104,                 //!< ファイル入力エラー
    TP_FILEOUTPUT_ERROR = 105,                //!< ファイル出力エラー
    TP_ENDOF_FILE_ERROR = 106,                //!< ファイルの終わりに達しました
    TP_ILLEGAL_TOKEN_ERROR = 107,             //!< トークンが正しくない
    TP_MISSING_LABEL_ERROR = 108,             //!< ラベルが見つからない
    TP_ILLEGAL_LABEL_ERROR = 109,             //!< ラベルが正しくない
    TP_ILLEGAL_ARRAY_LABEL_ERROR = 110,       //!< 配列型ラベルが正しくない
    TP_MISSING_ELEMENT_ERROR = 111,           //!< エレメントが見つからない
    TP_ILLEGAL_ELEMENT_ERROR = 112,           //!< エレメントが正しくない
    TP_NODE_END_ERROR = 113,             //!< ノードの終了文字が多い
    TP_NODE_END_MISSING_ERROR = 114,     //!< ノードの終了文字が無い
    TP_NODE_NOT_FOUND_ERROR = 115,       //!< ノードが見つからない
    TP_LABEL_ALREADY_USED_ERROR = 116,        //!< ラベルが既に使用されている
    TP_LABEL_ALREADY_USED_PATH_ERROR = 117,   //!< ラベルがパス内で既に使用されている
    TP_ILLEGAL_CURRENT_ELEMENT_ERROR = 118,   //!< カレントのエレメントが異常
    TP_ILLEGAL_PATH_ELEMENT_ERROR = 119,      //!< パスのエレメントが異常
    TP_MISSING_PATH_ELEMENT_ERROR = 120,      //!< パスのエレメントが見つからない
    TP_ILLEGAL_LABEL_PATH_ERROR = 121,        //!< パスのラベルが正しくない
    TP_UNKNOWN_ELEMENT_ERROR = 122,           //!< 不明のエレメント
    TP_MISSING_EQUAL_NOT_EQUAL_ERROR = 123,   //!< ==も!=も見つからない
    TP_MISSING_AND_OR_ERROR = 124,            //!< &&も||も見つからない
    TP_MISSING_CONDITION_EXPRESSION_ERROR = 125, //!< 条件式が見つからない
    TP_MISSING_CLOSED_BRANCKET_ERROR = 126,   //!< 条件式が見つからない
    TP_ILLEGAL_CONDITION_EXPRESSION_ERROR = 127,  //!< 条件式の記述が正しくない
    TP_ILLEGAL_DEPENDENCE_EXPRESSION_ERROR = 128, //!< 依存関係の記述が正しくない
    TP_MISSING_VALUE_ERROR = 129,             //!< 値が見つからない
    TP_ILLEGAL_VALUE_ERROR = 130,             //!< 値が正しくない
    TP_ILLEGAL_NUMERIC_VALUE_ERROR = 131,     //!< 数値が正しくない
    TP_ILLEGAL_VALUE_TYPE_ERROR = 132,        //!< ベクトルの値タイプが一致しない
    TP_MISSING_VECTOR_END_ERROR = 133,        //!< ベクトルの終了文字が無い
    TP_VALUE_CONVERSION_ERROR = 134,          //!< 値の変換エラー
    TP_MEMORY_ALLOCATION_ERROR = 135,         //!< メモリが確保できない
    TP_REMOVE_ELEMENT_ERROR = 136,            //!< エレメントの削除エラー
    TP_MISSING_COMMENT_END_ERROR = 137,       //!< コメントの終わりが見つからない
    TP_ID_OVER_ELEMENT_NUMBER_ERROR = 138,    //!< IDが要素数を超えている
    TP_GET_PARAMETER_ERROR = 139,             //!< パラメータ取得
    TP_UNSUPPORTED_ERROR = 199,               //!< サポートされていない
    TP_WARNING = 200,                         //!< エラー
    TP_UNDEFINED_VALUE_USED_WARNING = 201,    //!< 未定義のデータが使われている
    TP_UNRESOLVED_LABEL_USED_WARNING = 202,   //!< 未解決のラベルが使われている
} TextParserError;

#endif // __TEXTPARSER_COMMON_H__
