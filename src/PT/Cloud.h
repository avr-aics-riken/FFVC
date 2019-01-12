#ifndef _PT_CLOUD_H_
#define _PT_CLOUD_H_

//##################################################################################
//
// Copyright (c) 2019 Research Institute for Information Technology, Kyushu university
// All rights researved.
//
//##################################################################################

/**
 * @file   Cloud.h
 * @brief  Cloud class Header
 * @author riit
 * @note
 */

#include "DomainInfo.h"
#include "Chunk.h"
#include <vector>
#include <random>
#include "IntervalManager.h"
#include "Tracking.h"
#include "mydebug.h"
#include "PtDefine.h"
#include "PtComm.h"
#include "TextParser.h"
#include "EmitGroup.h"

using namespace Vec3class;
using std::vector;
using std::string;


/**
 * @brief 粒子群管理クラス
 * @note
 */
class Cloud : public DomainInfo {
protected:
  EmitGroup* Egrp;           ///< 粒子の開始点グループ
  vector<Chunk*> chunkList;  ///< チャンクリスト

  unsigned nParticle;        ///< 全粒子の数（ローカル）
  unsigned gParticle;        ///< 全粒子の数（グローバル）
  int nGrpEmit;              ///< 開始点のグループ数 set[@]の数
  int scheme;                ///< 積分方法
  bool flag_migration;       ///< マイグレーション発生 true
  REAL_TYPE dt;              ///< 時間積分幅
  unsigned buf_max_len;      ///< 送信用のバッファ長さの最大値 毎回異なる
  bool buf_updated;          ///< バッファ長さが更新されたときtrue

  int unit;                  ///< 指定座標の記述単位 {DIMENSIONAL | NONDIMENSIONAL}
  REAL_TYPE refLen;          ///< 代表長さ

  int* bcd;                  ///< BCindex B
  REAL_TYPE* vSource;        ///< 速度サンプリング元データ

  IntervalManager* Interval; ///< タイミング制御
  TextParser* tpCntl;        ///< TextParser
  Tracking* tr;              ///< Tracking

  PtComm PC;                 ///< 粒子通信クラス
  MPI_Request req[NDST*2];   ///< 非同期通信同期


public:
  /// デフォルトコンストラクタ
  Cloud() {}

  Cloud(int* m_bcd,
        REAL_TYPE* m_Vsrc,
        const REAL_TYPE m_refLen,
        const int m_unit,
        const REAL_TYPE dt,
        TextParser* m_tp)
  {
    nParticle = 0;
    gParticle = 0;
    nGrpEmit = 0;
    scheme = -1;
    buf_updated = false;
    flag_migration = false;

    this->dt         = dt;
    this->bcd        = m_bcd;
    this->vSource    = m_Vsrc;
    this->unit       = m_unit;
    this->refLen     = m_refLen;
    this->tpCntl     = m_tp;

    // 初期値として、BUF_UNIT*粒子分を確保 > buf_max_lenで100単位で更新
    buf_max_len = BUF_UNIT;
  }


  /// デストラクタ
  ~Cloud() {
    delete [] Interval;
  }


  // @brief 初期設定
  bool initCloud(FILE* fp);


  // @brief ランタイム
  bool tracking(const unsigned step, const double time);


  /// チャンク数を返す
  unsigned getChunkListSize() const {
    return (unsigned)chunkList.size();
  }


  // @brief マイグレーション処理
  bool migration();


  // @brief ファイル出力
  bool fileout();




protected:

  // @brief パラメータ表示
  void displayParam(FILE* fp);


  // @brief 受信粒子のアンパック
  void unpackParticle();


  // @brief 既に存在するグループIDの検索
  bool searchGrp(const int c);


  // @brief 開始点のユニークIDを割り振る
  bool determineUniqueID();


  /// @brief 粒子追跡情報を取得し，chunkに保持する
  bool setPTinfo();


  // @brief 粒子追跡パラメータ取得
  bool getTPparam(const string label_leaf,
                  const int odr);


  /// @brief 座標値を無次元化する
  /// @param [in,out] x  coordinate
  void normalizeCord(REAL_TYPE x[3])
  {
    x[0] /= refLen;
    x[1] /= refLen;
    x[2] /= refLen;
  }


  /// @brief 自領域内に存在するかどうかを判断
  /// @param [in] x  coordinate
  bool inOwnRegion(const REAL_TYPE x[3])
  {
    if ( x[0]<origin[0] || x[0]>origin[0]+region[0]
      || x[1]<origin[1] || x[1]>origin[1]+region[1]
      || x[2]<origin[2] || x[2]>origin[2]+region[2] ) return false;
    return true;
  }


  /// @brief 開始座標情報を取得し、chunkに保持(Line)
  bool setLine(const string label_base,
               const int odr);


  /// @brief 座標情報を取得し、chunkに保持(PointSet)
  bool setPointset(const string label_base,
                   const int odr);


  /// @brief 開始座標情報を取得し、chunkに保持(Disc)
  bool setDisc(const string label_base,
               const int odr);


  /// @brief メルセンヌツイスタ
  /// @param [in]     st  下限
  /// @param [in]     ed  上限
  /// @url http://vivi.dyndns.org/tech/cpp/random.html
  double mts(const double st, const double ed)
  {
    std::random_device rnd;   // 非決定的な乱数生成器を生成
    std::mt19937 mt(rnd());   // メルセンヌ・ツイスタの32ビット版、引数は初期シード値

    // [st, ed] 範囲の一様乱数
    std::uniform_real_distribution<> rand_s(st, ed);

    return rand_s(mt);
  }


  /// @brief 半径r内のサンプリング
  void samplingInCircle(const REAL_TYPE* cnt,
                        const REAL_TYPE* nv,
                        const REAL_TYPE radius,
                        const int nSample,
                        const int odr);


  /// @brief 指定法線nvがz軸の方向ベクトルに向かう回転角を計算する
  /// @param [in]  nv  指定法線
  /// @retval  回転角
  Vec3r getAngle(Vec3r nv);


  // @brief 回転ベクトルp(alpha, beta, gamma)でベクトルuを回転する
  Vec3r rotate(const Vec3r p, const Vec3r u);


  // @brief 回転ベクトルp(alpha, beta, gamma)に対して，-pでベクトルuを回転する
  Vec3r rotate_inv(const Vec3r p, const Vec3r u);
};

#endif // _PT_CLOUD_H_
