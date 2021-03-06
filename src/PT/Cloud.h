#ifndef _PT_CLOUD_H_
#define _PT_CLOUD_H_

//##################################################################################
//
// Copyright (c) 2019-2020 Research Institute for Information Technology, Kyushu university
// All rights researved.
//
//##################################################################################

/*
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
#include <stdlib.h>
#include <errno.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <iostream>
#include <fstream>
//#include <iomanip>


// FX10 profiler
#if defined __K_FPCOLL
#include "fjcoll.h"
#elif defined __FX_FAPP
#include "fj_tool/fjcoll.h"
#endif

#include "PerfMonitor.h"

using namespace pm_lib;
using namespace Vec3class;
using std::vector;
using std::string;


/**
 * @brief 粒子群管理クラス
 * @note
 */
class Cloud : public DomainInfo {
public:
  unsigned l_CurrentParticle; ///< 現在の全粒子の数（ローカル）
  unsigned g_CurrentParticle; ///< 現在の全粒子の数（グローバル）
  unsigned g_DeletedParticle; ///< 削除された全粒子数（グローバル）
  unsigned g_EmittedParticle; ///< 放出された全粒子数（グローバル）
  unsigned l_EmittedParticle; ///< 放出された全粒子数（ローカル）
  unsigned g_OutParticle;     ///< 領域外へ出た粒子数（グローバル）
  unsigned g_WallParticle;    ///< 壁を通過した粒子数（グローバル）>> 非アクティブ
  unsigned g_LimitParticle;   ///< 寿命が尽きた粒子数（グローバル）
  unsigned g_ClippedParticle; ///< クリップされた粒子数（グローバル）
  unsigned g_MigrateParticle; ///< マイグレーション時の送受信粒子数
  unsigned buf_max_particle;  ///< 送信用のバッファ長さ計算に使う粒子の最大数 毎回異なる
  int* Rmap;                  ///< PtCommクラスで作成した3x3x3のランクマップへのポインタ
  
  
protected:
  EmitGroup* Egrp;           ///< 粒子の開始点グループ
  vector<Chunk*> chunkList;  ///< チャンクリスト
  
  int nGrpEmit;              ///< 開始点のグループ数 set[@]の数
  unsigned nEmitParticle;    ///< 粒子放出点の数
  unsigned nEmission;        ///< 粒子放出の回数
  int EmitStart;             ///< 開始時刻
  int EmitInterval;          ///< 粒子放出インターバル
  int OutInterval;           ///< 出力インターバル
  int EmitLife;              ///< 寿命情報 -1 制御なし, 正数<MAX_LIFE未満
  int scheme;                ///< 積分方法
  REAL_TYPE dt;              ///< 時間積分幅
  bool buf_updated;          ///< バッファ長さが更新されたときtrue
  int clip_dir;              ///< クリップ処理の方向
  REAL_TYPE clip_val;        ///< 座標
  
  int unit;                  ///< 指定座標の記述単位 {DIMENSIONAL | NONDIMENSIONAL}
  REAL_TYPE refLen;          ///< 代表長さ
  REAL_TYPE refVel;          ///< 代表速度
  
  int* bcd;                  ///< BCindex B
  int* bid;                  ///< 境界条件ID
  REAL_TYPE* vSource;        ///< 速度サンプリング元データ
  
  IntervalManager Interval;  ///< タイミング制御
  TextParser* tpCntl;        ///< TextParser
  Tracking* tr;              ///< Tracking
  PtComm PC;                 ///< 粒子通信クラス
  PerfMonitor* PM;           ///< PerfMonitor class
  
  int out_format;            ///< 粒子出力ファイルフォーマット　0 - ascii, 1 - binary
  FILE* fpl1;                ///< ログ出力用のファイルポインタ
  FILE* fpl2;                ///< ログ出力用のファイルポインタ
  int restartRankFlag;       ///< リスタート時に粒子データを読むランクのみ > 1, else 0
  int restartFlag;           ///< リスタート指定時に1, else 0
  int restartStep;           ///< リスタートステップ
  int restartForm;           ///< リスタートファイルのフォーマット　0 - ascii, 1 - binary
  
  
  
public:
  /// デフォルトコンストラクタ
  Cloud() {
    Rmap = NULL;
  }
  
  Cloud(int* m_bcd,
        int* m_bid,
        REAL_TYPE* m_Vsrc,
        const REAL_TYPE dt,
        TextParser* m_tp,
        PerfMonitor* m_PM)
  {
    l_CurrentParticle = 0;
    g_CurrentParticle = 0;
    g_DeletedParticle = 0;
    g_OutParticle = 0;
    g_WallParticle = 0;
    g_LimitParticle = 0;
    g_ClippedParticle = 0;
    g_EmittedParticle = 0;
    l_EmittedParticle = 0;
    g_MigrateParticle = 0;
    
    nGrpEmit = 0;
    scheme = -1;
    buf_updated = false;
    out_format = -1;
    unit = NONDIMENSIONAL;
    refLen = 0.0;
    refVel = 0.0;
    nEmitParticle = 0;
    nEmission = 0;
    EmitStart = 0;
    EmitInterval = 0;
    EmitLife = -1;
    restartRankFlag = 0;
    restartFlag = 0;
    restartStep = -1;
    
    OutInterval = 0;
    clip_dir = -1;
    clip_val = 0.0;
    
    
    this->dt         = dt;
    this->bcd        = m_bcd;
    this->bid        = m_bid;
    this->vSource    = m_Vsrc;
    this->tpCntl     = m_tp;
    this->PM         = m_PM;
    
    Interval.setMode(IntervalManager::By_step);
    
    // 初期値として、BUF_UNIT*粒子分を確保 > buf_max_particleで100単位で更新
    buf_max_particle = BUF_UNIT;
  }
  
  
  
  /// デストラクタ
  ~Cloud() {
    if (fpl1) fclose(fpl1);
    if (fpl2) fclose(fpl2);
  }
  
  
  // ######################
  
  // @brief 初期設定
  bool initCloud(FILE* fp);
  
  
  // @brief ランタイム
  bool tracking(const unsigned step, const double time);
  
  
  
  
  // ######################
protected:
  
  // @brief chunkを登録
  // @param [in]  pos    座標ベクトル
  void registChunk(Vec3r pos);
  
  
  // @brieaf 粒子のchunkListへの追加
  void addParticle2ChunkList(particle p);
  
  
  // chunkListの登録粒子数を返す
  unsigned getNparticle() {
    unsigned tmp = 0;
    for(auto itr = chunkList.begin(); itr != chunkList.end(); ++itr) {
      tmp += (unsigned)(*itr)->getNpoints();
    }
    return tmp;
  }
  
  
  // @brieaf リスタート時の粒子データ読み込み(Binary)
  bool readRestartParticleBinary();
  
  
  // @brieaf リスタート時の粒子データ読み込み(Ascii)
  bool readRestartParticleAscii();
  
  
  // @brieaf リスタート時のランク番号情報
  bool readRestartRank();
  
  
  // @brief binary output
  bool write_binary(const unsigned step, const double time);
  
  
  // @brief meta file output
  bool write_filelist(const unsigned step);
  
  
  // @brief ascii output
  bool write_ascii(const unsigned step, const double time);
  
  
  // @brief ログ出力
  void logging(const unsigned step);
  
  
  // @brief マイグレーション処理
  bool migration();
  
  
  // @brief PMlib ラベル
  void set_timing_label();
  
  // @brief パラメータ表示
  void displayParam(FILE* fp);
  
  
  // @brief 受信粒子のアンパック
  void unpackParticle();
  
  
  // @brief 開始点のユニークIDを割り振る
  bool determineUniqueID();
  
  
  /// @brief 粒子追跡情報を取得し，chunkに保持する
  bool setPTinfo();
  
  
  
  /// @brief 座標値を無次元化する
  /// @param [in,out] x  coordinate
  void normalizeCord(REAL_TYPE* x)
  {
    x[0] /= refLen;
    x[1] /= refLen;
    x[2] /= refLen;
  }
  
  
  /// @brief 単位ベクトル化
  /// @param [in,out] x  coordinate
  void getUnitVector(REAL_TYPE* x)
  {
    REAL_TYPE s = sqrt( x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
    x[0] /= s;
    x[1] /= s;
    x[2] /= s;
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
  bool samplingInCircle(const REAL_TYPE* cnt,
                        const REAL_TYPE* nv,
                        const REAL_TYPE radius,
                        const int nSample);
  
  
  /// @brief 指定法線nvがz軸の方向ベクトルに向かう回転角を計算する
  /// @param [in]  nv  指定法線
  /// @retval  回転角
  Vec3r getAngle(Vec3r nv);
  
  
  // @brief 回転ベクトルp(alpha, beta, gamma)でベクトルuを回転する
  Vec3r rotate(const Vec3r p, const Vec3r u);
  
  
  // @brief 回転ベクトルp(alpha, beta, gamma)に対して，-pでベクトルuを回転する
  Vec3r rotate_inv(const Vec3r p, const Vec3r u);
  
  
  /**
   * @brief タイミング測定開始
   * @param [in] key ラベル
   */
  inline void TIMING_start(const string key)
  {
    // PMlib Intrinsic profiler
    PT_TIMING__ PM->start(key);
    
    const char* s_label = key.c_str();
    
    // Venus FX profiler
#if defined __K_FPCOLL
    start_collection( s_label );
#elif defined __FX_FAPP
    fapp_start( s_label, 0, 0);
#endif
  }
  
  
  /**
   * @brief タイミング測定終了
   * @param [in] key             ラベル
   * @param [in] flopPerTask    「タスク」あたりの計算量/通信量(バイト) (ディフォルト0)
   * @param [in] iterationCount  実行「タスク」数 (ディフォルト1)
   */
  inline void TIMING_stop(const string key, double flopPerTask=0.0, int iterationCount=1)
  {
    // Venus FX profiler
    const char* s_label = key.c_str();
    
#if defined __K_FPCOLL
    stop_collection( s_label );
#elif defined __FX_FAPP
    fapp_stop( s_label, 0, 0);
#endif
    
    // PMlib Intrinsic profiler
    PT_TIMING__ PM->stop(key, flopPerTask, (unsigned)iterationCount);
  }
  
  // ディレクトリがなければ作成、既存なら何もしない（単一ディレクトリ）
  bool c_mkdir(const char* path);
  
};

#endif // _PT_CLOUD_H_
