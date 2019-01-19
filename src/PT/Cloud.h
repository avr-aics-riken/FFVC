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
protected:
  EmitGroup* Egrp;           ///< 粒子の開始点グループ
  vector<Chunk*> chunkList;  ///< チャンクリスト

  int nGrpEmit;              ///< 開始点のグループ数 set[@]の数
	unsigned nEmitParticle;    ///< 粒子放出点の数
	unsigned nEmission;        ///< 粒子放出の回数
  int scheme;                ///< 積分方法
  bool flag_migration;       ///< マイグレーション発生 true
  REAL_TYPE dt;              ///< 時間積分幅
  unsigned buf_max_particle; ///< 送信用のバッファ長さ計算に使う粒子の最大数 毎回異なる
  bool buf_updated;          ///< バッファ長さが更新されたときtrue

  int unit;                  ///< 指定座標の記述単位 {DIMENSIONAL | NONDIMENSIONAL}
  REAL_TYPE refLen;          ///< 代表長さ

  int* bcd;                  ///< BCindex B
  REAL_TYPE* vSource;        ///< 速度サンプリング元データ

  IntervalManager* Interval; ///< タイミング制御
  TextParser* tpCntl;        ///< TextParser
  Tracking* tr;              ///< Tracking
  PtComm PC;                 ///< 粒子通信クラス
  PerfMonitor* PM;           ///< PerfMonitor class
  
  int log_interval;          ///< log file out interval
  int file_interval;         ///< log file out interval
  int file_format;           ///< 0 - ascii, 1 - binary
  FILE* fpl;                 ///< ログ出力用のファイルポインタ
  int nCommParticle;         ///< マイグレーション時の送受信粒子数
  unsigned nParticle;        ///< 全粒子の数（ローカル）
  unsigned gParticle;        ///< 全粒子の数（グローバル）
  
  int* Rmap;                 ///< PtCommクラスで作成した3x3x3のランクマップへのポインタ



public:
  /// デフォルトコンストラクタ
  Cloud() {
    Rmap = NULL;
  }

  Cloud(int* m_bcd,
        REAL_TYPE* m_Vsrc,
        const REAL_TYPE dt,
        TextParser* m_tp,
        PerfMonitor* m_PM)
  {
    nParticle = 0;
    gParticle = 0;
    nGrpEmit = 0;
    scheme = -1;
    buf_updated = false;
    flag_migration = false;
    nCommParticle = 0;
    log_interval = 0;
    file_interval = 0;
    file_format = -1;
    unit = NONDIMENSIONAL;
    refLen = 0.0;
		nEmitParticle = 0;
		nEmission = 0;

    this->dt         = dt;
    this->bcd        = m_bcd;
    this->vSource    = m_Vsrc;
    this->tpCntl     = m_tp;
    this->PM         = m_PM;

    // 初期値として、BUF_UNIT*粒子分を確保 > buf_max_particleで100単位で更新
    buf_max_particle = BUF_UNIT;
  }



  /// デストラクタ
  ~Cloud() {
    delete [] Interval;
    if (fpl) fclose(fpl);
  }


  // @brief 初期設定
  bool initCloud(FILE* fp);


  // @brief ランタイム
  bool tracking(const unsigned step, const double time);

  
  
  

    
protected:
  
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


  // @brief 既に存在するグループIDの検索
  bool searchGrp(const int c);


  // @brief 開始点のユニークIDを割り振る
  bool determineUniqueID();


  /// @brief 粒子追跡情報を取得し，chunkに保持する
  bool setPTinfo();


  // @brief 粒子追跡パラメータ取得
  bool getTPparam(const string label_leaf, int odr);


  /// @brief 座標値を無次元化する
  /// @param [in,out] x  coordinate
  void normalizeCord(REAL_TYPE x[3])
  {
    x[0] /= refLen;
    x[1] /= refLen;
    x[2] /= refLen;
  }


  /* // @brief 自領域内に存在するかどうかを判断
  /// @param [in] x  coordinate
  bool inOwnRegion(const Vec3r p)
  {
    if ( p.x < origin[0] || p.x >= origin[0]+region[0]
      || p.y < origin[1] || p.y >= origin[1]+region[1]
      || p.z < origin[2] || p.z >= origin[2]+region[2] ) return false;
    return true;
  }
	*/


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



};

#endif // _PT_CLOUD_H_
