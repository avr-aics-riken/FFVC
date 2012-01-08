/**@file
 * @brief 時間測定用クラス 宣言
 */

#ifndef CUTTIMING_H
#define CUTTIMING_H

#include <iostream>
#include <string>
#include <map>
#include <ctime>
#include <sys/time.h>
using namespace std;

namespace cutlib {

/// 時間測定ストップウオッチクラス
class Timing {
  string name_;     ///< ストップウオッチ名
  clock_t cStart_;  ///< スタートCPUクロック
  double wStart_;   ///< スタート経過時刻
  double cTime_;    ///< 累積CPU時間
  double wTime_;    ///< 累積経過時間
  int count_;       ///< 測定回数

  /// 時間測定ストップウオッチリスト
  static map<string, Timing> Timers_;

  /// コンストラクタ(private)
  /**
   * @param[in] name ストップウオッチ名
   */
  Timing(const string& name) : name_(name), cTime_(0.0), wTime_(0.0), count_(0) {}

  /// 経過時間取得(gettimeofdayシステムコール)
  double getWTime()
  {
    struct timeval tv;
    gettimeofday(&tv, 0);
    return (double)tv.tv_sec + (double)tv.tv_usec * 1.e-6;
  }

  /// 時間測定スタート
  void start() { cStart_ = clock(); wStart_ = getWTime(); }

  /// 時間測定ストップ
  void stop()
  {
    cTime_ += (double)(clock() - cStart_) / CLOCKS_PER_SEC;
    wTime_ += getWTime() - wStart_;
    ++count_;
  }

public:
  friend ostream& operator<<(ostream& os, const Timing& t);

  /// ストップウオッチnameをスタート
  /**
   * nameが未登録の場合は、新に作成
   * @param[in] name ストップウオッチ名
   */
  static void Start(const string& name);
 
  /// ストップウオッチnameをストップ
  /**
   * @param[in] name ストップウオッチ名
   */
  static void Stop(const string& name);

  /// ストップウオッチnameの内容を表示
  /**
   * @param[in] name ストップウオッチ名
   */
  static void Print(const string& name);

};

} /* namespace cutlib */

#endif /* CUTTIMING_H */
