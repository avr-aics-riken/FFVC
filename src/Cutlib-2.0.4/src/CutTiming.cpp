/**@file
 * @brief 時間測定用クラス 実装
 */

#include "CutTiming.h"

namespace cutlib {

ostream& operator<<(ostream& os, const Timing& t)
{
  os << t.name_ << ": count=" << t.count_ 
                << " cpu_time=" << t.cTime_ 
                << " elapsed_time="<< t.wTime_;
  return os;
}


/// 時間測定ストップウオッチリスト
map<string, Timing> Timing::Timers_;


/// ストップウオッチnameをスタート
/**
 * nameが未登録の場合は、新に作成
 * @param[in] name ストップウオッチ名
 */
void Timing::Start(const string& name)
{
  if (Timers_.find(name) == Timers_.end()) {
    Timers_.insert(pair<string, Timing>(name, Timing(name)));
  }
  map<string, Timing>::iterator itr = Timers_.find(name);
  itr->second.start();
}


/// ストップウオッチnameをストップ
/**
 * @param[in] name ストップウオッチ名
 */
void Timing::Stop(const string& name)
{
    map<string, Timing>::iterator itr = Timers_.find(name);
    if (itr != Timers_.end()) itr->second.stop();
}


/// ストップウオッチnameの内容を表示
/**
 * @param[in] name ストップウオッチ名
 */
void Timing::Print(const string& name)
{
    map<string, Timing>::iterator itr = Timers_.find(name);
    if (itr != Timers_.end()) {
      cout << itr->second << endl;
    }
}






} /* namespace cutlib */

