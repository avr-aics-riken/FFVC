#ifndef _PT_PARTICLE_H_
#define _PT_PARTICLE_H_

//##################################################################################
//
// Copyright (c) 2019 Research Institute for Information Technology, Kyushu university
// All rights researved.
//
//##################################################################################

/**
 * @file   Chunk.h
 * @brief  Chunk class Header
 * @author riit
 * @note 粒子単位の軌跡を管理
 */

#include <list>
#include <string.h>
#include "common/Vec3.h" // defined in Polylib
#include <stdio.h>
#include "FB_Define.h"
#include "PtDefine.h"
#include "EmitGroup.h"
#include "Tracking.h"

#define MAX_LIFE 33554431 // 2^{25}-1
#define MASK_26  0x7ffffff // 26 bit幅
#define MIGRATE_BIT 30
#define DST_BIT 25

#define IS_ACTIVE(a)  ( BIT_SHIFT(a, ACTIVE_BIT) )
#define IS_MIGRATE(a) ( BIT_SHIFT(a, MIGRATE_BIT) )


using namespace Vec3class;
using std::list;

enum IntegrationScheme
{
   pt_euler=0
  ,pt_rk2
  ,pt_rk4
};


// 粒子
typedef struct {
  Vec3r pos;             ///< 粒子座標
  int bf;                ///< ビットフィールド
} particle;              //   31 - active(1) / inactive(0)
                         //   30 - migrate (1) / stay (0)
                         //   25~29 移動先ランク(周囲の26方向)
                         //   24~0  放出後のライフカウント 0 ~ 33,554,431

/// 粒子群
class Chunk {
protected:
  int grp;                 ///< 管理用グループ番号Cloudから利用
  int uid;                 ///< 粒子の固有ID
  Vec3r origin;            ///< Chunkに割り振られた開始点（固定）
  int start_point;         ///< 初期の指定開始点(1) / マイグレート先(0)
  list<particle> pchunk;   ///< 粒子チャンク

public:
  Chunk() {
    uid = 0;
    grp = -1;
    start_point = 0;
  }

  /// コンストラクタ 作成時
  // @param [in] pv     粒子座標
  // @param [in] gp     グループID
  // @param [in] st     開始点のとき1
  Chunk(const Vec3r v, const int gp, const int st) {
    particle p;
    p.pos = v;
    p.bf = 0;
    Activate(p.bf);
    pchunk.push_back(p);
    grp = gp;
    origin = v;
    start_point = st;
  }

  /// コンストラクタ マイグレーション追加時
  // @param [in] pv     粒子座標
  // @param [in] gp     グループID
  // @param [in] pid    粒子ID
  Chunk(particle p, const int gp, const int pid) {
    Activate(p.bf);
    pchunk.push_back(p);
    grp = gp;
    uid = pid;
  }

  ~Chunk() {}



  // @brief pchunkに保持している粒子を積分し、マイグレーションと寿命判定
  // @retval マイグレーションが発生したら、trueを返す
  bool updatePosition(Tracking* tr,
                      const int scheme,
                      const int life,
                      const REAL_TYPE dt,
                      int& buf_len);


  // @brief マイグレーションフラグの立っている粒子をパックし、リストから削除
  void packParticle(REAL_TYPE* pbuf,
                    int* fbuf,
                    const int bsz,
                    int* pInfo);


  // @brief pchunkに開始点から粒子を追加
  void addParticle();


  // @brief pchunkに粒子を追加
  void addParticle(particle p);


  // @brief 粒子IDをセットする
  void setUid(const int m_id) {
    uid = m_id;
  }

  // @brief 粒子IDを返す
  int getUid() const {
    return uid;
  }


  // @brief 開始点かどうか
  bool isOrigin() const {
    return (start_point==1) ? true : false;
  }


  // @brief グループ番号をセットする
  void setGrp(const int m_grp) {
    grp = m_grp;
  }


  // @brief グループ番号を返す
  int getGrp() const {
    return grp;
  }


  // @brief pchunkのサイズを返す
  int getNpoints() const
  {
    return (int)pchunk.size();
  }


  // @brief ACTIVE_BITをOFF
  static inline void removeMigrate(int& idx)
  {
    idx & (~(0x1<<MIGRATE_BIT));
  }


protected:

  // @brief 送り先ランク方向をエンコードする
  // @param [in] 方向
  // @retval [0, 26]
  //
  //     LUT
  //          (k+1)                (k)               (k-1)
  //     +----+----+----+    +----+----+----+   +----+----+----+
  //     | 24 | 25 | 26 |    | 15 | 16 | 17 |   |  6 |  7 |  8 |
  //     +----+----+----+    +----+----+----+   +----+----+----+
  //     | 21 | 22 | 23 |    | 12 |    | 14 |   |  3 |  4 |  5 |
  //  j  +----+----+----+    +----+----+----+   +----+----+----+
  //  ^  | 18 | 19 | 20 |    |  9 | 10 | 11 |   |  0 |  1 |  2 |
  //  |  +----+----+----+    +----+----+----+   +----+----+----+
  //   --> i
  //
  inline int encMigrateDir(int& b, const Vec3i d) const {
     int m = d.z*3*3 + d.y*3 + d.x + 13;
     b &= (~(MASK_5 << DST_BIT) ); // 対象5bitをゼロにする
     b |= (m << DST_BIT);          // 書き込む
     return m;
  }


  // @brief 方向をデコード
  // @retval [0, 26]
  inline int decMigrateDir(const int b) const {
    return ( (b >> DST_BIT) & MASK_5 );
  }


  /*
   * @brief 26bitの正数をとりだす
   * @param [in] b  int variable
   */
  int getBit26 (const int b) const
  {
    return ( b & MASK_26 );
  }


  /*
   * @brief 26bit幅の値の設定
   * @param [in,out] b   int 変数
   * @param [in]     q   26-bit幅の数
   */
  void setBit26 (int& b, const int q)
  {
    if ( q > MAX_LIFE ) Exit(0);
    b &= ~MASK_26; // 対象26bitをゼロにする
    b |= q;        // 書き込む
  }


  // @brief ACTIVE_BITをON
  inline void Activate(int& idx)
  {
    idx | (0x1<<ACTIVE_BIT);
  }


  // @brief ACTIVE_BITをOFF
  inline void Inactivate(int& idx)
  {
    idx & (~(0x1<<ACTIVE_BIT));
  }


  // @brief MIGRATE_BITをON
  inline void stampMigrate(int& idx)
  {
    idx | (0x1<<MIGRATE_BIT);
  }

};

#endif // _PT_PARTICLE_H_
