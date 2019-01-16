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
#define MASK_25  0x1ffffff // 25 bit幅
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
  Vec3r vel;             ///< 粒子の並進速度成分
  int foo;               ///< padding とりあえずデバッグ用に使う
} particle;              //   他に必要な情報はここにいれる
                         //   ビットフィールド
                         //   31 - active(1) / inactive(0)
                         //   30 - migrate (1) / stay (0)
                         //   25~29 移動先ランク(周囲の26方向)
                         //   24~0  放出後のライフカウント 0 ~ 33,554,431
                         //         ライフカウントは放出回数、計算のステップ数ではない


/// 粒子塊　ストリークラインイメージ
class Chunk {
protected:
  int grp;                 ///< 管理用グループ番号Cloudから利用
  int uid;                 ///< 開始点の固有ID
  Vec3r EmitOrigin;        ///< Chunkに割り振られた開始点（固定）
  int startOrigin;         ///< 初期の指定開始点のチャンク(1) / マイグレート先(0)
  int EmitStep;            ///< 放出開始ステップ
  int myRank;              ///< ランク番号
  int interval;            ///< 放出間隔 step
  bool mig;                ///< マイグレーションフラグ true-あり, false-なし
  list<particle> pchunk;   ///< 粒子チャンク

public:
  Chunk() {
    uid = 0;
    grp = -1;
    startOrigin = 0;
    EmitStep = -1;
    interval = 0;
    mig = false;
  }

  /// コンストラクタ 作成時
  // @param [in] p       粒子座標
  // @param [in] gp      グループID
  // @param [in] st      オリジナル開始点のとき1
  // @param [in] stp     開始ステップ
  // @param [in] m_rank  ランク番号
  // @param [in] m_intvl 放出間隔
  Chunk(const Vec3r v,
        const int gp,
        const int st,
        const int stp,
        const int m_rank,
        const int m_intvl) {
    particle p;
    p.pos = v;
    p.bf  = 0;
    p.foo = 0;
    p.vel.assign(0.0, 0.0, 0.0);
    p.bf = Activate(p.bf);
    pchunk.push_back(p);
    grp = gp;
    EmitOrigin = v;
    startOrigin = st;
    EmitStep = stp;
    myRank = m_rank;
    interval = m_intvl;
    mig = false;
  }

  /// コンストラクタ マイグレーション追加時
  // @param [in] p       粒子座標
  // @param [in] gp      グループID
  // @param [in] pid     開始点ID
  // @param [in] stp     放出開始ステップ
  // @param [in] m_rank  ランク番号
  // @param [in] m_intvl 放出間隔
  Chunk(particle p,
        const int gp,
        int pid,
        const int stp,
        const int m_rank,
        const int m_intvl) {
    pchunk.push_back(p);
    grp = gp;
    uid = pid;
    EmitStep = stp;
    startOrigin = 0;
    myRank = m_rank;
    interval = m_intvl;
    mig = false;
  }

  ~Chunk() {}



  // @brief pchunkに保持している粒子を積分し、マイグレーションと寿命判定
  void updatePosition(Tracking* tr,
                      const int scheme,
                      const int life,
                      const REAL_TYPE dt,
                      int& max_particle,
                      const int* Rmap);


  // @brief マイグレーションフラグの立っている粒子をパックし、リストから削除
  void packParticle(REAL_TYPE* pbuf,
                    int* fbuf,
                    const int bsz,
                    int* pInfo);


  // @brief pchunkに開始点から粒子を追加
  void addParticleFromOrigin();


  // @brief pchunkの終端に粒子を追加
  // @note ライフカウントはparticleが持っている値を継承
  void addParticle(particle p)
  {
    pchunk.push_back(p);
  }


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
    return (startOrigin==1) ? true : false;
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
  unsigned getNpoints() const
  {
    return (unsigned)pchunk.size();
  }
  
  
  // @brief ascii output
  void write_ascii(FILE* fp);


  // @brief MIGRATE_BITと方向をOFF
  static inline int removeMigrate(const int b)
  {
    int c = b;
    c &= (~(0x2f << DST_BIT) ); // 対象6bitをゼロにする
    return c;
  }
  
  
  // @brief 25bit幅の正数をとりだす
  // @param [in] b  int variable
  static int getBit25 (const int b)
  {
    return ( b & MASK_25 );
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
  
  
  // @brief 方向 [0, 26] を返す
  // @note 自領域は13
  inline int getMigrateDir(const Vec3i d) const {
    int m = d.z*3*3 + d.y*3 + d.x + 13;
    return m;
  }


  // @brief 方向をデコード
  // @retval [0, 26]
  inline int decMigrateDir(const int b) const {
    return ( (b >> DST_BIT) & MASK_5 );
  }


  // @brief 25bit幅の値の設定
  // @param [in,out] b   int 変数
  // @param [in]     q   25-bit幅の数
  inline int setBit25 (int& b, const int q)
  {
    if ( q > MAX_LIFE ) exit(0);
    b &= ~MASK_25; // 対象25bitをゼロにする
    b |= q;        // 書き込む
  }


  // @brief ACTIVE_BITをON
  inline int Activate(int idx)
  {
    return ( idx | (0x1<<ACTIVE_BIT) );
  }


  // @brief ACTIVE_BITをOFF
  inline int Inactivate(int idx)
  {
    return ( idx & (~(0x1<<ACTIVE_BIT)) );
  }


  // @brief MIGRATE_BITをON
  inline int stampMigrate(int idx)
  {
    return ( idx | (0x1<<MIGRATE_BIT) );
  }

};

#endif // _PT_PARTICLE_H_
