/*
 * CPM - Cartesian Partition Manager
 *
 * Copyright (C) 2012 University of Tokyo.
 *
 */

/**
 * @file   cpm_ParaManager.cpp
 * パラレルマネージャクラスのソースファイル
 * @author University of Tokyo
 * @date   2012/05/31
 */
#include "cpm_ParaManager.h"

////////////////////////////////////////////////////////////////////////////////
// 唯一のインスタンスの取得
cpm_ParaManager*
cpm_ParaManager::get_instance()
{
  // 宣言
  static cpm_ParaManager instance;

  // ポインタ
  return &instance;
}

////////////////////////////////////////////////////////////////////////////////
// 唯一のインスタンスの取得(initialize処理も実行)
cpm_ParaManager*
cpm_ParaManager::get_instance(int &argc, char**& argv)
{
  // 宣言
  static cpm_ParaManager instance;

  // MPI_Init
  if( instance.Initialize(argc,argv) != CPM_SUCCESS )
  {
    return NULL;
  }

  // ポインタ
  return &instance;
}

////////////////////////////////////////////////////////////////////////////////
// コンストラクタ
cpm_ParaManager::cpm_ParaManager()
  : cpm_Base()
{
  // 並列数、ランク番号
  m_nRank  = 1;
  m_rankNo = 0;

  // プロセスグループリストをクリア
  m_procGrpList.clear();
  m_procGrpList.push_back(MPI_COMM_WORLD);

  // 領域管理情報マップのクリア
  m_voxelInfoMap.clear();
  m_rankNoMap.clear();

  // 袖通信バッファ情報のクリア
  m_bndCommInfoMap.clear();
}

////////////////////////////////////////////////////////////////////////////////
// デストラクタ
cpm_ParaManager::~cpm_ParaManager()
{
  // プロセスグループリストのクリア
  m_procGrpList.clear();

  // 領域管理情報マップの削除、クリア
  {
    VoxelInfoMap::iterator it  = m_voxelInfoMap.begin();
    VoxelInfoMap::iterator ite = m_voxelInfoMap.end();
    for( ; it!=ite; it++ )
    {
      if( it->second ) delete it->second;
    }
    m_voxelInfoMap.clear();
  }

  // ランク番号マップの削除、クリア
  {
    RankNoMap::iterator it  = m_rankNoMap.begin();
    RankNoMap::iterator ite = m_rankNoMap.end();
    for( ; it!=ite; it++ )
    {
      if( it->second ) delete [] it->second;
    }
    m_rankNoMap.clear();
  }

  // 袖通信バッファ情報の削除、クリア
  {
    BndCommInfoMap::iterator it  = m_bndCommInfoMap.begin();
    BndCommInfoMap::iterator ite = m_bndCommInfoMap.end();
    for( ; it!=ite; it++ )
    {
      if( it->second ) delete it->second;
    }
    m_bndCommInfoMap.clear();
  }

  // MPIのFinalize
  int flag1, flag2;
  MPI_Initialized(&flag1);
  MPI_Finalized(&flag2);
  if( flag1==true && flag2==false ) MPI_Finalize();
}

////////////////////////////////////////////////////////////////////////////////
// 初期化処理(MPI_Initは実行済みである必要がある)
cpm_ErrorCode
cpm_ParaManager::Initialize()
{
  m_nRank  = 1;
  m_rankNo = 0;

  // MPI_Init実行済みかチェック
  int flag1;
  MPI_Initialized(&flag1);
  if( flag1==false )
  {
    return CPM_ERROR_NO_MPI_INIT;
  }

  // プロセス並列数の取得
  if( MPI_Comm_size(MPI_COMM_WORLD, &m_nRank) != MPI_SUCCESS )
  {
    std::cerr << "MPI_Comm_size error." << std::endl;
    return CPM_ERROR_MPI;
  }
  if( m_nRank < 1 ) m_nRank = 1;

  // 自ランク番号の取得
  if( MPI_Comm_rank(MPI_COMM_WORLD, &m_rankNo) != MPI_SUCCESS )
  {
    std::cerr << "MPI_Comm_rank error." << std::endl;
    return CPM_ERROR_MPI;
  }
  if( m_rankNo < 0 ) m_rankNo = 0;

#ifdef _DEBUG
  flush(std::cout);
  for( int i=0;i<m_nRank;i++ )
  {
    if( i==m_rankNo )
      std::cout << "[" << m_rankNo << "] #rank=" << m_nRank << " parallel=" << IsParallel() << std::endl;
    flush(std::cout);
  }
  flush(std::cout);
#endif

  return CPM_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
// 初期化処理(MPI_Initも実行)
cpm_ErrorCode
cpm_ParaManager::Initialize( int &argc, char**& argv )
{
  m_nRank  = 1;
  m_rankNo = 0;

  // MPI_Init
  int flag1;
  MPI_Initialized(&flag1);
  if( flag1==false )
  {
    if( MPI_Init(&argc,&argv) != MPI_SUCCESS )
    {
      std::cerr << "MPI_Init error." << std::endl;
      return CPM_ERROR_MPI;
    }
  }

  // initialize
  return Initialize();
}

////////////////////////////////////////////////////////////////////////////////
// 並列実行であるかチェックする
bool
cpm_ParaManager::IsParallel()
{
  if( m_nRank <= 1 )
  {
    return false;
  }
  return true;
}

////////////////////////////////////////////////////////////////////////////////
// 並列実行であるかチェックする(const)
bool
cpm_ParaManager::IsParallel() const
{
  if( m_nRank <= 1 )
  {
    return false;
  }
  return true;
}

////////////////////////////////////////////////////////////////////////////////
// 領域分割
cpm_ErrorCode
cpm_ParaManager::VoxelInit( cpm_GlobalDomainInfo* domainInfo, size_t maxVC, size_t maxN, int procGrpNo )
{
  cpm_ErrorCode ret;
  cpm_VoxelInfo *voxelInfo = NULL;

  // 入力のチェック
  if( !domainInfo )
  {
    return CPM_ERROR_INVALID_PTR;
  }

  // コミュニケータを取得
  MPI_Comm comm = GetMPI_Comm( procGrpNo );
  if( IsCommNull( comm ) )
  {
    return CPM_ERROR_MPI_INVALID_COMM;
  }

  // 既に領域分割済みか
  if( m_voxelInfoMap.find(procGrpNo) != m_voxelInfoMap.end() )
  {
    return CPM_ERROR_ALREADY_VOXELINIIT;
  }

  // commの並列数を取得
  int nRank = 0;
  MPI_Comm_size( comm, &nRank );

  // プロセス数と活性サブドメイン数のチェック
  if( nRank != domainInfo->GetSubDomainNum() )
  {
    Abort(CPM_ERROR_MISMATCH_NP_SUBDOMAIN);
    return CPM_ERROR_MISMATCH_NP_SUBDOMAIN;
  }

  // インスタンス
  voxelInfo = new cpm_VoxelInfo();
  if( !voxelInfo )
  {
    Abort(CPM_ERROR_INVALID_PTR);
    return CPM_ERROR_INVALID_PTR;
  }

  // 領域分割情報の生成
  if( (ret = voxelInfo->Init( comm, domainInfo )) != CPM_SUCCESS )
  {
    delete voxelInfo;
    Abort(ret);
    return ret;
  }

  // VOXEL空間マップに登録
  if( !m_voxelInfoMap.insert(std::make_pair(procGrpNo, voxelInfo)).second )
  {
    m_procGrpList.pop_back();
    delete voxelInfo;
    return CPM_ERROR_INSERT_VOXELMAP;
  }

  // 袖通信バッファの設定
  SetBndCommBuffer( maxVC, maxN, procGrpNo );

  return CPM_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
// 領域分割(領域分割数を指定)
cpm_ErrorCode
cpm_ParaManager::VoxelInit( int div[3], int vox[3], REAL_TYPE origin[3], REAL_TYPE pitch[3]
                          , int obcid[6], size_t maxVC, size_t maxN, int procGrpNo )
{
  // DomainInfoを生成
  REAL_TYPE region[3];
  for( int i=0;i<3;i++ )
  {
    region[i] = pitch[i] * REAL_TYPE(vox[i]);
  }
  cpm_GlobalDomainInfo dInfo;
  dInfo.SetOrigin( origin );
  dInfo.SetPitch ( pitch );
  dInfo.SetRegion( region );
  dInfo.SetVoxNum( vox );
  dInfo.SetDivNum( div );

  // subdomain
  for( int k=0;k<div[2];k++ ){
  for( int j=0;j<div[1];j++ ){
  for( int i=0;i<div[0];i++ ){
    int pos[3] = {i, j, k};
    int bcid[6] = {0, 0, 0, 0, 0, 0};
    if( i==0 ) bcid[X_MINUS] = obcid[X_MINUS];
    if( j==0 ) bcid[Y_MINUS] = obcid[Y_MINUS];
    if( k==0 ) bcid[Z_MINUS] = obcid[Z_MINUS];
    if( i==div[0]-1 ) bcid[X_PLUS] = obcid[X_PLUS];
    if( j==div[1]-1 ) bcid[Y_PLUS] = obcid[Y_PLUS];
    if( k==div[2]-1 ) bcid[Z_PLUS] = obcid[Z_PLUS];

    cpm_ActiveSubDomainInfo dom( pos, bcid );
    dInfo.AddSubDomain( dom );
  }}}

  // 共通の処理
  return VoxelInit( &dInfo, maxVC, maxN, procGrpNo );
}

////////////////////////////////////////////////////////////////////////////////
// 領域分割(プロセスグループのランク数で自動領域分割)
cpm_ErrorCode
cpm_ParaManager::VoxelInit( int vox[3], REAL_TYPE origin[3], REAL_TYPE pitch[3]
                          , int obcid[6], size_t maxVC, size_t maxN, int procGrpNo )
{
  // ランク数
  int nrank = GetNumRank( procGrpNo );

  // 領域分割(とりあえずMPI_Dims_createで)
  int dims[3] = {0,0,0};
  if( MPI_Dims_create( nrank, 3, dims ) != MPI_SUCCESS )
  {
    return CPM_ERROR;
  }

  // 領域分割
  return VoxelInit( dims, vox, origin, pitch, obcid, maxVC, maxN, procGrpNo );
}

////////////////////////////////////////////////////////////////////////////////
// プロセスグループの作成
int
cpm_ParaManager::CreateProcessGroup( int nproc, int *proclist, int parentProcGrpNo )
{
  // 入力チェック
  if( parentProcGrpNo < 0 || nproc < 1 || !proclist )
  {
    return -1;
  }

  // 親のMPIコミュニケータを取得
  MPI_Comm parentComm = GetMPI_Comm( parentProcGrpNo );
  if( IsCommNull( parentComm ) )
  {
    return -1;
  }

  // 親のグループを取得
  MPI_Group parentGrp;
  MPI_Comm_group(parentComm, &parentGrp);

  // 新しいグループを生成
  MPI_Group newGrp;
  MPI_Group_incl(parentGrp, nproc, proclist, &newGrp);

  // 新しいコミュニケータを生成
  MPI_Comm  newComm;
  MPI_Comm_create(parentComm, newGrp, &newComm);
  if( IsCommNull(newComm) )
  {
    // 時ランクが含まれない
    return -1;
  }

  // コミュニケータを登録
  m_procGrpList.push_back(newComm);

  // プロセスグループの番号
  return m_procGrpList.size()-1;
}

////////////////////////////////////////////////////////////////////////////////
// VOXEL空間マップを検索
const cpm_VoxelInfo*
cpm_ParaManager::FindVoxelInfo( int procGrpNo )
{
  VoxelInfoMap::iterator it = m_voxelInfoMap.find(procGrpNo);
  if( it == m_voxelInfoMap.end() ) return NULL;
  return it->second;
}

////////////////////////////////////////////////////////////////////////////////
// 領域分割数を取得
const int*
cpm_ParaManager::GetDivNum( int procGrpNo )
{
  //VOXEL空間マップを検索
  const cpm_VoxelInfo *pVoxelInfo = FindVoxelInfo( procGrpNo );
  if( !pVoxelInfo ) return NULL;

  return pVoxelInfo->GetDivNum();
}

////////////////////////////////////////////////////////////////////////////////
// ピッチを取得
const REAL_TYPE*
cpm_ParaManager::GetPitch( int procGrpNo )
{
  //VOXEL空間マップを検索
  const cpm_VoxelInfo *pVoxelInfo = FindVoxelInfo( procGrpNo );
  if( !pVoxelInfo ) return NULL;

  return pVoxelInfo->GetPitch();
}

////////////////////////////////////////////////////////////////////////////////
// 全体ボクセル数を取得
const int*
cpm_ParaManager::GetGlobalVoxelSize( int procGrpNo )
{
  //VOXEL空間マップを検索
  const cpm_VoxelInfo *pVoxelInfo = FindVoxelInfo( procGrpNo );
  if( !pVoxelInfo ) return NULL;

  return pVoxelInfo->GetGlobalVoxelSize();
}

////////////////////////////////////////////////////////////////////////////////
// 全体空間の原点を取得
const REAL_TYPE*
cpm_ParaManager::GetGlobalOrigin( int procGrpNo )
{
  //VOXEL空間マップを検索
  const cpm_VoxelInfo *pVoxelInfo = FindVoxelInfo( procGrpNo );
  if( !pVoxelInfo ) return NULL;

  return pVoxelInfo->GetGlobalOrigin();
}

////////////////////////////////////////////////////////////////////////////////
// 全体空間サイズを取得
const REAL_TYPE*
cpm_ParaManager::GetGlobalRegion( int procGrpNo )
{
  //VOXEL空間マップを検索
  const cpm_VoxelInfo *pVoxelInfo = FindVoxelInfo( procGrpNo );
  if( !pVoxelInfo ) return NULL;

  return pVoxelInfo->GetGlobalRegion();
}

////////////////////////////////////////////////////////////////////////////////
// 自ランクのボクセル数を取得
const int*
cpm_ParaManager::GetLocalVoxelSize( int procGrpNo )
{
  //VOXEL空間マップを検索
  const cpm_VoxelInfo *pVoxelInfo = FindVoxelInfo( procGrpNo );
  if( !pVoxelInfo ) return NULL;

  return pVoxelInfo->GetLocalVoxelSize();
}

////////////////////////////////////////////////////////////////////////////////
// 自ランクの空間原点を取得
const REAL_TYPE*
cpm_ParaManager::GetLocalOrigin( int procGrpNo )
{
  //VOXEL空間マップを検索
  const cpm_VoxelInfo *pVoxelInfo = FindVoxelInfo( procGrpNo );
  if( !pVoxelInfo ) return NULL;

  return pVoxelInfo->GetLocalOrigin();
}

////////////////////////////////////////////////////////////////////////////////
// 自ランクの空間サイズを取得
const REAL_TYPE*
cpm_ParaManager::GetLocalRegion( int procGrpNo )
{
  //VOXEL空間マップを検索
  const cpm_VoxelInfo *pVoxelInfo = FindVoxelInfo( procGrpNo );
  if( !pVoxelInfo ) return NULL;

  return pVoxelInfo->GetLocalRegion();
}

////////////////////////////////////////////////////////////////////////////////
// 自ランクの領域分割位置を取得
const int*
cpm_ParaManager::GetDivPos( int procGrpNo )
{
  //VOXEL空間マップを検索
  const cpm_VoxelInfo *pVoxelInfo = FindVoxelInfo( procGrpNo );
  if( !pVoxelInfo ) return NULL;

  return pVoxelInfo->GetDivPos();
}

////////////////////////////////////////////////////////////////////////////////
// 自ランクの領域分割位置を取得
const int*
cpm_ParaManager::GetBCID( int procGrpNo )
{
  //VOXEL空間マップを検索
  const cpm_VoxelInfo *pVoxelInfo = FindVoxelInfo( procGrpNo );
  if( !pVoxelInfo ) return NULL;

  return pVoxelInfo->GetBCID();
}

////////////////////////////////////////////////////////////////////////////////
// 自ランクの始点VOXELの全体空間でのインデクスを取得
const int*
cpm_ParaManager::GetVoxelHeadIndex( int procGrpNo )
{
  //VOXEL空間マップを検索
  const cpm_VoxelInfo *pVoxelInfo = FindVoxelInfo( procGrpNo );
  if( !pVoxelInfo ) return NULL;

  return pVoxelInfo->GetVoxelHeadIndex();
}

////////////////////////////////////////////////////////////////////////////////
// 自ランクの終点VOXELの全体空間でのインデクスを取得
const int*
cpm_ParaManager::GetVoxelTailIndex( int procGrpNo )
{
  //VOXEL空間マップを検索
  const cpm_VoxelInfo *pVoxelInfo = FindVoxelInfo( procGrpNo );
  if( !pVoxelInfo ) return NULL;

  return pVoxelInfo->GetVoxelTailIndex();
}

////////////////////////////////////////////////////////////////////////////////
// 自ランクの隣接ランク番号を取得
const int*
cpm_ParaManager::GetNeighborRankID( int procGrpNo )
{
  //VOXEL空間マップを検索
  const cpm_VoxelInfo *pVoxelInfo = FindVoxelInfo( procGrpNo );
  if( !pVoxelInfo ) return NULL;

  return pVoxelInfo->GetNeighborRankID();
}

////////////////////////////////////////////////////////////////////////////////
// 自ランクの周期境界の隣接ランク番号を取得
const int*
cpm_ParaManager::GetPeriodicRankID( int procGrpNo )
{
  //VOXEL空間マップを検索
  const cpm_VoxelInfo *pVoxelInfo = FindVoxelInfo( procGrpNo );
  if( !pVoxelInfo ) return NULL;

  return pVoxelInfo->GetPeriodicRankID();
}

////////////////////////////////////////////////////////////////////////////////
// 袖通信バッファのセット
cpm_ErrorCode
cpm_ParaManager::SetBndCommBuffer( size_t maxVC, size_t maxN, int procGrpNo )
{
  if( maxVC==0 || maxN==0 )
  {
    return CPM_ERROR_BNDCOMM;
  }

  // local voxel size
  const int *sz = GetLocalVoxelSize(procGrpNo);
  if( !sz ) return CPM_ERROR_BNDCOMM_VOXELSIZE;

  // buffer size
  size_t nwX = size_t(sz[1]+2*maxVC) * size_t(sz[2]+2*maxVC) * maxVC * maxN;
  size_t nwY = size_t(sz[2]+2*maxVC) * size_t(sz[0]+2*maxVC) * maxVC * maxN;
  size_t nwZ = size_t(sz[0]+2*maxVC) * size_t(sz[1]+2*maxVC) * maxVC * maxN;

  // 送受信バッファ情報のインスタンス
  S_BNDCOMM_BUFFER *bufInfo = new S_BNDCOMM_BUFFER();
  if( !bufInfo )
  {
    return CPM_ERROR_INVALID_PTR;
  }

  bufInfo->m_maxVC = maxVC;
  bufInfo->m_maxN  = maxN;
  bufInfo->m_nwX   = nwX;
  bufInfo->m_nwY   = nwY;
  bufInfo->m_nwZ   = nwZ;

  // buffer
  for( int i=0;i<4;i++ )
  {
    bufInfo->m_bufX[i] = new REAL_BUF_TYPE[nwX];
    if( !bufInfo->m_bufX[i] )
    {
      delete bufInfo;
      return CPM_ERROR_INVALID_PTR;
    }
    bufInfo->m_bufY[i] = new REAL_BUF_TYPE[nwY];
    if( !bufInfo->m_bufY[i] )
    {
      delete bufInfo;
      return CPM_ERROR_INVALID_PTR;
    }
    bufInfo->m_bufZ[i] = new REAL_BUF_TYPE[nwZ];
    if( !bufInfo->m_bufZ[i] )
    {
      delete bufInfo;
      return CPM_ERROR_INVALID_PTR;
    }
  }

  //マップにセット
  BndCommInfoMap::iterator it = m_bndCommInfoMap.find(procGrpNo);
  if( it != m_bndCommInfoMap.end() )
  {
    delete it->second;
    m_bndCommInfoMap.erase( it );
  }
  m_bndCommInfoMap.insert( std::make_pair(procGrpNo, bufInfo) );

  return CPM_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
// 袖通信バッファサイズの取得
size_t
cpm_ParaManager::GetBndCommBufferSize( int procGrpNo )
{
  size_t mem = 0;

  // マップを検索
  if( procGrpNo < 0 )
  {
    //全体
    BndCommInfoMap::iterator it = m_bndCommInfoMap.begin();
    for( ; it!=m_bndCommInfoMap.end(); it++ )
    {
      if( !it->second ) continue;
      size_t nwX   = it->second->m_nwX;
      size_t nwY   = it->second->m_nwY;
      size_t nwZ   = it->second->m_nwZ;
      mem += ((nwX*4 + nwY*4 + nwZ*4) * sizeof(REAL_BUF_TYPE));
    }
  }
  else
  {
    // 指定プロセスグループを検索
    S_BNDCOMM_BUFFER *bInfo = GetBndCommBuffer( procGrpNo );
    if( !bInfo ) return 0;
    size_t nwX   = bInfo->m_nwX;
    size_t nwY   = bInfo->m_nwY;
    size_t nwZ   = bInfo->m_nwZ;
    mem = (nwX*4 + nwY*4 + nwZ*4) * sizeof(REAL_BUF_TYPE);
  }

  return mem;
}



////////////////////////////////////////////////////////////////////////////////
// 




////////////////////////////////////////////////////////////////////////////////
// flush
void
cpm_ParaManager::flush(std::ostream &out, int procGrpNo)
{
  for( int i=0;i<10;i++ )
  {
    Barrier(procGrpNo);
    out << std::flush;
    Barrier(procGrpNo);
  }
}

////////////////////////////////////////////////////////////////////////////////
// flush
void
cpm_ParaManager::flush(FILE *fp, int procGrpNo)
{
  if( !fp ) return;
  for( int i=0;i<10;i++ )
  {
    Barrier(procGrpNo);
    fflush(fp);
    Barrier(procGrpNo);
  }
}

#ifdef _DEBUG
#include <fstream>
// debug write
void
cpm_ParaManager::printVoxelInfo(int myrank)
{
  if( myrank < 0 )
  {
    flush(std::cout);
    for( int i=0;i<m_nRank;i++ )
    {
      if( i==m_rankNo )
      {
        std::cout << "####### VoxelInfo [" << i << " @ MPI_COMM_WORLD] #######" << std::endl;
        VoxelInfoMap::iterator itv = m_voxelInfoMap.begin();
        for( ;itv!=m_voxelInfoMap.end();itv++ )
        {
          // プロセスグループ番号、ボクセル情報、ランク番号
          int procGrpNo = itv->first;
          cpm_VoxelInfo *pV = itv->second;
          if( !pV ) continue;
          int rankNo = GetMyRankID(procGrpNo);
          std::cout << " *process group No [" << procGrpNo << "], rankNo[" << rankNo << "]" << std::endl;

          // 全体空間
          const int       *gdiv = pV->GetDivNum();
          const REAL_TYPE *gorg = pV->GetGlobalOrigin();
          const REAL_TYPE *gpch = pV->GetPitch();
          const REAL_TYPE *grgn = pV->GetGlobalRegion();
          const int       *gvox = pV->GetGlobalVoxelSize();
          std::cout << "  +----------------------------------" << std::endl;
          std::cout << "  global div = " << gdiv[0] << "," << gdiv[1] << "," << gdiv[2] << std::endl;
          std::cout << "  global org = " << gorg[0] << "," << gorg[1] << "," << gorg[2] << std::endl;
          std::cout << "  global pch = " << gpch[0] << "," << gpch[1] << "," << gpch[2] << std::endl;
          std::cout << "  global rgn = " << grgn[0] << "," << grgn[1] << "," << grgn[2] << std::endl;
          std::cout << "  global vox = " << gvox[0] << "," << gvox[1] << "," << gvox[2] << std::endl;

          // ローカル空間
          const REAL_TYPE *lorg = pV->GetLocalOrigin();
          const REAL_TYPE *lpch = pV->GetPitch();
          const REAL_TYPE *lrgn = pV->GetLocalRegion();
          const int       *lvox = pV->GetLocalVoxelSize();
          const int       *lpos = pV->GetDivPos();
          const int       *bcid = pV->GetBCID();
          const int       *head = pV->GetVoxelHeadIndex();
          const int       *neig = pV->GetNeighborRankID();
          const int       *peri = pV->GetPeriodicRankID();
          std::cout << "  +----------------------------------" << std::endl;
          std::cout << "  local  org = " << lorg[0] << "," << lorg[1] << "," << lorg[2] << std::endl;
          std::cout << "  local  pch = " << lpch[0] << "," << lpch[1] << "," << lpch[2] << std::endl;
          std::cout << "  local  rgn = " << lrgn[0] << "," << lrgn[1] << "," << lrgn[2] << std::endl;
          std::cout << "  local  vox = " << lvox[0] << "," << lvox[1] << "," << lvox[2] << std::endl;
          std::cout << "  local  pos = " << lpos[0] << "," << lpos[1] << "," << lpos[2] << std::endl;
          std::cout << "  local  head= " << head[0] << "," << head[1] << "," << head[2] << std::endl;
          std::cout << "  local  bcid= " << bcid[0] << "," << bcid[1] << "," << bcid[2] << ","
                                         << bcid[3] << "," << bcid[4] << "," << bcid[5] << std::endl;
          std::cout << "  local  neig= " << neig[0] << "," << neig[1] << "," << neig[2] << ","
                                         << neig[3] << "," << neig[4] << "," << neig[5] << std::endl;
          std::cout << "  local  peri= " << peri[0] << "," << peri[1] << "," << peri[2] << ","
                                         << peri[3] << "," << peri[4] << "," << peri[5] << std::endl;
        }
      }
      flush(std::cout);
    }
    flush(std::cout);
  }
  else
  {
    char fname[512];
    sprintf ( fname, "vinfo_%04d.log", myrank );
    std::ofstream ofs( fname );
    ofs << "####### VoxelInfo #######" << std::endl;
    VoxelInfoMap::iterator itv = m_voxelInfoMap.begin();
    for( ;itv!=m_voxelInfoMap.end();itv++ )
    {
      // プロセスグループ番号、ボクセル情報、ランク番号
      int procGrpNo = itv->first;
      cpm_VoxelInfo *pV = itv->second;
      if( !pV ) continue;
      int rankNo = GetMyRankID(procGrpNo);
      ofs << " *process group No [" << procGrpNo << "], rankNo[" << rankNo << "]" << std::endl;

      // 全体空間
      const int       *gdiv = pV->GetDivNum();
      const REAL_TYPE *gorg = pV->GetGlobalOrigin();
      const REAL_TYPE *gpch = pV->GetPitch();
      const REAL_TYPE *grgn = pV->GetGlobalRegion();
      const int       *gvox = pV->GetGlobalVoxelSize();
      ofs << "  +----------------------------------" << std::endl;
      ofs << "  global div = " << gdiv[0] << "," << gdiv[1] << "," << gdiv[2] << std::endl;
      ofs << "  global org = " << gorg[0] << "," << gorg[1] << "," << gorg[2] << std::endl;
      ofs << "  global pch = " << gpch[0] << "," << gpch[1] << "," << gpch[2] << std::endl;
      ofs << "  global rgn = " << grgn[0] << "," << grgn[1] << "," << grgn[2] << std::endl;
      ofs << "  global vox = " << gvox[0] << "," << gvox[1] << "," << gvox[2] << std::endl;

      // ローカル空間
      const REAL_TYPE *lorg = pV->GetLocalOrigin();
      const REAL_TYPE *lpch = pV->GetPitch();
      const REAL_TYPE *lrgn = pV->GetLocalRegion();
      const int       *lvox = pV->GetLocalVoxelSize();
      const int       *lpos = pV->GetDivPos();
      const int       *bcid = pV->GetBCID();
      const int       *head = pV->GetVoxelHeadIndex();
      const int       *neig = pV->GetNeighborRankID();
      const int       *peri = pV->GetPeriodicRankID();
      ofs << "  +----------------------------------" << std::endl;
      ofs << "  local  org = " << lorg[0] << "," << lorg[1] << "," << lorg[2] << std::endl;
      ofs << "  local  pch = " << lpch[0] << "," << lpch[1] << "," << lpch[2] << std::endl;
      ofs << "  local  rgn = " << lrgn[0] << "," << lrgn[1] << "," << lrgn[2] << std::endl;
      ofs << "  local  vox = " << lvox[0] << "," << lvox[1] << "," << lvox[2] << std::endl;
      ofs << "  local  pos = " << lpos[0] << "," << lpos[1] << "," << lpos[2] << std::endl;
      ofs << "  local  head= " << head[0] << "," << head[1] << "," << head[2] << std::endl;
      ofs << "  local  bcid= " << bcid[0] << "," << bcid[1] << "," << bcid[2] << ","
                               << bcid[3] << "," << bcid[4] << "," << bcid[5] << std::endl;
      ofs << "  local  neig= " << neig[0] << "," << neig[1] << "," << neig[2] << ","
                               << neig[3] << "," << neig[4] << "," << neig[5] << std::endl;
      ofs << "  local  peri= " << peri[0] << "," << peri[1] << "," << peri[2] << ","
                               << peri[3] << "," << peri[4] << "," << peri[5] << std::endl;
    }

    ofs << "####### rank map #######" << std::endl;
    itv = m_voxelInfoMap.begin();
    for( ;itv!=m_voxelInfoMap.end();itv++ )
    {
      // プロセスグループ番号、ボクセル情報、ランク番号
      int procGrpNo = itv->first;
      cpm_VoxelInfo *pV = itv->second;
      if( !pV ) continue;
      int rankNo = GetMyRankID(procGrpNo);
      ofs << " *process group No [" << procGrpNo << "], rankNo[" << rankNo << "]" << std::endl;
      // マップ、領域分割数
      int *map = pV->m_rankMap;
      const int *div = pV->GetDivNum();
      ofs << " *div = " << div[0] << " x " << div[1] << " x " << div[2] << std::endl;
      ofs << " *map-------------------------------" << std::endl;
      for( int k=div[2]-1;k>=0;k-- )
      {
        ofs << "  k=" << k << std::endl;
        for( int j=div[1]-1;j>=0;j-- )
        {
           ofs << "    ";
           for( int i=0;i<div[0];i++ )
           {
             ofs.width(3);
             ofs.setf(std::ios::left, std::ios::adjustfield);
             ofs << map[_IDX_S3D(i,j,k,div[0],div[1],div[2],0)];
             if( i < div[0]-1 ) ofs << ",";
           }
           ofs << std::endl;
        }
      }
    }
    ofs.close();
  }
}
#endif
