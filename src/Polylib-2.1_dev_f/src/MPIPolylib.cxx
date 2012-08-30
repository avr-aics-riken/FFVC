/*
 * Polylib - Polygon Management Library.
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2010-
 *
 */
#include <iostream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <map>
//#include "groups/PolygonGroup.h"
//#include "c_lang/CMPIPolylib.h"
#include "mpi.h"
#include "MPIPolylib.h"

// MPI通信用メッセージタグ
#define	MPITAG_NUM_CONFIG			1
#define	MPITAG_CONFIG				2
#define	MPITAG_NUM_TRIAS			3
#define MPITAG_TRIA_IDS				4
#define MPITAG_TRIAS				5


using namespace std;
using namespace PolylibNS;

////////////////////////////////////////////////////////////////////////////
/// 
/// クラス:MPIPolylib
/// Polylib継承したポリゴンを管理する為の並列版クラスライブラリです。
/// 
////////////////////////////////////////////////////////////////////////////


// public /////////////////////////////////////////////////////////////////////
MPIPolylib*
MPIPolylib::get_instance() {
	static MPIPolylib instance;
	return &instance;
}


// public /////////////////////////////////////////////////////////////////////
POLYLIB_STAT
MPIPolylib::init_parallel_info(
	MPI_Comm comm,
	float bpos[3], 
	unsigned int bbsize[3], 
	unsigned int gcsize[3], 
	float dx[3]
)
{
#ifdef DEBUG
	PL_DBGOSH << "MPIPolylib::init_parallel_info() in. " << endl;
#endif
	int i;

	// MPI情報の設定
	m_mycomm = comm;
	MPI_Comm_rank(comm, &m_myrank);
	MPI_Comm_size(comm, &m_numproc);

	// デバッグ出力用ランク番号文字列を設定
	std::ostringstream ostr;
	ostr << m_myrank;
	gs_rankno = "(rk:";
	gs_rankno += ostr.str();
	gs_rankno += ")";

#ifdef DEBUG
	PL_DBGOSH << "m_myrank: " << m_myrank << " m_numproc: " << m_numproc << endl;
#endif

	float bbsize_f[3], gcsize_f[3];
	for (i = 0; i < 3; i++) {
		bbsize_f[i] = (float)bbsize[i];
		gcsize_f[i] = (float)gcsize[i];
	}

	Vec3f v_bbsize(bbsize_f[0],bbsize_f[1],bbsize_f[2]);
	Vec3f v_gcsize(gcsize_f[0],gcsize_f[1],gcsize_f[2]);
	Vec3f v_bpos(bpos[0],bpos[1],bpos[2]);
	Vec3f v_dx(dx[0],dx[1],dx[2]);

	// 自PE領域情報を設定
	m_myproc.m_comm = comm;
	m_myproc.m_rank = m_myrank;
	m_myproc.m_area.m_bpos = v_bpos;
	m_myproc.m_area.m_bbsize = v_bbsize;
	m_myproc.m_area.m_gcsize = v_gcsize;
	m_myproc.m_area.m_dx = dx;
	m_myproc.m_area.m_gcell_min = v_bpos-( v_gcsize )*v_dx;
	m_myproc.m_area.m_gcell_max = v_bpos+( v_bbsize+v_gcsize )*v_dx;
	m_myproc.m_area.m_gcell_bbox.add(m_myproc.m_area.m_gcell_min);
	m_myproc.m_area.m_gcell_bbox.add(m_myproc.m_area.m_gcell_max);

#ifdef DEBUG
	PL_DBGOSH << "(my_rank:" << m_myrank << "):" <<"bpos      :" << v_bpos  << endl;
	PL_DBGOSH << "(my_rank:" << m_myrank << "):" <<"bbsize    :" << v_bbsize << endl;
	PL_DBGOSH << "(my_rank:" << m_myrank << "):" <<"gcsize    :" << v_gcsize << endl;
	PL_DBGOSH << "(my_rank:" << m_myrank << "):" <<"dx        :" << v_dx << endl;
	PL_DBGOSH << "(my_rank:" << m_myrank << "):" <<"gcell_min :"
		 << m_myproc.m_area.m_gcell_min << endl;
	PL_DBGOSH << "(my_rank:" << m_myrank << "):" <<"gcell_max :"
		 << m_myproc.m_area.m_gcell_max << endl;
#endif

	// 送信データ作成
	float send_buf[12];
	for (i = 0; i < 3; i++) {
		send_buf[i] = v_bpos[i];
	}
	for (i = 0; i < 3; i++) {
		send_buf[3+i] = v_bbsize[i];
	}
	for (i = 0; i < 3; i++) {
		send_buf[6+i] = v_gcsize[i];
	}
	for (i = 0; i < 3; i++) {
		send_buf[9+i] = v_dx[i];
	}

	// 受信領域確保
	float* recv_buf = new float[12 * m_numproc];

	// Allgather通信を行う
	if (MPI_Allgather(send_buf, 12, MPI_FLOAT, recv_buf, 12, MPI_FLOAT, comm) != MPI_SUCCESS) {
		PL_ERROSH << "[ERROR]MPIPolylib::init_parallel_info():MPI_Allgather "
				  << "faild." << endl;
		return PLSTAT_MPI_ERROR;
	}

	// 受信データの展開
	for (int irank = 0; irank < m_numproc; irank++) {
		// 自PE領域情報はスキップ
		if( irank == m_myrank ) continue;
		for (i = 0; i < 3; i++) {
			v_bpos[i] = recv_buf[i + 12*irank];
		}
		for (i = 0; i < 3; i++) {
			v_bbsize[i] = recv_buf[3 + i + 12*irank];
		} 
		for (i = 0; i < 3; i++) {
			v_gcsize[i] = recv_buf[6 + i + 12*irank];
		} 
		for (i = 0; i < 3; i++) {
			v_dx[i]	= recv_buf[9 + i + 12*irank];
		}
#ifdef DEBUG
		PL_DBGOSH << "(rank:" << irank << "):" <<"bpos  :" << v_bpos  << endl;
		PL_DBGOSH << "(rank:" << irank << "):" <<"bbsize:" << v_bbsize << endl;
		PL_DBGOSH << "(rank:" << irank << "):" <<"gcsize:" << v_gcsize << endl;
		PL_DBGOSH << "(rank:" << irank << "):" <<"dx    :" << v_dx << endl;
#endif
		ParallelInfo* proc = new (ParallelInfo);
		proc->m_comm = comm;
		proc->m_rank = irank;
		proc->m_area.m_bpos = v_bpos;
		proc->m_area.m_bbsize = v_bbsize;
		proc->m_area.m_gcsize = v_gcsize;
		proc->m_area.m_dx = v_dx;
		proc->m_area.m_gcell_min = v_bpos-( v_gcsize )*v_dx;
		proc->m_area.m_gcell_max = v_bpos+( v_bbsize+v_gcsize )*v_dx;
		proc->m_area.m_gcell_bbox.add(proc->m_area.m_gcell_min);
		proc->m_area.m_gcell_bbox.add(proc->m_area.m_gcell_max);

		// 全PE領域情報リストに追加
		m_other_procs.push_back(proc);

		// 自PE領域と隣接するPE領域情報はm_neibour_procsにも追加
		if( m_myproc.m_area.m_gcell_bbox.crossed(proc->m_area.m_gcell_bbox) ) {
			m_neibour_procs.push_back(proc);
#ifdef DEBUG
		PL_DBGOSH << m_myrank << ": " << "neighbour rank:" << proc->m_rank  << endl;
#endif
		}
	}
	// 受信領域あとしまつ
	delete[] recv_buf;

	return PLSTAT_OK;
}


#if 0 
// old tp version 
// public /////////////////////////////////////////////////////////////////////
POLYLIB_STAT
MPIPolylib::load_rank0(
	std::string config_filename
)
{
#ifdef DEBUG
	PL_DBGOSH << m_myrank << ": " << "MPIPolylib::load_rank0() in. " << endl;
#endif
	POLYLIB_STAT ret;
	string       config_contents;

	// for tp
	try {
		PolylibConfig base(config_filename);
		ret = make_group_tree(&base);
		if( ret != PLSTAT_OK ) return ret;
	}
	catch( POLYLIB_STAT e ){
		return e;
	}

	if( m_myrank == 0 ) {

	  // for tp
#if 0
	  // 設定ファイル読み込み。
	  if( (ret = Polylib::load_config_file( &config_contents, config_filename ))
	      != PLSTAT_OK ) {
	    PL_ERROSH << "[ERROR]MPIPolylib::load_rank0():Polylib::load_config() faild. returns:"
		      << PolylibStat2::String(ret) << endl;
	    return ret;
	  }

	  // グループ階層構造構築。
	  if( (ret = Polylib::make_group_tree( config_contents )) != PLSTAT_OK ) {
	    PL_ERROSH << "[ERROR]MPIPolylib::load_rank0():Polylib::make_group_tree() faild. returns:"
		      << PolylibStat2::String(ret) << endl;
	    return ret;
	  }

	  // 設定ファイルの内容を他PEへブロードキャストする
	  if( (ret = broadcast_config( config_contents )) != PLSTAT_OK ) {
	    PL_ERROSH << "[ERROR]MPIPolylib::load_rank0():broadcast_config() faild. returns:"
		      << PolylibStat2::String(ret) << endl;
	    return ret;
	  }
#endif

	  // ポリゴン情報を構築 (三角形IDファイルは不要なので、第二引数はダミー)
	  if( (ret = load_polygons(false, ID_BIN)) != PLSTAT_OK ) {
	    PL_ERROSH << "[ERROR]MPIPolylib::load_rank0():load_polygons() faild. returns:" << PolylibStat2::String(ret) << endl;
	    return ret;
	  }

	  // ポリゴン情報を他PEへ配信する。
	  if( (ret = send_polygons_to_all()) != PLSTAT_OK ) {
	    PL_ERROSH << "[ERROR]MPIPolylib::load_rank0():send_polygons_to_all() faild. returns:" << PolylibStat2::String(ret) << endl;
	    return ret;
	  }

	  // 他PE領域ポリゴン情報を削除して自領域分のみでデータ構造再構築
	  if( (ret = erase_outbounded_polygons()) != PLSTAT_OK ) {
	    PL_ERROSH << "[ERROR]MPIPolylib::load_rank0():erase_outbounded_polygons() faild. returns:" << PolylibStat2::String(ret) << endl;
	    return ret;
	  }
	}
	else {

	  // for tp
#if 0
	  // 設定ファイルの内容をrank0から受信する
	  if( (ret = broadcast_config_from_rank0()) != PLSTAT_OK ) {
	    PL_ERROSH << "[ERROR]MPIPolylib::load_rank0():broadcast_config_from_rank0() faild. returns:" << PolylibStat2::String(ret) << endl;
	    return ret;
	  }
#endif

	  // ポリゴン情報をrank0から受信する。
	  if( (ret = receive_polygons_from_rank0()) != PLSTAT_OK ) {
	    PL_ERROSH << "[ERROR]MPIPolylib::load_rank0():receive_polygons_from_rank0() faild. returns:" << PolylibStat2::String(ret) << endl;
	    return ret;
	  }
	}

	return PLSTAT_OK;

}
#endif



// new tp version  without  PolylibConfig
// public /////////////////////////////////////////////////////////////////////
POLYLIB_STAT
MPIPolylib::load_rank0(
	std::string config_filename
)
{
#ifdef DEBUG
	PL_DBGOSH << m_myrank << ": " << "MPIPolylib::load_rank0() in. " << endl;
#endif
	POLYLIB_STAT ret;
	string       config_contents;

	// for tp
	try {
	  //PolylibConfig base(config_filename);
		tp->read(config_filename);

		//ret = make_group_tree(&base);
		// only on rank0 ?
		ret = make_group_tree(tp);
		if( ret != PLSTAT_OK ) return ret;
	}
	catch( POLYLIB_STAT e ){
		return e;
	}

	if( m_myrank == 0 ) {

// for tp
#if 0
	  // 設定ファイル読み込み。
	  if( (ret = Polylib::load_config_file( &config_contents, config_filename ))
	      != PLSTAT_OK ) {
	    PL_ERROSH 
	      << "[ERROR]MPIPolylib::load_rank0():Polylib::load_config() faild."
	      <<" returns:" << PolylibStat2::String(ret) << endl;
	    return ret;
	  }

	  // グループ階層構造構築。
	  if( (ret = Polylib::make_group_tree( config_contents )) != PLSTAT_OK ) {
	    PL_ERROSH << "[ERROR]MPIPolylib::load_rank0():Polylib::make_group_tree() faild. returns:"
		      << PolylibStat2::String(ret) << endl;
	    return ret;
	  }

	  // 設定ファイルの内容を他PEへブロードキャストする
	  if( (ret = broadcast_config( config_contents )) != PLSTAT_OK ) {
	    PL_ERROSH << "[ERROR]MPIPolylib::load_rank0():broadcast_config() faild. returns:"
		      << PolylibStat2::String(ret) << endl;
	    return ret;
	  }
#endif

	  // ポリゴン情報を構築 (三角形IDファイルは不要なので、第二引数はダミー)
	  if( (ret = load_polygons(false, ID_BIN)) != PLSTAT_OK ) {
	    PL_ERROSH << "[ERROR]MPIPolylib::load_rank0():load_polygons() faild."
		      <<" returns:" << PolylibStat2::String(ret) << endl;
	    return ret;
	  }

	  // ポリゴン情報を他PEへ配信する。
	  if( (ret = send_polygons_to_all()) != PLSTAT_OK ) {
	    PL_ERROSH << "[ERROR]MPIPolylib::load_rank0():send_polygons_to_all()"
		      <<" faild. returns:" << PolylibStat2::String(ret) << endl;
	    return ret;
	  }

	  // 他PE領域ポリゴン情報を削除して自領域分のみでデータ構造再構築
	  if( (ret = erase_outbounded_polygons()) != PLSTAT_OK ) {
	    PL_ERROSH << "[ERROR]MPIPolylib::load_rank0():erase_outbounded_polygons()"
		      <<" faild. returns:" << PolylibStat2::String(ret) << endl;
	    return ret;
	  }
	} else {

// for tp
#if 0
	  // 設定ファイルの内容をrank0から受信する
	  if( (ret = broadcast_config_from_rank0()) != PLSTAT_OK ) {
	    PL_ERROSH << "[ERROR]MPIPolylib::load_rank0():broadcast_config_from_rank0()"
		      <<" faild. returns:" << PolylibStat2::String(ret) << endl;
	    return ret;
	  }
#endif

	    // ポリゴン情報をrank0から受信する。
	    if( (ret = receive_polygons_from_rank0()) != PLSTAT_OK ) {
	      PL_ERROSH << "[ERROR]MPIPolylib::load_rank0():receive_polygons_from_rank0()"
			<<" faild. returns:" << PolylibStat2::String(ret) << endl;
	      return ret;
	    }
	}
 
	return PLSTAT_OK;

}


// public /////////////////////////////////////////////////////////////////////
POLYLIB_STAT
MPIPolylib::load_parallel( 
	std::string config_filename,
	ID_FORMAT	id_format
)
{
#ifdef DEBUG
	PL_DBGOSH << "MPIPolylib::load_parallel() in. " << endl;
#endif
	POLYLIB_STAT ret;

	// 設定ファイル読み込み、グループ階層構造構築、ポリゴン情報構築。
	// MPIPolylib::save_parallel()で保存された設定ファイルであることが前提
	if( (ret = Polylib::load_with_idfile( config_filename, id_format )) != PLSTAT_OK ) {
		PL_ERROSH << "[ERROR]MPIPolylib::load_parallel():Polylib::load() faild. returns:" << PolylibStat2::String(ret) << endl;
		return ret;
	}

	return PLSTAT_OK;
}

#if 0 
// old version
// public /////////////////////////////////////////////////////////////////////
POLYLIB_STAT
MPIPolylib::save_rank0(
	std::string *p_config_filename,
	std::string stl_format,
	std::string extend
)
{
#ifdef DEBUG
	PL_DBGOSH << "MPIPolylib::save_rank0() in. " << endl;
#endif
	POLYLIB_STAT	ret;

	if( m_myrank == 0 ) {

		// 他rankからポリゴン情報を受信
		if( (ret = gather_polygons()) != PLSTAT_OK ) {
			PL_ERROSH << "[ERROR]MPIPolylib::save_rank0():gather_polygons() faild. returns:"
					  << PolylibStat2::String(ret) << endl;
			return ret;
		}

		// グループ階層構造、ポリゴン情報をファイルへ保存
// for tp
#if 0
		if( (ret = Polylib::save( p_config_filename, stl_format, extend )) != PLSTAT_OK ) {
			PL_ERROSH << "[ERROR]MPIPolylib::save_rank0():Polylib::save() faild. returns:"
					  << PolylibStat2::String(ret) << endl;
			return ret;
		}
#endif




	}
	else {

		// rank0へポリゴン情報を送信
		if( (ret = send_polygons_to_rank0()) != PLSTAT_OK ) {
			PL_ERROSH << "[ERROR]MPIPolylib::save_rank0():send_polygons_to_rank0() faild. returns:"
					  << PolylibStat2::String(ret) << endl;
			return ret;
		}
	}
	return PLSTAT_OK;
}
#endif

// new version 
// public /////////////////////////////////////////////////////////////////////
POLYLIB_STAT
MPIPolylib::save_rank0(
		       std::string *p_config_filename,
		       std::string stl_format,
		       std::string extend
		       )
{
#ifdef DEBUG
  PL_DBGOSH << "MPIPolylib::save_rank0() in. " << endl;
#endif
  POLYLIB_STAT	ret;

  if( m_myrank == 0 ) {
    // 他rankからポリゴン情報を受信
    if( (ret = gather_polygons()) != PLSTAT_OK ) {
      PL_ERROSH << "[ERROR]MPIPolylib::save_rank0():gather_polygons() faild. returns:"
		<< PolylibStat2::String(ret) << endl;
      return ret;
    }
    // グループ階層構造、ポリゴン情報をファイルへ保存

    // for tp

    if( (ret = Polylib::save( p_config_filename, stl_format, extend )) != PLSTAT_OK ) {
      PL_ERROSH << "[ERROR]MPIPolylib::save_rank0():Polylib::save() faild. returns:"
		<< PolylibStat2::String(ret) << endl;
      return ret;
    }
  }
  else {

    // rank0へポリゴン情報を送信
    if( (ret = send_polygons_to_rank0()) != PLSTAT_OK ) {
      PL_ERROSH << "[ERROR]MPIPolylib::save_rank0():send_polygons_to_rank0() faild."
		<<" returns:"<< PolylibStat2::String(ret) << endl;
      return ret;
    }
  }
  return PLSTAT_OK;
}


// public /////////////////////////////////////////////////////////////////////
POLYLIB_STAT
MPIPolylib::save_parallel(
	std::string *p_config_filename,
	std::string stl_format,
	std::string extend,
	ID_FORMAT	id_format
)
{
#ifdef DEBUG
	PL_DBGOSH << "MPIPolylib::save_parallel() in. " << endl;
#endif
	POLYLIB_STAT ret;

	// 各ランク毎に保存
	if( (ret = Polylib::save_with_rankno( p_config_filename, m_myrank, m_numproc-1, extend, stl_format, id_format)) != PLSTAT_OK ) {
		PL_ERROSH << "[ERROR]MPIPolylib::save_parallel():Polylib::save_with_rankno():failed. returns:" << PolylibStat2::String(ret) << endl;
		return ret;
	}
	return PLSTAT_OK;
}


// public ////////////////////////////////////////////////////////////////////
POLYLIB_STAT
MPIPolylib::move(
	PolylibMoveParams &params
)
{
#ifdef DEBUG
	PL_DBGOSH << "MPIPolylib::move() in. " << endl;
#endif
	POLYLIB_STAT ret;
	vector<PolygonGroup*>::iterator group_itr;
	PolygonGroup *p_pg;

	// 各ポリゴングループのmove()を実行
	// 全グループに対して
	for (group_itr = m_pg_list.begin(); group_itr != m_pg_list.end(); group_itr++) {
			p_pg = (*group_itr);

		// 移動する可能性のあるポリゴングループのみ対象
		if( p_pg->get_movable() ) {

			// move実行前から隣接PE領域に懸かっている三角形をmigrate対象除外リストに載せる
			if( (ret = select_excluded_trias( p_pg )) != PLSTAT_OK ) {
				PL_ERROSH << "[ERROR]MPIPolylib::move():select_exclude_trias() failed. returns:" << PolylibStat2::String(ret) << endl;
				return ret;
			}

			// move実行
			if( (ret = p_pg->move( params )) != PLSTAT_OK ) {
				PL_ERROSH << "[ERROR]MPIPolylib::move():(*group_itr)->move() failed. returns:" << PolylibStat2::String(ret) << endl;
				return ret;
			}

			// KD木を再構築 (三角形同士の位置関係が変化したため、再構築が必要)
			if( (ret = p_pg->rebuild_polygons()) != PLSTAT_OK ) {
				PL_ERROSH << "[ERROR]MPIPolylib::move():(*group_itr)->rebuild_polygons() failed. returns:" << PolylibStat2::String(ret) << endl;
				return ret;
			}
		}
	}
	return PLSTAT_OK;
}


// public /////////////////////////////////////////////////////////////////////
POLYLIB_STAT
MPIPolylib::migrate(
)
{
#ifdef DEBUG
	PL_DBGOSH << "MPIPolylib::migrate() in. " << endl;
#endif
	POLYLIB_STAT ret;
	unsigned int i, j;
	vector<PolygonGroup*>::iterator group_itr;
	PolygonGroup *p_pg;
	vector<ParallelInfo*>::iterator procs_itr;

	vector<PrivateTriangle*> const *p_trias;
	vector<int>   send_num_trias;
	int          *p_send_num_trias_array;
	vector<int>   send_tria_ids;
	int          *p_send_tria_ids_array;
	vector<float> send_trias;
	float        *p_send_trias_array;

	// 送信バッファ領域をMPI_WaitAll()後に纏めてdeleteするための
	// 配列アドレス記憶用ベクタ
	vector<int*>  send_num_trias_bufs;
	vector<int*>  send_tria_ids_bufs;
	vector<float*> send_trias_bufs;

	// 送信用MPI_Reqeust配列を確保
	MPI_Request *mpi_reqs = new MPI_Request[ m_neibour_procs.size() * 3 ]; // 隣接PEごとに3回Isendする
	MPI_Status  *mpi_stats = new MPI_Status[ m_neibour_procs.size() * 3 ];
	int reqs_pos = 0;

	//隣接PEごとに移動三角形情報を送信
	for (procs_itr = m_neibour_procs.begin(); procs_itr != m_neibour_procs.end(); procs_itr++) {

		// 送信用一時データ初期化
		send_num_trias.clear();
		send_trias.clear();
		send_tria_ids.clear();
		
		// 全ポリゴングループに対して
		for( group_itr=m_pg_list.begin(); group_itr!=m_pg_list.end(); group_itr++ ) {
			p_pg = (*group_itr);
			p_trias = NULL;

			// ポリゴン情報を持つグループだけ
			//if( p_pg->get_triangles() != NULL
			// && p_pg->get_triangles()->size() != 0 ) {

			// 移動する可能性のあるポリゴングループのみ対象
#ifdef DEBUG
PL_DBGOSH << "pg_name:" << p_pg->get_name() << " movable:" << p_pg->get_movable() << endl;
#endif
			if( p_pg->get_movable() ) {

				// 当該隣接PE領域への移動除外三角形IDリストを取得
				map< int, vector<int> >::iterator const itr =
					(*procs_itr)->m_exclusion_map.find( p_pg->get_internal_id() );

				// 当該隣接PE領域内にある移動フラグONの三角形を取得
				p_trias = p_pg->search_outbounded(
					(*procs_itr)->m_area.m_gcell_bbox, &((*itr).second) );
			}

			// グループIDと当該グループの三角形数の対を送信データに追加
			pack_num_trias( &send_num_trias, p_pg->get_internal_id(), p_trias );

			// 三角形情報を送信データに追加
			pack_trias( &send_trias, p_trias );

			// 三角形ID情報を送信データに追加
			pack_tria_ids( &send_tria_ids, p_trias );

			// search結果の後始末
			if( p_trias ) delete p_trias;
		}

		//-----  送信データをシリアライズ
		// 送信データ初期化
		p_send_num_trias_array = NULL;
		p_send_tria_ids_array  = NULL;
		p_send_trias_array     = NULL;

		// グループID,グループ毎三角形数リスト
		if( send_num_trias.size() > 0 ) {
			p_send_num_trias_array = new int[ send_num_trias.size() ];
		}
		for( i=0; i<send_num_trias.size(); i++ ) {
			p_send_num_trias_array[i] = send_num_trias[i];
		}
		// 三角形IDリスト
		if( send_tria_ids.size() > 0 ) {
			p_send_tria_ids_array  = new int[ send_tria_ids.size() ];
		}
		for( i=0; i<send_tria_ids.size(); i++ ) {
			p_send_tria_ids_array[i] = send_tria_ids[i];
		}
		// 三角形頂点リスト
		if( send_trias.size() > 0 ) {
			p_send_trias_array = new float[ send_trias.size() ];
		}
		for( i=0; i<send_trias.size(); i++ ) {
			p_send_trias_array[i] = send_trias[i];
		}

		// 送信データの先頭アドレスを記憶（MPI_Wait後にdeleteするため）
		send_num_trias_bufs.push_back( p_send_num_trias_array );
		send_tria_ids_bufs.push_back( p_send_tria_ids_array );
		send_trias_bufs.push_back( p_send_trias_array );

		// 当該PEへ非同期送信 (MPI_Wait()は後でまとめて行う)
#ifdef DEBUG
		PL_DBGOSH << "sending polygons rank:" << m_myrank <<  "->rank:"
				  << (*procs_itr)->m_rank << " ";
		for( i=0; i< send_num_trias.size(); i+=2 ) {
			PL_DBGOS << "(gid:" << send_num_trias[i] 
					 << ",num_tria:" << send_num_trias[i+1] << ")";
		}
		PL_DBGOS << endl;
#endif
		if (MPI_Isend( p_send_num_trias_array, send_num_trias.size(),
					MPI_INT, (*procs_itr)->m_rank, MPITAG_NUM_TRIAS,
					m_mycomm, &mpi_reqs[reqs_pos++] ) != MPI_SUCCESS) {
			PL_ERROSH << "[ERROR]MPIPolylib::migrate():MPI_Isend,"
					  << "MPITAG_NUM_TRIAS faild." << endl;
			return PLSTAT_MPI_ERROR;
		}
		if (MPI_Isend( p_send_tria_ids_array,  send_tria_ids.size(),
					MPI_INT, (*procs_itr)->m_rank, MPITAG_TRIA_IDS,
					m_mycomm, &mpi_reqs[reqs_pos++] ) != MPI_SUCCESS) {
			PL_ERROSH << "[ERROR]MPIPolylib::migrate():MPI_Isend,"
					  << " MPITAG_TRIA_IDS faild." << endl;
			return PLSTAT_MPI_ERROR;
		}
		if (MPI_Isend( p_send_trias_array,     send_trias.size(),
					MPI_FLOAT, (*procs_itr)->m_rank, MPITAG_TRIAS,
					m_mycomm, &mpi_reqs[reqs_pos++] ) != MPI_SUCCESS) {
			PL_ERROSH << "[ERROR]MPIPolylib::migrate():MPI_Isend,"
					  << " MPITAG_TRIAS faild." << endl;
			return PLSTAT_MPI_ERROR;
		}
	}

	//隣接PEごとに移動三角形情報を受信
	for (procs_itr = m_neibour_procs.begin(); procs_itr != m_neibour_procs.end(); procs_itr++) {
		int pos_id, pos_tria;
		MPI_Request mpi_req;
		MPI_Status  mpi_stat;

		// グループIDとグループ毎三角形数の対を非同期受信
		// グループ情報は各rank共有しているのでグループ数は予め分かっている
		int *p_intarray = new int[ m_pg_list.size()*2 ];
		if (MPI_Irecv( p_intarray, m_pg_list.size()*2, MPI_INT, (*procs_itr)->m_rank,
					MPITAG_NUM_TRIAS, m_mycomm, &mpi_req ) != MPI_SUCCESS) {
			PL_ERROSH << "[ERROR]MPIPolylib::migrate():MPI_Irecv,"
					  << ",MPITAG_NUM_TRIAS faild." << endl;
			return PLSTAT_MPI_ERROR;
		}
		if (MPI_Wait( &mpi_req, &mpi_stat ) != MPI_SUCCESS) {
			PL_ERROSH << "[ERROR]MPIPolylib::migrate():MPI_Wait,"
					  << "MPITAG_NUM_TRIAS faild." << endl;
			return PLSTAT_MPI_ERROR;
		}
#ifdef DEBUG
		PL_DBGOSH << "receiving polygons rank:" << (*procs_itr)->m_rank
				  <<  "->rank:" << m_myrank << " ";
		for( i=0; i< m_pg_list.size()*2-1; i+=2 ) {
			PL_DBGOS << "(gid:" << p_intarray[i] 
					 << ",num_tria:" << p_intarray[i+1] << ")";
		}
		PL_DBGOS << endl;
#endif

		// 受信三角形数を算出
		int total_tria_num = 0;
		for( i=1; i<m_pg_list.size() * 2; i+=2 ){
			total_tria_num += p_intarray[i];
		}

		// 三角形IDリストを非同期受信
		int *p_idarray = new int[ total_tria_num ];
		if (MPI_Irecv( p_idarray,  total_tria_num, MPI_INT, (*procs_itr)->m_rank,
					 MPITAG_TRIA_IDS, m_mycomm, &mpi_req ) != MPI_SUCCESS) {
			PL_ERROSH << "[ERROR]MPIPolylib::migrate():MPI_Irecv,"
					  << "MPI_INT faild." << endl;
			return PLSTAT_MPI_ERROR;
		}
		if (MPI_Wait( &mpi_req, &mpi_stat ) != MPI_SUCCESS) {
			PL_ERROSH << "[ERROR]MPIPolylib::migrate():MPI_Wait,"
					  << "MPI_INT faild." << endl;
			return PLSTAT_MPI_ERROR;
		}

		// 三角形リストを非同期受信
		float *p_triaarray = new float[ total_tria_num*3*3 ];
		if (MPI_Irecv( p_triaarray, total_tria_num*3*3, MPI_FLOAT, (*procs_itr)->m_rank,
					MPITAG_TRIAS, m_mycomm, &mpi_req ) != MPI_SUCCESS) {
			PL_ERROSH << "[ERROR]MPIPolylib::migrate():MPI_Irecv,MPI_FLOAT"
					  << " faild." << endl;
			return PLSTAT_MPI_ERROR;
		}
		if (MPI_Wait( &mpi_req, &mpi_stat ) != MPI_SUCCESS) {
			PL_ERROSH << "[ERROR]MPIPolylib::migrate():MPI_Wait,MPI_FLOAT"
					  << " faild." << endl;
			return PLSTAT_MPI_ERROR;
		}

		// 各ポリゴングループに対して三角形情報を追加
		pos_id = 0;
		pos_tria = 0;
		for( i=0; i<m_pg_list.size()*2-1; i+=2 ){

			// ポリゴングループID
			int pg_id = p_intarray[i];

			// 当該ポリゴングループの三角形数
			unsigned int num_trias = p_intarray[i+1];

			// グループIDのポリゴングループインスタンス取得
			PolygonGroup* p_pg = get_group( pg_id );
			if( p_pg == NULL ) {
				PL_ERROSH << "[ERROR]MPIPolylib::migrate():invalid pg_id:"
					  	<< pg_id << endl;
				return PLSTAT_NG;
			}

			// PrivateTriangleのベクタ - 受信データ配列からベクタへの変換用
			vector<PrivateTriangle*> tria_vec;

			// ベクタに受信データ内容をコピー
			for( j=0; j<num_trias; j++ ) {
				tria_vec.push_back(
						new PrivateTriangle(&p_triaarray[pos_tria], p_idarray[pos_id]) );
				pos_id++;
				pos_tria+=9;
			}

			// ポリゴングループに三角形リストを追加
			if( (ret = p_pg->add_triangles( &tria_vec )) != PLSTAT_OK ) {
				PL_ERROSH << "[ERROR]MPIPolylib::migrate():p_pg->add_triangles() failed. returns:"
						  << PolylibStat2::String(ret) << endl;
				return ret;
			}

			// ベクタの内容あとしまつ
			for( j=0; j<num_trias; j++ ) {
				delete tria_vec.at(j);
			}
		}

		// 受信領域あとしまつ
		delete[] p_intarray;
		delete[] p_idarray;
		delete[] p_triaarray;
	}

	// MPI_Isend()を纏めてアンロック
	if (MPI_Waitall( m_neibour_procs.size()*3, mpi_reqs, mpi_stats ) != MPI_SUCCESS) {
		PL_ERROSH << "[ERROR]MPIPolylib::migrate():MPI_Waitall failed." << endl;
		return PLSTAT_MPI_ERROR;
	}

	// 送信データ領域をdelete
	for( i=0; i<send_num_trias_bufs.size(); i++ ) {
		delete[] send_num_trias_bufs.at(i);
	}
	for( i=0; i<send_tria_ids_bufs.size(); i++ ) {
		delete[] send_tria_ids_bufs.at(i);
	}
	for( i=0; i<send_trias_bufs.size(); i++ ) {
		delete[] send_trias_bufs.at(i);
	}

	// 移動してきた三角形を含めたKD木を再構築
	for (group_itr = m_pg_list.begin(); group_itr != m_pg_list.end(); group_itr++) {
		p_pg = (*group_itr);

		// 移動可能グループだけ
		if( p_pg->get_movable() && p_pg->get_triangles() != NULL ) {

			// KD木を再構築
			if( (ret=p_pg->rebuild_polygons()) != PLSTAT_OK ) {
				PL_ERROSH << "[ERROR]MPIPolylib::migrate():p_pg->rebuild_polygons() failed. returns:"
						  << PolylibStat2::String(ret) << endl;
				return ret;
			}
		}
	}
	
	// 自PE領域外ポリゴン情報を消去
	if( erase_outbounded_polygons() != PLSTAT_OK ) {
		PL_ERROSH << "[ERROR]MPIPolylib::migrate():rebuild_polygons() failed." << endl;
	}

	// 後始末 2010.08.24
	delete[] mpi_reqs;
	delete[] mpi_stats;

#ifdef DEBUG
	PL_DBGOSH << "MPIPolylib::migrate() out normaly." << endl;
#endif
	return PLSTAT_OK;
}


// public /////////////////////////////////////////////////////////////////////
ParallelInfo* MPIPolylib::get_proc(int rank)
{
	vector<ParallelInfo*>::iterator itr;
	itr = m_other_procs.begin();
	for (; itr != m_other_procs.end(); itr++) {
		if ((*itr)->m_rank == rank) {
			return(*itr);
		}
	}
	return NULL;
}


// public /////////////////////////////////////////////////////////////////////
unsigned int MPIPolylib::used_memory_size()
{
	unsigned int								size;
	map< int, vector<int> >::iterator	ex;
	vector<ParallelInfo *>::iterator			pi;

	// Polylibクラスが管理している領域を取得
	size = Polylib::used_memory_size();

	// 自PE担当領域情報
	size += sizeof(ParallelInfo);
	size += m_myproc.m_exclusion_map.size() * (sizeof(int)+sizeof(vector<int>));
	for (ex = m_myproc.m_exclusion_map.begin(); 
							ex != m_myproc.m_exclusion_map.end(); ex++) {
		size += ex->second.size() * sizeof(int);
	}

	// 自PEを除く全PE担当領域情報リスト
	size += sizeof(vector<ParallelInfo *>);
	for (pi = m_other_procs.begin(); pi != m_other_procs.end(); pi++) {
		size += sizeof(ParallelInfo);
		size += (*pi)->m_exclusion_map.size() * (sizeof(int)+sizeof(vector<int>));
		for (ex = (*pi)->m_exclusion_map.begin(); 
									ex != (*pi)->m_exclusion_map.end(); ex++) {
			size += ex->second.size() * sizeof(int);
		}
	}

	// 隣接PE担当領域情報リスト
	size += sizeof(vector<ParallelInfo *>);
	size += m_neibour_procs.size() * sizeof(vector<ParallelInfo *>);

	// 自プロセスのランク数、全プロセス数
	size += sizeof(int) * 2;

	// 自プロセスが利用するコミュニケーター
	size += sizeof(MPI_Comm);

	return size;
}


// protected //////////////////////////////////////////////////////////////////
MPIPolylib::MPIPolylib() : Polylib()
{
}


// protected //////////////////////////////////////////////////////////////////
MPIPolylib::~MPIPolylib()
{
	vector<ParallelInfo*>::iterator itr;

	// 全PE領域情報リストを消去
	for ( itr=m_other_procs.begin(); itr != m_other_procs.end(); itr++ ) {
		delete *itr;
	}
	// 隣接PE領域情報リストについては、m_other_procsの内容と同じインスタンスを保持
	// するため、上記の処理で消去済。
}


// protected //////////////////////////////////////////////////////////////////
POLYLIB_STAT
MPIPolylib::broadcast_config(
	string config_contents
)
{
#ifdef DEBUG
	PL_DBGOSH << "MPIPolylib::broadcast_config() in. " << endl;
#endif

	// 文字数を送信(NULL文字も文字数にカウント）
	int data = config_contents.size()+1;
	if (MPI_Bcast( &data, 1, MPI_INT, 0, m_mycomm ) != MPI_SUCCESS) {
		PL_ERROSH << "[ERROR]MPIPolylib::broadcast_config():MPI_Bcast,"
				 << "MPI_INT faild." << endl;
		return PLSTAT_MPI_ERROR;
	}

	// ファイル内容文字列の送信
	if (MPI_Bcast( (void*)(config_contents.c_str()), config_contents.size()+1,
				MPI_CHAR, 0, m_mycomm ) != MPI_SUCCESS) {
		PL_ERROSH << "[ERROR]MPIPolylib::broadcast_config():MPI_Bcast,"
				  << "MPI_CHAR faild." << endl;
		return PLSTAT_MPI_ERROR;
	}

	return PLSTAT_OK;
}


// protected //////////////////////////////////////////////////////////////////
POLYLIB_STAT
MPIPolylib::send_polygons_to_all(
)
{
#ifdef DEBUG
	PL_DBGOSH << "MPIPolylib::send_polygons_to_all() in. " << endl;
#endif
	unsigned int i;
	vector<ParallelInfo*>::iterator	proc_itr;
	vector<PolygonGroup*>::iterator group_itr;
	PolygonGroup *p_pg;
	vector<PrivateTriangle*> const *p_trias;

	vector<int>   send_num_trias;
	int          *p_send_num_trias_array;
	vector<int>   send_tria_ids;
	int          *p_send_tria_ids_array;
	vector<float> send_trias;
	float        *p_send_trias_array;

	// 全PEに対して
	for (proc_itr = m_other_procs.begin(); proc_itr != m_other_procs.end(); proc_itr++) {

		// 送信用一時データ初期化
		p_trias = NULL;
		send_num_trias.clear();
		send_trias.clear();
		send_tria_ids.clear();
		
		// 全グループに対して
		for (group_itr = m_pg_list.begin(); group_itr != m_pg_list.end(); group_itr++) {
			p_pg = (*group_itr);
			p_trias = NULL;

			// ポリゴン情報を持つグループだけ
			if( p_pg->get_triangles() != NULL && p_pg->get_triangles()->size() != 0 ) {

				// 当該PE領域内に一部でも含まれるポリゴンを検索
				p_trias = p_pg->search( &((*proc_itr)->m_area.m_gcell_bbox), false );
			}

			// グループIDと当該グループの三角形数の対を送信データに追加
			pack_num_trias( &send_num_trias, p_pg->get_internal_id(), p_trias );

			// 三角形情報を送信データに追加
			pack_trias( &send_trias, p_trias );

			// 三角形ID情報を送信データに追加
			pack_tria_ids( &send_tria_ids, p_trias );

			// search結果の後始末
			if( p_trias ) delete p_trias;
		}

		//----  送信データをシリアライズ
		// 送信データ配列初期化
		p_send_num_trias_array = NULL;
		p_send_tria_ids_array  = NULL;
		p_send_trias_array     = NULL;
 
		// グループID,グループ毎三角形数リスト
		if( send_num_trias.size() > 0 ) {
			p_send_num_trias_array = new int[ send_num_trias.size() ];
		}
		for( i=0; i<send_num_trias.size(); i++ ) {
			p_send_num_trias_array[i] = send_num_trias[i];
		}
		// 三角形IDリスト
		if( send_tria_ids.size() > 0 ) {
			p_send_tria_ids_array  = new int[ send_tria_ids.size() ];
		}
		for( i=0; i<send_tria_ids.size(); i++ ) {
			p_send_tria_ids_array[i] = send_tria_ids[i];
		}
		// 三角形頂点座標リスト
		if( send_trias.size() > 0 ) {
			p_send_trias_array = new float[ send_trias.size() ];
		}
		for( i=0; i<send_trias.size(); i++ ) {
			p_send_trias_array[i] = send_trias[i];
		}

		// 当該PEへ送信
#ifdef DEBUG
		PL_DBGOSH << "sending polygons rank:0->rank:" << (*proc_itr)->m_rank << " ";
		for( i=0; i< send_num_trias.size(); i+=2 ) {
			PL_DBGOS << "(gid:" << send_num_trias[i] 
					 << ",num_tria:" << send_num_trias[i+1] << ")";
		}
		PL_DBGOS << endl;
#endif
		if (MPI_Send( p_send_num_trias_array, send_num_trias.size(), MPI_INT,
			(*proc_itr)->m_rank, MPITAG_NUM_TRIAS, m_mycomm ) != MPI_SUCCESS) {
			PL_ERROSH << "[ERROR]MPIPolylib::send_polygons_to_all():MPI_Send,"
					  << "MPITAG_NUM_TRIAS faild." << endl;
			return PLSTAT_MPI_ERROR;
		}
		if (MPI_Send( p_send_tria_ids_array,  send_tria_ids.size(),  MPI_INT,
			(*proc_itr)->m_rank, MPITAG_TRIA_IDS,  m_mycomm ) != MPI_SUCCESS) {
			PL_ERROSH << "[ERROR]MPIPolylib::send_polygons_to_all():MPI_Send,"
					  << "MPITAG_TRIA_IDS faild." << endl;
			return PLSTAT_MPI_ERROR;
		}
		if (MPI_Send( p_send_trias_array,     send_trias.size(),     MPI_FLOAT,
			(*proc_itr)->m_rank, MPITAG_TRIAS,     m_mycomm ) != MPI_SUCCESS) {
			PL_ERROSH << "[ERROR]MPIPolylib::send_polygons_to_all():MPI_Send,"
					  << "MPITAG_TRIAS faild." << endl;
			return PLSTAT_MPI_ERROR;
		}

		// 送信データの後始末
		if( p_send_num_trias_array ) delete[] p_send_num_trias_array;
		if( p_send_tria_ids_array )  delete[] p_send_tria_ids_array;
		if( p_send_trias_array )     delete[] p_send_trias_array;
	}
	return PLSTAT_OK;
}


// protected //////////////////////////////////////////////////////////////////
POLYLIB_STAT
MPIPolylib::pack_num_trias(
	vector<int>* p_vec,
	int group_id,
	const vector<PrivateTriangle*>* p_trias
)
{
#ifdef DEBUG
	PL_DBGOSH << "MPIPolylib::pack_num_trias() in. " << endl;
#endif
	// 出力配列に、グループID、グループ内三角形数、の順に追加
	int num = 0;
	if( p_trias ) num = p_trias->size();
	p_vec->push_back( group_id );
	p_vec->push_back( num );

	return PLSTAT_OK;
}


// protected //////////////////////////////////////////////////////////////////
POLYLIB_STAT
MPIPolylib::pack_trias(
	vector<float>* p_vec,
	const vector<PrivateTriangle*>* p_trias
)
{
#ifdef DEBUG
	PL_DBGOSH << "MPIPolylib::pack_trias() in. " << endl;
#endif
	if( p_trias == NULL ) return PLSTAT_OK;

	// 出力配列に、三角形頂点座標を順に追加
	for( unsigned int i=0; i<p_trias->size(); i++ ) {
		for( unsigned int j=0; j<3; j++ ) {
			for( unsigned int k=0; k<3; k++ ) {
				p_vec->push_back( p_trias->at(i)->get_vertex()[j][k] );
			}
		}
	}
	return PLSTAT_OK;
}


// protected //////////////////////////////////////////////////////////////////
POLYLIB_STAT
MPIPolylib::pack_tria_ids(
	vector<int>* p_vec,
	const vector<PrivateTriangle*>* p_trias
)
{
#ifdef DEBUG
	PL_DBGOSH << "MPIPolylib::pack_tria_ids() in. " << endl;
#endif
	if( p_trias == NULL ) return PLSTAT_OK;

	// 出力配列に、三角形IDを順に追加
	for( unsigned int i=0; i<p_trias->size(); i++ ) {
		p_vec->push_back( p_trias->at(i)->get_id() );
	}
	return PLSTAT_OK;
}


// protected //////////////////////////////////////////////////////////////////
POLYLIB_STAT
MPIPolylib::erase_outbounded_polygons(
)
{
#ifdef DEBUG
	PL_DBGOSH << "MPIPolylib::erase_outbounded_polygons() in. " << endl;
#endif
	POLYLIB_STAT ret;
	unsigned int i;
	vector<PolygonGroup*>::iterator group_itr;
	vector<PrivateTriangle*> const *p_trias;
	vector<PrivateTriangle*>  copy_trias;
	PolygonGroup *p_pg;

	// 各ポリゴングループのポリゴン情報を自領域分のみで再構築
	// 全グループに対して
	for (group_itr = m_pg_list.begin(); group_itr != m_pg_list.end(); group_itr++) {
		p_pg = (*group_itr);

		// ポリゴン情報を持つグループだけ
		if( p_pg->get_triangles() != NULL && p_pg->get_triangles()->size() != 0 ) {

			// 自領域内に一部でも含まれるポリゴンを検索
			p_trias = p_pg->search( &(m_myproc.m_area.m_gcell_bbox), false );

			// 検索結果のディープコピーを作成
			copy_trias.clear();
			if( p_trias ) {
				for( i=0; i<p_trias->size(); i++ ) {
					copy_trias.push_back( new PrivateTriangle(*(p_trias->at(i))) );
				}
			}

			// 検索結果でポリゴン情報を再構築
			if( (ret = p_pg->init( &copy_trias, true )) != PLSTAT_OK ) {
				PL_ERROSH << "[ERROR]MPIPolylib::erase_outbounded_polygons():p_pg->init() failed. returns:"
						  << PolylibStat2::String(ret) << endl;
				return ret;
			}

			// search結果の後始末
			if( p_trias ) delete p_trias;
			for( i=0; i<copy_trias.size(); i++ ) {
				delete copy_trias.at(i);
			}
		}
	}
	return PLSTAT_OK;
}


// protected //////////////////////////////////////////////////////////////////
POLYLIB_STAT
MPIPolylib::broadcast_config_from_rank0(
)
{
#ifdef DEBUG
	PL_DBGOSH << "MPIPolylib::broadcast_config_from_rank0() in. " << endl;
#endif
	POLYLIB_STAT ret;
	int str_len;
	char *p_str;
	string config_str;

	// 文字数をb_castで受信(NULL文字も文字数に含まれる)
	if (MPI_Bcast( &str_len, 1, MPI_INT, 0, m_mycomm ) != MPI_SUCCESS) {
		PL_ERROSH << "[ERROR]MPIPolylib::broadcast_config_from_rank0()"
				  << ":MPI_Bcast,MPI_INT faild." << endl;
		return PLSTAT_MPI_ERROR;
	}

	// XML文字列をb_castで受信
	p_str = new char[str_len];
	if (MPI_Bcast( p_str, str_len, MPI_CHAR, 0, m_mycomm ) != MPI_SUCCESS) {
		PL_ERROSH << "[ERROR]MPIPolylib::broadcast_config_from_rank0()"
				  << ":MPI_Bcast,MPI_CHAR faild." << endl;
		return PLSTAT_MPI_ERROR;
	}
#ifdef DEBUG
	PL_DBGOSH << "received config-string:" << p_str << endl;
#endif

	// XML文字列からグループ階層構造を構築
	config_str = p_str;
	if( (ret = make_group_tree( config_str )) != PLSTAT_OK ) {
		PL_ERROSH << "[ERROR]MPIPolylib::broadcast_config_from_rank0():make_group_tree() faild. returns:" << ret << endl;
		return ret;
	}

	return PLSTAT_OK;
}


// protected //////////////////////////////////////////////////////////////////
POLYLIB_STAT
MPIPolylib::receive_polygons_from_rank0(
)
{
#ifdef DEBUG
	PL_DBGOSH << "MPIPolylib::receive_polygons_from_rank0() in. " << endl;
#endif

	unsigned int i, j;
	unsigned int pos_id, pos_tria;
	MPI_Status mpi_stat;

	// グループIDとグループ毎三角形数の対をrank0から受信
	// グループ情報は配信済みなので、グループ数は予め分かっている
	int *p_intarray = new int[ m_pg_list.size() * 2 ];
	if (MPI_Recv( p_intarray, m_pg_list.size() * 2, MPI_INT, 0, 
				MPITAG_NUM_TRIAS, m_mycomm, &mpi_stat ) != MPI_SUCCESS) {
		PL_ERROSH << "[ERROR]MPIPolylib::receive_polygons_from_rank0()"
				  << ":MPI_Recv,MPITAG_NUM_TRIAS faild." << endl;
		return PLSTAT_MPI_ERROR;
	}
#ifdef DEBUG
	PL_DBGOSH << "    pintarray:(";
	for( i=0; i< m_pg_list.size()*2; i++ ) {
		PL_DBGOS << p_intarray[i] << ",";
	}
	PL_DBGOS << ")" << endl;
#endif

	// 自領域の全三角形数を算出
	unsigned int total_tria_num = 0;
	for( i=1; i<m_pg_list.size() * 2; i+=2 ){
		total_tria_num += p_intarray[i];
	}

	// 三角形IDリストをrank0から受信
	int *p_idarray = new int[ total_tria_num ];
	if (MPI_Recv( p_idarray,  total_tria_num, MPI_INT, 0, 
				MPITAG_TRIA_IDS, m_mycomm, &mpi_stat ) != MPI_SUCCESS) {
		PL_ERROSH << "[ERROR]MPIPolylib::receive_polygons_from_rank0()"
				  << ":MPI_Recv,MPITAG_TRIA_IDS faild." << endl;
		return PLSTAT_MPI_ERROR;
	}
#ifdef DEBUG
	PL_DBGOSH << "    pidarray:(";
	for( i=0; i<total_tria_num; i++ ) {
		PL_DBGOS << p_idarray[i] << ",";
	}
	PL_DBGOS << ")" << endl;
#endif

	// 三角形リストをrank0から受信
	float *p_triaarray = new float[ total_tria_num * 3 * 3 ];
	if (MPI_Recv( p_triaarray, total_tria_num * 3 * 3, MPI_FLOAT, 0, 
				MPITAG_TRIAS, m_mycomm, &mpi_stat ) != MPI_SUCCESS) {
		PL_ERROSH << "[ERROR]MPIPolylib::receive_polygons_from_rank0()"
				  << ":MPI_Recv,MPITAG_TRIAS faild." << endl;
		return PLSTAT_MPI_ERROR;
	}
#ifdef DEBUG
	PL_DBGOSH << "    ptriaarray:(";
	for( i=0; i<total_tria_num*3*3; i++ ) {
		PL_DBGOS << p_triaarray[i] << ",";
	}
	PL_DBGOS << ")" << endl;
#endif

	// 各ポリゴングループに対して三角形情報を設定＆KD木構築
	pos_id = 0;
	pos_tria = 0;
	for( i=0; i<m_pg_list.size()*2-1; i+=2 ){	// 偶数番目の値を処理

		// ポリゴングループID
		int pg_id = p_intarray[i];

		// 当該ポリゴングループの三角形数
		unsigned int num_trias = p_intarray[i+1];

		// グループIDのポリゴングループインスタンス取得
		PolygonGroup* p_pg = get_group( pg_id );
		if( p_pg == NULL ) {
			PL_ERROSH << "[ERROR]MPIPolylib::receive_polygons_from_rank0():invalid pg_id:"
					  << pg_id << endl;
			return PLSTAT_NG;
		}

		// PrivateTriangleのベクタ生成
		vector<PrivateTriangle*> tria_vec;

		// ベクタに受信データ内容を設定
		for( j=0; j<num_trias; j++ ) {
			tria_vec.push_back( new PrivateTriangle(&p_triaarray[pos_tria], p_idarray[pos_id]) );
			pos_id++;
			pos_tria+=9;
		}

		// ポリゴングループに三角形リストを設定、KD木構築
		if( p_pg->init( &tria_vec, true ) != PLSTAT_OK ) {
			PL_ERROSH << "[ERROR]MPIPolylib::receive_polygons_from_rank0():p_pg->init() failed:" << endl;
			return PLSTAT_NG;
		}

		// ベクタの内容あとしまつ
		for( j=0; j<num_trias; j++ ) {
			delete tria_vec.at(j);
		}
	}

	// 受信領域あとしまつ
	delete[] p_intarray;
	delete[] p_idarray;
	delete[] p_triaarray;

#ifdef DEBUG
	PL_DBGOSH << "MPIPolylib::receive_polygons_from_rank0() out. " << endl;
#endif
	return PLSTAT_OK;
}


// protected //////////////////////////////////////////////////////////////////
POLYLIB_STAT
MPIPolylib::gather_polygons(
)
{
#ifdef DEBUG
	PL_DBGOSH << "MPIPolylib::gather_polygons() in. " << endl;
#endif

	POLYLIB_STAT ret;
	unsigned int i, j;
	int rank;
	unsigned int pos_id, pos_tria;
	MPI_Status mpi_stat;

	// 全rankからポリゴン情報を受信
	for( rank=1; rank<m_numproc; rank++ ) {

		// グループIDとグループ毎三角形数の対を受信
		// グループ情報は全rank共通なので、グループ数は予め分かっている
		int *p_intarray = new int[ m_pg_list.size() * 2 ];
		if (MPI_Recv( p_intarray, m_pg_list.size() * 2, MPI_INT, rank, 
					MPITAG_NUM_TRIAS, m_mycomm, &mpi_stat ) != MPI_SUCCESS) {
			PL_ERROSH << "[ERROR]MPIPolylib::gather_polygons()"
					  << ":MPI_Recv,MPITAG_NUM_TRIAS faild.:rank=" << rank 
					  << endl;
			return PLSTAT_MPI_ERROR;
		}

		// 受信する全三角形数を算出
		unsigned int total_tria_num = 0;
		for( i=1; i<m_pg_list.size() * 2; i+=2 ){
			total_tria_num += p_intarray[i];
		}

		// 三角形IDリストを受信
		int *p_idarray = new int[ total_tria_num ];
		if (MPI_Recv( p_idarray,  total_tria_num, MPI_INT, rank, 
					MPITAG_TRIA_IDS, m_mycomm, &mpi_stat ) != MPI_SUCCESS) {
			PL_ERROSH << "[ERROR]MPIPolylib::gather_polygons()"
					  << ":MPI_Recv,MPITAG_TRIA_IDS faild.:rank=" << rank 
					  << endl;
			return PLSTAT_MPI_ERROR;
		}

		// 三角形リストを受信
		float *p_triaarray = new float[ total_tria_num * 3 * 3 ];
		if (MPI_Recv( p_triaarray, total_tria_num * 3 * 3, MPI_FLOAT, rank, 
					MPITAG_TRIAS, m_mycomm, &mpi_stat ) != MPI_SUCCESS) {
			PL_ERROSH << "[ERROR]MPIPolylib::gather_polygons()"
					  << ":MPI_Recv,MPITAG_TRIA_IDS faild.:rank=" << rank 
					  << endl;
			return PLSTAT_MPI_ERROR;
		}

		// 各ポリゴングループに対して受信した三角形情報を追加
		pos_id = 0;
		pos_tria = 0;
		for( i=0; i<m_pg_list.size() * 2; i++ ){

			// ポリゴングループID
			int pg_id = p_intarray[i];

			// 当該ポリゴングループの三角形数
			unsigned int num_trias = p_intarray[i+1];

			// グループIDのポリゴングループインスタンス取得
			PolygonGroup* p_pg = get_group( pg_id );
			if( p_pg == NULL ) {
				PL_ERROSH << "[ERROR]MPIPolylib::gather_polygons():invalid pg_id:"
						  << pg_id << endl;
				return PLSTAT_NG;
			}

			// PrivateTriangleのベクタ生成
			vector<PrivateTriangle*> tria_vec;

			// ベクタに受信データ内容を設定
			for( j=0; j<num_trias; j++ ) {
				tria_vec.push_back( new PrivateTriangle(&p_triaarray[pos_tria], p_idarray[pos_id]) );
				pos_id++;
				pos_tria+=9;
			}

			// ポリゴングループに三角形リストを追加
			if( (ret = p_pg->add_triangles( &tria_vec )) != PLSTAT_OK ) {
				PL_ERROSH << "[ERROR]MPIPolylib::gather_polygons():p_pg->init() failed. returns:" << PolylibStat2::String(ret) << endl;
				return ret;
			}

			// ベクタの内容あとしまつ
			for( j=0; j<num_trias; j++ ) {
				delete tria_vec.at(j);
			}

			// 次のポリゴングループID位置へ
			i++;
		}

		// 受信領域あとしまつ
		if( p_intarray != NULL )  delete[] p_intarray;
		if( p_idarray != NULL )   delete[] p_idarray;
		if( p_triaarray != NULL ) delete[] p_triaarray;
	}

	return PLSTAT_OK;
}


// protected //////////////////////////////////////////////////////////////////
POLYLIB_STAT
MPIPolylib::send_polygons_to_rank0(
)
{
#ifdef DEBUG
	PL_DBGOSH << "MPIPolylib::send_polygons_to_rank0() in. " << endl;
#endif
	unsigned int i;
	vector<PolygonGroup*>::iterator group_itr;
	PolygonGroup *p_pg;
	vector<PrivateTriangle*> const *p_trias;

	vector<int>   send_num_trias;
	int          *p_send_num_trias_array;
	vector<int>   send_tria_ids;
	int          *p_send_tria_ids_array;
	vector<float> send_trias;
	float        *p_send_trias_array;

	// 送信用一時データ初期化
	p_trias = NULL;
	send_num_trias.clear();
	send_trias.clear();
	send_tria_ids.clear();
	
	// 全グループに対して
	for (group_itr = m_pg_list.begin(); group_itr != m_pg_list.end(); group_itr++) {
		p_pg = (*group_itr);

		// ポリゴン情報を持つグループだけ
		if( p_pg->get_triangles() != NULL && p_pg->get_triangles()->size() != 0 ) {

			// 当該PE領域内に一部でも含まれるポリゴンを検索
			p_trias = p_pg->search( &(m_myproc.m_area.m_gcell_bbox), false );
		}

		// グループIDと当該グループの三角形数の対を送信データに追加
		pack_num_trias( &send_num_trias, p_pg->get_internal_id(), p_trias );

		// 三角形情報を送信データに追加
		pack_trias( &send_trias, p_trias );

		// 三角形ID情報を送信データに追加
		pack_tria_ids( &send_tria_ids, p_trias );

		// search結果の後始末
		if( p_trias ) delete p_trias;
		p_trias = NULL;
	}

	// 送信データをシリアライズ
	p_send_num_trias_array = new int[ send_num_trias.size() ];
	for( i=0; i<send_num_trias.size(); i++ )
		p_send_num_trias_array[i] = send_num_trias[i];

	p_send_tria_ids_array  = new int[ send_tria_ids.size() ];
	for( i=0; i<send_tria_ids.size(); i++ )
		p_send_tria_ids_array[i] = send_tria_ids[i];

	p_send_trias_array =     new float[ send_trias.size() ];
	for( i=0; i<send_trias.size(); i++ )
		p_send_trias_array[i] = send_trias[i];

	// rank0へ送信
#ifdef DEBUG
	PL_DBGOSH << "sending polygons rank:" << m_myrank << " -> rank:0 ";
	for( i=0; i< send_num_trias.size(); i+=2 ) {
		PL_DBGOS << "(gid:" << send_num_trias[i] 
				 << ",num_tria:" << send_num_trias[i+1] << ")";
	}
	PL_DBGOS << endl;
#endif
	if (MPI_Send( p_send_num_trias_array, send_num_trias.size(), MPI_INT, 0,
					MPITAG_NUM_TRIAS, m_mycomm ) != MPI_SUCCESS) {
		PL_ERROSH << "[ERROR]MPIPolylib::send_polygons_to_rank0()"
				  << ":MPI_Send,MPITAG_NUM_TRIAS faild." << endl;
		return PLSTAT_MPI_ERROR;
	}
	if (MPI_Send( p_send_tria_ids_array,  send_tria_ids.size(),  MPI_INT, 0,
					MPITAG_TRIA_IDS,  m_mycomm ) != MPI_SUCCESS) {
		PL_ERROSH << "[ERROR]MPIPolylib::send_polygons_to_rank0()"
				  << ":MPI_Send,MPITAG_TRIA_IDS faild." << endl;
		return PLSTAT_MPI_ERROR;
	}
	if (MPI_Send( p_send_trias_array,     send_trias.size(),   MPI_FLOAT, 0,
					MPITAG_TRIAS,     m_mycomm ) != MPI_SUCCESS) {
		PL_ERROSH << "[ERROR]MPIPolylib::send_polygons_to_rank0()"
				  << ":MPI_Send,MPITAG_TRIAS faild." << endl;
		return PLSTAT_MPI_ERROR;
	}

	// 送信データの後始末
	delete[] p_send_num_trias_array;
	delete[] p_send_tria_ids_array;
	delete[] p_send_trias_array;

	return PLSTAT_OK;
}

// protected //////////////////////////////////////////////////////////////////
POLYLIB_STAT
MPIPolylib::select_excluded_trias(
	PolygonGroup *p_pg
)
{
#ifdef DEBUG
	PL_DBGOSH << "MPIPolylib::select_excluded_trias() in. " << endl;
#endif

	unsigned int i, j;
	vector<PrivateTriangle*> const *p_trias;
	vector<int> ids;

	// 全隣接PEについて
	for( i=0; i<m_neibour_procs.size(); i++ ) {
		ids.clear();

		// 隣接PE領域(ガイドセル含)に懸かる三角形IDリストを作成
		p_trias = p_pg->search( &(m_neibour_procs.at(i)->m_area.m_gcell_bbox), false );
		for( j=0; j<p_trias->size(); j++ ) {
			ids.push_back( p_trias->at(j)->get_id() );
		}
#ifdef DEBUG
	PL_DBGOSH << "gid:" << p_pg->get_id() << " neibour_rank:" << m_neibour_procs.at(i)->m_rank
			  << " 除外三角形数:" << ids.size() << endl;
#endif

		// migrate除外三角形IDマップに追加
		m_neibour_procs.at(i)->m_exclusion_map[p_pg->get_internal_id()] = ids;

		// search結果あとしまつ
		delete p_trias;
	}
	return PLSTAT_OK;

}

// eof
