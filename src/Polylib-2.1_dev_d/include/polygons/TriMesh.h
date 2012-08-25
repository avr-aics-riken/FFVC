/*
 * Polylib - Polygon Management Library.
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2010-
 *
 */
#ifndef polylib_trimesh_h
#define polylib_trimesh_h

namespace PolylibNS {

class Polygons;
class VTree;
class PrivateTriangle;

////////////////////////////////////////////////////////////////////////////
///
/// クラス:TriMesh
/// 三角形ポリゴン集合を管理するクラス（KD木用に特化したクラス)。
///
////////////////////////////////////////////////////////////////////////////
class TriMesh : public Polygons {
public:
	///
	/// コンストラクタ。
	///
	TriMesh();

	///
	/// デストラクタ。
	///
	~TriMesh();

	///
	/// TriMeshクラスで管理する三角形ポリゴンリストを初期化し、引数で与えら
	/// れる三角形ポリゴンリストを設定する。
	/// 三角形ポリゴン用のメモリ領域は、Polylib内で新たに確保される。
	///
	///  @param[in] trias 設定する三角形ポリゴンリスト。
	///
	void init(
		const std::vector<PrivateTriangle*>	*trias
	);

	///
	/// 三角形ポリゴンリストに引数で与えられる三角形の複製を追加する。
	///
	/// @param[in] trias	設定する三角形ポリゴンリスト。
	/// @attention m_idが重複するインスタンスは追加されない。
	/// @attention KD木の再構築は行わない。
	///
	void add(
		const std::vector<PrivateTriangle*>  *trias
	);

	///
	/// ファイルからデータの初期化。
	///
	///  @param[in] fmap	ファイル名、ファイルフォーマット。
	///  @return PLSTAT_OK=成功/false=失敗
	///
	POLYLIB_STAT import(
		const std::map<std::string, std::string> fmap
	);

	///
	/// Polygonsクラスに含まれる全ポリゴン情報からKD木を作成する。
	///
	///  @return	POLYLIB_STATで定義される値が返る。
	///
	POLYLIB_STAT build();

	///
	/// TriMeshクラスが管理している三角形ポリゴン数を返す。
	///
	int triangles_num();

	///
	/// KD木探索により、指定矩形領域に含まれる三角形ポリゴンを抽出する。
	///
	///  @param[in] bbox	検索範囲を示す矩形領域。
	///  @param[in] every	true:3頂点が全て検索領域に含まれるものを抽出。
	///						false:1頂点でも検索領域に含まれるものを抽出。
	///  @return	抽出したポリゴンリストのポインタ。
	///	 @attention	三角形ポリゴンのメモリ領域は新たにPolylib内で確保される。
	///  @attention MPIPolylib内での利用が目的なので、ユーザは使用しないこと。
	///  @attention	オーバーロードメソッドあり。
	///
	const std::vector<PrivateTriangle*>* search(
		BBox	*bbox, 
		bool	every
	) const;

	///
	/// KD木探索により、指定矩形領域に含まれるポリゴンを抽出する。
	///
	///  @param[in]		bbox		検索範囲を示す矩形領域
	///  @param[in]		every		true:3頂点が全て検索領域に含まれるものを抽出。
	///								false:1頂点でも検索領域に含まれるものを抽出。
	///  @param[in,out] tri_list	抽出した三角形ポリゴンリストへのポインタ。
	///  @return	POLYLIB_STATで定義される値が返る。
	///  @attention	tri_listで戻される三角形ポリゴンのポインタは、Polylib内で
	///				保持されるアドレス値なので、ユーザはdeleteしないで下さい。
	///  @attention	オーバーロードメソッドあり。
	///
	POLYLIB_STAT search(
		BBox							*bbox,
		bool							every,
		std::vector<PrivateTriangle*>	*tri_list
	) const;

	///
	/// 線形探索により、指定矩形領域に含まれる三角形ポリゴンを抽出する。
	///
	///  @param[in] q_bbox	検索範囲を示す矩形領域。
	///  @param[in] every	true:3頂点が全て検索領域に含まれるものを抽出。
	///						false:1頂点でも検索領域に含まれるものを抽出。
	///  @return	抽出したポリゴンリストのポインタ。
	///	 @attention	三角形ポリゴンのメモリ領域は新たにPolylib内で確保される。
	///  @attention MPIPolylib内での利用が目的なので、ユーザは使用しないこと。
	///
	const std::vector<PrivateTriangle*>* linear_search(
		BBox	*q_bbox, 
		bool	every
	) const;

	///
	/// 線形探索により、指定矩形領域に含まれるポリゴンを抽出する。
	///
	///  @param[in]		q_bbox		検索範囲を示す矩形領域。
	///  @param[in]		every		true:3頂点が全て検索領域に含まれるものを抽出。
	///								false:1頂点でも検索領域に含まれるものを抽出。
	///  @param[in,out] tri_list	抽出した三角形ポリゴンリストへのポインタ。
	///  @return	POLYLIB_STATで定義される値が返る。
	///  @attention	tri_listで戻される三角形ポリゴンのポインタは、Polylib内で
	///				保持されるアドレス値なので、ユーザはdeleteしないで下さい。
	///  @attention	オーバーロードメソッドあり。
	///
	POLYLIB_STAT linear_search(
		BBox							*q_bbox, 
		bool							every,
		std::vector<PrivateTriangle*>	*tri_list
	) const;

	//=======================================================================
	// Setter/Getter
	//=======================================================================
	///
	/// TriMeshクラスが管理しているBoundingBoxを返す。
	///
	BBox get_bbox() const {
		return m_bbox;
	}

	///
	/// KD木クラスを取得。
	///
	/// @return KD木クラス。
	///
	VTree *get_vtree() const {
		return m_vtree;
	}

private:
	///
	/// 三角形ポリゴンリストの初期化。
	///
	void init_tri_list();

	//=======================================================================
	// クラス変数
	//=======================================================================
	/// 全三角形ポリゴンを外包するBoundingBox。
	BBox	m_bbox;

	/// KD木クラス。
	VTree	*m_vtree;

	/// MAX要素数。
	int		m_max_elements;
};

} //namespace PolylibNS

#endif  // polylib_trimesh_h
