/*
 * Polylib - Polygon Management Library.
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2010-
 *
 */
#ifndef polylib_triangle_h
#define polylib_triangle_h

#include "common/Vec3.h"

namespace PolylibNS{

////////////////////////////////////////////////////////////////////////////
///
/// クラス:Triangle
/// 入出力用インターフェースクラスであり、本ヘッダに対応する.cxxファイルは存在
/// しない。
///
////////////////////////////////////////////////////////////////////////////
class Triangle {
public:
	///
	/// コンストラクタ。
	///
	Triangle(){};

	///
	/// コンストラクタ。
	///
	/// @param[in] vertex ポリゴンの頂点。
	/// @attention 面積と法線はvertexを元に自動計算される。
	///
	Triangle(
		Vec3f	vertex[3]
	) {
		m_vertex[0] = vertex[0];
		m_vertex[1] = vertex[1];
		m_vertex[2] = vertex[2];
		calc_normal();
		calc_area();
	}

	///
	/// コンストラクタ。
	///
	/// @param[in] vertex	ポリゴンの頂点。
	/// @param[in] normal	法線。
	/// @attention 面積はvertexを元に自動計算される。
	///
	Triangle(
		Vec3f	vertex[3], 
		Vec3f	normal
	) {
		m_vertex[0] = vertex[0];
		m_vertex[1] = vertex[1];
		m_vertex[2] = vertex[2];
		m_normal = normal;
		calc_area();
	}

	///
	/// コンストラクタ。
	///
	/// @param[in] vertex	ポリゴンの頂点。
	/// @param[in] normal	法線。
	/// @param[in] area		ポリゴンの面積。
	///
	Triangle(
		Vec3f	vertex[3], 
		Vec3f	normal, 
		float	area
	) {
		m_vertex[0] = vertex[0];
		m_vertex[1] = vertex[1];
		m_vertex[2] = vertex[2];
		m_normal = normal;
		m_area = area;
	}

	//=======================================================================
	// Setter/Getter
	//=======================================================================
	///
	/// 頂点を設定。
	///
	/// @param[in] vertex		三角形の3頂点。
	/// @param[in] calc_normal	法線ベクトルを再計算するか？
	/// @param[in] calc_area	面積を再計算するか？
	///
	void set_vertexes(
		Vec3f	vertex[3], 
		bool	calc_normal, 
		bool	calc_area
	) {
		m_vertex[0] = vertex[0];
		m_vertex[1] = vertex[1];
		m_vertex[2] = vertex[2];
		if(calc_normal) this->calc_normal();
		if(calc_area) this->calc_area();
	}

	///
	/// vertexの配列を取得。
	///
	/// @return vertexの配列。
	///
	Vec3f* get_vertex() const {
		return const_cast<Vec3f*>(m_vertex);
	}

	///
	/// 法線ベクトルを取得。
	///
	/// @return 法線ベクトル。
	///
	Vec3f get_normal() const {
		return m_normal;
	}

	///
	/// 面積を取得。
	///
	/// @return 面積。
	///
	float get_area() const {
		return m_area;
	}

	///
	/// ユーザ定義IDを設定。
	///
	///
	void set_exid( int id ) {
		m_exid = id;
	}

	///
	/// ユーザ定義IDを取得。
	///
	/// @return ユーザ定義ID。
	///
	int get_exid() const {
		return m_exid;
	}

protected:
	///
	/// 法線ベクトル算出。
	///
	void calc_normal() {
		Vec3f a = m_vertex[1] - m_vertex[0];
		Vec3f b = m_vertex[2] - m_vertex[0];
		m_normal = (cross(a,b)).normalize();

	}

	///
	/// 面積算出。
	///
	void calc_area() {
		Vec3f a = m_vertex[1] - m_vertex[0];
		Vec3f b = m_vertex[2] - m_vertex[0];
		float al = a.length();
		float bl = b.length();
		float ab = dot(a,b);
		m_area = 0.5*sqrtf(al*al*bl*bl - ab*ab);
	}

	//=======================================================================
	// クラス変数
	//=======================================================================
	/// 三角形の頂点座標（反時計回りで並んでいる）。
	Vec3f	m_vertex[3];

	/// 三角形の法線ベクトル。
	Vec3f	m_normal;

	/// 三角形の面積。
	float	m_area;

	/// 三角形のユーザ定義ID
	int     m_exid;
};

////////////////////////////////////////////////////////////////////////////
///
/// クラス:PrivateTriangleクラス
/// Polylib内のデータ保存用の基本クラスです。
///
////////////////////////////////////////////////////////////////////////////
class PrivateTriangle : public Triangle {
public:
	///
	/// コンストラクタ。
	///
	/// @param[in] vertex	ポリゴンの頂点。
	/// @param[in] id		三角形ポリゴンID。
	///
	PrivateTriangle(
		Vec3f	vertex[3], 
		int		id
	) : Triangle(vertex) {
		m_id = id;
	}

	///
	/// コンストラクタ。
	///
	/// @param[in] vertex	ポリゴンの頂点。
	/// @param[in] normal	法線。
	/// @param[in] id		三角形ポリゴンID。
	///
	PrivateTriangle(
		Vec3f	vertex[3], 
		Vec3f	normal, 
		int		id
	) : Triangle(vertex, normal) {
		m_id = id;
	}

	///
	/// コンストラクタ。
	///
	/// @param[in] vertex	ポリゴンの頂点。
	/// @param[in] normal	法線。
	/// @param[in] area		ポリゴンの面積。
	/// @param[in] id		三角形ポリゴンID。
	///
	PrivateTriangle(
		Vec3f	vertex[3], 
		Vec3f	normal, 
		float	area, 
		int		id
	) : Triangle(vertex, normal, area) {
		m_id = id;
	}

	///
	/// コンストラクタ。
	///
	/// @param[in] tri		ポリゴン。
	/// @param[in] id		三角形ポリゴンID。
	///
	PrivateTriangle(
		Triangle	tri, 
		int			id
	) : Triangle(tri.get_vertex(), tri.get_normal()) {
		m_id = id;
	}

	///
	/// コンストラクタ。
	///
	/// @param[in] tri		ポリゴン。
	///
	PrivateTriangle(
		const PrivateTriangle	&tri 
	) : Triangle(tri.get_vertex(), tri.get_normal()) {
		m_id = tri.m_id;
	}

	///
	/// コンストラクタ。
	///
	/// @param[in] dim		ポリゴン頂点座標配列。
	/// @param[in] id		三角形ポリゴンID。
	///
	PrivateTriangle(
		float		*dim,
		int			id
	){
		for( int i=0; i<3; i++ ) {
			m_vertex[i].t[0] = *dim++;
			m_vertex[i].t[1] = *dim++;
			m_vertex[i].t[2] = *dim++;
		}
		m_id = id;
		calc_normal();
		calc_area();
	}

	//=======================================================================
	// Setter/Getter
	//=======================================================================
	///
	/// 三角形ポリゴンIDを設定。
	///
	///  @param[in] id	三角形ポリゴンID。
	///
	void set_id(int id)				{m_id = id;}

	///
	/// 三角形ポリゴンIDを返す。
	///
	///  @return 三角形ポリゴンID。
	///
	int get_id() const				{return m_id;}

protected:
	//=======================================================================
	// クラス変数
	//=======================================================================
	/// PolygonGroup内で一意となる三角形ポリゴンID。
	///
	int m_id;
};

} //namespace PolylibNS

#endif  // polylib_triangle_h

