/*
 * Polylib - Polygon Management Library.
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2010-
 *
 */
#ifndef bbox_h
#define bbox_h

#include <algorithm>
#include <list>
#include "common/Vec2.h"
#include "common/Vec3.h"
#include "common/axis.h"

namespace PolylibNS {

////////////////////////////////////////////////////////////////////////////
///
/// クラス:BBox
/// Bounding Boxを管理するクラス
///
////////////////////////////////////////////////////////////////////////////
class BBox {
public:
	BBox() { init(); }

	BBox(float _minx, float _miny, float _minz, float _maxx, float _maxy, float _maxz) {
		min[0] = _minx; min[1] = _miny; min[2] = _minz;
		max[0] = _maxx; max[1] = _maxy; max[2] = _maxz;
	}

	BBox(float _min[3], float _max[3]) {
		min[0] = _min[0]; min[1] = _min[1]; min[2] = _min[2];
		max[0] = _max[0]; max[1] = _max[1]; max[2] = _max[2];
	}

	BBox(const Vec3f& _min, const Vec3f& _max) {
		min = _min;
		max = _max;
	}

	void init() {
		float s = 1e10;
		min[0] = s;  min[1] = s;  min[2] = s;
		max[0] = -s; max[1] = -s; max[2] = -s;
	}

	void setMinMax(const Vec3f& _min, const Vec3f& _max) {
		min = _min;
		max = _max;
	}

	void add(const Vec3f& v) {
		min[0] = std::min(min[0], v[0]);
		min[1] = std::min(min[1], v[1]);
		min[2] = std::min(min[2], v[2]);
		max[0] = std::max(max[0], v[0]);
		max[1] = std::max(max[1], v[1]);
		max[2] = std::max(max[2], v[2]);
	}

	Vec3f getPoint(int idx) const {
		int x, y, z;
		Vec3f p = min;
		x = idx    & 0x1;
		y = idx>>1 & 0x1;
		z = idx>>2 & 0x1;
		if (x) p[0] = max[0];
		if (y) p[1] = max[1];
		if (z) p[2] = max[2];
		return p;
	}

	Vec3f center() const { return .5 * (min + max); }
	Vec3f size() const { return max - min; }
	float xsize() const { return max[0] - min[0]; }
	float ysize() const { return max[1] - min[1]; }
	float zsize() const { return max[2] - min[2]; }
	float length(const AxisEnum& axis) const { return max[axis] - min[axis]; }
	float diameter() const { return (max - min).length(); }
	AxisEnum getMaxAxis(float& length) const {
		Vec3f len = size();
		AxisEnum axis;
		if (len[0] > len[1]) {
			if (len[0] > len[2])	axis = AXIS_X;
			else					axis = AXIS_Z;
		}
		else {
			if (len[1] > len[2])	axis = AXIS_Y;
			else					axis = AXIS_Z;
		}
		length = len[axis];
		return axis;
	}

	///
	/// 引数で与えられた点が、このBBoxに含まれるかを判定する。
    /// @param[in] pos 試行する点
    /// @return 含まれる場合はtrue。他はfalse。
	///
	bool contain(const Vec3f& pos) const {
		if (min.t[0] <= pos.t[0] && pos.t[0] <= max.t[0] &&
			min.t[1] <= pos.t[1] && pos.t[1] <= max.t[1] &&
			min.t[2] <= pos.t[2] && pos.t[2] <= max.t[2]) {
			return true;
		}
		else {
			return false;
		}
	}

	///
	/// BBoxとBBoxの交差判定を行う。
	/// KD-Treeの交差判定と同じ。
    /// @param[in] bbox 試行するBBox
    /// @return 交差する場合はtrue。他はfalse。
	///
	bool crossed(const BBox& bbox) const {
		if (max.t[0] < bbox.min.t[0] || bbox.max.t[0] < min.t[0]) return false;
		if (max.t[1] < bbox.min.t[1] || bbox.max.t[1] < min.t[1]) return false;
		if (max.t[2] < bbox.min.t[2] || bbox.max.t[2] < min.t[2]) return false;
		return true;
	}

	///
	/// BBoxとBBoxの重複領域の抽出を行う。
	/// 自身の面と他方の辺との交差判定を行う。
    /// @param[in] other_bbox 試行するBBox
    /// @return 交差する場合はtrue。他はfalse。
	///
	BBox getCrossedRegion(BBox& other_bbox) const {

		// axis_idに垂直な軸配列
		// int axis[][2] = { {1,2}, {2,0}, {0,1} };

		// 点の包含判定
		std::list<Vec3f> inc_own_points;
		std::list<Vec3f> inc_other_points;
		for (int i=0; i<8; i++) {
			// 自分に相手の8頂点が含まれるかチェック
			if (contain(other_bbox.getPoint(i)) == true) {
				inc_other_points.push_back(other_bbox.getPoint(i));
			}
			// 相手に自分の8頂点が含まれるかチェック
			if (other_bbox.contain(getPoint(i)) == true) {
				inc_own_points.push_back(getPoint(i));
			}
		}

		// 自分のBBox面 と 他のBBox辺 の交差判定
		std::list<Vec3f> cross_own_face_with_other_side;
		for (int axis_id = 0; axis_id < 3; axis_id++) {
			Vec3f face[2][2];	//face[2面][長方形の対角点]
			Vec3f side[4][2];	//side[4辺][辺の両端の2点] :対点はmin-maxの順
			getFace(axis_id, face);
			other_bbox.getSide(axis_id, side);

			for (int face_id = 0; face_id < 2; face_id++) {
				for (int side_id = 0; side_id < 4; side_id++) {
					// 面と点が今見ている軸上でどこにあるかを検出。
					float z0 = (side[side_id][0])[axis_id];//辺id,点id,xyz
					float z1 = (side[side_id][1])[axis_id];//辺id,点id,xyz
					float z  = (face[face_id][0])[axis_id];//面id,点id,xyz

					if (z0 <= z && z <= z1) {
						// 2次元平面と点との包含判定
						Vec2f min = vec3to2(axis_id,face[face_id][0]);
						Vec2f max = vec3to2(axis_id,face[face_id][1]);
						Vec2f pos = vec3to2(axis_id,side[side_id][0]);

						std::cout << "min:" << min << std::endl;
						std::cout << "max:" << max << std::endl;
						std::cout << "pos:" << pos << std::endl;

						if (min[0] <= pos[0] && pos[0] <= max[0] &&
							min[1] <= pos[1] && pos[1] <= max[1]) {
							Vec3f vtx = side[side_id][0];
							vtx[axis_id] = z;
							cross_own_face_with_other_side.push_back(vtx);
						}
					}
				}
			}
		}

		// 他のBBox面 と 自分のBBox辺 の交差判定
		std::list<Vec3f> cross_other_face_with_own_side;
		for (int axis_id = 0; axis_id < 3; axis_id++) {
			Vec3f face[2][2];	//face[2面][長方形の対角点]
			Vec3f side[4][2];	//side[4辺][辺の両端の2点] :対点はmin-maxの順
			other_bbox.getFace(axis_id, face);
			getSide(axis_id, side);

			for (int face_id = 0; face_id < 2; face_id++) {
				for (int side_id = 0; side_id < 4; side_id++) {
					// 面と点が今見ている軸上でどこにあるかを検出。
					float z0 = (side[side_id][0])[axis_id];//辺id,点id,xyz
					float z1 = (side[side_id][1])[axis_id];//辺id,点id,xyz
					float z  = (face[face_id][0])[axis_id];//面id,点id,xyz

					if (z0 <= z && z <= z1) {
						// 2次元平面と点との包含判定
						Vec2f min = vec3to2(axis_id,face[face_id][0]);
						Vec2f max = vec3to2(axis_id,face[face_id][1]);
						Vec2f pos = vec3to2(axis_id,side[side_id][0]);

						if (min[0] <= pos[0] && pos[0] <= max[0] &&
							min[1] <= pos[1] && pos[1] <= max[1]) {
							Vec3f vtx = side[side_id][0];
							vtx[axis_id] = z;
							cross_other_face_with_own_side.push_back(vtx);
						}
					}
				}
			}
		}

		//返り値
		BBox cross_bbox;
		cross_bbox.init();
		std::list<Vec3f>::iterator itr;

		if (inc_own_points.size() == 1 && inc_other_points.size() == 1) {
			itr = inc_own_points.begin();
			cross_bbox.add(*itr);
			itr = inc_other_points.begin();
			cross_bbox.add(*itr);
		}
		else if ((inc_own_points.size() == 2 && inc_other_points.size() == 0) ||
				 (inc_own_points.size() == 4 && inc_other_points.size() == 0)) {
			itr = inc_own_points.begin();
			for (; itr != inc_own_points.end(); itr++) {
				cross_bbox.add(*itr);
			}
			itr = cross_other_face_with_own_side.begin();
			for (;itr != cross_other_face_with_own_side.end(); itr++) {
				cross_bbox.add(*itr);
			}
		}
		else if ((inc_own_points.size() == 0 && inc_other_points.size() == 2) ||
				 (inc_own_points.size() == 0 && inc_other_points.size() == 4)) {
			itr = inc_other_points.begin();
			for (; itr != inc_other_points.end(); itr++) {
				cross_bbox.add(*itr);
			}
			itr = cross_own_face_with_other_side.begin();
			for (; itr != cross_own_face_with_other_side.end(); itr++) {
				cross_bbox.add(*itr);
			}
		}
		else if (inc_own_points.size() == 8 && inc_other_points.size() == 0) {
			itr = inc_own_points.begin();
			for (; itr != inc_own_points.end(); itr++) {
				cross_bbox.add(*itr);
			}
		}
		else if (inc_own_points.size() == 0 && inc_other_points.size() == 8) {
			itr = inc_other_points.begin();
			for (; itr != inc_other_points.end(); itr++) {
				cross_bbox.add(*itr);
			}
		}
		else { // 貫通ケース
			itr = cross_other_face_with_own_side.begin();
			for (;itr != cross_other_face_with_own_side.end(); itr++) {
				cross_bbox.add(*itr);
			}
			itr = cross_own_face_with_other_side.begin();
			for (;itr != cross_own_face_with_other_side.end(); itr++) {
				cross_bbox.add(*itr);
			}
		}

		return cross_bbox;
	}

	///
	/// 引数axis_id(0=x,1=y,z=2)に垂直な成分を詰めて返す。
	///
	Vec2f vec3to2(int axis_id, Vec3f& v3) const {
		// axis perpendicular to the axis_id
		int axis[][2] = { {1,2}, {2,0}, {0,1} };
		Vec2f v2;
		v2[0] = v3[axis[axis_id][0]];
		v2[1] = v3[axis[axis_id][1]];
		return v2;
	}

	///
	/// 引数axis_id(0=x,1=y,z=2)に垂直な、このBBoxの面の対角点を返す。
    /// @param[in] axis_id 軸番号。0=x軸、1=y軸、2=z軸。
    /// @param[in] face BBoxの面の中で、軸に垂直な面の対角点。
	///
	// get 2 faces perpendicular to the axis_id
	void getFace(int axis_id, Vec3f face[2][2]) const {
		// min-max vertex indices of face
		int p_minmax[][2][2] = {
			{{0,6},{1,7}},
			{{0,5},{2,7}},
			{{0,3},{4,7}}
		};
		for (int i = 0; i < 2; i++) {
			for (int j = 0; j < 2; j++) {
				face[i][j] = getPoint(p_minmax[axis_id][i][j]);
#ifdef DEBUG
				std::cout.setf(std::ios::scientific, std::ios::floatfield);
				std::cout << "get_face[" 
					<< axis_id << "][" << i << "][" << j << "] "
					<< p_minmax[axis_id][i][j] << ": " 
					<< face[i][j][0] << " " << face[i][j][1] << " " 
					<< face[i][j][2] << std::endl;
				std::cout.unsetf(std::ios::scientific);
#endif
			}
		}
	}

	///
	/// 引数axis_id(0=x,1=y,z=2)に平行な、このBBoxの辺の端点を返す。
    /// @param[in] axis_id 軸番号。0=x軸、1=y軸、2=z軸。
    /// @param[in] side BBoxの辺の中で、軸に平行な辺の端点。
	///
	void getSide(int axis_id, Vec3f side[4][2]) const {
		int p_minmax[][4][2] = {
			{{0,1},{2,3}, {4,5},{6,7}},
			{{0,2},{1,3}, {4,6},{5,7}},
			{{0,4},{1,5}, {2,6},{3,7}}
		};
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 2; j++) {
				side[i][j] = getPoint(p_minmax[axis_id][i][j]);
#ifdef DEBUG
				std::cout.setf(std::ios::scientific, std::ios::floatfield);
				std::cout << "get_side[" 
					<< axis_id << "][" << i << "][" << j << "] "
					<< p_minmax[axis_id][i][j] << ": " 
					<< side[i][j][0] << " " << side[i][j][1] << " " 
					<< side[i][j][2] << std::endl;
				std::cout.unsetf(std::ios::scientific);
#endif
			}
		}
	}

	Vec3f min, max;
};

} // namespace PolylibNS

#endif  // bbox_h

