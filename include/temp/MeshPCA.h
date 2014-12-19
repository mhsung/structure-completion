/**
* file:	MeshPCA.h
* description:	Principal Component Analysis
*
* author:	Minhyuk Sung
* date:	(First) August 2009
		(Last) January 2014
*/

#ifndef _MESH_PCA_H_
#define _MESH_PCA_H_

#define FIXED_AXIS_BBOX

#include "pca.h"
#include "MyMesh.h"
#include <vector>

class MeshPCA
{
public:
	MeshPCA();
	static MeshPCA *constructor(std::vector<MyMesh::Point> &_points);
	static MeshPCA *constructor(std::vector<MyMesh::Point> &_points, std::vector<double> &_weights);
	static MeshPCA *constructor(MyMesh *_mesh);
	static MeshPCA *constructor(MyMesh *_mesh, std::vector<faceIndex> &_subset_fids);

private:
	MeshPCA(std::vector<MyMesh::Point> &_points);
	MeshPCA(std::vector<MyMesh::Point> &_points, std::vector<double> &_weights);
	MeshPCA(MyMesh *_mesh);
	MeshPCA(MyMesh *_mesh, std::vector<faceIndex> &_subset_fids);

public:
	~MeshPCA();

public:
	const MyMesh::Point get_pca_origin()const { return origin_; };
	const MyMesh::Normal get_pca_bounding_box_size()const { return bbox_size_; };
	const MyMesh::Normal *get_pca_axes()const { return axes_; };
	const MyMesh::Point *get_pca_bounding_box()const { return bbox_; };
	const float get_pca_diameter()const { return diameter_; };

	void set_pca_origin(MyMesh::Point _p) { origin_ = _p; };
	void set_pca_bounding_box_size(MyMesh::Normal _s) { bbox_size_ = _s; };

	// _axis: 0 - x, 1 - y, 2 - z axis
	void change_pca_axis_signs(short unsigned int _axis_order);

private:
	std::vector<MyMesh::Point> points_;
	std::vector<dReal> weights_;
	unsigned int num_data_;

	MyMesh::Point origin_;
	MyMesh::Normal bbox_size_;
	MyMesh::Normal axes_[3];
	std::vector<MyMesh::Point> normalized_points_;
	std::vector<MyMesh::Point> pca_points_;
	MyMesh::Point bbox_[8];
	double diameter_;

	void initialize();
	bool compute();
	void compute_origin();
	bool compute_axis();
	void compute_bounding_box();
	MyMesh::Point from_orig_to_pca_coord(MyMesh::Point _p);
	MyMesh::Point from_pca_to_orig_coord(MyMesh::Point _p);
};

#endif	// _MESH_PCA_H_