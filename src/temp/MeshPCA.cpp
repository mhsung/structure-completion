#include "MeshPCA.h"
#include "MeshPointSampling.h"

MeshPCA::MeshPCA()
{
	initialize();
}

MeshPCA *MeshPCA::constructor(std::vector<MyMesh::Point> &_points)
{
	if (_points.empty()) return NULL;

	MeshPCA *obj = new MeshPCA(_points);
	if (obj->num_data_ == 0)
	{
		delete obj;
		return NULL;
	}

	return obj;
}

MeshPCA *MeshPCA::constructor(std::vector<MyMesh::Point> &_points, std::vector<double> &_weights)
{
	if (_points.empty()
		|| _points.size() != _weights.size()) return NULL;

	MeshPCA *obj = new MeshPCA(_points, _weights);
	if (obj->num_data_ == 0)
	{
		delete obj;
		return NULL;
	}

	return obj;
}

MeshPCA *MeshPCA::constructor(MyMesh *_mesh)
{
	if (_mesh->n_faces() == 0) return NULL;

	MeshPCA *obj = new MeshPCA(_mesh);
	if (obj->num_data_ == 0)
	{
		delete obj;
		return NULL;
	}

	return obj;
}

MeshPCA *MeshPCA::constructor(MyMesh *_mesh, std::vector<faceIndex> &_subset_fids)
{
	if (_mesh->n_faces() == 0
		|| _subset_fids.empty()) return NULL;

	MeshPCA *obj = new MeshPCA(_mesh, _subset_fids);
	if (obj->num_data_ == 0)
	{
		delete obj;
		return NULL;
	}

	return obj;
}

MeshPCA::MeshPCA(std::vector<MyMesh::Point> &_points)
{
	initialize();

	points_ = _points;
	num_data_ = points_.size();
	weights_.resize(num_data_, 1.0);

	compute();
}

MeshPCA::MeshPCA(std::vector<MyMesh::Point> &_points,
	std::vector<double> &_weights)
{
	initialize();

	points_ = _points;
	weights_ = _weights;
	num_data_ = points_.size();
	assert(weights_.size() == num_data_);

	compute();
}

MeshPCA::MeshPCA(MyMesh *_mesh)
{
	initialize();

	//num_data_ = _mesh->n_vertices();
	//points_.resize(num_data_);
	//weights_.resize(num_data_);

	//// Vertex area weights
	//_mesh->request_vertex_areas();

	//for (unsigned int k = 0; k < num_data_; k++)
	//{
	//	MyMesh::VertexHandle vh = _mesh->vertex_handle(k);
	//	points_[k] = _mesh->point(vh);
	//	//weights_[k] = _mesh->property(_mesh->v_area_, vh);
	//	weights_[k] = 1;
	//}

	MeshPointSampling *sampling = new MeshPointSampling(_mesh);
	num_data_ = sampling->get_total_num_samples();

	points_.clear();
	points_.reserve(num_data_);
	weights_.resize(num_data_, 1.0);

	const std::vector<std::vector<MyMesh::Point>> &face_samples = sampling->get_face_samples();
	for (std::vector<std::vector<MyMesh::Point>>::const_iterator it = face_samples.begin();
		it != face_samples.end(); ++it)
	{
		points_.insert(points_.end(), (*it).begin(), (*it).end());
	}

	delete sampling;
	compute();
}

MeshPCA::MeshPCA(MyMesh *_mesh, std::vector<faceIndex> &_subset_fids)
{
	initialize();

	MeshPointSampling *sampling = new MeshPointSampling(_mesh, _subset_fids);
	num_data_ = sampling->get_total_num_samples();

	points_.clear();
	points_.reserve(num_data_);
	weights_.resize(num_data_, 1.0);

	const std::vector<std::vector<MyMesh::Point>> &face_samples = sampling->get_face_samples();
	for (std::vector<std::vector<MyMesh::Point>>::const_iterator it = face_samples.begin();
		it != face_samples.end(); ++it)
	{
		points_.insert(points_.end(), (*it).begin(), (*it).end());
	}

	delete sampling;
	compute();
}

MeshPCA::~MeshPCA()
{

}

void MeshPCA::initialize()
{
	num_data_ = 0;
	origin_ = MyMesh::Point(0.0);

	for (unsigned int i = 0; i < 3; i++)
	{
		// Default: x, y, z axes
		axes_[i] = MyMesh::Point(0.0);
		axes_[i][i] = 1.0;
		bbox_size_[i] = 0.0;
	}
		

	for (unsigned int i = 0; i < 8; i++)
	{
		bbox_[i] = MyMesh::Point(0.0);
	}

	points_.clear();
	weights_.clear();
	normalized_points_.clear();
	pca_points_.clear();
}

bool MeshPCA::compute()
{
	if (num_data_ == 0)
		return false;

	compute_origin();

#ifndef FIXED_AXIS_BBOX
	if (!compute_axis())
		return false;
#endif

	compute_bounding_box();
	return true;
}

void MeshPCA::compute_origin()
{
	normalized_points_.resize(num_data_);
	origin_ = MyMesh::Point(0.0);

	double total_weights = 0.0;
	for (unsigned int k = 0; k < num_data_; k++)
	{
		// weighted average position
		normalized_points_[k] = points_[k] * weights_[k];
		origin_ += normalized_points_[k];
		total_weights += weights_[k];
	}
	origin_ /= total_weights;
}

bool MeshPCA::compute_axis()
{
	const unsigned int num_dimension = 3;

	assert(num_data_ > 0);
	assert(points_.size() == num_data_);
	assert(weights_.size() == num_data_);


	CPCA *pca = new CPCA();
	Mat_DP input(num_data_, num_dimension);
	Vec_DP mean(0.0, num_dimension);
	Mat_DP output(num_dimension, num_dimension);


	//for (unsigned int k = 0; k < num_data_; k++)
	//	normalized_points_[k] -= mean_normalized_points;

	for (unsigned int k = 0; k < num_data_; k++)
		for (unsigned int j = 0; j < num_dimension; j++)
			input[k][j] = normalized_points_[k][j];


	if (!pca->Setup(input, 1.0))
	{
		std::cout << "error: PCA cannot be computed." << std::endl;
		return false;
	}
	pca->GetMeanMat(mean);
	pca->GetFVMat(output);
	delete pca;

	// transpose
	for (unsigned int k = 0; k < num_dimension; k++)
		for (unsigned int j = 0; j < num_dimension; j++)
			axes_[k][j] = (float)output[j][k];


	// debug
	//MyMesh::Point debug_mean = MyMesh::Point(0.0);
	//for (unsigned int k = 0; k < num_data_; k++)
	//	debug_mean += normalized_points_[k];
	//debug_mean /= num_data_;
	//printf("Mean = (%lf, %lf, %lf)\n", debug_mean[0], debug_mean[1], debug_mean[2]);


	//printf("Mean = (%lf, %lf, %lf)\n", mean[0], mean[1], mean[2]);
	//for (unsigned int i = 0; i < num_dimension; i++)
	//{
	//	double dot_prod = 0.0;
	//	double length = 0.0;
	//	for (unsigned int j = 0; j < num_dimension; j++)
	//	{
	//		length += output[i][j] * output[i][j];
	//		dot_prod += output[i][j] * output[(i + 1) % 3][j];
	//	}
	//	length = sqrt(length);
	//	printf(">> [%d] Length = %lf\n", i, length);
	//	printf(">> [%d - %d] Dot product = %lf\n", i, (i + 1) % 3, dot_prod);
	//	printf("\n");
	//}

	//for (unsigned int k = 0; k < num_data_; k++)
	//{
	//	MyMesh::Point pca_p = from_orig_to_pca_coord(points_[k]);
	//	MyMesh::Point orig_p = from_pca_to_orig_coord(pca_p);
	//	double distance = (orig_p - points_[k]).length();
	//	printf(">> %lf\n", distance);
	//	if (abs(distance) > 1.0E-3)
	//		system("pause");
	//}

	return true;
}

void MeshPCA::compute_bounding_box()
{
	const unsigned int num_dimension = 3;


	// bounding box
	pca_points_.clear();
	pca_points_.resize(num_data_);
	MyMesh::Point bb_min = MyMesh::Point(max<float>());
	MyMesh::Point bb_max = MyMesh::Point(-max<float>());

	for (unsigned int k = 0; k < num_data_; k++)
	{
		MyMesh::Point &pca_p = pca_points_[k];
		pca_p = from_orig_to_pca_coord(points_[k]);

		for (unsigned int i = 0; i < num_dimension; i++)
		{
			bb_min[i] = std::min(bb_min[i], pca_p[i]);
			bb_max[i] = std::max(bb_max[i], pca_p[i]);
		}
	}

	for (unsigned int i = 0; i < num_dimension; i++)
		bbox_size_[i] = bb_max[i] - bb_min[i];

	diameter_ = (bb_max - bb_min).length();

	bbox_[0] = MyMesh::Point(bb_min[0], bb_min[1], bb_min[2]);
	bbox_[1] = MyMesh::Point(bb_min[0], bb_min[1], bb_max[2]);
	bbox_[2] = MyMesh::Point(bb_min[0], bb_max[1], bb_max[2]);
	bbox_[3] = MyMesh::Point(bb_min[0], bb_max[1], bb_min[2]);

	bbox_[4] = MyMesh::Point(bb_max[0], bb_min[1], bb_min[2]);
	bbox_[5] = MyMesh::Point(bb_max[0], bb_min[1], bb_max[2]);
	bbox_[6] = MyMesh::Point(bb_max[0], bb_max[1], bb_max[2]);
	bbox_[7] = MyMesh::Point(bb_max[0], bb_max[1], bb_min[2]);

	for (unsigned int i = 0; i < 8; i++)
		bbox_[i] = from_pca_to_orig_coord(bbox_[i]);

	// Reassign the origin value as the center of the bounding box
	//
	bb_max = from_pca_to_orig_coord(bb_max);
	bb_min = from_pca_to_orig_coord(bb_min);
	origin_ = 0.5 * (bb_min + bb_max);
	//
}

MyMesh::Point MeshPCA::from_orig_to_pca_coord(MyMesh::Point _p)
{
	const unsigned int num_dimension = 3;

	_p = _p - origin_;
	MyMesh::Point pca_p = MyMesh::Point(0.0);

	for (unsigned int i = 0; i < num_dimension; i++)
		for (unsigned int j = 0; j < num_dimension; j++)
			pca_p[i] += (axes_[i][j] * _p[j]);

	return pca_p;
}

MyMesh::Point MeshPCA::from_pca_to_orig_coord(MyMesh::Point _p)
{
	const unsigned int num_dimension = 3;

	MyMesh::Point orig_p = MyMesh::Point(0.0);

	for (unsigned int i = 0; i < num_dimension; i++)
		for (unsigned int j = 0; j < num_dimension; j++)
			orig_p[i] += (axes_[j][i] * _p[j]);

	orig_p = orig_p + origin_;
	return orig_p;
}

void MeshPCA::change_pca_axis_signs(short unsigned int _axis_order)
{
	switch (_axis_order)
	{
	case 0:
	case 1:
	case 2:
		axes_[_axis_order] = -axes_[_axis_order];
		break;
	}
}
