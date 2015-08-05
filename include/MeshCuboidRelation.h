/**
* file:	MeshCuboidRelation.h
* description:	Define a mesh proxy part relations.
*
* author:	Minhyuk Sung
* date:	October 2014
*/

#ifndef _MESH_CUBOID_RELATION_H_
#define _MESH_CUBOID_RELATION_H_

#include "MyMesh.h"
#include "MeshCuboid.h"
#include "MeshCuboidStructure.h"

#include <string>
#include <vector>
#include <Eigen/Core>


class MeshCuboidAttributes {
public:
	MeshCuboidAttributes();
	MeshCuboidAttributes(const std::string _object_name);
	MeshCuboidAttributes(const MeshCuboidAttributes& _other);	// Copy constructor.
	virtual ~MeshCuboidAttributes();

	// [0 - 23] : 8 corner points.
	static const unsigned int k_num_attributes = 24;
	//static const unsigned int k_center_index = 0;
	static const unsigned int k_corner_index = 0;

	// Note:
	// Assume that z = 0 plane is the ground plane, and +y-axis is the normal direction.
	static const MyMesh::Normal k_up_direction;

	void initialize();

	void compute_attributes(const MeshCuboid *_cuboid);

	const Eigen::VectorXd get_attributes()const { return attributes_; }

	bool has_nan()const { return attributes_.hasNaN(); }

	static void get_attribute_collection_matrix(const std::list<MeshCuboidAttributes *>& _stats,
		Eigen::MatrixXd& _values);

	static bool save_attribute_collection(const std::list<MeshCuboidAttributes *>& _stats,
		const char* _filename);

private:
	std::string object_name_;
	Eigen::VectorXd attributes_;
};

class MeshCuboidFeatures {
public:
	MeshCuboidFeatures();
	MeshCuboidFeatures(const std::string _object_name);
	MeshCuboidFeatures(const MeshCuboidFeatures& _other);	// Copy constructor.
	virtual ~MeshCuboidFeatures();

	// [0 - 2] : center (local coordinate point).
	// [3 - 26] : 8 corner points  (local coordinate point).
	// [27] : center height.
	// [28 - 35] : corner heights.
	static const int k_num_features = 36;
	static const int k_corner_index = 3;
	// NOTE:
	// First 'k_num_local_coord_values' values are local coordinate values.
	static const int k_num_local_points = 9;
	static const unsigned int k_num_global_feature_values = k_num_features - 3 * k_num_local_points;

	void initialize();

	void compute_features(const MeshCuboid *_cuboid,
		Eigen::MatrixXd *_attributes_to_features_map = NULL);

	const Eigen::VectorXd get_features()const { return features_; }

	bool has_nan()const { return features_.hasNaN(); }

	static void get_feature_collection_matrix(const std::list<MeshCuboidFeatures *>& _stats,
		Eigen::MatrixXd& _values);

	static bool load_feature_collection(const char* _filename,
		std::list<MeshCuboidFeatures *>& _stats);

	static bool save_feature_collection(const char* _filename,
		const std::list<MeshCuboidFeatures *>& _stats);

private:
	std::string object_name_;
	Eigen::VectorXd features_;
};

class MeshCuboidTransformation{
public:
	MeshCuboidTransformation();
	MeshCuboidTransformation(const std::string _object_name);
	virtual ~MeshCuboidTransformation();

	void initialize();

	void compute_transformation(const MeshCuboid *_cuboid);
	Eigen::VectorXd get_transformed_features(const MeshCuboidFeatures& _other_features)const;
	Eigen::VectorXd get_transformed_features(const MeshCuboid *_other_cuboid)const;
	Eigen::VectorXd get_inverse_transformed_features(const MeshCuboidFeatures& _other_features)const;
	Eigen::VectorXd get_inverse_transformed_features(const MeshCuboid *_other_cuboid)const;

	void get_transformation(Eigen::Matrix3d &_rotation, Eigen::Vector3d &_translation) const;
	void get_inverse_transformation(Eigen::Matrix3d &_rotation, Eigen::Vector3d &_translation) const;

	void get_linear_map_transformation(Eigen::MatrixXd &_rotation, Eigen::MatrixXd &_translation)const;
	void get_linear_map_inverse_transformation(Eigen::MatrixXd &_rotation, Eigen::MatrixXd &_translation)const;
	
	static bool load_transformation_collection(const char* _filename,
		std::list<MeshCuboidTransformation *>& _stats);

	static bool save_transformation_collection(const char* _filename,
		const std::list<MeshCuboidTransformation *>& _stats);

private:
	std::string object_name_;
	Eigen::Vector3d first_translation_;
	Eigen::Matrix3d second_rotation_;
};

class MeshCuboidJointNormalRelations {
public:
	MeshCuboidJointNormalRelations();
	~MeshCuboidJointNormalRelations();

	static const int k_mat_size =
		2 * (2 * MeshCuboidFeatures::k_num_features - MeshCuboidFeatures::k_corner_index);

	bool load_joint_normal_csv(const char* _filename);
	bool save_joint_normal_csv(const char* _filename)const;

	bool load_joint_normal_dat(const char* _filename);

	static void get_pairwise_cuboid_features(
		const MeshCuboid *_cuboid_1, const MeshCuboid *_cuboid_2,
		const MeshCuboidTransformation *_transformation_1, const MeshCuboidTransformation *_transformation_2,
		Eigen::VectorXd &_pairwise_features_vec);

	static void get_pairwise_cuboid_features(
		const MeshCuboidFeatures &_features_1, const MeshCuboidFeatures &_features_2,
		const MeshCuboidTransformation *_transformation_1, const MeshCuboidTransformation *_transformation_2,
		Eigen::VectorXd &_pairwise_features_vec);

	double compute_error(const MeshCuboid *_cuboid_1, const MeshCuboid *_cuboid_2,
		const MeshCuboidTransformation *_transformation_1, const MeshCuboidTransformation *_transformation_2)const;

	// 1 is fixed and 2 is unknown.
	double compute_conditional_error(const MeshCuboid *_cuboid_1, const MeshCuboid *_cuboid_2,
		const MeshCuboidTransformation *_transformation_1)const;

	const Eigen::VectorXd &get_mean()const { return mean_; }
	const Eigen::MatrixXd &get_inv_cov()const { return inv_cov_; }

	void set_mean(const Eigen::VectorXd &_mean) { mean_ = _mean; }
	void set_inv_cov(const Eigen::MatrixXd &_inv_cov_) { inv_cov_ = _inv_cov_; }

private:
	Eigen::VectorXd mean_;
	Eigen::MatrixXd inv_cov_;
};

class MeshCuboidCondNormalRelations {
public:
	MeshCuboidCondNormalRelations();
	~MeshCuboidCondNormalRelations();

	bool load_cond_normal_csv(const char* _filename);
	bool save_cond_normal_csv(const char* _filename)const;

	bool load_cond_normal_dat(const char* _filename);

	static void get_pairwise_cuboid_features(
		const MeshCuboid *_cuboid_1, const MeshCuboid *_cuboid_2,
		const MeshCuboidTransformation *_transformation_1, const MeshCuboidTransformation *_transformation_2,
		Eigen::VectorXd &_global_features_vec_1, Eigen::VectorXd &_transformed_features_vec_12);

	static void get_pairwise_cuboid_features(
		const MeshCuboidFeatures &_features_1, const MeshCuboidFeatures &_features_2,
		const MeshCuboidTransformation *_transformation_1, const MeshCuboidTransformation *_transformation_2,
		Eigen::VectorXd &_global_features_vec_1, Eigen::VectorXd &_transformed_features_vec_12);

	double compute_error(const MeshCuboid *_cuboid_1, const MeshCuboid *_cuboid_2,
		const MeshCuboidTransformation *_transformation_1, const MeshCuboidTransformation *_transformation_2)const;

	const Eigen::MatrixXd &get_mean_A()const { return mean_A_; }
	const Eigen::VectorXd &get_mean_b()const { return mean_b_; }
	const Eigen::MatrixXd &get_inv_cov()const { return inv_cov_; }

	void set_mean_A(const Eigen::MatrixXd &_mean_A) { mean_A_ = _mean_A; }
	void set_mean_b(const Eigen::VectorXd &_mean_b) { mean_b_ = _mean_b; }
	void set_inv_cov(const Eigen::MatrixXd &_inv_cov_) { inv_cov_ = _inv_cov_; }

private:
	Eigen::MatrixXd mean_A_;
	Eigen::VectorXd mean_b_;
	Eigen::MatrixXd inv_cov_;
};

/*
class MeshCuboidPCARelations {
public:
	MeshCuboidPCARelations();
	~MeshCuboidPCARelations();

	static const unsigned int k_mat_size = 2 * (MeshCuboidFeatures::k_num_features);

	bool load_pca_csv(const char* _filename);
	bool save_pca_csv(const char* _filename)const;

	double compute_error(
		const Eigen::VectorXd &_features_vec_1,
		const Eigen::VectorXd &_features_vec_2)const;

	const Eigen::VectorXd &get_mean()const { return mean_; }
	const Eigen::MatrixXd &get_pca_bases()const { return pca_bases_; }

private:
	Eigen::VectorXd mean_;
	Eigen::MatrixXd pca_bases_;
};

class MeshCuboidCCARelations {
public:
	MeshCuboidCCARelations();
	~MeshCuboidCCARelations();

	bool load_cca_bases(const char* _filename);
	bool save_cca_bases(const char* _filename)const;

	double compute_error(
		const Eigen::VectorXd &_features_vec_1,
		const Eigen::VectorXd &_features_vec_2)const;

private:
	Eigen::VectorXd mean_1_;
	Eigen::VectorXd mean_2_;
	Eigen::VectorXd correlations_;
	Eigen::MatrixXd bases_12_;
	Eigen::MatrixXd bases_21_;
};
*/

#endif	// _MESH_CUBOID_RELATION_H_
