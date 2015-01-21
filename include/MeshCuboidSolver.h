#ifndef _MESH_CUBOID_SOLVER_H_
#define _MESH_CUBOID_SOLVER_H_

#include "MeshCuboid.h"
#include "MeshCuboidRelation.h"
#include "MeshCuboidPredictor.h"

#include <vector>
#include <string>
#include <QtOpenGL/qgl.h>

void update_cuboid_surface_points(
	MeshCuboidStructure &_cuboid_structure,
	const Real _modelview_matrix[16]);

void segment_sample_points(
	MeshCuboidStructure &_cuboid_structure);

void compute_labels_and_axes_configuration_potentials(
	const std::vector<Label>& _labels,
	const std::vector<MeshCuboid *>& _cuboids,
	const MeshCuboidPredictor &_predictor,
	Eigen::MatrixXd &_potential_mat,
	const std::vector< std::list<LabelIndex> > *_label_symmetries,
	bool _add_dummy_label = false);

void recognize_labels_and_axes_configurations(
	MeshCuboidStructure &_cuboid_structure,
	const MeshCuboidPredictor &_predictor,
	const std::string _log_filename,
	bool _use_symmetry_info = false,
	bool _add_dummy_label = false);

void test_recognize_labels_and_axes_configurations(
	const std::vector<Label> &_labels,
	const std::vector<MeshCuboid *>& _cuboids,
	const MeshCuboidPredictor &_predictor);

void get_optimization_formulation(
	const std::vector<MeshCuboid *>& _cuboids,
	const MeshCuboidPredictor &_predictor,
	Eigen::VectorXd &_init_values,
	Eigen::MatrixXd &_single_quadratic_term, Eigen::MatrixXd &_pair_quadratic_term,
	Eigen::VectorXd &_single_linear_term, Eigen::VectorXd &_pair_linear_term,
	double &_single_constant_term, double &_pair_constant_term,
	double &_single_total_energy, double &_pair_total_energy);

void get_optimization_error(
	const std::vector<MeshCuboid *>& _cuboids,
	const MeshCuboidPredictor &_predictor,
	double &_single_total_energy, double &_pair_total_energy);

void optimize_attributes_once(
	const std::vector<MeshCuboid *>& _cuboids,
	const MeshCuboidPredictor &_predictor,
	const double _quadprog_ratio);

void optimize_attributes(
	MeshCuboidStructure &_cuboid_structure,
	const Real _modelview_matrix[16],
	const MeshCuboidPredictor &_predictor,
	const double _quadprog_ratio,
	const std::string _log_filename,
	const unsigned int _max_num_iterations = 10,
	QGLWidget *_viewer = NULL);

void add_missing_cuboids_once(
	const std::vector<MeshCuboid *>& _given_cuboids,
	const std::list<LabelIndex> &_missing_label_indices,
	const MeshCuboidPredictor &_predictor,
	std::vector<MeshCuboid *>& _new_cuboids);

void add_missing_cuboids(
	MeshCuboidStructure &_cuboid_structure,
	const Real _modelview_matrix[16],
	const std::list<LabelIndex> &_missing_label_indices,
	const MeshCuboidPredictor &_predictor);

bool evaluate_segmentation(
	const MeshCuboidStructure &_cuboid_structure,
	std::vector<LabelIndex> _sample_point_label_indices,
	const std::string _mesh_name,
	const std::string _stats_filename);

/*
MeshCuboid *test_joint_normal_training(
	const MeshCuboid *_cuboid_1, const MeshCuboid *_cuboid_2,
	const LabelIndex _label_index_1, const LabelIndex _label_index_2,
	const std::vector< std::vector<MeshCuboidJointNormalRelations> >& _relations);
*/

#endif	// _MESH_CUBOID_SOLVER_H_