#ifndef _MESH_CUBOID_SOLVER_H_
#define _MESH_CUBOID_SOLVER_H_

#include "MeshCuboid.h"
#include "MeshCuboidRelation.h"
#include "MeshCuboidPredictor.h"

#include <vector>
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
	Eigen::MatrixXd &_potential_mat);

void recognize_labels_and_axes_configurations(
	MeshCuboidStructure &_cuboid_structure,
	const MeshCuboidPredictor &_predictor,
	const char* _log_filename);

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
	const char* _log_filename,
	const unsigned int _max_num_iterations = 30,
	QGLWidget *_viewer = NULL);

bool add_missing_cuboids(
	MeshCuboidStructure &_cuboid_structure,
	// NOTE:
	// 'MeshCuboidCondNormalRelationPredictor' Only.
	const MeshCuboidCondNormalRelationPredictor &_predictor);

MeshCuboid *test_joint_normal_training(
	const MeshCuboid *_cuboid_1, const MeshCuboid *_cuboid_2,
	const LabelIndex _label_index_1, const LabelIndex _label_index_2,
	const std::vector< std::vector<MeshCuboidJointNormalRelations> >& _relations);

#endif	// _MESH_CUBOID_SOLVER_H_