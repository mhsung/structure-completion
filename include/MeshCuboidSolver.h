#ifndef _MESH_CUBOID_SOLVER_H_
#define _MESH_CUBOID_SOLVER_H_

#include "GLViewerCore.h"
#include "MeshCuboid.h"
#include "MeshCuboidRelation.h"
#include "MeshCuboidPredictor.h"

#include <vector>
#include <string>
#include <Eigen/Core>


std::vector<int> solve_markov_random_field(
	const unsigned int _num_nodes,
	const unsigned int _num_labels,
	const Eigen::MatrixXd& _energy_mat);

Eigen::VectorXd solve_quadratic_programming(
	const Eigen::MatrixXd& _quadratic_term,
	const Eigen::VectorXd& _linear_term,
	const double _constant_term,
	Eigen::VectorXd* _init_values_vec = NULL);

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

void optimize_attributes_quadratic_once(
	const std::vector<MeshCuboid *>& _cuboids,
	const MeshCuboidPredictor &_predictor,
	const double _single_energy_term_weight);

void optimize_attributes_once(
	MeshCuboidStructure &_cuboid_structure,
	const MeshCuboidPredictor &_predictor,
	const double _single_energy_term_weight,
	const double _symmetry_energy_term_weight,
	bool _use_symmetry);

void optimize_attributes(
	MeshCuboidStructure &_cuboid_structure,
	const Real _modelview_matrix[16],
	const MeshCuboidPredictor &_predictor,
	const double _single_energy_term_weight,
	const double _symmetry_energy_term_weight,
	const unsigned int _max_num_iterations,
	const std::string _log_filename,
	GLViewerCore *_viewer,
	bool _use_symmetry);

void add_missing_cuboids_once(
	const std::vector<MeshCuboid *>& _given_cuboids,
	const std::list<LabelIndex> &_missing_label_indices,
	const MeshCuboidPredictor &_predictor,
	std::vector<MeshCuboid *>& _new_cuboids);

void add_missing_cuboids_once_simple(
	const MeshCuboid *_given_cuboids,
	const std::list<LabelIndex> &_missing_label_indices,
	const std::vector< std::vector<MeshCuboidJointNormalRelations *> > &_joint_normal_relations,
	std::vector<MeshCuboid *>& _new_cuboids);

bool add_missing_cuboids(
	MeshCuboidStructure &_cuboid_structure,
	const Real _modelview_matrix[16],
	const std::list<LabelIndex> &_missing_label_indices,
	const MeshCuboidPredictor &_predictor,
	//const std::vector< std::vector<MeshCuboidJointNormalRelations *> > &_joint_normal_relations,
	std::set<LabelIndex> &_ignored_label_indices);

//void symmetrize_cuboids(
//	MeshCuboidStructure &_cuboid_structure);

/*
MeshCuboid *test_joint_normal_training(
	const MeshCuboid *_cuboid_1, const MeshCuboid *_cuboid_2,
	const LabelIndex _label_index_1, const LabelIndex _label_index_2,
	const std::vector< std::vector<MeshCuboidJointNormalRelations> >& _relations);
*/

#endif	// _MESH_CUBOID_SOLVER_H_