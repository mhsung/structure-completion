#ifndef _MESH_CUBOID_NON_LINEAR_SOLVER_H_
#define _MESH_CUBOID_NON_LINEAR_SOLVER_H_

#include "MeshCuboid.h"
#include "MeshCuboidRelation.h"
#include "MeshCuboidPredictor.h"
#include "MeshCuboidSymmetryGroup.h"

// IPOPT.
#include "NLPFormulation.h"
#include "NLPEigenQuadFunction.h"
#include "NLPVectorExpression.h"
#include "IPOPTSolver.h"
#include "IpIpoptApplication.hpp"

#include "ANN/ANN.h"
#include <Eigen/Core>


class MeshCuboidNonLinearSolver
{
public:
	MeshCuboidNonLinearSolver(
		const std::vector<MeshCuboid *>& _cuboids,
		const std::vector<MeshCuboidReflectionSymmetryGroup *>& _reflection_symmetry_groups,
		const std::vector<MeshCuboidRotationSymmetryGroup *>& _rotation_symmetry_groups,
		const Real _neighbor_distance, const Real _symmetry_energy_term_weight);
	~MeshCuboidNonLinearSolver();


	inline unsigned int num_total_variables() const;
	inline unsigned int num_total_cuboid_corner_variables() const;
	inline unsigned int num_total_cuboid_axis_variables() const;
	inline unsigned int num_total_reflection_symmetry_group_variables() const;
	inline unsigned int num_total_rotation_symmetry_group_variables() const;

	static NLPVectorExpression create_vector_variable(
		const std::pair<Index, Index>& _index_size_pair);

	// Pair: (index, size)
	inline std::pair<Index, Index> get_cuboid_corner_variable_index_size(
		unsigned int _cuboid_index, unsigned int _corner_index) const;
	inline std::pair<Index, Index> get_cuboid_axis_variable_index_size(
		unsigned int _cuboid_index, unsigned int _axis_index) const;
	inline std::pair<Index, Index> get_reflection_symmetry_group_variable_n_index_size(
		unsigned int _symmetry_group_index) const;
	inline std::pair<Index, Index> get_reflection_symmetry_group_variable_t_index_size(
		unsigned int _symmetry_group_index) const;
	inline std::pair<Index, Index> get_rotation_symmetry_group_variable_n_index_size(
		unsigned int _symmetry_group_index) const;
	inline std::pair<Index, Index> get_rotation_symmetry_group_variable_t_index_size(
		unsigned int _symmetry_group_index) const;

	inline NLPVectorExpression create_cuboid_corner_variable(
		unsigned int _cuboid_index, unsigned int _corner_index) const;
	inline NLPVectorExpression create_cuboid_axis_variable(
		unsigned int _cuboid_index, unsigned int _axis_index) const;
	inline NLPVectorExpression create_reflection_symmetry_group_variable_n(
		unsigned int _symmetry_group_index) const;
	inline NLPVectorExpression create_reflection_symmetry_group_variable_t(
		unsigned int _symmetry_group_index) const;
	inline NLPVectorExpression create_rotation_symmetry_group_variable_n(
		unsigned int _symmetry_group_index) const;
	inline NLPVectorExpression create_rotation_symmetry_group_variable_t(
		unsigned int _symmetry_group_index) const;


	void optimize(
		const Eigen::MatrixXd& _cuboid_quadratic_term,
		const Eigen::VectorXd& _cuboid_linear_term,
		const double _cuboid_constant_term,
		Eigen::VectorXd* _init_values_vec = NULL);


private:
	// Core functions.
	void create_energy_functions(
		const Eigen::MatrixXd &_quadratic_term,
		const Eigen::VectorXd &_linear_term,
		const double _constant_term,
		std::vector<NLPFunction *> &_functions);

	void add_constraints(NLPFormulation &_formulation);

	bool compute_initial_values(const Eigen::VectorXd &_input, Eigen::VectorXd &_output);

	void update(const std::vector< Number >& _values);


	// Energy functions.
	NLPFunction *create_reflection_symmetry_group_energy_function();

	void create_reflection_symmetry_group_energy_function(
		const unsigned int _symmetry_group_index,
		const std::vector<ANNpointArray>& _cuboid_ann_points,
		const std::vector<ANNkd_tree *>& _cuboid_ann_kd_tree,
		Eigen::MatrixXd& _quadratic_term,
		Eigen::VectorXd& _linear_term,
		double &_constant_term);

	NLPFunction *create_rotation_symmetry_group_energy_function();

	void create_rotation_symmetry_group_energy_function(
		const unsigned int _symmetry_group_index,
		const std::vector<ANNpointArray>& _cuboid_ann_points,
		const std::vector<ANNkd_tree *>& _cuboid_ann_kd_tree,
		NLPExpression &_expression);

	void create_cuboid_sample_point_ann_trees(
		std::vector<ANNpointArray>& _cuboid_ann_points,
		std::vector<ANNkd_tree *>& _cuboid_ann_kd_tree) const;

	void delete_cuboid_sample_point_ann_trees(
		std::vector<ANNpointArray>& _cuboid_ann_points,
		std::vector<ANNkd_tree *>& _cuboid_ann_kd_tree) const;


	// Constraint functions.
	void add_cuboid_constraints(NLPFormulation &_formulation);
	void add_reflection_symmetry_group_constraints(NLPFormulation &_formulation);
	void add_rotation_symmetry_group_constraints(NLPFormulation &_formulation);

	void add_cuboid_constraints(
		const unsigned int _cuboid_index,
		NLPFormulation &_formulation);

	void add_reflection_symmetry_group_constraints(
		const unsigned int _symmetry_group_index,
		NLPFormulation &_formulation);

	void add_reflection_constraints(
		const NLPVectorExpression& _x_variable,
		const NLPVectorExpression& _y_variable,
		const NLPVectorExpression& _n_variable,
		const NLPVectorExpression& _t_variable,
		NLPFormulation &_formulation);

	void add_single_cuboid_reflection_constraints(
		const unsigned int _cuboid_index,
		const unsigned int _symmetry_group_index,
		const unsigned int _reflection_axis_index,
		NLPFormulation &_formulation);

	void add_pair_cuboid_reflection_constraints(
		const unsigned int _cuboid_index_1,
		const unsigned int _cuboid_index_2,
		const unsigned int _symmetry_group_index,
		const unsigned int _reflection_axis_index,
		NLPFormulation &_formulation);

	void add_rotation_symmetry_group_constraints(
		const unsigned int _symmetry_group_index,
		NLPFormulation &_formulation);


	// Initialization function.
	void compute_cuboid_axis_values(Eigen::VectorXd &_values);
	void compute_reflection_symmetry_group_values(Eigen::VectorXd &_values);
	void compute_rotation_symmetry_group_values(Eigen::VectorXd &_values);


	// Update function.
	void update_cuboids(const std::vector< Number >& _values);
	void update_reflection_symmetry_groups(const std::vector< Number >& _values);
	void update_rotation_symmetry_groups(const std::vector< Number >& _values);


private:
	const std::vector<MeshCuboid *>& cuboids_;
	const std::vector<MeshCuboidReflectionSymmetryGroup *>& reflection_symmetry_groups_;
	const std::vector<MeshCuboidRotationSymmetryGroup *>& rotation_symmetry_groups_;
	const Real neighbor_distance_;
	const double symmetry_energy_term_weight_;

	unsigned int num_cuboids_;
	unsigned int num_reflection_symmetry_groups_;
	unsigned int num_rotation_symmetry_groups_;

	unsigned int cuboid_corner_variable_start_index_;
	unsigned int cuboid_axis_variable_start_index_;
	unsigned int reflection_symmetry_group_variable_start_index_;
	unsigned int rotation_symmetry_group_variable_start_index_;

	const unsigned int num_cuboid_corner_variables_;
	const unsigned int num_cuboid_axis_variables_;
	const unsigned int num_reflection_symmetry_group_variables_;
	const unsigned int num_rotation_symmetry_group_variables_;
};

#endif	// _MESH_CUBOID_NON_LINEAR_SOLVER_H_