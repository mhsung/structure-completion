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

#include <Eigen/Core>


class MeshCuboidNonLinearSolver
{
public:
	MeshCuboidNonLinearSolver(
		const std::vector<MeshCuboid *>& _cuboids,
		const std::vector<MeshCuboidSymmetryGroup *>& _symmetry_groups);
	~MeshCuboidNonLinearSolver();

	inline unsigned int num_total_variables() const;
	inline unsigned int num_total_cuboid_corner_variables() const;
	inline unsigned int num_total_cuboid_axis_variables() const;
	inline unsigned int num_total_symmetry_gtoup_variables() const;

	static NLPVectorExpression create_vector_variable(
		const std::pair<Index, Index>& _index_size_pair);

	// Pair: (index, size)
	inline std::pair<Index, Index> get_cuboid_corner_variable_index_size(
		unsigned int _cuboid_index, unsigned int _corner_index) const;
	inline std::pair<Index, Index> get_cuboid_axis_variable_index_size(
		unsigned int _cuboid_index, unsigned int _axis_index) const;
	inline std::pair<Index, Index> get_symmetry_group_variable_n_index_size(
		unsigned int _symmetry_group_index) const;
	inline std::pair<Index, Index> get_symmetry_group_variable_t_index_size(
		unsigned int _symmetry_group_index) const;

	inline NLPVectorExpression get_cuboid_corner_variable(
		unsigned int _cuboid_index, unsigned int _corner_index) const;
	inline NLPVectorExpression get_cuboid_axis_variable(
		unsigned int _cuboid_index, unsigned int _axis_index) const;
	inline NLPVectorExpression get_symmetry_group_variable_n(
		unsigned int _symmetry_group_index) const;
	inline NLPVectorExpression get_symmetry_group_variable_t(
		unsigned int _symmetry_group_index) const;


	void optimize(
		const Eigen::MatrixXd& _quadratic_term,
		const Eigen::VectorXd& _linear_term,
		const double _constant_term,
		Eigen::VectorXd* _init_values_vec = NULL);


private:
	NLPEigenQuadFunction* create_quadratic_function(
		const Eigen::MatrixXd& _quadratic_term,
		const Eigen::VectorXd& _linear_term,
		const double _constant_term);

	// Constraint functions.
	void add_cuboid_constraints(NLPFormulation &_formulation);
	void add_symmetry_group_constraints(NLPFormulation &_formulation);

	void add_cuboid_constraints(
		const unsigned int _cuboid_index,
		NLPFormulation &_formulation);

	void add_symmetry_group_constraints(
		const unsigned int _symmetry_group_index,
		NLPFormulation &_formulation);

	void add_reflection_constraints(
		const NLPVectorExpression& _x_variable,
		const NLPVectorExpression& _y_variable,
		const NLPVectorExpression& _n_variable,
		const NLPVectorExpression& _t_variable,
		NLPFormulation &_formulation);

	void add_cuboid_reflection_constraints(
		const unsigned int _cuboid_index_1,
		const unsigned int _cuboid_index_2,
		const unsigned int _symmetry_group_index,
		const unsigned int _reflection_axis_index,
		NLPFormulation &_formulation);

	// Initialization function.
	bool compute_initial_values(const Eigen::VectorXd &_input, Eigen::VectorXd &_output);
	void compute_cuboid_axis_values(Eigen::VectorXd &_values);
	void compute_symmetry_group_values(Eigen::VectorXd &_values);

	// Update function.
	void update_cuboids(const std::vector< Number >& _values);
	void update_symmetry_groups(const std::vector< Number >& _values);


private:
	const std::vector<MeshCuboid *>& cuboids_;
	const std::vector<MeshCuboidSymmetryGroup *>& symmetry_groups_;

	unsigned int num_cuboids_;
	int num_symmetry_groups_;

	unsigned int cuboid_corner_variable_start_index_;
	unsigned int cuboid_axis_variable_start_index_;
	unsigned int symmetry_group_variable_start_index_;

	const unsigned int num_cuboid_corner_variables_;
	const unsigned int num_cuboid_axis_variables_;
	const unsigned int num_symmetry_group_variables_;

	Eigen::MatrixXd quadratic_term_;
	Eigen::VectorXd linear_term_;
	double constant_term_;
};

#endif	// _MESH_CUBOID_NON_LINEAR_SOLVER_H_