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

	inline unsigned int get_cuboid_corner_variable_start_index(
		unsigned int _cuboid_index) const;
	inline unsigned int get_cuboid_corner_variable_start_index(
		unsigned int _cuboid_index, unsigned int _corner_index) const;
	inline unsigned int get_cuboid_axis_variable_start_index(
		unsigned int _cuboid_index) const;
	inline unsigned int get_cuboid_axis_variable_start_index(
		unsigned int _cuboid_index, unsigned int axis_index) const;
	inline unsigned int get_symmetry_group_variable_start_index(
		unsigned int _symmetry_group_index) const;


	Eigen::VectorXd solve_quadratic_programming_with_constraints(
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
		const unsigned int _corner_variable_start_index,
		const unsigned int _axis_variable_start_index,
		NLPFormulation &_formulation);

	void add_symmetry_group_constraints(
		const unsigned int _symmetry_group_variable_index,
		const MeshCuboidSymmetryGroup *_symmetry_group,
		NLPFormulation &_formulation);

	void add_reflection_constraints(
		const unsigned int _x_variable_index,
		const unsigned int _y_variable_index,
		const unsigned int _symmetry_group_variable_index,
		NLPFormulation &_formulation);

	void add_cuboid_reflection_constraints(
		const unsigned int _cuboid_index_1,
		const unsigned int _cuboid_index_2,
		const unsigned int _symmetry_group_variable_index,
		const unsigned int _reflection_axis_index,
		NLPFormulation &_formulation);

	// Initialization function.
	bool compute_initial_values(const Eigen::VectorXd &_input, Eigen::VectorXd &_output);
	void compute_cuboid_axis_values(Eigen::VectorXd &_values);
	void compute_symmetry_group_values(Eigen::VectorXd &_values);


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