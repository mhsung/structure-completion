#ifndef _MESH_CUBOID_PREDICTOR_H_
#define _MESH_CUBOID_PREDICTOR_H_

#include "MeshCuboid.h"
#include "MeshCuboidRelation.h"
#include "MeshCuboidStructure.h"

#include <vector>

// Recognition:
// Recognize primitive labels and local coordinates.
class MeshCuboidPredictor {
public:
	MeshCuboidPredictor(unsigned int _num_labels);

	virtual void get_missing_label_indices(
		const std::list<LabelIndex> &_given_label_indices,
		std::list<LabelIndex> &_missing_label_indices)const;

	virtual Real get_single_potential(const MeshCuboid *_cuboid,
		const MeshCuboidAttributes *_attributes,
		const MeshCuboidTransformation *_transformation,
		const LabelIndex _label_index)const;

	virtual Real get_pair_potential(const MeshCuboid *_cuboid_1, const MeshCuboid *_cuboid_2,
		const MeshCuboidAttributes *_attributes_1, const MeshCuboidAttributes *_attributes_2,
		const MeshCuboidTransformation *_transformation_1, const MeshCuboidTransformation *_transformation_2,
		const LabelIndex _label_index_1, const LabelIndex _label_index_2)const;

	virtual void get_single_quadratic_form(MeshCuboid *_cuboid, const unsigned int _cuboid_index,
		Eigen::MatrixXd &_quadratic_term, Eigen::VectorXd &_linear_term, double& _constant_term)const;

	virtual Real get_pair_quadratic_form(const MeshCuboid *_cuboid_1, const MeshCuboid *_cuboid_2,
		const unsigned int _cuboid_index_1, const unsigned int _cuboid_index_2,
		const LabelIndex _label_index_1, const LabelIndex _label_index_2,
		Eigen::MatrixXd &_quadratic_term, Eigen::VectorXd &_linear_term, double& _constant_term)const;

	// 1 is fixed and 2 is unknown.
	virtual Real get_pair_conditional_quadratic_form(const MeshCuboid *_cuboid_1, const MeshCuboid *_cuboid_2,
		const unsigned int _cuboid_index_1, const unsigned int _cuboid_index_2,
		const LabelIndex _label_index_1, const LabelIndex _label_index_2,
		Eigen::MatrixXd &_quadratic_term, Eigen::VectorXd &_linear_term, double& _constant_term)const;

protected:
	const unsigned int num_labels_;
};

// Use joint normal relations for binary terms.
class MeshCuboidJointNormalRelationPredictor : public MeshCuboidPredictor{
public:
	MeshCuboidJointNormalRelationPredictor(
		const std::vector< std::vector<MeshCuboidJointNormalRelations *> > &_relations);

	virtual void get_missing_label_indices(
		const std::list<LabelIndex> &_given_label_indices,
		std::list<LabelIndex> &_missing_label_indices)const;

	virtual Real get_pair_potential(const MeshCuboid *_cuboid_1, const MeshCuboid *_cuboid_2,
		const MeshCuboidAttributes *_attributes_1, const MeshCuboidAttributes *_attributes_2,
		const MeshCuboidTransformation *_transformation_1, const MeshCuboidTransformation *_transformation_2,
		const LabelIndex _label_index_1, const LabelIndex _label_index_2)const;

	virtual Real get_pair_quadratic_form(const MeshCuboid *_cuboid_1, const MeshCuboid *_cuboid_2,
		const LabelIndex _label_index_1, const LabelIndex _label_index_2,
		const unsigned int _cuboid_index_1, const unsigned int _cuboid_index_2,
		Eigen::MatrixXd &_quadratic_term, Eigen::VectorXd &_linear_term, double& _constant_term)const;

	// 1 is fixed and 2 is unknown.
	virtual Real get_pair_conditional_quadratic_form(const MeshCuboid *_cuboid_1, const MeshCuboid *_cuboid_2,
		const LabelIndex _label_index_1, const LabelIndex _label_index_2,
		const unsigned int _cuboid_index_1, const unsigned int _cuboid_index_2,
		Eigen::MatrixXd &_quadratic_term, Eigen::VectorXd &_linear_term, double& _constant_term)const;

private:
	const std::vector< std::vector<MeshCuboidJointNormalRelations *> > &relations_;
};

// Use conditional normal relations for binary terms.
class MeshCuboidCondNormalRelationPredictor : public MeshCuboidPredictor{
public:
	MeshCuboidCondNormalRelationPredictor(
		const std::vector< std::vector<MeshCuboidCondNormalRelations *> > &_relations);

	virtual void get_missing_label_indices(
		const std::list<LabelIndex> &_given_label_indices,
		std::list<LabelIndex> &_missing_label_indices)const;

	virtual Real get_pair_potential(const MeshCuboid *_cuboid_1, const MeshCuboid *_cuboid_2,
		const MeshCuboidAttributes *_attributes_1, const MeshCuboidAttributes *_attributes_2,
		const MeshCuboidTransformation *_transformation_1, const MeshCuboidTransformation *_transformation_2,
		const LabelIndex _label_index_1, const LabelIndex _label_index_2)const;

	virtual Real get_pair_quadratic_form(const MeshCuboid *_cuboid_1, const MeshCuboid *_cuboid_2,
		const LabelIndex _label_index_1, const LabelIndex _label_index_2,
		const unsigned int _cuboid_index_1, const unsigned int _cuboid_index_2,
		Eigen::MatrixXd &_quadratic_term, Eigen::VectorXd &_linear_term, double& _constant_term)const;

private:
	const std::vector< std::vector<MeshCuboidCondNormalRelations *> > &relations_;
};

/*
// Use PCA relations for binary terms.
class MeshCuboidPCARelationPredictor : public MeshCuboidPredictor{
public:
	MeshCuboidPCARelationPredictor(
		const std::vector< std::vector<MeshCuboidPCARelations> > &_relations);

	virtual Real get_pair_potential(const MeshCuboid *_cuboid_1, const MeshCuboid *_cuboid_2,
		const MeshCuboidAttributes *_attributes_1, const MeshCuboidAttributes *_attributes_2,
		const MeshCuboidTransformation *_transformation_1, const MeshCuboidTransformation *_transformation_2,
		const LabelIndex _label_index_1, const LabelIndex _label_index_2)const;

	virtual Real get_pair_quadratic_form(const MeshCuboid *_cuboid_1, const MeshCuboid *_cuboid_2,
		const LabelIndex _label_index_1, const LabelIndex _label_index_2,
		const unsigned int _cuboid_index_1, const unsigned int _cuboid_index_2,
		Eigen::MatrixXd &_quadratic_term, Eigen::VectorXd &_linear_term, double& _constant_term)const;

private:
	const std::vector< std::vector<MeshCuboidPCARelations> > &relations_;
};

// Use CCA relations for binary terms.
class MeshCuboidCCARelationPredictor : public MeshCuboidPredictor{
public:
	MeshCuboidCCARelationPredictor(
		const std::vector< std::vector< std::vector<MeshCuboidCCARelations> > >& _relations);

	virtual Real get_pair_potential(const MeshCuboid *_cuboid_1, const MeshCuboid *_cuboid_2,
		const MeshCuboidAttributes *_attributes_1, const MeshCuboidAttributes *_attributes_2,
		const MeshCuboidTransformation *_transformation_1, const MeshCuboidTransformation *_transformation_2,
		const LabelIndex _label_index_1, const LabelIndex _label_index_2)const;

private:
	const std::vector< std::vector< std::vector<MeshCuboidCCARelations> > >& relations_;
};

// Use manual relations for both unary and binary terms.
class MeshCuboidManualRelationPredictor : public MeshCuboidPredictor{
public:
	MeshCuboidManualRelationPredictor(
		const std::vector<MeshCuboidStats> &_single_stats,
		const std::vector< std::vector<MeshCuboidStats> > &_pair_stats);

	virtual Real get_single_potential(const MeshCuboid *_cuboid,
		const MeshCuboidAttributes *_attributes,
		const MeshCuboidTransformation *_transformation,
		const LabelIndex _label_index)const;

	virtual Real get_pair_potential(const MeshCuboid *_cuboid_1, const MeshCuboid *_cuboid_2,
		const MeshCuboidAttributes *_attributes_1, const MeshCuboidAttributes *_attributes_2,
		const MeshCuboidTransformation *_transformation_1, const MeshCuboidTransformation *_transformation_2,
		const LabelIndex _label_index_1, const LabelIndex _label_index_2)const;

private:
	const std::vector<MeshCuboidStats> &single_stats_;
	const std::vector< std::vector<MeshCuboidStats> > &pair_stats_;
};
*/

#endif	// _MESH_CUBOID_PREDICTOR_H_