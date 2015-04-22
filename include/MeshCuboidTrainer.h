#ifndef _MESH_CUBOID_TRAINER_H_
#define _MESH_CUBOID_TRAINER_H_

#include "MeshCuboid.h"
#include "MeshCuboidRelation.h"
#include "MeshCuboidStructure.h"

#include <vector>
#include <string>


class MeshCuboidTrainer {
public:
	MeshCuboidTrainer();

	void clear();

	bool load_object_list(const std::string &_filename);
	bool load_features(const std::string &_filename_prefix);
	bool load_transformations(const std::string &_filename_prefix);

	void get_conflicted_labels(
		std::vector< std::list<LabelIndex> > &_cooccurrence_labels)const;

	void get_missing_label_index_groups(
		const std::list<LabelIndex> &_given_label_indices,
		std::list< std::list<LabelIndex> > &_missing_label_index_groups,
		const std::set<LabelIndex> *_ignored_label_indices = NULL)const;

	void get_joint_normal_relations(
		std::vector< std::vector<MeshCuboidJointNormalRelations *> > &_relations,
		const std::list<std::string> *_ignored_object_list = NULL)const;

	void get_cond_normal_relations(
		std::vector< std::vector<MeshCuboidCondNormalRelations *> > &_relations,
		const std::list<std::string> *_ignored_object_list = NULL)const;

	static void load_joint_normal_relations(
		const unsigned int _num_labels, const std::string _filename_prefix,
		std::vector< std::vector<MeshCuboidJointNormalRelations *> > &_relations);

	static void load_cond_normal_relations(
		const unsigned int _num_labels, const std::string _filename_prefix,
		std::vector< std::vector<MeshCuboidCondNormalRelations *> > &_relations);

protected:
	std::list<std::string> object_list_;
	std::vector< std::list<MeshCuboidFeatures *> > feature_list_;
	std::vector< std::list<MeshCuboidTransformation *> > transformation_list_;
};
#endif	// _MESH_CUBOID_TRAINER_H_