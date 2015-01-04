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

protected:
	std::list<std::string> object_list_;
	std::vector< std::list<MeshCuboidFeatures *> > feature_list_;
	std::vector< std::list<MeshCuboidTransformation *> > transformation_list_;
};

// Use joint normal relations for binary terms.
class MeshCuboidJointNormalRelationTrainer : public MeshCuboidTrainer{
public:
	MeshCuboidJointNormalRelationTrainer();

	void compute_relations(
		std::vector< std::vector<MeshCuboidJointNormalRelations *> > &_relations,
		const std::list<std::string> *_ignored_object_list = NULL)const;
};
#endif	// _MESH_CUBOID_TRAINER_H_