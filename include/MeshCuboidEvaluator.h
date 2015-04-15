#ifndef _MESH_CUBOID_EVALUATOR_H_
#define _MESH_CUBOID_EVALUATOR_H_

#include "MeshCuboid.h"
#include "MeshCuboidRelation.h"
#include "MeshCuboidStructure.h"

#include <map>
#include <vector>
#include <string>


class MeshCuboidEvaluator {
public:
	MeshCuboidEvaluator(
		const MeshCuboidStructure *_ground_truth_cuboid_structure);
		//const std::string _mesh_name,
		//const std::string _cuboid_structure_name);

	bool save_evaluate_results(
		const MeshCuboidStructure *_test_cuboid_structure,
		const std::string _filename, bool _verbose = true);
	
private:
	void evaluate_all(
		const MeshCuboidStructure *_test_cuboid_structure);

	//void evaluate_segmentation(
	//	const MeshCuboidStructure *_test_cuboid_structure);

	//void evaluate_cuboid_distance(
	//	const MeshCuboidStructure *_test_cuboid_structure);

	void evaluate_point_to_point_distances(
		const MeshCuboidStructure *_test_cuboid_structure,
		const std::string _filename);


protected:
	const MeshCuboidStructure *ground_truth_cuboid_structure_;
	const std::string mesh_name_;
	const std::string cuboid_structure_name_;

	std::map<std::string, Real> evaluation_results_;
};

#endif	// _MESH_CUBOID_EVALUATOR_H_