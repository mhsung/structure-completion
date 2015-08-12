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
	MeshCuboidEvaluator(MeshCuboidStructure *_ground_truth_cuboid_structure);

	void evaluate_point_to_point_distances(
		const MeshCuboidStructure *_test_cuboid_structure,
		const char *_filename);

	// NOTE: Should be called only when mesh label file is already loaded.
	void evaluate_point_labeling(
		const MeshCuboidStructure *_test_cuboid_structure,
		const char *_filename);

	// NOTE: Cuboid surface points are replaced.
	void evaluate_cuboid_distance(
		const MeshCuboidStructure *_test_cuboid_structure,
		const char *_filename);
	
private:
	void evaluate_point_to_point_distances(
		const std::vector<MeshSamplePoint *> _ground_truth_sample_points,
		const std::vector<MeshSamplePoint *> _test_sample_points,
		const char *_filename, bool _record_error = false);

protected:
	MeshCuboidStructure *ground_truth_cuboid_structure_;
	const std::string mesh_name_;
	const std::string cuboid_structure_name_;

	std::map<std::string, Real> evaluation_results_;
};

#endif	// _MESH_CUBOID_EVALUATOR_H_
