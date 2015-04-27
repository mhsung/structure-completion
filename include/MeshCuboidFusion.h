#ifndef _MESH_CUBOID_FUSION_H_
#define _MESH_CUBOID_FUSION_H_

#include "MeshCuboidStructure.h"

void reconstruct_fusion_simple(
	const MeshCuboidStructure &_symmetry_reconstruction,
	const MeshCuboidStructure &_database_reconstruction,
	MeshCuboidStructure &_output_cuboid_structure);

void reconstruct_fusion(const char *_mesh_filepath,
	const double *_snapshot_modelview_matrix,
	const double *_occlusion_modelview_matrix,
	const MeshCuboidStructure &_original_cuboid_structure,
	const MeshCuboidStructure &_symmetry_reconstruction,
	const MeshCuboidStructure &_database_reconstruction,
	MeshCuboidStructure &_output_cuboid_structure);

#endif	// _MESH_CUBOID_FUSION_H_