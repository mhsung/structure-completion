#include "MeshViewerCore.h"
#include "MeshCuboidEvaluator.h"
#include "MeshCuboidFusion.h"
#include "MeshCuboidParameters.h"
#include "simplerandom.h"

#include <sstream>
#include <QDir>
#include <QFileInfo>


void MeshViewerCore::reconstruct(
	const char *_mesh_filepath,
	const GLdouble *_snapshot_modelview_matrix,
	const GLdouble *_occlusion_modelview_matrix,
	const char *_output_file_prefix)
{
	bool ret;
	std::stringstream output_filename_sstr;

	QFileInfo file_info(_mesh_filepath);
	std::string mesh_name(file_info.baseName().toLocal8Bit());
	std::string dense_sample_filepath = FLAGS_data_root_path + FLAGS_dense_sample_path
		+ std::string("/") + mesh_name + std::string(".pts");


	updateGL();
	output_filename_sstr.clear(); output_filename_sstr.str("");
	output_filename_sstr << _output_file_prefix;
	snapshot(output_filename_sstr.str().c_str());

	output_filename_sstr.clear(); output_filename_sstr.str("");
	output_filename_sstr << _output_file_prefix << std::string(".arff");
	cuboid_structure_.save_cuboids(output_filename_sstr.str());


	//
	MyMesh ground_truth_mesh;
	MeshCuboidStructure ground_truth_cuboid_structure(&ground_truth_mesh);

	ret = true;
	ret = ret & ground_truth_cuboid_structure.load_labels((FLAGS_data_root_path +
		FLAGS_label_info_path + FLAGS_label_info_filename).c_str());
	ret = ret & ground_truth_cuboid_structure.load_label_symmetries((FLAGS_data_root_path +
		FLAGS_label_info_path + FLAGS_label_symmetry_info_filename).c_str());

	if (!ret)
	{
		do {
			std::cout << "Error: Cannot open label information files.";
			std::cout << '\n' << "Press the Enter key to continue.";
		} while (std::cin.get() != '\n');
	}

	ret = load_object_info(ground_truth_mesh, ground_truth_cuboid_structure,
		_mesh_filepath, LoadGroundTruthData);
	assert(ret);
	MeshCuboidEvaluator evaluator(&ground_truth_cuboid_structure);

	setDrawMode(COLORED_POINT_SAMPLES);
	MeshCuboidStructure cuboid_structure_copy(cuboid_structure_);


	// 1. Reconstruction using symmetry.
	reconstruct_symmetry_prior(_mesh_filepath, _snapshot_modelview_matrix, _occlusion_modelview_matrix);

	// NOTE:
	// The label of reconstructed points are recorded as confidence values.
	cuboid_structure_.set_sample_point_label_confidence_using_cuboids();

	output_filename_sstr.clear(); output_filename_sstr.str("");
	output_filename_sstr << _output_file_prefix << std::string("_symmetry");
	evaluator.evaluate_point_to_point_distances(&cuboid_structure_, output_filename_sstr.str().c_str());
	cuboid_structure_.save_sample_points_to_ply(output_filename_sstr.str().c_str());
	cuboid_structure_.save_sample_points((output_filename_sstr.str() + std::string(".pts")).c_str());
	cuboid_structure_.save_sample_point_labels((output_filename_sstr.str() + std::string("_label.arff")).c_str());
	MeshCuboidStructure symmetry_reconstruction(cuboid_structure_);
	
	updateGL();
	output_filename_sstr.clear(); output_filename_sstr.str("");
	output_filename_sstr << _output_file_prefix << std::string("_symmetry_accuracy");
	snapshot(output_filename_sstr.str().c_str());
	
	cuboid_structure_ = ground_truth_cuboid_structure;
	updateGL();
	output_filename_sstr.clear(); output_filename_sstr.str("");
	output_filename_sstr << _output_file_prefix << std::string("_symmetry_completeness");
	snapshot(output_filename_sstr.str().c_str());
	

	// 2. Reconstruction using database.
	cuboid_structure_ = cuboid_structure_copy;
	reconstruct_database_prior(_mesh_filepath);

	// NOTE:
	// The label of reconstructed points are recorded as confidence values.
	cuboid_structure_.set_sample_point_label_confidence_using_cuboids();

	output_filename_sstr.clear(); output_filename_sstr.str("");
	output_filename_sstr << _output_file_prefix << std::string("_database");
	evaluator.evaluate_point_to_point_distances(&cuboid_structure_, output_filename_sstr.str().c_str());
	cuboid_structure_.save_sample_points_to_ply(output_filename_sstr.str().c_str());
	cuboid_structure_.save_sample_points((output_filename_sstr.str() + std::string(".pts")).c_str());
	cuboid_structure_.save_sample_point_labels((output_filename_sstr.str() + std::string("_label.arff")).c_str());
	MeshCuboidStructure database_reconstruction(cuboid_structure_);

	updateGL();
	output_filename_sstr.clear(); output_filename_sstr.str("");
	output_filename_sstr << _output_file_prefix << std::string("_database_accuracy");
	snapshot(output_filename_sstr.str().c_str());

	cuboid_structure_ = ground_truth_cuboid_structure;
	updateGL();
	output_filename_sstr.clear(); output_filename_sstr.str("");
	output_filename_sstr << _output_file_prefix << std::string("_database_completeness");
	snapshot(output_filename_sstr.str().c_str());
		

	// 3. Fusion.
	cuboid_structure_ = cuboid_structure_copy;
	//MeshViewerCore::reconstruct_fusion_simple(symmetry_reconstruction, database_reconstruction);

	MyMesh original_mesh;
	MeshCuboidStructure original_cuboid_structure(&original_mesh);
	ret = load_object_info(original_mesh, original_cuboid_structure, _mesh_filepath, LoadSamplePoints);

	reconstruct_fusion(_mesh_filepath,
		_snapshot_modelview_matrix, _occlusion_modelview_matrix,
		original_cuboid_structure, symmetry_reconstruction, database_reconstruction, cuboid_structure_);
	

	// NOTE:
	// The label of reconstructed points are recorded as confidence values.
	cuboid_structure_.set_sample_point_label_confidence_using_cuboids();

	output_filename_sstr.clear(); output_filename_sstr.str("");
	output_filename_sstr << _output_file_prefix << std::string("_fusion");
	evaluator.evaluate_point_to_point_distances(&cuboid_structure_, output_filename_sstr.str().c_str());
	cuboid_structure_.save_sample_points_to_ply(output_filename_sstr.str().c_str());
	cuboid_structure_.save_sample_points((output_filename_sstr.str() + std::string(".pts")).c_str());
	cuboid_structure_.save_sample_point_labels((output_filename_sstr.str() + std::string("_label.arff")).c_str());

	updateGL();
	output_filename_sstr.clear(); output_filename_sstr.str("");
	output_filename_sstr << _output_file_prefix << std::string("_fusion_accuracy");
	snapshot(output_filename_sstr.str().c_str());

	cuboid_structure_ = ground_truth_cuboid_structure;
	updateGL();
	output_filename_sstr.clear(); output_filename_sstr.str("");
	output_filename_sstr << _output_file_prefix << std::string("_fusion_completeness");
	snapshot(output_filename_sstr.str().c_str());

	cuboid_structure_ = cuboid_structure_copy;

	setDrawMode(CUSTOM_VIEW);
}

void MeshViewerCore::reconstruct_symmetry_prior(const char *_mesh_filepath,
	const GLdouble *_snapshot_modelview_matrix,
	const GLdouble *_occlusion_modelview_matrix)
{
	QFileInfo file_info(_mesh_filepath);
	std::string mesh_name(file_info.baseName().toLocal8Bit());
	std::string dense_sample_filepath = FLAGS_data_root_path + FLAGS_dense_sample_path
		+ std::string("/") + mesh_name + std::string(".pts");

	bool ret = cuboid_structure_.load_dense_sample_points(dense_sample_filepath.c_str());
	assert(ret);
	set_modelview_matrix(_occlusion_modelview_matrix, false);
	remove_occluded_points();
	set_modelview_matrix(_snapshot_modelview_matrix);

	cuboid_structure_.copy_sample_points_to_symmetric_position();
}

void MeshViewerCore::reconstruct_database_prior(
	const char *_mesh_filepath,
	const std::vector<LabelIndex> *_reconstructed_label_indices)
{
	MyMesh example_mesh;
	MeshCuboidStructure example_cuboid_structure(&example_mesh);

	bool ret = true;
	ret = ret & example_cuboid_structure.load_labels((FLAGS_data_root_path +
		FLAGS_label_info_path + FLAGS_label_info_filename).c_str());
	ret = ret & example_cuboid_structure.load_label_symmetries((FLAGS_data_root_path +
		FLAGS_label_info_path + FLAGS_label_symmetry_info_filename).c_str());

	// Load symmetry groups.
	ret = ret & example_cuboid_structure.load_symmetry_groups((FLAGS_data_root_path +
		FLAGS_label_info_path + FLAGS_symmetry_group_info_filename).c_str());

	if (!ret)
	{
		do {
			std::cout << "Error: Cannot open label information files.";
			std::cout << '\n' << "Press the Enter key to continue.";
		} while (std::cin.get() != '\n');
	}

	unsigned int num_labels = cuboid_structure_.num_labels();
	assert(num_labels == example_cuboid_structure.num_labels());


	// For every file in the base path.
	QDir input_dir((FLAGS_data_root_path + FLAGS_mesh_path).c_str());
	assert(input_dir.exists());
	input_dir.setFilter(QDir::Files | QDir::Hidden | QDir::NoSymLinks);
	input_dir.setSorting(QDir::Name);


	// (Dissimilarity, mesh_filepath)
	std::vector< std::pair<Real, std::string> > label_matched_object(num_labels);
	for (unsigned int label_index = 0; label_index < num_labels; ++label_index)
		label_matched_object[label_index].first = std::numeric_limits<Real>::max();


	QFileInfoList dir_list = input_dir.entryInfoList();
	for (int i = 0; i < dir_list.size(); i++)
	{
		QFileInfo file_info = dir_list.at(i);
		if (file_info.exists() &&
			(file_info.suffix().compare("obj") == 0
			|| file_info.suffix().compare("off") == 0))
		{
			std::string mesh_filepath = std::string(file_info.filePath().toLocal8Bit());

			// Skip if the mesh is the same with the input mesh.
			if (mesh_filepath.compare(_mesh_filepath) == 0)
				continue;

			//QFileInfo file_info(mesh_filepath.c_str());
			std::string mesh_name(file_info.baseName().toLocal8Bit());
			std::string cuboid_filepath = FLAGS_training_dir + std::string("/") + mesh_name + std::string(".arff");

			bool ret = load_object_info(example_mesh, example_cuboid_structure,
				mesh_filepath.c_str(), LoadMesh, cuboid_filepath.c_str());
			if (!ret) continue;

			assert(example_cuboid_structure.num_labels() == num_labels);
			for (LabelIndex label_index = 0; label_index < num_labels; ++label_index)
			{
				assert(cuboid_structure_.label_cuboids_[label_index].size() <= 1);
				MeshCuboid *cuboid = NULL;
				if (!cuboid_structure_.label_cuboids_[label_index].empty())
					cuboid = cuboid_structure_.label_cuboids_[label_index].front();

				assert(example_cuboid_structure.label_cuboids_[label_index].size() <= 1);
				MeshCuboid *example_cuboid = NULL;
				if (!example_cuboid_structure.label_cuboids_[label_index].empty())
					example_cuboid = example_cuboid_structure.label_cuboids_[label_index].front();

				if (cuboid && example_cuboid)
				{
					// Measure similarity.
					MyMesh::Normal cuboid_size = cuboid->get_bbox_size();
					MyMesh::Normal example_cuboid_size = example_cuboid->get_bbox_size();
					MyMesh::Normal diff_size = (cuboid_size - example_cuboid_size);

					Real dissimilarity = 0.0;
					for (unsigned int i = 0; i < 3; ++i)
						dissimilarity += std::abs(diff_size[i]);

					if (dissimilarity < label_matched_object[label_index].first)
					{
						label_matched_object[label_index].first = dissimilarity;
						label_matched_object[label_index].second = mesh_filepath;
					}
				}
			}
		}
	}


	// NOTE:
	// Select the same 3D model for symmetric parts.
	for (std::vector< MeshCuboidSymmetryGroupInfo >::iterator jt = example_cuboid_structure.symmetry_group_info_.begin();
		jt != example_cuboid_structure.symmetry_group_info_.end(); ++jt)
	{
		MeshCuboidSymmetryGroupInfo &symmetry_group = (*jt);
		for (std::vector< std::pair<LabelIndex, LabelIndex> >::iterator kt = symmetry_group.pair_label_indices_.begin();
			kt != symmetry_group.pair_label_indices_.end(); ++kt)
		{
			LabelIndex label_index_1 = (*kt).first;
			LabelIndex label_index_2 = (*kt).second;

			if (label_matched_object[label_index_1].first <= label_matched_object[label_index_2].first)
				label_matched_object[label_index_2] = label_matched_object[label_index_1];
			else
				label_matched_object[label_index_1] = label_matched_object[label_index_2];
		}
	}


	if (_reconstructed_label_indices)
		cuboid_structure_.clear_label_sample_points(*_reconstructed_label_indices);
	else
		cuboid_structure_.clear_sample_points();

	for (LabelIndex label_index = 0; label_index < num_labels; ++label_index)
	{
		// Reconstruct designated labels only if '_reconstructed_label_indices' is provided.
		if (_reconstructed_label_indices)
		{
			bool exist = false;
			for (std::vector<LabelIndex>::const_iterator it = (*_reconstructed_label_indices).begin();
				it != (*_reconstructed_label_indices).end(); ++it)
			{
				if ((*it) == label_index)
				{
					exist = true;
					break;
				}
			}

			if (!exist) continue;
		}

		if (label_matched_object[label_index].first == std::numeric_limits<Real>::max())
			continue;

		std::string mesh_filepath = label_matched_object[label_index].second;
		QFileInfo file_info(mesh_filepath.c_str());
		std::string mesh_name(file_info.baseName().toLocal8Bit());
		std::string cuboid_filepath = FLAGS_training_dir + std::string("/") + mesh_name + std::string(".arff");

		bool ret = load_object_info(example_mesh, example_cuboid_structure,
			mesh_filepath.c_str(), LoadDenseSamplePoints, cuboid_filepath.c_str());
		assert(ret);

		assert(cuboid_structure_.label_cuboids_[label_index].size() <= 1);
		MeshCuboid *cuboid = NULL;
		if (!cuboid_structure_.label_cuboids_[label_index].empty())
			cuboid = cuboid_structure_.label_cuboids_[label_index].front();

		assert(example_cuboid_structure.label_cuboids_[label_index].size() <= 1);
		MeshCuboid *example_cuboid = NULL;
		if (!example_cuboid_structure.label_cuboids_[label_index].empty())
			example_cuboid = example_cuboid_structure.label_cuboids_[label_index].front();

		assert(cuboid);
		assert(example_cuboid);

		const int num_points = example_cuboid->num_sample_points();

		std::cout << "--------" << std::endl;
		std::cout << "Label (" << label_index << "):" << std::endl;
		std::cout << "Mesh name: " << mesh_name << std::endl;
		std::cout << "# of sample points: " << num_points << std::endl;
		std::cout << "Copying... ";

		for (int point_index = 0; point_index < num_points; ++point_index)
		{
			const MeshSamplePoint* sample_point = example_cuboid->get_sample_point(point_index);

			MyMesh::Point point = sample_point->point_;
			MyMesh::Normal normal = sample_point->normal_;

			//
			MyMesh::Point relative_point;
			MyMesh::Normal relative_normal;
			point = point - example_cuboid->get_bbox_center();

			for (unsigned int axis_index = 0; axis_index < 3; ++axis_index)
			{
				assert(example_cuboid->get_bbox_size()[axis_index] > 0);
				relative_point[axis_index] = dot(point, example_cuboid->get_bbox_axis(axis_index));
				relative_point[axis_index] *=
					(cuboid->get_bbox_size()[axis_index] / example_cuboid->get_bbox_size()[axis_index]);

				relative_normal[axis_index] = dot(normal, example_cuboid->get_bbox_axis(axis_index));
				relative_normal[axis_index] *=
					(cuboid->get_bbox_size()[axis_index] / example_cuboid->get_bbox_size()[axis_index]);
			}

			point = cuboid->get_bbox_center();
			normal = MyMesh::Normal(0.0);

			for (unsigned int axis_index = 0; axis_index < 3; ++axis_index)
			{
				point += (relative_point[axis_index] * cuboid->get_bbox_axis(axis_index));
				normal += (relative_normal[axis_index] * cuboid->get_bbox_axis(axis_index));
			}

			normal.normalize();
			//

			MeshSamplePoint *new_sample_point = cuboid_structure_.add_sample_point(point, normal);

			// Copy label confidence values.
			new_sample_point->label_index_confidence_ = sample_point->label_index_confidence_;

			cuboid->add_sample_point(new_sample_point);
		}

		std::cout << "Done." << std::endl;
	}
}

