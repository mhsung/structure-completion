#include "MeshViewerCore.h"
#include "MeshCuboidEvaluator.h"
#include "MeshCuboidFusion.h"
#include "MeshCuboidParameters.h"
#include "simplerandom.h"

#include <sstream>
#include <QDir>
#include <QFileInfo>


MyMesh::Point get_transformed_point(const MyMesh::Point _input_point,
	MeshCuboid *_input_cuboid, MeshCuboid *_target_cuboid)
{
	assert(_input_cuboid);
	assert(_target_cuboid);

	MyMesh::Point local_coord;
	MyMesh::Point point = _input_point - _input_cuboid->get_bbox_center();
	for (unsigned int axis_index = 0; axis_index < 3; ++axis_index)
	{
		if (_input_cuboid->get_bbox_size()[axis_index] <= 0)
		{
			local_coord[axis_index] = 0;
		}
		else
		{
			local_coord[axis_index] = dot(point, _input_cuboid->get_bbox_axis(axis_index));
			local_coord[axis_index] *= (_target_cuboid->get_bbox_size()[axis_index]
				/ _input_cuboid->get_bbox_size()[axis_index]);
		}
	}

	MyMesh::Point target_point = _target_cuboid->get_bbox_center();
	for (unsigned int axis_index = 0; axis_index < 3; ++axis_index)
		target_point += (local_coord[axis_index] * _target_cuboid->get_bbox_axis(axis_index));

	return target_point;
}

MyMesh::Normal get_transformed_normal(const MyMesh::Normal _input_normal,
	MeshCuboid *_input_cuboid, MeshCuboid *_target_cuboid)
{
	assert(_input_cuboid);
	assert(_target_cuboid);

	MyMesh::Point local_coord;
	for (unsigned int axis_index = 0; axis_index < 3; ++axis_index)
	{
		if (_input_cuboid->get_bbox_size()[axis_index] <= 0)
		{
			local_coord[axis_index] = 0;
		}
		else
		{
			local_coord[axis_index] = dot(_input_normal, _input_cuboid->get_bbox_axis(axis_index));
			local_coord[axis_index] *= (_target_cuboid->get_bbox_size()[axis_index]
				/ _input_cuboid->get_bbox_size()[axis_index]);
		}
	}

	MyMesh::Normal target_normal = MyMesh::Normal(0.0);
	for (unsigned int axis_index = 0; axis_index < 3; ++axis_index)
		target_normal += (local_coord[axis_index] * _target_cuboid->get_bbox_axis(axis_index));

	target_normal.normalize();
	return target_normal;
}

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


	output_filename_sstr.clear(); output_filename_sstr.str("");
	output_filename_sstr << _output_file_prefix << std::string(".arff");
	cuboid_structure_.save_cuboids(output_filename_sstr.str());

	//
	ret = cuboid_structure_.load_dense_sample_points(dense_sample_filepath.c_str());
	assert(ret);

	set_modelview_matrix(_occlusion_modelview_matrix, false);
	remove_occluded_points();
	set_modelview_matrix(_snapshot_modelview_matrix);

	MeshCuboidStructure cuboid_structure_copy(cuboid_structure_);
	//


	updateGL();
	output_filename_sstr.clear(); output_filename_sstr.str("");
	output_filename_sstr << _output_file_prefix;
	snapshot(output_filename_sstr.str().c_str());


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
		_mesh_filepath, LoadGroundTruthCuboids);
	assert(ret);
	MeshCuboidEvaluator evaluator(&ground_truth_cuboid_structure);
	//


	// 0. Evaluate point labeling and cuboid distance.
	// NOTE: Assume that ground truth mesh face labels are already loaded.
	output_filename_sstr.clear(); output_filename_sstr.str("");
	output_filename_sstr << _output_file_prefix << std::string("_labeling_stats.csv");
	evaluator.evaluate_point_labeling(&cuboid_structure_, output_filename_sstr.str().c_str());

	// NOTE: Cuboid surface points will be replaced.
	output_filename_sstr.clear(); output_filename_sstr.str("");
	output_filename_sstr << _output_file_prefix << std::string("_cuboid_distance.csv");
	evaluator.evaluate_cuboid_distance(&cuboid_structure_, output_filename_sstr.str().c_str());


	// 1. Reconstruction using symmetry.
	cuboid_structure_.copy_sample_points_to_symmetric_position();

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
	
	// Extract symmetry information.
	output_filename_sstr.clear(); output_filename_sstr.str("");
	output_filename_sstr << _output_file_prefix << std::string("_symmetry_info.txt");
	cuboid_structure_.save_symmetry_groups(output_filename_sstr.str());

	setDrawMode(POINT_SAMPLES);
	updateGL();
	output_filename_sstr.clear(); output_filename_sstr.str("");
	output_filename_sstr << _output_file_prefix << std::string("_symmetry_output");
	snapshot(output_filename_sstr.str().c_str());

	setDrawMode(COLORED_POINT_SAMPLES);
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

	setDrawMode(POINT_SAMPLES);
	updateGL();
	output_filename_sstr.clear(); output_filename_sstr.str("");
	output_filename_sstr << _output_file_prefix << std::string("_database_output");
	snapshot(output_filename_sstr.str().c_str());

	setDrawMode(COLORED_POINT_SAMPLES);
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

	setDrawMode(POINT_SAMPLES);
	updateGL();
	output_filename_sstr.clear(); output_filename_sstr.str("");
	output_filename_sstr << _output_file_prefix << std::string("_fusion_output");
	snapshot(output_filename_sstr.str().c_str());

	setDrawMode(COLORED_POINT_SAMPLES);
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

void MeshViewerCore::run_render_output()
{
	std::string mesh_filepath = FLAGS_data_root_path + FLAGS_mesh_path + std::string("/") + FLAGS_mesh_filename;
	QFileInfo mesh_file(mesh_filepath.c_str());
	std::string mesh_name = std::string(mesh_file.baseName().toLocal8Bit());


	//
	MyMesh ground_truth_mesh;
	MeshCuboidStructure ground_truth_cuboid_structure(&ground_truth_mesh);

	bool ret = true;
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
		mesh_filepath.c_str(), LoadGroundTruthCuboids);
	assert(ret);
	MeshCuboidEvaluator evaluator(&ground_truth_cuboid_structure);
	//

	std::string point_filepath;
	std::string label_filepath;
	std::stringstream result_file_prefix;
	std::stringstream output_filename_sstr;

	unsigned int candidate_index = 0;
	while (true)
	{
		result_file_prefix.clear(); result_file_prefix.str("");
		result_file_prefix << FLAGS_output_dir << "/" << mesh_name << "/" << mesh_name
			<< "_" << candidate_index;

		// 0. Render resulting cuboids.
		point_filepath = FLAGS_output_dir + std::string("/") + mesh_name + std::string("/")
			+ mesh_name + std::string("_input.pts");
		std::string cuboid_filename = result_file_prefix.str() + std::string(".arff");

		ret = load_result_info(mesh_, cuboid_structure_,
			mesh_filepath.c_str(), point_filepath.c_str(), NULL, cuboid_filename.c_str());
		if (!ret) break;
		open_modelview_matrix_file(FLAGS_pose_filename.c_str());

		//
		cuboid_structure_.compute_symmetry_groups();
		//

		setDrawMode(CUSTOM_VIEW);
		draw_cuboid_axes_ = false;
		output_filename_sstr.clear(); output_filename_sstr.str("");
		output_filename_sstr << result_file_prefix.str();
		snapshot(output_filename_sstr.str().c_str());

		cuboid_structure_.clear_cuboids();
		setDrawMode(POINT_SAMPLES);


		// 1. Reconstruction using symmetry.
		point_filepath = result_file_prefix.str() + std::string("_symmetry.pts");
		label_filepath = result_file_prefix.str() + std::string("_symmetry_label.arff");

		ret = load_result_info(mesh_, cuboid_structure_,
			mesh_filepath.c_str(), point_filepath.c_str(), label_filepath.c_str(), NULL);
		if (!ret) break;
		open_modelview_matrix_file(FLAGS_pose_filename.c_str());

		output_filename_sstr.clear(); output_filename_sstr.str("");
		output_filename_sstr << result_file_prefix.str() << std::string("_symmetry_output");
		snapshot(output_filename_sstr.str().c_str());
		cuboid_structure_ = ground_truth_cuboid_structure;


		// 2. Reconstruction using database.
		point_filepath = result_file_prefix.str() + std::string("_database.pts");
		label_filepath = result_file_prefix.str() + std::string("_database_label.arff");

		ret = load_result_info(mesh_, cuboid_structure_,
			mesh_filepath.c_str(), point_filepath.c_str(), label_filepath.c_str(), NULL);
		if (!ret) break;
		open_modelview_matrix_file(FLAGS_pose_filename.c_str());

		output_filename_sstr.clear(); output_filename_sstr.str("");
		output_filename_sstr << result_file_prefix.str() << std::string("_database_output");
		snapshot(output_filename_sstr.str().c_str());
		cuboid_structure_ = ground_truth_cuboid_structure;


		// 3. Fusion.
		point_filepath = result_file_prefix.str() + std::string("_fusion.pts");
		label_filepath = result_file_prefix.str() + std::string("_fusion_label.arff");

		ret = load_result_info(mesh_, cuboid_structure_,
			mesh_filepath.c_str(), point_filepath.c_str(), label_filepath.c_str(), NULL);
		if (!ret) break;
		open_modelview_matrix_file(FLAGS_pose_filename.c_str());

		output_filename_sstr.clear(); output_filename_sstr.str("");
		output_filename_sstr << result_file_prefix.str() << std::string("_fusion_output");
		snapshot(output_filename_sstr.str().c_str());
		cuboid_structure_ = ground_truth_cuboid_structure;

		++candidate_index;
	}

	setDrawMode(POINT_SAMPLES);

	// 4. Part Assembly.
	result_file_prefix.clear(); result_file_prefix.str("");
	result_file_prefix << FLAGS_part_assembly_dir << "/" << mesh_name << std::string("/") << mesh_name;

	point_filepath = result_file_prefix.str() + std::string("_assembly.pts");
	label_filepath = result_file_prefix.str() + std::string("_assembly_label.arff");

	ret = load_result_info(mesh_, cuboid_structure_,
		mesh_filepath.c_str(), point_filepath.c_str(), label_filepath.c_str(), NULL);
	if (ret)
	{
		open_modelview_matrix_file(FLAGS_pose_filename.c_str());
		output_filename_sstr.clear(); output_filename_sstr.str("");
		output_filename_sstr << result_file_prefix.str() << std::string("_assembly_output");
		snapshot(output_filename_sstr.str().c_str());
		cuboid_structure_ = ground_truth_cuboid_structure;
	}


	// 5. Symmetry Detection.
	result_file_prefix.clear(); result_file_prefix.str("");
	result_file_prefix << FLAGS_symmetry_detection_dir << "/" << mesh_name << std::string("/") << mesh_name;

	point_filepath = result_file_prefix.str() + std::string("_symm_detection.pts");
	label_filepath = result_file_prefix.str() + std::string("_symm_detection_label.arff");

	ret = load_result_info(mesh_, cuboid_structure_,
		mesh_filepath.c_str(), point_filepath.c_str(), label_filepath.c_str(), NULL);
	if (ret)
	{
		open_modelview_matrix_file(FLAGS_pose_filename.c_str());
		output_filename_sstr.clear(); output_filename_sstr.str("");
		output_filename_sstr << result_file_prefix.str() << std::string("_symm_detection_output");
		snapshot(output_filename_sstr.str().c_str());
		cuboid_structure_ = ground_truth_cuboid_structure;
	}
}

void MeshViewerCore::run_render_evaluation()
{
	std::string mesh_filepath = FLAGS_data_root_path + FLAGS_mesh_path + std::string("/") + FLAGS_mesh_filename;
	QFileInfo mesh_file(mesh_filepath.c_str());
	std::string mesh_name = std::string(mesh_file.baseName().toLocal8Bit());
	

	//
	MyMesh ground_truth_mesh;
	MeshCuboidStructure ground_truth_cuboid_structure(&ground_truth_mesh);

	bool ret = true;
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
		mesh_filepath.c_str(), LoadGroundTruthCuboids);
	assert(ret);
	MeshCuboidEvaluator evaluator(&ground_truth_cuboid_structure);

	setDrawMode(COLORED_POINT_SAMPLES);
	//

	std::string point_filepath;
	std::string label_filepath;
	std::stringstream result_file_prefix;
	std::stringstream output_filename_sstr;

	unsigned int candidate_index = 0;
	while (true)
	{
		result_file_prefix.clear(); result_file_prefix.str("");
		result_file_prefix << FLAGS_output_dir << "/" << mesh_name << std::string("/") << mesh_name
			<< "_" << candidate_index;


		// 1. Reconstruction using symmetry.
		point_filepath = result_file_prefix.str() + std::string("_symmetry.pts");
		label_filepath = result_file_prefix.str() + std::string("_symmetry_label.arff");

		ret = load_result_info(mesh_, cuboid_structure_,
			mesh_filepath.c_str(), point_filepath.c_str(), label_filepath.c_str(), NULL);
		if (!ret) break;
		open_modelview_matrix_file(FLAGS_pose_filename.c_str());
		cuboid_structure_.compute_label_cuboids();

		output_filename_sstr.clear(); output_filename_sstr.str("");
		output_filename_sstr << result_file_prefix.str() << std::string("_symmetry");
		evaluator.evaluate_point_to_point_distances(&cuboid_structure_, output_filename_sstr.str().c_str());

		updateGL();
		output_filename_sstr.clear(); output_filename_sstr.str("");
		output_filename_sstr << result_file_prefix.str() << std::string("_symmetry_accuracy");
		snapshot(output_filename_sstr.str().c_str());

		cuboid_structure_ = ground_truth_cuboid_structure;
		updateGL();
		output_filename_sstr.clear(); output_filename_sstr.str("");
		output_filename_sstr << result_file_prefix.str() << std::string("_symmetry_completeness");
		snapshot(output_filename_sstr.str().c_str());


		// 2. Reconstruction using database.
		point_filepath = result_file_prefix.str() + std::string("_database.pts");
		label_filepath = result_file_prefix.str() + std::string("_database_label.arff");

		ret = load_result_info(mesh_, cuboid_structure_,
			mesh_filepath.c_str(), point_filepath.c_str(), label_filepath.c_str(), NULL);
		if (!ret) break;
		open_modelview_matrix_file(FLAGS_pose_filename.c_str());
		cuboid_structure_.compute_label_cuboids();

		output_filename_sstr.clear(); output_filename_sstr.str("");
		output_filename_sstr << result_file_prefix.str() << std::string("_database");
		evaluator.evaluate_point_to_point_distances(&cuboid_structure_, output_filename_sstr.str().c_str());

		updateGL();
		output_filename_sstr.clear(); output_filename_sstr.str("");
		output_filename_sstr << result_file_prefix.str() << std::string("_database_accuracy");
		snapshot(output_filename_sstr.str().c_str());

		cuboid_structure_ = ground_truth_cuboid_structure;
		updateGL();
		output_filename_sstr.clear(); output_filename_sstr.str("");
		output_filename_sstr << result_file_prefix.str() << std::string("_database_completeness");
		snapshot(output_filename_sstr.str().c_str());


		// 3. Fusion.
		point_filepath = result_file_prefix.str() + std::string("_fusion.pts");
		label_filepath = result_file_prefix.str() + std::string("_fusion_label.arff");

		ret = load_result_info(mesh_, cuboid_structure_,
			mesh_filepath.c_str(), point_filepath.c_str(), label_filepath.c_str(), NULL);
		if (!ret) break;
		open_modelview_matrix_file(FLAGS_pose_filename.c_str());
		cuboid_structure_.compute_label_cuboids();

		output_filename_sstr.clear(); output_filename_sstr.str("");
		output_filename_sstr << result_file_prefix.str() << std::string("_fusion");
		evaluator.evaluate_point_to_point_distances(&cuboid_structure_, output_filename_sstr.str().c_str());

		updateGL();
		output_filename_sstr.clear(); output_filename_sstr.str("");
		output_filename_sstr << result_file_prefix.str() << std::string("_fusion_accuracy");
		snapshot(output_filename_sstr.str().c_str());

		cuboid_structure_ = ground_truth_cuboid_structure;
		updateGL();
		output_filename_sstr.clear(); output_filename_sstr.str("");
		output_filename_sstr << result_file_prefix.str() << std::string("_fusion_completeness");
		snapshot(output_filename_sstr.str().c_str());

		++candidate_index;
	}
	

	// 4. Part Assembly.
	result_file_prefix.clear(); result_file_prefix.str("");
	result_file_prefix << FLAGS_part_assembly_dir << "/" << mesh_name << std::string("/") << mesh_name;

	point_filepath = result_file_prefix.str() + std::string("_assembly.pts");
	label_filepath = result_file_prefix.str() + std::string("_assembly_label.arff");

	ret = load_result_info(mesh_, cuboid_structure_,
		mesh_filepath.c_str(), point_filepath.c_str(), label_filepath.c_str(), NULL);
	if (ret)
	{
		open_modelview_matrix_file(FLAGS_pose_filename.c_str());
		cuboid_structure_.compute_label_cuboids();

		output_filename_sstr.clear(); output_filename_sstr.str("");
		output_filename_sstr << result_file_prefix.str() << std::string("_assembly");
		evaluator.evaluate_point_to_point_distances(&cuboid_structure_, output_filename_sstr.str().c_str());

		updateGL();
		output_filename_sstr.clear(); output_filename_sstr.str("");
		output_filename_sstr << result_file_prefix.str() << std::string("_assembly_accuracy");
		snapshot(output_filename_sstr.str().c_str());

		cuboid_structure_ = ground_truth_cuboid_structure;
		updateGL();
		output_filename_sstr.clear(); output_filename_sstr.str("");
		output_filename_sstr << result_file_prefix.str() << std::string("_assembly_completeness");
		snapshot(output_filename_sstr.str().c_str());
	}


	// 5. Symmetry Detection.
	result_file_prefix.clear(); result_file_prefix.str("");
	result_file_prefix << FLAGS_symmetry_detection_dir << "/" << mesh_name << std::string("/") << mesh_name;

	point_filepath = result_file_prefix.str() + std::string("_symm_detection.pts");
	label_filepath = result_file_prefix.str() + std::string("_symm_detection_label.arff");

	ret = load_result_info(mesh_, cuboid_structure_,
		mesh_filepath.c_str(), point_filepath.c_str(), label_filepath.c_str(), NULL);
	if (ret)
	{
		open_modelview_matrix_file(FLAGS_pose_filename.c_str());
		cuboid_structure_.compute_label_cuboids();

		output_filename_sstr.clear(); output_filename_sstr.str("");
		output_filename_sstr << result_file_prefix.str() << std::string("_symm_detection_recon");
		evaluator.evaluate_point_to_point_distances(&cuboid_structure_, output_filename_sstr.str().c_str());

		updateGL();
		output_filename_sstr.clear(); output_filename_sstr.str("");
		output_filename_sstr << result_file_prefix.str() << std::string("_symm_detection_accuracy");
		snapshot(output_filename_sstr.str().c_str());

		cuboid_structure_ = ground_truth_cuboid_structure;
		updateGL();
		output_filename_sstr.clear(); output_filename_sstr.str("");
		output_filename_sstr << result_file_prefix.str() << std::string("_symm_detection_completeness");
		snapshot(output_filename_sstr.str().c_str());
	}
}

void MeshViewerCore::run_extract_symmetry_info()
{
	std::string mesh_filepath = FLAGS_data_root_path + FLAGS_mesh_path + std::string("/") + FLAGS_mesh_filename;
	QFileInfo mesh_file(mesh_filepath.c_str());
	std::string mesh_name = std::string(mesh_file.baseName().toLocal8Bit());


	//
	bool ret = true;
	ret = ret & cuboid_structure_.load_labels((FLAGS_data_root_path +
		FLAGS_label_info_path + FLAGS_label_info_filename).c_str());
	ret = ret & cuboid_structure_.load_label_symmetries((FLAGS_data_root_path +
		FLAGS_label_info_path + FLAGS_label_symmetry_info_filename).c_str());

	// Load symmetry groups.
	ret = ret & cuboid_structure_.load_symmetry_groups((FLAGS_data_root_path +
		FLAGS_label_info_path + FLAGS_symmetry_group_info_filename).c_str());

	if (!ret)
	{
		do {
			std::cout << "Error: Cannot open label information files.";
			std::cout << '\n' << "Press the Enter key to continue.";
		} while (std::cin.get() != '\n');
	}
	//

	std::string cuboid_filepath;
	std::stringstream result_file_prefix;
	std::stringstream output_filename_sstr;

	setDrawMode(CUSTOM_VIEW);
	open_modelview_matrix_file(FLAGS_pose_filename.c_str());

	unsigned int candidate_index = 0;
	while (true)
	{
		result_file_prefix.clear(); result_file_prefix.str("");
		result_file_prefix << FLAGS_output_dir << "/" << mesh_name << std::string("/") << mesh_name
			<< "_" << candidate_index;

		cuboid_filepath = result_file_prefix.str() + std::string(".arff");

		ret = load_object_info(mesh_, cuboid_structure_, mesh_filepath.c_str(), LoadMesh, cuboid_filepath.c_str());
		if (!ret) break;

		//
		cuboid_structure_.compute_symmetry_groups();
		updateGL();
		//
		output_filename_sstr.clear(); output_filename_sstr.str("");
		output_filename_sstr << result_file_prefix.str() << std::string("_symmetry_info.txt");
		cuboid_structure_.save_symmetry_groups(output_filename_sstr.str());

		++candidate_index;
	}
}

void MeshViewerCore::reconstruct_scan(
	const char *_mesh_filepath,
	const GLdouble *_snapshot_modelview_matrix,
	const GLdouble *_occlusion_modelview_matrix,
	const char *_output_file_prefix)
{
	std::cout << "Reconstructing scanned data (No evaluation)." << std::endl;

	bool ret;
	std::stringstream output_filename_sstr;

	QFileInfo file_info(_mesh_filepath);
	std::string mesh_name(file_info.baseName().toLocal8Bit());
	std::string sample_filepath = FLAGS_data_root_path + FLAGS_sample_path
		+ std::string("/") + mesh_name + std::string(".pts");
	std::string dense_sample_filepath = FLAGS_data_root_path + FLAGS_dense_sample_path
		+ std::string("/") + mesh_name + std::string(".pts");


	output_filename_sstr.clear(); output_filename_sstr.str("");
	output_filename_sstr << _output_file_prefix << std::string(".arff");
	cuboid_structure_.save_cuboids(output_filename_sstr.str());

	//
	ret = cuboid_structure_.load_dense_sample_points(dense_sample_filepath.c_str());
	assert(ret);

	set_modelview_matrix(_occlusion_modelview_matrix, false);
	remove_occluded_points();
	set_modelview_matrix(_snapshot_modelview_matrix);

	MeshCuboidStructure cuboid_structure_copy(cuboid_structure_);
	//


	updateGL();
	output_filename_sstr.clear(); output_filename_sstr.str("");
	output_filename_sstr << _output_file_prefix;
	snapshot(output_filename_sstr.str().c_str());

	setDrawMode(POINT_SAMPLES);


	// 1. Reconstruction using symmetry.
	cuboid_structure_.copy_sample_points_to_symmetric_position();

	// NOTE:
	// The label of reconstructed points are recorded as confidence values.
	cuboid_structure_.set_sample_point_label_confidence_using_cuboids();

	output_filename_sstr.clear(); output_filename_sstr.str("");
	output_filename_sstr << _output_file_prefix << std::string("_symmetry");
	cuboid_structure_.save_sample_points_to_ply(output_filename_sstr.str().c_str());
	cuboid_structure_.save_sample_points((output_filename_sstr.str() + std::string(".pts")).c_str());
	cuboid_structure_.save_sample_point_labels((output_filename_sstr.str() + std::string("_label.arff")).c_str());
	MeshCuboidStructure symmetry_reconstruction(cuboid_structure_);

	updateGL();
	output_filename_sstr.clear(); output_filename_sstr.str("");
	output_filename_sstr << _output_file_prefix << std::string("_symmetry");
	snapshot(output_filename_sstr.str().c_str());


	// 2. Reconstruction using database.
	cuboid_structure_ = cuboid_structure_copy;

	// FIXME: Fix this hack!
	std::string temp_label_info_path = FLAGS_label_info_path;
	std::string temp_mesh_path = FLAGS_mesh_path;
	std::string temp_sample_path = FLAGS_sample_path;
	std::string temp_dense_sample_path = FLAGS_dense_sample_path;
	std::string temp_mesh_label_path = FLAGS_mesh_label_path;
	std::string temp_sample_label_path = FLAGS_sample_label_path;

	FLAGS_label_info_path = FLAGS_retrieval_label_info_path;
	FLAGS_mesh_path = FLAGS_retrieval_mesh_path;
	FLAGS_sample_path = FLAGS_retrieval_sample_path;
	FLAGS_dense_sample_path = FLAGS_retrieval_dense_sample_path;
	FLAGS_mesh_label_path = FLAGS_retrieval_mesh_label_path;
	FLAGS_sample_label_path = FLAGS_retrieval_sample_label_path;
	//

	reconstruct_database_prior(_mesh_filepath);

	//
	FLAGS_label_info_path = temp_label_info_path;
	FLAGS_mesh_path = temp_mesh_path;
	FLAGS_sample_path = temp_sample_path;
	FLAGS_dense_sample_path = temp_dense_sample_path;
	FLAGS_mesh_label_path = temp_mesh_label_path;
	FLAGS_sample_label_path = temp_sample_label_path;
	//

	// NOTE:
	// The label of reconstructed points are recorded as confidence values.
	cuboid_structure_.set_sample_point_label_confidence_using_cuboids();

	output_filename_sstr.clear(); output_filename_sstr.str("");
	output_filename_sstr << _output_file_prefix << std::string("_database");
	cuboid_structure_.save_sample_points_to_ply(output_filename_sstr.str().c_str());
	cuboid_structure_.save_sample_points((output_filename_sstr.str() + std::string(".pts")).c_str());
	cuboid_structure_.save_sample_point_labels((output_filename_sstr.str() + std::string("_label.arff")).c_str());
	MeshCuboidStructure database_reconstruction(cuboid_structure_);

	updateGL();
	output_filename_sstr.clear(); output_filename_sstr.str("");
	output_filename_sstr << _output_file_prefix << std::string("_database");
	snapshot(output_filename_sstr.str().c_str());


	// 3. Fusion.
	cuboid_structure_ = cuboid_structure_copy;
	//MeshViewerCore::reconstruct_fusion_simple(symmetry_reconstruction, database_reconstruction);

	MyMesh original_mesh;
	MeshCuboidStructure original_cuboid_structure(&original_mesh);
	//ret = load_object_info(original_mesh, original_cuboid_structure, _mesh_filepath, LoadSamplePoints);

	ret = original_mesh.open_mesh(_mesh_filepath, false);
	assert(ret);

	std::cout << " - Load sample points." << std::endl;
	ret = original_cuboid_structure.load_sample_points(sample_filepath.c_str(), false);
	assert(ret);

	// NOTE: Do not add outlier samples.
	reconstruct_fusion(_mesh_filepath,
		_snapshot_modelview_matrix, _occlusion_modelview_matrix,
		original_cuboid_structure, symmetry_reconstruction, database_reconstruction, cuboid_structure_, false);


	// NOTE:
	// The label of reconstructed points are recorded as confidence values.
	cuboid_structure_.set_sample_point_label_confidence_using_cuboids();

	output_filename_sstr.clear(); output_filename_sstr.str("");
	output_filename_sstr << _output_file_prefix << std::string("_fusion");
	cuboid_structure_.save_sample_points_to_ply(output_filename_sstr.str().c_str());
	cuboid_structure_.save_sample_points((output_filename_sstr.str() + std::string(".pts")).c_str());
	cuboid_structure_.save_sample_point_labels((output_filename_sstr.str() + std::string("_label.arff")).c_str());

	updateGL();
	output_filename_sstr.clear(); output_filename_sstr.str("");
	output_filename_sstr << _output_file_prefix << std::string("_fusion");
	snapshot(output_filename_sstr.str().c_str());

	cuboid_structure_ = cuboid_structure_copy;

	setDrawMode(CUSTOM_VIEW);
}

void MeshViewerCore::reconstruct_database_prior(
	const char *_mesh_filepath,
	const std::vector<LabelIndex> *_reconstructed_label_indices)
{
	const Real part_assembly_voxel_size = FLAGS_param_part_assembly_voxel_size *
		cuboid_structure_.mesh_->get_object_diameter();


	QFileInfo mesh_file_info(_mesh_filepath);
	std::string mesh_name(mesh_file_info.baseName().toLocal8Bit());

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


	// (Similarity, mesh_filepath)
	std::vector< std::pair<Real, std::string> > label_matched_objects(num_labels);
	for (LabelIndex label_index = 0; label_index < num_labels; ++label_index)
		label_matched_objects[label_index].first = -1;


	QFileInfoList dir_list = input_dir.entryInfoList();
	for (int i = 0; i < dir_list.size(); i++)
	{
		QFileInfo example_file_info = dir_list.at(i);
		if (example_file_info.exists() &&
			(example_file_info.suffix().compare("obj") == 0
			|| example_file_info.suffix().compare("off") == 0))
		{
			std::string example_mesh_filepath = std::string(example_file_info.filePath().toLocal8Bit());
			std::string example_mesh_name(example_file_info.baseName().toLocal8Bit());
			std::string cuboid_filepath = FLAGS_training_dir + std::string("/") + example_mesh_name + std::string(".arff");

			// Skip if the mesh is the same with the input mesh.
			if (example_mesh_name.compare(mesh_name) == 0)
				continue;

			//bool ret = load_object_info(example_mesh, example_cuboid_structure,
			//	example_mesh_filepath.c_str(), LoadMesh, cuboid_filepath.c_str());

			bool ret = load_object_info(example_mesh, example_cuboid_structure,
				example_mesh_filepath.c_str(), LoadDenseSamplePoints, cuboid_filepath.c_str(), false);
			if (!ret) continue;

			assert(example_cuboid_structure.num_labels() == num_labels);
			for (LabelIndex label_index = 0; label_index < num_labels; ++label_index)
			{
				assert(cuboid_structure_.label_cuboids_[label_index].size() <= 1);
				MeshCuboid *input_cuboid = NULL;
				if (!cuboid_structure_.label_cuboids_[label_index].empty())
					input_cuboid = cuboid_structure_.label_cuboids_[label_index].front();

				assert(example_cuboid_structure.label_cuboids_[label_index].size() <= 1);
				MeshCuboid *example_cuboid = NULL;
				if (!example_cuboid_structure.label_cuboids_[label_index].empty())
					example_cuboid = example_cuboid_structure.label_cuboids_[label_index].front();

				if (input_cuboid && example_cuboid)
				{
					// Measure similarity.
					std::vector<MyMesh::Point> input_cuboid_sample_points;
					input_cuboid->get_sample_points(input_cuboid_sample_points);
					unsigned int num_cuboid_sample_points = input_cuboid_sample_points.size();

					std::vector<MyMesh::Point> example_cuboid_sample_points;
					example_cuboid->get_sample_points(example_cuboid_sample_points);
					unsigned int num_example_cuboid_sample_points = example_cuboid_sample_points.size();


					Real similarity = 0.0;

					if (input_cuboid_sample_points.empty())
					{
						// If the input cuboid has no sample point, measure similairty using cuboid size.
						MyMesh::Normal input_cuboid_size = input_cuboid->get_bbox_size();
						MyMesh::Normal example_cuboid_size = example_cuboid->get_bbox_size();
						MyMesh::Normal diff_size = (input_cuboid_size - example_cuboid_size);

						similarity = 1.0;
						for (unsigned int i = 0; i < 3; ++i)
							similarity *= std::exp(-diff_size[i] * diff_size[i] / (0.1));
					}
					else
					{
						if ( example_cuboid_sample_points.empty())
							continue;

						// Fit example cuboid to the input cuboid.
						for (unsigned int point_index = 0; point_index < num_example_cuboid_sample_points; ++point_index)
						{
							MyMesh::Point point = example_cuboid_sample_points[point_index];
							MyMesh::Point transformed_point = get_transformed_point(point, example_cuboid, input_cuboid);
							example_cuboid_sample_points[point_index] = transformed_point;
						}

						MyMesh::Point local_coord_bbox_min = input_cuboid_sample_points.front();
						MyMesh::Point local_coord_bbox_max = input_cuboid_sample_points.front();
						for (unsigned int point_index = 0; point_index < num_cuboid_sample_points; ++point_index)
						{
							for (unsigned int axis_index = 0; axis_index < 3; ++axis_index)
							{
								local_coord_bbox_min[axis_index] = std::min(local_coord_bbox_min[axis_index],
									input_cuboid_sample_points[point_index][axis_index]);
								local_coord_bbox_max[axis_index] = std::max(local_coord_bbox_max[axis_index],
									input_cuboid_sample_points[point_index][axis_index]);
							}
						}

						// Ignore if the input cuboid is too small.
						bool is_too_small_cuboid = false;
						for (unsigned int axis_index = 0; axis_index < 3; ++axis_index)
						{
							if ((local_coord_bbox_max[axis_index] - local_coord_bbox_min[axis_index]) <= 0)
							{
								is_too_small_cuboid = true;
								break;
							}
						}

						if (is_too_small_cuboid) continue;


						MeshCuboidVoxelGrid local_coord_voxels(local_coord_bbox_min, local_coord_bbox_max,
							part_assembly_voxel_size);

						Eigen::VectorXd input_voxel_occupancies;
						local_coord_voxels.get_voxel_occupancies(input_cuboid_sample_points, input_voxel_occupancies);
						assert(input_voxel_occupancies.rows() == local_coord_voxels.n_voxels());

						Eigen::VectorXd example_cuboid_voxel_occupancies;
						local_coord_voxels.get_voxel_occupancies(example_cuboid_sample_points, example_cuboid_voxel_occupancies);
						assert(example_cuboid_voxel_occupancies.rows() == local_coord_voxels.n_voxels());

						int num_input_occupied_voxels = (input_voxel_occupancies.array() >= 1).count();
						int num_example_occupied_voxels = (example_cuboid_voxel_occupancies.array() >= 1).count();
						if (num_input_occupied_voxels == 0 || num_example_occupied_voxels == 0)
							continue;

						Real score = input_voxel_occupancies.dot(example_cuboid_voxel_occupancies);
						similarity = (1 - 0.7) * (score / num_example_occupied_voxels)
							+ (0.7) * (score / num_input_occupied_voxels);
					}


					assert(similarity >= 0);

					if (similarity > label_matched_objects[label_index].first)
					{
						label_matched_objects[label_index].first = similarity;
						label_matched_objects[label_index].second = example_mesh_filepath;
					}
				}
			}
		}
	}


	// NOTE:
	// Select the same 3D model for symmetric parts.
	for (std::vector< MeshCuboidSymmetryGroupInfo >::iterator it = example_cuboid_structure.symmetry_group_info_.begin();
		it != example_cuboid_structure.symmetry_group_info_.end(); ++it)
	{
		MeshCuboidSymmetryGroupInfo &symmetry_group = (*it);
		for (std::vector< std::pair<LabelIndex, LabelIndex> >::iterator jt = symmetry_group.pair_label_indices_.begin();
			jt != symmetry_group.pair_label_indices_.end(); ++jt)
		{
			LabelIndex label_index_1 = (*jt).first;
			LabelIndex label_index_2 = (*jt).second;

			// If either one of symmetric cuboids exists.
			if (label_matched_objects[label_index_1].first >= 0
				|| label_matched_objects[label_index_2].first >= 0)
			{
				if (label_matched_objects[label_index_1].first >= label_matched_objects[label_index_2].first)
					label_matched_objects[label_index_2] = label_matched_objects[label_index_1];
				else
					label_matched_objects[label_index_1] = label_matched_objects[label_index_2];
			}
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

		if (label_matched_objects[label_index].first < 0)
			continue;

		std::string mesh_filepath = label_matched_objects[label_index].second;
		QFileInfo file_info(mesh_filepath.c_str());
		std::string mesh_name(file_info.baseName().toLocal8Bit());
		std::string cuboid_filepath = FLAGS_training_dir + std::string("/") + mesh_name + std::string(".arff");

		std::cout << "--------" << std::endl;
		std::cout << "Label (" << label_index << "):" << std::endl;
		std::cout << "Mesh name: " << mesh_name << std::endl;

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
		std::cout << "# of sample points: " << num_points << std::endl;
		std::cout << "Copying... ";

		for (int point_index = 0; point_index < num_points; ++point_index)
		{
			const MeshSamplePoint* sample_point = example_cuboid->get_sample_point(point_index);

			MyMesh::Point point = sample_point->point_;
			MyMesh::Normal normal = sample_point->normal_;

			MyMesh::Point transformed_point = get_transformed_point(point, example_cuboid, cuboid);
			MyMesh::Normal transformed_normal = get_transformed_normal(normal, example_cuboid, cuboid);

			MeshSamplePoint *new_sample_point = cuboid_structure_.add_sample_point(
				transformed_point, transformed_normal);

			// Copy label confidence values.
			new_sample_point->label_index_confidence_ = sample_point->label_index_confidence_;

			cuboid->add_sample_point(new_sample_point);
		}

		std::cout << "Done." << std::endl;
	}
}

