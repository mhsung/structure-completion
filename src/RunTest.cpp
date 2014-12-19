#include "MeshViewerWidget.h"
#include "MeshCuboidSolver.h"
#include "MeshCuboidNonLinearSolver.h"

#include <Eigen/Geometry>

void MeshViewerWidget::run_test_initialize()
{
	slotDrawMode(findAction(CUSTOM_VIEW));
	draw_cuboid_axes_ = true;

	cuboid_structure_.load_cuboids("cuboids.csv");

	unsigned int num_all_cuboid_labels = 2;

	joint_normal_relations_.clear();
	joint_normal_relations_.resize(num_all_cuboid_labels);
	for (LabelIndex label_index_1 = 0; label_index_1 < num_all_cuboid_labels; ++label_index_1)
	{
		joint_normal_relations_[label_index_1].resize(num_all_cuboid_labels);
		for (LabelIndex label_index_2 = 0; label_index_2 < num_all_cuboid_labels; ++label_index_2)
		{
			if (label_index_1 == label_index_2)
				continue;

			std::string relation_filename = "joint_normal_"
				+ std::to_string(label_index_1) + std::string("_")
				+ std::to_string(label_index_2) + std::string(".dat");
			bool ret = joint_normal_relations_[label_index_1][label_index_2].load_joint_normal_dat(
				relation_filename.c_str());
			if (!ret)
			{
				do {
					std::cout << '\n' << "Press the Enter key to continue.";
				} while (std::cin.get() != '\n');
			}
		}
	}

	delete test_joint_normal_predictor_;
	test_joint_normal_predictor_ = new MeshCuboidJointNormalRelationPredictor(joint_normal_relations_);

	all_cuboids_.clear();
	for (std::vector<std::vector<MeshCuboid *>>::iterator it = cuboid_structure_.label_cuboids_.begin();
		it != cuboid_structure_.label_cuboids_.end(); ++it)
		all_cuboids_.insert(all_cuboids_.end(), (*it).begin(), (*it).end());

	open_modelview_matrix_file("occlusion_pose.txt");
	set_view_direction();
	memcpy(test_occlusion_modelview_matrix_, modelview_matrix_, 16 * sizeof(double));
	open_modelview_matrix_file("pose.txt");

	update_cuboid_surface_points(cuboid_structure_, test_occlusion_modelview_matrix_);
	updateGL();
}

void MeshViewerWidget::run_test_translate(const MyMesh::Normal _translation)
{
	Eigen::Vector3d translation_vec;
	for (unsigned int i = 0; i < 3; ++i)
		translation_vec[i] = 0.001 * _translation[i];

	assert(!all_cuboids_.empty());
	all_cuboids_[0]->translate(translation_vec);

	update_cuboid_surface_points(cuboid_structure_, test_occlusion_modelview_matrix_);
	
	double single_total_energy = 0, pair_total_energy = 0;
	get_optimization_error(all_cuboids_, *test_joint_normal_predictor_, single_total_energy, pair_total_energy);
	std::cout << "Error: (pair = " << pair_total_energy
		<< ", single = " << single_total_energy << ")" << std::endl;

	updateGL();
}

void MeshViewerWidget::run_test_rotate(const Real _angle)
{
	Eigen::AngleAxisd rotation(_angle / 180.0 * M_PI, Eigen::Vector3d::UnitZ());

	assert(!all_cuboids_.empty());
	all_cuboids_[0]->rotate(rotation.toRotationMatrix(), false);
	all_cuboids_[0]->update_corner_points();

	Eigen::Matrix3d rotation_mat;
	for (unsigned int axis_index = 0; axis_index < 3; ++axis_index)
		for (unsigned int i = 0; i < 3; ++i)
			rotation_mat.row(axis_index)(i) = all_cuboids_[0]->get_bbox_axis(axis_index)[i];

	rotation.fromRotationMatrix(rotation_mat);
	std::cout << "Angle = " << rotation.angle() / M_PI * 180.0 << std::endl;

	update_cuboid_surface_points(cuboid_structure_, test_occlusion_modelview_matrix_);

	double single_total_energy = 0, pair_total_energy = 0;
	get_optimization_error(all_cuboids_, *test_joint_normal_predictor_, single_total_energy, pair_total_energy);
	std::cout << "Error: (pair = " << pair_total_energy
		<< ", single = " << single_total_energy << ")" << std::endl;

	updateGL();
}

void MeshViewerWidget::run_test_scale(const Real _scale_x, const Real _scale_y)
{
	assert(!all_cuboids_.empty());

	MyMesh::Point bbox_center = all_cuboids_[0]->get_bbox_center();
	MyMesh::Normal bbox_size = all_cuboids_[0]->get_bbox_size();

	bbox_center[0] += (0.001 * _scale_x);
	bbox_center[1] += (0.001 * _scale_y);

	bbox_size[0] += (2 * 0.001 * _scale_x);
	bbox_size[1] += (2 * 0.001 * _scale_y);

	all_cuboids_[0]->set_bbox_center(bbox_center);
	all_cuboids_[0]->set_bbox_size(bbox_size, true);

	update_cuboid_surface_points(cuboid_structure_, test_occlusion_modelview_matrix_);

	double single_total_energy = 0, pair_total_energy = 0;
	get_optimization_error(all_cuboids_, *test_joint_normal_predictor_, single_total_energy, pair_total_energy);
	std::cout << "Error: (pair = " << pair_total_energy
		<< ", single = " << single_total_energy << ")" << std::endl;

	updateGL();
}

void MeshViewerWidget::run_test_optimize()
{
	const double quadprog_ratio = 1E4;
	const unsigned int max_num_iterations = 30;

	optimize_attributes(cuboid_structure_, test_occlusion_modelview_matrix_,
		*test_joint_normal_predictor_, quadprog_ratio,
		"log.txt", max_num_iterations, this);
}