#include "MeshCuboidStructure.h"

#include <deque>
#include <fstream>
#include <iostream>

MeshCuboidStructure::MeshCuboidStructure(const MyMesh* _mesh)
	: mesh_(_mesh)
	, query_label_index_(0)
	, translation_(0.0)
	, scale_(1.0)
{
	assert(_mesh);
}

MeshCuboidStructure::~MeshCuboidStructure()
{
}

void MeshCuboidStructure::clear()
{
	clear_labels();

	for (std::vector<MeshSamplePoint *>::iterator it = sample_points_.begin();
		it != sample_points_.end(); ++it)
		delete (*it);
	sample_points_.clear();

	translation_ = MyMesh::Normal(0.0);
	scale_ = 1.0;
}

void MeshCuboidStructure::clear_labels()
{
	labels_.clear();

	for (std::vector< std::vector<MeshCuboid *> >::iterator it = label_cuboids_.begin();
		it != label_cuboids_.end(); ++it)
	{
		for (std::vector<MeshCuboid *>::iterator jt = (*it).begin(); jt != (*it).end(); ++jt)
			delete (*jt);
	}
	label_cuboids_.clear();

	query_label_index_ = 0;
}

void MeshCuboidStructure::apply_mesh_transformation()
{
	assert(mesh_);
	reset_transformation();
	scale(mesh_->get_scale());
	translate(mesh_->get_translation());
}

void MeshCuboidStructure::translate(const MyMesh::Normal _translate)
{
	for (std::vector<MeshSamplePoint *>::iterator it = sample_points_.begin();
		it != sample_points_.end(); ++it)
	{
		MyMesh::Point &p = (*it)->point_;
		p = p + _translate;
	}

	translation_ += _translate;
}

void MeshCuboidStructure::scale(const Real _scale)
{
	assert(_scale > 0);
	for (std::vector<MeshSamplePoint *>::iterator it = sample_points_.begin();
		it != sample_points_.end(); ++it)
	{
		MyMesh::Point &p = (*it)->point_;
		p = p * _scale;
	}

	scale_ *= _scale;
	translation_ *= _scale;
}

void MeshCuboidStructure::reset_transformation()
{
	if (translation_ != MyMesh::Point(0.0) || scale_ != 1.0)
	{
		scale(1.0 / scale_);
		translate(-translation_);
	}

	assert(translation_ == MyMesh::Point(0.0));
	assert(scale_ == 1.0);
}

bool MeshCuboidStructure::load_sample_points(const char *_filename, bool _verbose)
{
	std::ifstream file(_filename);
	if (!file)
	{
		std::cerr << "Can't open file: \"" << _filename << "\"" << std::endl;
		return false;
	}

	if (_verbose)
		std::cout << "Loading " << _filename << "..." << std::endl;

	clear();

	std::string buffer;

	for (unsigned int count = 0; !file.eof(); count++)
	{
		std::getline(file, buffer);

		std::stringstream strstr(buffer);
		std::string token;
		std::getline(strstr, token, ' ');
		FaceIndex corr_fid = atoi(token.c_str());
		assert(corr_fid >= 0);

		if (strstr.eof())
			continue;

		Real bx, by, bz;
		std::getline(strstr, token, ' ');
		bx = atof(token.c_str());
		std::getline(strstr, token, ' ');
		by = atof(token.c_str());
		std::getline(strstr, token, ' ');
		bz = atof(token.c_str());
		MyMesh::Point bary_coord = MyMesh::Point(bx, by, bz);

		Real px, py, pz;
		std::getline(strstr, token, ' ');
		px = atof(token.c_str());
		std::getline(strstr, token, ' ');
		py = atof(token.c_str());
		std::getline(strstr, token, ' ');
		pz = atof(token.c_str());
		MyMesh::Point pos = MyMesh::Point(px, py, pz);

		sample_points_.push_back(
			new MeshSamplePoint(corr_fid, bary_coord, pos));
	}

	file.close();

	apply_mesh_transformation();

	std::cout << "Done." << std::endl;

	return true;
}

bool MeshCuboidStructure::load_sample_point_labels(const char *_filename, bool _verbose)
{
	std::ifstream file(_filename);
	if (!file)
	{
		std::cerr << "Can't open file: \"" << _filename << "\"" << std::endl;
		return false;
	}

	if (_verbose)
		std::cout << "Loading " << _filename << "..." << std::endl;

	clear_labels();

	std::string buffer;
	Label new_label = 0;
	bool loading_alpha = true;

	SamplePointIndex sample_point_index = 0;
	while (!file.eof() && sample_point_index < num_sample_points())
	{
		std::getline(file, buffer);

		std::stringstream strstr(buffer);
		std::string token;
		if (strstr.eof())
		{
			std::cerr << "Error: Wrong file format: \"" << _filename << "\"" << std::endl;
			return false;
		}
		else std::getline(strstr, token, ' ');

		if (token.compare("@ATTRIBUTE") == 0)
		{
			// NOTE:
			// In this file format, labels are defined by the recorded order.
			labels_.push_back(new_label);
			++new_label;
		}
		else if (token[0] == '@')
		{
			// Skip.
		}
		else
		{
			if (loading_alpha)
			{
				// NOTE:
				// Now adding labels is finished.
				// The last label is garbage.
				if (labels_.size() <= 1)
				{
					std::cerr << "Error: Wrong file format: \"" << _filename << "\"" << std::endl;
					return false;
				}
				labels_.pop_back();
				loading_alpha = false;
			}

			std::stringstream strstr(buffer);
			std::string token;

			MeshSamplePoint* sample_point = sample_points_[sample_point_index];
			sample_point->label_index_confidence_.resize(num_labels());
			std::fill(sample_point->label_index_confidence_.begin(), sample_point->label_index_confidence_.end(), 0);

			// NOTE:
			// In this file format, labels are defined by the recorded order.
			for (LabelIndex label_index = 0; label_index < labels_.size(); ++label_index)
			{
				if (strstr.eof())
				{
					std::cerr << "Error: Wrong file format: \"" << _filename << "\"" << std::endl;
					return false;
				}
				else std::getline(strstr, token, ',');
				sample_point->label_index_confidence_[label_index] = atof(token.c_str());
			}

			++sample_point_index;
		}
	}
	file.close();

	assert(sample_point_index == num_sample_points());

	// NOTE:
	// Draws all points.
	query_label_index_ = num_labels();


	std::cout << "Done." << std::endl;
	return true;
}

bool MeshCuboidStructure::load_cuboids(const char *_filename, bool _verbose)
{
	std::ifstream file(_filename);
	if (!file)
	{
		std::cerr << "Can't open file: \"" << _filename << "\"" << std::endl;
		return false;
	}

	if (_verbose)
		std::cout << "Loading " << _filename << "..." << std::endl;

	clear();

	std::string buffer;

	for (unsigned int label_index = 0; !file.eof(); label_index++)
	{
		std::getline(file, buffer);
		std::stringstream strstr(buffer);
		std::string token;

		if (buffer == "" || strstr.eof())
			continue;

		MeshCuboid *new_cuboid = new MeshCuboid(label_index);
		MyMesh::Point bbox_center(0.0);
		std::array<MyMesh::Point, MeshCuboid::k_num_corners> bbox_corners;

		for (unsigned int corner_index = 0; corner_index < MeshCuboid::k_num_corners; ++corner_index)
		{
			Real px, py, pz;
			assert(!strstr.eof());
			std::getline(strstr, token, ',');
			px = atof(token.c_str());
			assert(!strstr.eof());
			std::getline(strstr, token, ',');
			py = atof(token.c_str());
			assert(!strstr.eof());
			std::getline(strstr, token, ',');
			pz = atof(token.c_str());
			MyMesh::Point pos = MyMesh::Point(px, py, pz);

			bbox_center = bbox_center + pos;
			bbox_corners[corner_index] = pos;
		}

		bbox_center = bbox_center / MeshCuboid::k_num_corners;
		new_cuboid->set_bbox_center(bbox_center);
		new_cuboid->set_bbox_corners(bbox_corners);
		new_cuboid->cuboidize();

		// A label is the same with its label index.
		labels_.push_back(label_index);

		std::vector<MeshCuboid *> new_label_cuboid;
		new_label_cuboid.push_back(new_cuboid);
		label_cuboids_.push_back(new_label_cuboid);
	}

	////
	//for (unsigned int label_index = 0; label_index < num_labels(); ++label_index)
	if (num_labels() > 0)
	{
		unsigned int label_index = 0;
	
		MeshCuboid *new_cuboid = label_cuboids_[label_index][0];

		MyMesh::Point corner_0 = new_cuboid->get_bbox_corner(0);
		MyMesh::Point corner_1 = new_cuboid->get_bbox_corner(1);
		MyMesh::Point corner_2 = new_cuboid->get_bbox_corner(2);

		Real pz = corner_0[2];

		Real min_x = corner_0[0];
		Real max_x = corner_1[0];

		Real min_y = corner_0[1];
		Real max_y = corner_2[1];

		unsigned int num_axis_points = 30;
		for (unsigned int i = 0; i < num_axis_points; ++i)
		{
			Real px = (max_x - min_x) / static_cast<Real>(num_axis_points - 1) * i + min_x;
			for (unsigned int j = 0; j < num_axis_points; ++j)
			{
				Real py = (max_y - min_y) / static_cast<Real>(num_axis_points - 1) * j + min_y;

				MyMesh::Point pos = MyMesh::Point(px, py, pz);
				MeshSamplePoint *sample_point = new MeshSamplePoint(0, MyMesh::Point(0.0), pos);
				sample_point->label_index_confidence_.resize(num_labels());

				for (unsigned int label_index = 0; label_index < num_labels(); label_index++)
					sample_point->label_index_confidence_[label_index] = 0;

				sample_point->label_index_confidence_[0] = 1.0;

				sample_points_.push_back(sample_point);
				new_cuboid->add_sample_point(sample_point);
			}
		}

		//MyMesh::Point corner_0 = new_cuboid->get_bbox_corner(0);
		//MyMesh::Point corner_2 = new_cuboid->get_bbox_corner(2);
		//MyMesh::Point corner_4 = new_cuboid->get_bbox_corner(4);

		//Real px = corner_0[0];

		//Real min_y = corner_0[1];
		//Real max_y = corner_2[1];

		//Real min_z = corner_0[2];
		//Real max_z = corner_4[2];

		//unsigned int num_axis_points = 11;
		//for (unsigned int i = 0; i < num_axis_points; ++i)
		//{
		//	Real py = (max_y - min_y) / static_cast<Real>(num_axis_points - 1) * i + min_y;
		//	for (unsigned int j = 0; j < num_axis_points; ++j)
		//	{
		//		Real pz = (max_z - min_z) / static_cast<Real>(num_axis_points - 1) * j + min_z;

		//		MyMesh::Point pos = MyMesh::Point(px, py, pz);
		//		MeshSamplePoint *sample_point = new MeshSamplePoint(0, MyMesh::Point(0.0), pos);
		//		sample_point->label_index_confidence_.resize(num_labels());

		//		for (unsigned int label_index = 0; label_index < num_labels(); label_index++)
		//			sample_point->label_index_confidence_[label_index] = 0;

		//		sample_point->label_index_confidence_[0] = 1.0;

		//		sample_points_.push_back(sample_point);
		//		new_cuboid->add_point(sample_point);
		//	}
		//}
	}
	////

	file.close();

	// NOTE:
	// Draws all points.
	query_label_index_ = num_labels();

	std::cout << "Done." << std::endl;

	return true;
}

std::vector<MeshCuboid *> MeshCuboidStructure::get_all_cuboids() const
{
	std::vector<MeshCuboid *> all_cuboids;

	for (std::vector< std::vector<MeshCuboid *> >::const_iterator it = label_cuboids_.begin();
		it != label_cuboids_.end(); ++it)
		all_cuboids.insert(all_cuboids.end(), (*it).begin(), (*it).end());

	return all_cuboids;
}

void MeshCuboidStructure::make_mesh_vertices_as_sample_points()
{
	assert(mesh_);
	clear();

	// NOTE:
	// Since now the sample points are generated from mesh, 
	// the transformation of mesh is applied first.
	apply_mesh_transformation();

	sample_points_.reserve(3 * mesh_->n_faces());

	for (MyMesh::FaceIter f_it = mesh_->faces_begin(); f_it != mesh_->faces_end(); ++f_it)
	{
		MyMesh::FaceHandle fh = f_it.handle();
		FaceIndex corr_fid = fh.idx();

		unsigned int i = 0;
		for (MyMesh::ConstFaceVertexIter fv_it = mesh_->cfv_iter(fh); fv_it; ++fv_it, ++i)
		{
			assert(i < 3);
			MyMesh::VertexHandle vh = fv_it.handle();

			MyMesh::Point bary_coord(0.0);
			bary_coord[i] = 1.0;

			MyMesh::Point pos = mesh_->point(vh);

			sample_points_.push_back(
				new MeshSamplePoint(corr_fid, bary_coord, pos));
		}
	}
}

bool MeshCuboidStructure::apply_point_cuboid_label_map(
	const std::vector<PointCuboidLabelMap>& _point_cuboid_label_maps,
	const std::vector<Label>& _all_cuboid_labels)
{
	const unsigned int num_point_labels = _point_cuboid_label_maps.size();

	assert(label_cuboids_.size() == num_labels());
	if (num_labels() > num_point_labels)
	{
		std::cerr << "Error: # of current labels is greater than 'num_point_labels' value." << std::endl;
		return false;
	}


	// Update labels.
	std::vector<Label> old_labels = labels_;
	labels_ = _all_cuboid_labels;
	

	// Update sample points.
	for (std::vector<MeshSamplePoint *>::iterator point_it = sample_points_.begin(); point_it != sample_points_.end(); ++point_it)
	{
		std::vector<Real> &label_index_confidence = (*point_it)->label_index_confidence_;
		std::vector<Real> new_label_index_confidence(num_labels(), 0);

		assert(label_index_confidence.size() <= num_point_labels);
		for (LabelIndex point_label_index = 0; point_label_index < label_index_confidence.size(); ++point_label_index)
		{
			for (std::list<Label>::const_iterator label_it = _point_cuboid_label_maps[point_label_index].mapped_cuboid_labels_.begin();
				label_it != _point_cuboid_label_maps[point_label_index].mapped_cuboid_labels_.end(); ++label_it)
			{
				// For all mapped cuboid labels.
				Label new_label = (*label_it);
				LabelIndex new_label_index = get_label_index(new_label);
				assert(new_label_index < num_labels());

				// Copy the same confidence value.
				new_label_index_confidence[new_label_index] = label_index_confidence[point_label_index];
			}
		}

		label_index_confidence.swap(new_label_index_confidence);
	}


	// Update cuboids.
	if (query_label_index_ == old_labels.size())
		query_label_index_ = num_labels();

	std::vector< std::vector<MeshCuboid *> > new_label_cuboids(num_labels());
	for (LabelIndex point_label_index = 0; point_label_index < label_cuboids_.size(); ++point_label_index)
	{
		std::vector<MeshCuboid *>& label_cuboids = label_cuboids_[point_label_index];
		if (label_cuboids.empty())
			continue;

		assert(!_point_cuboid_label_maps[point_label_index].mapped_cuboid_labels_.empty());

		// NOTE:
		// Assign anyone of mapped cuboid label.
		Label new_label = _point_cuboid_label_maps[point_label_index].mapped_cuboid_labels_.front();
		LabelIndex new_label_index = get_label_index(new_label);
		assert(new_label_index < num_labels());

		if (!_point_cuboid_label_maps[point_label_index].is_multiple_cuboids_)
		{
			// NOTE:
			// If there should be only one cuboid for this label, merge all of the existing cuboids.
			// When considering a scene, the only-one cuboid condition cannot be used.
			MeshCuboid *merged_cuboid = MeshCuboid::merge_cuboids(point_label_index, label_cuboids);
			assert(merged_cuboid);

			for (std::vector<MeshCuboid *>::iterator cuboid_it = label_cuboids.begin();
				cuboid_it != label_cuboids.end(); ++cuboid_it)
				delete (*cuboid_it);

			label_cuboids.clear();
			label_cuboids.push_back(merged_cuboid);
		}

		new_label_cuboids[new_label_index].insert(new_label_cuboids[new_label_index].end(),
			label_cuboids.begin(), label_cuboids.end());

		// Update label index for each cuboid.
		for (std::vector<MeshCuboid *>::iterator cuboid_it = label_cuboids.begin(); cuboid_it != label_cuboids.end(); ++cuboid_it)
		{
			MeshCuboid *cuboid = (*cuboid_it);
			assert(cuboid->get_label_index() == point_label_index);
			cuboid->set_label_index(new_label_index);
		}

		// Update query label index.
		if (query_label_index_ == point_label_index)
			query_label_index_ = new_label_index;
	}

	label_cuboids_.swap(new_label_cuboids);


	std::cout << "Done." << std::endl;
	return true;
}

void MeshCuboidStructure::apply_test()
{
	std::list<Label> duplicated_labels;

	//
	duplicated_labels.push_back(3);
	duplicated_labels.push_back(4);
	duplicated_labels.push_back(5);
	//

	for (std::list<Label>::iterator label_it = duplicated_labels.begin(); label_it != duplicated_labels.end(); ++label_it)
	{
		Label label = (*label_it);
		LabelIndex label_index = get_label_index(label);
		assert(label_index < label_cuboids_.size());

		for (std::vector<MeshCuboid *>::iterator cuboid_it = label_cuboids_[label_index].begin();
			cuboid_it != label_cuboids_[label_index].end(); ++cuboid_it)
			delete (*cuboid_it);
		label_cuboids_[label_index].clear();
	}
}

void MeshCuboidStructure::apply_mesh_face_labels_to_sample_points()
{
	assert(mesh_);
	labels_.clear();

	std::list< std::vector<FaceIndex> > all_label_faces;
	mesh_->get_all_label_faces(all_label_faces);

	if (all_label_faces.empty())
		return;

	for (std::list< std::vector<FaceIndex> >::iterator l_it = all_label_faces.begin();
		l_it != all_label_faces.end(); ++l_it)
	{
		std::vector<FaceIndex> &label_faces = (*l_it);
		assert(!label_faces.empty());

		FaceIndex fid = label_faces.front();
		MyMesh::FaceHandle fh = mesh_->face_handle(fid);
		Label label = mesh_->property(mesh_->face_label_, fh);

		labels_.push_back(label);
	}

	if (labels_.empty())
		return;

	for (SamplePointIndex sample_point_index = 0; sample_point_index < num_sample_points(); ++sample_point_index)
	{
		MeshSamplePoint* sample_point = sample_points_[sample_point_index];

		FaceIndex fid = sample_point->corr_fid_;
		assert(fid < mesh_->n_faces());

		MyMesh::FaceHandle fh = mesh_->face_handle(fid);
		Label label = mesh_->property(mesh_->face_label_, fh);
		LabelIndex label_index = get_label_index(label);
		assert(label_index < num_labels());

		sample_point->label_index_confidence_.resize(num_labels());
		std::fill(sample_point->label_index_confidence_.begin(), sample_point->label_index_confidence_.end(), 0);

		// Note:
		// The confidence of the given label becomes '1.0'.
		sample_point->label_index_confidence_[label_index] = 1.0;
	}
}

void MeshCuboidStructure::apply_mesh_face_labels_to_cuboids()
{
	// Apple mesh face labels to sample points.
	apply_mesh_face_labels_to_sample_points();

	std::list<MeshCuboid *> part_list;

	for (std::vector< std::vector<MeshCuboid *> >::iterator it = label_cuboids_.begin();
		it != label_cuboids_.end(); ++it)
	{
		for (std::vector<MeshCuboid *>::iterator jt = (*it).begin(); jt != (*it).end(); ++jt)
		{
			MeshCuboid *label_cuboid = (*jt);
			// Update the label based on the label confidence values of sample points.
			label_cuboid->update_label_using_sample_points();
			part_list.push_back(label_cuboid);
		}
	}

	label_cuboids_.clear();
	label_cuboids_.resize(num_labels());

	for (std::list<MeshCuboid *>::iterator it = part_list.begin(); it != part_list.end(); ++it)
	{
		MeshCuboid *cuboid = (*it);
		LabelIndex label_index = cuboid->get_label_index();
		label_cuboids_[label_index].push_back(cuboid);
	}

	// NOTE:
	// Draws all boxes.
	query_label_index_ = num_labels();
}

void MeshCuboidStructure::get_mesh_face_label_cuboids()
{
	// Apple mesh face labels to sample points.
	apply_mesh_face_labels_to_sample_points();

	compute_label_cuboids();
}

void MeshCuboidStructure::compute_label_cuboids()
{
	// Initialize.
	for (std::vector< std::vector<MeshCuboid *> >::iterator it = label_cuboids_.begin();
		it != label_cuboids_.end(); ++it)
	{
		for (std::vector<MeshCuboid *>::iterator jt = (*it).begin(); jt != (*it).end(); ++jt)
			delete (*jt);
	}
	label_cuboids_.clear();
	label_cuboids_.resize(num_labels());

	for (LabelIndex label_index = 0; label_index < num_labels(); ++label_index)
	{
		label_cuboids_[label_index].clear();

		MeshCuboid *cuboid = new MeshCuboid(label_index);
		std::vector<MeshSamplePoint *> label_sample_points;

		for (SamplePointIndex sample_point_index = 0; sample_point_index < num_sample_points();
			++sample_point_index)
		{
			MeshSamplePoint *sample_point = sample_points_[sample_point_index];
			assert(sample_point);

			// Select sample points which has sufficient confidence for the given label.
			if (sample_point->label_index_confidence_[label_index] >= param_confidence_tol)
				label_sample_points.push_back(sample_point);
		}

		cuboid->add_sample_points(label_sample_points);
		bool ret = cuboid->compute_bbox();
		if (!ret)
		{
			delete cuboid;
		}
		else
		{
			label_cuboids_[label_index].push_back(cuboid);
		}
	}

	split_label_cuboids();

	// NOTE:
	// Draws all boxes.
	query_label_index_ = num_labels();
}

void MeshCuboidStructure::find_the_largest_label_cuboids()
{
	// Find the largest part for each part.

	assert(label_cuboids_.size() == num_labels());

	for (LabelIndex label_index = 0; label_index < num_labels(); ++label_index)
	{
		if (label_cuboids_[label_index].size() <= 1)
			continue;

		Real max_volume = -1.0;
		MeshCuboid *max_volume_label_cuboid = NULL;

		for (std::vector<MeshCuboid *>::iterator it = label_cuboids_[label_index].begin();
			it != label_cuboids_[label_index].end(); ++it)
		{
			MeshCuboid *label_cuboid = (*it);
			assert(label_cuboid->get_label_index() == label_index);
			if (label_cuboid->get_bbox_volume() > max_volume)
			{
				max_volume = label_cuboid->get_bbox_volume();
				max_volume_label_cuboid = label_cuboid;
			}
		}
		assert(max_volume_label_cuboid);

		std::vector<MeshCuboid *> new_label_cuboids;
		new_label_cuboids.push_back(max_volume_label_cuboid);
		label_cuboids_[label_index].swap(new_label_cuboids);

		// Delete label parts except the largest one.
		for (std::vector<MeshCuboid *>::iterator it = new_label_cuboids.begin();
			it != new_label_cuboids.end(); ++it)
			if (*it != max_volume_label_cuboid)
				delete (*it);
	}
}

bool MeshCuboidStructure::set_new_label_indices(const std::vector<Label>& _labels)
{
	// If any existing label does not exist in the new set of labels, return false.
	for (LabelIndex old_label_index = 0; old_label_index < num_labels(); ++old_label_index)
	{
		Label label = labels_[old_label_index];
		bool found = false;
		for (LabelIndex new_label_index = 0; new_label_index < _labels.size(); ++new_label_index)
		{
			if (label == _labels[new_label_index])
			{
				found = true;
				break;
			}
		}
		if (!found)
			return false;
	}

	// Update labels.
	std::vector<Label> old_labels = labels_;
	labels_ = _labels;

	// Update sample points.
	for (std::vector<MeshSamplePoint *>::iterator point_it = sample_points_.begin(); point_it != sample_points_.end(); ++point_it)
	{
		std::vector<Real> &label_index_confidence = (*point_it)->label_index_confidence_;
		std::vector<Real> new_label_index_confidence(num_labels(), 0);

		for (LabelIndex old_label_index = 0; old_label_index < label_index_confidence.size(); ++old_label_index)
		{
			Label label = old_labels[old_label_index];
			LabelIndex new_label_index = get_label_index(label);
			new_label_index_confidence[new_label_index] = label_index_confidence[old_label_index];
		}

		label_index_confidence.swap(new_label_index_confidence);
	}

	// Update label cuboids.
	std::vector< std::vector<MeshCuboid *> > new_label_cuboids(num_labels());
	for (LabelIndex old_label_index = 0; old_label_index < label_cuboids_.size(); ++old_label_index)
	{
		std::vector<MeshCuboid *>& label_cuboids = label_cuboids_[old_label_index];
		Label label = old_labels[old_label_index];
		LabelIndex new_label_index = get_label_index(label);
		assert(new_label_index < num_labels());

		new_label_cuboids[new_label_index].insert(new_label_cuboids[new_label_index].end(),
			label_cuboids.begin(), label_cuboids.end());

		for (std::vector<MeshCuboid *>::iterator cuboid_it = label_cuboids.begin(); cuboid_it != label_cuboids.end(); ++cuboid_it)
		{
			MeshCuboid *cuboid = (*cuboid_it);
			assert(cuboid->get_label_index() == old_label_index);
			cuboid->set_label_index(new_label_index);
		}

		// Update query label index.
		if (query_label_index_ == old_label_index)
			query_label_index_ = new_label_index;
	}

	label_cuboids_.swap(new_label_cuboids);
	
	return true;
}

std::vector<LabelIndex> MeshCuboidStructure::get_sample_point_label_indices()
{
	// Get sample point labels from the label confidence values.
	std::vector<LabelIndex> sample_point_label_indices(num_sample_points());

	for (SamplePointIndex sample_point_index = 0; sample_point_index < num_sample_points(); ++sample_point_index)
	{
		MeshSamplePoint* sample_point = sample_points_[sample_point_index];
		assert(sample_point);
		assert(sample_point->label_index_confidence_.size() == labels_.size());

		Real max_confidence = -std::numeric_limits<Real>::max();
		Label max_confidence_label_index = 0;

		for (LabelIndex label_index = 0; label_index < num_labels(); ++label_index)
		{
			Real confidence = sample_point->label_index_confidence_[label_index];

			if (confidence > max_confidence)
			{
				max_confidence = confidence;
				max_confidence_label_index = label_index;
			}
		}

		sample_point_label_indices[sample_point_index] = max_confidence_label_index;
	}

	return sample_point_label_indices;
}

void MeshCuboidStructure::print_label_cuboids(const LabelIndex _label_index)const
{
	assert(label_cuboids_.size() == num_labels());
	assert(_label_index < num_labels());
	
	Label label = get_label(_label_index);
	std::vector<MeshCuboid *> cuboid = label_cuboids_[_label_index];
	std::cout << "Label (" << label << ")" << std::endl;

	unsigned int count_cuboids = 0;
	for (std::vector<MeshCuboid *>::const_iterator it = cuboid.begin(); it != cuboid.end();
		++it, ++count_cuboids)
	{
		std::cout << "[" << count_cuboids << "]" << std::endl;
		(*it)->print_cuboid();
	}
}

Label MeshCuboidStructure::get_label(const LabelIndex _label_index)const
{
	assert(_label_index < num_labels());
	return labels_[_label_index];
}

bool MeshCuboidStructure::exist_label(const Label _label,
	LabelIndex* _label_index)const
{
	for (LabelIndex label_index = 0; label_index < num_labels(); ++label_index)
	{
		if (labels_[label_index] == _label)
		{
			// NOTE:
			// Assign label index if the pointer is provided.
			if (_label_index) (*_label_index) = label_index;
			return true;
		}
	}

	return false;
}

LabelIndex MeshCuboidStructure::get_label_index(const Label _label)const
{
	for (LabelIndex label_index = 0; label_index < num_labels(); ++label_index)
		if (labels_[label_index] == _label)
			return label_index;

	assert(false);
	return 0;
}

Label MeshCuboidStructure::get_new_label()const
{
	Label max_label = 0;
	for (std::vector<Label>::const_iterator it = labels_.begin(); it != labels_.end(); ++it)
	{
		max_label = std::max(max_label, *it);
	}

	return max_label + 1;
}

void MeshCuboidStructure::split_label_cuboids()
{
	assert(mesh_);
	assert(label_cuboids_.size() == num_labels());
	Real object_diameter = mesh_->get_object_diameter();

	for (LabelIndex label_index = 0; label_index < num_labels(); ++label_index)
	{
		std::vector<MeshCuboid *> &cuboids = label_cuboids_[label_index];
		if (cuboids.empty()) continue;

		std::vector<MeshCuboid *> new_cuboids;

		for (std::vector<MeshCuboid *>::iterator it = cuboids.begin(); it != cuboids.end(); ++it)
		{
			MeshCuboid *cuboid = (*it);
			std::vector<MeshCuboid *> sub_cuboids = cuboid->split_cuboid(object_diameter);
			new_cuboids.insert(new_cuboids.end(), sub_cuboids.begin(), sub_cuboids.end());

			// Note:
			// Delete existing parts.
			delete cuboid;
		}

		cuboids.swap(new_cuboids);

		print_label_cuboids(label_index);
	}
}

void MeshCuboidStructure::remove_occluded_sample_points(
	const std::set<FaceIndex>& _visible_face_indices)
{
	assert(mesh_);

	unsigned int num_faces = mesh_->n_faces();
	bool *is_face_visible = new bool[num_faces];
	memset(is_face_visible, false, num_faces * sizeof(bool));

	for (std::set<FaceIndex>::const_iterator it = _visible_face_indices.begin();
		it != _visible_face_indices.end(); ++it)
	{
		FaceIndex fid = (*it);
		assert(fid < num_faces);
		is_face_visible[fid] = true;
	}

	for (std::vector<MeshSamplePoint *>::iterator it = sample_points_.begin();
		it != sample_points_.end();)	// No increment.
	{
		if (!is_face_visible[(*it)->corr_fid_])
			it = sample_points_.erase(it);
		else
			++it;
	}

	delete[] is_face_visible;
}
