#include "MeshCuboidStructure.h"

#include "MeshCuboidParameters.h"
#include "ICP.h"

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
	clear();
}

MeshCuboidStructure::MeshCuboidStructure(const MeshCuboidStructure& _other)
{
	deep_copy(_other);
}

MeshCuboidStructure& MeshCuboidStructure::operator=(const MeshCuboidStructure& _other)
{
	clear();
	deep_copy(_other);
	return (*this);
}

void MeshCuboidStructure::deep_copy(const MeshCuboidStructure& _other)
{
	// NOTE:
	// This instance should be cleaned prior to call this function.

	this->mesh_ = _other.mesh_;
	this->translation_ = _other.translation_;
	this->scale_ = _other.scale_;
	this->query_label_index_ = _other.query_label_index_;
	
	this->labels_ = _other.labels_;
	this->label_names_ = _other.label_names_;
	this->label_symmetries_ = _other.label_symmetries_;
	//this->label_children_ = _other.label_children_;
	this->symmetry_group_info_ = _other.symmetry_group_info_;


	// Deep copy sample points.
	assert(_other.sample_points_.size() == _other.num_sample_points());
	this->sample_points_.clear();
	this->sample_points_.reserve(_other.num_sample_points());

	for (std::vector<MeshSamplePoint *>::const_iterator it = _other.sample_points_.begin();
		it != _other.sample_points_.end(); ++it)
	{
		MeshSamplePoint *sample_point = new MeshSamplePoint(**it);
		this->sample_points_.push_back(sample_point);
	}

	// Deep copy label cuboids.
	assert(_other.label_cuboids_.size() == _other.num_labels());
	unsigned int num_labels = _other.num_labels();
	this->label_cuboids_.clear();
	this->label_cuboids_.resize(num_labels);

	for (LabelIndex label_index = 0; label_index < num_labels; ++label_index)
	{
		this->label_cuboids_[label_index].reserve(_other.label_cuboids_[label_index].size());
		for (std::vector<MeshCuboid *>::const_iterator it = _other.label_cuboids_[label_index].begin();
			it != _other.label_cuboids_[label_index].end(); ++it)
		{
			MeshCuboid *cuboid = new MeshCuboid(**it);

			const std::vector<MeshSamplePoint *> &cuboid_sample_points = cuboid->get_sample_points();
			std::vector<MeshSamplePoint *> new_cuboid_sample_points;

			for (std::vector<MeshSamplePoint *>::const_iterator it = cuboid_sample_points.begin();
				it != cuboid_sample_points.end(); ++it)
			{
				SamplePointIndex sample_point_index = (*it)->sample_point_index_;
				assert(sample_point_index < num_sample_points());
				assert(sample_points_[sample_point_index]->sample_point_index_ == sample_point_index);
				new_cuboid_sample_points.push_back(sample_points_[sample_point_index]);
			}

			// NOTE:
			// 'MeshCuboidStructure' class is a friend of 'MeshCuboid'.
			cuboid->sample_points_.swap(new_cuboid_sample_points);

			this->label_cuboids_[label_index].push_back(cuboid);
		}
	}

	// Deep copy symmetry groups.
	unsigned int num_symmetry_groups = _other.symmetry_groups_.size();
	this->symmetry_groups_.clear();
	this->symmetry_groups_.reserve(num_symmetry_groups);

	for (std::vector<MeshCuboidSymmetryGroup *>::const_iterator it = _other.symmetry_groups_.begin();
		it != _other.symmetry_groups_.end(); ++it)
	{
		assert(*it);
		MeshCuboidSymmetryGroup *symmetry_group = new MeshCuboidSymmetryGroup(**it);
		this->symmetry_groups_.push_back(symmetry_group);
	}
}

void MeshCuboidStructure::clear()
{
	clear_sample_points();
	clear_labels();
}

void MeshCuboidStructure::clear_sample_points()
{
	for (std::vector<MeshSamplePoint *>::iterator it = sample_points_.begin();
		it != sample_points_.end(); ++it)
		delete (*it);
	sample_points_.clear();

	//
	for (std::vector< std::vector<MeshCuboid *> >::iterator it = label_cuboids_.begin();
		it != label_cuboids_.end(); ++it)
	{
		for (std::vector<MeshCuboid *>::iterator jt = (*it).begin(); jt != (*it).end(); ++jt)
		{
			MeshCuboid* cuboid = (*jt);
			cuboid->clear_sample_points();
		}
	}
	//

	translation_ = MyMesh::Normal(0.0);
	scale_ = 1.0;
}

void MeshCuboidStructure::clear_label_sample_points(const std::vector<LabelIndex> &_label_indices)
{
	std::list<SamplePointIndex> deleted_sample_point_indices;

	// Collect sample points to be deleted.
	for (std::vector<LabelIndex>::const_iterator it = _label_indices.begin();
		it != _label_indices.end(); ++it)
	{
		LabelIndex label_index = (*it);
		for (std::vector<MeshCuboid *>::iterator jt = label_cuboids_[label_index].begin();
			jt != label_cuboids_[label_index].end(); ++jt)
		{
			MeshCuboid* cuboid = (*jt);

			const std::vector<MeshSamplePoint *> cuboid_sample_points = cuboid->get_sample_points();

			for (std::vector<MeshSamplePoint *>::const_iterator kt = cuboid_sample_points.begin();
				kt != cuboid_sample_points.end(); ++kt)
			{
				MeshSamplePoint* sample_point = (*kt);
				assert(sample_point);
				SamplePointIndex sample_point_index = sample_point->sample_point_index_;
				deleted_sample_point_indices.push_back(sample_point_index);
			}

			cuboid->clear_sample_points();
		}
	}

	// Delete sample points.
	for (std::list<SamplePointIndex>::iterator it = deleted_sample_point_indices.begin();
		it != deleted_sample_point_indices.end(); ++it)
	{
		SamplePointIndex sample_point_index = (*it);
		MeshSamplePoint* sample_point = this->sample_points_[sample_point_index];
		assert(sample_point);
		assert(sample_point->sample_point_index_ == sample_point_index);
		delete sample_point;
		this->sample_points_[sample_point_index] = NULL;
	}

	// Re-number sample points.
	SamplePointIndex new_sample_point_index = 0;
	for (std::vector<MeshSamplePoint *>::iterator it = sample_points_.begin();
		it != sample_points_.end(); )	// No increment.
	{
		if (!(*it))
		{
			it = sample_points_.erase(it);
		}
		else
		{
			(*it)->sample_point_index_ = new_sample_point_index;
			++it;
		}
	}
}

void MeshCuboidStructure::clear_cuboids()
{
	for (std::vector< std::vector<MeshCuboid *> >::iterator it = label_cuboids_.begin();
		it != label_cuboids_.end(); ++it)
	{
		for (std::vector<MeshCuboid *>::iterator jt = (*it).begin(); jt != (*it).end(); ++jt)
			delete (*jt);
	}

	label_cuboids_.clear();
	label_cuboids_.resize(num_labels());

	for (std::vector< MeshCuboidSymmetryGroup* >::iterator it = symmetry_groups_.begin();
		it != symmetry_groups_.end(); ++it)
		delete (*it);
	symmetry_groups_.clear();
}

void MeshCuboidStructure::clear_labels()
{
	clear_cuboids();

	labels_.clear();
	label_names_.clear();
	label_symmetries_.clear();
	//label_children_.clear();

	symmetry_group_info_.size();

	query_label_index_ = 0;
}

bool MeshCuboidStructure::load_cuboids(const std::string _filename, bool _verbose)
{
	if (labels_.empty())
	{
		std::cerr << "Error: Load label information first." << std::endl;
		return false;
	}

	std::ifstream file(_filename);
	if (!file)
	{
		std::cerr << "Can't open file: \"" << _filename << "\"" << std::endl;
		return false;
	}

	if (_verbose)
		std::cout << "Loading " << _filename << "..." << std::endl;


	clear_cuboids();

	std::string buffer;

	while (!file.eof())
	{
		std::getline(file, buffer);
		if (buffer == "") break;

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
			// Skip.
		}
		else if (token[0] == '@')
		{
			// Skip.
		}
		else
		{
			std::stringstream strstr(buffer);
			std::string token;
			bool is_failed = false;

			LabelIndex label_index;
			std::array<MyMesh::Normal, 3> bbox_axes;
			MyMesh::Point bbox_center;
			MyMesh::Normal bbox_size;

			if (std::getline(strstr, token, ',').fail()) is_failed = true;
			label_index = std::stoi(token);

			for (unsigned int axis_index = 0; axis_index < 3; ++axis_index)
			{
				for (unsigned int i = 0; i < 3; i++)
				{
					if (std::getline(strstr, token, ',').fail()) is_failed = true;
					bbox_axes[axis_index][i] = std::stof(token);
				}
			}

			for (unsigned int i = 0; i < 3; i++)
			{
				if (std::getline(strstr, token, ',').fail()) is_failed = true;
				bbox_center[i] = std::stof(token);
			}

			for (unsigned int i = 0; i < 3; i++)
			{
				if (std::getline(strstr, token, ',').fail()) is_failed = true;
				bbox_size[i] = std::stof(token);
			}

			if (is_failed)
			{
				std::cerr << "Error: Wrong file format: \"" << _filename << "\"" << std::endl;
				return false;
			}
			else if (label_index >= num_labels())
			{
				std::cerr << "Warning: The label index exceeds the number of labels"
					<< " (" << label_index << " >= " << num_labels() << ": \"" << _filename << "\"" << std::endl;
				continue;
			}

			MeshCuboid *cuboid = new MeshCuboid(label_index);
			cuboid->set_bbox_axes(bbox_axes);
			cuboid->set_bbox_center(bbox_center);
			cuboid->set_bbox_size(bbox_size);
			cuboid->update_corner_points();

			label_cuboids_[label_index].push_back(cuboid);
		}
	}
	file.close();


	// NOTE:
	// Draws all points.
	query_label_index_ = num_labels();


	std::cout << "Done." << std::endl;
	return true;
}

bool MeshCuboidStructure::save_cuboids(const std::string _filename, bool _verbose) const
{
	std::ofstream file(_filename);
	if (!file)
	{
		std::cerr << "Can't save file: \"" << _filename << "\"" << std::endl;
		return false;
	}

	if (_verbose)
		std::cout << "Saving " << _filename << "..." << std::endl;


	file << "@RELATION cuboid_structure" << std::endl;

	file << "@ATTRIBUTE label_index\tNUMERIC" << std::endl;

	file << "@ATTRIBUTE bbox_axis_0_x\tNUMERIC" << std::endl;
	file << "@ATTRIBUTE bbox_axis_0_y\tNUMERIC" << std::endl;
	file << "@ATTRIBUTE bbox_axis_0_z\tNUMERIC" << std::endl;

	file << "@ATTRIBUTE bbox_axis_1_x\tNUMERIC" << std::endl;
	file << "@ATTRIBUTE bbox_axis_1_y\tNUMERIC" << std::endl;
	file << "@ATTRIBUTE bbox_axis_1_z\tNUMERIC" << std::endl;

	file << "@ATTRIBUTE bbox_axis_2_x\tNUMERIC" << std::endl;
	file << "@ATTRIBUTE bbox_axis_2_y\tNUMERIC" << std::endl;
	file << "@ATTRIBUTE bbox_axis_2_z\tNUMERIC" << std::endl;

	file << "@ATTRIBUTE bbox_center_x\tNUMERIC" << std::endl;
	file << "@ATTRIBUTE bbox_center_y\tNUMERIC" << std::endl;
	file << "@ATTRIBUTE bbox_center_z\tNUMERIC" << std::endl;

	file << "@ATTRIBUTE bbox_size_x\tNUMERIC" << std::endl;
	file << "@ATTRIBUTE bbox_size_y\tNUMERIC" << std::endl;
	file << "@ATTRIBUTE bbox_size_z\tNUMERIC" << std::endl;

	file << "@DATA" << std::endl;
	

	const std::vector<MeshCuboid *> cuboids = get_all_cuboids();

	for (std::vector<MeshCuboid *>::const_iterator it = cuboids.begin();
		it != cuboids.end(); ++it)
	{
		const MeshCuboid *cuboid = (*it);
		assert(cuboid);
		bool first_item = true;

		if (first_item) first_item = false; else file << ",";
		file << cuboid->get_label_index();

		for (unsigned int axis_index = 0; axis_index < 3; ++axis_index)
		{
			for (unsigned int i = 0; i < 3; i++)
			{
				if (first_item) first_item = false; else file << ",";
				file << cuboid->get_bbox_axis(axis_index)[i];
			}
		}

		for (unsigned int i = 0; i < 3; i++)
		{
			if (first_item) first_item = false; else file << ",";
			file << cuboid->get_bbox_center()[i];
		}

		for (unsigned int i = 0; i < 3; i++)
		{
			if (first_item) first_item = false; else file << ",";
			file << cuboid->get_bbox_size()[i];
		}

		file << std::endl;
	}

	file.close();
	return true;
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

		scale_ = 1.0;
		translation_ = MyMesh::Point(0.0);
	}

	assert(translation_ == MyMesh::Point(0.0));
	assert(scale_ == 1.0);
}

bool MeshCuboidStructure::load_labels(const char *_filename, bool _verbose)
{
	std::ifstream file(_filename);
	if (!file)
	{
		std::cerr << "Can't open file: \"" << _filename << "\"" << std::endl;
		return false;
	}

	if (_verbose)
		std::cout << "Loading " << _filename << "..." << std::endl;


	// NOTE:
	// All cuboids are also deleted.
	clear_labels();


	std::string buffer;
	Label new_label = 0;

	while (!file.eof())
	{
		std::getline(file, buffer);
		if (buffer == "") break;

		std::stringstream sstr(buffer);
		
		const unsigned int num_tokens = 3;
		std::string tokens[num_tokens];

		for (unsigned int i = 0; i < num_tokens; ++i)
		{
			if (sstr.eof())
			{
				std::cerr << "Error: Wrong file format: \"" << _filename << "\"" << std::endl;
				return false;
			}

			std::getline(sstr, tokens[i], ' ');
		}

		if (tokens[1] != "pnts" || tokens[2] != "1")
		{
			std::cerr << "Error: Wrong file format: \"" << _filename << "\"" << std::endl;
			return false;
		}
		
		// NOTE:
		// In this file format, labels are defined by the recorded order.
		labels_.push_back(new_label);
		label_names_.push_back(tokens[0]);
		label_symmetries_.push_back(std::list<LabelIndex>());
		//label_children_.push_back(std::list<LabelIndex>());
		label_cuboids_.push_back(std::vector<MeshCuboid *>());
		++new_label;
	}

	file.close();


	// NOTE:
	// Draws all points.
	query_label_index_ = num_labels();


	std::cout << "Done." << std::endl;
	return true;
}

bool MeshCuboidStructure::load_label_symmetries(const char *_filename, bool _verbose)
{
	std::ifstream file(_filename);
	if (!file)
	{
		std::cerr << "Can't open file: \"" << _filename << "\"" << std::endl;
		return false;
	}

	if (_verbose)
		std::cout << "Loading " << _filename << "..." << std::endl;


	label_symmetries_.clear();
	label_symmetries_.resize(num_labels());


	std::string buffer;
	Label new_label = 0;

	while (!file.eof())
	{
		std::getline(file, buffer);
		if (buffer == "") break;

		std::stringstream sstr(buffer);
		std::list<LabelIndex> label_symmetry;

		while (!sstr.eof())
		{
			std::string token;
			std::getline(sstr, token, ' ');

			LabelIndex label_index = get_label_index(token);
			assert(label_index < num_labels());
			label_symmetry.push_back(label_index);
		}

		for (std::list<LabelIndex>::iterator it = label_symmetry.begin(); it != label_symmetry.end(); ++it)
		{
			LabelIndex label_index = (*it);

			for (std::list<LabelIndex>::iterator jt = label_symmetry.begin(); jt != label_symmetry.end(); ++jt)
			{
				LabelIndex n_label_index = (*jt);
				if (n_label_index != label_index)
				{
					label_symmetries_[label_index].push_back(n_label_index);
				}
			}
		}
	}

	file.close();


	std::cout << "Done." << std::endl;
	return true;
}

bool MeshCuboidStructure::load_symmetry_groups(const char *_filename, bool _verbose /*= true*/)
{
	std::ifstream file(_filename);
	if (!file)
	{
		std::cerr << "Can't open file: \"" << _filename << "\"" << std::endl;
		return false;
	}

	if (_verbose)
		std::cout << "Loading " << _filename << "..." << std::endl;


	for (std::vector< MeshCuboidSymmetryGroup* >::iterator it = symmetry_groups_.begin();
		it != symmetry_groups_.end(); ++it)
		delete (*it);
	symmetry_groups_.clear();
	symmetry_group_info_.clear();


	std::string buffer;
	MeshCuboidSymmetryGroupInfo new_symmetry_group;

	while (!file.eof())
	{
		std::getline(file, buffer);
		if (buffer == "") break;

		std::stringstream sstr(buffer);
		assert(!sstr.eof());
		std::string head_element;
		std::getline(sstr, head_element, ' ');

		std::vector<std::string> tokens;
		for (std::string each; std::getline(sstr, each, ' '); tokens.push_back(each));
		const unsigned int num_tokens = tokens.size();

		if (head_element == "symmetry_group")
		{
			//if (tokens.size() != 0)
			//{
			//	std::cerr << "Error: Wrong file format: \"" << _filename << "\"" << std::endl;
			//	return false;
			//}
			if (!new_symmetry_group.single_label_indices_.empty()
				|| !new_symmetry_group.pair_label_indices_.empty())
			{
				// Add the current symmetry group and create a new one.
				symmetry_group_info_.push_back(new_symmetry_group);
				new_symmetry_group = MeshCuboidSymmetryGroupInfo();
			}
		}
		else if (head_element == "axis_index")
		{
			if (tokens.size() != 1)
			{
				std::cerr << "Error: Wrong file format: \"" << _filename << "\"" << std::endl;
				return false;
			}
			else
			{
				unsigned int reflection_axis_index = atoi(tokens[0].c_str());
				assert(reflection_axis_index < 3);
				new_symmetry_group.reflection_axis_index_ = reflection_axis_index;
			}
		}
		else if (head_element == "single_label_indices")
		{
			for (unsigned int i = 0; i < num_tokens; ++i)
			{
				LabelIndex label_index = atoi(tokens[i].c_str());
				assert(label_index < labels_.size());
				new_symmetry_group.single_label_indices_.push_back(label_index);
			}
		}
		else if (head_element == "pair_label_indices")
		{
			if ((num_tokens % 2) != 0)
			{
				std::cerr << "Error: Wrong file format: \"" << _filename << "\"" << std::endl;
				return false;
			}
			else
			{
				for (unsigned int i = 0; i < num_tokens / 2; ++i)
				{
					LabelIndex label_index_1 = atoi(tokens[2 * i + 0].c_str());
					LabelIndex label_index_2 = atoi(tokens[2 * i + 1].c_str());
					assert(label_index_1 < labels_.size());
					assert(label_index_2 < labels_.size());
					new_symmetry_group.pair_label_indices_.push_back(
						std::make_pair(label_index_1, label_index_2));
				}
			}
		}
		else
		{
			std::cerr << "Error: Wrong file format: \"" << _filename << "\"" << std::endl;
			return false;
		}
	}

	if (!new_symmetry_group.single_label_indices_.empty()
		|| !new_symmetry_group.pair_label_indices_.empty())
	{
		symmetry_group_info_.push_back(new_symmetry_group);
		new_symmetry_group = MeshCuboidSymmetryGroupInfo();
	}

	file.close();

	std::cout << "Done." << std::endl;
	return true;
}

bool MeshCuboidStructure::load_sample_points(const char *_filename,
	bool _update_cuboid_memberships, bool _verbose)
{
	std::ifstream file(_filename);
	if (!file)
	{
		std::cerr << "Can't open file: \"" << _filename << "\"" << std::endl;
		return false;
	}

	if (_verbose)
		std::cout << "Loading " << _filename << "..." << std::endl;


	clear_sample_points();


	assert(mesh_);
	assert(mesh_->has_face_normals());

	std::string buffer;

	for (SamplePointIndex sample_point_index = 0; !file.eof(); ++sample_point_index)
	{
		std::getline(file, buffer);

		std::stringstream strstr(buffer);
		std::string token;
		std::getline(strstr, token, ' ');
		FaceIndex corr_fid = atoi(token.c_str());
		assert(corr_fid >= 0);
		assert(corr_fid < mesh_->n_faces());

		if (strstr.eof())
			continue;

		Real bx, by, bz;
		std::getline(strstr, token, ' ');
		bx = std::stof(token);
		std::getline(strstr, token, ' ');
		by = std::stof(token);
		std::getline(strstr, token, ' ');
		bz = std::stof(token);
		MyMesh::Point bary_coord = MyMesh::Point(bx, by, bz);

		Real px, py, pz;
		std::getline(strstr, token, ' ');
		px = std::stof(token);
		std::getline(strstr, token, ' ');
		py = std::stof(token);
		std::getline(strstr, token, ' ');
		pz = std::stof(token);
		MyMesh::Point point = MyMesh::Point(px, py, pz);
		
		MyMesh::Normal normal = mesh_->normal(mesh_->face_handle(corr_fid));

		MeshSamplePoint *sample_point = new MeshSamplePoint(sample_point_index, corr_fid, bary_coord, point, normal);
		sample_points_.push_back(sample_point);
		assert(sample_points_[sample_point_index] == sample_point);
	}

	file.close();

	apply_mesh_transformation();
	
	//
	if (_update_cuboid_memberships)
	{
		if (_verbose)
			std::cout << "Update cuboid memberships..." << std::endl;

		for (std::vector< std::vector<MeshCuboid *> >::iterator it = label_cuboids_.begin();
			it != label_cuboids_.end(); ++it)
		{
			for (std::vector<MeshCuboid *>::iterator jt = (*it).begin(); jt != (*it).end(); ++jt)
			{
				MeshCuboid* cuboid = (*jt);
				assert(cuboid->num_sample_points() == 0);

				for (std::vector<MeshSamplePoint *>::iterator kt = sample_points_.begin();
					kt != sample_points_.end(); ++kt)
				{
					MeshSamplePoint* sample_point = (*kt);
					if (cuboid->is_point_inside_cuboid(sample_point->point_))
					{
						cuboid->add_sample_point(sample_point);
					}
				}
			}
		}
	}
	//

	if (_verbose)
		std::cout << "Done." << std::endl;

	return true;
}

bool MeshCuboidStructure::load_dense_sample_points(const char *_filename, bool _verbose)
{
	// Compute segmentation of dense samples using segmented sparse samples.
	std::ifstream file(_filename);
	if (!file)
	{
		std::cerr << "Can't open file: \"" << _filename << "\"" << std::endl;
		return false;
	}

	if (_verbose)
		std::cout << "Loading " << _filename << "..." << std::endl;


	//
	std::vector< std::list<MeshCuboid *> > sparse_sample_to_cuboids(num_sample_points());

	for (std::vector< std::vector<MeshCuboid *> >::iterator it = label_cuboids_.begin();
		it != label_cuboids_.end(); ++it)
	{
		for (std::vector<MeshCuboid *>::iterator jt = (*it).begin(); jt != (*it).end(); ++jt)
		{
			MeshCuboid* cuboid = (*jt);
			const std::vector<MeshSamplePoint *> cuboid_sample_points = cuboid->get_sample_points();

			for (std::vector<MeshSamplePoint *>::const_iterator kt = cuboid_sample_points.begin();
				kt != cuboid_sample_points.end(); ++kt)
			{
				MeshSamplePoint* sample_point = (*kt);
				assert(sample_point);
				SamplePointIndex sample_point_index = sample_point->sample_point_index_;
				assert(sample_point_index < num_sample_points());
				sparse_sample_to_cuboids[sample_point_index].push_back(cuboid);
			}
		}
	}

	Eigen::MatrixXd sparse_sample_points(3, num_sample_points());

	// FIXME:
	// The type of indices should integer.
	// But, it causes compile errors in the 'ICP::get_closest_points' function.
	Eigen::MatrixXd sparse_sample_point_indices(1, num_sample_points());

	for (SamplePointIndex sample_point_index = 0; sample_point_index < num_sample_points();
		++sample_point_index)
	{
		assert(sample_points_[sample_point_index]);

		for (unsigned int i = 0; i < 3; ++i)
			sparse_sample_points.col(sample_point_index)(i) =
			sample_points_[sample_point_index]->point_[i];

		sparse_sample_point_indices.col(sample_point_index)(0) = 
			static_cast<double>(sample_point_index);
	}

	ANNpointArray sparse_sample_ann_points;
	ANNkd_tree *sparse_sample_ann_kd_tree = ICP::create_kd_tree(sparse_sample_points,
		sparse_sample_ann_points);
	assert(sparse_sample_ann_points);
	assert(sparse_sample_ann_kd_tree);
	//


	clear_sample_points();


	assert(mesh_);
	assert(mesh_->has_face_normals());

	std::string buffer;

	for (SamplePointIndex sample_point_index = 0; !file.eof(); ++sample_point_index)
	{
		std::getline(file, buffer);

		std::stringstream strstr(buffer);
		std::string token;
		std::getline(strstr, token, ' ');
		FaceIndex corr_fid = atoi(token.c_str());
		assert(corr_fid >= 0);
		assert(corr_fid < mesh_->n_faces());

		if (strstr.eof())
			continue;

		Real bx, by, bz;
		std::getline(strstr, token, ' ');
		bx = std::stof(token);
		std::getline(strstr, token, ' ');
		by = std::stof(token);
		std::getline(strstr, token, ' ');
		bz = std::stof(token);
		MyMesh::Point bary_coord = MyMesh::Point(bx, by, bz);

		Real px, py, pz;
		std::getline(strstr, token, ' ');
		px = std::stof(token);
		std::getline(strstr, token, ' ');
		py = std::stof(token);
		std::getline(strstr, token, ' ');
		pz = std::stof(token);
		MyMesh::Point point = MyMesh::Point(px, py, pz);

		MyMesh::Normal normal = mesh_->normal(mesh_->face_handle(corr_fid));

		MeshSamplePoint *sample_point = new MeshSamplePoint(sample_point_index, corr_fid, bary_coord, point, normal);
		sample_points_.push_back(sample_point);
		assert(sample_points_[sample_point_index] == sample_point);
	}

	file.close();

	apply_mesh_transformation();


	//
	Eigen::MatrixXd dense_sample_points(3, num_sample_points());
	Eigen::MatrixXd dense_sample_point_indices(1, num_sample_points());

	for (SamplePointIndex sample_point_index = 0; sample_point_index < num_sample_points();
		++sample_point_index)
	{
		assert(sample_points_[sample_point_index]);
		for (unsigned int i = 0; i < 3; ++i)
			dense_sample_points.col(sample_point_index)(i) =
			sample_points_[sample_point_index]->point_[i];
	}

	ICP::get_closest_points(sparse_sample_ann_kd_tree, dense_sample_points,
		sparse_sample_point_indices, dense_sample_point_indices);
	assert(dense_sample_point_indices.rows() == 1);
	assert(dense_sample_point_indices.cols() == num_sample_points());

	for (SamplePointIndex sample_point_index = 0; sample_point_index < num_sample_points();
		++sample_point_index)
	{
		SamplePointIndex sparse_sample_point_index =
			static_cast<SamplePointIndex>(dense_sample_point_indices.col(sample_point_index)(0));

		for (std::list<MeshCuboid *>::iterator it = sparse_sample_to_cuboids[sparse_sample_point_index].begin();
			it != sparse_sample_to_cuboids[sparse_sample_point_index].end(); ++it)
		{
			MeshCuboid *cuboid = (*it);
			cuboid->add_sample_point(sample_points_[sample_point_index]);
		}
	}


	annDeallocPts(sparse_sample_ann_points);
	delete sparse_sample_ann_kd_tree;
	//

	std::cout << "Done." << std::endl;

	return true;
}

bool MeshCuboidStructure::save_sample_points(const char *_filename, bool _verbose)
{
	std::ofstream file(_filename);
	if (!file)
	{
		std::cerr << "Can't save file: \"" << _filename << "\"" << std::endl;
		return false;
	}

	if (_verbose)
		std::cout << "Saving " << _filename << "..." << std::endl;

	// Note:
	// Save the original position of sample points.
	reset_transformation();

	for (SamplePointIndex sample_point_index = 0; sample_point_index < num_sample_points();
		++sample_point_index)
	{
		MeshSamplePoint *sample_point = sample_points_[sample_point_index];
		assert(sample_point);

		std::stringstream sstr;
		sstr << sample_point->corr_fid_ << " ";
		sstr << sample_point->bary_coord_[0] << " "
			<< sample_point->bary_coord_[1] << " "
			<< sample_point->bary_coord_[2] << " ";
		sstr << sample_point->point_[0] << " "
			<< sample_point->point_[1] << " "
			<< sample_point->point_[2] << " ";

		file << sstr.str() << std::endl;
	}

	file.close();

	apply_mesh_transformation();

	std::cout << "Done." << std::endl;

	return true;
}

bool MeshCuboidStructure::save_sample_points_to_ply(const char *_filename, bool _verbose)
{
	std::string ply_filename(_filename);
	ply_filename.append(".ply");
	std::ofstream file(ply_filename);
	if (!file)
	{
		std::cerr << "Can't save file: \"" << ply_filename << "\"" << std::endl;
		return false;
	}

	if (_verbose)
		std::cout << "Saving " << ply_filename << "..." << std::endl;

	// Note:
	// Save the original position of sample points.
	reset_transformation();

	// Print header.
	file << "ply" << std::endl;
	file << "format ascii 1.0" << std::endl;
	file << "element vertex " << num_sample_points() << std::endl;
	file << "property float x" << std::endl;
	file << "property float y" << std::endl;
	file << "property float z" << std::endl;
	file << "property float nx" << std::endl;
	file << "property float ny" << std::endl;
	file << "property float nz" << std::endl;
	file << "element face 0" << std::endl;
	file << "property list uchar int vertex_indices" << std::endl;
	file << "end_header" << std::endl;

	for (SamplePointIndex sample_point_index = 0; sample_point_index < num_sample_points();
		++sample_point_index)
	{
		MeshSamplePoint *sample_point = sample_points_[sample_point_index];
		assert(sample_point);

		std::stringstream sstr;
		sstr << sample_point->point_[0] << " "
			<< sample_point->point_[1] << " "
			<< sample_point->point_[2] << " "
			<< sample_point->normal_[0] << " "
			<< sample_point->normal_[1] << " "
			<< sample_point->normal_[2] << " ";

		file << sstr.str() << std::endl;
	}

	file.close();

	apply_mesh_transformation();

	std::cout << "Done." << std::endl;

	return true;
}

bool MeshCuboidStructure::load_sample_point_labels(const char *_filename, bool _verbose)
{
	if (labels_.empty())
	{
		std::cerr << "Error: Load label information first." << std::endl;
		return false;
	}
	else if (sample_points_.empty())
	{
		std::cerr << "Error: Load sample points first." << std::endl;
		return false;
	}

	std::ifstream file(_filename);
	if (!file)
	{
		std::cerr << "Can't open file: \"" << _filename << "\"" << std::endl;
		return false;
	}

	if (_verbose)
		std::cout << "Loading " << _filename << "..." << std::endl;


	std::string buffer;
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
			// Skip.
		}
		else if (token[0] == '@')
		{
			// Skip.
		}
		else
		{
			std::stringstream strstr(buffer);
			std::string token;

			MeshSamplePoint* sample_point = sample_points_[sample_point_index];
			assert(sample_point);
			sample_point->label_index_confidence_.clear();
			sample_point->label_index_confidence_.resize(num_labels());

			for (LabelIndex label_index = 0; label_index < num_labels(); ++label_index)
			{
				if (strstr.eof())
				{
					//std::cerr << "Error: Wrong file format: \"" << _filename << "\"" << std::endl;
					//return false;
					break;
				}
				else std::getline(strstr, token, ',');
				sample_point->label_index_confidence_[label_index] = std::stof(token);
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

bool MeshCuboidStructure::test_load_cuboids(const char *_filename, bool _verbose)
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
			px = std::stof(token);
			assert(!strstr.eof());
			std::getline(strstr, token, ',');
			py = std::stof(token);
			assert(!strstr.eof());
			std::getline(strstr, token, ',');
			pz = std::stof(token);
			MyMesh::Point pos = MyMesh::Point(px, py, pz);

			bbox_center = bbox_center + pos;
			bbox_corners[corner_index] = pos;
		}

		bbox_center = bbox_center / MeshCuboid::k_num_corners;
		new_cuboid->set_bbox_center(bbox_center);
		new_cuboid->set_bbox_corners(bbox_corners);
		new_cuboid->update_axes_center_size_corner_points();

		// A label is the same with its label index.
		labels_.push_back(label_index);

		std::vector<MeshCuboid *> new_label_cuboid;
		new_label_cuboid.push_back(new_cuboid);
		label_cuboids_.push_back(new_label_cuboid);
	}

	SamplePointIndex sample_point_index = 0;
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

				MyMesh::Point point = MyMesh::Point(px, py, pz);
				MeshSamplePoint *sample_point = new MeshSamplePoint(sample_point_index, 0, MyMesh::Point(0.0),
					point, MyMesh::Normal(0.0));
				sample_point->label_index_confidence_.resize(num_labels());
				++sample_point_index;

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

MeshSamplePoint *MeshCuboidStructure::add_sample_point(
	const MyMesh::Point& _point, const MyMesh::Normal& _normal)
{
	// NOTE:
	// Assume that the sample point index is the same with the index in the 'sample_points_' vector.
	SamplePointIndex new_sample_point_index = sample_points_.size();
	MeshSamplePoint *new_sample_point = new MeshSamplePoint(
		new_sample_point_index, 0, MyMesh::Point(0.0), _point, _normal);
	sample_points_.push_back(new_sample_point);
	return new_sample_point;
}

void MeshCuboidStructure::add_sample_points_from_mesh_vertices()
{
	assert(mesh_);
	assert(mesh_->has_vertex_normals());

	sample_points_.reserve(sample_points_.size() + 3 * mesh_->n_faces());
	SamplePointIndex sample_point_index = sample_points_.size();

	for (MyMesh::FaceIter f_it = mesh_->faces_begin(); f_it != mesh_->faces_end(); ++f_it)
	{
		MyMesh::FaceHandle fh = f_it.handle();
		FaceIndex corr_fid = fh.idx();
		Label label = mesh_->property(mesh_->face_label_, fh);

		LabelIndex label_index;
		bool ret = exist_label(label, &label_index);

		// NOTE:
		// If the mesh face label does not exist, ignore this mesh face.
		if (!ret) continue;
		assert(label_index < num_labels());

		unsigned int i = 0;
		for (MyMesh::ConstFaceVertexIter fv_it = mesh_->cfv_iter(fh); fv_it; ++fv_it, ++i)
		{
			assert(i < 3);
			MyMesh::VertexHandle vh = fv_it.handle();

			MyMesh::Point bary_coord(0.0);
			bary_coord[i] = 1.0;

			MyMesh::Point point = mesh_->point(vh);
			MyMesh::Normal normal = mesh_->normal(vh);

			MeshSamplePoint *sample_point = new MeshSamplePoint(sample_point_index, corr_fid, bary_coord, point, normal);

			sample_point->label_index_confidence_.clear();
			sample_point->label_index_confidence_.resize(num_labels(), 0.0);
			// Note:
			// The confidence of the given label becomes '1.0'.
			sample_point->label_index_confidence_[label_index] = 1.0;

			sample_points_.push_back(sample_point);
			++sample_point_index;
		}
	}

	assert(sample_point_index == num_sample_points());
}

void MeshCuboidStructure::remove_sample_points(const bool *is_sample_point_removed)
{
	assert(is_sample_point_removed);

	//
	for (std::vector< std::vector<MeshCuboid *> >::iterator it = label_cuboids_.begin();
		it != label_cuboids_.end(); ++it)
	{
		for (std::vector<MeshCuboid *>::iterator jt = (*it).begin(); jt != (*it).end(); ++jt)
		{
			MeshCuboid* cuboid = (*jt);
			cuboid->remove_sample_points(is_sample_point_removed);
		}
	}
	//

	SamplePointIndex new_sample_point_index = 0;
	for (std::vector<MeshSamplePoint *>::iterator it = sample_points_.begin();
		it != sample_points_.end();)	// No increment.
	{
		assert(*it);
		if (is_sample_point_removed[(*it)->sample_point_index_])
		{
			delete (*it);
			it = sample_points_.erase(it);
		}
		else
		{
			(*it)->sample_point_index_ = new_sample_point_index;
			++new_sample_point_index;
			++it;
		}
	}
}

void MeshCuboidStructure::apply_mesh_face_labels_to_sample_points()
{
	assert(mesh_);
	
	for (SamplePointIndex sample_point_index = 0; sample_point_index < num_sample_points(); ++sample_point_index)
	{
		MeshSamplePoint* sample_point = sample_points_[sample_point_index];

		FaceIndex fid = sample_point->corr_fid_;
		assert(fid < mesh_->n_faces());

		MyMesh::FaceHandle fh = mesh_->face_handle(fid);
		Label label = mesh_->property(mesh_->face_label_, fh);

		sample_point->label_index_confidence_.clear();
		sample_point->label_index_confidence_.resize(num_labels(), 0.0);

		LabelIndex label_index;
		bool ret = exist_label(label, &label_index);

		// NOTE:
		// If the mesh face label does not exist, ignore this sample point.
		if (ret)
		{
			assert(label_index < num_labels());

			// Note:
			// The confidence of the given label becomes '1.0'.
			sample_point->label_index_confidence_[label_index] = 1.0;
		}
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

void MeshCuboidStructure::get_mesh_face_label_cuboids(bool _add_sample_points_from_mesh_vertices)
{
	// Apple mesh face labels to sample points.
	apply_mesh_face_labels_to_sample_points();

	compute_label_cuboids();
}

void MeshCuboidStructure::compute_label_cuboids()
{
	clear_cuboids();

	for (LabelIndex label_index = 0; label_index < num_labels(); ++label_index)
	{
		assert(label_cuboids_[label_index].empty());
		MeshCuboid *cuboid = new MeshCuboid(label_index);
		std::vector<MeshSamplePoint *> label_sample_points;

		for (SamplePointIndex sample_point_index = 0; sample_point_index < num_sample_points();
			++sample_point_index)
		{
			MeshSamplePoint *sample_point = sample_points_[sample_point_index];
			assert(sample_point);

			// Select sample points which has sufficient confidence for the given label.
			if (sample_point->label_index_confidence_[label_index] >= FLAGS_param_min_sample_point_confidence)
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

LabelIndex MeshCuboidStructure::get_label_index(const std::string _label_name) const
{
	assert(label_names_.size() == num_labels());

	for (LabelIndex label_index = 0; label_index < num_labels(); ++label_index)
		if (label_names_[label_index] == _label_name)
			return label_index;

	assert(false);
	return 0;
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
	}
}

void MeshCuboidStructure::remove_symmetric_cuboids()
{
	// Remove cuboids in symmetric labels (when the same cuboids are duplicated for symmetric labels).

	assert(label_symmetries_.size() == num_labels());
	assert(label_cuboids_.size() == num_labels());

	bool *is_label_visited = new bool[num_labels()];
	memset(is_label_visited, false, num_labels()*sizeof(bool));

	for (LabelIndex label_index = 0; label_index < num_labels(); ++label_index)
	{
		if (is_label_visited[label_index])
		{
			for (std::vector<MeshCuboid *>::iterator jt = label_cuboids_[label_index].begin();
				jt != label_cuboids_[label_index].end(); ++jt)
				delete (*jt);
			label_cuboids_[label_index].clear();
			continue;
		}

		is_label_visited[label_index] = true;
		for (std::list<LabelIndex>::iterator it = label_symmetries_[label_index].begin();
			it != label_symmetries_[label_index].end(); ++it)
			is_label_visited[*it] = true;
	}

	delete[] is_label_visited;
}

void MeshCuboidStructure::compute_symmetry_groups()
{
	for (std::vector< MeshCuboidSymmetryGroup* >::iterator it = symmetry_groups_.begin();
		it != symmetry_groups_.end(); ++it)
		delete (*it);
	symmetry_groups_.clear();

	const std::vector<MeshCuboid*> cuboids = get_all_cuboids();

	for (std::vector< MeshCuboidSymmetryGroupInfo >::iterator it = symmetry_group_info_.begin();
		it != symmetry_group_info_.end(); ++it)
	{
		MeshCuboidSymmetryGroup* group = MeshCuboidSymmetryGroup::constructor(*it, cuboids);
		if (group) symmetry_groups_.push_back(group);
	}
}

void MeshCuboidStructure::copy_sample_points_to_symmetric_position()
{
	for (std::vector< MeshCuboidSymmetryGroup* >::const_iterator it = symmetry_groups_.begin();
		it != symmetry_groups_.end(); ++it)
	{
		copy_sample_points_to_symmetric_position(*it);
	}
}

void MeshCuboidStructure::copy_sample_points_to_symmetric_position(
	const MeshCuboidSymmetryGroup* _symmetry_group)
{
	assert(_symmetry_group);
	std::vector<MeshCuboid *> all_cuboids = get_all_cuboids();
	const unsigned int num_cuboids = all_cuboids.size();

	std::vector<unsigned int> single_cuboid_indices;
	_symmetry_group->get_single_cuboid_indices(all_cuboids, single_cuboid_indices);

	for (std::vector<unsigned int>::const_iterator it = single_cuboid_indices.begin();
		it != single_cuboid_indices.end(); ++it)
	{
		const unsigned int cuboid_index = (*it);
		copy_sample_points_to_symmetric_position(
			_symmetry_group, all_cuboids[cuboid_index], all_cuboids[cuboid_index]);
	}


	std::vector< std::pair<unsigned int, unsigned int> > pair_cuboid_indices;
	_symmetry_group->get_pair_cuboid_indices(all_cuboids, pair_cuboid_indices);

	for (std::vector< std::pair<unsigned int, unsigned int> >::const_iterator it = pair_cuboid_indices.begin();
		it != pair_cuboid_indices.end(); ++it)
	{
		const unsigned int cuboid_index_1 = (*it).first;
		const unsigned int cuboid_index_2 = (*it).second;

		// 1 -> 2.
		copy_sample_points_to_symmetric_position(
			_symmetry_group, all_cuboids[cuboid_index_1], all_cuboids[cuboid_index_2]);

		// 2 -> 1.
		copy_sample_points_to_symmetric_position(
			_symmetry_group, all_cuboids[cuboid_index_2], all_cuboids[cuboid_index_1]);
	}
}

void MeshCuboidStructure::copy_sample_points_to_symmetric_position(
	const MeshCuboidSymmetryGroup* _symmetry_group,
	const MeshCuboid *_cuboid_1, MeshCuboid *_cuboid_2)
{
	assert(_symmetry_group);
	if (!_cuboid_1 || !_cuboid_2)
		return;

	/*
	//
	Eigen::MatrixXd symmetric_sample_points_1(3, _cuboid_1->num_sample_points());
	Eigen::MatrixXd sample_points_2(3, _cuboid_2->num_sample_points());

	for (SamplePointIndex sample_point_index = 0; sample_point_index < _cuboid_1->num_sample_points();
		++sample_point_index)
	{
		assert(_cuboid_1->get_sample_point(sample_point_index));
		MyMesh::Point point_1 = _cuboid_1->get_sample_point(sample_point_index)->point_;
		//
		MyMesh::Point symmetric_point_1 = _symmetry_group->get_symmetric_point(point_1);
		//
		for (unsigned int i = 0; i < 3; ++i)
			symmetric_sample_points_1.col(sample_point_index)(i) = symmetric_point_1[i];
	}

	for (SamplePointIndex sample_point_index = 0; sample_point_index < _cuboid_2->num_sample_points();
		++sample_point_index)
	{
		assert(_cuboid_2->get_sample_point(sample_point_index));
		MyMesh::Point point_2 = _cuboid_2->get_sample_point(sample_point_index)->point_;
		for (unsigned int i = 0; i < 3; ++i)
			sample_points_2.col(sample_point_index)(i) = point_2[i];
	}

	const Real neighbor_distance = 10 * FLAGS_param_sample_point_neighbor_distance * mesh_->get_object_diameter();

	Eigen::Matrix3d rotation_mat;
	Eigen::Vector3d translation_vec;
	double icp_error = ICP::run_iterative_closest_points(symmetric_sample_points_1, sample_points_2,
		rotation_mat, translation_vec, &neighbor_distance);

	printf("%d -> %d: %lf\n", _cuboid_1->label_index_, _cuboid_2->label_index_, icp_error);
	//
	*/


	const int num_points_1 = _cuboid_1->num_sample_points();

	for (int point_index_1 = 0; point_index_1 < num_points_1; ++point_index_1)
	{
		const MeshSamplePoint* sample_point_1 = _cuboid_1->get_sample_point(point_index_1);
		assert(sample_point_1);
		MyMesh::Point point_1 = sample_point_1->point_;
		MyMesh::Normal normal_1 = sample_point_1->normal_;

		MyMesh::Point symmetric_point_1;
		//if (icp_error >= 0)
		//{
		//	// If ICP succeeded, copy the aligned point.
		//	for (unsigned int i = 0; i < 3; ++i)
		//		symmetric_point_1[i] = symmetric_sample_points_1.col(point_index_1)[i];
		//}
		//else
		{
			symmetric_point_1 = _symmetry_group->get_symmetric_point(point_1);
		}

		// FIXME:
		// IF ICP succeeded, how to compute the symmetric normal direction?
		MyMesh::Normal symmetric_normal_1 = _symmetry_group->get_symmetric_normal(normal_1);

		MeshSamplePoint *symmetric_sample_point = add_sample_point(symmetric_point_1, symmetric_normal_1);

		// Copy label confidence values.
		symmetric_sample_point->label_index_confidence_ = sample_point_1->label_index_confidence_;

		_cuboid_2->add_sample_point(symmetric_sample_point);
	}
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

/*
bool MeshCuboidStructure::is_label_group(LabelIndex _label_index)
{
	assert(_label_index < label_children_.size());
	bool ret = (!label_children_[_label_index].empty());
	return ret;
}

void MeshCuboidStructure::clear_symmetric_group_labels()
{
	std::vector<Label> new_labels;
	std::vector<std::string> new_label_names;
	std::vector< std::list<LabelIndex> > new_label_symmetries;
	std::vector< std::vector<MeshCuboid *> > new_label_cuboids;
	std::vector< std::list<LabelIndex> > new_label_children;

	new_labels.reserve(num_labels());
	new_label_names.reserve(num_labels());
	new_label_symmetries.reserve(num_labels());
	new_label_cuboids.reserve(num_labels());
	new_label_children.reserve(num_labels());

	for (unsigned int label_index = 0; label_index < num_labels(); ++label_index)
	{
		if (!is_label_group(label_index))
		{
			new_labels.push_back(labels_[label_index]);
			new_label_names.push_back(label_names_[label_index]);
			new_label_symmetries.push_back(label_symmetries_[label_index]);
			new_label_children.push_back(label_children_[label_index]);
			new_label_cuboids.push_back(label_cuboids_[label_index]);
		}
		else
		{
			for (std::vector<MeshCuboid *>::iterator it = label_cuboids_[label_index].begin();
				it != label_cuboids_[label_index].end(); ++it)
				delete (*it);
		}
	}

	labels_.swap(new_labels);
	label_names_.swap(new_label_names);
	label_symmetries_.swap(new_label_symmetries);
	label_cuboids_.swap(new_label_cuboids);
	label_children_.swap(new_label_children);
}

void MeshCuboidStructure::add_symmetric_group_labels()
{
	clear_symmetric_group_labels();

	const unsigned int num_given_labels = num_labels();
	assert(label_symmetries_.size() == num_given_labels);

	unsigned int num_new_labels = 0;

	bool *is_label_added = new bool[num_given_labels];
	memset(is_label_added, false, num_given_labels * sizeof(bool));

	for (LabelIndex label_index = 0; label_index < num_given_labels; ++label_index)
	{
		if (is_label_added[label_index])
			continue;

		std::list<LabelIndex> children_label_indices = label_symmetries_[label_index];
		if (children_label_indices.empty())
			continue;
		
		children_label_indices.push_front(label_index);

		for (std::list<LabelIndex>::const_iterator it = children_label_indices.begin();
			it != children_label_indices.end(); ++it)
		{
			LabelIndex n_label_index = (*it);
			assert(n_label_index < num_given_labels);
			is_label_added[n_label_index] = true;
		}


		std::stringstream new_label_name;
		new_label_name << "symmetry_group_" << num_new_labels;

		labels_.push_back(get_new_label());
		label_names_.push_back(new_label_name.str());
		label_symmetries_.push_back(std::list<LabelIndex>());
		label_cuboids_.push_back(std::vector<MeshCuboid*>());
		label_children_.push_back(children_label_indices);

		++num_new_labels;
	}

	delete[] is_label_added;

	// NOTE:
	// Draws all points.
	query_label_index_ = num_labels();
}

void MeshCuboidStructure::create_symmetric_group_cuboids()
{
	for (unsigned int label_index = 0; label_index < num_labels(); ++label_index)
	{
		if (!is_label_group(label_index))
			continue;

		std::vector<MeshCuboid *> children_cuboids;

		for (std::list<LabelIndex>::const_iterator it = label_children_[label_index].begin();
			it != label_children_[label_index].end(); ++it)
		{
			LabelIndex n_label_index = (*it);
			assert(n_label_index < num_labels());

			children_cuboids.insert(children_cuboids.end(),
				label_cuboids_[n_label_index].begin(), label_cuboids_[n_label_index].end());
		}

		if (!children_cuboids.empty())
		{
			MeshCuboid* merged_cuboid = MeshCuboid::merge_cuboids(label_index, children_cuboids);
			assert(merged_cuboid);
			merged_cuboid->set_group_cuboid(true);

			// Clear sample points.
			merged_cuboid->clear_sample_points();

			label_cuboids_[label_index].clear();
			label_cuboids_[label_index].push_back(merged_cuboid);
		}
	}
}

int MeshCuboidStructure::find_parent_label_index(
	const LabelIndex _label_index_1, const LabelIndex _label_index_2)
{
	assert(label_children_.size() == num_labels());

	for (LabelIndex label_index = 0; label_index < num_labels(); ++label_index)
	{
		bool is_child_label_index[2];
		memset(is_child_label_index, false, 2 * sizeof(bool));

		for (std::list<LabelIndex>::iterator it = label_children_[label_index].begin();
			it != label_children_[label_index].end(); ++it)
		{
			if ((*it) == _label_index_1) is_child_label_index[0] = true;
			if ((*it) == _label_index_2) is_child_label_index[1] = true;
		}

		if (is_child_label_index[0] && is_child_label_index[1])
			return label_index;
	}

	return -1;
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
*/