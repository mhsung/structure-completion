#include "MyMesh.h"

#include <fstream>
#include <iostream>
#include <limits>

//#include "ConvertFromOpenMesh.h"

MyMesh::MyMesh()
: prop_curvatures_(false)
, prop_face_areas_(false)
, prop_vertex_areas_(false)
, seed_vertex_index_(-1)
, query_vertex_index_(-1)
, bbox_center_(0.0)
, bbox_size_(0.0)
, object_diameter_(0.0)
, translation_(0.0)
, scale_(1.0)
, mesh_coloring_option_(FACE_COLOR)
{
	add_property(vertex_color_map_value_);
	add_property(face_color_map_value_);

	add_property(vertex_label_);
	add_property(face_label_);
};

MyMesh::~MyMesh()
{
	remove_property(vertex_color_map_value_);
	remove_property(face_color_map_value_);

	remove_property(vertex_label_);
	remove_property(face_label_);

	release_curvatures();
	release_face_areas();
	release_vertex_areas();
}

void MyMesh::clear()
{
	OpenMesh::TriMesh_ArrayKernelT<MyTraits>::clear();

	// Clear properties.
	MyMesh::ConstVertexIter v_it, v_start(vertices_begin()), v_end(vertices_end());
	MyMesh::ConstFaceIter f_it, f_start(faces_begin()), f_end(faces_end());

	// Clear vertex colors (set white).
	for (v_it = v_start; v_it != v_end; ++v_it)
		set_color(v_it.handle(), MyMesh::Color(255, 255, 255));

	// Clear face colors (set white).
	for (f_it = f_start; f_it != f_end; ++f_it)
		set_color(f_it.handle(), MyMesh::Color(255, 255, 255));

	// Clear curvatures.
	release_curvatures();

	// Clear color map values.
	for (v_it = v_start; v_it != v_end; ++v_it)
		property(vertex_color_map_value_, v_it) = 0.0;

	for (f_it = f_start; f_it != f_end; ++f_it)
		property(face_color_map_value_, f_it) = 0.0;

	// Clear labels.
	for (v_it = v_start; v_it != v_end; ++v_it)
		property(vertex_label_, v_it) = -1;

	for (f_it = f_start; f_it != f_end; ++f_it)
		property(face_label_, f_it) = -1;

	// Clear user-defined vertex indices.
	seed_vertex_index_ = -1;
	query_vertex_index_ = -1;
	feature_vertex_indices_.clear();

	// Clear bounding box.
	bbox_center_ = MyMesh::Point(0.0);
	bbox_size_ = MyMesh::Normal(0.0);
	object_diameter_ = 0.0;

	// Clear transformation.
	translation_ = MyMesh::Normal(0.0);
	scale_ = 1.0;
}

void MyMesh::initialize()
{
	// Check unreferenced vertices
	//bool *referenced_vertex = new bool[n_vertices()];
	//for (vertexIndex vid = 0; vid < n_vertices(); vid++)
	//	referenced_vertex[vid] = false;
	//for (e_it = edges_begin(); e_it != e_end; ++e_it)
	//{
	//	HalfedgeHandle heh = halfedge_handle(e_it.handle(), 0);	// 1st halfedge
	//	referenced_vertex[from_vertex_handle(heh).idx()] = true;
	//	referenced_vertex[to_vertex_handle(heh).idx()] = true;
	//}
	//for (vertexIndex vid = 0; vid < n_vertices(); vid++)
	//{
	//	if (!referenced_vertex[vid])
	//		std::cout << "Unreferenced vertex: " << vid << std::endl;
	//}
	//delete[] referenced_vertex;

	//make_face_normal_consistent();
	

	// Compute bounding box.
	MyMesh::Point bbox_min, bbox_max;
	bbox_min = bbox_max = OpenMesh::vector_cast<MyMesh::Point>(point(vertices_begin()));

	for (ConstVertexIter v_it = vertices_begin(); v_it != vertices_end(); ++v_it)
	{
		bbox_min.minimize(OpenMesh::vector_cast<MyMesh::Point>(point(v_it)));
		bbox_max.maximize(OpenMesh::vector_cast<MyMesh::Point>(point(v_it)));
	}

	bbox_center_ = 0.5 * (bbox_min + bbox_max);
	bbox_size_ = bbox_max - bbox_min;


	// Compute object diameter.
	object_diameter_ = 0.0;

	for (ConstVertexIter v_it = vertices_begin(); v_it != vertices_end(); ++v_it)
	{
		Real distance = (point(v_it) - bbox_center_).length();
		object_diameter_ = std::max(object_diameter_, distance);
	}

	object_diameter_ = 2 * object_diameter_;
}

void MyMesh::request_vertex_areas()
{
	if (!prop_vertex_areas_)
	{
		request_face_areas();
		add_property(vertex_area_);

		// Integrate adjacent face area into each vertex.
		for (MyMesh::ConstVertexIter v_it = vertices_begin(); v_it != vertices_end(); ++v_it)
		{
			MyMesh::VertexHandle vh = v_it.handle();
			property(vertex_area_, v_it) = 0.0f;
			unsigned int num_adjacent_faces = 0;

			for (MyMesh::ConstVertexFaceIter vf_it = cvf_iter(vh); vf_it; ++vf_it)
			{
				MyMesh::FaceHandle fh = vf_it.handle();
				property(vertex_area_, v_it) += property(face_area_, vf_it);
				++num_adjacent_faces;
			}

			property(vertex_area_, v_it) /= num_adjacent_faces;
		}
	}
}

void MyMesh::release_vertex_areas()
{
	if (prop_vertex_areas_)
	{
		remove_property(vertex_area_);
		prop_vertex_areas_ = false;
	}
}

void MyMesh::request_face_areas()
{
	if (!prop_face_areas_)
	{
		add_property(face_area_);

		// Compute each face area.
		for (MyMesh::ConstFaceIter f_it = faces_begin(); f_it != faces_end(); ++f_it)
		{
			MyMesh::FaceHandle fh = f_it.handle();
			MyMesh::VertexHandle faceVertices[3];

			MyMesh::ConstFaceVertexIter fv_it = cfv_iter(fh);
			for (unsigned int i = 0; fv_it && i < 3; ++fv_it, ++i)	// Trimesh
				faceVertices[i] = fv_it.handle();

			MyMesh::Normal face_edge_1 = point(faceVertices[1]) - point(faceVertices[0]);
			MyMesh::Normal face_edge_2 = point(faceVertices[2]) - point(faceVertices[0]);
			MyMesh::Normal face_normal = cross(face_edge_1, face_edge_2);

			property(face_area_, f_it) = (float)face_normal.norm() / 2.0f;
		}
	}
}

void MyMesh::release_face_areas()
{
	if (prop_face_areas_)
	{
		remove_property(face_area_);
		prop_face_areas_ = false;
	}
}

void MyMesh::request_curvatures()
{
	if (!prop_curvatures_)
	{
		add_property(principal_curvature_1_);
		add_property(principal_curvature_2_);
		add_property(principal_direction_1_);
		add_property(principal_curvature_2_);
		add_property(mean_curvature_);
		add_property(gaussian_curvature_);

		//TriMesh *tri_mesh = TriMeshPack<MyTraits>::need_curvatures(this, CURVATURE_SMOOTHING);

		//for(MyMesh::VertexIter v_it = this->vertices_begin(); v_it != this->vertices_end(); ++v_it)
		//{
		//	unsigned int vid = v_it.handle().idx();
		//	property(principal_curvature_1_, v_it) = tri_mesh->curv1[vid];
		//	property(principal_curvature_2_, v_it) = tri_mesh->curv2[vid];
		//	property(principal_direction_1_, v_it) = (Normal)tri_mesh->pdir1[vid];
		//	property(principal_direction_2_, v_it) = (Normal)tri_mesh->pdir2[vid];
		//	property(mean_curvature_, v_it) = (tri_mesh->curv1[vid] + tri_mesh->curv2[vid]) / 2.0f;
		//	property(gaussian_curvature_, v_it) = tri_mesh->curv1[vid] * tri_mesh->curv2[vid];
		//}
		//delete tri_mesh;

		prop_curvatures_ = true;
	}
}

void MyMesh::release_curvatures()
{
	if (prop_curvatures_)
	{
		remove_property(principal_curvature_1_);
		remove_property(principal_curvature_2_);
		remove_property(principal_direction_1_);
		remove_property(principal_curvature_2_);
		remove_property(mean_curvature_);
		remove_property(gaussian_curvature_);

		prop_curvatures_ = false;
	}
}

void MyMesh::translate(const Normal _translate)
{
	for (MyMesh::ConstVertexIter v_it = vertices_begin(); v_it != vertices_end(); ++v_it)
	{
		MyMesh::Point &p = point(v_it.handle());
		p = p + _translate;
	}

	translation_ += _translate;
	bbox_center_ += _translate;
}

void MyMesh::scale(const Real _scale)
{
	assert(_scale > 0);
	for (MyMesh::ConstVertexIter v_it = vertices_begin();
		v_it != vertices_end(); ++v_it)
	{
		MyMesh::Point &p = point(v_it.handle());
		p = p * _scale;
	}

	scale_ *= _scale;
	translation_ *= _scale;
	bbox_center_ *= _scale;
	bbox_size_ *= _scale;
	object_diameter_ *= _scale;
}

void MyMesh::reset_transformation()
{
	if (translation_ != MyMesh::Point(0.0) || scale_ != 1.0)
	{
		scale(1.0 / scale_);
		translate(-translation_);
	}

	assert(translation_ == MyMesh::Point(0.0));
	assert(scale_ == 1.0);
}

bool MyMesh::load_feature_vertices( const char *_filename, bool _verbose )
{
	std::ifstream file(_filename);
	if (!file)
	{
		std::cerr << "Can't open file: \"" << _filename << "\"" << std::endl;
		return false;
	}

	if (_verbose)
		std::cout << "Loading " << _filename << "..." << std::endl;

	feature_vertex_indices_.clear();

	char buffer[256];
	char text[256];
	VertexIndex vid = 0;

	for(unsigned int count = 0; !file.eof(); count++)
	{
		file.getline( buffer, 256 );

		int ret = sscanf(buffer, "%d, %s\n", &vid, text);
		if (ret < 1) continue;
		// NOTE:
		// Negative index indicates undefined.
		else if (vid < 0) continue;

		if (_verbose)
			std::cout << vid << ": " << text << std::endl;

		// Check the uniqueness.
		//bool is_unique_index = true;
		//for(VertexIndexArray::iterator it = feature_vertex_indices_.begin();
		//	it != feature_vertex_indices_.end(); ++it)
		//{
		//	if (*it == vid)
		//	{
		//		is_unique_index = false;
		//		break;
		//	}
		//}
		//if (is_unique_index)
		//	feature_vertex_indices_.push_back(vid);

		feature_vertex_indices_.push_back(vid);
	}
	file.close();

	if (_verbose)
		std::cout << "Done." << std::endl;

	// Eliminate duplicate features.
	//std::sort(feature_vertex_indices_.begin(), feature_vertex_indices_.end());
	//feature_vertex_indices_.erase(std::unique(feature_vertex_indices_.begin(),
	//	feature_vertex_indices_.end()), feature_vertex_indices_.end());

	return true;
}

bool MyMesh::save_feature_vertices( const char *_filename, bool _verbose ) const
{
	if(feature_vertex_indices_.empty())
	{
		std::cerr << "No feature vertex exists." << std::endl;
		return false;
	}

	std::ofstream file(_filename);
	if (!file)
	{
		std::cerr << "Can't save file: \"" << _filename << "\"" << std::endl;
		return false;
	}

	if (_verbose)
		std::cout << "Saving " << _filename << "..." << std::endl;

	for(VertexIndexArray::const_iterator it = feature_vertex_indices_.begin();
		it != feature_vertex_indices_.end(); ++it)
	{
		VertexIndex vid = (*it);
		// NOTE:
		// Negative index indicates undefined.
		assert(vid >= 0);

		if(_verbose)
			std::cout << vid << "," << property(vertex_color_map_value_)[*it] << std::endl;

		file << vid << "," << property(vertex_color_map_value_)[*it] << std::endl;
	}

	file.close();

	if (_verbose)
		std::cout << "Done." << std::endl;

	return true;
}

bool MyMesh::load_color_map(const char *_filename, bool _verbose)
{
	std::ifstream file(_filename);
	if (!file)
	{
		std::cerr << "Can't open file: \"" << _filename << "\"" << std::endl;
		return false;
	}

	if (_verbose)
		std::cout << "Loading " << _filename << "..." << std::endl;

	char buffer[256];
	unsigned int count = 0;
	float var = 0;

	unsigned int num_values = 0;
	if( mesh_coloring_option_ == VERTEX_COLOR ) num_values = n_vertices();
	else if( mesh_coloring_option_ == FACE_COLOR ) num_values = n_faces();

	for(; !file.eof(); count++)
	{
		file.getline( buffer, 256 );

		int ret = sscanf(buffer, "%f\n", &var);
		if(ret < 1) continue;

		if(count >= num_values)
		{
			std::cout << "Warning: Too many values." << std::endl;
			break;
		}

		if( mesh_coloring_option_ == VERTEX_COLOR ) property(vertex_color_map_value_)[count] = var;
		else if( mesh_coloring_option_ == FACE_COLOR ) property(face_color_map_value_)[count] = var;
	}

	if (count < num_values)
		std::cout << "Warning: Number of values are not enough." << std::endl;

	file.close();

	//
	if( mesh_coloring_option_ == VERTEX_COLOR ) set_vertex_color_map();
	else if( mesh_coloring_option_ == FACE_COLOR ) set_face_color_map();
	//

	if (_verbose)
		std::cout << "Done." << std::endl;

	return true;
}

bool MyMesh::save_color_map(const char *_filename, bool _verbose) const
{
	std::ofstream file(_filename);
	if (!file)
	{
		std::cerr << "Can't save file: \"" << _filename << "\"" << std::endl;
		return false;
	}

	if (_verbose)
		std::cout << "Saving " << _filename << "..." << std::endl;

	//
	if( mesh_coloring_option_ == VERTEX_COLOR )
	{
		for(MyMesh::VertexIter v_it = this->vertices_begin(); v_it != this->vertices_end(); ++v_it)
			file << property(vertex_color_map_value_, v_it) << std::endl;
	}
	else if( mesh_coloring_option_ == FACE_COLOR )
	{
		for(MyMesh::FaceIter f_it = this->faces_begin(); f_it != this->faces_end(); ++f_it)
			file << property(face_color_map_value_, f_it) << std::endl;
	}
	//

	file.close();

	if (_verbose)
		std::cout << "Done." << std::endl;

	return true;
}

bool MyMesh::load_vertex_color_map(const char *_filename, bool _verbose)
{
	mesh_coloring_option_ = VERTEX_COLOR;
	bool ret = load_color_map(_filename, _verbose);
	return ret;
}

bool MyMesh::save_vertex_color_map(const char *_filename, bool _verbose) const
{
	assert(mesh_coloring_option_ == VERTEX_COLOR);
	bool ret = save_color_map(_filename, _verbose);
	return ret;
}

bool MyMesh::load_face_color_map(const char *_filename, bool _verbose)
{
	mesh_coloring_option_ = FACE_COLOR;
	bool ret = load_color_map(_filename, _verbose);
	return ret;
}

bool MyMesh::save_face_color_map(const char *_filename, bool _verbose) const
{
	assert(mesh_coloring_option_ == FACE_COLOR);
	bool ret = save_color_map(_filename, _verbose);
	return ret;
}

bool MyMesh::load_face_label_simple(const char *_filename, bool _verbose)
{
	std::ifstream file(_filename);
	if (!file)
	{
		std::cerr << "Can't open file: \"" << _filename << "\"" << std::endl;
		return false;
	}

	if (_verbose)
		std::cout << "Loading " << _filename << "..." << std::endl;

	char buffer[256];
	unsigned int count = 0;
	Label label = 0;

	unsigned int num_faces = n_faces();
	for(; !file.eof(); count++)
	{
		file.getline( buffer, 256 );

		int ret = sscanf(buffer, "%d\n", &label);
		if (ret < 1) continue;
		// NOTE:
		// Negative index indicates undefined.
		else if (label < 0) continue;

		if (count >= num_faces)
		{
			std::cout << "Warning: Too many values." << std::endl;
			break;
		}

		property(face_label_)[count] = label;
	}

	if(count < num_faces)
		std::cout << "Warning: Number of values are not enough" << std::endl;

	file.close();

	//set_vertex_label_from_face_label();
	set_face_label_colors();

	std::cout << "Done." << std::endl;

	return true;
}

bool MyMesh::save_face_label_simple(const char *_filename, bool _verbose) const
{
	std::ofstream file(_filename);
	if (!file)
	{
		std::cerr << "Can't save file: \"" << _filename << "\"" << std::endl;
		return false;
	}

	if (_verbose)
		std::cout << "Saving " << _filename << "..." << std::endl;

	for (MyMesh::FaceIter f_it = this->faces_begin(); f_it != this->faces_end(); ++f_it)
	{
		Label label = property(face_label_, f_it);
		// NOTE:
		// Negative index indicates undefined.
		assert(label >= 0);

		file << label << std::endl;
	}

	file.close();

	if (_verbose)
		std::cout << "Done." << std::endl;

	return true;
}

bool MyMesh::load_face_label(const char *_filename, bool _verbose)
{
	std::ifstream file(_filename);
	if (!file)
	{
		std::cerr << "Can't open file: \"" << _filename << "\"" << std::endl;
		return false;
	}

	if (_verbose)
		std::cout << "Loading " << _filename << "..." << std::endl;

	int num_faces = n_faces();

	std::string buffer;
	Label label = 0;
	bool read_label = true;

	bool *is_face_assigned = new bool[num_faces];
	memset(is_face_assigned, false, num_faces * sizeof(bool));
	
	while (!file.eof())
	{
		std::getline(file, buffer);

		if (read_label)
		{
			label = atoi(buffer.c_str());
			// NOTE:
			// Negative index indicates undefined.
			if (label < 0)
			{
				std::cerr << "Wrong file format: \"" << _filename << "\"" << std::endl;
				return false;
			}

			read_label = false;
		}
		else
		{
			std::stringstream strstr(buffer);
			std::string token;

			for (int count = 0; std::getline(strstr, token, ' '); count++)
			{
				// NOTE:
				// The face index starts from '1' in OBJ format.
				FaceIndex fid = atoi(token.c_str()) - 1;
				assert(fid < static_cast<int>(n_faces()));
				property(face_label_)[fid] = label;
				is_face_assigned[fid] = true;
			}

			read_label = true;
		}
	}

	for (FaceIndex fid = 0; fid < num_faces; fid++)
	{
		if (!is_face_assigned[fid])
		{
			std::cout << "Warning: Label is not assigned to the face ID (" << fid << ")." << std::endl;
		}
	}


	delete [] is_face_assigned;
	file.close();

	//set_vertex_label_from_face_label();
	set_face_label_colors();

	if (_verbose)
		std::cout << "Done." << std::endl;

	return true;
}

bool MyMesh::save_face_label(const char *_filename, bool _verbose) const
{
	std::ofstream file(_filename);
	if (!file)
	{
		std::cerr << "Can't open saving file: \"" << _filename << "\"" << std::endl;
		return false;
	}

	if (_verbose)
		std::cout << "Saving " << _filename << "..." << std::endl;

	std::list<std::vector<FaceIndex>> all_label_faces;
	get_all_label_faces(all_label_faces);

	if (_verbose)
		std::cout << "Number of face clusters: " << all_label_faces.size() << std::endl;

	for (std::list<std::vector<FaceIndex>>::iterator l_it = all_label_faces.begin();
		l_it != all_label_faces.end(); ++l_it)
	{
		std::vector<FaceIndex> &label_faces = (*l_it);
		assert(!label_faces.empty());

		FaceIndex fid = label_faces.front();
		MyMesh::FaceHandle fh = face_handle(fid);
		Label label = property(face_label_, fh);

		file << label << std::endl;

		MyMesh::Point center = MyMesh::Point(0.0);

		for (std::vector<FaceIndex>::iterator f_it = label_faces.begin();
			f_it != label_faces.end(); ++f_it)
		{
			fid = (*f_it);
			// NOTE:
			// The face index starts from '1' in OBJ format.
			file << fid + 1 << " ";
		}

		file << std::endl;
	}

	file.close();

	if (_verbose)
		std::cout << "Done." << std::endl;

	return true;
}

void MyMesh::get_all_label_faces(std::list<std::vector<FaceIndex>> &_all_label_faces) const
{
	_all_label_faces.clear();

	for (MyMesh::FaceIter f_it = faces_begin(); f_it != faces_end(); ++f_it)
	{
		MyMesh::FaceHandle fh = f_it.handle();
		FaceIndex fid = fh.idx();
		int label = property(face_label_, fh);

		// NOTE:
		// Negative index indicates undefined.
		if (label < 0)
			continue;

		bool is_new_label = true;

		for (std::list<std::vector<FaceIndex>>::iterator c_it = _all_label_faces.begin();
			c_it != _all_label_faces.end(); ++c_it)
		{
			std::vector<FaceIndex> &face_cluster = (*c_it);
			assert(!face_cluster.empty());

			FaceIndex n_fid = face_cluster[0];
			MyMesh::FaceHandle n_fh = face_handle(n_fid);
			int n_label = property(face_label_, n_fh);
			if (n_label == label)
			{
				face_cluster.push_back(fid);
				is_new_label = false;
			}
		}

		if (is_new_label)
		{
			std::vector<FaceIndex> face_cluster;
			face_cluster.reserve(n_faces());
			face_cluster.push_back(fid);
			_all_label_faces.push_back(face_cluster);
		}
	}
}

void MyMesh::set_vertex_label_from_face_label()
{
	for (MyMesh::ConstFaceIter f_it = faces_begin(); f_it != faces_end(); ++f_it)
	{
		MyMesh::FaceHandle fh = f_it.handle();

		// NOTE:
		// When a vertex is on the boundary of two different label faces,
		// anyone of them is chosen by the traversing order.
		for (MyMesh::ConstFaceVertexIter fv_it = cfv_iter(fh); fv_it; ++fv_it)
		{
			MyMesh::VertexHandle vh = fv_it.handle();
			property(vertex_label_, vh) = property(face_label_, fh);
		}
	}
}

void MyMesh::get_edge_neighbor_vertex_indices(const unsigned int _eid,
	VertexIndex &_v1, VertexIndex &_v2) const
{
	MyMesh::HalfedgeHandle heh = halfedge_handle(edge_handle(_eid), 0);
	assert(edge_handle(heh).idx() == _eid);
	_v1 = from_vertex_handle(heh).idx();
	_v2 = to_vertex_handle(heh).idx();
}

MyMesh::Color MyMesh::get_label_color(const Label _label)
{
	// NOTE:
	// Negative index indicates undefined.
	assert(_label >= 0);

	srand(LABEL_COLORING_RANDOM_SEED * _label);

	Real r = rand() / (Real)RAND_MAX;
	Real g = rand() / (Real)RAND_MAX;
	Real b = rand() / (Real)RAND_MAX;

	MyMesh::Color label_color(
		(unsigned char)(r * 255), (unsigned char)(g * 255), (unsigned char)(b * 255));

	return label_color;
}

void MyMesh::gray_to_rgb_color(const Real _gray, Real &_r, Real &_g, Real &_b)
{
	// Assume that the input gray value and the output RGB value has range [0, 1].
	Real value = std::min(std::max(_gray, 0.0), 1.0);

	// Inverse gray value.
	value = 1.0 - value;

	// Red value.
	if (value >= 128 && value <= 192)	_r = (value - 128.0) / 65.0;
	else if (value > 192)	_r = 1.0;
	else	_r = 0.0;

	// Green value.
	if (value <= 64)	_g = value / 64.0;
	else if (value > 64 && value <= 192)	_g = 1.0;
	else	_g = 1.0 - ((value - 192.0) / 65.0);

	// Blue value.
	if (value <= 64)	_b = 1.0;
	else if (value > 64 && value <= 128)	_b = 1.0 - ((value - 64.0) / 65.0);
	else	_b = 0.0;
}

void MyMesh::clear_colors()
{
	MyMesh::ConstVertexIter v_it, v_start(vertices_begin()), v_end(vertices_end());
	MyMesh::ConstFaceIter f_it, f_start(faces_begin()), f_end(faces_end());

	// Clear vertex colors (set white).
	for (v_it = v_start; v_it != v_end; ++v_it)
		set_color(v_it.handle(), MyMesh::Color(255, 255, 255));

	// Clear face colors (set white).
	for (f_it = f_start; f_it != f_end; ++f_it)
		set_color(f_it.handle(), MyMesh::Color(255, 255, 255));
}

void MyMesh::set_vertex_color_map(const RealArray &_values)
{
	assert(_values.size() >= n_vertices());

	for(MyMesh::VertexIter v_it = this->vertices_begin(); v_it != this->vertices_end(); ++v_it)
		property(vertex_color_map_value_, v_it) = (float)_values[ v_it.handle().idx() ];
		
	set_vertex_color_map();
}

void MyMesh::set_vertex_color_map()
{
	mesh_coloring_option_ = VERTEX_COLOR;

	// Find minimum and maximum for given values.
	Real min_color_map_value = std::numeric_limits<Real>::max();
	Real max_color_map_value = -std::numeric_limits<Real>::max();
	for(MyMesh::VertexIter v_it = this->vertices_begin(); v_it != this->vertices_end(); ++v_it)
	{
		if (property(vertex_color_map_value_, v_it) > max_color_map_value)
			max_color_map_value = property(vertex_color_map_value_, v_it);

		if (property(vertex_color_map_value_, v_it) < min_color_map_value)
			min_color_map_value = property(vertex_color_map_value_, v_it);
	}

	for(MyMesh::VertexIter v_it = this->vertices_begin(); v_it != this->vertices_end(); ++v_it)
	{
		VertexIndex vid = v_it.handle().idx();
		Real color_map_value = (property(vertex_color_map_value_, v_it) - min_color_map_value)
			/ (max_color_map_value - min_color_map_value);

		Real red = 0.0, green = 0.0, blue = 0.0;
		gray_to_rgb_color(color_map_value, red, green, blue);

		this->set_color(vertex_handle(vid), MyMesh::Color(
			(unsigned char)(red * 255), (unsigned char)(green * 255), (unsigned char)(blue * 255)));
	}
}

void MyMesh::set_face_color_map(const RealArray &_values)
{
	assert(_values.size() >= n_faces());

	for(MyMesh::FaceIter f_it = this->faces_begin(); f_it != this->faces_end(); ++f_it)
		property(face_color_map_value_, f_it) = (float)_values[ f_it.handle().idx() ];

	set_face_color_map();
}

void MyMesh::set_face_color_map()
{
	mesh_coloring_option_ = FACE_COLOR;

	// Find minimum and maximum for given values.
	Real min_color_map_value = std::numeric_limits<Real>::max();
	Real max_color_map_value = -std::numeric_limits<Real>::max();
	for(MyMesh::FaceIter f_it = this->faces_begin(); f_it != this->faces_end(); ++f_it)
	{
		if (property(face_color_map_value_, f_it) > max_color_map_value)
			max_color_map_value = property(face_color_map_value_, f_it);
		if (property(face_color_map_value_, f_it) < min_color_map_value)
			min_color_map_value = property(face_color_map_value_, f_it);
	}
	

	for(MyMesh::FaceIter f_it = this->faces_begin(); f_it != this->faces_end(); ++f_it)
	{
		FaceIndex fid = f_it.handle().idx();
		Real color_map_value = (property(face_color_map_value_, f_it) - min_color_map_value)
			/ (max_color_map_value - min_color_map_value);

		Real red = 0.0, green = 0.0, blue = 0.0;
		gray_to_rgb_color(color_map_value, red, green, blue);

		this->set_color(face_handle(fid), MyMesh::Color(
			(unsigned char)(red * 255), (unsigned char)(green * 255), (unsigned char)(blue * 255)));
	}
}

void MyMesh::set_face_label_colors(bool _verbose)
{
	mesh_coloring_option_ = FACE_COLOR;

	std::list<std::vector<FaceIndex>> all_label_faces;
	get_all_label_faces(all_label_faces);

	if (_verbose)
		std::cout << "Number of labels: " << all_label_faces.size() << std::endl;

	for (std::list<std::vector<FaceIndex>>::iterator l_it = all_label_faces.begin();
		l_it != all_label_faces.end(); ++l_it)
	{
		std::vector<FaceIndex> &label_faces = (*l_it);
		assert(!label_faces.empty());

		FaceIndex fid = label_faces.front();
		MyMesh::FaceHandle fh = face_handle(fid);
		Label label = property(face_label_, fh);

		MyMesh::Color label_color = get_label_color(label);
		if (_verbose)
		{
			printf(" - label # = %d, color = (%d, %d, %d)\n",
				label, label_color[0], label_color[1], label_color[2]);
		}

		for (std::vector<FaceIndex>::iterator f_it = label_faces.begin();
			f_it != label_faces.end(); ++f_it)
		{
			set_color(face_handle((*f_it)), label_color);
		}
	}
}

void MyMesh::extract_local_min_max_feature_vertices(const RealArray &_values,
	const unsigned int _local_search_depth,
	const bool _extract_min, const bool _extract_max)
{
	assert(_values.size() >= n_vertices());
	feature_vertex_indices_.clear();

	for(MyMesh::VertexIter v_it = this->vertices_begin(); v_it != this->vertices_end(); ++v_it)
	{
		unsigned int vid = v_it.handle().idx();

		bool _is_local_min = true;
		bool _is_local_max = true;

		is_local_min_max(_values, vid, vid, _local_search_depth, _is_local_min, _is_local_max);

		if(_extract_min && _is_local_min)
		{
			feature_vertex_indices_.push_back(vid);
		}
		else if(_extract_max && _is_local_max)
		{
			feature_vertex_indices_.push_back(vid);
		}
	}
}

void MyMesh::is_local_min_max(const RealArray &_values,
	const VertexIndex _vid, const VertexIndex _n_vid,
	const unsigned int _local_search_depth,
	bool &_is_local_min, bool &_is_local_max) const
{
	// 'vid': To be tested.
	// 'n_vid': 'vid' itself or one of neighbors of 'vid'.
	// '_local_search_depth': Number of recursive call.

	unsigned int count_neighbors = 0;
	MyMesh::VertexHandle vh = vertex_handle(_n_vid);

	for(MyMesh::ConstVertexVertexIter vv_it = this->cvv_iter(vh); vv_it; ++vv_it, ++count_neighbors)
	{
		MyMesh::VertexHandle one_ring_vh = vv_it.handle();
		unsigned int nn_vid = one_ring_vh.idx();

		if(_vid == nn_vid) continue;
		else if( _values[_vid] >= _values[nn_vid] )	_is_local_min = false;
		else if( _values[_vid] <= _values[nn_vid] )	_is_local_max = false;

		if( !_is_local_min && !_is_local_max ) break;
	}

	if(_vid == _n_vid && count_neighbors == 0)
	{
		// Ignore isolated vertex.
		_is_local_min = _is_local_max = false;
	}

	if( (_is_local_min || _is_local_max) && _local_search_depth > 1 )
	{
		for(MyMesh::ConstVertexVertexIter vv_it = this->cvv_iter(vh); vv_it; ++vv_it)
		{
			MyMesh::VertexHandle one_ring_vh = vv_it.handle();
			unsigned int nn_vid = one_ring_vh.idx();

			// Call recursively.
			is_local_min_max(_values, _vid, nn_vid, _local_search_depth - 1,
				_is_local_min, _is_local_max);
		}
	}
}

void MyMesh::extract_zero_value_feature_vertices(const RealArray &_values,
	const Real _zero_threshold)
{
	// Less than zero_threshold portion of the absolute maximum value is defined as zero.
	assert(_values.size() >= n_vertices());
	feature_vertex_indices_.clear();

	// Find the absolute maximum.
	Real abs_max_value = 0.0, threshold = 0.0;
	for(MyMesh::VertexIter v_it = this->vertices_begin(); v_it != this->vertices_end(); ++v_it)
	{
		unsigned int vid = v_it.handle().idx();
		if(abs(_values[vid]) > abs_max_value)
			abs_max_value = abs(_values[vid]);
	}
	threshold = abs_max_value * _zero_threshold;

	for(MyMesh::VertexIter v_it = this->vertices_begin(); v_it != this->vertices_end(); ++v_it)
	{
		unsigned int vid = v_it.handle().idx();
		if(abs(_values[vid]) < threshold)
			feature_vertex_indices_.push_back(vid);
	}
}

void MyMesh::make_face_normal_consistent()
{
	request_face_normals();

	int num_edges = n_edges();
	int num_faces = n_faces();

	bool *is_edge_inverted = new bool[num_edges];
	memset(is_edge_inverted, false, num_edges * sizeof(bool));


	// Step #1: Check inverted adjacent face normals.
	for (EdgeIndex eid = 0; eid < num_edges; eid++)
	{
		EdgeHandle eh = edge_handle(eid);

		// Vertices on the adjacent edge.
		VertexHandle vh[2];	
		FaceHandle fh[2];

		// Face normals.
		Normal n[2];	

		// Vertices that are NOT on the adjacent edge, but on the opposite sides.
		// Original Positions.
		Point p[2];			


		// NOTE:
		// When assuming that two normal vectors are the same,
		// two faces are on the same plane.
		// Let assume that this plane is projected 2D plane,
		// and vh[0] is the origin, and (vh[1]-vh[0]) is y-axis direction.

		// Position with local coordinates.
		Point proj_p[2];	
		bool ret;

		for (unsigned int i = 0; i < 2; i++)
		{
			vh[i] = from_vertex_handle(halfedge_handle(eh, i));
			fh[i] = face_handle(halfedge_handle(eh, i));
			if (fh[i].idx() >= 0)	n[i] = normal(fh[i]);
		}

		if (vh[0].idx() < 0 || vh[1].idx() < 0
			|| fh[0].idx() < 0 || fh[1].idx() < 0)
			continue;

		Normal dir = (point(vh[1]) - point(vh[0])).normalize();

		for (unsigned int i = 0; i < 2; i++)
		{
			ret = false;
			for (ConstFaceVertexIter fv_it = cfv_iter(fh[i]); fv_it; ++fv_it)
			{
				if (fv_it.handle() != vh[0] && fv_it.handle() != vh[1])
				{
					p[i] = point(fv_it.handle());
					ret = true;
					break;
				}
			}
			assert(ret);

			Normal v = p[i] - point(vh[0]);

			// Projection (parallel) -> y
			Real y = dot(v, dir);
			Real x = sqrt(v.sqrnorm() - y*y);

			proj_p[i] = Point(x, y, 0);
		}

		Real dist = (p[1] - p[0]).norm();
		Real cos_angle = dot(n[0], n[1]);
		Real sin_angle = sqrt(1 - cos_angle*cos_angle);
		Point folded_p1;


		// Now, compute how p[1] is moved by two normal vectors.
		// The face f[1] is flipped by the normal difference (flipping axis is y-axis).
		folded_p1 = Point(
			cos_angle*proj_p[1][0],		// cos
			proj_p[1][1],
			sin_angle*proj_p[1][0]		// sin
			);

		// Estimated distance using current normal vectors.
		Real curr_normal_dist = (folded_p1 - proj_p[0]).norm();
		Real curr_normal_diff = curr_normal_dist - dist;

		// Consider the case that one of normal vector is inverted.
		// (angle) -> (-angle), (cos) -> (-cos), (sin) -> (sin)
		folded_p1 = Point(
			-cos_angle*proj_p[1][0],	// -cos
			proj_p[1][1],
			sin_angle*proj_p[1][0]		// sin
			);

		// Estimated distance using flipped current normal vectors.
		Real flip_normal_dist = (folded_p1 - proj_p[0]).norm();
		Real flip_normal_diff = flip_normal_dist - dist;

		if (curr_normal_diff > flip_normal_diff)
			is_edge_inverted[eid] = true;
	}


	// Step #2: Gather non-inverted adjacent faces and invert face normals of minority.
	typedef enum {
		SUBSET_A = 0,
		SUBSET_B = 1,
		VERIFIED = 2,
		NOT_VISITED = 3,
	} FaceTraversal;

	std::vector<FaceTraversal> face_subsets(num_faces, NOT_VISITED);
	unsigned int num_inverted_faces = 0;

	while (true)
	{
		FaceIndex curr_fid = 0;
		for (; curr_fid <= num_faces; curr_fid++)
			if (face_subsets[curr_fid] == NOT_VISITED)
				break;

		if (curr_fid >= num_faces)
			break;

		unsigned int num_subset_count[2] = { 0, 0 };
		
		FaceTraversal curr_subset_id = SUBSET_A;
		face_subsets[curr_fid] = curr_subset_id;
		num_subset_count[curr_subset_id]++;


		// Assume that some faces are not connected each other (multiple groups of connected faces)
		// for one connected faces group.
		bool all_neighbor_visited = false;
		while (!all_neighbor_visited)
		{
			FaceHandle fh = face_handle(curr_fid);
			all_neighbor_visited = true;

			for (ConstFaceHalfedgeIter fh_it = cfh_iter(fh); fh_it; ++fh_it)
			{
				assert(fh == face_handle(fh_it.handle()));
				HalfedgeHandle heh = opposite_halfedge_handle(fh_it.handle());

				EdgeIndex eid = edge_handle(heh).idx();
				FaceIndex n_fid = face_handle(heh).idx();

				if (face_subsets[n_fid] == NOT_VISITED)
				{
					// Flip subset ID (1 -> 2, 2 -> 1).
					if (is_edge_inverted[eid])
					{
						curr_subset_id = (curr_subset_id == SUBSET_A ? SUBSET_B : SUBSET_A);
					}

					assert(curr_subset_id == SUBSET_A || curr_subset_id == SUBSET_B);
					face_subsets[n_fid] = curr_subset_id;
					num_subset_count[curr_subset_id]++;
					curr_fid = n_fid;

					all_neighbor_visited = false;
					break;
				}
			}
		}

		// Assume that majority of normal vectors are correct.
		FaceTraversal minor_subset_id = (num_subset_count[0] > num_subset_count[1]
			? SUBSET_B : SUBSET_A);


		for (FaceIter f_it = faces_begin(); f_it != faces_end(); ++f_it)
		{
			FaceHandle fh = f_it.handle();
			FaceIndex fid = fh.idx();

			if (face_subsets[fid] == minor_subset_id)
			{
				// Invert face normal.
				set_normal(fh, -normal(fh));
				num_inverted_faces++;
			}

			if (face_subsets[fid] == SUBSET_A || face_subsets[fid] == SUBSET_B)
			{
				face_subsets[fid] = VERIFIED;
			}
		}
	}

	std::cout << num_inverted_faces << " face normals are inverted." << std::endl;

	delete [] is_edge_inverted;
}

void MyMesh::print_vertex_information(unsigned int _vid) const
{
	if(_vid > n_vertices()) return;

	std::cout << "- Index: " << _vid << std::endl;

	MyMesh::Point point = this->point(vertex_handle(_vid));
	std::cout << "- Position: " << point[0] << ", " << point[1] << ", " << point[2] << std::endl;

	MyMesh::Point normal = this->normal(vertex_handle(_vid));
	std::cout << "- Normal: " << normal[0] << ", " << normal[1] << ", " << normal[2] << std::endl;

	std::cout << "- Color map value: " << property(vertex_color_map_value_)[_vid] << std::endl;

	if (prop_curvatures_)
	{
		std::cout << "- Principal curvature 1: " << property(principal_curvature_1_)[_vid] << std::endl;
		std::cout << "- Principal curvature 2: " << property(principal_curvature_2_)[_vid] << std::endl;

		std::cout << "- Principal direction 1: " << property(principal_direction_1_)[_vid][0] << ", "
			<< property(principal_direction_1_)[_vid][1] << ", "
			<< property(principal_direction_1_)[_vid][2] << std::endl;

		std::cout << "- Principal direction 2: " << property(principal_direction_2_)[_vid][0] << ", "
			<< property(principal_direction_2_)[_vid][1] << ", "
			<< property(principal_direction_2_)[_vid][2] << std::endl;

		std::cout << "- Mean curvature: " << property(mean_curvature_)[_vid] << std::endl;

		std::cout << "- Gaussian curvature: " << property(gaussian_curvature_)[_vid] << std::endl;
	}

	std::cout << std::endl;
}

void MyMesh::print_user_defined_verices_information() const
{
	std::cout << "Seed vertex" << std::endl;
	print_vertex_information(seed_vertex_index_);
	std::cout << std::endl;

	std::cout << "Query vertex" << std::endl;
	print_vertex_information(query_vertex_index_);
	std::cout << std::endl;

	std::cout << "Feature vertices" << std::endl;
	for(VertexIndexArray::const_iterator it = feature_vertex_indices_.begin();
		it != feature_vertex_indices_.end(); ++it)
		std::cout << "[" << *it << "] : " << property(vertex_color_map_value_)[*it] << std::endl;
	std::cout << std::endl;

	std::cout << "Number of feature vertices: " << feature_vertex_indices_.size() << std::endl;
	std::cout << std::endl;
}
