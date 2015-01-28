#include "MeshViewerCore.h"

#include "glut_geometry.h"


MeshViewerCore::MeshViewerCore(GLViewerBase &_widget)
	: MeshViewerCoreT<MyMesh>(_widget)
	, point_size_(1.0)
	, selection_mode_(PICK_SEED)
	, cuboid_structure_(&mesh_)
	, draw_cuboid_axes_(true)
	, view_point_(0.0)
	, view_direction_(0.0)
{
	//// NOTE:
	//// Enable alpha channel.
	//QGLFormat fmt;
	//fmt.setAlpha(true);
	//occlusion_test_widget_ = new QGLOcculsionTestWidget(fmt);
	//draw_occlusion_test_points_ = false;

	// TEST.
	test_joint_normal_predictor_ = NULL;
	//
}

MeshViewerCore::~MeshViewerCore()
{
	//delete occlusion_test_widget_;

	// TEST.
	for (LabelIndex label_index_1 = 0; label_index_1 < joint_normal_relations_.size(); ++label_index_1)
		for (LabelIndex label_index_2 = 0; label_index_2 < joint_normal_relations_[label_index_1].size(); ++label_index_2)
			delete joint_normal_relations_[label_index_1][label_index_2];

	delete test_joint_normal_predictor_;
	//
}

bool MeshViewerCore::open_mesh(const char* _filename)
{
	bool ret = MeshViewerCoreT<MyMesh>::open_mesh(_filename);
	if (!ret) return false;

	mesh_.initialize();

	// NOTE:
	// Change scale so that the object diameter becomes 1.
	std::cout << "Object diameter: " << mesh_.get_object_diameter() << std::endl;
	std::cout << "Now scaled to become 1..." << std::endl;
	double scale = 1 / mesh_.get_object_diameter();
	mesh_.scale(scale);

	// Move mesh so that the object stands on the z = 0 plane.
	Vec3d translation = -mesh_.get_bbox_center();
	translation[2] += 0.5 * mesh_.get_bbox_size()[2];
	mesh_.translate(translation);

	Vec3f center(static_cast<float>(mesh_.get_bbox_center()[0]),
		static_cast<float>(mesh_.get_bbox_center()[1]),
		static_cast<float>(mesh_.get_bbox_center()[2]));
	float diagonal_size = static_cast<float>(mesh_.get_bbox_size().norm());
	float min_size = static_cast<float>(mesh_.get_bbox_size().min());

	// set center and radius
	set_scene_pos(center, diagonal_size * 0.5f);

	// for normal display
	normal_scale_ = min_size * 0.05f;

	updateGL();
	return true;
}

void MeshViewerCore::open_modelview_matrix_file(const char* _filename)
{
	assert(_filename);
	std::ifstream file(_filename);
	if (!file)
	{
		std::cerr << "Can't open file: \"" << _filename << "\"" << std::endl;
		return;
	}
	std::cout << "Loading " << _filename << "..." << std::endl;

	char buffer[256];
	unsigned int count = 0;
	float var = 0;

	GLdouble matrix[16];
	unsigned int num_values = 16;

	for (; !file.eof(); count++)
	{
		if (count >= num_values)
		{
			//std::cout << "Warning: Too many values." << std::endl;
			break;
		}

		file.getline(buffer, 256);

		int ret = sscanf(buffer, "%f\n", &var);
		if (ret < 1) continue;

		matrix[count] = var;
	}

	if (count < num_values)
	{
		std::cout << "Error: Number of values are not enough." << std::endl;
		return;
	}

	file.close();
	set_modelview_matrix(matrix);
	std::cout << "Done." << std::endl;
}

void MeshViewerCore::save_modelview_matrix_file(const char* _filename)
{
	assert(_filename);
	std::ofstream file(_filename);
	if(!file)
	{
		std::cerr << "Can't save file: \"" << _filename << "\"" << std::endl;
		return;
	}
	std::cout << "Saving " << _filename << "..." << std::endl;

	unsigned int num_values = 16;
	const GLdouble* matrix = modelview_matrix();
	for(unsigned int i = 0; i < num_values; i++)
	{
		file << matrix[i] << std::endl;
	}

	file.close();
	std::cout << "Done." << std::endl;
}

void MeshViewerCore::save_projection_matrix_file(const char* _filename)
{
	assert(_filename);
	std::ofstream file(_filename);
	if (!file)
	{
		std::cerr << "Can't save file: \"" << _filename << "\"" << std::endl;
		return;
	}
	std::cout << "Saving " << _filename << "..." << std::endl;

	unsigned int num_values = 16;
	const GLdouble* matrix = projection_matrix();
	for (unsigned int i = 0; i < num_values; i++)
	{
		file << matrix[i] << std::endl;
	}

	file.close();
	std::cout << "Done." << std::endl;
}

void MeshViewerCore::open_color_map_file(const char* _filename)
{
	assert(_filename);
	mesh_.load_vertex_color_map(_filename);
	updateGL();
}

void MeshViewerCore::save_color_map_file(const char* _filename)
{
	assert(_filename);
	mesh_.save_vertex_color_map(_filename);
	updateGL();
}

void MeshViewerCore::open_sample_point_file(const char* _filename)
{
	assert(_filename);
	bool ret = cuboid_structure_.load_sample_points(_filename);
	assert(ret);
	updateGL();
}

void MeshViewerCore::open_sample_point_label_file(const char* _filename)
{
	assert(_filename);
	bool ret = cuboid_structure_.load_sample_point_labels(_filename);
	assert(ret);
	updateGL();
}

void MeshViewerCore::open_mesh_face_label_file(const char* _filename)
{
	assert(_filename);
	bool ret = mesh_.load_face_label_simple(_filename);
	assert(ret);
	updateGL();
}

void MeshViewerCore::create_mesh_cuboids()
{
	cuboid_structure_.get_mesh_face_label_cuboids();
	updateGL();
}

void MeshViewerCore::apply_mesh_labels_to_cuboids()
{
	cuboid_structure_.apply_mesh_face_labels_to_cuboids();
	updateGL();
}

void MeshViewerCore::render_selection_mode()
{
	MyMesh::ConstVertexIter vit = mesh_.vertices_begin(), vEnd = mesh_.vertices_end();
	OpenMesh::VPropHandleT<MyMesh::Point> hPos;
	mesh_.get_property_handle(hPos, "v:points");

	glInitNames();
	glEnableClientState(GL_VERTEX_ARRAY);
	glVertexPointer(3, GL_DOUBLE, 0, &(mesh_.property(hPos, vit)));
	for(; vit != vEnd; ++vit)
	{
		glPushName(vit.handle().idx());
		glBegin(GL_POINTS);
		glArrayElement(vit.handle().idx());
		glEnd();
		glPopName();
	}
	glDisableClientState(GL_VERTEX_ARRAY);

	//// Face picking
	//Mesh::ConstFaceIter fIt(mesh().faces_begin()), fEnd(mesh().faces_end());
	//Mesh::ConstFaceVertexIter fvIt;
	//glInitNames();
	//for (; fIt!=fEnd; ++fIt)
	//{
	//	glPushName(fIt.handle().idx());
	//	glBegin(GL_TRIANGLES);
	//	glNormal3dv( &mesh().normal(fIt)[0] );

	//	fvIt = mesh().cfv_iter(fIt.handle()); 
	//	glVertex3dv( &mesh().point(fvIt)[0] );
	//	++fvIt;
	//	glVertex3dv( &mesh().point(fvIt)[0] );
	//	++fvIt;
	//	glVertex3dv( &mesh().point(fvIt)[0] );
	//	glEnd();
	//	glPopName();
	//}
}

VertexIndex MeshViewerCore::select_object(int x, int y)
{
	GLuint select_buffer[64];
	GLint hits, viewport[4];
	VertexIndex object_index = -1;

	glSelectBuffer(64, select_buffer);
	glRenderMode(GL_SELECT);

	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();
	glGetIntegerv(GL_VIEWPORT, viewport);

	// create 2x2 pixel picking region near cursor location 
	// 10 x 10 picking window
	gluPickMatrix((GLdouble) x, (GLdouble)(viewport[3] - y), 10, 10, viewport);
	gluPerspective(45.0, (GLdouble) width() / (GLdouble) height(),
		0.01*radius(), 100.0*radius());

	glMatrixMode( GL_MODELVIEW );
	glLoadMatrixd( modelview_matrix() );

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	render_selection_mode();

	swapBuffers();
	glMatrixMode(GL_PROJECTION);
	glPopMatrix();

	hits = glRenderMode(GL_RENDER);
	if(hits > 0)
	{
		GLuint* ptr = (GLuint*)select_buffer;
		object_index = ptr[3];
		std::cout << "Vertex: [" << object_index << "]" << std::endl;
	}
	return object_index;
}

void MeshViewerCore::mousePressEvent(const OpenMesh::Vec2f& _new_point_2d,
	bool _is_left_button, bool _is_mid_button, bool _is_right_button,
	bool _is_ctrl_pressed, bool _is_alt_pressed, bool _is_shift_pressed)
{
	if (_is_left_button && _is_shift_pressed )
	{
		// Selection mode.
		if (getDrawMode() == CUSTOM_VIEW || getDrawMode() == COLORED_RENDERING)
		{
			VertexIndex selected_vertex_index = -1;
			bool is_valid = false;

			// Set a seed vertex.
			if(selection_mode_ == PICK_SEED)
			{
				selected_vertex_index = select_object(_new_point_2d[0], _new_point_2d[1]);
				if(selected_vertex_index < mesh_.n_vertices())
				{
					mesh_.seed_vertex_index_ = selected_vertex_index;
					is_valid = true;
				}
			}

			// Set a query vertex.
			else if(selection_mode_ == PICK_QUERY)
			{
				selected_vertex_index = select_object(_new_point_2d[0], _new_point_2d[1]);
				if(selected_vertex_index < mesh_.n_vertices())
				{
					mesh_.query_vertex_index_ = selected_vertex_index;
					is_valid = true;
				}
			}

			/*
			// Add a feature vertex.
			else if (selection_mode_ == PICK_FEATURE)
			{
				QPoint newPoint2D = _event->pos();
				selected_vertex_index = select_object(newPoint2D.x(), newPoint2D.y());

				bool is_feature_vertex = false;
				for (VertexIndexArray::iterator it = mesh_.feature_vertex_indices_.begin();
					it != mesh_.feature_vertex_indices_.end(); it++)
				{
					if (*it == selected_vertex_index)
					{
						mesh_.feature_vertex_indices_.erase(it);
						is_feature_vertex = true;
						break;
					}
				}
				if (!is_feature_vertex)
				{
					mesh_.feature_vertex_indices_.push_back(selected_vertex_index);
					is_valid = true;
				}
			}
			*/
			
			if (is_valid)
			{
				mesh_.print_vertex_information(selected_vertex_index);
			}

			updateGL();
		}
	}
	else
	{
		GLViewerCore::mousePressEvent(_new_point_2d,
			_is_left_button, _is_mid_button, _is_right_button,
			_is_ctrl_pressed, _is_alt_pressed, _is_shift_pressed);
	}
}

void MeshViewerCore::keyPressEvent(const char _new_key,
	bool _is_ctrl_pressed, bool _is_alt_pressed, bool _is_shift_pressed)
{
	switch (_new_key)
	{
	//case Key_I:
	//	std::cout << std::endl;
	//	mesh_.show_picked_vertices_info();
	//	break;

	// FIXME.
	// Check whether the Ascii codes of special keys are correct.
	case 0x1B:	// Left
		if (cuboid_structure_.query_label_index_ - 1 < 0
			&& cuboid_structure_.num_labels() > 0)
		{
			cuboid_structure_.query_label_index_ = cuboid_structure_.num_labels() - 1;
		}
		else
		{
			--cuboid_structure_.query_label_index_;
		}
		std::cout << " - Label: " << cuboid_structure_.query_label_index_ << std::endl;
		updateGL();
		break;

	case 0x1A:	// Right
		// "mesh_.render_label_idx_ == mesh_.sample_point_num_labels_" draws all parts.
		if (cuboid_structure_.query_label_index_ + 1 > cuboid_structure_.num_labels())
		{
			cuboid_structure_.query_label_index_ = 0;
		}
		else
		{
			++cuboid_structure_.query_label_index_;
		}
		std::cout << " - Label: " << cuboid_structure_.query_label_index_ << std::endl;
		updateGL();
		break;

	/*
	// TEST.
	case 0x1B:	// Left
		if (_event->modifiers() & ControlModifier)
			run_test_rotate(-1.0);
		else if (_event->modifiers() & ShiftModifier)
			run_test_scale(-1.0, 0.0);
		else if (_event->modifiers() & AltModifier)
			run_test_translate(MyMesh::Normal(-1.0, 0.0, 0.0));
		break;

	case 0x1A:	// Right
		if (_event->modifiers() & ControlModifier)
			run_test_rotate(+1.0);
		else if (_event->modifiers() & ShiftModifier)
			run_test_scale(+1.0, 0.0);
		else if (_event->modifiers() & AltModifier)
			run_test_translate(MyMesh::Normal(+1.0, 0.0, 0.0));
		break;

	case 0x18:	// Up
		if (_event->modifiers() & ControlModifier);
		else if (_event->modifiers() & ShiftModifier)
			run_test_scale(0.0, +1.0);
		else if (_event->modifiers() & AltModifier)
			run_test_translate(MyMesh::Normal(0.0, +1.0, 0.0));
		break;

	case 0x19:	// Down
		if (_event->modifiers() & ControlModifier);
		else if (_event->modifiers() & ShiftModifier)
			run_test_scale(0.0, -1.0);
		else if (_event->modifiers() & AltModifier)
			run_test_translate(MyMesh::Normal(0.0, -1.0, 0.0));
		break;
	//
	*/

	case '1':
		std::cout << "Picking: [Seed]" << std::endl;
		selection_mode_ = PICK_SEED;
		break;

	case '2':
		std::cout << "Picking: [Query]" << std::endl;
		selection_mode_ = PICK_QUERY;
		break;

	case '3':
		std::cout << "Picking: [Feature]" << std::endl;
		selection_mode_ = PICK_FEATURE;
		break;

	default:
		GLViewerCore::keyPressEvent(_new_key,
			_is_ctrl_pressed, _is_alt_pressed, _is_shift_pressed);
		break;
	}
}

void draw_box(const std::array<MyMesh::Point, 8> &_corners)
{
	glBegin(GL_QUADS);

	// -X.
	glVertex3f(_corners[0][0], _corners[0][1], _corners[0][2]);
	glVertex3f(_corners[4][0], _corners[4][1], _corners[4][2]);
	glVertex3f(_corners[6][0], _corners[6][1], _corners[6][2]);
	glVertex3f(_corners[2][0], _corners[2][1], _corners[2][2]);

	// +X.
	glVertex3f(_corners[1][0], _corners[1][1], _corners[1][2]);
	glVertex3f(_corners[3][0], _corners[3][1], _corners[3][2]);
	glVertex3f(_corners[7][0], _corners[7][1], _corners[7][2]);
	glVertex3f(_corners[5][0], _corners[5][1], _corners[5][2]);
	
	// -Y.
	glVertex3f(_corners[0][0], _corners[0][1], _corners[0][2]);
	glVertex3f(_corners[1][0], _corners[1][1], _corners[1][2]);
	glVertex3f(_corners[5][0], _corners[5][1], _corners[5][2]);
	glVertex3f(_corners[4][0], _corners[4][1], _corners[4][2]);

	// +Y.
	glVertex3f(_corners[2][0], _corners[2][1], _corners[2][2]);
	glVertex3f(_corners[6][0], _corners[6][1], _corners[6][2]);
	glVertex3f(_corners[7][0], _corners[7][1], _corners[7][2]);
	glVertex3f(_corners[3][0], _corners[3][1], _corners[3][2]);

	// -Z.
	glVertex3f(_corners[0][0], _corners[0][1], _corners[0][2]);
	glVertex3f(_corners[2][0], _corners[2][1], _corners[2][2]);
	glVertex3f(_corners[3][0], _corners[3][1], _corners[3][2]);
	glVertex3f(_corners[1][0], _corners[1][1], _corners[1][2]);

	// +Z.
	glVertex3f(_corners[4][0], _corners[4][1], _corners[4][2]);
	glVertex3f(_corners[5][0], _corners[5][1], _corners[5][2]);
	glVertex3f(_corners[7][0], _corners[7][1], _corners[7][2]);
	glVertex3f(_corners[6][0], _corners[6][1], _corners[6][2]);

	glEnd();
}

void draw_arrow(MyMesh::Point _p1, MyMesh::Point _p2, GLdouble _arrow_size)
{
	MyMesh::Normal x_dir, y_dir, z_dir;
	z_dir = (_p2 - _p1).normalize();

	if (abs(dot(z_dir, MyMesh::Normal(0, 1, 0))) < 0.99)
	{
		x_dir = cross(MyMesh::Normal(0, 1, 0), z_dir);
		y_dir = cross(z_dir, x_dir);
	}
	else
	{
		y_dir = cross(z_dir, MyMesh::Normal(1, 0, 0));
		x_dir = cross(y_dir, z_dir);
	}

	GLdouble cylinder_matrix[16] = {
		x_dir[0], x_dir[1], x_dir[2], 0,
		y_dir[0], y_dir[1], y_dir[2], 0,
		z_dir[0], z_dir[1], z_dir[2], 0,
		_p1[0], _p1[1], _p1[2], 1
	};

	GLint polygon_mode[2];
	glGetIntegerv(GL_POLYGON_MODE, polygon_mode);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

	GLUquadric *cylinder = gluNewQuadric();
	gluQuadricDrawStyle(cylinder, GLU_FILL);
	gluQuadricNormals(cylinder, GLU_SMOOTH);
	gluQuadricOrientation(cylinder, GLU_OUTSIDE);

	glPushMatrix();
	glMultMatrixd(cylinder_matrix);

	GLdouble arrow_length = (_p2 - _p1).length() - (2 * _arrow_size);

	glTranslated(0, 0, _arrow_size);
	gluCylinder(cylinder, 0.3*_arrow_size, 0.3*_arrow_size, arrow_length, 20, 20);

	glRotated(180, 1, 0, 0);

	glRotated(-180, 1, 0, 0);
	glTranslated(0, 0, arrow_length);
	glutSolidCone(_arrow_size, _arrow_size, 20, 20);

	glPopMatrix();
	gluDeleteQuadric(cylinder);

	glPolygonMode(GL_FRONT, polygon_mode[0]);
	glPolygonMode(GL_BACK, polygon_mode[1]);
}

void MeshViewerCore::draw_scene(const std::string& _draw_mode)
{
	MeshViewerCoreT<MyMesh>::draw_scene(_draw_mode);

	if (_draw_mode == CUSTOM_VIEW)
	{
		glDisable(GL_LIGHTING);
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		draw_openmesh(_draw_mode);
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	}
	else if (_draw_mode == COLORED_RENDERING)
	{
		glEnable(GL_LIGHTING);
		glShadeModel(GL_SMOOTH);
		draw_openmesh(_draw_mode);
		setDefaultMaterial();
	}
	else if (_draw_mode == FACE_INDEX_RENDERING)
	{
		glDisable(GL_LIGHTING);
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		draw_openmesh(_draw_mode);
	}
}

void MeshViewerCore::draw_openmesh(const std::string& _drawmode)
{
	MeshViewerCoreT<MyMesh>::draw_openmesh(_drawmode);

	if (_drawmode == CUSTOM_VIEW) // -------------------------------------------
	{
		Mesh::ConstFaceIter f_it(mesh_.faces_begin()), f_end(mesh_.faces_end());
		Mesh::ConstFaceVertexIter fv_it;

		/*
		// Draw faces.
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		glDepthRange(0.01, 1.0);

		if(mesh_.mesh_coloring_option_ == mesh_.VERTEX_COLOR)
		{
			glBegin(GL_TRIANGLES);
			for (; f_it != f_end; ++f_it)
			{
				fv_it = mesh_.cfv_iter(f_it.handle()); 
				glColor(fv_it.handle());
				glVertex3dv( &mesh_.point(fv_it)[0] );
				++fv_it;
				glColor(fv_it.handle());
				glVertex3dv( &mesh_.point(fv_it)[0] );
				++fv_it;
				glColor(fv_it.handle());
				glVertex3dv( &mesh_.point(fv_it)[0] );
			}
			glEnd();
		}
		else if(mesh_.mesh_coloring_option_ == mesh_.FACE_COLOR)
		{
			glBegin(GL_TRIANGLES);
			for (; f_it != f_end; ++f_it)
			{
				glColor(f_it.handle());
				fv_it = mesh_.cfv_iter(f_it.handle()); 
				glVertex3dv( &mesh_.point(fv_it)[0] );
				++fv_it;
				glVertex3dv( &mesh_.point(fv_it)[0] );
				++fv_it;
				glVertex3dv( &mesh_.point(fv_it)[0] );
			}
			glEnd();
		}
		*/

		// Draw edges (Black).
		glPolygonMode( GL_FRONT_AND_BACK, GL_LINE);
		gray_color();
		glLineWidth(0.1f);
		glDepthRange( 0.0, 1.0 );
		draw_openmesh( "Wireframe" );

		/*
		// Draw a seed (Red).
		if( mesh_.seed_vertex_index_ < mesh_.n_vertices() )
		{
			GLdouble *point = &mesh_.point(mesh_.vertex_handle(mesh_.seed_vertex_index_))[0];
			glPushMatrix();
			glTranslatef(point[0], point[1], point[2]);
			red_color();
			glutSolidSphere((mesh_.get_object_diameter() * 0.01) * point_size_, 20, 20);
			glPopMatrix();
		}

		// Draw a query vertex (Blue).
		if( mesh_.query_vertex_index_ < mesh_.n_vertices() )
		{
			GLdouble *point = &mesh_.point(mesh_.vertex_handle(mesh_.query_vertex_index_))[0];
			glPushMatrix();
			glTranslatef(point[0], point[1], point[2]);
			blue_color();
			glutSolidSphere((mesh_.get_object_diameter() * 0.01) * point_size_, 20, 20);
			glPopMatrix();
		}
		*/

		/*
		// Draw occlusion test points.
		if (draw_occlusion_test_points_)
		{
			for (std::vector< std::pair<Mesh::Point, Real> >::iterator it = occlusion_test_points_.begin();
				it != occlusion_test_points_.end(); ++it)
			{
				GLdouble *point = &(*it).first[0];
				Real radius = (*it).second;
				radius = radius * (mesh_.get_object_diameter() * 0.0075) * point_size_;
				if (radius > 0)
				{
					glPushMatrix();
					glTranslatef(point[0], point[1], point[2]);
					red_color();

					glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
					glutSolidSphere(radius, 4, 4);
					glPopMatrix();
				}
			}
		}
		*/

		bool draw_all_labels = (cuboid_structure_.query_label_index_ == cuboid_structure_.num_labels());

		// Draw sample points (Black).
		if (draw_all_labels)
		{
			for (std::vector<MeshSamplePoint *>::iterator it = cuboid_structure_.sample_points_.begin();
				it != cuboid_structure_.sample_points_.end(); it++)
			{
				GLdouble *point = &(*it)->point_[0];

				Real radius = 1.0;
				if (cuboid_structure_.query_label_index_ < cuboid_structure_.num_labels())
				{
					assert((*it)->label_index_confidence_.size() == cuboid_structure_.num_labels());
					radius = (*it)->label_index_confidence_[cuboid_structure_.query_label_index_];
				}

				radius = radius * (mesh_.get_object_diameter() * 0.01) * point_size_;
				if (radius > 0)
				{
					glPushMatrix();
					glTranslatef(point[0], point[1], point[2]);
					black_color();

					glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
					glutSolidSphere(radius, 20, 20);
					glPopMatrix();
				}
			}
		}

		
		// For each label.
		for (LabelIndex label_index = 0; label_index < cuboid_structure_.label_cuboids_.size(); ++label_index)
		{
			// "cuboid_structure_.query_label_index_ == cuboid_structure_.label_cuboids_.size()" draws all parts.
			if (cuboid_structure_.query_label_index_ < cuboid_structure_.label_cuboids_.size())
			{
				if (label_index != cuboid_structure_.query_label_index_)
					continue;
			}

			Label label = cuboid_structure_.get_label(label_index);
			MyMesh::Color label_color = MyMesh::get_label_color(label);
			Real r = label_color[0] / 255.0;
			Real g = label_color[1] / 255.0;
			Real b = label_color[2] / 255.0;

			// For each cuboid.
			for (std::vector<MeshCuboid *>::iterator it = cuboid_structure_.label_cuboids_[label_index].begin();
				it != cuboid_structure_.label_cuboids_[label_index].end(); ++it)
			{
				MeshCuboid *cuboid = (*it);
				assert(cuboid);

				//// Draw sample points (Black).
				//if (!draw_all_labels)
				//{
				//	for (unsigned int point_index = 0; point_index < cuboid->num_sample_points(); ++point_index)
				//	{
				//		const MeshSamplePoint *sample_point = cuboid->get_sample_point(point_index);
				//		const GLdouble *point = &(sample_point->point_[0]);

				//		Real radius = (mesh_.get_object_diameter() * 0.005) * point_size_;
				//		if (radius > 0)
				//		{
				//			glPushMatrix();
				//			glTranslatef(point[0], point[1], point[2]);
				//			black_color();

				//			glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
				//			glutSolidSphere(radius, 20, 20);
				//			glPopMatrix();
				//		}
				//	}
				//}

				// Draw cuboid surface points (Red).
				for (unsigned int point_index = 0; point_index < cuboid->num_cuboid_surface_points(); ++point_index)
				{
					const MeshCuboidSurfacePoint *cuboid_surface_point =
						cuboid->get_cuboid_surface_point(point_index);
					MyMesh::Point point = cuboid_surface_point->point_;
					MyMesh::Normal normal = cuboid_surface_point->normal_;
					Real visibility = cuboid_surface_point->visibility_;
					Real radius = (visibility) * (mesh_.get_object_diameter() * 0.002) * point_size_;

					glPushMatrix();
					glTranslatef(point[0], point[1], point[2]);
					red_color();

					glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
					glutSolidSphere(radius, 20, 20);
					glPopMatrix();

					glBegin(GL_LINES);
					red_color();
					glLineWidth(4.0f);
					glVertex(point);
					glVertex(point + (2 * radius) * normal);
					glEnd();
				}

				// Draw point correspondences (sample point -> cuboid surface point).
				for (unsigned int sample_point_index = 0;
					sample_point_index < cuboid->num_sample_points(); ++sample_point_index)
				{
					int cuboid_surface_point_index =
						cuboid->get_sample_to_cuboid_surface_correspondences(sample_point_index);

					if (cuboid_surface_point_index < 0)
						continue;

					MyMesh::Point sample_point = cuboid->get_sample_point(sample_point_index)->point_;
					MyMesh::Point cuboid_surface_point = cuboid->get_cuboid_surface_point(
						cuboid_surface_point_index)->point_;

					glBegin(GL_LINES);
					black_color();
					glLineWidth(1.0f);
					glVertex(sample_point);
					glVertex(cuboid_surface_point);
					glEnd();
				}

				// Draw point correspondences (cuboid surface point -> sample point).
				for (unsigned int cuboid_surface_point_index = 0;
					cuboid_surface_point_index < cuboid->num_cuboid_surface_points();
					++cuboid_surface_point_index)
				{
					Real visibility = cuboid->get_cuboid_surface_point(
						cuboid_surface_point_index)->visibility_;
					if (visibility < 1.0)
						continue;

					int sample_point_index =
						cuboid->get_cuboid_surface_to_sample_correspondence(cuboid_surface_point_index);
					if (sample_point_index < 0)
						continue;

					MyMesh::Point sample_point = cuboid->get_sample_point(sample_point_index)->point_;
					MyMesh::Point cuboid_surface_point = cuboid->get_cuboid_surface_point(
						cuboid_surface_point_index)->point_;

					glBegin(GL_LINES);
					red_color();
					glLineWidth(1.0f);
					glVertex(sample_point);
					glVertex(cuboid_surface_point);
					glEnd();
				}

				// Draw cuboid axes.
				if (draw_cuboid_axes_)
				{
					glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
					red_color();
					draw_arrow(cuboid->get_bbox_center(), cuboid->get_bbox_center() +
						(0.5 * cuboid->get_bbox_size()[0] + 0.04 * mesh_.get_object_diameter()) * cuboid->get_bbox_axis(0),
						0.02 * mesh_.get_object_diameter());
					green_color();
					draw_arrow(cuboid->get_bbox_center(), cuboid->get_bbox_center() +
						(0.5 * cuboid->get_bbox_size()[1] + 0.04 * mesh_.get_object_diameter()) * cuboid->get_bbox_axis(1),
						0.02 * mesh_.get_object_diameter());
					blue_color();
					draw_arrow(cuboid->get_bbox_center(), cuboid->get_bbox_center() +
						(0.5 * cuboid->get_bbox_size()[2] + 0.04 * mesh_.get_object_diameter()) * cuboid->get_bbox_axis(2),
						0.02 * mesh_.get_object_diameter());
				}

				// Draw cuboid.
				glColor4f(r, g, b, 0.5f);
				glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
				glEnable(GL_BLEND);
				glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
				//draw_box(cuboid->get_bbox_corners());
				cuboid->draw_cuboid();
				glDisable(GL_BLEND);

				glColor4f(0.0f, 0.0f, 0.0f, 1.0f);
				glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
				glLineWidth(8.0f);
				//draw_box(cuboid->get_bbox_corners());
				cuboid->draw_cuboid();
				glLineWidth(1.0f);
			}
		}

		// Draw viewing direction.
		red_color();
		draw_arrow(view_point_ + (1.0 * mesh_.get_object_diameter()) * view_direction_,
			view_point_ + (1.25 * mesh_.get_object_diameter()) * view_direction_,
			0.02 * mesh_.get_object_diameter());
	}
	else if (_drawmode == COLORED_RENDERING) // -------------------------------------------
	{
		Mesh::ConstFaceIter f_it(mesh_.faces_begin()), f_end(mesh_.faces_end());
		Mesh::ConstFaceVertexIter fv_it;

		// Draw faces.
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		glDepthRange(0.01, 1.0);

		mesh_.update_normals();

		if(mesh_.mesh_coloring_option_ == mesh_.VERTEX_COLOR)
		{
			glBegin(GL_TRIANGLES);
			for (; f_it != f_end; ++f_it)
			{
				fv_it = mesh_.cfv_iter(f_it.handle()); 
				glMaterial(fv_it.handle());
				glNormal3dv( &mesh_.normal(fv_it)[0] );
				glVertex3dv( &mesh_.point(fv_it)[0] );
				++fv_it;
				glMaterial(fv_it.handle());
				glNormal3dv( &mesh_.normal(fv_it)[0] );
				glVertex3dv( &mesh_.point(fv_it)[0] );
				++fv_it;
				glMaterial(fv_it.handle());
				glNormal3dv( &mesh_.normal(fv_it)[0] );
				glVertex3dv( &mesh_.point(fv_it)[0] );
			}
			glEnd();

		}
		else if(mesh_.mesh_coloring_option_ == mesh_.FACE_COLOR)
		{
			glBegin(GL_TRIANGLES);
			for (; f_it != f_end; ++f_it)
			{
				glMaterial(f_it.handle());
				fv_it = mesh_.cfv_iter(f_it.handle()); 
				glNormal3dv( &mesh_.normal(fv_it)[0] );
				glVertex3dv( &mesh_.point(fv_it)[0] );
				++fv_it;
				glNormal3dv( &mesh_.normal(fv_it)[0] );
				glVertex3dv( &mesh_.point(fv_it)[0] );
				++fv_it;
				glNormal3dv( &mesh_.normal(fv_it)[0] );
				glVertex3dv( &mesh_.point(fv_it)[0] );
			}
			glEnd();
		}

		// Draw a seed (Black).
		if( mesh_.seed_vertex_index_ < mesh_.n_vertices() )
		{
			GLdouble *point = &mesh_.point(mesh_.vertex_handle(mesh_.seed_vertex_index_))[0];
			glPushMatrix();
			glTranslatef(point[0], point[1], point[2]);

			// Black
			OpenMesh::Vec4f m(  0.0f / 256.0f,  32.0f / 256.0f,  96.0f / 256.0f, 1.0f);
			glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, &m[0]); 

			glutSolidSphere((mesh_.get_object_diameter() * 0.01) * point_size_, 20, 20);
			glPopMatrix();
		}

		// Draw a query vertex (Red).
		if( mesh_.query_vertex_index_ < mesh_.n_vertices() )
		{
			GLdouble *point = &mesh_.point(mesh_.vertex_handle(mesh_.query_vertex_index_))[0];
			glPushMatrix();
			glTranslatef(point[0], point[1], point[2]);

			// Red
			OpenMesh::Vec4f m(192.0f / 256.0f,   0.0f / 256.0f,   0.0f / 256.0f, 1.0f);
			glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, &m[0]); 

			glutSolidSphere((mesh_.get_object_diameter() * 0.01) * point_size_, 20, 20);
			glPopMatrix();
		}

		// Draw feature points (Black).
		//if (mesh_.render_label_idx_ >= 0)
		//{
		//	assert(mesh_.render_label_idx_ < mesh_.num_labels_);

		//	for (std::vector<MeshSamplePoint *>::iterator it = mesh_.sample_points_.begin();
		//		it != mesh_.sample_points_.end(); it++)
		//	{
		//		GLdouble *point = &(*it)->pos_[0];
		//		dReal radius = (*it)->label_conf_[mesh_.render_label_idx_] * mesh_.avg_edge_length_ * pointSize;
		//		if (radius > 0)
		//		{
		//			glPushMatrix();
		//			glTranslatef(point[0], point[1], point[2]);

		//			// Black
		//			OpenMesh::Vec4f m(0.0f / 256.0f, 32.0f / 256.0f, 96.0f / 256.0f, 1.0f);
		//			glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, &m[0]);

		//			glutSolidSphere(radius, 20, 20);
		//			glPopMatrix();
		//		}
		//	}
		//}
	}
	else if (_drawmode == "Solid Smooth") // -------------------------------------------
	{
		srand(COLOR_RANDOM_SEED);

		// Draw feature points (Random).
		//if (mesh_.render_label_idx_ >= 0)
		//{
		//	assert(mesh_.render_label_idx_ < mesh_.num_labels_);

		//	for (std::vector<MeshSamplePoint *>::iterator it = mesh_.sample_points_.begin();
		//		it != mesh_.sample_points_.end(); it++)
		//	{
		//		GLdouble *point = &(*it)->pos_[0];
		//		dReal radius = (*it)->label_conf_[mesh_.render_label_idx_] * mesh_.avg_edge_length_ * pointSize;
		//		if (radius > 0)
		//		{
		//			glPushMatrix();
		//			glTranslatef(point[0], point[1], point[2]);

		//			// Random Color
		//			float r = rand() / (float)RAND_MAX;
		//			float g = rand() / (float)RAND_MAX;
		//			float b = rand() / (float)RAND_MAX;
		//			OpenMesh::Vec4f m(r, g, b, 1.0f);

		//			glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, &m[0]);
		//			glutSolidSphere(radius, 20, 20);
		//			setDefaultMaterial();

		//			glPopMatrix();
		//		}
		//	}
		//}

#ifdef DRAW_FEATURE_PAIR_LINE
		// Draw feature pair lines.
		if(mesh_.feature_idxs_.size() % 2 == 0)
		{
			for(unsigned int i = 0; i < mesh_.feature_idxs_.size(); i = i + 2)
			{
				VertexIndex vid_A = mesh_.feature_idxs_[i];
				VertexIndex vid_B = mesh_.feature_idxs_[i + 1];

				GLdouble *point_A = &(mesh_.point(mesh_.vertex_handle(vid_A))[0]);
				GLdouble *point_B = &(mesh_.point(mesh_.vertex_handle(vid_B))[0]);

				GLdouble coord[3];
				glMatrixMode( GL_MODELVIEW );
				glLoadIdentity();

				glBegin(GL_LINES);
				getTransformedCoord(modelview_matrix(), point_A, coord);
				glVertex3dv( coord );
				getTransformedCoord(modelview_matrix(), point_B, coord);
				glVertex3dv( coord );
				glEnd();
			}
		}
#endif
	}
	else if (_drawmode == FACE_INDEX_RENDERING)
	{
		// NOTE:
		// (255, 255, 255) is background color.
		if (mesh_.n_faces() >= static_cast<unsigned int>(std::pow(2, 8 * 3)))
		{
			std::cout << "Warning: Too many faces." << std::endl;
		}
		else
		{
			Mesh::ConstFaceIter f_it(mesh_.faces_begin()), f_end(mesh_.faces_end());
			Mesh::ConstFaceVertexIter fv_it;

			glBegin(GL_TRIANGLES);
			for (; f_it != f_end; ++f_it)
			{
				FaceIndex face_index = f_it->idx();
				GLubyte r = static_cast<GLubyte>(face_index % 256);
				GLubyte g = static_cast<GLubyte>((face_index >> 8) % 256);
				GLubyte b = static_cast<GLubyte>((face_index >> 16) % 256);
				glColor3ub(r, g, b);

				fv_it = mesh_.cfv_iter(f_it.handle());
				glVertex3dv(&mesh_.point(fv_it)[0]);
				++fv_it;
				glVertex3dv(&mesh_.point(fv_it)[0]);
				++fv_it;
				glVertex3dv(&mesh_.point(fv_it)[0]);
			}
			glEnd();
		}
	}
}

void MeshViewerCore::remove_occluded_points()
{
	std::string curr_draw_mode = getDrawMode();
	setDrawMode(FACE_INDEX_RENDERING);

	size_t w(width()), h(height());
	GLenum buffer(GL_BACK);

	std::vector<GLubyte> fbuffer(3 * w*h);
	std::set<FaceIndex> visible_face_indices;

	//qApp->processEvents();
	makeCurrent();
	updateGL();
	glFinish();

	glReadBuffer(buffer);
	glPixelStorei(GL_PACK_ALIGNMENT, 1);
	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

	updateGL();
	glReadPixels(0, 0, w, h, GL_RGB, GL_UNSIGNED_BYTE, &fbuffer[0]);

	unsigned int x, y, offset;

	for (y = 0; y < h; ++y) {
		for (x = 0; x < w; ++x) {
			offset = 3 * (y*w + x);
			unsigned int r = static_cast<unsigned int>(fbuffer[offset]);
			unsigned int g = static_cast<unsigned int>(fbuffer[offset + 1]);
			unsigned int b = static_cast<unsigned int>(fbuffer[offset + 2]);

			// NOTE:
			// (255, 255, 255) is background color.
			if (r >= 255 && g >= 255 && b >= 255)
				continue;

			FaceIndex face_index = r + 256 * g + 65536 * b;
			assert(face_index < mesh_.n_faces());
			visible_face_indices.insert(face_index);
		}
	}

	MyMesh::ConstFaceIter f_it, f_start(mesh_.faces_begin()), f_end(mesh_.faces_end());
	for (f_it = f_start; f_it != f_end; ++f_it)
		mesh_.set_color(f_it.handle(), MyMesh::Color(0, 255, 255));

	for (std::set<FaceIndex>::iterator it = visible_face_indices.begin();
		it != visible_face_indices.end(); ++it)
	{
		FaceIndex face_index = *it;
		MyMesh::FaceHandle fh = mesh_.face_handle(face_index);
		mesh_.set_color(fh, MyMesh::Color(255, 0, 0));
	}

	setDrawMode(curr_draw_mode);

	// Remove occluded sample points.
	cuboid_structure_.remove_occluded_sample_points(visible_face_indices);

	updateGL();
}

bool print_faces_info(const MyMesh &_mesh, const char *_filename)
{
	std::ofstream file(_filename);
	if( !file.is_open() )
	{
		std::cerr << "Error: [print_faces_info] failed to open file: " << _filename << std::endl;
		return false;
	}

	std::cout << "Saving face information... """ << _filename << """" << std::endl;
	std::cout << "Format:" << std::endl;
	std::cout << "Number of faces" << std::endl;
	std::cout << "fid,vid1,vid2,vid3," << std::endl;
	std::cout << "..." << std::endl;


	unsigned int num_faces = _mesh.n_faces();
	file << num_faces << std::endl;

	for(unsigned int fid = 0; fid < num_faces; fid++)
	{
		MyMesh::FaceHandle fh = _mesh.face_handle(fid);
		file << fid << ",";

		for(MyMesh::ConstFaceVertexIter fv_it = _mesh.cfv_iter(fh); fv_it; ++fv_it)
		{
			unsigned int vid = fv_it.handle().idx();
			file << vid << ",";
		}

		file << std::endl;
	}


	std::cout << "Done." << std::endl;
	std::cout << std::endl;
	
	file.close();
	return true;
}

bool print_vertices_info(const MyMesh &_mesh, const char *_filename)
{
	std::ofstream file(_filename);
	if( !file.is_open() )
	{
		std::cerr << "Error: Failed to open file: " << _filename << std::endl;
		return false;
	}

	std::cout << "Saving vertex information... """ << _filename << """" << std::endl;
	std::cout << "Format:" << std::endl;
	std::cout << "Number of vertices" << std::endl;
	std::cout << "vid,x,y,z," << std::endl;
	std::cout << "..." << std::endl;


	unsigned int num_vertices = _mesh.n_vertices();
	file << num_vertices << std::endl;

	for(unsigned int vid = 0; vid < num_vertices; vid++)
	{
		MyMesh::VertexHandle vh = _mesh.vertex_handle(vid);
		file << vid << ",";

		for(unsigned int i = 0; i < 3; i++)
			file << _mesh.point(vh)[i] << ",";

		file << std::endl;
	}


	std::cout << "Done." << std::endl;
	std::cout << std::endl;

	file.close();
	return true;
}

bool print_neighbors_info(const MyMesh &_mesh, const char *_filename)
{
	std::ofstream file(_filename);
	if( !file.is_open() )
	{
		std::cerr << "Error: Failed to open file: " << _filename << std::endl;
		return false;
	}

	std::cout << "Saving connectivity information... """ << _filename << """" << std::endl;
	std::cout << "Format:" << std::endl;
	std::cout << "Number of vertices" << std::endl;
	std::cout << "vid,Number of neighbors,n_vid1,n_vid2,..." << std::endl;
	std::cout << "..." << std::endl;


	unsigned int num_vertices = _mesh.n_vertices();
	file << num_vertices << std::endl;

	for(unsigned int vid = 0; vid < num_vertices; vid++)
	{
		MyMesh::VertexHandle vh = _mesh.vertex_handle(vid);
		file << vid << ",";

		unsigned int num_neighbors = 0;
		for(MyMesh::ConstVertexVertexIter vv_it = _mesh.cvv_iter(vh); vv_it; ++vv_it)
			num_neighbors++;

		file << num_neighbors << ",";

		for(MyMesh::ConstVertexVertexIter vv_it = _mesh.cvv_iter(vh); vv_it; ++vv_it)
		{
			unsigned int n_vid = vv_it.handle().idx();
			file << n_vid << ",";
		}

		file << std::endl;
	}


	std::cout << "Done." << std::endl;
	std::cout << std::endl;

	file.close();
	return true;
}

void MeshViewerCore::print_mesh_info()
{
	// Nov. 2012, Minhyuk Sung
	print_faces_info(mesh_, FACE_INFO_FILENAME);
	print_vertices_info(mesh_, VERTEX_INFO_FILENAME);
	print_neighbors_info(mesh_, NEIGHBOR_INFO_FILENAME);
}
