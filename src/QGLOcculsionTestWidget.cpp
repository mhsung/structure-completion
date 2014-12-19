#ifdef _MSC_VER
#  pragma warning(disable: 4267 4311 4305)
#endif

#include <assert.h>
#include <algorithm>
#include <iostream>
#include <limits>
#include <fstream>
#include <sstream>

#include "Angel.h"
#ifdef ARCH_DARWIN
#  include <glut.h>
#else
#  include <GL/glew.h>
#  include <GL/freeglut.h>
#endif

#include <qnamespace.h>
#include <qapplication.h>
#include <qtimer.h>
#include <qimage.h>

#include "QGLOcculsionTestWidget.h"


const char *sample_point_filename =
"D:/Development/stanford-projects/shape2pose/data/2_analysis/coseg_chairs/points/even1000/1.pts";
const char *color_vert_shader_filename = "glsl/OcclusionRendering.vert";
const char *color_frag_shader_filename = "glsl/OcclusionRendering.frag";


QGLOcculsionTestWidget::QGLOcculsionTestWidget(QWidget* _parent)
	: QGLWidget(_parent)
{
	initialize();
}

QGLOcculsionTestWidget::
QGLOcculsionTestWidget(QGLFormat& _fmt, QWidget* _parent)
: QGLWidget(_fmt, _parent)
{
	initialize();
}

void QGLOcculsionTestWidget::initialize(void)
{
	resolution_unit_ = DEFAULT_RESOLUTION_UNIT;
	resolution_unit_square_ = resolution_unit_ * resolution_unit_;
	resolution_unit_cubic_ = resolution_unit_ * resolution_unit_ * resolution_unit_;

	window_width_ = resolution_unit_cubic_;
	window_height_ = resolution_unit_cubic_;
	setFixedSize(window_width_, window_height_);

	shader_program_ = 0;

	vertex_array_object_ = 0;

	position_buffer_ = 0;
	position_attribute_ = 0;

	texcoord_buffer_ = 0;
	texcoord_attribute_ = 0;

	glsl_uniform_modelview_matrix_ = 0;
	modelview_matrix_ = Eigen::Matrix<float, 4, 4, Eigen::ColMajor>::Identity();

	glsl_uniform_sample_points_ = 0;
	sample_points_.clear();
	sample_points_start_index_ = 0;

	glsl_uniform_num_sample_points_ = 0;
	num_sample_points_ = 0;

	glsl_uniform_radius_ = 0;
	radius_ = DAFAULT_SAMPLE_POINT_RADIUS;

	glsl_uniform_bbox_center_ = 0;
	bbox_center_ = Eigen::Vector3f::Zero();

	glsl_uniform_bbox_size_ = 0;
	bbox_size_ = Eigen::Vector3f::Ones();

	num_vertices_ = 0;
	num_indices_ = 0;
	positions_mat_.resize(num_vertices_, 3);
	texcoords_mat_.resize(num_vertices_, 2);
	indices_vec_.resize(num_indices_);

	angle_ = 0;


	// qt stuff
	setAttribute(Qt::WA_NoSystemBackground, true);
	setFocusPolicy(Qt::StrongFocus);
	setAcceptDrops(false);
	setCursor(Qt::WaitCursor);
}

QGLOcculsionTestWidget::~QGLOcculsionTestWidget()
{
}

bool QGLOcculsionTestWidget::load_sample_points(const char *_filename, bool _verbose)
{
	std::ifstream file(_filename);
	if (!file)
	{
		std::cerr << "Can't open file: \"" << _filename << "\"" << std::endl;
		return false;
	}

	if (_verbose)
		std::cout << "Loading " << _filename << "..." << std::endl;

	std::string buffer;


	sample_points_.clear();

	while (!file.eof())
	{
		std::getline(file, buffer);
		if (buffer == "")
			continue;

		std::stringstream strstr(buffer);
		std::string token;
		std::getline(strstr, token, ' ');
		int corr_fid = atoi(token.c_str());
		assert(corr_fid >= 0);

		if (strstr.eof())
			continue;

		Eigen::Vector3f bicentric_point;
		std::getline(strstr, token, ' ');
		bicentric_point[0] = atof(token.c_str());
		std::getline(strstr, token, ' ');
		bicentric_point[1] = atof(token.c_str());
		std::getline(strstr, token, ' ');
		bicentric_point[2] = atof(token.c_str());

		Eigen::Vector3f point;
		std::getline(strstr, token, ' ');
		point[0] = atof(token.c_str());
		std::getline(strstr, token, ' ');
		point[1] = atof(token.c_str());
		std::getline(strstr, token, ' ');
		point[2] = atof(token.c_str());

		sample_points_.push_back(point);
	}

	num_sample_points_ = sample_points_.size();
	if (num_sample_points_ <= 0)
		return false;

	file.close();

	sample_points_start_index_ = 0;
	get_bounding_box();

	std::cout << "Done." << std::endl;
	return true;
}

void QGLOcculsionTestWidget::get_bounding_box()
{
	Eigen::Vector3f bbox_min_ = Eigen::Vector3f::Constant(std::numeric_limits<float>::max());
	Eigen::Vector3f bbox_max_ = Eigen::Vector3f::Constant(std::numeric_limits<float>::lowest());

	for (int point_index = 0; point_index < num_sample_points_; ++point_index)
	{
		for (unsigned int i = 0; i < 3; ++i)
		{
			bbox_min_[i] = std::min(sample_points_[point_index][i], bbox_min_[i]);
			bbox_max_[i] = std::max(sample_points_[point_index][i], bbox_max_[i]);
		}
	}

	bbox_center_ = 0.5 * (bbox_min_ + bbox_max_);
	bbox_size_ = bbox_max_ - bbox_min_;
}

void QGLOcculsionTestWidget::set_sample_points(const std::vector<Eigen::Vector3f> &_sample_points)
{
	assert(!_sample_points.empty());

	sample_points_.clear();
	sample_points_ = _sample_points;

	sample_points_start_index_ = 0;
	num_sample_points_ = sample_points_.size();

	sample_points_start_index_ = 0;
	get_bounding_box();
}

void QGLOcculsionTestWidget::set_modelview_matrix(const GLdouble _modelview_matrix[16])
{
	// NOTE:
	// '_modelview_matrix' is a column-major array.
	assert(_modelview_matrix);
	for (unsigned int col = 0; col < 4; ++col)
		for (unsigned int row = 0; row < 4; ++row)
			modelview_matrix_(row, col) = static_cast<float>(_modelview_matrix[4 * col + row]);
}

void QGLOcculsionTestWidget::generate_vertices()
{
	// NOTE:
	// Each square indicates one slice of voxels.
	// A multiple slices with different depth values will be rendered simultaneously.
	// The size of each axis of square is 'resolution_unit_square',
	// and the number of slices (number of depth values) is also 'resolution_unit_square'.
	// All slices are rendered in a 'resolution_unit' by 'resolution_unit' grid.

	// The position value [0, 1] which should be mapped to
	// [-0.55 * bbox_size_ + bbox_center_, +0.55 * bbox_size_ + bbox_center_].
	Eigen::Matrix<float, 4, 3, Eigen::RowMajor> square_position_mat;
	square_position_mat <<
		0.0, 0.0, 0.0,
		1.0, 0.0, 0.0,
		1.0, 1.0, 0.0,
		0.0, 1.0, 0.0;

	Eigen::Matrix<float, 4, 2, Eigen::RowMajor> square_texcoord_mat;
	square_texcoord_mat <<
		0.0, 0.0,
		resolution_unit_square_, 0.0,
		resolution_unit_square_, resolution_unit_square_,
		0.0, resolution_unit_square_;

	Eigen::Matrix<unsigned int, 6, 1> square_index_vec;
	square_index_vec <<
		0, 1, 2,
		2, 3, 0;

	num_vertices_ = 4 * resolution_unit_square_;
	num_indices_ = 6 * resolution_unit_square_;

	positions_mat_.resize(num_vertices_, 3);
	texcoords_mat_.resize(num_vertices_, 2);
	indices_vec_.resize(num_indices_);

	// All slices are rendered in a 'resolution_unit' by 'resolution_unit' grid.
	for (int row = 0; row < resolution_unit_; ++row)
	{
		for (int col = 0; col < resolution_unit_; ++col)
		{
			int depth_index = row * resolution_unit_ + col;
			assert(depth_index <= resolution_unit_square_ - 1);

			// Each voxel slice has a unique depth value.
			GLfloat depth = (GLfloat)depth_index / (resolution_unit_square_ - 1);
			assert(depth >= 0.0);
			assert(depth <= 1.0);

			for (unsigned int i = 0; i < 4; ++i)
			{
				positions_mat_.row(4 * depth_index + i) = square_position_mat.row(i);
				// Set depth.
				positions_mat_.row(4 * depth_index + i)[2] = depth;

				texcoords_mat_.row(4 * depth_index + i) = square_texcoord_mat.row(i)
					+ Eigen::RowVector2f(col * resolution_unit_square_, row * resolution_unit_square_);
				// Arrange in the grid.
				for (int k = 0; k < 2; ++k)
					texcoords_mat_.row(4 * depth_index + i)[k] = texcoords_mat_.row(4 * depth_index + i)[k]
						/ resolution_unit_cubic_ * 2.0 - 1.0;
			}

			for (unsigned int i = 0; i < 6; ++i)
				indices_vec_(6 * depth_index + i) = square_index_vec(i) + (4 * depth_index);
		}
	}
}

void QGLOcculsionTestWidget::load_buffer_data(){
	glGenVertexArrays(1, &vertex_array_object_);
	glBindVertexArray(vertex_array_object_);

	glGenBuffers(1, &position_buffer_);
	glBindBuffer(GL_ARRAY_BUFFER, position_buffer_);
	glBufferData(GL_ARRAY_BUFFER, num_vertices_ * 3 * sizeof(float), positions_mat_.data(), GL_STATIC_DRAW);

	glEnableVertexAttribArray(position_attribute_);
	glVertexAttribPointer(position_attribute_, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (const GLvoid *)0);

	glGenBuffers(1, &texcoord_buffer_);
	glBindBuffer(GL_ARRAY_BUFFER, texcoord_buffer_);
	glBufferData(GL_ARRAY_BUFFER, num_vertices_ * 2 * sizeof(float), texcoords_mat_.data(), GL_STATIC_DRAW);

	glEnableVertexAttribArray(texcoord_attribute_);
	glVertexAttribPointer(texcoord_attribute_, 2, GL_FLOAT, GL_FALSE, 2 * sizeof(float), (const GLvoid *)0);
}

GLuint QGLOcculsionTestWidget::load_shader(const char *vert_shader_file, const char *frag_shader_file)
{
	assert(vert_shader_file);
	assert(frag_shader_file);

	// "frag_color" is the name of output fragment color attributes.
	GLuint shader_program = InitShader(vert_shader_file, frag_shader_file, "frag_color");

	position_attribute_ = glGetAttribLocation(shader_program, "position");
	if (position_attribute_ < 0) {
		std::cerr << "Shader did not contain the 'position' attribute." << std::endl;
	}

	texcoord_attribute_ = glGetAttribLocation(shader_program, "texcoord");
	if (texcoord_attribute_ < 0) {
		std::cerr << "Shader did not contain the 'texcoord' attribute." << std::endl;
	}

	glsl_uniform_modelview_matrix_ = glGetUniformLocation(shader_program, "modelview_matrix");
	if (glsl_uniform_modelview_matrix_ < 0) {
		std::cerr << "Shader did not contain the 'modelview_matrix' uniform variable." << std::endl;
	}

	glsl_uniform_sample_points_ = glGetUniformLocation(shader_program, "sample_points");
	if (glsl_uniform_sample_points_ < 0) {
		std::cerr << "Shader did not contain the 'sample_points' uniform variable." << std::endl;
	}

	glsl_uniform_num_sample_points_ = glGetUniformLocation(shader_program, "num_sample_points");
	if (glsl_uniform_num_sample_points_ < 0) {
		std::cerr << "Shader did not contain the 'num_sample_points' uniform variable." << std::endl;
	}

	glsl_uniform_radius_ = glGetUniformLocation(shader_program, "radius");
	if (glsl_uniform_radius_ < 0) {
		std::cerr << "Shader did not contain the 'radius' uniform variable." << std::endl;
	}

	glsl_uniform_bbox_center_ = glGetUniformLocation(shader_program, "bbox_center");
	if (glsl_uniform_bbox_center_ < 0) {
		std::cerr << "Shader did not contain the 'bbox_center' uniform variable." << std::endl;
	}

	glsl_uniform_bbox_size_ = glGetUniformLocation(shader_program, "bbox_size");
	if (glsl_uniform_bbox_size_ < 0) {
		std::cerr << "Shader did not contain the 'bbox_size' uniform variable." << std::endl;
	}

	return shader_program;
}

void QGLOcculsionTestWidget::get_occlusion_test_result(const float _radius,
	std::array<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>, 4> &_occlusion_test_result)
{
	radius_ = _radius;

	for (unsigned int i = 0; i < 4; ++i)
		_occlusion_test_result[i].resize(window_width_, window_height_);

	// Reset the starting index to zero.
	sample_points_start_index_ = 0;

	GLenum buffer(GL_BACK);

	// Start to render.
	show();

	// Get pixel positions.
	qApp->processEvents();
	makeCurrent();
	updateGL();
	glFinish();

	glReadBuffer(buffer);
	glPixelStorei(GL_PACK_ALIGNMENT, 1);
	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
	paintGL();
	// RGB is mapped to XYZ.
	glReadPixels(0, 0, window_width_, window_height_, GL_RED, GL_FLOAT, _occlusion_test_result[0].data());
	glReadPixels(0, 0, window_width_, window_height_, GL_GREEN, GL_FLOAT, _occlusion_test_result[1].data());
	glReadPixels(0, 0, window_width_, window_height_, GL_BLUE, GL_FLOAT, _occlusion_test_result[2].data());
	assert(glGetError() == 0);
	
	// RGB has value [0, 1], which should be mapped to
	// [-0.55 * bbox_size_ + bbox_center_, +0.55 * bbox_size_ + bbox_center_].
	for (unsigned int i = 0; i < 3; ++i)
		_occlusion_test_result[i] =
		(_occlusion_test_result[i].array() - 0.5) * bbox_size_[i] * 1.1 + bbox_center_[i];


	// Get occlusion values.
	_occlusion_test_result[3].setConstant(1.0f);

	// NOTE:
	// At one time, only 'max_num_sample_points_' number of points can be rendered.
	// Thus, every 'max_num_sample_points_' points are rendered separately, and the results are merged.
	for (sample_points_start_index_ = 0; sample_points_start_index_ < num_sample_points_;
		sample_points_start_index_ += max_num_sample_points_)
	{
		Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> fbuffer(window_width_, window_height_);

		qApp->processEvents();
		makeCurrent();
		updateGL();
		glFinish();

		glReadBuffer(buffer);
		glPixelStorei(GL_PACK_ALIGNMENT, 1);
		glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
		paintGL();
		glReadPixels(0, 0, window_width_, window_height_, GL_ALPHA, GL_FLOAT, fbuffer.data());
		assert(glGetError() == 0);

		// Alpha = 1 means that the voxel is perfectly visible.
		// A voxel becomes invisible when it is occluded by any point.
		_occlusion_test_result[3] = _occlusion_test_result[3].array().min(fbuffer.array());
	}

	// Reset the starting index to zero.
	sample_points_start_index_ = 0;

	// Make this widget invisible.
	hide();
}

void QGLOcculsionTestWidget::initializeGL()
{
	// OpenGL Setup.
	glewExperimental = GL_TRUE;
	GLint GlewInitResult = glewInit();
	if (GlewInitResult != GLEW_OK)
		printf("ERROR: %s\n", glewGetErrorString(GlewInitResult));
	else if (GLEW_VERSION_3_3)
		std::cout << "Driver supports OpenGL 3.3\nDetails:" << std::endl;

	std::cout << "Using GLEW " << glewGetString(GLEW_VERSION) << std::endl;
	std::cout << "Vendor: " << glGetString(GL_VENDOR) << std::endl;
	std::cout << "Renderer: " << glGetString(GL_RENDERER) << std::endl;
	std::cout << "Version: " << glGetString(GL_VERSION) << std::endl;
	std::cout << "GLSL: " << glGetString(GL_SHADING_LANGUAGE_VERSION) << std::endl;

	GLint max_vertex_uniforms_vectors;
	glGetIntegerv(GL_MAX_VERTEX_UNIFORM_VECTORS, &max_vertex_uniforms_vectors);
	std::cout << "GL_MAX_VERTEX_UNIFORM_VECTORS: " << max_vertex_uniforms_vectors << std::endl;

	max_num_sample_points_ = static_cast<GLint>(max_vertex_uniforms_vectors / sizeof(GLfloat));


	generate_vertices();

	shader_program_ = load_shader(color_vert_shader_filename, color_frag_shader_filename);

	load_buffer_data();
}

void QGLOcculsionTestWidget::resizeGL(int _w, int _h)
{
	glViewport(0, 0, _w, _h);
	updateGL();
}

void QGLOcculsionTestWidget::paintGL()
{
	// NOTE:
	// At one time, only 'max_num_sample_points_' number of points can be rendered.
	// Thus, at most 'max_num_sample_points_' points are rendered from the given starting index.
	int num_partial_sample_points = num_sample_points_ - sample_points_start_index_;
	assert(num_partial_sample_points > 0);
	num_partial_sample_points = std::min(num_partial_sample_points, max_num_sample_points_);

	glClearColor(0.0, 0.0, 0.0, 1.0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glEnable(GL_ALPHA_TEST);

	// NOTE:
	// DO NOT USE OPENGL BLENDING FUNCTION HERE.

	glUseProgram(shader_program_);

	glUniformMatrix4fv(glsl_uniform_modelview_matrix_, 1, GL_FALSE, modelview_matrix_.data());
	glUniform3fv(glsl_uniform_sample_points_, num_partial_sample_points,
		static_cast<GLfloat *>(&sample_points_[sample_points_start_index_][0]));
	glUniform1i(glsl_uniform_num_sample_points_, num_partial_sample_points);
	glUniform1f(glsl_uniform_radius_, radius_);
	glUniform3fv(glsl_uniform_bbox_center_, 1, bbox_center_.data());
	glUniform3fv(glsl_uniform_bbox_size_, 1, bbox_size_.data());

	glDrawElements(GL_TRIANGLES, num_indices_, GL_UNSIGNED_INT, indices_vec_.data());
	assert(glGetError() == 0);
	glUseProgram(0);
}
