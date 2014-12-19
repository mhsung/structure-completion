#ifndef _QGL_OCCLUSION_TEST_WIDGET_H_
#define _QGL_OCCLUSION_TEST_WIDGET_H_

#define DEFAULT_RESOLUTION_UNIT		8
#define DAFAULT_SAMPLE_POINT_RADIUS	0.01

#include <array>

#include <QtOpenGL/qgl.h>
#include <Eigen/Core>
#include <Eigen/StdVector>

// NOTE:
// ALPHA CHANNEL SHOULD BE ENABLED BEFORE CONSTRUCTING AN INSTANCE.
// Example:
// QGLFormat fmt;
// fmt.setAlpha(true);
// occlusion_test_widget_ = new QGLOcculsionTestWidget(fmt);


class QGLOcculsionTestWidget : public QGLWidget
{
	Q_OBJECT

public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

	typedef QGLWidget Super;

	// Default constructor.
	QGLOcculsionTestWidget(QWidget* _parent = 0);

	// 
	QGLOcculsionTestWidget(QGLFormat& _fmt, QWidget* _parent = 0);

	// Destructor.
	virtual ~QGLOcculsionTestWidget();

	void set_sample_points(const std::vector<Eigen::Vector3f> &_sample_points);

	void set_modelview_matrix(const GLdouble _modelview_matrix[16]);

	// NOTE:
	// The value '1' means that the voxel is perfectly visible.
	void get_occlusion_test_result(const float _radius,
		std::array<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>, 4> &_occlusion_test_result);


private: // inherited
	// initialize OpenGL states (triggered by Qt)
	void initializeGL();

	// draw the scene (triggered by Qt)
	void paintGL();

	// handle resize events (triggered by Qt)
	void resizeGL(int w, int h);


private:
	void initialize(void);

	bool load_sample_points(const char *_filename, bool _verbose = true);

	void get_bounding_box();

	void generate_vertices();

	void load_buffer_data();

	GLuint load_shader(const char *vert_shader_file, const char *frag_shader_file);


private:
	int window_width_;
	int window_height_;

	int resolution_unit_;
	GLint max_num_sample_points_;

	int resolution_unit_square_;
	int resolution_unit_cubic_;
	
	GLuint shader_program_;

	GLuint vertex_array_object_;

	GLuint position_buffer_;
	GLint position_attribute_;

	GLuint texcoord_buffer_;
	GLint texcoord_attribute_;

	GLint glsl_uniform_modelview_matrix_;
	Eigen::Matrix<float, 4, 4, Eigen::ColMajor> modelview_matrix_;
	
	GLint glsl_uniform_sample_points_;
	std::vector<Eigen::Vector3f> sample_points_;
	int sample_points_start_index_;

	GLint glsl_uniform_num_sample_points_;
	int num_sample_points_;

	GLint glsl_uniform_radius_;
	GLfloat radius_;

	GLint glsl_uniform_bbox_center_;
	Eigen::Vector3f bbox_center_;

	GLint glsl_uniform_bbox_size_;
	Eigen::Vector3f bbox_size_;

	unsigned int num_vertices_;
	unsigned int num_indices_;

	Eigen::Matrix<float, Eigen::Dynamic, 3, Eigen::RowMajor> positions_mat_;
	Eigen::Matrix<float, Eigen::Dynamic, 2, Eigen::RowMajor> texcoords_mat_;
	Eigen::Matrix<unsigned int, Eigen::Dynamic, 1> indices_vec_;

	float angle_;
};

#endif	// _QGL_OCCLUSION_TEST_WIDGET_H_