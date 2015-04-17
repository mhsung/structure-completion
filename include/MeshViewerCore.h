#ifndef OPENMESHAPPS_VIEWERCORE_H
#define OPENMESHAPPS_VIEWERCORE_H

//== INCLUDES =================================================================

#include "MeshViewerCoreT.h"
#include "MyMesh.h"

#include "MeshCuboidStructure.h"
#include "MeshCuboidPredictor.h"
#include "MeshCuboidSymmetryGroup.h"


//== DEFINES ==================================================================

#define RECENT_LOAD				"recent_load.ini"
#define AUTO_LOAD_LB_OPERATOR	".laplacian.eigendata.dat"
#define CUSTOM_VIEW				"Custom view"
#define	COLORED_RENDERING		"Colored rendering"
#define	FACE_INDEX_RENDERING	"Face index rendering"
#define	COLORED_POINT_SAMPLES	"Colored point samples"
//#define COLOR_RANDOM_SEED		20091129
#define COLOR_RANDOM_SEED		20091130

#define FACE_INFO_FILENAME		"faces.csv"
#define VERTEX_INFO_FILENAME	"vertices.csv"
#define NEIGHBOR_INFO_FILENAME	"neighbors.csv"


//== CLASS DEFINITION =========================================================

class MeshViewerCore : public MeshViewerCoreT<MyMesh>
{
public:
    /// default constructor
	MeshViewerCore(GLViewerBase &_widget);
	virtual ~MeshViewerCore();

	void open_modelview_matrix_file(const char* _filename);
	void save_modelview_matrix_file(const char* _filename);
	void save_projection_matrix_file(const char* _filename);

	void open_color_map_file(const char* _filename);
	void save_color_map_file(const char* _filename);

	void open_sample_point_file(const char* _filename);
	void open_sample_point_label_file(const char* _filename);
	void open_mesh_face_label_file(const char* _filename);
	void open_cuboid_file(const char* _filename);
	void create_mesh_cuboids();
	void apply_mesh_labels_to_cuboids();

	void remove_occluded_points();

	// Experiment functions.
	void parse_arguments();
	void train();
	void train_file_files();
	void batch_predict();
	void predict();
	void batch_render_point_clusters();
	void batch_render_cuboids();
	void batch_render_points();
	void print_mesh_info();
	//void do_occlusion_test();
	void run_test();
	//

	// TEST.
	void test_initialize();
	void test_translate(const MyMesh::Normal _translation);
	void test_scale(const Real _scale_x, const Real _scale_y);
	void test_rotate(const Real _angle);
	void test_optimize();
	//


private:
	bool load_training_data(
		const char* _mesh_filepath,
		bool _load_dense_samples = false);

	bool load_training_data(
		const char* _mesh_filepath,
		MyMesh &_mesh,
		MeshCuboidStructure &_cuboid_structure,
		bool _load_dense_samples = false);

	void reconstruct(
		const char *_mesh_filepath,
		const char *_dense_sample_filepath,
		const GLdouble *_occlusion_modelview_matrix,
		const char *_output_file_prefix);

	void reconstruct_using_database(
		const std::vector<LabelIndex> *_reconstructed_label_indices = NULL);

	void set_view_direction();

	void set_random_view_direction(bool _set_modelview_matrix = false);


public:
	virtual void mousePressEvent(const OpenMesh::Vec2f& _new_point_2d,
		bool _is_left_button, bool _is_mid_button, bool _is_right_button,
		bool _is_ctrl_pressed, bool _is_alt_pressed, bool _is_shift_pressed);

	virtual void keyPressEvent(const char _new_key,
		bool _is_ctrl_pressed, bool _is_alt_pressed, bool _is_shift_pressed);


protected:
	//Added: Minhyuk Sung. 2009-10-22
	void render_selection_mode();
	VertexIndex select_object(int x, int y);

	//Added: Minhyuk Sung. 2009-08-11 (Last modified: 2015-04-16)
	void red_color()	{ glColor4ub(192, 0, 0, 255); }
	void green_color()	{ glColor4ub(0, 176, 80, 255); }
	void blue_color()	{ glColor4ub(0, 112, 192, 255); }
	void yellow_color()	{ glColor4ub(255, 192, 0, 255); }
	void violet_color()	{ glColor4ub(112, 48, 160, 255); }
	void black_color()	{ glColor4ub(0, 32, 96, 255); }
	void gray_color()	{ glColor4ub(96, 96, 96, 255); }

	/// draw the mesh
	virtual void draw_scene(const std::string& _draw_mode);
	virtual void draw_openmesh(const std::string& _drawmode);


private:
	double point_size_;
	enum
	{
		PICK_SEED,
		PICK_QUERY,
		PICK_FEATURE,
	}
	selection_mode_;

	MeshCuboidStructure cuboid_structure_;

	bool draw_cuboid_axes_;

	//QGLOcculsionTestWidget *occlusion_test_widget_;
	//std::vector< std::pair<MyMesh::Point, Real> > occlusion_test_points_;
	//bool draw_occlusion_test_points_;

	MyMesh::Point view_point_;
	MyMesh::Normal view_direction_;

	// TEST.
	std::vector< std::vector<MeshCuboidJointNormalRelations *> > joint_normal_relations_;
	MeshCuboidJointNormalRelationPredictor *test_joint_normal_predictor_;
	double test_occlusion_modelview_matrix_[16];
	std::vector<MeshCuboid *> all_cuboids_;
	//
};


#endif
