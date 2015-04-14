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
#define CUSTOM_VIEW				"Custom_View"
#define	COLORED_RENDERING		"Colored_Rendering"
#define	FACE_INDEX_RENDERING	"Face_Index_Rendering"
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

	virtual bool open_mesh(const char* _filename);

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
	void train();
	void train_file_files();
	void batch_predict();
	void predict();
	void batch_render_point_clusters();
	void batch_render_cuboids();
	void batch_render_points();
	void set_view_direction();
	void set_random_view_direction(bool _set_modelview_matrix = false);
	void parse_arguments();
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
	void reconstruct(
		const char *mesh_filepath,
		const char *dense_sample_filepath,
		const GLdouble *occlusion_modelview_matrix,
		const char *output_file_prefix);
	void reconstruct_using_database(
		const std::vector<LabelIndex> *_reconstructed_label_indices = NULL);


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

	//Added: Minhyuk Sung. 2009-08-11
	void red_color()	{ glColor4f(192.0f / 256.0f, 0.0f / 256.0f, 0.0f / 256.0f, 1.0f); }
	void green_color()	{ glColor4f(0.0f / 256.0f, 176.0f / 256.0f, 80.0f / 256.0f, 1.0f); }
	void blue_color()	{ glColor4f(0.0f / 256.0f, 112.0f / 256.0f, 192.0f / 256.0f, 1.0f); }
	void yellow_color()	{ glColor4f(255.0f / 256.0f, 192.0f / 256.0f, 0.0f / 256.0f, 1.0f); }
	void violet_color()	{ glColor4f(112.0f / 256.0f, 48.0f / 256.0f, 160.0f / 256.0f, 1.0f); }
	void black_color()	{ glColor4f(0.0f / 256.0f, 32.0f / 256.0f, 96.0f / 256.0f, 1.0f); }
	void gray_color()	{ glColor4f(127.0f / 256.0f, 127.0f / 256.0f, 127.0f / 256.0f, 1.0f); }

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
