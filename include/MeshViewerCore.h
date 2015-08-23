#ifndef OPENMESHAPPS_VIEWERCORE_H
#define OPENMESHAPPS_VIEWERCORE_H

//== INCLUDES =================================================================

#include "MeshViewerCoreT.h"
#include "MyMesh.h"

#include "MeshCuboidStructure.h"
#include "MeshCuboidPredictor.h"
#include "MeshCuboidSymmetryGroup.h"
#include "MeshCuboidTrainer.h"


//== DEFINES ==================================================================

#define RECENT_LOAD				"recent_load.ini"
#define AUTO_LOAD_LB_OPERATOR	".laplacian.eigendata.dat"
#define CUSTOM_VIEW				"Custom view"
#define	COLORED_RENDERING		"Colored rendering"
#define	FACE_INDEX_RENDERING	"Face index rendering"
#define	POINT_SAMPLES			"Point rendering"
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

	bool open_modelview_matrix_file(const char* _filename);
	bool save_modelview_matrix_file(const char* _filename);
	bool save_projection_matrix_file(const char* _filename);

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
	void compute_ground_truth_cuboids();
	void train();
	void batch_predict();
	void predict();
	void run_part_assembly();
	void run_symmetry_detection();
	void run_symmetry_detection_msh2pln();
	void run_baseline_stats();
	void run_render_output();
	void run_render_evaluation();
	void run_extract_symmetry_info();
	void batch_render_point_clusters();
	void batch_render_cuboids();
	void batch_render_points();
	void print_mesh_info();
	void run_fusion_test();
	//void run_occlusion_test();


	// TEST.
	void test_initialize();
	void test_translate(const MyMesh::Normal _translation);
	void test_scale(const Real _scale_x, const Real _scale_y);
	void test_rotate(const Real _angle);
	void test_optimize();

	void test_figure_1();
	//


private:
	typedef enum {
		LoadMesh,
		LoadSamplePoints,
		LoadDenseSamplePoints,
		LoadGroundTruthCuboids,
		LoadTestData,
		LoadDenseTestData
	} LoadObjectInfoOption;
	
	bool load_object_info(
		MyMesh &_mesh,
		MeshCuboidStructure &_cuboid_structure,
		const char* _mesh_filepath,
		const LoadObjectInfoOption _option,
		const char* _cuboid_filepath = NULL,
		bool _verbose = true);

	bool load_result_info(
		MyMesh &_mesh,
		MeshCuboidStructure &_cuboid_structure,
		const char* _mesh_filepath,
		const char* _sample_filepath,
		const char* _sample_label_filepath,
		const char* _cuboid_filepath);

	void evaluate_per_point_labeling(
		const char *_mesh_filepath, 
		const char *_output_file_prefix);

	void reconstruct(
		const char *_mesh_filepath,
		const GLdouble *_snapshot_modelview_matrix,
		const GLdouble *_occlusion_modelview_matrix,
		const char *_output_file_prefix);

	void reconstruct_scan(
		const char *_mesh_filepath,
		const GLdouble *_snapshot_modelview_matrix,
		const GLdouble *_occlusion_modelview_matrix,
		const char *_output_file_prefix);

	void reconstruct_database_prior(
		const char *_mesh_filepath,
		const std::vector<LabelIndex> *_reconstructed_label_indices = NULL);

	void run_part_assembly_align_database(const std::string _mesh_filepath,
		Real &_xy_size, Real &_z_size, Real &_angle);

	void run_part_assembly_render_alignment(const std::string _mesh_filepath,
		const Real _xy_size, const Real _z_size, const Real _angle, const std::string _output_filename);

	void run_part_assembly_match_parts(const std::string _mesh_filepath,
		const Real _xy_size, const Real _z_size, const Real _angle,
		const MeshCuboidTrainer &_trainer, std::vector<std::string> &_label_matched_objects);

	void run_part_assembly_reconstruction(const std::string _mesh_filepath,
		const Real _xy_size, const Real _z_size, const Real _angle,
		const std::vector<std::string> &_label_matched_objects);

	void render_part_assembly_cuboids();
		

	void set_view_direction();

	void set_random_view_direction(bool _set_modelview_matrix = false);

	void compute_view_plane_mask_range(const Real _modelview_matrix[16]);


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
	bool draw_point_correspondences_;

	//QGLOcculsionTestWidget *occlusion_test_widget_;
	//std::vector< std::pair<MyMesh::Point, Real> > occlusion_test_points_;
	//bool draw_occlusion_test_points_;

	MyMesh::Point view_point_;
	MyMesh::Normal view_direction_;

	// TEST.
	std::vector< std::vector<MeshCuboidJointNormalRelations *> > test_joint_normal_relations_;
	MeshCuboidJointNormalRelationPredictor *test_joint_normal_predictor_;
	double test_occlusion_modelview_matrix_[16];
	std::vector<MeshCuboid *> all_cuboids_;
	//
};


#endif
