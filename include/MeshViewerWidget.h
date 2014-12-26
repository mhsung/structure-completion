/*===========================================================================*\
 *                                                                           *
 *                               OpenMesh                                    *
 *      Copyright (C) 2001-2009 by Computer Graphics Group, RWTH Aachen      *
 *                           www.openmesh.org                                *
 *                                                                           *
 *---------------------------------------------------------------------------* 
 *  This file is part of OpenMesh.                                           *
 *                                                                           *
 *  OpenMesh is free software: you can redistribute it and/or modify         * 
 *  it under the terms of the GNU Lesser General Public License as           *
 *  published by the Free Software Foundation, either version 3 of           *
 *  the License, or (at your option) any later version with the              *
 *  following exceptions:                                                    *
 *                                                                           *
 *  If other files instantiate templates or use macros                       *
 *  or inline functions from this file, or you compile this file and         *
 *  link it with other files to produce an executable, this file does        *
 *  not by itself cause the resulting executable to be covered by the        *
 *  GNU Lesser General Public License. This exception does not however       *
 *  invalidate any other reasons why the executable file might be            *
 *  covered by the GNU Lesser General Public License.                        *
 *                                                                           *
 *  OpenMesh is distributed in the hope that it will be useful,              *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of           *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            *
 *  GNU Lesser General Public License for more details.                      *
 *                                                                           *
 *  You should have received a copy of the GNU LesserGeneral Public          *
 *  License along with OpenMesh.  If not,                                    *
 *  see <http://www.gnu.org/licenses/>.                                      *
 *                                                                           *
\*===========================================================================*/ 

/*===========================================================================*\
 *                                                                           *             
 *   $Revision: 137 $                                                         *
 *   $Date: 2009-06-04 10:46:29 +0200 (Do, 04. Jun 2009) $                   *
 *                                                                           *
\*===========================================================================*/

#ifndef OPENMESHAPPS_VIEWERWIDGET_HH
#define OPENMESHAPPS_VIEWERWIDGET_HH


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



//== INCLUDES =================================================================

#include <QApplication>
#include <QWidget>
#include <QString>
#include <QMessageBox>
#include <QFileDialog>
#include <OpenMesh/Tools/Utils/getopt.h>
#include <OpenMesh/Tools/Utils/Timer.hh>
#include "MeshViewerWidgetT.h"
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include "MyMesh.h"
#include "MeshCuboidStructure.h"
#include "MeshCuboidPredictor.h"
#include "QGLOcculsionTestWidget.h"


//== CLASS DEFINITION =========================================================

using namespace OpenMesh;  
using namespace OpenMesh::Attributes;

// Modified: Minhyuk Sung. 2009-08-09
/*
struct MyTraits : public OpenMesh::DefaultTraits
{
  HalfedgeAttributes(OpenMesh::Attributes::PrevHalfedge);
};

typedef OpenMesh::TriMesh_ArrayKernelT<MyTraits>  MyMesh;
*/


//== CLASS DEFINITION =========================================================

class MeshViewerWidget : public MeshViewerWidgetT<MyMesh>
{
    Q_OBJECT

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

	QGLOcculsionTestWidget *occlusion_test_widget_;
	std::vector< std::pair<MyMesh::Point, Real> > occlusion_test_points_;
	bool draw_occlusion_test_points_;

	MyMesh::Point view_point_;
	MyMesh::Normal view_direction_;

	// TEST
	std::vector< std::vector<MeshCuboidJointNormalRelations> > joint_normal_relations_;
	MeshCuboidJointNormalRelationPredictor *test_joint_normal_predictor_;
	double test_occlusion_modelview_matrix_[16];
	std::vector<MeshCuboid *> all_cuboids_;
	//


public:
    /// default constructor
    MeshViewerWidget(QWidget* parent = 0)
		: MeshViewerWidgetT<MyMesh>(parent)
		, point_size_(1.0)
		, selection_mode_(PICK_SEED)
		, cuboid_structure_(&mesh_)
		, draw_cuboid_axes_(true)
		, draw_occlusion_test_points_(false)
		, view_point_(0.0)
		, view_direction_(0.0)
    {
		// NOTE:
		// Enable alpha channel.
		QGLFormat fmt;
		fmt.setAlpha(true);
		occlusion_test_widget_ = new QGLOcculsionTestWidget(fmt);

		// TEST
		test_joint_normal_predictor_ = NULL;
		//
	}

	~MeshViewerWidget()
	{
		delete occlusion_test_widget_;

		// TEST
		delete test_joint_normal_predictor_;
		//
	}

	OpenMesh::IO::Options& options() { return _options; }
    const OpenMesh::IO::Options& options() const { return _options; }
    void setOptions(const OpenMesh::IO::Options& opts) { _options = opts; }

    void open_mesh_gui(QString fname)
    {
		OpenMesh::Utils::Timer t;
        t.start();
        if ( fname.isEmpty() || !open_mesh(fname.toLocal8Bit(), _options) )
        {
            QString msg = "Cannot read mesh from file:\n '";
            msg += fname;
            msg += "'";
            QMessageBox::critical( NULL, windowTitle(), msg);
        }

        t.stop();
        std::cout << "Loaded mesh in ~" << t.as_string() << std::endl;
    }
    void open_texture_gui(QString fname)
    {
        if ( fname.isEmpty() || !open_texture( fname.toLocal8Bit() ) )
        {
            QString msg = "Cannot load texture image from file:\n '";
            msg += fname;
            msg += "'\n\nPossible reasons:\n";
            msg += "- Mesh file didn't provide texture coordinates\n";
            msg += "- Texture file does not exist\n";
            msg += "- Texture file is not accessible.\n";
            QMessageBox::warning( NULL, windowTitle(), msg );
        }
    }


public slots:
	void query_open_mesh_file() {
		QString fileName = QFileDialog::getOpenFileName(this,
			tr("Open mesh file"),
			tr(""),
			tr("Meshes (*.obj *.off *.ply);;"
			"All Files (*)"));
		if (!fileName.isEmpty())
			open_mesh_gui(fileName);
	}

	void query_open_texture_file() {
		QString fileName = QFileDialog::getOpenFileName(this,
			tr("Open texture file"),
			tr(""),
			tr("Images (*.bmp *.jpg *.png *.tif);;"
			"All Files (*)"));
		if (!fileName.isEmpty())
			open_texture_gui(fileName);
	}

	void query_open_sample_point_file()
	{
		QString fileName = QFileDialog::getOpenFileName(this,
			tr("Open a sample point file"),
			tr(""),
			tr("Sample Points (*.pts);;"
			"All Files (*)"));
		if (!fileName.isEmpty())
			open_sample_point_file(fileName);
	}

	void query_open_sample_point_label_file()
	{
		QString fileName = QFileDialog::getOpenFileName(this,
			tr("Open a sample point label file"),
			tr(""),
			tr("Attributes (*.arff);;"
			"All Files (*)"));
		if (!fileName.isEmpty())
			open_sample_point_label_file(fileName);
	}

	void query_open_face_label_file()
	{
		QString fileName = QFileDialog::getOpenFileName(this,
			tr("Open a mesh face label file"),
			tr(""),
			tr("Labels (*.comp);;"
			"Labels (*.class);;"
			"All Files (*)"));
		if (!fileName.isEmpty())
			open_face_label_file(fileName);
	}

	void query_open_face_label_file_but_preserve_cuboids()
	{
		QString fileName = QFileDialog::getOpenFileName(this,
			tr("Open a mesh face label file"),
			tr(""),
			tr("Labels (*.comp);;"
			"Labels (*.class);;"
			"All Files (*)"));
		if (!fileName.isEmpty())
			open_face_label_file_but_preserve_cuboids(fileName);
	}

	void query_open_modelview_matrix_file()
	{
		QString fileName = QFileDialog::getOpenFileName(this,
			tr("Open a model view matrix file"),
			tr(""),
			tr("Text files (*.txt);;"
			"All Files (*)"));
		if (!fileName.isEmpty())
			open_modelview_matrix_file(fileName);
	}

	void query_save_modelview_matrix_file()
	{
		QString fileName = QFileDialog::getSaveFileName(this,
			tr("Save a model view matrix file"),
			tr(""),
			tr("Text files (*.txt);;"
			"All Files (*)"));
		if (!fileName.isEmpty())
			save_modelview_matrix_file(fileName);
	}

	void query_save_projection_matrix_file()
	{
		QString fileName = QFileDialog::getSaveFileName(this,
			tr("Save a projection matrix file"),
			tr(""),
			tr("Text files (*.txt);;"
			"All Files (*)"));
		if (!fileName.isEmpty())
			save_projection_matrix_file(fileName);
	}

	void open_sample_point_file(QString filename);
	void open_sample_point_label_file(QString filename);
	void open_face_label_file(QString filename);
	void open_face_label_file_but_preserve_cuboids(QString filename);

	void open_modelview_matrix_file(QString filename);
	void save_modelview_matrix_file(QString filename);
	void save_projection_matrix_file(QString filename);

	void query_open_color_map_file();
	void query_save_color_map_file();

	void remove_occluded_points();

	void run_quit()	{ qApp->quit(); };

	// Nov. 2012, Minhyuk Sung
	void run_print_mesh_info();
	void run_training();
	void run_prediction();
	void run_rendering_point_clusters();
	void run_occlusion_test();
	void set_view_direction();
	//void run_test_joint_normal_training();
	//void run_test_manual_relations();
	//void run_test_cca_relations();

	// TEST
	void run_test_initialize();
	void run_test_translate(const MyMesh::Normal _translation);
	void run_test_scale(const Real _scale_x, const Real _scale_y);
	void run_test_rotate(const Real _angle);
	void run_test_optimize();
	//


private:
    OpenMesh::IO::Options _options;


protected:
	//Added: Minhyuk Sung. 2009-10-22
	void render_selection_mode();
	VertexIndex select_object(int x, int y);
	virtual void mousePressEvent( QMouseEvent* _event );
	virtual void keyPressEvent( QKeyEvent* _event);

	//Added: Minhyuk Sung. 2009-08-11
	void red_color()	{ glColor4f(192.0f / 256.0f,   0.0f / 256.0f,   0.0f / 256.0f, 1.0f); }
	void green_color()	{ glColor4f(  0.0f / 256.0f, 176.0f / 256.0f,  80.0f / 256.0f, 1.0f); }
	void blue_color()	{ glColor4f(  0.0f / 256.0f, 112.0f / 256.0f, 192.0f / 256.0f, 1.0f); }
	void yellow_color()	{ glColor4f(255.0f / 256.0f, 192.0f / 256.0f,   0.0f / 256.0f, 1.0f); }
	void violet_color()	{ glColor4f(112.0f / 256.0f,  48.0f / 256.0f, 160.0f / 256.0f, 1.0f); }
	void black_color()	{ glColor4f(  0.0f / 256.0f,  32.0f / 256.0f,  96.0f / 256.0f, 1.0f); }

	/// draw the mesh
	virtual void draw_scene(const std::string& _draw_mode)
	{
		MeshViewerWidgetT::draw_scene(_draw_mode);
		if ( _draw_mode == CUSTOM_VIEW )
		{
			glDisable(GL_LIGHTING);
			glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
			draw_openmesh(_draw_mode);
			glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		}
		else if( _draw_mode == COLORED_RENDERING )
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

	virtual void draw_openmesh(const std::string& _drawmode);
};


#endif
