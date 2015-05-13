/*===========================================================================*\
 *                                                                           *
 *                               OpenMesh                                    *
 *      Copyright (C) 2001-2014 by Computer Graphics Group, RWTH Aachen      *
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
 *   $Revision: 990 $                                                         *
 *   $Date: 2014-02-05 10:01:07 +0100 (Mi, 05 Feb 2014) $                   *
 *                                                                           *
\*===========================================================================*/

#ifndef OPENMESHAPPS_VIEWERWIDGET_H
#define OPENMESHAPPS_VIEWERWIDGET_H

//== INCLUDES =================================================================

#include <QApplication>
#include <QWidget>
#include <QString>
#include <QMessageBox>
#include <QFileDialog>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Tools/Utils/getopt.h>
#include <OpenMesh/Tools/Utils/Timer.hh>

#include "QGLViewerWidget.h"
#include "MeshViewerCore.h"


//== CLASS DEFINITION =========================================================

class MeshViewerWidget : public QGLViewerWidget
{
    Q_OBJECT

public:
    /// default constructor
	MeshViewerWidget(QWidget* parent = 0);

	MeshViewerCore &get_mesh_viewer_core() { return mesh_viewer_core_; }

    void open_mesh_gui(QString fname)
    {
        OpenMesh::Utils::Timer t;
        t.start();
        if ( fname.isEmpty() || !mesh_viewer_core_.open_mesh(fname.toLocal8Bit()) )
        {
            QString msg = "Cannot read mesh from file:\n '";
            msg += fname;
            msg += "'";
            QMessageBox::critical( NULL, windowTitle(), msg );
        }
        t.stop();
        std::cout << "Loaded mesh in ~" << t.as_string() << std::endl;
    }

    //void open_texture_gui(QString fname)
    //{
    //    if ( fname.isEmpty() || !mesh_viewer_core_.open_texture( fname.toLocal8Bit() ) )
    //    {
    //        QString msg = "Cannot load texture image from file:\n '";
    //        msg += fname;
    //        msg += "'\n\nPossible reasons:\n";
    //        msg += "- Mesh file didn't provide texture coordinates\n";
    //        msg += "- Texture file does not exist\n";
    //        msg += "- Texture file is not accessible.\n";
    //        QMessageBox::warning( NULL, windowTitle(), msg );
    //    }
    //}

public slots:
    void query_open_mesh_file() {
        QString fileName = QFileDialog::getOpenFileName(this,
            tr("Open mesh file"),
            tr(""),
            tr("OFF Files (*.off);;"
			"OBJ Files (*.obj);;"
            "STL Files (*.stl);;"
            "All Files (*)"));
        if (!fileName.isEmpty())
            open_mesh_gui(fileName);
    }

    //void query_open_texture_file() {
    //    QString fileName = QFileDialog::getOpenFileName(this,
    //        tr("Open texture file"),
    //        tr(""),
    //        tr("PNG Files (*.png);;"
    //        "BMP Files (*.bmp);;"
    //        "GIF Files (*.gif);;"
    //        "JPEG Files (*.jpg);;"
    //        "TIFF Files (*.tif);;"
    //        "All Files (*)"));
    //    if (!fileName.isEmpty())
    //        open_texture_gui(fileName);
    //}

	void query_open_color_map_file() {
		QString fileName = QFileDialog::getOpenFileName(this,
			tr("Open a color map file"),
			tr(""),
			tr("Color maps (*.csv *.txt);;"
			"All Files (*)"));
		if (!fileName.isEmpty())
			mesh_viewer_core_.open_color_map_file(fileName.toLocal8Bit());
	}

	void query_save_color_map_file() {
		QString fileName = QFileDialog::getSaveFileName(this,
			tr("Open a color map file"),
			tr(""),
			tr("Color maps (*.csv *.txt);;"
			"All Files (*)"));
		if (!fileName.isEmpty())
			mesh_viewer_core_.save_color_map_file(fileName.toLocal8Bit());
	}

	void query_open_modelview_matrix_file()
	{
		QString fileName = QFileDialog::getOpenFileName(this,
			tr("Open a model view matrix file"),
			tr(""),
			tr("Text files (*.txt);;"
			"All Files (*)"));
		if (!fileName.isEmpty())
			mesh_viewer_core_.open_modelview_matrix_file(fileName.toLocal8Bit());
	}

	void query_save_modelview_matrix_file()
	{
		QString fileName = QFileDialog::getSaveFileName(this,
			tr("Save a model view matrix file"),
			tr(""),
			tr("Text files (*.txt);;"
			"All Files (*)"));
		if (!fileName.isEmpty())
			mesh_viewer_core_.save_modelview_matrix_file(fileName.toLocal8Bit());
	}

	void query_save_projection_matrix_file()
	{
		QString fileName = QFileDialog::getSaveFileName(this,
			tr("Save a projection matrix file"),
			tr(""),
			tr("Text files (*.txt);;"
			"All Files (*)"));
		if (!fileName.isEmpty())
			mesh_viewer_core_.save_projection_matrix_file(fileName.toLocal8Bit());
	}

	void query_open_sample_point_file()
	{
		QString fileName = QFileDialog::getOpenFileName(this,
			tr("Open a sample point file"),
			tr(""),
			tr("Sample Points (*.pts);;"
			"All Files (*)"));
		if (!fileName.isEmpty())
			mesh_viewer_core_.open_sample_point_file(fileName.toLocal8Bit());
	}

	void query_open_sample_point_label_file()
	{
		QString fileName = QFileDialog::getOpenFileName(this,
			tr("Open a sample point label file"),
			tr(""),
			tr("Attributes (*.arff);;"
			"All Files (*)"));
		if (!fileName.isEmpty())
			mesh_viewer_core_.open_sample_point_label_file(fileName.toLocal8Bit());
	}

	void query_open_mesh_face_label_file_and_create_cuboids()
	{
		QString fileName = QFileDialog::getOpenFileName(this,
			tr("Open a mesh face label file"),
			tr(""),
			tr("Labels (*.seg);;"
			"All Files (*)"));
		if (!fileName.isEmpty())
		{
			mesh_viewer_core_.open_mesh_face_label_file(fileName.toLocal8Bit());
			mesh_viewer_core_.create_mesh_cuboids();
		}
	}

	void query_open_face_label_file_and_apply_to_cuboids()
	{
		QString fileName = QFileDialog::getOpenFileName(this,
			tr("Open a mesh face label file"),
			tr(""),
			tr("Labels (*.seg);;"
			"All Files (*)"));
		if (!fileName.isEmpty())
		{
			mesh_viewer_core_.open_mesh_face_label_file(fileName.toLocal8Bit());
			mesh_viewer_core_.apply_mesh_labels_to_cuboids();
		}
	}

	void query_open_cuboid_file()
	{
		QString fileName = QFileDialog::getOpenFileName(this,
			tr("Open a cuboid file"),
			tr(""),
			tr("ARFF files (*.arff);;"
			"All Files (*)"));
		if (!fileName.isEmpty())
			mesh_viewer_core_.open_cuboid_file(fileName.toLocal8Bit());
	}

	void run_quit()	{ qApp->quit(); };


	void run_compute_ground_truth_cuboids() { mesh_viewer_core_.compute_ground_truth_cuboids(); }
	void run_train() { mesh_viewer_core_.train(); }
	void run_batch_predict() { mesh_viewer_core_.batch_predict(); }
	void run_predict() { mesh_viewer_core_.predict(); }
	void run_batch_render_cuboids() { mesh_viewer_core_.batch_render_cuboids(); }
	void run_batch_render_points() { mesh_viewer_core_.batch_render_points(); }
	void run_remove_occluded_points() { mesh_viewer_core_.remove_occluded_points(); }
	void run_print_arguments() { mesh_viewer_core_.parse_arguments(); }
	void run_print_mesh_info() { mesh_viewer_core_.print_mesh_info(); }
	void run_test() { mesh_viewer_core_.run_fusion_test(); }

	void run_test_initialize() { mesh_viewer_core_.test_initialize(); }
	void run_test_optimize() { mesh_viewer_core_.test_optimize(); }


protected: // inherited

	// initialize OpenGL states (triggered by Qt)
	virtual void initializeGL();

	// draw the scene (triggered by Qt)
	virtual void paintGL();

	// handle resize events (triggered by Qt)
	virtual void resizeGL(int w, int h);


	// Qt mouse events
	virtual void mousePressEvent(QMouseEvent* _event);
	virtual void mouseReleaseEvent(QMouseEvent* _event);
	virtual void mouseMoveEvent(QMouseEvent* _event);
	virtual void wheelEvent(QWheelEvent* _event);
	virtual void keyPressEvent(QKeyEvent* _event);


private:
	MeshViewerCore mesh_viewer_core_;
};


#endif
