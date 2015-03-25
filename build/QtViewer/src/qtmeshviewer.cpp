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

#ifdef _MSC_VER
#  pragma warning(disable: 4267 4311)
#endif

#define GLOG_NO_ABBREVIATED_SEVERITIES

#ifdef ARCH_DARWIN
	#include <glut.h>
#else
	#include <GL/glut.h>
	//#include <GL/glew.h>
	//#include <GL/freeglut.h>
#endif

#include <iostream>
#include <fstream>

#include <QApplication>
#include <QMessageBox>
#include <QMainWindow>
#include <QMenuBar>
#include <QFileDialog>

#include <gflags/gflags.h>
//#include <glog/logging.h>

#include "MeshViewerWidget.h"


void create_menu(QMainWindow &w);
void usage_and_exit(int xcode);

int main(int argc, char **argv)
{
	gflags::ParseCommandLineFlags(&argc, &argv, true);
	//google::InitGoogleLogging(argv[0]);

	// OpenGL check
	QApplication::setColorSpec(QApplication::CustomColor);
	QApplication app(argc, argv);

	glutInit(&argc, argv);

	//glutInitContextVersion(3, 3);
	//glutInitContextFlags(GLUT_FORWARD_COMPATIBLE);
	//glutInitContextProfile(GLUT_CORE_PROFILE);

	if (!QGLFormat::hasOpenGL()) {
		QString msg = "System has no OpenGL support!";
		QMessageBox::critical(0, QString("OpenGL"), msg + QString(argv[1]));
		return -1;
	}

	int c;
	OpenMesh::IO::Options opt;

	while ((c = getopt(argc, argv, "hbs")) != -1)
	{
		switch (c)
		{
		case 'b': opt += OpenMesh::IO::Options::Binary; break;
		case 'h':
			usage_and_exit(0);
		case 's': opt += OpenMesh::IO::Options::Swap; break;
		default:
			usage_and_exit(1);
		}
	}
	// create widget
	QMainWindow mainWin;
	MeshViewerWidget w(&mainWin);
	w.get_mesh_viewer_core().setOptions(opt);
	mainWin.setCentralWidget(&w);

	create_menu(mainWin);

	// static mesh, hence use strips
	// w.enable_strips();

	// Added: Minhyuk Sung. 2009-08-09
	//w.add_draw_mode("Colored Vertices");
	//w.add_draw_mode("Solid Colored Faces");
	QAction *a;
	a = w.add_draw_mode(CUSTOM_VIEW);
	a->setShortcut(QKeySequence(Qt::Key_D));
	a = w.add_draw_mode(COLORED_RENDERING);
	a = w.add_draw_mode(FACE_INDEX_RENDERING);

	mainWin.resize(800, 600 + 21);
	mainWin.show();

	//
	w.get_mesh_viewer_core().parse_arguments();
	//

	// load scene if specified on the command line
	if (optind < argc)
	{
		w.open_mesh_gui(argv[optind]);
	}

	//if (++optind < argc)
	//{
	//	w.open_texture_gui(argv[optind]);
	//}

	return app.exec();
}

void create_menu(QMainWindow &w)
{
	using namespace Qt;
	QAction* runAct;

	/* ------------------------	*/
	/* File Menu				*/
	/* ------------------------	*/

	QMenu *fileMenu = w.menuBar()->addMenu(w.tr("&File"));

	QAction* openAct = new QAction(w.tr("&Open mesh..."), &w);
	openAct->setShortcut(w.tr("Ctrl+O"));
	openAct->setStatusTip(w.tr("Open a mesh file"));
	QObject::connect(openAct, SIGNAL(triggered()), w.centralWidget(), SLOT(query_open_mesh_file()));
	fileMenu->addAction(openAct);

	//QAction* texAct = new QAction(w.tr("Open &texture..."), &w);
	//texAct->setShortcut(w.tr("Ctrl+T"));
	//texAct->setStatusTip(w.tr("Open a texture file"));
	//QObject::connect(texAct, SIGNAL(triggered()), w.centralWidget(), SLOT(query_open_texture_file()));
	//fileMenu->addAction(texAct);

	fileMenu->addSeparator();

	runAct = new QAction(w.tr("Open sample points..."), &w);
	runAct->setStatusTip(w.tr("Open a feature vertices file"));
	QObject::connect(runAct, SIGNAL(triggered()), w.centralWidget(), SLOT(query_open_sample_point_file()));
	fileMenu->addAction(runAct);

	runAct = new QAction(w.tr("Open sample point labels..."), &w);
	runAct->setStatusTip(w.tr("Save a feature vertices file"));
	QObject::connect(runAct, SIGNAL(triggered()), w.centralWidget(), SLOT(query_open_sample_point_label_file()));
	fileMenu->addAction(runAct);

	runAct = new QAction(w.tr("Open mesh face labels (Create cuboids)..."), &w);
	QObject::connect(runAct, SIGNAL(triggered()), w.centralWidget(),
		SLOT(query_open_mesh_face_label_file_and_create_cuboids()));
	fileMenu->addAction(runAct);

	runAct = new QAction(w.tr("Open mesh face labels (Apply to cuboids)..."), &w);
	QObject::connect(runAct, SIGNAL(triggered()), w.centralWidget(),
		SLOT(query_open_face_label_file_and_apply_to_cuboids()));
	fileMenu->addAction(runAct);

	runAct = new QAction(w.tr("Open cuboids..."), &w);
	QObject::connect(runAct, SIGNAL(triggered()), w.centralWidget(), SLOT(query_open_cuboid_file()));
	fileMenu->addAction(runAct);

	fileMenu->addSeparator();

	runAct = new QAction(w.tr("Open model view matrix..."), &w);
	runAct->setStatusTip(w.tr("Open a model view matrix file"));
	QObject::connect(runAct, SIGNAL(triggered()), w.centralWidget(), SLOT(query_open_modelview_matrix_file()));
	fileMenu->addAction(runAct);

	runAct = new QAction(w.tr("Save model view matrix..."), &w);
	runAct->setStatusTip(w.tr("Save a model view matrix file"));
	QObject::connect(runAct, SIGNAL(triggered()), w.centralWidget(), SLOT(query_save_modelview_matrix_file()));
	fileMenu->addAction(runAct);

	runAct = new QAction(w.tr("Save projection matrix..."), &w);
	runAct->setStatusTip(w.tr("Save a projection matrix file"));
	QObject::connect(runAct, SIGNAL(triggered()), w.centralWidget(), SLOT(query_save_projection_matrix_file()));
	fileMenu->addAction(runAct);

	fileMenu->addSeparator();

	runAct = new QAction(w.tr("&Quit"), &w);
	QObject::connect(runAct, SIGNAL(triggered()), w.centralWidget(), SLOT(run_quit()));
	fileMenu->addAction(runAct);



	/* ------------------------	*/
	/* Experiment Menu			*/
	/* ------------------------	*/

	runAct = new QAction(w.tr("Train"), &w);
	QObject::connect(runAct, SIGNAL(triggered()), w.centralWidget(), SLOT(run_train()));
	w.menuBar()->addAction(runAct);

	runAct = new QAction(w.tr("Train (From files)"), &w);
	QObject::connect(runAct, SIGNAL(triggered()), w.centralWidget(), SLOT(run_train_file_files()));
	w.menuBar()->addAction(runAct);

	runAct = new QAction(w.tr("Batch Predict"), &w);
	QObject::connect(runAct, SIGNAL(triggered()), w.centralWidget(), SLOT(run_batch_predict()));
	w.menuBar()->addAction(runAct);

	runAct = new QAction(w.tr("Render Cuboids"), &w);
	QObject::connect(runAct, SIGNAL(triggered()), w.centralWidget(), SLOT(run_batch_render_cuboids()));
	w.menuBar()->addAction(runAct);

	runAct = new QAction(w.tr("Render Points"), &w);
	QObject::connect(runAct, SIGNAL(triggered()), w.centralWidget(), SLOT(run_batch_render_points()));
	w.menuBar()->addAction(runAct);

	runAct = new QAction(w.tr("Set View Direction"), &w);
	QObject::connect(runAct, SIGNAL(triggered()), w.centralWidget(), SLOT(run_set_view_direction()));
	w.menuBar()->addAction(runAct);


	// TEST.
	runAct = new QAction(w.tr("[Test] Initialize"), &w);
	QObject::connect(runAct, SIGNAL(triggered()), w.centralWidget(), SLOT(run_test_initialize()));
	w.menuBar()->addAction(runAct);

	runAct = new QAction(w.tr("[Test] Optimize"), &w);
	QObject::connect(runAct, SIGNAL(triggered()), w.centralWidget(), SLOT(run_test_optimize()));
	w.menuBar()->addAction(runAct);
	//
}

void usage_and_exit(int xcode)
{
	std::cout << "Usage: meshviewer [-s] [mesh] [texture]\n" << std::endl;
	std::cout << "Options:\n"
		<< "  -b\n"
		<< "    Assume input to be binary.\n\n"
		<< "  -s\n"
		<< "    Reverse byte order, when reading binary files.\n"
		<< std::endl;
	exit(xcode);
}
