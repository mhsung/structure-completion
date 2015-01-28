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


#ifndef OPENMESHAPPS_QGLVIEWERWIDGET_HH
#define OPENMESHAPPS_QGLVIEWERWIDGET_HH


//== INCLUDES =================================================================


#include <OpenMesh/Core/Geometry/VectorT.hh>
#include <QtOpenGL/qgl.h>
#include <string>
#include <vector>
#include <map>

#include "GLViewerCore.h"


//== FORWARD DECLARATIONS =====================================================

class QMenu;
class QActionGroup;
class QAction;

//== CLASS DEFINITION =========================================================


class QGLViewerWidget : public QGLWidget, public GLViewerBase
{
	Q_OBJECT

public:
	typedef QGLWidget Super;

	// Default constructor.
	QGLViewerWidget(QWidget* _parent = 0);

	// 
	QGLViewerWidget(QGLFormat& _fmt, QWidget* _parent = 0);

	// Destructor.
	virtual ~QGLViewerWidget();

private:

	void init(void);

public:

	// GLViewerWidget functions
	virtual float GLViewerBaseWidth() { return width(); }

	virtual float GLViewerBaseHeight() { return height(); }

	virtual void GLViewerBasemakeCurrent() { makeCurrent(); }

	virtual void GLViewerBaseSwapBuffers() { swapBuffers(); }

	virtual void GLViewerBaseUpdateGL() { updateGL(); };

	virtual void GLViewerBaseSnapshot(const std::string &_filename) { slotSnapshot(_filename); }

	virtual void GLViewerBaseSetWindowTitle(const std::string &_title) { setWindowTitle(_title.c_str()); }

	virtual std::string GLViewerBaseGetDrawMode()
	{
		std::string draw_mode("");
		if (draw_mode_)
		{
			assert(draw_mode_ <= n_draw_modes_);
			draw_mode = draw_mode_names_[draw_mode_ - 1];
		}
		return draw_mode;
	}

	virtual void GLViewerBaseSetDrawMode(const std::string &_draw_mode)
	{
		slotDrawMode(findAction(_draw_mode.c_str()));
	}


	/// add draw mode to popup menu, and return the QAction created
	QAction *add_draw_mode(const std::string& _s);

	/// delete draw mode from popup menu
	void del_draw_mode(const std::string& _s);

	const std::string& current_draw_mode() const
	{
		return draw_mode_ ? draw_mode_names_[draw_mode_ - 1] : nomode_;
	}

	float fovy() const { return 45.0f; }

	QAction* findAction(const char *name);
	void addAction(QAction* action, const char* name);
	void removeAction(const char* name);
	void removeAction(QAction* action);


protected:

	//double performance(void);


	private slots:

	// popup menu clicked
	void slotDrawMode(QAction *_mode);
	void slotSnapshot(const std::string _filename);


protected: // inherited

	// initialize OpenGL states (triggered by Qt)
	virtual void initializeGL();

	// draw the scene (triggered by Qt)
	virtual void paintGL();

	// handle resize events (triggered by Qt)
	virtual void resizeGL(int w, int h);


	// Qt mouse events
	virtual void mousePressEvent(QMouseEvent*);
	virtual void mouseReleaseEvent(QMouseEvent*);
	virtual void mouseMoveEvent(QMouseEvent*);
	virtual void wheelEvent(QWheelEvent*);
	virtual void keyPressEvent(QKeyEvent*);


private:

	// popup menu for draw mode selection
	QMenu*               popup_menu_;
	QActionGroup*        draw_modes_group_;
	typedef std::map<QString, QAction*> ActionMap;
	ActionMap            names_to_actions;
	unsigned int              draw_mode_;
	unsigned int              n_draw_modes_;
	std::vector<std::string>  draw_mode_names_;
	static std::string        nomode_;

};


//=============================================================================
#endif // OPENMESHAPPS_QGLVIEWERWIDGET_HH
//=============================================================================

