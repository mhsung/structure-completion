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


#ifndef OPENMESHAPPS_QGLVIEWERWIDGET_HH
#define OPENMESHAPPS_QGLVIEWERWIDGET_HH

// If defined, both A, B objects move together
#define MOVE_BOTH_TOGETHER


//== INCLUDES =================================================================


#include <OpenMesh/Core/Geometry/VectorT.hh>
#include <QtOpenGL/qgl.h>
#include <string>
#include <vector>
#include <map>


//== FORWARD DECLARATIONS =====================================================

class QMenu;
class QActionGroup;
class QAction;

//== CLASS DEFINITION =========================================================


class QGLViewerWidget : public QGLWidget
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

	typedef enum
	{
		VIEW_A,
		VIEW_B
	} View;

	View get_current_view()
	{
		return view_;
	};

	virtual void set_current_view(View v)
	{
		view_ = v;
		if (view_ == VIEW_A) modelview_matrix_ = modelview_matrix_A;
		if (view_ == VIEW_B) modelview_matrix_ = modelview_matrix_B;
	};

	public slots:

	virtual void set_render_view_A()
	{
		render_modelview_matrix_ = modelview_matrix_A;
	};

	virtual void set_render_view_B()
	{
		render_modelview_matrix_ = modelview_matrix_B;
	};

	void toggle_view()
	{
		if (get_current_view() == VIEW_A)		set_current_view(VIEW_B);
		else if (get_current_view() == VIEW_B)	set_current_view(VIEW_A);
	};

	// Added: Minhyuk Sung. 2013-11-11
	void toggle_hide_A() { hide_A_ ^= 1; };
	void toggle_hide_B() { hide_B_ ^= 1; };


public:

	/* Sets the center and size of the whole scene.
	The _center is used as fixpoint for rotations and for adjusting
	the camera/viewer (see view_all()). */
	void set_scene_pos(const OpenMesh::Vec3f& _center, float _radius);

	/* view the whole scene: the eye point is moved far enough from the
	center so that the whole scene is visible. */
	void view_all();

	/// add draw mode to popup menu, and return the QAction created
	QAction *add_draw_mode(const std::string& _s);

	/// delete draw mode from popup menu
	void del_draw_mode(const std::string& _s);

	const std::string& current_draw_mode() const
	{
		return draw_mode_ ? draw_mode_names_[draw_mode_ - 1] : nomode_;
	}

	float radius() const { return radius_; }
	const OpenMesh::Vec3f& center() const { return center_; }

	const GLdouble* modelview_matrix() const  { return modelview_matrix_; }
	const GLdouble* projection_matrix() const { return projection_matrix_; }

	float fovy() const { return 45.0f; }

	QAction* findAction(const char *name);
	void addAction(QAction* action, const char* name);
	void removeAction(const char* name);
	void removeAction(QAction* action);

protected:

	// draw the scene: will be called by the painGL() method.
	virtual void draw_scene(const std::string& _draw_mode);

	double performance(void);

	void setDefaultMaterial(void);
	void setDefaultLight(void);

	// Modified: Minhyuk Sung. 2013-11-07
	// private -> protected
	// void slotSnapshot( void ) -> void slotSnapshot( const char *_filename = NULL )
	protected slots:

	// popup menu clicked
	void slotDrawMode(QAction *_mode);
	void slotSnapshot(const char *_filename = NULL);


private: // inherited

	// initialize OpenGL states (triggered by Qt)
	void initializeGL();

	// draw the scene (triggered by Qt)
	void paintGL();

	// handle resize events (triggered by Qt)
	void resizeGL(int w, int h);

protected:

	// Qt mouse events
	virtual void mousePressEvent(QMouseEvent*);
	virtual void mouseReleaseEvent(QMouseEvent*);
	virtual void mouseMoveEvent(QMouseEvent*);
	virtual void wheelEvent(QWheelEvent*);
	virtual void keyPressEvent(QKeyEvent*);

	// Added: Minhyuk Sung. 2009-11-24
	void set_modelview_matrix(GLdouble *matrix);

	static void getTransformedCoord(const GLdouble *matrix,
		const GLdouble *_coord, GLdouble output[3]);

	// Modified: Minhyuk Sung. 2013-11-07
	// private -> protected
protected:

	// updates projection matrix
	void update_projection_matrix();

	// translate the scene and update modelview matrix
	void translate(const OpenMesh::Vec3f& _trans);

	// rotate the scene (around its center) and update modelview matrix
	void rotate(const OpenMesh::Vec3f& _axis, float _angle);

	OpenMesh::Vec3f  center_;
	float            radius_;

	GLdouble    projection_matrix_[16];
	GLdouble    *modelview_matrix_;
	GLdouble    *render_modelview_matrix_;

protected:
	View view_;
	GLdouble    modelview_matrix_A[16];
	GLdouble    modelview_matrix_B[16];


	// popup menu for draw mode selection
	QMenu*               popup_menu_;
	QActionGroup*        draw_modes_group_;
	typedef std::map<QString, QAction*> ActionMap;
	ActionMap            names_to_actions;
	unsigned int              draw_mode_;
	unsigned int              n_draw_modes_;
	std::vector<std::string>  draw_mode_names_;
	static std::string        nomode_;

	// Added: Minhyuk Sung. 2013-11-11
	bool hide_A_;
	bool hide_B_;


	// virtual trackball: map 2D screen point to unit sphere
	bool map_to_sphere(const QPoint& _point, OpenMesh::Vec3f& _result);

	QPoint           last_point_2D_;
	OpenMesh::Vec3f  last_point_3D_;
	bool             last_point_ok_;

};


//=============================================================================
#endif // OPENMESHAPPS_QGLVIEWERWIDGET_HH
//=============================================================================

