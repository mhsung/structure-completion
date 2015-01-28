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

//== INCLUDES =================================================================

#ifdef _MSC_VER
#  pragma warning(disable: 4267 4311 4305)
#endif

#include <ctime>
#include <iomanip>
#include <sstream>
#include <algorithm>
// --------------------
#include <OpenMesh/Tools/Utils/Timer.hh>

#include "GLViewerCore.h"

#if !defined(M_PI)
#  define M_PI 3.1415926535897932
#endif

const double TRACKBALL_RADIUS = 0.6;


using namespace OpenMesh;


//== IMPLEMENTATION ========================================================== 

GLViewerCore::GLViewerCore(GLViewerBase &_base)
	: base_(_base)
{
	init();
}


//----------------------------------------------------------------------------

void
GLViewerCore::init(void)
{
}


//----------------------------------------------------------------------------

GLViewerCore::~GLViewerCore()
{
}


//----------------------------------------------------------------------------

void
GLViewerCore::setDefaultMaterial(void)
{
	GLfloat mat_a[] = { 0.1, 0.1, 0.1, 1.0 };
	GLfloat mat_d[] = { 0.7, 0.7, 0.5, 1.0 };
	GLfloat mat_s[] = { 1.0, 1.0, 1.0, 1.0 };
	GLfloat shine[] = { 120.0 };

	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, mat_a);
	glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, mat_d);
	glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, mat_s);
	glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, shine);
}


//----------------------------------------------------------------------------

void
GLViewerCore::setDefaultLight(void)
{
	GLfloat pos1[] = { 0.1, 0.1, -0.02, 0.0 };
	GLfloat pos2[] = { -0.1, 0.1, -0.02, 0.0 };
	GLfloat pos3[] = { 0.0, 0.0, 0.1, 0.0 };
	GLfloat col1[] = { 0.7, 0.7, 0.8, 1.0 };
	GLfloat col2[] = { 0.8, 0.7, 0.7, 1.0 };
	GLfloat col3[] = { 1.0, 1.0, 1.0, 1.0 };

	glEnable(GL_LIGHT0);
	glLightfv(GL_LIGHT0, GL_POSITION, pos1);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, col1);
	glLightfv(GL_LIGHT0, GL_SPECULAR, col1);

	glEnable(GL_LIGHT1);
	glLightfv(GL_LIGHT1, GL_POSITION, pos2);
	glLightfv(GL_LIGHT1, GL_DIFFUSE, col2);
	glLightfv(GL_LIGHT1, GL_SPECULAR, col2);

	glEnable(GL_LIGHT2);
	glLightfv(GL_LIGHT2, GL_POSITION, pos3);
	glLightfv(GL_LIGHT2, GL_DIFFUSE, col3);
	glLightfv(GL_LIGHT2, GL_SPECULAR, col3);
}


//----------------------------------------------------------------------------


void
GLViewerCore::initializeGL()
{
	// OpenGL state
	glClearColor(1.0, 1.0, 1.0, 0.0);
	glDisable(GL_DITHER);
	glEnable(GL_DEPTH_TEST);

	// Material
	setDefaultMaterial();

	// Lighting
	glLoadIdentity();
	setDefaultLight();

	// Fog
	GLfloat fogColor[4] = { 0.3, 0.3, 0.4, 1.0 };
	glFogi(GL_FOG_MODE, GL_LINEAR);
	glFogfv(GL_FOG_COLOR, fogColor);
	glFogf(GL_FOG_DENSITY, 0.35);
	glHint(GL_FOG_HINT, GL_DONT_CARE);
	glFogf(GL_FOG_START, 5.0f);
	glFogf(GL_FOG_END, 25.0f);

	// scene pos and size
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glGetDoublev(GL_MODELVIEW_MATRIX, modelview_matrix_);
	set_scene_pos(Vec3f(0.0, 0.0, 0.0), 1.0);
}


//----------------------------------------------------------------------------


void
GLViewerCore::resizeGL(int _w, int _h)
{
	update_projection_matrix();
	glViewport(0, 0, _w, _h);
	updateGL();
}


//----------------------------------------------------------------------------


void
GLViewerCore::paintGL()
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glMatrixMode(GL_PROJECTION);
	glLoadMatrixd(projection_matrix_);
	glMatrixMode(GL_MODELVIEW);
	glLoadMatrixd(modelview_matrix_);

	draw_scene(getDrawMode());
}


//----------------------------------------------------------------------------


void
GLViewerCore::draw_scene(const std::string& _draw_mode)
{
	if (_draw_mode == "Wireframe")
	{
		glDisable(GL_LIGHTING);
		glColor4f(0.0f, 0.0f, 0.0f, 1.0f);
		//glutWireTeapot(0.5);
	}

	else if (_draw_mode == "Solid Flat")
	{
		glEnable(GL_LIGHTING);
		glShadeModel(GL_FLAT);
		//glutSolidTeapot(0.5);
	}

	else if (_draw_mode == "Solid Smooth")
	{
		glEnable(GL_LIGHTING);
		glShadeModel(GL_SMOOTH);
		//glutSolidTeapot(0.5);
	}
}


//----------------------------------------------------------------------------


void
GLViewerCore::mousePressEvent(const OpenMesh::Vec2f& _new_point_2d,
bool _is_left_button, bool _is_mid_button, bool _is_right_button,
bool _is_ctrl_pressed, bool _is_alt_pressed, bool _is_shift_pressed)
{
	last_point_2D_ = _new_point_2d;
	last_point_ok_ = map_to_sphere(_new_point_2d, last_point_3D_);
}


//----------------------------------------------------------------------------


void
GLViewerCore::mouseMoveEvent(const OpenMesh::Vec2f& _new_point_2d,
bool _is_left_button, bool _is_mid_button, bool _is_right_button,
bool _is_ctrl_pressed, bool _is_alt_pressed, bool _is_shift_pressed)
{
	// Left button: rotate around center_
	// Middle button: translate object
	// Left & middle button: zoom in/out


	Vec3f  newPoint3D;
	bool   newPoint_hitSphere = map_to_sphere(_new_point_2d, newPoint3D);

	float dx = _new_point_2d[0] - last_point_2D_[0];
	float dy = _new_point_2d[1] - last_point_2D_[1];

	float w = width();
	float h = height();



	// enable GL context
	makeCurrent();


	// move in z direction
	if ((_is_left_button && _is_mid_button) ||
		(_is_left_button && _is_ctrl_pressed))
	{
		float value_y = radius_ * dy * 3.0 / h;
		translate(Vec3f(0.0, 0.0, value_y));
	}


	// move in x,y direction
	else if ((_is_mid_button) ||
		(_is_left_button && _is_alt_pressed))
	{
		float z = -(modelview_matrix_[2] * center_[0] +
			modelview_matrix_[6] * center_[1] +
			modelview_matrix_[10] * center_[2] +
			modelview_matrix_[14]) /
			(modelview_matrix_[3] * center_[0] +
			modelview_matrix_[7] * center_[1] +
			modelview_matrix_[11] * center_[2] +
			modelview_matrix_[15]);

		float aspect = w / h;
		float near_plane = 0.01 * radius_;
		float top = tan(fovy() / 2.0f*M_PI / 180.0f) * near_plane;
		float right = aspect*top;

		translate(Vec3f(2.0*dx / w*right / near_plane*z,
			-2.0*dy / h*top / near_plane*z,
			0.0f));
	}



	// rotate
	else if (_is_left_button) {

		if (last_point_ok_) {
			if ((newPoint_hitSphere = map_to_sphere(_new_point_2d, newPoint3D))) {
				Vec3f axis = last_point_3D_ % newPoint3D;
				if (axis.sqrnorm() < 1e-7) {
					axis = Vec3f(1, 0, 0);
				}
				else {
					axis.normalize();
				}
				// find the amount of rotation
				Vec3f d = last_point_3D_ - newPoint3D;
				float t = 0.5 * d.norm() / TRACKBALL_RADIUS;
				if (t < -1.0)
					t = -1.0;
				else if (t > 1.0)
					t = 1.0;
				float phi = 2.0 * asin(t);
				float angle = phi * 180.0 / M_PI;
				rotate(axis, angle);
			}
		}

	}


	// remember this point
	last_point_2D_ = _new_point_2d;
	last_point_3D_ = newPoint3D;
	last_point_ok_ = newPoint_hitSphere;

	// trigger redraw
	updateGL();
}

//----------------------------------------------------------------------------


void
GLViewerCore::mouseReleaseEvent(const OpenMesh::Vec2f& _new_point_2d,
bool _is_left_button, bool _is_mid_button, bool _is_right_button,
bool _is_ctrl_pressed, bool _is_alt_pressed, bool _is_shift_pressed)
{
	last_point_ok_ = false;
}


//-----------------------------------------------------------------------------


void GLViewerCore::wheelEvent(const float _new_delta,
	bool _is_ctrl_pressed, bool _is_alt_pressed, bool _is_shift_pressed)
{
	// Use the mouse wheel to zoom in/out

	float d = -_new_delta / 120.0 * 0.2 * radius_;
	translate(Vec3f(0.0, 0.0, d));
	updateGL();
}


//----------------------------------------------------------------------------


void GLViewerCore::keyPressEvent(const char _new_key,
	bool _is_ctrl_pressed, bool _is_alt_pressed, bool _is_shift_pressed)
{
	//std::time_t t = std::time(nullptr);
	std::stringstream time;
	//time << std::put_time(std::localtime(&t), "%Y_%b_%d_%H_%M_%S");
        time << "snapshot";

	switch (_new_key)
	{
	case 'p':
		snapshot(time.str());
		break;

	case 'h':
		std::cout << "Keys:\n";
		std::cout << "  Print\tMake snapshot\n";
		std::cout << "  C\tenable/disable back face culling\n";
		std::cout << "  F\tenable/disable fog\n";
		std::cout << "  I\tDisplay information\n";
		std::cout << "  N\tenable/disable display of vertex normals\n";
		std::cout << "  Shift N\tenable/disable display of face normals\n";
		std::cout << "  Shift P\tperformance check\n";
		break;

	case 'c':
		if (glIsEnabled(GL_CULL_FACE))
		{
			glDisable(GL_CULL_FACE);
			std::cout << "Back face culling: disabled\n";
		}
		else
		{
			glEnable(GL_CULL_FACE);
			std::cout << "Back face culling: enabled\n";
		}
		updateGL();
		break;

	case 'f':
		if (glIsEnabled(GL_FOG))
		{
			glDisable(GL_FOG);
			std::cout << "Fog: disabled\n";
		}
		else
		{
			glEnable(GL_FOG);
			std::cout << "Fog: enabled\n";
		}
		updateGL();
		break;

	case 'i':
		std::cout << "Scene radius: " << radius_ << std::endl;
		std::cout << "Scene center: " << center_ << std::endl;
		break;

	//case 'p':
	//	if (_is_shift_pressed)
	//	{
	//		double fps = performance();
	//		std::cout << "fps: "
	//			<< std::setiosflags(std::ios_base::fixed)
	//			<< fps << std::endl;
	//	}
		break;
	}
}

//----------------------------------------------------------------------------


void
GLViewerCore::translate(const OpenMesh::Vec3f& _trans)
{
	// Translate the object by _trans
	// Update modelview_matrix_
	makeCurrent();
	glLoadIdentity();
	glTranslated(_trans[0], _trans[1], _trans[2]);
	glMultMatrixd(modelview_matrix_);
	glGetDoublev(GL_MODELVIEW_MATRIX, modelview_matrix_);
}


//----------------------------------------------------------------------------


void
GLViewerCore::rotate(const OpenMesh::Vec3f& _axis, float _angle)
{
	// Rotate around center center_, axis _axis, by angle _angle
	// Update modelview_matrix_

	Vec3f t(modelview_matrix_[0] * center_[0] +
		modelview_matrix_[4] * center_[1] +
		modelview_matrix_[8] * center_[2] +
		modelview_matrix_[12],
		modelview_matrix_[1] * center_[0] +
		modelview_matrix_[5] * center_[1] +
		modelview_matrix_[9] * center_[2] +
		modelview_matrix_[13],
		modelview_matrix_[2] * center_[0] +
		modelview_matrix_[6] * center_[1] +
		modelview_matrix_[10] * center_[2] +
		modelview_matrix_[14]);

	makeCurrent();
	glLoadIdentity();
	glTranslatef(t[0], t[1], t[2]);
	glRotated(_angle, _axis[0], _axis[1], _axis[2]);
	glTranslatef(-t[0], -t[1], -t[2]);
	glMultMatrixd(modelview_matrix_);
	glGetDoublev(GL_MODELVIEW_MATRIX, modelview_matrix_);
}


//----------------------------------------------------------------------------


bool
GLViewerCore::map_to_sphere(const OpenMesh::Vec2f& _v2D, OpenMesh::Vec3f& _v3D)
{
	// This is actually doing the Sphere/Hyperbolic sheet hybrid thing,
	// based on Ken Shoemake's ArcBall in Graphics Gems IV, 1993.
	double x = (2.0*_v2D[0] - width()) / width();
	double y = -(2.0*_v2D[1] - height()) / height();
	double xval = x;
	double yval = y;
	double x2y2 = xval*xval + yval*yval;

	const double rsqr = TRACKBALL_RADIUS*TRACKBALL_RADIUS;
	_v3D[0] = xval;
	_v3D[1] = yval;
	if (x2y2 < 0.5*rsqr) {
		_v3D[2] = sqrt(rsqr - x2y2);
	}
	else {
		_v3D[2] = 0.5*rsqr / sqrt(x2y2);
	}

	return true;
}


//----------------------------------------------------------------------------


void
GLViewerCore::update_projection_matrix()
{
	makeCurrent();
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(45.0, (GLfloat)width() / (GLfloat)height(),
		0.01*radius_, 100.0*radius_);
	glGetDoublev(GL_PROJECTION_MATRIX, projection_matrix_);
	glMatrixMode(GL_MODELVIEW);
}


//----------------------------------------------------------------------------


void
GLViewerCore::view_all()
{
	translate(Vec3f(-(modelview_matrix_[0] * center_[0] +
		modelview_matrix_[4] * center_[1] +
		modelview_matrix_[8] * center_[2] +
		modelview_matrix_[12]),
		-(modelview_matrix_[1] * center_[0] +
		modelview_matrix_[5] * center_[1] +
		modelview_matrix_[9] * center_[2] +
		modelview_matrix_[13]),
		-(modelview_matrix_[2] * center_[0] +
		modelview_matrix_[6] * center_[1] +
		modelview_matrix_[10] * center_[2] +
		modelview_matrix_[14] +
		3.0*radius_)));
}


//----------------------------------------------------------------------------


void
GLViewerCore::set_scene_pos(const OpenMesh::Vec3f& _cog, float _radius)
{
	center_ = _cog;
	radius_ = _radius;
	glFogf(GL_FOG_START, 1.5*_radius);
	glFogf(GL_FOG_END, 3.0*_radius);

	update_projection_matrix();
	view_all();
}


//=============================================================================


void GLViewerCore::set_modelview_matrix(GLdouble *matrix)
{
	for (unsigned int i = 0; i < 16; i++)
		modelview_matrix_[i] = matrix[i];

	makeCurrent();
	glLoadIdentity();
	glMultMatrixd(modelview_matrix_);

	updateGL();
}


//=============================================================================
