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

#include <iomanip>
#include <sstream>
#include <algorithm>
// --------------------
#ifdef ARCH_DARWIN
#  include <glut.h>
#else
#  include <GL/glut.h>
#endif
// --------------------
#include <qnamespace.h>
#include <qapplication.h>
#include <qmenu.h>
#include <qcursor.h>
#include <qimage.h>
#include <qdatetime.h>
#include <QMouseEvent>
// --------------------
#include <OpenMesh/Tools/Utils/Timer.hh>

#include "QGLViewerWidget.h"

#if !defined(M_PI)
#  define M_PI 3.1415926535897932
#endif

const double TRACKBALL_RADIUS = 0.6;


using namespace Qt;
using namespace OpenMesh;


//== IMPLEMENTATION ========================================================== 

std::string QGLViewerWidget::nomode_ = "";

//----------------------------------------------------------------------------

QGLViewerWidget::QGLViewerWidget(QWidget* _parent)
	: QGLWidget(_parent)
{
	init();
}

//----------------------------------------------------------------------------

QGLViewerWidget::
QGLViewerWidget(QGLFormat& _fmt, QWidget* _parent)
: QGLWidget(_fmt, _parent)
{
	init();
}


//----------------------------------------------------------------------------

void
QGLViewerWidget::init(void)
{
	// qt stuff
	setAttribute(Qt::WA_NoSystemBackground, true);
	setFocusPolicy(Qt::StrongFocus);
	setAcceptDrops(true);
	setCursor(PointingHandCursor);


	// popup menu

	popup_menu_ = new QMenu(this);
	draw_modes_group_ = new QActionGroup(this);

	connect(draw_modes_group_, SIGNAL(triggered(QAction*)),
		this, SLOT(slotDrawMode(QAction*)));


	// init draw modes
	n_draw_modes_ = 0;
	//draw_mode_ = 3;
	QAction *a;
	a = add_draw_mode("Wireframe");
	a->setShortcut(QKeySequence(Key_W));
	add_draw_mode("Solid Flat");
	a = add_draw_mode("Solid Smooth");
	a->setShortcut(QKeySequence(Key_S));
	a->setChecked(true);

	slotDrawMode(a);
}


//----------------------------------------------------------------------------

QGLViewerWidget::~QGLViewerWidget()
{
}


//----------------------------------------------------------------------------


void
QGLViewerWidget::initializeGL()
{
}


//----------------------------------------------------------------------------


void
QGLViewerWidget::resizeGL(int _w, int _h)
{
}


//----------------------------------------------------------------------------


void
QGLViewerWidget::paintGL()
{
}


//----------------------------------------------------------------------------


void
QGLViewerWidget::mousePressEvent(QMouseEvent* _event)
{
	// popup menu
	if (_event->button() == RightButton)
	{
		popup_menu_->exec(QCursor::pos());
	}
}


//----------------------------------------------------------------------------


void
QGLViewerWidget::mouseMoveEvent(QMouseEvent* _event)
{
	// enable GL context
	makeCurrent();

	// trigger redraw
	updateGL();
}


//----------------------------------------------------------------------------


void
QGLViewerWidget::mouseReleaseEvent(QMouseEvent* /* _event */)
{
}


//-----------------------------------------------------------------------------


void QGLViewerWidget::wheelEvent(QWheelEvent* _event)
{
	_event->accept();
}


//----------------------------------------------------------------------------


void QGLViewerWidget::keyPressEvent(QKeyEvent* _event)
{
	switch (_event->key())
	{
	case Key_Q:
	case Key_Escape:
		qApp->quit();
	}

	_event->ignore();
}


//----------------------------------------------------------------------------


QAction*
QGLViewerWidget::add_draw_mode(const std::string& _s)
{
	++n_draw_modes_;
	draw_mode_names_.push_back(_s);

	QActionGroup *grp = draw_modes_group_;
	QAction* act = new QAction(tr(_s.c_str()), this);
	act->setCheckable(true);
	act->setData(n_draw_modes_);

	grp->addAction(act);
	popup_menu_->addAction(act);
	addAction(act, _s.c_str());

	return act;
}

void QGLViewerWidget::addAction(QAction* act, const char * name)
{
	names_to_actions[name] = act;
	Super::addAction(act);
}
void QGLViewerWidget::removeAction(QAction* act)
{
	ActionMap::iterator it = names_to_actions.begin(), e = names_to_actions.end();
	ActionMap::iterator found = e;
	for (; it != e; ++it) {
		if (it->second == act) {
			found = it;
			break;
		}
	}
	if (found != e) {
		names_to_actions.erase(found);
	}
	popup_menu_->removeAction(act);
	draw_modes_group_->removeAction(act);
	Super::removeAction(act);
}

void QGLViewerWidget::removeAction(const char* name)
{
	QString namestr = QString(name);
	ActionMap::iterator e = names_to_actions.end();

	ActionMap::iterator found = names_to_actions.find(namestr);
	if (found != e) {
		removeAction(found->second);
	}
}

QAction* QGLViewerWidget::findAction(const char* name)
{
	QString namestr = QString(name);
	ActionMap::iterator e = names_to_actions.end();

	ActionMap::iterator found = names_to_actions.find(namestr);
	if (found != e) {
		return found->second;
	}
	return 0;
}

//----------------------------------------------------------------------------


void
QGLViewerWidget::del_draw_mode(const std::string& _s)
{
	QString cmp = _s.c_str();
	QList<QAction*> actions_ = popup_menu_->actions();
	QList<QAction*>::iterator it = actions_.begin(), e = actions_.end();
	for (; it != e; ++it) {
		if ((*it)->text() == cmp) { break; }
	}

#if _DEBUG
	assert(it != e);
#else
	if (it == e)
		return;
#endif

	popup_menu_->removeAction(*it);
	//QActionGroup *grp = draw_modes_group_;

}


//----------------------------------------------------------------------------


void
QGLViewerWidget::slotDrawMode(QAction* _mode)
{
	// save draw mode
	draw_mode_ = _mode->data().toInt();
	updateGL();

	// check selected draw mode
	//popup_menu_->setItemChecked(draw_mode_, true);
}


//----------------------------------------------------------------------------

/*
double
QGLViewerWidget::performance()
{
	setCursor(Qt::WaitCursor);

	double fps(0.0);

	makeCurrent();
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();

	OpenMesh::Utils::Timer timer;

	unsigned int  frames = 60;
	const float   angle = 360.0 / (float)frames;
	unsigned int  i;
	Vec3f         axis;

	glFinish();

	timer.start();
	for (i = 0, axis = Vec3f(1, 0, 0); i<frames; ++i)
	{
		rotate(axis, angle); paintGL(); swapBuffers();
	}
	timer.stop();

	qApp->processEvents();

	timer.cont();
	for (i = 0, axis = Vec3f(0, 1, 0); i<frames; ++i)
	{
		rotate(axis, angle); paintGL(); swapBuffers();
	}
	timer.stop();

	qApp->processEvents();

	timer.cont();
	for (i = 0, axis = Vec3f(0, 0, 1); i<frames; ++i)
	{
		rotate(axis, angle); paintGL(); swapBuffers();
	}
	timer.stop();

	glFinish();
	timer.stop();

	glPopMatrix();
	updateGL();

	fps = ((3.0 * frames) / timer.seconds());

	setCursor(PointingHandCursor);

	return fps;
}
*/

void
QGLViewerWidget::slotSnapshot(const std::string _filename)
{
	QImage image;
	size_t w(width()), h(height());
	GLenum buffer(GL_BACK);

	try
	{
		image = QImage(w, h, QImage::Format_RGB32);

		std::vector<GLubyte> fbuffer(3 * w*h);

		qApp->processEvents();
		makeCurrent();
		updateGL();
		glFinish();

		glReadBuffer(buffer);
		glPixelStorei(GL_PACK_ALIGNMENT, 1);
		glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
		paintGL();
		glReadPixels(0, 0, w, h, GL_RGB, GL_UNSIGNED_BYTE, &fbuffer[0]);

		unsigned int x, y, offset;

		for (y = 0; y<h; ++y) {
			for (x = 0; x<w; ++x) {
				offset = 3 * (y*w + x);
				image.setPixel(x, h - y - 1, qRgb(fbuffer[offset],
					fbuffer[offset + 1],
					fbuffer[offset + 2]));
			}
		}


//		QString name = "snapshot-";
//#if defined(_MSC_VER)
//		{
//			std::stringstream s;
//			QDateTime         dt = QDateTime::currentDateTime();
//			s << dt.date().year()
//				<< std::setw(2) << std::setfill('0') << dt.date().month()
//				<< std::setw(2) << std::setfill('0') << dt.date().day()
//				<< std::setw(2) << std::setfill('0') << dt.time().hour()
//				<< std::setw(2) << std::setfill('0') << dt.time().minute()
//				<< std::setw(2) << std::setfill('0') << dt.time().second();
//			name += QString(s.str().c_str());
//		}
//#else
//		name += QDateTime::currentDateTime().toString("yyMMddhhmmss");
//#endif
		QString name = _filename.c_str();
		name += ".png";

		image.save(name, "PNG");
	}
	catch (std::bad_alloc&)
	{
		qWarning("Mem Alloc Error");
	}

}


//=============================================================================
