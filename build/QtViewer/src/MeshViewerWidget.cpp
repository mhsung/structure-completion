#include "MeshViewerWidget.h"

// --------------------
#include <qapplication.h>
#include <qnamespace.h>
#include <QMouseEvent>
// --------------------


MeshViewerWidget::MeshViewerWidget(QWidget* parent /*= 0*/)
	: QGLViewerWidget(parent)
	, mesh_viewer_core_(*this)
{
	initializeGL();

    add_draw_mode("Points");
    add_draw_mode("Hidden-Line");
#if defined(OM_USE_OSG) && OM_USE_OSG
    add_draw_mode("OpenSG Indices");
#endif
}

void MeshViewerWidget::initializeGL()
{
	mesh_viewer_core_.initializeGL();
}

void MeshViewerWidget::resizeGL(int w, int h)
{
	mesh_viewer_core_.resizeGL(w, h);
}

void MeshViewerWidget::paintGL()
{
	mesh_viewer_core_.paintGL();
}

void MeshViewerWidget::mousePressEvent(QMouseEvent* _event)
{
	//
	QGLViewerWidget::mousePressEvent(_event);
	//

	if (_event->button() != Qt::RightButton)
	{
		OpenMesh::Vec2f new_point_2d;
		new_point_2d[0] = _event->pos().x();
		new_point_2d[1] = _event->pos().y();

		bool is_left_button = (_event->buttons() == Qt::LeftButton);
		bool is_mid_button = (_event->buttons() == Qt::MidButton);
		bool is_right_button = (_event->buttons() == Qt::RightButton);

		bool is_ctrl_pressed = (_event->modifiers() & Qt::ControlModifier);
		bool is_alt_pressed = (_event->modifiers() & Qt::AltModifier);
		bool is_shift_pressed = (_event->modifiers() & Qt::ShiftModifier);

		mesh_viewer_core_.mousePressEvent(new_point_2d,
			is_left_button, is_mid_button, is_right_button,
			is_ctrl_pressed, is_alt_pressed, is_shift_pressed);
	}
}

void MeshViewerWidget::mouseReleaseEvent(QMouseEvent* _event)
{
	OpenMesh::Vec2f new_point_2d;
	new_point_2d[0] = _event->pos().x();
	new_point_2d[1] = _event->pos().y();

	bool is_left_button = (_event->buttons() == Qt::LeftButton);
	bool is_mid_button = (_event->buttons() == Qt::MidButton);
	bool is_right_button = (_event->buttons() == Qt::RightButton);

	bool is_ctrl_pressed = (_event->modifiers() & Qt::ControlModifier);
	bool is_alt_pressed = (_event->modifiers() & Qt::AltModifier);
	bool is_shift_pressed = (_event->modifiers() & Qt::ShiftModifier);

	mesh_viewer_core_.mouseReleaseEvent(new_point_2d,
		is_left_button, is_mid_button, is_right_button,
		is_ctrl_pressed, is_alt_pressed, is_shift_pressed);
}

void MeshViewerWidget::mouseMoveEvent(QMouseEvent* _event)
{
	OpenMesh::Vec2f new_point_2d;
	new_point_2d[0] = _event->pos().x();
	new_point_2d[1] = _event->pos().y();

	bool is_left_button = (_event->buttons() == Qt::LeftButton);
	bool is_mid_button = (_event->buttons() == Qt::MidButton);
	bool is_right_button = (_event->buttons() == Qt::RightButton);

	bool is_ctrl_pressed = (_event->modifiers() & Qt::ControlModifier);
	bool is_alt_pressed = (_event->modifiers() & Qt::AltModifier);
	bool is_shift_pressed = (_event->modifiers() & Qt::ShiftModifier);

	mesh_viewer_core_.mouseMoveEvent(new_point_2d,
		is_left_button, is_mid_button, is_right_button,
		is_ctrl_pressed, is_alt_pressed, is_shift_pressed);
}

void MeshViewerWidget::wheelEvent(QWheelEvent* _event)
{
	float new_delta;
	new_delta = _event->delta();

	bool is_ctrl_pressed = (_event->modifiers() & Qt::ControlModifier);
	bool is_alt_pressed = (_event->modifiers() & Qt::AltModifier);
	bool is_shift_pressed = (_event->modifiers() & Qt::ShiftModifier);

	mesh_viewer_core_.wheelEvent(new_delta,
		is_ctrl_pressed, is_alt_pressed, is_shift_pressed);

	//
	_event->accept();
	//
}

void MeshViewerWidget::keyPressEvent(QKeyEvent* _event)
{
	//
	switch (_event->key())
	{
	case Qt::Key_Q:
	case Qt::Key_Escape:
		qApp->quit();
	}
	//

	char new_key;
	new_key = _event->text().toStdString().c_str()[0];

	bool is_ctrl_pressed = (_event->modifiers() & Qt::ControlModifier);
	bool is_alt_pressed = (_event->modifiers() & Qt::AltModifier);
	bool is_shift_pressed = (_event->modifiers() & Qt::ShiftModifier);

	mesh_viewer_core_.keyPressEvent(new_key,
		is_ctrl_pressed, is_alt_pressed, is_shift_pressed);

	//
	_event->ignore();
	//
}
