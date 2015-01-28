// glutmeshviewer.cpp

//-----------------------------------------------------------------------------
// Includes
#include <GL/glut.h>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <string>

#include "MeshViewerCore.h"


//-----------------------------------------------------------------------------
// Global variables
static bool g_bLeftButtonDown = false;
unsigned int g_DrawMode = 0;
std::vector<std::string> g_DrawModeNames;


//-----------------------------------------------------------------------------
// Function definitions
void initializeGL(void);
void display(void);
void reshape(GLint width, GLint height);
void mouseButton(int button, int state, int x, int y);
void mouseMotion(int x, int y);
void Keyboard(unsigned char key, int x, int y);
void snapshot(const std::string filename);


//-----------------------------------------------------------------------------
// Class definition
class GlutMeshViewerBase : public GLViewerBase
{
	virtual float GLViewerBaseWidth() { return glutGet(GLUT_WINDOW_WIDTH); }

	virtual float GLViewerBaseHeight() { return glutGet(GLUT_WINDOW_HEIGHT); }

	virtual void GLViewerBasemakeCurrent() {
		int window = glutGetWindow();
		glutSetWindow(window);
	}

	virtual void GLViewerBaseSwapBuffers() { glutSwapBuffers(); }

	virtual void GLViewerBaseUpdateGL() { display(); }

	virtual void GLViewerBaseSnapshot(const std::string &filename) { snapshot(filename); }

	virtual void GLViewerBaseSetWindowTitle(const std::string &title) { glutSetWindowTitle(title.c_str()); }

	virtual std::string GLViewerBaseGetDrawMode()
	{
		assert(g_DrawMode < g_DrawModeNames.size());
		return g_DrawModeNames[g_DrawMode];
	};

	virtual void GLViewerBaseSetDrawMode(const std::string &_draw_mode)
	{
		for (unsigned int i = 0; i < g_DrawModeNames.size(); ++i)
			if (g_DrawModeNames[i] == _draw_mode)
			{
				g_DrawMode = i; return;
			}
		assert(false);
	}
};

GlutMeshViewerBase g_MeshViewerBase;
MeshViewerCore g_MeshViewer(g_MeshViewerBase);


//-----------------------------------------------------------------------------
// Implementation
void initializeGL(void)
{
	g_MeshViewer.initializeGL();
}

void display(void)
{
	g_MeshViewer.paintGL();
	glutSwapBuffers();
}

void reshape(GLint width, GLint height)
{
	g_MeshViewer.resizeGL(width, height);
}

void mouseButton(int button, int state, int x, int y)
{
	OpenMesh::Vec2f new_point_2d(x, y);

	bool is_left_button = (button == GLUT_LEFT_BUTTON);
	bool is_mid_button = (button == GLUT_MIDDLE_BUTTON);
	bool is_right_button = (button == GLUT_RIGHT_BUTTON);

	bool is_ctrl_pressed = (glutGetModifiers() == GLUT_ACTIVE_CTRL);
	bool is_alt_pressed = (glutGetModifiers() == GLUT_ACTIVE_ALT);
	bool is_shift_pressed = (glutGetModifiers() == GLUT_ACTIVE_SHIFT);

	if (state == GLUT_DOWN)
	{
		g_bLeftButtonDown = true;
		g_MeshViewer.mousePressEvent(new_point_2d,
			is_left_button, is_mid_button, is_right_button,
			is_ctrl_pressed, is_alt_pressed, is_shift_pressed);
	}
	else if (state == GLUT_UP)
	{
		g_bLeftButtonDown = false;
		g_MeshViewer.mouseReleaseEvent(new_point_2d,
			is_left_button, is_mid_button, is_right_button,
			is_ctrl_pressed, is_alt_pressed, is_shift_pressed);
	}
}

void mouseMotion(int x, int y)
{
	if (g_bLeftButtonDown)
	{
		OpenMesh::Vec2f new_point_2d(x, y);

		bool is_left_button = true;
		bool is_mid_button = false;
		bool is_right_button = false;

		bool is_ctrl_pressed = (glutGetModifiers() == GLUT_ACTIVE_CTRL);
		bool is_alt_pressed = (glutGetModifiers() == GLUT_ACTIVE_ALT);
		bool is_shift_pressed = (glutGetModifiers() == GLUT_ACTIVE_SHIFT);

		g_MeshViewer.mouseMoveEvent(new_point_2d,
			is_left_button, is_mid_button, is_right_button,
			is_ctrl_pressed, is_alt_pressed, is_shift_pressed);

		glutPostRedisplay();
	}
}

void Keyboard(unsigned char key, int x, int y)
{
	switch (key)
	{
	case 27:             // ESCAPE key
		exit(0);
	}

	bool is_ctrl_pressed = (glutGetModifiers() == GLUT_ACTIVE_CTRL);
	bool is_alt_pressed = (glutGetModifiers() == GLUT_ACTIVE_ALT);
	bool is_shift_pressed = (glutGetModifiers() == GLUT_ACTIVE_SHIFT);

	g_MeshViewer.keyPressEvent(key,
		is_ctrl_pressed, is_alt_pressed, is_shift_pressed);
}

void selectFromMenu(int command_id)
{
	g_DrawMode = command_id;
	assert(g_DrawMode < g_DrawModeNames.size());

	// Almost any menu selection requires a redraw
	glutPostRedisplay();
}

int buildPopupMenu(void)
{
	g_DrawModeNames.push_back("Wireframe");
	g_DrawModeNames.push_back("Solid Flat");
	g_DrawModeNames.push_back("Solid Smooth");

	g_DrawMode = g_DrawModeNames.size() - 1;

	int menu = glutCreateMenu(selectFromMenu);
	for (int i = 0; i < g_DrawModeNames.size(); ++i)
		glutAddMenuEntry(g_DrawModeNames[i].c_str(), i);

	return menu;
}

void snapshot(const std::string filename)
{
	GLenum buffer(GL_BACK);
	const int offset = 8;
	int w = glutGet(GLUT_WINDOW_WIDTH) - offset;
	int h = glutGet(GLUT_WINDOW_HEIGHT) - offset;

	std::vector<GLubyte> fbuffer(3 * w * h);

	glReadBuffer(buffer);
	glPixelStorei(GL_PACK_ALIGNMENT, 1);
	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
	display();
	glReadPixels(0, 0, w, h, GL_RGB, GL_UNSIGNED_BYTE, &fbuffer[0]);


	std::vector<GLubyte> image(3 * w*h);
	unsigned int x, y, fbuffer_offset, image_offset;

	for (y = 0; y < h; ++y) {
		for (x = 0; x < w; ++x) {
			fbuffer_offset = 3 * (y * w + x);
			image_offset = 3 * ((h - y - 1)* w + x);
			for (int i = 0; i < 3; ++i)
				image[image_offset + i] = fbuffer[fbuffer_offset + i];
		}
	}

	std::string name = filename;
	name.append(".ppm");
	FILE *ppmFile = fopen(name.c_str(), "wb");

	fprintf(ppmFile, "P6\n");
	fprintf(ppmFile, "%d %d\n", w, h);
	fprintf(ppmFile, "255\n");
	fwrite(&image[0], 1, w * h * 3, ppmFile);

	fclose(ppmFile);
}

int main(int argc, char** argv)
{
	// GLUT Window Initialization:
	glutInit(&argc, argv);
	glutInitWindowSize(640, 480);
	glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
	glutCreateWindow("GlutViewer");

	// Initialize OpenGL graphics state
	initializeGL();

	// Register callbacks
	glutDisplayFunc(display);
	glutReshapeFunc(reshape);
	glutMouseFunc(mouseButton);
	glutMotionFunc(mouseMotion);
	glutKeyboardFunc(Keyboard);

	// Create our popup menu
	buildPopupMenu();
	glutAttachMenu(GLUT_RIGHT_BUTTON);

	// To do
	//const char *mesh_filename = "mesh.off";
	//g_MeshViewer.open_mesh(mesh_filename);
	g_MeshViewer.parse_arguments();

	// Turn the flow of control over to GLUT
	glutMainLoop();
	return 0;
}
