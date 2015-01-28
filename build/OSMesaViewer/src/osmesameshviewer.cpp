// glutmeshviewer.cpp

//-----------------------------------------------------------------------------
// Includes
#include "GL/osmesa.h"
#include "GL/glu.h"
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <string>

#include <gflags/gflags.h>
#include <glog/logging.h>
#include <QImage>

#include "MeshViewerCore.h"


//-----------------------------------------------------------------------------
// Global variables
OSMesaContext g_Context;
void *g_Buffer;
static int g_Width = 640;
static int g_Height = 480;
std::string g_DrawMode = "Solid Smooth";


//-----------------------------------------------------------------------------
// Function definitions
void display(void);
void snapshot(const std::string filename);


//-----------------------------------------------------------------------------
// Class definition
class GlutMeshViewerBase : public GLViewerBase
{
	virtual float GLViewerBaseWidth() { return g_Width; }

	virtual float GLViewerBaseHeight() { return g_Height; }

	virtual void GLViewerBasemakeCurrent() { return; }

	virtual void GLViewerBaseSwapBuffers() { glFinish(); }

	virtual void GLViewerBaseUpdateGL() { display(); }

	virtual void GLViewerBaseSnapshot(const std::string &filename) { snapshot(filename); }

	virtual void GLViewerBaseSetWindowTitle(const std::string &title) { }

	virtual std::string GLViewerBaseGetDrawMode() { return g_DrawMode; };

	virtual void GLViewerBaseSetDrawMode(const std::string &_draw_mode) { g_DrawMode = _draw_mode; }
};

GlutMeshViewerBase g_MeshViewerBase;
MeshViewerCore g_MeshViewer(g_MeshViewerBase);


//-----------------------------------------------------------------------------
// Implementation

void Sphere(float radius, int slices, int stacks)
{
	GLUquadric *q = gluNewQuadric();
	gluQuadricNormals(q, GLU_SMOOTH);
	gluSphere(q, radius, slices, stacks);
	gluDeleteQuadric(q);
}


static void
Cone(float base, float height, int slices, int stacks)
{
	GLUquadric *q = gluNewQuadric();
	gluQuadricDrawStyle(q, GLU_FILL);
	gluQuadricNormals(q, GLU_SMOOTH);
	gluCylinder(q, base, 0.0, height, slices, stacks);
	gluDeleteQuadric(q);
}


static void
Torus(float innerRadius, float outerRadius, int sides, int rings)
{
	/* from GLUT... */
	int i, j;
	GLfloat theta, phi, theta1;
	GLfloat cosTheta, sinTheta;
	GLfloat cosTheta1, sinTheta1;
	const GLfloat ringDelta = 2.0 * M_PI / rings;
	const GLfloat sideDelta = 2.0 * M_PI / sides;

	theta = 0.0;
	cosTheta = 1.0;
	sinTheta = 0.0;
	for (i = rings - 1; i >= 0; i--) {
		theta1 = theta + ringDelta;
		cosTheta1 = cos(theta1);
		sinTheta1 = sin(theta1);
		glBegin(GL_QUAD_STRIP);
		phi = 0.0;
		for (j = sides; j >= 0; j--) {
			GLfloat cosPhi, sinPhi, dist;

			phi += sideDelta;
			cosPhi = cos(phi);
			sinPhi = sin(phi);
			dist = outerRadius + innerRadius * cosPhi;

			glNormal3f(cosTheta1 * cosPhi, -sinTheta1 * cosPhi, sinPhi);
			glVertex3f(cosTheta1 * dist, -sinTheta1 * dist, innerRadius * sinPhi);
			glNormal3f(cosTheta * cosPhi, -sinTheta * cosPhi, sinPhi);
			glVertex3f(cosTheta * dist, -sinTheta * dist, innerRadius * sinPhi);
		}
		glEnd();
		theta = theta1;
		cosTheta = cosTheta1;
		sinTheta = sinTheta1;
	}
}


void render_image(void)
{
	GLfloat light_ambient[] = { 0.0, 0.0, 0.0, 1.0 };
	GLfloat light_diffuse[] = { 1.0, 1.0, 1.0, 1.0 };
	GLfloat light_specular[] = { 1.0, 1.0, 1.0, 1.0 };
	GLfloat light_position[] = { 1.0, 1.0, 1.0, 0.0 };
	GLfloat red_mat[] = { 1.0, 0.2, 0.2, 1.0 };
	GLfloat green_mat[] = { 0.2, 1.0, 0.2, 1.0 };
	GLfloat blue_mat[] = { 0.2, 0.2, 1.0, 1.0 };


	glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse);
	glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular);
	glLightfv(GL_LIGHT0, GL_POSITION, light_position);

	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glEnable(GL_DEPTH_TEST);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(-2.5, 2.5, -2.5, 2.5, -10.0, 10.0);
	glMatrixMode(GL_MODELVIEW);

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glPushMatrix();
	glRotatef(20.0, 1.0, 0.0, 0.0);

	glPushMatrix();
	glTranslatef(-0.75, 0.5, 0.0);
	glRotatef(90.0, 1.0, 0.0, 0.0);
	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, red_mat);
	Torus(0.275, 0.85, 20, 20);
	glPopMatrix();

	glPushMatrix();
	glTranslatef(-0.75, -0.5, 0.0);
	glRotatef(270.0, 1.0, 0.0, 0.0);
	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, green_mat);
	Cone(1.0, 2.0, 16, 1);
	glPopMatrix();

	glPushMatrix();
	glTranslatef(0.75, 0.0, -1.0);
	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, blue_mat);
	Sphere(1.0, 20, 20);
	glPopMatrix();

	glPopMatrix();

	/* This is very important!!!
	*     * Make sure buffered commands are finished!!!
	*         */
	glFinish();
}

int osmesaInitWindowSize(int w, int h)
{
	g_Width = w;
	g_Height = h;

	g_Context = OSMesaCreateContextExt(OSMESA_RGBA, 16, 0, 0, NULL);
	if (!g_Context) {
		printf("OSMesaCreateContext failed!\n");
		return -1;
	}

	g_Buffer = malloc(g_Width * g_Height * 4 * sizeof(GLubyte));
	if (!g_Buffer) {
		printf("Alloc image buffer failed!\n");
		return -1;
	}

	if (!OSMesaMakeCurrent(g_Context, g_Buffer, GL_UNSIGNED_BYTE, g_Width, g_Height)) {
		printf("OSMesaMakeCurrent failed!\n");
		return -1;
	}


   {
	   int z, s, a;
	   glGetIntegerv(GL_DEPTH_BITS, &z);
	   glGetIntegerv(GL_STENCIL_BITS, &s);
	   glGetIntegerv(GL_ACCUM_RED_BITS, &a);
	   printf("Depth=%d Stencil=%d Accum=%d\n", z, s, a);
   }

   g_MeshViewer.initializeGL();
   return 0;
}

void osmesaFreeContext()
{
	free(g_Buffer);

	OSMesaDestroyContext(g_Context);
}

void display(void)
{
	g_MeshViewer.paintGL();
	glFinish();
}

void snapshot(const std::string filename)
{
	//// Save binary image
	//const int binary = 1;

	//FILE *f = fopen((filename + std::string(".ppm")).c_str(), "w");
	//if (f) {
	//	int i, x, y;
	//	const GLubyte *ptr = static_cast<GLubyte *>(g_Buffer);
	//	if (binary) {
	//		fprintf(f, "P6\n");
	//		fprintf(f, "# ppm-file created by osdemo.c\n");
	//		fprintf(f, "%i %i\n", g_Width, g_Height);
	//		fprintf(f, "255\n");
	//		fclose(f);
	//		f = fopen((filename + std::string(".ppm")).c_str(), "ab");  /* reopen in binary append mode */
	//		for (y = g_Height - 1; y >= 0; y--) {
	//			for (x = 0; x < g_Width; x++) {
	//				i = (y*g_Width + x) * 4;
	//				fputc(ptr[i], f);   /* write red */
	//				fputc(ptr[i + 1], f); /* write green */
	//				fputc(ptr[i + 2], f); /* write blue */
	//			}
	//		}
	//	}
	//	else {
	//		/*ASCII*/
	//		int counter = 0;
	//		fprintf(f, "P3\n");
	//		fprintf(f, "# ascii ppm file created by osdemo.c\n");
	//		fprintf(f, "%i %i\n", g_Width, g_Height);
	//		fprintf(f, "255\n");
	//		for (y = g_Height - 1; y >= 0; y--) {
	//			for (x = 0; x < g_Width; x++) {
	//				i = (y*g_Width + x) * 4;
	//				fprintf(f, " %3d %3d %3d", ptr[i], ptr[i + 1], ptr[i + 2]);
	//				counter++;
	//				if (counter % 5 == 0)
	//					fprintf(f, "\n");
	//			}
	//		}
	//	}
	//	fclose(f);
	//}

	// Save PNG images using Qt library
	try
	{
		QImage image;
		size_t w(g_Width), h(g_Height);
		image = QImage(w, h, QImage::Format_RGB32);

                std::vector<GLubyte> fbuffer(4 * w*h);
                memcpy(&fbuffer[0], g_Buffer, 4 * w*h * sizeof(GLubyte));

		unsigned int x, y, offset;
		for (y = 0; y < h; ++y) {
			for (x = 0; x < w; ++x) {
				offset = 4 * (y*w + x);
				image.setPixel(x, h - y - 1, qRgb(fbuffer[offset],
					fbuffer[offset + 1],
					fbuffer[offset + 2]));
			}
		}

		image.save((filename + std::string(".png")).c_str(), "PNG");
	}
	catch (std::bad_alloc&)
	{
		qWarning("Mem Alloc Error");
	}
}

int main(int argc, char** argv)
{
	gflags::ParseCommandLineFlags(&argc, &argv, true);
	google::InitGoogleLogging(argv[0]);

	// GLUT Window Initialization:
	osmesaInitWindowSize(640, 480);

	// To do
	//render_image();

	//const char *mesh_filename = "mesh.off";
	//g_MeshViewer.open_mesh(mesh_filename);
	//g_MeshViewer.updateGL();

	//snapshot("snapshot");
	//

	g_MeshViewer.parse_arguments();

	osmesaFreeContext();

	return 0;
}
