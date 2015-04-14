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


#ifndef OPENMESHAPPS_MESHVIEWERCORET_H
#define OPENMESHAPPS_MESHVIEWERCORET_H


//== INCLUDES =================================================================

#ifdef ARCH_DARWIN
#include <glut.h>
#else
#include <GL/glut.h>
#endif

#include <string>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/IO/Options.hh>
#include <OpenMesh/Core/Utils/GenProg.hh>
#include <OpenMesh/Core/Utils/color_cast.hh>
#include <OpenMesh/Core/Mesh/Attributes.hh>
#include <OpenMesh/Tools/Utils/StripifierT.hh>
#include <OpenMesh/Tools/Utils/Timer.hh>

#include "GLViewerCore.h"


//== FORWARDS =================================================================

//class QImage;


//== CLASS DEFINITION =========================================================


template <typename M>
class MeshViewerCoreT : public GLViewerCore
{
public:

	typedef M                             Mesh;
	typedef OpenMesh::StripifierT<Mesh>   MyStripifier;
public:

	/// default constructor
	MeshViewerCoreT(GLViewerBase &_widget)
		: GLViewerCore(_widget)
		, f_strips_(false)
		, tex_id_(0)
		, tex_mode_(GL_MODULATE)
		, strips_(mesh_)
		, use_color_(true)
		, show_vnormals_(false)
		, show_fnormals_(false)
	{
	}

	/// destructor
	~MeshViewerCoreT() {}

public:

	OpenMesh::IO::Options& options() { return mesh_.options(); }
	const OpenMesh::IO::Options& options() const { return mesh_.options(); }
	void setOptions(const OpenMesh::IO::Options& opts) { mesh_.setOptions(opts); }

	/// open mesh
	virtual bool open_mesh(const char* _filename);

	/// load texture
	//virtual bool open_texture( const char *_filename );
	//bool set_texture( QImage& _texsrc );

	//void enable_strips();
	//void disable_strips();  


	Mesh& mesh() { return mesh_; }
	const Mesh& mesh() const { return mesh_; }

	virtual void draw_scene(const std::string& _draw_mode);


protected:

	/// draw the mesh
	virtual void draw_openmesh(const std::string& _drawmode);


	void glVertex(const typename Mesh::VertexHandle _vh)
	{
		glVertex3dv(&mesh_.point(_vh)[0]);
	}

	void glVertex(const typename Mesh::Point& _p)
	{
		glVertex3dv(&_p[0]);
	}

	void glNormal(const typename Mesh::VertexHandle _vh)
	{
		glNormal3dv(&mesh_.normal(_vh)[0]);
	}

	void glTexCoord(const typename Mesh::VertexHandle _vh)
	{
		glTexCoord2dv(&mesh_.texcoord(_vh)[0]);
	}

	void glColor(const typename Mesh::VertexHandle _vh)
	{
		glColor3ubv(&mesh_.color(_vh)[0]);
	}

	// face properties

	void glNormal(const typename Mesh::FaceHandle _fh)
	{
		glNormal3dv(&mesh_.normal(_fh)[0]);
	}

	void glColor(const typename Mesh::FaceHandle _fh)
	{
		glColor3ubv(&mesh_.color(_fh)[0]);
	}

	void glMaterial(const typename Mesh::VertexHandle _vh,
		int _f = GL_FRONT_AND_BACK, int _m = GL_DIFFUSE)
	{
		OpenMesh::Vec3f c = OpenMesh::color_cast<OpenMesh::Vec3f>(mesh_.color(_vh));
		OpenMesh::Vec4f m(c[0], c[1], c[2], 1.0f);

		glMaterialfv(_f, _m, &m[0]);
	}

	void glMaterial(const typename Mesh::FaceHandle _fh,
		int _f = GL_FRONT_AND_BACK, int _m = GL_DIFFUSE)
	{
		OpenMesh::Vec3f c = OpenMesh::color_cast<OpenMesh::Vec3f>(mesh_.color(_fh));
		OpenMesh::Vec4f m(c[0], c[1], c[2], 1.0f);

		glMaterialfv(_f, _m, &m[0]);
	}


protected: // Strip support

	void compute_strips(void)
	{
		if (f_strips_)
		{
			strips_.clear();
			strips_.stripify();
		}
	}

protected: // inherited

	//virtual void keyPressEvent(QKeyEvent* _event);

protected:

	bool                   f_strips_; // enable/disable strip usage
	GLuint                 tex_id_;
	GLint                  tex_mode_;

	Mesh                   mesh_;
	MyStripifier           strips_;
	bool                   use_color_;
	bool                   show_vnormals_;
	bool                   show_fnormals_;
	double                 normal_scale_;
	OpenMesh::FPropHandleT< typename Mesh::Point > fp_normal_base_;
};


//=============================================================================
#if defined(OM_INCLUDE_TEMPLATES) && !defined(OPENMESHAPPS_MESHVIEWERCORET_CPP)
#  define OPENMESH_MESHVIEWERCORE_TEMPLATES
#  include "MeshViewerCoreT.cpp"
#endif
//=============================================================================
#endif // OPENMESHAPPS_MESHVIEWERCORET_H defined
//=============================================================================

