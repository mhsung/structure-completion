//////////////////////////////////////////////////////////////////////////////
//
//  --- Angel.h ---
//
//   The main header file for all examples from Angel 6th Edition
//
//////////////////////////////////////////////////////////////////////////////

#ifndef __ANGEL_H__
#define __ANGEL_H__

//----------------------------------------------------------------------------
// 
// --- Include system headers ---
//
#include <iostream>
#include <stdio.h>

//----------------------------------------------------------------------------
//
// --- Include OpenGL header files and helpers ---
//
//   The location of these files vary by operating system.  We've included
//     copies of open-source project headers in the "GL" directory local
//     this this "include" directory.
//

#ifdef __APPLE__  // include Mac OS X versions of headers
#  include <OpenGL/OpenGL.h>
#  include <GLUT/glut.h>
#else // non-Mac OS X operating systems
#  include <GL/glew.h>
#  include <GL/freeglut.h>
#  include <GL/freeglut_ext.h>
#endif  // __APPLE__

// Define a helpful macro for handling offsets into buffer objects
#define BUFFER_OFFSET( offset )   ((GLvoid*) (offset))

//----------------------------------------------------------------------------
//
//  --- Include our class libraries and constants ---
//

namespace Angel {

// Helper function to load a text file
char* readShaderSource(const char* shaderFile);

//  Helper function to load vertex and fragment shader files
GLuint InitShader( const char* vertexShaderFile,
		   const char* fragmentShaderFile ,
		   const char* outputAttributeName);
}  // namespace Angel


#define Print(x)  do { std::cerr << #x " = " << (x) << std::endl; } while(0)

//  Globally use our namespace in our example programs.
using namespace Angel;

#endif // __ANGEL_H__
