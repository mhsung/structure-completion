#include "Angel.h"

namespace Angel {

	// Create a NULL-terminated string by reading the provided file
	char* readShaderSource(const char* shaderFile)
	{
		FILE *filePointer;
		char *content = NULL;

		int count=0;

		if (shaderFile != NULL) {
			filePointer = fopen(shaderFile,"rt");

			if (filePointer != NULL) {

				fseek(filePointer, 0, SEEK_END);
				count = ftell(filePointer);
				rewind(filePointer);

				if (count > 0) {
					content = (char *)malloc(sizeof(char) * (count+1));
					count = fread(content,sizeof(char),count,filePointer);
					content[count] = '\0';
				}
				fclose(filePointer);
			}
		}
		return content;

	}


	// Create a GLSL program object from vertex and fragment shader files
	GLuint InitShader(const char* vShaderFile, const char* fShaderFile, const char* outputAttributeName)
	{
		struct Shader {
			const char*  filename;
			GLenum       type;
			GLchar*      source;
		}  shaders[2] = {
			{ vShaderFile, GL_VERTEX_SHADER, NULL },
			{ fShaderFile, GL_FRAGMENT_SHADER, NULL }
		};

		GLuint program = glCreateProgram();

		for ( int i = 0; i < 2; ++i ) {
			Shader& s = shaders[i];
			s.source = readShaderSource( s.filename );
			if ( shaders[i].source == NULL ) {
				std::cerr << "Failed to read " << s.filename << std::endl;
				exit( EXIT_FAILURE );
			}
			GLuint shader = glCreateShader( s.type );
			glShaderSource( shader, 1, (const GLchar**) &s.source, NULL );
			glCompileShader( shader );

			GLint  compiled;
			glGetShaderiv( shader, GL_COMPILE_STATUS, &compiled );
			if ( !compiled ) {
				std::cerr << s.filename << " failed to compile:" << std::endl;
				GLint  logSize;
				glGetShaderiv( shader, GL_INFO_LOG_LENGTH, &logSize );
				char* logMsg = new char[logSize];
				glGetShaderInfoLog( shader, logSize, NULL, logMsg );
				std::cerr << logMsg << std::endl;
				delete [] logMsg;

				exit( EXIT_FAILURE );
			}

			delete [] s.source;

			glAttachShader( program, shader );
		}

		/* Link output */
		glBindFragDataLocation(program, 0, outputAttributeName);

		/* link  and error check */
		glLinkProgram(program);

		GLint  linked;
		glGetProgramiv( program, GL_LINK_STATUS, &linked );
		if ( !linked ) {
			std::cerr << "Shader program failed to link" << std::endl;
			GLint  logSize;
			glGetProgramiv( program, GL_INFO_LOG_LENGTH, &logSize);
			char* logMsg = new char[logSize];
			glGetProgramInfoLog( program, logSize, NULL, logMsg );
			std::cerr << logMsg << std::endl;
			delete [] logMsg;

			exit( EXIT_FAILURE );
		}
		
		/* use program object */
		glUseProgram(program);

		return program;
	}

}  // Close namespace Angel block
