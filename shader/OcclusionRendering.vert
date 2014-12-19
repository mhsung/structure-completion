#version 150

in vec3 position;
in vec2 texcoord;
out vec3 frag_position;

void main (void)
{
    frag_position = position;

	gl_Position.xy = texcoord.xy;
	gl_Position.z = 0;
}