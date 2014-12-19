#version 150

uniform mat4 modelview_matrix;
uniform vec3 sample_points[256];
uniform int num_sample_points;
uniform float radius;
uniform vec3 bbox_center;
uniform vec3 bbox_size;

in vec3 frag_position;
out vec4 frag_color;

void main(void)
{
	frag_color.rgb = frag_position;

	// Alpha = 1 means that the voxel is perfectly visible.
	frag_color.a = 1.0;

	// 'frag_position' has value [0, 1], which should be mapped to
	// [-0.55 * bbox_size_ + bbox_center_, +0.55 * bbox_size_ + bbox_center_].
	vec3 voxel_point = (frag_position - 0.5) * bbox_size * 1.1 + bbox_center;

	// NOTE:
	// 'radius' should be greater than 0.

	// NOTE:
	// 'num_sample_points' should be less max_vertex_uniforms_vectors / sizeof(GLfloat).

	for(int i = 0; i < num_sample_points; ++i)
	{
		vec3 sample_point = sample_points[i];

		// Positions in the model view coordinates.
		vec4 lc_sample_point_4 = modelview_matrix * vec4(sample_point, 1.0);
		vec3 lc_sample_point = lc_sample_point_4.xyz / lc_sample_point_4.w;
		lc_sample_point.z -= radius;

		vec4 lc_voxel_point_4 = modelview_matrix * vec4(voxel_point, 1.0);
		vec3 lc_voxel_point = lc_voxel_point_4.xyz / lc_voxel_point_4.w;

		float lc_voxel_point_len = length(lc_voxel_point);
		float lc_sample_point_len = length(lc_sample_point);

		float dot_prod = dot(lc_voxel_point, lc_sample_point);

		// Ignore a sample point if either it or the voxel point is not visible
		// from this model view. 
		if (lc_sample_point.z >= 0 || lc_voxel_point.z >= 0
			|| lc_voxel_point_len == 0 || lc_sample_point_len == 0)
			continue;

		float cos_angle = dot_prod / (lc_voxel_point_len * lc_sample_point_len);
		float sin_angle = sqrt(1 - cos_angle * cos_angle);
		float tan_angle = sin_angle / cos_angle;
		
		if (dot_prod < lc_sample_point_len * lc_sample_point_len)
		{
			// NOTE:
			// If the following code is uncommented, a point is considered as a sphere.
			//float distance = length(voxel_point - sample_point) / radius;
			//frag_color.a = min(frag_color.a, distance);
		}
		else
		{
			float distance = (lc_sample_point_len * tan_angle) / radius;
			//frag_color.a = min(frag_color.a, distance);
			// NOTE:
			//
			if (distance < 1) frag_color.a = 0;
		}
	}
}