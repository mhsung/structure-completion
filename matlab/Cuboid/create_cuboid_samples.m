function[samples] = create_cuboid_samples(cuboid, start_face_index)

num_triangles = (size(cuboid, 2) / 3);

samples = [];

for face_index = 0:num_triangles-1
    v = cuboid(:, (3*face_index+1):(3*face_index+3));
    
    area = 0.5*norm(cross(v(:, 3)-v(:, 1), v(:, 3)-v(:, 2)));
    num_triangle_samples = max(round(1000 * area), 3);
    
    w = rand(3, num_triangle_samples);
    w = w ./ repmat(sum(w, 1), 3, 1);
    
    assert(num_triangle_samples >= 3);
    w(:, 1) = [1; 0; 0];
    w(:, 2) = [0; 1; 0];
    w(:, 3) = [0; 0; 1];
    
    p = v * w;
    
    f = (start_face_index + face_index) * ones(1, num_triangle_samples);
    
    samples = [samples, [f; w; p]];
end

