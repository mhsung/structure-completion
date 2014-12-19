function [] = save_cuboid_labels(filename, num_all_label_faces)

num_labels = length(num_all_label_faces);
accumulated_num_label_faces = 0;

file = fopen([filename, '.comp'], 'w');

for label = 1:num_labels
    fprintf(file, '%d\n', label - 1);
	
	num_label_samples = num_all_label_faces(label);
    
	for i = 1:num_label_samples
        accumulated_num_label_faces = accumulated_num_label_faces + 1;
        fprintf(file, '%d ', accumulated_num_label_faces);
    end
    
    fprintf(file, '\n');
end

fclose(file);