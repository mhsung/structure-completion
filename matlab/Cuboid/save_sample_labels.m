function [] = save_sample_labels(filename, num_all_label_samples)

num_labels = length(num_all_label_samples);

file = fopen([filename, '.arff'], 'w');

fprintf(file, '@RELATION pnts-featpnts\n');
for label = 1:num_labels + 1
    fprintf(file, '@ATTRIBUTE prediction-%d NUMERIC\n', label-1);
end
fprintf(file, '@DATA\n');

for label = 1:num_labels
	str = '';
	for i = 1:num_labels+1
		if label == i
			str = [str, '1'];
		else
			str = [str, '0'];
		end
		
		if i < num_labels+1
			str = [str, ','];
		end
	end
	
	num_label_samples = num_all_label_samples(label);
	for k = 1:num_label_samples
		fprintf(file, '%s\n', str);
	end
end

fclose(file);