function [] = save_samples_pts(filename, samples)

assert(size(samples, 1) == 7);

file = fopen([filename, '.pts'], 'w');

for i = 1:size(samples, 2)
    fprintf(file, '%d %f %f %f %f %f %f\n',...
        samples(1, i), samples(2, i), samples(3, i), samples(4, i),...
        samples(5, i), samples(6, i), samples(7, i));
end

fclose(file);