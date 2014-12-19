function [] = save_cuboid_obj(filename, cuboid)


assert(size(cuboid, 1) == 3);

file = fopen([filename, '.obj'], 'w');

for i = 1:size(cuboid, 2)
    fprintf(file, 'v %f %f %f\n',...
        cuboid(1, i), cuboid(2, i), cuboid(3, i));
end

for i = 1:(size(cuboid, 2) / 3)
    fprintf(file, 'f %d %d %d\n',...
        3*(i-1)+1, 3*(i-1)+2, 3*(i-1)+3);
end

fclose(file);