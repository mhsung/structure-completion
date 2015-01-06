close all;
clear all;

scale = 0.6;
input_directory = 'C:/project/app/cuboid-prediction/experiments/experiment_two_types_1/output/';
output_directory = 'C:/Users/Administrator/Dropbox/Public/web/two_types_chairs/';
listing = dir(input_directory);

for i = 1:length(listing)
    filename = listing(i).name;
    [pathstr,name,ext] = fileparts(filename);

    if strcmp(ext, '.png') == 1 && strcmp(name(end-1:end), '_0') == 1
        name_prefix = name(1:end-2);
        
        file0 = [input_directory, name_prefix, '_0.png'];
        file1 = [input_directory, name_prefix, '_1.png'];
        file2 = [input_directory, name_prefix, '_2.png'];
        %file3 = [input_directory, name_prefix, '_3.png'];
        file4 = [input_directory, name_prefix, '_4.png'];
        
        if exist(file0,'file') ~= 2; continue; end;
        if exist(file1,'file') ~= 2; continue; end;
        if exist(file2,'file') ~= 2; continue; end;
        %if exist(file3,'file') ~= 2; continue; end;
        if exist(file4,'file') ~= 2; continue; end;
        
        image0 = imread(file0);
        image1 = imread(file1);
        image2 = imread(file2);
        %image3 = imread(file3);
        image4 = imread(file4);
        
        H = vision.TextInserter('(a) Input');
        H.Font = 'Calibri';
        H.FontSize = 36;
        H.Location = [size(image0, 2)*0.2, size(image0, 1)*0.9];
        image0 = step(H, image0);
        
        H = vision.TextInserter('(a) Point Clustring');
        H.Font = 'Calibri';
        H.FontSize = 36;
        H.Location = [size(image1, 2)*0.2, size(image1, 1)*0.9];
        image1 = step(H, image1);
        
        H = vision.TextInserter('(b) Label & Axis Recognition');
        H.Font = 'Calibri';
        H.FontSize = 36;
        H.Location = [size(image2, 2)*0.2, size(image2, 1)*0.9];
        image2 = step(H, image2);
        
%         H = vision.TextInserter('(c) Corner Optimization');
%         H.Font = 'Calibri';
%         H.FontSize = 36;
%         H.Location = [size(image3, 2)*0.2, size(image3, 1)*0.9];
%         image3 = step(H, image3);
        
        H = vision.TextInserter('(d) Corner Optimization');
        H.Font = 'Calibri';
        H.FontSize = 36;
        H.Location = [size(image4, 2)*0.2, size(image4, 1)*0.9];
        image4 = step(H, image4);
        
        %image = [image1, image2; image3, image4];
        image = [image0, image1; image2, image4];
        image = imresize(image, scale);
        
        file = [output_directory, name_prefix, '.png']
        imwrite(image, file);
    end
end


