function get_colony_mask(path_folder)
%This code takes the bw images from the calibration (taken at the beginning
%of imaging) and extracts x, y coordinates, fits circle, saves these
%variables and an overlay image

conv = 0.1240; 
%0.1240 pixels = 1 um for 8x; 
%0.0775 pixels = 1 um for 5x;
%0.8681 pixels = 1 um for 56x;

%path_folder = '/Volumes/PORTKeioDat/20190910_2036_KQY4a';
colony_images_binary_folder = [path_folder filesep 'colony_images_binary'];
colony_images_folder = [path_folder filesep 'colony_images2'];
colony_images_edge_folder = [path_folder filesep 'colony_edge'];

if ~isdir([path_folder filesep 'colony_edge'])
    mkdir([path_folder filesep 'colony_edge']);
end
if ~isdir([path_folder filesep 'colony_mat'])
    mkdir([path_folder filesep 'colony_mat']);
end

files = dir2(colony_images_binary_folder);
nfiles = length(files);

se = strel('disk', 10);
for i = 1:nfiles
    image_name = files(i).name;
    k = strfind(image_name,'.tif');
    image_name_no_tif = image_name(1:k-1);
        
    I_filename = [colony_images_folder filesep image_name];
    Ibw_filename = [colony_images_binary_folder filesep image_name];
    
    I = imread(I_filename);
    Ibw = imread(Ibw_filename);
    Ifill = imfill(Ibw);
    Idil = imdilate(Ifill, se);
    Iperim = bwperim(Idil);
    Ioverlay = I.*uint16(imcomplement(Iperim));
    [y, x] = find(Iperim);
    
    if length(x)>1
            [a, b, R] = find_circles(x, y);
            output_image_filename = [path_folder filesep 'colony_edge' filesep 'bw' image_name];
            %plot_image_fitted_circle(output_image_filename, I, x, y, a, b, R, 'off');
            %saveas(f,[path_folder filesep 'colony_edge' filesep image_files(j).name 'bw'],'png');
            x = x/conv;
            y = y/conv;
            a = a/conv;
            b = b/conv;
            R = R/conv;
            save([path_folder filesep 'colony_mat' filesep image_name_no_tif '.mat'], 'x', 'y', 'a', 'b', 'R', 'conv');
            imwrite(Ioverlay, output_image_filename);
    end
        
end
