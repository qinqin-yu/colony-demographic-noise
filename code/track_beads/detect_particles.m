function [image_info, nfiles] = detect_particles(path_folder, masscut, radius)
%%% Parameters
%%%% initial feature detections
featsize=radius;%4         % radius of the particule (px), dected by a mask of size around 2 or 3 px (init=3)
%masscut=100000; 	      % cutting intensity for particule detection (init=100000)
%%%% minimal distance beetwen beads (first image)
%tolproxy=1;
%%%% other rejection processes
%barcc=1.0;          % excentricity (0=circular, 1=linear), usually around 0.1 or 0.2 (init=0.1)
%barrg=5;            % gyration radius (linked to the radius of particule in px) (init=5)
%IdivRg=0;      % minimal integrated intensity required to detect a particule (init=0)
%%%% cutoff on max displacement
%maxd=6;%4;%12;        % max displacement of a bead between 2 frames inside the mask (px)

%barcc_s=0.2;

%Nb_frames=1;

%Import files

%path=pwd;
%Ia=imread(fullfile(path,'0010.tif'));
%Ib=imread(fullfile(path,'0011.tif'));
%Ic=imread(fullfile(path,'0012.tif'));

%path_folder = ['/Users/qinqinyu/Documents/hallatschek_lab/data_changed/' ...
%    '20190129_keio_test/3595_focused/rfp1-4/'];

%path_folder = ['../../../data_changed/' ...
%    '20181118_ecoli_beads_high_mag/'];

oldpath = path;
%newpath = path_folder;
addpath(path_folder);

%Import cropped images (where best focus will be calculated) into a cell
imagefiles = dir2(strcat(path_folder, '/*.tif'));      
nfiles = length(imagefiles)    % Number of files found

image_info = cell(nfiles, 3);   
%image matix | t | z | c

%Add information about the images from the file names
%figure
colors = {'r', 'y', 'g', 'c', 'b', 'm'};
for ii=1:nfiles
%     if mod(ii, 10) == 0
%         ii
%     end
   ii
   currentfilename = imagefiles(ii).name;
   %[path_folder filesep currentfilename]
   currentimage = imread([path_folder filesep currentfilename]);

   image_info{ii, 2} = str2num(currentfilename(1:end-4));   %t
   %image_info{ii, 3} = str2num(currentfilename(8:10));  %z
   %image_info{ii, 4} = str2num(currentfilename(13:15)); %c
   %image_info{ii, 5} = std(im2double(currentimage(:)));
   
   %Increase contrast

   background = imopen(currentimage,strel('disk',12));%12));
   
   currentimage2=imsubtract(currentimage,background);
   
   contrastedimage=double(contrast(currentimage2,0.1,0.1));
   %contrastedimage=double(contrast(currentimage,0.1,0.1));
   
   image_info{ii, 1} = contrastedimage;    %image matrix
   
   M_PA = feature2D(contrastedimage,1,featsize,masscut);
   %M_PA = feature2D(currentimage,1,featsize,masscut);
   image_info{ii, 3} = M_PA;
%    if mod(ii,6) == 1
%         figure
%         imshow(contrastedimage,[])
%         hold on
%         plot(M_PA(:,1),M_PA(:,2),['+' colors{mod(ii,6)}],'MarkerSize',8)
%    elseif mod(ii,6) == 0
%         plot(M_PA(:,1),M_PA(:,2),['+' colors{6}],'MarkerSize',8)
%    else
%        %subplot(1,nfiles,ii), hold on, imshow(contrastedimage,[]), plot(M_PA(:,1),M_PA(:,2),'+r','MarkerSize',8)
%         plot(M_PA(:,1),M_PA(:,2),['+' colors{mod(ii, 6)}],'MarkerSize',8)
%    end

%     subplot(2, 2, 1);
%     imshow(currentimage, [63 121]);
%     subplot(2, 2, 2);
%     imshow(background, []);
%     subplot(2, 2, 3);
%     imshow(currentimage2, [63 121]);
%     %subplot(2, 2, 4);
%     f = figure;
%     ax = axes('Parent', f);
%     imshow(contrastedimage, []);
%     hold on
%     plot(ax,M_PA(:,1),M_PA(:,2),'.r','MarkerSize',10)
%     hold off

end

path(oldpath)
end