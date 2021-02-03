function []=main(path_folder, experiment, start_folder, masscut, maxd, radius, barrg)
% This function imports the focused images, detects the particles, does
% single particle tracking on subsequent images, then links particles in
% multiple images together. 

% INPUTS:

% Path folder: the folder containing the focused images for the experiment.

% masscut: the intensity cutoff for the detected particles

% maxd: the maximum distance (in pixels) that two particles are linked
% together in two consecutive timepoints

% radius: radius of the particle in pixels, dected by a mask of size around
%   2 or 3 px (suggest initial value = 3)

% barrg: gyration radius (linked to the radius of particule in px) (init=5)

% OUTPUTS: 

% mat file of detected particle positions and processed image in
%   particle_positions/[experiment_name]/p[strain_name].mat

% Saved image of tracks in 
%   tracks_images/[experiment_name]/t[strain_name].png, color coded by 
%   red: moving tracks, blue: static tracks

% mat file of tracks coordinates in 
%   tracks_mat/[experiment_name]/t[strain_name].png with:
%   xlinked, ylinked: x, y coordinates of moving tracks
%   xlinked_before, ylinked_before: x, y coordinates of all tracks
%   totdelx_all, totdely_all: total background drift
%   conv: conversion in microns/pixel

% Particle tracking parameters in tracking_parameters/[experiment_name].csv
% including:
%   experiment_all: name of experiment
%   strains_all: name of strain
%   radius_all
%   masscut_all
%   maxd_all
%   barrg_all
%   xmin_all (microns)
%   xmax_all (microns)
%   ymin_all (microns)
%   ymax_all (microns)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ADDITIONAL PARAMETERS

% Choose the right pixel to micron conversion for the objective
% magnification used. These are conversions for the Zeiss Axiozoom: 

conv = 0.8681; %for 56x: 0.8681 um/pixel
%conv = 4.0306; %for 16x: 0.2481 pixel/um
%conv = 0.4961; %for 32x: 0.4961 um/pixel

% Threshold number of particles as input to single particle tracking.
% Single particle tracking with more particles than this takes a long time
numpart_thresh = 800;

% Total (shortest path) displacement between beginning and end of track
% threshold. Below this distance, the track is classified as a static bead
% or noise.
delrthresh = 10;

% Can also choose to manually find the main folder of the experiment
%path_folder = uigetdir('/Volumes/PORTKeioDat/','Data changed folder'); 

% Import all subfolders (one for each colony position)
focused_images_folder = [path_folder filesep 'example_bead_images' filesep experiment];
files = dir2(focused_images_folder);
dirFlags = [files.isdir];
subFolders = files(dirFlags);

nfolders = length(subFolders);

ind = strfind(path_folder, '/');

tracking_par_filename = [path_folder filesep 'tracking_parameters' filesep experiment '.csv'];



% Some variables for storing the strains that are not computed because
% there were too many particles
not_computed = {};
index = 1;


% Loop through all subfolders
for i=start_folder:nfolders

    % Get the name of the strain without the 'f' at the end (denoting focused
    % images)
    strain_name = subFolders(i).name;
    k = strfind(strain_name,'f');
    strain_name = strain_name(1:k-1);

    disp(strcat('Detecting particles for:', strain_name))

    strain_focused_images_folder = [focused_images_folder filesep strain_name 'f/rfp'];

    % Detect the particles in the focused images
    [image_info, nfiles] = detect_particles(strain_focused_images_folder, masscut, radius);
    % image_info stores the information of all of the tracked particles
    
    % Check if the folder for saving particle positions already exists
    if ~isfolder([path_folder filesep 'particle_positions' filesep experiment])
        mkdir([path_folder filesep 'particle_positions' filesep experiment]);
    end

    % Save mat variables of particle positions
    save([path_folder filesep 'particle_positions' filesep experiment filesep 'p' strain_name '.mat'], 'image_info');
    
    % Single particle tracking
    % Calculate single step displacements and link multiple timepoints

    % Calculate the average number of particles detected per image
    all_particles = image_info(1:end,3);
    morethan_thresh = zeros(length(all_particles),1);
    all_numpart = zeros(length(all_particles),1);
    for part = 1:length(all_particles)
        particles = all_particles{part};
        numpart = size(particles,1);
        morethan_thresh(part) = (numpart>numpart_thresh);
        all_numpart(part) = numpart;
    end

    disp(strcat('Number of particles found:', num2str(mean(all_numpart)), '+/-', num2str(std(all_numpart))));

    % If the number of particles is more than the threshold, then save the
    % strain name and continue onto the next strain.

    if sum(morethan_thresh)>0
        not_computed{index} = strain_name;
        index = index + 1;
        continue
    end

    disp(strcat('Calculating single step displacements for:', strain_name))

    % Parameters for single particle tracking

    % Set the region used for correcting for background drift with static beads

    xmin = 0;
    % For the first row of strains, the static region is below the colony,
    % whereas for the subsequent rows, the static region is above the colony
    if str2num(strain_name)<12
        ymin = 900;%900;
        ymax = 1040;
    else
        ymin = 0;
        ymax = 100;%100;
    end
    xmax = 1388;%1388;

    % Threshold distance for calling a static particle
    rthresh = 0.5;%0.5;
    ythresh_min = -inf;
    ythresh_max = 0;

    % Initialize some parameters for saving the background drift information
    del_tot = cell(nfiles, 3); %by columns: delx, dely, delr
    xfinal_tot = cell(nfiles, 1);
    yfinal_tot = cell(nfiles, 1);

    prevdelx = 0;
    prevdely = 0;

    prevdelx_all = nan(1,nfiles-1);
    prevdely_all = nan(1,nfiles-1);

    totdelx_all = nan(1,nfiles-1);
    totdely_all = nan(1,nfiles-1);

    % Loop through all images
    for t = 1:nfiles-2%nfiles-1
        t
        % Run single particle tracking code for two consecutive timesteps
        [xplot, yplot] = spt_one_timestep_2(image_info, t, maxd, barrg);

        % Correct for background drift
        [xcorr, ycorr, totdelx, totdely] = back_drift_t1(xplot, yplot, xmin, xmax,...
            ymin, ymax, prevdelx, prevdely, 0);

        % Save the background drift
        prevdelx_all(t) = prevdelx;
        prevdely_all(t) = prevdely;
        totdelx_all(t) = totdelx;
        totdely_all(t) = totdely;

        % Save the total drifted value
        prevdelx = totdelx;
        prevdely = totdely; 

        % These are the corrected tracks
        xfinal = xcorr;
        yfinal = ycorr;

        %[xfinal, yfinal, delr] = stat_part_t1(xcorr, ycorr, rthresh, ythresh_min, ythresh_max, 0);
        %length = sum(~isnan(xfinal(1, :)))
        %figure
        % delx = diff(xfinal, 1);
        % dely = diff(yfinal, 1);
        % delr = sqrt(delx.^2 + dely.^2);
        % del_tot{t, 1} = delx'*conv;
        % del_tot{t, 2} = dely'*conv;
        % del_tot{t, 3} = delr'*conv;

        % Convert from pixels to microns
        xfinal_tot{t, 1} = xfinal'*conv;
        yfinal_tot{t, 1} = yfinal'*conv;
        % 
         %plot(xplot, yplot, 'g', 'LineWidth', 2)
         %plot(xcorr*conv, ycorr*conv, 'b', 'LineWidth', 1)
         %plot(xfinal*conv, yfinal*conv, 'r', 'LineWidth', 1)

    end

    % Link multiple time points
    disp(strcat('Linking multiple timepoints for: ', strain_name));

    [xlinked, ylinked] = link_multt(xfinal_tot, yfinal_tot);

    xlinked_before = xlinked;
    ylinked_before = ylinked;

    todelete = [];

    % Loop through each track
    for tracknum = 1:size(xlinked_before,2)
        xlinked_current = xlinked_before(:,tracknum);
        ylinked_current = ylinked_before(:,tracknum);
        xlinked_current_before = xlinked_current;

        % Get rid of NaN values
        xlinked_current = xlinked_current(~isnan(xlinked_current_before));
        ylinked_current = ylinked_current(~isnan(xlinked_current_before));

        % Calculate distance between first and last point
        delx = xlinked_current(end) - xlinked_current(1);
        dely = ylinked_current(end) - ylinked_current(1);
        delr = sqrt(delx^2 + dely^2);
        if delr<delrthresh
            %Set those tracks to NaN
            todelete = [todelete tracknum];
        end
    end

    % Set tracks that too short to NaN
    xlinked(:,todelete) = [];
    ylinked(:,todelete) = [];

    % Plot the tracks. Red are all kept tracks, blue are excluded tracks.
    f = figure('visible','off');
    hold on
    plot(xlinked_before, ylinked_before, 'b')
    plot(xlinked,ylinked, 'r')
    xlim([0 1388*conv]);
    ylim([0 1040*conv]);
    legend(strain_name);
    hold off

    % Save track files
    disp(strcat('Saving track files for: ', strain_name));

    % Check if the folder for saving tracks already exists
    if ~isfolder([path_folder filesep 'tracks_images' filesep experiment])
        mkdir([path_folder filesep 'tracks_images' filesep experiment]);
    end
    if ~isfolder([path_folder filesep 'tracks_mat' filesep experiment])
        mkdir([path_folder filesep 'tracks_mat' filesep experiment]);
    end

    % Save figure of tracks
    saveas(f,[path_folder filesep 'tracks_images' filesep experiment filesep 't' strain_name],'png')

    % Save mat variables of tracks
    save([path_folder filesep 'tracks_mat' filesep experiment filesep 't' strain_name '.mat'], 'xlinked', 'ylinked', 'xlinked_before', 'ylinked_before', 'totdelx_all', 'totdely_all', 'conv');
    
    % Close the figure
    close 
    
    % Check if there is already a file where the experiment parameters are
    % saved
    if isfile(tracking_par_filename)
        T = readtable(tracking_par_filename, 'ReadVariableNames', true);
        %T.experiment_all(1) = {'test'}
    else
        experiment_all = cell(nfolders,1);
        strains_all = nan(nfolders,1);
        radius_all = nan(nfolders,1);
        masscut_all = nan(nfolders,1);
        barrg_all = nan(nfolders,1);
        maxd_all = nan(nfolders,1);
        xmin_all = nan(nfolders,1);
        xmax_all = nan(nfolders,1);
        ymin_all = nan(nfolders,1);
        ymax_all = nan(nfolders,1);
        manually_changed = zeros(nfolders, 1);
        notes = cell(nfolders, 1);
        T = table(experiment_all, strains_all, radius_all, masscut_all, maxd_all, barrg_all, xmin_all, xmax_all, ymin_all, ymax_all, manually_changed, notes);
    end
    
    % Save the tracking parameters 
    ind = str2num(strain_name)+1;

    T(ind, 'experiment_all') = {experiment};
    T(ind, 'strains_all') = {str2num(strain_name)};
    T(ind, 'radius_all') = {radius};
    T(ind, 'masscut_all') = {masscut};
    T(ind, 'barrg_all') = {barrg};
    T(ind, 'maxd_all') = {maxd};
    T(ind, 'xmin_all') = {round(xmin*conv)};
    T(ind, 'xmax_all') = {round(xmax*conv)};
    T(ind, 'ymin_all') = {round(ymin*conv)};
    T(ind, 'ymax_all') = {round(ymax*conv)};

    % Save the tracking parameters in a csv file

    writetable(T, tracking_par_filename);
end


end