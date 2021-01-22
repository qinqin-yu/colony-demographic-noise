function [] = msd_continuous(path_name, folder)
%%%%CHANGE BEFORE STARTING%%%%
%Conditions on tracks
minlength = 30;%50; %in um
maxlength = inf;

%conv = 0.2481; %pixels/micron. Conversion for thresholding static particles (1 pixel)
conv = 0.8681; %for 56x: 0.8681 pixels/micron
%Window sizes (in pixels)
dx_pix = 5;
minL_pix = 5;
maxL_pix = 1000;

%Window sizes (in micron)
dx = dx_pix/conv;
minL = minL_pix/conv;
maxL = maxL_pix/conv;
all_L = minL:dx:maxL;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Get main directory for data
%path_folder = uigetdir('/Users/qinqinyu/Documents/hallatschek_lab/data_changed','Data changed folder'); 

%Get all files from the tracks mat folder
[path_name filesep 'tracks_mat' filesep folder]
files = dir2([path_name filesep 'tracks_mat' filesep folder]);

%Create folder for saving MSD images and mat files if it doesn't already
%extist
if ~isfolder([path_name filesep 'msd_images' filesep folder])
    mkdir([path_name filesep 'msd_images' filesep folder]);
end
if ~isfolder([path_name filesep 'msd_mat' filesep folder])
    mkdir([path_name filesep 'msd_mat' filesep folder]);
end

nfiles = length(files);

%Loop through all files
for m = 1:nfiles
    %Get strain name
    strain_name = files(m).name;
    k = strfind(strain_name,'.');
    strain_name = strain_name(2:k-1);
    if ~isempty(strain_name)
    
    %Display message to indicate which file is currently being processed
    disp(strcat('Calculating MSD for:', strain_name))
    
    %Import mat file
    [path_name filesep 'tracks_mat' filesep folder]
    filename = [path_name filesep 'tracks_mat' filesep folder filesep 't' strain_name '.mat'];
    load(filename);

    %Copy preconditioned track coordinates into new variable
    xlinked_ori = xlinked;
    ylinked_ori = ylinked;
    
    %Condition the track length based on input above
    [xlinked, ylinked] = condition_track_length(xlinked_ori, ylinked_ori, minlength, maxlength);
    
    vecx = diff(xlinked);
    vecy = diff(ylinked);
    vecr = sqrt(vecx.^2 + vecy.^2);

    xlinked_min_step = xlinked(1:end-1,:);
    ylinked_min_step = ylinked(1:end-1,:);
    xlinked_min_step(vecr<1/conv) = NaN;
    ylinked_min_step(vecr<1/conv) = NaN;
    
    %If there are still points in the conditioned tracks, then calculate
    %the MSD
    if ~isempty(xlinked)
        %Weight average by number of windows that fit into a track, minimum
        %of 5 tracks to do average
        %[msd_avg_weighted, msd_ste_weighted, ensemble_msd, ensemble_ste] = msd_continuous_avg_tracks_weighted(xlinked,ylinked,all_L,5,'N');
        
        [msd_avg_weighted, msd_ste_weighted, ensemble_msd, ensemble_ste] = msd_continuous_avg_tracks_weighted(xlinked_min_step,ylinked_min_step,all_L,5,'N');

        %Plot figure
        g = figure('visible', 'off');
        hold on
        
        plot(all_L, ensemble_msd, 'Color', [0.9 0.9 0.9]); 
        errorbar(all_L, msd_avg_weighted, msd_ste_weighted, 'ro-')
        %errorbar(all_L, msd_avg_weighted_ms, msd_ste_weighted_ms, 'bo-')


        xlim([5 3*10^3])
        ylim([10^-4 10^3])
        title(strain_name);
        xlabel('L [\mu m]');
        ylabel('MSD [\mu m^2]');
        set(gca, 'XScale', 'log', 'YScale', 'log', 'FontSize', 15)

        %Save figure
        saveas(g,[path_name filesep 'msd_images' filesep folder filesep 'mb' strain_name],'png')
        
        %Save mat file
        save([path_name filesep 'msd_mat' filesep folder filesep '/m' strain_name], 'all_L', 'msd_avg_weighted', 'msd_ste_weighted', 'ensemble_msd', 'ensemble_ste');
        close all
    end
    end
end
end