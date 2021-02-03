function local_front_roughness_fit(path_folder, experiment, start_folder)

colony_mat_folder = [path_folder filesep 'whole_colony_images/colony_mat' filesep experiment];
files = dir2(colony_mat_folder);
nfiles = length(files);

%Create folder for saving MSD images and mat files if it doesn't already
%extist
if ~isdir([path_folder filesep 'front_roughness' filesep 'images' filesep experiment])
    mkdir([path_folder filesep 'front_roughness' filesep 'images' filesep experiment]);
end
if ~isdir([path_folder filesep 'front_roughness' filesep 'mat' filesep experiment])
    mkdir([path_folder filesep 'front_roughness' filesep 'mat' filesep experiment]);
end

front_roughness = nan(96,1);

%Defining the window sizes (same as for msd measurement on tracks)
%conv = 0.2481; %pixels/micron. Conversion for thresholding static particles (1 pixel)
conv = 0.8681; %for 56x: 0.8681 um/pixel
%Window sizes (in pixels)
dx_pix = 10;%5;
minL_pix = 5;
maxL_pix = 1000;

%Window sizes (in micron)
dx = dx_pix/conv;
minL = minL_pix/conv;
maxL = maxL_pix/conv;
all_L = minL:dx:maxL;

%Loop through all files
for i=start_folder:nfiles
    i
    %Figure out the well position
    filename = files(i).name;
    k = strfind(filename, '_');
    strain_name = filename(1:k-1)
    full_strain_name = filename(1:k-1);
    k2 = well2num(strain_name);
    
    %Load in the front points and fitted circle
    load([colony_mat_folder filesep filename]);
    
    %Calculate the squared displacement for each point from the circle of best fit
    [d_all, sd_all, theta_all, radii_all] = calc_sd_circ_2(x, y, a, b, R);
    
    theta_all = theta_all+pi;
    %Sort the data so that it goes around the circle
    combined = [theta_all' radii_all'];
    sorted = sortrows(combined,1);
    theta_all_sorted = sorted(:,1);
    radii_all_sorted = sorted(:,2);

    %[all_L, rho_L] = calc_running_msd(theta_all, radii_all, dx);

    %Loop through windows and calculate the mean squared displacement
    rho_L = nan(1,length(all_L));
    rho_L_ste = nan(1,length(all_L));
    for j = 1:length(all_L)
        %The window size is the arcangle, convert to angle
        L = all_L(j);
        theta = L/R;
        
        %Create variable for saving squared displacements for each
        %definition of the window
        sd_window = nan(1,length(theta_all));
        
        %Loop over definitions of window
        for k = 1:length(theta_all_sorted)
            theta_all_sorted(k)+theta;
            
            %Find the first index that is just outside this window
            ind = find(theta_all_sorted>theta_all_sorted(k)+theta, 1);
            
            %Convert to arcangles and radii
            xvals = theta_all_sorted(k:ind-1)*R;
            yvals = radii_all_sorted(k:ind-1);
            
            %figure
            %hold on
            
            %If there are at least 3 points, then calculate the squared
            %displacement vector
            if length(xvals)>2
                %[sd, fit_params] = calc_sd(xvals, yvals);
                [sd, fit_params, xrot_shift, yrot_shift] = calc_sd_continuous(xvals, yvals, 'no');
                %sd_window(k) = sum(sd)/length(sd);
                sd_window(k) = sd/L;
            end
            
            %plot(xvals, yvals)
            %plot(xvals, xvals*fit_params(1)+fit_params(2))
%             %If such a point exists, then sum the squared displacements contained within these angles,
%             %divided by the number of values summed
%             if ~isempty(ind)
%                 sd_window(k) = sum(sd_all_sorted(k:ind-1))/(ind-k);
%             end
        end
        
        %Only keep non-nan values to get the number of defined windows
        windows = sum(~isnan(sd_window));  %number of defined windows
        %If there are at least 2 defined windows, then:
        if windows>1
            %Calculate the root mean squared displacement by summing and divding by
            %the number of definitions of window
            rho_L(j) = nansum(sqrt(sd_window))/sum(~isnan(sd_window));
            rho_L_ste(j) = nanstd(sqrt(sd_window))/sqrt(sum(~isnan(sd_window)));    %standard error of the mean
        end
    end
    
    %Can improve upon above by calculating ste for each rho_L and doing weighted average
    
    %Optional plotting
    g = figure('visible', 'off');
    errorbar(all_L, rho_L, rho_L_ste)
    xlim([10^1, 10^4])
    ylim([10^-1 10^2])
    set(gca, 'XScale', 'log', 'YScale', 'log', 'FontSize', 15)
    xlabel('Window size, [\mum]')
    ylabel('Mean squared front roughness [\mum^2]')
    
    %Save figure
    saveas(g,[path_folder filesep 'front_roughness/images' filesep experiment filesep strain_name],'png')

    %Save mat file
    save([path_folder filesep 'front_roughness/mat' filesep experiment filesep 'm' strain_name], 'all_L', 'rho_L', 'rho_L_ste');
    close all

end
end
