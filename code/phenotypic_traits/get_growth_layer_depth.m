function get_growth_layer_depth(path_folder)

    % Define parameters
    vye_min = 1; % To prevent beads that are randomly moving in back from contributing and that are stationary
    dxet_max = 50; % Neighborhood in x to look for

    % Conditioning total track length
    minlength = 0;  
    maxlength = inf;

    %Create folder for saving growth layer depth images and mat files if it doesn't already
    %extist
    if ~isdir([path_folder filesep 'growth_layer_depth_images'])
        mkdir([path_folder filesep 'growth_layer_depth_images']);
    end
    if ~isdir([path_folder filesep 'growth_layer_depth_mat'])
        mkdir([path_folder filesep 'growth_layer_depth_mat']);
    end

    d_all = nan(96, 1);
    derr_all = nan(96,1);

    files = dir2([path_folder filesep 'tracks_mat']);

    function_exp        = @(x,a) a(1)*exp(a(2)*x);
    a0 = [1,-1/100];

    for i = 1:length(files)

        filename = files(i).name;
        k = strfind(filename, '.');
        k2 = strfind(filename, 't');
        strain_name = filename(k2+1:k-1)

        load([path_folder filesep 'tracks_mat' filesep files(i).name]);

        if ~isempty(xlinked)

            [x, y, total_lengths] = condition_track_length(xlinked, ylinked, minlength, maxlength);

            % Calculate the displacements
            vx = diff(x);
            vy = diff(y);
            vr = sqrt(vx.^2 + vy.^2);

            vx(size(vx,1)+1,:) = nan(1,size(vx,2));
            vy(size(vy,1)+1,:) = nan(1,size(vy,2));
            vr(size(vr,1)+1,:) = nan(1,size(vr,2));

            % Calculate the overall direction of motion (unit vector e)
            Ex = nansum(nansum(vx));
            Ey = nansum(nansum(vy));
            Er = sqrt(Ex.^2+Ey.^2);

            ex = Ex/Er;
            ey = Ey/Er;
            e = [ex,ey];

            % Rotate coordinates
            theta = atan2(ex,ey); % to rotate 90 counterclockwise
            R = [cos(theta) -sin(theta); sin(theta) cos(theta)];

            xe = x*cos(theta) - y*sin(theta);
            ye = x*sin(theta) + y*cos(theta);

            vxe = vx*cos(theta) - vy*sin(theta);
            vye = vx*sin(theta) + vy*cos(theta);

            % Remove steps that are too small
            xe(vye<vye_min) = NaN;
            ye(vye<vye_min) = NaN;
            vxe(vye<vye_min) = NaN;
            vye(vye<vye_min) = NaN;

            n = size(x,2);
            nt = size(xe,1);

            dyet_final_all = nan(nt,n^2);
            ratio_final_all = nan(nt,n^2);

        %     all_struct = struct('dy',{},'ratio',{});

            % Loop through all timepoints
            for t = 1:nt

                % Get the points in that time point
                vxet = vxe(t,:);
                vyet = vye(t,:);
                xet = xe(t,:);
                yet = ye(t,:);

                % Calculate every pair's distance and velocity ratio in the y (rot) direction
                [xet1,xet2] = meshgrid(xet,xet);
                [yet1,yet2] = meshgrid(yet,yet);
                [vyet1,vyet2] = meshgrid(vyet,vyet);

                dxet = xet2-xet1;
                dyet = yet2-yet1;

                ratio = vyet2./vyet1;

                % Filter so that we're only looking at displacements in x that are less
                % than 200 away
                dyet_filt = dyet(abs(dxet)<dxet_max);
                ratio_filt = ratio(abs(dxet)<dxet_max);

                % Save into a new variable
                dyet_final_all(t,1:length(dyet_filt(:))) = -dyet_filt(:)';%-dyet(:)';
                ratio_final_all(t,1:length(ratio_filt(:))) = ratio_filt(:);%ratio(:)';

        %         % And save into the structure to preserve time information
        %         all_struct(t).dy = -dyet_filt(:)';
        %         all_struct(t).ratio = ratio_filt(:);
            end

        %     Bin ratio data
            if max(dyet_final_all(:))-min(dyet_final_all(:))<100
                bin_width = nanstd(dyet_final_all(:))*2/4;%10;
            else
                bin_width = 10;
            end
            %bin_width = 10;
            nbins = max(max(dyet_final_all));

            bins = -nbins:bin_width:nbins;
            bins_center = bins(1:end-1) + bin_width/2;
            mean = nan(1,length(bins)-1);
            std = nan(1,length(bins)-1);
            ste = nan(1,length(bins)-1);
            for j =2:length(bins)
                inds = dyet_final_all<bins(j) & dyet_final_all>bins(j-1);
                data = ratio_final_all(inds);
                nonan_len = length(data(~isnan(data)));
                if nonan_len>2
                    mean(j-1) = nanmean(data);
                    std(j-1) = nanstd(data);
                    ste(j-1) = std(j-1)/sqrt(nonan_len);
                end
            end

            bins_center_nonan = bins_center(~isnan(mean) & ste>0);
            mean_nonan = mean(~isnan(mean) & ste>0);
            ste_nonan = ste(~isnan(mean) & ste>0);

            if length(mean_nonan)>2


                [a,aerr,chisq,yfit] = gradsearch(bins_center_nonan(abs(bins_center_nonan)<dxet_max),mean_nonan(abs(bins_center_nonan)<dxet_max),ste_nonan(abs(bins_center_nonan)<dxet_max),function_exp,a0);
                d = -1/a(2);
                derr = d*abs(aerr(2)/a(2));

                d_all(str2num(strain_name)+1) = d;
                derr_all(str2num(strain_name)+1) = derr;

                h = zeros(1,5);

                g = figure('visible', 'off');
                hold on
                h(1) = scatter(dyet_final_all(:), ratio_final_all(:), 'MarkerEdgeColor', [0.8,0.8,0.8], 'DisplayName', 'raw data');
                h(2) = errorbar(bins_center_nonan, mean_nonan, ste_nonan, 'ko-','DisplayName', 'binned data');
                h(3) = plot(bins_center_nonan(abs(bins_center_nonan)<dxet_max),yfit, 'r-','DisplayName',['exponential fit, d = ' num2str(round(d,2)) '+/-' num2str(round(derr,2)) '\mum']);
                h(4) = plot([dxet_max,dxet_max],[min(ratio_final_all(:)), max(ratio_final_all(:))], 'k--');
                h(5) = plot([-dxet_max,-dxet_max],[min(ratio_final_all(:)), max(ratio_final_all(:))], 'k--');
                xlabel('-\Delta y (\mum)')
                ylabel('v_y(y_2)/v_y(y_1)')
                legend(h(1:3))
                xlim([-300,300])
                ylim([5*10^(-2), 0.5*10^2])
                set(gca,'YScale','log')

                %Save figure
                saveas(g,[path_folder filesep 'growth_layer_depth_images' filesep strain_name],'pdf')

                dy_binned = bins_center_nonan;
                ratio_binned = mean_nonan;
                ratio_err_binned = ste_nonan;
                dy = dyet_final_all;
                ratio = ratio_final_all;

                %Save mat file
                save([path_folder filesep 'growth_layer_depth_mat' filesep strain_name], 'dy_binned', 'ratio_binned', 'ratio_err_binned', 'dy', 'ratio', 'd', 'derr');
                close all
            end
        end
    end

    csvwrite([path_folder filesep 'summary_data/growth_layer_depth.csv'],[d_all, derr_all]);
end
