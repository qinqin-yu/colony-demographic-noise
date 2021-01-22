function [all_msd, all_ste, all_msd_runavg, all_ste_runavg, all_N_runavg, all_sd_all] = calc_all_L_msd(xraw, yraw, all_L, showplot)

% %Calculates the mean squared displacement for a range of window sizes
% %(all_L) on a trajectory. Two methods: non-overlapping windows (all_msd,
% %all_ste) and overlapping windows (all_msd_runavg, all_ste, runavg).

%Roughly rotate the track so that it lies on the x axis with the first
%point on the y axis and the last point on the positive x axis, and the
%line of best fit on the positive y axis.

%Create variable to save mean squared displacement in
all_msd = nan(1,length(all_L));
all_ste = nan(1,length(all_L));

all_msd_runavg = nan(1,length(all_L));
all_ste_runavg = nan(1,length(all_L));
all_N_runavg = nan(1,length(all_L));

all_sd_all = cell(1,length(all_L));

%Check that there are at least 3 points on the track
if length(xraw)<3
    return
else
    %First fit line of best fit through data points
    %fit_params_raw = polyfit(xraw, yraw, 1);
    xraw_old = xraw;
    yraw_old = yraw;
    [fit_params_noflip, S_noflip] = polyfit(xraw, yraw, 1);
    [fit_params_flip, S_flip] = polyfit(yraw, xraw, 1); %flipped across y = x
    [y_fitted_raw_noflip,delta_noflip] = polyval(fit_params_noflip,xraw,S_noflip);
    [y_fitted_raw_flip,delta_flip] = polyval(fit_params_flip,yraw,S_flip);
    if mean(delta_noflip)<=mean(delta_flip)
        %'no flip'
        y_fitted_raw = y_fitted_raw_noflip;
    else
        %'flip'
        xraw = yraw_old;
        yraw = xraw_old;
        y_fitted_raw = y_fitted_raw_flip;
    end

    %fit_params_raw=linortfit(xraw,yraw)
    %y_fitted_raw = fit_params_raw(1)*xraw + fit_params_raw(2); 

    %Rotate so that line of best fit is horizontal, with first point to the
    %left of the last point
    theta_raw = -atan2(y_fitted_raw(end)-y_fitted_raw(1), xraw(end)-xraw(1));
    %theta_raw = -atan(fit_params_raw(1));
    [xraw_rot, yraw_rot] = rotate_coords(xraw, yraw, xraw(1), yraw(1), theta_raw);
    if xraw_rot(end)<xraw_rot(1)
        xraw_rot = -xraw_rot;
    end

    %Shift values so that line of best fit lies on x axis and first point is at
    %x = 0 (, y = anything)
    yvals = yraw_rot - y_fitted_raw(1);
    xvals = xraw_rot - xraw_rot(1);

    %From now on, xvals and yvals are roughly shifted to be on the x axis, with
    %subsequent x extending in the positive x direction, and the first point
    %close to the orgin

    delx_fine = 0.01;

    %Loop through window lengths
    for iL = 1:length(all_L)

    [msd,ste,all_sd,msd_runavg,ste_runavg,all_sd_runavg,xvals_fine,yvals_fine] = calc_msd_continuous(xvals,yvals,delx_fine,all_L(iL));

    all_sd_all{iL} = all_sd_runavg(~isnan(all_sd_runavg))/all_L(iL);
    %all_sd_runavg
    %Calculate mean of all definitions of window and save to the mean square
    %displacement variable
    all_msd(iL) = msd;
    all_ste(iL) = ste;

    all_msd_runavg(iL) = msd_runavg;
    all_ste_runavg(iL) = ste_runavg;
    all_N_runavg(iL) = length(all_sd_runavg(~isnan(all_sd_runavg)));
    end

    if strcmp(showplot,'yes')
        %Plot to check that it looks right
        figure
        hold on
        %sgtitle('simulated gaussian random walk, x\in[1,100], dx = 1, \mu_x = 2, \mu_y = 0, \sigma_x = \sigma_y = 1\newline integral method', 'FontSize', 15)
        subplot(1,2,1)
        hold on
        plot(xraw_old, yraw_old, 'm')
        plot(xraw, yraw, 'k')
        plot(xraw_rot, yraw_rot, 'r')
        plot(xraw, y_fitted_raw, 'g')
        %plot(xvals, yvals)
        %plot(xvals, yvals)
        xlabel('x')
        ylabel('y')
        set(gca, 'FontSize', 15)
        hold off

        subplot(1,2,2)
        hold on
        errorbar(all_L,all_msd,all_ste)
        errorbar(all_L,all_msd_runavg,all_ste_runavg)
        %plot(all_L, all_L/10, 'k')
        set(gca, 'XScale', 'log', 'YScale', 'log', 'FontSize', 15)
        xlabel('L')
        ylabel('MSD')
        legend('non-overlapping windows', 'overlapping windows', 'MSD = L')
        hold off
    end
end
