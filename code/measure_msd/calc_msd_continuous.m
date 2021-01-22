function [msd,ste,all_sd,msd_runavg, ste_runavg, all_sd_runavg, xvals_fine,yvals_fine] = calc_msd_continuous(xvals,yvals,delx_fine,L)
%Loops through possible definitions of the window and calculates mean
%squared displacement. 

    if length(unique(xvals))<length(xvals)
        [xvals, ixa, ~] = unique(xvals);
        yvals = yvals(ixa);
    end
    
    %Create more fine x and y values for calculating window length
    xvals_fine = xvals(1):delx_fine:xvals(end);
    yvals_fine = interp1(xvals,yvals,xvals_fine);

    %Calculate the arc length between the 1st and nth point for the fine data
    delx_fine = diff(xvals_fine);
    dely_fine = diff(yvals_fine);
    delr_fine = sqrt(delx_fine.^2 + dely_fine.^2);
    totdelr_fine = [0 cumsum(delr_fine)];

    %Calculate number of windows
    numwin = floor(totdelr_fine(end)/(2*L));

    %Initiate variable for saving the squared displacement from each window
    all_sd = nan(1,numwin);

    %Loop through each definition of the window
    for k = 1:numwin

    %Find indices of x and y fine data that fall in window
    i = find(totdelr_fine>(k-1)*L,1);
    j = find(totdelr_fine>(k*L),1);

    %Save those fine data into a new variable
    xint = xvals_fine(i-1:j-1);
    yint = yvals_fine(i-1:j-1);

    %Also calculate range of values that are included
    xmin = xvals_fine(i-1);
    xmax = xvals_fine(j-1);
    ymin = yvals_fine(i-1);
    ymax = yvals_fine(j-1);

    %Pick out points from original data that fall in this range, and add the
    %fine end points (in order to get to target window length)
    x = [xint(1) xvals((xvals>=xmin) & (xvals<=xmax)) xint(end)];
    y = [yint(1) yvals((xvals>=xmin) & (xvals<=xmax)) yint(end)];

    [sd, fit_params, xrot_shift, yrot_shift] = calc_sd_continuous(x, y, 'no');

    %Divide by the length of the window and save to variable that contains all
    %definitions of the window
    all_sd(k) = sd;       %SHOULD USE ACTUAL LENGTH OF WINDOW?
    end

    if length(all_sd)>1
        msd = mean(all_sd)/L;
        ste = std(all_sd)/(sqrt(numwin)*L);
    else
        msd = NaN;
        ste = NaN;
    end

    all_sd_runavg = nan(1,length(xvals));

    for l = 1:length(xvals)-1

        xvals_fine_runavg = xvals(l):delx_fine:xvals(end);
        yvals_fine_runavg = interp1(xvals(l:end),yvals(l:end),xvals_fine_runavg);

        %Calculate the arc length between the 1st and nth point for the fine data
        delx_fine_runavg = diff(xvals_fine_runavg);
        dely_fine_runavg = diff(yvals_fine_runavg);
        delr_fine_runavg = sqrt(delx_fine_runavg.^2 + dely_fine_runavg.^2);
        totdelr_fine_runavg = [0 cumsum(delr_fine_runavg)];

        if isempty(find(totdelr_fine_runavg>L,1))
            break
        else
            %Find indices of x and y fine data that fall in window
            j = find(totdelr_fine_runavg>L,1);

            %Also calculate range of values that are included
            xmax = xvals_fine_runavg(j-1);

            %Pick out points from original data that fall in this range, and add the
            %fine end points (in order to get to target window length)
            x = [xvals((xvals>=xvals(l)) & (xvals<=xmax)) xvals_fine_runavg(j-1)];
            y = [yvals((xvals>=xvals(l)) & (xvals<=xmax)) yvals_fine_runavg(j-1)];

            [sd, fit_params, xrot_shift, yrot_shift] = calc_sd_continuous(x, y, 'no');

            all_sd_runavg(l) = sd;
        end
    end

    N = sum(~isnan(all_sd_runavg));
    if N>1
        msd_runavg = nanmean(all_sd_runavg)/L;
        ste_runavg = nanstd(all_sd_runavg)/(sqrt(N)*L);
    else
        msd_runavg = NaN;
        ste_runavg = NaN;
    end
end
