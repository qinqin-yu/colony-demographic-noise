%Extrapolate MSD to delta r = 100 um or 1000 um by using weighted least
%squares linear fitting in log-log space
function [] = get_msd_extrap_front(path_folder, experiment, deltar)
%Add path with curve fitting functions

%Get folders with data on MSD as a function of window length 
msd_folder = [path_folder filesep 'mat' filesep experiment];

%Make folder to store results if it doesn't already exist
if ~isdir([path_folder filesep 'extrapolated'])
    mkdir([path_folder filesep 'extrapolated']);
end

%Figure out number of files (trials) to work through
files = dir2(msd_folder);
nfiles = length(files);

%Natural log values for the window size that we're extrapolating to
logx = log(deltar);

%Iniaiate variables for saving the extrapolated values
yval_all = nan(96,1);

yval_all_err = nan(96, 1);

%Loop through all files
for p=1:nfiles
    filename = files(p).name;
    load([path_folder filesep 'mat' filesep experiment filesep filename]);
    msd_avg_weighted = rho_L;
    msd_ste_weighted = rho_L_ste;
    
    %Check that there are at least 5 MSD data points for proper
    %extrapolation
    numpoints = sum(~isnan(msd_avg_weighted));
    if numpoints > 4
        %Get index of file for saving order later
        indm = strfind(filename, 'm');
        inddot = strfind(filename, '.');
        ind = str2num(filename(indm+1:inddot-1))+1;

        %Convert to natural log space for doing linear fitting
        x = log(all_L(~isnan(msd_avg_weighted)));
        y = log(msd_avg_weighted(~isnan(msd_avg_weighted)));
        sig = msd_ste_weighted(~isnan(msd_avg_weighted))./msd_avg_weighted(~isnan(msd_avg_weighted));
        %Do linear fit in log log space using weighted least squares
        [a,aerr,chisq,yfit] = fitlin(x,y,sig); %Weird that the slopes are bigger than 2 - possibly bug in code?
        
        %Calculate extrapolated values
        log_yval = logx*a(2) + a(1);
        exp(log_yval);
        yval_all(ind) = exp(log_yval);
        
        %And errors
        
        log_yval_err_plus = logx*(a(2) + aerr(2)) + (a(1) + aerr(1));
        log_yval_err_minus = logx*(a(2) - aerr(2)) + (a(1) - aerr(1));
        log_yval_err = (log_yval_err_plus - log_yval_err_minus)/2;
        yval_all_err(ind) = yval_all(ind) * log_yval_err;
        
    end
end

%Write results to csv file
csvwrite([path_folder filesep 'extrapolated' filesep experiment '.csv'],[yval_all, yval_all_err])

end


