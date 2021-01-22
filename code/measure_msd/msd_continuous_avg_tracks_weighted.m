function [msd_avg_weighted, msd_ste_weighted, ensemble_msd, ensemble_ste] = msd_continuous_avg_tracks_weighted(xlinked,ylinked,all_L,mintracks,whichweight)
%Calculates the msd for each track, then averages them together, each
%weighted by 1/(standard error for each track)

%INPUTS
%xlinked, ylinked: coordinates for tracks (columns for different tracks,
%   rows for time)
%mintracks: minimum number of tracks for a window length in order to be
%   averaged

%OUTPUTS
%msd_avg: msd averaged across tracks for each window length
%msd_ste: standard error of the mean of the msd, calculated from all tracks
%   a given window length
%ensemble_msd: values for each track, rows are for each track, columns are 
%   window lengths
%ensemble_ste: standard error of the mean for each track, rows are for each 
%   track, columns are window lengths

xlinked_ori = xlinked;
ylinked_ori = ylinked;

ensemble_msd = nan(size(xlinked,2), length(all_L));
ensemble_ste = nan(size(xlinked,2), length(all_L));
ensemble_N = nan(size(xlinked,2), length(all_L));

%loop through all tracks
for i = 1:size(xlinked, 2)
        %Get non nan values of current track
        xraw = xlinked(:,i);
        yraw = ylinked(:,i);
        xraw = xraw(~isnan(xraw));
        yraw = yraw(~isnan(yraw));

        %If tracks is empty (only nan values) move onto next track
        if size(xraw)~=size(yraw) | isempty(xraw)
            continue
        else
        xraw = xraw(:)';
        yraw = yraw(:)';

        %Calculate MSD of track
        [all_msd, all_ste, all_msd_runavg, all_ste_runavg, all_N_runavg, all_sd_all] = calc_all_L_msd(xraw, yraw, all_L, 'no');

        ensemble_msd(i,:) = all_msd_runavg;
        ensemble_ste(i,:) = all_ste_runavg;
        ensemble_N(i,:) = all_N_runavg;
    end
end

%Average all tracks, weighted by the standard error of the mean of the
%track.
if strcmp(whichweight,'sigma')
    weights = 1./ensemble_ste;
elseif strcmp(whichweight,'N')
    weights = ensemble_N;
end

msd_avg_weighted = nansum((ensemble_msd.*weights))./nansum(weights);

%Only keep windows that have at least thresh number of tracks each
%thresh = 5;
N_runavg = sum(~isnan(ensemble_msd),1);
msd_avg_weighted(N_runavg<mintracks) = NaN;
N_tocount_runavg = N_runavg;
N_tocount_runavg(N_runavg<mintracks) = NaN;

msd_ste_weighted = sqrt((1./(N_tocount_runavg-1)).*(nansum(weights.*((ensemble_msd - msd_avg_weighted).^2),1)./nansum(weights,1)));
end
