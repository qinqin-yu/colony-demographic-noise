%% Calculate multiple timepoints from single timestep data
function [xlinked, ylinked] = link_multt(xfinal_tot, yfinal_tot)
%load('wo_antibiotic/replicate2/3815.mat');

%load(mat_file);

ti = 2;

xcurrent0 = xfinal_tot{1, 1}';
xcurrent0 = xcurrent0(~isnan(xcurrent0));
xcurrent0 = reshape(xcurrent0, [2,length(xcurrent0)/2]);

ycurrent0 = yfinal_tot{1, 1}';
ycurrent0 = ycurrent0(~isnan(ycurrent0));
ycurrent0 = reshape(ycurrent0, [2,length(ycurrent0)/2]);

xlinked = xcurrent0;
xlinked = [xlinked; NaN(1, size(xlinked, 2))];

ylinked = ycurrent0;
ylinked = [ylinked; NaN(1, size(ylinked, 2))];

%figure
%hold on
%plot(xcurrent0, ycurrent0, 'b', 'LineWidth', 2);

for t=ti:size(xfinal_tot, 1)
    xcurrent = xlinked(t, :);
    xnext = xfinal_tot{t,1}';
    
    ycurrent = ylinked(t, :);
    ynext = yfinal_tot{t,1}';
    
    %nansum(xnext)
    if nansum(nansum(xnext))>0
    %xcurrent = xcurrent(~isnan(xcurrent));
    %xcurrent = reshape(xcurrent, [2,length(xcurrent)/2]);
    %size(xnext)
    %xnext
    xnext = xnext(~isnan(xnext));
    %size(xnext)
    xnext = reshape(xnext, [2,length(xnext)/2]);
    
    %ycurrent = ycurrent(~isnan(ycurrent));
    %ycurrent = reshape(ycurrent, [2,length(ycurrent)/2]);
    ynext = ynext(~isnan(ynext));
    ynext = reshape(ynext, [2,length(ynext)/2]);

    %plot(xnext, ynext, 'b', 'LineWidth', 2);

    for i = 1:length(xcurrent)
        %i
        [I, J] = find(xnext(1, :) == xcurrent(i));
        [Iy, Jy] = find(ynext(1, :) == ycurrent(i));
        if isnan(xcurrent(i))
            xlinked(t+1, i) = NaN;
            ylinked(t+1, i) = NaN;
        elseif J~=Jy
            'could not match xy coordinates'
        elseif isempty(I)
            xlinked(t+1, i) = NaN;
            ylinked(t+1, i) = NaN;
        else
            %xcurrent(2, i)
            xlinked(t+1, i) = xnext(2,J);
            xnext(:,J) = NaN;
            
            ylinked(t+1, i) = ynext(2,J);
            ynext(:,J) = NaN;
            %ynext(:,J) = NaN;
        end
    end
    new_x = xnext(~isnan(xnext)); 
    new_x = reshape(new_x, [2,length(new_x)/2]);
    
    new_y = ynext(~isnan(ynext)); 
    new_y = reshape(new_y, [2,length(new_y)/2]);

    current_size = size(xlinked);
    %size(xlinked)
    xlinked = [xlinked NaN(current_size(1), size(new_x, 2))];
    ylinked = [ylinked NaN(current_size(1), size(new_y, 2))];

    xlinked(t:t+1, current_size(2) + 1:current_size(2) + size(new_x, 2)) = new_x;
    ylinked(t:t+1, current_size(2) + 1:current_size(2) + size(new_y, 2)) = new_y;
    
    else
        xlinked(t+1, 1:length(xcurrent)) = NaN;
        ylinked(t+1, 1:length(xcurrent)) = NaN;
    end
    
end

%plot(xlinked, ylinked, 'r', 'LineWidth', 1)
end

