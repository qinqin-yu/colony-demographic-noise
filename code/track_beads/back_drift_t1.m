function [xcorr, ycorr, totdelx, totdely] = back_drift_t1(xplot, yplot, xmin, xmax,...
    ymin, ymax, prevdelx, prevdely, plotfig)

%Subtract off background
t = size(xplot, 1);
delx_all = zeros(1, t-1);
dely_all = zeros(1, t-1);

xint = xplot>xmin & xplot<xmax;
yint = yplot>ymin & yplot<ymax;
% ymin
% ymax
% xint
% yint

xback = xplot;
yback = yplot;
xcorr = xplot;
ycorr = yplot;

%figure
%hold on
%plot(xback,  yback, 'ro');

xback(~xint | ~yint) = NaN;
yback(~xint | ~yint) = NaN;


%plot(xback, yback, 'kx');
%hold off

totdelx = prevdelx;
totdely = prevdely;

if nansum(xback(1,:),2)>0 && nansum(xback(2,:),2)>0
%for i = 2:t
    delx = nanmean(xback(2, :) - xback(1, :));
    dely = nanmean(yback(2, :) - yback(1, :));
    
    delx_all(2) = delx;
    dely_all(2) = dely;
    
    totdelx = totdelx + delx;
    totdely = totdely + dely;
    
    %disp(prevdelx)
    xcorr(1, :) = xplot(1, :) - prevdelx*ones(size(xplot(1, :)));
    ycorr(1, :) = yplot(1, :) - prevdely*ones(size(xplot(1, :)));
    
    xcorr(2, :) = xplot(2, :) - totdelx*ones(size(xplot(2, :)));
    ycorr(2, :) = yplot(2, :) - totdely*ones(size(xplot(2, :)));
    
%end
else 
    xcorr = xplot;
    ycorr = yplot;
end

if plotfig
    figure
    hold on
    plot(xplot, yplot, 'b', 'LineWidth', 2)
    %plot(xback, yback, 'r')
    plot(xcorr, ycorr, 'r', 'LineWidth', 2)
    hold off
end
end

% figure
% dely1 = yplot(2, :)-yplot(1, :);
% dely2 = yplot(3, :)-yplot(2, :);
% dely3 = yplot(4, :)-yplot(3, :);
% histogram(dely1);
% hold on
% histogram(dely2);
% histogram(dely3);