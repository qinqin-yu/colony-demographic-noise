function [sd, fit_params, xrot_shift, yrot_shift] = calc_sd_continuous(x, y, showplot)
% % Calculates the perpendicular squared displacement of continuous points 
% % x and y from its line of best fit, i.e. the area under the curve.
% % Outputs scalar sd which is the summed squared displacement, 
% % and fit vars y = ax+b.

%Fit the points using a first order polynomial
if length(x)<3
    sd = NaN;
    fit_params = nan(1,2);
    xrot_shift = x;
    yrot_shift = y;
else
fit_params = polyfit(x, y, 1);

%Calculate the fitted points
y_fitted = fit_params(1)*x + fit_params(2); 

%Rotate and shift so that the line of best fit lies on the x axis
theta = -atan2(y_fitted(end)-y_fitted(1), x(end)-x(1));
[xrot, yrot] = rotate_coords(x, y, x(1), y(1), theta);

[xrot_fitted, yrot_fitted] = rotate_coords(x, y_fitted, x(1), y(1), theta);

yrot_shift = yrot - yrot_fitted(1);
xrot_shift = xrot - xrot(1);
yrot_fitted_shift = yrot_fitted - yrot_fitted(1);

%Calculate the squared y values
yrot_shift_sq = yrot_shift.^2;

%Calculate the area under the curve using the trapezoid method
sd = trapz(xrot, yrot_shift_sq);
if sd<0
    sd = NaN;
end

if strcmp(showplot,'yes')
    figure
    hold on
    xverts = [xrot_shift(1:end-1); xrot_shift(1:end-1); xrot_shift(2:end); xrot_shift(2:end)];
    yverts = [zeros(1,size(yrot_shift,2)-1); yrot_shift_sq(1:end-1); yrot_shift_sq(2:end); zeros(1,size(yrot_shift,2)-1)];
    p = patch(xverts,yverts,'b');
    
    plot(x,y, 'b')
    plot(x,y_fitted,'k')
    plot(xrot_shift,yrot_shift_sq, 'g:')
    hold off
end

end
end