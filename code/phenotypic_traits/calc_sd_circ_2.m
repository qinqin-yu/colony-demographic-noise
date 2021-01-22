function [d_all, sd_all, theta_all, radii_all] = calc_sd_circ_2(xvals, yvals, a, b, R)
   % Calculates the perpendicular squared displacement of points xvals and
   % yvals from its line of best fit. Outputs vector sd which contains all
   % of the squared displacements, and fit vars y = ax+b.

  
    d_all = [];
    sd_all = [];
    theta_all = [];
    radii_all = [];
    for j = 1:length(xvals)
        theta = atan2(yvals(j)-b, xvals(j)-a);
        xcirc = R*cos(theta)+a;
        ycirc = R*sin(theta)+b;
        
        d = sqrt((xvals(j) - xcirc)^2 + (yvals(j)-ycirc)^2);
        d_all(j) = d;
        sd_all(j) = d^2;
        theta_all(j) = theta;
        radii_all(j) = sqrt((xvals(j) - a)^2 + (yvals(j) - b)^2);
%         if j == 500
%             figure
%             hold on
%             scatter(xvals(j), yvals(j));
%             scatter(xcirc, ycirc);
%             plot([xvals(j) a], [yvals(j) b]);
%             ang=0:0.01:2*pi; 
%             xp=R*cos(ang);
%             yp=R*sin(ang);
%             scatter(xvals, yvals, 1, 'r.')
%             plot(a+xp,b+yp, 'c');
%             hold off
%         end
    end
    
end