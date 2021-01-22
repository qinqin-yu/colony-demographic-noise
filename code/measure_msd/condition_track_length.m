function [xcond, ycond] = condition_track_length(x, y, minlength, maxlength)

todelete = [];

for tracknum = 1:size(x,2)
    xcur = x(:,tracknum);
    ycur = y(:,tracknum);
    xcur_ori = xcur;
    %Get rid of NaN values
    xcur = xcur(~isnan(xcur_ori));
    ycur = ycur(~isnan(xcur_ori));
    if isempty(xcur)
        todelete = [todelete tracknum];
        continue
    end
    %Calculate distance between first and last point for the track
    delx = xcur(end) - xcur(1);
    dely = ycur(end) - ycur(1);
    delr = sqrt(delx^2 + dely^2);
    %If the length of the track is outside of the bounds, remember index to
    %delete track
    if delr<minlength | delr>maxlength
        %Set those tracks to NaN
        todelete = [todelete tracknum];
    end
end

%Delete the tracks that are outside of the length bounds  
xcond = x;
ycond = y;
xcond(:,todelete) = [];
ycond(:,todelete) = [];

end