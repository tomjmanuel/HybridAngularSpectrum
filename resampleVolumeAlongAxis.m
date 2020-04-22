% Tom, 4/20
% function that resamples a transducer mask volume along the axis of the
% transducer to be used with HAS. Also applies sampling to shadow volume
% (CT data)
% Another way to describe what this function does is resample a volume with
% to marked points such that those marked points line up with the first
% dimension in a new volume 

% also returns rots which stores rotation about z and x axis respectively 

% transducer must have focus==1 and top of cap==2

function [nir, shadowr, rots] = resampleVolumeAlongAxis(im,shadow)

    % get starting angle
    I = find(im==1);
    [x y z] = ind2sub(size(im),I(1)); %focus position
    bt = [x y];

    I = find(im==2);
    [x y z] = ind2sub(size(im),I(1)); %top of cap
    tp = [x y];
    clear x y z I

    % this is initial vector
    vec = (bt-tp)./sqrt(sum((bt-tp).^2));

    % find angle btwn this vector and x-axis
    % convert vec components to degrees
    bar1 = acos(vec)*180/pi;
    ang = bar1(1);

    % rotate volume
    foor = imrotate3(im,-ang,[0 0 -1],'nearest','crop');

    %check new angle
    I = find(foor==1);
    [x y z] = ind2sub(size(foor),I(1)); %focus position
    bt = [x y];

    I = find(foor==2);
    [x y z] = ind2sub(size(foor),I(1)); %top of cap
    tp = [x y];
    vec2 = (bt-tp)./sqrt(sum((bt-tp).^2));

    bar2 = acos(vec2)*180/pi;
    'initial'
    round(bar1)
    'rotated'
    round(bar2)
    'angle about Z'
    round(ang)

    % now rotate about x
    %check new angle
    I = find(foor==1);
    [x y z] = ind2sub(size(foor),I(1)); %focus position
    bt = [x z];

    I = find(foor==2);
    [x y z] = ind2sub(size(foor),I(1)); %top of cap
    tp = [x z];
    clear x y z I
    vec3 = (bt-tp)./sqrt(sum((bt-tp).^2));
    bar3 = acos(vec3)*180/pi;
    ang2 = bar3(1);

    % rotate volume
    foor2 = imrotate3(foor,-ang2,[1 0 0],'nearest','crop');
    
    nir = foor2;
    
    
    % now rotate shadow volume
    shadow = imrotate3(shadow,-ang,[0 0 -1],'crop');
    shadowr = imrotate3(shadow,-ang2,[1 0 0],'crop');
    
    % store rotations
    rots = [-ang -ang2];
end