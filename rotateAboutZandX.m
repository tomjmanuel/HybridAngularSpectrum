% function that performs inverse operation following 
% resampleVolumeAlongAxis
% rots is an output of that function containing two rotation angles
% rots = [rotz rotx]
function [ni] = rotateAboutZandX(im,rots)

    foo = imrotate3(im,-rots(2),[1 0 0],'linear','crop');
    ni = imrotate3(foo,-rots(1),[0 0 -1],'linear','crop');

end