% function that pads the medium to improve the spectral support for HAS
% using this instead of padarray so that we can pad the medium with the
% average values of that slice

% inputs:   ps: pad size, this will be added to both ends of the 1st and 2nd
% dimension
%           sos and density: these are the medium matrices

% outputs: sos, density: padded 

function [sos2, den2] = padMedium(ps,sos,density)

    % first get vector which represents average slice property
    msos = squeeze(mean(mean(sos,1),2));
    mden = squeeze(mean(mean(density,1),2));
    
    if size(sos)~=size(density)
        error('medium arrays not same size')
    end
    
    % allocate memory for new media matrices
    dim = size(sos);
    dim(1)=dim(1)+2*ps;
    dim(2)=dim(2)+2*ps;
    sos2 = zeros(dim);
    den2 = zeros(dim);
    
    % next pad slices with their mean
    for i=1:length(msos)
        sos2(:,:,i) = padarray(sos(:,:,i),[ps ps],msos(i));
        den2(:,:,i) = padarray(density(:,:,i),[ps ps],mden(i));
    end
end