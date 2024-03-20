function saveMatrixAsGif(data, filename, delayTime)
% saveMatrixAsGif Saves an NxNxM matrix as a GIF file, with robust normalization.
%   data: The NxNxM input matrix.
%   filename: Name of the output GIF file.
%   delayTime: Delay time between frames in seconds.

[N, ~, M] = size(data); % Get the dimensions of the input matrix

% Calculate the 1st and 99th percentiles as the normalization range
pctLow = prctile(data(:), 1);
pctHigh = prctile(data(:), 99);

for i = 1:M
    % Extract the ith frame
    frame = squeeze(data(:,:,i));
    
    % Normalize the frame based on the percentile range
    % Any value below pctLow is set to 0, and above pctHigh is set to 1
    normalizedFrame = (frame - pctLow) / (pctHigh - pctLow);
    normalizedFrame(normalizedFrame < 0) = 0;
    normalizedFrame(normalizedFrame > 1) = 1;
    
    % Convert the normalized frame to uint8
    frameUint8 = im2uint8(normalizedFrame);
    
    % Convert the frame to an indexed image
    [ind, map] = gray2ind(frameUint8, 256);
    
    % Write to the GIF File
    if i == 1
        imwrite(ind, map, filename, 'gif', 'LoopCount', Inf, 'DelayTime', delayTime);
    else
        imwrite(ind, map, filename, 'gif', 'WriteMode', 'append', 'DelayTime', delayTime);
    end
end
end
