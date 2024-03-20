function plot_snr(images, sliceIdx, xCoordIdx)
    % Validate inputs
    if ~iscell(images) || isempty(images)
        error('Images must be provided in a non-empty cell array.');
    end
    numImages = numel(images);
    
    % Create a figure for the line graphs
    figure;
    hold on; % Hold on to superimpose line graphs
    colors = lines(numImages); % Get a set of distinct colors
    
    % Iterate through each image in the cell array
    for i = 1:numImages
        currentImage = images{i};
        if sliceIdx <= size(currentImage, 3) && xCoordIdx <= size(currentImage, 2)
            % Extract the line (intensity profile) for the given sliceIdx and xCoordIdx
            OCT_Line = squeeze(20*log10(abs(currentImage(:, xCoordIdx, sliceIdx))));
            % Plot the line graph for the current image
            plot(OCT_Line, 'LineWidth', 1, 'Color', colors(i,:));
        else
            warning(['Slice index or X-coordinate index out of bounds for image ', num2str(i)]);
        end
    end
    
    hold off; % Release hold
    title('Superimposed Line Graphs for Given Slice and X-Coordinate');
    xlabel('Y-coordinate');
    ylabel('Intensity');
    legend(arrayfun(@(x) ['Image ' num2str(x)], 1:numImages, 'UniformOutput', false), 'Location', 'Best');
end
