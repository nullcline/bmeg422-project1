% function PlotVolume(image)
%     for i = 1:size(image, 2)
% %         disp(i);
%         OCT_Image = squeeze(abs(image(:,:,i)).^2);
%         imagesc(imadjust(mat2gray(OCT_Image))); colormap(gray);
%         pause(0.01)
%     end
% end
function PlotVolume(image)
    % Create a figure
    fig = figure;
    % Initial indices
    sliceIdx = 1;
    xCoordIdx = 1;
    
    % Subplot for the original slice
    subplot(1, 2, 1);
    OCT_Image = squeeze(abs(image(:,:,sliceIdx)).^2);
    hImage1 = imagesc(imadjust(mat2gray(OCT_Image))); hold on;
    hLineIndicator = plot([xCoordIdx, xCoordIdx], ylim, 'r', 'LineWidth', 2); % Initial vertical line
    hold off;
    colormap(gray);
    title('Original Slice');
    
    % Subplot for the line plot
    subplot(1, 2, 2);
    OCT_Line = squeeze(20*log10(abs(image(:,xCoordIdx,sliceIdx))));
    hLine = plot(OCT_Line);
    title('Vertical Line Slice at X');
    xlabel('Y-coordinate');
    ylabel('Intensity');
    
    % Slider for the original slice
    slider1 = uicontrol('Parent', fig, 'Style', 'slider', ...
        'Units', 'normalized', 'Position', [0.1 0.01 0.35 0.05], ...
        'value', sliceIdx, 'min', 1, 'max', size(image, 3), ...
        'SliderStep', [1 1]./(size(image, 3)-1), ...
        'Callback', @updateSlice);
    
    % Slider for the X-coordinate
    slider2 = uicontrol('Parent', fig, 'Style', 'slider', ...
        'Units', 'normalized', 'Position', [0.55 0.01 0.35 0.05], ...
        'value', xCoordIdx, 'min', 1, 'max', size(image, 2), ...
        'SliderStep', [1 1]./(size(image, 2)-1), ...
        'Callback', @updateXCoord);
    
    % Callback function to update original slice plot
    function updateSlice(source, ~)
        sliceIdx = round(source.Value);
        disp(sliceIdx)
        OCT_Image = squeeze(abs(image(:,:,sliceIdx)).^2);
        set(hImage1, 'CData', imadjust(mat2gray(OCT_Image)));
        % Update the line plot for new slice
        OCT_Line = squeeze(20*log10(abs(image(:,xCoordIdx,sliceIdx))));
        set(hLine, 'YData', OCT_Line);
        % No need to move the vertical line indicator here since X-coordinate does not change
    end
    
    % Callback function to update X-coordinate and vertical line indicator
    function updateXCoord(source, ~)
        xCoordIdx = round(source.Value);
        disp(xCoordIdx)
        % Update the vertical line plot
        OCT_Line = squeeze(20*log10(abs(image(:,xCoordIdx,sliceIdx))));
        set(hLine, 'YData', OCT_Line);
        % Move the vertical line indicator to the new X-coordinate
        set(hLineIndicator, 'XData', [xCoordIdx, xCoordIdx]);
    end
end

