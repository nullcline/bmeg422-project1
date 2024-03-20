function PlotEnfaceAmp(image, StartScan, EndScan)
    depthROI = [StartScan, EndScan];

    image = squeeze(mean(image(depthROI(1):depthROI(2),:,:),1));
    
    figure();
    imagesc(imadjust(mat2gray(image))); colormap("gray"); title('avgOCT');
end