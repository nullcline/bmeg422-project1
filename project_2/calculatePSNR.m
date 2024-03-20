function psnrVal = calculatePSNR(imageA, imageB)
    % Ensure the images are double precision for calculations
    imageA = double(imageA);
    imageB = double(imageB);
    
    % Calculate MSE
    mse = mean((imageA(:) - imageB(:)).^2);
    
    % Assuming the images are 8-bit, MAX_I is 255
    % Adjust MAX_I if your images have a different range
    MAX_I = 255;
    
    % Calculate PSNR
    psnrVal = 20 * log10(MAX_I / sqrt(mse));
end