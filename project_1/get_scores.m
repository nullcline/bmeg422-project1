function [piqe_score, contrast_score, sharpness_score] = get_scores(fringe_data, ROI)
% image should be N x M array, Ref_Img should always work. assumes grayscale
    
    Ref_FFTData = fft(fringe_data);
    image = 20*log10(abs(Ref_FFTData(ROI(1):ROI(2),:)));

    piqe_score = piqe(image);
    
    contrast_score = (max(image(:)) - min(image(:))) / (max(image(:)) + min(image(:)));
    
    laplacianFilter = fspecial('laplacian');
    img_laplacian = imfilter(image, laplacianFilter, 'replicate');
    sharpness_score = var(img_laplacian(:));

    fprintf('PIQE Score: %f\nContrast Score: %f\nSharpness Score: %f\n', piqe_score, contrast_score, sharpness_score);

    
