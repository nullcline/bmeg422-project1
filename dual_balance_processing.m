loadloc = './data';

fn_A = 'RawOCT_A';
load(fullfile(loadloc,fn_A));
fn_B = 'RawOCT_B';
load(fullfile(loadloc,fn_B));

folder_LUT = './LUT';
fileID = fopen(fullfile(folder_LUT,'LUT_A.bin'),'r');
LUT_A = fread(fileID,'double'); fclose(fileID);
fileID = fopen(fullfile(folder_LUT,'LUT_B.bin'),'r');
LUT_B = fread(fileID,'double'); fclose(fileID);

numPoints = size(RawOCT_A, 1);
numAScans = size(RawOCT_A, 2);
numBscans = size(RawOCT_A, 3);

% dispersion
dispMaxOrder = 4;
coeffRange = 20;
depthROI = [100 350];

Ref_CplxRawOCT_A = hilbert(RawOCT_A(:,:,round(end/2)));
Ref_CplxRawOCT_B = hilbert(RawOCT_B(:,:,round(end/2)));

% check
% Ref_FFTData = fft(Ref_CplxRawOCT);
% Ref_Img = 20*log10(abs(Ref_FFTData(1:1024,:)));
% imagesc(Ref_Img); colormap("gray");

Ref_CplxRawOCT_Rescaled_A = reSampling_LUT(Ref_CplxRawOCT_A, LUT_A);
Ref_CplxRawOCT_Rescaled_B = reSampling_LUT(Ref_CplxRawOCT_B, LUT_B);

% subtract b from a
Ref_CplxRawOCT_Rescaled = Ref_CplxRawOCT_Rescaled_A - Ref_CplxRawOCT_Rescaled_B;

% check
% Ref_FFTData = fft(Ref_CplxRawOCT_Rescaled);
% Ref_Img = 20*log10(abs(Ref_FFTData(1:1024,:)));
% imagesc(Ref_Img); colormap("gray");

% remove median to remove dc
Ref_CplxRawOCT_DCSub = Ref_CplxRawOCT_Rescaled - (repmat(median(real(Ref_CplxRawOCT_Rescaled),2), [1,size(Ref_CplxRawOCT_Rescaled,2)]) ...
    +1j.*repmat(median(imag(Ref_CplxRawOCT_Rescaled),2), [1,size(Ref_CplxRawOCT_Rescaled,2)]));

% check
% Ref_FFTData = fft(Ref_CplxRawOCT_DCSub);
% Ref_Img = 20*log10(abs(Ref_FFTData(1:1024,:)));
% imagesc(Ref_Img); colormap("gray");

[dispCoeffs, ~, ~] = setDispCoeff(Ref_CplxRawOCT_DCSub, depthROI, dispMaxOrder, coeffRange);

Ref_CplxRawOCT_DisComp = compDisPhase(Ref_CplxRawOCT_DCSub, dispMaxOrder, dispCoeffs);

Ref_CplxRawOCT_HanWin = Ref_CplxRawOCT_DisComp.*repmat(hann(size(Ref_CplxRawOCT_DisComp, 1)), [1 size(Ref_CplxRawOCT_DisComp,2)]);

% check
Ref_FFTData = fft(Ref_CplxRawOCT_HanWin);
Ref_Img = 20*log10(abs(Ref_FFTData(21:500,:)));
imagesc(Ref_Img); colormap("gray"); title("HanWin with Dual Balance");

if true
    % going along each frame
    disp("initial processing, no tilt corr")

    for FrameNum = 1:numBscans
        RawOCT_A_Frame = RawOCT_A(:,:,FrameNum);
        RawOCT_B_Frame = RawOCT_B(:,:,FrameNum);
        
        rescaled_A = reSampling_LUT(RawOCT_A_Frame, LUT_A);
        rescaled_B = reSampling_LUT(RawOCT_B_Frame, LUT_B);
        
        rescaled = rescaled_A - rescaled_B;

        dc_sub = rescaled - (repmat(median(real(rescaled),2), [1,size(rescaled,2)]) ...
        +1j.*repmat(median(imag(rescaled),2), [1,size(rescaled,2)]));
        
        % no need to estimate, just apply
        dis_comp = compDisPhase(dc_sub, dispMaxOrder, dispCoeffs);
        
        han_win = dis_comp.*repmat(hann(size(dis_comp, 1)), [1 size(dis_comp,2)]);
    
        FFTData = fft(han_win);
    
        OCT(:,:,FrameNum) = flipud(FFTData(depthROI(1):depthROI(2),:));
    end
end


% % idk yet
usfac = 1;
numFrames = size(OCT, 3);
global_axial_motion = zeros([numFrames 1]);
OCT_mcorr = OCT;

for I = 1:numFrames
    %%% Every 'for' loop, reference frame will be the middle frame %%%
    [output, ~] = dftregistration(fft2(20.*log10(abs(OCT(:, :, round(numFrames./2))))),...
    fft2(20.*log10(abs(OCT(:, :, I)))), usfac);
    %%% Assign and save the shifting value for axial (yShift) %%%
    global_axial_motion(I) = round(output(3));
    OCT_mcorr(:, :, I) = circshift(OCT(:, :, I), [output(3), 0]);
end


% Save .tiff stack %
for i = 1:size(OCT_mcorr,3)
    img = imadjust(mat2gray(20.*log10(abs(OCT_mcorr(:,:,i)))));
    imwrite(img, 'OCT_dualbalance_mcorr.tiff', 'WriteMode', 'append', 'Compression','none');
end

numLines = size(OCT_mcorr, 2);
global_axial_tilt = zeros([numLines 1]);
OCT_tcorr = OCT_mcorr;
OCT_tcorr_fit = OCT_mcorr;

for I = 1:numLines
    %%% Every 'for' loop, reference frame will be the middle frame %%%
    [output, ~] = dftregistration(fft2(squeeze(OCT_mcorr(:, round(numLines./2), :))),...
    fft2(squeeze(OCT_mcorr(:, I, :))), usfac);
    %%% Assign and save the shifting value for lateral (xShift) and axial (yShift) %%%
    global_axial_tilt(I) = round(output(3));
    OCT_tcorr(:, I, :) = circshift(OCT_mcorr(:, I, :), [output(3), 0]);
end

x = [1:numLines]';
cx = polyfit(x,global_axial_tilt,2);
global_axial_tilt_fit = polyval(cx, x);

% using poly fit for tilt correction instead
for I = 1:numLines
    shift_amount = round(global_axial_tilt_fit(I));
        OCT_tcorr_fit(:, I, :) = circshift(OCT_mcorr(:, I, :), [shift_amount, 0]);
end


% Save .tiff stack %
for i = 1:size(OCT_tcorr,3)
    img = imadjust(mat2gray(20.*log10(abs(OCT_tcorr(:,:,i)))));
    imwrite(img, 'OCT_dualbalance_tcorr.tiff', 'WriteMode', 'append', 'Compression','none');
end

for i = 1:size(OCT_tcorr_fit,3)
    img = imadjust(mat2gray(20.*log10(abs(OCT_tcorr_fit(:,:,i)))));
    imwrite(img, 'OCT_dualbalance_tcorr_fit.tiff', 'WriteMode', 'append', 'Compression','none');
end


return
