cplxVol = 'OCT_BM_DBD';
load(cplxVol)

usfac = 1;
numAscans = size(OCT_BM_DBD,2);
numBscans = size(OCT_BM_DBD,3);
numMscans = numBscans/numAscans;
local_axial_motion = zeros([numBscans 1]);
local_lateral_motion = zeros([numBscans 1]);
OCT_BM_DBD_mcorr = OCT_BM_DBD;

% aligning chunks of M scans
for I = 1:numMscans:numBscans
    for J = 1:numMscans-1
        [output, ~] = dftregistration(fft2(20.*log10(abs(OCT_BM_DBD(:, :, I)))),...
        fft2(20.*log10(abs(OCT_BM_DBD(:, :, I+J)))), usfac);
        local_axial_motion(I+J) = round(output(3));
        local_lateral_motion(I+J) = round(output(4));
        OCT_BM_DBD_mcorr(:, :, I+J) = circshift(OCT_BM_DBD(:, :, I+J), [round(output(3)) round(output(4))]);
    end
end

for I = 1:numMscans:numBscans
    K = ((I-1)/numMscans)+1;
    BMscan = OCT_BM_DBD_mcorr(:,:,I:I+(numMscans-1));
    BMscan_sub = OCT_BM_DBD_mcorr(:,:,I:I+(numMscans-1)-1);

    % subtract the scans. 
    for J = 1:numMscans-1
        Xconj = BMscan(:,:,1+J).*conj(BMscan(:,:,1));
        BulkOff = repmat(angle(sum(Xconj,1)), [size(Xconj,1) 1]);
        
        BMscan(:,:,1+J) = BMscan(:,:,1+J).*exp(-1j*BulkOff);
        BMscan_sub(:,:,J) = BMscan(:,:,1) - BMscan(:,:,1+J);
    end

    AVG_1(:,:,K) = BMscan(:,:,1);
    VAR_1(:,:,K) = abs(BMscan(:,:,1));
    VAR_2(:,:,K) = abs(var(BMscan(:,:,1:2),0,3));
    VAR_3(:,:,K) = abs(var(BMscan(:,:,1:3),0,3));
    VAR_4(:,:,K) = abs(var(BMscan(:,:,1:4),0,3));
%   
    VAR_12(:,:,K) = abs(var(BMscan(:,:,[1 2]),0,3));
    VAR_13(:,:,K) = abs(var(BMscan(:,:,[1,3]),0,3));
    VAR_14(:,:,K) = abs(var(BMscan(:,:,[1,4]),0,3));
%     avgOCT_2(:,:,K) = mean(BMscan(:,:,1:2),3);
%     VAR_2(:,:,K) = abs(var(BMscan(:,:,1:2),0,3));
% 
%     avgOCT_3(:,:,K) = mean(BMscan(:,:,1:3),3);
%     VAR_3(:,:,K) = abs(var(BMscan(:,:,1:3),0,3));
% 
%     avgOCT_4(:,:,K) = mean(BMscan(:,:,1:4),3);
%     VAR_4(:,:,K) = abs(var(BMscan(:,:,1:4),0,3));
    % Substraction %
%     SUB(:,:,K) = mean(abs(BMscan_sub),3);
%     % Complex Differential Variance %
%     SUB_1(:,:,K) = mean(abs(BMscan_sub(:,:,1)), 3);
%     SUB_2(:,:,K) = mean(abs(BMscan_sub(:,:,1:2)), 3);
%     SUB_3(:,:,K) = mean(abs(BMscan_sub(:,:,1:3)), 3);
%     
%     SUB_1_2(:,:,K) = mean(abs(BMscan(:,:,1) - BMscan(:,:,2)), 3);
%     SUB_1_3(:,:,K) = mean(abs(BMscan(:,:,1) - BMscan(:,:,3)), 3);
%     SUB_1_4(:,:,K) = mean(abs(BMscan(:,:,1) - BMscan(:,:,4)), 3);

%     CDV_1(:,:,K) = abs(BMscan_sub(:,:,1));
%     CDV_2(:,:,K) = abs(var(BMscan_sub(:,:,1:2),0,3));
%     CDV_3(:,:,K) = abs(var(BMscan_sub(:,:,1:3),0,3));
end

correction_source = AVG_1;
usfac = 1;
numFrames = size(correction_source, 3);
global_axial_motion = zeros([numFrames 1]);
OCT_mcorr = correction_source;

for I = 1:numFrames
    %%% Every 'for' loop, reference frame will be the middle frame %%%
    [output, ~] = dftregistration(fft2(20.*log10(abs(correction_source(:, :, round(numFrames./2))))),...
    fft2(20.*log10(abs(correction_source(:, :, I)))), usfac);
    %%% Assign and save the shifting value for axial (yShift) %%%
    global_axial_motion(I) = round(output(3));
    OCT_mcorr(:, :, I) = circshift(correction_source(:, :, I), [output(3), 0]);
end

numLines = size(OCT_mcorr, 2);
global_axial_tilt = zeros([numLines 1]);
OCT_tcorr = OCT_mcorr;
OCT_tcorr_fit = OCT_mcorr;
for I = 1:numLines
    %%% Every 'for' loop, reference frame will be the middle frame %%%
    [output, ~] = dftregistration(fft2(squeeze(abs(OCT_mcorr(:, round(numLines./2), :)))),...
    fft2(squeeze(abs(OCT_mcorr(:, I, :)))), usfac);
    %%% Assign and save the shifting value for lateral (xShift) and axial (yShift) %%%
    global_axial_tilt(I) = round(output(3));
%     OCT_tcorr(:, I, :) = circshift(OCT_mcorr(:, I, :), [output(3), 0]);
end

x = [1:numLines]';
cx = polyfit(x,global_axial_tilt,2);
global_axial_tilt_fit = polyval(cx, x);

% apply to other images as needed 
uncorrected = {VAR_1, VAR_2, VAR_3, VAR_4, VAR_12, VAR_13, VAR_14};
corrections = {1, 2, 3, 4, 5, 6, 7};

for j = 1:length(uncorrected)
    % using poly fit for tilt correction instead

    target = uncorrected{j};
    for I = 1:numLines
        shift_amount = round(global_axial_tilt_fit(I));

        % global motion
        target(:, :, I) = circshift(target(:, :, I), [global_axial_motion(I), 0]);
        % tilt
        target(:, I, :) = circshift(target(:, I, :), [shift_amount, 0]);
    end
    corrections{j} = target;
end 

for i = 1:size(OCT_tcorr_fit, 3)
    disp(i);
    OCT_Image = squeeze((corrections{3}(i,:,:)));
    imagesc(imadjust(mat2gray(OCT_Image))); colormap(gray);
    pause(0.01)
end

OCT_Image = squeeze(VAR_14(240,:,:));
imagesc(imadjust(mat2gray(OCT_Image))); colormap(gray);
title("SUB_{(1,4)}")

% N_Mscans = [];
% I = 400; % choose slice
% 
% K = ((I-1)/numMscans)+1;
% BMscan = OCT_BM_DBD_mcorr(:,:,I:I+(numMscans-1));
% BMscan_sub = OCT_BM_DBD_mcorr(:,:,I:I+(numMscans-1)-1);
% 
% subtract the scans. 
% for J = 1:numMscans-1
%     Xconj = BMscan(:,:,1+J).*conj(BMscan(:,:,1));
%     BulkOff = repmat(angle(sum(Xconj,1)), [size(Xconj,1) 1]);
%     BMscan(:,:,1+J) = BMscan(:,:,1+J).*exp(-1j*BulkOff);
%     BMscan_sub(:,:,J) = BMscan(:,:,1) - BMscan(:,:,1+J);
% end
% 

figure;
for v = 350:360
    clf reset    
    hold on;
%     SUB_1_tcorr_slice = squeeze(CDV_1_tcorr(:,500,:));
%     SUB_2_tcorr_slice = squeeze(CDV_2_tcorr(:,500,:));
%     SUB_3_tcorr_slice = squeeze(CDV_3_tcorr(:,500,:));
    VAR_4_tcorr_slice = squeeze(VAR_4_tcorr(245,:,:));
%     
    vslice = 353
    plot(1:500,  20*log10(abs(SUB_1_tcorr_slice(:,vslice))))
    plot(1:500,  20*log10(abs(SUB_2_tcorr_slice(:,vslice))))
    plot(1:500,  20*log10(abs(SUB_3_tcorr_slice(:,vslice))))
%     plot(1:500,  20*log10(abs(VAR_4_tcorr_slice(:,vslice))))
    legend("(1)","(1,2)" , "(1,2,3)", "(1,2,3,4)")
    
    xlim([0 500])
    ylim([10 120])
    xlabel("Pixel")
    ylabel("Intensity (dB)")

    pause(0.5)
end
%%%

figure;
hold on
hslice = 1;
plot(1:500,  20*log10(abs(VAR_1_tcorr_slice(hslice,:))))
plot(1:500,  20*log10(abs(VAR_2_tcorr_slice(hslice,:))))
plot(1:500,  20*log10(abs(VAR_3_tcorr_slice(hslice,:))))
plot(1:500,  20*log10(abs(VAR_4_tcorr_slice(hslice,:))))
legend("(1)","(1,2)" , "(1,2,3)", "(1,2,3,4)")

xlim([0 500])
xlabel("Pixel")
ylabel("Intensity (dB)")

% Average OCT %
% avgOCT(:,:,K) = mean(BMscan,3);
% Variance %
% VAR(:,:,K) = abs(var(BMscan,0,1));
% Substraction %
% SUB(:,:,K) = mean(abs(BMscan_sub),3);
% Complex Differential Variance %
% CDV(:,:,K) = abs(var(BMscan_sub,0,3));
% 
% OCT_Image = squeeze(20.*log10(abs(octa_1)));
% imagesc(imadjust(mat2gray(OCT_Image))); colormap(gray);
% title("1")
% 
% 
% OCT_Image = squeeze(20.*log10(abs(octa_2)));
% imagesc(imadjust(mat2gray(OCT_Image))); colormap(gray);
% title("1,2")
% 
% OCT_Image = squeeze(20.*log10(abs(octa_3)));
% imagesc(imadjust(mat2gray(OCT_Image))); colormap(gray);
% title("1,2,3")
% 
% 
% OCT_Image = squeeze(20.*log10(abs(octa_4)));
% imagesc(imadjust(mat2gray(OCT_Image))); colormap(gray);
% title("1,2,3,4")





