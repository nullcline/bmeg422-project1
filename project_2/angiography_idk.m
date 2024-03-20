clear all; 
clc;
close all;
load(fullfile('OCT_BM_DBD.mat'));

OCT=OCT_BM_DBD;

%% Local Motion Correction 
usfac = 1; 
numMscans = 4;
numBscans = size(OCT,3); 
local_axial_motion = zeros([numBscans 1]); 
local_lateral_motion = zeros([numBscans 1]);

for I = 1:numMscans:numBscans
    for J = 1:numMscans-1
        [output, ~] = dftregistration(fft2(20.*log10(abs(OCT_BM_DBD(:, :, I)))),...
        fft2(20.*log10(abs(OCT_BM_DBD(:, :, I+J)))), usfac);
        local_axial_motion(I+J) = round(output(3));
        local_lateral_motion(I+J) = round(output(4));
        OCT(:, :, I+J) = circshift(OCT_BM_DBD(:, :, I+J), [round(output(3)) round(output(4))]);
    end
end

clear OCT_BM_DBD;

VAR_2 = zeros(500,500,500);
VAR_3 = zeros(500,500,500);
VAR_4 = zeros(500,500,500);

% VAR_4_no_bulk = zeros(500, 500 ,500);
%% OCT-A Processing
for I = 1:numMscans:numBscans
    K = ((I-1)/numMscans)+1;
    BMscan = OCT(:,:,I:I+(numMscans-1));
    BMscan_sub = OCT(:,:,I:I+(numMscans-1)-1);
%     BMscan_no_bulk = BMscan;
%     subtract the scans. 
    for J = 1:numMscans-1
        Xconj = BMscan(:,:,1+J).*conj(BMscan(:,:,1));
        BulkOff = repmat(angle(sum(Xconj,1)), [size(Xconj,1) 1]);
        
        BMscan(:,:,1+J) = BMscan(:,:,1+J).*exp(-1j*BulkOff);
        BMscan_sub(:,:,J) = BMscan(:,:,1) - BMscan(:,:,1+J);
    end
    
%   AVG_1_no_bulk_off(:,:,K) = BMscan(:,:,1);
    VAR_1(:,:,K) = abs(var(BMscan(:,:,1),0,3));
%     VAR_2(:,:,K) = abs(var(BMscan(:,:,1:2),0,3));
%     VAR_3(:,:,K) = abs(var(BMscan(:,:,1:3),0,3));
%     VAR_4(:,:,K) = abs(var(BMscan(:,:,1:4),0,3));
%     VAR_4_no_bulk(:,:,K) = abs(var(BMscan_no_bulk(:,:,1:4),0,3));


%     SUB_1_bulk(:,:,K) = mean(abs(BMscan_sub(:,:,1)), 3);
%     SUB_2(:,:,K) = mean(abs(BMscan_sub(:,:,1:2)), 3);
%     SUB_3(:,:,K) = mean(abs(BMscan_sub(:,:,1:3)), 3);
  
%     CDV_1(:,:,K) = abs(BMscan_sub(:,:,1));
%     CDV_2(:,:,K) = abs(var(BMscan_sub(:,:,1:2),0,3));
%     CDV_3(:,:,K) = abs(var(BMscan_sub(:,:,1:3),0,3));
end 
%% Global Motion correction

VAR_1 = abs(BMscan(:,:,1:4:end));

numFrames = size(AVG_1, 3); global_axial_motion = zeros([numFrames, 1]); 
for I = 1:numFrames 
    % Every 'for' loop, reference frame will be the middle frame
    [output, ~] = dftregistration(fft2(20.*log10(abs(AVG_1(:, :, round(numFrames./2))))),fft2(20.*log10(abs(AVG_1(:, :, I)))), usfac);
    % Assign and save the shifting value for axial (yShift) 
    global_axial_motion(I) = round(output(3)); 

%     AVG_mcorr(:, :, I) = circshift(AVG_1(:, :, I), [output(3), 0]); 
% %     VAR_1_mcorr(:, :, I) = circshift(VAR_1(:, :, I), [output(3), 0]);
%     VAR_2_mcorr(:, :, I) = circshift(VAR_2(:, :, I), [output(3), 0]); 
%     VAR_3_mcorr(:, :, I) = circshift(VAR_3(:, :, I), [output(3), 0]); 
%     VAR_4_mcorr(:, :, I) = circshift(VAR_4(:, :, I), [output(3), 0]); 
%     VAR_12_mcorr(:, :, I) = circshift(VAR_12(:, :, I), [output(3), 0]); 
%     VAR_13_mcorr(:, :, I) = circshift(VAR_13(:, :, I), [output(3), 0]); 
%     VAR_14_mcorr(:, :, I) = circshift(VAR_14(:, :, I), [output(3), 0]); 
    SUB_no_bulk_mcorr(:,:,I) = circshift(SUB_1(:, :, I), [output(3), 0]); 

end

%% Global Tilt Correction 
numLines = size(AVG_mcorr, 2); 
global_axial_tilt = zeros([numLines 1]); 
for I = 1:numLines 
    % Every 'for' loop, reference frame will be the middle frame %%% 
    [output, ~] = dftregistration(fft2(squeeze(AVG_mcorr(:, round(numLines./2), :))), fft2(squeeze(AVG_mcorr(:, I, :))), usfac); 
% Assign and save the shifting value for lateral (xShift) and axial (yShift) %%%
    global_axial_tilt(I) = round(output(3)); 
end

%% Polynomial fitting the tilt curve 
x = (1:numLines)'; cx = polyfit(x,global_axial_tilt,2); 
global_axial_tilt = polyval(cx, x); 
for I = 1:numLines 
%     AVG_tcorr(:, I, :) = circshift(AVG_mcorr(:, I, :), [round(global_axial_tilt(I), 0)]); 
% %     VAR__1_tcorr(:, I, :) = circshift(VAR_mcorr(:, I, :), [round(global_axial_tilt(I), 0)]); 
%     VAR__2_tcorr(:, I, :) = circshift(VAR_2_mcorr(:, I, :), [round(global_axial_tilt(I), 0)]); 
%     VAR__3_tcorr(:, I, :) = circshift(VAR_3_mcorr(:, I, :), [round(global_axial_tilt(I), 0)]); 
%     VAR__4_tcorr(:, I, :) = circshift(VAR_4_mcorr(:, I, :), [round(global_axial_tilt(I), 0)]); 
%     VAR__12_tcorr(:, I, :) = circshift(VAR_12_mcorr(:, I, :), [round(global_axial_tilt(I), 0)]); 
%     VAR__13_tcorr(:, I, :) = circshift(VAR_13_mcorr(:, I, :), [round(global_axial_tilt(I), 0)]); 
%     VAR__14_tcorr(:, I, :) = circshift(VAR_14_mcorr(:, I, :), [round(global_axial_tilt(I), 0)]); 
    no_bulk_tcorr(:, I, :) = circshift(no_bulk_mcorr(:, I, :), [round(global_axial_tilt(I), 0)]); 
end 


 % stuff to run 

var_1_corr = motion_correct(VAR_1, global_axial_motion, global_axial_tilt);
var_2_corr = motion_correct(VAR_2, global_axial_motion, global_axial_tilt);
var_3_corr = motion_correct(VAR_3, global_axial_motion, global_axial_tilt);
var_4_corr = motion_correct(VAR_4, global_axial_motion, global_axial_tilt);

 PlotEnfaceAmp(VAR__2_tcorr,VAR__3_tcorr,VAR__4_tcorr,VAR__12_tcorr,100,300);

 psnrVal = calculatePSNR(VAR__2_tcorr, VAR__3_tcorr);

% Save .tiff stack %
for i = 1:size(cdv_corr,3)
    img = imadjust(mat2gray((abs(var_corr(:,:,i)))));
    imwrite(img, 'OCT_var_corr.tiff', 'WriteMode', 'append', 'Compression','none');
end

 plot_snr({var_2_corr, var_3_corr, var_4_corr}, 173, 42)

