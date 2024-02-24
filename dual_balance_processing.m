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
