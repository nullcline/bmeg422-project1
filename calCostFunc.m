function cost = calCostFun(RawData, depthROI, maxDispOrders, dispCoeffs)

RawData_dispComp = compDisPhase(RawData,maxDispOrders,dispCoeffs);
FFTData_dispComp = fft(RawData_dispComp);

OCT      = abs(FFTData_dispComp(depthROI(1):depthROI(2),:)).^2;
NormOCT  = OCT./sum(OCT(:));

entropy  = -1*(NormOCT.*log(NormOCT));
cost     = sum(entropy(:));

end




