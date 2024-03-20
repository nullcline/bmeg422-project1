function [dispCoeffs, arrDispCoeff_output, arrCost_output] = setDispCoeff(RawData, depthROI, maxDispOrders, coeffRange)

dispCoeffs = zeros(1, maxDispOrders-1);

for i = 1 : length(dispCoeffs)
    arrDispCoeffRng = [coeffRange, -1 * coeffRange];
    arrCost = zeros(size(arrDispCoeffRng));

    for j = 1 : length(arrDispCoeffRng)
        dispCoeffs(i) = arrDispCoeffRng(j);
        arrCost(j)    = calCostFunc(RawData, depthROI, maxDispOrders, dispCoeffs);
    end
    
    for k = 1 : 50
        [~, idx]       = sort(arrCost);

        dispCoeff_min1 = arrDispCoeffRng(idx(1)); % index for the first minimum cost value
        dispCoeff_min2 = arrDispCoeffRng(idx(2)); % index for the second minimum cost value        

        dispCoeff_new   = (dispCoeff_min1 + dispCoeff_min2)/2;
        arrDispCoeffRng = [arrDispCoeffRng, dispCoeff_new];
        
        dispCoeffs(i) = dispCoeff_new;
        cost_new      = calCostFunc(RawData, depthROI, maxDispOrders, dispCoeffs);
        arrCost       = [arrCost, cost_new];
     end

    [~, argmin]   = min(arrCost);
    dispCoeffs(i) = arrDispCoeffRng(argmin);
    
    arrCost_output(i,:)      = arrCost;
    arrDispCoeff_output(i,:) = arrDispCoeffRng;
    
end

end
