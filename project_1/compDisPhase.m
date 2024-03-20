function RawData_dispComp = compDisPhase(RawData, maxDispOrders, dispCoeffs)

ScanPts         = size(RawData, 1);
LinePerFrame    = size(RawData, 2);
kLinear         = linspace(-1,1,ScanPts);
kaxis           = repmat(kLinear',1,LinePerFrame);

RawData_dispComp = RawData;

for i = 1:maxDispOrders-1
    RawData_dispComp = RawData_dispComp.*exp(1j.*(dispCoeffs(i)*(kaxis.^(i+1))));
end

end

