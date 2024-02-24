function RawData_rescaled = reSampling_LUT(RawData, LUT)

RawData_real        = real(RawData);
RawData_imag        = imag(RawData);
RawData_rescaled    = zeros(size(RawData));    
linSampIdx          = linspace(1,size(RawData,1),size(RawData,1));

for i = 1:size(RawData,2)
    RawData_rescaled(:,i) = (interp1(linSampIdx,RawData_real(:,i),LUT,'spline'))...
        +1j.*(interp1(linSampIdx,RawData_imag(:,i),LUT,'spline'));
end

end


