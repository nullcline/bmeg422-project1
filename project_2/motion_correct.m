function corrected = motion_correct(oct, global_axial_motion, global_axial_tilt)
    oct_mcorr = oct;
    corrected = oct;
    for I = 1:size(oct, 3) 
        oct_mcorr(:,:,I) = circshift(oct(:, :, I), [global_axial_motion(I), 0]);
    end

    for I = 1:size(oct, 2)
        corrected(:, I, :) = circshift(oct_mcorr(:, I, :), [round(global_axial_tilt(I), 0)]);
    end
end