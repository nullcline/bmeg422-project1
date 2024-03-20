function plot_oct(fringe_data, ROI, type, plot_title)
    % example usage: 
    %
    % plot_oct(Ref_CplxRawOCT_DCSub, [100 350], "FFT", "DC-Sub")
    %

    Ref_FFTData = fft(fringe_data);

    if type == "FFT"
        disp(size(fringe_data))
        Ref_FFTData = fft(fringe_data);
        Ref_Img = 20*log10(abs(Ref_FFTData));
        imagesc(Ref_Img); colormap("gray");
        xlabel("Position (µm)")
        ylabel("Depth (µm)")
        c = colorbar
        c.Title.String = "Intensity (dB)";
    else
        % plot fringe data, currently defaulted to first slice
        plot(1:2048, abs(Ref_FFTData(1:2048,1)))
        xlim([0 2048])
        xlabel("Wavenumber")
        ylabel("Intensity")
        title(plot_title)
    end