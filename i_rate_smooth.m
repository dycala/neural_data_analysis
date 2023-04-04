function [RateSmooth] = i_rate_smooth(Spikes)

    % get ISI inverse rate
    iRateFilled(:,1) = Spikes.SS_Bin1(:,1);

    [~,idx1] = min(abs(iRateFilled(:,1)-Spikes.SS_Raster(1,1)));

    for i = 2:length(Spikes.SS_Raster)
        idx2 = round((Spikes.SS_Raster(i,1)-Spikes.SS_Bin1(1,1))*1000)+1;
        iRateFilled(idx1:idx2,2) = 1/Spikes.SS_Raster(i,2);
        idx1 = idx2;
    end

    % gaussian
    sigma = .020; %Standard deviation of the kernel = 20
    edges=[-3*sigma:.001:3*sigma]; %Time ranges form -3*st. dev. to 3*st. dev.
    kernel = normpdf(edges,0,sigma); %Evaluate the Gaussian kernel
    kernel = kernel*.001; %Multiply by bin width so the probabilities sum to 1
    s=conv(iRateFilled(:,2),kernel); %Convolve spike data with the kernel
    center = ceil(length(edges)/2); %Find the index of the kernel center
    s=s(center:end-center+1); %Trim out the relevant portion of the spike density


    RateSmooth = [iRateFilled(:,1),s];

    % sample at 10 ms
    RateSmooth = downsample(RateSmooth,10);
end


