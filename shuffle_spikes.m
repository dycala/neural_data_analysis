
function [shuff] = shuffle_spikes(ReachS,raster,lagmax)

%% Get shuffled spikes in reg_kin epoch by sampling reg_kin ISIs%get ISIs from raster for each trial for later shuffle 

lag = (lagmax/10);

TrialRaster = [];
for i =1:length(ReachS)
    if ReachS(i).exclude == 0 
        % find data in which you will regress
        idx1 = knnsearch(raster(:,1),ReachS(i).kin_10ms(400,1));
        idx2 = knnsearch(raster(:,1),ReachS(i).kin_10ms(600,1));
        TrialRaster = vertcat(TrialRaster,raster(idx1-lag:idx2-lag,2));
    end
end

% make structures to later put spikes into
% make zero matrices 2010 because they will be downsampled by 10 to 201
for ii=1:100
   shuff(ii).binned = [];
   for i = 1:length(ReachS)
       if ReachS(i).exclude ==0
            shuff(ii).binned =  vertcat(shuff(ii).binned,zeros(2010,1));
       end
   end
      
    % randomly sample ISIs   
    idx=0;     
    while idx < length(shuff(ii).binned)
        rand = randi([1 length(TrialRaster)]);
        ISI = round(TrialRaster(rand,1)*1000); % so that each bin is a ms
        idx = (idx + ISI);
        if idx<length(shuff(ii).binned)
                shuff(ii).binned(idx,1) = 1;
        else 
            break
        end
    end

    % Bin at 10 ms
    a=1;
    iii=1;
    while iii<length(shuff(ii).binned)
        shuff(ii).binned10(a,1) = sum(shuff(ii).binned(iii:iii+9));
        a=a+1;
        iii=iii+10;
    end

    % get shuffled SS_Bin1 and raster for iRate Smooth
    Spikeshuff.SS_Bin1(:,1) = [0.001:0.001:length(shuff(ii).binned)/1000];
    Spikeshuff.SS_Bin1(:,2) = shuff(ii).binned;
    
    Spikeshuff = get_raster(Spikeshuff);

    shuff(ii).binned10smooth = i_rate_smooth(Spikeshuff);

end

% removed binned
shuff = rmfield(shuff,'binned');