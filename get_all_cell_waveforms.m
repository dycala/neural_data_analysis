function [distCutoff, corrCutoff, wfMAT] = get_all_cell_waveforms(path_loc,save_loc)
    
    %% Get waveforms for correlation across experiment

    cd(path)
    directory = dir('*_NPresults.mat'); 
    
    
    wfMAT = [];n = 1;
    for day  = 1:length(directory)
        load(directory(day).name)
    
        
        if ~isfield(cellData,'waveformCat1')
        myKsDir = convertStringsToChars(loc + "\continuous\Neuropix-PXI-100.0")
        % set parameters 
        gwfparams.dataDir = myKsDir;    % KiloSort/Phy output folder
        apD = 'continuous.dat';
        gwfparams.fileName = apD;         % .dat file containing the raw 
        gwfparams.dataType = 'int16';            % Data type of .dat file (this should be BP filtered)
        gwfparams.nCh = 384;                      % Number of channels that were streamed to disk in .dat file
        gwfparams.wfWin = [-40 41];          % Number of samples before and after spiketime to include in waveform
        gwfparams.nWf = 1000;                    % Number of waveforms per unit to pull out
    
        % load data
        sp = loadKSdir(myKsDir);
    
        numClusters = max(sp.clu);
        Tend = max(sp.st);
        Tstart = min(sp.st);
    
        for i =1 : numClusters
            rates(i,1) = length(sp.st(sp.clu==i))/Tend;
            rates(i,2) = i;
        end
    
        pre = [ReachS(:).stim] == 0;
        t1 = find(pre, 1, 'first');
        t2 = find(pre, 1, 'last');
        window1 = [ReachS(t1).out(1,1),ReachS(t2).out(end,1)];
    
        po = [ReachS(:).stim] == 2;
        t1 = find(po, 1, 'first');
        t2 = find(po, 1, 'last');
        window2 = [ReachS(t1).out(1,1),ReachS(t2).out(end,1)];
    
        for cell = 1:length(cellData)
        
            if cellData(cell).fr>10
                [~,idx]=min(abs(rates(:,1)-cellData(cell).fr));
                clusterID = idx;
                bestChan = cellData(cell).channel;
    
                % get number of channels (no bad channels)
                chNum = numel(readNPY(fullfile(gwfparams.dataDir, 'channel_map.npy'))+1);               % Order in which data was streamed to disk; must be 1-indexed for Matlab
    
                if bestChan<17
                    channels = 1:33;
                elseif bestChan>chNum-16
                    channels = chNum-32:chNum;
                else
                    channels = bestChan-16:bestChan+16;
                end
    
                % search baseline and washout epochs for spike templates
                for epoch = 1:2
    
                    allSpikeTimes = ceil(sp.st(sp.clu==clusterID)*30000); % Vector of cluster spike times (in samples) same length as .spikeClusters
                    allSpikeClusters = sp.clu(sp.clu==clusterID);
    
                    if epoch == 1
                        allSpikeTimesHalf = allSpikeTimes(allSpikeTimes>window1(1)*30000 & allSpikeTimes<window1(2)*30000 );
                        allSpikeClustHalf = allSpikeClusters(allSpikeTimes>window1(1)*30000 & allSpikeTimes<window1(2)*30000 );
                    else
                        allSpikeTimesHalf = allSpikeTimes(allSpikeTimes>window2(1)*30000 & allSpikeTimes<window2(2)*30000);
                        allSpikeClustHalf = allSpikeClusters(allSpikeTimes>window2(1)*30000 & allSpikeTimes<window2(2)*30000);
                    end
    
                    gwfparams.spikeTimes = allSpikeTimesHalf; % Vector of cluster spike times (in samples) same length as .spikeClusters
                    gwfparams.spikeClusters = allSpikeClustHalf;
    
                    wf = getWaveForms(gwfparams);
    
                    % get channel averages
                    onChan = squeeze(wf.waveForms(:,:,:,:));
                    chanCat = [];
                    for chanIDs = channels 
                        chanAvg = nanmean(squeeze(onChan(:,chanIDs,:)));
                        chanCat = [chanCat,chanAvg-mean(chanAvg)];
                    end
    
                    chanCatAll = [];
                    for chanIDs = 1:size( wf.waveForms,3)
                        chanAvg = nanmean(squeeze(onChan(:,chanIDs,:)));
                        chanCatAll = [chanCatAll,chanAvg-mean(chanAvg)];
                    end
    
                    if epoch == 1
                        chanCat1 = chanCat;
                        chanCatAll1 = chanCatAll;
                    else
                        chanCat2 = chanCat;
                        chanCatAll2 = chanCatAll;
                    end
                    clear allSpikeTimes allSpikeClust chanCat onChan chanCatAll
    
                end
    
            cellData(cell).waveformCat1 = chanCat1;
            cellData(cell).waveformCat2 = chanCat2;
            cellData(cell).whichChans = channels;
    
            % save waveforms in matrix for total correlations
            wfMAT(n).wf1 = chanCatAll1;
            wfMAT(n).wf2 = chanCatAll2;
            wfMAT(n).whichChans = channels;
            wfMAT(n).day = day;
            n = n+1;
    
            end
        end
    
        % get pearson correlation of waveforms
        for i = 1:length(cellData)
            if ~isempty(cellData(i).waveformCat1) &&length(cellData(i).waveformCat1)== 2706  && length(cellData(i).waveformCat2)== 2706
                cellData(i).pearson = corr(cellData(i).waveformCat1',cellData(i).waveformCat2');
            end
        end
    
        channelsXtot = [[-15,5,-5,15]];
        channelsXtotal = [];
        for i = 1:96
            channelsXtotal = [channelsXtotal,channelsXtot];
        end
    
        channelsYtot = [0,0];
        channelsYtotal = [];
        for i = 1:192
            channelsYtotal = [channelsYtotal,channelsYtot+(i*20)]
        end
    
        % get locations of start and end data
        for cell = 1:length(cellData)
            if ~isempty(cellData(cell).waveformCat1) && length(cellData(cell).waveformCat1)== 2706  && length(cellData(cell).waveformCat2)== 2706
                templateLen= 82;
                channelsY = channelsYtotal(cellData(cell).whichChans); % vertical loc in micrometers;
                channelsX = channelsXtotal(cellData(cell).whichChans); % lateral loc in micrometers;
    
                Xnum = 0; Xdenom = 0; Ynum = 0; Ydenom = 0;
                for channel = 0:32
                    start = (channel*templateLen)+1;
                    stop = start+templateLen-1;
    
                    % start data
                    amp = max(cellData(cell).waveformCat1(start:stop))-min(cellData(cell).waveformCat1(start:stop));
    
                    Xnum = Xnum + (amp^2)*channelsX(channel+1);
                    Xdenom = Xdenom + (amp^2);
    
                    Ynum = Ynum + (amp^2)*channelsY(channel+1);
                    Ydenom = Ydenom + (amp^2);
    
                    cellData(cell).startX = Xnum/Xdenom;
                    cellData(cell).startY = Ynum/Ydenom;
    
                    % for end data
                    amp = max(cellData(cell).waveformCat2(start:stop))-min(cellData(cell).waveformCat2(start:stop));
    
                    Xnum = Xnum + (amp^2)*channelsX(channel+1);
                    Xdenom = Xdenom + (amp^2);
    
                    Ynum = Ynum + (amp^2)*channelsY(channel+1);
                    Ydenom = Ydenom + (amp^2);
    
                    cellData(cell).endX = Xnum/Xdenom;
                    cellData(cell).endY = Ynum/Ydenom;
    
    
                end
            end
        end
    
        end
    end
    

    %% get correlations shuffle
    shuffcorr = [];
    templateLen = 82;
    for cell = 1:length(wfMAT)
        channelST = wfMAT(cell).whichChans(1)-1;
        channelEN = wfMAT(cell).whichChans(end)-1;
    
        start = (channelST*templateLen)+1;
        stop = ((channelEN*templateLen)+1)+templateLen-1;
        
        templates1 = wfMAT(cell).wf1(start:stop);
        
        for other = 1:length(wfMAT)
            if other ~= cell && wfMAT(cell).day == wfMAT(other).day && any(ismember(wfMAT(cell).whichChans,wfMAT(other).whichChans))
                templates2 = wfMAT(other).wf2(start:stop);
                
                [r,~] = corr(templates1',templates2');
                shuffcorr = [shuffcorr;r];
            end
        end
    end
    
    %% get correlations and distances of cells
    
    num = 1;
    empiricalCorr = []; empiricalCorr_PC = [];  empD_PC = []; empD = []; randD = [];
    numPC = zeros(length(directory),1);
    for day  = 1:length(directory)
        load(FolderInfo(day).name)
        disp(FolderInfo(day).name)
        for cell = 1:length(cellData)
            
            % get correlations
            empiricalCorr = [empiricalCorr;cellData(cell).pearson];
            
            if ~isempty(cellData(cell).cellID)
                empiricalCorr_PC = [empiricalCorr_PC;cellData(cell).pearson];
                numPC(day,1) = numPC(day,1)+1;
            end
    
            % get distances 
            if ~isempty(cellData(cell).startX)
                x1 = cellData(cell).startX;
                x2 = cellData(cell).endX;
                y1 = cellData(cell).startY;
                y2 = cellData(cell).endY;
                
                dist = sqrt((x1-x2)^2 + (y1-y2)^2);
                
                empD = [empD;dist];
                
                if ~isempty(cellData(cell).cellID)
                    empD_PC = [empD_PC;dist];
                end
                
                           
                for other = 1:length(cellData)
                    if cell~= other && any(ismember(cellData(cell).whichChans,cellData(other).whichChans))
                        x2 = cellData(other).endX;
                        y2 = cellData(other).endY;
                        dist = sqrt((x1-x2)^2 + (y1-y2)^2);
                        randD = [randD;dist];
                    end
                end
                
            end
        end
        clear cellData
    end
    
    
    distCutoff = prctile(randD,1);
    corrCutoff = prctile(shuffcorr,99);
end
