
function [] = cell_reach_analysis(path,save_loc)

    cd(path)
    
    directory = dir('*Compiled.mat');
    numreaches = 0
    
    % load individual kinematic files
    for day = 1:length(directory)
        
       %cd(pathloc)
       load(directory(day).name)
       disp(directory(day).name)
        
       ReachS = reach_out(ReachS);
    
       qm_thresh = 0.75;
       [ReachS] = qm_exclude(ReachS,qm_thresh,120);
       
       filter_data = 0;
       [ReachS] = reach_10ms_kinematics(ReachS,filter_data);
    
       Spikes = get_raster(Spikes);
       Bin10smooth = i_rate_smooth(Spikes);
       
        to_del = false(length(ReachS),1);
        for reach = 1:length(ReachS)
           if ReachS(reach).exclude == 1 
               to_del(reach) = 1;
           end
        end
        ReachS = ReachS(to_del == 0);
    
       % get smoothed FR aligned to Reaches
       for ii = 1:length(ReachS)
    
               numreaches = numreaches+1;
               % find threshold crossing
               cross = round(ReachS(ii).real_kin(602,1),3);
               [~,idx] = min(abs(Bin10smooth(:,1)-cross));
               % 2 sec on either side
               FR(day).reach(:,ii) = Bin10smooth(idx-200:idx+200,2);
               allReaches(day).reachxpos(:,ii) = ReachS(ii).kin_10ms(401-200:401+200,2);
               allReaches(day).reachspeed(:,ii) = ReachS(ii).kin_10ms(401-200:401+200,5);
               allReaches(day).reachxvel(:,ii) = ReachS(ii).kin_10ms(401-200:401+200,6);
               allReaches(day).reachyvel(:,ii) = ReachS(ii).kin_10ms(401-200:401+200,7);
               allReaches(day).reachzvel(:,ii) = ReachS(ii).kin_10ms(401-200:401+200,8);
               allReaches(day).reachxdecel(:,ii) = ReachS(ii).kin_10ms(401-200:401+200,10);
        end
    
       FR(day).average = mean((Spikes.SS_Bin10(:,2)*100)); % fr in Hz
       
    end
        
    %% get normalized FR of each cell sorted by time of modulation
    
    % determine whether cell is burster or pauser and time of max deviation
    data.frs = nan(length(FR),403);
    for i = 1:length(FR)
        
        data.frs(i,1:401) = nanmean(FR(i).reach,2)-FR(i).average;
        
        % find if cell increased or decreased FR in 1000 ms around threshold
        up = mean(data.frs(i,150:250));
        
        if up>0
           [~,peak(i)] = max(data.frs(i,150:250));
           nadir(i) = 0;
           data.frs(i,402) = peak(i);
        else 
           [~,nadir(i)] = min(data.frs(i,150:250));
           peak(i) = 0;
           data.frs(i,403) = nadir(i);
        end 
        
    end
    peaktimes = [peak(peak~=0)';nadir(nadir~=0)'];
    
    % make list of bursters and pausers
    for cell = 1:length(peak)
        if peak(cell) ~= 0
            data.bp_list{cell,1} = 'Burst';
        elseif nadir(cell) ~= 0
            data.bp_list{cell,1} = 'Pause';
        end
    end
    
    % first sort bursters
    data.frs = sortrows(data.frs,402);
    
    % find which are pausers and sort
    idx = sum(isnan(data.frs(:,403)))+1;
    data.frs(idx:end,:) = sortrows(data.frs(idx:end,:),403,'descend');
    
    data.burst_all = data.frs(1:idx-1,:);
    data.pause_all = data.frs(idx:end,:);
    
    data.burst_mean = nanmean(data.burst_all(:,1:401),1);
    data.pause_mean = nanmean(data.pause_all(:,1:401),1);
    
    % plot results
    hold on
    plot(data.burst_mean)
    plot(data.pause_mean)
    [val,idx]=max(data.burst_mean)
    [val,idx]=min(data.pause_mean)
    
    %% Firing rate heat map
    % normalize
    for i = 1:size(data.frs,1)
        trialMax = abs(max(data.frs(i,1:401)));
        trialMin = abs(min(data.frs(i,1:401)));
        if trialMax>=trialMin
            data.frs_norm(i,1:401) = data.frs(i,1:401)/trialMax;
        else 
            data.frs_norm(i,1:401) = data.frs(i,1:401)/abs(trialMin);
        end
    end
    
    % plot heat map
    [cmm]=cbrewer2('div', 'RdBu', 100, 'cubic');
    figure
    for num = 1:size(data.frs_norm,1)
    
        % flip
        row = size(data.frs_norm,1)+1-num;
        x = [-2000:10:2000];
        y = (ones(size(x))*row);
        z = zeros(size(x));
    
        col = data.frs_norm(num,:);  % This is the color, vary with x in this case.
        surface([x;x],[y;y],[z;z],[col;col],...
                'facecol','no',...
                'edgecol','interp',...
                'linew',1);
        colormap(flipud(cmm))
        
    end
    xlim([-2100 2100])
    ylim([-1 size(data.frs_norm,1)])
    colorbar()
    xlabel('Time from threshold crossing (ms)')
    ylabel('Cell')
    
    %% get reach averages
    for i = 1:length(allReaches)  
        data.grouped_reaches.speed(:,i) = nanmean(allReaches(i).reachspeed,2);
        data.grouped_reaches.xvel(:,i) = nanmean(allReaches(i).reachxvel,2);
        data.grouped_reaches.yvel(:,i) = nanmean(allReaches(i).reachyvel,2);
        data.grouped_reaches.zvel(:,i) = nanmean(allReaches(i).reachzvel,2);
    end
    
    % plot reaches
    hold on
    plot(nanmean(fRs(:,1:401)))
    plot(nanmean(data.grouped_reaches.speed,2))
    plot(nanmean(data.grouped_reaches.yvel,2))
    
    
    %% compute average FR
    
    baseline_fr = []; pause_baseline_fr = []; burst_baseline_fr = [];
    for cell = 1:length(FR)
        baseline_fr = [baseline_fr;FR(cell).average];
        if data.bp_list{cell,1} == 'Pause'
            pause_baseline_fr = [pause_baseline_fr;FR(cell).average];
        elseif data.bp_list{cell,1} == 'Burst'
            burst_baseline_fr = [burst_baseline_fr;FR(cell).average];
        end
    end
    
    data.mean_fr = nanmean(baseline_fr);
    data.mean_burst_fr = nanmean(burst_baseline_fr);
    data.mean_pause_fr = nanmean(pause_baseline_fr);
    
    %% bin reaches by reach velocity percentile (NNresub)
    
    percentiles = [0:20:100];
    
    clear binnedcells binnedreaches
    for prc = 1:length(percentiles)
        binnedcells(prc).all = [];
        binnedcells(prc).burst = [];
        binnedcells(prc).pause = [];
        
        binnedreaches(prc).allxv = [];
        binnedreaches(prc).allyv = [];
        binnedreaches(prc).allxp = [];
        binnedreaches(prc).allxa = [];
       
        binnedreaches(prc).burstxv = [];
        binnedreaches(prc).pausexv = [];
    end
    
    % put trials into bins
    for i = 1:length(FR)
        disp(i)
        clear val
        
        for ii = 1:size(FR(i).reach,2)
            [val(ii),~] = max(allReaches(i).reachxvel(195:205,ii));
        end
        
        % define velocity bins
        vel_quantile = prctile(val,percentiles);
        
        % sort FRs into bins
        for ii = 1:size(FR(i).reach,2)
            
            [maxvel,~] = max(allReaches(i).reachxvel(195:205,ii));
            
            for bin = 1:length(binnedcells)
                
    
                if maxvel<=vel_quantile(bin) 
                    
                   % bin all
                   binnedcells(bin).all = [binnedcells(bin).all,FR(i).reach(:,ii)];
                   
                   binnedreaches(bin).allxv = [binnedreaches(bin).allxv,allReaches(i).reachxvel(:,ii)];
                   binnedreaches(bin).allyv = [binnedreaches(bin).allyv,allReaches(i).reachyvel(:,ii)];
                   binnedreaches(bin).allxp = [binnedreaches(bin).allxp,allReaches(i).reachxpos(:,ii)];
                   binnedreaches(bin).allxa = [binnedreaches(bin).allxa,allReaches(i).reachxdecel(:,ii)];
                  
                   if data.bp_list{i,1} == 'Burst'
                       
                      binnedcells(bin).burst = [binnedcells(bin).burst,FR(i).reach(:,ii)];
                      binnedreaches(bin).burstxv = [binnedreaches(bin).burstxv,allReaches(i).reachxvel(:,ii)];
                   
                   elseif data.bp_list{i,1} == 'Pause'
                       
                      binnedcells(bin).pause = [binnedcells(bin).pause,FR(i).reach(:,ii)];
                      binnedreaches(bin).pausexv = [binnedreaches(bin).pausexv,allReaches(i).reachxvel(:,ii)];
    
                   end
                   
                   break
                   
                end
            end
    
        end
        
    end
    
    % combine first and second bins (first was 0th percentile value)
    binnedcells(2).all = [binnedcells(2).all ,binnedcells(1).all ];
    binnedcells(2).burst = [binnedcells(2).burst ,binnedcells(1).burst ];
    binnedcells(2).pause = [binnedcells(2).pause ,binnedcells(1).pause ];
    
    binnedreaches(2).allxv = [binnedreaches(2).allxv,binnedreaches(1).allxv];
    binnedreaches(2).allyv = [binnedreaches(2).allyv,binnedreaches(1).allyv];
    binnedreaches(2).allxa = [binnedreaches(2).allxa,binnedreaches(1).allxa];
    binnedreaches(2).allxp = [binnedreaches(2).allxp,binnedreaches(1).allxp];
    binnedreaches(2).burst = [binnedreaches(2).burstxv ,binnedreaches(1).burstxv ];
    binnedreaches(2).pause = [binnedreaches(2).pausexv ,binnedreaches(1).pausexv ];
    
    % drop unfilled first quantile
    binnedcells(1) = [];
    binnedreaches(1) = [];
    
    % plot velocity and FR
    figure
    [c1]=cbrewer2('seq', 'Purples', length(binnedcells)*2+5);
    t = [-2000:10:2000];
    cop_reaches = t'; cop_frs = t';
    figure
    for i = 1:length(binnedcells)
        subplot(2,1,1)
        title('Reach velocity curves')

        hold on
        plot(nanmean(binnedreaches(i).allxv,2),'Color',c1((i*2)+5,:))
        cop_reaches = [cop_reaches,nanmean(binnedreaches(i).allxv,2)];
        subplot(2,1,2)
        hold on
        plot(nanmean(binnedcells(i).all,2)-data.mean_fr,'Color',c1((i*2)+5,:))
        cop_frs = [cop_frs,nanmean(binnedcells(i).all,2)-data.mean_fr];
        
    end

    %% save data
    for i = 1:length(binnedreaches)
        data.population_quintiles.percentiles = percentiles  
        data.population_quintiles.reach_velocity(i).mean = nanmean(binnedreaches(i).allxv,2);
        data.population_quintiles.reach_velocity(i).sd = nanstd(binnedreaches(i).allxv');
        data.population_quintiles.reach_velocity(i).num = size(binnedreaches(i).allxv,2);
    
        data.population_quintiles.firing_rate(i).mean = nanmean(binnedcells(i).all,2);
        data.population_quintiles.firing_rate(i).sd = nanstd(binnedcells(i).all')';
        data.population_quintiles.firing_rate(i).num = size(binnedcells(i).all,2);
    end

    %% quantify fr suppression
    clear quantkeep
    % new quant
    clear drop_mean halfFRvalIdx reachStartIdx maxAccIdx maxVelIdx minAccIdx
    for bin = 1:length(binnedcells)
        
        % find max before
        indices = 175:225; 
        fr = nanmean(binnedcells(bin).all(indices,:),2);
        [val_max,idx_max] = max(fr);
        idx_max = idx_max+174;% adjust index to to match searched range
    
        % find min after
        indices = 175:225; 
        fr = nanmean(binnedcells(bin).all(indices,:),2);    
        [val_min,idx_min] = min(fr);
        idx_min = idx_min+174; 
        
        drop_mean(bin,1) = mean(binnedcells(bin).all(idx_min,:));
        drop_mean(bin,2) = std(binnedcells(bin).all(idx_min,:));
        drop_mean(bin,3) = numel(binnedcells(bin).all(idx_min,:));
    
        relative_drop_mean(bin,1) = mean(binnedcells(bin).all(idx_min,:))-data.mean_fr;
        relative_drop_mean(bin,2) = std(binnedcells(bin).all(idx_min,:));
        relative_drop_mean(bin,3) = numel(binnedcells(bin).all(idx_min,:));
    
        quantkeep(bin).which(:,1) = [(binnedcells(bin).all(idx_min,:))-data.mean_fr]';
       
        % find half drop time
        peak2trough_half = (val_max + val_min)/2;
        [~,half_idx] = min(abs(fr-peak2trough_half));
            
        halfFRvalIdx(bin) = floor((idx_min+idx_max)/2)-174 %half_idx;
        
        % get reach characteristics
        indices_reach = 175:225; % new
        [~,maxVelIdx(bin)] = max(nanmean(binnedreaches(bin).allxv(indices_reach,:),2));
        [~,maxAccIdx(bin)] = max(nanmean(binnedreaches(bin).allxa(indices_reach,:),2));
        [~,minAccIdx(bin)] = min(nanmean(binnedreaches(bin).allxa(indices_reach,:),2));
        
        % find reach start
        for t = 200:-1:175
            if nanmean(binnedreaches(bin).allxv(t,:),2)<2  && nanmean(binnedreaches(bin).allyv(t,:),2)<2  && nanmean(binnedreaches(bin).allxp(t,:),2)<0.5
                reachStartIdx(bin) = t-174;
                break
            end
        end
    
    end
    
    
    % plot fr drop
    scatter(percentiles(2:end),drop_mean(:,1))
    [h,p] = corr(percentiles(2:end)',drop_mean(:,1))
    
    % plot times
    cop = [halfFRvalIdx'-halfFRvalIdx',reachStartIdx'-halfFRvalIdx',maxAccIdx'-halfFRvalIdx',maxVelIdx'-halfFRvalIdx',minAccIdx'-halfFRvalIdx']*10;
    plot(nanmean(cop))
    
    % quantify number of bursters and pausers per bin
    clear Bnum Pnum    
    
    for i = 1:5
        hold on
        plot(data.population_quintiles.firing_rate(i).mean)
        plot(data.population_quintiles.reach_velocity(i).mean)
        scatter(out(2,i)+175,out(1,i))
    end
    
    for i = 1:length(binnedcells)
        Bnum(i) = size(binnedcells(i).burst,2);
        Pnum(i) = size(binnedcells(i).pause,2);
    end
    plot(Bnum(1,:)./(Bnum(1,:)+Pnum(1,:)))
    legend('number of bursters by bin','number of pausers by bin')
    cop2 = Bnum./(Bnum+Pnum)

end
