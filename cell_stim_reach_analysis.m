function [results] = cell_stim_reach_analysis(path_loc,type,plot_results)


    cd(path_loc)
    directory = dir('*results.mat'); 
    adda = 1;
    add = 1;
    add2 = 1;
    
    for day = 1:length(directory)
        
        % preprocess data
        load(directory(day).name,'ReachS','cellData')
        disp(directory(day).name)
        
        ReachS = reach_out(ReachS);
    
        qm_thresh = 0.75;
        [ReachS] = qm_exclude(ReachS,qm_thresh,120);
    
        % get kinmatic data binned at 10 ms
        filter_data = 0;
        [ReachS] = reach_10ms_kinematics(ReachS,filter_data);
        
     
        toDel = false(length(ReachS),1);
        for reach = 1:length(ReachS)
           if ReachS(reach).exclude == 1  || isempty(ReachS(reach).stim)
               toDel(reach) = 1;
           end
        end
        ReachS = ReachS(toDel == 0);
        
        for cell = 1:length(cellData)
        
            % don't consider non cspk cells if type is CS only
            if type == 'CS only'
                if isfield(cellData,'CS_Bin1')
                   if isempty(cellData(cell).CS_Bin1) 
                       continue
                   end
                else
                    continue
                end
            end

               allFR = []; 
               pre = []; 
               stim1 = []; 
               stimearly = []; 
               stimmid = [];
               stimend = []; 
               post1 = []; 
               postearly = []; 
               postmid = []; 
               postend = []; 
    
               preO = []; 
               stim1O = [];
               stimearlyO = [];
               stimmidO = []; 
               stimendO = [];
               post1O = [];
               postearlyO = []; 
               postmidO =[]; 
               postendO = []; 
    
               preall = []; 
               stimall = []; 
               postall = [];
    
               preallO = []; 
               stimallO = []; 
               postallO = []; 
               
               
               reach_vels = [];
    
               for reach = 1:length(ReachS)
    
                 
                    % find time of stim or midpoint to align to
                    if ReachS(reach).stim == 1
                        [~,idx] = min(abs(cellData(cell).Bin10smooth_G(:,1)-ReachS(reach).stimtime));
                    else
                        [~,idx] = min(abs(cellData(cell).Bin10smooth_G(:,1)-ReachS(reach).real_kin(602,1)));
                    end
                    
                    % rate across trial
                    trialFR = cellData(cell).Bin10smooth_G(idx-200:idx+200,2)-cellData(cell).fr;
                    allFR = [allFR,trialFR];
                    reach_vels = [reach_vels,ReachS(reach).real_kin(602-200:602+200,6)];
                    
                    % rate in stim window
                    outFR = (mean(cellData(cell).Bin10(idx:idx+5,2))*100)-cellData(cell).fr;  
    
                   
                    % length of blocks and midpoint
                    st = [ReachS(:).stim] == 1;
                    st1 = find(st, 1, 'first');
                    st2 = find(st, 1, 'last');
                    idxs = floor((st1+st2)/2);
                    
                    po = [ReachS(:).stim] == 2;
                    po1 = find(po, 1, 'first');
                    po2 = find(po, 1, 'last');
                    idxp = floor((po1+po2)/2);
                   
                    if ReachS(reach).stim == 0 && ReachS(reach+5).stim == 1
                        pre = [pre,trialFR];
                        preO = [preO,outFR];
                    end
                     
                    % stim 1
                    if ReachS(reach).stim == 1 && ReachS(reach-1).stim == 0
                        stim1 = [stim1,trialFR];
                        stim1O = [stim1O,outFR];
                    end
    
                    % first 5 stim 
                    if ReachS(reach).stim == 1 && ReachS(reach-5).stim == 0
                        stimearly = [stimearly,trialFR];
                        stimearlyO = [stimearlyO,outFR];
                    end
                    
                    % middle 5 stim 
                    if ReachS(reach).stim == 1 && reach>=idxs-2 && reach<=idxs+2
                        stimmid = [stimmid,trialFR];
                        stimmidO = [stimmidO,outFR];
                    end
    
                    % last 5 stim 
                    if ReachS(reach).stim == 1 && ReachS(reach+5).stim == 2
                        stimend = [stimend,trialFR];
                        stimendO = [stimendO,outFR];
                    end
    
                    %  wash 1
                    if ReachS(reach).stim == 2 && ReachS(reach-1).stim == 1
                        post1 = [post1,trialFR];
                        post1O = [post1O,outFR];
                    end
    
                    % first 5 wash
                    if ReachS(reach).stim == 2 && ReachS(reach-5).stim == 1
                        postearly = [postearly,trialFR];
                        postearlyO = [postearlyO,outFR];
                    end
                        
                    % middle 5 wash                
                    if ReachS(reach).stim == 2 && reach>=idxp-2 && reach<=idxp+2
                        postmid = [postmid,trialFR];
                        postmidO = [postmidO,outFR];
                    end
                        
                    % last 5 wash
                    if ReachS(reach).stim == 2 && reach+5>length(ReachS)
                        postend = [postend,trialFR];
                        postendO = [postendO,outFR];
                    end
                    
                    % all reaches from block
                    if ReachS(reach).stim == 0
                        preall = [preall,trialFR];
                        preallO = [preallO,outFR];
                    elseif ReachS(reach).stim == 1 
                        stimall = [stimall,trialFR];
                        stimallO = [stimallO,outFR];
                    elseif ReachS(reach).stim == 2
                        postall = [postall,trialFR];
                        postallO = [postallO,outFR];
                    end
                    
               end
               
                  
               % analyze reach modulated cells cells (if modified 1 SD with
               % reach)
               standy = nanstd(cellData(cell).Bin10smooth_G(:,2)-cellData(cell).fr);
               if max(abs(nanmean(allFR(150:250,:),2)))> standy
    
                    results.all.reachesvel(:,adda) = nanmean(reach_vels,2);
                   
                    results.all.pre(:,adda) = nanmean(pre,2);
                    results.all.stim1(:,adda) = nanmean(stim1,2);
                    results.all.stimearly(:,adda) = nanmean(stimearly,2);
                    results.all.stimmid(:,adda) = nanmean(stimmid,2);
                    results.all.stimend(:,adda) = nanmean(stimend,2);
                    results.all.post1(:,adda) = nanmean(post1,2);
                    results.all.postearly(:,adda) = nanmean(postearly,2);
                    results.all.postmid(:,adda) = nanmean(postmid,2);
                    results.all.postend(:,adda) = nanmean(postend,2);
    
                   % get effect size of stim
                   [h,p] = ttest2(preallO,stimallO);% was all of both

                   results.frchange(adda,1) = mean(stimallO)-mean(preallO);
                   results.frchange(adda,2) = p;
                   results.frchange(adda,3) = day;
                   results.frchange(adda,4) = cell;
    
                   adda = adda+1;
    
    
                   % stim up cells
                   if h == 1 && mean(preallO(end-4:end)) < mean(stimallO(1))
                       
                        results.stim.increase.list(add,1) = day;
                        results.stim.increase.list(add,2) = cell;
                        results.stim.increase.reachvel(:,add) = nanmean(reach_vels,2);
    
                        % add FR across reaches  
                        results.stim.increase.preall(:,add) = nanmean(preall,2);
                        results.stim.increase.stimall(:,add) = nanmean(stimall,2);
    
    
                        results.stim.increase.pre(:,add) = nanmean(pre,2);
                        results.stim.increase.stim1(:,add) = nanmean(stim1,2);
                        results.stim.increase.stimearly(:,add) = nanmean(stimearly,2);
                        results.stim.increase.stimmid(:,add) = nanmean(stimmid,2);
                        results.stim.increase.stimend(:,add) = nanmean(stimend,2);
                        results.stim.increase.post1(:,add) = nanmean(post1,2);
                        results.stim.increase.postearly(:,add) = nanmean(postearly,2);
                        results.stim.increase.postmid(:,add) = nanmean(postmid,2);
                        results.stim.increase.postend(:,add) = nanmean(postend,2);
    
                        % add FR during stim window
                        results.stim.increase.preO(add) = mean(preO)-mean(preO);
                        results.stim.increase.stim1O(add) = mean(stim1O)-mean(preO);
                        results.stim.increase.stimearlyO(add) = mean(stimearlyO,2)-mean(preO);
                        results.stim.increase.stimmidO(add) = mean(stimmidO)-mean(preO);
                        results.stim.increase.stimendO(add) = mean(stimendO)-mean(preO);
                        results.stim.increase.post1O(add) = mean(post1O)-mean(preO);
                        results.stim.increase.postearlyO(add) = mean(postearlyO)-mean(preO);
                        results.stim.increase.postmidO(add) = mean(postmidO)-mean(preO);
                        results.stim.increase.postendO(add) = mean(postendO)-mean(preO);
    
    
                        % get reach velocities
                        results.stim.increase.reach_vels(:,add) = nanmean(reach_vels,2) ;
    
                        add = add+1;
    
                   % stim down cells     
                   elseif h == 1 && mean(preallO(end-4:end)) > mean(stimallO(1))
    
                        results.stim.decrease.list(add2,1) = day;
                        results.stim.decrease.list(add2,2) = cell;
                        results.stim.decrease.reachvel(:,add2) = nanmean(reach_vels,2);
    
                        % add FR across reaches
                        results.stim.decrease.preall(:,add2) = nanmean(preall,2);
                        results.stim.decrease.stimall(:,add2) = nanmean(stimall,2);
    
                        results.stim.decrease.pre(:,add2) = nanmean(pre,2);
                        results.stim.decrease.stim1(:,add2) = nanmean(stim1,2);
                        results.stim.decrease.stimearly(:,add2) = nanmean(stimearly,2);
                        results.stim.decrease.stimmid(:,add2) = nanmean(stimmid,2);
                        results.stim.decrease.stimend(:,add2) = nanmean(stimend,2);
                        results.stim.decrease.post1(:,add2) = nanmean(post1,2);
                        results.stim.decrease.postearly(:,add2) = nanmean(postearly,2);
                        results.stim.decrease.postmid(:,add2) = nanmean(postmid,2);
                        results.stim.decrease.postend(:,add2) = nanmean(postend,2);
    
                        % add FR during stim window
                        results.stim.decrease.preO(add2) = mean(preO)-mean(preO);
                        results.stim.decrease.stim1O(add2) = mean(stim1O)-mean(preO);
                        results.stim.decrease.stimearlyO(add2) = mean(stimearlyO,2)-mean(preO);
                        results.stim.decrease.stimmidO(add2) = mean(stimmidO)-mean(preO);
                        results.stim.decrease.stimendO(add2) = mean(stimendO)-mean(preO);
                        results.stim.decrease.post1O(add2) = mean(post1O)-mean(preO);
                        results.stim.decrease.postearlyO(add2) = mean(postearlyO)-mean(preO);
                        results.stim.decrease.postmidO(add2) = mean(postmidO)-mean(preO);
                        results.stim.decrease.postendO(add2) = mean(postendO)-mean(preO);
                        
                        add2 = add2+1;
                   end
   
             
            end

        end
    
    end



    %% plot

    if plot_results == true
    
        results.comb_stim1 = [-(results.stim.decrease.stim1-results.stim.decrease.pre),results.stim.increase.stim1-results.stim.increase.pre]
        results.comb_stimearly = [-(results.stim.decrease.stimearly-results.stim.decrease.pre),results.stim.increase.stimearly-results.stim.increase.pre]
        results.comb_stimmid = [-(results.stim.decrease.stimmid-results.stim.decrease.pre),results.stim.increase.stimmid-results.stim.increase.pre]
        results.comb_stimlate = [-(results.stim.decrease.stimend-results.stim.decrease.pre),results.stim.increase.stimend-results.stim.increase.pre]
        
        results.comb_post1 = [-(results.stim.decrease.post1-results.stim.decrease.pre),(results.stim.increase.post1-results.stim.increase.pre)]
        results.comb_postearly = [-(results.stim.decrease.postearly-results.stim.decrease.pre),(results.stim.increase.postearly-results.stim.increase.pre)]
        results.comb_postmid = [-(results.stim.decrease.postmid-results.stim.decrease.pre),(results.stim.increase.postmid-results.stim.increase.pre)]
        results.comb_postlate = [-(results.stim.decrease.postend-results.stim.decrease.pre),(results.stim.increase.postend-results.stim.increase.pre)]
        
        net_stim_increasers = [results.stim.increase.preO',...
            results.stim.increase.stim1O',...
            results.stim.increase.stimearlyO',...
            results.stim.increase.stimmidO',...
            results.stim.increase.stimendO',...
            results.stim.increase.post1O',...
            results.stim.increase.postearlyO',...
            results.stim.increase.postmidO',...
            results.stim.increase.postendO'];
        
        net_stim_decreasers = [results.stim.decrease.preO',...
            results.stim.decrease.stim1O',...
            results.stim.decrease.stimearlyO',...
            results.stim.decrease.stimmidO',...
            results.stim.decrease.stimendO',...
            results.stim.decrease.post1O',...
            results.stim.decrease.postearlyO',...
            results.stim.decrease.postmidO',...
            results.stim.decrease.postendO'];
        
        results.net_stim_all = [net_stim_increasers;-(net_stim_decreasers)];



        figure
        plot(nanmean(results.net_stim_all,1))
        
        figure
        subplot(1,2,1)
        hold on
        plot(nanmean(results.comb_stim1,2))
        plot(nanmean(results.comb_stimearly,2))
        plot(nanmean(results.comb_stimmid,2))
        plot(nanmean(results.comb_stimlate,2))
        subplot(1,2,2)
        hold on
        plot(nanmean(results.comb_post1,2))
        plot(nanmean(results.comb_postearly,2))
        plot(nanmean(results.comb_postmid,2))
        plot(nanmean(results.comb_postlate,2))
        
        % stim increasers
        figure
        subplot(1,2,1)
        hold on
        plot(nanmean(results.stim.increase.pre,2))
        plot(nanmean(results.stim.increase.stim1,2)) 
        plot(nanmean(results.stim.increase.stimearly,2))
        plot(nanmean(results.stim.increase.stimmid,2))
        plot(nanmean(results.stim.increase.stimend,2))
        
        % stim decreasers
        subplot(1,2,2)
        hold on
        plot(nanmean(results.stim.decrease.pre,2))
        plot(nanmean(results.stim.decrease.stim1,2))
        plot(nanmean(results.stim.decrease.stimearly,2))
        plot(nanmean(results.stim.decrease.stimmid,2))
        plot(nanmean(results.stim.decrease.stimend,2))
        
        figure 
        subplot(2,1,1)
        ts = -200:200;
        hold on
        meanSEMplot(ts,results.all.stim1-results.all.pre,'c')
        meanSEMplot(ts,results.all.stimearly-results.all.pre,'c')
        meanSEMplot(ts,results.all.stimmid-results.all.pre,'g')
        meanSEMplot(ts,results.all.stimend-results.all.pre,'b')
        title('Stim block')
        subplot(2,1,2)
        hold on
        meanSEMplot(ts,results.all.post1-results.all.pre,'r')
        meanSEMplot(ts,results.all.postearly-results.all.pre,'r')
        meanSEMplot(ts,results.all.postmid-results.all.pre,'y')
        meanSEMplot(ts,results.all.postend-results.all.pre,'m')
        title('Washout block')
    end

end

