
function CS_results = cell_reach_cspk_analysis(path_loc, CS_type, window)

    %% analysis of cspks in all cells 
    
    % CS_type options: 'after' analyzes CSs that happen after reach start; 'before' analyzes CSs that happen before reach start
    % window is length in ms of time to analyze before or after each reach for Cspks
    cd(path_loc)
    directory = dir('*Compiled.mat');
    
    CS_results.totspikes_CStrials = [];
    CS_results.totspikes_noCStrials = [];
    CS_results.totspikes_CStrialsnext = [];
    
    CS_results.totspikes_CStrials_ts = [];
    CS_results.totspikes_noCStrials_ts = [];
    CS_results.totspikes_CStrialsnext_ts = [];
    
    CS_results.totCS_aligned_SS = [];
    CS_results.totnoCS_aligned_SS = [];
    CS_results.totCS_aligned_SSnext = [];
    CS_results.totCS_aligned_SSprev = [];
    
    CS_results.totCS_aligned_SSmeanbefore = [];
    CS_results.totCS_aligned_SSprevmeanbefore = [];
    CS_results.totCS_aligned_SSnextmeanbefore = [];
    
    CS_results.mean_tot_relendpoint_CSx = [];
    CS_results.mean_tot_relendpoint_CSy = [];
    CS_results.mean_tot_relendpoint_CSz = [];
    
    CS_results.mean_tot_relendpoint_noCSx = [];
    CS_results.mean_tot_relendpoint_noCSy = [];
    CS_results.mean_tot_relendpoint_noCSz = [];
    
    for day = 1:length(directory)
        
       % preprocess
       load(directory(day).name)
       
       ReachS = reach_out(ReachS);
    
       qm_thresh = 0.75;
       framerate = 120;
       [ReachS] = qm_exclude(ReachS,qm_thresh,framerate);
       
       [ReachS] = reach_10ms_kinematics(ReachS,1);
    
       Spikes = get_raster(Spikes);
       Bin10smooth = i_rate_smooth(Spikes);
       
       spikes_noCStrials = []; spikes_CStrials = []; spikes_CStrialsnext = [];
       spikes_CStrials_ts = []; spikes_noCStrials_ts = []; spikes_CStrialsnext_ts = [];
       CS_aligned_SS = []; CS_aligned_SSnext = []; CS_aligned_SSprev = [];
       CS_aligned_SSprevmeanbefore = []; CS_aligned_SSmeanbefore = []; CS_aligned_SSnextmeanbefore = [];
       endpoint_CSx = []; endpoint_CSy = []; endpoint_CSz = [];
       endpoint_noCSx = []; endpoint_noCSy = []; endpoint_noCSz = [];
       
       for ii = 1:length(ReachS)
           if ReachS(ii).exclude == 0 
                         
               baseline_fr = mean(Bin10smooth(:,2));
               
               % get index of midpoint, start, and end of outreach.
               [~,idx] = min(abs(Spikes.CS_Bin10(:,1)-ReachS(ii).real_kin(602,1)));
               [~,idx1] = min(abs(Bin10smooth(:,1)-(ReachS(ii).out(1,1))));
               [~,idx2] = min(abs(Bin10smooth(:,1)-(ReachS(ii).out(end,1))));
               
               % find where CS occurs in out up and lat
               if CS_type == "before"
                   CS_twindow = idx1-window:idx1;
               elseif CS_type == "after"
                   CS_twindow = idx1:idx2+window;
               end
                   
               for t = CS_twindow
                   if Spikes.CS_Bin10(t,2) == 1
                       
                       % align SS to kinematics
                       [~,CS_idx_SS] = min(abs(Bin10smooth(:,1)-Spikes.CS_Bin10(t,1)));
                       CS_aligned_SS = [CS_aligned_SS,Bin10smooth(CS_idx_SS-100:CS_idx_SS+100,2)];
                       CS_aligned_SSmeanbefore = [CS_aligned_SSmeanbefore,mean(Bin10smooth(CS_idx_SS-10:CS_idx_SS-1,2))];
                       
    
                       [~,idx_mid] = min(abs(Bin10smooth(:,1)-ReachS(ii).real_kin(602,1))); 
                       th_diff_idx = t-idx_mid;
                       
                       % align next reach if not last reach in exp
                       if ii~= length(ReachS)
                           [~,nidx_mid] = min(abs(Bin10smooth(:,1)-ReachS(ii+1).real_kin(602,1)));
                           nidx = nidx_mid+th_diff_idx;
                           CS_aligned_SSnext = [CS_aligned_SSnext,Bin10smooth(nidx-100:nidx+100,2)];
                           CS_aligned_SSnextmeanbefore = [CS_aligned_SSnextmeanbefore,mean(Bin10smooth(nidx-10:nidx-1,2))];
                         
                       end
                       
                       % align previous reach
                       if ii~= 1
                           [~,pidx_mid] = min(abs(Bin10smooth(:,1)-ReachS(ii-1).real_kin(602,1)));
                           pidx = pidx_mid+th_diff_idx;
                           CS_aligned_SSprev = [CS_aligned_SSprev,Bin10smooth(pidx-100:pidx+100,2)];
                           CS_aligned_SSprevmeanbefore = [CS_aligned_SSprevmeanbefore,mean(Bin10smooth(pidx-10:pidx-1,2))];
                          
                       end
                       
                   end
               end
    
               % detect CS in window
               if CS_type =="before"
                    numCS = sum(Spikes.CS_Bin10(idx1-window:idx1,2));
               elseif CS_type == "after"
                    numCS = sum(Spikes.CS_Bin10(idx1:idx2+window,2));
               end
              
               % kinematics of CS vs non CS reaches
               if numCS > 0
                    
                endpoint_CSx = [endpoint_CSx,ReachS(ii).out(end,2)];
                endpoint_CSy = [endpoint_CSy,ReachS(ii).out(end,3)];
                endpoint_CSz = [endpoint_CSz,ReachS(ii).out(end,4)];
                
                spikes_CStrials = [spikes_CStrials;mean(Bin10smooth(idx1:idx2+window,2))];
                spikes_CStrials_ts = [spikes_CStrials_ts,Bin10smooth(idx-100:idx+100,2)];
                   
                  % align trial after CS trial for analysis 
                  if ii~= length(ReachS)
                       [~,nidx1] = min(abs(Bin10smooth(:,1)-(ReachS(ii+1).out(1,1))));
                       [~,nidx2] = min(abs(Bin10smooth(:,1)-(ReachS(ii+1).out(end,1)))); 
                       [~,nidx_mid] = min(abs(Bin10smooth(:,1)-ReachS(ii+1).real_kin(602,1))); 
    
                       spikes_CStrialsnext = [spikes_CStrialsnext,mean(Bin10smooth(nidx1:nidx2+window,2))];
                       spikes_CStrialsnext_ts = [spikes_CStrialsnext_ts,Bin10smooth(nidx_mid-100:nidx_mid+100,2)];
                  end
                
               elseif numCS == 0
     
                   endpoint_noCSx = [endpoint_noCSx,ReachS(ii).out(end,2)];
                   endpoint_noCSy = [endpoint_noCSy,ReachS(ii).out(end,3)];
                   endpoint_noCSz = [endpoint_noCSz,ReachS(ii).out(end,4)];
    
                   spikes_noCStrials = [spikes_noCStrials;mean(Bin10smooth(idx1:idx2+window,2))];
                   spikes_noCStrials_ts = [spikes_noCStrials_ts,Bin10smooth(idx-100:idx+100,2)];
    
               end              
           end
       end
    
       % compare CS and no CS trials (in sessions that have both)
       if size(endpoint_CSx,2)>0 && size(endpoint_noCSx,2)>0
    
            % append mean spikes
    
            CS_results.totspikes_CStrials = [CS_results.totspikes_CStrials;nanmean(spikes_CStrials)];
            CS_results.totspikes_noCStrials = [CS_results.totspikes_noCStrials;nanmean(spikes_noCStrials)];
            CS_results.totspikes_CStrialsnext = [CS_results.totspikes_CStrialsnext,nanmean(spikes_CStrialsnext)];
            
            CS_results.totspikes_noCStrials_ts = [CS_results.totspikes_noCStrials_ts,nanmean(spikes_noCStrials_ts,2)-baseline_fr];
            CS_results.totspikes_CStrials_ts = [CS_results.totspikes_CStrials_ts,nanmean(spikes_CStrials_ts,2)-baseline_fr];
            CS_results.totspikes_CStrialsnext_ts = [CS_results.totspikes_CStrialsnext_ts,nanmean(spikes_CStrialsnext_ts,2)-baseline_fr];
            
            CS_results.totCS_aligned_SS = [CS_results.totCS_aligned_SS,nanmean(CS_aligned_SS,2)-baseline_fr];
            CS_results.totCS_aligned_SSnext = [CS_results.totCS_aligned_SSnext,nanmean(CS_aligned_SSnext,2)-baseline_fr];
            CS_results.totCS_aligned_SSprev = [CS_results.totCS_aligned_SSprev,nanmean(CS_aligned_SSprev,2)-baseline_fr];
    
    
            CS_results.totCS_aligned_SSmeanbefore = [CS_results.totCS_aligned_SSmeanbefore,nanmean(CS_aligned_SSmeanbefore)];        
            CS_results.totCS_aligned_SSprevmeanbefore = [CS_results.totCS_aligned_SSprevmeanbefore,nanmean(CS_aligned_SSprevmeanbefore)];
            CS_results.totCS_aligned_SSnextmeanbefore = [CS_results.totCS_aligned_SSnextmeanbefore,nanmean(CS_aligned_SSnextmeanbefore)];
    
            % endpoint relative to session median
            CS_results.mean_tot_relendpoint_CSx = [CS_results.mean_tot_relendpoint_CSx;mean(endpoint_CSx'-median([endpoint_CSx,endpoint_noCSx]))]; 
            CS_results.mean_tot_relendpoint_CSy = [CS_results.mean_tot_relendpoint_CSy;mean(endpoint_CSy'-median([endpoint_CSy,endpoint_noCSy]))];  
            CS_results.mean_tot_relendpoint_CSz = [CS_results.mean_tot_relendpoint_CSz;mean(endpoint_CSz'-median([endpoint_CSz,endpoint_noCSz]))]; 
    
            CS_results.mean_tot_relendpoint_noCSx = [CS_results.mean_tot_relendpoint_noCSx;mean(endpoint_noCSx'-median([endpoint_CSx,endpoint_noCSx]))]; 
            CS_results.mean_tot_relendpoint_noCSy = [CS_results.mean_tot_relendpoint_noCSy;mean(endpoint_noCSy'-median([endpoint_CSy,endpoint_noCSy]))];  
            CS_results.mean_tot_relendpoint_noCSz = [CS_results.mean_tot_relendpoint_noCSz;mean(endpoint_noCSz'-median([endpoint_CSz,endpoint_noCSz]))]; 
            
    
       end
    
       disp(day)
       %clear Spikes
       
    end
    
    % compute euclidean distance
    clear CS noCS data
    CS = [CS_results.mean_tot_relendpoint_CSx,CS_results.mean_tot_relendpoint_CSy,CS_results.mean_tot_relendpoint_CSz];
    noCS = [CS_results.mean_tot_relendpoint_noCSx,CS_results.mean_tot_relendpoint_noCSy,CS_results.mean_tot_relendpoint_noCSz];
    
    clear euc_CS euc_noCS
    for i = 1:length(CS)
        CS_results.euc_CS(i,1) = norm(CS(i,:));
        CS_results.euc_noCS(i,1) = norm(noCS(i,:));
    end
    
    mean(CS_results.euc_noCS)
    mean(CS_results.euc_CS)
    
    [h,p] = ttest(CS_results.euc_noCS,CS_results.euc_CS)
    cop = [CS_results.euc_noCS,CS_results.euc_CS]
    
    %% Figures
    % plot timeseries CS trials vs no CS trials
    figure
    meanSEMplot(-100:100,CS_results.totspikes_CStrials_ts,'r')
    meanSEMplot(-100:100,CS_results.totspikes_noCStrials_ts,'b')
    title('Sspk on CS trials vs no CS trials')
    
    % plot timeseries CS aligned trials vs previous trial vs next trial
    figure
    meanSEMplot(-100:100,CS_results.totCS_aligned_SS,'k')
    meanSEMplot(-100:100,CS_results.totCS_aligned_SSnext,'g')
    meanSEMplot(-100:100,CS_results.totCS_aligned_SSprev,'b')
    title('Sspk aligned to CS on CS trial, previous, and next trials')
    
    % analyze kinematic difference in CS trials and no CS trials
    % plot endpoints by CS vs no CS in out, up, and lat
    data_CS = [CS_results.mean_tot_relendpoint_CSx,CS_results.mean_tot_relendpoint_CSy,CS_results.mean_tot_relendpoint_CSz];
    data_noCS = [CS_results.mean_tot_relendpoint_noCSx,CS_results.mean_tot_relendpoint_noCSy,CS_results.mean_tot_relendpoint_noCSz];
    
    figure
    hold on
    for i = 1:length(CS_results.mean_tot_relendpoint_CSx)
        plot([data_CS(i,1);data_noCS(i,1)],[data_CS(i,2);data_noCS(i,2)],'k')
    end
    scatter(CS_results.mean_tot_relendpoint_CSx,CS_results.mean_tot_relendpoint_CSy,'filled','MarkerFaceColor','r','MarkerEdgeColor','k')
    scatter(CS_results.mean_tot_relendpoint_noCSx,CS_results.mean_tot_relendpoint_noCSy,'filled','MarkerFaceColor',[0.4,0.4,0.4],'MarkerEdgeColor','k')
    xlabel('Outward')
    ylabel('Upward')
    title('No CS vs CS trial endpoints')
    
    figure
    hold on
    for i = 1:length(CS_results.mean_tot_relendpoint_CSx)
        plot([data_CS(i,1);data_noCS(i,1)],[data_CS(i,3);data_noCS(i,3)],'k')
    end
    scatter(CS_results.mean_tot_relendpoint_CSx,CS_results.mean_tot_relendpoint_CSz,'filled','MarkerFaceColor','r','MarkerEdgeColor','k')
    scatter(CS_results.mean_tot_relendpoint_noCSx,CS_results.mean_tot_relendpoint_noCSz,'filled','MarkerFaceColor',[0.4,0.4,0.4],'MarkerEdgeColor','k')
    xlabel('Outward')
    ylabel('Lateral')
    title('No CS vs CS trial endpoints')


end
