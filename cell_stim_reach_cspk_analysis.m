function [results] = cell_stim_reach_cspk_analysis(path_loc)


    cd(path_loc)
    directory = dir('*.mat');    
    add = 1;
    
    for day = 1:length(directory)
        
        load(directory(day).name)
        
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
        
        if isfield(cellData,'CS_Bin1')
    
            for cell = 1:length(cellData)
                
                if ~isempty(cellData(cell).CS_Bin1)
    
                    
                    preOCS = []; 
                    stimOCS = []; 
                    postOCS = []; 
                    
                    outFR=[]; allFR = [];
    
                    for reach = 1:length(ReachS)
                       
                        if ReachS(reach).stim == 1
                            [~,idx] = min(abs(cellData(cell).CS_Bin1(:,1)-ReachS(reach).stimtime));
                        else
                            [~,idx] = min(abs(cellData(cell).CS_Bin1(:,1)-ReachS(reach).real_kin(602,1)));
                        end
    
    
                        % all SS and trial SS
                        [~,idxss] = min(abs(cellData(cell).Bin10smooth_G(:,1)-ReachS(reach).real_kin(602,1)));
                        allFR = [allFR,cellData(cell).Bin10smooth_G(idxss-200:idxss+200,2)];
                        trialFRSS = cellData(cell).Bin10smooth_G(idxss-200:idxss+200,2);
    
                        % CS rate during outreach
                        [~,idxss] = min(abs(cellData(cell).CS_Bin1(:,1)-ReachS(reach).real_kin(602,1)));
                        outFRCS = sum(cellData(cell).CS_Bin1(idxss:idxss+250,2)); % SS window + 250
                            
                        % collect all reaches from block
                        if ReachS(reach).stim == 0
                            preOCS = [preOCS,outFRCS];
                        elseif ReachS(reach).stim == 1
                            stimOCS = [stimOCS,outFRCS];
                        elseif ReachS(reach).stim == 2
                            postOCS = [postOCS,outFRCS];
                        end
                    end
                   
                    standy = nanstd(cellData(cell).Bin10smooth_G(:,2)-cellData(cell).fr);
                    if max(abs(nanmean(allFR(150:250,:),2)-cellData(cell).fr)) > standy
    
                        disp(directory(day).name)
    
                        results.preOCS(:,add) = mean(preOCS);
                        results.stimOCS(:,add) = mean(stimOCS);
                        results.postOCS(:,add) = mean(postOCS);
    
                        add = add+1;
    
                    end
                end
            end
        end
    end
    



end