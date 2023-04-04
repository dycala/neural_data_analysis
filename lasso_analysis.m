
function lasso_results = lasso_analysis(lasso_data)


    %% format data
    if isfield(lasso_data,'realdata_all')

        % get empirical data values
        for i = 1: length(lasso_data)
            % get empirical data
            lasso_results.empirical.rsq(i)= lasso_data(i).realdata_all.MSEmin.rsq;
            lasso_results.empirical.lagmax(i)= lasso_data(i).realdata_all.MSEmin.lagmax;
            lasso_results.empirical.predictors(:,i) = lasso_data(i).realdata_all.MSEmin.predictors;
        end
    end
    if isfield(lasso_data,'spikeshuff')
        % get spike shuff rsq avgs
        for i = 1: length(lasso_data)
            if ~isempty(lasso_data(i).spikeshuff)
            for ii = 1:100
                temp(ii) = lasso_data(i).spikeshuff.MSEmin(ii).rsq;
            end
            lasso_results.spike_shuff.ss_rsq(i) = mean(temp);
            clear temp
            end
        end
    end
    if isfield(lasso_data,'reachshuff')

        % get reach shuff rsq avgs
        for i = 1:length(lasso_data)
            if isempty(lasso_data(i).reachshuff) == 0
                for ii = 1:100
                    trial_rsq(ii) = lasso_data(i).reachshuff.MinMSE(ii).rsq;
                end
                lasso_results.reach_shuff.rs_rsq(i) = mean(trial_rsq);
                clear trial_rq
            end
        end
    end
    if isfield(lasso_data,'realdata_regshuff')
        for i = 1: length(lasso_data)
            for ii = 1:length(lasso_data(i).realdata_regshuff.MSEmin)
               if ~isempty(lasso_data(i).realdata_regshuff.MSEmin(ii).rsq)
                    regshuff_rsq(ii,i) = lasso_data(i).realdata_regshuff.MSEmin(ii).rsq;
               else
                   regshuff_rsq(ii,i) = lasso_results.empirical.rsq(i);
               end
            end
            % if values at end were not in lasso they were not included so fill
            % them in
            len = length(lasso_data(i).realdata_regshuff.MSEmin);
            if len<23
               num = 23-len ;
               for ii = 1:num
                   regshuff_rsq(len + ii,i) = lasso_results.empirical.rsq(i); 
               end
            end
        end
        
        % get fraction change r2 for shuffled regressors
        for i = 1:length(lasso_results.empirical.rsq)
            lasso_results.regressor_shuff.deltareg(:,i) = (lasso_results.empirical.rsq(i)-regshuff_rsq(:,i))/ lasso_results.empirical.rsq(i);
        end
        
        % which predictors were selected
        for i = 1:length(lasso_data)
            whichselected(:,i) = lasso_data(i).realdata_all.MSE1SE.predictors~=0;
        end
        lasso_results.regressor_shuff.whichselected = int8(whichselected);
    end


    for cell =1:length(lasso_data)
        cd(lasso_data(1).path);
        load(lasso_data(1).directory(cell).name);
    
        % get outward component
        ReachS = reach_out(ReachS);
    
        qm_thresh = 0.75;
        frame_rate = 120;
        [ReachS] = qm_exclude(ReachS,qm_thresh,frame_rate);
    
        % get kinmatic data binned at 10 ms
        filter_data = 0;
        [ReachS] = reach_10ms_kinematics(ReachS,filter_data);
        
        %% Get trial average rsq
        a=1;
        pred = []; 
        actu = [];
        % get reaches average
        for ii = 1:length(ReachS)
            if ReachS(ii).exclude == 0 
                idx1 = ((a-1)*201)+1;
                idx2 = idx1+200;
    
                pred = [pred;lasso_data(cell).realdata_all.MSEmin.predictedFR(idx1:idx2)'];
                actu = [actu;lasso_data(cell).realdata_all.MSEmin.actualFR(idx1:idx2)'];
    
                a=a+1;
            end
        end
        
        yhat = nanmean(pred);
        tResponses = nanmean(actu);
        
        ybar = mean(tResponses);
        SSresid = sum((tResponses-yhat).^2);
        SStotal = sum((tResponses-ybar).^2);
        rsq = 1 - (SSresid/SStotal);
        lasso_results.t_avg_rsq(cell) = rsq;
    
        %% bin error by position 
        Spikes = get_raster(Spikes);
        Bin10smooth = i_rate_smooth(Spikes);
    
        for kin = 2:4 %x,y,z position
                a=1;
        posBinDiff = [];
            % get reaches average
            for ii = 1:length(ReachS)
                if ReachS(ii).exclude == 0 
                    idx1 = ((a-1)*201)+1;
                    idx2 = idx1+200;
                    
        
                    reach = ReachS(ii).kin_10ms(401-100:401+100,kin);
                    diff = abs(lasso_data(cell).realdata_all.MSEmin.predictedFR(idx1:idx2)-lasso_data(cell).realdata_all.MSEmin.actualFR(idx1:idx2));
                    
                    % bin error by position
                    posBins = -2.05:0.1:2.05;
                    for i = 1:length(posBins)-1
                        posBinDiff(i,a) = nanmean(diff(reach>=posBins(i) & reach<posBins(i+1)));
                    end
                    
                    a=a+1;
                end
            end
            
            % get cell avg
            if kin ==2 
                totBinDiff_x(:,cell) = nanmean(posBinDiff,2);
            elseif kin == 3
                totBinDiff_y(:,cell) = nanmean(posBinDiff,2);
            elseif kin == 4
                totBinDiff_z(:,cell) = nanmean(posBinDiff,2);
            end    
    
        end
        
    end
    
    out_x(:,1) = nanmean(totBinDiff_x,2)
    out_x(:,2) = nanstd(totBinDiff_x')'
    out_y(:,1) = nanmean(totBinDiff_y,2)
    out_y(:,2) = nanstd(totBinDiff_y')'
    out_z(:,1) = nanmean(totBinDiff_z,2)
    out_z(:,2) = nanstd(totBinDiff_z')'
    for i = 1:size(totBinDiff_x)
        out_x(i,3) = length(totBinDiff_x(i,:)) - sum(isnan(totBinDiff_x(i,:)));
        out_y(i,3) = length(totBinDiff_y(i,:)) - sum(isnan(totBinDiff_y(i,:)));
        out_z(i,3) = length(totBinDiff_z(i,:)) - sum(isnan(totBinDiff_z(i,:)));
    end 
    positions = -2:0.1:2;
    lasso_results.position_bin_error.p = [positions'];
    lasso_results.position_bin_error.x= out_x;
    lasso_results.position_bin_error.y = out_y;1
    lasso_results.position_bin_error.z = out_z;
    
    figure 
    title('model error by position')
    hold on
    scatter(lasso_results.position_bin_error.p,lasso_results.position_bin_error.x(:,1))
    scatter(lasso_results.position_bin_error.p,lasso_results.position_bin_error.y(:,1))
    scatter(lasso_results.position_bin_error.p,lasso_results.position_bin_error.z(:,1))

end
