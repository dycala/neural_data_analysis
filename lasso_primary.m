
function [lasso_data] = lasso_primary(path,save_loc,program,lasso_data)
    %% Setup
    %program = 'empirical'; %'empirical' 'spike shuffle' 'reach shuffle' 'regressor shuffle'
    filter_data = 0;
    
    cd(path)
    directory = dir('*Compiled.mat'); 
    lasso_data(1).path = path;
    lasso_data(1).directory = directory;
    
    %% load data and Run LASSO for each cell
    loops_complete = 0;
    for i = 1:2%length(directory)
        cd(path);
        load(directory(i).name);
    
        % get outward component
        ReachS = reach_out(ReachS);
    
        qm_thresh = 0.75;
        frame_rate = 120;
        [ReachS] = qm_exclude(ReachS,qm_thresh,frame_rate);
    
        % get kinmatic data binned at 10 ms
        [ReachS] = reach_10ms_kinematics(ReachS,filter_data);
    
        % get binned and smoothed spikes for each reach
        Spikes = get_raster(Spikes);
        Bin10smooth = i_rate_smooth(Spikes);
    
        switch program
            case 'empirical'
                
                % empirical regression
                [realdata_all] = lasso_empirical(ReachS,Bin10smooth);
                lasso_data(i).name = directory(i).name(1:end-13);
                lasso_data(i).realdata_all = realdata_all;
                
            case 'spike shuffle'
    
                % get shuffled spikes (assumes a 200 point regression)
                raster = Spikes.SS_Raster;
                lagmax = lasso_data(i).realdata_all.MSEmin.lagmax;
                [shuff] = shuffle_spikes(ReachS,raster,lagmax);
    
                % LASSO spike shuffle
                [spikeshuffdata] = lasso_spike_shuffle(ReachS,shuff,lagmax);
                lasso_data(i).spikeshuff = spikeshuffdata;
    
            case 'reach shuffle'    
                
                % LASSO reach shuffle
                lagmax = lasso_data(i).realdata_all.MSEmin.lagmax;
                [reachshuffdata] = lasso_reach_shuffle(ReachS,Bin10smooth,lagmax);
                lasso_data(i).reachshuff = reachshuffdata;
    
            case 'regressor shuffle'  
                
                % shuffle by individual variable
                lagmax = lasso_data(i).realdata_all.MSEmin.lagmax;
                empirical_predictors = lasso_data(i).realdata_all.MSEmin.predictors;
                [realdata_regshuff] = lasso_regressor_shuffle(ReachS,Bin10smooth,lagmax,empirical_predictors);
                lasso_data(i).realdata_regshuff = realdata_regshuff;
        end
    
        loops_complete = loops_complete+1
    
    end

end