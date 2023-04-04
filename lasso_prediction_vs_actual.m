
function lasso_prediction_vs_actual(cell,lasso_data)

    %% Get FR prediction vs actual
    load(lasso_data(cell).directory(cell).name);

    % get outward component
    ReachS = reach_out(ReachS);

    qm_thresh = 0.75;
    framerate = 120;

    [ReachS] = qm_exclude(ReachS,qm_thresh,framerate);

    % get kinmatic data binned at 10 ms
    filter_data = 0;
    [ReachS] = reach_10ms_kinematics(ReachS,filter_data);
    
    a=1;
    reachp = []; pred = []; actu = []; reachv = [];
    % get reaches average
    for ii = 1:length(ReachS)
        if ReachS(ii).exclude == 0 
            idx1 = ((a-1)*201)+1;
            idx2 = idx1+200;

            reachv = [reachv;ReachS(ii).kin_10ms(401-100:401+100,6)'];
            reachp = [reachp;ReachS(ii).kin_10ms(401-100:401+100,2)'];
            pred = [pred;lasso_data(cell).realdata_all.MSEmin.predictedFR(idx1:idx2)'];
            actu = [actu;lasso_data(cell).realdata_all.MSEmin.actualFR(idx1:idx2)'];

            a=a+1;
        end
    end
    figure 

    hold on
    plot(nanmean(reachp),'k')
    plot(nanmean(reachp)+(nanstd(reachp)/sqrt(size(reachp,1))),'--k')
    plot(nanmean(reachp)-(nanstd(reachp)/sqrt(size(reachp,1))),'--k')
    
    hold on
    plot(nanmean(reachv),'k')
    plot(nanmean(reachv)+(nanstd(reachv)/sqrt(size(reachv,1))),'--k')
    plot(nanmean(reachv)-(nanstd(reachv)/sqrt(size(reachv,1))),'--k')

    plot(nanmean(pred),'b')
    plot(nanmean(pred)+(nanstd(pred)/sqrt(size(pred,1))),'--b')
    plot(nanmean(pred)-(nanstd(pred)/sqrt(size(pred,1))),'--b')

    plot(nanmean(actu),'r')
    plot(nanmean(actu)+(nanstd(actu)/sqrt(size(actu,1))),'--r')
    plot(nanmean(actu)-(nanstd(actu)/sqrt(size(actu,1))),'--r')
    
    title(cell)

end

