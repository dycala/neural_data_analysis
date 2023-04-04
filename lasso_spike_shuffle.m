%% Spike Shuffle Regression

function [spikeshuffdata] = lasso_spike_shuffle(ReachS,shuff,lagmax)
%% Set up predictors

% clip and concatenate reaches
ReachS(1).kin_cat_10ms = [];
for i=1:length(ReachS)
    if ReachS(i).exclude == 0 
        ReachS(1).kin_cat_10ms = vertcat(ReachS(1).kin_cat_10ms, ReachS(i).kin_10ms(300:500,:));
    end
end

% pull out predictors 
pos_x = ReachS(1).kin_cat_10ms(:,2);
pos_y = ReachS(1).kin_cat_10ms(:,3);
pos_z = ReachS(1).kin_cat_10ms(:,4);
vel = ReachS(1).kin_cat_10ms(:,5);

vel_x = ReachS(1).kin_cat_10ms(:,6);
vel_x_up_la = ReachS(1).kin_cat_10ms(:,6) >= 0;
vel_x_up = ReachS(1).kin_cat_10ms(:,6).*vel_x_up_la;
vel_x_down_la = ReachS(1).kin_cat_10ms(:,6) < 0;
vel_x_down = abs(ReachS(1).kin_cat_10ms(:,6).*vel_x_down_la);

clear vel_x_up_la
clear vel_x_down_la

vel_y = ReachS(1).kin_cat_10ms(:,7);
vel_y_up_la = ReachS(1).kin_cat_10ms(:,7) >= 0;
vel_y_up = ReachS(1).kin_cat_10ms(:,7).*vel_y_up_la;
vel_y_down_la = ReachS(1).kin_cat_10ms(:,7) < 0;
vel_y_down = abs(ReachS(1).kin_cat_10ms(:,7).*vel_y_down_la);

clear vel_y_up_la
clear vel_y_down_la

vel_z = ReachS(1).kin_cat_10ms(:,8);
vel_z_up_la = ReachS(1).kin_cat_10ms(:,8) >= 0;
vel_z_up = ReachS(1).kin_cat_10ms(:,8).*vel_z_up_la;
vel_z_down_la = ReachS(1).kin_cat_10ms(:,8) < 0;
vel_z_down = abs(ReachS(1).kin_cat_10ms(:,8).*vel_z_down_la);

clear vel_z_up_la
clear vel_z_down_la

acc = ReachS(1).kin_cat_10ms(:,9);

acc_x = ReachS(1).kin_cat_10ms(:,10);
acc_x_up_la = ReachS(1).kin_cat_10ms(:,10) >= 0;
acc_x_up = ReachS(1).kin_cat_10ms(:,10).*acc_x_up_la;
acc_x_down_la = ReachS(1).kin_cat_10ms(:,10) < 0;
acc_x_down = abs(ReachS(1).kin_cat_10ms(:,10).*acc_x_down_la);

clear acc_x_up_la
clear acc_x_down_la

acc_y = ReachS(1).kin_cat_10ms(:,11);
acc_y_up_la = ReachS(1).kin_cat_10ms(:,11) >= 0;
acc_y_up = ReachS(1).kin_cat_10ms(:,11).*acc_y_up_la;
acc_y_down_la = ReachS(1).kin_cat_10ms(:,11) < 0;
acc_y_down = abs(ReachS(1).kin_cat_10ms(:,11).*acc_y_down_la);

clear acc_y_up_la
clear acc_y_down_la

acc_z = ReachS(1).kin_cat_10ms(:,12);
acc_z_up_la = ReachS(1).kin_cat_10ms(:,12) >= 0;
acc_z_up = ReachS(1).kin_cat_10ms(:,12).*acc_z_up_la;
acc_z_down_la = ReachS(1).kin_cat_10ms(:,12) < 0;
acc_z_down = abs(ReachS(1).kin_cat_10ms(:,12).*acc_z_down_la);

clear acc_z_up_la
clear acc_z_down_la

Predictors = [pos_x pos_y pos_z vel vel_x vel_x_up vel_x_down vel_y vel_y_up vel_y_down vel_z vel_z_up vel_z_down acc acc_x acc_x_up acc_x_down acc_y acc_y_up acc_y_down acc_z acc_z_up acc_z_down];

clear pos_x pos_y pos_z vel vel_x vel_x_up vel_x_down vel_y vel_y_up vel_y_down vel_z vel_z_up vel_z_down acc acc_x acc_x_up acc_x_down acc_y acc_y_up acc_y_down acc_z acc_z_up acc_z_down

%% Regression

% concatenate FRs at lag
a=1;
for ii=1:100
    Responses = shuff(ii).binned10smooth(:,2);

    % regression
    [B,FitInfo] = lasso(Predictors,Responses,'Standardize',true,'CV',10,'PredictorNames',{'pos_x' 'pos_y' 'pos_z' 'vel' 'vel_x' 'vel_x_up' 'vel_x_down' 'vel_y' 'vel_y_up' 'vel_y_down' 'vel_z' 'vel_z_up' 'vel_z_down' 'acc' 'acc_x' 'acc_x_up' 'acc_x_down' 'acc_y' 'acc_y_up' 'acc_y_down' 'acc_z' 'acc_z_up' 'acc_z_down'});

    % put regression results into struct
    Regression(a).B=B;
    Regression(a).FitInfo=FitInfo;

    % find lambda min that minimized MSE
    idx4 = Regression(a).FitInfo.IndexMinMSE;

    % get coefficients at min lag min lambda
    X = Regression(a).B(:,idx4);
    coef0 = Regression(a).FitInfo.Intercept(1,idx4);

    yhat = Predictors(:,1)*X(1) + Predictors(:,2)*X(2) + Predictors(:,3)*X(3) + ...
         Predictors(:,4)*X(4) + Predictors(:,5)*X(5) + Predictors(:,6)*X(6) +...
         Predictors(:,7)*X(7) + Predictors(:,8)*X(8) + Predictors(:,9)*X(9) + ...
         Predictors(:,10)*X(10) + Predictors(:,11)*X(11) + Predictors(:,12)*X(12) + ...
         Predictors(:,13)*X(13) + Predictors(:,14)*X(14) + Predictors(:,15)*X(15) + ...
         Predictors(:,16)*X(16) + Predictors(:,17)*X(17) + Predictors(:,18)*X(18) + ...
         Predictors(:,19)*X(19) + Predictors(:,20)*X(20) + Predictors(:,21)*X(21) + ...
         Predictors(:,22)*X(22) + Predictors(:,23)*X(23) + coef0;

    % calculate R2 
    ybar = mean(Responses);
    SSresid = sum((Responses-yhat).^2);
    SStotal = sum((Responses-ybar).^2);
    rsq = 1 - (SSresid/SStotal);
    RSQ(a) = rsq;
    
    % report values
    spikeshuffdata.MSEmin(a).Lambda = Regression(a).FitInfo.Lambda;
    spikeshuffdata.MSEmin(a).DF = Regression(a).FitInfo.DF;
    spikeshuffdata.MSEmin(a).MSE = Regression(a).FitInfo.MSE;
    spikeshuffdata.MSEmin(a).rsq = rsq;
    spikeshuffdata.MSEmin(a).predictors = X;
    spikeshuffdata.MSEmin(1).predictornames = Regression(1).FitInfo.PredictorNames;
    spikeshuffdata.MSEmin(a).predictedFR = yhat;
    spikeshuffdata.MSEmin(a).actualFR = Responses;
    
    % find lambda min that at 1SE MSE
    idx4 = Regression(a).FitInfo.Index1SE;

    % get coefficients at min lag min lambda
    X = Regression(a).B(:,idx4);
    coef0 = Regression(a).FitInfo.Intercept(1,idx4);

    yhat = Predictors(:,1)*X(1) + Predictors(:,2)*X(2) + Predictors(:,3)*X(3) + ...
         Predictors(:,4)*X(4) + Predictors(:,5)*X(5) + Predictors(:,6)*X(6) +...
         Predictors(:,7)*X(7) + Predictors(:,8)*X(8) + Predictors(:,9)*X(9) + ...
         Predictors(:,10)*X(10) + Predictors(:,11)*X(11) + Predictors(:,12)*X(12) + ...
         Predictors(:,13)*X(13) + Predictors(:,14)*X(14) + Predictors(:,15)*X(15) + ...
         Predictors(:,16)*X(16) + Predictors(:,17)*X(17) + Predictors(:,18)*X(18) + ...
         Predictors(:,19)*X(19) + Predictors(:,20)*X(20) + Predictors(:,21)*X(21) + ...
         Predictors(:,22)*X(22) + Predictors(:,23)*X(23) + coef0;

    % calculate R2 
    ybar = mean(Responses);
    SSresid = sum((Responses-yhat).^2);
    SStotal = sum((Responses-ybar).^2);
    rsq = 1 - (SSresid/SStotal);
    RSQ(a) = rsq;
    
    % report values
    spikeshuffdata.MSE1SE(a).rsq = rsq;
    spikeshuffdata.MSE1SE(a).predictors = X;
    spikeshuffdata.MSE1SE(1).predictornames = Regression(1).FitInfo.PredictorNames;
    spikeshuffdata.MSE1SE(a).predictedFR = yhat;
    spikeshuffdata.MSE1SE(a).actualFR = Responses;
    
    % find lambda at DF3
    for iii = 1:length(Regression(a).FitInfo.DF)
        if Regression(a).FitInfo.DF(iii) == 3
            idx4 = iii;
            break
        end
    end
    
    % get coefficients at min lag min lambda
    X = Regression(a).B(:,idx4);
    coef0 = Regression(a).FitInfo.Intercept(1,idx4);

    yhat = Predictors(:,1)*X(1) + Predictors(:,2)*X(2) + Predictors(:,3)*X(3) + ...
         Predictors(:,4)*X(4) + Predictors(:,5)*X(5) + Predictors(:,6)*X(6) +...
         Predictors(:,7)*X(7) + Predictors(:,8)*X(8) + Predictors(:,9)*X(9) + ...
         Predictors(:,10)*X(10) + Predictors(:,11)*X(11) + Predictors(:,12)*X(12) + ...
         Predictors(:,13)*X(13) + Predictors(:,14)*X(14) + Predictors(:,15)*X(15) + ...
         Predictors(:,16)*X(16) + Predictors(:,17)*X(17) + Predictors(:,18)*X(18) + ...
         Predictors(:,19)*X(19) + Predictors(:,20)*X(20) + Predictors(:,21)*X(21) + ...
         Predictors(:,22)*X(22) + Predictors(:,23)*X(23) + coef0;

    % calculate R2 
    ybar = mean(Responses);
    SSresid = sum((Responses-yhat).^2);
    SStotal = sum((Responses-ybar).^2);
    rsq = 1 - (SSresid/SStotal);
    
    % report values
    spikeshuffdata.MSEDF3(a).rsq = rsq;
    spikeshuffdata.MSEDF3(a).predictors = X;
    spikeshuffdata.MSEDF3(1).predictornames = Regression(1).FitInfo.PredictorNames;
    spikeshuffdata.MSEDF3(a).predictedFR = yhat;
    spikeshuffdata.MSEDF3(a).actualFR = Responses;
    
    a=a+1;
end

