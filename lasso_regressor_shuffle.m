function [realdata_regshuff] = lasso_regressor_shuffle(ReachS,Bin10smooth,lagmax,empirical_predictors)

% lasso regression model

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
vel_x_down = ReachS(1).kin_cat_10ms(:,6).*vel_x_down_la;

clear vel_x_up_la
clear vel_x_down_la

vel_y = ReachS(1).kin_cat_10ms(:,7);
vel_y_up_la = ReachS(1).kin_cat_10ms(:,7) >= 0;
vel_y_up = ReachS(1).kin_cat_10ms(:,7).*vel_y_up_la;
vel_y_down_la = ReachS(1).kin_cat_10ms(:,7) < 0;
vel_y_down = ReachS(1).kin_cat_10ms(:,7).*vel_y_down_la;

clear vel_y_up_la
clear vel_y_down_la

vel_z = ReachS(1).kin_cat_10ms(:,8);
vel_z_up_la = ReachS(1).kin_cat_10ms(:,8) >= 0;
vel_z_up = ReachS(1).kin_cat_10ms(:,8).*vel_z_up_la;
vel_z_down_la = ReachS(1).kin_cat_10ms(:,8) < 0;
vel_z_down = ReachS(1).kin_cat_10ms(:,8).*vel_z_down_la;

clear vel_z_up_la
clear vel_z_down_la

acc = ReachS(1).kin_cat_10ms(:,9);

acc_x = ReachS(1).kin_cat_10ms(:,10);
acc_x_up_la = ReachS(1).kin_cat_10ms(:,10) >= 0;
acc_x_up = ReachS(1).kin_cat_10ms(:,10).*acc_x_up_la;
acc_x_down_la = ReachS(1).kin_cat_10ms(:,10) < 0;
acc_x_down = ReachS(1).kin_cat_10ms(:,10).*acc_x_down_la;

clear acc_x_up_la
clear acc_x_down_la

acc_y = ReachS(1).kin_cat_10ms(:,11);
acc_y_up_la = ReachS(1).kin_cat_10ms(:,11) >= 0;
acc_y_up = ReachS(1).kin_cat_10ms(:,11).*acc_y_up_la;
acc_y_down_la = ReachS(1).kin_cat_10ms(:,11) < 0;
acc_y_down = ReachS(1).kin_cat_10ms(:,11).*acc_y_down_la;

clear acc_y_up_la
clear acc_y_down_la

acc_z = ReachS(1).kin_cat_10ms(:,12);
acc_z_up_la = ReachS(1).kin_cat_10ms(:,12) >= 0;
acc_z_up = ReachS(1).kin_cat_10ms(:,12).*acc_z_up_la;
acc_z_down_la = ReachS(1).kin_cat_10ms(:,12) < 0;
acc_z_down = ReachS(1).kin_cat_10ms(:,12).*acc_z_down_la;

clear acc_z_up_la
clear acc_z_down_la

Predictors = [pos_x pos_y pos_z vel vel_x vel_x_up vel_x_down vel_y vel_y_up vel_y_down vel_z vel_z_up vel_z_down acc acc_x acc_x_up acc_x_down acc_y acc_y_up acc_y_down acc_z acc_z_up acc_z_down];

clear pos_x pos_y pos_z vel vel_x vel_x_up vel_x_down vel_y vel_y_up vel_y_down vel_z vel_z_up vel_z_down acc acc_x acc_x_up acc_x_down acc_y acc_y_up acc_y_down acc_z acc_z_up acc_z_down

for i = 1:length(empirical_predictors)
    if empirical_predictors(i,1) == 0
        Predictors(:,i) = Predictors(:,i)*0;
    end
end

%% Regression

% find aligned FR start and end points
for i=1:length(ReachS)
    if ReachS(i).exclude == 0 
        idx1(i,1)=knnsearch(Bin10smooth(:,1),ReachS(i).kin_10ms(300,1));
        idx2(i,1)=idx1(i,1)+200;
    end
end

a=1;
ii=lagmax/10;
step_size=1;
for n = 1:size(Predictors,2)
    % check if non-zero
    if sum(Predictors(:,n)) ~= 0
        % shuffle individual regressors
        Predictors_shuff = Predictors;
        col = Predictors(:,n);
        Predictors_shuff(:,n) = col(randperm(length(col)));

        % concatenate FRs at lag
        Responses = [];
        for i=1:length(ReachS)
            if ReachS(i).exclude == 0
                Responses = vertcat(Responses, Bin10smooth(idx1(i,1)+(ii*step_size):idx2(i,1)+(ii*step_size),2));
            end
        end

        % regression
        [B,FitInfo] = lasso(Predictors_shuff,Responses,'Standardize',true,'CV',10,'PredictorNames',{'pos_x' 'pos_y' 'pos_z' 'vel' 'vel_x' 'vel_x_up' 'vel_x_down' 'vel_y' 'vel_y_up' 'vel_y_down' 'vel_z' 'vel_z_up' 'vel_z_down' 'acc' 'acc_x' 'acc_x_up' 'acc_x_down' 'acc_y' 'acc_y_up' 'acc_y_down' 'acc_z' 'acc_z_up' 'acc_z_down'});

        % put regression results into struct
        Regression(n).B=B;
        Regression(n).FitInfo=FitInfo;

        % find lambda min that minimized MSE
        [~,idx4] = min(Regression(n).FitInfo.MSE);

        % get coefficients at min lag min lambda
        X = Regression(n).B(:,idx4);
        coef0 = Regression(n).FitInfo.Intercept(1,idx4);

        yhat = Predictors_shuff(:,1)*X(1) + Predictors_shuff(:,2)*X(2) + Predictors_shuff(:,3)*X(3) + ...
             Predictors_shuff(:,4)*X(4) + Predictors_shuff(:,5)*X(5) + Predictors_shuff(:,6)*X(6) +...
             Predictors_shuff(:,7)*X(7) + Predictors_shuff(:,8)*X(8) + Predictors_shuff(:,9)*X(9) + ...
             Predictors_shuff(:,10)*X(10) + Predictors_shuff(:,11)*X(11) + Predictors_shuff(:,12)*X(12) + ...
             Predictors_shuff(:,13)*X(13) + Predictors_shuff(:,14)*X(14) + Predictors_shuff(:,15)*X(15) + ...
             Predictors_shuff(:,16)*X(16) + Predictors_shuff(:,17)*X(17) + Predictors_shuff(:,18)*X(18) + ...
             Predictors_shuff(:,19)*X(19) + Predictors_shuff(:,20)*X(20) + Predictors_shuff(:,21)*X(21) + ...
             Predictors_shuff(:,22)*X(22) + Predictors_shuff(:,23)*X(23) + coef0;

        % calculate R2 
        ybar = mean(Responses);
        SSresid = sum((Responses-yhat).^2);
        SStotal = sum((Responses-ybar).^2);
        rsq = 1 - (SSresid/SStotal);

        % report data
        realdata_regshuff.MSEmin(n).Lambda = Regression(n).FitInfo.Lambda;
        realdata_regshuff.MSEmin(n).DF = Regression(n).FitInfo.DF;
        realdata_regshuff.MSEmin(n).MSE = Regression(n).FitInfo.MSE;
        realdata_regshuff.MSEmin(n).rsq = rsq;
        realdata_regshuff.MSEmin(n).predictors = X;
        realdata_regshuff.MSEmin(n).shuffledpredictor = Regression(n).FitInfo.PredictorNames(n);
        realdata_regshuff.MSEmin(n).predictedFR = yhat;
        realdata_regshuff.MSEmin(n).actualFR = Responses;

        % find lambda min that at 1SE MSE
        idx4 = Regression(n).FitInfo.Index1SE;

        % get coefficients at min lag min lambda
        X = Regression(n).B(:,idx4);
        coef0 = Regression(n).FitInfo.Intercept(1,idx4);

        yhat = Predictors_shuff(:,1)*X(1) + Predictors_shuff(:,2)*X(2) + Predictors_shuff(:,3)*X(3) + ...
             Predictors_shuff(:,4)*X(4) + Predictors_shuff(:,5)*X(5) + Predictors_shuff(:,6)*X(6) +...
             Predictors_shuff(:,7)*X(7) + Predictors_shuff(:,8)*X(8) + Predictors_shuff(:,9)*X(9) + ...
             Predictors_shuff(:,10)*X(10) + Predictors_shuff(:,11)*X(11) + Predictors_shuff(:,12)*X(12) + ...
             Predictors_shuff(:,13)*X(13) + Predictors_shuff(:,14)*X(14) + Predictors_shuff(:,15)*X(15) + ...
             Predictors_shuff(:,16)*X(16) + Predictors_shuff(:,17)*X(17) + Predictors_shuff(:,18)*X(18) + ...
             Predictors_shuff(:,19)*X(19) + Predictors_shuff(:,20)*X(20) + Predictors_shuff(:,21)*X(21) + ...
             Predictors_shuff(:,22)*X(22) + Predictors_shuff(:,23)*X(23) + coef0;

        % calculate R2 
        ybar = mean(Responses);
        SSresid = sum((Responses-yhat).^2);
        SStotal = sum((Responses-ybar).^2);
        rsq = 1 - (SSresid/SStotal);

        % report data
        realdata_regshuff.MSE1SE(n).rsq = rsq;
        realdata_regshuff.MSE1SE(n).predictors = X;
        realdata_regshuff.MSE1SE(n).shuffledpredictor = Regression(n).FitInfo.PredictorNames(n);
        realdata_regshuff.MSE1SE(n).predictedFR = yhat;
        realdata_regshuff.MSE1SE(n).actualFR = Responses;

        % find lambda min at MSE DF3
        for iii = 1:length(Regression(n).FitInfo.DF)
            if Regression(n).FitInfo.DF(iii) == 3
                idx4 = iii;
                break
            end
        end

        % get coefficients at min lag min lambda
        X = Regression(n).B(:,idx4);
        coef0 = Regression(n).FitInfo.Intercept(1,idx4);

        yhat = Predictors_shuff(:,1)*X(1) + Predictors_shuff(:,2)*X(2) + Predictors_shuff(:,3)*X(3) + ...
             Predictors_shuff(:,4)*X(4) + Predictors_shuff(:,5)*X(5) + Predictors_shuff(:,6)*X(6) +...
             Predictors_shuff(:,7)*X(7) + Predictors_shuff(:,8)*X(8) + Predictors_shuff(:,9)*X(9) + ...
             Predictors_shuff(:,10)*X(10) + Predictors_shuff(:,11)*X(11) + Predictors_shuff(:,12)*X(12) + ...
             Predictors_shuff(:,13)*X(13) + Predictors_shuff(:,14)*X(14) + Predictors_shuff(:,15)*X(15) + ...
             Predictors_shuff(:,16)*X(16) + Predictors_shuff(:,17)*X(17) + Predictors_shuff(:,18)*X(18) + ...
             Predictors_shuff(:,19)*X(19) + Predictors_shuff(:,20)*X(20) + Predictors_shuff(:,21)*X(21) + ...
             Predictors_shuff(:,22)*X(22) + Predictors_shuff(:,23)*X(23) + coef0;

        % calculate R2 
        ybar = mean(Responses);
        SSresid = sum((Responses-yhat).^2);
        SStotal = sum((Responses-ybar).^2);
        rsq = 1 - (SSresid/SStotal);

        % report data
        realdata_regshuff.MSEDF3(n).rsq = rsq;
        realdata_regshuff.MSEDF3(n).predictors = X;
        realdata_regshuff.MSEDF3(n).shuffledpredictor = Regression(n).FitInfo.PredictorNames(n);
        realdata_regshuff.MSEDF3(n).predictedFR = yhat;
        realdata_regshuff.MSEDF3(n).actualFR = Responses;
    end
end

