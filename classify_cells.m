
function [fr_summary] = classify_cells(path_loc) 

    cd(path_loc)
    directory = dir('*.mat');
    
    n = 1;
    for rec = 1:length(directory)
       
       % load data from NP files
       if contains(convertCharsToStrings(directory(rec).name),'NP')
    
           load(directory(rec).name,'cellData')
       
           for cell = 1:length(cellData)
    
               % only analyze cells marked as good
               if cellData(cell).quality == "good " 
    
                    Spikes.SS_Bin1 = cellData(cell).Bin1;
                    Spikes.SS_Bin10 = cellData(cell).Bin10;
                    Spikes.SS_Bin10smooth = cellData(cell).Bin10smooth;
    
                    % make raster
                    Spikes = get_raster(Spikes);
    
                    % calculate FR
                   fr = sum(Spikes.SS_Bin1(:,2))/(Spikes.SS_Bin1(end,1)-Spikes.SS_Bin1(1,1));
    
                   % calculate MAD ISI
                   medianISI = median(Spikes.SS_Raster(2:end,2));
                   ISIdiff = abs(Spikes.SS_Raster(2:end,2)-medianISI);
                   MAD = median(ISIdiff);
                   
                   % calculate CV
                   mISI = mean(Spikes.SS_Raster(2:end,2));
                   stdISI = std(Spikes.SS_Raster(2:end,2));
                   CV = stdISI/mISI;
    
                   % calculate CV2
                   for spike = 1:length(Spikes.SS_Raster)-1
                        CV2(spike) = abs(Spikes.SS_Raster(spike+1,2) - Spikes.SS_Raster(spike,2))/((Spikes.SS_Raster(spike+1,2) + Spikes.SS_Raster(spike,2))/2);
                   end
                   CV2 = mean(CV2);
    
                   % enter data
                   fr_summary(n).filename = directory(rec).name;
                   fr_summary(n).cellnum = cell;
                   fr_summary(n).MAD = MAD;
                   fr_summary(n).CV = CV;
                   fr_summary(n).CV2 = CV2;
                   fr_summary(n).fr = fr;
                   fr_summary(n).depth = cellData(cell).depth;
                   
                   % add PC layer depth if it was specified
                   if exist('PClayerdepth')
                        if cellData(cell).depth>=PClayerdepth
                            fr_summary(n).abovePClayerdepth = 1;
                        else
                            fr_summary(n).abovePClayerdepth = 0;
                        end
                   else
                       fr_summary(n).abovePClayerdepth = 1;
                   end
    
                   % check if cspk sorted
                   if isfield(cellData,'CS_Bin1')
                       if ~isempty(cellData(cell).CS_Bin1)
                           fr_summary(n).csSorted = 1;
                           fr_summary(n).csfr = sum(cellData(cell).CS_Bin1(:,2))/(cellData(cell).CS_Bin1(end,1)-cellData(cell).CS_Bin1(1,1));
                       else 
                           fr_summary(n).csSorted = 0; 
                       end
                   else
                       fr_summary(n).csSorted = 0;   
                   end
                   
                   clear CV2 CV fr csfr     
                   n = n+1;
               end
           end
    
       clear cellData
    
       % load data from psort files
       else 
           load(directory(rec).name)
       
           % calculate FR
           fr = sum(Spikes.SS_Bin1(:,2))/(Spikes.SS_Bin1(end,1)-Spikes.SS_Bin1(1,1));
    
           % calculate MAD ISI
           medianISI = median(Spikes.SS_Raster(2:end,2));
           ISIdiff = abs(Spikes.SS_Raster(2:end,2)-medianISI);
           MAD = median(ISIdiff);
                      
           % calculate CV
           mISI = mean(Spikes.SS_Raster(2:end,2));
           stdISI = std(Spikes.SS_Raster(2:end,2));
           CV = stdISI/mISI;
    
           % calculate CV2
           for spike = 1:length(Spikes.SS_Raster)-1
                CV2(spike) = abs(Spikes.SS_Raster(spike+1,2) - Spikes.SS_Raster(spike,2))/((Spikes.SS_Raster(spike+1,2) + Spikes.SS_Raster(spike,2))/2);
           end
           CV2 = mean(CV2);
    
           % enter data
           fr_summary(n).filename = directory(rec).name;
           fr_summary(n).cellnum = 1;
           fr_summary(n).MAD = MAD;
           fr_summary(n).CV = CV;
           fr_summary(n).CV2 = CV2;
           fr_summary(n).fr = fr;
           fr_summary(n).abovePClayerdepth = 1;
    
           if isfield(Spikes,'CS_Bin1')
               fr_summary(n).csSorted = 1;
               fr_summary(n).csfr = sum(Spikes.CS_Bin1(:,2))/(Spikes.CS_Bin1(end,1)-Spikes.CS_Bin1(1,1));
           else
               fr_summary(n).csSorted = 0;   
           end
    
           clear CV2 CV fr csfr Spikes MAD
           n = n+1;
           
       end              
    end
       
    %% classify each cell 
    n = 1; m = 1; o = 1;
    for i = 1:length(fr_summary)
    
        % classify as PC if Cspk is present
        if fr_summary(i).csSorted == 1
            fr_grouped.PC_cpsk(n,1) = fr_summary(i).MAD;
            fr_grouped.PC_cpsk(n,2) = fr_summary(i).CV2;
            fr_grouped.PC_cpsk(n,3) = fr_summary(i).fr;
            fr_summary(i).cellID = "PC";
            n = n+1;
    
        % classify as PC if fr criteria is met
        elseif fr_summary(i).CV2 > 0.2 && round(fr_summary(i).MAD,3)<0.008 && fr_summary(i).fr>=40 
            fr_grouped.PC_criteria(m,1) = fr_summary(i).MAD;
            fr_grouped.PC_criteria(m,2) = fr_summary(i).CV2;
            fr_grouped.PC_criteria(m,3) = fr_summary(i).fr;
            fr_summary(i).cellID = "PC";
            m= m+1;
    
        % othewise classify as other
        else
            fr_grouped.other(o,1) = fr_summary(i).MAD;
            fr_grouped.other(o,2) = fr_summary(i).CV2;
            fr_grouped.other(o,3) = fr_summary(i).fr;
            o= o+1;
        end
        
    end
    
    
    %% run tSNE

    % load initY for reproducibility
    load("tSNE_initY.mat")

    wCS_ID(1:length(fr_grouped.PC_cpsk),1) = {'PC-CS'};
    wCrit_ID(1:length(fr_grouped.PC_criteria),1) = {'PC-crit'};
    O_ID(1:length(fr_grouped.other),1) = {'other'};
    ID_type = [wCS_ID;wCrit_ID;O_ID];
    tot = [fr_grouped.PC_cpsk;fr_grouped.PC_criteria;fr_grouped.other];
    
    [Y,loss] = tsne(tot,'Algorithm','exact','Distance','euclidean','standardize',true,'LearnRate',5000,'InitialY',initY);
    figure
    gscatter(Y(:,1),Y(:,2),ID_type)
    
end
