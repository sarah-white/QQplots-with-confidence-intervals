% This code makes quantile-quantile (QQ) plots, which are a method to 
% compare two distributions. The code is oriented toward temperature
% datasets, but they could be anything (d18O, etc).
% It also makes normalized QQ plots, normalized to the 1:1 line 
% offset by the difference in mean between the two datasets.
% Confidence intervals are estimated by Monte Carlo bootstrapping, 
% and are plotted as a gray region.
% The code is designed to make many QQ and normalized QQ plots, to
% compare many datasets to a single reference dataset. Thus, the x axis
% dataset is the same on every plot.

% This code was written by Sarah White at UC Santa Cruz (smwhite@ucsc.edu), 
% for the following paper (please cite):
% White, S. M., Ravelo, A. C., & Polissar, P. J. (2018). 
% Dampened El Niño in the Early and Mid?Holocene Due To Insolation?Forced 
% Warming/Deepening of the Thermocline. Geophysical Research Letters, 
% 45(1), 316-326. https://doi.org/10.1002/2017GL075433

% The supplementary materials of the paper fully describe what the code
% does. There are also comments within the code.

% INPUT
% The code is written for paleoclimate data in which there are several 
% samples of different ages, each of which has an associated distribution 
% of temperature (or d18O, etc). 
% The code expects two Excel spreadsheets, one for the x axis dataset 
% (against which all other datsets will be plotted), and one for the y axis
% datasets. Y axis datasets should all be strung together in two columns,
% age and temperature. See example spreadsheet,
% "White_etal_2018_singleforamdata.xlsx".

% OUTPUT
% The code creates a QQ and a normalized QQ plot for every sample (actually,
% for every unique age in the y axis dataset). The plots all include
% confidence intervals, at a level set by the user (90%, 95%, etc).

clear all
close all

%%%%% User input
% What dataset do you want as the x axis? (coretop, modern temperatures, etc)?
% Set data pathnames below (filename, then tab name)
% The x axis dataset will be called "xaxis_temps" and the y axis dataset
% will be called "yaxis_alltemps"

%%%%% Get data
xaxis_temps=xlsread('White_etal_2018_singleforamdata.xlsx', 'coretop');
% What is your x axis dataset called?
x_axis_label='coretop {\circ}C';
% check which column your temperature data is in (1st, 2nd etc) and set below
xaxis_temps=xaxis_temps(:,2);
% Set filename for y axis data
yaxis_data=xlsread('White_etal_2018_singleforamdata.xlsx', 'downcore');
% check which column your temperature data is in
yaxis_alltemps=yaxis_data(:,2);

%%%%% User-set parameters
number_of_quantiles=50;                                         
% Number of Monte Carlo (bootstrap) datasets to create
nMC=10000;                                                      
confidence_interval1=95;
confidence_interval2=80;
% Go to plotting code to adjust x and y axis limits
% (lines 139-142 and 174-183)

%%%%% No more user adjustments below here!!
steps=(0 : 1/(number_of_quantiles) : 1);
nages=length(unique(yaxis_data(:,1)));
ages=unique(yaxis_data(:,1));

%%%%% Calculate x axis quantiles
% Sort the data, so that when plotted on the CDF, it will go from lowest to highest
xaxis_temps=sort(xaxis_temps); 
% In the next two lines of code, each new data point of xaxis_temps brings the 
% cumulative % of the data one step closer to 100%, with step sizes being 
% 1/length of the dataset. It is crucial that xaxis_temps is sorted;
% this way, when I plot it versus CDF_steps_xaxis (which is just equally 
% spaced numbers from 0 to 1), it comes out as a CDF 
CDF_steps_xaxis=[0 : 1/(length(xaxis_temps)-1) : 1]; 
% "xaxis_quantiles" are the quantile values of the data, 
% calculated by linear interpolation (default) between data points, 
% at intervals given by "number_of_quantiles" (e.g. 2%, 4%,... 98%)
xaxis_quantiles=interp1(CDF_steps_xaxis, xaxis_temps, steps);  

%%%%% Calculate quantiles, do Monte Carlo simulations, and make figures,
%%%%% looped through all downcore intervals
for a=1:nages
    % This finds the number of points with a given age
    i=find(yaxis_data(:,1)==ages(a));                            
    sacc_temps_interval=sort(yaxis_alltemps(i));
    % Create steps for the CDF
    CDF_steps=[0 : 1/(length(i)-1) : 1]; 
    % Here I calculate quantile values for data of a given age, 
    % in the same manner as for the x axis data
    quantiles=interp1(CDF_steps, sacc_temps_interval, steps);   
    randnums=rand(1,nMC*length(i));
    % Here I interpolate along the empirical CDF at random points, 
    % to generate random data points based on my data's distribution
    % (bootstrapping)
    MC_data=interp1(CDF_steps, sacc_temps_interval, randnums);
    % Breaks the big vector of bootstrapped data into nMC sets, with length
    % equal to the original dataset
    MC_data=reshape(MC_data, nMC, length(i));
    % Sorts each row, so we can calculate quantiles from each Monte Carlo bootstrapped
    % dataset
    MC_data=sort(MC_data,2);                                    
    MC_quantiles=zeros(nMC,number_of_quantiles+1);
    for b=1:nMC
        MC_quantiles(b,:)=interp1(CDF_steps, MC_data(b,:), steps);      
    end
    
    %%%%% Calculate confidence intervals
    % Sort each column
    MC_quantiles=sort(MC_quantiles,1);  
    % Choose confidence intervals from sorted Monte Carlo datasets. For
    % example, for 10,000 Monte Carlo datasets, 95% confidence interval
    % comes from 250th and 9750th values, for each quantile
    confidence_lowest(1,:)=MC_quantiles(round(((1-confidence_interval1/100)/2)*nMC),:);
    confidence_highest(1,:)=MC_quantiles(round((1-((1-confidence_interval1/100)/2))*nMC),:);
    confidence_midlo(1,:)=MC_quantiles(round(((1-confidence_interval2/100)/2)*nMC),:);
    confidence_midhi(1,:)=MC_quantiles(round((1-((1-confidence_interval2/100)/2))*nMC),:);
    
    %%%%% Make figures
    %%%%% QQ plots
    figure
    % scale the x axis for the normalized plots to the x axis dataset 
    norm_xmin=floor(min(xaxis_temps));
    norm_xmax=ceil(max(xaxis_temps));
    % scale the x and y axes for the QQ plots to the extremes of the downcore data
    QQ_xmin=floor(min(yaxis_alltemps));
    QQ_xmax=ceil(max(yaxis_alltemps));
    
    % Make plots
    % This makes the "out-and-back" x values for making a filled plot
    coretop_fillplot=[xaxis_quantiles(2:end-1), fliplr(xaxis_quantiles(2:end-1))];
    % This makes the "out-and-back" y values for making a filled plot
    confidence_95=[confidence_lowest(2:end-1), fliplr(confidence_highest(2:end-1))];      
    fill(coretop_fillplot, confidence_95, [0.75 0.75 0.75], 'EdgeColor', 'none')
    % Do you want axis limits defined by edges of dataset, or given values?
    % Comment out which of the following two lines of code you don't want
    xlim([QQ_xmin-1 QQ_xmax+1]), ylim([QQ_xmin-1 QQ_xmax+1])
    %xlim([15 29]), ylim([15 29])
    set(gca, 'fontsize', 20)
    hold on
    plot(xaxis_quantiles(2:end-1), quantiles(2:end-1), 'dk', 'MarkerFaceColor', 'k')
    xlabel(x_axis_label), ylabel([num2str(ages(a)), ' ka {\circ}C'])
    % Make 1:1 line (solid) and offset 1:1 line (dashed) (offset by difference in average
    % value between datasets)
    plot((QQ_xmin-1 :1: QQ_xmax+1),(QQ_xmin-1 :1: QQ_xmax+1), 'k','linewidth', 2)
    avg_offset=mean(sacc_temps_interval)-mean(xaxis_temps);
    plot((QQ_xmin-1 :1: QQ_xmax+1),(QQ_xmin-1+avg_offset :1: QQ_xmax+1+avg_offset), '--k')
    set(gca, 'linewidth', 2)
    hold off
    
    %%%%% normalized QQ plots
    figure
    % This gives the y value of the dashed line at every y point. 
    % It's equal to the x value of the dashed line (which is the value of 
    % the coretop dataset at that point) plus the mean offset between x and y datasets
    norm_1to1=xaxis_quantiles+avg_offset; 
    % Here we subtract the y value of the dashed line from the y dataset 
    norm_quantiles=quantiles-norm_1to1;         
    norm_conf_highest=confidence_highest-norm_1to1;
    norm_conf_lowest=confidence_lowest-norm_1to1;
    norm_conf95=[norm_conf_lowest(2:end-1), fliplr(norm_conf_highest(2:end-1))];
    % Make gray region for confidence interval
    fill(coretop_fillplot, norm_conf95, [0.75 0.75 0.75], 'EdgeColor', 'none')
    set(gca, 'fontsize', 24)
    hold on
    % Make the zero line
    plot(norm_xmin-1 :1: norm_xmax+1, zeros(1,(norm_xmax-norm_xmin)+3), '-k','linewidth', 2)
    plot(xaxis_quantiles(2:end-1), norm_quantiles(2:end-1), 'dk', 'MarkerFaceColor', 'k')
    xlabel(x_axis_label), ylabel(['normalized ', num2str(ages(a)), ' ka {\circ}C'])
    % Only use the next two lines of code if you want reversed axes, as for
    % d18O data
    %set(gca, 'xdir', 'reverse')
    %set(gca, 'ydir', 'reverse')
    % Do you want axis limits defined by edges of dataset, or given values?
    % Comment out which of the following two lines of code you don't want
    xlim([norm_xmin norm_xmax])
    %xlim([14 27])
    %xticks([14 16 18 20 22 24 26]);
    ylim([-1.5 1.5])
    set(gca, 'linewidth', 2)
end

