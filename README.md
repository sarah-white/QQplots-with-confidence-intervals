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

