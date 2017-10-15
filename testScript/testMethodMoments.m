% test method of moment
% run script
clc; close all; clear variables;
addpath(genpath(pwd))
%% 1) get data from Yahoo
freqIn = 'wk'; % time-series frequency 
freqOut = 12; % return frequency (quarterly) used in the model
Returns = getTimeSeries(freqIn,freqOut); % get time-series from Yahoo
% Returns = getTimeSeries(freqIn); % get time-series from Yahoo
minLength = min(cellfun('length',Returns)); % get the minimum length of each cell
Returns = cellfun(@(x) x(1:minLength),Returns,'UniformOutput',false);
Returns = cell2mat(Returns); % convert cell to matrix

[x] = MethodMomentsGM(Returns);