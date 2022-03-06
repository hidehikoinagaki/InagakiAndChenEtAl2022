function [] = plotNWBdata_PTsilencing()
% original code to plot NWB data @ https://neurodatawithoutborders.github.io/matnwb/tutorials/html/basicUsage.html
% modified the original code by Hidehiko to plot PT silencing data 2/23/2022 


% add matnwb to your path 
addpath('matnwb-2.4.0.0');

% example nwb file to read
nwb = nwbRead('out/SC009_112118.nwb');


%% read nwb file
unit_names = keys(nwb.analysis);

unit_ids = nwb.units.id.data.load(); % array of unit ids represented within this
% Initialize trials & times Map containers indexed by unit_ids
unit_trials = containers.Map('KeyType',class(unit_ids),'ValueType','any');
unit_times  = containers.Map('KeyType',class(unit_ids),'ValueType','any');
last_idx = 0;
for i = 1:length(unit_ids)
    unit_id = unit_ids(i);
    
    row = nwb.units.getRow(unit_id, 'useId', true, 'columns', {'spike_times', 'trialsID'});
    unit_trials(unit_id) = row.trialsID{1};
    unit_times(unit_id)  = row.spike_times{1};
end

sorted_ids = sort(unit_ids);
Photostim = struct(...
    'ind', true,... % mask into xs and ys for this photostim
    'name', 'none',...
    'stim_duration', 0,... % in seconds after the onset of normal go cue
    'stim_onset', 0,... % in seconds after the onset of normal go cue
    'gocue_onset', 0,...% in seconds after the onset of normal go cue
    'lickL_trials',0,... % sum of lick L trials
    'lickR_trials',0); % sum of lick R trials


% Initialize Map container of plotting data for each unit, stored as structure
Unit = containers.Map('KeyType',class(unit_ids),'ValueType','any');
unit_struct = struct(...
    'id', [],...
    'xs', [],...
    'ys', [],...
    'xlim', [-Inf Inf],...
    'sample', 0,...
    'delay', 0,...
    'response', 0,...
    'left_scatter', false,...
    'right_scatter', false,...
    'photostim', Photostim); % can have multiple photostim



% read data from indv units
for unit_id = unit_ids'
       
    unit_trial_id = unit_trials(unit_id);
    
    trial = nwb.intervals_trials.getRow(unit_trial_id, 'useId', true,...
        'columns', {'SampleOnset', 'DelayOnset', 'CueTime', 'GoodTrials', 'HitR', 'HitL',...
        'ErrR','ErrL','NoLickR','NoLickL','StimTrials', 'PhotostimulationType'});
    unit_sample      = trial.SampleOnset;
    unit_delay       = trial.DelayOnset;
    unit_response    = trial.CueTime;
    unit_good_trials = logical(trial.GoodTrials);
    
    % calculate mean sample/delay/response time based on good trials w.o.
    % early licks
    good_trials_for_timing = unit_good_trials &  (logical(trial.HitR) | logical(trial.HitL)) & ~logical(trial.StimTrials);
    avg_sample = mean(unit_sample(good_trials_for_timing));
    avg_delay = mean(unit_delay(good_trials_for_timing));
    avg_response = mean(unit_response(good_trials_for_timing));

    unit_trial_id = unit_trial_id(unit_good_trials);
    
    unit_spike_time = unit_times(unit_id);
    unit_spike_time = unit_spike_time(unit_good_trials);

    unit_stim_type = trial.PhotostimulationType;
       % stim type
%         '"0"--non-stimulation trials',...
%         '"7"--  PT silencing 0.5mW start at the Go cue for 1.0s',...
    
    unit_no_stim    = unit_good_trials & 0 == unit_stim_type;
    unit_stim       = unit_good_trials & 7 == unit_stim_type; % Perturbation PT lower 0.5mW
    
    
    % count number of trials per condition
    % we need to do this as there could be trial w.o. spikes
    first_trial = min(unit_trial_id); 
    last_trial  = max(unit_trial_id);
    
    trial_in_range = nwb.intervals_trials.getRow(first_trial:last_trial, 'useId', true,...
        'columns', {'GoodTrials', 'HitR', 'HitL',...
        'ErrR','ErrL','NoLickR','NoLickL','StimTrials', 'PhotostimulationType'});
    
    no_stim_L = sum((trial_in_range.HitL | trial_in_range.ErrL | trial_in_range.NoLickL) &...
        trial_in_range.GoodTrials & trial_in_range.PhotostimulationType == 0);
    no_stim_R = sum((trial_in_range.HitR | trial_in_range.ErrR | trial_in_range.NoLickR) &...
        trial_in_range.GoodTrials & trial_in_range.PhotostimulationType == 0);

    % Perturbation PT lower 0.5mW
    stim_L = sum((trial_in_range.HitL | trial_in_range.ErrL | trial_in_range.NoLickL) &...
        trial_in_range.GoodTrials & trial_in_range.PhotostimulationType == 7);
    stim_R = sum((trial_in_range.HitR | trial_in_range.ErrR | trial_in_range.NoLickR) &...
        trial_in_range.GoodTrials & trial_in_range.PhotostimulationType == 7);
    
    
    % summarize spike info for plotting
    
    xs = unit_spike_time - avg_response; % align to cue
    ys = unit_trial_id;
    
    curr_unit = unit_struct;
    
    curr_unit.xs = xs;
    curr_unit.ys = ys;
    curr_unit.left_scatter  = logical(trial.HitL | trial.ErrL | trial.NoLickL); % lick left trial spikes
    curr_unit.right_scatter = logical(trial.HitR | trial.ErrR | trial.NoLickR); % lick right trials spikes
    curr_unit.sample = avg_sample - avg_response;
    curr_unit.delay = avg_delay - avg_response;
    curr_unit.response = 0;

    % Photostim periods
    curr_unit.photostim.ind = unit_no_stim;
    curr_unit.photostim.lickL_trials = no_stim_L;
    curr_unit.photostim.lickR_trials = no_stim_R;
    
    % PT silencing
    if any(unit_stim)
        Stim = Photostim;
        Stim.ind = unit_stim;
        Stim.name = 'PT silencing';
        Stim.stim_duration = 1.2;
        Stim.stim_onset    = 0.0; % to go cue (s)
        Stim.gocue_onset   = 0.0; 
        Stim.lickL_trials = stim_L;
        Stim.lickR_trials = stim_R;
        curr_unit.photostim(end+1) = Stim;
    end


    Unit(unit_id) = curr_unit;
end
    

%% Plot Example Neurons  
neuron_ids   = 1:unit_id; 
num_conditions = 2; % photostim conditions: nostim, PPN perturbation
num_neurons = length(neuron_ids);
% Inititalize data structures for each summary plot of categorized neural spike data at specified stimulus condition
RasterPlot = struct(...
    'xs', 0,...
    'ys', 0);
ConditionPlot = struct(...
    'label', '',...
    'xlim', 0,...
    'sample', 0,...
    'delay', 0,...
    'response', 0,...
    'right_scatter', RasterPlot,...
    'left_scatter', RasterPlot,...
    'psth_bin_window', 0,...
    'stim_type', '');

    
% Plot neural spike data for each neuron and stimulus condition in a subplot array: num_neurons (rows) x num_conditions (columns)
for nn=1:num_neurons
    
    Neuron = Unit(neuron_ids(nn));
    
    % Initialize structure with neural + stimulus condition data
    CurrPlot = ConditionPlot;
    CurrPlot.xlim = [-4 3];
    CurrPlot.sample = Neuron.sample;
    CurrPlot.delay = Neuron.delay;   
    CurrPlot.psth_bin_window = 9;
    
    % Plot each neuron/condition
    figure;
    
    for cc=1:num_conditions        
        Stim = Neuron.photostim(cc);
        
        CurrPlot.stim_type = Stim.name;
        if strcmp(Stim.name, 'none')
            CurrPlot.label = sprintf('Neuron %d', neuron_ids(nn));    
        else
            CurrPlot.label = Stim.name;
        end
        
        CurrPlot.response   = Stim.gocue_onset;
        CurrPlot.stim_onset = Stim.stim_onset;
        CurrPlot.stim_dur   = Stim.stim_duration;
        
        
        stim_left_scatter_ind = Stim.ind & Neuron.left_scatter;
        stim_left_scatter_trials = Neuron.ys(stim_left_scatter_ind);
        
        CurrPlot.left_scatter.xs = Neuron.xs(stim_left_scatter_ind);
        [~,CurrPlot.left_scatter.ys] = ismember(stim_left_scatter_trials,unique(stim_left_scatter_trials));
        
        stim_right_scatter_ind = Stim.ind & Neuron.right_scatter;
        stim_right_scatter_trials = Neuron.ys(stim_right_scatter_ind);
        
        CurrPlot.right_scatter.xs = Neuron.xs(stim_right_scatter_ind);
        [~,CurrPlot.right_scatter.ys] = ismember(stim_right_scatter_trials,unique(stim_right_scatter_trials));
        
        CurrPlot.num_of_trials_right = Stim.lickR_trials;
        CurrPlot.num_of_trials_left  = Stim.lickL_trials;
        
        plot_condition(CurrPlot,cc,num_conditions);
    end
end

end



%% PSTH helper function

function [psth_xs, psth_ys] = calculate_psth(xs, bin_window, bin_width)
[bin_counts, edges]  = histcounts(xs, 'BinWidth', bin_width);
psth_xs = edges(1:end-1) + (bin_width / 2);
moving_avg_b = (1/bin_window) * ones(1,bin_window);
psth_ys = filter(moving_avg_b, 1, bin_counts)/bin_width;
end



function plot_condition(ConditionPlot,cc,num_conditions)
left_cdata = [1 0 0]; % red
right_cdata = [0 0 1]; % blue
% Calculate PSTH values
% moving average over 200 ms as per figure 1e
hist_bin_width = 0.2 / ConditionPlot.psth_bin_window;
[left_psth_xs, left_psth_ys] =...
    calculate_psth(ConditionPlot.left_scatter.xs, ConditionPlot.psth_bin_window, hist_bin_width);
[right_psth_xs, right_psth_ys] =...
    calculate_psth(ConditionPlot.right_scatter.xs, ConditionPlot.psth_bin_window, hist_bin_width);
right_scatter_offset = min(ConditionPlot.right_scatter.ys);
right_scatter_height = max(ConditionPlot.right_scatter.ys) - right_scatter_offset;
left_scatter_offset = min(ConditionPlot.left_scatter.ys);
left_scatter_height = max(ConditionPlot.left_scatter.ys) - left_scatter_offset;

left_psth_ys  = left_psth_ys/ConditionPlot.num_of_trials_left;
right_psth_ys = right_psth_ys/ConditionPlot.num_of_trials_right;

psth_height = max([left_psth_ys right_psth_ys]);
left_y_offset = left_scatter_offset;
right_y_offset = left_y_offset...
    + left_scatter_offset...
    + left_scatter_height...
    - right_scatter_offset;
subplot_height = right_y_offset...
    + right_scatter_offset...
    + right_scatter_height;



% PSTH
ax = subplot(2,num_conditions,cc+num_conditions);hold on
plot(ax, left_psth_xs, left_psth_ys, 'Color', left_cdata);
plot(ax, right_psth_xs, right_psth_ys, 'Color', right_cdata);
xline(ConditionPlot.sample,'k:');
xline(ConditionPlot.delay,'k:');
xline(ConditionPlot.response,'k:');
if strcmp(ConditionPlot.stim_type,'PT silencing') 
    fill([ConditionPlot.stim_onset ConditionPlot.stim_onset+ConditionPlot.stim_dur...
        ConditionPlot.stim_onset+ConditionPlot.stim_dur ConditionPlot.stim_onset],...
        [psth_height psth_height psth_height+1 psth_height+1],...
        'c','edgecolor','none');   
end
title(ax, ConditionPlot.label);
xlabel(ax, 'Time (Seconds)');
ylabel(ax, 'Spikes s^{-1}')
xticks(ax, [-2 0 2]);
ax.TickDir = 'out';
ax.XLim = ConditionPlot.xlim;
ax.YLim = [0 psth_height+1];
hold(ax, 'off');


% Scatter Plot
ax = subplot(2,num_conditions,cc);hold on
scatter(ax,...
    ConditionPlot.left_scatter.xs,...
    left_y_offset + ConditionPlot.left_scatter.ys,...
    'Marker', '.',...
    'CData', left_cdata,...
    'SizeData', 1);
scatter(ax,...
    ConditionPlot.right_scatter.xs,...
    right_y_offset + ConditionPlot.right_scatter.ys,...
    'Marker', '.',...
    'CData', right_cdata,...
    'SizeData', 1);
% sample, delay, response lines
xline(ConditionPlot.sample,'k:');
xline(ConditionPlot.delay,'k:');
xline(ConditionPlot.response,'k:');
if strcmp(ConditionPlot.stim_type, 'PT silencing')
    fill([ConditionPlot.stim_onset ConditionPlot.stim_onset+ConditionPlot.stim_dur...
        ConditionPlot.stim_onset+ConditionPlot.stim_dur ConditionPlot.stim_onset],...
        [subplot_height subplot_height subplot_height+2 subplot_height+2],...
        'c','edgecolor','none');   
end

title(ax, ConditionPlot.label);
xlabel(ax, 'Time (Seconds)');
ylabel(ax, 'Trials')
xticks(ax, [-2 0 2]);
ax.TickDir = 'out';
ax.XLim = ConditionPlot.xlim;
ax.YLim = [0 subplot_height+2];
hold(ax, 'off');

end

















