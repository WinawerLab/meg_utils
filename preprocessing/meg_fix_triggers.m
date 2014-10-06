function trigger = meg_fix_triggers(ts)
% Convert an analog signal representing triggers in binary, and return a
% vector of integer trigger values over time.
%
%   trigger = meg_fix_triggers(ts)
%
% Inputs:
%   ts: matrix of trigger values (time points by channel). Each column is
%       an analog representastion of a binary value (0 or 1). We thus convert
%       the columns to binary and then combine the values in a row from
%       binary to a base 10 integer. 
%
% Output: vector of base 10 trigger values over time (time points x 1)
%
% Example:
%       ts = [0 0 0; 0 0 0; 1 0 0; 0 0 0; 1 1 1];
%       trigger = meg_fix_triggers(ts)
%       This should return trigger = [0 0 1 0 7];

error('Not yet implemented')

% time
num_time_points = size(ts,1);

% t = 0:size(ts,1)-1; % samples

% number of trigger channels
num_channels = size(ts,2);

% rescale to [0 1]
ts = ts - min(ts(:));
ts = ts / max(ts(:));

% check whether triggers are indicated by a low value or a high value
if round(mean(ts)) == 0, trigger_is_high = true; 
else                     trigger_is_high = false; end 

% if triggers are indicated by a low value, then invert the signal
if ~trigger_is_high, ts = 1 - ts; end

% threshold, diff, and pad
ts = ts > 0.1;
ts = padarray(diff(ts), [1 0], 0, 'post');

% store the trigger time points here
trigger_onsets = zeros(num_time_points, num_channels);

% find trigger onsets in each trigger channel
for ii = 1:num_channels
    [~, inds] = findpeaks([0; ts(:,ii)],'MINPEAKDISTANCE', 20);
    trigger_onsets(inds,ii) = 1;
end

% check whether there are any cases of two triggers impossibly close
trigger_present      = sum(trigger_onsets,2) > 0;

trigger_present_inds = find(trigger_present);

trigger_problems     = [0; diff(trigger_present_inds) < 10];

trigger_problem_inds =find(trigger_problems);

bad_time_points = trigger_present_inds(trigger_problem_inds);

% Issue: most of the above indexes return all zeros, check on that 
% Solution: because of padding, the triggers are present -1 or -2 time
% points before the the time points listed in bad_time_points

%% go on from here (new general fix triggers algorithm WIP)


for ii=1:length(bad_time_points)
    ts(bad_time_points(ii),:) = ts(bad_time_points(ii)-1,:) + ts(bad_time_points(ii)-2,:);
    ts(bad_time_points(ii)-1,:) = zeros([1,size(ts,2)]);
    ts(bad_time_points(ii)-2,:) = zeros([1,size(ts,2)]);
end

%%
trigger = trigger_onsets * [1 2 4 8]';




% % Plot all the triggers
% figure(200); clf
% set(gca, 'ColorOrder', jet(4), 'Color', [.6 .6 .6]); hold all;
% plot(bsxfun(@plus, ts(:,161:164), (1:4)*40000), '-o',  'LineWidth', 2);
% legend(cellstr(num2str((161:164)')))
% set(gca, 'XTick', weird, 'XGrid', 'on');
% title('Trigger channels (indices and values)');
% xlabel('Time [s]');
% ylabel('Relative amplitude');
% legend('Trigger 1 / Chan 161', 'Trigger 2 / Chan 162', 'Trigger 3 / Chan 163', 'Trigger 4 / Chan 164');


%% Fix triggers when they have value 11 or higher
%

% Look whether this still holds!
for ii = 1:length(weird)
    
    % Get index of three triggers
    a = find(trig1.ind == weird(ii));
    b = find(trig2.ind == weird(ii));
    c = find(trig3.ind == weird(ii));
    c = find(trig4.ind == weird(ii));
    
    % If there is no element for trigger 1 and trigger 2, it means that
    % trigger 3 is 1 ms early. So we add +1
    if isempty(a) && isempty(b)
        trig3.ind(c) = trig3.ind(c) + 1;
        
        % If trigger 2 and trigger 3 have contain both elements, it means that
        % trigger 1 is 1 ms late. So we discard 1 ms of trigger 1.
    elseif ~isempty(b) && ~isempty(c)
        trig1.ind(find(trig1.ind==weird(ii)+1)) = trig1.ind(find(trig1.ind==weird(ii)+1)) - 1;
        
    end
end

%% Recompute triggers to check weird peaks again

trigger_onsets = zeros(length(t), 4);
trigger_onsets(trig1.ind,1) = 1;
trigger_onsets(trig2.ind,2) = 1;
trigger_onsets(trig3.ind,3) = 1;
trigger_onsets(trig4.ind,4) = 1;

trigger = trigger_onsets * [1 2 4 8]'; clear triggerBinary;

%% check that we have the right number of triggers
fprintf('%s\n','Look at sum of trigger to check misalignment')
fprintf('%s\t\t%s\n', 'Trigger nr', 'Sum')
for ii = 1:16; fprintf('%d\t%d\n', ii, sum(trigger == ii)); end

%% Fix other weird triggers (i.e. too much of a certain trigger)

trigger_times = find(trigger);
trigger_intervals = diff(trigger_times);
[sorted, inds] = sort(trigger_intervals);

weird_triggers = find(sorted==1);
nr_weird = length(weird_triggers);

weird2 = zeros(nr_weird,1);

for ii = 1:nr_weird
    weird2(ii) = trigger_times(inds(ii));
end

%% We know this code is ugly, but for now it's generic enough to handle all the weird trigger cases. (For the first data set).
for ii = 1:length(weird2);
    a = find(trig1.ind == weird2(ii));
    b = find(trig2.ind == weird2(ii));
    c = find(trig3.ind == weird2(ii));
    d = find(trig4.ind == weird2(ii));
    
    
    
    % If there is no element for trigger 1 and trigger 2, it means that
    % trigger 3 is 1 ms early. So we add +1
    if isempty(a) && isempty(d)
        
        timepoint = trig2.ind(b);
        if isempty(timepoint)
            timepoint = trig3.ind(c);
        end
        
        if ~isempty(find(trig1.ind == timepoint+1))
            nr_of_trigger_misaligned = (find(trig1.ind == timepoint+1));
            trig1.ind(nr_of_trigger_misaligned) = trig1.ind(nr_of_trigger_misaligned)-1;
            
            % If there are two late (or one early), then fix the second
            % late trigger
            if ~isempty(find(trig2.ind == timepoint+1))
                nr_of_trigger_misaligned = (find(trig2.ind == timepoint+1));
                trig2.ind(nr_of_trigger_misaligned) = trig2.ind(nr_of_trigger_misaligned)-1;
            end   
            
        elseif ~isempty(find(trig4.ind == timepoint+1))
            nr_of_trigger_misaligned = (find(trig4.ind == timepoint+1));
            trig4.ind(nr_of_trigger_misaligned) = trig4.ind(nr_of_trigger_misaligned)-1;
        end
        
        
    elseif isempty(b) && isempty(d)
        timepoint = trig1.ind(a);
        if isempty(timepoint)
            timepoint = trig4.ind(d);
        end
        
        if ~isempty(find(trig2.ind == timepoint+1))
            nr_of_trigger_misaligned = (find(trig2.ind == timepoint+1));
            trig2.ind(nr_of_trigger_misaligned) = trig2.ind(nr_of_trigger_misaligned)-1;
            
        elseif ~isempty(find(trig4.ind == timepoint+1))
            nr_of_trigger_misaligned = (find(trig4.ind == timepoint+1));
            trig4.ind(nr_of_trigger_misaligned) = trig4.ind(nr_of_trigger_misaligned)-1;
        end
        
    elseif isempty(a) && isempty(c)
        
        timepoint = trig4.ind(d);

        if ~isempty(find(trig1.ind == timepoint+1)) % 1
            nr_of_trigger_misaligned = (find(trig1.ind == timepoint+1));
            trig1.ind(nr_of_trigger_misaligned) = trig1.ind(nr_of_trigger_misaligned)-1;
        elseif ~isempty(find(trig3.ind == timepoint+1)) % 3
            nr_of_trigger_misaligned = (find(trig3.ind == timepoint+1));
            trig3.ind(nr_of_trigger_misaligned) = trig3.ind(nr_of_trigger_misaligned)-1;
        elseif ~isempty(find(trig2.ind == timepoint+1)) % 2
            nr_of_trigger_misaligned = (find(trig2.ind == timepoint+1));
            trig2.ind(nr_of_trigger_misaligned) = trig2.ind(nr_of_trigger_misaligned)-1;
            
            
        end
    end
    
end


%% Recompute triggers to check weird peaks again

trigger_onsets = zeros(length(t), 4);
trigger_onsets(trig1.ind,1) = 1;
trigger_onsets(trig2.ind,2) = 1;
trigger_onsets(trig3.ind,3) = 1;
trigger_onsets(trig4.ind,4) = 1;

trigger = trigger_onsets * [1 2 4 8]'; clear triggerBinary;

%% check that we have the right number of triggers
fprintf('%s\n','Look at sum of trigger to check misalignment')
fprintf('%s\t\t%s\n', 'Trigger nr', 'Sum')
for ii = 1:16; fprintf('%d\t%d\n', ii, sum(trigger == ii)); end

%% Do two double checks

% Check whether all conditions are shown 150 times.
for ii = 1:9
    assert(sum(trigger == ii) == 150);
end

% Check whether the sum of all the conditions, plus 30, is equal to the
% blanks
sum_of_condition = [];
for ii = 1:9
    sum_of_condition(ii) = sum(trigger == ii);
end

assert(sum(sum_of_condition)+30 == sum(trigger==10));









