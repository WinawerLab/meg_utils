function f = meg_check_triggers(fileName, trigChans, verbose)
% function meg_check_triggers(fileName, [trigChans], [verbose])

% INPUTS:
% fileName is the path to the sqd file
% trigChans (optional) is the trigger channels to plot. Default is 160:167.
% verbose (optional) [1 or 0] specify whether to print the trigger count 
%     for each channel
%
% OUTPUTS:
% f is a figure handle for the generated plot
%
% Rachel Denison
% 2014

%% deal with inputs
if nargin<2 || isempty(trigChans)
    trigChans = 160:167; % all trigger channels
end
if nargin<3 || isempty(verbose)
    verbose = 1;
end

%% plot triggers from each channel
f = figure;
hold all

for iTrig = 1:numel(trigChans)
    try
        triggers = all_trigger(fileName, trigChans(iTrig));
        trigTimes = triggers(:,1);
        
        plot(trigTimes, ones(size(trigTimes))+iTrig-1, '.')
        
        if verbose
            fprintf('Channel %d: %d triggers\n', trigChans(iTrig), numel(trigTimes))
        end
    catch err
        if strcmp(err.identifier, 'MATLAB:badsubscript')
            if verbose
                fprintf('Channel %d: 0 triggers\n', trigChans(iTrig))
            end
        end
    end 
end
ylim([0 10])
xlabel('time')
ylabel('triggers')
