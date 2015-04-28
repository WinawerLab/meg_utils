function [trl,Events] = mytrialfun_all(cfg, threshold, trigNum)
%MYTRIALFUN_ALL detects down-flank triggers from an MEG sqd file with
% specified trigger threshold and defines epochs of data that could be 
% used for further processing and analysis using Fieldtrip functions, 
% i.e., data that will be read in by FT_PREPROCESSING. 
%
%    INPUT variables:
%          cfg: a structure-array with the following fields to be specified
%                 cfg.dataset           - original .sqd file holding the  
%                                         trigger and MEG data
%                 cfg.trialdef.trig     - trigger channels used to encode
%                                         the trigger 'codes'( the 
%                                         numbering convention is that the 
%                                         first MEG channel is 0, not 1).
%                 cfg.trialdef.prestim  - pre-trigger time period to 
%                                         include in  the epoch (in 
%                                         seconds) 
%                 cfg.trialdef.poststim - post-trigger time period to
%                                         include in the epoch (in seconds) 
%    threshold: threshold for flank detection
%    trigNum  : total number of triggers actually used. The default is the 
%               number of triggers detected in the function, though one 
%               should specify this for error checking.
%               
%               
%
%    OUTPUT variables:
%          trl:    a Nx3 matrix (N is the number of trials) with trial 
%                  definition [begin end offset]. The first column contains 
%                  the sample-indices  of the begin of each trial relative 
%                  to the begin of the raw data, the second column contains
%                  the sample-indices of the end of each trial, and the 
%                  third column contains the offset of the trigger with
%                  respect to timepoint 0 of that trial.
% Events.trigger:  sample-indices of the trigger events relative to the 
%                  begin of the raw data for the defined trials. The
%                  trigger events correspond with the first sample of the 
%                  down-going flank.
% Events.channel:  trigger channels (as string vectors) correspond to each 
%                  trigger event.
%
%     by Sirui Liu

hdr   = ft_read_header(cfg.dataset);

trigger = ft_read_event(cfg.dataset,'trigindx',cfg.trialdef.trig,...
    'threshold',threshold,'detectflank','down');

% error check the triggers
if ~exist('trigNum','var')
    trigNum = [];
end
if isempty(trigNum)
    trigNum = length(trigger);
end

if length(trigger) ~= trigNum,
    warning(sprintf('wrong number of triggers, only found %d trigger\n',...
        length(trigger)));
end


pretrig  = - cfg.trialdef.prestim  * hdr.Fs;
posttrig =   cfg.trialdef.poststim * hdr.Fs;


for j = 1:length(trigger);  
      
    trlbegin = trigger(j).sample +  pretrig;     
    trlend   = trigger(j).sample + posttrig;       
%               begin       end   prestim 
    trl(j,:) = [trlbegin trlend pretrig];
    
    Events(j).trigger = trigger(j).sample;
    Events(j).channel = trigger(j).type;
%   Events(j).value   = trigger(j).value;
end


end
