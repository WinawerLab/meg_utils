function [ev_times] = eeg_read_evt_file(ev_pth)
% Read a NetStation events text file, extract event timing, and return a
% vector of event times in seconds

% ev_time = eeg_read_evt_file(pth)

% Reading in the evt file 'fopen(pth)' results in a 1 by many characters
% vector called fileID. We use textscan, with a size that should always be
% larger than the number of DIN's in any given session (i.e. 1000000).
                                                                         
% Since the .evt file is actually made completely of strings, we put ten of
% these '%s' in the 'formatSpec' field. Ten is the number of strings in a
% line of the .evt file, at least in the way we are currently exporting. 
  
% This outputs ten cells of size (# of DIN's x 1), if we assign a variable
% to the 7th cell in A we get a vector B of all the times in strings. 
% Then we run this through a datevec loop to create an actual vector
% 'ev_time' which contains numerical time values in seconds of the recorded
% DIN's 

% example inputs:
%       ev_pth = '/Volumes/server/Projects/EEG/SSEEG/Data/Timing_test_SSEEG_20150219/Timing_test_20150219.evt';
%       ev_pth = '/Volumes/server/Projects/EEG/SSEEG/Data/Pilot_SSEEG_20150129_wl_subj001/raw/Session_20150129_1007.evt';

fileID = fopen(ev_pth);

% These event files seem to have 10 columns. So read in 10 columns of cells
formatSpec = '%s%s%s%s%s%s%s%s%s%s';
sizeA = 1000000;
A = textscan(fileID, formatSpec, sizeA);
fclose(fileID);

% The 7th columnn contains time stamps
B = A{7};

% The first three rows are header rows
B = B(4:end);

[~, ~, ~, ~, M, S] = datevec(B); % Keep only minutes and seconds
ev_times = M*60+S; % Convert to only seconds


%% Garbage
% %% Convert a character vector of event onsets in time 
% % ('00:00.000' - Hours:Minutes:Seconds.Milliseconds) 
% 
% char_length = 14; % characters that define the onset time of each DIN
% 
% % the length of the initiation sequence should be n times the number of for
% % loops in the function flinitseq
% length_init_seq = 12;  
% 
% % ']' (marker) precedes the time values of the triggers in DIN file 
% marker = find(data==']');
% 
% % Get events from character vector
% events = [];
% 
% for ii = 1:size(marker,2)
%     events{ii} = data(marker(ii)+1:marker(ii)+char_length);
% end
% 
% clear data; clear marker;
% %% Transform these character times (HH:MM:SS.FFF) into a vector containing
% % only seconds
% 
% ev_time = [];
% for ii = 1:size(events,2)
%     [~, ~, ~, ~, M, S] = datevec(events{ii}); % Keep only minutes and seconds
%     ev_time(ii) = [M*60+S]; % Convert to only seconds
% end