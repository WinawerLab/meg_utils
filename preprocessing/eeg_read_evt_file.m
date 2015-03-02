function ev_time = eeg_read_evt_file(pth)
% Read a NetStation events text file, extract event timing, and return a
% vector of event times in seconds
%
%  ev_time = eeg_read_evt_file(pth)
%

% test
% pth = '/Volumes/server/Projects/EEG/SSEEG/Data/Pilot_SSEEG_20150129_wl_subj001/raw/Session_20150129_1007_pt1.evt';



% Reading in the evt file results in a 1 by many characters vector
data = fileread(pth);
%  

% fileID = fopen(pth);
% formatSpec = '%s\t\t%s\t%s\t%d\t%s\t%s\t%s\t%d';
% sizeA = [8 inf];
% A = textscan(fileID,formatSpec, sizeA);
% fclose(fileID);
%% Convert a character vector of event onsets in time 
% ('00:00.000' - Hours:Minutes:Seconds.Milliseconds) 

char_length = 14; % characters that define the onset time of each DIN

% the length of the initiation sequence should be n times the number of for
% loops in the function flinitseq
length_init_seq = 12;  

% ']' (marker) precedes the time values of the triggers in DIN file 
marker = find(data==']');

% Get events from character vector
events = [];

for ii = 1:size(marker,2)
    events{ii} = data(marker(ii)+1:marker(ii)+char_length);
end

clear data; clear marker;
%% Transform these character times (HH:MM:SS.FFF) into a vector containing
% only seconds

ev_time = [];
for ii = 1:size(events,2)
    [~, ~, ~, ~, M, S] = datevec(events{ii}); % Keep only minutes and seconds
    ev_time(ii) = [M*60+S]; % Convert to only seconds
end