function meg_add_fieldtrip_paths(ft_path, ft_subdirs)
% Add fieldtrip to the Matlab path.
%
% meg_add_fieldtrip_paths(ft_path, [ft_subdirs])
%
% Purpose:
%   Add fieldtrip, including certain subdirectories, to Matlab path.
% 
% Inputs
%   ft_path: path to fieldtrip
%   ft_subdirs: cell array of fieldtrip subdirectories to add to path or
%               string to indicate a set of subdirectories
%
% Outputs
%   none
%
% Example 1:
%   ft_path = '/Volumes/server/Projects/MEG/code/fieldtrip';
%   meg_add_fieldtrip_paths(ft_path);
%
% Example 2:
%   ft_path = '/Volumes/server/Projects/MEG/code/fieldtrip';
%   ft_subdirs = {'yokogawa', 'sqdproject'};
%   meg_add_fieldtrip_paths(ft_path, ft_subdirs);
%
% Example 3:
%   ft_path = '/Volumes/server/Projects/MEG/code/fieldtrip';
%   ft_subdirs = 'yokogawa_defaults';
%   meg_add_fieldtrip_paths(ft_path, ft_subdirs);

% Author: Winawer & Chua

% add the top level path to fieldtrip (no subdirectories)
addpath(ft_path);

% if requested, add fieldtrip external subdirectories
if exist('ft_subdirs', 'var')
    % if ft_subdirs is a string, then we interpret it
    if ~iscell(ft_subdirs)
        switch ft_subdirs
            case 'yokogawa_defaults'
                addpath(fullfile(ft_path, 'external', 'yokogawa'));
                addpath(fullfile(ft_path, 'external', 'sqdproject'));
                addpath(fullfile(ft_path, 'fieldtrip_private'));
                addpath(fullfile(ft_path, 'compat'));            
            otherwise
                warning('Unknown input %s', ft_subdirs);
        end
    else
        % if it is a cell, then add path to each entry in the cell array
        for ii = 1:length(ft_subdirs)
            addpath(fullfile(ft_path, 'external', ft_subdirs{ii}));
        end
    end
end

% Run the ft_defaults which will load the necessary paths:
ft_defaults

return

