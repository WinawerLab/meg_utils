function root_path = meg_utils_root_dir()
% Return the path to the root meg_utils directory
%
% root_path = meg_utils_root_dir()
%
% Purpose:
%   This function must reside in the directory at the base of the meg_utils
% directory structure.  It is used to determine the location of various
% sub-directories.
% 
% Inputs
%   none
% Outputs
%   root_path = full (absolute) path to meg_utils root directory
%
% Example:
%   fullfile(meg_utils_root_dir,'preprocessing')
%
% Author: Winawer & Chua

root_path=which('meg_utils_root_dir');

root_path=fileparts(root_path);

return

