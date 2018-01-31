function [fs_lh_overlay fs_rh_overlay] = tess_bst2fs(bst_subject, fs_subject, bst_overlay, varargin)
% [fs_lh_overlay fs_rh_overlay] = tess_bst2fs(bst_subject, fs_subject, bst_overlay)
% Converts the given Brainstorm subject cortical surface overlay into a pair of [LH RH] FreeSurfer
% subject cortical surface overlays.
% 
% The overlays are all 1D vectors consisting of a measurement or datum per vertex in the respective
% hemisphere; i.e., bst_overlay must have the same number of values as there are vertices in the
% bst_subject's down-sampled mesh surface. The resulting FreeSurfer overlays will have the same
% number of values as there are vertices in the respective hemisphere.
%
% Notes:
%  - The FreeSurfer subject must be an absolute path to the subject's FreeSurfer directory or a
%    subject id that can be found in the SUBJECTS_DIR environment variable.
%  - The BrainStorm subject must be any value for which bst_get('Subject', bst_subject) returns
%    correctly.
%  - Operates on the current protocol; make sure to switch protocols if necessary.
%

% Author: Noah C. Benson <nben@nyu.edu>, 2017

% Example:
% bs_dir = '/Volumes/server/MEG/brainstorm_db/';
% fs_dir = '/Volumes/server/Freesurfer_subjects/';
% [fs_lh_overlay fs_rh_overlay] = tess_bst2fs('wl_subj010', fullfile(fs_dir,'wl_subj010'), V1sources_coherent);

% Get the BrainStorm subject:
sub = bst_get('Subject', bst_subject);
% And the protocol
proto = bst_get('ProtocolInfo');
% Find the FreeSurfer subject directory:
if ~exist(fs_subject, 'dir')
    subsdir = getenv('SUBJECTS_DIR');
    if isempty(subsdir)
        error(sprintf('Could not find given freesurfer subject: %s', fs_subject));
    end
    fs_subject = fullfile(subsdir, fs_subject);
end
if ~exist(fs_subject, 'dir')
    error(sprintf('FreeSurfer subject directory not found: %s', fs_subject));
end

% Get the MRI structure
if ~isempty(sub.Anatomy)
    sMri = bst_memory('LoadMri', sub.Anatomy(sub.iAnatomy).FileName);
else
    sMri = [];
    warn('Anatomy field of subject not found; operating without MRI structure');
end

% Read in the pial matter surface
lh_tess_file = fullfile(fs_subject, 'surf', 'lh.pial');
rh_tess_file = fullfile(fs_subject, 'surf', 'rh.pial');
[lvertices, ~] = mne_read_surface(lh_tess_file);
[rvertices, ~] = mne_read_surface(rh_tess_file);
if size(lvertices, 1) == 3 && size(lvertices, 2) ~= 3, lvertices = lvertices'; end
if size(rvertices, 1) == 3 && size(rvertices, 2) ~= 3, rvertices = rvertices'; end

fs_vertices = [lvertices; rvertices];
% These operations were lifted from io/in_tess.m (around line 137)
fs_vertices = bst_bsxfun(@plus, fs_vertices, [128 129 128] / 1000);
% These operations were lifeted from io/in_tess.m (around line 208)
fs_vertices = cs_convert(sMri, 'mri', 'scs', fs_vertices);
% Okay, the FreeSurfervertices should now be oriented correctly!

% Next step: get the brainstorm pial surface
% First, find the white matter surface:
pial_surf_id = -1;
if nargin > 4
    % pial surface was provided
    wm_surf_name = varargin{5};
    for i = 1:numel(sub.Surface)
        if strcmp(sub.Surface(i).Comment, wm_surf_name)
            pial_surf_id = i;
            break;
        end
    end
else
    % we look for the surface named white_*V where the * is the lowest number of all surfs
    n = 10^10;
    for i = 1:numel(sub.Surface)
        if ~strcmp(sub.Surface(i).SurfaceType, 'Cortex'), continue; end
        nn = sscanf(sub.Surface(i).Comment, 'cortex_%dV');
        if isempty(nn), continue; end
        if nn < n
            n = nn;
            pial_surf_id = i;
        end
    end
end
if pial_surf_id < 0, error('No pial surface found for subject!'); end
% Next, get the surface data itself:
anat_dir = proto.SUBJECTS;
pial_path = fullfile(anat_dir, sub.Surface(pial_surf_id).FileName);
pial_surf = load(pial_path);
bst_vertices = pial_surf.Vertices;

% Next step: find the mapping between BrainStorm vertices and FreeSurfer vertices
idcs = dsearchn(bst_vertices, fs_vertices);
fs_overlay = bst_overlay(idcs);
fs_lh_overlay = fs_overlay(1:size(lvertices, 1));
fs_rh_overlay = fs_overlay((size(lvertices, 1) + 1):end);
% That's all!