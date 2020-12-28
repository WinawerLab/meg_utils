function bst_overlay = tess_fs2bst(bst_subject, fs_subject, fs_lh_overlay, fs_rh_overlay, varargin)
% bst_overlay = tess_fsb2bst(bst_subject, fs_subject, fs_lh_overlay, fs_rh_overlay)
% Converts the given FreeSurfer subject left and right hemisphere cortical surface overlays into a
% BrainStorm overlay for the pial surface of the given BrainStorm subject with the fewest vertices.
% 
% The overlays are all 1D vectors or 2D matrices consisting of a measurement or datum per vertex in
% the respective hemisphere; i.e., fs_lh_overlay must have the same number of rows as there are
% vertices in the fs_subject's left hemisphere surface. The result, bst_overlay, will have the same
% number of rows as there are vertices in the down-sampled pial surface of bst_subject. The left and
% right hemisphere overlays must have the same number of columns.
%
% Notes:
%  - The FreeSurfer subject must be an absolute path to the subject's FreeSurfer directory or a
%    subject id that can be found in the SUBJECTS_DIR environment variable.
%  - The BrainStorm subject must be any value for which bst_get('Subject', bst_subject) returns
%    correctly.
%  - Operates on the current protocol; make sure to switch protocols if necessary.
%
% Example:
%  %% (requires that the FreeSurfer toolbox is on the Matlab path)
%  >> fs_lh_curv = read_curv('/my/subjects/dir/bert/surf/lh.curv');
%  >> fs_rh_curv = read_curv('/my/subjects/dir/bert/surf/rh.curv');
%  %% (assumes that BrainStorm subject 'bert_bst' was imported from FreeSurfer's 'bert' subject)
%  >> bst_curv = tess_fs2bst('bert_bst', 'bert', fs_lh_curv, fs_rh_curv);

% Author: Noah C. Benson <nben@nyu.edu>, 2017


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

% Read in the pial matter surface (since Brainstorm uses pial surface for
% computing the forward model, and the vertices of the downsampled white, mid gray and
% pial surfaces are not identical)
lh_tess_file = fullfile(fs_subject, 'surf', 'lh.pial');
rh_tess_file = fullfile(fs_subject, 'surf', 'rh.pial');
[lvertices, ~] = mne_read_surface(lh_tess_file);
[rvertices, ~] = mne_read_surface(rh_tess_file);
if size(lvertices, 1) == 3 && size(lvertices, 2) ~= 3, lvertices = lvertices'; end
if size(rvertices, 1) == 3 && size(rvertices, 2) ~= 3, rvertices = rvertices'; end
if size(lvertices, 1) ~= size(fs_lh_overlay, 1)
    error('left hemisphere overlay has different number of vertices than surface');
elseif size(rvertices, 1) ~= size(fs_rh_overlay, 1)
    error('right hemisphere overlay has different number of vertices than surface');
end
fs_vertices = [lvertices; rvertices];
% These operations were lifted from io/in_tess.m (around line 137)
fs_vertices = bst_bsxfun(@plus, fs_vertices, [128 129 128] / 1000);
% These operations were lifeted from io/in_tess.m (around line 208)
fs_vertices = cs_convert(sMri, 'mri', 'scs', fs_vertices);
% Okay, the FreeSurfervertices should now be oriented correctly!

% Next step: get the brainstorm pial surface
% First, find the pial matter surface:
pial_surf_id = -1;
if nargin > 4
    % if other surface was provided
    alternative_surf_name = varargin{1};
    for i = 1:numel(sub.Surface)
        fn = strsplit(sub.Surface(i).FileName, '/'); % check file names
        if strcmp(fn{2}, [alternative_surf_name '.mat'])
            pial_surf_id = i;
            break;
        end
    end
else
    % we look for the surface named cortex_*V where the * is the lowest number of all surfs
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
idcs = dsearchn(fs_vertices, bst_vertices);
fs_overlay = [fs_lh_overlay; fs_rh_overlay];
bst_overlay = fs_overlay(idcs,:);
