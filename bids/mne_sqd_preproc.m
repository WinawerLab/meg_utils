function filenames = mne_sqd_preproc(sub, ses, raw_path, meta_path, out_path)
% filenames = mne_sqd_preproc(sub, ses, raw_path, meta_path, out_path)
% Yields a list of output sqd filenames (in the out_path directory) that
% result from importing the SQD files in raw_path, fixing their triggers,
% then splitting them at trigger=255.
%
% In order for this function to work properly, the following must be true:
%  * The raw_path must contain all the *.sqd files associated with the
%    session.
%  * The Marker sqd files should contain the string '_Marker' somewhere in
%    the filename.
%  * Any emptyroom recordings should contain the string '_Emptyroom' in the
%    filename (note that these are not case sensitive).
%  * There must be a *_HS.txt headshape file in the raw_path.
%  * There must be a *_Points.txt fiducials file in the raw_path.
%  * The meta_path must contain all the *.mat and *.tsv files that go with
%    the experiment and document the events.
%
% By default, if the out_path is not given, then the raw_path is used. The
% output files are named based on the matlab/tsv file names.

if ~exist('out_path', 'var'), out_path = raw_path; end

if isempty(which('sqdwrite'))
    ftwd = which('ft_write_data');
    if isempty(ftwd)
        error('Could not find sqdwrite or fieldtrip');
    end
    pth = fullfile(fileparts(fileparts(ftwd)), 'external', 'sqdproject');
    addpath(pth);
end

%% Step 1
%  Go through the raw directory and collect the files we know about.
hsp_files = {};
pts_files = {};
mrk_files = {};
fsn_files = {};
dat_files = {};
d = dir(raw_path);
for ii = 1:numel(d)
    name = d(ii).name;
    nm = lower(name);
    if nm(1) == '.', continue; end
    name = fullfile(raw_path, name);
    if     endsWith(nm, '_hs.txt'),        hsp_files{end+1} = name;
    elseif endsWith(nm, '_points.txt'),    pts_files{end+1} = name;
    elseif endsWith(nm, '.fsn'),           fsn_files{end+1} = name;
    elseif exist(name, 'dir'),             continue;
    elseif endsWith(nm, '.sqd')
        if     contains(nm, '_marker'),    mrk_files{end+1} = name;
        elseif contains(nm, '_emptyroom'), continue;
        elseif startsWith(nm, 'sub-'),     continue;
        else                               dat_files{end+1} = name;
            fprintf(' * %s\n', name);
        end
    else warning(sprintf('Unrecognized file: %s', name));
    end
end

%% Step 2
%  We scan the event/matlab files in the meta directory and sort them based
%  on their initial time (t0) so that we have everything in order.
matfiles = [];
mattrigs = {};
d = dir(meta_path);
for ii = 1:numel(d)
    name = d(ii).name;
    nm = lower(name);
    if nm(1) == '.' || ~endsWith(nm, '.mat'), continue; end
    fname = fullfile(meta_path, d(ii).name);
    tsvfl = [fname(1:end-4) '.tsv'];
    if ~exist(tsvfl, 'file'), continue; end
    bdat = bidsfile_split(name);
    bdat.filename = fname;
    bdat.basename = name;
    if strcmp(bdat.sub,sub) && strcmp(bdat.ses,ses)
        load(fname, 'time0');
        bdat.t0 = time0;
        matfiles = [matfiles, bdat];
        tmp = load(fname, 'stimulus');
        tmp = tmp.stimulus.trigSeq;
        tmp = tmp(tmp > 0);
        mattrigs{end+1} = tmp(:)';
    end
end
% sort the matfiles by time 0
[dummy ind] = sort([matfiles.t0]);
matfiles = matfiles(ind);
mattrigs = mattrigs(ind);

%% Step 3
%  Walk through the data files, fixing the triggers and splitting them into
%  individual runs from the 

if numel(dat_files) > 4
    warning(['There are more than 4 raw sqd data files in the raw ' ...
             'directory; because there is no clear way to deduce the ' ...
             'ordering of these files, it may take a long time to ' ...
             'test all orderings.']);
end

trigchs = 161:168;
% We'll test the orderings in this order...
all_perms = sortrows(perms(1:numel(dat_files)));

%% Step 3
%  Loop through the permutations until we find one that fits
for pp = 1:size(all_perms, 1)
    datfls = {dat_files{all_perms(pp,:)}};
    filenames = {};
    mti = 1;
    flno = 0;
    % okay, go through each data file to fix triggers then split them
    for ii = 1:numel(datfls)
        flnm = datfls{ii};
        sqd = sqdread(flnm);
        [ts,tnos,twh] = fix_triggers(sqd(:,trigchs));
        % fix the trigger channels to be digital:
        sqd(:,trigchs) = ts;
        % okay, now we line up triggers with the file contents...
        ntrigs = numel(tnos);
        % we now need to go through however many matfiles/mattrigs are in
        % the data file...
        while mti <= numel(mattrigs) && numel(tnos) > 0
            mts  = mattrigs{mti}(:);
            nmts = numel(mts);
            % These should start and end with 256
            if (mts(1) ~= 256 && mts(1) ~= 255) || (mts(end) ~= 256 && mts(end) ~= 255)
                error('file %s: trigSeq does not start/end with 256', ...
                      matfiles(mti).filename);
            elseif nmts > ntrigs
                warning('file %s: trigSeq size mismatch', ...
                        matfiles(mti).filename);
                mti = 0;
                break;
            elseif tnos(1) < 255 || tnos(nmts) < 255
                warning('file %s: start/stop mismatch', flnm);
                mti = 0;
                break;
            elseif ~isequal(tnos(2:nmts-1), mts(2:end-1))
                warning('file %s: value mismatch', flnm);
                mti = 0;
                break;
            end
            % okay, we have a file block
            sqdblock = sqd(twh(1):twh(nmts),:);
            sqdblock(1,   trigchs) = 0;
            sqdblock(end, trigchs) = 0;
            oflnm = matfiles(mti).basename;
            oflnm = [oflnm(1:end-3) 'sqd'];
            oflnm = fullfile(out_path, oflnm);
            if exist(oflnm, 'file'), delete(oflnm); end
            sqdw = sqdwrite(flnm, oflnm, sqdblock);
            if sqdw == -1
                error(sprintf('Complete failure to write file %s', oflnm));
            else
                flno = flno + 1;
                mti  = mti  + 1;
                tnos = tnos(nmts+1:end);
                twh  = twh(nmts+1:end);
                filenames{end+1} = oflnm;
            end
        end
        if mti == 0, break; end
    end
    if mti > 0, break; end
end

function [trigger,trigno,trigwh] = fix_triggers(ts)
   % This function is by Rachel Denison and was lifted from the WinawerLab
   % meg_utils repository. It has been lightly modified for the specific
   % use-case needed in this larger script then appended.

   % exclude outliers
   md = median(ts(:));
   ts(ts > 5*md) = NaN;
   % rescale to [0 1]
   ts = ts - min(ts(:));
   ts = ts / max(ts(:));
   % check whether triggers are indicated by a low value or a high value
   if round(mean(ts(:))) == 0, trigger_is_high = true; 
   else                     trigger_is_high = false; end 
   % if triggers are indicated by a low value, then invert the signal
   if ~trigger_is_high, ts = 1 - ts; end
   % threshold to binarize from analog recording 
   ts = ts > 0.1;
   % differentiate to isolate trigger onsets, and pad to preserve length
   ts = padarray(diff(ts), [1 0], 0, 'post');
   % rectify 
   ts(ts<0) = 0;
   % mark the samples as 1 when any trigger channel signaled 
   any_trigger      = sum(ts,2) > 0;
   % list the sample numbers (time ponts) when any trigger channel signalled
   any_trigger_inds = find(any_trigger);
   % check whether 2 triggers signalled in an impossibly short time period
   % (indicating that they were probably simultanuous, but slighly misaligned)
   too_close = diff(any_trigger_inds) < 10;
   too_close = any_trigger_inds(too_close);
   %% Realign bad triggers
   %   We find any case where two trigger channel onsets were within 2
   %   samples, and we re-align them to the first sample.
   for ii=1:length(too_close)
       ts(too_close(ii),:) = sum(ts(too_close(ii)+[-9:9],:));
       ts(too_close(ii) + [-9:-1 1:9],:) = 0;
   end
   % now check that there are no time points that are too close
   any_trigger                 = sum(ts,2) > 0;
   any_trigger_inds            = find(any_trigger);
   triggers_that_are_too_close = diff(any_trigger_inds) < 10;
   assert(sum(triggers_that_are_too_close) == 0);
   % return
   trigwh  = any_trigger_inds;
   trigger = ts * 5 * md;
   trigno  = ts * 2.^(0:size(ts,2)-1)';
   trigno  = trigno(trigwh);

function parts = bidsfile_split(filename)
   spl = strsplit(filename, '.');
   s = strjoin({spl{1:end-1}}, '.');
   spl = strsplit(s, '_');
   parts = [];
   for ii = 1:numel(spl)
       ss = strsplit(spl{ii}, '-');
       parts = setfield(parts, ss{1}, ss{2});
   end