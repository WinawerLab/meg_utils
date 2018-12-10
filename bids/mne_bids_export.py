# Some standard libraries:
import os, sys, six, h5py, time, warnings, pandas
import numpy       as np
import scipy       as sp

import mne, mne.io, mne_bids

def mne_read_raw_kit(sqd_files, mrk_files, elp_file, hsp_file,
                     subject=None, hsp_points=5000,
                     stim=[160,161,162,163,164,165,166,167],
                     stim_code='channel', slope='-',
                     preload=False, verbose=False):
    '''
    mne_read_raw_kit(sqf_files, mrk_files, elp_file, hsp_file) yields a list of mne
      Raw data objects corresponding to the given list of sqd files. The data are
      not pre-loaded. If a string is passed instead of a list for the filename, then
      only the single Raw object is returned.
    '''
    if not pimms.is_vector(sqd_files):
        tmp = mne_read_raw_kit([sqd_files], mrk_files, elp_file, hsp_file,
                               stim=stim, stim_code=stim_code, slope=slope,
                               preload=preload, verbose=verbose)
        return tmp[0]
    # We need to import the elp and hsp file points
    dat = []
    for flnm in [hsp_file, elp_file]:
        rows = []
        with open(flnm, 'r') as f:
            while True:
                r = f.readline().strip()
                if not r or r[0] == '\x00': break
                elif r.startswith('%'): continue
                else: rows.append(list(map(float, r.split())))
        dat.append(np.array(rows))
    (hsp_data,elp_data) = dat
    # Downsample hsp_points if necessary
    if hsp_points is not None:
        ii = np.random.choice(list(range(len(hsp_data))), hsp_points, replace=False)
        hsp_data = hsp_data[ii]
    raws = [
        mne.io.read_raw_kit(
            flnm,
            mrk=mrk_files,
            elp=elp_data,
            hsp=hsp_data,
            stim=stim,
            stim_code=stim_code,
            slope=slope,
            verbose=verbose,
            allow_unknown_format=True)
        for flnm in sqd_files]
    if subject:
        for raw in raws:
            raw.info['subject_id'] = subject
    # Annotate the trigger channels correctly
    chfix = {('MISC %03d'%chi):'stim' for chi in range(1,9)}
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        for raw in raws:
            raw.info['line_freq'] = 60;
            raw.set_channel_types(chfix)
    return raws

syntax_msg = '''
SYNTAX: python mne_bids_export.py <subject-id> <session-id> <raw-path> <meta-path> <bids-path>

The subject-id and session-id must be those used with BIDS.
The raw-path must be the directory containing all raw *.sqd files as well as the *_Marker?_*
  files and the *_Points.txt and *_HS.txt files.
The meta-path must contain all the *.tsv (events) and *.mat (stimuli) files.
The bids-path must be the output directory to which the BIDS data should be saved.
'''

if len(sys.argv) != 5:
    sys.stderr.write(syntax_msg)
    sys.exit(1)

(sub, ses, raws_path, meta_path, bids_path) = sys.argv

if not os.path.isdir(raws_path):
    sys.stderr.write('Could not find raw directory: %s' % raws_path)
    sys.exit(2)
if not os.path.isdir(meta_path):
    sys.stderr.write('Could not find meta directory: %s' % meta_path)
    sys.exit(2)

sqd_files = [os.path.join(raws_path, x)
             for x in os.listdir(raws_path)
             if x.startswith('sub-') and x.endswith('.sqd')]
mrk_files = [os.path.join(raws_path, x)
             for x in os.listdir(raws_path)
             if '_marker' in x.lower() and x.endswith('.sqd')]
elp_file = next((os.path.join(raws_path, x)
                 for x in os.listdir(raws_path)
                 if x.lower().endswith('_points.txt')),
                None)
hsp_file = next((os.path.join(raws_path, x)
                 for x in os.listdir(raws_path)
                 if x.lower().endswith('_hs.txt')),
                None)
if elp_file is None:
    sys.stderr.write('No _Points.txt file found in raw directory')
    sys.exit(3)
if hsp_file is None:
    sys.stderr.write('No _HS.txt file found in raw directory')
    sys.exit(3)

stim_dir = os.path.join(bids_path, 'stimuli')

raws = mne_read_raw_kit(sqd_files, mrk_files[0], elp_file, hsp_file, subject=sub)

for ii in range(len(raws)):
    #bname = mne_bids.make_bids_basename(subject='wlsubj051', session='NYUMEG01', task='bair')
    bname = os.path.split(raws[ii].filenames[0])[1];
    bname = '.'.join(bname.split('.')[:-1])
    bdata = {kk[0]:kk[1] for pp in bname.split('_') for kk in [pp.split('-')] if len(kk) == 2}
    mne_bids.write_raw_bids(raws[ii], bname, bids_path,
                            events_data=None, overwrite=True, verbose=False)
    # since mne_bids doesn't correctly write out our events data...
    pnm = mne_bids.make_bids_folders(subject=bdata['sub'], kind='meg', session=bdata['ses'],
                                     make_dir=False, output_path=bids_path)
    edata_in  = os.path.join(meta_path, bname + '.tsv')
    edata_out = os.path.join(pnm, bname + '_events.tsv')
    shutil.copyfile(edata_in, edata_out)
    # we need to make the stimulus directory also...

    if not os.path.isdir(stim_dir): os.makedirs(stim_dir)
    mdata_in = os.path.join(mata_path, bname + '.mat')
    mdata_out = os.path.join(stim_dir, bname + '.mat')
    shutil.copyfile(mdata_in, mdata_out)

# Done!
sys.exit(0)
