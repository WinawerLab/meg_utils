#! /bin/bash

# meg_to_bids.sh
# Script for processing MEG data into BIDS data
# Starts by running the mne_sqd_preproc.m script using matlab;
# then uses the meg-bids docker to run the mne_bids_export.py script.

function die {
    echo "$*"
    exit 1
}

SYNTAX="
SYNTAX: meg_to_bids.sh <subj> <sess> <raw_path> <meta_path> <out_path>

Input parameters:
  * subj: the subject ID
  * sess: the session ID
  * raw_path: the absolute or relative path to the directory containing the
    various raw *.sqd files from the MEG scanner; this should contain all of
    the *Marker*.sqd files, the *_HS.txt file, the *_Points.txt file, and the
    *_Emptyroom.sqd file (if any).
  * meta_path: the absolute or relative path to the directory containing the
    various *.tsv and *.mat files generated by the stimulus program.
  * out_path: the output BIDS project directory
"

[ "$#" -eq 5 ] || die "$SYNTAX"

subj="$1"
sess="$2"
rawp="$3"
metp="$4"
outp="$5"
time="`date '+%C%y-%m-%d %T'`"

# convert the above to absolute paths...
startpath="$PWD"
cd "$rawp"
rawp="$PWD"
cd "$startpath"
cd "$metp"
metp="$PWD"
cd "$startpath"
cd "$outp"
outp="$PWD"
cd "$startpath"

# Find this script's path (thanks to https://stackoverflow.com/questions/59895/get-the-source-directory-of-a-bash-script-from-within-the-script-itself)
SOURCE="${BASH_SOURCE[0]}"
while [ -h "$SOURCE" ]; do # resolve $SOURCE until the file is no longer a symlink
    DIR="$( cd -P "$( dirname "$SOURCE" )" >/dev/null 2>&1 && pwd )"
    SOURCE="$(readlink "$SOURCE")"
    [[ $SOURCE != /* ]] && SOURCE="$DIR/$SOURCE" # if $SOURCE was a relative symlink, we need to resolve it relative to the path where the symlink file was located
done
spth="$( cd -P "$( dirname "$SOURCE" )" >/dev/null 2>&1 && pwd )"


cat <<EOF
================================================================================
MEG to BIDS script

Author:      Noah C. Benson <nben@nyu.edu>
Date:        $time
Run Path:    '$PWD'
Script Path: '$spth'

Parameters:
  * Subject:   '$subj'
  * Session:   '$sess'
  * Raw Path:  '$rawp'
  * Meta Path: '$metp'
  * Output:    '$outp'

EOF


################################################################################
# Run the Matlab Preprocessing...

# make a temp dir
tmpp=`mktemp -d "$PWD/temp_XXXXXX" 2>/dev/null || mktemp -d -t 'mytmpdir'`
#tmpp="$rawp"
cat <<EOF
--------------------------------------------------------------------------------
Step 1: Preprocessing Data (Matlab)

  * Temporary directory: '$tmpp'
EOF
echo "  * Running Matlab script..."
echo ""
echo "________________________________________________________________________________"
matlab -nodesktop -nodisplay -nosplash <<EOF
tbUse fieldtrip;
addpath('$spth');
flnms = mne_sqd_preproc('$subj', '$sess', '$rawp', '$metp', '$tmpp');
for ii = 1:numel(flnms)
  fprintf('  * %s\n', flnms{ii});
end
EOF
echo "________________________________________________________________________________"
echo ""

# copy the points/hs txt files over
if [ "$rawp" != "$tmpp" ]
then cp "$rawp"/*.txt "$tmpp"
fi
# and the marker files
cp "$rawp"/*Marker*.sqd "$tmpp"


################################################################################
# Run the meg-to-bids docker...
cat <<EOF
--------------------------------------------------------------------------------
Step 2: Running MNE-to-BIDS script...

  * Temporary directory: '$tmpp'
EOF

echo "  * Checking Python docker... "
lns=`docker images winawerlab/meg-bids | wc -l`
#lns="0"
if [ "$lns" -lt 2 ]
then echo "  * Building Python docker... "
     curpth="$PWD"
     cd "$spth"/docker
     docker build --no-cache --tag winawerlab/meg-bids "$spth"/docker
     cd "$curpth"
fi

echo "  * Running Python script..."
echo "________________________________________________________________________________"
echo ""
docker run --rm -it -v "$tmpp":"/raw" -v "$metp":/meta -v "$outp":/out \
       winawerlab/meg-bids "$subj" "$sess" /raw /meta /out
echo "________________________________________________________________________________"
echo ""


echo "Done!"
exit 0
