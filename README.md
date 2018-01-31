meg_utils
=========
This is a set of Matlab tools for analyzing MEG or EEG data. 
The data we work with is acquired at NYU's MEG lab http://psych.nyu.edu/nellab/meglab.html. 
Or at the NYU's EEG lab usinga 128 channel EGI cap & NetStation software.

In this repository you can find the following folders:

* experiments
	* gamma_meg
		* experimentTools:  Scripts and their subfunctions to create a visual gamma experiment
		* HPC:				Scripts and functions to run analyses on the High Performance Cluster (HPC)
		* Scripts:			
			* Analysis: 	Scripts and their subfunctions to run a full analysis trees
			* Visualization: Scripts and functions to run full visualization analyses
	* sseeg
	* ssmeg
	
* external
	* field_trip:			Functions that depend on or support fieldtrip analyses
	* Other functions written by collaborators
	
* ft_development: 			Functions that support fieldtrip analysis

* preprocessing:			General functions that are used in preprocessing of MEG or EEG data

* tutorials:				Scripts that show example analysis

* utils:					Functions that make our like easier

* Visualization:			General functions and scripts specifically for visualizing MEG or EEG data