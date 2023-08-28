README Attune Cruise Underway Processing

The main Script for processing a cruise is process_wrapper_2021.m
This function takes as input a basepath, but also needs to be adjusted in the first section based on the attune settings for that cruise and what you would like to do. 

basepath is usually something like this: 
basepath = '\\sosiknas1\Lab_data\Attune\cruise_data\20201013_EN657\'

Make sure that if you are starting a new cruise, Step 1 = 1. In theory you should be able to runs step 1-7 all at once, by setting them all to 1, 
but frequently you have to redo parts, so they can all be turned on or off. You just need to make sure the previous steps have been done already for any new steps you are going to run. 


In order to perform class assignments or test_class_assignments.m, step2 of the wrapper, you need a classifier function. 
We have been making cruise-specific classifier functions such as "assign_class_EN657"

It may make sense to start with a copy of a classifier from a different cruise with similar settings or at a similar time of year. 
These functions take fcs data as input and output class assingments like Syn, Euk, etc. 
They start by drawing first approximation boxes around clusters "main_gates" 
Then those edges are adjusted based on the distribution of particles within them. 
Sometimes I have different thresholds for different parts of the cruise, and I use the variable, "phase = 1,2 etc." to define those sections. 

I like using the script, test_class_assignments.m, to adjust the classifier for the whole cruise before creating and editing any of the class files. 

Other things to pay attention to in the wrapper script are the "OD2setting" and the "moviechannels" which will depend on the cruise settings as recorded in 
Settings_configuration_history.xlxs
moviechannels can be set to "early" for all cruises that we did not have scattering measured on GL1 (TN368 and prior). 
OD2setting is 'SSC' or 'GL1'. 
Since our settings have been relatively steady for a while, it typically should be: moviechannels = 'late'; OD2setting = 'GL1'. 

flags like "appendonly", "dont_overwrite_volumes", and "makemovieasyougo" can be adjusted based on whether you are redoing any of the steps for a cruise you have already processed or partially processed.

filetype2exclude contains a list of strings that might be in filenames we don't want to include. You can always add more to this list. 

stepsize can be adjusted to make a movie that doesn't include every single file, but now that we have test_class_assignents.m, this isn't as necessary. 

SSCDIM refers to the side scattering measurement we want to use for size calibration. It can be A or H for area or heigh. We have used Area for all the cruises in my (bethany's) 2021 processing, but we used to use H. 

Beadtype will always be "FCB". That variabile is a relic of me trying to use the Performance test beads when FCB beads weren't available. 


*** 

Once the wrapper has been run for steps 1-7 (or 1-6 with makemovieasyougo), all particles have been assigned a class and volume. 

Step 8 is match attune samples to underway environmental meansurements.
For new cruises, and new ships especially, you should check that underway data is being imported correctly into matlab before proceeding. 

Usually the restapi address works well as the input: 
e.g.  uw_fullname = 'https://nes-lter-data.whoi.edu/api/underway/en608.csv'
But sometimes you need to specificy the path to a local csv 
e.g.  uw_fullname = '\\sosiknas1\Lab_data\LTER\20200201_EN649\scs\proc\cruise\Data60Sec_Cruise_20200201-004500.csv';


Step 9 creates AttuneVolTable.mat 
This step discretizes the particles by volumes into the volume bins used for model fitting and division rate estiamtion. 
This is also the function where we parse the underway data names to get single standaradized columns for measurements of latitude, longitude, salinity and other underway data of interest based on Stace's recommendations of which variabiles were best. 

Her recommendations are here: 

SAMOS database has Armstrong intermediate (not research) products in netcdf with quality-assessment:
    https://samos.coaps.fsu.edu/html/cruise_data_availability.php?data=0&ships[]=WARL&cruises[]=AR61-B&sort=overall&submit=search

Stace's underway REST API assessment of which sensor to use in ultimate LTER product (useful for Endeavor cruises):
    https://docs.google.com/spreadsheets/d/1R0vQKF7I-KzXqVeSWXhESREcIgr_WD6-VAj_F71tqWs/edit?usp=sharing

R2R post-navigation products available for EN and AR (expect about a year or so after cruise), e.g.,
    https://www.rvdata.us/search/cruise/EN649

Last is to check each cruise ship-provided data README for notes.


When new ships are used for future cruises, we may need to add cases to these parsing steps. For example, HRS2303 has some weird variable names. 


***************


Finally, I ran generate_attune_table_EDI.m to cut AttuneVolTable down to just the files and variables we share with EDI. 

This script has the depth values chosen for each vessel, so we will need to add values if new ships are added 


Emily's notes on steps 8-9 August 2023:

Step 8 in the process wrapper 

With a new cruise,
load in the csv and see if you have all the right columns and whether they have names
The two things that can happen are you have maybe 2 columns smooshed into one column, wrong delimiters
Or, for some reason no headers, or headers ended up in first row of data.

Once you figure out what commands you need in matlab to read the table correctly, you can put them in another else if in match_Attune_underway_LTER.m

Some of them might be redundant, but just made a new case instead of exploring.

Goal is to have any of the measurements int he table imported into matlab and to have the dates in a readable format. 

Date function in matlab now seems pretty good and you can decide what format you want things in etc, 
Our BEthany Attune script currently does both.

For cruises that don't have a regular API read it saves an underway.mat in the Attune data area for that cruise.

In Step 9, 
get_cruise_voldists_fromEDItable2.m
Starting at line 54 parsing names from tables.
Stace gave input on which variables to use, for example furuno...


-----------------------------------------



Questions from Emily 7 Dec 2022:

-What to look for on beads steps. What might tip me off about a problem?

-How to do phases of the gating if needed:

Example for different phases in EN655 simpler, TN368 big mess example. The phase is only within the asign class fiile.
Basically put time with matlab dates at the top, and then phases around sections of the thresholds.
-> the assign_class_EN687 that I am currently using is one that floats the bottom of the Syn, not all of them do. 


-what else do you look for when making the movies? I was mostly looking at Syn. The green Euks semed fine. The other Euks I didn’t pay much attention, yet. Maybe need a different plot?

-Diagonal line – what to look for?

-How to make the movies

-Keeping track of what’s done

-Final steps after process wrapper? Should I do all of them. Need the EDI table!
-Going to revisit this after current EDI submission - this was resolved in steps 8-9 added by Bethany


---------------------------------


Notes from Summer 2023 on Grazer processing special case:

Values that we use in calibration are in the individual class files

line 180 ish in use calibration stats linear, we actually
use the bead value and could adjust manually there. 

use calibration looks at the hv setting in the datafile, and uses bead
value from FCB bead runs with same hv setting. Looks like it
is already working for Ron Brown.

Skip process wrapper step three for AR29

scatter_value comes from what SSC channel and dimension,

Take home messages:
use calibration stats linear already acounts for if there was
a change in hv settings, AND we added a line around line 182
for if you just want to set the bead value to a specific
value. 

Problem with RB1904 was that the settings in step 2 are 
also read in the other steps and maybe I had not left
OD2 setting as 'none'

