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

filetype2exclude contains a list of strings that might be in filenames we don't want to include. You can always ad more to this list. 

stepsize can be adjusted to make a movie that doesn't include every single file, but now that we have test_class_assignents.m, this isn't as necessary. 

SSCDIM refers to the side scattering measurement we want to use for size calibration. It can be A or H for area or heigh. We have used Area for all the cruises in my (bethany's) 2021 processing, but we used to use H. 

Beadtype will always be "FCB". That variabile is a relic of me trying to use the Performance test beads when FCB beads weren't available. 





Once the wrapper has been run for steps 1-7 (or 1-6 with makemovieasyougo), 
run match_Attune_underway_LTER.m  to combine Attune table with underway data from REST API 

***************
***Emily can probably skip this part for now... Bethany used it as described here
Because I was interested in having size distributions for Syn and Euks, I then ran get_cruise_voldists_fromEDItable2.m to generate AttuneVolTable. 
This script also includes parsing of the underway data names to get a single measurements of latitude, longitude, salinity and other underway data of interest based on Stace's recommendations of which variabiles were best. 
Her recommendations are here: 

SAMOS database has Armstrong intermediate (not research) products in netcdf with quality-assessment:
    https://samos.coaps.fsu.edu/html/cruise_data_availability.php?data=0&ships[]=WARL&cruises[]=AR61-B&sort=overall&submit=search

Stace's underway REST API assessment of which sensor to use in ultimate LTER product (useful for Endeavor cruises):
    https://docs.google.com/spreadsheets/d/1R0vQKF7I-KzXqVeSWXhESREcIgr_WD6-VAj_F71tqWs/edit?usp=sharing

R2R post-navigation products available for EN and AR (expect about a year or so after cruise), e.g.,
    https://www.rvdata.us/search/cruise/EN649

Last is to check each cruise ship-provided data README for notes.
***************


Finally, I ran generate_attune_table_EDI.m to cut AttuneVolTable down to just the files and variables we share with EDI. 


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
-Going to revisit this after current EDI submission



