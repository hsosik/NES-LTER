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


************

Finally, I ran generate_attune_table_EDI.m to cut AttuneVolTable down to just the files and variables we share with EDI. 

This script has the depth values chosen for each vessel, so we will need to add values if new ships are added 



**** Model Fitting for Estimating Division Rates of Syn and Euks ******

Step 1. 

Setup_Cruise_Sliding_Model_Days_Syn_and_Euks.m 

   This script will create folders for the cruise here: 

  '\\sosiknas1\Lab_data\Attune\cruise_data\Division-rate-model\SlidingWindow\';
  '\\sosiknas1\Lab_data\Attune\cruise_data\Division-rate-model\eukaryotes\SlidingWindow\';

And within each of those folders save the input files for each day (starting at each hour) in the model fitting path. 
The only input to this function is the full path of the AttuneVolTable.mat file

This script requires input as it runs to see if light data looks ok and dawn is correctly identified. 
plot will appear with sunlightd data over time, make sure that big red dot
looks to be at approximately dawn. This can be messed up if there is
artificial peak in sunlight during the night or some other light noise at night. 
Note that dawn is only ever needed rounded to the hour, so don't worry about precision beyond that. 


Step 2. 
ApplyModel_Cruise.m 

This script will find the input files generated above and fit the model to all of them
saving a corresponding output file for each input file. 

The input to this function is the path for the folder that was generated above (either for euks or syn) 
Will have to run call the function twice to fit to syn and then euks. 

e.g. inputspath = ''\\sosiknas1\Lab_data\Attune\cruise_data\Division-rate-model\SlidingWindow\HRS2303'

Right now, the model scripts that are used to fit are those on the nas in the MVCO/FCB directory
 but if the model is ever updated, either those scripts or the path in line 12 of this function will need to change

The main roles of this function are to deal with the input file structure, check to see if its a euk or syn input, 
interpolate the light data (which for some reason we never incorporated into our input files and always do at this step)
And then call the optimization problem. 

This script will use the "FastStart" model fitting procedure, where the first parameter guess is the best fit parameters of the last day fit. 
There are still 39 other well dispersed parameter guesses used to start optimizations and the results are compared. 
But we just make sure that at least one of the initial guesses is the same as the best fit from "yesterday" (or "tomorrow" depending which order we're going through the list)


Step 3. 
Summarize_Model_Outputs.m 

This script goes through the output files generated by the model fitting procedure. 
Some user input is good for quality control. 

Requires two inputs: 
The first is the path where the model input and output files are stored, same as inputspath above
inputspath = '\\sosiknas1\Lab_data\Attune\cruise_data\Division-rate-model\SlidingWindow\HRS2303';

The other is the path where the underway attune data is and where we will save the summary table: 
outpath = '\\sosiknas1\Lab_data\Attune\cruise_data\20230429_HRS2303\bead_calibrated'; 

This script creates a file ModelOuputs.mat in the outpath directory. In this mat file are three tables

EukModelOutputs and SynModelOutputs

I would encourage future users to look at these closely for errors and to use these to change the values
in the QC_Check column to 0 if they look suspicious. 

Some columns in the table contain values for every attune underway sample within the 24 hour period. 
For example sampledates has all the samples included, and lat, lon, and temperature also contain the values
from the AttuneTable_uw_match for each of those samples. 

the lastdawn column contains the time that we have identified as the previous dawn. 
This is something that can easily be wrong and a day fit might need to be rerun (only for syn. Doesn't affect the picoeuk division rates since we don't limit division after dawn for them) 
Also, the dawnstart column equals 1 if the lastdawn and the daysstarttime are the same hour of the day
but sometimes the lastdawn changes from e.g. 12 to 11 just as the daystart changes from 11 to 12 and you "miss" a dawn in the dawnstart column. 
So some dawnsstarts may need to be chosen manually as a nearest neighbor. Maybe we could improve this in the future. 

Modelresults column stores best fit modelfits same as for MVCO model fitting. 
     modelresults - 1 by 23 vector with results of optimization process. 
     The first entry is the day in Matlab datenum form. 
     Entries 2:15 are best-fit model parameters described below: 
     modelresults(16) is the negative log likelihood. 
     modelresults(17) is the estimated division rate for the assemblage, while modelresults(18:19) are the estimated division rates for each of the two subpopulations. 
     modelresults(20:21) are the relative proportions of the two subpopulations at the end of the simulated day (as opposed to the starting proportions which is one of the parameters), 
     and modelresults(22) is the ExitFlag from createOptimProblem for that best run. 
     Modelresults(23) is the length of modelfits, or the number of model runs the optimization process tried.

gmax1=theta(1); %max fraction of cells growing into next size class, subpopn 1 
b1=theta(2);  %shape parameter for division function, subpopn 1 
E_star1=theta(3); %shape parameter of growth function (point where function switches from linear to constant), subpopn 1 
dmax1=theta(4); %max fraction of cells able to divide in a given size class, subpopn 1 
gmax2=theta(5); %max fraction of cells growing into next size class, subpopn 2 
b2=theta(6); %shape parameter for division function, subpopn 2 
E_star2=theta(7); %shape parameter of growth function (point where function switches from linear to constant), subpopn 2 
dmax2=theta(8); %max fraction of cells able to divide in a given size class, subpopn 2 
f=theta(9); %proportion parameter, specifies starting fraction of subpopn 1 
m1=theta(10); %mean volume for starting cell size distribution, subpopn 1 
m2=theta(11); %mean volume for starting cell size distribution, subpopn 2  
sigma1=theta(12); %variance parameter for starting cell size distributions for popn 1  
sigma2=theta(13); %variance parameter for starting cell size distributions for popn 2  
s=theta(14); %overdispersion parameter for the Dirichlet-multinomial distribution 



Also in ModelOuputs.mat is a merged table with Syn and Euk results, useful for plotting but with many fewer variables for quality control. 





****



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

