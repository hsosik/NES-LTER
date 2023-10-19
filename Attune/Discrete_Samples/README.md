# Attune_Cruise_Process
Scripts for processing Sosik Lab attune flow cytometry data from discrete niskin samples

Temporarilly, scripts are replicated in two repos: 
https://github.com/beefowler/Attune_Cruise_Process
   and
https://github.com/hsosik/NES-LTER/tree/master/Attune/Discrete_sample_processing



Below is the protocol for processing discrete samples. 

# Attune Discrete Sample Protocol 

1. The main script for processing a cruise of discrete samples is "Process_Preserved_AllSteps.m"

To run this script, you will need to adjust the top for a specific cruise: 
    basepath = where to look for preserved fcm sample fcs fles
    cruisename 
    restpath = where to look for CTD data as metadata for niskin files 
    elogpath = where to look for metadata associated with underway discrete files if there are any

    you can also adjust the flags Step1 and Step5only if you don't need to redo file set up (step 1) or if you only want to redo size calibration (step 5)
    if you are processing a cruise from the beginning, set Step1 ==1 and Step5only == 0. 

    hierarchical_gates can be set to either 'True' or 'False'. This has to do with whether gates were drawn within other gates. 
        early cruises did not have hierarchical gating, but OTZ cruises that Alexi gated do, and subsequent cruises usually will have hierarchical gates. 
        As of 6/29/23, code will only look for a single parent gate. No double parents. Grandparents are fine. 

Then run the script. Fingers crossed it goes smoothly through all the steps. 

This script generates: 
    - an outputs folder within the preserved folder where the outputs will be
    - Processing_variables.mat which saves the paths and hierarchical gate setting used for processing
    - FCSList.mat which has all the fcs files, identified cast and niskins, as well as dates processed, trigger channels and hv, and "Vol_analyzed_ml" which is 80% of the total sample volume
        80% is because this is the proportion of the sample used in quantification of the cells. 
    - Gated_Table.mat which again has a row for every fcs file, this time with gate names and logic for every gate in the corresponding aws files.
        gate_assignments has a 1 or 0 for every particle within every gate within the aws file. This is what will be used to generate final class assingments for each particle. 
        polygon columns provide data relating to the actual polygons used for gating. For example, polygon_names might include an entry "syn", which has a corresponding cell in 
            the "polygon_vars" (variables) entry which indicates the channels used for drawing the polygon ("SSC-H" and "GL2-H"), and also a corresponding cell within the 
            entry for "polygon_vals" (values) with the coordinates of the X number of points used to draw this polygon (e.g. [3e3, 1e6], [3e4, 1e6], ...) 
      After Step3 metadata from REST API is added to the gated_table so you'll have lat, lon, temperature, salinity etc. 
    - A folder within outputs labelled "class" which  has class files for every fcs file. Like for underway processing, class files include final class assingments for every particle, notes about which assignment is which,
         volume estimates for every particle and the bead values and filename used for volume estimates. 
    - A folder within outputs labelled "figs" with figures for gate assingments, grouped according to run type. 
    - SummaryTable.mat which translates Gated_Table into a table with a row for each Niskin, rather than each fsc file. 
        The fcs file used to report Syn, Euks, Prochlorococcus, and Bacteria are listed. 
    - EDI_table.mat, a version of the summary table with better column names, biovolume estimates etc and without prochlorococcus counts. This table will be remade in the next step. 


2. Once AllSteps function has been run, run "format_table_for_EDI_group"

Make sure that list of directories at the top includes all those cruises you want included in the merged table. 
Then run. 
This script will overwrite the EDI tables for each of the curises made in the previous step, so that they are up to date. 

This script generates: 
      the EDI_Table.mat files within each cruise output folder, 
        and a merged table for all the cruises saved as a mat file and a csv: 
      '\\sosiknas1\Lab_data\Attune\cruise_data\Attune_Discrete_Table.mat'
      '\\sosiknas1\Lab_data\Attune\EDI_data_packages\Attune_transect_FCMdiscrete\attune-transect-discrete-samples.csv'

3. quality control step: "compare_filenames_to_samplelog.m"