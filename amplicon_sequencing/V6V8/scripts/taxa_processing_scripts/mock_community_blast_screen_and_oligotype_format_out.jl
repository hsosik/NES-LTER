#For each dataset (mock communities and mvco timeseries):
#import tsv taxbyseq file, create unique sequence identifier, create fasta file for blast
#run blast from shell
#import results, screen and parse, write out fasta file of Syn sequences for next steps

#Kristen R. Hunter-Cevera, Marine Biological Laboratory, Spring 2020

using DataFrames, CSV
using BioSequences
using Gadfly
using JLD2


##
mock_df = DataFrame(CSV.File("/home/kristen/Documents/V6V8_analysis/data/vamps_downloads/taxbyseq_tables/mock_community_taxbyseq.tsv"; delim = '\t',header = 2))

#for now, remove dock samples:
delete!(mock_df,Symbol("syn_mock_communities--DockA_run7"))
select!(mock_df, Not( Symbol("syn_mock_communities--DockB_run7")))

rename!(mock_df,Symbol("syn_mock_communities--SynCommunityA_run7") => :SynA)
rename!(mock_df,Symbol("syn_mock_communities--SynCommunityB_run7") => :SynB)
rename!(mock_df,Symbol("syn_mock_communities--SynCommunityC_run7") => :SynC)
rename!(mock_df,Symbol("syn_mock_communities--SynCommunityD_run7") => :SynD)
rename!(mock_df,Symbol("syn_mock_communities--M1_run3")  => :M1)
rename!(mock_df,Symbol("syn_mock_communities--MVM2_run2") => :M2)

#remove rows that have totals of zeros after dock samples removed:
[sum(col) for col = eachcol(mock_df[!,[:SynA, :SynB]])]
zeroind=findall([sum(row) for row = eachrow(mock_df[!,Between(:M1,:SynD)])] .== 0)
deleterows!(mock_df, zeroind)
synind = findall(map(x -> occursin(r"Synecho",x), mock_df[:,:Taxonomy]))
#colwise(sum,mock_df[[ Symbol("syn_mock_communities--SynCommunityA_run7"),  Symbol("syn_mock_communities--SynCommunityB_run7")]])
mock_df.uniq_ID=collect(1:size(mock_df,1))

## write FASTA file of unqiue sequences with unique identifiers:
fastafile="/home/kristen/Documents/V6V8_analysis/analysis_products/blast/mock_communities.fasta"
w = open(FASTA.Writer, fastafile)
for q=1:size(mock_df,1)
    rec = FASTA.Record(string(mock_df[q,:uniq_ID]),mock_df[q,:Sequence])
    write(w, rec)
end
close(w)

## run blast:

outfile="/home/kristen/Documents/V6V8_analysis/analysis_products/blast/mock_community_blast_results"
run(`bash /home/kristen/Documents/V6V8_analysis/scripts/taxa_processing_scripts/run_blast.sh $fastafile $outfile`)


## Import the blast output file and merge with database:

#convert output file from blastn search to DataFrame
blast_df = DataFrame(CSV.File("/home/kristen/Documents/V6V8_analysis/analysis_products/blast/mock_community_blast_results"; delim = '\t', header = false))

#put headers on output file
blast_col = ["uniq_ID", "sseqid" ,"pident" ,"length", "mismatch","gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"];
names!(blast_df, Symbol.(blast_col))

#merge dataframes:
test_df=join(mock_df,blast_df,on = :uniq_ID, kind=:left)
test_df.query_length=abs.(test_df[!,:qstart] .- test_df[!,:qend])
test_df.subject_length=abs.(test_df[!,:sstart] .- test_df[!,:send])

synind=findall(occursin.(r"Synecho",test_df[!,:Taxonomy]))
noblast_hitsyn = findall(ismissing.(test_df[synind,:pident]))
#good that this is empty - all gast syn had blast syn hits
ii=setdiff(1:length(synind),noblast_hitsyn)
synind=synind[ii]


## some logic tests to see how pident, bitscore, etc relate....plus some plots

# 2 QUESTIONS:
#Of the sequences classified as Syn by GAST - how many are not?
#Of the sequences not classified as Syn by GAST - how many might be?


#🌵 Question 1:

#histograms to show distribution of blast results:
Gadfly.set_default_plot_size(20cm,20cm)
p1=plot(x=test_df[synind,:pident], Geom.histogram,Guide.title("Percent ID")) #most should be > 95
p2=plot(x=test_df[synind,:mismatch], Geom.histogram,Guide.title("Mismatch")) #most should be < 10
p3=plot(x=test_df[synind,:gapopen], Geom.histogram,Guide.title("Gaps Open")) #most should be < 3
p4=plot(x=test_df[synind,:bitscore], Geom.histogram,Guide.title("BitScore")) #should be above 800

gridstack([p1 p2; p3 p4])

##
tt=findall(test_df[synind,:pident] .<= 90) #so not too many...and usually low abundance
[sum(col) for col = eachcol(mock_df[synind[tt],Between(:M1,:SynD)])]

## if needed, write these out to blast against full nucleotide database:
w = open(FASTA.Writer, "/home/kristen/Documents/V6V8_analysis/mock_syngast_lowPID.fasta")
for q=1:size(tt,1)
    rec = FASTA.Record(string(mock_df[synind[tt[q]],:uniq_ID]),mock_df[synind[tt[q]],:Sequence])
    write(w, rec)
end
close(w)

# so, all in all for the mock syn - no worries for question 1! ✅



#🌵 Question 2:
noblast = findall(ismissing.(test_df[!,:pident])) #no blast results to syn
synblast = setdiff(1:size(test_df,1),noblast) #total blast results that had match to syn
nonsyn = setdiff(synblast,synind) #remove gast matches to syn, what came back as syn?

plot(x=test_df[nonsyn,:pident], Geom.histogram) #most should be > 95
#most of these sequences are below 90...which is good, but what about the ones that arent't?

#look at the sequences that really show a good chance of being a coastal syn:
#tt=findall((test_df[nonsyn,:pident] .>= 90) .& (test_df[nonsyn,:query_length] .>= 440) .& (test_df[nonsyn,:subject_length] .>= 440)
#    .& (test_df[nonsyn,:gapopen] .<= 3) .& (test_df[nonsyn,:mismatch] .<= 10) .& (test_df[nonsyn,:bitscore] .>= 800)) #only consider full length queries

#-or-:

#because these are mock communities - what are the other cyanos?
ss=findall(map(x -> occursin(r"Cyano",x),test_df[nonsyn,:Taxonomy]))
[sum(col) for col = eachcol(test_df[nonsyn[ss],Between(:M1,:SynD)])] #still very few

tt=findall(test_df[nonsyn,:bitscore] .>= 700) #only consider full length queries
plot(x=test_df[nonsyn[ss],:mismatch], Geom.histogram)

taxa_to_check=unique(test_df[nonsyn[tt],:Taxonomy])

[sum(col) for col = eachcol(test_df[nonsyn[tt],Between(:M1,:SynD)])] #still very few

rowsum=[sum(row) for row = eachrow(test_df[nonsyn[tt],Between(:M1,:SynD)])] #still very few
sum(rowsum) #this is also a low amount...

# ii=findall(map(x -> occursin(r"Prochlorococcus",x),test_df[nonsyn[tt],:Taxonomy]))

## hmmm...so some cyanos that are matching to Pro, Cyanobium, etc. Parse these for blast nucleotide check!
wcyano = open(FASTA.Writer, "/home/kristen/Documents/V6V8_analysis/mock_nonsyngast_highPID_cyano.fasta")
wother = open(FASTA.Writer, "/home/kristen/Documents/V6V8_analysis/mock_nonsyngast_highPID_other.fasta")

for q=1:size(tt,1)
    if occursin(r"Cyano",mock_df[nonsyn[tt[q]],:Taxonomy])
        rec = FASTA.Record(string(mock_df[nonsyn[tt[q]],:uniq_ID]),mock_df[nonsyn[tt[q]],:Sequence])
        write(wcyano, rec)
    else
        rec = FASTA.Record(string(mock_df[nonsyn[tt[q]],:uniq_ID]),mock_df[nonsyn[tt[q]],:Sequence])
        write(wother, rec)
    end
end
close(wcyano)
close(wother)
## hmmmm...so it does seem that these sequences could be considered Syn..despite being labeled as Pro or a Microcytis...probably should include them!
# Write out sequences in oliogytping format!

#index to use:
sind=union(synind,nonsyn[tt])
df=test_df[sind,:]
# Compile for oligotyping! woot-woot!

#wA = open(FASTA.Writer, "/home/kristen/Documents/V6V8_analysis/analysis_products/syn_sequences/mock_synseqs_2016exp.fasta")
w = open(FASTA.Writer, "/home/kristen/Documents/V6V8_analysis/analysis_products/syn_sequences/mock_synseqs_2018exp.fasta")

##

for n in names(df[!,Between(:SynA,:SynD)]) #names(df[!,Between(:M1,:M2)])#names(df[!,Between(:SynA,:SynD)]) #names(df[!,Between(:MV308_rd2,:MV400Ster)])  #names(df[!,Between(:MV263_rd2,:MV375)])

        ind=findall(map(x -> x != 0, df[:,n])) #returns the indexes that are not zero
        #df[ind,[n, 68]]
        @show n
        global test=sum(df[ind,n])
        m=match(r"(?<sample>(M[1,2]{1}|Syn[A,B,C,D]{1}))",string(n))

        #m=match(r"(?<sample>MV\d{3}[A,B]{0,1}\w*)",string(n)) #for MVCO time series
        #@show m[:smaple]

        global count=0 #important for oligotyping file - adds number to each sequence
        for j=1:length(ind) #different Syn seqs in this sample
            for i=1:df[ind[j],n] #for each syn seq, duplicate to create fasta file for oligotyping
                #@show
                global count += 1
                rec = FASTA.Record(string("Sample-",m[:sample],"_Read",count),df[ind[j],:Sequence])
                write(w, rec)
            end
        end
        @show test
        @show count
        if count != test
            print("Nooooo - sum does not match sequences printed?")
        end
    #end
end

close(w)
