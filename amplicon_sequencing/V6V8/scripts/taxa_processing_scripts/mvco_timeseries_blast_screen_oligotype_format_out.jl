#Import .tsv taxbyseq file, create unique sequence identifier, create fasta file for blast
#run blast from shell
#import results, screen and parse, write out fasta file of Syn sequences for next steps

#Kristen R. Hunter-Cevera, Marine Biological Laboratory, Spring 2020

using DataFrames, CSV
using BioSequences
using Gadfly, Compose
##

# ind2 = findall(map(x -> occursin(r"Synecho",x), taxbyseq_df2[:,:Taxonomy])) #make a Syn only database
# taxbyseq_df2=taxbyseq_df2[ind2,:]

#import taxbyseq file:
mvco_df = DataFrame(CSV.File("/home/kristen/Documents/V6V8_analysis/data/vamps_downloads/taxbyseq_tables/mvco_timeseries_taxbyseq.tsv"; delim = '\t',header=2))
mvco_df.uniq_ID=collect(1:size(mvco_df,1)) #just use numbers for now?

#easier handling:
for n in names(mvco_df[!,Between(Symbol("MVCO_2010_2018_timeseries--MV254_run2"),Symbol("MVCO_2010_2018_timeseries--MV400Ster_run7"))])
     rename!(mvco_df, n => replace(string(n),"MVCO_2010_2018_timeseries--" => ""))
end

## setup for blast check:
fastafile="/home/kristen/Documents/V6V8_analysis/analysis_products/mvco_timeseries_uniq.fasta"

w = open(FASTA.Writer, fastafile)
for q=1:size(mvco_df,1)
    rec = FASTA.Record(string(mvco_df[q,:uniq_ID]),mvco_df[q,:Sequence])
    write(w, rec)
end
close(w)

## Run blast:

outfile="/home/kristen/Documents/V6V8_analysis/analysis_products/blast/mvco_timerseries_blast"
run(`bash /home/kristen/Documents/V6V8_analysis/scripts/taxa_processing_scripts/run_blast.sh $fastafile $outfile`)


## once have run blast, reimport results:

#convert output file from blastn search to DataFrame
blast_df = DataFrame(CSV.File("/home/kristen/Documents/V6V8_analysis/analysis_products/blast/mvco_timeseries_blast"; delim = '\t', header = false))

#put headers on output file
blast_col = ["uniq_ID", "sseqid" ,"pident" ,"length", "mismatch","gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"];
names!(blast_df, Symbol.(blast_col))

## merge dataframes:

test_df=join(mvco_df,blast_df,on = :uniq_ID, kind=:left)

test_df.query_length=abs.(test_df[!,:qstart] .- test_df[!,:qend])
test_df.subject_length=abs.(test_df[!,:sstart] .- test_df[!,:send])

## Sometimes blast will come back with two answers for a sequence match -
#not sure if that's a threading issue or something being mismanaged
#this is despite the the fact that only one answer was requested from the program!

#Upon merging, if the dataframe is presented with two options - I think it jsut duplicates the line to add both matches

t=convert(Array,test_df[!,:uniq_ID]) #merged dataframe is longer than mvco_df
g=setdiff(1:size(mvco_df,1),indexin(unique(t), t)) #find the indices of duplicate unqiue IDs
gg=unique(test_df[g,:uniq_ID])

global rec=zeros(size(gg,1),1)

for q=1:length(gg)

    ii=findall(test_df[!,:uniq_ID] .== gg[q]) #find matching indexes
    println(size(ii))
    k=argmin(test_df[ii,:bitscore])
    deleterows!(test_df,ii[k]) #remove the result that had the lower bit score, but keep record of lines removed!
    rec[q]=ii[k]

end

#size of rec should match size of gg!

## great, now check how the syn sequences look:

#find all the Syn based on GAST matches:
synind=findall(occursin.("Synecho",test_df[!,:Taxonomy]))

## prelim plot:

Gadfly.set_default_plot_size(20cm,20cm)
p1=plot(x=test_df[synind,:pident], Geom.histogram,Guide.title("Percent ID")) #most should be > 95
p2=plot(x=test_df[synind,:mismatch], Geom.histogram,Guide.title("Mismatch")) #most should be < 10
p3=plot(x=test_df[synind,:gapopen], Geom.histogram,Guide.title("Gaps Open")) #most should be < 3
p4=plot(x=test_df[synind,:bitscore], Geom.histogram,Guide.title("BitScore")) #should be above 800

gridstack([p1 p2; p3 p4])

##

synind=findall(occursin.(r"Synecho",test_df[!,:Taxonomy]))
noblast_hitsyn = findall(ismissing.(test_df[synind,:pident])) #any syn (classified by gast) that did not match to our blast database
#that should be empty or no hits - all gast syn basically had blast syn hits
ii=setdiff(1:length(synind),noblast_hitsyn)
synind=synind[ii]

#what do min values look like?
jj=findall(test_df[synind,:bitscore] .<= 700)
minimum(test_df[synind,:bitscore])

findall(test_df[synind,:pident] .<= 90)
minimum(test_df[synind,:pident])

findall(test_df[synind,:mismatch] .>= 10)
maximum(test_df[synind,:mismatch])

findall(test_df[synind,:gapopen] .>= 3)
maximum(test_df[synind,:gapopen])

# Some surprisingly low values, but I think for this, we will keep these sequences in and let oligotyping sort it out!

## 🌵 Now, do any sequences need to be added?
# I.e. are there Syn sequences not tagged by GAST but captured by blast?

noblast = findall(ismissing.(test_df[!,:pident])) #no blast results to syn
synblast = setdiff(1:size(test_df,1),noblast) #total blast results that had match to syn
nonsyn = setdiff(synblast,synind) #remove gast matches to syn, what came back as syn?

plot(x=test_df[nonsyn,:bitscore], Geom.histogram)

kk=findall(test_df[nonsyn,:bitscore] .>= 800)
jj=findall(test_df[nonsyn,:bitscore] .>= 700)
unique(test_df[nonsyn[kk],:Taxonomy])
#minimum(test_df[synind,:bitscore])

sum([sum(row) for row = eachrow(test_df[nonsyn[jj],Between(:MV254_run2,:MV400Ster_run7)])])
#So, only a few hundred sequences come back as not-Syn by GAST, but a very good Syn Blast match to our database

#check in general these categories:
ii=findall(test_df[!,:Taxonomy] .== "Bacteria;Cyanobacteria;empty_class;Chroococcales;empty_family")
ii=findall(test_df[!,:Taxonomy] .== "Bacteria;Cyanobacteria;empty_class")
ii=findall(occursin.("Prochloro",test_df[!,:Taxonomy]))
sum([sum(row) for row = eachrow(test_df[ii,Between(:MV254_run2,:MV400Ster_run7)])]) #low pro count! good!

## with the exception of matches to Pro, we're gonna inlcude these extra blast hits :D
#adjusted synind:

ii=findall(occursin.("Prochloro",test_df[!,:Taxonomy]))
ind=setdiff(nonsyn[jj],ii) #jj here refers to bitscore threshold matches
unique(test_df[ind,:Taxonomy])

newind=union(synind,ind)


## Write sequences for Syn to fasta file, in oligotype format

df=test_df[newind,:] #not a deep copy

w = open(FASTA.Writer, "/home/kristen/Documents/V6V8_analysis/analysis_products/syn_sequences/mvco_syn_timeseries_for_oligotyping.fasta")

##
global syncount=DataFrame(event=Union{String,Missing}[],count=Union{Float64,Missing}[])
for n in names(df[!,Between(:MV254_run2,:MV400Ster_run7)])  #names(df[!,Between(:SynA,:SynD)]) # #names(df[!,Between(:MV263_rd2,:MV375)])

    # ## We need to split and clean up dataframe a little bit as there are repeat samples, OOI cruise samples and mock samples...
    # if occursin(r"--\d{3}",string(n)) #mvco sample

        ind=findall(map(x -> x != 0, df[:,n])) #returns the indexes that are not zero
        #df[ind,[n, 68]]
        @show n
        global test=sum(df[ind,n])
        #m=match(r"(?<sample>(M[1,2]{1}|Syn[A,B,C,D]{1}))",string(n))

        m=match(r"(?<sample>MV\d{3}[A,B]{0,1}\w*)",string(n)) #for MVCO time series
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

        push!(syncount,[m[:sample],count])
        @show test
        @show count
        if count != test
            print("Nooooo - sum does not match sequences printed?")
        end
    #end
end

##
close(w)

#and can check sums here - compare with grep commands on just made fasta file:
sum([sum(row) for row = eachrow(test_df[newind,Between(:MV254_run2,:MV400Ster_run7)])]) #Excellent, numbers should match!

#and if need total sequence counts per sample:
syncount.totals=[sum(cols) for cols = eachcol(test_df[!,Between(:MV254_run2,:MV400Ster_run7)])]
syncount |> CSV.write("/home/kristen/Documents/V6V8_analysis/analysis_products/syn_sequences/mvco_total_and_syn_counts.csv",delim=",")


## OLDER, PREVIOUS CODE:

# m=match(r"(?<event>(\d{3}[A,B]{0,1}_tr|\d{3}[A,B]{0,1}))","MV263A_tr")

#remodel this database to screen out non-mvco timeseries samples: NO LONGER NEEDED SINCE REUPLOAD TO VAMPS!
# redo5=("263","325","336","337","339","357","358")
# redo7=("308","309","310","338")
#
# for name in names(taxbyseq_df)
#         if occursin(r"(--C\d{1,2}N\d{1,2})|(MVBUCKET)|(MVM\d{1})|(Dock)|(SynCommunity)",string(name)) #for now, don't include these samples in analysis
#                 println(name)
#                 delete!(taxbyseq_df,name)
#         elseif occursin(r"MVCO_V6V8_Feb2016--",string(name)) #remove datasetname
#                 rename!(taxbyseq_df, name => replace(string(name),"MVCO_V6V8_Feb2016--" => ""))
#         elseif occursin(r"MVCO_ALLREADS_20190108--",string(name)) #remove dataset name and add MV
#                 m=match(r"MVCO_ALLREADS_20190108--(?<event>\d{3})",string(name))
#                 if m[:event] in redo5
#                         rename!(taxbyseq_df, name => string("MV",m[:event],"_rd2"))
#                 else
#                         rename!(taxbyseq_df, name => string("MV",m[:event]))
#                 end
#         elseif occursin(r"MVCO_TestRun2_August2018--",string(name)) #remove dataset name and add MV
#                 m=match(r"MVCO_TestRun2_August2018--(?<event>MV363A|MV\d{3})",string(name))
#                 rename!(taxbyseq_df, name => string(m[:event],"_tr2"))
#         elseif occursin(r"MVCO_timeseries_run7--",string(name)) #remove dataset name and add MV
#                 m=match(r"MVCO_timeseries_run7--(?<event>\d{3}(Milli|Ster){0,1})",string(name))
#                 if m[:event] in redo7
#                     rename!(taxbyseq_df, name => string("MV",m[:event],"_rd2"))
#                 else
#                     rename!(taxbyseq_df, name => string("MV",m[:event]))
#                 end
#         end
# end

## alter database to only keep mock community samples:
# taxbyseq_df=taxbyseq_df2
# for name in names(taxbyseq_df)
#         if occursin(r"(SynCommunity|MVM\d{1})|Sequence",string(name)) #for now, don't include these samples in analysis
#             println(name)
#             if occursin(r"MVCO_V6V8_Feb2016--",string(name)) #remove datasetname
#                 rename!(taxbyseq_df, name => replace(string(name),"MVCO_V6V8_Feb2016--" => ""))
#             elseif occursin(r"MVCO_timeseries_run7--",string(name)) #remove dataset name
#                 rename!(taxbyseq_df, name => replace(string(name),"MVCO_timeseries_run7--" => ""))
#             end
#         else
#             delete!(taxbyseq_df,name)
#         end
# end
