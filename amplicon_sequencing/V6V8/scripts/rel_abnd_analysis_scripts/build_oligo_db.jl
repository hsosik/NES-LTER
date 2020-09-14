#Import oligotype counts and MVCO data - combine to form one DataFrame
#Kristen Hunter-Cevera, Marine Biological Laboratory, Spring 2020

using DataFrames, CSV
using Statistics
using JLD2
using Dates

## IMPORT OLIGO AND SYN COUNTS:

filename="/home/kristen/Documents/V6V8_analysis/analysis_products/syn_oligotyping/mvco_syn/mvco_syn_timeseries_for_oligotyping_pynast_aligned-c19-s1-a0.0-A0-M200/MATRIX-COUNT.txt"  #taxbyseq_nov2018_only.tsv"
#filename="/home/kristen/Documents/V6V8_analysis/data/synseqs/test-c25-s1-a0.0-A100-M40/MATRIX-COUNT.txt"  #taxbyseq_nov2018_only.tsv"
oligo_counts_df = CSV.File(filename, header=1) |> DataFrame!

for j=2:size(oligo_counts_df,2) #rename columns
    rename!(oligo_counts_df,:($j) => "Oligo$(j-1)")
end

#detect and add event number as separate column:
nums=map(x-> match(r"(363[A,B])|(\d{3})",x),oligo_counts_df[:,1]) #special handling for 363A and 363B - event numbers
nums2=string.(collect(m.match for m in nums))
oligo_counts_df[:event]=nums2 #make it into a column!

oligo_counts_df.samples=map(x -> replace(x,"Sample-" => ""),oligo_counts_df[!,:samples])

# ugh=string.(oligo_counts_df[!,:samples])
# [replace(row,"Sample-" => "") for row in ugh]

oligo_counts_df=oligo_counts_df[!,[1;size(oligo_counts_df,2);2:(size(oligo_counts_df,2)-1)]] #reshuffle column names

# Import total syn reads and total seq reads....syn counts for the moment comes from o-get-sample...but should have this as out when prints oligo seqs!
filename="/home/kristen/Documents/V6V8_analysis/analysis_products/syn_sequences/mvco_total_and_syn_counts.csv"  #taxbyseq_nov2018_only.tsv"
total_syncounts_df = CSV.File(filename, header=1) |> DataFrame!
rename!(total_syncounts_df,:event => :samples)
rename!(total_syncounts_df,:count => :syncounts)
rename!(total_syncounts_df,:totals => :totalcount)
oligo_counts_df=join(oligo_counts_df,total_syncounts_df, on = :samples)

# Finished importing oligo and syn counts!

## Combine, tr (trial run) and rd (redo) samples:

eventlist=unique(oligo_counts_df[:,:event])
ind_rec=[]

for q=1:size(eventlist,1)
    ind=findall(map(x->occursin(Regex(eventlist[q]),x),oligo_counts_df[:,:event]))

    if size(ind,1) > 1
        @show q
        append!(ind_rec,ind)
        #@show B[q,:event]
        combo=sum.(eachcol(oligo_counts_df[ind,3:size(oligo_counts_df,2)]))
        A=[string("Sample-MV", oligo_counts_df[ind[1],:event],"_combined"),oligo_counts_df[ind[1],:event],combo...]
        push!(oligo_counts_df,A)
        deleterows!(oligo_counts_df,ind)
        #can push rows as NamedTuples
    end
end

##FINALLY..................................................................






########## MATCH TO MVCO DATES, TEMPERATURE, LIGHT, ETC #############################
#from MVCO_nut_reps, match event number to date, from date match temp, light, nutrients

#load in MVCO environmental and syn data:
@load "/home/kristen/Documents/V6V8_analysis/analysis_products/mvco_nut_data.jld2"


#There is a real different between NaN and missing; NaN is for calc like Num/0 or something
#Missing truly refers to a missing value....
## replace NaN's with missing values:
for q=7:19
    @. MVCO_nut_reps[isnan(MVCO_nut_reps[:,q]),q] = missing
end
#database entry mistake?
ii=findall(map(x -> occursin(r"273",x),MVCO_nut_reps[:,1]))
MVCO_nut_reps[ii[1],7]=0.0

# a[isnan.(a)] .= 0
# or you can use the @. macro to dot everything for you:
# @. a[isnan(a)] = 0


nutrient_indexes=collect(findall(map(x->occursin(Regex(n),x),MVCO_nut_reps[:,1])) for n in oligo_counts_df[:,:event])
missing_a_date=findall(map(x->size(x,1) == 0,nutrient_indexes))
#seems 272 does not have an entry in MVCO_nut_reps!

## DATES AND NUTRIENTS:
 #event_dates=Dict() #as a dictionary, yes...but want in dataframe:
datevec=Array{DateTime}(undef,size(nutrient_indexes,1))
NO3vec=Array{Union{Float64,Missing}}(undef,size(nutrient_indexes,1))
NH4vec=Array{Union{Float64,Missing}}(undef,size(nutrient_indexes,1))
PO4vec=Array{Union{Float64,Missing}}(undef,size(nutrient_indexes,1))
SiOHvec=Array{Union{Float64,Missing}}(undef,size(nutrient_indexes,1))
NO3vecsurf=Array{Union{Float64,Missing}}(undef,size(nutrient_indexes,1))
NH4vecsurf=Array{Union{Float64,Missing}}(undef,size(nutrient_indexes,1))
PO4vecsurf=Array{Union{Float64,Missing}}(undef,size(nutrient_indexes,1))
SiOHvecsurf=Array{Union{Float64,Missing}}(undef,size(nutrient_indexes,1))

for j =1:size(nutrient_indexes,1)
    if !isempty(nutrient_indexes[j])
        #@show j
        datevec[j]=Date(split(MVCO_nut_reps[nutrient_indexes[j][1],3]," ")[1],"yyyy-mm-dd")
        #event_dates[nums2[j]]=Date(split(MVCO_nut_reps[t[j][1],3]," ")[1],"yyyy-mm-dd") #for dictionary

        #average of all depths:
        NO3vec[j]=mean(skipmissing(MVCO_nut_reps[nutrient_indexes[j],8:10]))
        NH4vec[j]=mean(skipmissing(MVCO_nut_reps[nutrient_indexes[j],11:13]))
        PO4vec[j]=mean(skipmissing(MVCO_nut_reps[nutrient_indexes[j],17:19]))
        SiOHvec[j]=mean(skipmissing(MVCO_nut_reps[nutrient_indexes[j],14:16]))

        #average of 'surface' measurements < 6m deep
        temp_ind=nutrient_indexes[j][findall(map(x -> x<=6, MVCO_nut_reps[nutrient_indexes[j],7]))]
        @show j
        @show temp_ind
        NO3vecsurf[j]=mean(skipmissing(MVCO_nut_reps[temp_ind,8:10]))
        NH4vecsurf[j]=mean(skipmissing(MVCO_nut_reps[temp_ind,11:13]))
        PO4vecsurf[j]=mean(skipmissing(MVCO_nut_reps[temp_ind,17:19]))
        SiOHvecsurf[j]=mean(skipmissing(MVCO_nut_reps[temp_ind,14:16]))

    else
         @show j
         println("Special case for event 272")
         datevec[j]=Date("2011-08-06","yyyy-mm-dd")
         NO3vec[j]=missing
         NH4vec[j]=missing
         PO4vec[j]=missing
         SiOHvec[j]=missing
         NO3vecsurf[j]=missing
         NH4vecsurf[j]=missing
         PO4vecsurf[j]=missing
         SiOHvecsurf[j]=missing
         #event_dates[nums2[j]]=Date("2011-08-06","yyyy-mm-dd") #for dictionary
    end
end

# collect(d[1] for d in t) #would have worked if not for empty vector....
oligo_counts_df[:date]=datevec
oligo_counts_df[:NO3]=NO3vec
oligo_counts_df[:NH4]=NH4vec
oligo_counts_df[:PO4]=PO4vec
oligo_counts_df[:SiOH]=SiOHvec
oligo_counts_df[:NO3surf]=NO3vecsurf
oligo_counts_df[:NH4surf]=NH4vecsurf
oligo_counts_df[:PO4surf]=PO4vecsurf
oligo_counts_df[:SiOHsurf]=SiOHvecsurf
## woo-hoo!


## MATCH TO SYN ABND, LIGHT AND TEMPERATURE!
@load "/home/kristen/Documents/V6V8_analysis/analysis_products/mvco_envdata.jld2" Tbeam_corr light_int
@load "/home/kristen/Documents/V6V8_analysis/analysis_products/mvco_syndata.jld2" time_syn daily_syn
# ind=findall(map(x -> x == datevec[1],time_syn))
# daily_syn[ind] #yes! this works!
##
ind=collect(findall(map(x -> x == dd,time_syn)) for dd in datevec) #use date as the index....

synvec=Array{Float64}(undef,size(nutrient_indexes,1))
tempvec=Array{Float64}(undef,size(nutrient_indexes,1))
lightvec=Array{Float64}(undef,size(nutrient_indexes,1))
lightvec_wkavg=Array{Float64}(undef,size(nutrient_indexes,1))

Tbeam_corr=recode(Tbeam_corr, NaN => missing)
light_int=recode(light_int, NaN => missing)

light_avg=vcat((mean(skipmissing(x)) for x in eachrow(light_int))...) #thanks stackoverflow!
temp_avg=vcat((mean(skipmissing(x)) for x in eachrow(Tbeam_corr))...) #thanks stackoverflow!



for j=1:size(ind,1)

    i0=LinearIndices(light_int)[ind[j]]
    i0=i0[1]
    i1=ind[j][1][1] #row index
    i2=ind[j][1][2] #col index

    synvec[j]=daily_syn[i0]

    if ismissing(Tbeam_corr[i0])
        @show j
        tempvec[j]=mean(skipmissing(Tbeam_corr[i0-1:i0+1]))
        #println("missing temp")
        #println(tempvec[j])

        if ismissing(tempvec[j]) | isnan(tempvec[j]) #and in case do not have even surrounding data:
            @show j
            tempvec[j]=temp_avg[i1]
            println("too much missing for temp!")
            println(tempvec[j])
        end

    else
        tempvec[j]=Tbeam_corr[i0]
    end

    if ismissing(light_int[i0])

        lightvec[j]=light_avg[i1]
        #println("missing light")
        #println(lightvec[j])
    else
        lightvec[j]=light_int[i0]
    end

    lightvec_wkavg[j] = mean(skipmissing(light_int[i0-7:i0]))  #averages the week before the sample time point

end

##
oligo_counts_df[:synconc]=synvec
oligo_counts_df[:temp]=tempvec
oligo_counts_df[:light]=lightvec
oligo_counts_df[:lightweek]=lightvec_wkavg
#Got it! Now syn oligos with matching env data is complete!


sort!(oligo_counts_df,:date)
oligo_counts_df.month=@. Dates.month(oligo_counts_df[:,:date])


@save "/home/kristen/Documents/V6V8_analysis/analysis_products/mvco_syn_oligo_counts_plus_env_data.jld2" oligo_counts_df









#plotting..............................
# gr()
# scatter([1;2],[3;4])
# ugh=plot(B[:,:date], B[:,:Oligo1])
# savefig(ugh,"test.png")
# scatter!(B[:date], B[:synconc],yscale = :log10)
# xaxis!("Year day",(0, 366),0:100:360,font(12, "Times"))
