#Script that matches oligotype sequences to database of full length Synechococcus sequences

using BioSequences
using BioAlignments
using JLD2
using DataFrames, CSV

using Gadfly, Colors, Compose, ColorSchemes, Cairo, Fontconfig
include("alignment_functions.jl")
## load in that reference database!

#unclear exactly how this database was made in the first place....oops!
@load "/home/kristen/Documents/V6V8_analysis/analysis_products/marine_syn_updated_db.jld2" synref

#must change string sequences back to BioSequence types!
synref.sequence = DNASequence.(synref[!,:sequence])
synref.Length = length.(synref[!,:sequence])


## we need to do some clean up before can use it!

#exclude sequences that are less than full length: goes from 232 seqs to 229
ii=findall(synref[!,:Length] .> 700)
synref=synref[ii,:]

#to find out how many have unambiguous clades;
#findall((synref[!,:clade] .== "?") .| (synref[!,:clade] .== "-"))

## first shorten to V6-V8 region!
synref.V6V8=Array{Any,1}(missing,size(synref,1)) #keep the alignments for future reference!
for q=1:size(synref,1)
    synref[q,:V6V8]=find_V6V8_region(synref[q,:sequence])
end

synref.V6V8_length=map(x -> length(x),synref[!,:V6V8])

ii=findall(synref[!,:V6V8_length] .>= 300) #exclude shorter sequences - goes down to 225 sequences
synref=synref[ii,:]
#the double check:
#kk=findall(isempty.(synref[!,:V6V8]))



## alright, soooo, we need to show two types of aligment relationships:

# how doe strain reps in the database differ from each other?
# how are these connected to the oligotype sequences?

#matrix appears to be the easiest way to do this!

#we could just add the oligotype sequences in and designate as new 'clade'....

# ii=findall(synrefB[!,:clade] .== clade)
# x=synrefB[ii,:V6V8]
# idx = unique(z -> x[z], 1:length(x))


## per clade, find the unique sequences, including short ones:
clade_list=unique(synref[!,:clade]) #find unqiue sequences and reorganize!
# clade_list=sort(setdiff(clade_list,["-";"?"]))

#indexin(["5.3"],clade_list0)
## for this, it's just easier to do by hand!
clade_list=["5.3";
"5.2MV1";
 "5.2MV2";
 "5.2MV3";
 "CB4"   ;
 "CB5"   ;
 "I"    ;
 "II"    ;
 "III"   ;
 "IV"    ;
 "V"     ;
 "VI"    ;
 "VII"   ;
 "VIII";
  "IX" ;
 "CRD1" ;
 "WPC1"  ;
 "WPC2"  ;
 "XV"    ;
 "XVI"]
##
unqrefs=DataFrame(ID=String[],strain_reps=Set[],clade=String[],V6V8=BioSequence[])

for clade in clade_list

    println(string("analyzing sequences within clade: ",clade))
    jj=findall((synref[!,:clade] .== clade) .& (synref[!,:V6V8_length] .>300))
    withinclade_array=deepcopy(synref[jj,[:strain]])
    scoremodel = AffineGapScoreModel(match=5, mismatch=-4, gap_open=-10,gap_extend=-10)

    for q=1:size(jj,1)
        #println(q)
        temp_df=align_to_database2!(synref[jj[q],:V6V8],synref[jj,:],scoremodel,:V6V8,:strain) #align to unique sequences only here
        #synrefD[j,:alignments]=temp_df #record alignments
        tempID=synref[jj[q],:strain] #column name
        withinclade_array[Symbol(tempID)] = temp_df[!,:mismatch]
    end

    # now, sort through this identity matrix!
    global counter=1
    for q=1:size(jj,1)

        #first, check if this strain has already been accounted for:
        if !any(map(x -> withinclade_array[q,:strain] in x,unqrefs[!,:strain_reps])) #if not, then proceed!

            qq=findall(convert(Array{Any},withinclade_array[q,:]) .== 0)
            qq=qq .- 1 #to match columns

            if length(qq) == 1 #then just this strain has this sequence and it should be unique!
                push!(unqrefs,[string(clade,"-",counter),Set([synref[jj[q],:strain]]),clade,synref[jj[q],:V6V8]])
                global counter += 1
            else
                strains=Set(withinclade_array[qq,:strain])

                #indexes back into synref!
                qqind=indexin(withinclade_array[qq,:strain],synref[!,:strain])
                y=argmax(synref[qqind,:V6V8_length]) #index of longest v6v8 read in qqind

                push!(unqrefs,[string(clade,"-",counter),strains,clade,synref[qqind[y],:V6V8]])
                global counter += 1
            end
        end
    end
end

## sweet! and we can double check this:

findall(map(x -> "HB1133" in x,unqrefs[!,:strain_reps])) #this should just be one!
findall(map(x -> "RCC1026" in x,unqrefs[!,:strain_reps])) #this should just be one!

## add the lengths!
unqrefs.V6V8_length  = length.(unqrefs[!,:V6V8])


## and now, add in the oligotypes!

oligo_file = open(FASTA.Reader, "/home/kristen/Documents/V6V8_analysis/analysis_products/syn_oligotyping/mvco_syn/mvco_syn_timeseries_for_oligotyping_pynast_aligned-c19-s1-a0.0-A0-M200/oligo_repseqs.fasta") #test-c25-s1-a0.0-A100-M40/oligo_representative_no_gaps.fasta") #open the fasta file

for record in oligo_file
    println(string(FASTA.identifier(record)))
    push!(unqrefs, [string(FASTA.identifier(record)), Set([]),"",FASTA.sequence(record),length(FASTA.sequence(record))])
end
close(oligo_file)

## cool, now we can do alignments for each sequence to see how they match up!

clade_array=deepcopy(unqrefs[!,[:clade,:ID]])

scoremodel = AffineGapScoreModel(match=5, mismatch=-4, gap_open=-10,gap_extend=-10)
#scoremodel = AffineGapScoreModel(match=5, mismatch=-4, gap_open=-4,gap_extend=-1)

for q=1:size(unqrefs,1)

        println(q)
        temp_df=align_to_database2!(unqrefs[q,:V6V8],unqrefs,scoremodel,:V6V8,:ID) #align to unique sequences only here
        #synrefD[j,:alignments]=temp_df #record alignments
        tempID=unqrefs[q,:ID]
        clade_array[Symbol(tempID)] = temp_df[!,:mismatch]

end

## a double check on the oligotypes:

k=findall(clade_array[!,:ID] .== "O3")
kk=findall(convert(Array{Any},clade_array[k,:]) .== 0)
qq= [ i[2] for i in kk ]
qq=qq .- 2
println(clade_array[qq,:ID])
## sweet! and the plotting!


spy(convert(Array,clade_array[:,3:end]))

xlabel=unqrefs[!,:ID]

xticks = collect(1:size(unqrefs,1))
labels2 = Dict(zip(xticks, xlabel))

spy(convert(Array,clade_array[:,3:end]),
    Scale.x_discrete(labels = x -> labels2[x]),
    Scale.y_discrete(labels = x -> labels2[x]))

## save these matches

@save "/home/kristen/Documents/V6V8_analysis/analysis_products/oligo_clade_match_dataframes.jld2" unqrefs clade_array

## and export for fig in matlab, if needed:

CSV.write("/home/kristen/Documents/V6V8_analysis/scripts/rel_abnd_analysis_scripts/matlab_scripts/V6V8_matrix_match.csv", clade_array, delim = ',')


# notes; ah-ya! O3 has additional matches to M11.1 and M11.2, but clade designation is ambiguous from full length 16S
#O14 also matches M12.1, which is clade II from 16S...

## SI table to show sequences used:

synref_for_table=deepcopy(synref[!,[:clade,:strain,:accession,:Source,:V6V8_length]])

## some housekeeping:

#exclude the strains that don't have clade designation:
ii=findall((synref[!,:clade] .== "?") .| (synref[!,:clade] .== "-"))
jj=setdiff(collect(1:size(synref,1)),ii)
synref_for_table=synref_for_table[jj,:]

#replace!(synref_for_table[!,:Source],"Hunter-Cevera et al. 2016" => "Hunter-Cevera 2014")


## now match each strain to a unqiue V6V8 sequence yo!

v6v8_match=Array{Union{Missing, String}}(missing, size(synref_for_table,1))
for q=1:size(synref_for_table,1)
    println(q)
    qq=findall(map(x -> synref_for_table[q,:strain] in x,unqrefs[!,:strain_reps]))
    v6v8_match[q] = unqrefs[qq,:ID][1]
end

synref_for_table.V6V8_unq_match = v6v8_match

##

sort!(synref_for_table,:V6V8_unq_match)
CSV.write("/home/kristen/Documents/V6V8_analysis/paper_V6V8/PNAS_style/syn_ref_table.csv", synref_for_table, delim = ",")


## If need to write out these sequences in a fasta format:

synref_for_table_seqs=deepcopy(synref[!,[:clade,:strain,:accession,:Source,:sequence]])

#exclude the strains that don't have clade designation:
ii=findall((synref[!,:clade] .== "?") .| (synref[!,:clade] .== "-"))
jj=setdiff(collect(1:size(synref,1)),ii)
synref_for_table_seqs=synref_for_table_seqs[jj,:]

fastafile="/home/kristen/Documents/V6V8_analysis/analysis_products/marine_syn_db.fasta"

w = open(FASTA.Writer, fastafile)
for q=1:size(synref_for_table_seqs,1)
    rec = FASTA.Record(string(synref_for_table_seqs[q,:strain]),string(synref_for_table_seqs[q,:accession]),synref_for_table_seqs[q,:sequence])
    write(w, rec)
end
close(w)


## or if want to put this table in latex:
t=Array{Union{Missing, String}}(missing, size(synref_for_table,1))
t .= "&"

n=Array{Union{Missing, String}}(missing, size(synref_for_table,1))
n .= "\\"

insert!(synref_for_table, 2, t, :L1)
insert!(synref_for_table, 4, t, :L2)
insert!(synref_for_table, 6, t, :L3)
insert!(synref_for_table, 8, t, :L4)
insert!(synref_for_table, 10, t, :L5)
insert!(synref_for_table, 12, n, :E)

##

CSV.write("/home/kristen/Documents/V6V8_analysis/scripts/rel_abnd_analysis_scripts/matlab_scripts/syn_ref_table.csv", synref_for_table, delim = "   ")



## if need to double check alignments directly:

scoremodel = AffineGapScoreModel(match=5, mismatch=-4, gap_open=-10,gap_extend=-10)

# s1 = synrefB[ii[idx[1]],:V6V8];
# s2 = synrefB[ii[idx[2]],:V6V8];
s1 = unqrefs[12,:V6V8]#unqseq[6,:V6V8];
s2 = unqrefs[18,:V6V8];
res = pairalign(LocalAlignment(), s1, s2, scoremodel)
aln=alignment(res)
