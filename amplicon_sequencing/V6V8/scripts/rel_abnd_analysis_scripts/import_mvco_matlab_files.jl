#Import matlab files and saving them as .JLD2 files!

#Kristen R. Hunter-Cevera, Marine Biological Laboratory, Spring 2020

using MAT
using Dates
using JLD2
include("/home/kristen/Documents/useful_julia_scripts/save_by_regex.jl")

file = matopen("/mnt/Lab_data/MVCO/FCB/Syn_and_MVCO_packaged_data/syndata_current.mat")
varnames = names(file)


#assign variables one by one:
#mu_avg=read(file, "mu_avg") # note that this does NOT introduce a variable ``varname`` into scope

# @eval $(Symbol(string)) = $value or @eval $(Symbol(string)) = value.

## Success! Finally!
#to_add=[] #does this need to allocated first?
for var in varnames

    if occursin("time",var) #need to deal with time differently as these are matlab dates/times
        @show var

        ex = :(temp = read(file,$var)) #import variable as temporary
        eval(ex)

        global to_add = temp .- temp[1,1] #reset to 0 days first entry should be Jan-1-03!
        @show typeof(to_add)
        v=Symbol(var)
        ex=:($v = DateTime(2003,1,1) + map(x->Dates.Day(x), to_add)) #sweet! Success!
        eval(ex)

    else
        v=Symbol(var)
        ex = :($v = read(file,$var))
        #@show ex
        eval(ex)
    end
end

close(file)
##

#saving with regex! thanks to Dehann!
save("/home/kristen/Documents/V6V8_analysis/analysis_products/mvco_syndata.jld2", [r"syn",r"mu",r"net",r"loss",r"PE",r"SSC",r"vol",r"CHL",r"allgrowthrates",r"allMR",r"modelallmatdate",r"allmatdat",r"FCB"])
#*syn* *mu* *net* *loss* *PE* *SSC* *vol* *CHL* allgrowthrates allMR modelallmatdate allmatdat* FCBnum
#older, default save:
#@save("/home/kristen/Documents/V6V8_analysis/scripts/analysis_products/mvco_syndata.jld2")



## and the environmental data!

file = matopen("/mnt/Lab_data/MVCO/FCB/Syn_and_MVCO_packaged_data/mvco_envdata_current.mat")
varnames = names(file)

#to_add=[] #does this need to allocated first?
for var in varnames

    if occursin("time",var)#need to deal with time differently as these are matlab dates/times
        @show var

        ex = :(temp = read(file,$var)) #import variable as temporary
        eval(ex)
        global to_add = temp .- temp[1,1] #reset to 0 days first entry should be Jan-1-03!
        # @show typeof(to_add)
        v=Symbol(var)

        uu=findall(map(x->isnan(x),temp))
        if isempty(uu)

            to_add = map(x->Dates.Day(x), to_add)
            @show typeof(to_add)
            ex=:($v = DateTime(2003,1,1) + to_add) #sweet! Success!
            @show ex
            eval(ex)

        end
    else
        v=Symbol(var)
        ex = :($v = read(file,$var))
        #@show ex
        eval(ex)
    end
end


close(file)
##
save("/home/kristen/Documents/V6V8_analysis/analysis_products/mvco_envdata.jld2",[r"light", r"Tcorr", r"corr_years", r"Tbeam_corr", r"dy_", r"_month nutyears", r"nutyrdy", r"time_nut", r"SiOH", r"NH4", r"NO3", r"PO4"])
#*light* *Tcorr* corr_years Tbeam_corr dy_*  *_month nutyears nutyrdy time_nut* SiOH_wk NH4_wk NO3_wk PO4_wk SiOH_*_mn NH4_*_mn NO3_*_mn PO4_*_mn'
#@save "/home/kristen/Documents/V6V8_analysis/scripts/analysis_products/mvco_envdata.jld2"

#EVAL AND EXPRESSIONS:
# var="mu_avg"
# v=Symbol(var)
# ex = :($v = read(file,$var))
# eval(ex)

# a=1
# b=3
# ex=:(a+b)
# eval(ex)

#we need to work with time....so, easy way?
#731582.0 = Jan 1 2003 > work it!
# to_add = time_syn .- time_syn[1,1]
#
# map(x->Dates.Day(x), to_add)
# a= DateTime(2003,1,1) + map(x->Dates.Day(x), to_add) #sweet! Success!

## to match mvco event numbers:
file = matopen("/home/kristen/Documents/seasons_of_syn_paper/matlab_scripts/nut_data_reps.mat")
varnames = names(file)
for var in varnames
    v=Symbol(var)
    ex = :($v = read(file,$var))
    #@show ex
    eval(ex)
end
close(file)
##

@save "/home/kristen/Documents/V6V8_analysis/scripts/analysis_products/mvco_nut_data.jld2" MVCO_nut_reps header_nut




### plotting....
# gr()
# scatter(ydmu, mu_avg)
# xaxis!("Year day",(0, 366),0:100:360,font(12, "Times"))
# # xlabel = "my label",
# #    xlims = (0,10),
# #    xticks = 0:0.5:10,
# #    xscale = :log,
# #    xflip = true,
# #    xtickfont = font(20, "Courier"))
#
# plot(time_syn,daily_syn, yscale = :log10, size=(600,200), framestyle =:box, legend=:outertopright)
