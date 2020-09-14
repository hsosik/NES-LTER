
#usage
#varinfo(r"oligo").content[1].rows[2:end][1][1]

#Dehann Fourie, Massachusetts Institute of Technology, Spring 2020

function varinfo_dict(rgx::Regex, dict::Dict=Dict{String, Any}())
  # get Markdown tables of variables matching expression
  vartables = varinfo(rgx).content[1].rows[2:end]
  # extract variable names as strings
  varnames = map(x->x[1], vartables)

  # combine names with values in a dictionary for JLD2/FileIO interface
  for v in varnames
    dict[v] = getfield(Main, Symbol(v))
  end

  return dict
end

#usage:
#varinfo_dict(r"time")

import FileIO: save

function save(filename::AbstractString, rgx::Regex)
  save(filename, varinfo_dict(rgx))
end

function save(filename::AbstractString, rgxs::Union{Vector{Regex}, NTuple{N, Regex}}) where N
  vardict = Dict{String, Any}()
  for rx in rgxs
    varinfo_dict(rx, vardict)
  end

  save(filename, vardict)
end


#usage:
#save(<file>, [r"PE",r"mu",r"light"])
#save(<file>, (r"PE",r"mu",r"light"))
