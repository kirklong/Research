#Routines for reading in various (mostly delimited) astronomical data files

##############################################################################
#Load all data from a HEASARC tdat file into a Julia Dataframe.
#requires initial open to find params then uses CSV read to dataframe
#DJM, Jan 2019

using DelimitedFiles, CSV, FITSIO, DataFrames, Statistics

function tdat2df(tdatin,fields=[])
    #tdatin="/data/heasarc/heasarc_cxoxassist.tdat"
    #tdat = readlines(tdatin)
    fin = open(tdatin)

    nline=1
    namearr=[]
    endhd=0
    for inln in eachline(fin)
#find the tag that denotes the start of data
        if occursin(r"^<DATA>",inln)
          #println(inln)
          endhd = nline
          break
        end
        nline=nline+1

#find the line containing number of data rows and parse
        if occursin(r"TOTAL ROWS",inln)
          totstrt = findfirst("TOTAL ROWS",inln)[1]
          nrows=parse(Int64,inln[totstrt+12:end])
          #println("Number of Rows: ",nrows)
        end

#Parse the line containing the column names
        if occursin(r"^line\[1\]",inln)
          eqrange = findfirst(" = ",inln)[end]
          namearr = split(inln[eqrange+1:end])
          #println(namearr[1])
        end

    end
    close(fin)

#tdat files have trailing delimiter, for now pad with empty column
    push!(namearr,"null0")

#Went through header and got main info, now read data
#Not elegant, but it works
    tdatdf = CSV.read(tdatin,delim='|',header=namearr,datarow=endhd+1,footerskip=1)
    return(tdatdf)

end
