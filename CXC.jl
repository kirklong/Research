#Routines for reading in various Chandra Observatory Data Products

using DelimitedFiles, CSV, FITSIO, DataFrames, Statistics

##############################################################################
function cxctimein(infile)

#hri names include a _1rh delimiter
#hri photons: time xdet ydet
#_1rh ends srcname and next . ends obs name

#pspc names include an _rp or _wp
#pspc files header: ntimes tstart tend ?????
#pspc photons: time xdet ydet pha
#_rp or _wp delimiter ends srcname and next _ ends obs name
#srcb=CSV.read(infile; header=false,datarow=8,delim=' ',ignorerepeated=true)

#println(srcb[1:10])
  #return srcb[:Column3]
  #println(size(srcb)[1])

  fsize=filesize(infile)
  if fsize < 500
    return([])
  end

  srcb=readdlm(infile, skipstart=7)
   if size(srcb)[1] > 0   #test for empty list
     #println("Hi there\n")
     #return srcb[:Column3]
     return srcb[1:end,2]
   else
      return([])
   end

end

##############################################################################
#Creates a delimited file containing information on all CXC source in field.   
#DJM, Mar 2019
function cxc_srcinfo(obsid,catfile)

    mydf = DataFrame()


    f = FITS(catfile)
    header = read_header(f[2])
    ra=read(f[2], "RA"; case_sensitive=true)
    dec=read(f[2], "DEC"; case_sensitive=true)
    raerr=read(f[2], "RA_ERR"; case_sensitive=true)
    decerr=read(f[2], "DEC_ERR"; case_sensitive=true)
    netcts=read(f[2], "NET_COUNTS"; case_sensitive=true)
    netctserr=read(f[2], "NET_COUNTS_ERR"; case_sensitive=true)
    bkgcts=read(f[2], "BKG_COUNTS"; case_sensitive=true)
    bkgctserr=read(f[2], "BKG_COUNTS_ERR"; case_sensitive=true)
    exposure=read(f[2], "EXPTIME"; case_sensitive=true)
    psfsize=read(f[2], "PSF_SIZE"; case_sensitive=true)
    shape  =read(f[2], "SHAPE"; case_sensitive=true)
    psfratio=read(f[2], "PSFRATIO"; case_sensitive=true)

    exptime=read(f[2], "EXPTIME"; case_sensitive=true)
    curmjdref =header["MJDREF"] 
    curexptime=header["ONTIME"] 
    curtstart =header["TSTART"]  #seconds since ref MJD
    curtstop  =header["TSTOP"]  #seconds since ref MJD
    curmission=header["TELESCOP"]  #seconds since ref MJD
    curinstr  =header["INSTRUME"]  #seconds since ref MJD
    curtimeres=header["TIMEDEL"]  #seconds since ref MJD

    #dateend=header["DATE-END"] 
    #exptime=header["DATE-OBS"] 

    tstart = curmjdref + curtstart/60.0/60.0/24.0
    tstop  = curmjdref + curtstop/60.0/60.0/24.0

    nlen = length(ra)
    num=collect(1:nlen)
    sdate=[ tstart for n in 1:nlen]
    edate=[ tstop for n in 1:nlen]
    expo=[ exptime for n in 1:nlen]
    mission=[ curmission for n in 1:nlen]
    instr=[ curinstr for n in 1:nlen]
    timeres=[ curtimeres for n in 1:nlen]
    mjdref =[ curmjdref for n in 1:nlen]
    mission =[ curmission for n in 1:nlen]
    instr   =[ curinstr for n in 1:nlen]

    mydf.mission=mission
    mydf.instr=instr
    mydf.num=num 
    mydf.mjdref=mjdref
    mydf.sdate=sdate
    mydf.edate=edate
    mydf.expo=expo
    mydf.timeres=timeres
    mydf.ra = ra
    mydf.raerr = raerr
    mydf.dec = dec 
    mydf.decerr = decerr 
    mydf.netcts = netcts
    mydf.netctserr = netctserr
    mydf.bkgcts = bkgcts
    mydf.bkgctserr = bkgctserr
    mydf.exptime = exptime
    mydf.psfsize = psfsize
    mydf.shape = shape
    mydf.psfratio = psfratio

    close(f)

    outname = string(obsid,"_srclist.csv")
    CSV.write(outname,mydf,delim=',')
    println(mydf)
    return(mydf)

end

##############################################################################
#Creates a file containing all significant powers from all fft files   
#---DJM, Mar 2019
#---function cxc_fftinfo(obsid,datadir,srcinfo,suff;plimit=10.0,fperiods=[])
#---should consider adding obsid to srcinfo file rather than as parameter

function cxc_fftinfo(obsid,datadir,srcinfo,suff;plimit=10.0)

    cd(datadir)
    srcf=string(obsid,"_srclist.csv")
    srcdf=CSV.read(srcf)

    fnames = readdir(datadir)
    tnames = fnames[occursin.(suff,fnames)]  #get files with suffix
    nnames= length(tnames)
    nsuff = length(suff)

    allobs = Array{Float32}[]
    allnum = Array{Float32}[]
    allsdate = Array{Float32}[]
    alledate = Array{Float32}[]
    allobsdur =Array{Float32}[]
    allexpo =Array{Float32}[]
    alltres =Array{Float32}[]
    allra = Array{Float32}[]
    alldec = Array{Float32}[]
    allfreq = Array{Float32}[]
    allpow = Array{Float32}[]
    #dfall = DataFrame()
    allinstr = Array()
    allmission  = Array()

    #for i in 1:nnames    #loop over fft files 
    for curfile in tnames    #loop over fft files 

        #curroot = tnames[i][1:end-nsuff]
        curroot = curfile[1:end-nsuff]
        #fname = tnames[i]
        fname = curfile

        #read in fft information
        curfft = CSV.read(fname,delim="\t",header=false)
        #normalize fft
        meanfft = Statistics.mean(curfft[2][2:end-1]) #avoid dc level
        avefft=curfft[2]/meanfft
        npows=length(avefft)
       
        #find normalized powers > plimit
        #avefft[avefft .>= plimit] 
        mid = convert(Int64,floor(npows/2.)) + 1
        highp=findall(x -> x>=plimit, avefft[2:mid])
        highp = highp .+ 1
        nhighp = length(highp)

#skip the rest if no high powers
        if nhighp == 0 
           continue
        end

#get vectors of high powers and associated frequencies
        freqarr=curfft[1][highp]
        powarr=avefft[highp]

        ##println(fname)
        #fend=findlast(".",fname)[1] - 1
        #filehd=fname[1:fend]

#tokenize the fft file name
        mystr=replace(fname,r"[_.]" => " ")
        smatch = eachmatch(r"\w+",mystr)
        smatcharr=collect(smatch)

#replicate obsid and source number as integers from standard fft file name
        obsid=parse(Int64,smatcharr[1].match)
        srcnum=parse(Int64,smatcharr[3].match)

#create vectors giving the obsid and source number
        obsarr = fill(obsid,nhighp)
        numarr = fill(srcnum,nhighp)


#MUST BE CAREFUL HERE - get the source information not from the
#assuming the index of the dataframe is the same as the source number,
#but instead get the dataframe row where the :num field is the same
#as srcnum.  Necessary since not every .times file has an fft

###BUT srcdf should have an entry for each source, and be ordered
#by source number since it comes from the fits file, so can just
#use the srcnum row
        ##strtarr = fill(srcdf[:start][i],nhighp)
        #strtarr = fill(srcdf[:sdate][i],nhighp)
        #expoarr = fill(srcdf[:expo][i],nhighp)
        #raarr = fill(srcdf[:ra][i],nhighp)
        #decarr = fill(srcdf[:dec][i],nhighp)

        strtarr = fill(srcdf[:sdate][srcnum],nhighp)
        expoarr = fill(srcdf[:expo][srcnum],nhighp)
        raarr = fill(srcdf[:ra][srcnum],nhighp)
        decarr = fill(srcdf[:dec][srcnum],nhighp)

#    println(powarr,"\n")
        curdf = DataFrame(obsid=obsid,srcnum=srcnum,startdate=strtarr,exptime=expoarr,ra=raarr,dec=decarr,freq=freqarr,power=powarr)

        #append!(dfall, curdf)
#        println(curdf,"\n")

        append!(allobs,obsarr)  #write the new summed flux
        append!(allnum, numarr)
        append!(allstrt, strtarr)
        append!(allexpo, expoarr)
        append!(allra, raarr)
        append!(alldec, decarr)
        append!(allfreq, freqarr)
        append!(allpow, powarr)

        #supplement the output arrays

    end

    dfall = DataFrame(obsid=allobs,srcnum=allnum,startdate=allstrt,exptime=allexpo,ra=allra,dec=alldec,freq=allfreq,power=allpow)
    fout = string(obsid,"_fftinfo.csv") 
    CSV.write(fout,dfall)

end
