using Dates
#define various parameter columns from the MATLAB created csv files
#and Swift/BAT definitions.  Note: reduced form csv files from newcsv
#are more self-descriptive because they have headers
elo8=[14,20,24,35,50,75,100,150];
ehi8=[20,24,35,50,75,100,150,195];
elo4=[14,24,50,100];
ehi4=[24,50,100,195];
elo=elo8
ehi=ehi8
snrind=collect(155:162);
fluxind=collect(24:31);
fluxerrind=collect(32:39);
tstartind=4;
tstopind=5;
expoind=6;
maxsnrind=164;
totsnrind=163;
det2area=0.16;  #BAT area in cm^2 per detector cell; divide rates by this
secperday=86400.0;
dayperyear=365.25  #Days per Julian Year
tzero=DateTime(2001, 1, 1, 0, 0, 0);  #BAT Swift Mission Elapsed Time is Seconds since 1 Jan 2001
jdzero=Dates.datetime2julian(tzero);  #JD of zero time
mjdzero=jdzero-2400000.5; #Modified Julian Date Definition

#read in from FITS file from directory "root", make a new dataframe and
#write a lovely new jld file.  Only saves the columns noted below.
#to save all, use "fits2csv"
##############################################################################
#root: directory containing the fits file
#oldfits: name of the fits file
function fits2dict(root,oldfits,out=0)
    #print(oldfits)
    secperday=86400.0;

    inname = oldfits

    #for now, keeping all data in the current directory
    #csvin = readtable(oldcsv,header=false)
    fin = FITS(string(root,"/",oldfits))
    #println(fin)
    #curexpo = read(f[2], "EXPOSURE")

    srcname = read(fin[2], "NAME")
    obsid = read(fin[2],"OBS_ID")

    #maxsnr=csvin[:,maxsnrind]
    ra = read(fin[2],"RA_OBJ")
    dec = read(fin[2],"DEC_OBJ")
    #maxsnr=csvin[:,maxsnrind]
    rapnt = read(fin[2],"RA_PNT")
    decpnt = read(fin[2],"DEC_PNT")

    #tstart=mjdzero+csvin[:,tstartind]/secperday
    tstartin = read(fin[2], "TIME")
    tstart = mjdzero.+tstartin/secperday
    tstopin=read(fin[2], "TIME_STOP")
    tstop= mjdzero.+tstopin/secperday
    #texpo=csvin[:,expoind]
    texpo = read(fin[2],"EXPOSURE")

    #flux=csvin[:,fluxind]
    rate = read(fin[2],"RATE")
    #fluxerr=csvin[:,fluxerrind]
    rateerr = read(fin[2],"RATE_ERR")

    #flux=csvin[:,fluxind]
    bkg = read(fin[2],"BKG")
    #fluxerr=csv in[:,fluxerrind]
    bkgerr = read(fin[2],"BKG_ERR")


    #flux=csvin[:,fluxind]
    chi2 = read(fin[2],"CHI2")
    #fluxerr=csvin[:,fluxerrind]
    chi2dof = read(fin[2],"DOF")
    detect_status = read(fin[2],"DETECT_STATUS")

    #snr=csvin[:,snrind]
    vectsnr = read(fin[2],"VECTSNR")
    #maximum of vect/tot snr
    maxsnr = read(fin[2],"SNR")
    #S/N for all bands combined
    totsnr = read(fin[2],"TOTSNR")
    ngood = findall(x->x>0., totsnr)  #for now, fin d positive totsnr

    #Number of enabled detectors
    ngoodpix = read(fin[2],"NGOODPIX")
    #Partial Coding Fraction, Fraction of pixels active
    pcodefr=read(fin[2],"PCODEFR")

    #Processing flags, should see which need to be monitored
    backapp = read(fin[2],"BACKAPP")
    acolapp = read(fin[2],"ACOLAPP")
    pcodeapp = read(fin[2],"PCODEAPP")
    ffapp = read(fin[2],"FFAPP")
    ngpixapp = read(fin[2],"NGPIXAPP")
    bdistapp = read(fin[2],"BDISTAPP")
    mskwtapp = read(fin[2],"MSKWTAPP")
    occapp = read(fin[2],"OCCAPP")

    close(fin)

    #HERE IS WHERE WE MIGHT PUT DATA FILTERING BASED UPON VALUES

   df = Dict(:srcname=>srcname[1],:ra=>ra[1],:dec=>dec[1],:obsid=>obsid,:rapnt=>rapnt,:decpnt=>decpnt,
        :tstart=>tstart,:tstop=>tstop,:texpo=>texpo,:rate=>rate,:rateerr=>rateerr,:bkg=>bkg,bkgerr=>bkgerr,
        :chi2=>chi2,:chi2dof=>chi2dof,:vectsnr=>vectsnr,:maxsnr=>maxsnr,:totsnr=>totsnr,:ngoodpix=>ngoodpix,
        :pcodefr=>pcodefr,:backapp=>backapp,:acolapp=>acolapp,:pcodeapp=>pcodeapp,ffapp=>ffapp,:ngpixapp=>ngpixapp,
        :bdistapp=>bdistapp,:mskwtapp=>mskwtapp,:occapp=>occapp)   #This dataframe holds only selected columns from original

    #this section checks for a non-zero third argument meaning a jld file is written
    if out != 0
        strend = findfirst(".fits",oldfits)
        #println(strend)
        short = oldfits[1:strend[1]-1]
        save(string(root,"/",short,".jld"),"batdat",df)
    end

    return df
    end

    ##############################################################################
    #rebin pointing flux, fluxerror, according to some time length in days
#using the dictionaries defined for BAT/Swift per Lister analysis
#
function binflux(indict,bindays,lowind,hiind)

###INPUTS
    #indict: format as per fits2jld
    #bindays: The number of days in a time bin, i
    #lowind: the index of the lowest energy bin
    #hiind: the index of the highest energy bin

###OUTPUT:
    #A dataframe containing the flux, flux errors, time spans, and number of segments

    flux=[]
    fluxerr=[]
    exposure=[]
    midtime=[]    #note midtime is for splitting the time-span, not mid- exposure
    npoints=[]

    #unlikely the input dictionary is sorted in any way, collect the exposure bin indices (by start time)
    sdate=minimum(indict[:tstart])
    binnum=trunc.(Int,(indict[:tstart].-sdate) / bindays) #get time bins
    maxbin=maximum(binnum)

    for i = 0:maxbin
        good=findall(x->x==i, binnum)
        ngood=size(good)[1]
        if ngood==0
            continue
        end
    #
        fluxtmp=0
        fluxerrtmp=0
        expotmp=0

        ###Note that bat survey data is Gaussian statistics, not Poisson
        for j in good  #convert things to counts for calcs before going back to flux
            fluxtmp+=sum(Array(indict[:rate][lowind:hiind,j]))*indict[:texpo][j]
            #fluxerrtmp+=sum(Array(trial[j,2+lowind:2+hiind]))*trial[j,2] #energy columns start at 3
            #add the sum of the squared errors times exposure
            fluxerrtmp+=sum(Array(indict[:rateerr][lowind:hiind,j]).*Array(indict[:rateerr][lowind:hiind,j]))*indict[:texpo][j]^2
            #sqrt(sum(Array(trial[4,2+1:2+8]).*Array(trial[4,2+1:2+8])))
            expotmp+=indict[:texpo][j]
        end

        push!(flux,fluxtmp/expotmp)  #write the new summed flux
        push!(fluxerr,sqrt(fluxerrtmp)/expotmp)  #write the new summed flux
        push!(exposure,expotmp)  #write the new summed flux
        push!(midtime,sdate+bindays*(i+0.5))  #write the new summed flux
        push!(npoints,ngood)  #write the new summed flux

    end

    return DataFrame(flux=flux,fluxe=fluxerr,expo=exposure,tmid=midtime,npoints=npoints)

end

################################################################################################

###creates a histogram from "data" based upon min/max/nbins construct.
#Can probably be replaced with a package function like fit from StatsBase

function myhist(data, min, max, nbins)
  N = length(data)             # How many elements in the input vector 'data' ?
  delta = (max-min)/nbins      # Bin size is inferred here from the maximal, minimal, and bin number
    out = zeros(nbins)           # Let's initialize ty6yhy7 nujihe output data structures for the bin count
  bin = zeros(nbins)           # and for the bin centres...

  start = min                  # Left edge
  for k=1:nbins
    stop   = start + delta   # Right edge
    out[k] = length(findall((data .>= start) .& (data .< stop))) # Count how many elements are between left and right
    bin[k] = start + delta/2. # Centre of the bin
    start  = stop            # New left edge
   end
   return out, bin
  end
  #########################################################################################
  
