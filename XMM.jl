#Routines for processing various XMM SAS Data Products

using DelimitedFiles, CSV, FITSIO, DataFrames, Statistics

##############################################################################
#Creates a set of lightcurves for detected sources given an XMM observation
#directory
#
#expects full path to directory with no trailing / and assumes SAS setup
#including definition of SAS_CCF and SAS_ODF! could do this in script 
#including seeing if  odfingest and/or cifbuild need running but do later
#
#---DJM, Sep 2019

function procxmm(datadir::String)

    ppsdir=string(datadir,"/PPS")
    odfdir=string(datadir,"/ODF")
    ppsnames = readdir(ppsdir)

#Assume for now the stupid XMM .FTZ files have been unzipped as .FIT files
    f8str="IMAGE_8000.FIT"
    imnames = fnames[occursin.(f8str,ppsnames)]  #get files with suffix
    nimages= length(imnames)

#now have the set of broad-spectrum image names, use to ID the proper event
#files and extract the photon times for each source
    for i in 1:nimages
        fname = imnames[i]
        obstrt=findlast(f8str,fname)[1] - 6
        obsname = fname[obstrt:obstrt+5]

#now have obsname (like which mos or pn observation) corresponds to this
#image to further refine commands and locate the event file
        curobsfiles=ppsnames[occursin.(obsname,ppsnames)]
##again making the file search two steps, can probably do a single step
#if I use regular expressions
        curevtfile=curobsfile[occursin.("EVL",curobsfile)]


#example MOS filtering command
#evselect table=mos1.fits withfilteredset=yes \
#$   $ expression='(PATTERN $<=$ 12)&&(PI in [200:12000])&&#XMMEA_EM' \
#$   $ filteredset=mos1_filt.fits filtertype=expression keepfilteroutput=yes \
#$   $ updateexposure=yes filterexposure=yes




#example MOS lightcurve generation command
#evselect table=mos1.fits withrateset=yes rateset=mos1_ltcrv.fits \
#$   $ maketimecolumn=yes timecolumn=TIME timebinsize=100 makeratecolumn=yes [] fv mos1_ltcrv.fits &






#in particular, we can use the current image file to get a filtered gti file

        fend=findlast(,fname)[1] - 1
        filehd=fname[1:fend]
    #println(filehd,"\n")

        println(fname," ",size(intimes)[1])
        if size(intimes)[1] <= 1
            continue
        end

        fftoutstr=string(filehd,".fft")
        fftpltstr=string(filehd,".png")
        fftjldstr=string(filehd,".jld")

        writedlm(fftoutstr,hcat(myfreqs,mypow))
        #save(fftjldstr, "freqs", myfreqs, "powers", mypow)

        pyplot()
        titlestr=string("FFT of ",filehd)
        myplot=plot(myfreqs[2:end],mypow[2:end]/avepow,title=titlestr,xlabel="freq",ylabel="Power")
        png(myplot,filehd)
    end
end
