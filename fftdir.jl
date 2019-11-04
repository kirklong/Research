#perform FFTs on times file in a particular directory.
using DelimitedFiles, Plots, Statistics, JLD, DataFrames

include("/home/djm/Programming/Julia/DJM.jl")

#function fftdir(datadir,suff=suff,inst="",mmult=mmult,srchlim=srchlim,nlimit=50,by2=by2,noplot=noplot)
function fftdir(datadir::String,suff::String,delt::Float64)

    fnames = readdir(datadir)
    tnames = fnames[occursin.(suff,fnames)]  #get files with suffix
    nnames= length(tnames)
    #delt=3.24104

    for i in 1:nnames
        fname = string(datadir,tnames[i])
        #println(fname)
        fend=findlast(".",fname)[1] - 1
        filehd=fname[1:fend]
    #println(filehd,"\n")

        intimes=DJM.cxctimein(fname)
        println(fname," ",size(intimes)[1])
        if size(intimes)[1] <= 1
            continue
        end

        myfft = DJM.ffttimes(intimes,delt)
        mypow=abs2.(myfft)
        avepow=mean(mypow[2:end])

        nfreqs=length(myfft)
        myfreqs = DJM.freqarr(delt,nfreqs)

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
# pro fft_dir_124,datadir,limit,suff=suff,mmult=mmult,inst=inst,srchlim=srchlim,nlimit=nlimit,by2=by2,noplot=noplot
#
# ;t1=systime(/seconds)
# ;set the amount of memory to use and adjust the maximum size
# ;of the fft array.  Will need to print a report on how much
# ;of the data is used, or loop over performing independent ffts
# ;until the entire datafile is checked
# if not(keyword_set(mmult)) then mmult=1
# if not(keyword_set(nlimit)) then nlimit=250
# if not(keyword_set(inst)) then inst="chandra"
# if not(keyword_set(by2)) then by2=0
# if not(keyword_set(noplot)) then noplot=0
# if not(keyword_set(srchlim)) then srchlim=0
# if not(keyword_set(suff)) then suff="times_bc"

#
# ;limit the array to less than the available memory by allowing
# ;the user to but it the desired size in 128 Megabyte units
# memsize=mmult*(2L^27)
# binlimit=long(memsize/8L) ;max size of complex array
# ;print, "Memory used: ",memsize, " # of complex bins: ",binlimit
#
# ;process the name of the file to get the root directory
# thedir=datadir
# srcfiles=file_search(thedir + "/source_*." + suff)
# fitsfile=file_search(thedir + "/*evt2_bary.fits")
# srcfile=file_search(thedir + "/*src2.fits")
#
# ;make certain you can read the times file using the
# ;appropriate template
# if(n_elements(chandra_template) eq 0) then restore, "chandra_template.sav"
# nfiles=n_elements(srcfiles)
# ntimes=fltarr(nfiles)
#
# ; get the necessary fits file info, keeping in mind that the
# ; appropriate event times have already been extracted
# inhead=headfits(fitsfile[0],EXTEN=1)
# dateend=fxpar(inhead,"DATE-END")
# dateobs=fxpar(inhead,"DATE-OBS")
# timeres=float(fxpar(inhead,"TIMEDEL"))
# tstart = fxpar(inhead,"TSTART")
# tstop = fxpar(inhead,"TSTOP")
# livetime=tstop-tstart
# exposure=sxpar(inhead,"EXPOSURE")
# print, dateend, timeres, tstop - tstart
#
# ;autosrc=mrdfits(datadir+"/auto_src.fits",1,srchd)
# sources=mrdfits(srcfile[0],1,srchd)
# nfiles=n_elements(sources)
#
# nbins=exposure/timeres
# ;print, nbins, livetime, timeres
# nind=findgen(30)
# binby2=2L^nind
# above=where(binby2 ge nbins)
# print, above,nbins
# maxbins=binlimit
# if (above[0] ge 0) then maxbins = binby2[above[0]]
# if( (maxbins gt binlimit) or (above[0] eq -1) ) then maxbins=binlimit
# cfrac = float(maxbins)/float(nbins)     ;fraction of observation covered
#
# set_plot, 'z'
#
# highp1=fltarr(nfiles)
# highf1=fltarr(nfiles)
# highp2=fltarr(nfiles)
# highf2=fltarr(nfiles)
# highp4=fltarr(nfiles)
# highf4=fltarr(nfiles)
#
# rootdir = file_dirname(srcfiles[0])
# endid=strpos(rootdir,"primary")
# strtid=strpos(rootdir,path_sep(),endid-2,/reverse_search)
# rootdir=strmid(rootdir,strtid+1,endid-strtid-2)
# maxsrch=(maxbins/2) - srchlim
# maxfreq=1.0/timeres/2.0
# minsrch=srchlim
# sind=strpos(srcfiles,"_")
# eind=strpos(srcfiles,"."+suff)
#
# ;
# for i = 0,nfiles-1 do begin
#
#     times=read_ascii(srcfiles[i],template=chandra_template,count=reccts)
#     if(reccts eq 0) then continue
#     nevents=long(n_elements(times.t))
#     if(nevents lt nlimit) then continue
#     tactual=times.t - times.t[0]
#     under1=strpos(srcfiles[i],"_")
#     dot1=strpos(srcfiles[i],".")
#
# ;get the number from the times file format note that srcnum != i
#     ;srcnum=strmid(srcfiles[i],under1+1,dot1-under1-1)
#     srcnum=fix(strmid(srcfiles[i],sind[i]+1,eind[i]-sind[i]-1))
#     ntimes[srcnum]=nevents
#     bintimes=complexarr(maxbins)
#     ntimed=0L
#
#     ;bintimes=histogram(tactual,binsize=timeres,nbins=maxbins)
#     for j = 0L, nevents-1 do begin
#        tint=long(tactual[j]/timeres)
#        if(tint lt maxbins) then bintimes[tint]=bintimes[tint] + complex(1,0)
#        ntimed=ntimed+1
#     endfor
#
#     bintimes=temporary(fft(bintimes,/overwrite))
#     bintimes=abs(temporary(bintimes))
#     bintimes=(temporary(bintimes))^2
#     undefine, tactual
#     undefine, times
#
#     pave = (mean(bintimes))
#     highp1[srcnum] =max(bintimes[minsrch:maxsrch],maxbin)/pave
#     highf1[srcnum] = maxbin/timeres/maxbins
#     ;highp[i] =max(bintimes[minsrch:maxsrch],maxbin)/pave
#     ;highf[i] = maxbin/timeres/maxbins
#     ;memstr=memory()
#     ;print, "Memory after calculating highp: ",memstr
#     delf=1.0/timeres
#
#     titlestr="Inst: " + inst + " Obsid: " + rootdir + " Src #: " + string(srcnum)
#
#     ;substr=" Ra: " + string(autosrc[srcnum-1].ra) + " Dec: " + string(autosrc[srcnum-1].dec)
#     substr=" Ra: " + string(sources[srcnum].ra) + " Dec: " + string(sources[srcnum].dec)
#     ;substr=" Ra: " + string(autosrc[srcnum-1].ra) + " Dec: " + string(autosrc[srcnum-1].dec)
#     ;ymax=1.1*highp[i]
#     ymax=1.1*highp1[srcnum]
#
#     plot, bintimes[0:maxbins/2]/pave,$
#     xtitle="Bin #",ytitle="Normalized Power",title=titlestr,xticks=8,$
#     xrange=[0,maxbins/2-1],xstyle=1,yrange=[0,ymax],ystyle=1
#     xyouts, 0.6,0.9,"Start Time: "+dateobs,/normal
#     xyouts, 0.6,0.85,"End Time: "+dateend,/normal
#     xyouts, 0.6,0.80,"Time Resolution: " + string(timeres),/normal
#     xyouts, 0.6, 0.75,"RA: " + string(sources[srcnum-1].ra),/normal
#     xyouts, 0.6, 0.70,"DEC: " + string(sources[srcnum-1].dec),/normal
#
#     nyquist = 0.5/timeres
#     xyouts, 0.6, 0.65,"Nyquist: " + string(nyquist),/normal
#     xyouts, 0.6, 0.60,"# of Events: " + string(nevents),/normal
#     ;memstr=memory()
#     ;print, "Memory after plotting: ",memstr
#
#     basefile=file_basename(srcfiles[i],".times_bc")
#     write_jpeg, datadir + basefile + ".jpg",tvrd(),quality=95
#
#     ndot=strpos(srcfiles[i],".")
#     tname=srcfiles[i]
#     basesrc=strmid(tname,0,ndot)
#     savpow=basesrc+"_powers124.sav"
#
#     nhigh1=where(bintimes ge limit*pave,ncount)
#     if(ncount gt 0) then begin
#         powers1=bintimes[nhigh1]/pave
# ;        save, timeres,pave,nhigh,powers,filename=savpow
#     endif
#
#     nhigh2=lonarr(1)
#     powers2=fltarr(1)
#     k=0L
#     for k = long(srchlim),long(maxbins/4) do begin
#        sum2=(bintimes[k]+bintimes[2*k])/pave
#        if(sum2 gt 1.41421*limit) then begin
#          nhigh2=[nhigh2,k]
#          powers2=[powers2,sum2]
#        endif
#     endfor
#     highp2[srcnum] =max(powers2,maxbin)
#     highf2[srcnum] = nhigh2[maxbin]/timeres/maxbins
# ;    bintimes[0:maxbins/4-1]=bintimes[0:maxbins/4-1]+bintimes[maxbins/4:maxbins/2-1]
# ;    nhigh2=where(bintimes[0:maxbins/4-1] ge 1.414*limit*pave,ncount)
# ;    if(ncount gt 0) then begin
# ;        powers2=bintimes[nhigh2]/pave
# ;        save, timeres,pave,nhigh,powers,filename=savpow
# ;    endif
#     nhigh4=lonarr(1)
#     powers4=fltarr(1)
#     for k = long(srchlim),long(maxbins/8) do begin
#        sum4=(bintimes[k]+bintimes[2*k]+bintimes[3*k]+bintimes[4*k])/pave
#        if(sum4 gt 2.0*limit) then begin
#          nhigh4=[nhigh4,k]
#          powers4=[powers4,sum4]
#        endif
#     endfor
# ;    highp4[srcnum] =max(npowers4,maxbin)
# ;    highf4[srcnum] = maxbin/timeres/maxbins
#     highp4[srcnum] =max(powers4,maxbin)
#     highf4[srcnum] = nhigh4[maxbin]/timeres/maxbins
# ;    bintimes[0:maxbins/8-1]=bintimes[0:maxbins/8-1]+bintimes[maxbins/8:maxbins/4-1]
# ;    nhigh4=where(bintimes[0:maxbins/8-1] ge 2.*limit*pave,ncount)
# ;    if(ncount gt 0) then begin
# ;        powers4=bintimes[nhigh4]/pave
# ;        save, timeres,pave,nhigh,powers,filename=savpow
# ;    endif
#
#     save, timeres,pave,nhigh1,nhigh2,nhigh4,powers1,powers2,powers4,filename=savpow
#
#     undefine, bintimes
#
# ;    memstr=memory()
# ;    print, "Memory after end of source file: ",memstr
#     ;t2=systime(/seconds)
#     ;print, "Elapsed time: ",t2-t1
#     undefine, bintimes
#
# endfor
#
# savsum=datadir+"fft_summary_124.sav"
# save,srcfiles,ntimes,highp1,highf1,highp2,highf2,highp4,highf4,tstart,tstop,exposure,timeres,maxbins,cfrac,inhead,autosrc, filename=savsum
#
# end
