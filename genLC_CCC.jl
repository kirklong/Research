#!/usr/bin/env julia
using Plots, CSV, HypothesisTests, Printf, DataFrames

function genCCC(a,b) #implement circular correlation coefficient from Fisher & Lee (1983)
    #reference: https://cnx.org/contents/rQ-KKFjj@3/Circular-Correlation-Coefficient
    n=length(a)
    num,den1,den2=0,0,0
    for i=1:(n-1)
        for j=(i+1):n
            num+=sin(a[i]-a[j])*sin(b[i]-b[j])
            den1+=sin(a[i]-a[j])^2
            den2+=sin(b[i]-b[j])^2
        end
    end
    den=sqrt(den1*den2)
    r=num/den
    return r
end

function getFloatTimes(dataList)
    timeList=zeros(length(dataList)-7) #-7 because 1st 7 are bs
    for i=8:length(dataList) #real data starts at 8
        dat=dataList[i]
        time=parse.(Float64,dat[8:27]) #this is always where time values are
        timeList[i-7]=time #-7 to get back to 1
    end
    return timeList
end

function rescaleTimes(times)
    dur=maximum(times)-minimum(times)
    rescaledTimes=zeros(length(times))
    for i=1:length(rescaledTimes)
        rescaledTimes[i]=times[i]-minimum(times)
    end
    return rescaledTimes
end

function getInd(rescaledTimes,T; binNum=11) #T=period of source
    indList=zeros(length(rescaledTimes))
    for i=1:length(indList)
        mod=rescaledTimes[i]%T #value between 0 and T
        ind=floor(binNum*mod/T)+1 #ie if T is 110 and we want 11 bins and mod value was 9 the index is floor(11*9/110)+1=1
        indList[i]=ind
    end
    return indList
end

function makeLC(indList; binNum=11,srcinfo="",shifted=false,norm=false) #pick odd number so there's a middle
    bins=zeros(binNum)
    for i=1:length(indList)
        ind=Int(indList[i])
        bins[ind]+=1
    end
    GR.inline("png")
    if shifted==false
        if norm==false
            p=plot(bins,label="",title="LC for $srcinfo \nCOUNTS: $(sum(bins))",yerr=sqrt.(bins),dpi=300)
        else
            p=plot(bins./sum(bins),label="",title="Normalized LC for $srcinfo \nCOUNTS: $(sum(bins))",yerr=sqrt.(bins)./sum(bins),dpi=300)
        end
    else
        #this should be generalized with KS test at some point to accomadate shifting other than just to peak
        half=Int((length(bins)+1)/2) #assuming odd bin num
        binShift=zeros(length(bins))
        maxInd=findmax(bins)[2]
        delta=abs(maxInd-half)
        if maxInd>half #need to move things "backwards"
            for i=1:length(binShift)
                if i<=(length(bins)-delta)
                    binShift[i]=bins[i+delta]
                else
                    binShift[i]=bins[length(bins)+1-i] #+1 because always need >=1 for ind
                end
            end
            if norm==false
                p=plot(binShift,label="",title="Shifted LC for $srcinfo \nCOUNTS: $(sum(binShift))",yerr=sqrt.(binShift),dpi=300)
            else
                p=plot(binShift./sum(binShift),label="",title="Normalized shifted LC for $srcinfo \nCOUNTS: $(sum(binShift))",yerr=sqrt.(binShift)./sum(binShift),dpi=300)
            end
        elseif maxInd<half #need to move things "forwards"
            for i=1:length(binShift)
                if i>=(1+delta)
                    binShift[i]=bins[i-delta]
                else
                    binShift[i]=bins[length(bins)-delta+i]
                end
            end
            if norm==false
                p=plot(binShift,label="",title="Shifted LC for $srcinfo \nCOUNTS: $(sum(binShift))",yerr=sqrt.(binShift),dpi=300)
            else
                p=plot(binShift./sum(binShift),label="",title="Normalized shifted LC for $srcinfo \nCOUNTS: $(sum(binShift))",yerr=sqrt.(binShift)./sum(binShift),dpi=300)
            end
        else
            if norm==false
                p=plot(bins,label="",title="Shifted LC for $srcinfo \nCOUNTS: $(sum(bins))",yerr=sqrt.(bins),dpi=300)
            else
                p=plot(bins./sum(bins),label="",title="Normalized shifted LC for $srcinfo \nCOUNTS: $(sum(bins))",yerr=sqrt.(bins)./sum(bins),dpi=300)
            end
        end
    return p
    end
end

function genAllLC(; fileListPath="obsid_src_test.txt",shiftedVar=false,normalized=false)
    #println("got here")
    #srcs=CSV.read(fileListPath,header=false) #read in from txt file generated elsewhere
    srcs=DataFrame!(CSV.File(fileListPath; header=false))
    obsid=Int.(srcs.Column1[:])
    srcNum=Int.(srcs.Column2[:])
    srcT=srcs.Column3[:]
    for i=1:length(obsid)
        GR.inline("png")
        print(i/length(obsid)*100," % complete\r")
        currentObsid=obsid[i]
        currentSrcNum=srcNum[i]
        currentT=srcT[i]
        srcinfoString="OBSID: $(currentObsid) SRC: $(currentSrcNum)"
        filePath="$(currentObsid)_src_$(currentSrcNum).times" #assume run from directory with files
        data=readlines(filePath)
        times=getFloatTimes(data)
        rescaledTimes=rescaleTimes(times)
        ind=getInd(rescaledTimes,currentT)
        lc=makeLC(ind;shifted=shiftedVar,srcinfo=srcinfoString,norm=normalized)
        if shiftedVar==false
            if normalized==false
                png(lc,"plots/$(currentObsid)_src_$(currentSrcNum).png")
            else
                png(lc,"plots/$(currentObsid)_src_$(currentSrcNum)_NORMALIZED.png")
            end
        else
            if normalized==false
                png(lc,"plots/$(currentObsid)_src_$(currentSrcNum)_SHIFTED.png")
            else
                png(lc,"plots/$(currentObsid)_src_$(currentSrcNum)_SHIFTED_NORMALIZED")
            end
        end
        closeall()
    end
end

function makeLCPair(i,indList1,indList2,f1,f2; binNum=11,srcinfo1="",srcinfo2="",expected=0) #pick odd number so there's a middle
    bins1,bins2=zeros(binNum),zeros(binNum)
    for i=1:length(indList1)
        ind1=Int(indList1[i])
        bins1[ind1]+=1
    end
    for i=1:length(indList2)
        ind2=Int(indList2[i])
        bins2[ind2]+=1
    end

    half=Int((length(bins1)+1)/2) #assuming odd bin num
    maxInd1,maxInd2=findmax(bins1)[2],findmax(bins2)[2]
    delta1,delta2=abs(maxInd1-half),abs(maxInd2-half)

    function shiftBinsMid(maxInd,half,delta,bins) #shifts bins so peak is midpoint
        binShift=zeros(length(bins))
        if maxInd>half #need to move things "backwards"
            for i=1:length(binShift)
                if i<=(length(bins)-delta)
                    binShift[i]=bins[i+delta]
                else
                    binShift[i]=bins[length(bins)+1-i] #+1 because always need >=1 for ind
                end
            end
        elseif maxInd<half #need to move things "forwards"
            for i=1:length(binShift)
                if i>=(1+delta)
                    binShift[i]=bins[i-delta]
                else
                    binShift[i]=bins[length(bins)-delta+i]
                end
            end
        else
            binShift=bins
        end
        return binShift
    end

    function cumulativeBins(bins) #creates cumulative sum binning out of existing bins, not needed when using CCC
        cumulative=zeros(length(bins))
        for i=1:length(bins)
            cumulative[i]=sum(bins[1:i])
        end
        return cumulative
    end

    function findBestFit(a,b) #similar to fx in old genLC but now for CCC, in future could combine and add another param to control test
        #question -- do a and b need to scaled in some way to go better with the sine fx?
        #in paper a and b are just "angular variables"
        rList=zeros(length(a))
        counter=1
        arrangeList=[]
        reArrangement=b
        for i=1:length(a)
            test=genCCC(a,reArrangement)
            rList[i]=test
            tempArr=zeros(length(reArrangement))
            last=reArrangement[end]
            tempArr[1]=last #shift last bin in first bin
            tempArr[2:end]=reArrangement[1:(end-1)]
            reArrangement=tempArr
            push!(arrangeList,reArrangement)
        end
        maxrVal,ind=findmax(rList)
        if maximum(rList)==minimum(rList)
            println("no variance when shifting bins at i = $i")
            println("info1: $srcinfo1 info2: $srcinfo2")
        end
        return maxrVal,arrangeList[ind]
    end

    shift1=shiftBinsMid(maxInd1,half,delta1,bins1) #shit bins1 so peak is in middle
    fitStat,shift2=findBestFit(shift1,bins2) #get best match to bins1 shift
    yerr1,yerr2=sqrt.(shift1),sqrt.(shift2)
    norm1,norm2=shift1./sum(shift1),shift2./sum(shift2)
    #p=plot(norm1,label="$srcinfo1 COUNTS: $(sum(shift1))",ribbon=yerr1./sum(shift1),fillalpha=0.5)#yerr=yerr1./sum(shift1)) #VERBOSE WAY
    #p=plot!(norm2,label="$srcinfo2 COUNTS: $(sum(shift2))",ribbon=yerr2./sum(shift2),fillalpha=0.5)#yerr=yerr2./sum(shift2))
    p=plot(norm1,label="$srcinfo1:$(sum(shift1))",ribbon=yerr1./sum(shift1),fillalpha=0.5)#yerr=yerr1./sum(shift1)) #NON VERBOSE
    p=plot!(norm2,label="$srcinfo2:$(sum(shift2))",ribbon=yerr2./sum(shift2),fillalpha=0.5)#yerr=yerr2./sum(shift2))

    function genAnnotations(binShift,left;kwargs...)
        annotations=[]
        if left==true #for some reason could not get this to work in kwargs...
            for i=1:length(binShift)
                push!(annotations,(i,binShift[i]/sum(binShift),Plots.text("$(binShift[i])",8,:left;kwargs...)))
            end
        else
            for i=1:length(binShift)
                push!(annotations,(i,binShift[i]/sum(binShift),Plots.text("$(binShift[i])",8,:right;kwargs...)))
            end
        end
        return annotations
    end

    a1,a2=genAnnotations(shift1,true;color=:blue),genAnnotations(shift2,false;color=:red)
    for i=1:length(a1)
        p=annotate!(a1[i][1],a1[i][2],a1[i][3]) #add annotations to data points
        p=annotate!(a2[i][1],a2[i][2],a2[i][3])
    end

    yMax1,yMax2=maximum(norm1),maximum(norm2)
    if yMax1>yMax2
        pMax=yMax1
    else
        pMax=yMax2
    end
    GR.inline("png")
    p=plot!(title="Shifted and normalized LC pair \nr = $(@sprintf("%.3f",fitStat)) expected: $(@sprintf("%.1f",expected))\nf1 = $(@sprintf("%.6f",f1))    f2 = $(@sprintf("%.6f",f2)) (Hz)",legend=:best,dpi=300)#ylims=(0,pMax+0.08),dpi=300) #manually set ylim so legend clears data
    return p
end


function genPairLC(; fileListPath="obsid_src_test.txt") #need to both shift and normalize for this to be relevant
    #srcs=CSV.read(fileListPath,header=false) #read in from txt file generated elsewhere
    srcs=DataFrame!(CSV.File(fileListPath; header=false))
    obsid=Int.(srcs.Column1[:])
    srcNum=Int.(srcs.Column2[:])
    srcT=srcs.Column3[:]
    expected=srcs.Column4[:]
    power=srcs.Column5[:]
    alternateFileCounter=1 #some files are same obsid/same source so need to append something to file name
    for i=1:2:length(obsid)
        GR.inline("png") #this shouldn't be required, but without this the plots won't save so I added it everywhere it might be needed...new issue when upgrading to 1.4.2
        print("$(@sprintf("%.1f",i/length(obsid)*100)) % complete\r")
        obsid1=obsid[i]
        srcNum1=srcNum[i]
        obsid2=obsid[i+1]
        srcNum2=srcNum[i+1]
        T1,T2=srcT[i],srcT[i+1]
        f1,f2=1/T1,1/T2
        p1,p2=power[i],power[i+1]
        expVal=expected[i]
        #srcinfoString1="OBSID: $(obsid1) SRC: $(srcNum1)" #verbose way
        #srcinfoString2="OBSID: $(obsid2) SRC: $(srcNum2)"
        srcinfoString1="$(obsid1):$(srcNum1):$(@sprintf("%.1f",p1))" #non verbose
        srcinfoString2="$(obsid2):$(srcNum2):$(@sprintf("%.1f",p2))"
        file1="$(obsid1)_src_$(srcNum1).times" #assume run from directory with these files
        file2="$(obsid2)_src_$(srcNum2).times"
        data1,data2=readlines(file1),readlines(file2)
        times1,times2=getFloatTimes(data1),getFloatTimes(data2)
        rescaledTimes1,rescaledTimes2=rescaleTimes(times1),rescaleTimes(times2)
        ind1,ind2=getInd(rescaledTimes1,T1),getInd(rescaledTimes2,T2)
        lc=makeLCPair(i,ind1,ind2,f1,f2;srcinfo1=srcinfoString1,srcinfo2=srcinfoString2,expected=expVal)
        fileName="plots/pair_$(obsid1)_src_$(srcNum1)_w_$(obsid2)_src_$(srcNum2)"
        if isfile("$(fileName).png")==false
            #println("I got here")
            png(lc,"$(fileName).png")
            alternateFileCounter=1 #reset file counter for use if there is a duplicate
        else
            #println("I got here")
            png(lc,"$(fileName)_$(alternateFileCounter).png")
            alternateFileCounter+=1 #so next alternate can have a valid name
        end
        #png(lc,"plots/pair_$(obsid1)_src_$(srcNum1)_w_$(obsid2)_src_$(srcNum2).png") old way of doing things before if statement above
        closeall()
    end
end

#genAllLC()
#genAllLC(normalized=true)
#genAllLC(shiftedVar=true)
#genAllLC(shiftedVar=true,normalized=true)
genPairLC()
