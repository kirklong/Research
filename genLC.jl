#!/usr/bin/env julia
using Plots, CSV
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

function getInd(rescaledTimes)
    indList=zeros(length(rescaledTimes))
    for i=1:length(indList)
        indList[i]=floor(rescaledTimes[i]%10)+1
    end
    return indList
end

function makeLC(indList; binNum=11,srcinfo="",shifted=false,norm=false) #pick odd number so there's a middle
    bins=zeros(binNum)
    for i=1:length(indList)
        ind=Int(indList[i])
        bins[ind]+=1
    end
    if shifted==false
        if norm==false
            p=plot(bins,label="",title="LC for $srcinfo \nCOUNTS: $(sum(bins))",yerr=sqrt.(bins))
        else
            p=plot(bins./sum(bins),label="",title="Normalized LC for $srcinfo \nCOUNTS: $(sum(bins))",yerr=sqrt.(bins)./sum(bins))
        end
    else
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
                p=plot(binShift,label="",title="Shifted LC for $srcinfo \nCOUNTS: $(sum(binShift))",yerr=sqrt.(binShift))
            else
                p=plot(binShift./sum(binShift),label="",title="Normalized shifted LC for $srcinfo \nCOUNTS: $(sum(binShift))",yerr=sqrt.(binShift)./sum(binShift))
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
                p=plot(binShift,label="",title="Shifted LC for $srcinfo \nCOUNTS: $(sum(binShift))",yerr=sqrt.(binShift))
            else
                p=plot(binShift./sum(binShift),label="",title="Normalized shifted LC for $srcinfo \nCOUNTS: $(sum(binShift))",yerr=sqrt.(binShift)./sum(binShift))
            end
        else
            if norm==false
                p=plot(bins,label="",title="Shifted LC for $srcinfo \nCOUNTS: $(sum(bins))",yerr=sqrt.(bins))
            else
                p=plot(bins./sum(bins),label="",title="Normalized shifted LC for $srcinfo \nCOUNTS: $(sum(bins))",yerr=sqrt.(bins)./sum(bins))
            end
        end
    return p
    end
end

function genAllLC(; fileListPath="obsid_src_test.txt",shiftedVar=false,normalized=false)
    #println("got here")
    srcs=CSV.read(fileListPath,header=false) #read in from txt file generated elsewhere
    obsid=srcs.Column1[:]
    srcNum=srcs.Column2[:]
    for i=1:length(obsid)
        print(i/length(obsid)*100," % complete\r")
        currentObsid=obsid[i]
        currentSrcNum=srcNum[i]
        srcinfoString="OBSID: $(currentObsid) SRC: $(currentSrcNum)"
        filePath="$(currentObsid)_src_$(currentSrcNum).times" #assume run from directory with files
        data=readlines(filePath)
        times=getFloatTimes(data)
        rescaledTimes=rescaleTimes(times)
        ind=getInd(rescaledTimes)
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

function makeLCPair(indList1,indList2; binNum=11,srcinfo1="",srcinfo2="") #pick odd number so there's a middle
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

    function shiftBins(maxInd,half,delta,bins)
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

    shift1,shift2=shiftBins(maxInd1,half,delta1,bins1),shiftBins(maxInd2,half,delta2,bins2)
    yerr1,yerr2=sqrt.(shift1),sqrt.(shift2)
    norm1,norm2=shift1./sum(shift1),shift2./sum(shift2)
    p=plot(norm1,label="$srcinfo1 COUNTS: $(sum(shift1))",ribbon=yerr1./sum(shift1),fillalpha=0.5)#yerr=yerr1./sum(shift1))
    p=plot!(norm2,label="$srcinfo2 COUNTS: $(sum(shift2))",ribbon=yerr2./sum(shift2),fillalpha=0.5)#yerr=yerr2./sum(shift2))
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
    p=plot!(title="Shifted and normalized LC pair",legend=:best,ylims=(0,pMax+0.08)) #manually set ylim so legend clears data
    return p
end


function genPairLC(; fileListPath="obsid_src_test.txt") #need to both shift and normalize for this to be relevant
    srcs=CSV.read(fileListPath,header=false) #read in from txt file generated elsewhere
    obsid=srcs.Column1[:]
    srcNum=srcs.Column2[:]
    for i=1:2:length(obsid)
        print(i/length(obsid)*100," % complete\r")
        obsid1=obsid[i]
        srcNum1=srcNum[i]
        obsid2=obsid[i+1]
        srcNum2=srcNum[i+1]
        srcinfoString1="OBSID: $(obsid1) SRC: $(srcNum1)"
        srcinfoString2="OBSID: $(obsid2) SRC: $(srcNum2)"
        #filePath="pair$(i)_$(obsid1)_src_$(srcNum1)_w_$(obsid2)_src_$(srcNum2).times" #assume run from directory with files
        file1="$(obsid1)_src_$(srcNum1).times"
        file2="$(obsid2)_src_$(srcNum2).times"
        data1,data2=readlines(file1),readlines(file2)
        times1,times2=getFloatTimes(data1),getFloatTimes(data2)
        rescaledTimes1,rescaledTimes2=rescaleTimes(times1),rescaleTimes(times2)
        ind1,ind2=getInd(rescaledTimes1),getInd(rescaledTimes2)
        lc=makeLCPair(ind1,ind2;srcinfo1=srcinfoString1,srcinfo2=srcinfoString2)
        png(lc,"plots/pair_$(obsid1)_src_$(srcNum1)_w_$(obsid2)_src_$(srcNum2).png")
        closeall()
    end
end

#genAllLC()
#genAllLC(normalized=true)
#genAllLC(shiftedVar=true)
#genAllLC(shiftedVar=true,normalized=true)
#genPairLC()
