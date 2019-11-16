#this file contains all the background matrix generation tasks

#pdot required for some calculations later
function pdot(f) #change in period due to accretion, serious estimating involved
    Pdot2 = -2  #seconds/year
    Pdot1 = -5
    P1 = 0
    P2 = 3
    plog=log10(1/f)

    dpdot = (Pdot2 - Pdot1)/(P2-P1)
    logpdot = Pdot1 + (plog-P1)*dpdot
    return(10^logpdot)
end

#this function creates massive background power matrix, generalized now for both overlapping and non-overlapping cases
function createPowMatrix(powList,overLapBool) #overLapBool: put 1 for True (overlapping) or 0 for false (non-overlapping)
    n=ceil(maximum(powList.power[:]))
    df=0.0001
    nk=ceil(maximum(powList.freq[:])/df)
    powMatrix=zeros(101,101,Int(nk)+1) #any power greater than 100 will go in 101,nk for freq z direction, +1 because 1 is for normal matrix
    justPow=powList.power[:]
    srcInd=[]
    sameCounter=0 #used to keep track of # of pairs excluded as result of same obsid
    for i in 1:length(justPow)
        print(format(i/length(justPow)*100,precision=2),"% complete\r") #output % tracker
        obsid1=powList.obsid[i] #used for checking later--we want to exclude pairs in same observations
        if ceil(powList.power[i])>100
            currentRow=101 #power is too big
        else
            currentRow=ceil(powList.power[i]) #row index for use later
        end

        for j in (i+1):length(justPow) #i+1 accounts for indistinguishability
            obsid2=powList.obsid[j]
            if obsid1!=obsid2

                if ceil(powList.power[j])>100
                    currentCol=101 #other power too big
                else
                   currentCol=ceil(powList.power[j]) #column index for use later
                end
                if overLapBool==0 #non-overlapping
                    rdist = DJM.gcdist(powList.ra[j],powList.dec[j],powList.ra[i],powList.dec[i])
                    if (rdist*60.0 > 0.01666) #not in the same place
                        powMatrix[Int(currentRow),Int(currentCol),1]+=1 #found a pair, doesn't account for overcounting (indistinguishability)
                        #stuff below for freq z part of matrix
                        fAvg=(powList.freq[i]+powList.freq[j])/2
                        pAvg=1/fAvg
                        pDotAvg=pdot(fAvg) #pdot takes freq but returns period
                        dt=abs(powList.startdate[i]-powList.startdate[j])*24*60*60 #days to seconds
                        pMax=pAvg+dt*pDotAvg #upper bound for period from accretion model
                        pMin=pAvg-dt*pDotAvg
                        fRangeMax=1/pMin #upper freq bound corresponds to lower period bound
                        fRangeMin=1/pMax
                        if fAvg>=fRangeMin && fAvg<=fRangeMax
                            kInd=ceil(fAvg/df)
                            powMatrix[Int(currentRow),Int(currentCol),Int(kInd+1)]+=1 #kInd+1 because 1 is base matrix
                        end
                    end
                else #overlapping case
                    rdist = DJM.gcdist(powList.ra[j],powList.dec[j],powList.ra[i],powList.dec[i])
                    if (rdist*60.0 < 0.01666) #MODIFIED NOW IN THE SAME PLACE
                        powMatrix[Int(currentRow),Int(currentCol),1]+=1 #found a pair, doesn't account for overcounting (indistinguishability)
                        #stuff below for freq z part of matrix
                        fAvg=(powList.freq[i]+powList.freq[j])/2
                        pAvg=1/fAvg
                        pDotAvg=pdot(fAvg) #pdot takes freq but returns period
                        dt=abs(powList.startdate[i]-powList.startdate[j])*24*60*60 #days to seconds
                        pMax=pAvg+dt*pDotAvg #upper bound for period from accretion model
                        pMin=pAvg-dt*pDotAvg
                        fRangeMax=1/pMin #upper freq bound corresponds to lower period bound
                        fRangeMin=1/pMax
                            if fAvg>=fRangeMin && fAvg<=fRangeMax
                                push!(srcInd,[(i,j)]) #i and j are location of power pairs in original data AND frequencies close
                                kInd=ceil(fAvg/df)
                                powMatrix[Int(currentRow),Int(currentCol),Int(kInd+1)]+=1 #kInd+1 because 1 is base matrix
                            end
                    end
                end
            else
                sameCounter+=1
            end
        end
    end
    println(sameCounter," potential matches had the same obsid and were excluded.")
    if overLapBool==0 #for non-overlapping case we don't care about specific sources
        return powMatrix
    else
        return powMatrix,srcInd
    end
end
