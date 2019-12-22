using Random
function createBootStrapMatrix(powList,nTrials,iterNum)  #nTrials for random points, iterNum for number of times this function will run
    n=ceil(maximum(powList.power[:]))
    df=0.0001
    nk=ceil(maximum(powList.freq[:])/df)
    powMatrix=zeros(101,101,Int(nk)+1) #any power greater than 100 will go in 101,nk for freq z direction, +1 because 1 is for normal matrix
    justPow=powList.power[:]
    counter=0
    while counter<nTrials
        print("Iteration ",iterNum," ",format(counter/nTrials*100,precision=2),"% complete\r") #output % tracker
        powInd=rand(1:length(justPow),2) #two random numbers between 1 and n
        pow1=powInd[1]
        pow2=powInd[2]
        obsid1=powList.obsid[pow1] #used for checking later--we want to exclude pairs in same observations
        if ceil(powList.power[pow1])>100
            currentRow=101 #power is too big
        else
            currentRow=ceil(powList.power[pow1]) #row index for use later
        end
        obsid2=powList.obsid[pow2]
        if obsid1!=obsid2
            if ceil(powList.power[pow2])>100
                currentCol=101 #other power too big
            else
               currentCol=ceil(powList.power[pow2]) #column index for use later
            end
            rdist = DJM.gcdist(powList.ra[pow2],powList.dec[pow2],powList.ra[pow1],powList.dec[pow1])
            if (rdist*60.0 > 0.01666) #not in the same place
                counter+=1 #only count it as a trial if the obsids are different and they're in different places
                powMatrix[Int(currentRow),Int(currentCol),1]+=1 #found a pair, doesn't account for overcounting (indistinguishability)
                #stuff below for freq z part of matrix
                fAvg=(powList.freq[pow1]+powList.freq[pow2])/2
                pAvg=1/fAvg
                pDotAvg=pdot(fAvg) #pdot takes freq but returns period
                dt=abs(powList.startdate[pow1]-powList.startdate[pow2])*24*60*60 #days to seconds
                pMax=pAvg+dt*pDotAvg #upper bound for period from accretion model
                pMin=pAvg-dt*pDotAvg
                fRangeMax=1/pMin #upper freq bound corresponds to lower period bound
                fRangeMin=1/pMax
                if fAvg>=fRangeMin && fAvg<=fRangeMax
                    kInd=ceil(fAvg/df)
                    powMatrix[Int(currentRow),Int(currentCol),Int(kInd+1)]+=1 #kInd+1 because 1 is base matrix
                end
            end
        end
    end
    return powMatrix
end

function genStdMatrix(iterNum,powList,nTrials)
    allMatrices=zeros(iterNum)
    for i=1:iterNum
        allMatrices[i]=createBootStrapMatrix(powList,nTrials,i)
    end
    df=0.0001
    nk=ceil(maximum(powList.freq[:])/df)
    stdMatrix=zeros(101,101,Int(nk)+1)
    matrixCounter,freqCounter,rowCounter=1,1,1
    for f=2:(nk+1) #current frequency, +1 shift because 1st slice is not freq dependent
        print(format(f/(nK+1)*100,precision=2),"% complete\r") #output % tracker
        for i=1:101 #row counter
            for j=1:101 #col counter
                stop=false
                m=1
                data=zeros(length(allMatrices)) #placeholder to get accumulate accross matrices
                while stop==false
                    currentMatrix=allMatrices[m]
                    currentEntry=currentMatrix[i,j,f]
                    data[m]+=currentEntry
                    if m<length(allMatrices)
                        m+=1
                    else #reached the end of loop, find std and place in appropriate cell
                        currentStd=std(data)
                        stdMatrix[i,j,f]=currentStd
                        stop=true
                    end
                end
            end
        end
    end
    return stdMatrix
end
