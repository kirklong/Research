#This file contains various fxs for obtaining proabilities from background bgrndMap

#function get_probability returns probability of specific power pair at specific freq bin
function get_probability(bgrndMap,pair,freqBin)
    pow1,pow2=pair[1],pair[2]
    if freqBin==0
        kInd=1 #first slice is general map #THIS IS REDUNDANT BULLSHIT
    else
        kInd=freqBin+1 #frequency specific bins start 1 slice back
    end
    omegaTotal=sum(bgrndMap[:,:,kInd]) #total # of microstates in slice
    omegaPow=bgrndMap[pow1,pow2,kInd]
    if omegaTotal==0 #sometimes there are no
        p=0
    else
        p=omegaPow/omegaTotal #probability equals # in microstate/total # of possibilities
    end
    return p
end

#function generateProbList4Pow generates a probability list for a single pow with all other pow entries at specific freq bin
function generateProbList4Pow(bgrndMap,pow,freqbin)
    probPowWOthers=zeros(size(bgrndMap)[1])
    for i in 1:length(probPowWOthers)
        probPowWOthers[i]=get_probability(bgrndMap,[pow,i],freqbin)
    end
    return probPowWOthers
end

#function generateProbALLPow does the same thing as above function but for all possible pow pairs at specific freq bin
function generateProbALLPow(bgrndMap,freqbin)
    powList=range(1,stop=101,length=101) #pow bins are 1-100 with 101 representing all pow >100
    probMatrix=zeros(length(powList),length(powList))
    for pow in powList
        probMatrix[:,Int(pow)]=generateProbList4Pow(bgrndMap,Int(pow),freqbin) #column = results of prob function
    end
    return probMatrix
end

#function generateAllProb builds on generateProbALLPow by generalizing over all frequency bins
function generateAllProb(bgrndMap)
    nk=size(bgrndMap)[3] #size of k
    nij=size(bgrndMap)[1] #size of i,j
    probMatrix=zeros(nij,nij,nk)
    for freqbin=0:(nk-1) #start at 0 because of how I defined kInd in get_probability fx
        kInd=freqbin+1
        print(format(kInd/nk*100,precision=2),"% complete\r")
        probMatrix[:,:,kInd]=generateProbALLPow(bgrndMap,freqbin)
    end
    return probMatrix
end
