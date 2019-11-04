#deltaT=20.  #years
dflen = size(pows9)[1]
ndiff=0
nmatch = 0
fi=[]
fj=[]
pi=[]
pj=[]
obsi=[]
obsj=[]
srci=[]
srcj=[]
rai=[]
raj=[]
dci=[]
dcj=[]

for i = 1:dflen
    for j = i:dflen
        global ndiff
        global nmatch
        if pows9[:obsid][i] != pows9[:obsid][j]
            rdist = DJM.gcdist(pows9[:ra][j],pows9[:dec][j],pows9[:ra][i],pows9[:dec][i])
            fave = (pows9[:freq][i] + pows9[:freq][j])/2.0
            Pdot = pdot(fave)
            deltaT = abs(pows9[:startdate][i] - pows9[:startdate][j])/365.25
            deltaf = fave*fave*Pdot*deltaT  #likely limit of frequency change

            df1=1/pows9[:exptime][i]
            df2=1/pows9[:exptime][j]
            df = maximum([df1,df2])
            #frange = df1+df2
            numdf = ceil(deltaf/df) #how many of the largest frequency spaces could it have drifted?

            fdiff = abs(pows9[:freq][i] - pows9[:freq][j])
            frange = df*numdf

            if (rdist*60.0 > 3*0.01666)
            #if (rdist*60.0 < 0.01666)
                ndiff = ndiff + 1
            end

            if ((rdist*60.0 > 3*0.01666) && (fdiff < frange))
            #if ((rdist*60.0 < 0.01666) && (fdiff < frange))
                println(pows9[:obsid][i]," ",pows9[:srcnum][i]," ",pows9[:freq][i]," ",pows9[:power][i]," ",pows9[:obsid][j]," ",pows9[:srcnum][j]," ",pows9[:freq][j]," ",pows9[:power][j])
                nmatch = nmatch + 1
                push!(fi,pows9[:freq][i])
                push!(fj,pows9[:freq][j])
                push!(pi,pows9[:power][i])
                push!(pj,pows9[:power][j])
                push!(obsi,pows9[:obsid][i])
                push!(obsj,pows9[:obsid][j])
                push!(srci,pows9[:srcnum][i])
                push!(srcj,pows9[:srcnum][j])
                push!(rai,pows9[:ra][i])
                push!(raj,pows9[:ra][j])
                push!(dci,pows9[:dec][i])
                push!(dcj,pows9[:dec][j])
            end

        end
    end
end
notsame=[obsi srci rai dci fi pi obsj srcj raj dcj fj pj]
#return(notsame)
