function pdot(f)
    Pdot2 = log10(0.01)  #seconds/year
    Pdot1 = log10(0.00001)
    P1 = log10(1.0)
    P2 = log10(1000.0)
    plog=log10(1/f)

    dpdot = (Pdot2 - Pdot1)/(P2-P1)
    logpdot = Pdot1 + (plog-P1)*dpdot
    return(10^logpdot)
end
