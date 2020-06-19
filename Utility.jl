

#for now assume decimal degrees
function gcdist(lon1,lat1,lon2,lat2) #Kirk changed

    d2r = pi/180.0
    lon1 = lon1*d2r
    lon2 = lon2*d2r
    lat1 = lat1*d2r
    lat2 = lat2*d2r

    longdif =  abs(lon2-lon1)
    acosArg = sin(lat1)*sin(lat1) + cos(lat1)*cos(lat1)*cos(longdif)
    if acosArg > 1
        acosArg = 1
    end
    rdist = acos(acosArg)
    return rdist/d2r
    # try
    #     rdist = acos(sin(lat1)*sin(lat1) + cos(lat1)*cos(lat1)*cos(longdif))
    #     return rdist/d2r #return in degrees
    # catch
    #     return 0
    # end
     #return 0 if acos arg is 1.0000000000000002 (result of consistent overflow problem)

end

function gcdist(lon1::Vector,lat1::Vector,lon2,lat2)
  d2r = pi/180.0
  lon1 = lon1*d2r
  lon2 = lon2*d2r
  lat1 = lat1*d2r
  lat2 = lat2*d2r

  latdif =  lat1.-lat2
  rdist = acos.(sin.(lon1).*sin(lon2) + cos.(lon1).*cos(lon2).*cos.(latdif))

  return(rdist/d2r)  #return decimal degrees

end

function overlap(lon1::Vector,lat1::Vector,lon2,lat2)
  d2r = pi/180.0
  lon1 = lon1*d2r
  lon2 = lon2*d2r
  lat1 = lat1*d2r
  lat2 = lat2*d2r

  latdif =  lat1.-lat2
  rdist = acos.(sin.(lon1).*sin(lon2) + cos.(lon1).*cos(lon2).*cos.(latdif))

  return(rdist/d2r)  #return decimal degrees

end
