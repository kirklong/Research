#!/usr/bin/env python
import csv,subprocess,os,glob,shutil,time
import astropy.io.fits as fits
import pandas as pd
import numpy as np

# Large parts of this code made possible thanks to the examples/work done by
# Evandro Martinez Ribeiro: https://evandromr.github.io/

# This function taken from https://github.com/evandromr/python_xmmscripts
def rpcdata(odffolder, process='yes'):
    """creates a new Summary File for the observation and Reprocess XMM data"""

    if not os.path.isdir('rpcdata'):
        os.mkdir('rpcdata/')

    os.chdir('rpcdata/')
    rpcdata_dir = os.getcwd()
    # Point to Raw Observation Data directory
    os.environ['SAS_ODF'] = odffolder
    # Build calibration files index
    if process == 'yes':
        subprocess.call(['cifbuild'])
        # Point to the calibration index
        os.environ['SAS_CCF'] = os.path.abspath('ccf.cif')
        # Create observation summary
        subprocess.call(['odfingest'])
        # Point to the Summary file
        os.environ['SAS_ODF'] = os.path.abspath(glob.glob('*SUM.SAS')[0])
        # Reprocess data
        subprocess.call(['epproc'])  # PN camera
        subprocess.call(['emproc'])  # MOS cameras
    elif process == 'no':
        os.environ['SAS_CCF'] = os.path.abspath('ccf.cif')
        os.environ['SAS_ODF'] = os.path.abspath(glob.glob('*SUM.SAS')[0])
    else:
        print ("Something wrong: choose process = 'yes' or 'no'")
        raw_input("Please enter <Ctrl-C> to terminate, and check for errors")

    # (one can use "(ep/em)chain" instead, see the SAS documentation)
    os.chdir('../')

    return rpcdata_dir

#this function mostly taken from same source as above, modified slightly to remove extraneous things
#we don't want + apply barycentric correction
def clearevents(camera, rpcdata_dir='../rpcdata'):
    """Clear event file from times with high background flares"""

    curdir = os.getcwd()
    cam = camera.lower()

    if cam == 'pn':
        exp1 = "expression=#XMMEA_EP && (PI>10000&&PI<12000) && (PATTERN==0)"
        expgti = "expression=RATE<=0.4"
        exp2 = "expression=#XMMEA_EP && gti(pn_gti.ds, TIME) && (PI > 150)"
        exp3 = "expression=#XMMEA_EP && (PI>300 && PI<12000) && (PATTERN<=4)\
 && FLAG==0"
    elif cam == 'mos1':
        exp1 = "expression=#XMMEA_EM && (PI>10000) && (PATTERN==0)"
        expgti = "expression=RATE<=0.35"
        exp2 = "expression=#XMMEA_EM && gti(mos1_gti.ds, TIME) && (PI > 150)"
        exp3 = "expression=#XMMEA_EM && (PI>150 && PI<10000) && (PATTERN<=12)\
 && FLAG==0"
    elif cam == 'mos2':
        exp1 = "expression=#XMMEA_EM && (PI>10000) && (PATTERN==0)"
        expgti = "expression=RATE<=0.35"
        exp2 = "expression=#XMMEA_EM && gti(mos2_gti.ds, TIME) && (PI > 150)"
        exp3 = "expression=#XMMEA_EM && (PI>150 && PI<10000) && (PATTERN<=12)\
 && FLAG==0"
    else:
        print ("Something is wrong, the camera doesn't exist")
        raw_input("Please press 'Ctrl-C' to terminate and check errors")

    if not os.path.isdir(cam):
        os.mkdir(cam)

    os.chdir(cam)

    shutil.copyfile(os.path.abspath(
        glob.glob(rpcdata_dir+'/'+'*{0}*Evts.ds'.format(camera.upper()))[0]),
        'pnevents.ds')
    evtfile = 'pnevents.ds'

    # Creates a GTI (good time interval)
    subprocess.call(
        ['tabgtigen', 'table={0}_rate.ds'.format(cam),
         'gtiset={0}_gti.ds'.format(cam), expgti])

    # Creates a clean Events File with the events on the GTI
    subprocess.call(
        ['evselect', 'table={0}:EVENTS'.format(evtfile), 'withfilteredset=yes',
         'keepfilteroutput=yes', 'filteredset={0}_clean.ds'.format(cam), exp2])

    shutil.copy(evtfile, 'evts_barycen.ds') #create copy to apply barycentric correction to
    subprocess.call(['barycen', 'table=evts_barycen.ds:EVENTS'.format(cam)]) #call barycentric correction on clean event file

    os.chdir(curdir)

    return True

#I wrote this one!
def obsidIterator(obsidFolder, listDestFolder):
    """iterate through all obsids in a folder containing raw PPS products"""
    """return all obsids of interest and all source region file names"""
    curdir = os.getcwd()

    if os.path.isdir(obsidFolder):
        os.chdir(obsidFolder)
        subprocess.call(['ls > {0}/obsids.txt'.format(listDestFolder)],shell=True)
        with open('{0}/obsids.txt'.format(listDestFolder)) as f:
            reader = csv.reader(f)
            for row in reader:
                currObsid = row[0]
                ppsPath = obsidFolder + '/' + currObsid + '/PPS'
                os.chdir(ppsPath)
                subprocess.call(['ls *SRCREG* > {0}/{1}_srcreglist.txt'.format(listDestFolder,currObsid)],shell=True)
    else:
        print("OBSID folder does not exist as entered")
    os.chdir(curdir)

#I wrote this one!
def getSrcReg(obsid,obsidFolder):
    with open('{0}_srcreglist.txt'.format(obsid)) as f1:
        reader1 = csv.reader(f1)
        for row in reader1:
            cam = row[0][11:13] #either M1, M2, PN (assuming file name size is static)
            if cam == 'M1' or cam == 'M2' or cam == 'PN':
                filename=obsidFolder + '/' + obsid + '/PPS/' + row[0] #raw file name
                srcString = row[0][17:-4] #ie SRCREG0001
                lines = open(filename).read().splitlines()
                if not os.path.isdir('regions/{0}'.format(obsid)):
                    os.mkdir('regions/{0}'.format(obsid))
                with open('regions/{0}/{1}_{2}_{3}_PHYS.reg'.format(obsid,obsid,cam,srcString),'w') as outFile:
                    outFile.write("# Region file format: DS9 version 4.1\n")
                    outFile.write('global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
                    outFile.write("physical\n")
                with open('regions/{0}/{1}_{2}_{3}_FK5.reg'.format(obsid,obsid,cam,srcString),'w') as outFile:
                    outFile.write("# Region file format: DS9 version 4.1\n")
                    outFile.write('global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
                    outFile.write("FK5\n")
                for line in lines:
                    #print(line[:8])
                    if line[:8] == 'detector':
                        for i in range(8,len(line)):
                            if line[i] == 'c' and line[i+1]=='i':
                                with open('regions/{0}/{1}_{2}_{3}_PHYS.reg'.format(obsid,obsid,cam,srcString),'a') as outFile:
                                    outFile.write("{0}\n".format(line[i:]))
                                break
                    elif line[:3] == 'fk5':
                        for i in range(3,len(line)):
                            if line[i] == 'c' and line[i+1]=='i':
                                with open('regions/{0}/{1}_{2}_{3}_FK5.reg'.format(obsid,obsid,cam,srcString),'a') as outFile:
                                    outFile.write("{0}\n".format(line[i:]))
                                break
            else:
                print("something is wrong -- file size not static?")
    return srcString #not relevant for this function, but this is a hack way to use it later

#I wrote this one!
def getCircles(regListFile):
    src=open(regListFile)
    srcregions=src.readlines()
    src.close()
    srcregions=[src.strip() for src in srcregions[3:]] #start at 3 to discard header
    circles = []
    for i in range(len(srcregions)):
        for j in range(len(srcregions[i])):
            if srcregions[i][j] == ")": #stop at second paranthetical
                circles.append(srcregions[i][:j+1]) #append just circle portion of text file
                break
    return circles

#taken from source above but modified somewhat substantially
def lcextraction(abin, arange, emin, emax, table, srcregion, camera,
                 tstart, tstop, OBSID, srcString, iter):
    ''' Extract a barycentric corrected lightcurve
        in the energy range['emin':'emax'], with time bins of 'bin'
    '''
    #srclc = "0810880201_M1_lc_SRCREG0001_1_0310keV_bin1.ds"
    srclc = "{0}_{1}_lc_{2}_{3}_{4}keV_bin{5}.ds".format(OBSID, cam, srcString, iter, arange, abin)
    srcimg = "{0}_{1}_img_{2}_{3}_{4}keV_bin{5}.ds".format(OBSID, cam, srcString, iter, arange, abin)
    psimg = "{0}_{1}_lcPlot_{2}_{3}_{4}keV_bin{5}.ps".format(OBSID, cam, srcString, iter, arange, abin)

    if cam == 'PN':
        srcexp = 'expression=#XMMEA_EP && (PI IN [{0}:{1}]) && PATTERN <=4\
 && FLAG==0 && ((X,Y) IN {2})'.format(emin, emax, srcregion)

    elif cam == 'M1':
        srcexp = 'expression=#XMMEA_EM && (PI IN [{0}:{1}]) && PATTERN <=12\
 && FLAG==0 && ((X,Y) IN {2})'.format(emin, emax, srcregion)

    elif cam == 'M2':
        srcexp = 'expression=#XMMEA_EM && (PI IN [{0}:{1}]) && PATTERN <=12\
 && FLAG==0 && ((X,Y) IN {2})'.format(emin, emax, srcregion)

    else:
        print ("Something is wrong, the camera doesn't exist")
        raw_input("Please press 'Ctrl-C' to terminate and check for errors")
    currdir = os.getcwd()
    if not os.path.isdir('lc/{0}'.format(OBSID)):
        os.mkdir('lc/{0}'.format(OBSID))
    os.chdir('lc/{0}'.format(OBSID)) #change into lc directory so we don't spam pipeline directory
    subprocess.call(
        ['evselect', 'table={0}'.format(table),'energycolumn=PI','rateset={0}'.format(srclc),
        'timebinsize={0}'.format(abin),'maketimecolumn=yes','makeratecolumn=yes',
        'imageset={0}'.format(srcimg),'xcolumn=X', 'ycolumn=Y','ximagebinsize=80','yimagebinsize=80',
        'timemin={0}'.format(tstart), 'timemax={0}'.format(tstop),srcexp] #testing simple version
    )

    subprocess.call(
        ['dsplot', 'table={0}'.format(srclc), 'withx=yes', 'withy=yes',
         'x=TIME', 'y=RATE',
         'plotter=xmgrace -hardcopy -printfile {0}'.format(psimg)])
    os.chdir(currdir)
    return True

#taken from source above, relatively unchanged
def findtimes(rpcdata_dir, camera):
    """Find the initial and final times of observation for the given camera"""
    currentDir = os.getcwd()
    os.chdir(rpcdata_dir)
    #print(os.getcwd())
    if camera == 'M1':
        camera = 'MOS1'
    elif camera == 'M2':
        camera = 'MOS2'
    elif camera == 'PN':
        camera = 'PN'
    else:
        print('camera not entered correctly')
    evtfile = glob.glob('*{0}*ImagingEvts.ds'.format(camera.upper()))
    #print(evtfile
    tstart = fits.getval(evtfile[0], 'TSTART', ext=1)
    tstop = fits.getval(evtfile[0], 'TSTOP', ext=1)
    os.chdir(currentDir)
    return tstart, tstop

#I wrote this one!
def lc2CSV(infile,obsid):
    currdir = os.getcwd()
    os.chdir('lc/{0}'.format(obsid)) #change into lc directory where files are located
    outfile = infile[:-3] + '.csv' #change the extension from .ds to .csv
    subprocess.call(['fdump {0}[1] {1} rows=- columns=- prhead=no, fldsep=","'.format(infile, outfile)],shell=True)
    os.chdir(currdir) #back out to pipeline
    outFilePath = "{0}/lc/{1}/{2}".format(currdir,obsid,outfile)
    return outFilePath #don't print the header, print all the rows and columns from the rate index, return file name for later use

#converting Julia code I wrote to Python...stupid.
#I wrote this one!
def csv2DF(csvFile):
    rate =[]
    time =[]
    error =[]
    with open(csvFile) as f1:
        reader1 = csv.reader(f1)
        rowCount = 1
        for row in reader1:
            if rowCount > 4:
                rate.append(float(row[1]))
                error.append(float(row[2]))
                time.append(float(row[3]))
            rowCount+=1
    df = pd.DataFrame({'rate' : np.array(rate), 'time' : np.array(time), 'error' : np.array(error)})
    return df

#I wrote this one!
def genFFTDF(df):
    fftOut=np.fft.fft(df.rate[:]) #rate is really counts (counts/s recorded every sec)
    fftOut=fftOut[1:] #get rid of Nyqhuist
    fftOut=fftOut[:len(fftOut)//2] #take half (rounding midpoint down if needed)
    fftOut=np.abs(fftOut)**2 #square magnitude
    avgPow=np.mean(fftOut)
    fftOut=fftOut/avgPow #normalize
    dt=df.time[len(df.time[:])-1]-df.time[0] #for some reason -1 indexing wasn't working here whatever
    freq=[i/dt for i in range(len(fftOut))] #freq(i) = bin#/totalTime
    fftDF=pd.DataFrame({'freq' : freq, 'power' : fftOut})
    return fftDF

#I wrote this one!
def saveFFT(obsid,cam,srcReg,srcNum,fftDF):
    currdir = os.getcwd()
    os.chdir('FFT')
    if not os.path.isdir('{0}'.format(obsid)):
        os.mkdir('{0}'.format(obsid))
    os.chdir('{0}'.format(obsid))
    exportFile="{0}_{1}_{2}_{3}_FFT.csv".format(obsid,cam,srcReg,srcNum)
    fftDF.to_csv(exportFile,index=False)
    os.chdir(currdir)

#I wrote this one!
def createAllInfo(powCutoff,obsid,cam,srcReg,srcNum,startdate,ra,dec,fitsDF,fftInfo):
    outFile = "xmm.fftinfo.all{0}".format(powCutoff)
    exptime = fitsDF.time[len(fitsDF.time[:])-1]-fitsDF.time[0] #-1 returns key error for some reason? whatever
    obsidCam = "{0}_{1}".format(obsid,cam)
    srcRegNum = "{0}_{1}".format(srcReg,srcNum)
    with open(outFile,'a') as f:
        for i in range(len(fftInfo.power[:])):
            if fftInfo.power[i] >= powCutoff:
                try: #for some reason I can't write to the file?
                    freq = fftInfo.freq[i]
                    power = fftInfo.power[i]
                    f.write("{0},{1},{2},{3},{4},{5},{6},{7}\n".format(obsidCam,srcRegNum,startdate,exptime,ra,dec,freq,power))
                except:
                    print("error in writing to file -- logging")
                    errFile="/home/kirk/Documents/research/XMM_Newton/python_xmmscripts/pipeline/err.txt"
                    with open(errFile,'a') as f1:
                        f1.write("{0},{1},index={2}\n{3}\n".format(obsidCam,srcRegNum,i,fftInfo))
                        return 1 #if this stupid ass bug pops up add to the errCount
    return 0 #if no problems, don't add to errCount

# OVERALL SCRIPTING PROCESS
# BEFORE RUNNING: HEAsoft and SAS need to be initialized, and the SAS variable SAS_CCFPATH needs to be defined (script will do the other two)

#1 get obsids -- started 12:13
print("getting obsids -- STEP 1")
obsidFolder = "/home/kirk/Documents/research/XMM_Newton/obsids" #my machine
#obsidFolder = "/Users/djm/Public/xmm/smc/obsids" #Daryl's computer
listDestFolder = '/home/kirk/Documents/research/XMM_Newton/python_xmmscripts/pipeline'
#listDestFolder = "/Users/djm/Public/xmm/smc/pipeline" #Daryl's computer
obsidIterator(obsidFolder, listDestFolder)
#2 get srcreg files LOOP (i)
f = open('obsids.txt')
obsids = f.readlines()
f.close()
errCount=0
clean=False
for obsid in obsids:
    odfFolder = "{0}/{1}/ODF".format(obsidFolder,obsid.strip())
    #time to run set up nonsense is roughly 10 min (informal test)
    #0 SAS set up from other guy's XMM Python stuff (not really modified)
    #NOTE: this takes a substantial amount of time to run (~5-10 min per obsid) and some/all of this may be unneccessary?
    #print ("Reprocessing data for OBSID: {0}".format(obsid))
    #rpcdata_dir = rpcdata(odfFolder) #create the rpcdata directory if it doesn't exist, run processing tasks like odfingest and cifbuild
    #print ("Finished data processing\nCleaning events...")
    #clearevents(camera='pn', rpcdata_dir=rpcdata_dir)
    #clearevents(camera='mos1', rpcdata_dir=rpcdata_dir)
    #clearevents(camera='mos2', rpcdata_dir=rpcdata_dir)
    print ("DONE")
    obsid = obsid[:-1].strip() #get rid of \n character
    print("getting srcreg files -- STEP 2")
    srcString = getSrcReg(obsid,obsidFolder)
    currdir = os.getcwd()
    if not os.path.isdir('regions/{0}'.format(obsid)):
        os.mkdir('regions/{0}'.format(obsid))
    os.chdir('regions/{0}'.format(obsid))
    subprocess.call(['ls > tempList.txt'],shell=True) #hack way to loop over file names
    f = open('tempList.txt')
    files = f.readlines()
    f.close()
    os.chdir(currdir)
    #3 get detector circles for SAS commands LOOP (j)
    for name in files:
        name = name[:-1] #remove \n
        print("getting detector circles -- STEP 3")
        if name[-8:-4] == 'PHYS':
            #set up for lc extraction
            emin=200 #eV
            emax=12000 #eV
            cam = name.split("_")[1]
            if cam == 'M1':
                camFolder = 'mos1'
            elif cam == 'M2':
                camFolder = 'mos2'
            elif cam == 'PN':
                camFolder = 'pn'
            else:
                print("camera option error")
            table='{0}/{1}/evts_barycen.ds'.format(listDestFolder,camFolder)
            abin=1 #1second lc
            arange = '{0}{1}'.format(emin/100,emax/100)
            tstart, tstop = findtimes('{0}/rpcdata/'.format(listDestFolder), cam)
            circles = getCircles("regions/{0}".format(name))
            fk5File = name[:-8]+'FK5.reg'
            f = open("regions/{0}/{1}".format(obsid,fk5File))
            raDec = f.readlines()
            f.close()
            print("getting RA and DEC info -- STEP 4")
            raAvg = 0
            decAvg = 0
            for i in range(3,len(raDec)): #3 offsets because of header
                raDecSplit = raDec[i].split(",")
                ra = float(raDecSplit[0].split()[1]) #split twice to get rid of preceeding "circle( "
                raAvg += ra/(len(raDec)-3) #average ra
                dec = float(raDecSplit[1])
                decAvg += dec/(len(raDec)-3) #average dec

            #5 extract LCs
            for i in range(len(circles)):
                print("extracting LCs -- STEP 5")
                srcregion = circles[i]
                lcextraction(abin,arange,emin,emax,table,srcregion,cam,tstart,tstop,obsid,srcString,i)
                #5 send lc to CSV LOOP (i)
                print("making LC CSV -- STEP 6")
                infile="{0}_{1}_lc_{2}_{3}_{4}keV_bin{5}.ds".format(obsid, cam, srcString, i, arange, abin)
                csvFile = lc2CSV(infile,obsid)
                #6 get RA and DEC info
                # print("getting RA and DEC -- STEP 6")
                # raDec = raDec[i+3].split(",") #+3 because need to skip header, then raDec[i] should correspond to physical circle referenced
                # ra = raDec[0].split()[1] #second split needed to remove leading "circle( " business
                # dec = raDec[1]
                ra = raAvg #taking average because for some reason there are less RA/DEC circles than PHYS circles but this is probs close enough
                dec = decAvg
                #7 make and save FFT
                print("making FFTs -- STEP 7")
                csvDF = csv2DF(csvFile)
                fftDF = genFFTDF(csvDF)
                saveFFT(obsid,cam,srcString,i,fftDF)
                #8 write to all info file
                startdate = tstart
                powCutoff = 8
                print("writing to all info file -- STEP 8")
                errCount+=createAllInfo(powCutoff,obsid,cam,srcString,i,startdate,ra,dec,csvDF,fftDF)
                #7 call Julia program to generate fft and allInfo list -- sadly not using becauase can't get julia initialization time down
                #subprocess.call(['fftFileGen.jl {0} {1} {2} {3} {4} {5} {6} {7} {8} {9}'.format(obsid,cam,srcString,i,ra,dec,startdate,abin,arange,powCutoff)],shell=True)
    #clear out old files for next obsid
    if clean == True:
        print("cleaning up files")
        currdir = os.getcwd()
        print(currdir)
        try:
            os.chdir('rpcdata')
            subprocess.call(['rm *'],shell=True)
            os.chdir(currdir)
        except:
            print("could not access rpcdata directory")
        try:
            os.chdir('mos1')
            subprocess.call(['rm *'],shell=True)
            os.chdir(currdir)
        except:
            print("could not access mos1 directory")
        try:
            os.chdir('mos2')
            subprocess.call(['rm *'],shell=True)
            os.chdir(currdir)
        except:
            print("could not access mos1 directory")
        try:
            os.chdir('pn')
            subprocess.call(['rm *'],shell=True)
            os.chdir(currdir)
        except:
            print("could not access mos1 directory")

    print("finished cleaning, starting next obsid")
print("congrats -- you made it! \n total errors: {0}".format(errCount))
#finished ~
