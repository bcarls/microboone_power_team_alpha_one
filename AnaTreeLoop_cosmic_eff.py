#!/usr/local/bin/python

import ROOT
import math
import rootplot.root2matplotlib as r2m
import matplotlib.pyplot as plt
import pickle
import sys

class PyTree(ROOT.TTree):
    def Hi(self):
        return 'Hi'

def getTreeVariable(t,trackModule,varName):
    return getattr(t,varName+"_"+trackModule)


def main(argv):

    inputFile = ''
    inputFileList = ''
    nevts = 0
    # Parse arguments
    args = argv[1:]
    while len(args) > 0:
        if args[0] == '-s' or args[0] == '--source' :
            if len(args) > 1:
                inputFile = args[1]
                del args[0:2]
        elif args[0] == '-S' or args[0] == '--source-list' :
            if len(args) > 1:
                inputFileList = args[1]
                del args[0:2]
        elif args[0] == '-n' or args[0] == '--nevts' :
            if len(args) > 1:
                nevts = int(args[1])
                del args[0:2]

    if inputFile == '' and inputFileList == '':
        print 'No input file(s) specified. Use "-s" or "--source" to specify one. Additionally, a file list can be supplied with "-S" or "--source-list".'
        return 0

    # Track module
    #trackModule = "trackkalmanhitpandoraneutrino"
    #trackModule = "trackkalmanhitpandoracosmic"
    trackModule = "trackkalmanhit"
    #trackModule = "trackkalmanhitcc"
    #trackModule = "trackkalmanhitpandora"

    # The stitching track modules
    #trackModule = "stitchcc"

    plt.rc('text', usetex=True)
    #plt.rc('font', family='serif')
   


    # Read in file(s) and create the TChain c
    inputFiles = []
    c = ROOT.TChain("analysistree/anatree") 
    if inputFile != '':
        inputFiles.append(inputFile)
    elif inputFileList != '':
        inputFiles = open(inputFileList).read().splitlines()

    for f in inputFiles:
        print 'Adding input file: '+f
        c.Add(f)

    nEntries = c.GetEntries()
    if nevts == 0:
        nevts = nEntries
    print 'Found '+str(nEntries)+' events.'
    print 'Looping over '+str(nevts)+' of them.'

    hNumTracks = []
    hTrackPhi = []
    hLongestTrackPhi = []
    hTrackTheta = []
    hLongestTrackTheta = []
    hTrackThetaXZ = []
    hLongestTrackThetaXZ = []
    hTrackThetaYZ = []
    hLongestTrackThetaYZ = []
    hTrackMom = []
    hTrackLen = []
    hLongestTrackLen = []
    hTrackCosmicScore = []
    hLongestTrackCosmicScore = []
    hTrackPIDA = []
    hLongestTrackPIDA = []
    hTrackNumHits = []
    hLongestTrackNumHits = []
  
    # Unfortunately, we have to use ROOT to divide histograms, use it to handle the efficiency histograms
    # Look at the efficiency of cosmic tagging as a function of the number of track hits
    hNumHitsTracksTrueCosmics = ROOT.TH1F("hNumHitsTracksTrueCosmics ","",10,0,2000)
    hNumHitsTracksTrueCosmics.Sumw2()
    hNumHitsTracksFoundCosmics = ROOT.TH1F("hNumHitsTracksFoundCosmics ","",10,0,2000)
    hNumHitsTracksFoundCosmics.Sumw2()
    # Look at the mistag rate
    hNumHitsTracksTrueNotCosmics = ROOT.TH1F("hNumHitsTracksTrueNotCosmics ","",10,0,2000)
    hNumHitsTracksTrueNotCosmics.Sumw2()
    hNumHitsTracksNotCosmicFoundCosmics = ROOT.TH1F("hNumHitsTracksNotCosmicFoundCosmics ","",10,0,2000)
    hNumHitsTracksNotCosmicFoundCosmics.Sumw2()



    # Loop over events 
    for ientry in xrange(nEntries):
        nb = c.GetEntry(ientry)
        if nb <= 0:
           continue
     
        if ientry >= nevts:
            break

        longestTrackLen = 0
        longestTrackPhi = 0
        longestTrackTheta = 0
        longestTrackThetaXZ = 0
        longestTrackThetaYZ = 0
        longestTrackCosmicScore = 0
        longestTrackPIDA = 0
        longestTrackNumHits = 0
        foundLongestTrack = False

        hNumTracks.append(getattr(c,"ntracks_"+trackModule))

        # Loop over all tracks in event
        for kentry in xrange(getattr(c,"ntracks_"+trackModule)):

            # If the cosmic score makes no sense, continue
            if getattr(c,"trkcosmicscore_tagger_"+trackModule)[kentry] < -10:
                continue


            # Was this track a true cosmic?
            # Have to access this a little differently, the matrix a[i][j] translates to a[i*n_cols+j]
            # We also have to check the origin over all three planes and take the one with the most votes
            countsCosmic = 0
            isCosmic = False
            if getattr(c,"trkorigin_"+trackModule)[kentry*3+0] == 2:
                countsCosmic+=1
            if getattr(c,"trkorigin_"+trackModule)[kentry*3+1] == 2:
                countsCosmic+=1
            if getattr(c,"trkorigin_"+trackModule)[kentry*3+2] == 2:
                countsCosmic+=1
            # If isCosmic is > 1, then we have two votes for a cosmic, good enough for me
            if countsCosmic > 1:
                isCosmic = True


            # Does the track have a cosmic tag?
            hasCosmicTag = False
            if getattr(c,"trkcosmicscore_tagger_"+trackModule)[kentry] > 0.5 or ( getattr(c,"trkcosmicscore_tagger_"+trackModule)[kentry] > 0.3 and  getattr(c,"trkcosmicscore_flashmatch_"+trackModule)[kentry]):
            #if getattr(c,"trkcosmicscore_tagger_"+trackModule)[kentry] > 0.5:
            #if getattr(c,"trkcosmicscore_tagger_"+trackModule)[kentry] > 0.3:
                hasCosmicTag = True



            if isCosmic:
                hNumHitsTracksTrueCosmics.Fill(getattr(c,"ntrkhits_"+trackModule)[kentry])
                # Was this tagged as a cosmic?
                if hasCosmicTag:
                    hNumHitsTracksFoundCosmics.Fill(getattr(c,"ntrkhits_"+trackModule)[kentry])
            else:
                hNumHitsTracksTrueNotCosmics.Fill(getattr(c,"ntrkhits_"+trackModule)[kentry])
                # Was this tagged as a cosmic?
                if hasCosmicTag:
                    hNumHitsTracksNotCosmicFoundCosmics.Fill(getattr(c,"ntrkhits_"+trackModule)[kentry])

            # Veto the track if it has a cosmic tag
            #if getattr(c,"trkcosmicscore_tagger_"+trackModule)[kentry] > 0:
            #if getattr(c,"trkCosmicScore_"+trackModule)[kentry] > 0:
            if hasCosmicTag:
                continue













            # Veto the track if it has a pida > 10
            if getattr(c,"trkpidpida_"+trackModule)[kentry] > 10:
                continue

            if getattr(c,"trklen_"+trackModule)[kentry] > longestTrackLen:
                longestTrackLen = getattr(c,"trklen_"+trackModule)[kentry]
                longestTrackPhi = getattr(c,"trkphi_"+trackModule)[kentry]
                longestTrackTheta = getattr(c,"trktheta_"+trackModule)[kentry]
                longestTrackThetaXZ = getattr(c,"trkthetaxz_"+trackModule)[kentry]
                longestTrackThetaYZ = getattr(c,"trkthetayz_"+trackModule)[kentry]
                longestTrackCosmicScore = getattr(c,"trkcosmicscore_tagger_"+trackModule)[kentry]
                #longestTrackCosmicScore = getattr(c,"trkCosmicScore_"+trackModule)[kentry]
                longestTrackPIDA = getattr(c,"trkpidpida_"+trackModule)[kentry]
                longestTrackNumHits = getattr(c,"ntrkhits_"+trackModule)[kentry]
                foundLongestTrack = True

            hTrackPhi.append(getattr(c,"trkphi_"+trackModule)[kentry])
            hTrackTheta.append(getattr(c,"trktheta_"+trackModule)[kentry])
            hTrackThetaXZ.append(getattr(c,"trkthetaxz_"+trackModule)[kentry])
            hTrackThetaYZ.append(getattr(c,"trkthetayz_"+trackModule)[kentry])
            hTrackMom.append(getattr(c,"trkmom_"+trackModule)[kentry])
            hTrackLen.append(getattr(c,"trklen_"+trackModule)[kentry])
            hTrackCosmicScore.append(getattr(c,"trkcosmicscore_tagger_"+trackModule)[kentry])
            #hTrackCosmicScore.append(getattr(c,"trkCosmicScore_"+trackModule)[kentry])
            hTrackPIDA.append(getattr(c,"trkpidpida_"+trackModule)[kentry])
            hTrackNumHits.append(getattr(c,"ntrkhits_"+trackModule)[kentry])
       
        if foundLongestTrack:
            hLongestTrackLen.append(longestTrackLen)
            hLongestTrackPhi.append(longestTrackPhi)
            hLongestTrackTheta.append(longestTrackTheta)
            hLongestTrackThetaXZ.append(longestTrackThetaXZ)
            hLongestTrackThetaYZ.append(longestTrackThetaYZ)
            hLongestTrackCosmicScore.append(longestTrackCosmicScore)
            hLongestTrackPIDA.append(longestTrackPIDA)
            hLongestTrackNumHits.append(longestTrackNumHits)
   
    plt.figure()
    plt.hist(hNumTracks)
    plt.xlabel(r'Number of Tracks')
    plt.margins(y=0.05)
    plt.ylim(ymin=0)
    #pickle.dump(plt.gcf(),file('NumTracks.pickle','w'))
    plt.savefig('NumTracks.png')

    plt.figure()
    plt.hist(hTrackPhi, bins = 100, range = (-3.14,3.14))
    plt.xlim([-3.14,3.14])
    plt.xlabel(r'track $\phi$')
    plt.ylabel(r'tracks/(0.0628 radians)')
    plt.margins(y=0.05)
    plt.ylim(ymin=0)
    #pickle.dump(plt.gcf(),file('TrackPhi.pickle','w'))
    plt.savefig('TrackPhi.png')
    
    plt.figure()
    plt.hist(hLongestTrackPhi, bins = 100, range = (-3.14,3.14))
    plt.xlim([-3.14,3.14])
    plt.xlabel(r'longest track $\phi$')
    plt.ylabel(r'tracks/(0.0628 radians)')
    plt.margins(y=0.05)
    plt.ylim(ymin=0)
    #pickle.dump(plt.gcf(),file('LongestTrackPhi.pickle','w'))
    plt.savefig('LongestTrackPhi.png')
    
    plt.figure()
    plt.hist(hTrackTheta, bins = 100, range = (0,3.14))
    plt.xlim([0,3.14])
    plt.xlabel(r'track $\theta$')
    plt.ylabel(r'tracks/(0.0314 radians)')
    plt.margins(y=0.05)
    plt.ylim(ymin=0)
    #pickle.dump(plt.gcf(),file('TrackTheta.pickle','w'))
    plt.savefig('TrackTheta.png')
    
    plt.figure()
    plt.hist(hLongestTrackTheta, bins = 100, range = (0,3.14))
    plt.xlim([0,3.14])
    plt.xlabel(r'longest track $\theta$')
    plt.ylabel(r'tracks/(0.0314 radians)')
    plt.margins(y=0.05)
    plt.ylim(ymin=0)
    #pickle.dump(plt.gcf(),file('LongestTrackTheta.pickle','w'))
    plt.savefig('LongestTrackTheta.png')
    
    plt.figure()
    plt.hist(hTrackThetaXZ, bins = 100, range = (-3.14,3.14))
    plt.xlim([-3.14,3.14])
    plt.xlabel(r'track $\theta$ in the XZ plane')
    plt.margins(y=0.05)
    plt.ylim(ymin=0)
    #pickle.dump(plt.gcf(),file('TrackThetaXZ.pickle','w'))
    plt.savefig('TrackThetaXZ.png')
    
    plt.figure()
    plt.hist(hLongestTrackThetaXZ, bins = 100, range = (-3.14,3.14))
    plt.xlim([-3.14,3.14])
    plt.xlabel(r'longest track $\theta$ in the XZ plane')
    plt.margins(y=0.05)
    plt.ylim(ymin=0)
    #pickle.dump(plt.gcf(),file('LongestTrackThetaXZ.pickle','w'))
    plt.savefig('LongestTrackThetaXZ.png')
    
    plt.figure()
    plt.hist(hTrackThetaYZ, bins = 100, range = (0,3.14))
    plt.xlim([0,3.14])
    plt.xlabel(r'track $\theta$ in the YZ plane')
    plt.margins(y=0.05)
    plt.ylim(ymin=0)
    #pickle.dump(plt.gcf(),file('TrackThetaYZ.pickle','w'))
    plt.savefig('TrackThetaYZ.png')
    
    plt.figure()
    plt.hist(hLongestTrackThetaYZ, bins = 100, range = (0,3.14))
    plt.xlim([0,3.14])
    plt.xlabel(r'longest track $\theta$ in the YZ plane')
    plt.margins(y=0.05)
    plt.ylim(ymin=0)
    #pickle.dump(plt.gcf(),file('LongestTrackThetaYZ.pickle','w'))
    plt.savefig('LongestTrackThetaYZ.png')
    
    plt.figure()
    plt.hist2d(hLongestTrackPhi,hLongestTrackTheta,100)
    plt.colorbar()
    plt.xlim([-3.14,3.14])
    plt.ylim([0,3.14])
    plt.xlabel(r'longest track $\phi$')
    plt.ylabel(r'longest track $\theta$')
    #pickle.dump(plt.gcf(),file('LongestTrackPhiTheta.pickle','w'))
    plt.savefig('LongestTrackPhiTheta.png')
    
    plt.figure()
    plt.hist2d(hTrackPhi,hTrackTheta,100)
    plt.colorbar()
    plt.xlim([-3.14,3.14])
    plt.ylim([0,3.14])
    plt.xlabel(r'track $\phi$')
    plt.ylabel(r'track $\theta$')
    #pickle.dump(plt.gcf(),file('TrackPhiTheta.pickle','w'))
    plt.savefig('TrackPhiTheta.png')
    
    plt.figure()
    plt.hist(hTrackMom, bins = 100, range = (0,1))
    plt.xlim([0,1])
    plt.margins(y=0.05)
    plt.ylim(ymin=0)
    #pickle.dump(plt.gcf(),file('TrackMom.pickle','w'))
    plt.savefig('TrackMom.png')
    
    plt.figure()
    plt.hist(hTrackLen, bins = 100, range = (0,2000))
    plt.xlim([0,2000])
    plt.xlabel(r'track length (cm)')
    plt.margins(y=0.05)
    plt.ylim(ymin=0)
    #pickle.dump(plt.gcf(),file('TrackLen.pickle','w'))
    plt.savefig('TrackLen.png')
   
    plt.figure()
    plt.hist(hLongestTrackLen, bins = 100, range = (0,2000))
    plt.xlim([0,2000])
    plt.xlabel(r'longest track length (cm)')
    plt.margins(y=0.05)
    plt.ylim(ymin=0)
    #pickle.dump(plt.gcf(),file('LongestTrackLen.pickle','w'))
    plt.savefig('LongestTrackLen.png')

    plt.figure()
    plt.hist(hTrackCosmicScore, bins = 100, range = (-1.5,1.5))
    plt.xlim([-1.5,1.5])
    plt.xlabel(r'track cosmic score')
    plt.margins(y=0.05)
    plt.ylim(ymin=0)
    #pickle.dump(plt.gcf(),file('TrackCosmicScore.pickle','w'))
    plt.savefig('TrackCosmicScore.png')

    plt.figure()
    plt.hist(hLongestTrackCosmicScore, bins = 100, range = (-1.5,1.5))
    plt.xlim([-1.5,1.5])
    plt.xlabel(r'longest track cosmic score')
    plt.margins(y=0.05)
    plt.ylim(ymin=0)
    #pickle.dump(plt.gcf(),file('LongestTrackCosmicScore.pickle','w'))
    plt.savefig('LongestTrackCosmicScore.png')

    plt.figure()
    plt.hist(hTrackPIDA, bins = 100, range = (0.0,100.0))
    plt.xlim([0.0,100.0])
    plt.xlabel(r'track PIDA')
    plt.margins(y=0.05)
    plt.ylim(ymin=0)
    #pickle.dump(plt.gcf(),file('TrackPIDA.pickle','w'))
    plt.savefig('TrackPIDA.png')

    plt.figure()
    plt.hist(hLongestTrackPIDA, bins = 100, range = (0.0,100.0))
    plt.xlim([0.0,100.0])
    plt.xlabel(r'longest track PIDA')
    plt.margins(y=0.05)
    plt.ylim(ymin=0)
    #pickle.dump(plt.gcf(),file('LongestTrackPIDA.pickle','w'))
    plt.savefig('LongestTrackPIDA.png')

    plt.figure()
    plt.hist(hTrackNumHits, bins = 100, range = (0,2000))
    plt.xlim([0,2000])
    plt.xlabel(r'track number of hits')
    plt.margins(y=0.05)
    plt.ylim(ymin=0)
    #pickle.dump(plt.gcf(),file('TrackNumHits.pickle','w'))
    plt.savefig('TrackNumHits.png')
   
    plt.figure()
    plt.hist(hLongestTrackNumHits, bins = 100, range = (0,2000))
    plt.xlim([0,2000])
    plt.xlabel(r'longest track number of hits')
    plt.margins(y=0.05)
    plt.ylim(ymin=0)
    #pickle.dump(plt.gcf(),file('LongestTrackNumHits.pickle','w'))
    plt.savefig('LongestTrackNumHits.png')


    plt.figure()
    hNumHitsTracksFoundCosmicsR2M = r2m.Hist(hNumHitsTracksFoundCosmics)
    hNumHitsTracksFoundCosmicsR2M.hist()
    plt.margins(y=0.05)
    plt.ylim(ymin=0)
    plt.xlabel(r'track number of hits')
    #pickle.dump(plt.gcf(),file('NumHitsTracksFoundCosmics.pickle','w'))
    plt.savefig('NumHitsTracksFoundCosmics.png')

    plt.figure()
    hNumHitsTracksTrueCosmicsR2M = r2m.Hist(hNumHitsTracksTrueCosmics)
    hNumHitsTracksTrueCosmicsR2M.hist()
    plt.margins(y=0.05)
    plt.ylim(ymin=0)
    plt.xlabel(r'track number of hits')
    #pickle.dump(plt.gcf(),file('NumHitsTracksTrueCosmics.pickle','w'))
    plt.savefig('NumHitsTracksTrueCosmics.png')

    plt.figure()
    hCosmicTagEfficiencyNumHits = hNumHitsTracksFoundCosmics
    # Divide the histograms and compute the binomial errors for each bin (that's what the "B" is for)
    hCosmicTagEfficiencyNumHits.Divide(hNumHitsTracksFoundCosmics,hNumHitsTracksTrueCosmics,1.0,1.0,"B")
    hCosmicTagEfficiencyNumHitsR2M = r2m.Hist(hCosmicTagEfficiencyNumHits)
    hCosmicTagEfficiencyNumHitsR2M.errorbar(fmt='o', yerr=True, xerr=True, color='blue')
    plt.margins(y=0.05)
    plt.ylim(ymin=0)
    plt.xlabel(r'track number of hits')
    #pickle.dump(plt.gcf(),file('CosmicTagEfficiencyNumHits.pickle','w'))
    plt.savefig('CosmicTagEfficiencyNumHits.png')




    plt.figure()
    hNumHitsTracksNotCosmicFoundCosmicsR2M = r2m.Hist(hNumHitsTracksNotCosmicFoundCosmics)
    hNumHitsTracksNotCosmicFoundCosmicsR2M.hist()
    plt.margins(y=0.05)
    plt.ylim(ymin=0)
    plt.xlabel(r'track number of hits')
    #pickle.dump(plt.gcf(),file('NumHitsTracksNotCosmicFoundCosmics.pickle','w'))
    plt.savefig('NumHitsTracksNotCosmicFoundCosmics.png')
    
    plt.figure()
    hNumHitsTracksTrueNotCosmicsR2M = r2m.Hist(hNumHitsTracksTrueNotCosmics)
    hNumHitsTracksTrueNotCosmicsR2M.hist()
    plt.margins(y=0.05)
    plt.ylim(ymin=0)
    plt.xlabel(r'track number of hits')
    #pickle.dump(plt.gcf(),file('NumHitsTracksTrueNotCosmics.pickle','w'))
    plt.savefig('NumHitsTracksTrueNotCosmics.png')

    plt.figure()
    hCosmicMistagRateNumHits = hNumHitsTracksNotCosmicFoundCosmics
    # Divide the histograms and compute the binomial errors for each bin (that's what the "B" is for)
    hCosmicMistagRateNumHits.Divide(hCosmicMistagRateNumHits,hNumHitsTracksTrueNotCosmics,1.0,1.0,"B")
    hCosmicMistagRateNumHitsR2M = r2m.Hist(hCosmicMistagRateNumHits)
    hCosmicMistagRateNumHitsR2M.errorbar(fmt='o', yerr=True, xerr=True, color='blue')
    plt.margins(y=0.05)
    plt.ylim(ymin=0)
    plt.xlabel(r'track number of hits')
    #pickle.dump(plt.gcf(),file('CosmicMistagRateNumHits.pickle','w'))
    plt.savefig('CosmicMistagRateNumHits.png')




    totalEfficiency = 0
    totalEfficiencyError = 0
    numNonZeroBins = 0
    for binx in xrange(hCosmicTagEfficiencyNumHits.GetNbinsX()):
        if hNumHitsTracksTrueNotCosmics.GetBinContent(binx) == 0:
            continue
        numNonZeroBins+=1
        totalEfficiency+= hCosmicTagEfficiencyNumHits.GetBinContent(binx)
        totalEfficiencyError+=hCosmicTagEfficiencyNumHits.GetBinError(binx)*hCosmicTagEfficiencyNumHits.GetBinError(binx)

    print 'The total cosmic tagging efficiency is '+str(totalEfficiency/numNonZeroBins)
    print 'The total cosmic tagging efficiency error is '+str(math.sqrt(totalEfficiencyError)/numNonZeroBins)
    


    totalMistag = 0
    totalMistagError= 0
    numNonZeroBins = 0
    totalTracksNotCosmicFoundCosmics = 0
    totalTracksTrueNotCosmics = 0
    for binx in xrange(hCosmicMistagRateNumHits.GetNbinsX()):
        if hNumHitsTracksTrueNotCosmics.GetBinContent(binx) == 0:
            continue
        numNonZeroBins+= 1
        totalMistag+=hCosmicMistagRateNumHits.GetBinContent(binx)
        totalMistagError+=hCosmicMistagRateNumHits.GetBinError(binx)*hCosmicMistagRateNumHits.GetBinError(binx)


    print 'The total cosmic tagging mistag rate is '+str(totalMistag/numNonZeroBins)
    print 'The total cosmic tagging mistag error is '+str(math.sqrt(totalMistagError)/numNonZeroBins)





if __name__ == '__main__':
    sys.exit(main(sys.argv))









