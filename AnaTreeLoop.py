#!/usr/local/bin/python

import ROOT
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

    # Track module name
    trackModule = "trackkalmanhitpandoraneutrino"
    #trackModule = "trackkalmanhitpandora"
    #trackModule = "trackkalmanhit"

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
        foundLongestTrack = False

        hNumTracks.append(getattr(c,"ntracks_"+trackModule))

        # Loop over tracks in event
        for kentry in xrange(getattr(c,"ntracks_"+trackModule)):

            # Veto the track if it has a cosmic tag
            if getattr(c,"trkcosmicscore_tagger_"+trackModule)[kentry] > 0:
            #if getattr(c,"trkCosmicScore_"+trackModule)[kentry] > 0:
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
       
        if foundLongestTrack:
            hLongestTrackLen.append(longestTrackLen)
            hLongestTrackPhi.append(longestTrackPhi)
            hLongestTrackTheta.append(longestTrackTheta)
            hLongestTrackThetaXZ.append(longestTrackThetaXZ)
            hLongestTrackThetaYZ.append(longestTrackThetaYZ)
            hLongestTrackCosmicScore.append(longestTrackCosmicScore)
            hLongestTrackPIDA.append(longestTrackPIDA)
   
    plt.figure()
    plt.hist(hNumTracks)
    plt.xlabel(r'Number of Tracks')
    pickle.dump(plt.gcf(),file('NumTracks.pickle','w'))
    plt.savefig('NumTracks.png')

    plt.figure()
    plt.hist(hTrackPhi, bins = 100, range = (-3.14,3.14))
    plt.xlim([-3.14,3.14])
    plt.xlabel(r'track $\phi$')
    plt.ylabel(r'tracks/(0.0628 radians)')
    pickle.dump(plt.gcf(),file('TrackPhi.pickle','w'))
    plt.savefig('TrackPhi.png')
    
    plt.figure()
    plt.hist(hLongestTrackPhi, bins = 100, range = (-3.14,3.14))
    plt.xlim([-3.14,3.14])
    plt.xlabel(r'longest track $\phi$')
    plt.ylabel(r'tracks/(0.0628 radians)')
    pickle.dump(plt.gcf(),file('LongestTrackPhi.pickle','w'))
    plt.savefig('LongestTrackPhi.png')
    
    plt.figure()
    plt.hist(hTrackTheta, bins = 100, range = (0,3.14))
    plt.xlim([0,3.14])
    plt.xlabel(r'track $\theta$')
    plt.ylabel(r'tracks/(0.0314 radians)')
    pickle.dump(plt.gcf(),file('TrackTheta.pickle','w'))
    plt.savefig('TrackTheta.png')
    
    plt.figure()
    plt.hist(hLongestTrackTheta, bins = 100, range = (0,3.14))
    plt.xlim([0,3.14])
    plt.xlabel(r'longest track $\theta$')
    plt.ylabel(r'tracks/(0.0314 radians)')
    pickle.dump(plt.gcf(),file('LongestTrackTheta.pickle','w'))
    plt.savefig('LongestTrackTheta.png')
    
    plt.figure()
    plt.hist(hTrackThetaXZ, bins = 100, range = (-3.14,3.14))
    plt.xlim([-3.14,3.14])
    plt.xlabel(r'track $\theta$ in the XZ plane')
    pickle.dump(plt.gcf(),file('TrackThetaXZ.pickle','w'))
    plt.savefig('TrackThetaXZ.png')
    
    plt.figure()
    plt.hist(hLongestTrackThetaXZ, bins = 100, range = (-3.14,3.14))
    plt.xlim([-3.14,3.14])
    plt.xlabel(r'longest track $\theta$ in the XZ plane')
    pickle.dump(plt.gcf(),file('LongestTrackThetaXZ.pickle','w'))
    plt.savefig('LongestTrackThetaXZ.png')
    
    plt.figure()
    plt.hist(hTrackThetaYZ, bins = 100, range = (0,3.14))
    plt.xlim([0,3.14])
    plt.xlabel(r'track $\theta$ in the YZ plane')
    pickle.dump(plt.gcf(),file('TrackThetaYZ.pickle','w'))
    plt.savefig('TrackThetaYZ.png')
    
    plt.figure()
    plt.hist(hLongestTrackThetaYZ, bins = 100, range = (0,3.14))
    plt.xlim([0,3.14])
    plt.xlabel(r'longest track $\theta$ in the YZ plane')
    pickle.dump(plt.gcf(),file('LongestTrackThetaYZ.pickle','w'))
    plt.savefig('LongestTrackThetaYZ.png')
    
    plt.figure()
    plt.hist2d(hLongestTrackPhi,hLongestTrackTheta,100)
    plt.colorbar()
    plt.xlim([-3.14,3.14])
    plt.ylim([0,3.14])
    plt.xlabel(r'longest track $\phi$')
    plt.ylabel(r'longest track $\theta$')
    pickle.dump(plt.gcf(),file('LongestTrackPhiTheta.pickle','w'))
    plt.savefig('LongestTrackPhiTheta.png')
    
    plt.figure()
    plt.hist2d(hTrackPhi,hTrackTheta,100)
    plt.colorbar()
    plt.xlim([-3.14,3.14])
    plt.ylim([0,3.14])
    plt.xlabel(r'track $\phi$')
    plt.ylabel(r'track $\theta$')
    pickle.dump(plt.gcf(),file('TrackPhiTheta.pickle','w'))
    plt.savefig('TrackPhiTheta.png')
    
    plt.figure()
    plt.hist(hTrackMom, bins = 100, range = (0,1))
    plt.xlim([0,1])
    pickle.dump(plt.gcf(),file('TrackMom.pickle','w'))
    plt.savefig('TrackMom.png')
    
    plt.figure()
    plt.hist(hTrackLen, bins = 100, range = (0,2000))
    plt.xlim([0,2000])
    plt.xlabel(r'track length (cm)')
    pickle.dump(plt.gcf(),file('TrackLen.pickle','w'))
    plt.savefig('TrackLen.png')
   
    plt.figure()
    plt.hist(hLongestTrackLen, bins = 100, range = (0,2000))
    plt.xlim([0,2000])
    plt.xlabel(r'longest track length (cm)')
    pickle.dump(plt.gcf(),file('LongestTrackLen.pickle','w'))
    plt.savefig('LongestTrackLen.png')

    plt.figure()
    plt.hist(hTrackCosmicScore, bins = 100, range = (-1.5,1.5))
    plt.xlim([-1.5,1.5])
    plt.xlabel(r'track cosmic score')
    pickle.dump(plt.gcf(),file('TrackCosmicScore.pickle','w'))
    plt.savefig('TrackCosmicScore.png')

    plt.figure()
    plt.hist(hLongestTrackCosmicScore, bins = 100, range = (-1.5,1.5))
    plt.xlim([-1.5,1.5])
    plt.xlabel(r'longest track cosmic score')
    pickle.dump(plt.gcf(),file('LongestTrackCosmicScore.pickle','w'))
    plt.savefig('LongestTrackCosmicScore.png')

    plt.figure()
    plt.hist(hTrackPIDA, bins = 100, range = (0.0,100.0))
    plt.xlim([0.0,100.0])
    plt.xlabel(r'track PIDA')
    pickle.dump(plt.gcf(),file('TrackPIDA.pickle','w'))
    plt.savefig('TrackPIDA.png')

    plt.figure()
    plt.hist(hLongestTrackPIDA, bins = 100, range = (0.0,100.0))
    plt.xlim([0.0,100.0])
    plt.xlabel(r'longest track PIDA')
    pickle.dump(plt.gcf(),file('LongestTrackPIDA.pickle','w'))
    plt.savefig('LongestTrackPIDA.png')


if __name__ == '__main__':
    sys.exit(main(sys.argv))









