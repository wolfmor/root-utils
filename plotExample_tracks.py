import ROOT
from utils import PlotFactory, HistoSample, TreeSample, Variable,  Ratio

pf = PlotFactory(
    outputpath='.',  # change this to a path that suits you
    outputpattern='compareStopDirect_VARIABLE',  # change this to a name that suits you (VARIABLE will be replaced by the variable name)
    outputformat='pdf',

    normalize=True,  # for now, we only want to look at differences in shape
    ylabel='Fraction of Tracks',
    ylabelratio='Stop/Direct',

    yaxisrangeratio=(0.0001, 2.9999),
    ratiohlines=[1],

    text='2016',
)

path = '/pnfs/desy.de/cms/tier2/store/user/mowolf/NTuples/NTuplesV10/'  # the path where the ntuples are stored
eventselection = None  # don't apply any additional event selection cuts
trackselection = 'track_isSignalTrack==1'  # we only want to look at our signal track

ntestfiles = 0  # if it takes too long to process all files, use only a few to test

pf.add_samples([
    # add two direct higgsino production samples
    TreeSample(ntestfiles=ntestfiles, category='line', name='direct_mChi115_dm0p57',
               tree='tEvent', files=path + 'SignalV4_16/step3_higgsinoDm0eDmpm_RunIISpring21UL16FS_susyall_mChipm115GeV_dm0p568GeV_part*_NTuple_job*.root',
               eventselection=eventselection, vectorselection=trackselection,
               color=ROOT.kAzure+1),
    TreeSample(ntestfiles=ntestfiles, category='line', name='direct_mChi115_dm0p97',
               tree='tEvent', files=path + 'SignalV4_16/step3_higgsinoDm0eDmpm_RunIISpring21UL16FS_susyall_mChipm115GeV_dm0p968GeV_part*_NTuple_job*.root',
               eventselection=eventselection, vectorselection=trackselection,
               color=ROOT.kRed),

    # add the corresponding stop production samples
    TreeSample(ntestfiles=ntestfiles, category='marker', name='stop_mStop800_mChi115_dm0p6',
               tree='tEvent', files=path + 'SignalStopV3_16/higgsino_Summer16_stopstop_800GeV_mChipm115GeV_dm0p6GeV_pu35_part*_NTuple_job*.root',
               eventselection=eventselection, vectorselection=trackselection,
               color=ROOT.kAzure+1),
    TreeSample(ntestfiles=ntestfiles, category='marker', name='stop_mStop800_mChi115_dm1p0',
               tree='tEvent', files=path + 'SignalStopV3_16/higgsino_Summer16_stopstop_800GeV_mChipm115GeV_dm1p0GeV_pu35_part*_NTuple_job*.root',
               eventselection=eventselection, vectorselection=trackselection,
               color=ROOT.kRed),

    # to compare the signal tracks with all other tracks in the event (note the different "vectorselection")
    TreeSample(ntestfiles=ntestfiles, category='line', name='direct_insignalBkg',
               tree='tEvent', files=path + 'SignalV4_16/step3_higgsinoDm0eDmpm_RunIISpring21UL16FS_susyall_mChipm*GeV_dm*GeV_part1of25_NTuple_job*.root',  # it's enough to only use the *part1* file of each model point
               eventselection=eventselection, vectorselection='track_isSignalTrack==0&&track_isSusyTrack==0',
               color=ROOT.kBlack),

    TreeSample(ntestfiles=ntestfiles, category='marker', name='stop_insignalBkg',
               tree='tEvent', files=path + 'SignalStopV3_16/higgsino_Summer16_stopstop_*GeV_mChipm*GeV_dm*GeV_pu35_part1of25_NTuple_job*.root',  # it's enough to only use the *part1* file of each model point
               eventselection=eventselection, vectorselection='track_isSignalTrack==0&&track_isSusyTrack==0',
               color=ROOT.kBlack),
])

pf.add_variables([
    # some track observables
    Variable(name='track_pt', axisrange=(0, 5), nbins=100),
    Variable(name='track_log10_dxy_', axisrange=(-5, 1), nbins=60),
    Variable(name='track_log10_dz_', axisrange=(-5, 1), nbins=60),
    Variable(name='track_dphiMet', axisrange=(-3.5, 3.5), nbins=70),

    Variable(name='track_drminJet30', axisrange=(0., 10), nbins=50),

    # add more! e.g. isolation-related variables, impact parameter (IP),...
])

pf.add_ratios([
    # to visualize the difference between the production modes
    Ratio(category='ratio', name='stop_mStop800_mChi115_dm0p6:direct_mChi115_dm0p57'),
    Ratio(category='ratio', name='stop_mStop800_mChi115_dm1p0:direct_mChi115_dm0p97'),
])

# make the plots!
# (might take some time depending on the size of the input trees)
pf.process()
