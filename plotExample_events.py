import ROOT
from utils import PlotFactory, HistoSample, TreeSample, Variable,  Ratio

pf = PlotFactory(
    outputpath='.',  # change this to a path that suits you
    outputpattern='compareStopDirect_VARIABLE',  # change this to a name that suits you (VARIABLE will be replaced by the variable name)
    outputformat='pdf',

    normalize=True,  # for now, we only want to look at differences in shape
    ylabel='Fraction of Events',
    ylabelratio='Stop/Direct',

    yaxisrangeratio=(0.0001, 2.9999),
    ratiohlines=[1],

    text='2016',
)

path = '/pnfs/desy.de/cms/tier2/store/user/mowolf/NTuples/NTuplesV10/'  # the path where the ntuples are stored
eventselection = None  # don't apply any additional event selection cuts

ntestfiles = 0  # if it takes too long to process all files, use only a few to test

pf.add_samples([
    # add two direct higgsino production samples
    TreeSample(ntestfiles=ntestfiles, category='line', name='direct_mChi115_dm0p57',
               tree='tEvent', files=path + 'SignalV4_16/step3_higgsinoDm0eDmpm_RunIISpring21UL16FS_susyall_mChipm115GeV_dm0p568GeV_part*_NTuple_job*.root',
               eventselection=eventselection,
               color=ROOT.kAzure+1),
    TreeSample(ntestfiles=ntestfiles, category='line', name='direct_mChi115_dm0p97',
               tree='tEvent', files=path + 'SignalV4_16/step3_higgsinoDm0eDmpm_RunIISpring21UL16FS_susyall_mChipm115GeV_dm0p968GeV_part*_NTuple_job*.root',
               eventselection=eventselection,
               color=ROOT.kRed),

    # add the corresponding stop production samples
    TreeSample(ntestfiles=ntestfiles, category='marker', name='stop_mStop800_mChi115_dm0p6',
               tree='tEvent', files=path + 'SignalStopV3_16/higgsino_Summer16_stopstop_800GeV_mChipm115GeV_dm0p6GeV_pu35_part*_NTuple_job*.root',
               eventselection=eventselection,
               color=ROOT.kAzure+1),
    TreeSample(ntestfiles=ntestfiles, category='marker', name='stop_mStop800_mChi115_dm1p0',
               tree='tEvent', files=path + 'SignalStopV3_16/higgsino_Summer16_stopstop_800GeV_mChipm115GeV_dm1p0GeV_pu35_part*_NTuple_job*.root',
               eventselection=eventselection,
               color=ROOT.kRed),
])

pf.add_variables([
    # first some GEN information
    Variable(name='n_chiC1', axisrange=(0, 3), nbins=3),
    Variable(name='n_chiN2', axisrange=(0, 3), nbins=3),
    Variable(name='chiC1_log10_decaylengthXY_', axisrange=(-5, 1), nbins=60),
    Variable(name='chiC1_log10_chidecaylengthZ_', axisrange=(-5, 1), nbins=60),
    Variable(name='chiC1_pionPt', axisrange=(-1, 9), nbins=100),
    Variable(name='chiC1_hasMatchedTrackPion', axisrange=(0, 2), nbins=2),

    # then some event observables
    Variable(name='met_pt', axisrange=(0, 1000), nbins=100),
    Variable(name='n_jet_100', axisrange=(0, 5), nbins=5),
    Variable(name='n_jet_30_btagloose', axisrange=(0, 10), nbins=10),
    Variable(name='n_jet_30_btagmedium', axisrange=(0, 10), nbins=10),
    Variable(name='n_jet_30_btagtight', axisrange=(0, 10), nbins=10),

    # add more! e.g. number of (isolated) leptons, HT(miss),...
])

pf.add_ratios([
    # to visualize the difference between the production modes
    Ratio(category='ratio', name='stop_mStop800_mChi115_dm0p6:direct_mChi115_dm0p57'),
    Ratio(category='ratio', name='stop_mStop800_mChi115_dm1p0:direct_mChi115_dm0p97'),
])

# make the plots!
# (might take some time depending on the size of the input trees)
pf.process()
