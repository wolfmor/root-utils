import ROOT

import re
from glob import glob
import sys
sys.path.insert(0, '/afs/desy.de/user/w/wolfmor/cmssw/CMSSW_10_2_18/src/SoftDisplacedTrack/root-utils')
sys.path.insert(0, '/afs/desy.de/user/w/wolfmor/cmssw/CMSSW_10_6_34/src/SoftDisplacedPion/ntuplizer')
sys.path.insert(0, '/data/dust/user/wolfmor/NTupleFileLists')

from utils import PlotFactory, HistoSample, TreeSample, Variable,  Ratio
from myvariables import myvariablesV14, myvariables_V15
from mydatasets import mydatasets_V14, mydatasets_V15


# setuplcg102centos9

print(sys.argv)
args = {}
if len(sys.argv) > 1:
    args['era'] = sys.argv[1]
    args['eventclass'] = sys.argv[2]
    args['dataset'] = sys.argv[3]
    args['outputsubsubfolder'] = sys.argv[4]

# see https://cms-analysis.docs.cern.ch/guidelines/plotting/colors
colors = [
    ROOT.TColor.GetColor('#3f90da'),
    ROOT.TColor.GetColor('#ffa90e'),
    ROOT.TColor.GetColor('#bd1f01'),
    ROOT.TColor.GetColor('#94a4a2'),
    ROOT.TColor.GetColor('#832db6'),
    ROOT.TColor.GetColor('#a96b59'),
    ROOT.TColor.GetColor('#e76300'),
    ROOT.TColor.GetColor('#b9ac70'),
    ROOT.TColor.GetColor('#717581'),
    ROOT.TColor.GetColor('#92dadd'),
]

# TODO: what ntuple/friends version and era(s)?
# if making signal histos, make 16APV and 16 and 17/18 separately due to different PU weights
tag = 'NTuplesV15'
friends = 'NTuplesV15p15'  # None  #
variables = myvariables_V15
datasets = mydatasets_V15
eras = [
    'era16_UL_APV',
    'era16_UL',
    'era17_UL',
    'era18_UL',
]
if len(args) > 0:
    eras = [args['era']]

emptyntupletxts = '/data/dust/user/wolfmor/EmptyNTupleTxts/' + tag + '/ERA/DATASET.txt'

# TODO: refiner?
# refiner = '_refined20241130'  # ''  #
refiners = [
    '',
    # '_refined20241221_2',
    # '_refined20250102',
    # '_refined20250102_1',
    # '_refined20250102_2',
    # '_refined20250102_3',
    '_refined20241221_2_ensemble',
]
for refiner in refiners:

    nn_names = [
        # 'PyKeras_V14_20231025_1_multiclass',
        # 'PyKeras_V14_20231121_multiclass',
        # 'PyKeras_V14_20240228_1_multiclass',
        # 'PyKeras_V14_20240319_4_multiclass',
        # 'PyKeras_V14_20240327_multiclass',
        # 'PyKeras_V14_20240328_multiclass',
        # 'PyKeras_V14_20240323_multiclass',
        # 'PyKeras_V14_20240403_multiclass',
        # 'PyKeras_V14_20240403_1_multiclass',
        # 'PyKeras_V14_20240403_2_multiclass',
        # 'PyKeras_V14_20240405_1_multiclass',
        # 'PyKeras_V14_20240408_multiclass',
        # 'PyKeras_V14_20240408_1_multiclass',
        # 'PyKeras_V14_20240409_multiclass',
        # 'PyKeras_V14_20240409_1_multiclass',
        # 'PyKeras_V14_20240411_multiclass',
        # 'PyKeras_V14_20240411_1_multiclass',
        'PyKeras_V15_20240711_multiclass',
    ]
    notparametrized = [
        'PyKeras_V14_20240328_multiclass',
    ]

    nn_names = [n + refiner for n in nn_names]

    nn_name = nn_names[0]  # take the first one for the selection cuts

    # SRs for V15_20240711: 0.3: >10.5, 0.6: >9.5, 1.0: >8.0

    nn_vars = {
        'PyKeras_V15_20240711_multiclass': ['track_pt', 'track_abs_eta_', 'track_distPVAssPVxy', 'track_distPVAssPVz', 'track_distPVAssSVxy', 'track_distPVAssSVz', 'track_log10_IPsigXY_', 'track_log10_IPsigZ_', 'track_log10_IPxy_', 'track_log10_IPz_', 'track_log10_IPsigXYPU_', 'track_log10_IPsigZPU_', 'track_log10_IPxyPU_', 'track_log10_IPzPU_', 'track_log10_IPsigXYAssPV_', 'track_log10_IPsigZAssPV_', 'track_log10_IPxyAssPV_', 'track_log10_IPzAssPV_', 'track_log10_IPsigXYPUAssPV_', 'track_log10_IPsigZPUAssPV_', 'track_log10_IPxyPUAssPV_', 'track_log10_IPzPUAssPV_', 'track_log10_dxy_', 'track_log10_dz_', 'track_log10_dxyPU_', 'track_log10_dzPU_', 'track_log10_dxyError_', 'track_log10_dzError_', 'track_pfAbsIso', 'track_drminTrack5', 'track_drmin2ndTrack5', 'track_drminJet15', 'track_drminJet30', 'track_abs_detaLeadingJet_', 'track_abs_dphiLeadingJet_', 'track_abs_dphiMet_', 'met_pt']
    }

    refiner_vars = {
        '_refined20241101': ['track_pt', 'track_drminJet15', 'track_log10_dzError_', 'track_log10_dxyError_', 'track_log10_IPsigZPUAssPV_'],
        '_refined20241111_1': ['track_pt', 'track_drminJet15', 'track_log10_dzError_', 'track_log10_dxyError_', 'track_log10_IPsigZPUAssPV_'],
        '_refined20241117_3': ['track_pt', 'track_drminJet15', 'track_log10_dzError_', 'track_log10_dxyError_', 'track_log10_IPsigZPUAssPV_'],
        '_refined20241117_4': ['track_pt', 'track_drminJet15', 'track_log10_dzError_', 'track_log10_dxyError_', 'track_log10_IPsigZPUAssPV_'],
        '_refined20241117_5': ['track_pt', 'track_drminJet15', 'track_log10_dzError_', 'track_log10_dxyError_', 'track_log10_IPsigZPUAssPV_'],
        '_refined20241120_1': ['track_pt', 'track_log10_dzError_', 'track_log10_dxyError_'],
        '_refined20241120_2': ['track_log10_dzError_', 'track_log10_dxyError_'],
        '_refined20241121': ['track_log10_dzError_', 'track_log10_dxyError_',
                             'track_log10_IPsigXY_', 'track_log10_IPsigZ_', 'track_log10_IPsigXYPU_', 'track_log10_IPsigZPU_',
                             'track_log10_IPsigXYAssPV_', 'track_log10_IPsigZAssPV_', 'track_log10_IPsigXYPUAssPV_', 'track_log10_IPsigZPUAssPV_'],
        '_refined20241121_1': ['track_log10_dzError_', 'track_log10_dxyError_',
                               'track_log10_IPsigXY_', 'track_log10_IPsigZ_', 'track_log10_IPsigXYPU_', 'track_log10_IPsigZPU_',
                               'track_log10_IPsigXYAssPV_', 'track_log10_IPsigZAssPV_', 'track_log10_IPsigXYPUAssPV_', 'track_log10_IPsigZPUAssPV_'],
        '_refined20241124': ['track_log10_dzError_', 'track_log10_dxyError_',
                             'track_log10_IPsigZPUAssPV_', 'track_log10_IPsigXY_', 'track_log10_IPsigZ_'],
        '_refined20241124_1': ['track_log10_dzError_', 'track_log10_dxyError_',
                               'track_log10_IPsigZPUAssPV_', 'track_log10_IPsigXY_', 'track_log10_IPsigZ_'],
        '_refined20241126': ['track_pt', 'track_log10_dzError_', 'track_log10_dxyError_'],
        '_refined20241128': ['track_pt', 'track_log10_dzError_', 'track_log10_dxyError_'],
        '_refined20241130': ['track_pt', 'track_log10_dzError_', 'track_log10_dxyError_'],
        '_refined20241204_3': ['track_pt', 'track_log10_dzError_', 'track_log10_dxyError_'],
        '_refined20241205_3': ['track_pt', 'track_log10_dzError_', 'track_log10_dxyError_'],
        '_refined20241205_4': ['track_pt', 'track_log10_dzError_', 'track_log10_dxyError_'],
        '_refined20241205_5': ['track_pt', 'track_log10_dzError_', 'track_log10_dxyError_'],
        '_refined20241206': ['track_pt', 'track_log10_dzError_', 'track_log10_dxyError_', 'track_log10_IPsigZPUAssPV_'],
        '_refined20241206_1': ['track_pt', 'track_log10_dzError_', 'track_log10_dxyError_', 'track_log10_IPsigZPUAssPV_'],
        '_refined20241206_2': ['track_pt', 'track_log10_dzError_', 'track_log10_dxyError_', 'track_log10_IPsigZPUAssPV_'],
        '_refined20241209': ['track_pt', 'track_log10_dzError_', 'track_log10_dxyError_'],
        '_refined20241209_1': ['track_pt', 'track_log10_dzError_', 'track_log10_dxyError_'],
        '_refined20241210': ['track_pt', 'track_log10_dzError_', 'track_log10_dxyError_', 'track_log10_IPsigZPUAssPV_'],
        '_refined20241210_1': ['track_pt', 'track_log10_dzError_', 'track_log10_dxyError_', 'track_log10_IPsigZPUAssPV_'],
        '_refined20241221': ['track_pt', 'track_log10_dzError_', 'track_log10_dxyError_',
                             'track_log10_IPsigXY_', 'track_log10_IPsigZ_', 'track_log10_IPsigXYPU_', 'track_log10_IPsigZPU_',
                             'track_log10_IPsigXYAssPV_', 'track_log10_IPsigZAssPV_', 'track_log10_IPsigXYPUAssPV_', 'track_log10_IPsigZPUAssPV_'],
        '_refined20241221_1': ['track_pt', 'track_log10_dzError_', 'track_log10_dxyError_', 'track_log10_IPsigZPUAssPV_'],
        '_refined20241221_2': ['track_pt', 'track_log10_dzError_', 'track_log10_dxyError_'],
        '_refined20241221_3': ['track_pt', 'track_log10_dzError_', 'track_log10_dxyError_', 'track_log10_IPsigZPUAssPV_'],
        '_refined20241229': ['track_pt', 'track_log10_dzError_', 'track_log10_dxyError_', 'track_log10_IPsigZPUAssPV_'],
        '_refined20241229_1': ['track_pt', 'track_log10_dzError_', 'track_log10_dxyError_'],
        '_refined20250102': ['track_pt', 'track_log10_dzError_', 'track_log10_dxyError_'],
        '_refined20250102_1': ['track_pt', 'track_log10_dzError_', 'track_log10_dxyError_'],
        '_refined20250102_2': ['track_pt', 'track_log10_dzError_', 'track_log10_dxyError_'],
        '_refined20250102_3': ['track_pt', 'track_log10_dzError_', 'track_log10_dxyError_'],
        '_refined20241221_2_ensemble': ['track_pt', 'track_log10_dzError_', 'track_log10_dxyError_'],
    }

    # UL (from AN2021_197_v10)
    lumi_DataMET_16Bver1 = 0.
    lumi_DataMET_16Bver2 = 5.826
    lumi_DataMET_16C_HIPM = 2.621
    lumi_DataMET_16D_HIPM = 4.286
    lumi_DataMET_16E_HIPM = 4.066
    lumi_DataMET_16F_HIPM = 2.865
    lumi_DataMET_16F = 0.584
    lumi_DataMET_16G = 7.653
    lumi_DataMET_16H = 8.740

    # from https://twiki.cern.ch/twiki/bin/viewauth/CMS/SUSRecommendationsRun2UltraLegacy
    # 2016: 35.9/fb, 2017: 41.5/fb, 2018: 59.7/fb, total: 137/fb
    lumis = {
        # unfortunately those don't add up to 35.9 but to 19.664+16.977=36.641, anyway...
        'era16_UL_APV': lumi_DataMET_16Bver1 + lumi_DataMET_16Bver2 + lumi_DataMET_16C_HIPM + lumi_DataMET_16D_HIPM + lumi_DataMET_16E_HIPM + lumi_DataMET_16F_HIPM,
        'era16_UL': lumi_DataMET_16F + lumi_DataMET_16G + lumi_DataMET_16H,
        'era17_UL': 41.5,
        'era18_UL': 59.8,
    }
    lumi = sum([lumis[era] for era in lumis if era in eras])


    # ####
    # event-level
    # ####

    noselection = '1'
    basicselection = 'met_pt>250&&ht>250&&htMiss>250&&ht5/ht<2&&n_jet_100>0&&dphiminMetJets>0.5&&n_photon_iso<1&&n_lepton_iso<1&&n_tau_20_tight<1&&n_jet_HEM1516veto<1&&n_jet_30<5'
    basicselection300 = 'met_pt>300&&htMiss>300&&ht5/ht<2&&ht>100&&n_jet_100>0&&dphiminMetJets>0.5&&n_photon_iso<1&&n_lepton_iso<1&&n_tau_20_tight<1&&n_jet_HEM1516veto<1&&n_jet_30<5'
    nobtagselection = basicselection + '&&n_jet_30_btagDeepCSVmedium<1'
    nobtagselection300 = basicselection300 + '&&n_jet_30_btagDeepCSVmedium<1'

    controlselection0p9 = nobtagselection300 + '&&maxscore_' + nn_name + '_Signal_deltam0p3<0.95' + '&&maxscore_' + nn_name + '_Signal_deltam0p6<0.95' + '&&maxscore_' + nn_name + '_Signal_deltam1p0<0.95'

    diagonalcut0p3vs1p0selection = nobtagselection300 + '&&maxscore_' + nn_name + '_Signal_deltam1p0<maxscore_' + nn_name + '_Signal_deltam0p3'
    diagonalcut0p6vs1p0selection = nobtagselection300 + '&&maxscore_' + nn_name + '_Signal_deltam1p0<maxscore_' + nn_name + '_Signal_deltam0p6'
    diagonalcut1p0vs0p3selection = nobtagselection300 + '&&maxscore_' + nn_name + '_Signal_deltam0p3<maxscore_' + nn_name + '_Signal_deltam1p0'
    diagonalcut1p0vs0p6selection = nobtagselection300 + '&&maxscore_' + nn_name + '_Signal_deltam0p6<maxscore_' + nn_name + '_Signal_deltam1p0'

    diagonalcut0p3vs0p6and1p0selection = nobtagselection300 + '&&maxscore_' + nn_name + '_Signal_deltam0p6<maxscore_' + nn_name + '_Signal_deltam0p3' + '&&maxscore_' + nn_name + '_Signal_deltam1p0<maxscore_' + nn_name + '_Signal_deltam0p3'
    diagonalcut0p6vs0p3and1p0selection = nobtagselection300 + '&&maxscore_' + nn_name + '_Signal_deltam0p3<maxscore_' + nn_name + '_Signal_deltam0p6' + '&&maxscore_' + nn_name + '_Signal_deltam1p0<maxscore_' + nn_name + '_Signal_deltam0p6'
    diagonalcut1p0vs0p3and0p6selection = nobtagselection300 + '&&maxscore_' + nn_name + '_Signal_deltam0p3<maxscore_' + nn_name + '_Signal_deltam1p0' + '&&maxscore_' + nn_name + '_Signal_deltam0p6<maxscore_' + nn_name + '_Signal_deltam1p0'

    srlikecut0p6 = diagonalcut0p6vs0p3and1p0selection + '&&maxscore_' + nn_name + '_Signal_deltam0p6>(1/(1+exp(-7)))'
    srcut0p6 = diagonalcut0p6vs0p3and1p0selection + '&&maxscore_' + nn_name + '_Signal_deltam0p6>(1/(1+exp(-8.75)))'

    # TODO: event-level selection
    thesrdm = ''  # '0p3' '0p6' '1p0'
    selection = nobtagselection300  # noselection  # nobtagselection300  #
    requiremctrigger = False

    # TODO: group events by classes?
    groupeventsbyclasses = False
    # TODO: uncomment event classes
    # TODO: adapt deltam in cutstring!
    eventclasses = [
        # TODO: use this instead?
        ('fromtruetau', 'Sec. from W(#tau)', '(maxscore_' + nn_name + '_Signal_deltam0p3_genMatchMotherIsTheTau==1)&&', colors[0]),

        # ('fromtruetaulep', 'Sec. from W(#tau #rightarrow e,mu)', '(maxscore_' + nn_name + '_Signal_deltam0p3_genMatchMotherIsTheTau==1)&&(maxscore_' + nn_name + '_Signal_deltam0p3_genMatchMotherTauDecay>19)&&', ROOT.kViolet),
        # ('fromtruetau3prong', 'Sec. from W(#tau #rightarrow 3 pr.)', '(maxscore_' + nn_name + '_Signal_deltam0p3_genMatchMotherIsTheTau==1)&&(maxscore_' + nn_name + '_Signal_deltam0p3_genMatchMotherTauDecay>9)&&(maxscore_' + nn_name + '_Signal_deltam0p3_genMatchMotherTauDecay<20)&&', ROOT.kAzure-6),
        # ('fromtruetau1prong', 'Sec. from W(#tau #rightarrow 1 pr.)', '(maxscore_' + nn_name + '_Signal_deltam0p3_genMatchMotherIsTheTau==1)&&(maxscore_' + nn_name + '_Signal_deltam0p3_genMatchMotherTauDecay<10)&&', ROOT.kAzure+6),

        # TODO: use this instead?
        ('notau', 'PV-associated', '(maxscore_' + nn_name + '_Signal_deltam0p3_hasGenMatch>0)&&(maxscore_' + nn_name + '_Signal_deltam0p3_genMatchMotherIsTheTau==0)&&', colors[1]),

        # ('prompt', 'Prompt', '(maxscore_' + nn_name + '_Signal_deltam0p3_hasGenMatch>0)&&(maxscore_' + nn_name + '_Signal_deltam0p3_genMatchIsPrompt>0)&&(maxscore_' + nn_name + '_Signal_deltam0p3_genMatchMotherIsTheTau==0)&&', ROOT.kGray),
        # ('secondary', 'Secondary', '(maxscore_' + nn_name + '_Signal_deltam0p3_hasGenMatch>0)&&(maxscore_' + nn_name + '_Signal_deltam0p3_genMatchIsPrompt==0)&&(maxscore_' + nn_name + '_Signal_deltam0p3_genMatchMotherIsTheTau==0)&&', ROOT.kGray+2),

        ('nogenmatch', 'Spurious', '(maxscore_' + nn_name + '_Signal_deltam0p3_hasGenMatch==0)&&', colors[3]),
    ]
    if len(args) > 0:
        eventclasses = [ec for ec in eventclasses if ec[0] == args['eventclass']]


    # ####
    # track-level
    # ####

    trackqualityselection = 'track_quality>2&&abs(track_eta)<2.4'
    trackdz1selection = trackqualityselection + '&&track_dz<1'
    tracksanityselection = trackqualityselection + '&&track_pt<50'
    trackfinalselection = trackqualityselection + '&&track_pt<20&&track_isPfCand>0'

    trackmaxscore0p3selection = trackqualityselection + '&&track_' + nn_name + '_Signal_deltam0p3_isMaxscore==1'
    trackmaxscore0p6selection = trackqualityselection + '&&track_' + nn_name + '_Signal_deltam0p6_isMaxscore==1'
    trackmaxscore1p0selection = trackqualityselection + '&&track_' + nn_name + '_Signal_deltam1p0_isMaxscore==1'

    # TODO: track-level?
    tracklevel = True
    # TODO: track-level selection
    trackselection = trackfinalselection + '&&track_' + nn_name + '_Signal_deltam0p3>(1/(1+exp(5)))' + '&&track_' + nn_name + '_Background_fromtruetau_deltam0p3>(1/(1+exp(-2)))'
    trackselection_refined = trackselection.replace('track_pt', 'track_pt' + refiner)

    # TODO: group tracks by classes?
    grouptracksbyclasses = True
    # TODO: uncomment track classes
    trackclasses = [
        ('fromtruetau', 'Sec. from W(#tau)', '(track_genMatchMotherIsTheTau==1)&&', colors[0], 1.),  # ROOT.kViolet

        # ('fromtruetaulep', 'Sec. from W(#tau #rightarrow e,mu)', '(track_genMatchMotherIsTheTau==1)&&(track_genMatchMotherTauDecay>19)&&', colors[4]),  # ROOT.kViolet
        # ('fromtruetau3prong', 'Sec. from W(#tau #rightarrow 3 pr.)', '(track_genMatchMotherIsTheTau==1)&&(track_genMatchMotherTauDecay>9)&&(track_genMatchMotherTauDecay<20)&&', colors[9]),  # ROOT.kAzure-6
        # ('fromtruetau1prong', 'Sec. from W(#tau #rightarrow 1 pr.)', '(track_genMatchMotherIsTheTau==1)&&(track_genMatchMotherTauDecay<10)&&', colors[0]),  # ROOT.kAzure+6

        ('notau', 'PV-associated', '(track_hasGenMatch>0)&&(track_genMatchMotherIsTheTau==0)&&', colors[1], 1.),  # ROOT.kGray

        # ('prompt', 'PV-associated', '(track_hasGenMatch>0)&&(track_genMatchIsPrompt>0)&&(track_genMatchMotherIsTheTau==0)&&', colors[1]),  # ROOT.kGray
        # ('secondary', '', '(track_hasGenMatch>0)&&((track_genMatchIsDirectHadronDecayProduct>0)||(track_genMatchIsDirectTauDecayProduct>0))&&(track_genMatchMotherIsTheTau==0)&&', colors[1]),  # ROOT.kGray+2

        ('nogenmatch', 'Spurious', '(track_hasGenMatch==0)&&', colors[3], 1.),  # ROOT.kGray+3
    ]
    # trackclasses = [
    #     # ('fromtruetau', 'Sec. from W(#tau)', '(track_genMatchMotherIsTheTau==1)&&', ROOT.kAzure-4),
    #
    #     ('fromtruetau1prong', 'Sec. from W(#tau #rightarrow 1 pr.)', '(track_genMatchMotherIsTheTau==1)&&(track_genMatchMotherTauDecay<10)&&', ROOT.kAzure+6),
    #     ('fromtruetau3prong', 'Sec. from W(#tau #rightarrow 3 pr.)', '(track_genMatchMotherIsTheTau==1)&&(track_genMatchMotherTauDecay>9)&&(track_genMatchMotherTauDecay<20)&&', ROOT.kAzure-6),
    #     ('fromtruetaulep', 'Sec. from W(#tau #rightarrow e,mu)', '(track_genMatchMotherIsTheTau==1)&&(track_genMatchMotherTauDecay>19)&&', ROOT.kViolet),
    #
    #     ('hard', 'Hard', '(track_genMatchMotherIsTheTau==0)&&(track_genMatchIsFromHardProcess==1)&&', ROOT.kGray),
    #     ('ue', 'UE', '(track_genMatchMotherIsTheTau==0)&&(track_genMatchIsFromHardProcess==0)&&', ROOT.kGray+2),
    #     ('nogenmatch', 'No GEN Match', '(track_genMatchIsFromHardProcess==-1)&&', ROOT.kGray+3),
    # ]
    if len(args) > 0:
        trackclasses = [tc for tc in trackclasses if tc[0] == args['eventclass']]


    # TODO: test?
    ntestfiles = 0
    ntestfilesdata = 0  # 1  # 100

    # TODO: tree/histo options
    savehistos = True  # better set to False if using histo cache
    usehistocache = False
    checkcounter = False  # only happening if not using histo cache

    # TODO: what to include?
    includeallsignals = False
    includesignal = False
    includeinsignalbkg = False  # only if track-level
    includecleaneddymc = False
    includecleaneddydata = False
    includedata = False
    if len(args) > 0:
        includeallsignals = args['dataset'] == 'allsignals'
        includesignal = args['dataset'] == 'signal'
        includeinsignalbkg = args['dataset'] == 'insignalbkg'  # only if track-level
        includecleaneddymc = args['dataset'] == 'cleaneddymc' or (args['dataset'].startswith('DYJetsToLL_M-50_HT-') and args['dataset'].endswith('_cleaned'))
        includecleaneddydata = args['dataset'] == 'cleaneddydata'
        includedata = args['dataset'] == 'data'

    # standardweight = str(lumi) + '*weight_PU_data:=[TFile f("/data/dust/user/wolfmor/NTupleStuff/PUweights_data16EF.root"); TH1D *h = (TH1D*)f.Get("n_pvdata:STACK"); Int_t bin = h->FindBin(n_pv); return h->GetBinContent(bin);]' \
    #                              '*XSEC/NSIM'
    standardweight = 'LUMI*weight_PU_MCData*crossSection/COUNTER'

    cleaneddyzptweight = {
        'era16_UL_APV': 'weight_zpt:=[TFile f("/data/dust/user/wolfmor/NTupleStuff/CleanedDY_ZptWeight_era16APV_UL/histos.root"); TH1D *h = (TH1D*)f.Get("zGamma_ptratioSTACKVScleanedDYmc"); Int_t bin = h->FindBin(zGamma_pt); return h->GetBinContent(bin);]',
        'era16_UL': 'weight_zpt:=[TFile f("/data/dust/user/wolfmor/NTupleStuff/CleanedDY_ZptWeight_era16_UL/histos.root"); TH1D *h = (TH1D*)f.Get("zGamma_ptratioSTACKVScleanedDYmc"); Int_t bin = h->FindBin(zGamma_pt); return h->GetBinContent(bin);]',
        'era17_UL': 'weight_zpt:=[TFile f("/data/dust/user/wolfmor/NTupleStuff/CleanedDY_ZptWeight_era17_UL/histos.root"); TH1D *h = (TH1D*)f.Get("zGamma_ptratioSTACKVScleanedDYmc"); Int_t bin = h->FindBin(zGamma_pt); return h->GetBinContent(bin);]',
        'era18_UL': 'weight_zpt:=[TFile f("/data/dust/user/wolfmor/NTupleStuff/CleanedDY_ZptWeight_era18_UL/histos.root"); TH1D *h = (TH1D*)f.Get("zGamma_ptratioSTACKVScleanedDYmc"); Int_t bin = h->FindBin(zGamma_pt); return h->GetBinContent(bin);]',
    }

    # mcdzerrorweight_tracklevel = 'weight_dzerrorTRACKLEVEL:=[get_mcdzerrorweight_vector(track_log10_dzError_);]'
    # mcdzerrorweight_tracklevel_cpp_code = '''
    # auto get_mcdzerrorweight_vector = [](const ROOT::VecOps::RVec<float>& track_log10_dzError_) {
    #     TFile f("/data/dust/user/wolfmor/NTupleStuff/MCdzErrorWeight/histos.root");
    #     TH1D *h = (TH1D*)f.Get("track_log10_dzError_ratiocleanedDYdataVSSTACK");
    #
    #     const auto size = track_log10_dzError_.size();
    #     ROOT::VecOps::RVec<double> weight_vec(size);
    #
    #     for (size_t i = 0; i < size; ++i) {
    #         Float_t theDzError = track_log10_dzError_[i];
    #         Int_t bin = h->FindBin(theDzError);
    #         weight_vec[i] = h->GetBinContent(bin);
    #     }
    #
    #     return weight_vec;
    # };
    # '''
    # ROOT.gInterpreter.Declare(mcdzerrorweight_tracklevel_cpp_code)

    signalweight = standardweight.replace('COUNTER', 'numSimEvents')
    if eras == ['era16_UL_APV']:
        signalweight = signalweight.replace('weight_PU_MCData', 'weight_PU_SignalData16pre')
    elif eras == ['era16_UL']:
        signalweight = signalweight.replace('weight_PU_MCData', 'weight_PU_SignalData16post')
    elif eras == ['era16_UL_APV', 'era16_UL'] or eras == ['era16_UL', 'era16_UL_APV']:
        signalweight = signalweight.replace('weight_PU_MCData', 'weight_PU_SignalData16full')

    # standardweight += '*' + mcdzerrorweight_tracklevel  # TODO: do this?

    # TODO: don't replace refiner name (but otherwise file name too long)
    outputpath = '/afs/desy.de/user/w/wolfmor/Plots/SoftDisplacedTrack/' \
                 + (tag if friends is None else 'FrieNdTuples_' + friends) + '/' \
                 + selection\
                     .replace(diagonalcut0p3vs0p6and1p0selection.replace(nobtagselection300, ''), '&&diagonalCut0p3vs0p6and1p0_' + nn_name.replace(refiner, '')) \
                     .replace(diagonalcut0p6vs0p3and1p0selection.replace(nobtagselection300, ''), '&&diagonalCut0p6vs0p3and1p0_' + nn_name.replace(refiner, '')) \
                     .replace(diagonalcut1p0vs0p3and0p6selection.replace(nobtagselection300, ''), '&&diagonalCut1p0vs0p3and0p6_' + nn_name.replace(refiner, '')) \
                     .replace('&&maxscore_' + nn_name + '_Signal_deltam1p0<maxscore_' + nn_name + '_Signal_deltam0p3', '&&diagonalCut0p3vs1p0_' + nn_name) \
                     .replace('&&maxscore_' + nn_name + '_Signal_deltam1p0<maxscore_' + nn_name + '_Signal_deltam0p6', '&&diagonalCut0p6vs1p0_' + nn_name) \
                     .replace('&&maxscore_' + nn_name + '_Signal_deltam0p3<maxscore_' + nn_name + '_Signal_deltam1p0', '&&diagonalCut1p0vs0p3_' + nn_name) \
                     .replace('&&maxscore_' + nn_name + '_Signal_deltam0p6<maxscore_' + nn_name + '_Signal_deltam1p0', '&&diagonalCut1p0vs0p6_' + nn_name) \
                 # + ''
    if requiremctrigger:
        outputpath += '_withMCtrigger'
    if tracklevel:
        outputpath += '/tracks/' \
                      + trackselection\
                          .replace('&&track_' + nn_name + '_Signal_deltam0p3_isMaxscore==1', '&&isMaxscore0p3_' + nn_name) \
                          .replace('&&track_' + nn_name + '_Signal_deltam0p6_isMaxscore==1', '&&isMaxscore0p6_' + nn_name) \
                          .replace('&&track_' + nn_name + '_Signal_deltam1p0_isMaxscore==1', '&&isMaxscore1p0_' + nn_name) \
                      # + '_trackclasses'  # '_alternativetrackclasses'

    outputsubfolder = 'era' + '-'.join(eras).replace('era', '').replace('UL', '').replace('_', '')
    # TODO: output subfolder?
    outputsubfolder += '/' + args['outputsubsubfolder'] + refiner if len(args) > 0 else ''  # '/SRbinning_allsignal'  # '/partCleanedDYMC_zptweighted'  # '/zptWeight'  # '/partData'  # '_signalPurity1p0_deltam0p3'  # '_cleaningTestOldVsNew'

    pf = PlotFactory(
        outputpath=outputpath,
        outputsubfolder=outputsubfolder,

        # TODO: output pattern/format?
        outputpattern='VARIABLE',  # _onlyClDYunscaled_kfactor_zptweighted _wClDY_zptweighted_onlyZinv _SRbinning
        outputformat=[],  # 'pdf',  # ['pdf', 'png'],  #

        # TODO: histo options
        uoflowbins=True,
        normalize=False,  # don't normalize when dividing stack into multiple parts otherwise can't combine anymore!
        ylabel='Tracks' if tracklevel else 'Events',  # 'A.U.',  # 'Simulated Events',  # Fraction of

        # TODO: ratio panel options

        yaxisrangeratio=(0.0001, 1.99999),
        # yaxisrangeratio=(0.0001, 2.99999),

        # ylabelratio='Purity',
        # ylabelratio='#scale[0.75]{#splitline{Cut Sig.}{#color[622]{In-sig. B. / Stack}}}' if friends is not None else 'Data / MC',
        # ylabelratio='#scale[0.85]{#splitline{Cut Sig.}{#color[418]{Cl. DY / Stack}}}' if friends is not None else 'Data / MC',
        # ylabelratio='#scale[0.85]{#splitline{Bin Sig.}{#color[418]{Cl. DY / Stack}}}' if friends is not None else 'Data / MC',
        # ylabelratio='#scale[0.55]{#splitline{#splitline{Data / MC}{#color[632]{Cut Sig.}}}{#color[418]{Cl. DY / MC}}}' if friends is not None else 'Data / MC',
        # ylabelratio='#scale[0.85]{#splitline{In-sig. B / B}{Cl. DY / B}}',  # #color[418]{}
        # ylabelratio='#scale[0.85]{#splitline{#color[418]{Cl. DY / MC}}{#color[835]{Z(inv) / MC}}}',
        # ylabelratio='#scale[0.5]{#splitline{#color[418]{Cl. DY / Stack}}{#splitline{#color[835]{Z(inv) / Stack}}{Cut Sig.}}}',
        # ylabelratio='#scale[0.8]{Cl. DY / Z(inv)}',  # #color[418]{}
        # ylabelratio='#scale[0.8]{Cl. DY / Stack}',  # #color[418]{}
        # ylabelratio='#scale[0.8]{Z(inv) / Cl. DY}',
        ylabelratio='Data / MC',
        # ylabelratio='Cut Sig.',
        # ylabelratio='Bin Sig.',
        # ylabelratio='#scale[0.85]{#splitline{#color[632]{Bin Sig.}}{Data / MC}}',
        # ylabelratio='In-sig. B / B',
        # ylabelratio='S/(S+B)',

        # ratiohlines=None,
        # ratiohlines=[0.15],
        # ratiohlines=[1, 2],

        ncolumnslegend=2,

        text=str(round(lumi, 1)) + ' fb^{-1} (13 TeV)',
        extratext='Work in progress' if includedata or includecleaneddydata else '#splitline{Work in progress}{Simulation}',

        # for "private work" label
        # cmstext='Private work',
        # extratext='#scale[0.8]{(CMS data/simulation)}' if includedata or includecleaneddydata else '#scale[0.8]{(CMS simulation)}',

        # poslegend=(0.45, 0.55, 0.93, 0.9),  # use if no ratio panel
    )

    # TODO: what backgrounds to include?
    backgrounds = [

        # ('DY(lep)Jets', ROOT.kCyan+2, '', 1., [
        #     'DYJetsToLL_M-50_Zpt-200toInf',
        # ]),
        #
        # ('TTZ/TTW', ROOT.kOrange+1, '', 1., [
        #     'TTZToLLNuNu_M-10',
        #     'TTZToQQ',
        #     'TTWJetsToLNu',
        #     'TTWJetsToQQ',
        # ]),

        # TODO: use this for track-level categories (otherwise tau labels not shown)
        ('WJets', ROOT.kViolet+1, '', 1.21, [
            'WJetsToLNu_HT-100To200',  # <0.1% for nobtagselection
            'WJetsToLNu_HT-200To400',
            'WJetsToLNu_HT-400To600',
            'WJetsToLNu_HT-600To800',
            'WJetsToLNu_HT-800To1200',
            'WJetsToLNu_HT-1200To2500',
            'WJetsToLNu_HT-2500ToInf',
        ]),

        # ('QCD', ROOT.kGray+1, '', 1., [
        #     'QCD_HT100to200',
        #     'QCD_HT200to300',
        #     'QCD_HT300to500',
        #     'QCD_HT500to700',
        #     'QCD_HT700to1000',
        #     'QCD_HT1000to1500',
        #     'QCD_HT1500to2000',
        #     'QCD_HT2000toInf',
        # ]),

        ('Diboson', ROOT.kCyan, '', 1., [
            'WW',
            'WZ',
            'ZZ',
        ]),

        ('SingleTop', ROOT.kOrange+2, '', 1., [
            'ST_t-channel_antitop',
            'ST_t-channel_top',
            'ST_tW_antitop',
            'ST_tW_top',
        ]),

        ('TT(lep)Jets', ROOT.kOrange+3, '', 1., [
            'TTJets_DiLept',
            'TTJets_SingleLeptFromT',
            'TTJets_SingleLeptFromTbar',
        ]),

        # TODO: use this for event-level plots (when NOT using eventclasses)
        # ('W(e,mu)Jets', colors[4], 'wBoson_tauDecayMode<0&&', 1.21, [  # ROOT.kViolet+1
        #     # 'WJetsToLNu_HT-100To200',  # <0.1% for nobtagselection
        #     'WJetsToLNu_HT-200To400',
        #     'WJetsToLNu_HT-400To600',
        #     'WJetsToLNu_HT-600To800',
        #     'WJetsToLNu_HT-800To1200',
        #     'WJetsToLNu_HT-1200To2500',
        #     'WJetsToLNu_HT-2500ToInf',
        # ]),
        # ('W(tau)Jets', colors[0], 'wBoson_tauDecayMode>-1&&', 1.21, [  # ROOT.kAzure-4
        #     # 'WJetsToLNu_HT-100To200',  # <0.1% for nobtagselection
        #     'WJetsToLNu_HT-200To400',
        #     'WJetsToLNu_HT-400To600',
        #     'WJetsToLNu_HT-600To800',
        #     'WJetsToLNu_HT-800To1200',
        #     'WJetsToLNu_HT-1200To2500',
        #     'WJetsToLNu_HT-2500ToInf',
        # ]),

        # TODO: zinv?
        ('Z(inv)Jets', colors[2], '', 1.23, [  # ROOT.kTeal-5
            'ZJetsToNuNu_HT-100To200',  # <0.1% for nobtagselection
            'ZJetsToNuNu_HT-200To400',
            'ZJetsToNuNu_HT-400To600',
            'ZJetsToNuNu_HT-600To800',
            'ZJetsToNuNu_HT-800To1200',
            'ZJetsToNuNu_HT-1200To2500',
            'ZJetsToNuNu_HT-2500ToInf',
        ]),
    ]

    if tracklevel and grouptracksbyclasses:

        for icls, cls in enumerate(trackclasses):
            for ibkg, bkg in enumerate(backgrounds):

                if 'fromtruetau' in cls[0] and not bkg[0] in ['WJets', 'W(tau)Jets']: continue  # TODO: adapt?

                if 'COUNTER' in standardweight:
                    pf.add_samples([
                        TreeSample(ntestfiles=ntestfiles, category='stack', name=cls[0]+era+b+bkg[2], title=cls[1] if (iera == 0 and ib == 0 and ibkg == 0) else '',
                                   tree='tEvent', files=datasets[era][b],
                                   friendfileslambda=None if friends is None else lambda x: x.replace('/NTuples/', '/FrieNdTuples/').replace('/' + tag + '/', '/' + friends + '/').replace('.root', '_friend.root'),
                                   eventselection=bkg[2]+selection+('&&triggerfired_met>0' if requiremctrigger else ''),
                                   weight=str(bkg[3])+'*'+standardweight.replace('LUMI', str(lumis[era])), checkcounter=checkcounter,
                                   vectorselection=cls[2]+trackselection_refined if tracklevel else None,
                                   skipemptyfiles=emptyntupletxts.replace('ERA', era).replace('DATASET', b),
                                   color=cls[3], usehistocache=usehistocache,
                                   usedataframeandweightfrom=False if icls == 0 else trackclasses[0][0]+era+b+bkg[2])
                        for iera, era in enumerate(eras)
                        for ib, b in enumerate(bkg[4]) if len(args) == 0 or b == args['dataset']
                    ])
                else:
                    pf.add_samples([
                        TreeSample(ntestfiles=ntestfiles, category='stack', name=cls[0]+era+bkg[0].replace(' ', '').replace('/', ''), title=cls[1] if (iera == 0 and ibkg == 0) else '',
                                   tree='tEvent', files=[datasets[era][b] for ib, b in enumerate(bkg[4]) if len(args) == 0 or b == args['dataset']],
                                   friendfileslambda=None if friends is None else lambda x: x.replace('/NTuples/', '/FrieNdTuples/').replace('/' + tag + '/', '/' + friends + '/').replace('.root', '_friend.root'),
                                   eventselection=bkg[2]+selection+('&&triggerfired_met>0' if requiremctrigger else ''),
                                   weight=str(bkg[3])+'*'+standardweight.replace('LUMI', str(lumis[era])), checkcounter=checkcounter,
                                   vectorselection=cls[2]+trackselection_refined if tracklevel else None,
                                   color=cls[3], usehistocache=usehistocache)
                        for iera, era in enumerate(eras)
                    ])

    elif groupeventsbyclasses:

        for icls, cls in enumerate(eventclasses):
            for ibkg, bkg in enumerate(backgrounds):

                if 'fromtruetau' in cls[0] and not bkg[0] in ['WJets', 'W(tau)Jets']: continue  # TODO: adapt?

                if 'COUNTER' in standardweight:
                    pf.add_samples([
                        TreeSample(ntestfiles=ntestfiles, category='stack', name=cls[0]+era+b+bkg[2], title=cls[1] if (iera == 0 and ib == 0 and ibkg == 0) else '',
                                   tree='tEvent', files=datasets[era][b],
                                   friendfileslambda=None if friends is None else lambda x: x.replace('/NTuples/', '/FrieNdTuples/').replace('/' + tag + '/', '/' + friends + '/').replace('.root', '_friend.root'),
                                   eventselection=cls[2]+bkg[2]+selection+('&&triggerfired_met>0' if requiremctrigger else ''),
                                   weight=str(bkg[3])+'*'+standardweight.replace('LUMI', str(lumis[era])), checkcounter=checkcounter,
                                   skipemptyfiles=emptyntupletxts.replace('ERA', era).replace('DATASET', b),
                                   color=cls[3], usehistocache=usehistocache)
                        for iera, era in enumerate(eras)
                        for ib, b in enumerate(bkg[4]) if len(args) == 0 or b == args['dataset']
                    ])
                else:
                    pf.add_samples([
                        TreeSample(ntestfiles=ntestfiles, category='stack', name=cls[0]+era+bkg[0].replace(' ', '').replace('/', ''), title=cls[1] if (iera == 0 and ibkg == 0) else '',
                                   tree='tEvent', files=[datasets[era][b] for ib, b in enumerate(bkg[4]) if len(args) == 0 or b == args['dataset']],
                                   friendfileslambda=None if friends is None else lambda x: x.replace('/NTuples/', '/FrieNdTuples/').replace('/' + tag + '/', '/' + friends + '/').replace('.root', '_friend.root'),
                                   eventselection=cls[2]+bkg[2]+selection+('&&triggerfired_met>0' if requiremctrigger else ''),
                                   weight=str(bkg[3])+'*'+standardweight.replace('LUMI', str(lumis[era])), checkcounter=checkcounter,
                                   color=cls[3], usehistocache=usehistocache)
                        for iera, era in enumerate(eras)
                    ])

    else:

        for bkg in backgrounds:
            if 'COUNTER' in standardweight:
                pf.add_samples([
                    TreeSample(ntestfiles=ntestfiles, category='stack', name=era+b+bkg[2], title=bkg[0] if (iera == 0 and ib == 0) else '',
                               tree='tEvent', files=datasets[era][b],
                               friendfileslambda=None if friends is None else lambda x: x.replace('/NTuples/', '/FrieNdTuples/').replace('/' + tag + '/', '/' + friends + '/').replace('.root', '_friend.root'),
                               eventselection=bkg[2]+selection+('&&triggerfired_met>0' if requiremctrigger else ''),
                               weight=str(bkg[3])+'*'+standardweight.replace('LUMI', str(lumis[era])), checkcounter=checkcounter,
                               vectorselection=trackselection_refined if tracklevel else None,
                               skipemptyfiles=emptyntupletxts.replace('ERA', era).replace('DATASET', b),
                               color=bkg[1], usehistocache=usehistocache)
                    for iera, era in enumerate(eras)
                    for ib, b in enumerate(bkg[4]) if len(args) == 0 or b == args['dataset']
                ])
            else:
                pf.add_samples([
                    TreeSample(ntestfiles=ntestfiles, category='stack', name=era+bkg[0].replace(' ', '').replace('/', ''), title=bkg[0] if iera == 0 else '',
                               tree='tEvent', files=[datasets[era][b] for ib, b in enumerate(bkg[4]) if len(args) == 0 or b == args['dataset']],
                               friendfileslambda=None if friends is None else lambda x: x.replace('/NTuples/', '/FrieNdTuples/').replace('/' + tag + '/', '/' + friends + '/').replace('.root', '_friend.root'),
                               eventselection=bkg[2]+selection+('&&triggerfired_met>0' if requiremctrigger else ''),
                               weight=str(bkg[3])+'*'+standardweight.replace('LUMI', str(lumis[era])), checkcounter=checkcounter,
                               vectorselection=trackselection_refined if tracklevel else None,
                               color=bkg[1], usehistocache=usehistocache)
                    for iera, era in enumerate(eras)
                ])


    if includecleaneddymc:

        cldymc = ('Cleaned DY MC', ROOT.kGreen+2, '', 1.23, [
            'DYJetsToLL_M-50_HT-100to200_cleaned',  # <0.1% for nobtagselection
            'DYJetsToLL_M-50_HT-200to400_cleaned',
            'DYJetsToLL_M-50_HT-400to600_cleaned',
            'DYJetsToLL_M-50_HT-600to800_cleaned',
            'DYJetsToLL_M-50_HT-800to1200_cleaned',
            'DYJetsToLL_M-50_HT-1200to2500_cleaned',
            'DYJetsToLL_M-50_HT-2500toInf_cleaned',
        ])

        if tracklevel and grouptracksbyclasses:
            for icls, cls in enumerate(trackclasses):

                if 'fromtruetau' in cls[0]: continue

                # TODO: tighter Z mass window
                pf.add_samples([
                    TreeSample(ntestfiles=ntestfiles, category='stack', name=cls[0]+era+b+cldymc[2],  # TODO: group='cleanedDYmc',
                               title=cldymc[0] + ' (' + cls[1] + ')' if (iera == 0 and ib == 0) else '',
                               tree='tEvent', files=datasets[era][b],
                               friendfileslambda=None if friends is None else lambda x: x.replace('/NTuples/', '/FrieNdTuples/').replace('/' + tag + '/', '/' + friends + '/').replace('.root', '_friend.root'),  # .replace('V15p0', 'V15p14')
                               eventselection='(cleaning_muonsCleaned==1)&&(cleaning_l1dBetaRelIso<0.2)&&(cleaning_l2dBetaRelIso<0.2)&&(cleaning_zPt>200)&&'+cldymc[2]+selection,
                               weight=str(cldymc[3])+'*'+standardweight.replace('LUMI', str(lumis[era])) + '*' + cleaneddyzptweight[era],   # + '*' + cleaneddymetweight,   # TODO: cleaneddy...weight?,
                               checkcounter=checkcounter,
                               vectorselection=cls[2]+trackselection_refined if tracklevel else None,
                               skipemptyfiles=emptyntupletxts.replace('ERA', era).replace('DATASET', b),
                               # TODO: scale cleanedDYmc?
                               # TODO: don't scale to full range?
                               scaleto=None,  # 'STACK:200:END',  # 'cleanedDYdata:START:END',  # 'STACK:START:END',  #
                               # scaleby=cls[4],  # TODO: MC correction factors?
                               color=cls[3],  # TODO: cldymc[1],
                               linestyle=ROOT.kSolid,
                               usehistocache=usehistocache)
                    for iera, era in enumerate(eras)
                    for ib, b in enumerate(cldymc[4]) if len(args) == 0 or b == args['dataset']
                ])

        else:

            # TODO: tighter Z mass window
            pf.add_samples([
                TreeSample(ntestfiles=ntestfiles, category='line', name=era+b+cldymc[2], group='cleanedDYmc', title=cldymc[0] if (iera == 0 and ib == 0) else '',
                           tree='tEvent', files=datasets[era][b],
                           friendfileslambda=None if friends is None else lambda x: x.replace('/NTuples/', '/FrieNdTuples/').replace('/' + tag + '/', '/' + friends + '/').replace('.root', '_friend.root'),
                           eventselection='(cleaning_muonsCleaned==1)&&(cleaning_l1dBetaRelIso<0.2)&&(cleaning_l2dBetaRelIso<0.2)&&(cleaning_zPt>200)&&'+cldymc[2]+selection,
                           weight=str(cldymc[3])+'*'+standardweight.replace('LUMI', str(lumis[era])) + '*' + cleaneddyzptweight[era],   # + '*' + cleaneddymetweight,   # TODO: cleaneddy...weight?,
                           checkcounter=checkcounter,
                           vectorselection=trackselection_refined if tracklevel else None,
                           skipemptyfiles=emptyntupletxts.replace('ERA', era).replace('DATASET', b),
                           # TODO: scale cleanedDYmc?
                           # TODO: don't scale to full range?
                           scaleto=None,  # 'STACK:200:END',  # 'cleanedDYdata:START:END',  # 'STACK:START:END',  #
                           # scaleby=1.23,  # TODO: add k-factor of 1.23?
                           color=cldymc[1], linestyle=ROOT.kSolid,
                           usehistocache=usehistocache)
                for iera, era in enumerate(eras)
                for ib, b in enumerate(cldymc[4]) if len(args) == 0 or b == args['dataset']
            ])

        # # TODO: tighter Z mass window
        # pf.add_samples([
        #     TreeSample(ntestfiles=ntestfiles, category='line', name=era+'cleanedDYmc', group='cleanedDYmc', title='Cleaned DY MC',
        #                tree='tEvent', files=[datasets[era][d] for d in datasets[era] if 'DYJetsToLL_M-50_HT' in d and d.endswith('_cleaned')],
        #                friendfileslambda=None if friends is None else lambda x: x.replace('/NTuples/', '/FrieNdTuples/').replace('/' + tag + '/', '/' + friends + '/').replace('.root', '_friend.root'),
        #                eventselection='(cleaning_muonsCleaned==1)&&(cleaning_l1dBetaRelIso<0.2)&&(cleaning_l2dBetaRelIso<0.2)&&(cleaning_zPt>200)&&'+selection,
        #                weight=standardweight.replace('COUNTER', 'numSimEvents').replace('LUMI', str(lumis[era])) + '*' + cleaneddyzptweight[era],   # + '*' + cleaneddymetweight,   # TODO: cleaneddy...weight?
        #                checkcounter=checkcounter,
        #                vectorselection=trackselection_refined if tracklevel else None,
        #                skipemptyfiles=emptyntupletxts.replace('ERA', era).replace('DATASET', b),
        #                # TODO: scale cleanedDYmc?
        #                # TODO: don't scale to full range?
        #                scaleto=None,  # 'STACK:200:END',  # 'cleanedDYdata:START:END',  # 'STACK:START:END',  #
        #                # scaleby=1.23,  # TODO: add k-factor of 1.23?
        #                color=ROOT.kGreen+2, linestyle=ROOT.kSolid,
        #                usehistocache=usehistocache,
        #                # histocachefile=outputpath.replace('ht5/ht', 'ht5overht') + '/' + {'era16_UL_APV': 'era16APV', 'era16_UL': 'era16', 'era17_UL': 'era17', 'era18_UL': 'era18'}[era] + '/zptweighted_partCleanedDYMC/histos.root'  # TODO: specify histo cache file?
        #                )
        #     for iera, era in enumerate(eras)
        # ])

        # TODO: add cleanedDYmc ratio?
        # pf.add_ratios([
        #     Ratio(category='ratio', name='cleanedDYmc:STACK', drawoptions=['hist', 'e2:=3002']),
        #     # Ratio(category='ratio', name='STACK:cleanedDYmc', drawoptions=['hist', 'e2:=3002'], savehistos=True, variablestosave=['zGamma_pt']),
        #     # Ratio(category='ratio', name='cleanedDYmc:cleanedDYdata', drawoptions=['hist', 'e2:=3002']),
        # ])


    if tracklevel and includeinsignalbkg:
        pf.add_samples([
            TreeSample(ntestfiles=ntestfiles, category='line', name=era+'insignal_bkg', group='insignal_bkg', title='In-signal Bkg.',
                       tree='tEvent', files=[datasets[era][d].replace('_part*of25_', '_part1of25_') for d in datasets[era] if 'SignalV4' in d],
                       friendfileslambda=None if friends is None else lambda x: x.replace('/NTuples/', '/FrieNdTuples/').replace('/' + tag + '/', '/' + friends + '/').replace('.root', '_friend.root'),
                       eventselection=selection, weight=signalweight.replace('LUMI', str(lumis[era])), checkcounter=checkcounter,
                       vectorselection='track_isSusyTrack==0&&'+trackselection if tracklevel else None,
                       # TODO: scale in-signal background?
                       scaleto=None,  # 'STACK:START:END',  #
                       color=ROOT.kBlack, linestyle=ROOT.kDashed,  # ROOT.kRed-10
                       usehistocache=usehistocache)
            for iera, era in enumerate(eras)
        ])

        # TODO: add in-signal background ratio?
        # pf.add_ratios([
        #     Ratio(category='ratio', name='insignal_bkg:STACK'),
        # ])


    if includeallsignals:

        if len(refiner) > 0: continue

        signalmodelpoints = {}
        for path in [datasets[era][d] for era in eras for d in datasets[era] if 'SignalV4' in d]:
            for f in glob(path):

                mChi = re.search(r'mChipm(.*?)GeV', f).group(1)
                dm = re.search(r'dm(.*?)GeV', f).group(1)

                signalmodelpoint = 'mChi' + mChi + '_dm' + dm

                if signalmodelpoint not in signalmodelpoints.keys():
                    signalmodelpoints[signalmodelpoint] = [f]
                else:
                    signalmodelpoints[signalmodelpoint].append(f)

        for signalmodelpoint in signalmodelpoints:

            if 'mChi225' in signalmodelpoint or 'mChi400' in signalmodelpoint: continue  # not present in 17/18

            pf.add_samples([
                TreeSample(ntestfiles=ntestfiles, category='line', name=era+'signal_direct_' + signalmodelpoint, group='signal_direct_' + signalmodelpoint, title='#scale[0.8]{(#chi,#Deltam^{#pm})=(' + signalmodelpoint.split('_')[0].replace('mChi', '') + ',' + signalmodelpoint.split('_')[1].replace('dm', '').replace('p', '.') + ')}',
                           tree='tEvent', files=signalmodelpoints[signalmodelpoint],
                           friendfileslambda=None if friends is None else lambda x: x.replace('/NTuples/', '/FrieNdTuples/').replace('/' + tag + '/', '/' + friends + '/').replace('.root', '_friend.root'),
                           eventselection=selection,  weight=signalweight.replace('LUMI', str(lumis[era])),
                           vectorselection='track_isSignalTrack==1&&'+trackselection if tracklevel else None,
                           color=ROOT.kRed, linestyle=ROOT.kSolid,
                           usehistocache=usehistocache, checkcounter=checkcounter)
                for iera, era in enumerate(eras)
            ])

            # TODO: add significance scans?
            # if friends is not None:
            #     pf.add_ratios([
            #         # Ratio(category='cutsig', name='signal_direct_' + signalmodelpoint + ':STACK:0.2'),
            #         Ratio(category='binsig', name='signal_direct_' + signalmodelpoint + ':STACK:0.2'),
            #     ])

    elif includesignal:

        if len(refiner) > 0: continue

        pf.add_samples([
            TreeSample(ntestfiles=ntestfiles, category='line', name=era+'signal_direct_mChi115_dm0p27', group='signal_direct_mChi115_dm0p27', title='#scale[1.0]{(#chi,#Deltam^{#pm})=(115,0.3)}',
                       tree='tEvent', files=[datasets[era][d].replace('mChipm*GeV', 'mChipm115GeV').replace('dm*GeV', 'dm0p268GeV') for d in datasets[era] if 'SignalV4' in d],
                       friendfileslambda=None if friends is None else lambda x: x.replace('/NTuples/', '/FrieNdTuples/').replace('/' + tag + '/', '/' + friends + '/').replace('.root', '_friend.root'),
                       eventselection=selection,  weight=signalweight.replace('LUMI', str(lumis[era])),
                       vectorselection='track_isSignalTrack==1&&'+trackselection if tracklevel else None,
                       color=colors[7],  # ROOT.kOrange,
                       linestyle=ROOT.kSolid,
                       usehistocache=usehistocache, checkcounter=checkcounter)
            for iera, era in enumerate(eras)
        ])

        pf.add_samples([
            TreeSample(ntestfiles=ntestfiles, category='line', name=era+'signal_direct_mChi115_dm0p57', group='signal_direct_mChi115_dm0p57', title='#scale[1.0]{(#chi,#Deltam^{#pm})=(115,0.6)}',
                       tree='tEvent', files=[datasets[era][d].replace('mChipm*GeV', 'mChipm115GeV').replace('dm*GeV', 'dm0p568GeV') for d in datasets[era] if 'SignalV4' in d],
                       friendfileslambda=None if friends is None else lambda x: x.replace('/NTuples/', '/FrieNdTuples/').replace('/' + tag + '/', '/' + friends + '/').replace('.root', '_friend.root'),
                       eventselection=selection,  weight=signalweight.replace('LUMI', str(lumis[era])),
                       vectorselection='track_isSignalTrack==1&&'+trackselection if tracklevel else None,
                       color=colors[6],  # ROOT.kMagenta,
                       linestyle=ROOT.kSolid,
                       usehistocache=usehistocache, checkcounter=checkcounter)
            for iera, era in enumerate(eras)
        ])

        pf.add_samples([
            TreeSample(ntestfiles=ntestfiles, category='line', name=era+'signal_direct_mChi115_dm0p97', group='signal_direct_mChi115_dm0p97', title='#scale[1.0]{(#chi,#Deltam^{#pm})=(115,1.0)}',
                       tree='tEvent', files=[datasets[era][d].replace('mChipm*GeV', 'mChipm115GeV').replace('dm*GeV', 'dm0p968GeV') for d in datasets[era] if 'SignalV4' in d],
                       friendfileslambda=None if friends is None else lambda x: x.replace('/NTuples/', '/FrieNdTuples/').replace('/' + tag + '/', '/' + friends + '/').replace('.root', '_friend.root'),
                       eventselection=selection,  weight=signalweight.replace('LUMI', str(lumis[era])),
                       vectorselection='track_isSignalTrack==1&&'+trackselection if tracklevel else None,
                       color=colors[5],  # ROOT.kRed,
                       linestyle=ROOT.kSolid,
                       usehistocache=usehistocache, checkcounter=checkcounter)
            for iera, era in enumerate(eras)
        ])

        pf.add_samples([
            TreeSample(ntestfiles=ntestfiles, category='line', name=era+'signal_direct_mChi160_dm0p59', group='signal_direct_mChi160_dm0p59', title='#scale[1.0]{(#chi,#Deltam^{#pm})=(160,0.6)}',
                       tree='tEvent', files=[datasets[era][d].replace('mChipm*GeV', 'mChipm160GeV').replace('dm*GeV', 'dm0p587GeV') for d in datasets[era] if 'SignalV4' in d],
                       friendfileslambda=None if friends is None else lambda x: x.replace('/NTuples/', '/FrieNdTuples/').replace('/' + tag + '/', '/' + friends + '/').replace('.root', '_friend.root'),
                       eventselection=selection,  weight=signalweight.replace('LUMI', str(lumis[era])),
                       vectorselection='track_isSignalTrack==1&&'+trackselection if tracklevel else None,
                       color=colors[6],  # ROOT.kMagenta+2,
                       linestyle=ROOT.kDotted,  # ROOT.kSolid,
                       usehistocache=usehistocache, checkcounter=checkcounter)
            for iera, era in enumerate(eras)
        ])

        # TODO: add significance scans?
        # if friends is not None:
        #     pf.add_ratios([
        #         Ratio(category='binsig', name='signal_direct_mChi115_dm0p97:STACK:0.2', drawoptions=['hist', 'e2:=3002']),
        #         # Ratio(category='cutsig', name='signal_direct_mChi115_dm0p27:STACK:0.2', drawoptions=['hist', 'e2:=3002']),
        #         # Ratio(category='cutsig', name='signal_direct_mChi115_dm0p57:STACK:0.2', drawoptions=['hist', 'e2:=3002']),
        #         # Ratio(category='cutsig', name='signal_direct_mChi115_dm0p97:STACK:0.2', drawoptions=['hist', 'e2:=3002']),
        #         # Ratio(category='cutsig', name='signal_direct_mChi160_dm0p59:STACK:0.2', drawoptions=['hist', 'e2:=3002']),
        #     ])

        # for purity plots
        if False:
            pf.add_samples([

                # TreeSample(ntestfiles=ntestfiles, category='line', name='signal_direct_mChi115_dm0p27_susy', title='#scale[0.8]{SUSY matched}',
                #            tree='tEvent', files=path_signal + '/SignalV4_16/step3_higgsinoDm0eDmpm_RunIISpring21UL16FS_susyall_mChipm115GeV_dm0p268GeV_part*_NTuple_job*.root',
                #            friendfileslambda=None if friends is None else lambda x: x.replace('/NTuples/', '/FrieNdTuples/').replace('/' + tag + '/', '/' + friends + '/').replace('.root', '_friend.root'),
                #            eventselection=selection+'&&maxscore_' + nn_name + '_Signal_deltam1p0_isSusyTrack==1',
                #            color=ROOT.kOrange+1, linestyle=ROOT.kDotted,
                #            usehistocache=usehistocache),

                TreeSample(ntestfiles=ntestfiles, category='line', name='signal_direct_mChi115_dm0p97_susy', title='#scale[0.8]{SUSY matched}',
                           tree='tEvent', files=path_signal + '/SignalV4_16/step3_higgsinoDm0eDmpm_RunIISpring21UL16FS_susyall_mChipm115GeV_dm0p968GeV_part*_NTuple_job*.root',
                           friendfileslambda=None if friends is None else lambda x: x.replace('/NTuples/', '/FrieNdTuples/').replace('/' + tag + '/', '/' + friends + '/').replace('.root', '_friend.root'),
                           eventselection=selection+'&&maxscore_' + nn_name + '_Signal_deltam0p3_isSusyTrack==1',
                           color=ROOT.kRed+2, linestyle=ROOT.kDotted,
                           usehistocache=usehistocache),

                # TreeSample(ntestfiles=ntestfiles, category='line', name='signal_direct_mChi115_dm0p27_signal', title='#scale[0.8]{Signal pion matched}',
                #            tree='tEvent', files=path_signal + '/SignalV4_16/step3_higgsinoDm0eDmpm_RunIISpring21UL16FS_susyall_mChipm115GeV_dm0p268GeV_part*_NTuple_job*.root',
                #            friendfileslambda=None if friends is None else lambda x: x.replace('/NTuples/', '/FrieNdTuples/').replace('/' + tag + '/', '/' + friends + '/').replace('.root', '_friend.root'),
                #            eventselection=selection+'&&maxscore_' + nn_name + '_Signal_deltam1p0_isSignalTrack==1',
                #            color=ROOT.kOrange, linestyle=ROOT.kDashed,
                #            usehistocache=usehistocache),

                TreeSample(ntestfiles=ntestfiles, category='line', name='signal_direct_mChi115_dm0p97_signal', title='#scale[0.8]{Signal pion matched}',
                           tree='tEvent', files=path_signal + '/SignalV4_16/step3_higgsinoDm0eDmpm_RunIISpring21UL16FS_susyall_mChipm115GeV_dm0p968GeV_part*_NTuple_job*.root',
                           friendfileslambda=None if friends is None else lambda x: x.replace('/NTuples/', '/FrieNdTuples/').replace('/' + tag + '/', '/' + friends + '/').replace('.root', '_friend.root'),
                           eventselection=selection+'&&maxscore_' + nn_name + '_Signal_deltam0p3_isSignalTrack==1',
                           color=ROOT.kRed, linestyle=ROOT.kDashed,
                           usehistocache=usehistocache),
            ])

            pf.add_ratios([
                # Ratio(category='efficiency', name='signal_direct_mChi115_dm0p27_signal:signal_direct_mChi115_dm0p27'),
                Ratio(category='efficiency', name='signal_direct_mChi115_dm0p97_signal:signal_direct_mChi115_dm0p97'),
                # Ratio(category='efficiency', name='signal_direct_mChi115_dm0p27_susy:signal_direct_mChi115_dm0p27'),
                Ratio(category='efficiency', name='signal_direct_mChi115_dm0p97_susy:signal_direct_mChi115_dm0p97'),
            ])


    if includecleaneddydata:

        if len(refiner) > 0: continue

        # TODO: tighter Z mass window
        pf.add_samples([
            TreeSample(ntestfiles=ntestfilesdata, category='marker', name=era+'cleanedDYdata', group='cleanedDYdata', title='Cleaned DY Data',
                       tree='tEvent', files=[datasets[era][d] for d in datasets[era] if 'DataMuon' in d],
                       friendfileslambda=None if friends is None else lambda x: x.replace('/NTuples/', '/FrieNdTuples/').replace('/' + tag + '/', '/' + friends + '/').replace('.root', '_friend.root'),
                       eventselection='triggerfired_met>0&&(cleaning_muonsCleaned==1)&&(cleaning_l1dBetaRelIso<0.2)&&(cleaning_l2dBetaRelIso<0.2)&&(cleaning_zPt>200)&&'+selection,  # TODO: what trigger? 18.7.2024 change to triggerfired_met>0
                       # weight=cleaneddymetweight,  # TODO: cleaneddymetweight?
                       vectorselection=trackselection if tracklevel else None,
                       skipemptyfiles=emptyntupletxts.replace('ERA', era).replace('DATASET', 'DataMuon'),
                       # TODO: scale cleanedDYdata?
                       # TODO: don't scale to full range?
                       scaleto=None,  # 'cleanedDYmc:START:END',  # 'data:START:END',  # 'STACK:START:END',  #
                       color=ROOT.kGreen+2, usehistocache=usehistocache)
            for iera, era in enumerate(eras)
        ])

        # TODO: add cleanedDYdata ratio?
        # pf.add_ratios([
        #     Ratio(category='ratio', name='cleanedDYdata:cleanedDYmc'),
        #     # Ratio(category='ratio', name='cleanedDYdata:STACK'),
        #     # Ratio(category='ratio', name='cleanedDYdata:data'),
        # ])


    if includedata:

        if len(refiner) > 0: continue

        pf.add_samples([
            TreeSample(ntestfiles=ntestfilesdata, category='marker', name=era+'data', group='data', title='Data',
                       tree='tEvent', files=[datasets[era][d] for d in datasets[era] if 'DataMET' in d],
                       friendfileslambda=None if friends is None else lambda x: x.replace('/NTuples/', '/FrieNdTuples/').replace('/' + tag + '/', '/' + friends + '/').replace('.root', '_friend.root'),
                       eventselection='triggerfired_met>0&&'+selection,
                       vectorselection=trackselection if tracklevel else None,
                       skipemptyfiles=emptyntupletxts.replace('ERA', era).replace('DATASET', 'DataMET'),
                       color=ROOT.kBlack, usehistocache=usehistocache)
            for iera, era in enumerate(eras)
        ])

        # TODO: add data ratio?
        # pf.add_ratios([
        #     Ratio(category='ratio', name='data:STACK'),
        # ])


    if tracklevel:

        pf.add_variables([
            Variable.fromlist(variables[v], savehistos=savehistos) for v in variables if v.startswith('track_')
        ])

        if refiner in refiner_vars:
            pf.add_variables([
                Variable.fromlist(variables[v], savehistos=savehistos, append=refiner)
                for v in refiner_vars[refiner] + ['track_pt_widerange']
            ])

    else:

        # for GEN Z pT weight
        # import numpy as np
        # pf.add_variables([
        #     Variable.fromlist(['zGamma_pt', 'p_{T}(Z_{GEN})', (-25, 2500),
        #                        list(np.arange(-25, 500, 25)) + list(np.arange(500, 1000, 50)) + list(np.arange(1000, 2000, 100)) + list(np.arange(2000, 2500, 250)) + [2500],
        #                        101])
        # ])

        pf.add_variables([
            Variable.fromlist(variables[v], savehistos=savehistos) for v in variables if not v.startswith('track_') and not (v == 'deltam' or v == 'deltam_true')
        ])

        # if includedata:
        #     for v in pf.variables:
        #         if v.name in ['met_pt', 'htMiss', 'ht', 'ht5']:
        #             v.blind = ['data:500:END']


    if friends is not None:

        if tracklevel:
            nodes = [
                ('Signal', 'Signal', (-15, 15), 120),
                ('Background_nogenmatch', 'No GEN Match', (-15, 15), 120),
                ('Background_prompt', 'Prompt', (-15, 15), 120),
                ('Background_secondary', 'Secondary', (-15, 15), 120),
                ('Background_fromtruetau', 'W(#tau)', (-15, 15), 120)
            ]
        else:
            nodes = [
                ('Signal', 'Signal', (-5, 15), 80),
                ('Background_fromtruetau', 'W(#tau)', (-15, 15), 60)
            ]

        deltams = [
            ('deltam0p3', '0.3'),
            ('deltam0p6', '0.6'),
            ('deltam1p0', '1.0'),
        ]
        deltams_by_nn = {}
        for nn in nn_names:
            if nn in notparametrized:
                deltams_by_nn[nn] = [('', 'X')]
            else:
                deltams_by_nn[nn] = deltams

        if tracklevel:

            pf.add_variable(
                Variable('track_sanitycheck', title='This should be 0', vartoplot='track_random-track_random_sanitycheck', axisrange=(-10, 11), nbins=21, savehistos=savehistos)
            )

            for nn in nn_names:
                for name, title, axisrange, nbins in nodes:
                    for dmname, dmtitle in deltams_by_nn[nn]:

                        pf.add_variables([
                            Variable('track_' + nn + '_' + name + '_' + dmname + '_logit_raw', title='#tilde{P}(' + title + ' | #Deltam=' + dmtitle + ' GeV)', savehistos=savehistos,
                                     axisrange=axisrange, nbins=nbins,  # blind='data:4:END' if includedata and name == 'Signal' else None,  # nbins='flat:signal_direct_mChi115_dm0p27:0.2',  #
                                     vartoplot='log(track_' + nn + '_' + name + '_' + dmname + '/(1-track_' + nn + '_' + name + '_' + dmname + '))'),

                            # Variable('track_' + nn + '_' + name + '_' + dmname + '_logit_abs', title='#tilde{P}(' + title + ' | #Deltam=' + dmtitle + ' GeV)', savehistos=savehistos,
                            #          axisrange=axisrange, nbins=nbins,  blind='data:START:END' if includedata and name == 'Signal' else None,  # nbins='flat:signal_direct_mChi115_dm0p27:0.2',  #
                            #          vartoplot='log((abs(track_' + nn + '_' + name + '_' + dmname + ')+1e-6)/(1-abs(track_' + nn + '_' + name + '_' + dmname + ')+1e-6))'),

                            Variable('track_' + nn + '_' + name + '_' + dmname + '_logit', title='#tilde{P}(' + title + ' | #Deltam=' + dmtitle + ' GeV)', savehistos=savehistos,
                                     axisrange=axisrange, nbins=nbins,  # blind='data:4:END' if includedata and name == 'Signal' else None,  # nbins='flat:signal_direct_mChi115_dm0p27:0.2',  #
                                     vartoplot='log((track_' + nn + '_' + name + '_' + dmname + '+1e-6)/(1-track_' + nn + '_' + name + '_' + dmname + '+1e-6))'),

                            Variable('track_' + nn + '_' + name + '_' + dmname, title='P(' + title + ' | #Deltam=' + dmtitle + ' GeV)', savehistos=savehistos,
                                     axisrange=(0, 1), nbins=50),  # blind='data:0.9:END' if includedata and name == 'Signal' else None),
                        ])

                        if name == 'Signal':
                            pf.add_variables([
                                Variable('track_' + nn + '_' + name + '_' + dmname + '_ranking', title='In-event Ranking for P(' + title + ' | #Deltam=' + dmtitle + ' GeV)', savehistos=savehistos,
                                         axisrange=(0, 20), nbins=20),  # ,  blind='data:START:END' if includedata and name == 'Signal' else None),
                            ])

        else:

            pf.add_variable(
                Variable('sanitycheck', title='This should be 0', vartoplot='random-random_sanitycheck', axisrange=(-10, 11), nbins=21, savehistos=savehistos)
            )

            for nn in nn_names:
                for name, title, axisrange, nbins in nodes:
                    for dmname, dmtitle in deltams_by_nn[nn]:
                        pf.add_variables([
                            Variable('maxscore_' + nn + '_' + name + '_' + dmname + '_logit_raw', title='#tilde{P}_{max}(' + title + ' | #Deltam=' + dmtitle + ' GeV)', savehistos=savehistos,
                                     axisrange=axisrange, nbins=nbins,  # blind='data:4:END' if includedata else None,  # and name == 'Signal'
                                     vartoplot='log(maxscore_' + nn + '_' + name + '_' + dmname + '/(1-maxscore_' + nn + '_' + name + '_' + dmname + '))'),

                            Variable('maxscore_' + nn + '_' + name + '_' + dmname + '_logit', title='#tilde{P}_{max}(' + title + ' | #Deltam=' + dmtitle + ' GeV)', savehistos=savehistos,
                                     axisrange=(-5, 15), nbins=80,  # ,nbins='flat:signal_direct_mChi115_dm' + {'0.3': '0p268' if includeallsignals else '0p27', '0.6': '0p568' if includeallsignals else '0p57', '1.0': '0p968' if includeallsignals else '0p97'}[dmtitle] + ':111',  # TODO: flat S binning (adapt nbins)?
                                     # blind='data:4:END' if includedata else None,  # and name == 'Signal'
                                     vartoplot='log((maxscore_' + nn + '_' + name + '_' + dmname + '+1e-6)/(1-maxscore_' + nn + '_' + name + '_' + dmname + '+1e-6))', verbosesignificance=name == 'Signal'),

                            Variable('maxscore_' + nn + '_' + name + '_' + dmname + '_logit_fullrange', title='#tilde{P}_{max}(' + title + ' | #Deltam=' + dmtitle + ' GeV)', savehistos=savehistos,
                                     axisrange=(-15, 15), nbins=120,  # blind='data:4:END' if includedata else None,  # and name == 'Signal'
                                     vartoplot='log((maxscore_' + nn + '_' + name + '_' + dmname + '+1e-6)/(1-maxscore_' + nn + '_' + name + '_' + dmname + '+1e-6))'),

                            Variable('maxscore_' + nn + '_' + name + '_' + dmname, title='P_{max}(' + title + ' | #Deltam=' + dmtitle + ' GeV)', savehistos=savehistos,
                                     axisrange=(0, 1), nbins=50),  # ,  blind='data:0.9:END' if includedata else None),  # and name == 'Signal'

                            Variable('maxscore_' + nn + '_' + name + '_' + dmname + '_zoom', title='P_{max}(' + title + ' | #Deltam=' + dmtitle + ' GeV)', savehistos=savehistos,
                                     axisrange=(0.99, 1), nbins=50,  # blind='data:0.991:END' if includedata else None,  # and name == 'Signal'
                                     vartoplot='maxscore_' + nn + '_' + name + '_' + dmname),
                        ])

    # SR binning
    # 0.3: max cut sig for cut at 10.75 --> signal efficiency 0.3% --> 333 flat bins --> changed to 35 (tighest bin not at 10.75 but at 11.0)
    # 0.6: max cut sig for cut at 9.5 --> signal efficiency 1.6% --> 63 flat bins --> changed to 28
    # 1.0: max cut sig for cut at 8.0 --> signal efficiency 0.9% --> 111 flat bins --> changed to 24

    pf.process(dryrun=False)
