# list of variables: name, title, (vartoplot,) axisrange, rebin, nbins
myvariables = {

    # '': ['', '', (0, 5), 1, 50],

    'triggerfired': ['triggerfired', 'Trigger Fired', (0, 2), 1, 2],
    'weight_PU_MCData': ['weight_PU_MCData', 'weight_PU_MCData', (0, 10), 1, 50],
    'weight_PU_SignalMC': ['weight_PU_SignalMC', 'weight_PU_SignalMC', (0, 10), 1, 50],
    'weight_PU_SignalData': ['weight_PU_SignalData', 'weight_PU_SignalData', (0, 10), 1, 50],

    'n_pv': ['n_pv', 'N PV', (0, 100), 1, 50],
    'n_trueInteractions': ['n_trueInteractions', 'N True Int.', (-2, 100), 1, 51],

    'pv0_numTracks': ['pv0_numTracks', 'PV0 N Tracks', 'pv_numTracks[0]', (0, 250), 1, 50],
    'pv0_normalizedChi2': ['pv0_normalizedChi2', 'PV0 Normalized Chi2', 'pv_normalizedChi2[0]', (0, 10), 1, 50],
    'pv0_z': ['pv0_z', 'PV0 z (cm)', 'pv_z[0]', (-30, 30), 1, 60],

    'met_pt': ['met_pt', 'p_{T}^{miss} (GeV)', (0, 2500), 1, 50],  # [0, 250, 300, 350, 400, 450, 500, 600, 700, 1000, 1500, 2500]
    'met_phi': ['met_phi', '#phi(p_{T}^{miss})', (-4, 4), 1, 40],

    'ht': ['ht', 'H_{T} (GeV)', (0, 5000), 1, 50],  # (jets w/ p_{T} > 30 GeV, |#eta| < 2.4)
    'ht5': ['ht5', 'H_{T}^{5} (GeV)', (0, 5000), 1, 50],  # (jets w/ p_{T} > 30 GeV, |#eta| < 5)
    'ht5overht': ['ht5overht', 'H_{T}(|#eta| < 5) / H_{T}(|#eta| < 2.4) ', 'ht5/ht', (0, 10), 1, 50],
    'htMiss': ['htMiss', 'H_{T}^{miss} (GeV)', (0, 2500), 1, 50],  # (jets w/ p_{T} > 30 GeV, |#eta| < 5)

    'n_jet_15': ['n_jet_15', 'N Jet (p_{T} > 15 GeV)', (0, 100), 1, 50],
    'n_jet_30': ['n_jet_30', 'N Jet (p_{T} > 30 GeV)', (0, 20), 1, 20],
    'n_jet_50': ['n_jet_50', 'N Jet (p_{T} > 50 GeV)', (0, 20), 1, 20],
    'n_jet_100': ['n_jet_100', 'N Jet (p_{T} > 100 GeV)', (0, 10), 1, 10],
    'n_jet_200': ['n_jet_200', 'N Jet (p_{T} > 200 GeV)', (0, 10), 1, 10],

    'n_jet_15_btagloose': ['n_jet_15_btagloose', 'N Jet (p_{T} > 15 GeV, b-Tagged Loose WP)', (0, 10), 1, 10],
    'n_jet_15_btagmedium': ['n_jet_15_btagmedium', 'N Jet (p_{T} > 15 GeV, b-Tagged Medium WP)', (0, 10), 1, 10],
    'n_jet_15_btagtight': ['n_jet_15_btagtight', 'N Jet (p_{T} > 15 GeV, b-Tagged Tight WP)', (0, 10), 1, 10],
    
    'n_jet_30_btagloose': ['n_jet_30_btagloose', 'N Jet (p_{T} > 30 GeV, b-Tagged Loose WP)', (0, 10), 1, 10],
    'n_jet_30_btagmedium': ['n_jet_30_btagmedium', 'N Jet (p_{T} > 30 GeV, b-Tagged Medium WP)', (0, 10), 1, 10],
    'n_jet_30_btagtight': ['n_jet_30_btagtight', 'N Jet (p_{T} > 30 GeV, b-Tagged Tight WP)', (0, 10), 1, 10],

    'mtMetLeadingJet': ['mtMetLeadingJet', 'm_{T}(p_{T}^{miss}, Leading Jet) (GeV)', (0, 2500), 1, 50],
    'dphiminMetJets': ['dphiminMetJets', '#Delta#phi_{min}(p_{T}^{miss}, Jet_{1,2,3,4})', (0, 3.5), 1, 35],

    'n_photon': ['n_photon', 'N Photon', (0, 10), 1, 10],
    'n_photon_iso': ['n_photon_iso', 'N Isolated Photon', (0, 10), 1, 10],
    
    'n_electron': ['n_electron', 'N Electron', (0, 10), 1, 10],
    'n_electron_iso': ['n_electron_iso', 'N Isolated Electron', (0, 10), 1, 10],
    
    'n_muon': ['n_muon', 'N Muon', (0, 10), 1, 10],
    'n_muon_iso': ['n_muon_iso', 'N Isolated Muon', (0, 10), 1, 10],
    
    'n_lepton': ['n_lepton', 'N Lepton', (0, 10), 1, 10],
    'n_lepton_iso': ['n_lepton_iso', 'N Isolated Lepton', (0, 10), 1, 10],

    'n_tau_vloose': ['n_tau_vloose', 'N Tau VLoose WP (No p_{T} Cut)', (0, 20), 1, 20],
    'n_tau_loose': ['n_tau_loose', 'N Tau Loose WP (No p_{T} Cut)', (0, 20), 1, 20],
    'n_tau_medium': ['n_tau_medium', 'N Tau Medium WP (No p_{T} Cut)', (0, 20), 1, 20],
    'n_tau_tight': ['n_tau_tight', 'N Tau Tight WP (No p_{T} Cut)', (0, 20), 1, 20],
    'n_tau_vtight': ['n_tau_vtight', 'N Tau VTight WP (No p_{T} Cut)', (0, 20), 1, 20],
    'n_tau_vvtight': ['n_tau_vvtight', 'N Tau VVTight WP (No p_{T} Cut)', (0, 20), 1, 20],

    'n_tau_20_vloose': ['n_tau_20_vloose', 'N Tau VLoose WP', (0, 10), 1, 10],
    'n_tau_20_loose': ['n_tau_20_loose', 'N Tau Loose WP', (0, 10), 1, 10],
    'n_tau_20_medium': ['n_tau_20_medium', 'N Tau Medium WP', (0, 10), 1, 10],
    'n_tau_20_tight': ['n_tau_20_tight', 'N Tau Tight WP', (0, 10), 1, 10],
    'n_tau_20_vtight': ['n_tau_20_vtight', 'N Tau VTight WP', (0, 10), 1, 10],
    'n_tau_20_vvtight': ['n_tau_20_vvtight', 'N Tau VVTight WP', (0, 10), 1, 10],

    'n_track': ['n_track', 'N Track (Preselection)', (0, 1000), 1, 50],

    'n_track_total_per_pv': ['n_track_total_per_pv', 'N Track / PV', 'n_track_total/n_pv', (0, 200), 1, 50],


    'deltam': ['deltam', '', (0, 1.5), 1, 150],
    'deltam_true': ['deltam_true', '', (0, 1.5), 1, 150],


    'track_charge': ['track_charge', 'Charge', (-1, 2), 1, 3],
    'track_pt': ['track_pt', 'p_{T} (GeV)', (0, 6), 1, 60],
    'track_pt_widerange': ['track_pt_widerange', 'p_{T} (GeV)', 'track_pt', (0, 30), 1, 60],
    'track_ptError_D_pt': ['track_ptError_D_pt', 'p^{Error}_{T} / p_{T}', (0, 0.5), 1, 50],
    'track_log10_ptError_D_pt_': ['track_log10_ptError_D_pt_', 'log_{10}(p^{Error}_{T} / p_{T})', (-3, 0), 1, 60],
    'track_eta': ['track_eta', '#eta', (-4, 4), 1, 40],
    'track_abs_eta_': ['track_abs_eta_', '|#eta|', 'abs(track_eta)', (0, 4), 1, 40],
    'track_phi': ['track_phi', '#phi', (-4, 4), 1, 40],
    'track_etaError': ['track_etaError', '#eta^{Error}', (0, 0.1), 1, 50],
    'track_phiError': ['track_phiError', '#phi^{Error}', (0, 0.1), 1, 50],

    'track_isPfCand': ['track_isPfCand', 'Track is PFC', (0, 2), 1, 2],
    'track_abs_pfCandPdgId_': ['track_abs_pfCandPdgId_', 'PdgId PFC', 'abs(track_pfCandPdgId)', (0, 250), 1, 250],
    'track_pfCandPt': ['track_pfCandPt', 'p_{T}^{PFC}', (-1, 10), 1, 60],
    'track_pfCandEnergy': ['track_pfCandEnergy', 'E^{PFC}', (-1, 20), 1, 42],
    'track_pfCandEcalEnergy': ['track_pfCandEcalEnergy', 'E_{ECAL}^{PFC}', (-1, 20), 1, 42],
    'track_pfCandHcalEnergy': ['track_pfCandHcalEnergy', 'E_{HCAL}^{PFC}', (-1, 20), 1, 42],

    # min(...) doesn't work with RDF "vectorized"
    # 'track_min_1_associatedPV_': ['track_min_1_associatedPV_', 'Associated Vertex (None, PV, PU)', 'min(1.,track_associatedPV)', (-1, 2), 1, 3],
    'track_associatedPV': ['track_associatedPV', 'Associated Vertex', (-1, 20), 1, 21],

    'track_IPsig': ['track_IPsig', 'IP Significance w.r.t. PV', (0, 150), 1, 30],
    'track_IPxyz': ['track_IPxyz', 'IPxyz w.r.t. PV (cm)', (0, 2), 1, 40],
    'track_IPxy': ['track_IPxy', 'IPxy w.r.t. PV (cm)', (0, 2), 1, 40],
    'track_IPz': ['track_IPz', 'IPz w.r.t. PV (cm)', (0, 2), 1, 40],

    'track_log10_IPsig_': ['track_log10_IPsig_', 'log_{10}(IP Significance w.r.t. PV)', (-2, 4), 1, 30],
    'track_log10_IPxyz_': ['track_log10_IPxyz_', 'log_{10}(IPxyz w.r.t. PV (cm))', (-6, 2), 1, 40],
    'track_log10_IPxy_': ['track_log10_IPxy_', 'log_{10}(IPxy w.r.t. PV (cm))', (-6, 2), 1, 40],
    'track_log10_IPz_': ['track_log10_IPz_', 'log_{10}(IPz w.r.t. PV (cm))', (-6, 2), 1, 40],
    
    'track_IPsigPU': ['track_IPsigPU', 'IP Significance w.r.t. PU', (0, 150), 1, 30],
    'track_IPxyzPU': ['track_IPxyzPU', 'IPxyz w.r.t. PU (cm)', (0, 2), 1, 40],
    'track_IPxyPU': ['track_IPxyPU', 'IPxy w.r.t. PU (cm)', (0, 2), 1, 40],
    'track_IPzPU': ['track_IPzPU', 'IPz w.r.t. PU (cm)', (0, 2), 1, 40],

    'track_log10_IPsigPU_': ['track_log10_IPsigPU_', 'log_{10}(IP Significance w.r.t. PU)', (-2, 4), 1, 30],
    'track_log10_IPxyzPU_': ['track_log10_IPxyzPU_', 'log_{10}(IPxyz w.r.t. PU (cm))', (-6, 2), 1, 40],
    'track_log10_IPxyPU_': ['track_log10_IPxyPU_', 'log_{10}(IPxy w.r.t. PU (cm))', (-6, 2), 1, 40],
    'track_log10_IPzPU_': ['track_log10_IPzPU_', 'log_{10}(IPz w.r.t. PU (cm))', (-6, 2), 1, 40],

    'track_dxy': ['track_dxy', 'dxy w.r.t. PV (cm)', (0, 2), 1, 40],
    'track_dz': ['track_dz', 'dz w.r.t. PV (cm)', (0, 2), 1, 40],
    'track_log10_dxy_': ['track_log10_dxy_', 'log_{10}(dxy w.r.t. PV (cm))', (-6, 2), 1, 40],
    'track_log10_dz_': ['track_log10_dz_', 'log_{10}(dz w.r.t. PV (cm))', (-6, 2), 1, 40],  # [-6, -3., -0.2, 0, 2]

    'track_dxyPU': ['track_dxyPU', 'dxy w.r.t. PU (cm)', (0, 2), 1, 40],
    'track_dzPU': ['track_dzPU', 'dz w.r.t. PU (cm)', (0, 2), 1, 40],
    'track_log10_dxyPU_': ['track_log10_dxyPU_', 'log_{10}(dxy w.r.t. PU (cm))', (-6, 2), 1, 40],
    'track_log10_dzPU_': ['track_log10_dzPU_', 'log_{10}(dz w.r.t. PU (cm))', (-6, 2), 1, 40],

    'track_log10_dxyError_': ['track_log10_dxyError_', 'log_{10}(dxy^{Error} (cm))', (-3, 1), 1, 40],
    'track_log10_dzError_': ['track_log10_dzError_', 'log_{10}(dz^{Error} (cm))', (-3, 1), 1, 40],

    'track_pfAbsIso': ['track_pfAbsIso', 'Abs. Iso PF', (-1, 30), 1, 31],
    'track_pfRelIso': ['track_pfRelIso', 'Rel. Iso PF', (-1, 50), 1, 51],
    'track_drminPf': ['track_drminPf', '#DeltaR_{min}^{PF}', (0, 1), 1, 50],
    'track_numneighboursPf': ['track_numneighboursPf', 'N Neighbours PF', (0, 30), 1, 30],
    
    'track_chPfAbsIso': ['track_chPfAbsIso', 'Abs. Iso Ch. PF', (-1, 30), 1, 31],
    'track_chPfRelIso': ['track_chPfRelIso', 'Rel. Iso Ch. PF', (-1, 50), 1, 51],
    'track_drminChPf': ['track_drminChPf', '#DeltaR_{min}^{Ch. PF}', (0, 10), 1, 50],
    'track_numneighboursChPf': ['track_numneighboursChPf', 'N Neighbours Ch. PF', (0, 30), 1, 30],
    
    'track_tkAbsIso0': ['track_tkAbsIso0', 'Abs. Iso Tracks p_{T} > 0 GeV', (-1, 30), 1, 31],
    'track_tkRelIso0': ['track_tkRelIso0', 'Rel. Iso Tracks p_{T} > 0 GeV', (-1, 50), 1, 51],
    'track_drminTrack0': ['track_drminTrack0', '#DeltaR_{min}^{Tracks pT > 0 GeV}', (0, 10), 1, 50],
    'track_numneighboursTrack0': ['track_numneighboursTrack0', 'N Neighbours Tracks p_{T} > 0 GeV', (0, 30), 1, 30],
    
    'track_tkAbsIso1': ['track_tkAbsIso1', 'Abs. Iso Tracks p_{T} > 1 GeV', (-1, 30), 1, 31],
    'track_tkRelIso1': ['track_tkRelIso1', 'Rel. Iso Tracks p_{T} > 1 GeV', (-1, 50), 1, 51],
    'track_drminTrack1': ['track_drminTrack1', '#DeltaR_{min}^{Tracks pT > 1 GeV}', (0, 10), 1, 50],
    'track_numneighboursTrack1': ['track_numneighboursTrack1', 'N Neighbours Tracks p_{T} > 1 GeV', (0, 30), 1, 30],
    
    'track_tkAbsIso10': ['track_tkAbsIso10', 'Abs. Iso Tracks p_{T} > 10 GeV', (-1, 30), 1, 31],
    'track_tkRelIso10': ['track_tkRelIso10', 'Rel. Iso Tracks p_{T} > 10 GeV', (-1, 50), 1, 51],
    'track_drminTrack10': ['track_drminTrack10', '#DeltaR_{min}^{Tracks pT > 10 GeV}', (0, 10), 1, 50],
    'track_numneighboursTrack10': ['track_numneighboursTrack10', 'N Neighbours Tracks p_{T} > 10 GeV', (0, 30), 1, 30],

    # 'track_jetIso0': ['track_jetIso0', 'Jet Iso p_{T} > 0 GeV', (0, 200), 1, 50],
    # 'track_jetIsoMulti0': ['track_jetIsoMulti0', 'Jet Iso Multi p_{T} > 0 GeV', (0, 100), 1, 50],
    # 'track_drminJet0': ['track_drminJet0', '#DeltaR_{min}^{jets pT > 0 GeV}', (0, 10), 1, 50],
    # 'track_btagJet0': ['track_btagJet0', 'b-tag closest jet p_{T} > 0 GeV', (-2, 1), 1, 30],
    
    # 'track_jetIso10': ['track_jetIso10', 'Jet Iso p_{T} > 10 GeV', (0, 200), 1, 50],
    # 'track_jetIsoMulti10': ['track_jetIsoMulti10', 'Jet Iso Multi p_{T} > 10 GeV', (0, 100), 1, 50],
    # 'track_drminJet10': ['track_drminJet10', '#DeltaR_{min}^{jets pT > 10 GeV}', (0, 10), 1, 50],
    # 'track_btagJet10': ['track_btagJet10', 'b-tag closest jet p_{T} > 10 GeV', (-2, 1), 1, 30],
    
    'track_jetIso15': ['track_jetIso15', 'Jet Iso p_{T} > 15 GeV', (0, 200), 1, 50],
    'track_jetIsoMulti15': ['track_jetIsoMulti15', 'Jet Iso Multi p_{T} > 15 GeV', (0, 100), 1, 50],
    'track_drminJet15': ['track_drminJet15', '#DeltaR_{min}^{jets pT > 15 GeV}', (0, 10), 1, 50],
    'track_btagJet15': ['track_btagJet15', 'b-tag closest jet p_{T} > 15 GeV', (-2, 1), 1, 30],
    
    'track_jetIso20': ['track_jetIso20', 'Jet Iso p_{T} > 20 GeV', (0, 200), 1, 50],
    'track_jetIsoMulti20': ['track_jetIsoMulti20', 'Jet Iso Multi p_{T} > 20 GeV', (0, 100), 1, 50],
    'track_drminJet20': ['track_drminJet20', '#DeltaR_{min}^{jets pT > 20 GeV}', (0, 10), 1, 50],
    'track_btagJet20': ['track_btagJet20', 'b-tag closest jet p_{T} > 20 GeV', (-2, 1), 1, 30],
    
    'track_jetIso30': ['track_jetIso30', 'Jet Iso p_{T} > 30 GeV', (0, 200), 1, 50],
    'track_jetIsoMulti30': ['track_jetIsoMulti30', 'Jet Iso Multi p_{T} > 30 GeV', (0, 100), 1, 50],
    'track_drminJet30': ['track_drminJet30', '#DeltaR_{min}^{jets pT > 30 GeV}', (0, 10), 1, 50],
    'track_btagJet30': ['track_btagJet30', 'b-tag closest jet p_{T} > 30 GeV', (-2, 1), 1, 30],

    'track_jetIsoNoLepton15': ['track_jetIsoNoLepton15', 'Jet Iso p_{T} > 15 GeV no lepton', (0, 200), 1, 50],
    'track_jetIsoMultiNoLepton15': ['track_jetIsoMultiNoLepton15', 'Jet Iso Multi p_{T} > 15 GeV no lepton', (0, 100), 1, 50],
    'track_drminJetNoLepton15': ['track_drminJetNoLepton15', '#DeltaR_{min}^{jets pT > 15 GeV} no lepton', (0, 10), 1, 50],
    'track_btagJetNoLepton15': ['track_btagJetNoLepton15', 'b-tag closest jet p_{T} > 15 GeV no lepton', (-2, 1), 1, 30],

    'track_neHadAbsIso0': ['track_neHadAbsIso0', 'Abs. Iso ne. had. p_{T} > 0 GeV', (0, 50), 1, 50],
    'track_drminNeHad0': ['track_drminNeHad0', '#DeltaR_{min}^{ne. had. pT > 0 GeV}', (0, 10), 1, 50],
    # 'track_invmNeHad0': ['track_invmNeHad0', 'Inv. Mass(track, ne. had. p_{T} > 0 GeV)', (0, 50), 1, 50],
    
    'track_neHadAbsIso1': ['track_neHadAbsIso1', 'Abs. Iso ne. had. p_{T} > 1 GeV', (0, 50), 1, 50],
    'track_drminNeHad1': ['track_drminNeHad1', '#DeltaR_{min}^{ne. had. pT > 1 GeV}', (0, 10), 1, 50],
    # 'track_invmNeHad1': ['track_invmNeHad1', 'Inv. Mass(track, ne. had. p_{T} > 1 GeV) (GeV)', (0, 50), 1, 50],
    
    'track_neHadAbsIso10': ['track_neHadAbsIso10', 'Abs. Iso ne. had. p_{T} > 10 GeV', (0, 50), 1, 50],
    'track_drminNeHad10': ['track_drminNeHad10', '#DeltaR_{min}^{ne. had. pT > 10 GeV}', (0, 10), 1, 50],
    # 'track_invmNeHad10': ['track_invmNeHad10', 'Inv. Mass(track, ne. had. p_{T} > 10 GeV) (GeV)', (0, 50), 1, 50],

    'track_drminPhoton': ['track_drminPhoton', '#DeltaR_{min}^{photon}', (0, 10), 1, 50],
    'track_drminElectron': ['track_drminElectron', '#DeltaR_{min}^{electron}', (0, 10), 1, 50],
    'track_drminMuon': ['track_drminMuon', '#DeltaR_{min}^{muon}', (0, 10), 1, 50],

    'track_isTauLeadPfChHadCand0': ['track_isTauLeadPfChHadCand0', 'isLeadPfChHadCand Tau p_{T} > 0 GeV', (0, 2), 1, 2],
    'track_drminTau0': ['track_drminTau0', '#DeltaR_{min}^{tau pT > 0 GeV}', (0, 10), 1, 50],
    'track_mvaDiscrTau0': ['track_mvaDiscrTau0', 'MVA Discr. Closest Tau p_{T} > 0 GeV', (-1, 1), 1, 50],
    'track_decayModeTau0': ['track_decayModeTau0', 'DM Closest Tau p_{T} > 0 GeV', (-1, 25), 1, 26],
    'track_dr3highestWpTau0': ['track_dr3highestWpTau0', 'Highest WP #DeltaR<0.3 Tau p_{T} > 0 GeV', (0, 10), 1, 10],
    'track_dr4highestWpTau0': ['track_dr4highestWpTau0', 'Highest WP #DeltaR<0.4 Tau p_{T} > 0 GeV', (0, 10), 1, 10],
    'track_dr5highestWpTau0': ['track_dr5highestWpTau0', 'Highest WP #DeltaR<0.5 Tau p_{T} > 0 GeV', (0, 10), 1, 10],
    
    'track_isTauLeadPfChHadCand20': ['track_isTauLeadPfChHadCand20', 'isLeadPfChHadCand Tau p_{T} > 20 GeV', (0, 2), 1, 2],
    'track_drminTau20': ['track_drminTau20', '#DeltaR_{min}^{tau pT > 20 GeV}', (0, 10), 1, 50],
    'track_mvaDiscrTau20': ['track_mvaDiscrTau20', 'MVA Discr. Closest Tau p_{T} > 20 GeV', (-1, 1), 1, 50],
    'track_decayModeTau20': ['track_decayModeTau20', 'DM Closest Tau p_{T} > 20 GeV', (-1, 25), 1, 26],
    'track_dr3highestWpTau20': ['track_dr3highestWpTau20', 'Highest WP #DeltaR<0.3 Tau p_{T} > 20 GeV', (0, 10), 1, 10],
    'track_dr4highestWpTau20': ['track_dr4highestWpTau20', 'Highest WP #DeltaR<0.4 Tau p_{T} > 20 GeV', (0, 10), 1, 10],
    'track_dr5highestWpTau20': ['track_dr5highestWpTau20', 'Highest WP #DeltaR<0.5 Tau p_{T} > 20 GeV', (0, 10), 1, 10],

    # 'track_detaLeadingJet': ['track_detaLeadingJet', '#Delta#eta(Track, Leading Jet)', (0, 6), 1, 30],
    'track_abs_detaLeadingJet_': ['track_abs_detaLeadingJet_', '|#Delta#eta(Track, Leading Jet)|', 'abs(track_detaLeadingJet)', (0, 6), 1, 30],
    # 'track_dphiLeadingJet': ['track_dphiLeadingJet', '#Delta#phi(Track, Leading Jet)', (-4, 4), 1, 40],
    'track_abs_dphiLeadingJet_': ['track_abs_dphiLeadingJet_', '|#Delta#phi(Track, Leading Jet)|', 'abs(track_dphiLeadingJet)', (0, 4), 1, 40],
    'track_dphiMet': ['track_dphiMet', '#Delta#phi(Track, p_{T}^{miss})', (-4, 4), 1, 40],
    'track_abs_dphiMet_': ['track_abs_dphiMet_', '|#Delta#phi(Track, p_{T}^{miss})|', 'abs(track_dphiMet)', (0, 4), 1, 40],
    'track_dphiMetPca': ['track_dphiMetPca', '#Delta#phi^{PCA}(Track, p_{T}^{miss})', (-4, 4), 1, 40],
    'track_abs_dphiMetPca_': ['track_abs_dphiMetPca_', '|#Delta#phi^{PCA}(Track, p_{T}^{miss})|', 'abs(track_dphiMetPca)', (0, 4), 1, 40],

    'track_chi2': ['track_chi2', '#chi^{2}', (0, 10), 1, 50],
    'track_quality': ['track_quality', 'Quality', (1, 4), 1, 3],
    'track_numValidHits': ['track_numValidHits', 'N Valid Hits', (0, 40), 1, 40],
    'track_numLostHits': ['track_numLostHits', 'N Lost Hits', (0, 10), 1, 10],

    'track_hasGenMatch': ['track_hasGenMatch', 'Has GEN Match', (0, 3), 1, 3],
    # 'track_abs_genMatchPdgId_': ['track_abs_genMatchPdgId_', '|GEN Match PdgId|', 'abs(track_genMatchPdgId)', (0, 5000), 1, 100],
    'track_genMatchPt': ['track_genMatchPt', 'GEN Match p_{T} (GeV)', (-1, 50), 1, 51],
    'track_genMatchStatus': ['track_genMatchStatus', 'GEN Match Status', (-1, 60), 1, 61],
    'track_genMatchIsHardProcess': ['track_genMatchIsHardProcess', 'GEN Match isHardProcess', (-1, 2), 1, 3],
    'track_genMatchIsFromHardProcess': ['track_genMatchIsFromHardProcess', 'GEN Match isFromHardProcess', (-1, 2), 1, 3],
    'track_genMatchIsPrompt': ['track_genMatchIsPrompt', 'GEN Match isPrompt', (-1, 2), 1, 3],
    'track_genMatchIsDirectHadronDecayProduct': ['track_genMatchIsDirectHadronDecayProduct', 'GEN Match isDirectHadronDecayProduct', (-1, 2), 1, 3],
    'track_genMatchIsDirectTauDecayProduct': ['track_genMatchIsDirectTauDecayProduct', 'GEN Match isDirectTauDecayProduct', (-1, 2), 1, 3],
    # 'track_abs_genMatchMotherPdgId_': ['track_abs_genMatchMotherPdgId_', '|GEN Match Mother PdgId|', 'abs(track_genMatchMotherPdgId)', (0, 110000), 1, 110],
    'track_genMatchMotherPt': ['track_genMatchMotherPt', 'GEN Match Mother p_{T} (GeV)', (-5, 300), 1, 61],
    'track_genMatchMotherStatus': ['track_genMatchMotherStatus', 'GEN Match Mother Status', (-1, 80), 1, 81],
    'track_genMatchMotherIsHardProcess': ['track_genMatchMotherIsHardProcess', 'GEN Match Mother isHardProcess', (-1, 2), 1, 3],
    'track_genMatchMotherIsTheTau': ['track_genMatchMotherIsTheTau', 'GEN Match Mother isTheTau', (-1, 2), 1, 3],
    'track_genMatchMotherTauDecay': ['track_genMatchMotherTauDecay', 'GEN Match Mother TauDecay', (-1, 25), 1, 26],

    'track_drminGenTauJet': ['track_drminGenTauJet', '#DeltaR_{min}^{GEN Tau}', (0, 10), 1, 50],
}

myvariablesV10 = myvariables.copy()
myvariablesV10.update({

    'track_IPsigAssPV': ['track_IPsigAssPV', 'IP Significance w.r.t. Ass. PV', (0, 150), 1, 30],
    'track_IPxyzAssPV': ['track_IPxyzAssPV', 'IPxyz w.r.t. Ass. PV (cm)', (0, 2), 1, 40],
    'track_IPxyAssPV': ['track_IPxyAssPV', 'IPxy w.r.t. Ass. PV (cm)', (0, 2), 1, 40],
    'track_IPzAssPV': ['track_IPzAssPV', 'IPz w.r.t. Ass. PV (cm)', (0, 2), 1, 40],

    'track_log10_IPsigAssPV_': ['track_log10_IPsigAssPV_', 'log_{10}(IP Significance w.r.t. Ass. PV)', (-2, 4), 1, 30],
    'track_log10_IPxyzAssPV_': ['track_log10_IPxyzAssPV_', 'log_{10}(IPxyz w.r.t. Ass. PV (cm))', (-6, 2), 1, 40],
    'track_log10_IPxyAssPV_': ['track_log10_IPxyAssPV_', 'log_{10}(IPxy w.r.t. Ass. PV (cm))', (-6, 2), 1, 40],
    'track_log10_IPzAssPV_': ['track_log10_IPzAssPV_', 'log_{10}(IPz w.r.t. Ass. PV (cm))', (-6, 2), 1, 40],

    'track_IPsigPUAssPV': ['track_IPsigPUAssPV', 'IP Significance w.r.t. PU (excl. ass. PV)', (0, 150), 1, 30],
    'track_IPxyzPUAssPV': ['track_IPxyzPUAssPV', 'IPxyz w.r.t. PU (excl. ass. PV) (cm)', (0, 2), 1, 40],
    'track_IPxyPUAssPV': ['track_IPxyPUAssPV', 'IPxy w.r.t. PU (excl. ass. PV) (cm)', (0, 2), 1, 40],
    'track_IPzPUAssPV': ['track_IPzPUAssPV', 'IPz w.r.t. PU (excl. ass. PV) (cm)', (0, 2), 1, 40],

    'track_log10_IPsigPUAssPV_': ['track_log10_IPsigPUAssPV_', 'log_{10}(IP Significance w.r.t. PU (excl. ass. PV))', (-2, 4), 1, 30],
    'track_log10_IPxyzPUAssPV_': ['track_log10_IPxyzPUAssPV_', 'log_{10}(IPxyz w.r.t. PU (excl. ass. PV) (cm))', (-6, 2), 1, 40],
    'track_log10_IPxyPUAssPV_': ['track_log10_IPxyPUAssPV_', 'log_{10}(IPxy w.r.t. PU (excl. ass. PV) (cm))', (-6, 2), 1, 40],
    'track_log10_IPzPUAssPV_': ['track_log10_IPzPUAssPV_', 'log_{10}(IPz w.r.t. PU (excl. ass. PV) (cm))', (-6, 2), 1, 40],
    
    'track_dxyAssPV': ['track_dxyAssPV', 'dxy w.r.t. Ass. PV (cm)', (0, 2), 1, 40],
    'track_dzAssPV': ['track_dzAssPV', 'dz w.r.t. Ass. PV (cm)', (0, 2), 1, 40],
    'track_log10_dxyAssPV_': ['track_log10_dxyAssPV_', 'log_{10}(dxy w.r.t. Ass. PV (cm))', (-6, 2), 1, 40],
    'track_log10_dzAssPV_': ['track_log10_dzAssPV_', 'log_{10}(dz w.r.t. Ass. PV (cm))', (-6, 2), 1, 40],

    'track_dxyPUAssPV': ['track_dxyPUAssPV', 'dxy w.r.t. PU (excl. ass. PV) (cm)', (0, 2), 1, 40],
    'track_dzPUAssPV': ['track_dzPUAssPV', 'dz w.r.t. PU (excl. ass. PV) (cm)', (0, 2), 1, 40],
    'track_log10_dxyPUAssPV_': ['track_log10_dxyPUAssPV_', 'log_{10}(dxy w.r.t. PU (excl. ass. PV) (cm))', (-6, 2), 1, 40],
    'track_log10_dzPUAssPV_': ['track_log10_dzPUAssPV_', 'log_{10}(dz w.r.t. PU (excl. ass. PV) (cm))', (-6, 2), 1, 40],

    'track_tkAbsIso5': ['track_tkAbsIso5', 'Abs. Iso Tracks p_{T} > 5 GeV', (-1, 30), 1, 31],
    'track_tkRelIso5': ['track_tkRelIso5', 'Rel. Iso Tracks p_{T} > 5 GeV', (-1, 50), 1, 51],
    'track_drminTrack5': ['track_drminTrack5', '#DeltaR_{min}^{tracks pT > 5 GeV}', (0, 10), 1, 50],
    'track_numneighboursTrack5': ['track_numneighboursTrack5', 'N neighbours tracks p_{T} > 5 GeV', (0, 30), 1, 30],

    # 'track_minvJet0': ['track_minvJet0', 'Inv. Mass(track, jet p_{T} > 0 GeV) (GeV)', (-5, 250), 1, 51],
    # 'track_minvJet10': ['track_minvJet10', 'Inv. Mass(track, jet p_{T} > 10 GeV) (GeV)', (-5, 250), 1, 51],
    # 'track_minvJet15': ['track_minvJet15', 'Inv. Mass(track, jet p_{T} > 15 GeV) (GeV)', (-5, 250), 1, 51],
    # 'track_minvJet20': ['track_minvJet20', 'Inv. Mass(track, jet p_{T} > 20 GeV) (GeV)', (-5, 250), 1, 51],
    # 'track_minvJet30': ['track_minvJet30', 'Inv. Mass(track, jet p_{T} > 30 GeV) (GeV)', (-5, 250), 1, 51],
    
    # 'track_minvJetNoLepton15': ['track_minvJetNoLepton15', 'Inv. Mass(track, jet p_{T} > 15 GeV) no lepton (GeV)', (-5, 250), 1, 51],

    'track_bjetLooseIso15': ['track_bjetLooseIso15', 'B-Jet Loose Iso p_{T} > 15 GeV', (0, 200), 1, 50],
    'track_bjetLooseIsoMulti15': ['track_bjetLooseIsoMulti15', 'B-Jet Loose Iso Multi p_{T} > 15 GeV', (0, 100), 1, 50],
    'track_drminBjetLoose15': ['track_drminBjetLoose15', '#DeltaR_{min}^{b-jets loose pT > 15 GeV}', (0, 10), 1, 50],
    'track_btagBjetLoose15': ['track_btagBjetLoose15', 'b-tag closest b-jet loose p_{T} > 15 GeV', (-2, 1), 1, 30],
    # 'track_minvBjetLoose15': ['track_minvBjetLoose15', 'Inv. Mass(track, b-jet loose p_{T} > 15 GeV) (GeV)', (-5, 250), 1, 51],
    
    'track_bjetLooseIso30': ['track_bjetLooseIso30', 'B-Jet Loose Iso p_{T} > 30 GeV', (0, 200), 1, 50],
    'track_bjetLooseIsoMulti30': ['track_bjetLooseIsoMulti30', 'B-Jet Loose Iso Multi p_{T} > 30 GeV', (0, 100), 1, 50],
    'track_drminBjetLoose30': ['track_drminBjetLoose30', '#DeltaR_{min}^{b-jets loose pT > 30 GeV}', (0, 10), 1, 50],
    'track_btagBjetLoose30': ['track_btagBjetLoose30', 'b-tag closest b-jet loose p_{T} > 30 GeV', (-2, 1), 1, 30],
    # 'track_minvBjetLoose30': ['track_minvBjetLoose30', 'Inv. Mass(track, b-jet loose p_{T} > 30 GeV) (GeV)', (-5, 250), 1, 51],
    
    'track_bjetMediumIso15': ['track_bjetMediumIso15', 'B-Jet Medium Iso p_{T} > 15 GeV', (0, 200), 1, 50],
    'track_bjetMediumIsoMulti15': ['track_bjetMediumIsoMulti15', 'B-Jet Medium Iso Multi p_{T} > 15 GeV', (0, 100), 1, 50],
    'track_drminBjetMedium15': ['track_drminBjetMedium15', '#DeltaR_{min}^{b-jets Medium pT > 15 GeV}', (0, 10), 1, 50],
    'track_btagBjetMedium15': ['track_btagBjetMedium15', 'b-tag closest b-jet medium p_{T} > 15 GeV', (-2, 1), 1, 30],
    # 'track_minvBjetMedium15': ['track_minvBjetMedium15', 'Inv. Mass(track, b-jet medium p_{T} > 15 GeV) (GeV)', (-5, 250), 1, 51],
    
    'track_bjetMediumIso30': ['track_bjetMediumIso30', 'B-Jet Medium Iso p_{T} > 30 GeV', (0, 200), 1, 50],
    'track_bjetMediumIsoMulti30': ['track_bjetMediumIsoMulti30', 'B-Jet Medium Iso Multi p_{T} > 30 GeV', (0, 100), 1, 50],
    'track_drminBjetMedium30': ['track_drminBjetMedium30', '#DeltaR_{min}^{b-jets Medium pT > 30 GeV}', (0, 10), 1, 50],
    'track_btagBjetMedium30': ['track_btagBjetMedium30', 'b-tag closest b-jet medium p_{T} > 30 GeV', (-2, 1), 1, 30],
    # 'track_minvBjetMedium30': ['track_minvBjetMedium30', 'Inv. Mass(track, b-jet medium p_{T} > 30 GeV) (GeV)', (-5, 250), 1, 51],
    
    'track_bjetTightIso15': ['track_bjetTightIso15', 'B-Jet Tight Iso p_{T} > 15 GeV', (0, 200), 1, 50],
    'track_bjetTightIsoMulti15': ['track_bjetTightIsoMulti15', 'B-Jet Tight Iso Multi p_{T} > 15 GeV', (0, 100), 1, 50],
    'track_drminBjetTight15': ['track_drminBjetTight15', '#DeltaR_{min}^{b-jets Tight pT > 15 GeV}', (0, 10), 1, 50],
    'track_btagBjetTight15': ['track_btagBjetTight15', 'b-tag closest b-jet tight p_{T} > 15 GeV', (-2, 1), 1, 30],
    # 'track_minvBjetTight15': ['track_minvBjetTight15', 'Inv. Mass(track, b-jet tight p_{T} > 15 GeV) (GeV)', (-5, 250), 1, 51],
    
    'track_bjetTightIso30': ['track_bjetTightIso30', 'B-Jet Tight Iso p_{T} > 30 GeV', (0, 200), 1, 50],
    'track_bjetTightIsoMulti30': ['track_bjetTightIsoMulti30', 'B-Jet Tight Iso Multi p_{T} > 30 GeV', (0, 100), 1, 50],
    'track_drminBjetTight30': ['track_drminBjetTight30', '#DeltaR_{min}^{b-jets Tight pT > 30 GeV}', (0, 10), 1, 50],
    'track_btagBjetTight30': ['track_btagBjetTight30', 'b-tag closest b-jet tight p_{T} > 30 GeV', (-2, 1), 1, 30],
    # 'track_minvBjetTight30': ['track_minvBjetTight30', 'Inv. Mass(track, b-jet tight p_{T} > 30 GeV) (GeV)', (-5, 250), 1, 51],
    
    'track_neHadAbsIso5': ['track_neHadAbsIso5', 'Abs. Iso ne. had. p_{T} > 5 GeV', (0, 50), 1, 50],
    'track_drminNeHad5': ['track_drminNeHad5', '#DeltaR_{min}^{ne. had. pT > 5 GeV}', (0, 10), 1, 50],
    # 'track_invmNeHad5': ['track_invmNeHad5', 'Inv. Mass(track, ne. had. p_{T} > 5 GeV) (GeV)', (0, 50), 1, 50],
})

myvariablesV14 = myvariablesV10.copy()
myvariablesV14.update({

    'zGamma_pt': ['zGamma_pt', 'p_{T}(Z_{GEN}) (GeV)', (-50, 2500), 1, 51],
    'zGamma_eta': ['zGamma_eta', '#eta(Z_{GEN})', (-4, 4), 1, 40],
    'wBoson_pt': ['wBoson_pt', 'p_{T}(W_{GEN}) (GeV)', (-20, 1000), 1, 51],
    'wBoson_eta': ['wBoson_eta', '#eta(W_{GEN})', (-4, 4), 1, 40],
    'wBoson_tauPt': ['wBoson_tauPt', 'p_{T}(#tau_{GEN}) (GeV)', (-20, 1000), 1, 51],
    'wBoson_tauEta': ['wBoson_tauEta', '#eta(#tau_{GEN})', (-4, 4), 1, 40],
    'wBoson_tauPtVis': ['wBoson_tauPtVis', 'p_{T}(#tau_{vis,GEN}) (GeV)', (-20, 1000), 1, 51],
    'wBoson_tauEtaVis': ['wBoson_tauEtaVis', '#eta(#tau_{vis,GEN})', (-4, 4), 1, 40],

    'cleaning_invm': ['cleaning_invm', 'DY Cleaning m(Z) (GeV)', (-3, 150), 1, 51],
    'cleaning_zPt': ['cleaning_zPt', 'DY Cleaning p_{T}(Z) (GeV)', (-50, 2500), 1, 51],
    'cleaning_metPtBeforeCleaning': ['cleaning_metPtBeforeCleaning', 'p_{T}^{miss} before DY Cleaning (GeV)', (-20, 1000), 1, 51],
    'cleaning_l1Pt': ['cleaning_l1Pt', 'DY Cleaning p_{T}(Lepton_{1}) (GeV)', (-20, 1000), 1, 51],
    'cleaning_l2Pt': ['cleaning_l2Pt', 'DY Cleaning p_{T}(Lepton_{2}) (GeV)', (-20, 1000), 1, 51],
    'cleaning_l1Eta': ['cleaning_l1Eta', 'DY Cleaning #eta(Lepton_{1})', (-4, 4), 1, 40],
    'cleaning_l2Eta': ['cleaning_l2Eta', 'DY Cleaning #eta(Lepton_{2})', (-4, 4), 1, 40],
    'cleaning_l1dBetaRelIso': ['cleaning_l1dBetaRelIso', 'DY Cleaning Rel. Iso. Lepton_{1}', (-1, 3), 1, 40],
    'cleaning_l2dBetaRelIso': ['cleaning_l2dBetaRelIso', 'DY Cleaning Rel. Iso. Lepton_{2}', (-1, 3), 1, 40],

    'maxbtagJet_pt': ['maxbtagJet_pt', 'Max. b-Tag Jet p_{T} (GeV)', (0, 1000), 1, 50],
    'maxbtagJet_eta': ['maxbtagJet_eta', 'Max. b-Tag Jet #eta', (-4, 4), 1, 40],
    'maxbtagJet_phi': ['maxbtagJet_phi', 'Max. b-Tag Jet #phi', (-4, 4), 1, 40],
    'maxbtagJet_btag': ['maxbtagJet_btag', 'Max. b-Tag Jet b-Tag CSVv2', (-2, 1), 1, 30],
    'maxbtagJet_btagDeepCSV': ['maxbtagJet_btagDeepCSV', 'Max. b-Tag Jet b-Tag DeepCSV', (-2, 1), 1, 30],

    'leadingJet_pt': ['leadingJet_pt', 'Leading Jet p_{T} (GeV)', (0, 1000), 1, 50],
    'leadingJet_eta': ['leadingJet_eta', 'Leading Jet #eta', (-4, 4), 1, 40],
    'leadingJet_phi': ['leadingJet_phi', 'Leading Jet #phi', (-4, 4), 1, 40],

    'n_jet_15_btagDeepCSVloose': ['n_jet_15_btagDeepCSVloose', 'N Jets b-Tag Loose WP (jets w/ p_{T} > 15 GeV)', (0, 10), 1, 10],
    'n_jet_15_btagDeepCSVmedium': ['n_jet_15_btagDeepCSVmedium', 'N Jets b-Tag Medium WP (jets w/ p_{T} > 15 GeV)', (0, 10), 1, 10],
    'n_jet_15_btagDeepCSVtight': ['n_jet_15_btagDeepCSVtight', 'N Jets b-Tag Tight WP (jets w/ p_{T} > 15 GeV)', (0, 10), 1, 10],
    
    'n_jet_30_btagDeepCSVloose': ['n_jet_30_btagDeepCSVloose', 'N Jets b-Tag Loose WP (jets w/ p_{T} > 30 GeV)', (0, 10), 1, 10],
    'n_jet_30_btagDeepCSVmedium': ['n_jet_30_btagDeepCSVmedium', 'N Jets b-Tag Medium WP (jets w/ p_{T} > 30 GeV)', (0, 10), 1, 10],
    'n_jet_30_btagDeepCSVtight': ['n_jet_30_btagDeepCSVtight', 'N Jets b-Tag Tight WP (jets w/ p_{T} > 30 GeV)', (0, 10), 1, 10],

    'track_drmin2ndTrack0': ['track_drmin2ndTrack0', '#DeltaR_{2nd min}^{tracks pT > 0 GeV}', (-1, 10), 1, 55],
    'track_drmin2ndTrack1': ['track_drmin2ndTrack1', '#DeltaR_{2nd min}^{tracks pT > 1 GeV}', (-1, 10), 1, 55],
    'track_drmin2ndTrack5': ['track_drmin2ndTrack5', '#DeltaR_{2nd min}^{tracks pT > 5 GeV}', (-1, 10), 1, 55],
    'track_drmin2ndTrack10': ['track_drmin2ndTrack10', '#DeltaR_{2nd min}^{tracks pT > 10 GeV}', (-1, 10), 1, 55],

    'track_distPVAssPVxy': ['track_distPVAssPVxy', '#Deltaxy(PV, associated PV) (cm)', (-0.1, 1), 1, 55],
    'track_distPVAssPVz': ['track_distPVAssPVz', '#Deltaz(PV, associated PV) (cm)', (-1, 11), 1, 60],

    'log10_crossSection_D_numSimEvents_': ['log10_crossSection_D_numSimEvents_', 'log_{10}(X-sec. / N sim. events)', 'log10(crossSection/numSimEvents)', (-5, 5), 1, 50],
})

myvariables_V15 = myvariablesV14.copy()
myvariables_V15.pop('triggerfired')
myvariables_V15.update({
    'triggerfired_met': ['triggerfired_met', 'MET/MHT Trigger Fired', (0, 2), 1, 2],
    'triggerfired_muon': ['triggerfired_muon', 'Muon Trigger Fired', (0, 2), 1, 2],

    'wBoson_tauDecayMode': ['wBoson_tauDecayMode', '#tau Decay Mode', (-1, 25), 1, 26],
    'wBoson_tauDecayMode_alt': ['wBoson_tauDecayMode_alt', '#tau Decay Mode', (-1, 25), 1, 26],

    'weight_genInfo': ['weight_genInfo', 'genInfo Weight', (0, 10), 1, 50],

    'era': ['era', 'Era', (-1, 4), 1, 5],

    'n_inclusivesv': ['n_inclusivesv', 'N Inclusive SVs', (0, 50), 1, 50],

    'n_jet_HEM1516veto': ['n_jet_HEM1516veto', 'N Jets HEM1516 Veto', (0, 10), 1, 10],
    'n_track_disaptrackcandidate': ['n_track_disaptrackcandidate', 'N Disap. Track Candidates', (0, 10), 1, 10],

    'track_genMatchEta': ['track_genMatchEta', 'GEN Match #eta', (-4, 4), 1, 40],
    'track_genMatchPhi': ['track_genMatchPhi', 'GEN Match #phi', (-4, 4), 1, 40],

    'track_genMatchDrmin': ['track_genMatchDrmin', '#DeltaR(Track, GEN Match)', (-1, 10), 1, 55],
    'track_genMatchDxyzmin': ['track_genMatchDxyzmin', '#DeltaXYZ(Track, GEN Match)', (-1, 10), 1, 55],
    'track_genMatchDrminold': ['track_genMatchDrminold', '#DeltaR(Track, GEN Match)', (-1, 10), 1, 55],

    'track_genMatchMotherTauDecay_alt': ['track_genMatchMotherTauDecay_alt', 'GEN Match Mother TauDecay', (-1, 25), 1, 26],

    'track_isDisaptrackcandidate': ['track_isDisaptrackcandidate', 'Track is Disap. Track Cand.', (0, 2), 1, 2],

    'track_hasAssociatedPV': ['track_hasAssociatedPV', 'Track has associated PV', (0, 2), 1, 2],

    'track_log10_distPVAssPVxy_': ['track_log10_distPVAssPVxy_', 'log_{10}(#Deltaxy(PV, associated PV) (cm))', 'log10(track_distPVAssPVxy)', (-5, 1), 1, 60],
    'track_log10_distPVAssPVz_': ['track_log10_distPVAssPVz_', 'log_{10}(#Deltaz(PV, associated PV) (cm))', 'log10(track_distPVAssPVz)', (-5, 1), 1, 60],

    'track_hasAssociatedSV': ['track_hasAssociatedSV', 'Track has associated SV', (0, 2), 1, 2],
    'track_distPVAssSVxy': ['track_distPVAssSVxy', '#Deltaxy(PV, associated SV) (cm)', (-0.1, 1), 1, 55],
    'track_distPVAssSVz': ['track_distPVAssSVz', '#Deltaz(PV, associated SV) (cm)', (-1, 11), 1, 60],
    'track_log10_distPVAssSVxy_': ['track_log10_distPVAssSVxy_', 'log_{10}(#Deltaxy(PV, associated SV) (cm))', 'log10(track_distPVAssSVxy)', (-5, 1), 1, 60],
    'track_log10_distPVAssSVz_': ['track_log10_distPVAssSVz_', 'log_{10}(#Deltaz(PV, associated SV) (cm))', 'log10(track_distPVAssSVz)', (-5, 1), 1, 60],

    'track_IPsigSV': ['track_IPsigSV', 'IP Significance w.r.t. SV', (0, 150), 1, 30],
    'track_IPxyzSV': ['track_IPxyzSV', 'IPxyz w.r.t. SV (cm)', (0, 2), 1, 40],
    'track_IPxySV': ['track_IPxySV', 'IPxy w.r.t. SV (cm)', (0, 2), 1, 40],
    'track_IPzSV': ['track_IPzSV', 'IPz w.r.t. SV (cm)', (0, 2), 1, 40],

    'track_log10_IPsigSV_': ['track_log10_IPsigSV_', 'log_{10}(IP Significance w.r.t. SV)', (-2, 4), 1, 30],
    'track_log10_IPxyzSV_': ['track_log10_IPxyzSV_', 'log_{10}(IPxyz w.r.t. SV (cm))', (-6, 2), 1, 40],
    'track_log10_IPxySV_': ['track_log10_IPxySV_', 'log_{10}(IPxy w.r.t. SV (cm))', (-6, 2), 1, 40],
    'track_log10_IPzSV_': ['track_log10_IPzSV_', 'log_{10}(IPz w.r.t. SV (cm))', (-6, 2), 1, 40],

    'track_IPsigXY': ['track_IPsigXY', 'IPxy Significance w.r.t. PV', (0, 150), 1, 30],
    'track_log10_IPsigXY_': ['track_log10_IPsigXY_', 'log_{10}(IPxy Significance w.r.t. PV)', (-2, 4), 1, 30],
    'track_IPsigZ': ['track_IPsigZ', 'IPz Significance w.r.t. PV', (0, 150), 1, 30],
    'track_log10_IPsigZ_': ['track_log10_IPsigZ_', 'log_{10}(IPz Significance w.r.t. PV)', (-2, 4), 1, 30],

    'track_IPsigXYPU': ['track_IPsigXYPU', 'IPxy Significance w.r.t. PU', (0, 150), 1, 30],
    'track_log10_IPsigXYPU_': ['track_log10_IPsigXYPU_', 'log_{10}(IPxy Significance w.r.t. PU)', (-2, 4), 1, 30],
    'track_IPsigZPU': ['track_IPsigZPU', 'IPz Significance w.r.t. PU', (0, 150), 1, 30],
    'track_log10_IPsigZPU_': ['track_log10_IPsigZPU_', 'log_{10}(IPz Significance w.r.t. PU)', (-2, 4), 1, 30],
    
    'track_IPsigXYAssPV': ['track_IPsigXYAssPV', 'IPxy Significance w.r.t. Ass. PV', (0, 150), 1, 30],
    'track_log10_IPsigXYAssPV_': ['track_log10_IPsigXYAssPV_', 'log_{10}(IPxy Significance w.r.t. Ass. PV)', (-2, 4), 1, 30],
    'track_IPsigZAssPV': ['track_IPsigZAssPV', 'IPz Significance w.r.t. Ass. PV', (0, 150), 1, 30],
    'track_log10_IPsigZAssPV_': ['track_log10_IPsigZAssPV_', 'log_{10}(IPz Significance w.r.t. Ass. PV)', (-2, 4), 1, 30],
    
    'track_IPsigXYPUAssPV': ['track_IPsigXYPUAssPV', 'IPxy Significance w.r.t. PU (excl. ass. PV)', (0, 150), 1, 30],
    'track_log10_IPsigXYPUAssPV_': ['track_log10_IPsigXYPUAssPV_', 'log_{10}(IPxy Significance w.r.t. PU (excl. ass. PV))', (-2, 4), 1, 30],
    'track_IPsigZPUAssPV': ['track_IPsigZPUAssPV', 'IPz Significance w.r.t. PU (excl. ass. PV)', (0, 150), 1, 30],
    'track_log10_IPsigZPUAssPV_': ['track_log10_IPsigZPUAssPV_', 'log_{10}(IPz Significance w.r.t. PU (excl. ass. PV))', (-2, 4), 1, 30],

    'cleaning_tracksRemoved': ['cleaning_tracksRemoved', 'DY Cleaning N Tracks Removed', (-1, 4), 1, 5],
    'cleaning_pfcandsRemoved': ['cleaning_pfcandsRemoved', 'DY Cleaning N PFCs Removed', (-1, 4), 1, 5],
    'cleaning_photonsRemoved': ['cleaning_photonsRemoved', 'DY Cleaning N Photons Removed', (-1, 4), 1, 5],
    'cleaning_jetsRemoved': ['cleaning_jetsRemoved', 'DY Cleaning N Jets Removed', (-1, 4), 1, 5],
    'cleaning_svsRemoved': ['cleaning_svsRemoved', 'DY Cleaning N SVs Removed', (-1, 4), 1, 5],
    
    'HLT_PFMET90_PFMHT90_IDTight_v': ['HLT_PFMET90_PFMHT90_IDTight_v', 'HLT_PFMET90_PFMHT90_IDTight_v', (0, 2), 1, 2],
    'HLT_PFMET100_PFMHT100_IDTight_v': ['HLT_PFMET100_PFMHT100_IDTight_v', 'HLT_PFMET100_PFMHT100_IDTight_v', (0, 2), 1, 2],
    'HLT_PFMET110_PFMHT110_IDTight_v': ['HLT_PFMET110_PFMHT110_IDTight_v', 'HLT_PFMET110_PFMHT110_IDTight_v', (0, 2), 1, 2],
    'HLT_PFMET120_PFMHT120_IDTight_v': ['HLT_PFMET120_PFMHT120_IDTight_v', 'HLT_PFMET120_PFMHT120_IDTight_v', (0, 2), 1, 2],
    'HLT_PFMET130_PFMHT130_IDTight_v': ['HLT_PFMET130_PFMHT130_IDTight_v', 'HLT_PFMET130_PFMHT130_IDTight_v', (0, 2), 1, 2],
    'HLT_PFMET140_PFMHT140_IDTight_v': ['HLT_PFMET140_PFMHT140_IDTight_v', 'HLT_PFMET140_PFMHT140_IDTight_v', (0, 2), 1, 2],
    'HLT_PFMET100_PFMHT100_IDTight_PFHT60_v': ['HLT_PFMET100_PFMHT100_IDTight_PFHT60_v', 'HLT_PFMET100_PFMHT100_IDTight_PFHT60_v', (0, 2), 1, 2],
    'HLT_PFMET110_PFMHT110_IDTight_PFHT60_v': ['HLT_PFMET110_PFMHT110_IDTight_PFHT60_v', 'HLT_PFMET110_PFMHT110_IDTight_PFHT60_v', (0, 2), 1, 2],
    'HLT_PFMET120_PFMHT120_IDTight_PFHT60_v': ['HLT_PFMET120_PFMHT120_IDTight_PFHT60_v', 'HLT_PFMET120_PFMHT120_IDTight_PFHT60_v', (0, 2), 1, 2],
    'HLT_PFMET130_PFMHT130_IDTight_PFHT60_v': ['HLT_PFMET130_PFMHT130_IDTight_PFHT60_v', 'HLT_PFMET130_PFMHT130_IDTight_PFHT60_v', (0, 2), 1, 2],
    'HLT_PFMET140_PFMHT140_IDTight_PFHT60_v': ['HLT_PFMET140_PFMHT140_IDTight_PFHT60_v', 'HLT_PFMET140_PFMHT140_IDTight_PFHT60_v', (0, 2), 1, 2],
    'HLT_PFMETNoMu90_PFMHTNoMu90_IDTight_v': ['HLT_PFMETNoMu90_PFMHTNoMu90_IDTight_v', 'HLT_PFMETNoMu90_PFMHTNoMu90_IDTight_v', (0, 2), 1, 2],
    'HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_v': ['HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_v', 'HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_v', (0, 2), 1, 2],
    'HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_v': ['HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_v', 'HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_v', (0, 2), 1, 2],
    'HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v': ['HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v', 'HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v', (0, 2), 1, 2],
    'HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_v': ['HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_v', 'HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_v', (0, 2), 1, 2],
    'HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_v': ['HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_v', 'HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_v', (0, 2), 1, 2],
    'HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60_v': ['HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60_v', 'HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60_v', (0, 2), 1, 2],
    'HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_PFHT60_v': ['HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_PFHT60_v', 'HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_PFHT60_v', (0, 2), 1, 2],
    'HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v': ['HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v', 'HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v', (0, 2), 1, 2],
    'HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_PFHT60_v': ['HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_PFHT60_v', 'HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_PFHT60_v', (0, 2), 1, 2],
    'HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_PFHT60_v': ['HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_PFHT60_v', 'HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_PFHT60_v', (0, 2), 1, 2],

    'htoverhtMiss': ['htoverhtMiss', 'H_{T} / H_{T}^{miss}', 'ht/htMiss', (0, 10), 1, 50],
    'ht5overhtMiss': ['ht5overhtMiss', 'H_{T}^{5} / H_{T}^{miss}', 'ht5/htMiss', (0, 10), 1, 50],

    # TODO: add PFCand related variables (energy/pT,...)

})

myvariables_V15_variations = myvariables_V15.copy()
myvariables_V15_variations.update({

    'met_ptJECup': ['met_ptJECup', 'p_{T}^{miss} (GeV)', (0, 2500), 1, 50],
    'met_phiJECup': ['met_phiJECup', '#phi(p_{T}^{miss})', (-4, 4), 1, 40],
    'htJECup': ['htJECup', 'H_{T} (GeV)', (0, 5000), 1, 50],
    'ht5JECup': ['ht5JECup', 'H_{T}^{5} (GeV)', (0, 5000), 1, 50],
    'htMissJECup': ['htMissJECup', 'H_{T}^{miss} (GeV)', (0, 2500), 1, 50],
    'n_jet_100JECup': ['n_jet_100JECup', 'N Jet (p_{T} > 100 GeV)', (0, 10), 1, 10],
    'n_jet_30_btagDeepCSVmediumJECup': ['n_jet_30_btagDeepCSVmediumJECup', 'N Jets b-Tag Medium WP (jets w/ p_{T} > 30 GeV)', (0, 10), 1, 10],
    'dphiminMetJetsJECup': ['dphiminMetJetsJECup', '#Delta#phi_{min}(p_{T}^{miss}, Jet_{1,2,3,4})', (0, 3.5), 1, 35],

    'met_ptJECdn': ['met_ptJECdn', 'p_{T}^{miss} (GeV)', (0, 2500), 1, 50],
    'met_phiJECdn': ['met_phiJECdn', '#phi(p_{T}^{miss})', (-4, 4), 1, 40],
    'htJECdn': ['htJECdn', 'H_{T} (GeV)', (0, 5000), 1, 50],
    'ht5JECdn': ['ht5JECdn', 'H_{T}^{5} (GeV)', (0, 5000), 1, 50],
    'htMissJECdn': ['htMissJECdn', 'H_{T}^{miss} (GeV)', (0, 2500), 1, 50],
    'n_jet_100JECdn': ['n_jet_100JECdn', 'N Jet (p_{T} > 100 GeV)', (0, 10), 1, 10],
    'n_jet_30_btagDeepCSVmediumJECdn': ['n_jet_30_btagDeepCSVmediumJECdn', 'N Jets b-Tag Medium WP (jets w/ p_{T} > 30 GeV)', (0, 10), 1, 10],
    'dphiminMetJetsJECdn': ['dphiminMetJetsJECdn', '#Delta#phi_{min}(p_{T}^{miss}, Jet_{1,2,3,4})', (0, 3.5), 1, 35],

    'met_ptJERup': ['met_ptJERup', 'p_{T}^{miss} (GeV)', (0, 2500), 1, 50],
    'met_phiJERup': ['met_phiJERup', '#phi(p_{T}^{miss})', (-4, 4), 1, 40],
    'htJERup': ['htJERup', 'H_{T} (GeV)', (0, 5000), 1, 50],
    'ht5JERup': ['ht5JERup', 'H_{T}^{5} (GeV)', (0, 5000), 1, 50],
    'htMissJERup': ['htMissJERup', 'H_{T}^{miss} (GeV)', (0, 2500), 1, 50],
    'n_jet_100JERup': ['n_jet_100JERup', 'N Jet (p_{T} > 100 GeV)', (0, 10), 1, 10],
    'n_jet_30_btagDeepCSVmediumJERup': ['n_jet_30_btagDeepCSVmediumJERup', 'N Jets b-Tag Medium WP (jets w/ p_{T} > 30 GeV)', (0, 10), 1, 10],
    'dphiminMetJetsJERup': ['dphiminMetJetsJERup', '#Delta#phi_{min}(p_{T}^{miss}, Jet_{1,2,3,4})', (0, 3.5), 1, 35],

    'met_ptJERdn': ['met_ptJERdn', 'p_{T}^{miss} (GeV)', (0, 2500), 1, 50],
    'met_phiJERdn': ['met_phiJERdn', '#phi(p_{T}^{miss})', (-4, 4), 1, 40],
    'htJERdn': ['htJERdn', 'H_{T} (GeV)', (0, 5000), 1, 50],
    'ht5JERdn': ['ht5JERdn', 'H_{T}^{5} (GeV)', (0, 5000), 1, 50],
    'htMissJERdn': ['htMissJERdn', 'H_{T}^{miss} (GeV)', (0, 2500), 1, 50],
    'n_jet_100JERdn': ['n_jet_100JERdn', 'N Jet (p_{T} > 100 GeV)', (0, 10), 1, 10],
    'n_jet_30_btagDeepCSVmediumJERdn': ['n_jet_30_btagDeepCSVmediumJERdn', 'N Jets b-Tag Medium WP (jets w/ p_{T} > 30 GeV)', (0, 10), 1, 10],
    'dphiminMetJetsJERdn': ['dphiminMetJetsJERdn', '#Delta#phi_{min}(p_{T}^{miss}, Jet_{1,2,3,4})', (0, 3.5), 1, 35],

})

myvariables_V15_leptoncr = myvariables_V15.copy()
myvariables_V15_leptoncr.update({
    'track_IPxy_widerange': ['track_IPxy_widerange', 'IPxy w.r.t. PV (cm)', 'track_IPxy', (0, 20), 1, 200],
    'track_dxy_widerange': ['track_dxy_widerange', 'dxy w.r.t. PV (cm)', 'track_dxy', (0, 20), 1, 200],
})

