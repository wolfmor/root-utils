"""Module to plot, create, and edit ROOT histograms.

To use ROOT.RDataFrame to process TreeSamples a rather recent ROOT version is required, e.g.,
source /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.24.06/x86_64-centos7-gcc48-opt/bin/thisroot.sh
"""

import os
import ctypes
import ROOT
from glob import glob
from array import array
from plotstyle import Plotstyle

ROOT.gEnv.SetValue('RooFit.Banner', 0)

# enable multi-threading to speed up processing of ROOT.RDataFrame (if supported by ROOT version)
try:
    ROOT.EnableImplicitMT()  # doesn't work with ROOT.RDataFrame().Range()
except AttributeError:
    pass


class PlotFactory:
    """Main class to make plots and/or create histograms.

    Methods
    -------
    process()
        Create plots and/or save histograms according to the added variables, samples, and ratios.
    add_variables(variablelist)
        Add a list of variables to be processed.
    add_variable(variable)
        Add a variable to be processed.
    add_samples(samplelist)
        Add a list of samples to be processed.
    add_sample(sample)
        Add a sample to be processed.
    add_ratios(ratiolist)
        Add a list of ratios to be processed.
    add_ratio(ratio)
        Add a ratio to be processed.
    """

    def __init__(self,
                 inputfiles=None,
                 inputpattern='VARIABLESAMPLE',
                 outputpath='.',
                 outputpattern='VARIABLE',
                 outputformat='pdf',
                 normalize=False,
                 axes='linlog',
                 ylabel='Events',
                 ylabelratio='Data / Prediction',
                 yaxisrangeratio=(0.0001, 1.9999),
                 ratiohlines=None,
                 linewidth=1,
                 markersize=2,
                 poslegend=(0.4, 0.55, 0.93, 0.9),
                 ncolumnslegend=1,
                 boldlegend=False,
                 text='36 fb^{-1} (13 TeV)',
                 extratext='Simulation Private Work',
                 height=1280,
                 width=None,
                 ipos=0,
                 uoflowbins=False):
        """
        Parameters
        ----------
        inputfiles : list[str], optional
            List of input ROOT files containing histograms.
        inputpattern : str, optional
            Naming scheme for histograms stored in the inputfiles,
            VARIABLE is replaced by the variable name, SAMPLE is replaced by the sample name.
        outputpath : str, default '.'
            Where to save the output plots and histograms, directory is created if it doesn't exist.
        outputpattern : str, default 'VARIABLE'
            Naming scheme for output plots, VARIABLE is replaced by the variable name.
        outputformat : str or list[str], default 'pdf'
            Format(s) to save output plots, has to be supported by ROOT.TCanvas.Print()
        normalize : bool, default False
            Scale all histograms to unity.
        axes : {'linlog', 'lin', 'log'}
            Draw plots with a linear y-axis, logarithmic y-axis, or both.
        ylabel : str, default 'Events'
            Label for y-axis of main plot (upper panel).
        ylabelratio : str, default 'Data / Prediction'
            Label for y-axis of lower panel.
        yaxisrangeratio : tuple[float], default (0.0001, 1.9999)
            Range of y-axis of lower panel, choose values just above/below a round number
            to avoid the respective tick label.
        ratiohlines : list[float or tuple[float]], default [1.]
            Straight line(s) to draw in lower panel, if element is float a horizontal line is drawn at the given value,
            if element is tuple of the form (x1, y1, x2, y2) a line is drawn from (x1,y1) to (x2,y2).
        linewidth : int, default 1
            Linewidth used for samples of category marker or line and for ratios,
            corresponding to ROOT.TAttLine.
        markersize : int, default 2
            Markersize used for samples of category marker,
            corresponding to ROOT.TAttMarker.
        poslegend : tuple[float], default (0.4, 0.55, 0.93, 0.9)
            Position of legend, specified by tuple (x1, y1, x2, y2).
        ncolumnslegend : int, default 1
            Number of columns used in legend.
        boldlegend : bool, default False
            Use bold font in legend.
        text : str, default '36 fb^{-1} (13 TeV)'
            Text to write on top right corner.
        extratext : str, default '#splitline{Work in progress}{Simulation}'
            Text to write next to CMS stamp.
        height : int, default 1280
            Height of plot in pixels.
        width : int or None, default None
            Width of plot in pixels, if None then width is chosen proportionally to height.
        ipos : int, default 0
            Position of CMS stamp + extratext, e.g., 0 for top left out of frame, 11 for top left inside frame.
            For details see CMS_lumi.py.
        uoflowbins : bool, default False
            Whether to plot the under- and overflow bins

        """

        if inputfiles is None:
            inputfiles = []

        if ratiohlines is None:
            ratiohlines = [1.]

        self.variables = []

        self.stacksamples = []
        self.markersamples = []
        self.linesamples = []

        self.ratios = []

        self.inputfiles = [ROOT.TFile(inputfile) for inputfile in inputfiles]
        self.inputpattern = inputpattern

        if outputpath[-1] == '/':
            self.outputpath = outputpath[:-1]
        else:
            self.outputpath = outputpath

        if not os.path.exists(outputpath):
            os.makedirs(outputpath)

        self.outputpattern = outputpattern
        self.outputformat = outputformat

        self.normalize = normalize
        if axes not in ['linlog', 'lin', 'log']:
            raise NotImplementedError('axes must be one of: lin, log, linlog')
        self.axes = axes

        self.ylabel = ylabel
        self.ylabelratio = ylabelratio

        self.yaxisrangeratio = yaxisrangeratio
        self.ratiohlines = ratiohlines

        self.linewidth = linewidth
        self.markersize = markersize

        self.p = None
        self.text = text
        self.extratext = extratext
        self.height = height
        self.width = width
        self.ipos = ipos

        self.legend = None
        self.poslegend = poslegend
        self.ncolumnslegend = ncolumnslegend
        self.boldlegend = boldlegend
        self.uoflowbins = uoflowbins

        self.histos = {}
        self.stacks = {}
        self.sums = {}
        self.ratiohistos = {}
        self.emptylinhistos = {}
        self.emptyloghistos = {}

    def process(self):

        if self.width is None:
            if len(self.ratios) > 0:
                self.width = int(1. * self.height)
            else:
                self.width = int(1.3 * self.height)

        self.p = Plotstyle(text=self.text, extratext=self.extratext, H_ref=self.height, W_ref=self.width, iPos=self.ipos)
        mystyle = self.p.setStyle()
        mystyle.cd()

        self.legend = ROOT.TLegend(self.poslegend[0], self.poslegend[1], self.poslegend[2], self.poslegend[3])
        self.legend.SetNColumns(self.ncolumnslegend)
        self.legend.SetFillStyle(0)
        self.legend.SetBorderSize(0)
        if self.boldlegend: self.legend.SetTextFont(62)
        else: self.legend.SetTextFont(42)

        self._get_histos()
        self._make_ratios()
        self._style_histos()
        self._draw_plots()
        self._save_histos()

    def add_variables(self, variablelist):
        """

        Parameters
        ----------
        variablelist :
        """
        for variable in variablelist:
            self.add_variable(variable)

    def add_variable(self, variable):

        if str(variable) in [str(v) for v in self.variables]:
            raise AssertionError(str(variable) + ' is a duplicate (names must be unique)')

        self.variables.append(variable)

    def add_samples(self, samplelist):
        for sample in samplelist:
            self.add_sample(sample)

    def add_sample(self, sample):

        if str(sample) in [str(s) for s in self.stacksamples + self.markersamples + self.linesamples]:
            raise AssertionError(str(sample) + ' is a duplicate (names must be unique)')
        
        if sample.category == 'stack':
            self.stacksamples.append(sample)
        elif sample.category == 'marker':
            self.markersamples.append(sample)
        elif sample.category == 'line':
            self.linesamples.append(sample)
        else:
            raise NotImplementedError('unknown sample category')

    def add_ratios(self, ratiolist):
        for ratio in ratiolist:
            self.add_ratio(ratio)

    def add_ratio(self, ratio):
        self.ratios.append(ratio)

    def _get_histos(self):

        print('\n# Get Histos')

        makeflat = {}
        for s in self.stacksamples + self.markersamples + self.linesamples:

            print(s)

            if isinstance(s, TreeSample):

                # first book all histograms
                for v in self.variables:

                    if not s.modifyvarname(v.vartoplot) == s.modifyvarname(v.name) and s.modifyvarname(v.name) not in s.df.GetColumnNames():
                        s.df = s.df.Define(
                            s.modifyvarname(v.name), s.modifyvarname(v.vartoplot)
                        )

                    if type(v.nbins) == str and 'flat' in v.nbins:
                        makeflat[v] = v.nbins
                        v.nbins = (v.axisrange[1] - v.axisrange[0]) * 100  # TODO: might be problematic for large ranges

                    if s.vectorselection is None:
                        if s.weight is None:
                            self.histos[v + s] = s.df.Histo1D(
                                ('h' + str(v + s), '', v.nbins, v.axisrange[0], v.axisrange[1]), s.modifyvarname(v.name)
                            )
                        else:
                            self.histos[v + s] = s.df.Define(
                                '_w', s.weight
                            ).Histo1D(
                                ('h' + str(v + s), '', v.nbins, v.axisrange[0], v.axisrange[1]), s.modifyvarname(v.name), '_w'
                            )
                    else:
                        if s.weight is None:
                            self.histos[v + s] = s.df.Define(
                                s.modifyvarname(v.name) + '_pass', s.modifyvarname(v.name) + '[' + s.vectorselection + ']'
                            ).Histo1D(
                                ('h' + str(v + s), '', v.nbins, v.axisrange[0], v.axisrange[1]), s.modifyvarname(v.name) + '_pass'
                            )
                        else:
                            self.histos[v + s] = s.df.Define(
                                s.modifyvarname(v.name) + '_pass', s.modifyvarname(v.name) + '[' + s.vectorselection + ']'
                            ).Define(
                                '_w', s.weight
                            ).Histo1D(
                                ('h' + str(v + s), '', v.nbins, v.axisrange[0], v.axisrange[1]), s.modifyvarname(v.name) + '_pass', '_w'
                            )

                # then trigger the filling of the histos
                s.df.Report().Print()

                # get the pointer to the real histogram instead of the RResultPtr
                for v in self.variables:
                    self.histos[v + s] = self.histos[v + s].GetPtr()
                    self.histos[v + s].UseCurrentStyle()

            else:

                for v in self.variables:

                    if s.file is None:
                        self.histos[v + s] = self.inputfiles[s.inputfileindex].Get(
                            self.inputpattern.replace('VARIABLE', v.vartoplot).replace('SAMPLE', s.name)
                        )
                    else:
                        self.histos[v + s] = s.file.Get(
                            self.inputpattern.replace('VARIABLE', v.vartoplot).replace('SAMPLE', s.name)
                        )

                    if type(v.rebin) == list:
                        self.histos[v + s] = self.histos[v + s].Rebin(
                            len(v.rebin)-1, self.histos[v + s].GetName() + 'rebinned', array('d', v.rebin)
                        )
                    else:
                        self.histos[v + s].Rebin(v.rebin)

                    self.histos[v + s].GetXaxis().SetRangeUser(v.axisrange[0], v.axisrange[1])
                    self.histos[v + s].UseCurrentStyle()


        for iv, v in enumerate(self.variables):
            
            if v in makeflat:
                
                flatsample = makeflat[v].split(':')[1]
                flatnbins = float(makeflat[v].split(':')[2])

                if flatnbins < 1:  # interpret as relative stat error in each bin
                    flatnperbin = 1. / flatnbins**2.
                    flatnbins = int(self.histos[v + flatsample].Integral() / flatnperbin)
                else:
                    flatnbins = int(flatnbins)

                p = array('d', [i/float(flatnbins) for i in range(1, flatnbins)])
                quantiles = array('d', (flatnbins-1)*[0.])
                self.histos[v + flatsample].GetQuantiles(flatnbins-1, quantiles, p)

                # in order to make sure to use existing bin edges
                roundedquantiles = []
                for q in quantiles:
                    roundedquantiles.append(round(q, 2))
                quantiles = array('d', roundedquantiles)

                quantiles.insert(0, v.axisrange[0])
                quantiles.append(v.axisrange[1])

                for s in self.stacksamples + self.markersamples + self.linesamples:
                    self.histos[v + s] = self.histos[v + s].Rebin(
                        flatnbins, self.histos[v + s].GetName() + 'rebinned', quantiles
                    )

            if v.blind is not None:
                for blindthis in v.blind:
                    if not len(blindthis.split(':')) == 3:
                        raise NotImplementedError('cannot interpret blind: ' + blindthis)

                    blindsample, blindstart, blindend = blindthis.split(':')

                    blindstart = blindstart.replace('START', str(v.axisrange[0])).replace('END', str(v.axisrange[1]))
                    blindend = blindend.replace('START', str(v.axisrange[0])).replace('END', str(v.axisrange[1]))

                    blindstartbin = self.histos[v + blindsample].FindBin(float(blindstart))
                    blindendbin = self.histos[v + blindsample].FindBin(float(blindend))

                    for blindbin in range(blindstartbin, blindendbin+1):

                        self.histos[v + blindsample].SetBinContent(blindbin, 0.)

                        # also blind the under- and overflow bins
                        if blindbin == 1:
                            self.histos[v + blindsample].SetBinContent(0, 0.)
                        if blindbin == self.histos[v + blindsample].GetNbinsX():
                            self.histos[v + blindsample].SetBinContent(self.histos[v + blindsample].GetNbinsX() + 1, 0.)

            for s in self.stacksamples + self.markersamples + self.linesamples:
                if not self.histos[v + s].GetSumw2N():
                    self.histos[v + s].Sumw2()

                if self.uoflowbins:
                    self.histos[v + s].GetXaxis().SetRange(0, self.histos[v + s].GetNbinsX() + 1)


            if len(self.stacksamples) > 0:
                self.sums[v] = self.histos[v + self.stacksamples[0]].Clone('sum' + v.name)
                self.sums[v].Reset('ICESM')
                for s in self.stacksamples:
                    self.sums[v].Add(self.histos[v + s])

            if self.normalize:
                for s in self.markersamples + self.linesamples:
                    if self.histos[v + s].Integral() > 0:
                        self.histos[v + s].Scale(1. / self.histos[v + s].Integral())
                if len(self.stacksamples) > 0 and self.sums[v].Integral() > 0:
                    for s in self.stacksamples:
                        self.histos[v + s].Scale(1. / self.sums[v].Integral())
                    self.sums[v].Scale(1. / self.sums[v].Integral())

            for s in self.markersamples + self.linesamples:
                if s.scaleto is None: continue
                scaletoargs = s.scaleto.split(':')
                if len(scaletoargs) == 3 and scaletoargs[0] in [str(_s) for _s in self.stacksamples + self.markersamples + self.linesamples] or scaletoargs[0] == 'STACK':

                    scaletostart = float(scaletoargs[1].replace('START', str(v.axisrange[0])))
                    scaletoend = float(scaletoargs[2].replace('END', str(v.axisrange[1])))

                    scaletostartbin = self.histos[v + s].GetXaxis().FindBin(scaletostart)
                    scaletoendbin = self.histos[v + s].GetXaxis().FindBin(scaletoend)

                    if scaletoargs[0] == 'STACK':
                        scaletotarget = self.sums[v].Integral(scaletostartbin, scaletoendbin)
                    else:
                        scaletotarget = self.histos[str(v) + scaletoargs[0]].Integral(scaletostartbin, scaletoendbin)

                    if self.histos[v + s].Integral() > 0:
                        self.histos[v + s].Scale(scaletotarget / self.histos[v + s].Integral())

                elif len(scaletoargs) == 1 and scaletoargs[0].replace('.', '', 1).isdigit():

                    if self.histos[v + s].Integral() > 0:
                        self.histos[v + s].Scale(float(scaletoargs[0]) / self.histos[v + s].Integral())

                else:
                    raise NotImplementedError('cannot interpret scaleto: ' + s.scaleto)


            if len(self.stacksamples) > 0:
                self.stacks[v] = ROOT.THStack('stack' + v.name, '')
                for s in self.stacksamples:
                    self.stacks[v].Add(self.histos[v + s])

            if iv == 0:
                for s in self.markersamples:
                    if len(s.title) > 0:
                        self.legend.AddEntry(self.histos[v + s], s.title, 'pe')

                for s in self.stacksamples[::-1]:
                    if len(s.title) > 0:
                        self.legend.AddEntry(self.histos[v + s], s.title, 'f')

                for s in self.linesamples:
                    if len(s.title) > 0:
                        self.legend.AddEntry(self.histos[v + s], s.title, 'l')

    def _make_ratios(self):

        print('\n# Make Ratios')

        for v in self.variables:

            print(v)

            for r in self.ratios:

                if r.category == 'ratio':

                    if not len(r.name.split(':')) == 2:
                        raise NotImplementedError('cannot interpret ratio')

                    numerator = r.name.split(':')[0]
                    denominator = r.name.split(':')[1]

                    if numerator == 'STACK':
                        self.ratiohistos[v + r] = self.sums[v].Clone(v + r)
                    else:
                        self.ratiohistos[v + r] = self.histos[v.name + numerator].Clone(v + r)

                    if denominator == 'STACK':
                        self.ratiohistos[v + r].Divide(self.sums[v])
                    else:
                        self.ratiohistos[v + r].Divide(self.histos[v.name + denominator])

                elif r.category in ['cutsig', 'binsig']:

                    if not len(r.name.split(':')) == 3:
                        raise NotImplementedError('cannot interpret ratio')

                    signal = r.name.split(':')[0]
                    background = r.name.split(':')[1]
                    syserrB = float(r.name.split(':')[2])

                    for s in self.stacksamples + self.markersamples + self.linesamples:
                        self.ratiohistos[v + r] = self.histos[v + s].Clone(v + r)
                        self.ratiohistos[v + r].Reset()
                        break

                    nbins = self.ratiohistos[v + r].GetNbinsX()
                    for b in range(nbins):

                        b += 1

                        if r.category == 'binsig': nbins = b

                        staterrS = ctypes.c_double(0.)
                        staterrB = ctypes.c_double(0.)

                        if signal == 'STACK':
                            S = self.sums[v].IntegralAndError(b, nbins, staterrS)
                        else:
                            S = self.histos[v.name + signal].IntegralAndError(b, nbins, staterrS)

                        if background == 'STACK':
                            B = self.sums[v].IntegralAndError(b, nbins, staterrB)
                        else:
                            B = self.histos[v.name + background].IntegralAndError(b, nbins, staterrB)

                        # TODO: implement error like this?
                        if S > 0 and B > 0:

                            sig = ROOT.RooStats.AsimovSignificance(S, B, syserrB*B)

                            sigUp = ROOT.RooStats.AsimovSignificance(S+ROOT.TMath.Sqrt(S), B-ROOT.TMath.Sqrt(B), syserrB*(B-ROOT.TMath.Sqrt(B)))
                            sigDown = ROOT.RooStats.AsimovSignificance(S-ROOT.TMath.Sqrt(S), B+ROOT.TMath.Sqrt(B), syserrB*(B+ROOT.TMath.Sqrt(B)))

                            self.ratiohistos[v + r].SetBinContent(b, sig)
                            self.ratiohistos[v + r].SetBinError(b, 0.5*(sigUp-sigDown))

                elif r.category == 'calibration':

                    if not len(r.name.split(':')) == 2:
                        raise NotImplementedError('cannot interpret calibration')

                    signal = r.name.split(':')[0]
                    background = r.name.split(':')[1]

                    if signal == 'STACK':
                        signalhisto = self.sums[v].Clone('calibration_signalhisto')
                    else:
                        signalhisto = self.histos[v.name + signal].Clone('calibration_signalhisto')

                    if background == 'STACK':
                        backgroundhisto = self.sums[v].Clone('calibration_backgroundhisto')
                    else:
                        backgroundhisto = self.histos[v.name + background].Clone('calibration_backgroundhisto')

                    if signalhisto.Integral() == 0 or backgroundhisto.Integral() == 0:

                        self.ratiohistos[v + r] = None
                        print('skipping calibration histogram because integral is zero')

                    else:

                        signalhisto.Scale(1. / signalhisto.Integral())
                        backgroundhisto.Scale(1. / backgroundhisto.Integral())

                        backgroundhisto.Add(signalhisto)

                        self.ratiohistos[v + r] = signalhisto.Clone(v + r)
                        self.ratiohistos[v + r].Divide(backgroundhisto)

                elif r.category == 'efficiency':

                    if not len(r.name.split(':')) == 2:
                        raise NotImplementedError('cannot interpret efficiency')

                    name_pass = r.name.split(':')[0]
                    name_total = r.name.split(':')[1]

                    if name_pass == 'STACK':
                        histo_pass = self.sums[v].Clone('efficiency_passhisto')
                    else:
                        histo_pass = self.histos[v.name + name_pass].Clone('efficiency_passhisto')

                    if name_total == 'STACK':
                        histo_total = self.sums[v].Clone('efficiency_totalhisto')
                    else:
                        histo_total = self.histos[v.name + name_total].Clone('efficiency_totalhisto')

                    if ROOT.TEfficiency.CheckConsistency(histo_pass, histo_total):

                        self.ratiohistos[v + r + 'eff'] = ROOT.TEfficiency(histo_pass, histo_total)
                        self.ratiohistos[v + r + 'eff'].SetName(r.name + 'eff')

                        self.ratiohistos[v + r] = histo_pass.Clone(r.name + 'empty')
                        self.ratiohistos[v + r].Reset()
                    else:
                        self.ratiohistos[v + r] = None
                        print('skipping efficiency because histograms are not compatible')

                else:
                    raise NotImplementedError('unknown ratio type')

    def _style_histos(self):

        print('\n# Style Histos')

        for v in self.variables:

            print(v)

            for s in self.stacksamples + self.markersamples + self.linesamples:
                self.emptylinhistos[v] = self.histos[v + s].Clone('emptylin' + v.name)
                self.emptyloghistos[v] = self.histos[v + s].Clone('emptylog' + v.name)
                break

            for emptyhisto in [self.emptylinhistos[v], self.emptyloghistos[v]]:
                emptyhisto.Reset()
                emptyhisto.GetYaxis().SetTitle(self.ylabel)
                if len(self.ratios) > 0:
                    emptyhisto.GetXaxis().SetLabelSize(0)
                else:
                    emptyhisto.GetXaxis().SetTitle(v.title)

            globalmin = min([self.sums[v].GetMinimum(0) if len(self.stacksamples) > 0 else float('inf')] +
                            [self.histos[v + s].GetMinimum(0) for s in self.markersamples + self.linesamples])

            globalmax = max([self.sums[v].GetMaximum() if len(self.stacksamples) > 0 else 0] +
                            [self.histos[v + s].GetMaximum() for s in self.markersamples + self.linesamples])

            if v.ymaxlog is None: logmax = globalmax
            else: logmax = v.ymaxlog
            if v.yminlog is None: logmin = globalmin
            else: logmin = v.yminlog
            logrange = ROOT.TMath.Log10(logmax) - ROOT.TMath.Log10(logmin)

            if v.yminlin is None:
                self.emptylinhistos[v].SetMinimum(0.)
            else:
                self.emptylinhistos[v].SetMinimum(v.yminlin)

            if v.ymaxlin is None:
                self.emptylinhistos[v].SetMaximum(2. * globalmax)
            else:
                self.emptylinhistos[v].SetMaximum(v.ymaxlin)

            if v.yminlog is None:
                self.emptyloghistos[v].SetMinimum(0.5 * globalmin)
            else:
                self.emptyloghistos[v].SetMinimum(v.yminlog)

            if v.ymaxlog is None:
                self.emptyloghistos[v].SetMaximum(globalmax * 10 ** max(1, logrange))
            else:
                self.emptyloghistos[v].SetMaximum(v.ymaxlog)

            for s in self.stacksamples:
                self.histos[v + s].SetLineWidth(0)
                self.histos[v + s].SetLineColor(s.color)

                self.histos[v + s].SetMarkerSize(0)

                self.histos[v + s].SetFillStyle(s.fillstyle)
                self.histos[v + s].SetFillColor(s.color)

            for s in self.markersamples:
                self.histos[v + s].SetLineWidth(self.linewidth)
                self.histos[v + s].SetLineColor(s.color)

                self.histos[v + s].SetMarkerSize(self.markersize)
                self.histos[v + s].SetMarkerColor(s.color)

            for s in self.linesamples:
                self.histos[v + s].SetLineWidth(self.linewidth)
                self.histos[v + s].SetLineStyle(s.linestyle)
                self.histos[v + s].SetLineColor(s.color)

                self.histos[v + s].SetMarkerSize(0)

            isfirstratiohisto = True
            for r in self.ratios:

                if self.ratiohistos[v + r] is None: continue

                self.ratiohistos[v + r].GetXaxis().SetTitle(v.title)
                self.ratiohistos[v + r].GetYaxis().SetTitle(self.ylabelratio)

                self.ratiohistos[v + r].SetLineWidth(self.linewidth)
                self.ratiohistos[v + r].SetLineStyle(self.histos[v.name + r.name.split(':')[0]].GetLineStyle())
                self.ratiohistos[v + r].SetLineColor(self.histos[v.name + r.name.split(':')[0]].GetLineColor())

                self.ratiohistos[v + r].SetMarkerSize(self.histos[v.name + r.name.split(':')[0]].GetMarkerSize())
                self.ratiohistos[v + r].SetMarkerColor(self.histos[v.name + r.name.split(':')[0]].GetMarkerColor())

                if r.category == 'efficiency':

                    self.ratiohistos[str(v + r) + 'eff'].SetLineWidth(self.linewidth)
                    self.ratiohistos[str(v + r) + 'eff'].SetLineStyle(self.histos[v.name + r.name.split(':')[0]].GetLineStyle())
                    self.ratiohistos[str(v + r) + 'eff'].SetLineColor(self.histos[v.name + r.name.split(':')[0]].GetLineColor())

                    self.ratiohistos[str(v + r) + 'eff'].SetMarkerSize(self.histos[v.name + r.name.split(':')[0]].GetMarkerSize())
                    self.ratiohistos[str(v + r) + 'eff'].SetMarkerColor(self.histos[v.name + r.name.split(':')[0]].GetMarkerColor())


                if isfirstratiohisto:

                    self.ratiohistos[v + r].SetMinimum(self.yaxisrangeratio[0])
                    self.ratiohistos[v + r].SetMaximum(self.yaxisrangeratio[1])

                    isfirstratiohisto = False

    def _draw_plots(self):

        print('\n# Draw Plots')

        for v in self.variables:

            print(v)

            xlow = v.axisrange[0]
            xhigh = v.axisrange[1]
            if self.uoflowbins:
                xlow = self.emptylinhistos[v].GetXaxis().GetBinLowEdge(0)
                xhigh = self.emptylinhistos[v].GetXaxis().GetBinUpEdge(self.emptylinhistos[v].GetNbinsX() + 1)

            lines = []
            for ratiohline in self.ratiohlines:

                if type(ratiohline) == tuple:

                    if not len(ratiohline) == 4:
                        raise NotImplementedError('cannot interpret ratiohline')

                    ratiohline = [float(rhl.replace('START', str(xlow)).replace('END', str(xhigh))) if type(rhl) == str else rhl for rhl in ratiohline]

                    lines.append(ROOT.TLine(ratiohline[0], ratiohline[1], ratiohline[2], ratiohline[3]))
                    lines[-1].SetLineWidth(1)
                    lines[-1].SetLineColor(ROOT.kBlack)

                else:

                    lines.append(ROOT.TLine(xlow, ratiohline, xhigh, ratiohline))
                    lines[-1].SetLineWidth(1)
                    lines[-1].SetLineColor(ROOT.kBlack)

            if self.axes == 'lin' or self.axes == 'log':

                canvas = ROOT.TCanvas('c' + v.name, 'c' + v.name, self.width, self.height)
                canvas.Divide(1, 1)

                if len(self.ratios) > 0:

                    subpads = self.p.addLowerPads(canvas)

                    subpads['1_1'].cd()

                    if self.axes == 'lin':
                        self.emptylinhistos[v].Draw('axis')
                        self._draw_histos(v, log=False)
                    else:
                        self.emptyloghistos[v].Draw('axis')
                        self._draw_histos(v, log=True)

                    subpads['1_2'].cd()

                    for r in self.ratios:

                        if self.ratiohistos[v + r] is None: continue

                        if r.name.split(':')[0] in [s.name for s in self.markersamples]:
                            self.ratiohistos[v + r].Draw('e x0 same')
                            if r.category == 'efficiency':
                                self.ratiohistos[str(v + r) + 'eff'].Draw('same')
                        # TODO: errors on ratio histos?
                        # elif r.name.split(':')[0] in [s.name for s in self.linesamples]:
                        #     self.ratiohistos[v + r].Draw('hist same')
                        #     self.ratiohistos[v + r].Draw('e x0 same')
                        else:
                            self.ratiohistos[v + r].Draw('hist same')
                            if r.category == 'efficiency':
                                self.ratiohistos[str(v + r) + 'eff'].Draw('same')

                        self.p.adjustLowerHisto(self.ratiohistos[v + r])

                    for line in lines:
                        line.Draw('same')

                else:

                    self.p.preparePad()

                    if self.axes == 'lin':
                        self.emptylinhistos[v].Draw('axis')
                        self._draw_histos(v, log=False)
                    else:
                        self.emptyloghistos[v].Draw('axis')
                        self._draw_histos(v, log=True)

                    self.p.postparePad()

            else:

                canvas = ROOT.TCanvas('c' + v.name, 'c' + v.name, 2 * self.width, self.height)
                canvas.Divide(2, 1)

                if len(self.ratios) > 0:

                    subpads = self.p.addLowerPads(canvas)

                    subpads['1_1'].cd()

                    self.emptylinhistos[v].Draw()
                    self._draw_histos(v, log=False)

                    subpads['2_1'].cd()

                    self.emptyloghistos[v].Draw()
                    self._draw_histos(v, log=True)

                    for lowerpad in ['1_2', '2_2']:

                        subpads[lowerpad].cd()
                        # TODO: draw small arrows if ratio out of range
                        for r in self.ratios:

                            if self.ratiohistos[v + r] is None: continue

                            if r.name.split(':')[0] in [s.name for s in self.markersamples]:
                                self.ratiohistos[v + r].Draw('e x0 same')
                                if r.category == 'efficiency':
                                    self.ratiohistos[str(v + r) + 'eff'].Draw('same')
                            # TODO: errors on ratio histos?
                            # elif r.name.split(':')[0] in [s.name for s in self.linesamples]:
                            #     self.ratiohistos[v + r].Draw('hist same')
                            #     self.ratiohistos[v + r].Draw('e x0 same')
                            else:
                                self.ratiohistos[v + r].Draw('hist same')
                                if r.category == 'efficiency':
                                    self.ratiohistos[str(v + r) + 'eff'].Draw('same')

                            self.p.adjustLowerHisto(self.ratiohistos[v + r])

                        for line in lines:
                            line.Draw('same')

                else:

                    canvas.cd(1)

                    self.p.preparePad()

                    self.emptylinhistos[v].Draw('axis')
                    self._draw_histos(v, log=False)

                    self.p.postparePad()

                    canvas.cd(2)

                    self.p.preparePad()

                    self.emptyloghistos[v].Draw('axis')
                    self._draw_histos(v, log=True)

                    self.p.postparePad()

            canvas.Update()
            canvas.Draw()
            if type(self.outputformat) == list:
                for outputfmt in self.outputformat:
                    canvas.Print(self.outputpath
                                 + '/' + self.outputpattern.replace('VARIABLE', v.name)
                                 + '.' + outputfmt)
            else:
                canvas.Print(self.outputpath
                             + '/' + self.outputpattern.replace('VARIABLE', v.name)
                             + '.' + self.outputformat)

    def _draw_histos(self, v, log):

        ROOT.gPad.Update()

        if len(self.stacksamples) > 0:
            self.stacks[v].Draw('hist same noclear')

        for s in self.linesamples:
            self.histos[v + s].Draw('hist same')

        for s in self.markersamples:
            self.histos[v + s].Draw('e x0 same')

        self.legend.Draw('same')

        if log: ROOT.gPad.SetLogy()

        self.p.postparePad()

    def _save_histos(self):


        if len([obj for obj in self.ratios + self.variables if obj.savehistos]) > 0:

            print('\n# Save Histos')

            fout = ROOT.TFile(self.outputpath + '/histos.root', 'recreate')
            dirout = {}

            for r in self.ratios:
                if r.savehistos and type(r.variablestosave) == list:
                    for v in r.variablestosave:
                        if self.ratiohistos[v + str(r)] is None: continue
                        self.ratiohistos[v + str(r)].Write()

            for v in self.variables:
                if v.savehistos:
                    if v.samplestosave == 'ALL':

                        for s in self.stacksamples + self.markersamples + self.linesamples:

                            if s.name not in dirout:
                                dirout[s.name] = fout.mkdir(s.name)
                            dirout[s.name].cd()

                            self.histos[v + s].Write()

                        if len(self.stacksamples) > 0:

                            if 'STACK' not in dirout:
                                dirout['STACK'] = fout.mkdir('STACK')
                            dirout['STACK'].cd()

                            self.sums[v].Write()

                    elif type(v.samplestosave) == list:

                        for s in v.samplestosave:

                            if s not in dirout:
                                dirout[s] = fout.mkdir(s)
                            dirout[s].cd()

                            if s == 'STACK':
                                self.sums[v].Write()
                            else:
                                self.histos[v + str(s)].Write()

            print('just created ' + str(fout.GetName()))


class Sample:
    def __init__(self, category, name,

                 title=None,

                 color=ROOT.kBlack,
                 linestyle=ROOT.kSolid,
                 fillstyle=1001,

                 scaleto=None):

        if category not in ['stack', 'marker', 'line']:
            raise NotImplementedError('unknown sample category')
        self.category = category

        self.name = name
        if title is None:
            self.title = name
        else:
            self.title = title

        self.color = color
        self.linestyle = linestyle
        self.fillstyle = fillstyle

        self.scaleto = scaleto

    def __str__(self):
        return self.name

    def __add__(self, other):
        return str(self) + str(other)


class HistoSample(Sample):
    def __init__(self, category, name,

                 title=None,

                 inputfileindex=0,
                 histofile=None,

                 color=ROOT.kBlack,
                 linestyle=ROOT.kSolid,
                 fillstyle=1001,

                 scaleto=None):

        Sample.__init__(self, category, name, title, color, linestyle, fillstyle, scaleto)

        self.inputfileindex = inputfileindex
        if histofile is None:
            self.file = None
        else:
            self.file = ROOT.TFile(histofile)


class TreeSample(Sample):
    def __init__(self, category, name, tree, files,

                 title=None,

                 eventselection=None,
                 vectorselection=None,
                 weight=None,
                 modifyvarname=lambda varname: varname,

                 color=ROOT.kBlack,
                 linestyle=ROOT.kSolid,
                 fillstyle=1001,

                 scaleto=None,

                 ntestfiles=0
                 ):

        print('Initializing ' + name)

        Sample.__init__(self, category, name, title, color, linestyle, fillstyle, scaleto)

        if ntestfiles:
            if type(files) == list:
                files = glob(files[0])[:ntestfiles]
            else:
                files = glob(files)[:ntestfiles]

        if eventselection is None or eventselection == '' or eventselection == '1':
            self.df = ROOT.RDataFrame(tree, files)
        else:
            self.df = ROOT.RDataFrame(tree, files).Filter(eventselection, ' selection for ' + self.name)

        if vectorselection is None or vectorselection == '' or vectorselection == '1':
            self.vectorselection = None
        else:
            self.vectorselection = vectorselection


        if weight is not None:

            weighttosplit = ''
            oktoreplace = True
            for ichar, character in enumerate(weight):
                if character in ['/', '*'] and oktoreplace:
                    weighttosplit += 'SPLITHERE'
                else:
                    weighttosplit += character

                if ichar > 1 and weight[ichar-2:ichar+1] == ':=[':
                    oktoreplace = False

                if not oktoreplace and character == ']':
                    oktoreplace = True

            weights = weighttosplit.split('SPLITHERE')
            for w in weights:
                if w in self.df.GetColumnNames() or w in ['XSEC', 'NSIM', 'COUNTER'] or w.replace('.', '', 1).isdigit(): continue
                wargs = w.split(':=[')
                if len(wargs) == 2:
                    self.df = self.df.Define(wargs[0], wargs[1].replace(']', ''))
                    weight = weight.replace(':=[' + wargs[1], '')
                else:
                    raise NotImplementedError('cannot interpret weight: ', w)


            self.smallchain = None
            if 'XSEC' in weight or 'NSIM' in weight:
                print(' get small chain')
                self.smallchain = ROOT.TChain(tree)
                if type(files) == list:
                    self.smallchain.Add(files[0], 1)
                else:
                    self.smallchain.Add(files, 1)

                if self.smallchain.GetEntries() > 0:

                    self.smallchain.GetEntry(0)

                    if 'XSEC' in weight:
                        print(' get XSEC')
                        if self.smallchain is not None:
                            weight = weight.replace('XSEC', str(self.smallchain.crossSection))
                        else:
                            weight = weight.replace('XSEC', '1.')

                    if 'NSIM' in weight:
                        print(' get NSIM')
                        if self.smallchain is not None:
                            weight = weight.replace('NSIM', str(self.smallchain.numSimEvents))
                        else:
                            weight = weight.replace('NSIM', '1.')
                else:
                    print(' small chain is empty')

            if 'COUNTER' in weight:
                print(' get COUNTER')
                self.counter = ROOT.RDataFrame('tCounter', files)
                weight = weight.replace('COUNTER',
                                        str(round(self.counter.Count().GetValue(), 1)))

            print(' weight ' + weight)
        self.weight = weight

        self.modifyvarname = modifyvarname

        # TODO: implement usage of TreeSamples without RDataFrame with:
        # deactivation of unused branches via setbranchstatus
        # entrylist used correctly?
        #
        # self.tree = tree
        # self.files = files
        #
        # self.selection = selection
        # print('# Prepare Sample ' + self.name)
        #
        # self.chain = ROOT.TChain()
        # for t in self.trees:
        #     for f in self.files:
        #         self.chain.Add(f + '/' + t)
        #
        # self.chain.Draw('>>entrylist' + self.name, selection, 'entrylist')
        # self.entrylist = ROOT.gDirectory.Get('entrylist' + self.name)
        #
        # self.chain.SetEntryList(self.entrylist)


class Variable:
    def __init__(self, name,
                 title=None,
                 vartoplot=None,
                 blind=None,
                 axisrange=(0, 1),
                 rebin=1,
                 nbins=100,
                 yminlin=None, ymaxlin=None, yminlog=None, ymaxlog=None,
                 savehistos=False,
                 samplestosave='ALL'):

        self.name = name
        if title is None:
            self.title = name
        else:
            self.title = title

        if vartoplot is None:
            self.vartoplot = name
        else:
            self.vartoplot = vartoplot

        if blind is None or blind == '':
            self.blind = None
        else:
            if type(blind) == list:
                self.blind = blind
            else:
                self.blind = [blind]


        self.axisrange = axisrange
        self.rebin = rebin
        self.nbins = nbins

        self.yminlin = yminlin
        self.ymaxlin = ymaxlin

        self.yminlog = yminlog
        self.ymaxlog = ymaxlog

        self.savehistos = savehistos
        self.samplestosave = samplestosave

    def __str__(self):
        return self.name

    def __add__(self, other):
        return str(self) + str(other)

    @classmethod
    def fromlist(cls, lst):
        if len(lst) == 5:
            return cls(name=lst[0], title=lst[1], axisrange=lst[2], rebin=lst[3], nbins=lst[4])
        elif len(lst) == 6:
            return cls(name=lst[0], title=lst[1], vartoplot=lst[2], axisrange=lst[3], rebin=lst[4], nbins=lst[5])
        else:
            raise NotImplementedError('cannot interpret lists of length ' + str(len(lst)))


class Ratio:
    def __init__(self, category, name, savehistos=False, variablestosave=None):

        if variablestosave is None:
            variablestosave = []

        if category not in ['ratio', 'cutsig', 'binsig', 'calibration', 'efficiency']:
            raise NotImplementedError('unknown ratio category')
        self.category = category

        self.name = name

        self.savehistos = savehistos
        self.variablestosave = variablestosave

    def __str__(self):
        return self.category + self.name

    def __add__(self, other):
        return str(self) + str(other)


if __name__ == '__main__':
    pass
