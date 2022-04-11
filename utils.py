import ROOT
from plotstyle import Plotstyle
from array import array
import ctypes

ROOT.gEnv.SetValue('RooFit.Banner', 0)

try:
    ROOT.EnableImplicitMT()
except AttributeError:
    pass


"""
To use RDataFrames to process TreeSamples a rather recent ROOT version is required, e.g.
source /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.24.06/x86_64-centos7-gcc48-opt/bin/thisroot.sh
"""


class PlotFactory:

    def __init__(self, inputfiles=None, inputpattern='VARIABLESAMPLE',
                 outputpath='.', outputpattern='VARIABLE', outputformat='pdf',
                 normalize=False, axes='linlog',
                 ylabel='Events', ylabelratio='Data / Prediction',
                 yaxisrangeratio=(0.0001, 1.9999),
                 ratiohlines=None,
                 linewidth=1, markersize=2,
                 poslegend=(0.4, 0.55, 0.95, 0.9), ncolumnslegend=1, boldlegend=False,
                 text='36 fb^{-1} (13 TeV)', extratext='#splitline{Work in progress}{Simulation}',
                 height=1280, width=None, ipos=11):

        if inputfiles is None:
            inputfiles = []

        if ratiohlines is None:
            ratiohlines = [1]

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

        self.outputpattern = outputpattern
        self.outputformat = outputformat

        self.normalize = normalize
        if axes not in ['lin', 'log', 'linlog']:
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

        # TODO: implement exporting histos/trees
        # fout = ROOT.TFile('/nfs/dust/cms/user/wolfmor/NTupleStuff/PUweights.root', 'recreate')
        # for ratio in self.ratiohistos:
        #     self.ratiohistos[ratio].Write()
        # fout.Close()

    def add_variables(self, variablelist):
        for variable in variablelist:
            self.add_variable(variable)

    def add_variable(self, variable):
        self.variables.append(variable)

    def add_samples(self, samplelist):
        for sample in samplelist:
            self.add_sample(sample)

    def add_sample(self, sample):
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

        for s in self.stacksamples + self.markersamples + self.linesamples:

            print s

            if isinstance(s, TreeSample):

                # TODO: implement histo weights with tCounter.GetEntries()

                # first book all histograms
                for v in self.variables:

                    if s.vectorselection is None:
                        if s.weight is None:
                            self.histos[v + s] = s.df.Histo1D(
                                ('h' + str(v + s), '', v.nbins, v.axisrange[0], v.axisrange[1]), v.vartoplot
                            )
                        else:
                            self.histos[v + s] = s.df.Define(
                                '_w', s.weight
                            ).Histo1D(
                                ('h' + str(v + s), '', v.nbins, v.axisrange[0], v.axisrange[1]), v.vartoplot, '_w'
                            )
                    else:
                        if s.weight is None:
                            self.histos[v + s] = s.df.Define(
                                v.name + '_pass', v.vartoplot + '[' + s.vectorselection + ']'
                            ).Histo1D(
                                ('h' + str(v + s), '', v.nbins, v.axisrange[0], v.axisrange[1]), v.name + '_pass'
                            )
                        else:
                            self.histos[v + s] = s.df.Define(
                                v.name + '_pass', v.vartoplot + '[' + s.vectorselection + ']'
                            ).Define(
                                '_w', s.weight
                            ).Histo1D(
                                ('h' + str(v + s), '', v.nbins, v.axisrange[0], v.axisrange[1]), v.name + '_pass', '_w'
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

            if len(self.stacksamples) > 0:
                self.sums[v] = self.histos[v + self.stacksamples[0]].Clone('sum' + v.name)
                self.sums[v].Reset('ICESM')
                for s in self.stacksamples:
                    self.sums[v].Add(self.histos[v + s])

            if self.normalize:
                for s in self.markersamples + self.linesamples:
                    self.histos[v + s].Scale(1. / self.histos[v + s].Integral())
                if len(self.stacksamples) > 0:
                    for s in self.stacksamples:
                        self.histos[v + s].Scale(1. / self.sums[v].Integral())
                    self.sums[v].Scale(1. / self.sums[v].Integral())

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

                elif r.category == 'cutsig':

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

                        if S > 0 and B > 0:
                            self.ratiohistos[v + r].SetBinContent(b, ROOT.RooStats.AsimovSignificance(S, B, syserrB*B))

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

            for ir, r in enumerate(self.ratios):
                self.ratiohistos[v + r].GetXaxis().SetTitle(v.title)
                self.ratiohistos[v + r].GetYaxis().SetTitle(self.ylabelratio)

                self.ratiohistos[v + r].SetLineWidth(self.linewidth)
                self.ratiohistos[v + r].SetLineStyle(self.histos[v.name + r.name.split(':')[0]].GetLineStyle())
                self.ratiohistos[v + r].SetLineColor(self.histos[v.name + r.name.split(':')[0]].GetLineColor())

                self.ratiohistos[v + r].SetMarkerSize(self.histos[v.name + r.name.split(':')[0]].GetMarkerSize())
                self.ratiohistos[v + r].SetMarkerColor(self.histos[v.name + r.name.split(':')[0]].GetMarkerColor())

                if ir == 0:
                    self.ratiohistos[v + r].SetMinimum(self.yaxisrangeratio[0])
                    self.ratiohistos[v + r].SetMaximum(self.yaxisrangeratio[1])

    def _draw_plots(self):

        print('\n# Draw Plots')

        for v in self.variables:

            print(v)

            lines = []
            for ratiohline in self.ratiohlines:
                lines.append(ROOT.TLine(v.axisrange[0], ratiohline, v.axisrange[1], ratiohline))
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

                        if r.name.split(':')[0] in [s.name for s in self.markersamples]:
                            self.ratiohistos[v + r].Draw('e x0 same')
                        else:
                            self.ratiohistos[v + r].Draw('hist same')

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

                        for r in self.ratios:

                            if r.name.split(':')[0] in [s.name for s in self.markersamples]:
                                self.ratiohistos[v + r].Draw('e x0 same')
                            else:
                                self.ratiohistos[v + r].Draw('hist same')

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


class Sample:
    def __init__(self, category, name,

                 title=None,

                 color=ROOT.kBlack,
                 linestyle=ROOT.kSolid,
                 fillstyle=1001):

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
                 fillstyle=1001):

        Sample.__init__(self, category, name, title, color, linestyle, fillstyle)

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

                 color=ROOT.kBlack,
                 linestyle=ROOT.kSolid,
                 fillstyle=1001):

        Sample.__init__(self, category, name, title, color, linestyle, fillstyle)

        self.vectorselection = vectorselection
        self.weight = weight

        if eventselection is None or eventselection == '':
            self.df = ROOT.RDataFrame(tree, files)
        else:
            self.df = ROOT.RDataFrame(tree, files).Filter(eventselection, ' selection for ' + self.name)

        # TODO: implement usage of TreeSamples without RDataFrame with:
        # TODO: deactivation of unused branches via setbranchstatus
        # TODO: entrylist used correctly?
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
                 axisrange=(0, 1),
                 rebin=1,
                 nbins=100,
                 yminlin=None, ymaxlin=None, yminlog=None, ymaxlog=None):

        self.name = name
        if title is None:
            self.title = name
        else:
            self.title = title

        if vartoplot is None:
            self.vartoplot = name
        else:
            self.vartoplot = vartoplot

        self.axisrange = axisrange
        self.rebin = rebin
        self.nbins = nbins

        self.yminlin = yminlin
        self.ymaxlin = ymaxlin

        self.yminlog = yminlog
        self.ymaxlog = ymaxlog

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
            return None


class Ratio:
    def __init__(self, category, name):

        if category not in ['ratio', 'cutsig']:
            raise NotImplementedError('unknown ratio category')
        self.category = category

        self.name = name

    def __str__(self):
        return self.name

    def __add__(self, other):
        return str(self) + str(other)


if __name__ == '__main__':

    pass
