#!python
# -*- coding: UTF-8 -*-
# %%
import os,sys

#%%
try:
    os.chdir(os.path.join(os.getcwd(), '../../../tmp'))
    print(os.getcwd())
except:
    pass
# %% [markdown]
# # JetFake Transformation for mt QCD 2017

# %%
import numpy as np
import matplotlib.pyplot as plt
import re, yaml
import copy
import itertools

# %%
import ROOT
ROOT.ROOT.EnableImplicitMT(15)

# %%
import logging

# try:
#     logger
# except NameError:
logger = logging.getLogger()
logger.setLevel(logging.INFO)
formatter = logging.Formatter("%(name)s - %(levelname)s - %(message)s")

handler = logging.StreamHandler()
handler.setFormatter(formatter)
logger.addHandler(handler)

# %%


trainingVarL = {
    "bpt_1", "dijetpt", "eta_1", "jdeta", "jpt_1", "jpt_2", "ME_q2v1",
    "ME_q2v2", "mjj", "m_sv", "m_sv_puppi", "mTdileptonMET",
    "mTdileptonMET_puppi", "m_vis", "nbtag", "njets", "pt_1", "pt_2", "ptvis"
}
ffVarL = {"decayMode_2", "pt_2", "njets", "m_vis", "mt_1", "iso_1"}
relVarL = trainingVarL.union(ffVarL).union({"eventWeight"})


# %%
from fake_factor_derivation.regions import ParSpaceRegion,plotvarInRegions,plotvar

# %%

jfQCD = ParSpaceRegion("2016", "mt", "QCD")
jfWjet = ParSpaceRegion("2016", "mt", "W+jets")
jfttbar = ParSpaceRegion("2016", "mt", "ttbar")

# %%
lastvals = {}
class DynBins(object):
    def __init__(self, slarr, blarr, predictionVars):
        self.slarr = slarr
        self.blarr = blarr
        print(len(self.slarr[predictionVars[0]]))
        #print(len(self.blarr[predictionVars[0]]))

        self.predictionVars = predictionVars
        ### histogram both regions in the generated bins
        ### use dynamic binning, from the bl region -> no 1/0 bins
        for var in self.predictionVars:
            print(var)
            print(np.histogram(self.blarr[var]))
        self.binbordersD = {
            var: self.dyn1dBins(self.blarr[var], self.blarr["eventWeight"])
            for var in self.predictionVars
        }
        self.blh, self.edges = np.histogramdd(
            [self.blarr[var] for var in self.predictionVars],
            weights=self.blarr["eventWeight"],
            bins=np.array(
                [self.binbordersD[var] for var in self.predictionVars]))
        self.slh, _ = np.histogramdd(
            [self.slarr[var] for var in self.predictionVars],
            weights=self.slarr["eventWeight"],
            bins=np.array(
                [self.binbordersD[var] for var in self.predictionVars]))
        self.bincentersL = [
            self.binbordersToCenters(self.binbordersD[var])
            for var in self.predictionVars
        ]

        ## currently, the merge of the 1d bins into a nbinning lev
        ## join bins  a var that do not have sufficient statistic due to the "multiplication" of the 1d binnings
        ## the binning should be reduced in the variable with most bins first:
        predictionVarsByLengthL = list(
            enumerate(copy.deepcopy(self.predictionVars)))
        predictionVarsByLengthL.sort(key=lambda k: len(self.binbordersD[k[1]]),
                                     reverse=True)
        for ivar, var in predictionVarsByLengthL:
            bb = self.binbordersD[var]
            bc = self.bincentersL[ivar]
            ibin = len(bc) - 1
            while (ibin >= 1):
                slplane = self.slh.take(indices=ibin, axis=ivar)
                blplane = self.blh.take(indices=ibin, axis=ivar)
                ### continue joining bins of the selected variables as long as 50% of the events at the selected index dont have enough events
                while (any([
                        np.sum(plane < 10.) > 0.5 * len(plane)
                        for plane in (slplane, blplane)
                ])):
                    logger.debug(
                        "for variable {} deleting merging bin {} with previous bin"
                        .format(var, ibin))
                    ##  the new 'row' is the element wise sum of the already existing row + the new one
                    self.slh = self.replaceWithSumToPreviousIndex(
                        self.slh, ivar, ibin)
                    self.blh = self.replaceWithSumToPreviousIndex(
                        self.blh, ivar, ibin)
                    # set the upper bin border to the next bin
                    bc[ibin - 1] = bb[ibin]

                    ## hopefully this does del bc[ibin], bb[ibin]
                    bc = np.delete(bc, ibin)
                    bc = np.delete(bb, ibin)
                    ### count bin down
                    ibin -= 1
                    ### set the planes again, so the while can check them again
                    slplane = self.slh.take(indices=ibin, axis=ivar)
                    blplane = self.blh.take(indices=ibin, axis=ivar)
                ibin -= 1

        print(self.slh.sum())
        #print(self.blh.sum())
        ## filter indices of bins with less than 10 events
        ## generate a list of index tuples for slh, blh
        self.idxs = list(np.ndindex(*self.slh.shape))
        self.valididxs = [
            idx for idx in self.idxs if self.slh[idx] > 5 and self.blh[idx] > 5
        ]
        ## select the bins with mor then 10 events in signal and bkg region as x values for our regression
        self.xv = np.array(list(
            itertools.product(*np.array(self.bincentersL))))[self.valididxs]
        self.yv = np.array(
            map(lambda idx: self.slh[idx] / self.blh[idx], self.valididxs))
        self.plotyMat()

    def plotyMat(self):
        self.ymat = np.full(self.slh.shape, np.nan)
        for idx, y in zip(self.valididxs, self.yv):
            self.ymat[idx] = y
        import seaborn as sns
        import matplotlib.pyplot as plt
        plt.close()
        plt.cla()
        plt.clf()
        fit, axes = plt.subplots(3, 1)
        plt.gcf().set_size_inches((6, 8))

        g = sns.heatmap(self.ymat, ax=axes[0])
        g.set
        g.set_xticklabels(self.bincentersL[0].round(0))
        g.set_yticklabels(self.bincentersL[1].round(0))
        g = sns.heatmap(self.slh, ax=axes[1])
        g.set_xticklabels(self.bincentersL[0].round(0))
        g.set_yticklabels(self.bincentersL[1].round(0))
        g = sns.heatmap(self.blh, ax=axes[2])
        g.set_xticklabels(self.bincentersL[0].round(0))
        g.set_yticklabels(self.bincentersL[1].round(0))
        plt.xlabel(self.predictionVars[0])
        plt.ylabel(self.predictionVars[1])
        plt.show()

        # # x values for the fit
        # self.xsc = np.arange(
        #     self.binborders[0], self.binborders[-1] + 200,
        #     200. / np.abs(self.binborders[-1] - self.binborders[0]))

    def dyn1dBins(self, arr, weightsarr):
        import numpy.lib.recfunctions as rfn
        varS = rfn.merge_arrays([arr, weightsarr])
        #sort the values first
        varS.sort(order="f0")

        marker = len(varS) - 1
        markerval = varS["f0"][marker]
        stdbinwith = 10
        # we start from the maximium so it is out first binborder
        binborders = [varS["f0"][-1]]
        varmin = varS["f0"][0]

        ##
        while (0 != marker):
            curbinsum = 0
            #make sure the sum of the event weights in the bin is at least 100
            while (curbinsum < 100 and marker>0):
                #add the weight of the current event
                curbinsum += varS["f1"][marker]
                marker -= 1
                markerval = varS["f0"][marker]
            traveledDistance = binborders[-1] - markerval

            ##  make sure the size of the bin is at letzt stdbinsize
            while (traveledDistance < stdbinwith and marker != 0):
                traveledDistance += np.abs(markerval - varS["f0"][marker - 1])
                marker -= 1
                markerval = varS["f0"][marker]
            binborders.append(markerval)

        binborders.reverse()
        return binborders

    def binbordersToCenters(self, arr):
        return np.array([(arr[i] + arr[i + 1]) / 2.0
                         for i in range(len(arr) - 1)])

    ## return a tuple indexing the axis , to with naxis=3, axis=1, fix=3
    ## this returns (:,3,:)
    def getIndexTuple(self, naxis, axis, fix):
        return (tuple([[slice(None), fix][axi == axis]
                       for axi in range(naxis)]))

    ### takes and array, the selected axis and an index
    ### replaces the values along that axis for index idx-1 with the sum of values along this axis with index idx-1 and idx
    ###
    def replaceWithSumToPreviousIndex(self, arr, ivar, idx):
        plane = arr.take(indices=idx, axis=ivar) + arr.take(indices=idx - 1,
                                                            axis=ivar)
        s = self.getIndexTuple(len(arr.shape), axis=ivar, fix=idx - 1)
        arr[s] = plane
        return (np.delete(arr, obj=idx, axis=ivar))


# %%
class fakefactor(object):
    def __init__(self, jfRS, filterstring, predictionVars):
        ## get the numpy arrays from the jfObject
        self.jfRS = jfRS
        self.filterstring = filterstring
        self.predictionVars = predictionVars
        self.meta = jfRS.meta
        self.blarr = jfRS.DR_bl.RDF.Filter(filterstring).AsNumpy(
            predictionVars + ["eventWeight"])
        self.slarr = jfRS.DR_sl.RDF.Filter(filterstring).AsNumpy(
            predictionVars + ["eventWeight"])

        binningObj = DynBins(self.slarr, self.blarr, self.predictionVars)
        self.slh = binningObj.slh
        self.blh = binningObj.blh
        self.xv = binningObj.xv
        self.yv = binningObj.yv
        self.valididxs = binningObj.valididxs
        ### borders for comparing the old fakefaktors

        # if jfRS.meta["bkg"]=="QCD":
        #     self.binborders=np.array([16.8, 20.0, 23.2, 26.4, 29.599999999999998, 32.800000000000004, 36.0, 39.2, 42.400000000000006, 45.60000000000001, 48.800000000000004, 52.0, 55.2, 58.400000000000006, 61.60000000000001, 64.8, 68.00000000000001, 71.2, 74.4, 77.60000000000001, 80.8, 84.00000000000001, 87.2, 90.4, 93.60000000000001, 96.80000000000001], np.float64)
        # #elif jfRS.meta["bkg"]=="W+jets": self.binborders=np.array([29,32.5,36,42,65,500], np.float32)
        # else: raise Exception

        # #self.applyCorrections()
        self.yerror = self.calcyerror()
        # self.npars, self.pars = self.getFitPars()
        # ## number of standar ddeviations used for calculating the error bands
        # self.nstd = 1
        # ## number of parameters in the fit function
        # self.fitsettings = {
        #     "f": self.fitfunction,
        #     "xdata": self.xv,
        #     "sigma": self.yerror,
        #     "bounds":
        #     np.array([self.pars[key][0] for key in self.pars.keys()]).T,
        #     "p0": [self.pars[key][1] for key in self.pars.keys()],
        #     "maxfev": 100000,
        #     "loss": "soft_l1",
        #     "method": "trf"
        # }
        # self.popt = self.bestfit()
        # self.ff = self.fitfunction(self.xv, *popt)
        # self.fitresy = self.fitfunction(self.xsc, *popt)
        # self.rooFitPercV = self.roofitError()
        # self.plotFF()

    def calcyerror(self):
        ### Error calculation
        ## σ_y= abs(y)sqrt((∂_s y)^2*σ_s^2+(∂_b y)^2*σ_b^2+2(∂_s y)(∂_b y)σ_s^2*σ_b^2)
        # y,s,b>0 y=s/b
        # = y*sqrt((σ_s/s)^2+(σ_b/b)^2+σ_s*σ_b/(s*b))
        # poisson error
        # = y*sqrt(sqrt(1/s)^2+sqrt(1/b)^2-2*σ_sb/(s*b))
        # = y*sqrt(sqrt(1/s)^2+sqrt(1/b)^2)
        yerror = np.array([
            self.yv[i] * np.sqrt(1 / self.slh[self.valididxs[i]] +
                                 1 / self.blh[self.valididxs[i]])
            for i in range(len(self.yv))
        ])
        return yerror

    def getFitPars(self):
        center = np.mean(self.binborders)
        width = float(self.binborders[-1] - self.binborders[1])
        heigth = np.max(self.yv) - np.min(self.yv)
        from collections import OrderedDict
        pars = OrderedDict()
        pars["c1"] = [(0, 1), min(1, np.mean(self.yv))]
        pars["c2"] = [(np.mean(self.yv) - 2 * heigth,
                       np.mean(self.yv) + 2 * heigth),
                      np.mean(self.yv)]
        #pars["l1"]=[(0,np.max(self.xv)),center]
        pars["l1"] = [(0, np.inf), center]
        pars["s1"] = [(0, np.inf), 1 / width]
        pars["b1"] = [(-np.inf, np.inf), 1 / width]
        #pars["a1"]=[(-(np.mean(self.yv)+2*heigth),np.mean(self.yv)+2*heigth),np.mean(self.yv)]
        pars["a1"] = [(-np.inf, np.inf), np.mean(self.yv)]
        npars = len(pars.keys())
        return npars, pars

    def fitfunction(self, xe, c1, c2, l1, s1, b1, a1):
        nx = xe.size
        xs = 4 / (self.xv[-1] - self.xv[-2]) * (xe - self.xv[-2])
        lastvals["xs"] = xs
        sig = 1 / (1 + np.exp(-xs))
        lastvals["sig"] = sig
        const = np.array([c1 for x in range(nx)])
        lastvals["sig"] = const
        x1 = s1 * (xe - l1)
        lastvals["x1"] = x1
        #e1=np.abs(2)*2
        #poly=a1*x1**2/(x1**2+1)+xe*b1
        poly = a1 + xe * b1
        lastvals["poly"] = poly
        p1 = const + poly

        p2 = np.array([c2 for x in range(nx)])

        funcsum = p1 * (1 - sig) + p2 * sig
        return funcsum

    def bestfit(self):
        from scipy.optimize import curve_fit
        import seaborn as sns
        ### covariance matrix based error calculation
        popt, pcov = curve_fit(ydata=self.yv, **self.fitsettings)
        perr = np.sqrt(np.diag(pcov))
        fity, fityu, fityd = [
            self.fitfunction(self.xsc, *par) for par in
            [popt, popt + self.nstd * perr, popt - self.nstd * perr]
        ]
        print "pot pars"
        for i, key in enumerate(self.pars.keys()):
            print [key, popt[i], self.pars[self.pars.keys()[i]]]
        # print "Cov matrix"
        # sns.heatmap(pcov)
        # plt.show()
        print "diag error"
        print [(key, perr[i]) for i, key in enumerate(self.pars.keys())
               if perr[i] != 0]

        # plt.cla()
        # plt.clf()
        # plt.gcf().set_size_inches((6,5))
        # plt.style.use('seaborn-darkgrid')
        # #plt.xscale("log")
        # plt.xlim(self.binborders[0],self.binborders[-1]*2)
        # plt.ylim(0,np.max(self.yv)+.1)
        # plt.errorbar(self.xv,self.yv,yerr=self.yerror,fmt="o", capsize=3, marker="x")
        # plt.plot(self.xsc,fity, label="fit")
        # plt.show()
        return (popt)

    def roofitError(self):
        from scipy.optimize import curve_fit
        from scipy.stats import norm
        ### Roofit - like error calutation run the fit multiple times with y values sampled
        rep = 20
        ysampled = (self.yerror * np.random.randn(rep, len(self.yv))) + self.yv
        parMat = np.ndarray(shape=[rep, self.npars])
        for i in range(rep):
            popt, pcov = curve_fit(ydata=ysampled[i], **self.fitsettings)
            parMat[i] = popt

        ## init array for collecting the collecting the upper/median/lower limit for each xs value
        fitres = np.ndarray(shape=[len(self.xsc), 3], dtype=float)
        yfitvals = np.array(
            [self.fitfunction(self.xsc, *parD_) for parD_ in parMat])
        for xi, segment in enumerate(yfitvals.T):
            fitres[xi] = np.percentile(segment, [
                norm.cdf(-self.nstd, 0, 1) * 100, 50,
                norm.cdf(self.nstd, 0, 1) * 100
            ])
        return fitres

    def plotFF(self):
        matplotlib.rcdefaults()
        plt.close()
        plt.cla()
        plt.clf()
        plt.gcf().set_size_inches((6, 5))
        plt.style.use('seaborn-darkgrid')
        plt.xscale("log")
        plt.xlim(self.binborders[0], self.binborders[-1] + 200)
        plt.ylim(0, np.max(self.fitresy) + .2)
        plt.errorbar(self.xv,
                     self.yv,
                     yerr=self.yerror,
                     fmt="o",
                     capsize=3,
                     marker="x")

        plt.plot(self.xsc, self.fitresy, label="fit")
        # oldffxv, oldffyv=self.getOldFF()
        # print oldffxv
        # print oldffyv
        # plt.plot(oldffxv,oldffyv,label="oldFF", marker='o', markersize=3)
        plt.gca().fill_between(self.xsc,
                               self.rooFitPercV[:, 2],
                               self.rooFitPercV[:, 0],
                               alpha=.25,
                               label=str(self.nstd) + "-sigma roofitMethod")
        plt.legend()
        plt.title(self.var + " for " + str(self.meta) + " with " +
                  self.filterstring)
        for yvar in self.binborders:
            plt.axvline(yvar, ls=":", color="darkgrey")
        plt.tight_layout()
        plt.savefig(
            "/home/mscham/fakeFaktors/plots/" +
            "-".join([self.meta[key] for key in ["era", "channel", "bkg"]]) +
            "-" + self.filterstring)
        plt.show()

    # def getOldFF(self):
    #     f=ROOT.TFile("/portal/ekpbms2/home/jbechtel/fakefactors/new/CMSSW_8_0_25/src/ViennaTool/sim/mt/CR_QCD_pt_data.root")
    #     slh=f.Get("hh_t_pt")
    #     blh=f.Get("hh_l_pt")
    #     # %%
    #     bincenters=[slh.GetBinCenter(i) for i in range(slh.GetNbinsX())]
    #     if bincenters!=[blh.GetBinCenter(i) for i in range(blh.GetNbinsX())]:
    #         raise Exception
    #     slr=np.array([slh.GetBinContent(i) for i in range(slh.GetNbinsX())])
    #     blr=np.array([blh.GetBinContent(i) for i in range(blh.GetNbinsX())])
    #     xv=[]
    #     ff=[]
    #     for i in range(len(bincenters)):
    #         if np.abs(blr[i])>0.0:
    #             xv.append(bincenters[i])
    #             ff.append(slr[i]/blr[i])
    #     return np.array(xv),np.array(ff)


#%%
a = fakefactor(jfQCD, "njets==0", ["pt_2", "m_vis"])
#a = fakefactor(jfWjet, "1==1", ["pt_2", "m_vis"])


# %%
