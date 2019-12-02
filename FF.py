#!python
# -*- coding: UTF-8 -*-
# %%
import os
try:
    os.chdir(os.path.join(os.getcwd(), '../../../tmp'))
    print(os.getcwd())
except:
    pass
# %% [markdown]
# # JetFake Transformation for mt QCD 2017

# %%
#import uproot
import numpy as np
#import pandas as pd
import matplotlib.pyplot as plt
import re, yaml
import copy
import itertools

# %%
import ROOT
ROOT.ROOT.EnableImplicitMT(15)
from shape_producer.cutstring import Cut, Cuts, Weights, Weight
from shape_producer.channel import *

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
from shape_producer.estimation_methods_2017 import DataEstimation, ZTTEstimation, ZJEstimation, ZLEstimation, TTLEstimation, TTJEstimation, TTTEstimation, VVTEstimation, VVJEstimation, VVLEstimation, WEstimation, ggHEstimation, qqHEstimation, EWKZEstimation, ZTTEmbeddedEstimation, NewFakeEstimationTT, NewFakeEstimationLT
from shape_producer.era import Run2016, Run2017, Run2018

# %%
database = "datasets/datasets.json"
base_path = "/ceph/htautau/2017/ntuples/"
channelDict = {
    "2016": {
        "mt": MTSM2016(),
        "et": ETSM2016(),
        "tt": TTSM2016(),
        "em": EMSM2016()
    },
    "2017": {
        "mt": MTSM2017(),
        "et": ETSM2017(),
        "tt": TTSM2017(),
        "em": EMSM2017()
    },
    "2018": {
        "mt": MTSM2018(),
        "et": ETSM2018(),
        "tt": TTSM2018(),
        "em": EMSM2018()
    }
}
eraD = {
    "2016": Run2016(database),
    "2017": Run2017(database),
    #"2018":Run2018(database)
}


# %%
def cutDB(meta, cutName):
    if cutName == "tauDisc":
        if meta["channel"] in ["et", "mt"]:
            #c=Cut("byTightIsolationMVArun2017v2DBoldDMwLT2017_2>0.5", "tauDisc")
            c = Cut("byTightDeepTau2017v2p1VSjet_2>0.5", "tauDisc")
            if not meta["signalLikeRegion"]:
                c.invert()
            return c
        elif meta["channel"] == "tt":
            #### tight -> medium ??? like in AN
            # c=Cut("(byMediumIsolationMVArun2017v2DBoldDMwLT2017_2>0.5&&byMediumIsolationMVArun2017v2DBoldDMwLT2017_1<0.5)||(byMediumIsolationMVArun2017v2DBoldDMwLT2017_1>0.5&&byMediumIsolationMVArun2017v2DBoldDMwLT2017_2<0.5)", "tauDisc")
            c = Cut(
                "(byMediumDeepTau2017v2p1VSjet_2>0.5&&byMediumDeepTau2017v2p1VSjet_1<0.5)||(byMediumDeepTau2017v2p1VSjet_1>0.5&&byMediumDeepTau2017v2p1VSjet_2<0.5)",
                "tauDisc")
            if meta["signalLikeRegion"]:
                c.invert()
            return c
        else:
            print cutName
            print meta
            raise Exception
    elif cutName == "signalEnrichment":
        if meta["bkg"] == "QCD":
            if not meta["determinationRegion"]:
                ## opposite sign
                return Cut("q_1*q_2<0", "signalEnrichment")
            else:
                return Cut("q_1*q_2>0", "signalEnrichment")
        elif meta["bkg"] == "W+jets":
            if meta["determinationRegion"]:
                return Cut("mt_2<50", "signalEnrichment")
            else:
                return Cut("mt_2>70", "signalEnrichment")
        elif meta["bkg"] == "ttbar":
            return None
        else:
            print cutName
            print meta
            raise Exception
    elif cutName == "DR":
        if meta["bkg"] == "QCD":
            return (Cut("(iso_1>0.05)*(mt_1<40)", "QCDffpresel"))
        elif meta["bkg"] == "W+jets":
            return (Cut("nbtag==0", "W+jetsffpresel"))
    else:
        logger.fatal("No cut with name {} from {}".format(cutName, str(meta)))
        raise Exception
    """elif cutName=="baseline":
        #### untested, not in the current ntuples
        cutL=[]
        if eraName=="2017":
            if channelName=="mt":
                ###μ
                anyTriggerL=["trg_singlemuon_24", "trg_singlemuon_27","trg_crossmuon_mu20tau27"]
                cutL.append(Cut("||".join([x+">0.5" for x in  anyTriggerL]), "anyFilterFF"))
                ### add muon id here
                ##guess:
                cutL.append(Cut("id_m_medium_1>.5","muon_id"))
                ######verified:
                cutL.append(Cut("pt_1>25", "muon_pt"))
                cutL.append(Cut("-2.1<eta_1&&eta_1<2.1", "muon_eta"))
                cutL.append(Cut("iso_1<0.15", "muon_iso"))
                cutL.append(Cut("dZ_1<.2&&d0_1<.045", "muon_impact"))
                #cutL.append(Cut("dilepton_veto<0.5", "dilepton_veto"))

                ##τ
                #this migth not be the right one: it's always 1 in our analysis
                #cutL.append(Cut("decayModeFinding_2>0.5","bydecay"))
                ####
                cutL.append(Cut("pt_2>23", "tau_pt"))
                cutL.append(Cut("-2.3<eta_2&&eta_2<2.3", "tau_eta"))
                ### Laut AN:
                #cutL.append(Cut("byTightIsolationMVArun2v1DBoldDMwLT_2>.5", "tau_iso"))
                ## laut Janek
                cutL.append(Cut("byTightIsolationMVArun2017v2DBoldDMwLT2017_2>0.5","tau_iso213"))
                cutL.append(Cut("dZ_2<.2","tau_impact"))



                ### after this???
                #cutL.append(Cut("againstMuonTight3_2>0.5", "againstMuonDiscriminator"))
                #cutL.append(Cut("againstElectronVLooseMVA6_2>0.5", "againstElectronDiscriminator"))


                #cutL.append(Cut("byVLooseIsolationMVArun2017v2DBoldDMwLT2017_2>.5","VLooseTauIso"))
                Cut("extraelec_veto<0.5", "extraelec_veto"),
            Cut("extramuon_veto<0.5", "extramuon_veto"),
            Cut("byTightIsolationMVArun2v1DBoldDMwLT_2>0.5", "tau_iso"),
            Cut("trg_singlemuon==1", "trg_singlemuon"))

            # elif channelName=="et":
            #     anyTriggerL=["Ele27_WPTight_Gsf", "Ele32_WPTight_Gsf", "Ele35_WPTight_Gsf", "Ele24_eta2p1", "_WPTight_Gsf_Loose", "ChargedIsoPFTau30", "_eta2p1 CrossL1"]
            # elif channelName=="tt":
            #     anyTriggerL=["DoubleMediumChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg","DoubleTightChargedIsoPFTau40_Trk1_eta2p1_Reg","DoubleTightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg"]
            #### trigger var not translated
            #cutL.append(Cut("("+"||".join(anyTriggerL)+")","baseline")
        return cutL
        """


trainingVarL = {
    "bpt_1", "dijetpt", "eta_1", "jdeta", "jpt_1", "jpt_2", "ME_q2v1",
    "ME_q2v2", "mjj", "m_sv", "m_sv_puppi", "mTdileptonMET",
    "mTdileptonMET_puppi", "m_vis", "nbtag", "njets", "pt_1", "pt_2", "ptvis"
}
ffVarL = {"decayMode_2", "pt_2", "njets", "m_vis", "mt_1", "iso_1"}
relVarL = trainingVarL.union(ffVarL).union({"eventWeight"})


# %%
class jetFakeEstimation(object):
    def __init__(self, eraName, channelName, bkgName):
        self.meta = {"era": eraName, "channel": channelName, "bkg": bkgName}
        self.era = eraD[eraName]
        self.channel = copy.deepcopy(channelDict[eraName][channelName])

        if self.meta["era"] not in ["2016", "2017"]: raise Exception

        ### remove old isolation cuts and add the very Loose isolation, that is needed to to exclude other backgrounds from all regions
        # for cut_ in self.channel.cuts.names:
        #    self.channel.cuts.remove(cut_)
        # for cut_ in cutDB(channelName,"baseline", eraName="2017"):
        #    self.channel.cuts.add(cut_)
        if self.meta["channel"] in ["et", "mt"]:
            self.channel.cuts.remove("tau_iso")
            self.channel.cuts.remove("os")
            ##switch to deeptauID
            self.channel.cuts.get("againstMuonDiscriminator"
                                  ).variable = "byTightDeepTau2017v2p1VSmu_2"
            self.channel.cuts.get("againstElectronDiscriminator"
                                  ).variable = "byVLooseDeepTau2017v2p1VSe_2"
            self.channel.cuts.remove("trg_selection")
            self.channel.cuts.add(
                Cut(
                    "(pt_2>30 && ((trg_singlemuon == 1) || (trg_mutaucross == 1)))",
                    "trg_selection"))
            #self.channel.cuts.add(Cut("byVLooseIsolationMVArun2017v2DBoldDMwLT2017_2>.5","VLooseTauIso"))
            self.channel.cuts.add(
                Cut("byVLooseDeepTau2017v2p1VSjet_2>.5", "VLooseTauIso"))

        if self.meta["channel"] == "tt":
            self.channel.cuts.remove("tau_1_iso")
            self.channel.cuts.remove("tau_2_iso")
            self.channel.cuts.remove("os")
            ##switch to deeptauid
            raise Exception
            # self.channel.cuts.add(Cut("byVLooseIsolationMVArun2017v2DBoldDMwLT2017_1>.5 &&byVLooseIsolationMVArun2017v2DBoldDMwLT2017_2>.5","VLooseTauIso"))
            self.channel.cuts.add(
                Cut(
                    "byVLooseDeepTau2017v2p1VSjet_1>.5 &&byVLooseDeepTau2017v2p1VSjet_2>.5",
                    "VLooseTauIso"))

        # ###OldFFWeights
        # if self.meta["channel"] in ["et","mt"]:
        #     self.fakeWeightstring="ff2_nom"
        # elif self.meta["channel"]=="tt":
        #     self.fakeWeightstring="(0.5*ffcutDB(channelName,cutName,bkgName,signalLikeRegion, eraName="2017")1_nom*(byTightIsolationMVArun2017v2DBoldDMwLT2017_1<0.5)+0.5*ff2_nom*(byTightIsolationMVArun2017v2DBoldDMwLT2017_2<0.5))"

        self.CutString = self.channel.cuts.expand()
        print self.CutString
        self.estimation = DataEstimation(
            era=eraD[eraName],
            directory="/ceph/htautau/" + eraName + "/ntuples/",
            friend_directory=[
                "/ceph/htautau/" + eraName + "/" + ft
                for ft in ["ff_friends/", "mela_friends/", "svfit_friends/"]
            ],
            channel=self.channel)
        self.createRDF()
        if self.meta["bkg"] == "ttbar":
            self.createttbarRDF()
        self.createClosureRDF()
        self.SR = jFRegion(self.channel,
                           self.meta,
                           self.RDF,
                           signalLikeRegion=True,
                           determinationRegion=False)

        self.AR = jFRegion(self.channel,
                           self.meta,
                           self.RDF,
                           signalLikeRegion=False,
                           determinationRegion=False)

        ## Define the Background and Signal like Determination Region
        if self.meta["bkg"] != "ttbar":
            self.DR_sl = jFRegion(self.channel,
                                  self.meta,
                                  self.RDF,
                                  signalLikeRegion=True,
                                  determinationRegion=True)
            self.DR_bl = jFRegion(self.channel,
                                  self.meta,
                                  self.RDF,
                                  signalLikeRegion=False,
                                  determinationRegion=True)
        else:
            self.DR_sl = jFRegion(self.channel,
                                  self.meta,
                                  self.ttbarRDF,
                                  signalLikeRegion=True,
                                  determinationRegion=True)
            self.DR_bl = jFRegion(self.channel,
                                  self.meta,
                                  self.ttbarRDF,
                                  signalLikeRegion=False,
                                  determinationRegion=True)

        ## Define RDFs for the closure correction
        self.Closure_sl = jFRegion(self.channel,
                                   self.meta,
                                   self.ClosureRDF,
                                   signalLikeRegion=True,
                                   determinationRegion=True)
        self.Closure_bl = jFRegion(self.channel,
                                   self.meta,
                                   self.ClosureRDF,
                                   signalLikeRegion=False,
                                   determinationRegion=True)

    def createRDF(self):
        tree_path = self.channel.name + "_nominal/ntuple"
        self.rdfFilePath = "/ceph/mscham/data/fakefaktorRDFs/{}-{}.root".format(
            self.meta["era"], self.meta["channel"])
        if os.path.isfile(self.rdfFilePath) and os.path.exists(
                self.rdfFilePath):
            self.RDF = ROOT.RDataFrame(tree_path, self.rdfFilePath)
        else:
            logger.info("Creating new RDF!")
            #exit(1)
            self.chain = ROOT.TChain(tree_path)  ##
            self.chainsD = {}
            dontDelMyChainsL = []
            for i, ntupleFilename in enumerate(self.estimation.get_files()):
                #if i!=0: continue
                #### get the file basename
                filename = os.path.basename(os.path.normpath(ntupleFilename))
                ### instance the chain with selector
                if "ntuple" not in self.chainsD:
                    self.chainsD["ntuple"] = ROOT.TChain(tree_path)
                    logger.info("creating ntuple chain")

                if not os.path.exists(ntupleFilename):
                    logger.fatal("File does not exist: {}".format(path))
                    raise Exception
                logger.info("Adding ntuple:{}".format(ntupleFilename))
                self.chainsD["ntuple"].AddFile(ntupleFilename)

                j = 0
                # Make sure, that friend files are put in the same order together
                for friend in self.estimation.get_friend_files(
                ):  # friend in mela, svfit, ...
                    friendFileName = friend[i]
                    if not os.path.exists(friendFileName):
                        logger.fatal(
                            "File does not exist: {}".format(friendFileName))
                        raise Exception

                    friendVarL = ROOT.RDataFrame(
                        tree_path, str(friendFileName)).GetColumnNames()
                    if "ff2_nom" in friendVarL: friendName = "ff"
                    elif "ME_phi" in friendVarL: friendName = "mela"
                    elif "m_sv" in friendVarL: friendName = "svfit"
                    else:
                        friendName = str(j)
                        j = j + 1

                    if friendName not in self.chainsD:
                        self.chainsD[friendName] = ROOT.TChain(tree_path)
                        logger.info("creating friend tree chain for " +
                                    friendName)
                    logger.info("Adding {} as {} friend.".format(
                        friendFileName, friendName))
                    self.chainsD[friendName].AddFile(friendFileName)

            ### Collect the friend Chains
            self.chain.Add(self.chainsD["ntuple"])
            for key in self.chainsD.keys():
                if key == "ntuple": continue
                logger.info("Add " + key + " to chain")
                self.chain.AddFriend(self.chainsD[key], key)
            chain_numentries = self.chain.GetEntries()
            if chain_numentries == 0:
                logger.fatal(
                    "Chain (before skimming) does not contain any events.")
                raise Exception
            logger.debug("Found {} events.".format(chain_numentries))
            ### Convert Chain to an RDF
            self.RDF = ROOT.RDataFrame(self.chain)

            opt = ROOT.ROOT.RDF.RSnapshotOptions()
            opt.fMode = "RECREATE"
            self.RDF.Snapshot(tree_path, self.rdfFilePath, ".*", opt)

        logger.debug("Skim events with cut string: {}".format(self.CutString))
        self.RDF = self.RDF.Filter(self.CutString)
        self.eventCount = self.RDF.Count().GetValue()

        logger.debug("RDF after basic cuts contains {} events.".format(
            self.eventCount))
        self.varL = list(self.RDF.GetColumnNames())

    def createttbarRDF(self):
        tree_path = self.channel.name + "_nominal/ntuple"
        self.rdfttbarFilePath = "/ceph/mscham/data/fakefaktorRDFs/{}-{}-ttbar.root".format(
            self.meta["era"], self.meta["channel"])
        if os.path.isfile(self.rdfttbarFilePath) and os.path.exists(
                self.rdfttbarFilePath):
            self.ttbarRDF = ROOT.RDataFrame(tree_path, self.rdfttbarFilePath)
        else:
            raise Exception

        logger.debug("Skim ttbar with cut string: {}".format(self.CutString))
        self.ttbarRDF = self.RDF.Filter(self.CutString)
        self.ttbareventCount = self.ttbarRDF.Count().GetValue()

        logger.debug("ttbarRDF after basic cuts contains {} events.".format(
            self.ttbareventCount))
        self.ttbarvarL = list(self.ttbarRDF.GetColumnNames())

    def createClosureRDF(self):
        tree_path = self.channel.name + "_nominal/ntuple"
        FilePath = "/ceph/mscham/data/fakefaktorRDFs/{}-{}-Closure.root".format(
            self.meta["era"], self.meta["channel"])
        if os.path.isfile(FilePath) and os.path.exists(FilePath):
            self.ClosureRDF = ROOT.RDataFrame(tree_path, FilePath)
        else:
            raise Exception

        logger.debug("Skim Closure with cut string: {}".format(self.CutString))
        self.ClosureRDF = self.RDF.Filter(self.CutString)
        self.ClosureEventCount = self.ClosureRDF.Count().GetValue()

        logger.debug("ClosureRDF after basic cuts contains {} events.".format(
            self.ClosureEventCount))
        self.ClosureVarL = list(self.ClosureRDF.GetColumnNames())


class jFRegion(jetFakeEstimation):
    def __init__(self, channel, meta, RDF, signalLikeRegion,
                 determinationRegion):
        self.channel = copy.deepcopy(channel)
        self.meta = copy.deepcopy(meta)
        self.meta["signalLikeRegion"] = signalLikeRegion
        self.meta["determinationRegion"] = determinationRegion

        cutNames = []
        if determinationRegion:
            cutNames.append("DR")
        cutNames.append("tauDisc")
        cutNames.append("signalEnrichment")
        cutsToAdd = [cutDB(self.meta, cutName) for cutName in cutNames]
        for cut_ in cutsToAdd:
            # ttbar + signalEnrichment returns none
            if cut_ != None:
                self.channel.cuts.add(cut_)

        self.CutString = self.channel.cuts.expand()
        self.RDF = RDF.Filter(self.CutString)


def plotvarInRegions(jfE, var):
    npmatAR = jfE.AR.RDF.AsNumpy([var, "eventWeight"])
    npmatSR = jfE.SR.RDF.AsNumpy([var, "eventWeight"])
    npmatDR_sl = jfE.DR_sl.RDF.AsNumpy([var, "eventWeight"])
    npmatDR_bl = jfE.DR_bl.RDF.AsNumpy([var, "eventWeight"])
    plt.gcf().set_size_inches((10, 6))
    plt.title(var)

    nameL = ["DR_sl", "SR", "DR_bl", "AR"]
    for i, ds in enumerate([npmatDR_sl, npmatSR, npmatDR_bl, npmatAR]):
        plt.subplot(2, 2, i + 1)
        plt.hist(ds[var], range=(0, 500), weights=ds["eventWeight"], bins=100)
        plt.yscale('log')
        plt.title(nameL[i])

    plt.tight_layout()
    plt.show()


def plotvar(jfE, var):
    plt.gcf().set_size_inches((10, 6))
    plt.title(var)
    if var in jfE.RDF.GetColumnNames():
        npmat = jfE.RDF.AsNumpy([var, "eventWeight"])
        plt.hist(npmat[var], weights=npmat["eventWeight"], bins=20)
    else:
        npmat = jfE.RDF.Define("tmpplotvar",
                               var).AsNumpy(["tmpplotvar", "eventWeight"])
        plt.hist(npmat["tmpplotvar"], weights=npmat["eventWeight"], bins=20)
    plt.yscale('log')
    plt.title(var + " from " + str(jfE.meta))
    plt.tight_layout()
    plt.show()


# %%

jfQCD = jetFakeEstimation("2016", "mt", "QCD")
jfWjet = jetFakeEstimation("2016", "mt", "W+jets")
jfttbar = jetFakeEstimation("2016", "mt", "ttbar")

# %%
lastvals = {}


class DynBins(object):
    def __init__(self, slarr, blarr, predictionVars):
        self.slarr = slarr
        self.blarr = blarr
        self.predictionVars = predictionVars
        ### histogram both regions in the generated bins
        ### use dynamic binning, from the bl region -> no 1/0 bins
        self.binbordersD = {
            var: self.dynbins(self.blarr[var], self.blarr["eventWeight"])
            for var in predictionVars
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
        print(self.slh.sum())

        ## currently, the merge of the 1d bins into a nbinning lev
        ## join bins  a var that do not have sufficient statistic due to the "multiplication" of the 1d binnings
        ## the binning should be reduced in the variable with most bins first:
        predictionVarsByLengthL = list(
            enumerate(copy.deepcopy(self.predictionVars)))
        predictionVarsByLengthL.sort(key=lambda k: len(self.binbordersD[k[1]]),
                                     reverse=True)
        print(predictionVarsByLengthL)
        print(self.predictionVars)
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
                    print(
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
        ## generate a list of index tuples for slh, blh
        # self.idxs = [
        #     tuple(x) for x in itertools.product(
        #         *[list(range(len(l))) for l in self.bincentersL])
        # ]
        ## filter indices of bins with less than 10 events
        self.idxs = list(np.ndindex(*self.slh.shape))
        self.valididxs = [
            idx for idx in self.idxs if self.slh[idx] > 5 and self.blh[idx] > 5
        ]
        print(self.slh.sum())
        ## select the bins with mor then 10 events in signal and bkg region as x values for our regression
        self.xv = np.array(list(
            itertools.product(*np.array(self.bincentersL))))[self.valididxs]
        self.yv = np.array(
            map(lambda idx: self.slh[idx] / self.blh[idx], self.valididxs))
        self.ymat = np.full(self.slh.shape, np.nan)
        for idx, y in zip(self.valididxs, self.yv):
            self.ymat[idx] = y
        import seaborn as sns
        sns.heatmap(self.ymat)

        # # x values for the fit
        # self.xsc = np.arange(
        #     self.binborders[0], self.binborders[-1] + 200,
        #     200. / np.abs(self.binborders[-1] - self.binborders[0]))

    def dynbins(self, arr, weightsarr):
        import numpy.lib.recfunctions as rfn
        # if weights == None:
        #     weightsarr = [np.random.uniform(0, .7) for i in range(len(arr))]
        #     weightsarr = [1 for i in range(len(arr))]
        varS = rfn.merge_arrays([arr, weightsarr])
        #print varS
        varS.sort(order="f0")

        marker = len(varS) - 1
        markerval = varS["f0"][marker]
        stdbinwith = 10
        binborders = [varS["f0"][-1]]
        varmin = varS["f0"][0]

        while (0 != marker):
            curbinsum = 0
            #make sure the sum of the event weights in the bin is at least 100
            while (curbinsum < 100):
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

    def calcyerror(self):
        ### Error calculation
        ## σ_y= abs(y)sqrt((∂_s y)^2*σ_s^2+(∂_b y)^2*σ_b^2+2(∂_s y)(∂_b y)σ_s^2*σ_b^2)
        # y,s,b>0 y=s/b
        # = y*sqrt((σ_s/s)^2+(σ_b/b)^2+σ_s*σ_b/(s*b))
        # poisson error
        # = y*sqrt(sqrt(1/s)^2+sqrt(1/b)^2-2*σ_sb/(s*b))
        # = y*sqrt(sqrt(1/s)^2+sqrt(1/b)^2)
        yerror = np.array([
            self.yv[i] * np.sqrt(1 / self.slh[i] + 1 / self.blh[i])
            for i in range(len(self.yv))
        ])
        return yerror

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
        self.slh=binningObj.slh
        self.blh=binningObj.blh
        ### borders for comparing the old fakefaktors

        # if jfRS.meta["bkg"]=="QCD":
        #     self.binborders=np.array([16.8, 20.0, 23.2, 26.4, 29.599999999999998, 32.800000000000004, 36.0, 39.2, 42.400000000000006, 45.60000000000001, 48.800000000000004, 52.0, 55.2, 58.400000000000006, 61.60000000000001, 64.8, 68.00000000000001, 71.2, 74.4, 77.60000000000001, 80.8, 84.00000000000001, 87.2, 90.4, 93.60000000000001, 96.80000000000001], np.float64)
        # #elif jfRS.meta["bkg"]=="W+jets": self.binborders=np.array([29,32.5,36,42,65,500], np.float32)
        # else: raise Exception

        #self.applyCorrections()
        # self.yerror = self.calcyerror()
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


class ClosureCorrection(fakefactor):
    def __init__(self, ffObject):
        ## get the numpy arrays from the jfObject
        self.jfRS = ffObject.jfRS
        self.filterstring = ffObject.filterstring
        var = "m_vis"
        self.var = var
        self.meta = ffObject.meta
        self.blarr = self.jfRS.Closure_bl.RDF.Filter(
            self.filterstring).AsNumpy([var, "eventWeight"])
        self.slarr = self.jfRS.Closure_sl.RDF.Filter(
            self.filterstring).AsNumpy([var, "eventWeight"])
        self.binborders = super(ClosureCorrection,
                                self).dynbins(self.blarr[var])


#%%
a = fakefactor(jfQCD, "1==1", ["pt_2", "m_vis"])

# l=[]
# for jf in [jfQCD,jfWjet]:
#     for filter in ["njets==0","njets>0"]:
#         l.append(fakefactor(jf,filter,"pt_2"))

# # %% ### compare old FF's
# f=ROOT.TFile("/portal/ekpbms2/home/jbechtel/fakefactors/new/CMSSW_8_0_25/src/ViennaTool/sim/mt/CR_QCD_pt_data.root")
# slh=f.Get("hh_t_pt")
# blh=f.Get("hh_l_pt")

# binborder=[slh.GetBinLowEdge(0)]
# for i in range(slh.GetNbinsX()):
#     binborder.append(slh.GetBinLowEdge(i)+slh.GetBinWidth(i))
# print binborder
# bincenters=np.array([slh.GetBinCenter(i) for i in range(slh.GetNbinsX())])
# if all(bincenters!=[blh.GetBinCenter(i) for i in range(blh.GetNbinsX())]):
#     raise Exception
# slbins=np.array([slh.GetBinContent(i) for i in range(slh.GetNbinsX())])
# blbins=np.array([blh.GetBinContent(i) for i in range(blh.GetNbinsX())])

# #plt.plot(a.binbordersToCenters(a.binborders),a.blh,label="new")
# plt.plot(bincenters,a.blh/blbins,label="new/old")
# plt.title("bl ff")
# plt.legend(loc="best")
# plt.savefig("/home/mscham/bl.png")
# plt.show()

# #plt.plot(a.binbordersToCenters(a.binborders),a.slh,label="new")
# plt.plot(bincenters,a.slh/slbins,label="new/old")
# plt.title("signal-like ff")
# plt.legend(loc="best")
# plt.savefig("/home/mscham/sl.png")
# plt.show()

# %%
