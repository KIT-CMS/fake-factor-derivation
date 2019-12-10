#!python
# -*- coding: UTF-8 -*-
#%%
import ROOT
ROOT.ROOT.EnableImplicitMT(15)
from shape_producer.channel import *
from shape_producer.cutstring import Cut, Cuts, Weights, Weight
from shape_producer.era import Run2016, Run2017, Run2018
import copy,os

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
    "2018":Run2018(database)
}

from shape_producer.estimation_methods_2017 import DataEstimation, ZTTEstimation, ZJEstimation, ZLEstimation, TTLEstimation, TTJEstimation, TTTEstimation, VVTEstimation, VVJEstimation, VVLEstimation, WEstimation, ggHEstimation, qqHEstimation, EWKZEstimation, ZTTEmbeddedEstimation, NewFakeEstimationTT, NewFakeEstimationLT
from fake_factor_derivation.cuts import cutDB

class ParSpaceRegion(object):
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
        logger.debug(self.CutString)
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
        self.SR = ParSpaceCrop(self.channel,
                               self.meta,
                               self.RDF,
                               signalLikeRegion=True,
                               determinationRegion=False)

        self.AR = ParSpaceCrop(self.channel,
                               self.meta,
                               self.RDF,
                               signalLikeRegion=False,
                               determinationRegion=False)

        ## Define the Background and Signal like Determination Region
        if self.meta["bkg"] != "ttbar":
            self.DR_sl = ParSpaceCrop(self.channel,
                                      self.meta,
                                      self.RDF,
                                      signalLikeRegion=True,
                                      determinationRegion=True)
            self.DR_bl = ParSpaceCrop(self.channel,
                                      self.meta,
                                      self.RDF,
                                      signalLikeRegion=False,
                                      determinationRegion=True)
        else:
            self.DR_sl = ParSpaceCrop(self.channel,
                                      self.meta,
                                      self.ttbarRDF,
                                      signalLikeRegion=True,
                                      determinationRegion=True)
            self.DR_bl = ParSpaceCrop(self.channel,
                                      self.meta,
                                      self.ttbarRDF,
                                      signalLikeRegion=False,
                                      determinationRegion=True)

        ## Define RDFs for the closure correction
        self.Closure_sl = ParSpaceCrop(self.channel,
                                       self.meta,
                                       self.ClosureRDF,
                                       signalLikeRegion=True,
                                       determinationRegion=True)
        self.Closure_bl = ParSpaceCrop(self.channel,
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


class ParSpaceCrop(ParSpaceRegion):
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
