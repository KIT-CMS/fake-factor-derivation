#!python
# -*- coding: UTF-8 -*-
#%%
from shape_producer.cutstring import Cut, Cuts, Weights, Weight
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
