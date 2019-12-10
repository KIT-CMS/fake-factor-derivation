# %%
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