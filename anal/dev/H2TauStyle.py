from ROOT import TColor, kViolet, kBlue, kRed, kMagenta

class Style:

    def __init__(self,
                 markerStyle=8,
                 markerColor=1,
                 markerSize=1,
                 lineStyle=1,
                 lineColor=1,
                 lineWidth=2,
                 fillColor=None,
                 fillStyle=1001):
        self.markerStyle = markerStyle
        self.markerColor = markerColor
        self.markerSize = markerSize
        self.lineStyle = lineStyle
        self.lineColor = lineColor
        self.lineWidth = lineWidth
        if fillColor is None:
            self.fillColor = lineColor
        else:
            self.fillColor = fillColor
        self.fillStyle = fillStyle

    def formatHisto(self, hist, title=None):
        hist.SetMarkerStyle(self.markerStyle)
        hist.SetMarkerColor(self.markerColor)
        hist.SetMarkerSize(self.markerSize)
        hist.SetLineStyle(self.lineStyle)
        hist.SetLineColor(self.lineColor)
        hist.SetLineWidth(self.lineWidth)
        hist.SetFillColor(self.fillColor)
        hist.SetFillStyle(self.fillStyle)
        if title != None:
            hist.SetTitle(title)
        return hist

# the following standard files are defined and ready to be used.
# more standard styles can be added on demand.
# user defined styles can be created in the same way in any python module

sBlack = Style()
sData = Style(fillStyle=0, markerSize=1.3)
sBlue = Style(markerColor=4, fillColor=4)
sGreen = Style(markerColor=8, fillColor=8)
sRed = Style(markerColor=2, fillColor=2)
sYellow = Style(lineColor=1, markerColor=5, fillColor=5)
sViolet = Style(lineColor=1, markerColor=kViolet, fillColor=kViolet)

qcdcol = TColor.GetColor(250,202,255)
sHTT_QCD = Style(lineColor=1, markerColor=qcdcol, fillColor=qcdcol)
#dycol =  TColor.GetColor(248,206,104)
dycol =  TColor.GetColor(248,186,102)
orange_1 =  TColor.GetColor(255,204,153)
sHTT_DYJets = Style(lineColor=1, markerColor=dycol, fillColor=dycol)
wcol = TColor.GetColor(222,90,106)
sHTT_WJets = Style(lineColor=1, markerColor=wcol, fillColor=wcol)

taucol = TColor.GetColor(51,102,255)
s_taucol = Style(lineColor=1, markerColor=taucol, fillColor=taucol)
taucol_1 = TColor.GetColor(153,153,204)
s_taucol_1 = Style(lineColor=1, markerColor=taucol_1, fillColor=taucol_1)
taucol_2 = TColor.GetColor(153,153,255)
s_taucol_2 = Style(lineColor=1, markerColor=taucol_2, fillColor=taucol_2)
taucol_3 = TColor.GetColor(204,204,255) 
s_taucol_3 = Style(lineColor=1, markerColor=taucol_3, fillColor=taucol_3)
taucol_4 = TColor.GetColor(51,102,153)     
s_taucol_4 = Style(lineColor=1, markerColor=taucol_4, fillColor=taucol_4)

bc_mu_1 = TColor.GetColor(250,202,255)
s_bc_mu_1=Style(lineColor=1, markerColor=bc_mu_1, fillColor=bc_mu_1)
bc_mu_2 = TColor.GetColor(219,166,225)
s_bc_mu_2=Style(lineColor=1, markerColor=bc_mu_2, fillColor=bc_mu_2)
bc_mu_3 = TColor.GetColor(232,135,243)
s_bc_mu_3=Style(lineColor=1, markerColor=bc_mu_3, fillColor=bc_mu_3)
bc_mu_4 = TColor.GetColor(243,217,246)
s_bc_mu_4=Style(lineColor=1, markerColor=bc_mu_4, fillColor=bc_mu_4)
bc_mu_5 = TColor.GetColor(203,83,216)
s_bc_mu_5=Style(lineColor=1, markerColor=bc_mu_5, fillColor=bc_mu_5)
jpsi_1 = TColor.GetColor(81,213,218)
s_jpsi_1=Style(lineColor=1, markerColor=jpsi_1, fillColor=jpsi_1)
jpsi_2 = TColor.GetColor(32,115,118)
s_jpsi_2=Style(lineColor=1, markerColor=jpsi_2, fillColor=jpsi_2)
jpsi_3 = TColor.GetColor(175,231,233)
s_jpsi_3=Style(lineColor=1, markerColor=jpsi_3, fillColor=jpsi_3)
jpsi_4 = TColor.GetColor(124,187,189)
s_jpsi_4=Style(lineColor=1, markerColor=jpsi_4, fillColor=jpsi_4)
jpsi_5 = TColor.GetColor(1,127,131)
s_jpsi_5=Style(lineColor=1, markerColor=jpsi_5, fillColor=jpsi_5)
others = TColor.GetColor(204,255,204)
s_others = Style(lineColor=1, markerColor=others, fillColor=others)
tau_1 = TColor.GetColor(247,243,126)
s_tau_1 = Style(lineColor=1, markerColor=tau_1,fillColor=tau_1)

lob_1 = TColor.GetColor(255,154,154)
s_lob_1=Style(lineColor=1, markerColor=lob_1, fillColor=lob_1)
pink_1 = TColor.GetColor(255,204,255)                          
s_pink_1=Style(lineColor=1, markerColor=pink_1, fillColor=pink_1)   

ttcol = TColor.GetColor(155,152,204)
sHTT_TTJets = Style(lineColor=1, markerColor=ttcol, fillColor=ttcol)

ttcol_v2 = TColor.GetColor(135,206,250)
sHTT_TTJets_v2 = Style(lineColor=1, markerColor=ttcol_v2, fillColor=ttcol_v2)

sHTT_Higgs_1 = Style(lineColor=2, markerSize=0, markerColor=2, lineWidth=3, lineStyle=1, fillColor=0)
sHTT_Higgs_2 = Style(lineColor=3, markerSize=0, markerColor=4, lineWidth=3, lineStyle=1, fillColor=0)
sHTT_Higgs_3 = Style(lineColor=4, markerSize=0, markerColor=4, lineWidth=3, lineStyle=1, fillColor=0)
sHTT_Higgs_4 = Style(lineColor=6, markerSize=0, markerColor=6, lineWidth=3, lineStyle=1, fillColor=0)
sHTT_Higgs_5 = Style(lineColor=7, markerSize=0, markerColor=7, lineWidth=3, lineStyle=1, fillColor=0)
sHTT_Higgs_6 = Style(lineColor=8, markerSize=0, markerColor=8, lineWidth=3, lineStyle=1, fillColor=0)
sHTT_Higgs_7 = Style(lineColor=9, markerSize=0, markerColor=9, lineWidth=3, lineStyle=1, fillColor=0)
sHTT_Higgs_8 = Style(lineColor=11, markerSize=0, markerColor=11, lineWidth=3, lineStyle=1, fillColor=0)
sHTT_Higgs_9 = Style(lineColor=12, markerSize=0, markerColor=12, lineWidth=3, lineStyle=1, fillColor=0)
sHTT_Higgs_10 = Style(lineColor=46, markerSize=0, markerColor=46, lineWidth=3, lineStyle=1, fillColor=0)
sHTT_Higgs_11 = Style(lineColor=47, markerSize=0, markerColor=47, lineWidth=3, lineStyle=1, fillColor=0)
sHTT_Higgs_12 = Style(lineColor=48, markerSize=0, markerColor=48, lineWidth=3, lineStyle=1, fillColor=0)


tHTT_Higgs_1 = Style(lineColor=2, markerSize=0, markerColor=2, lineWidth=3, lineStyle=1, fillColor=0)
tHTT_Higgs_2 = Style(lineColor=3, markerSize=0, markerColor=4, lineWidth=3, lineStyle=1, fillColor=0)
tHTT_Higgs_3 = Style(lineColor=4, markerSize=0, markerColor=4, lineWidth=3, lineStyle=1, fillColor=0)
tHTT_Higgs_4 = Style(lineColor=6, markerSize=0, markerColor=6, lineWidth=3, lineStyle=1, fillColor=0)
tHTT_Higgs_5 = Style(lineColor=7, markerSize=0, markerColor=7, lineWidth=3, lineStyle=1, fillColor=0)
tHTT_Higgs_6 = Style(lineColor=8, markerSize=0, markerColor=8, lineWidth=3, lineStyle=1, fillColor=0)
tHTT_Higgs_7 = Style(lineColor=9, markerSize=0, markerColor=9, lineWidth=3, lineStyle=1, fillColor=0)
tHTT_Higgs_8 = Style(lineColor=11, markerSize=0, markerColor=11, lineWidth=3, lineStyle=1, fillColor=0)
tHTT_Higgs_9 = Style(lineColor=12, markerSize=0, markerColor=12, lineWidth=3, lineStyle=1, fillColor=0)
tHTT_Higgs_10 = Style(lineColor=46, markerSize=0, markerColor=46, lineWidth=3, lineStyle=1, fillColor=0)
tHTT_Higgs_11 = Style(lineColor=47, markerSize=0, markerColor=47, lineWidth=3, lineStyle=1, fillColor=0)
tHTT_Higgs_12 = Style(lineColor=48, markerSize=0, markerColor=48, lineWidth=3, lineStyle=1, fillColor=0)

pHTT_Higgs_1 = Style(lineColor=2, markerSize=0, markerColor=2, lineWidth=3, lineStyle=1, fillColor=0)
pHTT_Higgs_2 = Style(lineColor=3, markerSize=0, markerColor=4, lineWidth=3, lineStyle=1, fillColor=0)
pHTT_Higgs_3 = Style(lineColor=2, markerSize=0, markerColor=2, lineWidth=3, lineStyle=1, fillColor=0)
pHTT_Higgs_4 = Style(lineColor=6, markerSize=0, markerColor=6, lineWidth=5, lineStyle=1, fillColor=0)
pHTT_Higgs_5 = Style(lineColor=7, markerSize=0, markerColor=7, lineWidth=3, lineStyle=1, fillColor=0)
pHTT_Higgs_6 = Style(lineColor=8, markerSize=0, markerColor=8, lineWidth=3, lineStyle=1, fillColor=0)
pHTT_Higgs_7 = Style(lineColor=9, markerSize=0, markerColor=9, lineWidth=3, lineStyle=1, fillColor=0)
pHTT_Higgs_8 = Style(lineColor=11, markerSize=0, markerColor=11, lineWidth=3, lineStyle=1, fillColor=0)
pHTT_Higgs_9 = Style(lineColor=12, markerSize=0, markerColor=12, lineWidth=3, lineStyle=1, fillColor=0)
pHTT_Higgs_10 = Style(lineColor=46, markerSize=0, markerColor=46, lineWidth=3, lineStyle=1, fillColor=0)
pHTT_Higgs_11 = Style(lineColor=47, markerSize=0, markerColor=47, lineWidth=3, lineStyle=1, fillColor=0)
pHTT_Higgs_12 = Style(lineColor=48, markerSize=0, markerColor=48, lineWidth=3, lineStyle=1, fillColor=0)



sHTT_Higgs = Style(lineColor=kBlue, markerColor=2, lineStyle=2, fillColor=0)
zlcol = TColor.GetColor(100,182,232)
sHTT_ZL = Style(lineColor=1, markerColor=zlcol, fillColor=zlcol)
dibosoncol = TColor.GetColor(222,140,106)
sHTT_VV = Style(lineColor=1, markerColor=dibosoncol, fillColor=dibosoncol)

jtfake = TColor.GetColor(100,222,106)
sHTT_jtfake = Style(lineColor=1, markerColor=jtfake, fillColor=jtfake)

lowmass = TColor.GetColor(240,175,60)
sHTT_lowmass = Style(lineColor=1, markerColor=lowmass, fillColor=lowmass)


sBlackSquares = Style(markerStyle=21)
sBlueSquares = Style(lineColor=4, markerStyle=21, markerColor=4)
sGreenSquares = Style(lineColor=8, markerStyle=21, markerColor=8)
sRedSquares = Style(lineColor=2, markerStyle=21, markerColor=2)


styleSet = [sBlue, sGreen, sRed, sYellow, sViolet, sBlackSquares, sBlueSquares, sGreenSquares, sRedSquares]
iStyle = 0


def nextStyle():
    global iStyle
    style = styleSet[iStyle]
    iStyle = iStyle+1
    if iStyle >= len(styleSet):
        iStyle = 0
    return style

histPref = {}
histPref['Data'] = {'style':sData, 'layer':2999, 'legend':'Observed'}
histPref['data_*'] = {'style':sData, 'layer':2999, 'legend':'Observed'}
#histPref['ZTT'] = {'style':sHTT_DYJets, 'layer':4, 'legend':'Z#rightarrow#tau#tau'}
#histPref['DY50'] = {'style':sHTT_DYJets, 'layer':4, 'legend':'Z#rightarrow#tau#tau (M > 50)'}
#histPref['DY10to50'] = {'style':sHTT_ZL, 'layer':5, 'legend':'Z#rightarrow#tau#tau (10 < M < 50)'}
histPref['bg_ul*'] = {'style':sHTT_DYJets, 'layer':10, 'legend':'J/#psi X bkg'}
histPref['bg_bc'] = {'style':sHTT_QCD, 'layer':1, 'legend':'Bc BG'}
histPref['bg_norm'] = {'style':sHTT_jtfake, 'layer':1, 'legend':'Bg norm.'}
histPref['sig_had'] = {'style':sHTT_TTJets_v2, 'layer':4, 'legend':'Signal (3prong)'}
#ALL possible Bc decays
histPref['bc_jpsi_tau_3p'] = {'style':s_taucol_3, 'layer':4, 'legend':'Bc#rightarrowJ/#psi#tau_{h}#nu'}
histPref['bc_jpsi_tau_N3p'] = {'style':s_lob_1, 'layer':4, 'legend':'Bc#rightarrowJ/#psi#tau_{oth}#nu'}
histPref['bc_jpsi_tau_mu'] = {'style':s_taucol_1, 'layer':4, 'legend':'Bc#rightarrowJ/#psi#tau_{#mu}#nu'}
histPref['bc_jpsi_mu'] = {'style':s_pink_1, 'layer':4, 'legend':'Bc#rightarrowJ/#psi#mu#nu'}
histPref['bc_charmonium_mu'] = {'style':s_bc_mu_1, 'layer':4, 'legend':'Bc#rightarrow(c#bar{c})_{other}#mu#nu'}
histPref['bc_psi2s_mu'] = {'style':s_bc_mu_1, 'layer':4, 'legend':'Bc#rightarrow#psi(2s)#mu#nu'}
histPref['bc_chic0_mu'] = {'style':s_bc_mu_2, 'layer':4, 'legend':'Bc#rightarrow#chi_{c0}#mu#nu'}
histPref['bc_chic1_mu'] = {'style':s_bc_mu_3, 'layer':4, 'legend':'Bc#rightarrow#chi_{c1}#mu#nu'}
histPref['bc_chic2_mu'] = {'style':s_bc_mu_4, 'layer':4, 'legend':'Bc#rightarrow#chi_{c2}#mu#nu'}
histPref['bc_hc_mu'] = {'style':s_bc_mu_5, 'layer':4, 'legend':'Bc#rightarrowh_{c}#mu#nu'}
histPref['bc_psi2s_tau'] = {'style':s_tau_1, 'layer':4, 'legend':'Bc#rightarrow#psi(2s)#tau#nu'}
histPref['bc_chic0_tau'] = {'style':s_bc_mu_2, 'layer':4, 'legend':'Bc#rightarrow#chi_{c0}#tau#nu'}
histPref['bc_chic1_tau'] = {'style':s_bc_mu_3, 'layer':4, 'legend':'Bc#rightarrow#chi_{c1}#tau#nu'}
histPref['bc_chic2_tau'] = {'style':s_bc_mu_4, 'layer':4, 'legend':'Bc#rightarrow#chi_{c2}#tau#nu'}
histPref['bc_hc_tau'] = {'style':s_bc_mu_5, 'layer':4, 'legend':'Bc#rightarrowh_{c}#tau#nu'}
histPref['bc_jpsi_pi'] = {'style':s_jpsi_1, 'layer':4, 'legend':'Bc#rightarrowJ/#psi#pi'}
histPref['bc_jpsi_3pi'] = {'style':s_jpsi_2, 'layer':4, 'legend':'Bc#rightarrowJ/#psi3#pi'}
histPref['bc_jpsi_ds'] = {'style':s_jpsi_3, 'layer':4, 'legend':'Bc#rightarrowJ/#psi D^{(*)+}_{s}'}
histPref['bc_jpsi_5pi'] = {'style':s_jpsi_4, 'layer':4, 'legend':'Bc#rightarrowJ/#psi5#pi'}
histPref['bc_jpsi_pions'] = {'style':s_lob_1, 'layer':4, 'legend':'Bc#rightarrowJ/#psi pions'}
histPref['bc_others'] = {'style':s_others, 'layer':4, 'legend':'Bc#rightarrow others'}

histPref['sig_3pp*'] = {'style':sHTT_TTJets, 'layer':3, 'legend':'Signal (3prong + #pi^{0})'}
histPref['sig_others*'] = {'style':sHTT_WJets, 'layer':2, 'legend':'Signal (others)'} 
histPref['signal*'] = {'style':sHTT_TTJets_v2, 'layer':1, 'legend':'Signal'} 
histPref['TTT'] = {'style':sHTT_TTJets, 'layer':10, 'legend':'t#bar{t} with real #tau_{h}'} 
histPref['TTJ'] = {'style':sHTT_TTJets, 'layer':10, 'legend':'t#bar{t} with jets'} 
histPref['T*tW*'] = {'style':sHTT_TTJets, 'layer':1, 'legend':'Single t'}
histPref['electroweak*'] = {'style':sHTT_WJets, 'layer':0.2, 'legend':'Electroweak'} 
histPref['EWK'] = {'style':sHTT_DYJets, 'layer':0.2, 'legend':'Electroweak'} 
histPref['TTo*'] = {'style':sHTT_TTJets, 'layer':1, 'legend':'Single t'} 
histPref['Single t'] = {'style':sHTT_TTJets, 'layer':1, 'legend':'Single t'} 
histPref['ST*'] = {'style':sHTT_TTJets, 'layer':1, 'legend':'Single top quark'} 
histPref['WW*'] = {'style':sHTT_VV, 'layer':0.9, 'legend':'Diboson'} 
histPref['WZ*'] = {'style':sHTT_VV, 'layer':0.8, 'legend':'Diboson'} 
histPref['ZZ*'] = {'style':sHTT_VV, 'layer':0.7, 'legend':'Diboson'} 
histPref['Diboson'] = {'style':sHTT_VV, 'layer':0.7, 'legend':'Diboson'} 
histPref['VV*'] = {'style':sHTT_VV, 'layer':0.7, 'legend':'Diboson'} 
histPref['QCD_'] = {'style':sHTT_QCD, 'layer':0.1, 'legend':'Jet#rightarrow#tau_{h} fakes'}
histPref['QCD'] = {'style':sHTT_QCD, 'layer':0.1, 'legend':'QCD multijet'}
histPref['W'] = {'style':sHTT_WJets, 'layer':3, 'legend':'W+jets'}  
histPref['Electroweak'] = {'style':sHTT_DYJets, 'layer':3, 'legend':'Electroweak'}  
histPref['W1'] = {'style':sHTT_WJets, 'layer':3, 'legend':'W+jets'}  
histPref['W2'] = {'style':sHTT_WJets, 'layer':3, 'legend':'W+jets'}  
histPref['W3'] = {'style':sHTT_WJets, 'layer':3, 'legend':'W+jets'}  
histPref['W4'] = {'style':sHTT_WJets, 'layer':3, 'legend':'W+jets'}  
histPref['WJets'] = {'style':sHTT_WJets, 'layer':3, 'legend':'W+jets'}  
histPref['ZJ*'] = {'style':sHTT_jtfake, 'layer':3.1, 'legend':'Z#rightarrow#tau#tau/Z#rightarrow ll, j#rightarrow#tau'}
histPref['ZL*'] = {'style':sHTT_ZL, 'layer':3.2, 'legend':'Z#rightarrow ll'}
histPref['ZLL'] = {'style':sHTT_ZL, 'layer':3.2, 'legend':'Z#rightarrow ll'}
histPref['Ztt_TL'] = {'style':sViolet, 'layer':4.1, 'legend':'Z#rightarrow#tau#tau/Z#rightarrow ll, j#rightarrow#tau'}
histPref['Higgs*'] = {'style':sHTT_Higgs, 'layer':1001, 'legend':None}
#histPref['Signal*'] = {'style':sHTT_Higgs, 'layer':1001, 'legend':None}
#histPref['TotalSig*'] = {'style':sHTT_Higgs, 'layer':1001, 'legend':'SM H(125)'}

histPref['Signal_M200'] = {'style':sHTT_Higgs_1, 'layer':1001, 'legend':'LQ s-channel: 200 GeV'}
histPref['Signal_M300'] = {'style':sHTT_Higgs_1, 'layer':1001, 'legend':'LQ s-channel: 300 GeV'}
histPref['Signal_M400'] = {'style':sHTT_Higgs_1, 'layer':1001, 'legend':'LQ (400 GeV)'}
histPref['Signal_M500'] = {'style':sHTT_Higgs_12, 'layer':1001, 'legend':'LQ s-channel: 500 GeV'}
histPref['Signal_M600'] = {'style':sHTT_Higgs_3, 'layer':1001, 'legend':'LQ s-channel: 600 GeV'}
histPref['Signal_M700'] = {'style':sHTT_Higgs_4, 'layer':1001, 'legend':'LQ s-channel: 700 GeV'}
histPref['Signal_M800'] = {'style':pHTT_Higgs_3, 'layer':1001, 'legend':'LQ (800 GeV)'}
#histPref['TotalSig'] = {'style':pHTT_Higgs_3, 'layer':1001, 'legend':'#splitline{scalar LQ}{(700 GeV, #lambda = 1, #beta = 1)}'}
histPref['TotalSig'] = {'style':pHTT_Higgs_3, 'layer':1001, 'legend':'Scalar LQ 700 GeV'}
histPref['Signal_M900'] = {'style':sHTT_Higgs_6, 'layer':1001, 'legend':'LQ s-channel: 900 GeV'}
histPref['Signal_M1000'] = {'style':sHTT_Higgs_7, 'layer':1001, 'legend':'LQ s-channel: 1000 GeV'}
histPref['Signal_M1200'] = {'style':sHTT_Higgs_8, 'layer':1001, 'legend':'LQ s-channel: 1200 GeV'}
histPref['Signal_M1500'] = {'style':sHTT_Higgs_9, 'layer':1001, 'legend':'LQ s-channel: 1500 GeV'}
histPref['Signal_M2000'] = {'style':sHTT_Higgs_10, 'layer':1001, 'legend':'LQ s-channel: 2000 GeV'}

histPref['Signal_t_M200'] = {'style':tHTT_Higgs_1, 'layer':1001, 'legend':'LQ t-channel: 200 GeV'}
histPref['Signal_t_M300'] = {'style':tHTT_Higgs_2, 'layer':1001, 'legend':'LQ t-channel: 300 GeV'}
histPref['Signal_t_M400'] = {'style':tHTT_Higgs_2, 'layer':1001, 'legend':'LQ t-channel: 400 GeV'}
histPref['Signal_t_M500'] = {'style':tHTT_Higgs_12, 'layer':1001, 'legend':'LQ t-channel: 500 GeV'}
histPref['Signal_t_M600'] = {'style':tHTT_Higgs_3, 'layer':1001, 'legend':'LQ t-channel: 600 GeV'}
histPref['Signal_t_M700'] = {'style':tHTT_Higgs_4, 'layer':1001, 'legend':'LQ t-channel: 700 GeV'}
histPref['Signal_t_M800'] = {'style':tHTT_Higgs_5, 'layer':1001, 'legend':'LQ t-channel: 800 GeV'}
histPref['Signal_t_M900'] = {'style':tHTT_Higgs_6, 'layer':1001, 'legend':'LQ t-channel: 900 GeV'}
histPref['Signal_t_M1000'] = {'style':tHTT_Higgs_7, 'layer':1001, 'legend':'LQ t-channel: 1000 GeV'}
histPref['Signal_t_M1200'] = {'style':tHTT_Higgs_8, 'layer':1001, 'legend':'LQ t-channel: 1200 GeV'}
histPref['Signal_t_M1500'] = {'style':tHTT_Higgs_9, 'layer':1001, 'legend':'LQ t-channel: 1500 GeV'}
histPref['Signal_t_M2000'] = {'style':tHTT_Higgs_10, 'layer':1001, 'legend':'LQ t-channel: 2000 GeV'}


histPref['Signal_pair_M200'] = {'style':pHTT_Higgs_1, 'layer':1001, 'legend':'LQ pair: 200 GeV'}
histPref['Signal_pair_M300'] = {'style':pHTT_Higgs_3, 'layer':1001, 'legend':'LQ pair: 300 GeV'}
histPref['Signal_pair_M400'] = {'style':pHTT_Higgs_3, 'layer':1001, 'legend':'LQ pair: 400 GeV'}
histPref['Signal_pair_M500'] = {'style':pHTT_Higgs_12, 'layer':1001, 'legend':'LQ pair: 500 GeV'}
histPref['Signal_pair_M600'] = {'style':pHTT_Higgs_3, 'layer':1001, 'legend':'LQ pair: 600 GeV'}
histPref['Signal_pair_M700'] = {'style':pHTT_Higgs_4, 'layer':1001, 'legend':'LQ pair: 700 GeV'}
histPref['Signal_pair_M800'] = {'style':pHTT_Higgs_5, 'layer':1001, 'legend':'LQ pair: 800 GeV'}
histPref['Signal_pair_M900'] = {'style':pHTT_Higgs_6, 'layer':1001, 'legend':'LQ pair: 900 GeV'}
histPref['Signal_pair_M1000'] = {'style':pHTT_Higgs_7, 'layer':1001, 'legend':'LQ pair: 1000 GeV'}
histPref['Signal_pair_M1200'] = {'style':pHTT_Higgs_8, 'layer':1001, 'legend':'LQ pair: 1200 GeV'}
histPref['Signal_pair_M1500'] = {'style':pHTT_Higgs_9, 'layer':1001, 'legend':'LQ pair: 1500 GeV'}
histPref['Signal_pair_M2000'] = {'style':pHTT_Higgs_10, 'layer':1001, 'legend':'LQ pair: 2000 GeV'}




