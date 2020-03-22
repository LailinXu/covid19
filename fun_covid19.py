#!/usr/bin/env python

# import ROOT
import sys
# path to root lib, can be found via `root-config  --libdir`
ROOTSYS = '/Applications/root/root_v6.16.00/lib'
sys.path.append(ROOTSYS)

import ROOT
from ROOT import gROOT, gStyle, TFile, TH1F, TCanvas, TLegend, TF1, TPad, TPaveText
from math import sqrt, fabs

from datetime import date,timedelta, datetime

## fit COVID-19 cases 
## a better model can be found here: https://www.thelancet.com/action/showPdf?pii=S0140-6736%2820%2930260-9

def read_data(infile="", dname=""):
  """
  read data and fill a histogram
  data format: 
    day1 Ncases1
    day2 Ncases2
    ...
  """

  hname="hist_"+dname
  nbin=100
  hist=TH1F(hname, hname, nbin, 0.5, nbin+0.5)
  hist.Sumw2()

  list_label=[]
  list_data=[]
  with open(infile, "r") as fin:
    for line in fin:
      line=line.rstrip()
      sline=line.split() 
      if len(sline)>=2:
        nc=int(sline[1])
        list_label.append(sline[0])
        list_data.append(nc)

  # extend data/label
  if len(list_data) < nbin:
    llast=list_label[-1]
    for ib in range(nbin):
      if ib < len(list_data): continue
      else:
        list_label.append( str(int(llast) + ib-len(list_data)+1) ) 
        list_data.append(0)

  # convert to real date
  list_date=convert_date(list_label)

  # fill hist
  for il, d in enumerate(list_data):
    hist.SetBinContent(il+1, d)
    hist.SetBinError(il+1, sqrt(d)) ## Gaussian error
    #hist.GetXaxis().SetBinLabel(il+1, list_label[il])
    hist.GetXaxis().SetBinLabel(il+1, list_date[il])

  hist.GetXaxis().LabelsOption("v")
  return hist

def convert_date(list_label=[], ref_date="1:20/02/15"):
  """
  convert index to real date
  ref_date:  the reference date, for example, "1:20/02/15": "1" means the index "1" from the raw input data; and "20/02/15" is the
             corresponding actual date
  """

  ref_ind=ref_date.split(":")[0]
  ref_d=ref_date.split(":")[1]

  datetime_object = datetime.strptime(ref_d, '%y/%m/%d')

  list_date=[]
  for il, la in enumerate(list_label):
    diff=int(la)-int(ref_ind)

    rdate= (datetime_object + timedelta(days=diff)).date().strftime('%y-%m-%d')
    
    list_date.append(rdate)


  return list_date

def get_diff(hist=None, func=None):
  """
  compare the predicted, func, with the observed, hist
  and get the difference
  """
  
  hres=hist.Clone()
  hname=hist.GetName()
  hres.SetName(hname+"_diff")
  hres.SetTitle(hname+"_diff")

  nb=hist.GetNbinsX()
  for ib in range(nb):
    y=hist.GetBinContent(ib+1)
    yerr=hist.GetBinError(ib+1)
    if y==0: continue

    yp=func.Eval(ib+1.)

    diff=y-yp

    hres.SetBinContent(ib+1, diff)
    hres.SetBinError(ib+1, yerr)

  return hres

def get_pred(hist=None, func=None):
  """
  get the prediction from the fitted function
  """
  
  hres=hist.Clone()
  hname=hist.GetName()
  hres.SetName(hname+"_pred")
  hres.SetTitle(hname+"_pred")

  nb=hist.GetNbinsX()
  for ib in range(nb):
    y=hist.GetBinContent(ib+1)
    yerr=hist.GetBinError(ib+1)
    if y>0: continue

    yp=func.Eval(ib+1.)
    yerr=sqrt(yp)

    hres.SetBinContent(ib+1, int(yp))
    hres.SetBinError(ib+1, yerr)

  return hres


def fit(infile="", label="", output=""):

  tfout=TFile.Open("%s.root" % output, "RECREATE")

  ## read data and fill a histogram
  hist=read_data(infile)

  xoffset=0.5
  xmin=0+xoffset
  xmax=45+xoffset
  ymin=0
  ymax=hist.GetMaximum()
 
  ## do fit
  fun1=TF1("fun1", "expo", xmin, xmax)
  fun1.SetLineColor(ROOT.kRed)
  hist.Fit("fun1")

  fun2=TF1("fun2", "[0]*exp([1]*x+[2]*1./x)", xmin, xmax)
  #fun2=TF1("fun2", "[0]*exp([1]*x+[2]*1./x+[3]*x*x)", xmin, xmax)
  fun2.SetLineColor(ROOT.kBlue)
  fun2.SetParameter(0, fun1.GetParameter(0))
  fun2.SetParameter(1, fun1.GetParameter(1))
  hist.Fit("fun2")

  fun2.SetParameter(0, fun2.GetParameter(0))
  fun2.SetParameter(1, fun2.GetParameter(1))
  fun2.SetParameter(2, fun2.GetParameter(2))
  #fun2.SetParameter(3, fun2.GetParameter(3))
  hist.Fit("fun2")

  prob=fun2.GetProb()
  chi2=fun2.GetChisquare()
  print "prob:", prob
  print "chi2:", chi2
  
  gROOT.SetStyle("ATLAS")

  yscale=2.5

  doRatio=1
  ## plot
  canvas_w,canvas_h=800,600
  fraction=1.
  ratio=0.2
  if doRatio:
    canvas_w,canvas_h=800,int(600*(1+ratio))
    fraction=ratio+0.2

  gStyle.SetPadLeftMargin(0.15)
  gStyle.SetPadRightMargin(0.10)
  gStyle.SetPadTopMargin(0.05)
  gStyle.SetPadBottomMargin(0.15)
  MyC=TCanvas('MyC','MyC',canvas_w,canvas_h)
  MyC.SetTicks(1,1)
  ymax*=yscale

  if doRatio:
    Pad1 = TPad("p1","p1",0,fraction*1.0/(fraction+1),1,1,0,0) # x1,y1,x2,y2
    Pad1.SetMargin(0.15,0.10,0.03,0.05)
    Pad2 = TPad("p2","p2",0,0,1,fraction*1.0/(fraction+1),0,0)
    Pad2.SetMargin(0.15,0.10,0.15/fraction,0.08)
    Pad1.Draw()
    Pad2.Draw()

  ## draw data and fitted lines
  hist.SetMarkerStyle(20)
  hist.SetMarkerSize(1.8)
  hist.SetLineWidth(2)
  hist.GetXaxis().SetTitle("Days")
  hist.GetYaxis().SetTitle("Cases")
  hist.GetYaxis().SetTitleOffset(1.1*hist.GetYaxis().GetTitleOffset())
  hist.GetXaxis().SetTitleOffset(1.1*hist.GetXaxis().GetTitleOffset())
  hist.GetXaxis().SetRangeUser(xmin, xmax)
  hist.SetMaximum(ymax)

  if doRatio:
    Pad1.cd()
  else:
    MyC.cd()

  gStyle.SetPaintTextFormat("      g")
  hist.Draw("petext")
  
  fun1.Draw("same")
  fun2.Draw("same")

  ## plot predictions
  hpred=get_pred(hist, fun2)
  hpred.SetMarkerColor(ROOT.kBlue)
  hpred.Draw("petextsame")
  hist.Draw("petextsame")

  ## Legend
  l = TLegend(0.20, 0.67, 0.35, 0.92)
  l.SetFillColor(10)
  l.SetBorderSize(0)
  l.SetTextSize(0.04)

  l.AddEntry(hist, label, "p")
  l.AddEntry(hpred, "Predicted", "p")
  flabel1="Fit: %.2fexp( %.2f*t)" % (fun1.GetParameter(0), fun1.GetParameter(1))
  l.AddEntry(fun1, flabel1, "l")
  flabel2="Fit: %.2fexp( %.2f*t + %.2f/t)" % (fun2.GetParameter(0), fun2.GetParameter(1), fun2.GetParameter(2))
  #flabel2="Fit: %.2fexp( %.2f*t + %.2f/t + %.2f/t^{2})" % (fun2.GetParameter(0), fun2.GetParameter(1), fun2.GetParameter(2), fun2.GetParameter(3))
  l.AddEntry(fun2, flabel2, "l")
  
  l.Draw()

  ## show fit quality
  """
  tl=TPaveText(0.30, 0.57, 0.50, 0.62, "NDC")
  tl.SetFillColor(10)
  tl.SetBorderSize(0)
  tl.SetTextSize(0.04)
  tl.AddText("p-value: %.3f" % prob)
  tl.Draw()
  """

  ## compare predicted and observed
  if doRatio:
    hres=get_diff(hist, fun2)
    hres.GetXaxis().SetTitle(hist.GetXaxis().GetTitle())
    hres.GetYaxis().SetTitle("Data-Pred.")
    hres.GetXaxis().SetTitleOffset(1.0)
    hres.GetYaxis().SetTitleOffset(0.4)
    hres.GetXaxis().SetLabelOffset( hres.GetXaxis().GetLabelOffset()*2.0 )
    hres.GetXaxis().SetTitleSize(hist.GetXaxis().GetTitleSize()*2.5)
    hres.GetYaxis().SetTitleSize(hist.GetYaxis().GetTitleSize()*2.5)
    hres.GetXaxis().SetLabelSize(hist.GetXaxis().GetLabelSize()*2.5)
    hres.GetYaxis().SetLabelSize(hist.GetYaxis().GetLabelSize()*2)
    hres.GetXaxis().SetTitleOffset(2.*hist.GetXaxis().GetTitleOffset())

    hist.GetXaxis().SetTitle("")
    hist.GetXaxis().SetLabelSize(0)

    Pad2.cd()
    hres.GetFunction("fun2").SetBit(ROOT.TF1.kNotDraw)
    dymin=hres.GetMinimum()
    dymax=hres.GetMaximum()
    dy=max(fabs(dymin), fabs(dymax))*1.1
    dy=min(dy, 2000)
    hres.GetYaxis().SetRangeUser(-dy, dy)
    hres.Draw()

  MyC.SaveAs("%s_linear.png" % output)
  hist.SetMinimum(11)
  hist.SetMaximum( hist.GetMaximum() * 2.)
  if doRatio:
    Pad1.SetLogy(1)
  else:
    MyC.SetLogy(1)
  MyC.SaveAs("%s_log.png" % output)
  
  tfout.cd()
  hist.Write()
  fun1.Write()
  fun2.Write()
  tfout.Close()

def fit_IT():

  infile="data/IT_20200320.txt"
  label="Italy"
  output="data/"+label

  fit(infile=infile, label=label, output=output)


if __name__ == "__main__":

  fit_IT()
