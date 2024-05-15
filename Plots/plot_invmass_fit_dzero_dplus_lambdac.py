#!/usr/bin/env python3

# Copyright 2019-2020 CERN and copyright holders of ALICE O2.
# See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
# All rights not expressly granted are reserved.

# This software is distributed under the terms of the GNU General Public
# License v3 (GPL Version 3), copied verbatim in the file "COPYING".

# In applying this license CERN does not waive the privileges and immunities
# granted to it by virtue of its status as an Intergovernmental Organization
# or submit itself to any jurisdiction.

"""
file: plot_invmass_fit_dzero_dplus_lambdac.py
brief: script to produce invariant mass fit plot for article, the CFG file should be the one used for inv mass fit
usage: python3 plot_invmass_fit_dzero_dplus_lambdac.py CFG
author: Alexandre Bigot <alexandre.bigot@cern.ch>, Strasbourg University
"""

import argparse

import yaml
from ROOT import (TF1, TCanvas, TDatabasePDG, TFile, TLatex, TLegend, TMath,
                  TMCParticle, gROOT, kAzure, kBlack, kBlue, kGreen, kFullCircle, kRed)

from StyleFormatter import SetGlobalStyle, SetObjectStyle

# enumerator
D0, DPLUS, LAMBDAC_TO_PKPI, LAMBDAC_TO_PK0S = 0, 1, 2, 3

# colours
RED = kRed + 1
BLUE = kBlue + 1
AZURE = kAzure + 4
Green = kGreen + 1

# conversion
GEV2MEV = 1000

# canvas dimensions
WIDTH = 520
HEIGHT = 500

# text size
SIZE_TEXT_LAT_ALICE = 26
SIZE_TEXT_LAT_LABEL_FOR_COLL_SYSTEM = 22
SIZE_TEXT_LAT_LABEL = 18
SIZE_TEXT_LEGEND = 18


def get_name_infile(particle):
    """
    Helper method to get the name of the input file according to the particle

    Parameters
    ----------
    - particle (int): particle ID

    Returns
    ----------
    - name_infile (string): name of the input file
    """

    name_infile = ""
    if particle == D0:
        name_infile = "../Results/D0/RawYieldsData_D0_pPb5TeV_FD_pos10.root"
    elif particle == DPLUS:
        name_infile = "../Results/Dplus/rawYield_Dplus_nonprompt_enhanced.root"
    elif particle == LAMBDAC_TO_PKPI:
        name_infile = ""
    elif particle == LAMBDAC_TO_PK0S:
        name_infile = ""

    return name_infile


def get_title_xaxis(particle):
    """
    Helper method to get the title of x axis according to the particle

    Parameters
    ----------
    - particle (int): particle ID

    Returns
    ----------
    - title_xaxis (string): title of x axis
    """

    title_xaxis = ""
    if particle == D0:
        title_xaxis = "#it{M}(K#pi) (GeV/#it{c}^{2})"
    elif particle == DPLUS:
        title_xaxis = "#it{M}(#piK#pi) (GeV/#it{c}^{2})"
    elif particle == LAMBDAC_TO_PKPI:
        title_xaxis = "#it{M}(pK#pi) (GeV/#it{c}^{2})"
    elif particle == LAMBDAC_TO_PK0S:
        title_xaxis = "#it{M}(pK^{0}_{S}) (GeV/#it{c}^{2})"

    return title_xaxis


def get_h_value_err(h, i_bin, convert_to_mev=False):
    """
    Helper method to get bin content and error of an histogram

    Parameters
    ----------
    - h (TH1): histogram
    - i_bin (int): bin number
    - convert_to_mev (int): apply conversion from GeV to MeV

    Returns
    ----------
    - value (float): bin content of h
    - error (float): bin error of h
    """

    value = h.GetBinContent(i_bin)
    error = h.GetBinError(i_bin)

    if convert_to_mev:
        value *= GEV2MEV
        error *= GEV2MEV

    return value, error


def draw_info(lat_label, particle):
    """
    Helper method to draw particle-dependent information on canvas

    Parameters
    ----------
    - lat_label (TLatex): TLatex instance
    - particle (int): particle ID
    """

    info = ""
    fnonprompt = ""
    if particle == D0:
        info = "D^{0} #rightarrow K^{#font[122]{-}}#pi^{+} and charge conj."
        # fnonprompt = "#it{f}_{ non-prompt}^{ raw} = 0.750 #pm  0.016 (stat.) #pm 0.008 (syst.)"
        fnonprompt = "#it{f}_{ non-prompt}^{ raw} = 0.531"
    elif particle == DPLUS:
        info = "D^{+} #rightarrow #pi^{+}K^{#font[122]{-}}#pi^{+} and charge conj."
        fnonprompt = "#it{f}_{ non-prompt}^{ raw} = 0.787 #pm  0.022 (stat.) #pm 0.016 (syst.)"
    elif particle == LAMBDAC_TO_PKPI:
        info = "#Lambda_{c}^{+}  #rightarrow pK^{#font[122]{-}}#pi^{+} and charge conj."
        fnonprompt = "#it{f}_{ non-prompt}^{ raw} = 0.630 #pm  0.056 (stat.) #pm 0.050 (syst.)"
    elif particle == LAMBDAC_TO_PK0S:
        info = "#Lambda_{c}^{+}  #rightarrow pK^{0}_{S} and charge conj."
        fnonprompt = "#it{f}_{ non-prompt}^{ raw} = 0.549 #pm  0.138 (stat.) #pm 0.055 (syst.)"

    lat_label.DrawLatex(0.19, 0.72, info)
    lat_label.DrawLatex(0.19, 0.16, fnonprompt)


def save_canvas(canvas, particle, pt_mins, pt_maxs, i_pt):
    """
    Helper method to save canvas according to particle

    Parameters
    ----------
    - canvas (TCanvas): a canvas
    - particle (int): particle ID
    """

    out_dir = "./"
    name = ""
    if particle == D0:
        name = "Dzero"
    elif particle == DPLUS:
        name = "Dplus"
    elif particle == LAMBDAC_TO_PKPI:
        name = "LambdacToPKPi"
    elif particle == LAMBDAC_TO_PK0S:
        name = "LambdacToPKzeroShort"

    for ext in ["pdf", "png", "eps"]:
        canvas.SaveAs(out_dir + "InvMassFit" + name + f"Pt_{pt_mins[i_pt]:.0f}_{pt_maxs[i_pt]:.0f}." + ext)


# pylint: disable=too-many-locals,too-many-statements
def main(particle, i_pt, cfg, batch):
    """
    Main method for a single bin (for article plots)

    Parameters
    ----------
    - particle (int): particle ID
    - i_pt (int): pT bin number
    """

    SetGlobalStyle(padtopmargin=0.07, padleftmargin=0.14, padbottommargin=0.125, titleoffsety=1.3, titleoffsetx=1., maxdigits=3)

    # import confiruables
    pt_mins = cfg["pPb5TeVFD"]["PtMin"]
    pt_maxs = cfg["pPb5TeVFD"]["PtMax"]
    mass_mins = cfg["pPb5TeVFD"]["MassMin"]
    mass_maxs = cfg["pPb5TeVFD"]["MassMax"]
    rebin = cfg["pPb5TeVFD"]["Rebin"]

    name_infile = get_name_infile(particle)

    file = TFile.Open(name_infile)

    hmean = file.Get("hRawYieldsMean")
    hsigma = file.Get("hRawYieldsSigma")
    if particle == D0:
        hsignal = file.Get("hRawYields")
    else:
        hsignal = file.Get("hRawYieldsSignal")

    name_hmass = f"hMass_{10*pt_mins[i_pt]:.0f}_{10*pt_maxs[i_pt]:.0f}"
    hmass = file.Get(name_hmass)
    hmass.Rebin(rebin[i_pt])

    title_xaxis = get_title_xaxis(particle)
    width_bin = hmass.GetBinWidth(i_pt+1)
    bin_max = hmass.GetMaximumBin()
    if particle == D0:
        ymin, ymax = 0., 1.15*(hmass.GetMaximum() + hmass.GetBinError(bin_max))
    else:
        ymin, ymax = 40, 1.3*(hmass.GetMaximum() + hmass.GetBinError(bin_max))
    title = f"{pt_mins[i_pt]:.0f} < #it{{p}}_{{T}} < {pt_maxs[i_pt]:.0f} GeV/#it{{c}};{title_xaxis};" \
        f"Counts per {width_bin*GEV2MEV:.0f} MeV/#it{{c}}^{{2}}"

    if particle == D0:
        fit_tot = file.Get(f"fTot_{pt_mins[i_pt]:.0f}.0_{pt_maxs[i_pt]:.0f}.0")
#    print(f"fTot_{pt_mins[i_pt]:.0f}_{pt_maxs[i_pt]:.0f}")
        fit_bkg = file.Get(f"fBkg_{pt_mins[i_pt]:.0f}.0_{pt_maxs[i_pt]:.0f}.0")
        fit_refl = file.Get(f"freflect")
    else:
        fit_tot = file.Get(f"fTot_{pt_mins[i_pt]:.0f}_{pt_maxs[i_pt]:.0f}")
        fit_bkg = file.Get(f"fBkg_{pt_mins[i_pt]:.0f}_{pt_maxs[i_pt]:.0f}")

    mean, err_mean = get_h_value_err(hmean, i_pt+1, True)
    sigma, _ = get_h_value_err(hsigma, i_pt+1, True)
    signal, err_signal = get_h_value_err(hsignal, i_pt+1)

    lat_alice = TLatex()
    lat_alice.SetNDC()
    lat_alice.SetTextSize(SIZE_TEXT_LAT_ALICE)
    lat_alice.SetTextFont(43)
    lat_alice.SetTextColor(kBlack)

    lat_label = TLatex()
    lat_label.SetNDC()
    lat_label.SetTextFont(43)
    lat_label.SetTextColor(kBlack)

    # lat_label = TLatex()
    # lat_label.SetNDC()
    # lat_label.SetTextFont(43)
    # lat_label.SetTextColor(kBlack)

    # str_mu = f"#it{{#mu}} = ({mean:.0f} #pm {err_mean:.0f}) MeV/#it{{c}}^{{2}}"
    # str_sigma = f"#it{{#sigma}} = {sigma:.0f} MeV/#it{{c}}^{{2}}"
    str_sig = f'#it{{S}} = {signal:.0f} #pm {err_signal:.0f}'

    if particle == D0:
        legend = TLegend(0.6, 0.51, 0.87, 0.72)
    else:
        legend = TLegend(0.62, 0.58, 0.85, 0.72)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextFont(43)
    legend.SetTextSize(SIZE_TEXT_LEGEND)
    legend.AddEntry(fit_tot, 'Total fit function', 'l')
    legend.AddEntry(fit_bkg, '#splitline{Combinatorial}{background}', 'l')
    if particle == D0:
        legend.AddEntry(fit_refl, 'K#minus#pi reflected', 'l')

    c = TCanvas("c", "", WIDTH, HEIGHT)
    frame = c.DrawFrame(mass_mins[i_pt], ymin, mass_maxs[i_pt], ymax, title)
    frame.GetYaxis().SetDecimals()

    SetObjectStyle(hmass, linewidth=3, linecolor=kBlack)
    SetObjectStyle(fit_tot, linewidth=3, linecolor=kBlue)
    SetObjectStyle(fit_bkg, linewidth=3, linecolor=kRed, linestyle=2)
    if particle == D0:
        SetObjectStyle(fit_refl, linewidth=3, linecolor=kGreen+2, linestyle=9)

    hmass.Draw("same")
    fit_bkg.Draw("same")
    fit_tot.Draw("same")
    if particle == D0:
        fit_refl.Draw("same")

    lat_alice.DrawLatex(0.19, 0.85, 'ALICE Preliminary')
    lat_label.SetTextSize(SIZE_TEXT_LAT_LABEL_FOR_COLL_SYSTEM)
    lat_label.DrawLatex(0.19, 0.79, 'pp, #sqrt{#it{s}} = 13.6 TeV')
    lat_label.SetTextSize(SIZE_TEXT_LAT_LABEL)
    draw_info(lat_label, particle)
    lat_label.DrawLatex(0.19, 0.66, f'{pt_mins[i_pt]:.0f} < #it{{p}}_{{T}} < {pt_maxs[i_pt]:.0f} GeV/#it{{c}}')
    lat_label.DrawLatex(0.7, 0.85, '#font[122]{-}0.5 < #it{y} < 0.5')
    # lat_label.DrawLatex(0.19, 0.64, str_mu)
    # lat_label.DrawLatex(0.19, 0.58, str_sigma)
    lat_label.DrawLatex(0.19, 0.6, str_sig)

    legend.Draw()

    c.Update()

    save_canvas(c, particle, pt_mins, pt_maxs, i_pt)

    if not batch:
        input("Press enter to exit")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Arguments")
    parser.add_argument("config", metavar="text", default="config.yml", help="config file name for ml")
    parser.add_argument("--batch", help="suppress video output", action="store_true")
    args = parser.parse_args()

    print("Loading analysis configuration: ...", end="\r")
    with open(args.config, "r", encoding="utf-8") as yml_cfg:
        configuration = yaml.load(yml_cfg, yaml.FullLoader)
    print("Loading analysis configuration: Done!")

    main(particle=D0, i_pt=0, cfg=configuration, batch=args.batch)
    # main(particle=DPLUS, i_pt=3, cfg=configuration, batch=args.batch)
