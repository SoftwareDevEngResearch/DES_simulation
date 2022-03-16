#!/usr/bin/env python
# coding: utf-8

import argparse
import sys
from bs4 import BeautifulSoup
import requests
import numpy as np
import pandas as pd
import math

#############################################################################################
#                           COMMAND LINE INTERFACE SETUP
#############################################################################################

parser = argparse.ArgumentParser(description='Determine peaks in DES spectrum given isotope, activity, and absorber dimensions')
parser.add_argument('-i', '--isotope', required=True, type=str, help='isotope in DES spectrum, example: Am-241')
parser.add_argument('-a', '--activity', required=True, type=float, help='activity of sample in Bq')
parser.add_argument('-x', '--width', required=True, type=float, help='absorber width in cm')
parser.add_argument('-y', '--length', required=True, type=float, help='absorber length in cm')
parser.add_argument('-z', '--thickness', required=True, type=float, help='absorber thickness in cm')

args = parser.parse_args(sys.argv[1:])

isotope = args.isotope
activity = args.activity
x = args.width
y = args.length
z = args.thickness

t_vals = [x/2, y/2, z/2]    #places the sample in the center of absorber

#############################################################################################
#                           FUNCTION DEFINITIONS
#############################################################################################

def get_url(ISOTOPE):
    return f'https://www.nndc.bnl.gov/nudat2/decaysearchdirect.jsp?nuc={ISOTOPE}&unc=nds'

def isfloat(value):
    try:
        float(value)
        return True
    except ValueError:
        return False

def get_gammas(page):
    """Grabs Q-value and creates arrays of gamma ray energies and intensities from tables found at nndc.bnl.gov"""
    soup = BeautifulSoup(page.content,'html.parser')
    table = soup.find_all('table')[-2] #Gammas table
    tr = table.find_all('tr')
    energies = []
    intensities = []
    for row in tr:
        energy = row.find_all('td')[1].text.strip().split(' ')[0]
        if isfloat(energy):
            energies.append(float(energy))  #from gammas table, create array of energies
        intensity = row.find_all('td')[2].text.strip().split(' ')[0]
        if isfloat(intensity):
            intensities.append(float(intensity)/100)  #from gammas table, create array of intensities
    QTable = soup.find_all('table')[0]
    tr = QTable.find_all('tr')
    QValue = float(tr[1].find_all('td')[-4].text.replace(u'\xa0','').strip())
    escape_energies = QValue - np.asarray(energies)
    return np.asarray(energies), escape_energies, np.asarray(intensities), QValue

def get_des_energies(isotope):
    url = get_url(isotope)
    isotope_url = requests.get(url)
    gamma_energies, escape_energies, intensities, qval = get_gammas(isotope_url)
    return pd.DataFrame({'Gamma Energy':gamma_energies, 'Escape Energy':escape_energies, 'Intensity':intensities}), qval

def del_rows(df):
    df.drop(df.index[df['Intensity'] < 0.1], inplace = True) #drop rows with intensity < 1%
    gamma_energies = df['Gamma Energy']
    return df, np.asarray(gamma_energies)

def attenuation(activity, mu_vals, t_vals):
    """Performs gamma attenuation through matter calculation to determine which gammas escape"""
    escape_mu = []
    for x in mu_vals:
        for y in t_vals:
            I = activity * math.exp(-x * y * 19.3)
            if I > activity*0.001:      #adds gamma energy to list of escaped gammas if signifcant activity remains
                [escape_mu.append(float(x)) for x in mu_vals if x not in escape_mu]
    return np.asarray(escape_mu)

def recall_energies(escape_mu, df):
    escaped_gammas = []
    for x in escape_mu:
        energy = df.loc[df['Mu'] == x]['Gamma Energy']
        escaped_gammas.append(float(energy))
    return np.asarray(escaped_gammas)

def recall_escape_energies(escaped_gammas, df):
    escape_energies = []
    for x in escaped_gammas:
        energy = df.loc[df['Gamma Energy'] == x]['Escape Energy']
        escape_energies.append(float(energy))
    return np.asarray(escape_energies)

def print_DES_spectrum(isotope, escape_energies, qval):
    print(isotope + " will have escape peaks at: ", escape_energies)
    print(isotope + " has a full energy peak at: ", qval)

nuclear_data = pd.read_csv("Au_attencoeff.csv",index_col='mu/rho')

def get_mu(gamma_energies, df, colname):
    """From csv file, grabs appropriate attenuation coefficient based on gamma energy"""
    mu_vals = []
    for x in gamma_energies:
        exactmatch = df[df[colname] == x]
        if not exactmatch.empty:
            mu_vals.append(float(exactmatch.index))
        else:
            lowerneighbor_ind = df[df[colname] < x][colname].idxmax()
            upperneighbor_ind = df[df[colname] > x][colname].idxmin()
            mu = min([lowerneighbor_ind, upperneighbor_ind])
            mu_vals.append(float(mu))
    return np.asarray(mu_vals), pd.DataFrame({'Mu':mu_vals, 'Gamma Energy':gamma_energies})

#############################################################################################
#                           FUNCTION CALLS
#############################################################################################

full_df, qval = get_des_energies(isotope)
red_df, gamma_energies = del_rows(full_df)
mu_vals, muE_df = get_mu(gamma_energies, nuclear_data, 'Energy.1')
escape_mu = attenuation(activity, mu_vals, t_vals)
escaped_gammas = recall_energies(escape_mu, muE_df)
escape_energies = recall_escape_energies(escaped_gammas, red_df)
print_DES_spectrum(isotope, escape_energies, qval)
