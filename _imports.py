# Imports here
import numpy as np
import matplotlib.pyplot as plt
import os

# Conversions and constants

radius_e_nmi = 3443.922786 # earth radius nmi
radius_e_km  = 6378.145

du_nmi = 3443.922786 # n.mi. per DU

nmi_ft = 2.092567257e7 / 3443.922786 # ft / nmi
ft_nmi = nmi_ft**-1 # nmi / ft

ft_mt = 1 / 0.3048
mt_ft = ft_mt**-1

nmi_km = 1/1.852
km_nmi = nmi_km**-1

sec_day = (60*60*24)**-1
sec_hour = (60*60)**-1