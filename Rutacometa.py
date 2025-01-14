#!/usr/bin/env python
# coding: utf-8

# In[22]:


import numpy as np
from matplotlib import pyplot as plt
from matplotlib.collections import LineCollection
from skyfield.api import Star, load
from skyfield.data import hipparcos, mpc, stellarium
from skyfield.projections import build_stereographic_projection
from astroquery.simbad import Simbad
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.wcs import WCS
from astropy.visualization.wcsaxes import WCSAxes
from skyfield.constants import GM_SUN_Pitjeva_2005_km3_s2 as GM_SUN


# In[23]:


def mapa_cometa(inicio, final):
    #Tiempo en que el se grafica
    ts = load.timescale()
    t_comet = ts.utc(2025, 1, range(inicio, final+1,1)) #se le suma 1 para que quedé en la fecha final ingresada
    t = t_comet[len(t_comet) // 2]  # fecha media
    
    # Diseño de la grafica
    plt.style.use("dark_background")
    plt.rcParams['font.family']='DeJavu Serif'
    plt.rcParams['font.serif'] = ['Times New Roman']
    #Para que el diseño salga se centró en M72, para cometas distintos se debe ver el cielo y ver en que se centra
    OBJECT = 'M72'
    #campo de visión
    FOV = 30.0
    #magnitud máxima de las estrellas que se veran
    MAG = 5.0

    TABLE = Simbad.query_object(OBJECT)
    RA = TABLE['RA'][0]
    DEC = TABLE['DEC'][0]
    COORD = SkyCoord(f"{RA} {DEC}", unit=(u.hourangle, u.deg), frame='fk5')
    
    eph=load('de440s.bsp') #esto puedo fallar, ojo
    earth = eph['earth']
    sun = eph['sun']
    
    #descarga los datos de la posición del cometa
    
    with load.open(mpc.COMET_URL) as f:
        comets = mpc.load_comets_dataframe(f)

        comets = (comets.sort_values('reference')
              .groupby('designation', as_index=False).last()
              .set_index('designation', drop=False))

    row = comets.loc['C/2024 G3 (ATLAS)']
    comet = sun + mpc.comet_orbit(row, ts, GM_SUN)
    
    url=('https://raw.githubusercontent.com/Stellarium/stellarium/master'
       '/skycultures/modern_st/constellationship.fab') #esto también puede fallar
    with load.open(url) as f:
        constellations = stellarium.parse_constellations(f)

    edges = [edge for name, edges in constellations for edge in edges]
    edges_star1 = [star1 for star1, star2 in edges]
    edges_star2 = [star2 for star1, star2 in edges]

    # Se saca los datos de las estrellas del Hipparcos
    with load.open(hipparcos.URL) as f:
        stars = hipparcos.load_dataframe(f)
        
    center = earth.at(t).observe(comet)
    projection = build_stereographic_projection(center)
    
    
    star_positions = earth.at(t).observe(Star.from_dataframe(stars))
    stars['x'], stars['y'] = projection(star_positions)

    # Para que aparezcan las estrellas.
    bright_stars = (stars.magnitude <= MAG)
    magnitude = stars['magnitude'][bright_stars]
    marker_size = (0.5 + MAG - magnitude) ** 2.0
    comet_x, comet_y = projection(earth.at(t_comet).observe(comet))
    
    xy1 = stars[['x', 'y']].loc[edges_star1].values
    xy2 = stars[['x', 'y']].loc[edges_star2].values
    lines_xy = np.rollaxis(np.array([xy1, xy2]), 1)

    
    #Aquí se hace la gráfica
    fig = plt.figure(figsize=[10, 10])
    wcs = WCS(naxis=2)
    wcs.wcs.crpix = [1, 1]
    wcs.wcs.cdelt = np.array([-FOV / 360, FOV / 360])
    wcs.wcs.crval = [COORD.ra.deg, COORD.dec.deg]
    wcs.wcs.ctype = ["RA---STG", "DEC--STG"]

    ax = fig.add_subplot(111, projection=wcs)

    # Líneas de las constelaciones
    ax.add_collection(LineCollection(lines_xy, colors='cyan', linewidths=1, linestyle='-'))

    # Las estrellas
    ax.scatter(stars['x'][bright_stars], stars['y'][bright_stars],
               s=marker_size, color='white', zorder=2)

    comet_color = 'red' #color de la marca del cometa
    offset = 0.009 #que tan lejos está la fecha de la marca del cometa

    ax.plot(comet_x, comet_y, '+', c=comet_color, zorder=3)

    for xi, yi, tstr in zip(comet_x, comet_y, t_comet.utc_strftime('%d')): #si se ´pone m%/%d sale mes y el días como mes/dia
        tstr = tstr.lstrip('0')
        text = ax.text(xi + offset, yi - offset, tstr, color=comet_color,
                   ha='left', va='top', fontsize=9, weight='bold', zorder=-1)
        text.set_alpha(0.5)

    angle = np.pi - FOV / 360.0 * np.pi
    limit = np.sin(angle) / (1.0 - np.cos(angle))

    # limites de la gráfica, se pusó tres porque queda mejor
    ax.set_xlim(-3*limit, 3*limit)
    ax.set_ylim(-3*limit, 3*limit)
    ax.set_aspect('equal')

    # cuadricula de las coordenadas ecuatoriales
    ax.coords.grid(True, color='white', linestyle='dotted')

    
    ax.coords[0].set_axislabel('RA (horas)')
    ax.coords[1].set_axislabel('Dec (grados)')
    ax.coords[0].set_major_formatter('hh:mm:ss')
    ax.coords[1].set_major_formatter('dd:mm:ss')


    ax.set_title('Trayectoria del Cometa C/2024 G3 del {} al {}'.format(
        t_comet[0].utc_strftime('%Y %B %d'),
        t_comet[-1].utc_strftime('%Y %B %d'),))


# In[ ]:





# In[24]:


def hacer_mapa():
    print('Ingrese entre que días de enero quiere ver el cometa:')
    inicio = int(input('Fecha inicial: '))
    final = int(input('Fecha final: '))
    if inicio < 1: 
        print('Fecha fuera de rango')
    elif final > 31:
        print('Fecha fuera de rango')
    else: 
        mapa_cometa(inicio,final)


# In[52]:


hacer_mapa()


# In[53]:


# todo lo que esta abajo de aqui, ignorelo


# In[ ]:





# In[55]:


2024-1054


# In[31]:


def Hsp(lat,delta):
    Hsp=np.arccos((-9.89e-3)-np.sin(lat)*np.sin(delta)/(np.cos(lat)*np.cos(delta)))
    Hsp1=360-Hsp
    return Hsp, Hsp1


# In[33]:


def Asp(lat, delta):
    Asp=np.arccos(np.sin(delta)/np.cos(lat))
    Asp1=360-Asp
    return Asp, Asp1


# In[35]:


def TSLsp(alpha,Hs, Hp):
    TSLs=alpha+Hs
    TSLp=alpha+Hp
    return TSLs, TSLp


# In[37]:


def TSGsp_w(TSLs, TSLp, long):
    #long debe estar en grados
    TSGs_w=TSLs+(long/15)
    TSGp_w=TSLp+(long/15)
    return TSGs_w, TSGp_w


# In[38]:


def TSGsp_e(TSLs, TSLp, long):
    #long debe estar en grados
    TSGs_e=TSLs-(long/15)
    TSGp_e=TSLp-(long/15)
    return TSGs_e, TSGp_e


# In[46]:


def TSG0(FJ):
    T=(FJ-2451545.0)/36525
    TSG0=6.697374558333333+2400.0513369055557*T+(2.586111111111111e-05)*T**2
    if TSG0>24:
        TSG0=TSG0-int(TSG0/24)*24
    else:
        TSG0=TSG0
    return TSG0


# In[47]:


def TUsp(TSGs, TSGp, TSG0):
    TUs=(TSGs-TSG0)/1.0027379
    TUp=(TSGp-TSG0)/1.0027379
    return TUs, TUp


# In[48]:


def TLsp(TUs, TUp, HH):
    TLs=TUs+HH
    TLp=TUp+HH
    return TUs, TUp


# In[49]:


###########################################################################################################################
#PRUEBA
###########################################################################################################################


# In[51]:


alpha_p=(5+(55/60))*15 #que en grados



# In[ ]:





# In[ ]:





# In[50]:


###########################################################################################################################
#FINAL DE LA PRUEBA
###########################################################################################################################


# In[43]:


0+ (0/60) + (0.09310/3600)


# In[40]:


60*60


# In[45]:


int(45.99)


# In[ ]:




