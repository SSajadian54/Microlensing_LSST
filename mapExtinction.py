from __future__ import print_function
from astropy.coordinates import SkyCoord
from dustmaps.bayestar import BayestarQuery
from dustmaps.sfd import SFDQuery
import astropy.units as units
from dustmaps.bayestar import BayestarQuery
from dustmaps.planck import PlanckQuery
#from dustmaps.config import config
#config.reset()


bayestar = BayestarQuery(max_samples=2, version='bayestar2019')


#coords = SkyCoord(180., 0., unit='deg', frame='galactic')
#bayestar = BayestarQuery(max_samples=2, version='bayestar2019')
#E = bayestar(coords, mode='random_sample')


coords = SkyCoord(10.0*units.deg, 2.0*units.deg, distance=(5+1.0)*1000.*units.pc, frame='galactic')
sfd = SFDQuery()
ebv = sfd(coords)

print('E(B-V) = {:.3f} mag'.format(ebv))



for i in range(8): 
    coords = SkyCoord(10.0*units.deg, 2.0*units.deg, distance=(i+1.0)*1000.*units.pc, frame='galactic')
    E = bayestar(coords, mode='mean'   )
    print ("Reddening: ",  E,   3.1*E  ,  (i+1.0)*1000. )
