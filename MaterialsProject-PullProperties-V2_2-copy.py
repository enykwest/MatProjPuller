#!/usr/bin/env python3
"""
Created on Fri Sep 11 15:37:14 2020
Code for pulling data from Materials Project.
When handed a chemical space and a set of properties, this code will search
the materials database for matching materials, retreive the properties,
and store them in a csv file.

@author: Erik Nykwest
"""
#%% Define Variables (Inputs)

#Ease of Use
Ele1 = "Re"
Ele2 = "O"

# Where should everything be saved? (required)
path = r'FullPathToDir\{}x{}y'.format(Ele1, Ele2) # to catch odd behavior, path must already exist

# Database name (required)
data = r"{}x{}y-V2_2.csv".format(Ele1, Ele2) # we will start by using CSV files, but this is infeasible once the database gets large.

# Materials Project Key (required)
MPKey = r'PutYourMPKeyHere'

# Chemical Space (required)
# Follow MongoDB format
# https://pymatgen.org/usage.html#the-query-method
# https://github.com/materialsproject/mapidoc

ChemSpc = {
            "elements":{"$in":[Ele1], "$all": [Ele2]},
            "nelements":2
            } # Dictionary, Only consider binary oxides that include [ Fe,
            # Ti, Ni,Cr, Os, Pt, Si, V, W, Hf, Rh, Ir ]

### Properties of interest
props = dict() # Setting up an empty dictionary which we will use later

# List of properties we can pull using MPRester().query
# Defaults
query_props = ['task_id', 'structure' ] # Minimum requirements
                         
query_props = ['task_id', 'full_formula', 'pretty_formula',
               'structure', 'e_above_hull', 'formation_energy_per_atom'] # 'spacegroup'?


# How to retreive Cohesive Energy
# Why can't I pull this with a query?!?!?
def CE(S, ID):
    with MPRester(MPKey) as MPR:
        ce = BePatient(lambda: MPR.get_cohesive_energy(ID, per_atom=True))
        return ce  

# A Dict() of properties pulled via a method on the structure object
method_props = {
                "cohesive_energy": CE,
                }

# Atomistic properties, properties that are atomic site dependent
atomic_method_props = {
                "cohesive_energy": CE,
                }

# Wrap it all back up into our dictionary
props['query_props'] = query_props
props['method_props'] = method_props

#%% Import Modules
import os
import time

#Pymatgen materials package
import pymatgen as pmg
import pymatgen.io.cif
#from pymatgen import MPRester # depreciated
from pymatgen.ext.matproj import MPRester # also depreciated, but still works as of 10/2022

#Data science package to handle data
import pandas as pd
import numpy as np
import re

# Modules for LocalGeometry Finder
from pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder import LocalGeometryFinder
from pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies import SimplestChemenvStrategy
from pymatgen.analysis.chemenv.coordination_environments.structure_environments import LightStructureEnvironments
# Possible work around for bad pymatgen oxidation state guesses
from pymatgen.analysis.bond_valence import BVAnalyzer


#%% Define Functions
def BePatient(func):
    # MPREST requests can be buggy and sometimes time out
    # This decorator allows for multiple tries to connect before giving up
    x = None
    for i in range(10): # try no more than 10 times
        try:
            x = func()
            break
        except: 
            print("Time out, trying again")
            time.sleep(1) # wait 1 second and try again
            continue
        
    return x


# Because the authors are EVIL, the .get_shannon_radius method
# only takes in ROMAN NUMERALS for coordination number...
def roman(CN):
    
    RN = "None" # Number as roman Numeral
    
    numerals = { 1:"I", 2:"II", 3:"III", 4:"IV", 5:"V",
                 6:"VI", 7:"VII", 8:"VIII", 9:"IX", 10:"X",
                 -1:"None"}
    
    try:
        RN = numerals[CN]
    except:
        pass # Leave RN as "None"
        
    return RN


def UpdateDatabase(MPKey, data, ChemSpc, props, replace = False, write_cif=False):
    """
    Inputs:
    MPKey = Materials Project Database API Key
    data = .csv file to store properties
    ChemSpc = a string used to query the materials project database for a list of materials
    props is a dictionary of material properties to retrieve about each material.
        - props['query_props'] is a list of of properties to be requested from the MP
        database along with the initial query.
        - props['method_props'] is a dictionary of functions that act that act on the 
        retrieved structure to generate an output. Typically a method function.
        
    Variables
    FullData = The Pandas DataFrame that holds ALL of the data for the
        materials and gets written to a .csv file at each iteration.
    pds = a Pandas DataSeries that holds all of the INTRINSIC material properties
        (like cohesive energy).
    pdf = a Pandas DataFrame that holds all of the atom-level properties
        (like coordination enviornment).
    """
    
    # Try to add data to existing database
    try:
        # If there is an existing csv file, open it, append new data, and save it
        FullData = pd.read_csv(data, index_col=0)
    except(FileNotFoundError):
        # If these is no csvfile, make one
        FullData = pd.DataFrame(columns=['task_id'],dtype=str)
        
    #Connect to Materials Project Database and pull material properties
    with MPRester(MPKey) as MPR:
        # list of dictionaries of material properties
        mats = BePatient(lambda: MPR.query(ChemSpc, props['query_props']))
    
    # For each material, extract propeties from materials dictionary
    for mdict in mats:
        
        # Do we skip materials already in database?
        if replace:
            # Replace Existing Data
            FullData = FullData[ FullData['task_id'] != mdict['task_id'] ]

        elif FullData["task_id"].isin([mdict['task_id']]).any():
            continue # skip materials already in database
            
        
        # Extract structure
        MatStruct = mdict.pop('structure') # remove and store structure
        
        # Add Oxidation state
        try:
            # BVAnalyzer gives better answers, but sometimes can't be solved
            BV=BVAnalyzer()
            MatStruct=BV.get_oxi_state_decorated_structure(structure=MatStruct)
        except:
            # Returns 0 when guess fails
            # Returns decimal oxidation states sometimes
            MatStruct.add_oxidation_state_by_guess() # WARNING: This changes the
                # species in both the structure and the cif file to add the oxidation
                # state. For example: Fe -> Fe2.05+
                
        # Add query properties to a pandas data series
        pds = pd.Series(mdict)
        
        # Run and store method properties
        for mkey in props['method_props']:
            pds[mkey] = props['method_props'][mkey](MatStruct,pds['task_id'])
            
        
        
        ### Extract atomic site properties and assemble dataframe
        pdf = pd.DataFrame(dtype=str)
        MatStructDF = MatStruct.as_dataframe()
        
        '''
        The "Species" is currently being stored as a 
        'pymatgen.core.composition.Composition', but when you 
        write to a CSV file it becomes a string, which will read:      
            
        Species:OxidationState:Count
        
        For example: Fe2.05+1 means there is 1 Fe atom with oxidation state +2.05
        
        To avoid these issues the next couple steps extract the element 
        and oxidation state from the composition.
        '''
        MatStructDF["Species"] = MatStructDF.Species.apply(lambda x: str(x)[:-1])
        
        # Seperate element from oxidation state
        EleOxiCount = MatStructDF['Species'].apply(
            lambda x: re.match(r'([A-Z,a-z]*)([0-9]?\.?[0-9]*[+-])', x))
        Ele = EleOxiCount.apply(lambda x: x.group(1))
        Oxi = EleOxiCount.apply(lambda x: x.group(2))
        # BVAnalyzer is stupid and writes O1- as O- so to fix that..
        Oxi[ Oxi.isin(['+']) ] = '1+'
        Oxi[ Oxi.isin(['-']) ] = '1-'
        
        #Switch Order
        Oxi = Oxi.apply( lambda x: x[-1] + x[:-1] )
        
        # Update DateFrame
        pdf['Species'] = MatStructDF['Species']
        pdf['Element'] = Ele
        pdf['Oxidation_State'] = Oxi
        pdf['magmom'] = MatStructDF['magmom']
        
        
        
        ### Calculate Structure Enviroments (Coordination Enviornment)
        # lse.coordination_environments is a list of dictionaries of properties
        lse = LGF(MatStruct)
        
        #filter coordination environments to replace None values with indexable values
        filtered=[
                  [{'ce_symbol': 'None:-1', 'ce_fraction': 0, 'csm': 1, 'permutation': [0, 0, 0, 0, 0, 0, 0]}]
                  if v is None else v for v in lse.coordination_environments
                  ]
        
        # The above object doesn't play nice with list slicing for some weird reason...
        # Add ce_symbol to pdf
        filtered = np.array(filtered).T
        ce = list()
        for i in filtered[0]:
            ce.append(i['ce_symbol'])     
        pdf['ce_symbol'] = ce
        
        # Add CSM to pdf
        csm = list()
        for i in filtered[0]:
            csm.append(i['csm'])     
        pdf['csm'] = csm
        
        # Calculate Average Coordination
        c = [ int(x.rsplit(":")[1]) for x in ce ]
        c = [ x for x in c if x >= 0 ] # remove "None" values from our average
        # AVERAGE is a gobal property so we add it to pds NOT pdf
        pds["Average Coordination"] = np.average(c)
        pds["STD of Coordination"] = np.std(c)
        
        
        
        ### Get Shannon Radius
        # We need a species object to invoke .get_shannon_radius method
        CreateSpecies = pymatgen.core.periodic_table.Species
        
        Species = [ [CreateSpecies(S[0], round(float(S[1]))), round(float(S[2].rsplit(":")[1])) ] for S in 
                   zip(pdf.Element, pdf.Oxidation_State, pdf.ce_symbol)]
        # Convert Coordination Number to Roman Numerals
        Species = [ [S, roman(CN) ] for (S, CN) in Species ]
        
        ShannonIR = [] # Shannon Ionic Radius
        ShannonCR = [] # Shannon Ionic Crystal Radius
        IonicR = [] # Ionic Radius from MatsPrj
        for S, CN in Species:
            IonicR.append(S.ionic_radius)
            
            try:
                ShannonIR.append(S.get_shannon_radius(CN, radius_type="ionic"))
                ShannonCR.append(S.get_shannon_radius(CN, radius_type="crystal"))
            except:
                ShannonIR.append(np.nan)
                ShannonCR.append(np.nan)
                
        # Add Ionic Radii to pdf
        pdf["ionic_radius"] = IonicR
        pdf["Shannon_IR"] = ShannonIR
        pdf["Shannon_CR"] = ShannonCR
    
    
    
        ### If wanted, write the structure to a cif file
        if write_cif:
            cif_file = pds["full_formula"].replace(r" ", r"")+"_"+pds["task_id"]+'.cif' # ChemicalFormula_MaterialID.cif
            pmg.io.cif.CifWriter(MatStruct).write_file(cif_file)
            pds["Cif File"] = cif_file #Location of the cif file
        else:
            pds["Cif File"] = "None"
        
        # Add "Global" properties to atom-decomposed dataframe
        # (e.g. every atom in the same material has the same "Cohesive Energy")
        for i in pds.index:
            pdf[i] = pds[i]
            
        # Add data to existing database, May create duplicate entries
        FullData = FullData.append(pdf, ignore_index=True)

        #debug
        #break
        
        # Save after each pass? or only at end? add or remove this block from the loop!
        # Save File
        FullData.to_csv(data) #WARNING: The "Species" is currently being stored
        # as a 'pymatgen.core.composition.Composition', but when you write to
        # a CSV file it becomes a string which will read
        # Species:OxidationState:Count
        # For example: Fe2.05+1 means there is 1 Fe atom with oxidation state +2.05
        print("Saving")
   
    return

            
def LGF(structure):
    # Local Geometry Finder
    
    ### Setting up LocalGeometryFinder
    lgf = LocalGeometryFinder() 

    lgf.setup_parameters(centering_type='centroid', include_central_site_in_centroid=True,
                         structure_refinement=lgf.STRUCTURE_REFINEMENT_NONE,
                         bva_distance_scale_factor=1.015,
                         spg_analyzer_options={'symprec':0.1, 'angle_tolerance': 5}
                        )
    
    #logging.basicConfig(filename='chemenv_structure_environments' + str(i) + '.log', format='%(levelname)s:%(module)s:%(funcName)s:%(message)s', level=logging.DEBUG)
    lgf.setup_structure(structure)
    
    #parameters suggested by David Waroquiers in email
    strategy = SimplestChemenvStrategy()
    
    
    ### Compute Enviorments and filter
    try: #Some structures throw errors, but I'm not about to debug LSE
        se = lgf.compute_structure_environments(maximum_distance_factor=1.4,
                                            minimum_angle_factor=0.3,
                                            ) # removed only_atoms="U"
        
        # From all structure enviorments, use chosen strategy to downselect to enviorments of interest
        lse = LightStructureEnvironments.from_structure_environments(strategy=strategy, structure_environments=se)
        
    except:
        #Some structures throw errors, but I'm not about to debug LSE
        def lse():
            # Dummy function to attach properties to.
            pass
        
        lse.coordination_environments = [ None ]*len(structure.as_dataframe())
        pass
    
    #debug
    #coord = lse.coordination_environments
    #print(structure[0])
    #print(coord)
    
    
    
    return lse           
#%% Move to folder, Pull IDs from Materials Project

os.chdir(path) # Move to working Directory
UpdateDatabase(MPKey, data, ChemSpc, props, replace=False, write_cif=True) # Pull data, overwrite CSV file


print("Done")