"""
Tools to get the optimal threshold and Bhattacharyya distance for a separability between two classes.
The field name of the classes have to be named CLASS.

S.Laventure / S.Alleaume (IRSTEA)    
"""

##Shapefile_en_entree=vector
##Classe_1=string 1
##Classe_2=string 2
##Nom_du_champs_des_statistiques=string minmean
##Nombre_echantillons=number 30
##th=output string
##Fichier_txt_de_sortie=output file

import sys, os
import random
import numpy as np # sort, ...
from collections import * # defaultdict
try :
    import ogr, gdal
except :
    from osgeo import ogr, gdal

def extract_data(Shapefile_en_entree, Classe_1, Classe_2, Nom_du_champs_des_statistiques):
    
    """
    To extract data from a shapefile
    
    """
    # Initial variable
    data_1 = []
    data_2 = []
        
    # Open the input vector
    data_source = ogr.GetDriverByName('ESRI Shapefile').Open(Shapefile_en_entree, 0)

    if data_source is None:
        print('Could not open file')
        sys.exit(1)

    print('Shapefile opening : ' + data_source.GetLayer().GetName())

    shp_ogr = data_source.GetLayer()
    # Projection
    # Import input shapefile projection
    srsObj = shp_ogr.GetSpatialRef()
    # Conversion to syntax ESRI
    srsObj.MorphToESRI()

    in_feature = shp_ogr.SetNextByIndex(0) # Polygons initialisation
    in_feature = shp_ogr.GetNextFeature()
    # Loop on input polygons to create a output polygons
    while in_feature:
        if str(in_feature.GetField('class')) in Classe_1:
            data_1.append(in_feature.GetField(str(Nom_du_champs_des_statistiques)))
        elif str(in_feature.GetField('class')) in Classe_2:
            data_2.append(in_feature.GetField(str(Nom_du_champs_des_statistiques)))
            
        in_feature.Destroy()
        # Next polygon
        in_feature = shp_ogr.GetNextFeature()

    return data_1, data_2

def seath(tab_1, tab_2, Classe_1, Classe_2, Nom_du_champs_des_statistiques, Nombre_echantillons):
    
    """
    To extract the optimal threshold for a separability between two classes
    
    Source article : SEaTH A new tool for automated feature extraction in the context of object-based image analysis S. Nussbaum et al.
    
    Source Info : Kenji Ose (IRSTEA) et Nathalie St-Geours (IRSTEA)
    """
    
    # Nb_samples random samples
    tab_1 = random.sample(tab_1, Nombre_echantillons)
    tab_2 = random.sample(tab_2, Nombre_echantillons)
    
    field_class = [tab_1, tab_2]

    # ## Compute Bhattacharyya distance ###
    # ##############################

    # Compute mean and variance
    m = defaultdict(list) # Average
    v = defaultdict(list) # Variance
    p = defaultdict(list) # Likelihood
    C = defaultdict(list) # Classes

    B = [] # Bhattacharyya distance
    # Optimal threshold
    seuil1 = []
    seuil2 = []
    seuil = []

    J = [] # Jeffries Matusita distance (measure separability between 2 classes on a 0 to 2 scale)
    threshold = '' # Optimal threshold
    th = '' #ÃÂ Optimal threshold for a output data

    for mark in range(len(field_class)):
        C[mark] = np.array(field_class[mark])
        m[mark].append(np.mean(C[mark]))
        v[mark].append(np.var(C[mark]))
        p[mark].append(1 / float(len(C[mark])))

    print "m : ", m
    print "v : ", v

    # Mean, standard deviation and likelihood initialisation phase for 2 classes 
    m1 = m[0]
    m2 = m[1]
    v1 = v[0]
    v2 = v[1]
    p1 = p[0]
    p2 = p[1]

    for i in range(len(m[0])):
        B.append(( (1/float(8)) * ( (m1[i] - m2[i])**2 ) * (2 / ( v1[i] + v2[i] )) ) + ( (1/float(2)) * np.log( ( v1[i] + v2[i] ) / ( 2 * np.sqrt(v1[i] * v2[i]) ) ) ))
        J.append(2 * ( 1 - np.exp( -B[i] ) ))
        
        ### Optimal threshold calculation ###
        ######################
        # Bayes theorem solution
        A = np.log( np.sqrt( v1[i] / v2[i] ) * ( p2[i] / p1[i] ))
        D = np.sqrt( v1[i] * v2[i] ) * np.sqrt( ( m1[i] - m2[i] )**2 + 2 * A * ( v1[i] -  v2[i] ) )
        seuil1.append(( 1 / ( v1[i] - v2[i] ) ) * ( m2[i] * v1[i] - m1[i] * v2[i] + D ))
        seuil2.append(( 1 / ( v1[i] - v2[i] ) ) * ( m2[i] * v1[i] - m1[i] * v2[i] - D ))
        
        # Optimal threshold
        # Logical condition depending on article figure 2
        if ( seuil1[i] > m2[i] and seuil1[i] < m1[i] ) or ( seuil1[i] > m1[i] and seuil1[i] < m2[i] ) :
            print "Valid  threshold !"
        else:
            seuil1[i] = ""
        
        if ( seuil2[i] > m2[i] and seuil2[i] < m1[i] ) or ( seuil2[i] > m1[i] and seuil2[i] < m2[i] ) :
            print "Valid  threshold !"
        else:
            seuil2[i] = ""

        # Final condition
        if ( seuil1[i] == "" and seuil2[i] == "" ) or ( seuil1[i] != "" and seuil2[i] != "" ):
            seuil.append("")
        elif ( seuil1[i] != "" and seuil2[i] == "" ):
            seuil.append(seuil1[i])
        elif ( seuil1[i] == "" and seuil2[i] != "" ):
            seuil.append(seuil2[i])

    print("Bhattacharyya distance ", B)
    print("J : ", J)
    print("Threshold 1 : ", seuil1)
    print("Threshold 2 : ", seuil2)
    print("Optimal threshold :", seuil)

    for i in range(len(seuil)):
        if seuil[i] != "" and m1[i] > m2[i]:
            print('To extract the class 1, we need to this index > ' + str(seuil[i]))
            threshold = 'To extract the class 1, we need to this index ' + str(Nom_du_champs_des_statistiques) + ' > ' + str(seuil[i])
            th = 'CASE WHEN "' + str(Nom_du_champs_des_statistiques) + '" ' + '>' + str(seuil[i]) + ' THEN ' + str(Classe_1) + ' ELSE ' + str(Classe_2) + ' END'
        elif seuil[i] != "" and m1[i] < m2[i]:
            print('To extract the class 1, we need to this index < ' + str(seuil[i]))
            threshold = 'To extract the class 1, we need to this index ' + str(Nom_du_champs_des_statistiques) + ' < ' + str(seuil[i])
            th = 'CASE WHEN "' + str(Nom_du_champs_des_statistiques) + '" ' + '<' + str(seuil[i]) + ' THEN ' + str(Classe_1) + ' ELSE ' + str(Classe_2) + ' END'
        else:
            print('Not discrimination for this index!')
            threshold = 'Not discrimination for this index!'
            th=''
    #        sys.exit(1)
    return threshold, th

# Main process
# Convert input classe name to list. To permit to use several variables.
if ',' in Classe_1:
    Classe_1 = Classe_1.split(',')
else:
    Classe_1 =[Classe_1]
if ',' in Classe_2:
    Classe_2 = Classe_2.split(',')
else:
    Classe_2 =[Classe_2]
# Call 2 main functions
tab_1, tab_2 = extract_data(Shapefile_en_entree, Classe_1, Classe_2, Nom_du_champs_des_statistiques)
# Tricks to launch again Seath proccessing if there isn't discrimination. 11 times max.
p = 0
th = ''
while p < 11 :
    if th == '' :
        threshold, th = seath(tab_1, tab_2, Classe_1[0], Classe_2[0], Nom_du_champs_des_statistiques, Nombre_echantillons)
        p = p + 1
    else:
        p = 11

# Extract results in a file
if os.path.exists(Fichier_txt_de_sortie):
    os.remove(Fichier_txt_de_sortie)
f_out = open(Fichier_txt_de_sortie, "wb")
f_out.write(str(threshold))
f_out.close()
