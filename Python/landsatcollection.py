#----------------------------------------------------------------------------------------#
#---------------------------------------//GEESEBAL//-------------------------------------#
#GEESEBAL - GOOGLE EARTH ENGINE APP FOR SURFACE ENERGY BALANCE ALGORITHM FOR LAND (SEBAL)
#CREATE BY: LEONARDO LAIPELT, RAFAEL KAYSER, ANDERSON RUHOFF AND AYAN FLEISCHMANN 
#PROJECT - ET BRASIL https://etbrasil.org/
#LAB - HIDROLOGIA DE GRANDE ESCALA [HGE] website: https://www.ufrgs.br/hge/author/hge/
#UNIVERSITY - UNIVERSIDADE FEDERAL DO RIO GRANDE DO SUL - UFRGS
#RIO GRANDE DO SUL, BRAZIL

#DOI
#VERSION 0.1
#CONTATC US: leonardo.laipelt@ufrgs.br

#----------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------#

#PYTHON PACKAGES
#Call EE
import ee

#SURFACE REFLECTANCE
#ATMOSPHERICALLY CORRECTED

#GET LANDSAT 8 COLLECTIONS BY PATH ROW
def fexp_landsat_8PathRow(start_date,end_date,n_path, n_row,th_cloud_cover):
    col_SR_L8 =(ee.ImageCollection('LANDSAT/LC08/C01/T1_SR')
                        .filterDate(start_date, end_date)
                        .filterMetadata('WRS_PATH', 'equals', n_path)
                        .filterMetadata('WRS_ROW', 'equals', n_row)
                        .select([0,1,2,3,4,5,6,7,10],["UB","B","GR","R","NIR","SWIR_1","SWIR_2","BRT","pixel_qa"])
                        .filterMetadata('CLOUD_COVER', 'less_than', th_cloud_cover));
    return col_SR_L8;

#GET LANDSAT 7 COLLECTIONS BY PATH ROW
def fexp_landsat_7PathRow(start_date,end_date,n_path, n_row,th_cloud_cover):
    
    col_SR_L7 =(ee.ImageCollection('LANDSAT/LE07/C01/T1_SR')
                        .filterDate(start_date, end_date)
                        .filterMetadata('WRS_PATH', 'equals', n_path)
                        .filterMetadata('WRS_ROW', 'equals', n_row)                        
                        .select([0,1,2,3,4,5,6,9], ["B","GR","R","NIR","SWIR_1","BRT","SWIR_2", "pixel_qa"])
                        .filterMetadata('CLOUD_COVER', 'less_than', th_cloud_cover));
        

    return col_SR_L7; 

#GET LANDSAT 5 COLLECTIONS BY PATH ROW
def fexp_landsat_5PathRow(start_date,end_date,n_path, n_row,th_cloud_cover):
    col_SR_L5 =(ee.ImageCollection('LANDSAT/LT05/C01/T1_SR')
                        .filterDate(start_date, end_date)
                        .filterMetadata('WRS_PATH', 'equals', n_path)
                        .filterMetadata('WRS_ROW', 'equals', n_row)
                        .select([0,1,2,3,4,5,6,9], ["B","GR","R","NIR","SWIR_1","BRT","SWIR_2", "pixel_qa"])
                        .filterMetadata('CLOUD_COVER', 'less_than', th_cloud_cover));

    return col_SR_L5;

#GET LANDSAT 7 COLLECTIONS BY COORDINATE
def fexp_landsat_7Coordinate(start_date,end_date,coordinate,th_cloud_cover):
    
    col_SR_L7 =(ee.ImageCollection('LANDSAT/LE07/C01/T1_SR')
                        .filterDate(start_date, end_date)
                        .filterBounds(coordinate)           
                        .select([0,1,2,3,4,5,6,9], ["B","GR","R","NIR","SWIR_1","BRT","SWIR_2", "pixel_qa"])
                        .filterMetadata('CLOUD_COVER', 'less_than', th_cloud_cover));
        

    return col_SR_L7;  

#GET LANDSAT 8 COLLECTIONS BY COORDINATE
def fexp_landsat_8Coordinate(start_date,end_date,coordinate,th_cloud_cover):
    col_SR_L8 =(ee.ImageCollection('LANDSAT/LC08/C01/T1_SR')
                        .filterDate(start_date, end_date)
                        .filterBounds(coordinate)
                        .select([0,1,2,3,4,5,6,7,10],["UB","B","GR","R","NIR","SWIR_1","SWIR_2","BRT","pixel_qa"])
                        .filterMetadata('CLOUD_COVER', 'less_than', th_cloud_cover));
    return col_SR_L8;

#GET LANDSAT 5 COLLECTIONS BY COORDINATE
def fexp_landsat_5Coordinate(start_date,end_date,coordinate,th_cloud_cover):
    col_SR_L5 =(ee.ImageCollection('LANDSAT/LT05/C01/T1_SR')
                        .filterDate(start_date, end_date)
                        .filterBounds(coordinate)
                        .select([0,1,2,3,4,5,6,9], ["B","GR","R","NIR","SWIR_1","BRT","SWIR_2", "pixel_qa"])
                        .filterMetadata('CLOUD_COVER', 'less_than', th_cloud_cover));

    return col_SR_L5;

