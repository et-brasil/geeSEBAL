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

#GLOBAL LAND DATA ASSIMILATION SYSTEM (GLDAS)
#1984 TO 1999-12-31 - GLDAS 2.0
#2000 TO PRESENT - GLDAS 2.1
#3h, 6h, 9h, 12h, 15h, 18h, 21h, 00h

def get_GLDAS(image,time_start):

    GLDAS20=ee.ImageCollection("NASA/GLDAS/V20/NOAH/G025/T3H").filterDate('1984-01-01','1999-12-31')

    DATASET=GLDAS20.merge(ee.ImageCollection('NASA/GLDAS/V021/NOAH/G025/T3H'))

    #LINEAR INTERPOLATION
    TIME_START_NUM=ee.Number(time_start)
    PREVIOUS_TIME=TIME_START_NUM.subtract(3*60*60*1000)
    NEXT_TIME=TIME_START_NUM.add(3*60*60*1000)

    
    PREVIOUS_IMAGE=(DATASET.filter(ee.Filter.date(PREVIOUS_TIME,TIME_START_NUM))
                          .limit(1, 'system:time_start', False).first())
    
    NEXT_IMAGE=(DATASET.filter(ee.Filter.date(TIME_START_NUM,NEXT_TIME))
                          .limit(1, 'system:time_start', False).first())

    IMAGE_PREVIOUS_TIME= ee.Number(PREVIOUS_IMAGE.get('system:time_start'))

    IMAGE_NEXT_TIME=ee.Number(NEXT_IMAGE.get('system:time_start'))

    DELTA_TIME=(TIME_START_NUM.subtract(IMAGE_PREVIOUS_TIME)).divide(IMAGE_NEXT_TIME.subtract(IMAGE_PREVIOUS_TIME))

    
    #DAILY MEAN VARIABLES
    dataset=DATASET.filter(ee.Filter.date(ee.Date(time_start),ee.Date(time_start).advance(1,'day'))).mean()
    
    #DAY OF THE YEAR
    dateStr = ee.Date(time_start);
    doy = dateStr.getRelative('day', 'year');
    Pi=ee.Number(3.14);

    #INVERSE RELATIVE DISTANCE EARTH-SUN 
    #ALLEN ET AL.(1998)
    d1 =  ee.Number(2).multiply(ee.Number(Pi)).divide(ee.Number(365));
    d2 = d1.multiply(doy);
    d3 = d2.cos();
    dr = ee.Number(1).add(ee.Number(0.033).multiply(d3));

    #SOLAR DECLINATION [RADIANS]
    #ASCE REPORT (2005)
    e1 =  ee.Number(2).multiply(ee.Number(Pi)).multiply(doy);
    e2 = e1.divide(ee.Number(365));
    e3 = e2.subtract(ee.Number(1.39));
    e4 = e3.sin();
    solar_dec = ee.Number(0.409).multiply(e4);

    #GET COORDINATES
    i_Rn24_coord =dataset.addBands([ee.Image.pixelLonLat()]);

    #SUNSET  HOUR ANGLE [RADIANS]
    #ASCE REPORT (2005)
    i_lat_rad = (i_Rn24_coord.select('latitude').multiply(ee.Number(Pi))).divide(ee.Number(180));
    i_sun_hour = i_lat_rad.expression(
    'acos(- tan(lat)* tan(solar_dec))', {
          'lat' : i_lat_rad,
          'solar_dec' : solar_dec}).rename('sun_hour');

    #SOLAR CONSTANT
    gsc = ee.Number(4.92); #[MJ M-2 H-1] 

    #EXTRATERRESTRIAL RADIATION 24H  [MJ M-2 D-1]
    #ASCE REPORT (2005)
    i_Ra_24h = i_sun_hour.expression(
    '(24/pi)*Gcs * dr * ( (omega * sin(lat_rad)* sin(solar_dec)) +  (cos(lat_rad) * cos(solar_dec) * sin(omega)))*11.5740',{
          'pi' : ee.Number(Pi), 
          'Gcs' : gsc,
          'dr': dr,
          'omega': i_sun_hour,
          'solar_dec': solar_dec,
          'lat_rad': i_lat_rad}).rename('Ra_24h');
    
    i_Ra_24h=i_Ra_24h.select('Ra_24h').reduce(ee.Reducer.mean());
   
    #INCOMING SHORT-WAVE RADIATION DAILY EAN [W M-2] 
    i_Rs_24h = dataset.select('SWdown_f_tavg').rename('SW_Down')

    #ALBEDO
    #IF GLDAS
    #i_albedo_ls = dataset.select('Albedo_inst').reduce(ee.Reducer.mean()).divide(100);

    # TASUMI
    i_albedo_ls =image.select('ALFA').first()

    #NET RADIATION 24H [W M-2]
    #BRUIN (1982)
    i_Rn_24h = i_Ra_24h.expression(
    '((1 - albedo) * i_Rs_24h) - (Cs * (i_Rs_24h / i_Ra_24h))',{
           'albedo' : i_albedo_ls,
           'i_Rs_24h' : i_Rs_24h,
           'Cs': ee.Number(110), #CONSTANT
           'i_Ra_24h': i_Ra_24h}).rename('Rn24h_G');


    #GAP FILLING
    d_Rn_24h_med = i_Rn_24h.reduceRegion(
    reducer= ee.Reducer.mean(),
    geometry= image.geometry(),
    scale= 20000,
    maxPixels=10e14);   
    n_Rn_24h_med = ee.Number(d_Rn_24h_med.get('Rn24h_G'));
    i_Rn_24h = i_Rn_24h.unmask(-9999);
    i_Rn_24h =i_Rn_24h.where(i_Rn_24h.eq(-9999), n_Rn_24h_med);

    #RESAMPLE
    i_Rn_24h=i_Rn_24h.clip(image.geometry().bounds()).reduceNeighborhood(
      reducer= ee.Reducer.mean(),
      kernel= ee.Kernel.square(18000,"meters")
      ).reproject(i_Rn_24h.projection().atScale(1000)).resample('bilinear').clip(image.geometry()).rename('Rn24h_G')
    

    #SPECIFIC HUMIDITY [KG KG-1]
    i_q_med =NEXT_IMAGE.select('Qair_f_inst')\
              .subtract(PREVIOUS_IMAGE.select('Qair_f_inst'))\
              .multiply(DELTA_TIME).add(PREVIOUS_IMAGE.select('Qair_f_inst'))

    #AIR TEMPERATURE [K] 
    i_air_temp_K = NEXT_IMAGE.select('Tair_f_inst')\
              .subtract(PREVIOUS_IMAGE.select('Tair_f_inst'))\
              .multiply(DELTA_TIME).add(PREVIOUS_IMAGE.select('Tair_f_inst'))

    #WIND SPEED [M S-1]
    i_ux_med = NEXT_IMAGE.select('Wind_f_inst')\
              .subtract(PREVIOUS_IMAGE.select('Wind_f_inst'))\
              .multiply(DELTA_TIME).add(PREVIOUS_IMAGE.select('Wind_f_inst')).rename('ux_G')

    #PRESSURE [PA] CONVERTED TO KPA
    i_P_med= NEXT_IMAGE.select('Psurf_f_inst')\
              .subtract(PREVIOUS_IMAGE.select('Psurf_f_inst'))\
              .multiply(DELTA_TIME).add(PREVIOUS_IMAGE.select('Psurf_f_inst')).divide(ee.Number(1000))

    #GAP FILLING
    d_ux_med = i_ux_med.reduceRegion(
        reducer= ee.Reducer.mean(),
        geometry= image.geometry(),
        scale=20000,
        maxPixels=10e14);
    n_ux_med = ee.Number(d_ux_med.get('ux_G'));
    i_ux_med = i_ux_med.unmask(-9999);
    i_ux_med = i_ux_med.where(i_ux_med.eq(-9999), n_ux_med);

    #RESAMPLE
    i_ux_med=i_ux_med.clip(image.geometry().bounds()).reduceNeighborhood(
        reducer= ee.Reducer.mean(),
        kernel= ee.Kernel.square(18000,"meters")
      ).reproject(i_ux_med.projection().atScale(1000)).resample('bilinear').clip(image.geometry()).rename('ux_G')
       
    #AIR TEMPERATURE [C]
    i_air_temp_C = i_air_temp_K.subtract(ee.Number(273.15)).rename('AirT_G');

    #GAP FILLING
    d_air_temp_med = i_air_temp_C.reduceRegion(
        reducer= ee.Reducer.mean(),
        geometry= image.geometry(),
        scale= 20000,
        maxPixels=10e14);
    n_air_temp_med = ee.Number(d_air_temp_med.get('AirT_G'));
    i_air_temp_C = i_air_temp_C.unmask(-9999);
    i_air_temp_C =i_air_temp_C.where(i_air_temp_C.eq(-9999), n_air_temp_med);  
    
    #RESAMPLE
    i_air_temp_C=i_air_temp_C.clip(image.geometry().bounds()).reduceNeighborhood(
        reducer= ee.Reducer.mean(),
        kernel= ee.Kernel.square(18000,"meters")
      ).reproject(i_air_temp_C.projection().atScale(1000)).resample('bilinear').clip(image.geometry()).rename('AirT_G')
    
    #ACTUAL VAPOR PRESSURE [KPA]
    i_ea=i_P_med.expression('(1/0.622)*Q*P',{
            'Q': i_q_med,
            'P':i_P_med}).rename('i_ea');

    #SATURATED VAPOR PRESSURE [KPA] 
    i_esat=i_air_temp_C.expression('0.6108*(exp((17.27*T_air)/(T_air+237.3)))',{
        'T_air':i_air_temp_C}).rename('e_sat');

    #RELATIVE HUMIDITY (%)
    #i_RH = (i_q_med.divide(i_q_sat)).multiply(100).rename('RH_G');
    i_RH=i_ea.divide(i_esat).multiply(100).rename('RH_G');

    #GAP FILLING
    d_RH_med = i_RH.reduceRegion(
        reducer= ee.Reducer.mean(),
        geometry= image.geometry(),
        scale= 20000,
        maxPixels=10e14);
    n_RH_med = ee.Number(d_RH_med.get('RH_G'));
    i_RH = i_RH.unmask(-9999);
    i_RH =i_RH.where(i_RH.eq(-9999), n_RH_med);

    #RESAMPLE
    i_RH=i_RH.clip(image.geometry().bounds()).reduceNeighborhood(
        reducer= ee.Reducer.mean(),
        kernel= ee.Kernel.square(18000,"meters")
      ).reproject(i_RH.projection().atScale(1000)).resample('bilinear').clip(image.geometry()).rename('RH_G')
    

    #CLIP 
    #i_Rn_24h=i_Rn_24h.mask(camada_clip);
    #i_air_temp_C=i_air_temp_C.mask(camada_clip);
    #i_RH=i_RH.mask(camada_clip);
    #i_ux_med=i_ux_med.mask(camada_clip);    

    #CONCATENATES IMAGES
    col_GLDAS = ee.Image.cat(i_Rn_24h, i_air_temp_C, i_RH, i_ux_med,i_Rs_24h);
    
    return col_GLDAS;

if __name__ == "__main__":
    get_GLDAS()
