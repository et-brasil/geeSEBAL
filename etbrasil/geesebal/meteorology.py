#----------------------------------------------------------------------------------------#
#---------------------------------------//GEESEBAL//-------------------------------------#
#GEESEBAL - GOOGLE EARTH ENGINE APP FOR SURFACE ENERGY BALANCE ALGORITHM FOR LAND (SEBAL)
#CREATE BY: LEONARDO LAIPELT, RAFAEL KAYSER, ANDERSON RUHOFF AND AYAN FLEISCHMANN
#PROJECT - ET BRASIL https://etbrasil.org/
#LAB - HIDROLOGIA DE GRANDE ESCALA [HGE] website: https://www.ufrgs.br/hge/author/hge/
#UNIVERSITY - UNIVERSIDADE FEDERAL DO RIO GRANDE DO SUL - UFRGS
#RIO GRANDE DO SUL, BRAZIL

#DOI
#VERSION 0.1.1
#CONTACT US: leonardo.laipelt@ufrgs.br

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

def get_meteorology(image,time_start):

    meteo_inst_source = 'ECMWF/ERA5_LAND/HOURLY'

    DATASET = ee.ImageCollection(meteo_inst_source)

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
    i_Rn24_coord =DATASET.first().addBands([ee.Image.pixelLonLat()]);

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
    i_Rs_24h = ee.ImageCollection(meteo_inst_source)\
                .filterDate(ee.Date(time_start).advance(-11,'hour'),ee.Date(time_start).advance(13,'hour'))\
                .select("surface_solar_radiation_downwards_hourly")\
                .sum()\
                .divide(86400).rename('SW_Down')

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

    # AIR TEMPERATURE [K]
    tair_c = NEXT_IMAGE.select('temperature_2m')\
        .subtract(PREVIOUS_IMAGE.select('temperature_2m'))\
        .multiply(DELTA_TIME).add(PREVIOUS_IMAGE.select('temperature_2m'))\
        .rename('AirT_G')

    # WIND SPEED [M S-1]
    wind_u = NEXT_IMAGE.select('u_component_of_wind_10m')\
        .subtract(PREVIOUS_IMAGE.select('u_component_of_wind_10m'))\
        .multiply(DELTA_TIME).add(PREVIOUS_IMAGE.select('u_component_of_wind_10m'))

    wind_v = NEXT_IMAGE.select('v_component_of_wind_10m')\
        .subtract(PREVIOUS_IMAGE.select('v_component_of_wind_10m'))\
        .multiply(DELTA_TIME).add(PREVIOUS_IMAGE.select('v_component_of_wind_10m'))

    # TODO: CGM check if the select calls are needed
    wind_med = wind_u.expression(
        'sqrt(ux_u ** 2 + ux_v ** 2)', {'ux_u': wind_u, 'ux_v': wind_v},
    ).rename('ux_G')

    wind_med = wind_med.expression(
        'ux * (4.87) / log(67.8 * z - 5.42)', {'ux': wind_med, 'z': 10.0}).rename('ux_G')

    # PRESSURE [PA] CONVERTED TO KPA
    tdp = NEXT_IMAGE.select('dewpoint_temperature_2m')\
        .subtract(PREVIOUS_IMAGE.select('dewpoint_temperature_2m'))\
        .multiply(DELTA_TIME).add(PREVIOUS_IMAGE.select('dewpoint_temperature_2m'))\
        .rename('tdp')

    # ACTUAL VAPOR PRESSURE [KPA]
    ea = tdp.expression(
        '0.6108 * (exp((17.27 * T_air) / (T_air + 237.3)))',{
        'T_air': tdp.subtract(273.15)})

    # SATURATED VAPOR PRESSURE [KPA]
    esat = tair_c.expression(
        '0.6108 * (exp((17.27 * T_air) / (T_air + 237.3)))', {'T_air': tair_c.subtract(273.15)})

    # RELATIVE HUMIDITY (%)
    rh = ea.divide(esat).multiply(100).rename('RH_G')

    # Resample
    tair_c = tair_c.subtract(273.15).resample('bilinear')
    wind_med = wind_med.resample('bilinear')
    rh = rh.resample('bilinear')
    swdown24h = i_Rs_24h.resample('bilinear')
    rn24h = i_Rn_24h.resample('bilinear')


    #CONCATENATES IMAGES
    col_meteorology = ee.Image.cat(rn24h, tair_c, rh, wind_med, swdown24h)

    return col_meteorology

if __name__ == "__main__":
    get_meteorology()
