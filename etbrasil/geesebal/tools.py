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

#SPECTRAL INDICES MODULE
def fexp_spec_ind(image):

    #NORMALIZED DIFFERENCE VEGETATION INDEX (NDVI)
    ndvi =  image.normalizedDifference(['NIR', 'R']).rename('NDVI');

    #ENHANCED VEGETATION INDEX (EVI)
    evi = image.expression('2.5 * ((N - R) / (N + (6 * R) - (7.5 * B) + 1))', {
        'N': image.select('NIR').divide(10000),
        'R': image.select('R').divide(10000),
        'B': image.select('B').divide(10000),}).rename('EVI');

    #SOIL ADHUSTED VEGETATION INDEX (SAVI)
    savi = image.expression(
      '((1 + 0.5)*(B5 - B4)) / (0.5 + (B5 + B4))', {
        'B4': image.select('R').divide(10000),
        'B5': image.select('NIR').divide(10000),
    }).rename('SAVI');

    #NORMALIZED DIFFERENCE WATER INDEX (NDWI)
    ndwi =  image.normalizedDifference(['GR', 'NIR']).rename('NDWI');
    savi1 = savi.where(savi.gt(0.689), 0.689);

    #LEAF AREA INDEX (LAI)
    lai = image.expression(
      '-(log(( 0.69-SAVI)/0.59 ) / 0.91)', {'SAVI': savi1}).rename('LAI');

    NDVI_adjust=(ndvi.clamp(0.0, 1.00))
    fipar = (NDVI_adjust.multiply(1).subtract(ee.Number(0.05))).rename('fipar')
    fipar = fipar.clamp(0,1)

    #BROAD-BAND SURFACE EMISSIVITY (e_0)
    e_0 = image.expression(
      '0.95 + 0.01 * LAI',{
        'LAI': lai});
    e_0 = e_0.where(lai.gt(3), 0.98).rename('e_0');

    #NARROW BAND TRANSMISSIVITY (e_NB)
    e_NB = image.expression(
      '0.97 + (0.0033 * LAI)',{'LAI': lai});
    e_NB = e_NB.where(lai.gt(3), 0.98).rename('e_NB');
    log_eNB = e_NB.log();

    #LAND SURFACE TEMPERATURE (LST) [K]
    comp_onda = ee.Number(1.115e-05);
    lst = image.expression(
      'Tb / ( 1+ ( ( comp_onda * Tb / fator) * log_eNB))',{
        'Tb': image.select('BRT').divide(10),
        'comp_onda': comp_onda,
        'log_eNB': log_eNB,
        'fator': ee.Number(1.438e-02),
      }).rename('T_LST');

    #RESCALED BRIGHTNESS TEMPERATURE
    brt_r = image.select('BRT').divide(10).rename('BRT_R');
    proj = image.select('B').projection();
    latlon = ee.Image.pixelLonLat().reproject(proj);
    coords = latlon.select(['longitude', 'latitude']);

    #FOR FUTHER USE
    pos_ndvi = ndvi.updateMask(ndvi.gt(0)).rename('pos_NDVI');
    ndvi_neg =  pos_ndvi.multiply(-1).rename('NDVI_neg');
    #lst_neg = lst.multiply(-1).rename('LST_neg');
    int_ = ee.Image(1).rename('int');
    #lst_nw = lst.updateMask(ndwi.lte(0)).rename('LST_NW');
    sd_ndvi = ee.Image(1).rename('sd_ndvi');

    #ADD BANDS
    image = image.addBands([ndvi,evi,savi,lst,lai,e_0,e_NB,coords,ndvi_neg,pos_ndvi,int_, sd_ndvi, ndwi, brt_r]);
    return image

#LAND SURFACE TEMPERATURE CORRECTION
#JIMENEZ-MUNOZ ET AL. (2009)
#NOT WORKING PROPERLY
def fexp_lst_export(img_main,img_main_RAD,landsat_version,refpoly):
    #NCEP ATMOSPHERIC DATA
    bdate = ee.Date(img_main.get('system:time_start')).format('YYYY-MM-dd');
    edate = ee.Date(bdate).advance(1, 'day');

    #SURFACE WATER VAPOUR VARIABLE FROM NCEP
    #NOTE THAT EACH OBSERVATION DURING DAY
    #IS A COLLECTION THERE ARE 4 OBSERVATION PER DAY
    tairColl_wv = (ee.ImageCollection('NCEP_RE/surface_wv')
                        .filterDate(bdate, edate));

    #CONVERT EACH TIME OF OBSERVATIONS FROM A COLLECTION TO BANDS
    size_w = tairColl_wv.size();
    list_w = tairColl_wv.toList(size_w);

    #SELECT THE IMAGE CORRESPONDING 12H
    # WTR IN NCEP DATA [kg/m^2 ]
    # CONVERT TO  g/cm^2: 1 kg/m2 = 0.1 g/cm^2
    wv = ee.Image(list_w.get(2)).select([0], ['SRWVAP12']).multiply(0.1);

    #CREATING A MEAN OF THE NCEP PRODUCT (PIXEL 2.5 DEGREES)
    d_wv_med = wv.reduceRegion(
    reducer= ee.Reducer.mean(),
    geometry= refpoly,
    scale= 25000,
    maxPixels=9000000000000);
    n_wv_med = ee.Number(d_wv_med.get('SRWVAP12'));
    wv = wv.unmask(-9999);
    wv =wv.where(wv, n_wv_med);

    #THERMAL RADIANCE-AT-THE-SENSOR [W sr−1 m−2 µm−1]
    radtemp = img_main_RAD.select('T_RAD');

    #BRIGHTNESS TEMPERATURE [K]
    brightemp = img_main.select('BRT_R');

    # EMISSIVITY
    e = img_main.select('e_NB');

    #PLANCK'S CONSTANT VALUE [W µm4 m−2 sr−1]
    c1 = ee.Image(1.19104*1e8)

    #CONSTANT (K)
    c2 = (ee.Image(14387.7));

    #CENTRAL WAVELENGHT OF THE THERMAL BAND OF
    #THE LANDSAT SENSOR
    lambda1 = ee.Image(11.457);

    #GAMMA PARASTATIDIS (2017)
    gamma = (radtemp.multiply(c2).divide(brightemp.pow(2))
              .multiply(radtemp.multiply(lambda1.pow(4)).divide(c1)
                        .add(lambda1.pow(-1))).pow(-1));

    # DELTA PARASTATIDIS (2017)
    delta = brightemp.subtract(radtemp.multiply(gamma));

    #PSIX - ATMOSPHERIC FUNCTIONS
    #JIMENEZ-MUNOZ (2009)
    #COEFFICIENT TABLES FOR ATMOSPHERIC PARAMETERIZATION
    #COEFFICIENTS FOR LANDSAT 5 - TIGR61 JIMENEZ-MUNOZ (2009)
    if (landsat_version == 'LANDSAT_5'):
        c11 = 0.08735; c12 = -0.09553; c13 = 1.10188;
        c21 = -0.69188;c22 = -0.58185; c23 = -0.29887;
        c31 = -0.03724; c32 = 1.53065; c33 = -0.45476;

    if (landsat_version == 'LANDSAT_7'):
        c11 = 0.07593; c12 = -0.07132; c13 = 1.08565;
        c21 = -0.61438; c22 = -0.70916; c23 = -0.19379;
        c31 = -0.02892; c32 = 1.46051; c33 = -0.43199;

    # COEFFICIENTS FOR LANDSAT 8 JIMMENEZ-MUNOZ (2014)
    if (landsat_version == 'LANDSAT_8'):
        c11 = 0.04019; c12 = 0.02916; c13 = 1.01523;
        c21 = -0.38333; c22 = -1.50294; c23 = 0.20324;
        c31 = 0.00918; c32 = 1.36072; c33 = -0.27514;

    #CALC PSIX
    psi1 = (ee.Image(c11).multiply(wv.pow(ee.Image(2)))
            .add(ee.Image(c12).multiply(wv))
            .add(ee.Image(c13)));
    psi2 = (ee.Image(c21).multiply(wv.pow(ee.Image(2)))
            .add(ee.Image(c22).multiply(wv))
            .add(ee.Image(c23)));
    psi3 = (ee.Image(c31).multiply(wv.pow(ee.Image(2)))
            .add(ee.Image(c32).multiply(wv))
            .add(ee.Image(c33)));

    # LAND SURFACE TEMPERATURE CORRECTION
    # (Eq. 3 - Jimenez-Munoz, 2009)
    LStemp = (gamma.multiply(((psi1.multiply(radtemp)).add(psi2).divide(e)).add(psi3))
                .add(delta).rename('T_LST'))

    #OTHER TEMPERATURE BANDS
    ndwi = img_main.select('NDWI');
    lst_nw = LStemp.updateMask(ndwi.lte(0)).rename('LST_NW');
    lst_neg = lst_nw.multiply(-1).rename('LST_neg');

    #ADD BANDS
    img_main = img_main.addBands([LStemp,lst_nw,lst_neg]);
    return img_main

#LAND SURFACE TEMPERATURE WITH DEM CORRECTION AND ASPECT/SLOPE
#JAAFAR AND AHMAD (2020)
#PYSEBAL (BASTIAANSSEN) Reference?

def LST_DEM_correction(image, z_alt, T_air, UR,SUN_ELEVATION,hour,minuts):

    #SOLAR CONSTANT [W M-2]
    gsc = ee.Number(1367)

    #DAY OF YEAR
    dateStr = image.date()
    doy = dateStr.getRelative('day', 'year')
    Pi=ee.Number(3.14)

    #INVERSE RELATIVE  DISTANCE EARTH-SUN
    d1 =  ee.Number(2).multiply(Pi).divide(ee.Number(365))
    d2 = d1.multiply(doy)
    d3 = d2.cos()
    dr = ee.Number(1).add(ee.Number(0.033).multiply(d3))

    #ATMOSPHERIC PRESSURE [KPA]
    #SHUTTLEWORTH (2012)
    pres = image.expression(
       '101.3 * ((293 - (0.0065 * Z))/ 293) ** 5.26 ', {
       'Z' : z_alt}).rename('P_ATM')

    #SATURATION VAPOR PRESSURE (es) [KPA]
    es = image.expression(
       ' 0.6108 *(exp( (17.27 * T_air) / (T_air + 237.3)))', {
       'T_air': T_air}).rename('ES')

    #ACTUAL VAPOR PRESSURE (ea) [KPA]
    ea = es.multiply(UR).divide(100).rename('EA')

    #WATER IN THE ATMOSPHERE [mm]
    #Garrison and Adler (1990)
    W = image.expression(
       '(0.14 * EA * PATM) + 2.1', {
       'PATM' : pres,
       'EA' : ea}).rename('W_ATM')

    #SOLAR ZENITH ANGLE OVER A HORZONTAL SURFACE
    solar_zenith = ee.Number(90).subtract(SUN_ELEVATION)
    degree2radian = 0.01745
    solar_zenith_radians = solar_zenith.multiply(degree2radian)
    cos_theta = solar_zenith_radians.cos()

    #BROAD-BAND ATMOSPHERIC TRANSMISSIVITY (tao_sw)
    #ASCE-EWRI (2005)
    tao_sw = image.expression(
       '0.35 + 0.627 * exp(((-0.00146 * P)/(Kt * cos_theta)) - (0.075 * (W / cos_theta)**0.4))', {
       'P' : pres,
       'W': W,
       'Kt' : ee.Number(1),
       'cos_theta' : cos_theta}).rename('Tao_sw')

    #AIR DENSITY [KG M-3]
    air_dens = image.expression(
        '(1000* Pair)/(1.01*LST*287)',{
        'Pair': pres,
        'LST': image.select('T_LST')})

    #TEMPERATURE LAPSE RATE (0.0065)
    Temp_lapse_rate= ee.Number(0.0065)

    #LAND SURFACE TEMPERATURE CORRECTION DEM [K]
    Temp_corr= image.select('T_LST').add(z_alt.select('elevation').multiply(Temp_lapse_rate))

    #COS ZENITH ANGLE SUN ELEVATION #ALLEN ET AL. (2006)
    slope_aspect = ee.Terrain.products(z_alt);

    B=(ee.Number(360).divide(ee.Number(365))).multiply(doy.subtract(ee.Number(81)));
    delta = ee.Image(ee.Number(23.45).multiply(degree2radian).sin().asin().multiply(B.multiply(degree2radian).sin()))
    s =slope_aspect.select('slope').multiply(degree2radian)
    gamma = (slope_aspect.select('aspect').subtract(180)).multiply(degree2radian)
    phi = image.select('latitude').multiply(degree2radian)

    #CONSTANTS ALLEN ET AL. (2006)
    a = ee.Image((delta.sin().multiply(phi.cos()).multiply(s.sin()).multiply(gamma.cos())).subtract(delta.sin().multiply(phi.sin().multiply(s.cos()))));
    b = (delta.cos().multiply(phi.cos()).multiply(s.cos())).add(delta.cos().multiply(phi.sin().multiply(s.sin()).multiply(gamma.cos())));
    c= (delta.cos().multiply(s.sin()).multiply(gamma.sin()))

    #GET IMAGE CENTROID
    image_center =image.geometry().centroid();
    longitude_center=ee.Number(image_center.coordinates().get(0));

    #DELTA GTM
    DELTA_GTM =longitude_center.divide(15).int();

    min_to_hour=minuts.divide(60);

    #LOCAL HOUR TIME
    Local_hour_time = hour.add(DELTA_GTM).add(min_to_hour);

    HOUR_A = (Local_hour_time.subtract(12)).multiply(15);

    w = HOUR_A.multiply(degree2radian) ;

    cos_zn =image.expression(
        '-a +b*w_cos +c*w_sin',{
        'a': a,
        'b': b,
        'c': c,
        'w_cos': w.cos(),
        'w_sin': w.sin()})

    #LAND SURFACE TEMPERATURE WITH ASPECT/SLOPE CORRECTION [K]
    TS_DEM=image.expression(

        '(Temp_corr + (Gsc * dr * Transm_corr * cos_zn -Gsc * dr * Transm_corr * cos_zenith_flat) / (air_dens * 1004 * 0.050))',{
        'Temp_corr':Temp_corr,
        'Gsc':gsc,
        'dr':dr,
        'Transm_corr':tao_sw,
        'cos_zenith_flat':cos_theta,
        'cos_zn':cos_zn,
        'air_dens':air_dens}).rename('T_LST_DEM')

    #MASKS FOR SELECT PRE-CANDIDATES PIXELS
    lst_neg = TS_DEM.multiply(-1).rename('LST_neg');
    ndwi = image.select('NDWI');
    lst_nw = TS_DEM.updateMask(ndwi.lte(0)).rename('LST_NW');

    #ADD BANDS
    image=image.addBands([TS_DEM, lst_neg, lst_nw])

    return image

#INSTANTANEOUS OUTGOING LONG-WAVE RADIATION (Rl_up) [W M-2]
def fexp_radlong_up(image):
    #BROAD-BAND SURFACE THERMAL EMISSIVITY
    #TASUMI ET AL. (2003)
    #ALLEN ET AL. (2007)

   emi = image.expression('0.95 + (0.01 * LAI)', {
           'LAI' : image.select('LAI')})
   #LAI
   lai = image.select('LAI');
   emi = emi.where(lai.gt(3), 0.98)
   stefBol = ee.Image(5.67e-8)

   Rl_up = image.expression(
      'emi * stefBol * (LST ** 4)', {
      'emi' : emi,
      'stefBol': stefBol,
      'LST': image.select('T_LST')}).rename('Rl_up')

   #ADD BANDS
   image = image.addBands([Rl_up])
   return image

#INSTANTANEOUS INCOMING SHORT-WAVE RADIATION (Rs_down) [W M-2]
def fexp_radshort_down(image, z_alt, T_air, UR,SUN_ELEVATION):

    #SOLAR CONSTANT
    gsc = ee.Number(1367) #[W M-2]

    #DAY OF THE YEAR
    dateStr = image.date()
    doy = dateStr.getRelative('day', 'year')
    Pi=ee.Number(3.14)

    #INVERSE RELATIVE  DISTANCE EARTH-SUN
    d1 =  ee.Number(2).multiply(Pi).divide(ee.Number(365));
    d2 = d1.multiply(doy);
    d3 = d2.cos();
    dr = ee.Number(1).add(ee.Number(0.033).multiply(d3));

    #ATMOSPHERIC PRESSURE [KPA]
    #SHUTTLEWORTH (2012)
    pres = image.expression(
       '101.3 * ((293 - (0.0065 * Z))/ 293) ** 5.26 ', {
       'Z' : z_alt}).rename('P_ATM')

    #SATURATION VAPOR PRESSURE (es) [KPA]
    es = image.expression(
       ' 0.6108 *(exp( (17.27 * T_air) / (T_air + 237.3)))', {
       'T_air': T_air}).rename('ES')

    #ACTUAL VAPOR PRESSURE (ea) [KPA]
    ea = es.multiply(UR).divide(100).rename('EA')

    #WATER IN THE ATMOSPHERE [mm]
    #GARRISON AND ADLER (1990)
    W = image.expression(
       '(0.14 * EA * PATM) + 2.1', {
       'PATM' : pres,
       'EA' : ea, }).rename('W_ATM')

    #SOLAR ZENITH ANGLE OVER A HORIZONTAL SURFACE
    solar_zenith = ee.Number(90).subtract(SUN_ELEVATION)
    degree2radian = 0.01745
    solar_zenith_radians = solar_zenith.multiply(degree2radian)
    cos_theta = solar_zenith_radians.cos()

    #BROAD-BAND ATMOSPHERIC TRANSMISSIVITY (tao_sw)
    #ASCE-EWRI (2005)
    tao_sw = image.expression(
       '0.35 + 0.627 * exp(((-0.00146 * P)/(Kt * cos_theta)) - (0.075 * (W / cos_theta)**0.4))', {
       'P' : pres,
       'W': W,
       'Kt' : ee.Number(1),
       'cos_theta' : cos_theta}).rename('Tao_sw')

    #INSTANTANEOUS SHORT-WAVE RADIATION (Rs_down) [W M-2]
    Rs_down = image.expression(
       '(gsc * cos_theta * tao_sw * dr)', {
       'gsc' : gsc,
       'tao_sw': tao_sw,
       'dr' : dr,
       'cos_theta' : cos_theta}).rename('Rs_down')

    #ADD BANDS
    image =  image.addBands([Rs_down, tao_sw,es,ea]);
    return image

    #INSTANTANEOUS INCOMING LONGWAVE RADIATION (Rl_down) [W M-2]
    #ALLEN ET AL (2007)
def fexp_radlong_down(image, n_Ts_cold):

    log_taosw = image.select('Tao_sw').log()
    Rl_down = image.expression(
       '(0.85 * (- log_taosw) ** 0.09) * stefBol * (n_Ts_cold ** 4)  ',{
       'log_taosw' : log_taosw,
       'stefBol': ee.Number(5.67e-8),
       'n_Ts_cold' : n_Ts_cold}).rename('Rl_down')

    #ADD BANDS
    image = image.addBands([Rl_down])
    return image

    #INSTANTANEOUS NET RADIATON BALANCE (Rn) [W M-2]
def fexp_radbalance(image):

    Rn = image.expression(
       '((1-alfa) * Rs_down) + Rl_down - Rl_up - ((1 - e_0) * Rl_down) ', {
       'alfa' : image.select('ALFA'),
       'Rs_down': image.select('Rs_down'),
       'Rl_down' : image.select('Rl_down'),
       'Rl_up' : image.select('Rl_up'),
       'e_0' : image.select('e_0')}).rename('Rn')

    #ADD BANDS
    image = image.addBands(Rn)
    return image

    #SOIL HEAT FLUX (G) [W M-2]
    #BASTIAANSSEN (2000)
def fexp_soil_heat(image):

    G = image.expression(
        'Rn * (T_LST - 273.15) * ( 0.0038 + (0.0074 * ALFA)) *  (1 - 0.98 * (NDVI ** 4)) ', {
         'Rn' : image.select('Rn'),
         'NDVI': image.select('NDVI'),
         'ALFA': image.select('ALFA'),
         'T_LST': image.select('T_LST_DEM')}).rename('G');

    #ADD BANDS
    image = image.addBands([G]);
    return image

    #SENSIBLE HEAT FLUX (H) [W M-2]
def fexp_sensible_heat_flux(image, ux, UR, Rn24hobs, n_Ts_cold, d_hot_pixel, date_string, refpoly):

    #VEGETATION HEIGHTS  [M]
    n_veg_hight = ee.Number(3)

    #WIND SPEED AT HEIGHT Zx [M]
    n_zx = ee.Number(2)

    #BLENDING HEIGHT [M]
    n_hight = ee.Number(200)

    #AIR SPECIFIC HEAT [J kg-1/K-1]
    n_Cp = ee.Number(1004)

    #VON KARMAN'S CONSTANT
    n_K = ee.Number(0.41)

    #TS COLD PIXEL
    n_Ts_cold = ee.Number(n_Ts_cold)
    #TS HOT PIXEL
    n_Ts_hot = ee.Number(d_hot_pixel.get('temp'))
    #G HOT PIXEL
    n_G_hot = ee.Number(d_hot_pixel.get('G').getInfo())
    #RN HOT PIXEL
    n_Rn_hot = ee.Number(d_hot_pixel.get('Rn').getInfo())
    #LAT AND LON HOT PIXEL
    n_long_hot = ee.Number(d_hot_pixel.get('x').getInfo())
    n_lat_hot = ee.Number(d_hot_pixel.get('y').getInfo())
    #POINT GEOMETRY
    p_hot_pix =  ee.Geometry.Point([n_long_hot, n_lat_hot])

    #SAVI
    i_savi = image.select('SAVI')

    #MOMENTUM ROUGHNESS LENGHT (ZOM) AT THE WEATHER STATION [M]
    #BRUTSAERT (1982)
    n_zom = n_veg_hight.multiply(0.12)

    #FRICTION VELOCITY AT WEATHER STATION [M S-1]
    i_ufric_ws = i_savi.expression(
      '(n_K * ux)/ log(n_zx /n_zom)', {
              'n_K': n_K,
              'n_zx': n_zx,
              'n_zom': n_zom,
              'ux': ux }).rename('ufric_ws')

    #WIND SPEED AT BLENDING HEIGHT AT THE WEATHER STATION [M S-1]
    i_u200 = i_savi.expression(
      'i_ufric_ws *  (log(n_hight/n_zom)/n_K)', {
              'i_ufric_ws' : i_ufric_ws,
              'n_hight' : n_hight,
              'n_zom' : n_zom,
              'n_K' : n_K}).rename('i_u200')

    #MOMENTUM ROUGHNESS LENGHT (ZOM) FOR EACH PIXEL [M]
    i_zom = i_savi.expression(
     'exp((5.62 * (SAVI))-5.809)', {
              'SAVI' : i_savi,}).rename('zom')
    #ADD BAND
    image=image.addBands(i_zom.select('zom'))

    #FRICTION VELOCITY FOR EACH PIXEL  [M S-1]
    i_ufric = i_savi.expression(
      '(n_K *u200) /(log(hight/i_zom))', {
              'u200' : i_u200,
              'hight': n_hight,
              'i_zom':n_zom,
              'n_K': n_K }).rename('u_fr')
    #ADD BAND
    image=image.addBands(i_ufric.select('u_fr'))

    #AERODYNAMIC RESISTANCE TO HEAT TRANSPORT (rah) [S M-1]

    #Z1 AND Z2 ARE HEIGHTS [M] ABOVE THE ZERO PLANE DISPLACEMENT
    #OF THE VEGETATION
    z1= ee.Number(0.1);
    z2= ee.Number(2);
    i_rah = i_ufric.expression(
      '(log(z2/z1))/(i_ufric*0.41)', {
              'z2' : z2,
              'z1': z1,
              'i_ufric':i_ufric }).rename('rah')

    i_rah_first = i_rah.rename('rah_first')

    #AIR DENSITY HOT PIXEL
    n_ro_hot= (ee.Number(-0.0046).multiply(n_Ts_hot)).add(ee.Number(2.5538))

    #========ITERATIVE PROCESS=========#

    #SENSIBLE HEAT FLUX AT THE HOT PIXEL (H_hot)
    n_H_hot = ee.Number(n_Rn_hot).subtract(ee.Number(n_G_hot))

    #ITERATIVE VARIABLES
    n= ee.Number(1)
    n_dif= ee.Number(1)
    n_dif_min = ee.Number(0.1)
    list_dif = ee.List([])
    list_dT_hot = ee.List([])
    list_rah_hot = ee.List([])
    list_coef_a = ee.List([])
    list_coef_b = ee.List([])

    #NUMBER OF ITERATIVE STEPS: 15
    #CAN BE CHANGED, BUT BE AWARE THAT
    #A MINIMUM NUMBER OF ITERATIVE PROCESSES
    #IS NECESSARY TO ACHIEVE RAH AND H ESTIMATIONS

    #========INIT ITERATION========#
    for n in range(15):

    #AERODYNAMIC RESISTANCE TO HEAT TRANSPORT
    #IN HOT PIXEL
        d_rah_hot = i_rah.reduceRegion(
            reducer= ee.Reducer.first(),
            geometry= p_hot_pix,
            scale= 30,
            maxPixels=9000000000)

        n_rah_hot =   ee.Number(d_rah_hot.get('rah'))

    #NEAR SURFACE TEMPERATURE DIFFERENCE IN HOT PIXEL (dT= Tz1-Tz2)  [K]
    # dThot= Hhot*rah/(ρCp)
        n_dT_hot = (n_H_hot.multiply(n_rah_hot)).divide(n_ro_hot.multiply(n_Cp))

    #NEAR SURFACE TEMPERATURE DIFFERENCE IN COLD PIXEL (dT= tZ1-tZ2)
        n_dT_cold = ee.Number(0)
    # dT =  aTs + b
    #ANGULAR COEFFICIENT
        n_coef_a = (n_dT_cold.subtract(n_dT_hot)).divide(n_Ts_cold.subtract(n_Ts_hot))

    #LINEAR COEFFICIENT
        n_coef_b = n_dT_hot.subtract(n_coef_a.multiply(n_Ts_hot))

    #dT FOR EACH PIXEL [K]
        i_lst_med = image.select('T_LST_DEM')
        i_dT_int = ee.Image(0).clip(refpoly).expression(
            '(n_coef_a * i_lst_med) + n_coef_b', {
            'n_coef_a' : n_coef_a,
            'n_coef_b': n_coef_b,
            'i_lst_med':i_lst_med }).rename('dT')

    #AIR TEMPERATURE (TA) FOR EACH PIXEL (TA=TS-dT) [K]
        i_Ta = i_lst_med.expression(
            'i_lst_med - i_dT_int', {
            'i_lst_med' : i_lst_med,
            'i_dT_int': i_dT_int})

    #AIR DENSITY (ro) [KM M-3]
        i_ro = i_Ta.expression(
    '(-0.0046 * i_Ta) + 2.5538', {
            'i_Ta' : i_Ta}).rename('ro')

    #SENSIBLE HEAT FLUX (H) FOR EACH PIXEL  [W M-2]
        i_H_int = i_dT_int.expression(
      '(i_ro*n_Cp*i_dT_int)/i_rah', {
              'i_ro' : i_ro,
              'n_Cp': n_Cp,
              'i_dT_int':i_dT_int,
              'i_rah':i_rah }).rename('H')
    #GET VALUE
        d_H_int = i_H_int.reduceRegion(
            reducer= ee.Reducer.first(),
           geometry= p_hot_pix,
            scale= 30,
            maxPixels=9000000000)
        n_H_int =   ee.Number(d_H_int.get('H'))

    #MONIN-OBUKHOV LENGTH (L)
    #FOR STABILITY CONDITIONS OF THE ATMOSPHERE IN THE ITERATIVE PROCESS
        i_L_int = i_dT_int.expression(
                '-(i_ro*n_Cp*(i_ufric**3)*i_lst_med)/(0.41*9.81*i_H_int)',
                {'i_ro' : i_ro,
                 'n_Cp': n_Cp,
                 'i_ufric':i_ufric,
                 'i_lst_med':i_lst_med,
                 'i_H_int':i_H_int }).rename('L')

    #STABILITY CORRECTIONS FOR MOMENTUM AND HEAT TRANSPORT
    #PAULSON (1970)
    #WEBB (1970)
        img = ee.Image(0).clip(refpoly);

    #STABILITY CORRECTIONS FOR STABLE CONDITIONS
        i_psim_200 = img.expression(
                '-5*(hight/i_L_int)', {'hight' : ee.Number(200),'i_L_int': i_L_int}).rename('psim_200')
        i_psih_2 = img.expression(
                '-5*(hight/i_L_int)',{'hight' : ee.Number(2),'i_L_int': i_L_int}).rename('psih_2')
        i_psih_01 = img.expression(
                '-5*(hight/i_L_int)',{'hight' : ee.Number(0.1),'i_L_int': i_L_int}).rename('psih_01')

    #FOR DIFFERENT HEIGHT
        i_x200 = i_L_int.expression(
                '(1-(16*(hight/i_L_int)))**0.25',
                {'hight' : ee.Number(200),'i_L_int': i_L_int}).rename('i_x200')
        i_x2 = i_L_int.expression(
                '(1-(16*(hight/i_L_int)))**0.25',
                {'hight' : ee.Number(2),'i_L_int': i_L_int}).rename('i_x2')
        i_x01 = i_L_int.expression(
                '(1-(16*(hight/i_L_int)))**0.25',
                {'hight' : ee.Number(0.1),'i_L_int': i_L_int})

    #STABILITY CORRECTIONS FOR UNSTABLE CONDITIONS
        i_psimu_200 = i_x200.expression(
                '2*log((1+i_x200)/2)+log((1+i_x200**2)/2)-2*atan(i_x200)+0.5*pi',
                {'i_x200' : i_x200,'pi': ee.Number(3.14159265)})
        i_psihu_2 = i_x2.expression(
                '2*log((1+i_x2**2)/2)',
                {'i_x2' : i_x2})
        i_psihu_01 = i_x01.expression(
                '2*log((1+i_x01**2)/2)',
                {'i_x01' : i_x01})

    #FOR EACH PIXEL
        i_psim_200 = i_psim_200.where(i_L_int.lt(0), i_psimu_200)
        i_psih_2 = i_psih_2.where(i_L_int.lt(0), i_psihu_2)
        i_psih_01 = i_psih_01.where(i_L_int.lt(0), i_psihu_01)
        i_psim_200 = i_psim_200.where(i_L_int.eq(0), 0)
        i_psih_2 = i_psih_2.where(i_L_int.eq(0), 0);
        i_psih_01 = i_psih_01.where(i_L_int.eq(0), 0)

        if n==1:
            i_psim_200_exp = i_psim_200
            i_psih_2_exp = i_psih_2
            i_psih_01_exp = i_psih_01
            i_L_int_exp = i_L_int
            i_H_int_exp = i_H_int
            i_dT_int_exp = i_dT_int
            i_rah_exp = i_rah

    #CORRECTED VALUE FOR THE FRICTION VELOCITY (i_ufric) [M S-1]
        i_ufric = i_ufric.expression(
                '(u200*0.41)/(log(hight/i_zom)-i_psim_200)',{
                 'u200' : i_u200,
                 'hight': n_hight,
                 'i_zom':i_zom,
                 'i_psim_200': i_psim_200}).rename('ufric_star')

    #CORRECTED VALUE FOR THE AERODYNAMIC RESISTANCE TO THE HEAT TRANSPORT (rah) [S M-1]
        i_rah = i_rah.expression(
                '(log(z2/z1)-psi_h2+psi_h01)/(i_ufric*0.41)',
                {'z2' : z2,'z1': z1, 'i_ufric':i_ufric, 'psi_h2':i_psih_2, 'psi_h01':i_psih_01}).rename('rah')
        if n==1:
            n_dT_hot_old = n_dT_hot
            n_rah_hot_old = n_rah_hot
            n_dif = ee.Number(1)

        if n > 1:
            n_dT_hot_abs = n_dT_hot.abs()
            n_dT_hot_old_abs = n_dT_hot_old.abs()
            n_rah_hot_abs = n_rah_hot.abs()
            n_rah_hot_old_abs = n_rah_hot_old.abs()
            n_dif=(n_dT_hot_abs.subtract(n_dT_hot_old_abs).add(n_rah_hot_abs).subtract(n_rah_hot_old_abs)).abs()
            n_dT_hot_old = n_dT_hot
            n_rah_hot_old = n_rah_hot

        #INSERT EACH ITERATION VALUE INTO A LIST
        list_dif = list_dif.add(n_dif);
        list_coef_a = list_coef_a.add(n_coef_a)
        list_coef_b = list_coef_b.add(n_coef_b)
        list_dT_hot = list_dT_hot.add(n_dT_hot)
        list_rah_hot = list_rah_hot.add(n_rah_hot)

    #=========END ITERATION =========#

    #GET FINAL rah, dT AND H
    i_rah_final = i_rah.rename('rah') #[SM-1]
    i_dT_final = i_dT_int.rename('dT') #[K]
    i_H_final = i_H_int.expression(  #[W M-2]
            '(i_ro*n_Cp*i_dT_int)/i_rah',{
             'i_ro' : i_ro,
             'n_Cp': n_Cp,
             'i_dT_int':i_dT_final,
             'i_rah':i_rah_final }).rename('H')

    #ADD BANDS
    image = image.addBands([i_H_final, i_rah_final, i_dT_final,
                            i_rah_first,image.select('zom'),image.select('u_fr'),i_ufric])
    return image
