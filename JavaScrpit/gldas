//-------------------------------------------------------------------------------------------//
//----------------------------------------//GEESEBAL//---------------------------------------//
//GEESEBAL - GOOGLE EARTH ENGINE APP FOR SURFACE ENERGY BALANCE ALGORITHM FOR LAND (SEBAL)
//CREATE BY: LEONARDO LAIPELT, RAFAEL KAYSER, ANDERSON RUHOFF AND AYAN FLEISCHMANN 
//PROJECT - ET BRASIL https://etbrasil.org/
//LAB - HIDROLOGIA DE GRANDE ESCALA [HGE] website: https://www.ufrgs.br/hge/author/hge/
//UNIVERSITY - UNIVERSIDADE FEDERAL DO RIO GRANDE DO SUL - UFRGS
//BRAZIL, RIO GRANDE DO SUL
//
//DOI
//VERSION 1.0
//CONTATC US: leonardo.laipelt@ufrgs.br

exports.gldas = function(image,time_start) {

    
    var GLDAS20=ee.ImageCollection("NASA/GLDAS/V20/NOAH/G025/T3H").filterDate('1984-01-01','1999-12-31');
    
    var DATASET=GLDAS20.merge(ee.ImageCollection('NASA/GLDAS/V021/NOAH/G025/T3H'));
    
    //LINEAR INTERPOLATION
    var TIME_START_NUM=ee.Number(time_start);
    var PREVIOUS_TIME=TIME_START_NUM.subtract(3*60*60*1000);
    var NEXT_TIME=TIME_START_NUM.add(3*60*60*1000);    
    
        
    var PREVIOUS_IMAGE=(DATASET.filter(ee.Filter.date(PREVIOUS_TIME,TIME_START_NUM))
                          .limit(1, 'system:time_start', false).first());
    
    var NEXT_IMAGE=(DATASET.filter(ee.Filter.date(TIME_START_NUM,NEXT_TIME))
                          .limit(1, 'system:time_start', false).first());

    var IMAGE_PREVIOUS_TIME= ee.Number(PREVIOUS_IMAGE.get('system:time_start'));

    var IMAGE_NEXT_TIME=ee.Number(NEXT_IMAGE.get('system:time_start'));

    var DELTA_TIME=(TIME_START_NUM.subtract(IMAGE_PREVIOUS_TIME)).divide(IMAGE_NEXT_TIME.subtract(IMAGE_PREVIOUS_TIME));

    //DAILY MEAN VARIABLES
    var dataset=DATASET.filter(ee.Filter.date(ee.Date(time_start),ee.Date(time_start).advance(1,'day'))).mean();
    
    //DAY OF THE YEAR
    var dateStr = ee.Date(time_start);
    var doy = dateStr.getRelative('day', 'year');
    
    //INVERSE RELATIVE DISTANCE EARTH-SUN
    //ALLEN ET AL.(1998)
    var d1 =  ee.Number(2).multiply(Math.PI).divide(ee.Number(365));
    var d2 = d1.multiply(doy);
    var d3 = d2.cos();
    var dr = ee.Number(1).add(ee.Number(0.033).multiply(d3));
  
    //SOLAR DECLINATION [RADIANS]
    //ASCE REPORT (2005)
    var e1 =  ee.Number(2).multiply(Math.PI).multiply(doy);
    var e2 = e1.divide(ee.Number(365));
    var e3 = e2.subtract(ee.Number(1.39));
    var e4 = e3.sin();
    var solar_dec = ee.Number(0.409).multiply(e4);
  
    //GET COORDINATES
    var i_Rn24_coord = dataset.addBands([ee.Image.pixelLonLat()]);
    
    //SUNSET  HOUR ANGLE [RADIANS]
    //ASCE REPORT (2005)    
    var i_lat_rad = (i_Rn24_coord.select('latitude').multiply(ee.Number(Math.PI))).divide(ee.Number(180));
    var i_sun_hour = i_lat_rad.expression(
    ' acos(- tan(lat)* tan(solar_dec))', {
      'lat' : i_lat_rad,
      'solar_dec' : solar_dec}).rename('sun_hour');
  
    //SOLAR CONSTANT
    var gsc = ee.Number(4.92); //[MJ M-2 H-1] 
  
    //EXTRATERRESTRIAL RADIATION 24H  [MJ M-2 D-1]
    //ASCE REPORT (2005)
    var i_Ra_24h = i_sun_hour.expression(
    '(24/pi) * Gcs * dr * ( (omega * sin(lat_rad)* sin(solar_dec)) +  (cos(lat_rad) * cos(solar_dec) * sin(omega))) * 11.5740',{
      'pi' : ee.Number(Math.PI),
      'Gcs' : gsc,
      'dr': dr, 'omega': i_sun_hour, 
      'solar_dec': solar_dec, 
      'lat_rad': i_lat_rad}).rename('Ra_24h');
  
    // INCOMING SHORT-WAVE RADIATION DAILY EAN [W M-2] 
    var i_Rs_24h = dataset.select('SWdown_f_tavg').reduce(ee.Reducer.mean());
    
    // ALBEDO
    //IF GLDAS
    //var i_albedo = dataset.select('Albedo_inst').reduce(ee.Reducer.mean()).divide(100);
    //TASUMI
    var i_albedo =image.select('ALFA');
    
    //NET RADIATION 24H [W M-2]
    //BRUIN (1982)

    var i_Rn_24h = i_Ra_24h.expression(
    '((1 - albedo) * i_Rs_24h) - (Cs * (i_Rs_24h / i_Ra_24h))',{
      'albedo' : i_albedo, 
      'i_Rs_24h' : i_Rs_24h, 
      'Cs': ee.Number(110), 
      'i_Ra_24h': i_Ra_24h}).rename('Rn24h_G');
    
     //GAP FILLING
    var d_Rn_24h_med = i_Rn_24h.select('Rn24h_G').reduceRegion({reducer: ee.Reducer.mean(),geometry: image.geometry().bounds(),scale: 20000,maxPixels: 10e14});
    var n_Rn_24h_med = ee.Number(d_Rn_24h_med.get('Rn24h_G'));
    var i_Rn_24h = i_Rn_24h.unmask(-9999);
    var i_Rn_24h =i_Rn_24h.where(i_Rn_24h.eq(-9999), n_Rn_24h_med).clip(image.select('R').geometry().bounds());
    
    //RESAMPLE
    var i_Rn_24h=i_Rn_24h.clip(image.geometry().bounds()).reduceNeighborhood({
      reducer: ee.Reducer.mean(),
      kernel: ee.Kernel.square(18000,"meters"),
    }).reproject(i_Rn_24h.projection().atScale(1000)).resample('bilinear').clip(image.geometry()).rename('Rn24h_G')


    //SPECIFIC HUMIDITY [KG KG-1]
    var i_q_med =NEXT_IMAGE.select('Qair_f_inst')
              .subtract(PREVIOUS_IMAGE.select('Qair_f_inst'))
              .multiply(DELTA_TIME).add(PREVIOUS_IMAGE.select('Qair_f_inst'));

    //AIR TEMPERATURE [K] 
    var i_air_temp_K = NEXT_IMAGE.select('Tair_f_inst')
              .subtract(PREVIOUS_IMAGE.select('Tair_f_inst'))
              .multiply(DELTA_TIME).add(PREVIOUS_IMAGE.select('Tair_f_inst'));

    //WIND SPEED [M S-1]
    var i_ux_med = NEXT_IMAGE.select('Wind_f_inst')
              .subtract(PREVIOUS_IMAGE.select('Wind_f_inst'))
              .multiply(DELTA_TIME).add(PREVIOUS_IMAGE.select('Wind_f_inst')).rename('ux_G');

    //PRESSURE [PA] CONVERTED TO KPA
    var i_P_med= NEXT_IMAGE.select('Psurf_f_inst')
              .subtract(PREVIOUS_IMAGE.select('Psurf_f_inst'))
              .multiply(DELTA_TIME).add(PREVIOUS_IMAGE.select('Psurf_f_inst')).divide(ee.Number(1000));


    //GAP FILLING
    var d_ux_med = i_ux_med.reduceRegion({
        reducer: ee.Reducer.mean(),
        geometry: image.geometry(),
        scale:20000,
        maxPixels:10e14});
    var n_ux_med = ee.Number(d_ux_med.get('ux_G'));
    var i_ux_med = i_ux_med.unmask(-9999);
    var i_ux_med = i_ux_med.where(i_ux_med.eq(-9999), n_ux_med);
    
    //RESAMPLE
    var i_ux_med=i_ux_med.clip(image.geometry().bounds()).reduceNeighborhood({
      reducer: ee.Reducer.mean(),
      kernel: ee.Kernel.square(18000,"meters"),
    }).reproject(i_ux_med.projection().atScale(1000)).resample('bilinear').clip(image.geometry()).rename('ux_G');
    
    //AIR TEMPERATURE [C]
    var i_air_temp_C = i_air_temp_K.subtract(ee.Number(273.15)).rename('AirT_G');

    //GAP FILLING
    var d_air_temp_med = i_air_temp_C.reduceRegion({
        reducer: ee.Reducer.mean(),
        geometry: image.geometry(),
        scale: 20000,
        maxPixels:10e14});
    var n_air_temp_med = ee.Number(d_air_temp_med.get('AirT_G'));
    var i_air_temp_C = i_air_temp_C.unmask(-9999);
    var i_air_temp_C =i_air_temp_C.where(i_air_temp_C.eq(-9999), n_air_temp_med);  
    
    //RESAMPLE
    var i_air_temp_C=i_air_temp_C.clip(image.geometry().bounds()).reduceNeighborhood({
      reducer: ee.Reducer.mean(),
      kernel: ee.Kernel.square(18000,"meters"),
    }).reproject(i_air_temp_C.projection().atScale(1000)).resample('bilinear').clip(image.geometry()).rename('AirT_G');
    
    //ACTUAL VAPOR PRESSURE [KPA]
    var i_ea=i_P_med.expression('(1/0.622)*Q*P',{
            'Q': i_q_med,
            'P':i_P_med}).rename('i_ea');

    //SATURATED VAPOR PRESSURE [KPA] 
    var i_esat=i_air_temp_C.expression('0.6108*(exp((17.27*T_air)/(T_air+237.3)))',{
        'T_air':i_air_temp_C}).rename('e_sat');

    //RELATIVE HUMIDITY (%)
    //i_RH = (i_q_med.divide(i_q_sat)).multiply(100).rename('RH_G');
    var i_RH=i_ea.divide(i_esat).multiply(100).rename('RH_G');

    //GAP FILLING
    var d_RH_med = i_RH.reduceRegion({
        reducer: ee.Reducer.mean(),
        geometry: image.geometry(),
        scale: 20000,
        maxPixels:10e14});
    var n_RH_med = ee.Number(d_RH_med.get('RH_G'));
    var i_RH = i_RH.unmask(-9999);
    var i_RH =i_RH.where(i_RH.eq(-9999), n_RH_med);
    
    //RESAMPLE
    var i_RH=i_RH.clip(image.geometry().bounds()).reduceNeighborhood({
      reducer: ee.Reducer.mean(),
      kernel: ee.Kernel.square(18000,"meters"),
    }).reproject(i_RH.projection().atScale(1000)).resample('bilinear').clip(image.geometry()).rename('RH_G');
    
    // FINAL COLLECTION
    var col_GLDAS = ee.Image.cat(i_Rn_24h, i_air_temp_C, i_RH, i_ux_med);
    
return col_GLDAS};
