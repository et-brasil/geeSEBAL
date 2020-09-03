//GEESEBAL - GOOGLE EARTH ENGINE APP FOR SURFACE ENERGY BALANCE ALGORITHM FOR LAND (SEBAL)
//CREATE BY: LEONARDO LAIPELT, RAFAEL KAYSER, ANDERSON RUHOFF AND AYAN FLEISCHMANN 
//PROJECT - ET BRASIL
//LAB - HIDROLOGIA DE GRANDE ESCALA [HGE] website: https://www.ufrgs.br/hge/author/hge/
//UNIVERSITY - UNIVERSIDADE FEDERAL DO RIO GRANDE DO SUL - UFRGS
//BRAZIL, RIO GRANDE DO SUL
//
//DOI
//VERSION 1.0
//CONTATC US: leonardo.laipelt@ufrgs.br
exports.get_meteorological_data = function(image, date_init,flag_meteorological) {

//GLDAS
if (flag_meteorological ===0){
    var Rn_Cs= ee.Number(110);
    var start_date=ee.Date(date_init);
    var edate = ee.Date(start_date).advance(1, 'day');
    
    var GLDAS20=ee.ImageCollection("NASA/GLDAS/V20/NOAH/G025/T3H").filterDate('1984-01-01','1999-12-31');
    
    var dataset =GLDAS20.merge(ee.ImageCollection('NASA/GLDAS/V021/NOAH/G025/T3H'));
    
    //get GLDAS collection
    var dataset =dataset
      .filterDate(start_date, edate)
      .filterBounds(image.geometry().bounds());
      
    //albedo from landsat image  
    var i_albedo_ls = image.select('ALFA');
  
    // ---------------------------------------------------  NET DAILY RADIATION -------------------------------------------------------------------
  
    // *** Extraterrestrial Radiation - 24-Hour
    //Day of year
    var i_first = dataset.select('Albedo_inst').first();
    var dateStr = i_first.date();
    var doy = dateStr.getRelative('day', 'year');
    
    //Relative Earth–Sun distance 
    var d1 =  ee.Number(2).multiply(Math.PI).divide(ee.Number(365));
    var d2 = d1.multiply(doy);
    var d3 = d2.cos();
    var dr = ee.Number(1).add(ee.Number(0.033).multiply(d3));
  
    // Solar declination (sun_dec) [radians] - ref.  ASCE Report (2005)
    var e1 =  ee.Number(2).multiply(Math.PI).multiply(doy);
    var e2 = e1.divide(ee.Number(365));
    var e3 = e2.subtract(ee.Number(1.39));
    var e4 = e3.sin();
    var solar_dec = ee.Number(0.409).multiply(e4);
  
    // Get coordinates
   // var i_Rn24_coord = ee.Image.pixelLonLat().reproject(dataset.select('Albedo_inst').first().mask(image.select('R')).projection());
  
    var i_Rn24_coord = i_first.addBands([ee.Image.pixelLonLat()]);
    // Sunset hour angle (i_sun_hour) [radians] - ref.  ASCE Report (2005)
    var i_lat_rad = (i_Rn24_coord.select('latitude').multiply(ee.Number(Math.PI))).divide(ee.Number(180));
    var i_sun_hour = i_lat_rad.expression(
    ' acos(- tan(lat)* tan(solar_dec))', {'lat' : i_lat_rad, 'solar_dec' : solar_dec}).rename('sun_hour');
  
    //Solar constant
    var gsc = ee.Number(4.92); //  ( MJ m-2 h-1) is the solar constant
  
    // Extraterrestrial Radiation for 24-Hour Periods (Ra_24h) - [MJ m-2 d-1] - ref.  ASCE Report (2005)
    var i_Ra_24h = i_sun_hour.expression(
    '(24/pi) * Gcs * dr * ( (omega * sin(lat_rad)* sin(solar_dec)) +  (cos(lat_rad) * cos(solar_dec) * sin(omega))) * 11.5740',
    {'pi' : ee.Number(Math.PI), 'Gcs' : gsc, 'dr': dr, 'omega': i_sun_hour, 'solar_dec': solar_dec, 'lat_rad': i_lat_rad}).rename('Ra_24h');
  
    // INCOMING SHORTWAVE RADIATION DAILY MEAN (W/ m²)
    var i_Rs_24h = dataset.select('SWdown_f_tavg').reduce(ee.Reducer.mean());
    // ALBEDO
    var i_albedo = dataset.select('Albedo_inst').reduce(ee.Reducer.mean()).divide(100);
  
    // *** NET DAILY RADIATION (Rn_24h) - [W/m²]
    var i_Rn_24h = i_Ra_24h.expression(
    '((1 - albedo) * i_Rs_24h) - (Cs * (i_Rs_24h / i_Ra_24h))',
    {'albedo' : i_albedo_ls, 'i_Rs_24h' : i_Rs_24h, 'Cs': ee.Number(Rn_Cs), 'i_Ra_24h': i_Ra_24h}).rename('Rn24h_G');
    
     //Apply a reduce region function
    //var d_Rn_24h_med = i_Rn_24h_init.select('Rn24h').reduceRegion({reducer: ee.Reducer.mean(),geometry: image.geometry().bounds(),scale: 20000,maxPixels: 10e14});
    //var n_Rn_24h_med = ee.Number(d_Rn_24h_med.get('Rn24h'));
  
    //i_Rn_24h = i_Rn_24h_init.unmask(-9999);
   //i_Rn_24h =i_Rn_24h_init.where(i_Rn_24h_init, n_Rn_24h_med).clip(image.select('R').geometry().bounds()).rename('Rn24h_G');
    //var i_Rn_24h = i_Rn_24h.mask(image.select('R')).resample('bilinear').reproject({crs:image.projection().crs(), scale: 250});
   // var i_Rn_24h=i_Rn_24h_init.clip(image.geometry().bounds()).reduceNeighborhood({
   //   reducer: ee.Reducer.mean(),
   //   kernel: ee.Kernel.square(18000,"meters"),
   // }).reproject(i_Rn_24h_init.projection().atScale(1000)).resample('bilinear').clip(image.geometry()).rename('Rn24h_G')


  
    // Image list
    var size_t = dataset.size();
    var list_t = dataset.toList(size_t);
   
    //**** Specific humidity
    var i_q_09 = ee.Image(list_t.get(4)).select('Qair_f_inst'); // 9h
    var i_q_12 = ee.Image(list_t.get(5)).select('Qair_f_inst');  //12h
    var i_q_med = (i_q_09.add(i_q_12)).divide(ee.Number(2));
    
    
    // *** Air temperature
    var i_air_temp_09 = ee.Image(list_t.get(4)).select('Tair_f_inst');   //9h
    var i_air_temp_12 = ee.Image(list_t.get(5)).select('Tair_f_inst');  //12h
    var i_air_temp_K = (i_air_temp_09.add(i_air_temp_12)).divide(ee.Number(2));
  
    // *** Wind speed
    var i_ux_09 = ee.Image(list_t.get(4)).select('Wind_f_inst');  //9h
    var i_ux_12 = ee.Image(list_t.get(5)).select('Wind_f_inst');  //12h
    var i_ux_med = ((i_ux_09.add(i_ux_12)).divide(ee.Number(2))).rename('ux_G');
    //Pressure
    var i_P_09 = ee.Image(list_t.get(4)).select('Psurf_f_inst');  //9h
    var i_P_12 = ee.Image(list_t.get(5)).select('Psurf_f_inst');  //12h
    var i_P_med = ((i_P_09.add(i_P_12)).divide(ee.Number(2)));
}
//
if (flag_meteorological ===1){
  
  
    
    var time_start=ee.Date(date_init);
    var time_start_num=ee.Number(date_init);
    var previous_time_CFSV2 =time_start_num.subtract(6*60*60*1000);
    var next_time_CFSV2 =time_start_num.add(6*60*60*1000);

    var CSFV2_previous_image=ee.Image(ee.ImageCollection('NOAA/CFSV2/FOR6H')
                          .filter(ee.Filter.date(previous_time_CFSV2,time_start_num))
                          .limit(1, 'system:time_start', false).first());

    var CSFV2_next_image=ee.Image(ee.ImageCollection('NOAA/CFSV2/FOR6H')
                          .filter(ee.Filter.date(time_start_num,next_time_CFSV2))
                          .limit(1, 'system:time_start', false).first());             
    
    var CFSV2_image_previous_time= ee.Number(CSFV2_previous_image.get('system:time_start')) ;

    var CFSV2_image_next_time=ee.Number(CSFV2_next_image.get('system:time_start')) ;    
    
    var Delta_Time=(time_start_num.subtract(CFSV2_image_previous_time)).divide(CFSV2_image_next_time.subtract(CFSV2_image_previous_time));
    //Net Radiation
    var Radiation_bands=(ee.ImageCollection('NOAA/CFSV2/FOR6H')
               .filter(ee.Filter.date(ee.Date(time_start),ee.Date(time_start).advance(1,'day')))).mean();
    
   // var size_t = Radiation_bands.size();
  // / var list_t = Radiation_bands.toList(size_t);
   
   // var Radiation_bands=(ee.Image(list_t.get(1)).add(ee.Image(list_t.get(2))).add(ee.Image(list_t.get(3)))).divide(ee.Number(3));
    

    var doy = time_start.getRelative('day', 'year');
    //Relative Earth–Sun distance 
    var d1 =  ee.Number(2).multiply(Math.PI).divide(ee.Number(365));
    var d2 = d1.multiply(doy);
    var d3 = d2.cos();
    var dr = ee.Number(1).add(ee.Number(0.033).multiply(d3));
  
    // Solar declination (sun_dec) [radians] - ref.  ASCE Report (2005)
    var e1 =  ee.Number(2).multiply(Math.PI).multiply(doy);
    var e2 = e1.divide(ee.Number(365));
    var e3 = e2.subtract(ee.Number(1.39));
    var e4 = e3.sin();
    var solar_dec = ee.Number(0.409).multiply(e4);
    
    var i_Rn24_coord = Radiation_bands.addBands([ee.Image.pixelLonLat()]);
    // Sunset hour angle (i_sun_hour) [radians] - ref.  ASCE Report (2005)
    var i_lat_rad = (i_Rn24_coord.select('latitude').multiply(ee.Number(Math.PI))).divide(ee.Number(180));
    var i_sun_hour = i_lat_rad.expression(
    ' acos(- tan(lat)* tan(solar_dec))', {'lat' : i_lat_rad, 'solar_dec' : solar_dec}).rename('sun_hour');
  
    //Solar constant
    var gsc = ee.Number(4.92); //  ( MJ m-2 h-1) is the solar constant
  
    // Extraterrestrial Radiation for 24-Hour Periods (Ra_24h) - [MJ m-2 d-1] - ref.  ASCE Report (2005)
    var i_Ra_24h = i_sun_hour.expression(
    '(24/pi) * Gcs * dr * ( (omega * sin(lat_rad)* sin(solar_dec)) +  (cos(lat_rad) * cos(solar_dec) * sin(omega))) * 11.5740',
    {'pi' : ee.Number(Math.PI), 'Gcs' : gsc, 'dr': dr, 'omega': i_sun_hour, 'solar_dec': solar_dec, 'lat_rad': i_lat_rad}).rename('Ra_24h');
  
    // INCOMING SHORTWAVE RADIATION DAILY MEAN (W/ m²)
    var i_Rs_24h = Radiation_bands.select('Downward_Short-Wave_Radiation_Flux_surface_6_Hour_Average');
    var ALBEDO=image.select("ALFA").reproject(Radiation_bands.projection().atScale(1000)).resample('bilinear').clip(image.geometry());

    // *** NET DAILY RADIATION (Rn_24h) - [W/m²]
    //var i_Rn_24h_init = i_Ra_24h.expression(
    //'((1 - albedo) * i_Rs_24h) - (Cs * (i_Rs_24h / i_Ra_24h))',
    //{'albedo' : ALBEDO, 'i_Rs_24h' : i_Rs_24h, 'Cs': ee.Number(110), 'i_Ra_24h': i_Ra_24h}).rename('Rn24h');   
    
    //var i_Rn_24h_init=Radiation_bands.expression(
    //    '(SWdown - SWup) + (Longdown - Longup)',{
    //        'SWdown':Radiation_bands.select('Downward_Short-Wave_Radiation_Flux_surface_6_Hour_Average'),
    //        'SWup':Radiation_bands.select('Upward_Short-Wave_Radiation_Flux_surface_6_Hour_Average'),
    //        'Longdown':Radiation_bands.select('Downward_Long-Wave_Radp_Flux_surface_6_Hour_Average'),
    //        'Longup':Radiation_bands.select('Upward_Long-Wave_Radp_Flux_surface_6_Hour_Average')
            
    //        }
    //    ).rename('Rn24h');
        
  var Rs=i_Rs_24h.multiply(0.0864).rename('Rs')
    
    var dateStr = image.date()
    var doy = dateStr.getRelative('day', 'year')
    var Pi=ee.Number(3.14159)
    var dr = image.expression(
    '1 + (0.033 * cos((2 * pi/365) * J) )',
    {'J': doy, 'pi': Pi})
     
    var Sd = image.expression(
    "0.40928 * sin( (((2 * pi) / 365) * J) - 1.39 )", 
    {'J': doy, 'pi': Pi})
    
    var convert=Pi.divide(ee.Number((180)))
    
    var lat = i_lat_rad
    
    var Ws = image.expression(
    'acos(-tan(Lat) * tan(Sd))', 
    {'Lat':lat, 'Sd':Sd})
    
    var Gsc = 0.0820 
     
    var Raa = image.expression(
    'Ws * sin(Lat) * sin(Sd) + cos(Lat) * cos(Sd) * sin(Ws)', 
    {'Ws':Ws, 'Lat':lat, 'Sd':Sd}).rename('Raa')
    
    var Ra = image.expression(
    '((24 * 60) / pi) * Gsc * Dr * Raa', 
    {'pi':Pi,'Gsc': Gsc, 'Dr':dr, 'Raa':Raa}).rename('Ra')
  
    var SRTM_ELEVATION ='USGS/SRTMGL1_003';
  var srtm = ee.Image(SRTM_ELEVATION).clip(image.geometry().bounds());
  var elev = srtm.select('elevation');
    
    var Rso = image.expression(
      '(0.75 + 2E-5 * z) * Ra',
      {'z':elev, 'Ra':Ra}
      ).rename('Rso')  
    
    var a = 0.23
    
    var kRs = 0.175  
    
    //Rs = image.expression(
   // 'kRs * sqrt((Tmax - Tmin)) * Ra',
    //{'Tmax': Tmax, 'Tmin': Tmin, 'Ra': Ra, 'kRs': kRs}
    //).rename('Rs')
    
    
    var Rns = image.expression(
    '(1 - a) * Rs',
    {'Rs': Rs, 'a':ALBEDO}).rename('Rns')
  
    var Tmax=(ee.ImageCollection('NOAA/CFSV2/FOR6H')
               .filter(ee.Filter.date(ee.Date(time_start),ee.Date(time_start).advance(1,'day')))).select(5).mean();
    
    var Tmin=(ee.ImageCollection('NOAA/CFSV2/FOR6H')
               .filter(ee.Filter.date(ee.Date(time_start),ee.Date(time_start).advance(1,'day')))).select(7).mean();    
    
    var o = 4.901E-9 
    
    var ea = image.expression(
     ' 0.6108 *(exp( (17.27 * T_air) / (T_air + 237.3)))', {
       'T_air': Tmin.subtract(273.15),  
       }).rename('EA_Min')  
    

               
    var nMaxT=Tmax.reduceRegion({
    reducer:ee.Reducer.mean(),
    geometry:image.geometry(),
    scale:500,  
    maxPixels:1e12});
    
    var Tmax=Tmax.where(Tmax,ee.Number(nMaxT.get('Maximum_temperature_height_above_ground_6_Hour_Interval')));

    var nMinT=Tmin.reduceRegion({
    reducer:ee.Reducer.mean(),
    geometry:image.geometry(),
    scale:500,  
    maxPixels:1e12});
    
    var Tmin=Tmin.where(Tmin,ee.Number(nMinT.get('Minimum_temperature_height_above_ground_6_Hour_Interval')));
    
    var Rnl = image.expression(
    'o * ((Tmax ** 4 + Tmin ** 4)/2.0) * (0.34 - 0.14 * sqrt(ea)) * ((Rs/Rso))', 
    {'o': o, 'Tmax':Tmax, 'Tmin':Tmin, 'ea':ea, 'Rs': Rs, 'Rso': Rso} ).rename('Rnl')
    

    var Rn = image.expression(
    'Rns - Rnl',
    {'Rns': Rns, 'Rnl': Rnl}).rename('Rn_Allen')


    var i_Rn_24h_init=Rn.multiply(11.5).rename('Rn24h');
       
        
    //  'Downward_Short-Wave_Radiation_Flux_surface_6_Hour_Average_mean' - 'Upward_Short-Wave_Radiation_Flux_surface_6_Hour_Average_mean'- "+ 
     //     "(b('Downward_Long-Wave_Radp_Flux_surface_6_Hour_Average_mean') - b('Upward_Long-Wave_Radp_Flux_surface_6_Hour_Average_mean'  
        
        
    //var i_Rn_24h=Radiation_bands.expression(
    //    '((1-ALFA)*SWdown) + (Longdown - Longup)',{
      //      'SWdown':Radiation_bands.select('Downward_Short-Wave_Radiation_Flux_surface_6_Hour_Average'),
      //      'ALFA':ALBEDO,
        //    'Longdown':Radiation_bands.select('Downward_Long-Wave_Radp_Flux_surface_6_Hour_Average'),
        //    'Longup':Radiation_bands.select('Upward_Long-Wave_Radp_Flux_surface_6_Hour_Average')
            
          //  }
       // ).rename('Rn24h_G');
   var i_Rn_24h=i_Rn_24h_init.clip(image.geometry().bounds()).reduceNeighborhood({
   reducer: ee.Reducer.mean(),
    kernel: ee.Kernel.square(15000,"meters"),
   }).reproject(i_Rn_24h_init.projection().atScale(500)).resample('bilinear').clip(image.geometry()).rename('Rn24h_G')      
        
   // var i_Rn_24h=i_Rn_24h.mask(image.select('R')).reproject({crs:image.projection().crs(), scale: 1000})
    //var d_Rn_24h_med = i_Rn_24h.select('Rn24h_G').reduceRegion({reducer: ee.Reducer.mean(),geometry: image.geometry().bounds(),scale: 20000,maxPixels: 10e14});
    //var n_Rn_24h_med = ee.Number(d_Rn_24h_med.get('Rn24h_G'));
    //var i_Rn_24h = i_Rn_24h.mask(image.select('R')).float().resample('bilinear').reproject({crs:image.projection().crs(), scale: 100});
    
   /// i_Rn_24h =i_Rn_24h.where(i_Rn_24h, n_Rn_24h_med).clip(image.select('R').geometry().bounds());
    var i_q_med= CSFV2_next_image.select(12).subtract(CSFV2_previous_image.select(12)).multiply(Delta_Time).add(CSFV2_previous_image.select(12));
    //Air Temeprature K
    var i_air_temp_K= ee.Image(CSFV2_next_image.select(13).subtract(CSFV2_previous_image.select(13)).multiply(Delta_Time).add(CSFV2_previous_image.select(13)));
    //Wind Speed
    var i_ux_u=CSFV2_next_image.select(14).subtract(CSFV2_previous_image.select(14)).multiply(Delta_Time).add(CSFV2_previous_image.select(14)).rename('ux_u');
    var i_ux_v=CSFV2_next_image.select(17).subtract(CSFV2_previous_image.select(17)).multiply(Delta_Time).add(CSFV2_previous_image.select(17)).rename('ux_v');


    var i_ux_med=i_ux_u.expression(
        'sqrt((ux_u)**2 + (ux_v) ** 2) * (4.87/log((67.8 * 10) - 5.42))',{
            'ux_u':i_ux_u.select('ux_u'),
            'ux_v':i_ux_v.select('ux_v')
            }).rename('ux_G');
    //Pressure
    var i_P_med=CSFV2_next_image.select(10).subtract(CSFV2_previous_image.select(10)).multiply(Delta_Time).add(CSFV2_previous_image.select(10));

}
    //Apply a reduce region function
   // var d_ux_med = i_ux_med.select('ux_G').reduceRegion({reducer: ee.Reducer.mean(),geometry: image.geometry().bounds(),scale: 20000,maxPixels: 10e14});
   // var n_ux_med = ee.Number(d_ux_med.get('ux_G'));
    //i_ux_med = i_ux_med.unmask(-9999);
    //i_ux_med = i_ux_med.where(i_ux_med, n_ux_med).clip(image.select('R').geometry().bounds());
    //var i_ux_med_resample = i_ux_med.select('ux_G').float().resample('bilinear').reproject({crs:'EPSG:3857', scale: 100}).rename('ux_G');
    // *** Relative humidity
  
    // Air temperature (celsius)
    var i_air_temp_C = i_air_temp_K.subtract(ee.Number(273.15)).rename('AirT');
    
    var i_air_temp_C=i_air_temp_C.clip(image.geometry().bounds()).reduceNeighborhood({
     reducer: ee.Reducer.mean(),
      kernel: ee.Kernel.square(15000,"meters"),
    }).reproject(i_air_temp_C.projection().atScale(500)).resample('bilinear').clip(image.geometry()).rename('AirT_G')          
    //Apply a reduce region function
    var d_air_temp_med = i_air_temp_C.select('AirT_G').reduceRegion({reducer: ee.Reducer.mean(),geometry: image.geometry().bounds(),scale: 20000,maxPixels: 10e14});
    var n_air_temp_med = ee.Number(d_air_temp_med.get('AirT_G'));
  
    //i_air_temp_C = i_air_temp_C.unmask(-9999);
    i_air_temp_C =i_air_temp_C.where(i_air_temp_C, n_air_temp_med).clip(image.select('R').geometry().bounds());
    //var i_air_temp_C_resample = i_air_temp_C.select('AirT_G').float().resample('bilinear').reproject({crs:'EPSG:3857', scale: 100}).rename('AirT_G');
    // Saturated vapor pressure (esat)
    var i_esat = i_air_temp_C.expression(
    '0.6108 *(exp( (17.27 * T_air) / (T_air + 237.3)))', {'T_air' : i_air_temp_C}).rename('e_sat');
  
    var i_ea=i_P_med.expression('(1/0.622)*Q*P',{
            'Q': i_q_med,
            'P':i_P_med}).rename('i_ea');

    var i_RH=i_ea.divide(i_esat).multiply(100).rename('RH_G');  
  
  
  
    // Density of water vapor in saturated air (ro_sat)
    //var i_ro_sat = i_air_temp_K.expression(
   // '(i_esat * 1000) / (RV * T) ', {'i_esat' : i_esat, 'T' : i_air_temp_K, 'RV': ee.Number(461.51)}).rename('ro_sat');  // RV -  gas constant of water vapor
  
    //  Saturation specific humidity
    //var i_q_sat = i_air_temp_K.expression(
    //'i_ro_sat / (i_ro_sat + ro_dry) ', {'i_ro_sat' : i_ro_sat, 'ro_dry' : ee.Number(1.225)}).rename('q_sat');
  
    // Relative humidity
    //var i_RH = (i_q_med.divide(i_q_sat)).multiply(100).rename('RH_G');
    
    //Apply a reduce region function
    //var d_RH_med = i_RH.select('RH_G').reduceRegion({reducer: ee.Reducer.mean(),geometry: image.geometry().bounds(),scale: 20000,maxPixels: 10e14});
   // var n_RH_med = ee.Number(d_RH_med.get('RH_G'));
  
    //i_RH = i_RH.unmask(-9999);
   //  i_RH =i_RH.where(i_RH, n_RH_med).clip(image.select('R').geometry().bounds());
   // var i_RH_resample = i_RH.select('RH_G').float().resample('bilinear').reproject({crs:'EPSG:3857', scale: 100}).rename('RH_G');
    
    // FINAL COLLECTION
    var col_GLDAS = ee.Image.cat(i_Rn_24h, i_air_temp_C, i_RH, i_ux_med);
    
    var d_col_GLDAS=col_GLDAS.reduceRegion(ee.Reducer.mean(),image.geometry().bounds(),20000)
return col_GLDAS};
