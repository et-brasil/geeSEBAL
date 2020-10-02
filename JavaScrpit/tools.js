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
//VERSION 0.1
//CONTATC US: leonardo.laipelt@ufrgs.br
//----------------------------------------------------------------------------------------//
//----------------------------------------------------------------------------------------//
//----------------------------------------------------------------------------------------//

exports.spec_ind = function(image) {

    //NORMALIZED DIFFERENCE VEGETATION INDEX (NDVI)
    var ndvi =  image.normalizedDifference(['NIR', 'R']).rename('NDVI');
    
    //ENHANCED VEGETATION INDEX (EVI)
    var evi = image.expression('2.5 * ((N - R) / (N + (6 * R) - (7.5 * B) + 1))', {
        'N': image.select('NIR').divide(10000),
        'R': image.select('R').divide(10000),
        'B': image.select('B').divide(10000),}).rename('EVI');
    
    //SOIL ADHUSTED VEGETATION INDEX (SAVI)
    var savi = image.expression(
        '((1 + 0.5)*(B5 - B4)) / (0.5 + (B5 + B4))', {
        'B4': image.select('R').multiply(0.0001),
        'B5': image.select('NIR').multiply(0.0001),
    }).rename('SAVI');
 
    var savi1 = savi.where(savi.gt(0.689), 0.689); 
 
    //NORMALIZED DIFFERENCE WATER INDEX (NDWI)
    var ndwi =  image.normalizedDifference(['GR', 'NIR']).rename('NDWI');
  
   //LEAF AREA INDEX (LAI)
    var lai = image.expression(
      '-(log(( 0.69-SAVI)/0.59 ) / 0.91)',
      {'SAVI': savi1}).rename('LAI');
 
    //BROAD-BAND SURFACE EMISSIVITY (e_0)
    var e_0 = image.expression(
      '0.95 + 0.01 * LAI',{
        'LAI': lai});
      e_0 = e_0.where(lai.gt(3), 0.98).rename('e_0');
  
    //NARROW BAND TRANSMISSIVITY (e_NB)
    var e_NB = image.expression(
      '0.97 + (0.0033 * LAI)',{'LAI': lai});
    e_NB = e_NB.where(lai.gt(3), 0.98).rename('e_NB');
    var log_eNB = e_NB.log();
    
  
    //LAND SURFACE TEMPERATURE (LST) [K]
    var comp_onda = ee.Number(1.122e-05); 
    
    var lst = image.expression(
      'Tb / ( 1+ ( ( comp_onda * Tb / fator) * log_eNB))',{
        'Tb': image.select('BRT').multiply(0.1),
        'comp_onda': comp_onda,
        'log_eNB': log_eNB,
        'fator': ee.Number(1.438e-02),
      }).rename('T_LST');
      
    //GET COORDINATES
    var proj = image.select('B').projection();
    var latlon=image.select('B').addBands(ee.Image.pixelLonLat());
    //var latlon = ee.Image.pixelLonLat().reproject(proj);
    var coords = latlon.select(['longitude', 'latitude']);
 
    //FOR FUTHER USE
    var pos_ndvi = ndvi.updateMask(ndvi.gt(0)).rename('pos_NDVI');  
    var ndvi_neg =  pos_ndvi.multiply(-1).rename('NDVI_neg');
    // var lst_neg = lst.multiply(-1).rename('T_LST_neg');
    var int = ee.Image(1).rename('int');
  
  //var no_water = ndwi.updateMask(ndwi.lte(0)).rename('NOWAT');  
 // var lst_nw = lst.updateMask(ndwi.lte(0)).rename('LST_NW');  
  //var lst_neg = lst.multiply(-1).rename('LST_neg');

  image = image.addBands([ndvi,evi,savi,lai,lst,e_0,e_NB,coords,ndvi_neg,pos_ndvi,int,ndwi]);
  return image;

};

//LAND SURFACE TEMPERATURE WITH DEM CORRECTION AND ASPECT/SLOPE
//JAAFAR AND AHMAD (2020)
//PYSEBAL (BASTIAANSSEN) Reference?
exports.LST_correction =function(image,z_alt,T_air,UR,SUN_ELEVATION, hour, min){
  
  //SOLAR CONSTANT [W M-2]
  var    gsc = ee.Number(1367);
  //DAY OF YEAR
  var    dateStr = image.date();
  var    doy = dateStr.getRelative('day', 'year');
  var    Pi=ee.Number(3.14);

  //INVERSE RELATIVE  DISTANCE EARTH-SUN
  var    d1 =  ee.Number(2).multiply(Pi).divide(ee.Number(365));
  var    d2 = d1.multiply(doy);
  var    d3 = d2.cos();
  var    dr = ee.Number(1).add(ee.Number(0.033).multiply(d3));
  
  //ATMOSPHERIC PRESSURE [KPA] 
  //SHUTTLEWORTH (2012)
  var    pres = image.expression(
    '101.3 * ((293 - (0.0065 * Z))/ 293) ** 5.26 ', {
       'Z' : z_alt}).rename('P_ATM');
  
  //SATURATION VAPOR PRESSURE (es) [KPA]
  var    es = image.expression(
     ' 0.6108 *(exp( (17.27 * T_air) / (T_air + 237.3)))', {
       'T_air': T_air}).rename('ES');
  
  //ACTUAL VAPOR PRESSURE (ea) [KPA]
  var    ea = es.multiply(UR).divide(100).rename('EA');
  
  //WATER IN THE ATMOSPHERE [mm]
  //Garrison and Adler (1990)
  var    W = image.expression(  
    '(0.14 * EA * PATM) + 2.1', {
       'PATM' : pres,
       'EA' : ea}).rename('W_ATM');

  //SOLAR ZENITH ANGLE OVER A HORZONTAL SURFACE
  var    solar_zenith = ee.Number(90).subtract(SUN_ELEVATION);
  var    degree2radian = ee.Number(0.01745);
  var    solar_zenith_radians = solar_zenith.multiply(degree2radian);
  var    cos_theta = solar_zenith_radians.cos();  

  //BROAD-BAND ATMOSPHERIC TRANSMISSIVITY (tao_sw)
  //ASCE-EWRI (2005)    
  var    tao_sw = image.expression( 
       '0.35 + 0.627 * exp(((-0.00146 * P)/(Kt * cos_theta)) - (0.075 * (W / cos_theta)**0.4))', {
       'P' : pres,
       'W': W,
       'Kt' : ee.Number(1),
      'cos_theta' : cos_theta,
   }).rename('Tao_sw');

  //AIR DENSITY [KG M-3]    
  var    air_dens = image.expression(
      '(1000* Pair)/(1.01*LST*287)',{
      'Pair': pres,
      'LST': image.select('T_LST')});

  //TEMPERATURE LAPSE RATE (0.0065)
  var Temp_lapse_rate= ee.Number(0.0065) ;

  //LAND SURFACE TEMPERATURE CORRECTION DEM [K]
  var Temp_corr= image.select('T_LST').add(z_alt.select('elevation').multiply(Temp_lapse_rate));
    
  //COS ZENITH ANGLE SUN ELEVATION #ALLEN ET AL. (2006)
  var    slope_aspect = ee.Terrain.products(z_alt);
      
  var    B=(ee.Number(360).divide(ee.Number(365))).multiply(doy.subtract(ee.Number(81)));
  var    delta = ee.Image(ee.Number(23.45).multiply(degree2radian).sin().asin().multiply(B.multiply(degree2radian).sin()));
  var    s =slope_aspect.select('slope').multiply(degree2radian);
  var    gamma = (slope_aspect.select('aspect').subtract(180)).multiply(degree2radian);
  var image_coord = image.addBands(ee.Image.pixelLonLat());
  var    phi = image.select('latitude').multiply(degree2radian);
    
  //CONSTANTS ALLEN ET AL. (2006)
  var    a = ee.Image((delta.sin().multiply(phi.cos()).multiply(s.sin()).multiply(gamma.cos())).subtract(delta.sin().multiply(phi.sin().multiply(s.cos()))));
  var    b = (delta.cos().multiply(phi.cos()).multiply(s.cos())).add(delta.cos().multiply(phi.sin().multiply(s.sin()).multiply(gamma.cos())));
  var    c= (delta.cos().multiply(s.sin()).multiply(gamma.sin()));
  
  //GET IMAGE CENTROID
  var image_center =image.geometry().centroid();
  var longitude_center=ee.Number(image_center.coordinates().get(0));
  
  //DELTA GTM
  var delta_GTM =longitude_center.divide(15).int();
  var min_to_hour=min.divide(60);
  
  //LAND SURFACE TEMPERATURE WITH ASPECT/SLOPE CORRECTION [K]
  var Local_hour_time = hour.add(delta_GTM).add(min_to_hour);
  
  var HRA = (Local_hour_time.subtract(12)).multiply(15);
  
  var w = HRA.multiply(degree2radian) ;
  
  var cos_zn =image.expression(
          '-a +b*w_cos +c*w_sin',{
          'a': a,
          'b': b,
          'c': c,
          'w_cos': w.cos(),
          'w_sin': w.sin()
      }) ;
      
  var    TS_DEM=image.expression(
  
        '(Temp_corr + (Gsc * dr * Transm_corr * cos_zn -Gsc * dr * Transm_corr * cos_zenith_flat) / (air_dens * 1004 * 0.050))',{
        'Temp_corr':Temp_corr,
        'Gsc':gsc,
        'dr':dr,
        'Transm_corr':tao_sw,
        'cos_zenith_flat':cos_theta,
        'cos_zn':cos_zn,
        'air_dens':air_dens
    
        }).rename('T_LST_DEM');

  //MASKS FOR SELECT PRE-CANDIDATES PIXELS   
   var ndwi = image.select('NDWI');
   var lst_nw = image.select('T_LST').updateMask(ndwi.lte(0).and(cos_zn.gt(0.6))).rename('LST_NW'); 
   var lst_neg = image.select('T_LST').multiply(-1).updateMask(cos_zn.gt(0.6)).rename('LST_neg');
    
  image=image.addBands([image.select('T_LST').rename('T_LST_DEM'),lst_nw,lst_neg,cos_zn.rename('Solar_angle_cos')]);
    
  return image;
  
};

  //INSTANTANEOUS OUTGOING LONG-WAVE RADIATION (Rl_up) [W M-2]
exports.fexp_radlong_up = function(image) {
  
  //BROAD-BAND SURFACE THERMAL EMISSIVITY
  //TASUMI ET AL. (2003)
  //ALLEN ET AL. (2007)
  var emi = image.expression(
    '0.95 + (0.01 * LAI)', {'LAI' : image.select('LAI')});
  //LAI
  var lai = image.select('LAI');
  emi = emi.where(lai.gt(3), 0.98);
  var stefBol = ee.Image(5.67e-8);
  
   var Rl_up = image.expression(
      'emi * stefBol * (LST ** 4)', {
      'emi' : emi,
      'stefBol': stefBol,
      'LST': image.select('T_LST'),
    }).rename('Rl_up');
  
  image = image.addBands([Rl_up]);

return image;
};

  //INSTANTANEOUS INCOMING SHORT-WAVE RADIATION (Rs_down) [W M-2]
exports.fexp_radshort_down = function(image, z_alt, T_air, UR,sun_elevation) {
  
  //SOLAR CONSTANT
  var gsc = ee.Number(1367); //[W M-2] 
  
  //DAY OF THE YEAR 
  var dateStr = image.date();
  var doy = dateStr.getRelative('day', 'year');
  
 //INVERSE RELATIVE  DISTANCE EARTH-SUN
  var d1 =  ee.Number(2).multiply(Math.PI).divide(ee.Number(365));
  var d2 = d1.multiply(doy);
  var d3 = d2.cos();
  var dr = ee.Number(1).add(ee.Number(0.033).multiply(d3));

  //ATMOSPHERIC PRESSURE [KPA]
  //SHUTTLEWORTH (2012)
  var pres = image.expression(
       '101.3 * ((293 - (0.0065 * Z))/ 293) ** 5.26 ', {
       'Z' : z_alt,
  }).rename('P_ATM');
    
  //SATURATION VAPOR PRESSURE (es) [KPA]
  var es = image.expression(
       ' 0.6108 *(exp( (17.27 * T_air) / (T_air + 237.3)))', {
       'T_air': T_air,  
  }).rename('ES');
     
  //ACTUAL VAPOR PRESSURE (ea) [KPA]
  var ea = es.multiply(UR).divide(100);
  
  //WATER IN THE ATMOSPHERE [mm]
  //GARRISON AND ADLER (1990)
  var W = image.expression(  
       '(0.14 * EA * PATM) + 2.1', {
       'PATM' : pres,
       'EA' : ea,  
  }).rename('W_ATM');
    
  //SOLAR ZENITH ANGLE OVER A HORIZONTAL SURFACE
  var solar_zenith = ee.Number(90).subtract(sun_elevation);
  var degree2radian = 0.01745;
  var solar_zenith_radians = solar_zenith.multiply(degree2radian);
  var cos_theta = solar_zenith_radians.cos();
 
  //BROAD-BAND ATMOSPHERIC TRANSMISSIVITY (tao_sw) 
  //ASCE-EWRI (2005)
   var tao_sw = image.expression( 
       '0.35 + 0.627 * exp(((-0.00146 * P)/(Kt * cos_theta)) - (0.075 * (W / cos_theta)**0.4))', {
       'P' : pres,
       'W': W,
       'Kt' : ee.Number(1),
       'cos_theta' : image.select('Solar_angle_cos'), 
   }).rename('Tao_sw');

  //INSTANTANEOUS SHORT-WAVE RADIATION (Rs_down) [W M-2]   
  var Rs_down = image.expression( 
   //'(gsc * cos_theta * tao_sw) / (dr * tao_sw)', {
       '(gsc * cos_theta * tao_sw * dr)', {
       'gsc' : gsc,
       'tao_sw': tao_sw,
       'dr' : dr,
       'cos_theta' : image.select('Solar_angle_cos'),  
       
    }).rename('Rs_down');
    
  image =  image.addBands([Rs_down, tao_sw]);
return image;
};

  //INSTANTANEOUS INCOMING LONGWAVE RADIATION (Rl_down) [W M-2]
  //ALLEN ET AL (2007)
exports.fexp_radlong_down = function(image, n_Ts_cold) {
  
  var log_taosw = image.select('Tao_sw').log();
  // (Rl_down=((0.85*(-log(tao_SW))^0.09).*(Ts_cold.^4).*5.67e-8))
  // Incoming longwave radiation 
  var Rl_down = image.expression( 
       '(0.85 * (- log_taosw) ** 0.09) * stefBol * (n_Ts_cold ** 4)  ', {
       'log_taosw' : log_taosw,
       'stefBol': ee.Number(5.67e-8),
       'n_Ts_cold' : n_Ts_cold,  
       
    }).rename('Rl_down');
 
  image = image.addBands([Rl_down]);

return image; 
};

//INSTANTANEOUS NET RADIATON BALANCE (Rn) [W M-2]
exports.fexp_radbalance = function(image) {
  
  // Rn= ((1-alfa)*Rs_down) + RL_down - RL_up  -   (1-e_0).*RL_down;
   var Rn = image.expression( 
       '((1-alfa) * Rs_down) + Rl_down - Rl_up - ((1 - e_0) * Rl_down) ', {
       'alfa' : image.select('ALFA'),
       'Rs_down': image.select('Rs_down'),
       'Rl_down' : image.select('Rl_down'),
       'Rl_up' : image.select('Rl_up'), 
       'e_0' : image.select('e_0'), 
       
    }).rename('Rn');
 
  //return  image.addBands([Rn]);
  image = image.addBands([Rn]);
   
return image;
};
 
  //SOIL HEAT FLUX (G) [W M-2]
  //BASTIAANSSEN (2000)
exports.fexp_soil_heat = function(image) {
  
      var G = image.expression( 
      'Rn * (T_LST - 273.15) * ( 0.0038 + (0.0074 * ALFA)) *  (1 - 0.98 * (NDVI ** 4)) ', {
         'Rn' : image.select('Rn'),
         'NDVI': image.select('NDVI'),
        'ALFA': image.select('ALFA'),
        'T_LST': image.select('T_LST_DEM')}).rename('G');
       
 // var G = image.expression( 
 //   'Rn * (0.05 + 0.18 * exp(-0.521 * LAI))', {
 //      'Rn' : image.select('Rn'),
 //      'LAI': image.select('LAI'),
 //       }).rename('G');     

 // var G1 = image.expression( 
  //  '(0.084 * Rn) + (1.80 * (T_LST - 273.15))', {
 //      'Rn' : image.select('Rn'),
  //     'T_LST': image.select('T_LST'),
  //      }).rename('G1');
  
 // var lai = image.select('LAI');
 //  G = G.where(lai.lt(0.5),G1).rename('G');
 
 // var G = image.expression( 
 //    '(0.3 * Rn) * (1 - 0.98 * (NDVI ** 4)) ', {
 //       'Rn' : image.select('Rn'),
 //       'NDVI': image.select('NDVI')}).rename('G');
       
 //  var G1 = image.expression( 
 //    '(0.5 * Rn) ', {'Rn' : image.select('Rn')});      
       
 //  var ndvi = image.select('NDVI');
 // G = G.where(ndvi.lt(0),G1).rename('G');

 image = image.addBands([G]);

return image;
};

//SENSIBLE HEAT FLUX (H) [W M-2]
  exports.fexp_sensible_heat_flux = function(image, ux, UR, Rn24hobs, n_Ts_cold, d_hot_pixel) {
  
  //VEGETATION HEIGHTS  [M]
  var n_veg_hight = ee.Number(3);
  
  //WIND SPEED AT HEIGHT Zx [M]
  var n_zx = ee.Number(2);
  
  //BLENDING HEIGHT [M]
  var n_hight = ee.Number(200);
  
  //AIR SPECIFIC HEAT [J kg-1/K-1]
  var n_Cp = ee.Number(1004); 
  
  //VON KARMAN'S CONSTANT
  var n_K = ee.Number(0.41);
  
  //TS COLD PIXEL
  n_Ts_cold = ee.Number(n_Ts_cold);
  
  //TS HOT PIXEL
  var n_Ts_hot = ee.Number(d_hot_pixel.get('temp'));
  
  //G HOT PIXEL
  var n_G_hot = ee.Number(d_hot_pixel.get('G'));
  
  //RN HOT PIXEL
  var n_Rn_hot = ee.Number(d_hot_pixel.get('Rn'));
  
  //LAT AND LON HOT PIXEL
  var n_long_hot = ee.Number(d_hot_pixel.get('x'));
  var n_lat_hot = ee.Number(d_hot_pixel.get('y'));
  
  //POINT GEOMETRY
  var p_hot_pix =  ee.Geometry.Point([n_long_hot, n_lat_hot]);

  //SAVI
  var i_savi = image.select('SAVI');
  
  //MOMENTUM ROUGHNESS LENGHT (ZOM) AT THE WEATHER STATION [M]
  //BRUTSAERT (1982)
  var  n_zom = n_veg_hight.multiply(0.12); 
  
  //FRICTION VELOCITY AT WEATHER STATION [M S-1]
  var i_ufric_ws = i_savi.expression( 
      '(n_K * ux)/ log(n_zx /n_zom)', {'n_K': n_K, 'n_zx': n_zx, 'n_zom': n_zom, 'ux': ux }); 

  //WIND SPEED AT BLENDING HEIGHT AT THE WEATHER STATION [M S-1]
  var i_u200 = i_savi.expression( 
      'i_ufric_ws *  (log(n_hight/n_zom)/n_K)', {'i_ufric_ws' : i_ufric_ws, 'n_hight' : n_hight, 'n_zom' : n_zom, 'n_K' : n_K});
      
  //MOMENTUM ROUGHNESS LENGHT (ZOM) FOR EACH PIXEL [M]
  var i_zom = i_savi.expression( 
      'exp((5.62 * (SAVI))-5.809)', {'SAVI' : i_savi,}); 
  
  //FRICTION VELOCITY FOR EACH PIXEL  [M S-1]  
  var i_ufric = i_zom.expression( 
      '(n_K *u200) /(log(hight/i_zom))', {'u200' : i_u200,'hight': n_hight, 'i_zom':n_zom, 'n_K': n_K }).rename('u_fr');
  
  //AERODYNAMIC RESISTANCE TO HEAT TRANSPORT (rah) [S M-1]
  //Z1 AND Z2 ARE HEIGHTS [M] ABOVE THE ZERO PLANE DISPLACEMENT
  //OF THE VEGETATION  
  var z1= ee.Number(0.1);
  var z2= ee.Number(2);   
  
  var i_rah = i_ufric.expression( 
      '(log(z2/z1))/(i_ufric*0.41)', {'z2' : z2,'z1': z1, 'i_ufric':i_ufric }).rename('rah'); 
      
  var i_rah_first = i_rah.rename('rah_first');
  
  //AIR DENSITY HOT PIXEL
  var n_ro_hot= (ee.Number(-0.0046).multiply(n_Ts_hot)).add(ee.Number(2.5538));

  //========ITERATIVE PROCESS=========//

  //SENSIBLE HEAT FLUX AT THE HOT PIXEL (H_hot) 
  var n_H_hot = ee.Number(n_Rn_hot).subtract(ee.Number(n_G_hot));
  
  //ITERATIVE VARIABLES
  var n= ee.Number(1);
  var n_dif= ee.Number(1);
  var n_dif_min = ee.Number(0.1);
  var list_dif = ee.List([]);
  var list_dT_hot = ee.List([]);
  var list_rah_hot = ee.List([]);
  var list_coef_a = ee.List([]);
  var list_coef_b = ee.List([]);
  
  //NUMBER OF ITERATIVE STEPS: 15
  //CAN BE CHANGED, BUT BE AWARE THAT
  //A MINIMUM NUMBER OF ITERATIVE PROCESSES
  //IS NECESSARY TO ACHIEVE RAH AND H ESTIMATIONS

  //========INIT ITERATION========//
  
  for (n = 1; n < 16; n++) {
 
    //AERODYNAMIC RESISTANCE TO HEAT TRANSPORT
    //IN HOT PIXEL
    var d_rah_hot = i_rah.select('rah').reduceRegion({reducer: ee.Reducer.first(), geometry: p_hot_pix, scale: 30,maxPixels: 10e14});

    var n_rah_hot =   ee.Number(d_rah_hot.get('rah'));    
    
    //NEAR SURFACE TEMPERATURE DIFFERENCE IN HOT PIXEL (dT= Tz1-Tz2)  [K]
    //dThot= Hhot*rah/(ρCp)
    var n_dT_hot = (n_H_hot.multiply(n_rah_hot)).divide(n_ro_hot.multiply(n_Cp));
    
    //NEAR SURFACE TEMPERATURE DIFFERENCE IN COLD PIXEL (dT= tZ1-tZ2)
    var n_dT_cold = ee.Number(0);
  
    // dT =  aTs + b
    //ANGULAR COEFFICIENT 
    var n_coef_a = (n_dT_cold.subtract(n_dT_hot)).divide(n_Ts_cold.subtract(n_Ts_hot));
  
    //LINEAR COEFFICIENT 
    var n_coef_b = n_dT_hot.subtract(n_coef_a.multiply(n_Ts_hot));
   
    //dT FOR EACH PIXEL [K]
    var i_lst_med = image.select('T_LST_DEM');
    var i_dT_int = ee.Image(0).clip( image.geometry().bounds()).expression( 
        '(n_coef_a * i_lst_med) + n_coef_b', {
        'n_coef_a' : n_coef_a,
        'n_coef_b': n_coef_b,
        'i_lst_med':i_lst_med }).rename('dT'); 
  
    //AIR TEMPERATURE (TA) FOR EACH PIXEL (TA=TS-dT) [K]
    var i_Ta = i_lst_med.expression( 
      'i_lst_med - i_dT_int', {
      'i_lst_med' : i_lst_med,
      'i_dT_int': i_dT_int});   
    
    //AIR DENSITY (ro) [KM M-3]
    var i_ro = i_Ta.expression( 
      '(-0.0046 * i_Ta) + 2.5538', {
      'i_Ta' : i_Ta}).rename('ro');  
 
    //SENSIBLE HEAT FLUX (H) FOR EACH PIXEL  [W M-2]
    var i_H_int = i_dT_int.expression( 
      '(i_ro*n_Cp*i_dT_int)/i_rah', {
      'i_ro' : i_ro,
      'n_Cp': n_Cp,
      'i_dT_int':i_dT_int,
      'i_rah':i_rah }).rename('H');
      
      
    //MONIN-OBUKHOV LENGTH (L)
    //FOR STABILITY CONDITIONS OF THE ATMOSPHERE IN THE ITERATIVE PROCESS
    var i_L_int = i_dT_int.expression( 
      '-(i_ro*n_Cp*(i_ufric**3)*i_lst_med)/(0.41*9.81*i_H_int)',{
      'i_ro' : i_ro,
      'n_Cp': n_Cp,
      'i_ufric':i_ufric,
      'i_lst_med':i_lst_med,
      'i_H_int':i_H_int }).rename('L');
  
    //STABILITY CORRECTIONS FOR MOMENTUM AND HEAT TRANSPORT
    //PAULSON (1970)
    //WEBB (1970)    
    var img = ee.Image(0).clip( image.geometry().bounds());
    
    //STABILITY CORRECTIONS FOR STABLE CONDITIONS
    var i_psim_200 = img.expression( 
       '-5*(hight/i_L_int)', {
       'hight' : ee.Number(200),
       'i_L_int': i_L_int}).rename('psim_200');
    var i_psih_2 = img.expression( 
      '-5*(hight/i_L_int)',{
      'hight' : ee.Number(2),
      'i_L_int': i_L_int}).rename('psih_2');
  
   var i_psih_01 = img.expression( 
       '-5*(hight/i_L_int)',{
       'hight' : ee.Number(0.1),
       'i_L_int': i_L_int}).rename('psih_01');
  
  //FOR DIFFERENT HEIGHT
    var i_x200 = i_L_int.expression( 
        '(1-(16*(hight/i_L_int)))**0.25',{
        'hight' : ee.Number(200),
        'i_L_int': i_L_int}).rename('i_x200');
    var i_x2 = i_L_int.expression( 
      '(1-(16*(hight/i_L_int)))**0.25',
    {'hight' : ee.Number(2),'i_L_int': i_L_int}).rename('i_x2');
    var i_x01 = i_L_int.expression( 
      '(1-(16*(hight/i_L_int)))**0.25',
    {'hight' : ee.Number(0.1),'i_L_int': i_L_int}); 
  
  //STABILITY CORRECTIONS FOR UNSTABLE CONDITIONS
    var i_psimu_200 = i_x200.expression( 
     '2*log((1+i_x200)/2)+log((1+i_x200**2)/2)-2*atan(i_x200)+0.5*pi',
    {'i_x200' : i_x200,'pi': ee.Number(3.14159265)});
    var i_psihu_2 = i_x2.expression( 
      '2*log((1+i_x2**2)/2)',
    {'i_x2' : i_x2});
    var i_psihu_01 = i_x01.expression( 
     '2*log((1+i_x01**2)/2)',
   {'i_x01' : i_x01});
  
  //FOR EACH PIXEL
    i_psim_200 = i_psim_200.where(i_L_int.lt(0), i_psimu_200);
    i_psih_2 = i_psih_2.where(i_L_int.lt(0), i_psihu_2);
    i_psih_01 = i_psih_01.where(i_L_int.lt(0), i_psihu_01);
  
    i_psim_200 = i_psim_200.where(i_L_int.eq(0), 0);
    i_psih_2 = i_psih_2.where(i_L_int.eq(0), 0);
    i_psih_01 = i_psih_01.where(i_L_int.eq(0), 0);
    
    
    if (n === 1) {
      var i_psim_200_exp = i_psim_200;
      var i_psih_2_exp = i_psih_2;
      var i_psih_01_exp = i_psih_01;
      var i_L_int_exp = i_L_int;
      var i_H_int_exp = i_H_int;
      var i_dT_int_exp = i_dT_int;
      var i_rah_exp = i_rah;
  }
  
    //CORRECTED VALUE FOR THE FRICTION VELOCITY (i_ufric) [M S-1]
    i_ufric = i_ufric.expression( 
      '(u200*0.41)/(log(hight/i_zom)-i_psim_200)',
      {'u200' : i_u200,'hight': n_hight, 'i_zom':i_zom,'i_psim_200': i_psim_200});
      
    //CORRECTED VALUE FOR THE AERODYNAMIC RESISTANCE TO THE HEAT TRANSPORT (rah) [S M-1]
    i_rah = i_rah.expression( 
      '(log(z2/z1)-psi_h2+psi_h01)/(i_ufric*0.41)', 
      {'z2' : z2,'z1': z1, 'i_ufric':i_ufric, 'psi_h2':i_psih_2, 'psi_h01':i_psih_01}).rename('rah');
  
  if (n === 1) {
    var n_dT_hot_old = n_dT_hot;
    var n_rah_hot_old = n_rah_hot;
    n_dif = ee.Number(1);
  }
 
  if (n > 1) {

    var n_dT_hot_abs = n_dT_hot.abs();
    var n_dT_hot_old_abs = n_dT_hot_old.abs();
    var n_rah_hot_abs = n_rah_hot.abs();
    var n_rah_hot_old_abs = n_rah_hot_old.abs();
    n_dif = (n_dT_hot_abs.subtract(n_dT_hot_old_abs).add(n_rah_hot_abs).subtract(n_rah_hot_old_abs)).abs();
    n_dT_hot_old = n_dT_hot;
    n_rah_hot_old = n_rah_hot;
  }
  
  //INSERT EACH ITERATION VALUE INTO A LIST
  list_dif = list_dif.add(n_dif);
  list_coef_a = list_coef_a.add(n_coef_a);
  list_coef_b = list_coef_b.add(n_coef_b);
  list_dT_hot = list_dT_hot.add(n_dT_hot);
  list_rah_hot = list_rah_hot.add(n_rah_hot);

} 
  //=========END ITERATION =========//
  
  //GET FINAL rah, dT AND H
  var i_rah_final = i_rah.rename('rah'); //[SM-1]
  var i_dT_final = i_dT_int.rename('dT'); //[K]
  var i_H_final = i_H_int.expression( //[W M-2]
      '(i_ro*n_Cp*i_dT_int)/i_rah', {
      'i_ro' : i_ro,
      'n_Cp': n_Cp,
      'i_dT_int':i_dT_final,
      'i_rah':i_rah_final }).rename('H');

  image = image.addBands([i_H_final, i_rah_final, i_dT_final, i_rah_first]);
  return image;
  };




