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
exports.LST_correction =function(image,z_alt,T_air,UR,SUN_ELEVATION, hour, min){
  
var    gsc = ee.Number(1367) //%(W/m²) is the solar constant
var    dateStr = image.date()
var    doy = dateStr.getRelative('day', 'year')
var    Pi=ee.Number(3.14)
    //Relative Earth–Sun distance 
var    d1 =  ee.Number(2).multiply(Pi).divide(ee.Number(365));
var    d2 = d1.multiply(doy);
var    d3 = d2.cos();
var    dr = ee.Number(1).add(ee.Number(0.033).multiply(d3));
    //Atmospheric pressure (kPa) - ASCE-EWRI (2005)
var    pres = image.expression(
    '101.3 * ((293 - (0.0065 * Z))/ 293) ** 5.26 ', {
       'Z' : z_alt,
  }).rename('P_ATM');
    //Saturation Vapor Pressure (es)
var    es = image.expression(
     ' 0.6108 *(exp( (17.27 * T_air) / (T_air + 237.3)))', {
       'T_air': T_air,  
  }).rename('ES')  ;
    //Actual Vapor Pressure (ea) - ea=es*UR/100
var    ea = es.multiply(UR).divide(100).rename('EA');
    //Water in the atmosphere
    // (Garrison and Adler (1990)
var    W = image.expression(  
    '(0.14 * EA * PATM) + 2.1', {
       'PATM' : pres,
       'EA' : ea, // near-surface vapor pressure
  }).rename('W_ATM')      ;
//Solar zenith angle over a horizontal surface
var    solar_zenith = ee.Number(90).subtract(SUN_ELEVATION);
var    degree2radian = ee.Number(0.01745);
var    solar_zenith_radians = solar_zenith.multiply(degree2radian);
var    cos_theta = solar_zenith_radians.cos()   ;  
    //Broad-band atmospheric transmissivity (tao_sw)
    //ASCE-EWRI (2005)
    
var    tao_sw = image.expression( 
     '0.35 + 0.627 * exp(((-0.00146 * P)/(Kt * cos_theta)) - (0.075 * (W / cos_theta)**0.4))', {
       'P' : pres,
       'W': W,
       'Kt' : ee.Number(1),
      'cos_theta' : cos_theta,
   }).rename('Tao_sw');
    
var    air_dens = image.expression(
        '(1000* Pair)/(1.01*LST*287)',{
            'Pair': pres,
            'LST': image.select('T_LST'),
            });

var Temp_lapse_rate= ee.Number(0.0065) ;

var    Temp_corr= image.select('T_LST').add(z_alt.select('elevation').multiply(Temp_lapse_rate));
    

var    slope_aspect = ee.Terrain.products(z_alt);
    
var    B=(ee.Number(360).divide(ee.Number(365))).multiply(doy.subtract(ee.Number(81)));
var    delta = ee.Image(ee.Number(23.45).multiply(degree2radian).sin().asin().multiply(B.multiply(degree2radian).sin()));
var    s =slope_aspect.select('slope').multiply(degree2radian);
var    gamma = (slope_aspect.select('aspect').subtract(180)).multiply(degree2radian);

var image_coord = image.addBands(ee.Image.pixelLonLat());

var    phi = image.select('latitude').multiply(degree2radian);

var    a = ee.Image((delta.sin().multiply(phi.cos()).multiply(s.sin()).multiply(gamma.cos())).subtract(delta.sin().multiply(phi.sin().multiply(s.cos()))));

var    b = (delta.cos().multiply(phi.cos()).multiply(s.cos())).add(delta.cos().multiply(phi.sin().multiply(s.sin()).multiply(gamma.cos())));
var    c= (delta.cos().multiply(s.sin()).multiply(gamma.sin()));

var centro_imagem =image.geometry().centroid();

var longitude_center=ee.Number(centro_imagem.coordinates().get(0));

var delta_GTM =longitude_center.divide(15).int();

var min_to_hour=min.divide(60);

var Local_hour_time = hour.add(delta_GTM).add(min_to_hour);

var HRA = (Local_hour_time.subtract(12)).multiply(15);

var w = HRA.multiply(degree2radian) ;

var    cos_zn =image.expression(
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
   
   var ndwi = image.select('NDWI');
   var lst_nw = image.select('T_LST').updateMask(ndwi.lte(0).and(cos_zn.gt(0.6))).rename('LST_NW'); 
   var lst_neg = image.select('T_LST').multiply(-1).updateMask(cos_zn.gt(0.6)).rename('LST_neg');
    
    image=image.addBands([image.select('T_LST').rename('T_LST_DEM'),lst_nw,lst_neg,cos_zn.rename('Solar_angle_cos')]);
    
    return image;
  
};

// -----------------------------------------------------------------------------------------------
// OUTGOING LONGWAVE RADIATION (Rl_up) 
exports.fexp_radlong_up = function(image) {
  
  // Broad-band surface thermal emissivity (e_0 = 0.95 + 0.01 LAI for LAI lte 3)
  // Tasumi (2003) apud Allen (2007)
  var emi = image.expression(
    '0.95 + (0.01 * LAI)', {'LAI' : image.select('LAI')});
    
  var lai = image.select('LAI');
  emi = emi.where(lai.gt(3), 0.98);
 
 
  var stefBol = ee.Image(5.67e-8);
  
   var Rl_up = image.expression(
    'emi * stefBol * (LST ** 4)', {
      'emi' : emi,
      'stefBol': stefBol,
      'LST': image.select('T_LST'),
    }).rename('Rl_up');
  
  //return image.addBands([Rl_up]);
  image = image.addBands([Rl_up]);
  
  //}  // end function

return image;
};

// -------------------------------------------------------------------------------------------------------
// INCOMING SHORTWAVE RADIATION (Rs_down)(W/m²)
exports.fexp_radshort_down = function(image, z_alt, T_air, UR,sun_elevation) {
  
  //Solar constant
  var gsc = ee.Number(1367); //%(W/m²) is the solar constant
  
  //Day of year
  var dateStr = image.date();
  var doy = dateStr.getRelative('day', 'year');
  
 //Relative Earth–Sun distance 
  var d1 =  ee.Number(2).multiply(Math.PI).divide(ee.Number(365));
  var d2 = d1.multiply(doy);
  var d3 = d2.cos();
  var dr = ee.Number(1).add(ee.Number(0.033).multiply(d3));

  
  //Atmospheric pressure (kPa) - ASCE-EWRI (2005)
  var pres = image.expression(
    '101.3 * ((293 - (0.0065 * Z))/ 293) ** 5.26 ', {
       'Z' : z_alt,
  }).rename('P_ATM');
    

  // Saturation Vapor Pressure (es)
  var es = image.expression(
     ' 0.6108 *(exp( (17.27 * T_air) / (T_air + 237.3)))', {
       'T_air': T_air,  
  }).rename('ES');
     
  
  // Actual Vapor Pressure (ea) - ea=es*UR/100
  var ea = es.multiply(UR).divide(100);
  
  // Water in the atmosphere
  // [W unit - milimeters] (Garrison and Adler (1990)
  var W = image.expression(  
    '(0.14 * EA * PATM) + 2.1', {
       'PATM' : pres,
       'EA' : ea,  //near-surface vapor pressure
  }).rename('W_ATM');
    
    
  // Solar zenith angle over a horizontal surface
  var solar_zenith = ee.Number(90).subtract(sun_elevation);
  var degree2radian = 0.01745;
  var solar_zenith_radians = solar_zenith.multiply(degree2radian);
  var cos_theta = solar_zenith_radians.cos();
 
  // Broad-band atmospheric transmissivity (tao_sw)
  // ASCE-EWRI (2005)
   var tao_sw = image.expression(  //
     '0.35 + 0.627 * exp(((-0.00146 * P)/(Kt * cos_theta)) - (0.075 * (W / cos_theta)**0.4))', {
       'P' : pres,
       'W': W,
       'Kt' : ee.Number(1),
      'cos_theta' : image.select('Solar_angle_cos'),  //near-surface vapor pressure
   }).rename('Tao_sw');

  
  // Incoming shortwave radiation 
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

//---------------------------------------------------------------------------------------------------------
// INCOMING LONGWAVE RADIATION (Rl_down) (W/m²)
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
 
 //return  image.addBands([Rl_down]);
 image = image.addBands([Rl_down]);
  //}  // end function
  
return image; 
};

// --------------------------------------------------------------------------------------------------------
// NET RADIATON BALANCE (Rn) (W/m²)
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
 

// -------------------------------------------------------------------------------------------------------
// SOIL HEAT FLUX (G) (W/m²)
exports.fexp_soil_heat = function(image) {
  
  //Bastianssen (2000)
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

// ----------------------------------------- COLD PIXEL -----------------------------------------------------------------------
exports.fexp_cold_pixel = function(image, p_top_NDVI, p_coldest_Ts) {

  var d_perc_top_NDVI = image.select('NDVI_neg').reduceRegion({
    reducer: ee.Reducer.percentile([p_top_NDVI]), 
    geometry: image.geometry().bounds(), 
    scale: 30,
    bestEffort: true,
    maxPixels: 8e+12,
  });
  var n_perc_top_NDVI = ee.Number(d_perc_top_NDVI.get('NDVI_neg'));
  var i_top_NDVI = image.updateMask(image.select('NDVI_neg').lte(n_perc_top_NDVI));
  var d_perc_low_LST = i_top_NDVI.select('LST_NW').reduceRegion({
    reducer: ee.Reducer.percentile([p_coldest_Ts]), 
    geometry: image.geometry().bounds(), 
    scale: 30,
    maxPixels: 8e+12,
    });
  var n_perc_low_LST = ee.Number(d_perc_low_LST.get('LST_NW'));
  var i_cold_lst = i_top_NDVI.updateMask(i_top_NDVI.select('LST_NW').lte(n_perc_low_LST));
  var c_lst_cold20 =  i_cold_lst.updateMask(image.select('LST_NW').gte(200));
  //Creates a reducer that outputs the minimum value of its (first) input.
  var med_lst_cold20 = c_lst_cold20.select('LST_NW')
    .reduceRegion({reducer:  ee.Reducer.median(), geometry: image.geometry().bounds(), scale: 30, maxPixels: 1e9,});
   var n_med_lst_cold20 = ee.Number(med_lst_cold20.get('LST_NW'));  
    var sum_final_cold_pix = c_lst_cold20.select('int')
   .reduceRegion({reducer:  ee.Reducer.sum(), geometry: image.geometry().bounds(), scale: 30, maxPixels: 1e9,});
    var n_sum_final_cold_pix = ee.Number(sum_final_cold_pix.get('int'));  

   var dif_temp = c_lst_cold20.expression( 
      'abs(LST - LST_med_20cold)', {
        'LST' : c_lst_cold20.select('LST_NW'),
        'LST_med_20cold': n_med_lst_cold20,
      }).rename('dif_temp'); 

   c_lst_cold20 = c_lst_cold20.addBands(dif_temp);
   var d_red_min = c_lst_cold20.select('dif_temp', 'LST_NW', 'longitude','latitude', 'NDVI').reduceRegion({
    reducer: ee.Reducer.min(5), 
    geometry: image.geometry().bounds(), 
    scale: 30,
    maxPixels: 1e+09,
    });
    
  var n_Ts_cold = ee.Number(d_red_min.get('min1'));  
  var n_long_cold = ee.Number(d_red_min.get('min2'));
  var n_lat_cold = ee.Number(d_red_min.get('min3'));
  var n_ndvi_cold = ee.Number(d_red_min.get('min4'));
 
  // Make a Dictionary on the server.
  var d_cold_pixel = ee.Dictionary({
  temp: n_Ts_cold,
  ndvi: n_ndvi_cold,
  x: n_long_cold,
  y: n_lat_cold,
  sum: n_sum_final_cold_pix
  });
  
return d_cold_pixel;
};

 // --------------------------------------------- HOT PIXEL -----------------------------------------------------------
exports.fexp_hot_pixel = function(image,  p_lowest_NDVI, p_hottest_Ts) {

  // 1 - Identify the down xx% NDVI pixels ****************************************************************************************
   // Calculate percentile from  ndvi
   var d_perc_down_ndvi = image.select('pos_NDVI').reduceRegion({
      reducer: ee.Reducer.percentile([p_lowest_NDVI]), 
      geometry: image.geometry().bounds(), 
      scale: 30,
      maxPixels:10e14,
    });
    var n_perc_low_NDVI = ee.Number(d_perc_down_ndvi.get('pos_NDVI'));
 
  var i_low_NDVI = image.updateMask(image.select('pos_NDVI').lte(n_perc_low_NDVI));
 
  // 2 - Identify the hottest pixels  *********************************************************************************************
   var d_perc_top_lst = i_low_NDVI.select('LST_neg').reduceRegion({
      reducer: ee.Reducer.percentile([p_hottest_Ts]), 
      geometry: image.geometry().bounds(), 
      scale: 30,
      maxPixels: 10e14,
    });
  var n_perc_top_lst = ee.Number(d_perc_top_lst.get('LST_neg'));
  var i_top_LST = i_low_NDVI.updateMask(i_low_NDVI.select('LST_neg').lte(n_perc_top_lst));
  var c_lst_hotpix =  i_top_LST;
 
  var med_lst_hotpix = c_lst_hotpix.select('LST_NW')
    .reduceRegion({reducer: ee.Reducer.median(), geometry: image.geometry().bounds(), scale: 30, maxPixels: 10e14,});
  var n_med_lst_hotpix = ee.Number(med_lst_hotpix.get('LST_NW'));  //transforma de objeto para numero
 
  var sum_final_hot_pix = c_lst_hotpix.select('int')
    .reduceRegion({reducer:  ee.Reducer.sum(), geometry: image.geometry().bounds(), scale: 30, maxPixels: 10e14,});
   var n_sum_final_hot_pix = ee.Number(sum_final_hot_pix.get('int'));  //transforma de objeto para numero
  
   var dif_temp = c_lst_hotpix.expression( 
      'abs(LST - LST_med_hotpix)', {
        'LST' : c_lst_hotpix.select('LST_NW'),
        'LST_med_hotpix': ee.Number(n_med_lst_hotpix),
      }).rename('dif_temp'); 

  c_lst_hotpix = c_lst_hotpix.addBands(dif_temp);
  
    var d_min_diftemp_hot = c_lst_hotpix.select('dif_temp', 'LST_NW', 'Rn', 'G','SAVI','NDVI','longitude','latitude').reduceRegion({
    reducer: ee.Reducer.min(8), 
    geometry: image.geometry().bounds(), 
    scale: 30,
    maxPixels: 10e14,
    });
  
  var n_Ts_hot = d_min_diftemp_hot.get('min1');  
  var n_Rn_hot = d_min_diftemp_hot.get('min2');
  var n_G_hot = d_min_diftemp_hot.get('min3');
  var n_savi_hot = d_min_diftemp_hot.get('min4');
  var n_ndvi_hot = d_min_diftemp_hot.get('min5');
  var n_long_hot = d_min_diftemp_hot.get('min6');
  var n_lat_hot = d_min_diftemp_hot.get('min7');
  
  // Make a Dictionary on the server.
  var d_hot_pixel = ee.Dictionary({
    temp: n_Ts_hot,
    x: n_long_hot,
    y: n_lat_hot,
    Rn: n_Rn_hot,
    G: n_G_hot,
    ndvi: n_ndvi_hot,
    sum: n_sum_final_hot_pix,
  });

  return d_hot_pixel;
};

 // ------------------------------------------------------------------------------------------------
  // Sensible Heat Flux (H)
  exports.fexp_sensible_heat_flux = function(image, ux, UR, Rn24hobs, n_Ts_cold, d_hot_pixel) {


  var n_veg_hight = ee.Number(3);
  var n_zx = ee.Number(2);
  var n_hight = ee.Number(100); 
  var n_Cp = ee.Number(1004); //Air specific heat (J/kg/K)
  var n_K = ee.Number(0.41); // k is von Karman’s constant
  n_Ts_cold = ee.Number(n_Ts_cold);
  var n_Ts_hot = ee.Number(d_hot_pixel.get('temp'));
  var n_G_hot = ee.Number(d_hot_pixel.get('G'));
  var n_Rn_hot = ee.Number(d_hot_pixel.get('Rn'));
  var n_long_hot = ee.Number(d_hot_pixel.get('x'));
  var n_lat_hot = ee.Number(d_hot_pixel.get('y'));
  var p_hot_pix =  ee.Geometry.Point([n_long_hot, n_lat_hot]);
  
  //------------------------------------------------------------------------------------------------------------------------------------
  
  // savi index 
  var i_savi = image.select('SAVI');
  
  //Momentum roughness length (zom) at the weather station:
  var  n_zom = n_veg_hight.multiply(0.12);  //% at the weather station
  
  //Friction velocity at the weather station
    var i_ufric_ws = i_savi.expression( 
      '(n_K * ux)/ log(n_zx /n_zom)', {'n_K': n_K, 'n_zx': n_zx, 'n_zom': n_zom, 'ux': ux }); 

  // Wind speed at blending height at the weather station
  var i_u200 = i_savi.expression( 
      'i_ufric_ws *  (log(n_hight/n_zom)/n_K)', {'i_ufric_ws' : i_ufric_ws, 'n_hight' : n_hight, 'n_zom' : n_zom, 'n_K' : n_K});
      
  // Momentum roughness length (zom) for each pixel:
  var i_zom = i_savi.expression( 
      'exp((5.62 * (SAVI))-5.809)', {'SAVI' : i_savi,}); 
  
      
  //Friction velocity for each pixel:  
  var i_ufric = i_zom.expression( 
      '(n_K *u200) /(log(hight/i_zom))', {'u200' : i_u200,'hight': n_hight, 'i_zom':n_zom, 'n_K': n_K }).rename('u_fr');
  
  
  // Aerodynamic resistance to heat transport (rah)
  var z1= ee.Number(0.1);
  var z2= ee.Number(2);   
  
  
  var i_rah = i_ufric.expression( 
      '(log(z2/z1))/(i_ufric*0.41)', {'z2' : z2,'z1': z1, 'i_ufric':i_ufric }).rename('rah'); 
      
      
  var i_rah_first = i_rah.rename('rah_first');
  
    var n_ro_hot= (ee.Number(-0.0046).multiply(n_Ts_hot)).add(ee.Number(2.5538));

  //% #### Here begins the interactive process: ####
  
  // Sensible heat flux at the hot pixel
  var n_H_hot = ee.Number(n_Rn_hot).subtract(ee.Number(n_G_hot));
  
  var n= ee.Number(1);
  var n_dif= ee.Number(1);
  var n_dif_min = ee.Number(0.1);
  
  
  var list_dif = ee.List([]);
  var list_dT_hot = ee.List([]);
  var list_rah_hot = ee.List([]);
  var list_coef_a = ee.List([]);
  var list_coef_b = ee.List([]);
  
  
  // -------------------------------------------------------------------------------------------------------------

  
  for (n = 1; n < 16; n++) {
 
    // aerodynamic resistance to heat transport in hot pixel
    var d_rah_hot = i_rah.select('rah').reduceRegion({reducer: ee.Reducer.first(), geometry: p_hot_pix, scale: 30,maxPixels: 10e14});
    // Convert to Number for further use
    var n_rah_hot =   ee.Number(d_rah_hot.get('rah'));    
    
    // Near surface temperature difference in hot pixel (dT = Tz1 – Tz2)
    // dThot= Hhot*rah/(ρCp)
    var n_dT_hot = (n_H_hot.multiply(n_rah_hot)).divide(n_ro_hot.multiply(n_Cp));
    
    // Near surface temperature difference in cold pixel (dT = Tz1 – Tz2)
    var n_dT_cold = ee.Number(0);
  
    // Angular coefficient
    var n_coef_a = (n_dT_cold.subtract(n_dT_hot)).divide(n_Ts_cold.subtract(n_Ts_hot));
  
    // Linear coefficient
    var n_coef_b = n_dT_hot.subtract(n_coef_a.multiply(n_Ts_hot));
   
    // dT for each pixel
    var i_lst_med = image.select('T_LST_DEM');
    var i_dT_int = ee.Image(0).clip( image.geometry().bounds()).expression( 
      '(n_coef_a * i_lst_med) + n_coef_b', {'n_coef_a' : n_coef_a,'n_coef_b': n_coef_b, 'i_lst_med':i_lst_med }).rename('dT'); 
  
    // Air temperature (Ta) for each pixel (Ta = Ts-dT)
    var i_Ta = i_lst_med.expression( 
    'i_lst_med - i_dT_int', {'i_lst_med' : i_lst_med,'i_dT_int': i_dT_int});  //// 
    
    // ro (ρ) - air density (kg/m3)
    var i_ro = i_Ta.expression( 
    '(-0.0046 * i_Ta) + 2.5538', {'i_Ta' : i_Ta}).rename('ro');  // ro=-0.0046.*Ta+2.5538;
 
    // Sensible heat flux (H) for each pixel - iteration
    var i_H_int = i_dT_int.expression( 
      '(i_ro*n_Cp*i_dT_int)/i_rah', {'i_ro' : i_ro,'n_Cp': n_Cp, 'i_dT_int':i_dT_int, 'i_rah':i_rah }).rename('H');
      
      
    // Monin-Obukhov length (L) - iteration
    var i_L_int = i_dT_int.expression( 
      '-(i_ro*n_Cp*(i_ufric**3)*i_lst_med)/(0.41*9.81*i_H_int)',
      {'i_ro' : i_ro,'n_Cp': n_Cp, 'i_ufric':i_ufric, 'i_lst_med':i_lst_med, 'i_H_int':i_H_int }).rename('L');
  
    // Stability corrections for momentum and heat transport
    
    var img = ee.Image(0).clip( image.geometry().bounds());
    
     // stability corrections for stable conditions
    var i_psim_200 = img.expression( 
     '-5*(hight/i_L_int)', {'hight' : ee.Number(200),'i_L_int': i_L_int}).rename('psim_200');
    
    
    var i_psih_2 = img.expression( 
     '-5*(hight/i_L_int)',{'hight' : ee.Number(2),'i_L_int': i_L_int}).rename('psih_2');
  
   var i_psih_01 = img.expression( 
     '-5*(hight/i_L_int)',{'hight' : ee.Number(0.1),'i_L_int': i_L_int}).rename('psih_01');
  
  //unstable conditions
  
    // x for different height
    var i_x200 = i_L_int.expression( 
      '(1-(16*(hight/i_L_int)))**0.25',
    {'hight' : ee.Number(200),'i_L_int': i_L_int}).rename('i_x200');
    
    var i_x2 = i_L_int.expression( 
      '(1-(16*(hight/i_L_int)))**0.25',
    {'hight' : ee.Number(2),'i_L_int': i_L_int}).rename('i_x2');
    
    var i_x01 = i_L_int.expression( 
      '(1-(16*(hight/i_L_int)))**0.25',
    {'hight' : ee.Number(0.1),'i_L_int': i_L_int}); 
  
    // stability corrections for unstable conditions
   var i_psimu_200 = i_x200.expression( 
     '2*log((1+i_x200)/2)+log((1+i_x200**2)/2)-2*atan(i_x200)+0.5*pi',
    {'i_x200' : i_x200,'pi': ee.Number(3.14159265)});
    
    var i_psihu_2 = i_x2.expression( 
      '2*log((1+i_x2**2)/2)',
    {'i_x2' : i_x2});
  
    var i_psihu_01 = i_x01.expression( 
     '2*log((1+i_x01**2)/2)',
   {'i_x01' : i_x01});
  
  
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
  
    // Corrected value for the friction velocity (u_asterisk)
    // u_asterisk=(u200.*0.41)./(log(hight./zom_pixel)-psi_m200);
    i_ufric = i_ufric.expression( 
      '(u200*0.41)/(log(hight/i_zom)-i_psim_200)',
      {'u200' : i_u200,'hight': n_hight, 'i_zom':i_zom,'i_psim_200': i_psim_200});
      
    
    // Corrected value for the aerodinamic resistance to the heat transport (rah)
    //rah=(log(z2/z1)-psi_h2+psi_h01)./(u_asterisk*0.41);
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
  // insert each iteration value into a list
  list_dif = list_dif.add(n_dif);
  list_coef_a = list_coef_a.add(n_coef_a);
  list_coef_b = list_coef_b.add(n_coef_b);
  list_dT_hot = list_dT_hot.add(n_dT_hot);
  list_rah_hot = list_rah_hot.add(n_rah_hot);

}   // end

  var i_rah_final = i_rah.rename('rah');
  var i_dT_final = i_dT_int.rename('dT');
  var i_H_final = i_H_int.expression( 
    '(i_ro*n_Cp*i_dT_int)/i_rah', {'i_ro' : i_ro,'n_Cp': n_Cp, 'i_dT_int':i_dT_final, 'i_rah':i_rah_final }).rename('H');

  image = image.addBands([i_H_final, i_rah_final, i_dT_final, i_rah_first]);
  return image;
  };
exports.fexp_inst_et = function(image, Rn24hobs) {
  
  var i_Rn = image.select('Rn');
  var i_G = image.select('G');
  var i_lst = image.select('T_LST_DEM');
  var i_H_final = image.select('H');
   //Map.addLayer(i_H_final,{min:0,max: 300},'H') ;
  // Map.addLayer(i_G,{min:0,max: 100},'G')  
  // Map.addLayer(i_Rn,{min:0,max: 600},'Rn')  
  i_G=i_G.where(i_G.lt(0),0);
  i_H_final=i_H_final.where(i_H_final.lt(0),0);
  //var d_Rn_24h_med = Rn24hobs.reduceRegion({reducer: ee.Reducer.mean(),geometry: image.geometry().bounds(),scale: 20000,maxPixels: 900000000});
  //var n_Rn_24h_med = ee.Number(d_Rn_24h_med.get('Rn24h_G'));
  //var Rn24hobs =Rn24hobs.where(Rn24hobs, n_Rn_24h_med);
  // Latent heat flux - instantaneous value for the time of the satellite overpass (W/m2). (LE)
  var i_lambda_ET = i_H_final.expression( 
    '(i_Rn-i_G-i_H_fim)', {'i_Rn' : i_Rn,'i_G': i_G, 'i_H_fim':i_H_final }).rename('LE');
    
 // Map.addLayer(i_lambda_ET,{min:0,max: 700},'LE')  
  i_lambda_ET=i_lambda_ET.where(i_lambda_ET.lt(0),0);
  // Latent heat of vaporization or the heat absorbed when a kilogram of water evaporates - lambda (J/kg).
    var i_lambda = i_H_final.expression( 
    '(2.501-0.002361*(Ts-273.15))', {'Ts' : i_lst });
  
  // Instantaneous value of ET (mm/H)
  var i_ET_inst = i_H_final.expression( 
    '0.0036 * (i_lambda_ET/i_lambda)', {'i_lambda_ET' : i_lambda_ET, 'i_lambda' : i_lambda  }).rename('ET_inst');
  
 // Evaporative fraction 
   var i_FE = i_H_final.expression( 
    'i_lambda_ET/(i_Rn-i_G)', {'i_lambda_ET' : i_lambda_ET, 'i_Rn' : i_Rn, 'i_G' : i_G }).rename('FE');
   i_FE=i_FE.where(i_lambda_ET.lt(0),0);
 // FE = lambda_ET./(Rn-G);
  // ET24h (mm/day) - density (kg m-3)
  var i_ET24h_calc = i_H_final.expression( 
  '(86.4 *i_FE * Rn24hobs)/(i_lambda * dens)', {'i_FE' : i_FE, 'i_lambda' : i_lambda, 'Rn24hobs' : Rn24hobs, 'dens': ee.Number(1000) }).rename('ET_24h');
  //i_ET24h_calc=i_ET24h_calc.where(i_ET24h_calc.lt(0),0)
  image = image.addBands([i_ET_inst, i_ET24h_calc, i_lambda_ET,i_FE]);
  return image;

};  
  



