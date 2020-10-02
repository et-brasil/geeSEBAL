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
