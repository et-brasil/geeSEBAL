//GEESEBAL - GOOGLE EARTH ENGINE APP FOR SURFACE ENERGY BALANCE ALGORITHM FOR LAND (SEBAL)
//CREATE BY: LEONARDO LAIPELT, RAFAEL KAYSER, ANDERSON RUHOFF AND AYAN FLEISCHMANN 
//PROJECT - ET BRASIL
//LAB - HIDROLOGIA DE GRANDE ESCALA [HGE] website: https://www.ufrgs.br/hge/author/hge/
//UNIVERSITY - UNIVERSIDADE FEDERAL DO RIO GRANDE DO SUL - UFRGS
//BRAZIL, RIO GRANDE DO SUL
//
//DOI
//VERSION 0.1
//CONTATC US: leonardo.laipelt@ufrgs.br

//EARTH ENGINE APP
//DEVELOPED BY: LEONARDO LAIPELT
//YEAR: 2020

//DEBUG
var debug = ui.url.get('debug', true);

//Null
var layersLists = null;
var NDVIlayer =null;
var currentImage = null;
var Area_Download =null;
var inspectorNDVI =null;
var imageSelect = null;


//REQUIERES
var Masks = require('users/leolaipelt/geeSEBAL:masks');
var Spectral_Indices = require('users/leolaipelt/geeSEBAL:tools');
var ET_24h = require('users/leolaipelt/geeSEBAL:image');

//BARS
function makeBARS(visParams) {
 var bar= {
     bbox: [0, 0, 1, 0.1], dimensions: '70x10',format: 'png',min: 0,max: 1,palette: visParams.palette,};
 return bar;}
 
//CREATE BARS AND LEGENDS
function BarImages(visParams,legendTitle){
  
  var colorBar = ui.Thumbnail({
  image: ee.Image.pixelLonLat().select(0),params: makeBARS(visParams),style: {stretch: 'horizontal', margin: '0px 4px', width: '200px', height: '24px'}}); 
  
  var legendLabels = ui.Panel({
  widgets: [ ui.Label(visParams.min),ui.Label(((visParams.max +visParams.min) / 2),{margin: '4px 4px', textAlign: 'center', stretch: 'horizontal'}),
    ui.Label(visParams.max, {margin: '4px 4px'})],layout: ui.Panel.Layout.flow('horizontal')});   

  var TitleLegend = ui.Label({
  value: legendTitle,style: {fontWeight: 'bold'}});
  
  var legendPanel = ui.Panel({ layout: ui.Panel.Layout.flow('vertical'), widgets:[TitleLegend, colorBar, legendLabels],style: {position: 'bottom-right', padding: '0px'}});
  
  return legendPanel;
}
//MIN and MAX
var minBox_NDVI =ui.Textbox({value: 0.0 , style:{color: '#000000', width: '97px'}});
var maxBox_NDVI =ui.Textbox({value: 1.0 , style:{color: '#000000', width: '97px'}});
var maxBox_FE =ui.Textbox({value: 1.0 , style:{color: '#000000', width: '97px'}});
var minBox_FE =ui.Textbox({value: 0.0, style:{color: '#000000', width: '97px'}});
var maxBox_LST =ui.Textbox({value: 320.0, style:{color: '#000000', width: '97px'}});
var minBox_LST =ui.Textbox({value: 285.0 , style:{color: '#000000', width: '97px'}});
var maxBox_ET =ui.Textbox({value: 7.0 , style:{color: '#000000', width: '97px'}});
var minBox_ET =ui.Textbox({value: 0.0 , style:{color: '#000000', width: '97px'}});

//PALETTES
//ET
var ET_Palette ={ min: 0, max: 7, 
 palette: ['deac9c', 'EDD9A6', 'f2dc8d', 'fff199', 'b5e684', '3BB369', '20998F', '25b1c1', '16678A', '114982', '0B2C7A']
 };
 
//NDVI
var NDVI_Palette ={ 
 min:0,  max: 1,palette: ['FFFFFF', 'CE7E45', 'DF923D', 'F1B555', 'FCD163', '99B718',
 '74A901', '66A000', '529400', '3E8601','023B01', '012E01', '011D01', '011301']};

//LST
var LST_Palette = { 
 min: 285,max: 320,
 palette: ['blue', 'limegreen', 'yellow', 'darkorange', 'red']};
var FE_Palette = { 
 min: 0,max: 1,palette: ["ffffff",'aec3d4']};
 
//COLLECTIONS
var LANDSAT_5=ee.ImageCollection("LANDSAT/LT05/C01/T1_SR").select([0,1,2,3,4,5,6,9],
["B","GR","R","NIR","SWIR_1","BRT","SWIR_2", "pixel_qa"]);
var LANDSAT_7=ee.ImageCollection("LANDSAT/LE07/C01/T1_SR").select([0,1,2,3,4,5,6,9],
["B","GR","R","NIR","SWIR_1","BRT","SWIR_2", "pixel_qa"]);
var LANDSAT_8=ee.ImageCollection("LANDSAT/LC08/C01/T1_SR").select([0,1,2,3,4,5,6,7,10],
["UB","B","GR","R","NIR","SWIR_1","SWIR_2","BRT","pixel_qa"]);

//LOGOS
var logo_eeSEBAL=ee.Image('projects/et-brasil/assets/eeSEBAL/geeSEBAL_logo_white').resample('bicubic');
var logo_ETBRASIL = ee.Image("projects/et-brasil/assets/LOGO_ET_BRASIL").resample('bicubic');
var logo_UFRGS = ee.Image("projects/et-brasil/assets/eeSEBAL/LOGOS_UFRGS_HGE_2").resample('bicubic');

//STYLE BUTTON
var styleButton ={color: '#ccc', fontWeight: 'bold', textAlign: 'center', border: '0.5px solid #dcdcdc', backgroundColor: '#f3f3f3', padding: '0px',
margin: 'auto', fontSize: '15px', width: '60px', height: '30px'};

//PRODUCTS
//RGB
var getRGB = function(){ 
 var RGB_image=currentImage.select(['R','GR','B']);
 return RGB_image;};
//ET
var getET24h = function(){ 
 var ET_calculation = ET_24h.ET_estimation(currentImage,parseFloat(topNDVI.getValue()),parseFloat(coldestTs.getValue()),
                                                 parseFloat(lowestNDVI.getValue()),parseFloat(hottestTs.getValue()));
 return ET_calculation ; };
//NDVI
var getndvi = function(){ 
 var ndvi = currentImage.normalizedDifference(['NIR', 'R']).rename('NDVI');
 return ndvi.select('NDVI') ; } ;
//LST
var getlst = function(){ 
  var lst = Spectral_Indices.spec_ind(currentImage);
  return lst.select('T_LST') ;
 };

//EF
var getFE = function(){ 
 var FE =ET_24h.ET_estimation(currentImage,parseFloat(topNDVI.getValue()),parseFloat(coldestTs.getValue()),
                                                 parseFloat(lowestNDVI.getValue()),parseFloat(hottestTs.getValue()));
 return FE.select('FE'); 
 } ;

//GET NAME IMAGE
var getName = function(){
  var time_start=currentImage.get('system:time_start');
  var date=ee.Date(time_start);
  var year=ee.Number(date.get('year'));
  var month=ee.Number(date.get('month'));
  var day=ee.Number(date.get('day'));
  var Name_image=ee.String(year).cat(ee.String('_')).cat(ee.String(month)).cat(ee.String('_')).cat(ee.String(day));
  return  Name_image;
 };
 
 //GET URL DOWNLOADS
 //STAND BY
var getURLNDVI= function(){ 
 var NDVI =getndvi();
 var Name=getName();
 var NameMap=Name.cat(ee.String('_NDVI'));
 var NameDownload=NameMap.getInfo();
 var downloadUrl = (NDVI.select(['NDVI']).multiply(10000).int()).divide(10000).getDownloadURL({ 
 name: NameDownload, scale: 100,geometry:NDVI.geometry().bounds(), maxPixels: 10e14 });
 return downloadUrl};

var getURLLST= function(){ 
 var lst = Spectral_Indices.spec_ind(currentImage);
 var Name=getName();
 var NameMap=Name.cat(ee.String('_LST'));
 var NameDownload=NameMap.getInfo();
 var downloadUrl = (lst.select(['T_LST']).multiply(10000).int()).divide(10000).getDownloadURL({ 
 name: NameDownload, scale: 100,geometry:lst.geometry().bounds(),maxPixels: 10e14 });
 return downloadUrl}  ;

var getURLRGB= function(){ 
 var RGB_image=(currentImage.select(['R','GR','B']).multiply(10000).int()).divide(10000);
 var Name=getName();
 var NameMap=Name.cat(ee.String('_RGB'));
 var NameDownload=NameMap.getInfo();
 var downloadUrl = RGB_image.getDownloadURL({ 
 name: NameDownload,geometry:RGB_image.geometry().bounds(), filePerBand: true, scale: 150,maxPixels: 10e14 });
 return downloadUrl};
 
var getURLFE= function(){ 
 var FE =ET_24h.ET_estimation(currentImage,NDVI_top_number,coldestTS_number,lowestNDVI_number,hottestTs_number);
 var Name=getName();
 var NameMap=Name.cat(ee.String('_EF'));
 var NameDownload=NameMap.getInfo();
 var downloadUrl = (FE.select(['FE']).multiply(10000).int()).divide(10000).getDownloadURL({ 
 name: NameDownload, scale: 100,maxPixels: 10e14,geometry:FE.geometry().bounds() });
 return downloadUrl}; 
 
var getURLET= function(){ 
 var ET_calculation = ET_24h.ET_estimation(currentImage,parseFloat(topNDVI.getValue()),parseFloat(coldestTs.getValue()),
                                                 parseFloat(lowestNDVI.getValue()),parseFloat(hottestTs.getValue()));
 var Name=getName();
 var NameMap=Name.cat(ee.String('ET_24h'));
 var NameDownload=NameMap.getInfo();
  var URL=ET_calculation.evaluate(function(result) {
   
   var downloadUrl = ((result.select(['ET_24h']).multiply(10000).int()).divide(10000)).getDownloadURL({ 
    name: NameDownload, scale: 100,geometry:result.geometry().bounds(),maxPixels: 10e18 });
    return downloadUrl;
    });
 
 return URL;

};

var panel = ui.Panel({
  layout: ui.Panel.Layout.flow('vertical'),
  style: {maxWidth: '430px'}
});

var map =ui.Map({
  center: {lon: -50, lat: -10, zoom: 4}, 
  });
map.style().set('cursor', 'crosshair');

var logo_EESEBAL=ui.Thumbnail({image:logo_eeSEBAL,params:{bands:['b1','b2','b3'],min:0,max:255},
style:{width:'50%', padding: '40px 10px 10px 10px', margin: 'auto'}});

var Title=ui.Label({ 
  value:'Version 0.1 - last updated: 18/09/2020',
    style: { fontSize: '14px', textAlign: 'center',  margin: 'auto', padding: '0px 30px'}
  });
 
var Subtitle=ui.Label({ 
  value:'GEE Surface Energy Balance Algorithm for Land (SEBAL)',
    style: { fontSize: '18px', textAlign: 'center',  margin: '10px 10px 10px 10px',}
  });
  
var WritingBy=ui.Label({ 
  value:'Written by Leonardo Laipelt, Rafael Kayser, Anderson Ruhoff, Institute of Hydraulic Research, Federal University of Rio Grande do Sul',
    style: { fontSize: '13px', textAlign: 'center',  margin: '10px 10px 10px 10px',}
  });
  
var CiteUs= ui.Label({value: 'If you use this tool, please cite: Laipelt et al. (2020), eeSEBAL: '+ 
    'A Google Earth Engine application for long term evapotranspiration estimation.',
  style: {fontSize: '12px', textAlign: 'center',padding: '5px 5px', margin: 'auto'}
});
var DOI=    ui.Label({value: 'doi:',
  style: {fontSize: '12px', textAlign: 'center',padding: '5px 5px', margin: 'auto'}
}) ;

var contact = ui.Label({
  value: 'Contact - leonardo.laipelt@ufrgs.br',
  style: {fontSize: '12px', textAlign: 'center',padding: '5px 5px', margin: 'auto'}
});
var ButtonFAQ = ui.Button({
  label: 'How to use this tool?',
   style: { fontSize: '26px', textAlign: 'top-center', color: '#000000',padding: '10px', margin: 'auto', width: '170px'},
   onClick: function(){map.add(FAQ_PANEL)}
});
var ButtonAboutSEBAL = ui.Button({
  label: 'About SEBAL',
  style: { fontSize: '26px', textAlign: 'top-center', color: '#000000',padding: '10px', margin: 'auto', width: '170px'},
  onClick: function(){map.add(ABOUT_SEBAL_Panel)}
  
});
var Panels_Button= ui.Panel({
   layout: ui.Panel.Layout.flow('horizontal'),
   widgets: [ButtonFAQ,ButtonAboutSEBAL]
});

var FAQ_PANEL=ui.Panel({
    layout: ui.Panel.Layout.flow('vertical'),
    widgets:[
    ui.Label({ value:'How to use this app?',
    style: { fontSize: '30px', fontWeight: 'bold', textAlign: 'left',padding: '10px',  margin: 'auto'}
  }), ui.Label({ value:'1 – Select data range to search for images from Landsat 5,7 and'+
  '8 datasets. It is important that data is writing as ‘YYYY-MM-DD’ (example: 2000-01-01).',
     style: {textAlign: 'left',padding: '10px',  margin: 'auto'}
  }),
  ui.Label({ value:'2 – Click on the map to select a location and then click in the search button to obtain a list of images.',
     style: {textAlign: 'left',padding: '10px',  margin: 'auto'}
  }),
  ui.Label({ value:'3 – Select a scene to show NDVI, LST, RGB and evapotranspiration products.',
     style: {textAlign: 'left',padding: '10px',  margin: 'auto'}
  }),
  ui.Label({ value:'4– Anchor Pixels percentiles: geeSEBAL has an automatic calibration as explained in our paper (Laipelt et al.(2020)) based on '+
  'Allen et al. (2013). Its possible to change those percentiles to obtained others results of Fraction Evaporative and daily Evapotranspiration.'+
  'Percentiles range between 0-100 (%) and are set as default 5%,20%,10% and 20% (Allen et al. 2013). Evapotranspiration estimation may be not possible'+
  'for some groups of percentiles due to not be possible to select hot and cold pixels with those percentiles.  See more in Allen et al. (2013), Laipelt. et al. (2020).',
     style: {textAlign: 'left',padding: '10px',  margin: 'auto'}
  }), 
  ui.Label({
    value: 'Allen, Richard G., Burnett, Boyd, Kramber, William, Huntington, Justin, Kjaersgaard, Jeppe, Kilic, Ayse, Kelly, Carlos, and Trezza, Ricardo, 2013.'+
    'Automated Calibration of the METRIC‐Landsat Evapotranspiration Process. Journal of the American Water Resources Association (JAWRA) . 49( 3): 563– 576',
     style: {textAlign: 'left',padding: '10px',  margin: 'auto',fontSize:'10px'},
     targetUrl: 'https://doi.org/10.1111/jawr.12056'
  }), ui.Label({
    value: 'Laipelt, L.; Ruhoff, A.L.; Fleischmann, A.S.; Kayser, R.H.B.; Kich, E.M.; da Rocha, H.R.; Neale, C.M.U. Assessment of an Automated Calibration of the SEBAL'+
    'Algorithm to Estimate Dry-Season Surface-Energy Partitioning in a Forest–Savanna Transition in Brazil. Remote Sens. 2020, 12, 1108.',
    style: {textAlign: 'left',padding: '10px',  margin: 'auto',fontSize:'10px'},
    targetUrl: 'https://doi.org/10.3390/rs12071108'
  })
        ],
    
    style: {
    position: 'top-center',
    shown: true,
    width: '40%',
    height: 'auto',
    padding: '5px',
    margin: '10px',
    }
  });

var ABOUT_SEBAL_Panel=ui.Panel({
    layout: ui.Panel.Layout.flow('vertical'),
    widgets:[
    ui.Label({ value:'Surface Energy Balance Algorithm for Land (SEBAL)',
    style: { fontSize: '30px', fontWeight: 'bold', textAlign: 'left',padding: '10px',  margin: 'auto'}
  }),
  ui.Label({ value:'Surface Energy Balance Algorithm for Land (SEBAL) was developed and validated by Bastiaanssen (Bastiaanssen, 1995; Bastiaanssen et al., 1998a, 1998b)'+
  'to estimate evapotranspiration (ET) from energy balance equation (Rn – G = LE + H), where LE, Rn, G and H are Latent Heat Flux, Net Radiation, Soil Heat Flux and Sensible Heat Flux, respectivelt.'+
  'SEBAL estimate LE as a residual of others energy fluxes (LE = Rn - LE - G).',
     style: {textAlign: 'left',padding: '10px',  margin: 'auto'}
  }),
   ui.Label({ value:'SEBAL algorithm has an internal calibration, assuming a linear relationship between dT and LST across domain area, where dT is designed as a vertical air temperature (Ta)'+
   ' floating over the land surface, considering two extreme conditions. At the hot and dry extreme condition,'+ 
   'LE is zero and H is equal to the available energy, whereas at the cold and wet extreme condition, H is zero and LE is equal to the available energy.',
     style: {textAlign: 'left',padding: '10px',  margin: 'auto'}
  }),
   ui.Label({ value:'Allen et al. (2013) suggested an automatic calibration procedure to represent extreme conditions to use in METRIC (adaptable for SEBAL),'+
   'using endmembers candidates from pre-defined percentiles of LST and NDVI. Allen et al. (2013) defined a subset of endmembers within the highest 5% of NDVI and the lowest 20% of LST to select cold extreme condition,'+
   'whereas endmembers within the lowest 10% NDVI and within the highest 20% of LST are used to select the hot extreme conditions.',
     style: {textAlign: 'left',padding: '10px',  margin: 'auto'}
  }),ui.Label({
    value: 'Bastiaanssen, W.G.M., 1995. Regionalization of surface flux densities and moisture indicators in composite terrain: a remote sensing approach under clear skies in Mediterranean climates. Dr. thesis,'+
    'Wageningen Agric. Univ. Wageningen Netherlands. SC-DLO, Wageningen. ',
     style: {textAlign: 'left',padding: '10px',  margin: 'auto',fontSize:'10px'},
     targetUrl: 'https://doi.org/90-5485-465-0'
  }),ui.Label({
    value: 'Bastiaanssen, W.G.M., Menenti, M., Feddes, R.A., Holtslag, A.A.M., 1998a. A remote sensing surface energy balance algorithm for land (SEBAL): 1. Formulation. J. Hydrol. 212–213, 198–212.' ,
     style: {textAlign: 'left',padding: '10px',  margin: 'auto',fontSize:'10px'},
     targetUrl: 'https://doi.org/10.1016/S0022-1694(98)00253-4'
  }),ui.Label({
    value: 'Bastiaanssen, W.G.M., Pelgrum, H., Wang, J., Ma, Y., Moreno, J.F., Roerink, G.J., van der Wal, T., 1998c. A remote sensing surface energy balance algorithm for land (SEBAL): 2. Validation. J. Hydrol. 212–213, 213–229.',
     style: {textAlign: 'left',padding: '10px',  margin: 'auto',fontSize:'10px'},
     targetUrl: 'https://doi.org/10.1016/S0022-1694(98)00254-6'
  }),
   ui.Label({
    value: 'Allen, Richard G., Burnett, Boyd, Kramber, William, Huntington, Justin, Kjaersgaard, Jeppe, Kilic, Ayse, Kelly, Carlos, and Trezza, Ricardo, 2013.'+
    'Automated Calibration of the METRIC‐Landsat Evapotranspiration Process. Journal of the American Water Resources Association (JAWRA) . 49( 3): 563– 576',
     style: {textAlign: 'left',padding: '10px',  margin: 'auto',fontSize:'10px'},
     targetUrl: 'https://doi.org/10.1111/jawr.12056'
  })
        ],
    style: { position: 'top-center',shown: true,width: '40%',height: 'auto',padding: '0%',margin: '0%',}});  
  
var CloseButton = ui.Button({
  label: 'Close',style: { fontSize: '26px', position: 'top-left', color: '#000000',padding: '5px', margin: '5px', width: '80px'},
  onClick: function(){
    map.remove(FAQ_PANEL);
  }
});
var CloseButtonSEBAL = ui.Button({
  label: 'Close',style: { fontSize: '26px', position: 'top-left', color: '#000000',padding: '5px', margin: '5px', width: '80px'},
  onClick: function(){
    map.remove(ABOUT_SEBAL_Panel);
  }
});

panel.add(logo_EESEBAL);
panel.add(Title);
panel.add(Subtitle);
panel.add(WritingBy);
//panel.add(CiteUs);
//panel.add(DOI);
panel.add(Panels_Button);
panel.add(contact);

//FAQ PANEL
FAQ_PANEL.add(CloseButton);
ABOUT_SEBAL_Panel.add(CloseButtonSEBAL);

var DataSearchPanel= ui.Panel({
  layout: ui.Panel.Layout.flow("horizontal"), 
  style: {position: 'top-center',}
});

//START DATE
var startDate = ui.Textbox({ 
 value: "2019-01-01", 
 placeholder:"yyyy-mm-dd",
style: { color: "#000000"}
});

//END DATE
var endDate = ui.Textbox({ 
 value: "2019-12-31", 
 placeholder:"yyyy-mm-dd" ,
 style: { color: "#000000"}});

var Box_Data= ui.Panel({
  layout: ui.Panel.Layout.flow('vertical'),
  style: {backgroundColor: "#367bf0", color: "#ffffff", margin: '2px 5px'},//border: '0.5px solid #000000' 
  widgets:
  
  [ui.Label({ value:'Data search',
    style: { fontSize: '18px', textAlign: 'center',backgroundColor: "#367bf0",  margin: '10px 10px 10px 10px'}
  }),
  
  ui.Panel({
    layout:ui.Panel.Layout.flow("horizontal"), 
    style: {position: 'top-center',margin: 'auto', width: '100%'},
    widgets: [ui.Panel({ layout: ui.Panel.Layout.flow("horizontal"), 
    style: {position: 'top-center', margin: 'auto'},
    widgets: [
      startDate,endDate]})]})]});

panel.add(Box_Data);

//var Meteorological_data = ui.Select({ 
// items: [ 
// { label:"Global Land Data Assimilation System - GLDAS", value: 0 }, 
// { label:"NCEP Climate Forecast System Version 2 - CFSV2", value: 1}, 

// ], 
// value: 0,
 //style:{color:'#000000'},
 //onChange: function(value) {
   
 //  flag_meteorological = value;
 //}, 
//});
/*
var Box_Meteorological= ui.Panel({
  layout: ui.Panel.Layout.flow('vertical'),
  style: {backgroundColor: "#367bf0", color: "#ffffff", margin: '2px 5px'},//border: '0.5px solid #000000' 
  widgets:
  
  [ui.Label({ value:'Meteorological data',
    style: { fontSize: '18px', textAlign: 'center',backgroundColor: "#367bf0",  margin: '10px 10px 10px 10px'}
  }),
  
  ui.Panel({
    layout:ui.Panel.Layout.flow("horizontal"), 
      style: {position: 'top-center',margin: 'auto', width: '100%'},
    widgets: [ui.Panel({ layout: ui.Panel.Layout.flow("horizontal"), 
    style: {position: 'top-center', margin: 'auto'},
    widgets: [
      Meteorological_data]})]})]});
*/

//BOX LOCATION
var Box_Location= ui.Panel({
  layout: ui.Panel.Layout.flow('vertical'),
  style: {backgroundColor: "#367bf0", color: "#ffffff", margin: '2px 5px'},//border: '0.5px solid #000000' 
  widgets:
  
  [ui.Label({ value:'Location Information',
    style: { fontSize: '18px', textAlign: 'center',backgroundColor: "#367bf0",  margin: '10px 10px 10px 10px'}
  }),
  
  ui.Panel({
    layout:ui.Panel.Layout.flow("vertical"), 
    style: {position: 'top-center'},
    widgets: [ui.Label({
      value:"Click on the map to select a location and then click in the search button to obtain a list of images.",
      style:{color: '#000000', margin: 'auto', padding: '10px'}})]
    
})] });

//ENDMEMBERS
var topNDVI = ui.Textbox({
  value: 5.0,
  onChange: function(value) {
    // set value with a dedicated method
    topNDVI.setValue(value);
  },
  style: { margin: 'auto', width: '80px', padding:'10px',color: "#000000"}
  });
  
var coldestTs = ui.Textbox({
  value: 20.0,
  onChange: function(value) {
    coldestTs.setValue(value);
  },
  style: { margin: 'auto', width: '80px', padding:'10px',color: "#000000"}});
  
var lowestNDVI = ui.Textbox({
  value: 10.0,
  onChange: function(value) {
    lowestNDVI.setValue(value);

  },
  style: { margin: 'auto', width: '80px', padding:'10px',color: "#000000"}});

var hottestTs = ui.Textbox({
  value: 20.0,
  onChange: function(value) {
    hottestTs.setValue(value);

  },
  style: { margin: 'auto', width: '80px', padding:'10px',color: "#000000"}});   
  
var AnchorPixelsPercentile = ui.Panel({
  layout: ui.Panel.Layout.flow('vertical'),
   style: {backgroundColor: "#367bf0", color: "#ffffff", margin: 'auto', width: '100%'},
  widgets: [ui.Label({ value:'Percentiles for automated calibration',
    style: { fontSize: '18px', textAlign: 'center',backgroundColor: "#367bf0",  margin: '10px 10px 10px 10px'}
  }),
  
  ui.Panel({
    layout:ui.Panel.Layout.flow("horizontal"), 
    style: {position: 'top-center', margin:'auto', width: '100%', padding:' 0px 30px'},
    widgets: [ ui.Panel({layout:ui.Panel.Layout.flow("vertical"), 
    style: {position: 'top-center'}, widgets:[ui.Label('Top NDVI',{'color':'#000000'}),topNDVI]}),
     ui.Panel({layout:ui.Panel.Layout.flow("vertical"), 
    style: {position: 'top-center'}, widgets:[ui.Label('Coldest LST',{'color':'#000000'}),coldestTs]}),
     ui.Panel({layout:ui.Panel.Layout.flow("vertical"), 
    style: {position: 'top-center'}, widgets:[ui.Label('Lowest NDVI',{'color':'#000000'}),lowestNDVI]}),
     ui.Panel({layout:ui.Panel.Layout.flow("vertical"), 
    style: {position: 'top-center'}, widgets:[ui.Label('Hottest LST',{'color':'#000000'}),hottestTs]}),
   ]})]
  
});

var Value_NDVI=ui.Label({
  value:'0',
  style:{textAlign: 'center', margin:'auto'}
  });
var Value_LST=ui.Label({
  value:'0',
  style:{textAlign: 'center', margin:'auto'}
  });
var Value_FE=ui.Label({
  value:'0',
  style:{textAlign: 'center', margin:'auto'}
  });

var Value_ET24h=ui.Label({
  value:'0',
  style:{textAlign: 'center', margin:'auto'}
  }) ; 
  
var PixelsValue_onclick =  ui.Panel({
    layout:ui.Panel.Layout.flow("horizontal"), 
    style: {position: 'bottom-center', margin:'auto', width: 'auto', padding:' 0px 0px'},
    widgets: [ ui.Panel({layout:ui.Panel.Layout.flow("vertical"), 
    style: {textAlign: 'center', margin:'auto'}, widgets:[ui.Label('NDVI',{'color':'#000000'}),Value_NDVI]}),
     ui.Panel({layout:ui.Panel.Layout.flow("vertical"), 
    style: {textAlign: 'center', margin:'auto'}, widgets:[ui.Label('EF',{'color':'#000000'}),Value_FE]}),
     ui.Panel({layout:ui.Panel.Layout.flow("vertical"), 
    style: {textAlign: 'center', margin:'auto'}, widgets:[ui.Label('LST',{'color':'#000000'}),Value_LST]}),
     ui.Panel({layout:ui.Panel.Layout.flow("vertical"), 
    style: {textAlign: 'center', margin:'auto'}, widgets:[ui.Label('ET 24h',{'color':'#000000'}),Value_ET24h]}),
   ]});

var lon_text =ui.Textbox({ 
 value: "Longitude", 
 style: { color: "#000000"}}); 

var lat_text =ui.Textbox({ 
 value: "Latitude", 
 style: { color: "#000000"}}); 

var latlon_search= ui.Button({
  label: 'Search',
   style: { fontSize: '26px', textAlign: 'top-center', color: '#000000',padding: '10px', margin: 'auto', width: '170px'}});

var LatLon_location=  ui.Panel({
    layout:ui.Panel.Layout.flow("vertical"), 
    style: {position: 'top-center', margin: 'auto'},
    widgets: [ui.Panel({ layout: ui.Panel.Layout.flow("horizontal"), 
    style: {position: 'top-center', margin: 'auto'},
    widgets: [lon_text,lat_text]}),
    latlon_search]});


var Image_Panel = ui.Panel({
  layout: ui.Panel.Layout.flow('vertical'),
  style: {backgroundColor: "#367bf0", color: "#ffffff", margin: '2px 5px'},//border: '0.5px solid #000000' 
  widgets:
  
  [ui.Label({ value:'Select scene',
    style: { fontSize: '18px', textAlign: 'center',backgroundColor: "#367bf0",  margin: '10px 10px 10px 10px'}
  })]});

//NDVI VALUES TEXT
var inspectorNDVI = ui.Panel([ui.Label('Click to get NDVI value')]);
var inspectorLST = ui.Panel([ui.Label('Click to get LST value')]);
var inspectorFE = ui.Panel([ui.Label('Click to get EF value')]);
var inspectorET = ui.Panel([ui.Label('Click to get Daily ET value')]);

 //map.add(inspectorNDVI);
var RGB_Button= ui.Button("RGB",
function(){ 
 var RGB =getRGB();
 
 var Name=getName();
 var NameMap=Name.cat(ee.String('_RGB'));
 map.addLayer(RGB, { 
 min: 0, 
 max: 3000,
 gamma: 1.4,
 bands: ['R','GR','B']
 },NameMap.getInfo()); 
 
 },false,{color: '#000000', width: '75%'});
var NDVI_Button= ui.Button("NDVI",
function(){ 
 map.remove(inspectorNDVI);
 map.remove(NDVIBARS);
 map.add(NDVIBARS);
 inspectorNDVI.widgets().set(0, ui.Label({
    value: 'Click to get NDVI value'}));
 var NDVI =(getndvi().multiply(100).int()).divide(100);
 
 var Name=getName();
 var NameMap=Name.cat(ee.String('_NDVI'));
 
 NDVIlayer = ui.Map.Layer( NDVI,NDVI_Palette ,NameMap.getInfo());

 map.layers().add(NDVIlayer); 
 
 minBox_NDVI.onChange(function(value) {
    map.remove(NDVIBARS); map.add(NDVIBARS); minBox_NDVI.setValue(value);
    NDVI_Palette.min =parseFloat(value);
    NDVIlayer.setVisParams(NDVI_Palette);
    NDVIBARS.widgets().reset();
    NDVIBARS.widgets().set(0,BarImages(NDVI_Palette,'NDVI'));

});

 maxBox_NDVI.onChange(function(value) {
    map.remove(NDVIBARS); map.add(NDVIBARS); maxBox_NDVI.setValue(value);
    NDVI_Palette.max =parseFloat(value);
    NDVIlayer.setVisParams(NDVI_Palette);
    NDVIBARS.widgets().reset();
    NDVIBARS.widgets().set(0,BarImages(NDVI_Palette,'NDVI'));

});
 
 map.add(inspectorNDVI);
 map.onClick(function(coords){
 
  var location = 'lon: ' + coords.lon.toFixed(4) + ' ' +
                 'lat: ' + coords.lat.toFixed(4);
  var click_point = ee.Geometry.Point(coords.lon, coords.lat);
  lon_text.setValue(coords.lon);
  lat_text.setValue(coords.lat);

  //inspectorNDVI.widgets().reset()
  inspectorNDVI.widgets().set(0, ui.Label({
    value: 'Loading...',
    style: {color: 'gray'}
  }));
  var point = ee.Geometry.Point(coords.lon, coords.lat);
  var sample = NDVI.sample(point, 30);
  var computedValue = sample.first().get('NDVI');
  
  computedValue.evaluate(function(result) {
    // When the server returns the value, show it.
    inspectorNDVI.widgets().set(0, ui.Label({
      value: 'NDVI: ' + result.toFixed(2),
      style: {width: '80px', textAlign: 'center'}
    }));
  
  });
 });
 },
 false,{color: '#000000', width: '75%'});

var LST_Button= ui.Button("LST",
function(){
 map.remove(inspectorLST);
 map.remove(LSTBARS);
 map.add(LSTBARS);
 
 inspectorLST.widgets().set(0, ui.Label({value: 'Click to get LST value'}));
 var LST =(getlst().multiply(100).int()).divide(100);

 var Name=getName();
 var NameMap=Name.cat(ee.String('_LST'));

 var LSTlayer = ui.Map.Layer( LST, LST_Palette,NameMap.getInfo());
 map.layers().add(LSTlayer); 
 
 maxBox_LST.onChange( function(value) {
    map.remove(LSTBARS); map.add(LSTBARS); maxBox_LST.setValue(value);
    LST_Palette.max =parseFloat(value);
    LSTlayer.setVisParams(LST_Palette);
    LSTBARS.widgets().reset();
    LSTBARS.widgets().set(0,BarImages(LST_Palette,'LST'));

});
 minBox_LST.onChange( function(value) {
    map.remove(LSTBARS); map.add(LSTBARS); minBox_LST.setValue(value);
    LST_Palette.min =parseFloat(value);
    LSTlayer.setVisParams(LST_Palette);
    LSTBARS.widgets().reset();
    LSTBARS.widgets().set(0,BarImages(LST_Palette,'LST'));

 });
 map.add(inspectorLST);
 map.onClick(function(coords){
  
  var location = 'lon: ' + coords.lon.toFixed(4) + ' ' +
                 'lat: ' + coords.lat.toFixed(4);
  var click_point = ee.Geometry.Point(coords.lon, coords.lat);
  lon_text.setValue(coords.lon);
  lat_text.setValue(coords.lat);
  
  //CLICK 
  //inspectorLST.widgets().reset()
  inspectorLST.widgets().set(0, ui.Label({
    value: 'Loading...',
    style: {color: 'gray'}
  }));
  var point = ee.Geometry.Point(coords.lon, coords.lat);
  var sample = LST.sample(point, 30);
  var computedValue = sample.first().get('T_LST');
  
  computedValue.evaluate(function(result) {
    inspectorLST.widgets().set(0, ui.Label({
      value: 'LST: ' + result.toFixed(2) ,
      style: {width: '80px', textAlign: 'center'}
    }) 
    );
  });
 });
 },
 false,{color: '#000000', width: '75%'});
 
var FE_Button= ui.Button("EF",
function(){
 map.remove(inspectorFE);
 map.remove(FEBARS);
 map.add(FEBARS);
 inspectorFE.widgets().set(0, ui.Label({value: 'Click to get EF value'})); 
 var FE =(getFE().multiply(100).int()).divide(100);

 var Name=getName();
 var NameMap=Name.cat(ee.String('_EF'));
 var FElayer = ui.Map.Layer(FE,FE_Palette ,NameMap.getInfo());
 map.layers().add(FElayer);
 
 maxBox_FE.onChange( function(value) {
    map.remove(FEBARS); map.add(FEBARS); maxBox_FE.setValue(value);
    FE_Palette.max =parseFloat(value);
    FElayer.setVisParams(FE_Palette);
    FEBARS.widgets().reset();
    FEBARS.widgets().set(0,BarImages(FE_Palette,'FE'));

});
 minBox_FE.onChange( function(value) {
    map.remove(FEBARS); map.add(FEBARS); minBox_FE.setValue(value);
    FE_Palette.min =parseFloat(value);
    FElayer.setVisParams(FE_Palette);
    FEBARS.widgets().reset();
    FEBARS.widgets().set(0,BarImages(FE_Palette,'EF'));
 });

 map.add(inspectorFE);
 map.onClick(function(coords){
  var location = 'lon: ' + coords.lon.toFixed(4) + ' ' +
                 'lat: ' + coords.lat.toFixed(4);
  var click_point = ee.Geometry.Point(coords.lon, coords.lat);
  lon_text.setValue(coords.lon);
  lat_text.setValue(coords.lat);
  
  //CLICK 
  inspectorFE.widgets().reset();
  inspectorFE.widgets().set(0, ui.Label({
    value: 'Loading...',
    style: {color: 'gray'}
  }));
  var point = ee.Geometry.Point(coords.lon, coords.lat);
  var sample = FE.sample(point, 30);
  var computedValue = sample.first().get('FE');
  
  computedValue.evaluate(function(result) {
    // When the server returns the value, show it.
    inspectorFE.widgets().set(0, ui.Label({
      value: 'EF: ' + result.toFixed(2),
      style: {width: '80px', textAlign: 'center'}
    }));
  });
 
 });
},
false,{color: '#000000', width: '75%'});
  
var ET_Button= ui.Button("Daily ET",
function(){
 map.remove(inspectorET);
 map.remove(ETBARS);
 map.add(ETBARS);
 inspectorET.widgets().set(0, ui.Label({value: 'Click to get daily ET value'})); 
 
 var et24h =(getET24h().select('ET_24h').multiply(100).int()).divide(100);
 var LE =(getET24h().select('LE').multiply(100).int()).divide(100);
 var G =(getET24h().select('G').multiply(100).int()).divide(100);
 var H =(getET24h().select('H').multiply(100).int()).divide(100);
 var Rn =(getET24h().select('Rn').multiply(100).int()).divide(100);

 var Name=getName();
 var NameMap=Name.cat(ee.String('_ET'));

 var ETlayer = ui.Map.Layer(et24h,ET_Palette,NameMap.getInfo());

 map.layers().add(ETlayer);

 maxBox_ET.onChange( function(value) {
    map.remove(ETBARS); map.add(ETBARS); maxBox_ET.setValue(value);
    ET_Palette.max =parseFloat(value);
    ETlayer.setVisParams(ET_Palette);
    ETBARS.widgets().reset();
    ETBARS.widgets().set(0,BarImages(ET_Palette,'Daily ET mm day-1'));

});
 minBox_ET.onChange( function(value) {
    map.remove(ETBARS); map.add(ETBARS); minBox_ET.setValue(value);
    ET_Palette.min =parseFloat(value);
    ETlayer.setVisParams(ET_Palette);
    ETBARS.widgets().reset();
    ETBARS.widgets().set(0,BarImages(ET_Palette,'Daily ET mm day-1'));

});
 map.add(inspectorET);
 map.onClick(function(coords){
  
  var location = 'lon: ' + coords.lon.toFixed(4) + ' ' +
                 'lat: ' + coords.lat.toFixed(4);
  var click_point = ee.Geometry.Point(coords.lon, coords.lat);
  lon_text.setValue(coords.lon);
  lat_text.setValue(coords.lat);
  //CLICK 
  inspectorET.widgets().reset();
  inspectorET.widgets().set(0, ui.Label({
    value: 'Loading...',
    style: {color: 'gray'}
  }));
  var point = ee.Geometry.Point(coords.lon, coords.lat);
  var sample = et24h.sample(point, 30);
  var computedValue = sample.first().get('ET_24h');
  
  computedValue.evaluate(function(result) {
    
    inspectorET.widgets().set(0, ui.Label({
      value: 'ET: ' + result.toFixed(2),
      style: {width: '80px', textAlign: 'center'}
    }));
  });
 });
 },
 false,{color: '#000000', width: '75%'});

var urlLabel = ui.Label('⇓',{ fontSize: '18px', textAlign: 'center',backgroundColor: "#ffffff",
border: '0.5px solid #000000' ,color: "#000000", padding: '5px', width: '50px', height: '30px',  margin: '10px', shown: false});

//DOWNLOAD BUTTONS
//STAND BY
var DownloadNDVI_Button = ui.Label({value:"⇓", style: styleButton});
var CheckBoxNDVI = ui.Checkbox('', false);
var DownloadRGB_Button = ui.Label({value:"⇓", style: styleButton});
var CheckBoxRGB = ui.Checkbox('', false);  
var DownloadLST_Button = ui.Label({value:"⇓", style: styleButton});  
var CheckBoxLST = ui.Checkbox('', false); 
var DownloadFE_Button = ui.Label({value:"⇓", style: styleButton}); 
var CheckBoxFE = ui.Checkbox('', false);  
var DownloadET_Button = ui.Label({value:"⇓", style: styleButton});  
var CheckBoxET = ui.Checkbox('', false);

//CheckBoxRGB.onChange(function() {DownloadRGB_Button.setUrl(getURLRGB())})
//CheckBoxNDVI.onChange(function() {DownloadNDVI_Button.setUrl(getURLNDVI())})
//CheckBoxLST.onChange(function() {DownloadLST_Button.setUrl(getURLLST())})
//CheckBoxFE.onChange(function() {DownloadFE_Button.setUrl(getURLFE())})
//CheckBoxET.onChange(function() {DownloadET_Button.setUrl(getURLET())})

var Products_Panel = ui.Panel({
  layout: ui.Panel.Layout.flow('vertical'),
  style: {backgroundColor: "#367bf0", color: "#ffffff", margin: '2px 5px'},//border: '0.5px solid #000000' 
  widgets:
  [ui.Label({ value:'Products',
    style: { fontSize: '18px', textAlign: 'center',backgroundColor: "#367bf0",  margin: '10px 10px 10px 10px'}
  }),
   ui.Panel({
      layout:ui.Panel.Layout.flow("vertical"), 
      style: {position: 'top-center'},
      widgets: [
        ui.Panel({
          layout: ui.Panel.Layout.flow("horizontal"),
                style: {position: 'top-center'},
                widgets: [RGB_Button, DownloadRGB_Button]
        })]})
        ,ui.Panel({
      layout:ui.Panel.Layout.flow("vertical"), 
      style: {position: 'top-center'},
      widgets: [
        ui.Panel({
          layout: ui.Panel.Layout.flow("horizontal"),
                style: {position: 'top-center'},
                widgets: [NDVI_Button, DownloadNDVI_Button]
        }),
          ui.Panel({
          layout: ui.Panel.Layout.flow("horizontal"),
                style: {position: 'top-center'},
                widgets: [ui.Label('min:',{color: '#000000'}),minBox_NDVI,ui.Label('max:',{color: '#000000'}), maxBox_NDVI]
          
        })]}),  ui.Panel({
      layout:ui.Panel.Layout.flow("vertical"), 
      style: {position: 'top-center'},
      widgets: [
        ui.Panel({
          layout: ui.Panel.Layout.flow("horizontal"),
                style: {position: 'top-center'},
                widgets: [LST_Button, DownloadLST_Button]
          
        }),
          ui.Panel({
          layout: ui.Panel.Layout.flow("horizontal"),
                style: {position: 'top-center'},
                widgets: [ui.Label('min:',{color: '#000000'}),minBox_LST,ui.Label('max:',{color: '#000000'}), maxBox_LST]
          
        })]}),
          ui.Panel({
      layout:ui.Panel.Layout.flow("vertical"), 
      style: {position: 'top-center'},
      widgets: [
        ui.Panel({
          layout: ui.Panel.Layout.flow("horizontal"),
                style: {position: 'top-center'},
                widgets: [FE_Button, DownloadFE_Button]
        }),
          ui.Panel({
          layout: ui.Panel.Layout.flow("horizontal"),
                style: {position: 'top-center'},
                widgets: [ui.Label('min:',{color: '#000000'}),minBox_FE,ui.Label('max:',{color: '#000000'}), maxBox_FE]
        })]}),
          ui.Panel({
      layout:ui.Panel.Layout.flow("vertical"), 
      style: {position: 'top-center'},
      widgets: [
        ui.Panel({
          layout: ui.Panel.Layout.flow("horizontal"),
                style: {position: 'top-center'},
                widgets: [ET_Button, DownloadET_Button]
        }),
          ui.Panel({
          layout: ui.Panel.Layout.flow("horizontal"),
                style: {position: 'top-center'},
                widgets: [ui.Label('min:',{color: '#000000'}),minBox_ET,ui.Label('max:',{color: '#000000'}), maxBox_ET]
          })]}),
]});

function buildImageSelect(items){ 
 var imageSelect = ui.Select({ 
 items: items, 
 placeholder: "List of scenes", 
 style: { 
 width: '95%', 
 margin: '10px 10px 10px 10px',
 position: 'top-center'}});
 
 return imageSelect; 
}

//IMAGES 
map.onClick(function(coords){
  //map.layers().reset();
  var location = 'lon: ' + coords.lon.toFixed(4) + ' ' +
                 'lat: ' + coords.lat.toFixed(4);
  var click_point = ee.Geometry.Point(coords.lon, coords.lat);
  //map.layers().set(1, ui.Map.Layer(click_point, {color: '#000000'},'Point'));
  lon_text.setValue(coords.lon);
  lat_text.setValue(coords.lat)});

latlon_search.onClick(function() {
  var click_point=ee.Geometry.Point([lon_text.getValue(),lat_text.getValue()]);
  map.layers().reset();
   map.layers().set(1, ui.Map.Layer(click_point, {color: '#000000'},'Point'));
  map.setCenter(lon_text.getValue(),lat_text.getValue(),7);
  
  panel.remove(imageSelect);
  panel.remove(Image_Panel);
  panel.remove(Products_Panel);
  panel.remove(AnchorPixelsPercentile);
  
  map.remove(ETBARS);
  map.remove(LSTBARS);
  map.remove(NDVIBARS);
  map.remove(FEBARS);

  var IMAGES_AVAILABLES=LANDSAT_5.filterBounds(click_point).filterDate(startDate.getValue(), endDate.getValue())
  .map(Masks.f_cloudMaskL457_SR).map(Masks.f_albedoL5L7).merge(
    LANDSAT_7.filterBounds(click_point).filterDate(startDate.getValue(), endDate.getValue())
  .map(Masks.f_cloudMaskL457_SR).map(Masks.f_albedoL5L7)).merge(
    LANDSAT_8.filterBounds(click_point).filterDate(startDate.getValue(), endDate.getValue())
  .map(Masks.f_cloudMaskL8_SR).map(Masks.f_albedoL8));
  
    IMAGES_AVAILABLES.sort("system:time_start").evaluate(function(col){ 
    var items = []; 
    col.features.forEach(function(feature){ 
    
    var label = "LANDSAT ID"+ feature.properties.LANDSAT_ID + " / " 
     + "Cloud " + feature.properties.CLOUD_COVER + "% / " ;
    var value = feature.properties;
     
     items.push({label: label, value: value});}); 
  
  imageSelect = buildImageSelect(items);
  imageSelect.onChange(function(properties){
  
  panel.remove(Products_Panel);
  panel.remove(AnchorPixelsPercentile);
  
  map.remove(NDVIBARS);
  map.remove(ETBARS);
  map.remove(FEBARS); 
  map.remove(LSTBARS);
  map.remove(inspectorNDVI);
  map.remove(inspectorLST); 
  map.remove(inspectorFE);
  map.remove(inspectorET);
  
  currentImage = ee.Image(IMAGES_AVAILABLES
  
 .filterMetadata("LANDSAT_ID", "equals", properties.LANDSAT_ID) 
 .first() 
 );
 //DownloadNDVI_Button.setUrl(getURLNDVI())
// DownloadRGB_Button.setUrl(getURLRGB())
// DownloadLST_Button.setUrl(getURLLST())
 //DownloadFE_Button.setUrl(getURLFE())
// DownloadET_Button.setUrl(getURLET())
 
 panel.add(Products_Panel);
 panel.add(AnchorPixelsPercentile);
 
 map.layers().reset();
 map.remove(inspectorNDVI);
 map.remove(inspectorLST); 
 map.remove(inspectorFE);
 map.remove(inspectorET);
 }); 
  
  //ADD PANEL
  panel.add(Image_Panel);
  panel.add(imageSelect);
  
  //map.add(ETBARS);
  //map.add(LSTBARS);
 // map.add(NDVIBARS);
  //map.add(FEBARS);
  
  //map.add(getValuePixels());
  ///IMAGES
});

});
//function getValuePixels(){
//   map.add(PixelsValue_onclick);
//   var ET_calculation = ET_24h.ET_estimation(currentImage,parseFloat(topNDVI.getValue()),parseFloat(coldestTs.getValue()),
//                                                 parseFloat(lowestNDVI.getValue()),parseFloat(hottestTs.getValue()));
//  
//  var max_value = ET_calculation.select(['NDVI','T_LST','FE','ET_24h']).reduceRegion({
//  reducer: ee.Reducer.max(),
//  geometry: ee.Geometry.Point([lon_text.getValue(),lat_text.getValue()]),
//  scale : 30});
//  print(max_value);
//}

//BARS
var ETBARS=BarImages(ET_Palette,'Daily ET mm day-1');
var LSTBARS=BarImages(LST_Palette,'Land Surface Temperature K');
var NDVIBARS=BarImages(NDVI_Palette,'NDVI');
var FEBARS=BarImages(FE_Palette,'Evaporative Fraction');

//PANEL
//panel.add(Box_Meteorological);
panel.add(Box_Location);
panel.add(LatLon_location);

//LOGOS PANEL
var Logos_PANEL=ui.Panel({
    style: {
    width: '205px',
    height: 'auto',
    padding: '10px',
    position: 'bottom-right'
    },
    widgets:[ui.Thumbnail({image:logo_ETBRASIL,params:{bands:['b1','b2','b3'],min:0,max:255},style:{width:'190px',height:'auto', margin: 'auto'}}),
    ui.Thumbnail({image:logo_UFRGS,params:{bands:['b1','b2','b3'],min:0,max:255},style:{width:'190px',height:'auto', margin: 'auto'}})
    ]
  });


ui.root.clear();
ui.root.add(panel);
map.add(Logos_PANEL);

ui.root.add(map);
