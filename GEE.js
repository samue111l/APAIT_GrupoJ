// APAIT · Proyecto GEE · Grupo J
// Incendio del sur de Ávila (valle del Alberche), julio–agosto de 2026

// ===== 0. PARÁMETROS ========================================================
var aoi = ee.Geometry.Rectangle([-4.85, 40.20, -4.15, 40.55]);
var PRE  = ['2026-07-01', '2026-07-16'];   // fin exclusivo: antes del primer foco
var POST = ['2026-08-16', '2026-09-01'];   // tras el control oficial (15 ago)
var ESCALA = 20;                           // m
var MIN_PIXELES = 25;                      // 25 píxeles de 20 m = 1 ha
var UMBRAL_BIBLIO = 0.10;                  // USGS

Map.centerObject(aoi, 10);
Map.addLayer(ee.Image().paint(aoi, 0, 2), {palette: ['yellow']}, 'AOI');

// ===== 1. COMPOSITES PRE Y POST (Sentinel-2 L2A + Cloud Score+) ==============
var S2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED');
var CS = ee.ImageCollection('GOOGLE/CLOUD_SCORE_PLUS/V1/S2_HARMONIZED');

function preparar(img) {
  return img.select(['B8', 'B12']).divide(10000)
    .updateMask(img.select('cs_cdf').gte(0.60));
}
function composite(fechas) {
  return S2.filterBounds(aoi).filterDate(fechas[0], fechas[1])
    .linkCollection(CS, ['cs_cdf']).map(preparar).median().clip(aoi);
}
// Nº de imágenes y % de nubes
[PRE, POST].forEach(function (f) {
  var col = S2.filterBounds(aoi).filterDate(f[0], f[1]);
  print('Imágenes ' + f[0] + ' a ' + f[1], col.size());
  print('% de nubes', col.aggregate_array('CLOUDY_PIXEL_PERCENTAGE'));
});

var pre  = composite(PRE);
var post = composite(POST);
Map.addLayer(post, {bands: ['B12', 'B8', 'B8'], min: 0, max: 0.4}, 'POST (SWIR)', false);

// ===== 2. dNBR ==============================================================
var dNBR = pre.normalizedDifference(['B8', 'B12'])
  .subtract(post.normalizedDifference(['B8', 'B12'])).rename('dNBR');
Map.addLayer(dNBR, {min: -0.2, max: 0.8,
  palette: ['1a9850', 'ffffbf', 'fdae61', 'd73027', '7b3294']}, 'dNBR');
