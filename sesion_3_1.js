/*******************************************************************************
 * APAIT · Proyecto GEE · Sesión 3 · Demo 3.1
 * CÁLCULO DE SUPERFICIES: TOTAL, POR CLASE Y POR UNIDAD ESPACIAL
 *
 * Ejemplo: incendio de la Sierra de la Culebra (Sarracín–Ferreras, Zamora),
 *          15–24 de junio de 2022.
 * La superficie se calcula sumando el área real de cada píxel (pixelArea),
 * nunca contando píxeles, y siempre con una escala explícita.
 ******************************************************************************/

// ===== 0. ZONA DE ESTUDIO Y PARÁMETROS ======================================
// Polígono que engloba el área quemada (de Mahíde/Villardeciervos al valle del Tera)
var aoi = ee.Geometry.Polygon([[
  [-6.45, 41.78], [-5.95, 41.78], [-5.95, 42.03], [-6.45, 42.03], [-6.45, 41.78]
]]);
var PRE    = ['2022-05-15', '2022-06-14'];
var POST   = ['2022-06-25', '2022-07-15'];   // antes del incendio de Losacio (17 jul), al este
var ESCALA = 20;                             // m
var MIN_PIXELES = 25;                        // 25 píxeles de 20 m = 1 ha

Map.centerObject(aoi, 10);
Map.addLayer(ee.Image().paint(aoi, 0, 2), {palette: ['yellow']}, '0. AOI');

// ===== 1. RESULTADO DE LA SESIÓN 2 ==========================================
// Composites → dNBR → umbral de Otsu → limpieza → máscaras auxiliares.
// Es la misma cadena de los scripts S2_3 y S2_4a, condensada.
var S2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED');
var CS = ee.ImageCollection('GOOGLE/CLOUD_SCORE_PLUS/V1/S2_HARMONIZED');
var BANDAS = ['B2', 'B3', 'B4', 'B8', 'B11', 'B12'];

function mascaraS2(img) {
  return img.select(BANDAS).divide(10000)
    .updateMask(img.select('cs_cdf').gte(0.60))
    .set('system:time_start', img.get('system:time_start'));
}
function compositeS2(fechas) {
  return S2.filterBounds(aoi).filterDate(fechas[0], fechas[1])
    .linkCollection(CS, ['cs_cdf']).map(mascaraS2).median().clip(aoi);
}
function otsu(histograma) {
  histograma = ee.Dictionary(histograma);
  var cuentas = ee.Array(histograma.get('histogram'));
  var medias  = ee.Array(histograma.get('bucketMeans'));
  var n       = medias.length().get([0]);
  var total   = cuentas.reduce(ee.Reducer.sum(), [0]).get([0]);
  var suma    = medias.multiply(cuentas).reduce(ee.Reducer.sum(), [0]).get([0]);
  var media   = suma.divide(total);
  var varEntre = ee.List.sequence(1, n.subtract(1)).map(function (i) {
    var cA = cuentas.slice(0, 0, i);
    var nA = cA.reduce(ee.Reducer.sum(), [0]).get([0]);
    var mA = medias.slice(0, 0, i).multiply(cA)
      .reduce(ee.Reducer.sum(), [0]).get([0]).divide(nA);
    var nB = total.subtract(nA);
    var mB = suma.subtract(nA.multiply(mA)).divide(nB);
    return nA.multiply(mA.subtract(media).pow(2))
      .add(nB.multiply(mB.subtract(media).pow(2)));
  });
  return medias.slice(0, 0, n.subtract(1)).sort(varEntre).get([-1]);
}

var pre  = compositeS2(PRE);
var post = compositeS2(POST);
var dNBR = pre.normalizedDifference(['B8', 'B12'])
  .subtract(post.normalizedDifference(['B8', 'B12'])).rename('dNBR');
var umbral = ee.Number(otsu(dNBR.reduceRegion({
  reducer: ee.Reducer.histogram(255, 0.01),
  geometry: aoi, scale: ESCALA, maxPixels: 1e10, bestEffort: true
}).get('dNBR')));

var moda = dNBR.gt(umbral).unmask(0).clip(aoi)
  .focalMode({radius: 1, kernelType: 'square', units: 'pixels'});
var tam = moda.selfMask().connectedPixelCount({maxSize: 100, eightConnected: true});
var limpia = moda.updateMask(tam.gte(MIN_PIXELES)).unmask(0);

var aguaPerm = ee.Image('JRC/GSW1_4/GlobalSurfaceWater').select('occurrence').gt(80).unmask(0);
var wc = ee.ImageCollection('ESA/WorldCover/v200').first().select('Map').clip(aoi);
var combustible = wc.eq(10).or(wc.eq(20)).or(wc.eq(30)).or(wc.eq(40));

// Máscara final 0/1 (1 = afectado) y clases de severidad USGS
var afectado = limpia.and(aguaPerm.not()).and(combustible).rename('afectado').clip(aoi);
var severidad = ee.Image(0)
  .where(dNBR.gte(0.10), 1).where(dNBR.gte(0.27), 2)
  .where(dNBR.gte(0.44), 3).where(dNBR.gte(0.66), 4)
  .updateMask(afectado).unmask(0).rename('severidad').toByte().clip(aoi);

var PALETA_SEV  = ['ffffff', 'fee08b', 'fdae61', 'f46d43', 'a50026'];
var NOMBRES_SEV = ['Sin afectar', 'Baja', 'Moderada-baja', 'Moderada-alta', 'Alta'];
print('1. Umbral de Otsu (dNBR):', umbral);
Map.addLayer(post, {bands: ['B12', 'B8', 'B4'], min: 0, max: 0.4}, '1a. Composite POST (SWIR)', false);
Map.addLayer(dNBR, {min: -0.2, max: 0.8,
  palette: ['1a9850', 'ffffbf', 'fdae61', 'd73027', '7b3294']}, '1b. dNBR', false);
Map.addLayer(afectado.selfMask(), {palette: ['a50026']}, '1c. Afectado (resultado de la sesión 2)');

// ===== 2. QUÉ SE SUMA: EL ÁREA DE CADA PÍXEL ================================
print('2. Área media de un píxel a ' + ESCALA + ' m (m²):', ee.Image.pixelArea().reduceRegion({
  reducer: ee.Reducer.mean(), geometry: aoi, scale: ESCALA, maxPixels: 1e10
}).get('area'));

// Imagen de hectáreas por píxel: se multiplica por cualquier máscara 0/1
var HA = ee.Image.pixelArea().divide(1e4).rename('ha');

// ===== 3. SUPERFICIE TOTAL ==================================================
// tileScale: 4 reparte el cálculo en teselas más pequeñas (evita errores de memoria)
var haTotal = afectado.multiply(HA).rename('ha').reduceRegion({
  reducer: ee.Reducer.sum(), geometry: aoi, scale: ESCALA, maxPixels: 1e10, tileScale: 4
}).get('ha');
print('3. Superficie afectada (ha):', haTotal);

// La escala importa: la misma máscara medida a otras escalas
[20, 60, 200].forEach(function (s) {
  print('3. Superficie a ' + s + ' m (ha):', afectado.multiply(HA).rename('ha').reduceRegion({
    reducer: ee.Reducer.sum(), geometry: aoi, scale: s, maxPixels: 1e10, tileScale: 4
  }).get('ha'));
});

// ===== 4. FUNCIÓN DE SUPERFICIE POR CLASE (reductor agrupado) ================
// Devuelve una FeatureCollection {clase, nombre, ha} lista para tabla y gráfico
function areaPorClase(imgClases, banda, nombres) {
  var grupos = HA.addBands(imgClases.select(banda)).reduceRegion({
    reducer: ee.Reducer.sum().group({groupField: 1, groupName: 'clase'}),
    geometry: aoi, scale: ESCALA, maxPixels: 1e10, tileScale: 4
  }).get('groups');
  return ee.FeatureCollection(ee.List(grupos).map(function (g) {
    g = ee.Dictionary(g);
    var clase = ee.Number(g.get('clase'));
    return ee.Feature(null, {
      clase: clase,
      nombre: ee.Dictionary(nombres).get(clase.format('%d'), clase.format('%d')),
      ha: ee.Number(g.get('sum')).round()
    });
  }));
}

// ===== 5. POR CLASE DE SEVERIDAD ============================================
Map.addLayer(severidad.selfMask(), {min: 0, max: 4, palette: PALETA_SEV}, '5. Severidad');
var nombresSev = {'1': 'Baja', '2': 'Moderada-baja', '3': 'Moderada-alta', '4': 'Alta'};
var haSev = areaPorClase(severidad.selfMask(), 'severidad', nombresSev);
print('5. Hectáreas por severidad:', haSev);
print(ui.Chart.feature.byFeature({features: haSev, xProperty: 'nombre', yProperties: ['ha']})
  .setChartType('ColumnChart')
  .setOptions({title: '5. Hectáreas por clase de severidad', legend: {position: 'none'},
               colors: ['f46d43'], vAxis: {title: 'ha'}}));

// ===== 6. POR COBERTURA DEL SUELO (antes del incendio) =======================
var nombresWC = {
  '10': 'Árboles', '20': 'Matorral', '30': 'Pastizal', '40': 'Cultivos', '50': 'Urbano',
  '60': 'Suelo desnudo', '70': 'Nieve', '80': 'Agua', '90': 'Humedal', '95': 'Manglar', '100': 'Musgo'
};
var wcAfectado = wc.updateMask(afectado);
Map.addLayer(wcAfectado, {}, '6. Cobertura del área afectada (WorldCover)');
var haWC = areaPorClase(wcAfectado, 'Map', nombresWC);
print('6. Hectáreas por cobertura:', haWC);
print(ui.Chart.feature.byFeature({features: haWC, xProperty: 'nombre', yProperties: ['ha']})
  .setChartType('PieChart')
  .setOptions({title: '6. Superficie afectada por cobertura del suelo'}));

// ===== 7. POR UNIDAD ADMINISTRATIVA (FAO GAUL) ===============================
// En España, GAUL nivel 2 = provincias. Para municipios, subid como asset los
// límites del IGN (Assets → New → Shape files) y usadlos en lugar de GAUL.
var unidades = ee.FeatureCollection('FAO/GAUL/2015/level2').filterBounds(aoi);
var porUnidad = afectado.multiply(HA).rename('ha').reduceRegions({
  collection: unidades, reducer: ee.Reducer.sum(), scale: ESCALA, tileScale: 4
});
print('7. Hectáreas por provincia:', porUnidad.select(['ADM2_NAME', 'sum']));
Map.addLayer(ee.Image().paint(unidades, 0, 2), {palette: ['000000']}, '7. Límites provinciales', false);

// ===== 8. POR CELDAS DE UNA MALLA (distribución espacial) ===================
// Malla de 2 km en ETRS89 / UTM 30N; cada celda guarda sus hectáreas afectadas
var malla = aoi.coveringGrid(ee.Projection('EPSG:25830'), 2000);
var porCelda = afectado.multiply(HA).rename('ha').reduceRegions({
  collection: malla, reducer: ee.Reducer.sum(), scale: ESCALA, tileScale: 4
});
print('8. Nº de celdas de 2 km:', malla.size());
// Coropletas: cada celda pintada con su valor (máximo posible: 400 ha)
Map.addLayer(ee.Image().float().paint(porCelda, 'sum').clip(aoi),
  {min: 0, max: 400, palette: ['fff5f0', 'fcbba1', 'fb6a4a', 'cb181d', '67000d']},
  '8a. Hectáreas afectadas por celda de 2 km');
Map.addLayer(ee.Image().paint(malla, 0, 1), {palette: ['808080']}, '8b. Malla de 2 km', false);

// Histograma: ¿cuántas celdas están poco, medio o muy afectadas?
print(ui.Chart.feature.histogram({features: porCelda.filter(ee.Filter.gt('sum', 0)),
  property: 'sum', maxBuckets: 20})
  .setOptions({title: '8. Celdas según hectáreas afectadas', hAxis: {title: 'ha por celda'},
               legend: {position: 'none'}}));

// ===== 9. EXPORTAR LAS TABLAS ===============================================
Export.table.toDrive({collection: haSev, description: 'ha_por_severidad',
  folder: 'APAIT', fileFormat: 'CSV', selectors: ['nombre', 'ha']});
Export.table.toDrive({collection: haWC, description: 'ha_por_cobertura',
  folder: 'APAIT', fileFormat: 'CSV', selectors: ['nombre', 'ha']});
