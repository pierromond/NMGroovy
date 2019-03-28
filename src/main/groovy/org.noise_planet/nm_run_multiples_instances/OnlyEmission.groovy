package org.noise_planet.nm_run_multiples_instances

import com.esotericsoftware.kryo.Kryo
import com.esotericsoftware.kryo.serializers.MapSerializer
import groovy.sql.Sql

/** author : Aumond Pierre

 Ce code permet de faire une unique de simulation

 Le résultat se trouve dans un dossier zippé
 Il s'agit d'une version Postgre/Postgis (lecture/ecriture plus lent qu'avec h2gis mais jointures beaucoup plus rapides)

 **/

// Importation des librairies nécessaire au code
import groovy.time.TimeCategory
import groovy.transform.SourceURI
import org.apache.commons.io.FileUtils
import org.cts.CRSFactory
import org.cts.crs.CoordinateReferenceSystem
import org.cts.crs.GeodeticCRS
import org.cts.op.CoordinateOperationFactory
import org.cts.registry.EPSGRegistry
import org.cts.registry.RegistryManager
import org.h2gis.api.EmptyProgressVisitor
import org.h2gis.functions.io.csv.CSVDriverFunction
import org.h2gis.functions.spatial.crs.ST_Transform
import org.h2gis.functions.spatial.volume.GeometryExtrude
import org.h2gis.utilities.JDBCUtilities
import org.h2gis.utilities.SFSUtilities
import org.h2gis.utilities.SpatialResultSet
import org.h2gis.utilities.TableLocation
import org.h2gis.utilities.wrapper.ConnectionWrapper
import org.locationtech.jts.algorithm.Angle
import org.locationtech.jts.geom.*
import org.locationtech.jts.geom.util.GeometryEditor
import org.orbisgis.noisemap.core.*

import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.Paths
import java.sql.Connection
import java.sql.DriverManager
import java.sql.PreparedStatement
import java.sql.SQLException
import java.util.zip.ZipEntry
import java.util.zip.ZipOutputStream

class OnlyEmission {
    static void main(String[] args) {
        OnlyEmission oneRun = new OnlyEmission()
        oneRun.run()
    }

    void run() {
/**
 ///////////////////////////////////////////
 // Paramètres d'entrée et initialisations //
 ///////////////////////////////////////////

 */
        String workspace_output = "D:\\aumond\\Documents\\Recherche_hors_projet\\2019_04_PasseportRecherche\\Traitement"
        // le workspace output doit respecter une arborescence specifique
        //output
        //|- data
        //|- config

        String input_file = "D:\\aumond\\Documents\\Recherche_hors_projet\\2019_04_PasseportRecherche\\Donnees\\Eleves6.csv"

        @SourceURI
        URI sourceUri
        Path scriptLocation = Paths.get(sourceUri)
        String rootPath = scriptLocation.getParent().getParent().getParent().getParent().toString()

  /**
 ///////////////////////////////////////////
 // FIN Paramètres d'entrée et initialisations //
 ///////////////////////////////////////////

 */


        try {
            PrintWriter pw = new PrintWriter(new File("D:\\sortieModelFlow2.csv"))
            System.out.println("Init")
            CSVDriverFunction csv = new CSVDriverFunction()

            CSVDriverFunction exp = new CSVDriverFunction()
// init du vecteur de fréquences

            /*VL;ML;PL;Vit;Vitest;NbVoies;Dist;Haute*/

            int nvar = 8

            ArrayList<Double> VL = new ArrayList<Double>()
            ArrayList<Double> ML = new ArrayList<Double>()
            ArrayList<Double> PL = new ArrayList<Double>()
            ArrayList<Double> Vit = new ArrayList<Double>()
            ArrayList<Double> Vitest = new ArrayList<Double>()
            ArrayList<Double> NbVoies = new ArrayList<Double>()
            ArrayList<Double> Dist = new ArrayList<Double>()
            ArrayList<Double> Haute = new ArrayList<Double>()
            ArrayList<Double> Temps = new ArrayList<Double>()
            ArrayList<Double> Vlh = new ArrayList<Double>()
            ArrayList<Double> MlH = new ArrayList<Double>()
            ArrayList<Double> Plh = new ArrayList<Double>()
            ArrayList<Double> PK = new ArrayList<Double>()
            ArrayList<Double> Result = new ArrayList<Double>()
            def list = [63, 125, 250, 500, 1000, 2000, 4000, 8000]


            new File(input_file).splitEachLine(";") { fields ->
                def res = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
                def kk = 0
                VL.add(fields[0].toFloat() )
                ML.add(fields[1].toFloat() )
                PL.add(fields[2].toFloat() )
                Vit.add(fields[3].toFloat() )
                Vitest.add(fields[4].toFloat() )
                NbVoies.add(fields[5].toFloat() )
                Dist.add(fields[6].toFloat() )
                Haute.add(fields[7].toFloat())
                Temps.add(fields[8].toFloat())
                Vlh.add(fields[9].toFloat())
                MlH.add(fields[10].toFloat())
                Plh.add(fields[11].toFloat())
                PK.add(fields[12].toFloat())

                for (f in list) {
                    //RSParametersCnossos srcParameters = new RSParametersCnossos(fields[4].toFloat(), fields[4].toFloat(), fields[4].toFloat(), 0, 0,
                    //        fields[0].toFloat(), fields[1].toFloat(), fields[2].toFloat(), 0, 0, f, 20, "NL01", 0, 0, 250, 1)
                    //RSParametersCnossos srcParameters = new RSParametersCnossos(fields[4].toFloat(), fields[4].toFloat(), fields[4].toFloat(), 0, 0,
                    //        3600, 0, 0, 0, 0, f, 20, "NL01", 0, 0, 250, 1)
                    RSParametersCnossos srcParameters = new RSParametersCnossos(fields[4].toFloat(), fields[4].toFloat(), fields[4].toFloat(), 0, 0,
                            fields[9].toFloat(), fields[10].toFloat(), fields[11].toFloat(), 0, 0, f, 20, "NL01", 0, 0, 250, 1)


                    srcParameters.setSlopePercentage(RSParametersCnossos.computeSlope(0, 0, 10))
                    res[kk] = EvaluateRoadSourceCnossos.evaluate(srcParameters)
                    res[kk] = 10 * Math.log10(Math.pow(10, EvaluateRoadSourceCnossos.evaluate(srcParameters) / 10) + Math.pow(10, res[kk] / 10))
                    kk++

                }

                res[8]=wToDba(DbaToW(res[0]-26.2)+DbaToW(res[1]-16.1)+DbaToW(res[2]-8.6)+DbaToW(res[3]-3.2)+DbaToW(res[4])+DbaToW(res[5]+1.2)+DbaToW(res[6]+1)+DbaToW(res[7]-1.1))


                Result.add(res)
                pw.write(fields[12].toFloat() + ";" + res[8].toString() + "\n")

            }



            pw.close()


int a=1
        } finally {

        }
    }
/**
 * Combine l'ensemble des rayons en prenant en compte favorable defavorable
 * @param aGlobal
 * @param propData
 * @param propDataOut
 * @param p
 */
    private void computeWithMeteo(HashMap<Integer, double[]> aGlobal, PropagationProcessPathData propData, ComputeRaysOut propDataOut, double[] p) {

        EvaluateAttenuationCnossos evaluateAttenuationCnossos = new EvaluateAttenuationCnossos()
        double[] aGlobalMeteo
        double[] aGlobalPathFav
        double[] aGlobalPathHom

        for (PropagationPath propath : propDataOut.propagationPaths) {
            double angleRad = Angle.angle(propath.pointList.get(0).coordinate, propath.pointList.get(propath.pointList.size() - 1).coordinate)
            double rad2rose = (-angleRad + Math.PI / 2)
            int roseindex = Math.round(rad2rose / (2 * Math.PI / p.length))

            propath.setFavorable(false)
            evaluateAttenuationCnossos.evaluate(propath, propData)
            aGlobalPathFav = evaluateAttenuationCnossos.getaGlobal()
            // todo bug meteo
            propath.setFavorable(true)
            evaluateAttenuationCnossos.evaluate(propath, propData)
            aGlobalPathHom = evaluateAttenuationCnossos.getaGlobal()

            aGlobalMeteo = sumArrayWithPonderation(aGlobalPathFav, aGlobalPathHom, p[roseindex])

            if (aGlobal.containsKey(propath.idSource)) {
                aGlobalMeteo = sumArray(aGlobalMeteo, aGlobal.get(propath.idSource))
                aGlobal.replace(propath.idSource, aGlobalMeteo)
            } else {
                aGlobal.put(propath.idSource, aGlobalMeteo)
            }


        }
    }

/**
 * Permet de passer de dBA en pression
 * @param dBA
 * @return
 */
    private static double DbaToW(double dBA) {
        return Math.pow(10.0, dBA / 10.0)
    }
/**
 * Passe de pression en dBA
 * @param w
 * @return
 */
    private static double wToDba(double w) {
        return 10 * Math.log10(w)
    }

/**
 * Fais la somme de deux arrays
 * @param array1
 * @param array2
 * @return
 */
    private double[] sumArray(double[] array1, double[] array2) {
        double[] sum = new double[array1.length]
        for (int i = 0; i < array1.length; i++) {
            sum[i] = wToDba(DbaToW(array1[i]) + DbaToW(array2[i]))
        }
        return sum
    }

/**
 * Fais la moyenne ponderee par p de deux arrays
 * @param array1
 * @param array2
 * @param p
 * @return
 */
    private double[] sumArrayWithPonderation(double[] array1, double[] array2, double p) {
        double[] sum = new double[array1.length]
        for (int i = 0; i < array1.length; i++) {
            sum[i] = wToDba(p * DbaToW(array1[i]) + (1 - p) * DbaToW(array2[i]))
        }
        return sum
    }

/**
 * Permet d'exploser les multipolygones en poylgones
 * @param geom
 * @param polygon
 */
    private static void addGeometry(List<Geometry> geom, Geometry polygon) {
        if (polygon instanceof Polygon) {
            geom.add((Polygon) polygon)
        } else {
            for (int i = 0; i < polygon.getNumGeometries(); i++) {
                addGeometry(geom, polygon.getGeometryN(i))
            }
        }

    }

    /**
     * Permet de tracer les rayons dans un vtk lisible dans paraview
     * @param filename
     * @param propDataOut
     * @throws IOException
     */
    private void writeKML_Line(String filename, ComputeRaysOut propDataOut, FastObstructionTest mesh) throws IOException {

        FileWriter fileWriter = new FileWriter(filename)
        fileWriter.write('<?xml version="1.0" encoding="UTF-8"?>\n')
        fileWriter.write('<kml xmlns="http://www.opengis.net/kml/2.2">\n')
        fileWriter.write('<Document>\n')
        fileWriter.write('\t<name>Paths</name>\n')
        fileWriter.write('<Style id="yellowLineGreenPoly">\n')
        fileWriter.write('<LineStyle>\n')
        fileWriter.write('<color>7f00ffff</color>\n')
        fileWriter.write('<width>1</width>\n')
        fileWriter.write('</LineStyle>\n')
        fileWriter.write('<PolyStyle>\n')
        fileWriter.write('<color>7f00ff00</color>\n')
        fileWriter.write('</PolyStyle>\n')
        fileWriter.write('</Style>\n')

        GeometryFactory gf = new GeometryFactory()
        CRSFactory cRSFactory = new CRSFactory()

        RegistryManager registryManager = cRSFactory.getRegistryManager()
        registryManager.addRegistry(new EPSGRegistry())

        CoordinateReferenceSystem inputCRS = cRSFactory.getCRS("EPSG:2154")
        CoordinateReferenceSystem targetCRS = cRSFactory.getCRS("EPSG:4326")
        List<GeometryEditor.CoordinateOperation> op = CoordinateOperationFactory.createCoordinateOperations((GeodeticCRS) inputCRS, (GeodeticCRS) targetCRS)

        for (int j = 0; j < propDataOut.propagationPaths.size(); j++) {
            fileWriter.write('\t<Placemark>\n')
            fileWriter.write('\t\t<name>Absolute Extruded</name>\n')
            fileWriter.write('\t\t<styleUrl>#yellowLineGreenPoly</styleUrl>\n')
            fileWriter.write('\t\t\t<LineString>\n')
            fileWriter.write('\t\t\t<altitudeMode>relativeToGround</altitudeMode>\n')
            fileWriter.write('\t\t\t<coordinates>\n')
            for (PropagationPath.PointPath p : propDataOut.propagationPaths.get(j).getPointList()) {
                //Geometry pointLine = GeometryFactory.createPointFromInternalCoord(p.coordinate)
                double z = mesh.getHeightAtPosition(p.coordinate)

                Geometry outPutGeom = (Point) gf.createPoint(p.coordinate).copy()
                outPutGeom.geometryChanged()
                outPutGeom.apply(new ST_Transform.CRSTransformFilter(op.get(0)))
                outPutGeom.setSRID(4326)
                Coordinate coordinate = outPutGeom.getCoordinate()


                fileWriter.write("\t\t\t" + coordinate.x.toString() + "," + coordinate.y.toString() + "," + (p.coordinate.z - z).toString() + "\n")

            }
            fileWriter.write('\t\t\t</coordinates>\n')
            fileWriter.write('\t\t\t</LineString>\n')
            fileWriter.write('\t</Placemark>\n')
        }
        fileWriter.write('\t</Document>\n')
        fileWriter.write('</kml>')
        fileWriter.close()
    }
/**
 * Permet de tracer les rayons dans un vtk lisible dans paraview
 * @param filename
 * @param propDataOut
 * @throws IOException
 */
    private void writeKML_Result(String filename, FastObstructionTest mesh, List<Coordinate> receivers, HashMap<Integer, Double> Result_by_receivers) throws IOException {

        FileWriter fileWriter = new FileWriter(filename)
        fileWriter.write('<?xml version="1.0" encoding="UTF-8"?>\n')
        fileWriter.write('<kml xmlns="http://www.opengis.net/kml/2.2" xmlns:gx="http://www.google.com/kml/ext/2.2" xmlns:kml="http://www.opengis.net/kml/2.2" xmlns:atom="http://www.w3.org/2005/Atom">\n')
        fileWriter.write('<Document>\n')
        fileWriter.write('<name>Results</name>\n')


        fileWriter.write('<StyleMap id="msn_shaded_dot1">\n')
        fileWriter.write('	<Pair>\n')
        fileWriter.write('	<key>normal</key>\n')
        fileWriter.write('  <styleUrl>#sn_shaded_dot1</styleUrl>\n')
        fileWriter.write('	</Pair>\n')
        fileWriter.write('   <Pair>\n')
        fileWriter.write('   <key>highlight</key>\n')
        fileWriter.write('		<styleUrl>#sh_shaded_dot1</styleUrl>\n')
        fileWriter.write('    </Pair>\n')
        fileWriter.write('</StyleMap>\n')
        fileWriter.write('     <Style id="sn_shaded_dot1">\n')
        fileWriter.write('      <IconStyle>\n')
        fileWriter.write('     <color>6414F000</color>\n')
        fileWriter.write('		<scale>1.2</scale>\n')
        fileWriter.write('     <Icon>\n')
        fileWriter.write('     <href>http://maps.google.com/mapfiles/kml/shapes/shaded_dot.png</href>\n')
        fileWriter.write('    </Icon>\n')
        fileWriter.write('	</IconStyle>\n')
        fileWriter.write('<LabelStyle>\n')
        fileWriter.write('<color>0FFFFFF</color>\n')
        fileWriter.write('</LabelStyle>\n')
        fileWriter.write('<BalloonStyle>\n')
        fileWriter.write('</BalloonStyle>\n')
        fileWriter.write('<ListStyle>\n')
        fileWriter.write('</ListStyle>\n')
        fileWriter.write('</Style>\n')
        fileWriter.write('<Style id="sh_shaded_dot1">\n')
        fileWriter.write('<IconStyle>\n')
        fileWriter.write('	<color>6414F000</color>\n')
        fileWriter.write(' <scale>1.4</scale>\n')
        fileWriter.write('	<Icon>\n')
        fileWriter.write('		<href>http://maps.google.com/mapfiles/kml/shapes/shaded_dot.png</href>\n')
        fileWriter.write('	</Icon>\n')
        fileWriter.write(' </IconStyle>\n')
        fileWriter.write('<LabelStyle>\n')
        fileWriter.write('	<color>0FFFFFF</color>\n')
        fileWriter.write('</LabelStyle>\n')
        fileWriter.write('<BalloonStyle>\n')
        fileWriter.write('</BalloonStyle>\n')
        fileWriter.write('<ListStyle>\n')
        fileWriter.write('   </ListStyle>\n')
        fileWriter.write('</Style>\n')

        fileWriter.write('<StyleMap id="msn_shaded_dot2">\n')
        fileWriter.write('	<Pair>\n')
        fileWriter.write('	<key>normal</key>\n')
        fileWriter.write('  <styleUrl>#sn_shaded_dot2</styleUrl>\n')
        fileWriter.write('	</Pair>\n')
        fileWriter.write('   <Pair>\n')
        fileWriter.write('   <key>highlight</key>\n')
        fileWriter.write('		<styleUrl>#sh_shaded_dot2</styleUrl>\n')
        fileWriter.write('    </Pair>\n')
        fileWriter.write('</StyleMap>\n')
        fileWriter.write('     <Style id="sn_shaded_dot2">\n')
        fileWriter.write('      <IconStyle>\n')
        fileWriter.write('     <color>6414F078</color>\n')
        fileWriter.write('		<scale>1.2</scale>\n')
        fileWriter.write('     <Icon>\n')
        fileWriter.write('     <href>http://maps.google.com/mapfiles/kml/shapes/shaded_dot.png</href>\n')
        fileWriter.write('    </Icon>\n')
        fileWriter.write('	</IconStyle>\n')
        fileWriter.write('<LabelStyle>\n')
        fileWriter.write('<color>0FFFFFF</color>\n')
        fileWriter.write('</LabelStyle>\n')
        fileWriter.write('<BalloonStyle>\n')
        fileWriter.write('</BalloonStyle>\n')
        fileWriter.write('<ListStyle>\n')
        fileWriter.write('</ListStyle>\n')
        fileWriter.write('</Style>\n')
        fileWriter.write('<Style id="sh_shaded_dot2">\n')
        fileWriter.write('<IconStyle>\n')
        fileWriter.write('	<color>6414F078</color>\n')
        fileWriter.write(' <scale>1.4</scale>\n')
        fileWriter.write('	<Icon>\n')
        fileWriter.write('		<href>http://maps.google.com/mapfiles/kml/shapes/shaded_dot.png</href>\n')
        fileWriter.write('	</Icon>\n')
        fileWriter.write(' </IconStyle>\n')
        fileWriter.write('<LabelStyle>\n')
        fileWriter.write('	<color>0FFFFFF</color>\n')
        fileWriter.write('</LabelStyle>\n')
        fileWriter.write('<BalloonStyle>\n')
        fileWriter.write('</BalloonStyle>\n')
        fileWriter.write('<ListStyle>\n')
        fileWriter.write('   </ListStyle>\n')
        fileWriter.write('</Style>\n')

        fileWriter.write('<StyleMap id="msn_shaded_dot3">\n')
        fileWriter.write('	<Pair>\n')
        fileWriter.write('	<key>normal</key>\n')
        fileWriter.write('  <styleUrl>#sn_shaded_dot3</styleUrl>\n')
        fileWriter.write('	</Pair>\n')
        fileWriter.write('   <Pair>\n')
        fileWriter.write('   <key>highlight</key>\n')
        fileWriter.write('		<styleUrl>#sh_shaded_dot3</styleUrl>\n')
        fileWriter.write('    </Pair>\n')
        fileWriter.write('</StyleMap>\n')
        fileWriter.write('     <Style id="sn_shaded_dot3">\n')
        fileWriter.write('      <IconStyle>\n')
        fileWriter.write('     <color>6414F0FF</color>\n')
        fileWriter.write('		<scale>1.2</scale>\n')
        fileWriter.write('     <Icon>\n')
        fileWriter.write('     <href>http://maps.google.com/mapfiles/kml/shapes/shaded_dot.png</href>\n')
        fileWriter.write('    </Icon>\n')
        fileWriter.write('	</IconStyle>\n')
        fileWriter.write('<LabelStyle>\n')
        fileWriter.write('<color>0FFFFFF</color>\n')
        fileWriter.write('</LabelStyle>\n')
        fileWriter.write('<BalloonStyle>\n')
        fileWriter.write('</BalloonStyle>\n')
        fileWriter.write('<ListStyle>\n')
        fileWriter.write('</ListStyle>\n')
        fileWriter.write('</Style>\n')
        fileWriter.write('<Style id="sh_shaded_dot3">\n')
        fileWriter.write('<IconStyle>\n')
        fileWriter.write('	<color>6414F0FF</color>\n')
        fileWriter.write(' <scale>1.4</scale>\n')
        fileWriter.write('	<Icon>\n')
        fileWriter.write('		<href>http://maps.google.com/mapfiles/kml/shapes/shaded_dot.png</href>\n')
        fileWriter.write('	</Icon>\n')
        fileWriter.write(' </IconStyle>\n')
        fileWriter.write('<LabelStyle>\n')
        fileWriter.write('	<color>0FFFFFF</color>\n')
        fileWriter.write('</LabelStyle>\n')
        fileWriter.write('<BalloonStyle>\n')
        fileWriter.write('</BalloonStyle>\n')
        fileWriter.write('<ListStyle>\n')
        fileWriter.write('   </ListStyle>\n')
        fileWriter.write('</Style>\n')

        fileWriter.write('<StyleMap id="msn_shaded_dot4">\n')
        fileWriter.write('	<Pair>\n')
        fileWriter.write('	<key>normal</key>\n')
        fileWriter.write('  <styleUrl>#sn_shaded_dot4</styleUrl>\n')
        fileWriter.write('	</Pair>\n')
        fileWriter.write('   <Pair>\n')
        fileWriter.write('   <key>highlight</key>\n')
        fileWriter.write('		<styleUrl>#sh_shaded_dot4</styleUrl>\n')
        fileWriter.write('    </Pair>\n')
        fileWriter.write('</StyleMap>\n')
        fileWriter.write('     <Style id="sn_shaded_dot4">\n')
        fileWriter.write('      <IconStyle>\n')
        fileWriter.write('     <color>6414B4FF</color>\n')
        fileWriter.write('		<scale>1.2</scale>\n')
        fileWriter.write('     <Icon>\n')
        fileWriter.write('     <href>http://maps.google.com/mapfiles/kml/shapes/shaded_dot.png</href>\n')
        fileWriter.write('    </Icon>\n')
        fileWriter.write('	</IconStyle>\n')
        fileWriter.write('<LabelStyle>\n')
        fileWriter.write('<color>0FFFFFF</color>\n')
        fileWriter.write('</LabelStyle>\n')
        fileWriter.write('<BalloonStyle>\n')
        fileWriter.write('</BalloonStyle>\n')
        fileWriter.write('<ListStyle>\n')
        fileWriter.write('</ListStyle>\n')
        fileWriter.write('</Style>\n')
        fileWriter.write('<Style id="sh_shaded_dot4">\n')
        fileWriter.write('<IconStyle>\n')
        fileWriter.write('	<color>6414B4FF</color>\n')
        fileWriter.write(' <scale>1.4</scale>\n')
        fileWriter.write('	<Icon>\n')
        fileWriter.write('		<href>http://maps.google.com/mapfiles/kml/shapes/shaded_dot.png</href>\n')
        fileWriter.write('	</Icon>\n')
        fileWriter.write(' </IconStyle>\n')
        fileWriter.write('<LabelStyle>\n')
        fileWriter.write('	<color>0FFFFFF</color>\n')
        fileWriter.write('</LabelStyle>\n')
        fileWriter.write('<BalloonStyle>\n')
        fileWriter.write('</BalloonStyle>\n')
        fileWriter.write('<ListStyle>\n')
        fileWriter.write('   </ListStyle>\n')
        fileWriter.write('</Style>\n')

        fileWriter.write('<StyleMap id="msn_shaded_dot5">\n')
        fileWriter.write('	<Pair>\n')
        fileWriter.write('	<key>normal</key>\n')
        fileWriter.write('  <styleUrl>#sn_shaded_dot5</styleUrl>\n')
        fileWriter.write('	</Pair>\n')
        fileWriter.write('   <Pair>\n')
        fileWriter.write('   <key>highlight</key>\n')
        fileWriter.write('		<styleUrl>#sh_shaded_dot5</styleUrl>\n')
        fileWriter.write('    </Pair>\n')
        fileWriter.write('</StyleMap>\n')
        fileWriter.write('     <Style id="sn_shaded_dot5">\n')
        fileWriter.write('      <IconStyle>\n')
        fileWriter.write('     <color>641478FF</color>\n')
        fileWriter.write('		<scale>1.2</scale>\n')
        fileWriter.write('     <Icon>\n')
        fileWriter.write('     <href>http://maps.google.com/mapfiles/kml/shapes/shaded_dot.png</href>\n')
        fileWriter.write('    </Icon>\n')
        fileWriter.write('	</IconStyle>\n')
        fileWriter.write('<LabelStyle>\n')
        fileWriter.write('<color>0FFFFFF</color>\n')
        fileWriter.write('</LabelStyle>\n')
        fileWriter.write('<BalloonStyle>\n')
        fileWriter.write('</BalloonStyle>\n')
        fileWriter.write('<ListStyle>\n')
        fileWriter.write('</ListStyle>\n')
        fileWriter.write('</Style>\n')
        fileWriter.write('<Style id="sh_shaded_dot5">\n')
        fileWriter.write('<IconStyle>\n')
        fileWriter.write('	<color>641478FF</color>\n')
        fileWriter.write(' <scale>1.4</scale>\n')
        fileWriter.write('	<Icon>\n')
        fileWriter.write('		<href>http://maps.google.com/mapfiles/kml/shapes/shaded_dot.png</href>\n')
        fileWriter.write('	</Icon>\n')
        fileWriter.write(' </IconStyle>\n')
        fileWriter.write('<LabelStyle>\n')
        fileWriter.write('	<color>0FFFFFF</color>\n')
        fileWriter.write('</LabelStyle>\n')
        fileWriter.write('<BalloonStyle>\n')
        fileWriter.write('</BalloonStyle>\n')
        fileWriter.write('<ListStyle>\n')
        fileWriter.write('   </ListStyle>\n')
        fileWriter.write('</Style>\n')

        fileWriter.write('<StyleMap id="msn_shaded_dot6">\n')
        fileWriter.write('	<Pair>\n')
        fileWriter.write('	<key>normal</key>\n')
        fileWriter.write('  <styleUrl>#sn_shaded_dot6</styleUrl>\n')
        fileWriter.write('	</Pair>\n')
        fileWriter.write('   <Pair>\n')
        fileWriter.write('   <key>highlight</key>\n')
        fileWriter.write('		<styleUrl>#sh_shaded_dot6</styleUrl>\n')
        fileWriter.write('    </Pair>\n')
        fileWriter.write('</StyleMap>\n')
        fileWriter.write('     <Style id="sn_shaded_dot6">\n')
        fileWriter.write('      <IconStyle>\n')
        fileWriter.write('     <color>64143CFF</color>\n')
        fileWriter.write('		<scale>1.2</scale>\n')
        fileWriter.write('     <Icon>\n')
        fileWriter.write('     <href>http://maps.google.com/mapfiles/kml/shapes/shaded_dot.png</href>\n')
        fileWriter.write('    </Icon>\n')
        fileWriter.write('	</IconStyle>\n')
        fileWriter.write('<LabelStyle>\n')
        fileWriter.write('<color>0FFFFFF</color>\n')
        fileWriter.write('</LabelStyle>\n')
        fileWriter.write('<BalloonStyle>\n')
        fileWriter.write('</BalloonStyle>\n')
        fileWriter.write('<ListStyle>\n')
        fileWriter.write('</ListStyle>\n')
        fileWriter.write('</Style>\n')
        fileWriter.write('<Style id="sh_shaded_dot6">\n')
        fileWriter.write('<IconStyle>\n')
        fileWriter.write('	<color>64143CFF</color>\n')
        fileWriter.write(' <scale>1.4</scale>\n')
        fileWriter.write('	<Icon>\n')
        fileWriter.write('		<href>http://maps.google.com/mapfiles/kml/shapes/shaded_dot.png</href>\n')
        fileWriter.write('	</Icon>\n')
        fileWriter.write(' </IconStyle>\n')
        fileWriter.write('<LabelStyle>\n')
        fileWriter.write('	<color>0FFFFFF</color>\n')
        fileWriter.write('</LabelStyle>\n')
        fileWriter.write('<BalloonStyle>\n')
        fileWriter.write('</BalloonStyle>\n')
        fileWriter.write('<ListStyle>\n')
        fileWriter.write('   </ListStyle>\n')
        fileWriter.write('</Style>\n')

        fileWriter.write('<StyleMap id="msn_shaded_dot7">\n')
        fileWriter.write('	<Pair>\n')
        fileWriter.write('	<key>normal</key>\n')
        fileWriter.write('  <styleUrl>#sn_shaded_dot7</styleUrl>\n')
        fileWriter.write('	</Pair>\n')
        fileWriter.write('   <Pair>\n')
        fileWriter.write('   <key>highlight</key>\n')
        fileWriter.write('		<styleUrl>#sh_shaded_dot7</styleUrl>\n')
        fileWriter.write('    </Pair>\n')
        fileWriter.write('</StyleMap>\n')
        fileWriter.write('     <Style id="sn_shaded_dot7">\n')
        fileWriter.write('      <IconStyle>\n')
        fileWriter.write('     <color>641400FF</color>\n')
        fileWriter.write('		<scale>1.2</scale>\n')
        fileWriter.write('     <Icon>\n')
        fileWriter.write('     <href>http://maps.google.com/mapfiles/kml/shapes/shaded_dot.png</href>\n')
        fileWriter.write('    </Icon>\n')
        fileWriter.write('	</IconStyle>\n')
        fileWriter.write('<LabelStyle>\n')
        fileWriter.write('<color>0FFFFFF</color>\n')
        fileWriter.write('</LabelStyle>\n')
        fileWriter.write('<BalloonStyle>\n')
        fileWriter.write('</BalloonStyle>\n')
        fileWriter.write('<ListStyle>\n')
        fileWriter.write('</ListStyle>\n')
        fileWriter.write('</Style>\n')
        fileWriter.write('<Style id="sh_shaded_dot7">\n')
        fileWriter.write('<IconStyle>\n')
        fileWriter.write('	<color>641400FF</color>\n')
        fileWriter.write(' <scale>1.4</scale>\n')
        fileWriter.write('	<Icon>\n')
        fileWriter.write('		<href>http://maps.google.com/mapfiles/kml/shapes/shaded_dot.png</href>\n')
        fileWriter.write('	</Icon>\n')
        fileWriter.write(' </IconStyle>\n')
        fileWriter.write('<LabelStyle>\n')
        fileWriter.write('	<color>0FFFFFF</color>\n')
        fileWriter.write('</LabelStyle>\n')
        fileWriter.write('<BalloonStyle>\n')
        fileWriter.write('</BalloonStyle>\n')
        fileWriter.write('<ListStyle>\n')
        fileWriter.write('   </ListStyle>\n')
        fileWriter.write('</Style>\n')


        GeometryFactory gf = new GeometryFactory()
        CRSFactory cRSFactory = new CRSFactory()

        RegistryManager registryManager = cRSFactory.getRegistryManager()
        registryManager.addRegistry(new EPSGRegistry())

        CoordinateReferenceSystem inputCRS = cRSFactory.getCRS("EPSG:2154")
        CoordinateReferenceSystem targetCRS = cRSFactory.getCRS("EPSG:4326")
        List<GeometryEditor.CoordinateOperation> op = CoordinateOperationFactory.createCoordinateOperations((GeodeticCRS) inputCRS, (GeodeticCRS) targetCRS)

        Iterator it = Result_by_receivers.entrySet().iterator()
        while (it.hasNext()) {
            Map.Entry pair = (Map.Entry) it.next()
            Integer idReceiver = pair.getKey()
            Double Res = pair.getValue()
            Coordinate coordinate = receivers.get(idReceiver)

            fileWriter.write('\t<Placemark>\n')
            fileWriter.write('\t<name>Result</name>')
            fileWriter.write('\t\t<description>'+ Res.toString() +'</description>\n')
            fileWriter.write('\t\t<styleUrl>#msn_shaded_dot'+ (Math.round((Res-50)*15/100)).toString() +'</styleUrl>\n')
            fileWriter.write('\t\t<Point>\n')
            fileWriter.write('\t\t	<altitudeMode>relativeToGround</altitudeMode>\n')
            fileWriter.write('\t\t<gx:drawOrder>1</gx:drawOrder>\n')
            fileWriter.write('\t\t<coordinates>')

            double z = mesh.getHeightAtPosition(coordinate)

            Geometry outPutGeom = (Point) gf.createPoint(coordinate).copy()
            outPutGeom.geometryChanged()
            outPutGeom.apply(new ST_Transform.CRSTransformFilter(op.get(0)))
            outPutGeom.setSRID(4326)
            Coordinate coordinate2 = outPutGeom.getCoordinate()
            fileWriter.write("\t\t\t" + coordinate2.x.toString() + "," + coordinate2.y.toString() + "," + (coordinate.z - z).toString() + "\n")


            fileWriter.write('\t\t\t</coordinates>\n')
            fileWriter.write('\t\t\t</Point>\n')
            fileWriter.write('\t</Placemark>\n')
        }
        it.remove()
        fileWriter.write('\t</Document>\n')
        fileWriter.write('</kml>')
        fileWriter.close()
}

/**
 * Permet de tracer les rayons dans un vtk lisible dans paraview
 * @param filename
 * @param propDataOut
 * @throws IOException
 */
    private void writeVTK(String filename, ComputeRaysOut propDataOut) throws IOException {


        FileWriter fileWriter = new FileWriter(filename)
        fileWriter.write("# vtk DataFile Version 2.0\n")
        fileWriter.write("PropagationPath\n")
        fileWriter.write("ASCII\n")
        fileWriter.write("DATASET POLYDATA\n")
        int nbPoints = 0
        for (int j = 0; j < propDataOut.propagationPaths.size(); j++) {
            nbPoints = nbPoints + propDataOut.propagationPaths.get(j).getPointList().size()
        }
        fileWriter.write("\n")
        fileWriter.write("POINTS " + String.valueOf(nbPoints) + " float\n")

        GeometryFactory geometryFactory = new GeometryFactory()
        List<Coordinate> coordinates = new ArrayList<>()
        for (int j = 0; j < propDataOut.propagationPaths.size(); j++) {
            for (PropagationPath.PointPath p : propDataOut.propagationPaths.get(j).getPointList()) {
                coordinates.add(p.coordinate)
                fileWriter.write(String.valueOf(p.coordinate.x) + " " + String.valueOf(p.coordinate.y) + " " + String.valueOf(p.coordinate.z) + "\n")
            }
        }

        fileWriter.write("\n")
        fileWriter.write("LINES " + String.valueOf(propDataOut.propagationPaths.size()) + " " + String.valueOf(nbPoints + propDataOut.propagationPaths.size()) + "\n")
        int i = 0
        for (int j = 0; j < propDataOut.propagationPaths.size(); j++) {
            fileWriter.write(String.valueOf(propDataOut.propagationPaths.get(j).getPointList().size()))

            for (PropagationPath.PointPath p : propDataOut.propagationPaths.get(j).getPointList()) {
                fileWriter.write(" " + String.valueOf(i))
                i++
            }
            fileWriter.write("\n")
        }


        fileWriter.close()
    }

/**
 * permet de faire un fichier ply avec la topographie pour lecture dans paraview
 * @param filename
 * @param mesh
 * @throws IOException
 * @throws LayerDelaunayError
 */
    private void writePLY(String filename, MeshBuilder mesh) throws IOException, LayerDelaunayError {
        PointsMerge pointsMerge = new PointsMerge(0.01)
        List<Geometry> triVertices2 = new ArrayList<>()
        Map<String, Integer> vertices2 = new HashMap<>()
        List<Coordinate> vertices3 = new ArrayList<>()
        GeometryFactory geometryFactory = new GeometryFactory()
        int k = 0
        for (MeshBuilder.PolygonWithHeight polygon : mesh.getPolygonWithHeight()) {
            GeometryCollection buildingExtruded = GeometryExtrude.extrudePolygonAsGeometry((Polygon) polygon.getGeometry(), polygon.getHeight())
            addGeometry(triVertices2, buildingExtruded)
            for (Coordinate coordinate : buildingExtruded.getCoordinates()) {
                vertices2.put(coordinate.toString(), k)
                vertices3.add(coordinate)
                k++
            }

        }
        int vertexCountG = mesh.getVertices().size()
        int vertexCountB = vertices3.size()
        int faceCountG = mesh.getTriangles().size()
        int faceCountB = triVertices2.size()
        int vertexCount = vertexCountG + vertexCountB
        int faceCount = faceCountG + faceCountB
        FileWriter fileWriter = new FileWriter(filename)
        fileWriter.write("ply\n")
        fileWriter.write("format ascii 1.0\n")
        fileWriter.write("element vertex " + vertexCount + "\n")
        fileWriter.write("property float x\n")
        fileWriter.write("property float y\n")
        fileWriter.write("property float z\n")
        fileWriter.write("property uchar green\n")
        fileWriter.write("property uchar red\n")
        fileWriter.write("property uchar blue\n")
        fileWriter.write("element face " + faceCount + "\n")
        fileWriter.write("property list uchar int vertex_index\n")
        fileWriter.write("end_header\n")

        for (int i = 0; i < vertexCountG; i++) {
            fileWriter.write(mesh.getVertices().get(i).x + " " + mesh.getVertices().get(i).y + " " + mesh.getVertices().get(i).z + " " + "255 0 0\n")
        }
        // Iterating over values only
        for (Coordinate vertice : vertices3) {
            fileWriter.write(vertice.x + " " + vertice.y + " " + vertice.z + " " + "0 0 255\n")
        }

        for (int i = 0; i < faceCountG; i++) {
            fileWriter.write("3 " + mesh.getTriangles().get(i).getA() + " " + mesh.getTriangles().get(i).getB() + " " + (mesh.getTriangles().get(i).getC()) + "\n")
        }
        for (int i = 0; i < faceCountB; i++) {
            Coordinate[] coordinates = triVertices2.get(i).getCoordinates()
            fileWriter.write(coordinates.length + " ")
            for (int j = 0; j < coordinates.length; j++) {
                fileWriter.write((vertexCountG + vertices2.get(coordinates[j].toString())) + " ")
            }
            fileWriter.write("\n")
        }
        fileWriter.close()
    }

/**
 * permet de faire un fichier KML avec la topographie pour lecture dans GoogleEarth
 * @param filename
 * @param mesh
 * @throws IOException
 * @throws LayerDelaunayError
 */

    private void writeKML(String filename, MeshBuilder mesh, FastObstructionTest manager) throws IOException, LayerDelaunayError {

        PointsMerge pointsMerge = new PointsMerge(0.01)
        List<Geometry> triVertices2 = new ArrayList<>()
        Map<String, Integer> vertices2 = new HashMap<>()
        List<Coordinate> vertices3 = new ArrayList<>()
        GeometryFactory geometryFactory = new GeometryFactory()

        int vertexCountG = mesh.getVertices().size()
        int vertexCountB = vertices3.size()
        int faceCountG = mesh.getTriangles().size()
        int faceCountB = triVertices2.size()
        int vertexCount = vertexCountG + vertexCountB
        int faceCount = faceCountG + faceCountB
        FileWriter fileWriter = new FileWriter(filename)
        fileWriter.write('<?xml version="1.0" encoding="UTF-8"?>\n')
        fileWriter.write('<kml xmlns="http://www.opengis.net/kml/2.2">\n')
        fileWriter.write('<Document>\n')
        CRSFactory cRSFactory = new CRSFactory()
        GeometryFactory gf = new GeometryFactory()
        RegistryManager registryManager = cRSFactory.getRegistryManager()
        registryManager.addRegistry(new EPSGRegistry())

        CoordinateReferenceSystem inputCRS = cRSFactory.getCRS("EPSG:2154")
        CoordinateReferenceSystem targetCRS = cRSFactory.getCRS("EPSG:4326")
        List<GeometryEditor.CoordinateOperation> op = CoordinateOperationFactory.createCoordinateOperations((GeodeticCRS) inputCRS, (GeodeticCRS) targetCRS)



        int k = 0
        for (MeshBuilder.PolygonWithHeight polygon : mesh.getPolygonWithHeight()) {
            fileWriter.write('\t<Placemark>\n')
            fileWriter.write('\t\t<name>The Pentagon</name>\n')
            fileWriter.write('\t\t<Polygon>\n')
            fileWriter.write('\t\t\t<extrude>1</extrude>\n')
            fileWriter.write('\t\t\t<altitudeMode>relativeToGround</altitudeMode>\n')
            fileWriter.write('\t\t\t<outerBoundaryIs>\n')
            fileWriter.write('\t\t\t<LinearRing>\n')
            fileWriter.write('\t\t\t<coordinates>\n')
            double height = polygon.getHeight()

            //Geometry building = (Polygon) polygon.getGeometry()
            //Geometry outPutGeom = building.copy()
            //outPutGeom.geometryChanged()
            //outPutGeom.apply(new ST_Transform.CRSTransformFilter(op.get(0)))
            //outPutGeom.setSRID(4326)

            // addGeometry(triVertices2, buildingExtruded)

            for (Coordinate coordinate : polygon.getGeometry().getCoordinates()) {

                double z = manager.getHeightAtPosition(coordinate)
                Geometry outPutGeom = (Point) gf.createPoint(coordinate).copy()
                outPutGeom.geometryChanged()
                outPutGeom.apply(new ST_Transform.CRSTransformFilter(op.get(0)))
                outPutGeom.setSRID(4326)
                Coordinate coordinate2 = outPutGeom.getCoordinate()

                fileWriter.write("\t\t\t" + coordinate2.x.toString() + "," + coordinate2.y.toString() + "," + (height).toString() + "\n")
            }

            fileWriter.write('\t\t\t</coordinates>\n')
            fileWriter.write('\t\t\t</LinearRing>\n')
            fileWriter.write('\t\t\t</outerBoundaryIs>\n')
            fileWriter.write('\t\t</Polygon>\n')
            fileWriter.write('\t</Placemark>\n')

        }

        fileWriter.write('\t</Document>\n')
        fileWriter.write('</kml>')
        fileWriter.close()
    }
/**
 * Permet d'importer la topo depuis la base postgre
 * @param connection
 * @param mesh
 * @param buildingsTableName
 * @param zField
 * @throws SQLException
 */
    protected void fetchCellDem(Connection connection, MeshBuilder mesh, String buildingsTableName, String zField) throws SQLException {

        PreparedStatement st = connection.prepareStatement(
                "SELECT the_geom, " + zField + " FROM " +
                        buildingsTableName)

        SpatialResultSet rs = st.executeQuery().unwrap(SpatialResultSet.class)
        while (rs.next()) {
            Geometry pt = rs.getGeometry()
            if (pt != null) {
                mesh.addTopographicPoint(new Coordinate(pt.getCoordinate().x, pt.getCoordinate().y, rs.getDouble(zField)))
            }
        }

    }

/**
 * permet d'importe les recepteurs depuis la base postgre
 * @param connection
 * @param receiverTableName
 * @param receivers
 * @param skipReceivers
 */
    protected List<Long> fetchCellReceiver_withindex(Connection connection, String receiverTableName, List<Coordinate> receivers, Set<Long> skipReceivers) {
        // Fetch receivers
        //int intPk = JDBCUtilities.getIntegerPrimaryKey(connection, receiverTableName)
        List<Long> receiversPk = new ArrayList<>()
        String pkSelect
        //if (intPk >= 1) {
        //    pkSelect = ", " + JDBCUtilities.getFieldName(connection.getMetaData(), receiverTableName, intPk)
        //} else {
        pkSelect = ", id "
        //}
        //List<String> geomFields = SFSUtilities.getGeometryFields(connection, TableLocation.parse(receiverTableName))

        PreparedStatement st = connection.prepareStatement("SELECT the_geom" + pkSelect + " FROM " + receiverTableName)
        SpatialResultSet rs = st.executeQuery().unwrap(SpatialResultSet.class)
        while (rs.next()) {
            if (!pkSelect.isEmpty()) {
                long receiverPk2 = rs.getLong(2)
                if (skipReceivers.contains(receiverPk2)) {
                    continue
                }
                receiversPk.add(receiverPk2)
            }
            Geometry pt = rs.getGeometry()
            if (pt != null) {
                receivers.add(pt.getCoordinate())
            }


        }
        return receiversPk
    }

/**
 * permet d'importer les sources depuis la base postgre
 * @param connection
 * @param allSourceGeometries
 * @param sourceGeometries
 * @param sourcesTableName
 * @param sourcePk
 * @param wj_sources
 * @param sourcesIndex
 * @throws SQLException
 */
    protected void fetchCellSource_withindex(Connection connection, List<Geometry> allSourceGeometries,
                                             List<Geometry> sourceGeometries, String sourcesTableName, List<Long> sourcePk, List<ArrayList<Double>> wj_sources, QueryGeometryStructure sourcesIndex) throws SQLException {
        double[] db_field_ids = [3, 4, 5, 6, 7, 8, 9, 10]
        int idSource = 0
        PreparedStatement st = connection.prepareStatement("SELECT * FROM " + sourcesTableName)
        SpatialResultSet rs = st.executeQuery().unwrap(SpatialResultSet.class)
        while (rs.next()) {
            Geometry geo = rs.getGeometry()
            // todo get primarykey as for receiver, for the moment take only the first column
            sourcePk.add(rs.getLong(1))
            if (geo != null) {
                ArrayList<Double> wj_spectrum = new ArrayList<>()

                double sumPow = 0
                for (Integer idcol : db_field_ids) {
                    double wj = DbaToW(rs.getDouble(idcol))
                    wj_spectrum.add(wj)
                    sumPow += wj
                }
                if (allSourceGeometries != null) {
                    allSourceGeometries.add(geo)
                }
                if (sumPow > 0) {
                    wj_sources.add(wj_spectrum)
                    sourcesIndex.appendGeometry(geo, idSource)
                    sourceGeometries.add(geo)
                    idSource++
                }
            }
        }
    }

/**
 * Permet d'importer les types de sol depuis la base postgre

 * @param connection
 * @param geoWithSoil
 * @param soilTableName
 * @param gField
 * @throws SQLException
 */
    protected void fetchCellSoilAreas(Connection connection, List<GeoWithSoilType> geoWithSoil, String soilTableName, String gField) throws SQLException {
        PreparedStatement st = connection.prepareStatement("SELECT the_geom, " + gField + " FROM " + soilTableName)
        SpatialResultSet rs = st.executeQuery().unwrap(SpatialResultSet.class)
        while (rs.next()) {
            Geometry poly = rs.getGeometry()
            if (poly != null) {
                geoWithSoil.add(new GeoWithSoilType(poly, rs.getDouble(gField)))
            }
        }

    }

/**
 * permet d'importer les bâtiment depuis la base postgre
 * @param connection
 * @param mesh
 * @param buildingsTableName
 * @param ALPHA_FIELD_NAME
 * @param heightField
 * @throws SQLException
 */
    protected void fetchCellBuildings(Connection connection, MeshBuilder mesh, String buildingsTableName, String ALPHA_FIELD_NAME, String heightField) throws SQLException {
        boolean fetchAlpha = JDBCUtilities.hasField(connection, buildingsTableName, ALPHA_FIELD_NAME)

        List<String> geomFields = SFSUtilities.getGeometryFields(connection, TableLocation.parse(buildingsTableName))
        String topoGeomName = geomFields.get(0)

        double wallAbsorption = 0.1
        String additionalQuery
        if (!heightField.isEmpty()) {
            additionalQuery = ", " + TableLocation.quoteIdentifier(heightField)
        }
        if (fetchAlpha) {
            additionalQuery = ", " + ALPHA_FIELD_NAME
        }
        PreparedStatement st = connection.prepareStatement(
                "SELECT " + TableLocation.quoteIdentifier(topoGeomName) + ", " + heightField + " FROM " +
                        buildingsTableName)

        SpatialResultSet rs = st.executeQuery().unwrap(SpatialResultSet.class)
        while (rs.next()) {
            //if we don't have height of building
            Geometry building = rs.getGeometry()
            if (building != null) {
                if (building instanceof Polygon || building instanceof MultiPolygon) {
                    mesh.addGeometry(building,
                            heightField.isEmpty() ? Double.MAX_VALUE : rs.getDouble(heightField),
                            fetchAlpha ? rs.getDouble(ALPHA_FIELD_NAME) : wallAbsorption)
                }
            }
        }
    }

/**
 * permet de copier coller les fichiers resultats
 * @param source
 * @param dest
 * @throws IOException
 */
    private static void copyFileUsingStream(File source, File dest) throws IOException {
        InputStream is = null
        OutputStream os = null
        try {
            is = new FileInputStream(source)
            os = new FileOutputStream(dest)
            byte[] buffer = new byte[1024]
            int length
            while ((length = is.read(buffer)) > 0) {
                os.write(buffer, 0, length)
            }
        } finally {
            is.close()
            os.close()
        }
    }

}
