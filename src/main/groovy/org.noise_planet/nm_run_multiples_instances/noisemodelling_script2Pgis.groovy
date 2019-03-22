// author : Aumond Pierre

// Ce code permet de faire une série de simulation à partir :
// - d'un plan d'expérience appelé Exp_comp.csv
// - et d'un fichier de config Config.csv
// Le résultat se trouve dans un dossier zippé
// Il s'agit d'une version Postgre/Postgis (lecture/ecriture plus lent qu'avec h2gis mais jointures beaucoup plus rapides)

// Importation des librairies nécessaire au code
import groovy.sql.Sql
import groovy.transform.Field
import org.apache.commons.io.FileUtils
import org.h2gis.functions.io.shp.SHPDriverFunction
import org.h2gis.functions.spatial.volume.GeometryExtrude
import org.h2gis.utilities.JDBCUtilities
import org.h2gis.utilities.SFSUtilities
import org.h2gis.utilities.SpatialResultSet
import org.h2gis.utilities.TableLocation
import org.locationtech.jts.geom.Coordinate
import org.locationtech.jts.geom.Envelope
import org.locationtech.jts.geom.Geometry
import org.locationtech.jts.geom.GeometryCollection
import org.locationtech.jts.geom.GeometryFactory
import org.locationtech.jts.geom.LineString
import org.locationtech.jts.geom.MultiPolygon
import org.locationtech.jts.geom.Polygon
import org.locationtech.jts.io.WKTWriter
import org.orbisgis.noisemap.core.ComputeRays
import org.orbisgis.noisemap.core.ComputeRaysOut
import org.orbisgis.noisemap.core.EvaluateAttenuationCnossos
import org.orbisgis.noisemap.core.FastObstructionTest
import org.orbisgis.noisemap.core.GeoWithSoilType
import org.orbisgis.noisemap.core.LayerDelaunayError
import org.orbisgis.noisemap.core.MeshBuilder
import org.orbisgis.noisemap.core.PointsMerge
import org.orbisgis.noisemap.core.PropagationDebugInfo
import org.orbisgis.noisemap.core.PropagationPath
import org.orbisgis.noisemap.core.PropagationProcessData
import org.orbisgis.noisemap.core.PropagationProcessOut_Att_f
import org.orbisgis.noisemap.core.PropagationProcessPathData
import org.orbisgis.noisemap.core.QueryGeometryStructure
import org.orbisgis.noisemap.core.QueryQuadTree
import org.h2gis.functions.io.csv.CSVDriverFunction
import org.h2gis.utilities.wrapper.ConnectionWrapper

import java.sql.Connection
import java.sql.DriverManager
import java.sql.PreparedStatement
import java.sql.SQLException
import java.sql.Statement

import org.orbisgis.noisemap.core.EvaluateRoadSourceCnossos
import org.orbisgis.noisemap.core.RSParametersCnossos

import groovy.time.*
import java.nio.file.Files
import java.util.zip.ZipEntry
import java.util.zip.ZipOutputStream
import org.h2gis.api.EmptyProgressVisitor



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
    double[] sum = new double[array1.length];
    for (int i = 0; i < array1.length; i++) {
        sum[i] = wToDba(DbaToW(array1[i]) + DbaToW(array2[i]));
    }
    return sum;
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
        geom.add((Polygon) polygon);
    } else {
        for (int i = 0; i < polygon.getNumGeometries(); i++) {
            addGeometry(geom, polygon.getGeometryN(i));
        }
    }

}

/**
 * Combine l'ensemble des rayons en prenant en compte favorable defavorable
 * todo pour le moment on prend juste p, mais en fait pour chacun des rayons il faut regarder la direction principale puis combiner avec les favroses
 * @param aGlobal
 * @param propData
 * @param propDataOut
 * @param p
 */
private void computeWithMeteo(HashMap<Integer, double[]> aGlobal, PropagationProcessPathData propData, ComputeRaysOut propDataOut, double p) {
    EvaluateAttenuationCnossos evaluateAttenuationCnossos = new EvaluateAttenuationCnossos()
    for (PropagationPath propath : propDataOut.propagationPaths) {
        aGlobalMeteo = aGlobal.get(propath.idSource)
        propath.setFavorable(false)
        evaluateAttenuationCnossos.evaluate(propath, propData)
        if (aGlobalMeteo != null) {
            aGlobalMeteo = sumArray(evaluateAttenuationCnossos.getaGlobal(), aGlobalMeteo);
        } else {
            aGlobalMeteo = evaluateAttenuationCnossos.getaGlobal()
        }

        propath.setFavorable(true);
        evaluateAttenuationCnossos.evaluate(propath, propData);
        aGlobalMeteo = sumArrayWithPonderation(aGlobalMeteo, evaluateAttenuationCnossos.getaGlobal(), p);
        if (aGlobal.containsKey(propath.idSource)) {
            aGlobalMeteo = sumArray(aGlobalMeteo, aGlobal.get(propath.idSource))
            aGlobal.replace(propath.idSource, aGlobalMeteo)
        } else {
            aGlobal.put(propath.idSource, aGlobalMeteo)
        }

    }
}

/**
 * Permet de tracer les rayons dans un vtk lisible dans paraview
 * @param filename
 * @param propDataOut
 * @throws IOException
 */
private void writeVTK(String filename, ComputeRaysOut propDataOut) throws IOException {


    FileWriter fileWriter = new FileWriter(filename);
    fileWriter.write("# vtk DataFile Version 2.0\n");
    fileWriter.write("PropagationPath\n");
    fileWriter.write("ASCII\n");
    fileWriter.write("DATASET POLYDATA\n");
    int nbPoints = 0;
    for (int j = 0; j < propDataOut.propagationPaths.size(); j++) {
        nbPoints = nbPoints + propDataOut.propagationPaths.get(j).getPointList().size();
    }
    fileWriter.write("\n");
    fileWriter.write("POINTS " + String.valueOf(nbPoints) + " float\n");

    GeometryFactory geometryFactory = new GeometryFactory();
    List<Coordinate> coordinates = new ArrayList<>();
    for (int j = 0; j < propDataOut.propagationPaths.size(); j++) {
        for (PropagationPath.PointPath p : propDataOut.propagationPaths.get(j).getPointList()) {
            coordinates.add(p.coordinate);
            fileWriter.write(String.valueOf(p.coordinate.x) + " " + String.valueOf(p.coordinate.y) + " " + String.valueOf(p.coordinate.z) + "\n");
        }
    }

    fileWriter.write("\n");
    fileWriter.write("LINES " + String.valueOf(propDataOut.propagationPaths.size()) + " " + String.valueOf(nbPoints + propDataOut.propagationPaths.size()) + "\n");
    int i = 0;
    for (int j = 0; j < propDataOut.propagationPaths.size(); j++) {
        fileWriter.write(String.valueOf(propDataOut.propagationPaths.get(j).getPointList().size()));

        for (PropagationPath.PointPath p : propDataOut.propagationPaths.get(j).getPointList()) {
            fileWriter.write(" " + String.valueOf(i));
            i++;
        }
        fileWriter.write("\n");
    }


    fileWriter.close();
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
    int k = 0;
    for (MeshBuilder.PolygonWithHeight polygon : mesh.getPolygonWithHeight()) {
        GeometryCollection buildingExtruded = GeometryExtrude.extrudePolygonAsGeometry((Polygon) polygon.getGeometry(), polygon.getHeight())
        addGeometry(triVertices2, buildingExtruded)
        for (Coordinate coordinate : buildingExtruded.getCoordinates()) {
            vertices2.put(coordinate.toString(), k)
            vertices3.add(coordinate)
            k++
        }

    }
    int vertexCountG = mesh.getVertices().size();
    int vertexCountB = vertices3.size();
    int faceCountG = mesh.getTriangles().size();
    int faceCountB = triVertices2.size();
    int vertexCount = vertexCountG + vertexCountB;
    int faceCount = faceCountG + faceCountB;
    FileWriter fileWriter = new FileWriter(filename);
    fileWriter.write("ply\n");
    fileWriter.write("format ascii 1.0\n");
    fileWriter.write("element vertex " + vertexCount + "\n");
    fileWriter.write("property float x\n");
    fileWriter.write("property float y\n");
    fileWriter.write("property float z\n");
    fileWriter.write("property uchar green\n");
    fileWriter.write("property uchar red\n");
    fileWriter.write("property uchar blue\n");
    fileWriter.write("element face " + faceCount + "\n");
    fileWriter.write("property list uchar int vertex_index\n");
    fileWriter.write("end_header\n");

    for (int i = 0; i < vertexCountG; i++) {
        fileWriter.write(mesh.getVertices().get(i).x + " " + mesh.getVertices().get(i).y + " " + mesh.getVertices().get(i).z + " " + "255 0 0\n")
    }
    // Iterating over values only
    for (Coordinate vertice : vertices3) {
        fileWriter.write(vertice.x + " " + vertice.y + " " + vertice.z + " " + "0 0 255\n")
    }

    for (int i = 0; i < faceCountG; i++) {
        fileWriter.write("3 " + mesh.getTriangles().get(i).getA() + " " + mesh.getTriangles().get(i).getB() + " " + (mesh.getTriangles().get(i).getC()) + "\n");
    }
    for (int i = 0; i < faceCountB; i++) {
        Coordinate[] coordinates = triVertices2.get(i).getCoordinates();
        fileWriter.write(coordinates.length + " ");
        for (int j = 0; j < coordinates.length; j++) {
            fileWriter.write((vertexCountG + vertices2.get(coordinates[j].toString())) + " ");
        }
        fileWriter.write("\n");
    }
    fileWriter.close();
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
protected void fetchCellReceiver_withindex(Connection connection, String receiverTableName, List<Coordinate> receivers, Set<Long> skipReceivers) {
    // Fetch receivers

    int intPk = JDBCUtilities.getIntegerPrimaryKey(connection, receiverTableName)
    String pkSelect = ""
    if (intPk >= 1) {
        pkSelect = ", " + JDBCUtilities.getFieldName(connection.getMetaData(), receiverTableName, intPk);
    } else {
        pkSelect = ", id "
    }
    PreparedStatement st = connection.prepareStatement("SELECT the_geom" + pkSelect + " FROM " + receiverTableName)
    SpatialResultSet rs = st.executeQuery().unwrap(SpatialResultSet.class)
    while (rs.next()) {
        if (!pkSelect.isEmpty()) {
            long receiverPk = rs.getLong(2);
            if (skipReceivers.contains(receiverPk)) {
                continue;
            }
            receiversPk.add(receiverPk)
        }
        Geometry pt = rs.getGeometry();
        if (pt != null) {
            receivers.add(pt.getCoordinate());
        }


    }
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
    db_field_ids = [3, 4, 5, 6, 7, 8, 9, 10]
    int idSource = 0;
    PreparedStatement st = connection.prepareStatement("SELECT * FROM " + sourcesTableName)
    SpatialResultSet rs = st.executeQuery().unwrap(SpatialResultSet.class)
    while (rs.next()) {
        Geometry geo = rs.getGeometry();
        // todo get primarykey as for receiver, for the moment take only the first column
        sourcePk.add(rs.getLong(1));
        if (geo != null) {
            ArrayList<Double> wj_spectrum = new ArrayList<>();

            double sumPow = 0;
            for (Integer idcol : db_field_ids) {
                double wj = DbaToW(rs.getDouble(idcol));
                wj_spectrum.add(wj);
                sumPow += wj;
            }
            if (allSourceGeometries != null) {
                allSourceGeometries.add(geo);
            }
            if (sumPow > 0) {
                wj_sources.add(wj_spectrum);
                sourcesIndex.appendGeometry(geo, idSource);
                sourceGeometries.add(geo);
                idSource++;
            }
        }
    }
}

/**
 * Permet d'importer les types de sol depuis la base postgre
 * todo pour le moment marche pas bien
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
    boolean fetchAlpha = JDBCUtilities.hasField(connection, buildingsTableName, ALPHA_FIELD_NAME);

    List<String> geomFields = SFSUtilities.getGeometryFields(connection, TableLocation.parse(buildingsTableName))
    String topoGeomName = geomFields.get(0)

    double wallAbsorption = 0.1;
    String additionalQuery = "";
    if (!heightField.isEmpty()) {
        additionalQuery = ", " + TableLocation.quoteIdentifier(heightField);
    }
    if (fetchAlpha) {
        additionalQuery = ", " + ALPHA_FIELD_NAME;
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

/**
 * La, c'est vraiment le début du code, attention, attention...
 */

////////////////////////////////////////////
// Paramètres d'entrée et initialisations //
////////////////////////////////////////////

//
String workspace_input = "D:\\aumond\\Documents\\CENSE\\WP2\\Analyses\\Incertitudes\\NM_Pierre\\input"
String workspace_output = "D:\\aumond\\Documents\\CENSE\\WP2\\Analyses\\Incertitudes\\NM_Pierre\\output"
// le workspace input et output doivent respecter une arborescence specifique
// data
// config

boolean clean = false // reset the postgis database
boolean printply = false // sort un ply de la zone
int nvar = 0
int nr = 0
int nSimu = 0   // Nombre de simulations
int batrec = 3  // Localisation des récepteurs :  (1) : bâtiments, (2) : récepteurs, (3) : bâtiments + récepteurs
int n_comp = 0
int i_read = 1  // Nombre de lignes d'entête
String h2_url = "jdbc:postgresql_h2://127.0.0.1:5433/noisemodelling"    // Adresse de la base de données PostGreSQL
user_name = "postgres"
user_password = "12345678"

String zip_filename_root = "results_C"      // Racine du nom de l'archive des résultats (.zip)
String config_filename = "\\config\\Config_GG.csv"    // Nom du fichier de configuration
String exp_filename =  "\\config\\Exp_comp1m_GG.csv"

// Noms des tables en entrée et des attributs
String zone_name = "zone_cense_2km"
String buildings_table_name = "cerema_ppbe_lorient_batiment_pop"
String alpha_field = ""
String height_field_name = "height"
String receivers_table_name="receivers"
String sources_table_name ="roads_src_zone"

// Paramètres de propagation
int reflexion_order = 0
int diffraction_order = 0
double max_src_dist = 500
double max_ref_dist = 500
double min_ref_dist = 1.0
double wall_alpha = 0.2
double forget_source = 0.1
int cell_id = 0
boolean compute_vertical_diffraction = false


// Ouverture de la base de donnees PostGreSQL qui recueillera les infos
ArrayList<String> urls = new ArrayList<>()
ArrayList<Connection> connections = new ArrayList<>()

// Base PostGIS
Class.forName("org.orbisgis.postgis_jts.Driver")
Connection connection = new ConnectionWrapper(DriverManager.getConnection(h2_url, user_name, user_password))
connections.add(connection)

def url2 = h2_url
def driver = 'org.orbisgis.postgis_jts.Driver'
def user = user_name
def password = user_password
def sql = Sql.newInstance(url2, user, password, driver)

// Dossier des résultats
File newFile = new File(workspace_output + "/config/" )
FileUtils.cleanDirectory(newFile)
newFile = new File(workspace_output + "/data/" )
FileUtils.cleanDirectory(newFile)

try {
    System.out.println("Init")
    CSVDriverFunction csv = new CSVDriverFunction()
    SHPDriverFunction shp = new SHPDriverFunction()
    // Fichier de config
    File src_config = new File(workspace_input + config_filename)
    File dest_config = new File(workspace_output + config_filename)
    copyFileUsingStream(src_config, dest_config)
    // Fichier de plan d'expérience
    File src_exp = new File(workspace_input + exp_filename)
    File dest_exp = new File(workspace_output + exp_filename)
    copyFileUsingStream(src_exp, dest_exp)

    // Lecture du fichier de config
    new File(workspace_output + config_filename).splitEachLine(",") {
        fields ->
            switch (i_read) {
                case 1:
                    nvar = fields[0].toInteger()
                case 2:
                    nr = fields[0].toInteger()
                case 3:
                    nSimu = fields[0].toInteger()
                case 4:
                    n_comp = fields[0].toInteger()
                default: break
            }
            i_read = i_read + 1
    }

    CSVDriverFunction exp = new CSVDriverFunction()
    //////////////////////
    // Import file text
    //////////////////////

    ArrayList<Integer> Simu = new ArrayList<Integer>()
    ArrayList<Double> CalcTime = new ArrayList<Double>()

    // Init des variables
    ArrayList<Double> Deb_d_lv = new ArrayList<Double>()
    ArrayList<Double> Deb_d_mv = new ArrayList<Double>()
    ArrayList<Double> Deb_d_hv = new ArrayList<Double>()
    ArrayList<Double> Deb_d_wav = new ArrayList<Double>()
    ArrayList<Double> Deb_d_wbv = new ArrayList<Double>()

    ArrayList<Double> Deb_e_lv = new ArrayList<Double>()
    ArrayList<Double> Deb_e_mv = new ArrayList<Double>()
    ArrayList<Double> Deb_e_hv = new ArrayList<Double>()
    ArrayList<Double> Deb_e_wav = new ArrayList<Double>()
    ArrayList<Double> Deb_e_wbv = new ArrayList<Double>()

    ArrayList<Double> Deb_n_lv = new ArrayList<Double>()
    ArrayList<Double> Deb_n_mv = new ArrayList<Double>()
    ArrayList<Double> Deb_n_hv = new ArrayList<Double>()
    ArrayList<Double> Deb_n_wav = new ArrayList<Double>()
    ArrayList<Double> Deb_n_wbv = new ArrayList<Double>()

    ArrayList<Double> Speed_d_lv = new ArrayList<Double>()
    ArrayList<Double> Speed_d_mv = new ArrayList<Double>()
    ArrayList<Double> Speed_d_hv = new ArrayList<Double>()
    ArrayList<Double> Speed_d_wav = new ArrayList<Double>()
    ArrayList<Double> Speed_d_wbv = new ArrayList<Double>()

    ArrayList<Double> Speed_e_lv = new ArrayList<Double>()
    ArrayList<Double> Speed_e_mv = new ArrayList<Double>()
    ArrayList<Double> Speed_e_hv = new ArrayList<Double>()
    ArrayList<Double> Speed_e_wav = new ArrayList<Double>()
    ArrayList<Double> Speed_e_wbv = new ArrayList<Double>()

    ArrayList<Double> Speed_n_lv = new ArrayList<Double>()
    ArrayList<Double> Speed_n_mv = new ArrayList<Double>()
    ArrayList<Double> Speed_n_hv = new ArrayList<Double>()
    ArrayList<Double> Speed_n_wav = new ArrayList<Double>()
    ArrayList<Double> Speed_n_wbv = new ArrayList<Double>()

    ArrayList<Double> slopeZ = new ArrayList<Double>()

    ArrayList<Double> Tempd_year = new ArrayList<Double>()
    ArrayList<Double> Tempe_year = new ArrayList<Double>()
    ArrayList<Double> Tempn_year = new ArrayList<Double>()

    ArrayList<Integer> R_Road = new ArrayList<Double>()

    ArrayList<Double> Junc_dist = new ArrayList<Double>()
    ArrayList<Integer> R_Junc = new ArrayList<Double>()

    // juste des petites variables pour contrôler le temps de calcul
    double TimeE = 0
    double TimeE2 = 0
    double TimeA = 0
    double TimeR = 0

    // Remplissage des variables avec le contenu du fichier plan d'exp
    // todo a completer avec les variables de la propagation
    new File(workspace_output + exp_filename).splitEachLine(",") { fields ->

        Deb_d_lv.add(fields[0].toFloat())
        Deb_d_mv.add(fields[1].toFloat())
        Deb_d_hv.add(fields[2].toFloat())
        Deb_d_wav.add(fields[3].toFloat())
        Deb_d_wbv.add(fields[4].toFloat())

        Deb_e_lv.add(fields[5].toFloat())
        Deb_e_mv.add(fields[6].toFloat())
        Deb_e_hv.add(fields[7].toFloat())
        Deb_e_wav.add(fields[8].toFloat())
        Deb_e_wbv.add(fields[9].toFloat())

        Deb_n_lv.add(fields[10].toFloat())
        Deb_n_mv.add(fields[11].toFloat())
        Deb_n_hv.add(fields[12].toFloat())
        Deb_n_wav.add(fields[13].toFloat())
        Deb_n_wbv.add(fields[14].toFloat())

        Speed_d_lv.add(fields[15].toFloat())
        Speed_d_mv.add(fields[16].toFloat())
        Speed_d_hv.add(fields[17].toFloat())
        Speed_d_wav.add(fields[18].toFloat())
        Speed_d_wbv.add(fields[19].toFloat())

        Speed_e_lv.add(fields[20].toFloat())
        Speed_e_mv.add(fields[21].toFloat())
        Speed_e_hv.add(fields[22].toFloat())
        Speed_e_wav.add(fields[23].toFloat())
        Speed_e_wbv.add(fields[24].toFloat())

        Speed_n_lv.add(fields[25].toFloat())
        Speed_n_mv.add(fields[26].toFloat())
        Speed_n_hv.add(fields[27].toFloat())
        Speed_n_wav.add(fields[28].toFloat())
        Speed_n_wbv.add(fields[29].toFloat())

        slopeZ.add(fields[30].toFloat())

        Tempd_year.add(fields[31].toFloat())
        Tempe_year.add(fields[32].toFloat())
        Tempn_year.add(fields[33].toFloat())

        R_Road.add(fields[34].toInteger())

        Junc_dist.add(fields[35].toFloat())
        R_Junc.add(fields[36].toInteger())
        Simu.add(fields[37].toInteger())
        i_read = i_read + 1
    }

    //////////////////////
    // Ici on debute les calculs
    //////////////////////

    // Écriture d'un fichier avec le temps de calcul
    PrintWriter CalcTime_record = new PrintWriter(new File(workspace_output + "/data/ComputationTime.txt"))
    StringBuilder sb2 = new StringBuilder()

    // If clean = true, recompute input tables
    if (clean) {
        System.out.println("Clean database...")
        // Nettoyage de la base de données
        sql.execute(new File("../../../sql/Cleaning_database.sql").text)
        // Préparation de la zone d'étude
        sql.execute("drop table if exists zone;")
        def zone_sql_select = "create table zone as select * from " + zone_name + ";"
        sql.execute(zone_sql_select)
        //sql.execute("alter table zone rename column geom to the_geom;")
        sql.execute("alter table zone alter column the_geom type geometry using ST_SetSRID(the_geom, 2154);")
        sql.execute("create index zone_the_geom_gist on zone using GIST (the_geom);")
        // Préparation des tables d'entée types batiments, routes, etc., sur la zone (Need postgis extension sfcgal)
        sql.execute(new File("../../../sql/LoadTables2Pgis.sql").text)
    }
    // ------------------------------------------------------------ //
    // ----------- Initialisation et import des données ----------- //
    // -----------        (sol et des bâtiments)        ----------- //
    // ------------------------------------------------------------ //

    System.out.println("Import Tables...")
    // Import GeometrySoilType
    List<GeoWithSoilType> geoWithSoilTypeList = new ArrayList<>()
    System.out.println("LandUse...")
    fetchCellSoilAreas(connection, geoWithSoilTypeList, "land_use_zone_capteur2", "g")

    MeshBuilder mesh = new MeshBuilder()

    // Import Buildings
    System.out.println("Buildings...")
    fetchCellBuildings(connection, mesh, "buildings_zone", "", "height")

    // Enveloppe de la scène
    Envelope cellEnvelope = mesh.getEnvelope()

    // Importer le Sol
    System.out.println("Altitudes...")
    //sql.execute("alter table full_dem_lite2 rename column geom to the_geom;")
    fetchCellDem(connection, mesh, "dem_lite2", "contour")

    mesh.finishPolygonFeeding(cellEnvelope)
    FastObstructionTest manager = new FastObstructionTest(mesh.getPolygonWithHeight(), mesh.getTriangles(),
            mesh.getTriNeighbors(), mesh.getVertices())

    // Ecriture de la topo dans un fichier polygon ouvrable avec paraview par exemple
    if (printply) {
        System.out.println("Impression en cours...")
        String filename2 = workspace_output + "\\" + zone_name + ".ply"
        try {
            writePLY(filename2, mesh)
        } catch (IOException e) {
            e.printStackTrace()
        }
    }

    // import sources
    List<Long> sourcesPk = new ArrayList<>()
    ArrayList<Geometry> sourceGeometries = new ArrayList<>()
    ArrayList<ArrayList<Double>> wj_sources = new ArrayList<>()
    QueryGeometryStructure sourcesIndex = new QueryQuadTree()

    Set<Long> computedReceivers = new HashSet<Long>()
    List<Coordinate> receivers = new ArrayList<>()
    fetchCellReceiver_withindex(connection, "receivers", receivers, computedReceivers)

    // receiver index
    @Field List<Long> receiversPk = new ArrayList<>()

    // init du vecteur de fréquences
    List<Integer> db_field_freq = [63, 125, 250, 500, 1000, 2000, 4000, 8000]

    // rose of favourable conditions
    // todo pour le moment pas pris en compte
    double[] favrose = [0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25]
    double[] favrose_6_18 = [0.41, 0.39, 0.38, 0.37, 0.36, 0.35, 0.34, 0.34, 0.34, 0.38, 0.44, 0.47, 0.48, 0.49, 0.49, 0.47, 0.45, 0.43];
    double[] favrose_18_22 = [0.53, 0.45, 0.41, 0.39, 0.37, 0.36, 0.36, 0.38, 0.41, 0.53, 0.62, 0.65, 0.67, 0.67, 0.68, 0.69, 0.68, 0.64];
    double[] favrose_22_6 = [0.64, 0.57, 0.51, 0.48, 0.46, 0.43, 0.41, 0.38, 0.36, 0.38, 0.44, 0.49, 0.52, 0.55, 0.58, 0.61, 0.64, 0.67];
    double[] favrose2 = [0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25];

    // ca c'est la table avec tous les rayons, attention gros espace memoire !
    HashMap<Integer,ComputeRaysOut> propaMap = new HashMap<>()

    // ----------------------------------
    // Et la on commence la boucle sur les simus
    // ----------------------------------
    System.out.println("Run Simu...")
    for (int r = 0; r < nSimu; ++r) {
        System.out.println(Math.round(100 * r / nSimu) + "%")
        rowNum = 0
        def timeStart = new Date()

        // Quand c'est le premier tout de simu, il y a des petites initialisations des tables postgresql nécessaires
        if (r == 0) {
            sql.execute(new File("../../../sql/LoadSrcsFromPGis.sql").text)
        }

        // Ici on rempli les tables pour permettre le calcul
        sql.execute("delete from receiver_lvl_day_zone;")
        sql.execute("delete from receiver_lvl_evening_zone;")
        sql.execute("delete from receiver_lvl_night_zone;")

        def timeStart2 = new Date()

        sql.execute("truncate roads_src_zone;")
        System.out.println("Compute Noise Emission...")
        def qry2 = 'INSERT INTO roads_src_zone (id, the_geom, ' +
                'db_m_d63, db_m_d125, db_m_d250, db_m_d500, db_m_d1000,db_m_d2000, db_m_d4000, db_m_d8000,' +
                'db_m_e63, db_m_e125, db_m_e250, db_m_e500, db_m_e1000,db_m_e2000, db_m_e4000, db_m_e8000,' +
                'db_m_n63, db_m_n125, db_m_n250, db_m_n500, db_m_n1000,db_m_n2000, db_m_n4000, db_m_n8000) ' +
                'VALUES (?,?, ' +
                '?, ?, ?, ?, ?,?, ?, ?,' +
                '?, ?, ?, ?, ?,?, ?, ?,' +
                '?, ?, ?, ?, ?,?, ?, ?)'
        sql.withBatch(100, qry2) { ps ->
            sql.eachRow('SELECT id, the_geom,\n' +
                    'lv_d_speed,mv_d_speed,hv_d_speed,wav_d_speed,wbv_d_speed,\n' +
                    'lv_e_speed,mv_e_speed,hv_e_speed,wav_e_speed,wbv_e_speed,\n' +
                    'lv_n_speed,mv_n_speed,hv_n_speed,wav_n_speed,wbv_n_speed,\n' +
                    'vl_d_per_hour,ml_d_per_hour,pl_d_per_hour,wa_d_per_hour,wb_d_per_hour,\n' +
                    'vl_e_per_hour,ml_e_per_hour,pl_e_per_hour,wa_e_per_hour,wb_e_per_hour,\n' +
                    'vl_n_per_hour,ml_n_per_hour,pl_n_per_hour,wa_n_per_hour,wb_n_per_hour,\n' +
                    'Zstart,Zend, Juncdist, Junc_type,road_pav FROM ROADS_TRAFFIC_ZONE_CAPTEUR_format2') { row ->

                //connect
                def list = [63, 125, 250, 500, 1000, 2000, 4000, 8000]
                def kk = 0
                def res_d = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
                def res_e = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
                def res_n = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

                int id = row[0]
                def the_geom = row[1]
                double lv_d_speed = row[2]
                double mv_d_speed = row[3]
                double hv_d_speed = row[4]
                double wav_d_speed = row[5]
                double wbv_d_speed = row[6]
                double lv_e_speed = row[7]
                double mv_e_speed = row[8]
                double hv_e_speed = row[9]
                double wav_e_speed = row[10]
                double wbv_e_speed = row[11]
                double lv_n_speed = row[12]
                double mv_n_speed = row[13]
                double hv_n_speed = row[14]
                double wav_n_speed = row[15]
                double wbv_n_speed = row[16]
                double vl_d_per_hour = row[17]
                double ml_d_per_hour = row[18]
                double pl_d_per_hour = row[19]
                double wa_d_per_hour = row[20]
                double wb_d_per_hour = row[21]
                double vl_e_per_hour = row[22]
                double ml_e_per_hour = row[23]
                double pl_e_per_hour = row[24]
                double wa_e_per_hour = row[25]
                double wb_e_per_hour = row[26]
                double vl_n_per_hour = row[27]
                double ml_n_per_hour = row[28]
                double pl_n_per_hour = row[29]
                double wa_n_per_hour = row[30]
                double wb_n_per_hour = row[31]
                double Zstart = row[32]
                double Zend = row[33]
                double Juncdist = row[34]
                int Junc_type = row[35]
                int road_pav = row[36]

                if (R_Junc[r] == 1) {
                    Juncdist = 200
                    Junc_dist[r] = 1
                }

                // todo corriger ça, maintenant c'est plus des nombres mais le nom du revêtement directement
                switch (R_Road[r]) {
                    case 1:
                        //15:FR1 17:FR2 19:FR3
                        if (road_pav >= 2 && road_pav <= 4) {
                            road_pav = 15
                        } else {
                            road_pav = 17
                        }
                        break
                    case 2:
                        road_pav = 17
                        break
                    case 3:
                        if (road_pav == 5 && Math.random() < 0.2) {
                            road_pav = 19
                        } else {
                            road_pav = 17
                        }
                        break
                }


                // Ici on calcule les valeurs d'emission par tronçons et par fréquence
                for (f in list) {
                    // fois 0.5 car moitié dans un sens et moitié dans l'autre
                    RSParametersCnossos srcParameters_d = new RSParametersCnossos(lv_d_speed * Speed_d_lv[r], mv_d_speed * Speed_d_mv[r], hv_d_speed * Speed_d_hv[r], wav_d_speed * Speed_d_wav[r], wbv_d_speed * Speed_d_wbv[r],
                            vl_d_per_hour * Deb_d_lv[r] * 0.5, ml_d_per_hour * Deb_d_mv[r] * 0.5, pl_d_per_hour * Deb_d_hv[r] * 0.5, wa_d_per_hour * Deb_d_wav[r] * 0.5, wb_d_per_hour * Deb_d_wbv[r] * 0.5,
                            f, Tempd_year[r], "NL01", 0, 0, Juncdist * Junc_dist[r], Junc_type)
                    RSParametersCnossos srcParameters_e = new RSParametersCnossos(lv_e_speed * Speed_e_lv[r], mv_e_speed * Speed_e_mv[r], hv_e_speed * Speed_e_hv[r], wav_e_speed * Speed_e_wav[r], wbv_e_speed * Speed_e_wbv[r],
                            vl_e_per_hour * Deb_e_lv[r] * 0.5, ml_e_per_hour * Deb_e_mv[r] * 0.5, pl_e_per_hour * Deb_e_hv[r] * 0.5, wa_e_per_hour * Deb_e_wav[r] * 0.5, wb_e_per_hour * Deb_e_wbv[r] * 0.5,
                            f, Tempe_year[r], "NL01", 0, 0, Juncdist * Junc_dist[r], Junc_type)
                    RSParametersCnossos srcParameters_n = new RSParametersCnossos(lv_n_speed * Speed_n_lv[r], mv_n_speed * Speed_n_mv[r], hv_n_speed * Speed_n_hv[r], wav_n_speed * Speed_n_wav[r], wbv_n_speed * Speed_n_wbv[r],
                            vl_n_per_hour * Deb_n_lv[r] * 0.5, ml_n_per_hour * Deb_n_mv[r] * 0.5, pl_n_per_hour * Deb_n_hv[r] * 0.5, wa_n_per_hour * Deb_n_wav[r] * 0.5, wb_n_per_hour * Deb_n_wbv[r] * 0.5,
                            f, Tempn_year[r], "NL01", 0, 0, Juncdist * Junc_dist[r], Junc_type)

                    srcParameters_d.setSlopePercentage(RSParametersCnossos.computeSlope(Zstart, Zend, the_geom.getLength()))
                    srcParameters_e.setSlopePercentage(RSParametersCnossos.computeSlope(Zstart, Zend, the_geom.getLength()))
                    srcParameters_n.setSlopePercentage(RSParametersCnossos.computeSlope(Zstart, Zend, the_geom.getLength()))
                    res_d[kk] = EvaluateRoadSourceCnossos.evaluate(srcParameters_d)
                    res_e[kk] = EvaluateRoadSourceCnossos.evaluate(srcParameters_e)
                    res_n[kk] = EvaluateRoadSourceCnossos.evaluate(srcParameters_n)
                    srcParameters_d.setSlopePercentage(RSParametersCnossos.computeSlope(Zend, Zstart, the_geom.getLength()))
                    srcParameters_e.setSlopePercentage(RSParametersCnossos.computeSlope(Zend, Zstart, the_geom.getLength()))
                    srcParameters_n.setSlopePercentage(RSParametersCnossos.computeSlope(Zend, Zstart, the_geom.getLength()))
                    res_d[kk] = 10 * Math.log10(Math.pow(10, EvaluateRoadSourceCnossos.evaluate(srcParameters_d) / 10) + Math.pow(10, res_d[kk] / 10))
                    res_e[kk] = 10 * Math.log10(Math.pow(10, EvaluateRoadSourceCnossos.evaluate(srcParameters_e) / 10) + Math.pow(10, res_e[kk] / 10))
                    res_n[kk] = 10 * Math.log10(Math.pow(10, EvaluateRoadSourceCnossos.evaluate(srcParameters_n) / 10) + Math.pow(10, res_n[kk] / 10))


                    kk++
                }

                //def timeStart3 = System.currentTimeMillis()
                // On rempli la table correspondante des valeurs des emissions en dB

                ps.addBatch(id, the_geom,
                        res_d[0], res_d[1], res_d[2], res_d[3], res_d[4], res_d[5], res_d[6], res_d[7],
                        res_e[0], res_e[1], res_e[2], res_e[3], res_e[4], res_e[5], res_e[6], res_e[7],
                        res_n[0], res_n[1], res_n[2], res_n[3], res_n[4], res_n[5], res_n[6], res_n[7])

                /*def params = [id, the_geom,
                              res_d[0], res_d[1], res_d[2], res_d[3], res_d[4], res_d[5], res_d[6], res_d[7],
                              res_e[0], res_e[1], res_e[2], res_e[3], res_e[4], res_e[5], res_e[6], res_e[7],
                              res_n[0], res_n[1], res_n[2], res_n[3], res_n[4], res_n[5], res_n[6], res_n[7]]
                sql.execute(qry, params)*/
                //TimeE2 = System.currentTimeMillis() - timeStart3

            }
        }
        TimeE = TimeCategory.minus(new Date(), timeStart2).toMilliseconds()

        // Ici on importe les sources
        // -> on est d'accord que ecrire dans
        // la base de données pour pouvoir lire ça sert pas à grand chose
        // et ça prend beaucoup de temps de calcul, donc c'est optimisable
        System.out.println("Load Sources in Pgis...")
        fetchCellSource_withindex(connection, null, sourceGeometries, "roads_src_zone", sourcesPk, wj_sources, sourcesIndex);

        // Je crois que ça ça sert pas à grand chose mais il faudra vérifier
        sql.execute(new File("../../../sql/Lorient_LoadTables3Pgis.sql").text)


        if (r == 0) {
            System.out.println("Compute Rays...")
            // Configure noisemap with specified receivers
            timeStart2 = new Date()
            //-----------------------------------------------------------------
            // ----------- ICI On calcul les rayons entre sources et recepteurs (c est pour ça r=0, on le fait qu'une fois)
            //-----------------------------------------------------------------

            // on a pas vraiment besoin de energeticsum, enfin on en discutera
            double[] energeticSum = new double[8]

            // ici on configure les data pour tirer les rayons
            PropagationProcessData rayData = new PropagationProcessData(new ArrayList<Coordinate>(), manager, sourcesIndex, sourceGeometries, wj_sources, db_field_freq, 1, 0, 500, 500, 1.0, 0.2, favrose2, 0.1, 0, null, geoWithSoilTypeList, true)
            // phase d'init
            rayData.makeRelativeZToAbsoluteOnlySources()


            // et la pour chacun des receptuers, on va chercher les rayons vers les sources
            for (int pk = 0; pk < receivers.size(); pk++) {
                System.out.println(100*pk/receivers.size()+ " %")
                ComputeRaysOut propDataOut = new ComputeRaysOut()
                ComputeRays propaProcess = new ComputeRays(rayData, propDataOut)
                propaProcess.initStructures()

                receivers.get(pk).setCoordinate(new Coordinate(receivers.get(pk).x ,receivers.get(pk).y ,receivers.get(pk).z + manager.getHeightAtPosition(receivers.get(pk))))
                // Et ici on calcule les rayons
                propaProcess.computeRaysAtPosition(receivers.get(pk), receiversPk.get(pk).toInteger(), energeticSum, null)

                // on stocke dans propaMap
                if (!propDataOut.getPropagationPaths().isEmpty()) {
                    propaMap.put(pk, propDataOut)
                    // Ca c'est pour voir les rayons, je met en commentaire mais ça peut s'activer
                    String filename = "D:/aumond/Desktop/Results/Incertitudes" + pk.toString() + ".vtk";
                    try {
                        writeVTK(filename, propDataOut)
                    } catch (IOException e) {
                        e.printStackTrace()
                    }
                }
            }
            TimeR = TimeCategory.minus(new Date(), timeStart2).toMilliseconds()


        }

        timeStart2 = new Date()
        // Ici on rentre dans la phase calcul de la matrice de transfer
        ComputeRaysOut output = new ComputeRaysOut()
        System.out.println("Compute Attenuation...")

        Iterator it = propaMap.entrySet().iterator()
        while (it.hasNext()) {
            Map.Entry pair = (Map.Entry)it.next()
            int idReceiver = pair.getKey()

            HashMap<Integer, double[]> aGlobal = new HashMap<>()

            //todo ça c'est des parametres a faire bouger (parmi d'autres)
            double p = 0.5 // probability favourable conditions
            PropagationProcessPathData propData = new PropagationProcessPathData()
            propData.setTemperature(10)
            propData.setHumidity(70)
            propData.setPrime2520(true)

            // et c'est cette ligne le coeur de calcul
            computeWithMeteo(aGlobal, propData, propaMap.get(idReceiver), p)

            // ca c'est pour ecrire dans la table postgre
            // puisque la combinaison sources + transfer matrix
            // se fait encore en postgis
            // mais dans le futur pourquoi pas rester en groovy
            // ca aura plus de sens
            Iterator it2 = aGlobal.entrySet().iterator()
            def qry = 'INSERT INTO RECEIVER_LVL_DAY_ZONE (IDRECEPTEUR, IDSOURCE,' +
                    'ATT63, ATT125, ATT250, ATT500, ATT1000,ATT2000, ATT4000, ATT8000) ' +
                    'VALUES (?,?,?,?,?,?,?,?,?,?);'
            sql.withBatch(100, qry) { ps ->
                while (it2.hasNext()) {
                    Map.Entry pair2 = (Map.Entry) it2.next()
                    double[] att = (double[]) pair2.getValue()
                    int idRecv = receiversPk.get(idReceiver).toInteger()
                    int idSrc = sourcesPk.get((Integer) pair2.getKey()).toInteger()
                    output.addVerticeSoundLevel(idRecv, (Integer) pair2.getKey(), (double[]) pair2.getValue())




                    ps.addBatch(idRecv, idSrc,
                            att[0], att[1], att[2], att[3], att[4], att[5], att[6], att[7])


                    it2.remove() // avoids a ConcurrentModificationException
                }
                //ps.executeBatch()
            }






            it.remove()
        }




        TimeA = TimeCategory.minus(new Date(), timeStart2).toMilliseconds()

        // A partir de la jusqua la fin c'est de la jointure de table,
        // pour faire emetteur + matrice de trasnfer + recepteurs
        // et de l'ecrtirue des fichiers de sorties

        // todo compute for evening and night
        sql.execute '''drop table if exists receiver_lvl_evening_zone, receiver_lvl_night_zone;'''
        sql.execute '''create table receiver_lvl_evening_zone as select * from receiver_lvl_day_zone;'''
        sql.execute '''create table receiver_lvl_night_zone as select * from receiver_lvl_day_zone;'''

        csv.exportTable(connection, "receiver_lvl_day_zone", new File(workspace_output + "/data/receiver_lvl_day_zone.csv"), new EmptyProgressVisitor())
        csv.exportTable(connection, "receiver_lvl_evening_zone", new File(workspace_output + "/data/receiver_lvl_evening_zone.csv"), new EmptyProgressVisitor())
        csv.exportTable(connection, "receiver_lvl_night_zone", new File(workspace_output + "/data/receiver_lvl_night_zone.csv"), new EmptyProgressVisitor())

         if (r == 0) {
            sql.execute(new File("../../../sql/Reception_PrimaryPgis.sql").text)
        }

        sql.execute(new File("../../../sql/ReceptionPgis.sql").text)
        switch (batrec) {
            case 1:
                csvFile = new File(workspace_output + "/data/Simu_" + r + "_pop_lvl_n.csv")
                exp.exportTable(connection, "pop_lvl_n", csvFile, new EmptyProgressVisitor())
                csvFile = new File(workspace_output + "/data/Simu_" + r + "_pop_lvl_den.csv")
                exp.exportTable(connection, "pop_lvl_den", csvFile, new EmptyProgressVisitor())
                break
            case 2:
                csvFile = new File(workspace_output + "/data/Simu_" + r + "_Lden.csv")
                exp.exportTable(connection, "Lden", csvFile, new EmptyProgressVisitor())
                sql.execute '''drop table if exists receivers_lden_zone;'''
                sql.execute '''create table receivers_lden_zone as 
                                        select cast(s.db as float) as Lden, r.the_geom as the_geom 
                                        from lvl_receiver_lvl_day_zone as s, receivers as r 
                                        where s.id = r.id;'''
                break
            case 3:
                csvFile = new File(workspace_output + "/data/Simu_" + r + "_pop_lvl_n.csv")
                exp.exportTable(connection, "pop_lvl_n", csvFile, new EmptyProgressVisitor())
                csvFile = new File(workspace_output + "/data/Simu_" + r + "_pop_lvl_den.csv")
                exp.exportTable(connection, "pop_lvl_den", csvFile, new EmptyProgressVisitor())
                csvFile = new File(workspace_output + "/data/Simu_" + r + "_Lden.csv")
                exp.exportTable(connection, "Lden", csvFile, new EmptyProgressVisitor())
                break
        }

        def timeStop = new Date()
        CalcTime[r] = TimeCategory.minus(timeStop, timeStart).toMilliseconds()
        sb2.append(Simu[r].toString() + "\t")
        sb2.append(CalcTime[r].toString() + "\n")

    }

    CalcTime_record.write(sb2.toString())
    CalcTime_record.close()

    //la on zip les results
    String zipFileName = "results_C" + n_comp.toString() + ".zip"
    String inputDir = workspace_output

    ZipOutputStream output = new ZipOutputStream(new FileOutputStream(workspace_output + zipFileName))

    new File(inputDir).eachFile() { file ->
        if (!file.isFile()) {
            return
        }
        println file.name.toString()
        println file.toString()

        output.putNextEntry(new ZipEntry(file.name.toString())) // Create the name of the entry in the ZIP

        InputStream input = new FileInputStream(file)

        // Stream the document data to the ZIP
        Files.copy(input, output)
        output.closeEntry() // End of current document in ZIP
        input.close()
    }
    output.close() // End of all documents - ZIP is complete


} finally {

    sql.close()
    connection.close()
}
