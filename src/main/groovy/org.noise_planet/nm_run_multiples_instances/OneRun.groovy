package org.noise_planet.nm_run_multiples_instances

import com.esotericsoftware.kryo.Kryo
import com.esotericsoftware.kryo.io.Input
import com.esotericsoftware.kryo.io.Output
import com.esotericsoftware.kryo.serializers.MapSerializer

/** author : Aumond Pierre

 Ce code permet de faire une unique de simulation

 Le résultat se trouve dans un dossier zippé
 Il s'agit d'une version Postgre/Postgis (lecture/ecriture plus lent qu'avec h2gis mais jointures beaucoup plus rapides)

 **/

// Importation des librairies nécessaire au code
import groovy.sql.Sql

import org.apache.commons.io.FileUtils
import org.cts.CRSFactory
import org.cts.crs.CoordinateReferenceSystem
import org.cts.crs.GeodeticCRS
import org.cts.registry.EPSGRegistry
import org.cts.registry.RegistryManager
import org.h2gis.functions.spatial.crs.EPSGTuple
import org.h2gis.functions.spatial.crs.ST_Transform
import org.h2gis.functions.spatial.volume.GeometryExtrude
import org.h2gis.utilities.JDBCUtilities
import org.h2gis.utilities.SFSUtilities
import org.h2gis.utilities.SpatialResultSet
import org.h2gis.utilities.TableLocation
import org.locationtech.jts.algorithm.Angle
import org.locationtech.jts.geom.Coordinate
import org.locationtech.jts.geom.Envelope
import org.locationtech.jts.geom.Geometry
import org.locationtech.jts.geom.GeometryCollection
import org.locationtech.jts.geom.GeometryFactory
import org.locationtech.jts.geom.MultiPolygon
import org.locationtech.jts.geom.Point
import org.locationtech.jts.geom.Polygon
import org.locationtech.jts.geom.impl.CoordinateArraySequence
import org.locationtech.jts.geom.util.GeometryEditor
import org.orbisgis.noisemap.core.ComputeRays
import org.orbisgis.noisemap.core.ComputeRaysOut
import org.orbisgis.noisemap.core.EvaluateAttenuationCnossos
import org.orbisgis.noisemap.core.FastObstructionTest
import org.orbisgis.noisemap.core.GeoWithSoilType
import org.orbisgis.noisemap.core.LayerDelaunayError
import org.orbisgis.noisemap.core.MeshBuilder
import org.orbisgis.noisemap.core.PointsMerge
import org.orbisgis.noisemap.core.PropagationPath
import org.orbisgis.noisemap.core.PropagationProcessData
import org.orbisgis.noisemap.core.PropagationProcessPathData
import org.orbisgis.noisemap.core.QueryGeometryStructure
import org.orbisgis.noisemap.core.QueryQuadTree
import org.h2gis.functions.io.csv.CSVDriverFunction
import org.h2gis.utilities.wrapper.ConnectionWrapper
import org.cts.op.CoordinateOperationFactory
import java.sql.Connection
import java.sql.DriverManager
import java.sql.PreparedStatement
import java.sql.SQLException

import org.orbisgis.noisemap.core.EvaluateRoadSourceCnossos
import org.orbisgis.noisemap.core.RSParametersCnossos

import groovy.time.*
import java.nio.file.Files
import java.util.zip.ZipEntry
import java.util.zip.ZipOutputStream
import org.h2gis.api.EmptyProgressVisitor

import groovy.transform.SourceURI
import java.nio.file.Path
import java.nio.file.Paths




class OneRun {
    static void main(String[] args) {
        OneRun oneRun = new OneRun()
        oneRun.run()
    }

    void run() {
/**
 ///////////////////////////////////////////
 // Paramètres d'entrée et initialisations //
 ///////////////////////////////////////////

 */
        String workspace_output = "D:\\output"
        // le workspace output doit respecter une arborescence specifique
        //output
        //|- data
        //|- config

        boolean clean = true // reset the postgis database
        boolean printply = true // sort un ply de la zone
        boolean saveRays = true
        boolean loadRays = false
        int batrec = 2
        // Localisation des récepteurs :  (1) : bâtiments, (2) : récepteurs, (3) : bâtiments + récepteurs

        boolean H2 = true

        String h2_url
        String user_name
        String user_password
        if (!H2) {
            h2_url = "jdbc:postgresql_h2://127.0.0.1:5433/noisemodelling"
            user_name = "postgres"
            user_password = "12345678"
        } else {
            // Ouverture de la base de données H2 qui recueillera les infos
            h2_url = "jdbc:h2:/D:/db/database;DB_CLOSE_DELAY=30;DEFRAG_ALWAYS=TRUE"
            user_name = "sa"
            user_password = ""
        }

        String zip_filename_root = "results_C"      // Racine du nom de l'archive des résultats (.zip)

        // Noms des tables en entrée et des attributs
        String zone_name = "zone_cense_2km"
        String buildings_table_name = "buildings_zone"
        String alpha_field = ""
        String height_field_name = "hauteur"
        String g_field_name = "g"
        String receivers_table_name = "receivers2"
        String sources_table_name = "roads_src_zone"
        String sources_table_name2 = "ROADS_TRAFFIC_ZONE_CAPTEUR_format2"
        String soil_table_name = "land_use_zone_capteur"
        String topo_table_name = "dem_lite2"

        @SourceURI
        URI sourceUri
        Path scriptLocation = Paths.get(sourceUri)
        String rootPath = scriptLocation.getParent().getParent().getParent().getParent().toString()

        // Paramètres de propagation
        int reflexion_order = 1
        int diffraction_order = 10
        double max_src_dist = 500
        double max_ref_dist = 500
        double min_ref_dist = 1.0
        double wall_alpha = 0.1 // todo pour le moment cette valeur ne peut pas être changé
        double forget_source = 0.1 // todo pour le moment cette valeur est inutile
        boolean compute_vertical_diffraction = true

        // rose of favourable conditions
        double[] favrose = [0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25]
        double[] favrose_6_18 = [0.41, 0.39, 0.38, 0.37, 0.36, 0.35, 0.34, 0.34, 0.34, 0.38, 0.44, 0.47, 0.48, 0.49, 0.49, 0.47, 0.45, 0.43]
        double[] favrose_18_22 = [0.53, 0.45, 0.41, 0.39, 0.37, 0.36, 0.36, 0.38, 0.41, 0.53, 0.62, 0.65, 0.67, 0.67, 0.68, 0.69, 0.68, 0.64]
        double[] favrose_22_6 = [0.64, 0.57, 0.51, 0.48, 0.46, 0.43, 0.41, 0.38, 0.36, 0.38, 0.44, 0.49, 0.52, 0.55, 0.58, 0.61, 0.64, 0.67]
        double[] favrose2 = [0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25]

/**
 ///////////////////////////////////////////
 // FIN Paramètres d'entrée et initialisations //
 ///////////////////////////////////////////

 */

// Ouverture de la base de donnees PostGreSQL qui recueillera les infos
        ArrayList<String> urls = new ArrayList<>()
        ArrayList<Connection> connections = new ArrayList<>()

        // Connexion à la base
        def driver
        if (!H2) {
            driver = 'org.orbisgis.postgis_jts.Driver'
        } else {
            driver = 'org.h2.Driver'
        }
        Class.forName(driver)
        Connection connection = new ConnectionWrapper(DriverManager.getConnection(h2_url, user_name, user_password))
        connections.add(connection)
        def sql = Sql.newInstance(h2_url, user_name, user_password, driver)

        // Dossier des résultats
        File newFile = new File(workspace_output + "/config/")
        FileUtils.cleanDirectory(newFile)
        newFile = new File(workspace_output + "/data/")
        FileUtils.cleanDirectory(newFile)

        try {
            System.out.println("Init")
            CSVDriverFunction csv = new CSVDriverFunction()

            CSVDriverFunction exp = new CSVDriverFunction()
            //////////////////////
            // Import file text
            //////////////////////

            ArrayList<Integer> Simu = new ArrayList<Integer>()
            ArrayList<Double> CalcTime = new ArrayList<Double>()

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
                sql.execute(new File(rootPath + "/sql/Cleaning_database.sql").text)
                // import des shapefiles ici
                sql.execute('drop table ZONE_CENSE_2KM if exists ')
                sql.execute('drop table buildings_zone if exists ')
                sql.execute('drop table receivers2 if exists ')
                sql.execute('drop table land_use_zone_capteur if exists ')
                sql.execute('drop table dem_lite2 if exists ')
                sql.execute('drop table roads_src_zone if exists ')
                sql.execute('drop table ROADS_TRAFFIC_ZONE_CAPTEUR_format2 if exists ')

                sql.execute("CALL File_table('D:/Documents/Stage/LETG-IFSTTAR/IFSTTAR/couches_clean/emprise.shp','zone_cense_2km')")
                sql.execute("CALL File_table('D:/Documents/Stage/LETG-IFSTTAR/IFSTTAR/couches_clean/buildings.shp','buildings_zone')")
                sql.execute("CALL File_table('D:/Documents/Stage/LETG-IFSTTAR/IFSTTAR/couches_clean/recepteurs1.shp','receivers2')")
                sql.execute("CALL File_table('D:/Documents/Stage/LETG-IFSTTAR/IFSTTAR/couches_clean/ground_g2.shp','land_use_zone_capteur')")
                sql.execute("CALL File_table('D:/Documents/Stage/LETG-IFSTTAR/IFSTTAR/couches_clean/MNT1_lite.shp','dem_lite2')")
                sql.execute("CALL File_table('D:/Documents/Stage/LETG-IFSTTAR/IFSTTAR/couches_clean/vide.shp','roads_src_zone')")
                sql.execute("CALL File_table('D:/Documents/Stage/LETG-IFSTTAR/IFSTTAR/couches_clean/routes_adaptees.shp','ROADS_TRAFFIC_ZONE_CAPTEUR_format2')")

                sql.execute("create spatial index on zone_cense_2km(the_geom)")
                sql.execute("create spatial index on buildings_zone(the_geom)")
                sql.execute("create spatial index on receivers2(the_geom)")
                sql.execute("create spatial index on land_use_zone_capteur(the_geom)")
                sql.execute("create spatial index on dem_lite2(the_geom)")
                sql.execute("create spatial index on roads_src_zone(the_geom)")
                sql.execute("create spatial index on ROADS_TRAFFIC_ZONE_CAPTEUR_format2(the_geom)")

                // Préparation de la zone d'étude

                //sql.execute("drop table if exists zone;")
                //def zone_sql_select = "create table zone as select * from " + zone_name + ";"
                //sql.execute(zone_sql_select)
                //sql.execute("alter table zone rename column geom to the_geom;")
                //sql.execute("alter table zone alter column the_geom type geometry using ST_SetSRID(the_geom, 2154);")
                //sql.execute("create index zone_the_geom_gist on zone using GIST (the_geom);")
                // Préparation des tables d'entée types batiments, routes, etc., sur la zone (Need postgis extension sfcgal)
                //sql.execute(new File(rootPath + "/sql/LoadTables2Pgis2.sql").text)
            }
            // ------------------------------------------------------------ //
            // ----------- Initialisation et import des données ----------- //
            // -----------        (sol et des bâtiments)        ----------- //
            // ------------------------------------------------------------ //

            System.out.println("Import Tables...")
            // Import GeometrySoilType
            List<GeoWithSoilType> geoWithSoilTypeList = new ArrayList<>()
            System.out.println("LandUse...")
            fetchCellSoilAreas(connection, geoWithSoilTypeList, soil_table_name, g_field_name)

            MeshBuilder mesh = new MeshBuilder()

            // Import Buildings
            System.out.println("Buildings...")
            fetchCellBuildings(connection, mesh, buildings_table_name, alpha_field, height_field_name)

            // Enveloppe de la scène
            Envelope cellEnvelope = mesh.getEnvelope()

            // Importer le Sol
            System.out.println("Altitudes...")
            //sql.execute("alter table full_dem_lite2 rename column geom to the_geom;")
            fetchCellDem(connection, mesh, topo_table_name, "contour")

            mesh.finishPolygonFeeding(cellEnvelope)
            FastObstructionTest manager = new FastObstructionTest(mesh.getPolygonWithHeight(), mesh.getTriangles(),
                    mesh.getTriNeighbors(), mesh.getVertices())

            // Ecriture de la topo dans un fichier polygon ouvrable avec paraview par exemple
            if (printply) {
                System.out.println("Impression en cours...")
                String filename2 = workspace_output + "\\" + zone_name + ".ply"
                String filename3 = workspace_output + "\\" + zone_name + ".kml"
                try {
                   // writePLY(filename2, mesh)
                    writeKML(filename3, mesh, manager)
                } catch (IOException e) {
                    e.printStackTrace()
                }
            }

            // import sources
            List<Long> sourcesPk = new ArrayList<>()
            ArrayList<Geometry> sourceGeometries = new ArrayList<>()
            ArrayList<ArrayList<Double>> wj_sources = new ArrayList<>()
            QueryGeometryStructure sourcesIndex = new QueryQuadTree()

            System.out.println("Import receivers...")
            // receiver index
            List<Long> receiversPk = new ArrayList<>()
            Set<Long> computedReceivers = new HashSet<Long>()
            List<Coordinate> receivers = new ArrayList<>()
            receiversPk = fetchCellReceiver_withindex(connection, receivers_table_name, receivers, computedReceivers)

            // init du vecteur de fréquences
            List<Integer> db_field_freq = [63, 125, 250, 500, 1000, 2000, 4000, 8000]

            // ca c'est la table avec tous les rayons, attention gros espace memoire !
            HashMap<Integer, ComputeRaysOut> propaMap = new HashMap<>()

            // ----------------------------------
            // Et la on commence la boucle sur les simus
            // ----------------------------------
            System.out.println("Run Simu...")

            def timeStart = new Date()

            sql.execute(new File(rootPath + "/sql/LoadSrcsFromPGis.sql").text)

            // Ici on rempli les tables pour permettre le calcul
            sql.execute("delete from receiver_lvl_day_zone;")
            sql.execute("delete from receiver_lvl_evening_zone;")
            sql.execute("delete from receiver_lvl_night_zone;")

            def timeStart2 = new Date()
            if (!H2) {
                sql.execute("truncate " + sources_table_name + ";")
            } else {
                sql.execute("truncate table " + sources_table_name + ";")
            }

            System.out.println("Compute Noise Emission...")
            def qry2 = 'INSERT INTO ' + sources_table_name + '(id, the_geom, ' +
                    'db_m_d63, db_m_d125, db_m_d250, db_m_d500, db_m_d1000,db_m_d2000, db_m_d4000, db_m_d8000,' +
                    'db_m_e63, db_m_e125, db_m_e250, db_m_e500, db_m_e1000,db_m_e2000, db_m_e4000, db_m_e8000,' +
                    'db_m_n63, db_m_n125, db_m_n250, db_m_n500, db_m_n1000,db_m_n2000, db_m_n4000, db_m_n8000) ' +
                    'VALUES (?,?, ' +
                    '?, ?, ?, ?, ?,?, ?, ?,' +
                    '?, ?, ?, ?, ?,?, ?, ?,' +
                    '?, ?, ?, ?, ?,?, ?, ?)'
            sql.withBatch(100, qry2) { ps ->
                // memes vleurs d e et n
                sql.eachRow('SELECT id, the_geom,\n' +
                        'lv_d_speed,mv_d_speed,hv_d_speed,wav_d_spee,wbv_d_spee,\n' +
                        'lv_e_speed,mv_e_speed,hv_e_speed,wav_e_spee,wbv_e_spee,\n' +
                        'lv_n_speed,mv_n_speed,hv_n_speed,wav_n_spee,wbv_n_spee,\n' +
                        'vl_d_per_h,ml_d_per_h,pl_d_per_h,wa_d_per_h,wb_d_per_h,\n' +
                        'vl_e_per_h,ml_e_per_h,pl_e_per_h,wa_e_per_h,wb_e_per_h,\n' +
                        'vl_n_per_h,ml_n_per_h,pl_n_per_h,wa_n_per_h,wb_n_per_h,\n' +
                        'Zstart,Zend, Juncdist, Junc_type,road_pav FROM ' + sources_table_name2 + ' ;') { row ->

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

                    // Ici on calcule les valeurs d'emission par tronçons et par fréquence
                    for (f in list) {
                        // fois 0.5 car moitié dans un sens et moitié dans l'autre
                        RSParametersCnossos srcParameters_d = new RSParametersCnossos(lv_d_speed, mv_d_speed, hv_d_speed, wav_d_speed, wbv_d_speed,
                                vl_d_per_hour * 0.5, ml_d_per_hour * 0.5, pl_d_per_hour * 0.5, wa_d_per_hour * 0.5, wb_d_per_hour * 0.5,
                                f, 20, "NL01", 0, 0, Juncdist, Junc_type)
                        RSParametersCnossos srcParameters_e = new RSParametersCnossos(lv_e_speed, mv_e_speed, hv_e_speed, wav_e_speed, wbv_e_speed,
                                vl_e_per_hour * 0.5, ml_e_per_hour * 0.5, pl_e_per_hour * 0.5, wa_e_per_hour * 0.5, wb_e_per_hour * 0.5,
                                f, 20, "NL01", 0, 0, 200, Junc_type)
                        RSParametersCnossos srcParameters_n = new RSParametersCnossos(lv_n_speed, mv_n_speed, hv_n_speed, wav_n_speed, wbv_n_speed,
                                vl_n_per_hour * 0.5, ml_n_per_hour * 0.5, pl_n_per_hour * 0.5, wa_n_per_hour * 0.5, wb_n_per_hour * 0.5,
                                f, 20, "NL01", 0, 0, 200, Junc_type)

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

                    // On rempli la table correspondante des valeurs des emissions en dB
                    ps.addBatch(id, the_geom,
                            res_d[0], res_d[1], res_d[2], res_d[3], res_d[4], res_d[5], res_d[6], res_d[7],
                            res_e[0], res_e[1], res_e[2], res_e[3], res_e[4], res_e[5], res_e[6], res_e[7],
                            res_n[0], res_n[1], res_n[2], res_n[3], res_n[4], res_n[5], res_n[6], res_n[7])

                }
            }

            sql.execute('ALTER TABLE ' + sources_table_name + ' ALTER COLUMN ID SET NOT NULL')
            sql.execute('ALTER TABLE ' + sources_table_name + ' ADD PRIMARY KEY (ID)')
            // Ici on importe les sources
            // -> on est d'accord que ecrire dans
            // la base de données pour pouvoir lire ça sert pas à grand chose
            // et ça prend beaucoup de temps de calcul, donc c'est optimisable
            System.out.println("Load Sources in Pgis...")

            fetchCellSource_withindex(connection, null, sourceGeometries, sources_table_name, sourcesPk, wj_sources, sourcesIndex)

            if (!loadRays){
                System.out.println("Compute Rays...")
                // Configure noisemap with specified receivers
                //-----------------------------------------------------------------
                // ----------- ICI On calcul les rayons entre sources et recepteurs (c est pour ça r=0, on le fait qu'une fois)
                //-----------------------------------------------------------------

                // on a pas vraiment besoin de energeticsum, enfin on en discutera
                double[] energeticSum = new double[8]

                // ici on configure les data pour tirer les rayons
                PropagationProcessData rayData = new PropagationProcessData(new ArrayList<Coordinate>(), manager, sourcesIndex, sourceGeometries, wj_sources, db_field_freq, reflexion_order, diffraction_order, max_src_dist, max_ref_dist, min_ref_dist, wall_alpha, favrose2, forget_source, 0, null, geoWithSoilTypeList, true)
                // phase d'init
                rayData.makeRelativeZToAbsoluteOnlySources()

                // et la pour chacun des receptuers, on va chercher les rayons vers les sources
                //
                for (int pk = 0; pk < receivers.size(); pk++) {
                    System.out.println(100 * pk / receivers.size() + " %")
                    ComputeRaysOut propDataOut = new ComputeRaysOut()
                    ComputeRays propaProcess = new ComputeRays(rayData, propDataOut)
                    propaProcess.initStructures()

                    receivers.get(pk).setCoordinate(new Coordinate(receivers.get(pk).x, receivers.get(pk).y, receivers.get(pk).z + manager.getHeightAtPosition(receivers.get(pk))))
                    // Et ici on calcule les rayons
                    propaProcess.computeRaysAtPosition(receivers.get(pk), receiversPk.get(pk).toInteger(), energeticSum, null)

                    // on stocke dans propaMap
                    if (!propDataOut.getPropagationPaths().isEmpty()) {
                        propaMap.put(pk, propDataOut)
                        // Ca c'est pour voir les rayons, je met en commentaire mais ça peut s'activer
                        String filename = "D:/results/coucou" + pk.toString() + ".vtk"
                        String filenameKML = "D:/results/coucou" + pk.toString() + ".kml"
                        try {
                            //writeVTK(filename, propDataOut)
                            writeKML_Line(filenameKML, propDataOut, manager)
                        } catch (IOException e) {
                            e.printStackTrace()
                        }
                    }
                   /* if (propaMap.size()>10){
                        break
                    }*/
                }
            }


            // Ici on rentre dans la phase calcul de la matrice de transfer
            ComputeRaysOut output = new ComputeRaysOut()
            System.out.println("Compute Attenuation...")

            Kryo kryo = new Kryo()
            MapSerializer serializer = new MapSerializer()
            kryo.register(HashMap.class, serializer)
            kryo.register(LinkedHashMap.class, serializer)
            serializer.setKeyClass(String.class, kryo.getSerializer(String.class))
            serializer.setKeysCanBeNull(false)
            serializer.setKeyClass(String.class, kryo.getSerializer(String.class))

           String filenamebin = workspace_output + "\\Rays.bin"
          /*  if (saveRays && !loadRays){
                Output outputBin = new Output(new FileOutputStream(filenamebin))
                kryo.writeObject(outputBin, propaMap)
                outputBin.close()
            }

            if (loadRays){
                Input input = new Input(new FileInputStream(filenamebin))
                propaMap = kryo.readObject(input, HashMap.class)
                input.close()
            }
*/

            Iterator it = propaMap.entrySet().iterator()
            while (it.hasNext()) {
                Map.Entry pair = (Map.Entry) it.next()
                int idReceiver = pair.getKey()

                HashMap<Integer, double[]> aGlobal = new HashMap<>()

                PropagationProcessPathData propData = new PropagationProcessPathData()
                propData.setTemperature(10)
                propData.setHumidity(70)
                propData.setPrime2520(true)

                // et c'est cette ligne le coeur de calcul
                computeWithMeteo(aGlobal, propData, propaMap.get(idReceiver), favrose_22_6)

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
                }

                it.remove()
            }

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

            sql.execute(new File(rootPath + "/sql/Reception_Primary.sql").text)
            sql.execute(new File(rootPath + "/sql/Reception.sql").text)
            switch (batrec) {
                case 1:
                    File csvFile = new File(workspace_output + "/data/Simu_" + 0 + "_pop_lvl_n.csv")
                    exp.exportTable(connection, "pop_lvl_n", csvFile, new EmptyProgressVisitor())
                    csvFile = new File(workspace_output + "/data/Simu_" + 0 + "_pop_lvl_den.csv")
                    exp.exportTable(connection, "pop_lvl_den", csvFile, new EmptyProgressVisitor())
                    break
                case 2:
                    File csvFile = new File(workspace_output + "/data/Simu_" + 0 + "_Lden.csv")
                    exp.exportTable(connection, "Lden", csvFile, new EmptyProgressVisitor())
                    sql.execute '''drop table if exists receivers_lden_zone;'''
                    sql.execute '''create table receivers_lden_zone as 
                                        select cast(s.db as float) as Lden, r.the_geom as the_geom 
                                        from lvl_receiver_lvl_day_zone as s, receivers as r 
                                        where s.id = r.id;'''
                    break
                case 3:
                    File csvFile = new File(workspace_output + "/data/Simu_" + 0 + "_pop_lvl_n.csv")
                    exp.exportTable(connection, "pop_lvl_n", csvFile, new EmptyProgressVisitor())
                    csvFile = new File(workspace_output + "/data/Simu_" + 0 + "_pop_lvl_den.csv")
                    exp.exportTable(connection, "pop_lvl_den", csvFile, new EmptyProgressVisitor())
                    csvFile = new File(workspace_output + "/data/Simu_" + 0 + "_Lden.csv")
                    exp.exportTable(connection, "Lden", csvFile, new EmptyProgressVisitor())
                    break
            }

            def timeStop = new Date()
            CalcTime[0] = TimeCategory.minus(timeStop, timeStart).toMilliseconds()
            sb2.append(Simu[0].toString() + "\t")
            sb2.append(CalcTime[0].toString() + "\n")

            CalcTime_record.write(sb2.toString())
            CalcTime_record.close()

            //la on zip les results
            String zipFileName = "results.zip"
            String inputDir = workspace_output

            ZipOutputStream output_ = new ZipOutputStream(new FileOutputStream(workspace_output + zipFileName))

            new File(inputDir).eachFile() { file ->
                if (!file.isFile()) {
                    return
                }
                println file.name.toString()
                println file.toString()

                output_.putNextEntry(new ZipEntry(file.name.toString())) // Create the name of the entry in the ZIP

                InputStream input = new FileInputStream(file)

                // Stream the document data to the ZIP
                Files.copy(input, output_)
                output_.closeEntry() // End of current document in ZIP
                input.close()
            }
            output_.close() // End of all documents - ZIP is complete


        } finally {

            sql.close()
            connection.close()
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
            propath.setFavorable(false)
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


                fileWriter.write("\t\t\t"+coordinate.x.toString()+","+ coordinate.y.toString()+","+ (p.coordinate.z - z).toString() +"\n")

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

                fileWriter.write("\t\t\t"+coordinate2.x.toString()+","+ coordinate2.y.toString()+","+ (height).toString() +"\n")
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
