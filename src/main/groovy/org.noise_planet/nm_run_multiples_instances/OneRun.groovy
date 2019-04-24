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

import org.h2gis.functions.io.csv.CSVDriverFunction
import org.h2gis.utilities.wrapper.ConnectionWrapper
import org.cts.op.CoordinateOperationFactory

import java.sql.Connection
import java.sql.DriverManager
import java.sql.PreparedStatement
import java.sql.SQLException

import org.noise_planet.noisemodelling.propagation.ComputeRays
import org.noise_planet.noisemodelling.propagation.ComputeRaysOut
import org.noise_planet.noisemodelling.propagation.FastObstructionTest
import org.noise_planet.noisemodelling.propagation.IComputeRaysOut
import org.noise_planet.noisemodelling.propagation.PropagationPath
import org.noise_planet.noisemodelling.propagation.PropagationProcessData
import org.noise_planet.noisemodelling.propagation.PropagationProcessPathData
import org.noise_planet.noisemodelling.propagation.jdbc.PointNoiseMap
import org.noise_planet.noisemodelling.emission.RSParametersCnossos
import org.noise_planet.noisemodelling.emission.EvaluateRoadSourceCnossos

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

        String workspace_output = "D:\\aumond\\Documents\\PROJETS\\NOISEMODELLING\\debug\\output"
        // le workspace output doit respecter une arborescence specifique
        //output
        //|- data
        //|- config
        System.out.println("Total memory (bytes): " +
                Runtime.getRuntime().totalMemory())

        boolean clean = true //true // reset the postgis database

        boolean loadRays = false

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
            def database_path = "D:/db3/database"
            def database_file = new File(database_path+".mv.db")
            if (database_file.exists()) {
                database_file.delete()
            }
            h2_url = "jdbc:h2:/"+database_path+";DB_CLOSE_DELAY=30;DEFRAG_ALWAYS=TRUE"
            user_name = "sa"
            user_password = ""
        }


        // Noms des tables en entrée et des attributs

        String sources_table_name = "roads_src_zone" // ne pas modifier
        String sources_table_name2 = "ROADS_TRAFFIC_ZONE_CAPTEUR_format2"

        @SourceURI
        URI sourceUri
        Path scriptLocation = Paths.get(sourceUri)
        String rootPath = scriptLocation.getParent().getParent().getParent().getParent().toString()

        // Paramètres de propagation
        int reflexion_order = 0
        int diffraction_order = 0
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

        if(H2) {
            sql.execute("CREATE ALIAS IF NOT EXISTS H2GIS_SPATIAL FOR \"org.h2gis.functions.factory.H2GISFunctions.load\";")
            sql.execute("CALL H2GIS_SPATIAL();")
        }

        // Dossier des résultats
        File newFile = new File(workspace_output + "/config/")
        if(newFile.exists()) {
            FileUtils.cleanDirectory(newFile)
        } else {
            newFile.mkdir()
        }
        newFile = new File(workspace_output + "/data/")
        if(newFile.exists()) {
            FileUtils.cleanDirectory(newFile)
        } else {
            newFile.mkdir()
        }


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
                sql.execute('drop table land_use_zone_capteur2 if exists ')
                sql.execute('drop table dem_lite2 if exists ')
                sql.execute('drop table roads_src_zone if exists ')
                sql.execute('drop table ROADS_TRAFFIC_ZONE_CAPTEUR_format2 if exists ')

                sql.execute("CALL File_table('D:\\aumond\\Documents\\Recherche_hors_projet\\2019_03_GCorbeau_Oiseaux\\LastRun\\couches_clean1\\zone.shp','zone_cense_2km')")
                sql.execute("CALL File_table('D:\\aumond\\Documents\\Recherche_hors_projet\\2019_03_GCorbeau_Oiseaux\\LastRun\\couches_clean1\\batiments2.shp','buildings_zone')")
                sql.execute("CALL File_table('D:\\aumond\\Documents\\Recherche_hors_projet\\2019_03_GCorbeau_Oiseaux\\LastRun\\couches_clean1\\recv.shp','receivers2')")
                sql.execute("CALL File_table('D:\\aumond\\Documents\\Recherche_hors_projet\\2019_03_GCorbeau_Oiseaux\\LastRun\\couches_clean1\\occsol.shp','land_use_zone_capteur2')")
                sql.execute("CALL File_table('D:\\aumond\\Documents\\Recherche_hors_projet\\2019_03_GCorbeau_Oiseaux\\LastRun\\couches_clean1\\DEM.shp','DEM_LITE2')")
                sql.execute("CALL File_table('D:\\aumond\\Documents\\Recherche_hors_projet\\2019_03_GCorbeau_Oiseaux\\LastRun\\couches_clean1\\vide.shp','roads_src_zone')")
                sql.execute("CALL File_table('D:\\aumond\\Documents\\Recherche_hors_projet\\2019_03_GCorbeau_Oiseaux\\LastRun\\couches_clean1\\route_adapt2.shp','ROADS_TRAFFIC_ZONE_CAPTEUR_format2')")

                sql.execute("create spatial index on zone_cense_2km(the_geom)")
                sql.execute("create spatial index on buildings_zone(the_geom)")
                sql.execute("create spatial index on receivers2(the_geom)")
                sql.execute("create spatial index on land_use_zone_capteur2(the_geom)")
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

            System.out.println("Import Tables all...")
            // Import GeometrySoilType

            //System.out.println("LandUse...")
            //fetchCellSoilAreas(connection, geoWithSoilTypeList, soil_table_name, g_field_name)

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


            if (!loadRays) {
                System.out.println("Compute Rays...")
                // Configure noisemap with specified receivers
                //-----------------------------------------------------------------
                // ----------- ICI On calcul les rayons entre sources et recepteurs (c est pour ça r=0, on le fait qu'une fois)
                //-----------------------------------------------------------------

                PointNoiseMap pointNoiseMap = new PointNoiseMap("BUILDINGS_ZONE", "ROADS_SRC_ZONE", "RECEIVERS2")
                pointNoiseMap.setComputeHorizontalDiffraction(false)
                pointNoiseMap.setComputeVerticalDiffraction(false)
                pointNoiseMap.setSoundReflectionOrder(0)
                pointNoiseMap.setHeightField("HAUTEUR")
                pointNoiseMap.setDemTable("DEM_LITE2")
                pointNoiseMap.setMaximumPropagationDistance(5000)
                pointNoiseMap.setMaximumReflectionDistance(50)
                pointNoiseMap.setWallAbsorption(0.1)
                pointNoiseMap.setSoilTableName("LAND_USE_ZONE_CAPTEUR2")
                pointNoiseMap.setThreadCount(8)


                pointNoiseMap.initialize(connection, new EmptyProgressVisitor())
                pointNoiseMap.setComputeRaysOutFactory(new JDBCComputeRaysOut())

                List<ComputeRaysOut.verticeSL> allLevels = new ArrayList<>()
                Set<Long> receivers_ = new HashSet<>()
                for (int i = 0; i < pointNoiseMap.getGridDim(); i++) {
                    for (int j = 0; j < pointNoiseMap.getGridDim(); j++) {
                        IComputeRaysOut out = pointNoiseMap.evaluateCell(connection, i, j, new EmptyProgressVisitor(), receivers_)
                        if (out instanceof ComputeRaysOut) {
                            allLevels.addAll(((ComputeRaysOut) out).getVerticesSoundLevel())
                        }
                    }
                }

                allLevels



            }

             // Ici on rentre dans la phase calcul de la matrice de transfer

            System.out.println("Compute Attenuation...")

            Kryo kryo = new Kryo()
            MapSerializer serializer = new MapSerializer()
            kryo.register(HashMap.class, serializer)
            kryo.register(LinkedHashMap.class, serializer)
            serializer.setKeyClass(String.class, kryo.getSerializer(String.class))
            serializer.setKeysCanBeNull(false)
            serializer.setKeyClass(String.class, kryo.getSerializer(String.class))

           String filenamebin = workspace_output + "\\Rays.bin"

            if (saveRays){
                Output outputBin = new Output(new FileOutputStream(filenamebin))
                kryo.writeObject(outputBin, out)
                outputBin.close()
            }
/*
            if (loadRays){
                Input input = new Input(new FileInputStream(filenamebin))
                propaMap = kryo.readObject(input, HashMap.class)
                input.close()
            }


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
            CalcTime_record.close() */

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

    private static class JDBCPropagationData implements PointNoiseMap.PropagationProcessDataFactory {
        @Override
        PropagationProcessData create(FastObstructionTest freeFieldFinder) {
            return new DirectPropagationProcessData(freeFieldFinder);
        }
    }

    private static class JDBCComputeRaysOut implements PointNoiseMap.IComputeRaysOutFactory {
        @Override
        IComputeRaysOut create(PropagationProcessData threadData, PropagationProcessPathData pathData) {
            return new RayOut(false, pathData, (DirectPropagationProcessData)threadData)
        }
    }

    private static class RayOut extends ComputeRaysOut {
        private DirectPropagationProcessData processData

        RayOut(boolean keepRays, PropagationProcessPathData pathData, DirectPropagationProcessData processData) {
            super(keepRays, pathData, processData)
            this.processData = processData
        }

        @Override
        double[] computeAttenuation(PropagationProcessPathData pathData, long sourceId, double sourceLi, long receiverId, List<PropagationPath> propagationPath) {
            double[] attenuation = super.computeAttenuation(pathData, sourceId, sourceLi, receiverId, propagationPath);
            double[] soundLevel = ComputeRays.wToDba(ComputeRays.multArray(processData.wjSources.get((int)sourceId), ComputeRays.dbaToW(attenuation)));
            return soundLevel
        }
    }

    private static class DirectPropagationProcessData extends PropagationProcessData {
        private List<double[]> wjSources = new ArrayList<>();
        private final static String[] powerColumns = ["db_m63", "db_m125","db_m250", "db_m500", "db_m1000", "db_m2000", "db_m4000", "db_m8000"]

        DirectPropagationProcessData(FastObstructionTest freeFieldFinder) {
            super(freeFieldFinder);
        }


        @Override
        void addSource(Long pk, Geometry geom, SpatialResultSet rs)  throws SQLException {
            super.addSource(pk, geom, rs)
            def sl = new double[powerColumns.length]
            int i = 0
            for(String columnName : powerColumns) {
                sl[i++] = ComputeRays.dbaToW(rs.getDouble(columnName))
            }
            wjSources.add(sl)
        }

        @Override
        double[] getMaximalSourcePower(int sourceId) {
            return wjSources.get(sourceId)
        }
    }
}