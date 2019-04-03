----------------------
-- INPUT PARAMETERS --
----------------------

--________________--
-- ZONE CENSE 2KM --
drop table if exists zone_capteur;
create table  zone_capteur as select * from zone_cense_2km;
ALTER TABLE zone_capteur ALTER COLUMN the_geom TYPE geometry USING ST_SetSRID(the_geom,2154);
CREATE INDEX zone_capteur_the_geom_gist ON zone_capteur USING GIST (the_geom);
-- ZONE TEST --
--drop table if exists zone_capteur;
--create table zone_capteur as select pk pk, ST_SetSRID(the_geom, 2154) the_geom from FENCE_TEST;
--ALTER TABLE zone_capteur ALTER COLUMN the_geom TYPE geometry USING ST_SetSRID(the_geom,2154);
--CREATE INDEX zone_capteur_the_geom_gist ON zone_capteur USING GIST (the_geom);
-- create linked table zone_cense_2km ('org.orbisgis.postgis_jts.Driver', 
--                                     'jdbc:postgresql_h2://137.121.123.1:5432/cense', 
--                                     'cense_read', '9vugV6to25u9', 
--                                     'geo_results', 'zone_cense_2km');

--___________--
-- BUIDLINGS --

drop table if exists buildings;
create table  buildings as select * from cerema_ppbe_lorient_batiment_pop;
ALTER TABLE buildings ALTER COLUMN the_geom TYPE geometry USING ST_SetSRID(the_geom,2154);
-- create linked table buildings ('org.orbisgis.postgis_jts.Driver', 
--                                'jdbc:postgresql_h2://137.121.123.1:5432/cense', 
--                                'cense_read', '9vugV6to25u9', 
--                                'geo_reference', 'cerema_ppbe_lorient_batiment_pop');

-- Extract buildings for the area zone_capteur
drop table if exists buildings_zone_capteur;
create table buildings_zone_capteur as
    select b.*
    from buildings as b, 
         zone_capteur as z
    where ST_intersects(b.the_geom, z.the_geom);


-- Simplify the geometry
drop table if exists buildings_simp_merge_0;
create table buildings_simp_merge_0 as 
    select ST_Union(ST_Accum(ST_SimplifyPreserveTopology(ST_Buffer(the_geom,-0.1,'endcap=square join=bevel'), 0.1))) as the_geom, 
           hauteur, pop_bati
    from buildings_zone_capteur
    group by hauteur, pop_bati;

drop table if exists buildings_simp_merge;
CREATE TABLE buildings_simp_merge AS SELECT pop_bati, hauteur, (ST_DUMP(the_geom)).geom::geometry(Polygon,0) AS the_geom FROM buildings_simp_merge_0 where ST_IsValid(the_geom); 

-- Explode the geometry
drop table if exists buildings_zone_capteur;
create table buildings_zone_capteur AS select ST_FORCE3DZ(the_geom) the_geom, hauteur height, pop_bati pop_bati from buildings_simp_merge;
ALTER TABLE buildings_zone_capteur ADD COLUMN id_build SERIAL PRIMARY KEY;
CREATE INDEX buildings_zone_capteur_the_geom_gist ON buildings_zone_capteur USING GIST (the_geom);

-- Explode the geometry
drop table if exists buildings_zone_capteur2D;
create table buildings_zone_capteur2D as 
    select ST_FORCE2D(the_geom) the_geom, hauteur height , pop_bati pop_bati
    from buildings_simp_merge;
ALTER TABLE buildings_zone_capteur2D ADD COLUMN id_build SERIAL PRIMARY KEY;


--CREATE INDEX buildings_zone_capteur_the_geom_gist ON buildings_zone_capteur USING GIST (the_geom);

--_______--
-- ROADS --

-- Load roads
drop table if exists roads;
create table  roads as select * from cerema_ppbe_lorient_route;


-- create linked table roads ('org.orbisgis.postgis_jts.Driver', 
--                            'jdbc:postgresql_h2://137.121.123.1:5432/cense', 
--                            'cense_read', '9vugV6to25u9', 
--                            'geo_reference', 'cerema_ppbe_lorient_route');




-- Projection from Lambert CC48
drop table if exists roads_proj;
create table roads_proj as
    select pk, ST_Transform(ST_SetSRID(the_geom, 2154), 2154) as the_geom, 
    CASE WHEN importance ='NC' THEN 0
            ELSE importance::float
       END  importance,
           prec_plani, prec_alti, nature, largeur, nb_voies, z_ini, z_fin, 
           vitesse, tmja, annee, tmja_lorie, src_tmja, pourc_pl, src_pr_pl, vit_vl, vit_pl, pk2,id,
           mt as flow_d, me as flow_e, mn as flow_n, pt as pourc_pl_d, pe as pourc_pl_e, pn as pourc_pl_n
    from roads;
drop table if exists roads;



-- Extract roads for the area zone_capteur
drop table if exists roads_zone_capteur;
create table roads_zone_capteur as
    select r.pk, r.the_geom,
    prec_plani, prec_alti, nature, largeur, nb_voies, z_ini, z_fin, 
           vitesse, tmja, annee, tmja_lorie, src_tmja, pourc_pl, src_pr_pl, vit_vl, vit_pl, pk2,id, 
            flow_d, flow_e, flow_n, pourc_pl_d, pourc_pl_e, pourc_pl_n, importance
    from roads_proj as r, 
         zone_capteur as z
    where ST_intersects(r.the_geom, z.the_geom);


-- Simplify the geometry
drop table if exists roads_simp_merge_0;
create table roads_simp_merge_0 as 
    select ST_FORCE2D(ST_SimplifyPreserveTopology(the_geom, 0.1)) as the_geom, 
           prec_plani, prec_alti, nature, largeur, nb_voies, z_ini, z_fin, 
           vitesse, tmja, annee, tmja_lorie, src_tmja, pourc_pl, src_pr_pl, vit_vl, vit_pl, pk2,id,  
           flow_d, flow_e, flow_n, pourc_pl_d, pourc_pl_e, pourc_pl_n, importance
    from roads_zone_capteur;

drop table if exists roads_simp_merge;
CREATE TABLE roads_simp_merge AS SELECT prec_plani, prec_alti, nature, largeur, nb_voies, z_ini, z_fin, importance,
           vitesse, tmja, annee, tmja_lorie, src_tmja, pourc_pl, src_pr_pl, vit_vl, vit_pl, pk2,id,  
           flow_d, flow_e, flow_n, pourc_pl_d, pourc_pl_e, pourc_pl_n, (ST_DUMP(the_geom)).geom::geometry(LineString,0) AS the_geom FROM roads_simp_merge_0; 


-- Merge the geometries
drop table if exists roads_zone_capteur;
create table roads_zone_capteur as 
    select ST_SubDivide(ST_Segmentize(the_geom::geometry,5),8) the_geom, prec_plani, prec_alti, largeur, nb_voies, z_ini, z_fin, 
           vitesse, tmja, annee, tmja_lorie, src_tmja, pourc_pl, src_pr_pl, flow_d, flow_e, flow_n, 
           pourc_pl_d, pourc_pl_e, pourc_pl_n, vit_vl, vit_pl, pk2,id as idst, importance
    from roads_simp_merge;
drop table if exists roads_simp_merge;
ALTER TABLE roads_zone_capteur ADD COLUMN id SERIAL PRIMARY KEY;

-- Compute road traffic parameters
drop table if exists roads_traffic_zone_capteur;
create table roads_traffic_zone_capteur as 
    select id, the_geom, vit_vl as speed_lv, vit_pl as speed_hv, 
           z_ini as beginZ, z_fin as endZ, 
           -- Day
           cast(flow_d*(100-pourc_pl_d)/100. as float) as flow_lv_hour_d, 
           cast(flow_d*pourc_pl_d/100. as float) as flow_hv_hour_d,
           -- Evening
           cast(flow_e*(100-pourc_pl_e)/100. as float) as flow_lv_hour_e, 
           cast(flow_e*pourc_pl_e/100. as float) as flow_hv_hour_e, 
           -- Night
           cast(flow_n*(100-pourc_pl_n)/100. as float) as flow_lv_hour_n, 
           cast(flow_n*pourc_pl_n/100. as float) as flow_hv_hour_n, importance
    from roads_zone_capteur;
-- Compute road traffic parameters

ALTER TABLE roads_traffic_zone_capteur ALTER COLUMN the_geom TYPE geometry USING ST_SetSRID(the_geom,2154);

--________--
-- Junctions --
-- Load junctions

drop table if exists junction;
--create table junction as
--SELECT a.pk2 ida, b.pk2 idb, ST_Intersection(a.the_geom, b.the_geom)
--FROM roads_zone_capteur AS a, roads_zone_capteur AS b 
--WHERE ST_Intersects(a.the_geom,b.the_geom) 
--  AND ST_Touches(a.the_geom, b.the_geom) 
--  AND a.pk2 < b.pk2 
--  AND a.idst <> b.idst
--GROUP BY  a.pk2, b.pk2, ST_Intersection(a.the_geom, b.the_geom);
--
--
--CREATE INDEX roads_proj_the_geom_gist ON roads_proj USING GIST (the_geom);
--drop table if exists junction2;
--CREATE TABLE junction2 as SELECT ST_Intersection(a.the_geom, b.the_geom), Count(Distinct a.id) FROM roads_proj as a, roads_proj as b WHERE ST_Touches(a.the_geom, b.the_geom) AND a.id != b.id GROUP BY ST_Intersection(a.the_geom, b.the_geom);
--
--
--drop table if exists Junc_Line;
--create table Junc_Line as select * from Junction_line;
--ALTER TABLE Junc_Line ALTER COLUMN the_geom TYPE geometry USING ST_SetSRID(the_geom,2154);
--
--drop table if exists Junc_Point;
--create table Junc_Point as select * from Junction_point;
--ALTER TABLE Junc_Point ALTER COLUMN the_geom TYPE geometry USING ST_SetSRID(the_geom,2154);
--
--drop table if exists  t1,t2,t3;
--create table t1 as SELECT a.id id,  0 t,ST_DISTANCE(a.the_geom,ST_ClosestPoint(b.the_geom,a.the_geom)) dist from roads_traffic_zone_capteur a, JUNC_POINT b where ST_DWithin(a.the_geom,b.the_geom,100) and b.highway='mini_roundabout';
--create table t2 as SELECT a.id id ,  1 t,ST_DISTANCE(a.the_geom,ST_ClosestPoint(b.the_geom,a.the_geom)) dist from roads_traffic_zone_capteur a, JUNC_POINT b where ST_DWithin(a.the_geom,b.the_geom,100) and (b.highway='stop' or b.highway='traffic_signals');
--create table t3 as SELECT a.id id,  0 t,ST_DISTANCE(a.the_geom,ST_ClosestPoint(b.the_geom,a.the_geom)) dist from roads_traffic_zone_capteur a, JUNC_LINE b where ST_DWithin(a.the_geom,b.the_geom,100);
--
--drop table  if exists dist_junc;
--create table dist_junc as select t1.id, t1.dist, t1.t from t1 UNION SELECT t2.id, t2.dist,t2.t from t2 UNION SELECT t3.id, t3.dist,t3.t from t3 order by dist;
--
--drop table  if exists dist_junc2;
--create table dist_junc2 as select id, min(dist) from dist_junc d1 group by id;
--
--drop table  if exists dist_junc3;
--create table dist_junc3 as select d1.id, d1.dist, d1.t from dist_junc2 d2, dist_junc d1  where d1.id=d2.id and d1.dist = d2.dist;
--
--drop table if exists roads_traffic_zone_capteur_junc;
--create table roads_traffic_zone_capteur_junc as 
--    select a.*,
--            COALESCE(b.dist,'200') dist,
--            COALESCE(b.t,'0') t     
--    from roads_traffic_zone_capteur a LEFT JOIN dist_junc2 b ON a.id = b.id;
--
--drop table if exists ROADS_TRAFFIC_ZONE_CAPTEUR_format;
--create table ROADS_TRAFFIC_ZONE_CAPTEUR_format as select id id, ST_Translate(ST_Force3DZ(the_geom),0,0,0.05) the_geom,
--speed_lv lv_d_speed,speed_lv mv_d_speed,speed_hv hv_d_speed,speed_lv wav_d_speed,speed_lv wbv_d_speed,
--speed_lv lv_e_speed,speed_lv mv_e_speed,speed_hv hv_e_speed,speed_lv wav_e_speed,speed_lv wbv_e_speed,
--speed_lv lv_n_speed,speed_lv mv_n_speed,speed_hv hv_n_speed,speed_lv wav_n_speed,speed_lv wbv_n_speed,
--flow_lv_hour_d vl_d_per_hour,0 ml_d_per_hour,flow_hv_hour_d pl_d_per_hour,0 wa_d_per_hour,0 wb_d_per_hour,
--flow_lv_hour_e vl_e_per_hour,0 ml_e_per_hour,flow_hv_hour_e pl_e_per_hour,0 wa_e_per_hour,0 wb_e_per_hour,
--flow_lv_hour_n vl_n_per_hour,0 ml_n_per_hour,flow_hv_hour_n pl_n_per_hour,0 wa_n_per_hour,0 wb_n_per_hour,
--beginZ Zstart,endZ Zend,dist Juncdist, t Junc_type, importance road_pav from roads_traffic_zone_capteur_junc;
--
--ALTER TABLE ROADS_TRAFFIC_ZONE_CAPTEUR_format ADD PRIMARY KEY(ID);

--________--
-- DEM --
-- Load dem
drop table if exists dem2;
create table  dem2 as select ST_SetSRID(a.the_geom,2154) the_geom, contour contour from cerema_ppbe_lorient_crb_niv_1m a, zone_capteur b where st_intersects(ST_SetSRID(a.the_geom,2154), b.the_geom) ;
CREATE INDEX dem2_the_geom_gist ON dem2 USING GIST (the_geom);

-- create linked table ground ('org.orbisgis.postgis_jts.Driver', 
--                             'jdbc:postgresql_h2://137.121.123.1:5432/cense', 
--                             'cense_read', '9vugV6to25u9', 
--                             'geo_reference', 'cerema_ppbe_lorient_absorption');
drop table if exists dem2_points;
create table dem2_points as SELECT (dp).geom As the_geom, z as Z FROM (select ST_DumpPoints(the_geom) AS dp, contour as z from dem2) as foo;
CREATE INDEX dem2_points_the_geom_gist ON dem2 USING GIST (the_geom);

--________--
-- GROUND --

-- Load ground
drop table if exists ground;
create table  ground as select * from cerema_ppbe_lorient_absorption;
ALTER TABLE ground ALTER COLUMN the_geom TYPE geometry USING ST_SetSRID(the_geom,2154);
CREATE INDEX ground_the_geom_gist ON ground USING GIST (the_geom);

-- create linked table ground ('org.orbisgis.postgis_jts.Driver', 
--                             'jdbc:postgresql_h2://137.121.123.1:5432/cense', 
--                             'cense_read', '9vugV6to25u9', 
--                             'geo_reference', 'cerema_ppbe_lorient_absorption');

-- Extract ground for the area zone_capteur
drop table if exists land_use_zone_capteur;
create table land_use_zone_capteur as
    select g.pk, g.the_geom, g.babs as G
    from ground as g, 
         zone_capteur as z
    where ST_intersects(g.the_geom, z.the_geom);
CREATE INDEX land_use_zone_capteur_the_geom_gist ON land_use_zone_capteur USING GIST (the_geom);

-- Simplify the geometry
drop table if exists land_use_simp_merge_0;
create table land_use_simp_merge_0 as 
    select ST_SimplifyPreserveTopology(a.the_geom,0.1) as the_geom, a.G 
    from land_use_zone_capteur a  ;

drop table if exists land_use_simp_merge_1;
create table land_use_simp_merge_1 as select a.the_geom, a.g from land_use_simp_merge_0 a where ST_IsValid(the_geom)=true;

drop table if exists land_use_simp_merge;
CREATE TABLE land_use_simp_merge AS SELECT G , (ST_DUMP(the_geom)).geom::geometry AS the_geom FROM land_use_simp_merge_1; 

-- Explode the geometry
drop table if exists land_use_zone_capteur_tmp;
create table land_use_zone_capteur_tmp as 
    select  the_geom, G 
    from land_use_simp_merge;
drop table if exists land_use_simp_merge;
ALTER TABLE land_use_zone_capteur_tmp ADD COLUMN id SERIAL PRIMARY KEY;

--CREATE EXTENSION postgis_sfcgal;
drop table if exists land_use_zone_capteur_0;
create table land_use_zone_capteur_0 as select ST_DelaunayTriangles(the_geom,1,0) as the_geom, G 
    from land_use_zone_capteur_tmp where ST_IsValid(the_geom)=true;


-- Optimise land use table for faster computation
drop table if exists land_use_zone_capteur;
create table land_use_zone_capteur as 
    select (ST_DUMP(the_geom)).geom::geometry AS the_geom, G 
    from land_use_zone_capteur_0;


-- Export data to shapefiles
--CALL SHPWrite('D:\aumond\Documents\CENSE\WP2\Analyses\Incertitudes\Incertitudes\data\output\land_use_zone_capteur.shp', 'LAND_USE_ZONE_CAPTEUR');

---------------------------
-- CREATE RECEIVERS GRID --
---------------------------

drop table if exists receivers_build_0;
create table receivers_build_0 as SELECT ID_build build_id, ST_ExteriorRing(ST_Buffer(ST_SimplifyPreserveTopology(b.the_geom,2), 2, 'quad_segs=0 endcap=flat')) the_geom  from BUILDINGS_ZONE_CAPTEUR b ;
ALTER TABLE receivers_build_0 ADD COLUMN id SERIAL PRIMARY KEY;

drop table if exists receivers_build;
CREATE TABLE receivers_build AS SELECT build_ID, (ST_DUMP(the_geom)).geom::geometry(LineString,0) AS the_geom FROM receivers_build_0; 

--drop table if exists indexed_points;
--CREATE TABLE indexed_points AS SELECT build_ID,(ST_DUMP(the_geom)).geom::geometry(Polygon,0) AS the_geom FROM receivers_build_0; 
--ALTER TABLE indexed_points ADD COLUMN id SERIAL PRIMARY KEY;

--drop table if exists receivers_delete;
--create table receivers_delete as SELECT r.ID, r.the_geom,r.build_ID from indexed_points r,  BUILDINGS_ZONE_CAPTEUR b where st_intersects(b.the_geom, r.the_geom) and b.the_geom && r.the_geom;
--delete from indexed_points r where exists (select 1 from receivers_delete rd where r.ID=rd.ID);
--drop table if exists receivers_delete;


drop table if exists bb;
create table bb as select ST_EXPAND(ST_Collect(ST_ACCUM(b.the_geom)),-2000,-2000)  the_geom from BUILDINGS_ZONE_CAPTEUR2D b;

drop table if exists receivers_build_ratio; 
create table receivers_build_ratio as select a.* from receivers_build_0 a, bb b where st_intersects(b.the_geom, a.the_geom) ORDER BY random() LIMIT 10;


drop table if exists indexed_points;
--create table indexed_points(the_geom geometry, gid integer);
create table indexed_points(old_edge_id integer, the_geom geometry, number_on_line integer, gid integer);

CREATE OR REPLACE FUNCTION VBEB_CLASSE2(wanted_len double precision)
RETURNS integer AS
$$
DECLARE
    current_fractional double precision := 0.0;
    current_number_of_point integer := 1;
    i record;
BEGIN
    FOR i IN 
        SELECT id as id_column, st_transform(the_geom, 2154) as the_geom, st_length(st_transform(the_geom, 2154)) as line_length, build_id as build_id FROM receivers_build_ratio 
        LOOP
        current_fractional := 0.0;
        WHILE current_fractional <= (1.0)::double precision 
            LOOP
                INSERT INTO indexed_points(old_edge_id, the_geom, number_on_line, gid) VALUES (i.id_column, ST_LineInterpolatePoint(i.the_geom, current_fractional), current_number_of_point, i.build_id);
                current_fractional := current_fractional + (wanted_len/ i.line_length); 
                current_number_of_point := current_number_of_point + 1;
            END LOOP;
       END LOOP;
       return i.id_column;
    END;
$$ LANGUAGE PLPGSQL ;

select VBEB_CLASSE2(5);


ALTER TABLE indexed_points ADD COLUMN id SERIAL PRIMARY KEY;
CREATE INDEX indexed_points_the_geom_gist ON indexed_points USING GIST (the_geom);

drop table if exists receivers_delete;
create table receivers_delete as SELECT r.ID, r.the_geom,r.gid build_id from indexed_points r,  BUILDINGS_ZONE_CAPTEUR b where st_intersects(b.the_geom, r.the_geom);
delete from indexed_points r where exists (select 1 from receivers_delete rd where r.ID=rd.ID);
drop table if exists receivers_delete;

alter table indexed_points DROP column id ;
alter table indexed_points add column id serial ;
ALTER TABLE indexed_points ADD PRIMARY KEY (ID);

drop table if exists receivers;
create table receivers as select id id, gid gid, ST_Translate(ST_force3d(the_geom),0,0,4) the_geom from indexed_points;
ALTER TABLE receivers ADD PRIMARY KEY (ID);
CREATE INDEX receivers_the_geom_gist ON receivers USING GIST (the_geom);


-- delete all sources that are not necessary
drop table if exists cc;
create table cc as select ST_Union(ST_EXPAND(b.the_geom,250,250))  the_geom from receivers b;
drop table if exists source_delete;
create table source_delete as SELECT r.* from ROADS_TRAFFIC_ZONE_CAPTEUR_format r,  cc b where not st_intersects(b.the_geom, r.the_geom);
delete from ROADS_TRAFFIC_ZONE_CAPTEUR_format r where exists (select 1 from source_delete rd where r.ID=rd.ID);
drop table if exists source_delete;





--drop table if exists receivers_build;
--CREATE TABLE receivers_build AS SELECT build_ID, (ST_DUMP(the_geom)).geom::geometry(LineString,0) AS the_geom FROM receivers_build_0; 
-- Export data to shapefiles
--CALL SHPWrite('D:\aumond\Documents\CENSE\WP2\Analyses\Incertitudes\Incertitudes\data\output\indexed_points.shp', 'INDEXED_POINTS');




