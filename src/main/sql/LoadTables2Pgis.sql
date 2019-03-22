-----------------------------
-- LOAD AND PROCESS TABLES --
-----------------------------

-- ____________ --
-- Study region --
/*
drop table if exists zone;
create table zone as select * from zone_cense_2km;
alter table zone rename column geom to the_geom;
alter table zone alter column the_geom type geometry using ST_SetSRID(the_geom, 2154);
create index zone_the_geom_gist on zone using GIST (the_geom);
*/
/*
---- Connection to the Cense PostGreSQL database --
-- drop table if exists zone;
-- create table zone as select pk as pk, ST_SetSRID(the_geom, 2154) as the_geom from FENCE_TEST;
-- alter table zone alter column the_geom type geometry using ST_SetSRID(the_geom, 2154);
-- create index zone_the_geom_gist on zone using GIST (the_geom);
-- create linked table zone_cense_2km ('org.orbisgis.postgis_jts.Driver',
--                                     'jdbc:postgresql_h2://137.121.123.1:5432/cense', 
--                                     'cense_read', '9vugV6to25u9', 
--                                     'geo_results', 'zone_cense_2km');
*/

-- _________ --
-- Buildings --

drop table if exists buildings;
create table buildings as select * from cerema_ppbe_lorient_batiment_pop;
-- alter table buildings rename column geom to the_geom;
alter table buildings alter column the_geom type geometry using ST_SetSRID(the_geom, 2154);
/*
---- Connection to the Cense PostGreSQL database --
-- create linked table buildings ('org.orbisgis.postgis_jts.Driver',
--                                'jdbc:postgresql_h2://137.121.123.1:5432/cense', 
--                                'cense_read', '9vugV6to25u9', 
--                                'geo_reference', 'cerema_ppbe_lorient_batiment_pop');
*/

-- Extract buildings for the study region
drop table if exists buildings_zone;
create table buildings_zone as
    select b.*
    from buildings as b, 
         zone as z
    where ST_intersects(b.the_geom, z.the_geom);

-- Simplify the geometry
drop table if exists buildings_simp_merge_0;
create table buildings_simp_merge_0 as 
    select ST_Union(ST_Accum(ST_SimplifyPreserveTopology(ST_Buffer(the_geom,-0.1,'endcap=square join=bevel'), 0.1)))
        as the_geom, hauteur, pop_bati
    from buildings_zone
    group by hauteur, pop_bati;

drop table if exists buildings_simp_merge;
create table buildings_simp_merge as 
    select pop_bati, hauteur, (ST_Dump(the_geom)).geom::geometry(Polygon,0) as the_geom 
    from buildings_simp_merge_0 
    where ST_IsValid(the_geom);
drop table if exists buildings_simp_merge_0;

-- Explode the geometry
drop table if exists buildings_zone;
create table buildings_zone as
    select ST_Force3DZ(the_geom) the_geom, hauteur as height, pop_bati as pop_bati
    from buildings_simp_merge;
alter table buildings_zone add column id_build serial primary key;
create index buildings_zone_the_geom_gist on buildings_zone using gist (the_geom);

-- Explode the geometry
drop table if exists buildings_2D_zone;
create table buildings_2D_zone as 
    select ST_Force2D(the_geom) as the_geom, hauteur as height, pop_bati as pop_bati
    from buildings_simp_merge;
alter table buildings_2D_zone add column id_build serial primary key;
--create index buildings_zone_the_geom_gist on buildings_zone using gist (the_geom);

--_______--
-- Roads --

-- Load roads
drop table if exists roads;
create table roads as select * from cerema_ppbe_lorient_route;
--alter table roads rename column geom to the_geom;
/*
---- Connection to the Cense PostGreSQL database --
-- create linked table roads ('org.orbisgis.postgis_jts.Driver', 
--                            'jdbc:postgresql_h2://137.121.123.1:5432/cense', 
--                            'cense_read', '9vugV6to25u9', 
--                            'geo_reference', 'cerema_ppbe_lorient_route');
*/

-- Projection from Lambert CC48
drop table if exists roads_proj;
create table roads_proj as
    select pk, ST_Transform(ST_SetSRID(the_geom, 2154), 2154) as the_geom, 
    case when importance ='NC' then 0
            else importance::float
       end importance,
           prec_plani, prec_alti, nature, largeur, nb_voies, z_ini, z_fin, 
           vitesse, tmja, annee, tmja_lorie, src_tmja, pourc_pl, src_pr_pl, vit_vl, vit_pl, id,
           mt as flow_d, me as flow_e, mn as flow_n, pt as pourc_pl_d, pe as pourc_pl_e, pn as pourc_pl_n
    from roads;
drop table if exists roads;

-- Extract roads for the study region
drop table if exists roads_zone;
create table roads_zone as
    select r.pk, r.the_geom, prec_plani, prec_alti, nature, largeur, nb_voies, z_ini, z_fin, 
           vitesse, tmja, annee, tmja_lorie, src_tmja, pourc_pl, src_pr_pl, vit_vl, vit_pl, id, 
           flow_d, flow_e, flow_n, pourc_pl_d, pourc_pl_e, pourc_pl_n, importance
    from roads_proj as r, zone as z
    where ST_intersects(r.the_geom, z.the_geom);

-- Simplify the geometry
drop table if exists roads_simp_merge_0;
create table roads_simp_merge_0 as 
    select ST_Force2D(ST_SimplifyPreserveTopology(the_geom, 0.1)) as the_geom, 
           prec_plani, prec_alti, nature, largeur, nb_voies, z_ini, z_fin, 
           vitesse, tmja, annee, tmja_lorie, src_tmja, pourc_pl, src_pr_pl, vit_vl, vit_pl, pk, id,  
           flow_d, flow_e, flow_n, pourc_pl_d, pourc_pl_e, pourc_pl_n, importance
    from roads_zone;

drop table if exists roads_simp_merge;
create table roads_simp_merge as 
    select prec_plani, prec_alti, nature, largeur, nb_voies, z_ini, z_fin, importance,
           vitesse, tmja, annee, tmja_lorie, src_tmja, pourc_pl, src_pr_pl, vit_vl, vit_pl, pk, id,  
           flow_d, flow_e, flow_n, pourc_pl_d, pourc_pl_e, pourc_pl_n,
           (ST_Dump(the_geom)).geom::geometry(LineString,0) as the_geom
    from roads_simp_merge_0;
drop table if exists roads_simp_merge_0;

-- Merge the geometries
drop table if exists roads_zone;
create table roads_zone as 
    select ST_SubDivide(ST_Segmentize(the_geom::geometry,5),8) as the_geom, prec_plani, prec_alti, largeur, nb_voies,
           z_ini, z_fin, vitesse, tmja, annee, tmja_lorie, src_tmja, pourc_pl, src_pr_pl, flow_d, flow_e, flow_n,
           pourc_pl_d, pourc_pl_e, pourc_pl_n, vit_vl, vit_pl, pk, id as idst, importance
    from roads_simp_merge;
drop table if exists roads_simp_merge;
alter table roads_zone add column id serial primary key;

-- Compute road traffic parameters
drop table if exists roads_traffic_zone;
create table roads_traffic_zone as 
    select id, the_geom, vit_vl as speed_lv, vit_pl as speed_hv, z_ini as beginZ, z_fin as endZ,
           -- Day
           cast(flow_d*(100-pourc_pl_d)/100. as float) as flow_lv_hour_d, 
           cast(flow_d*pourc_pl_d/100. as float) as flow_hv_hour_d,
           -- Evening
           cast(flow_e*(100-pourc_pl_e)/100. as float) as flow_lv_hour_e, 
           cast(flow_e*pourc_pl_e/100. as float) as flow_hv_hour_e, 
           -- Night
           cast(flow_n*(100-pourc_pl_n)/100. as float) as flow_lv_hour_n, 
           cast(flow_n*pourc_pl_n/100. as float) as flow_hv_hour_n, importance
    from roads_zone;
alter table roads_traffic_zone 
    alter column the_geom type geometry using ST_SetSRID(the_geom, 2154);

drop table if exists roads_traffic_zone_format;
create table roads_traffic_zone_format as 
    select id id, ST_Translate(ST_Force3DZ(the_geom),0,0,0.05) as the_geom,
           speed_lv lv_d_speed, speed_lv mv_d_speed, speed_hv hv_d_speed, speed_lv wav_d_speed, speed_lv wbv_d_speed,
           speed_lv lv_e_speed, speed_lv mv_e_speed, speed_hv hv_e_speed, speed_lv wav_e_speed, speed_lv wbv_e_speed,
           speed_lv lv_n_speed, speed_lv mv_n_speed, speed_hv hv_n_speed, speed_lv wav_n_speed, speed_lv wbv_n_speed,
           flow_lv_hour_d vl_d_per_hour, 0 ml_d_per_hour,flow_hv_hour_d pl_d_per_hour, 0 wa_d_per_hour, 0 wb_d_per_hour,
           flow_lv_hour_e vl_e_per_hour, 0 ml_e_per_hour,flow_hv_hour_e pl_e_per_hour, 0 wa_e_per_hour, 0 wb_e_per_hour,
           flow_lv_hour_n vl_n_per_hour, 0 ml_n_per_hour,flow_hv_hour_n pl_n_per_hour, 0 wa_n_per_hour, 0 wb_n_per_hour,
           beginZ as Zstart, endZ as Zend, importance as road_pav
    from roads_traffic_zone;
alter table roads_traffic_zone_format add primary key(id);

-- _______________________ --
-- Digital Elevation Model --

-- Extract dem for the study region
drop table if exists dem1;
create table dem1 as select * from cerema_ppbe_lorient_crb_niv_1m;
--alter table dem1 rename column geom to the_geom;
drop table if exists dem2;
create table dem2 as 
    select ST_SetSRID(a.the_geom, 2154) as the_geom, contour as contour
    from dem1 as a, zone as b 
    where ST_Intersects(ST_SetSRID(a.the_geom, 2154), b.the_geom) ;
create index dem2_the_geom_gist on dem2 using GIST (the_geom);
/*
---- Connection to the Cense PostGreSQL database --
-- create linked table ground ('org.orbisgis.postgis_jts.Driver', 
--                             'jdbc:postgresql_h2://137.121.123.1:5432/cense', 
--                             'cense_read', '9vugV6to25u9', 
--                             'geo_reference', 'cerema_ppbe_lorient_absorption');
*/

drop table if exists dem2_points;
create table dem2_points as
    select (dp).geom as the_geom, z as Z
    from (select ST_DumpPoints(the_geom) as dp, contour as z from dem2) as foo;
create index dem2_points_the_geom_gist on dem2 using gist (the_geom);

-- ______ --
-- Ground --

-- Load ground
drop table if exists ground;
create table ground as select * from cerema_ppbe_lorient_absorption;
--alter table ground rename column geom to the_geom;
alter table ground alter column the_geom type geometry using ST_SetSRID(the_geom, 2154);
create index ground_the_geom_gist on ground using GIST (the_geom);
/*
---- Connection to the Cense PostGreSQL database --
-- create linked table ground ('org.orbisgis.postgis_jts.Driver', 
--                             'jdbc:postgresql_h2://137.121.123.1:5432/cense', 
--                             'cense_read', '9vugV6to25u9', 
--                             'geo_reference', 'cerema_ppbe_lorient_absorption');
*/

-- Extract ground for the study region
drop table if exists land_use_zone;
create table land_use_zone as
    select g.the_geom, g.babs as G	-- g.pk, 
    from ground as g, 
         zone as z
    where ST_Intersects(g.the_geom, z.the_geom);
create index land_use_zone_the_geom_gist on land_use_zone using gist (the_geom);

-- Simplify the geometry
drop table if exists land_use_simp_merge_0;
create table land_use_simp_merge_0 as 
    select ST_SimplifyPreserveTopology(a.the_geom,0.1) as the_geom, a.G 
    from land_use_zone a  ;

drop table if exists land_use_simp_merge_1;
create table land_use_simp_merge_1 as 
    select a.the_geom, a.g 
    from land_use_simp_merge_0 a 
    where ST_IsValid(the_geom)=true;
drop table if exists land_use_simp_merge_0;

drop table if exists land_use_simp_merge;
create table land_use_simp_merge as 
    select G, (ST_Dump(the_geom)).geom::geometry as the_geom
    from land_use_simp_merge_1; 

-- Optimise land use table for faster computation
drop table if exists land_use_zone;
create table land_use_zone as 
    select ST_FORCE2D((ST_Dump(the_geom)).geom::geometry) as the_geom, G 
    from land_use_simp_merge_1;	-- land_use_zone_0
drop table if exists land_use_simp_merge_1;

-- Export data to shapefiles
--CALL SHPWrite('D:\aumond\Documents\CENSE\WP2\Analyses\Incertitudes\Incertitudes\data\output\land_use_zone.shp', 'Land_USE_zone');

---------------------------
-- CREATE RECEIVERS GRID --
---------------------------

drop table if exists receivers_build_0;
create table receivers_build_0 as
    select id_build as build_id,
           ST_ExteriorRing(ST_Buffer(ST_SimplifyPreserveTopology(b.the_geom,2),2,'quad_segs=0 endcap=flat')) as the_geom
    from buildings_zone as b;
alter table receivers_build_0 add column id serial primary key;

drop table if exists receivers_build;
create table receivers_build as
    select build_id, (ST_Dump(the_geom)).geom::geometry(LineString,0) as the_geom
    from receivers_build_0;
--drop table if exists indexed_points;
--create table indexed_points as select build_ID,(ST_Dump(the_geom)).geom::geometry(Polygon,0) as the_geom from receivers_build_0; 
--alter table indexed_points add column id SERIAL primary key;

--drop table if exists receivers_delete;
--create table receivers_delete as select r.ID, r.the_geom,r.build_ID from indexed_points r,  BUILDINGS_zone b where ST_Intersects(b.the_geom, r.the_geom) and b.the_geom && r.the_geom;
--delete from indexed_points r where exists (select 1 from receivers_delete rd where r.ID=rd.ID);
--drop table if exists receivers_delete;

drop table if exists bb;
create table bb as 
    select ST_Expand(ST_Collect(ST_Accum(b.the_geom)),-2000,-2000) as the_geom
    from buildings_2D_zone as b;

drop table if exists receivers_build_ratio; 
create table receivers_build_ratio as
    select a.*
    from receivers_build_0 as a, bb as b
    where ST_Intersects(b.the_geom, a.the_geom)
    order by random() limit 10;
drop table if exists receivers_build_0;

drop table if exists indexed_points;
create table indexed_points(old_edge_id integer, the_geom geometry, number_on_line integer, gid integer);

create or replace function VBEB_CLasSE2(wanted_len double precision)
returns integer as
$$
declare
    current_fractional double precision := 0.0;
    current_number_of_point integer := 1;
    i record;
begin
    for i in
        select id as id_column,
               ST_Transform(the_geom, 2154) as the_geom,
               ST_Length(ST_Transform(the_geom, 2154)) as line_length,
               build_id as build_id
            from receivers_build_ratio
        loop
        current_fractional := 0.0;
        while current_fractional <= (1.0)::double precision 
            loop
                insert into indexed_points(old_edge_id, the_geom, number_on_line, gid)
                    values (i.id_column, ST_LineInterpolatePoint(i.the_geom, current_fractional), current_number_of_point, i.build_id);
                current_fractional := current_fractional + (wanted_len / i.line_length);
                current_number_of_point := current_number_of_point + 1;
            end loop;
       end loop;
       return i.id_column;
    end;
$$ language PLPGSQL ;

select VBEB_CLasSE2(5);

alter table indexed_points add column id serial primary key;
create index indexed_points_the_geom_gist on indexed_points using gist (the_geom);

drop table if exists receivers_delete;
create table receivers_delete as
    select r.id, r.the_geom, r.gid build_id
    from indexed_points as r, buildings_zone as b
    where ST_Intersects(b.the_geom, r.the_geom);
delete from indexed_points as r where exists (select 1 from receivers_delete as rd where r.id=rd.id);
drop table if exists receivers_delete;

alter table indexed_points drop column id ;
alter table indexed_points add column id serial ;
alter table indexed_points add primary key (id);

drop table if exists receivers;
create table receivers as
    select id as id, gid as gid, ST_Translate(ST_Force3DZ(the_geom),0,0,4) as the_geom
    from indexed_points;
alter table receivers add primary key (id);
create index receivers_the_geom_gist on receivers using gist (the_geom);

-- delete all sources that are not necessary
drop table if exists cc;
create table cc as
    select ST_Union(ST_Expand(b.the_geom,250,250)) as the_geom
    from receivers b;
drop table if exists source_delete;
create table source_delete as
    select r.*
    from roads_traffic_zone_format as r, cc as b
    where not ST_Intersects(b.the_geom, r.the_geom);
delete from roads_traffic_zone_format as r 
    where exists (select 1 from source_delete as rd where r.id=rd.id);
drop table if exists source_delete;
drop table if exists cc;

-- Export data to shapefiles
--CALL SHPWrite('D:\aumond\Documents\CENSE\WP2\Analyses\Incertitudes\Incertitudes\data\output\indexed_points.shp', 'INDEXED_POINTS');

