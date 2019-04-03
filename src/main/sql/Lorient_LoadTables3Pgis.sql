
-- Fix the height of sources at 0.05 m
drop table if exists roads_src_zone_capteur;
create table roads_src_zone_capteur as 
    select id, the_geom, 
           db_m_d63, db_m_d125, db_m_d250, db_m_d500, db_m_d1000,db_m_d2000, db_m_d4000, db_m_d8000, 
           db_m_e63, db_m_e125, db_m_e250, db_m_e500, db_m_e1000,db_m_e2000, db_m_e4000, db_m_e8000, 
           db_m_n63, db_m_n125, db_m_n250, db_m_n500, db_m_n1000,db_m_n2000, db_m_n4000, db_m_n8000
    from roads_src_fqs_zone_capteur;
    
--alter table roads_src_zone_capteur DROP column id ;
--alter table roads_src_zone_capteur add column id serial ;
ALTER TABLE roads_src_zone_capteur ADD PRIMARY KEY (ID);

drop table if exists roads_src_zone_capteur100;
create table roads_src_zone_capteur100 as select id, the_geom, 
           100 db_m63, 100 db_m125, 100 db_m250, 100 db_m500, 100 db_m1000,100 db_m2000, 100 db_m4000, 100 db_m8000
            from roads_src_zone_capteur;
ALTER TABLE roads_src_zone_capteur100 ADD PRIMARY KEY (ID);
