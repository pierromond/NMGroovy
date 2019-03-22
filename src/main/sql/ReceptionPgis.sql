----------------------------------------------
-- Compute sound level on receivers --
----------------------------------------------

drop table if exists lvl_receiver_lvl_day_zone;
create table lvl_receiver_lvl_day_zone as select 
        r.the_geom the_geom,
        r.id id,
        r.gid build_id,
           10.*log(
       power(10,10.*log(sum(power(10,(s.db_m_d63+l.att63)/10)))/10)
        +power(10,10.*log(sum(power(10,(s.db_m_d125+l.att125)/10)))/10)
        +power(10,10.*log(sum(power(10,(s.db_m_d250+l.att250)/10)))/10)
        +power(10,10.*log(sum(power(10,(s.db_m_d500+l.att500)/10)))/10)
        +power(10,10.*log(sum(power(10,(s.db_m_d1000+l.att1000)/10)))/10)
        +power(10,10.*log(sum(power(10,(s.db_m_d2000+l.att2000)/10)))/10)
        +power(10,10.*log(sum(power(10,(s.db_m_d4000+l.att4000)/10)))/10)
        +power(10,10.*log(sum(power(10,(s.db_m_d8000+l.att8000)/10)))/10)
       ) db      
        from receiver_lvl_day_zone l, receivers r,  roads_src_zone s
        where l.idsource=s.id and r.id = l.idrecepteur and s.db_m_d63+l.att63>0
        group by r.id, r.the_geom, r.gid;

drop table if exists lvl_receiver_lvl_evening_zone;
create table lvl_receiver_lvl_evening_zone as select 
        r.the_geom the_geom,
        r.id id,
        r.gid build_id,
           10.*log(
       power(10,10.*log(sum(power(10,(s.db_m_e63+l.att63)/10)))/10)
        +power(10,10.*log(sum(power(10,(s.db_m_e125+l.att125)/10)))/10)
        +power(10,10.*log(sum(power(10,(s.db_m_e250+l.att250)/10)))/10)
        +power(10,10.*log(sum(power(10,(s.db_m_e500+l.att500)/10)))/10)
        +power(10,10.*log(sum(power(10,(s.db_m_e1000+l.att1000)/10)))/10)
        +power(10,10.*log(sum(power(10,(s.db_m_e2000+l.att2000)/10)))/10)
        +power(10,10.*log(sum(power(10,(s.db_m_e4000+l.att4000)/10)))/10)
        +power(10,10.*log(sum(power(10,(s.db_m_e8000+l.att8000)/10)))/10)
       ) db      
        from receiver_lvl_evening_zone l, receivers r,  roads_src_zone s
        where l.idsource=s.id and r.id = l.idrecepteur and s.db_m_e63+l.att63>0
        group by r.id, r.the_geom, r.gid;
            
        
drop table if exists lvl_receiver_lvl_night_zone;
create table lvl_receiver_lvl_night_zone as select 
        r.the_geom the_geom,
        r.id id,
        r.gid build_id,
           10.*log(
       power(10,10.*log(sum(power(10,(s.db_m_n63+l.att63)/10)))/10)
        +power(10,10.*log(sum(power(10,(s.db_m_n125+l.att125)/10)))/10)
        +power(10,10.*log(sum(power(10,(s.db_m_n250+l.att250)/10)))/10)
        +power(10,10.*log(sum(power(10,(s.db_m_n500+l.att500)/10)))/10)
        +power(10,10.*log(sum(power(10,(s.db_m_n1000+l.att1000)/10)))/10)
        +power(10,10.*log(sum(power(10,(s.db_m_n2000+l.att2000)/10)))/10)
        +power(10,10.*log(sum(power(10,(s.db_m_n4000+l.att4000)/10)))/10)
        +power(10,10.*log(sum(power(10,(s.db_m_n8000+l.att8000)/10)))/10)
       ) db      
        from receiver_lvl_night_zone l, receivers r,  roads_src_zone s
        where l.idsource=s.id and r.id = l.idrecepteur and s.db_m_n63+l.att63>0
        group by r.id, r.the_geom, r.gid;



--_______________________________
-- Compute the Lden

drop table if exists lday, levening, lnight;
create table lday as select id, db, build_id from lvl_receiver_lvl_day_zone;
create table levening as select id, db, build_id from lvl_receiver_lvl_evening_zone;
create table lnight as select id, db, build_id from lvl_receiver_lvl_night_zone;


drop table if exists lden;
create table lden as
    select ld.id id, 
           10.*log((12./24.)*power(10, ld.db/10.) + (4./24.)*power(10, (le.db+5.)/10.) + (8./24.)*power(10, (ln.db+10.)/10.)) lden, 
           ld.build_id build_id
    from lday ld,
         levening le,
         lnight ln 
    where ld.id=le.id 
    and le.id=ln.id;

alter table lden add column pk serial ;
alter table lden alter column id set not null;
alter table lden add primary key (id);

--------------------------
-- Exposed population --
--------------------------

--drop table  if exists pop_lvl_den;
--create table pop_lvl_den as select b.id_build ,  b.pop_bati , max(a.lden) db from lden a, buildings2 b where a.build_id = b.id_build group by b.id_build;

--drop table  if exists pop_lvl_n;
--create table pop_lvl_n as select b.id_build ,  b.pop_bati , max(a.db)  db from lvl_receiver_lvl_night_zone a, buildings_zone2d b where a.build_id = b.id_build group by b.id_build;
