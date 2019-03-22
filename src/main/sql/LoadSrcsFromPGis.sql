drop table if exists roads_src_zone;

create table roads_src_zone (id integer, the_geom  geometry, 
				   db_m_d63 double precision, db_m_d125 double precision, db_m_d250 double precision,  db_m_d500 double precision, 
				   db_m_d1000 double precision, db_m_d2000 double precision, db_m_d4000 double precision, db_m_d8000 double precision, 
				   db_m_e63 double precision, db_m_e125 double precision,db_m_e250 double precision,db_m_e500 double precision, 
				   db_m_e1000 double precision,db_m_e2000 double precision, db_m_e4000 double precision, db_m_e8000 double precision,
				   db_m_n63 double precision, db_m_n125 double precision, db_m_n250 double precision,db_m_n500 double precision, 
				   db_m_n1000 double precision, db_m_n2000 double precision, db_m_n4000 double precision, db_m_n8000 double precision);

drop table if exists receiver_lvl_day_zone, receiver_lvl_evening_zone, receiver_lvl_night_zone;

create table receiver_lvl_day_zone (idrecepteur integer, idsource integer, 
					  att63 double precision, att125 double precision, att250 double precision, att500 double precision, 
					  att1000 double precision, att2000 double precision, att4000 double precision, att8000 double precision);

create table receiver_lvl_evening_zone as select * from receiver_lvl_day_zone;

create table receiver_lvl_night_zone as select * from receiver_lvl_day_zone;
