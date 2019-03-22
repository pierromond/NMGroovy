alter table receiver_lvl_day_zone alter column idrecepteur set not null;
alter table receiver_lvl_day_zone alter column idsource set not null;
alter table receiver_lvl_day_zone add primary key (idrecepteur, idsource);


alter table receiver_lvl_evening_zone alter column idrecepteur set not null;
alter table receiver_lvl_evening_zone alter column idsource set not null;
alter table receiver_lvl_evening_zone add primary key (idrecepteur, idsource);


alter table receiver_lvl_night_zone alter column idrecepteur set not null;
alter table receiver_lvl_night_zone alter column idsource set not null;
alter table receiver_lvl_night_zone add primary key (idrecepteur, idsource);
