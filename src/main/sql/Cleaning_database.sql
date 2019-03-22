------------------------------
-- CLEANING OF THE DATABASE --
------------------------------
drop table if exists
    bb, buildings, buildings_2d_zone, buildings_simp_merge, buildings_simp_merge_0, buildings_zone, dem1,
    dem2, dem2_points, ground, indexed_points, land_use_simp_merge, land_use_zone, receiver_lvl_day_zone,
    receiver_lvl_evening_zone, receiver_lvl_night_zone, receivers, receivers_build, receivers_build_ratio,
    roads_proj, roads_simp_merge_0, roads_src_zone, roads_traffic_zone, roads_traffic_zone_format,
    roads_zone, zone;
-- Removed from the list : cerema_ppbe_lorient_absorption, cerema_ppbe_lorient_batiment_pop, cerema_ppbe_lorient_crb_niv_1m, cerema_ppbe_lorient_route, zone_2km
