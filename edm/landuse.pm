#!/usr/bin/perl -wT

##############################################################################
# PURPOSE: pm - landuse coverage metadata used by landuseserver.
# NOTES:   Usage:
#          use strict;
#          BEGIN { push @INC, '/data/Apache/cgi-bin'; }
#          use landuse; # Use my $coverage_metadata_reference = landuse::get_coverage_metadata_reference();
#          pop @INC; # Remove local directory from path.
#          package main;
#          ...
# HISTORY: 2026-06-209 plessel.todd@epa.gov
# STATUS:  unreviewed tested
##############################################################################

use strict;
package landuse;

##############################################################################

my @coverage_metadata = (

  {
    name        => 'sediment_gridded_mud',
    units       => '%',
    description => 'Gulf Sediment Model gridded sediment mud'
  },
  {
    name        => 'sediment_gridded_mudu',
    units       => '%',
    description => 'Gulf Sediment Model gridded sediment mudu'
  },
  {
    name        => 'sediment_gridded_sand',
    units       => '%',
    description => 'Gulf Sediment Model gridded sediment sand'
  },
  {
    name        => 'sediment_gridded_gravel',
    units       => '%',
    description => 'Gulf Sediment Model gridded sediment gravel'
  },
  {
    name        => 'sediment_gridded_rock',
    units       => '%',
    description => 'Gulf Sediment Model gridded sediment rock'
  },
  {
    name        => 'nlcd_2011',
    units       => '-',
    description => 'USGS National Land Coverage Dataset land use 2011 categories'
  },
  {
    name        => 'nlcd_2006',
    units       => '-',
    description => 'USGS National Land Coverage Dataset land use 2006 categories'
  },
  {
    name        => 'nlcd_2001',
    units       => '-',
    description => 'USGS National Land Coverage Dataset land use 2001 categories'
  },
  {
    name        => 'nlcd_1992',
    units       => '-',
    description => 'USGS National Land Coverage Dataset land use 1992 categories'
  },
  {
    name        => 'nlcd_change_1992_2001',
    units       => '-',
    description => 'USGS National Land Coverage Dataset changed land use categories 1992-2001'
  },
  {
    name        => 'nlcd_change_2001_2006',
    units       => '-',
    description => 'USGS National Land Coverage Dataset changed land use categories 2001-2006'
  },
  {
    name        => 'ccap_2006',
    units       => '-',
    description => 'NOAA Coverage Change land use 2006 categories'
  },
  {
    name        => 'ccap_2001',
    units       => '-',
    description => 'NOAA Coverage Change land use 2001 categories'
  },
  {
    name        => 'decadal_iclus_land_use_2000_2100',
    units       => '-',
    description => 'ICLUSv2 100m gridded middle-of-the-road scenario decadal land use categories 2000-2100'
  },
  {
    name        => 'impervious_2006',
    units       => '%',
    description => 'USGS percent imperervious land surface 2006'
  },
  {
    name        => 'impervious_2001',
    units       => '%',
    description => 'USGS percent impervious land surface 2001'
  },
  {
    name        => 'impervious_change_2001_2006',
    units       => '%',
    description => 'USGS percent impervious land surface change 2001-2006'
  },
  {
    name        => 'decadal_iclus_a1_impervious_2010_2100',
    units       => '%',
    description => 'ICLUS Scenario A1 decadal impervious land surface 2010-2100'
  },
  {
    name        => 'decadal_iclus_a2_impervious_2010_2100',
    units       => '%',
    description => 'ICLUS Scenario A2 decadal impervious land surface 2010-2100'
  },
  {
    name        => 'decadal_iclus_b1_impervious_2010_2100',
    units       => '%',
    description => 'ICLUS Scenario B1 decadal impervious land surface 2010-2100'
  },
  {
    name        => 'decadal_iclus_b2_impervious_2010_2100',
    units       => '%',
    description => 'ICLUS Scenario B2 decadal impervious land surface 2010-2100'
  },
  {
    name        => 'decadal_iclus_bc_impervious_2010_2100',
    units       => '%',
    description => 'ICLUS Scenario BC decadal impervious land surface 2010-2100'
  },
  {
    name        => 'elevation',
    units       => 'm',
    description => 'NOAA topography/bathymetry elevation'
  },
  {
    name        => 'elevation_adak_ak_mhw',
    units       => 'm',
    description => 'NOAA Tsunami Inundation topography/bathymetry elevation for Adak, AK (MHW)'
  },
  {
    name        => 'elevation_agua_hedionda_creek_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Agua Hedionda Creek, CA (MHW)'
  },
  {
    name        => 'elevation_akutan_ak_mhhw',
    units       => 'm',
    description => 'NOAA Tsunami Inundation topography/bathymetry elevation for Akutan, AK (MHW)'
  },
  {
    name        => 'elevation_albemarle_sound_nc_mllw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Albemarle Sound NC (MLLW)'
  },
  {
    name        => 'elevation_albion_river_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Albion River, CA (MHW)'
  },
  {
    name        => 'elevation_alder_creek_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Alder Creek, CA (MHW)'
  },
  {
    name        => 'elevation_aliso_canyon_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Aliso Canyon, CA (MHW)'
  },
  {
    name        => 'elevation_aliso_creek_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Aliso Creek, CA (MHW)'
  },
  {
    name        => 'elevation_alligator_cove_fl_navd88',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Alligator Cove, FL (NAVD88)'
  },
  {
    name        => 'elevation_alsea_bay_or_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Alsea Bay, OR (MHW)'
  },
  {
    name        => 'elevation_altahama_sound_ga_navd88',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Altahama Sound, GA (NAVD88)'
  },
  {
    name        => 'elevation_alva_paul_creek_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Alva Paul Creek, CA (MHW)'
  },
  {
    name        => 'elevation_andree_clark_bird_refuge_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Andree Clark Bird Refuge, CA (MHW)'
  },
  {
    name        => 'elevation_annisquam_river_ma_mllw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Annisquam River, MA (MLLW)'
  },
  {
    name        => 'elevation_apalachee_bay_fl_mllw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Apalachicola Bay, FL (MLLW)'
  },
  {
    name        => 'elevation_apalachicola_bay_fl_mllw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Apalachee Bay, FL (MLLW)'
  },
  {
    name        => 'elevation_aransas_bay_1_tx_mllw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Aransas Bay 1, TX (MLLW)'
  },
  {
    name        => 'elevation_aransas_bay_2_tx_mllw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Aransas Bay 2, TX (MLLW)'
  },
  {
    name        => 'elevation_arena_cove_ca_mhw',
    units       => 'm',
    description => 'NOAA Tsunami Inundation topography/bathymetry elevation for Arena Cove, CA (MHW)'
  },
  {
    name        => 'elevation_arroyo_burro_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Arroyo Burro, CA (MHW)'
  },
  {
    name        => 'elevation_arroyo_de_la_cruz_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Arroyo De La Cruz, CA (MHW)'
  },
  {
    name        => 'elevation_arroyo_del_corral_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Arroyo Del Corral, CA (MHW)'
  },
  {
    name        => 'elevation_arroyo_del_puerto_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Arroyo Del Puerto, CA (MHW)'
  },
  {
    name        => 'elevation_arroyo_grande_creek_lagoon_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Arroyo Grande Creek Lagoon, CA (MHW)'
  },
  {
    name        => 'elevation_arroyo_laguna_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Arroyo Laguna, CA (MHW)'
  },
  {
    name        => 'elevation_arthur_kill_nj_mllw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Arthur Kill, NJ (MLLW)'
  },
  {
    name        => 'elevation_arthur_kill_ny_mllw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Arthur Kill, NY (MLLW)'
  },
  {
    name        => 'elevation_assawoman_bay_de_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Assawoman Bay, DE (MHW)'
  },
  {
    name        => 'elevation_assawoman_bay_md_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Assawoman Bay, MD (MHW)'
  },
  {
    name        => 'elevation_astoria_or_mhw',
    units       => 'm',
    description => 'NOAA Tsunami Inundation topography/bathymetry elevation for Astoria, OR (MHW)'
  },
  {
    name        => 'elevation_atka_ak_mhw',
    units       => 'm',
    description => 'NOAA Tsunami Inundation topography/bathymetry elevation for Atka, AK (MHW)'
  },
  {
    name        => 'elevation_atlantic_city_nj_mhw',
    units       => 'm',
    description => 'NOAA Tsunami Inundation topography/bathymetry elevation for Atlantic City, NJ (MHW)'
  },
  {
    name        => 'elevation_baffin_bay_tx_mllw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Baffin Bay, TX (MLLW)'
  },
  {
    name        => 'elevation_baldwin_creek_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Baldwin Creek, CA (MHW)'
  },
  {
    name        => 'elevation_bar_harbor_me_mhw',
    units       => 'm',
    description => 'NOAA Tsunami Inundation topography/bathymetry elevation for Bar Harbor, ME (MHW)'
  },
  {
    name        => 'elevation_bar_harbor_me_navd88',
    units       => 'm',
    description => 'NOAA Tsunami Inundation topography/bathymetry elevation for Bar Harbor, ME (NAVD88)'
  },
  {
    name        => 'elevation_barnegat_bay_nj_mllw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Barnegat Bay, NJ (MLLW)'
  },
  {
    name        => 'elevation_barnstable_harbor_ma_navd88',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Barnstable Bay, MA (NAVD88)'
  },
  {
    name        => 'elevation_bc_msl',
    units       => 'm',
    description => 'NOAA Tsunami Inundation topography/bathymetry elevation for British Columbia (MSL)'
  },
  {
    name        => 'elevation_bear_river_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Bear River, CA (MHW)'
  },
  {
    name        => 'elevation_beaver_creek_or_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Beaver Creek, OR (MHW)'
  },
  {
    name        => 'elevation_bell_canyon_creek_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Bell Canyon Creek, CA (MHW)'
  },
  {
    name        => 'elevation_bellingham_bay_wa_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Bellingham Bay, WA (MHW)'
  },
  {
    name        => 'elevation_berry_creek_or_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Berry Creek, OR (MHW)'
  },
  {
    name        => 'elevation_big_creek_1_or_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Big Creek 1, OR (MHW)'
  },
  {
    name        => 'elevation_big_creek_or_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Big Creek, OR (MHW)'
  },
  {
    name        => 'elevation_big_lagoon_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Big Lagoon, CA (MHW)'
  },
  {
    name        => 'elevation_big_river_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Big River, CA (MHW)'
  },
  {
    name        => 'elevation_big_spring_creek_fl_mllw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Big Spring Creek, FL (MLLW)'
  },
  {
    name        => 'elevation_big_sur_river_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Big Sur River, CA (MHW)'
  },
  {
    name        => 'elevation_biloxi_ms_mhw',
    units       => 'm',
    description => 'NOAA Tsunami Inundation topography/bathymetry elevation for Biloxi, MS (MHW)'
  },
  {
    name        => 'elevation_birch_bay_wa_mllw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Birch Bay, WA (MLLW)'
  },
  {
    name        => 'elevation_biscayne_bay_fl_navd88',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Biscayne Bay, FL (NAVD88)'
  },
  {
    name        => 'elevation_blackwater_sound_fl_navd88',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Blackwater Sound, FL (NAVD88)'
  },
  {
    name        => 'elevation_blue_hill_bay_me_navd88',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Blue Hill Bay, ME (NAVD88)'
  },
  {
    name        => 'elevation_bob_creek_or_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Bob Creek, OR (MHW)'
  },
  {
    name        => 'elevation_bodega_bay_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Bodega Bay, CA (MHW)'
  },
  {
    name        => 'elevation_boggy_bay_fl_mllw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Boggy Bay, FL (MLLW)'
  },
  {
    name        => 'elevation_boggy_creek_fl_mllw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Boggy Creek, FL (MLLW)'
  },
  {
    name        => 'elevation_bogue_sound_nc_navd88',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Bogue Sound, NC (NAVD88)'
  },
  {
    name        => 'elevation_bolinas_lagoon_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Bolinas Lagoon, CA (MHW)'
  },
  {
    name        => 'elevation_bolsa_chica_ecological_reserve_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Bolsa Chica Ecological Reserve, CA (MHW)'
  },
  {
    name        => 'elevation_boone_creek_wa_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Boone Creek, WA (MHW)'
  },
  {
    name        => 'elevation_boston_harbor_ma_mllw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Boston Harbor, MA (MLLW)'
  },
  {
    name        => 'elevation_brave_boat_harbor_me_navd88',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Brave Boat Harbor, ME (NAVD88)'
  },
  {
    name        => 'elevation_brazos_river_tx_navd88',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Brazos River, TX (NAVD88)'
  },
  {
    name        => 'elevation_broad_river_rodgers_river_fl_navd88',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Broad River/Port Royal Sound, SC (NAVD88)'
  },
  {
    name        => 'elevation_broad_river_sc_navd88',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Broad River/Rodgers River, FL (NAVD88)'
  },
  {
    name        => 'elevation_broad_sound_ma_mllw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Broad Sound, MA (MLLW)'
  },
  {
    name        => 'elevation_brush_creek_1_or_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Brush Creek 1, OR (MHW)'
  },
  {
    name        => 'elevation_brush_creek_2_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Brush Creek 2, CA (MHW)'
  },
  {
    name        => 'elevation_buck_creek_fl_mllw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Buck Creek, FL (MLLW)'
  },
  {
    name        => 'elevation_buena_vista_lagoon_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Buena Vista Lagoon, CA (MHW)'
  },
  {
    name        => 'elevation_bulls_bay_sc_navd88',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Bulls Bay, SC (NAVD88)'
  },
  {
    name        => 'elevation_bumblebee_creek_fl_mllw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Bumblebee Creek, FL (MLLW)'
  },
  {
    name        => 'elevation_bush_river_md_navd88',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Bush River, MD (NAVD88)'
  },
  {
    name        => 'elevation_buzzards_bay_ma_mllw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Buzzards Bay, MA (MLLW)'
  },
  {
    name        => 'elevation_calera_creek_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Calera Creek, CA (MHW)'
  },
  {
    name        => 'elevation_camp_creek_wa_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Camp Creek, WA (MHW)'
  },
  {
    name        => 'elevation_canada_de_la_gaviota_creek_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Canada De La Gaviota Creek, CA (MHW)'
  },
  {
    name        => 'elevation_canada_del_refugio_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Canada Del Refugio, CA (MHW)'
  },
  {
    name        => 'elevation_cape_creek_or_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Cape Creek, OR (MHW)'
  },
  {
    name        => 'elevation_cape_fear_river_nc_navd88',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Cape Fear River, NC (NAVD88)'
  },
  {
    name        => 'elevation_cape_hatteras_nc_mhw',
    units       => 'm',
    description => 'NOAA Tsunami Inundation topography/bathymetry elevation for Cape Hatteras, NC (MHW)'
  },
  {
    name        => 'elevation_cape_romaine_refuge_sc_navd88',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Cape Romaine Refuge, SC (NAVD88)'
  },
  {
    name        => 'elevation_capers_creek_sc_navd88',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Capers Creek, SC (NAVD88)'
  },
  {
    name        => 'elevation_carmel_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Carmel, CA (MHW)'
  },
  {
    name        => 'elevation_carpinteria_creek_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Carpinteria Creek, CA (MHW)'
  },
  {
    name        => 'elevation_carpinteria_salt_marsh_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Carpinteria Salt Marsh, CA (MHW)'
  },
  {
    name        => 'elevation_casco_bay_me_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Casco Bay, ME (MHW)'
  },
  {
    name        => 'elevation_casper_creek_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Casper Creek, CA (MHW)'
  },
  {
    name        => 'elevation_cayucos_creek_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Cayucos Creek, CA (MHW)'
  },
  {
    name        => 'elevation_cedar_creek_wa_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Cedar Creek, WA (MHW)'
  },
  {
    name        => 'elevation_central_ca_mhw',
    units       => 'm',
    description => 'NOAA Tsunami Inundation topography/bathymetry elevation for Central CA (MHW)'
  },
  {
    name        => 'elevation_central_ca_navd88',
    units       => 'm',
    description => 'NOAA Tsunami Inundation topography/bathymetry elevation for Central CA (NAVD88)'
  },
  {
    name        => 'elevation_central_fl_mhw',
    units       => 'm',
    description => 'NOAA Tsunami Inundation topography/bathymetry elevation for Central FL (MHW)'
  },
  {
    name        => 'elevation_central_fl_navd88',
    units       => 'm',
    description => 'NOAA Tsunami Inundation topography/bathymetry elevation for Central FL (NAVD88)'
  },
  {
    name        => 'elevation_central_gulf_msl',
    units       => 'm',
    description => 'NOAA Tsunami Inundation topography/bathymetry elevation for Central Gulf (MSL)'
  },
  {
    name        => 'elevation_central_or_mhw',
    units       => 'm',
    description => 'NOAA Tsunami Inundation topography/bathymetry elevation for Central OR (MHW)'
  },
  {
    name        => 'elevation_central_pacific_msl',
    units       => 'm',
    description => 'NOAA Tsunami Inundation topography/bathymetry elevation for Central Pacific (MSL)'
  },
  {
    name        => 'elevation_chadbourne_gulch_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Chadbourne Gulch, CA (MHW)'
  },
  {
    name        => 'elevation_channel_islands_harbor_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Channel Islands Harbor, CA (MHW)'
  },
  {
    name        => 'elevation_charleston_harbor_sc_navd88',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Charleston Harbor, SC (NAVD88)'
  },
  {
    name        => 'elevation_charlotte_harbor_fl_mllw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Charlotte Harbor, FL (MLLW)'
  },
  {
    name        => 'elevation_chenega_ak_mhhw',
    units       => 'm',
    description => 'NOAA Tsunami Inundation topography/bathymetry elevation for Chenega, AK (MHHW)'
  },
  {
    name        => 'elevation_chesapeake_bay_mainstem_1_md_mllw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Chesapeake Bay Mainstem 1, MD (MLLW)'
  },
  {
    name        => 'elevation_chesapeake_bay_mainstem_1_va_mllw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Chesapeake Bay Mainstem 1, VA (MLLW)'
  },
  {
    name        => 'elevation_chesapeake_bay_mainstem_2_md_mllw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Chesapeake Bay Mainstem 2, MD (MLLW)'
  },
  {
    name        => 'elevation_chesapeake_bay_mainstem_2_va_mllw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Chesapeake Bay Mainstem 2, VA (MLLW)'
  },
  {
    name        => 'elevation_chester_river_md_navd88',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Chester River, MD (NAVD88)'
  },
  {
    name        => 'elevation_chetco_river_or_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Chetco River, OR (MHW)'
  },
  {
    name        => 'elevation_chignik_ak_mhw',
    units       => 'm',
    description => 'NOAA Tsunami Inundation topography/bathymetry elevation for Chignik, AK (MHW)'
  },
  {
    name        => 'elevation_china_creek_or_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for China Creek, OR (MHW)'
  },
  {
    name        => 'elevation_chincoteague_bay_md_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Chincoteague Bay, MD (MHW)'
  },
  {
    name        => 'elevation_chincoteague_bay_va_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Chincoteague Bay, VA (MHW)'
  },
  {
    name        => 'elevation_choctawhatchee_bay_fl_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Choctawhatchee Bay, FL (MHW)'
  },
  {
    name        => 'elevation_choptank_river_md_navd88',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Choptank River, MD (NAVD88)'
  },
  {
    name        => 'elevation_clam_beach_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Clam Beach, CA (MHW)'
  },
  {
    name        => 'elevation_clay_creek_fl_mllw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Clay Creek, FL (MLLW)'
  },
  {
    name        => 'elevation_clearwater_saint_joseph_sound_fl_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Clearwater/St. Joseph Sound, FL (MHW)'
  },
  {
    name        => 'elevation_cobscook_me_mllw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Cobscook, ME (MLLW)'
  },
  {
    name        => 'elevation_cockelburr_canyon_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Cockelburr Canyon, CA (MHW)'
  },
  {
    name        => 'elevation_cocohatchee_river_fl_navd88',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Cocohatchee River, FL (NAVD88)'
  },
  {
    name        => 'elevation_columbia_river_or_mllw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Columbia River, OR (MLLW)'
  },
  {
    name        => 'elevation_columbia_river_wa_mllw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Columbia River, WA (MLLW)'
  },
  {
    name        => 'elevation_conner_creek_wa_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Conner Creek, WA (MHW)'
  },
  {
    name        => 'elevation_coos_bay_or_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Coos Bay, OR (MHW)'
  },
  {
    name        => 'elevation_copalis_river_wa_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Copalis River, WA (MHW)'
  },
  {
    name        => 'elevation_coquille_river_or_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Coquille River, OR (MHW)'
  },
  {
    name        => 'elevation_corpus_christi_bay_tx_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Corpus Christi Bay, TX (MHW)'
  },
  {
    name        => 'elevation_corpus_christi_tx_mhw',
    units       => 'm',
    description => 'NOAA Tsunami Inundation topography/bathymetry elevation for Corpus Christi, TX (MHW)'
  },
  {
    name        => 'elevation_cottaneva_creek_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Cottaneva Creek, CA (MHW)'
  },
  {
    name        => 'elevation_craig_ak_mhw',
    units       => 'm',
    description => 'NOAA Tsunami Inundation topography/bathymetry elevation for Craig, AK (MHW)'
  },
  {
    name        => 'elevation_crescent_city_ca_mhw',
    units       => 'm',
    description => 'NOAA Tsunami Inundation topography/bathymetry elevation for Crescent City, CA (MHW)'
  },
  {
    name        => 'elevation_crescent_city_ca_navd88',
    units       => 'm',
    description => 'NOAA Tsunami Inundation topography/bathymetry elevation for Crescent City, CA (NAVD88)'
  },
  {
    name        => 'elevation_crescent_city_harbor_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Crescent City Harbor, CA (MHW)'
  },
  {
    name        => 'elevation_crooked_creek_fl_mllw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Crooked Creek, FL (MLLW)'
  },
  {
    name        => 'elevation_crystal_river_1_fl_mllw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Crystal River 1, FL (MLLW)'
  },
  {
    name        => 'elevation_crystal_river_2_fl_mllw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Crystal River 2, FL (MLLW)'
  },
  {
    name        => 'elevation_crystal_river_3_fl_mllw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Crystal River 3, FL (MLLW)'
  },
  {
    name        => 'elevation_cumberland_sound_fl_navd88',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Cumberland Sound/Saint Marys River, FL (NAVD88)'
  },
  {
    name        => 'elevation_cumberland_sound_ga_navd88',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Cumberland Sound/Saint Marys River, GA (NAVD88)'
  },
  {
    name        => 'elevation_cummins_creek_or_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Cummins Creek, OR (MHW)'
  },
  {
    name        => 'elevation_currituck_sound_nc_navd88',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Currituck Sound, NC (NAVD88)'
  },
  {
    name        => 'elevation_currituck_sound_va_navd88',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Currituck Sound, VA (NAVD88)'
  },
  {
    name        => 'elevation_cut_creek_or_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Cut Creek, OR (MHW)'
  },
  {
    name        => 'elevation_dairy_gulch_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Dairy Gulch, CA (MHW)'
  },
  {
    name        => 'elevation_daley_lake_or_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Daley Lake, OR (MHW)'
  },
  {
    name        => 'elevation_dallus_creek_fl_mllw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Dallus Creek, FL (MLLW)'
  },
  {
    name        => 'elevation_damariscaotta_me_mllw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Damariscaotta, ME (MLLW)'
  },
  {
    name        => 'elevation_dana_point_harbor_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Dana Point Harbor, CA (MHW)'
  },
  {
    name        => 'elevation_days_creek_fl_mllw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Days Creek, FL (MLLW)'
  },
  {
    name        => 'elevation_daytona_beach_fl_mhw',
    units       => 'm',
    description => 'NOAA Tsunami Inundation topography/bathymetry elevation for Daytona Beach, FL (MHW)'
  },
  {
    name        => 'elevation_deadman_bay_fl_mllw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Deadman Bay, FL (MLLW)'
  },
  {
    name        => 'elevation_deboy_sound_ga_navd88',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Deboy Sound, GA (NAVD88)'
  },
  {
    name        => 'elevation_delaware_bay_1_de_mllw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Delaware Bay 1, DE (MLLW)'
  },
  {
    name        => 'elevation_delaware_bay_1_nj_mllw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Delaware Bay 1, NJ (MLLW)'
  },
  {
    name        => 'elevation_delaware_bay_2_de_mllw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Delaware Bay 2, DE (MLLW)'
  },
  {
    name        => 'elevation_delaware_bay_2_nj_mllw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Delaware Bay 2, NJ (MLLW)'
  },
  {
    name        => 'elevation_delaware_bay_3_de_mllw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Delaware Bay 3, DE (MLLW)'
  },
  {
    name        => 'elevation_delaware_bay_3_nj_mllw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Delaware Bay 3, NJ (MLLW)'
  },
  {
    name        => 'elevation_depoe_bay_or_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Depoe Bay, OR (MHW)'
  },
  {
    name        => 'elevation_devereux_lagoon_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Devereux Lagoon, CA (MHW)'
  },
  {
    name        => 'elevation_devils_lake_or_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Devils Lake, OR (MHW)'
  },
  {
    name        => 'elevation_dewes_hamlin_creek_sc_navd88',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Dewes/Hamlin Creek, SC (NAVD88)'
  },
  {
    name        => 'elevation_drakes_estero_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Drakes Estero, CA (MHW)'
  },
  {
    name        => 'elevation_drayton_harbor_wa_mllw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Drayton Harbor, WA (MLLW)'
  },
  {
    name        => 'elevation_duck_creek_wa_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Duck Creek, WA (MHW)'
  },
  {
    name        => 'elevation_dutch_harbor_ak_mhw',
    units       => 'm',
    description => 'NOAA Tsunami Inundation topography/bathymetry elevation for Dutch Harbor, AK (MHW)'
  },
  {
    name        => 'elevation_duxbury_bay_ma_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Duxbury Bay/Plymouth Harbor, MA (MHW)'
  },
  {
    name        => 'elevation_eaglenest_creek_fl_mllw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Eaglenest Creek, FL (MLLW)'
  },
  {
    name        => 'elevation_eastern_bay_md_navd88',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Eastern Bay, MD (NAVD88)'
  },
  {
    name        => 'elevation_eastern_gulf_msl',
    units       => 'm',
    description => 'NOAA Tsunami Inundation topography/bathymetry elevation for Eastern Gulf (MSL)'
  },
  {
    name        => 'elevation_ecola_creek_or_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Ecola Creek, OR (MHW)'
  },
  {
    name        => 'elevation_eel_river_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Eel River, CA (MHW)'
  },
  {
    name        => 'elevation_eleven_prong_bay_fl_mllw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Eleven Prong Bay, FL (MLLW)'
  },
  {
    name        => 'elevation_elfin_cove_ak_mhw',
    units       => 'm',
    description => 'NOAA Tsunami Inundation topography/bathymetry elevation for Elfin Cove, AK (MHW)'
  },
  {
    name        => 'elevation_elk_creek_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Elk Creek, CA (MHW)'
  },
  {
    name        => 'elevation_elk_creek_wa_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Elk Creek, WA (MHW)'
  },
  {
    name        => 'elevation_elk_river_2_wa_mllw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Elk River 2, WA (MLLW)'
  },
  {
    name        => 'elevation_elk_river_3_or_mllw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Elk River 3, OR (MLLW)'
  },
  {
    name        => 'elevation_elk_river_md_navd88',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Elk River, MD (NAVD88)'
  },
  {
    name        => 'elevation_espa_lagoon_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Espa Lagoon, CA (MHW)'
  },
  {
    name        => 'elevation_essex_bay_ma_mllw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Essex Bay, MA (MLLW)'
  },
  {
    name        => 'elevation_estero_bay_fl_navd88',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Estero Bay, FL (NAVD88)'
  },
  {
    name        => 'elevation_euchre_creek_or_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Euchre Creek, OR (MHW)'
  },
  {
    name        => 'elevation_eureka_ca_mhw',
    units       => 'm',
    description => 'NOAA Tsunami Inundation topography/bathymetry elevation for Eureka, CA (MHW)'
  },
  {
    name        => 'elevation_eureka_ca_navd88',
    units       => 'm',
    description => 'NOAA Tsunami Inundation topography/bathymetry elevation for Eureka, CA (NAVD88)'
  },
  {
    name        => 'elevation_fern_canyon_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Fern Canyon, CA (MHW)'
  },
  {
    name        => 'elevation_fillman_bayou_fl_mllw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Fillman Bayou, FL (MLLW)'
  },
  {
    name        => 'elevation_first_bay_fl_navd88',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for First Bay, FL (NAVD88)'
  },
  {
    name        => 'elevation_fish_creek_fl_mllw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Fish Creek, FL (MLLW)'
  },
  {
    name        => 'elevation_floras_creek_2_or_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Floras Creek 2, OR (MHW)'
  },
  {
    name        => 'elevation_florida_bay_fl_mllw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Florida Bay, FL (MLLW)'
  },
  {
    name        => 'elevation_fort_bragg_ca_mhw',
    units       => 'm',
    description => 'NOAA Tsunami Inundation topography/bathymetry elevation for Fort Bragg, CA (MHW)'
  },
  {
    name        => 'elevation_fort_bragg_ca_navd88',
    units       => 'm',
    description => 'NOAA Tsunami Inundation topography/bathymetry elevation for Fort Bragg, CA (NAVD88)'
  },
  {
    name        => 'elevation_french_canyon_2_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for French Canyon 2, CA (MHW)'
  },
  {
    name        => 'elevation_frenchmans_bay_me_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Frenchman\'s Bay, ME (MHW)'
  },
  {
    name        => 'elevation_frenchmans_creek_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Frenchman\'s Creek, CA (MHW)'
  },
  {
    name        => 'elevation_galvaston_bay_1_tx_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Galvaston Bay 1, TX (MHW)'
  },
  {
    name        => 'elevation_galvaston_bay_2_tx_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Galvaston Bay 2, TX (MHW)'
  },
  {
    name        => 'elevation_galveston_tx_mhw',
    units       => 'm',
    description => 'NOAA Tsunami Inundation topography/bathymetry elevation for Galveston, TX (MHW)'
  },
  {
    name        => 'elevation_garcia_river_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Garcia River, CA (MHW)'
  },
  {
    name        => 'elevation_gardiner_bay_peconic_bay_1_ny_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Gardiner Bay/Peconic Bay 1, NY (MHW)'
  },
  {
    name        => 'elevation_gardiner_bay_peconic_bay_2_ny_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Gardiner Bay/Peconic Bay 2, NY (MHW)'
  },
  {
    name        => 'elevation_garibaldi_or_mhw',
    units       => 'm',
    description => 'NOAA Tsunami Inundation topography/bathymetry elevation for Garibaldi, OR (MHW)'
  },
  {
    name        => 'elevation_garibaldi_or_navd88',
    units       => 'm',
    description => 'NOAA Tsunami Inundation topography/bathymetry elevation for Garibaldi, OR (NAVD88)'
  },
  {
    name        => 'elevation_garrapata_creek_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Garrapata Creek, CA (MHW)'
  },
  {
    name        => 'elevation_gazos_creek_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Gazos Creek, CA (MHW)'
  },
  {
    name        => 'elevation_glouchester_harbor_ma_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Glouchester Harbor, MA, (MHW)'
  },
  {
    name        => 'elevation_goleta_slough_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Goleta Slough, CA (MHW)'
  },
  {
    name        => 'elevation_goodman_creek_wa_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Goodman Creek, WA (MHW)'
  },
  {
    name        => 'elevation_goosefare_bay_me_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Goosefare Bay, ME (MHW)'
  },
  {
    name        => 'elevation_gouldsboro_bay_me_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Gouldsboro Bay, ME (MHW)'
  },
  {
    name        => 'elevation_grays_harbor_wa_mllw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Grays Harbor, WA (MLLW)'
  },
  {
    name        => 'elevation_great_bay_me_mllw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Great Bay, ME (MLLW)'
  },
  {
    name        => 'elevation_great_bay_nh_mllw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Great Bay, NH (MLLW)'
  },
  {
    name        => 'elevation_great_bay_nj_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Great Bay, NJ (MHW)'
  },
  {
    name        => 'elevation_great_south_bay_ny_mllw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Great South Bay, NY (MLLW)'
  },
  {
    name        => 'elevation_great_wilmico_fleets_bay_va_navd88',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Great Wilmico/Fleets Bay, VA (NAVD88)'
  },
  {
    name        => 'elevation_greenwood_creek_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Greenwood Creek, CA (MHW)'
  },
  {
    name        => 'elevation_gualala_river_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Gualala River, CA (MHW)'
  },
  {
    name        => 'elevation_gunpowder_river_md_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Gunpowder River, MD (MHW)'
  },
  {
    name        => 'elevation_guthrie_creek_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Guthrie Creek, CA (MHW)'
  },
  {
    name        => 'elevation_half_moon_bay_airport_lagoon_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Half Moon Bay Airport Lagoon, CA (MHW)'
  },
  {
    name        => 'elevation_hampton_harbor_nh_mllw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Hampton Harbor, NH (MLLW)'
  },
  {
    name        => 'elevation_hanalei_hi_mhw',
    units       => 'm',
    description => 'NOAA Tsunami Inundation topography/bathymetry elevation for Hanalei, HI (MHW)'
  },
  {
    name        => 'elevation_hare_creek_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Hare Creek, CA (MHW)'
  },
  {
    name        => 'elevation_haycock_harbor_baileys_mistake_me_mllw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Haycock Harbor/Baileys Mistake, ME (MLLW)'
  },
  {
    name        => 'elevation_hi_msl',
    units       => 'm',
    description => 'NOAA Tsunami Inundation topography/bathymetry elevation for HI (MSL)'
  },
  {
    name        => 'elevation_hilo_hi_mhw',
    units       => 'm',
    description => 'NOAA Tsunami Inundation topography/bathymetry elevation for Hilo, HI (MHW)'
  },
  {
    name        => 'elevation_hoh_river_wa_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Hoh River, WA (MHW)'
  },
  {
    name        => 'elevation_honda_creek_lagoon_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Honda Creek Lagoon, CA (MHW)'
  },
  {
    name        => 'elevation_honga_river_md_navd88',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Honga River, MD (NAVD88)'
  },
  {
    name        => 'elevation_horseshoe_cove_fl_mllw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Horseshoe Cove, FL (MLLW)'
  },
  {
    name        => 'elevation_hubbard_creek_or_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Hubbard Creek, OR (MHW)'
  },
  {
    name        => 'elevation_hudson_river_nj_mllw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Hudson River, NJ (MLLW)'
  },
  {
    name        => 'elevation_hudson_river_ny_mllw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Hudson River, NY (MLLW)'
  },
  {
    name        => 'elevation_humboldt_bay_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Humboldt Bay, CA (MHW)'
  },
  {
    name        => 'elevation_humboldt_ca_mhw',
    units       => 'm',
    description => 'NOAA Tsunami Inundation topography/bathymetry elevation for Humboldt, CA (MHW)'
  },
  {
    name        => 'elevation_huntington_channel_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Huntington Channel, CA (MHW)'
  },
  {
    name        => 'elevation_icw4_fl_navd88',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Icw4, FL (NAVD88)'
  },
  {
    name        => 'elevation_icw_wolf_mobile_al_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Icw Wolf Mobile, AL (MHW)'
  },
  {
    name        => 'elevation_indian_bay_fl_mllw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Indian Bay, FL (MLLW)'
  },
  {
    name        => 'elevation_indian_river_bay_de_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Indian River Bay, DE (MHW)'
  },
  {
    name        => 'elevation_indian_river_bay_fl_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Indian River Bay, FL (MHW)'
  },
  {
    name        => 'elevation_inglenook_fen_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Inglenook Fen, CA (MHW)'
  },
  {
    name        => 'elevation_irish_gulch_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Irish Gulch, CA (MHW)'
  },
  {
    name        => 'elevation_jackass_creek_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Jackass Creek, CA (MHW)'
  },
  {
    name        => 'elevation_jalama_creek_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Jalama Creek, CA (MHW)'
  },
  {
    name        => 'elevation_jamaica_bay_ny_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Jamaica Bay, NY (MHW)'
  },
  {
    name        => 'elevation_james_river_va_navd88',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for James River, VA (NAVD88)'
  },
  {
    name        => 'elevation_jekyll_sound_saint_andrews_sound_ga_navd88',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Jekyll Sound/Saint Andrews Sound, GA (NAVD88)'
  },
  {
    name        => 'elevation_joe_creek_wa_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Joe Creek, WA (MHW)'
  },
  {
    name        => 'elevation_johnson_creek_3_or_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Johnson Creek 3, OR (MHW)'
  },
  {
    name        => 'elevation_johnson_creek_fl_mllw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Johnson Creek, FL (MLLW)'
  },
  {
    name        => 'elevation_juneau_ak_mhhw',
    units       => 'm',
    description => 'NOAA Tsunami Inundation topography/bathymetry elevation for Juneau, AK (MHHW)'
  },
  {
    name        => 'elevation_kachemak_bay_ak_mhw',
    units       => 'm',
    description => 'NOAA Tsunami Inundation topography/bathymetry elevation for Kachemak Bay, AK (MHW)'
  },
  {
    name        => 'elevation_kahului_hi_mllw',
    units       => 'm',
    description => 'NOAA Tsunami Inundation topography/bathymetry elevation for Kahului, HI (MLLW)'
  },
  {
    name        => 'elevation_kailua_kona_hi_mhw',
    units       => 'm',
    description => 'NOAA Tsunami Inundation topography/bathymetry elevation for Kailua-Kona, HI (MHW)'
  },
  {
    name        => 'elevation_kalaloch_creek_wa_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Kalaloch Creek, WA (MHW)'
  },
  {
    name        => 'elevation_katama_bay_ma_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Katama Bay, MA (MHW)'
  },
  {
    name        => 'elevation_kauai_hi_mhw',
    units       => 'm',
    description => 'NOAA Tsunami Inundation topography/bathymetry elevation for Kauai, HI (MHW)'
  },
  {
    name        => 'elevation_kawaihae_hi_mhw',
    units       => 'm',
    description => 'NOAA Tsunami Inundation topography/bathymetry elevation for Kawaihae, HI (MHW)'
  },
  {
    name        => 'elevation_keauhou_hi_mhw',
    units       => 'm',
    description => 'NOAA Tsunami Inundation topography/bathymetry elevation for Keauhou, HI (MHW)'
  },
  {
    name        => 'elevation_kennebec_river_estuary_me_mllw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Kennebec River Estuary, ME (MLLW)'
  },
  {
    name        => 'elevation_key_west_fl_mhw',
    units       => 'm',
    description => 'NOAA Tsunami Inundation topography/bathymetry elevation for Key West, FL (MHW)'
  },
  {
    name        => 'elevation_key_west_fl_navd88',
    units       => 'm',
    description => 'NOAA Tsunami Inundation topography/bathymetry elevation for Key West, FL (NAVD88)'
  },
  {
    name        => 'elevation_king_cove_ak_mhw',
    units       => 'm',
    description => 'NOAA Tsunami Inundation topography/bathymetry elevation for King Cove, AK (MHW)'
  },
  {
    name        => 'elevation_klamath_river_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Klamath River, CA (MHW)'
  },
  {
    name        => 'elevation_kodiak_ak_mhw',
    units       => 'm',
    description => 'NOAA Tsunami Inundation topography/bathymetry elevation for Kodiak, AK (MHW)'
  },
  {
    name        => 'elevation_la_push_wa_mhw',
    units       => 'm',
    description => 'NOAA Tsunami Inundation topography/bathymetry elevation for La Push, WA (MHW)'
  },
  {
    name        => 'elevation_lagoon_creek_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Lagoon Creek, CA (MHW)'
  },
  {
    name        => 'elevation_laguna_creek_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Laguna Creek, CA (MHW)'
  },
  {
    name        => 'elevation_lahaina_hi_mhw',
    units       => 'm',
    description => 'NOAA Tsunami Inundation topography/bathymetry elevation for Lahaina, HI (MHW)'
  },
  {
    name        => 'elevation_lake_cleone_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Lake Cleone, CA (MHW)'
  },
  {
    name        => 'elevation_lake_earl_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Lake Earl, CA (MHW)'
  },
  {
    name        => 'elevation_lake_lucerne_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Lake Lucerne, CA (MHW)'
  },
  {
    name        => 'elevation_lake_lytle_or_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Lake Lytle, OR (MHW)'
  },
  {
    name        => 'elevation_lake_washington_wa_llw',
    units       => 'm',
    description => 'NOAA Tsunami Inundation topography/bathymetry elevation for Lake Washington, WA (LLW)'
  },
  {
    name        => 'elevation_lake_worth_fl_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Lake Worth, FL (MHW)'
  },
  {
    name        => 'elevation_lakes_reeds_bay_nj_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Lakes Reeds Bay, NJ (MHW)'
  },
  {
    name        => 'elevation_las_flores_creek_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Las Flores Creek, CA (MHW)'
  },
  {
    name        => 'elevation_las_pulgas_creek_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Las Pulgas Creek, CA (MHW)'
  },
  {
    name        => 'elevation_little_creek_1_or_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Little Creek 1, OR (MHW)'
  },
  {
    name        => 'elevation_little_creek_3_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Little Creek 3, CA (MHW)'
  },
  {
    name        => 'elevation_little_kennebec_englishman_bay_me_mllw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Little Kennebec/Englishman Bay, ME (MLLW)'
  },
  {
    name        => 'elevation_little_lagoon_al_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Little Lagoon, AL (MHW)'
  },
  {
    name        => 'elevation_little_narragansett_bay_ct_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Little Narragansett Bay, CT (MHW)'
  },
  {
    name        => 'elevation_little_pico_creek_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Little Pico Creek, CA (MHW)'
  },
  {
    name        => 'elevation_little_pine_1_bay_fl_mllw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Little Pine 1 Bay, FL (MLLW)'
  },
  {
    name        => 'elevation_little_river_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Little River, CA (MHW)'
  },
  {
    name        => 'elevation_little_river_sc_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Little River, SC (MHW)'
  },
  {
    name        => 'elevation_little_rocky_creek_fl_mllw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Little Rocky Creek, FL (MLLW)'
  },
  {
    name        => 'elevation_little_salmon_creek_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Little Salmon Creek, CA (MHW)'
  },
  {
    name        => 'elevation_little_spring_creek_fl_mllw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Little Spring Creek, FL (MLLW)'
  },
  {
    name        => 'elevation_little_sur_lagoon_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Little Sur Lagoon, CA (MHW)'
  },
  {
    name        => 'elevation_lockwoods_folly_river_nc_navd88',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Lockwoods Folly River, NC (NAVD88)'
  },
  {
    name        => 'elevation_loma_alta_marsh_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Loma Alta Marsh, CA (MHW)'
  },
  {
    name        => 'elevation_lombardi_creek_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Lombardi Creek, CA (MHW)'
  },
  {
    name        => 'elevation_long_island_sound_ct_mllw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Long Island Sound, CT (MLLW)'
  },
  {
    name        => 'elevation_long_island_sound_ny_mllw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Long Island Sound, NY (MLLW)'
  },
  {
    name        => 'elevation_los_penasquitos_lagoon_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Los Penasquitos Lagoon, CA (MHW)'
  },
  {
    name        => 'elevation_lows_bay_fl_mllw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Lows Bay, FL (MLLW)'
  },
  {
    name        => 'elevation_loxahatchee_river_fl_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Loxahatchee River, FL (MHW)'
  },
  {
    name        => 'elevation_lunada_bay_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Lunada Bay, CA (MHW)'
  },
  {
    name        => 'elevation_lynnhaven_river_va_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Lynnhaven River, VA (MHW)'
  },
  {
    name        => 'elevation_machias_bay_me_mllw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Machias Bay, ME (MLLW)'
  },
  {
    name        => 'elevation_mad_river_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Mad River, CA (MHW)'
  },
  {
    name        => 'elevation_maddaket_harbor_ma_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Maddaket Harbor, MA (MHW)'
  },
  {
    name        => 'elevation_magothy_river_md_navd88',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Magothy River, MD (NAVD88)'
  },
  {
    name        => 'elevation_makah_bay_wa_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Magu Lagoon, CA (MHW)'
  },
  {
    name        => 'elevation_malibu_lagoon_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Makah Bay, WA (MHW)'
  },
  {
    name        => 'elevation_malpaso_creek_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Malibu Lagoon, CA (MHW)'
  },
  {
    name        => 'elevation_manasquan_river_nj_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Malpaso Creek, CA (MHW)'
  },
  {
    name        => 'elevation_manchester_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Manasquan River, NJ (MHW)'
  },
  {
    name        => 'elevation_matagorda_bay_tx_mllw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Manchester, CA (MHW)'
  },
  {
    name        => 'elevation_mattole_river_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Matagorda Bay, TX (MLLW)'
  },
  {
    name        => 'elevation_mcnutt_gulch_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Mattole River, CA (MHW)'
  },
  {
    name        => 'elevation_merrimack_river_ma_mllw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Mcnutt Gulch, CA (MHW)'
  },
  {
    name        => 'elevation_miami_fl_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Merrimack River, MA (MLLW)'
  },
  {
    name        => 'elevation_miami_fl_navd88',
    units       => 'm',
    description => 'NOAA Tsunami Inundation topography/bathymetry elevation for Miami, FL (MHW)'
  },
  {
    name        => 'elevation_middle_and_black_rivers_md_mhw',
    units       => 'm',
    description => 'NOAA Tsunami Inundation topography/bathymetry elevation for Miami, FL (NAVD88)'
  },
  {
    name        => 'elevation_mississippi_sound_al_mllw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Middle and Black Rivers, MD (MHW)'
  },
  {
    name        => 'elevation_mississippi_sound_la_mllw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Mississippi Sound, AL (MLLW)'
  },
  {
    name        => 'elevation_mississippi_sound_ms_mllw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Mississippi Sound, LA (MLLW)'
  },
  {
    name        => 'elevation_misson_bay_lagoon_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Mississippi Sound, MS (MLLW)'
  },
  {
    name        => 'elevation_misson_creek_lagoon_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Misson Bay Lagoon, CA (MHW)'
  },
  {
    name        => 'elevation_mitchell_creek_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Misson Creek Lagoon, CA (MHW)'
  },
  {
    name        => 'elevation_mobile_al_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Mitchell Creek, CA (MHW)'
  },
  {
    name        => 'elevation_mobile_al_navd88',
    units       => 'm',
    description => 'NOAA Tsunami Inundation topography/bathymetry elevation for Mobile, AL (MHW)'
  },
  {
    name        => 'elevation_mobile_bay_al_mhw',
    units       => 'm',
    description => 'NOAA Tsunami Inundation topography/bathymetry elevation for Mobile, AL (NAVD88)'
  },
  {
    name        => 'elevation_moclips_river_wa_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Mobile Bay, AL (MHW)'
  },
  {
    name        => 'elevation_mojack_bay_va_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Moclips River, WA (MHW)'
  },
  {
    name        => 'elevation_montauk_ny_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Mojack Bay, VA (MHW)'
  },
  {
    name        => 'elevation_monterey_bay_ca_mhw',
    units       => 'm',
    description => 'NOAA Tsunami Inundation topography/bathymetry elevation for Montauk, NY (MHW)'
  },
  {
    name        => 'elevation_monterey_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Monterey Bay, CA (MHW)'
  },
  {
    name        => 'elevation_monterey_ca_navd88',
    units       => 'm',
    description => 'NOAA Tsunami Inundation topography/bathymetry elevation for Monterey, CA (MHW)'
  },
  {
    name        => 'elevation_moose_cove_me_mllw',
    units       => 'm',
    description => 'NOAA Tsunami Inundation topography/bathymetry elevation for Monterey, CA (NAVD88)'
  },
  {
    name        => 'elevation_moran_lake_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Moose Cove, ME (MLLW)'
  },
  {
    name        => 'elevation_morehead_city_nc_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Moran Lake, CA (MHW)'
  },
  {
    name        => 'elevation_morehead_city_nc_navd88',
    units       => 'm',
    description => 'NOAA Tsunami Inundation topography/bathymetry elevation for Morehead City, NC (MHW)'
  },
  {
    name        => 'elevation_morro_bay_ca_mhw',
    units       => 'm',
    description => 'NOAA Tsunami Inundation topography/bathymetry elevation for Morehead City, NC (NAVD88)'
  },
  {
    name        => 'elevation_morro_creek_lagoon_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Morro Bay, CA (MHW)'
  },
  {
    name        => 'elevation_mosquito_creek_wa_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Morro Creek Lagoon, CA (MHW)'
  },
  {
    name        => 'elevation_mugu_lagoon_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Mosquito Creek, WA (MHW)'
  },
  {
    name        => 'elevation_murrells_inlet_sc_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Murrells Inlet, SC (MHW)'
  },
  {
    name        => 'elevation_muscongus_bay_me_mllw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Muscongus Bay, ME (MLLW)'
  },
  {
    name        => 'elevation_mussel_creek_or_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Mussel Creek, OR (MHW)'
  },
  {
    name        => 'elevation_myers_creek_or_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Myers Creek, OR (MHW)'
  },
  {
    name        => 'elevation_myrtle_beach_sc_mhw',
    units       => 'm',
    description => 'NOAA Tsunami Inundation topography/bathymetry elevation for Myrtle Beach, SC (MHW)'
  },
  {
    name        => 'elevation_nantucket_harbor_ma_mhw',
    units       => 'm',
    description => 'NOAA Tsunami Inundation topography/bathymetry elevation for Nantucket, MA (MHW)'
  },
  {
    name        => 'elevation_nantucket_ma_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Nantucket Harbor, MA (MHW)'
  },
  {
    name        => 'elevation_naples_bay_fl_navd88',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Naples Bay, FL (NAVD88)'
  },
  {
    name        => 'elevation_narragansett_bay_1_ma_navd88',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Narragansett Bay 1, MA (NAVD88)'
  },
  {
    name        => 'elevation_narragansett_bay_1_ri_navd88',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Narragansett Bay 1, RI (NAVD88)'
  },
  {
    name        => 'elevation_narragansett_bay_2_ma_navd88',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Narragansett Bay 2, MA (NAVD88)'
  },
  {
    name        => 'elevation_narragansett_bay_2_ri_navd88',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Narragansett Bay 2, RI (NAVD88)'
  },
  {
    name        => 'elevation_narragansett_bay_3_ma_navd88',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Narragansett Bay 3, MA (NAVD88)'
  },
  {
    name        => 'elevation_narragansett_bay_3_ri_navd88',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Narragansett Bay 3, RI (NAVD88)'
  },
  {
    name        => 'elevation_narraguagas_bay_me_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Narraguagas Bay, ME (MHW)'
  },
  {
    name        => 'elevation_nassau_sound_fl_navd88',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Nassau Sound, FL (NAVD88)'
  },
  {
    name        => 'elevation_navarro_river_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Navarro River, CA (MHW)'
  },
  {
    name        => 'elevation_navesink_nj_navd88',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Navesink, NJ (NAVD88)'
  },
  {
    name        => 'elevation_nawiliwili_hi_mlw',
    units       => 'm',
    description => 'NOAA Tsunami Inundation topography/bathymetry elevation for Nawiliwili, HI (MLW)'
  },
  {
    name        => 'elevation_nc9_crescent_beach_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for NC9 Crescent Beach, CA (MHW)'
  },
  {
    name        => 'elevation_nc_intercoastal_waterway_nc_navd88',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for NC Intercoastal Waterway, NC (NAVD88)'
  },
  {
    name        => 'elevation_necanicum_river_estuary_or_mllw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Necanicum River Estuary, OR (MLLW)'
  },
  {
    name        => 'elevation_nehalem_bay_or_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Nehalem Bay, OR (MHW)'
  },
  {
    name        => 'elevation_neskowin_creek_or_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Neskowin Creek, OR (MHW)'
  },
  {
    name        => 'elevation_nestucca_bay_or_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Nestucca Bay, OR (MHW)'
  },
  {
    name        => 'elevation_netarts_bay_or_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Netarts Bay, OR (MHW)'
  },
  {
    name        => 'elevation_new_orleans_la_mhw',
    units       => 'm',
    description => 'NOAA Tsunami Inundation topography/bathymetry elevation for New Orleans, LA (MHW)'
  },
  {
    name        => 'elevation_new_orleans_la_mllw',
    units       => 'm',
    description => 'NOAA Tsunami Inundation topography/bathymetry elevation for New Orleans, LA (MLLW)'
  },
  {
    name        => 'elevation_new_orleans_la_navd88',
    units       => 'm',
    description => 'NOAA Tsunami Inundation topography/bathymetry elevation for New Orleans, LA (NAVD88)'
  },
  {
    name        => 'elevation_new_river_2_or_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for New River 2, OR (MHW)'
  },
  {
    name        => 'elevation_new_river_nc_navd88',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for New River, NC (NAVD88)'
  },
  {
    name        => 'elevation_new_river_sc_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for New River/Calibogue Sound, SC (MHW)'
  },
  {
    name        => 'elevation_new_york_harbor_1_nj_navd88',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for New York Harbor 1, NJ (NAVD88)'
  },
  {
    name        => 'elevation_new_york_harbor_1_ny_navd88',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for New York Harbor 1, NY (NAVD88)'
  },
  {
    name        => 'elevation_new_york_harbor_2_nj_navd88',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for New York Harbor 2, NJ (NAVD88)'
  },
  {
    name        => 'elevation_new_york_harbor_2_ny_navd88',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for New York Harbor 2, NY (NAVD88)'
  },
  {
    name        => 'elevation_newark_bay_kill_van_kull_nj_navd88',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Newark Bay Kill Van Kull, NJ (NAVD88)'
  },
  {
    name        => 'elevation_newark_bay_kill_van_kull_ny_navd88',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Newark Bay Kill Van Kull, NY (NAVD88)'
  },
  {
    name        => 'elevation_newport_harbor_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Newport Harbor, CA (MHW)'
  },
  {
    name        => 'elevation_newport_or_navd88',
    units       => 'm',
    description => 'NOAA Tsunami Inundation topography/bathymetry elevation for Newport, OR (NAVD88)'
  },
  {
    name        => 'elevation_nikolski_ak_mhw',
    units       => 'm',
    description => 'NOAA Tsunami Inundation topography/bathymetry elevation for Nikolski, AK (MHW)'
  },
  {
    name        => 'elevation_nj_coastal_bays_great_egg_harbor_nj_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for NJ Coastal Bays Great Egg Harbor, NJ (MHW)'
  },
  {
    name        => 'elevation_north_edisto_river_sc_navd88',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for North Edisto River, SC (NAVD88)'
  },
  {
    name        => 'elevation_north_south_rivers_ma_mllw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for North South Rivers, MA (MLLW)'
  },
  {
    name        => 'elevation_north_south_santee_rivers_sc_navd88',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for North South Santee Rivers, SC (NAVD88)'
  },
  {
    name        => 'elevation_northeast_atlantic_msl',
    units       => 'm',
    description => 'NOAA Tsunami Inundation topography/bathymetry elevation for Northeast Atlantic (MSL)'
  },
  {
    name        => 'elevation_northeast_river_md_navd88',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Northeast River, MD (NAVD88)'
  },
  {
    name        => 'elevation_northern_ca_msl',
    units       => 'm',
    description => 'NOAA Tsunami Inundation topography/bathymetry elevation for Northern CA (MSL)'
  },
  {
    name        => 'elevation_northern_gulf_navd88',
    units       => 'm',
    description => 'NOAA Tsunami Inundation topography/bathymetry elevation for Northern Gulf (NAVD88)'
  },
  {
    name        => 'elevation_northern_or_mhw',
    units       => 'm',
    description => 'NOAA Tsunami Inundation topography/bathymetry elevation for Northern OR (MHW)'
  },
  {
    name        => 'elevation_northwest_pacific_msl',
    units       => 'm',
    description => 'NOAA Tsunami Inundation topography/bathymetry elevation for Northwest Pacific (MSL)'
  },
  {
    name        => 'elevation_noyo_river_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Noyo River, CA (MHW)'
  },
  {
    name        => 'elevation_oahu_hi_mhw',
    units       => 'm',
    description => 'NOAA Tsunami Inundation topography/bathymetry elevation for Oahu, HI (MHW)'
  },
  {
    name        => 'elevation_ocean_city_md_mhw',
    units       => 'm',
    description => 'NOAA Tsunami Inundation topography/bathymetry elevation for Ocean City, MD (MHW)'
  },
  {
    name        => 'elevation_oceanside_harbor_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Oceanside Harbor, CA (MHW)'
  },
  {
    name        => 'elevation_old_creek_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Old Creek, CA (MHW)'
  },
  {
    name        => 'elevation_ormond_beach_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Ormond Beach, CA (MHW)'
  },
  {
    name        => 'elevation_ossabaw_sound_ga_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Ossabaw Sound, GA (MHW)'
  },
  {
    name        => 'elevation_ossagon_creek_2_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Ossagon Creek 2, CA (MHW)'
  },
  {
    name        => 'elevation_otter_cove_me_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Otter Cove, ME (MHW)'
  },
  {
    name        => 'elevation_outer_doctors_bay_fl_navd88',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Outer Doctors Bay, FL (NAVD88)'
  },
  {
    name        => 'elevation_oyster_creek_fl_mllw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Oyster Creek, FL (MLLW)'
  },
  {
    name        => 'elevation_oyster_creek_tx_navd88',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Oyster Creek, TX (NAVD88)'
  },
  {
    name        => 'elevation_ozette_river_wa_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Ozette River, WA (MHW)'
  },
  {
    name        => 'elevation_padilla_bay_wa_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Padilla Bay, WA (MHW)'
  },
  {
    name        => 'elevation_palm_beach_fl_mhw',
    units       => 'm',
    description => 'NOAA Tsunami Inundation topography/bathymetry elevation for Palm Beach, FL (MHW)'
  },
  {
    name        => 'elevation_palm_beach_fl_navd88',
    units       => 'm',
    description => 'NOAA Tsunami Inundation topography/bathymetry elevation for Palm Beach, FL (NAVD88)'
  },
  {
    name        => 'elevation_pamlico_sound_nc_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Pamlico Sound, NC (MHW)'
  },
  {
    name        => 'elevation_panama_city_fl_mhw',
    units       => 'm',
    description => 'NOAA Tsunami Inundation topography/bathymetry elevation for Panama City, FL (MHW)'
  },
  {
    name        => 'elevation_panama_city_fl_navd88',
    units       => 'm',
    description => 'NOAA Tsunami Inundation topography/bathymetry elevation for Panama City, FL (NAVD88)'
  },
  {
    name        => 'elevation_passamaquoddy_bay_st_croix_river_me_mllw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Passamaquoddy Bay St Croix River, ME (MLLWW)'
  },
  {
    name        => 'elevation_patapsco_river_md_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Patapsco River, MD (MHW)'
  },
  {
    name        => 'elevation_patuxent_river_md_navd88',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Patuxent River, MD (NAVD88)'
  },
  {
    name        => 'elevation_pearl_harbor_hi_mhw',
    units       => 'm',
    description => 'NOAA Tsunami Inundation topography/bathymetry elevation for Pearl Harbor, HI (MHW)'
  },
  {
    name        => 'elevation_pebble_beach_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Pebble Beach, CA (MHW)'
  },
  {
    name        => 'elevation_penobscot_bay_me_mllw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Penobscot Bay, ME (MLLW)'
  },
  {
    name        => 'elevation_perdido_bay_al_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Perdido Bay, AL (MHW)'
  },
  {
    name        => 'elevation_perdido_bay_fl_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Perdido Bay, FL (MHW)'
  },
  {
    name        => 'elevation_pescadero_creek_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Pescadero Creek, CA (MHW)'
  },
  {
    name        => 'elevation_piankatank_river_va_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Piankatank River, VA (MHW)'
  },
  {
    name        => 'elevation_pico_creek_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Pico Creek, CA (MHW)'
  },
  {
    name        => 'elevation_pilarcitos_creek_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Pilcitos Creek, CA (MHW)'
  },
  {
    name        => 'elevation_pillar_point_harbor_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Pillar Point Harbor, CA (MHW)'
  },
  {
    name        => 'elevation_pismo_creek_lagoon_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Pismo Creek Lagoon, CA (MHW)'
  },
  {
    name        => 'elevation_pistol_river_or_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Pistol River, OR (MHW)'
  },
  {
    name        => 'elevation_pithlachascotee_river_fl_mllw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Pithlachascotee River, FL (MLLW)'
  },
  {
    name        => 'elevation_pleasant_bay_ma_mllw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Pleasant Bay, ME (MHW)'
  },
  {
    name        => 'elevation_pleasant_bay_me_mhw',
    units       => 'm',
    description => 'NOAA Tsunami Inundation topography/bathymetry elevation for Pleasant Bay/Chatham Harbor, MA (MLLW)'
  },
  {
    name        => 'elevation_plum_island_sound_ma_mllw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Plum Island Sound, MA (MLLW)'
  },
  {
    name        => 'elevation_pomponio_creek_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Pomponio Creek, CA (MHW)'
  },
  {
    name        => 'elevation_ponce_de_leon_inlet_fl_navd88',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Ponce De Leon Inlet, FL (NAVD88)'
  },
  {
    name        => 'elevation_poquoson_back_rivers_va_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Poquoson Back Rivers, VA (MHW)'
  },
  {
    name        => 'elevation_port_alexander_ak_mhw',
    units       => 'm',
    description => 'NOAA Tsunami Inundation topography/bathymetry elevation for Port Alexander, AK (MHW)'
  },
  {
    name        => 'elevation_port_hueneme_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Port Hueneme, CA (MHW)'
  },
  {
    name        => 'elevation_port_orford_head_or_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Port Orford Head, OR (MHW)'
  },
  {
    name        => 'elevation_port_orford_or_mhw',
    units       => 'm',
    description => 'NOAA Tsunami Inundation topography/bathymetry elevation for Port Orford, OR (MHW)'
  },
  {
    name        => 'elevation_port_san_luis_ca_mhw',
    units       => 'm',
    description => 'NOAA Tsunami Inundation topography/bathymetry elevation for Port San Luis, CA (MHW)'
  },
  {
    name        => 'elevation_port_townsend_wa_mhw',
    units       => 'm',
    description => 'NOAA Tsunami Inundation topography/bathymetry elevation for Port Townsend, WA (MHW)'
  },
  {
    name        => 'elevation_port_townsend_wa_navd88',
    units       => 'm',
    description => 'NOAA Tsunami Inundation topography/bathymetry elevation for Port Townsend, WA (NAVD88)'
  },
  {
    name        => 'elevation_portland_me_mhw',
    units       => 'm',
    description => 'NOAA Tsunami Inundation topography/bathymetry elevation for Portland, ME (MHW)'
  },
  {
    name        => 'elevation_potomac_river_md_mllw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Potomac River, MD (MLLW)'
  },
  {
    name        => 'elevation_potomac_river_va_mllw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Potomac River, VA (MLLW)'
  },
  {
    name        => 'elevation_prince_william_sound_ak_mhhw',
    units       => 'm',
    description => 'NOAA Tsunami Inundation topography/bathymetry elevation for Prince William Sound, AK (MHHW)'
  },
  {
    name        => 'elevation_prospect_harbor_me_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Prospect Harbor, ME (MHW)'
  },
  {
    name        => 'elevation_pudding_creek_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Pudding Creek, CA (MHW)'
  },
  {
    name        => 'elevation_puget_sound_wa_mhw',
    units       => 'm',
    description => 'NOAA Tsunami Inundation topography/bathymetry elevation for Puget Sound, WA (MHW)'
  },
  {
    name        => 'elevation_puget_sound_wa_mllw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Puget Sound, WA (MLLW)'
  },
  {
    name        => 'elevation_puget_sound_wa_navd88',
    units       => 'm',
    description => 'NOAA Tsunami Inundation topography/bathymetry elevation for Puget Sound, WA (NAVD88)'
  },
  {
    name        => 'elevation_queets_river_wa_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Queets River, WA (MHW)'
  },
  {
    name        => 'elevation_quillayute_river_wa_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Quillayute River, WA (MHW)'
  },
  {
    name        => 'elevation_quinault_river_wa_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Quinault River, WA (MHW)'
  },
  {
    name        => 'elevation_quoddy_narrows_me_mllw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Quoddy Narrows, ME (MLLW)'
  },
  {
    name        => 'elevation_raft_river_wa_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Raft River, WA (MHW)'
  },
  {
    name        => 'elevation_rappahannock_river_va_navd88',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Rappahannock River, VA (NAVD88)'
  },
  {
    name        => 'elevation_raritan_bay_nj_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Raritan Bay, NJ (MHW)'
  },
  {
    name        => 'elevation_raritan_bay_ny_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Raritan Bay, NY (MHW)'
  },
  {
    name        => 'elevation_redwood_creek_big_lagoon_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Redwood Creek Big Lagoon, CA (MHW)'
  },
  {
    name        => 'elevation_redwood_creek_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Redwood Creek, CA (MHW)'
  },
  {
    name        => 'elevation_rehoboth_bay_de_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Rehoboth Bay, DE (MHW)'
  },
  {
    name        => 'elevation_rincon_creek_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Rincon Creek, CA (MHW)'
  },
  {
    name        => 'elevation_rock_1_bay_fl_mllw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Rock 1 Bay, FL (MLLW)'
  },
  {
    name        => 'elevation_rock_creek_or_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Rock Creek, OR (MHW)'
  },
  {
    name        => 'elevation_rodeo_lagoon_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Rodeo Lagoon, CA (MHW)'
  },
  {
    name        => 'elevation_rogue_river_or_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Rogue River, OR (MHW)'
  },
  {
    name        => 'elevation_russian_gulch_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Russian Gulch, CA (MHW)'
  },
  {
    name        => 'elevation_russian_gulch_creek_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Russian Gulch Creek, CA (MHW)'
  },
  {
    name        => 'elevation_russian_river_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Russian River, CA (MHW)'
  },
  {
    name        => 'elevation_sabine_lake_la_mllw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Sabine Lake, LA (MLLW)'
  },
  {
    name        => 'elevation_sabine_lake_tx_mllw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Sabine Lake, TX (MLLW)'
  },
  {
    name        => 'elevation_saco_bay_me_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Saco Bay, ME (MHW)'
  },
  {
    name        => 'elevation_saint_andrew_bay_fl_navd88',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Saint Andrew Bay, FL (NAVD88)'
  },
  {
    name        => 'elevation_saint_augustine_fl_navd88',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Saint Augustine, FL (NAVD88)'
  },
  {
    name        => 'elevation_saint_johns_river_fl_navd88',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Saint Johns River, FL (NAVD88)'
  },
  {
    name        => 'elevation_saint_joseph_bay_fl_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Saint Joseph Bay, FL (MHW)'
  },
  {
    name        => 'elevation_salem_sound_ma_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Salem Sound, MA (MHW)'
  },
  {
    name        => 'elevation_salmon_creek_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Salmon Creek, CA (MHW)'
  },
  {
    name        => 'elevation_salmon_river_estuary_or_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Salmon River Estuary, OR (MHW)'
  },
  {
    name        => 'elevation_salt_creek_2_fl_mllw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Salt Creek 2, FL (MLLW)'
  },
  {
    name        => 'elevation_salt_creek_fl_mllw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Salt Creek, FL (MLLW)'
  },
  {
    name        => 'elevation_samish_bay_wa_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Samish Bay, WA (MHW)'
  },
  {
    name        => 'elevation_san_antonio_bay_tx_mllw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for San Antonio Bay, TX (MLLW)'
  },
  {
    name        => 'elevation_san_antonio_creek_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for San Antonio Creek, CA (MHW)'
  },
  {
    name        => 'elevation_san_carpoforo_creek_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for San Carpoforo Creek, CA (MHW)'
  },
  {
    name        => 'elevation_san_diego_bay_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for San Diego Bay, CA (MHW)'
  },
  {
    name        => 'elevation_san_diego_ca_mhw',
    units       => 'm',
    description => 'NOAA Tsunami Inundation topography/bathymetry elevation for San Diego, CA (MHW)'
  },
  {
    name        => 'elevation_san_diego_ca_navd88',
    units       => 'm',
    description => 'NOAA Tsunami Inundation topography/bathymetry elevation for San Diego, CA (NAVD88)'
  },
  {
    name        => 'elevation_san_dieguito_river_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for San Dieguito River, CA (MHW)'
  },
  {
    name        => 'elevation_san_elijo_lagoon_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for San Elijo Lagoon, CA (MHW)'
  },
  {
    name        => 'elevation_san_francisco_bay_1_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for San Francisco Bay 1, CA (MHW)'
  },
  {
    name        => 'elevation_san_francisco_bay_2_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for San Francisco Bay 2, CA (MHW)'
  },
  {
    name        => 'elevation_san_francisco_bay_3_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for San Francisco Bay 3, CA (MHW)'
  },
  {
    name        => 'elevation_san_francisco_bay_4_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for San Francisco Bay 4, CA (MHW)'
  },
  {
    name        => 'elevation_san_francisco_bay_5_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for San Francisco Bay 5, CA (MHW)'
  },
  {
    name        => 'elevation_san_francisco_ca_mhw',
    units       => 'm',
    description => 'NOAA Tsunami Inundation topography/bathymetry elevation for San Francisco, CA (MHW)'
  },
  {
    name        => 'elevation_san_francisco_ca_navd88',
    units       => 'm',
    description => 'NOAA Tsunami Inundation topography/bathymetry elevation for San Francisco, CA (NAVD88)'
  },
  {
    name        => 'elevation_san_gregorio_creek_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for San Gregorio Creek, CA (MHW)'
  },
  {
    name        => 'elevation_san_juan_creek_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for San Juan Creek, CA (MHW)'
  },
  {
    name        => 'elevation_san_luiz_obispo_bay_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for San Luiz Obispo Bay, CA (MHW)'
  },
  {
    name        => 'elevation_san_luiz_rey_river_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for San Luiz Rey River, CA (MHW)'
  },
  {
    name        => 'elevation_san_mateo_lagoon_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for San Mateo Lagoon, CA (MHW)'
  },
  {
    name        => 'elevation_san_onofre_creek_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for San Onofre Creek, CA (MHW)'
  },
  {
    name        => 'elevation_san_pedro_bay_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for San Pedro Bay, CA (MHW)'
  },
  {
    name        => 'elevation_san_pedro_creek_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for San Pedro Creek, CA (MHW)'
  },
  {
    name        => 'elevation_san_simeon_creek_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for San Simeon Creek, CA (MHW)'
  },
  {
    name        => 'elevation_sand_lake_estuary_or_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Sand Lake Estuary, OR (MHW)'
  },
  {
    name        => 'elevation_sand_point_ak_mhw',
    units       => 'm',
    description => 'NOAA Tsunami Inundation topography/bathymetry elevation for Sand Point, AK (MHW)'
  },
  {
    name        => 'elevation_sanders_creek_fl_mllw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Sanders Creek, FL (MLLW)'
  },
  {
    name        => 'elevation_santa_ana_river_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Santa Ana River, CA (MHW)'
  },
  {
    name        => 'elevation_santa_barbara_ca_mhw',
    units       => 'm',
    description => 'NOAA Tsunami Inundation topography/bathymetry elevation for Santa Barbara, CA (MHW)'
  },
  {
    name        => 'elevation_santa_clara_river_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Santa Clara River, CA (MHW)'
  },
  {
    name        => 'elevation_santa_margarita_river_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Santa Margarita River, CA (MHW)'
  },
  {
    name        => 'elevation_santa_maria_river_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Santa Maria River, CA (MHW)'
  },
  {
    name        => 'elevation_santa_monica_bay_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Santa Monica Bay, CA (MHW)'
  },
  {
    name        => 'elevation_santa_monica_ca_mhw',
    units       => 'm',
    description => 'NOAA Tsunami Inundation topography/bathymetry elevation for Santa Monica, CA (MHW)'
  },
  {
    name        => 'elevation_santa_monica_ca_navd88',
    units       => 'm',
    description => 'NOAA Tsunami Inundation topography/bathymetry elevation for Santa Monica, CA (NAVD88)'
  },
  {
    name        => 'elevation_santa_monica_canyon_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Santa Monica Canyon, CA (MHW)'
  },
  {
    name        => 'elevation_santa_rosa_creek_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Santa Rosa Creek, CA (MHW)'
  },
  {
    name        => 'elevation_santa_ynez_river_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Santa Ynez River, CA (MHW)'
  },
  {
    name        => 'elevation_sapelo_sound_ga_navd88',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Sapelo Sound, GA (NAVD88)'
  },
  {
    name        => 'elevation_sarasota_bay_fl_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Sarasota Bay, FL (MHW)'
  },
  {
    name        => 'elevation_sassafrass_river_md_navd88',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Sassafrass River, MD (NAVD88)'
  },
  {
    name        => 'elevation_savanna_river_ga_mhw',
    units       => 'm',
    description => 'NOAA Tsunami Inundation topography/bathymetry elevation for Savannah, GA (MHW)'
  },
  {
    name        => 'elevation_savannah_ga_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Savanna River, GA (MHW)'
  },
  {
    name        => 'elevation_schwan_lagoon_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Schwan Lagoon, CA (MHW)'
  },
  {
    name        => 'elevation_scott_creek_lagoon_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Scott Creek Lagoon, CA (MHW)'
  },
  {
    name        => 'elevation_seaside_creek_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Seaside Creek, CA (MHW)'
  },
  {
    name        => 'elevation_seaside_or_mhw',
    units       => 'm',
    description => 'NOAA Tsunami Inundation topography/bathymetry elevation for Seaside, OR (MHW)'
  },
  {
    name        => 'elevation_seldovia_ak_mhw',
    units       => 'm',
    description => 'NOAA Tsunami Inundation topography/bathymetry elevation for Seldovia, AK (MHW)'
  },
  {
    name        => 'elevation_sengekontacket_pond_ma_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Sengekontacket Pond, MA (MHW)'
  },
  {
    name        => 'elevation_severn_river_md_navd88',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Severn River, MD (NAVD88)'
  },
  {
    name        => 'elevation_seward_ak_mhw',
    units       => 'm',
    description => 'NOAA Tsunami Inundation topography/bathymetry elevation for Seward, AK (MHW)'
  },
  {
    name        => 'elevation_sewee_bay_sc_navd88',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Sewee Bay, SC (NAVD88)'
  },
  {
    name        => 'elevation_shallotte_river_nc_navd88',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Shallotte River, NC (NAVD88)'
  },
  {
    name        => 'elevation_shark_river_nj_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Shark River, NJ (MHW)'
  },
  {
    name        => 'elevation_sheepscot_bay_me_mllw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Sheepscot Bay, ME (MLLW)'
  },
  {
    name        => 'elevation_shemya_ak_mhw',
    units       => 'm',
    description => 'NOAA Tsunami Inundation topography/bathymetry elevation for Shemya, AK (MHW)'
  },
  {
    name        => 'elevation_shinnecock_bay_ny_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Shinnecock Bay, NY (MHW)'
  },
  {
    name        => 'elevation_shrewsbury_river_nj_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Shrewsbury River, NJ (MHW)'
  },
  {
    name        => 'elevation_siletz_bay_or_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Siletz Bay, OR (MHW)'
  },
  {
    name        => 'elevation_siltcoos_river_or_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Siltcoos River, OR (MHW)'
  },
  {
    name        => 'elevation_sink_creek_fl_mllw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Sink Creek, FL (MLLW)'
  },
  {
    name        => 'elevation_sitka_ak_mhw',
    units       => 'm',
    description => 'NOAA Tsunami Inundation topography/bathymetry elevation for Sitka, AK (MHW)'
  },
  {
    name        => 'elevation_sixes_river_or_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Sixes River, OR (MHW)'
  },
  {
    name        => 'elevation_smith_lake_or_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Smith Lake, OR (MHW)'
  },
  {
    name        => 'elevation_smith_river_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Smith River, CA (MHW)'
  },
  {
    name        => 'elevation_somes_sound_sw_bay_me_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Somes Sound Sw Bay, ME (MHW)'
  },
  {
    name        => 'elevation_south_edisto_river_sc_navd88',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for South Edisto River, SC (NAVD88)'
  },
  {
    name        => 'elevation_south_padre_island_tx_mhw',
    units       => 'm',
    description => 'NOAA Tsunami Inundation topography/bathymetry elevation for South Padre Island, TX (MHW)'
  },
  {
    name        => 'elevation_south_padre_island_tx_navd88',
    units       => 'm',
    description => 'NOAA Tsunami Inundation topography/bathymetry elevation for South Padre Island, TX (NAVD88)'
  },
  {
    name        => 'elevation_southcentral_ak_mhw',
    units       => 'm',
    description => 'NOAA Tsunami Inundation topography/bathymetry elevation for Southcentral, AK (MHW)'
  },
  {
    name        => 'elevation_southeast_ak_mhhw',
    units       => 'm',
    description => 'NOAA Tsunami Inundation topography/bathymetry elevation for Southeast, AK (MHHW)'
  },
  {
    name        => 'elevation_southeast_atlantic_msl',
    units       => 'm',
    description => 'NOAA Tsunami Inundation topography/bathymetry elevation for Southeast Atlantic (MSL)'
  },
  {
    name        => 'elevation_southern_ca_msl',
    units       => 'm',
    description => 'NOAA Tsunami Inundation topography/bathymetry elevation for Southern CA (MSL)'
  },
  {
    name        => 'elevation_southern_cap_coastal_ponds_ma_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Southern Cap Coastal Ponds, MA (MHW)'
  },
  {
    name        => 'elevation_spring_warrior_creek_fl_mllw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Spring Warrior Creek, FL (MLLW)'
  },
  {
    name        => 'elevation_squashan_creek_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Squashan Creek, CA (MHW)'
  },
  {
    name        => 'elevation_squaw_creek_or_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Squaw Creek, OR (MHW)'
  },
  {
    name        => 'elevation_steamboat_creek_wa_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Steamboat Creek, WA (MHW)'
  },
  {
    name        => 'elevation_stone_lagoon_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Stone Lagoon, CA (MHW)'
  },
  {
    name        => 'elevation_storter_bay_fl_mllw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Storter Bay, FL (MLLW)'
  },
  {
    name        => 'elevation_strait_of_juan_de_fuca_wa_mllw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Strait Of Juan De Fuca, WA (MLLW)'
  },
  {
    name        => 'elevation_stump_beach_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Stump Beach, CA (MHW)'
  },
  {
    name        => 'elevation_suislaw_river_or_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Suislaw River, OR (MHW)'
  },
  {
    name        => 'elevation_sunset_bay_or_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Sunset Bay, OR (MHW)'
  },
  {
    name        => 'elevation_sutton_creek_or_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Sutton Creek, OR (MHW)'
  },
  {
    name        => 'elevation_suwanee_river_fl_mllw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Suwanee River, FL (MLLW)'
  },
  {
    name        => 'elevation_sycamore_creek_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Sycamore Creek, CA (MHW)'
  },
  {
    name        => 'elevation_tahkenitch_creek_or_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Tahkenitch Creek, OR (MHW)'
  },
  {
    name        => 'elevation_tahola_wa_mhw',
    units       => 'm',
    description => 'NOAA Tsunami Inundation topography/bathymetry elevation for Tahola, WA (MHW)'
  },
  {
    name        => 'elevation_tajiguas_creek_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Tajiguas Creek, CA (MHW)'
  },
  {
    name        => 'elevation_tampa_bay_1_fl_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Tampa Bay 1, FL (MHW)'
  },
  {
    name        => 'elevation_tampa_bay_2_fl_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Tampa Bay 2, FL (MHW)'
  },
  {
    name        => 'elevation_tampa_bay_3_fl_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Tampa Bay 3, FL (MHW)'
  },
  {
    name        => 'elevation_tampa_bay_fl_mhw',
    units       => 'm',
    description => 'NOAA Tsunami Inundation topography/bathymetry elevation for Tampa Bay, FL (MHW)'
  },
  {
    name        => 'elevation_tampa_bay_fl_navd88',
    units       => 'm',
    description => 'NOAA Tsunami Inundation topography/bathymetry elevation for Tampa Bay, FL (NAVD88)'
  },
  {
    name        => 'elevation_tatitlek_ak_mhhw',
    units       => 'm',
    description => 'NOAA Tsunami Inundation topography/bathymetry elevation for Tatitlek, AK (MHHW)'
  },
  {
    name        => 'elevation_ten_mile_river_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Ten Mile River, CA (MHW)'
  },
  {
    name        => 'elevation_ten_thousand_islands_marco_bay_fl_navd88',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Ten Thousand Islands Marco Bay, FL (NAVD88)'
  },
  {
    name        => 'elevation_tenmile_creek_1_or_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Tenmile Creek 1, OR (MHW)'
  },
  {
    name        => 'elevation_tenmile_creek_2_or_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Tenmile Creek 2, OR (MHW)'
  },
  {
    name        => 'elevation_tennessee_valley_lagoon_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Tennessee Valley Lagoon, CA (MHW)'
  },
  {
    name        => 'elevation_thomas_creek_or_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Thomas Creek, OR (MHW)'
  },
  {
    name        => 'elevation_threemile_creek_or_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Threemile Creek, OR (MHW)'
  },
  {
    name        => 'elevation_tijuana_river_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Tijuana River, CA (MHW)'
  },
  {
    name        => 'elevation_tillamook_bay_or_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Tillamook Bay, OR (MHW)'
  },
  {
    name        => 'elevation_tomales_bay_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Tomales Bay, CA (MHW)'
  },
  {
    name        => 'elevation_trancas_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Trancas, CA (MHW)'
  },
  {
    name        => 'elevation_tunitas_creek_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Tunitas Creek, CA (MHW)'
  },
  {
    name        => 'elevation_twomile_creek_2_or_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Twomile Creek 2, OR (MHW)'
  },
  {
    name        => 'elevation_umpqua_river_or_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Umpqua River, OR (MHW)'
  },
  {
    name        => 'elevation_unalaska_ak_mhhw',
    units       => 'm',
    description => 'NOAA Tsunami Inundation topography/bathymetry elevation for Unalaska, AK (MHHW)'
  },
  {
    name        => 'elevation_upper_laguna_madre_tx_mllw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Upper Laguna Madre, TX (MLLW)'
  },
  {
    name        => 'elevation_usal_creek_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Usal Creek, CA (MHW)'
  },
  {
    name        => 'elevation_ventura_harbor_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Ventura Harbor, CA (MHW)'
  },
  {
    name        => 'elevation_ventura_river_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Ventura River, CA (MHW)'
  },
  {
    name        => 'elevation_villa_creek_lagoon_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Villa Creek Lagoon, CA (MHW)'
  },
  {
    name        => 'elevation_virgin_creek_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Virgin Creek, CA (MHW)'
  },
  {
    name        => 'elevation_virginia_beach_va_mhw',
    units       => 'm',
    description => 'NOAA Tsunami Inundation topography/bathymetry elevation for Virginia Beach, VA (MHW)'
  },
  {
    name        => 'elevation_waccasassa_river_fl_mllw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Waccasassa River, FL (MLLW)'
  },
  {
    name        => 'elevation_waddell_creek_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Waddell Creek, CA (MHW)'
  },
  {
    name        => 'elevation_wages_creek_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Wages Creek, CA (MHW)'
  },
  {
    name        => 'elevation_weeki_wachee_bay_fl_mllw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Weeki Wachee Bay, FL (MLLW)'
  },
  {
    name        => 'elevation_western_gulf_msl',
    units       => 'm',
    description => 'NOAA Tsunami Inundation topography/bathymetry elevation for Western Gulf (MSL)'
  },
  {
    name        => 'elevation_whackup_creek_fl_mllw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Whackup Creek, FL (MLLW)'
  },
  {
    name        => 'elevation_whale_creek_wa_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Whale Creek, WA (MHW)'
  },
  {
    name        => 'elevation_whaleshead_creek_or_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Whaleshead Creek, OR (MHW)'
  },
  {
    name        => 'elevation_whisky_run_or_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Whisky Run, OR (MHW)'
  },
  {
    name        => 'elevation_whitewater_bay_fl_mllw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Whitewater Bay, FL (MLLW)'
  },
  {
    name        => 'elevation_whittier_ak_mhhw',
    units       => 'm',
    description => 'NOAA Tsunami Inundation topography/bathymetry elevation for Whittier, AK (MHHW)'
  },
  {
    name        => 'elevation_wilder_creek_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Wilder Creek, CA (MHW)'
  },
  {
    name        => 'elevation_willapa_bay_wa_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Willapa Bay, WA (MHW)'
  },
  {
    name        => 'elevation_wilson_creek_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Wilson Creek, CA (MHW)'
  },
  {
    name        => 'elevation_winchuck_river_or_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Winchuck River, OR (MHW)'
  },
  {
    name        => 'elevation_withlacoochee_bay_fl_mllw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Withlacoochee Bay, FL (MLLW)'
  },
  {
    name        => 'elevation_wreck_creek_wa_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Wreck Creek, WA (MHW)'
  },
  {
    name        => 'elevation_yachats_river_or_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Yachats River, OR (MHW)'
  },
  {
    name        => 'elevation_yakutat_ak_mhhw',
    units       => 'm',
    description => 'NOAA Tsunami Inundation topography/bathymetry elevation for Yakutat, AK (MHHW)'
  },
  {
    name        => 'elevation_yaquina_bay_or_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Yaquina Bay, OR (MHW)'
  },
  {
    name        => 'elevation_yates_creek_fl_mllw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Yates Creek, FL (MLLW)'
  },
  {
    name        => 'elevation_younger_lagoon_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Younger Lagoon, CA (MHW)'
  },
  {
    name        => 'elevation_zuma_canyon_ca_mhw',
    units       => 'm',
    description => 'NOAA Estuarine Bathymetry elevation for Zuma Canyon, CA (MHW)'
  },
  {
    name        => 'annual_nadp_wet_nh4',
    units       => 'kgN/ha/yr',
    description => 'NADP annual wet deposition of ammonium'
  },
  {
    name        => 'annual_nadp_wet_no3',
    units       => 'kgN/ha/yr',
    description => 'NADP annual wet deposition of nitrate'
  },
  {
    name        => 'annual_nadp_wet_totn',
    units       => 'kgN/ha/yr',
    description => 'NADP annual wet deposition of total nitrogen'
  },
  {
    name        => 'annual_nadpcmaq_dry_nh3',
    units       => 'kgN/ha/yr',
    description => 'NADP+CMAQ annual dry deposition of ammonia'
  },
  {
    name        => 'annual_nadpcmaq_dry_nh4',
    units       => 'kgN/ha/yr',
    description => 'NADP+CMAQ annual dry deposition of ammonium'
  },
  {
    name        => 'annual_nadpcmaq_dry_no3',
    units       => 'kgN/ha/yr',
    description => 'NADP+CMAQ annual dry deposition of nitrate'
  },
  {
    name        => 'annual_nadpcmaq_dry_totn',
    units       => 'kgN/ha/yr',
    description => 'NADP+CMAQ annual dry deposition of total nitrogen'
  },
  {
    name        => 'annual_nadpcmaq_dry_totnp',
    units       => 'kgN/ha/yr',
    description => 'NADP+CMAQ annual dry deposition of total nitrogen as a percent of dry+wet'
  },
  {
    name        => 'annual_nadpcmaq_dry_totun',
    units       => 'kgN/ha/yr',
    description => 'NADP+CMAQ annual dry deposition of total unmeasured nitrogen species'
  },
  {
    name        => 'annual_nadpcmaq_dry_totunp',
    units       => 'kgN/ha/yr',
    description => 'NADP+CMAQ annual dry deposition of total unmeasured nitrogen species as a percent of dry+wet'
  },
  {
    name        => 'annual_nadpcmaq_dry_net_nh3',
    units       => 'kgN/ha/yr',
    description => 'NADP+CMAQ annual dry deposition of net flux (deposition - emission) ammonia'
  },
  {
    name        => 'annual_nadpcmaq_dry_oxidized_totn',
    units       => 'kgN/ha/yr',
    description => 'NADP+CMAQ annual dry deposition of oxidized total nitrogen'
  },
  {
    name        => 'annual_nadpcmaq_dry_oxidized_totnp',
    units       => 'kgN/ha/yr',
    description => 'NADP+CMAQ annual dry deposition of oxidized total nitrogen as a percent of dry+wet'
  },
  {
    name        => 'annual_nadpcmaq_dry_reduced_totn',
    units       => 'kgN/ha/yr',
    description => 'NADP+CMAQ annual dry deposition of reduced total nitrogen'
  },
  {
    name        => 'annual_nadpcmaq_dry_reduced_totnp',
    units       => 'kgN/ha/yr',
    description => 'NADP+CMAQ annual dry deposition of reduced total nitrogen as a percent of dry+wet'
  },
  {
    name        => 'annual_nadpcmaq_wet_nh4',
    units       => 'kgN/ha/yr',
    description => 'NADP+CMAQ annual wet deposition of ammonium'
  },
  {
    name        => 'annual_nadpcmaq_wet_no3',
    units       => 'kgN/ha/yr',
    description => 'NADP+CMAQ annual wet deposition of nitrate'
  },
  {
    name        => 'annual_nadpcmaq_wet_totn',
    units       => 'kgN/ha/yr',
    description => 'NADP+CMAQ annual wet deposition of total nitrogen'
  },
  {
    name        => 'annual_nadpcmaq_oxidized_totn',
    units       => 'kgN/ha/yr',
    description => 'NADP+CMAQ annual deposition of oxidized total nitrogen'
  },
  {
    name        => 'annual_nadpcmaq_oxidized_totnp',
    units       => 'kgN/ha/yr',
    description => 'NADP+CMAQ annual deposition of oxidized total nitrogen as a percent of dry+wet'
  },
  {
    name        => 'annual_nadpcmaq_reduced_totn',
    units       => 'kgN/ha/yr',
    description => 'NADP+CMAQ annual deposition of reduced total nitrogen'
  },
  {
    name        => 'annual_nadpcmaq_reduced_totnp',
    units       => 'kgN/ha/yr',
    description => 'NADP+CMAQ annual deposition of reduced total nitrogen as a percent of dry+wet'
  },
  {
    name        => 'annual_nadpcmaq_total_totn',
    units       => 'kgN/ha/yr',
    description => 'NADP+CMAQ annual deposition of total nitrogen'
  },
  {
    name        => 'monthly_temperature_minimum',
    units       => 'degC',
    description => 'PRISM 30-year (1971-2000) average monthly low temperature '
  },
  {
    name        => 'monthly_temperature_maximum',
    units       => 'degC',
    description => 'PRISM 30-year (1971-2000) average monthly high temperature'
  },
  {
    name        => 'monthly_precipitation_mean',
    units       => 'mm',
    description => 'PRISM 30-year (1971-2000) average monthly precipitation'
  },
  {
    name        => 'annual_temperature_minimum',
    units       => 'degC',
    description => 'PRISM 30-year (1971-2000) average annual low temperature '
  },
  {
    name        => 'annual_temperature_maximum',
    units       => 'degC',
    description => 'PRISM 30-year (1971-2000) average annual high temperature'
  },
  {
    name        => 'annual_precipitation_mean',
    units       => 'mm',
    description => 'PRISM 30-year (1971-2000) average annual precipitation'
  },
  {
    name        => 'seagrass',
    units       => '-',
    description => 'USGS seagrass density'
  },
  {
    name        => 'nitrogen.source_gridded.crops',
    units       => 'kgN/ha/yr',
    description => 'EPA NADP nitrogen gridded source crops'
  },
  {
    name        => 'nitrogen.source_gridded.fertilizer',
    units       => 'kgN/ha/yr',
    description => 'EPA NADP nitrogen gridded source fertilizer'
  },
  {
    name        => 'nitrogen.source_gridded.manure',
    units       => 'kgN/ha/yr',
    description => 'EPA NADP nitrogen gridded source manure'
  },
  {
    name        => 'impervious_ball_mountain_vt',
    units       => '-',
    description => 'EPA one meter impervious grid cell Ball Mountain, VT watershed'
  },
  {
    name        => 'impervious_black_river_vt',
    units       => '-',
    description => 'EPA one meter impervious grid cell Black River, VT watershed'
  },
  {
    name        => 'impervious_bloody_brook_vt',
    units       => '-',
    description => 'EPA one meter impervious grid cell Bloody Brook , VT watershed'
  },
  {
    name        => 'impervious_branch_brook_vt',
    units       => '-',
    description => 'EPA one meter impervious grid cell Branch Brook, VT watershed'
  },
  {
    name        => 'impervious_calendar_brook_vt',
    units       => '-',
    description => 'EPA one meter impervious grid cell Calendar Brook, VT watershed'
  },
  {
    name        => 'impervious_east_creek_vt',
    units       => '-',
    description => 'EPA one meter impervious grid cell East Creek, VT watershed'
  },
  {
    name        => 'impervious_jay_branch_vt',
    units       => '-',
    description => 'EPA one meter impervious grid cell Jay Branch, VT watershed'
  },
  {
    name        => 'impervious_jewett_brook_vt',
    units       => '-',
    description => 'EPA one meter impervious grid cell Jewett Brook, VT watershed'
  },
  {
    name        => 'impervious_lulls_brook_vt',
    units       => '-',
    description => 'EPA one meter impervious grid cell Lulls Brook, VT watershed'
  },
  {
    name        => 'impervious_mallets_bay_vt',
    units       => '-',
    description => 'EPA one meter impervious grid cell Mallets Bay, VT watershed'
  },
  {
    name        => 'impervious_mallets_bay10_vt',
    units       => '-',
    description => 'EPA one meter impervious grid cell Mallets Bay (HUC 10), VT watershed'
  },
  {
    name        => 'impervious_mallets_creek_vt',
    units       => '-',
    description => 'EPA one meter impervious grid cell Mallets Creek, VT watershed'
  },
  {
    name        => 'impervious_mill_river_vt',
    units       => '-',
    description => 'EPA one meter impervious grid cell Mill Brook, VT watershed'
  },
  {
    name        => 'impervious_mill_brook_vt',
    units       => '-',
    description => 'EPA one meter impervious grid cell Mill River, VT watershed'
  },
  {
    name        => 'impervious_moon_brook_vt',
    units       => '-',
    description => 'EPA one meter impervious grid cell Moon Brook, VT watershed'
  },
  {
    name        => 'impervious_muddy_brook_vt',
    units       => '-',
    description => 'EPA one meter impervious grid cell Muddy Brook, VT watershed'
  },
  {
    name        => 'impervious_north_branch_deerfield_river_vt',
    units       => '-',
    description => 'EPA one meter impervious grid cell North Branch Deerfield River, VT watershed'
  },
  {
    name        => 'impervious_ottauquechee_river_vt',
    units       => '-',
    description => 'EPA one meter impervious grid cell Ottauquechee River, VT watershed'
  },
  {
    name        => 'impervious_ottauquechee_river_headwaters_vt',
    units       => '-',
    description => 'EPA one meter impervious grid cell Ottauquechee River Headwaters, VT watershed'
  },
  {
    name        => 'impervious_shelburne_river_vt',
    units       => '-',
    description => 'EPA one meter impervious grid cell Shelburne River, VT watershed'
  },
  {
    name        => 'impervious_snipe_island_vt',
    units       => '-',
    description => 'EPA one meter impervious grid cell Snipe Island, VT watershed'
  },
  {
    name        => 'impervious_winooski_river_vt',
    units       => '-',
    description => 'EPA one meter impervious grid cell Winooski River, VT watershed'
  },
  {
    name        => 'meris_weekly_cyanobacteria_al',
    units       => 'cells/mL',
    description => 'MERIS satellite 300 meter weekly composite cyanobacteria Alabama 2002-05-21 to 2012-04-01'
  },
  {
    name        => 'meris_weekly_cyanobacteria_ar',
    units       => 'cells/mL',
    description => 'MERIS satellite 300 meter weekly composite cyanobacteria Arkansas 2002-05-21 to 2012-04-01'
  },
  {
    name        => 'meris_weekly_cyanobacteria_az',
    units       => 'cells/mL',
    description => 'MERIS satellite 300 meter weekly composite cyanobacteria Arizona 2002-05-21 to 2012-04-01'
  },
  {
    name        => 'meris_weekly_cyanobacteria_ca',
    units       => 'cells/mL',
    description => 'MERIS+MERIS satellite 300 meter weekly composite cyanobacteria California 2002-05-21 to 2012-04-01'
  },
  {
    name        => 'meris_weekly_cyanobacteria_co',
    units       => 'cells/mL',
    description => 'MERIS satellite 300 meter weekly composite cyanobacteria Colorado 2002-05-21 to 2012-04-01'
  },
  {
    name        => 'meris_weekly_cyanobacteria_ct',
    units       => 'cells/mL',
    description => 'MERIS satellite 300 meter weekly composite cyanobacteria Connecticut 2002-05-21 to 2012-04-01'
  },
  {
    name        => 'meris_weekly_cyanobacteria_de',
    units       => 'cells/mL',
    description => 'MERIS satellite 300 meter weekly composite cyanobacteria Delaware 2002-05-21 to 2012-04-01'
  },
  {
    name        => 'meris_weekly_cyanobacteria_fl',
    units       => 'cells/mL',
    description => 'MERIS+MERIS satellite 300 meter weekly composite cyanobacteria Florida 2002-05-21 to 2012-04-01'
  },
  {
    name        => 'meris_weekly_cyanobacteria_ga',
    units       => 'cells/mL',
    description => 'MERIS satellite 300 meter weekly composite cyanobacteria Georgia 2002-05-21 to 2012-04-01'
  },
  {
    name        => 'meris_weekly_cyanobacteria_ia',
    units       => 'cells/mL',
    description => 'MERIS satellite 300 meter weekly composite cyanobacteria Iowa 2002-05-21 to 2012-04-01'
  },
  {
    name        => 'meris_weekly_cyanobacteria_id',
    units       => 'cells/mL',
    description => 'MERIS satellite 300 meter weekly composite cyanobacteria Idaho 2002-05-21 to 2012-04-01'
  },
  {
    name        => 'meris_weekly_cyanobacteria_il',
    units       => 'cells/mL',
    description => 'MERIS satellite 300 meter weekly composite cyanobacteria Illinois 2002-05-21 to 2012-04-01'
  },
  {
    name        => 'meris_weekly_cyanobacteria_in',
    units       => 'cells/mL',
    description => 'MERIS satellite 300 meter weekly composite cyanobacteria Indiana 2002-05-21 to 2012-04-01'
  },
  {
    name        => 'meris_weekly_cyanobacteria_ks',
    units       => 'cells/mL',
    description => 'MERIS satellite 300 meter weekly composite cyanobacteria Kansas 2002-05-21 to 2012-04-01'
  },
  {
    name        => 'meris_weekly_cyanobacteria_ky',
    units       => 'cells/mL',
    description => 'MERIS satellite 300 meter weekly composite cyanobacteria Kentucky 2002-05-21 to 2012-04-01'
  },
  {
    name        => 'meris_weekly_cyanobacteria_la',
    units       => 'cells/mL',
    description => 'MERIS satellite 300 meter weekly composite cyanobacteria Louisiana 2002-05-21 to 2012-04-01'
  },
  {
    name        => 'meris_weekly_cyanobacteria_ma',
    units       => 'cells/mL',
    description => 'MERIS satellite 300 meter weekly composite cyanobacteria Massachusetts 2002-05-21 to 2012-04-01'
  },
  {
    name        => 'meris_weekly_cyanobacteria_md',
    units       => 'cells/mL',
    description => 'MERIS satellite 300 meter weekly composite cyanobacteria Maryland 2002-05-21 to 2012-04-01'
  },
  {
    name        => 'meris_weekly_cyanobacteria_me',
    units       => 'cells/mL',
    description => 'MERIS satellite 300 meter weekly composite cyanobacteria Maine 2002-05-21 to 2012-04-01'
  },
  {
    name        => 'meris_weekly_cyanobacteria_mi',
    units       => 'cells/mL',
    description => 'MERIS satellite 300 meter weekly composite cyanobacteria Michigan 2002-05-21 to 2012-04-01'
  },
  {
    name        => 'meris_weekly_cyanobacteria_mn',
    units       => 'cells/mL',
    description => 'MERIS satellite 300 meter weekly composite cyanobacteria Minnesota 2002-05-21 to 2012-04-01'
  },
  {
    name        => 'meris_weekly_cyanobacteria_mo',
    units       => 'cells/mL',
    description => 'MERIS satellite 300 meter weekly composite cyanobacteria Missouri 2002-05-21 to 2012-04-01'
  },
  {
    name        => 'meris_weekly_cyanobacteria_ms',
    units       => 'cells/mL',
    description => 'MERIS satellite 300 meter weekly composite cyanobacteria Mississippi 2002-05-21 to 2012-04-01'
  },
  {
    name        => 'meris_weekly_cyanobacteria_mt',
    units       => 'cells/mL',
    description => 'MERIS satellite 300 meter weekly composite cyanobacteria Montana 2002-05-21 to 2012-04-01'
  },
  {
    name        => 'meris_weekly_cyanobacteria_nc',
    units       => 'cells/mL',
    description => 'MERIS satellite 300 meter weekly composite cyanobacteria North Carolina 2002-05-21 to 2012-04-01'
  },
  {
    name        => 'meris_weekly_cyanobacteria_nd',
    units       => 'cells/mL',
    description => 'MERIS satellite 300 meter weekly composite cyanobacteria North Dakota 2002-05-21 to 2012-04-01'
  },
  {
    name        => 'meris_weekly_cyanobacteria_ne',
    units       => 'cells/mL',
    description => 'MERIS satellite 300 meter weekly composite cyanobacteria Nebraska 2002-05-21 to 2012-04-01'
  },
  {
    name        => 'meris_weekly_cyanobacteria_nh',
    units       => 'cells/mL',
    description => 'MERIS satellite 300 meter weekly composite cyanobacteria New Hampshire 2002-05-21 to 2012-04-01'
  },
  {
    name        => 'meris_weekly_cyanobacteria_nj',
    units       => 'cells/mL',
    description => 'MERIS satellite 300 meter weekly composite cyanobacteria New Jersey 2002-05-21 to 2012-04-01'
  },
  {
    name        => 'meris_weekly_cyanobacteria_nm',
    units       => 'cells/mL',
    description => 'MERIS satellite 300 meter weekly composite cyanobacteria New Mexico 2002-05-21 to 2012-04-01'
  },
  {
    name        => 'meris_weekly_cyanobacteria_nv',
    units       => 'cells/mL',
    description => 'MERIS satellite 300 meter weekly composite cyanobacteria Nevada 2002-05-21 to 2012-04-01'
  },
  {
    name        => 'meris_weekly_cyanobacteria_ny',
    units       => 'cells/mL',
    description => 'MERIS satellite 300 meter weekly composite cyanobacteria New York 2002-05-21 to 2012-04-01'
  },
  {
    name        => 'meris_weekly_cyanobacteria_oh',
    units       => 'cells/mL',
    description => 'MERIS satellite 300 meter weekly composite cyanobacteria Ohio 2002-05-21 to 2012-04-01'
  },
  {
    name        => 'meris_weekly_cyanobacteria_ok',
    units       => 'cells/mL',
    description => 'MERIS satellite 300 meter weekly composite cyanobacteria Oklahoma 2002-05-21 to 2012-04-01'
  },
  {
    name        => 'meris_weekly_cyanobacteria_or',
    units       => 'cells/mL',
    description => 'MERIS satellite 300 meter weekly composite cyanobacteria Oregon 2002-05-21 to 2012-04-01'
  },
  {
    name        => 'meris_weekly_cyanobacteria_pa',
    units       => 'cells/mL',
    description => 'MERIS satellite 300 meter weekly composite cyanobacteria Pennsylvania 2002-05-21 to 2012-04-01'
  },
  {
    name        => 'meris_weekly_cyanobacteria_ri',
    units       => 'cells/mL',
    description => 'MERIS satellite 300 meter weekly composite cyanobacteria Rhode Island 2002-05-21 to 2012-04-01'
  },
  {
    name        => 'meris_weekly_cyanobacteria_sc',
    units       => 'cells/mL',
    description => 'MERIS satellite 300 meter weekly composite cyanobacteria South Carolina 2002-05-21 to 2012-04-01'
  },
  {
    name        => 'meris_weekly_cyanobacteria_sd',
    units       => 'cells/mL',
    description => 'MERIS satellite 300 meter weekly composite cyanobacteria South Dakota 2002-05-21 to 2012-04-01'
  },
  {
    name        => 'meris_weekly_cyanobacteria_tn',
    units       => 'cells/mL',
    description => 'MERIS satellite 300 meter weekly composite cyanobacteria Tennessee 2002-05-21 to 2012-04-01'
  },
  {
    name        => 'meris_weekly_cyanobacteria_tx',
    units       => 'cells/mL',
    description => 'MERIS satellite 300 meter weekly composite cyanobacteria Texas 2002-05-21 to 2012-04-01'
  },
  {
    name        => 'meris_weekly_cyanobacteria_ut',
    units       => 'cells/mL',
    description => 'MERIS satellite 300 meter weekly composite cyanobacteria Utah 2002-05-21 to 2012-04-01'
  },
  {
    name        => 'meris_weekly_cyanobacteria_va',
    units       => 'cells/mL',
    description => 'MERIS satellite 300 meter weekly composite cyanobacteria Virginia 2002-05-21 to 2012-04-01'
  },
  {
    name        => 'meris_weekly_cyanobacteria_vt',
    units       => 'cells/mL',
    description => 'MERIS satellite 300 meter weekly composite cyanobacteria Vermont 2002-05-21 to 2012-04-01'
  },
  {
    name        => 'meris_weekly_cyanobacteria_wa',
    units       => 'cells/mL',
    description => 'MERIS satellite 300 meter weekly composite cyanobacteria Washington 2002-05-21 to 2012-04-01'
  },
  {
    name        => 'meris_weekly_cyanobacteria_wi',
    units       => 'cells/mL',
    description => 'MERIS satellite 300 meter weekly composite cyanobacteria Wisconsin 2002-05-21 to 2012-04-01'
  },
  {
    name        => 'meris_weekly_cyanobacteria_wv',
    units       => 'cells/mL',
    description => 'MERIS satellite 300 meter weekly composite cyanobacteria West Virginia 2002-05-21 to 2012-04-01'
  },
  {
    name        => 'meris_weekly_cyanobacteria_wy',
    units       => 'cells/mL',
    description => 'MERIS satellite 300 meter weekly composite cyanobacteria Wyoming 2002-05-21 to 2012-04-01'
  },
  {
    name        => 'sentinel_weekly_cyanobacteria_al',
    units       => 'cells/mL',
    description => 'SENTINEL satellite 300 meter weekly composite cyanobacteria Alabama 2017-01-01 to 2018-08-12'
  },
  {
    name        => 'sentinel_weekly_cyanobacteria_ar',
    units       => 'cells/mL',
    description => 'SENTINEL satellite 300 meter weekly composite cyanobacteria Arkansas 2017-01-01 to 2018-08-12'
  },
  {
    name        => 'sentinel_weekly_cyanobacteria_az',
    units       => 'cells/mL',
    description => 'SENTINEL satellite 300 meter weekly composite cyanobacteria Arizona 2017-01-01 to 2018-08-12'
  },
  {
    name        => 'sentinel_weekly_cyanobacteria_ca',
    units       => 'cells/mL',
    description => 'MERIS+SENTINEL satellite 300 meter weekly composite cyanobacteria California 2017-01-01 to 2018-06-03'
  },
  {
    name        => 'sentinel_weekly_cyanobacteria_co',
    units       => 'cells/mL',
    description => 'SENTINEL satellite 300 meter weekly composite cyanobacteria Colorado 2017-01-01 to 2018-08-12'
  },
  {
    name        => 'sentinel_weekly_cyanobacteria_ct',
    units       => 'cells/mL',
    description => 'SENTINEL satellite 300 meter weekly composite cyanobacteria Connecticut 2017-01-01 to 2018-08-12'
  },
  {
    name        => 'sentinel_weekly_cyanobacteria_de',
    units       => 'cells/mL',
    description => 'SENTINEL satellite 300 meter weekly composite cyanobacteria Delaware 2017-01-01 to 2018-08-12'
  },
  {
    name        => 'sentinel_weekly_cyanobacteria_fl',
    units       => 'cells/mL',
    description => 'MERIS+SENTINEL satellite 300 meter weekly composite cyanobacteria Florida 2002-05-21 to 2018-08-12'
  },
  {
    name        => 'sentinel_weekly_cyanobacteria_ga',
    units       => 'cells/mL',
    description => 'SENTINEL satellite 300 meter weekly composite cyanobacteria Georgia 2017-01-01 to 2018-08-12'
  },
  {
    name        => 'sentinel_weekly_cyanobacteria_ia',
    units       => 'cells/mL',
    description => 'SENTINEL satellite 300 meter weekly composite cyanobacteria Iowa 2017-01-01 to 2018-08-12'
  },
  {
    name        => 'sentinel_weekly_cyanobacteria_id',
    units       => 'cells/mL',
    description => 'SENTINEL satellite 300 meter weekly composite cyanobacteria Idaho 2017-01-01 to 2018-08-12'
  },
  {
    name        => 'sentinel_weekly_cyanobacteria_il',
    units       => 'cells/mL',
    description => 'SENTINEL satellite 300 meter weekly composite cyanobacteria Illinois 2017-01-01 to 2018-08-12'
  },
  {
    name        => 'sentinel_weekly_cyanobacteria_in',
    units       => 'cells/mL',
    description => 'SENTINEL satellite 300 meter weekly composite cyanobacteria Indiana 2017-01-01 to 2018-08-12'
  },
  {
    name        => 'sentinel_weekly_cyanobacteria_ks',
    units       => 'cells/mL',
    description => 'SENTINEL satellite 300 meter weekly composite cyanobacteria Kansas 2017-01-01 to 2018-08-12'
  },
  {
    name        => 'sentinel_weekly_cyanobacteria_ky',
    units       => 'cells/mL',
    description => 'SENTINEL satellite 300 meter weekly composite cyanobacteria Kentucky 2017-01-01 to 2018-08-12'
  },
  {
    name        => 'sentinel_weekly_cyanobacteria_la',
    units       => 'cells/mL',
    description => 'SENTINEL satellite 300 meter weekly composite cyanobacteria Louisiana 2017-01-01 to 2018-08-12'
  },
  {
    name        => 'sentinel_weekly_cyanobacteria_ma',
    units       => 'cells/mL',
    description => 'SENTINEL satellite 300 meter weekly composite cyanobacteria Massachusetts 2017-01-01 to 2018-08-12'
  },
  {
    name        => 'sentinel_weekly_cyanobacteria_md',
    units       => 'cells/mL',
    description => 'SENTINEL satellite 300 meter weekly composite cyanobacteria Maryland 2017-01-01 to 2018-08-12'
  },
  {
    name        => 'sentinel_weekly_cyanobacteria_me',
    units       => 'cells/mL',
    description => 'SENTINEL satellite 300 meter weekly composite cyanobacteria Maine 2017-01-01 to 2018-08-12'
  },
  {
    name        => 'sentinel_weekly_cyanobacteria_mi',
    units       => 'cells/mL',
    description => 'SENTINEL satellite 300 meter weekly composite cyanobacteria Michigan 2017-01-01 to 2018-08-12'
  },
  {
    name        => 'sentinel_weekly_cyanobacteria_mn',
    units       => 'cells/mL',
    description => 'SENTINEL satellite 300 meter weekly composite cyanobacteria Minnesota 2017-01-01 to 2018-08-12'
  },
  {
    name        => 'sentinel_weekly_cyanobacteria_mo',
    units       => 'cells/mL',
    description => 'SENTINEL satellite 300 meter weekly composite cyanobacteria Missouri 2017-01-01 to 2018-08-12'
  },
  {
    name        => 'sentinel_weekly_cyanobacteria_ms',
    units       => 'cells/mL',
    description => 'SENTINEL satellite 300 meter weekly composite cyanobacteria Mississippi 2017-01-01 to 2018-08-12'
  },
  {
    name        => 'sentinel_weekly_cyanobacteria_mt',
    units       => 'cells/mL',
    description => 'SENTINEL satellite 300 meter weekly composite cyanobacteria Montana 2017-01-01 to 2018-08-12'
  },
  {
    name        => 'sentinel_weekly_cyanobacteria_nc',
    units       => 'cells/mL',
    description => 'SENTINEL satellite 300 meter weekly composite cyanobacteria North Carolina 2017-01-01 to 2017-12-31'
  },
  {
    name        => 'sentinel_weekly_cyanobacteria_nd',
    units       => 'cells/mL',
    description => 'SENTINEL satellite 300 meter weekly composite cyanobacteria North Dakota 2017-01-01 to 2018-08-12'
  },
  {
    name        => 'sentinel_weekly_cyanobacteria_ne',
    units       => 'cells/mL',
    description => 'SENTINEL satellite 300 meter weekly composite cyanobacteria Nebraska 2017-01-01 to 2018-08-12'
  },
  {
    name        => 'sentinel_weekly_cyanobacteria_nh',
    units       => 'cells/mL',
    description => 'SENTINEL satellite 300 meter weekly composite cyanobacteria New Hampshire 2017-01-01 to 2018-08-12'
  },
  {
    name        => 'sentinel_weekly_cyanobacteria_nj',
    units       => 'cells/mL',
    description => 'SENTINEL satellite 300 meter weekly composite cyanobacteria New Jersey 2017-01-01 to 2018-08-12'
  },
  {
    name        => 'sentinel_weekly_cyanobacteria_nm',
    units       => 'cells/mL',
    description => 'SENTINEL satellite 300 meter weekly composite cyanobacteria New Mexico 2017-01-01 to 2018-08-12'
  },
  {
    name        => 'sentinel_weekly_cyanobacteria_nv',
    units       => 'cells/mL',
    description => 'SENTINEL satellite 300 meter weekly composite cyanobacteria Nevada 2017-01-01 to 2018-08-12'
  },
  {
    name        => 'sentinel_weekly_cyanobacteria_ny',
    units       => 'cells/mL',
    description => 'SENTINEL satellite 300 meter weekly composite cyanobacteria New York 2017-01-01 to 2018-06-03'
  },
  {
    name        => 'sentinel_weekly_cyanobacteria_oh',
    units       => 'cells/mL',
    description => 'SENTINEL satellite 300 meter weekly composite cyanobacteria Ohio 2017-01-01 to 2018-08-12'
  },
  {
    name        => 'sentinel_weekly_cyanobacteria_ok',
    units       => 'cells/mL',
    description => 'SENTINEL satellite 300 meter weekly composite cyanobacteria Oklahoma 2017-01-01 to 2018-08-12'
  },
  {
    name        => 'sentinel_weekly_cyanobacteria_or',
    units       => 'cells/mL',
    description => 'SENTINEL satellite 300 meter weekly composite cyanobacteria Oregon 2017-01-01 to 2018-06-03'
  },
  {
    name        => 'sentinel_weekly_cyanobacteria_pa',
    units       => 'cells/mL',
    description => 'SENTINEL satellite 300 meter weekly composite cyanobacteria Pennsylvania 2017-01-01 to 2018-08-12'
  },
  {
    name        => 'sentinel_weekly_cyanobacteria_ri',
    units       => 'cells/mL',
    description => 'SENTINEL satellite 300 meter weekly composite cyanobacteria Rhode Island 2017-01-01 to 2018-08-12'
  },
  {
    name        => 'sentinel_weekly_cyanobacteria_sc',
    units       => 'cells/mL',
    description => 'SENTINEL satellite 300 meter weekly composite cyanobacteria South Carolina 2017-01-01 to 2018-08-12'
  },
  {
    name        => 'sentinel_weekly_cyanobacteria_sd',
    units       => 'cells/mL',
    description => 'SENTINEL satellite 300 meter weekly composite cyanobacteria South Dakota 2017-01-01 to 2018-08-12'
  },
  {
    name        => 'sentinel_weekly_cyanobacteria_tn',
    units       => 'cells/mL',
    description => 'SENTINEL satellite 300 meter weekly composite cyanobacteria Tennessee 2017-01-01 to 2018-08-12'
  },
  {
    name        => 'sentinel_weekly_cyanobacteria_tx',
    units       => 'cells/mL',
    description => 'SENTINEL satellite 300 meter weekly composite cyanobacteria Texas 2017-01-01 to 2018-08-12'
  },
  {
    name        => 'sentinel_weekly_cyanobacteria_ut',
    units       => 'cells/mL',
    description => 'SENTINEL satellite 300 meter weekly composite cyanobacteria Utah 2017-01-01 to 2018-08-12'
  },
  {
    name        => 'sentinel_weekly_cyanobacteria_va',
    units       => 'cells/mL',
    description => 'SENTINEL satellite 300 meter weekly composite cyanobacteria Virginia 2017-01-01 to 2018-08-12'
  },
  {
    name        => 'sentinel_weekly_cyanobacteria_vt',
    units       => 'cells/mL',
    description => 'SENTINEL satellite 300 meter weekly composite cyanobacteria Vermont 2017-01-01 to 2018-08-12'
  },
  {
    name        => 'sentinel_weekly_cyanobacteria_wa',
    units       => 'cells/mL',
    description => 'SENTINEL satellite 300 meter weekly composite cyanobacteria Washington 2017-01-01 to 2018-08-12'
  },
  {
    name        => 'sentinel_weekly_cyanobacteria_wi',
    units       => 'cells/mL',
    description => 'SENTINEL satellite 300 meter weekly composite cyanobacteria Wisconsin 2017-01-01 to 2018-08-12'
  },
  {
    name        => 'sentinel_weekly_cyanobacteria_wv',
    units       => 'cells/mL',
    description => 'SENTINEL satellite 300 meter weekly composite cyanobacteria West Virginia 2017-01-01 to 2018-08-12'
  },
  {
    name        => 'sentinel_weekly_cyanobacteria_wy',
    units       => 'cells/mL',
    description => 'SENTINEL satellite 300 meter weekly composite cyanobacteria Wyoming 2017-01-01 to 2018-08-12'
  },
  {
    name        => 'states',
    units       => '-',
    description => 'US states'
  },
  {
    name        => 'provinces',
    units       => '-',
    description => 'Canadian Provinces'
  },
  {
    name        => 'counties',
    units       => '-',
    description => 'US counties'
  },
  {
    name        => 'cities',
    units       => '-',
    description => 'US cities'
  },
  {
    name        => 'roads',
    units       => '-',
    description => 'US roads'
  },
  {
    name        => 'tributaries',
    units       => '-',
    description => 'NHDFlowline stream level 1 major tributaries'
  },
  {
    name        => 'hucs',
    units       => '-',
    description => 'USGS Hydrologic Units Catalog'
  },
  {
    name        => 'watersheds',
    units       => '-',
    description => 'NHD/EPA watersheds'
  },
  {
    name        => 'estuaries',
    units       => '-',
    description => 'USGS/EPA estuaries'
  },
  {
    name        => 'lakes',
    units       => '-',
    description => 'MERIS/EPA lakes'
  },
  {
    name        => 'coastal_zone_management_areas',
    units       => '-',
    description => 'NOAA Coastal Zone Management Areas'
  },
  {
    name        => 'legislative_districts',
    units       => '-',
    description => 'US Census Bureau Congressional Legislative Districts'
  },
  {
    name        => 'wetlands',
    units       => '-',
    description => 'National wetlands classification'
  },
  {
    name        => 'population_iclus1',
    units       => '-',
    description => 'Integrated Climate and Land Use System version 1 projected population'
  },
  {
    name        => 'population_iclus2',
    units       => '-',
    description => 'Integrated Climate and Land Use System version 2 projected population'
  },
  {
    name        => 'land_use_iclus_watershed',
    units       => '%',
    description => 'Integrated Climate and Land Use System version 2 projected decadal land use 2000-2100 percent of each category per watershed'
  },
  {
    name        => 'land_use_nlcd_2011_watershed',
    units       => '%',
    description => 'USGS National Land Cover Dataset 2011 percent of each category per watershed'
  },
  {
    name        => 'land_use_nlcd_2006_watershed',
    units       => '%',
    description => 'USGS National Land Cover Dataset 2006 percent of each category per watershed'
  },
  {
    name        => 'soil',
    units       => '-',
    description => 'USGS STATSGO soil characterization'
  },
  {
    name        => 'nitrogen.estuary_cmaq',
    units       => '-',
    description => 'USGS/EPA CMAQ modeled nitrogen in estuaries'
  },
  {
    name        => 'dry_nh3.estuary_cmaq',
    units       => '-',
    description => 'USGS/EPA CMAQ modeled dry deposition of NH3 in estuaries'
  },
  {
    name        => 'dry_natural_nh3.estuary_cmaq',
    units       => '-',
    description => 'USGS/EPA CMAQ modeled dry deposition of natural NH3 in estuaries'
  },
  {
    name        => 'dry_nhx.estuary_cmaq',
    units       => '-',
    description => 'USGS/EPA CMAQ modeled dry deposition of NHx in estuaries'
  },
  {
    name        => 'dry_so2.estuary_cmaq',
    units       => '-',
    description => 'USGS/EPA CMAQ modeled dry deposition of SO2 in estuaries'
  },
  {
    name        => 'dry_oxidized_nox.estuary_cmaq',
    units       => '-',
    description => 'USGS/EPA CMAQ modeled dry deposition of oxidized NOx in estuaries'
  },
  {
    name        => 'dry_oxidized_organic_nitrogen.estuary_cmaq',
    units       => '-',
    description => 'USGS/EPA CMAQ modeled dry deposition of oxidized organic nitrogen in estuaries'
  },
  {
    name        => 'dry_oxidized_no3.estuary_cmaq',
    units       => '-',
    description => 'USGS/EPA CMAQ modeled dry deposition of oxidized NO3 in estuaries'
  },
  {
    name        => 'dry_oxidized_nitrogen.estuary_cmaq',
    units       => '-',
    description => 'USGS/EPA CMAQ modeled dry deposition of oxidized nitrogen in estuaries'
  },
  {
    name        => 'dry_reduced_nitrogen.estuary_cmaq',
    units       => '-',
    description => 'USGS/EPA CMAQ modeled dry deposition of reduced nitrogen in estuaries'
  },
  {
    name        => 'dry_sulfur.estuary_cmaq',
    units       => '-',
    description => 'USGS/EPA CMAQ modeled dry deposition of sulfur in estuaries'
  },
  {
    name        => 'chloride.estuary_cmaq',
    units       => '-',
    description => 'USGS/EPA CMAQ modeled deposition of chloride in estuaries'
  },
  {
    name        => 'total_nitrogen.estuary_cmaq',
    units       => '-',
    description => 'USGS/EPA CMAQ modeled deposition of total nitrogen in estuaries'
  },
  {
    name        => 'oxidized_nitrogen.estuary_cmaq',
    units       => '-',
    description => 'USGS/EPA CMAQ modeled deposition of oxidized nitrogen in estuaries'
  },
  {
    name        => 'reduced_nitrogen.estuary_cmaq',
    units       => '-',
    description => 'USGS/EPA CMAQ modeled deposition of reduced nitrogen in estuaries'
  },
  {
    name        => 'sulfur.estuary_cmaq',
    units       => '-',
    description => 'USGS/EPA CMAQ modeled deposition of sulfur in estuaries'
  },
  {
    name        => 'wet_nh3.estuary_cmaq',
    units       => '-',
    description => 'USGS/EPA CMAQ modeled wet deposition of NH3 in estuaries'
  },
  {
    name        => 'wet_chloride.estuary_cmaq',
    units       => '-',
    description => 'USGS/EPA CMAQ modeled wet deposition of chloride in estuaries'
  },
  {
    name        => 'wet_so4.estuary_cmaq',
    units       => '-',
    description => 'USGS/EPA CMAQ modeled wet deposition of SO4 in estuaries'
  },
  {
    name        => 'wet_oxidized_no3.estuary_cmaq',
    units       => '-',
    description => 'USGS/EPA CMAQ modeled wet deposition of oxidized NO3 in estuaries'
  },
  {
    name        => 'wet_oxidized_nox.estuary_cmaq',
    units       => '-',
    description => 'USGS/EPA CMAQ modeled wet deposition of oxidized NOx in estuaries'
  },
  {
    name        => 'wet_organic_nitrogen.estuary_cmaq',
    units       => '-',
    description => 'USGS/EPA CMAQ modeled wet deposition of organic nitrogen in estuaries'
  },
  {
    name        => 'wet_oxidized_nitrogen.estuary_cmaq',
    units       => '-',
    description => 'USGS/EPA CMAQ modeled wet deposition of oxidized in estuaries'
  },
  {
    name        => 'wet_reduced_nitrogen.estuary_cmaq',
    units       => '-',
    description => 'USGS/EPA CMAQ modeled wet deposition of reduced nitrogen in estuaries'
  },
  {
    name        => 'wet_sulfur.estuary_cmaq',
    units       => '-',
    description => 'USGS/EPA CMAQ modeled wet deposition of sulfur in estuaries'
  },
  {
    name        => 'nitrogen.subestuary_cmaq',
    units       => '-',
    description => 'USGS/EPA CMAQ modeled nitrogen in subestuaries'
  },
  {
    name        => 'dry_nh3.subestuary_cmaq',
    units       => '-',
    description => 'USGS/EPA CMAQ modeled dry deposition of NH3 in subestuaries'
  },
  {
    name        => 'dry_natural_nh3.subestuary_cmaq',
    units       => '-',
    description => 'USGS/EPA CMAQ modeled dry deposition of natural NH3 in subestuaries'
  },
  {
    name        => 'dry_nhx.subestuary_cmaq',
    units       => '-',
    description => 'USGS/EPA CMAQ modeled dry deposition of NHx in subestuaries'
  },
  {
    name        => 'dry_so2.subestuary_cmaq',
    units       => '-',
    description => 'USGS/EPA CMAQ modeled dry deposition of SO2 in subestuaries'
  },
  {
    name        => 'dry_oxidized_nox.subestuary_cmaq',
    units       => '-',
    description => 'USGS/EPA CMAQ modeled dry deposition of oxidized NOx in subestuaries'
  },
  {
    name        => 'dry_oxidized_organic_nitrogen.subestuary_cmaq',
    units       => '-',
    description => 'USGS/EPA CMAQ modeled dry deposition of oxidized organic nitrogen in subestuaries'
  },
  {
    name        => 'dry_oxidized_no3.subestuary_cmaq',
    units       => '-',
    description => 'USGS/EPA CMAQ modeled dry deposition of oxidized NO3 in subestuaries'
  },
  {
    name        => 'dry_oxidized_nitrogen.subestuary_cmaq',
    units       => '-',
    description => 'USGS/EPA CMAQ modeled dry deposition of oxidized nitrogen in subestuaries'
  },
  {
    name        => 'dry_reduced_nitrogen.subestuary_cmaq',
    units       => '-',
    description => 'USGS/EPA CMAQ modeled dry deposition of reduced nitrogen in subestuaries'
  },
  {
    name        => 'dry_sulfur.subestuary_cmaq',
    units       => '-',
    description => 'USGS/EPA CMAQ modeled dry deposition of sulfur in subestuaries'
  },
  {
    name        => 'chloride.subestuary_cmaq',
    units       => '-',
    description => 'USGS/EPA CMAQ modeled deposition of chloride in subestuaries'
  },
  {
    name        => 'total_nitrogen.subestuary_cmaq',
    units       => '-',
    description => 'USGS/EPA CMAQ modeled deposition of total nitrogen in subestuaries'
  },
  {
    name        => 'oxidized_nitrogen.subestuary_cmaq',
    units       => '-',
    description => 'USGS/EPA CMAQ modeled deposition of oxidized nitrogen in subestuaries'
  },
  {
    name        => 'reduced_nitrogen.subestuary_cmaq',
    units       => '-',
    description => 'USGS/EPA CMAQ modeled deposition of reduced nitrogen in subestuaries'
  },
  {
    name        => 'sulfur.subestuary_cmaq',
    units       => '-',
    description => 'USGS/EPA CMAQ modeled deposition of sulfur in subestuaries'
  },
  {
    name        => 'wet_nh3.subestuary_cmaq',
    units       => '-',
    description => 'USGS/EPA CMAQ modeled wet deposition of NH3 in subestuaries'
  },
  {
    name        => 'wet_chloride.subestuary_cmaq',
    units       => '-',
    description => 'USGS/EPA CMAQ modeled wet deposition of chloride in subestuaries'
  },
  {
    name        => 'wet_so4.subestuary_cmaq',
    units       => '-',
    description => 'USGS/EPA CMAQ modeled wet deposition of SO4 in subestuaries'
  },
  {
    name        => 'wet_oxidized_no3.subestuary_cmaq',
    units       => '-',
    description => 'USGS/EPA CMAQ modeled wet deposition of oxidized NO3 in subestuaries'
  },
  {
    name        => 'wet_oxidized_nox.subestuary_cmaq',
    units       => '-',
    description => 'USGS/EPA CMAQ modeled wet deposition of oxidized NOx in subestuaries'
  },
  {
    name        => 'wet_organic_nitrogen.subestuary_cmaq',
    units       => '-',
    description => 'USGS/EPA CMAQ modeled wet deposition of organic nitrogen in subestuaries'
  },
  {
    name        => 'wet_oxidized_nitrogen.subestuary_cmaq',
    units       => '-',
    description => 'USGS/EPA CMAQ modeled wet deposition of oxidized in subestuaries'
  },
  {
    name        => 'wet_reduced_nitrogen.subestuary_cmaq',
    units       => '-',
    description => 'USGS/EPA CMAQ modeled wet deposition of reduced nitrogen in subestuaries'
  },
  {
    name        => 'wet_sulfur.subestuary_cmaq',
    units       => '-',
    description => 'USGS/EPA CMAQ modeled wet deposition of sulfur in subestuaries'
  },
  {
    name        => 'nitrogen.load_estuary_sparrow_1992',
    units       => '-',
    description => 'USGS SPARROW 1992 model estuary nitrogen loads'
  },
  {
    name        => 'nitrogen.load_estuary_sparrow_2002_mrb1',
    units       => '-',
    description => 'USGS SPARROW 2002 model estuary major river basin 1 nitrogen loads'
  },
  {
    name        => 'nitrogen.load_estuary_sparrow_2002_mrb2',
    units       => '-',
    description => 'USGS SPARROW 2002 model estuary major river basin 2 nitrogen loads'
  },
  {
    name        => 'nitrogen.load_estuary_sparrow_2002_mrb5',
    units       => '-',
    description => 'USGS SPARROW 2002 model estuary major river basin 5 nitrogen loads'
  },
  {
    name        => 'nitrogen.load_estuary_sparrow_2002_mrb7',
    units       => '-',
    description => 'USGS SPARROW 2002 model estuary major river basin 7 nitrogen loads'
  },
  {
    name        => 'nitrogen.load_estuary_sparrow_2011_mrb1mrb2',
    units       => '-',
    description => 'USGS SPARROW 2011 model estuary major river basin 1 and 2 nitrogen loads'
  },
  {
    name        => 'nitrogen.load_estuary_sparrow_2011_mrb5',
    units       => '-',
    description => 'USGS SPARROW 2011 model estuary major river basin 5 nitrogen loads'
  },
  {
    name        => 'nitrogen.load_estuary_sparrow_2011_mrb7',
    units       => '-',
    description => 'USGS SPARROW 2011 model estuary major river basin 7 nitrogen loads'
  },
  {
    name        => 'nitrogen.load_estuary_sparrow_2011_mrb8',
    units       => '-',
    description => 'USGS SPARROW 2011 model estuary major river basin 8 nitrogen loads'
  },
  {
    name        => 'nitrogen.load_subestuary_sparrow_1992',
    units       => '-',
    description => 'USGS SPARROW 1992 model subestuary nitrogen loads'
  },
  {
    name        => 'nitrogen.load_subestuary_sparrow_2002_mrb1',
    units       => '-',
    description => 'USGS SPARROW 2002 model subestuary major river basin 1 nitrogen loads'
  },
  {
    name        => 'nitrogen.load_subestuary_sparrow_2002_mrb2',
    units       => '-',
    description => 'USGS SPARROW 2002 model subestuary major river basin 2 nitrogen loads'
  },
  {
    name        => 'nitrogen.load_subestuary_sparrow_2002_mrb5',
    units       => '-',
    description => 'USGS SPARROW 2002 model subestuary major river basin 5 nitrogen loads'
  },
  {
    name        => 'nitrogen.load_subestuary_sparrow_2002_mrb7',
    units       => '-',
    description => 'USGS SPARROW 2002 model subestuary major river basin 7 nitrogen loads'
  },
  {
    name        => 'total_nitrogen.source_watershed_nadp',
    units       => '-',
    description => 'USGS/EPA NADP modeled deposition of nitrogen in watersheds'
  },
  {
    name        => 'wet_nitrogen.source_watershed_nadp',
    units       => '-',
    description => 'USGS/EPA NADP modeled wet deposition of nitrogen in watersheds'
  },
  {
    name        => 'dry_nitrogen.source_watershed_nadp',
    units       => '-',
    description => 'USGS/EPA NADP modeled dry deposition of nitrogen in watersheds'
  },
  {
    name        => 'oxidized_nitrogen.source_watershed_nadp',
    units       => '-',
    description => 'USGS/EPA NADP modeled deposition of oxidized nitrogen in watersheds'
  },
  {
    name        => 'reduced_nitrogen.source_watershed_nadp',
    units       => '-',
    description => 'USGS/EPA NADP modeled deposition of reduced nitrogen in watersheds'
  },
  {
    name        => 'dry_oxidized_nitrogen.source_watershed_nadp',
    units       => '-',
    description => 'USGS/EPA NADP modeled deposition of dry oxidized nitrogen in watersheds'
  },
  {
    name        => 'dry_reduced_nitrogen.source_watershed_nadp',
    units       => '-',
    description => 'USGS/EPA NADP modeled deposition of dry reduced nitrogen in watersheds'
  },
  {
    name        => 'wet_no3.source_watershed_nadp',
    units       => '-',
    description => 'USGS/EPA NADP modeled wet deposition of nitrate in watersheds'
  },
  {
    name        => 'wet_nh4.source_watershed_nadp',
    units       => '-',
    description => 'USGS/EPA NADP modeled wet deposition of ammonium in watersheds'
  },
  {
    name        => 'dry_no3.source_watershed_nadp',
    units       => '-',
    description => 'USGS/EPA NADP modeled dry deposition nitrate in watersheds'
  },
  {
    name        => 'dry_nh4.source_watershed_nadp',
    units       => '-',
    description => 'USGS/EPA NADP modeled dry deposition of ammonium in watersheds'
  },
  {
    name        => 'nitrogen.source_watershed_cmaq',
    units       => '-',
    description => 'USGS/EPA CMAQ modeled nitrogen in watersheds'
  },
  {
    name        => 'dry_nh3.source_watershed_cmaq',
    units       => '-',
    description => 'USGS/EPA CMAQ modeled dry deposition of NH3 in watersheds'
  },
  {
    name        => 'dry_natural_nh3.source_watershed_cmaq',
    units       => '-',
    description => 'USGS/EPA CMAQ modeled dry deposition of natural NH3 in watersheds'
  },
  {
    name        => 'dry_nhx.source_watershed_cmaq',
    units       => '-',
    description => 'USGS/EPA CMAQ modeled dry deposition of NHx in watersheds'
  },
  {
    name        => 'dry_so2.source_watershed_cmaq',
    units       => '-',
    description => 'USGS/EPA CMAQ modeled dry deposition of SO2 in watersheds'
  },
  {
    name        => 'dry_oxidized_nox.source_watershed_cmaq',
    units       => '-',
    description => 'USGS/EPA CMAQ modeled dry deposition of oxidized NOx in watersheds'
  },
  {
    name        => 'dry_oxidized_organic_nitrogen.source_watershed_cmaq',
    units       => '-',
    description => 'USGS/EPA CMAQ modeled dry deposition of oxidized organic nitrogen in watersheds'
  },
  {
    name        => 'dry_oxidized_no3.source_watershed_cmaq',
    units       => '-',
    description => 'USGS/EPA CMAQ modeled dry deposition of oxidized NO3 in watersheds'
  },
  {
    name        => 'dry_oxidized_nitrogen.source_watershed_cmaq',
    units       => '-',
    description => 'USGS/EPA CMAQ modeled dry deposition of oxidized nitrogen in watersheds'
  },
  {
    name        => 'dry_reduced_nitrogen.source_watershed_cmaq',
    units       => '-',
    description => 'USGS/EPA CMAQ modeled dry deposition of reduced nitrogen in watersheds'
  },
  {
    name        => 'dry_sulfur.source_watershed_cmaq',
    units       => '-',
    description => 'USGS/EPA CMAQ modeled dry deposition of sulfur in watersheds'
  },
  {
    name        => 'chloride.source_watershed_cmaq',
    units       => '-',
    description => 'USGS/EPA CMAQ modeled deposition of chloride in watersheds'
  },
  {
    name        => 'total_nitrogen.source_watershed_cmaq',
    units       => '-',
    description => 'USGS/EPA CMAQ modeled deposition of total nitrogen in watersheds'
  },
  {
    name        => 'oxidized_nitrogen.source_watershed_cmaq',
    units       => '-',
    description => 'USGS/EPA CMAQ modeled deposition of oxidized nitrogen in watersheds'
  },
  {
    name        => 'reduced_nitrogen.source_watershed_cmaq',
    units       => '-',
    description => 'USGS/EPA CMAQ modeled deposition of reduced nitrogen in watersheds'
  },
  {
    name        => 'sulfur.source_watershed_cmaq',
    units       => '-',
    description => 'USGS/EPA CMAQ modeled deposition of sulfur in watersheds'
  },
  {
    name        => 'wet_nh3.source_watershed_cmaq',
    units       => '-',
    description => 'USGS/EPA CMAQ modeled wet deposition of NH3 in watersheds'
  },
  {
    name        => 'wet_chloride.source_watershed_cmaq',
    units       => '-',
    description => 'USGS/EPA CMAQ modeled wet deposition of chloride in watersheds'
  },
  {
    name        => 'wet_so4.source_watershed_cmaq',
    units       => '-',
    description => 'USGS/EPA CMAQ modeled wet deposition of SO4 in watersheds'
  },
  {
    name        => 'wet_oxidized_no3.source_watershed_cmaq',
    units       => '-',
    description => 'USGS/EPA CMAQ modeled wet deposition of oxidized NO3 in watersheds'
  },
  {
    name        => 'wet_oxidized_nox.source_watershed_cmaq',
    units       => '-',
    description => 'USGS/EPA CMAQ modeled wet deposition of oxidized NOx in watersheds'
  },
  {
    name        => 'wet_organic_nitrogen.source_watershed_cmaq',
    units       => '-',
    description => 'USGS/EPA CMAQ modeled wet deposition of organic nitrogen in watersheds'
  },
  {
    name        => 'wet_oxidized_nitrogen.source_watershed_cmaq',
    units       => '-',
    description => 'USGS/EPA CMAQ modeled wet deposition of oxidized in watersheds'
  },
  {
    name        => 'wet_reduced_nitrogen.source_watershed_cmaq',
    units       => '-',
    description => 'USGS/EPA CMAQ modeled wet deposition of reduced nitrogen in watersheds'
  },
  {
    name        => 'wet_sulfur.source_watershed_cmaq',
    units       => '-',
    description => 'USGS/EPA CMAQ modeled wet deposition of sulfur in watersheds'
  },
  {
    name        => 'nitrogen.source_watershed_point',
    units       => '-',
    description => 'USGS/EPA Point-source nitrogen in watersheds'
  },
  {
    name        => 'nitrogen.source_watershed_nonpoint',
    units       => '-',
    description => 'USGS/EPA Non-point source nitrogen in watersheds'
  },
  {
    name        => 'phosphorus.source_watershed_point',
    units       => '-',
    description => 'USGS/EPA Point-source phosphorus in watersheds'
  },
  {
    name        => 'sediment_nca',
    units       => '-',
    description => 'EPA-NCA seabed sediment point multi-data 1990-2006'
  },
  {
    name        => 'sediment_nca_1990-2006',
    units       => '-',
    description => 'EPA-NCA seabed sediment point multi-data 1990-2006'
  },
  {
    name        => 'sediment_nca_2015',
    units       => '-',
    description => 'EPA-NCA seabed sediment point multi-data 2015'
  },
  {
    name        => 'sediment_calculated',
    units       => '-',
    description => 'USGS seabed sediment point multi-data (calculated)'
  },
  {
    name        => 'sediment_extracted',
    units       => '-',
    description => 'USGS seabed sediment point multi-data (extracted)'
  },
  {
    name        => 'sediment_parsed',
    units       => '-',
    description => 'USGS seabed sediment point multi-data (parsed)'
  },
  {
    name        => 'sediment_kriged',
    units       => '-',
    description => 'The Nature Conservancy kriged sediment'
  },
  {
    name        => 'coastal_vulnerability',
    units       => '-',
    description => 'Coastal vulnerability measures'
  },
  {
    name        => 'sensitivity_estuary',
    units       => '-',
    description => 'Estuary nutrient sensitivity measures'
  },
  {
    name        => 'groundwater_contact_catchment',
    units       => '-',
    description => 'Groundwater contact time (days) in catchments.'
  },
  {
    name        => 'groundwater_contact_watershed',
    units       => 'days',
    description => 'Groundwater contact time (days) in watersheds.'
  },
  {
    name        => 'stream_discharge_nhd_line',
    units       => 'days',
    description => 'NHD stream discharge measures'
  },
  {
    name        => 'stream_discharge_usgs_point',
    units       => '-',
    description => 'USGS gauge discharge measures'
  },
  {
    name        => 'tide_point',
    units       => '-',
    description => 'NOAA gauge tide height and volume measures'
  },
  {
    name        => 'tide_current',
    units       => '-',
    description => 'NOAA gauge tide current measures'
  },
  {
    name        => 'longshore_current',
    units       => '-',
    description => 'UCSD high-frequency gauge network longshore current measures'
  },
  {
    name        => 'estuary_flushing',
    units       => '-',
    description => 'NOAA gauge estuary flushing measures'
  },
  {
    name        => 'salinity_point',
    units       => '-',
    description => 'NOAA gauge estuary salinity measures'
  },
  {
    name        => 'flowlines_upper_ct_river',
    units       => '-',
    description => 'USGS SPARROW Modeled nutrients in flowlines of Upper Connecticut River'
  },
  {
    name        => 'flowlines_puget_sound_watershed',
    units       => '-',
    description => 'USGS SPARROW Modeled nutrients in _point of Puget Sound watershed'
  },
  {
    name        => 'discharge_points_federal_puget_sound_watershed',
    units       => '-',
    description => 'Federal-licensed discharge points in Puget Sound watershed'
  },
  {
    name        => 'discharge_points_state_puget_sound_watershed',
    units       => '-',
    description => 'State-licensed discharge points in Puget Sound watershed'
  },
  {
    name        => 'gi_bmp_installations_new_england',
    units       => '-',
    description => 'Green Infrastructure Best Management Practices Installations - New England'
  },
  {
    name        => 'impervious_nhdplus_new_england',
    units       => '-',
    description => 'Impervious Cover NHDPlus - New England'
  },
  {
    name        => 'land_change_2020_current_chesapeake',
    units       => '-',
    description => 'Chesapeake Bay Land Change Model 2020 prediction of current policy'
  },
  {
    name        => 'land_change_2030_current_chesapeake',
    units       => '-',
    description => 'Chesapeake Bay Land Change Model 2030 prediction of current policy'
  },
  {
    name        => 'land_change_2040_current_chesapeake',
    units       => '-',
    description => 'Chesapeake Bay Land Change Model 2040 prediction of current policy'
  },
  {
    name        => 'land_change_2020_deregulation_chesapeake',
    units       => '-',
    description => 'Chesapeake Bay Land Change Model 2020 prediction of deregulation policy'
  },
  {
    name        => 'land_change_2030_deregulation_chesapeake',
    units       => '-',
    description => 'Chesapeake Bay Land Change Model 2030 prediction of deregulation policy'
  },
  {
    name        => 'land_change_2040_deregulation_chesapeake',
    units       => '-',
    description => 'Chesapeake Bay Land Change Model 2040 prediction of deregulation policy'
  },
  {
    name        => 'land_change_2020_redevelopment_chesapeake',
    units       => '-',
    description => 'Chesapeake Bay Land Change Model 2020 prediction of redevelopment policy'
  },
  {
    name        => 'land_change_2030_redevelopment_chesapeake',
    units       => '-',
    description => 'Chesapeake Bay Land Change Model 2030 prediction of redevelopment policy'
  },
  {
    name        => 'land_change_2040_redevelopment_chesapeake',
    units       => '-',
    description => 'Chesapeake Bay Land Change Model 2040 prediction of redevelopment policy'
  },
  {
    name        => 'land_change_2020_conservation_chesapeake',
    units       => '-',
    description => 'Chesapeake Bay Land Change Model 2020 prediction of conservation policy'
  },
  {
    name        => 'land_change_2030_conservation_chesapeake',
    units       => '-',
    description => 'Chesapeake Bay Land Change Model 2030 prediction of conservation policy'
  },
  {
    name        => 'land_change_2040_conservation_2040_chesapeake',
    units       => '-',
    description => 'Chesapeake Bay Land Change Model 2040 prediction of conservation policy'
  },
  {
    name        => 'stream_temperature_point_08_lower_columbia_river',
    units       => '-',
    description => 'Stream temperature point August Lower Columbia River'
  },
  {
    name        => 'stream_temperature_point_08_middle_columbia_river',
    units       => '-',
    description => 'Stream temperature point August Middle Columbia River'
  },
  {
    name        => 'stream_temperature_point_median_07_new_england',
    units       => '-',
    description => 'Stream temperature point median July New England'
  },
  {
    name        => 'stream_temperature_point_median_08_new_england',
    units       => '-',
    description => 'Stream temperature point median August New England'
  },
  {
    name        => 'stream_temperature_point_mean_daily_range_07_new_england',
    units       => '-',
    description => 'Stream temperature point mean daily_range July New England'
  },
  {
    name        => 'stream_temperature_point_mean_daily_range_08_new_england',
    units       => '-',
    description => 'Stream temperature point mean daily range August New England'
  },
  {
    name        => 'stream_temperature_point_maximum_daily_increase_07_new_england',
    units       => '-',
    description => 'Stream temperature point maximum daily increase July New England'
  },
  {
    name        => 'stream_temperature_point_maximum_daily_increase_08_new_england',
    units       => '-',
    description => 'Stream temperature point maximum daily increase August New England'
  },
  {
    name        => 'stream_temperature_point_maximum_daily_decrease_07_new_england',
    units       => '-',
    description => 'Stream temperature point maximum daily decrease July New England'
  },
  {
    name        => 'stream_temperature_point_maximum_daily_decrease_08_new_england',
    units       => '-',
    description => 'Stream temperature point maximum daily decrease August New England'
  },
  {
    name        => 'stream_temperature_point_maximum_new_england',
    units       => '-',
    description => 'Stream temperature point maximum New England'
  },
  {
    name        => 'stream_temperature_point_day_of_maximum_new_england',
    units       => '-',
    description => 'Stream temperature point day of maximum New England'
  },
  {
    name        => 'stream_temperature_line_08_lower_columbia_river',
    units       => '-',
    description => 'Stream temperature line August Lower Columbia River'
  },
  {
    name        => 'stream_temperature_line_08_middle_columbia_river',
    units       => '-',
    description => 'Stream temperature line August Middle Columbia River'
  },
  {
    name        => 'stream_temperature_line_current_08_upper_rogue_river',
    units       => '-',
    description => 'Stream temperature line current August Upper Rogue River'
  },
  {
    name        => 'stream_temperature_line_future_08_upper_rogue_river',
    units       => '-',
    description => 'Stream temperature line future August Upper Rogue River'
  },
  {
    name        => 'stream_temperature_line_07_new_england',
    units       => '-',
    description => 'Stream temperature line July New England'
  },
  {
    name        => 'stream_temperature_line_meduxnekeag_river',
    units       => '-',
    description => 'Stream temperature line Meduxnekeag River'
  },
  {
    name        => 'stream_shade_line_08_lower_columbia_river',
    units       => '-',
    description => 'Stream shade line August Lower Columbia River'
  },
  {
    name        => 'stream_shade_line_08_middle_columbia_river',
    units       => '-',
    description => 'Stream shade line August Middle Columbia River'
  },
  {
    name        => 'seagrass_point',
    units       => '-',
    description => 'Seagrass point data'
  },
  {
    name        => 'seagrass_point2',
    units       => '-',
    description => 'Seagrass point representation of tiny polygon area data'
  },
  {
    name        => 'seagrass_polygon',
    units       => '-',
    description => 'Seagrass polygon data'
  },
  {
    name        => 'esat1_water_quality',
    units       => '-',
    description => 'Estuary satellite-matched buoy water quality data.'
  },
  {
    name        => 'esat2_water_quality',
    units       => '-',
    description => 'Estuary satellite-matched buoy water quality data'
  },
  {
    name        => 'hspf_blackstone10_huc10_point.hydrology',
    units       => '-',
    description => 'HSPF Blackstone10 domain HUC10 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_blackstone10_huc10_polygon.hydrology',
    units       => '-',
    description => 'HSPF Blackstone10 domain HUC10 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_blackstone10_huc12_point.hydrology',
    units       => '-',
    description => 'HSPF Blackstone10 domain HUC12 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_blackstone10_huc12_polygon.hydrology',
    units       => '-',
    description => 'HSPF Blackstone10 domain HUC12 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_blackstone1_huc10_point.hydrology',
    units       => '-',
    description => 'HSPF Blackstone1 domain HUC10 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_blackstone1_huc10_polygon.hydrology',
    units       => '-',
    description => 'HSPF Blackstone1 domain HUC10 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_blackstone1_huc12_point.hydrology',
    units       => '-',
    description => 'HSPF Blackstone1 domain HUC12 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_blackstone1_huc12_polygon.hydrology',
    units       => '-',
    description => 'HSPF Blackstone1 domain HUC12 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_blackstone2_huc10_point.hydrology',
    units       => '-',
    description => 'HSPF Blackstone2 domain HUC10 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_blackstone2_huc10_polygon.hydrology',
    units       => '-',
    description => 'HSPF Blackstone2 domain HUC10 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_blackstone2_huc12_point.hydrology',
    units       => '-',
    description => 'HSPF Blackstone2 domain HUC12 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_blackstone2_huc12_polygon.hydrology',
    units       => '-',
    description => 'HSPF Blackstone2 domain HUC12 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_blackstone34_huc10_point.hydrology',
    units       => '-',
    description => 'HSPF Blackstone34 domain HUC10 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_blackstone34_huc10_polygon.hydrology',
    units       => '-',
    description => 'HSPF Blackstone34 domain HUC10 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_blackstone34_huc12_point.hydrology',
    units       => '-',
    description => 'HSPF Blackstone34 domain HUC12 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_blackstone34_huc12_polygon.hydrology',
    units       => '-',
    description => 'HSPF Blackstone34 domain HUC12 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_blackstone35_huc10_point.hydrology',
    units       => '-',
    description => 'HSPF Blackstone35 domain HUC10 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_blackstone35_huc10_polygon.hydrology',
    units       => '-',
    description => 'HSPF Blackstone35 domain HUC10 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_blackstone35_huc12_point.hydrology',
    units       => '-',
    description => 'HSPF Blackstone35 domain HUC12 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_blackstone35_huc12_polygon.hydrology',
    units       => '-',
    description => 'HSPF Blackstone35 domain HUC12 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_blackstone36_huc10_point.hydrology',
    units       => '-',
    description => 'HSPF Blackstone36 domain HUC10 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_blackstone36_huc10_polygon.hydrology',
    units       => '-',
    description => 'HSPF Blackstone36 domain HUC10 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_blackstone36_huc12_point.hydrology',
    units       => '-',
    description => 'HSPF Blackstone36 domain HUC12 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_blackstone36_huc12_polygon.hydrology',
    units       => '-',
    description => 'HSPF Blackstone36 domain HUC12 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_blackstone37_huc10_point.hydrology',
    units       => '-',
    description => 'HSPF Blackstone37 domain HUC10 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_blackstone37_huc10_polygon.hydrology',
    units       => '-',
    description => 'HSPF Blackstone37 domain HUC10 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_blackstone37_huc12_point.hydrology',
    units       => '-',
    description => 'HSPF Blackstone37 domain HUC12 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_blackstone37_huc12_polygon.hydrology',
    units       => '-',
    description => 'HSPF Blackstone37 domain HUC12 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_blackstone38_huc10_point.hydrology',
    units       => '-',
    description => 'HSPF Blackstone38 domain HUC10 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_blackstone38_huc10_polygon.hydrology',
    units       => '-',
    description => 'HSPF Blackstone38 domain HUC10 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_blackstone38_huc12_point.hydrology',
    units       => '-',
    description => 'HSPF Blackstone38 domain HUC12 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_blackstone38_huc12_polygon.hydrology',
    units       => '-',
    description => 'HSPF Blackstone38 domain HUC12 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_blackstone39_huc10_point.hydrology',
    units       => '-',
    description => 'HSPF Blackstone39 domain HUC10 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_blackstone39_huc10_polygon.hydrology',
    units       => '-',
    description => 'HSPF Blackstone39 domain HUC10 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_blackstone39_huc12_point.hydrology',
    units       => '-',
    description => 'HSPF Blackstone39 domain HUC12 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_blackstone39_huc12_polygon.hydrology',
    units       => '-',
    description => 'HSPF Blackstone39 domain HUC12 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_blackstone40_huc10_point.hydrology',
    units       => '-',
    description => 'HSPF Blackstone40 domain HUC10 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_blackstone40_huc10_polygon.hydrology',
    units       => '-',
    description => 'HSPF Blackstone40 domain HUC10 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_blackstone40_huc12_point.hydrology',
    units       => '-',
    description => 'HSPF Blackstone40 domain HUC12 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_blackstone40_huc12_polygon.hydrology',
    units       => '-',
    description => 'HSPF Blackstone40 domain HUC12 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_blackstone41_huc10_point.hydrology',
    units       => '-',
    description => 'HSPF Blackstone41 domain HUC10 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_blackstone41_huc10_polygon.hydrology',
    units       => '-',
    description => 'HSPF Blackstone41 domain HUC10 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_blackstone41_huc12_point.hydrology',
    units       => '-',
    description => 'HSPF Blackstone41 domain HUC12 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_blackstone41_huc12_polygon.hydrology',
    units       => '-',
    description => 'HSPF Blackstone41 domain HUC12 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_blackstone42_huc10_point.hydrology',
    units       => '-',
    description => 'HSPF Blackstone42 domain HUC10 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_blackstone42_huc10_polygon.hydrology',
    units       => '-',
    description => 'HSPF Blackstone42 domain HUC10 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_blackstone42_huc12_point.hydrology',
    units       => '-',
    description => 'HSPF Blackstone42 domain HUC12 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_blackstone42_huc12_polygon.hydrology',
    units       => '-',
    description => 'HSPF Blackstone42 domain HUC12 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_blackstone43_huc10_point.hydrology',
    units       => '-',
    description => 'HSPF Blackstone43 domain HUC10 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_blackstone43_huc10_polygon.hydrology',
    units       => '-',
    description => 'HSPF Blackstone43 domain HUC10 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_blackstone43_huc12_point.hydrology',
    units       => '-',
    description => 'HSPF Blackstone43 domain HUC12 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_blackstone43_huc12_polygon.hydrology',
    units       => '-',
    description => 'HSPF Blackstone43 domain HUC12 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_blackstone9_huc10_point.hydrology',
    units       => '-',
    description => 'HSPF Blackstone9 domain HUC10 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_blackstone9_huc10_polygon.hydrology',
    units       => '-',
    description => 'HSPF Blackstone9 domain HUC10 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_blackstone9_huc12_point.hydrology',
    units       => '-',
    description => 'HSPF Blackstone9 domain HUC12 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_blackstone9_huc12_polygon.hydrology',
    units       => '-',
    description => 'HSPF Blackstone9 domain HUC12 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_charles11_huc10_point.hydrology',
    units       => '-',
    description => 'HSPF Charles11 domain HUC10 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_charles11_huc10_point.recharge_adjustment',
    units       => '-',
    description => 'HSPF Charles11 domain HUC10 point recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_charles11_huc10_polygon.hydrology',
    units       => '-',
    description => 'HSPF Charles11 domain HUC10 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_charles11_huc10_polygon.recharge_adjustment',
    units       => '-',
    description => 'HSPF Charles11 domain HUC10 polygon recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_charles11_huc12_point.hydrology',
    units       => '-',
    description => 'HSPF Charles11 domain HUC12 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_charles11_huc12_point.recharge_adjustment',
    units       => '-',
    description => 'HSPF Charles11 domain HUC12 point recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_charles11_huc12_polygon.hydrology',
    units       => '-',
    description => 'HSPF Charles11 domain HUC12 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_charles11_huc12_polygon.recharge_adjustment',
    units       => '-',
    description => 'HSPF Charles11 domain HUC12 polygon recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_charles3_huc10_point.hydrology',
    units       => '-',
    description => 'HSPF Charles3 domain HUC10 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_charles3_huc10_point.recharge_adjustment',
    units       => '-',
    description => 'HSPF Charles3 domain HUC10 point recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_charles3_huc10_polygon.hydrology',
    units       => '-',
    description => 'HSPF Charles3 domain HUC10 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_charles3_huc10_polygon.recharge_adjustment',
    units       => '-',
    description => 'HSPF Charles3 domain HUC10 polygon recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_charles3_huc12_point.hydrology',
    units       => '-',
    description => 'HSPF Charles3 domain HUC12 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_charles3_huc12_point.recharge_adjustment',
    units       => '-',
    description => 'HSPF Charles3 domain HUC12 point recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_charles3_huc12_polygon.hydrology',
    units       => '-',
    description => 'HSPF Charles3 domain HUC12 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_charles3_huc12_polygon.recharge_adjustment',
    units       => '-',
    description => 'HSPF Charles3 domain HUC12 polygon recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_charles44_huc10_point.hydrology',
    units       => '-',
    description => 'HSPF Charles44 domain HUC10 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_charles44_huc10_point.recharge_adjustment',
    units       => '-',
    description => 'HSPF Charles44 domain HUC10 point recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_charles44_huc10_polygon.hydrology',
    units       => '-',
    description => 'HSPF Charles44 domain HUC10 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_charles44_huc10_polygon.recharge_adjustment',
    units       => '-',
    description => 'HSPF Charles44 domain HUC10 polygon recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_charles44_huc12_point.hydrology',
    units       => '-',
    description => 'HSPF Charles44 domain HUC12 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_charles44_huc12_point.recharge_adjustment',
    units       => '-',
    description => 'HSPF Charles44 domain HUC12 point recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_charles44_huc12_polygon.hydrology',
    units       => '-',
    description => 'HSPF Charles44 domain HUC12 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_charles44_huc12_polygon.recharge_adjustment',
    units       => '-',
    description => 'HSPF Charles44 domain HUC12 polygon recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_charles45_huc10_point.hydrology',
    units       => '-',
    description => 'HSPF Charles45 domain HUC10 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_charles45_huc10_point.recharge_adjustment',
    units       => '-',
    description => 'HSPF Charles45 domain HUC10 point recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_charles45_huc10_polygon.hydrology',
    units       => '-',
    description => 'HSPF Charles45 domain HUC10 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_charles45_huc10_polygon.recharge_adjustment',
    units       => '-',
    description => 'HSPF Charles45 domain HUC10 polygon recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_charles45_huc12_point.hydrology',
    units       => '-',
    description => 'HSPF Charles45 domain HUC12 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_charles45_huc12_point.recharge_adjustment',
    units       => '-',
    description => 'HSPF Charles45 domain HUC12 point recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_charles45_huc12_polygon.hydrology',
    units       => '-',
    description => 'HSPF Charles45 domain HUC12 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_charles45_huc12_polygon.recharge_adjustment',
    units       => '-',
    description => 'HSPF Charles45 domain HUC12 polygon recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_charles46_huc10_point.hydrology',
    units       => '-',
    description => 'HSPF Charles46 domain HUC10 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_charles46_huc10_point.recharge_adjustment',
    units       => '-',
    description => 'HSPF Charles46 domain HUC10 point recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_charles46_huc10_polygon.hydrology',
    units       => '-',
    description => 'HSPF Charles46 domain HUC10 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_charles46_huc10_polygon.recharge_adjustment',
    units       => '-',
    description => 'HSPF Charles46 domain HUC10 polygon recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_charles46_huc12_point.hydrology',
    units       => '-',
    description => 'HSPF Charles46 domain HUC12 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_charles46_huc12_point.recharge_adjustment',
    units       => '-',
    description => 'HSPF Charles46 domain HUC12 point recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_charles46_huc12_polygon.hydrology',
    units       => '-',
    description => 'HSPF Charles46 domain HUC12 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_charles46_huc12_polygon.recharge_adjustment',
    units       => '-',
    description => 'HSPF Charles46 domain HUC12 polygon recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_charles47_huc10_point.hydrology',
    units       => '-',
    description => 'HSPF Charles47 domain HUC10 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_charles47_huc10_point.recharge_adjustment',
    units       => '-',
    description => 'HSPF Charles47 domain HUC10 point recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_charles47_huc10_polygon.hydrology',
    units       => '-',
    description => 'HSPF Charles47 domain HUC10 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_charles47_huc10_polygon.recharge_adjustment',
    units       => '-',
    description => 'HSPF Charles47 domain HUC10 polygon recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_charles47_huc12_point.hydrology',
    units       => '-',
    description => 'HSPF Charles47 domain HUC12 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_charles47_huc12_point.recharge_adjustment',
    units       => '-',
    description => 'HSPF Charles47 domain HUC12 point recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_charles47_huc12_polygon.hydrology',
    units       => '-',
    description => 'HSPF Charles47 domain HUC12 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_charles47_huc12_polygon.recharge_adjustment',
    units       => '-',
    description => 'HSPF Charles47 domain HUC12 polygon recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_charles48_huc10_point.hydrology',
    units       => '-',
    description => 'HSPF Charles48 domain HUC10 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_charles48_huc10_point.recharge_adjustment',
    units       => '-',
    description => 'HSPF Charles48 domain HUC10 point recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_charles48_huc10_polygon.hydrology',
    units       => '-',
    description => 'HSPF Charles48 domain HUC10 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_charles48_huc10_polygon.recharge_adjustment',
    units       => '-',
    description => 'HSPF Charles48 domain HUC10 polygon recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_charles48_huc12_point.hydrology',
    units       => '-',
    description => 'HSPF Charles48 domain HUC12 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_charles48_huc12_point.recharge_adjustment',
    units       => '-',
    description => 'HSPF Charles48 domain HUC12 point recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_charles48_huc12_polygon.hydrology',
    units       => '-',
    description => 'HSPF Charles48 domain HUC12 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_charles48_huc12_polygon.recharge_adjustment',
    units       => '-',
    description => 'HSPF Charles48 domain HUC12 polygon recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_farmington12_huc10_point.hydrology',
    units       => '-',
    description => 'HSPF Farmington12 domain HUC10 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_farmington12_huc10_point.recharge_adjustment',
    units       => '-',
    description => 'HSPF Farmington12 domain HUC10 point recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_farmington12_huc10_polygon.hydrology',
    units       => '-',
    description => 'HSPF Farmington12 domain HUC10 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_farmington12_huc10_polygon.recharge_adjustment',
    units       => '-',
    description => 'HSPF Farmington12 domain HUC10 polygon recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_farmington12_huc12_point.hydrology',
    units       => '-',
    description => 'HSPF Farmington12 domain HUC12 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_farmington12_huc12_point.recharge_adjustment',
    units       => '-',
    description => 'HSPF Farmington12 domain HUC12 point recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_farmington12_huc12_polygon.hydrology',
    units       => '-',
    description => 'HSPF Farmington12 domain HUC12 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_farmington12_huc12_polygon.recharge_adjustment',
    units       => '-',
    description => 'HSPF Farmington12 domain HUC12 polygon recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_farmington49_huc10_point.hydrology',
    units       => '-',
    description => 'HSPF Farmington49 domain HUC10 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_farmington49_huc10_point.recharge_adjustment',
    units       => '-',
    description => 'HSPF Farmington49 domain HUC10 point recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_farmington49_huc10_polygon.hydrology',
    units       => '-',
    description => 'HSPF Farmington49 domain HUC10 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_farmington49_huc10_polygon.recharge_adjustment',
    units       => '-',
    description => 'HSPF Farmington49 domain HUC10 polygon recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_farmington49_huc12_point.hydrology',
    units       => '-',
    description => 'HSPF Farmington49 domain HUC12 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_farmington49_huc12_point.recharge_adjustment',
    units       => '-',
    description => 'HSPF Farmington49 domain HUC12 point recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_farmington49_huc12_polygon.hydrology',
    units       => '-',
    description => 'HSPF Farmington49 domain HUC12 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_farmington49_huc12_polygon.recharge_adjustment',
    units       => '-',
    description => 'HSPF Farmington49 domain HUC12 polygon recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_farmington4_huc10_point.hydrology',
    units       => '-',
    description => 'HSPF Farmington4 domain HUC10 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_farmington4_huc10_point.recharge_adjustment',
    units       => '-',
    description => 'HSPF Farmington4 domain HUC10 point recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_farmington4_huc10_polygon.hydrology',
    units       => '-',
    description => 'HSPF Farmington4 domain HUC10 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_farmington4_huc10_polygon.recharge_adjustment',
    units       => '-',
    description => 'HSPF Farmington4 domain HUC10 polygon recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_farmington4_huc12_point.hydrology',
    units       => '-',
    description => 'HSPF Farmington4 domain HUC12 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_farmington4_huc12_point.recharge_adjustment',
    units       => '-',
    description => 'HSPF Farmington4 domain HUC12 point recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_farmington4_huc12_polygon.hydrology',
    units       => '-',
    description => 'HSPF Farmington4 domain HUC12 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_farmington4_huc12_polygon.recharge_adjustment',
    units       => '-',
    description => 'HSPF Farmington4 domain HUC12 polygon recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_farmington50_huc10_point.hydrology',
    units       => '-',
    description => 'HSPF Farmington50 domain HUC10 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_farmington50_huc10_point.recharge_adjustment',
    units       => '-',
    description => 'HSPF Farmington50 domain HUC10 point recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_farmington50_huc10_polygon.hydrology',
    units       => '-',
    description => 'HSPF Farmington50 domain HUC10 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_farmington50_huc10_polygon.recharge_adjustment',
    units       => '-',
    description => 'HSPF Farmington50 domain HUC10 polygon recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_farmington50_huc12_point.hydrology',
    units       => '-',
    description => 'HSPF Farmington50 domain HUC12 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_farmington50_huc12_point.recharge_adjustment',
    units       => '-',
    description => 'HSPF Farmington50 domain HUC12 point recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_farmington50_huc12_polygon.hydrology',
    units       => '-',
    description => 'HSPF Farmington50 domain HUC12 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_farmington50_huc12_polygon.recharge_adjustment',
    units       => '-',
    description => 'HSPF Farmington50 domain HUC12 polygon recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_farmington51_huc10_point.hydrology',
    units       => '-',
    description => 'HSPF Farmington51 domain HUC10 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_farmington51_huc10_point.recharge_adjustment',
    units       => '-',
    description => 'HSPF Farmington51 domain HUC10 point recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_farmington51_huc10_polygon.hydrology',
    units       => '-',
    description => 'HSPF Farmington51 domain HUC10 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_farmington51_huc10_polygon.recharge_adjustment',
    units       => '-',
    description => 'HSPF Farmington51 domain HUC10 polygon recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_farmington51_huc12_point.hydrology',
    units       => '-',
    description => 'HSPF Farmington51 domain HUC12 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_farmington51_huc12_point.recharge_adjustment',
    units       => '-',
    description => 'HSPF Farmington51 domain HUC12 point recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_farmington51_huc12_polygon.hydrology',
    units       => '-',
    description => 'HSPF Farmington51 domain HUC12 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_farmington51_huc12_polygon.recharge_adjustment',
    units       => '-',
    description => 'HSPF Farmington51 domain HUC12 polygon recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_farmington52_huc10_point.hydrology',
    units       => '-',
    description => 'HSPF Farmington52 domain HUC10 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_farmington52_huc10_point.recharge_adjustment',
    units       => '-',
    description => 'HSPF Farmington52 domain HUC10 point recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_farmington52_huc10_polygon.hydrology',
    units       => '-',
    description => 'HSPF Farmington52 domain HUC10 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_farmington52_huc10_polygon.recharge_adjustment',
    units       => '-',
    description => 'HSPF Farmington52 domain HUC10 polygon recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_farmington52_huc12_point.hydrology',
    units       => '-',
    description => 'HSPF Farmington52 domain HUC12 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_farmington52_huc12_point.recharge_adjustment',
    units       => '-',
    description => 'HSPF Farmington52 domain HUC12 point recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_farmington52_huc12_polygon.hydrology',
    units       => '-',
    description => 'HSPF Farmington52 domain HUC12 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_farmington52_huc12_polygon.recharge_adjustment',
    units       => '-',
    description => 'HSPF Farmington52 domain HUC12 polygon recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_farmington53_huc10_point.hydrology',
    units       => '-',
    description => 'HSPF Farmington53 domain HUC10 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_farmington53_huc10_point.recharge_adjustment',
    units       => '-',
    description => 'HSPF Farmington53 domain HUC10 point recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_farmington53_huc10_polygon.hydrology',
    units       => '-',
    description => 'HSPF Farmington53 domain HUC10 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_farmington53_huc10_polygon.recharge_adjustment',
    units       => '-',
    description => 'HSPF Farmington53 domain HUC10 polygon recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_farmington53_huc12_point.hydrology',
    units       => '-',
    description => 'HSPF Farmington53 domain HUC12 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_farmington53_huc12_point.recharge_adjustment',
    units       => '-',
    description => 'HSPF Farmington53 domain HUC12 point recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_farmington53_huc12_polygon.hydrology',
    units       => '-',
    description => 'HSPF Farmington53 domain HUC12 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_farmington53_huc12_polygon.recharge_adjustment',
    units       => '-',
    description => 'HSPF Farmington53 domain HUC12 polygon recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_ipswich13_huc10_point.hydrology',
    units       => '-',
    description => 'HSPF Ipswich13 domain HUC10 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_ipswich13_huc10_point.recharge_adjustment',
    units       => '-',
    description => 'HSPF Ipswich13 domain HUC10 point recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_ipswich13_huc10_polygon.hydrology',
    units       => '-',
    description => 'HSPF Ipswich13 domain HUC10 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_ipswich13_huc10_polygon.recharge_adjustment',
    units       => '-',
    description => 'HSPF Ipswich13 domain HUC10 polygon recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_ipswich13_huc12_point.hydrology',
    units       => '-',
    description => 'HSPF Ipswich13 domain HUC12 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_ipswich13_huc12_point.recharge_adjustment',
    units       => '-',
    description => 'HSPF Ipswich13 domain HUC12 point recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_ipswich13_huc12_polygon.hydrology',
    units       => '-',
    description => 'HSPF Ipswich13 domain HUC12 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_ipswich13_huc12_polygon.recharge_adjustment',
    units       => '-',
    description => 'HSPF Ipswich13 domain HUC12 polygon recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_ipswich54_huc10_point.hydrology',
    units       => '-',
    description => 'HSPF Ipswich54 domain HUC10 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_ipswich54_huc10_point.recharge_adjustment',
    units       => '-',
    description => 'HSPF Ipswich54 domain HUC10 point recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_ipswich54_huc10_polygon.hydrology',
    units       => '-',
    description => 'HSPF Ipswich54 domain HUC10 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_ipswich54_huc10_polygon.recharge_adjustment',
    units       => '-',
    description => 'HSPF Ipswich54 domain HUC10 polygon recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_ipswich54_huc12_point.hydrology',
    units       => '-',
    description => 'HSPF Ipswich54 domain HUC12 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_ipswich54_huc12_point.recharge_adjustment',
    units       => '-',
    description => 'HSPF Ipswich54 domain HUC12 point recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_ipswich54_huc12_polygon.hydrology',
    units       => '-',
    description => 'HSPF Ipswich54 domain HUC12 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_ipswich54_huc12_polygon.recharge_adjustment',
    units       => '-',
    description => 'HSPF Ipswich54 domain HUC12 polygon recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_ipswich55_huc10_point.hydrology',
    units       => '-',
    description => 'HSPF Ipswich55 domain HUC10 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_ipswich55_huc10_point.recharge_adjustment',
    units       => '-',
    description => 'HSPF Ipswich55 domain HUC10 point recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_ipswich55_huc10_polygon.hydrology',
    units       => '-',
    description => 'HSPF Ipswich55 domain HUC10 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_ipswich55_huc10_polygon.recharge_adjustment',
    units       => '-',
    description => 'HSPF Ipswich55 domain HUC10 polygon recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_ipswich55_huc12_point.hydrology',
    units       => '-',
    description => 'HSPF Ipswich55 domain HUC12 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_ipswich55_huc12_point.recharge_adjustment',
    units       => '-',
    description => 'HSPF Ipswich55 domain HUC12 point recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_ipswich55_huc12_polygon.hydrology',
    units       => '-',
    description => 'HSPF Ipswich55 domain HUC12 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_ipswich55_huc12_polygon.recharge_adjustment',
    units       => '-',
    description => 'HSPF Ipswich55 domain HUC12 polygon recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_ipswich56_huc10_point.hydrology',
    units       => '-',
    description => 'HSPF Ipswich56 domain HUC10 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_ipswich56_huc10_point.recharge_adjustment',
    units       => '-',
    description => 'HSPF Ipswich56 domain HUC10 point recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_ipswich56_huc10_polygon.hydrology',
    units       => '-',
    description => 'HSPF Ipswich56 domain HUC10 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_ipswich56_huc10_polygon.recharge_adjustment',
    units       => '-',
    description => 'HSPF Ipswich56 domain HUC10 polygon recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_ipswich56_huc12_point.hydrology',
    units       => '-',
    description => 'HSPF Ipswich56 domain HUC12 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_ipswich56_huc12_point.recharge_adjustment',
    units       => '-',
    description => 'HSPF Ipswich56 domain HUC12 point recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_ipswich56_huc12_polygon.hydrology',
    units       => '-',
    description => 'HSPF Ipswich56 domain HUC12 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_ipswich56_huc12_polygon.recharge_adjustment',
    units       => '-',
    description => 'HSPF Ipswich56 domain HUC12 polygon recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_ipswich57_huc10_point.hydrology',
    units       => '-',
    description => 'HSPF Ipswich57 domain HUC10 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_ipswich57_huc10_point.recharge_adjustment',
    units       => '-',
    description => 'HSPF Ipswich57 domain HUC10 point recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_ipswich57_huc10_polygon.hydrology',
    units       => '-',
    description => 'HSPF Ipswich57 domain HUC10 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_ipswich57_huc10_polygon.recharge_adjustment',
    units       => '-',
    description => 'HSPF Ipswich57 domain HUC10 polygon recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_ipswich57_huc12_point.hydrology',
    units       => '-',
    description => 'HSPF Ipswich57 domain HUC12 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_ipswich57_huc12_point.recharge_adjustment',
    units       => '-',
    description => 'HSPF Ipswich57 domain HUC12 point recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_ipswich57_huc12_polygon.hydrology',
    units       => '-',
    description => 'HSPF Ipswich57 domain HUC12 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_ipswich57_huc12_polygon.recharge_adjustment',
    units       => '-',
    description => 'HSPF Ipswich57 domain HUC12 polygon recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_ipswich58_huc10_point.hydrology',
    units       => '-',
    description => 'HSPF Ipswich58 domain HUC10 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_ipswich58_huc10_point.recharge_adjustment',
    units       => '-',
    description => 'HSPF Ipswich58 domain HUC10 point recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_ipswich58_huc10_polygon.hydrology',
    units       => '-',
    description => 'HSPF Ipswich58 domain HUC10 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_ipswich58_huc10_polygon.recharge_adjustment',
    units       => '-',
    description => 'HSPF Ipswich58 domain HUC10 polygon recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_ipswich58_huc12_point.hydrology',
    units       => '-',
    description => 'HSPF Ipswich58 domain HUC12 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_ipswich58_huc12_point.recharge_adjustment',
    units       => '-',
    description => 'HSPF Ipswich58 domain HUC12 point recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_ipswich58_huc12_polygon.hydrology',
    units       => '-',
    description => 'HSPF Ipswich58 domain HUC12 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_ipswich58_huc12_polygon.recharge_adjustment',
    units       => '-',
    description => 'HSPF Ipswich58 domain HUC12 polygon recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_ipswich5_huc10_point.hydrology',
    units       => '-',
    description => 'HSPF Ipswich5 domain HUC10 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_ipswich5_huc10_point.recharge_adjustment',
    units       => '-',
    description => 'HSPF Ipswich5 domain HUC10 point recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_ipswich5_huc10_polygon.hydrology',
    units       => '-',
    description => 'HSPF Ipswich5 domain HUC10 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_ipswich5_huc10_polygon.recharge_adjustment',
    units       => '-',
    description => 'HSPF Ipswich5 domain HUC10 polygon recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_ipswich5_huc12_point.hydrology',
    units       => '-',
    description => 'HSPF Ipswich5 domain HUC12 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_ipswich5_huc12_point.recharge_adjustment',
    units       => '-',
    description => 'HSPF Ipswich5 domain HUC12 point recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_ipswich5_huc12_polygon.hydrology',
    units       => '-',
    description => 'HSPF Ipswich5 domain HUC12 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_ipswich5_huc12_polygon.recharge_adjustment',
    units       => '-',
    description => 'HSPF Ipswich5 domain HUC12 polygon recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_pawcatuck14_huc10_point.hydrology',
    units       => '-',
    description => 'HSPF Pawcatuck14 domain HUC10 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_pawcatuck14_huc10_point.recharge_adjustment',
    units       => '-',
    description => 'HSPF Pawcatuck14 domain HUC10 point recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_pawcatuck14_huc10_polygon.hydrology',
    units       => '-',
    description => 'HSPF Pawcatuck14 domain HUC10 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_pawcatuck14_huc10_polygon.recharge_adjustment',
    units       => '-',
    description => 'HSPF Pawcatuck14 domain HUC10 polygon recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_pawcatuck14_huc12_point.hydrology',
    units       => '-',
    description => 'HSPF Pawcatuck14 domain HUC12 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_pawcatuck14_huc12_point.recharge_adjustment',
    units       => '-',
    description => 'HSPF Pawcatuck14 domain HUC12 point recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_pawcatuck14_huc12_polygon.hydrology',
    units       => '-',
    description => 'HSPF Pawcatuck14 domain HUC12 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_pawcatuck14_huc12_polygon.recharge_adjustment',
    units       => '-',
    description => 'HSPF Pawcatuck14 domain HUC12 polygon recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_pawcatuck59_huc10_point.hydrology',
    units       => '-',
    description => 'HSPF Pawcatuck59 domain HUC10 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_pawcatuck59_huc10_point.recharge_adjustment',
    units       => '-',
    description => 'HSPF Pawcatuck59 domain HUC10 point recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_pawcatuck59_huc10_polygon.hydrology',
    units       => '-',
    description => 'HSPF Pawcatuck59 domain HUC10 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_pawcatuck59_huc10_polygon.recharge_adjustment',
    units       => '-',
    description => 'HSPF Pawcatuck59 domain HUC10 polygon recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_pawcatuck59_huc12_point.hydrology',
    units       => '-',
    description => 'HSPF Pawcatuck59 domain HUC12 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_pawcatuck59_huc12_point.recharge_adjustment',
    units       => '-',
    description => 'HSPF Pawcatuck59 domain HUC12 point recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_pawcatuck59_huc12_polygon.hydrology',
    units       => '-',
    description => 'HSPF Pawcatuck59 domain HUC12 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_pawcatuck59_huc12_polygon.recharge_adjustment',
    units       => '-',
    description => 'HSPF Pawcatuck59 domain HUC12 polygon recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_pawcatuck60_huc10_point.hydrology',
    units       => '-',
    description => 'HSPF Pawcatuck60 domain HUC10 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_pawcatuck60_huc10_point.recharge_adjustment',
    units       => '-',
    description => 'HSPF Pawcatuck60 domain HUC10 point recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_pawcatuck60_huc10_polygon.hydrology',
    units       => '-',
    description => 'HSPF Pawcatuck60 domain HUC10 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_pawcatuck60_huc10_polygon.recharge_adjustment',
    units       => '-',
    description => 'HSPF Pawcatuck60 domain HUC10 polygon recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_pawcatuck60_huc12_point.hydrology',
    units       => '-',
    description => 'HSPF Pawcatuck60 domain HUC12 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_pawcatuck60_huc12_point.recharge_adjustment',
    units       => '-',
    description => 'HSPF Pawcatuck60 domain HUC12 point recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_pawcatuck60_huc12_polygon.hydrology',
    units       => '-',
    description => 'HSPF Pawcatuck60 domain HUC12 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_pawcatuck60_huc12_polygon.recharge_adjustment',
    units       => '-',
    description => 'HSPF Pawcatuck60 domain HUC12 polygon recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_pawcatuck61_huc10_point.hydrology',
    units       => '-',
    description => 'HSPF Pawcatuck61 domain HUC10 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_pawcatuck61_huc10_point.recharge_adjustment',
    units       => '-',
    description => 'HSPF Pawcatuck61 domain HUC10 point recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_pawcatuck61_huc10_polygon.hydrology',
    units       => '-',
    description => 'HSPF Pawcatuck61 domain HUC10 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_pawcatuck61_huc10_polygon.recharge_adjustment',
    units       => '-',
    description => 'HSPF Pawcatuck61 domain HUC10 polygon recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_pawcatuck61_huc12_point.hydrology',
    units       => '-',
    description => 'HSPF Pawcatuck61 domain HUC12 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_pawcatuck61_huc12_point.recharge_adjustment',
    units       => '-',
    description => 'HSPF Pawcatuck61 domain HUC12 point recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_pawcatuck61_huc12_polygon.hydrology',
    units       => '-',
    description => 'HSPF Pawcatuck61 domain HUC12 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_pawcatuck61_huc12_polygon.recharge_adjustment',
    units       => '-',
    description => 'HSPF Pawcatuck61 domain HUC12 polygon recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_pawcatuck62_huc10_point.hydrology',
    units       => '-',
    description => 'HSPF Pawcatuck62 domain HUC10 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_pawcatuck62_huc10_point.recharge_adjustment',
    units       => '-',
    description => 'HSPF Pawcatuck62 domain HUC10 point recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_pawcatuck62_huc10_polygon.hydrology',
    units       => '-',
    description => 'HSPF Pawcatuck62 domain HUC10 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_pawcatuck62_huc10_polygon.recharge_adjustment',
    units       => '-',
    description => 'HSPF Pawcatuck62 domain HUC10 polygon recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_pawcatuck62_huc12_point.hydrology',
    units       => '-',
    description => 'HSPF Pawcatuck62 domain HUC12 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_pawcatuck62_huc12_point.recharge_adjustment',
    units       => '-',
    description => 'HSPF Pawcatuck62 domain HUC12 point recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_pawcatuck62_huc12_polygon.hydrology',
    units       => '-',
    description => 'HSPF Pawcatuck62 domain HUC12 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_pawcatuck62_huc12_polygon.recharge_adjustment',
    units       => '-',
    description => 'HSPF Pawcatuck62 domain HUC12 polygon recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_pawcatuck63_huc10_point.hydrology',
    units       => '-',
    description => 'HSPF Pawcatuck63 domain HUC10 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_pawcatuck63_huc10_point.recharge_adjustment',
    units       => '-',
    description => 'HSPF Pawcatuck63 domain HUC10 point recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_pawcatuck63_huc10_polygon.hydrology',
    units       => '-',
    description => 'HSPF Pawcatuck63 domain HUC10 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_pawcatuck63_huc10_polygon.recharge_adjustment',
    units       => '-',
    description => 'HSPF Pawcatuck63 domain HUC10 polygon recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_pawcatuck63_huc12_point.hydrology',
    units       => '-',
    description => 'HSPF Pawcatuck63 domain HUC12 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_pawcatuck63_huc12_point.recharge_adjustment',
    units       => '-',
    description => 'HSPF Pawcatuck63 domain HUC12 point recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_pawcatuck63_huc12_polygon.hydrology',
    units       => '-',
    description => 'HSPF Pawcatuck63 domain HUC12 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_pawcatuck63_huc12_polygon.recharge_adjustment',
    units       => '-',
    description => 'HSPF Pawcatuck63 domain HUC12 polygon recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_pawcatuck6_huc10_point.hydrology',
    units       => '-',
    description => 'HSPF Pawcatuck6 domain HUC10 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_pawcatuck6_huc10_point.recharge_adjustment',
    units       => '-',
    description => 'HSPF Pawcatuck6 domain HUC10 point recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_pawcatuck6_huc10_polygon.hydrology',
    units       => '-',
    description => 'HSPF Pawcatuck6 domain HUC10 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_pawcatuck6_huc10_polygon.recharge_adjustment',
    units       => '-',
    description => 'HSPF Pawcatuck6 domain HUC10 polygon recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_pawcatuck6_huc12_point.hydrology',
    units       => '-',
    description => 'HSPF Pawcatuck6 domain HUC12 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_pawcatuck6_huc12_point.recharge_adjustment',
    units       => '-',
    description => 'HSPF Pawcatuck6 domain HUC12 point recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_pawcatuck6_huc12_polygon.hydrology',
    units       => '-',
    description => 'HSPF Pawcatuck6 domain HUC12 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_pawcatuck6_huc12_polygon.recharge_adjustment',
    units       => '-',
    description => 'HSPF Pawcatuck6 domain HUC12 polygon recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_sudbury15_huc10_point.hydrology',
    units       => '-',
    description => 'HSPF Sudbury15 domain HUC10 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_sudbury15_huc10_point.recharge_adjustment',
    units       => '-',
    description => 'HSPF Sudbury15 domain HUC10 point recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_sudbury15_huc10_polygon.hydrology',
    units       => '-',
    description => 'HSPF Sudbury15 domain HUC10 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_sudbury15_huc10_polygon.recharge_adjustment',
    units       => '-',
    description => 'HSPF Sudbury15 domain HUC10 polygon recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_sudbury15_huc12_point.hydrology',
    units       => '-',
    description => 'HSPF Sudbury15 domain HUC12 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_sudbury15_huc12_point.recharge_adjustment',
    units       => '-',
    description => 'HSPF Sudbury15 domain HUC12 point recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_sudbury15_huc12_polygon.hydrology',
    units       => '-',
    description => 'HSPF Sudbury15 domain HUC12 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_sudbury15_huc12_polygon.recharge_adjustment',
    units       => '-',
    description => 'HSPF Sudbury15 domain HUC12 polygon recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_sudbury64_huc10_point.hydrology',
    units       => '-',
    description => 'HSPF Sudbury64 domain HUC10 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_sudbury64_huc10_point.recharge_adjustment',
    units       => '-',
    description => 'HSPF Sudbury64 domain HUC10 point recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_sudbury64_huc10_polygon.hydrology',
    units       => '-',
    description => 'HSPF Sudbury64 domain HUC10 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_sudbury64_huc10_polygon.recharge_adjustment',
    units       => '-',
    description => 'HSPF Sudbury64 domain HUC10 polygon recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_sudbury64_huc12_point.hydrology',
    units       => '-',
    description => 'HSPF Sudbury64 domain HUC12 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_sudbury64_huc12_point.recharge_adjustment',
    units       => '-',
    description => 'HSPF Sudbury64 domain HUC12 point recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_sudbury64_huc12_polygon.hydrology',
    units       => '-',
    description => 'HSPF Sudbury64 domain HUC12 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_sudbury64_huc12_polygon.recharge_adjustment',
    units       => '-',
    description => 'HSPF Sudbury64 domain HUC12 polygon recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_sudbury65_huc10_point.hydrology',
    units       => '-',
    description => 'HSPF Sudbury65 domain HUC10 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_sudbury65_huc10_point.recharge_adjustment',
    units       => '-',
    description => 'HSPF Sudbury65 domain HUC10 point recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_sudbury65_huc10_polygon.hydrology',
    units       => '-',
    description => 'HSPF Sudbury65 domain HUC10 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_sudbury65_huc10_polygon.recharge_adjustment',
    units       => '-',
    description => 'HSPF Sudbury65 domain HUC10 polygon recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_sudbury65_huc12_point.hydrology',
    units       => '-',
    description => 'HSPF Sudbury65 domain HUC12 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_sudbury65_huc12_point.recharge_adjustment',
    units       => '-',
    description => 'HSPF Sudbury65 domain HUC12 point recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_sudbury65_huc12_polygon.hydrology',
    units       => '-',
    description => 'HSPF Sudbury65 domain HUC12 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_sudbury65_huc12_polygon.recharge_adjustment',
    units       => '-',
    description => 'HSPF Sudbury65 domain HUC12 polygon recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_sudbury66_huc10_point.hydrology',
    units       => '-',
    description => 'HSPF Sudbury66 domain HUC10 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_sudbury66_huc10_point.recharge_adjustment',
    units       => '-',
    description => 'HSPF Sudbury66 domain HUC10 point recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_sudbury66_huc10_polygon.hydrology',
    units       => '-',
    description => 'HSPF Sudbury66 domain HUC10 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_sudbury66_huc10_polygon.recharge_adjustment',
    units       => '-',
    description => 'HSPF Sudbury66 domain HUC10 polygon recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_sudbury66_huc12_point.hydrology',
    units       => '-',
    description => 'HSPF Sudbury66 domain HUC12 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_sudbury66_huc12_point.recharge_adjustment',
    units       => '-',
    description => 'HSPF Sudbury66 domain HUC12 point recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_sudbury66_huc12_polygon.hydrology',
    units       => '-',
    description => 'HSPF Sudbury66 domain HUC12 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_sudbury66_huc12_polygon.recharge_adjustment',
    units       => '-',
    description => 'HSPF Sudbury66 domain HUC12 polygon recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_sudbury67_huc10_point.hydrology',
    units       => '-',
    description => 'HSPF Sudbury67 domain HUC10 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_sudbury67_huc10_point.recharge_adjustment',
    units       => '-',
    description => 'HSPF Sudbury67 domain HUC10 point recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_sudbury67_huc10_polygon.hydrology',
    units       => '-',
    description => 'HSPF Sudbury67 domain HUC10 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_sudbury67_huc10_polygon.recharge_adjustment',
    units       => '-',
    description => 'HSPF Sudbury67 domain HUC10 polygon recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_sudbury67_huc12_point.hydrology',
    units       => '-',
    description => 'HSPF Sudbury67 domain HUC12 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_sudbury67_huc12_point.recharge_adjustment',
    units       => '-',
    description => 'HSPF Sudbury67 domain HUC12 point recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_sudbury67_huc12_polygon.hydrology',
    units       => '-',
    description => 'HSPF Sudbury67 domain HUC12 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_sudbury67_huc12_polygon.recharge_adjustment',
    units       => '-',
    description => 'HSPF Sudbury67 domain HUC12 polygon recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_sudbury68_huc10_point.hydrology',
    units       => '-',
    description => 'HSPF Sudbury68 domain HUC10 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_sudbury68_huc10_point.recharge_adjustment',
    units       => '-',
    description => 'HSPF Sudbury68 domain HUC10 point recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_sudbury68_huc10_polygon.hydrology',
    units       => '-',
    description => 'HSPF Sudbury68 domain HUC10 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_sudbury68_huc10_polygon.recharge_adjustment',
    units       => '-',
    description => 'HSPF Sudbury68 domain HUC10 polygon recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_sudbury68_huc12_point.hydrology',
    units       => '-',
    description => 'HSPF Sudbury68 domain HUC12 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_sudbury68_huc12_point.recharge_adjustment',
    units       => '-',
    description => 'HSPF Sudbury68 domain HUC12 point recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_sudbury68_huc12_polygon.hydrology',
    units       => '-',
    description => 'HSPF Sudbury68 domain HUC12 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_sudbury68_huc12_polygon.recharge_adjustment',
    units       => '-',
    description => 'HSPF Sudbury68 domain HUC12 polygon recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_sudbury7_huc10_point.hydrology',
    units       => '-',
    description => 'HSPF Sudbury7 domain HUC10 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_sudbury7_huc10_point.recharge_adjustment',
    units       => '-',
    description => 'HSPF Sudbury7 domain HUC10 point recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_sudbury7_huc10_polygon.hydrology',
    units       => '-',
    description => 'HSPF Sudbury7 domain HUC10 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_sudbury7_huc10_polygon.recharge_adjustment',
    units       => '-',
    description => 'HSPF Sudbury7 domain HUC10 polygon recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_sudbury7_huc12_point.hydrology',
    units       => '-',
    description => 'HSPF Sudbury7 domain HUC12 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_sudbury7_huc12_point.recharge_adjustment',
    units       => '-',
    description => 'HSPF Sudbury7 domain HUC12 point recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_sudbury7_huc12_polygon.hydrology',
    units       => '-',
    description => 'HSPF Sudbury7 domain HUC12 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_sudbury7_huc12_polygon.recharge_adjustment',
    units       => '-',
    description => 'HSPF Sudbury7 domain HUC12 polygon recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_taunton16_huc10_point.hydrology',
    units       => '-',
    description => 'HSPF Taunton16 domain HUC10 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_taunton16_huc10_point.recharge_adjustment',
    units       => '-',
    description => 'HSPF Taunton16 domain HUC10 point recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_taunton16_huc10_polygon.hydrology',
    units       => '-',
    description => 'HSPF Taunton16 domain HUC10 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_taunton16_huc10_polygon.recharge_adjustment',
    units       => '-',
    description => 'HSPF Taunton16 domain HUC10 polygon recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_taunton16_huc12_point.hydrology',
    units       => '-',
    description => 'HSPF Taunton16 domain HUC12 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_taunton16_huc12_point.recharge_adjustment',
    units       => '-',
    description => 'HSPF Taunton16 domain HUC12 point recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_taunton16_huc12_polygon.hydrology',
    units       => '-',
    description => 'HSPF Taunton16 domain HUC12 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_taunton16_huc12_polygon.recharge_adjustment',
    units       => '-',
    description => 'HSPF Taunton16 domain HUC12 polygon recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_taunton69_huc10_point.hydrology',
    units       => '-',
    description => 'HSPF Taunton69 domain HUC10 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_taunton69_huc10_point.recharge_adjustment',
    units       => '-',
    description => 'HSPF Taunton69 domain HUC10 point recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_taunton69_huc10_polygon.hydrology',
    units       => '-',
    description => 'HSPF Taunton69 domain HUC10 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_taunton69_huc10_polygon.recharge_adjustment',
    units       => '-',
    description => 'HSPF Taunton69 domain HUC10 polygon recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_taunton69_huc12_point.hydrology',
    units       => '-',
    description => 'HSPF Taunton69 domain HUC12 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_taunton69_huc12_point.recharge_adjustment',
    units       => '-',
    description => 'HSPF Taunton69 domain HUC12 point recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_taunton69_huc12_polygon.hydrology',
    units       => '-',
    description => 'HSPF Taunton69 domain HUC12 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_taunton69_huc12_polygon.recharge_adjustment',
    units       => '-',
    description => 'HSPF Taunton69 domain HUC12 polygon recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_taunton70_huc10_point.hydrology',
    units       => '-',
    description => 'HSPF Taunton70 domain HUC10 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_taunton70_huc10_point.recharge_adjustment',
    units       => '-',
    description => 'HSPF Taunton70 domain HUC10 point recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_taunton70_huc10_polygon.hydrology',
    units       => '-',
    description => 'HSPF Taunton70 domain HUC10 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_taunton70_huc10_polygon.recharge_adjustment',
    units       => '-',
    description => 'HSPF Taunton70 domain HUC10 polygon recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_taunton70_huc12_point.hydrology',
    units       => '-',
    description => 'HSPF Taunton70 domain HUC12 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_taunton70_huc12_point.recharge_adjustment',
    units       => '-',
    description => 'HSPF Taunton70 domain HUC12 point recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_taunton70_huc12_polygon.hydrology',
    units       => '-',
    description => 'HSPF Taunton70 domain HUC12 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_taunton70_huc12_polygon.recharge_adjustment',
    units       => '-',
    description => 'HSPF Taunton70 domain HUC12 polygon recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_taunton71_huc10_point.hydrology',
    units       => '-',
    description => 'HSPF Taunton71 domain HUC10 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_taunton71_huc10_point.recharge_adjustment',
    units       => '-',
    description => 'HSPF Taunton71 domain HUC10 point recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_taunton71_huc10_polygon.hydrology',
    units       => '-',
    description => 'HSPF Taunton71 domain HUC10 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_taunton71_huc10_polygon.recharge_adjustment',
    units       => '-',
    description => 'HSPF Taunton71 domain HUC10 polygon recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_taunton71_huc12_point.hydrology',
    units       => '-',
    description => 'HSPF Taunton71 domain HUC12 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_taunton71_huc12_point.recharge_adjustment',
    units       => '-',
    description => 'HSPF Taunton71 domain HUC12 point recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_taunton71_huc12_polygon.hydrology',
    units       => '-',
    description => 'HSPF Taunton71 domain HUC12 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_taunton71_huc12_polygon.recharge_adjustment',
    units       => '-',
    description => 'HSPF Taunton71 domain HUC12 polygon recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_taunton72_huc10_point.hydrology',
    units       => '-',
    description => 'HSPF Taunton72 domain HUC10 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_taunton72_huc10_point.recharge_adjustment',
    units       => '-',
    description => 'HSPF Taunton72 domain HUC10 point recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_taunton72_huc10_polygon.hydrology',
    units       => '-',
    description => 'HSPF Taunton72 domain HUC10 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_taunton72_huc10_polygon.recharge_adjustment',
    units       => '-',
    description => 'HSPF Taunton72 domain HUC10 polygon recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_taunton72_huc12_point.hydrology',
    units       => '-',
    description => 'HSPF Taunton72 domain HUC12 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_taunton72_huc12_point.recharge_adjustment',
    units       => '-',
    description => 'HSPF Taunton72 domain HUC12 point recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_taunton72_huc12_polygon.hydrology',
    units       => '-',
    description => 'HSPF Taunton72 domain HUC12 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_taunton72_huc12_polygon.recharge_adjustment',
    units       => '-',
    description => 'HSPF Taunton72 domain HUC12 polygon recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_taunton73_huc10_point.hydrology',
    units       => '-',
    description => 'HSPF Taunton73 domain HUC10 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_taunton73_huc10_point.recharge_adjustment',
    units       => '-',
    description => 'HSPF Taunton73 domain HUC10 point recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_taunton73_huc10_polygon.hydrology',
    units       => '-',
    description => 'HSPF Taunton73 domain HUC10 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_taunton73_huc10_polygon.recharge_adjustment',
    units       => '-',
    description => 'HSPF Taunton73 domain HUC10 polygon recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_taunton73_huc12_point.hydrology',
    units       => '-',
    description => 'HSPF Taunton73 domain HUC12 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_taunton73_huc12_point.recharge_adjustment',
    units       => '-',
    description => 'HSPF Taunton73 domain HUC12 point recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_taunton73_huc12_polygon.hydrology',
    units       => '-',
    description => 'HSPF Taunton73 domain HUC12 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_taunton73_huc12_polygon.recharge_adjustment',
    units       => '-',
    description => 'HSPF Taunton73 domain HUC12 polygon recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_taunton8_huc10_point.hydrology',
    units       => '-',
    description => 'HSPF Taunton8 domain HUC10 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_taunton8_huc10_point.recharge_adjustment',
    units       => '-',
    description => 'HSPF Taunton8 domain HUC10 point recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_taunton8_huc10_polygon.hydrology',
    units       => '-',
    description => 'HSPF Taunton8 domain HUC10 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_taunton8_huc10_polygon.recharge_adjustment',
    units       => '-',
    description => 'HSPF Taunton8 domain HUC10 polygon recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_taunton8_huc12_point.hydrology',
    units       => '-',
    description => 'HSPF Taunton8 domain HUC12 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_taunton8_huc12_point.recharge_adjustment',
    units       => '-',
    description => 'HSPF Taunton8 domain HUC12 point recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_taunton8_huc12_polygon.hydrology',
    units       => '-',
    description => 'HSPF Taunton8 domain HUC12 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'hspf_taunton8_huc12_polygon.recharge_adjustment',
    units       => '-',
    description => 'HSPF Taunton8 domain HUC12 polygon recharge adjustment runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'swmm_blackstone17_huc10_point.hydrology',
    units       => '-',
    description => 'SWMM Blackstone17 domain HUC10 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'swmm_blackstone17_huc10_point.loadings_tn',
    units       => '-',
    description => 'SWMM Blackstone17 domain HUC10 point loadings of total nitrogen runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_blackstone17_huc10_point.loadings_tp',
    units       => '-',
    description => 'SWMM Blackstone17 domain HUC10 point loadings of total phosphorus runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_blackstone17_huc10_point.loadings_tss',
    units       => '-',
    description => 'SWMM Blackstone17 domain HUC10 point loadings of total suspended solids runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_blackstone17_huc10_point.loadings_zn',
    units       => '-',
    description => 'SWMM Blackstone17 domain HUC10 point loadings of total zinc runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_blackstone17_huc10_polygon.hydrology',
    units       => '-',
    description => 'SWMM Blackstone17 domain HUC10 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'swmm_blackstone17_huc10_polygon.loadings_tn',
    units       => '-',
    description => 'SWMM Blackstone17 domain HUC10 polygon loadings of total nitrogen runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_blackstone17_huc10_polygon.loadings_tp',
    units       => '-',
    description => 'SWMM Blackstone17 domain HUC10 polygon loadings of total phosphorus runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_blackstone17_huc10_polygon.loadings_tss',
    units       => '-',
    description => 'SWMM Blackstone17 domain HUC10 polygon loadings of total suspended solids runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_blackstone17_huc10_polygon.loadings_zn',
    units       => '-',
    description => 'SWMM Blackstone17 domain HUC10 polygon loadings of total zinc runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_blackstone17_huc12_point.hydrology',
    units       => '-',
    description => 'SWMM Blackstone17 domain HUC12 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'swmm_blackstone17_huc12_point.loadings_tn',
    units       => '-',
    description => 'SWMM Blackstone17 domain HUC12 point loadings of total nitrogen runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_blackstone17_huc12_point.loadings_tp',
    units       => '-',
    description => 'SWMM Blackstone17 domain HUC12 point loadings of total phosphorus runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_blackstone17_huc12_point.loadings_tss',
    units       => '-',
    description => 'SWMM Blackstone17 domain HUC12 point loadings of total suspended solids runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_blackstone17_huc12_point.loadings_zn',
    units       => '-',
    description => 'SWMM Blackstone17 domain HUC12 point loadings of total zinc runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_blackstone17_huc12_polygon.hydrology',
    units       => '-',
    description => 'SWMM Blackstone17 domain HUC12 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'swmm_blackstone17_huc12_polygon.loadings_tn',
    units       => '-',
    description => 'SWMM Blackstone17 domain HUC12 polygon loadings of total nitrogen runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_blackstone17_huc12_polygon.loadings_tp',
    units       => '-',
    description => 'SWMM Blackstone17 domain HUC12 polygon loadings of total phosphorus runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_blackstone17_huc12_polygon.loadings_tss',
    units       => '-',
    description => 'SWMM Blackstone17 domain HUC12 polygon loadings of total suspended solids runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_blackstone17_huc12_polygon.loadings_zn',
    units       => '-',
    description => 'SWMM Blackstone17 domain HUC12 polygon loadings of total zinc runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_blackstone18_huc10_point.hydrology',
    units       => '-',
    description => 'SWMM Blackstone18 domain HUC10 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'swmm_blackstone18_huc10_point.loadings_tn',
    units       => '-',
    description => 'SWMM Blackstone18 domain HUC10 point loadings of total nitrogen runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_blackstone18_huc10_point.loadings_tp',
    units       => '-',
    description => 'SWMM Blackstone18 domain HUC10 point loadings of total phosphorus runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_blackstone18_huc10_point.loadings_tss',
    units       => '-',
    description => 'SWMM Blackstone18 domain HUC10 point loadings of total suspended solids runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_blackstone18_huc10_point.loadings_zn',
    units       => '-',
    description => 'SWMM Blackstone18 domain HUC10 point loadings of total zinc runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_blackstone18_huc10_polygon.hydrology',
    units       => '-',
    description => 'SWMM Blackstone18 domain HUC10 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'swmm_blackstone18_huc10_polygon.loadings_tn',
    units       => '-',
    description => 'SWMM Blackstone18 domain HUC10 polygon loadings of total nitrogen runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_blackstone18_huc10_polygon.loadings_tp',
    units       => '-',
    description => 'SWMM Blackstone18 domain HUC10 polygon loadings of total phosphorus runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_blackstone18_huc10_polygon.loadings_tss',
    units       => '-',
    description => 'SWMM Blackstone18 domain HUC10 polygon loadings of total suspended solids runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_blackstone18_huc10_polygon.loadings_zn',
    units       => '-',
    description => 'SWMM Blackstone18 domain HUC10 polygon loadings of total zinc runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_blackstone18_huc12_point.hydrology',
    units       => '-',
    description => 'SWMM Blackstone18 domain HUC12 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'swmm_blackstone18_huc12_point.loadings_tn',
    units       => '-',
    description => 'SWMM Blackstone18 domain HUC12 point loadings of total nitrogen runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_blackstone18_huc12_point.loadings_tp',
    units       => '-',
    description => 'SWMM Blackstone18 domain HUC12 point loadings of total phosphorus runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_blackstone18_huc12_point.loadings_tss',
    units       => '-',
    description => 'SWMM Blackstone18 domain HUC12 point loadings of total suspended solids runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_blackstone18_huc12_point.loadings_zn',
    units       => '-',
    description => 'SWMM Blackstone18 domain HUC12 point loadings of total zinc runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_blackstone18_huc12_polygon.hydrology',
    units       => '-',
    description => 'SWMM Blackstone18 domain HUC12 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'swmm_blackstone18_huc12_polygon.loadings_tn',
    units       => '-',
    description => 'SWMM Blackstone18 domain HUC12 polygon loadings of total nitrogen runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_blackstone18_huc12_polygon.loadings_tp',
    units       => '-',
    description => 'SWMM Blackstone18 domain HUC12 polygon loadings of total phosphorus runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_blackstone18_huc12_polygon.loadings_tss',
    units       => '-',
    description => 'SWMM Blackstone18 domain HUC12 polygon loadings of total suspended solids runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_blackstone18_huc12_polygon.loadings_zn',
    units       => '-',
    description => 'SWMM Blackstone18 domain HUC12 polygon loadings of total zinc runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_blackstone25_huc10_point.hydrology',
    units       => '-',
    description => 'SWMM Blackstone25 domain HUC10 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'swmm_blackstone25_huc10_point.loadings_tn',
    units       => '-',
    description => 'SWMM Blackstone25 domain HUC10 point loadings of total nitrogen runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_blackstone25_huc10_point.loadings_tp',
    units       => '-',
    description => 'SWMM Blackstone25 domain HUC10 point loadings of total phosphorus runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_blackstone25_huc10_point.loadings_tss',
    units       => '-',
    description => 'SWMM Blackstone25 domain HUC10 point loadings of total suspended solids runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_blackstone25_huc10_point.loadings_zn',
    units       => '-',
    description => 'SWMM Blackstone25 domain HUC10 point loadings of total zinc runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_blackstone25_huc10_polygon.hydrology',
    units       => '-',
    description => 'SWMM Blackstone25 domain HUC10 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'swmm_blackstone25_huc10_polygon.loadings_tn',
    units       => '-',
    description => 'SWMM Blackstone25 domain HUC10 polygon loadings of total nitrogen runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_blackstone25_huc10_polygon.loadings_tp',
    units       => '-',
    description => 'SWMM Blackstone25 domain HUC10 polygon loadings of total phosphorus runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_blackstone25_huc10_polygon.loadings_tss',
    units       => '-',
    description => 'SWMM Blackstone25 domain HUC10 polygon loadings of total suspended solids runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_blackstone25_huc10_polygon.loadings_zn',
    units       => '-',
    description => 'SWMM Blackstone25 domain HUC10 polygon loadings of total zinc runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_blackstone25_huc12_point.hydrology',
    units       => '-',
    description => 'SWMM Blackstone25 domain HUC12 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'swmm_blackstone25_huc12_point.loadings_tn',
    units       => '-',
    description => 'SWMM Blackstone25 domain HUC12 point loadings of total nitrogen runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_blackstone25_huc12_point.loadings_tp',
    units       => '-',
    description => 'SWMM Blackstone25 domain HUC12 point loadings of total phosphorus runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_blackstone25_huc12_point.loadings_tss',
    units       => '-',
    description => 'SWMM Blackstone25 domain HUC12 point loadings of total suspended solids runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_blackstone25_huc12_point.loadings_zn',
    units       => '-',
    description => 'SWMM Blackstone25 domain HUC12 point loadings of total zinc runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_blackstone25_huc12_polygon.hydrology',
    units       => '-',
    description => 'SWMM Blackstone25 domain HUC12 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'swmm_blackstone25_huc12_polygon.loadings_tn',
    units       => '-',
    description => 'SWMM Blackstone25 domain HUC12 polygon loadings of total nitrogen runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_blackstone25_huc12_polygon.loadings_tp',
    units       => '-',
    description => 'SWMM Blackstone25 domain HUC12 polygon loadings of total phosphorus runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_blackstone25_huc12_polygon.loadings_tss',
    units       => '-',
    description => 'SWMM Blackstone25 domain HUC12 polygon loadings of total suspended solids runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_blackstone25_huc12_polygon.loadings_zn',
    units       => '-',
    description => 'SWMM Blackstone25 domain HUC12 polygon loadings of total zinc runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_blackstone26_huc10_point.hydrology',
    units       => '-',
    description => 'SWMM Blackstone26 domain HUC10 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'swmm_blackstone26_huc10_point.loadings_tn',
    units       => '-',
    description => 'SWMM Blackstone26 domain HUC10 point loadings of total nitrogen runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_blackstone26_huc10_point.loadings_tp',
    units       => '-',
    description => 'SWMM Blackstone26 domain HUC10 point loadings of total phosphorus runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_blackstone26_huc10_point.loadings_tss',
    units       => '-',
    description => 'SWMM Blackstone26 domain HUC10 point loadings of total suspended solids runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_blackstone26_huc10_point.loadings_zn',
    units       => '-',
    description => 'SWMM Blackstone26 domain HUC10 point loadings of total zinc runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_blackstone26_huc10_polygon.hydrology',
    units       => '-',
    description => 'SWMM Blackstone26 domain HUC10 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'swmm_blackstone26_huc10_polygon.loadings_tn',
    units       => '-',
    description => 'SWMM Blackstone26 domain HUC10 polygon loadings of total nitrogen runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_blackstone26_huc10_polygon.loadings_tp',
    units       => '-',
    description => 'SWMM Blackstone26 domain HUC10 polygon loadings of total phosphorus runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_blackstone26_huc10_polygon.loadings_tss',
    units       => '-',
    description => 'SWMM Blackstone26 domain HUC10 polygon loadings of total suspended solids runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_blackstone26_huc10_polygon.loadings_zn',
    units       => '-',
    description => 'SWMM Blackstone26 domain HUC10 polygon loadings of total zinc runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_blackstone26_huc12_point.hydrology',
    units       => '-',
    description => 'SWMM Blackstone26 domain HUC12 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'swmm_blackstone26_huc12_point.loadings_tn',
    units       => '-',
    description => 'SWMM Blackstone26 domain HUC12 point loadings of total nitrogen runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_blackstone26_huc12_point.loadings_tp',
    units       => '-',
    description => 'SWMM Blackstone26 domain HUC12 point loadings of total phosphorus runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_blackstone26_huc12_point.loadings_tss',
    units       => '-',
    description => 'SWMM Blackstone26 domain HUC12 point loadings of total suspended solids runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_blackstone26_huc12_point.loadings_zn',
    units       => '-',
    description => 'SWMM Blackstone26 domain HUC12 point loadings of total zinc runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_blackstone26_huc12_polygon.hydrology',
    units       => '-',
    description => 'SWMM Blackstone26 domain HUC12 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'swmm_blackstone26_huc12_polygon.loadings_tn',
    units       => '-',
    description => 'SWMM Blackstone26 domain HUC12 polygon loadings of total nitrogen runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_blackstone26_huc12_polygon.loadings_tp',
    units       => '-',
    description => 'SWMM Blackstone26 domain HUC12 polygon loadings of total phosphorus runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_blackstone26_huc12_polygon.loadings_tss',
    units       => '-',
    description => 'SWMM Blackstone26 domain HUC12 polygon loadings of total suspended solids runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_blackstone26_huc12_polygon.loadings_zn',
    units       => '-',
    description => 'SWMM Blackstone26 domain HUC12 polygon loadings of total zinc runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_charles19_huc10_point.hydrology',
    units       => '-',
    description => 'SWMM Charles19 domain HUC10 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'swmm_charles19_huc10_point.loadings_tn',
    units       => '-',
    description => 'SWMM Charles19 domain HUC10 point loadings of total nitrogen runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_charles19_huc10_point.loadings_tp',
    units       => '-',
    description => 'SWMM Charles19 domain HUC10 point loadings of total phosphorus runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_charles19_huc10_point.loadings_tss',
    units       => '-',
    description => 'SWMM Charles19 domain HUC10 point loadings of total suspended solids runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_charles19_huc10_point.loadings_zn',
    units       => '-',
    description => 'SWMM Charles19 domain HUC10 point loadings of total zinc runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_charles19_huc10_polygon.hydrology',
    units       => '-',
    description => 'SWMM Charles19 domain HUC10 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'swmm_charles19_huc10_polygon.loadings_tn',
    units       => '-',
    description => 'SWMM Charles19 domain HUC10 polygon loadings of total nitrogen runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_charles19_huc10_polygon.loadings_tp',
    units       => '-',
    description => 'SWMM Charles19 domain HUC10 polygon loadings of total phosphorus runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_charles19_huc10_polygon.loadings_tss',
    units       => '-',
    description => 'SWMM Charles19 domain HUC10 polygon loadings of total suspended solids runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_charles19_huc10_polygon.loadings_zn',
    units       => '-',
    description => 'SWMM Charles19 domain HUC10 polygon loadings of total zinc runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_charles19_huc12_point.hydrology',
    units       => '-',
    description => 'SWMM Charles19 domain HUC12 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'swmm_charles19_huc12_point.loadings_tn',
    units       => '-',
    description => 'SWMM Charles19 domain HUC12 point loadings of total nitrogen runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_charles19_huc12_point.loadings_tp',
    units       => '-',
    description => 'SWMM Charles19 domain HUC12 point loadings of total phosphorus runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_charles19_huc12_point.loadings_tss',
    units       => '-',
    description => 'SWMM Charles19 domain HUC12 point loadings of total suspended solids runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_charles19_huc12_point.loadings_zn',
    units       => '-',
    description => 'SWMM Charles19 domain HUC12 point loadings of total zinc runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_charles19_huc12_polygon.hydrology',
    units       => '-',
    description => 'SWMM Charles19 domain HUC12 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'swmm_charles19_huc12_polygon.loadings_tn',
    units       => '-',
    description => 'SWMM Charles19 domain HUC12 polygon loadings of total nitrogen runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_charles19_huc12_polygon.loadings_tp',
    units       => '-',
    description => 'SWMM Charles19 domain HUC12 polygon loadings of total phosphorus runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_charles19_huc12_polygon.loadings_tss',
    units       => '-',
    description => 'SWMM Charles19 domain HUC12 polygon loadings of total suspended solids runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_charles19_huc12_polygon.loadings_zn',
    units       => '-',
    description => 'SWMM Charles19 domain HUC12 polygon loadings of total zinc runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_charles27_huc10_point.hydrology',
    units       => '-',
    description => 'SWMM Charles27 domain HUC10 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'swmm_charles27_huc10_point.loadings_tn',
    units       => '-',
    description => 'SWMM Charles27 domain HUC10 point loadings of total nitrogen runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_charles27_huc10_point.loadings_tp',
    units       => '-',
    description => 'SWMM Charles27 domain HUC10 point loadings of total phosphorus runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_charles27_huc10_point.loadings_tss',
    units       => '-',
    description => 'SWMM Charles27 domain HUC10 point loadings of total suspended solids runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_charles27_huc10_point.loadings_zn',
    units       => '-',
    description => 'SWMM Charles27 domain HUC10 point loadings of total zinc runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_charles27_huc10_polygon.hydrology',
    units       => '-',
    description => 'SWMM Charles27 domain HUC10 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'swmm_charles27_huc10_polygon.loadings_tn',
    units       => '-',
    description => 'SWMM Charles27 domain HUC10 polygon loadings of total nitrogen runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_charles27_huc10_polygon.loadings_tp',
    units       => '-',
    description => 'SWMM Charles27 domain HUC10 polygon loadings of total phosphorus runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_charles27_huc10_polygon.loadings_tss',
    units       => '-',
    description => 'SWMM Charles27 domain HUC10 polygon loadings of total suspended solids runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_charles27_huc10_polygon.loadings_zn',
    units       => '-',
    description => 'SWMM Charles27 domain HUC10 polygon loadings of total zinc runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_charles27_huc12_point.hydrology',
    units       => '-',
    description => 'SWMM Charles27 domain HUC12 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'swmm_charles27_huc12_point.loadings_tn',
    units       => '-',
    description => 'SWMM Charles27 domain HUC12 point loadings of total nitrogen runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_charles27_huc12_point.loadings_tp',
    units       => '-',
    description => 'SWMM Charles27 domain HUC12 point loadings of total phosphorus runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_charles27_huc12_point.loadings_tss',
    units       => '-',
    description => 'SWMM Charles27 domain HUC12 point loadings of total suspended solids runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_charles27_huc12_point.loadings_zn',
    units       => '-',
    description => 'SWMM Charles27 domain HUC12 point loadings of total zinc runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_charles27_huc12_polygon.hydrology',
    units       => '-',
    description => 'SWMM Charles27 domain HUC12 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'swmm_charles27_huc12_polygon.loadings_tn',
    units       => '-',
    description => 'SWMM Charles27 domain HUC12 polygon loadings of total nitrogen runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_charles27_huc12_polygon.loadings_tp',
    units       => '-',
    description => 'SWMM Charles27 domain HUC12 polygon loadings of total phosphorus runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_charles27_huc12_polygon.loadings_tss',
    units       => '-',
    description => 'SWMM Charles27 domain HUC12 polygon loadings of total suspended solids runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_charles27_huc12_polygon.loadings_zn',
    units       => '-',
    description => 'SWMM Charles27 domain HUC12 polygon loadings of total zinc runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_farmington20_huc10_point.hydrology',
    units       => '-',
    description => 'SWMM Farmington20 domain HUC10 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'swmm_farmington20_huc10_point.loadings_tn',
    units       => '-',
    description => 'SWMM Farmington20 domain HUC10 point loadings of total nitrogen runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_farmington20_huc10_point.loadings_tp',
    units       => '-',
    description => 'SWMM Farmington20 domain HUC10 point loadings of total phosphorus runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_farmington20_huc10_point.loadings_tss',
    units       => '-',
    description => 'SWMM Farmington20 domain HUC10 point loadings of total suspended solids runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_farmington20_huc10_point.loadings_zn',
    units       => '-',
    description => 'SWMM Farmington20 domain HUC10 point loadings of total zinc runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_farmington20_huc10_polygon.hydrology',
    units       => '-',
    description => 'SWMM Farmington20 domain HUC10 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'swmm_farmington20_huc10_polygon.loadings_tn',
    units       => '-',
    description => 'SWMM Farmington20 domain HUC10 polygon loadings of total nitrogen runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_farmington20_huc10_polygon.loadings_tp',
    units       => '-',
    description => 'SWMM Farmington20 domain HUC10 polygon loadings of total phosphorus runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_farmington20_huc10_polygon.loadings_tss',
    units       => '-',
    description => 'SWMM Farmington20 domain HUC10 polygon loadings of total suspended solids runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_farmington20_huc10_polygon.loadings_zn',
    units       => '-',
    description => 'SWMM Farmington20 domain HUC10 polygon loadings of total zinc runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_farmington20_huc12_point.hydrology',
    units       => '-',
    description => 'SWMM Farmington20 domain HUC12 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'swmm_farmington20_huc12_point.loadings_tn',
    units       => '-',
    description => 'SWMM Farmington20 domain HUC12 point loadings of total nitrogen runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_farmington20_huc12_point.loadings_tp',
    units       => '-',
    description => 'SWMM Farmington20 domain HUC12 point loadings of total phosphorus runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_farmington20_huc12_point.loadings_tss',
    units       => '-',
    description => 'SWMM Farmington20 domain HUC12 point loadings of total suspended solids runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_farmington20_huc12_point.loadings_zn',
    units       => '-',
    description => 'SWMM Farmington20 domain HUC12 point loadings of total zinc runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_farmington20_huc12_polygon.hydrology',
    units       => '-',
    description => 'SWMM Farmington20 domain HUC12 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'swmm_farmington20_huc12_polygon.loadings_tn',
    units       => '-',
    description => 'SWMM Farmington20 domain HUC12 polygon loadings of total nitrogen runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_farmington20_huc12_polygon.loadings_tp',
    units       => '-',
    description => 'SWMM Farmington20 domain HUC12 polygon loadings of total phosphorus runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_farmington20_huc12_polygon.loadings_tss',
    units       => '-',
    description => 'SWMM Farmington20 domain HUC12 polygon loadings of total suspended solids runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_farmington20_huc12_polygon.loadings_zn',
    units       => '-',
    description => 'SWMM Farmington20 domain HUC12 polygon loadings of total zinc runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_farmington28_huc10_point.hydrology',
    units       => '-',
    description => 'SWMM Farmington28 domain HUC10 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'swmm_farmington28_huc10_point.loadings_tn',
    units       => '-',
    description => 'SWMM Farmington28 domain HUC10 point loadings of total nitrogen runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_farmington28_huc10_point.loadings_tp',
    units       => '-',
    description => 'SWMM Farmington28 domain HUC10 point loadings of total phosphorus runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_farmington28_huc10_point.loadings_tss',
    units       => '-',
    description => 'SWMM Farmington28 domain HUC10 point loadings of total suspended solids runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_farmington28_huc10_point.loadings_zn',
    units       => '-',
    description => 'SWMM Farmington28 domain HUC10 point loadings of total zinc runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_farmington28_huc10_polygon.hydrology',
    units       => '-',
    description => 'SWMM Farmington28 domain HUC10 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'swmm_farmington28_huc10_polygon.loadings_tn',
    units       => '-',
    description => 'SWMM Farmington28 domain HUC10 polygon loadings of total nitrogen runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_farmington28_huc10_polygon.loadings_tp',
    units       => '-',
    description => 'SWMM Farmington28 domain HUC10 polygon loadings of total phosphorus runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_farmington28_huc10_polygon.loadings_tss',
    units       => '-',
    description => 'SWMM Farmington28 domain HUC10 polygon loadings of total suspended solids runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_farmington28_huc10_polygon.loadings_zn',
    units       => '-',
    description => 'SWMM Farmington28 domain HUC10 polygon loadings of total zinc runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_farmington28_huc12_point.hydrology',
    units       => '-',
    description => 'SWMM Farmington28 domain HUC12 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'swmm_farmington28_huc12_point.loadings_tn',
    units       => '-',
    description => 'SWMM Farmington28 domain HUC12 point loadings of total nitrogen runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_farmington28_huc12_point.loadings_tp',
    units       => '-',
    description => 'SWMM Farmington28 domain HUC12 point loadings of total phosphorus runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_farmington28_huc12_point.loadings_tss',
    units       => '-',
    description => 'SWMM Farmington28 domain HUC12 point loadings of total suspended solids runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_farmington28_huc12_point.loadings_zn',
    units       => '-',
    description => 'SWMM Farmington28 domain HUC12 point loadings of total zinc runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_farmington28_huc12_polygon.hydrology',
    units       => '-',
    description => 'SWMM Farmington28 domain HUC12 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'swmm_farmington28_huc12_polygon.loadings_tn',
    units       => '-',
    description => 'SWMM Farmington28 domain HUC12 polygon loadings of total nitrogen runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_farmington28_huc12_polygon.loadings_tp',
    units       => '-',
    description => 'SWMM Farmington28 domain HUC12 polygon loadings of total phosphorus runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_farmington28_huc12_polygon.loadings_tss',
    units       => '-',
    description => 'SWMM Farmington28 domain HUC12 polygon loadings of total suspended solids runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_farmington28_huc12_polygon.loadings_zn',
    units       => '-',
    description => 'SWMM Farmington28 domain HUC12 polygon loadings of total zinc runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_ipswich21_huc10_point.hydrology',
    units       => '-',
    description => 'SWMM Ipswich21 domain HUC10 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'swmm_ipswich21_huc10_point.loadings_tn',
    units       => '-',
    description => 'SWMM Ipswich21 domain HUC10 point loadings of total nitrogen runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_ipswich21_huc10_point.loadings_tp',
    units       => '-',
    description => 'SWMM Ipswich21 domain HUC10 point loadings of total phosphorus runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_ipswich21_huc10_point.loadings_tss',
    units       => '-',
    description => 'SWMM Ipswich21 domain HUC10 point loadings of total suspended solids runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_ipswich21_huc10_point.loadings_zn',
    units       => '-',
    description => 'SWMM Ipswich21 domain HUC10 point loadings of total zinc runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_ipswich21_huc10_polygon.hydrology',
    units       => '-',
    description => 'SWMM Ipswich21 domain HUC10 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'swmm_ipswich21_huc10_polygon.loadings_tn',
    units       => '-',
    description => 'SWMM Ipswich21 domain HUC10 polygon loadings of total nitrogen runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_ipswich21_huc10_polygon.loadings_tp',
    units       => '-',
    description => 'SWMM Ipswich21 domain HUC10 polygon loadings of total phosphorus runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_ipswich21_huc10_polygon.loadings_tss',
    units       => '-',
    description => 'SWMM Ipswich21 domain HUC10 polygon loadings of total suspended solids runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_ipswich21_huc10_polygon.loadings_zn',
    units       => '-',
    description => 'SWMM Ipswich21 domain HUC10 polygon loadings of total zinc runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_ipswich21_huc12_point.hydrology',
    units       => '-',
    description => 'SWMM Ipswich21 domain HUC12 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'swmm_ipswich21_huc12_point.loadings_tn',
    units       => '-',
    description => 'SWMM Ipswich21 domain HUC12 point loadings of total nitrogen runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_ipswich21_huc12_point.loadings_tp',
    units       => '-',
    description => 'SWMM Ipswich21 domain HUC12 point loadings of total phosphorus runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_ipswich21_huc12_point.loadings_tss',
    units       => '-',
    description => 'SWMM Ipswich21 domain HUC12 point loadings of total suspended solids runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_ipswich21_huc12_point.loadings_zn',
    units       => '-',
    description => 'SWMM Ipswich21 domain HUC12 point loadings of total zinc runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_ipswich21_huc12_polygon.hydrology',
    units       => '-',
    description => 'SWMM Ipswich21 domain HUC12 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'swmm_ipswich21_huc12_polygon.loadings_tn',
    units       => '-',
    description => 'SWMM Ipswich21 domain HUC12 polygon loadings of total nitrogen runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_ipswich21_huc12_polygon.loadings_tp',
    units       => '-',
    description => 'SWMM Ipswich21 domain HUC12 polygon loadings of total phosphorus runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_ipswich21_huc12_polygon.loadings_tss',
    units       => '-',
    description => 'SWMM Ipswich21 domain HUC12 polygon loadings of total suspended solids runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_ipswich21_huc12_polygon.loadings_zn',
    units       => '-',
    description => 'SWMM Ipswich21 domain HUC12 polygon loadings of total zinc runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_ipswich29_huc10_point.hydrology',
    units       => '-',
    description => 'SWMM Ipswich29 domain HUC10 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'swmm_ipswich29_huc10_point.loadings_tn',
    units       => '-',
    description => 'SWMM Ipswich29 domain HUC10 point loadings of total nitrogen runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_ipswich29_huc10_point.loadings_tp',
    units       => '-',
    description => 'SWMM Ipswich29 domain HUC10 point loadings of total phosphorus runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_ipswich29_huc10_point.loadings_tss',
    units       => '-',
    description => 'SWMM Ipswich29 domain HUC10 point loadings of total suspended solids runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_ipswich29_huc10_point.loadings_zn',
    units       => '-',
    description => 'SWMM Ipswich29 domain HUC10 point loadings of total zinc runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_ipswich29_huc10_polygon.hydrology',
    units       => '-',
    description => 'SWMM Ipswich29 domain HUC10 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'swmm_ipswich29_huc10_polygon.loadings_tn',
    units       => '-',
    description => 'SWMM Ipswich29 domain HUC10 polygon loadings of total nitrogen runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_ipswich29_huc10_polygon.loadings_tp',
    units       => '-',
    description => 'SWMM Ipswich29 domain HUC10 polygon loadings of total phosphorus runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_ipswich29_huc10_polygon.loadings_tss',
    units       => '-',
    description => 'SWMM Ipswich29 domain HUC10 polygon loadings of total suspended solids runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_ipswich29_huc10_polygon.loadings_zn',
    units       => '-',
    description => 'SWMM Ipswich29 domain HUC10 polygon loadings of total zinc runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_ipswich29_huc12_point.hydrology',
    units       => '-',
    description => 'SWMM Ipswich29 domain HUC12 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'swmm_ipswich29_huc12_point.loadings_tn',
    units       => '-',
    description => 'SWMM Ipswich29 domain HUC12 point loadings of total nitrogen runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_ipswich29_huc12_point.loadings_tp',
    units       => '-',
    description => 'SWMM Ipswich29 domain HUC12 point loadings of total phosphorus runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_ipswich29_huc12_point.loadings_tss',
    units       => '-',
    description => 'SWMM Ipswich29 domain HUC12 point loadings of total suspended solids runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_ipswich29_huc12_point.loadings_zn',
    units       => '-',
    description => 'SWMM Ipswich29 domain HUC12 point loadings of total zinc runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_ipswich29_huc12_polygon.hydrology',
    units       => '-',
    description => 'SWMM Ipswich29 domain HUC12 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'swmm_ipswich29_huc12_polygon.loadings_tn',
    units       => '-',
    description => 'SWMM Ipswich29 domain HUC12 polygon loadings of total nitrogen runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_ipswich29_huc12_polygon.loadings_tp',
    units       => '-',
    description => 'SWMM Ipswich29 domain HUC12 polygon loadings of total phosphorus runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_ipswich29_huc12_polygon.loadings_tss',
    units       => '-',
    description => 'SWMM Ipswich29 domain HUC12 polygon loadings of total suspended solids runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_ipswich29_huc12_polygon.loadings_zn',
    units       => '-',
    description => 'SWMM Ipswich29 domain HUC12 polygon loadings of total zinc runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_pawcatuck22_huc10_point.hydrology',
    units       => '-',
    description => 'SWMM Pawcatuck22 domain HUC10 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'swmm_pawcatuck22_huc10_point.loadings_tn',
    units       => '-',
    description => 'SWMM Pawcatuck22 domain HUC10 point loadings of total nitrogen runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_pawcatuck22_huc10_point.loadings_tp',
    units       => '-',
    description => 'SWMM Pawcatuck22 domain HUC10 point loadings of total phosphorus runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_pawcatuck22_huc10_point.loadings_tss',
    units       => '-',
    description => 'SWMM Pawcatuck22 domain HUC10 point loadings of total suspended solids runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_pawcatuck22_huc10_point.loadings_zn',
    units       => '-',
    description => 'SWMM Pawcatuck22 domain HUC10 point loadings of total zinc runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_pawcatuck22_huc10_polygon.hydrology',
    units       => '-',
    description => 'SWMM Pawcatuck22 domain HUC10 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'swmm_pawcatuck22_huc10_polygon.loadings_tn',
    units       => '-',
    description => 'SWMM Pawcatuck22 domain HUC10 polygon loadings of total nitrogen runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_pawcatuck22_huc10_polygon.loadings_tp',
    units       => '-',
    description => 'SWMM Pawcatuck22 domain HUC10 polygon loadings of total phosphorus runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_pawcatuck22_huc10_polygon.loadings_tss',
    units       => '-',
    description => 'SWMM Pawcatuck22 domain HUC10 polygon loadings of total suspended solids runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_pawcatuck22_huc10_polygon.loadings_zn',
    units       => '-',
    description => 'SWMM Pawcatuck22 domain HUC10 polygon loadings of total zinc runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_pawcatuck22_huc12_point.hydrology',
    units       => '-',
    description => 'SWMM Pawcatuck22 domain HUC12 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'swmm_pawcatuck22_huc12_point.loadings_tn',
    units       => '-',
    description => 'SWMM Pawcatuck22 domain HUC12 point loadings of total nitrogen runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_pawcatuck22_huc12_point.loadings_tp',
    units       => '-',
    description => 'SWMM Pawcatuck22 domain HUC12 point loadings of total phosphorus runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_pawcatuck22_huc12_point.loadings_tss',
    units       => '-',
    description => 'SWMM Pawcatuck22 domain HUC12 point loadings of total suspended solids runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_pawcatuck22_huc12_point.loadings_zn',
    units       => '-',
    description => 'SWMM Pawcatuck22 domain HUC12 point loadings of total zinc runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_pawcatuck22_huc12_polygon.hydrology',
    units       => '-',
    description => 'SWMM Pawcatuck22 domain HUC12 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'swmm_pawcatuck22_huc12_polygon.loadings_tn',
    units       => '-',
    description => 'SWMM Pawcatuck22 domain HUC12 polygon loadings of total nitrogen runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_pawcatuck22_huc12_polygon.loadings_tp',
    units       => '-',
    description => 'SWMM Pawcatuck22 domain HUC12 polygon loadings of total phosphorus runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_pawcatuck22_huc12_polygon.loadings_tss',
    units       => '-',
    description => 'SWMM Pawcatuck22 domain HUC12 polygon loadings of total suspended solids runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_pawcatuck22_huc12_polygon.loadings_zn',
    units       => '-',
    description => 'SWMM Pawcatuck22 domain HUC12 polygon loadings of total zinc runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_pawcatuck30_huc10_point.hydrology',
    units       => '-',
    description => 'SWMM Pawcatuck30 domain HUC10 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'swmm_pawcatuck30_huc10_point.loadings_tn',
    units       => '-',
    description => 'SWMM Pawcatuck30 domain HUC10 point loadings of total nitrogen runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_pawcatuck30_huc10_point.loadings_tp',
    units       => '-',
    description => 'SWMM Pawcatuck30 domain HUC10 point loadings of total phosphorus runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_pawcatuck30_huc10_point.loadings_tss',
    units       => '-',
    description => 'SWMM Pawcatuck30 domain HUC10 point loadings of total suspended solids runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_pawcatuck30_huc10_point.loadings_zn',
    units       => '-',
    description => 'SWMM Pawcatuck30 domain HUC10 point loadings of total zinc runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_pawcatuck30_huc10_polygon.hydrology',
    units       => '-',
    description => 'SWMM Pawcatuck30 domain HUC10 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'swmm_pawcatuck30_huc10_polygon.loadings_tn',
    units       => '-',
    description => 'SWMM Pawcatuck30 domain HUC10 polygon loadings of total nitrogen runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_pawcatuck30_huc10_polygon.loadings_tp',
    units       => '-',
    description => 'SWMM Pawcatuck30 domain HUC10 polygon loadings of total phosphorus runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_pawcatuck30_huc10_polygon.loadings_tss',
    units       => '-',
    description => 'SWMM Pawcatuck30 domain HUC10 polygon loadings of total suspended solids runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_pawcatuck30_huc10_polygon.loadings_zn',
    units       => '-',
    description => 'SWMM Pawcatuck30 domain HUC10 polygon loadings of total zinc runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_pawcatuck30_huc12_point.hydrology',
    units       => '-',
    description => 'SWMM Pawcatuck30 domain HUC12 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'swmm_pawcatuck30_huc12_point.loadings_tn',
    units       => '-',
    description => 'SWMM Pawcatuck30 domain HUC12 point loadings of total nitrogen runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_pawcatuck30_huc12_point.loadings_tp',
    units       => '-',
    description => 'SWMM Pawcatuck30 domain HUC12 point loadings of total phosphorus runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_pawcatuck30_huc12_point.loadings_tss',
    units       => '-',
    description => 'SWMM Pawcatuck30 domain HUC12 point loadings of total suspended solids runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_pawcatuck30_huc12_point.loadings_zn',
    units       => '-',
    description => 'SWMM Pawcatuck30 domain HUC12 point loadings of total zinc runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_pawcatuck30_huc12_polygon.hydrology',
    units       => '-',
    description => 'SWMM Pawcatuck30 domain HUC12 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'swmm_pawcatuck30_huc12_polygon.loadings_tn',
    units       => '-',
    description => 'SWMM Pawcatuck30 domain HUC12 polygon loadings of total nitrogen runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_pawcatuck30_huc12_polygon.loadings_tp',
    units       => '-',
    description => 'SWMM Pawcatuck30 domain HUC12 polygon loadings of total phosphorus runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_pawcatuck30_huc12_polygon.loadings_tss',
    units       => '-',
    description => 'SWMM Pawcatuck30 domain HUC12 polygon loadings of total suspended solids runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_pawcatuck30_huc12_polygon.loadings_zn',
    units       => '-',
    description => 'SWMM Pawcatuck30 domain HUC12 polygon loadings of total zinc runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_sudbury23_huc10_point.hydrology',
    units       => '-',
    description => 'SWMM Sudbury23 domain HUC10 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'swmm_sudbury23_huc10_point.loadings_tn',
    units       => '-',
    description => 'SWMM Sudbury23 domain HUC10 point loadings of total nitrogen runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_sudbury23_huc10_point.loadings_tp',
    units       => '-',
    description => 'SWMM Sudbury23 domain HUC10 point loadings of total phosphorus runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_sudbury23_huc10_point.loadings_tss',
    units       => '-',
    description => 'SWMM Sudbury23 domain HUC10 point loadings of total suspended solids runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_sudbury23_huc10_point.loadings_zn',
    units       => '-',
    description => 'SWMM Sudbury23 domain HUC10 point loadings of total zinc runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_sudbury23_huc10_polygon.hydrology',
    units       => '-',
    description => 'SWMM Sudbury23 domain HUC10 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'swmm_sudbury23_huc10_polygon.loadings_tn',
    units       => '-',
    description => 'SWMM Sudbury23 domain HUC10 polygon loadings of total nitrogen runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_sudbury23_huc10_polygon.loadings_tp',
    units       => '-',
    description => 'SWMM Sudbury23 domain HUC10 polygon loadings of total phosphorus runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_sudbury23_huc10_polygon.loadings_tss',
    units       => '-',
    description => 'SWMM Sudbury23 domain HUC10 polygon loadings of total suspended solids runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_sudbury23_huc10_polygon.loadings_zn',
    units       => '-',
    description => 'SWMM Sudbury23 domain HUC10 polygon loadings of total zinc runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_sudbury23_huc12_point.hydrology',
    units       => '-',
    description => 'SWMM Sudbury23 domain HUC12 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'swmm_sudbury23_huc12_point.loadings_tn',
    units       => '-',
    description => 'SWMM Sudbury23 domain HUC12 point loadings of total nitrogen runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_sudbury23_huc12_point.loadings_tp',
    units       => '-',
    description => 'SWMM Sudbury23 domain HUC12 point loadings of total phosphorus runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_sudbury23_huc12_point.loadings_tss',
    units       => '-',
    description => 'SWMM Sudbury23 domain HUC12 point loadings of total suspended solids runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_sudbury23_huc12_point.loadings_zn',
    units       => '-',
    description => 'SWMM Sudbury23 domain HUC12 point loadings of total zinc runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_sudbury23_huc12_polygon.hydrology',
    units       => '-',
    description => 'SWMM Sudbury23 domain HUC12 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'swmm_sudbury23_huc12_polygon.loadings_tn',
    units       => '-',
    description => 'SWMM Sudbury23 domain HUC12 polygon loadings of total nitrogen runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_sudbury23_huc12_polygon.loadings_tp',
    units       => '-',
    description => 'SWMM Sudbury23 domain HUC12 polygon loadings of total phosphorus runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_sudbury23_huc12_polygon.loadings_tss',
    units       => '-',
    description => 'SWMM Sudbury23 domain HUC12 polygon loadings of total suspended solids runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_sudbury23_huc12_polygon.loadings_zn',
    units       => '-',
    description => 'SWMM Sudbury23 domain HUC12 polygon loadings of total zinc runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_sudbury31_huc10_point.hydrology',
    units       => '-',
    description => 'SWMM Sudbury31 domain HUC10 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'swmm_sudbury31_huc10_point.loadings_tn',
    units       => '-',
    description => 'SWMM Sudbury31 domain HUC10 point loadings of total nitrogen runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_sudbury31_huc10_point.loadings_tp',
    units       => '-',
    description => 'SWMM Sudbury31 domain HUC10 point loadings of total phosphorus runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_sudbury31_huc10_point.loadings_tss',
    units       => '-',
    description => 'SWMM Sudbury31 domain HUC10 point loadings of total suspended solids runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_sudbury31_huc10_point.loadings_zn',
    units       => '-',
    description => 'SWMM Sudbury31 domain HUC10 point loadings of total zinc runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_sudbury31_huc10_polygon.hydrology',
    units       => '-',
    description => 'SWMM Sudbury31 domain HUC10 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'swmm_sudbury31_huc10_polygon.loadings_tn',
    units       => '-',
    description => 'SWMM Sudbury31 domain HUC10 polygon loadings of total nitrogen runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_sudbury31_huc10_polygon.loadings_tp',
    units       => '-',
    description => 'SWMM Sudbury31 domain HUC10 polygon loadings of total phosphorus runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_sudbury31_huc10_polygon.loadings_tss',
    units       => '-',
    description => 'SWMM Sudbury31 domain HUC10 polygon loadings of total suspended solids runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_sudbury31_huc10_polygon.loadings_zn',
    units       => '-',
    description => 'SWMM Sudbury31 domain HUC10 polygon loadings of total zinc runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_sudbury31_huc12_point.hydrology',
    units       => '-',
    description => 'SWMM Sudbury31 domain HUC12 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'swmm_sudbury31_huc12_point.loadings_tn',
    units       => '-',
    description => 'SWMM Sudbury31 domain HUC12 point loadings of total nitrogen runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_sudbury31_huc12_point.loadings_tp',
    units       => '-',
    description => 'SWMM Sudbury31 domain HUC12 point loadings of total phosphorus runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_sudbury31_huc12_point.loadings_tss',
    units       => '-',
    description => 'SWMM Sudbury31 domain HUC12 point loadings of total suspended solids runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_sudbury31_huc12_point.loadings_zn',
    units       => '-',
    description => 'SWMM Sudbury31 domain HUC12 point loadings of total zinc runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_sudbury31_huc12_polygon.hydrology',
    units       => '-',
    description => 'SWMM Sudbury31 domain HUC12 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'swmm_sudbury31_huc12_polygon.loadings_tn',
    units       => '-',
    description => 'SWMM Sudbury31 domain HUC12 polygon loadings of total nitrogen runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_sudbury31_huc12_polygon.loadings_tp',
    units       => '-',
    description => 'SWMM Sudbury31 domain HUC12 polygon loadings of total phosphorus runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_sudbury31_huc12_polygon.loadings_tss',
    units       => '-',
    description => 'SWMM Sudbury31 domain HUC12 polygon loadings of total suspended solids runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_sudbury31_huc12_polygon.loadings_zn',
    units       => '-',
    description => 'SWMM Sudbury31 domain HUC12 polygon loadings of total zinc runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_taunton24_huc10_point.hydrology',
    units       => '-',
    description => 'SWMM Taunton24 domain HUC10 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'swmm_taunton24_huc10_point.loadings_tn',
    units       => '-',
    description => 'SWMM Taunton24 domain HUC10 point loadings of total nitrogen runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_taunton24_huc10_point.loadings_tp',
    units       => '-',
    description => 'SWMM Taunton24 domain HUC10 point loadings of total phosphorus runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_taunton24_huc10_point.loadings_tss',
    units       => '-',
    description => 'SWMM Taunton24 domain HUC10 point loadings of total suspended solids runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_taunton24_huc10_point.loadings_zn',
    units       => '-',
    description => 'SWMM Taunton24 domain HUC10 point loadings of total zinc runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_taunton24_huc10_polygon.hydrology',
    units       => '-',
    description => 'SWMM Taunton24 domain HUC10 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'swmm_taunton24_huc10_polygon.loadings_tn',
    units       => '-',
    description => 'SWMM Taunton24 domain HUC10 polygon loadings of total nitrogen runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_taunton24_huc10_polygon.loadings_tp',
    units       => '-',
    description => 'SWMM Taunton24 domain HUC10 polygon loadings of total phosphorus runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_taunton24_huc10_polygon.loadings_tss',
    units       => '-',
    description => 'SWMM Taunton24 domain HUC10 polygon loadings of total suspended solids runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_taunton24_huc10_polygon.loadings_zn',
    units       => '-',
    description => 'SWMM Taunton24 domain HUC10 polygon loadings of total zinc runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_taunton24_huc12_point.hydrology',
    units       => '-',
    description => 'SWMM Taunton24 domain HUC12 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'swmm_taunton24_huc12_point.loadings_tn',
    units       => '-',
    description => 'SWMM Taunton24 domain HUC12 point loadings of total nitrogen runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_taunton24_huc12_point.loadings_tp',
    units       => '-',
    description => 'SWMM Taunton24 domain HUC12 point loadings of total phosphorus runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_taunton24_huc12_point.loadings_tss',
    units       => '-',
    description => 'SWMM Taunton24 domain HUC12 point loadings of total suspended solids runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_taunton24_huc12_point.loadings_zn',
    units       => '-',
    description => 'SWMM Taunton24 domain HUC12 point loadings of total zinc runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_taunton24_huc12_polygon.hydrology',
    units       => '-',
    description => 'SWMM Taunton24 domain HUC12 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'swmm_taunton24_huc12_polygon.loadings_tn',
    units       => '-',
    description => 'SWMM Taunton24 domain HUC12 polygon loadings of total nitrogen runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_taunton24_huc12_polygon.loadings_tp',
    units       => '-',
    description => 'SWMM Taunton24 domain HUC12 polygon loadings of total phosphorus runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_taunton24_huc12_polygon.loadings_tss',
    units       => '-',
    description => 'SWMM Taunton24 domain HUC12 polygon loadings of total suspended solids runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_taunton24_huc12_polygon.loadings_zn',
    units       => '-',
    description => 'SWMM Taunton24 domain HUC12 polygon loadings of total zinc runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_taunton32_huc10_point.hydrology',
    units       => '-',
    description => 'SWMM Taunton32 domain HUC10 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'swmm_taunton32_huc10_point.loadings_tn',
    units       => '-',
    description => 'SWMM Taunton32 domain HUC10 point loadings of total nitrogen runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_taunton32_huc10_point.loadings_tp',
    units       => '-',
    description => 'SWMM Taunton32 domain HUC10 point loadings of total phosphorus runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_taunton32_huc10_point.loadings_tss',
    units       => '-',
    description => 'SWMM Taunton32 domain HUC10 point loadings of total suspended solids runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_taunton32_huc10_point.loadings_zn',
    units       => '-',
    description => 'SWMM Taunton32 domain HUC10 point loadings of total zinc runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_taunton32_huc10_polygon.hydrology',
    units       => '-',
    description => 'SWMM Taunton32 domain HUC10 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'swmm_taunton32_huc10_polygon.loadings_tn',
    units       => '-',
    description => 'SWMM Taunton32 domain HUC10 polygon loadings of total nitrogen runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_taunton32_huc10_polygon.loadings_tp',
    units       => '-',
    description => 'SWMM Taunton32 domain HUC10 polygon loadings of total phosphorus runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_taunton32_huc10_polygon.loadings_tss',
    units       => '-',
    description => 'SWMM Taunton32 domain HUC10 polygon loadings of total suspended solids runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_taunton32_huc10_polygon.loadings_zn',
    units       => '-',
    description => 'SWMM Taunton32 domain HUC10 polygon loadings of total zinc runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_taunton32_huc12_point.hydrology',
    units       => '-',
    description => 'SWMM Taunton32 domain HUC12 point hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'swmm_taunton32_huc12_point.loadings_tn',
    units       => '-',
    description => 'SWMM Taunton32 domain HUC12 point loadings of total nitrogen runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_taunton32_huc12_point.loadings_tp',
    units       => '-',
    description => 'SWMM Taunton32 domain HUC12 point loadings of total phosphorus runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_taunton32_huc12_point.loadings_tss',
    units       => '-',
    description => 'SWMM Taunton32 domain HUC12 point loadings of total suspended solids runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_taunton32_huc12_point.loadings_zn',
    units       => '-',
    description => 'SWMM Taunton32 domain HUC12 point loadings of total zinc runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_taunton32_huc12_polygon.hydrology',
    units       => '-',
    description => 'SWMM Taunton32 domain HUC12 polygon hydrology runoff (in/acre/hour), recharge (in/acre/hour), precipitation (in/hour), temperature (degF) and potential evapotranspiration (in/hour).'
  },
  {
    name        => 'swmm_taunton32_huc12_polygon.loadings_tn',
    units       => '-',
    description => 'SWMM Taunton32 domain HUC12 polygon loadings of total nitrogen runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_taunton32_huc12_polygon.loadings_tp',
    units       => '-',
    description => 'SWMM Taunton32 domain HUC12 polygon loadings of total phosphorus runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_taunton32_huc12_polygon.loadings_tss',
    units       => '-',
    description => 'SWMM Taunton32 domain HUC12 polygon loadings of total suspended solids runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'swmm_taunton32_huc12_polygon.loadings_zn',
    units       => '-',
    description => 'SWMM Taunton32 domain HUC12 polygon loadings of total zinc runoff (in/acre/hour), recharge (in/acre/hour).'
  },
  {
    name        => 'greenspace_housing',
    units       => '-',
    description => 'Modeled greenspace area housing development statistics'
  },

);



my %region_metadata = (

  'atlantic' => {
    'states' => {
      begin_date  => '2018-04-01',
      end_date    => '2018-04-01',
      west_bound  =>  -87.63490,
      east_bound  =>  -66.88544,
      south_bound =>   24.39631,
      north_bound =>   47.45985
    },
    'provinces' => {
      begin_date  => '2018-04-01',
      end_date    => '2018-04-01',
      west_bound  =>  -79.76106,
      east_bound  =>  -57.10549,
      south_bound =>   44.57031,
      north_bound =>   62.58286
    },
    'counties' => {
      begin_date  => '2018-04-01',
      end_date    => '2018-04-01',
      west_bound  =>  -87.63490,
      east_bound  =>  -66.88544,
      south_bound =>   24.39631,
      north_bound =>   47.45985
    },
    'cities' => {
      begin_date  => '2018-04-01',
      end_date    => '2018-04-01',
      west_bound  =>  -87.63300,
      east_bound  =>  -67.23016,
      south_bound =>   24.54409,
      north_bound =>   46.88095
    },
    'roads' => {
      begin_date  => '2018-04-01',
      end_date    => '2018-04-01',
      west_bound  =>  -87.41240,
      east_bound  =>  -67.78126,
      south_bound =>   25.45512,
      north_bound =>   46.14542
    },
    'tributaries' => {
      begin_date  => '2018-04-01',
      end_date    => '2018-04-01',
      west_bound  =>  -87.44991,
      east_bound  =>  -66.98840,
      south_bound =>   25.10267,
      north_bound =>   45.93816
    },
    'hucs' => {
      begin_date  => '2018-04-01',
      end_date    => '2018-04-01',
      west_bound  =>  -84.44240,
      east_bound  =>  -66.01875,
      south_bound =>   24.39533,
      north_bound =>   48.09970
    },
    'watersheds' => {
      begin_date  => '2018-04-01',
      end_date    => '2018-04-01',
      west_bound  =>  -84.44409,
      east_bound  =>  -66.94960,
      south_bound =>   24.95160,
      north_bound =>   46.46890
    },
    'estuaries' => {
      begin_date  => '2018-04-01',
      end_date    => '2018-04-01',
      west_bound  =>  -87.51210,
      east_bound  =>  -66.93271,
      south_bound =>   24.55000,
      north_bound =>   45.50007
    },
    'lakes' => {
      begin_date  => '2018-04-01',
      end_date    => '2018-04-01',
      west_bound  =>  -85.75558,
      east_bound  =>  -67.12168,
      south_bound =>   25.13734,
      north_bound =>   47.27269
    },
    'legislative_districts' => {
      begin_date  => '2018-04-01',
      end_date    => '2018-04-01',
      west_bound  =>  -87.63490,
      east_bound  =>  -66.88544,
      south_bound =>   24.39631,
      north_bound =>   47.45985
    },

    'nlcd_2011' => {
      begin_date  => '2011-01-01',
      end_date    => '2011-12-31',
      west_bound  => -86,
      east_bound  => -65,
      south_bound => 23,
      north_bound => 50
    },
    'nlcd_2006' => {
      begin_date  => '2006-01-01',
      end_date    => '2006-12-31',
      west_bound  => -86,
      east_bound  => -65,
      south_bound => 23,
      north_bound => 50
    },
    'nlcd_2001' => {
      begin_date  => '2001-01-01',
      end_date    => '2001-12-31',
      west_bound  => -86,
      east_bound  => -65,
      south_bound => 23,
      north_bound => 50
    },
    'nlcd_1992' => {
      begin_date  => '1992-01-01',
      end_date    => '1992-12-31',
      west_bound  => -86,
      east_bound  => -65,
      south_bound => 23,
      north_bound => 50
    },
    'nlcd_change_1992_2001' => {
      begin_date  => '1992-01-01',
      end_date    => '2001-12-31',
      west_bound  => -86,
      east_bound  => -65,
      south_bound => 23,
      north_bound => 50
    },
    'nlcd_change_2001_2006' => {
      begin_date  => '2001-01-01',
      end_date    => '2006-12-31',
      west_bound  => -86,
      east_bound  => -65,
      south_bound => 23,
      north_bound => 50
    },
    'ccap_2006' => {
      begin_date  => '2006-01-01',
      end_date    => '2006-12-31',
      west_bound  => -89,
      east_bound  => -65,
      south_bound => 22,
      north_bound => 50
    },
    'ccap_2001' => {
      begin_date  => '2001-01-01',
      end_date    => '2001-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'decadal_iclus_land_use_2000_2100' => {
      begin_date  => '2000-01-01',
      end_date    => '2100-12-31',
      west_bound  => -86,
      east_bound  => -65,
      south_bound => 23,
      north_bound => 50
    },
    'iclus_land_use_watershed' => {
      begin_date  => '2000-01-01',
      end_date    => '2100-12-31',
      west_bound  => -86,
      east_bound  => -65,
      south_bound => 23,
      north_bound => 50
    },
    'nlcd_2011_land_use_watershed' => {
      begin_date  => '2011-01-01',
      end_date    => '2011-12-31',
      west_bound  => -86,
      east_bound  => -65,
      south_bound => 23,
      north_bound => 50
    },
    'nlcd_2006_land_use_watershed' => {
      begin_date  => '2006-01-01',
      end_date    => '2006-12-31',
      west_bound  => -86,
      east_bound  => -65,
      south_bound => 23,
      north_bound => 50
    },
    'impervious_2006' => {
      begin_date  => '2006-01-01',
      end_date    => '2006-12-31',
      west_bound  => -86,
      east_bound  => -65,
      south_bound => 23,
      north_bound => 50
    },
    'impervious_2001' => {
      begin_date  => '2001-01-01',
      end_date    => '2001-12-31',
      west_bound  => -86,
      east_bound  => -65,
      south_bound => 23,
      north_bound => 50
    },
    'impervious_change_2001_2006' => {
      begin_date  => '2001-01-01',
      end_date    => '2006-12-31',
      west_bound  => -86,
      east_bound  => -65,
      south_bound => 23,
      north_bound => 50
    },
    'decadal_iclus_a1_impervious_2010_2100' => {
      begin_date  => '2010-01-01',
      end_date    => '2100-12-31',
      west_bound  => -86,
      east_bound  => -65,
      south_bound => 23,
      north_bound => 50
    },
    'decadal_iclus_a2_impervious_2010_2100' => {
      begin_date  => '2010-01-01',
      end_date    => '2100-12-31',
      west_bound  => -86,
      east_bound  => -65,
      south_bound => 23,
      north_bound => 50
    },
    'decadal_iclus_b1_impervious_2010_2100' => {
      begin_date  => '2010-01-01',
      end_date    => '2100-12-31',
      west_bound  => -86,
      east_bound  => -65,
      south_bound => 23,
      north_bound => 50
    },
    'decadal_iclus_b2_impervious_2010_2100' => {
      begin_date  => '2010-01-01',
      end_date    => '2100-12-31',
      west_bound  => -86,
      east_bound  => -65,
      south_bound => 23,
      north_bound => 50
    },
    'decadal_iclus_bc_impervious_2010_2100' => {
      begin_date  => '2010-01-01',
      end_date    => '2100-12-31',
      west_bound  => -86,
      east_bound  => -65,
      south_bound => 23,
      north_bound => 50
    },

    'elevation' => {
      begin_date  => '2004-01-01',
      end_date    => '2004-12-31',
      west_bound  => -180,
      east_bound  => -66,
      south_bound => 18,
      north_bound => 66
    },

    # From: https://www.ngdc.noaa.gov/dem/squareCellGrid/search

    'elevation_albemarle_sound_nc_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -77.282,
      east_bound  => -75.306,
      south_bound => 35.346,
      north_bound => 36.779
    },
    'elevation_alligator_cove_fl_navd88' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -81.246,
      east_bound  => -81.225,
      south_bound => 25.598,
      north_bound => 25.612
    },
    'elevation_altahama_sound_ga_navd88' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -81.509,
      east_bound  => -81.258,
      south_bound => 31.13,
      north_bound => 31.384
    },
    'elevation_annisquam_river_ma_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -70.75,
      east_bound  => -70.625,
      south_bound => 42.589,
      north_bound => 42.695
    },
    'elevation_apalachicola_bay_fl_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -85.31,
      east_bound  => -84.378,
      south_bound => 29.51,
      north_bound => 30.034
    },
    'elevation_apalachee_bay_fl_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -84.664,
      east_bound  => -83.704,
      south_bound => 29.802,
      north_bound => 30.385
    },
    'elevation_atlantic_city_nj_mhw' => {
      begin_date  => '2007-10-10',
      end_date    => '2007-10-10',
      west_bound  => -75.05,
      east_bound  => -74,
      south_bound => 38.85,
      north_bound => 39.75
    },
    'elevation_bar_harbor_me_mhw' => {
      begin_date  => '2011-10-17',
      end_date    => '2011-10-17',
      west_bound  => -68.7,
      east_bound  => -67.62,
      south_bound => 43.69,
      north_bound => 44.69
    },
    'elevation_bar_harbor_me_navd88' => {
      begin_date  => '2011-10-17',
      end_date    => '2011-10-17',
      west_bound  => -68.7,
      east_bound  => -67.62,
      south_bound => 43.69,
      north_bound => 44.69
    },
    'elevation_barnegat_bay_nj_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -74.509,
      east_bound  => -73.877,
      south_bound => 39.461,
      north_bound => 40.134
    },
    'elevation_barnstable_bay_ma_navd88' => {
      begin_date  => '2020-06-17',
      end_date    => '2020-06-17',
      west_bound  => -70.444,
      east_bound  => -70.183,
      south_bound => 41.659,
      north_bound => 41.784
    },
    'elevation_big_spring_creek_fl_mllw' => {
      begin_date  => '2020-06-19',
      end_date    => '2020-06-19',
      west_bound  => -83.754,
      east_bound  => -83.728,
      south_bound => 29.957,
      north_bound => 29.976
    },
    'elevation_biscayne_bay_fl_navd88' => {
      begin_date  => '2020-06-19',
      end_date    => '2020-06-19',
      west_bound  => -80.542,
      east_bound  => -79.981,
      south_bound => 25.225,
      north_bound => 26.074
    },
    'elevation_blackwater_sound_fl_navd88' => {
      begin_date  => '2020-06-19',
      end_date    => '2020-06-19',
      west_bound  => -80.425,
      east_bound  => -80.356,
      south_bound => 25.089,
      north_bound => 25.174
    },
    'elevation_blue_hill_bay_me_navd88' => {
      begin_date  => '2020-06-19',
      end_date    => '2020-06-19',
      west_bound  => -68.746,
      east_bound  => -68.167,
      south_bound => 44.082,
      north_bound => 44.595
    },
    'elevation_boggy_bay_fl_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -82.748,
      east_bound  => -82.71,
      south_bound => 28.271,
      north_bound => 28.306
    },
    'elevation_boggy_creek_fl_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -83.337,
      east_bound  => -83.267,
      south_bound => 29.446,
      north_bound => 29.514
    },
    'elevation_bogue_sound_navd88' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -77.375,
      east_bound  => -76.04,
      south_bound => 34.394,
      north_bound => 35.196
    },
    'elevation_boston_harbor_ma_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -71.179,
      east_bound  => -70.761,
      south_bound => 42.168,
      north_bound => 42.472
    },
    'elevation_brave_boat_harbor_me_navd88' => {
      begin_date  => '2020-06-22',
      end_date    => '2020-06-22',
      west_bound  => -70.719,
      east_bound  => -70.611,
      south_bound => 43.053,
      north_bound => 43.14
    },
    'elevation_broad_river_rodgers_river_fl_navd88' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -81.206,
      east_bound  => -80.996,
      south_bound => 25.423,
      north_bound => 25.547
    },
    'elevation_broad_river_sc_navd88' => {
      begin_date  => '2020-06-22',
      end_date    => '2020-06-22',
      west_bound  => -81.012,
      east_bound  => -80.443,
      south_bound => 32.197,
      north_bound => 32.666
    },
    'elevation_broad_sound_ma_mllw' => {
      begin_date  => '2020-06-22',
      end_date    => '2020-06-22',
      west_bound  => -71.075,
      east_bound  => -70.876,
      south_bound => 42.36,
      north_bound => 42.489
    },
    'elevation_buck_creek_fl_mllw' => {
      begin_date  => '2020-06-19',
      end_date    => '2020-06-19',
      west_bound  => -83.402,
      east_bound  => -83.388,
      south_bound => 29.542,
      north_bound => 29.555
    },
    'elevation_bulls_bay_sc_navd88' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -79.647,
      east_bound  => -79.406,
      south_bound => 32.877,
      north_bound => 33.123
    },
    'elevation_bumblebee_creek_fl_mllw' => {
      begin_date  => '2020-06-19',
      end_date    => '2020-06-19',
      west_bound  => -83.197,
      east_bound  => -83.147,
      south_bound => 29.342,
      north_bound => 29.388
    },
    'elevation_bush_river_md_navd88' => {
      begin_date  => '2020-06-19',
      end_date    => '2020-06-19',
      west_bound  => -76.339,
      east_bound  => -76.164,
      south_bound => 39.321,
      north_bound => 39.495
    },
    'elevation_buzzards_bay_ma_mllw' => {
      begin_date  => '2020-06-19',
      end_date    => '2020-06-19',
      west_bound  => -71.254,
      east_bound  => -70.368,
      south_bound => 41.292,
      north_bound => 41.897
    },
    'elevation_cape_fear_river_nc_navd88' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -78.097,
      east_bound  => -77.854,
      south_bound => 33.848,
      north_bound => 34.274
    },
    'elevation_cape_hatteras_nc_mhw' => {
      begin_date  => '2006-07-06',
      end_date    => '2006-07-06',
      west_bound  => -76.05,
      east_bound  => -75.05,
      south_bound => 34.75,
      north_bound => 35.8
    },
    'elevation_cape_romaine_refuge_sc_navd88' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -79.468,
      east_bound  => -79.293,
      south_bound => 32.993,
      north_bound => 33.142
    },
    'elevation_capers_creek_sc_navd88' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -79.739,
      east_bound  => -78.627,
      south_bound => 32.825,
      north_bound => 32.94
    },
    'elevation_casco_bay_me_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -70.462,
      east_bound  => -69.688,
      south_bound => 43.465,
      north_bound => 44.041
    },
    'elevation_central_fl_mhw' => {
      begin_date  => '2014-11-03',
      end_date    => '2014-11-03',
      west_bound  => -80.88,
      east_bound  => -79.95,
      south_bound => 27.27,
      north_bound => 28.88
    },
    'elevation_central_fl_navd88' => {
      begin_date  => '2014-11-03',
      end_date    => '2014-11-03',
      west_bound  => -80.88,
      east_bound  => -79.95,
      south_bound => 27.27,
      north_bound => 28.88
    },
    'elevation_charleston_harbor_sc_navd88' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -80.099,
      east_bound  => -79.708,
      south_bound => 32.671,
      north_bound => 33.009
    },
    'elevation_charlotte_harbor_fl_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -82.5112,
      east_bound  => -81.643,
      south_bound => 26.345,
      north_bound => 27.157
    },
    'elevation_chesapeake_bay_mainstem_1_md_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -76.732,
      east_bound  => -75.467,
      south_bound => 36.785,
      north_bound => 38.06
    },
    'elevation_chesapeake_bay_mainstem_1_va_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -76.732,
      east_bound  => -75.467,
      south_bound => 36.785,
      north_bound => 38.06
    },
    'elevation_chesapeake_bay_mainstem_2_md_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -76.987,
      east_bound  => -75.552,
      south_bound => 37.768,
      north_bound => 39.576
    },
    'elevation_chesapeake_bay_mainstem_2_va_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -76.987,
      east_bound  => -75.552,
      south_bound => 37.768,
      north_bound => 39.576
    },
    'elevation_chester_river_md_navd88' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -76.309,
      east_bound  => -75.834,
      south_bound => 38.918,
      north_bound => 39.321
    },
    'elevation_chincoteague_bay_md_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -74.94,
      east_bound  => -74.936,
      south_bound => 37.775,
      north_bound => 38.436
    },
    'elevation_chincoteague_bay_va_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -74.94,
      east_bound  => -74.936,
      south_bound => 37.775,
      north_bound => 38.436
    },
    'elevation_choctawhatchee_bay_fl_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -86.663,
      east_bound  => -85.831,
      south_bound => 30.242,
      north_bound => 30.591
    },
    'elevation_choptank_river_md_navd88' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -76.704,
      east_bound  => -75.782,
      south_bound => 38.374,
      north_bound => 39.053
    },
    'elevation_clay_creek_fl_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -83.521,
      east_bound  => -83.471,
      south_bound => 29.711,
      north_bound => 29.763
    },
    'elevation_clearwater_saint_joseph_sound_fl_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -82.907,
      east_bound  => -82.675,
      south_bound => 27.863,
      north_bound => 28.258
    },
    'elevation_cobscook_me_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -67.303,
      east_bound  => -66.894,
      south_bound => 44.737,
      north_bound => 45.027
    },
    'elevation_cocohatchee_river_fl_navd88' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -81.852,
      east_bound  => -81.794,
      south_bound => 26.258,
      north_bound => 26.332
    },
    'elevation_crooked_creek_fl_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -83.56,
      east_bound  => -83.513,
      south_bound => 29.738,
      north_bound => 29.761
    },
    'elevation_crystal_river_1_fl_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -82.799,
      east_bound  => -82.546,
      south_bound => 28.807,
      north_bound => 29.007
    },
    'elevation_crystal_river_2_fl_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -82.77,
      east_bound  => -82.617,
      south_bound => 28.707,
      north_bound => 28.851
    },
    'elevation_crystal_river_3_fl_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -82.73,
      east_bound  => -82.552,
      south_bound => 28.611,
      north_bound => 28.766
    },
    'elevation_cumberland_sound_fl_navd88' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -81.666,
      east_bound  => -81.378,
      south_bound => 30.609,
      north_bound => 31.0
    },
    'elevation_cumberland_sound_ga_navd88' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -81.666,
      east_bound  => -81.378,
      south_bound => 30.609,
      north_bound => 31.0
    },
    'elevation_currituck_sound_nc_navd88' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -76.429,
      east_bound  => -75.563,
      south_bound => 36.055,
      north_bound => 36.796
    },
    'elevation_currituck_sound_va_navd88' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -76.429,
      east_bound  => -75.563,
      south_bound => 36.055,
      north_bound => 36.796
    },
    'elevation_dallus_creek_fl_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -83.504,
      east_bound  => -83.456,
      south_bound => 29.707,
      north_bound => 29.729
    },
    'elevation_damariscaotta_me_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -69.781,
      east_bound  => -69.424,
      south_bound => 43.748,
      north_bound => 44.072
    },
    'elevation_days_creek_fl_mllw' => {
      begin_date  => '2007-10-31',
      end_date    => '2007-10-31',
      west_bound  => -83.41,
      east_bound  => -83.386,
      south_bound => 29.562,
      north_bound => 29.577
    },
    'elevation_daytona_beach_fl_mhw' => {
      begin_date  => '2007-10-31',
      end_date    => '2007-10-31',
      west_bound  => -81.4,
      east_bound  => -80.5,
      south_bound => 28.85,
      north_bound => 29.7
    },
    'elevation_deadman_bay_fl_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -83.475,
      east_bound  => -83.332,
      south_bound => 29.576,
      north_bound => 29.72
    },
    'elevation_deboy_sound_ga_navd88' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -81.513,
      east_bound  => -81.238,
      south_bound => 31.322,
      north_bound => 31.524
    },
    'elevation_delaware_bay_1_de_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -75.59,
      east_bound  => -74.672,
      south_bound => 38.666,
      north_bound => 39.379
    },
    'elevation_delaware_bay_1_nj_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -75.59,
      east_bound  => -74.672,
      south_bound => 38.666,
      north_bound => 39.379
    },
    'elevation_delaware_bay_2_de_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -75.786,
      east_bound  => -74.972,
      south_bound => 39.152,
      north_bound => 40.054
    },
    'elevation_delaware_bay_2_nj_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -75.786,
      east_bound  => -74.972,
      south_bound => 39.152,
      north_bound => 40.054
    },
    'elevation_delaware_bay_3_de_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -75.156,
      east_bound  => -74.678,
      south_bound => 39.837,
      north_bound => 40.272
    },
    'elevation_delaware_bay_3_nj_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -75.156,
      east_bound  => -74.678,
      south_bound => 39.837,
      north_bound => 40.272
    },
    'elevation_dewes_hamlin_creek_sc_navd88' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -79.854,
      east_bound  => -79.706,
      south_bound => 32.751,
      north_bound => 32.895
    },
    'elevation_duxbury_bay_ma_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -70.762,
      east_bound  => -70.538,
      south_bound => 41.904,
      north_bound => 42.089
    },
    'elevation_eaglenest_creek_fl_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -83.734,
      east_bound  => -83.72,
      south_bound => 29.94,
      north_bound => 29.954
    },
    'elevation_eastern_bay_md_navd88' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -76.394,
      east_bound  => -76.038,
      south_bound => 38.731,
      north_bound => 39.011
    },
    'elevation_eleven_prong_bay_fl_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -82.823,
      east_bound  => -82.774,
      south_bound => 29.071,
      north_bound => 29.108
    },
    'elevation_elk_river_md_navd88' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -76.016,
      east_bound  => -75.559,
      south_bound => 39.365,
      north_bound => 39.64
    },
    'elevation_essex_bay_ma_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -70.807,
      east_bound  => -70.66,
      south_bound => 42.608,
      north_bound => 42.704
    },
    'elevation_estero_bay_fl_navd88' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -81.947,
      east_bound  => -81.779,
      south_bound => 26.314,
      north_bound => 26.524
    },
    'elevation_fillman_bayou_fl_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -82.703,
      east_bound  => -82.674,
      south_bound => 28.401,
      north_bound => 28.421
    },
    'elevation_first_bay_fl_navd88' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -81.244,
      east_bound  => -80.98,
      south_bound => 25.474,
      north_bound => 25.666
    },
    'elevation_fish_creek_fl_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -83.586,
      east_bound  => -83.496,
      south_bound => 29.758,
      north_bound => 29.812
    },
    'elevation_florida_bay_fl_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -81.807,
      east_bound  => -80.189,
      south_bound => 24.399,
      north_bound => 25.503
    },
    'elevation_frenchmans_bay_me_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -68.505,
      east_bound  => -67.926,
      south_bound => 44.236,
      north_bound => 44.662
    },
    'elevation_gardiner_bay_peconic_bay_1_ny_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -72.672,
      east_bound  => -72.242,
      south_bound => 40.837,
      north_bound => 41.166
    },
    'elevation_gardiner_bay_peconic_bay_2_ny_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -72.448,
      east_bound  => -71.879,
      south_bound => 40.907,
      north_bound => 41.267
    },
    'elevation_glouchester_harbor_ma_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -70.741,
      east_bound  => -70.611,
      south_bound => 42.556,
      north_bound => 42.652
    },
    'elevation_goosefare_bay_me_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -70.467,
      east_bound  => -70.354,
      south_bound => 43.359,
      north_bound => 43.444
    },
    'elevation_gouldsboro_bay_me_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -68.081,
      east_bound  => -67.882,
      south_bound => 44.374,
      north_bound => 44.53
    },
    'elevation_great_bay_me_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -71.04,
      east_bound  => -70.567,
      south_bound => 42.919,
      north_bound => 43.279
    },
    'elevation_great_bay_nh_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -71.04,
      east_bound  => -70.567,
      south_bound => 42.919,
      north_bound => 43.279
    },
    'elevation_great_bay_nj_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -74.7,
      east_bound  => -74.243,
      south_bound => 39.41,
      north_bound => 39.697
    },
    'elevation_great_south_bay_ny_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -74.025,
      east_bound  => -72.497,
      south_bound => 40.31,
      north_bound => 41.061
    },
    'elevation_great_wilmico_fleets_bay_va_navd88' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -76.472,
      east_bound  => -76.259,
      south_bound => 37.773,
      north_bound => 37.888
    },
    'elevation_gunpowder_river_md_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -76.461,
      east_bound  => -76.239,
      south_bound => 39.274,
      north_bound => 39.448
    },
    'elevation_hampton_harbor_nh_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -70.909,
      east_bound  => -70.731,
      south_bound => 42.807,
      north_bound => 42.972
    },
    'elevation_haycock_harbor_baileys_mistake_me_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -67.103,
      east_bound  => -67.007,
      south_bound => 44.729,
      north_bound => 44.79
    },
    'elevation_honga_river_md_navd88' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -76.306,
      east_bound  => -76.034,
      south_bound => 38.225,
      north_bound => 38.426
    },
    'elevation_horseshoe_cove_fl_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -83.314,
      east_bound  => -83.164,
      south_bound => 29.371,
      north_bound => 29.517
    },
    'elevation_hudson_river_nj_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -74.682,
      east_bound  => -73.038,
      south_bound => 40.631,
      north_bound => 42.824
    },
    'elevation_hudson_river_ny_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -74.682,
      east_bound  => -73.038,
      south_bound => 40.631,
      north_bound => 42.824
    },
    'elevation_icw4_fl_navd88' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -82.456,
      east_bound  => -82.395,
      south_bound => 27.011,
      north_bound => 27.104
    },
    'elevation_indian_bay_fl_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -81.158,
      east_bound  => -79.836,
      south_bound => 27.009,
      north_bound => 28.912
    },
    'elevation_indian_bay_fl_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -82.688,
      east_bound  => -82.657,
      south_bound => 28.433,
      north_bound => 28.467
    },
    'elevation_indian_river_bay_de_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -75.313,
      east_bound  => -75.034,
      south_bound => 38.504,
      north_bound => 38.667
    },
    'elevation_indian_river_bay_fl_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -81.158,
      east_bound  => -79.836,
      south_bound => 27.009,
      north_bound => 28.912
    },
    'elevation_jamaica_bay_ny_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -73.97,
      east_bound  => -73.701,
      south_bound => 40.523,
      north_bound => 40.706
    },
    'elevation_james_river_va_navd88' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -77.634,
      east_bound  => -76.035,
      south_bound => 36.675,
      north_bound => 37.539
    },
    'elevation_jekyll_sound_saint_andrews_sound_ga_navd88' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -81.72,
      east_bound  => -81.384,
      south_bound => 30.919,
      north_bound => 31.148
    },
    'elevation_johnson_creek_fl_navd88' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -83.213,
      east_bound  => -83.175,
      south_bound => 29.378,
      north_bound => 29.403
    },
    'elevation_katama_bay_ma_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -70.533,
      east_bound  => -70.443,
      south_bound => 41.337,
      north_bound => 41.407
    },
    'elevation_kennebec_river_estuary_me_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -70.18,
      east_bound  => -69.532,
      south_bound => 43.697,
      north_bound => 44.359
    },
    'elevation_key_west_fl_mhw' => {
      begin_date  => '2011-09-01',
      end_date    => '2011-09-01',
      west_bound  => -82.18,
      east_bound  => -81.27,
      south_bound => 23.98,
      north_bound => 25.02
    },
    'elevation_key_west_fl_navd88' => {
      begin_date  => '2011-09-01',
      end_date    => '2011-09-01',
      west_bound  => -82.18,
      east_bound  => -81.27,
      south_bound => 23.98,
      north_bound => 25.02
    },
    'elevation_lake_worth_fl_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -80.223,
      east_bound  => -80.009,
      south_bound => 25.925,
      north_bound => 26.331
    },
    'elevation_lakes_reeds_bay_nj_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -74.565,
      east_bound  => -74.281,
      south_bound => 39.309,
      north_bound => 39.532
    },
    'elevation_little_kennebec_englishman_bay_me_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -67.754,
      east_bound  => -67.308,
      south_bound => 44.4,
      north_bound => 44.727
    },
    'elevation_little_narragansett_bay_ct_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -72.063,
      east_bound  => -71.812,
      south_bound => 41.229,
      north_bound => 41.419
    },
    'elevation_little_pine_1_bay_fl_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -82.691,
      east_bound  => -82.64,
      south_bound => 28.459,
      north_bound => 28.509
    },
    'elevation_little_river_sc_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -78.51,
      east_bound  => -78.5,
      south_bound => 33.821,
      north_bound => 33.915
    },
    'elevation_little_rocky_river_fl_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -83.358,
      east_bound  => -83.317,
      south_bound => 29.471,
      north_bound => 29.504
    },
    'elevation_little_spring_creek_fl_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -83.663,
      east_bound  => -83.639,
      south_bound => 29.897,
      north_bound => 29.916
    },
    'elevation_lockwoods_folly_river_nc_navd88' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -78.278,
      east_bound  => -78.139,
      south_bound => 33.897,
      north_bound => 33.983
    },
    'elevation_long_island_sound_ct_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -74.128,
      east_bound  => -71.672,
      south_bound => 40.446,
      north_bound => 42.085
    },
    'elevation_long_island_sound_ny_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -74.128,
      east_bound  => -71.672,
      south_bound => 40.446,
      north_bound => 42.085
    },
    'elevation_lows_bay_fl_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -82.823,
      east_bound  => -82.76,
      south_bound => 29.05,
      north_bound => 29.095
    },
    'elevation_loxahatchee_river_fl_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -80.319,
      east_bound  => -79.883,
      south_bound => 26.321,
      north_bound => 27.18
    },
    'elevation_lynnhaven_river_va_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -76.16,
      east_bound  => -75.955,
      south_bound => 36.802,
      north_bound => 36.941
    },
    'elevation_machias_bay_me_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -67.525,
      east_bound  => -67.162,
      south_bound => 44.521,
      north_bound => 44.782
    },
    'elevation_maddaket_harbor_ma_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -70.272,
      east_bound  => -70.17,
      south_bound => 41.258,
      north_bound => 41.319
    },
    'elevation_magothy_river_md_navd88' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -76.567,
      east_bound  => -76.402,
      south_bound => 39.028,
      north_bound => 39.121
    },
    'elevation_manasquan_river_nj_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -74.125,
      east_bound  => -74.02,
      south_bound => 40.071,
      north_bound => 40.14
    },
    'elevation_merrimack_river_ma_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -71.211,
      east_bound  => -70.759,
      south_bound => 42.637,
      north_bound => 42.908
    },
    'elevation_miami_fl_mhw' => {
      begin_date  => '2015-01-21',
      end_date    => '2015-01-21',
      west_bound  => -80.41,
      east_bound  => -79.35,
      south_bound => 25.28,
      north_bound => 26.32
    },
    'elevation_miami_fl_navd88' => {
      begin_date  => '2015-01-21',
      end_date    => '2015-01-21',
      west_bound  => -80.41,
      east_bound  => -79.35,
      south_bound => 25.28,
      north_bound => 26.32
    },
    'elevation_middle_and_black_rivers_md_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -76.565,
      east_bound  => -76.305,
      south_bound => 39.196,
      north_bound => 39.371
    },
    'elevation_mojack_bay_va_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -76.567,
      east_bound  => -76.219,
      south_bound => 37.225,
      north_bound => 37.516
    },
    'elevation_montauk_ny_mhw' => {
      begin_date  => '2007-09-28',
      end_date    => '2007-09-28',
      west_bound  => -72.6,
      east_bound  => -71.5,
      south_bound => 40.6,
      north_bound => 41.4
    },
    'elevation_moose_cove_me_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -67.129,
      east_bound  => -67.077,
      south_bound => 44.721,
      north_bound => 44.752
    },
    'elevation_morehead_city_nc_mhw' => {
      begin_date  => '2011-04-15',
      end_date    => '2011-04-15',
      west_bound  => -77.27,
      east_bound  => -76,
      south_bound => 34.37,
      north_bound => 35.57
    },
    'elevation_morehead_city_nc_navd88' => {
      begin_date  => '2011-04-15',
      end_date    => '2011-04-15',
      west_bound  => -77.27,
      east_bound  => -76,
      south_bound => 34.37,
      north_bound => 35.57
    },
    'elevation_myrtle_beach_sc_mhw' => {
      begin_date  => '2006-05-26',
      end_date    => '2006-05-26',
      west_bound  => -79.2,
      east_bound  => -78.4,
      south_bound => 33.25,
      north_bound => 33.95
    },
    'elevation_murrells_inlet_sc_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -79.088,
      east_bound  => -78.982,
      south_bound => 33.49,
      north_bound => 33.592
    },
    'elevation_muscongus_bay_me_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -69.612,
      east_bound  => -69.064,
      south_bound => 43.764,
      north_bound => 44.162
    },
    'elevation_nantucket_ma_mhw' => {
      begin_date  => '2008-10-10',
      end_date    => '2008-10-10',
      west_bound  => -70.67,
      east_bound  => -69.49,
      south_bound => 40.81,
      north_bound => 41.71
    },
    'elevation_nantucket_harbor_ma_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -70.137,
      east_bound  => -69.966,
      south_bound => 41.252,
      north_bound => 41.388
    },
    'elevation_naples_bay_fl_navd88' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -81.813,
      east_bound  => -81.739,
      south_bound => 26.05,
      north_bound => 26.173
    },
    'elevation_narragansett_bay_1_ma_navd88' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -71.441,
      east_bound  => -71.284,
      south_bound => 41.712,
      north_bound => 41.877
    },
    'elevation_narragansett_bay_1_ri_navd88' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -71.441,
      east_bound  => -71.284,
      south_bound => 41.712,
      north_bound => 41.877
    },
    'elevation_narragansett_bay_2_ma_navd88' => {
      begin_date  => '2020-06-18',
      end_date  => '2020-06-18',
      west_bound  => -71.197,
      east_bound  => -71.041,
      south_bound => 41.684,
      north_bound => 41.911
    },
    'elevation_narragansett_bay_2_ri_navd88' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -71.197,
      east_bound  => -71.041,
      south_bound => 41.684,
      north_bound => 41.911
    },
    'elevation_narragansett_bay_3_ma_navd88' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -71.573,
      east_bound  => -71.033,
      south_bound => 41.389,
      north_bound => 41.839
    },
    'elevation_narragansett_bay_3_ri_navd88' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -71.573,
      east_bound  => -71.033,
      south_bound => 41.389,
      north_bound => 41.839
    },
    'elevation_narraguagas_bay_me_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -68.028,
      east_bound  => -67.682,
      south_bound => 44.360,
      north_bound => 44.654
    },
    'elevation_nassau_sound_fl_navd88' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -81.732,
      east_bound  => -81.367,
      south_bound => 30.406,
      north_bound => 30.660
    },
    'elevation_navesink_nj_navd88' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -74.111,
      east_bound  => -73.963,
      south_bound => 40.296,
      north_bound => 40.422
    },
    'elevation_nc_intercoastal_waterway_nc_navd88' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -77.946,
      east_bound  => -77.366,
      south_bound => 33.952,
      north_bound => 34.596
    },
    'elevation_new_river_nc_navd88' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -77.508,
      east_bound  => -77.247,
      south_bound => 34.507,
      north_bound => 34.767
    },
    'elevation_new_river_sc_mhw' => {
      begin_date  => '2020-06-22',
      end_date    => '2020-06-122',
      west_bound  => -81.07,
      east_bound  => -80.62,
      south_bound => 31.968,
      north_bound => 32.325
    },
    'elevation_new_york_harbor_1_nj_navd88' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -74.124,
      east_bound  => -73.761,
      south_bound => 40.568,
      north_bound => 40.902
    },
    'elevation_new_york_harbor_1_ny_navd88' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -74.124,
      east_bound  => -73.761,
      south_bound => 40.568,
      north_bound => 40.902
    },
    'elevation_new_york_harbor_2_nj_navd88' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -74.253,
      east_bound  => -73.880,
      south_bound => 40.370,
      north_bound => 40.653
    },
    'elevation_new_york_harbor_2_ny_navd88' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -74.253,
      east_bound  => -73.880,
      south_bound => 40.370,
      north_bound => 40.653
    },
    'elevation_newark_bay_kill_van_kull_nj_navd88' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -74.204,
      east_bound  => -73.976,
      south_bound => 40.622,
      north_bound => 40.951
    },
    'elevation_newark_bay_kill_van_kull_ny_navd88' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -74.204,
      east_bound  => -73.976,
      south_bound => 40.622,
      north_bound => 40.951
    },
    'elevation_nj_coastal_bays_great_egg_harbor_nj_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -75.038,
      east_bound  => -74.342,
      south_bound => 38.895,
      north_bound => 39.522
    },
    'elevation_north_edisto_river_sc_navd88' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -80.416,
      east_bound  => -80.057,
      south_bound => 32.509,
      north_bound => 32.775
    },
    'elevation_north_south_rivers_ma_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -70.802,
      east_bound  => -70.596,
      south_bound => 42.056,
      north_bound => 42.224
    },
    'elevation_north_south_santee_rivers_sc_navd88' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -79.539,
      east_bound  => -79.198,
      south_bound => 33.088,
      north_bound => 33.269
    },
    'elevation_northeast_atlantic_msl' => {
      begin_date  => '1998-11-16',
      end_date    => '1998-11-16',
      west_bound  => -80,
      east_bound  => -66,
      south_bound => 40,
      north_bound => 48
    },
    'elevation_northeast_river_md_navd88' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -76.012,
      east_bound  => -75.940,
      south_bound => 39.525,
      north_bound => 39.608
    },
    'elevation_ocean_city_md_mhw' => {
      begin_date  => '2009-10-16',
      end_date    => '2009-10-16',
      west_bound  => -75.58,
      east_bound  => -74.71,
      south_bound => 37.68,
      north_bound => 38.87
    },
    'elevation_ossabaw_sound_ga_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -81.408,
      east_bound  => -80.900,
      south_bound => 31.702,
      north_bound => 32.129
    },
    'elevation_otter_cove_me_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -68.221,
      east_bound  => -68.168,
      south_bound => 44.292,
      north_bound => 44.336
    },
    'elevation_outer_doctors_bay_fl_navd88' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -81.821,
      east_bound  => -81.800,
      south_bound => 26.163,
      north_bound => 26.208
    },
    'elevation_oyster_creek_fl_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -83.577,
      east_bound  => -83.564,
      south_bound => 29.765,
      north_bound => 29.776
    },
    'elevation_palm_beach_fl_mhw' => {
      begin_date  => '2010-09-16',
      end_date    => '2010-09-16',
      west_bound  => -80.36,
      east_bound  => -79.42,
      south_bound => 26.29,
      north_bound => 27.31
    },
    'elevation_palm_beach_fl_navd88' => {
      begin_date  => '2010-09-16',
      end_date    => '2010-09-16',
      west_bound  => -80.36,
      east_bound  => -79.42,
      south_bound => 26.29,
      north_bound => 27.31
    },
    'elevation_pamlico_sound_nc_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -77.545,
      east_bound  => -75.173,
      south_bound => 34.537,
      north_bound => 36.231
    },
    'elevation_panama_city_fl_mhw' => {
      begin_date  => '2010-07-01',
      end_date    => '2010-07-01',
      west_bound  => -86.1,
      east_bound  => -85.2,
      south_bound => 29.55,
      north_bound => 30.5
    },
    'elevation_panama_city_fl_navd88' => {
      begin_date  => '2010-07-01',
      end_date    => '2010-07-01',
      west_bound  => -86.1,
      east_bound  => -85.2,
      south_bound => 29.55,
      north_bound => 30.5
    },
    'elevation_passamaquoddy_bay_st_croix_river_me_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -67.377,
      east_bound  => -66.938,
      south_bound => 44.909,
      north_bound => 45.244
    },
    'elevation_patapsco_river_md_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -76.689,
      east_bound  => -76.390,
      south_bound => 39.105,
      north_bound => 39.326
    },
    'elevation_patuxent_river_md_navd88' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -76.810,
      east_bound  => -76.319,
      south_bound => 38.283,
      north_bound => 38.745
    },
    'elevation_penobscot_bay_me_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -69.454,
      east_bound  => -68.106,
      south_bound => 43.830,
      north_bound => 44.936
    },
    'elevation_perdido_bay_al_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -87.655,
      east_bound  => -87.308,
      south_bound => 30.247,
      north_bound => 30.538
    },
    'elevation_perdido_bay_fl_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -87.655,
      east_bound  => -87.308,
      south_bound => 30.247,
      north_bound => 30.538
    },
    'elevation_piankatank_river_va_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -76.616,
      east_bound  => -76.245,
      south_bound => 37.419,
      north_bound => 37.633
    },
    'elevation_pithlachascotee_river_fl_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -82.752,
      east_bound  => -82.704,
      south_bound => 28.231,
      north_bound => 28.289
    },
    'elevation_pleasant_bay_me_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -67.875,
      east_bound  => -67.642,
      south_bound => 44.459,
      north_bound => 44.672
    },
    'elevation_pleasant_bay_ma_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -70.051,
      east_bound  => -69.878,
      south_bound => 41.626,
      north_bound => 41.797
    },
    'elevation_plum_island_sound_ma_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -70.917,
      east_bound  => -70.709,
      south_bound => 42.651,
      north_bound => 42.814
    },
    'elevation_ponce_de_leon_inlet_fl_navd88' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -81.290,
      east_bound  => -80.488,
      south_bound => 28.641,
      north_bound => 29.506
    },
    'elevation_poquoson_back_rivers_va_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -76.504,
      east_bound  => -76.227,
      south_bound => 37.001,
      north_bound => 37.229
    },
    'elevation_portland_me_mhw' => {
      begin_date  => '2008-12-31',
      end_date    => '2008-12-31',
      west_bound  => -70.74,
      east_bound  => -69.63,
      south_bound => 43,
      north_bound => 43.99
    },
    'elevation_potomac_river_md_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -77.627,
      east_bound  => -75.962,
      south_bound => 37.717,
      north_bound => 39.101
    },
    'elevation_potomac_river_va_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -77.627,
      east_bound  => -75.962,
      south_bound => 37.717,
      north_bound => 39.101
    },
    'elevation_prospect_harbor_me_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -68.043,
      east_bound  => -67.957,
      south_bound => 44.372,
      north_bound => 44.429
    },
    'elevation_quoddy_narrows_me_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -67.032,
      east_bound  => -66.910,
      south_bound => 44.791,
      north_bound => 44.878
    },
    'elevation_rappahannock_river_va_navd88' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -77.614,
      east_bound  => -76.155,
      south_bound => 37.548,
      north_bound => 38.321
    },
    'elevation_raritan_bay_nj_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -74.516,
      east_bound  => -74.098,
      south_bound => 40.372,
      north_bound => 40.581
    },
    'elevation_raritan_bay_ny_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -74.516,
      east_bound  => -74.098,
      south_bound => 40.372,
      north_bound => 40.581
    },
    'elevation_rehoboth_bay_de_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -75.246,
      east_bound  => -75.015,
      south_bound => 38.588,
      north_bound => 38.809
    },
    'elevation_rock_1_bay_fl_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -82.669,
      east_bound  => -82.635,
      south_bound => 28.538,
      north_bound => 28.585
    },
    'elevation_saco_bay_me_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -70.510,
      east_bound  => -70.139,
      south_bound => 43.382,
      north_bound => 43.637
    },
    'elevation_saint_andrew_bay_fl_navd88' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -85.971,
      east_bound  => -85.016,
      south_bound => 29.745,
      north_bound => 30.390
    },
    'elevation_saint_augustine_fl_navd88' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -81.543,
      east_bound  => -81.039,
      south_bound => 29.544,
      north_bound => 30.271
    },
    'elevation_saint_johns_river_fl_navd88' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -81.879,
      east_bound  => -81.339,
      south_bound => 29.946,
      north_bound => 30.522
    },
    'elevation_saint_joseph_bay_fl_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -85.442,
      east_bound  => -85.217,
      south_bound => 29.661,
      north_bound => 29.931
    },
    'elevation_salt_creek_fl_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -83.595,
      east_bound  => -83.576,
      south_bound => 29.817,
      north_bound => 29.836
    },
    'elevation_salt_creek_2_fl_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -83.480,
      east_bound  => -83.416,
      south_bound => 29.678,
      north_bound => 29.735
    },
    'elevation_salem_sound_ma_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -70.965,
      east_bound  => -70.723,
      south_bound => 42.457,
      north_bound => 42.608
    },
    'elevation_sanders_creek_fl_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -83.208,
      east_bound  => -83.167,
      south_bound => 29.359,
      north_bound => 29.402
    },
    'elevation_sapelo_sound_ga_navd88' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -81.444,
      east_bound  => -81.135,
      south_bound => 31.460,
      north_bound => 31.702
    },
    'elevation_sarasota_bay_fl_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -82.753,
      east_bound  => -82.376,
      south_bound => 27.069,
      north_bound => 27.504
    },
    'elevation_sassafrass_river_md_navd88' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -76.115,
      east_bound  => -75.805,
      south_bound => 39.308,
      north_bound => 39.428
    },
    'elevation_savannah_ga_mhw' => {
      begin_date  => '2006-12-05',
      end_date    => '2006-12-05',
      west_bound  => -81.35,
      east_bound  => -80.35,
      south_bound => 31.5,
      north_bound => 32.45
    },
    'elevation_savanna_river_ga_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -81.293,
      east_bound  => -80.843,
      south_bound => 31.966,
      north_bound => 32.267
    },
    'elevation_sengekontacket_pond_ma_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -70.592,
      east_bound  => -70.516,
      south_bound => 41.395,
      north_bound => 41.463
    },
    'elevation_severn_river_md_navd88' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -76.651,
      east_bound  => -76.426,
      south_bound => 38.940,
      north_bound => 39.086
    },
    'elevation_sewee_bay_sc_navd88' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -79.678,
      east_bound  => -79.592,
      south_bound => 32.884,
      north_bound => 32.975
    },
    'elevation_shallotte_river_nc_navd88' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -78.423,
      east_bound  => -78.304,
      south_bound => 33.879,
      north_bound => 33.985
    },
    'elevation_shark_river_nj_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -74.068,
      east_bound  => -74.000,
      south_bound => 40.166,
      north_bound => 40.212
    },
    'elevation_sheepscot_bay_me_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -69.975,
      east_bound  => -69.459,
      south_bound => 43.651,
      north_bound => 44.100
    },
    'elevation_shinnecock_bay_ny_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -72.625,
      east_bound  => -72.383,
      south_bound => 40.775,
      north_bound => 40.922
    },
    'elevation_shrewsbury_river_nj_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -74.074,
      east_bound  => -73.954,
      south_bound => 40.286,
      north_bound => 40.381
    },
    'elevation_sink_creek_fl_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -83.405,
      east_bound  => -83.389,
      south_bound => 29.553,
      north_bound => 29.567
    },
    'elevation_somes_sound_sw_bay_me_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -68.390,
      east_bound  => -68.158,
      south_bound => 44.216,
      north_bound => 44.402
    },
    'elevation_south_edisto_river_sc_navd88' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -80.491,
      east_bound  => -80.235,
      south_bound => 32.464,
      north_bound => 32.709
    },
    'elevation_southeast_atlantic_msl' => {
      begin_date  => '1998-11-16',
      end_date    => '1998-11-16',
      west_bound  => -85,
      east_bound  => -68,
      south_bound => 31,
      north_bound => 40
    },
    'elevation_southern_cap_coastal_ponds_ma_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -70.583,
      east_bound  => -70.194,
      south_bound => 41.489,
      north_bound => 41.723
    },
      'elevation_spring_warrior_creek_fl_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -83.688,
      east_bound  => -83.647,
      south_bound => 29.910,
      north_bound => 29.949
    },
    'elevation_storter_bay_fl_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -81.379,
      east_bound  => -81.078,
      south_bound => 25.581,
      north_bound => 25.844
    },
    'elevation_suwanee_river_fl_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -83.227,
      east_bound  => -82.842,
      south_bound => 29.263,
      north_bound => 29.717
    },
    'elevation_tampa_bay_1_fl_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -82.544,
      east_bound  => -82.301,
      south_bound => 27.770,
      north_bound => 28.038
    },
    'elevation_tampa_bay_2_fl_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -82.779,
      east_bound  => -82.476,
      south_bound => 27.784,
      north_bound => 28.077
    },
    'elevation_tampa_bay_3_fl_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -82.903,
      east_bound  => -82.327,
      south_bound => 27.372,
      north_bound => 27.914
    },
    'elevation_tampa_bay_fl_mhw' => {
      begin_date  => '2014-09-11',
      end_date    => '2014-09-1',
      west_bound  => -83.41,
      east_bound  => -82.37,
      south_bound => 27.22,
      north_bound => 28.26
    },
    'elevation_tampa_bay_fl_navd88' => {
      begin_date  => '2014-09-11',
      end_date    => '2014-09-1',
      west_bound  => -83.41,
      east_bound  => -82.37,
      south_bound => 27.22,
      north_bound => 28.26
    },
    'elevation_ten_thousand_islands_marco_bay_fl_navd88' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -81.834,
      east_bound  => -81.241,
      south_bound => 25.730,
      north_bound => 26.082
    },
    'elevation_virginia_beach_va_mhw' => {
      begin_date  => '2007-01-11',
      end_date    => '2007-01-11',
      west_bound  => -76.6,
      east_bound  => -75.4,
      south_bound => 36.45,
      north_bound => 37.5
    },
    'elevation_waccasassa_river_fl_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -82.938,
      east_bound  => -82.642,
      south_bound => 29.052,
      north_bound => 29.252
    },
    'elevation_weeki_wachee_bay_fl_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -82.675,
      east_bound  => -82.617,
      south_bound => 28.511,
      north_bound => 28.560
    },
    'elevation_whackup_creek_fl_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -83.344,
      east_bound  => -83.314,
      south_bound => 29.472,
      north_bound => 29.496
    },
    'elevation_whitewater_bay_fl_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -81.206,
      east_bound  => -80.793,
      south_bound => 25.086,
      north_bound => 25.502
    },
    'elevation_withlacoochee_bay_fl_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -82.818,
      east_bound  => -82.704,
      south_bound => 28.989,
      north_bound => 29.087
    },
    'elevation_yates_creek_fl_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -83.656,
      east_bound  => -83.647,
      south_bound => 29.894,
      north_bound => 29.903
    },

    'nadp.gridded.wet_nh4' => {
      begin_date  => '1985-01-01',
      end_date    => '2014-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'nadp.gridded.wet_no3' => {
      begin_date  => '1985-01-01',
      end_date    => '2014-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'nadp.gridded.wet_totn' => {
      begin_date  => '1985-01-01',
      end_date    => '2014-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'nadpcmaq.gridded.dry_nh3' => {
      begin_date  => '2000-01-01',
      end_date    => '2014-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'nadpcmaq.gridded.dry_nh4' => {
      begin_date  => '2000-01-01',
      end_date    => '2014-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'nadpcmaq.gridded.dry_no3' => {
      begin_date  => '2000-01-01',
      end_date    => '2014-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'nadpcmaq.gridded.dry_totn' => {
      begin_date  => '2000-01-01',
      end_date    => '2014-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'nadpcmaq.gridded.dry_totnp' => {
      begin_date  => '2000-01-01',
      end_date    => '2014-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'nadpcmaq.gridded.dry_totun' => {
      begin_date  => '2000-01-01',
      end_date    => '2014-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'nadpcmaq.gridded.dry_totunp' => {
      begin_date  => '2000-01-01',
      end_date    => '2014-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'nadpcmaq.gridded.dry_net_nh3' => {
      begin_date  => '2000-01-01',
      end_date    => '2014-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'nadpcmaq.gridded.dry_oxidized_totn' => {
      begin_date  => '2000-01-01',
      end_date    => '2014-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'nadpcmaq.gridded.dry_oxidized_totnp' => {
      begin_date  => '2000-01-01',
      end_date    => '2014-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'nadpcmaq.gridded.dry_reduced_totn' => {
      begin_date  => '2000-01-01',
      end_date    => '2014-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'nadpcmaq.gridded.dry_reduced_totnp' => {
      begin_date  => '2000-01-01',
      end_date    => '2014-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'nadpcmaq.gridded.wet_nh4' => {
      begin_date  => '2000-01-01',
      end_date    => '2014-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'nadpcmaq.gridded.wet_no3' => {
      begin_date  => '2000-01-01',
      end_date    => '2014-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'nadpcmaq.gridded.wet_totn' => {
      begin_date  => '2000-01-01',
      end_date    => '2014-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'nadpcmaq.gridded.oxidized_totn' => {
      begin_date  => '2000-01-01',
      end_date    => '2014-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'nadpcmaq.gridded.oxidized_totnp' => {
      begin_date  => '2000-01-01',
      end_date    => '2014-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'nadpcmaq.gridded.reduced_totn' => {
      begin_date  => '2000-01-01',
      end_date    => '2014-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'nadpcmaq.gridded.reduced_totnp' => {
      begin_date  => '2000-01-01',
      end_date    => '2014-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'nadpcmaq.gridded.total_totn' => {
      begin_date  => '2000-01-01',
      end_date    => '2014-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },

    'temperature_minimum' => {
      begin_date  => '1971-01-01',
      end_date    => '2000-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'temperature_maximum' => {
      begin_date  => '1971-01-01',
      end_date    => '2000-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'precipitation_mean' => {
      begin_date  => '1971-01-01',
      end_date    => '2000-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'monthly_temperature_minimum' => {
      begin_date  => '1971-01-01',
      end_date    => '2000-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'monthly_temperature_maximum' => {
      begin_date  => '1971-01-01',
      end_date    => '2000-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'monthly_precipitation_mean' => {
      begin_date  => '1971-01-01',
      end_date    => '2000-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'annual_temperature_minimum' => {
      begin_date  => '1971-01-01',
      end_date    => '2000-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'annual_temperature_maximum' => {
      begin_date  => '1971-01-01',
      end_date    => '2000-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'annual_precipitation_mean' => {
      begin_date  => '1971-01-01',
      end_date    => '2000-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'seagrass' => {
      begin_date  => '1979-01-01',
      end_date    => '2009-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'nitrogen.source_gridded.crops' => {
      begin_date  => '2006-01-01',
      end_date    => '2006-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'nitrogen.source_gridded.fertilizer' => {
      begin_date  => '2006-01-01',
      end_date    => '2006-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'nitrogen.source_gridded.manure' => {
      begin_date  => '2006-01-01',
      end_date    => '2006-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'wetlands' => {
      begin_date  => '1981-01-01',
      end_date    => '2011-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'population_iclus1' => {
      begin_date  => '2010-01-01',
      end_date    => '2090-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'population_iclus2' => {
      begin_date  => '2020-01-01',
      end_date    => '2100-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'land_use_iclus_watershed' => {
      begin_date  => '2000-01-01',
      end_date    => '2100-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'land_use_nlcd_2011_watershed' => {
      begin_date  => '2011-01-01',
      end_date    => '2011-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'land_use_nlcd_2006_watershed' => {
      begin_date  => '2006-01-01',
      end_date    => '2006-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'soil' => {
      begin_date  => '2012-01-01',
      end_date    => '2012-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'nitrogen.estuary_cmaq' => { # Old
      begin_date  => '2002-01-01',
      end_date    => '2008-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'dry_nh3.estuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'dry_natural_nh3.estuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'dry_nhx.estuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'dry_so2.estuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'dry_oxidized_nox.estuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'dry_oxidized_organic_nitrogen.estuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'dry_oxidized_no3.estuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'dry_oxidized_nitrogen.estuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'dry_reduced_nitrogen.estuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'dry_sulfur.estuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'chloride.estuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'total_nitrogen.estuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'oxidized_nitrogen.estuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'reduced_nitrogen.estuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'sulfur.estuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'wet_nh3.estuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'wet_chloride.estuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'wet_so4.estuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'wet_oxidized_no3.estuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'wet_oxidized_nox.estuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'wet_organic_nitrogen.estuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'wet_oxidized_nitrogen.estuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'wet_reduced_nitrogen.estuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'wet_sulfur.estuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'nitrogen.subestuary_cmaq' => { # Old
      begin_date  => '2002-01-01',
      end_date    => '2008-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'dry_nh3.subestuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'dry_natural_nh3.subestuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'dry_nhx.subestuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'dry_so2.subestuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'dry_oxidized_nox.subestuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'dry_oxidized_organic_nitrogen.subestuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'dry_oxidized_no3.subestuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'dry_oxidized_nitrogen.subestuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'dry_reduced_nitrogen.subestuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'dry_sulfur.subestuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'chloride.subestuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'total_nitrogen.subestuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'oxidized_nitrogen.subestuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'reduced_nitrogen.subestuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'sulfur.subestuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'wet_nh3.subestuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'wet_chloride.subestuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'wet_so4.subestuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'wet_oxidized_no3.subestuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'wet_oxidized_nox.subestuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'wet_organic_nitrogen.subestuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'wet_oxidized_nitrogen.subestuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'wet_reduced_nitrogen.subestuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'wet_sulfur.subestuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'nitrogen.load_estuary_sparrow_1992' => {
      begin_date  => '1992-01-01',
      end_date    => '1992-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'nitrogen.load_estuary_sparrow_2002_mrb1' => {
      begin_date  => '1992-01-01',
      end_date    => '1992-12-31',
      west_bound  => -77.5,
      east_bound  => -67,
      south_bound => 36.7,
      north_bound => 45.5
    },
    'nitrogen.load_estuary_sparrow_2002_mrb2' => {
      begin_date  => '1992-01-01',
      end_date    => '1992-12-31',
      west_bound  => -77.6,
      east_bound  => -75.5,
      south_bound => 35.6,
      north_bound => 36.6
    },
    'nitrogen.load_estuary_sparrow_2011_mrb1mrb2' => {
      begin_date  => '2011-01-01',
      end_date    => '2011-12-31',
      west_bound  => -77.6,
      east_bound  => -67,
      south_bound => 35.6,
      north_bound => 45.5
    },
    'nitrogen.load_subestuary_sparrow_1992' => {
      begin_date  => '1992-01-01',
      end_date    => '1992-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'nitrogen.load_subestuary_sparrow_2002_mrb1' => {
      begin_date  => '1992-01-01',
      end_date    => '1992-12-31',
      west_bound  => -77.5,
      east_bound  => -67,
      south_bound => 36.7,
      north_bound => 45.5
    },
    'nitrogen.load_subestuary_sparrow_2002_mrb2' => {
      begin_date  => '1992-01-01',
      end_date    => '1992-12-31',
      west_bound  => -77.6,
      east_bound  => -75.5,
      south_bound => 35.6,
      north_bound => 36.6
    },
    'total_nitrogen.source_watershed_nadp' => {
      begin_date  => '2000-01-01',
      end_date    => '2014-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'wet_nitrogen.source_watershed_nadp' => {
      begin_date  => '2000-01-01',
      end_date    => '2014-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'dry_nitrogen.source_watershed_nadp' => {
      begin_date  => '2000-01-01',
      end_date    => '2014-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'oxidized_nitrogen.source_watershed_nadp' => {
      begin_date  => '2000-01-01',
      end_date    => '2014-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'reduced_nitrogen.source_watershed_nadp' => {
      begin_date  => '2000-01-01',
      end_date    => '2014-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'dry_oxidized_nitrogen.source_watershed_nadp' => {
      begin_date  => '2000-01-01',
      end_date    => '2014-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'dry_reduced_nitrogen.source_watershed_nadp' => {
      begin_date  => '2000-01-01',
      end_date    => '2014-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'wet_no3.source_watershed_nadp' => {
      begin_date  => '1994-01-01',
      end_date    => '2014-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'wet_nh4.source_watershed_nadp' => {
      begin_date  => '1994-01-01',
      end_date    => '2014-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'dry_no3.source_watershed_nadp' => {
      begin_date  => '2000-01-01',
      end_date    => '2014-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'dry_nh4.source_watershed_nadp' => {
      begin_date  => '2000-01-01',
      end_date    => '2014-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'nitrogen.source_watershed_cmaq' => { # Old
      begin_date  => '2008-01-01',
      end_date    => '2011-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'dry_nh3.source_watershed_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'dry_natural_nh3.source_watershed_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'dry_nhx.source_watershed_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'dry_so2.source_watershed_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'dry_oxidized_nox.source_watershed_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'dry_oxidized_organic_nitrogen.source_watershed_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'dry_oxidized_no3.source_watershed_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'dry_oxidized_nitrogen.source_watershed_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'dry_reduced_nitrogen.source_watershed_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'dry_sulfur.source_watershed_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'chloride.source_watershed_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'total_nitrogen.source_watershed_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'oxidized_nitrogen.source_watershed_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'reduced_nitrogen.source_watershed_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'sulfur.source_watershed_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'wet_nh3.source_watershed_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'wet_chloride.source_watershed_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'wet_so4.source_watershed_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'wet_oxidized_no3.source_watershed_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'wet_oxidized_nox.source_watershed_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'wet_organic_nitrogen.source_watershed_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'wet_oxidized_nitrogen.source_watershed_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'wet_reduced_nitrogen.source_watershed_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'wet_sulfur.source_watershed_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'nitrogen.source_watershed_point' => {
      begin_date  => '2007-01-01',
      end_date    => '2011-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'nitrogen.source_watershed_nonpoint' => {
      begin_date  => '2006-01-01',
      end_date    => '2006-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'phosphorus.source_watershed_point' => {
      begin_date  => '2007-01-01',
      end_date    => '2015-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'sediment_nca' => {
      begin_date  => '1990-01-01',
      end_date    => '2006-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'sediment_nca_1990-2006' => {
      begin_date  => '1990-01-01',
      end_date    => '2006-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'sediment_nca_2015' => {
      begin_date  => '2015-01-01',
      end_date    => '2015-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'sediment_calculated' => {
      begin_date  => '1990-01-01',
      end_date    => '2006-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'sediment_extracted' => {
      begin_date  => '1990-01-01',
      end_date    => '2006-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'sediment_parsed' => {
      begin_date  => '1990-01-01',
      end_date    => '2006-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'sediment_kriged' => {
      begin_date  => '2005-01-01',
      end_date    => '2014-12-31',
      west_bound  => -77,
      east_bound  => -64,
      south_bound => 35,
      north_bound => 46
    },
    'coastal_vulnerability' => {
      begin_date  => '1976-01-01',
      end_date    => '1995-12-31',
      west_bound  => -82,
      east_bound  => -67,
      south_bound => 24,
      north_bound => 46
    },
    'sensitivity_estuary' => {
      begin_date  => '2016-10-05',
      end_date    => '2016-10-05',
      west_bound  => -82,
      east_bound  => -67,
      south_bound => 24,
      north_bound => 46
    },
    'groundwater_contact_catchment' => {
      begin_date  => '2002-01-01',
      end_date    => '2002-12-31',
      west_bound  => -82,
      east_bound  => -67,
      south_bound => 24,
      north_bound => 46
    },
    'groundwater_contact_watershed' => {
      begin_date  => '2002-01-01',
      end_date    => '2002-12-31',
      west_bound  => -82,
      east_bound  => -67,
      south_bound => 24,
      north_bound => 46
    },
    'stream_discharge_nhd_line' => {
      begin_date  => '1952-01-01',
      end_date    => '2017-12-31',
      west_bound  => -85,
      east_bound  => -67,
      south_bound => 25,
      north_bound => 45
    },
    'stream_discharge_usgs_point' => {
      begin_date  => '1952-01-01',
      end_date    => '2017-12-31',
      west_bound  => -88,
      east_bound  => -67,
      south_bound => 27,
      north_bound => 45
    },
    'tide_point' => {
      begin_date  => '1950-01-01',
      end_date    => '2018-12-31',
      west_bound  => -88,
      east_bound  => -67,
      south_bound => 27,
      north_bound => 45
    },
    'tide_current' => {
      begin_date  => '2017-01-01',
      end_date    => '2019-12-31',
      west_bound  => -88,
      east_bound  => -67,
      south_bound => 24,
      north_bound => 45
    },
    'longshore_current' => {
      begin_date  => '2012-01-01',
      end_date    => '2020-12-31',
      west_bound  => -88,
      east_bound  => -67,
      south_bound => 24,
      north_bound => 45
    },
    'estuary_flushing' => {
      begin_date  => '1950-01-01',
      end_date    => '2020-12-31',
      west_bound  => -88,
      east_bound  => -67,
      south_bound => 24,
      north_bound => 45
    },
    'salinity_point' => {
      begin_date  => '1990-01-01',
      end_date    => '2018-12-31',
      west_bound  => -88,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 46
    },
    'flowlines_upper_ct_river' => {
      begin_date  => '2012-01-01',
      end_date    => '2018-12-31',
      west_bound  => -72.993,
      east_bound  => -71.102,
      south_bound => 42.76,
      north_bound => 45.29
    },
   'gi_bmp_installations_new_england' => {
      begin_date  => '2012-12-31',
      end_date    => '2015-06-12',
      west_bound  => -74.72,
      east_bound  => -67.36,
      south_bound => 40.44,
      north_bound => 47.54
    },
    'impervious_nhdplus_new_england' => {
      begin_date  => '2014-01-01',
      end_date    => '2014-01-01',
      west_bound  => -74.68,
      east_bound  => -65.7,
      south_bound => 40.04,
      north_bound => 48.782
    },
    'impervious_ball_mountain_vt' => {
      begin_date  => '2014-11-30',
      end_date    => '2014-11-30',
      west_bound  => -72.96927,
      east_bound  => -72.74163,
      south_bound => 43.00982,
      north_bound => 43.16081
    },
    'impervious_black_river_vt' => {
      begin_date  => '2014-11-30',
      end_date    => '2014-11-30',
      west_bound  => -72.89554,
      east_bound  => -72.34710,
      south_bound => 43.23006,
      north_bound => 43.57497
    },
    'impervious_bloody_brook_vt' => {
      begin_date  => '2014-11-30',
      end_date    => '2014-11-30',
      west_bound  => -72.43119,
      east_bound  => -72.05692,
      south_bound => 43.60100,
      north_bound => 43.84380
    },
    'impervious_branch_brook_vt' => {
      begin_date  => '2014-11-30',
      end_date    => '2014-11-30',
      west_bound  => -72.84215,
      east_bound  => -72.55687,
      south_bound => 43.30715,
      north_bound => 43.50186
    },
    'impervious_calendar_brook_vt' => {
      begin_date  => '2014-11-30',
      end_date    => '2014-11-30',
      west_bound  => -72.15850,
      east_bound  => -71.94029,
      south_bound => 44.54922,
      north_bound => 44.70891
    },
    'impervious_east_creek_vt' => {
      begin_date  => '2014-11-30',
      end_date    => '2014-11-30',
      west_bound  => -73.07816,
      east_bound  => -72.76167,
      south_bound => 43.56533,
      north_bound => 43.80843
    },
    'impervious_jay_branch_vt' => {
      begin_date  => '2014-11-30',
      end_date    => '2014-11-30',
      west_bound  => -72.59749,
      east_bound  => -72.25492,
      south_bound => 44.81928,
      north_bound => 45.06633
    },
    'impervious_jewett_brook_vt' => {
      begin_date  => '2014-11-30',
      end_date    => '2014-11-30',
      west_bound  => -73.20868,
      east_bound  => -73.00478,
      south_bound => 44.76437,
      north_bound => 44.91525
    },
    'impervious_lulls_brook_vt' => {
      begin_date  => '2014-11-30',
      end_date    => '2014-11-30',
      west_bound  => -72.59237,
      east_bound  => -72.23276,
      south_bound => 43.45014,
      north_bound => 43.69695
    },
    'impervious_mallets_bay10_vt' => {
      begin_date  => '2014-11-30',
      end_date    => '2014-11-30',
      west_bound  => -73.37203,
      east_bound  => -72.99797,
      south_bound => 44.44312,
      north_bound => 44.70178
    },
    'impervious_mallets_bay_vt' => {
      begin_date  => '2014-11-30',
      end_date    => '2014-11-30',
      west_bound  => -73.36774,
      east_bound  => -73.02573,
      south_bound => 44.44670,
      north_bound => 44.68213
    },
    'impervious_mallets_creek_vt' => {
      begin_date  => '2014-11-30',
      end_date    => '2014-11-30',
      west_bound  => -73.21811,
      east_bound  => -73.01123,
      south_bound => 44.50673,
      north_bound => 44.67909
    },
    'impervious_mill_brook_vt' => {
      begin_date  => '2014-11-30',
      end_date    => '2014-11-30',
      west_bound  => -72.97613,
      east_bound  => -72.73228,
      south_bound => 44.08476,
      north_bound => 44.26059
    },
    'impervious_mill_river_vt' => {
      begin_date  => '2014-11-30',
      end_date    => '2014-11-30',
      west_bound  => -73.20595,
      east_bound  => -73.00983,
      south_bound => 44.68966,
      north_bound => 44.83256
    },
    'impervious_moon_brook_vt' => {
      begin_date  => '2014-11-30',
      end_date    => '2014-11-30',
      west_bound  => -73.12229,
      east_bound  => -72.84459,
      south_bound => 43.47189,
      north_bound => 43.72002
    },
    'impervious_muddy_brook_vt' => {
      begin_date  => '2014-11-30',
      end_date    => '2014-11-30',
      west_bound  => -73.22225,
      east_bound  => -73.04496,
      south_bound => 44.34886,
      north_bound => 44.49121
    },
    'impervious_north_branch_deerfield_river_vt' => {
      begin_date  => '2014-11-30',
      end_date    => '2014-11-30',
      west_bound  => -72.98122,
      east_bound  => -72.72832,
      south_bound => 42.80714,
      north_bound => 43.02503
    },
    'impervious_ottauquechee_river_headwaters_vt' => {
      begin_date  => '2014-11-30',
      end_date    => '2014-11-30',
      west_bound  => -72.89036,
      east_bound  => -72.60851,
      south_bound => 43.52262,
      north_bound => 43.74325
    },
    'impervious_ottauquechee_river_vt' => {
      begin_date  => '2014-11-30',
      end_date    => '2014-11-30',
      west_bound  => -72.98809,
      east_bound  => -72.20105,
      south_bound => 43.34998,
      north_bound => 43.87024
    },
    'impervious_shelburne_river_vt' => {
      begin_date  => '2014-11-30',
      end_date    => '2014-11-30',
      west_bound  => -73.32527,
      east_bound  => -72.98207,
      south_bound => 44.24329,
      north_bound => 44.50658
    },
    'impervious_snipe_island_vt' => {
      begin_date  => '2014-11-30',
      end_date    => '2014-11-30',
      west_bound  => -73.12803,
      east_bound  => -72.78950,
      south_bound => 44.31307,
      north_bound => 44.54340
    },
    'impervious_winooski_river_vt' => {
      begin_date  => '2014-11-30',
      end_date    => '2014-11-30',
      west_bound  => -73.34566,
      east_bound  => -72.97283,
      south_bound => 44.34826,
      north_bound => 44.63423
    },

    # weekly_cyanobacteria are not prefixed by domain so list them all here:

    'meris_weekly_cyanobacteria_al' => {
      begin_date  => '2017-01-01',
      end_date    => '2018-08-12',
      west_bound  =>  -88.4771,
      east_bound  =>  -84.8907,
      south_bound =>   30.2299,
      north_bound =>   35.0186
    },
    'meris_weekly_cyanobacteria_ar' => {
      begin_date  => '2017-01-01',
      end_date    => '2018-08-12',
      west_bound  =>  -94.6216,
      east_bound  =>  -89.6405,
      south_bound =>   33.0064,
      north_bound =>   36.4972
    },
    'meris_weekly_cyanobacteria_az' => {
      begin_date  => '2017-01-01',
      end_date    => '2018-08-12',
      west_bound  => -114.8260,
      east_bound  => -109.0428,
      south_bound =>   31.3319,
      north_bound =>   37.0076
    },
    'meris_weekly_cyanobacteria_ca' => {
      begin_date  => '2002-05-21',
      end_date    => '2018-06-03',
      west_bound  => -124.3972,
      east_bound  => -114.1215,
      south_bound =>   32.5321,
      north_bound =>   42.0072
    },
    'meris_weekly_cyanobacteria_co' => {
      begin_date  => '2017-01-01',
      end_date    => '2018-08-12',
      west_bound  => -109.0603,
      east_bound  => -102.0319,
      south_bound =>   36.9853,
      north_bound =>   41.0058
    },
    'meris_weekly_cyanobacteria_ct' => {
      begin_date  => '2017-01-01',
      end_date    => '2018-08-12',
      west_bound  =>  -73.7298,
      east_bound  =>  -71.7821,
      south_bound =>   40.9947,
      north_bound =>   42.0523
    },
    'meris_weekly_cyanobacteria_de' => {
      begin_date  => '2017-01-01',
      end_date    => '2018-08-12',
      west_bound  =>  -75.7956,
      east_bound  =>  -75.0429,
      south_bound =>   38.4459,
      north_bound =>   39.8437
    },
    'meris_weekly_cyanobacteria_fl' => {
      begin_date  => '2002-05-21',
      end_date    => '2018-08-12',
      west_bound  =>  -87.6298,
      east_bound  =>  -80.0482,
      south_bound =>   24.9527,
      north_bound =>   31.0073
    },
    'meris_weekly_cyanobacteria_ga' => {
      begin_date  => '2017-01-01',
      end_date    => '2018-08-12',
      west_bound  =>  -85.6132,
      east_bound  =>  -80.8916,
      south_bound =>   30.3576,
      north_bound =>   35.0056
    },
    'meris_weekly_cyanobacteria_ia' => {
      begin_date  => '2017-01-01',
      end_date    => '2018-08-12',
      west_bound  =>  -96.6454,
      east_bound  =>  -90.1393,
      south_bound =>   40.3682,
      north_bound =>   43.5037
    },
    'meris_weekly_cyanobacteria_id' => {
      begin_date  => '2017-01-01',
      end_date    => '2018-08-12',
      west_bound  => -117.2417,
      east_bound  => -111.0430,
      south_bound =>   41.9909,
      north_bound =>   49.0031
    },
    'meris_weekly_cyanobacteria_il' => {
      begin_date  => '2017-01-01',
      end_date    => '2018-08-12',
      west_bound  =>  -91.5208,
      east_bound  =>  -87.5024,
      south_bound =>   36.9831,
      north_bound =>   42.5139
    },
    'meris_weekly_cyanobacteria_in' => {
      begin_date  => '2017-01-01',
      end_date    => '2018-08-12',
      west_bound  =>  -88.1059,
      east_bound  =>  -84.7832,
      south_bound =>   37.7725,
      north_bound =>   41.7678
    },
    'meris_weekly_cyanobacteria_ks' => {
      begin_date  => '2017-01-01',
      end_date    => '2018-08-12',
      west_bound  => -102.0560,
      east_bound  =>  -94.5968,
      south_bound =>   36.9852,
      north_bound =>   40.0082
    },
    'meris_weekly_cyanobacteria_ky' => {
      begin_date  => '2017-01-01',
      end_date    => '2018-08-12',
      west_bound  =>  -89.5726,
      east_bound  =>  -81.9538,
      south_bound =>   36.4929,
      north_bound =>   39.1446
    },
    'meris_weekly_cyanobacteria_la' => {
      begin_date  => '2017-01-01',
      end_date    => '2018-08-12',
      west_bound  =>  -94.0460,
      east_bound  =>  -89.0192,
      south_bound =>   28.9359,
      north_bound =>   33.0259
    },
    'meris_weekly_cyanobacteria_ma' => {
      begin_date  => '2017-01-01',
      end_date    => '2018-08-12',
      west_bound  =>  -73.5035,
      east_bound  =>  -69.9141,
      south_bound =>   41.2346,
      north_bound =>   42.8887
    },
    'meris_weekly_cyanobacteria_md' => {
      begin_date  => '2017-01-01',
      end_date    => '2018-08-12',
      west_bound  =>  -79.4943,
      east_bound  =>  -75.0414,
      south_bound =>   37.9666,
      north_bound =>   39.7285
    },
    'meris_weekly_cyanobacteria_me' => {
      begin_date  => '2017-01-01',
      end_date    => '2018-08-12',
      west_bound  =>  -71.0923,
      east_bound  =>  -66.9648,
      south_bound =>   43.0874,
      north_bound =>   47.4576
    },
    'meris_weekly_cyanobacteria_mi' => {
      begin_date  => '2017-01-01',
      end_date    => '2018-08-12',
      west_bound  =>  -90.4131,
      east_bound  =>  -82.4151,
      south_bound =>   41.6938,
      north_bound =>   48.1761
    },
    'meris_weekly_cyanobacteria_mn' => {
      begin_date  => '2017-01-01',
      end_date    => '2018-08-12',
      west_bound  =>  -97.2345,
      east_bound  =>  -89.5252,
      south_bound =>   43.4944,
      north_bound =>   49.3741
    },
    'meris_weekly_cyanobacteria_mo' => {
      begin_date  => '2017-01-01',
      end_date    => '2018-08-12',
      west_bound  =>  -95.7720,
      east_bound  =>  -89.1003,
      south_bound =>   35.9860,
      north_bound =>   40.6142
    },
    'meris_weekly_cyanobacteria_ms' => {
      begin_date  => '2017-01-01',
      end_date    => '2018-08-12',
      west_bound  =>  -91.6478,
      east_bound  =>  -88.0871,
      south_bound =>   30.1912,
      north_bound =>   35.0081
    },
    'meris_weekly_cyanobacteria_mt' => {
      begin_date  => '2017-01-01',
      end_date    => '2018-08-12',
      west_bound  => -116.0686,
      east_bound  => -104.0401,
      south_bound =>   44.3499,
      north_bound =>   49.0055
    },
    'meris_weekly_cyanobacteria_nc' => {
      begin_date  => '2017-01-01',
      end_date    => '2017-12-31',
      west_bound  =>  -84.3281,
      east_bound  =>  -75.4515,
      south_bound =>   33.8785,
      north_bound =>   36.5936
    },
    'meris_weekly_cyanobacteria_nd' => {
      begin_date  => '2017-01-01',
      end_date    => '2018-08-12',
      west_bound  => -104.0671,
      east_bound  =>  -96.5475,
      south_bound =>   45.9282,
      north_bound =>   48.9961
    },
    'meris_weekly_cyanobacteria_ne' => {
      begin_date  => '2017-01-01',
      end_date    => '2018-08-12',
      west_bound  => -104.0609,
      east_bound  =>  -95.3058,
      south_bound =>   39.9889,
      north_bound =>   43.0079
    },
    'meris_weekly_cyanobacteria_nh' => {
      begin_date  => '2017-01-01',
      end_date    => '2018-08-12',
      west_bound  =>  -72.5581,
      east_bound  =>  -70.7279,
      south_bound =>   42.6949,
      north_bound =>   45.3042
    },
    'meris_weekly_cyanobacteria_nj' => {
      begin_date  => '2017-01-01',
      end_date    => '2018-08-12',
      west_bound  =>  -75.5747,
      east_bound  =>  -73.8900,
      south_bound =>   38.9530,
      north_bound =>   41.3529
    },
    'meris_weekly_cyanobacteria_nm' => {
      begin_date  => '2017-01-01',
      end_date    => '2018-08-12',
      west_bound  => -109.0555,
      east_bound  => -102.9932,
      south_bound =>   31.3397,
      north_bound =>   37.0047
    },
    'meris_weekly_cyanobacteria_nv' => {
      begin_date  => '2017-01-01',
      end_date    => '2018-08-12',
      west_bound  => -120.0009,
      east_bound  => -114.0328,
      south_bound =>   34.9952,
      north_bound =>   42.0017
    },
    'meris_weekly_cyanobacteria_ny' => {
      begin_date  => '2017-01-01',
      end_date    => '2018-06-03',
      west_bound  =>  -79.7679,
      east_bound  =>  -71.8668,
      south_bound =>   40.5023,
      north_bound =>   45.0088
    },
    'meris_weekly_cyanobacteria_oh' => {
      begin_date  => '2017-01-01',
      end_date    => '2018-08-12',
      west_bound  =>  -84.8165,
      east_bound  =>  -80.5147,
      south_bound =>   38.3968,
      north_bound =>   41.9913
    },
    'meris_weekly_cyanobacteria_ok' => {
      begin_date  => '2017-01-01',
      end_date    => '2018-08-12',
      west_bound  => -103.0021,
      east_bound  =>  -94.4258,
      south_bound =>   33.6174,
      north_bound =>   37.0038
    },
    'meris_weekly_cyanobacteria_or' => {
      begin_date  => '2017-01-01',
      end_date    => '2018-06-03',
      west_bound  => -124.5643,
      east_bound  => -116.4661,
      south_bound =>   41.9840,
      north_bound =>   46.2385
    },
    'meris_weekly_cyanobacteria_pa' => {
      begin_date  => '2017-01-01',
      end_date    => '2018-08-12',
      west_bound  =>  -80.5306,
      east_bound  =>  -74.6965,
      south_bound =>   39.7156,
      north_bound =>   42.2701
    },
    'meris_weekly_cyanobacteria_ri' => {
      begin_date  => '2017-01-01',
      end_date    => '2018-08-12',
      west_bound  =>  -71.8713,
      east_bound  =>  -71.1120,
      south_bound =>   41.3191,
      north_bound =>   42.0171
    },
    'meris_weekly_cyanobacteria_sc' => {
      begin_date  => '2017-01-01',
      end_date    => '2018-08-12',
      west_bound  => -83.354,
      east_bound  => -78.576,
      south_bound =>  32.065,
      north_bound =>  35.2112
    },
    'meris_weekly_cyanobacteria_sd' => {
      begin_date  => '2017-01-01',
      end_date    => '2018-08-12',
      west_bound  => -104.0658,
      east_bound  =>  -96.4349,
      south_bound =>   42.4848,
      north_bound =>   45.9474
    },
    'meris_weekly_cyanobacteria_tn' => {
      begin_date  => '2017-01-01',
      end_date    => '2018-08-12',
      west_bound  =>  -90.3097,
      east_bound  =>  -81.6471,
      south_bound =>   34.9851,
      north_bound =>   36.6826
    },
    'meris_weekly_cyanobacteria_tx' => {
      begin_date  => '2017-01-01',
      end_date    => '2018-08-12',
      west_bound  => -106.6542,
      east_bound  =>  -93.5026,
      south_bound =>   25.8418,
      north_bound =>   36.4989
    },
    'meris_weekly_cyanobacteria_ut' => {
      begin_date  => '2017-01-01',
      end_date    => '2018-08-12',
      west_bound  => -114.0517,
      east_bound  => -109.0391,
      south_bound =>   36.9880,
      north_bound =>   42.0069
    },
    'meris_weekly_cyanobacteria_va' => {
      begin_date  => '2017-01-01',
      end_date    => '2018-08-12',
      west_bound  =>  -83.6795,
      east_bound  =>  -75.2371,
      south_bound =>   36.5379,
      north_bound =>   39.4596
    },
    'meris_weekly_cyanobacteria_vt' => {
      begin_date  => '2017-01-01',
      end_date    => '2018-08-12',
      west_bound  => -73.4398,
      east_bound  => -71.5,
      south_bound =>  42.7231,
      north_bound =>  45.016
    },
    'meris_weekly_cyanobacteria_wa' => {
      begin_date  => '2017-01-01',
      end_date    => '2018-08-12',
      west_bound  => -124.7378,
      east_bound  => -116.9137,
      south_bound =>   45.5394,
      north_bound =>   49.0025
    },
    'meris_weekly_cyanobacteria_wi' => {
      begin_date  => '2017-01-01',
      end_date    => '2018-08-12',
      west_bound  =>  -92.8902,
      east_bound  =>  -86.9641,
      south_bound =>   42.4855,
      north_bound =>   46.9553
    },
    'meris_weekly_cyanobacteria_wv' => {
      begin_date  => '2017-01-01',
      end_date    => '2018-08-12',
      west_bound  =>  -82.6516,
      east_bound  =>  -77.7239,
      south_bound =>   37.2012,
      north_bound =>   40.6410
    },
    'meris_weekly_cyanobacteria_wy' => {
      begin_date  => '2017-01-01',
      end_date    => '2018-08-12',
      west_bound  => -111.0582,
      east_bound  => -104.0463,
      south_bound =>   40.9906,
      north_bound =>   45.0080
    },

    'sentinel_weekly_cyanobacteria_al' => {
      begin_date  => '2017-01-01',
      end_date    => '2018-08-12',
      west_bound  =>  -88.4771,
      east_bound  =>  -84.8907,
      south_bound =>   30.2299,
      north_bound =>   35.0186
    },
    'sentinel_weekly_cyanobacteria_ar' => {
      begin_date  => '2017-01-01',
      end_date    => '2018-08-12',
      west_bound  =>  -94.6216,
      east_bound  =>  -89.6405,
      south_bound =>   33.0064,
      north_bound =>   36.4972
    },
    'sentinel_weekly_cyanobacteria_az' => {
      begin_date  => '2017-01-01',
      end_date    => '2018-08-12',
      west_bound  => -114.8260,
      east_bound  => -109.0428,
      south_bound =>   31.3319,
      north_bound =>   37.0076
    },
    'sentinel_weekly_cyanobacteria_ca' => {
      begin_date  => '2002-05-21',
      end_date    => '2018-06-03',
      west_bound  => -124.3972,
      east_bound  => -114.1215,
      south_bound =>   32.5321,
      north_bound =>   42.0072
    },
    'sentinel_weekly_cyanobacteria_co' => {
      begin_date  => '2017-01-01',
      end_date    => '2018-08-12',
      west_bound  => -109.0603,
      east_bound  => -102.0319,
      south_bound =>   36.9853,
      north_bound =>   41.0058
    },
    'sentinel_weekly_cyanobacteria_ct' => {
      begin_date  => '2017-01-01',
      end_date    => '2018-08-12',
      west_bound  =>  -73.7298,
      east_bound  =>  -71.7821,
      south_bound =>   40.9947,
      north_bound =>   42.0523
    },
    'sentinel_weekly_cyanobacteria_de' => {
      begin_date  => '2017-01-01',
      end_date    => '2018-08-12',
      west_bound  =>  -75.7956,
      east_bound  =>  -75.0429,
      south_bound =>   38.4459,
      north_bound =>   39.8437
    },
    'sentinel_weekly_cyanobacteria_fl' => {
      begin_date  => '2002-05-21',
      end_date    => '2018-08-12',
      west_bound  =>  -87.6298,
      east_bound  =>  -80.0482,
      south_bound =>   24.9527,
      north_bound =>   31.0073
    },
    'sentinel_weekly_cyanobacteria_ga' => {
      begin_date  => '2017-01-01',
      end_date    => '2018-08-12',
      west_bound  =>  -85.6132,
      east_bound  =>  -80.8916,
      south_bound =>   30.3576,
      north_bound =>   35.0056
    },
    'sentinel_weekly_cyanobacteria_ia' => {
      begin_date  => '2017-01-01',
      end_date    => '2018-08-12',
      west_bound  =>  -96.6454,
      east_bound  =>  -90.1393,
      south_bound =>   40.3682,
      north_bound =>   43.5037
    },
    'sentinel_weekly_cyanobacteria_id' => {
      begin_date  => '2017-01-01',
      end_date    => '2018-08-12',
      west_bound  => -117.2417,
      east_bound  => -111.0430,
      south_bound =>   41.9909,
      north_bound =>   49.0031
    },
    'sentinel_weekly_cyanobacteria_il' => {
      begin_date  => '2017-01-01',
      end_date    => '2018-08-12',
      west_bound  =>  -91.5208,
      east_bound  =>  -87.5024,
      south_bound =>   36.9831,
      north_bound =>   42.5139
    },
    'sentinel_weekly_cyanobacteria_in' => {
      begin_date  => '2017-01-01',
      end_date    => '2018-08-12',
      west_bound  =>  -88.1059,
      east_bound  =>  -84.7832,
      south_bound =>   37.7725,
      north_bound =>   41.7678
    },
    'sentinel_weekly_cyanobacteria_ks' => {
      begin_date  => '2017-01-01',
      end_date    => '2018-08-12',
      west_bound  => -102.0560,
      east_bound  =>  -94.5968,
      south_bound =>   36.9852,
      north_bound =>   40.0082
    },
    'sentinel_weekly_cyanobacteria_ky' => {
      begin_date  => '2017-01-01',
      end_date    => '2018-08-12',
      west_bound  =>  -89.5726,
      east_bound  =>  -81.9538,
      south_bound =>   36.4929,
      north_bound =>   39.1446
    },
    'sentinel_weekly_cyanobacteria_la' => {
      begin_date  => '2017-01-01',
      end_date    => '2018-08-12',
      west_bound  =>  -94.0460,
      east_bound  =>  -89.0192,
      south_bound =>   28.9359,
      north_bound =>   33.0259
    },
    'sentinel_weekly_cyanobacteria_ma' => {
      begin_date  => '2017-01-01',
      end_date    => '2018-08-12',
      west_bound  =>  -73.5035,
      east_bound  =>  -69.9141,
      south_bound =>   41.2346,
      north_bound =>   42.8887
    },
    'sentinel_weekly_cyanobacteria_md' => {
      begin_date  => '2017-01-01',
      end_date    => '2018-08-12',
      west_bound  =>  -79.4943,
      east_bound  =>  -75.0414,
      south_bound =>   37.9666,
      north_bound =>   39.7285
    },
    'sentinel_weekly_cyanobacteria_me' => {
      begin_date  => '2017-01-01',
      end_date    => '2018-08-12',
      west_bound  =>  -71.0923,
      east_bound  =>  -66.9648,
      south_bound =>   43.0874,
      north_bound =>   47.4576
    },
    'sentinel_weekly_cyanobacteria_mi' => {
      begin_date  => '2017-01-01',
      end_date    => '2018-08-12',
      west_bound  =>  -90.4131,
      east_bound  =>  -82.4151,
      south_bound =>   41.6938,
      north_bound =>   48.1761
    },
    'sentinel_weekly_cyanobacteria_mn' => {
      begin_date  => '2017-01-01',
      end_date    => '2018-08-12',
      west_bound  =>  -97.2345,
      east_bound  =>  -89.5252,
      south_bound =>   43.4944,
      north_bound =>   49.3741
    },
    'sentinel_weekly_cyanobacteria_mo' => {
      begin_date  => '2017-01-01',
      end_date    => '2018-08-12',
      west_bound  =>  -95.7720,
      east_bound  =>  -89.1003,
      south_bound =>   35.9860,
      north_bound =>   40.6142
    },
    'sentinel_weekly_cyanobacteria_ms' => {
      begin_date  => '2017-01-01',
      end_date    => '2018-08-12',
      west_bound  =>  -91.6478,
      east_bound  =>  -88.0871,
      south_bound =>   30.1912,
      north_bound =>   35.0081
    },
    'sentinel_weekly_cyanobacteria_mt' => {
      begin_date  => '2017-01-01',
      end_date    => '2018-08-12',
      west_bound  => -116.0686,
      east_bound  => -104.0401,
      south_bound =>   44.3499,
      north_bound =>   49.0055
    },
    'sentinel_weekly_cyanobacteria_nc' => {
      begin_date  => '2017-01-01',
      end_date    => '2017-12-31',
      west_bound  =>  -84.3281,
      east_bound  =>  -75.4515,
      south_bound =>   33.8785,
      north_bound =>   36.5936
    },
    'sentinel_weekly_cyanobacteria_nd' => {
      begin_date  => '2017-01-01',
      end_date    => '2018-08-12',
      west_bound  => -104.0671,
      east_bound  =>  -96.5475,
      south_bound =>   45.9282,
      north_bound =>   48.9961
    },
    'sentinel_weekly_cyanobacteria_ne' => {
      begin_date  => '2017-01-01',
      end_date    => '2018-08-12',
      west_bound  => -104.0609,
      east_bound  =>  -95.3058,
      south_bound =>   39.9889,
      north_bound =>   43.0079
    },
    'sentinel_weekly_cyanobacteria_nh' => {
      begin_date  => '2017-01-01',
      end_date    => '2018-08-12',
      west_bound  =>  -72.5581,
      east_bound  =>  -70.7279,
      south_bound =>   42.6949,
      north_bound =>   45.3042
    },
    'sentinel_weekly_cyanobacteria_nj' => {
      begin_date  => '2017-01-01',
      end_date    => '2018-08-12',
      west_bound  =>  -75.5747,
      east_bound  =>  -73.8900,
      south_bound =>   38.9530,
      north_bound =>   41.3529
    },
    'sentinel_weekly_cyanobacteria_nm' => {
      begin_date  => '2017-01-01',
      end_date    => '2018-08-12',
      west_bound  => -109.0555,
      east_bound  => -102.9932,
      south_bound =>   31.3397,
      north_bound =>   37.0047
    },
    'sentinel_weekly_cyanobacteria_nv' => {
      begin_date  => '2017-01-01',
      end_date    => '2018-08-12',
      west_bound  => -120.0009,
      east_bound  => -114.0328,
      south_bound =>   34.9952,
      north_bound =>   42.0017
    },
    'sentinel_weekly_cyanobacteria_ny' => {
      begin_date  => '2017-01-01',
      end_date    => '2018-06-03',
      west_bound  =>  -79.7679,
      east_bound  =>  -71.8668,
      south_bound =>   40.5023,
      north_bound =>   45.0088
    },
    'sentinel_weekly_cyanobacteria_oh' => {
      begin_date  => '2017-01-01',
      end_date    => '2018-08-12',
      west_bound  =>  -84.8165,
      east_bound  =>  -80.5147,
      south_bound =>   38.3968,
      north_bound =>   41.9913
    },
    'sentinel_weekly_cyanobacteria_ok' => {
      begin_date  => '2017-01-01',
      end_date    => '2018-08-12',
      west_bound  => -103.0021,
      east_bound  =>  -94.4258,
      south_bound =>   33.6174,
      north_bound =>   37.0038
    },
    'sentinel_weekly_cyanobacteria_or' => {
      begin_date  => '2017-01-01',
      end_date    => '2018-06-03',
      west_bound  => -124.5643,
      east_bound  => -116.4661,
      south_bound =>   41.9840,
      north_bound =>   46.2385
    },
    'sentinel_weekly_cyanobacteria_pa' => {
      begin_date  => '2017-01-01',
      end_date    => '2018-08-12',
      west_bound  =>  -80.5306,
      east_bound  =>  -74.6965,
      south_bound =>   39.7156,
      north_bound =>   42.2701
    },
    'sentinel_weekly_cyanobacteria_ri' => {
      begin_date  => '2017-01-01',
      end_date    => '2018-08-12',
      west_bound  =>  -71.8713,
      east_bound  =>  -71.1120,
      south_bound =>   41.3191,
      north_bound =>   42.0171
    },
    'sentinel_weekly_cyanobacteria_sc' => {
      begin_date  => '2017-01-01',
      end_date    => '2018-08-12',
      west_bound  => -83.354,
      east_bound  => -78.576,
      south_bound =>  32.065,
      north_bound =>  35.2112
    },
    'sentinel_weekly_cyanobacteria_sd' => {
      begin_date  => '2017-01-01',
      end_date    => '2018-08-12',
      west_bound  => -104.0658,
      east_bound  =>  -96.4349,
      south_bound =>   42.4848,
      north_bound =>   45.9474
    },
    'sentinel_weekly_cyanobacteria_tn' => {
      begin_date  => '2017-01-01',
      end_date    => '2018-08-12',
      west_bound  =>  -90.3097,
      east_bound  =>  -81.6471,
      south_bound =>   34.9851,
      north_bound =>   36.6826
    },
    'sentinel_weekly_cyanobacteria_tx' => {
      begin_date  => '2017-01-01',
      end_date    => '2018-08-12',
      west_bound  => -106.6542,
      east_bound  =>  -93.5026,
      south_bound =>   25.8418,
      north_bound =>   36.4989
    },
    'sentinel_weekly_cyanobacteria_ut' => {
      begin_date  => '2017-01-01',
      end_date    => '2018-08-12',
      west_bound  => -114.0517,
      east_bound  => -109.0391,
      south_bound =>   36.9880,
      north_bound =>   42.0069
    },
    'sentinel_weekly_cyanobacteria_va' => {
      begin_date  => '2017-01-01',
      end_date    => '2018-08-12',
      west_bound  =>  -83.6795,
      east_bound  =>  -75.2371,
      south_bound =>   36.5379,
      north_bound =>   39.4596
    },
    'sentinel_weekly_cyanobacteria_vt' => {
      begin_date  => '2017-01-01',
      end_date    => '2018-08-12',
      west_bound  => -73.4398,
      east_bound  => -71.5,
      south_bound =>  42.7231,
      north_bound =>  45.016
    },
    'sentinel_weekly_cyanobacteria_wa' => {
      begin_date  => '2017-01-01',
      end_date    => '2018-08-12',
      west_bound  => -124.7378,
      east_bound  => -116.9137,
      south_bound =>   45.5394,
      north_bound =>   49.0025
    },
    'sentinel_weekly_cyanobacteria_wi' => {
      begin_date  => '2017-01-01',
      end_date    => '2018-08-12',
      west_bound  =>  -92.8902,
      east_bound  =>  -86.9641,
      south_bound =>   42.4855,
      north_bound =>   46.9553
    },
    'sentinel_weekly_cyanobacteria_wv' => {
      begin_date  => '2017-01-01',
      end_date    => '2018-08-12',
      west_bound  =>  -82.6516,
      east_bound  =>  -77.7239,
      south_bound =>   37.2012,
      north_bound =>   40.6410
    },
    'sentinel_weekly_cyanobacteria_wy' => {
      begin_date  => '2017-01-01',
      end_date    => '2018-08-12',
      west_bound  => -111.0582,
      east_bound  => -104.0463,
      south_bound =>   40.9906,
      north_bound =>   45.0080
    },

    'land_change_2020_current_chesapeake' => {
      begin_date  => '2020-01-01',
      end_date    => '2020-12-31',
      west_bound  => -79.92,
      east_bound  => -75.048,
      south_bound => 37.116,
      north_bound => 40.344
    },
    'land_change_2030_current_chesapeake' => {
      begin_date  => '2030-01-01',
      end_date    => '2030-12-31',
      west_bound  => -79.92,
      east_bound  => -75.048,
      south_bound => 37.116,
      north_bound => 40.344
    },
    'land_change_2040_current_chesapeake' => {
      begin_date  => '2040-01-01',
      end_date    => '2040-12-31',
      west_bound  => -79.92,
      east_bound  => -75.048,
      south_bound => 37.116,
      north_bound => 40.344
    },
    'land_change_2020_deregulation_chesapeake' => {
      begin_date  => '2020-01-01',
      end_date    => '2020-12-31',
      west_bound  => -79.92,
      east_bound  => -75.048,
      south_bound => 37.116,
      north_bound => 40.344
    },
    'land_change_2030_deregulation_chesapeake' => {
      begin_date  => '2030-01-01',
      end_date    => '2030-12-31',
      west_bound  => -79.92,
      east_bound  => -75.048,
      south_bound => 37.116,
      north_bound => 40.344
    },
    'land_change_2040_deregulation_chesapeake' => {
      begin_date  => '2040-01-01',
      end_date    => '2040-12-31',
      west_bound  => -79.92,
      east_bound  => -75.048,
      south_bound => 37.116,
      north_bound => 40.344
    },
    'land_change_2020_redevelopment_chesapeake' => {
      begin_date  => '2020-01-01',
      end_date    => '2020-12-31',
      west_bound  => -79.92,
      east_bound  => -75.048,
      south_bound => 37.116,
      north_bound => 40.344
    },
    'land_change_2030_redevelopment_chesapeake' => {
      begin_date  => '2030-01-01',
      end_date    => '2030-12-31',
      west_bound  => -79.92,
      east_bound  => -75.048,
      south_bound => 37.116,
      north_bound => 40.344
    },
    'land_change_2040_redevelopment_chesapeake' => {
      begin_date  => '2040-01-01',
      end_date    => '2040-12-31',
      west_bound  => -79.92,
      east_bound  => -75.048,
      south_bound => 37.116,
      north_bound => 40.344
    },
    'land_change_2020_conservation_chesapeake' => {
      begin_date  => '2020-01-01',
      end_date    => '2020-12-31',
      west_bound  => -79.92,
      east_bound  => -75.048,
      south_bound => 37.116,
      north_bound => 40.344
    },
    'land_change_2030_conservation_chesapeake' => {
      begin_date  => '2030-01-01',
      end_date    => '2030-12-31',
      west_bound  => -79.92,
      east_bound  => -75.048,
      south_bound => 37.116,
      north_bound => 40.344
    },
    'land_change_2040_conservation_chesapeake' => {
      begin_date  => '2040-01-01',
      end_date    => '2040-12-31',
      west_bound  => -79.92,
      east_bound  => -75.048,
      south_bound => 37.116,
      north_bound => 40.344
    },
    'stream_temperature_point_median_07_new_england' => {
      begin_date  => '1995-01-01',
      end_date    => '2010-12-31',
      west_bound  => -73.511939,
      east_bound  => -67.088797,
      south_bound => 41.169369,
      north_bound => 46.144257
    },
    'stream_temperature_point_median_08_new_england' => {
      begin_date  => '1995-01-01',
      end_date    => '2010-12-31',
      west_bound  => -73.511939,
      east_bound  => -67.088797,
      south_bound => 41.169369,
      north_bound => 46.144257
    },
    'stream_temperature_point_mean_daily_range_07_new_england' => {
      begin_date  => '1995-01-01',
      end_date    => '2010-12-31',
      west_bound  => -73.511939,
      east_bound  => -67.088797,
      south_bound => 41.169369,
      north_bound => 46.144257
    },
    'stream_temperature_point_mean_daily_range_08_new_england' => {
      begin_date  => '1995-01-01',
      end_date    => '2010-12-31',
      west_bound  => -73.511939,
      east_bound  => -67.088797,
      south_bound => 41.169369,
      north_bound => 46.144257
    },
    'stream_temperature_point_maximum_daily_increase_07_new_england' => {
      begin_date  => '1995-01-01',
      end_date    => '2010-12-31',
      west_bound  => -73.511939,
      east_bound  => -67.088797,
      south_bound => 41.169369,
      north_bound => 46.144257
    },
    'stream_temperature_point_maximum_daily_increase_08_new_england' => {
      begin_date  => '1995-01-01',
      end_date    => '2010-12-31',
      west_bound  => -73.511939,
      east_bound  => -67.088797,
      south_bound => 41.169369,
      north_bound => 46.144257
    },
    'stream_temperature_point_maximum_daily_decrease_07_new_england' => {
      begin_date  => '1995-01-01',
      end_date    => '2010-12-31',
      west_bound  => -73.511939,
      east_bound  => -67.088797,
      south_bound => 41.169369,
      north_bound => 46.144257
    },
    'stream_temperature_point_maximum_daily_decrease_08_new_england' => {
      begin_date  => '1995-01-01',
      end_date    => '2010-12-31',
      west_bound  => -73.511939,
      east_bound  => -67.088797,
      south_bound => 41.169369,
      north_bound => 46.144257
    },
    'stream_temperature_point_maximum_new_england' => {
      begin_date  => '1995-01-01',
      end_date    => '2010-12-31',
      west_bound  => -73.511939,
      east_bound  => -67.088797,
      south_bound => 41.169369,
      north_bound => 46.144257
    },
    'stream_temperature_point_day_of_maximum_new_england' => {
      begin_date  => '1995-01-01',
      end_date    => '2010-12-31',
      west_bound  => -73.511939,
      east_bound  => -67.088797,
      south_bound => 41.169369,
      north_bound => 46.144257
    },
    'stream_temperature_line_07_new_england' => {
      begin_date  => '1995-01-01',
      end_date    => '2010-12-31',
      west_bound  => -73.713,
      east_bound  => -66.988,
      south_bound => 41.001,
      north_bound => 46.462
    },

    'stream_temperature_line_meduxnekeag_river' => {
      begin_date  => '2010-07-01',
      end_date    => '2011-08-31',
      west_bound  => -68.124,
      east_bound  => -67.572,
      south_bound => 45.911,
      north_bound => 46.452
    },

    # WSM:

    'hspf_blackstone10_huc10_point.hydrology' => {
      begin_date  => '1960-01-01',
      end_date    => '2009-12-31',
      west_bound  =>  -71.57708,
      east_bound  =>  -71.57708,
      south_bound =>   42.01962,
      north_bound =>   42.01962
    },
    'hspf_blackstone10_huc10_polygon.hydrology' => {
      begin_date  => '1960-01-01',
      end_date    => '2009-12-31',
      west_bound  =>  -71.78048,
      east_bound  =>  -71.35340,
      south_bound =>   41.85489,
      north_bound =>   42.23843
    },
    'hspf_blackstone10_huc12_point.hydrology' => {
      begin_date  => '1960-01-01',
      end_date    => '2009-12-31',
      west_bound  =>  -71.71881,
      east_bound  =>  -71.39789,
      south_bound =>   41.91013,
      north_bound =>   42.15551
    },
    'hspf_blackstone10_huc12_polygon.hydrology' => {
      begin_date  => '1960-01-01',
      end_date    => '2009-12-31',
      west_bound  =>  -71.78048,
      east_bound  =>  -71.35340,
      south_bound =>   41.85489,
      north_bound =>   42.23843
    },
    'hspf_blackstone1_huc10_point.hydrology' => {
      begin_date  => '1960-01-01',
      end_date    => '2009-12-31',
      west_bound  =>  -71.76534,
      east_bound  =>  -71.76534,
      south_bound =>   42.19291,
      north_bound =>   42.19291
    },
    'hspf_blackstone1_huc10_polygon.hydrology' => {
      begin_date  => '1960-01-01',
      end_date    => '2009-12-31',
      west_bound  =>  -71.92799,
      east_bound  =>  -71.60905,
      south_bound =>   42.03787,
      north_bound =>   42.35075
    },
    'hspf_blackstone1_huc12_point.hydrology' => {
      begin_date  => '1960-01-01',
      end_date    => '2009-12-31',
      west_bound  =>  -71.85803,
      east_bound  =>  -71.71627,
      south_bound =>   42.09478,
      north_bound =>   42.26947
    },
    'hspf_blackstone1_huc12_polygon.hydrology' => {
      begin_date  => '1960-01-01',
      end_date    => '2009-12-31',
      west_bound  =>  -71.92799,
      east_bound  =>  -71.60905,
      south_bound =>   42.03787,
      north_bound =>   42.35075
    },
    'hspf_blackstone2_huc10_point.hydrology' => {
      begin_date  => '1960-01-01',
      end_date    => '2009-12-31',
      west_bound  =>  -71.57708,
      east_bound  =>  -71.57708,
      south_bound =>   42.01962,
      north_bound =>   42.01962
    },
    'hspf_blackstone2_huc10_polygon.hydrology' => {
      begin_date  => '1960-01-01',
      end_date    => '2009-12-31',
      west_bound  =>  -71.78048,
      east_bound  =>  -71.35340,
      south_bound =>   41.85489,
      north_bound =>   42.23843
    },
    'hspf_blackstone2_huc12_point.hydrology' => {
      begin_date  => '1960-01-01',
      end_date    => '2009-12-31',
      west_bound  =>  -71.71881,
      east_bound  =>  -71.39789,
      south_bound =>   41.91013,
      north_bound =>   42.15551
    },
    'hspf_blackstone2_huc12_polygon.hydrology' => {
      begin_date  => '1960-01-01',
      end_date    => '2009-12-31',
      west_bound  =>  -71.78048,
      east_bound  =>  -71.35340,
      south_bound =>   41.85489,
      north_bound =>   42.23843
    },
    'hspf_blackstone34_huc10_point.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.76534,
      east_bound  =>  -71.76534,
      south_bound =>   42.19291,
      north_bound =>   42.19291
    },
    'hspf_blackstone34_huc10_polygon.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.92799,
      east_bound  =>  -71.60905,
      south_bound =>   42.03787,
      north_bound =>   42.35075
    },
    'hspf_blackstone34_huc12_point.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.85803,
      east_bound  =>  -71.71627,
      south_bound =>   42.09478,
      north_bound =>   42.26947
    },
    'hspf_blackstone34_huc12_polygon.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.92799,
      east_bound  =>  -71.60905,
      south_bound =>   42.03787,
      north_bound =>   42.35075
    },
    'hspf_blackstone35_huc10_point.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.76534,
      east_bound  =>  -71.76534,
      south_bound =>   42.19291,
      north_bound =>   42.19291
    },
    'hspf_blackstone35_huc10_polygon.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.92799,
      east_bound  =>  -71.60905,
      south_bound =>   42.03787,
      north_bound =>   42.35075
    },
    'hspf_blackstone35_huc12_point.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.85803,
      east_bound  =>  -71.71627,
      south_bound =>   42.09478,
      north_bound =>   42.26947
    },
    'hspf_blackstone35_huc12_polygon.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.92799,
      east_bound  =>  -71.60905,
      south_bound =>   42.03787,
      north_bound =>   42.35075
    },
    'hspf_blackstone36_huc10_point.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.76534,
      east_bound  =>  -71.76534,
      south_bound =>   42.19291,
      north_bound =>   42.19291
    },
    'hspf_blackstone36_huc10_polygon.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.92799,
      east_bound  =>  -71.60905,
      south_bound =>   42.03787,
      north_bound =>   42.35075
    },
    'hspf_blackstone36_huc12_point.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.85803,
      east_bound  =>  -71.71627,
      south_bound =>   42.09478,
      north_bound =>   42.26947
    },
    'hspf_blackstone36_huc12_polygon.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.92799,
      east_bound  =>  -71.60905,
      south_bound =>   42.03787,
      north_bound =>   42.35075
    },
    'hspf_blackstone37_huc10_point.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.76534,
      east_bound  =>  -71.76534,
      south_bound =>   42.19291,
      north_bound =>   42.19291
    },
    'hspf_blackstone37_huc10_polygon.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.92799,
      east_bound  =>  -71.60905,
      south_bound =>   42.03787,
      north_bound =>   42.35075
    },
    'hspf_blackstone37_huc12_point.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.85803,
      east_bound  =>  -71.71627,
      south_bound =>   42.09478,
      north_bound =>   42.26947
    },
    'hspf_blackstone37_huc12_polygon.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.92799,
      east_bound  =>  -71.60905,
      south_bound =>   42.03787,
      north_bound =>   42.35075
    },
    'hspf_blackstone38_huc10_point.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.76534,
      east_bound  =>  -71.76534,
      south_bound =>   42.19291,
      north_bound =>   42.19291
    },
    'hspf_blackstone38_huc10_polygon.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.92799,
      east_bound  =>  -71.60905,
      south_bound =>   42.03787,
      north_bound =>   42.35075
    },
    'hspf_blackstone38_huc12_point.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.85803,
      east_bound  =>  -71.71627,
      south_bound =>   42.09478,
      north_bound =>   42.26947
    },
    'hspf_blackstone38_huc12_polygon.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.92799,
      east_bound  =>  -71.60905,
      south_bound =>   42.03787,
      north_bound =>   42.35075
    },
    'hspf_blackstone39_huc10_point.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.57708,
      east_bound  =>  -71.57708,
      south_bound =>   42.01962,
      north_bound =>   42.01962
    },
    'hspf_blackstone39_huc10_polygon.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.78048,
      east_bound  =>  -71.35340,
      south_bound =>   41.85489,
      north_bound =>   42.23843
    },
    'hspf_blackstone39_huc12_point.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.71881,
      east_bound  =>  -71.39789,
      south_bound =>   41.91013,
      north_bound =>   42.15551
    },
    'hspf_blackstone39_huc12_polygon.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.78048,
      east_bound  =>  -71.35340,
      south_bound =>   41.85489,
      north_bound =>   42.23843
    },
    'hspf_blackstone40_huc10_point.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.57708,
      east_bound  =>  -71.57708,
      south_bound =>   42.01962,
      north_bound =>   42.01962
    },
    'hspf_blackstone40_huc10_polygon.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.78048,
      east_bound  =>  -71.35340,
      south_bound =>   41.85489,
      north_bound =>   42.23843
    },
    'hspf_blackstone40_huc12_point.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.71881,
      east_bound  =>  -71.39789,
      south_bound =>   41.91013,
      north_bound =>   42.15551
    },
    'hspf_blackstone40_huc12_polygon.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.78048,
      east_bound  =>  -71.35340,
      south_bound =>   41.85489,
      north_bound =>   42.23843
    },
    'hspf_blackstone41_huc10_point.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.57708,
      east_bound  =>  -71.57708,
      south_bound =>   42.01962,
      north_bound =>   42.01962
    },
    'hspf_blackstone41_huc10_polygon.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.78048,
      east_bound  =>  -71.35340,
      south_bound =>   41.85489,
      north_bound =>   42.23843
    },
    'hspf_blackstone41_huc12_point.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.71881,
      east_bound  =>  -71.39789,
      south_bound =>   41.91013,
      north_bound =>   42.15551
    },
    'hspf_blackstone41_huc12_polygon.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.78048,
      east_bound  =>  -71.35340,
      south_bound =>   41.85489,
      north_bound =>   42.23843
    },
    'hspf_blackstone42_huc10_point.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.57708,
      east_bound  =>  -71.57708,
      south_bound =>   42.01962,
      north_bound =>   42.01962
    },
    'hspf_blackstone42_huc10_polygon.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.78048,
      east_bound  =>  -71.35340,
      south_bound =>   41.85489,
      north_bound =>   42.23843
    },
    'hspf_blackstone42_huc12_point.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.71881,
      east_bound  =>  -71.39789,
      south_bound =>   41.91013,
      north_bound =>   42.15551
    },
    'hspf_blackstone42_huc12_polygon.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.78048,
      east_bound  =>  -71.35340,
      south_bound =>   41.85489,
      north_bound =>   42.23843
    },
    'hspf_blackstone43_huc10_point.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.57708,
      east_bound  =>  -71.57708,
      south_bound =>   42.01962,
      north_bound =>   42.01962
    },
    'hspf_blackstone43_huc10_polygon.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.78048,
      east_bound  =>  -71.35340,
      south_bound =>   41.85489,
      north_bound =>   42.23843
    },
    'hspf_blackstone43_huc12_point.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.71881,
      east_bound  =>  -71.39789,
      south_bound =>   41.91013,
      north_bound =>   42.15551
    },
    'hspf_blackstone43_huc12_polygon.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.78048,
      east_bound  =>  -71.35340,
      south_bound =>   41.85489,
      north_bound =>   42.23843
    },
    'hspf_blackstone9_huc10_point.hydrology' => {
      begin_date  => '1960-01-01',
      end_date    => '2009-12-31',
      west_bound  =>  -71.76534,
      east_bound  =>  -71.76534,
      south_bound =>   42.19291,
      north_bound =>   42.19291
    },
    'hspf_blackstone9_huc10_polygon.hydrology' => {
      begin_date  => '1960-01-01',
      end_date    => '2009-12-31',
      west_bound  =>  -71.92799,
      east_bound  =>  -71.60905,
      south_bound =>   42.03787,
      north_bound =>   42.35075
    },
    'hspf_blackstone9_huc12_point.hydrology' => {
      begin_date  => '1960-01-01',
      end_date    => '2009-12-31',
      west_bound  =>  -71.85803,
      east_bound  =>  -71.71627,
      south_bound =>   42.09478,
      north_bound =>   42.26947
    },
    'hspf_blackstone9_huc12_polygon.hydrology' => {
      begin_date  => '1960-01-01',
      end_date    => '2009-12-31',
      west_bound  =>  -71.92799,
      east_bound  =>  -71.60905,
      south_bound =>   42.03787,
      north_bound =>   42.35075
    },
    'hspf_charles11_huc10_point.hydrology' => {
      begin_date  => '1961-01-01',
      end_date    => '2009-12-31',
      west_bound  =>  -71.40021,
      east_bound  =>  -71.21325,
      south_bound =>   42.15069,
      north_bound =>   42.32819
    },
    'hspf_charles11_huc10_point.recharge_adjustment' => {
      begin_date  => '1961-01-01',
      end_date    => '2009-12-31',
      west_bound  =>  -71.40021,
      east_bound  =>  -71.21325,
      south_bound =>   42.15069,
      north_bound =>   42.32819
    },
    'hspf_charles11_huc10_polygon.hydrology' => {
      begin_date  => '1961-01-01',
      end_date    => '2009-12-31',
      west_bound  =>  -71.54945,
      east_bound  =>  -71.00937,
      south_bound =>   42.02564,
      north_bound =>   42.45390
    },
    'hspf_charles11_huc10_polygon.recharge_adjustment' => {
      begin_date  => '1961-01-01',
      end_date    => '2009-12-31',
      west_bound  =>  -71.54945,
      east_bound  =>  -71.00937,
      south_bound =>   42.02564,
      north_bound =>   42.45390
    },
    'hspf_charles11_huc12_point.hydrology' => {
      begin_date  => '1961-01-01',
      end_date    => '2009-12-31',
      west_bound  =>  -71.47076,
      east_bound  =>  -71.12661,
      south_bound =>   42.12082,
      north_bound =>   42.40220
    },
    'hspf_charles11_huc12_point.recharge_adjustment' => {
      begin_date  => '1961-01-01',
      end_date    => '2009-12-31',
      west_bound  =>  -71.47076,
      east_bound  =>  -71.12661,
      south_bound =>   42.12082,
      north_bound =>   42.40220
    },
    'hspf_charles11_huc12_polygon.hydrology' => {
      begin_date  => '1961-01-01',
      end_date    => '2009-12-31',
      west_bound  =>  -71.54945,
      east_bound  =>  -71.00937,
      south_bound =>   42.02564,
      north_bound =>   42.45390
    },
    'hspf_charles11_huc12_polygon.recharge_adjustment' => {
      begin_date  => '1961-01-01',
      end_date    => '2009-12-31',
      west_bound  =>  -71.54945,
      east_bound  =>  -71.00937,
      south_bound =>   42.02564,
      north_bound =>   42.45390
    },
    'hspf_charles3_huc10_point.hydrology' => {
      begin_date  => '1961-01-01',
      end_date    => '2009-12-31',
      west_bound  =>  -71.40021,
      east_bound  =>  -71.21325,
      south_bound =>   42.15069,
      north_bound =>   42.32819
    },
    'hspf_charles3_huc10_point.recharge_adjustment' => {
      begin_date  => '1961-01-01',
      end_date    => '2009-12-31',
      west_bound  =>  -71.40021,
      east_bound  =>  -71.21325,
      south_bound =>   42.15069,
      north_bound =>   42.32819
    },
    'hspf_charles3_huc10_polygon.hydrology' => {
      begin_date  => '1961-01-01',
      end_date    => '2009-12-31',
      west_bound  =>  -71.54945,
      east_bound  =>  -71.00937,
      south_bound =>   42.02564,
      north_bound =>   42.45390
    },
    'hspf_charles3_huc10_polygon.recharge_adjustment' => {
      begin_date  => '1961-01-01',
      end_date    => '2009-12-31',
      west_bound  =>  -71.54945,
      east_bound  =>  -71.00937,
      south_bound =>   42.02564,
      north_bound =>   42.45390
    },
    'hspf_charles3_huc12_point.hydrology' => {
      begin_date  => '1961-01-01',
      end_date    => '2009-12-31',
      west_bound  =>  -71.47076,
      east_bound  =>  -71.12661,
      south_bound =>   42.12082,
      north_bound =>   42.40220
    },
    'hspf_charles3_huc12_point.recharge_adjustment' => {
      begin_date  => '1961-01-01',
      end_date    => '2009-12-31',
      west_bound  =>  -71.47076,
      east_bound  =>  -71.12661,
      south_bound =>   42.12082,
      north_bound =>   42.40220
    },
    'hspf_charles3_huc12_polygon.hydrology' => {
      begin_date  => '1961-01-01',
      end_date    => '2009-12-31',
      west_bound  =>  -71.54945,
      east_bound  =>  -71.00937,
      south_bound =>   42.02564,
      north_bound =>   42.45390
    },
    'hspf_charles3_huc12_polygon.recharge_adjustment' => {
      begin_date  => '1961-01-01',
      end_date    => '2009-12-31',
      west_bound  =>  -71.54945,
      east_bound  =>  -71.00937,
      south_bound =>   42.02564,
      north_bound =>   42.45390
    },
    'hspf_charles44_huc10_point.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.40021,
      east_bound  =>  -71.21325,
      south_bound =>   42.15069,
      north_bound =>   42.32819
    },
    'hspf_charles44_huc10_point.recharge_adjustment' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.40021,
      east_bound  =>  -71.21325,
      south_bound =>   42.15069,
      north_bound =>   42.32819
    },
    'hspf_charles44_huc10_polygon.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.54945,
      east_bound  =>  -71.00937,
      south_bound =>   42.02564,
      north_bound =>   42.45390
    },
    'hspf_charles44_huc10_polygon.recharge_adjustment' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.54945,
      east_bound  =>  -71.00937,
      south_bound =>   42.02564,
      north_bound =>   42.45390
    },
    'hspf_charles44_huc12_point.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.47076,
      east_bound  =>  -71.12661,
      south_bound =>   42.12082,
      north_bound =>   42.40220
    },
    'hspf_charles44_huc12_point.recharge_adjustment' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.47076,
      east_bound  =>  -71.12661,
      south_bound =>   42.12082,
      north_bound =>   42.40220
    },
    'hspf_charles44_huc12_polygon.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.54945,
      east_bound  =>  -71.00937,
      south_bound =>   42.02564,
      north_bound =>   42.45390
    },
    'hspf_charles44_huc12_polygon.recharge_adjustment' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.54945,
      east_bound  =>  -71.00937,
      south_bound =>   42.02564,
      north_bound =>   42.45390
    },
    'hspf_charles45_huc10_point.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.40021,
      east_bound  =>  -71.21325,
      south_bound =>   42.15069,
      north_bound =>   42.32819
    },
    'hspf_charles45_huc10_point.recharge_adjustment' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.40021,
      east_bound  =>  -71.21325,
      south_bound =>   42.15069,
      north_bound =>   42.32819
    },
    'hspf_charles45_huc10_polygon.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.54945,
      east_bound  =>  -71.00937,
      south_bound =>   42.02564,
      north_bound =>   42.45390
    },
    'hspf_charles45_huc10_polygon.recharge_adjustment' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.54945,
      east_bound  =>  -71.00937,
      south_bound =>   42.02564,
      north_bound =>   42.45390
    },
    'hspf_charles45_huc12_point.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.47076,
      east_bound  =>  -71.12661,
      south_bound =>   42.12082,
      north_bound =>   42.40220
    },
    'hspf_charles45_huc12_point.recharge_adjustment' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.47076,
      east_bound  =>  -71.12661,
      south_bound =>   42.12082,
      north_bound =>   42.40220
    },
    'hspf_charles45_huc12_polygon.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.54945,
      east_bound  =>  -71.00937,
      south_bound =>   42.02564,
      north_bound =>   42.45390
    },
    'hspf_charles45_huc12_polygon.recharge_adjustment' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.54945,
      east_bound  =>  -71.00937,
      south_bound =>   42.02564,
      north_bound =>   42.45390
    },
    'hspf_charles46_huc10_point.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.40021,
      east_bound  =>  -71.21325,
      south_bound =>   42.15069,
      north_bound =>   42.32819
    },
    'hspf_charles46_huc10_point.recharge_adjustment' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.40021,
      east_bound  =>  -71.21325,
      south_bound =>   42.15069,
      north_bound =>   42.32819
    },
    'hspf_charles46_huc10_polygon.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.54945,
      east_bound  =>  -71.00937,
      south_bound =>   42.02564,
      north_bound =>   42.45390
    },
    'hspf_charles46_huc10_polygon.recharge_adjustment' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.54945,
      east_bound  =>  -71.00937,
      south_bound =>   42.02564,
      north_bound =>   42.45390
    },
    'hspf_charles46_huc12_point.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.47076,
      east_bound  =>  -71.12661,
      south_bound =>   42.12082,
      north_bound =>   42.40220
    },
    'hspf_charles46_huc12_point.recharge_adjustment' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.47076,
      east_bound  =>  -71.12661,
      south_bound =>   42.12082,
      north_bound =>   42.40220
    },
    'hspf_charles46_huc12_polygon.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.54945,
      east_bound  =>  -71.00937,
      south_bound =>   42.02564,
      north_bound =>   42.45390
    },
    'hspf_charles46_huc12_polygon.recharge_adjustment' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.54945,
      east_bound  =>  -71.00937,
      south_bound =>   42.02564,
      north_bound =>   42.45390
    },
    'hspf_charles47_huc10_point.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.40021,
      east_bound  =>  -71.21325,
      south_bound =>   42.15069,
      north_bound =>   42.32819
    },
    'hspf_charles47_huc10_point.recharge_adjustment' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.40021,
      east_bound  =>  -71.21325,
      south_bound =>   42.15069,
      north_bound =>   42.32819
    },
    'hspf_charles47_huc10_polygon.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.54945,
      east_bound  =>  -71.00937,
      south_bound =>   42.02564,
      north_bound =>   42.45390
    },
    'hspf_charles47_huc10_polygon.recharge_adjustment' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.54945,
      east_bound  =>  -71.00937,
      south_bound =>   42.02564,
      north_bound =>   42.45390
    },
    'hspf_charles47_huc12_point.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.47076,
      east_bound  =>  -71.12661,
      south_bound =>   42.12082,
      north_bound =>   42.40220
    },
    'hspf_charles47_huc12_point.recharge_adjustment' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.47076,
      east_bound  =>  -71.12661,
      south_bound =>   42.12082,
      north_bound =>   42.40220
    },
    'hspf_charles47_huc12_polygon.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.54945,
      east_bound  =>  -71.00937,
      south_bound =>   42.02564,
      north_bound =>   42.45390
    },
    'hspf_charles47_huc12_polygon.recharge_adjustment' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.54945,
      east_bound  =>  -71.00937,
      south_bound =>   42.02564,
      north_bound =>   42.45390
    },
    'hspf_charles48_huc10_point.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.40021,
      east_bound  =>  -71.21325,
      south_bound =>   42.15069,
      north_bound =>   42.32819
    },
    'hspf_charles48_huc10_point.recharge_adjustment' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.40021,
      east_bound  =>  -71.21325,
      south_bound =>   42.15069,
      north_bound =>   42.32819
    },
    'hspf_charles48_huc10_polygon.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.54945,
      east_bound  =>  -71.00937,
      south_bound =>   42.02564,
      north_bound =>   42.45390
    },
    'hspf_charles48_huc10_polygon.recharge_adjustment' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.54945,
      east_bound  =>  -71.00937,
      south_bound =>   42.02564,
      north_bound =>   42.45390
    },
    'hspf_charles48_huc12_point.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.47076,
      east_bound  =>  -71.12661,
      south_bound =>   42.12082,
      north_bound =>   42.40220
    },
    'hspf_charles48_huc12_point.recharge_adjustment' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.47076,
      east_bound  =>  -71.12661,
      south_bound =>   42.12082,
      north_bound =>   42.40220
    },
    'hspf_charles48_huc12_polygon.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.54945,
      east_bound  =>  -71.00937,
      south_bound =>   42.02564,
      north_bound =>   42.45390
    },
    'hspf_charles48_huc12_polygon.recharge_adjustment' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.54945,
      east_bound  =>  -71.00937,
      south_bound =>   42.02564,
      north_bound =>   42.45390
    },
    'hspf_farmington12_huc10_point.hydrology' => {
      begin_date  => '1970-01-01',
      end_date    => '2005-12-31',
      west_bound  =>  -73.11394,
      east_bound  =>  -72.77703,
      south_bound =>   41.77390,
      north_bound =>   42.10068
    },
    'hspf_farmington12_huc10_point.recharge_adjustment' => {
      begin_date  => '1970-01-01',
      end_date    => '2005-12-31',
      west_bound  =>  -73.11394,
      east_bound  =>  -72.77703,
      south_bound =>   41.77390,
      north_bound =>   42.10068
    },
    'hspf_farmington12_huc10_polygon.hydrology' => {
      begin_date  => '1970-01-01',
      end_date    => '2005-12-31',
      west_bound  =>  -73.20258,
      east_bound  =>  -72.62316,
      south_bound =>   41.63613,
      north_bound =>   42.30952
    },
    'hspf_farmington12_huc10_polygon.recharge_adjustment' => {
      begin_date  => '1970-01-01',
      end_date    => '2005-12-31',
      west_bound  =>  -73.20258,
      east_bound  =>  -72.62316,
      south_bound =>   41.63613,
      north_bound =>   42.30952
    },
    'hspf_farmington12_huc12_point.hydrology' => {
      begin_date  => '1970-01-01',
      end_date    => '2005-12-31',
      west_bound  =>  -73.12987,
      east_bound  =>  -72.75031,
      south_bound =>   41.69660,
      north_bound =>   42.19990
    },
    'hspf_farmington12_huc12_point.recharge_adjustment' => {
      begin_date  => '1970-01-01',
      end_date    => '2005-12-31',
      west_bound  =>  -73.12987,
      east_bound  =>  -72.75031,
      south_bound =>   41.69660,
      north_bound =>   42.19990
    },
    'hspf_farmington12_huc12_polygon.hydrology' => {
      begin_date  => '1970-01-01',
      end_date    => '2005-12-31',
      west_bound  =>  -73.20258,
      east_bound  =>  -72.62316,
      south_bound =>   41.63613,
      north_bound =>   42.30952
    },
    'hspf_farmington12_huc12_polygon.recharge_adjustment' => {
      begin_date  => '1970-01-01',
      end_date    => '2005-12-31',
      west_bound  =>  -73.20258,
      east_bound  =>  -72.62316,
      south_bound =>   41.63613,
      north_bound =>   42.30952
    },
    'hspf_farmington49_huc10_point.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -73.11394,
      east_bound  =>  -72.77703,
      south_bound =>   41.77390,
      north_bound =>   42.10068
    },
    'hspf_farmington49_huc10_point.recharge_adjustment' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -73.11394,
      east_bound  =>  -72.77703,
      south_bound =>   41.77390,
      north_bound =>   42.10068
    },
    'hspf_farmington49_huc10_polygon.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -73.20258,
      east_bound  =>  -72.62316,
      south_bound =>   41.63613,
      north_bound =>   42.30952
    },
    'hspf_farmington49_huc10_polygon.recharge_adjustment' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -73.20258,
      east_bound  =>  -72.62316,
      south_bound =>   41.63613,
      north_bound =>   42.30952
    },
    'hspf_farmington49_huc12_point.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -73.12987,
      east_bound  =>  -72.75031,
      south_bound =>   41.69660,
      north_bound =>   42.19990
    },
    'hspf_farmington49_huc12_point.recharge_adjustment' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -73.12987,
      east_bound  =>  -72.75031,
      south_bound =>   41.69660,
      north_bound =>   42.19990
    },
    'hspf_farmington49_huc12_polygon.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -73.20258,
      east_bound  =>  -72.62316,
      south_bound =>   41.63613,
      north_bound =>   42.30952
    },
    'hspf_farmington49_huc12_polygon.recharge_adjustment' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -73.20258,
      east_bound  =>  -72.62316,
      south_bound =>   41.63613,
      north_bound =>   42.30952
    },
    'hspf_farmington4_huc10_point.hydrology' => {
      begin_date  => '1970-01-01',
      end_date    => '2005-12-31',
      west_bound  =>  -73.11394,
      east_bound  =>  -72.77703,
      south_bound =>   41.77390,
      north_bound =>   42.10068
    },
    'hspf_farmington4_huc10_point.recharge_adjustment' => {
      begin_date  => '1970-01-01',
      end_date    => '2005-12-31',
      west_bound  =>  -73.11394,
      east_bound  =>  -72.77703,
      south_bound =>   41.77390,
      north_bound =>   42.10068
    },
    'hspf_farmington4_huc10_polygon.hydrology' => {
      begin_date  => '1970-01-01',
      end_date    => '2005-12-31',
      west_bound  =>  -73.20258,
      east_bound  =>  -72.62316,
      south_bound =>   41.63613,
      north_bound =>   42.30952
    },
    'hspf_farmington4_huc10_polygon.recharge_adjustment' => {
      begin_date  => '1970-01-01',
      end_date    => '2005-12-31',
      west_bound  =>  -73.20258,
      east_bound  =>  -72.62316,
      south_bound =>   41.63613,
      north_bound =>   42.30952
    },
    'hspf_farmington4_huc12_point.hydrology' => {
      begin_date  => '1970-01-01',
      end_date    => '2005-12-31',
      west_bound  =>  -73.12987,
      east_bound  =>  -72.75031,
      south_bound =>   41.69660,
      north_bound =>   42.19990
    },
    'hspf_farmington4_huc12_point.recharge_adjustment' => {
      begin_date  => '1970-01-01',
      end_date    => '2005-12-31',
      west_bound  =>  -73.12987,
      east_bound  =>  -72.75031,
      south_bound =>   41.69660,
      north_bound =>   42.19990
    },
    'hspf_farmington4_huc12_polygon.hydrology' => {
      begin_date  => '1970-01-01',
      end_date    => '2005-12-31',
      west_bound  =>  -73.20258,
      east_bound  =>  -72.62316,
      south_bound =>   41.63613,
      north_bound =>   42.30952
    },
    'hspf_farmington4_huc12_polygon.recharge_adjustment' => {
      begin_date  => '1970-01-01',
      end_date    => '2005-12-31',
      west_bound  =>  -73.20258,
      east_bound  =>  -72.62316,
      south_bound =>   41.63613,
      north_bound =>   42.30952
    },
    'hspf_farmington50_huc10_point.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -73.11394,
      east_bound  =>  -72.77703,
      south_bound =>   41.77390,
      north_bound =>   42.10068
    },
    'hspf_farmington50_huc10_point.recharge_adjustment' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -73.11394,
      east_bound  =>  -72.77703,
      south_bound =>   41.77390,
      north_bound =>   42.10068
    },
    'hspf_farmington50_huc10_polygon.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -73.20258,
      east_bound  =>  -72.62316,
      south_bound =>   41.63613,
      north_bound =>   42.30952
    },
    'hspf_farmington50_huc10_polygon.recharge_adjustment' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -73.20258,
      east_bound  =>  -72.62316,
      south_bound =>   41.63613,
      north_bound =>   42.30952
    },
    'hspf_farmington50_huc12_point.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -73.12987,
      east_bound  =>  -72.75031,
      south_bound =>   41.69660,
      north_bound =>   42.19990
    },
    'hspf_farmington50_huc12_point.recharge_adjustment' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -73.12987,
      east_bound  =>  -72.75031,
      south_bound =>   41.69660,
      north_bound =>   42.19990
    },
    'hspf_farmington50_huc12_polygon.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -73.20258,
      east_bound  =>  -72.62316,
      south_bound =>   41.63613,
      north_bound =>   42.30952
    },
    'hspf_farmington50_huc12_polygon.recharge_adjustment' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -73.20258,
      east_bound  =>  -72.62316,
      south_bound =>   41.63613,
      north_bound =>   42.30952
    },
    'hspf_farmington51_huc10_point.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -73.11394,
      east_bound  =>  -72.77703,
      south_bound =>   41.77390,
      north_bound =>   42.10068
    },
    'hspf_farmington51_huc10_point.recharge_adjustment' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -73.11394,
      east_bound  =>  -72.77703,
      south_bound =>   41.77390,
      north_bound =>   42.10068
    },
    'hspf_farmington51_huc10_polygon.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -73.20258,
      east_bound  =>  -72.62316,
      south_bound =>   41.63613,
      north_bound =>   42.30952
    },
    'hspf_farmington51_huc10_polygon.recharge_adjustment' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -73.20258,
      east_bound  =>  -72.62316,
      south_bound =>   41.63613,
      north_bound =>   42.30952
    },
    'hspf_farmington51_huc12_point.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -73.12987,
      east_bound  =>  -72.75031,
      south_bound =>   41.69660,
      north_bound =>   42.19990
    },
    'hspf_farmington51_huc12_point.recharge_adjustment' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -73.12987,
      east_bound  =>  -72.75031,
      south_bound =>   41.69660,
      north_bound =>   42.19990
    },
    'hspf_farmington51_huc12_polygon.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -73.20258,
      east_bound  =>  -72.62316,
      south_bound =>   41.63613,
      north_bound =>   42.30952
    },
    'hspf_farmington51_huc12_polygon.recharge_adjustment' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -73.20258,
      east_bound  =>  -72.62316,
      south_bound =>   41.63613,
      north_bound =>   42.30952
    },
    'hspf_farmington52_huc10_point.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -73.11394,
      east_bound  =>  -72.77703,
      south_bound =>   41.77390,
      north_bound =>   42.10068
    },
    'hspf_farmington52_huc10_point.recharge_adjustment' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -73.11394,
      east_bound  =>  -72.77703,
      south_bound =>   41.77390,
      north_bound =>   42.10068
    },
    'hspf_farmington52_huc10_polygon.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -73.20258,
      east_bound  =>  -72.62316,
      south_bound =>   41.63613,
      north_bound =>   42.30952
    },
    'hspf_farmington52_huc10_polygon.recharge_adjustment' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -73.20258,
      east_bound  =>  -72.62316,
      south_bound =>   41.63613,
      north_bound =>   42.30952
    },
    'hspf_farmington52_huc12_point.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -73.12987,
      east_bound  =>  -72.75031,
      south_bound =>   41.69660,
      north_bound =>   42.19990
    },
    'hspf_farmington52_huc12_point.recharge_adjustment' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -73.12987,
      east_bound  =>  -72.75031,
      south_bound =>   41.69660,
      north_bound =>   42.19990
    },
    'hspf_farmington52_huc12_polygon.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -73.20258,
      east_bound  =>  -72.62316,
      south_bound =>   41.63613,
      north_bound =>   42.30952
    },
    'hspf_farmington52_huc12_polygon.recharge_adjustment' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -73.20258,
      east_bound  =>  -72.62316,
      south_bound =>   41.63613,
      north_bound =>   42.30952
    },
    'hspf_farmington53_huc10_point.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -73.11394,
      east_bound  =>  -72.77703,
      south_bound =>   41.77390,
      north_bound =>   42.10068
    },
    'hspf_farmington53_huc10_point.recharge_adjustment' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -73.11394,
      east_bound  =>  -72.77703,
      south_bound =>   41.77390,
      north_bound =>   42.10068
    },
    'hspf_farmington53_huc10_polygon.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -73.20258,
      east_bound  =>  -72.62316,
      south_bound =>   41.63613,
      north_bound =>   42.30952
    },
    'hspf_farmington53_huc10_polygon.recharge_adjustment' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -73.20258,
      east_bound  =>  -72.62316,
      south_bound =>   41.63613,
      north_bound =>   42.30952
    },
    'hspf_farmington53_huc12_point.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -73.12987,
      east_bound  =>  -72.75031,
      south_bound =>   41.69660,
      north_bound =>   42.19990
    },
    'hspf_farmington53_huc12_point.recharge_adjustment' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -73.12987,
      east_bound  =>  -72.75031,
      south_bound =>   41.69660,
      north_bound =>   42.19990
    },
    'hspf_farmington53_huc12_polygon.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -73.20258,
      east_bound  =>  -72.62316,
      south_bound =>   41.63613,
      north_bound =>   42.30952
    },
    'hspf_farmington53_huc12_polygon.recharge_adjustment' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -73.20258,
      east_bound  =>  -72.62316,
      south_bound =>   41.63613,
      north_bound =>   42.30952
    },
    'hspf_ipswich13_huc10_point.hydrology' => {
      begin_date  => '1961-01-01',
      end_date    => '2009-12-31',
      west_bound  =>  -71.13159,
      east_bound  =>  -71.13159,
      south_bound =>   42.56714,
      north_bound =>   42.56714
    },
    'hspf_ipswich13_huc10_point.recharge_adjustment' => {
      begin_date  => '1961-01-01',
      end_date    => '2009-12-31',
      west_bound  =>  -71.13159,
      east_bound  =>  -71.13159,
      south_bound =>   42.56714,
      north_bound =>   42.56714
    },
    'hspf_ipswich13_huc10_polygon.hydrology' => {
      begin_date  => '1961-01-01',
      end_date    => '2009-12-31',
      west_bound  =>  -71.21963,
      east_bound  =>  -70.81332,
      south_bound =>   42.49842,
      north_bound =>   42.70367
    },
    'hspf_ipswich13_huc10_polygon.recharge_adjustment' => {
      begin_date  => '1961-01-01',
      end_date    => '2009-12-31',
      west_bound  =>  -71.21963,
      east_bound  =>  -70.81332,
      south_bound =>   42.49842,
      north_bound =>   42.70367
    },
    'hspf_ipswich13_huc12_point.hydrology' => {
      begin_date  => '1961-01-01',
      end_date    => '2009-12-31',
      west_bound  =>  -71.13159,
      east_bound  =>  -70.89173,
      south_bound =>   42.56714,
      north_bound =>   42.67464
    },
    'hspf_ipswich13_huc12_point.recharge_adjustment' => {
      begin_date  => '1961-01-01',
      end_date    => '2009-12-31',
      west_bound  =>  -71.13159,
      east_bound  =>  -70.89173,
      south_bound =>   42.56714,
      north_bound =>   42.67464
    },
    'hspf_ipswich13_huc12_polygon.hydrology' => {
      begin_date  => '1961-01-01',
      end_date    => '2009-12-31',
      west_bound  =>  -71.21963,
      east_bound  =>  -70.81332,
      south_bound =>   42.49842,
      north_bound =>   42.70367
    },
    'hspf_ipswich13_huc12_polygon.recharge_adjustment' => {
      begin_date  => '1961-01-01',
      end_date    => '2009-12-31',
      west_bound  =>  -71.21963,
      east_bound  =>  -70.81332,
      south_bound =>   42.49842,
      north_bound =>   42.70367
    },
    'hspf_ipswich54_huc10_point.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.13159,
      east_bound  =>  -71.13159,
      south_bound =>   42.56714,
      north_bound =>   42.56714
    },
    'hspf_ipswich54_huc10_point.recharge_adjustment' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.13159,
      east_bound  =>  -71.13159,
      south_bound =>   42.56714,
      north_bound =>   42.56714
    },
    'hspf_ipswich54_huc10_polygon.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.21963,
      east_bound  =>  -70.81332,
      south_bound =>   42.49842,
      north_bound =>   42.70367
    },
    'hspf_ipswich54_huc10_polygon.recharge_adjustment' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.21963,
      east_bound  =>  -70.81332,
      south_bound =>   42.49842,
      north_bound =>   42.70367
    },
    'hspf_ipswich54_huc12_point.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.13159,
      east_bound  =>  -70.89173,
      south_bound =>   42.56714,
      north_bound =>   42.67464
    },
    'hspf_ipswich54_huc12_point.recharge_adjustment' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.13159,
      east_bound  =>  -70.89173,
      south_bound =>   42.56714,
      north_bound =>   42.67464
    },
    'hspf_ipswich54_huc12_polygon.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.21963,
      east_bound  =>  -70.81332,
      south_bound =>   42.49842,
      north_bound =>   42.70367
    },
    'hspf_ipswich54_huc12_polygon.recharge_adjustment' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.21963,
      east_bound  =>  -70.81332,
      south_bound =>   42.49842,
      north_bound =>   42.70367
    },
    'hspf_ipswich55_huc10_point.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.13159,
      east_bound  =>  -71.13159,
      south_bound =>   42.56714,
      north_bound =>   42.56714
    },
    'hspf_ipswich55_huc10_point.recharge_adjustment' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.13159,
      east_bound  =>  -71.13159,
      south_bound =>   42.56714,
      north_bound =>   42.56714
    },
    'hspf_ipswich55_huc10_polygon.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.21963,
      east_bound  =>  -70.81332,
      south_bound =>   42.49842,
      north_bound =>   42.70367
    },
    'hspf_ipswich55_huc10_polygon.recharge_adjustment' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.21963,
      east_bound  =>  -70.81332,
      south_bound =>   42.49842,
      north_bound =>   42.70367
    },
    'hspf_ipswich55_huc12_point.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.13159,
      east_bound  =>  -70.89173,
      south_bound =>   42.56714,
      north_bound =>   42.67464
    },
    'hspf_ipswich55_huc12_point.recharge_adjustment' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.13159,
      east_bound  =>  -70.89173,
      south_bound =>   42.56714,
      north_bound =>   42.67464
    },
    'hspf_ipswich55_huc12_polygon.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.21963,
      east_bound  =>  -70.81332,
      south_bound =>   42.49842,
      north_bound =>   42.70367
    },
    'hspf_ipswich55_huc12_polygon.recharge_adjustment' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.21963,
      east_bound  =>  -70.81332,
      south_bound =>   42.49842,
      north_bound =>   42.70367
    },
    'hspf_ipswich56_huc10_point.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.13159,
      east_bound  =>  -71.13159,
      south_bound =>   42.56714,
      north_bound =>   42.56714
    },
    'hspf_ipswich56_huc10_point.recharge_adjustment' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.13159,
      east_bound  =>  -71.13159,
      south_bound =>   42.56714,
      north_bound =>   42.56714
    },
    'hspf_ipswich56_huc10_polygon.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.21963,
      east_bound  =>  -70.81332,
      south_bound =>   42.49842,
      north_bound =>   42.70367
    },
    'hspf_ipswich56_huc10_polygon.recharge_adjustment' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.21963,
      east_bound  =>  -70.81332,
      south_bound =>   42.49842,
      north_bound =>   42.70367
    },
    'hspf_ipswich56_huc12_point.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.13159,
      east_bound  =>  -70.89173,
      south_bound =>   42.56714,
      north_bound =>   42.67464
    },
    'hspf_ipswich56_huc12_point.recharge_adjustment' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.13159,
      east_bound  =>  -70.89173,
      south_bound =>   42.56714,
      north_bound =>   42.67464
    },
    'hspf_ipswich56_huc12_polygon.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.21963,
      east_bound  =>  -70.81332,
      south_bound =>   42.49842,
      north_bound =>   42.70367
    },
    'hspf_ipswich56_huc12_polygon.recharge_adjustment' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.21963,
      east_bound  =>  -70.81332,
      south_bound =>   42.49842,
      north_bound =>   42.70367
    },
    'hspf_ipswich57_huc10_point.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.13159,
      east_bound  =>  -71.13159,
      south_bound =>   42.56714,
      north_bound =>   42.56714
    },
    'hspf_ipswich57_huc10_point.recharge_adjustment' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.13159,
      east_bound  =>  -71.13159,
      south_bound =>   42.56714,
      north_bound =>   42.56714
    },
    'hspf_ipswich57_huc10_polygon.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.21963,
      east_bound  =>  -70.81332,
      south_bound =>   42.49842,
      north_bound =>   42.70367
    },
    'hspf_ipswich57_huc10_polygon.recharge_adjustment' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.21963,
      east_bound  =>  -70.81332,
      south_bound =>   42.49842,
      north_bound =>   42.70367
    },
    'hspf_ipswich57_huc12_point.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.13159,
      east_bound  =>  -70.89173,
      south_bound =>   42.56714,
      north_bound =>   42.67464
    },
    'hspf_ipswich57_huc12_point.recharge_adjustment' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.13159,
      east_bound  =>  -70.89173,
      south_bound =>   42.56714,
      north_bound =>   42.67464
    },
    'hspf_ipswich57_huc12_polygon.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.21963,
      east_bound  =>  -70.81332,
      south_bound =>   42.49842,
      north_bound =>   42.70367
    },
    'hspf_ipswich57_huc12_polygon.recharge_adjustment' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.21963,
      east_bound  =>  -70.81332,
      south_bound =>   42.49842,
      north_bound =>   42.70367
    },
    'hspf_ipswich58_huc10_point.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.13159,
      east_bound  =>  -71.13159,
      south_bound =>   42.56714,
      north_bound =>   42.56714
    },
    'hspf_ipswich58_huc10_point.recharge_adjustment' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.13159,
      east_bound  =>  -71.13159,
      south_bound =>   42.56714,
      north_bound =>   42.56714
    },
    'hspf_ipswich58_huc10_polygon.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.21963,
      east_bound  =>  -70.81332,
      south_bound =>   42.49842,
      north_bound =>   42.70367
    },
    'hspf_ipswich58_huc10_polygon.recharge_adjustment' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.21963,
      east_bound  =>  -70.81332,
      south_bound =>   42.49842,
      north_bound =>   42.70367
    },
    'hspf_ipswich58_huc12_point.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.13159,
      east_bound  =>  -70.89173,
      south_bound =>   42.56714,
      north_bound =>   42.67464
    },
    'hspf_ipswich58_huc12_point.recharge_adjustment' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.13159,
      east_bound  =>  -70.89173,
      south_bound =>   42.56714,
      north_bound =>   42.67464
    },
    'hspf_ipswich58_huc12_polygon.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.21963,
      east_bound  =>  -70.81332,
      south_bound =>   42.49842,
      north_bound =>   42.70367
    },
    'hspf_ipswich58_huc12_polygon.recharge_adjustment' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.21963,
      east_bound  =>  -70.81332,
      south_bound =>   42.49842,
      north_bound =>   42.70367
    },
    'hspf_ipswich5_huc10_point.hydrology' => {
      begin_date  => '1961-01-01',
      end_date    => '2009-12-31',
      west_bound  =>  -71.13159,
      east_bound  =>  -71.13159,
      south_bound =>   42.56714,
      north_bound =>   42.56714
    },
    'hspf_ipswich5_huc10_point.recharge_adjustment' => {
      begin_date  => '1961-01-01',
      end_date    => '2009-12-31',
      west_bound  =>  -71.13159,
      east_bound  =>  -71.13159,
      south_bound =>   42.56714,
      north_bound =>   42.56714
    },
    'hspf_ipswich5_huc10_polygon.hydrology' => {
      begin_date  => '1961-01-01',
      end_date    => '2009-12-31',
      west_bound  =>  -71.21963,
      east_bound  =>  -70.81332,
      south_bound =>   42.49842,
      north_bound =>   42.70367
    },
    'hspf_ipswich5_huc10_polygon.recharge_adjustment' => {
      begin_date  => '1961-01-01',
      end_date    => '2009-12-31',
      west_bound  =>  -71.21963,
      east_bound  =>  -70.81332,
      south_bound =>   42.49842,
      north_bound =>   42.70367
    },
    'hspf_ipswich5_huc12_point.hydrology' => {
      begin_date  => '1961-01-01',
      end_date    => '2009-12-31',
      west_bound  =>  -71.13159,
      east_bound  =>  -70.89173,
      south_bound =>   42.56714,
      north_bound =>   42.67464
    },
    'hspf_ipswich5_huc12_point.recharge_adjustment' => {
      begin_date  => '1961-01-01',
      end_date    => '2009-12-31',
      west_bound  =>  -71.13159,
      east_bound  =>  -70.89173,
      south_bound =>   42.56714,
      north_bound =>   42.67464
    },
    'hspf_ipswich5_huc12_polygon.hydrology' => {
      begin_date  => '1961-01-01',
      end_date    => '2009-12-31',
      west_bound  =>  -71.21963,
      east_bound  =>  -70.81332,
      south_bound =>   42.49842,
      north_bound =>   42.70367
    },
    'hspf_ipswich5_huc12_polygon.recharge_adjustment' => {
      begin_date  => '1961-01-01',
      end_date    => '2009-12-31',
      west_bound  =>  -71.21963,
      east_bound  =>  -70.81332,
      south_bound =>   42.49842,
      north_bound =>   42.70367
    },
    'hspf_pawcatuck14_huc10_point.hydrology' => {
      begin_date  => '1960-01-01',
      end_date    => '1999-12-31',
      west_bound  =>  -71.84246,
      east_bound  =>  -71.63932,
      south_bound =>   41.44597,
      north_bound =>   41.56840
    },
    'hspf_pawcatuck14_huc10_point.recharge_adjustment' => {
      begin_date  => '1960-01-01',
      end_date    => '2009-12-31',
      west_bound  =>  -71.84246,
      east_bound  =>  -71.63932,
      south_bound =>   41.44597,
      north_bound =>   41.56840
    },
    'hspf_pawcatuck14_huc10_polygon.hydrology' => {
      begin_date  => '1960-01-01',
      end_date    => '1999-12-31',
      west_bound  =>  -71.93323,
      east_bound  =>  -71.48661,
      south_bound =>   41.31577,
      north_bound =>   41.68141
    },
    'hspf_pawcatuck14_huc10_polygon.recharge_adjustment' => {
      begin_date  => '1960-01-01',
      end_date    => '2009-12-31',
      west_bound  =>  -71.93323,
      east_bound  =>  -71.48661,
      south_bound =>   41.31577,
      north_bound =>   41.68141
    },
    'hspf_pawcatuck14_huc12_point.hydrology' => {
      begin_date  => '1960-01-01',
      end_date    => '1999-12-31',
      west_bound  =>  -71.88957,
      east_bound  =>  -71.54656,
      south_bound =>   41.38083,
      north_bound =>   41.59852
    },
    'hspf_pawcatuck14_huc12_point.recharge_adjustment' => {
      begin_date  => '1960-01-01',
      end_date    => '2009-12-31',
      west_bound  =>  -71.88957,
      east_bound  =>  -71.54656,
      south_bound =>   41.38083,
      north_bound =>   41.59852
    },
    'hspf_pawcatuck14_huc12_polygon.hydrology' => {
      begin_date  => '1960-01-01',
      end_date    => '1999-12-31',
      west_bound  =>  -71.93323,
      east_bound  =>  -71.48661,
      south_bound =>   41.31577,
      north_bound =>   41.68141
    },
    'hspf_pawcatuck14_huc12_polygon.recharge_adjustment' => {
      begin_date  => '1960-01-01',
      end_date    => '2009-12-31',
      west_bound  =>  -71.93323,
      east_bound  =>  -71.48661,
      south_bound =>   41.31577,
      north_bound =>   41.68141
    },
    'hspf_pawcatuck59_huc10_point.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.84246,
      east_bound  =>  -71.63932,
      south_bound =>   41.44597,
      north_bound =>   41.56840
    },
    'hspf_pawcatuck59_huc10_point.recharge_adjustment' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.84246,
      east_bound  =>  -71.63932,
      south_bound =>   41.44597,
      north_bound =>   41.56840
    },
    'hspf_pawcatuck59_huc10_polygon.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.93323,
      east_bound  =>  -71.48661,
      south_bound =>   41.31577,
      north_bound =>   41.68141
    },
    'hspf_pawcatuck59_huc10_polygon.recharge_adjustment' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.93323,
      east_bound  =>  -71.48661,
      south_bound =>   41.31577,
      north_bound =>   41.68141
    },
    'hspf_pawcatuck59_huc12_point.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.88957,
      east_bound  =>  -71.54656,
      south_bound =>   41.38083,
      north_bound =>   41.59852
    },
    'hspf_pawcatuck59_huc12_point.recharge_adjustment' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.88957,
      east_bound  =>  -71.54656,
      south_bound =>   41.38083,
      north_bound =>   41.59852
    },
    'hspf_pawcatuck59_huc12_polygon.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.93323,
      east_bound  =>  -71.48661,
      south_bound =>   41.31577,
      north_bound =>   41.68141
    },
    'hspf_pawcatuck59_huc12_polygon.recharge_adjustment' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.93323,
      east_bound  =>  -71.48661,
      south_bound =>   41.31577,
      north_bound =>   41.68141
    },
    'hspf_pawcatuck60_huc10_point.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.84246,
      east_bound  =>  -71.63932,
      south_bound =>   41.44597,
      north_bound =>   41.56840
    },
    'hspf_pawcatuck60_huc10_point.recharge_adjustment' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.84246,
      east_bound  =>  -71.63932,
      south_bound =>   41.44597,
      north_bound =>   41.56840
    },
    'hspf_pawcatuck60_huc10_polygon.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.93323,
      east_bound  =>  -71.48661,
      south_bound =>   41.31577,
      north_bound =>   41.68141
    },
    'hspf_pawcatuck60_huc10_polygon.recharge_adjustment' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.93323,
      east_bound  =>  -71.48661,
      south_bound =>   41.31577,
      north_bound =>   41.68141
    },
    'hspf_pawcatuck60_huc12_point.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.88957,
      east_bound  =>  -71.54656,
      south_bound =>   41.38083,
      north_bound =>   41.59852
    },
    'hspf_pawcatuck60_huc12_point.recharge_adjustment' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.88957,
      east_bound  =>  -71.54656,
      south_bound =>   41.38083,
      north_bound =>   41.59852
    },
    'hspf_pawcatuck60_huc12_polygon.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.93323,
      east_bound  =>  -71.48661,
      south_bound =>   41.31577,
      north_bound =>   41.68141
    },
    'hspf_pawcatuck60_huc12_polygon.recharge_adjustment' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.93323,
      east_bound  =>  -71.48661,
      south_bound =>   41.31577,
      north_bound =>   41.68141
    },
    'hspf_pawcatuck61_huc10_point.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.84246,
      east_bound  =>  -71.63932,
      south_bound =>   41.44597,
      north_bound =>   41.56840
    },
    'hspf_pawcatuck61_huc10_point.recharge_adjustment' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.84246,
      east_bound  =>  -71.63932,
      south_bound =>   41.44597,
      north_bound =>   41.56840
    },
    'hspf_pawcatuck61_huc10_polygon.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.93323,
      east_bound  =>  -71.48661,
      south_bound =>   41.31577,
      north_bound =>   41.68141
    },
    'hspf_pawcatuck61_huc10_polygon.recharge_adjustment' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.93323,
      east_bound  =>  -71.48661,
      south_bound =>   41.31577,
      north_bound =>   41.68141
    },
    'hspf_pawcatuck61_huc12_point.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.88957,
      east_bound  =>  -71.54656,
      south_bound =>   41.38083,
      north_bound =>   41.59852
    },
    'hspf_pawcatuck61_huc12_point.recharge_adjustment' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.88957,
      east_bound  =>  -71.54656,
      south_bound =>   41.38083,
      north_bound =>   41.59852
    },
    'hspf_pawcatuck61_huc12_polygon.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.93323,
      east_bound  =>  -71.48661,
      south_bound =>   41.31577,
      north_bound =>   41.68141
    },
    'hspf_pawcatuck61_huc12_polygon.recharge_adjustment' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.93323,
      east_bound  =>  -71.48661,
      south_bound =>   41.31577,
      north_bound =>   41.68141
    },
    'hspf_pawcatuck62_huc10_point.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.84246,
      east_bound  =>  -71.63932,
      south_bound =>   41.44597,
      north_bound =>   41.56840
    },
    'hspf_pawcatuck62_huc10_point.recharge_adjustment' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.84246,
      east_bound  =>  -71.63932,
      south_bound =>   41.44597,
      north_bound =>   41.56840
    },
    'hspf_pawcatuck62_huc10_polygon.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.93323,
      east_bound  =>  -71.48661,
      south_bound =>   41.31577,
      north_bound =>   41.68141
    },
    'hspf_pawcatuck62_huc10_polygon.recharge_adjustment' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.93323,
      east_bound  =>  -71.48661,
      south_bound =>   41.31577,
      north_bound =>   41.68141
    },
    'hspf_pawcatuck62_huc12_point.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.88957,
      east_bound  =>  -71.54656,
      south_bound =>   41.38083,
      north_bound =>   41.59852
    },
    'hspf_pawcatuck62_huc12_point.recharge_adjustment' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.88957,
      east_bound  =>  -71.54656,
      south_bound =>   41.38083,
      north_bound =>   41.59852
    },
    'hspf_pawcatuck62_huc12_polygon.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.93323,
      east_bound  =>  -71.48661,
      south_bound =>   41.31577,
      north_bound =>   41.68141
    },
    'hspf_pawcatuck62_huc12_polygon.recharge_adjustment' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.93323,
      east_bound  =>  -71.48661,
      south_bound =>   41.31577,
      north_bound =>   41.68141
    },
    'hspf_pawcatuck63_huc10_point.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.84246,
      east_bound  =>  -71.63932,
      south_bound =>   41.44597,
      north_bound =>   41.56840
    },
    'hspf_pawcatuck63_huc10_point.recharge_adjustment' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.84246,
      east_bound  =>  -71.63932,
      south_bound =>   41.44597,
      north_bound =>   41.56840
    },
    'hspf_pawcatuck63_huc10_polygon.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.93323,
      east_bound  =>  -71.48661,
      south_bound =>   41.31577,
      north_bound =>   41.68141
    },
    'hspf_pawcatuck63_huc10_polygon.recharge_adjustment' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.93323,
      east_bound  =>  -71.48661,
      south_bound =>   41.31577,
      north_bound =>   41.68141
    },
    'hspf_pawcatuck63_huc12_point.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.88957,
      east_bound  =>  -71.54656,
      south_bound =>   41.38083,
      north_bound =>   41.59852
    },
    'hspf_pawcatuck63_huc12_point.recharge_adjustment' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.88957,
      east_bound  =>  -71.54656,
      south_bound =>   41.38083,
      north_bound =>   41.59852
    },
    'hspf_pawcatuck63_huc12_polygon.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.93323,
      east_bound  =>  -71.48661,
      south_bound =>   41.31577,
      north_bound =>   41.68141
    },
    'hspf_pawcatuck63_huc12_polygon.recharge_adjustment' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.93323,
      east_bound  =>  -71.48661,
      south_bound =>   41.31577,
      north_bound =>   41.68141
    },
    'hspf_pawcatuck6_huc10_point.hydrology' => {
      begin_date  => '1960-01-01',
      end_date    => '1999-12-31',
      west_bound  =>  -71.84246,
      east_bound  =>  -71.63932,
      south_bound =>   41.44597,
      north_bound =>   41.56840
    },
    'hspf_pawcatuck6_huc10_point.recharge_adjustment' => {
      begin_date  => '1960-01-01',
      end_date    => '2009-12-31',
      west_bound  =>  -71.84246,
      east_bound  =>  -71.63932,
      south_bound =>   41.44597,
      north_bound =>   41.56840
    },
    'hspf_pawcatuck6_huc10_polygon.hydrology' => {
      begin_date  => '1960-01-01',
      end_date    => '1999-12-31',
      west_bound  =>  -71.93323,
      east_bound  =>  -71.48661,
      south_bound =>   41.31577,
      north_bound =>   41.68141
    },
    'hspf_pawcatuck6_huc10_polygon.recharge_adjustment' => {
      begin_date  => '1960-01-01',
      end_date    => '2009-12-31',
      west_bound  =>  -71.93323,
      east_bound  =>  -71.48661,
      south_bound =>   41.31577,
      north_bound =>   41.68141
    },
    'hspf_pawcatuck6_huc12_point.hydrology' => {
      begin_date  => '1960-01-01',
      end_date    => '1999-12-31',
      west_bound  =>  -71.88957,
      east_bound  =>  -71.54656,
      south_bound =>   41.38083,
      north_bound =>   41.59852
    },
    'hspf_pawcatuck6_huc12_point.recharge_adjustment' => {
      begin_date  => '1960-01-01',
      end_date    => '2009-12-31',
      west_bound  =>  -71.88957,
      east_bound  =>  -71.54656,
      south_bound =>   41.38083,
      north_bound =>   41.59852
    },
    'hspf_pawcatuck6_huc12_polygon.hydrology' => {
      begin_date  => '1960-01-01',
      end_date    => '1999-12-31',
      west_bound  =>  -71.93323,
      east_bound  =>  -71.48661,
      south_bound =>   41.31577,
      north_bound =>   41.68141
    },
    'hspf_pawcatuck6_huc12_polygon.recharge_adjustment' => {
      begin_date  => '1960-01-01',
      end_date    => '2009-12-31',
      west_bound  =>  -71.93323,
      east_bound  =>  -71.48661,
      south_bound =>   41.31577,
      north_bound =>   41.68141
    },
    'hspf_sudbury15_huc10_point.hydrology' => {
      begin_date  => '1960-01-01',
      end_date    => '2009-12-31',
      west_bound  =>  -71.48964,
      east_bound  =>  -71.46364,
      south_bound =>   42.31622,
      north_bound =>   42.44976
    },
    'hspf_sudbury15_huc10_point.recharge_adjustment' => {
      begin_date  => '1960-01-01',
      end_date    => '2009-12-31',
      west_bound  =>  -71.48964,
      east_bound  =>  -71.46364,
      south_bound =>   42.31622,
      north_bound =>   42.44976
    },
    'hspf_sudbury15_huc10_polygon.hydrology' => {
      begin_date  => '1960-01-01',
      end_date    => '2009-12-31',
      west_bound  =>  -71.72177,
      east_bound  =>  -71.25150,
      south_bound =>   42.19754,
      north_bound =>   42.64674
    },
    'hspf_sudbury15_huc10_polygon.recharge_adjustment' => {
      begin_date  => '1960-01-01',
      end_date    => '2009-12-31',
      west_bound  =>  -71.72177,
      east_bound  =>  -71.25150,
      south_bound =>   42.19754,
      north_bound =>   42.64674
    },
    'hspf_sudbury15_huc12_point.hydrology' => {
      begin_date  => '1960-01-01',
      end_date    => '2009-12-31',
      west_bound  =>  -71.65187,
      east_bound  =>  -71.30854,
      south_bound =>   42.25174,
      north_bound =>   42.57971
    },
    'hspf_sudbury15_huc12_point.recharge_adjustment' => {
      begin_date  => '1960-01-01',
      end_date    => '2009-12-31',
      west_bound  =>  -71.65187,
      east_bound  =>  -71.30854,
      south_bound =>   42.25174,
      north_bound =>   42.57971
    },
    'hspf_sudbury15_huc12_polygon.hydrology' => {
      begin_date  => '1960-01-01',
      end_date    => '2009-12-31',
      west_bound  =>  -71.72177,
      east_bound  =>  -71.25150,
      south_bound =>   42.19754,
      north_bound =>   42.64674
    },
    'hspf_sudbury15_huc12_polygon.recharge_adjustment' => {
      begin_date  => '1960-01-01',
      end_date    => '2009-12-31',
      west_bound  =>  -71.72177,
      east_bound  =>  -71.25150,
      south_bound =>   42.19754,
      north_bound =>   42.64674
    },
    'hspf_sudbury64_huc10_point.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.48964,
      east_bound  =>  -71.46364,
      south_bound =>   42.31622,
      north_bound =>   42.44976
    },
    'hspf_sudbury64_huc10_point.recharge_adjustment' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.48964,
      east_bound  =>  -71.46364,
      south_bound =>   42.31622,
      north_bound =>   42.44976
    },
    'hspf_sudbury64_huc10_polygon.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.72177,
      east_bound  =>  -71.25150,
      south_bound =>   42.19754,
      north_bound =>   42.64674
    },
    'hspf_sudbury64_huc10_polygon.recharge_adjustment' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.72177,
      east_bound  =>  -71.25150,
      south_bound =>   42.19754,
      north_bound =>   42.64674
    },
    'hspf_sudbury64_huc12_point.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.65187,
      east_bound  =>  -71.30854,
      south_bound =>   42.25174,
      north_bound =>   42.57971
    },
    'hspf_sudbury64_huc12_point.recharge_adjustment' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.65187,
      east_bound  =>  -71.30854,
      south_bound =>   42.25174,
      north_bound =>   42.57971
    },
    'hspf_sudbury64_huc12_polygon.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.72177,
      east_bound  =>  -71.25150,
      south_bound =>   42.19754,
      north_bound =>   42.64674
    },
    'hspf_sudbury64_huc12_polygon.recharge_adjustment' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.72177,
      east_bound  =>  -71.25150,
      south_bound =>   42.19754,
      north_bound =>   42.64674
    },
    'hspf_sudbury65_huc10_point.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.48964,
      east_bound  =>  -71.46364,
      south_bound =>   42.31622,
      north_bound =>   42.44976
    },
    'hspf_sudbury65_huc10_point.recharge_adjustment' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.48964,
      east_bound  =>  -71.46364,
      south_bound =>   42.31622,
      north_bound =>   42.44976
    },
    'hspf_sudbury65_huc10_polygon.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.72177,
      east_bound  =>  -71.25150,
      south_bound =>   42.19754,
      north_bound =>   42.64674
    },
    'hspf_sudbury65_huc10_polygon.recharge_adjustment' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.72177,
      east_bound  =>  -71.25150,
      south_bound =>   42.19754,
      north_bound =>   42.64674
    },
    'hspf_sudbury65_huc12_point.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.65187,
      east_bound  =>  -71.30854,
      south_bound =>   42.25174,
      north_bound =>   42.57971
    },
    'hspf_sudbury65_huc12_point.recharge_adjustment' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.65187,
      east_bound  =>  -71.30854,
      south_bound =>   42.25174,
      north_bound =>   42.57971
    },
    'hspf_sudbury65_huc12_polygon.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.72177,
      east_bound  =>  -71.25150,
      south_bound =>   42.19754,
      north_bound =>   42.64674
    },
    'hspf_sudbury65_huc12_polygon.recharge_adjustment' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.72177,
      east_bound  =>  -71.25150,
      south_bound =>   42.19754,
      north_bound =>   42.64674
    },
    'hspf_sudbury66_huc10_point.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.48964,
      east_bound  =>  -71.46364,
      south_bound =>   42.31622,
      north_bound =>   42.44976
    },
    'hspf_sudbury66_huc10_point.recharge_adjustment' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.48964,
      east_bound  =>  -71.46364,
      south_bound =>   42.31622,
      north_bound =>   42.44976
    },
    'hspf_sudbury66_huc10_polygon.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.72177,
      east_bound  =>  -71.25150,
      south_bound =>   42.19754,
      north_bound =>   42.64674
    },
    'hspf_sudbury66_huc10_polygon.recharge_adjustment' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.72177,
      east_bound  =>  -71.25150,
      south_bound =>   42.19754,
      north_bound =>   42.64674
    },
    'hspf_sudbury66_huc12_point.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.65187,
      east_bound  =>  -71.30854,
      south_bound =>   42.25174,
      north_bound =>   42.57971
    },
    'hspf_sudbury66_huc12_point.recharge_adjustment' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.65187,
      east_bound  =>  -71.30854,
      south_bound =>   42.25174,
      north_bound =>   42.57971
    },
    'hspf_sudbury66_huc12_polygon.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.72177,
      east_bound  =>  -71.25150,
      south_bound =>   42.19754,
      north_bound =>   42.64674
    },
    'hspf_sudbury66_huc12_polygon.recharge_adjustment' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.72177,
      east_bound  =>  -71.25150,
      south_bound =>   42.19754,
      north_bound =>   42.64674
    },
    'hspf_sudbury67_huc10_point.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.48964,
      east_bound  =>  -71.46364,
      south_bound =>   42.31622,
      north_bound =>   42.44976
    },
    'hspf_sudbury67_huc10_point.recharge_adjustment' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.48964,
      east_bound  =>  -71.46364,
      south_bound =>   42.31622,
      north_bound =>   42.44976
    },
    'hspf_sudbury67_huc10_polygon.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.72177,
      east_bound  =>  -71.25150,
      south_bound =>   42.19754,
      north_bound =>   42.64674
    },
    'hspf_sudbury67_huc10_polygon.recharge_adjustment' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.72177,
      east_bound  =>  -71.25150,
      south_bound =>   42.19754,
      north_bound =>   42.64674
    },
    'hspf_sudbury67_huc12_point.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.65187,
      east_bound  =>  -71.30854,
      south_bound =>   42.25174,
      north_bound =>   42.57971
    },
    'hspf_sudbury67_huc12_point.recharge_adjustment' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.65187,
      east_bound  =>  -71.30854,
      south_bound =>   42.25174,
      north_bound =>   42.57971
    },
    'hspf_sudbury67_huc12_polygon.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.72177,
      east_bound  =>  -71.25150,
      south_bound =>   42.19754,
      north_bound =>   42.64674
    },
    'hspf_sudbury67_huc12_polygon.recharge_adjustment' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.72177,
      east_bound  =>  -71.25150,
      south_bound =>   42.19754,
      north_bound =>   42.64674
    },
    'hspf_sudbury68_huc10_point.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.48964,
      east_bound  =>  -71.46364,
      south_bound =>   42.31622,
      north_bound =>   42.44976
    },
    'hspf_sudbury68_huc10_point.recharge_adjustment' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.48964,
      east_bound  =>  -71.46364,
      south_bound =>   42.31622,
      north_bound =>   42.44976
    },
    'hspf_sudbury68_huc10_polygon.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.72177,
      east_bound  =>  -71.25150,
      south_bound =>   42.19754,
      north_bound =>   42.64674
    },
    'hspf_sudbury68_huc10_polygon.recharge_adjustment' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.72177,
      east_bound  =>  -71.25150,
      south_bound =>   42.19754,
      north_bound =>   42.64674
    },
    'hspf_sudbury68_huc12_point.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.65187,
      east_bound  =>  -71.30854,
      south_bound =>   42.25174,
      north_bound =>   42.57971
    },
    'hspf_sudbury68_huc12_point.recharge_adjustment' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.65187,
      east_bound  =>  -71.30854,
      south_bound =>   42.25174,
      north_bound =>   42.57971
    },
    'hspf_sudbury68_huc12_polygon.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.72177,
      east_bound  =>  -71.25150,
      south_bound =>   42.19754,
      north_bound =>   42.64674
    },
    'hspf_sudbury68_huc12_polygon.recharge_adjustment' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.72177,
      east_bound  =>  -71.25150,
      south_bound =>   42.19754,
      north_bound =>   42.64674
    },
    'hspf_sudbury7_huc10_point.hydrology' => {
      begin_date  => '1960-01-01',
      end_date    => '2009-12-31',
      west_bound  =>  -71.48964,
      east_bound  =>  -71.46364,
      south_bound =>   42.31622,
      north_bound =>   42.44976
    },
    'hspf_sudbury7_huc10_point.recharge_adjustment' => {
      begin_date  => '1960-01-01',
      end_date    => '2009-12-31',
      west_bound  =>  -71.48964,
      east_bound  =>  -71.46364,
      south_bound =>   42.31622,
      north_bound =>   42.44976
    },
    'hspf_sudbury7_huc10_polygon.hydrology' => {
      begin_date  => '1960-01-01',
      end_date    => '2009-12-31',
      west_bound  =>  -71.72177,
      east_bound  =>  -71.25150,
      south_bound =>   42.19754,
      north_bound =>   42.64674
    },
    'hspf_sudbury7_huc10_polygon.recharge_adjustment' => {
      begin_date  => '1960-01-01',
      end_date    => '2009-12-31',
      west_bound  =>  -71.72177,
      east_bound  =>  -71.25150,
      south_bound =>   42.19754,
      north_bound =>   42.64674
    },
    'hspf_sudbury7_huc12_point.hydrology' => {
      begin_date  => '1960-01-01',
      end_date    => '2009-12-31',
      west_bound  =>  -71.65187,
      east_bound  =>  -71.30854,
      south_bound =>   42.25174,
      north_bound =>   42.57971
    },
    'hspf_sudbury7_huc12_point.recharge_adjustment' => {
      begin_date  => '1960-01-01',
      end_date    => '2009-12-31',
      west_bound  =>  -71.65187,
      east_bound  =>  -71.30854,
      south_bound =>   42.25174,
      north_bound =>   42.57971
    },
    'hspf_sudbury7_huc12_polygon.hydrology' => {
      begin_date  => '1960-01-01',
      end_date    => '2009-12-31',
      west_bound  =>  -71.72177,
      east_bound  =>  -71.25150,
      south_bound =>   42.19754,
      north_bound =>   42.64674
    },
    'hspf_sudbury7_huc12_polygon.recharge_adjustment' => {
      begin_date  => '1960-01-01',
      end_date    => '2009-12-31',
      west_bound  =>  -71.72177,
      east_bound  =>  -71.25150,
      south_bound =>   42.19754,
      north_bound =>   42.64674
    },
    'hspf_taunton16_huc10_point.hydrology' => {
      begin_date  => '1960-01-01',
      end_date    => '2009-12-31',
      west_bound  =>  -71.22375,
      east_bound  =>  -70.92933,
      south_bound =>   41.87677,
      north_bound =>   42.05023
    },
    'hspf_taunton16_huc10_point.recharge_adjustment' => {
      begin_date  => '1960-01-01',
      end_date    => '2009-12-31',
      west_bound  =>  -71.22375,
      east_bound  =>  -70.92933,
      south_bound =>   41.87677,
      north_bound =>   42.05023
    },
    'hspf_taunton16_huc10_polygon.hydrology' => {
      begin_date  => '1960-01-01',
      end_date    => '2009-12-31',
      west_bound  =>  -71.34374,
      east_bound  =>  -70.71200,
      south_bound =>   41.70641,
      north_bound =>   42.16041
    },
    'hspf_taunton16_huc10_polygon.recharge_adjustment' => {
      begin_date  => '1960-01-01',
      end_date    => '2009-12-31',
      west_bound  =>  -71.34374,
      east_bound  =>  -70.71200,
      south_bound =>   41.70641,
      north_bound =>   42.16041
    },
    'hspf_taunton16_huc12_point.hydrology' => {
      begin_date  => '1960-01-01',
      end_date    => '2009-12-31',
      west_bound  =>  -71.26141,
      east_bound  =>  -70.84823,
      south_bound =>   41.79535,
      north_bound =>   42.08117
    },
    'hspf_taunton16_huc12_point.recharge_adjustment' => {
      begin_date  => '1960-01-01',
      end_date    => '2009-12-31',
      west_bound  =>  -71.26141,
      east_bound  =>  -70.84823,
      south_bound =>   41.79535,
      north_bound =>   42.08117
    },
    'hspf_taunton16_huc12_polygon.hydrology' => {
      begin_date  => '1960-01-01',
      end_date    => '2009-12-31',
      west_bound  =>  -71.34374,
      east_bound  =>  -70.71200,
      south_bound =>   41.70641,
      north_bound =>   42.16041
    },
    'hspf_taunton16_huc12_polygon.recharge_adjustment' => {
      begin_date  => '1960-01-01',
      end_date    => '2009-12-31',
      west_bound  =>  -71.34374,
      east_bound  =>  -70.71200,
      south_bound =>   41.70641,
      north_bound =>   42.16041
    },
    'hspf_taunton69_huc10_point.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.22375,
      east_bound  =>  -70.92933,
      south_bound =>   41.87677,
      north_bound =>   42.05023
    },
    'hspf_taunton69_huc10_point.recharge_adjustment' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.22375,
      east_bound  =>  -70.92933,
      south_bound =>   41.87677,
      north_bound =>   42.05023
    },
    'hspf_taunton69_huc10_polygon.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.34374,
      east_bound  =>  -70.71200,
      south_bound =>   41.70641,
      north_bound =>   42.16041
    },
    'hspf_taunton69_huc10_polygon.recharge_adjustment' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.34374,
      east_bound  =>  -70.71200,
      south_bound =>   41.70641,
      north_bound =>   42.16041
    },
    'hspf_taunton69_huc12_point.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.26141,
      east_bound  =>  -70.84823,
      south_bound =>   41.79535,
      north_bound =>   42.08117
    },
    'hspf_taunton69_huc12_point.recharge_adjustment' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.26141,
      east_bound  =>  -70.84823,
      south_bound =>   41.79535,
      north_bound =>   42.08117
    },
    'hspf_taunton69_huc12_polygon.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.34374,
      east_bound  =>  -70.71200,
      south_bound =>   41.70641,
      north_bound =>   42.16041
    },
    'hspf_taunton69_huc12_polygon.recharge_adjustment' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.34374,
      east_bound  =>  -70.71200,
      south_bound =>   41.70641,
      north_bound =>   42.16041
    },
    'hspf_taunton70_huc10_point.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.22375,
      east_bound  =>  -70.92933,
      south_bound =>   41.87677,
      north_bound =>   42.05023
    },
    'hspf_taunton70_huc10_point.recharge_adjustment' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.22375,
      east_bound  =>  -70.92933,
      south_bound =>   41.87677,
      north_bound =>   42.05023
    },
    'hspf_taunton70_huc10_polygon.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.34374,
      east_bound  =>  -70.71200,
      south_bound =>   41.70641,
      north_bound =>   42.16041
    },
    'hspf_taunton70_huc10_polygon.recharge_adjustment' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.34374,
      east_bound  =>  -70.71200,
      south_bound =>   41.70641,
      north_bound =>   42.16041
    },
    'hspf_taunton70_huc12_point.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.26141,
      east_bound  =>  -70.84823,
      south_bound =>   41.79535,
      north_bound =>   42.08117
    },
    'hspf_taunton70_huc12_point.recharge_adjustment' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.26141,
      east_bound  =>  -70.84823,
      south_bound =>   41.79535,
      north_bound =>   42.08117
    },
    'hspf_taunton70_huc12_polygon.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.34374,
      east_bound  =>  -70.71200,
      south_bound =>   41.70641,
      north_bound =>   42.16041
    },
    'hspf_taunton70_huc12_polygon.recharge_adjustment' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.34374,
      east_bound  =>  -70.71200,
      south_bound =>   41.70641,
      north_bound =>   42.16041
    },
    'hspf_taunton71_huc10_point.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.22375,
      east_bound  =>  -70.92933,
      south_bound =>   41.87677,
      north_bound =>   42.05023
    },
    'hspf_taunton71_huc10_point.recharge_adjustment' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.22375,
      east_bound  =>  -70.92933,
      south_bound =>   41.87677,
      north_bound =>   42.05023
    },
    'hspf_taunton71_huc10_polygon.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.34374,
      east_bound  =>  -70.71200,
      south_bound =>   41.70641,
      north_bound =>   42.16041
    },
    'hspf_taunton71_huc10_polygon.recharge_adjustment' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.34374,
      east_bound  =>  -70.71200,
      south_bound =>   41.70641,
      north_bound =>   42.16041
    },
    'hspf_taunton71_huc12_point.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.26141,
      east_bound  =>  -70.84823,
      south_bound =>   41.79535,
      north_bound =>   42.08117
    },
    'hspf_taunton71_huc12_point.recharge_adjustment' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.26141,
      east_bound  =>  -70.84823,
      south_bound =>   41.79535,
      north_bound =>   42.08117
    },
    'hspf_taunton71_huc12_polygon.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.34374,
      east_bound  =>  -70.71200,
      south_bound =>   41.70641,
      north_bound =>   42.16041
    },
    'hspf_taunton71_huc12_polygon.recharge_adjustment' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.34374,
      east_bound  =>  -70.71200,
      south_bound =>   41.70641,
      north_bound =>   42.16041
    },
    'hspf_taunton72_huc10_point.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.22375,
      east_bound  =>  -70.92933,
      south_bound =>   41.87677,
      north_bound =>   42.05023
    },
    'hspf_taunton72_huc10_point.recharge_adjustment' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.22375,
      east_bound  =>  -70.92933,
      south_bound =>   41.87677,
      north_bound =>   42.05023
    },
    'hspf_taunton72_huc10_polygon.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.34374,
      east_bound  =>  -70.71200,
      south_bound =>   41.70641,
      north_bound =>   42.16041
    },
    'hspf_taunton72_huc10_polygon.recharge_adjustment' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.34374,
      east_bound  =>  -70.71200,
      south_bound =>   41.70641,
      north_bound =>   42.16041
    },
    'hspf_taunton72_huc12_point.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.26141,
      east_bound  =>  -70.84823,
      south_bound =>   41.79535,
      north_bound =>   42.08117
    },
    'hspf_taunton72_huc12_point.recharge_adjustment' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.26141,
      east_bound  =>  -70.84823,
      south_bound =>   41.79535,
      north_bound =>   42.08117
    },
    'hspf_taunton72_huc12_polygon.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.34374,
      east_bound  =>  -70.71200,
      south_bound =>   41.70641,
      north_bound =>   42.16041
    },
    'hspf_taunton72_huc12_polygon.recharge_adjustment' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.34374,
      east_bound  =>  -70.71200,
      south_bound =>   41.70641,
      north_bound =>   42.16041
    },
    'hspf_taunton73_huc10_point.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.22375,
      east_bound  =>  -70.92933,
      south_bound =>   41.87677,
      north_bound =>   42.05023
    },
    'hspf_taunton73_huc10_point.recharge_adjustment' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.22375,
      east_bound  =>  -70.92933,
      south_bound =>   41.87677,
      north_bound =>   42.05023
    },
    'hspf_taunton73_huc10_polygon.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.34374,
      east_bound  =>  -70.71200,
      south_bound =>   41.70641,
      north_bound =>   42.16041
    },
    'hspf_taunton73_huc10_polygon.recharge_adjustment' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.34374,
      east_bound  =>  -70.71200,
      south_bound =>   41.70641,
      north_bound =>   42.16041
    },
    'hspf_taunton73_huc12_point.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.26141,
      east_bound  =>  -70.84823,
      south_bound =>   41.79535,
      north_bound =>   42.08117
    },
    'hspf_taunton73_huc12_point.recharge_adjustment' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.26141,
      east_bound  =>  -70.84823,
      south_bound =>   41.79535,
      north_bound =>   42.08117
    },
    'hspf_taunton73_huc12_polygon.hydrology' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.34374,
      east_bound  =>  -70.71200,
      south_bound =>   41.70641,
      north_bound =>   42.16041
    },
    'hspf_taunton73_huc12_polygon.recharge_adjustment' => {
      begin_date  => '1975-01-01',
      end_date    => '2004-12-31',
      west_bound  =>  -71.34374,
      east_bound  =>  -70.71200,
      south_bound =>   41.70641,
      north_bound =>   42.16041
    },
    'hspf_taunton8_huc10_point.hydrology' => {
      begin_date  => '1960-01-01',
      end_date    => '2009-12-31',
      west_bound  =>  -71.22375,
      east_bound  =>  -70.92933,
      south_bound =>   41.87677,
      north_bound =>   42.05023
    },
    'hspf_taunton8_huc10_point.recharge_adjustment' => {
      begin_date  => '1960-01-01',
      end_date    => '2009-12-31',
      west_bound  =>  -71.22375,
      east_bound  =>  -70.92933,
      south_bound =>   41.87677,
      north_bound =>   42.05023
    },
    'hspf_taunton8_huc10_polygon.hydrology' => {
      begin_date  => '1960-01-01',
      end_date    => '2009-12-31',
      west_bound  =>  -71.34374,
      east_bound  =>  -70.71200,
      south_bound =>   41.70641,
      north_bound =>   42.16041
    },
    'hspf_taunton8_huc10_polygon.recharge_adjustment' => {
      begin_date  => '1960-01-01',
      end_date    => '2009-12-31',
      west_bound  =>  -71.34374,
      east_bound  =>  -70.71200,
      south_bound =>   41.70641,
      north_bound =>   42.16041
    },
    'hspf_taunton8_huc12_point.hydrology' => {
      begin_date  => '1960-01-01',
      end_date    => '2009-12-31',
      west_bound  =>  -71.26141,
      east_bound  =>  -70.84823,
      south_bound =>   41.79535,
      north_bound =>   42.08117
    },
    'hspf_taunton8_huc12_point.recharge_adjustment' => {
      begin_date  => '1960-01-01',
      end_date    => '2009-12-31',
      west_bound  =>  -71.26141,
      east_bound  =>  -70.84823,
      south_bound =>   41.79535,
      north_bound =>   42.08117
    },
    'hspf_taunton8_huc12_polygon.hydrology' => {
      begin_date  => '1960-01-01',
      end_date    => '2009-12-31',
      west_bound  =>  -71.34374,
      east_bound  =>  -70.71200,
      south_bound =>   41.70641,
      north_bound =>   42.16041
    },
    'hspf_taunton8_huc12_polygon.recharge_adjustment' => {
      begin_date  => '1960-01-01',
      end_date    => '2009-12-31',
      west_bound  =>  -71.34374,
      east_bound  =>  -70.71200,
      south_bound =>   41.70641,
      north_bound =>   42.16041
    },
    'swmm_blackstone17_huc10_point.hydrology' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.76534,
      east_bound  =>  -71.76534,
      south_bound =>   42.19291,
      north_bound =>   42.19291
    },
    'swmm_blackstone17_huc10_point.loadings_tn' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.76534,
      east_bound  =>  -71.76534,
      south_bound =>   42.19291,
      north_bound =>   42.19291
    },
    'swmm_blackstone17_huc10_point.loadings_tp' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.76534,
      east_bound  =>  -71.76534,
      south_bound =>   42.19291,
      north_bound =>   42.19291
    },
    'swmm_blackstone17_huc10_point.loadings_tss' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.76534,
      east_bound  =>  -71.76534,
      south_bound =>   42.19291,
      north_bound =>   42.19291
    },
    'swmm_blackstone17_huc10_point.loadings_zn' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.76534,
      east_bound  =>  -71.76534,
      south_bound =>   42.19291,
      north_bound =>   42.19291
    },
    'swmm_blackstone17_huc10_polygon.hydrology' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.92799,
      east_bound  =>  -71.60905,
      south_bound =>   42.03787,
      north_bound =>   42.35075
    },
    'swmm_blackstone17_huc10_polygon.loadings_tn' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.92799,
      east_bound  =>  -71.60905,
      south_bound =>   42.03787,
      north_bound =>   42.35075
    },
    'swmm_blackstone17_huc10_polygon.loadings_tp' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.92799,
      east_bound  =>  -71.60905,
      south_bound =>   42.03787,
      north_bound =>   42.35075
    },
    'swmm_blackstone17_huc10_polygon.loadings_tss' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.92799,
      east_bound  =>  -71.60905,
      south_bound =>   42.03787,
      north_bound =>   42.35075
    },
    'swmm_blackstone17_huc10_polygon.loadings_zn' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.92799,
      east_bound  =>  -71.60905,
      south_bound =>   42.03787,
      north_bound =>   42.35075
    },
    'swmm_blackstone17_huc12_point.hydrology' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.85803,
      east_bound  =>  -71.71627,
      south_bound =>   42.09478,
      north_bound =>   42.26947
    },
    'swmm_blackstone17_huc12_point.loadings_tn' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.85803,
      east_bound  =>  -71.71627,
      south_bound =>   42.09478,
      north_bound =>   42.26947
    },
    'swmm_blackstone17_huc12_point.loadings_tp' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.85803,
      east_bound  =>  -71.71627,
      south_bound =>   42.09478,
      north_bound =>   42.26947
    },
    'swmm_blackstone17_huc12_point.loadings_tss' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.85803,
      east_bound  =>  -71.71627,
      south_bound =>   42.09478,
      north_bound =>   42.26947
    },
    'swmm_blackstone17_huc12_point.loadings_zn' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.85803,
      east_bound  =>  -71.71627,
      south_bound =>   42.09478,
      north_bound =>   42.26947
    },
    'swmm_blackstone17_huc12_polygon.hydrology' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.92799,
      east_bound  =>  -71.60905,
      south_bound =>   42.03787,
      north_bound =>   42.35075
    },
    'swmm_blackstone17_huc12_polygon.loadings_tn' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.92799,
      east_bound  =>  -71.60905,
      south_bound =>   42.03787,
      north_bound =>   42.35075
    },
    'swmm_blackstone17_huc12_polygon.loadings_tp' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.92799,
      east_bound  =>  -71.60905,
      south_bound =>   42.03787,
      north_bound =>   42.35075
    },
    'swmm_blackstone17_huc12_polygon.loadings_tss' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.92799,
      east_bound  =>  -71.60905,
      south_bound =>   42.03787,
      north_bound =>   42.35075
    },
    'swmm_blackstone17_huc12_polygon.loadings_zn' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.92799,
      east_bound  =>  -71.60905,
      south_bound =>   42.03787,
      north_bound =>   42.35075
    },
    'swmm_blackstone18_huc10_point.hydrology' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.57708,
      east_bound  =>  -71.57708,
      south_bound =>   42.01962,
      north_bound =>   42.01962
    },
    'swmm_blackstone18_huc10_point.loadings_tn' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.57708,
      east_bound  =>  -71.57708,
      south_bound =>   42.01962,
      north_bound =>   42.01962
    },
    'swmm_blackstone18_huc10_point.loadings_tp' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.57708,
      east_bound  =>  -71.57708,
      south_bound =>   42.01962,
      north_bound =>   42.01962
    },
    'swmm_blackstone18_huc10_point.loadings_tss' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.57708,
      east_bound  =>  -71.57708,
      south_bound =>   42.01962,
      north_bound =>   42.01962
    },
    'swmm_blackstone18_huc10_point.loadings_zn' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.57708,
      east_bound  =>  -71.57708,
      south_bound =>   42.01962,
      north_bound =>   42.01962
    },
    'swmm_blackstone18_huc10_polygon.hydrology' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.78048,
      east_bound  =>  -71.35340,
      south_bound =>   41.85489,
      north_bound =>   42.23843
    },
    'swmm_blackstone18_huc10_polygon.loadings_tn' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.78048,
      east_bound  =>  -71.35340,
      south_bound =>   41.85489,
      north_bound =>   42.23843
    },
    'swmm_blackstone18_huc10_polygon.loadings_tp' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.78048,
      east_bound  =>  -71.35340,
      south_bound =>   41.85489,
      north_bound =>   42.23843
    },
    'swmm_blackstone18_huc10_polygon.loadings_tss' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.78048,
      east_bound  =>  -71.35340,
      south_bound =>   41.85489,
      north_bound =>   42.23843
    },
    'swmm_blackstone18_huc10_polygon.loadings_zn' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.78048,
      east_bound  =>  -71.35340,
      south_bound =>   41.85489,
      north_bound =>   42.23843
    },
    'swmm_blackstone18_huc12_point.hydrology' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.71881,
      east_bound  =>  -71.39789,
      south_bound =>   41.91013,
      north_bound =>   42.15551
    },
    'swmm_blackstone18_huc12_point.loadings_tn' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.71881,
      east_bound  =>  -71.39789,
      south_bound =>   41.91013,
      north_bound =>   42.15551
    },
    'swmm_blackstone18_huc12_point.loadings_tp' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.71881,
      east_bound  =>  -71.39789,
      south_bound =>   41.91013,
      north_bound =>   42.15551
    },
    'swmm_blackstone18_huc12_point.loadings_tss' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.71881,
      east_bound  =>  -71.39789,
      south_bound =>   41.91013,
      north_bound =>   42.15551
    },
    'swmm_blackstone18_huc12_point.loadings_zn' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.71881,
      east_bound  =>  -71.39789,
      south_bound =>   41.91013,
      north_bound =>   42.15551
    },
    'swmm_blackstone18_huc12_polygon.hydrology' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.78048,
      east_bound  =>  -71.35340,
      south_bound =>   41.85489,
      north_bound =>   42.23843
    },
    'swmm_blackstone18_huc12_polygon.loadings_tn' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.78048,
      east_bound  =>  -71.35340,
      south_bound =>   41.85489,
      north_bound =>   42.23843
    },
    'swmm_blackstone18_huc12_polygon.loadings_tp' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.78048,
      east_bound  =>  -71.35340,
      south_bound =>   41.85489,
      north_bound =>   42.23843
    },
    'swmm_blackstone18_huc12_polygon.loadings_tss' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.78048,
      east_bound  =>  -71.35340,
      south_bound =>   41.85489,
      north_bound =>   42.23843
    },
    'swmm_blackstone18_huc12_polygon.loadings_zn' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.78048,
      east_bound  =>  -71.35340,
      south_bound =>   41.85489,
      north_bound =>   42.23843
    },
    'swmm_blackstone25_huc10_point.hydrology' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.76534,
      east_bound  =>  -71.76534,
      south_bound =>   42.19291,
      north_bound =>   42.19291
    },
    'swmm_blackstone25_huc10_point.loadings_tn' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.76534,
      east_bound  =>  -71.76534,
      south_bound =>   42.19291,
      north_bound =>   42.19291
    },
    'swmm_blackstone25_huc10_point.loadings_tp' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.76534,
      east_bound  =>  -71.76534,
      south_bound =>   42.19291,
      north_bound =>   42.19291
    },
    'swmm_blackstone25_huc10_point.loadings_tss' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.76534,
      east_bound  =>  -71.76534,
      south_bound =>   42.19291,
      north_bound =>   42.19291
    },
    'swmm_blackstone25_huc10_point.loadings_zn' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.76534,
      east_bound  =>  -71.76534,
      south_bound =>   42.19291,
      north_bound =>   42.19291
    },
    'swmm_blackstone25_huc10_polygon.hydrology' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.92799,
      east_bound  =>  -71.60905,
      south_bound =>   42.03787,
      north_bound =>   42.35075
    },
    'swmm_blackstone25_huc10_polygon.loadings_tn' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.92799,
      east_bound  =>  -71.60905,
      south_bound =>   42.03787,
      north_bound =>   42.35075
    },
    'swmm_blackstone25_huc10_polygon.loadings_tp' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.92799,
      east_bound  =>  -71.60905,
      south_bound =>   42.03787,
      north_bound =>   42.35075
    },
    'swmm_blackstone25_huc10_polygon.loadings_tss' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.92799,
      east_bound  =>  -71.60905,
      south_bound =>   42.03787,
      north_bound =>   42.35075
    },
    'swmm_blackstone25_huc10_polygon.loadings_zn' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.92799,
      east_bound  =>  -71.60905,
      south_bound =>   42.03787,
      north_bound =>   42.35075
    },
    'swmm_blackstone25_huc12_point.hydrology' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.85803,
      east_bound  =>  -71.71627,
      south_bound =>   42.09478,
      north_bound =>   42.26947
    },
    'swmm_blackstone25_huc12_point.loadings_tn' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.85803,
      east_bound  =>  -71.71627,
      south_bound =>   42.09478,
      north_bound =>   42.26947
    },
    'swmm_blackstone25_huc12_point.loadings_tp' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.85803,
      east_bound  =>  -71.71627,
      south_bound =>   42.09478,
      north_bound =>   42.26947
    },
    'swmm_blackstone25_huc12_point.loadings_tss' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.85803,
      east_bound  =>  -71.71627,
      south_bound =>   42.09478,
      north_bound =>   42.26947
    },
    'swmm_blackstone25_huc12_point.loadings_zn' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.85803,
      east_bound  =>  -71.71627,
      south_bound =>   42.09478,
      north_bound =>   42.26947
    },
    'swmm_blackstone25_huc12_polygon.hydrology' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.92799,
      east_bound  =>  -71.60905,
      south_bound =>   42.03787,
      north_bound =>   42.35075
    },
    'swmm_blackstone25_huc12_polygon.loadings_tn' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.92799,
      east_bound  =>  -71.60905,
      south_bound =>   42.03787,
      north_bound =>   42.35075
    },
    'swmm_blackstone25_huc12_polygon.loadings_tp' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.92799,
      east_bound  =>  -71.60905,
      south_bound =>   42.03787,
      north_bound =>   42.35075
    },
    'swmm_blackstone25_huc12_polygon.loadings_tss' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.92799,
      east_bound  =>  -71.60905,
      south_bound =>   42.03787,
      north_bound =>   42.35075
    },
    'swmm_blackstone25_huc12_polygon.loadings_zn' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.92799,
      east_bound  =>  -71.60905,
      south_bound =>   42.03787,
      north_bound =>   42.35075
    },
    'swmm_blackstone26_huc10_point.hydrology' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.57708,
      east_bound  =>  -71.57708,
      south_bound =>   42.01962,
      north_bound =>   42.01962
    },
    'swmm_blackstone26_huc10_point.loadings_tn' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.57708,
      east_bound  =>  -71.57708,
      south_bound =>   42.01962,
      north_bound =>   42.01962
    },
    'swmm_blackstone26_huc10_point.loadings_tp' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.57708,
      east_bound  =>  -71.57708,
      south_bound =>   42.01962,
      north_bound =>   42.01962
    },
    'swmm_blackstone26_huc10_point.loadings_tss' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.57708,
      east_bound  =>  -71.57708,
      south_bound =>   42.01962,
      north_bound =>   42.01962
    },
    'swmm_blackstone26_huc10_point.loadings_zn' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.57708,
      east_bound  =>  -71.57708,
      south_bound =>   42.01962,
      north_bound =>   42.01962
    },
    'swmm_blackstone26_huc10_polygon.hydrology' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.78048,
      east_bound  =>  -71.35340,
      south_bound =>   41.85489,
      north_bound =>   42.23843
    },
    'swmm_blackstone26_huc10_polygon.loadings_tn' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.78048,
      east_bound  =>  -71.35340,
      south_bound =>   41.85489,
      north_bound =>   42.23843
    },
    'swmm_blackstone26_huc10_polygon.loadings_tp' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.78048,
      east_bound  =>  -71.35340,
      south_bound =>   41.85489,
      north_bound =>   42.23843
    },
    'swmm_blackstone26_huc10_polygon.loadings_tss' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.78048,
      east_bound  =>  -71.35340,
      south_bound =>   41.85489,
      north_bound =>   42.23843
    },
    'swmm_blackstone26_huc10_polygon.loadings_zn' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.78048,
      east_bound  =>  -71.35340,
      south_bound =>   41.85489,
      north_bound =>   42.23843
    },
    'swmm_blackstone26_huc12_point.hydrology' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.71881,
      east_bound  =>  -71.39789,
      south_bound =>   41.91013,
      north_bound =>   42.15551
    },
    'swmm_blackstone26_huc12_point.loadings_tn' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.71881,
      east_bound  =>  -71.39789,
      south_bound =>   41.91013,
      north_bound =>   42.15551
    },
    'swmm_blackstone26_huc12_point.loadings_tp' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.71881,
      east_bound  =>  -71.39789,
      south_bound =>   41.91013,
      north_bound =>   42.15551
    },
    'swmm_blackstone26_huc12_point.loadings_tss' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.71881,
      east_bound  =>  -71.39789,
      south_bound =>   41.91013,
      north_bound =>   42.15551
    },
    'swmm_blackstone26_huc12_point.loadings_zn' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.71881,
      east_bound  =>  -71.39789,
      south_bound =>   41.91013,
      north_bound =>   42.15551
    },
    'swmm_blackstone26_huc12_polygon.hydrology' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.78048,
      east_bound  =>  -71.35340,
      south_bound =>   41.85489,
      north_bound =>   42.23843
    },
    'swmm_blackstone26_huc12_polygon.loadings_tn' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.78048,
      east_bound  =>  -71.35340,
      south_bound =>   41.85489,
      north_bound =>   42.23843
    },
    'swmm_blackstone26_huc12_polygon.loadings_tp' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.78048,
      east_bound  =>  -71.35340,
      south_bound =>   41.85489,
      north_bound =>   42.23843
    },
    'swmm_blackstone26_huc12_polygon.loadings_tss' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.78048,
      east_bound  =>  -71.35340,
      south_bound =>   41.85489,
      north_bound =>   42.23843
    },
    'swmm_blackstone26_huc12_polygon.loadings_zn' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.78048,
      east_bound  =>  -71.35340,
      south_bound =>   41.85489,
      north_bound =>   42.23843
    },
    'swmm_charles19_huc10_point.hydrology' => {
      begin_date  => '1961-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.40021,
      east_bound  =>  -71.21325,
      south_bound =>   42.15069,
      north_bound =>   42.32819
    },
    'swmm_charles19_huc10_point.loadings_tn' => {
      begin_date  => '1961-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.40021,
      east_bound  =>  -71.21325,
      south_bound =>   42.15069,
      north_bound =>   42.32819
    },
    'swmm_charles19_huc10_point.loadings_tp' => {
      begin_date  => '1961-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.40021,
      east_bound  =>  -71.21325,
      south_bound =>   42.15069,
      north_bound =>   42.32819
    },
    'swmm_charles19_huc10_point.loadings_tss' => {
      begin_date  => '1961-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.40021,
      east_bound  =>  -71.21325,
      south_bound =>   42.15069,
      north_bound =>   42.32819
    },
    'swmm_charles19_huc10_point.loadings_zn' => {
      begin_date  => '1961-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.40021,
      east_bound  =>  -71.21325,
      south_bound =>   42.15069,
      north_bound =>   42.32819
    },
    'swmm_charles19_huc10_polygon.hydrology' => {
      begin_date  => '1961-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.54945,
      east_bound  =>  -71.00937,
      south_bound =>   42.02564,
      north_bound =>   42.45390
    },
    'swmm_charles19_huc10_polygon.loadings_tn' => {
      begin_date  => '1961-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.54945,
      east_bound  =>  -71.00937,
      south_bound =>   42.02564,
      north_bound =>   42.45390
    },
    'swmm_charles19_huc10_polygon.loadings_tp' => {
      begin_date  => '1961-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.54945,
      east_bound  =>  -71.00937,
      south_bound =>   42.02564,
      north_bound =>   42.45390
    },
    'swmm_charles19_huc10_polygon.loadings_tss' => {
      begin_date  => '1961-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.54945,
      east_bound  =>  -71.00937,
      south_bound =>   42.02564,
      north_bound =>   42.45390
    },
    'swmm_charles19_huc10_polygon.loadings_zn' => {
      begin_date  => '1961-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.54945,
      east_bound  =>  -71.00937,
      south_bound =>   42.02564,
      north_bound =>   42.45390
    },
    'swmm_charles19_huc12_point.hydrology' => {
      begin_date  => '1961-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.47076,
      east_bound  =>  -71.12661,
      south_bound =>   42.12082,
      north_bound =>   42.40220
    },
    'swmm_charles19_huc12_point.loadings_tn' => {
      begin_date  => '1961-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.47076,
      east_bound  =>  -71.12661,
      south_bound =>   42.12082,
      north_bound =>   42.40220
    },
    'swmm_charles19_huc12_point.loadings_tp' => {
      begin_date  => '1961-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.47076,
      east_bound  =>  -71.12661,
      south_bound =>   42.12082,
      north_bound =>   42.40220
    },
    'swmm_charles19_huc12_point.loadings_tss' => {
      begin_date  => '1961-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.47076,
      east_bound  =>  -71.12661,
      south_bound =>   42.12082,
      north_bound =>   42.40220
    },
    'swmm_charles19_huc12_point.loadings_zn' => {
      begin_date  => '1961-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.47076,
      east_bound  =>  -71.12661,
      south_bound =>   42.12082,
      north_bound =>   42.40220
    },
    'swmm_charles19_huc12_polygon.hydrology' => {
      begin_date  => '1961-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.54945,
      east_bound  =>  -71.00937,
      south_bound =>   42.02564,
      north_bound =>   42.45390
    },
    'swmm_charles19_huc12_polygon.loadings_tn' => {
      begin_date  => '1961-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.54945,
      east_bound  =>  -71.00937,
      south_bound =>   42.02564,
      north_bound =>   42.45390
    },
    'swmm_charles19_huc12_polygon.loadings_tp' => {
      begin_date  => '1961-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.54945,
      east_bound  =>  -71.00937,
      south_bound =>   42.02564,
      north_bound =>   42.45390
    },
    'swmm_charles19_huc12_polygon.loadings_tss' => {
      begin_date  => '1961-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.54945,
      east_bound  =>  -71.00937,
      south_bound =>   42.02564,
      north_bound =>   42.45390
    },
    'swmm_charles19_huc12_polygon.loadings_zn' => {
      begin_date  => '1961-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.54945,
      east_bound  =>  -71.00937,
      south_bound =>   42.02564,
      north_bound =>   42.45390
    },
    'swmm_charles27_huc10_point.hydrology' => {
      begin_date  => '1961-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.40021,
      east_bound  =>  -71.21325,
      south_bound =>   42.15069,
      north_bound =>   42.32819
    },
    'swmm_charles27_huc10_point.loadings_tn' => {
      begin_date  => '1961-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.40021,
      east_bound  =>  -71.21325,
      south_bound =>   42.15069,
      north_bound =>   42.32819
    },
    'swmm_charles27_huc10_point.loadings_tp' => {
      begin_date  => '1961-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.40021,
      east_bound  =>  -71.21325,
      south_bound =>   42.15069,
      north_bound =>   42.32819
    },
    'swmm_charles27_huc10_point.loadings_tss' => {
      begin_date  => '1961-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.40021,
      east_bound  =>  -71.21325,
      south_bound =>   42.15069,
      north_bound =>   42.32819
    },
    'swmm_charles27_huc10_point.loadings_zn' => {
      begin_date  => '1961-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.40021,
      east_bound  =>  -71.21325,
      south_bound =>   42.15069,
      north_bound =>   42.32819
    },
    'swmm_charles27_huc10_polygon.hydrology' => {
      begin_date  => '1961-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.54945,
      east_bound  =>  -71.00937,
      south_bound =>   42.02564,
      north_bound =>   42.45390
    },
    'swmm_charles27_huc10_polygon.loadings_tn' => {
      begin_date  => '1961-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.54945,
      east_bound  =>  -71.00937,
      south_bound =>   42.02564,
      north_bound =>   42.45390
    },
    'swmm_charles27_huc10_polygon.loadings_tp' => {
      begin_date  => '1961-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.54945,
      east_bound  =>  -71.00937,
      south_bound =>   42.02564,
      north_bound =>   42.45390
    },
    'swmm_charles27_huc10_polygon.loadings_tss' => {
      begin_date  => '1961-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.54945,
      east_bound  =>  -71.00937,
      south_bound =>   42.02564,
      north_bound =>   42.45390
    },
    'swmm_charles27_huc10_polygon.loadings_zn' => {
      begin_date  => '1961-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.54945,
      east_bound  =>  -71.00937,
      south_bound =>   42.02564,
      north_bound =>   42.45390
    },
    'swmm_charles27_huc12_point.hydrology' => {
      begin_date  => '1961-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.47076,
      east_bound  =>  -71.12661,
      south_bound =>   42.12082,
      north_bound =>   42.40220
    },
    'swmm_charles27_huc12_point.loadings_tn' => {
      begin_date  => '1961-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.47076,
      east_bound  =>  -71.12661,
      south_bound =>   42.12082,
      north_bound =>   42.40220
    },
    'swmm_charles27_huc12_point.loadings_tp' => {
      begin_date  => '1961-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.47076,
      east_bound  =>  -71.12661,
      south_bound =>   42.12082,
      north_bound =>   42.40220
    },
    'swmm_charles27_huc12_point.loadings_tss' => {
      begin_date  => '1961-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.47076,
      east_bound  =>  -71.12661,
      south_bound =>   42.12082,
      north_bound =>   42.40220
    },
    'swmm_charles27_huc12_point.loadings_zn' => {
      begin_date  => '1961-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.47076,
      east_bound  =>  -71.12661,
      south_bound =>   42.12082,
      north_bound =>   42.40220
    },
    'swmm_charles27_huc12_polygon.hydrology' => {
      begin_date  => '1961-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.54945,
      east_bound  =>  -71.00937,
      south_bound =>   42.02564,
      north_bound =>   42.45390
    },
    'swmm_charles27_huc12_polygon.loadings_tn' => {
      begin_date  => '1961-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.54945,
      east_bound  =>  -71.00937,
      south_bound =>   42.02564,
      north_bound =>   42.45390
    },
    'swmm_charles27_huc12_polygon.loadings_tp' => {
      begin_date  => '1961-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.54945,
      east_bound  =>  -71.00937,
      south_bound =>   42.02564,
      north_bound =>   42.45390
    },
    'swmm_charles27_huc12_polygon.loadings_tss' => {
      begin_date  => '1961-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.54945,
      east_bound  =>  -71.00937,
      south_bound =>   42.02564,
      north_bound =>   42.45390
    },
    'swmm_charles27_huc12_polygon.loadings_zn' => {
      begin_date  => '1961-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.54945,
      east_bound  =>  -71.00937,
      south_bound =>   42.02564,
      north_bound =>   42.45390
    },
    'swmm_farmington20_huc10_point.hydrology' => {
      begin_date  => '1970-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -73.11394,
      east_bound  =>  -72.77703,
      south_bound =>   41.77390,
      north_bound =>   42.10068
    },
    'swmm_farmington20_huc10_point.loadings_tn' => {
      begin_date  => '1970-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -73.11394,
      east_bound  =>  -72.77703,
      south_bound =>   41.77390,
      north_bound =>   42.10068
    },
    'swmm_farmington20_huc10_point.loadings_tp' => {
      begin_date  => '1970-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -73.11394,
      east_bound  =>  -72.77703,
      south_bound =>   41.77390,
      north_bound =>   42.10068
    },
    'swmm_farmington20_huc10_point.loadings_tss' => {
      begin_date  => '1970-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -73.11394,
      east_bound  =>  -72.77703,
      south_bound =>   41.77390,
      north_bound =>   42.10068
    },
    'swmm_farmington20_huc10_point.loadings_zn' => {
      begin_date  => '1970-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -73.11394,
      east_bound  =>  -72.77703,
      south_bound =>   41.77390,
      north_bound =>   42.10068
    },
    'swmm_farmington20_huc10_polygon.hydrology' => {
      begin_date  => '1970-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -73.20258,
      east_bound  =>  -72.62316,
      south_bound =>   41.63613,
      north_bound =>   42.30952
    },
    'swmm_farmington20_huc10_polygon.loadings_tn' => {
      begin_date  => '1970-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -73.20258,
      east_bound  =>  -72.62316,
      south_bound =>   41.63613,
      north_bound =>   42.30952
    },
    'swmm_farmington20_huc10_polygon.loadings_tp' => {
      begin_date  => '1970-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -73.20258,
      east_bound  =>  -72.62316,
      south_bound =>   41.63613,
      north_bound =>   42.30952
    },
    'swmm_farmington20_huc10_polygon.loadings_tss' => {
      begin_date  => '1970-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -73.20258,
      east_bound  =>  -72.62316,
      south_bound =>   41.63613,
      north_bound =>   42.30952
    },
    'swmm_farmington20_huc10_polygon.loadings_zn' => {
      begin_date  => '1970-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -73.20258,
      east_bound  =>  -72.62316,
      south_bound =>   41.63613,
      north_bound =>   42.30952
    },
    'swmm_farmington20_huc12_point.hydrology' => {
      begin_date  => '1970-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -73.12987,
      east_bound  =>  -72.75031,
      south_bound =>   41.69660,
      north_bound =>   42.19990
    },
    'swmm_farmington20_huc12_point.loadings_tn' => {
      begin_date  => '1970-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -73.12987,
      east_bound  =>  -72.75031,
      south_bound =>   41.69660,
      north_bound =>   42.19990
    },
    'swmm_farmington20_huc12_point.loadings_tp' => {
      begin_date  => '1970-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -73.12987,
      east_bound  =>  -72.75031,
      south_bound =>   41.69660,
      north_bound =>   42.19990
    },
    'swmm_farmington20_huc12_point.loadings_tss' => {
      begin_date  => '1970-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -73.12987,
      east_bound  =>  -72.75031,
      south_bound =>   41.69660,
      north_bound =>   42.19990
    },
    'swmm_farmington20_huc12_point.loadings_zn' => {
      begin_date  => '1970-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -73.12987,
      east_bound  =>  -72.75031,
      south_bound =>   41.69660,
      north_bound =>   42.19990
    },
    'swmm_farmington20_huc12_polygon.hydrology' => {
      begin_date  => '1970-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -73.20258,
      east_bound  =>  -72.62316,
      south_bound =>   41.63613,
      north_bound =>   42.30952
    },
    'swmm_farmington20_huc12_polygon.loadings_tn' => {
      begin_date  => '1970-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -73.20258,
      east_bound  =>  -72.62316,
      south_bound =>   41.63613,
      north_bound =>   42.30952
    },
    'swmm_farmington20_huc12_polygon.loadings_tp' => {
      begin_date  => '1970-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -73.20258,
      east_bound  =>  -72.62316,
      south_bound =>   41.63613,
      north_bound =>   42.30952
    },
    'swmm_farmington20_huc12_polygon.loadings_tss' => {
      begin_date  => '1970-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -73.20258,
      east_bound  =>  -72.62316,
      south_bound =>   41.63613,
      north_bound =>   42.30952
    },
    'swmm_farmington20_huc12_polygon.loadings_zn' => {
      begin_date  => '1970-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -73.20258,
      east_bound  =>  -72.62316,
      south_bound =>   41.63613,
      north_bound =>   42.30952
    },
    'swmm_farmington28_huc10_point.hydrology' => {
      begin_date  => '1970-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -73.11394,
      east_bound  =>  -72.77703,
      south_bound =>   41.77390,
      north_bound =>   42.10068
    },
    'swmm_farmington28_huc10_point.loadings_tn' => {
      begin_date  => '1970-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -73.11394,
      east_bound  =>  -72.77703,
      south_bound =>   41.77390,
      north_bound =>   42.10068
    },
    'swmm_farmington28_huc10_point.loadings_tp' => {
      begin_date  => '1970-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -73.11394,
      east_bound  =>  -72.77703,
      south_bound =>   41.77390,
      north_bound =>   42.10068
    },
    'swmm_farmington28_huc10_point.loadings_tss' => {
      begin_date  => '1970-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -73.11394,
      east_bound  =>  -72.77703,
      south_bound =>   41.77390,
      north_bound =>   42.10068
    },
    'swmm_farmington28_huc10_point.loadings_zn' => {
      begin_date  => '1970-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -73.11394,
      east_bound  =>  -72.77703,
      south_bound =>   41.77390,
      north_bound =>   42.10068
    },
    'swmm_farmington28_huc10_polygon.hydrology' => {
      begin_date  => '1970-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -73.20258,
      east_bound  =>  -72.62316,
      south_bound =>   41.63613,
      north_bound =>   42.30952
    },
    'swmm_farmington28_huc10_polygon.loadings_tn' => {
      begin_date  => '1970-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -73.20258,
      east_bound  =>  -72.62316,
      south_bound =>   41.63613,
      north_bound =>   42.30952
    },
    'swmm_farmington28_huc10_polygon.loadings_tp' => {
      begin_date  => '1970-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -73.20258,
      east_bound  =>  -72.62316,
      south_bound =>   41.63613,
      north_bound =>   42.30952
    },
    'swmm_farmington28_huc10_polygon.loadings_tss' => {
      begin_date  => '1970-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -73.20258,
      east_bound  =>  -72.62316,
      south_bound =>   41.63613,
      north_bound =>   42.30952
    },
    'swmm_farmington28_huc10_polygon.loadings_zn' => {
      begin_date  => '1970-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -73.20258,
      east_bound  =>  -72.62316,
      south_bound =>   41.63613,
      north_bound =>   42.30952
    },
    'swmm_farmington28_huc12_point.hydrology' => {
      begin_date  => '1970-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -73.12987,
      east_bound  =>  -72.75031,
      south_bound =>   41.69660,
      north_bound =>   42.19990
    },
    'swmm_farmington28_huc12_point.loadings_tn' => {
      begin_date  => '1970-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -73.12987,
      east_bound  =>  -72.75031,
      south_bound =>   41.69660,
      north_bound =>   42.19990
    },
    'swmm_farmington28_huc12_point.loadings_tp' => {
      begin_date  => '1970-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -73.12987,
      east_bound  =>  -72.75031,
      south_bound =>   41.69660,
      north_bound =>   42.19990
    },
    'swmm_farmington28_huc12_point.loadings_tss' => {
      begin_date  => '1970-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -73.12987,
      east_bound  =>  -72.75031,
      south_bound =>   41.69660,
      north_bound =>   42.19990
    },
    'swmm_farmington28_huc12_point.loadings_zn' => {
      begin_date  => '1970-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -73.12987,
      east_bound  =>  -72.75031,
      south_bound =>   41.69660,
      north_bound =>   42.19990
    },
    'swmm_farmington28_huc12_polygon.hydrology' => {
      begin_date  => '1970-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -73.20258,
      east_bound  =>  -72.62316,
      south_bound =>   41.63613,
      north_bound =>   42.30952
    },
    'swmm_farmington28_huc12_polygon.loadings_tn' => {
      begin_date  => '1970-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -73.20258,
      east_bound  =>  -72.62316,
      south_bound =>   41.63613,
      north_bound =>   42.30952
    },
    'swmm_farmington28_huc12_polygon.loadings_tp' => {
      begin_date  => '1970-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -73.20258,
      east_bound  =>  -72.62316,
      south_bound =>   41.63613,
      north_bound =>   42.30952
    },
    'swmm_farmington28_huc12_polygon.loadings_tss' => {
      begin_date  => '1970-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -73.20258,
      east_bound  =>  -72.62316,
      south_bound =>   41.63613,
      north_bound =>   42.30952
    },
    'swmm_farmington28_huc12_polygon.loadings_zn' => {
      begin_date  => '1970-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -73.20258,
      east_bound  =>  -72.62316,
      south_bound =>   41.63613,
      north_bound =>   42.30952
    },
    'swmm_ipswich21_huc10_point.hydrology' => {
      begin_date  => '1961-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.13159,
      east_bound  =>  -71.13159,
      south_bound =>   42.56714,
      north_bound =>   42.56714
    },
    'swmm_ipswich21_huc10_point.loadings_tn' => {
      begin_date  => '1961-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.13159,
      east_bound  =>  -71.13159,
      south_bound =>   42.56714,
      north_bound =>   42.56714
    },
    'swmm_ipswich21_huc10_point.loadings_tp' => {
      begin_date  => '1961-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.13159,
      east_bound  =>  -71.13159,
      south_bound =>   42.56714,
      north_bound =>   42.56714
    },
    'swmm_ipswich21_huc10_point.loadings_tss' => {
      begin_date  => '1961-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.13159,
      east_bound  =>  -71.13159,
      south_bound =>   42.56714,
      north_bound =>   42.56714
    },
    'swmm_ipswich21_huc10_point.loadings_zn' => {
      begin_date  => '1961-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.13159,
      east_bound  =>  -71.13159,
      south_bound =>   42.56714,
      north_bound =>   42.56714
    },
    'swmm_ipswich21_huc10_polygon.hydrology' => {
      begin_date  => '1961-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.21963,
      east_bound  =>  -70.81332,
      south_bound =>   42.49842,
      north_bound =>   42.70367
    },
    'swmm_ipswich21_huc10_polygon.loadings_tn' => {
      begin_date  => '1961-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.21963,
      east_bound  =>  -70.81332,
      south_bound =>   42.49842,
      north_bound =>   42.70367
    },
    'swmm_ipswich21_huc10_polygon.loadings_tp' => {
      begin_date  => '1961-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.21963,
      east_bound  =>  -70.81332,
      south_bound =>   42.49842,
      north_bound =>   42.70367
    },
    'swmm_ipswich21_huc10_polygon.loadings_tss' => {
      begin_date  => '1961-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.21963,
      east_bound  =>  -70.81332,
      south_bound =>   42.49842,
      north_bound =>   42.70367
    },
    'swmm_ipswich21_huc10_polygon.loadings_zn' => {
      begin_date  => '1961-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.21963,
      east_bound  =>  -70.81332,
      south_bound =>   42.49842,
      north_bound =>   42.70367
    },
    'swmm_ipswich21_huc12_point.hydrology' => {
      begin_date  => '1961-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.13159,
      east_bound  =>  -70.89173,
      south_bound =>   42.56714,
      north_bound =>   42.67464
    },
    'swmm_ipswich21_huc12_point.loadings_tn' => {
      begin_date  => '1961-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.13159,
      east_bound  =>  -70.89173,
      south_bound =>   42.56714,
      north_bound =>   42.67464
    },
    'swmm_ipswich21_huc12_point.loadings_tp' => {
      begin_date  => '1961-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.13159,
      east_bound  =>  -70.89173,
      south_bound =>   42.56714,
      north_bound =>   42.67464
    },
    'swmm_ipswich21_huc12_point.loadings_tss' => {
      begin_date  => '1961-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.13159,
      east_bound  =>  -70.89173,
      south_bound =>   42.56714,
      north_bound =>   42.67464
    },
    'swmm_ipswich21_huc12_point.loadings_zn' => {
      begin_date  => '1961-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.13159,
      east_bound  =>  -70.89173,
      south_bound =>   42.56714,
      north_bound =>   42.67464
    },
    'swmm_ipswich21_huc12_polygon.hydrology' => {
      begin_date  => '1961-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.21963,
      east_bound  =>  -70.81332,
      south_bound =>   42.49842,
      north_bound =>   42.70367
    },
    'swmm_ipswich21_huc12_polygon.loadings_tn' => {
      begin_date  => '1961-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.21963,
      east_bound  =>  -70.81332,
      south_bound =>   42.49842,
      north_bound =>   42.70367
    },
    'swmm_ipswich21_huc12_polygon.loadings_tp' => {
      begin_date  => '1961-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.21963,
      east_bound  =>  -70.81332,
      south_bound =>   42.49842,
      north_bound =>   42.70367
    },
    'swmm_ipswich21_huc12_polygon.loadings_tss' => {
      begin_date  => '1961-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.21963,
      east_bound  =>  -70.81332,
      south_bound =>   42.49842,
      north_bound =>   42.70367
    },
    'swmm_ipswich21_huc12_polygon.loadings_zn' => {
      begin_date  => '1961-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.21963,
      east_bound  =>  -70.81332,
      south_bound =>   42.49842,
      north_bound =>   42.70367
    },
    'swmm_ipswich29_huc10_point.hydrology' => {
      begin_date  => '1961-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.13159,
      east_bound  =>  -71.13159,
      south_bound =>   42.56714,
      north_bound =>   42.56714
    },
    'swmm_ipswich29_huc10_point.loadings_tn' => {
      begin_date  => '1961-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.13159,
      east_bound  =>  -71.13159,
      south_bound =>   42.56714,
      north_bound =>   42.56714
    },
    'swmm_ipswich29_huc10_point.loadings_tp' => {
      begin_date  => '1961-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.13159,
      east_bound  =>  -71.13159,
      south_bound =>   42.56714,
      north_bound =>   42.56714
    },
    'swmm_ipswich29_huc10_point.loadings_tss' => {
      begin_date  => '1961-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.13159,
      east_bound  =>  -71.13159,
      south_bound =>   42.56714,
      north_bound =>   42.56714
    },
    'swmm_ipswich29_huc10_point.loadings_zn' => {
      begin_date  => '1961-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.13159,
      east_bound  =>  -71.13159,
      south_bound =>   42.56714,
      north_bound =>   42.56714
    },
    'swmm_ipswich29_huc10_polygon.hydrology' => {
      begin_date  => '1961-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.21963,
      east_bound  =>  -70.81332,
      south_bound =>   42.49842,
      north_bound =>   42.70367
    },
    'swmm_ipswich29_huc10_polygon.loadings_tn' => {
      begin_date  => '1961-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.21963,
      east_bound  =>  -70.81332,
      south_bound =>   42.49842,
      north_bound =>   42.70367
    },
    'swmm_ipswich29_huc10_polygon.loadings_tp' => {
      begin_date  => '1961-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.21963,
      east_bound  =>  -70.81332,
      south_bound =>   42.49842,
      north_bound =>   42.70367
    },
    'swmm_ipswich29_huc10_polygon.loadings_tss' => {
      begin_date  => '1961-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.21963,
      east_bound  =>  -70.81332,
      south_bound =>   42.49842,
      north_bound =>   42.70367
    },
    'swmm_ipswich29_huc10_polygon.loadings_zn' => {
      begin_date  => '1961-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.21963,
      east_bound  =>  -70.81332,
      south_bound =>   42.49842,
      north_bound =>   42.70367
    },
    'swmm_ipswich29_huc12_point.hydrology' => {
      begin_date  => '1961-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.13159,
      east_bound  =>  -70.89173,
      south_bound =>   42.56714,
      north_bound =>   42.67464
    },
    'swmm_ipswich29_huc12_point.loadings_tn' => {
      begin_date  => '1961-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.13159,
      east_bound  =>  -70.89173,
      south_bound =>   42.56714,
      north_bound =>   42.67464
    },
    'swmm_ipswich29_huc12_point.loadings_tp' => {
      begin_date  => '1961-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.13159,
      east_bound  =>  -70.89173,
      south_bound =>   42.56714,
      north_bound =>   42.67464
    },
    'swmm_ipswich29_huc12_point.loadings_tss' => {
      begin_date  => '1961-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.13159,
      east_bound  =>  -70.89173,
      south_bound =>   42.56714,
      north_bound =>   42.67464
    },
    'swmm_ipswich29_huc12_point.loadings_zn' => {
      begin_date  => '1961-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.13159,
      east_bound  =>  -70.89173,
      south_bound =>   42.56714,
      north_bound =>   42.67464
    },
    'swmm_ipswich29_huc12_polygon.hydrology' => {
      begin_date  => '1961-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.21963,
      east_bound  =>  -70.81332,
      south_bound =>   42.49842,
      north_bound =>   42.70367
    },
    'swmm_ipswich29_huc12_polygon.loadings_tn' => {
      begin_date  => '1961-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.21963,
      east_bound  =>  -70.81332,
      south_bound =>   42.49842,
      north_bound =>   42.70367
    },
    'swmm_ipswich29_huc12_polygon.loadings_tp' => {
      begin_date  => '1961-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.21963,
      east_bound  =>  -70.81332,
      south_bound =>   42.49842,
      north_bound =>   42.70367
    },
    'swmm_ipswich29_huc12_polygon.loadings_tss' => {
      begin_date  => '1961-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.21963,
      east_bound  =>  -70.81332,
      south_bound =>   42.49842,
      north_bound =>   42.70367
    },
    'swmm_ipswich29_huc12_polygon.loadings_zn' => {
      begin_date  => '1961-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.21963,
      east_bound  =>  -70.81332,
      south_bound =>   42.49842,
      north_bound =>   42.70367
    },
    'swmm_pawcatuck22_huc10_point.hydrology' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.84246,
      east_bound  =>  -71.63932,
      south_bound =>   41.44597,
      north_bound =>   41.56840
    },
    'swmm_pawcatuck22_huc10_point.loadings_tn' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.84246,
      east_bound  =>  -71.63932,
      south_bound =>   41.44597,
      north_bound =>   41.56840
    },
    'swmm_pawcatuck22_huc10_point.loadings_tp' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.84246,
      east_bound  =>  -71.63932,
      south_bound =>   41.44597,
      north_bound =>   41.56840
    },
    'swmm_pawcatuck22_huc10_point.loadings_tss' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.84246,
      east_bound  =>  -71.63932,
      south_bound =>   41.44597,
      north_bound =>   41.56840
    },
    'swmm_pawcatuck22_huc10_point.loadings_zn' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.84246,
      east_bound  =>  -71.63932,
      south_bound =>   41.44597,
      north_bound =>   41.56840
    },
    'swmm_pawcatuck22_huc10_polygon.hydrology' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.93323,
      east_bound  =>  -71.48661,
      south_bound =>   41.31577,
      north_bound =>   41.68141
    },
    'swmm_pawcatuck22_huc10_polygon.loadings_tn' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.93323,
      east_bound  =>  -71.48661,
      south_bound =>   41.31577,
      north_bound =>   41.68141
    },
    'swmm_pawcatuck22_huc10_polygon.loadings_tp' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.93323,
      east_bound  =>  -71.48661,
      south_bound =>   41.31577,
      north_bound =>   41.68141
    },
    'swmm_pawcatuck22_huc10_polygon.loadings_tss' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.93323,
      east_bound  =>  -71.48661,
      south_bound =>   41.31577,
      north_bound =>   41.68141
    },
    'swmm_pawcatuck22_huc10_polygon.loadings_zn' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.93323,
      east_bound  =>  -71.48661,
      south_bound =>   41.31577,
      north_bound =>   41.68141
    },
    'swmm_pawcatuck22_huc12_point.hydrology' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.88957,
      east_bound  =>  -71.54656,
      south_bound =>   41.38083,
      north_bound =>   41.59852
    },
    'swmm_pawcatuck22_huc12_point.loadings_tn' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.88957,
      east_bound  =>  -71.54656,
      south_bound =>   41.38083,
      north_bound =>   41.59852
    },
    'swmm_pawcatuck22_huc12_point.loadings_tp' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.88957,
      east_bound  =>  -71.54656,
      south_bound =>   41.38083,
      north_bound =>   41.59852
    },
    'swmm_pawcatuck22_huc12_point.loadings_tss' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.88957,
      east_bound  =>  -71.54656,
      south_bound =>   41.38083,
      north_bound =>   41.59852
    },
    'swmm_pawcatuck22_huc12_point.loadings_zn' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.88957,
      east_bound  =>  -71.54656,
      south_bound =>   41.38083,
      north_bound =>   41.59852
    },
    'swmm_pawcatuck22_huc12_polygon.hydrology' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.93323,
      east_bound  =>  -71.48661,
      south_bound =>   41.31577,
      north_bound =>   41.68141
    },
    'swmm_pawcatuck22_huc12_polygon.loadings_tn' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.93323,
      east_bound  =>  -71.48661,
      south_bound =>   41.31577,
      north_bound =>   41.68141
    },
    'swmm_pawcatuck22_huc12_polygon.loadings_tp' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.93323,
      east_bound  =>  -71.48661,
      south_bound =>   41.31577,
      north_bound =>   41.68141
    },
    'swmm_pawcatuck22_huc12_polygon.loadings_tss' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.93323,
      east_bound  =>  -71.48661,
      south_bound =>   41.31577,
      north_bound =>   41.68141
    },
    'swmm_pawcatuck22_huc12_polygon.loadings_zn' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.93323,
      east_bound  =>  -71.48661,
      south_bound =>   41.31577,
      north_bound =>   41.68141
    },
    'swmm_pawcatuck30_huc10_point.hydrology' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.84246,
      east_bound  =>  -71.63932,
      south_bound =>   41.44597,
      north_bound =>   41.56840
    },
    'swmm_pawcatuck30_huc10_point.loadings_tn' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.84246,
      east_bound  =>  -71.63932,
      south_bound =>   41.44597,
      north_bound =>   41.56840
    },
    'swmm_pawcatuck30_huc10_point.loadings_tp' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.84246,
      east_bound  =>  -71.63932,
      south_bound =>   41.44597,
      north_bound =>   41.56840
    },
    'swmm_pawcatuck30_huc10_point.loadings_tss' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.84246,
      east_bound  =>  -71.63932,
      south_bound =>   41.44597,
      north_bound =>   41.56840
    },
    'swmm_pawcatuck30_huc10_point.loadings_zn' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.84246,
      east_bound  =>  -71.63932,
      south_bound =>   41.44597,
      north_bound =>   41.56840
    },
    'swmm_pawcatuck30_huc10_polygon.hydrology' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.93323,
      east_bound  =>  -71.48661,
      south_bound =>   41.31577,
      north_bound =>   41.68141
    },
    'swmm_pawcatuck30_huc10_polygon.loadings_tn' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.93323,
      east_bound  =>  -71.48661,
      south_bound =>   41.31577,
      north_bound =>   41.68141
    },
    'swmm_pawcatuck30_huc10_polygon.loadings_tp' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.93323,
      east_bound  =>  -71.48661,
      south_bound =>   41.31577,
      north_bound =>   41.68141
    },
    'swmm_pawcatuck30_huc10_polygon.loadings_tss' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.93323,
      east_bound  =>  -71.48661,
      south_bound =>   41.31577,
      north_bound =>   41.68141
    },
    'swmm_pawcatuck30_huc10_polygon.loadings_zn' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.93323,
      east_bound  =>  -71.48661,
      south_bound =>   41.31577,
      north_bound =>   41.68141
    },
    'swmm_pawcatuck30_huc12_point.hydrology' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.88957,
      east_bound  =>  -71.54656,
      south_bound =>   41.38083,
      north_bound =>   41.59852
    },
    'swmm_pawcatuck30_huc12_point.loadings_tn' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.88957,
      east_bound  =>  -71.54656,
      south_bound =>   41.38083,
      north_bound =>   41.59852
    },
    'swmm_pawcatuck30_huc12_point.loadings_tp' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.88957,
      east_bound  =>  -71.54656,
      south_bound =>   41.38083,
      north_bound =>   41.59852
    },
    'swmm_pawcatuck30_huc12_point.loadings_tss' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.88957,
      east_bound  =>  -71.54656,
      south_bound =>   41.38083,
      north_bound =>   41.59852
    },
    'swmm_pawcatuck30_huc12_point.loadings_zn' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.88957,
      east_bound  =>  -71.54656,
      south_bound =>   41.38083,
      north_bound =>   41.59852
    },
    'swmm_pawcatuck30_huc12_polygon.hydrology' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.93323,
      east_bound  =>  -71.48661,
      south_bound =>   41.31577,
      north_bound =>   41.68141
    },
    'swmm_pawcatuck30_huc12_polygon.loadings_tn' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.93323,
      east_bound  =>  -71.48661,
      south_bound =>   41.31577,
      north_bound =>   41.68141
    },
    'swmm_pawcatuck30_huc12_polygon.loadings_tp' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.93323,
      east_bound  =>  -71.48661,
      south_bound =>   41.31577,
      north_bound =>   41.68141
    },
    'swmm_pawcatuck30_huc12_polygon.loadings_tss' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.93323,
      east_bound  =>  -71.48661,
      south_bound =>   41.31577,
      north_bound =>   41.68141
    },
    'swmm_pawcatuck30_huc12_polygon.loadings_zn' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.93323,
      east_bound  =>  -71.48661,
      south_bound =>   41.31577,
      north_bound =>   41.68141
    },
    'swmm_sudbury23_huc10_point.hydrology' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.48964,
      east_bound  =>  -71.46364,
      south_bound =>   42.31622,
      north_bound =>   42.44976
    },
    'swmm_sudbury23_huc10_point.loadings_tn' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.48964,
      east_bound  =>  -71.46364,
      south_bound =>   42.31622,
      north_bound =>   42.44976
    },
    'swmm_sudbury23_huc10_point.loadings_tp' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.48964,
      east_bound  =>  -71.46364,
      south_bound =>   42.31622,
      north_bound =>   42.44976
    },
    'swmm_sudbury23_huc10_point.loadings_tss' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.48964,
      east_bound  =>  -71.46364,
      south_bound =>   42.31622,
      north_bound =>   42.44976
    },
    'swmm_sudbury23_huc10_point.loadings_zn' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.48964,
      east_bound  =>  -71.46364,
      south_bound =>   42.31622,
      north_bound =>   42.44976
    },
    'swmm_sudbury23_huc10_polygon.hydrology' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.72177,
      east_bound  =>  -71.25150,
      south_bound =>   42.19754,
      north_bound =>   42.64674
    },
    'swmm_sudbury23_huc10_polygon.loadings_tn' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.72177,
      east_bound  =>  -71.25150,
      south_bound =>   42.19754,
      north_bound =>   42.64674
    },
    'swmm_sudbury23_huc10_polygon.loadings_tp' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.72177,
      east_bound  =>  -71.25150,
      south_bound =>   42.19754,
      north_bound =>   42.64674
    },
    'swmm_sudbury23_huc10_polygon.loadings_tss' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.72177,
      east_bound  =>  -71.25150,
      south_bound =>   42.19754,
      north_bound =>   42.64674
    },
    'swmm_sudbury23_huc10_polygon.loadings_zn' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.72177,
      east_bound  =>  -71.25150,
      south_bound =>   42.19754,
      north_bound =>   42.64674
    },
    'swmm_sudbury23_huc12_point.hydrology' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.65187,
      east_bound  =>  -71.30854,
      south_bound =>   42.25174,
      north_bound =>   42.57971
    },
    'swmm_sudbury23_huc12_point.loadings_tn' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.65187,
      east_bound  =>  -71.30854,
      south_bound =>   42.25174,
      north_bound =>   42.57971
    },
    'swmm_sudbury23_huc12_point.loadings_tp' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.65187,
      east_bound  =>  -71.30854,
      south_bound =>   42.25174,
      north_bound =>   42.57971
    },
    'swmm_sudbury23_huc12_point.loadings_tss' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.65187,
      east_bound  =>  -71.30854,
      south_bound =>   42.25174,
      north_bound =>   42.57971
    },
    'swmm_sudbury23_huc12_point.loadings_zn' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.65187,
      east_bound  =>  -71.30854,
      south_bound =>   42.25174,
      north_bound =>   42.57971
    },
    'swmm_sudbury23_huc12_polygon.hydrology' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.72177,
      east_bound  =>  -71.25150,
      south_bound =>   42.19754,
      north_bound =>   42.64674
    },
    'swmm_sudbury23_huc12_polygon.loadings_tn' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.72177,
      east_bound  =>  -71.25150,
      south_bound =>   42.19754,
      north_bound =>   42.64674
    },
    'swmm_sudbury23_huc12_polygon.loadings_tp' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.72177,
      east_bound  =>  -71.25150,
      south_bound =>   42.19754,
      north_bound =>   42.64674
    },
    'swmm_sudbury23_huc12_polygon.loadings_tss' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.72177,
      east_bound  =>  -71.25150,
      south_bound =>   42.19754,
      north_bound =>   42.64674
    },
    'swmm_sudbury23_huc12_polygon.loadings_zn' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.72177,
      east_bound  =>  -71.25150,
      south_bound =>   42.19754,
      north_bound =>   42.64674
    },
    'swmm_sudbury31_huc10_point.hydrology' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.48964,
      east_bound  =>  -71.46364,
      south_bound =>   42.31622,
      north_bound =>   42.44976
    },
    'swmm_sudbury31_huc10_point.loadings_tn' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.48964,
      east_bound  =>  -71.46364,
      south_bound =>   42.31622,
      north_bound =>   42.44976
    },
    'swmm_sudbury31_huc10_point.loadings_tp' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.48964,
      east_bound  =>  -71.46364,
      south_bound =>   42.31622,
      north_bound =>   42.44976
    },
    'swmm_sudbury31_huc10_point.loadings_tss' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.48964,
      east_bound  =>  -71.46364,
      south_bound =>   42.31622,
      north_bound =>   42.44976
    },
    'swmm_sudbury31_huc10_point.loadings_zn' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.48964,
      east_bound  =>  -71.46364,
      south_bound =>   42.31622,
      north_bound =>   42.44976
    },
    'swmm_sudbury31_huc10_polygon.hydrology' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.72177,
      east_bound  =>  -71.25150,
      south_bound =>   42.19754,
      north_bound =>   42.64674
    },
    'swmm_sudbury31_huc10_polygon.loadings_tn' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.72177,
      east_bound  =>  -71.25150,
      south_bound =>   42.19754,
      north_bound =>   42.64674
    },
    'swmm_sudbury31_huc10_polygon.loadings_tp' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.72177,
      east_bound  =>  -71.25150,
      south_bound =>   42.19754,
      north_bound =>   42.64674
    },
    'swmm_sudbury31_huc10_polygon.loadings_tss' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.72177,
      east_bound  =>  -71.25150,
      south_bound =>   42.19754,
      north_bound =>   42.64674
    },
    'swmm_sudbury31_huc10_polygon.loadings_zn' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.72177,
      east_bound  =>  -71.25150,
      south_bound =>   42.19754,
      north_bound =>   42.64674
    },
    'swmm_sudbury31_huc12_point.hydrology' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.65187,
      east_bound  =>  -71.30854,
      south_bound =>   42.25174,
      north_bound =>   42.57971
    },
    'swmm_sudbury31_huc12_point.loadings_tn' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.65187,
      east_bound  =>  -71.30854,
      south_bound =>   42.25174,
      north_bound =>   42.57971
    },
    'swmm_sudbury31_huc12_point.loadings_tp' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.65187,
      east_bound  =>  -71.30854,
      south_bound =>   42.25174,
      north_bound =>   42.57971
    },
    'swmm_sudbury31_huc12_point.loadings_tss' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.65187,
      east_bound  =>  -71.30854,
      south_bound =>   42.25174,
      north_bound =>   42.57971
    },
    'swmm_sudbury31_huc12_point.loadings_zn' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.65187,
      east_bound  =>  -71.30854,
      south_bound =>   42.25174,
      north_bound =>   42.57971
    },
    'swmm_sudbury31_huc12_polygon.hydrology' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.72177,
      east_bound  =>  -71.25150,
      south_bound =>   42.19754,
      north_bound =>   42.64674
    },
    'swmm_sudbury31_huc12_polygon.loadings_tn' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.72177,
      east_bound  =>  -71.25150,
      south_bound =>   42.19754,
      north_bound =>   42.64674
    },
    'swmm_sudbury31_huc12_polygon.loadings_tp' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.72177,
      east_bound  =>  -71.25150,
      south_bound =>   42.19754,
      north_bound =>   42.64674
    },
    'swmm_sudbury31_huc12_polygon.loadings_tss' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.72177,
      east_bound  =>  -71.25150,
      south_bound =>   42.19754,
      north_bound =>   42.64674
    },
    'swmm_sudbury31_huc12_polygon.loadings_zn' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.72177,
      east_bound  =>  -71.25150,
      south_bound =>   42.19754,
      north_bound =>   42.64674
    },
    'swmm_taunton24_huc10_point.hydrology' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.22375,
      east_bound  =>  -70.92933,
      south_bound =>   41.87677,
      north_bound =>   42.05023
    },
    'swmm_taunton24_huc10_point.loadings_tn' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.22375,
      east_bound  =>  -70.92933,
      south_bound =>   41.87677,
      north_bound =>   42.05023
    },
    'swmm_taunton24_huc10_point.loadings_tp' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.22375,
      east_bound  =>  -70.92933,
      south_bound =>   41.87677,
      north_bound =>   42.05023
    },
    'swmm_taunton24_huc10_point.loadings_tss' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.22375,
      east_bound  =>  -70.92933,
      south_bound =>   41.87677,
      north_bound =>   42.05023
    },
    'swmm_taunton24_huc10_point.loadings_zn' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.22375,
      east_bound  =>  -70.92933,
      south_bound =>   41.87677,
      north_bound =>   42.05023
    },
    'swmm_taunton24_huc10_polygon.hydrology' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.34374,
      east_bound  =>  -70.71200,
      south_bound =>   41.70641,
      north_bound =>   42.16041
    },
    'swmm_taunton24_huc10_polygon.loadings_tn' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.34374,
      east_bound  =>  -70.71200,
      south_bound =>   41.70641,
      north_bound =>   42.16041
    },
    'swmm_taunton24_huc10_polygon.loadings_tp' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.34374,
      east_bound  =>  -70.71200,
      south_bound =>   41.70641,
      north_bound =>   42.16041
    },
    'swmm_taunton24_huc10_polygon.loadings_tss' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.34374,
      east_bound  =>  -70.71200,
      south_bound =>   41.70641,
      north_bound =>   42.16041
    },
    'swmm_taunton24_huc10_polygon.loadings_zn' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.34374,
      east_bound  =>  -70.71200,
      south_bound =>   41.70641,
      north_bound =>   42.16041
    },
    'swmm_taunton24_huc12_point.hydrology' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.26141,
      east_bound  =>  -70.84823,
      south_bound =>   41.79535,
      north_bound =>   42.08117
    },
    'swmm_taunton24_huc12_point.loadings_tn' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.26141,
      east_bound  =>  -70.84823,
      south_bound =>   41.79535,
      north_bound =>   42.08117
    },
    'swmm_taunton24_huc12_point.loadings_tp' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.26141,
      east_bound  =>  -70.84823,
      south_bound =>   41.79535,
      north_bound =>   42.08117
    },
    'swmm_taunton24_huc12_point.loadings_tss' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.26141,
      east_bound  =>  -70.84823,
      south_bound =>   41.79535,
      north_bound =>   42.08117
    },
    'swmm_taunton24_huc12_point.loadings_zn' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.26141,
      east_bound  =>  -70.84823,
      south_bound =>   41.79535,
      north_bound =>   42.08117
    },
    'swmm_taunton24_huc12_polygon.hydrology' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.34374,
      east_bound  =>  -70.71200,
      south_bound =>   41.70641,
      north_bound =>   42.16041
    },
    'swmm_taunton24_huc12_polygon.loadings_tn' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.34374,
      east_bound  =>  -70.71200,
      south_bound =>   41.70641,
      north_bound =>   42.16041
    },
    'swmm_taunton24_huc12_polygon.loadings_tp' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.34374,
      east_bound  =>  -70.71200,
      south_bound =>   41.70641,
      north_bound =>   42.16041
    },
    'swmm_taunton24_huc12_polygon.loadings_tss' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.34374,
      east_bound  =>  -70.71200,
      south_bound =>   41.70641,
      north_bound =>   42.16041
    },
    'swmm_taunton24_huc12_polygon.loadings_zn' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.34374,
      east_bound  =>  -70.71200,
      south_bound =>   41.70641,
      north_bound =>   42.16041
    },
    'swmm_taunton32_huc10_point.hydrology' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.22375,
      east_bound  =>  -70.92933,
      south_bound =>   41.87677,
      north_bound =>   42.05023
    },
    'swmm_taunton32_huc10_point.loadings_tn' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.22375,
      east_bound  =>  -70.92933,
      south_bound =>   41.87677,
      north_bound =>   42.05023
    },
    'swmm_taunton32_huc10_point.loadings_tp' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.22375,
      east_bound  =>  -70.92933,
      south_bound =>   41.87677,
      north_bound =>   42.05023
    },
    'swmm_taunton32_huc10_point.loadings_tss' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.22375,
      east_bound  =>  -70.92933,
      south_bound =>   41.87677,
      north_bound =>   42.05023
    },
    'swmm_taunton32_huc10_point.loadings_zn' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.22375,
      east_bound  =>  -70.92933,
      south_bound =>   41.87677,
      north_bound =>   42.05023
    },
    'swmm_taunton32_huc10_polygon.hydrology' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.34374,
      east_bound  =>  -70.71200,
      south_bound =>   41.70641,
      north_bound =>   42.16041
    },
    'swmm_taunton32_huc10_polygon.loadings_tn' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.34374,
      east_bound  =>  -70.71200,
      south_bound =>   41.70641,
      north_bound =>   42.16041
    },
    'swmm_taunton32_huc10_polygon.loadings_tp' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.34374,
      east_bound  =>  -70.71200,
      south_bound =>   41.70641,
      north_bound =>   42.16041
    },
    'swmm_taunton32_huc10_polygon.loadings_tss' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.34374,
      east_bound  =>  -70.71200,
      south_bound =>   41.70641,
      north_bound =>   42.16041
    },
    'swmm_taunton32_huc10_polygon.loadings_zn' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.34374,
      east_bound  =>  -70.71200,
      south_bound =>   41.70641,
      north_bound =>   42.16041
    },
    'swmm_taunton32_huc12_point.hydrology' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.26141,
      east_bound  =>  -70.84823,
      south_bound =>   41.79535,
      north_bound =>   42.08117
    },
    'swmm_taunton32_huc12_point.loadings_tn' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.26141,
      east_bound  =>  -70.84823,
      south_bound =>   41.79535,
      north_bound =>   42.08117
    },
    'swmm_taunton32_huc12_point.loadings_tp' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.26141,
      east_bound  =>  -70.84823,
      south_bound =>   41.79535,
      north_bound =>   42.08117
    },
    'swmm_taunton32_huc12_point.loadings_tss' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.26141,
      east_bound  =>  -70.84823,
      south_bound =>   41.79535,
      north_bound =>   42.08117
    },
    'swmm_taunton32_huc12_point.loadings_zn' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.26141,
      east_bound  =>  -70.84823,
      south_bound =>   41.79535,
      north_bound =>   42.08117
    },
    'swmm_taunton32_huc12_polygon.hydrology' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.34374,
      east_bound  =>  -70.71200,
      south_bound =>   41.70641,
      north_bound =>   42.16041
    },
    'swmm_taunton32_huc12_polygon.loadings_tn' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.34374,
      east_bound  =>  -70.71200,
      south_bound =>   41.70641,
      north_bound =>   42.16041
    },
    'swmm_taunton32_huc12_polygon.loadings_tp' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.34374,
      east_bound  =>  -70.71200,
      south_bound =>   41.70641,
      north_bound =>   42.16041
    },
    'swmm_taunton32_huc12_polygon.loadings_tss' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.34374,
      east_bound  =>  -70.71200,
      south_bound =>   41.70641,
      north_bound =>   42.16041
    },
    'swmm_taunton32_huc12_polygon.loadings_zn' => {
      begin_date  => '1960-01-01',
      end_date    => '2015-12-31',
      west_bound  =>  -71.34374,
      east_bound  =>  -70.71200,
      south_bound =>   41.70641,
      north_bound =>   42.16041
    },
    'greenspace_housing' => {
      begin_date  => '2014-01-01',
      end_date    => '2014-12-31',
      west_bound  =>  -124.902,
      east_bound  =>   -66.019,
      south_bound =>   24.395,
      north_bound =>   49.796
    },
    'seagrass_point' => {
      begin_date  => '1978-01-01',
      end_date    => '2017-12-31',
      west_bound  =>  -124.71,
      east_bound  =>  -117.12,
      south_bound =>   32.55,
      north_bound =>   48.98
    },
    'seagrass_point2' => {
      begin_date  => '1978-01-01',
      end_date    => '2017-12-31',
      west_bound  =>  -124.71,
      east_bound  =>  -117.12,
      south_bound =>   32.55,
      north_bound =>   48.98
    },
    'seagrass_polygon' => {
      begin_date  => '1978-01-01',
      end_date    => '2017-12-31',
      west_bound  =>  -124.51,
      east_bound  =>  -117.98,
      south_bound =>   32.59,
      north_bound =>   48.99
    },
    'esat1_water_quality' => {
      begin_date  => '1984-01-01',
      end_date    => '2021-12-31',
      west_bound  =>  -124.02,
      east_bound  =>  -70.52,
      south_bound =>   24.91,
      north_bound =>   48.56
    },
    'esat2_water_quality' => {
      begin_date  => '2018-01-01',
      end_date    => '2020-12-31',
      west_bound  =>  -124.02,
      east_bound  =>  -70.52,
      south_bound =>   24.91,
      north_bound =>   48.56
    }
  },

  'gulf' => {
    'states' => {
      begin_date  => '2018-04-01',
      end_date    => '2018-04-01',
      west_bound  => -106.64565,
      east_bound  =>  -79.97431,
      south_bound =>   24.39631,
      north_bound =>   36.50070
    },
    'counties' => {
      begin_date  => '2018-04-01',
      end_date    => '2018-04-01',
      west_bound  => -106.64565,
      east_bound  =>  -79.97431,
      south_bound =>   24.39631,
      north_bound =>   36.50070
    },
    'cities' => {
      begin_date  => '2018-04-01',
      end_date    => '2018-04-01',
      west_bound  => -106.17019,
      east_bound  =>  -80.03116,
      south_bound =>   24.54409,
      north_bound =>   36.40657
    },
    'roads' => {
      begin_date  => '2018-04-01',
      end_date    => '2018-04-01',
      west_bound  => -106.58294,
      east_bound  =>  -80.06662,
      south_bound =>   25.45512,
      north_bound =>   35.32735
    },
    'tributaries' => {
      begin_date  => '2018-04-01',
      end_date    => '2018-04-01',
      west_bound  => -101.42156,
      east_bound  =>  -80.05212,
      south_bound =>   25.10267,
      north_bound =>   33.60087
    },
    'hucs' => {
      begin_date  => '2018-04-01',
      end_date    => '2018-04-01',
      west_bound  => -103.87052,
      east_bound  =>  -79.97402,
      south_bound =>   24.39533,
      north_bound =>   35.15195
    },
    'watersheds' => {
      begin_date  => '2018-04-01',
      end_date    => '2018-04-01',
      west_bound  => -103.87144,
      east_bound  =>  -80.03321,
      south_bound =>   24.95149,
      north_bound =>   35.15012
    },
    'estuaries' => {
      begin_date  => '2018-04-01',
      end_date    => '2018-04-01',
      west_bound  =>  -97.77243,
      east_bound  =>  -80.03195,
      south_bound =>   24.55000,
      north_bound =>   31.00019
    },
    'lakes' => {
      begin_date  => '2018-04-01',
      end_date    => '2018-04-01',
      west_bound  => -103.99124,
      east_bound  =>  -80.61124,
      south_bound =>   25.13734,
      north_bound =>   35.72883
    },
    'legislative_districts' => {
      begin_date  => '2018-04-01',
      end_date    => '2018-04-01',
      west_bound  => -106.64565,
      east_bound  =>  -79.97431,
      south_bound =>   24.39631,
      north_bound =>   36.50070
    },

    'nlcd_2011' => {
      begin_date  => '2011-01-01',
      end_date    => '2011-12-31',
      west_bound  => -105,
      east_bound  => -77,
      south_bound => 24,
      north_bound => 38
    },
    'nlcd_2006' => {
      begin_date  => '2006-01-01',
      end_date    => '2006-12-31',
      west_bound  => -105,
      east_bound  => -77,
      south_bound => 24,
      north_bound => 38
    },
    'nlcd_2001' => {
      begin_date  => '2001-01-01',
      end_date    => '2001-12-31',
      west_bound  => -105,
      east_bound  => -77,
      south_bound => 24,
      north_bound => 38
    },
    'nlcd_1992' => {
      begin_date  => '1992-01-01',
      end_date    => '1992-12-31',
      west_bound  => -105,
      east_bound  => -77,
      south_bound => 24,
      north_bound => 38
    },
    'nlcd_change_1992_2001' => {
      begin_date  => '1992-01-01',
      end_date    => '2001-12-31',
      west_bound  => -105,
      east_bound  => -77,
      south_bound => 24,
      north_bound => 38
    },
    'nlcd_change_2001_2006' => {
      begin_date  => '2001-01-01',
      end_date    => '2006-12-31',
      west_bound  => -105,
      east_bound  => -77,
      south_bound => 24,
      north_bound => 38
    },
    'ccap_2006' => {
      begin_date  => '2006-01-01',
      end_date    => '2006-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'ccap_2001' => {
      begin_date  => '2001-01-01',
      end_date    => '2001-12-31',
      west_bound  => -102,
      east_bound  => -77,
      south_bound => 24,
      north_bound => 36
    },
    'decadal_iclus_land_use_2000_2100' => {
      begin_date  => '2000-01-01',
      end_date    => '2100-12-31',
      west_bound  => -105,
      east_bound  => -77,
      south_bound => 24,
      north_bound => 38
    },
    'iclus_land_use_watershed' => {
      begin_date  => '2000-01-01',
      end_date    => '2100-12-31',
      west_bound  => -105,
      east_bound  => -77,
      south_bound => 24,
      north_bound => 38
    },
    'nlcd_2011_land_use_watershed' => {
      begin_date  => '2011-01-01',
      end_date    => '2011-12-31',
      west_bound  => -105,
      east_bound  => -77,
      south_bound => 24,
      north_bound => 38
    },
    'nlcd_2006_land_use_watershed' => {
      begin_date  => '2006-01-01',
      end_date    => '2006-12-31',
      west_bound  => -105,
      east_bound  => -77,
      south_bound => 24,
      north_bound => 38
    },
    'impervious_2006' => {
      begin_date  => '2006-01-01',
      end_date    => '2006-12-31',
      west_bound  => -105,
      east_bound  => -77,
      south_bound => 24,
      north_bound => 38
    },
    'impervious_2001' => {
      begin_date  => '2001-01-01',
      end_date    => '2001-12-31',
      west_bound  => -105,
      east_bound  => -77,
      south_bound => 24,
      north_bound => 38
    },
    'impervious_change_2001_2006' => {
      begin_date  => '2001-01-01',
      end_date    => '2006-12-31',
      west_bound  => -105,
      east_bound  => -77,
      south_bound => 24,
      north_bound => 38
    },
    'decadal_iclus_a1_impervious_2010_2100' => {
      begin_date  => '2010-01-01',
      end_date    => '2100-12-31',
      west_bound  => -105,
      east_bound  => -77,
      south_bound => 24,
      north_bound => 38
    },
    'decadal_iclus_a2_impervious_2010_2100' => {
      begin_date  => '2010-01-01',
      end_date    => '2100-12-31',
      west_bound  => -105,
      east_bound  => -77,
      south_bound => 24,
      north_bound => 38
    },
    'decadal_iclus_b1_impervious_2010_2100' => {
      begin_date  => '2010-01-01',
      end_date    => '2100-12-31',
      west_bound  => -105,
      east_bound  => -77,
      south_bound => 24,
      north_bound => 38
    },
    'decadal_iclus_b2_impervious_2010_2100' => {
      begin_date  => '2010-01-01',
      end_date    => '2100-12-31',
      west_bound  => -105,
      east_bound  => -77,
      south_bound => 24,
      north_bound => 38
    },
    'decadal_iclus_bc_impervious_2010_2100' => {
      begin_date  => '2010-01-01',
      end_date    => '2100-12-31',
      west_bound  => -105,
      east_bound  => -77,
      south_bound => 24,
      north_bound => 38
    },

    'elevation' => {
      begin_date  => '2004-01-01',
      end_date    => '2004-12-31',
      west_bound  => -180,
      east_bound  => -66,
      south_bound => 18,
      north_bound => 66
    },

    # From: https://www.ngdc.noaa.gov/dem/squareCellGrid/search

    'elevation_alligator_cove_fl_navd88' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -81.246,
      east_bound  => -81.225,
      south_bound => 25.598,
      north_bound => 25.612
    },
    'elevation_apalachicola_bay_fl_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -85.31,
      east_bound  => -84.378,
      south_bound => 29.51,
      north_bound => 30.034
    },
    'elevation_apalachee_bay_fl_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -84.664,
      east_bound  => -83.704,
      south_bound => 29.802,
      north_bound => 30.385
    },
    'elevation_aransas_bay_1_tx_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -97.198,
      east_bound  => -96.853,
      south_bound => 27.817,
      north_bound => 28.274
    },
    'elevation_aransas_bay_2_tx_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -97.238,
      east_bound  => -97.004,
      south_bound => 27.983,
      north_bound => 28.210
    },
    'elevation_baffin_bay_tx_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -97.751,
      east_bound  => -97.395,
      south_bound => 27.216,
      north_bound => 27.482
    },
    'elevation_big_spring_creek_fl_mllw' => {
      begin_date  => '2020-06-19',
      end_date    => '2020-06-19',
      west_bound  => -83.754,
      east_bound  => -83.728,
      south_bound => 29.957,
      north_bound => 29.976
    },
    'elevation_biloxi_ms_mhw' => {
      begin_date  => '2007-03-29',
      end_date    => '2007-03-29',
      west_bound  => -89.3,
      east_bound  => -88.3,
      south_bound => 29.7,
      north_bound => 30.6
    },
    'elevation_biscayne_bay_fl_navd88' => {
      begin_date  => '2020-06-19',
      end_date    => '2020-06-19',
      west_bound  => -80.542,
      east_bound  => -79.981,
      south_bound => 25.225,
      north_bound => 26.074
    },
    'elevation_blackwater_sound_fl_navd88' => {
      begin_date  => '2020-06-19',
      end_date    => '2020-06-19',
      west_bound  => -80.425,
      east_bound  => -80.356,
      south_bound => 25.089,
      north_bound => 25.174
    },
    'elevation_boggy_bay_fl_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -82.748,
      east_bound  => -82.71,
      south_bound => 28.271,
      north_bound => 28.306
    },
    'elevation_boggy_creek_fl_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -83.337,
      east_bound  => -82.267,
      south_bound => 29.446,
      north_bound => 29.514
    },
    'elevation_brazos_river_tx_navd88' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -95.589,
      east_bound  => -95.279,
      south_bound => 28.828,
      north_bound => 29.039
    },
    'elevation_broad_river_rodgers_river_fl_navd88' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -81.206,
      east_bound  => -80.996,
      south_bound => 25.423,
      north_bound => 25.547
    },
    'elevation_buck_creek_fl_mllw' => {
      begin_date  => '2020-06-19',
      end_date    => '2020-06-19',
      west_bound  => -83.402,
      east_bound  => -83.388,
      south_bound => 29.542,
      north_bound => 29.555
    },
    'elevation_bumblebee_creek_fl_mllw' => {
      begin_date  => '2020-06-19',
      end_date    => '2020-06-19',
      west_bound  => -83.197,
      east_bound  => -83.147,
      south_bound => 29.342,
      north_bound => 29.388
    },
    'elevation_central_fl_mhw' => {
      begin_date  => '2014-11-03',
      end_date    => '2014-11-03',
      west_bound  => -80.88,
      east_bound  => -79.95,
      south_bound => 27.27,
      north_bound => 28.88
    },
    'elevation_central_fl_navd88' => {
      begin_date  => '2014-11-03',
      end_date    => '2014-11-03',
      west_bound  => -80.88,
      east_bound  => -79.95,
      south_bound => 27.27,
      north_bound => 28.88
    },
    'elevation_central_gulf_msl' => {
      begin_date  => '2001-07-20',
      end_date    => '2001-07-20',
      west_bound  => -94,
      east_bound  => -87,
      south_bound => 24,
      north_bound => 35
    },
    'elevation_charlotte_harbor_fl_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -82.5112,
      east_bound  => -81.643,
      south_bound => 26.345,
      north_bound => 27.157
    },
    'elevation_choctawhatchee_bay_fl_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -86.663,
      east_bound  => -85.831,
      south_bound => 30.242,
      north_bound => 30.591
    },
    'elevation_clay_creek_fl_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -83.521,
      east_bound  => -83.471,
      south_bound => 29.711,
      north_bound => 29.763
    },
    'elevation_clearwater_saint_joseph_sound_fl_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -82.907,
      east_bound  => -82.675,
      south_bound => 27.8631,
      north_bound => 28.258
    },
    'elevation_cocohatchee_river_fl_navd88' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -81.852,
      east_bound  => -81.794,
      south_bound => 26.258,
      north_bound => 26.332
    },
    'elevation_corpus_christi_bay_tx_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -97.604,
      east_bound  => -97.042,
      south_bound => 27.599,
      north_bound => 27.912
    },
    'elevation_corpus_christi_tx_mhw' => {
      begin_date  => '2007-05-04',
      end_date    => '2007-05-04',
      west_bound  => -97.7,
      east_bound  => -96.65,
      south_bound => 27.3,
      north_bound => 28.15
    },
    'elevation_crooked_creek_fl_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -83.56,
      east_bound  => -83.513,
      south_bound => 29.738,
      north_bound => 29.761
    },
    'elevation_crystal_river_1_fl_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -82.799,
      east_bound  => -82.546,
      south_bound => 28.807,
      north_bound => 29.007
    },
    'elevation_crystal_river_2_fl_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -82.77,
      east_bound  => -82.617,
      south_bound => 28.707,
      north_bound => 28.851
    },
    'elevation_crystal_river_3_fl_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -82.73,
      east_bound  => -82.552,
      south_bound => 28.611,
      north_bound => 28.766
    },
    'elevation_cumberland_sound_fl_navd88' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -81.666,
      east_bound  => -81.378,
      south_bound => 30.609,
      north_bound => 31.0
    },
    'elevation_dallus_creek_fl_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -83.504,
      east_bound  => -83.456,
      south_bound => 29.707,
      north_bound => 29.729
    },
    'elevation_days_creek_fl_mllw' => {
      begin_date  => '2007-10-31',
      end_date    => '2007-10-31',
      west_bound  => -83.41,
      east_bound  => -83.386,
      south_bound => 29.562,
      north_bound => 29.577
    },
    'elevation_daytona_beach_fl_mhw' => {
      begin_date  => '2007-10-31',
      end_date    => '2007-10-31',
      west_bound  => -81.4,
      east_bound  => -80.5,
      south_bound => 28.85,
      north_bound => 29.7
    },
    'elevation_deadman_bay_fl_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -83.475,
      east_bound  => -83.332,
      south_bound => 29.576,
      north_bound => 29.72
    },
    'elevation_eaglenest_creek_fl_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -83.734,
      east_bound  => -83.72,
      south_bound => 29.94,
      north_bound => 29.954
    },
    'elevation_eastern_gulf_msl' => {
      begin_date  => '2001-07-20',
      end_date    => '2001-07-20',
      west_bound  => -87,
      east_bound  => -78,
      south_bound => 24,
      north_bound => 35
    },
    'elevation_eleven_prong_bay_fl_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -82.823,
      east_bound  => -82.774,
      south_bound => 29.071,
      north_bound => 29.108
    },
    'elevation_estero_bay_fl_navd88' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -81.947,
      east_bound  => -81.779,
      south_bound => 26.314,
      north_bound => 26.524
    },
    'elevation_fillman_bayou_fl_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -82.703,
      east_bound  => -82.674,
      south_bound => 28.401,
      north_bound => 28.421
    },
    'elevation_first_bay_fl_navd88' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -81.244,
      east_bound  => -80.98,
      south_bound => 25.474,
      north_bound => 25.666
    },
    'elevation_fish_creek_fl_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -83.586,
      east_bound  => -83.496,
      south_bound => 29.758,
      north_bound => 29.812
    },
    'elevation_florida_bay_fl_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -81.807,
      east_bound  => -80.189,
      south_bound => 24.399,
      north_bound => 25.503
    },
    'elevation_galvaston_bay_1_tx_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -95.253,
      east_bound  => -94.654,
      south_bound => 29.487,
      north_bound => 29.806
    },
    'elevation_galvaston_bay_2_tx_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -95.257,
      east_bound  => -94.446,
      south_bound => 28.981,
      north_bound => 29.593
    },
    'elevation_galveston_tx_mhw' => {
      begin_date  => '2007-05-14',
      end_date    => '2007-05-14',
      west_bound  => -95.25,
      east_bound  => -94.3,
      south_bound => 28.85,
      north_bound => 29.8
    },
    'elevation_horseshoe_cove_fl_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -83.314,
      east_bound  => -83.164,
      south_bound => 29.371,
      north_bound => 29.517
    },
    'elevation_icw4_fl_navd88' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -82.456,
      east_bound  => -82.395,
      south_bound => 27.011,
      north_bound => 27.104
    },
    'elevation_icw_wolf_mobile_al_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -87.735,
      east_bound  => -87.608,
      south_bound => 30.265,
      north_bound => 30.315
    },
    'elevation_indian_bay_fl_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -81.158,
      east_bound  => -79.836,
      south_bound => 27.009,
      north_bound => 28.912
    },
    'elevation_indian_bay_fl_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -82.688,
      east_bound  => -82.657,
      south_bound => 28.433,
      north_bound => 28.467
    },
    'elevation_indian_river_bay_fl_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -81.158,
      east_bound  => -79.836,
      south_bound => 27.009,
      north_bound => 28.912
    },
    'elevation_johnson_creek_fl_navd88' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -83.213,
      east_bound  => -83.175,
      south_bound => 29.378,
      north_bound => 29.403
    },
    'elevation_key_west_fl_mhw' => {
      begin_date  => '2011-09-01',
      end_date    => '2011-09-01',
      west_bound  => -82.18,
      east_bound  => -81.27,
      south_bound => 23.98,
      north_bound => 25.02
    },
    'elevation_key_west_fl_navd88' => {
      begin_date  => '2011-09-01',
      end_date    => '2011-09-01',
      west_bound  => -82.18,
      east_bound  => -81.27,
      south_bound => 23.98,
      north_bound => 25.02
    },
    'elevation_lake_worth_fl_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -80.223,
      east_bound  => -80.009,
      south_bound => 25.925,
      north_bound => 26.331
    },
    'elevation_little_lagoon_al_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -87.834,
      east_bound  => -87.686,
      south_bound => 30.221,
      north_bound => 30.272
    },
    'elevation_little_pine_1_bay_fl_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -82.691,
      east_bound  => -82.64,
      south_bound => 28.459,
      north_bound => 28.509
    },
    'elevation_little_rocky_river_fl_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -83.358,
      east_bound  => -83.317,
      south_bound => 29.471,
      north_bound => 29.504
    },
    'elevation_little_spring_creek_fl_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -83.663,
      east_bound  => -83.639,
      south_bound => 29.897,
      north_bound => 29.916
    },
    'elevation_lows_bay_fl_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -82.823,
      east_bound  => -82.76,
      south_bound => 29.05,
      north_bound => 29.095
    },
    'elevation_loxahatchee_river_fl_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -80.319,
      east_bound  => -79.883,
      south_bound => 26.321,
      north_bound => 27.18
    },
    'elevation_matagorda_bay_tx_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -96.669,
      east_bound  => -95.512,
      south_bound => 28.290,
      north_bound => 28.898
    },
    'elevation_miami_fl_mhw' => {
      begin_date  => '2015-01-21',
      end_date    => '2015-01-21',
      west_bound  => -80.41,
      east_bound  => -79.35,
      south_bound => 25.28,
      north_bound => 26.32
    },
    'elevation_miami_fl_navd88' => {
      begin_date  => '2015-01-21',
      end_date    => '2015-01-21',
      west_bound  => -80.41,
      east_bound  => -79.35,
      south_bound => 25.28,
      north_bound => 26.32
    },
    'elevation_mississippi_sound_al_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -89.524,
      east_bound  => -88.043,
      south_bound => 29.917,
      north_bound => 30.558
    },
    'elevation_mississippi_sound_la_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -89.524,
      east_bound  => -88.043,
      south_bound => 29.917,
      north_bound => 30.558
    },
    'elevation_mississippi_sound_ms_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -89.524,
      east_bound  => -88.043,
      south_bound => 29.917,
      north_bound => 30.558
    },
    'elevation_mobile_al_mhw' => {
      begin_date  => '2009-11-30',
      end_date    => '2009-11-30',
      west_bound  => -88.3,
      east_bound  => -87.65,
      south_bound => 30,
      north_bound => 31
    },
    'elevation_mobile_bay_al_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -88.275,
      east_bound  => -87.594,
      south_bound => 30.157,
      north_bound => 31.038
    },
    'elevation_mobile_al_navd88' => {
      begin_date  => '2009-11-30',
      end_date    => '2009-11-30',
      west_bound  => -88.3,
      east_bound  => -87.65,
      south_bound => 30,
      north_bound => 31
    },
    'elevation_naples_bay_fl_navd88' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -81.813,
      east_bound  => -81.739,
      south_bound => 26.05,
      north_bound => 26.173
    },
    'elevation_nassau_sound_fl_navd88' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -81.732,
      east_bound  => -81.367,
      south_bound => 30.406,
      north_bound => 30.660
    },
    'elevation_new_orleans_la_mhw' => {
      begin_date  => '2010-04-30',
      end_date    => '2010-04-30',
      west_bound  => -90.65,
      east_bound  => -89.3,
      south_bound => 29.7,
      north_bound => 30.5
    },
    'elevation_new_orleans_la_mllw' => {
      begin_date  => '2010-04-30',
      end_date    => '2010-04-30',
      west_bound  => -90.65,
      east_bound  => -89.3,
      south_bound => 29.7,
      north_bound => 30.5
    },
    'elevation_new_orleans_la_navd88' => {
      begin_date  => '2010-04-30',
      end_date    => '2010-04-30',
      west_bound  => -90.65,
      east_bound  => -89.3,
      south_bound => 29.7,
      north_bound => 30.5
    },
    'elevation_northern_gulf_navd88' => {
      begin_date  => '2010-12-31',
      end_date    => '2010-12-31',
      west_bound  => -90.75,
      east_bound  => -85,
      south_bound => 28.5,
      north_bound => 31.25
    },
    'elevation_outer_doctors_bay_fl_navd88' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -81.821,
      east_bound  => -81.800,
      south_bound => 26.163,
      north_bound => 26.208
    },
    'elevation_oyster_creek_fl_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -83.577,
      east_bound  => -83.564,
      south_bound => 29.765,
      north_bound => 29.776
    },
    'elevation_oyster_creek_tx_navd88' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -95.353,
      east_bound  => -95.250,
      south_bound => 28.952,
      north_bound => 29.028
    },
    'elevation_palm_beach_fl_mhw' => {
      begin_date  => '2010-09-16',
      end_date    => '2010-09-16',
      west_bound  => -80.36,
      east_bound  => -79.42,
      south_bound => 26.29,
      north_bound => 27.31
    },
    'elevation_palm_beach_fl_navd88' => {
      begin_date  => '2010-09-16',
      end_date    => '2010-09-16',
      west_bound  => -80.36,
      east_bound  => -79.42,
      south_bound => 26.29,
      north_bound => 27.31
    },
    'elevation_panama_city_fl_mhw' => {
      begin_date  => '2010-07-01',
      end_date    => '2010-07-01',
      west_bound  => -86.1,
      east_bound  => -85.2,
      south_bound => 29.55,
      north_bound => 30.5
    },
    'elevation_panama_city_fl_navd88' => {
      begin_date  => '2010-07-01',
      end_date    => '2010-07-01',
      west_bound  => -86.1,
      east_bound  => -85.2,
      south_bound => 29.55,
      north_bound => 30.5
    },
    'elevation_perdido_bay_fl_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -87.655,
      east_bound  => -87.308,
      south_bound => 30.247,
      north_bound => 30.538
    },
    'elevation_south_padre_island_tx_mhw' => {
      begin_date  => '2012-05-04',
      end_date    => '2012-05-04',
      west_bound  => -97.65,
      east_bound  => -96.8,
      south_bound => 25.82,
      north_bound => 26.86
    },
    'elevation_pithlachascotee_river_fl_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -82.752,
      east_bound  => -82.704,
      south_bound => 28.231,
      north_bound => 28.289
    },
    'elevation_ponce_de_leon_inlet_fl_navd88' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -81.290,
      east_bound  => -80.488,
      south_bound => 28.641,
      north_bound => 29.506
    },
    'elevation_rock_1_bay_fl_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -82.669,
      east_bound  => -82.635,
      south_bound => 28.538,
      north_bound => 28.585
    },
    'elevation_sabine_lake_la_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -94.164,
      east_bound  => -93.632,
      south_bound => 29.671,
      north_bound => 30.212
    },
    'elevation_sabine_lake_tx_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -94.164,
      east_bound  => -93.632,
      south_bound => 29.671,
      north_bound => 30.212
    },
    'elevation_saint_andrew_bay_fl_navd88' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -85.971,
      east_bound  => -85.016,
      south_bound => 29.745,
      north_bound => 30.390
    },
    'elevation_saint_augustine_fl_navd88' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -81.543,
      east_bound  => -81.039,
      south_bound => 29.544,
      north_bound => 30.271
    },
    'elevation_saint_johns_river_fl_navd88' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -81.879,
      east_bound  => -81.339,
      south_bound => 29.946,
      north_bound => 30.522
    },
    'elevation_saint_joseph_bay_fl_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -85.442,
      east_bound  => -85.217,
      south_bound => 29.661,
      north_bound => 29.931
    },
    'elevation_salt_creek_fl_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -83.595,
      east_bound  => -83.576,
      south_bound => 29.817,
      north_bound => 29.836
    },
    'elevation_salt_creek_2_fl_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -83.480,
      east_bound  => -83.416,
      south_bound => 29.678,
      north_bound => 29.735
    },
    'elevation_san_antonio_bay_tx_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -96.914,
      east_bound  => -96.396,
      south_bound => 28.097,
      north_bound => 28.490
    },
    'elevation_sanders_creek_fl_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -83.208,
      east_bound  => -83.167,
      south_bound => 29.359,
      north_bound => 29.402
    },
    'elevation_sarasota_bay_fl_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -82.753,
      east_bound  => -82.376,
      south_bound => 27.069,
      north_bound => 27.504
    },
    'elevation_sink_creek_fl_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -83.405,
      east_bound  => -83.389,
      south_bound => 29.553,
      north_bound => 29.567
    },
    'elevation_south_padre_island_tx_navd88' => {
      begin_date  => '2012-05-04',
      end_date    => '2012-05-04',
      west_bound  => -97.65,
      east_bound  => -96.8,
      south_bound => 25.82,
      north_bound => 26.86
    },
    'elevation_spring_warrior_creek_fl_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -83.688,
      east_bound  => -83.647,
      south_bound => 29.910,
      north_bound => 29.949
    },
    'elevation_storter_bay_fl_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -81.379,
      east_bound  => -81.078,
      south_bound => 25.581,
      north_bound => 25.844
    },
    'elevation_suwanee_river_fl_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -83.227,
      east_bound  => -82.842,
      south_bound => 29.263,
      north_bound => 29.717
    },
    'elevation_tampa_bay_1_fl_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -82.544,
      east_bound  => -82.301,
      south_bound => 27.770,
      north_bound => 28.038
    },
    'elevation_tampa_bay_2_fl_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -82.779,
      east_bound  => -82.476,
      south_bound => 27.784,
      north_bound => 28.077
    },
    'elevation_tampa_bay_3_fl_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -82.903,
      east_bound  => -82.327,
      south_bound => 27.372,
      north_bound => 27.914
    },
    'elevation_tampa_bay_fl_mhw' => {
      begin_date  => '2014-09-11',
      end_date    => '2014-09-1',
      west_bound  => -83.41,
      east_bound  => -82.37,
      south_bound => 27.22,
      north_bound => 28.26
    },
    'elevation_tampa_bay_fl_navd88' => {
      begin_date  => '2014-09-11',
      end_date    => '2014-09-1',
      west_bound  => -83.41,
      east_bound  => -82.37,
      south_bound => 27.22,
      north_bound => 28.26
    },
    'elevation_ten_thousand_islands_marco_bay_fl_navd88' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -81.834,
      east_bound  => -81.241,
      south_bound => 25.730,
      north_bound => 26.082
    },
    'elevation_upper_laguna_madre_tx_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -97.611,
      east_bound  => -97.122,
      south_bound => 25.944,
      north_bound => 27.701
    },
    'elevation_waccasassa_river_fl_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -82.938,
      east_bound  => -82.642,
      south_bound => 29.052,
      north_bound => 29.252
    },
    'elevation_weeki_wachee_bay_fl_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -82.675,
      east_bound  => -82.617,
      south_bound => 28.511,
      north_bound => 28.560
    },
    'elevation_western_gulf_msl' => {
      begin_date  => '2001-07-20',
      end_date    => '2001-07-20',
      west_bound  => -107,
      east_bound  => -94,
      south_bound => 25,
      north_bound => 37
    },
    'elevation_whackup_creek_fl_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -83.344,
      east_bound  => -83.314,
      south_bound => 29.472,
      north_bound => 29.496
    },
    'elevation_whitewater_bay_fl_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -81.206,
      east_bound  => -80.793,
      south_bound => 25.086,
      north_bound => 25.502
    },
    'elevation_withlacoochee_bay_fl_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -82.818,
      east_bound  => -82.704,
      south_bound => 28.989,
      north_bound => 29.087
    },
    'elevation_yates_creek_fl_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -83.656,
      east_bound  => -83.647,
      south_bound => 29.894,
      north_bound => 29.903
    },

    'nadp.gridded.wet_nh4' => {
      begin_date  => '1985-01-01',
      end_date    => '2014-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'nadp.gridded.wet_no3' => {
      begin_date  => '1985-01-01',
      end_date    => '2014-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'nadp.gridded.wet_totn' => {
      begin_date  => '1985-01-01',
      end_date    => '2014-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'nadpcmaq.gridded.dry_nh3' => {
      begin_date  => '2000-01-01',
      end_date    => '2014-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'nadpcmaq.gridded.dry_nh4' => {
      begin_date  => '2000-01-01',
      end_date    => '2014-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'nadpcmaq.gridded.dry_no3' => {
      begin_date  => '2000-01-01',
      end_date    => '2014-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'nadpcmaq.gridded.dry_totn' => {
      begin_date  => '2000-01-01',
      end_date    => '2014-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'nadpcmaq.gridded.dry_totnp' => {
      begin_date  => '2000-01-01',
      end_date    => '2014-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'nadpcmaq.gridded.dry_totun' => {
      begin_date  => '2000-01-01',
      end_date    => '2014-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'nadpcmaq.gridded.dry_totunp' => {
      begin_date  => '2000-01-01',
      end_date    => '2014-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'nadpcmaq.gridded.dry_net_nh3' => {
      begin_date  => '2000-01-01',
      end_date    => '2014-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'nadpcmaq.gridded.dry_oxidized_totn' => {
      begin_date  => '2000-01-01',
      end_date    => '2014-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'nadpcmaq.gridded.dry_oxidized_totnp' => {
      begin_date  => '2000-01-01',
      end_date    => '2014-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'nadpcmaq.gridded.dry_reduced_totn' => {
      begin_date  => '2000-01-01',
      end_date    => '2014-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'nadpcmaq.gridded.dry_reduced_totnp' => {
      begin_date  => '2000-01-01',
      end_date    => '2014-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'nadpcmaq.gridded.wet_nh4' => {
      begin_date  => '2000-01-01',
      end_date    => '2014-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'nadpcmaq.gridded.wet_no3' => {
      begin_date  => '2000-01-01',
      end_date    => '2014-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'nadpcmaq.gridded.wet_totn' => {
      begin_date  => '2000-01-01',
      end_date    => '2014-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'nadpcmaq.gridded.oxidized_totn' => {
      begin_date  => '2000-01-01',
      end_date    => '2014-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'nadpcmaq.gridded.oxidized_totnp' => {
      begin_date  => '2000-01-01',
      end_date    => '2014-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'nadpcmaq.gridded.reduced_totn' => {
      begin_date  => '2000-01-01',
      end_date    => '2014-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'nadpcmaq.gridded.reduced_totnp' => {
      begin_date  => '2000-01-01',
      end_date    => '2014-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'nadpcmaq.gridded.total_totn' => {
      begin_date  => '2000-01-01',
      end_date    => '2014-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },

    'temperature_minimum' => {
      begin_date  => '1971-01-01',
      end_date    => '2000-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'temperature_maximum' => {
      begin_date  => '1971-01-01',
      end_date    => '2000-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'precipitation_mean' => {
      begin_date  => '1971-01-01',
      end_date    => '2000-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'monthly_temperature_minimum' => {
      begin_date  => '1971-01-01',
      end_date    => '2000-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'monthly_temperature_maximum' => {
      begin_date  => '1971-01-01',
      end_date    => '2000-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'monthly_precipitation_mean' => {
      begin_date  => '1971-01-01',
      end_date    => '2000-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'annual_temperature_minimum' => {
      begin_date  => '1971-01-01',
      end_date    => '2000-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'annual_temperature_maximum' => {
      begin_date  => '1971-01-01',
      end_date    => '2000-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'annual_precipitation_mean' => {
      begin_date  => '1971-01-01',
      end_date    => '2000-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'seagrass' => {
      begin_date  => '1979-01-01',
      end_date    => '2009-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'nitrogen.source_gridded.crops' => {
      begin_date  => '2006-01-01',
      end_date    => '2006-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'nitrogen.source_gridded.fertilizer' => {
      begin_date  => '2006-01-01',
      end_date    => '2006-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'nitrogen.source_gridded.manure' => {
      begin_date  => '2006-01-01',
      end_date    => '2006-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'wetlands' => {
      begin_date  => '1981-01-01',
      end_date    => '2011-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'population_iclus1' => {
      begin_date  => '2010-01-01',
      end_date    => '2090-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'population_iclus2' => {
      begin_date  => '2020-01-01',
      end_date    => '2090-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'land_use_iclus_watershed' => {
      begin_date  => '2000-01-01',
      end_date    => '2100-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'land_use_nlcd_2011_watershed' => {
      begin_date  => '2011-01-01',
      end_date    => '2011-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'land_use_nlcd_2006_watershed' => {
      begin_date  => '2006-01-01',
      end_date    => '2006-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'soil' => {
      begin_date  => '2012-01-01',
      end_date    => '2012-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'nitrogen.estuary_cmaq' => { # Old
      begin_date  => '2002-01-01',
      end_date    => '2008-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'dry_nh3.estuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'dry_natural_nh3.estuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'dry_nhx.estuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'dry_so2.estuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'dry_oxidized_nox.estuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'dry_oxidized_organic_nitrogen.estuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'dry_oxidized_no3.estuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'dry_oxidized_nitrogen.estuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'dry_reduced_nitrogen.estuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'dry_sulfur.estuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'chloride.estuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'total_nitrogen.estuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'oxidized_nitrogen.estuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'reduced_nitrogen.estuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'sulfur.estuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'wet_nh3.estuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'wet_chloride.estuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'wet_so4.estuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'wet_oxidized_no3.estuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'wet_oxidized_nox.estuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'wet_organic_nitrogen.estuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'wet_oxidized_nitrogen.estuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'wet_reduced_nitrogen.estuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'wet_sulfur.estuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'nitrogen.subestuary_cmaq' => { # Old
      begin_date  => '2002-01-01',
      end_date    => '2008-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'dry_nh3.subestuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'dry_natural_nh3.subestuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'dry_nhx.subestuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'dry_so2.subestuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'dry_oxidized_nox.subestuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'dry_oxidized_organic_nitrogen.subestuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'dry_oxidized_no3.subestuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'dry_oxidized_nitrogen.subestuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'dry_reduced_nitrogen.subestuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'dry_sulfur.subestuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'chloride.subestuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'total_nitrogen.subestuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'oxidized_nitrogen.subestuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'reduced_nitrogen.subestuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'sulfur.subestuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'wet_nh3.subestuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'wet_chloride.subestuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'wet_so4.subestuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'wet_oxidized_no3.subestuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'wet_oxidized_nox.subestuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'wet_organic_nitrogen.subestuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'wet_oxidized_nitrogen.subestuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'wet_reduced_nitrogen.subestuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'wet_sulfur.subestuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -85,
      east_bound  => -66,
      south_bound => 24,
      north_bound => 48
    },
    'nitrogen.load_estuary_sparrow_1992' => {
      begin_date  => '1992-01-01',
      end_date    => '1992-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'nitrogen.load_estuary_sparrow_2002_mrb5' => {
      begin_date  => '1992-01-01',
      end_date    => '1992-12-31',
      west_bound  => -97.3,
      east_bound  => -96.8,
      south_bound => 25.9,
      north_bound => 30.2
    },
    'nitrogen.load_estuary_sparrow_2011_mrb5' => {
      begin_date  => '2011-01-01',
      end_date    => '2011-12-31',
      west_bound  => -97.3,
      east_bound  => -96.8,
      south_bound => 25.9,
      north_bound => 30.2
    },
    'nitrogen.load_subestuary_sparrow_1992' => {
      begin_date  => '1992-01-01',
      end_date    => '1992-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'nitrogen.load_subestuary_sparrow_2002_mrb5' => {
      begin_date  => '1992-01-01',
      end_date    => '1992-12-31',
      west_bound  => -97.3,
      east_bound  => -96.8,
      south_bound => 25.9,
      north_bound => 30.2
    },
    'total_nitrogen.source_watershed_nadp' => {
      begin_date  => '2000-01-01',
      end_date    => '2014-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'wet_nitrogen.source_watershed_nadp' => {
      begin_date  => '2000-01-01',
      end_date    => '2014-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'dry_nitrogen.source_watershed_nadp' => {
      begin_date  => '2000-01-01',
      end_date    => '2014-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'oxidized_nitrogen.source_watershed_nadp' => {
      begin_date  => '2000-01-01',
      end_date    => '2014-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'reduced_nitrogen.source_watershed_nadp' => {
      begin_date  => '2000-01-01',
      end_date    => '2014-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'dry_oxidized_nitrogen.source_watershed_nadp' => {
      begin_date  => '2000-01-01',
      end_date    => '2014-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'dry_reduced_nitrogen.source_watershed_nadp' => {
      begin_date  => '2000-01-01',
      end_date    => '2014-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'wet_no3.source_watershed_nadp' => {
      begin_date  => '1994-01-01',
      end_date    => '2014-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'wet_nh4.source_watershed_nadp' => {
      begin_date  => '1994-01-01',
      end_date    => '2014-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'dry_no3.source_watershed_nadp' => {
      begin_date  => '2000-01-01',
      end_date    => '2014-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'dry_nh4.source_watershed_nadp' => {
      begin_date  => '2000-01-01',
      end_date    => '2014-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'nitrogen.source_watershed_cmaq' => { # Old
      begin_date  => '2002-01-01',
      end_date    => '2008-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'dry_nh3.source_watershed_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'dry_natural_nh3.source_watershed_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'dry_nhx.source_watershed_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'dry_so2.source_watershed_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'dry_oxidized_nox.source_watershed_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'dry_oxidized_organic_nitrogen.source_watershed_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'dry_oxidized_no3.source_watershed_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'dry_oxidized_nitrogen.source_watershed_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'dry_reduced_nitrogen.source_watershed_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'dry_sulfur.source_watershed_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'chloride.source_watershed_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'total_nitrogen.source_watershed_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'oxidized_nitrogen.source_watershed_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'reduced_nitrogen.source_watershed_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'sulfur.source_watershed_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'wet_nh3.source_watershed_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'wet_chloride.source_watershed_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'wet_so4.source_watershed_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'wet_oxidized_no3.source_watershed_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'wet_oxidized_nox.source_watershed_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'wet_organic_nitrogen.source_watershed_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'wet_oxidized_nitrogen.source_watershed_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'wet_reduced_nitrogen.source_watershed_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'wet_sulfur.source_watershed_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'nitrogen.source_watershed_point' => {
      begin_date  => '2007-01-01',
      end_date    => '2011-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'nitrogen.source_watershed_nonpoint' => {
      begin_date  => '2006-01-01',
      end_date    => '2006-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'phosphorus.source_watershed_point' => {
      begin_date  => '2007-01-01',
      end_date    => '2015-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'sediment_nca' => {
      begin_date  => '1990-01-01',
      end_date    => '2006-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'sediment_nca_1990-2006' => {
      begin_date  => '1990-01-01',
      end_date    => '2006-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'sediment_nca_2015' => {
      begin_date  => '2015-01-01',
      end_date    => '2015-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'sediment_calculated' => {
      begin_date  => '1990-01-01',
      end_date    => '2006-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'sediment_extracted' => {
      begin_date  => '1990-01-01',
      end_date    => '2006-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'sediment_parsed' => {
      begin_date  => '1990-01-01',
      end_date    => '2006-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'sediment_gridded' => {
      begin_date  => '2005-01-01',
      end_date    => '2014-12-31',
      west_bound  => -104,
      east_bound  => -79,
      south_bound => 24,
      north_bound => 38
    },
    'coastal_vulnerability' => {
      begin_date  => '1976-01-01',
      end_date    => '1995-12-31',
      west_bound  => -98,
      east_bound  => -80,
      south_bound => 25,
      north_bound => 31
    },
    'sensitivity_estuary' => {
      begin_date  => '2016-10-05',
      end_date    => '2016-10-05',
      west_bound  => -98,
      east_bound  => -80,
      south_bound => 25,
      north_bound => 31
    },
    'groundwater_contact_catchment' => {
      begin_date  => '2002-01-01',
      end_date    => '2002-12-31',
      west_bound  => -98,
      east_bound  => -80,
      south_bound => 25,
      north_bound => 31
    },
    'groundwater_contact_watershed' => {
      begin_date  => '2002-01-01',
      end_date    => '2002-12-31',
      west_bound  => -98,
      east_bound  => -80,
      south_bound => 25,
      north_bound => 31
    },
    'stream_discharge_nhd_line' => {
      begin_date  => '1952-01-01',
      end_date    => '2017-12-31',
      west_bound  => -85,
      east_bound  => -80,
      south_bound => 25,
      north_bound => 30
    },
    'stream_discharge_usgs_point' => {
      begin_date  => '1952-01-01',
      end_date    => '2017-12-31',
      west_bound  => -98,
      east_bound  => -81,
      south_bound => 27,
      north_bound => 32
    },
    'tide_point' => {
      begin_date  => '1950-01-01',
      end_date    => '2018-12-31',
      west_bound  => -98,
      east_bound  => -81,
      south_bound => 27,
      north_bound => 32
    },
    'tide_current' => {
      begin_date  => '2017-01-01',
      end_date    => '2019-12-31',
      west_bound  => -98,
      east_bound  => -80,
      south_bound => 25,
      north_bound => 30
    },
    'longshore_current' => {
      begin_date  => '2012-01-01',
      end_date    => '2020-12-31',
      west_bound  => -98,
      east_bound  => -80,
      south_bound => 25,
      north_bound => 30
    },
    'estuary_flushing' => {
      begin_date  => '1950-01-01',
      end_date    => '2020-12-31',
      west_bound  => -98,
      east_bound  => -80,
      south_bound => 25,
      north_bound => 30
    },
    'salinity_point' => {
      begin_date  => '1990-01-01',
      end_date    => '2018-12-31',
      west_bound  => -98,
      east_bound  => -80,
      south_bound => 24,
      north_bound => 31
    },
    'seagrass_point' => {
      begin_date  => '1978-01-01',
      end_date    => '2017-12-31',
      west_bound  => -98,
      east_bound  => -80,
      south_bound => 24,
      north_bound => 31
    },
    'seagrass_point2' => {
      begin_date  => '1978-01-01',
      end_date    => '2017-12-31',
      west_bound  => -98,
      east_bound  => -80,
      south_bound => 24,
      north_bound => 31
    },
    'seagrass_polygon' => {
      begin_date  => '1978-01-01',
      end_date    => '2017-12-31',
      west_bound  => -98,
      east_bound  => -80,
      south_bound => 24,
      north_bound => 31
    },
    'esat1_water_quality' => {
      begin_date  => '1984-01-01',
      end_date    => '2021-12-31',
      west_bound  =>  -124.02,
      east_bound  =>  -70.52,
      south_bound =>   24.91,
      north_bound =>   48.56
    },
    'esat2_water_quality' => {
      begin_date  => '2018-01-01',
      end_date    => '2020-12-31',
      west_bound  =>  -124.02,
      east_bound  =>  -70.52,
      south_bound =>   24.91,
      north_bound =>   48.56
    }
  },

  'pacific' => {
    'states' => {
      begin_date  => '2018-04-01',
      end_date    => '2018-04-01',
      west_bound  => -124.84897,
      east_bound  => -104.03969,
      south_bound =>   32.52883,
      north_bound =>   49.00243
    },
    'provinces' => {
      begin_date  => '2018-04-01',
      end_date    => '2018-04-01',
      west_bound  => -139.05220,
      east_bound  => -114.05376,
      south_bound =>   48.30892,
      north_bound =>   60.00006
    },
    'counties' => {
      begin_date  => '2018-04-01',
      end_date    => '2018-04-01',
      west_bound  => -124.84897,
      east_bound  => -104.03969,
      south_bound =>   32.52883,
      north_bound =>   49.00243
    },
    'cities' => {
      begin_date  => '2018-04-01',
      end_date    => '2018-04-01',
      west_bound  => -124.43730,
      east_bound  => -104.13282,
      south_bound =>   32.53973,
      north_bound =>   49.00209
    },
    'roads' => {
      begin_date  => '2018-04-01',
      end_date    => '2018-04-01',
      west_bound  => -124.21652,
      east_bound  => -104.04659,
      south_bound =>   32.54218,
      north_bound =>   49.00209
    },
    'tributaries' => {
      begin_date  => '2018-04-01',
      end_date    => '2018-04-01',
      west_bound  => -124.72442,
      east_bound  => -116.46121,
      south_bound =>   32.54973,
      north_bound =>   48.99995
    },
    'hucs' => {
      begin_date  => '2018-04-01',
      end_date    => '2018-04-01',
      west_bound  => -124.90219,
      east_bound  => -109.76248,
      south_bound =>   32.13300,
      north_bound =>   50.63471
    },
    'watersheds' => {
      begin_date  => '2018-04-01',
      end_date    => '2018-04-01',
      west_bound  => -124.71877,
      east_bound  => -109.76215,
      south_bound =>   32.16560,
      north_bound =>   52.88004
    },
    'estuaries' => {
      begin_date  => '2018-04-01',
      end_date    => '2018-04-01',
      west_bound  => -124.74974,
      east_bound  => -117.04971,
      south_bound =>   32.54969,
      north_bound =>   49.00225
    },
    'lakes' => {
      begin_date  => '2018-04-01',
      end_date    => '2018-04-01',
      west_bound  => -124.66945,
      east_bound  => -104.17377,
      south_bound =>   33.09148,
      north_bound =>   49.00093
    },
    'legislative_districts' => {
      begin_date  => '2018-04-01',
      end_date    => '2018-04-01',
      west_bound  => -124.84897,
      east_bound  => -104.03969,
      south_bound =>   32.52883,
      north_bound =>   49.00243
    },

    'nlcd_2011' => {
      begin_date  => '2011-01-01',
      end_date    => '2011-12-31',
      west_bound  => -128,
      east_bound  => -107,
      south_bound => 31,
      north_bound => 51
    },
    'nlcd_2006' => {
      begin_date  => '2006-01-01',
      end_date    => '2006-12-31',
      west_bound  => -128,
      east_bound  => -107,
      south_bound => 31,
      north_bound => 51
    },
    'nlcd_2001' => {
      begin_date  => '2001-01-01',
      end_date    => '2001-12-31',
      west_bound  => -128,
      east_bound  => -107,
      south_bound => 31,
      north_bound => 51
    },
    'nlcd_1992' => {
      begin_date  => '1992-01-01',
      end_date    => '1992-12-31',
      west_bound  => -128,
      east_bound  => -107,
      south_bound => 31,
      north_bound => 51
    },
    'nlcd_change_1992_2001' => {
      begin_date  => '1992-01-01',
      end_date    => '2001-12-31',
      west_bound  => -128,
      east_bound  => -107,
      south_bound => 31,
      north_bound => 51
    },
    'nlcd_change_2001_2006' => {
      begin_date  => '2001-01-01',
      end_date    => '2006-12-31',
      west_bound  => -128,
      east_bound  => -107,
      south_bound => 31,
      north_bound => 51
    },
    'ccap_2006' => {
      begin_date  => '2006-01-01',
      end_date    => '2006-12-31',
      west_bound  => -128,
      east_bound  => -114,
      south_bound => 31,
      north_bound => 50
    },
    'ccap_2001' => {
      begin_date  => '2001-01-01',
      end_date    => '2001-12-31',
      west_bound  => -128,
      east_bound  => -114,
      south_bound => 31,
      north_bound => 50
    },
    'decadal_iclus_land_use_2000_2100' => {
      begin_date  => '2000-01-01',
      end_date    => '2100-12-31',
      west_bound  => -128,
      east_bound  => -107,
      south_bound => 31,
      north_bound => 51
    },
    'iclus_land_use_watershed' => {
      begin_date  => '2000-01-01',
      end_date    => '2100-12-31',
      west_bound  => -128,
      east_bound  => -107,
      south_bound => 31,
      north_bound => 51
    },
    'nlcd_2011_land_use_watershed' => {
      begin_date  => '2011-01-01',
      end_date    => '2011-12-31',
      west_bound  => -128,
      east_bound  => -107,
      south_bound => 31,
      north_bound => 51
    },
    'nlcd_2006_land_use_watershed' => {
      begin_date  => '2006-01-01',
      end_date    => '2006-12-31',
      west_bound  => -128,
      east_bound  => -107,
      south_bound => 31,
      north_bound => 51
    },
    'impervious_2006' => {
      begin_date  => '2001-01-01',
      end_date    => '2001-12-31',
      west_bound  => -128,
      east_bound  => -107,
      south_bound => 31,
      north_bound => 51
    },
    'impervious_2001' => {
      begin_date  => '2001-01-01',
      end_date    => '2001-12-31',
      west_bound  => -128,
      east_bound  => -107,
      south_bound => 31,
      north_bound => 51
    },
    'impervious_change_2001_2006' => {
      begin_date  => '2001-01-01',
      end_date    => '2006-12-31',
      west_bound  => -128,
      east_bound  => -107,
      south_bound => 31,
      north_bound => 51
    },
    'decadal_iclus_a1_impervious_2010_2100' => {
      begin_date  => '2010-01-01',
      end_date    => '2100-12-31',
      west_bound  => -128,
      east_bound  => -107,
      south_bound => 31,
      north_bound => 51
    },
    'decadal_iclus_a2_impervious_2010_2100' => {
      begin_date  => '2010-01-01',
      end_date    => '2100-12-31',
      west_bound  => -128,
      east_bound  => -107,
      south_bound => 31,
      north_bound => 51
    },
    'decadal_iclus_b1_impervious_2010_2100' => {
      begin_date  => '2010-01-01',
      end_date    => '2100-12-31',
      west_bound  => -128,
      east_bound  => -107,
      south_bound => 31,
      north_bound => 51
    },
    'decadal_iclus_b2_impervious_2010_2100' => {
      begin_date  => '2010-01-01',
      end_date    => '2100-12-31',
      west_bound  => -128,
      east_bound  => -107,
      south_bound => 31,
      north_bound => 51
    },
    'decadal_iclus_bc_impervious_2010_2100' => {
      begin_date  => '2010-01-01',
      end_date    => '2100-12-31',
      west_bound  => -128,
      east_bound  => -107,
      south_bound => 31,
      north_bound => 51
    },

    'elevation' => {
      begin_date  => '2004-01-01',
      end_date    => '2004-12-31',
      west_bound  => -180,
      east_bound  => -66,
      south_bound => 18,
      north_bound => 66
    },

    # From: https://www.ngdc.noaa.gov/dem/squareCellGrid/search

    'elevation_agua_hedionda_creek_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -117.348,
      east_bound  => -117.305,
      south_bound => 33.131,
      north_bound => 33.155
    },
    'elevation_albion_river_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -123.791,
      east_bound  => -123.702,
      south_bound => 39.201,
      north_bound => 39.263
    },
    'elevation_alder_creek_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -123.699,
      east_bound  => -123.689,
      south_bound => 38.998,
      north_bound => 39.012
    },
    'elevation_aliso_canyon_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -117.444,
      east_bound  => -117.439,
      south_bound => 33.264,
      north_bound => 33.268
    },
    'elevation_aliso_creek_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -117.755,
      east_bound  => -117.747,
      south_bound => 33.510,
      north_bound => 33.513
    },
    'elevation_alsea_bay_or_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -124.115,
      east_bound  => -123.872,
      south_bound => 44.351,
      north_bound => 44.494
    },
    'elevation_alva_paul_creek_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -120.869,
      east_bound  => -120.866,
      south_bound => 35.397,
      north_bound => 35.400
    },
    'elevation_andree_clark_bird_refuge_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -119.665,
      east_bound  => -119.656,
      south_bound => 34.417,
      north_bound => 34.423
    },
    'elevation_arena_cove_ca_mhw' => {
      begin_date  => '2009-11-30',
      end_date    => '2009-11-30',
      west_bound  => -124.43,
      east_bound  => -123.43,
      south_bound => 38.4,
      north_bound => 39.4
    },
    'elevation_arroyo_del_corral_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -121.288,
      east_bound  => -121.284,
      south_bound => 35.682,
      north_bound => 35.685
    },
    'elevation_arroyo_del_puerto_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -121.190,
      east_bound  => -121.188,
      south_bound => 35.643,
      north_bound => 35.645
    },
    'elevation_arroyo_burro_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -119.744,
      east_bound  => -119.740,
      south_bound => 34.402,
      north_bound => 34.405
    },
    'elevation_arroyo_de_la_cruz_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -121.312,
      east_bound  => -121.301,
      south_bound => 35.706,
      north_bound => 35.713
    },
    'elevation_arroyo_grande_creek_lagoon_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -120.636,
      east_bound  => -120.621,
      south_bound => 35.094,
      north_bound => 35.112
    },
    'elevation_arroyo_laguna_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -121.221,
      east_bound  => -121.218,
      south_bound => 35.650,
      north_bound => 35.655
    },
    'elevation_arthur_kill_nj_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -74.324,
      east_bound  => -74.118,
      south_bound => 40.478,
      north_bound => 40.676
    },
    'elevation_arthur_kill_ny_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -74.324,
      east_bound  => -74.118,
      south_bound => 40.478,
      north_bound => 40.676
    },
    'elevation_assawoman_bay_de_mhw' => {
      begin_date  => '2020-06-22',
      end_date    => '2020-06-22',
      west_bound  => -75.991,
      east_bound  => -74.991,
      south_bound => 38.297,
      north_bound => 38.562
    },
    'elevation_assawoman_bay_md_mhw' => {
      begin_date  => '2020-06-22',
      end_date    => '2020-06-22',
      west_bound  => -75.991,
      east_bound  => -74.991,
      south_bound => 38.297,
      north_bound => 38.562
    },
    'elevation_astoria_or_mhw' => {
      begin_date  => '2016-06-18',
      end_date    => '2016-06-18',
      west_bound  => -124.59,
      east_bound  => -123.6,
      south_bound => 45.9,
      north_bound => 47.16
    },
    'elevation_baldwin_creek_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -122.126,
      east_bound  => -122.121,
      south_bound => 36.966,
      north_bound => 36.971
    },
    'elevation_bc_msl' => {
      begin_date  => '2013-07-09',
      end_date    => '2013-07-09',
      west_bound  => -68.7,
      east_bound  => -67.62,
      south_bound => 43.69,
      north_bound => 44.69
    },
    'elevation_bear_river_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -124.397,
      east_bound  => -124.382,
      south_bound => 40.468,
      north_bound => 40.484
    },
    'elevation_beaver_creek_or_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -124.089,
      east_bound  => -124.024,
      south_bound => 44.497,
      north_bound => 44.534
    },
    'elevation_bell_canyon_creek_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -119.914,
      east_bound  => -119.912,
      south_bound => 34.428,
      north_bound => 34.431
    },
    'elevation_bellingham_bay_wa_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -122.755,
      east_bound  => -122.404,
      south_bound => 48.550,
      north_bound => 48.824
    },
    'elevation_berry_creek_or_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -124.126,
      east_bound  => -124.124,
      south_bound => 44.090,
      north_bound => 44.092
    },
    'elevation_big_creek_1_or_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -124.093,
      east_bound  => -124.080,
      south_bound => 44.367,
      north_bound => 44.374
    },
    'elevation_big_creek_or_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -124.117,
      east_bound  => -124.115,
      south_bound => 44.175,
      north_bound => 44.177
    },
    'elevation_big_lagoon_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -124.153,
      east_bound  => -124.082,
      south_bound => 41.150,
      north_bound => 41.218
    },
    'elevation_big_river_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -123.798,
      east_bound  => -123.723,
      south_bound => 39.278,
      north_bound => 39.320
    },
    'elevation_big_sur_river_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -121.861,
      east_bound  => -121.855,
      south_bound => 36.279,
      north_bound => 36.284
    },
    'elevation_birch_bay_wa_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -122.840,
      east_bound  => -122.723,
      south_bound => 48.880,
      north_bound => 48.958
    },
    'elevation_bob_creek_or_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -124.114,
      east_bound  => -124.003,
      south_bound => 44.213,
      north_bound => 44.267
    },
    'elevation_bodega_bay_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -123.102,
      east_bound  => -122.923,
      south_bound => 38.201,
      north_bound => 38.358
    },
    'elevation_bolinas_lagoon_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -122.699,
      east_bound  => -122.644,
      south_bound => 37.892,
      north_bound => 37.943
    },
    'elevation_bolsa_chica_ecological_reserve_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -118.073,
      east_bound  => -117.985,
      south_bound => 33.670,
      north_bound => 33.737
    },
    'elevation_boone_creek_wa_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -124.194,
      east_bound  => -124.184,
      south_bound => 47.157,
      north_bound => 47.162
    },
    'elevation_brush_creek_1_or_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -124.449,
      east_bound  => -124.444,
      south_bound => 42.685,
      north_bound => 42.688
    },
    'elevation_brush_creek_2_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -123.720,
      east_bound  => -123.697,
      south_bound => 38.961,
      north_bound => 38.988
    },
    'elevation_buena_vista_lagoon_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -117.363,
      east_bound  => -117.337,
      south_bound => 33.164,
      north_bound => 33.181
    },
    'elevation_calera_creek_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -122.497,
      east_bound  => -122.492,
      south_bound => 37.611,
      north_bound => 37.614
    },
    'elevation_camp_creek_wa_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -124.330,
      east_bound  => -124.324,
      south_bound => 47.397,
      north_bound => 47.401
    },
    'elevation_canada_de_la_gaviota_creek_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -120.229,
      east_bound  => -120.225,
      south_bound => 34.470,
      north_bound => 34.474
    },
    'elevation_canada_del_refugio_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -120.081,
      east_bound  => -120.055,
      south_bound => 34.462,
      north_bound => 34.500
    },
    'elevation_cape_creek_or_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -124.124,
      east_bound  => -124.121,
      south_bound => 44.133,
      north_bound => 44.135
    },
    'elevation_carmel_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -121.988,
      east_bound  => -121.901,
      south_bound => 36.508,
      north_bound => 36.580
    },
    'elevation_carpinteria_creek_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -119.521,
      east_bound  => -119.517,
      south_bound => 34.390,
      north_bound => 34.392
    },
    'elevation_carpinteria_salt_marsh_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -119.548,
      east_bound  => -119.518,
      south_bound => 34.391,
      north_bound => 34.410
    },
    'elevation_casper_creek_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -123.819,
      east_bound  => -123.812,
      south_bound => 39.358,
      north_bound => 39.363
    },
    'elevation_cayucos_creek_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -120.912,
      east_bound  => -120.903,
      south_bound => 35.446,
      north_bound => 35.451
    },
    'elevation_cedar_creek_wa_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -124.417,
      east_bound  => -124.412,
      south_bound => 47.710,
      north_bound => 47.713
    },
    'elevation_northern_or_mhw' => {
      begin_date  => '2004-10-06',
      end_date    => '2004-10-06',
      west_bound  => -124.5,
      east_bound  => -123.5,
      south_bound => 45.89,
      north_bound => 47.36
    },
    'elevation_central_ca_mhw' => {
      begin_date  => '2010-07-20',
      end_date    => '2010-07-20',
      west_bound  => -124,
      east_bound  => -121.85,
      south_bound => 37,
      north_bound => 38.75
    },
    'elevation_central_ca_navd88' => {
      begin_date  => '2010-07-20',
      end_date    => '2010-07-20',
      west_bound  => -124,
      east_bound  => -121.85,
      south_bound => 37,
      north_bound => 38.75
    },
    'elevation_central_or_mhw' => {
      begin_date  => '2008-09-03',
      end_date    => '2008-09-03',
      west_bound  => -124.63,
      east_bound  => -123.81,
      south_bound => 43.7,
      north_bound => 45.14
    },
    'elevation_central_pacific_msl' => {
      begin_date  => '2003-04-21',
      end_date    => '2003-04-21',
      west_bound  => -127,
      east_bound  => -117,
      south_bound => 37,
      north_bound => 44
    },
    'elevation_chadbourne_gulch_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -123.783,
      east_bound  => -123.781,
      south_bound => 39.613,
      north_bound => 39.615
    },
    'elevation_channel_islands_harbor_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -119.264,
      east_bound  => -119.208,
      south_bound => 34.151,
      north_bound => 34.212
    },
    'elevation_chetco_river_or_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -124.283,
      east_bound  => -124.218,
      south_bound => 42.042,
      north_bound => 42.079
    },
    'elevation_china_creek_or_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -124.118,
      east_bound  => -124.083,
      south_bound => 44.138,
      north_bound => 44.174
    },
    'elevation_clam_beach_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -124.121,
      east_bound  => -124.118,
      south_bound => 40.986,
      north_bound => 40.988
    },
    'elevation_cockelburr_canyon_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -117.433,
      east_bound  => -117.430,
      south_bound => 33.250,
      north_bound => 33.252
    },
    'elevation_columbia_river_or_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -124.406,
      east_bound  => -121.637,
      south_bound => 45.134,
      north_bound => 46.751
    },
    'elevation_columbia_river_wa_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -124.406,
      east_bound  => -121.637,
      south_bound => 45.134,
      north_bound => 46.751
    },
    'elevation_conner_creek_wa_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -124.195,
      east_bound  => -124.153,
      south_bound => 47.054,
      north_bound => 47.123
    },
    'elevation_coos_bay_or_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -124.436,
      east_bound  => -124.046,
      south_bound => 43.217,
      north_bound => 43.523
    },
    'elevation_copalis_river_wa_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -124.221,
      east_bound  => -124.132,
      south_bound => 47.083,
      north_bound => 47.166
    },
    'elevation_coquille_river_or_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -124.469,
      east_bound  => -124.180,
      south_bound => 43.064,
      north_bound => 43.252
    },
    'elevation_cottaneva_creek_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -123.829,
      east_bound  => -123.824,
      south_bound => 39.733,
      north_bound => 39.737
    },
    'elevation_crescent_city_harbor_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -124.202,
      east_bound  => -124.178,
      south_bound => 41.735,
      north_bound => 41.755
    },
    'elevation_crescent_city_ca_mhw' => {
      begin_date  => '2010-07-14',
      end_date    => '2010-07-14',
      west_bound  => -125.25,
      east_bound  => -123.88,
      south_bound => 41.42,
      north_bound => 42.53
    },
    'elevation_crescent_city_ca_navd88' => {
      begin_date  => '2010-07-14',
      end_date    => '2010-07-14',
      west_bound  => -125.25,
      east_bound  => -123.88,
      south_bound => 41.42,
      north_bound => 42.53
    },
    'elevation_cummins_creek_or_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -124.107,
      east_bound  => -124.020,
      south_bound => 44.241,
      north_bound => 44.283
    },
    'elevation_cut_creek_or_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -124.402,
      east_bound  => -124.386,
      south_bound => 43.183,
      north_bound => 43.189
    },
    'elevation_dairy_gulch_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -122.092,
      east_bound  => -122.090,
      south_bound => 36.954,
      north_bound => 36.956
    },
    'elevation_daley_lake_or_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -123.977,
      east_bound  => -123.969,
      south_bound => 45.140,
      north_bound => 45.146
    },
    'elevation_dana_point_harbor_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -117.709,
      east_bound  => -117.688,
      south_bound => 33.451,
      north_bound => 33.466
    },
    'elevation_depoe_bay_or_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -124.064,
      east_bound  => -124.057,
      south_bound => 44.806,
      north_bound => 44.811
    },
    'elevation_devereux_lagoon_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -119.884,
      east_bound  => -119.871,
      south_bound => 34.408,
      north_bound => 34.418
    },
    'elevation_devils_lake_or_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -124.019,
      east_bound  => -124.016,
      south_bound => 44.967,
      north_bound => 44.969
    },
    'elevation_drakes_estero_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -123.021,
      east_bound  => -122.774,
      south_bound => 37.949,
      north_bound => 38.128
    },
    'elevation_drayton_harbor_wa_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -122.814,
      east_bound  => -122.687,
      south_bound => 48.932,
      north_bound => 49.014
    },
    'elevation_duck_creek_wa_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -124.327,
      east_bound  => -124.323,
      south_bound => 47.387,
      north_bound => 47.390
    },
    'elevation_ecola_creek_or_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -123.969,
      east_bound  => -123.949,
      south_bound => 45.892,
      north_bound => 45.906
    },
    'elevation_eel_river_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -124.341,
      east_bound  => -124.227,
      south_bound => 40.598,
      north_bound => 40.672
    },
    'elevation_elk_creek_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -123.708,
      east_bound  => -123.676,
      south_bound => 39.085,
      north_bound => 39.109
    },
    'elevation_elk_creek_wa_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -124.203,
      east_bound  => -124.191,
      south_bound => 47.188,
      north_bound => 47.194
    },
    'elevation_elk_river_2_wa_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -124.026,
      east_bound  => -123.967,
      south_bound => 46.823,
      north_bound => 46.856
    },
    'elevation_elk_river_3_or_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -124.559,
      east_bound  => -124.471,
      south_bound => 42.735,
      north_bound => 42.818
    },
    'elevation_espa_lagoon_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -124.078,
      east_bound  => -124.070,
      south_bound => 41.355,
      north_bound => 41.360
    },
    'elevation_euchre_creek_or_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -124.396,
      east_bound  => -124.382,
      south_bound => 42.555,
      north_bound => 42.564
    },
    'elevation_eureka_ca_mhw' => {
      begin_date  => '2009-11-30',
      end_date    => '2009-11-30',
      west_bound  => -125.01,
      east_bound  => -124,
      south_bound => 40.27,
      north_bound => 41.42
    },
    'elevation_eureka_ca_navd88' => {
      begin_date  => '2009-11-30',
      end_date    => '2009-11-30',
      west_bound  => -125.01,
      east_bound  => -124,
      south_bound => 40.27,
      north_bound => 41.42
    },
    'elevation_fern_canyon_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -124.076,
      east_bound  => -124.056,
      south_bound => 41.397,
      north_bound => 41.415
    },
    'elevation_floras_creek_2_or_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -124.505,
      east_bound  => -124.490,
      south_bound => 42.910,
      north_bound => 42.921
    },
    'elevation_fort_bragg_ca_mhw' => {
      begin_date  => '2012-01-15',
      end_date    => '2012-01-15',
      west_bound  => -125.05,
      east_bound  => -123.72,
      south_bound => 39.38,
      north_bound => 40.28
    },
    'elevation_fort_bragg_ca_navd88' => {
      begin_date  => '2012-01-15',
      end_date    => '2012-01-15',
      west_bound  => -125.05,
      east_bound  => -123.72,
      south_bound => 39.38,
      north_bound => 40.28
    },
      'elevation_frenchmans_creek_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -122.452,
      east_bound  => -122.445,
      south_bound => 37.479,
      north_bound => 37.484
    },
    'elevation_french_canyon_2_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -117.442,
      east_bound  => -117.434,
      south_bound => 33.262,
      north_bound => 33.266
    },
    'elevation_garibaldi_or_mhw' => {
      begin_date  => '2007-08-22',
      end_date    => '2007-08-22',
      west_bound  => -124.5,
      east_bound  => -123.7,
      south_bound => 45.1,
      north_bound => 45.95
    },
    'elevation_garibaldi_or_navd88' => {
      begin_date  => '2007-08-22',
      end_date    => '2007-08-22',
      west_bound  => -124.5,
      east_bound  => -123.7,
      south_bound => 45.1,
      north_bound => 45.95
    },
    'elevation_garcia_river_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -123.734,
      east_bound  => -123.718,
      south_bound => 38.941,
      north_bound => 38.958
    },
    'elevation_garrapata_creek_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -121.916,
      east_bound  => -121.914,
      south_bound => 36.416,
      north_bound => 36.418
    },
    'elevation_gazos_creek_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -122.363,
      east_bound  => -122.358,
      south_bound => 37.160,
      north_bound => 37.166
    },
    'elevation_goleta_slough_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -119.866,
      east_bound  => -119.818,
      south_bound => 34.409,
      north_bound => 34.438
    },
    'elevation_goodman_creek_wa_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -124.516,
      east_bound  => -124.504,
      south_bound => 47.822,
      north_bound => 47.830
    },
    'elevation_grays_harbor_wa_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -124.317,
      east_bound  => -123.397,
      south_bound => 46.595,
      north_bound => 47.180
    },
    'elevation_greenwood_creek_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -123.719,
      east_bound  => -123.713,
      south_bound => 39.124,
      north_bound => 39.128
    },
    'elevation_gualala_river_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -123.539,
      east_bound  => -123.501,
      south_bound => 38.752,
      north_bound => 38.779
    },
    'elevation_guthrie_creek_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -124.366,
      east_bound  => -124.362,
      south_bound => 40.540,
      north_bound => 40.543
    },
    'elevation_half_moon_bay_airport_lagoon_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -122.498,
      east_bound  => -122.493,
      south_bound => 37.500,
      north_bound => 37.505
    },
    'elevation_hare_creek_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -123.813,
      east_bound  => -123.803,
      south_bound => 39.416,
      north_bound => 39.420
    },
    'elevation_hoh_river_wa_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -124.498,
      east_bound  => -124.396,
      south_bound => 47.716,
      north_bound => 47.783
    },
    'elevation_honda_creek_lagoon_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -120.637,
      east_bound  => -120.635,
      south_bound => 34.607,
      north_bound => 34.609
    },
    'elevation_hubbard_creek_or_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -124.479,
      east_bound  => -124.475,
      south_bound => 42.733,
      north_bound => 42.737
    },
    'elevation_humboldt_bay_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -124.346,
      east_bound  => -124.001,
      south_bound => 40.633,
      north_bound => 40.916
    },
    'elevation_humboldt_ca_mhw' => {
      begin_date  => '2004-07-02',
      end_date    => '2004-07-02',
      west_bound  => -124.4,
      east_bound  => -124,
      south_bound => 40.55,
      north_bound => 41
    },
    'elevation_huntington_channel_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -118.009,
      east_bound  => -117.936,
      south_bound => 33.622,
      north_bound => 33.692
    },
    'elevation_inglenook_fen_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -123.781,
      east_bound  => -123.772,
      south_bound => 39.516,
      north_bound => 39.527
    },
    'elevation_irish_gulch_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -123.692,
      east_bound  => -123.687,
      south_bound => 39.019,
      north_bound => 39.022
    },
    'elevation_jackass_creek_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -123.917,
      east_bound  => -123.913,
      south_bound => 39.880,
      north_bound => 39.883
    },
    'elevation_jalama_creek_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -120.504,
      east_bound  => -120.499,
      south_bound => 34.511,
      north_bound => 34.513
    },
    'elevation_joe_creek_wa_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -124.206,
      east_bound  => -124.182,
      south_bound => 47.202,
      north_bound => 47.212
    },
    'elevation_johnson_creek_3_or_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -124.433,
      east_bound  => -124.360,
      south_bound => 43.044,
      north_bound => 43.105
    },
    'elevation_kalaloch_creek_wa_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -124.377,
      east_bound  => -124.368,
      south_bound => 47.605,
      north_bound => 47.610
    },
    'elevation_klamath_river_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -124.094,
      east_bound  => -123.989,
      south_bound => 41.500,
      north_bound => 41.565
    },
    'elevation_lagoon_creek_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -124.105,
      east_bound  => -124.091,
      south_bound => 41.583,
      north_bound => 41.598
    },
    'elevation_laguna_creek_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -122.157,
      east_bound  => -122.151,
      south_bound => 36.981,
      north_bound => 36.986
    },
    'elevation_lake_cleone_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -123.800,
      east_bound  => -123.789,
      south_bound => 39.485,
      north_bound => 39.491
    },
    'elevation_lake_earl_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -124.240,
      east_bound  => -124.152,
      south_bound => 41.782,
      north_bound => 41.867
    },
    'elevation_lake_lucerne_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -122.410,
      east_bound  => -122.392,
      south_bound => 37.220,
      north_bound => 37.229
    },
    'elevation_lake_lytle_or_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -123.944,
      east_bound  => -123.934,
      south_bound => 45.629,
      north_bound => 45.638
    },
    'elevation_lake_washington_wa_llw' => {
      begin_date  => '2010-01-05',
      end_date    => '2010-01-05',
      west_bound  => -122.3,
      east_bound  => -122.16,
      south_bound => 47.45,
      north_bound => 47.78
    },
    'elevation_la_push_wa_mhw' => {
      begin_date  => '2007-07-26',
      end_date    => '2007-07-26',
      west_bound  => -125.05,
      east_bound  => -124.35,
      south_bound => 47.55,
      north_bound => 48.45
    },
    'elevation_las_flores_creek_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -117.466,
      east_bound  => -117.460,
      south_bound => 33.289,
      north_bound => 33.292
    },
    'elevation_las_pulgas_creek_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -117.468,
      east_bound  => -117.464,
      south_bound => 33.292,
      north_bound => 33.295
    },
    'elevation_little_creek_1_or_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -124.062,
      east_bound  => -124.057,
      south_bound => 44.662,
      north_bound => 44.670
    },
    'elevation_little_creek_3_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -123.794,
      east_bound  => -123.783,
      south_bound => 39.273,
      north_bound => 39.281
    },
    'elevation_little_pico_creek_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -121.165,
      east_bound  => -121.162,
      south_bound => 35.633,
      north_bound => 35.635
    },
    'elevation_little_river_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -124.118,
      east_bound  => -124.100,
      south_bound => 41.013,
      north_bound => 41.033
    },
    'elevation_little_salmon_creek_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -123.770,
      east_bound  => -123.765,
      south_bound => 39.214,
      north_bound => 39.217
    },
    'elevation_little_sur_lagoon_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -121.896,
      east_bound  => -121.885,
      south_bound => 36.330,
      north_bound => 36.337
    },
    'elevation_loma_alta_marsh_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -117.370,
      east_bound  => -117.365,
      south_bound => 33.176,
      north_bound => 33.180
    },
    'elevation_lombardi_creek_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -122.114,
      east_bound  => -122.111,
      south_bound => 36.962,
      north_bound => 36.964
    },
    'elevation_los_penasquitos_lagoon_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -117.265,
      east_bound  => -117.240,
      south_bound => 32.920,
      north_bound => 32.941
    },
    'elevation_lunada_bay_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -118.431,
      east_bound  => -118.420,
      south_bound => 33.765,
      north_bound => 33.775
    },
    'elevation_mad_river_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -124.148,
      east_bound  => -124.047,
      south_bound => 40.892,
      north_bound => 40.985
    },
    'elevation_mugu_lagoon_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -119.148,
      east_bound  => -119.057,
      south_bound => 34.079,
      north_bound => 34.147
    },
    'elevation_makah_bay_wa_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -124.714,
      east_bound  => -124.601,
      south_bound => 48.283,
      north_bound => 48.370
    },
    'elevation_malibu_lagoon_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -118.687,
      east_bound  => -118.678,
      south_bound => 34.031,
      north_bound => 34.040
    },
    'elevation_malpaso_creek_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -121.939,
      east_bound  => -121.873,
      south_bound => 36.451,
      north_bound => 36.494
    },
    'elevation_manchester_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -123.706,
      east_bound  => -123.693,
      south_bound => 38.986,
      north_bound => 38.999
    },
    'elevation_mattole_river_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -124.360,
      east_bound  => -124.307,
      south_bound => 40.285,
      north_bound => 40.312
    },
    'elevation_mcnutt_gulch_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -124.364,
      east_bound  => -124.360,
      south_bound => 40.351,
      north_bound => 40.354
    },
    'elevation_misson_bay_lagoon_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -117.273,
      east_bound  => -117.182,
      south_bound => 32.743,
      north_bound => 32.817
    },
    'elevation_misson_creek_lagoon_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -119.693,
      east_bound  => -119.684,
      south_bound => 34.411,
      north_bound => 34.419
    },
    'elevation_mitchell_creek_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -123.818,
      east_bound  => -123.813,
      south_bound => 39.390,
      north_bound => 39.394
    },
    'elevation_moclips_river_wa_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -124.222,
      east_bound  => -124.198,
      south_bound => 47.235,
      north_bound => 47.252
    },
    'elevation_monterey_bay_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -122.188,
      east_bound  => -121.606,
      south_bound => 36.533,
      north_bound => 37.052
    },
    'elevation_monterey_ca_mhw' => {
      begin_date  => '2012-01-13',
      end_date    => '2012-01-13',
      west_bound  => -122.42,
      east_bound  => -121.28,
      south_bound => 35.7,
      north_bound => 37.2
    },
    'elevation_monterey_ca_navd88' => {
      begin_date  => '2012-01-13',
      end_date    => '2012-01-13',
      west_bound  => -122.42,
      east_bound  => -121.28,
      south_bound => 35.7,
      north_bound => 37.2
    },
    'elevation_moran_lake_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -121.979,
      east_bound  => -121.974,
      south_bound => 36.955,
      north_bound => 36.960
    },
    'elevation_morro_bay_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -120.890,
      east_bound  => -120.794,
      south_bound => 35.299,
      north_bound => 35.384
    },
    'elevation_morro_creek_lagoon_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -120.865,
      east_bound  => -120.860,
      south_bound => 35.375,
      north_bound => 35.380
    },
    'elevation_mosquito_creek_wa_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -124.484,
      east_bound  => -124.478,
      south_bound => 47.798,
      north_bound => 47.802
    },
    'elevation_mussel_creek_or_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -124.401,
      east_bound  => -124.398,
      south_bound => 42.614,
      north_bound => 42.617
    },
    'elevation_myers_creek_or_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -124.413,
      east_bound  => -124.409,
      south_bound => 42.306,
      north_bound => 42.308
    },
    'elevation_navarro_river_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -123.771,
      east_bound  => -123.701,
      south_bound => 39.166,
      north_bound => 39.211
    },
    'elevation_nc9_crescent_beach_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -124.153,
      east_bound  => -124.147,
      south_bound => 41.722,
      north_bound => 41.728
    },
    'elevation_necanicum_river_estuary_or_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -124.011,
      east_bound  => -123.852,
      south_bound => 45.975,
      north_bound => 46.137
    },
    'elevation_nehalem_bay_or_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -123.973,
      east_bound  => -123.820,
      south_bound => 45.639,
      north_bound => 45.758
    },
    'elevation_neskowin_creek_or_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -123.988,
      east_bound  => -123.977,
      south_bound => 45.095,
      north_bound => 45.107
    },
    'elevation_nestucca_bay_or_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -124.004,
      east_bound  => -123.893,
      south_bound => 45.132,
      north_bound => 45.225
    },
    'elevation_netarts_bay_or_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -124.000,
      east_bound  => -123.900,
      south_bound => 45.358,
      north_bound => 45.451
    },
    'elevation_new_river_2_or_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -124.495,
      east_bound  => -124.404,
      south_bound => 42.975,
      north_bound => 43.058
    },
    'elevation_newport_harbor_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -117.959,
      east_bound  => -117.832,
      south_bound => 33.572,
      north_bound => 33.677
    },
    'elevation_newport_or_navd88' => {
      begin_date  => '2004-09-28',
      end_date    => '2004-09-28',
      west_bound  => -124.15,
      east_bound  => -124,
      south_bound => 44.54,
      north_bound => 44.7
    },
    'elevation_northern_ca_msl' => {
      begin_date  => '2005-11-30',
      end_date    => '2005-11-30',
      west_bound  => -125,
      east_bound  => -124,
      south_bound => 41,
      north_bound => 44.5
    },
    'elevation_northwest_pacific_msl' => {
      begin_date  => '2003-04-21',
      end_date    => '2003-04-21',
      west_bound  => -127,
      east_bound  => -116,
      south_bound => 44,
      north_bound => 49
    },
    'elevation_noyo_river_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -123.812,
      east_bound  => -123.750,
      south_bound => 39.421,
      north_bound => 39.446
    },
    'elevation_oceanside_harbor_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -117.410,
      east_bound  => -117.385,
      south_bound => 33.201,
      north_bound => 33.223
    },
    'elevation_old_creek_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -120.888,
      east_bound  => -120.886,
      south_bound => 35.433,
      north_bound => 35.436
    },
    'elevation_ormond_beach_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -119.194,
      east_bound  => -119.165,
      south_bound => 34.128,
      north_bound => 34.155
    },
    'elevation_ossagon_creek_2_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -124.075,
      east_bound  => -124.051,
      south_bound => 41.415,
      north_bound => 41.448
    },
    'elevation_ozette_river_wa_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -124.710,
      east_bound  => -124.701,
      south_bound => 48.177,
      north_bound => 48.183
    },
    'elevation_padilla_bay_wa_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -122.663,
      east_bound  => -122.345,
      south_bound => 48.410,
      north_bound => 48.620
    },
    'elevation_pebble_beach_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -124.233,
      east_bound  => -124.221,
      south_bound => 41.762,
      north_bound => 41.771
    },
    'elevation_pescadero_creek_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -122.415,
      east_bound  => -122.395,
      south_bound => 37.248,
      north_bound => 37.275
    },
    'elevation_pico_creek_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -121.150,
      east_bound  => -121.146,
      south_bound => 35.615,
      north_bound => 35.618
    },
    'elevation_pilarcitos_creek_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -122.449,
      east_bound  => -122.438,
      south_bound => 37.467,
      north_bound => 37.477
    },
    'elevation_pillar_point_harbor_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -122.5,
      east_bound  => -122.470,
      south_bound => 37.448,
      north_bound => 37.509
    },
    'elevation_pismo_creek_lagoon_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -120.642,
      east_bound  => -120.633,
      south_bound => 35.130,
      north_bound => 35.140
    },
    'elevation_pistol_river_or_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -124.415,
      east_bound  => -124.397,
      south_bound => 42.268,
      north_bound => 42.292
    },
    'elevation_pomponio_creek_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -122.408,
      east_bound  => -122.404,
      south_bound => 37.296,
      north_bound => 37.300
    },
    'elevation_port_hueneme_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -119.219,
      east_bound  => -119.198,
      south_bound => 34.140,
      north_bound => 34.157
    },
    'elevation_port_orford_head_or_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -124.521,
      east_bound  => -124.497,
      south_bound => 42.742,
      north_bound => 42.762
    },
    'elevation_port_orford_or_mhw' => {
      begin_date  => '2008-03-14',
      end_date    => '2008-03-14',
      west_bound  => -124.97,
      east_bound  => -124.13,
      south_bound => 42.53,
      north_bound => 43.73
    },
    'elevation_port_san_luis_ca_mhw' => {
      begin_date  => '2011-05-20',
      end_date    => '2011-05-20',
      west_bound  => -121.3,
      east_bound  => -120.45,
      south_bound => 34.6,
      north_bound => 35.7
    },
    'elevation_port_townsend_wa_mhw' => {
      begin_date  => '2011-02-07',
      end_date    => '2011-02-07',
      west_bound  => -123.62,
      east_bound  => -122.2,
      south_bound => 47.91,
      north_bound => 48.79
    },
    'elevation_port_townsend_wa_navd88' => {
      begin_date  => '2011-02-07',
      end_date    => '2011-02-07',
      west_bound  => -123.62,
      east_bound  => -122.2,
      south_bound => 47.91,
      north_bound => 48.79
    },
    'elevation_pudding_creek_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -123.810,
      east_bound  => -123.781,
      south_bound => 39.443,
      north_bound => 39.465
    },
    'elevation_puget_sound_wa_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -124.114,
      east_bound  => -121.283,
      south_bound => 46.796,
      north_bound => 49.234
    },
    'elevation_puget_sound_wa_mhw' => {
      begin_date  => '2014-06-24',
      end_date    => '2014-06-24',
      west_bound  => -123.18,
      east_bound  => -122.16,
      south_bound => 47.01,
      north_bound => 48.19
    },
    'elevation_puget_sound_wa_navd88' => {
      begin_date  => '2014-06-24',
      end_date    => '2014-06-24',
      west_bound  => -123.18,
      east_bound  => -122.16,
      south_bound => 47.01,
      north_bound => 48.19
    },
    'elevation_queets_river_wa_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -124.359,
      east_bound  => -124.304,
      south_bound => 47.521,
      north_bound => 47.556
    },
    'elevation_quillayute_river_wa_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -124.652,
      east_bound  => -124.594,
      south_bound => 47.902,
      north_bound => 47.934
    },
    'elevation_quinault_river_wa_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -124.301,
      east_bound  => -124.159,
      south_bound => 47.320,
      north_bound => 47.377
    },
    'elevation_raft_river_wa_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -124.343,
      east_bound  => -124.318,
      south_bound => 47.450,
      north_bound => 47.468
    },
    'elevation_redwood_creek_big_lagoon_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -122.579,
      east_bound  => -122.572,
      south_bound => 37.859,
      north_bound => 37.863
    },
    'elevation_redwood_creek_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -124.095,
      east_bound  => -124.071,
      south_bound => 41.282,
      north_bound => 41.302
    },
    'elevation_rincon_creek_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -119.485,
      east_bound  => -119.437,
      south_bound => 34.373,
      north_bound => 34.407
    },
    'elevation_rock_creek_or_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -124.117,
      east_bound  => -124.011,
      south_bound => 44.171,
      north_bound => 44.207
    },
    'elevation_rodeo_lagoon_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -122.539,
      east_bound  => -122.523,
      south_bound => 37.828,
      north_bound => 37.835
    },
    'elevation_rogue_river_or_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -124.438,
      east_bound  => -124.385,
      south_bound => 42.412,
      north_bound => 42.448
    },
    'elevation_russian_gulch_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -123.806,
      east_bound  => -123.795,
      south_bound => 39.328,
      north_bound => 39.332
    },
    'elevation_russian_gulch_creek_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -123.159,
      east_bound  => -123.152,
      south_bound => 38.466,
      north_bound => 38.471
    },
    'elevation_russian_river_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -123.138,
      east_bound  => -123.051,
      south_bound => 38.418,
      north_bound => 38.468
    },
    'elevation_salmon_creek_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -123.069,
      east_bound  => -123.049,
      south_bound => 38.348,
      north_bound => 38.361
    },
    'elevation_salmon_river_estuary_or_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -124.021,
      east_bound  => -123.938,
      south_bound => 45.007,
      north_bound => 45.061
    },
    'elevation_samish_bay_wa_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -122.627,
      east_bound  => -122.338,
      south_bound => 48.473,
      north_bound => 48.698
    },
    'elevation_san_antonio_creek_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -120.623,
      east_bound  => -120.615,
      south_bound => 34.792,
      north_bound => 34.803
    },
    'elevation_san_carpoforo_creek_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -121.325,
      east_bound  => -121.317,
      south_bound => 35.762,
      north_bound => 35.769
    },
    'elevation_san_diego_bay_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -117.277,
      east_bound  => -117.057,
      south_bound => 32.573,
      north_bound => 32.768
    },
    'elevation_san_diego_ca_mhw' => {
      begin_date  => '2012-03-07',
      end_date    => '2012-03-07',
      west_bound  => -117.83,
      east_bound  => -117,
      south_bound => 32.45,
      north_bound => 33.6
    },
    'elevation_san_diego_ca_navd88' => {
      begin_date  => '2012-03-07',
      end_date    => '2012-03-07',
      west_bound  => -117.83,
      east_bound  => -117,
      south_bound => 32.45,
      north_bound => 33.6
    },
    'elevation_san_dieguito_river_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -117.275,
      east_bound  => -117.228,
      south_bound => 32.956,
      north_bound => 32.987
    },
    'elevation_san_elijo_lagoon_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -117.286,
      east_bound  => -117.243,
      south_bound => 32.996,
      north_bound => 33.023
    },
    'elevation_san_francisco_bay_1_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -122.571,
      east_bound  => -122.373,
      south_bound => 37.755,
      north_bound => 37.923
    },
    'elevation_san_francisco_bay_2_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -122.728,
      east_bound  => -121.876,
      south_bound => 37.373,
      north_bound => 38.100
    },
    'elevation_san_francisco_bay_3_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -122.713,
      east_bound  => -121.996,
      south_bound => 37.856,
      north_bound => 38.412
    },
    'elevation_san_francisco_bay_4_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -122.174,
      east_bound  => -121.804,
      south_bound => 37.975,
      north_bound => 38.224
    },
    'elevation_san_francisco_bay_5_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -122.187,
      east_bound  => -121.905,
      south_bound => 37.360,
      north_bound => 37.571
    },
    'elevation_san_francisco_ca_mhw' => {
      begin_date  => '2010-07-20',
      end_date    => '2010-07-20',
      west_bound  => -123.3,
      east_bound  => -121.85,
      south_bound => 37.32,
      north_bound => 38.48
    },
    'elevation_san_francisco_ca_navd88' => {
      begin_date  => '2010-07-20',
      end_date    => '2010-07-20',
      west_bound  => -123.3,
      east_bound  => -121.85,
      south_bound => 37.32,
      north_bound => 38.48
    },
    'elevation_san_gregorio_creek_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -122.405,
      east_bound  => -122.393,
      south_bound => 37.319,
      north_bound => 37.325
    },
    'elevation_san_juan_creek_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -117.688,
      east_bound  => -117.675,
      south_bound => 33.461,
      north_bound => 33.475
    },
    'elevation_san_luiz_obispo_bay_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -120.768,
      east_bound  => -120.672,
      south_bound => 35.140,
      north_bound => 35.204
    },
    'elevation_san_luiz_rey_river_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -117.396,
      east_bound  => -117.378,
      south_bound => 33.201,
      north_bound => 33.211
    },
    'elevation_san_mateo_lagoon_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -117.596,
      east_bound  => -117.590,
      south_bound => 33.385,
      north_bound => 33.390
    },
    'elevation_san_onofre_creek_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -117.580,
      east_bound  => -117.576,
      south_bound => 33.380,
      north_bound => 33.384
    },
    'elevation_san_pedro_bay_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -118.328,
      east_bound  => -118.004,
      south_bound => 33.657,
      north_bound => 33.887
    },
    'elevation_san_pedro_creek_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -122.506,
      east_bound  => -122.504,
      south_bound => 37.595,
      north_bound => 37.597
    },
    'elevation_san_simeon_creek_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -121.129,
      east_bound  => -121.121,
      south_bound => 35.593,
      north_bound => 35.598
    },
    'elevation_sand_lake_estuary_or_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -123.978,
      east_bound  => -123.910,
      south_bound => 45.254,
      north_bound => 45.308
    },
    'elevation_santa_ana_river_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -117.976,
      east_bound  => -117.921,
      south_bound => 33.621,
      north_bound => 33.690
    },
    'elevation_santa_barbara_ca_mhw' => {
      begin_date  => '2008-12-15',
      end_date    => '2008-12-15',
      west_bound  => -120.51,
      east_bound  => -119.14,
      south_bound => 33.77,
      north_bound => 34.62
    },
    'elevation_santa_clara_river_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -119.268,
      east_bound  => -119.237,
      south_bound => 34.216,
      north_bound => 34.243
    },
    'elevation_santa_margarita_river_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -117.422,
      east_bound  => -117.378,
      south_bound => 33.218,
      north_bound => 33.248
    },
    'elevation_santa_maria_river_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -120.653,
      east_bound  => -120.635,
      south_bound => 34.958,
      north_bound => 34.978
    },
    'elevation_santa_monica_bay_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -118.482,
      east_bound  => -118.407,
      south_bound => 33.948,
      north_bound => 34.002
    },
    'elevation_santa_monica_ca_mhw' => {
      begin_date  => '2010-03-12',
      end_date    => '2010-03-12',
      west_bound  => -119.14,
      east_bound  => -117.8,
      south_bound => 33.2,
      north_bound => 34.2
    },
    'elevation_santa_monica_ca_navd88' => {
      begin_date  => '2010-03-12',
      end_date    => '2010-03-12',
      west_bound  => -119.14,
      east_bound  => -117.8,
      south_bound => 33.2,
      north_bound => 34.2
    },
    'elevation_santa_monica_canyon_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -118.521,
      east_bound  => -118.518,
      south_bound => 34.025,
      north_bound => 34.028
    },
    'elevation_santa_rosa_creek_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -121.111,
      east_bound  => -121.103,
      south_bound => 35.566,
      north_bound => 35.571
    },
    'elevation_santa_ynez_river_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -120.605,
      east_bound  => -120.578,
      south_bound => 34.681,
      north_bound => 34.698
    },
    'elevation_schwan_lagoon_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -122.000,
      east_bound  => -121.988,
      south_bound => 36.961,
      north_bound => 36.970
    },
    'elevation_scott_creek_lagoon_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -122.232,
      east_bound  => -122.224,
      south_bound => 37.039,
      north_bound => 37.045
    },
    'elevation_seaside_creek_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -123.768,
      east_bound  => -123.762,
      south_bound => 39.556,
      north_bound => 39.561
    },
    'elevation_seaside_or_mhw' => {
      begin_date  => '2004-04-01',
      end_date    => '2004-04-01',
      west_bound  => -124.04,
      east_bound  => -123.89,
      south_bound => 45.9,
      north_bound => 46.08
    },
   'elevation_siletz_bay_or_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -124.085,
      east_bound  => -123.864,
      south_bound => 44.779,
      north_bound => 44.958
    },
    'elevation_siltcoos_river_or_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -124.160,
      east_bound  => -124.130,
      south_bound => 43.870,
      north_bound => 43.885
    },
    'elevation_sixes_river_or_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -124.583,
      east_bound  => -124.509,
      south_bound => 42.806,
      north_bound => 42.864
    },
    'elevation_smith_lake_or_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -123.951,
      east_bound  => -123.947,
      south_bound => 45.587,
      north_bound => 45.592
    },
    'elevation_smith_river_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -124.231,
      east_bound  => -124.140,
      south_bound => 41.892,
      north_bound => 41.967
    },
    'elevation_squashan_creek_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -124.072,
      east_bound  => -124.067,
      south_bound => 41.388,
      north_bound => 41.392
    },
    'elevation_squaw_creek_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -124.115,
      east_bound  => -124.069,
      south_bound => 44.202,
      north_bound => 44.209
    },
    'elevation_steamboat_creek_wa_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -124.407,
      east_bound  => -124.339,
      south_bound => 47.669,
      north_bound => 47.697
    },
    'elevation_stone_lagoon_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -124.109,
      east_bound  => -124.081,
      south_bound => 41.227,
      north_bound => 41.260
    },
    'elevation_strait_of_juan_de_fuca_wa_mllw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -124.912,
      east_bound  => -122.948,
      south_bound => 47.803,
      north_bound => 48.811
    },
    'elevation_stump_beach_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -123.337,
      east_bound  => -123.331,
      south_bound => 38.577,
      north_bound => 38.582
    },
    'elevation_suislaw_river_or_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -124.169,
      east_bound  => -123.821,
      south_bound => 43.904,
      north_bound => 44.102
    },
    'elevation_sunset_bay_or_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -124.376,
      east_bound  => -124.368,
      south_bound => 43.330,
      north_bound => 43.334
    },
    'elevation_sutton_creek_or_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -124.135,
      east_bound  => -124.124,
      south_bound => 44.049,
      north_bound => 44.059
    },
    'elevation_sycamore_creek_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -119.668,
      east_bound  => -119.666,
      south_bound => 34.416,
      north_bound => 34.419
    },
    'elevation_tahkenitch_creek_or_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -124.183,
      east_bound  => -124.147,
      south_bound => 43.790,
      north_bound => 43.819
    },
    'elevation_tahola_wa_mhw' => {
      begin_date  => '2009-07-24',
      end_date    => '2009-07-24',
      west_bound  => -125.3,
      east_bound  => -123.7,
      south_bound => 46.99,
      north_bound => 47.65
    },
    'elevation_tajiguas_creek_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -120.112,
      east_bound  => -120.083,
      south_bound => 34.464,
      north_bound => 34.500
    },
    'elevation_ten_mile_river_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -123.770,
      east_bound  => -123.733,
      south_bound => 39.523,
      north_bound => 39.560
    },
    'elevation_tenmile_creek_1_or_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -124.268,
      east_bound  => -124.138,
      south_bound => 43.547,
      north_bound => 43.663
    },
    'elevation_tenmile_creek_2_or_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -124.113,
      east_bound  => -124.106,
      south_bound => 44.222,
      north_bound => 44.228
    },
    'elevation_tennessee_valley_lagoon_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -122.553,
      east_bound  => -122.550,
      south_bound => 37.841,
      north_bound => 37.844
    },
    'elevation_thomas_creek_or_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -124.362,
      east_bound  => -124.359,
      south_bound => 42.165,
      north_bound => 42.166
    },
    'elevation_threemile_creek_or_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -124.183,
      east_bound  => -124.157,
      south_bound => 43.743,
      north_bound => 43.754
    },
    'elevation_tijuana_river_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -117.142,
      east_bound  => -117.092,
      south_bound => 32.533,
      north_bound => 32.582
    },
    'elevation_tillamook_bay_or_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -124.016,
      east_bound  => -123.733,
      south_bound => 45.386,
      north_bound => 45.613
    },
    'elevation_tomales_bay_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -123.036,
      east_bound  => -122.769,
      south_bound => 38.049,
      north_bound => 38.272
    },
    'elevation_trancas_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -118.843,
      east_bound  => -118.840,
      south_bound => 34.028,
      north_bound => 34.031
    },
    'elevation_tunitas_creek_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -122.402,
      east_bound  => -122.398,
      south_bound => 37.356,
      north_bound => 37.360
    },
    'elevation_twomile_creek_2_or_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -124.392,
      east_bound  => -124.368,
      south_bound => 43.219,
      north_bound => 43.240
    },
    'elevation_umpqua_river_or_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -124.262,
      east_bound  => -123.802,
      south_bound => 43.580,
      north_bound => 43.864
    },
    'elevation_usal_creek_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -123.853,
      east_bound  => -123.845,
      south_bound => 39.829,
      north_bound => 39.835
    },
    'elevation_ventura_harbor_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -119.279,
      east_bound  => -119.253,
      south_bound => 34.238,
      north_bound => 34.264
    },
    'elevation_ventura_river_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -119.312,
      east_bound  => -119.305,
      south_bound => 34.273,
      north_bound => 34.281
    },
    'elevation_villa_creek_lagoon_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -120.974,
      east_bound  => -120.969,
      south_bound => 35.459,
      north_bound => 35.464
    },
    'elevation_virgin_creek_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -123.806,
      east_bound  => -123.798,
      south_bound => 39.469,
      north_bound => 39.475
    },
    'elevation_waddell_creek_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -122.281,
      east_bound  => -122.272,
      south_bound => 37.093,
      north_bound => 37.101
    },
    'elevation_wages_creek_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -123.786,
      east_bound  => -123.782,
      south_bound => 39.649,
      north_bound => 39.652
    },
    'elevation_whale_creek_wa_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -124.346,
      east_bound  => -124.341,
      south_bound => 47.487,
      north_bound => 47.491
    },
    'elevation_whaleshead_creek_or_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -124.358,
      east_bound  => -124.353,
      south_bound => 42.142,
      north_bound => 42.145
    },
    'elevation_whisky_run_or_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -124.422,
      east_bound  => -124.371,
      south_bound => 43.183,
      north_bound => 43.239
    },
    'elevation_wilder_creek_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -122.083,
      east_bound  => -122.075,
      south_bound => 36.952,
      north_bound => 36.958
    },
    'elevation_willapa_bay_wa_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -124.356,
      east_bound  => -123.455,
      south_bound => 46.274,
      north_bound => 46.995
    },
    'elevation_wilson_creek_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -124.102,
      east_bound  => -124.100,
      south_bound => 41.603,
      north_bound => 41.605
    },
    'elevation_winchuck_river_or_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -124.217,
      east_bound  => -124.193,
      south_bound => 41.993,
      north_bound => 42.010
    },
    'elevation_wreck_creek_wa_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -124.238,
      east_bound  => -124.231,
      south_bound => 47.283,
      north_bound => 47.287
    },
    'elevation_yachats_river_or_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -124.109,
      east_bound  => -124.073,
      south_bound => 44.299,
      north_bound => 44.316
    },
    'elevation_yaquina_bay_or_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -124.112,
      east_bound  => -123.854,
      south_bound => 44.525,
      north_bound => 44.674
    },
    'elevation_younger_lagoon_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -122.069,
      east_bound  => -122.065,
      south_bound => 36.949,
      north_bound => 36.954
    },
    'elevation_zuma_canyon_ca_mhw' => {
      begin_date  => '2020-06-18',
      end_date    => '2020-06-18',
      west_bound  => -118.822,
      east_bound  => -118.819,
      south_bound => 34.013,
      north_bound => 34.016
    },

    'nadp.gridded.wet_nh4' => {
      begin_date  => '1985-01-01',
      end_date    => '2014-12-31',
      west_bound  => -125,
      east_bound  => -109,
      south_bound => 32,
      north_bound => 49
    },
    'nadp.gridded.wet_no3' => {
      begin_date  => '1985-01-01',
      end_date    => '2014-12-31',
      west_bound  => -125,
      east_bound  => -109,
      south_bound => 32,
      north_bound => 49
    },
    'nadp.gridded.wet_totn' => {
      begin_date  => '1985-01-01',
      end_date    => '2014-12-31',
      west_bound  => -125,
      east_bound  => -109,
      south_bound => 32,
      north_bound => 49
    },
    'nadpcmaq.gridded.dry_nh3' => {
      begin_date  => '2000-01-01',
      end_date    => '2014-12-31',
      west_bound  => -125,
      east_bound  => -109,
      south_bound => 32,
      north_bound => 49
    },
    'nadpcmaq.gridded.dry_nh4' => {
      begin_date  => '2000-01-01',
      end_date    => '2014-12-31',
      west_bound  => -125,
      east_bound  => -109,
      south_bound => 32,
      north_bound => 49
    },
    'nadpcmaq.gridded.dry_no3' => {
      begin_date  => '2000-01-01',
      end_date    => '2014-12-31',
      west_bound  => -125,
      east_bound  => -109,
      south_bound => 32,
      north_bound => 49
    },
    'nadpcmaq.gridded.dry_totn' => {
      begin_date  => '2000-01-01',
      end_date    => '2014-12-31',
      west_bound  => -125,
      east_bound  => -109,
      south_bound => 32,
      north_bound => 49
    },
    'nadpcmaq.gridded.dry_totnp' => {
      begin_date  => '2000-01-01',
      end_date    => '2014-12-31',
      west_bound  => -125,
      east_bound  => -109,
      south_bound => 32,
      north_bound => 49
    },
    'nadpcmaq.gridded.dry_totun' => {
      begin_date  => '2000-01-01',
      end_date    => '2014-12-31',
      west_bound  => -125,
      east_bound  => -109,
      south_bound => 32,
      north_bound => 49
    },
    'nadpcmaq.gridded.dry_totunp' => {
      begin_date  => '2000-01-01',
      end_date    => '2014-12-31',
      west_bound  => -125,
      east_bound  => -109,
      south_bound => 32,
      north_bound => 49
    },
    'nadpcmaq.gridded.dry_net_nh3' => {
      begin_date  => '2000-01-01',
      end_date    => '2014-12-31',
      west_bound  => -125,
      east_bound  => -109,
      south_bound => 32,
      north_bound => 49
    },
    'nadpcmaq.gridded.dry_oxidized_totn' => {
      begin_date  => '2000-01-01',
      end_date    => '2014-12-31',
      west_bound  => -125,
      east_bound  => -109,
      south_bound => 32,
      north_bound => 49
    },
    'nadpcmaq.gridded.dry_oxidized_totnp' => {
      begin_date  => '2000-01-01',
      end_date    => '2014-12-31',
      west_bound  => -125,
      east_bound  => -109,
      south_bound => 32,
      north_bound => 49
    },
    'nadpcmaq.gridded.dry_reduced_totn' => {
      begin_date  => '2000-01-01',
      end_date    => '2014-12-31',
      west_bound  => -125,
      east_bound  => -109,
      south_bound => 32,
      north_bound => 49
    },
    'nadpcmaq.gridded.dry_reduced_totnp' => {
      begin_date  => '2000-01-01',
      end_date    => '2014-12-31',
      west_bound  => -125,
      east_bound  => -109,
      south_bound => 32,
      north_bound => 49
    },
    'nadpcmaq.gridded.wet_nh4' => {
      begin_date  => '2000-01-01',
      end_date    => '2014-12-31',
      west_bound  => -125,
      east_bound  => -109,
      south_bound => 32,
      north_bound => 49
    },
    'nadpcmaq.gridded.wet_no3' => {
      begin_date  => '2000-01-01',
      end_date    => '2014-12-31',
      west_bound  => -125,
      east_bound  => -109,
      south_bound => 32,
      north_bound => 49
    },
    'nadpcmaq.gridded.wet_totn' => {
      begin_date  => '2000-01-01',
      end_date    => '2014-12-31',
      west_bound  => -125,
      east_bound  => -109,
      south_bound => 32,
      north_bound => 49
    },
    'nadpcmaq.gridded.oxidized_totn' => {
      begin_date  => '2000-01-01',
      end_date    => '2014-12-31',
      west_bound  => -125,
      east_bound  => -109,
      south_bound => 32,
      north_bound => 49
    },
    'nadpcmaq.gridded.oxidized_totnp' => {
      begin_date  => '2000-01-01',
      end_date    => '2014-12-31',
      west_bound  => -125,
      east_bound  => -109,
      south_bound => 32,
      north_bound => 49
    },
    'nadpcmaq.gridded.reduced_totn' => {
      begin_date  => '2000-01-01',
      end_date    => '2014-12-31',
      west_bound  => -125,
      east_bound  => -109,
      south_bound => 32,
      north_bound => 49
    },
    'nadpcmaq.gridded.reduced_totnp' => {
      begin_date  => '2000-01-01',
      end_date    => '2014-12-31',
      west_bound  => -125,
      east_bound  => -109,
      south_bound => 32,
      north_bound => 49
    },
    'nadpcmaq.gridded.total_totn' => {
      begin_date  => '2000-01-01',
      end_date    => '2014-12-31',
      west_bound  => -125,
      east_bound  => -109,
      south_bound => 32,
      north_bound => 49
    },

    'temperature_minimum' => {
      begin_date  => '1971-01-01',
      end_date    => '2000-12-31',
      west_bound  => -125,
      east_bound  => -109,
      south_bound => 32,
      north_bound => 49
    },
    'temperature_maximum' => {
      begin_date  => '1971-01-01',
      end_date    => '2000-12-31',
      west_bound  => -125,
      east_bound  => -109,
      south_bound => 32,
      north_bound => 49
    },
    'precipitation_mean' => {
      begin_date  => '1971-01-01',
      end_date    => '2000-12-31',
      west_bound  => -125,
      east_bound  => -109,
      south_bound => 32,
      north_bound => 49
    },
    'monthly_temperature_minimum' => {
      begin_date  => '1971-01-01',
      end_date    => '2000-12-31',
      west_bound  => -125,
      east_bound  => -109,
      south_bound => 32,
      north_bound => 49
    },
    'monthly_temperature_maximum' => {
      begin_date  => '1971-01-01',
      end_date    => '2000-12-31',
      west_bound  => -125,
      east_bound  => -109,
      south_bound => 32,
      north_bound => 49
    },
    'monthly_precipitation_mean' => {
      begin_date  => '1971-01-01',
      end_date    => '2000-12-31',
      west_bound  => -125,
      east_bound  => -109,
      south_bound => 32,
      north_bound => 49
    },
    'annual_temperature_minimum' => {
      begin_date  => '1971-01-01',
      end_date    => '2000-12-31',
      west_bound  => -125,
      east_bound  => -109,
      south_bound => 32,
      north_bound => 49
    },
    'annual_temperature_maximum' => {
      begin_date  => '1971-01-01',
      end_date    => '2000-12-31',
      west_bound  => -125,
      east_bound  => -109,
      south_bound => 32,
      north_bound => 49
    },
    'annual_precipitation_mean' => {
      begin_date  => '1971-01-01',
      end_date    => '2000-12-31',
      west_bound  => -125,
      east_bound  => -109,
      south_bound => 32,
      north_bound => 49
    },
    'nitrogen.source_gridded.crops' => {
      begin_date  => '2006-01-01',
      end_date    => '2006-12-31',
      west_bound  => -125,
      east_bound  => -109,
      south_bound => 32,
      north_bound => 49
    },
    'nitrogen.source_gridded.fertilizer' => {
      begin_date  => '2006-01-01',
      end_date    => '2006-12-31',
      west_bound  => -125,
      east_bound  => -109,
      south_bound => 32,
      north_bound => 49
    },
    'nitrogen.source_gridded.manure' => {
      begin_date  => '2006-01-01',
      end_date    => '2006-12-31',
      west_bound  => -125,
      east_bound  => -109,
      south_bound => 32,
      north_bound => 49
    },
    'seagrass' => {
      begin_date  => '1979-01-01',
      end_date    => '2009-12-31',
      west_bound  => -125,
      east_bound  => -109,
      south_bound => 32,
      north_bound => 49
    },
    'wetlands' => {
      begin_date  => '1981-01-01',
      end_date    => '2011-12-31',
      west_bound  => -125,
      east_bound  => -109,
      south_bound => 32,
      north_bound => 49
    },
    'population_iclus1' => {
      begin_date  => '2010-01-01',
      end_date    => '2090-12-31',
      west_bound  => -125,
      east_bound  => -109,
      south_bound => 32,
      north_bound => 49
    },
    'population_iclus2' => {
      begin_date  => '2020-01-01',
      end_date    => '2100-12-31',
      west_bound  => -125,
      east_bound  => -109,
      south_bound => 32,
      north_bound => 49
    },
    'land_use_iclus_watershed' => {
      begin_date  => '2000-01-01',
      end_date    => '2100-12-31',
      west_bound  => -125,
      east_bound  => -109,
      south_bound => 32,
      north_bound => 49
    },
    'land_use_nlcd_2011_watershed' => {
      begin_date  => '2011-01-01',
      end_date    => '2011-12-31',
      west_bound  => -125,
      east_bound  => -109,
      south_bound => 32,
      north_bound => 49
    },
    'land_use_nlcd_2006_watershed' => {
      begin_date  => '2006-01-01',
      end_date    => '2006-12-31',
      west_bound  => -125,
      east_bound  => -109,
      south_bound => 32,
      north_bound => 49
    },
    'soil' => {
      begin_date  => '2012-01-01',
      end_date    => '2012-12-31',
      west_bound  => -125,
      east_bound  => -109,
      south_bound => 32,
      north_bound => 49
    },
    'nitrogen.estuary_cmaq' => { # Old
      begin_date  => '2002-01-01',
      end_date    => '2008-12-31',
      west_bound  => -125,
      east_bound  => -109,
      south_bound => 32,
      north_bound => 49
    },
    'dry_nh3.estuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -125,
      east_bound  => -123,
      south_bound => 42,
      north_bound => 49
    },
    'dry_natural_nh3.estuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -125,
      east_bound  => -123,
      south_bound => 42,
      north_bound => 49
    },
    'dry_nhx.estuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -125,
      east_bound  => -123,
      south_bound => 42,
      north_bound => 49
    },
    'dry_so2.estuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -125,
      east_bound  => -123,
      south_bound => 42,
      north_bound => 49
    },
    'dry_oxidized_nox.estuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -125,
      east_bound  => -123,
      south_bound => 42,
      north_bound => 49
    },
    'dry_oxidized_organic_nitrogen.estuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -125,
      east_bound  => -123,
      south_bound => 42,
      north_bound => 49
    },
    'dry_oxidized_no3.estuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -125,
      east_bound  => -123,
      south_bound => 42,
      north_bound => 49
    },
    'dry_oxidized_nitrogen.estuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -125,
      east_bound  => -123,
      south_bound => 42,
      north_bound => 49
    },
    'dry_reduced_nitrogen.estuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -125,
      east_bound  => -123,
      south_bound => 42,
      north_bound => 49
    },
    'dry_sulfur.estuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -125,
      east_bound  => -123,
      south_bound => 42,
      north_bound => 49
    },
    'chloride.estuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -125,
      east_bound  => -123,
      south_bound => 42,
      north_bound => 49
    },
    'total_nitrogen.estuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -125,
      east_bound  => -123,
      south_bound => 42,
      north_bound => 49
    },
    'oxidized_nitrogen.estuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -125,
      east_bound  => -123,
      south_bound => 42,
      north_bound => 49
    },
    'reduced_nitrogen.estuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -125,
      east_bound  => -123,
      south_bound => 42,
      north_bound => 49
    },
    'sulfur.estuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -125,
      east_bound  => -123,
      south_bound => 42,
      north_bound => 49
    },
    'wet_nh3.estuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -125,
      east_bound  => -123,
      south_bound => 42,
      north_bound => 49
    },
    'wet_chloride.estuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -125,
      east_bound  => -123,
      south_bound => 42,
      north_bound => 49
    },
    'wet_so4.estuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -125,
      east_bound  => -123,
      south_bound => 42,
      north_bound => 49
    },
    'wet_oxidized_no3.estuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -125,
      east_bound  => -123,
      south_bound => 42,
      north_bound => 49
    },
    'wet_oxidized_nox.estuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -125,
      east_bound  => -123,
      south_bound => 42,
      north_bound => 49
    },
    'wet_organic_nitrogen.estuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -125,
      east_bound  => -123,
      south_bound => 42,
      north_bound => 49
    },
    'wet_oxidized_nitrogen.estuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -125,
      east_bound  => -123,
      south_bound => 42,
      north_bound => 49
    },
    'wet_reduced_nitrogen.estuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -125,
      east_bound  => -123,
      south_bound => 42,
      north_bound => 49
    },
    'wet_sulfur.estuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -125,
      east_bound  => -123,
      south_bound => 42,
      north_bound => 49
    },
    'nitrogen.subestuary_cmaq' => { # Old
      begin_date  => '2002-01-01',
      end_date    => '2008-12-31',
      west_bound  => -125,
      east_bound  => -123,
      south_bound => 42,
      north_bound => 49
    },
    'dry_nh3.subestuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -125,
      east_bound  => -123,
      south_bound => 42,
      north_bound => 49
    },
    'dry_natural_nh3.subestuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -125,
      east_bound  => -123,
      south_bound => 42,
      north_bound => 49
    },
    'dry_nhx.subestuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -125,
      east_bound  => -123,
      south_bound => 42,
      north_bound => 49
    },
    'dry_so2.subestuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -125,
      east_bound  => -123,
      south_bound => 42,
      north_bound => 49
    },
    'dry_oxidized_nox.subestuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -125,
      east_bound  => -123,
      south_bound => 42,
      north_bound => 49
    },
    'dry_oxidized_organic_nitrogen.subestuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -125,
      east_bound  => -123,
      south_bound => 42,
      north_bound => 49
    },
    'dry_oxidized_no3.subestuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -125,
      east_bound  => -123,
      south_bound => 42,
      north_bound => 49
    },
    'dry_oxidized_nitrogen.subestuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -125,
      east_bound  => -123,
      south_bound => 42,
      north_bound => 49
    },
    'dry_reduced_nitrogen.subestuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -125,
      east_bound  => -123,
      south_bound => 42,
      north_bound => 49
    },
    'dry_sulfur.subestuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -125,
      east_bound  => -123,
      south_bound => 42,
      north_bound => 49
    },
    'chloride.subestuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -125,
      east_bound  => -123,
      south_bound => 42,
      north_bound => 49
    },
    'total_nitrogen.subestuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -125,
      east_bound  => -123,
      south_bound => 42,
      north_bound => 49
    },
    'oxidized_nitrogen.subestuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -125,
      east_bound  => -123,
      south_bound => 42,
      north_bound => 49
    },
    'reduced_nitrogen.subestuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -125,
      east_bound  => -123,
      south_bound => 42,
      north_bound => 49
    },
    'sulfur.subestuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -125,
      east_bound  => -123,
      south_bound => 42,
      north_bound => 49
    },
    'wet_nh3.subestuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -125,
      east_bound  => -123,
      south_bound => 42,
      north_bound => 49
    },
    'wet_chloride.subestuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -125,
      east_bound  => -123,
      south_bound => 42,
      north_bound => 49
    },
    'wet_so4.subestuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -125,
      east_bound  => -123,
      south_bound => 42,
      north_bound => 49
    },
    'wet_oxidized_no3.subestuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -125,
      east_bound  => -123,
      south_bound => 42,
      north_bound => 49
    },
    'wet_oxidized_nox.subestuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -125,
      east_bound  => -123,
      south_bound => 42,
      north_bound => 49
    },
    'wet_organic_nitrogen.subestuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -125,
      east_bound  => -123,
      south_bound => 42,
      north_bound => 49
    },
    'wet_oxidized_nitrogen.subestuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -125,
      east_bound  => -123,
      south_bound => 42,
      north_bound => 49
    },
    'wet_reduced_nitrogen.subestuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -125,
      east_bound  => -123,
      south_bound => 42,
      north_bound => 49
    },
    'wet_sulfur.subestuary_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2012-12-31',
      west_bound  => -125,
      east_bound  => -123,
      south_bound => 42,
      north_bound => 49
    },
    'nitrogen.load_estuary_sparrow_1992' => {
      begin_date  => '1992-01-01',
      end_date    => '1992-12-31',
      west_bound  => -125,
      east_bound  => -109,
      south_bound => 32,
      north_bound => 49
    },
    'nitrogen.load_estuary_sparrow_2002_mrb7' => {
      begin_date  => '1992-01-01',
      end_date    => '1992-12-31',
      west_bound  => -125,
      east_bound  => -123,
      south_bound => 42,
      north_bound => 49
    },
    'nitrogen.load_estuary_sparrow_2011_mrb7' => {
      begin_date  => '2011-01-01',
      end_date    => '2011-12-31',
      west_bound  => -125,
      east_bound  => -121,
      south_bound => 42,
      north_bound => 49
    },
    'nitrogen.load_estuary_sparrow_2002_mrb8' => {
      begin_date  => '2011-01-01',
      end_date    => '2011-12-31',
      west_bound  => -125,
      east_bound  => -117,
      south_bound => 32,
      north_bound => 42
    },
    'nitrogen.load_subestuary_sparrow_1992' => {
      begin_date  => '1992-01-01',
      end_date    => '1992-12-31',
      west_bound  => -125,
      east_bound  => -109,
      south_bound => 32,
      north_bound => 49
    },
    'nitrogen.load_subestuary_sparrow_2002_mrb7' => {
      begin_date  => '1992-01-01',
      end_date    => '1992-12-31',
      west_bound  => -125,
      east_bound  => -123,
      south_bound => 42,
      north_bound => 49
    },
    'total_nitrogen.source_watershed_nadp' => {
      begin_date  => '2000-01-01',
      end_date    => '2014-12-31',
      west_bound  => -125,
      east_bound  => -109,
      south_bound => 32,
      north_bound => 49
    },
    'wet_nitrogen.source_watershed_nadp' => {
      begin_date  => '2000-01-01',
      end_date    => '2014-12-31',
      west_bound  => -125,
      east_bound  => -109,
      south_bound => 32,
      north_bound => 49
    },
    'dry_nitrogen.source_watershed_nadp' => {
      begin_date  => '2000-01-01',
      end_date    => '2014-12-31',
      west_bound  => -125,
      east_bound  => -109,
      south_bound => 32,
      north_bound => 49
    },
    'oxidized_nitrogen.source_watershed_nadp' => {
      begin_date  => '2000-01-01',
      end_date    => '2014-12-31',
      west_bound  => -125,
      east_bound  => -109,
      south_bound => 32,
      north_bound => 49
    },
    'reduced_nitrogen.source_watershed_nadp' => {
      begin_date  => '2000-01-01',
      end_date    => '2014-12-31',
      west_bound  => -125,
      east_bound  => -109,
      south_bound => 32,
      north_bound => 49
    },
    'dry_oxidized_nitrogen.source_watershed_nadp' => {
      begin_date  => '2000-01-01',
      end_date    => '2014-12-31',
      west_bound  => -125,
      east_bound  => -109,
      south_bound => 32,
      north_bound => 49
    },
    'dry_reduced_nitrogen.source_watershed_nadp' => {
      begin_date  => '2000-01-01',
      end_date    => '2014-12-31',
      west_bound  => -125,
      east_bound  => -109,
      south_bound => 32,
      north_bound => 49
    },
    'wet_no3.source_watershed_nadp' => {
      begin_date  => '1994-01-01',
      end_date    => '2014-12-31',
      west_bound  => -125,
      east_bound  => -109,
      south_bound => 32,
      north_bound => 49
    },
    'wet_nh4.source_watershed_nadp' => {
      begin_date  => '1994-01-01',
      end_date    => '2014-12-31',
      west_bound  => -125,
      east_bound  => -109,
      south_bound => 32,
      north_bound => 49
    },
    'dry_no3.source_watershed_nadp' => {
      begin_date  => '2000-01-01',
      end_date    => '2014-12-31',
      west_bound  => -125,
      east_bound  => -109,
      south_bound => 32,
      north_bound => 49
    },
    'dry_nh4.source_watershed_nadp' => {
      begin_date  => '2000-01-01',
      end_date    => '2014-12-31',
      west_bound  => -125,
      east_bound  => -109,
      south_bound => 32,
      north_bound => 49
    },
    'nitrogen.source_watershed_cmaq' => {
      begin_date  => '2002-01-01',
      end_date    => '2008-12-31',
      west_bound  => -125,
      east_bound  => -109,
      south_bound => 32,
      north_bound => 49
    },
    'nitrogen.source_watershed_point' => {
      begin_date  => '2007-01-01',
      end_date    => '2011-12-31',
      west_bound  => -125,
      east_bound  => -109,
      south_bound => 32,
      north_bound => 49
    },
    'nitrogen.source_watershed_nonpoint' => {
      begin_date  => '2006-01-01',
      end_date    => '2006-12-31',
      west_bound  => -125,
      east_bound  => -109,
      south_bound => 32,
      north_bound => 49
    },
    'phosphorus.source_watershed_point' => {
      begin_date  => '2007-01-01',
      end_date    => '2016-12-31',
      west_bound  => -125,
      east_bound  => -109,
      south_bound => 32,
      north_bound => 49
    },
    'sediment_nca' => {
      begin_date  => '1990-01-01',
      end_date    => '2006-12-31',
      west_bound  => -125,
      east_bound  => -109,
      south_bound => 32,
      north_bound => 49
    },
    'sediment_nca_1990-2006' => {
      begin_date  => '1990-01-01',
      end_date    => '2006-12-31',
      west_bound  => -125,
      east_bound  => -109,
      south_bound => 32,
      north_bound => 49
    },
    'sediment_nca_2015' => {
      begin_date  => '2015-01-01',
      end_date    => '2015-12-31',
      west_bound  => -125,
      east_bound  => -109,
      south_bound => 32,
      north_bound => 49
    },
    'sediment_nca' => {
      begin_date  => '1990-01-01',
      end_date    => '2006-12-31',
      west_bound  => -125,
      east_bound  => -109,
      south_bound => 32,
      north_bound => 49
    },
    'sediment_calculated' => {
      begin_date  => '1990-01-01',
      end_date    => '2006-12-31',
      west_bound  => -125,
      east_bound  => -109,
      south_bound => 32,
      north_bound => 49
    },
    'sediment_extracted' => {
      begin_date  => '1990-01-01',
      end_date    => '2006-12-31',
      west_bound  => -125,
      east_bound  => -109,
      south_bound => 32,
      north_bound => 49
    },
    'sediment_parsed' => {
      begin_date  => '1990-01-01',
      end_date    => '2006-12-31',
      west_bound  => -125,
      east_bound  => -109,
      south_bound => 32,
      north_bound => 49
    },
    'coastal_vulnerability' => {
      begin_date  => '1976-01-01',
      end_date    => '1995-12-31',
      west_bound  => -82,
      east_bound  => -67,
      south_bound => 24,
      north_bound => 46
    },
    'sensitivity_estuary' => {
      begin_date  => '2016-10-05',
      end_date    => '2016-10-05',
      west_bound  => -82,
      east_bound  => -67,
      south_bound => 24,
      north_bound => 46
    },
    'groundwater_contact_catchment' => {
      begin_date  => '2002-01-01',
      end_date    => '2002-12-31',
      west_bound  => -82,
      east_bound  => -67,
      south_bound => 24,
      north_bound => 46
    },
    'groundwater_contact_watershed' => {
      begin_date  => '2002-01-01',
      end_date    => '2002-12-31',
      west_bound  => -82,
      east_bound  => -67,
      south_bound => 24,
      north_bound => 46
    },
    'stream_discharge_nhd_line' => {
      begin_date  => '1952-01-01',
      end_date    => '2017-12-31',
      west_bound  => -125,
      east_bound  => -117,
      south_bound => 33,
      north_bound => 49
    },
    'stream_discharge_usgs_point' => {
      begin_date  => '1952-01-01',
      end_date    => '2017-12-31',
      west_bound  => -125,
      east_bound  => -117,
      south_bound => 33,
      north_bound => 49
    },
    'tide_point' => {
      begin_date  => '1952-01-01',
      end_date    => '2017-12-31',
      west_bound  => -125,
      east_bound  => -117,
      south_bound => 33,
      north_bound => 49
    },
    'tide_current' => {
      begin_date  => '2017-01-01',
      end_date    => '2019-12-31',
      west_bound  => -125,
      east_bound  => -117,
      south_bound => 32,
      north_bound => 49
    },
    'longshore_current' => {
      begin_date  => '2012-01-01',
      end_date    => '2020-12-31',
      west_bound  => -125,
      east_bound  => -117,
      south_bound => 32,
      north_bound => 49
    },
    'estuary_flushing' => {
      begin_date  => '1950-01-01',
      end_date    => '2020-12-31',
      west_bound  => -125,
      east_bound  => -117,
      south_bound => 32,
      north_bound => 49
    },
    'salinity_point' => {
      begin_date  => '1990-01-01',
      end_date    => '2018-12-31',
      west_bound  => -125,
      east_bound  => -117,
      south_bound => 32,
      north_bound => 49
    },
    'seagrass_point' => {
      begin_date  => '1978-01-01',
      end_date    => '2017-12-31',
      west_bound  =>  -124.71,
      east_bound  =>  -117.12,
      south_bound =>   32.55,
      north_bound =>   48.98
    },
    'seagrass_polygon' => {
      begin_date  => '1978-01-01',
      end_date    => '2017-12-31',
      west_bound  =>  -124.51,
      east_bound  =>  -117.98,
      south_bound =>   32.59,
      north_bound =>   48.99
    },
    'esat1_water_quality' => {
      begin_date  => '1984-01-01',
      end_date    => '2021-12-31',
      west_bound  =>  -124.02,
      east_bound  =>  -70.52,
      south_bound =>   24.91,
      north_bound =>   48.56
    },
    'esat2_water_quality' => {
      begin_date  => '2018-01-01',
      end_date    => '2020-12-31',
      west_bound  =>  -124.02,
      east_bound  =>  -70.52,
      south_bound =>   24.91,
      north_bound =>   48.56
    },
    'stream_temperature_point_08_lower_columbia_river' => {
      begin_date  => '2000-01-01',
      end_date    => '2080-12-31',
      west_bound  => -125.0,
      east_bound  => -116.0,
      south_bound => 42.0,
      north_bound => 46.0
    },
    'stream_temperature_point_08_middle_columbia_river' => {
      begin_date  => '2000-01-01',
      end_date    => '2080-12-31',
      west_bound  => -125.0,
      east_bound  => -116.0,
      south_bound => 42.0,
      north_bound => 46.0
    },
    'stream_temperature_line_08_lower_columbia_river' => {
      begin_date  => '2000-01-01',
      end_date    => '2080-12-31',
      west_bound  => -125.0,
      east_bound  => -116.0,
      south_bound => 42.0,
      north_bound => 46.0
    },
    'stream_temperature_line_08_middle_columbia_river' => {
      begin_date  => '2000-01-01',
      end_date    => '2080-12-31',
      west_bound  => -125.0,
      east_bound  => -116.0,
      south_bound => 42.0,
      north_bound => 46.0
    },
    'stream_temperature_line_current_08_upper_rogue_river' => {
      begin_date  => '1990-08-01',
      end_date    => '2018-08-31',
      west_bound  => -122.994,
      east_bound  => -122.146,
      south_bound => 42.221,
      north_bound => 43.128
    },
    'stream_temperature_line_future_08_upper_rogue_river' => {
      begin_date  => '1990-08-01',
      end_date    => '2018-08-31',
      west_bound  => -122.994,
      east_bound  => -122.146,
      south_bound => 42.221,
      north_bound => 43.128
    },
    'stream_shade_line_08_lower_columbia_river' => {
      begin_date  => '2000-01-01',
      end_date    => '2080-12-31',
      west_bound  => -125.0,
      east_bound  => -116.0,
      south_bound => 42.0,
      north_bound => 46.0
    },
    'stream_shade_line_08_middle_columbia_river' => {
      begin_date  => '2000-01-01',
      end_date    => '2080-12-31',
      west_bound  => -125.0,
      east_bound  => -116.0,
      south_bound => 42.0,
      north_bound => 46.0
    }
  },

  'alaska' => {
    'states' => {
      begin_date  => '2018-04-01',
      end_date    => '2018-04-01',
      west_bound  => -179.23109,
      east_bound  =>  179.85968,
      south_bound =>   51.17509,
      north_bound =>   71.43979
    },
    'provinces' => {
      begin_date  => '2018-04-01',
      end_date    => '2018-04-01',
      west_bound  => -141.01807,
      east_bound  => -114.05376,
      south_bound =>   48.30892,
      north_bound =>   69.64746
    },
    'counties' => {
      begin_date  => '2018-04-01',
      end_date    => '2018-04-01',
      west_bound  => -179.23109,
      east_bound  =>  179.85968,
      south_bound =>   51.17509,
      north_bound =>   71.43979
    },
    'cities' => {
      begin_date  => '2018-04-01',
      end_date    => '2018-04-01',
      west_bound  => -165.43036,
      east_bound  => -122.70919,
      south_bound =>   48.89251,
      north_bound =>   71.30684
    },
    'roads' => {
      begin_date  => '2018-04-01',
      end_date    => '2018-04-01',
      west_bound  => -149.91720,
      east_bound  => -149.24358,
      south_bound =>   61.08498,
      north_bound =>   61.57123
    },
    'hucs' => {
      begin_date  => '2018-04-01',
      end_date    => '2018-04-01',
      west_bound  => -179.22964,
      east_bound  =>  179.85669,
      south_bound =>   51.15683,
      north_bound =>   71.43957
    },
    'legislative_districts' => {
      begin_date  => '2018-04-01',
      end_date    => '2018-04-01',
      west_bound  => -179.23109,
      east_bound  =>  179.85968,
      south_bound =>   51.17509,
      north_bound =>   71.43979
    },

    'elevation' => {
      begin_date  => '2004-01-01',
      end_date    => '2004-12-31',
      west_bound  => -180,
      east_bound  => -66,
      south_bound => 18,
      north_bound => 66
    },

    # From: https://www.ngdc.noaa.gov/dem/squareCellGrid/search

    'elevation_adak_ak_mhw' => {
      begin_date  => '2009-04-03',
      end_date    => '2009-04-03',
      west_bound  => -177.24,
      east_bound  => -175.78,
      south_bound => 51.25,
      north_bound => 52.16
    },
    'elevation_akutan_ak_mhhw' => {
      begin_date  => '2009-04-20',
      end_date    => '2009-04-20',
      west_bound  => -165.85,
      east_bound  => -165.67,
      south_bound => 54.09,
      north_bound => 54.17
    },
    'elevation_atka_ak_mhw' => {
      begin_date  => '2010-03-30',
      end_date    => '2010-03-30',
      west_bound  => -175.78,
      east_bound  => -172.9,
      south_bound => 51.59,
      north_bound => 52.49
    },
    'elevation_chenega_ak_mhhw' => {
      begin_date  => '2011-01-24',
      end_date    => '2011-01-24',
      west_bound  => -148.17,
      east_bound  => -147.89,
      south_bound => 59.97,
      north_bound => 60.19
    },
    'elevation_chignik_ak_mhw' => {
      begin_date  => '2008-10-24',
      end_date    => '2008-10-24',
      west_bound  => -158.54,
      east_bound  => -158.21,
      south_bound => 56.23,
      north_bound => 56.49
    },
    'elevation_craig_ak_mhw' => {
      begin_date  => '2012-06-29',
      end_date    => '2012-06-29',
      west_bound  => -133.91,
      east_bound  => -132.88,
      south_bound => 55.04,
      north_bound => 55.81
    },
    'elevation_dutch_harbor_ak_mhw' => {
      begin_date  => '2006-08-22',
      end_date    => '2006-08-22',
      west_bound  => -167.2,
      east_bound  => -165.9,
      south_bound => 53.5,
      north_bound => 54.35
    },
    'elevation_elfin_cove_ak_mhw' => {
      begin_date  => '2011-03-31',
      end_date    => '2011-03-31',
      west_bound  => -137.27,
      east_bound  => -135.97,
      south_bound => 57.53,
      north_bound => 58.67
    },
    'elevation_juneau_ak_mhhw' => {
      begin_date  => '2010-09-30',
      end_date    => '2010-09-30',
      west_bound  => -135.23,
      east_bound  => -133.91,
      south_bound => 58.03,
      north_bound => 58.91
    },
    'elevation_kachemak_bay_ak_mhw' => {
      begin_date  => '2010-12-07',
      end_date    => '2010-12-07',
      west_bound  => -152.1,
      east_bound  => -150.89,
      south_bound => 59.32,
      north_bound => 59.82
    },
    'elevation_king_cove_ak_mhw' => {
      begin_date  => '2013-12-31',
      end_date    => '2013-12-31',
      west_bound  => -162.8,
      east_bound  => -161.84,
      south_bound => 54.58,
      north_bound => 55.34
    },
    'elevation_kodiak_ak_mhw' => {
      begin_date  => '2013-11-04',
      end_date    => '2013-11-04',
      west_bound  => -152.8,
      east_bound  => -152,
      south_bound => 57.58,
      north_bound => 58
    },
    'elevation_nikolski_ak_mhw' => {
      begin_date  => '2010-03-15',
      end_date    => '2010-30-15',
      west_bound  => -169.66,
      east_bound  => -168.09,
      south_bound => 52.29,
      north_bound => 53.19
    },
    'elevation_prince_william_sound_ak_mhhw' => {
      begin_date  => '2009-04-20',
      end_date    => '2009-04-20',
      west_bound  => -148.77,
      east_bound  => -145.33,
      south_bound => 59.65,
      north_bound => 61.35
    },
    'elevation_sand_point_ak_mhw' => {
      begin_date  => '2006-09-29',
      end_date    => '2006-09-29',
      west_bound  => -161,
      east_bound  => -159.8,
      south_bound => 55.05,
      north_bound => 55.7
    },
    'elevation_seldovia_ak_mhw' => {
      begin_date  => '2004-06-22',
      end_date    => '2004-06-22',
      west_bound  => -151.78,
      east_bound  => -151.69,
      south_bound => 59.39,
      north_bound => 59.47
    },
    'elevation_seward_ak_mhw' => {
      begin_date  => '2005-06-22',
      end_date    => '2005-06-22',
      west_bound  => -149.47,
      east_bound  => -149.31,
      south_bound => 60.06,
      north_bound => 60.16
    },
    'elevation_shemya_ak_mhw' => {
      begin_date  => '2009-07-31',
      end_date    => '2009-07-31',
      west_bound  => -172.45,
      east_bound  => -174.21,
      south_bound => 52.27,
      north_bound => 53.06
    },
    'elevation_sitka_ak_mhw' => {
      begin_date  => '2011-07-26',
      end_date    => '2011-07-26',
      west_bound  => -135.83,
      east_bound  => -135.08,
      south_bound => 56.78,
      north_bound => 57.27
    },
    'elevation_southcentral_ak_mhw' => {
      begin_date  => '2006-08-23',
      end_date    => '2006-08-23',
      west_bound  => -156,
      east_bound  => -147,
      south_bound => 55,
      north_bound => 62
    },
    'elevation_southcentral_ak_mhw' => {
      begin_date  => '2006-08-23',
      end_date    => '2006-08-23',
      west_bound  => -156,
      east_bound  => -147,
      south_bound => 55,
      north_bound => 62
    },
    'elevation_southeast_ak_mhhw' => {
      begin_date  => '2010-08-30',
      end_date    => '2010-08-30',
      west_bound  => -137.51,
      east_bound  => -129.49,
      south_bound => 54.59,
      north_bound => 59.61
    },
    'elevation_tatitlek_ak_mhhw' => {
      begin_date  => '2014-01-20',
      end_date    => '2014-01-20',
      west_bound  => -146.79,
      east_bound  => -146.61,
      south_bound => 60.81,
      north_bound => 60.95
    },
    'elevation_unalaska_ak_mhhw' => {
      begin_date  => '2012-10-16',
      end_date    => '2012-10-16',
      west_bound  => -166.66,
      east_bound  => -166.42,
      south_bound => 53.8,
      north_bound => 53.94
    },
    'elevation_whittier_ak_mhhw' => {
      begin_date  => '2009-04-20',
      end_date    => '2009-04-20',
      west_bound  => -148.75,
      east_bound  => -148.49,
      south_bound => 60.75,
      north_bound => 60.85
    },
    'elevation_yakutat_ak_mhhw' => {
      begin_date  => '2009-05-27',
      end_date    => '2009-05-27',
      west_bound  => -139.93,
      east_bound  => -139.57,
      south_bound => 59.41,
      north_bound => 59.67
    }
  },

  'hawaii' => {
    'states' => {
      begin_date  => '2018-04-01',
      end_date    => '2018-04-01',
      west_bound  => -178.44359,
      east_bound  => -154.75579,
      south_bound =>   18.86546,
      north_bound =>   28.51727
    },
    'counties' => {
      begin_date  => '2018-04-01',
      end_date    => '2018-04-01',
      west_bound  => -178.44359,
      east_bound  => -154.75579,
      south_bound =>   18.86546,
      north_bound =>   28.51727
    },
    'cities' => {
      begin_date  => '2018-04-01',
      end_date    => '2018-04-01',
      west_bound  => -159.73181,
      east_bound  => -154.87984,
      south_bound =>   19.44399,
      north_bound =>   22.22975
    },
    'roads' => {
      begin_date  => '2018-04-01',
      end_date    => '2018-04-01',
      west_bound  => -158.10418,
      east_bound  => -157.75702,
      south_bound =>   21.27780,
      north_bound =>   21.49046
    },
    'hucs' => {
      begin_date  => '2018-04-01',
      end_date    => '2018-04-01',
      west_bound  => -160.30307,
      east_bound  => -154.75308,
      south_bound =>   18.85975,
      north_bound =>   22.28533
    },
    'legislative_districts' => {
      begin_date  => '2018-04-01',
      end_date    => '2018-04-01',
      west_bound  => -178.44359,
      east_bound  => -154.75579,
      south_bound =>   18.86546,
      north_bound =>   28.51727
    },
    'elevation' => {
      begin_date  => '2004-01-01',
      end_date    => '2004-12-31',
      west_bound  => -125,
      east_bound  => -117,
      south_bound => 32,
      north_bound => 49
    },

    # From: https://www.ngdc.noaa.gov/dem/squareCellGrid/search

    'elevation_hi_msl' => {
      begin_date  => '2004-09-01',
      end_date    => '2004-09-01',
      west_bound  => -161,
      east_bound  => -153,
      south_bound => 18,
      north_bound => 23
    },
    'elevation_hanalei_hi_mhw' => {
      begin_date  => '2011-03-21',
      end_date    => '2011-03-21',
      west_bound  => -159.91,
      east_bound  => -519.09,
      south_bound => 22.09,
      north_bound => 22.66
    },
    'elevation_hilo_hi_mhw' => {
      begin_date  => '2011-08-12',
      end_date    => '2011-08-12',
      west_bound  => -155.43,
      east_bound  => -154.5,
      south_bound => 19.21,
      north_bound => 20.22
    },
    'elevation_kahului_hi_mllw' => {
      begin_date  => '2005-09-09',
      end_date    => '2005-09-09',
      west_bound  => -156.9,
      east_bound  => -155.9,
      south_bound => 20.4,
      north_bound => 21.4
    },
    'elevation_kailua_kona_hi_mhw' => {
      begin_date  => '2006-11-21',
      end_date    => '2006-11-21',
      west_bound  => -156.2,
      east_bound  => -155.8,
      south_bound => 19.55,
      north_bound => 20.3
    },
    'elevation_kauai_hi_mhw' => {
      begin_date  => '2012-05-31',
      end_date    => '2012-05-31',
      west_bound  => -160.4,
      east_bound  => -159.09,
      south_bound => 21.55,
      north_bound => 22.4
    },
    'elevation_kawaihae_hi_mhw' => {
      begin_date  => '2011-06-24',
      end_date    => '2011-06-24',
      west_bound  => -156.37,
      east_bound  => -155.37,
      south_bound => 19.87,
      north_bound => 20.45
    },
    'elevation_keauhou_hi_mhw' => {
      begin_date  => '2011-06-24',
      end_date    => '2011-06-24',
      west_bound  => -156.31,
      east_bound  => -155.75,
      south_bound => 19.22,
      north_bound => 19.92
    },
    'elevation_lahaina_hi_mhw' => {
      begin_date  => '2007-03-30',
      end_date    => '2007-03-30',
      west_bound  => -156.9,
      east_bound  => -156.55,
      south_bound => 20.7,
      north_bound => 21.1
    },
    'elevation_nawiliwili_hi_mlw' => {
      begin_date  => '2006-01-12',
      end_date    => '2006-01-12',
      west_bound  => -159.5,
      east_bound  => -159.1,
      south_bound => 21.7,
      north_bound => 22.1
    },
    'elevation_oahu_hi_mhw' => {
      begin_date  => '2011-10-01',
      end_date    => '2011-10-01',
      west_bound  => -158.45,
      east_bound  => -157.45,
      south_bound => 21.02,
      north_bound => 22.06
    },
    'elevation_pearl_harbor_hi_mhw' => {
      begin_date  => '2006-06-29',
      end_date    => '2006-06-29',
      west_bound  => -158.07,
      east_bound  => -157.81,
      south_bound => 21.25,
      north_bound => 21.39
    }
  },

  'interior_sw' => {
    'states' => {
      begin_date  => '2018-04-01',
      end_date    => '2018-04-01',
      west_bound  => -124.70354,
      east_bound  =>  -93.50804,
      south_bound =>   25.83716,
      north_bound =>   49.00091
    },
    'counties' => {
      begin_date  => '2018-04-01',
      end_date    => '2018-04-01',
      west_bound  => -124.70354,
      east_bound  =>  -93.50804,
      south_bound =>   25.83716,
      north_bound =>   49.00091
    },
    'cities' => {
      begin_date  => '2018-04-01',
      end_date    => '2018-04-01',
      west_bound  => -124.43730,
      east_bound  =>  -93.71540,
      south_bound =>   25.85936,
      north_bound =>   48.70987
    },
    'roads' => {
      begin_date  => '2018-04-01',
      end_date    => '2018-04-01',
      west_bound  => -124.21652,
      east_bound  =>  -93.70136,
      south_bound =>   25.89511,
      north_bound =>   47.71432
    },
    'hucs' => {
      begin_date  => '2018-04-01',
      end_date    => '2018-04-01',
      west_bound  => -120.46492,
      east_bound  =>  -97.14611,
      south_bound =>   25.83707,
      north_bound =>   43.45218
    },
    'lakes' => {
      begin_date  => '2018-04-01',
      end_date    => '2018-04-01',
      west_bound  => -124.22568,
      east_bound  =>  -93.56366,
      south_bound =>   26.01253,
      north_bound =>   48.80421
    },
    'legislative_districts' => {
      begin_date  => '2018-04-01',
      end_date    => '2018-04-01',
      west_bound  => -124.70354,
      east_bound  =>  -93.50804,
      south_bound =>   25.83716,
      north_bound =>   49.00091
    },
    'flowlines_puget_sound_watershed' => {
      begin_date  => '2024-01-01',
      end_date    => '2024-12-31',
      west_bound  => -124.7272,
      east_bound  => -124.7177,
      south_bound => 48.3826,
      north_bound => 49.0394
    },
    'discharge_points_federal_puget_sound_watershed' => {
      begin_date  => '2005-01-01',
      end_date    => '2020-12-31',
      west_bound  => -124.662,
      east_bound  => -122.274,
      south_bound => 47.134,
      north_bound => 48.814
    },
    'discharge_points_state_puget_sound_watershed' => {
      begin_date  => '2005-01-01',
      end_date    => '2020-12-31',
      west_bound  => -124.297,
      east_bound  => -121.142,
      south_bound => 46.86,
      north_bound => 48.979
    }
  },

  'great_lakes' => {
    'states' => {
      begin_date  => '2018-04-01',
      end_date    => '2018-04-01',
      west_bound  => -104.05788,
      east_bound  =>  -71.77749,
      south_bound =>   36.97030,
      north_bound =>   49.38436
    },
    'provinces' => {
      begin_date  => '2018-04-01',
      end_date    => '2018-04-01',
      west_bound  => -110.00802,
      east_bound  =>  -57.10549,
      south_bound =>   41.68144,
      north_bound =>   62.58286
    },
    'counties' => {
      begin_date  => '2018-04-01',
      end_date    => '2018-04-01',
      west_bound  => -104.05788,
      east_bound  =>  -71.77749,
      south_bound =>   36.97030,
      north_bound =>   49.38436
    },
    'cities' => {
      begin_date  => '2018-04-01',
      end_date    => '2018-04-01',
      west_bound  => -103.91718,
      east_bound  =>  -71.90675,
      south_bound =>   36.96654,
      north_bound =>   48.85811
    },
    'roads' => {
      begin_date  => '2018-04-01',
      end_date    => '2018-04-01',
      west_bound  => -104.05593,
      east_bound  =>  -72.46947,
      south_bound =>   37.02361,
      north_bound =>   49.00067
    },
    'tributaries' => {
      begin_date  => '2018-04-01',
      end_date    => '2018-04-01',
      west_bound  =>  -92.31052,
      east_bound  =>  -74.37351,
      south_bound =>   41.32767,
      north_bound =>   49.04478
    },
    'hucs' => {
      begin_date  => '2018-04-01',
      end_date    => '2018-04-01',
      west_bound  => -104.79320,
      east_bound  =>  -73.39095,
      south_bound =>   40.39467,
      north_bound =>   50.16802
    },
    'watersheds' => {
      begin_date  => '2018-04-01',
      end_date    => '2018-04-01',
      west_bound  =>  -93.21372,
      east_bound  =>  -74.02793,
      south_bound =>   40.39475,
      north_bound =>   51.12923
    },
    'lakes' => {
      begin_date  => '2018-04-01',
      end_date    => '2018-04-01',
      west_bound  => -115.75729,
      east_bound  =>  -90.11437,
      south_bound =>   36.24090,
      north_bound =>   48.99995
    },
    'legislative_districts' => {
      begin_date  => '2018-04-01',
      end_date    => '2018-04-01',
      west_bound  => -104.05788,
      east_bound  =>  -71.77749,
      south_bound =>   36.97030,
      north_bound =>   49.38436
    }
  },

  'lower_mississippi' => {
    'states' => {
      begin_date  => '2018-04-01',
      end_date    => '2018-04-01',
      west_bound  => -109.06020,
      east_bound  =>  -75.16644,
      south_bound =>   25.83716,
      north_bound =>   41.00344
    },
    'counties' => {
      begin_date  => '2018-04-01',
      end_date    => '2018-04-01',
      west_bound  => -109.06020,
      east_bound  =>  -75.16644,
      south_bound =>   25.83716,
      north_bound =>   41.00344
    },
    'cities' => {
      begin_date  => '2018-04-01',
      end_date    => '2018-04-01',
      west_bound  => -109.10694,
      east_bound  =>  -75.56276,
      south_bound =>   25.85936,
      north_bound =>   40.75509
    },
    'roads' => {
      begin_date  => '2018-04-01',
      end_date    => '2018-04-01',
      west_bound  => -109.05128,
      east_bound  =>  -75.96892,
      south_bound =>   25.89511,
      north_bound =>   41.00232
    },
    'tributaries' => {
      begin_date  => '2018-04-01',
      end_date    => '2018-04-01',
      west_bound  =>  -93.83979,
      east_bound  =>  -89.12869,
      south_bound =>   29.09775,
      north_bound =>   34.99561
    },
    'hucs' => {
      begin_date  => '2018-04-01',
      end_date    => '2018-04-01',
      west_bound  => -106.59976,
      east_bound  =>  -81.26112,
      south_bound =>   28.85430,
      north_bound =>   39.38325
    },
    'watersheds' => {
      begin_date  => '2018-04-01',
      end_date    => '2018-04-01',
      west_bound  =>  -93.66833,
      east_bound  =>  -88.78032,
      south_bound =>   28.99238,
      north_bound =>   33.39162
    },
    'estuaries' => {
      begin_date  => '2018-04-01',
      end_date    => '2018-04-01',
      west_bound  =>  -93.69902,
      east_bound  =>  -88.81840,
      south_bound =>   29.00024,
      north_bound =>   30.47649
    },
    'lakes' => {
      begin_date  => '2018-04-01',
      end_date    => '2018-04-01',
      west_bound  => -108.48937,
      east_bound  =>  -75.81044,
      south_bound =>   26.01253,
      north_bound =>   40.94529
    },
    'legislative_districts' => {
      begin_date  => '2018-04-01',
      end_date    => '2018-04-01',
      west_bound  => -109.06020,
      east_bound  =>  -75.16644,
      south_bound =>   25.83716,
      north_bound =>   41.00344
    },
    'sediment_nca_2015' => {
      begin_date  => '2015-01-01',
      end_date    => '2015-12-31',
      west_bound  => -93.3,
      east_bound  => -88.9,
      south_bound => 29,
      north_bound => 30.3
    },
    'steam_discharge_nhd_line' => {
      begin_date  => '1952-01-01',
      end_date    => '2017-12-31',
      west_bound  => -93,
      east_bound  => -90,
      south_bound => 30,
      north_bound => 31
    },
    'steam_discharge_usgs_point' => {
      begin_date  => '1952-01-01',
      end_date    => '2017-12-31',
      west_bound  => -93,
      east_bound  => -90,
      south_bound => 30,
      north_bound => 31
    },
    'tide_point' => {
      begin_date  => '1950-01-01',
      end_date    => '2018-12-31',
      west_bound  => -93,
      east_bound  => -90,
      south_bound => 30,
      north_bound => 31
    },
    'tide_current' => {
      begin_date  => '2017-01-01',
      end_date    => '2019-12-31',
      west_bound  => -93,
      east_bound  => -89,
      south_bound => 29,
      north_bound => 30
    },
    'longshore_current' => {
      begin_date  => '2012-01-01',
      end_date    => '2020-12-31',
      west_bound  => -93,
      east_bound  => -89,
      south_bound => 29,
      north_bound => 30
    },
    'estuary_flushing' => {
      begin_date  => '1950-01-01',
      end_date    => '2020-12-31',
      west_bound  => -93,
      east_bound  => -89,
      south_bound => 29,
      north_bound => 30
    },
    'salinity_point' => {
      begin_date  => '1990-01-01',
      end_date    => '2018-12-31',
      west_bound  => -94,
      east_bound  => -88,
      south_bound => 29,
      north_bound => 31
    }
  },

  'upper_mississippi' => {
    'states' => {
      begin_date  => '2018-04-01',
      end_date    => '2018-04-01',
      west_bound  => -104.05788,
      east_bound  =>  -71.77749,
      south_bound =>   33.75288,
      north_bound =>   49.38436
    },
    'counties' => {
      begin_date  => '2018-04-01',
      end_date    => '2018-04-01',
      west_bound  => -104.05788,
      east_bound  =>  -71.77749,
      south_bound =>   33.75288,
      north_bound =>   49.38436
    },
    'cities' => {
      begin_date  => '2018-04-01',
      end_date    => '2018-04-01',
      west_bound  => -103.91718,
      east_bound  =>  -71.90675,
      south_bound =>   33.87338,
      north_bound =>   48.85811
    },
    'roads' => {
      begin_date  => '2018-04-01',
      end_date    => '2018-04-01',
      west_bound  => -104.05593,
      east_bound  =>  -72.46947,
      south_bound =>   34.19382,
      north_bound =>   49.00067
    },
    'hucs' => {
      begin_date  => '2018-04-01',
      end_date    => '2018-04-01',
      west_bound  =>  -97.36394,
      east_bound  =>  -77.83937,
      south_bound =>   35.31332,
      north_bound =>   47.77165
    },
    'lakes' => {
      begin_date  => '2018-04-01',
      end_date    => '2018-04-01',
      west_bound  => -103.76586,
      east_bound  =>  -71.90773,
      south_bound =>   34.25690,
      north_bound =>   49.37019
    },
    'legislative_districts' => {
      begin_date  => '2018-04-01',
      end_date    => '2018-04-01',
      west_bound  => -104.05788,
      east_bound  =>  -71.77749,
      south_bound =>   33.75288,
      north_bound =>   49.38436
    }
  },

  'missouri_river' => {
    'states' => {
      begin_date  => '2018-04-01',
      end_date    => '2018-04-01',
      west_bound  => -116.04915,
      east_bound  =>  -89.09897,
      south_bound =>   35.99568,
      north_bound =>   49.38436
    },
    'provinces' => {
      begin_date  => '2018-04-01',
      end_date    => '2018-04-01',
      west_bound  => -120.00270,
      east_bound  => -101.36230,
      south_bound =>   48.99674,
      north_bound =>   60.00007
    },
    'counties' => {
      begin_date  => '2018-04-01',
      end_date    => '2018-04-01',
      west_bound  => -116.04915,
      east_bound  =>  -89.09897,
      south_bound =>   35.99568,
      north_bound =>   49.38436
    },
    'cities' => {
      begin_date  => '2018-04-01',
      end_date    => '2018-04-01',
      west_bound  => -115.59728,
      east_bound  =>  -89.32832,
      south_bound =>   36.16184,
      north_bound =>   48.85811
    },
    'roads' => {
      begin_date  => '2018-04-01',
      end_date    => '2018-04-01',
      west_bound  => -115.69476,
      east_bound  =>  -89.21127,
      south_bound =>   35.99957,
      north_bound =>   49.00067
    },
    'hucs' => {
      begin_date  => '2018-04-01',
      end_date    => '2018-04-01',
      west_bound  => -114.16663,
      east_bound  =>  -90.11421,
      south_bound =>   37.02460,
      north_bound =>   50.00105
    },
    'lakes' => {
      begin_date  => '2018-04-01',
      end_date    => '2018-04-01',
      west_bound  => -115.75729,
      east_bound  =>  -90.11437,
      south_bound =>   36.24090,
      north_bound =>   48.99995
    },
    'legislative_districts' => {
      begin_date  => '2018-04-01',
      end_date    => '2018-04-01',
      west_bound  => -116.04915,
      east_bound  =>  -89.09897,
      south_bound =>   35.99568,
      north_bound =>   49.38436
    }
  }
);



# @coverages is used by parse_coverage() whereas the other arrays are used by
# print_coverage_offering_brief() and print_coverage_offering().
# These should match and somehow be merged. UGLY.

my @coverages = (

  'atlantic.states',
  'gulf.states',
  'pacific.states',
  'alaska.states',
  'hawaii.states',
  'interior_sw.states',
  'great_lakes.states',
  'missouri_river.states',
  'lower_mississippi.states',
  'upper_mississippi.states',

  'atlantic.provinces',
  'pacific.provinces',
  'great_lakes.provinces',
  'missouri_river.provinces',

  'atlantic.counties',
  'gulf.counties',
  'pacific.counties',
  'alaska.counties',
  'hawaii.counties',
  'interior_sw.counties',
  'great_lakes.counties',
  'missouri_river.counties',
  'lower_mississippi.counties',
  'upper_mississippi.counties',

  'atlantic.cities',
  'gulf.cities',
  'pacific.cities',
  'alaska.cities',
  'hawaii.cities',
  'interior_sw.cities',
  'great_lakes.cities',
  'missouri_river.cities',
  'lower_mississippi.cities',
  'upper_mississippi.cities',

  'atlantic.roads',
  'gulf.roads',
  'pacific.roads',
  'alaska.roads',
  'hawaii.roads',
  'interior_sw.roads',
  'great_lakes.roads',
  'missouri_river.roads',
  'lower_mississippi.roads',
  'upper_mississippi.roads',

  'atlantic.tributaries',
  'gulf.tributaries',
  'pacific.tributaries',
  'great_lakes.tributaries',
  'lower_mississippi.tributaries',

  'atlantic.hucs',
  'gulf.hucs',
  'pacific.hucs',
  'alaska.hucs',
  'hawaii.hucs',
  'interior_sw.hucs',
  'great_lakes.hucs',
  'missouri_river.hucs',
  'lower_mississippi.hucs',
  'upper_mississippi.hucs',

  'atlantic.watersheds',
  'gulf.watersheds',
  'pacific.watersheds',
  'great_lakes.watersheds',
  'lower_mississippi.watersheds',

  'atlantic.estuaries',
  'gulf.estuaries',
  'pacific.estuaries',
  'lower_mississippi.estuaries',

  'atlantic.lakes',
  'gulf.lakes',
  'pacific.lakes',
  'interior_sw.lakes',
  'great_lakes.lakes',
  'missouri_river.lakes',
  'lower_mississippi.lakes',
  'upper_mississippi.lakes',

  'atlantic.coastal_zone_management_areas',
  'gulf.coastal_zone_management_areas',
  'pacific.coastal_zone_management_areas',
  'alaska.coastal_zone_management_areas',
  'hawaii.coastal_zone_management_areas',
  'interior_sw.coastal_zone_management_areas',
  'lower_mississippi.coastal_zone_management_areas',
  'upper_mississippi.coastal_zone_management_areas',

  'atlantic.legislative_districts',
  'gulf.legislative_districts',
  'pacific.legislative_districts',
  'alaska.legislative_districts',
  'hawaii.legislative_districts',
  'interior_sw.legislative_districts',
  'great_lakes.legislative_districts',
  'missouri_river.legislative_districts',
  'lower_mississippi.legislative_districts',
  'upper_mississippi.legislative_districts',

  # Data:

  'atlantic.nlcd_2011',
  'atlantic.nlcd_2006',
  'atlantic.nlcd_2001',
  'atlantic.nlcd_1992',
  'atlantic.nlcd_change_1992_2001',
  'atlantic.nlcd_change_2001_2006',
  'atlantic.ccap_2006',
  'atlantic.ccap_2001',
  'atlantic.decadal_iclus_land_use_2000_2100',
  'atlantic.land_use_iclus_watershed',
  'atlantic.land_use_nlcd_2011_watershed',
  'atlantic.land_use_nlcd_2006_watershed',
  'atlantic.impervious_2006',
  'atlantic.impervious_2001',
  'atlantic.impervious_change_2001_2006',
  'atlantic.decadal_iclus_a1_impervious_2010_2100',
  'atlantic.decadal_iclus_a2_impervious_2010_2100',
  'atlantic.decadal_iclus_b1_impervious_2010_2100',
  'atlantic.decadal_iclus_b2_impervious_2010_2100',
  'atlantic.decadal_iclus_bc_impervious_2010_2100',

  'gulf.nlcd_2011',
  'gulf.nlcd_2006',
  'gulf.nlcd_2001',
  'gulf.nlcd_1992',
  'gulf.nlcd_change_1992_2001',
  'gulf.nlcd_change_2001_2006',
  'gulf.ccap_2006',
  'gulf.ccap_2001',
  'gulf.decadal_iclus_land_use_2000_2100',
  'gulf.land_use_iclus_watershed',
  'gulf.land_use_nlcd_2011_watershed',
  'gulf.land_use_nlcd_2006_watershed',
  'gulf.impervious_2006',
  'gulf.impervious_2001',
  'gulf.impervious_change_2001_2006',
  'gulf.decadal_iclus_a1_impervious_2010_2100',
  'gulf.decadal_iclus_a2_impervious_2010_2100',
  'gulf.decadal_iclus_b1_impervious_2010_2100',
  'gulf.decadal_iclus_b2_impervious_2010_2100',
  'gulf.decadal_iclus_bc_impervious_2010_2100',

  'pacific.nlcd_2011',
  'pacific.nlcd_2006',
  'pacific.nlcd_2001',
  'pacific.nlcd_1992',
  'pacific.nlcd_change_1992_2001',
  'pacific.nlcd_change_2001_2006',
  'pacific.ccap_2006',
  'pacific.ccap_2001',
  'pacific.decadal_iclus_land_use_2000_2100',
  'pacific.land_use_iclus_watershed',
  'pacific.land_use_nlcd_2011_watershed',
  'pacific.land_use_nlcd_2006_watershed',
  'pacific.impervious_2006',
  'pacific.impervious_2001',
  'pacific.impervious_change_2001_2006',
  'pacific.decadal_iclus_a1_impervious_2010_2100',
  'pacific.decadal_iclus_a2_impervious_2010_2100',
  'pacific.decadal_iclus_b1_impervious_2010_2100',
  'pacific.decadal_iclus_b2_impervious_2010_2100',
  'pacific.decadal_iclus_bc_impervious_2010_2100',

  'atlantic.elevation',
  'gulf.elevation',
  'pacific.elevation',
  'alaska.elevation',
  'hawaii.elevation',

  'elevation_adak_ak_mhw',
  'elevation_agua_hedionda_creek_ca_mhw',
  'elevation_akutan_ak_mhhw',
  'elevation_albemarle_sound_nc_mllw',
  'elevation_albion_river_ca_mhw',
  'elevation_alder_creek_ca_mhw',
  'elevation_aliso_canyon_ca_mhw',
  'elevation_aliso_creek_ca_mhw',
  'elevation_alligator_cove_fl_navd88',
  'elevation_alsea_bay_or_mhw',
  'elevation_altahama_sound_ga_navd88',
  'elevation_alva_paul_creek_ca_mhw',
  'elevation_andree_clark_bird_refuge_ca_mhw',
  'elevation_annisquam_river_ma_mllw',
  'elevation_apalachee_bay_fl_mllw',
  'elevation_apalachicola_bay_fl_mllw',
  'elevation_aransas_bay_1_tx_mllw',
  'elevation_aransas_bay_2_tx_mllw',
  'elevation_arena_cove_ca_mhw',
  'elevation_arroyo_burro_ca_mhw',
  'elevation_arroyo_de_la_cruz_ca_mhw',
  'elevation_arroyo_del_corral_ca_mhw',
  'elevation_arroyo_del_puerto_ca_mhw',
  'elevation_arroyo_grande_creek_lagoon_ca_mhw',
  'elevation_arroyo_laguna_ca_mhw',
  'elevation_arthur_kill_nj_mllw',
  'elevation_arthur_kill_ny_mllw',
  'elevation_assawoman_bay_de_mhw',
  'elevation_assawoman_bay_md_mhw',
  'elevation_astoria_or_mhw',
  'elevation_atka_ak_mhw',
  'elevation_atlantic_city_nj_mhw',
  'elevation_baffin_bay_tx_mllw',
  'elevation_baldwin_creek_ca_mhw',
  'elevation_bar_harbor_me_mhw',
  'elevation_bar_harbor_me_navd88',
  'elevation_barnegat_bay_nj_mllw',
  'elevation_barnstable_bay_ma_navd88',
  'elevation_bc_msl',
  'elevation_bear_river_ca_mhw',
  'elevation_beaver_creek_or_mhw',
  'elevation_bell_canyon_creek_ca_mhw',
  'elevation_bellingham_bay_wa_mhw',
  'elevation_berry_creek_or_mhw',
  'elevation_big_creek_1_or_mhw',
  'elevation_big_creek_or_mhw',
  'elevation_big_lagoon_ca_mhw',
  'elevation_big_river_ca_mhw',
  'elevation_big_spring_creek_fl_mllw',
  'elevation_big_sur_river_ca_mhw',
  'elevation_biloxi_ms_mhw',
  'elevation_birch_bay_wa_mllw',
  'elevation_biscayne_bay_fl_navd88',
  'elevation_blackwater_sound_fl_navd88',
  'elevation_blue_hill_bay_me_navd88',
  'elevation_bob_creek_or_mhw',
  'elevation_bodega_bay_ca_mhw',
  'elevation_boggy_bay_fl_mllw',
  'elevation_boggy_creek_fl_mllw',
  'elevation_bogue_sound_nc_navd88',
  'elevation_bolinas_lagoon_ca_mhw',
  'elevation_bolsa_chica_ecological_reserve_ca_mhw',
  'elevation_boone_creek_wa_mhw',
  'elevation_boston_harbor_ma_mllw',
  'elevation_brave_boat_harbor_me_navd88',
  'elevation_brazos_river_tx_navd88',
  'elevation_broad_river_rodgers_river_fl_navd88',
  'elevation_broad_river_sc_navd88',
  'elevation_broad_sound_ma_mllw',
  'elevation_brush_creek_1_or_mhw',
  'elevation_brush_creek_2_ca_mhw',
  'elevation_buck_creek_fl_mllw',
  'elevation_buena_vista_lagoon_ca_mhw',
  'elevation_bulls_bay_sc_navd88',
  'elevation_bumblebee_creek_fl_mllw',
  'elevation_bush_river_md_navd88',
  'elevation_buzzards_bay_ma_mllw',
  'elevation_calera_creek_ca_mhw',
  'elevation_camp_creek_wa_mhw',
  'elevation_canada_de_la_gaviota_creek_ca_mhw',
  'elevation_canada_del_refugio_ca_mhw',
  'elevation_cape_creek_or_mhw',
  'elevation_cape_fear_river_nc_navd88',
  'elevation_cape_hatteras_nc_mhw',
  'elevation_cape_romaine_refuge_sc_navd88',
  'elevation_capers_creek_sc_navd88',
  'elevation_carmel_ca_mhw',
  'elevation_carpinteria_creek_ca_mhw',
  'elevation_carpinteria_salt_marsh_ca_mhw',
  'elevation_casco_bay_me_mhw',
  'elevation_casper_creek_ca_mhw',
  'elevation_cayucos_creek_ca_mhw',
  'elevation_cedar_creek_wa_mhw',
  'elevation_central_ca_mhw',
  'elevation_central_ca_navd88',
  'elevation_central_fl_mhw',
  'elevation_central_fl_navd88',
  'elevation_central_gulf_msl',
  'elevation_central_or_mhw',
  'elevation_central_pacific_msl',
  'elevation_chadbourne_gulch_ca_mhw',
  'elevation_channel_islands_harbor_ca_mhw',
  'elevation_charleston_harbor_sc_navd88',
  'elevation_charlotte_harbor_fl_mllw',
  'elevation_chenega_ak_mhhw',
  'elevation_chesapeake_bay_mainstem_1_md_mllw',
  'elevation_chesapeake_bay_mainstem_1_va_mllw',
  'elevation_chesapeake_bay_mainstem_2_md_mllw',
  'elevation_chesapeake_bay_mainstem_2_va_mllw',
  'elevation_chester_river_md_navd88',
  'elevation_chetco_river_or_mhw',
  'elevation_chignik_ak_mhw',
  'elevation_china_creek_or_mhw',
  'elevation_chincoteague_bay_md_mhw',
  'elevation_chincoteague_bay_va_mhw',
  'elevation_choctawhatchee_bay_fl_mhw',
  'elevation_choctawhatchee_bay_fl_mhw',
  'elevation_choptank_river_md_navd88',
  'elevation_clam_beach_ca_mhw',
  'elevation_clay_creek_fl_mllw',
  'elevation_clearwater_saint_joseph_sound_fl_mhw',
  'elevation_cobscook_me_mllw',
  'elevation_cockelburr_canyon_ca_mhw',
  'elevation_cocohatchee_river_fl_navd88',
  'elevation_columbia_river_or_mllw',
  'elevation_columbia_river_wa_mllw',
  'elevation_conner_creek_wa_mhw',
  'elevation_coos_bay_or_mhw',
  'elevation_copalis_river_wa_mhw',
  'elevation_coquille_river_or_mhw',
  'elevation_corpus_christi_bay_tx_mhw',
  'elevation_corpus_christi_tx_mhw',
  'elevation_cottaneva_creek_ca_mhw',
  'elevation_craig_ak_mhw',
  'elevation_crescent_city_ca_mhw',
  'elevation_crescent_city_ca_navd88',
  'elevation_crescent_city_harbor_ca_mhw',
  'elevation_crooked_creek_fl_mllw',
  'elevation_crystal_river_1_fl_mllw',
  'elevation_crystal_river_2_fl_mllw',
  'elevation_crystal_river_3_fl_mllw',
  'elevation_cumberland_sound_fl_navd88',
  'elevation_cumberland_sound_ga_navd88',
  'elevation_cummins_creek_or_mhw',
  'elevation_currituck_sound_nc_navd88',
  'elevation_currituck_sound_va_navd88',
  'elevation_cut_creek_or_mhw',
  'elevation_dairy_gulch_ca_mhw',
  'elevation_daley_lake_or_mhw',
  'elevation_dallus_creek_fl_mllw',
  'elevation_damariscaotta_me_mllw',
  'elevation_dana_point_harbor_ca_mhw',
  'elevation_days_creek_fl_mllw',
  'elevation_daytona_beach_fl_mhw',
  'elevation_deadman_bay_fl_mllw',
  'elevation_deboy_sound_ga_navd88',
  'elevation_delaware_bay_1_de_mllw',
  'elevation_delaware_bay_1_nj_mllw',
  'elevation_delaware_bay_2_de_mllw',
  'elevation_delaware_bay_2_nj_mllw',
  'elevation_delaware_bay_3_de_mllw',
  'elevation_delaware_bay_3_nj_mllw',
  'elevation_depoe_bay_or_mhw',
  'elevation_devereux_lagoon_ca_mhw',
  'elevation_devils_lake_or_mhw',
  'elevation_dewes_hamlin_creek_sc_navd88',
  'elevation_drakes_estero_ca_mhw',
  'elevation_drayton_harbor_wa_mllw',
  'elevation_duck_creek_wa_mhw',
  'elevation_dutch_harbor_ak_mhw',
  'elevation_duxbury_bay_ma_mhw',
  'elevation_eaglenest_creek_fl_mllw',
  'elevation_eastern_bay_md_navd88',
  'elevation_eastern_gulf_msl',
  'elevation_ecola_creek_or_mhw',
  'elevation_eel_river_ca_mhw',
  'elevation_eleven_prong_bay_fl_mllw',
  'elevation_elfin_cove_ak_mhw',
  'elevation_elk_creek_ca_mhw',
  'elevation_elk_creek_wa_mhw',
  'elevation_elk_river_2_wa_mllw',
  'elevation_elk_river_3_or_mllw',
  'elevation_elk_river_md_navd88',
  'elevation_espa_lagoon_ca_mhw',
  'elevation_essex_bay_ma_mllw',
  'elevation_estero_bay_fl_navd88',
  'elevation_euchre_creek_or_mhw',
  'elevation_eureka_ca_mhw',
  'elevation_eureka_ca_navd88',
  'elevation_fern_canyon_ca_mhw',
  'elevation_fillman_bayou_fl_mllw',
  'elevation_first_bay_fl_navd88',
  'elevation_fish_creek_fl_mllw',
  'elevation_floras_creek_2_or_mhw',
  'elevation_florida_bay_fl_mllw',
  'elevation_fort_bragg_ca_mhw',
  'elevation_fort_bragg_ca_navd88',
  'elevation_french_canyon_2_ca_mhw',
  'elevation_frenchmans_bay_me_mhw',
  'elevation_frenchmans_creek_ca_mhw',
  'elevation_galvaston_bay_1_tx_mhw',
  'elevation_galvaston_bay_2_tx_mhw',
  'elevation_galveston_tx_mhw',
  'elevation_garcia_river_ca_mhw',
  'elevation_gardiner_bay_peconic_bay_1_ny_mhw',
  'elevation_gardiner_bay_peconic_bay_2_ny_mhw',
  'elevation_garibaldi_or_mhw',
  'elevation_garibaldi_or_navd88',
  'elevation_garrapata_creek_ca_mhw',
  'elevation_gazos_creek_ca_mhw',
  'elevation_glouchester_harbor_ma_mhw',
  'elevation_goleta_slough_ca_mhw',
  'elevation_goodman_creek_wa_mhw',
  'elevation_goosefare_bay_me_mhw',
  'elevation_gouldsboro_bay_me_mhw',
  'elevation_grays_harbor_wa_mllw',
  'elevation_great_bay_me_mllw',
  'elevation_great_bay_nh_mllw',
  'elevation_great_bay_nj_mhw',
  'elevation_great_south_bay_ny_mllw',
  'elevation_great_wilmico_fleets_bay_va_navd88',
  'elevation_greenwood_creek_ca_mhw',
  'elevation_gualala_river_ca_mhw',
  'elevation_gunpowder_river_md_mhw',
  'elevation_guthrie_creek_ca_mhw',
  'elevation_half_moon_bay_airport_lagoon_ca_mhw',
  'elevation_hampton_harbor_nh_mllw',
  'elevation_hanalei_hi_mhw',
  'elevation_hare_creek_ca_mhw',
  'elevation_haycock_harbor_baileys_mistake_me_mllw',
  'elevation_hi_msl',
  'elevation_hilo_hi_mhw',
  'elevation_hoh_river_wa_mhw',
  'elevation_honda_creek_lagoon_ca_mhw',
  'elevation_honga_river_md_navd88',
  'elevation_horseshoe_cove_fl_mllw',
  'elevation_hubbard_creek_or_mhw',
  'elevation_hudson_river_nj_mllw',
  'elevation_hudson_river_ny_mllw',
  'elevation_humboldt_bay_ca_mhw',
  'elevation_humboldt_ca_mhw',
  'elevation_huntington_channel_ca_mhw',
  'elevation_icw4_fl_navd88',
  'elevation_icw_wolf_mobile_al_mhw',
  'elevation_indian_bay_fl_mllw',
  'elevation_indian_river_bay_de_mhw',
  'elevation_indian_river_bay_fl_mhw',
  'elevation_inglenook_fen_ca_mhw',
  'elevation_irish_gulch_ca_mhw',
  'elevation_jackass_creek_ca_mhw',
  'elevation_jalama_creek_ca_mhw',
  'elevation_jamaica_bay_ny_mhw',
  'elevation_james_river_va_navd88',
  'elevation_jekyll_sound_saint_andrews_sound_ga_navd88',
  'elevation_joe_creek_wa_mhw',
  'elevation_johnson_creek_3_or_mhw',
  'elevation_johnson_creek_fl_mllw',
  'elevation_juneau_ak_mhhw',
  'elevation_kachemak_bay_ak_mhw',
  'elevation_kahului_hi_mllw',
  'elevation_kailua_kona_hi_mhw',
  'elevation_kalaloch_creek_wa_mhw',
  'elevation_katama_bay_ma_mhw',
  'elevation_kauai_hi_mhw',
  'elevation_kawaihae_hi_mhw',
  'elevation_keauhou_hi_mhw',
  'elevation_kennebec_river_estuary_me_mllw',
  'elevation_key_west_fl_mhw',
  'elevation_key_west_fl_navd88',
  'elevation_king_cove_ak_mhw',
  'elevation_klamath_river_ca_mhw',
  'elevation_kodiak_ak_mhw',
  'elevation_la_push_wa_mhw',
  'elevation_lagoon_creek_ca_mhw',
  'elevation_laguna_creek_ca_mhw',
  'elevation_lahaina_hi_mhw',
  'elevation_lake_cleone_ca_mhw',
  'elevation_lake_earl_ca_mhw',
  'elevation_lake_lucerne_ca_mhw',
  'elevation_lake_lytle_or_mhw',
  'elevation_lake_washington_wa_llw',
  'elevation_lake_worth_fl_mhw',
  'elevation_lakes_reeds_bay_nj_mhw',
  'elevation_las_flores_creek_ca_mhw',
  'elevation_las_pulgas_creek_ca_mhw',
  'elevation_little_creek_1_or_mhw',
  'elevation_little_creek_3_ca_mhw',
  'elevation_little_kennebec_englishman_bay_me_mllw',
  'elevation_little_lagoon_al_mhw',
  'elevation_little_narragansett_bay_ct_mhw',
  'elevation_little_pico_creek_ca_mhw',
  'elevation_little_pine_1_bay_fl_mllw',
  'elevation_little_river_ca_mhw',
  'elevation_little_river_sc_mhw',
  'elevation_little_rocky_creek_fl_mllw',
  'elevation_little_salmon_creek_ca_mhw',
  'elevation_little_spring_creek_fl_mllw',
  'elevation_little_sur_lagoon_ca_mhw',
  'elevation_lockwoods_folly_river_nc_navd88',
  'elevation_loma_alta_marsh_ca_mhw',
  'elevation_lombardi_creek_ca_mhw',
  'elevation_long_island_sound_ct_mllw',
  'elevation_long_island_sound_ny_mllw',
  'elevation_los_penasquitos_lagoon_ca_mhw',
  'elevation_lows_bay_fl_mllw',
  'elevation_loxahatchee_river_fl_mhw',
  'elevation_lunada_bay_ca_mhw',
  'elevation_lynnhaven_river_va_mhw',
  'elevation_machias_bay_me_mllw',
  'elevation_mad_river_ca_mhw',
  'elevation_maddaket_harbor_ma_mhw',
  'elevation_magothy_river_md_navd88',
  'elevation_makah_bay_wa_mhw',
  'elevation_malibu_lagoon_ca_mhw',
  'elevation_malpaso_creek_ca_mhw',
  'elevation_manasquan_river_nj_mhw',
  'elevation_manchester_ca_mhw',
  'elevation_matagorda_bay_tx_mllw',
  'elevation_mattole_river_ca_mhw',
  'elevation_mcnutt_gulch_ca_mhw',
  'elevation_merrimack_river_ma_mllw',
  'elevation_miami_fl_mhw',
  'elevation_miami_fl_navd88',
  'elevation_middle_and_black_rivers_md_mhw',
  'elevation_mississippi_sound_al_mllw',
  'elevation_mississippi_sound_la_mllw',
  'elevation_mississippi_sound_ms_mllw',
  'elevation_misson_bay_lagoon_ca_mhw',
  'elevation_misson_creek_lagoon_ca_mhw',
  'elevation_mitchell_creek_ca_mhw',
  'elevation_mobile_al_mhw',
  'elevation_mobile_al_navd88',
  'elevation_mobile_bay_al_mhw',
  'elevation_moclips_river_wa_mhw',
  'elevation_mojack_bay_va_mhw',
  'elevation_montauk_ny_mhw',
  'elevation_monterey_bay_ca_mhw',
  'elevation_monterey_ca_mhw',
  'elevation_monterey_ca_navd88',
  'elevation_moose_cove_me_mllw',
  'elevation_moran_lake_ca_mhw',
  'elevation_morehead_city_nc_mhw',
  'elevation_morehead_city_nc_navd88',
  'elevation_morro_bay_ca_mhw',
  'elevation_morro_creek_lagoon_ca_mhw',
  'elevation_mosquito_creek_wa_mhw',
  'elevation_mugu_lagoon_ca_mhw',
  'elevation_murrells_inlet_sc_mhw',
  'elevation_muscongus_bay_me_mllw',
  'elevation_mussel_creek_or_mhw',
  'elevation_myers_creek_or_mhw',
  'elevation_myrtle_beach_sc_mhw',
  'elevation_nantucket_harbor_ma_mhw',
  'elevation_nantucket_ma_mhw',
  'elevation_naples_bay_fl_navd88',
  'elevation_narragansett_bay_1_ma_navd88',
  'elevation_narragansett_bay_1_ri_navd88',
  'elevation_narragansett_bay_2_ma_navd88',
  'elevation_narragansett_bay_2_ri_navd88',
  'elevation_narragansett_bay_3_ma_navd88',
  'elevation_narragansett_bay_3_ri_navd88',
  'elevation_narraguagas_bay_me_mhw',
  'elevation_nassau_sound_fl_navd88',
  'elevation_navarro_river_ca_mhw',
  'elevation_navesink_nj_navd88',
  'elevation_nawiliwili_hi_mlw',
  'elevation_nc9_crescent_beach_ca_mhw',
  'elevation_nc_intercoastal_waterway_nc_navd88',
  'elevation_necanicum_river_estuary_or_mllw',
  'elevation_nehalem_bay_or_mhw',
  'elevation_neskowin_creek_or_mhw',
  'elevation_nestucca_bay_or_mhw',
  'elevation_netarts_bay_or_mhw',
  'elevation_new_orleans_la_mhw',
  'elevation_new_orleans_la_mllw',
  'elevation_new_orleans_la_navd88',
  'elevation_new_river_2_or_mhw',
  'elevation_new_river_nc_nad88',
  'elevation_new_river_nc_navd88',
  'elevation_new_river_sc_mhw',
  'elevation_new_river_sc_mhw',
  'elevation_new_york_harbor_1_nj_navd88',
  'elevation_new_york_harbor_1_ny_navd88',
  'elevation_new_york_harbor_2_nj_navd88',
  'elevation_new_york_harbor_2_ny_navd88',
  'elevation_newark_bay_kill_van_kull_nj_navd88',
  'elevation_newark_bay_kill_van_kull_ny_navd88',
  'elevation_newport_harbor_ca_mhw',
  'elevation_newport_or_navd88',
  'elevation_nikolski_ak_mhw',
  'elevation_nj_coastal_bays_great_egg_harbor_nj_mhw',
  'elevation_north_edisto_river_sc_navd88',
  'elevation_north_south_rivers_ma_mllw',
  'elevation_north_south_santee_rivers_sc_navd88',
  'elevation_northeast_atlantic_msl',
  'elevation_northeast_river_md_navd88',
  'elevation_northern_ca_msl',
  'elevation_northern_gulf_navd88',
  'elevation_northern_or_mhw',
  'elevation_northwest_pacific_msl',
  'elevation_noyo_river_ca_mhw',
  'elevation_oahu_hi_mhw',
  'elevation_ocean_city_md_mhw',
  'elevation_oceanside_harbor_ca_mhw',
  'elevation_old_creek_ca_mhw',
  'elevation_ormond_beach_ca_mhw',
  'elevation_ossabaw_sound_ga_mhw',
  'elevation_ossagon_creek_2_ca_mhw',
  'elevation_otter_cove_me_mhw',
  'elevation_outer_doctors_bay_fl_navd88',
  'elevation_oyster_creek_fl_mllw',
  'elevation_oyster_creek_tx_navd88',
  'elevation_ozette_river_wa_mhw',
  'elevation_padilla_bay_wa_mhw',
  'elevation_palm_beach_fl_mhw',
  'elevation_palm_beach_fl_navd88',
  'elevation_pamlico_sound_nc_mhw',
  'elevation_panama_city_fl_mhw',
  'elevation_panama_city_fl_navd88',
  'elevation_passamaquoddy_bay_st_croix_river_me_mllw',
  'elevation_patapsco_river_md_mhw',
  'elevation_patuxent_river_md_navd88',
  'elevation_pearl_harbor_hi_mhw',
  'elevation_pebble_beach_ca_mhw',
  'elevation_penobscot_bay_me_mllw',
  'elevation_perdido_bay_al_mhw',
  'elevation_perdido_bay_fl_mhw',
  'elevation_pescadero_creek_ca_mhw',
  'elevation_piankatank_river_va_mhw',
  'elevation_pico_creek_ca_mhw',
  'elevation_pilarcitos_creek_ca_mhw',
  'elevation_pillar_point_harbor_ca_mhw',
  'elevation_pismo_creek_lagoon_ca_mhw',
  'elevation_pistol_river_or_mhw',
  'elevation_pithlachascotee_river_fl_mllw',
  'elevation_pleasant_bay_ma_mllw',
  'elevation_pleasant_bay_me_mhw',
  'elevation_plum_island_sound_ma_mllw',
  'elevation_pomponio_creek_ca_mhw',
  'elevation_ponce_de_leon_inlet_fl_navd88',
  'elevation_poquoson_back_rivers_va_mhw',
  'elevation_port_alexander_ak_mhw',
  'elevation_port_hueneme_ca_mhw',
  'elevation_port_orford_head_or_mhw',
  'elevation_port_orford_or_mhw',
  'elevation_port_san_luis_ca_mhw',
  'elevation_port_townsend_wa_mhw',
  'elevation_port_townsend_wa_navd88',
  'elevation_portland_me_mhw',
  'elevation_potomac_river_md_mllw',
  'elevation_potomac_river_va_mllw',
  'elevation_prince_william_sound_ak_mhhw',
  'elevation_prospect_harbor_me_mhw',
  'elevation_pudding_creek_ca_mhw',
  'elevation_puget_sound_wa_mhw',
  'elevation_puget_sound_wa_mllw',
  'elevation_puget_sound_wa_navd88',
  'elevation_queets_river_wa_mhw',
  'elevation_quillayute_river_wa_mhw',
  'elevation_quinault_river_wa_mhw',
  'elevation_quoddy_narrows_me_mllw',
  'elevation_raft_river_wa_mhw',
  'elevation_rappahannock_river_va_navd88',
  'elevation_raritan_bay_nj_mhw',
  'elevation_raritan_bay_ny_mhw',
  'elevation_redwood_creek_big_lagoon_ca_mhw',
  'elevation_redwood_creek_ca_mhw',
  'elevation_rehoboth_bay_de_mhw',
  'elevation_rincon_creek_ca_mhw',
  'elevation_rock_1_bay_fl_mllw',
  'elevation_rock_creek_or_mhw',
  'elevation_rodeo_lagoon_ca_mhw',
  'elevation_rogue_river_or_mhw',
  'elevation_russian_gulch_ca_mhw',
  'elevation_russian_gulch_creek_ca_mhw',
  'elevation_russian_river_ca_mhw',
  'elevation_sabine_lake_la_mllw',
  'elevation_sabine_lake_tx_mllw',
  'elevation_saco_bay_me_mhw',
  'elevation_saint_andrew_bay_fl_navd88',
  'elevation_saint_augustine_fl_navd88',
  'elevation_saint_johns_river_fl_navd88',
  'elevation_saint_joseph_bay_fl_mhw',
  'elevation_salem_sound_ma_mhw',
  'elevation_salmon_creek_ca_mhw',
  'elevation_salmon_river_estuary_or_mhw',
  'elevation_salt_creek_2_fl_mllw',
  'elevation_salt_creek_fl_mllw',
  'elevation_samish_bay_wa_mhw',
  'elevation_san_antonio_bay_tx_mllw',
  'elevation_san_antonio_creek_ca_mhw',
  'elevation_san_carpoforo_creek_ca_mhw',
  'elevation_san_diego_bay_ca_mhw',
  'elevation_san_diego_ca_mhw',
  'elevation_san_diego_ca_navd88',
  'elevation_san_dieguito_river_ca_mhw',
  'elevation_san_elijo_lagoon_ca_mhw',
  'elevation_san_francisco_bay_1_ca_mhw',
  'elevation_san_francisco_bay_2_ca_mhw',
  'elevation_san_francisco_bay_3_ca_mhw',
  'elevation_san_francisco_bay_4_ca_mhw',
  'elevation_san_francisco_bay_5_ca_mhw',
  'elevation_san_francisco_ca_mhw',
  'elevation_san_francisco_ca_navd88',
  'elevation_san_gregorio_creek_ca_mhw',
  'elevation_san_juan_creek_ca_mhw',
  'elevation_san_luiz_obispo_bay_ca_mhw',
  'elevation_san_luiz_rey_river_ca_mhw',
  'elevation_san_mateo_lagoon_ca_mhw',
  'elevation_san_onofre_creek_ca_mhw',
  'elevation_san_pedro_bay_ca_mhw',
  'elevation_san_pedro_creek_ca_mhw',
  'elevation_san_simeon_creek_ca_mhw',
  'elevation_sand_lake_estuary_or_mhw',
  'elevation_sand_point_ak_mhw',
  'elevation_sanders_creek_fl_mllw',
  'elevation_santa_ana_river_ca_mhw',
  'elevation_santa_barbara_ca_mhw',
  'elevation_santa_clara_river_ca_mhw',
  'elevation_santa_margarita_river_ca_mhw',
  'elevation_santa_maria_river_ca_mhw',
  'elevation_santa_monica_bay_ca_mhw',
  'elevation_santa_monica_ca_mhw',
  'elevation_santa_monica_ca_mhw',
  'elevation_santa_monica_canyon_ca_mhw',
  'elevation_santa_rosa_creek_ca_mhw',
  'elevation_santa_ynez_river_ca_mhw',
  'elevation_sapelo_sound_ga_navd88',
  'elevation_sarasota_bay_fl_mhw',
  'elevation_sassafrass_river_md_navd88',
  'elevation_savanna_river_ga_mhw',
  'elevation_savannah_ga_mhw',
  'elevation_schwan_lagoon_ca_mhw',
  'elevation_scott_creek_lagoon_ca_mhw',
  'elevation_seaside_creek_ca_mhw',
  'elevation_seaside_or_mhw',
  'elevation_seldovia_ak_mhw',
  'elevation_sengekontacket_pond_ma_mhw',
  'elevation_severn_river_md_navd88',
  'elevation_seward_ak_mhw',
  'elevation_sewee_bay_sc_navd88',
  'elevation_shallotte_river_nc_navd88',
  'elevation_shark_river_nj_mhw',
  'elevation_sheepscot_bay_me_mllw',
  'elevation_shemya_ak_mhw',
  'elevation_shinnecock_bay_ny_mhw',
  'elevation_shrewsbury_river_nj_mhw',
  'elevation_siletz_bay_or_mhw',
  'elevation_siltcoos_river_or_mhw',
  'elevation_sink_creek_fl_mllw',
  'elevation_sitka_ak_mhw',
  'elevation_sixes_river_or_mhw',
  'elevation_smith_lake_or_mhw',
  'elevation_smith_river_ca_mhw',
  'elevation_somes_sound_sw_bay_me_mhw',
  'elevation_south_edisto_river_sc_navd88',
  'elevation_south_padre_island_tx_mhw',
  'elevation_south_padre_island_tx_navd88',
  'elevation_southcentral_ak_mhw',
  'elevation_southeast_ak_mhhw',
  'elevation_southeast_atlantic_msl',
  'elevation_southern_ca_msl',
  'elevation_southern_cap_coastal_ponds_ma_mhw',
  'elevation_spring_warrior_creek_fl_mllw',
  'elevation_squashan_creek_ca_mhw',
  'elevation_squaw_creek_or_mhw',
  'elevation_steamboat_creek_wa_mhw',
  'elevation_stone_lagoon_ca_mhw',
  'elevation_storter_bay_fl_mllw',
  'elevation_strait_of_juan_de_fuca_wa_mllw',
  'elevation_stump_beach_ca_mhw',
  'elevation_suislaw_river_or_mhw',
  'elevation_sunset_bay_or_mhw',
  'elevation_sutton_creek_or_mhw',
  'elevation_suwanee_river_fl_mllw',
  'elevation_sycamore_creek_ca_mhw',
  'elevation_tahkenitch_creek_or_mhw',
  'elevation_tahola_wa_mhw',
  'elevation_tajiguas_creek_ca_mhw',
  'elevation_tampa_bay_1_fl_mhw',
  'elevation_tampa_bay_2_fl_mhw',
  'elevation_tampa_bay_3_fl_mhw',
  'elevation_tampa_bay_fl_mhw',
  'elevation_tampa_bay_fl_navd88',
  'elevation_tatitlek_ak_mhhw',
  'elevation_ten_mile_river_ca_mhw',
  'elevation_ten_thousand_islands_marco_bay_fl_navd88',
  'elevation_tenmile_creek_1_or_mhw',
  'elevation_tenmile_creek_2_or_mhw',
  'elevation_tennessee_valley_lagoon_ca_mhw',
  'elevation_thomas_creek_or_mhw',
  'elevation_threemile_creek_or_mhw',
  'elevation_tijuana_river_ca_mhw',
  'elevation_tillamook_bay_or_mhw',
  'elevation_tomales_bay_ca_mhw',
  'elevation_trancas_ca_mhw',
  'elevation_tunitas_creek_ca_mhw',
  'elevation_twomile_creek_2_or_mhw',
  'elevation_umpqua_river_or_mhw',
  'elevation_unalaska_ak_mhhw',
  'elevation_upper_laguna_madre_tx_mllw',
  'elevation_usal_creek_ca_mhw',
  'elevation_ventura_harbor_ca_mhw',
  'elevation_ventura_river_ca_mhw',
  'elevation_villa_creek_lagoon_ca_mhw',
  'elevation_virgin_creek_ca_mhw',
  'elevation_virginia_beach_va_mhw',
  'elevation_waccasassa_river_fl_mllw',
  'elevation_waddell_creek_ca_mhw',
  'elevation_wages_creek_ca_mhw',
  'elevation_weeki_wachee_bay_fl_mllw',
  'elevation_western_gulf_msl',
  'elevation_whackup_creek_fl_mllw',
  'elevation_whale_creek_wa_mhw',
  'elevation_whaleshead_creek_or_mhw',
  'elevation_whisky_run_or_mhw',
  'elevation_whitewater_bay_fl_mllw',
  'elevation_whittier_ak_mhhw',
  'elevation_wilder_creek_ca_mhw',
  'elevation_willapa_bay_wa_mhw',
  'elevation_wilson_creek_ca_mhw',
  'elevation_winchuck_river_or_mhw',
  'elevation_withlacoochee_bay_fl_mllw',
  'elevation_wreck_creek_wa_mhw',
  'elevation_yachats_river_or_mhw',
  'elevation_yakutat_ak_mhhw',
  'elevation_yaquina_bay_or_mhw',
  'elevation_yates_creek_fl_mllw',
  'elevation_younger_lagoon_ca_mhw',
  'elevation_zuma_canyon_ca_mhw',

  'atlantic.nadp.gridded.annual_wet_nh4',
  'atlantic.nadp.gridded.annual_wet_no3',
  'atlantic.nadp.gridded.annual_wet_totn',
  'atlantic.nadpcmaq.gridded.annual_dry_nh3',
  'atlantic.nadpcmaq.gridded.annual_dry_nh4',
  'atlantic.nadpcmaq.gridded.annual_dry_no3',
  'atlantic.nadpcmaq.gridded.annual_dry_totn',
  'atlantic.nadpcmaq.gridded.annual_dry_totnp',
  'atlantic.nadpcmaq.gridded.annual_dry_totun',
  'atlantic.nadpcmaq.gridded.annual_dry_totunp',
  'atlantic.nadpcmaq.gridded.annual_dry_net_nh3',
  'atlantic.nadpcmaq.gridded.annual_dry_oxidized_totn',
  'atlantic.nadpcmaq.gridded.annual_dry_oxidized_totnp',
  'atlantic.nadpcmaq.gridded.annual_dry_reduced_totn',
  'atlantic.nadpcmaq.gridded.annual_dry_reduced_totnp',
  'atlantic.nadpcmaq.gridded.annual_wet_nh4',
  'atlantic.nadpcmaq.gridded.annual_wet_no3',
  'atlantic.nadpcmaq.gridded.annual_wet_totn',
  'atlantic.nadpcmaq.gridded.annual_oxidized_totn',
  'atlantic.nadpcmaq.gridded.annual_oxidized_totnp',
  'atlantic.nadpcmaq.gridded.annual_reduced_totn',
  'atlantic.nadpcmaq.gridded.annual_reduced_totnp',
  'atlantic.nadpcmaq.gridded.annual_total_totn',

  'gulf.nadp.gridded.annual_wet_nh4',
  'gulf.nadp.gridded.annual_wet_no3',
  'gulf.nadp.gridded.annual_wet_totn',
  'gulf.nadpcmaq.gridded.annual_dry_nh3',
  'gulf.nadpcmaq.gridded.annual_dry_nh4',
  'gulf.nadpcmaq.gridded.annual_dry_no3',
  'gulf.nadpcmaq.gridded.annual_dry_totn',
  'gulf.nadpcmaq.gridded.annual_dry_totnp',
  'gulf.nadpcmaq.gridded.annual_dry_totun',
  'gulf.nadpcmaq.gridded.annual_dry_totunp',
  'gulf.nadpcmaq.gridded.annual_dry_net_nh3',
  'gulf.nadpcmaq.gridded.annual_dry_oxidized_totn',
  'gulf.nadpcmaq.gridded.annual_dry_oxidized_totnp',
  'gulf.nadpcmaq.gridded.annual_dry_reduced_totn',
  'gulf.nadpcmaq.gridded.annual_dry_reduced_totnp',
  'gulf.nadpcmaq.gridded.annual_wet_nh4',
  'gulf.nadpcmaq.gridded.annual_wet_no3',
  'gulf.nadpcmaq.gridded.annual_wet_totn',
  'gulf.nadpcmaq.gridded.annual_oxidized_totn',
  'gulf.nadpcmaq.gridded.annual_oxidized_totnp',
  'gulf.nadpcmaq.gridded.annual_reduced_totn',
  'gulf.nadpcmaq.gridded.annual_reduced_totnp',
  'gulf.nadpcmaq.gridded.annual_total_totn',
    
  'pacific.nadp.gridded.annual_wet_nh4',
  'pacific.nadp.gridded.annual_wet_no3',
  'pacific.nadp.gridded.annual_wet_totn',
  'pacific.nadpcmaq.gridded.annual_dry_nh3',
  'pacific.nadpcmaq.gridded.annual_dry_nh4',
  'pacific.nadpcmaq.gridded.annual_dry_no3',
  'pacific.nadpcmaq.gridded.annual_dry_totn',
  'pacific.nadpcmaq.gridded.annual_dry_totnp',
  'pacific.nadpcmaq.gridded.annual_dry_totun',
  'pacific.nadpcmaq.gridded.annual_dry_totunp',
  'pacific.nadpcmaq.gridded.annual_dry_net_nh3',
  'pacific.nadpcmaq.gridded.annual_dry_oxidized_totn',
  'pacific.nadpcmaq.gridded.annual_dry_oxidized_totnp',
  'pacific.nadpcmaq.gridded.annual_dry_reduced_totn',
  'pacific.nadpcmaq.gridded.annual_dry_reduced_totnp',
  'pacific.nadpcmaq.gridded.annual_wet_nh4',
  'pacific.nadpcmaq.gridded.annual_wet_no3',
  'pacific.nadpcmaq.gridded.annual_wet_totn',
  'pacific.nadpcmaq.gridded.annual_oxidized_totn',
  'pacific.nadpcmaq.gridded.annual_oxidized_totnp',
  'pacific.nadpcmaq.gridded.annual_reduced_totn',
  'pacific.nadpcmaq.gridded.annual_reduced_totnp',
  'pacific.nadpcmaq.gridded.annual_total_totn',

  'atlantic.annual_temperature_minimum',
  'gulf.annual_temperature_minimum',
  'pacific.annual_temperature_minimum',
  'atlantic.annual_temperature_maximum',
  'gulf.annual_temperature_maximum',
  'pacific.annual_temperature_maximum',

  'atlantic.monthly_temperature_minimum',
  'gulf.monthly_temperature_minimum',
  'pacific.monthly_temperature_minimum',
  'atlantic.monthly_temperature_maximum',
  'gulf.monthly_temperature_maximum',
  'pacific.monthly_temperature_maximum',

  'atlantic.annual_precipitation_mean',
  'gulf.annual_precipitation_mean',
  'pacific.annual_precipitation_mean',

  'atlantic.monthly_precipitation_mean',
  'gulf.monthly_precipitation_mean',
  'pacific.monthly_precipitation_mean',

  'atlantic.annual_precipitation',
  'gulf.annual_precipitation',
  'pacific.annual_precipitation',

  'atlantic.wetlands',
  'gulf.wetlands',
  'pacific.wetlands',
  'atlantic.seagrass',
  'gulf.seagrass',
  'pacific.seagrass',

  'atlantic.sediment_nca',
  'atlantic.sediment_nca_1990-2006',
  'atlantic.sediment_nca_2015',
  'atlantic.sediment_calculated',
  'atlantic.sediment_extracted',
  'atlantic.sediment_parsed',
  'atlantic.sediment_kriged',
  'gulf.sediment_nca',
  'gulf.sediment_nca_1990-2006',
  'gulf.sediment_nca_2015',
  'gulf.sediment_calculated',
  'gulf.sediment_extracted',
  'gulf.sediment_parsed',
  'gulf.sediment_gridded_mud',
  'gulf.sediment_gridded_mudu',
  'gulf.sediment_gridded_sand',
  'gulf.sediment_gridded_gravel',
  'gulf.sediment_gridded_rock',
  'pacific.sediment_nca',
  'pacific.sediment_nca_1990-2006',
  'pacific.sediment_nca_2015',
  'pacific.sediment_calculated',
  'pacific.sediment_extracted',
  'pacific.sediment_parsed',
  'lower_mississippi.sediment_nca_2015',

  'atlantic.soil',
  'gulf.soil',
  'pacific.soil',

  'atlantic.coastal_vulnerability',
  'gulf.coastal_vulnerability',
  'pacific.coastal_vulnerability',

  'atlantic.sensitivity_estuary',
  'gulf.sensitivity_estuary',
  'pacific.sensitivity_estuary',

  'atlantic.groundwater_contact_catchment',
  'gulf.groundwater_contact_catchment',
  'pacific.groundwater_contact_catchment',

  'atlantic.groundwater_contact_watershed',
  'gulf.groundwater_contact_watershed',
  'pacific.groundwater_contact_watershed',

  'atlantic.monthly.stream_discharge_nhd_line',
  'gulf.monthly.stream_discharge_nhd_line',
  'pacific.monthly.stream_discharge_nhd_line',
  'lower_mississippi.monthly.stream_discharge_nhd_line',

  'atlantic.monthly.stream_discharge_usgs_point',
  'gulf.monthly.stream_discharge_usgs_point',
  'pacific.monthly.stream_discharge_usgs_point',
  'lower_mississippi.monthly.stream_discharge_usgs_point',

  'atlantic.yearly.stream_discharge_nhd_line',
  'gulf.yearly.stream_discharge_nhd_line',
  'pacific.yearly.stream_discharge_nhd_line',
  'lower_mississippi.yearly.stream_discharge_nhd_line',

  'atlantic.yearly.stream_discharge_usgs_point',
  'gulf.yearly.stream_discharge_usgs_point',
  'pacific.yearly.stream_discharge_usgs_point',
  'lower_mississippi.yearly.stream_discharge_usgs_point',

  'atlantic.monthly.tide_point',
  'gulf.monthly.tide_point',
  'pacific.monthly.tide_point',
  'lower_mississippi.monthly.tide_point',

  'atlantic.yearly.tide_point',
  'gulf.yearly.tide_point',
  'pacific.yearly.tide_point',
  'lower_mississippi.yearly.tide_point',

  'atlantic.hourly.tide_current',
  'gulf.hourly.tide_current',
  'pacific.hourly.tide_current',
  'lower_mississippi.hourly.tide_current',

  'atlantic.monthly.tide_current',
  'gulf.monthly.tide_current',
  'pacific.monthly.tide_current',
  'lower_mississippi.monthly.tide_current',

  'atlantic.monthly.longshore_current',
  'gulf.monthly.longshore_current',
  'pacific.monthly.longshore_current',
  'lower_mississippi.monthly.longshore_current',

  'atlantic.yearly.estuary_flushing_ffm',
  'atlantic.yearly.estuary_flushing_tpm',
  'atlantic.monthly.estuary_flushing_ffm',
  'atlantic.monthly.estuary_flushing_tpm',
  'atlantic.monthly.estuary_flushing_froude',
  'atlantic.monthly.estuary_flushing_qv',

  'gulf.yearly.estuary_flushing_ffm',
  'gulf.yearly.estuary_flushing_tpm',
  'gulf.monthly.estuary_flushing_ffm',
  'gulf.monthly.estuary_flushing_tpm',
  'gulf.monthly.estuary_flushing_froude',
  'gulf.monthly.estuary_flushing_qv',


  'lower_mississippi.yearly.estuary_flushing_ffm',
  'lower_mississippi.yearly.estuary_flushing_tpm',
  'lower_mississippi.monthly.estuary_flushing_ffm',
  'lower_mississippi.monthly.estuary_flushing_tpm',
  'lower_mississippi.monthly.estuary_flushing_froude',
  'lower_mississippi.monthly.estuary_flushing_qv',

  'pacific.yearly.estuary_flushing_ffm',
  'pacific.yearly.estuary_flushing_tpm',
  'pacific.monthly.estuary_flushing_ffm',
  'pacific.monthly.estuary_flushing_tpm',
  'pacific.monthly.estuary_flushing_froude',
  'pacific.monthly.estuary_flushing_qv',

  'atlantic.salinity_point',
  'gulf.salinity_point',
  'lower_mississippi.salinity_point',
  'pacific.salinity_point',

  'atlantic.seagrass_point',
  'gulf.seagrass_point',
  'pacific.seagrass_point',

  'atlantic.seagrass_point2',
  'gulf.seagrass_point2',

  'atlantic.seagrass_polygon',
  'gulf.seagrass_polygon',
  'pacific.seagrass_polygon',

  'fl.seagrass_polygon',

  'eastfl.seagrass_polygon',
  'westfl.seagrass_polygon',
  'tampafl.seagrass_polygon',

  'atlantic.esat1_water_quality',
  'gulf.esat1_water_quality',
  'pacific.esat1_water_quality',

  'atlantic.esat2_water_quality',
  'gulf.esat2_water_quality',
  'pacific.esat2_water_quality',

  'flowlines_upper_ct_river',
  'flowlines_puget_sound_watershed',
  'discharge_points_federal_puget_sound_watershed',
  'discharge_points_state_puget_sound_watershed',

  'atlantic.nitrogen.estuary_cmaq', # Old
  'gulf.nitrogen.estuary_cmaq', # Old
  'pacific.nitrogen.estuary_cmaq', # Old

  'atlantic.dry_nh3.estuary_cmaq',
  'atlantic.dry_natural_nh3.estuary_cmaq',
  'atlantic.dry_nhx.estuary_cmaq',
  'atlantic.dry_so2.estuary_cmaq',
  'atlantic.dry_oxidized_nox.estuary_cmaq',
  'atlantic.dry_oxidized_organic_nitrogen.estuary_cmaq',
  'atlantic.dry_oxidized_no3.estuary_cmaq',
  'atlantic.dry_oxidized_nitrogen.estuary_cmaq',
  'atlantic.dry_reduced_nitrogen.estuary_cmaq',
  'atlantic.dry_sulfur.estuary_cmaq',
  'atlantic.chloride.estuary_cmaq',
  'atlantic.total_nitrogen.estuary_cmaq',
  'atlantic.oxidized_nitrogen.estuary_cmaq',
  'atlantic.reduced_nitrogen.estuary_cmaq',
  'atlantic.sulfur.estuary_cmaq',
  'atlantic.wet_nh3.estuary_cmaq',
  'atlantic.wet_chloride.estuary_cmaq',
  'atlantic.wet_so4.estuary_cmaq',
  'atlantic.wet_oxidized_no3.estuary_cmaq',
  'atlantic.wet_oxidized_nox.estuary_cmaq',
  'atlantic.wet_organic_nitrogen.estuary_cmaq',
  'atlantic.wet_oxidized_nitrogen.estuary_cmaq',
  'atlantic.wet_reduced_nitrogen.estuary_cmaq',
  'atlantic.wet_sulfur.estuary_cmaq',

  'gulf.dry_nh3.estuary_cmaq',
  'gulf.dry_natural_nh3.estuary_cmaq',
  'gulf.dry_nhx.estuary_cmaq',
  'gulf.dry_so2.estuary_cmaq',
  'gulf.dry_oxidized_nox.estuary_cmaq',
  'gulf.dry_oxidized_organic_nitrogen.estuary_cmaq',
  'gulf.dry_oxidized_no3.estuary_cmaq',
  'gulf.dry_oxidized_nitrogen.estuary_cmaq',
  'gulf.dry_reduced_nitrogen.estuary_cmaq',
  'gulf.dry_sulfur.estuary_cmaq',
  'gulf.chloride.estuary_cmaq',
  'gulf.total_nitrogen.estuary_cmaq',
  'gulf.oxidized_nitrogen.estuary_cmaq',
  'gulf.reduced_nitrogen.estuary_cmaq',
  'gulf.sulfur.estuary_cmaq',
  'gulf.wet_nh3.estuary_cmaq',
  'gulf.wet_chloride.estuary_cmaq',
  'gulf.wet_so4.estuary_cmaq',
  'gulf.wet_oxidized_no3.estuary_cmaq',
  'gulf.wet_oxidized_nox.estuary_cmaq',
  'gulf.wet_organic_nitrogen.estuary_cmaq',
  'gulf.wet_oxidized_nitrogen.estuary_cmaq',
  'gulf.wet_reduced_nitrogen.estuary_cmaq',
  'gulf.wet_sulfur.estuary_cmaq',

  'pacific.dry_nh3.estuary_cmaq',
  'pacific.dry_natural_nh3.estuary_cmaq',
  'pacific.dry_nhx.estuary_cmaq',
  'pacific.dry_so2.estuary_cmaq',
  'pacific.dry_oxidized_nox.estuary_cmaq',
  'pacific.dry_oxidized_organic_nitrogen.estuary_cmaq',
  'pacific.dry_oxidized_no3.estuary_cmaq',
  'pacific.dry_oxidized_nitrogen.estuary_cmaq',
  'pacific.dry_reduced_nitrogen.estuary_cmaq',
  'pacific.dry_sulfur.estuary_cmaq',
  'pacific.chloride.estuary_cmaq',
  'pacific.total_nitrogen.estuary_cmaq',
  'pacific.oxidized_nitrogen.estuary_cmaq',
  'pacific.reduced_nitrogen.estuary_cmaq',
  'pacific.sulfur.estuary_cmaq',
  'pacific.wet_nh3.estuary_cmaq',
  'pacific.wet_chloride.estuary_cmaq',
  'pacific.wet_so4.estuary_cmaq',
  'pacific.wet_oxidized_no3.estuary_cmaq',
  'pacific.wet_oxidized_nox.estuary_cmaq',
  'pacific.wet_organic_nitrogen.estuary_cmaq',
  'pacific.wet_oxidized_nitrogen.estuary_cmaq',
  'pacific.wet_reduced_nitrogen.estuary_cmaq',
  'pacific.wet_sulfur.estuary_cmaq',

  'atlantic.nitrogen.subestuary_cmaq', # Old
  'gulf.nitrogen.subestuary_cmaq', # Old
  'pacific.nitrogen.subestuary_cmaq', # Old

  'atlantic.dry_nh3.subestuary_cmaq',
  'atlantic.dry_natural_nh3.subestuary_cmaq',
  'atlantic.dry_nhx.subestuary_cmaq',
  'atlantic.dry_so2.subestuary_cmaq',
  'atlantic.dry_oxidized_nox.subestuary_cmaq',
  'atlantic.dry_oxidized_organic_nitrogen.subestuary_cmaq',
  'atlantic.dry_oxidized_no3.subestuary_cmaq',
  'atlantic.dry_oxidized_nitrogen.subestuary_cmaq',
  'atlantic.dry_reduced_nitrogen.subestuary_cmaq',
  'atlantic.dry_sulfur.subestuary_cmaq',
  'atlantic.chloride.subestuary_cmaq',
  'atlantic.total_nitrogen.subestuary_cmaq',
  'atlantic.oxidized_nitrogen.subestuary_cmaq',
  'atlantic.reduced_nitrogen.subestuary_cmaq',
  'atlantic.sulfur.subestuary_cmaq',
  'atlantic.wet_nh3.subestuary_cmaq',
  'atlantic.wet_chloride.subestuary_cmaq',
  'atlantic.wet_so4.subestuary_cmaq',
  'atlantic.wet_oxidized_no3.subestuary_cmaq',
  'atlantic.wet_oxidized_nox.subestuary_cmaq',
  'atlantic.wet_organic_nitrogen.subestuary_cmaq',
  'atlantic.wet_oxidized_nitrogen.subestuary_cmaq',
  'atlantic.wet_reduced_nitrogen.subestuary_cmaq',
  'atlantic.wet_sulfur.subestuary_cmaq',

  'gulf.dry_nh3.subestuary_cmaq',
  'gulf.dry_natural_nh3.subestuary_cmaq',
  'gulf.dry_nhx.subestuary_cmaq',
  'gulf.dry_so2.subestuary_cmaq',
  'gulf.dry_oxidized_nox.subestuary_cmaq',
  'gulf.dry_oxidized_organic_nitrogen.subestuary_cmaq',
  'gulf.dry_oxidized_no3.subestuary_cmaq',
  'gulf.dry_oxidized_nitrogen.subestuary_cmaq',
  'gulf.dry_reduced_nitrogen.subestuary_cmaq',
  'gulf.dry_sulfur.subestuary_cmaq',
  'gulf.chloride.subestuary_cmaq',
  'gulf.total_nitrogen.subestuary_cmaq',
  'gulf.oxidized_nitrogen.subestuary_cmaq',
  'gulf.reduced_nitrogen.subestuary_cmaq',
  'gulf.sulfur.subestuary_cmaq',
  'gulf.wet_nh3.subestuary_cmaq',
  'gulf.wet_chloride.subestuary_cmaq',
  'gulf.wet_so4.subestuary_cmaq',
  'gulf.wet_oxidized_no3.subestuary_cmaq',
  'gulf.wet_oxidized_nox.subestuary_cmaq',
  'gulf.wet_organic_nitrogen.subestuary_cmaq',
  'gulf.wet_oxidized_nitrogen.subestuary_cmaq',
  'gulf.wet_reduced_nitrogen.subestuary_cmaq',
  'gulf.wet_sulfur.subestuary_cmaq',

  'pacific.dry_nh3.subestuary_cmaq',
  'pacific.dry_natural_nh3.subestuary_cmaq',
  'pacific.dry_nhx.subestuary_cmaq',
  'pacific.dry_so2.subestuary_cmaq',
  'pacific.dry_oxidized_nox.subestuary_cmaq',
  'pacific.dry_oxidized_organic_nitrogen.subestuary_cmaq',
  'pacific.dry_oxidized_no3.subestuary_cmaq',
  'pacific.dry_oxidized_nitrogen.subestuary_cmaq',
  'pacific.dry_reduced_nitrogen.subestuary_cmaq',
  'pacific.dry_sulfur.subestuary_cmaq',
  'pacific.chloride.subestuary_cmaq',
  'pacific.total_nitrogen.subestuary_cmaq',
  'pacific.oxidized_nitrogen.subestuary_cmaq',
  'pacific.reduced_nitrogen.subestuary_cmaq',
  'pacific.sulfur.subestuary_cmaq',
  'pacific.wet_nh3.subestuary_cmaq',
  'pacific.wet_chloride.subestuary_cmaq',
  'pacific.wet_so4.subestuary_cmaq',
  'pacific.wet_oxidized_no3.subestuary_cmaq',
  'pacific.wet_oxidized_nox.subestuary_cmaq',
  'pacific.wet_organic_nitrogen.subestuary_cmaq',
  'pacific.wet_oxidized_nitrogen.subestuary_cmaq',
  'pacific.wet_reduced_nitrogen.subestuary_cmaq',
  'pacific.wet_sulfur.subestuary_cmaq',

  'atlantic.nitrogen.source_gridded.crops',
  'atlantic.nitrogen.source_gridded.fertilizer',
  'atlantic.nitrogen.source_gridded.manure',
  'gulf.nitrogen.source_gridded.crops',
  'gulf.nitrogen.source_gridded.fertilizer',
  'gulf.nitrogen.source_gridded.manure',
  'pacific.nitrogen.source_gridded.crops',
  'pacific.nitrogen.source_gridded.fertilizer',
  'pacific.nitrogen.source_gridded.manure',

  'atlantic.total_nitrogen.source_watershed_nadp',
  'gulf.total_nitrogen.source_watershed_nadp',
  'pacific.total_nitrogen.source_watershed_nadp',

  'atlantic.wet_nitrogen.source_watershed_nadp',
  'gulf.wet_nitrogen.source_watershed_nadp',
  'pacific.wet_nitrogen.source_watershed_nadp',

  'atlantic.dry_nitrogen.source_watershed_nadp',
  'gulf.dry_nitrogen.source_watershed_nadp',
  'pacific.dry_nitrogen.source_watershed_nadp',

  'atlantic.oxidized_nitrogen.source_watershed_nadp',
  'gulf.oxidized_nitrogen.source_watershed_nadp',
  'pacific.oxidized_nitrogen.source_watershed_nadp',

  'atlantic.reduced_nitrogen.source_watershed_nadp',
  'gulf.reduced_nitrogen.source_watershed_nadp',
  'pacific.reduced_nitrogen.source_watershed_nadp',

  'atlantic.dry_oxidized_nitrogen.source_watershed_nadp',
  'gulf.dry_oxidized_nitrogen.source_watershed_nadp',
  'pacific.dry_oxidized_nitrogen.source_watershed_nadp',

  'atlantic.dry_reduced_nitrogen.source_watershed_nadp',
  'gulf.dry_reduced_nitrogen.source_watershed_nadp',
  'pacific.dry_reduced_nitrogen.source_watershed_nadp',

  'atlantic.wet_no3.source_watershed_nadp',
  'gulf.wet_no3.source_watershed_nadp',
  'pacific.wet_no3.source_watershed_nadp',

  'atlantic.wet_nh4.source_watershed_nadp',
  'gulf.wet_nh4.source_watershed_nadp',
  'pacific.wet_nh4.source_watershed_nadp',

  'atlantic.dry_no3.source_watershed_nadp',
  'gulf.dry_no3.source_watershed_nadp',
  'pacific.dry_no3.source_watershed_nadp',

  'atlantic.dry_nh4.source_watershed_nadp',
  'gulf.dry_nh4.source_watershed_nadp',
  'pacific.dry_nh4.source_watershed_nadp',

  'atlantic.nitrogen.source_watershed_cmaq', # Old
  'gulf.nitrogen.source_watershed_cmaq', # Old
  'pacific.nitrogen.source_watershed_cmaq', # Old

  'atlantic.dry_nh3.source_watershed_cmaq',
  'atlantic.dry_natural_nh3.source_watershed_cmaq',
  'atlantic.dry_nhx.source_watershed_cmaq',
  'atlantic.dry_so2.source_watershed_cmaq',
  'atlantic.dry_oxidized_nox.source_watershed_cmaq',
  'atlantic.dry_oxidized_organic_nitrogen.source_watershed_cmaq',
  'atlantic.dry_oxidized_no3.source_watershed_cmaq',
  'atlantic.dry_oxidized_nitrogen.source_watershed_cmaq',
  'atlantic.dry_reduced_nitrogen.source_watershed_cmaq',
  'atlantic.dry_sulfur.source_watershed_cmaq',
  'atlantic.chloride.source_watershed_cmaq',
  'atlantic.total_nitrogen.source_watershed_cmaq',
  'atlantic.oxidized_nitrogen.source_watershed_cmaq',
  'atlantic.reduced_nitrogen.source_watershed_cmaq',
  'atlantic.sulfur.source_watershed_cmaq',
  'atlantic.wet_nh3.source_watershed_cmaq',
  'atlantic.wet_chloride.source_watershed_cmaq',
  'atlantic.wet_so4.source_watershed_cmaq',
  'atlantic.wet_oxidized_no3.source_watershed_cmaq',
  'atlantic.wet_oxidized_nox.source_watershed_cmaq',
  'atlantic.wet_organic_nitrogen.source_watershed_cmaq',
  'atlantic.wet_oxidized_nitrogen.source_watershed_cmaq',
  'atlantic.wet_reduced_nitrogen.source_watershed_cmaq',
  'atlantic.wet_sulfur.source_watershed_cmaq',

  'gulf.dry_nh3.source_watershed_cmaq',
  'gulf.dry_natural_nh3.source_watershed_cmaq',
  'gulf.dry_nhx.source_watershed_cmaq',
  'gulf.dry_so2.source_watershed_cmaq',
  'gulf.dry_oxidized_nox.source_watershed_cmaq',
  'gulf.dry_oxidized_organic_nitrogen.source_watershed_cmaq',
  'gulf.dry_oxidized_no3.source_watershed_cmaq',
  'gulf.dry_oxidized_nitrogen.source_watershed_cmaq',
  'gulf.dry_reduced_nitrogen.source_watershed_cmaq',
  'gulf.dry_sulfur.source_watershed_cmaq',
  'gulf.chloride.source_watershed_cmaq',
  'gulf.total_nitrogen.source_watershed_cmaq',
  'gulf.oxidized_nitrogen.source_watershed_cmaq',
  'gulf.reduced_nitrogen.source_watershed_cmaq',
  'gulf.sulfur.source_watershed_cmaq',
  'gulf.wet_nh3.source_watershed_cmaq',
  'gulf.wet_chloride.source_watershed_cmaq',
  'gulf.wet_so4.source_watershed_cmaq',
  'gulf.wet_oxidized_no3.source_watershed_cmaq',
  'gulf.wet_oxidized_nox.source_watershed_cmaq',
  'gulf.wet_organic_nitrogen.source_watershed_cmaq',
  'gulf.wet_oxidized_nitrogen.source_watershed_cmaq',
  'gulf.wet_reduced_nitrogen.source_watershed_cmaq',
  'gulf.wet_sulfur.source_watershed_cmaq',

  'pacific.dry_nh3.source_watershed_cmaq',
  'pacific.dry_natural_nh3.source_watershed_cmaq',
  'pacific.dry_nhx.source_watershed_cmaq',
  'pacific.dry_so2.source_watershed_cmaq',
  'pacific.dry_oxidized_nox.source_watershed_cmaq',
  'pacific.dry_oxidized_organic_nitrogen.source_watershed_cmaq',
  'pacific.dry_oxidized_no3.source_watershed_cmaq',
  'pacific.dry_oxidized_nitrogen.source_watershed_cmaq',
  'pacific.dry_reduced_nitrogen.source_watershed_cmaq',
  'pacific.dry_sulfur.source_watershed_cmaq',
  'pacific.chloride.source_watershed_cmaq',
  'pacific.total_nitrogen.source_watershed_cmaq',
  'pacific.oxidized_nitrogen.source_watershed_cmaq',
  'pacific.reduced_nitrogen.source_watershed_cmaq',
  'pacific.sulfur.source_watershed_cmaq',
  'pacific.wet_nh3.source_watershed_cmaq',
  'pacific.wet_chloride.source_watershed_cmaq',
  'pacific.wet_so4.source_watershed_cmaq',
  'pacific.wet_oxidized_no3.source_watershed_cmaq',
  'pacific.wet_oxidized_nox.source_watershed_cmaq',
  'pacific.wet_organic_nitrogen.source_watershed_cmaq',
  'pacific.wet_oxidized_nitrogen.source_watershed_cmaq',
  'pacific.wet_reduced_nitrogen.source_watershed_cmaq',
  'pacific.wet_sulfur.source_watershed_cmaq',

  'atlantic.nitrogen.source_watershed_point',
  'gulf.nitrogen.source_watershed_point',
  'pacific.nitrogen.source_watershed_point',

  'atlantic.nitrogen.source_watershed_nonpoint',
  'gulf.nitrogen.source_watershed_nonpoint',
  'pacific.nitrogen.source_watershed_nonpoint',

  'atlantic.nitrogen.load_estuary_sparrow_1992',
  'gulf.nitrogen.load_estuary_sparrow_1992',
  'pacific.nitrogen.load_estuary_sparrow_1992',

  'atlantic.nitrogen.load_estuary_sparrow_2002_mrb1',
  'atlantic.nitrogen.load_estuary_sparrow_2002_mrb2',
  'atlantic.nitrogen.load_estuary_sparrow_2011_mrb1mrb2',
  'gulf.nitrogen.load_estuary_sparrow_2002_mrb2',
  'gulf.nitrogen.load_estuary_sparrow_2002_mrb5',
  'gulf.nitrogen.load_estuary_sparrow_2011_mrb5',
  'pacific.nitrogen.load_estuary_sparrow_2002_mrb7',
  'pacific.nitrogen.load_estuary_sparrow_2011_mrb7',
  'pacific.nitrogen.load_estuary_sparrow_2011_mrb8',

  'atlantic.nitrogen.load_subestuary_sparrow_1992',
  'gulf.nitrogen.load_subestuary_sparrow_1992',
  'pacific.nitrogen.load_subestuary_sparrow_1992',

  'atlantic.nitrogen.load_subestuary_sparrow_2002_mrb1',
  'atlantic.nitrogen.load_subestuary_sparrow_2002_mrb2',
  'gulf.nitrogen.load_subestuary_sparrow_2002_mrb2',
  'gulf.nitrogen.load_subestuary_sparrow_2002_mrb5',
  'pacific.nitrogen.load_subestuary_sparrow_2002_mrb7',

  'atlantic.phosphorus.source_watershed_point',
  'gulf.phosphorus.source_watershed_point',
  'pacific.phosphorus.source_watershed_point',

  'atlantic.population_iclus1',
  'gulf.population_iclus1',
  'pacific.population_iclus1',

  'atlantic.population_iclus2',
  'gulf.population_iclus2',
  'pacific.population_iclus2',

  'gi_bmp_installations_new_england',
  'impervious_nhdplus_new_england',
  'impervious_ball_mountain_vt',
  'impervious_black_river_vt',
  'impervious_bloody_brook_vt',
  'impervious_branch_brook_vt',
  'impervious_calendar_brook_vt',
  'impervious_east_creek_vt',
  'impervious_jay_branch_vt',
  'impervious_jewett_brook_vt',
  'impervious_lulls_brook_vt',
  'impervious_mallets_bay_vt',
  'impervious_mallets_bay10_vt',
  'impervious_mallets_creek_vt',
  'impervious_mill_brook_vt',
  'impervious_mill_river_vt',
  'impervious_moon_brook_vt',
  'impervious_muddy_brook_vt',
  'impervious_north_branch_deerfield_river_vt',
  'impervious_ottauquechee_river_vt',
  'impervious_ottauquechee_river_headwaters_vt',
  'impervious_shelburne_river_vt',
  'impervious_snipe_island_vt',
  'impervious_winooski_river_vt',

  'meris_weekly_cyanobacteria_al',
  'meris_weekly_cyanobacteria_ar',
  'meris_weekly_cyanobacteria_az',
  'meris_weekly_cyanobacteria_ca',
  'meris_weekly_cyanobacteria_co',
  'meris_weekly_cyanobacteria_ct',
  'meris_weekly_cyanobacteria_de',
  'meris_weekly_cyanobacteria_fl',
  'meris_weekly_cyanobacteria_ga',
  'meris_weekly_cyanobacteria_ia',
  'meris_weekly_cyanobacteria_id',
  'meris_weekly_cyanobacteria_il',
  'meris_weekly_cyanobacteria_in',
  'meris_weekly_cyanobacteria_ks',
  'meris_weekly_cyanobacteria_ky',
  'meris_weekly_cyanobacteria_la',
  'meris_weekly_cyanobacteria_ma',
  'meris_weekly_cyanobacteria_md',
  'meris_weekly_cyanobacteria_me',
  'meris_weekly_cyanobacteria_mi',
  'meris_weekly_cyanobacteria_mn',
  'meris_weekly_cyanobacteria_mo',
  'meris_weekly_cyanobacteria_ms',
  'meris_weekly_cyanobacteria_mt',
  'meris_weekly_cyanobacteria_nc',
  'meris_weekly_cyanobacteria_nd',
  'meris_weekly_cyanobacteria_ne',
  'meris_weekly_cyanobacteria_nh',
  'meris_weekly_cyanobacteria_nj',
  'meris_weekly_cyanobacteria_nm',
  'meris_weekly_cyanobacteria_nv',
  'meris_weekly_cyanobacteria_ny',
  'meris_weekly_cyanobacteria_oh',
  'meris_weekly_cyanobacteria_ok',
  'meris_weekly_cyanobacteria_or',
  'meris_weekly_cyanobacteria_pa',
  'meris_weekly_cyanobacteria_ri',
  'meris_weekly_cyanobacteria_sc',
  'meris_weekly_cyanobacteria_sd',
  'meris_weekly_cyanobacteria_tn',
  'meris_weekly_cyanobacteria_tx',
  'meris_weekly_cyanobacteria_ut',
  'meris_weekly_cyanobacteria_va',
  'meris_weekly_cyanobacteria_vt',
  'meris_weekly_cyanobacteria_wa',
  'meris_weekly_cyanobacteria_wi',
  'meris_weekly_cyanobacteria_wv',
  'meris_weekly_cyanobacteria_wy',

  'sentinel_weekly_cyanobacteria_al',
  'sentinel_weekly_cyanobacteria_ar',
  'sentinel_weekly_cyanobacteria_az',
  'sentinel_weekly_cyanobacteria_ca',
  'sentinel_weekly_cyanobacteria_co',
  'sentinel_weekly_cyanobacteria_ct',
  'sentinel_weekly_cyanobacteria_de',
  'sentinel_weekly_cyanobacteria_fl',
  'sentinel_weekly_cyanobacteria_ga',
  'sentinel_weekly_cyanobacteria_ia',
  'sentinel_weekly_cyanobacteria_id',
  'sentinel_weekly_cyanobacteria_il',
  'sentinel_weekly_cyanobacteria_in',
  'sentinel_weekly_cyanobacteria_ks',
  'sentinel_weekly_cyanobacteria_ky',
  'sentinel_weekly_cyanobacteria_la',
  'sentinel_weekly_cyanobacteria_ma',
  'sentinel_weekly_cyanobacteria_md',
  'sentinel_weekly_cyanobacteria_me',
  'sentinel_weekly_cyanobacteria_mi',
  'sentinel_weekly_cyanobacteria_mn',
  'sentinel_weekly_cyanobacteria_mo',
  'sentinel_weekly_cyanobacteria_ms',
  'sentinel_weekly_cyanobacteria_mt',
  'sentinel_weekly_cyanobacteria_nc',
  'sentinel_weekly_cyanobacteria_nd',
  'sentinel_weekly_cyanobacteria_ne',
  'sentinel_weekly_cyanobacteria_nh',
  'sentinel_weekly_cyanobacteria_nj',
  'sentinel_weekly_cyanobacteria_nm',
  'sentinel_weekly_cyanobacteria_nv',
  'sentinel_weekly_cyanobacteria_ny',
  'sentinel_weekly_cyanobacteria_oh',
  'sentinel_weekly_cyanobacteria_ok',
  'sentinel_weekly_cyanobacteria_or',
  'sentinel_weekly_cyanobacteria_pa',
  'sentinel_weekly_cyanobacteria_ri',
  'sentinel_weekly_cyanobacteria_sc',
  'sentinel_weekly_cyanobacteria_sd',
  'sentinel_weekly_cyanobacteria_tn',
  'sentinel_weekly_cyanobacteria_tx',
  'sentinel_weekly_cyanobacteria_ut',
  'sentinel_weekly_cyanobacteria_va',
  'sentinel_weekly_cyanobacteria_vt',
  'sentinel_weekly_cyanobacteria_wa',
  'sentinel_weekly_cyanobacteria_wi',
  'sentinel_weekly_cyanobacteria_wv',
  'sentinel_weekly_cyanobacteria_wy',

  'land_change_2020_current_chesapeake',
  'land_change_2030_current_chesapeake',
  'land_change_2040_current_chesapeake',
  'land_change_2020_deregulation_chesapeake',
  'land_change_2030_deregulation_chesapeake',
  'land_change_2040_deregulation_chesapeake',
  'land_change_2020_redevelopment_chesapeake',
  'land_change_2030_redevelopment_chesapeake',
  'land_change_2040_redevelopment_chesapeake',
  'land_change_2020_conservation_chesapeake',
  'land_change_2030_conservation_chesapeake',
  'land_change_2040_conservation_chesapeake',

  'stream_temperature_point_08_lower_columbia_river',
  'stream_temperature_point_08_middle_columbia_river',
  'stream_temperature_point_median_07_new_england',
  'stream_temperature_point_median_08_new_england',
  'stream_temperature_point_mean_daily_range_07_new_england',
  'stream_temperature_point_mean_daily_range_08_new_england',
  'stream_temperature_point_maximum_daily_increase_07_new_england',
  'stream_temperature_point_maximum_daily_increase_08_new_england',
  'stream_temperature_point_maximum_daily_decrease_07_new_england',
  'stream_temperature_point_maximum_daily_decrease_08_new_england',
  'stream_temperature_point_maximum_new_england',
  'stream_temperature_point_day_of_maximum_new_england',
  'stream_temperature_line_08_lower_columbia_river',
  'stream_temperature_line_08_middle_columbia_river',
  'stream_temperature_line_current_08_upper_rogue_river',
  'stream_temperature_line_future_08_upper_rogue_river',
  'stream_temperature_line_07_new_england',
  'stream_temperature_line_meduxnekeag_river',
  'stream_shade_line_08_lower_columbia_river',
  'stream_shade_line_08_middle_columbia_river',

  'hspf_blackstone10_huc10_point.hydrology',
  'hspf_blackstone10_huc10_polygon.hydrology',
  'hspf_blackstone10_huc12_point.hydrology',
  'hspf_blackstone10_huc12_polygon.hydrology',
  'hspf_blackstone1_huc10_point.hydrology',
  'hspf_blackstone1_huc10_polygon.hydrology',
  'hspf_blackstone1_huc12_point.hydrology',
  'hspf_blackstone1_huc12_polygon.hydrology',
  'hspf_blackstone2_huc10_point.hydrology',
  'hspf_blackstone2_huc10_polygon.hydrology',
  'hspf_blackstone2_huc12_point.hydrology',
  'hspf_blackstone2_huc12_polygon.hydrology',
  'hspf_blackstone34_huc10_point.hydrology',
  'hspf_blackstone34_huc10_polygon.hydrology',
  'hspf_blackstone34_huc12_point.hydrology',
  'hspf_blackstone34_huc12_polygon.hydrology',
  'hspf_blackstone35_huc10_point.hydrology',
  'hspf_blackstone35_huc10_polygon.hydrology',
  'hspf_blackstone35_huc12_point.hydrology',
  'hspf_blackstone35_huc12_polygon.hydrology',
  'hspf_blackstone36_huc10_point.hydrology',
  'hspf_blackstone36_huc10_polygon.hydrology',
  'hspf_blackstone36_huc12_point.hydrology',
  'hspf_blackstone36_huc12_polygon.hydrology',
  'hspf_blackstone37_huc10_point.hydrology',
  'hspf_blackstone37_huc10_polygon.hydrology',
  'hspf_blackstone37_huc12_point.hydrology',
  'hspf_blackstone37_huc12_polygon.hydrology',
  'hspf_blackstone38_huc10_point.hydrology',
  'hspf_blackstone38_huc10_polygon.hydrology',
  'hspf_blackstone38_huc12_point.hydrology',
  'hspf_blackstone38_huc12_polygon.hydrology',
  'hspf_blackstone39_huc10_point.hydrology',
  'hspf_blackstone39_huc10_polygon.hydrology',
  'hspf_blackstone39_huc12_point.hydrology',
  'hspf_blackstone39_huc12_polygon.hydrology',
  'hspf_blackstone40_huc10_point.hydrology',
  'hspf_blackstone40_huc10_polygon.hydrology',
  'hspf_blackstone40_huc12_point.hydrology',
  'hspf_blackstone40_huc12_polygon.hydrology',
  'hspf_blackstone41_huc10_point.hydrology',
  'hspf_blackstone41_huc10_polygon.hydrology',
  'hspf_blackstone41_huc12_point.hydrology',
  'hspf_blackstone41_huc12_polygon.hydrology',
  'hspf_blackstone42_huc10_point.hydrology',
  'hspf_blackstone42_huc10_polygon.hydrology',
  'hspf_blackstone42_huc12_point.hydrology',
  'hspf_blackstone42_huc12_polygon.hydrology',
  'hspf_blackstone43_huc10_point.hydrology',
  'hspf_blackstone43_huc10_polygon.hydrology',
  'hspf_blackstone43_huc12_point.hydrology',
  'hspf_blackstone43_huc12_polygon.hydrology',
  'hspf_blackstone9_huc10_point.hydrology',
  'hspf_blackstone9_huc10_polygon.hydrology',
  'hspf_blackstone9_huc12_point.hydrology',
  'hspf_blackstone9_huc12_polygon.hydrology',
  'hspf_charles11_huc10_point.hydrology',
  'hspf_charles11_huc10_point.recharge_adjustment',
  'hspf_charles11_huc10_polygon.hydrology',
  'hspf_charles11_huc10_polygon.recharge_adjustment',
  'hspf_charles11_huc12_point.hydrology',
  'hspf_charles11_huc12_point.recharge_adjustment',
  'hspf_charles11_huc12_polygon.hydrology',
  'hspf_charles11_huc12_polygon.recharge_adjustment',
  'hspf_charles3_huc10_point.hydrology',
  'hspf_charles3_huc10_point.recharge_adjustment',
  'hspf_charles3_huc10_polygon.hydrology',
  'hspf_charles3_huc10_polygon.recharge_adjustment',
  'hspf_charles3_huc12_point.hydrology',
  'hspf_charles3_huc12_point.recharge_adjustment',
  'hspf_charles3_huc12_polygon.hydrology',
  'hspf_charles3_huc12_polygon.recharge_adjustment',
  'hspf_charles44_huc10_point.hydrology',
  'hspf_charles44_huc10_point.recharge_adjustment',
  'hspf_charles44_huc10_polygon.hydrology',
  'hspf_charles44_huc10_polygon.recharge_adjustment',
  'hspf_charles44_huc12_point.hydrology',
  'hspf_charles44_huc12_point.recharge_adjustment',
  'hspf_charles44_huc12_polygon.hydrology',
  'hspf_charles44_huc12_polygon.recharge_adjustment',
  'hspf_charles45_huc10_point.hydrology',
  'hspf_charles45_huc10_point.recharge_adjustment',
  'hspf_charles45_huc10_polygon.hydrology',
  'hspf_charles45_huc10_polygon.recharge_adjustment',
  'hspf_charles45_huc12_point.hydrology',
  'hspf_charles45_huc12_point.recharge_adjustment',
  'hspf_charles45_huc12_polygon.hydrology',
  'hspf_charles45_huc12_polygon.recharge_adjustment',
  'hspf_charles46_huc10_point.hydrology',
  'hspf_charles46_huc10_point.recharge_adjustment',
  'hspf_charles46_huc10_polygon.hydrology',
  'hspf_charles46_huc10_polygon.recharge_adjustment',
  'hspf_charles46_huc12_point.hydrology',
  'hspf_charles46_huc12_point.recharge_adjustment',
  'hspf_charles46_huc12_polygon.hydrology',
  'hspf_charles46_huc12_polygon.recharge_adjustment',
  'hspf_charles47_huc10_point.hydrology',
  'hspf_charles47_huc10_point.recharge_adjustment',
  'hspf_charles47_huc10_polygon.hydrology',
  'hspf_charles47_huc10_polygon.recharge_adjustment',
  'hspf_charles47_huc12_point.hydrology',
  'hspf_charles47_huc12_point.recharge_adjustment',
  'hspf_charles47_huc12_polygon.hydrology',
  'hspf_charles47_huc12_polygon.recharge_adjustment',
  'hspf_charles48_huc10_point.hydrology',
  'hspf_charles48_huc10_point.recharge_adjustment',
  'hspf_charles48_huc10_polygon.hydrology',
  'hspf_charles48_huc10_polygon.recharge_adjustment',
  'hspf_charles48_huc12_point.hydrology',
  'hspf_charles48_huc12_point.recharge_adjustment',
  'hspf_charles48_huc12_polygon.hydrology',
  'hspf_charles48_huc12_polygon.recharge_adjustment',
  'hspf_farmington12_huc10_point.hydrology',
  'hspf_farmington12_huc10_point.recharge_adjustment',
  'hspf_farmington12_huc10_polygon.hydrology',
  'hspf_farmington12_huc10_polygon.recharge_adjustment',
  'hspf_farmington12_huc12_point.hydrology',
  'hspf_farmington12_huc12_point.recharge_adjustment',
  'hspf_farmington12_huc12_polygon.hydrology',
  'hspf_farmington12_huc12_polygon.recharge_adjustment',
  'hspf_farmington49_huc10_point.hydrology',
  'hspf_farmington49_huc10_point.recharge_adjustment',
  'hspf_farmington49_huc10_polygon.hydrology',
  'hspf_farmington49_huc10_polygon.recharge_adjustment',
  'hspf_farmington49_huc12_point.hydrology',
  'hspf_farmington49_huc12_point.recharge_adjustment',
  'hspf_farmington49_huc12_polygon.hydrology',
  'hspf_farmington49_huc12_polygon.recharge_adjustment',
  'hspf_farmington4_huc10_point.hydrology',
  'hspf_farmington4_huc10_point.recharge_adjustment',
  'hspf_farmington4_huc10_polygon.hydrology',
  'hspf_farmington4_huc10_polygon.recharge_adjustment',
  'hspf_farmington4_huc12_point.hydrology',
  'hspf_farmington4_huc12_point.recharge_adjustment',
  'hspf_farmington4_huc12_polygon.hydrology',
  'hspf_farmington4_huc12_polygon.recharge_adjustment',
  'hspf_farmington50_huc10_point.hydrology',
  'hspf_farmington50_huc10_point.recharge_adjustment',
  'hspf_farmington50_huc10_polygon.hydrology',
  'hspf_farmington50_huc10_polygon.recharge_adjustment',
  'hspf_farmington50_huc12_point.hydrology',
  'hspf_farmington50_huc12_point.recharge_adjustment',
  'hspf_farmington50_huc12_polygon.hydrology',
  'hspf_farmington50_huc12_polygon.recharge_adjustment',
  'hspf_farmington51_huc10_point.hydrology',
  'hspf_farmington51_huc10_point.recharge_adjustment',
  'hspf_farmington51_huc10_polygon.hydrology',
  'hspf_farmington51_huc10_polygon.recharge_adjustment',
  'hspf_farmington51_huc12_point.hydrology',
  'hspf_farmington51_huc12_point.recharge_adjustment',
  'hspf_farmington51_huc12_polygon.hydrology',
  'hspf_farmington51_huc12_polygon.recharge_adjustment',
  'hspf_farmington52_huc10_point.hydrology',
  'hspf_farmington52_huc10_point.recharge_adjustment',
  'hspf_farmington52_huc10_polygon.hydrology',
  'hspf_farmington52_huc10_polygon.recharge_adjustment',
  'hspf_farmington52_huc12_point.hydrology',
  'hspf_farmington52_huc12_point.recharge_adjustment',
  'hspf_farmington52_huc12_polygon.hydrology',
  'hspf_farmington52_huc12_polygon.recharge_adjustment',
  'hspf_farmington53_huc10_point.hydrology',
  'hspf_farmington53_huc10_point.recharge_adjustment',
  'hspf_farmington53_huc10_polygon.hydrology',
  'hspf_farmington53_huc10_polygon.recharge_adjustment',
  'hspf_farmington53_huc12_point.hydrology',
  'hspf_farmington53_huc12_point.recharge_adjustment',
  'hspf_farmington53_huc12_polygon.hydrology',
  'hspf_farmington53_huc12_polygon.recharge_adjustment',
  'hspf_ipswich13_huc10_point.hydrology',
  'hspf_ipswich13_huc10_point.recharge_adjustment',
  'hspf_ipswich13_huc10_polygon.hydrology',
  'hspf_ipswich13_huc10_polygon.recharge_adjustment',
  'hspf_ipswich13_huc12_point.hydrology',
  'hspf_ipswich13_huc12_point.recharge_adjustment',
  'hspf_ipswich13_huc12_polygon.hydrology',
  'hspf_ipswich13_huc12_polygon.recharge_adjustment',
  'hspf_ipswich54_huc10_point.hydrology',
  'hspf_ipswich54_huc10_point.recharge_adjustment',
  'hspf_ipswich54_huc10_polygon.hydrology',
  'hspf_ipswich54_huc10_polygon.recharge_adjustment',
  'hspf_ipswich54_huc12_point.hydrology',
  'hspf_ipswich54_huc12_point.recharge_adjustment',
  'hspf_ipswich54_huc12_polygon.hydrology',
  'hspf_ipswich54_huc12_polygon.recharge_adjustment',
  'hspf_ipswich55_huc10_point.hydrology',
  'hspf_ipswich55_huc10_point.recharge_adjustment',
  'hspf_ipswich55_huc10_polygon.hydrology',
  'hspf_ipswich55_huc10_polygon.recharge_adjustment',
  'hspf_ipswich55_huc12_point.hydrology',
  'hspf_ipswich55_huc12_point.recharge_adjustment',
  'hspf_ipswich55_huc12_polygon.hydrology',
  'hspf_ipswich55_huc12_polygon.recharge_adjustment',
  'hspf_ipswich56_huc10_point.hydrology',
  'hspf_ipswich56_huc10_point.recharge_adjustment',
  'hspf_ipswich56_huc10_polygon.hydrology',
  'hspf_ipswich56_huc10_polygon.recharge_adjustment',
  'hspf_ipswich56_huc12_point.hydrology',
  'hspf_ipswich56_huc12_point.recharge_adjustment',
  'hspf_ipswich56_huc12_polygon.hydrology',
  'hspf_ipswich56_huc12_polygon.recharge_adjustment',
  'hspf_ipswich57_huc10_point.hydrology',
  'hspf_ipswich57_huc10_point.recharge_adjustment',
  'hspf_ipswich57_huc10_polygon.hydrology',
  'hspf_ipswich57_huc10_polygon.recharge_adjustment',
  'hspf_ipswich57_huc12_point.hydrology',
  'hspf_ipswich57_huc12_point.recharge_adjustment',
  'hspf_ipswich57_huc12_polygon.hydrology',
  'hspf_ipswich57_huc12_polygon.recharge_adjustment',
  'hspf_ipswich58_huc10_point.hydrology',
  'hspf_ipswich58_huc10_point.recharge_adjustment',
  'hspf_ipswich58_huc10_polygon.hydrology',
  'hspf_ipswich58_huc10_polygon.recharge_adjustment',
  'hspf_ipswich58_huc12_point.hydrology',
  'hspf_ipswich58_huc12_point.recharge_adjustment',
  'hspf_ipswich58_huc12_polygon.hydrology',
  'hspf_ipswich58_huc12_polygon.recharge_adjustment',
  'hspf_ipswich5_huc10_point.hydrology',
  'hspf_ipswich5_huc10_point.recharge_adjustment',
  'hspf_ipswich5_huc10_polygon.hydrology',
  'hspf_ipswich5_huc10_polygon.recharge_adjustment',
  'hspf_ipswich5_huc12_point.hydrology',
  'hspf_ipswich5_huc12_point.recharge_adjustment',
  'hspf_ipswich5_huc12_polygon.hydrology',
  'hspf_ipswich5_huc12_polygon.recharge_adjustment',
  'hspf_pawcatuck14_huc10_point.hydrology',
  'hspf_pawcatuck14_huc10_point.recharge_adjustment',
  'hspf_pawcatuck14_huc10_polygon.hydrology',
  'hspf_pawcatuck14_huc10_polygon.recharge_adjustment',
  'hspf_pawcatuck14_huc12_point.hydrology',
  'hspf_pawcatuck14_huc12_point.recharge_adjustment',
  'hspf_pawcatuck14_huc12_polygon.hydrology',
  'hspf_pawcatuck14_huc12_polygon.recharge_adjustment',
  'hspf_pawcatuck59_huc10_point.hydrology',
  'hspf_pawcatuck59_huc10_point.recharge_adjustment',
  'hspf_pawcatuck59_huc10_polygon.hydrology',
  'hspf_pawcatuck59_huc10_polygon.recharge_adjustment',
  'hspf_pawcatuck59_huc12_point.hydrology',
  'hspf_pawcatuck59_huc12_point.recharge_adjustment',
  'hspf_pawcatuck59_huc12_polygon.hydrology',
  'hspf_pawcatuck59_huc12_polygon.recharge_adjustment',
  'hspf_pawcatuck60_huc10_point.hydrology',
  'hspf_pawcatuck60_huc10_point.recharge_adjustment',
  'hspf_pawcatuck60_huc10_polygon.hydrology',
  'hspf_pawcatuck60_huc10_polygon.recharge_adjustment',
  'hspf_pawcatuck60_huc12_point.hydrology',
  'hspf_pawcatuck60_huc12_point.recharge_adjustment',
  'hspf_pawcatuck60_huc12_polygon.hydrology',
  'hspf_pawcatuck60_huc12_polygon.recharge_adjustment',
  'hspf_pawcatuck61_huc10_point.hydrology',
  'hspf_pawcatuck61_huc10_point.recharge_adjustment',
  'hspf_pawcatuck61_huc10_polygon.hydrology',
  'hspf_pawcatuck61_huc10_polygon.recharge_adjustment',
  'hspf_pawcatuck61_huc12_point.hydrology',
  'hspf_pawcatuck61_huc12_point.recharge_adjustment',
  'hspf_pawcatuck61_huc12_polygon.hydrology',
  'hspf_pawcatuck61_huc12_polygon.recharge_adjustment',
  'hspf_pawcatuck62_huc10_point.hydrology',
  'hspf_pawcatuck62_huc10_point.recharge_adjustment',
  'hspf_pawcatuck62_huc10_polygon.hydrology',
  'hspf_pawcatuck62_huc10_polygon.recharge_adjustment',
  'hspf_pawcatuck62_huc12_point.hydrology',
  'hspf_pawcatuck62_huc12_point.recharge_adjustment',
  'hspf_pawcatuck62_huc12_polygon.hydrology',
  'hspf_pawcatuck62_huc12_polygon.recharge_adjustment',
  'hspf_pawcatuck63_huc10_point.hydrology',
  'hspf_pawcatuck63_huc10_point.recharge_adjustment',
  'hspf_pawcatuck63_huc10_polygon.hydrology',
  'hspf_pawcatuck63_huc10_polygon.recharge_adjustment',
  'hspf_pawcatuck63_huc12_point.hydrology',
  'hspf_pawcatuck63_huc12_point.recharge_adjustment',
  'hspf_pawcatuck63_huc12_polygon.hydrology',
  'hspf_pawcatuck63_huc12_polygon.recharge_adjustment',
  'hspf_pawcatuck6_huc10_point.hydrology',
  'hspf_pawcatuck6_huc10_point.recharge_adjustment',
  'hspf_pawcatuck6_huc10_polygon.hydrology',
  'hspf_pawcatuck6_huc10_polygon.recharge_adjustment',
  'hspf_pawcatuck6_huc12_point.hydrology',
  'hspf_pawcatuck6_huc12_point.recharge_adjustment',
  'hspf_pawcatuck6_huc12_polygon.hydrology',
  'hspf_pawcatuck6_huc12_polygon.recharge_adjustment',
  'hspf_sudbury15_huc10_point.hydrology',
  'hspf_sudbury15_huc10_point.recharge_adjustment',
  'hspf_sudbury15_huc10_polygon.hydrology',
  'hspf_sudbury15_huc10_polygon.recharge_adjustment',
  'hspf_sudbury15_huc12_point.hydrology',
  'hspf_sudbury15_huc12_point.recharge_adjustment',
  'hspf_sudbury15_huc12_polygon.hydrology',
  'hspf_sudbury15_huc12_polygon.recharge_adjustment',
  'hspf_sudbury64_huc10_point.hydrology',
  'hspf_sudbury64_huc10_point.recharge_adjustment',
  'hspf_sudbury64_huc10_polygon.hydrology',
  'hspf_sudbury64_huc10_polygon.recharge_adjustment',
  'hspf_sudbury64_huc12_point.hydrology',
  'hspf_sudbury64_huc12_point.recharge_adjustment',
  'hspf_sudbury64_huc12_polygon.hydrology',
  'hspf_sudbury64_huc12_polygon.recharge_adjustment',
  'hspf_sudbury65_huc10_point.hydrology',
  'hspf_sudbury65_huc10_point.recharge_adjustment',
  'hspf_sudbury65_huc10_polygon.hydrology',
  'hspf_sudbury65_huc10_polygon.recharge_adjustment',
  'hspf_sudbury65_huc12_point.hydrology',
  'hspf_sudbury65_huc12_point.recharge_adjustment',
  'hspf_sudbury65_huc12_polygon.hydrology',
  'hspf_sudbury65_huc12_polygon.recharge_adjustment',
  'hspf_sudbury66_huc10_point.hydrology',
  'hspf_sudbury66_huc10_point.recharge_adjustment',
  'hspf_sudbury66_huc10_polygon.hydrology',
  'hspf_sudbury66_huc10_polygon.recharge_adjustment',
  'hspf_sudbury66_huc12_point.hydrology',
  'hspf_sudbury66_huc12_point.recharge_adjustment',
  'hspf_sudbury66_huc12_polygon.hydrology',
  'hspf_sudbury66_huc12_polygon.recharge_adjustment',
  'hspf_sudbury67_huc10_point.hydrology',
  'hspf_sudbury67_huc10_point.recharge_adjustment',
  'hspf_sudbury67_huc10_polygon.hydrology',
  'hspf_sudbury67_huc10_polygon.recharge_adjustment',
  'hspf_sudbury67_huc12_point.hydrology',
  'hspf_sudbury67_huc12_point.recharge_adjustment',
  'hspf_sudbury67_huc12_polygon.hydrology',
  'hspf_sudbury67_huc12_polygon.recharge_adjustment',
  'hspf_sudbury68_huc10_point.hydrology',
  'hspf_sudbury68_huc10_point.recharge_adjustment',
  'hspf_sudbury68_huc10_polygon.hydrology',
  'hspf_sudbury68_huc10_polygon.recharge_adjustment',
  'hspf_sudbury68_huc12_point.hydrology',
  'hspf_sudbury68_huc12_point.recharge_adjustment',
  'hspf_sudbury68_huc12_polygon.hydrology',
  'hspf_sudbury68_huc12_polygon.recharge_adjustment',
  'hspf_sudbury7_huc10_point.hydrology',
  'hspf_sudbury7_huc10_point.recharge_adjustment',
  'hspf_sudbury7_huc10_polygon.hydrology',
  'hspf_sudbury7_huc10_polygon.recharge_adjustment',
  'hspf_sudbury7_huc12_point.hydrology',
  'hspf_sudbury7_huc12_point.recharge_adjustment',
  'hspf_sudbury7_huc12_polygon.hydrology',
  'hspf_sudbury7_huc12_polygon.recharge_adjustment',
  'hspf_taunton16_huc10_point.hydrology',
  'hspf_taunton16_huc10_point.recharge_adjustment',
  'hspf_taunton16_huc10_polygon.hydrology',
  'hspf_taunton16_huc10_polygon.recharge_adjustment',
  'hspf_taunton16_huc12_point.hydrology',
  'hspf_taunton16_huc12_point.recharge_adjustment',
  'hspf_taunton16_huc12_polygon.hydrology',
  'hspf_taunton16_huc12_polygon.recharge_adjustment',
  'hspf_taunton69_huc10_point.hydrology',
  'hspf_taunton69_huc10_point.recharge_adjustment',
  'hspf_taunton69_huc10_polygon.hydrology',
  'hspf_taunton69_huc10_polygon.recharge_adjustment',
  'hspf_taunton69_huc12_point.hydrology',
  'hspf_taunton69_huc12_point.recharge_adjustment',
  'hspf_taunton69_huc12_polygon.hydrology',
  'hspf_taunton69_huc12_polygon.recharge_adjustment',
  'hspf_taunton70_huc10_point.hydrology',
  'hspf_taunton70_huc10_point.recharge_adjustment',
  'hspf_taunton70_huc10_polygon.hydrology',
  'hspf_taunton70_huc10_polygon.recharge_adjustment',
  'hspf_taunton70_huc12_point.hydrology',
  'hspf_taunton70_huc12_point.recharge_adjustment',
  'hspf_taunton70_huc12_polygon.hydrology',
  'hspf_taunton70_huc12_polygon.recharge_adjustment',
  'hspf_taunton71_huc10_point.hydrology',
  'hspf_taunton71_huc10_point.recharge_adjustment',
  'hspf_taunton71_huc10_polygon.hydrology',
  'hspf_taunton71_huc10_polygon.recharge_adjustment',
  'hspf_taunton71_huc12_point.hydrology',
  'hspf_taunton71_huc12_point.recharge_adjustment',
  'hspf_taunton71_huc12_polygon.hydrology',
  'hspf_taunton71_huc12_polygon.recharge_adjustment',
  'hspf_taunton72_huc10_point.hydrology',
  'hspf_taunton72_huc10_point.recharge_adjustment',
  'hspf_taunton72_huc10_polygon.hydrology',
  'hspf_taunton72_huc10_polygon.recharge_adjustment',
  'hspf_taunton72_huc12_point.hydrology',
  'hspf_taunton72_huc12_point.recharge_adjustment',
  'hspf_taunton72_huc12_polygon.hydrology',
  'hspf_taunton72_huc12_polygon.recharge_adjustment',
  'hspf_taunton73_huc10_point.hydrology',
  'hspf_taunton73_huc10_point.recharge_adjustment',
  'hspf_taunton73_huc10_polygon.hydrology',
  'hspf_taunton73_huc10_polygon.recharge_adjustment',
  'hspf_taunton73_huc12_point.hydrology',
  'hspf_taunton73_huc12_point.recharge_adjustment',
  'hspf_taunton73_huc12_polygon.hydrology',
  'hspf_taunton73_huc12_polygon.recharge_adjustment',
  'hspf_taunton8_huc10_point.hydrology',
  'hspf_taunton8_huc10_point.recharge_adjustment',
  'hspf_taunton8_huc10_polygon.hydrology',
  'hspf_taunton8_huc10_polygon.recharge_adjustment',
  'hspf_taunton8_huc12_point.hydrology',
  'hspf_taunton8_huc12_point.recharge_adjustment',
  'hspf_taunton8_huc12_polygon.hydrology',
  'hspf_taunton8_huc12_polygon.recharge_adjustment',
  'swmm_blackstone17_huc10_point.hydrology',
  'swmm_blackstone17_huc10_point.loadings_tn',
  'swmm_blackstone17_huc10_point.loadings_tp',
  'swmm_blackstone17_huc10_point.loadings_tss',
  'swmm_blackstone17_huc10_point.loadings_zn',
  'swmm_blackstone17_huc10_polygon.hydrology',
  'swmm_blackstone17_huc10_polygon.loadings_tn',
  'swmm_blackstone17_huc10_polygon.loadings_tp',
  'swmm_blackstone17_huc10_polygon.loadings_tss',
  'swmm_blackstone17_huc10_polygon.loadings_zn',
  'swmm_blackstone17_huc12_point.hydrology',
  'swmm_blackstone17_huc12_point.loadings_tn',
  'swmm_blackstone17_huc12_point.loadings_tp',
  'swmm_blackstone17_huc12_point.loadings_tss',
  'swmm_blackstone17_huc12_point.loadings_zn',
  'swmm_blackstone17_huc12_polygon.hydrology',
  'swmm_blackstone17_huc12_polygon.loadings_tn',
  'swmm_blackstone17_huc12_polygon.loadings_tp',
  'swmm_blackstone17_huc12_polygon.loadings_tss',
  'swmm_blackstone17_huc12_polygon.loadings_zn',
  'swmm_blackstone18_huc10_point.hydrology',
  'swmm_blackstone18_huc10_point.loadings_tn',
  'swmm_blackstone18_huc10_point.loadings_tp',
  'swmm_blackstone18_huc10_point.loadings_tss',
  'swmm_blackstone18_huc10_point.loadings_zn',
  'swmm_blackstone18_huc10_polygon.hydrology',
  'swmm_blackstone18_huc10_polygon.loadings_tn',
  'swmm_blackstone18_huc10_polygon.loadings_tp',
  'swmm_blackstone18_huc10_polygon.loadings_tss',
  'swmm_blackstone18_huc10_polygon.loadings_zn',
  'swmm_blackstone18_huc12_point.hydrology',
  'swmm_blackstone18_huc12_point.loadings_tn',
  'swmm_blackstone18_huc12_point.loadings_tp',
  'swmm_blackstone18_huc12_point.loadings_tss',
  'swmm_blackstone18_huc12_point.loadings_zn',
  'swmm_blackstone18_huc12_polygon.hydrology',
  'swmm_blackstone18_huc12_polygon.loadings_tn',
  'swmm_blackstone18_huc12_polygon.loadings_tp',
  'swmm_blackstone18_huc12_polygon.loadings_tss',
  'swmm_blackstone18_huc12_polygon.loadings_zn',
  'swmm_blackstone25_huc10_point.hydrology',
  'swmm_blackstone25_huc10_point.loadings_tn',
  'swmm_blackstone25_huc10_point.loadings_tp',
  'swmm_blackstone25_huc10_point.loadings_tss',
  'swmm_blackstone25_huc10_point.loadings_zn',
  'swmm_blackstone25_huc10_polygon.hydrology',
  'swmm_blackstone25_huc10_polygon.loadings_tn',
  'swmm_blackstone25_huc10_polygon.loadings_tp',
  'swmm_blackstone25_huc10_polygon.loadings_tss',
  'swmm_blackstone25_huc10_polygon.loadings_zn',
  'swmm_blackstone25_huc12_point.hydrology',
  'swmm_blackstone25_huc12_point.loadings_tn',
  'swmm_blackstone25_huc12_point.loadings_tp',
  'swmm_blackstone25_huc12_point.loadings_tss',
  'swmm_blackstone25_huc12_point.loadings_zn',
  'swmm_blackstone25_huc12_polygon.hydrology',
  'swmm_blackstone25_huc12_polygon.loadings_tn',
  'swmm_blackstone25_huc12_polygon.loadings_tp',
  'swmm_blackstone25_huc12_polygon.loadings_tss',
  'swmm_blackstone25_huc12_polygon.loadings_zn',
  'swmm_blackstone26_huc10_point.hydrology',
  'swmm_blackstone26_huc10_point.loadings_tn',
  'swmm_blackstone26_huc10_point.loadings_tp',
  'swmm_blackstone26_huc10_point.loadings_tss',
  'swmm_blackstone26_huc10_point.loadings_zn',
  'swmm_blackstone26_huc10_polygon.hydrology',
  'swmm_blackstone26_huc10_polygon.loadings_tn',
  'swmm_blackstone26_huc10_polygon.loadings_tp',
  'swmm_blackstone26_huc10_polygon.loadings_tss',
  'swmm_blackstone26_huc10_polygon.loadings_zn',
  'swmm_blackstone26_huc12_point.hydrology',
  'swmm_blackstone26_huc12_point.loadings_tn',
  'swmm_blackstone26_huc12_point.loadings_tp',
  'swmm_blackstone26_huc12_point.loadings_tss',
  'swmm_blackstone26_huc12_point.loadings_zn',
  'swmm_blackstone26_huc12_polygon.hydrology',
  'swmm_blackstone26_huc12_polygon.loadings_tn',
  'swmm_blackstone26_huc12_polygon.loadings_tp',
  'swmm_blackstone26_huc12_polygon.loadings_tss',
  'swmm_blackstone26_huc12_polygon.loadings_zn',
  'swmm_charles19_huc10_point.hydrology',
  'swmm_charles19_huc10_point.loadings_tn',
  'swmm_charles19_huc10_point.loadings_tp',
  'swmm_charles19_huc10_point.loadings_tss',
  'swmm_charles19_huc10_point.loadings_zn',
  'swmm_charles19_huc10_polygon.hydrology',
  'swmm_charles19_huc10_polygon.loadings_tn',
  'swmm_charles19_huc10_polygon.loadings_tp',
  'swmm_charles19_huc10_polygon.loadings_tss',
  'swmm_charles19_huc10_polygon.loadings_zn',
  'swmm_charles19_huc12_point.hydrology',
  'swmm_charles19_huc12_point.loadings_tn',
  'swmm_charles19_huc12_point.loadings_tp',
  'swmm_charles19_huc12_point.loadings_tss',
  'swmm_charles19_huc12_point.loadings_zn',
  'swmm_charles19_huc12_polygon.hydrology',
  'swmm_charles19_huc12_polygon.loadings_tn',
  'swmm_charles19_huc12_polygon.loadings_tp',
  'swmm_charles19_huc12_polygon.loadings_tss',
  'swmm_charles19_huc12_polygon.loadings_zn',
  'swmm_charles27_huc10_point.hydrology',
  'swmm_charles27_huc10_point.loadings_tn',
  'swmm_charles27_huc10_point.loadings_tp',
  'swmm_charles27_huc10_point.loadings_tss',
  'swmm_charles27_huc10_point.loadings_zn',
  'swmm_charles27_huc10_polygon.hydrology',
  'swmm_charles27_huc10_polygon.loadings_tn',
  'swmm_charles27_huc10_polygon.loadings_tp',
  'swmm_charles27_huc10_polygon.loadings_tss',
  'swmm_charles27_huc10_polygon.loadings_zn',
  'swmm_charles27_huc12_point.hydrology',
  'swmm_charles27_huc12_point.loadings_tn',
  'swmm_charles27_huc12_point.loadings_tp',
  'swmm_charles27_huc12_point.loadings_tss',
  'swmm_charles27_huc12_point.loadings_zn',
  'swmm_charles27_huc12_polygon.hydrology',
  'swmm_charles27_huc12_polygon.loadings_tn',
  'swmm_charles27_huc12_polygon.loadings_tp',
  'swmm_charles27_huc12_polygon.loadings_tss',
  'swmm_charles27_huc12_polygon.loadings_zn',
  'swmm_farmington20_huc10_point.hydrology',
  'swmm_farmington20_huc10_point.loadings_tn',
  'swmm_farmington20_huc10_point.loadings_tp',
  'swmm_farmington20_huc10_point.loadings_tss',
  'swmm_farmington20_huc10_point.loadings_zn',
  'swmm_farmington20_huc10_polygon.hydrology',
  'swmm_farmington20_huc10_polygon.loadings_tn',
  'swmm_farmington20_huc10_polygon.loadings_tp',
  'swmm_farmington20_huc10_polygon.loadings_tss',
  'swmm_farmington20_huc10_polygon.loadings_zn',
  'swmm_farmington20_huc12_point.hydrology',
  'swmm_farmington20_huc12_point.loadings_tn',
  'swmm_farmington20_huc12_point.loadings_tp',
  'swmm_farmington20_huc12_point.loadings_tss',
  'swmm_farmington20_huc12_point.loadings_zn',
  'swmm_farmington20_huc12_polygon.hydrology',
  'swmm_farmington20_huc12_polygon.loadings_tn',
  'swmm_farmington20_huc12_polygon.loadings_tp',
  'swmm_farmington20_huc12_polygon.loadings_tss',
  'swmm_farmington20_huc12_polygon.loadings_zn',
  'swmm_farmington28_huc10_point.hydrology',
  'swmm_farmington28_huc10_point.loadings_tn',
  'swmm_farmington28_huc10_point.loadings_tp',
  'swmm_farmington28_huc10_point.loadings_tss',
  'swmm_farmington28_huc10_point.loadings_zn',
  'swmm_farmington28_huc10_polygon.hydrology',
  'swmm_farmington28_huc10_polygon.loadings_tn',
  'swmm_farmington28_huc10_polygon.loadings_tp',
  'swmm_farmington28_huc10_polygon.loadings_tss',
  'swmm_farmington28_huc10_polygon.loadings_zn',
  'swmm_farmington28_huc12_point.hydrology',
  'swmm_farmington28_huc12_point.loadings_tn',
  'swmm_farmington28_huc12_point.loadings_tp',
  'swmm_farmington28_huc12_point.loadings_tss',
  'swmm_farmington28_huc12_point.loadings_zn',
  'swmm_farmington28_huc12_polygon.hydrology',
  'swmm_farmington28_huc12_polygon.loadings_tn',
  'swmm_farmington28_huc12_polygon.loadings_tp',
  'swmm_farmington28_huc12_polygon.loadings_tss',
  'swmm_farmington28_huc12_polygon.loadings_zn',
  'swmm_ipswich21_huc10_point.hydrology',
  'swmm_ipswich21_huc10_point.loadings_tn',
  'swmm_ipswich21_huc10_point.loadings_tp',
  'swmm_ipswich21_huc10_point.loadings_tss',
  'swmm_ipswich21_huc10_point.loadings_zn',
  'swmm_ipswich21_huc10_polygon.hydrology',
  'swmm_ipswich21_huc10_polygon.loadings_tn',
  'swmm_ipswich21_huc10_polygon.loadings_tp',
  'swmm_ipswich21_huc10_polygon.loadings_tss',
  'swmm_ipswich21_huc10_polygon.loadings_zn',
  'swmm_ipswich21_huc12_point.hydrology',
  'swmm_ipswich21_huc12_point.loadings_tn',
  'swmm_ipswich21_huc12_point.loadings_tp',
  'swmm_ipswich21_huc12_point.loadings_tss',
  'swmm_ipswich21_huc12_point.loadings_zn',
  'swmm_ipswich21_huc12_polygon.hydrology',
  'swmm_ipswich21_huc12_polygon.loadings_tn',
  'swmm_ipswich21_huc12_polygon.loadings_tp',
  'swmm_ipswich21_huc12_polygon.loadings_tss',
  'swmm_ipswich21_huc12_polygon.loadings_zn',
  'swmm_ipswich29_huc10_point.hydrology',
  'swmm_ipswich29_huc10_point.loadings_tn',
  'swmm_ipswich29_huc10_point.loadings_tp',
  'swmm_ipswich29_huc10_point.loadings_tss',
  'swmm_ipswich29_huc10_point.loadings_zn',
  'swmm_ipswich29_huc10_polygon.hydrology',
  'swmm_ipswich29_huc10_polygon.loadings_tn',
  'swmm_ipswich29_huc10_polygon.loadings_tp',
  'swmm_ipswich29_huc10_polygon.loadings_tss',
  'swmm_ipswich29_huc10_polygon.loadings_zn',
  'swmm_ipswich29_huc12_point.hydrology',
  'swmm_ipswich29_huc12_point.loadings_tn',
  'swmm_ipswich29_huc12_point.loadings_tp',
  'swmm_ipswich29_huc12_point.loadings_tss',
  'swmm_ipswich29_huc12_point.loadings_zn',
  'swmm_ipswich29_huc12_polygon.hydrology',
  'swmm_ipswich29_huc12_polygon.loadings_tn',
  'swmm_ipswich29_huc12_polygon.loadings_tp',
  'swmm_ipswich29_huc12_polygon.loadings_tss',
  'swmm_ipswich29_huc12_polygon.loadings_zn',
  'swmm_pawcatuck22_huc10_point.hydrology',
  'swmm_pawcatuck22_huc10_point.loadings_tn',
  'swmm_pawcatuck22_huc10_point.loadings_tp',
  'swmm_pawcatuck22_huc10_point.loadings_tss',
  'swmm_pawcatuck22_huc10_point.loadings_zn',
  'swmm_pawcatuck22_huc10_polygon.hydrology',
  'swmm_pawcatuck22_huc10_polygon.loadings_tn',
  'swmm_pawcatuck22_huc10_polygon.loadings_tp',
  'swmm_pawcatuck22_huc10_polygon.loadings_tss',
  'swmm_pawcatuck22_huc10_polygon.loadings_zn',
  'swmm_pawcatuck22_huc12_point.hydrology',
  'swmm_pawcatuck22_huc12_point.loadings_tn',
  'swmm_pawcatuck22_huc12_point.loadings_tp',
  'swmm_pawcatuck22_huc12_point.loadings_tss',
  'swmm_pawcatuck22_huc12_point.loadings_zn',
  'swmm_pawcatuck22_huc12_polygon.hydrology',
  'swmm_pawcatuck22_huc12_polygon.loadings_tn',
  'swmm_pawcatuck22_huc12_polygon.loadings_tp',
  'swmm_pawcatuck22_huc12_polygon.loadings_tss',
  'swmm_pawcatuck22_huc12_polygon.loadings_zn',
  'swmm_pawcatuck30_huc10_point.hydrology',
  'swmm_pawcatuck30_huc10_point.loadings_tn',
  'swmm_pawcatuck30_huc10_point.loadings_tp',
  'swmm_pawcatuck30_huc10_point.loadings_tss',
  'swmm_pawcatuck30_huc10_point.loadings_zn',
  'swmm_pawcatuck30_huc10_polygon.hydrology',
  'swmm_pawcatuck30_huc10_polygon.loadings_tn',
  'swmm_pawcatuck30_huc10_polygon.loadings_tp',
  'swmm_pawcatuck30_huc10_polygon.loadings_tss',
  'swmm_pawcatuck30_huc10_polygon.loadings_zn',
  'swmm_pawcatuck30_huc12_point.hydrology',
  'swmm_pawcatuck30_huc12_point.loadings_tn',
  'swmm_pawcatuck30_huc12_point.loadings_tp',
  'swmm_pawcatuck30_huc12_point.loadings_tss',
  'swmm_pawcatuck30_huc12_point.loadings_zn',
  'swmm_pawcatuck30_huc12_polygon.hydrology',
  'swmm_pawcatuck30_huc12_polygon.loadings_tn',
  'swmm_pawcatuck30_huc12_polygon.loadings_tp',
  'swmm_pawcatuck30_huc12_polygon.loadings_tss',
  'swmm_pawcatuck30_huc12_polygon.loadings_zn',
  'swmm_sudbury23_huc10_point.hydrology',
  'swmm_sudbury23_huc10_point.loadings_tn',
  'swmm_sudbury23_huc10_point.loadings_tp',
  'swmm_sudbury23_huc10_point.loadings_tss',
  'swmm_sudbury23_huc10_point.loadings_zn',
  'swmm_sudbury23_huc10_polygon.hydrology',
  'swmm_sudbury23_huc10_polygon.loadings_tn',
  'swmm_sudbury23_huc10_polygon.loadings_tp',
  'swmm_sudbury23_huc10_polygon.loadings_tss',
  'swmm_sudbury23_huc10_polygon.loadings_zn',
  'swmm_sudbury23_huc12_point.hydrology',
  'swmm_sudbury23_huc12_point.loadings_tn',
  'swmm_sudbury23_huc12_point.loadings_tp',
  'swmm_sudbury23_huc12_point.loadings_tss',
  'swmm_sudbury23_huc12_point.loadings_zn',
  'swmm_sudbury23_huc12_polygon.hydrology',
  'swmm_sudbury23_huc12_polygon.loadings_tn',
  'swmm_sudbury23_huc12_polygon.loadings_tp',
  'swmm_sudbury23_huc12_polygon.loadings_tss',
  'swmm_sudbury23_huc12_polygon.loadings_zn',
  'swmm_sudbury31_huc10_point.hydrology',
  'swmm_sudbury31_huc10_point.loadings_tn',
  'swmm_sudbury31_huc10_point.loadings_tp',
  'swmm_sudbury31_huc10_point.loadings_tss',
  'swmm_sudbury31_huc10_point.loadings_zn',
  'swmm_sudbury31_huc10_polygon.hydrology',
  'swmm_sudbury31_huc10_polygon.loadings_tn',
  'swmm_sudbury31_huc10_polygon.loadings_tp',
  'swmm_sudbury31_huc10_polygon.loadings_tss',
  'swmm_sudbury31_huc10_polygon.loadings_zn',
  'swmm_sudbury31_huc12_point.hydrology',
  'swmm_sudbury31_huc12_point.loadings_tn',
  'swmm_sudbury31_huc12_point.loadings_tp',
  'swmm_sudbury31_huc12_point.loadings_tss',
  'swmm_sudbury31_huc12_point.loadings_zn',
  'swmm_sudbury31_huc12_polygon.hydrology',
  'swmm_sudbury31_huc12_polygon.loadings_tn',
  'swmm_sudbury31_huc12_polygon.loadings_tp',
  'swmm_sudbury31_huc12_polygon.loadings_tss',
  'swmm_sudbury31_huc12_polygon.loadings_zn',
  'swmm_taunton24_huc10_point.hydrology',
  'swmm_taunton24_huc10_point.loadings_tn',
  'swmm_taunton24_huc10_point.loadings_tp',
  'swmm_taunton24_huc10_point.loadings_tss',
  'swmm_taunton24_huc10_point.loadings_zn',
  'swmm_taunton24_huc10_polygon.hydrology',
  'swmm_taunton24_huc10_polygon.loadings_tn',
  'swmm_taunton24_huc10_polygon.loadings_tp',
  'swmm_taunton24_huc10_polygon.loadings_tss',
  'swmm_taunton24_huc10_polygon.loadings_zn',
  'swmm_taunton24_huc12_point.hydrology',
  'swmm_taunton24_huc12_point.loadings_tn',
  'swmm_taunton24_huc12_point.loadings_tp',
  'swmm_taunton24_huc12_point.loadings_tss',
  'swmm_taunton24_huc12_point.loadings_zn',
  'swmm_taunton24_huc12_polygon.hydrology',
  'swmm_taunton24_huc12_polygon.loadings_tn',
  'swmm_taunton24_huc12_polygon.loadings_tp',
  'swmm_taunton24_huc12_polygon.loadings_tss',
  'swmm_taunton24_huc12_polygon.loadings_zn',
  'swmm_taunton32_huc10_point.hydrology',
  'swmm_taunton32_huc10_point.loadings_tn',
  'swmm_taunton32_huc10_point.loadings_tp',
  'swmm_taunton32_huc10_point.loadings_tss',
  'swmm_taunton32_huc10_point.loadings_zn',
  'swmm_taunton32_huc10_polygon.hydrology',
  'swmm_taunton32_huc10_polygon.loadings_tn',
  'swmm_taunton32_huc10_polygon.loadings_tp',
  'swmm_taunton32_huc10_polygon.loadings_tss',
  'swmm_taunton32_huc10_polygon.loadings_zn',
  'swmm_taunton32_huc12_point.hydrology',
  'swmm_taunton32_huc12_point.loadings_tn',
  'swmm_taunton32_huc12_point.loadings_tp',
  'swmm_taunton32_huc12_point.loadings_tss',
  'swmm_taunton32_huc12_point.loadings_zn',
  'swmm_taunton32_huc12_polygon.hydrology',
  'swmm_taunton32_huc12_polygon.loadings_tn',
  'swmm_taunton32_huc12_polygon.loadings_tp',
  'swmm_taunton32_huc12_polygon.loadings_tss',
  'swmm_taunton32_huc12_polygon.loadings_zn',
  'greenspace_housing',

  # Downloadable files:

  'wwtp_baseline_removal.csv',
  'urban_bmp_data.csv',
  'road_bmps.csv',
  'mun_state_all_usgs.csv',
  'mun_fed_all_usgs.csv'

);



my @regions = (
  'atlantic',
  'gulf',
  'pacific',
  'alaska',
  'hawaii',
  'interior_sw',
  'great_lakes',
  'missouri_river',
  'lower_mississippi',
  'upper_mississippi'
);

###############################################################################


sub get_coverage_metadata_reference { return \@coverage_metadata; }
sub get_region_metadata_reference { return \%region_metadata; }
sub get_regions { return @regions; }
sub get_coverages_reference { return \@coverages; }

1; # Required last executed statement of module. Do not add code below here!
