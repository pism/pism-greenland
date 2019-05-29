import PISM

ctx = PISM.Context()
ctx.log.set_threshold(3)

geometry_filename = ctx.config.get_string("input.file")

if geometry_filename == "":
    geometry_filename = "pism_Greenland_ismip6_1000m_mcb_jpl_v4.nc"

grid = PISM.IceGrid.FromFile(ctx.ctx, geometry_filename, ["bed"], PISM.NOT_PERIODIC)

geometry = PISM.Geometry(grid)
geometry.ice_thickness.regrid(geometry_filename)
geometry.ice_area_specific_volume.set(0.0)
geometry.bed_elevation.regrid(geometry_filename)
geometry.sea_level_elevation.set(0.0)
geometry.ensure_consistency(0)

model = PISM.CalvingHayhurstCalving(grid)

model.init()

model.update(geometry.cell_type, geometry.ice_thickness, geometry.sea_level_elevation, geometry.bed_elevation)

model.calving_rate().dump(ctx.config.get_string("output.file_name"))
