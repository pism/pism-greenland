import PISM
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

# set up the option parser
parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
parser.description = "Generating scripts for warming experiments."
parser.add_argument("FILE", nargs=1, help="Bootstrap file", default=None)
parser.add_argument("--routing_file", dest="routing_file", help="routing file", default="JIB_DMI-HIRHAM5_GL2_ERAI_1980_2014_MRROS_YMM_EPSG3413_3600m_0.nc")
parser.add_argument("--th_file", dest="th_file", help="Thermal forcing file", default=None)

options = parser.parse_args()
input_file = options.FILE[0]
routing_file = options.routing_file
th_file = options.th_file

context = PISM.Context()
ctx = context.ctx
config = context.config

registration = PISM.CELL_CENTER


grid = PISM.IceGrid.FromFile(ctx, input_file, ("bed", "thickness"), registration)
geometry = PISM.Geometry(grid)
geometry.ice_thickness.regrid(input_file, critical=True)
geometry.bed_elevation.regrid(input_file, critical=True)
min_thickness = config.get_double("geometry.ice_free_thickness_standard")
geometry.ensure_consistency(min_thickness)

potential_temperature = 0.0

if th_file is None:
    th_file = 'th_given.nc'
    PISM.util.prepare_output(th_file)

    Th = PISM.IceModelVec2S(grid, "theta_ocean", PISM.WITHOUT_GHOSTS)
    Th.set_attrs("climate", "potential temperature", "Kelvin", "")
    Th.set(potential_temperature)
    Th.write(th_file)

config.set_string("frontal_melt.routing.file", th_file)
config.set_string("hydrology.surface_input_file", routing_file)
config.set_double("hydrology.tillwat_max", 0.0)

inputs = PISM.FrontalMeltInputs()
cell_area = grid.dx() * grid.dy()
water_density = config.get_double("constants.fresh_water.density")

Wtill = PISM.IceModelVec2S(grid, "tillwat", PISM.WITHOUT_GHOSTS)
Wtill.set_attrs("model_state", "effective thickness of subglacial water stored in till", "m", "")
Wtill.set(0.0)

P = PISM.IceModelVec2S(grid, "overburden_pressure", PISM.WITHOUT_GHOSTS);
P.set_attrs("internal", "overburden pressure", "Pa", "");
P.set(0.0);

hydrology = PISM.RoutingHydrology(grid)
hydrology.initialize(Wtill, Wtill, P)

fmelt = PISM.FrontalMeltDischargeRouting(grid)
fmelt.init(geometry)
