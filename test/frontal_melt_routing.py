import PISM
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

context = PISM.Context()
ctx = context.ctx
config = context.config
log = ctx.log()

def geometry(grid, input_file):
    "Allocate storage for ice geometry and initialize it from an input file."

    geometry = PISM.Geometry(grid)
    geometry.ice_thickness.regrid(input_file, critical=True)
    geometry.bed_elevation.regrid(input_file, critical=True)
    min_thickness = config.get_double("geometry.ice_free_thickness_standard")
    geometry.ensure_consistency(min_thickness)

    return geometry

def create_potential_temperature(grid, file_name, theta=274.15):

    PISM.util.prepare_output(file_name)

    Th = PISM.IceModelVec2S(grid, "theta_ocean", PISM.WITHOUT_GHOSTS)
    Th.set_attrs("climate", "potential temperature", "Kelvin", "")
    Th.set(theta)
    Th.write(file_name)

def hydrology(grid):
    """Allocate and initialize the hydrology model, starting with zero
    till and transportable water amounts. Also, note that the water
    pressure is not used in the initialization of the routing model,
    so using zero pressure is OK here.

    """

    Wtill = PISM.IceModelVec2S(grid, "tillwat", PISM.WITHOUT_GHOSTS)
    Wtill.set_attrs("model_state", "effective thickness of subglacial water stored in till", "m", "")
    Wtill.set(0.0)

    P = PISM.IceModelVec2S(grid, "overburden_pressure", PISM.WITHOUT_GHOSTS);
    P.set_attrs("internal", "overburden pressure", "Pa", "");
    P.set(0.0);

    hydrology = PISM.RoutingHydrology(grid)
    hydrology.init(Wtill, Wtill, P)

    return hydrology

def frontal_melt(grid):
    "Allocate and initialize the frontal melt model."
    fmelt = PISM.FrontalMeltDischargeRouting(grid)
    fmelt.init(geometry)

    return fmelt

if __name__ == "__main__":
    # set up the option parser
    parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
    parser.description = """Uses PISM's routing hydrology and frontal melt models to compute
frontal melt corresponding to provided surface water input rates."""
    parser.add_argument("BOOTSTRAP_FILE", nargs=1, help="Bootstrap file", default=None)
    parser.add_argument("--routing_file", dest="routing_file", help="routing file",
                        default="JIB_DMI-HIRHAM5_GL2_ERAI_1980_2014_MRROS_YMM_EPSG3413_3600m_0.nc")
    parser.add_argument("--th_file", dest="th_file", help="Thermal forcing file", default=None)
    parser.add_argument("--theta", dest="theta", help="Thermal forcing Default=274.15K", type=float, default=274.15)
    parser.add_argument("-o", dest="output_file", help="output file name",
                        default="frontal_melt.nc")

    options = parser.parse_args()
    input_file = options.BOOTSTRAP_FILE[0]
    routing_file = options.routing_file
    th_file = options.th_file
    output_file = options.output_file
    theta = options.theta
    
    # create the grid
    registration = PISM.CELL_CORNER
    grid = PISM.IceGrid.FromFile(ctx, input_file, ("bed", "thickness"), registration)

    if th_file is None:
        th_file = "th_file.nc"
        create_potential_temperature(grid, th_file, theta=theta)

    # initialize ice geometry
    geometry = geometry(grid, input_file)

    # create the hydrology model
    config.set_double("hydrology.tillwat_max", 0.0)
    config.set_string("hydrology.surface_input_file", routing_file)
    hydrology = hydrology(grid)

    # create the frontal melt model
    config.set_string("frontal_melt.routing.file", th_file)
    frontal_melt = frontal_melt(grid)

    # initialize model time; note that the calendar ("standard") will
    # be reset by init_from_file()
    time = PISM.Time_Calendar(ctx.com(), config, "standard", ctx.unit_system())
    time.init_from_file(routing_file, log, True) # do set start time

    # initialize water input rate forcing
    f = PISM.PIO(ctx.com(), "netcdf3", routing_file, PISM.PISM_READONLY)
    water_input_rate = PISM.IceModelVec2T_ForcingField(grid, f, "water_input_rate", "",
                                                       12, 12, False)
    f.close()
    water_input_rate.set_attrs("climate", "water input rate", "m s-1", "")

    # not periodic
    period = 0
    # reference time is irrelevant
    reference_time = 0
    water_input_rate.init(routing_file, period, reference_time)

    basal_melt_rate = PISM.IceModelVec2S(grid, "basal_melt_rate", PISM.WITHOUT_GHOSTS)
    basal_melt_rate.set(0.0)

    discharge = PISM.IceModelVec2S(grid, "discharge", PISM.WITHOUT_GHOSTS)
    discharge.set_attrs("internal", "discharge", "kg", "")
    discharge.set(0.0)

    hydro_inputs                    = PISM.HydrologyInputs()
    hydro_inputs.cell_type          = geometry.cell_type
    hydro_inputs.ice_thickness      = geometry.ice_thickness
    hydro_inputs.bed_elevation      = geometry.bed_elevation
    hydro_inputs.surface_input_rate = water_input_rate
    hydro_inputs.basal_melt_rate    = basal_melt_rate

    frontal_melt_inputs = PISM.FrontalMeltInputs()
    frontal_melt_inputs.geometry = geometry
    frontal_melt_inputs.subglacial_discharge_at_grounding_line = discharge

    time_name = config.get_string("time.dimension_name")
    output = PISM.PIO(ctx.com(), config.get_string("output.format"),
                      output_file, PISM.PISM_READWRITE_MOVE)
    PISM.define_time(output, time_name, time.calendar(), time.units_string(), ctx.unit_system())
    output.put_att_text(time_name, "bounds", "time_bounds")

    bounds = PISM.TimeBoundsMetadata("time_bounds", time_name, ctx.unit_system())
    bounds.set_string("units", time.units_string())

    # run models, stepping through time one record of forcing at at time
    output_record = 0

    while time.current() < time.end():
        t = time.current()
        dt = water_input_rate.max_timestep(t).value()

        water_input_rate.update(t, dt)
        water_input_rate.average(t, dt)

        log.message(2, "{}, dt={} days".format(time.date(), dt / 86400.0))

        hydrology.update(t, dt, hydro_inputs)

        # convert mass change at the grounding line (positive
        # corresponds to mass gain) to discharge (positive corresponds
        # to mass loss)
        discharge.copy_from(hydrology.mass_change_at_grounding_line())
        discharge.scale(-1.0)

        frontal_melt.update(frontal_melt_inputs, t, dt)

        time.step(dt)

        PISM.append_time(output, time_name, time.current())
        PISM.write_time_bounds(output, bounds, output_record, [t, t + dt], PISM.PISM_DOUBLE)
        output_record += 1

        frontal_melt.frontal_melt_rate().write(output)
        hydrology.subglacial_water_thickness().write(output)
        hydrology.total_input_rate().write(output)
        discharge.write(output)

    output.close()
