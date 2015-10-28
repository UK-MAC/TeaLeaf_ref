#!/usr/bin/python

import random
import sys
import re

class Params(object):
    num_states = 0
    x_cells = 10
    y_cells = 10

    xmin = 0.0
    ymin = 0.0
    xmax = 100.0
    ymax = 100.0

    tl_max_iters=10000
    tl_coefficient=1
    eps=10e-15

    profiler_on=0

    dtinit=0.004
    dtmax=1.0
    dtmin=0.0000001
    dtrise=1.5

    visit_frequency=0
    summary_frequency=10
    g_ibig=100
    g_small=1.0e-16
    g_big  =1.0e+21
    end_time=10.0
    end_step=g_ibig

    g_xdir=1
    g_ydir=2

    geoms = {"rectangle":1, "circle":2, "point":3}
    controls = {1:"sound", 2:"xvel", 3:"yvel", 4:"div"}

    num_fields = 6
    fields = {"density":0,
              "energy0":0, "energy1":0,
              "u":0,
              "p":0,
              "sd":0}
    correct_field_order = ["density",
                           "energy0",
                           "energy1",
                           "u",
                           "p",
                           "sd"]

    packing_field_numbers = {
        "density"   : 1,
        "energy0"   : 2,
        "energy1"   : 3,
        "u"         : 4,
        "p"         : 5,
        "sd"        : 6,
    }

    sizes = {
        "density" : (0, 0),
        "energy0" : (0, 0),
        "energy1" : (0, 0),
        "u" : (0, 0),
        "p" : (0, 0),
        "sd" : (0, 0),
    }

    CHUNK_LEFT = 1
    CHUNK_RIGHT = 2
    CHUNK_BOTTOM = 3
    CHUNK_TOP = 4
    EXTERNAL_FACE = -1

    def __init__(self, INFILE):
        """ Reads input from clover.in into params """

        input_params = []
        with open(INFILE) as in_file:
            input_params = in_file.read()

        print
        print "Running from input:"
        print 20*"-"
        print input_params,
        print 20*"-"
        print

        # gets a parameter from the input file in the form 'x=2.4'
        def getval(param, string):
            #pattern = param + r'[+-]?(\d+(\.\d*)?|\.\d+)([eE][+-]?\d+)?'
            #pattern = param + "([+-]?(?=\d*[.eE])(?=\.?\d)\d*\.?\d*(?:[eE][+-]?\d+)?|\d+)"
            #pattern = param + "=([-]?(\d*\.\d+|\d+))"

            pattern = param + r"=([+-]?(\d*\.)?\d+(e[+-]?\d+)?|\w+)"

            match = re.search(pattern, string)
            if not match:
                return None
            else:
                return match.group(1)

        # match state information
        states_read = [int(i) for i in re.findall("[^#]state (\d)", input_params)]
        if not len(states_read):
            raise Exception("No states defined")
        else:
            self.num_states = max(states_read)

        self.states = {}

        for ss in xrange(1, 1+self.num_states):

            # check there is info
            state_match = re.search("state %d (.*)" % ss, input_params)
            if not state_match or len(state_match.groups()) < 1:
                # if there was no info for state, raise exception
                raise Exception("State information not defined")

            # specification for state values
            state_info = state_match.group(1)

            # intiialise as empty
            self.states[ss] = {
                "energy":0.0,
                "density":0.0,
                "xvel":0.0,
                "yvel":0.0,
                "xmin":0.0,
                "ymin":0.0,
                "xmax":0.0,
                "ymax":0.0,
                "radius":0.0,
                "geometry":0
            }

            # gets match, sets in dict - none are needed?
            def match_state_param(param, type_func):
                x = getval(param, state_info)
                if not x:
                    return
                else:
                    self.states[ss][param] = type_func(x)

            # geometry
            geom_convert = lambda x: self.geoms[x]
            geom_pattern = "geometry=(\w+)"
            match = re.search(geom_pattern, state_info)
            if match:
                self.states[ss]["geometry"] = geom_convert(match.group(1))

            # stuff
            match_state_param("xmin", float)
            match_state_param("ymin", float)
            match_state_param("xmax", float)
            match_state_param("ymax", float)
            match_state_param("radius", float)

            # fluid self?
            match_state_param("density", float)
            match_state_param("energy", float)

        # removes all lines with "state" in them
        # so state info isnt matched by regex
        no_states = re.sub("state.*", " ", input_params)
        no_states = re.sub("!.*", " ", no_states)

        # match other paramaters - float or int
        def match_param(param, type_func, needed=False, attr=None, default=None):
            # if a default value is set, use it
            if default != None:
                setattr(self, param, type_func(default))

            x = getval(param, no_states)
            if not x and needed:
                raise Exception(param + " not defined but needed")
            elif not x:
                return
            else:
                if attr:
                    setattr(self, attr, type_func(x))
                else:
                    setattr(self, param, type_func(x))

        def check_param_enabled(param):
            match = re.search(param, no_states)
            if not match:
                setattr(self, param, False)
            else:
                setattr(self, param, True)

        # size
        match_param("x_cells", int)
        match_param("y_cells", int)

        # stuff
        match_param("xmin", float)
        match_param("ymin", float)
        match_param("xmax", float)
        match_param("ymax", float)

        # tealeaf
        match_param("tl_max_iters", int, False)
        match_param("tl_coefficient", int, False)
        match_param("tl_eps", float, needed=True)
        match_param("tiles_per_task", int, False, default=1)

        check_param_enabled("tl_use_cg")
        check_param_enabled("tl_use_dpcg")
        check_param_enabled("reflective_boundary")

        if not (self.tl_use_cg or self.tl_use_dpcg):
            raise Exception("Must specify tl_use_cg or tl_use_dpcg in input file")

        # preconditioners
        match_param("tl_preconditioner_type", str, default="none")

        # timestep things
        match_param("initial_timestep", float, attr="dtinit")
        match_param("timestep_rise", float, attr="dtrise")
        match_param("max_timestep", float, attr="dtmax")
        match_param("end_time", float)
        match_param("end_step", int)

        # correction thing
        dx = (self.xmax - self.xmin) / self.x_cells
        dy = (self.ymax - self.ymin) / self.y_cells

        for ss in xrange(2, 1+self.num_states):
            self.states[ss]["xmin"] += dx/100.0
            self.states[ss]["ymin"] += dy/100.0
            self.states[ss]["xmax"] -= dx/100.0
            self.states[ss]["ymax"] -= dy/100.0

