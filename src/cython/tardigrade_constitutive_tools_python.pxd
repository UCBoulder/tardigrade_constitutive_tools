from libcpp.vector cimport vector

import numpy as np
cimport numpy as np


cdef extern from "tardigrade_constitutive_tools.h" namespace "tardigradeConstitutiveTools":

    void decomposeGreenLagrangeStrain(const vector[double] &, vector[double] &, double &)

    void midpointEvolution(const double &, const vector[double] &,\
                           const vector[double] &, const vector[double] &,\
                           vector[double] &, vector[double] &, const vector[double] &)

    void midpointEvolution(const double &, const vector[double] &,\
                           const vector[double] &, const vector[double] &,\
                           vector[double] &, vector[double] &, vector[vector[double]] &,\
                           const vector[double] &)
