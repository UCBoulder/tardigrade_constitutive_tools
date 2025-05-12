from libcpp.vector cimport vector

import numpy as np
cimport numpy as np

cimport tardigrade_error_tools_python


cdef extern from "tardigrade_constitutive_tools.h" namespace "tardigradeConstitutiveTools":

    void decomposeGreenLagrangeStrain(const vector[double] &, vector[double] &, double &)

    void midpointEvolution(const double &, const vector[double] &,\
                           const vector[double] &, const vector[double] &,\
                           vector[double] &, vector[double] &, const vector[double] &)

    void midpointEvolution(const double &, const vector[double] &,\
                           const vector[double] &, const vector[double] &,\
                           vector[double] &, vector[double] &, vector[vector[double]] &,\
                           const vector[double] &)
