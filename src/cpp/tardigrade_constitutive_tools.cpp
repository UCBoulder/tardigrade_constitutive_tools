/**
 *****************************************************************************
 * \file tardigrade_constitutive_tools.cpp
 *****************************************************************************
 * A collection of tools useful for constitutive models. These tools are
 * intended to be generalized expressions which perform operations commonly
 * encountered in the development of constitutive models. This will enable
 * users to develop new models quickly and (in principle) with less errors
 * resulting in delays.
 *
 * Developers should use tardigrade_vector_tools to perform vector multiplications and
 * matrix solves since this library can be independently checked. Also using
 * templates and typedef for the data is strongly encouraged so that single
 * and double precision values (in addition to other types) can be used
 * simply without significant reworking of the code. Passing by reference is
 * also encouraged if appropriate to allow users of FORTRAN to interface with
 * the code without extensive modifications.
 *****************************************************************************
 */

#include <tardigrade_constitutive_tools.h>

#include <algorithm>

namespace tardigradeConstitutiveTools {

    floatType deltaDirac(const unsigned int i, const unsigned int j) {
        /*!
         * The delta dirac function \f$\delta\f$
         *
         * if i==j return 1
         * if i!=j return 0
         *
         * \param i: The first index
         * \param j: The second index
         */

        return (floatType)(i == j);
    }

    void rotateMatrix(const floatVector &A, const floatVector &Q, floatVector &rotatedA) {
        /*!
         * Rotate a matrix \f$A\f$ using the orthogonal matrix \f$Q\f$ with the form
         *
         * \f$A'_{ij} = Q_{Ii} A_{IJ} Q_{Jj}\f$
         *
         * TODO: Generalize to non square matrices
         *
         * \param &A: The matrix to be rotated ( \f$A\f$ )
         * \param &Q: The rotation matrix ( \f$Q\f$Q )
         * \param &rotatedA: The rotated matrix ( \f$A'\f$ )
         */

        rotatedA = floatVector(A.size(), 0);
        rotateMatrix(std::cbegin(A), std::cend(A), std::cbegin(Q), std::cend(Q), std::sqrt(A.size()),
                     std::begin(rotatedA), std::end(rotatedA));

        return;
    }

    void computeDeformationGradient(const floatVector &displacementGradient, floatVector &F, const bool isCurrent) {
        /*!
         * Compute the deformation gradient from the gradient of the displacement
         *
         * If isCurrent = false
         *
         * \f$ \bf{F} = \frac{\partial \bf{u}}{\partial \bf{X} } u_i + \bf{I} \f$
         *
         * else if isCurrent = true
         *
         * \f$ \bf{F} = \left(\bf{I} - \frac{\partial \bf{u}}{\partial \bf{x}}\right)^{-1} \f$
         *
         * \param &displacementGradient: The gradient of the displacement with respect to either the
         *     current or previous position.
         * \param &F: The deformation gradient
         * \param &isCurrent: Boolean indicating whether the gradient is taken w.r.t. the current (true)
         *     or reference (false) position.
         */

        const unsigned int dim     = (unsigned int)std::pow(displacementGradient.size(), 0.5);
        const unsigned int sot_dim = dim * dim;

        F = floatVector(sot_dim, 0);

        TARDIGRADE_ERROR_TOOLS_CHECK((dim == 1) || (dim == 2) || (dim == 3),
                                     "The dimension of the displacement gradient is " + std::to_string(dim) +
                                         " but must be 1, 2, or 3");

        if (dim == 3) {
            computeDeformationGradient<3>(std::cbegin(displacementGradient), std::cend(displacementGradient),
                                          std::begin(F), std::end(F), isCurrent);
        } else if (dim == 2) {
            computeDeformationGradient<2>(std::cbegin(displacementGradient), std::cend(displacementGradient),
                                          std::begin(F), std::end(F), isCurrent);
        } else if (dim == 1) {
            computeDeformationGradient<1>(std::cbegin(displacementGradient), std::cend(displacementGradient),
                                          std::begin(F), std::end(F), isCurrent);
        }

        return;
    }

    void computeDeformationGradient(const floatVector &displacementGradient, floatVector &F, floatVector &dFdGradU,
                                    const bool isCurrent) {
        /*!
         * Compute the deformation gradient from the gradient of the displacement
         *
         * If isCurrent = false
         *
         * \f$ \bf{F} = \frac{\partial \bf{u}}{\partial \bf{X} } u_i + \bf{I} \f$
         *
         * else if isCurrent = true
         *
         * \f$ \bf{F} = \left(\bf{I} - \frac{\partial \bf{u}}{\partial \bf{x}}\right)^{-1} \f$
         *
         * \param &displacementGradient: The gradient of the displacement with respect to either the
         *     current or previous position.
         * \param &F: The deformation gradient
         * \param &dFdGradU: The derivative of the deformation gradient w.r.t. the displacement gradient
         * \param &isCurrent: Boolean indicating whether the gradient is taken w.r.t. the current (true)
         *     or reference (false) position.
         */

        const unsigned int dim     = (unsigned int)std::pow(displacementGradient.size(), 0.5);
        const unsigned int sot_dim = dim * dim;

        TARDIGRADE_ERROR_TOOLS_CHECK((dim == 1) || (dim == 2) || (dim == 3),
                                     "The dimension of the displacement gradient is " + std::to_string(dim) +
                                         " but must be 1, 2, or 3");

        F        = floatVector(sot_dim, 0);
        dFdGradU = floatVector(sot_dim * sot_dim, 0);

        if (dim == 3) {
            computeDeformationGradient<3>(std::begin(displacementGradient), std::end(displacementGradient),
                                          std::begin(F), std::end(F), std::begin(dFdGradU), std::end(dFdGradU),
                                          isCurrent);

        } else if (dim == 2) {
            computeDeformationGradient<2>(std::begin(displacementGradient), std::end(displacementGradient),
                                          std::begin(F), std::end(F), std::begin(dFdGradU), std::end(dFdGradU),
                                          isCurrent);

        } else if (dim == 1) {
            computeDeformationGradient<1>(std::begin(displacementGradient), std::end(displacementGradient),
                                          std::begin(F), std::end(F), std::begin(dFdGradU), std::end(dFdGradU),
                                          isCurrent);
        }

        return;
    }

    void computeRightCauchyGreen(const floatVector &deformationGradient, floatVector &C) {
        /*!
         * Compute the Right Cauchy-Green deformation tensor ( \f$C\f$ )
         *
         * \f$C_{IJ} = F_{iI} F_{iJ}\f$
         *
         * \param &deformationGradient: A reference to the deformation gradient ( \f$F\f$ )
         * \param &C: The resulting Right Cauchy-Green deformation tensor ( \f$C\f$ )
         *
         * The deformation gradient is organized as F11, F12, F13, F21, F22, F23, F31, F32, F33
         *
         * The Right Cauchy-Green deformation tensor is organized as C11, C12, C13, C21, C22, C23, C31, C32, C33
         */

        const unsigned int dim     = (unsigned int)std::pow(deformationGradient.size(), 0.5);
        const unsigned int sot_dim = dim * dim;

        TARDIGRADE_ERROR_TOOLS_CHECK((dim == 1) || (dim == 2) || (dim == 3),
                                     "The dimension of the deformation gradient is " + std::to_string(dim) +
                                         " but must be 1, 2, or 3");

        C = floatVector(sot_dim, 0);

        if (dim == 3) {
            computeRightCauchyGreen<3>(std::begin(deformationGradient), std::end(deformationGradient), std::begin(C),
                                       std::end(C));

        } else if (dim == 2) {
            computeRightCauchyGreen<2>(std::begin(deformationGradient), std::end(deformationGradient), std::begin(C),
                                       std::end(C));

        } else if (dim == 1) {
            computeRightCauchyGreen<1>(std::begin(deformationGradient), std::end(deformationGradient), std::begin(C),
                                       std::end(C));
        }

        return;
    }

    void computeRightCauchyGreen(const floatVector &deformationGradient, floatVector &C, floatMatrix &dCdF) {
        /*!
         * Compute the Right Cauchy-Green deformation tensor ( \f$C\f$ ) from the deformation gradient ( \f$F\f$ )
         *
         * \f$C_{IJ} = F_{iI} F_{iJ}\f$
         *
         * \param &deformationGradient: A reference to the deformation gradient ( \f$F\f$ )
         * \param &C: The resulting Right Cauchy-Green deformation tensor ( \f$C\f$ )
         * \param &dCdF: The Jacobian of the Right Cauchy-Green deformation tensor
         *     with regards to the deformation gradient ( \f$\frac{\partial C}{\partial F}\f$ ).
         *
         * The deformation gradient is organized as F11, F12, F13, F21, F22, F23, F31, F32, F33
         *
         * The Right Cauchy-Green deformation tensor is organized as C11, C12, C13, C21, C22, C23, C31, C32, C33
         */

        constexpr unsigned int dim     = 3;
        constexpr unsigned int sot_dim = dim * dim;

        floatVector _dCdF;

        TARDIGRADE_ERROR_TOOLS_CATCH(computeRightCauchyGreen(deformationGradient, C, _dCdF));

        dCdF = tardigradeVectorTools::inflate(_dCdF, sot_dim, sot_dim);

        return;
    }

    void computeRightCauchyGreen(const floatVector &deformationGradient, floatVector &C, floatVector &dCdF) {
        /*!
         * Compute the Right Cauchy-Green deformation tensor ( \f$C\f$ ) from the deformation gradient ( \f$F\f$ )
         *
         * \f$C_{IJ} = F_{iI} F_{iJ}\f$
         *
         * \param &deformationGradient: A reference to the deformation gradient ( \f$F\f$ )
         * \param &C: The resulting Right Cauchy-Green deformation tensor ( \f$C\f$ )
         * \param &dCdF: The Jacobian of the Right Cauchy-Green deformation tensor
         *     with regards to the deformation gradient ( \f$\frac{\partial C}{\partial F}\f$ ).
         *
         * The deformation gradient is organized as F11, F12, F13, F21, F22, F23, F31, F32, F33
         *
         * The Right Cauchy-Green deformation tensor is organized as C11, C12, C13, C21, C22, C23, C31, C32, C33
         */

        const unsigned int dim = (unsigned int)std::pow(deformationGradient.size(), 0.5);

        TARDIGRADE_ERROR_TOOLS_CHECK((dim == 3) || (dim == 2) || (dim == 1),
                                     "The dimension of the deformation gradient is " + std::to_string(dim) +
                                         " but must be 1, 2, or 3");

        C    = floatVector(dim * dim, 0);
        dCdF = floatVector(dim * dim * dim * dim, 0);

        if (dim == 3) {
            computeRightCauchyGreen<3>(std::begin(deformationGradient), std::end(deformationGradient), std::begin(C),
                                       std::end(C), std::begin(dCdF), std::end(dCdF));

        } else if (dim == 2) {
            computeRightCauchyGreen<2>(std::begin(deformationGradient), std::end(deformationGradient), std::begin(C),
                                       std::end(C), std::begin(dCdF), std::end(dCdF));

        } else if (dim == 1) {
            computeRightCauchyGreen<1>(std::begin(deformationGradient), std::end(deformationGradient), std::begin(C),
                                       std::end(C), std::begin(dCdF), std::end(dCdF));
        }

        return;
    }

    void computeGreenLagrangeStrain(const floatVector &deformationGradient, floatVector &E) {
        /*!
         * Compute the Green-Lagrange strain ( \f$E\f$ ) from the deformation gradient ( \f$F\f$ ). The operation is:
         *
         * \f$E = 0.5 (F_{iI} F_{iJ} - \delta_{IJ})\f$
         *
         * Where \f$F\f$ is the deformation gradient and \f$\delta\f$ is the kronecker delta.
         *
         * \param &deformationGradient: A reference to the deformation gradient ( \f$F\f$ ).
         * \param &E: The resulting Green-Lagrange strain ( \f$E\f$ ).
         *
         * The deformation gradient is organized as  F11, F12, F13, F21, F22, F23, F31, F32, F33
         *
         * The Green-Lagrange strain is organized as E11, E12, E13, E21, E22, E23, E31, E32, E33
         */

        const unsigned int dim = (unsigned int)std::pow(deformationGradient.size(), 0.5);

        TARDIGRADE_ERROR_TOOLS_CHECK((dim == 3) || (dim == 2) || (dim == 1),
                                     "The dimension of the deformation gradient is " + std::to_string(dim) +
                                         " but must be 1, 2, or 3");

        E = floatVector(dim * dim, 0);

        if (dim == 3) {
            computeGreenLagrangeStrain<3>(std::begin(deformationGradient), std::end(deformationGradient), std::begin(E),
                                          std::end(E));

        } else if (dim == 2) {
            computeGreenLagrangeStrain<2>(std::begin(deformationGradient), std::end(deformationGradient), std::begin(E),
                                          std::end(E));

        } else if (dim == 1) {
            computeGreenLagrangeStrain<1>(std::begin(deformationGradient), std::end(deformationGradient), std::begin(E),
                                          std::end(E));
        }

        return;
    }

    void computeGreenLagrangeStrain(const floatVector &deformationGradient, floatVector &E, floatMatrix &dEdF) {
        /*!
         * Compute the Green-Lagrange strain ( \f$E\f$ ) from the deformation gradient ( \f$F\f$ ) and it's jacobian.
         *
         * \param &deformationGradient: A reference to the deformation gradient ( \f$F\f$ ).
         * \param &E: The resulting Green-Lagrange strain ( \f$E\f$ ).
         * \param &dEdF: The jacobian of the Green-Lagrange strain w.r.t. the
         *     deformation gradient ( \f$\frac{\partial E}{\partial F}\f$ ).
         *
         * The deformation gradient is organized as  F11, F12, F13, F21, F22, F23, F31, F32, F33
         *
         * The Green-Lagrange strain is organized as E11, E12, E13, E21, E22, E23, E31, E32, E33
         */

        floatVector _dEdF;

        TARDIGRADE_ERROR_TOOLS_CATCH(computeGreenLagrangeStrain(deformationGradient, E, _dEdF));

        dEdF = tardigradeVectorTools::inflate(_dEdF, deformationGradient.size(), deformationGradient.size());

        return;
    }

    void computeGreenLagrangeStrain(const floatVector &deformationGradient, floatVector &E, floatVector &dEdF) {
        /*!
         * Compute the Green-Lagrange strain ( \f$E\f$ ) from the deformation gradient ( \f$F\f$ ) and it's jacobian.
         *
         * \param &deformationGradient: A reference to the deformation gradient ( \f$F\f$ ).
         * \param &E: The resulting Green-Lagrange strain ( \f$E\f$ ).
         * \param &dEdF: The jacobian of the Green-Lagrange strain w.r.t. the
         *     deformation gradient ( \f$\frac{\partial E}{\partial F}\f$ ).
         *
         * The deformation gradient is organized as  F11, F12, F13, F21, F22, F23, F31, F32, F33
         *
         * The Green-Lagrange strain is organized as E11, E12, E13, E21, E22, E23, E31, E32, E33
         */

        const unsigned int dim = (unsigned int)std::pow(deformationGradient.size(), 0.5);

        TARDIGRADE_ERROR_TOOLS_CHECK((dim == 3) || (dim == 2) || (dim == 1),
                                     "The dimension of the deformation gradient is " + std::to_string(dim) +
                                         " but must be 1, 2, or 3");

        E    = floatVector(dim * dim, 0);
        dEdF = floatVector(dim * dim * dim * dim, 0);

        if (dim == 3) {
            computeGreenLagrangeStrain<3>(std::begin(deformationGradient), std::end(deformationGradient), std::begin(E),
                                          std::end(E), std::begin(dEdF), std::end(dEdF));

        } else if (dim == 2) {
            computeGreenLagrangeStrain<2>(std::begin(deformationGradient), std::end(deformationGradient), std::begin(E),
                                          std::end(E), std::begin(dEdF), std::end(dEdF));

        } else if (dim == 1) {
            computeGreenLagrangeStrain<1>(std::begin(deformationGradient), std::end(deformationGradient), std::begin(E),
                                          std::end(E), std::begin(dEdF), std::end(dEdF));
        }

        return;
    }

    void computeDGreenLagrangeStrainDF(const floatVector &deformationGradient, floatMatrix &dEdF) {
        /*!
         * Compute the derivative of the Green-Lagrange strain ( \f$E\f$ )w.r.t. the deformation gradient ( \f$F\f$ ).
         *
         * \f$\frac{\partial E_{IJ}}{\partial F_{kK}} = 0.5 ( \delta_{IK} F_{kJ} + F_{kI} \delta_{JK})\f$
         *
         * Where \f$F\f$ is the deformation gradient and \f$\delta\f$ is the kronecker delta.
         *
         * \param &deformationGradient: A reference to the deformation gradient ( \f$F\f$ ).
         * \param &dEdF: The resulting gradient ( \f$\frac{\partial E}{\partial F}\f$ ).
         *
         * The deformation gradient is organized as  F11, F12, F13, F21, F22, F23, F31, F32, F33
         */

        floatVector _dEdF;

        TARDIGRADE_ERROR_TOOLS_CATCH(computeDGreenLagrangeStrainDF(deformationGradient, _dEdF));

        dEdF = tardigradeVectorTools::inflate(_dEdF, deformationGradient.size(), deformationGradient.size());

        return;
    }

    void computeDGreenLagrangeStrainDF(const floatVector &deformationGradient, floatVector &dEdF) {
        /*!
         * Compute the derivative of the Green-Lagrange strain ( \f$E\f$ )w.r.t. the deformation gradient ( \f$F\f$ ).
         *
         * \f$\frac{\partial E_{IJ}}{\partial F_{kK}} = 0.5 ( \delta_{IK} F_{kJ} + F_{kI} \delta_{JK})\f$
         *
         * Where \f$F\f$ is the deformation gradient and \f$\delta\f$ is the kronecker delta.
         *
         * \param &deformationGradient: A reference to the deformation gradient ( \f$F\f$ ).
         * \param &dEdF: The resulting gradient ( \f$\frac{\partial E}{\partial F}\f$ ).
         *
         * The deformation gradient is organized as  F11, F12, F13, F21, F22, F23, F31, F32, F33
         */

        const unsigned int dim = (unsigned int)std::pow(deformationGradient.size(), 0.5);

        TARDIGRADE_ERROR_TOOLS_CHECK((dim == 3) || (dim == 2) || (dim == 1),
                                     "The dimension of the deformation gradient is " + std::to_string(dim) +
                                         " but must be 1, 2, or 3");

        dEdF = floatVector(dim * dim * dim * dim, 0);

        if (dim == 3) {
            computeDGreenLagrangeStrainDF<3>(std::begin(deformationGradient), std::end(deformationGradient),
                                             std::begin(dEdF), std::end(dEdF));

        } else if (dim == 2) {
            computeDGreenLagrangeStrainDF<2>(std::begin(deformationGradient), std::end(deformationGradient),
                                             std::begin(dEdF), std::end(dEdF));

        } else if (dim == 1) {
            computeDGreenLagrangeStrainDF<1>(std::begin(deformationGradient), std::end(deformationGradient),
                                             std::begin(dEdF), std::end(dEdF));
        }

        return;
    }

    void decomposeGreenLagrangeStrain(const floatVector &E, floatVector &Ebar, floatType &J) {
        /*!
         * Decompose the Green-Lagrange strain tensor ( \f$E\f$ ) into isochoric ( \f$\bar{E}\f$ ) and volumetric (
         * \f$J\f$ ) parts where
         *
         * \f$J = det(F) = sqrt(det(2*E + I))\f$
         *
         * \f$\bar{E}_{IJ} = 0.5*((1/(J**(2/3))) F_{iI} F_{iJ} - I_{IJ}) = (1/(J**(2/3)))*E_{IJ} + 0.5(1/(J**(2/3)) -
         * 1)*I_{IJ}\f$
         *
         * \param &E: The Green-Lagrange strain tensor ( \f$E\f$ )
         * \param &Ebar: The isochoric Green-Lagrange strain tensor ( \f$\bar{E}\f$ ).
         *     format = E11, E12, E13, E21, E22, E23, E31, E32, E33
         * \param &J: The Jacobian of deformation ( \f$J\f$ )
         */

        const unsigned int dim = (unsigned int)std::pow(E.size(), 0.5);

        TARDIGRADE_ERROR_TOOLS_CHECK((dim == 3) || (dim == 2) || (dim == 1),
                                     "The dimension of the deformation gradient is " + std::to_string(dim) +
                                         " but must be 1, 2, or 3");

        Ebar = floatVector(dim * dim, 0);

        if (dim == 3) {
            decomposeGreenLagrangeStrain<3>(std::begin(E), std::end(E), std::begin(Ebar), std::end(Ebar), J);

        } else if (dim == 2) {
            decomposeGreenLagrangeStrain<2>(std::begin(E), std::end(E), std::begin(Ebar), std::end(Ebar), J);

        } else if (dim == 1) {
            decomposeGreenLagrangeStrain<1>(std::begin(E), std::end(E), std::begin(Ebar), std::end(Ebar), J);
        }

        return;
    }

    void decomposeGreenLagrangeStrain(const floatVector &E, floatVector &Ebar, floatType &J, floatVector &dEbardE,
                                      floatVector &dJdE) {
        /*!
         * Decompose the Green-Lagrange strain tensor ( \f$E\f$ ) into isochoric ( \f$\bar{E}\f$ ) and volumetric (
         * \f$J\f$ ) parts where
         *
         * \f$J = det(F) = sqrt(det(2*E + I))\f$
         *
         * \f$\bar{E}_{IJ} = 0.5*((1/(J**(2/3))) F_{iI} F_{iJ} - I_{IJ}) = (1/(J**(2/3)))*E_{IJ} + 0.5(1/(J**(2/3)) -
         * 1)*I_{IJ}\f$
         *
         * \param &E: The Green-Lagrange strain tensor ( \f$E\f$ )
         * \param &Ebar: The isochoric Green-Lagrange strain tensor ( \f$\bar{E}\f$ ).
         *     format = E11, E12, E13, E21, E22, E23, E31, E32, E33
         * \param &J: The Jacobian of deformation ( \f$J\f$ )
         * \param &dEbardE: The derivative of the isochoric Green-Lagrange strain
         *     tensor w.r.t. the total strain tensor ( \f$\frac{\partial \bar{E}}{\partial E}\f$ ).
         * \param &dJdE: The derivative of the jacobian of deformation w.r.t. the
         *     Green-Lagrange strain tensor ( \f$\frac{\partial J}{\partial E}\f$ ).
         */

        const unsigned int dim = (unsigned int)std::pow(E.size(), 0.5);

        TARDIGRADE_ERROR_TOOLS_CHECK((dim == 3) || (dim == 2) || (dim == 1),
                                     "The dimension of the deformation gradient is " + std::to_string(dim) +
                                         " but must be 1, 2, or 3");

        Ebar    = floatVector(dim * dim, 0);
        dEbardE = floatVector(dim * dim * dim * dim, 0);
        dJdE    = floatVector(dim * dim, 0);

        if (dim == 3) {
            decomposeGreenLagrangeStrain<3>(std::begin(E), std::end(E), std::begin(Ebar), std::end(Ebar), J,
                                            std::begin(dEbardE), std::end(dEbardE), std::begin(dJdE), std::end(dJdE));

        } else if (dim == 2) {
            decomposeGreenLagrangeStrain<2>(std::begin(E), std::end(E), std::begin(Ebar), std::end(Ebar), J,
                                            std::begin(dEbardE), std::end(dEbardE), std::begin(dJdE), std::end(dJdE));

        } else if (dim == 1) {
            decomposeGreenLagrangeStrain<1>(std::begin(E), std::end(E), std::begin(Ebar), std::end(Ebar), J,
                                            std::begin(dEbardE), std::end(dEbardE), std::begin(dJdE), std::end(dJdE));
        }

        return;
    }

    void decomposeGreenLagrangeStrain(const floatVector &E, floatVector &Ebar, floatType &J, floatMatrix &dEbardE,
                                      floatVector &dJdE) {
        /*!
         * Decompose the Green-Lagrange strain tensor ( \f$E\f$ ) into isochoric ( \f$\bar{E}\f$ ) and volumetric (
         * \f$J\f$ ) parts where
         *
         * \f$J = det(F) = sqrt(det(2*E + I))\f$
         *
         * \f$\bar{E}_{IJ} = 0.5*((1/(J**(2/3))) F_{iI} F_{iJ} - I_{IJ}) = (1/(J**(2/3)))*E_{IJ} + 0.5(1/(J**(2/3)) -
         * 1)*I_{IJ}\f$
         *
         * \param &E: The Green-Lagrange strain tensor ( \f$E\f$ )
         * \param &Ebar: The isochoric Green-Lagrange strain tensor ( \f$\bar{E}\f$ ).
         *     format = E11, E12, E13, E21, E22, E23, E31, E32, E33
         * \param &J: The Jacobian of deformation ( \f$J\f$ )
         * \param &dEbardE: The derivative of the isochoric Green-Lagrange strain
         *     tensor w.r.t. the total strain tensor ( \f$\frac{\partial \bar{E}}{\partial E}\f$ ).
         * \param &dJdE: The derivative of the jacobian of deformation w.r.t. the
         *     Green-Lagrange strain tensor ( \f$\frac{\partial J}{\partial E}\f$ ).
         */

        constexpr unsigned int dim     = 3;
        constexpr unsigned int sot_dim = dim * dim;

        floatVector _dEbardE;

        TARDIGRADE_ERROR_TOOLS_CATCH(decomposeGreenLagrangeStrain(E, Ebar, J, _dEbardE, dJdE));

        dEbardE = tardigradeVectorTools::inflate(_dEbardE, sot_dim, sot_dim);

        return;
    }

    void mapPK2toCauchy(const floatVector &PK2Stress, const floatVector &deformationGradient,
                        floatVector &cauchyStress) {
        /*!
         * Map the PK2 stress ( \f$P^{II}\f$ ) to the current configuration resulting in the Cauchy stress (
         * \f$\sigma\f$ ).
         *
         * \f$\sigma_{ij} = (1/det(F)) F_{iI} P^{II}_{IJ} F_{jJ}\f$
         *
         * where \f$F\f$ is the deformation gradient
         *
         * \param &PK2Stress: The Second Piola-Kirchoff stress ( \f$P^{II}\f$ )
         * \param &deformationGradient: The total deformation gradient ( \f$F\f$ ).
         * \param &cauchyStress: The Cauchy stress (\f$\sigma\f$ ).
         */

        const unsigned int dim = (unsigned int)std::pow(deformationGradient.size(), 0.5);

        TARDIGRADE_ERROR_TOOLS_CHECK((dim == 3) || (dim == 2) || (dim == 1),
                                     "The dimension of the deformation gradient is " + std::to_string(dim) +
                                         " but must be 1, 2, or 3");

        cauchyStress = floatVector(dim * dim, 0);

        if (dim == 3) {
            mapPK2toCauchy<3>(std::begin(PK2Stress), std::end(PK2Stress), std::begin(deformationGradient),
                              std::end(deformationGradient), std::begin(cauchyStress), std::end(cauchyStress));

        } else if (dim == 2) {
            mapPK2toCauchy<2>(std::begin(PK2Stress), std::end(PK2Stress), std::begin(deformationGradient),
                              std::end(deformationGradient), std::begin(cauchyStress), std::end(cauchyStress));

        } else if (dim == 1) {
            mapPK2toCauchy<1>(std::begin(PK2Stress), std::end(PK2Stress), std::begin(deformationGradient),
                              std::end(deformationGradient), std::begin(cauchyStress), std::end(cauchyStress));
        }

        return;
    }

    void WLF(const floatType &temperature, const floatVector &WLFParameters, floatType &factor) {
        /*!
         * An implementation of the Williams-Landel-Ferry equation.
         *
         * \f$factor = 10**((-C_1*(T - T_r))/(C_2 + T - T_r))\f$
         *
         * where \f$T\f$ is the temperature, \f$T_r\f$ is the reference temperature, and \f$C_1\f$ and \f$C_2\f$ are
         * parameters
         *
         * \param &temperature: The temperature \f$T\f$
         * \param &WLFParameters: The parameters for the function [\f$T_r\f$, \f$C_1\f$, \f$C_2\f$]
         * \param &factor: The shift factor
         */

        TARDIGRADE_ERROR_TOOLS_CATCH(WLF(temperature, std::begin(WLFParameters), std::end(WLFParameters), factor));

        return;
    }

    void WLF(const floatType &temperature, const floatVector &WLFParameters, floatType &factor, floatType &dFactordT) {
        /*!
         * An implementation of the Williams-Landel-Ferry equation that also returns the gradient w.r.t. \f$T\f$
         *
         * \param &temperature: The temperature ( \f$T\f$ )
         * \param &WLFParameters: The parameters for the function [\f$T_r\f$, \f$C_1\f$, \f$C_2\f$]
         * \param &factor: The shift factor
         * \param &dFactordT: The derivative of the shift factor w.r.t. the temperature ( \f$\frac{\partial
         * factor}{\partial T}\f$ )
         */

        TARDIGRADE_ERROR_TOOLS_CATCH(WLF(temperature, std::begin(WLFParameters), std::end(WLFParameters), factor,
                                         dFactordT));

        return;
    }

    void computeDFDt(const floatVector &velocityGradient, const floatVector &deformationGradient, floatVector &DFDt) {
        /*!
         * Compute the total time derivative of the deformation gradient.
         *
         * \f$\dot{F}_{iI} = L_{ij} F_{jI}\f$
         *
         * \param &velocityGradient: The velocity gradient \f$L_{ij}\f$
         * \param &deformationGradient: The deformation gradient \f$F_{iI}\f$
         * \param &DFDt: The total time derivative of the deformation gradient
         */

        const unsigned int dim = (unsigned int)std::pow(deformationGradient.size(), 0.5);

        TARDIGRADE_ERROR_TOOLS_CHECK((dim == 3) || (dim == 2) || (dim == 1),
                                     "The dimension of the deformation gradient is " + std::to_string(dim) +
                                         " but must be 1, 2, or 3");

        DFDt = floatVector(dim * dim, 0);

        if (dim == 3) {
            computeDFDt<3>(std::begin(velocityGradient), std::end(velocityGradient), std::begin(deformationGradient),
                           std::end(deformationGradient), std::begin(DFDt), std::end(DFDt));

        } else if (dim == 2) {
            computeDFDt<2>(std::begin(velocityGradient), std::end(velocityGradient), std::begin(deformationGradient),
                           std::end(deformationGradient), std::begin(DFDt), std::end(DFDt));

        } else if (dim == 1) {
            computeDFDt<1>(std::begin(velocityGradient), std::end(velocityGradient), std::begin(deformationGradient),
                           std::end(deformationGradient), std::begin(DFDt), std::end(DFDt));
        }

        return;
    }

    void computeDFDt(const floatVector &velocityGradient, const floatVector &deformationGradient, floatVector &DFDt,
                     floatVector &dDFDtdL, floatVector &dDFDtdF) {
        /*!
         * Compute the total time derivative of the deformation gradient
         * and return the partial derivatives w.r.t. L and F.
         *
         * \f$\dot{F}_{iI} = L_{ij} F_{jI}\f$
         * \f$\frac{\partial \dot{F}_{iI}}{\partial L_{kl}} = \delta_{ik} F{lI}\f$
         * \f$\frac{\partial \dot{F}_{iI}}{\partial F_{kK}} = L_{ik} \delta{IK}\f$
         *
         * \param &velocityGradient: The velocity gradient \f$L_{ij}\f$
         * \param &deformationGradient: The deformation gradient \f$F_{iI}\f$
         * \param &DFDt: The total time derivative of the deformation gradient
         * \param &dDFDtdL: The derivative of the total time derivative of the deformation gradient
         *     with respect to the velocity gradient.
         * \param &dDFDtdF: The derivative of the total time derivative of the deformation gradient
         *     with respect to the deformation gradient.
         */

        const unsigned int dim = (unsigned int)std::pow(deformationGradient.size(), 0.5);

        TARDIGRADE_ERROR_TOOLS_CHECK((dim == 3) || (dim == 2) || (dim == 1),
                                     "The dimension of the deformation gradient is " + std::to_string(dim) +
                                         " but must be 1, 2, or 3");

        DFDt    = floatVector(dim * dim, 0);
        dDFDtdL = floatVector(dim * dim * dim * dim, 0);
        dDFDtdF = floatVector(dim * dim * dim * dim, 0);

        if (dim == 3) {
            computeDFDt<3>(std::begin(velocityGradient), std::end(velocityGradient), std::begin(deformationGradient),
                           std::end(deformationGradient), std::begin(DFDt), std::end(DFDt), std::begin(dDFDtdL),
                           std::end(dDFDtdL), std::begin(dDFDtdF), std::end(dDFDtdF));

        } else if (dim == 2) {
            computeDFDt<2>(std::begin(velocityGradient), std::end(velocityGradient), std::begin(deformationGradient),
                           std::end(deformationGradient), std::begin(DFDt), std::end(DFDt), std::begin(dDFDtdL),
                           std::end(dDFDtdL), std::begin(dDFDtdF), std::end(dDFDtdF));

        } else if (dim == 1) {
            computeDFDt<1>(std::begin(velocityGradient), std::end(velocityGradient), std::begin(deformationGradient),
                           std::end(deformationGradient), std::begin(DFDt), std::end(DFDt), std::begin(dDFDtdL),
                           std::end(dDFDtdL), std::begin(dDFDtdF), std::end(dDFDtdF));
        }

        return;
    }

    void computeDFDt(const floatVector &velocityGradient, const floatVector &deformationGradient, floatVector &DFDt,
                     floatMatrix &dDFDtdL, floatMatrix &dDFDtdF) {
        /*!
         * Compute the total time derivative of the deformation gradient
         * and return the partial derivatives w.r.t. L and F.
         *
         * \f$\dot{F}_{iI} = L_{ij} F_{jI}\f$
         * \f$\frac{\partial \dot{F}_{iI}}{\partial L_{kl}} = \delta_{ik} F{lI}\f$
         * \f$\frac{\partial \dot{F}_{iI}}{\partial F_{kK}} = L_{ik} \delta{IK}\f$
         *
         * \param &velocityGradient: The velocity gradient \f$L_{ij}\f$
         * \param &deformationGradient: The deformation gradient \f$F_{iI}\f$
         * \param &DFDt: The total time derivative of the deformation gradient
         * \param &dDFDtdL: The derivative of the total time derivative of the deformation gradient
         *     with respect to the velocity gradient.
         * \param &dDFDtdF: The derivative of the total time derivative of the deformation gradient
         *     with respect to the deformation gradient.
         */

        // Assume 3D
        constexpr unsigned int dim     = 3;
        constexpr unsigned int sot_dim = dim * dim;

        floatVector _dDFDtdL, _dDFDtdF;

        TARDIGRADE_ERROR_TOOLS_CATCH(computeDFDt(velocityGradient, deformationGradient, DFDt, _dDFDtdL, _dDFDtdF));

        dDFDtdL = tardigradeVectorTools::inflate(_dDFDtdL, sot_dim, sot_dim);

        dDFDtdF = tardigradeVectorTools::inflate(_dDFDtdF, sot_dim, sot_dim);

        return;
    }

    void midpointEvolution(const floatType &Dt, const floatVector &Ap, const floatVector &DApDt,
                           const floatVector &DADt, floatVector &dA, floatVector &A, const floatVector &alpha) {
        /*!
         * Perform midpoint rule based evolution of a vector.
         *
         * alpha=0 (implicit)
         *
         * alpha=1 (explicit)
         *
         * \param &Dt: The change in time.
         * \param &Ap: The previous value of the vector
         * \param &DApDt: The previous time rate of change of the vector.
         * \param &DADt: The current time rate of change of the vector.
         * \param &dA: The change in the vector
         * \param &A: The current value of the vector.
         * \param &alpha: The integration parameter.
         */

        TARDIGRADE_ERROR_TOOLS_CHECK((Ap.size() == DApDt.size()) && (Ap.size() == DADt.size()),
                                     "The size of the previous value of the vector and the two rates are not equal");

        TARDIGRADE_ERROR_TOOLS_CHECK(Ap.size() == alpha.size(),
                                     "The size of the alpha vector is not the same size as the previous vector value");

        dA = floatVector(Ap.size(), 0);

        A = floatVector(Ap.size(), 0);

        midpointEvolution(Dt, std::begin(Ap), std::end(Ap), std::begin(DApDt), std::end(DApDt), std::begin(DADt),
                          std::end(DADt), std::begin(dA), std::end(dA), std::begin(A), std::end(A), std::begin(alpha),
                          std::end(alpha));

        return;
    }

    void midpointEvolutionFlatJ(const floatType &Dt, const floatVector &Ap, const floatVector &DApDt,
                                const floatVector &DADt, floatVector &dA, floatVector &A, floatVector &DADADt,
                                const floatVector &alpha) {
        /*!
         * Perform midpoint rule based evolution of a vector and return the jacobian.
         *
         * alpha=0 (implicit)
         *
         * alpha=1 (explicit)
         *
         * \param &Dt: The change in time.
         * \param &Ap: The previous value of the vector
         * \param &DApDt: The previous time rate of change of the vector.
         * \param &DADt: The current time rate of change of the vector.
         * \param &dA: The change in value of the vector.
         * \param &A: The current value of the vector.
         * \param &DADADt: The gradient of A w.r.t. the current rate of change.
         * \param &alpha: The integration parameter.
         */

        const unsigned int A_size = Ap.size();

        dA     = floatVector(A_size, 0);
        A      = floatVector(A_size, 0);
        DADADt = floatVector(A_size * A_size, 0);

        TARDIGRADE_ERROR_TOOLS_CATCH(midpointEvolution(Dt, std::begin(Ap), std::end(Ap), std::begin(DApDt),
                                                       std::end(DApDt), std::begin(DADt), std::end(DADt),
                                                       std::begin(dA), std::end(dA), std::begin(A), std::end(A),
                                                       std::begin(DADADt), std::end(DADADt), std::begin(alpha),
                                                       std::end(alpha)))

        return;
    }

    void midpointEvolution(const floatType &Dt, const floatVector &Ap, const floatVector &DApDt,
                           const floatVector &DADt, floatVector &dA, floatVector &A, floatMatrix &DADADt,
                           const floatVector &alpha) {
        /*!
         * Perform midpoint rule based evolution of a vector and return the jacobian.
         *
         * alpha=0 (implicit)
         *
         * alpha=1 (explicit)
         *
         * \param &Dt: The change in time.
         * \param &Ap: The previous value of the vector
         * \param &DApDt: The previous time rate of change of the vector.
         * \param &DADt: The current time rate of change of the vector.
         * \param &dA: The change in value of the vector.
         * \param &A: The current value of the vector.
         * \param &DADADt: The gradient of A w.r.t. the current rate of change.
         * \param &alpha: The integration parameter.
         */

        floatVector _DADADt;

        TARDIGRADE_ERROR_TOOLS_CATCH(midpointEvolutionFlatJ(Dt, Ap, DApDt, DADt, dA, A, _DADADt, alpha))

        DADADt = tardigradeVectorTools::inflate(_DADADt, A.size(), A.size());

        return;
    }

    void midpointEvolutionFlatJ(const floatType &Dt, const floatVector &Ap, const floatVector &DApDt,
                                const floatVector &DADt, floatVector &dA, floatVector &A, floatVector &DADADt,
                                floatVector &DADADtp, const floatVector &alpha) {
        /*!
         * Perform midpoint rule based evolution of a vector and return the jacobian.
         *
         * alpha=0 (implicit)
         *
         * alpha=1 (explicit)
         *
         * Note that the gradient of A w.r.t. Ap is identity and the gradient of dA w.r.t. Ap is zero
         *
         * \param &Dt: The change in time.
         * \param &Ap: The previous value of the vector
         * \param &DApDt: The previous time rate of change of the vector.
         * \param &DADt: The current time rate of change of the vector.
         * \param &dA: The change in value of the vector.
         * \param &A: The current value of the vector.
         * \param &DADADt: The gradient of A w.r.t. the current rate of change.
         * \param &DADADtp: The gradient of A w.r.t. the previous rate of change.
         * \param &alpha: The integration parameter.
         */

        const unsigned int A_size = Ap.size();

        dA      = floatVector(A_size, 0);
        A       = floatVector(A_size, 0);
        DADADt  = floatVector(A_size * A_size, 0);
        DADADtp = floatVector(A_size * A_size, 0);

        TARDIGRADE_ERROR_TOOLS_CATCH(midpointEvolution(Dt, std::begin(Ap), std::end(Ap), std::begin(DApDt),
                                                       std::end(DApDt), std::begin(DADt), std::end(DADt),
                                                       std::begin(dA), std::end(dA), std::begin(A), std::end(A),
                                                       std::begin(DADADt), std::end(DADADt), std::begin(DADADtp),
                                                       std::end(DADADtp), std::begin(alpha), std::end(alpha)))

        return;
    }

    void midpointEvolution(const floatType &Dt, const floatVector &Ap, const floatVector &DApDt,
                           const floatVector &DADt, floatVector &dA, floatVector &A, floatMatrix &DADADt,
                           floatMatrix &DADADtp, const floatVector &alpha) {
        /*!
         * Perform midpoint rule based evolution of a vector and return the jacobian.
         *
         * alpha=0 (implicit)
         *
         * alpha=1 (explicit)
         *
         * Note that the gradient of A w.r.t. Ap is identity and the gradient of dA w.r.t. Ap is zero
         *
         * \param &Dt: The change in time.
         * \param &Ap: The previous value of the vector
         * \param &DApDt: The previous time rate of change of the vector.
         * \param &DADt: The current time rate of change of the vector.
         * \param &dA: The change in value of the vector.
         * \param &A: The current value of the vector.
         * \param &DADADt: The gradient of A w.r.t. the current rate of change.
         * \param &DADADtp: The gradient of A w.r.t. the previous rate of change.
         * \param &alpha: The integration parameter.
         */

        floatVector _DADADt, _DADADtp;

        TARDIGRADE_ERROR_TOOLS_CATCH(midpointEvolutionFlatJ(Dt, Ap, DApDt, DADt, dA, A, _DADADt, _DADADtp, alpha));

        DADADt = tardigradeVectorTools::inflate(_DADADt, A.size(), A.size());

        DADADtp = tardigradeVectorTools::inflate(_DADADtp, A.size(), A.size());

        return;
    }

    void midpointEvolution(const floatType &Dt, const floatVector &Ap, const floatVector &DApDt,
                           const floatVector &DADt, floatVector &dA, floatVector &A, const floatType alpha) {
        /*!
         * Perform midpoint rule based evolution of a vector. Defaults to the trapezoidal rule.
         *
         * alpha=0 (implicit)
         *
         * alpha=1 (explicit)
         *
         * \param &Dt: The change in time.
         * \param &Ap: The previous value of the vector
         * \param &DApDt: The previous time rate of change of the vector.
         * \param &DADt: The current time rate of change of the vector.
         * \param &dA: The change in the value of the vector.
         * \param &A: The current value of the vector.
         * \param alpha: The integration parameter.
         */

        dA = floatVector(Ap.size(), 0);
        A  = floatVector(Ap.size(), 0);

        TARDIGRADE_ERROR_TOOLS_CATCH(midpointEvolution(Dt, std::begin(Ap), std::end(Ap), std::begin(DApDt),
                                                       std::end(DApDt), std::begin(DADt), std::end(DADt),
                                                       std::begin(dA), std::end(dA), std::begin(A), std::end(A), alpha))

        return;
    }

    void midpointEvolutionFlatJ(const floatType &Dt, const floatVector &Ap, const floatVector &DApDt,
                                const floatVector &DADt, floatVector &dA, floatVector &A, floatVector &DADADt,
                                const floatType alpha) {
        /*!
         * Perform midpoint rule based evolution of a vector. Defaults to the trapezoidal rule.
         *
         * alpha=0 (implicit)
         *
         * alpha=1 (explicit)
         *
         * \param &Dt: The change in time.
         * \param &Ap: The previous value of the vector
         * \param &DApDt: The previous time rate of change of the vector.
         * \param *DADt: The current time rate of change of the vector.
         * \param &dA: The change in the vector
         * \param &A: The current value of the vector.
         * \param &DADADt: The derivative of the vector w.r.t. the rate of change of the vector.
         * \param alpha: The integration parameter.
         */

        dA     = floatVector(Ap.size(), 0);
        A      = floatVector(Ap.size(), 0);
        DADADt = floatVector(Ap.size() * Ap.size(), 0);

        TARDIGRADE_ERROR_TOOLS_CATCH(midpointEvolution(Dt, std::begin(Ap), std::end(Ap), std::begin(DApDt),
                                                       std::end(DApDt), std::begin(DADt), std::end(DADt),
                                                       std::begin(dA), std::end(dA), std::begin(A), std::end(A),
                                                       std::begin(DADADt), std::end(DADADt), alpha))

        return;
    }

    void midpointEvolutionFlatJ(const floatType &Dt, const floatVector &Ap, const floatVector &DApDt,
                                const floatVector &DADt, floatVector &dA, floatVector &A, floatVector &DADADt,
                                floatVector &DADADtp, const floatType alpha) {
        /*!
         * Perform midpoint rule based evolution of a vector. Defaults to the trapezoidal rule.
         *
         * alpha=0 (implicit)
         *
         * alpha=1 (explicit)
         *
         * Note that the gradient of A w.r.t. Ap is identity and the gradient of dA w.r.t. Ap is zero
         *
         * \param &Dt: The change in time.
         * \param &Ap: The previous value of the vector
         * \param &DApDt: The previous time rate of change of the vector.
         * \param *DADt: The current time rate of change of the vector.
         * \param &dA: The change in the vector
         * \param &A: The current value of the vector.
         * \param &DADADt: The derivative of the vector w.r.t. the rate of change of the vector.
         * \param &DADADtp: The derivative of the vector w.r.t. the previous rate of change of the vector.
         * \param alpha: The integration parameter.
         */

        dA      = floatVector(Ap.size(), 0);
        A       = floatVector(Ap.size(), 0);
        DADADt  = floatVector(Ap.size() * Ap.size(), 0);
        DADADtp = floatVector(Ap.size() * Ap.size(), 0);

        TARDIGRADE_ERROR_TOOLS_CATCH(midpointEvolution(Dt, std::begin(Ap), std::end(Ap), std::begin(DApDt),
                                                       std::end(DApDt), std::begin(DADt), std::end(DADt),
                                                       std::begin(dA), std::end(dA), std::begin(A), std::end(A),
                                                       std::begin(DADADt), std::end(DADADt), std::begin(DADADtp),
                                                       std::end(DADADtp), alpha))

        return;
    }

    void midpointEvolution(const floatType &Dt, const floatVector &Ap, const floatVector &DApDt,
                           const floatVector &DADt, floatVector &dA, floatVector &A, floatMatrix &DADADt,
                           const floatType alpha) {
        /*!
         * Perform midpoint rule based evolution of a vector. Defaults to the trapezoidal rule.
         *
         * alpha=0 (implicit)
         *
         * alpha=1 (explicit)
         *
         * \param &Dt: The change in time.
         * \param &Ap: The previous value of the vector
         * \param &DApDt: The previous time rate of change of the vector.
         * \param *DADt: The current time rate of change of the vector.
         * \param &dA: The change in the vector
         * \param &A: The current value of the vector.
         * \param &DADADt: The derivative of the vector w.r.t. the rate of change of the vector.
         * \param alpha: The integration parameter.
         */

        return midpointEvolution(Dt, Ap, DApDt, DADt, dA, A, DADADt, alpha * floatVector(Ap.size(), 1));
    }

    void midpointEvolution(const floatType &Dt, const floatVector &Ap, const floatVector &DApDt,
                           const floatVector &DADt, floatVector &dA, floatVector &A, floatMatrix &DADADt,
                           floatMatrix &DADADtp, const floatType alpha) {
        /*!
         * Perform midpoint rule based evolution of a vector. Defaults to the trapezoidal rule.
         *
         * alpha=0 (implicit)
         *
         * alpha=1 (explicit)
         *
         * Note that the gradient of A w.r.t. Ap is identity and the gradient of dA w.r.t. Ap is zero
         *
         * \param &Dt: The change in time.
         * \param &Ap: The previous value of the vector
         * \param &DApDt: The previous time rate of change of the vector.
         * \param *DADt: The current time rate of change of the vector.
         * \param &dA: The change in the vector
         * \param &A: The current value of the vector.
         * \param &DADADt: The derivative of the vector w.r.t. the rate of change of the vector.
         * \param &DADADtp: The derivative of the vector w.r.t. the previous rate of change of the vector.
         * \param alpha: The integration parameter.
         */

        return midpointEvolution(Dt, Ap, DApDt, DADt, dA, A, DADADt, DADADtp, alpha * floatVector(Ap.size(), 1));
    }

    void evolveF(const floatType &Dt, const floatVector &previousDeformationGradient, const floatVector &Lp,
                 const floatVector &L, floatVector &dF, floatVector &deformationGradient, const floatType alpha,
                 const unsigned int mode) {
        /*!
         * Evolve the deformation gradient ( F ) using the midpoint integration method.
         *
         * mode 1:
         * \f$F_{iI}^{t + 1} = \left[\delta_{ij} - \Delta t \left(1 - \alpha \right) L_{ij}^{t+1} \right]^{-1}
         * \left[F_{iI}^{t} + \Delta t \alpha \dot{F}_{iI}^{t} \right]\f$
         *
         * mode 2:
         * \f$F_{iI}^{t + 1} = \left[F_{iJ}^{t} + \Delta t \alpha \dot{F}_{iJ}^{t} \right] \left[\delta_{IJ} - \Delta T
         * \left( 1- \alpha \right) L_{IJ}^{t+1} \right]^{-1}\f$
         *
         * \param &Dt: The change in time.
         * \param &previousDeformationGradient: The previous value of the deformation gradient
         * \param &Lp: The previous velocity gradient in the current configuration (mode 1) or
         *     reference configuration (mode 2).
         * \param &L: The current velocity gradient in the current configuration (mode 1) or
         *     reference configuration (mode 2).
         * \param &dF: The change in the deformation gradient \f$\Delta \bf{F}\f$ such that \f$F_{iI}^{t+1} = F_{iI}^t +
         * \Delta F_{iI}\f$ \param &deformationGradient: The computed current deformation gradient. \param alpha: The
         * integration parameter. \param mode: The mode of the ODE. Whether the velocity gradient is known in the
         *     current (mode 1) or reference (mode 2) configuration.
         */

        const unsigned int dim = (unsigned int)std::pow(previousDeformationGradient.size(), 0.5);

        dF                  = floatVector(dim * dim, 0);
        deformationGradient = floatVector(dim * dim, 0);

        TARDIGRADE_ERROR_TOOLS_CHECK((dim == 3) || (dim == 2) || (dim == 1),
                                     "The dimension of the deformation gradient is " + std::to_string(dim) +
                                         " but must be 1, 2, or 3");

        TARDIGRADE_ERROR_TOOLS_CHECK((mode == 1) || (mode == 2), "The mode of evolution " + std::to_string(mode) +
                                                                     " is not recognized. It must be 1 or 2.");

        if (mode == 1) {
            if (dim == 3) {
                evolveF<3, 1>(Dt, std::begin(previousDeformationGradient), std::end(previousDeformationGradient),
                              std::begin(Lp), std::end(Lp), std::begin(L), std::end(L), std::begin(dF), std::end(dF),
                              std::begin(deformationGradient), std::end(deformationGradient), alpha);
            } else if (dim == 2) {
                evolveF<2, 1>(Dt, std::begin(previousDeformationGradient), std::end(previousDeformationGradient),
                              std::begin(Lp), std::end(Lp), std::begin(L), std::end(L), std::begin(dF), std::end(dF),
                              std::begin(deformationGradient), std::end(deformationGradient), alpha);
            } else if (dim == 1) {
                evolveF<1, 1>(Dt, std::begin(previousDeformationGradient), std::end(previousDeformationGradient),
                              std::begin(Lp), std::end(Lp), std::begin(L), std::end(L), std::begin(dF), std::end(dF),
                              std::begin(deformationGradient), std::end(deformationGradient), alpha);
            }

        } else if (mode == 2) {
            if (dim == 3) {
                evolveF<3, 2>(Dt, std::begin(previousDeformationGradient), std::end(previousDeformationGradient),
                              std::begin(Lp), std::end(Lp), std::begin(L), std::end(L), std::begin(dF), std::end(dF),
                              std::begin(deformationGradient), std::end(deformationGradient), alpha);
            } else if (dim == 2) {
                evolveF<2, 2>(Dt, std::begin(previousDeformationGradient), std::end(previousDeformationGradient),
                              std::begin(Lp), std::end(Lp), std::begin(L), std::end(L), std::begin(dF), std::end(dF),
                              std::begin(deformationGradient), std::end(deformationGradient), alpha);
            } else if (dim == 1) {
                evolveF<1, 2>(Dt, std::begin(previousDeformationGradient), std::end(previousDeformationGradient),
                              std::begin(Lp), std::end(Lp), std::begin(L), std::end(L), std::begin(dF), std::end(dF),
                              std::begin(deformationGradient), std::end(deformationGradient), alpha);
            }
        }

        return;
    }

    void evolveF(const floatType &Dt, const floatVector &previousDeformationGradient, const floatVector &Lp,
                 const floatVector &L, floatVector &deformationGradient, const floatType alpha,
                 const unsigned int mode) {
        /*!
         * Evolve the deformation gradient ( F ) using the midpoint integration method.
         *
         * mode 1:
         * \f$F_{iI}^{t + 1} = \left[\delta_{ij} - \Delta t \left(1 - \alpha \right) L_{ij}^{t+1} \right]^{-1}
         * \left[F_{iI}^{t} + \Delta t \alpha \dot{F}_{iI}^{t} \right]\f$
         *
         * mode 2:
         * \f$F_{iI}^{t + 1} = \left[F_{iJ}^{t} + \Delta t \alpha \dot{F}_{iJ}^{t} \right] \left[\delta_{IJ} - \Delta T
         * \left( 1- \alpha \right) L_{IJ}^{t+1} \right]^{-1}\f$
         *
         * \param &Dt: The change in time.
         * \param &previousDeformationGradient: The previous value of the deformation gradient
         * \param &Lp: The previous velocity gradient in the current configuration (mode 1) or
         *     reference configuration (mode 2).
         * \param &L: The current velocity gradient in the current configuration (mode 1) or
         *     reference configuration (mode 2).
         * \param &deformationGradient: The computed current deformation gradient.
         * \param alpha: The integration parameter.
         * \param mode: The mode of the ODE. Whether the velocity gradient is known in the
         *     current (mode 1) or reference (mode 2) configuration.
         */

        const unsigned int dim = (unsigned int)std::pow(previousDeformationGradient.size(), 0.5);

        floatVector dF(dim * dim, 0);
        deformationGradient = floatVector(dim * dim, 0);

        TARDIGRADE_ERROR_TOOLS_CHECK((dim == 3) || (dim == 2) || (dim == 1),
                                     "The dimension of the deformation gradient is " + std::to_string(dim) +
                                         " but must be 1, 2, or 3");

        TARDIGRADE_ERROR_TOOLS_CHECK((mode == 1) || (mode == 2), "The mode of evolution " + std::to_string(mode) +
                                                                     " is not recognized. It must be 1 or 2.");

        if (mode == 1) {
            if (dim == 3) {
                evolveF<3, 1>(Dt, std::begin(previousDeformationGradient), std::end(previousDeformationGradient),
                              std::begin(Lp), std::end(Lp), std::begin(L), std::end(L), std::begin(dF), std::end(dF),
                              std::begin(deformationGradient), std::end(deformationGradient), alpha);
            } else if (dim == 2) {
                evolveF<2, 1>(Dt, std::begin(previousDeformationGradient), std::end(previousDeformationGradient),
                              std::begin(Lp), std::end(Lp), std::begin(L), std::end(L), std::begin(dF), std::end(dF),
                              std::begin(deformationGradient), std::end(deformationGradient), alpha);
            } else if (dim == 1) {
                evolveF<1, 1>(Dt, std::begin(previousDeformationGradient), std::end(previousDeformationGradient),
                              std::begin(Lp), std::end(Lp), std::begin(L), std::end(L), std::begin(dF), std::end(dF),
                              std::begin(deformationGradient), std::end(deformationGradient), alpha);
            }

        } else if (mode == 2) {
            if (dim == 3) {
                evolveF<3, 2>(Dt, std::begin(previousDeformationGradient), std::end(previousDeformationGradient),
                              std::begin(Lp), std::end(Lp), std::begin(L), std::end(L), std::begin(dF), std::end(dF),
                              std::begin(deformationGradient), std::end(deformationGradient), alpha);
            } else if (dim == 2) {
                evolveF<2, 2>(Dt, std::begin(previousDeformationGradient), std::end(previousDeformationGradient),
                              std::begin(Lp), std::end(Lp), std::begin(L), std::end(L), std::begin(dF), std::end(dF),
                              std::begin(deformationGradient), std::end(deformationGradient), alpha);
            } else if (dim == 1) {
                evolveF<1, 2>(Dt, std::begin(previousDeformationGradient), std::end(previousDeformationGradient),
                              std::begin(Lp), std::end(Lp), std::begin(L), std::end(L), std::begin(dF), std::end(dF),
                              std::begin(deformationGradient), std::end(deformationGradient), alpha);
            }
        }

        return;
    }

    void evolveF(const floatType &Dt, const floatVector &previousDeformationGradient, const floatVector &Lp,
                 const floatVector &L, floatVector &dF, floatVector &deformationGradient, floatMatrix &dFdL,
                 const floatType alpha, const unsigned int mode) {
        /*!
         * Evolve the deformation gradient ( F ) using the midpoint integration method and return the jacobian w.r.t. L.
         *
         * mode 1:
         * \f$F_{iI}^{t + 1} = \left[\delta_{ij} - \Delta t \left(1 - \alpha \right) L_{ij}^{t+1} \right]^{-1}
         * \left[F_{iI}^{t} + \Delta t \alpha \dot{F}_{iI}^{t} \right]\f$ \f$\frac{\partial F_{jI}^{t + 1}}{\partial
         * L_{kl}^{t+1}} = \left[\delta_{kj} - \Delta t \left(1 - \alpha\right) L_{kj}\right]^{-1} \Delta t \left(1 -
         * \alpha\right) F_{lI}^{t + 1}\f$
         *
         * mode 2:
         * \f$F_{iI}^{t + 1} = \left[F_{iJ}^{t} + \Delta t \alpha \dot{F}_{iJ}^{t} \right] \left[\delta_{IJ} - \Delta T
         * \left( 1- \alpha \right) L_{IJ}^{t+1} \right]^{-1}\f$ \f$\frac{\partial F_{iJ}^{t + 1}}{\partial L_{KL}} =
         * \Delta t (1 - \alpha) F_{iK}^{t + 1} \left[\delta_{JL} - \right ]\f$
         *
         * \param &Dt: The change in time.
         * \param &previousDeformationGradient: The previous value of the deformation gradient
         * \param &Lp: The previous velocity gradient.
         * \param &L: The current velocity gradient.
         * \param &dF: The change in the deformation gradient \f$\Delta \bf{F}\f$ such that \f$F_{iI}^{t+1} = F_{iI}^t +
         * \Delta F_{iI}\f$ \param &deformationGradient: The computed current deformation gradient. \param &dFdL: The
         * derivative of the deformation gradient w.r.t. the velocity gradient \param alpha: The integration parameter (
         * 0 for implicit, 1 for explicit ) \param mode: The form of the ODE. See above for details.
         */

        floatVector _dFdL;

        TARDIGRADE_ERROR_TOOLS_CATCH(evolveFFlatJ(Dt, previousDeformationGradient, Lp, L, dF, deformationGradient,
                                                  _dFdL, alpha, mode));

        dFdL = tardigradeVectorTools::inflate(_dFdL, deformationGradient.size(), deformationGradient.size());

        return;
    }

    void evolveFFlatJ(const floatType &Dt, const floatVector &previousDeformationGradient, const floatVector &Lp,
                      const floatVector &L, floatVector &dF, floatVector &deformationGradient, floatVector &dFdL,
                      const floatType alpha, const unsigned int mode) {
        /*!
         * Evolve the deformation gradient ( F ) using the midpoint integration method and return the jacobian w.r.t. L.
         *
         * mode 1:
         * \f$F_{iI}^{t + 1} = \left[\delta_{ij} - \Delta t \left(1 - \alpha \right) L_{ij}^{t+1} \right]^{-1}
         * \left[F_{iI}^{t} + \Delta t \alpha \dot{F}_{iI}^{t} \right]\f$ \f$\frac{\partial F_{jI}^{t + 1}}{\partial
         * L_{kl}^{t+1}} = \left[\delta_{kj} - \Delta t \left(1 - \alpha\right) L_{kj}\right]^{-1} \Delta t \left(1 -
         * \alpha\right) F_{lI}^{t + 1}\f$
         *
         * mode 2:
         * \f$F_{iI}^{t + 1} = \left[F_{iJ}^{t} + \Delta t \alpha \dot{F}_{iJ}^{t} \right] \left[\delta_{IJ} - \Delta T
         * \left( 1- \alpha \right) L_{IJ}^{t+1} \right]^{-1}\f$ \f$\frac{\partial F_{iJ}^{t + 1}}{\partial L_{KL}} =
         * \Delta t (1 - \alpha) F_{iK}^{t + 1} \left[\delta_{JL} - \right ]\f$
         *
         * \param &Dt: The change in time.
         * \param &previousDeformationGradient: The previous value of the deformation gradient
         * \param &Lp: The previous velocity gradient.
         * \param &L: The current velocity gradient.
         * \param &dF: The change in the deformation gradient \f$\Delta \bf{F}\f$ such that \f$F_{iI}^{t+1} = F_{iI}^t +
         * \Delta F_{iI}\f$ \param &deformationGradient: The computed current deformation gradient. \param &dFdL: The
         * derivative of the deformation gradient w.r.t. the velocity gradient \param alpha: The integration parameter (
         * 0 for implicit, 1 for explicit ) \param mode: The form of the ODE. See above for details.
         */

        const unsigned int dim = (unsigned int)std::pow(previousDeformationGradient.size(), 0.5);

        dF                  = floatVector(dim * dim, 0);
        deformationGradient = floatVector(dim * dim, 0);
        dFdL                = floatVector(dim * dim * dim * dim, 0);

        TARDIGRADE_ERROR_TOOLS_CHECK((dim == 3) || (dim == 2) || (dim == 1),
                                     "The dimension of the deformation gradient is " + std::to_string(dim) +
                                         " but must be 1, 2, or 3");

        TARDIGRADE_ERROR_TOOLS_CHECK((mode == 1) || (mode == 2), "The mode of evolution " + std::to_string(mode) +
                                                                     " is not recognized. It must be 1 or 2.");

        if (mode == 1) {
            if (dim == 3) {
                evolveF<3, 1>(Dt, std::begin(previousDeformationGradient), std::end(previousDeformationGradient),
                              std::begin(Lp), std::end(Lp), std::begin(L), std::end(L), std::begin(dF), std::end(dF),
                              std::begin(deformationGradient), std::end(deformationGradient), std::begin(dFdL),
                              std::end(dFdL), alpha);
            } else if (dim == 2) {
                evolveF<2, 1>(Dt, std::begin(previousDeformationGradient), std::end(previousDeformationGradient),
                              std::begin(Lp), std::end(Lp), std::begin(L), std::end(L), std::begin(dF), std::end(dF),
                              std::begin(deformationGradient), std::end(deformationGradient), std::begin(dFdL),
                              std::end(dFdL), alpha);
            } else if (dim == 1) {
                evolveF<1, 1>(Dt, std::begin(previousDeformationGradient), std::end(previousDeformationGradient),
                              std::begin(Lp), std::end(Lp), std::begin(L), std::end(L), std::begin(dF), std::end(dF),
                              std::begin(deformationGradient), std::end(deformationGradient), std::begin(dFdL),
                              std::end(dFdL), alpha);
            }

        } else if (mode == 2) {
            if (dim == 3) {
                evolveF<3, 2>(Dt, std::begin(previousDeformationGradient), std::end(previousDeformationGradient),
                              std::begin(Lp), std::end(Lp), std::begin(L), std::end(L), std::begin(dF), std::end(dF),
                              std::begin(deformationGradient), std::end(deformationGradient), std::begin(dFdL),
                              std::end(dFdL), alpha);
            } else if (dim == 2) {
                evolveF<2, 2>(Dt, std::begin(previousDeformationGradient), std::end(previousDeformationGradient),
                              std::begin(Lp), std::end(Lp), std::begin(L), std::end(L), std::begin(dF), std::end(dF),
                              std::begin(deformationGradient), std::end(deformationGradient), std::begin(dFdL),
                              std::end(dFdL), alpha);
            } else if (dim == 1) {
                evolveF<1, 2>(Dt, std::begin(previousDeformationGradient), std::end(previousDeformationGradient),
                              std::begin(Lp), std::end(Lp), std::begin(L), std::end(L), std::begin(dF), std::end(dF),
                              std::begin(deformationGradient), std::end(deformationGradient), std::begin(dFdL),
                              std::end(dFdL), alpha);
            }
        }
    }

    void evolveFFlatJ(const floatType &Dt, const floatVector &previousDeformationGradient, const floatVector &Lp,
                      const floatVector &L, floatVector &deformationGradient, floatVector &dFdL, const floatType alpha,
                      const unsigned int mode) {
        /*!
         * Evolve the deformation gradient ( F ) using the midpoint integration method and return the jacobian w.r.t. L.
         *
         * mode 1:
         * \f$F_{iI}^{t + 1} = \left[\delta_{ij} - \Delta t \left(1 - \alpha \right) L_{ij}^{t+1} \right]^{-1}
         * \left[F_{iI}^{t} + \Delta t \alpha \dot{F}_{iI}^{t} \right]\f$ \f$\frac{\partial F_{jI}^{t + 1}}{\partial
         * L_{kl}^{t+1}} = \left[\delta_{kj} - \Delta t \left(1 - \alpha\right) L_{kj}\right]^{-1} \Delta t \left(1 -
         * \alpha\right) F_{lI}^{t + 1}\f$
         *
         * mode 2:
         * \f$F_{iI}^{t + 1} = \left[F_{iJ}^{t} + \Delta t \alpha \dot{F}_{iJ}^{t} \right] \left[\delta_{IJ} - \Delta T
         * \left( 1- \alpha \right) L_{IJ}^{t+1} \right]^{-1}\f$ \f$\frac{\partial F_{iJ}^{t + 1}}{\partial L_{KL}} =
         * \Delta t (1 - \alpha) F_{iK}^{t + 1} \left[\delta_{JL} - \right ]\f$
         *
         * \param &Dt: The change in time.
         * \param &previousDeformationGradient: The previous value of the deformation gradient
         * \param &Lp: The previous velocity gradient.
         * \param &L: The current velocity gradient.
         * \param &deformationGradient: The computed current deformation gradient.
         * \param &dFdL: The derivative of the deformation gradient w.r.t. the velocity gradient
         * \param alpha: The integration parameter ( 0 for implicit, 1 for explicit )
         * \param mode: The form of the ODE. See above for details.
         */

        floatVector dF;

        return evolveFFlatJ(Dt, previousDeformationGradient, Lp, L, dF, deformationGradient, dFdL, alpha, mode);
    }

    void evolveF(const floatType &Dt, const floatVector &previousDeformationGradient, const floatVector &Lp,
                 const floatVector &L, floatVector &deformationGradient, floatMatrix &dFdL, const floatType alpha,
                 const unsigned int mode) {
        /*!
         * Evolve the deformation gradient ( F ) using the midpoint integration method and return the jacobian w.r.t. L.
         *
         * mode 1:
         * \f$F_{iI}^{t + 1} = \left[\delta_{ij} - \Delta t \left(1 - \alpha \right) L_{ij}^{t+1} \right]^{-1}
         * \left[F_{iI}^{t} + \Delta t \alpha \dot{F}_{iI}^{t} \right]\f$ \f$\frac{\partial F_{jI}^{t + 1}}{\partial
         * L_{kl}^{t+1}} = \left[\delta_{kj} - \Delta t \left(1 - \alpha\right) L_{kj}\right]^{-1} \Delta t \left(1 -
         * \alpha\right) F_{lI}^{t + 1}\f$
         *
         * mode 2:
         * \f$F_{iI}^{t + 1} = \left[F_{iJ}^{t} + \Delta t \alpha \dot{F}_{iJ}^{t} \right] \left[\delta_{IJ} - \Delta T
         * \left( 1- \alpha \right) L_{IJ}^{t+1} \right]^{-1}\f$ \f$\frac{\partial F_{iJ}^{t + 1}}{\partial L_{KL}} =
         * \Delta t (1 - \alpha) F_{iK}^{t + 1} \left[\delta_{JL} - \right ]\f$
         *
         * \param &Dt: The change in time.
         * \param &previousDeformationGradient: The previous value of the deformation gradient
         * \param &Lp: The previous velocity gradient.
         * \param &L: The current velocity gradient.
         * \param &deformationGradient: The computed current deformation gradient.
         * \param &dFdL: The derivative of the deformation gradient w.r.t. the velocity gradient
         * \param alpha: The integration parameter ( 0 for implicit, 1 for explicit )
         * \param mode: The form of the ODE. See above for details.
         */

        floatVector dF;

        return evolveF(Dt, previousDeformationGradient, Lp, L, dF, deformationGradient, dFdL, alpha, mode);
    }

    void evolveF(const floatType &Dt, const floatVector &previousDeformationGradient, const floatVector &Lp,
                 const floatVector &L, floatVector &dF, floatVector &deformationGradient, floatMatrix &dFdL,
                 floatMatrix &ddFdFp, floatMatrix &dFdFp, floatMatrix &dFdLp, const floatType alpha,
                 const unsigned int mode) {
        /*!
         * Evolve the deformation gradient ( F ) using the midpoint integration method and return the jacobian w.r.t. L.
         *
         * mode 1:
         * \f$F_{iI}^{t + 1} = \left[\delta_{ij} - \Delta t \left(1 - \alpha \right) L_{ij}^{t+1} \right]^{-1}
         * \left[F_{iI}^{t} + \Delta t \alpha \dot{F}_{iI}^{t} \right]\f$ \f$\frac{\partial F_{jI}^{t + 1}}{\partial
         * L_{kl}^{t+1}} = \left[\delta_{kj} - \Delta t \left(1 - \alpha\right) L_{kj}\right]^{-1} \Delta t \left(1 -
         * \alpha\right) F_{lI}^{t + 1}\f$
         *
         * mode 2:
         * \f$F_{iI}^{t + 1} = \left[F_{iJ}^{t} + \Delta t \alpha \dot{F}_{iJ}^{t} \right] \left[\delta_{IJ} - \Delta T
         * \left( 1- \alpha \right) L_{IJ}^{t+1} \right]^{-1}\f$ \f$\frac{\partial F_{iJ}^{t + 1}}{\partial L_{KL}} =
         * \Delta t (1 - \alpha) F_{iK}^{t + 1} \left[\delta_{JL} - \right ]\f$
         *
         * \param &Dt: The change in time.
         * \param &previousDeformationGradient: The previous value of the deformation gradient
         * \param &Lp: The previous velocity gradient.
         * \param &L: The current velocity gradient.
         * \param &dF: The change in the deformation gradient \f$\Delta \bf{F}\f$ such that \f$F_{iI}^{t+1} = F_{iI}^t +
         * \Delta F_{iI}\f$ \param &deformationGradient: The computed current deformation gradient. \param &dFdL: The
         * derivative of the deformation gradient w.r.t. the velocity gradient \param &ddFdFp: The derivative of the
         * change in the deformation gradient w.r.t. the previous deformation gradient \param &dFdFp: The derivative of
         * the deformation gradient w.r.t. the previous deformation gradient \param &dFdLp: The derivative of the
         * deformation gradient w.r.t. the previous velocity gradient \param alpha: The integration parameter ( 0 for
         * implicit, 1 for explicit ) \param mode: The form of the ODE. See above for details.
         */

        floatVector _dFdL;
        floatVector _ddFdFp;
        floatVector _dFdFp;
        floatVector _dFdLp;

        TARDIGRADE_ERROR_TOOLS_CATCH(evolveFFlatJ(Dt, previousDeformationGradient, Lp, L, dF, deformationGradient,
                                                  _dFdL, _ddFdFp, _dFdFp, _dFdLp, alpha, mode));

        dFdL   = tardigradeVectorTools::inflate(_dFdL, deformationGradient.size(), deformationGradient.size());
        ddFdFp = tardigradeVectorTools::inflate(_ddFdFp, deformationGradient.size(), deformationGradient.size());
        dFdFp  = tardigradeVectorTools::inflate(_dFdFp, deformationGradient.size(), deformationGradient.size());
        dFdLp  = tardigradeVectorTools::inflate(_dFdLp, deformationGradient.size(), deformationGradient.size());

        return;
    }

    void evolveFFlatJ(const floatType &Dt, const floatVector &previousDeformationGradient, const floatVector &Lp,
                      const floatVector &L, floatVector &dF, floatVector &deformationGradient, floatVector &dFdL,
                      floatVector &ddFdFp, floatVector &dFdFp, floatVector &dFdLp, const floatType alpha,
                      const unsigned int mode) {
        /*!
         * Evolve the deformation gradient ( F ) using the midpoint integration method and return the jacobian w.r.t. L.
         *
         * mode 1:
         * \f$F_{iI}^{t + 1} = \left[\delta_{ij} - \Delta t \left(1 - \alpha \right) L_{ij}^{t+1} \right]^{-1}
         * \left[F_{iI}^{t} + \Delta t \alpha \dot{F}_{iI}^{t} \right]\f$ \f$\frac{\partial F_{jI}^{t + 1}}{\partial
         * L_{kl}^{t+1}} = \left[\delta_{kj} - \Delta t \left(1 - \alpha\right) L_{kj}\right]^{-1} \Delta t \left(1 -
         * \alpha\right) F_{lI}^{t + 1}\f$
         *
         * mode 2:
         * \f$F_{iI}^{t + 1} = \left[F_{iJ}^{t} + \Delta t \alpha \dot{F}_{iJ}^{t} \right] \left[\delta_{IJ} - \Delta T
         * \left( 1- \alpha \right) L_{IJ}^{t+1} \right]^{-1}\f$ \f$\frac{\partial F_{iJ}^{t + 1}}{\partial L_{KL}} =
         * \Delta t (1 - \alpha) F_{iK}^{t + 1} \left[\delta_{JL} - \right ]\f$
         *
         * \param &Dt: The change in time.
         * \param &previousDeformationGradient: The previous value of the deformation gradient
         * \param &Lp: The previous velocity gradient.
         * \param &L: The current velocity gradient.
         * \param &dF: The change in the deformation gradient \f$\Delta \bf{F}\f$ such that \f$F_{iI}^{t+1} = F_{iI}^t +
         * \Delta F_{iI}\f$ \param &deformationGradient: The computed current deformation gradient. \param &dFdL: The
         * derivative of the deformation gradient w.r.t. the velocity gradient \param &ddFdFp: The derivative of the
         * change in the deformation gradient w.r.t. the previous deformation gradient \param &dFdFp: The derivative of
         * the deformation gradient w.r.t. the previous deformation gradient \param &dFdLp: The derivative of the
         * deformation gradient w.r.t. the previous velocity gradient \param alpha: The integration parameter ( 0 for
         * implicit, 1 for explicit ) \param mode: The form of the ODE. See above for details.
         */

        const unsigned int dim = (unsigned int)std::pow(previousDeformationGradient.size(), 0.5);

        dF                  = floatVector(dim * dim, 0);
        deformationGradient = floatVector(dim * dim, 0);
        dFdL                = floatVector(dim * dim * dim * dim, 0);
        ddFdFp              = floatVector(dim * dim * dim * dim, 0);
        dFdFp               = floatVector(dim * dim * dim * dim, 0);
        dFdLp               = floatVector(dim * dim * dim * dim, 0);

        TARDIGRADE_ERROR_TOOLS_CHECK((dim == 3) || (dim == 2) || (dim == 1),
                                     "The dimension of the deformation gradient is " + std::to_string(dim) +
                                         " but must be 1, 2, or 3");

        TARDIGRADE_ERROR_TOOLS_CHECK((mode == 1) || (mode == 2), "The mode of evolution " + std::to_string(mode) +
                                                                     " is not recognized. It must be 1 or 2.");

        if (mode == 1) {
            if (dim == 3) {
                evolveF<3, 1>(Dt, std::begin(previousDeformationGradient), std::end(previousDeformationGradient),
                              std::begin(Lp), std::end(Lp), std::begin(L), std::end(L), std::begin(dF), std::end(dF),
                              std::begin(deformationGradient), std::end(deformationGradient), std::begin(dFdL),
                              std::end(dFdL), std::begin(ddFdFp), std::end(ddFdFp), std::begin(dFdFp), std::end(dFdFp),
                              std::begin(dFdLp), std::end(dFdLp), alpha);
            } else if (dim == 2) {
                evolveF<2, 1>(Dt, std::begin(previousDeformationGradient), std::end(previousDeformationGradient),
                              std::begin(Lp), std::end(Lp), std::begin(L), std::end(L), std::begin(dF), std::end(dF),
                              std::begin(deformationGradient), std::end(deformationGradient), std::begin(dFdL),
                              std::end(dFdL), std::begin(ddFdFp), std::end(ddFdFp), std::begin(dFdFp), std::end(dFdFp),
                              std::begin(dFdLp), std::end(dFdLp), alpha);
            } else if (dim == 1) {
                evolveF<1, 1>(Dt, std::begin(previousDeformationGradient), std::end(previousDeformationGradient),
                              std::begin(Lp), std::end(Lp), std::begin(L), std::end(L), std::begin(dF), std::end(dF),
                              std::begin(deformationGradient), std::end(deformationGradient), std::begin(dFdL),
                              std::end(dFdL), std::begin(ddFdFp), std::end(ddFdFp), std::begin(dFdFp), std::end(dFdFp),
                              std::begin(dFdLp), std::end(dFdLp), alpha);
            }

        } else if (mode == 2) {
            if (dim == 3) {
                evolveF<3, 2>(Dt, std::begin(previousDeformationGradient), std::end(previousDeformationGradient),
                              std::begin(Lp), std::end(Lp), std::begin(L), std::end(L), std::begin(dF), std::end(dF),
                              std::begin(deformationGradient), std::end(deformationGradient), std::begin(dFdL),
                              std::end(dFdL), std::begin(ddFdFp), std::end(ddFdFp), std::begin(dFdFp), std::end(dFdFp),
                              std::begin(dFdLp), std::end(dFdLp), alpha);
            } else if (dim == 2) {
                evolveF<2, 2>(Dt, std::begin(previousDeformationGradient), std::end(previousDeformationGradient),
                              std::begin(Lp), std::end(Lp), std::begin(L), std::end(L), std::begin(dF), std::end(dF),
                              std::begin(deformationGradient), std::end(deformationGradient), std::begin(dFdL),
                              std::end(dFdL), std::begin(ddFdFp), std::end(ddFdFp), std::begin(dFdFp), std::end(dFdFp),
                              std::begin(dFdLp), std::end(dFdLp), alpha);
            } else if (dim == 1) {
                evolveF<1, 2>(Dt, std::begin(previousDeformationGradient), std::end(previousDeformationGradient),
                              std::begin(Lp), std::end(Lp), std::begin(L), std::end(L), std::begin(dF), std::end(dF),
                              std::begin(deformationGradient), std::end(deformationGradient), std::begin(dFdL),
                              std::end(dFdL), std::begin(ddFdFp), std::end(ddFdFp), std::begin(dFdFp), std::end(dFdFp),
                              std::begin(dFdLp), std::end(dFdLp), alpha);
            }
        }

        return;
    }

    void evolveFFlatJ(const floatType &Dt, const floatVector &previousDeformationGradient, const floatVector &Lp,
                      const floatVector &L, floatVector &deformationGradient, floatVector &dFdL, floatVector &dFdFp,
                      floatVector &dFdLp, const floatType alpha, const unsigned int mode) {
        /*!
         * Evolve the deformation gradient ( F ) using the midpoint integration method and return the jacobian w.r.t. L.
         *
         * mode 1:
         * \f$F_{iI}^{t + 1} = \left[\delta_{ij} - \Delta t \left(1 - \alpha \right) L_{ij}^{t+1} \right]^{-1}
         * \left[F_{iI}^{t} + \Delta t \alpha \dot{F}_{iI}^{t} \right]\f$ \f$\frac{\partial F_{jI}^{t + 1}}{\partial
         * L_{kl}^{t+1}} = \left[\delta_{kj} - \Delta t \left(1 - \alpha\right) L_{kj}\right]^{-1} \Delta t \left(1 -
         * \alpha\right) F_{lI}^{t + 1}\f$
         *
         * mode 2:
         * \f$F_{iI}^{t + 1} = \left[F_{iJ}^{t} + \Delta t \alpha \dot{F}_{iJ}^{t} \right] \left[\delta_{IJ} - \Delta T
         * \left( 1- \alpha \right) L_{IJ}^{t+1} \right]^{-1}\f$ \f$\frac{\partial F_{iJ}^{t + 1}}{\partial L_{KL}} =
         * \Delta t (1 - \alpha) F_{iK}^{t + 1} \left[\delta_{JL} - \right ]\f$
         *
         * \param &Dt: The change in time.
         * \param &previousDeformationGradient: The previous value of the deformation gradient
         * \param &Lp: The previous velocity gradient.
         * \param &L: The current velocity gradient.
         * \param &deformationGradient: The computed current deformation gradient.
         * \param &dFdL: The derivative of the deformation gradient w.r.t. the velocity gradient
         * \param &dFdFp: The derivative of the deformation gradient w.r.t. the previous deformation gradient
         * \param &dFdLp: The derivative of the deformation gradient w.r.t. the previous velocity gradient
         * \param alpha: The integration parameter ( 0 for implicit, 1 for explicit )
         * \param mode: The form of the ODE. See above for details.
         */

        floatVector dF;
        floatVector ddFdFp;

        return evolveFFlatJ(Dt, previousDeformationGradient, Lp, L, dF, deformationGradient, dFdL, ddFdFp, dFdFp, dFdLp,
                            alpha, mode);
    }

    void evolveF(const floatType &Dt, const floatVector &previousDeformationGradient, const floatVector &Lp,
                 const floatVector &L, floatVector &deformationGradient, floatMatrix &dFdL, floatMatrix &dFdFp,
                 floatMatrix &dFdLp, const floatType alpha, const unsigned int mode) {
        /*!
         * Evolve the deformation gradient ( F ) using the midpoint integration method and return the jacobian w.r.t. L.
         *
         * mode 1:
         * \f$F_{iI}^{t + 1} = \left[\delta_{ij} - \Delta t \left(1 - \alpha \right) L_{ij}^{t+1} \right]^{-1}
         * \left[F_{iI}^{t} + \Delta t \alpha \dot{F}_{iI}^{t} \right]\f$ \f$\frac{\partial F_{jI}^{t + 1}}{\partial
         * L_{kl}^{t+1}} = \left[\delta_{kj} - \Delta t \left(1 - \alpha\right) L_{kj}\right]^{-1} \Delta t \left(1 -
         * \alpha\right) F_{lI}^{t + 1}\f$
         *
         * mode 2:
         * \f$F_{iI}^{t + 1} = \left[F_{iJ}^{t} + \Delta t \alpha \dot{F}_{iJ}^{t} \right] \left[\delta_{IJ} - \Delta T
         * \left( 1- \alpha \right) L_{IJ}^{t+1} \right]^{-1}\f$ \f$\frac{\partial F_{iJ}^{t + 1}}{\partial L_{KL}} =
         * \Delta t (1 - \alpha) F_{iK}^{t + 1} \left[\delta_{JL} - \right ]\f$
         *
         * \param &Dt: The change in time.
         * \param &previousDeformationGradient: The previous value of the deformation gradient
         * \param &Lp: The previous velocity gradient.
         * \param &L: The current velocity gradient.
         * \param &deformationGradient: The computed current deformation gradient.
         * \param &dFdL: The derivative of the deformation gradient w.r.t. the velocity gradient
         * \param &dFdFp: The derivative of the deformation gradient w.r.t. the previous deformation gradient
         * \param &dFdLp: The derivative of the deformation gradient w.r.t. the previous velocity gradient
         * \param alpha: The integration parameter ( 0 for implicit, 1 for explicit )
         * \param mode: The form of the ODE. See above for details.
         */

        floatVector dF;
        floatMatrix ddFdFp;

        return evolveF(Dt, previousDeformationGradient, Lp, L, dF, deformationGradient, dFdL, ddFdFp, dFdFp, dFdLp,
                       alpha, mode);
    }

    void computeUnitNormal(const floatVector &A, floatVector &Anorm) {
        /*!
         * Compute the unit normal of a second order tensor (or strictly speaking
         * any tensor).
         *
         * \param &A: The second order tensor
         * \param &Anorm: The unit normal in the direction of A
         */

        Anorm = floatVector(A.size(), 0);

        computeUnitNormal(std::begin(A), std::end(A), std::begin(Anorm), std::end(Anorm));

        return;
    }

    void computeUnitNormal(const floatVector &A, floatVector &Anorm, floatVector &dAnormdA) {
        /*!
         * Compute the unit normal of a second order tensor (or strictly speaking any
         * tensor) and the gradient of that unit normal w.r.t. the tensor.
         *
         * \param &A: The second order tensor
         * \param &Anorm: The unit normal in the direction of A
         * \param &dAnormdA: The gradient of the unit normal w.r.t. A
         */

        Anorm    = floatVector(A.size(), 0);
        dAnormdA = floatVector(A.size() * A.size(), 0);

        computeUnitNormal(std::begin(A), std::end(A), std::begin(Anorm), std::end(Anorm), std::begin(dAnormdA),
                          std::end(dAnormdA));

        return;
    }

    void computeUnitNormal(const floatVector &A, floatVector &Anorm, floatMatrix &dAnormdA) {
        /*!
         * Compute the unit normal of a second order tensor (or strictly speaking any
         * tensor) and the gradient of that unit normal w.r.t. the tensor.
         *
         * \param &A: The second order tensor
         * \param &Anorm: The unit normal in the direction of A
         * \param &dAnormdA: The gradient of the unit normal w.r.t. A
         */

        const unsigned int A_size = A.size();

        floatVector _dAnormdA;

        TARDIGRADE_ERROR_TOOLS_CATCH(computeUnitNormal(A, Anorm, _dAnormdA));

        dAnormdA = tardigradeVectorTools::inflate(_dAnormdA, A_size, A_size);

        return;
    }

    void pullBackVelocityGradient(const floatVector &velocityGradient, const floatVector &deformationGradient,
                                  floatVector &pulledBackVelocityGradient) {
        /*!
         * Pull back the velocity gradient to the configuration indicated by deformationGradient, i.e.
         *
         * \f$totalDeformationGradient_{iI} = deformationGradient_{i \bar{I}} remainingDeformationGradient_{\bar{I}I}\f$
         *
         * This is done via
         *
         * \f$L_{\bar{I} \bar{J}} = deformationGradient_{\bar{I} i}^{-1} velocityGradient_{ij}
         * deformationGradient_{j\bar{J}}\f$
         *
         * \param &velocityGradient: The velocity gradient in the current configuration.
         * \param &deformationGradient: The deformation gradient between the desired configuration
         *     and the current configuration.
         * \param &pulledBackVelocityGradient: The pulled back velocity gradient.
         */

        const unsigned int dim = (unsigned int)std::pow(velocityGradient.size(), 0.5);

        pulledBackVelocityGradient = floatVector(dim * dim, 0);

        TARDIGRADE_ERROR_TOOLS_CHECK((dim == 3) || (dim == 2) || (dim == 1),
                                     "The dimension of the deformation gradient is " + std::to_string(dim) +
                                         " but must be 1, 2, or 3");

        if (dim == 3) {
            pullBackVelocityGradient<3>(std::begin(velocityGradient), std::end(velocityGradient),
                                        std::begin(deformationGradient), std::end(deformationGradient),
                                        std::begin(pulledBackVelocityGradient), std::end(pulledBackVelocityGradient));
        } else if (dim == 2) {
            pullBackVelocityGradient<2>(std::begin(velocityGradient), std::end(velocityGradient),
                                        std::begin(deformationGradient), std::end(deformationGradient),
                                        std::begin(pulledBackVelocityGradient), std::end(pulledBackVelocityGradient));
        } else if (dim == 1) {
            pullBackVelocityGradient<1>(std::begin(velocityGradient), std::end(velocityGradient),
                                        std::begin(deformationGradient), std::end(deformationGradient),
                                        std::begin(pulledBackVelocityGradient), std::end(pulledBackVelocityGradient));
        }

        return;
    }

    void pullBackVelocityGradient(const floatVector &velocityGradient, const floatVector &deformationGradient,
                                  floatVector &pulledBackVelocityGradient, floatVector &dPullBackLdL,
                                  floatVector &dPullBackLdF) {
        /*!
         * Pull back the velocity gradient to the configuration indicated by deformationGradient, i.e.
         *
         * \f$totalDeformationGradient_{iI} = deformationGradient_{i \bar{I}} remainingDeformationGradient_{\bar{I}I}\f$
         *
         * This is done via
         *
         * \f$L_{\bar{I} \bar{J}} = deformationGradient_{\bar{I} i}^{-1} velocityGradient_{ij}
         * deformationGradient_{j\bar{J}}\f$
         *
         * \param &velocityGradient: The velocity gradient in the current configuration.
         * \param &deformationGradient: The deformation gradient between the desired configuration
         *     and the current configuration.
         * \param &pulledBackVelocityGradient: The pulled back velocity gradient.
         * \param &dPullBackLdL: The gradient of the pulled back velocity gradient
         *     w.r.t. the velocity gradient.
         * \param &dPullBackLdF: The gradient of the pulled back velocity gradient
         *     w.r.t. the deformation gradient.
         */

        const unsigned int dim = (unsigned int)std::pow(velocityGradient.size(), 0.5);

        pulledBackVelocityGradient = floatVector(dim * dim, 0);
        dPullBackLdL               = floatVector(dim * dim * dim * dim, 0);
        dPullBackLdF               = floatVector(dim * dim * dim * dim, 0);

        TARDIGRADE_ERROR_TOOLS_CHECK((dim == 3) || (dim == 2) || (dim == 1),
                                     "The dimension of the deformation gradient is " + std::to_string(dim) +
                                         " but must be 1, 2, or 3");

        if (dim == 3) {
            pullBackVelocityGradient<3>(std::begin(velocityGradient), std::end(velocityGradient),
                                        std::begin(deformationGradient), std::end(deformationGradient),
                                        std::begin(pulledBackVelocityGradient), std::end(pulledBackVelocityGradient),
                                        std::begin(dPullBackLdL), std::end(dPullBackLdL), std::begin(dPullBackLdF),
                                        std::end(dPullBackLdF));
        } else if (dim == 2) {
            pullBackVelocityGradient<2>(std::begin(velocityGradient), std::end(velocityGradient),
                                        std::begin(deformationGradient), std::end(deformationGradient),
                                        std::begin(pulledBackVelocityGradient), std::end(pulledBackVelocityGradient),
                                        std::begin(dPullBackLdL), std::end(dPullBackLdL), std::begin(dPullBackLdF),
                                        std::end(dPullBackLdF));
        } else if (dim == 1) {
            pullBackVelocityGradient<1>(std::begin(velocityGradient), std::end(velocityGradient),
                                        std::begin(deformationGradient), std::end(deformationGradient),
                                        std::begin(pulledBackVelocityGradient), std::end(pulledBackVelocityGradient),
                                        std::begin(dPullBackLdL), std::end(dPullBackLdL), std::begin(dPullBackLdF),
                                        std::end(dPullBackLdF));
        }

        return;
    }

    void pullBackVelocityGradient(const floatVector &velocityGradient, const floatVector &deformationGradient,
                                  floatVector &pulledBackVelocityGradient, floatMatrix &dPullBackLdL,
                                  floatMatrix &dPullBackLdF) {
        /*!
         * Pull back the velocity gradient to the configuration indicated by deformationGradient, i.e.
         *
         * \f$totalDeformationGradient_{iI} = deformationGradient_{i \bar{I}} remainingDeformationGradient_{\bar{I}I}\f$
         *
         * This is done via
         *
         * \f$L_{\bar{I} \bar{J}} = deformationGradient_{\bar{I} i}^{-1} velocityGradient_{ij}
         * deformationGradient_{j\bar{J}}\f$
         *
         * \param &velocityGradient: The velocity gradient in the current configuration.
         * \param &deformationGradient: The deformation gradient between the desired configuration
         *     and the current configuration.
         * \param &pulledBackVelocityGradient: The pulled back velocity gradient.
         * \param &dPullBackLdL: The gradient of the pulled back velocity gradient
         *     w.r.t. the velocity gradient.
         * \param &dPullBackLdF: The gradient of the pulled back velocity gradient
         *     w.r.t. the deformation gradient.
         */

        // Assume 3D
        constexpr unsigned int dim     = 3;
        constexpr unsigned int sot_dim = dim * dim;

        floatVector _dPullBackLdL, _dPullBackLdF;

        TARDIGRADE_ERROR_TOOLS_CATCH(pullBackVelocityGradient(velocityGradient, deformationGradient,
                                                              pulledBackVelocityGradient, _dPullBackLdL,
                                                              _dPullBackLdF));

        dPullBackLdL = tardigradeVectorTools::inflate(_dPullBackLdL, sot_dim, sot_dim);

        dPullBackLdF = tardigradeVectorTools::inflate(_dPullBackLdF, sot_dim, sot_dim);

        return;
    }

    void quadraticThermalExpansion(const floatType &temperature, const floatType &referenceTemperature,
                                   const floatVector &linearParameters, const floatVector &quadraticParameters,
                                   floatVector &thermalExpansion) {
        /*!
         * Define a quadratic equation for the thermal expansion. This could be the
         * thermal strain or the value of the stretch tensor.
         *
         * \f$ e^{\theta}_{ij} = a_{ij} \left(\theta - \theta_0\right) + b_{ij} \left(\theta^2 - \theta_0^2\right)\f$
         *
         * Where \f$e^{\theta}_{ij}\f$ is the thermal expansion, \f$a_{ij}\f$ are the linear parameters,
         * \f$b_{ij}\f$ are the quadratic parameters, \f$\theta\f$ is the current temperature, and \f$\theta_0\f$
         * is the reference temperature.
         *
         * \param &temperature: The temperature
         * \param &referenceTemperature: The reference temperature
         * \param &linearParameters: The linear thermal expansion parameters.
         * \param &quadraticParameters: The quadratic thermal expansion parameters.
         * \param &thermalExpansion: The resulting thermal expansion.
         */

        thermalExpansion = floatVector(linearParameters.size(), 0);

        quadraticThermalExpansion(temperature, referenceTemperature, std::begin(linearParameters),
                                  std::end(linearParameters), std::begin(quadraticParameters),
                                  std::end(quadraticParameters), std::begin(thermalExpansion),
                                  std::end(thermalExpansion));

        return;
    }

    void quadraticThermalExpansion(const floatType &temperature, const floatType &referenceTemperature,
                                   const floatVector &linearParameters, const floatVector &quadraticParameters,
                                   floatVector &thermalExpansion, floatVector &thermalExpansionJacobian) {
        /*!
         * Define a quadratic equation for the thermal expansion. This could be the
         * thermal strain or the value of the stretch tensor.
         *
         * \f$ e^{\theta}_{ij} = a_{ij} \left(\theta - \theta_0\right) + b_{ij} \left(\theta^2 - \theta_0^2\right)\f$
         *
         * Where \f$e^{\theta}_{ij}\f$ is the thermal expansion, \f$a_{ij}\f$ are the linear parameters,
         * \f$b_{ij}\f$ are the quadratic parameters, \f$\theta\f$ is the current temperature, and \f$\theta_0\f$
         * is the reference temperature.
         *
         * \param &temperature: The temperature
         * \param &referenceTemperature: The reference temperature
         * \param &linearParameters: The linear thermal expansion parameters.
         * \param &quadraticParameters: The quadratic thermal expansion parameters.
         * \param &thermalExpansion: The resulting thermal expansion.
         * \param &thermalExpansionJacobian: The gradient of the thermal expansion w.r.t.
         *     the temperature.
         */

        thermalExpansion         = floatVector(linearParameters.size(), 0);
        thermalExpansionJacobian = floatVector(linearParameters.size(), 0);

        quadraticThermalExpansion(temperature, referenceTemperature, std::begin(linearParameters),
                                  std::end(linearParameters), std::begin(quadraticParameters),
                                  std::end(quadraticParameters), std::begin(thermalExpansion),
                                  std::end(thermalExpansion), std::begin(thermalExpansionJacobian),
                                  std::end(thermalExpansionJacobian));

        return;
    }

    void pushForwardGreenLagrangeStrain(const floatVector &greenLagrangeStrain, const floatVector &deformationGradient,
                                        floatVector &almansiStrain) {
        /*!
         * Push forward the Green-Lagrange strain to the current configuration.
         *
         * \f$e_{ij} = F_{Ii}^{-1} E_{IJ} F_{Jj}^{-1}\f$
         *
         * where \f$e_{ij}\f$ is the Almansi strain (the strain in the current configuration, \f$F_{iI}^{-1}\f$ is the
         * inverse of the deformation gradient, and \f$E_{IJ}\f$ is the Green-Lagrange strain.
         *
         * \param &greenLagrangeStrain: The Green-Lagrange strain.
         * \param &deformationGradient: The deformation gradient mapping between configurations.
         * \param &almansiStrain: The strain in the current configuration indicated by the deformation gradient.
         */

        const unsigned int dim = (unsigned int)std::pow(greenLagrangeStrain.size(), 0.5);

        almansiStrain = floatVector(dim * dim, 0);

        TARDIGRADE_ERROR_TOOLS_CHECK((dim == 3) || (dim == 2) || (dim == 1),
                                     "The dimension of the deformation gradient is " + std::to_string(dim) +
                                         " but must be 1, 2, or 3");

        if (dim == 3) {
            pushForwardGreenLagrangeStrain<3>(std::begin(greenLagrangeStrain), std::end(greenLagrangeStrain),
                                              std::begin(deformationGradient), std::end(deformationGradient),
                                              std::begin(almansiStrain), std::end(almansiStrain));
        } else if (dim == 2) {
            pushForwardGreenLagrangeStrain<2>(std::begin(greenLagrangeStrain), std::end(greenLagrangeStrain),
                                              std::begin(deformationGradient), std::end(deformationGradient),
                                              std::begin(almansiStrain), std::end(almansiStrain));
        } else if (dim == 1) {
            pushForwardGreenLagrangeStrain<1>(std::begin(greenLagrangeStrain), std::end(greenLagrangeStrain),
                                              std::begin(deformationGradient), std::end(deformationGradient),
                                              std::begin(almansiStrain), std::end(almansiStrain));
        }

        return;
    }

    void pushForwardGreenLagrangeStrain(const floatVector &greenLagrangeStrain, const floatVector &deformationGradient,
                                        floatVector &almansiStrain, floatVector &dAlmansiStraindE,
                                        floatVector &dAlmansiStraindF) {
        /*!
         * Push forward the Green-Lagrange strain to the current configuration
         * and return the jacobians.
         *
         * \f$e_{ij} = F_{Ii}^{-1} E_{IJ} F_{Jj}^{-1}\f$
         *
         * \f$\frac{\partial e_{ij}}{\partial E_{KL}} = F_{Ki}^{-1} F_{Kj}^{-1}\f$
         *
         * \f$\frac{\partial e_{ij}}{\partial F_{kK}} = -F_{Ik}^{-1} F_{Ki}^{-1} E_{IJ} F_{J j}^{-1} - F_{Ii}^{-1}
         * E_{IJ} F_{Jk}^{-1} F_{Kj}^{-1}\f$
         *
         * where \f$e_{ij}\f$ is the Almansi strain (the strain in the current configuration, \f$F_{iI}^{-1}\f$ is the
         * inverse of the deformation gradient, and \f$E_{IJ}\f$ is the Green-Lagrange strain.
         *
         * \param &greenLagrangeStrain: The Green-Lagrange strain.
         * \param &deformationGradient: The deformation gradient mapping between configurations.
         * \param &almansiStrain: The strain in the current configuration indicated by the deformation gradient.
         * \param &dAlmansiStraindE: Compute the derivative of the Almansi strain w.r.t. the Green-Lagrange strain.
         * \param &dAlmansiStraindF: Compute the derivative of the Almansi strain w.r.t. the deformation gradient.
         */

        const unsigned int dim = (unsigned int)std::pow(greenLagrangeStrain.size(), 0.5);

        almansiStrain    = floatVector(dim * dim, 0);
        dAlmansiStraindE = floatVector(dim * dim * dim * dim, 0);
        dAlmansiStraindF = floatVector(dim * dim * dim * dim, 0);

        TARDIGRADE_ERROR_TOOLS_CHECK((dim == 3) || (dim == 2) || (dim == 1),
                                     "The dimension of the deformation gradient is " + std::to_string(dim) +
                                         " but must be 1, 2, or 3");

        if (dim == 3) {
            pushForwardGreenLagrangeStrain<3>(std::begin(greenLagrangeStrain), std::end(greenLagrangeStrain),
                                              std::begin(deformationGradient), std::end(deformationGradient),
                                              std::begin(almansiStrain), std::end(almansiStrain),
                                              std::begin(dAlmansiStraindE), std::end(dAlmansiStraindE),
                                              std::begin(dAlmansiStraindF), std::end(dAlmansiStraindF));
        } else if (dim == 2) {
            pushForwardGreenLagrangeStrain<2>(std::begin(greenLagrangeStrain), std::end(greenLagrangeStrain),
                                              std::begin(deformationGradient), std::end(deformationGradient),
                                              std::begin(almansiStrain), std::end(almansiStrain),
                                              std::begin(dAlmansiStraindE), std::end(dAlmansiStraindE),
                                              std::begin(dAlmansiStraindF), std::end(dAlmansiStraindF));
        } else if (dim == 1) {
            pushForwardGreenLagrangeStrain<1>(std::begin(greenLagrangeStrain), std::end(greenLagrangeStrain),
                                              std::begin(deformationGradient), std::end(deformationGradient),
                                              std::begin(almansiStrain), std::end(almansiStrain),
                                              std::begin(dAlmansiStraindE), std::end(dAlmansiStraindE),
                                              std::begin(dAlmansiStraindF), std::end(dAlmansiStraindF));
        }

        return;
    }

    void pushForwardGreenLagrangeStrain(const floatVector &greenLagrangeStrain, const floatVector &deformationGradient,
                                        floatVector &almansiStrain, floatMatrix &dAlmansiStraindE,
                                        floatMatrix &dAlmansiStraindF) {
        /*!
         * Push forward the Green-Lagrange strain to the current configuration
         * and return the jacobians.
         *
         * \f$e_{ij} = F_{Ii}^{-1} E_{IJ} F_{Jj}^{-1}\f$
         *
         * \f$\frac{\partial e_{ij}}{\partial E_{KL}} = F_{Ki}^{-1} F_{Kj}^{-1}\f$
         *
         * \f$\frac{\partial e_{ij}}{\partial F_{kK}} = -F_{Ik}^{-1} F_{Ki}^{-1} E_{IJ} F_{J j}^{-1} - F_{Ii}^{-1}
         * E_{IJ} F_{Jk}^{-1} F_{Kj}^{-1}\f$
         *
         * where \f$e_{ij}\f$ is the Almansi strain (the strain in the current configuration, \f$F_{iI}^{-1}\f$ is the
         * inverse of the deformation gradient, and \f$E_{IJ}\f$ is the Green-Lagrange strain.
         *
         * \param &greenLagrangeStrain: The Green-Lagrange strain.
         * \param &deformationGradient: The deformation gradient mapping between configurations.
         * \param &almansiStrain: The strain in the current configuration indicated by the deformation gradient.
         * \param &dAlmansiStraindE: Compute the derivative of the Almansi strain w.r.t. the Green-Lagrange strain.
         * \param &dAlmansiStraindF: Compute the derivative of the Almansi strain w.r.t. the deformation gradient.
         */

        // Assume 3D
        constexpr unsigned int dim     = 3;
        constexpr unsigned int sot_dim = dim * dim;

        floatVector _dAlmansiStraindE, _dAlmansiStraindF;

        TARDIGRADE_ERROR_TOOLS_CATCH(pushForwardGreenLagrangeStrain(greenLagrangeStrain, deformationGradient,
                                                                    almansiStrain, _dAlmansiStraindE,
                                                                    _dAlmansiStraindF));

        dAlmansiStraindE = tardigradeVectorTools::inflate(_dAlmansiStraindE, sot_dim, sot_dim);
        dAlmansiStraindF = tardigradeVectorTools::inflate(_dAlmansiStraindF, sot_dim, sot_dim);

        return;
    }

    void pullBackAlmansiStrain(const floatVector &almansiStrain, const floatVector &deformationGradient,
                               floatVector &greenLagrangeStrain) {
        /*!
         * Pull back the almansi strain to the configuration indicated by the deformation gradient.
         *
         * \param &almansiStrain: The strain in the deformation gradient's current configuration.
         * \param &deformationGradient: The deformation gradient between configurations.
         * \param &greenLagrangeStrain: The Green-Lagrange strain which corresponds to the reference
         *     configuration of the deformation gradient.
         */

        const unsigned int dim = (unsigned int)std::pow(almansiStrain.size(), 0.5);

        greenLagrangeStrain = floatVector(dim * dim, 0);

        TARDIGRADE_ERROR_TOOLS_CHECK((dim == 3) || (dim == 2) || (dim == 1),
                                     "The dimension of the deformation gradient is " + std::to_string(dim) +
                                         " but must be 1, 2, or 3");

        if (dim == 3) {
            pullBackAlmansiStrain<3>(std::begin(almansiStrain), std::end(almansiStrain),
                                     std::begin(deformationGradient), std::end(deformationGradient),
                                     std::begin(greenLagrangeStrain), std::end(greenLagrangeStrain));
        } else if (dim == 2) {
            pullBackAlmansiStrain<2>(std::begin(almansiStrain), std::end(almansiStrain),
                                     std::begin(deformationGradient), std::end(deformationGradient),
                                     std::begin(greenLagrangeStrain), std::end(greenLagrangeStrain));
        } else if (dim == 1) {
            pullBackAlmansiStrain<1>(std::begin(almansiStrain), std::end(almansiStrain),
                                     std::begin(deformationGradient), std::end(deformationGradient),
                                     std::begin(greenLagrangeStrain), std::end(greenLagrangeStrain));
        }

        return;
    }

    void pullBackAlmansiStrain(const floatVector &almansiStrain, const floatVector &deformationGradient,
                               floatVector &greenLagrangeStrain, floatVector &dEde, floatVector &dEdF) {
        /*!
         * Pull back the almansi strain to the configuration indicated by the deformation gradient.
         *
         * Also return the Jacobians.
         *
         * \param &almansiStrain: The strain in the deformation gradient's current configuration.
         * \param &deformationGradient: The deformation gradient between configurations.
         * \param &greenLagrangeStrain: The Green-Lagrange strain which corresponds to the reference
         *     configuration of the deformation gradient.
         * \param &dEde: The derivative of the Green-Lagrange strain w.r.t. the Almansi strain.
         * \param &dEdF: The derivative of the Green-Lagrange strain w.r.t. the deformation gradient
         */

        const unsigned int dim = (unsigned int)std::pow(almansiStrain.size(), 0.5);

        greenLagrangeStrain = floatVector(dim * dim, 0);
        dEde                = floatVector(dim * dim * dim * dim, 0);
        dEdF                = floatVector(dim * dim * dim * dim, 0);

        TARDIGRADE_ERROR_TOOLS_CHECK((dim == 3) || (dim == 2) || (dim == 1),
                                     "The dimension of the deformation gradient is " + std::to_string(dim) +
                                         " but must be 1, 2, or 3");

        if (dim == 3) {
            pullBackAlmansiStrain<3>(std::begin(almansiStrain), std::end(almansiStrain),
                                     std::begin(deformationGradient), std::end(deformationGradient),
                                     std::begin(greenLagrangeStrain), std::end(greenLagrangeStrain), std::begin(dEde),
                                     std::end(dEde), std::begin(dEdF), std::end(dEdF));
        } else if (dim == 2) {
            pullBackAlmansiStrain<2>(std::begin(almansiStrain), std::end(almansiStrain),
                                     std::begin(deformationGradient), std::end(deformationGradient),
                                     std::begin(greenLagrangeStrain), std::end(greenLagrangeStrain), std::begin(dEde),
                                     std::end(dEde), std::begin(dEdF), std::end(dEdF));
        } else if (dim == 1) {
            pullBackAlmansiStrain<1>(std::begin(almansiStrain), std::end(almansiStrain),
                                     std::begin(deformationGradient), std::end(deformationGradient),
                                     std::begin(greenLagrangeStrain), std::end(greenLagrangeStrain), std::begin(dEde),
                                     std::end(dEde), std::begin(dEdF), std::end(dEdF));
        }

        return;
    }

    void pullBackAlmansiStrain(const floatVector &almansiStrain, const floatVector &deformationGradient,
                               floatVector &greenLagrangeStrain, floatMatrix &dEde, floatMatrix &dEdF) {
        /*!
         * Pull back the almansi strain to the configuration indicated by the deformation gradient.
         *
         * Also return the Jacobians.
         *
         * \param &almansiStrain: The strain in the deformation gradient's current configuration.
         * \param &deformationGradient: The deformation gradient between configurations.
         * \param &greenLagrangeStrain: The Green-Lagrange strain which corresponds to the reference
         *     configuration of the deformation gradient.
         * \param &dEde: The derivative of the Green-Lagrange strain w.r.t. the Almansi strain.
         * \param &dEdF: The derivative of the Green-Lagrange strain w.r.t. the deformation gradient
         */

        // Assume 3d
        constexpr unsigned int dim     = 3;
        constexpr unsigned int sot_dim = dim * dim;

        floatVector _dEde, _dEdF;

        TARDIGRADE_ERROR_TOOLS_CATCH(pullBackAlmansiStrain(almansiStrain, deformationGradient, greenLagrangeStrain,
                                                           _dEde, _dEdF))

        dEde = tardigradeVectorTools::inflate(_dEde, sot_dim, sot_dim);
        dEdF = tardigradeVectorTools::inflate(_dEdF, sot_dim, sot_dim);

        return;
    }

    void computeSymmetricPart(const floatVector &A, floatVector &symmA, unsigned int &dim) {
        /*!
         * Compute the symmetric part of a second order tensor ( \f$A\f$ ) and return it.
         *
         * \f$symm( A )_ij = \frac{1}{2}\left(A_{ij} + A_{ji}\right)\f$
         *
         * \param &A: A constant reference to the second order tensor to process ( \f$A\f$ )
         * \param &symmA: The symmetric part of A ( \f$A^{symm}\f$ )
         * \param &dim: The dimension of A. Note that this is an output used for help
         *     with computing the Jacobian. If you don't need dim as an output use the
         *     version of this function without it.
         */

        symmA = floatVector(A.size());

        computeSymmetricPart(std::begin(A), std::end(A), std::begin(symmA), std::end(symmA), dim);

        return;
    }

    void computeSymmetricPart(const floatVector &A, floatVector &symmA) {
        /*!
         * Compute the symmetric part of a second order tensor ( \f$A\f$ ) and return it.
         *
         * \f$symm( A )_ij = \frac{1}{2}\left(A_{ij} + A_{ji}\right)\f$
         *
         * \param &A: A constant reference to the second order tensor to process ( \f$A\f$ )
         * \param &symmA: The symmetric part of A ( \f$A^{symm}\f$ )
         */

        unsigned int dim;
        symmA = floatVector(A.size());
        return computeSymmetricPart(std::begin(A), std::end(A), std::begin(symmA), std::end(symmA), dim);
    }

    void computeSymmetricPart(const floatVector &A, floatVector &symmA, floatVector &dSymmAdA) {
        /*!
         * Compute the symmetric part of a second order tensor ( \f$A\f$ ) and return it.
         *
         * \f$( A )^{symm}_{ij} = \frac{1}{2}\left(A_{ij} + A_{ji}\right)\f$
         *
         * Also computes the jacobian
         *
         * \f$\frac{\partial A^{symm}_{ij}}{\partial A_{kl}} = \frac{1}{2}\left( \delta_{ik} \delta_{jl} +
         * \delta_{jk}\delta_{il} \right)\f$
         *
         * \param &A: A constant reference to the second order tensor to process ( \f$A\f$ )
         * \param &symmA: The symmetric part of A ( \f$A^{symm}\f$ )
         * \param &dSymmAdA: The Jacobian of the symmetric part of A w.r.t. A ( \f$\frac{\partial A^{symm}}{\partial
         * A}\f$ )
         */

        symmA    = floatVector(A.size(), 0);
        dSymmAdA = floatVector(A.size() * A.size(), 0);

        computeSymmetricPart(std::begin(A), std::end(A), std::begin(symmA), std::end(symmA), std::begin(dSymmAdA),
                             std::end(dSymmAdA));

        return;
    }

    void computeSymmetricPart(const floatVector &A, floatVector &symmA, floatMatrix &dSymmAdA) {
        /*!
         * Compute the symmetric part of a second order tensor ( \f$A\f$ ) and return it.
         *
         * \f$( A )^{symm}_{ij} = \frac{1}{2}\left(A_{ij} + A_{ji}\right)\f$
         *
         * Also computes the jacobian
         *
         * \f$\frac{\partial A^{symm}_{ij}}{\partial A_{kl}} = \frac{1}{2}\left( \delta_{ik} \delta_{jl} +
         * \delta_{jk}\delta_{il} \right) \f$
         *
         * \param &A: A constant reference to the second order tensor to process ( \f$A\f$ )
         * \param &symmA: The symmetric part of A ( \f$A^{symm}\f$ )
         * \param &dSymmAdA: The Jacobian of the symmetric part of A w.r.t. A ( \f$\frac{\partial A^{symm}}{\partial
         * A}\f$ )
         */

        symmA = floatVector(A.size());
        floatVector _dSymmAdA(A.size() * A.size());

        TARDIGRADE_ERROR_TOOLS_CATCH(computeSymmetricPart(std::begin(A), std::end(A), std::begin(symmA),
                                                          std::end(symmA), std::begin(_dSymmAdA), std::end(_dSymmAdA)));

        dSymmAdA = tardigradeVectorTools::inflate(_dSymmAdA, A.size(), A.size());

        return;
    }

    void pushForwardPK2Stress(const floatVector &PK2, const floatVector &F, floatVector &cauchyStress) {
        /*!
         * Push the Second Piola-Kirchhoff stress forward to the current configuration resulting in the Cauchy stress
         *
         * \f$ \sigma_{ij} = \frac{1}{J} F_{iI} S_{IJ} F_{jJ} \f$
         *
         * \param &PK2: The Second Piola-Kirchhoff stress \f$ S_{IJ} \f$
         * \param &F: The deformation gradient \f$ F_{iI} \f$
         * \param &cauchyStress: The Cauchy stress \f$ \sigma_{ij} \f$
         */

        const unsigned int dim = (unsigned int)std::pow(PK2.size(), 0.5);

        cauchyStress = floatVector(dim * dim, 0);

        TARDIGRADE_ERROR_TOOLS_CHECK((dim == 3) || (dim == 2) || (dim == 1), "The dimension of the PK2 stress is " +
                                                                                 std::to_string(dim) +
                                                                                 " but must be 1, 2, or 3");

        if (dim == 3) {
            pushForwardPK2Stress<3>(std::begin(PK2), std::end(PK2), std::begin(F), std::end(F),
                                    std::begin(cauchyStress), std::end(cauchyStress));
        } else if (dim == 2) {
            pushForwardPK2Stress<2>(std::begin(PK2), std::end(PK2), std::begin(F), std::end(F),
                                    std::begin(cauchyStress), std::end(cauchyStress));
        } else if (dim == 1) {
            pushForwardPK2Stress<1>(std::begin(PK2), std::end(PK2), std::begin(F), std::end(F),
                                    std::begin(cauchyStress), std::end(cauchyStress));
        }

        return;
    }

    void pushForwardPK2Stress(const floatVector &PK2, const floatVector &F, floatVector &cauchyStress,
                              floatVector &dCauchyStressdPK2, floatVector &dCauchyStressdF) {
        /*!
         * Push the Second Piola-Kirchhoff stress forward to the current configuration resulting in the Cauchy stress
         *
         * \f$ \sigma_{ij} = \frac{1}{J} F_{iI} S_{IJ} F_{jJ} \f$
         *
         * \param &PK2: The Second Piola-Kirchhoff stress \f$ S_{IJ} \f$
         * \param &F: The deformation gradient \f$ F_{iI} \f$
         * \param &cauchyStress: The Cauchy stress \f$ \sigma_{ij} \f$
         * \param &dCauchyStressdPK2: The gradient of the Cauchy stress w.r.t. the PK2 stress
         * \param &dCauchyStressdF: The gradient of the Cauchy stress w.r.t. the deformation gradient
         */

        const unsigned int dim = (unsigned int)std::pow(PK2.size(), 0.5);

        cauchyStress      = floatVector(dim * dim, 0);
        dCauchyStressdPK2 = floatVector(dim * dim * dim * dim, 0);
        dCauchyStressdF   = floatVector(dim * dim * dim * dim, 0);

        TARDIGRADE_ERROR_TOOLS_CHECK((dim == 3) || (dim == 2) || (dim == 1), "The dimension of the PK2 stress is " +
                                                                                 std::to_string(dim) +
                                                                                 " but must be 1, 2, or 3");

        if (dim == 3) {
            pushForwardPK2Stress<3>(std::begin(PK2), std::end(PK2), std::begin(F), std::end(F),
                                    std::begin(cauchyStress), std::end(cauchyStress), std::begin(dCauchyStressdPK2),
                                    std::end(dCauchyStressdPK2), std::begin(dCauchyStressdF),
                                    std::end(dCauchyStressdF));
        } else if (dim == 2) {
            pushForwardPK2Stress<2>(std::begin(PK2), std::end(PK2), std::begin(F), std::end(F),
                                    std::begin(cauchyStress), std::end(cauchyStress), std::begin(dCauchyStressdPK2),
                                    std::end(dCauchyStressdPK2), std::begin(dCauchyStressdF),
                                    std::end(dCauchyStressdF));
        } else if (dim == 1) {
            pushForwardPK2Stress<1>(std::begin(PK2), std::end(PK2), std::begin(F), std::end(F),
                                    std::begin(cauchyStress), std::end(cauchyStress), std::begin(dCauchyStressdPK2),
                                    std::end(dCauchyStressdPK2), std::begin(dCauchyStressdF),
                                    std::end(dCauchyStressdF));
        }

        return;
    }

    void pushForwardPK2Stress(const floatVector &PK2, const floatVector &F, floatVector &cauchyStress,
                              floatMatrix &dCauchyStressdPK2, floatMatrix &dCauchyStressdF) {
        /*!
         * Push the Second Piola-Kirchhoff stress forward to the current configuration resulting in the Cauchy stress
         *
         * \f$ \sigma_{ij} = \frac{1}{J} F_{iI} S_{IJ} F_{jJ} \f$
         *
         * \param &PK2: The Second Piola-Kirchhoff stress \f$ S_{IJ} \f$
         * \param &F: The deformation gradient \f$ F_{iI} \f$
         * \param &cauchyStress: The Cauchy stress \f$ \sigma_{ij} \f$
         * \param &dCauchyStressdPK2: The gradient of the Cauchy stress w.r.t. the PK2 stress
         * \param &dCauchyStressdF: The gradient of the Cauchy stress w.r.t. the deformation gradient
         */

        floatVector _dCauchyStressdPK2, _dCauchyStressdF;

        TARDIGRADE_ERROR_TOOLS_CATCH(pushForwardPK2Stress(PK2, F, cauchyStress, _dCauchyStressdPK2, _dCauchyStressdF))

        dCauchyStressdPK2 = tardigradeVectorTools::inflate(_dCauchyStressdPK2, cauchyStress.size(), PK2.size());

        dCauchyStressdF = tardigradeVectorTools::inflate(_dCauchyStressdF, cauchyStress.size(), F.size());

        return;
    }

    void pullBackCauchyStress(const floatVector &cauchyStress, const floatVector &F, floatVector &PK2) {
        /*!
         * Pull back the Cauchy stress to an earlier configuration resulting in the second Piola-Kirchhoff stress
         *
         * \f$ S_{IJ} = J F^{-1}_{Ii} \sigma_{ij} F^{-1}_{Jj} \f$
         *
         * where \f$S_{IJ}\f$ are the components of the second Piola-Kirchhoff stress tensor, \f$J \f$ is the
         * determinant of the deformation gradient \f$\bf{F}\f$ which has components \f$F_{iI}\f$, and
         * \f$ \sigma_{ij} \f$ are the components of the Cauchy stress.
         *
         * \param &cauchyStress: The cauchy stress tensor in row-major form (all nine components)
         * \param &F: The deformation gradient
         * \param &PK2: The resulting second Piola-Kirchhoff stress
         */

        const unsigned int dim = (unsigned int)std::pow(cauchyStress.size(), 0.5);

        PK2 = floatVector(dim * dim, 0);

        TARDIGRADE_ERROR_TOOLS_CHECK((dim == 3) || (dim == 2) || (dim == 1), "The dimension of the Cauchy stress is " +
                                                                                 std::to_string(dim) +
                                                                                 " but must be 1, 2, or 3");

        if (dim == 3) {
            pullBackCauchyStress<3>(std::begin(cauchyStress), std::end(cauchyStress), std::begin(F), std::end(F),
                                    std::begin(PK2), std::end(PK2));
        } else if (dim == 2) {
            pullBackCauchyStress<2>(std::begin(cauchyStress), std::end(cauchyStress), std::begin(F), std::end(F),
                                    std::begin(PK2), std::end(PK2));
        } else if (dim == 1) {
            pullBackCauchyStress<1>(std::begin(cauchyStress), std::end(cauchyStress), std::begin(F), std::end(F),
                                    std::begin(PK2), std::end(PK2));
        }

        return;
    }

    void pullBackCauchyStress(const floatVector &cauchyStress, const floatVector &F, floatVector &PK2,
                              floatVector &dPK2dCauchyStress, floatVector &dPK2dF) {
        /*!
         * Pull back the Cauchy stress to an earlier configuration resulting in the second Piola-Kirchhoff stress
         *
         * \f$ S_{IJ} = J F^{-1}_{Ii} \sigma_{ij} F^{-1}_{Jj} \f$
         *
         * where \f$S_{IJ}\f$ are the components of the second Piola-Kirchhoff stress tensor, \f$J \f$ is the
         * determinant of the deformation gradient \f$\bf{F}\f$ which has components \f$F_{iI}\f$, and
         * \f$ \sigma_{ij} \f$ are the components of the Cauchy stress.
         *
         * \param &cauchyStress: The cauchy stress tensor in row-major form (all nine components)
         * \param &F: The deformation gradient
         * \param &PK2: The resulting second Piola-Kirchhoff stress
         * \param &dPK2dCauchyStress: The directional derivative of the second Piola-Kirchhoff stress tensor w.r.t.
         *     the Cauchy stress
         * \param &dPK2dF: The directional derivative of the second Piola-Kirchhoff stress tensor w.r.t. the
         *     deformation gradient
         */

        const unsigned int dim = (unsigned int)std::pow(cauchyStress.size(), 0.5);

        PK2               = floatVector(dim * dim, 0);
        dPK2dCauchyStress = floatVector(dim * dim * dim * dim, 0);
        dPK2dF            = floatVector(dim * dim * dim * dim, 0);

        TARDIGRADE_ERROR_TOOLS_CHECK((dim == 3) || (dim == 2) || (dim == 1), "The dimension of the Cauchy stress is " +
                                                                                 std::to_string(dim) +
                                                                                 " but must be 1, 2, or 3");

        if (dim == 3) {
            pullBackCauchyStress<3>(std::begin(cauchyStress), std::end(cauchyStress), std::begin(F), std::end(F),
                                    std::begin(PK2), std::end(PK2), std::begin(dPK2dCauchyStress),
                                    std::end(dPK2dCauchyStress), std::begin(dPK2dF), std::end(dPK2dF));
        } else if (dim == 2) {
            pullBackCauchyStress<2>(std::begin(cauchyStress), std::end(cauchyStress), std::begin(F), std::end(F),
                                    std::begin(PK2), std::end(PK2), std::begin(dPK2dCauchyStress),
                                    std::end(dPK2dCauchyStress), std::begin(dPK2dF), std::end(dPK2dF));
        } else if (dim == 1) {
            pullBackCauchyStress<1>(std::begin(cauchyStress), std::end(cauchyStress), std::begin(F), std::end(F),
                                    std::begin(PK2), std::end(PK2), std::begin(dPK2dCauchyStress),
                                    std::end(dPK2dCauchyStress), std::begin(dPK2dF), std::end(dPK2dF));
        }

        return;
    }

    void pullBackCauchyStress(const floatVector &cauchyStress, const floatVector &F, floatVector &PK2,
                              floatMatrix &dPK2dCauchyStress, floatMatrix &dPK2dF) {
        /*!
         * Pull back the Cauchy stress to an earlier configuration resulting in the second Piola-Kirchhoff stress
         *
         * \f$ S_{IJ} = J F^{-1}_{Ii} \sigma_{ij} F^{-1}_{Jj} \f$
         *
         * where \f$S_{IJ}\f$ are the components of the second Piola-Kirchhoff stress tensor, \f$J \f$ is the
         * determinant of the deformation gradient \f$\bf{F}\f$ which has components \f$F_{iI}\f$, and
         * \f$ \sigma_{ij} \f$ are the components of the Cauchy stress.
         *
         * \param &cauchyStress: The cauchy stress tensor in row-major form (all nine components)
         * \param &F: The deformation gradient
         * \param &PK2: The resulting second Piola-Kirchhoff stress
         * \param &dPK2dCauchyStress: The directional derivative of the second Piola-Kirchhoff stress tensor w.r.t.
         *     the Cauchy stress
         * \param &dPK2dF: The directional derivative of the second Piola-Kirchhoff stress tensor w.r.t. the
         *     deformation gradient
         */

        floatVector _dPK2dCauchyStress, _dPK2dF;

        TARDIGRADE_ERROR_TOOLS_CATCH(pullBackCauchyStress(cauchyStress, F, PK2, _dPK2dCauchyStress, _dPK2dF))

        dPK2dCauchyStress = tardigradeVectorTools::inflate(_dPK2dCauchyStress, PK2.size(), cauchyStress.size());

        dPK2dF = tardigradeVectorTools::inflate(_dPK2dF, PK2.size(), F.size());

        return;
    }

    void evolveFExponentialMap(const floatType &Dt, const floatVector &previousDeformationGradient,
                               const floatVector &Lp, const floatVector &L, floatVector &deformationGradient,
                               const floatType alpha) {
        /*!
         * Evolve the deformation gradient using the exponential map. Assumes the evolution equation is of the form
         *
         * \f$ \dot{F}_{iI} = \ell_{ij} F_{jI} \f$
         *
         * \param &Dt: The change in time
         * \param &previousDeformationGradient: The previous value of the deformation gradient
         * \param &Lp: The previous value of the velocity gradient
         * \param &L: The current value of the velocity gradient
         * \param &deformationGradient: The computed value of the deformation gradient
         * \param &alpha: The integration parameter (0 is explicit and 1 is implicit)
         */

        const unsigned int dim = (unsigned int)std::pow(previousDeformationGradient.size(), 0.5);

        deformationGradient = floatVector(dim * dim, 0);

        TARDIGRADE_ERROR_TOOLS_CHECK((dim == 3) || (dim == 2) || (dim == 1),
                                     "The dimension of the deformation gradient is " + std::to_string(dim) +
                                         " but must be 1, 2, or 3");

        if (dim == 3) {
            evolveFExponentialMap<3>(Dt, std::begin(previousDeformationGradient), std::end(previousDeformationGradient),
                                     std::begin(Lp), std::end(Lp), std::begin(L), std::end(L),
                                     std::begin(deformationGradient), std::end(deformationGradient), alpha);
        } else if (dim == 2) {
            evolveFExponentialMap<2>(Dt, std::begin(previousDeformationGradient), std::end(previousDeformationGradient),
                                     std::begin(Lp), std::end(Lp), std::begin(L), std::end(L),
                                     std::begin(deformationGradient), std::end(deformationGradient), alpha);
        } else if (dim == 1) {
            evolveFExponentialMap<1>(Dt, std::begin(previousDeformationGradient), std::end(previousDeformationGradient),
                                     std::begin(Lp), std::end(Lp), std::begin(L), std::end(L),
                                     std::begin(deformationGradient), std::end(deformationGradient), alpha);
        }
    }

    void evolveFExponentialMap(const floatType &Dt, const floatVector &previousDeformationGradient,
                               const floatVector &Lp, const floatVector &L, floatVector &deformationGradient,
                               floatVector &dFdL, const floatType alpha) {
        /*!
         * Evolve the deformation gradient using the exponential map. Assumes the evolution equation is of the form
         *
         * \f$ \dot{F}_{iI} = \ell_{ij} F_{jI} \f$
         *
         * \param &Dt: The change in time
         * \param &previousDeformationGradient: The previous value of the deformation gradient
         * \param &Lp: The previous value of the velocity gradient
         * \param &L: The current value of the velocity gradient
         * \param &deformationGradient: The computed value of the deformation gradient
         * \param &dFdL: The derivative of the deformation gradient w.r.t. the velocity gradient
         * \param &alpha: The integration parameter (0 is explicit and 1 is implicit)
         */

        const unsigned int dim = (unsigned int)std::pow(previousDeformationGradient.size(), 0.5);

        deformationGradient = floatVector(dim * dim, 0);
        dFdL                = floatVector(dim * dim * dim * dim, 0);

        TARDIGRADE_ERROR_TOOLS_CHECK((dim == 3) || (dim == 2) || (dim == 1),
                                     "The dimension of the deformation gradient is " + std::to_string(dim) +
                                         " but must be 1, 2, or 3");

        if (dim == 3) {
            evolveFExponentialMap<3>(Dt, std::begin(previousDeformationGradient), std::end(previousDeformationGradient),
                                     std::begin(Lp), std::end(Lp), std::begin(L), std::end(L),
                                     std::begin(deformationGradient), std::end(deformationGradient), std::begin(dFdL),
                                     std::end(dFdL), alpha);
        } else if (dim == 2) {
            evolveFExponentialMap<2>(Dt, std::begin(previousDeformationGradient), std::end(previousDeformationGradient),
                                     std::begin(Lp), std::end(Lp), std::begin(L), std::end(L),
                                     std::begin(deformationGradient), std::end(deformationGradient), std::begin(dFdL),
                                     std::end(dFdL), alpha);
        } else if (dim == 1) {
            evolveFExponentialMap<1>(Dt, std::begin(previousDeformationGradient), std::end(previousDeformationGradient),
                                     std::begin(Lp), std::end(Lp), std::begin(L), std::end(L),
                                     std::begin(deformationGradient), std::end(deformationGradient), std::begin(dFdL),
                                     std::end(dFdL), alpha);
        }
    }

    void evolveFExponentialMap(const floatType &Dt, const floatVector &previousDeformationGradient,
                               const floatVector &Lp, const floatVector &L, floatVector &deformationGradient,
                               floatVector &dFdL, floatVector &dFdFp, floatVector &dFdLp, const floatType alpha) {
        /*!
         * Evolve the deformation gradient using the exponential map. Assumes the evolution equation is of the form
         *
         * \f$ \dot{F}_{iI} = \ell_{ij} F_{jI} \f$
         *
         * \param &Dt: The change in time
         * \param &previousDeformationGradient: The previous value of the deformation gradient
         * \param &Lp: The previous value of the velocity gradient
         * \param &L: The current value of the velocity gradient
         * \param &deformationGradient: The computed value of the deformation gradient
         * \param &dFdL: The derivative of the deformation gradient w.r.t. the velocity gradient
         * \param &dFdFp: The derivative of the deformation gradient w.r.t. the previous deformation gradient
         * \param &dFdLp: The derivative of the deformation gradient w.r.t. the previous velocity gradient
         * \param &alpha: The integration parameter (0 is explicit and 1 is implicit)
         */

        const unsigned int dim = (unsigned int)std::pow(previousDeformationGradient.size(), 0.5);

        deformationGradient = floatVector(dim * dim, 0);
        dFdL                = floatVector(dim * dim * dim * dim, 0);
        dFdFp               = floatVector(dim * dim * dim * dim, 0);
        dFdLp               = floatVector(dim * dim * dim * dim, 0);

        TARDIGRADE_ERROR_TOOLS_CHECK((dim == 3) || (dim == 2) || (dim == 1),
                                     "The dimension of the deformation gradient is " + std::to_string(dim) +
                                         " but must be 1, 2, or 3");

        if (dim == 3) {
            evolveFExponentialMap<3>(Dt, std::begin(previousDeformationGradient), std::end(previousDeformationGradient),
                                     std::begin(Lp), std::end(Lp), std::begin(L), std::end(L),
                                     std::begin(deformationGradient), std::end(deformationGradient), std::begin(dFdL),
                                     std::end(dFdL), std::begin(dFdFp), std::end(dFdFp), std::begin(dFdLp),
                                     std::end(dFdLp), alpha);
        } else if (dim == 2) {
            evolveFExponentialMap<2>(Dt, std::begin(previousDeformationGradient), std::end(previousDeformationGradient),
                                     std::begin(Lp), std::end(Lp), std::begin(L), std::end(L),
                                     std::begin(deformationGradient), std::end(deformationGradient), std::begin(dFdL),
                                     std::end(dFdL), std::begin(dFdFp), std::end(dFdFp), std::begin(dFdLp),
                                     std::end(dFdLp), alpha);
        } else if (dim == 1) {
            evolveFExponentialMap<1>(Dt, std::begin(previousDeformationGradient), std::end(previousDeformationGradient),
                                     std::begin(Lp), std::end(Lp), std::begin(L), std::end(L),
                                     std::begin(deformationGradient), std::end(deformationGradient), std::begin(dFdL),
                                     std::end(dFdL), std::begin(dFdFp), std::end(dFdFp), std::begin(dFdLp),
                                     std::end(dFdLp), alpha);
        }
    }

    void computeDCurrentNormalVectorDF(const floatVector &normalVector, const floatVector &F,
                                       floatVector &dNormalVectordF) {
        /*!
         * Compute the derivative of the normal vector in the current configuration w.r.t. the deformation gradient
         *
         * \param &normalVector: The unit normal vector in the current configuration
         * \param &F: The deformation gradient
         * \param &dNormalVectordF: The derivative of the normal vector w.r.t. the deformation gradient
         */

        const unsigned int dim = normalVector.size();

        dNormalVectordF = floatVector(dim * dim * dim, 0);

        TARDIGRADE_ERROR_TOOLS_CHECK((dim == 3) || (dim == 2) || (dim == 1),
                                     "The dimension of the deformation gradient is " + std::to_string(dim) +
                                         " but must be 1, 2, or 3");

        if (dim == 3) {
            computeDCurrentNormalVectorDF<3>(std::begin(normalVector), std::end(normalVector), std::begin(F),
                                             std::end(F), std::begin(dNormalVectordF), std::end(dNormalVectordF));
        } else if (dim == 2) {
            computeDCurrentNormalVectorDF<2>(std::begin(normalVector), std::end(normalVector), std::begin(F),
                                             std::end(F), std::begin(dNormalVectordF), std::end(dNormalVectordF));
        } else if (dim == 1) {
            computeDCurrentNormalVectorDF<1>(std::begin(normalVector), std::end(normalVector), std::begin(F),
                                             std::end(F), std::begin(dNormalVectordF), std::end(dNormalVectordF));
        }

        return;
    }

    void computeDCurrentAreaWeightedNormalVectorDF(const floatVector &normalVector, const floatVector &F,
                                                   floatVector &dAreaWeightedNormalVectordF) {
        /*!
         * Compute the derivative of the area weighted normal vector w.r.t. the deformation gradient i.e.
         *
         * \f$ \frac{\partial}{\partial F_{bB}} \left( n_i da \right) \f$
         *
         * Note that if the user passes in the unit normal vector, then the result will be more convenient for the
         * construction of the jacobian of a surface integral in the current configuration.
         *
         * \param &normalVector: The normal vector (a unit vector is likely what is desired)
         * \param &F: The deformation gradient
         * \param &dAreaWeightedNormalVectordF: The derivative of the area weighted normal vector w.r.t. the deformation
         * gradient
         */

        const unsigned int dim = normalVector.size();

        dAreaWeightedNormalVectordF = floatVector(dim * dim * dim, 0);

        TARDIGRADE_ERROR_TOOLS_CHECK((dim == 3) || (dim == 2) || (dim == 1),
                                     "The dimension of the deformation gradient is " + std::to_string(dim) +
                                         " but must be 1, 2, or 3");

        if (dim == 3) {
            computeDCurrentAreaWeightedNormalVectorDF<3>(std::begin(normalVector), std::end(normalVector),
                                                         std::begin(F), std::end(F),
                                                         std::begin(dAreaWeightedNormalVectordF),
                                                         std::end(dAreaWeightedNormalVectordF));
        } else if (dim == 2) {
            computeDCurrentAreaWeightedNormalVectorDF<2>(std::begin(normalVector), std::end(normalVector),
                                                         std::begin(F), std::end(F),
                                                         std::begin(dAreaWeightedNormalVectordF),
                                                         std::end(dAreaWeightedNormalVectordF));
        } else if (dim == 1) {
            computeDCurrentAreaWeightedNormalVectorDF<1>(std::begin(normalVector), std::end(normalVector),
                                                         std::begin(F), std::end(F),
                                                         std::begin(dAreaWeightedNormalVectordF),
                                                         std::end(dAreaWeightedNormalVectordF));
        }

        return;
    }

    void computeDCurrentAreaDF(const floatVector &normalVector, const floatVector &F, floatVector &dCurrentAreadF) {
        /*!
         * Compute the derivative of the current area w.r.t. the deformation gradient
         *
         * \param &normalVector: The current unit normal vector
         * \param &F: The deformation gradient
         * \param &dCurrentAreadF: The derivative of the current surface area w.r.t. F
         */

        const unsigned int dim = normalVector.size();

        dCurrentAreadF = floatVector(dim * dim, 0);

        TARDIGRADE_ERROR_TOOLS_CHECK((dim == 3) || (dim == 2) || (dim == 1),
                                     "The dimension of the deformation gradient is " + std::to_string(dim) +
                                         " but must be 1, 2, or 3");

        if (dim == 3) {
            computeDCurrentAreaDF<3>(std::begin(normalVector), std::end(normalVector), std::begin(F), std::end(F),
                                     std::begin(dCurrentAreadF), std::end(dCurrentAreadF));
        } else if (dim == 2) {
            computeDCurrentAreaDF<2>(std::begin(normalVector), std::end(normalVector), std::begin(F), std::end(F),
                                     std::begin(dCurrentAreadF), std::end(dCurrentAreadF));
        } else if (dim == 1) {
            computeDCurrentAreaDF<1>(std::begin(normalVector), std::end(normalVector), std::begin(F), std::end(F),
                                     std::begin(dCurrentAreadF), std::end(dCurrentAreadF));
        }

        return;
    }

    void computeDCurrentNormalVectorDGradU(const floatVector &normalVector, const floatVector &gradU,
                                           floatVector &dNormalVectordGradU, const bool isCurrent) {
        /*!
         * Compute the derivative of the normal vector in the current configuration w.r.t. the displacement gradient
         *
         * \param &normalVector: The unit normal vector in the current configuration
         * \param &gradU: The displacement gradient
         * \param &dNormalVectordGradU: The derivative of the normal vector w.r.t. the displacement gradient
         * \param &isCurrent: Whether the displacement gradient is with respect to the reference or current
         * configuration
         */

        const unsigned int dim = normalVector.size();

        dNormalVectordGradU = floatVector(dim * dim * dim, 0);

        TARDIGRADE_ERROR_TOOLS_CHECK((dim == 3) || (dim == 2) || (dim == 1),
                                     "The dimension of the deformation gradient is " + std::to_string(dim) +
                                         " but must be 1, 2, or 3");

        if (dim == 3) {
            computeDCurrentNormalVectorDGradU<3>(std::begin(normalVector), std::end(normalVector), std::begin(gradU),
                                                 std::end(gradU), std::begin(dNormalVectordGradU),
                                                 std::end(dNormalVectordGradU));
        } else if (dim == 2) {
            computeDCurrentNormalVectorDGradU<2>(std::begin(normalVector), std::end(normalVector), std::begin(gradU),
                                                 std::end(gradU), std::begin(dNormalVectordGradU),
                                                 std::end(dNormalVectordGradU));
        } else if (dim == 1) {
            computeDCurrentNormalVectorDGradU<1>(std::begin(normalVector), std::end(normalVector), std::begin(gradU),
                                                 std::end(gradU), std::begin(dNormalVectordGradU),
                                                 std::end(dNormalVectordGradU));
        }

        return;
    }

    void computeDCurrentAreaWeightedNormalVectorDGradU(const floatVector &normalVector, const floatVector &gradU,
                                                       floatVector &dAreaWeightedNormalVectordGradU,
                                                       const bool   isCurrent) {
        /*!
         * Compute the derivative of the area weighted normal vector w.r.t. the displacement gradient
         *
         * \f$ \frac{\partial}{\partial u_{i,j}} \left( n_i da \right) \f$
         *
         * Note that if the user passes in the unit normal vector, then the result will be more convenient for the
         * construction of the jacobian of a surface integral in the current configuration.
         *
         * \param &normalVector: The normal vector (a unit vector is likely what is desired)
         * \param &gradU: The displacement gradient
         * \param &dAreaWeightedNormalVectordGradU: The derivative of the area weighted normal vector w.r.t. the
         * displacement gradient \param &isCurrent: Whether the displacement gradient is with respect to the reference
         * or current configuration
         */

        const unsigned int dim = normalVector.size();

        dAreaWeightedNormalVectordGradU = floatVector(dim * dim * dim, 0);

        TARDIGRADE_ERROR_TOOLS_CHECK((dim == 3) || (dim == 2) || (dim == 1),
                                     "The dimension of the deformation gradient is " + std::to_string(dim) +
                                         " but must be 1, 2, or 3");

        if (dim == 3) {
            computeDCurrentAreaWeightedNormalVectorDGradU<3>(std::begin(normalVector), std::end(normalVector),
                                                             std::begin(gradU), std::end(gradU),
                                                             std::begin(dAreaWeightedNormalVectordGradU),
                                                             std::end(dAreaWeightedNormalVectordGradU));
        } else if (dim == 2) {
            computeDCurrentAreaWeightedNormalVectorDGradU<2>(std::begin(normalVector), std::end(normalVector),
                                                             std::begin(gradU), std::end(gradU),
                                                             std::begin(dAreaWeightedNormalVectordGradU),
                                                             std::end(dAreaWeightedNormalVectordGradU));
        } else if (dim == 1) {
            computeDCurrentAreaWeightedNormalVectorDGradU<1>(std::begin(normalVector), std::end(normalVector),
                                                             std::begin(gradU), std::end(gradU),
                                                             std::begin(dAreaWeightedNormalVectordGradU),
                                                             std::end(dAreaWeightedNormalVectordGradU));
        }

        return;
    }

    void computeDCurrentAreaDGradU(const floatVector &normalVector, const floatVector &gradU,
                                   floatVector &dCurrentAreadGradU, const bool isCurrent) {
        /*!
         * Compute the derivative of the current area w.r.t. the displacement gradient
         *
         * \param &normalVector: The current unit normal vector
         * \param &gradU: The displacement gradient
         * \param &dCurrentAreadGradU: The derivative of the current surface area w.r.t. the displacement gradient
         * \param &isCurrent: Whether the displacement gradient is with respect to the reference or current
         * configuration
         */

        const unsigned int dim = normalVector.size();

        dCurrentAreadGradU = floatVector(dim * dim, 0);

        TARDIGRADE_ERROR_TOOLS_CHECK((dim == 3) || (dim == 2) || (dim == 1),
                                     "The dimension of the deformation gradient is " + std::to_string(dim) +
                                         " but must be 1, 2, or 3");

        if (dim == 3) {
            computeDCurrentAreaDGradU<3>(std::begin(normalVector), std::end(normalVector), std::begin(gradU),
                                         std::end(gradU), std::begin(dCurrentAreadGradU), std::end(dCurrentAreadGradU));
        } else if (dim == 2) {
            computeDCurrentAreaDGradU<2>(std::begin(normalVector), std::end(normalVector), std::begin(gradU),
                                         std::end(gradU), std::begin(dCurrentAreadGradU), std::end(dCurrentAreadGradU));
        } else if (dim == 1) {
            computeDCurrentAreaDGradU<1>(std::begin(normalVector), std::end(normalVector), std::begin(gradU),
                                         std::end(gradU), std::begin(dCurrentAreadGradU), std::end(dCurrentAreadGradU));
        }

        return;
        //
        //
        //
        //
        //
        //
        //
        //        constexpr unsigned int dim = 3;
        //        constexpr unsigned int sot_dim = dim * dim;
        //
        //        floatVector F;
        //
        //        floatVector dFdGradU;
        //
        //        computeDeformationGradient( gradU, F, dFdGradU, isCurrent );
        //
        //        floatVector dCurrentAreadF;
        //
        //        computeDCurrentAreaDF( normalVector, F, dCurrentAreadF );
        //
        //        dCurrentAreadGradU = floatVector( sot_dim, 0 );
        //
        //        for ( unsigned int i = 0; i < sot_dim; i++ ){
        //
        //            for ( unsigned int j = 0; j < sot_dim; j++ ){
        //
        //                dCurrentAreadGradU[ j ] += dCurrentAreadF[ i ] * dFdGradU[ sot_dim * i + j ];
        //
        //            }
        //
        //        }
    }

}  // namespace tardigradeConstitutiveTools
