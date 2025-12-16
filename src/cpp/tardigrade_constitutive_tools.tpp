/**
 *****************************************************************************
 * \file tardigrade_constitutive_tools.h
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
 *
 * Error handling is taken care of using the tardigrade_error_tools library. All
 * functions should return tardigrade_error_tools::node* objects in the case of an error
 * and NULL otherwise. Wrapping functions in try, except can be used to do
 * this efficiently.
 *****************************************************************************
 */

namespace tardigradeConstitutiveTools {

    template <unsigned int dim, class A_iterator, class Q_iterator, class rotatedA_iterator>
    void rotateMatrix(const A_iterator &A_begin, const A_iterator &A_end, const Q_iterator &Q_begin,
                      const Q_iterator &Q_end, rotatedA_iterator rotatedA_begin, rotatedA_iterator rotatedA_end) {
        /*!
         * Rotate a matrix \f$A\f$ using the orthogonal matrix \f$Q\f$ with the form
         *
         * \f$A'_{ij} = Q_{Ii} A_{IJ} Q_{Jj}\f$
         *
         * TODO: Generalize to non square matrices
         *
         * \param &A_begin: The starting iterator of the matrix to be rotated ( \f$A\f$ )
         * \param &A_end: The stopping iterator matrix to be rotated ( \f$A\f$ )
         * \param &Q_begin: The starting iterator of the rotation matrix ( \f$Q\f$Q )
         * \param &Q_end: The stopping iterator of the rotation matrix ( \f$Q\f$Q )
         * \param &rotatedA_begin: The rotated matrix ( \f$A'\f$ )
         * \param &rotatedA_end: The rotated matrix ( \f$A'\f$ )
         */

        TARDIGRADE_ERROR_TOOLS_CHECK((unsigned int)(A_end - A_begin) == dim * dim,
                                     "A has a size of " + std::to_string((unsigned int)(A_end - A_begin)) +
                                         " and must be a square matrix of size " + std::to_string(dim * dim));

        TARDIGRADE_ERROR_TOOLS_CHECK((unsigned int)(Q_end - Q_begin) == dim * dim,
                                     "Q has a size of " + std::to_string((unsigned int)(Q_end - Q_begin)) +
                                         " and must be a square matrix of size " + std::to_string(dim * dim));

        TARDIGRADE_ERROR_TOOLS_CHECK((unsigned int)(rotatedA_end - rotatedA_begin) == dim * dim,
                                     "rotatedA has a size of " +
                                         std::to_string((unsigned int)(rotatedA_end - rotatedA_begin)) +
                                         " and must be a square matrix of size " + std::to_string(dim * dim));

        using rotatedA_type = typename std::iterator_traits<rotatedA_iterator>::value_type;
        std::array<rotatedA_type, dim * dim> temp;
        std::fill(std::begin(temp), std::end(temp), 0);

        for (unsigned int i = 0; i < dim; ++i) {
            for (unsigned int j = 0; j < dim; ++j) {
                for (unsigned int k = 0; k < dim; ++k) {
                    temp[dim * j + k] += *(Q_begin + dim * i + j) * (*(A_begin + dim * i + k));
                }
            }
        }

        for (unsigned int i = 0; i < dim; ++i) {
            for (unsigned int j = 0; j < dim; ++j) {
                for (unsigned int k = 0; k < dim; ++k) {
                    *(rotatedA_begin + dim * i + k) += temp[dim * i + j] * (*(Q_begin + dim * j + k));
                }
            }
        }

        return;
    }

    template <class A_iterator, class Q_iterator, class rotatedA_iterator>
    void rotateMatrix(const A_iterator &A_begin, const A_iterator &A_end, const Q_iterator &Q_begin,
                      const Q_iterator &Q_end, const unsigned int dim, rotatedA_iterator rotatedA_begin,
                      rotatedA_iterator rotatedA_end) {
        /*!
         * Rotate a matrix \f$A\f$ using the orthogonal matrix \f$Q\f$ with the form
         *
         * \f$A'_{ij} = Q_{Ii} A_{IJ} Q_{Jj}\f$
         *
         * TODO: Generalize to non square matrices
         *
         * \param &A_begin: The starting iterator of the matrix to be rotated ( \f$A\f$ )
         * \param &A_end: The stopping iterator matrix to be rotated ( \f$A\f$ )
         * \param &Q_begin: The starting iterator of the rotation matrix ( \f$Q\f$Q )
         * \param &Q_end: The stopping iterator of the rotation matrix ( \f$Q\f$Q )
         * \param dim: The number of rows/columns in the matrices
         * \param &rotatedA_begin: The rotated matrix ( \f$A'\f$ )
         * \param &rotatedA_end: The rotated matrix ( \f$A'\f$ )
         */

        TARDIGRADE_ERROR_TOOLS_CHECK((unsigned int)(A_end - A_begin) == dim * dim,
                                     "A has a size of " + std::to_string((unsigned int)(A_end - A_begin)) +
                                         " and must be a square matrix of size " + std::to_string(dim * dim));

        TARDIGRADE_ERROR_TOOLS_CHECK((unsigned int)(Q_end - Q_begin) == dim * dim,
                                     "Q has a size of " + std::to_string((unsigned int)(Q_end - Q_begin)) +
                                         " and must be a square matrix of size " + std::to_string(dim * dim));

        TARDIGRADE_ERROR_TOOLS_CHECK((unsigned int)(rotatedA_end - rotatedA_begin) == dim * dim,
                                     "rotatedA has a size of " +
                                         std::to_string((unsigned int)(rotatedA_end - rotatedA_begin)) +
                                         " and must be a square matrix of size " + std::to_string(dim * dim));

        using rotatedA_type = typename std::iterator_traits<rotatedA_iterator>::value_type;
        std::vector<rotatedA_type> temp(dim * dim, rotatedA_type());

        for (unsigned int i = 0; i < dim; ++i) {
            for (unsigned int j = 0; j < dim; ++j) {
                for (unsigned int k = 0; k < dim; ++k) {
                    temp[dim * j + k] += *(Q_begin + dim * i + j) * (*(A_begin + dim * i + k));
                }
            }
        }

        for (unsigned int i = 0; i < dim; ++i) {
            for (unsigned int j = 0; j < dim; ++j) {
                for (unsigned int k = 0; k < dim; ++k) {
                    *(rotatedA_begin + dim * i + k) += temp[dim * i + j] * (*(Q_begin + dim * j + k));
                }
            }
        }

        return;
    }

    template <unsigned int dim, class displacementGradient_iterator, class deformationGradient_iterator>
    void computeDeformationGradient(const displacementGradient_iterator &displacementGradient_begin,
                                    const displacementGradient_iterator &displacementGradient_end,
                                    deformationGradient_iterator         deformationGradient_begin,
                                    deformationGradient_iterator deformationGradient_end, const bool isCurrent) {
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
         * \param &displacementGradient_begin: The starting iterator of the gradient of the displacement with respect to
         * either the current or previous position. \param &displacementGradient_end: The stopping iterator of the
         * gradient of the displacement with respect to either the current or previous position. \param
         * &deformationGradient_begin: The starting iterator of the deformation gradient \param
         * &deformationGradient_end: The stopping iterator of the deformation gradient \param &isCurrent: Boolean
         * indicating whether the gradient is taken w.r.t. the current (true) or reference (false) position.
         */

        using deformationGradient_type = typename std::iterator_traits<deformationGradient_iterator>::value_type;

        constexpr unsigned int sot_dim = dim * dim;

        TARDIGRADE_ERROR_TOOLS_CHECK((unsigned int)(displacementGradient_end - displacementGradient_begin) == sot_dim,
                                     "The displacement gradient has " +
                                         std::to_string((unsigned int)(displacementGradient_end -
                                                                       displacementGradient_begin)) +
                                         " elements but must have " + std::to_string(sot_dim));

        TARDIGRADE_ERROR_TOOLS_CHECK((unsigned int)(deformationGradient_end - deformationGradient_begin) == sot_dim,
                                     "The deformation gradient has " +
                                         std::to_string((unsigned int)(deformationGradient_end -
                                                                       deformationGradient_begin)) +
                                         " elements but must have " + std::to_string(sot_dim));

        std::copy(displacementGradient_begin, displacementGradient_end, deformationGradient_begin);

        if (isCurrent) {
            std::transform(deformationGradient_begin, deformationGradient_end, deformationGradient_begin,
                           std::negate<>());

            for (unsigned int i = 0; i < dim; ++i) {
                *(deformationGradient_begin + dim * i + i) += 1;
            }

            Eigen::Map<Eigen::Matrix<deformationGradient_type, dim, dim> > map(&(*deformationGradient_begin));
            map = map.inverse().eval();

        } else {
            for (unsigned int i = 0; i < dim; ++i) {
                *(deformationGradient_begin + dim * i + i) += 1.;
            }
        }
    }

    template <unsigned int dim, class displacementGradient_iterator, class deformationGradient_iterator,
              class dFdGradU_iterator>
    void computeDeformationGradient(const displacementGradient_iterator &displacementGradient_begin,
                                    const displacementGradient_iterator &displacementGradient_end,
                                    deformationGradient_iterator         deformationGradient_begin,
                                    deformationGradient_iterator         deformationGradient_end,
                                    dFdGradU_iterator dFdGradU_begin, dFdGradU_iterator dFdGradU_end,
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
         * \param &displacementGradient_begin: The starting iterator of the gradient of the displacement with respect to
         * either the current or previous position. \param &displacementGradient_end: The stopping iterator of the
         * gradient of the displacement with respect to either the current or previous position. \param
         * &deformationGradient_begin: The starting iterator of the deformation gradient \param
         * &deformationGradient_end: The stopping iterator of the deformation gradient \param &dFdGradU_begin: The
         * starting iterator of the derivative of the deformation gradient w.r.t. the displacement gradient \param
         * &dFdGradU_end: The stopping iterator of the derivative of the deformation gradient w.r.t. the displacement
         * gradient \param &isCurrent: Boolean indicating whether the gradient is taken w.r.t. the current (true) or
         * reference (false) position.
         */

        using deformationGradient_type = typename std::iterator_traits<deformationGradient_iterator>::value_type;

        constexpr unsigned int sot_dim = dim * dim;
        constexpr unsigned int fot_dim = sot_dim * sot_dim;

        TARDIGRADE_ERROR_TOOLS_CHECK((unsigned int)(displacementGradient_end - displacementGradient_begin) == sot_dim,
                                     "The displacement gradient has " +
                                         std::to_string((unsigned int)(displacementGradient_end -
                                                                       displacementGradient_begin)) +
                                         " elements but must have " + std::to_string(sot_dim));

        TARDIGRADE_ERROR_TOOLS_CHECK((unsigned int)(deformationGradient_end - deformationGradient_begin) == sot_dim,
                                     "The deformation gradient has " +
                                         std::to_string((unsigned int)(deformationGradient_end -
                                                                       deformationGradient_begin)) +
                                         " elements but must have " + std::to_string(sot_dim));

        TARDIGRADE_ERROR_TOOLS_CHECK(
            (unsigned int)(dFdGradU_end - dFdGradU_begin) == fot_dim,
            "The gradient of the deformation gradient with respect to the displacement gradient has a size of " +
                std::to_string((unsigned int)(dFdGradU_end - dFdGradU_begin)) + " but should have a size of " +
                std::to_string(fot_dim));

        std::copy(displacementGradient_begin, displacementGradient_end, deformationGradient_begin);

        std::fill(dFdGradU_begin, dFdGradU_end, deformationGradient_type());

        if (isCurrent) {
            std::transform(deformationGradient_begin, deformationGradient_end, deformationGradient_begin,
                           std::negate<>());

            for (unsigned int i = 0; i < dim; ++i) {
                *(deformationGradient_begin + dim * i + i) += 1;
            }

            Eigen::Map<Eigen::Matrix<deformationGradient_type, dim, dim> > map(&(*deformationGradient_begin));
            map = map.inverse().eval();

            for (unsigned int i = 0; i < dim; ++i) {
                for (unsigned int j = 0; j < dim; ++j) {
                    for (unsigned int k = 0; k < dim; ++k) {
                        for (unsigned int l = 0; l < dim; ++l) {
                            *(dFdGradU_begin + dim * sot_dim * i + sot_dim * j + dim * k + l) =
                                (*(deformationGradient_begin + dim * i + k)) *
                                (*(deformationGradient_begin + dim * l + j));
                        }
                    }
                }
            }

        } else {
            for (unsigned int i = 0; i < dim; ++i) {
                *(deformationGradient_begin + dim * i + i) += 1.;
            }

            for (unsigned int i = 0; i < sot_dim; ++i) {
                *(dFdGradU_begin + sot_dim * i + i) += 1;
            }
        }
    }

    template <unsigned int dim, class deformationGradient_iterator, class C_iterator>
    void computeRightCauchyGreen(const deformationGradient_iterator &deformationGradient_begin,
                                 const deformationGradient_iterator &deformationGradient_end, C_iterator C_begin,
                                 C_iterator C_end) {
        /*!
         * Compute the Right Cauchy-Green deformation tensor ( \f$C\f$ )
         *
         * \f$C_{IJ} = F_{iI} F_{iJ}\f$
         *
         * \param &deformationGradient_begin: The starting iterator of the deformation gradient ( \f$F\f$ )
         * \param &deformationGradient_end: The stopping iterator of the deformation gradient ( \f$F\f$ )
         * \param C_begin: The starting iterator of the resulting Right Cauchy-Green deformation tensor ( \f$C\f$ )
         * \param C_end: The stopping iterator of the resulting Right Cauchy-Green deformation tensor ( \f$C\f$ )
         *
         * The deformation gradient is organized as F11, F12, F13, F21, F22, F23, F31, F32, F33
         *
         * The Right Cauchy-Green deformation tensor is organized as C11, C12, C13, C21, C22, C23, C31, C32, C33
         */

        using deformationGradient_type = typename std::iterator_traits<deformationGradient_iterator>::value_type;
        using C_type                   = typename std::iterator_traits<C_iterator>::value_type;

        constexpr unsigned int sot_dim = dim * dim;

        TARDIGRADE_ERROR_TOOLS_CHECK((unsigned int)(deformationGradient_end - deformationGradient_begin) == sot_dim,
                                     "The deformation gradient has a size of " +
                                         std::to_string((unsigned int)(deformationGradient_end -
                                                                       deformationGradient_begin)) +
                                         " but must have a size of " + std::to_string(sot_dim))

        TARDIGRADE_ERROR_TOOLS_CHECK((unsigned int)(C_end - C_begin) == sot_dim,
                                     "The Right Cauchy Green has a size of " +
                                         std::to_string((unsigned int)(C_end - C_begin)) + " but must have a size of " +
                                         std::to_string(sot_dim))

        std::fill(C_begin, C_end, C_type());

        Eigen::Map<const Eigen::Matrix<deformationGradient_type, dim, dim, Eigen::RowMajor> > F(
            &(*deformationGradient_begin), dim, dim);
        Eigen::Map<Eigen::Matrix<C_type, dim, dim, Eigen::RowMajor> > C_map(&(*C_begin), dim, dim);

        C_map = (F.transpose() * F).eval();
    }

    template <unsigned int dim, class deformationGradient_iterator, class C_iterator, class dCdF_iterator>
    void computeRightCauchyGreen(const deformationGradient_iterator &deformationGradient_begin,
                                 const deformationGradient_iterator &deformationGradient_end, C_iterator C_begin,
                                 C_iterator C_end, dCdF_iterator dCdF_begin, dCdF_iterator dCdF_end) {
        /*!
         * Compute the Right Cauchy-Green deformation tensor ( \f$C\f$ ) from the deformation gradient ( \f$F\f$ )
         *
         * \f$C_{IJ} = F_{iI} F_{iJ}\f$
         *
         * \param &deformationGradient_begin: The starting iterator of the deformation gradient ( \f$F\f$ )
         * \param &deformationGradient_end: The stopping iterator of the deformation gradient ( \f$F\f$ )
         * \param &C_begin: The starting iterator of the resulting Right Cauchy-Green deformation tensor ( \f$C\f$ )
         * \param &C_end: The starting iterator of the resulting Right Cauchy-Green deformation tensor ( \f$C\f$ )
         * \param &dCdF_begin: The starting iterator of the Jacobian of the Right Cauchy-Green deformation tensor
         *     with regards to the deformation gradient ( \f$\frac{\partial C}{\partial F}\f$ ).
         * \param &dCdF_end: The stopping iterator of the Jacobian of the Right Cauchy-Green deformation tensor
         *     with regards to the deformation gradient ( \f$\frac{\partial C}{\partial F}\f$ ).
         *
         * The deformation gradient is organized as F11, F12, F13, F21, F22, F23, F31, F32, F33
         *
         * The Right Cauchy-Green deformation tensor is organized as C11, C12, C13, C21, C22, C23, C31, C32, C33
         */

        using dCdF_type = typename std::iterator_traits<dCdF_iterator>::value_type;

        TARDIGRADE_ERROR_TOOLS_CATCH(computeRightCauchyGreen<dim>(deformationGradient_begin, deformationGradient_end,
                                                                  C_begin, C_end));

        TARDIGRADE_ERROR_TOOLS_CHECK((unsigned int)(dCdF_end - dCdF_begin),
                                     "The derivative of the right Cauchy-Green deformation tensor with respect to the "
                                     "deformation gradient has a size of " +
                                         std::to_string((unsigned int)(dCdF_end - dCdF_begin)) +
                                         " but must have a size of " + std::to_string(dim * dim * dim * dim))

        // Assemble the Jacobian
        std::fill(dCdF_begin, dCdF_end, dCdF_type());

        for (unsigned int I = 0; I < dim; ++I) {
            for (unsigned int J = 0; J < dim; ++J) {
                for (unsigned int k = 0; k < dim; ++k) {
                    *(dCdF_begin + dim * dim * dim * I + dim * dim * J + dim * k + I) +=
                        *(deformationGradient_begin + dim * k + J);
                    *(dCdF_begin + dim * dim * dim * I + dim * dim * J + dim * k + J) +=
                        *(deformationGradient_begin + dim * k + I);
                }
            }
        }
    }

    template <unsigned int dim, class deformationGradient_iterator, class E_iterator>
    void computeGreenLagrangeStrain(const deformationGradient_iterator &deformationGradient_begin,
                                    const deformationGradient_iterator &deformationGradient_end, E_iterator E_begin,
                                    E_iterator E_end) {
        /*!
         * Compute the Green-Lagrange strain ( \f$E\f$ ) from the deformation gradient ( \f$F\f$ ). The operation is:
         *
         * \f$E = 0.5 (F_{iI} F_{iJ} - \delta_{IJ})\f$
         *
         * Where \f$F\f$ is the deformation gradient and \f$\delta\f$ is the kronecker delta.
         *
         * \param &deformationGradient_begin: A starting iterator of the deformation gradient ( \f$F\f$ ).
         * \param &deformationGradient_end: A stopping iterator of the deformation gradient ( \f$F\f$ ).
         * \param &E_begin: The starting iterator of the resulting Green-Lagrange strain ( \f$E\f$ ).
         * \param &E_end: The stopping iterator of the resulting Green-Lagrange strain ( \f$E\f$ ).
         *
         * The deformation gradient is organized as  F11, F12, F13, F21, F22, F23, F31, F32, F33
         *
         * The Green-Lagrange strain is organized as E11, E12, E13, E21, E22, E23, E31, E32, E33
         */

        TARDIGRADE_ERROR_TOOLS_CATCH(computeRightCauchyGreen<dim>(deformationGradient_begin, deformationGradient_end,
                                                                  E_begin, E_end););

        for (unsigned int i = 0; i < dim; ++i) {
            *(E_begin + dim * i + i) -= 1;
        }

        std::transform(E_begin, E_end, E_begin, std::bind(std::multiplies<>(), std::placeholders::_1, 0.5));
    }

    template <unsigned int dim, class deformationGradient_iterator, class E_iterator, class dEdF_iterator>
    void computeGreenLagrangeStrain(const deformationGradient_iterator &deformationGradient_begin,
                                    const deformationGradient_iterator &deformationGradient_end, E_iterator E_begin,
                                    E_iterator E_end, dEdF_iterator dEdF_begin, dEdF_iterator dEdF_end) {
        /*!
         * Compute the Green-Lagrange strain ( \f$E\f$ ) from the deformation gradient ( \f$F\f$ ). The operation is:
         *
         * \f$E = 0.5 (F_{iI} F_{iJ} - \delta_{IJ})\f$
         *
         * Where \f$F\f$ is the deformation gradient and \f$\delta\f$ is the kronecker delta.
         *
         * \param &deformationGradient_begin: A starting iterator of the deformation gradient ( \f$F\f$ ).
         * \param &deformationGradient_end: A stopping iterator of the deformation gradient ( \f$F\f$ ).
         * \param &E_begin: The starting iterator of the resulting Green-Lagrange strain ( \f$E\f$ ).
         * \param &E_end: The stopping iterator of the resulting Green-Lagrange strain ( \f$E\f$ ).
         * \param &dEdF_begin: The starting iterator of the jacobian of the Green-Lagrange strain w.r.t. the
         *     deformation gradient ( \f$\frac{\partial E}{\partial F}\f$ )
         * \param &dEdF_end: The stopping iterator of the jacobian of the Green-Lagrange strain w.r.t. the
         *     deformation gradient ( \f$\frac{\partial E}{\partial F}\f$ )
         *
         * The deformation gradient is organized as  F11, F12, F13, F21, F22, F23, F31, F32, F33
         *
         * The Green-Lagrange strain is organized as E11, E12, E13, E21, E22, E23, E31, E32, E33
         */

        TARDIGRADE_ERROR_TOOLS_CATCH(computeRightCauchyGreen<dim>(deformationGradient_begin, deformationGradient_end,
                                                                  E_begin, E_end, dEdF_begin, dEdF_end));

        for (unsigned int i = 0; i < dim; ++i) {
            *(E_begin + dim * i + i) -= 1;
        }

        std::transform(E_begin, E_end, E_begin, std::bind(std::multiplies<>(), std::placeholders::_1, 0.5));

        std::transform(dEdF_begin, dEdF_end, dEdF_begin, std::bind(std::multiplies<>(), std::placeholders::_1, 0.5));
    }

    template <unsigned int dim, class deformationGradient_iterator, class dEdF_iterator>
    void computeDGreenLagrangeStrainDF(const deformationGradient_iterator &deformationGradient_begin,
                                       const deformationGradient_iterator &deformationGradient_end,
                                       dEdF_iterator dEdF_begin, dEdF_iterator dEdF_end) {
        /*!
         * Compute the derivative of the Green-Lagrange strain ( \f$E\f$ )w.r.t. the deformation gradient ( \f$F\f$ ).
         *
         * \f$\frac{\partial E_{IJ}}{\partial F_{kK}} = 0.5 ( \delta_{IK} F_{kJ} + F_{kI} \delta_{JK})\f$
         *
         * Where \f$F\f$ is the deformation gradient and \f$\delta\f$ is the kronecker delta.
         *
         * \param &deformationGradient_begin: The starting iterator of the deformation gradient ( \f$F\f$ ).
         * \param &deformationGradient_end: The stopping iterator of the deformation gradient ( \f$F\f$ ).
         * \param &dEdF_begin: The starting iterator of the resulting gradient ( \f$\frac{\partial E}{\partial F}\f$ ).
         * \param &dEdF_end: The stopping iterator resulting gradient ( \f$\frac{\partial E}{\partial F}\f$ ).
         *
         * The deformation gradient is organized as  F11, F12, F13, F21, F22, F23, F31, F32, F33
         */

        using dEdF_type = typename std::iterator_traits<dEdF_iterator>::value_type;
        std::fill(dEdF_begin, dEdF_end, dEdF_type());

        for (unsigned int I = 0; I < dim; ++I) {
            for (unsigned int J = 0; J < dim; ++J) {
                for (unsigned int k = 0; k < dim; ++k) {
                    *(dEdF_begin + dim * dim * dim * I + dim * dim * J + dim * k + I) +=
                        *(deformationGradient_begin + dim * k + J);
                    *(dEdF_begin + dim * dim * dim * I + dim * dim * J + dim * k + J) +=
                        *(deformationGradient_begin + dim * k + I);
                }
            }
        }

        std::transform(dEdF_begin, dEdF_end, dEdF_begin, std::bind(std::multiplies<>(), std::placeholders::_1, 0.5));
    }

    template <unsigned int dim, class E_iterator, class Ebar_iterator, typename J_type>
    void decomposeGreenLagrangeStrain(const E_iterator &E_begin, const E_iterator &E_end, Ebar_iterator Ebar_begin,
                                      Ebar_iterator Ebar_end, J_type &J) {
        /*!
         * Decompose the Green-Lagrange strain tensor ( \f$E\f$ ) into isochoric ( \f$\bar{E}\f$ ) and volumetric (
         * \f$J\f$ ) parts where
         *
         * \f$J = det(F) = sqrt(det(2*E + I))\f$
         *
         * \f$\bar{E}_{IJ} = 0.5*((1/(J**(2/3))) F_{iI} F_{iJ} - I_{IJ}) = (1/(J**(2/3)))*E_{IJ} + 0.5(1/(J**(2/3)) -
         * 1)*I_{IJ}\f$
         *
         * \param &E_begin: The starting iterator of the Green-Lagrange strain tensor ( \f$E\f$ )
         * \param &E_end: The stopping iterator of the Green-Lagrange strain tensor ( \f$E\f$ )
         * \param &Ebar_begin: The starting iterator of the isochoric Green-Lagrange strain tensor ( \f$\bar{E}\f$ ).
         *     format = E11, E12, E13, E21, E22, E23, E31, E32, E33
         * \param &Ebar_end: The stopping iterator of the isochoric Green-Lagrange strain tensor ( \f$\bar{E}\f$ ).
         *     format = E11, E12, E13, E21, E22, E23, E31, E32, E33
         * \param &J: The Jacobian of deformation ( \f$J\f$ )
         */

        using E_type = typename std::iterator_traits<E_iterator>::value_type;

        TARDIGRADE_ERROR_TOOLS_CHECK((unsigned int)(E_end - E_begin),
                                     "The Green-Lagrange strain tensor has a size of " +
                                         std::to_string((unsigned int)(E_end - E_begin)) +
                                         " but should have a size of " + std::to_string(dim * dim));

        TARDIGRADE_ERROR_TOOLS_CHECK((unsigned int)(Ebar_end - Ebar_begin),
                                     "The isochoric Green-Lagrange strain tensor has a size of " +
                                         std::to_string((unsigned int)(Ebar_end - Ebar_begin)) +
                                         " but should have a size of " + std::to_string(dim * dim));

        std::array<E_type, dim * dim> F_squared = {E_type()};

        std::transform(E_begin, E_end, std::begin(F_squared),
                       std::bind(std::multiplies<>(), std::placeholders::_1, 2.0));

        for (unsigned int i = 0; i < dim; ++i) {
            F_squared[dim * i + i] += 1;
        }

        Eigen::Map<Eigen::Matrix<E_type, dim, dim, Eigen::RowMajor> > F_squared_map(F_squared.data());
        J_type                                                        Jsq = F_squared_map.determinant();

        TARDIGRADE_ERROR_TOOLS_CHECK(Jsq > 0, "the determinant of the Green-Lagrange strain is negative")

        J = std::sqrt(Jsq);

        std::transform(E_begin, E_end, Ebar_begin,
                       std::bind(std::divides<>(), std::placeholders::_1, std::pow(J, 2. / 3)));

        for (unsigned int i = 0; i < dim; ++i) {
            *(Ebar_begin + dim * i + i) += 0.5 * (1 / std::pow(J, 2. / 3) - 1);
        }
    }

    template <unsigned int dim, class E_iterator, class Ebar_iterator, typename J_type, class dEbardE_iterator,
              class dJdE_iterator>
    void decomposeGreenLagrangeStrain(const E_iterator &E_begin, const E_iterator &E_end, Ebar_iterator Ebar_begin,
                                      Ebar_iterator Ebar_end, J_type &J, dEbardE_iterator dEbardE_begin,
                                      dEbardE_iterator dEbardE_end, dJdE_iterator dJdE_begin, dJdE_iterator dJdE_end) {
        /*!
         * Decompose the Green-Lagrange strain tensor ( \f$E\f$ ) into isochoric ( \f$\bar{E}\f$ ) and volumetric (
         * \f$J\f$ ) parts where
         *
         * \f$J = det(F) = sqrt(det(2*E + I))\f$
         *
         * \f$\bar{E}_{IJ} = 0.5*((1/(J**(2/3))) F_{iI} F_{iJ} - I_{IJ}) = (1/(J**(2/3)))*E_{IJ} + 0.5(1/(J**(2/3)) -
         * 1)*I_{IJ}\f$
         *
         * \param &E_begin: The starting iterator of the Green-Lagrange strain tensor ( \f$E\f$ )
         * \param &E_end: The stopping iterator of the Green-Lagrange strain tensor ( \f$E\f$ )
         * \param &Ebar_begin: The starting iterator of the isochoric Green-Lagrange strain tensor ( \f$\bar{E}\f$ ).
         *     format = E11, E12, E13, E21, E22, E23, E31, E32, E33
         * \param &Ebar_end: The stopping iterator of the isochoric Green-Lagrange strain tensor ( \f$\bar{E}\f$ ).
         *     format = E11, E12, E13, E21, E22, E23, E31, E32, E33
         * \param &J: The Jacobian of deformation ( \f$J\f$ )
         * \param &dEbardE_begin: The starting iterator derivative of the isochoric Green-Lagrange strain
         *     tensor w.r.t. the total strain tensor ( \f$\frac{\partial \bar{E}}{\partial E}\f$ ).
         * \param &dEbardE_end: The stopping iterator derivative of the isochoric Green-Lagrange strain
         *     tensor w.r.t. the total strain tensor ( \f$\frac{\partial \bar{E}}{\partial E}\f$ ).
         * \param &dJdE_begin: The starting iterator derivative of the jacobian of deformation w.r.t. the
         *     Green-Lagrange strain tensor ( \f$\frac{\partial J}{\partial E}\f$ ).
         * \param &dJdE_end: The stopping iterator derivative of the jacobian of deformation w.r.t. the
         *     Green-Lagrange strain tensor ( \f$\frac{\partial J}{\partial E}\f$ ).
         */

        using E_type       = typename std::iterator_traits<E_iterator>::value_type;
        using dJdE_type    = typename std::iterator_traits<dJdE_iterator>::value_type;
        using dEbardE_type = typename std::iterator_traits<dEbardE_iterator>::value_type;

        TARDIGRADE_ERROR_TOOLS_CHECK((unsigned int)(dEbardE_end - dEbardE_begin) == dim * dim * dim * dim,
                                     "The derivative of the isochoric Green-Lagrange strain with respect to the "
                                     "Green-Lagrange strain has a size of " +
                                         std::to_string((unsigned int)(dEbardE_end - dEbardE_begin)) +
                                         " but should have a size of " + std::to_string(dim * dim * dim * dim));

        TARDIGRADE_ERROR_TOOLS_CHECK(
            (unsigned int)(dJdE_end - dJdE_begin) == dim * dim,
            "The derivative of the Jacobian with respect to the Green-Lagrange strain has a size of " +
                std::to_string((unsigned int)(dJdE_end - dJdE_begin)) + " but should have a size of " +
                std::to_string(dim * dim));

        TARDIGRADE_ERROR_TOOLS_CATCH(decomposeGreenLagrangeStrain<dim>(E_begin, E_end, Ebar_begin, Ebar_end, J));

        // Compute the derivative of the jacobian of deformation w.r.t. the Green-Lagrange strain
        std::array<dJdE_type, dim * dim> dJdE_inv;

        std::transform(E_begin, E_end, std::begin(dJdE_inv), std::bind(std::multiplies<>(), std::placeholders::_1, 2));

        for (unsigned int i = 0; i < dim; ++i) {
            dJdE_inv[dim * i + i] += 1;
        }

        Eigen::Map<Eigen::Matrix<E_type, dim, dim, Eigen::RowMajor> >    dJdE_inv_map(dJdE_inv.data());
        Eigen::Map<Eigen::Matrix<dJdE_type, dim, dim, Eigen::RowMajor> > dJdE_map(&(*dJdE_begin));

        dJdE_map = (J * dJdE_inv_map.inverse()).eval();

        // Compute the derivative of the isochoric part of the Green-Lagrange strain w.r.t. the Green-Lagrange strain
        J_type invJ23 = std::pow(J, -2. / 3);
        J_type invJ53 = std::pow(J, -5. / 3);

        std::fill(dEbardE_begin, dEbardE_end, dEbardE_type());

        for (unsigned int i = 0; i < dim * dim; ++i) {
            *(dEbardE_begin + dim * dim * i + i) += invJ23;
        }

        for (unsigned int i = 0; i < dim; i++) {
            for (unsigned int j = 0; j < dim; j++) {
                for (unsigned int k = 0; k < dim; k++) {
                    *(dEbardE_begin + dim * dim * dim * i + dim * dim * i + dim * j + k) -=
                        (1. / 3) * invJ53 * (*(dJdE_begin + dim * j + k));

                    for (unsigned int l = 0; l < dim; l++) {
                        *(dEbardE_begin + dim * dim * dim * i + dim * dim * j + dim * k + l) -=
                            (2. / 3) * invJ53 * (*(E_begin + dim * i + j)) * (*(dJdE_begin + dim * k + l));
                    }
                }
            }
        }
    }

    template <unsigned int dim, class PK2Stress_iterator, class deformationGradient_iterator,
              class cauchyStress_iterator>
    void mapPK2toCauchy(const PK2Stress_iterator &PK2Stress_begin, const PK2Stress_iterator &PK2Stress_end,
                        const deformationGradient_iterator &deformationGradient_begin,
                        const deformationGradient_iterator &deformationGradient_end,
                        cauchyStress_iterator cauchyStress_begin, cauchyStress_iterator cauchyStress_end) {
        /*!
         * Map the PK2 stress ( \f$P^{II}\f$ ) to the current configuration resulting in the Cauchy stress (
         * \f$\sigma\f$ ).
         *
         * \f$\sigma_{ij} = (1/det(F)) F_{iI} P^{II}_{IJ} F_{jJ}\f$
         *
         * where \f$F\f$ is the deformation gradient
         *
         * \param &PK2Stress_begin: The starting iterator of the Second Piola-Kirchoff stress ( \f$P^{II}\f$ )
         * \param &PK2Stress_end: The stopping iterator of the Second Piola-Kirchoff stress ( \f$P^{II}\f$ )
         * \param &deformationGradient_begin: The starting iterator total deformation gradient ( \f$F\f$ ).
         * \param &deformationGradient_end: The stopping iterator total deformation gradient ( \f$F\f$ ).
         * \param cauchyStress_begin: The starting iterator of the Cauchy stress (\f$\sigma\f$ ).
         * \param cauchyStress_end: The stopping iterator of the Cauchy stress (\f$\sigma\f$ ).
         */

        using F_type            = typename std::iterator_traits<deformationGradient_iterator>::value_type;
        using cauchyStress_type = typename std::iterator_traits<cauchyStress_iterator>::value_type;

        TARDIGRADE_ERROR_TOOLS_CHECK((unsigned int)(PK2Stress_end - PK2Stress_begin) == dim * dim,
                                     "The PK2 stress has a size of " +
                                         std::to_string((unsigned int)(PK2Stress_end - PK2Stress_begin)) +
                                         " but must have a size of " + std::to_string(dim * dim));

        TARDIGRADE_ERROR_TOOLS_CHECK((unsigned int)(deformationGradient_end - deformationGradient_begin) == dim * dim,
                                     "The deformation gradient has a size of " +
                                         std::to_string((unsigned int)(deformationGradient_end -
                                                                       deformationGradient_begin)) +
                                         " but must have a size of " + std::to_string(dim * dim));

        // Compute the determinant of the deformation gradient
        Eigen::Map<const Eigen::Matrix<F_type, dim, dim, Eigen::RowMajor> > map(&(*deformationGradient_begin));
        F_type                                                              detF = map.determinant();

        // Initialize the Cauchy stress
        std::array<cauchyStress_type, dim * dim> temp_sot;
        std::fill(std::begin(temp_sot), std::end(temp_sot), cauchyStress_type());
        std::fill(cauchyStress_begin, cauchyStress_end, cauchyStress_type());

        for (unsigned int i = 0; i < dim; i++) {
            for (unsigned int I = 0; I < dim; I++) {
                for (unsigned int j = 0; j < dim; j++) {
                    temp_sot[dim * i + j] +=
                        (*(deformationGradient_begin + dim * i + I)) * (*(PK2Stress_begin + dim * I + j));
                }
            }
        }

        for (unsigned int i = 0; i < dim; i++) {
            for (unsigned int j = 0; j < dim; j++) {
                for (unsigned int I = 0; I < dim; I++) {
                    *(cauchyStress_begin + dim * i + j) +=
                        temp_sot[dim * i + I] * (*(deformationGradient_begin + dim * j + I));
                }
            }
        }

        std::transform(cauchyStress_begin, cauchyStress_end, cauchyStress_begin,
                       std::bind(std::divides<>(), std::placeholders::_1, detF));

        return;
    }

    template <typename temperature_type, class WLFParameters_iterator, typename factor_type>
    void WLF(const temperature_type &temperature, const WLFParameters_iterator &WLFParameters_begin,
             const WLFParameters_iterator &WLFParameters_end, factor_type &factor) {
        /*!
         * An implementation of the Williams-Landel-Ferry equation.
         *
         * \f$factor = 10**((-C_1*(T - T_r))/(C_2 + T - T_r))\f$
         *
         * where \f$T\f$ is the temperature, \f$T_r\f$ is the reference temperature, and \f$C_1\f$ and \f$C_2\f$ are
         * parameters
         *
         * \param &temperature: The temperature \f$T\f$
         * \param &WLFParameters_begin: The starting iterator of the parameters for the function [\f$T_r\f$, \f$C_1\f$,
         * \f$C_2\f$] \param &WLFParameters_end: The stopping iterator of the parameters for the function [\f$T_r\f$,
         * \f$C_1\f$, \f$C_2\f$] \param &factor: The shift factor
         */

        TARDIGRADE_ERROR_TOOLS_CHECK((unsigned int)(WLFParameters_end - WLFParameters_begin) == 3,
                                     "The parameter array has a size of " +
                                         std::to_string((unsigned int)(WLFParameters_end - WLFParameters_begin)) +
                                         " but should have a size of 3");

        factor_type denominator = *(WLFParameters_begin + 2) + (temperature - *(WLFParameters_begin + 0));

        TARDIGRADE_ERROR_TOOLS_CHECK(std::fabs(denominator) > 1e-9, "Zero in the denominator");

        factor =
            std::pow(10., -(*(WLFParameters_begin + 1)) * (temperature - *(WLFParameters_begin + 0)) / denominator);

        return;
    }

    template <typename temperature_type, class WLFParameters_iterator, typename factor_type, typename dFactordT_type>
    void WLF(const temperature_type &temperature, const WLFParameters_iterator &WLFParameters_begin,
             const WLFParameters_iterator &WLFParameters_end, factor_type &factor, dFactordT_type &dFactordT) {
        /*!
         * An implementation of the Williams-Landel-Ferry equation.
         *
         * \f$factor = 10**((-C_1*(T - T_r))/(C_2 + T - T_r))\f$
         *
         * where \f$T\f$ is the temperature, \f$T_r\f$ is the reference temperature, and \f$C_1\f$ and \f$C_2\f$ are
         * parameters
         *
         * \param &temperature: The temperature \f$T\f$
         * \param &WLFParameters_begin: The starting iterator of the parameters for the function [\f$T_r\f$, \f$C_1\f$,
         * \f$C_2\f$] \param &WLFParameters_end: The stopping iterator of the parameters for the function [\f$T_r\f$,
         * \f$C_1\f$, \f$C_2\f$] \param &factor: The shift factor \param &dFactordT: The derivative of the shift factor
         * w.r.t. the temperature
         */

        TARDIGRADE_ERROR_TOOLS_CATCH(WLF(temperature, WLFParameters_begin, WLFParameters_end, factor));

        dFactordT =
            std::log(10) * factor *
            (-*(WLFParameters_begin + 1) / (*(WLFParameters_begin + 2) + temperature - *(WLFParameters_begin + 0)) +
             *(WLFParameters_begin + 1) * (temperature - *(WLFParameters_begin + 0)) /
                 std::pow(*(WLFParameters_begin + 2) + temperature - *(WLFParameters_begin + 0), 2.));
    }

    template <unsigned int dim, class velocityGradient_iterator, class deformationGradient_iterator,
              class DFDt_iterator>
    void computeDFDt(const velocityGradient_iterator    &velocityGradient_begin,
                     const velocityGradient_iterator    &velocityGradient_end,
                     const deformationGradient_iterator &deformationGradient_begin,
                     const deformationGradient_iterator &deformationGradient_end, DFDt_iterator DFDt_begin,
                     DFDt_iterator DFDt_end) {
        /*!
         * Compute the total time derivative of the deformation gradient.
         *
         * \f$\dot{F}_{iI} = L_{ij} F_{jI}\f$
         *
         * \param &velocityGradient_begin: The starting iterator of the velocity gradient \f$L_{ij}\f$
         * \param &velocityGradient_end: The stopping iterator of the velocity gradient \f$L_{ij}\f$
         * \param &deformationGradient_begin: The starting iterator of the deformation gradient \f$F_{iI}\f$
         * \param &deformationGradient_end: The stopping iterator of the deformation gradient \f$F_{iI}\f$
         * \param &DFDt_begin: The starting iterator of the total time derivative of the deformation gradient
         * \param &DFDt_end: The stopping iterator of the total time derivative of the deformation gradient
         */

        using DFDt_type = typename std::iterator_traits<DFDt_iterator>::value_type;

        TARDIGRADE_ERROR_TOOLS_CHECK((unsigned int)(velocityGradient_end - velocityGradient_begin) == dim * dim,
                                     "The velocity gradient has a size of " +
                                         std::to_string((unsigned int)(velocityGradient_end - velocityGradient_begin)) +
                                         " and must have a size of " + std::to_string(dim * dim));

        TARDIGRADE_ERROR_TOOLS_CHECK((unsigned int)(deformationGradient_end - deformationGradient_begin) == dim * dim,
                                     "The deformation gradient has a size of " +
                                         std::to_string((unsigned int)(deformationGradient_end -
                                                                       deformationGradient_begin)) +
                                         " and must have a size of " + std::to_string(dim * dim));

        TARDIGRADE_ERROR_TOOLS_CHECK((unsigned int)(DFDt_end - DFDt_begin) == dim * dim,
                                     "The total temporal derivative of the deformation gradient has a size of " +
                                         std::to_string((unsigned int)(DFDt_end - DFDt_begin)) +
                                         " and must have a size of " + std::to_string(dim * dim));

        std::fill(DFDt_begin, DFDt_end, DFDt_type());

        for (unsigned int i = 0; i < dim; ++i) {
            for (unsigned int I = 0; I < dim; ++I) {
                for (unsigned int j = 0; j < dim; ++j) {
                    *(DFDt_begin + dim * i + j) +=
                        (*(velocityGradient_begin + dim * i + I)) * (*(deformationGradient_begin + dim * I + j));
                }
            }
        }

        return;
    }

    template <unsigned int dim, class velocityGradient_iterator, class deformationGradient_iterator,
              class DFDt_iterator, class dDFDtdL_iterator, class dDFDtdF_iterator>
    void computeDFDt(const velocityGradient_iterator    &velocityGradient_begin,
                     const velocityGradient_iterator    &velocityGradient_end,
                     const deformationGradient_iterator &deformationGradient_begin,
                     const deformationGradient_iterator &deformationGradient_end, DFDt_iterator DFDt_begin,
                     DFDt_iterator DFDt_end, dDFDtdL_iterator dDFDtdL_begin, dDFDtdL_iterator dDFDtdL_end,
                     dDFDtdF_iterator dDFDtdF_begin, dDFDtdF_iterator dDFDtdF_end) {
        /*!
         * Compute the total time derivative of the deformation gradient.
         *
         * \f$\dot{F}_{iI} = L_{ij} F_{jI}\f$
         *
         * \param &velocityGradient_begin: The starting iterator of the velocity gradient \f$L_{ij}\f$
         * \param &velocityGradient_end: The stopping iterator of the velocity gradient \f$L_{ij}\f$
         * \param &deformationGradient_begin: The starting iterator of the deformation gradient \f$F_{iI}\f$
         * \param &deformationGradient_end: The stopping iterator of the deformation gradient \f$F_{iI}\f$
         * \param &DFDt_begin: The starting iterator of the total time derivative of the deformation gradient
         * \param &DFDt_end: The stopping iterator of the total time derivative of the deformation gradient
         * \param &dDFDtdL_begin: The starting iterator of the derivative of the total time derivative of the
         * deformation gradient with respect to the velocity gradient. \param &dDFDtdL_end: The stopping iterator of the
         * derivative of the total time derivative of the deformation gradient with respect to the velocity gradient.
         * \param &dDFDtdF_begin: The starting iterator of the derivative of the total time derivative of the
         * deformation gradient with respect to the deformation gradient. \param &dDFDtdF_end: The stopping iterator of
         * the derivative of the total time derivative of the deformation gradient with respect to the deformation
         * gradient.
         */

        using dDFDtdL_type = typename std::iterator_traits<dDFDtdL_iterator>::value_type;
        using dDFDtdF_type = typename std::iterator_traits<dDFDtdF_iterator>::value_type;

        TARDIGRADE_ERROR_TOOLS_CHECK((unsigned int)(dDFDtdL_end - dDFDtdL_begin) == dim * dim * dim * dim,
                                     "The derivative of the total temporal derivative of the deformation gradient with "
                                     "respect to the velocity gradient has a size of " +
                                         std::to_string((unsigned int)(dDFDtdL_end - dDFDtdL_begin)) +
                                         " and must have a size of " + std::to_string(dim * dim * dim * dim));

        TARDIGRADE_ERROR_TOOLS_CHECK((unsigned int)(dDFDtdF_end - dDFDtdF_begin) == dim * dim * dim * dim,
                                     "The derivative of the total temporal derivative of the deformation gradient with "
                                     "respect to the deformation gradient has a size of " +
                                         std::to_string((unsigned int)(dDFDtdF_end - dDFDtdF_begin)) +
                                         " and must have a size of " + std::to_string(dim * dim * dim * dim));

        TARDIGRADE_ERROR_TOOLS_CATCH(computeDFDt<dim>(velocityGradient_begin, velocityGradient_end,
                                                      deformationGradient_begin, deformationGradient_end, DFDt_begin,
                                                      DFDt_end));

        std::fill(dDFDtdL_begin, dDFDtdL_end, dDFDtdL_type());
        std::fill(dDFDtdF_begin, dDFDtdF_end, dDFDtdF_type());

        for (unsigned int i = 0; i < dim; ++i) {
            for (unsigned int I = 0; I < dim; ++I) {
                for (unsigned int k = 0; k < dim; ++k) {
                    *(dDFDtdL_begin + dim * dim * dim * i + dim * dim * I + dim * i + k) =
                        *(deformationGradient_begin + dim * k + I);
                    *(dDFDtdF_begin + dim * dim * dim * i + dim * dim * I + dim * k + I) =
                        *(velocityGradient_begin + dim * i + k);
                }
            }
        }

        return;
    }

    template <typename Dt_type, class Ap_iterator, class DApDt_iterator, class DADt_iterator, class dA_iterator,
              class A_iterator, class alpha_iterator>
    void midpointEvolution(const Dt_type &Dt, const Ap_iterator &Ap_begin, const Ap_iterator &Ap_end,
                           const DApDt_iterator &DApDt_begin, const DApDt_iterator &DApDt_end,
                           const DADt_iterator &DADt_begin, const DADt_iterator &DADt_end, dA_iterator dA_begin,
                           dA_iterator dA_end, A_iterator A_begin, A_iterator A_end, alpha_iterator alpha_begin,
                           alpha_iterator alpha_end) {
        /*!
         * Perform midpoint rule based evolution of a vector.
         *
         * alpha=0 (implicit)
         *
         * alpha=1 (explicit)
         *
         * \param &Dt: The change in time.
         * \param &Ap_begin: The starting iterator of the previous value of the vector
         * \param &Ap_end: The stopping iterator of the previous value of the vector
         * \param &DApDt_begin: The starting iterator of the previous time rate of change of the vector.
         * \param &DApDt_end: The stopping iterator of the previous time rate of change of the vector.
         * \param &DADt_begin: The starting iterator of the current time rate of change of the vector.
         * \param &DADt_end: The stopping iterator of the current time rate of change of the vector.
         * \param &dA_begin: The starting iterator of the change in the vector
         * \param &dA_end: The stopping iterator of the change in the vector
         * \param &A_begin: The starting iterator of the current value of the vector.
         * \param &A_end: The stopping iterator of the current value of the vector.
         * \param &alpha_begin: The starting iterator of the integration parameter.
         * \param &alpha_end: The stopping iterator of the integration parameter.
         */

        using dA_type = typename std::iterator_traits<dA_iterator>::value_type;
        using A_type  = typename std::iterator_traits<A_iterator>::value_type;

        TARDIGRADE_ERROR_TOOLS_EVAL(unsigned int Ap_size = (unsigned int)(Ap_end - Ap_begin);)

        TARDIGRADE_ERROR_TOOLS_CHECK(Ap_size == (unsigned int)(DApDt_end - DApDt_begin),
                                     "DApDt has a size of " + std::to_string((unsigned int)(DApDt_end - DApDt_begin)) +
                                         " and must be consistent with Ap which has a size of " +
                                         std::to_string(Ap_size));

        TARDIGRADE_ERROR_TOOLS_CHECK(Ap_size == (unsigned int)(DADt_end - DADt_begin),
                                     "DADt has a size of " + std::to_string((unsigned int)(DADt_end - DADt_begin)) +
                                         " and must be consistent with Ap which has a size of " +
                                         std::to_string(Ap_size));

        TARDIGRADE_ERROR_TOOLS_CHECK(Ap_size == (unsigned int)(dA_end - dA_begin),
                                     "dA has a size of " + std::to_string((unsigned int)(dA_end - dA_begin)) +
                                         " and must be consistent with Ap which has a size of " +
                                         std::to_string(Ap_size));

        TARDIGRADE_ERROR_TOOLS_CHECK(Ap_size == (unsigned int)(A_end - A_begin),
                                     "A has a size of " + std::to_string((unsigned int)(A_end - A_begin)) +
                                         " and must be consistent with Ap which has a size of " +
                                         std::to_string(Ap_size));

        TARDIGRADE_ERROR_TOOLS_CHECK(Ap_size == (unsigned int)(alpha_end - alpha_begin),
                                     "alpha has a size of " + std::to_string((unsigned int)(alpha_end - alpha_begin)) +
                                         " and must be consistent with Ap which has a size of " +
                                         std::to_string(Ap_size));

        std::fill(dA_begin, dA_end, dA_type());

        std::fill(A_begin, A_end, A_type());

        for (auto v = std::pair<unsigned int, alpha_iterator>(0, alpha_begin); v.second != alpha_end;
             ++v.first, ++v.second) {
            TARDIGRADE_ERROR_TOOLS_CHECK(((*v.second) >= 0) && ((*v.second) <= 1), "Alpha must be between 0 and 1");

            *(dA_begin + v.first) =
                Dt * (*v.second * (*(DApDt_begin + v.first)) + (1 - *v.second) * (*(DADt_begin + v.first)));

            *(A_begin + v.first) = *(Ap_begin + v.first) + *(dA_begin + v.first);
        }

        return;
    }

    template <typename Dt_type, class Ap_iterator, class DApDt_iterator, class DADt_iterator, class dA_iterator,
              class A_iterator, class DADADt_iterator, class alpha_iterator>
    void midpointEvolution(const Dt_type &Dt, const Ap_iterator &Ap_begin, const Ap_iterator &Ap_end,
                           const DApDt_iterator &DApDt_begin, const DApDt_iterator &DApDt_end,
                           const DADt_iterator &DADt_begin, const DADt_iterator &DADt_end, dA_iterator dA_begin,
                           dA_iterator dA_end, A_iterator A_begin, A_iterator A_end, DADADt_iterator DADADt_begin,
                           DADADt_iterator DADADt_end, alpha_iterator alpha_begin, alpha_iterator alpha_end) {
        /*!
         * Perform midpoint rule based evolution of a vector.
         *
         * alpha=0 (implicit)
         *
         * alpha=1 (explicit)
         *
         * \param &Dt: The change in time.
         * \param &Ap_begin: The starting iterator of the previous value of the vector
         * \param &Ap_end: The stopping iterator of the previous value of the vector
         * \param &DApDt_begin: The starting iterator of the previous time rate of change of the vector.
         * \param &DApDt_end: The stopping iterator of the previous time rate of change of the vector.
         * \param &DADt_begin: The starting iterator of the current time rate of change of the vector.
         * \param &DADt_end: The stopping iterator of the current time rate of change of the vector.
         * \param &dA_begin: The starting iterator of the change in the vector
         * \param &dA_end: The stopping iterator of the change in the vector
         * \param &A_begin: The starting iterator of the current value of the vector.
         * \param &A_end: The stopping iterator of the current value of the vector.
         * \param &DADADt_begin: The starting iterator of the gradient of A w.r.t. the current rate of change.
         * \param &DADADt_end: The stopping iterator of the gradient of A w.r.t. the current rate of change.
         * \param &alpha_begin: The starting iterator of the integration parameter.
         * \param &alpha_end: The stopping iterator of the integration parameter.
         */

        using DADADt_type = typename std::iterator_traits<DADADt_iterator>::value_type;

        const unsigned int Ap_size = (unsigned int)(Ap_end - Ap_begin);

        TARDIGRADE_ERROR_TOOLS_CHECK((unsigned int)(DADADt_end - DADADt_begin) == (Ap_size * Ap_size),
                                     "DADADt has a size of " +
                                         std::to_string((unsigned int)(DADADt_end - DADADt_begin)) +
                                         " but must have a size of " + std::to_string(Ap_size * Ap_size));

        TARDIGRADE_ERROR_TOOLS_CATCH(midpointEvolution(Dt, Ap_begin, Ap_end, DApDt_begin, DApDt_end, DADt_begin,
                                                       DADt_end, dA_begin, dA_end, A_begin, A_end, alpha_begin,
                                                       alpha_end))

        std::fill(DADADt_begin, DADADt_end, DADADt_type());

        for (auto v = std::pair<unsigned int, alpha_iterator>(0, alpha_begin); v.second != alpha_end;
             ++v.first, ++v.second) {
            *(DADADt_begin + Ap_size * v.first + v.first) = Dt * (1 - *v.second);
        }
    }

    template <typename Dt_type, class Ap_iterator, class DApDt_iterator, class DADt_iterator, class dA_iterator,
              class A_iterator, class DADADt_iterator, class DADApDt_iterator, class alpha_iterator>
    void midpointEvolution(const Dt_type &Dt, const Ap_iterator &Ap_begin, const Ap_iterator &Ap_end,
                           const DApDt_iterator &DApDt_begin, const DApDt_iterator &DApDt_end,
                           const DADt_iterator &DADt_begin, const DADt_iterator &DADt_end, dA_iterator dA_begin,
                           dA_iterator dA_end, A_iterator A_begin, A_iterator A_end, DADADt_iterator DADADt_begin,
                           DADADt_iterator DADADt_end, DADApDt_iterator DADApDt_begin, DADApDt_iterator DADApDt_end,
                           alpha_iterator alpha_begin, alpha_iterator alpha_end) {
        /*!
         * Perform midpoint rule based evolution of a vector.
         *
         * alpha=0 (implicit)
         *
         * alpha=1 (explicit)
         *
         * \param &Dt: The change in time.
         * \param &Ap_begin: The starting iterator of the previous value of the vector
         * \param &Ap_end: The stopping iterator of the previous value of the vector
         * \param &DApDt_begin: The starting iterator of the previous time rate of change of the vector.
         * \param &DApDt_end: The stopping iterator of the previous time rate of change of the vector.
         * \param &DADt_begin: The starting iterator of the current time rate of change of the vector.
         * \param &DADt_end: The stopping iterator of the current time rate of change of the vector.
         * \param &dA_begin: The starting iterator of the change in the vector
         * \param &dA_end: The stopping iterator of the change in the vector
         * \param &A_begin: The starting iterator of the current value of the vector.
         * \param &A_end: The stopping iterator of the current value of the vector.
         * \param &DADADt_begin: The starting iterator of the gradient of A w.r.t. the current rate of change.
         * \param &DADADt_end: The stopping iterator of the gradient of A w.r.t. the current rate of change.
         * \param &DADApDt_begin: The starting iterator of the gradient of A w.r.t. the previous rate of change.
         * \param &DADApDt_end: The stopping iterator of the gradient of A w.r.t. the previous rate of change.
         * \param &alpha_begin: The starting iterator of the integration parameter.
         * \param &alpha_end: The stopping iterator of the integration parameter.
         */

        using DADApDt_type = typename std::iterator_traits<DADApDt_iterator>::value_type;

        const unsigned int Ap_size = (unsigned int)(Ap_end - Ap_begin);

        TARDIGRADE_ERROR_TOOLS_CHECK((unsigned int)(DADApDt_end - DADApDt_begin) == (Ap_size * Ap_size),
                                     "DADApDt has a size of " +
                                         std::to_string((unsigned int)(DADApDt_end - DADApDt_begin)) +
                                         " but must have a size of " + std::to_string(Ap_size * Ap_size));

        TARDIGRADE_ERROR_TOOLS_CATCH(midpointEvolution(Dt, Ap_begin, Ap_end, DApDt_begin, DApDt_end, DADt_begin,
                                                       DADt_end, dA_begin, dA_end, A_begin, A_end, DADADt_begin,
                                                       DADADt_end, alpha_begin, alpha_end))

        std::fill(DADApDt_begin, DADApDt_end, DADApDt_type());

        for (auto v = std::pair<unsigned int, alpha_iterator>(0, alpha_begin); v.second != alpha_end;
             ++v.first, ++v.second) {
            *(DADApDt_begin + Ap_size * v.first + v.first) = Dt * (*v.second);
        }
    }

    template <typename Dt_type, class Ap_iterator, class DApDt_iterator, class DADt_iterator, class dA_iterator,
              class A_iterator, typename alpha_type>
    void midpointEvolution(const Dt_type &Dt, const Ap_iterator &Ap_begin, const Ap_iterator &Ap_end,
                           const DApDt_iterator &DApDt_begin, const DApDt_iterator &DApDt_end,
                           const DADt_iterator &DADt_begin, const DADt_iterator &DADt_end, dA_iterator dA_begin,
                           dA_iterator dA_end, A_iterator A_begin, A_iterator A_end, alpha_type alpha) {
        /*!
         * Perform midpoint rule based evolution of a vector.
         *
         * alpha=0 (implicit)
         *
         * alpha=1 (explicit)
         *
         * \param &Dt: The change in time.
         * \param &Ap_begin: The starting iterator of the previous value of the vector
         * \param &Ap_end: The stopping iterator of the previous value of the vector
         * \param &DApDt_begin: The starting iterator of the previous time rate of change of the vector.
         * \param &DApDt_end: The stopping iterator of the previous time rate of change of the vector.
         * \param &DADt_begin: The starting iterator of the current time rate of change of the vector.
         * \param &DADt_end: The stopping iterator of the current time rate of change of the vector.
         * \param &dA_begin: The starting iterator of the change in the vector
         * \param &dA_end: The stopping iterator of the change in the vector
         * \param &A_begin: The starting iterator of the current value of the vector.
         * \param &A_end: The stopping iterator of the current value of the vector.
         * \param &alpha: The integration parameter (defaults to 0.5)
         */

        using dA_type = typename std::iterator_traits<dA_iterator>::value_type;
        using A_type  = typename std::iterator_traits<A_iterator>::value_type;

        TARDIGRADE_ERROR_TOOLS_EVAL(unsigned int Ap_size = (unsigned int)(Ap_end - Ap_begin);)

        TARDIGRADE_ERROR_TOOLS_CHECK(Ap_size == (unsigned int)(DApDt_end - DApDt_begin),
                                     "DApDt has a size of " + std::to_string((unsigned int)(DApDt_end - DApDt_begin)) +
                                         " and must be consistent with Ap which has a size of " +
                                         std::to_string(Ap_size));

        TARDIGRADE_ERROR_TOOLS_CHECK(Ap_size == (unsigned int)(DADt_end - DADt_begin),
                                     "DADt has a size of " + std::to_string((unsigned int)(DADt_end - DADt_begin)) +
                                         " and must be consistent with Ap which has a size of " +
                                         std::to_string(Ap_size));

        TARDIGRADE_ERROR_TOOLS_CHECK(Ap_size == (unsigned int)(dA_end - dA_begin),
                                     "dA has a size of " + std::to_string((unsigned int)(dA_end - dA_begin)) +
                                         " and must be consistent with Ap which has a size of " +
                                         std::to_string(Ap_size));

        TARDIGRADE_ERROR_TOOLS_CHECK(Ap_size == (unsigned int)(A_end - A_begin),
                                     "A has a size of " + std::to_string((unsigned int)(A_end - A_begin)) +
                                         " and must be consistent with Ap which has a size of " +
                                         std::to_string(Ap_size));

        TARDIGRADE_ERROR_TOOLS_CHECK(((alpha >= 0) && (alpha <= 1)), "alpha has a value of " + std::to_string(alpha) +
                                                                         " but it must be between 0 and 1");

        std::fill(dA_begin, dA_end, dA_type());

        std::fill(A_begin, A_end, A_type());

        for (auto v = std::pair<unsigned int, Ap_iterator>(0, Ap_begin); v.second != Ap_end; ++v.first, ++v.second) {
            *(dA_begin + v.first) = Dt * (alpha * (*(DApDt_begin + v.first)) + (1 - alpha) * (*(DADt_begin + v.first)));

            *(A_begin + v.first) = *v.second + *(dA_begin + v.first);
        }

        return;
    }

    template <typename Dt_type, class Ap_iterator, class DApDt_iterator, class DADt_iterator, class dA_iterator,
              class A_iterator, class DADADt_iterator, typename alpha_type>
    void midpointEvolution(const Dt_type &Dt, const Ap_iterator &Ap_begin, const Ap_iterator &Ap_end,
                           const DApDt_iterator &DApDt_begin, const DApDt_iterator &DApDt_end,
                           const DADt_iterator &DADt_begin, const DADt_iterator &DADt_end, dA_iterator dA_begin,
                           dA_iterator dA_end, A_iterator A_begin, A_iterator A_end, DADADt_iterator DADADt_begin,
                           DADADt_iterator DADADt_end, alpha_type alpha) {
        /*!
         * Perform midpoint rule based evolution of a vector.
         *
         * alpha=0 (implicit)
         *
         * alpha=1 (explicit)
         *
         * \param &Dt: The change in time.
         * \param &Ap_begin: The starting iterator of the previous value of the vector
         * \param &Ap_end: The stopping iterator of the previous value of the vector
         * \param &DApDt_begin: The starting iterator of the previous time rate of change of the vector.
         * \param &DApDt_end: The stopping iterator of the previous time rate of change of the vector.
         * \param &DADt_begin: The starting iterator of the current time rate of change of the vector.
         * \param &DADt_end: The stopping iterator of the current time rate of change of the vector.
         * \param &dA_begin: The starting iterator of the change in the vector
         * \param &dA_end: The stopping iterator of the change in the vector
         * \param &A_begin: The starting iterator of the current value of the vector.
         * \param &A_end: The stopping iterator of the current value of the vector.
         * \param &DADADt_begin: The starting iterator of the gradient of A w.r.t. the current rate of change.
         * \param &DADADt_end: The stopping iterator of the gradient of A w.r.t. the current rate of change.
         * \param &alpha: The integration parameter (defaults to 0.5)
         */

        using DADADt_type = typename std::iterator_traits<DADADt_iterator>::value_type;

        const unsigned int Ap_size = (unsigned int)(Ap_end - Ap_begin);

        TARDIGRADE_ERROR_TOOLS_CHECK((unsigned int)(DADADt_end - DADADt_begin) == (Ap_size * Ap_size),
                                     "DADADt has a size of " +
                                         std::to_string((unsigned int)(DADADt_end - DADADt_begin)) +
                                         " but must have a size of " + std::to_string(Ap_size * Ap_size));

        TARDIGRADE_ERROR_TOOLS_CATCH(midpointEvolution(Dt, Ap_begin, Ap_end, DApDt_begin, DApDt_end, DADt_begin,
                                                       DADt_end, dA_begin, dA_end, A_begin, A_end, alpha))

        std::fill(DADADt_begin, DADADt_end, DADADt_type());

        for (auto v = std::pair<unsigned int, Ap_iterator>(0, Ap_begin); v.second != Ap_end; ++v.first, ++v.second) {
            *(DADADt_begin + Ap_size * v.first + v.first) = Dt * (1 - alpha);
        }
    }

    template <typename Dt_type, class Ap_iterator, class DApDt_iterator, class DADt_iterator, class dA_iterator,
              class A_iterator, class DADADt_iterator, class DADApDt_iterator, typename alpha_type>
    void midpointEvolution(const Dt_type &Dt, const Ap_iterator &Ap_begin, const Ap_iterator &Ap_end,
                           const DApDt_iterator &DApDt_begin, const DApDt_iterator &DApDt_end,
                           const DADt_iterator &DADt_begin, const DADt_iterator &DADt_end, dA_iterator dA_begin,
                           dA_iterator dA_end, A_iterator A_begin, A_iterator A_end, DADADt_iterator DADADt_begin,
                           DADADt_iterator DADADt_end, DADApDt_iterator DADApDt_begin, DADApDt_iterator DADApDt_end,
                           alpha_type alpha) {
        /*!
         * Perform midpoint rule based evolution of a vector.
         *
         * alpha=0 (implicit)
         *
         * alpha=1 (explicit)
         *
         * \param &Dt: The change in time.
         * \param &Ap_begin: The starting iterator of the previous value of the vector
         * \param &Ap_end: The stopping iterator of the previous value of the vector
         * \param &DApDt_begin: The starting iterator of the previous time rate of change of the vector.
         * \param &DApDt_end: The stopping iterator of the previous time rate of change of the vector.
         * \param &DADt_begin: The starting iterator of the current time rate of change of the vector.
         * \param &DADt_end: The stopping iterator of the current time rate of change of the vector.
         * \param &dA_begin: The starting iterator of the change in the vector
         * \param &dA_end: The stopping iterator of the change in the vector
         * \param &A_begin: The starting iterator of the current value of the vector.
         * \param &A_end: The stopping iterator of the current value of the vector.
         * \param &DADADt_begin: The starting iterator of the gradient of A w.r.t. the current rate of change.
         * \param &DADADt_end: The stopping iterator of the gradient of A w.r.t. the current rate of change.
         * \param &DADApDt_begin: The starting iterator of the gradient of A w.r.t. the previous rate of change.
         * \param &DADApDt_end: The stopping iterator of the gradient of A w.r.t. the previous rate of change.
         * \param &alpha: The integration parameter (defaults to 0.5)
         */

        using DADApDt_type = typename std::iterator_traits<DADApDt_iterator>::value_type;

        const unsigned int Ap_size = (unsigned int)(Ap_end - Ap_begin);

        TARDIGRADE_ERROR_TOOLS_CHECK((unsigned int)(DADApDt_end - DADApDt_begin) == (Ap_size * Ap_size),
                                     "DADApDt has a size of " +
                                         std::to_string((unsigned int)(DADApDt_end - DADApDt_begin)) +
                                         " but must have a size of " + std::to_string(Ap_size * Ap_size));

        TARDIGRADE_ERROR_TOOLS_CATCH(midpointEvolution(Dt, Ap_begin, Ap_end, DApDt_begin, DApDt_end, DADt_begin,
                                                       DADt_end, dA_begin, dA_end, A_begin, A_end, DADADt_begin,
                                                       DADADt_end, alpha))

        std::fill(DADApDt_begin, DADApDt_end, DADApDt_type());

        for (auto v = std::pair<unsigned int, Ap_iterator>(0, Ap_begin); v.second != Ap_end; ++v.first, ++v.second) {
            *(DADApDt_begin + Ap_size * v.first + v.first) = Dt * alpha;
        }
    }

    template <unsigned int dim, unsigned int mode, typename Dt_type, class previousDeformationGradient_iterator,
              class Lp_iterator, class L_iterator, class dF_iterator, class deformationGradient_iterator,
              typename alpha_type>
    void evolveF(const Dt_type &Dt, const previousDeformationGradient_iterator &previousDeformationGradient_begin,
                 const previousDeformationGradient_iterator &previousDeformationGradient_end,
                 const Lp_iterator &Lp_begin, const Lp_iterator &Lp_end, const L_iterator &L_begin,
                 const L_iterator &L_end, dF_iterator dF_begin, dF_iterator dF_end,
                 deformationGradient_iterator deformationGradient_begin,
                 deformationGradient_iterator deformationGradient_end, const alpha_type alpha) {
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
         * \param &previousDeformationGradient_begin: The starting iterator of the previous value of the deformation
         * gradient \param &previousDeformationGradient_end: The stopping iterator of the previous value of the
         * deformation gradient \param &Lp_begin: The starting iterator previous velocity gradient in the current
         * configuration (mode 1) or reference configuration (mode 2). \param &Lp_end: The stopping iterator previous
         * velocity gradient in the current configuration (mode 1) or reference configuration (mode 2). \param &L_begin:
         * The starting iterator of the current velocity gradient in the current configuration (mode 1) or reference
         * configuration (mode 2). \param &L_end: The stopping iterator of the current velocity gradient in the current
         * configuration (mode 1) or reference configuration (mode 2). \param &dF_begin: The starting iterator of the
         * change in the deformation gradient \f$\Delta \bf{F}\f$ such that \f$F_{iI}^{t+1} = F_{iI}^t + \Delta
         * F_{iI}\f$ \param &dF_end: The stopping iterator of the change in the deformation gradient \f$\Delta \bf{F}\f$
         * such that \f$F_{iI}^{t+1} = F_{iI}^t + \Delta F_{iI}\f$ \param &deformationGradient_begin: The starting
         * iterator of the computed current deformation gradient. \param &deformationGradient_end: The starting iterator
         * of the computed current deformation gradient. \param alpha: The integration parameter.
         */

        using deformationGradient_type = typename std::iterator_traits<deformationGradient_iterator>::value_type;

        TARDIGRADE_ERROR_TOOLS_CHECK((unsigned int)(deformationGradient_end - deformationGradient_begin) == dim * dim,
                                     "The deformation gradient has a size of " +
                                         std::to_string((unsigned int)(deformationGradient_end -
                                                                       deformationGradient_begin)) +
                                         " but it should have a size of " + std::to_string(dim * dim))

        TARDIGRADE_ERROR_TOOLS_CHECK(
            (unsigned int)(previousDeformationGradient_end - previousDeformationGradient_begin) == dim * dim,
            "The previous deformation gradient has a size of " +
                std::to_string((unsigned int)(previousDeformationGradient_end - previousDeformationGradient_begin)) +
                " but it should have a size of " + std::to_string(dim * dim))

        TARDIGRADE_ERROR_TOOLS_CHECK((unsigned int)(Lp_end - Lp_begin) == dim * dim,
                                     "The previous velocity gradient has a size of " +
                                         std::to_string((unsigned int)(Lp_end - Lp_begin)) +
                                         " but it should have a size of " + std::to_string(dim * dim))

        TARDIGRADE_ERROR_TOOLS_CHECK((unsigned int)(L_end - L_begin) == dim * dim,
                                     "The velocity gradient has a size of " +
                                         std::to_string((unsigned int)(L_end - L_begin)) +
                                         " but it should have a size of " + std::to_string(dim * dim))

        TARDIGRADE_ERROR_TOOLS_CHECK((mode == 1) || (mode == 2), "The mode of evolution " + std::to_string(mode) +
                                                                     " is not recognized. It must be 1 or 2.");

        // Construct the velocity gradient at time t + alpha
        std::array<deformationGradient_type, dim * dim> LtpAlpha;
        for (unsigned int i = 0; i < dim * dim; ++i) {
            LtpAlpha[i] = alpha * (*(Lp_begin + i)) + (1 - alpha) * (*(L_begin + i));
        }

        // Compute the right hand side
        std::array<deformationGradient_type, dim * dim> RHS;

        Eigen::Map<const Eigen::Matrix<deformationGradient_type, dim, dim, Eigen::RowMajor> > Fp(
            &(*previousDeformationGradient_begin));
        Eigen::Map<const Eigen::Matrix<deformationGradient_type, dim, dim, Eigen::RowMajor> > Lt(LtpAlpha.data());
        Eigen::Map<Eigen::Matrix<deformationGradient_type, dim, dim, Eigen::RowMajor> > RHS_map(RHS.data(), dim, dim);

        if (mode == 1) {
            RHS_map = (Lt * Fp * Dt).eval();
        }
        if (mode == 2) {
            RHS_map = (Fp * Lt * Dt).eval();
        }

        // Compute the left-hand side
        std::array<deformationGradient_type, dim * dim> LHS;
        std::transform(L_begin, L_end, std::begin(LHS),
                       std::bind(std::multiplies<>(), std::placeholders::_1, -Dt * (1 - alpha)));

        for (unsigned int i = 0; i < dim; i++) {
            LHS[dim * i + i] += 1;
        }

        std::fill(dF_begin, dF_end, deformationGradient_type());

        Eigen::Map<Eigen::Matrix<deformationGradient_type, dim, dim, Eigen::RowMajor> > LHS_map(LHS.data());
        Eigen::Map<Eigen::Matrix<deformationGradient_type, dim, dim, Eigen::RowMajor> > dF_map(&(*dF_begin));

        if (mode == 1) {
            dF_map = (LHS_map.inverse() * RHS_map).eval();
        }
        if (mode == 2) {
            dF_map = (RHS_map * LHS_map.inverse()).eval();
        }

        std::transform(previousDeformationGradient_begin, previousDeformationGradient_end, dF_begin,
                       deformationGradient_begin, std::plus<>());

        return;
    }

    template <unsigned int dim, unsigned int mode, typename Dt_type, class previousDeformationGradient_iterator,
              class Lp_iterator, class L_iterator, class dF_iterator, class deformationGradient_iterator,
              class dFdL_iterator, typename alpha_type>
    void evolveF(const Dt_type &Dt, const previousDeformationGradient_iterator &previousDeformationGradient_begin,
                 const previousDeformationGradient_iterator &previousDeformationGradient_end,
                 const Lp_iterator &Lp_begin, const Lp_iterator &Lp_end, const L_iterator &L_begin,
                 const L_iterator &L_end, dF_iterator dF_begin, dF_iterator dF_end,
                 deformationGradient_iterator deformationGradient_begin,
                 deformationGradient_iterator deformationGradient_end, dFdL_iterator dFdL_begin, dFdL_iterator dFdL_end,
                 const alpha_type alpha) {
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
         * \param &previousDeformationGradient_begin: The starting iterator of the previous value of the deformation
         * gradient \param &previousDeformationGradient_end: The stopping iterator of the previous value of the
         * deformation gradient \param &Lp_begin: The starting iterator previous velocity gradient in the current
         * configuration (mode 1) or reference configuration (mode 2). \param &Lp_end: The stopping iterator previous
         * velocity gradient in the current configuration (mode 1) or reference configuration (mode 2). \param &L_begin:
         * The starting iterator of the current velocity gradient in the current configuration (mode 1) or reference
         * configuration (mode 2). \param &L_end: The stopping iterator of the current velocity gradient in the current
         * configuration (mode 1) or reference configuration (mode 2). \param &dF_begin: The starting iterator of the
         * change in the deformation gradient \f$\Delta \bf{F}\f$ such that \f$F_{iI}^{t+1} = F_{iI}^t + \Delta
         * F_{iI}\f$ \param &dF_end: The stopping iterator of the change in the deformation gradient \f$\Delta \bf{F}\f$
         * such that \f$F_{iI}^{t+1} = F_{iI}^t + \Delta F_{iI}\f$ \param &deformationGradient_begin: The starting
         * iterator of the computed current deformation gradient. \param &deformationGradient_end: The starting iterator
         * of the computed current deformation gradient. \param &dFdL_begin: The starting iterator of the derivative of
         * the deformation gradient w.r.t. the velocity gradient \param &dFdL_end: The stopping iterator of the
         * derivative of the deformation gradient w.r.t. the velocity gradient \param alpha: The integration parameter.
         */

        using deformationGradient_type = typename std::iterator_traits<deformationGradient_iterator>::value_type;

        TARDIGRADE_ERROR_TOOLS_CHECK(
            (unsigned int)(dFdL_end - dFdL_begin) == dim * dim * dim * dim,
            "The derivative of the deformation gradient w.r.t. the velocity gradient has a size of " +
                std::to_string((unsigned int)(dFdL_end - dFdL_begin)) + " but must have a size of " +
                std::to_string(dim * dim * dim * dim));

        evolveF<dim, mode>(Dt, previousDeformationGradient_begin, previousDeformationGradient_end, Lp_begin, Lp_end,
                           L_begin, L_end, dF_begin, dF_end, deformationGradient_begin, deformationGradient_end, alpha);

        // Compute the left-hand side
        std::array<deformationGradient_type, dim * dim> LHS, inv_LHS;
        std::transform(L_begin, L_end, std::begin(LHS),
                       std::bind(std::multiplies<>(), std::placeholders::_1, -Dt * (1 - alpha)));

        for (unsigned int i = 0; i < dim; i++) {
            LHS[dim * i + i] += 1;
        }

        Eigen::Map<Eigen::Matrix<deformationGradient_type, dim, dim, Eigen::RowMajor> > LHS_map(LHS.data());
        Eigen::Map<Eigen::Matrix<deformationGradient_type, dim, dim, Eigen::RowMajor> > inv_LHS_map(inv_LHS.data());

        inv_LHS_map = LHS_map.inverse().eval();

        std::fill(dFdL_begin, dFdL_end, deformationGradient_type());

        // Compute the jacobian
        if (mode == 1) {
            for (unsigned int j = 0; j < dim; ++j) {
                for (unsigned int I = 0; I < dim; ++I) {
                    for (unsigned int k = 0; k < dim; ++k) {
                        for (unsigned int l = 0; l < dim; ++l) {
                            *(dFdL_begin + dim * dim * dim * j + dim * dim * I + dim * k + l) +=
                                inv_LHS[dim * j + k] * (*(deformationGradient_begin + dim * l + I));
                        }
                    }
                }
            }
        } else if (mode == 2) {
            for (unsigned int j = 0; j < dim; ++j) {
                for (unsigned int I = 0; I < dim; ++I) {
                    for (unsigned int K = 0; K < dim; ++K) {
                        for (unsigned int _L = 0; _L < dim; ++_L) {
                            *(dFdL_begin + dim * dim * dim * j + dim * dim * I + dim * K + _L) +=
                                inv_LHS[dim * _L + I] * (*(deformationGradient_begin + dim * j + K));
                        }
                    }
                }
            }
        }

        std::transform(dFdL_begin, dFdL_end, dFdL_begin,
                       std::bind(std::multiplies<>(), std::placeholders::_1, Dt * (1 - alpha)));

        return;
    }

    template <unsigned int dim, unsigned int mode, typename Dt_type, class previousDeformationGradient_iterator,
              class Lp_iterator, class L_iterator, class dF_iterator, class deformationGradient_iterator,
              class dFdL_iterator, class ddFdFp_iterator, class dFdFp_iterator, class dFdLp_iterator,
              typename alpha_type>
    void evolveF(const Dt_type &Dt, const previousDeformationGradient_iterator &previousDeformationGradient_begin,
                 const previousDeformationGradient_iterator &previousDeformationGradient_end,
                 const Lp_iterator &Lp_begin, const Lp_iterator &Lp_end, const L_iterator &L_begin,
                 const L_iterator &L_end, dF_iterator dF_begin, dF_iterator dF_end,
                 deformationGradient_iterator deformationGradient_begin,
                 deformationGradient_iterator deformationGradient_end, dFdL_iterator dFdL_begin, dFdL_iterator dFdL_end,
                 dFdFp_iterator ddFdFp_begin, ddFdFp_iterator ddFdFp_end, dFdFp_iterator dFdFp_begin,
                 dFdFp_iterator dFdFp_end, dFdLp_iterator dFdLp_begin, dFdLp_iterator dFdLp_end,
                 const alpha_type alpha) {
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
         * \param &previousDeformationGradient_begin: The starting iterator of the previous value of the deformation
         * gradient \param &previousDeformationGradient_end: The stopping iterator of the previous value of the
         * deformation gradient \param &Lp_begin: The starting iterator previous velocity gradient in the current
         * configuration (mode 1) or reference configuration (mode 2). \param &Lp_end: The stopping iterator previous
         * velocity gradient in the current configuration (mode 1) or reference configuration (mode 2). \param &L_begin:
         * The starting iterator of the current velocity gradient in the current configuration (mode 1) or reference
         * configuration (mode 2). \param &L_end: The stopping iterator of the current velocity gradient in the current
         * configuration (mode 1) or reference configuration (mode 2). \param &dF_begin: The starting iterator of the
         * change in the deformation gradient \f$\Delta \bf{F}\f$ such that \f$F_{iI}^{t+1} = F_{iI}^t + \Delta
         * F_{iI}\f$ \param &dF_end: The stopping iterator of the change in the deformation gradient \f$\Delta \bf{F}\f$
         * such that \f$F_{iI}^{t+1} = F_{iI}^t + \Delta F_{iI}\f$ \param &deformationGradient_begin: The starting
         * iterator of the computed current deformation gradient. \param &deformationGradient_end: The starting iterator
         * of the computed current deformation gradient. \param &dFdL_begin: The starting iterator of the derivative of
         * the deformation gradient w.r.t. the velocity gradient \param &dFdL_end: The stopping iterator of the
         * derivative of the deformation gradient w.r.t. the velocity gradient \param &ddFdFp_begin: The starting
         * iterator of the derivative of the change in the deformation gradient w.r.t. the previous deformation gradient
         * \param &ddFdFp_end: The stopping iterator of the derivative of the change in the deformation gradient w.r.t.
         * the previous deformation gradient \param &dFdFp_begin: The starting iterator of the derivative of the
         * deformation gradient w.r.t. the previous deformation gradient \param &dFdFp_end: The stopping iterator of the
         * derivative of the deformation gradient w.r.t. the previous deformation gradient \param &dFdLp_begin: The
         * starting iterator of the derivative of the deformation gradient w.r.t. the previous velocity gradient \param
         * &dFdLp_end: The stopping iterator of the derivative of the deformation gradient w.r.t. the previous velocity
         * gradient \param alpha: The integration parameter.
         */

        using deformationGradient_type = typename std::iterator_traits<deformationGradient_iterator>::value_type;

        TARDIGRADE_ERROR_TOOLS_CHECK(
            (unsigned int)(dFdL_end - dFdL_begin) == dim * dim * dim * dim,
            "The derivative of the deformation gradient w.r.t. the velocity gradient has a size of " +
                std::to_string((unsigned int)(dFdL_end - dFdL_begin)) + " but must have a size of " +
                std::to_string(dim * dim * dim * dim));

        TARDIGRADE_ERROR_TOOLS_CHECK((unsigned int)(ddFdFp_end - ddFdFp_begin) == dim * dim * dim * dim,
                                     "The derivative of the change in the deformation gradient w.r.t. the previous "
                                     "deformation gradient has a size of " +
                                         std::to_string((unsigned int)(ddFdFp_end - ddFdFp_begin)) +
                                         " but must have a size of " + std::to_string(dim * dim * dim * dim));

        TARDIGRADE_ERROR_TOOLS_CHECK(
            (unsigned int)(dFdFp_end - dFdFp_begin) == dim * dim * dim * dim,
            "The derivative of the deformation gradient w.r.t. the previous deformation gradient has a size of " +
                std::to_string((unsigned int)(dFdFp_end - dFdFp_begin)) + " but must have a size of " +
                std::to_string(dim * dim * dim * dim));

        TARDIGRADE_ERROR_TOOLS_CHECK(
            (unsigned int)(dFdLp_end - dFdLp_begin) == dim * dim * dim * dim,
            "The derivative of the deformation gradient w.r.t. the previous velocity gradient has a size of " +
                std::to_string((unsigned int)(dFdLp_end - dFdLp_begin)) + " but must have a size of " +
                std::to_string(dim * dim * dim * dim));

        evolveF<dim, mode>(Dt, previousDeformationGradient_begin, previousDeformationGradient_end, Lp_begin, Lp_end,
                           L_begin, L_end, dF_begin, dF_end, deformationGradient_begin, deformationGradient_end, alpha);

        // Construct the velocity gradient at time t + alpha
        std::array<deformationGradient_type, dim * dim> LtpAlpha;
        for (unsigned int i = 0; i < dim * dim; ++i) {
            LtpAlpha[i] = alpha * (*(Lp_begin + i)) + (1 - alpha) * (*(L_begin + i));
        }

        // Compute the left-hand side
        std::array<deformationGradient_type, dim * dim> LHS, inv_LHS;
        std::transform(L_begin, L_end, std::begin(LHS), std::begin(LHS),
                       std::bind(std::multiplies<>(), std::placeholders::_1, -Dt * (1 - alpha)));

        for (unsigned int i = 0; i < dim; i++) {
            LHS[dim * i + i] += 1;
        }

        Eigen::Map<Eigen::Matrix<deformationGradient_type, dim, dim, Eigen::RowMajor> > LHS_map(LHS.data());
        Eigen::Map<Eigen::Matrix<deformationGradient_type, dim, dim, Eigen::RowMajor> > inv_LHS_map(inv_LHS.data());

        inv_LHS_map = LHS_map.inverse().eval();

        std::fill(dFdL_begin, dFdL_end, deformationGradient_type());
        std::fill(ddFdFp_begin, ddFdFp_end, deformationGradient_type());
        std::fill(dFdLp_begin, dFdLp_end, deformationGradient_type());

        if (mode == 1) {
            for (unsigned int j = 0; j < dim; ++j) {
                for (unsigned int I = 0; I < dim; ++I) {
                    for (unsigned int k = 0; k < dim; ++k) {
                        for (unsigned int l = 0; l < dim; ++l) {
                            *(dFdL_begin + dim * dim * dim * j + dim * dim * I + dim * k + l) +=
                                inv_LHS[dim * j + k] * (*(deformationGradient_begin + dim * l + I));
                            *(dFdLp_begin + dim * dim * dim * j + dim * dim * I + dim * k + l) +=
                                inv_LHS[dim * j + k] * (*(previousDeformationGradient_begin + dim * l + I));
                            *(ddFdFp_begin + dim * dim * dim * j + dim * dim * I + dim * k + I) +=
                                inv_LHS[dim * j + l] * LtpAlpha[dim * l + k];
                        }
                    }
                }
            }
        }
        if (mode == 2) {
            for (unsigned int j = 0; j < dim; j++) {
                for (unsigned int I = 0; I < dim; I++) {
                    for (unsigned int K = 0; K < dim; K++) {
                        for (unsigned int _L = 0; _L < dim; _L++) {
                            *(dFdL_begin + dim * dim * dim * j + dim * dim * I + dim * K + _L) +=
                                inv_LHS[dim * _L + I] * (*(deformationGradient_begin + dim * j + K));
                            *(dFdLp_begin + dim * dim * dim * j + dim * dim * I + dim * K + _L) +=
                                inv_LHS[dim * _L + I] * (*(previousDeformationGradient_begin + dim * j + K));
                            *(ddFdFp_begin + dim * dim * dim * j + dim * dim * I + dim * j + K) +=
                                LtpAlpha[dim * K + _L] * inv_LHS[dim * _L + I];
                        }
                    }
                }
            }
        }

        std::transform(dFdL_begin, dFdL_end, dFdL_begin,
                       std::bind(std::multiplies<>(), std::placeholders::_1, Dt * (1 - alpha)));

        std::transform(dFdLp_begin, dFdLp_end, dFdLp_begin,
                       std::bind(std::multiplies<>(), std::placeholders::_1, Dt * alpha));

        std::transform(ddFdFp_begin, ddFdFp_end, ddFdFp_begin,
                       std::bind(std::multiplies<>(), std::placeholders::_1, Dt));

        std::copy(ddFdFp_begin, ddFdFp_end, dFdFp_begin);
        for (unsigned int i = 0; i < dim * dim; ++i) {
            *(dFdFp_begin + dim * dim * i + i) += 1;
        }
    }

    template <class A_iterator, class Anorm_iterator>
    void computeUnitNormal(const A_iterator &A_begin, const A_iterator &A_end, Anorm_iterator Anorm_begin,
                           Anorm_iterator Anorm_end) {
        /*!
         * Compute the unit normal of a second order tensor (or strictly speaking
         * any tensor).
         *
         * \param &A_begin: The starting iterator for the incoming vector
         * \param &A_end: The stopping iterator for the incoming vector
         * \param &Anorm_begin: The starting iterator of the unit normal in the direction of A
         * \param &Anorm_end: The stopping iterator of the unit normal in the direction of A
         */

        using A_type     = typename std::iterator_traits<A_iterator>::value_type;
        using Anorm_type = typename std::iterator_traits<Anorm_iterator>::value_type;

        TARDIGRADE_ERROR_TOOLS_CHECK((unsigned int)(A_end - A_begin) == (unsigned int)(Anorm_end - Anorm_begin),
                                     "A has a size of " + std::to_string((unsigned int)(A_end - A_begin)) +
                                         " but Anorm has a size of " +
                                         std::to_string((unsigned int)(Anorm_end - Anorm_begin)));

        auto norm = tardigradeVectorTools::l2norm<A_type>(A_begin, A_end);

        if (tardigradeVectorTools::fuzzyEquals(norm, 0.)) {
            std::fill(Anorm_begin, Anorm_end, Anorm_type());
        } else {
            std::transform(A_begin, A_end, Anorm_begin, std::bind(std::divides<>(), std::placeholders::_1, norm));
        }
    }

    template <class A_iterator, class Anorm_iterator, class dAnormdA_iterator>
    void computeUnitNormal(const A_iterator &A_begin, const A_iterator &A_end, Anorm_iterator Anorm_begin,
                           Anorm_iterator Anorm_end, dAnormdA_iterator dAnormdA_begin, dAnormdA_iterator dAnormdA_end) {
        /*!
         * Compute the unit normal of a second order tensor (or strictly speaking
         * any tensor).
         *
         * \param &A_begin: The starting iterator for the incoming vector
         * \param &A_end: The stopping iterator for the incoming vector
         * \param &Anorm_begin: The starting iterator of the unit normal in the direction of A
         * \param &Anorm_end: The stopping iterator of the unit normal in the direction of A
         * \param &dAnormdA_begin: The starting iterator of the gradient of the unit normal in the direction of A w.r.t.
         * A \param &dAnormdA_end: The stopping iterator of the gradient of the unit normal in the direction of A w.r.t.
         * A
         */

        using A_type     = typename std::iterator_traits<A_iterator>::value_type;
        using Anorm_type = typename std::iterator_traits<Anorm_iterator>::value_type;

        auto size = (unsigned int)(A_end - A_begin);

        TARDIGRADE_ERROR_TOOLS_CHECK(size == (unsigned int)(Anorm_end - Anorm_begin),
                                     "A has a size of " + std::to_string(size) + " but Anorm has a size of " +
                                         std::to_string((unsigned int)(Anorm_end - Anorm_begin)));

        TARDIGRADE_ERROR_TOOLS_CHECK(size * size == (unsigned int)(dAnormdA_end - dAnormdA_begin),
                                     "dAnormdA should have a size of " + std::to_string(size * size) +
                                         " but it has a size of " +
                                         std::to_string((unsigned int)(dAnormdA_end - dAnormdA_begin)));

        auto norm = tardigradeVectorTools::l2norm<A_type>(A_begin, A_end);

        if (tardigradeVectorTools::fuzzyEquals(norm, 0.)) {
            std::fill(Anorm_begin, Anorm_end, Anorm_type());
        } else {
            std::transform(A_begin, A_end, Anorm_begin, std::bind(std::divides<>(), std::placeholders::_1, norm));
        }

        std::fill(dAnormdA_begin, dAnormdA_end, Anorm_type());

        for (auto v = std::pair<unsigned int, Anorm_iterator>(0, Anorm_begin); v.second != Anorm_end;
             ++v.first, ++v.second) {
            *(dAnormdA_begin + size * v.first + v.first) += 1;

            for (auto w = std::pair<unsigned int, Anorm_iterator>(0, Anorm_begin); w.second != Anorm_end;
                 ++w.first, ++w.second) {
                *(dAnormdA_begin + size * v.first + w.first) -= (*v.second) * (*w.second);
            }
        }

        std::transform(dAnormdA_begin, dAnormdA_end, dAnormdA_begin,
                       std::bind(std::divides<>(), std::placeholders::_1, norm));
    }

    template <unsigned int dim, class velocityGradient_iterator, class deformationGradient_iterator,
              class pulledBackVelocityGradient_iterator>
    void pullBackVelocityGradient(const velocityGradient_iterator    &velocityGradient_begin,
                                  const velocityGradient_iterator    &velocityGradient_end,
                                  const deformationGradient_iterator &deformationGradient_begin,
                                  const deformationGradient_iterator &deformationGradient_end,
                                  pulledBackVelocityGradient_iterator pulledBackVelocityGradient_begin,
                                  pulledBackVelocityGradient_iterator pulledBackVelocityGradient_end) {
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
         * \param &velocityGradient_begin: The starting iterator of the velocity gradient in the current configuration.
         * \param &velocityGradient_end: The stopping iterator of the velocity gradient in the current configuration.
         * \param &deformationGradient_begin: The starting iterator of the deformation gradient between the desired
         * configuration and the current configuration. \param &deformationGradient_end: The stopping iterator of the
         * deformation gradient between the desired configuration and the current configuration. \param
         * &pulledBackVelocityGradient_begin: The starting iterator of the pulled back velocity gradient. \param
         * &pulledBackVelocityGradient_end: The stopping iterator of the pulled back velocity gradient.
         */

        using velocityGradient_type    = typename std::iterator_traits<velocityGradient_iterator>::value_type;
        using deformationGradient_type = typename std::iterator_traits<deformationGradient_iterator>::value_type;
        using pulledBackVelocityGradient_type =
            typename std::iterator_traits<pulledBackVelocityGradient_iterator>::value_type;

        TARDIGRADE_ERROR_TOOLS_CHECK(dim * dim == (unsigned int)(velocityGradient_end - velocityGradient_begin),
                                     "The velocity gradient has a size of " +
                                         std::to_string((unsigned int)(velocityGradient_end - velocityGradient_begin)) +
                                         " but it shold be " + std::to_string(dim * dim));

        TARDIGRADE_ERROR_TOOLS_CHECK(dim * dim == (unsigned int)(deformationGradient_end - deformationGradient_begin),
                                     "The deformation gradient has a size of " +
                                         std::to_string((unsigned int)(deformationGradient_end -
                                                                       deformationGradient_begin)) +
                                         " but it shold be " + std::to_string(dim * dim));

        TARDIGRADE_ERROR_TOOLS_CHECK(
            dim * dim == (unsigned int)(pulledBackVelocityGradient_end - pulledBackVelocityGradient_begin),
            "The pulled back velocity gradient has a size of " +
                std::to_string((unsigned int)(pulledBackVelocityGradient_end - pulledBackVelocityGradient_begin)) +
                " but it shold be " + std::to_string(dim * dim));

        Eigen::Map<const Eigen::Matrix<deformationGradient_type, dim, dim, Eigen::RowMajor> > F(
            &(*deformationGradient_begin));
        Eigen::Map<const Eigen::Matrix<velocityGradient_type, dim, dim, Eigen::RowMajor> > L(
            &(*velocityGradient_begin));
        Eigen::Map<Eigen::Matrix<pulledBackVelocityGradient_type, dim, dim, Eigen::RowMajor> > pullBackL(
            &(*pulledBackVelocityGradient_begin));

        pullBackL = (F.inverse() * L * F).eval();

        return;
    }

    template <unsigned int dim, class velocityGradient_iterator, class deformationGradient_iterator,
              class pulledBackVelocityGradient_iterator, class dPullBackLdL_iterator, class dPullBackLdF_iterator>
    void pullBackVelocityGradient(const velocityGradient_iterator    &velocityGradient_begin,
                                  const velocityGradient_iterator    &velocityGradient_end,
                                  const deformationGradient_iterator &deformationGradient_begin,
                                  const deformationGradient_iterator &deformationGradient_end,
                                  pulledBackVelocityGradient_iterator pulledBackVelocityGradient_begin,
                                  pulledBackVelocityGradient_iterator pulledBackVelocityGradient_end,
                                  dPullBackLdL_iterator dPullBackLdL_begin, dPullBackLdL_iterator dPullBackLdL_end,
                                  dPullBackLdF_iterator dPullBackLdF_begin, dPullBackLdF_iterator dPullBackLdF_end) {
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
         * \param &velocityGradient_begin: The starting iterator of the velocity gradient in the current configuration.
         * \param &velocityGradient_end: The stopping iterator of the velocity gradient in the current configuration.
         * \param &deformationGradient_begin: The starting iterator of the deformation gradient between the desired
         * configuration and the current configuration. \param &deformationGradient_end: The stopping iterator of the
         * deformation gradient between the desired configuration and the current configuration. \param
         * &pulledBackVelocityGradient_begin: The starting iterator of the pulled back velocity gradient. \param
         * &pulledBackVelocityGradient_end: The stopping iterator of the pulled back velocity gradient. \param
         * &dPullBackLdL_begin: The starting iterator of the gradient of the pulled back velocity gradient w.r.t. the
         * velocity gradient. \param &dPullBackLdL_end: The stopping iterator of the gradient of the pulled back
         * velocity gradient w.r.t. the velocity gradient. \param &dPullBackLdF_begin: The starting iterator of the
         * gradient of the pulled back velocity gradient w.r.t. the deformation gradient. \param &dPullBackLdF_end: The
         * stopping iterator of the gradient of the pulled back velocity gradient w.r.t. the deformation gradient.
         */

        using velocityGradient_type    = typename std::iterator_traits<velocityGradient_iterator>::value_type;
        using deformationGradient_type = typename std::iterator_traits<deformationGradient_iterator>::value_type;
        using pulledBackVelocityGradient_type =
            typename std::iterator_traits<pulledBackVelocityGradient_iterator>::value_type;
        using dPullBackLdL_type = typename std::iterator_traits<dPullBackLdL_iterator>::value_type;
        using dPullBackLdF_type = typename std::iterator_traits<dPullBackLdF_iterator>::value_type;

        TARDIGRADE_ERROR_TOOLS_CHECK(dim * dim == (unsigned int)(velocityGradient_end - velocityGradient_begin),
                                     "The velocity gradient has a size of " +
                                         std::to_string((unsigned int)(velocityGradient_end - velocityGradient_begin)) +
                                         " but it shold be " + std::to_string(dim * dim));

        TARDIGRADE_ERROR_TOOLS_CHECK(dim * dim == (unsigned int)(deformationGradient_end - deformationGradient_begin),
                                     "The deformation gradient has a size of " +
                                         std::to_string((unsigned int)(deformationGradient_end -
                                                                       deformationGradient_begin)) +
                                         " but it shold be " + std::to_string(dim * dim));

        TARDIGRADE_ERROR_TOOLS_CHECK(
            dim * dim == (unsigned int)(pulledBackVelocityGradient_end - pulledBackVelocityGradient_begin),
            "The pulled back velocity gradient has a size of " +
                std::to_string((unsigned int)(pulledBackVelocityGradient_end - pulledBackVelocityGradient_begin)) +
                " but it shold be " + std::to_string(dim * dim));

        TARDIGRADE_ERROR_TOOLS_CHECK(
            dim * dim * dim * dim == (unsigned int)(dPullBackLdL_end - dPullBackLdL_begin),
            "The derivative of the pulled back velocity gradient with respect to the velocity gradient has a size of " +
                std::to_string((unsigned int)(dPullBackLdL_end - dPullBackLdL_begin)) + " but it shold be " +
                std::to_string(dim * dim * dim * dim));

        TARDIGRADE_ERROR_TOOLS_CHECK(dim * dim * dim * dim == (unsigned int)(dPullBackLdF_end - dPullBackLdF_begin),
                                     "The derivative of the pulled back velocity gradient with respect to the "
                                     "deformation gradient has a size of " +
                                         std::to_string((unsigned int)(dPullBackLdF_end - dPullBackLdF_begin)) +
                                         " but it shold be " + std::to_string(dim * dim * dim * dim));

        Eigen::Map<const Eigen::Matrix<deformationGradient_type, dim, dim, Eigen::RowMajor> > F(
            &(*deformationGradient_begin));
        Eigen::Map<const Eigen::Matrix<velocityGradient_type, dim, dim, Eigen::RowMajor> > L(
            &(*velocityGradient_begin));
        Eigen::Map<Eigen::Matrix<pulledBackVelocityGradient_type, dim, dim, Eigen::RowMajor> > pullBackL(
            &(*pulledBackVelocityGradient_begin));

        std::array<deformationGradient_type, dim * dim>        inverseDeformationGradient;
        std::array<pulledBackVelocityGradient_type, dim * dim> term2;

        Eigen::Map<Eigen::Matrix<deformationGradient_type, dim, dim, Eigen::RowMajor> > invF_map(
            inverseDeformationGradient.data());
        invF_map = F.inverse().eval();

        Eigen::Map<Eigen::Matrix<pulledBackVelocityGradient_type, dim, dim, Eigen::RowMajor> > term2_map(term2.data());
        term2_map = (invF_map * L).eval();

        // Pull back the velocity gradient
        pullBackL = (term2_map * F).eval();

        // Construct the gradients
        std::fill(dPullBackLdL_begin, dPullBackLdL_end, dPullBackLdL_type());
        std::fill(dPullBackLdF_begin, dPullBackLdF_end, dPullBackLdF_type());

        for (unsigned int I = 0; I < dim; ++I) {
            for (unsigned int J = 0; J < dim; ++J) {
                for (unsigned int k = 0; k < dim; ++k) {
                    for (unsigned int l = 0; l < dim; ++l) {
                        *(dPullBackLdL_begin + dim * dim * dim * I + dim * dim * J + dim * k + l) =
                            inverseDeformationGradient[dim * I + k] * (*(deformationGradient_begin + dim * l + J));
                    }

                    *(dPullBackLdF_begin + dim * dim * dim * I + dim * dim * J + dim * k + J) += term2[dim * I + k];

                    for (unsigned int K = 0; K < dim; ++K) {
                        *(dPullBackLdF_begin + dim * dim * dim * I + dim * dim * J + dim * k + K) -=
                            inverseDeformationGradient[dim * I + k] *
                            (*(pulledBackVelocityGradient_begin + dim * K + J));
                    }
                }
            }
        }

        return;
    }

    template <typename temperature_type, typename referenceTemperature_type, class linearParameters_iterator,
              class quadraticParameters_iterator, class thermalExpansion_iterator>
    void quadraticThermalExpansion(const temperature_type             &temperature,
                                   const referenceTemperature_type    &referenceTemperature,
                                   const linearParameters_iterator    &linearParameters_begin,
                                   const linearParameters_iterator    &linearParameters_end,
                                   const quadraticParameters_iterator &quadraticParameters_begin,
                                   const quadraticParameters_iterator &quadraticParameters_end,
                                   thermalExpansion_iterator           thermalExpansion_begin,
                                   thermalExpansion_iterator           thermalExpansion_end) {
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
         * \param &linearParameters_begin: The starting iterator of the linear thermal expansion parameters.
         * \param &linearParameters_end: The stopping iterator of the linear thermal expansion parameters.
         * \param &quadraticParameters_begin: The starting iterator of the quadratic thermal expansion parameters.
         * \param &quadraticParameters_end: The stopping iterator of the quadratic thermal expansion parameters.
         * \param &thermalExpansion_begin: The starting iterator of the resulting thermal expansion.
         * \param &thermalExpansion_end: The stopping iterator of the resulting thermal expansion.
         */

        TARDIGRADE_ERROR_TOOLS_EVAL(const unsigned int parameters_dim =
                                        (unsigned int)(linearParameters_end - linearParameters_begin);)

        TARDIGRADE_ERROR_TOOLS_CHECK(
            parameters_dim == (unsigned int)(quadraticParameters_end - quadraticParameters_begin),
            "The quadratic parameters array has a size of " +
                std::to_string((unsigned int)(quadraticParameters_end - quadraticParameters_begin)) +
                " and should have a size of " + std::to_string(parameters_dim))

        TARDIGRADE_ERROR_TOOLS_CHECK(parameters_dim == (unsigned int)(thermalExpansion_end - thermalExpansion_begin),
                                     "The thermal expansion array has a size of " +
                                         std::to_string((unsigned int)(thermalExpansion_end - thermalExpansion_begin)) +
                                         " and should have a size of " + std::to_string(parameters_dim))

        for (auto v = std::pair<unsigned int, thermalExpansion_iterator>(0, thermalExpansion_begin);
             v.second != thermalExpansion_end; ++v.first, ++v.second) {
            *v.second = (*(linearParameters_begin + v.first)) * temperature +
                        (*(quadraticParameters_begin + v.first)) * temperature * temperature -
                        (*(linearParameters_begin + v.first)) * referenceTemperature -
                        (*(quadraticParameters_begin + v.first)) * referenceTemperature * referenceTemperature;
        }
    }

    template <typename temperature_type, typename referenceTemperature_type, class linearParameters_iterator,
              class quadraticParameters_iterator, class thermalExpansion_iterator,
              class thermalExpansionJacobian_iterator>
    void quadraticThermalExpansion(const temperature_type             &temperature,
                                   const referenceTemperature_type    &referenceTemperature,
                                   const linearParameters_iterator    &linearParameters_begin,
                                   const linearParameters_iterator    &linearParameters_end,
                                   const quadraticParameters_iterator &quadraticParameters_begin,
                                   const quadraticParameters_iterator &quadraticParameters_end,
                                   thermalExpansion_iterator           thermalExpansion_begin,
                                   thermalExpansion_iterator           thermalExpansion_end,
                                   thermalExpansionJacobian_iterator   thermalExpansionJacobian_begin,
                                   thermalExpansionJacobian_iterator   thermalExpansionJacobian_end) {
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
         * \param &linearParameters_begin: The starting iterator of the linear thermal expansion parameters.
         * \param &linearParameters_end: The stopping iterator of the linear thermal expansion parameters.
         * \param &quadraticParameters_begin: The starting iterator of the quadratic thermal expansion parameters.
         * \param &quadraticParameters_end: The stopping iterator of the quadratic thermal expansion parameters.
         * \param thermalExpansion_begin: The starting iterator of the resulting thermal expansion.
         * \param thermalExpansion_end: The stopping iterator of the resulting thermal expansion.
         * \param thermalExpansionJacobian_begin: The starting iterator of the Jacobian of the thermal expansion w.r.t.
         * temperature \param thermalExpansionJacobian_end: The stopping iterator of the Jacobian of the thermal
         * expansion w.r.t. temperature
         */

        TARDIGRADE_ERROR_TOOLS_EVAL(const unsigned int parameters_dim =
                                        (unsigned int)(linearParameters_end - linearParameters_begin);)

        TARDIGRADE_ERROR_TOOLS_CHECK(
            parameters_dim == (unsigned int)(quadraticParameters_end - quadraticParameters_begin),
            "The quadratic parameters array has a size of " +
                std::to_string((unsigned int)(quadraticParameters_end - quadraticParameters_begin)) +
                " and should have a size of " + std::to_string(parameters_dim))

        TARDIGRADE_ERROR_TOOLS_CHECK(parameters_dim == (unsigned int)(thermalExpansion_end - thermalExpansion_begin),
                                     "The thermal expansion array has a size of " +
                                         std::to_string((unsigned int)(thermalExpansion_end - thermalExpansion_begin)) +
                                         " and should have a size of " + std::to_string(parameters_dim))

        TARDIGRADE_ERROR_TOOLS_CHECK(parameters_dim ==
                                         (unsigned int)(thermalExpansionJacobian_end - thermalExpansionJacobian_begin),
                                     "The thermal expansion Jacobian array has a size of " +
                                         std::to_string((unsigned int)(thermalExpansion_end - thermalExpansion_begin)) +
                                         " and should have a size of " + std::to_string(parameters_dim))

        for (auto v = std::pair<unsigned int, thermalExpansion_iterator>(0, thermalExpansion_begin);
             v.second != thermalExpansion_end; ++v.first, ++v.second) {
            *v.second = (*(linearParameters_begin + v.first)) * temperature +
                        (*(quadraticParameters_begin + v.first)) * temperature * temperature -
                        (*(linearParameters_begin + v.first)) * referenceTemperature -
                        (*(quadraticParameters_begin + v.first)) * referenceTemperature * referenceTemperature;

            *(thermalExpansionJacobian_begin + v.first) =
                (*(linearParameters_begin + v.first)) + 2 * (*(quadraticParameters_begin + v.first)) * temperature;
        }
    }

    template <unsigned int dim, class greenLagrangeStrain_iterator, class deformationGradient_iterator,
              class almansiStrain_iterator>
    void pushForwardGreenLagrangeStrain(const greenLagrangeStrain_iterator &greenLagrangeStrain_begin,
                                        const greenLagrangeStrain_iterator &greenLagrangeStrain_end,
                                        const deformationGradient_iterator &deformationGradient_begin,
                                        const deformationGradient_iterator &deformationGradient_end,
                                        almansiStrain_iterator              almansiStrain_begin,
                                        almansiStrain_iterator              almansiStrain_end) {
        /*!
         * Push forward the Green-Lagrange strain to the current configuration.
         *
         * \f$e_{ij} = F_{Ii}^{-1} E_{IJ} F_{Jj}^{-1}\f$
         *
         * where \f$e_{ij}\f$ is the Almansi strain (the strain in the current configuration, \f$F_{iI}^{-1}\f$ is the
         * inverse of the deformation gradient, and \f$E_{IJ}\f$ is the Green-Lagrange strain.
         *
         * \param &greenLagrangeStrain_begin: The starting iterator of the Green-Lagrange strain.
         * \param &greenLagrangeStrain_end: The stopping iterator of the Green-Lagrange strain.
         * \param &deformationGradient_begin: The starting iterator of the deformation gradient mapping between
         * configurations. \param &deformationGradient_end: The stopping iterator of the deformation gradient mapping
         * between configurations. \param &almansiStrain_begin: The starting iterator of the strain in the current
         * configuration indicated by the deformation gradient. \param &almansiStrain_end: The stopping iterator of the
         * strain in the current configuration indicated by the deformation gradient.
         */

        using greenLagrangeStrain_type = typename std::iterator_traits<greenLagrangeStrain_iterator>::value_type;
        using deformationGradient_type = typename std::iterator_traits<deformationGradient_iterator>::value_type;
        using almansiStrain_type       = typename std::iterator_traits<almansiStrain_iterator>::value_type;

        TARDIGRADE_ERROR_TOOLS_CHECK(dim * dim == (unsigned int)(greenLagrangeStrain_end - greenLagrangeStrain_begin),
                                     "The Green-Lagrange strain has a size of " +
                                         std::to_string((unsigned int)(greenLagrangeStrain_end -
                                                                       greenLagrangeStrain_begin)) +
                                         " but should have a size of " + std::to_string(dim * dim))

        TARDIGRADE_ERROR_TOOLS_CHECK(dim * dim == (unsigned int)(deformationGradient_end - deformationGradient_begin),
                                     "The deformation gradient has a size of " +
                                         std::to_string((unsigned int)(deformationGradient_end -
                                                                       deformationGradient_begin)) +
                                         " but should have a size of " + std::to_string(dim * dim))

        TARDIGRADE_ERROR_TOOLS_CHECK(dim * dim == (unsigned int)(almansiStrain_end - almansiStrain_begin),
                                     "The Almansi strain has a size of " +
                                         std::to_string((unsigned int)(almansiStrain_end - almansiStrain_begin)) +
                                         " but should have a size of " + std::to_string(dim * dim))

        std::array<deformationGradient_type, dim * dim> inverseDeformationGradient;
        Eigen::Map<const Eigen::Matrix<deformationGradient_type, dim, dim, Eigen::RowMajor> > F(
            &(*deformationGradient_begin));
        Eigen::Map<Eigen::Matrix<deformationGradient_type, dim, dim, Eigen::RowMajor> > inv_F(
            inverseDeformationGradient.data());

        inv_F = F.inverse().eval();

        Eigen::Map<const Eigen::Matrix<greenLagrangeStrain_type, dim, dim, Eigen::RowMajor> > E(
            &(*greenLagrangeStrain_begin));
        Eigen::Map<Eigen::Matrix<almansiStrain_type, dim, dim, Eigen::RowMajor> > e(&(*almansiStrain_begin));

        // Map the Green-Lagrange strain to the current configuration
        e = (inv_F.transpose() * E * inv_F).eval();

        return;
    }

    template <unsigned int dim, class greenLagrangeStrain_iterator, class deformationGradient_iterator,
              class almansiStrain_iterator, class dAlmansiStraindE_iterator, class dAlmansiStraindF_iterator>
    void pushForwardGreenLagrangeStrain(const greenLagrangeStrain_iterator &greenLagrangeStrain_begin,
                                        const greenLagrangeStrain_iterator &greenLagrangeStrain_end,
                                        const deformationGradient_iterator &deformationGradient_begin,
                                        const deformationGradient_iterator &deformationGradient_end,
                                        almansiStrain_iterator              almansiStrain_begin,
                                        almansiStrain_iterator              almansiStrain_end,
                                        dAlmansiStraindE_iterator           dAlmansiStraindE_begin,
                                        dAlmansiStraindE_iterator           dAlmansiStraindE_end,
                                        dAlmansiStraindF_iterator           dAlmansiStraindF_begin,
                                        dAlmansiStraindF_iterator           dAlmansiStraindF_end) {
        /*!
         * Push forward the Green-Lagrange strain to the current configuration.
         *
         * \f$e_{ij} = F_{Ii}^{-1} E_{IJ} F_{Jj}^{-1}\f$
         *
         * where \f$e_{ij}\f$ is the Almansi strain (the strain in the current configuration, \f$F_{iI}^{-1}\f$ is the
         * inverse of the deformation gradient, and \f$E_{IJ}\f$ is the Green-Lagrange strain.
         *
         * \param &greenLagrangeStrain_begin: The starting iterator of the Green-Lagrange strain.
         * \param &greenLagrangeStrain_end: The stopping iterator of the Green-Lagrange strain.
         * \param &deformationGradient_begin: The starting iterator of the deformation gradient mapping between
         * configurations. \param &deformationGradient_end: The stopping iterator of the deformation gradient mapping
         * between configurations. \param &almansiStrain_begin: The starting iterator of the strain in the current
         * configuration indicated by the deformation gradient. \param &almansiStrain_end: The stopping iterator of the
         * strain in the current configuration indicated by the deformation gradient. \param &dAlmansiStraindE_begin:
         * The starting iterator for the Jacobian of the Almansi strain w.r.t. the Green-Lagrange strain. \param
         * &dAlmansiStraindE_end: The stopping iterator for the Jacobian of the Almansi strain w.r.t. the Green-Lagrange
         * strain. \param &dAlmansiStraindF_begin: The starting iterator for the Jacobian of the Almansi strain w.r.t.
         * the deformation gradient. \param &dAlmansiStraindF_end: The stopping iterator for the Jacobian of the Almansi
         * strain w.r.t. the deformation gradient.
         */

        using greenLagrangeStrain_type = typename std::iterator_traits<greenLagrangeStrain_iterator>::value_type;
        using deformationGradient_type = typename std::iterator_traits<deformationGradient_iterator>::value_type;
        using almansiStrain_type       = typename std::iterator_traits<almansiStrain_iterator>::value_type;

        TARDIGRADE_ERROR_TOOLS_CHECK(dim * dim == (unsigned int)(greenLagrangeStrain_end - greenLagrangeStrain_begin),
                                     "The Green-Lagrange strain has a size of " +
                                         std::to_string((unsigned int)(greenLagrangeStrain_end -
                                                                       greenLagrangeStrain_begin)) +
                                         " but should have a size of " + std::to_string(dim * dim))

        TARDIGRADE_ERROR_TOOLS_CHECK(dim * dim == (unsigned int)(deformationGradient_end - deformationGradient_begin),
                                     "The deformation gradient has a size of " +
                                         std::to_string((unsigned int)(deformationGradient_end -
                                                                       deformationGradient_begin)) +
                                         " but should have a size of " + std::to_string(dim * dim))

        TARDIGRADE_ERROR_TOOLS_CHECK(dim * dim == (unsigned int)(almansiStrain_end - almansiStrain_begin),
                                     "The Almansi strain has a size of " +
                                         std::to_string((unsigned int)(almansiStrain_end - almansiStrain_begin)) +
                                         " but should have a size of " + std::to_string(dim * dim))

        TARDIGRADE_ERROR_TOOLS_CHECK(
            dim * dim * dim * dim == (unsigned int)(dAlmansiStraindE_end - dAlmansiStraindE_begin),
            "The derivative of the Almansi strain with respect to the Green-Lagrange strain has a size of " +
                std::to_string((unsigned int)(dAlmansiStraindE_end - dAlmansiStraindE_begin)) +
                " but should have a size of " + std::to_string(dim * dim * dim * dim))

        TARDIGRADE_ERROR_TOOLS_CHECK(
            dim * dim * dim * dim == (unsigned int)(dAlmansiStraindF_end - dAlmansiStraindF_begin),
            "The derivative of the Almansi strain with respect to the deformation gradient has a size of " +
                std::to_string((unsigned int)(dAlmansiStraindF_end - dAlmansiStraindF_begin)) +
                " but should have a size of " + std::to_string(dim * dim * dim * dim))

        std::array<deformationGradient_type, dim * dim> inverseDeformationGradient;
        Eigen::Map<const Eigen::Matrix<deformationGradient_type, dim, dim, Eigen::RowMajor> > F(
            &(*deformationGradient_begin));
        Eigen::Map<Eigen::Matrix<deformationGradient_type, dim, dim, Eigen::RowMajor> > inv_F(
            inverseDeformationGradient.data());

        inv_F = F.inverse().eval();

        Eigen::Map<const Eigen::Matrix<greenLagrangeStrain_type, dim, dim, Eigen::RowMajor> > E(
            &(*greenLagrangeStrain_begin));
        Eigen::Map<Eigen::Matrix<almansiStrain_type, dim, dim, Eigen::RowMajor> > e(&(*almansiStrain_begin));

        // Map the Green-Lagrange strain to the current configuration
        e = (inv_F.transpose() * E * inv_F).eval();

        // Compute the jacobians
        for (unsigned int i = 0; i < dim; ++i) {
            for (unsigned int j = 0; j < dim; ++j) {
                for (unsigned int K = 0; K < dim; ++K) {
                    for (unsigned int L = 0; L < dim; ++L) {
                        *(dAlmansiStraindE_begin + dim * dim * dim * i + dim * dim * j + dim * K + L) =
                            inverseDeformationGradient[dim * K + i] * inverseDeformationGradient[dim * L + j];
                        *(dAlmansiStraindF_begin + dim * dim * dim * i + dim * dim * j + dim * K + L) =
                            -inverseDeformationGradient[dim * L + i] * (*(almansiStrain_begin + dim * K + j)) -
                            inverseDeformationGradient[dim * L + j] * (*(almansiStrain_begin + dim * i + K));
                    }
                }
            }
        }

        return;
    }

    template <unsigned int dim, class almansiStrain_iterator, class deformationGradient_iterator,
              class greenLagrangeStrain_iterator>
    void pullBackAlmansiStrain(const almansiStrain_iterator       &almansiStrain_begin,
                               const almansiStrain_iterator       &almansiStrain_end,
                               const deformationGradient_iterator &deformationGradient_begin,
                               const deformationGradient_iterator &deformationGradient_end,
                               greenLagrangeStrain_iterator        greenLagrangeStrain_begin,
                               greenLagrangeStrain_iterator        greenLagrangeStrain_end) {
        /*!
         * Pull back the almansi strain to the configuration indicated by the deformation gradient.
         *
         * \param &almansiStrain_begin: The starting iterator of the strain in the deformation gradient's current
         * configuration. \param &almansiStrain_end: The stopping iterator of the strain in the deformation gradient's
         * current configuration. \param &deformationGradient_begin: The starting iterator of the deformation gradient
         * between configurations. \param &deformationGradient_end: The stopping iterator of the deformation gradient
         * between configurations. \param &greenLagrangeStrain_begin: The starting iterator of the Green-Lagrange strain
         * which corresponds to the reference configuration of the deformation gradient. \param
         * &greenLagrangeStrain_end: The stopping iterator of the Green-Lagrange strain which corresponds to the
         * reference configuration of the deformation gradient.
         */

        using almansiStrain_type       = typename std::iterator_traits<almansiStrain_iterator>::value_type;
        using deformationGradient_type = typename std::iterator_traits<deformationGradient_iterator>::value_type;
        using greenLagrangeStrain_type = typename std::iterator_traits<greenLagrangeStrain_iterator>::value_type;

        TARDIGRADE_ERROR_TOOLS_CHECK(dim * dim == (unsigned int)(almansiStrain_end - almansiStrain_begin),
                                     "The Almansi strain has a size of " +
                                         std::to_string((unsigned int)(almansiStrain_end - almansiStrain_begin)) +
                                         " but should have a size of " + std::to_string(dim * dim))

        TARDIGRADE_ERROR_TOOLS_CHECK(dim * dim == (unsigned int)(deformationGradient_end - deformationGradient_begin),
                                     "The deformation gradient has a size of " +
                                         std::to_string((unsigned int)(deformationGradient_end -
                                                                       deformationGradient_begin)) +
                                         " but should have a size of " + std::to_string(dim * dim))

        TARDIGRADE_ERROR_TOOLS_CHECK(dim * dim == (unsigned int)(greenLagrangeStrain_end - greenLagrangeStrain_begin),
                                     "The Green-Lagrange strain has a size of " +
                                         std::to_string((unsigned int)(greenLagrangeStrain_end -
                                                                       greenLagrangeStrain_begin)) +
                                         " but should have a size of " + std::to_string(dim * dim))

        Eigen::Map<const Eigen::Matrix<deformationGradient_type, dim, dim, Eigen::RowMajor> > F(
            &(*deformationGradient_begin));

        Eigen::Map<Eigen::Matrix<greenLagrangeStrain_type, dim, dim, Eigen::RowMajor> > E(
            &(*greenLagrangeStrain_begin));
        Eigen::Map<const Eigen::Matrix<almansiStrain_type, dim, dim, Eigen::RowMajor> > e(&(*almansiStrain_begin));

        E = (F.transpose() * e * F).eval();

        return;
    }

    template <unsigned int dim, class almansiStrain_iterator, class deformationGradient_iterator,
              class greenLagrangeStrain_iterator, class dEde_iterator, class dEdF_iterator>
    void pullBackAlmansiStrain(const almansiStrain_iterator       &almansiStrain_begin,
                               const almansiStrain_iterator       &almansiStrain_end,
                               const deformationGradient_iterator &deformationGradient_begin,
                               const deformationGradient_iterator &deformationGradient_end,
                               greenLagrangeStrain_iterator        greenLagrangeStrain_begin,
                               greenLagrangeStrain_iterator greenLagrangeStrain_end, dEde_iterator dEde_begin,
                               dEde_iterator dEde_end, dEdF_iterator dEdF_begin, dEdF_iterator dEdF_end) {
        /*!
         * Pull back the almansi strain to the configuration indicated by the deformation gradient.
         *
         * \param &almansiStrain_begin: The starting iterator of the strain in the deformation gradient's current
         * configuration. \param &almansiStrain_end: The stopping iterator of the strain in the deformation gradient's
         * current configuration. \param &deformationGradient_begin: The starting iterator of the deformation gradient
         * between configurations. \param &deformationGradient_end: The stopping iterator of the deformation gradient
         * between configurations. \param greenLagrangeStrain_begin: The starting iterator of the Green-Lagrange strain
         * which corresponds to the reference configuration of the deformation gradient. \param greenLagrangeStrain_end:
         * The stopping iterator of the Green-Lagrange strain which corresponds to the reference configuration of the
         * deformation gradient. \param dEde_begin: The starting iterator of the Jacobian of the Green-Lagrange strain
         * w.r.t. the Almansi strain \param dEde_end: The stopping iterator of the Jacobian of the Green-Lagrange strain
         * w.r.t. the Almansi strain \param dEdF_begin: The starting iterator of the Jacobian of the Green-Lagrange
         * strain w.r.t. the deformation gradient \param dEdF_end: The stopping iterator of the Jacobian of the
         * Green-Lagrange strain w.r.t. the deformation gradient
         */

        using almansiStrain_type       = typename std::iterator_traits<almansiStrain_iterator>::value_type;
        using deformationGradient_type = typename std::iterator_traits<deformationGradient_iterator>::value_type;
        using greenLagrangeStrain_type = typename std::iterator_traits<greenLagrangeStrain_iterator>::value_type;
        using dEde_type                = typename std::iterator_traits<dEde_iterator>::value_type;
        using dEdF_type                = typename std::iterator_traits<dEdF_iterator>::value_type;

        TARDIGRADE_ERROR_TOOLS_CHECK(dim * dim == (unsigned int)(almansiStrain_end - almansiStrain_begin),
                                     "The Almansi strain has a size of " +
                                         std::to_string((unsigned int)(almansiStrain_end - almansiStrain_begin)) +
                                         " but should have a size of " + std::to_string(dim * dim))

        TARDIGRADE_ERROR_TOOLS_CHECK(dim * dim == (unsigned int)(deformationGradient_end - deformationGradient_begin),
                                     "The deformation gradient has a size of " +
                                         std::to_string((unsigned int)(deformationGradient_end -
                                                                       deformationGradient_begin)) +
                                         " but should have a size of " + std::to_string(dim * dim))

        TARDIGRADE_ERROR_TOOLS_CHECK(dim * dim == (unsigned int)(greenLagrangeStrain_end - greenLagrangeStrain_begin),
                                     "The Green-Lagrange strain has a size of " +
                                         std::to_string((unsigned int)(greenLagrangeStrain_end -
                                                                       greenLagrangeStrain_begin)) +
                                         " but should have a size of " + std::to_string(dim * dim))

        Eigen::Map<const Eigen::Matrix<deformationGradient_type, dim, dim, Eigen::RowMajor> > F(
            &(*deformationGradient_begin));

        Eigen::Map<Eigen::Matrix<greenLagrangeStrain_type, dim, dim, Eigen::RowMajor> > E(
            &(*greenLagrangeStrain_begin));
        Eigen::Map<const Eigen::Matrix<almansiStrain_type, dim, dim, Eigen::RowMajor> > e(&(*almansiStrain_begin));

        E = (F.transpose() * e * F).eval();
        std::fill(dEde_begin, dEde_end, dEde_type());
        std::fill(dEdF_begin, dEdF_end, dEdF_type());

        for (unsigned int I = 0; I < dim; ++I) {
            for (unsigned int J = 0; J < dim; ++J) {
                for (unsigned int K = 0; K < dim; ++K) {
                    for (unsigned int L = 0; L < dim; ++L) {
                        *(dEde_begin + dim * dim * dim * I + dim * dim * J + dim * K + L) =
                            (*(deformationGradient_begin + dim * K + I)) * (*(deformationGradient_begin + dim * L + J));
                        *(dEdF_begin + dim * dim * dim * I + dim * dim * J + dim * K + I) +=
                            (*(almansiStrain_begin + dim * K + L)) * (*(deformationGradient_begin + dim * L + J));
                        *(dEdF_begin + dim * dim * dim * I + dim * dim * J + dim * K + J) +=
                            (*(deformationGradient_begin + dim * L + I)) * (*(almansiStrain_begin + dim * L + K));
                    }
                }
            }
        }
    }

    template <unsigned int dim, class A_iterator, class symmA_iterator>
    void computeSymmetricPart(const A_iterator &A_begin, const A_iterator &A_end, symmA_iterator symmA_begin,
                              symmA_iterator symmA_end) {
        /*!
         * Compute the symmetric part of a second order tensor ( \f$A\f$ ) and return it.
         *
         * \f$symm( A )_ij = \frac{1}{2}\left(A_{ij} + A_{ji}\right)\f$
         *
         * \param &A_begin: A constant reference to the starting iterator of the second order tensor to process (
         * \f$A\f$ ) \param &A_end: A constant reference to the stopping iterator of the second order tensor to process
         * ( \f$A\f$ ) \param symmA_begin: The starting iterator of the symmetric part of A ( \f$A^{symm}\f$ ) \param
         * symmA_end: The stopping iterator of the symmetric part of A ( \f$A^{symm}\f$ )
         */

        TARDIGRADE_ERROR_TOOLS_CHECK((unsigned int)(A_end - A_begin) == dim * dim,
                                     "A has a size of " + std::to_string((unsigned int)(A_end - A_begin)) +
                                         " but the nearest square matrix has a row and column size of " +
                                         std::to_string(dim))

        TARDIGRADE_ERROR_TOOLS_CHECK((unsigned int)(symmA_end - symmA_begin) == dim * dim,
                                     "symmA has a size of " + std::to_string((unsigned int)(symmA_end - symmA_begin)) +
                                         " but should have a size of " + std::to_string(dim * dim))

        for (unsigned int i = 0; i < dim; ++i) {
            for (unsigned int j = 0; j < dim; ++j) {
                *(symmA_begin + dim * i + j) = 0.5 * ((*(A_begin + dim * i + j)) + (*(A_begin + dim * j + i)));
            }
        }
    }

    template <unsigned int dim, class A_iterator, class symmA_iterator, class dSymmAdA_iterator>
    void computeSymmetricPart(const A_iterator &A_begin, const A_iterator &A_end, symmA_iterator symmA_begin,
                              symmA_iterator symmA_end, dSymmAdA_iterator dSymmAdA_begin,
                              dSymmAdA_iterator dSymmAdA_end) {
        /*!
         * Compute the symmetric part of a second order tensor ( \f$A\f$ ) and return it.
         *
         * \f$symm( A )_ij = \frac{1}{2}\left(A_{ij} + A_{ji}\right)\f$
         *
         * \param &A_begin: A constant reference to the starting iterator of the second order tensor to process (
         * \f$A\f$ ) \param &A_end: A constant reference to the stopping iterator of the second order tensor to process
         * ( \f$A\f$ ) \param symmA_begin: The starting iterator of the symmetric part of A ( \f$A^{symm}\f$ ) \param
         * symmA_end: The stopping iterator of the symmetric part of A ( \f$A^{symm}\f$ ) \param dSymmAdA_begin: The
         * starting iterator of the symmetric part of the derivative of the symmetric part of A with respect to A (
         * \f$\frac{\partial A^{symm}}{\partial A}\f$ ) \param dSymmAdA_end: The stopping iterator of the symmetric part
         * of the derivative of the symmetric part of A with respect to A ( \f$\frac{\partial A^{symm}}{\partial A}\f$ )
         */

        using dSymmAdA_type = typename std::iterator_traits<dSymmAdA_iterator>::value_type;

        TARDIGRADE_ERROR_TOOLS_CHECK((unsigned int)(dSymmAdA_end - dSymmAdA_begin) == dim * dim * dim * dim,
                                     "dSymmAdA has a size of " +
                                         std::to_string((unsigned int)(dSymmAdA_end - dSymmAdA_begin)) +
                                         " but must have a size of " + std::to_string(dim * dim * dim * dim))

        TARDIGRADE_ERROR_TOOLS_CATCH(computeSymmetricPart<dim>(A_begin, A_end, symmA_begin, symmA_end));

        std::fill(dSymmAdA_begin, dSymmAdA_end, dSymmAdA_type());

        for (unsigned int i = 0; i < dim; ++i) {
            for (unsigned int j = 0; j < dim; ++j) {
                *(dSymmAdA_begin + dim * dim * dim * i + dim * dim * j + dim * i + j) += 0.5;
                *(dSymmAdA_begin + dim * dim * dim * i + dim * dim * j + dim * j + i) += 0.5;
            }
        }

        return;
    }

    template <class A_iterator, class symmA_iterator>
    void computeSymmetricPart(const A_iterator &A_begin, const A_iterator &A_end, symmA_iterator symmA_begin,
                              symmA_iterator symmA_end, unsigned int &dim) {
        /*!
         * Compute the symmetric part of a second order tensor ( \f$A\f$ ) and return it.
         *
         * \f$symm( A )_ij = \frac{1}{2}\left(A_{ij} + A_{ji}\right)\f$
         *
         * \param &A_begin: A constant reference to the starting iterator of the second order tensor to process (
         * \f$A\f$ ) \param &A_end: A constant reference to the stopping iterator of the second order tensor to process
         * ( \f$A\f$ ) \param symmA_begin: The starting iterator of the symmetric part of A ( \f$A^{symm}\f$ ) \param
         * symmA_end: The stopping iterator of the symmetric part of A ( \f$A^{symm}\f$ ) \param &dim: The dimension of
         * A. Note that this is an output used for help with computing the Jacobian. If you don't need dim as an output
         * use the version of this function without it.
         */

        // Get the dimension of A
        dim = (unsigned int)std::pow((unsigned int)(A_end - A_begin), 0.5);

        TARDIGRADE_ERROR_TOOLS_CHECK((unsigned int)(A_end - A_begin) == dim * dim,
                                     "A has a size of " + std::to_string((unsigned int)(A_end - A_begin)) +
                                         " but the nearest square matrix has a row and column size of " +
                                         std::to_string(dim))

        TARDIGRADE_ERROR_TOOLS_CHECK((unsigned int)(symmA_end - symmA_begin) == dim * dim,
                                     "symmA has a size of " + std::to_string((unsigned int)(symmA_end - symmA_begin)) +
                                         " but should have a size of " + std::to_string(dim * dim))

        for (unsigned int i = 0; i < dim; ++i) {
            for (unsigned int j = 0; j < dim; ++j) {
                *(symmA_begin + dim * i + j) = 0.5 * ((*(A_begin + dim * i + j)) + (*(A_begin + dim * j + i)));
            }
        }
    }

    template <class A_iterator, class symmA_iterator, class dSymmAdA_iterator>
    void computeSymmetricPart(const A_iterator &A_begin, const A_iterator &A_end, symmA_iterator symmA_begin,
                              symmA_iterator symmA_end, dSymmAdA_iterator dSymmAdA_begin,
                              dSymmAdA_iterator dSymmAdA_end) {
        /*!
         * Compute the symmetric part of a second order tensor ( \f$A\f$ ) and return it.
         *
         * \f$symm( A )_ij = \frac{1}{2}\left(A_{ij} + A_{ji}\right)\f$
         *
         * \param &A_begin: A constant reference to the starting iterator of the second order tensor to process (
         * \f$A\f$ ) \param &A_end: A constant reference to the stopping iterator of the second order tensor to process
         * ( \f$A\f$ ) \param symmA_begin: The starting iterator of the symmetric part of A ( \f$A^{symm}\f$ ) \param
         * symmA_end: The stopping iterator of the symmetric part of A ( \f$A^{symm}\f$ ) \param dSymmAdA_begin: The
         * starting iterator of the symmetric part of the derivative of the symmetric part of A with respect to A (
         * \f$\frac{\partial A^{symm}}{\partial A}\f$ ) \param dSymmAdA_end: The stopping iterator of the symmetric part
         * of the derivative of the symmetric part of A with respect to A ( \f$\frac{\partial A^{symm}}{\partial A}\f$ )
         */

        using dSymmAdA_type = typename std::iterator_traits<dSymmAdA_iterator>::value_type;

        TARDIGRADE_ERROR_TOOLS_EVAL(const unsigned int Asize = (unsigned int)(A_end - A_begin);)

        TARDIGRADE_ERROR_TOOLS_CHECK((unsigned int)(dSymmAdA_end - dSymmAdA_begin) == Asize * Asize,
                                     "dSymmAdA has a size of " +
                                         std::to_string((unsigned int)(dSymmAdA_end - dSymmAdA_begin)) +
                                         " but must have a size of " + std::to_string(Asize * Asize))

        unsigned int dim;
        TARDIGRADE_ERROR_TOOLS_CATCH(computeSymmetricPart(A_begin, A_end, symmA_begin, symmA_end, dim));

        std::fill(dSymmAdA_begin, dSymmAdA_end, dSymmAdA_type());

        for (unsigned int i = 0; i < dim; ++i) {
            for (unsigned int j = 0; j < dim; ++j) {
                *(dSymmAdA_begin + dim * dim * dim * i + dim * dim * j + dim * i + j) += 0.5;
                *(dSymmAdA_begin + dim * dim * dim * i + dim * dim * j + dim * j + i) += 0.5;
            }
        }

        return;
    }

    template <unsigned int dim, class PK2_iterator, class F_iterator, class cauchyStress_iterator>
    void pushForwardPK2Stress(const PK2_iterator &PK2_begin, const PK2_iterator &PK2_end, const F_iterator &F_begin,
                              const F_iterator &F_end, cauchyStress_iterator cauchyStress_begin,
                              cauchyStress_iterator cauchyStress_end) {
        /*!
         * Push the Second Piola-Kirchhoff stress forward to the current configuration resulting in the Cauchy stress
         *
         * \f$ \sigma_{ij} = \frac{1}{J} F_{iI} S_{IJ} F_{jJ} \f$
         *
         * \param &PK2_begin: The starting iterator of the Second Piola-Kirchhoff stress \f$ S_{IJ} \f$
         * \param &PK2_end: The stopping iterator of the Second Piola-Kirchhoff stress \f$ S_{IJ} \f$
         * \param &F_begin: The starting iterator of the deformation gradient \f$ F_{iI} \f$
         * \param &F_end: The stopping iterator of the deformation gradient \f$ F_{iI} \f$
         * \param &cauchyStress_begin: The starting iterator of the Cauchy stress \f$ \sigma_{ij} \f$
         * \param &cauchyStress_end: The stopping iterator of the Cauchy stress \f$ \sigma_{ij} \f$
         */

        using PK2_type          = typename std::iterator_traits<PK2_iterator>::value_type;
        using F_type            = typename std::iterator_traits<F_iterator>::value_type;
        using cauchyStress_type = typename std::iterator_traits<cauchyStress_iterator>::value_type;

        TARDIGRADE_ERROR_TOOLS_CHECK(dim * dim == (unsigned int)(PK2_end - PK2_begin),
                                     "The PK2 stress has a size of " +
                                         std::to_string((unsigned int)(PK2_end - PK2_begin)) +
                                         " but must have a size of " + std::to_string(dim * dim))

        TARDIGRADE_ERROR_TOOLS_CHECK(dim * dim == (unsigned int)(F_end - F_begin),
                                     "The deformation gradient has a size of " +
                                         std::to_string((unsigned int)(F_end - F_begin)) + " but must have a size of " +
                                         std::to_string(dim * dim))

        TARDIGRADE_ERROR_TOOLS_CHECK(dim * dim == (unsigned int)(cauchyStress_end - cauchyStress_begin),
                                     "The Cauchy stress has a size of " +
                                         std::to_string((unsigned int)(cauchyStress_end - cauchyStress_begin)) +
                                         " but must have a size of " + std::to_string(dim * dim))

        Eigen::Map<const Eigen::Matrix<PK2_type, dim, dim, Eigen::RowMajor> >    PK2_map(&(*PK2_begin));
        Eigen::Map<const Eigen::Matrix<F_type, dim, dim, Eigen::RowMajor> >      F_map(&(*F_begin));
        Eigen::Map<Eigen::Matrix<cauchyStress_type, dim, dim, Eigen::RowMajor> > cauchyStress_map(
            &(*cauchyStress_begin));

        F_type J = F_map.determinant();

        cauchyStress_map = (F_map * PK2_map * F_map.transpose() / J).eval();

        return;
    }

    template <unsigned int dim, class PK2_iterator, class F_iterator, class cauchyStress_iterator,
              class dCauchyStressdPK2_iterator, class dCauchyStressdF_iterator>
    void pushForwardPK2Stress(const PK2_iterator &PK2_begin, const PK2_iterator &PK2_end, const F_iterator &F_begin,
                              const F_iterator &F_end, cauchyStress_iterator cauchyStress_begin,
                              cauchyStress_iterator      cauchyStress_end,
                              dCauchyStressdPK2_iterator dCauchyStressdPK2_begin,
                              dCauchyStressdPK2_iterator dCauchyStressdPK2_end,
                              dCauchyStressdF_iterator   dCauchyStressdF_begin,
                              dCauchyStressdF_iterator   dCauchyStressdF_end) {
        /*!
         * Push the Second Piola-Kirchhoff stress forward to the current configuration resulting in the Cauchy stress
         *
         * \f$ \sigma_{ij} = \frac{1}{J} F_{iI} S_{IJ} F_{jJ} \f$
         *
         * \param &PK2_begin: The starting iterator of the Second Piola-Kirchhoff stress \f$ S_{IJ} \f$
         * \param &PK2_end: The stopping iterator of the Second Piola-Kirchhoff stress \f$ S_{IJ} \f$
         * \param &F_begin: The starting iterator of the deformation gradient \f$ F_{iI} \f$
         * \param &F_end: The stopping iterator of the deformation gradient \f$ F_{iI} \f$
         * \param &cauchyStress_begin: The starting iterator of the Cauchy stress \f$ \sigma_{ij} \f$
         * \param &cauchyStress_end: The stopping iterator of the Cauchy stress \f$ \sigma_{ij} \f$
         * \param &dCauchyStressdPK2_begin: The starting iterator of the gradient of the Cauchy stress w.r.t. the PK2
         * stress \param &dCauchyStressdPK2_end: The stopping iterator of the gradient of the Cauchy stress w.r.t. the
         * PK2 stress \param &dCauchyStressdF_begin: The starting iterator of the gradient of the Cauchy stress w.r.t.
         * the deformation gradient \param &dCauchyStressdF_end: The stopping iterator of the gradient of the Cauchy
         * stress w.r.t. the deformation gradient
         */

        using PK2_type               = typename std::iterator_traits<PK2_iterator>::value_type;
        using F_type                 = typename std::iterator_traits<F_iterator>::value_type;
        using cauchyStress_type      = typename std::iterator_traits<cauchyStress_iterator>::value_type;
        using dCauchyStressdPK2_type = typename std::iterator_traits<dCauchyStressdPK2_iterator>::value_type;
        using dCauchyStressdF_type   = typename std::iterator_traits<dCauchyStressdF_iterator>::value_type;

        TARDIGRADE_ERROR_TOOLS_CHECK(dim * dim == (unsigned int)(PK2_end - PK2_begin),
                                     "The PK2 stress has a size of " +
                                         std::to_string((unsigned int)(PK2_end - PK2_begin)) +
                                         " but must have a size of " + std::to_string(dim * dim))

        TARDIGRADE_ERROR_TOOLS_CHECK(dim * dim == (unsigned int)(F_end - F_begin),
                                     "The deformation gradient has a size of " +
                                         std::to_string((unsigned int)(F_end - F_begin)) + " but must have a size of " +
                                         std::to_string(dim * dim))

        TARDIGRADE_ERROR_TOOLS_CHECK(dim * dim == (unsigned int)(cauchyStress_end - cauchyStress_begin),
                                     "The Cauchy stress has a size of " +
                                         std::to_string((unsigned int)(cauchyStress_end - cauchyStress_begin)) +
                                         " but must have a size of " + std::to_string(dim * dim))

        TARDIGRADE_ERROR_TOOLS_CHECK(
            dim * dim * dim * dim == (unsigned int)(dCauchyStressdPK2_end - dCauchyStressdPK2_begin),
            "The derivative of the Cauchy stress w.r.t. the PK2 stress has a size of " +
                std::to_string((unsigned int)(dCauchyStressdPK2_end - dCauchyStressdPK2_begin)) +
                " but must have a size of " + std::to_string(dim * dim * dim * dim))

        TARDIGRADE_ERROR_TOOLS_CHECK(
            dim * dim * dim * dim == (unsigned int)(dCauchyStressdF_end - dCauchyStressdF_begin),
            "The derivative of the Cauchy stress w.r.t. the deformation gradient has a size of " +
                std::to_string((unsigned int)(dCauchyStressdF_end - dCauchyStressdF_begin)) +
                " but must have a size of " + std::to_string(dim * dim * dim * dim))

        std::array<F_type, dim * dim> dJdF;

        Eigen::Map<const Eigen::Matrix<PK2_type, dim, dim, Eigen::RowMajor> >    PK2_map(&(*PK2_begin));
        Eigen::Map<const Eigen::Matrix<F_type, dim, dim, Eigen::RowMajor> >      F_map(&(*F_begin));
        Eigen::Map<Eigen::Matrix<cauchyStress_type, dim, dim, Eigen::RowMajor> > cauchyStress_map(
            &(*cauchyStress_begin));
        Eigen::Map<Eigen::Matrix<F_type, dim, dim, Eigen::RowMajor> > dJdF_map(dJdF.data());

        F_type J = F_map.determinant();

        cauchyStress_map = (F_map * PK2_map * F_map.transpose() / J).eval();
        dJdF_map         = (J * F_map.inverse().transpose()).eval();

        std::fill(dCauchyStressdPK2_begin, dCauchyStressdPK2_end, dCauchyStressdPK2_type());
        std::fill(dCauchyStressdF_begin, dCauchyStressdF_end, dCauchyStressdF_type());

        for (unsigned int i = 0; i < dim; ++i) {
            for (unsigned int j = 0; j < dim; ++j) {
                for (unsigned int A = 0; A < dim; ++A) {
                    for (unsigned int B = 0; B < dim; ++B) {
                        *(dCauchyStressdPK2_begin + dim * dim * dim * i + dim * dim * j + dim * A + B) +=
                            (*(F_begin + dim * i + A)) * (*(F_begin + dim * j + B));

                        *(dCauchyStressdF_begin + dim * dim * dim * i + dim * dim * j + dim * A + B) -=
                            (*(cauchyStress_begin + dim * i + j)) * dJdF[dim * A + B];

                        *(dCauchyStressdF_begin + dim * dim * dim * i + dim * dim * j + dim * i + A) +=
                            (*(PK2_begin + dim * A + B)) * (*(F_begin + dim * j + B));

                        *(dCauchyStressdF_begin + dim * dim * dim * i + dim * dim * j + dim * j + A) +=
                            (*(F_begin + dim * i + B)) * (*(PK2_begin + dim * B + A));
                    }
                }
            }
        }

        std::transform(dCauchyStressdPK2_begin, dCauchyStressdPK2_end, dCauchyStressdPK2_begin,
                       std::bind(std::divides<>(), std::placeholders::_1, J));

        std::transform(dCauchyStressdF_begin, dCauchyStressdF_end, dCauchyStressdF_begin,
                       std::bind(std::divides<>(), std::placeholders::_1, J));

        return;
    }

    template <unsigned int dim, class cauchyStress_iterator, class F_iterator, class PK2_iterator>
    void pullBackCauchyStress(const cauchyStress_iterator &cauchyStress_begin,
                              const cauchyStress_iterator &cauchyStress_end, const F_iterator &F_begin,
                              const F_iterator &F_end, PK2_iterator PK2_begin, PK2_iterator PK2_end) {
        /*!
         * Pull back the Cauchy stress to an earlier configuration resulting in the second Piola-Kirchhoff stress
         *
         * \f$ S_{IJ} = J F^{-1}_{Ii} \sigma_{ij} F^{-1}_{Jj} \f$
         *
         * where \f$S_{IJ}\f$ are the components of the second Piola-Kirchhoff stress tensor, \f$J \f$ is the
         * determinant of the deformation gradient \f$\bf{F}\f$ which has components \f$F_{iI}\f$, and
         * \f$ \sigma_{ij} \f$ are the components of the Cauchy stress.
         *
         * \param &cauchyStress_begin: The starting iterator of the cauchy stress tensor in row-major form (all nine
         * components) \param &cauchyStress_end: The stopping iterator of the cauchy stress tensor in row-major form
         * (all nine components) \param &F_begin: The starting iterator of the deformation gradient \param &F_end: The
         * stopping iterator of the deformation gradient \param &PK2_begin: The starting iterator of the resulting
         * second Piola-Kirchhoff stress \param &PK2_end: The stopping iterator of the resulting second Piola-Kirchhoff
         * stress
         */

        using cauchyStress_type = typename std::iterator_traits<cauchyStress_iterator>::value_type;
        using F_type            = typename std::iterator_traits<F_iterator>::value_type;
        using PK2_type          = typename std::iterator_traits<PK2_iterator>::value_type;

        TARDIGRADE_ERROR_TOOLS_CHECK(dim * dim == (unsigned int)(cauchyStress_end - cauchyStress_begin),
                                     "The Cauchy stress has a size of " +
                                         std::to_string((unsigned int)(cauchyStress_end - cauchyStress_begin)) +
                                         " but must have a size of " + std::to_string(dim * dim))

        TARDIGRADE_ERROR_TOOLS_CHECK(dim * dim == (unsigned int)(F_end - F_begin),
                                     "The deformation gradient has a size of " +
                                         std::to_string((unsigned int)(F_end - F_begin)) + " but must have a size of " +
                                         std::to_string(dim * dim))

        TARDIGRADE_ERROR_TOOLS_CHECK(dim * dim == (unsigned int)(PK2_end - PK2_begin),
                                     "The PK2 stress has a size of " +
                                         std::to_string((unsigned int)(PK2_end - PK2_begin)) +
                                         " but must have a size of " + std::to_string(dim * dim))

        std::array<F_type, dim * dim> F_inv;

        Eigen::Map<const Eigen::Matrix<cauchyStress_type, dim, dim, Eigen::RowMajor> > cauchyStress_map(
            &(*cauchyStress_begin));
        Eigen::Map<const Eigen::Matrix<F_type, dim, dim, Eigen::RowMajor> > F_map(&(*F_begin));
        Eigen::Map<Eigen::Matrix<F_type, dim, dim, Eigen::RowMajor> >       F_inv_map(F_inv.data());
        Eigen::Map<Eigen::Matrix<PK2_type, dim, dim, Eigen::RowMajor> >     PK2_map(&(*PK2_begin));

        F_inv_map = F_map.inverse();
        F_type J  = F_map.determinant();

        PK2_map = (J * F_inv_map * cauchyStress_map * F_inv_map.transpose()).eval();

        return;
    }

    template <unsigned int dim, class cauchyStress_iterator, class F_iterator, class PK2_iterator,
              class dPK2dCauchyStress_iterator, class dPK2dF_iterator>
    void pullBackCauchyStress(const cauchyStress_iterator &cauchyStress_begin,
                              const cauchyStress_iterator &cauchyStress_end, const F_iterator &F_begin,
                              const F_iterator &F_end, PK2_iterator PK2_begin, PK2_iterator PK2_end,
                              dPK2dCauchyStress_iterator dPK2dCauchyStress_begin,
                              dPK2dCauchyStress_iterator dPK2dCauchyStress_end, dPK2dF_iterator dPK2dF_begin,
                              dPK2dF_iterator dPK2dF_end) {
        /*!
         * Pull back the Cauchy stress to an earlier configuration resulting in the second Piola-Kirchhoff stress
         *
         * \f$ S_{IJ} = J F^{-1}_{Ii} \sigma_{ij} F^{-1}_{Jj} \f$
         *
         * where \f$S_{IJ}\f$ are the components of the second Piola-Kirchhoff stress tensor, \f$J \f$ is the
         * determinant of the deformation gradient \f$\bf{F}\f$ which has components \f$F_{iI}\f$, and
         * \f$ \sigma_{ij} \f$ are the components of the Cauchy stress.
         *
         * \param &cauchyStress_begin: The starting iterator of the cauchy stress tensor in row-major form (all nine
         * components) \param &cauchyStress_end: The stopping iterator of the cauchy stress tensor in row-major form
         * (all nine components) \param &F_begin: The starting iterator of the deformation gradient \param &F_end: The
         * stopping iterator of the deformation gradient \param &PK2_begin: The starting iterator of the resulting
         * second Piola-Kirchhoff stress \param &PK2_end: The stopping iterator of the resulting second Piola-Kirchhoff
         * stress \param &dPK2dCauchyStress_begin: The starting iterator of the directional derivative of the second
         * Piola-Kirchhoff stress tensor w.r.t. the Cauchy stress \param &dPK2dCauchyStress_end: The stopping iterator
         * of the directional derivative of the second Piola-Kirchhoff stress tensor w.r.t. the Cauchy stress \param
         * &dPK2dF_begin: The starting iterator of the directional derivative of the second Piola-Kirchhoff stress
         * tensor w.r.t. the deformation gradient \param &dPK2dF_end: The stopping iterator of the directional
         * derivative of the second Piola-Kirchhoff stress tensor w.r.t. the deformation gradient
         */

        using cauchyStress_type      = typename std::iterator_traits<cauchyStress_iterator>::value_type;
        using F_type                 = typename std::iterator_traits<F_iterator>::value_type;
        using PK2_type               = typename std::iterator_traits<PK2_iterator>::value_type;
        using dPK2dCauchyStress_type = typename std::iterator_traits<dPK2dCauchyStress_iterator>::value_type;
        using dPK2dF_type            = typename std::iterator_traits<dPK2dF_iterator>::value_type;

        TARDIGRADE_ERROR_TOOLS_CHECK(dim * dim == (unsigned int)(cauchyStress_end - cauchyStress_begin),
                                     "The Cauchy stress has a size of " +
                                         std::to_string((unsigned int)(cauchyStress_end - cauchyStress_begin)) +
                                         " but must have a size of " + std::to_string(dim * dim))

        TARDIGRADE_ERROR_TOOLS_CHECK(dim * dim == (unsigned int)(F_end - F_begin),
                                     "The deformation gradient has a size of " +
                                         std::to_string((unsigned int)(F_end - F_begin)) + " but must have a size of " +
                                         std::to_string(dim * dim))

        TARDIGRADE_ERROR_TOOLS_CHECK(dim * dim == (unsigned int)(PK2_end - PK2_begin),
                                     "The PK2 stress has a size of " +
                                         std::to_string((unsigned int)(PK2_end - PK2_begin)) +
                                         " but must have a size of " + std::to_string(dim * dim))

        TARDIGRADE_ERROR_TOOLS_CHECK(
            dim * dim * dim * dim == (unsigned int)(dPK2dCauchyStress_end - dPK2dCauchyStress_begin),
            "The derivative of the PK2 stress w.r.t. the Cauchy stress has a size of " +
                std::to_string((unsigned int)(dPK2dCauchyStress_end - dPK2dCauchyStress_begin)) +
                " but must have a size of " + std::to_string(dim * dim * dim * dim))

        TARDIGRADE_ERROR_TOOLS_CHECK(dim * dim * dim * dim == (unsigned int)(dPK2dF_end - dPK2dF_begin),
                                     "The derivative of the PK2 stress w.r.t. the deformation gradient has a size of " +
                                         std::to_string((unsigned int)(dPK2dF_end - dPK2dF_begin)) +
                                         " but must have a size of " + std::to_string(dim * dim * dim * dim))

        std::array<F_type, dim * dim> F_inv;

        Eigen::Map<const Eigen::Matrix<cauchyStress_type, dim, dim, Eigen::RowMajor> > cauchyStress_map(
            &(*cauchyStress_begin));
        Eigen::Map<const Eigen::Matrix<F_type, dim, dim, Eigen::RowMajor> > F_map(&(*F_begin));
        Eigen::Map<Eigen::Matrix<F_type, dim, dim, Eigen::RowMajor> >       F_inv_map(F_inv.data());
        Eigen::Map<Eigen::Matrix<PK2_type, dim, dim, Eigen::RowMajor> >     PK2_map(&(*PK2_begin));

        F_inv_map = F_map.inverse();
        F_type J  = F_map.determinant();

        PK2_map = (J * F_inv_map * cauchyStress_map * F_inv_map.transpose()).eval();

        std::fill(dPK2dCauchyStress_begin, dPK2dCauchyStress_end, dPK2dCauchyStress_type());
        std::fill(dPK2dF_begin, dPK2dF_end, dPK2dF_type());

        for (unsigned int A = 0; A < dim; A++) {
            for (unsigned int B = 0; B < dim; B++) {
                for (unsigned int k = 0; k < dim; k++) {
                    for (unsigned int l = 0; l < dim; l++) {
                        *(dPK2dCauchyStress_begin + dim * dim * dim * A + dim * dim * B + dim * k + l) =
                            J * F_inv[dim * A + k] * F_inv[dim * B + l];

                        *(dPK2dF_begin + dim * dim * dim * A + dim * dim * B + dim * k + l) =
                            F_inv[dim * l + k] * (*(PK2_begin + dim * A + B)) -
                            F_inv[dim * A + k] * (*(PK2_begin + dim * l + B)) -
                            F_inv[dim * B + k] * (*(PK2_begin + dim * A + l));
                    }
                }
            }
        }

        return;
    }

    template <unsigned int dim, typename Dt_type, class previousDeformationGradient_iterator, class Lp_iterator,
              class L_iterator, class deformationGradient_iterator, typename alpha_type>
    void evolveFExponentialMap(const Dt_type                              &Dt,
                               const previousDeformationGradient_iterator &previousDeformationGradient_begin,
                               const previousDeformationGradient_iterator &previousDeformationGradient_end,
                               const Lp_iterator &Lp_begin, const Lp_iterator &Lp_end, const L_iterator &L_begin,
                               const L_iterator &L_end, deformationGradient_iterator deformationGradient_begin,
                               deformationGradient_iterator deformationGradient_end, const alpha_type alpha) {
        /*!
         * Evolve the deformation gradient using the exponential map. Assumes the evolution equation is of the form
         *
         * \f$ \dot{F}_{iI} = \ell_{ij} F_{jI} \f$
         *
         * \param &Dt: The change in time
         * \param &previousDeformationGradient_begin: The starting iterator of the previous value of the deformation
         * gradient \param &previousDeformationGradient_end: The stopping iterator of the previous value of the
         * deformation gradient \param &Lp_begin: The starting iterator of the previous value of the velocity gradient
         * \param &Lp_end: The stopping iterator of the previous value of the velocity gradient
         * \param &L_begin: The starting iterator of the current value of the velocity gradient
         * \param &L_end: The stopping iterator of the current value of the velocity gradient
         * \param &deformationGradient_begin: The starting iterator of the computed value of the deformation gradient
         * \param &deformationGradient_end: The stopping iterator of the computed value of the deformation gradient
         * \param &alpha: The integration parameter (0 is explicit and 1 is implicit)
         */

        using deformationGradient_type = typename std::iterator_traits<deformationGradient_iterator>::value_type;

        TARDIGRADE_ERROR_TOOLS_CHECK(
            (unsigned int)(previousDeformationGradient_end - previousDeformationGradient_begin),
            "The previous deformation gradient has a size of " +
                std::to_string((unsigned int)(previousDeformationGradient_end - previousDeformationGradient_begin)) +
                " but must have a size of " + std::to_string(dim * dim))

        TARDIGRADE_ERROR_TOOLS_CHECK((unsigned int)(Lp_end - Lp_begin),
                                     "The previous velocity gradient has a size of " +
                                         std::to_string((unsigned int)(Lp_end - Lp_begin)) +
                                         " but must have a size of " + std::to_string(dim * dim))

        TARDIGRADE_ERROR_TOOLS_CHECK((unsigned int)(L_end - L_begin),
                                     "The velocity gradient has a size of " +
                                         std::to_string((unsigned int)(L_end - L_begin)) + " but must have a size of " +
                                         std::to_string(dim * dim))

        TARDIGRADE_ERROR_TOOLS_CHECK((unsigned int)(deformationGradient_end - deformationGradient_begin),
                                     "The velocity gradient has a size of " +
                                         std::to_string((unsigned int)(L_end - L_begin)) + " but must have a size of " +
                                         std::to_string(dim * dim))

        std::array<deformationGradient_type, dim * dim> DtLalpha, expDtLalpha, tempVector1, tempVector2, tempVector3;

        for (unsigned int i = 0; i < dim * dim; ++i) {
            DtLalpha[i] = Dt * ((1 - alpha) * (*(Lp_begin + i)) + alpha * (*(L_begin + i)));
        }

        TARDIGRADE_ERROR_TOOLS_CATCH(tardigradeVectorTools::computeMatrixExponentialScalingAndSquaring(
            std::begin(DtLalpha), std::end(DtLalpha), dim, std::begin(tempVector1), std::end(tempVector1),
            std::begin(tempVector2), std::end(tempVector2), std::begin(tempVector3), std::end(tempVector3),
            std::begin(expDtLalpha), std::end(expDtLalpha)))

        std::fill(deformationGradient_begin, deformationGradient_end, deformationGradient_type());

        for (unsigned int i = 0; i < dim; ++i) {
            for (unsigned int j = 0; j < dim; ++j) {
                for (unsigned int k = 0; k < dim; ++k) {
                    *(deformationGradient_begin + dim * i + k) +=
                        expDtLalpha[dim * i + j] * (*(previousDeformationGradient_begin + dim * j + k));
                }
            }
        }
    }

    template <unsigned int dim, typename Dt_type, class previousDeformationGradient_iterator, class Lp_iterator,
              class L_iterator, class deformationGradient_iterator, class dFdL_iterator, typename alpha_type>
    void evolveFExponentialMap(const Dt_type                              &Dt,
                               const previousDeformationGradient_iterator &previousDeformationGradient_begin,
                               const previousDeformationGradient_iterator &previousDeformationGradient_end,
                               const Lp_iterator &Lp_begin, const Lp_iterator &Lp_end, const L_iterator &L_begin,
                               const L_iterator &L_end, deformationGradient_iterator deformationGradient_begin,
                               deformationGradient_iterator deformationGradient_end, dFdL_iterator dFdL_begin,
                               dFdL_iterator dFdL_end, const alpha_type alpha) {
        /*!
         * Evolve the deformation gradient using the exponential map. Assumes the evolution equation is of the form
         *
         * \f$ \dot{F}_{iI} = \ell_{ij} F_{jI} \f$
         *
         * \param &Dt: The change in time
         * \param &previousDeformationGradient_begin: The starting iterator of the previous value of the deformation
         * gradient \param &previousDeformationGradient_end: The stopping iterator of the previous value of the
         * deformation gradient \param &Lp_begin: The starting iterator of the previous value of the velocity gradient
         * \param &Lp_end: The stopping iterator of the previous value of the velocity gradient
         * \param &L_begin: The starting iterator of the current value of the velocity gradient
         * \param &L_end: The stopping iterator of the current value of the velocity gradient
         * \param &deformationGradient_begin: The starting iterator of the computed value of the deformation gradient
         * \param &deformationGradient_end: The stopping iterator of the computed value of the deformation gradient
         * \param &dFdL_begin: The starting iterator of the derivative of the deformation gradient w.r.t. the velocity
         * gradient \param &dFdL_end: The stopping iterator of the derivative of the deformation gradient w.r.t. the
         * velocity gradient \param &alpha: The integration parameter (0 is explicit and 1 is implicit)
         */

        using deformationGradient_type = typename std::iterator_traits<deformationGradient_iterator>::value_type;

        TARDIGRADE_ERROR_TOOLS_CHECK(
            (unsigned int)(previousDeformationGradient_end - previousDeformationGradient_begin),
            "The previous deformation gradient has a size of " +
                std::to_string((unsigned int)(previousDeformationGradient_end - previousDeformationGradient_begin)) +
                " but must have a size of " + std::to_string(dim * dim))

        TARDIGRADE_ERROR_TOOLS_CHECK((unsigned int)(Lp_end - Lp_begin),
                                     "The previous velocity gradient has a size of " +
                                         std::to_string((unsigned int)(Lp_end - Lp_begin)) +
                                         " but must have a size of " + std::to_string(dim * dim))

        TARDIGRADE_ERROR_TOOLS_CHECK((unsigned int)(L_end - L_begin),
                                     "The velocity gradient has a size of " +
                                         std::to_string((unsigned int)(L_end - L_begin)) + " but must have a size of " +
                                         std::to_string(dim * dim))

        TARDIGRADE_ERROR_TOOLS_CHECK((unsigned int)(deformationGradient_end - deformationGradient_begin),
                                     "The velocity gradient has a size of " +
                                         std::to_string((unsigned int)(L_end - L_begin)) + " but must have a size of " +
                                         std::to_string(dim * dim))

        std::array<deformationGradient_type, dim * dim> DtLalpha, expDtLalpha, tempVector1, tempVector2, tempVector3;

        for (unsigned int i = 0; i < dim * dim; ++i) {
            DtLalpha[i] = Dt * ((1 - alpha) * (*(Lp_begin + i)) + alpha * (*(L_begin + i)));
        }

        std::array<deformationGradient_type, dim * dim * dim * dim> dExpDtLalphadL, tempMatrix1, tempMatrix2;

        TARDIGRADE_ERROR_TOOLS_CATCH(tardigradeVectorTools::computeMatrixExponentialScalingAndSquaring(
            std::begin(DtLalpha), std::end(DtLalpha), dim, std::begin(tempVector1), std::end(tempVector1),
            std::begin(tempVector2), std::end(tempVector2), std::begin(tempVector3), std::end(tempVector3),
            std::begin(tempMatrix1), std::end(tempMatrix1), std::begin(tempMatrix2), std::end(tempMatrix2),
            std::begin(expDtLalpha), std::end(expDtLalpha), std::begin(dExpDtLalphadL), std::end(dExpDtLalphadL)))

        std::transform(std::begin(dExpDtLalphadL), std::end(dExpDtLalphadL), std::begin(dExpDtLalphadL),
                       std::bind(std::multiplies<>(), std::placeholders::_1, Dt * alpha));

        std::fill(deformationGradient_begin, deformationGradient_end, deformationGradient_type());

        std::fill(dFdL_begin, dFdL_end, deformationGradient_type());

        for (unsigned int i = 0; i < dim; ++i) {
            for (unsigned int j = 0; j < dim; ++j) {
                for (unsigned int k = 0; k < dim; ++k) {
                    *(deformationGradient_begin + dim * i + k) +=
                        expDtLalpha[dim * i + j] * (*(previousDeformationGradient_begin + dim * j + k));

                    for (unsigned int ab = 0; ab < dim * dim; ++ab) {
                        *(dFdL_begin + dim * dim * dim * i + dim * dim * k + ab) +=
                            dExpDtLalphadL[dim * dim * dim * i + dim * dim * j + ab] *
                            (*(previousDeformationGradient_begin + dim * j + k));
                    }
                }
            }
        }
    }

    template <unsigned int dim, typename Dt_type, class previousDeformationGradient_iterator, class Lp_iterator,
              class L_iterator, class deformationGradient_iterator, class dFdL_iterator, class dFdFp_iterator,
              class dFdLp_iterator, typename alpha_type>
    void evolveFExponentialMap(const Dt_type                              &Dt,
                               const previousDeformationGradient_iterator &previousDeformationGradient_begin,
                               const previousDeformationGradient_iterator &previousDeformationGradient_end,
                               const Lp_iterator &Lp_begin, const Lp_iterator &Lp_end, const L_iterator &L_begin,
                               const L_iterator &L_end, deformationGradient_iterator deformationGradient_begin,
                               deformationGradient_iterator deformationGradient_end, dFdL_iterator dFdL_begin,
                               dFdL_iterator dFdL_end, dFdFp_iterator dFdFp_begin, dFdFp_iterator dFdFp_end,
                               dFdLp_iterator dFdLp_begin, dFdLp_iterator dFdLp_end, const alpha_type alpha) {
        /*!
         * Evolve the deformation gradient using the exponential map. Assumes the evolution equation is of the form
         *
         * \f$ \dot{F}_{iI} = \ell_{ij} F_{jI} \f$
         *
         * \param &Dt: The change in time
         * \param &previousDeformationGradient_begin: The starting iterator of the previous value of the deformation
         * gradient \param &previousDeformationGradient_end: The stopping iterator of the previous value of the
         * deformation gradient \param &Lp_begin: The starting iterator of the previous value of the velocity gradient
         * \param &Lp_end: The stopping iterator of the previous value of the velocity gradient
         * \param &L_begin: The starting iterator of the current value of the velocity gradient
         * \param &L_end: The stopping iterator of the current value of the velocity gradient
         * \param &deformationGradient_begin: The starting iterator of the computed value of the deformation gradient
         * \param &deformationGradient_end: The stopping iterator of the computed value of the deformation gradient
         * \param &dFdL_begin: The starting iterator of the derivative of the deformation gradient w.r.t. the velocity
         * gradient \param &dFdL_end: The stopping iterator of the derivative of the deformation gradient w.r.t. the
         * velocity gradient \param &dFdFp_begin: The starting iterator of the derivative of the deformation gradient
         * w.r.t. the previous deformation gradient \param &dFdFp_end: The stopping iterator of the derivative of the
         * deformation gradient w.r.t. the previous deformation gradient \param &dFdLp_begin: The starting iterator of
         * the derivative of the deformation gradient w.r.t. the previous velocity gradient \param &dFdLp_end: The
         * stopping iterator of the derivative of the deformation gradient w.r.t. the previous velocity gradient \param
         * &alpha: The integration parameter (0 is explicit and 1 is implicit)
         */

        using deformationGradient_type = typename std::iterator_traits<deformationGradient_iterator>::value_type;

        TARDIGRADE_ERROR_TOOLS_CHECK(
            (unsigned int)(previousDeformationGradient_end - previousDeformationGradient_begin),
            "The previous deformation gradient has a size of " +
                std::to_string((unsigned int)(previousDeformationGradient_end - previousDeformationGradient_begin)) +
                " but must have a size of " + std::to_string(dim * dim))

        TARDIGRADE_ERROR_TOOLS_CHECK((unsigned int)(Lp_end - Lp_begin),
                                     "The previous velocity gradient has a size of " +
                                         std::to_string((unsigned int)(Lp_end - Lp_begin)) +
                                         " but must have a size of " + std::to_string(dim * dim))

        TARDIGRADE_ERROR_TOOLS_CHECK((unsigned int)(L_end - L_begin),
                                     "The velocity gradient has a size of " +
                                         std::to_string((unsigned int)(L_end - L_begin)) + " but must have a size of " +
                                         std::to_string(dim * dim))

        TARDIGRADE_ERROR_TOOLS_CHECK((unsigned int)(deformationGradient_end - deformationGradient_begin),
                                     "The velocity gradient has a size of " +
                                         std::to_string((unsigned int)(L_end - L_begin)) + " but must have a size of " +
                                         std::to_string(dim * dim))

        std::array<deformationGradient_type, dim * dim> DtLalpha, expDtLalpha, tempVector1, tempVector2, tempVector3;

        for (unsigned int i = 0; i < dim * dim; ++i) {
            DtLalpha[i] = Dt * ((1 - alpha) * (*(Lp_begin + i)) + alpha * (*(L_begin + i)));
        }

        std::array<deformationGradient_type, dim * dim * dim * dim> dExpDtLalphadL, dExpDtLalphadLp, tempMatrix1,
            tempMatrix2;

        TARDIGRADE_ERROR_TOOLS_CATCH(tardigradeVectorTools::computeMatrixExponentialScalingAndSquaring(
            std::begin(DtLalpha), std::end(DtLalpha), dim, std::begin(tempVector1), std::end(tempVector1),
            std::begin(tempVector2), std::end(tempVector2), std::begin(tempVector3), std::end(tempVector3),
            std::begin(tempMatrix1), std::end(tempMatrix1), std::begin(tempMatrix2), std::end(tempMatrix2),
            std::begin(expDtLalpha), std::end(expDtLalpha), std::begin(dExpDtLalphadL), std::end(dExpDtLalphadL)))

        std::transform(std::begin(dExpDtLalphadL), std::end(dExpDtLalphadL), std::begin(dExpDtLalphadLp),
                       std::bind(std::multiplies<>(), std::placeholders::_1, Dt * (1 - alpha)));

        std::transform(std::begin(dExpDtLalphadL), std::end(dExpDtLalphadL), std::begin(dExpDtLalphadL),
                       std::bind(std::multiplies<>(), std::placeholders::_1, Dt * alpha));

        std::fill(deformationGradient_begin, deformationGradient_end, deformationGradient_type());

        std::fill(dFdL_begin, dFdL_end, deformationGradient_type());

        std::fill(dFdFp_begin, dFdFp_end, deformationGradient_type());

        std::fill(dFdLp_begin, dFdLp_end, deformationGradient_type());

        for (unsigned int i = 0; i < dim; ++i) {
            for (unsigned int j = 0; j < dim; ++j) {
                for (unsigned int k = 0; k < dim; ++k) {
                    *(deformationGradient_begin + dim * i + k) +=
                        expDtLalpha[dim * i + j] * (*(previousDeformationGradient_begin + dim * j + k));

                    *(dFdFp_begin + dim * dim * dim * i + dim * dim * j + dim * k + j) += expDtLalpha[dim * i + k];

                    for (unsigned int ab = 0; ab < dim * dim; ++ab) {
                        *(dFdL_begin + dim * dim * dim * i + dim * dim * k + ab) +=
                            dExpDtLalphadL[dim * dim * dim * i + dim * dim * j + ab] *
                            (*(previousDeformationGradient_begin + dim * j + k));

                        *(dFdLp_begin + dim * dim * dim * i + dim * dim * k + ab) +=
                            dExpDtLalphadLp[dim * dim * dim * i + dim * dim * j + ab] *
                            (*(previousDeformationGradient_begin + dim * j + k));
                    }
                }
            }
        }
    }

}  // namespace tardigradeConstitutiveTools
