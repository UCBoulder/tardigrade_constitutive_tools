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

}  // namespace tardigradeConstitutiveTools
