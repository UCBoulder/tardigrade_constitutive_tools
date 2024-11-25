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

#include<tardigrade_constitutive_tools.h>

#include<algorithm>

namespace tardigradeConstitutiveTools{

    floatType deltaDirac(const unsigned int i, const unsigned int j){
        /*!
         * The delta dirac function \f$\delta\f$
         *
         * if i==j return 1
         * if i!=j return 0
         *
         * \param i: The first index
         * \param j: The second index
         */

        if (i==j){
            return 1.;
        }
        return 0;
    }

    errorOut rotateMatrix(const floatVector &A, const floatVector &Q, floatVector &rotatedA){
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

        //Check the size of A
        TARDIGRADE_ERROR_TOOLS_CHECK( A.size( ) == Q.size( ), "A and Q must have the same number of values" );

        //Set the dimension to be the square-root of the size of A
        const unsigned int dim = std::sqrt(A.size());
        TARDIGRADE_ERROR_TOOLS_CHECK( ( A.size() % dim ) == 0, "A must be square");

        //Resize rotated A
        rotatedA.resize(A.size());

        for (unsigned int i=0; i<dim; i++){
            for (unsigned int j=0; j<dim; j++){
                for (unsigned int I=0; I<dim; I++){
                    for (unsigned int J=0; J<dim; J++){
                        rotatedA[i*dim + j] += Q[I*dim + i] * A[I*dim + J] * Q[J*dim + j];
                    }
                }
            }
        }

        return NULL;
    }

    void computeDeformationGradient( const floatVector &displacementGradient, floatVector &F, const bool isCurrent ){
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

        const unsigned int dim = ( unsigned int )std::pow( displacementGradient.size( ), 0.5 );
        const unsigned int sot_dim = dim * dim;

        TARDIGRADE_ERROR_TOOLS_CHECK( displacementGradient.size( ) == sot_dim, "The displacement gradienthas " + std::to_string( displacementGradient.size( ) ) + " values but the dimension has been determined to be " + std::to_string( dim ) + "." );

        F = floatVector( sot_dim, 0 );

        std::copy( displacementGradient.begin( ),
                   displacementGradient.end( ),
                   F.begin( ) );

        if ( isCurrent ){

            std::transform( F.cbegin( ), F.cend( ), F.begin( ), std::negate< floatType >( ) );

            for ( unsigned int i = 0; i < dim; i++ ){ F[ dim * i + i ] += 1; }

            Eigen::Map< Eigen::Matrix< floatType, -1, -1 > > map( F.data( ), dim, dim );
            map = map.inverse( ).eval( );

        }
        else{

            for ( unsigned int i = 0; i < dim; i++ ){ F[ dim * i + i ] += 1.; }

        }

        return;

    }

    void computeDeformationGradient( const floatVector &displacementGradient, floatVector &F, floatVector &dFdGradU, const bool isCurrent ){
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

        const unsigned int dim = ( unsigned int )std::pow( displacementGradient.size( ), 0.5 );
        const unsigned int sot_dim = dim * dim;

        TARDIGRADE_ERROR_TOOLS_CHECK( displacementGradient.size( ) == sot_dim, "The displacement gradienthas " + std::to_string( displacementGradient.size( ) ) + " values but the dimension has been determined to be " + std::to_string( dim ) + "." );

        F = floatVector( sot_dim, 0 );

        dFdGradU = floatVector( sot_dim * sot_dim, 0 );

        std::copy( displacementGradient.begin( ),
                   displacementGradient.end( ),
                   F.begin( ) );

        if ( isCurrent ){

            std::transform( F.cbegin( ), F.cend( ), F.begin( ), std::negate< floatType >( ) );

            for ( unsigned int i = 0; i < dim; i++ ){ F[ dim * i + i ] += 1; }

            Eigen::Map< Eigen::Matrix< floatType, -1, -1 > > map( F.data( ), dim, dim );
            map = map.inverse( ).eval( );

            for ( unsigned int i = 0; i < dim; i++ ){

                for ( unsigned int j = 0; j < dim; j++ ){

                    for ( unsigned int k = 0; k < dim; k++ ){

                        for ( unsigned int l = 0; l < dim; l++ ){

                            dFdGradU[ dim * sot_dim * i + sot_dim * j + dim * k + l ]
                                = F[ dim * i + k ] * F[ dim * l + j ];

                        }

                    }

                }

            }

        }
        else{

            for ( unsigned int i = 0; i < dim; i++ ){ F[ dim * i + i ] += 1.; }

            for ( unsigned int i = 0; i < sot_dim; i++ ){ dFdGradU[ sot_dim * i + i ] += 1; }

        }

        return;

    }

    errorOut computeRightCauchyGreen( const floatVector &deformationGradient, floatVector &C ){
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

        //Assume 3D
        constexpr unsigned int dim = 3;
        constexpr unsigned int sot_dim = dim * dim;

        TARDIGRADE_ERROR_TOOLS_CHECK( deformationGradient.size( ) == sot_dim, "The deformation gradient must be 3D" );

        C = floatVector( sot_dim, 0 );

        Eigen::Map< const Eigen::Matrix< floatType, dim, dim, Eigen::RowMajor > > F( deformationGradient.data( ), dim, dim );
        Eigen::Map< Eigen::Matrix< floatType, dim, dim, Eigen::RowMajor > > C_map( C.data( ), dim, dim );

        C_map = ( F.transpose( ) * F ).eval( );

        return NULL;
    }

    errorOut computeRightCauchyGreen( const floatVector &deformationGradient, floatVector &C, floatMatrix &dCdF ){
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

        constexpr unsigned int dim = 3;
        constexpr unsigned int sot_dim = dim * dim;

        floatVector _dCdF;

        TARDIGRADE_ERROR_TOOLS_CATCH( computeRightCauchyGreen( deformationGradient, C, _dCdF ) );

        dCdF = tardigradeVectorTools::inflate( _dCdF, sot_dim, sot_dim );

        return NULL;

    }

    errorOut computeRightCauchyGreen( const floatVector &deformationGradient, floatVector &C, floatVector &dCdF ){
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

        //Assume 3D
        constexpr unsigned int dim = 3;
        constexpr unsigned int sot_dim = dim * dim;

        TARDIGRADE_ERROR_TOOLS_CATCH( computeRightCauchyGreen( deformationGradient, C ) );
        
        //Assemble the Jacobian
        
        dCdF = floatVector( sot_dim * sot_dim, 0 );
        
        for ( unsigned int I = 0; I < dim; I++ ){
            for ( unsigned int J = 0; J < dim; J++ ){
                for ( unsigned int k = 0; k < dim; k++ ){
                    dCdF[ dim * sot_dim * I + sot_dim * J + dim * k + I ] += deformationGradient[ dim * k + J ];
                    dCdF[ dim * sot_dim * I + sot_dim * J + dim * k + J ] += deformationGradient[ dim * k + I ];
                }
            }
        }

        return NULL;

    }

    errorOut computeGreenLagrangeStrain( const floatVector &deformationGradient,
                                         floatVector &E ){
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

        TARDIGRADE_ERROR_TOOLS_CHECK( deformationGradient.size( ) == 9, "The deformation gradient must be 3D." );

        constexpr unsigned int dim=3;
        E.resize( dim * dim );

        for ( unsigned int I = 0; I < dim; I++ ){
            E[ dim * I + I ] -= 1;
            for ( unsigned int J = 0; J < dim; J++ ){
                for ( unsigned int i = 0; i < dim; i++ ){
                    E[ dim * I + J ] += deformationGradient[ dim * i + I ] * deformationGradient[ dim * i + J ];
                }
                E[ dim * I + J ] *= 0.5;
            }
        }
        return NULL;
    }

    errorOut computeGreenLagrangeStrain( const floatVector &deformationGradient, floatVector &E, floatMatrix &dEdF){
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

        TARDIGRADE_ERROR_TOOLS_CATCH( computeGreenLagrangeStrain( deformationGradient, E, _dEdF ) );

        dEdF = tardigradeVectorTools::inflate( _dEdF, deformationGradient.size( ), deformationGradient.size( ) );

        return NULL;

    }

    errorOut computeGreenLagrangeStrain( const floatVector &deformationGradient, floatVector &E, floatVector &dEdF){
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
        
        TARDIGRADE_ERROR_TOOLS_CATCH( computeGreenLagrangeStrain( deformationGradient, E ) );

        TARDIGRADE_ERROR_TOOLS_CATCH( computeDGreenLagrangeStrainDF( deformationGradient, dEdF ) );

        return NULL;
    }

    errorOut computeDGreenLagrangeStrainDF(const floatVector &deformationGradient, floatMatrix &dEdF){
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

        TARDIGRADE_ERROR_TOOLS_CATCH( computeDGreenLagrangeStrainDF( deformationGradient, _dEdF ) );

        dEdF = tardigradeVectorTools::inflate( _dEdF, deformationGradient.size( ), deformationGradient.size( ) );

        return NULL;

    }

    errorOut computeDGreenLagrangeStrainDF(const floatVector &deformationGradient, floatVector &dEdF){
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

        TARDIGRADE_ERROR_TOOLS_CHECK( deformationGradient.size( ) == 9, "the Green-Lagrange strain must be 3D" );

        dEdF = floatVector( 81, 0 );
        for ( unsigned int I = 0; I < 3; I++ ){
            for ( unsigned int J = 0; J < 3; J++ ){
                for ( unsigned int k = 0; k < 3; k++ ){
                    dEdF[ 3 * 9 * I + 9 * J + 3 * k + I ] += 0.5 * deformationGradient[ 3 * k + J ];
                    dEdF[ 3 * 9 * I + 9 * J + 3 * k + J ] += 0.5 * deformationGradient[ 3 * k + I ];
                }
            }
        }

        return NULL;
    }

    errorOut decomposeGreenLagrangeStrain( const floatVector &E, floatVector &Ebar, floatType &J ){
        /*!
         * Decompose the Green-Lagrange strain tensor ( \f$E\f$ ) into isochoric ( \f$\bar{E}\f$ ) and volumetric ( \f$J\f$ ) parts where
         *
         * \f$J = det(F) = sqrt(det(2*E + I))\f$
         *
         * \f$\bar{E}_{IJ} = 0.5*((1/(J**(2/3))) F_{iI} F_{iJ} - I_{IJ}) = (1/(J**(2/3)))*E_{IJ} + 0.5(1/(J**(2/3)) - 1)*I_{IJ}\f$
         *
         * \param &E: The Green-Lagrange strain tensor ( \f$E\f$ )
         * \param &Ebar: The isochoric Green-Lagrange strain tensor ( \f$\bar{E}\f$ ).
         *     format = E11, E12, E13, E21, E22, E23, E31, E32, E33
         * \param &J: The Jacobian of deformation ( \f$J\f$ )
         */

        constexpr unsigned int dim = 3;
        constexpr unsigned int sot_dim = dim * dim;

        TARDIGRADE_ERROR_TOOLS_CHECK( E.size() == sot_dim, "the Green-Lagrange strain must be 3D");

        //Construct the identity tensor
        floatVector F_squared = 2 * E;
        for ( unsigned int i = 0; i < dim; i++ ){ F_squared[ dim * i + i ] += 1; }

        Eigen::Map< Eigen::Matrix< floatType, dim, dim, Eigen::RowMajor > > F_squared_map( F_squared.data( ), dim, dim );

        floatType Jsq = F_squared_map.determinant( );

        TARDIGRADE_ERROR_TOOLS_CHECK( Jsq > 0, "the determinant of the Green-Lagrange strain is negative");

        J = sqrt(Jsq);
        Ebar = E/(pow(J, 2./3));
       
        for ( unsigned int i = 0; i < dim; i++ ){ Ebar[ dim * i + i ] += 0.5*(1/pow(J, 2./3) - 1); }

        return NULL;
    }

    errorOut decomposeGreenLagrangeStrain(const floatVector &E, floatVector &Ebar, floatType &J,
                                          floatVector &dEbardE, floatVector &dJdE){
        /*!
         * Decompose the Green-Lagrange strain tensor ( \f$E\f$ ) into isochoric ( \f$\bar{E}\f$ ) and volumetric ( \f$J\f$ ) parts where
         *
         * \f$J = det(F) = sqrt(det(2*E + I))\f$
         *
         * \f$\bar{E}_{IJ} = 0.5*((1/(J**(2/3))) F_{iI} F_{iJ} - I_{IJ}) = (1/(J**(2/3)))*E_{IJ} + 0.5(1/(J**(2/3)) - 1)*I_{IJ}\f$
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

        constexpr unsigned int dim = 3;
        constexpr unsigned int sot_dim = dim * dim;

        TARDIGRADE_ERROR_TOOLS_CATCH( decomposeGreenLagrangeStrain(E, Ebar, J) );

        //Compute the derivative of the jacobian of deformation w.r.t. the Green-Lagrange strain
        dJdE = 2 * E;
        for ( unsigned int i = 0; i < 3; i++ ){ dJdE[ dim * i + i ] += 1; };

        Eigen::Map< Eigen::Matrix< floatType, dim, dim, Eigen::RowMajor > > dJdE_map( dJdE.data( ), dim, dim );

        dJdE_map = ( J * dJdE_map.inverse( ) ).eval( );

        //Compute the derivative of the isochoric part of the Green-Lagrange strain w.r.t. the Green-Lagrange strain
        floatVector eye( sot_dim );
        for ( unsigned int i = 0; i < dim; i++ ){ eye[ dim * i + i ] = 1.; };
        floatMatrix EYE = tardigradeVectorTools::eye<floatType>(9);

        floatType invJ23 = 1./pow(J, 2./3);
        floatType invJ53 = 1./pow(J, 5./3);

        dEbardE = floatVector( sot_dim * sot_dim, 0 );

        for ( unsigned int i = 0; i < sot_dim; i++ ){
            dEbardE[ sot_dim * i + i ] += invJ23;
        }

        for ( unsigned int i = 0; i < dim; i++ ){

            for ( unsigned int j = 0; j < dim; j++ ){

                for ( unsigned int k = 0; k < dim; k++ ){

                    dEbardE[ sot_dim * dim * i + sot_dim * i + dim * j + k ] -= (1./3) * invJ53 * dJdE[ dim * j + k ];

                    for ( unsigned int l = 0; l < dim; l++ ){

                        dEbardE[ sot_dim * dim * i + sot_dim * j + dim * k + l ] -= (2./3) * invJ53 * E[ dim * i + j ] * dJdE[ dim * k + l ];

                    }

                }

            }

        }

        return NULL;
    }

    errorOut decomposeGreenLagrangeStrain(const floatVector &E, floatVector &Ebar, floatType &J,
                                          floatMatrix &dEbardE, floatVector &dJdE){
        /*!
         * Decompose the Green-Lagrange strain tensor ( \f$E\f$ ) into isochoric ( \f$\bar{E}\f$ ) and volumetric ( \f$J\f$ ) parts where
         *
         * \f$J = det(F) = sqrt(det(2*E + I))\f$
         *
         * \f$\bar{E}_{IJ} = 0.5*((1/(J**(2/3))) F_{iI} F_{iJ} - I_{IJ}) = (1/(J**(2/3)))*E_{IJ} + 0.5(1/(J**(2/3)) - 1)*I_{IJ}\f$
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

        constexpr unsigned int dim = 3;
        constexpr unsigned int sot_dim = dim * dim;

        floatVector _dEbardE;

        TARDIGRADE_ERROR_TOOLS_CATCH( decomposeGreenLagrangeStrain( E, Ebar, J, _dEbardE, dJdE ) );

        dEbardE = tardigradeVectorTools::inflate( _dEbardE, sot_dim, sot_dim );

        return NULL;
    }

    errorOut mapPK2toCauchy(const floatVector &PK2Stress, const floatVector &deformationGradient, floatVector &cauchyStress){
        /*!
         * Map the PK2 stress ( \f$P^{II}\f$ ) to the current configuration resulting in the Cauchy stress ( \f$\sigma\f$ ).
         *
         * \f$\sigma_{ij} = (1/det(F)) F_{iI} P^{II}_{IJ} F_{jJ}\f$
         *
         * where \f$F\f$ is the deformation gradient
         *
         * \param &PK2Stress: The Second Piola-Kirchoff stress ( \f$P^{II}\f$ )
         * \param &deformationGradient: The total deformation gradient ( \f$F\f$ ).
         * \param &cauchyStress: The Cauchy stress (\f$\sigma\f$ ).
         */

        constexpr unsigned int dim = 3;
        constexpr unsigned int sot_dim = dim * dim;

        TARDIGRADE_ERROR_TOOLS_CHECK( PK2Stress.size( ) == sot_dim, "The cauchy stress must have nine components (3D)");

        TARDIGRADE_ERROR_TOOLS_CHECK( deformationGradient.size() == PK2Stress.size(), "The deformation gradient and the PK2 stress don't have the same size");

        //Compute the determinant of the deformation gradient
        Eigen::Map< const Eigen::Matrix< floatType, dim, dim, Eigen::RowMajor > > map( deformationGradient.data( ), dim, dim );
        floatType detF = map.determinant( );

        //Initialize the Cauchy stress
        floatVector temp_sot( sot_dim, 0 );
        cauchyStress = floatVector( sot_dim, 0);

        for (unsigned int i=0; i<dim; i++){
            for (unsigned int I=0; I<dim; I++){
                for (unsigned int j=0; j<dim; j++){
                    temp_sot[ dim * i + j ] += deformationGradient[ 3 * i + I ] * PK2Stress[ 3 * I + j ];
                }
            }
        }

        temp_sot /= detF;

        for ( unsigned int i = 0; i < dim; i++ ){
            for ( unsigned int j = 0; j < dim; j++ ){
                for ( unsigned int I = 0; I < dim; I++ ){
                        cauchyStress[ dim * i + j ] += temp_sot[ dim * i + I ]*deformationGradient[ 3 * j + I ];
                }
            }
        }
        return NULL;
    }

    errorOut WLF(const floatType &temperature, const floatVector &WLFParameters, floatType &factor){
        /*!
         * An implementation of the Williams-Landel-Ferry equation.
         *
         * \f$factor = 10**((-C_1*(T - T_r))/(C_2 + T - T_r))\f$
         *
         * where \f$T\f$ is the temperature, \f$T_r\f$ is the reference temperature, and \f$C_1\f$ and \f$C_2\f$ are parameters
         *
         * \param &temperature: The temperature \f$T\f$
         * \param &WLFParameters: The parameters for the function [\f$T_r\f$, \f$C_1\f$, \f$C_2\f$]
         * \param &factor: The shift factor
         */

        TARDIGRADE_ERROR_TOOLS_CHECK( WLFParameters.size() == 3, "The parameters have the wrong number of terms");

        floatType Tr = WLFParameters[0];
        floatType C1 = WLFParameters[1];
        floatType C2 = WLFParameters[2];

        TARDIGRADE_ERROR_TOOLS_CHECK( !tardigradeVectorTools::fuzzyEquals(C2 + (temperature - Tr), 0.), "Zero in the denominator");

        factor = pow(10., -C1*(temperature - Tr)/(C2 + (temperature - Tr)));

        return NULL;
    }

    errorOut WLF(const floatType &temperature, const floatVector &WLFParameters, floatType &factor, floatType &dfactordT){
        /*!
         * An implementation of the Williams-Landel-Ferry equation that also returns the gradient w.r.t. \f$T\f$
         *
         * \param &temperature: The temperature ( \f$T\f$ )
         * \param &WLFParameters: The parameters for the function [\f$T_r\f$, \f$C_1\f$, \f$C_2\f$]
         * \param &factor: The shift factor
         * \param &dfactordT: The derivative of the shift factor w.r.t. the temperature ( \f$\frac{\partial factor}{\partial T}\f$ )
         */

        TARDIGRADE_ERROR_TOOLS_CATCH( WLF(temperature, WLFParameters, factor) );

        floatType Tr = WLFParameters[0];
        floatType C1 = WLFParameters[1];
        floatType C2 = WLFParameters[2];

        dfactordT = std::log(10)*factor*(-C1/(C2 + temperature - Tr) + (C1*(temperature - Tr)/pow(C2 + temperature - Tr, 2)));

        return NULL;
    }

    errorOut computeDFDt(const floatVector &velocityGradient, const floatVector &deformationGradient, floatVector &DFDt){
        /*!
         * Compute the total time derivative of the deformation gradient.
         *
         * \f$\dot{F}_{iI} = L_{ij} F_{jI}\f$
         *
         * \param &velocityGradient: The velocity gradient \f$L_{ij}\f$
         * \param &deformationGradient: The deformation gradient \f$F_{iI}\f$
         * \param &DFDt: The total time derivative of the deformation gradient
         */

        //Assume 3D
        constexpr unsigned int dim = 3;
        constexpr unsigned int sot_dim = dim * dim;

        TARDIGRADE_ERROR_TOOLS_CHECK( velocityGradient.size() == deformationGradient.size(), "The velocity gradient and deformation gradient must have the same size");

        TARDIGRADE_ERROR_TOOLS_CHECK( velocityGradient.size() == sot_dim, "The velocity gradient doesn't have enough entries");

        DFDt = floatVector(velocityGradient.size(), 0);

        for (unsigned int i=0; i<dim; i++){
            for (unsigned int I=0; I<dim; I++){
                for (unsigned int j=0; j<dim; j++){
                    DFDt[dim*i + I] += velocityGradient[dim*i + j] * deformationGradient[dim*j + I];
                }
            }
        }
        return NULL;
    }

    errorOut computeDFDt(const floatVector &velocityGradient, const floatVector &deformationGradient, floatVector &DFDt,
                         floatVector &dDFDtdL, floatVector &dDFDtdF){
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

        //Assume 3D
        constexpr unsigned int dim = 3;
        constexpr unsigned int sot_dim = dim * dim;

        TARDIGRADE_ERROR_TOOLS_CATCH( computeDFDt(velocityGradient, deformationGradient, DFDt) );

        //Form the identity tensor
        floatVector eye(dim*dim, 0);
        tardigradeVectorTools::eye(eye);

        //Form the partial w.r.t. L and F
        dDFDtdL = floatVector( sot_dim * sot_dim, 0 );
        dDFDtdF = floatVector( sot_dim * sot_dim, 0 );;

        for (unsigned int i=0; i<dim; i++){
            for (unsigned int I=0; I<dim; I++){
                for (unsigned int k=0; k<dim; k++){
                    dDFDtdL[ sot_dim * dim * i + sot_dim * I + dim * i + k] = deformationGradient[ dim * k + I ];
                    dDFDtdF[ sot_dim * dim * i + sot_dim * I + dim * k + I] = velocityGradient[ dim * i + k ];
                }
            }
        }

        return NULL;
    }

    errorOut computeDFDt(const floatVector &velocityGradient, const floatVector &deformationGradient, floatVector &DFDt,
                         floatMatrix &dDFDtdL, floatMatrix &dDFDtdF){
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

        //Assume 3D
        constexpr unsigned int dim = 3;
        constexpr unsigned int sot_dim = dim * dim;

        floatVector _dDFDtdL, _dDFDtdF;

        TARDIGRADE_ERROR_TOOLS_CATCH( computeDFDt( velocityGradient, deformationGradient, DFDt, _dDFDtdL, _dDFDtdF ) );

        dDFDtdL = tardigradeVectorTools::inflate( _dDFDtdL, sot_dim, sot_dim );

        dDFDtdF = tardigradeVectorTools::inflate( _dDFDtdF, sot_dim, sot_dim );

        return NULL;
    }

    errorOut midpointEvolution( const floatType &Dt, const floatVector &Ap, const floatVector &DApDt, const floatVector &DADt,
                                floatVector &dA, floatVector &A, const floatVector &alpha ){
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

        TARDIGRADE_ERROR_TOOLS_CHECK( ( Ap.size( ) == DApDt.size( ) ) && ( Ap.size( ) == DADt.size( ) ), "The size of the previous value of the vector and the two rates are not equal" );

        TARDIGRADE_ERROR_TOOLS_CHECK( Ap.size( ) == alpha.size( ), "The size of the alpha vector is not the same size as the previous vector value" );

        dA = floatVector( Ap.size( ), 0 );

        A = floatVector( Ap.size( ), 0 );

        unsigned int i = 0;

        for ( auto ai = alpha.begin( ); ai != alpha.end( ); ai++, i++ ){

            TARDIGRADE_ERROR_TOOLS_CHECK( ( ( *ai ) >= 0) && ( ( *ai ) <= 1 ), "Alpha must be between 0 and 1" );

            dA[ i ] = Dt * ( *ai * DApDt[ i ] + ( 1 - *ai ) * DADt[ i ] );

            A[ i ]  = Ap[ i ] + dA[ i ];

        }

        return NULL;

    }

    errorOut midpointEvolutionFlatJ( const floatType &Dt, const floatVector &Ap, const floatVector &DApDt, const floatVector &DADt,
                                     floatVector &dA, floatVector &A, floatVector &DADADt, const floatVector &alpha ){
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

        TARDIGRADE_ERROR_TOOLS_CATCH( midpointEvolution( Dt, Ap, DApDt, DADt, dA, A, alpha ) )

        const unsigned int A_size = A.size( );

        DADADt = floatVector( A_size * A_size, 0 );

        unsigned int i = 0;

        for ( auto ai = alpha.begin( ); ai != alpha.end( ); ai++, i++ ){

            DADADt[ A_size * i + i ] = Dt * ( 1 - *ai );

        }

        return NULL;

    }

    errorOut midpointEvolution( const floatType &Dt, const floatVector &Ap, const floatVector &DApDt, const floatVector &DADt,
                                floatVector &dA, floatVector &A, floatMatrix &DADADt, const floatVector &alpha ){
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

        TARDIGRADE_ERROR_TOOLS_CATCH( midpointEvolutionFlatJ( Dt, Ap, DApDt, DADt, dA, A, _DADADt, alpha ) )

        DADADt = tardigradeVectorTools::inflate( _DADADt, A.size( ), A.size( ) );

        return NULL;

    }


    errorOut midpointEvolutionFlatJ( const floatType &Dt, const floatVector &Ap, const floatVector &DApDt, const floatVector &DADt,
                                     floatVector &dA, floatVector &A, floatVector &DADADt, floatVector &DADADtp,
                                     const floatVector &alpha ){
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

        TARDIGRADE_ERROR_TOOLS_CATCH( midpointEvolutionFlatJ( Dt, Ap, DApDt, DADt, dA, A, DADADt, alpha ) );

        const unsigned int A_size = A.size( );

        DADADtp = floatVector( A_size * A_size, 0 );

        unsigned int i = 0;

        for ( auto ai = alpha.begin( ); ai != alpha.end( ); ai++, i++ ){

            DADADtp[ A_size * i + i ] = Dt * ( *ai );

        }

        return NULL;

    }

    errorOut midpointEvolution( const floatType &Dt, const floatVector &Ap, const floatVector &DApDt, const floatVector &DADt,
                                floatVector &dA, floatVector &A, floatMatrix &DADADt, floatMatrix &DADADtp,
                                const floatVector &alpha ){
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

        TARDIGRADE_ERROR_TOOLS_CATCH( midpointEvolutionFlatJ( Dt, Ap, DApDt, DADt, dA, A, _DADADt, _DADADtp, alpha ) );

        DADADt  = tardigradeVectorTools::inflate( _DADADt,  A.size( ), A.size( ) );

        DADADtp = tardigradeVectorTools::inflate( _DADADtp, A.size( ), A.size( ) );

        return NULL;

    }

    errorOut midpointEvolution( const floatType &Dt, const floatVector &Ap, const floatVector &DApDt, const floatVector &DADt,
                                floatVector &dA, floatVector &A, const floatType alpha ){
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

        return midpointEvolution( Dt, Ap, DApDt, DADt, dA, A, alpha * floatVector( Ap.size( ), 1 ) );

    }

    errorOut midpointEvolutionFlatJ( const floatType &Dt, const floatVector &Ap, const floatVector &DApDt, const floatVector &DADt,
                                     floatVector &dA, floatVector &A, floatVector &DADADt, const floatType alpha ){
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

        return midpointEvolutionFlatJ( Dt, Ap, DApDt, DADt, dA, A, DADADt, alpha * floatVector( Ap.size( ), 1 ) );

    }

    errorOut midpointEvolutionFlatJ( const floatType &Dt, const floatVector &Ap, const floatVector &DApDt, const floatVector &DADt,
                                     floatVector &dA, floatVector &A, floatVector &DADADt, floatVector &DADADtp,
                                     const floatType alpha ){
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

        return midpointEvolutionFlatJ( Dt, Ap, DApDt, DADt, dA, A, DADADt, DADADtp, alpha * floatVector( Ap.size( ), 1 ) );

    }

    errorOut midpointEvolution( const floatType &Dt, const floatVector &Ap, const floatVector &DApDt, const floatVector &DADt,
                                floatVector &dA, floatVector &A, floatMatrix &DADADt, const floatType alpha ){
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

        return midpointEvolution( Dt, Ap, DApDt, DADt, dA, A, DADADt, alpha * floatVector( Ap.size( ), 1 ) );

    }

    errorOut midpointEvolution( const floatType &Dt, const floatVector &Ap, const floatVector &DApDt, const floatVector &DADt,
                                floatVector &dA, floatVector &A, floatMatrix &DADADt, floatMatrix &DADADtp,
                                const floatType alpha ){
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

        return midpointEvolution( Dt, Ap, DApDt, DADt, dA, A, DADADt, DADADtp, alpha * floatVector( Ap.size( ), 1 ) );

    }

    errorOut evolveF(const floatType &Dt, const floatVector &previousDeformationGradient, const floatVector &Lp, const floatVector &L,
                     floatVector &dF, floatVector &deformationGradient, const floatType alpha, const unsigned int mode){
        /*!
         * Evolve the deformation gradient ( F ) using the midpoint integration method.
         *
         * mode 1:
         * \f$F_{iI}^{t + 1} = \left[\delta_{ij} - \Delta t \left(1 - \alpha \right) L_{ij}^{t+1} \right]^{-1} \left[F_{iI}^{t} + \Delta t \alpha \dot{F}_{iI}^{t} \right]\f$
         *
         * mode 2:
         * \f$F_{iI}^{t + 1} = \left[F_{iJ}^{t} + \Delta t \alpha \dot{F}_{iJ}^{t} \right] \left[\delta_{IJ} - \Delta T \left( 1- \alpha \right) L_{IJ}^{t+1} \right]^{-1}\f$
         *
         * \param &Dt: The change in time.
         * \param &previousDeformationGradient: The previous value of the deformation gradient
         * \param &Lp: The previous velocity gradient in the current configuration (mode 1) or
         *     reference configuration (mode 2).
         * \param &L: The current velocity gradient in the current configuration (mode 1) or
         *     reference configuration (mode 2).
         * \param &dF: The change in the deformation gradient \f$\Delta \bf{F}\f$ such that \f$F_{iI}^{t+1} = F_{iI}^t + \Delta F_{iI}\f$
         * \param &deformationGradient: The computed current deformation gradient.
         * \param alpha: The integration parameter.
         * \param mode: The mode of the ODE. Whether the velocity gradient is known in the
         *     current (mode 1) or reference (mode 2) configuration.
         */

        //Assumes 3D
        constexpr unsigned int dim = 3;
        constexpr unsigned int sot_dim = dim * dim;

        TARDIGRADE_ERROR_TOOLS_CHECK( previousDeformationGradient.size( ) == sot_dim, "The deformation gradient doesn't have enough terms (require 9 for 3D)" );

        TARDIGRADE_ERROR_TOOLS_CHECK( Lp.size( ) == previousDeformationGradient.size( ), "The previous velocity gradient and deformation gradient aren't the same size" );

        TARDIGRADE_ERROR_TOOLS_CHECK( previousDeformationGradient.size( ) == L.size( ), "The previous deformation gradient and the current velocity gradient aren't the same size" );

        TARDIGRADE_ERROR_TOOLS_CHECK( ( mode == 1 ) || ( mode == 2 ), "The mode of evolution is not recognized" );

        //Compute L^{t + \alpha}
        floatVector LtpAlpha = alpha * Lp + ( 1 - alpha ) * L;

        //Compute the right hand side

        floatVector RHS( sot_dim, 0 );
        Eigen::Map< const Eigen::Matrix< floatType, dim, dim, Eigen::RowMajor > > Fp( previousDeformationGradient.data( ), dim, dim );
        Eigen::Map< const Eigen::Matrix< floatType, dim, dim, Eigen::RowMajor > > Lt( LtpAlpha.data( ), dim, dim );
        Eigen::Map< Eigen::Matrix< floatType, dim, dim, Eigen::RowMajor > > RHS_map( RHS.data( ), dim, dim );

        if ( mode == 1 ){
            RHS_map = ( Lt * Fp ).eval( );
        }
        if ( mode == 2 ){
            RHS_map = ( Fp * Lt ).eval( );
        }

        RHS *= Dt;

        //Compute the left-hand side
        floatVector eye( dim * dim );

        tardigradeVectorTools::eye( eye );

        floatVector invLHS = -Dt * ( 1 - alpha ) * L;
        for ( unsigned int i = 0; i < dim; i++ ){ invLHS[ dim * i + i ] += 1; }

        dF = floatVector( sot_dim, 0 );

        Eigen::Map< Eigen::Matrix< floatType, dim, dim, Eigen::RowMajor > > invLHS_map( invLHS.data( ), dim, dim );
        Eigen::Map< Eigen::Matrix< floatType, dim, dim, Eigen::RowMajor > > dF_map( dF.data( ), dim, dim );

        invLHS_map = invLHS_map.inverse( ).eval( );

        if ( mode == 1 ){

            dF_map = ( invLHS_map * RHS_map ).eval( );

        }
        if ( mode == 2 ){

            dF_map = ( RHS_map * invLHS_map ).eval( );

        }

        deformationGradient = previousDeformationGradient + dF;

        return NULL;

    }

    errorOut evolveF(const floatType &Dt, const floatVector &previousDeformationGradient, const floatVector &Lp, const floatVector &L,
                     floatVector &deformationGradient, const floatType alpha, const unsigned int mode){
        /*!
         * Evolve the deformation gradient ( F ) using the midpoint integration method.
         *
         * mode 1:
         * \f$F_{iI}^{t + 1} = \left[\delta_{ij} - \Delta t \left(1 - \alpha \right) L_{ij}^{t+1} \right]^{-1} \left[F_{iI}^{t} + \Delta t \alpha \dot{F}_{iI}^{t} \right]\f$
         *
         * mode 2:
         * \f$F_{iI}^{t + 1} = \left[F_{iJ}^{t} + \Delta t \alpha \dot{F}_{iJ}^{t} \right] \left[\delta_{IJ} - \Delta T \left( 1- \alpha \right) L_{IJ}^{t+1} \right]^{-1}\f$
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

        floatVector dF;

        return evolveF( Dt, previousDeformationGradient, Lp, L, dF, deformationGradient, alpha, mode );

    }

    errorOut evolveF( const floatType &Dt, const floatVector &previousDeformationGradient, const floatVector &Lp, const floatVector &L,
                      floatVector &dF, floatVector &deformationGradient, floatMatrix &dFdL, const floatType alpha, const unsigned int mode ){
        /*!
         * Evolve the deformation gradient ( F ) using the midpoint integration method and return the jacobian w.r.t. L.
         *
         * mode 1:
         * \f$F_{iI}^{t + 1} = \left[\delta_{ij} - \Delta t \left(1 - \alpha \right) L_{ij}^{t+1} \right]^{-1} \left[F_{iI}^{t} + \Delta t \alpha \dot{F}_{iI}^{t} \right]\f$
         * \f$\frac{\partial F_{jI}^{t + 1}}{\partial L_{kl}^{t+1}} = \left[\delta_{kj} - \Delta t \left(1 - \alpha\right) L_{kj}\right]^{-1} \Delta t \left(1 - \alpha\right) F_{lI}^{t + 1}\f$
         *
         * mode 2:
         * \f$F_{iI}^{t + 1} = \left[F_{iJ}^{t} + \Delta t \alpha \dot{F}_{iJ}^{t} \right] \left[\delta_{IJ} - \Delta T \left( 1- \alpha \right) L_{IJ}^{t+1} \right]^{-1}\f$
         * \f$\frac{\partial F_{iJ}^{t + 1}}{\partial L_{KL}} = \Delta t (1 - \alpha) F_{iK}^{t + 1} \left[\delta_{JL} - \right ]\f$
         *
         * \param &Dt: The change in time.
         * \param &previousDeformationGradient: The previous value of the deformation gradient
         * \param &Lp: The previous velocity gradient.
         * \param &L: The current velocity gradient.
         * \param &dF: The change in the deformation gradient \f$\Delta \bf{F}\f$ such that \f$F_{iI}^{t+1} = F_{iI}^t + \Delta F_{iI}\f$
         * \param &deformationGradient: The computed current deformation gradient.
         * \param &dFdL: The derivative of the deformation gradient w.r.t. the velocity gradient
         * \param alpha: The integration parameter ( 0 for implicit, 1 for explicit )
         * \param mode: The form of the ODE. See above for details.
         */

        constexpr unsigned int dim = 3;
        const unsigned int sot_dim = dim * dim;

        floatVector _dFdL;

        TARDIGRADE_ERROR_TOOLS_CATCH( evolveFFlatJ( Dt, previousDeformationGradient, Lp, L, dF, deformationGradient, _dFdL, alpha, mode ) );

        dFdL = tardigradeVectorTools::inflate( _dFdL, sot_dim, sot_dim );

        return NULL;

    }

    errorOut evolveFFlatJ( const floatType &Dt, const floatVector &previousDeformationGradient, const floatVector &Lp, const floatVector &L,
                           floatVector &dF, floatVector &deformationGradient, floatVector &dFdL, const floatType alpha, const unsigned int mode ){
        /*!
         * Evolve the deformation gradient ( F ) using the midpoint integration method and return the jacobian w.r.t. L.
         *
         * mode 1:
         * \f$F_{iI}^{t + 1} = \left[\delta_{ij} - \Delta t \left(1 - \alpha \right) L_{ij}^{t+1} \right]^{-1} \left[F_{iI}^{t} + \Delta t \alpha \dot{F}_{iI}^{t} \right]\f$
         * \f$\frac{\partial F_{jI}^{t + 1}}{\partial L_{kl}^{t+1}} = \left[\delta_{kj} - \Delta t \left(1 - \alpha\right) L_{kj}\right]^{-1} \Delta t \left(1 - \alpha\right) F_{lI}^{t + 1}\f$
         *
         * mode 2:
         * \f$F_{iI}^{t + 1} = \left[F_{iJ}^{t} + \Delta t \alpha \dot{F}_{iJ}^{t} \right] \left[\delta_{IJ} - \Delta T \left( 1- \alpha \right) L_{IJ}^{t+1} \right]^{-1}\f$
         * \f$\frac{\partial F_{iJ}^{t + 1}}{\partial L_{KL}} = \Delta t (1 - \alpha) F_{iK}^{t + 1} \left[\delta_{JL} - \right ]\f$
         *
         * \param &Dt: The change in time.
         * \param &previousDeformationGradient: The previous value of the deformation gradient
         * \param &Lp: The previous velocity gradient.
         * \param &L: The current velocity gradient.
         * \param &dF: The change in the deformation gradient \f$\Delta \bf{F}\f$ such that \f$F_{iI}^{t+1} = F_{iI}^t + \Delta F_{iI}\f$
         * \param &deformationGradient: The computed current deformation gradient.
         * \param &dFdL: The derivative of the deformation gradient w.r.t. the velocity gradient
         * \param alpha: The integration parameter ( 0 for implicit, 1 for explicit )
         * \param mode: The form of the ODE. See above for details.
         */

        //Assumes 3D
        constexpr unsigned int dim = 3;
        const unsigned int sot_dim = dim * dim;

        TARDIGRADE_ERROR_TOOLS_CATCH( evolveF( Dt, previousDeformationGradient, Lp, L, dF, deformationGradient, alpha, mode) );

        //Compute the left hand side
        floatVector invLHS = -Dt * ( 1 - alpha ) * L;
        for ( unsigned int i = 0; i < dim; i++ ){ invLHS[ dim * i + i ] += 1; }

        Eigen::Map< Eigen::Matrix< floatType, dim, dim, Eigen::RowMajor > > invLHS_map( invLHS.data( ), dim, dim );

        invLHS_map = invLHS_map.inverse( ).eval( );

        //Compute the jacobian
        dFdL = floatVector( sot_dim * sot_dim, 0 );
        if ( mode == 1 ){
            for ( unsigned int j = 0; j < dim; j++ ){
                for ( unsigned int I = 0; I < dim; I++ ){
                    for ( unsigned int k = 0; k < dim; k++ ){
                        for ( unsigned int l = 0; l < dim; l++ ){
                            dFdL[ dim * sot_dim * j + sot_dim * I + dim * k + l ] += Dt * ( 1 - alpha ) * invLHS[ dim * j + k ] * deformationGradient[ dim * l + I ];
                        }
                    }
                }
            }
        }
        if ( mode == 2 ){
            for ( unsigned int j = 0; j < dim; j++ ){
                for ( unsigned int I = 0; I < dim; I++ ){
                    for ( unsigned int K = 0; K < dim; K++ ){
                        for ( unsigned int _L = 0; _L < dim; _L++ ){
                            dFdL[ dim * sot_dim * j + sot_dim * I + dim * K + _L ] += Dt * ( 1 - alpha ) * invLHS[ dim * _L + I ] * deformationGradient[ dim * j + K ];
                        }
                    }
                }
            }
        }
        return NULL;
    }

    errorOut evolveFFlatJ( const floatType &Dt, const floatVector &previousDeformationGradient, const floatVector &Lp, const floatVector &L,
                           floatVector &deformationGradient, floatVector &dFdL, const floatType alpha, const unsigned int mode ){
        /*!
         * Evolve the deformation gradient ( F ) using the midpoint integration method and return the jacobian w.r.t. L.
         *
         * mode 1:
         * \f$F_{iI}^{t + 1} = \left[\delta_{ij} - \Delta t \left(1 - \alpha \right) L_{ij}^{t+1} \right]^{-1} \left[F_{iI}^{t} + \Delta t \alpha \dot{F}_{iI}^{t} \right]\f$
         * \f$\frac{\partial F_{jI}^{t + 1}}{\partial L_{kl}^{t+1}} = \left[\delta_{kj} - \Delta t \left(1 - \alpha\right) L_{kj}\right]^{-1} \Delta t \left(1 - \alpha\right) F_{lI}^{t + 1}\f$
         *
         * mode 2:
         * \f$F_{iI}^{t + 1} = \left[F_{iJ}^{t} + \Delta t \alpha \dot{F}_{iJ}^{t} \right] \left[\delta_{IJ} - \Delta T \left( 1- \alpha \right) L_{IJ}^{t+1} \right]^{-1}\f$
         * \f$\frac{\partial F_{iJ}^{t + 1}}{\partial L_{KL}} = \Delta t (1 - \alpha) F_{iK}^{t + 1} \left[\delta_{JL} - \right ]\f$
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

        return evolveFFlatJ( Dt, previousDeformationGradient, Lp, L, dF, deformationGradient, dFdL, alpha, mode );

    }

    errorOut evolveF( const floatType &Dt, const floatVector &previousDeformationGradient, const floatVector &Lp, const floatVector &L,
                      floatVector &deformationGradient, floatMatrix &dFdL, const floatType alpha, const unsigned int mode ){
        /*!
         * Evolve the deformation gradient ( F ) using the midpoint integration method and return the jacobian w.r.t. L.
         *
         * mode 1:
         * \f$F_{iI}^{t + 1} = \left[\delta_{ij} - \Delta t \left(1 - \alpha \right) L_{ij}^{t+1} \right]^{-1} \left[F_{iI}^{t} + \Delta t \alpha \dot{F}_{iI}^{t} \right]\f$
         * \f$\frac{\partial F_{jI}^{t + 1}}{\partial L_{kl}^{t+1}} = \left[\delta_{kj} - \Delta t \left(1 - \alpha\right) L_{kj}\right]^{-1} \Delta t \left(1 - \alpha\right) F_{lI}^{t + 1}\f$
         *
         * mode 2:
         * \f$F_{iI}^{t + 1} = \left[F_{iJ}^{t} + \Delta t \alpha \dot{F}_{iJ}^{t} \right] \left[\delta_{IJ} - \Delta T \left( 1- \alpha \right) L_{IJ}^{t+1} \right]^{-1}\f$
         * \f$\frac{\partial F_{iJ}^{t + 1}}{\partial L_{KL}} = \Delta t (1 - \alpha) F_{iK}^{t + 1} \left[\delta_{JL} - \right ]\f$
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

        return evolveF( Dt, previousDeformationGradient, Lp, L, dF, deformationGradient, dFdL, alpha, mode );

    }

    errorOut evolveF( const floatType &Dt, const floatVector &previousDeformationGradient, const floatVector &Lp, const floatVector &L,
                      floatVector &dF, floatVector &deformationGradient, floatMatrix &dFdL, floatMatrix &ddFdFp, floatMatrix &dFdFp, floatMatrix &dFdLp, const floatType alpha, const unsigned int mode ){
        /*!
         * Evolve the deformation gradient ( F ) using the midpoint integration method and return the jacobian w.r.t. L.
         *
         * mode 1:
         * \f$F_{iI}^{t + 1} = \left[\delta_{ij} - \Delta t \left(1 - \alpha \right) L_{ij}^{t+1} \right]^{-1} \left[F_{iI}^{t} + \Delta t \alpha \dot{F}_{iI}^{t} \right]\f$
         * \f$\frac{\partial F_{jI}^{t + 1}}{\partial L_{kl}^{t+1}} = \left[\delta_{kj} - \Delta t \left(1 - \alpha\right) L_{kj}\right]^{-1} \Delta t \left(1 - \alpha\right) F_{lI}^{t + 1}\f$
         *
         * mode 2:
         * \f$F_{iI}^{t + 1} = \left[F_{iJ}^{t} + \Delta t \alpha \dot{F}_{iJ}^{t} \right] \left[\delta_{IJ} - \Delta T \left( 1- \alpha \right) L_{IJ}^{t+1} \right]^{-1}\f$
         * \f$\frac{\partial F_{iJ}^{t + 1}}{\partial L_{KL}} = \Delta t (1 - \alpha) F_{iK}^{t + 1} \left[\delta_{JL} - \right ]\f$
         *
         * \param &Dt: The change in time.
         * \param &previousDeformationGradient: The previous value of the deformation gradient
         * \param &Lp: The previous velocity gradient.
         * \param &L: The current velocity gradient.
         * \param &dF: The change in the deformation gradient \f$\Delta \bf{F}\f$ such that \f$F_{iI}^{t+1} = F_{iI}^t + \Delta F_{iI}\f$
         * \param &deformationGradient: The computed current deformation gradient.
         * \param &dFdL: The derivative of the deformation gradient w.r.t. the velocity gradient
         * \param &ddFdFp: The derivative of the change in the deformation gradient w.r.t. the previous deformation gradient
         * \param &dFdFp: The derivative of the deformation gradient w.r.t. the previous deformation gradient
         * \param &dFdLp: The derivative of the deformation gradient w.r.t. the previous velocity gradient
         * \param alpha: The integration parameter ( 0 for implicit, 1 for explicit )
         * \param mode: The form of the ODE. See above for details.
         */

        //Assumes 3D
        constexpr unsigned int dim = 3;
        const unsigned int sot_dim = dim * dim;

        floatVector _dFdL;
        floatVector _ddFdFp;
        floatVector _dFdFp;
        floatVector _dFdLp;

        TARDIGRADE_ERROR_TOOLS_CATCH( evolveFFlatJ( Dt, previousDeformationGradient, Lp, L,
                                                    dF, deformationGradient, _dFdL, _ddFdFp, _dFdFp, _dFdLp, alpha, mode ) );

        dFdL   = tardigradeVectorTools::inflate( _dFdL,   sot_dim, sot_dim );
        ddFdFp = tardigradeVectorTools::inflate( _ddFdFp, sot_dim, sot_dim );
        dFdFp  = tardigradeVectorTools::inflate( _dFdFp,  sot_dim, sot_dim );
        dFdLp  = tardigradeVectorTools::inflate( _dFdLp,  sot_dim, sot_dim );

        return NULL;

    }
    errorOut evolveFFlatJ( const floatType &Dt, const floatVector &previousDeformationGradient, const floatVector &Lp, const floatVector &L,
                           floatVector &dF, floatVector &deformationGradient, floatVector &dFdL, floatVector &ddFdFp, floatVector &dFdFp, floatVector &dFdLp, const floatType alpha, const unsigned int mode ){
        /*!
         * Evolve the deformation gradient ( F ) using the midpoint integration method and return the jacobian w.r.t. L.
         *
         * mode 1:
         * \f$F_{iI}^{t + 1} = \left[\delta_{ij} - \Delta t \left(1 - \alpha \right) L_{ij}^{t+1} \right]^{-1} \left[F_{iI}^{t} + \Delta t \alpha \dot{F}_{iI}^{t} \right]\f$
         * \f$\frac{\partial F_{jI}^{t + 1}}{\partial L_{kl}^{t+1}} = \left[\delta_{kj} - \Delta t \left(1 - \alpha\right) L_{kj}\right]^{-1} \Delta t \left(1 - \alpha\right) F_{lI}^{t + 1}\f$
         *
         * mode 2:
         * \f$F_{iI}^{t + 1} = \left[F_{iJ}^{t} + \Delta t \alpha \dot{F}_{iJ}^{t} \right] \left[\delta_{IJ} - \Delta T \left( 1- \alpha \right) L_{IJ}^{t+1} \right]^{-1}\f$
         * \f$\frac{\partial F_{iJ}^{t + 1}}{\partial L_{KL}} = \Delta t (1 - \alpha) F_{iK}^{t + 1} \left[\delta_{JL} - \right ]\f$
         *
         * \param &Dt: The change in time.
         * \param &previousDeformationGradient: The previous value of the deformation gradient
         * \param &Lp: The previous velocity gradient.
         * \param &L: The current velocity gradient.
         * \param &dF: The change in the deformation gradient \f$\Delta \bf{F}\f$ such that \f$F_{iI}^{t+1} = F_{iI}^t + \Delta F_{iI}\f$
         * \param &deformationGradient: The computed current deformation gradient.
         * \param &dFdL: The derivative of the deformation gradient w.r.t. the velocity gradient
         * \param &ddFdFp: The derivative of the change in the deformation gradient w.r.t. the previous deformation gradient
         * \param &dFdFp: The derivative of the deformation gradient w.r.t. the previous deformation gradient
         * \param &dFdLp: The derivative of the deformation gradient w.r.t. the previous velocity gradient
         * \param alpha: The integration parameter ( 0 for implicit, 1 for explicit )
         * \param mode: The form of the ODE. See above for details.
         */

        //Assumes 3D
        constexpr unsigned int dim = 3;
        const unsigned int sot_dim = dim * dim;

        TARDIGRADE_ERROR_TOOLS_CATCH( evolveF( Dt, previousDeformationGradient, Lp, L, dF, deformationGradient, alpha, mode) );

        //Compute L^{t + \alpha}
        floatVector LtpAlpha = alpha * Lp + ( 1 - alpha ) * L;

        //Compute the left hand side
        floatVector eye( dim * dim );

        tardigradeVectorTools::eye( eye );

        floatVector invLHS = -Dt * ( 1 - alpha ) * L;
        for ( unsigned int i = 0; i < dim; i++ ){ invLHS[ dim * i + i ] += 1; }

        Eigen::Map< Eigen::Matrix< floatType, dim, dim, Eigen::RowMajor > > invLHS_map( invLHS.data( ), dim, dim );

        invLHS_map = invLHS_map.inverse( ).eval( );

        //Compute the jacobian
        dFdL   = floatVector( sot_dim * sot_dim, 0 );
        ddFdFp = floatVector( sot_dim * sot_dim, 0 );
        dFdLp  = floatVector( sot_dim * sot_dim, 0 );
        if ( mode == 1 ){
            for ( unsigned int j = 0; j < dim; j++ ){
                for ( unsigned int I = 0; I < dim; I++ ){
                    for ( unsigned int k = 0; k < dim; k++ ){
                        for ( unsigned int l = 0; l < dim; l++ ){
                            dFdL[ dim * sot_dim * j + sot_dim * I + dim * k + l ]  += Dt * ( 1 - alpha ) * invLHS[ dim * j + k ] * deformationGradient[ dim * l + I ];
                            dFdLp[ dim * sot_dim * j + sot_dim * I + dim * k + l ] += invLHS[ dim * j + k ] * Dt * alpha * previousDeformationGradient[ dim * l + I ];
                            ddFdFp[ dim * sot_dim * j + sot_dim * I + dim * k + I ] += Dt * invLHS[ dim * j + l ] * LtpAlpha[ dim * l + k ];
                        }
                    }
                }
            }
        }
        if ( mode == 2 ){
            for ( unsigned int j = 0; j < dim; j++ ){
                for ( unsigned int I = 0; I < dim; I++ ){
                    for ( unsigned int K = 0; K < dim; K++ ){
                        for ( unsigned int _L = 0; _L < dim; _L++ ){
                            dFdL[ dim * sot_dim * j + sot_dim * I + dim * K + _L ]  += Dt * ( 1 - alpha ) * invLHS[ dim * _L + I ] * deformationGradient[ dim * j + K ];
                            dFdLp[ dim * sot_dim * j + sot_dim * I + dim * K + _L ] += Dt * alpha * invLHS[ dim * _L + I ] * previousDeformationGradient[ dim * j + K ];
                            ddFdFp[ dim * sot_dim * j + sot_dim * I + dim * j + K ] += Dt * LtpAlpha[ dim * K + _L ] * invLHS[ dim * _L + I ];
                        }
                    }
                }
            }
        }

        dFdFp = ddFdFp;
        for ( unsigned int i = 0; i < sot_dim; i++ ){ dFdFp[ sot_dim * i + i ] += 1; }

        return NULL;
    }

    errorOut evolveFFlatJ( const floatType &Dt, const floatVector &previousDeformationGradient, const floatVector &Lp, const floatVector &L,
                           floatVector &deformationGradient, floatVector &dFdL, floatVector &dFdFp, floatVector &dFdLp, const floatType alpha, const unsigned int mode ){
        /*!
         * Evolve the deformation gradient ( F ) using the midpoint integration method and return the jacobian w.r.t. L.
         *
         * mode 1:
         * \f$F_{iI}^{t + 1} = \left[\delta_{ij} - \Delta t \left(1 - \alpha \right) L_{ij}^{t+1} \right]^{-1} \left[F_{iI}^{t} + \Delta t \alpha \dot{F}_{iI}^{t} \right]\f$
         * \f$\frac{\partial F_{jI}^{t + 1}}{\partial L_{kl}^{t+1}} = \left[\delta_{kj} - \Delta t \left(1 - \alpha\right) L_{kj}\right]^{-1} \Delta t \left(1 - \alpha\right) F_{lI}^{t + 1}\f$
         *
         * mode 2:
         * \f$F_{iI}^{t + 1} = \left[F_{iJ}^{t} + \Delta t \alpha \dot{F}_{iJ}^{t} \right] \left[\delta_{IJ} - \Delta T \left( 1- \alpha \right) L_{IJ}^{t+1} \right]^{-1}\f$
         * \f$\frac{\partial F_{iJ}^{t + 1}}{\partial L_{KL}} = \Delta t (1 - \alpha) F_{iK}^{t + 1} \left[\delta_{JL} - \right ]\f$
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

        return evolveFFlatJ( Dt, previousDeformationGradient, Lp, L, dF, deformationGradient, dFdL, ddFdFp, dFdFp, dFdLp, alpha, mode );

    }

    errorOut evolveF( const floatType &Dt, const floatVector &previousDeformationGradient, const floatVector &Lp, const floatVector &L,
                      floatVector &deformationGradient, floatMatrix &dFdL, floatMatrix &dFdFp, floatMatrix &dFdLp, const floatType alpha, const unsigned int mode ){
        /*!
         * Evolve the deformation gradient ( F ) using the midpoint integration method and return the jacobian w.r.t. L.
         *
         * mode 1:
         * \f$F_{iI}^{t + 1} = \left[\delta_{ij} - \Delta t \left(1 - \alpha \right) L_{ij}^{t+1} \right]^{-1} \left[F_{iI}^{t} + \Delta t \alpha \dot{F}_{iI}^{t} \right]\f$
         * \f$\frac{\partial F_{jI}^{t + 1}}{\partial L_{kl}^{t+1}} = \left[\delta_{kj} - \Delta t \left(1 - \alpha\right) L_{kj}\right]^{-1} \Delta t \left(1 - \alpha\right) F_{lI}^{t + 1}\f$
         *
         * mode 2:
         * \f$F_{iI}^{t + 1} = \left[F_{iJ}^{t} + \Delta t \alpha \dot{F}_{iJ}^{t} \right] \left[\delta_{IJ} - \Delta T \left( 1- \alpha \right) L_{IJ}^{t+1} \right]^{-1}\f$
         * \f$\frac{\partial F_{iJ}^{t + 1}}{\partial L_{KL}} = \Delta t (1 - \alpha) F_{iK}^{t + 1} \left[\delta_{JL} - \right ]\f$
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

        return evolveF( Dt, previousDeformationGradient, Lp, L, dF, deformationGradient, dFdL, ddFdFp, dFdFp, dFdLp, alpha, mode );

    }

    floatType mac(const floatType &x){
        /*!
         * Compute the Macaulay brackets of a scalar x
         *
         * returns x if x>0, 0 otherwise
         *
         * \param &x: The incoming scalar.
         */

        return 0.5 * ( fabs( x ) + x );
    }

    floatType mac( const floatType &x, floatType &dmacdx ){
        /*!
         * Compute the Macaulay brackets of the scalar x and
         * return the jacobian as well.
         * 
         * returns x if x>0, 0 otherwise
         * 
         * The Jacobian is the Heaviside function
         * 
         * \param &x: The incoming scalar
         * \param &dmacdx: The returned jacobian
         */

        dmacdx = 0;
        if ( x >= 0 ){ dmacdx = 1; }
        return mac( x );
    }

    errorOut computeUnitNormal(const floatVector &A, floatVector &Anorm){
        /*!
         * Compute the unit normal of a second order tensor (or strictly speaking
         * any tensor).
         *
         * \param &A: The second order tensor
         * \param &Anorm: The unit normal in the direction of A
         */

        const unsigned int A_size = A.size( );

        floatType norm = sqrt(tardigradeVectorTools::inner(A, A));

        if ( tardigradeVectorTools::fuzzyEquals( norm, 0. ) ){
            Anorm = floatVector( A_size, 0 );
        }
        else {
            Anorm = A/norm;
        }

        return NULL;
    }

    errorOut computeUnitNormal(const floatVector &A, floatVector &Anorm, floatVector &dAnormdA){
        /*!
         * Compute the unit normal of a second order tensor (or strictly speaking any
         * tensor) and the gradient of that unit normal w.r.t. the tensor.
         *
         * \param &A: The second order tensor
         * \param &Anorm: The unit normal in the direction of A
         * \param &dAnormdA: The gradient of the unit normal w.r.t. A
         */

        const unsigned int A_size = A.size( );

        floatType norm = sqrt(tardigradeVectorTools::inner(A, A));

        if ( tardigradeVectorTools::fuzzyEquals( norm, 0. ) ){
            Anorm = floatVector( A_size, 0 );
        }
        else {
            Anorm = A/norm;
        }

        dAnormdA = floatVector( A_size * A_size, 0 );

        for ( unsigned int i = 0; i < A_size; i++ ){

            dAnormdA[ A_size * i + i ] += 1;

            for ( unsigned int j = 0; j < A_size; j++ ){

                dAnormdA[ A_size * i + j ] -= Anorm[ i ] * Anorm[ j ];

            }

        }

        dAnormdA /= norm;

        return NULL;
    }

    errorOut computeUnitNormal(const floatVector &A, floatVector &Anorm, floatMatrix &dAnormdA){
        /*!
         * Compute the unit normal of a second order tensor (or strictly speaking any
         * tensor) and the gradient of that unit normal w.r.t. the tensor.
         *
         * \param &A: The second order tensor
         * \param &Anorm: The unit normal in the direction of A
         * \param &dAnormdA: The gradient of the unit normal w.r.t. A
         */

        const unsigned int A_size = A.size( );

        floatVector _dAnormdA;

        TARDIGRADE_ERROR_TOOLS_CATCH( computeUnitNormal( A, Anorm, _dAnormdA ) );

        dAnormdA = tardigradeVectorTools::inflate( _dAnormdA, A_size, A_size );

        return NULL;
    }

    errorOut pullBackVelocityGradient(const floatVector &velocityGradient, const floatVector &deformationGradient,
                                      floatVector &pulledBackVelocityGradient){
        /*!
         * Pull back the velocity gradient to the configuration indicated by deformationGradient, i.e.
         *
         * \f$totalDeformationGradient_{iI} = deformationGradient_{i \bar{I}} remainingDeformationGradient_{\bar{I}I}\f$
         *
         * This is done via
         *
         * \f$L_{\bar{I} \bar{J}} = deformationGradient_{\bar{I} i}^{-1} velocityGradient_{ij} deformationGradient_{j\bar{J}}\f$
         *
         * \param &velocityGradient: The velocity gradient in the current configuration.
         * \param &deformationGradient: The deformation gradient between the desired configuration
         *     and the current configuration.
         * \param &pulledBackVelocityGradient: The pulled back velocity gradient.
         */

        //Assume 3D
        constexpr unsigned int dim = 3;
        constexpr unsigned int sot_dim = dim * dim;

        //Invert the deformation gradient
        pulledBackVelocityGradient = floatVector( sot_dim, 0 );
        Eigen::Map< const Eigen::Matrix< floatType, dim, dim, Eigen::RowMajor > > F( deformationGradient.data( ), dim, dim );
        Eigen::Map< const Eigen::Matrix< floatType, dim, dim, Eigen::RowMajor > > L( velocityGradient.data( ), dim, dim );
        Eigen::Map< Eigen::Matrix< floatType, dim, dim, Eigen::RowMajor > > pullBackL( pulledBackVelocityGradient.data( ), dim, dim );

        pullBackL = ( F.inverse( ) * L * F ).eval( );

        return NULL;
    }

    errorOut pullBackVelocityGradient(const floatVector &velocityGradient, const floatVector &deformationGradient,
                                      floatVector &pulledBackVelocityGradient, floatVector &dPullBackLdL,
                                      floatVector &dPullBackLdF){
        /*!
         * Pull back the velocity gradient to the configuration indicated by deformationGradient, i.e.
         *
         * \f$totalDeformationGradient_{iI} = deformationGradient_{i \bar{I}} remainingDeformationGradient_{\bar{I}I}\f$
         *
         * This is done via
         *
         * \f$L_{\bar{I} \bar{J}} = deformationGradient_{\bar{I} i}^{-1} velocityGradient_{ij} deformationGradient_{j\bar{J}}\f$
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

        //Assume 3D
        constexpr unsigned int dim = 3;
        constexpr unsigned int sot_dim = dim * dim;

        pulledBackVelocityGradient = floatVector( sot_dim, 0 );
        Eigen::Map< const Eigen::Matrix< floatType, dim, dim, Eigen::RowMajor > > F( deformationGradient.data( ), dim, dim );
        Eigen::Map< const Eigen::Matrix< floatType, dim, dim, Eigen::RowMajor > > L( velocityGradient.data( ), dim, dim );
        Eigen::Map< Eigen::Matrix< floatType, dim, dim, Eigen::RowMajor > > pullBackL( pulledBackVelocityGradient.data( ), dim, dim );

        floatVector inverseDeformationGradient = deformationGradient;
        Eigen::Map< Eigen::Matrix< floatType, dim, dim, Eigen::RowMajor > > invF_map( inverseDeformationGradient.data( ), dim, dim );
        invF_map = invF_map.inverse( ).eval( );

        floatVector term2( sot_dim, 0 );
        Eigen::Map< Eigen::Matrix< floatType, dim, dim, Eigen::RowMajor > > term2_map( term2.data( ), dim, dim );

        term2_map = ( invF_map * L ).eval( );

        //Pull back the velocity gradient
        pullBackL = ( term2_map * F ).eval( );

        //Construct the gradients
        dPullBackLdL = floatVector( sot_dim * sot_dim, 0 );
        dPullBackLdF = floatVector( sot_dim * sot_dim, 0 );

        for (unsigned int I=0; I<dim; I++){
            for (unsigned int J=0; J<dim; J++){
                for (unsigned int k=0; k<dim; k++){
                    for (unsigned int l=0; l<dim; l++){
                        dPullBackLdL[sot_dim * dim * I + sot_dim * J + dim * k + l ] = inverseDeformationGradient[dim*I + k] * deformationGradient[dim*l + J];
                    }

                    dPullBackLdF[sot_dim * dim * I + sot_dim * J + dim * k + J ] += term2[dim*I + k];

                    for ( unsigned int K = 0; K < dim; K++ ){
                        dPullBackLdF[sot_dim * dim * I + sot_dim * J + dim * k + K ] -= inverseDeformationGradient[dim*I + k] * pulledBackVelocityGradient[dim*K + J];
                    }
                }
            }
        }

        return NULL;
    }

    errorOut pullBackVelocityGradient(const floatVector &velocityGradient, const floatVector &deformationGradient,
                                      floatVector &pulledBackVelocityGradient, floatMatrix &dPullBackLdL,
                                      floatMatrix &dPullBackLdF){
        /*!
         * Pull back the velocity gradient to the configuration indicated by deformationGradient, i.e.
         *
         * \f$totalDeformationGradient_{iI} = deformationGradient_{i \bar{I}} remainingDeformationGradient_{\bar{I}I}\f$
         *
         * This is done via
         *
         * \f$L_{\bar{I} \bar{J}} = deformationGradient_{\bar{I} i}^{-1} velocityGradient_{ij} deformationGradient_{j\bar{J}}\f$
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

        //Assume 3D
        constexpr unsigned int dim = 3;
        constexpr unsigned int sot_dim = dim * dim;

        floatVector _dPullBackLdL, _dPullBackLdF;

        TARDIGRADE_ERROR_TOOLS_CATCH( pullBackVelocityGradient( velocityGradient, deformationGradient, pulledBackVelocityGradient, _dPullBackLdL, _dPullBackLdF ) );

        dPullBackLdL = tardigradeVectorTools::inflate( _dPullBackLdL, sot_dim, sot_dim );

        dPullBackLdF = tardigradeVectorTools::inflate( _dPullBackLdF, sot_dim, sot_dim );

        return NULL;
    }

    errorOut quadraticThermalExpansion(const floatType &temperature, const floatType &referenceTemperature,
                                       const floatVector &linearParameters, const floatVector &quadraticParameters,
                                       floatVector &thermalExpansion){
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

        TARDIGRADE_ERROR_TOOLS_CHECK( linearParameters.size() == quadraticParameters.size(), "The linear and quadratic parameters must have the same length");

        thermalExpansion = linearParameters * temperature          + quadraticParameters * temperature * temperature
                         - linearParameters * referenceTemperature - quadraticParameters * referenceTemperature * referenceTemperature;

        return NULL;
    }

    errorOut quadraticThermalExpansion(const floatType &temperature, const floatType &referenceTemperature,
                                       const floatVector &linearParameters, const floatVector &quadraticParameters,
                                       floatVector &thermalExpansion, floatVector &thermalExpansionJacobian){
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

        TARDIGRADE_ERROR_TOOLS_CATCH( quadraticThermalExpansion(temperature, referenceTemperature, linearParameters, quadraticParameters,
                                                                thermalExpansion) )

        thermalExpansionJacobian = linearParameters + 2 * quadraticParameters * temperature;

        return NULL;
    }

    errorOut pushForwardGreenLagrangeStrain(const floatVector &greenLagrangeStrain, const floatVector &deformationGradient,
                                            floatVector &almansiStrain){
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

        //Assume 3D
        constexpr unsigned int dim = 3;
        constexpr unsigned int sot_dim = dim * dim;

        floatVector inverseDeformationGradient = deformationGradient;
        Eigen::Map< Eigen::Matrix< floatType, dim, dim, Eigen::RowMajor > > invF( inverseDeformationGradient.data( ), dim, dim );
        invF = invF.inverse( ).eval( );

        almansiStrain = floatVector( sot_dim, 0 );

        Eigen::Map< const Eigen::Matrix< floatType, dim, dim, Eigen::RowMajor > > E( greenLagrangeStrain.data( ), dim, dim );
        Eigen::Map< Eigen::Matrix< floatType, dim, dim, Eigen::RowMajor > > e( almansiStrain.data( ), dim, dim );

        //Map the Green-Lagrange strain to the current configuration
        e = ( invF.transpose( ) * E * invF ).eval( );

        return NULL;
    }

    errorOut pushForwardGreenLagrangeStrain(const floatVector &greenLagrangeStrain, const floatVector &deformationGradient,
                                            floatVector &almansiStrain, floatVector &dAlmansiStraindE, floatVector &dAlmansiStraindF){
        /*!
         * Push forward the Green-Lagrange strain to the current configuration
         * and return the jacobians.
         *
         * \f$e_{ij} = F_{Ii}^{-1} E_{IJ} F_{Jj}^{-1}\f$
         *
         * \f$\frac{\partial e_{ij}}{\partial E_{KL}} = F_{Ki}^{-1} F_{Kj}^{-1}\f$
         *
         * \f$\frac{\partial e_{ij}}{\partial F_{kK}} = -F_{Ik}^{-1} F_{Ki}^{-1} E_{IJ} F_{J j}^{-1} - F_{Ii}^{-1} E_{IJ} F_{Jk}^{-1} F_{Kj}^{-1}\f$
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

        //Assume 3D
        constexpr unsigned int dim = 3;
        constexpr unsigned int sot_dim = dim * dim;

        floatVector inverseDeformationGradient = deformationGradient;
        Eigen::Map< Eigen::Matrix< floatType, dim, dim, Eigen::RowMajor > > invF( inverseDeformationGradient.data( ), dim, dim );
        invF = invF.inverse( ).eval( );

        almansiStrain = floatVector( sot_dim, 0 );

        Eigen::Map< const Eigen::Matrix< floatType, dim, dim, Eigen::RowMajor > > E( greenLagrangeStrain.data( ), dim, dim );
        Eigen::Map< Eigen::Matrix< floatType, dim, dim, Eigen::RowMajor > > e( almansiStrain.data( ), dim, dim );

        //Map the Green-Lagrange strain to the current configuration
        e = ( invF.transpose( ) * E * invF ).eval( );

        //Compute the jacobians
        dAlmansiStraindE = floatVector( sot_dim * sot_dim, 0 );
        dAlmansiStraindF = floatVector( sot_dim * sot_dim, 0 );
        for (unsigned int i=0; i<dim; i++){
            for (unsigned int j=0; j<dim; j++){
                for (unsigned int K=0; K<dim; K++){
                    for (unsigned int L=0; L<dim; L++){
                        dAlmansiStraindE[sot_dim * dim * i + sot_dim * j + dim * K + L ] = inverseDeformationGradient[dim*K + i] *
                                                                                           inverseDeformationGradient[dim*L + j];
                        dAlmansiStraindF[sot_dim * dim * i + sot_dim * j + dim * K + L ] = -inverseDeformationGradient[dim*L + i ] * almansiStrain[dim*K+j]
                                                                                           -inverseDeformationGradient[dim*L + j ] * almansiStrain[dim*i+K];
                    }
                }
            }
        }

        return NULL;
    }

    errorOut pushForwardGreenLagrangeStrain(const floatVector &greenLagrangeStrain, const floatVector &deformationGradient,
                                            floatVector &almansiStrain, floatMatrix &dAlmansiStraindE, floatMatrix &dAlmansiStraindF){
        /*!
         * Push forward the Green-Lagrange strain to the current configuration
         * and return the jacobians.
         *
         * \f$e_{ij} = F_{Ii}^{-1} E_{IJ} F_{Jj}^{-1}\f$
         *
         * \f$\frac{\partial e_{ij}}{\partial E_{KL}} = F_{Ki}^{-1} F_{Kj}^{-1}\f$
         *
         * \f$\frac{\partial e_{ij}}{\partial F_{kK}} = -F_{Ik}^{-1} F_{Ki}^{-1} E_{IJ} F_{J j}^{-1} - F_{Ii}^{-1} E_{IJ} F_{Jk}^{-1} F_{Kj}^{-1}\f$
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

        //Assume 3D
        constexpr unsigned int dim = 3;
        constexpr unsigned int sot_dim = dim * dim;

        floatVector _dAlmansiStraindE, _dAlmansiStraindF;

        TARDIGRADE_ERROR_TOOLS_CATCH( pushForwardGreenLagrangeStrain( greenLagrangeStrain, deformationGradient, almansiStrain, _dAlmansiStraindE, _dAlmansiStraindF ) );

        dAlmansiStraindE = tardigradeVectorTools::inflate( _dAlmansiStraindE, sot_dim, sot_dim );
        dAlmansiStraindF = tardigradeVectorTools::inflate( _dAlmansiStraindF, sot_dim, sot_dim );

        return NULL;

    } 

    errorOut pullBackAlmansiStrain( const floatVector &almansiStrain, const floatVector &deformationGradient,
                                    floatVector &greenLagrangeStrain ){
        /*!
         * Pull back the almansi strain to the configuration indicated by the deformation gradient.
         *
         * \param &almansiStrain: The strain in the deformation gradient's current configuration.
         * \param &deformationGradient: The deformation gradient between configurations.
         * \param &greenLagrangeStrain: The Green-Lagrange strain which corresponds to the reference
         *     configuration of the deformation gradient.
         */

        //Assume 3d
        constexpr unsigned int dim = 3;
        constexpr unsigned int sot_dim = dim * dim;

        Eigen::Map< const Eigen::Matrix< floatType, dim, dim, Eigen::RowMajor > > F( deformationGradient.data( ), dim, dim );

        greenLagrangeStrain = floatVector( sot_dim, 0 );

        Eigen::Map< Eigen::Matrix< floatType, dim, dim, Eigen::RowMajor > > E( greenLagrangeStrain.data( ), dim, dim );
        Eigen::Map< const Eigen::Matrix< floatType, dim, dim, Eigen::RowMajor > > e( almansiStrain.data( ), dim, dim );

        E = ( F.transpose( ) * e * F ).eval( );

        return NULL;
    }

    errorOut pullBackAlmansiStrain( const floatVector &almansiStrain, const floatVector &deformationGradient,
                                    floatVector &greenLagrangeStrain, floatVector &dEde, floatVector &dEdF ){
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

        //Assume 3d
        constexpr unsigned int dim = 3;
        constexpr unsigned int sot_dim = dim * dim;

        TARDIGRADE_ERROR_TOOLS_CATCH( pullBackAlmansiStrain( almansiStrain, deformationGradient, greenLagrangeStrain ) )

        floatVector eye( dim * dim );
        tardigradeVectorTools::eye( eye );

        dEde = floatVector( sot_dim * sot_dim, 0 );
        dEdF = floatVector( sot_dim * sot_dim, 0 );

        for ( unsigned int I = 0; I < dim; I++ ){
            for ( unsigned int J = 0; J < dim; J++ ){
                for ( unsigned int K = 0; K < dim; K++ ){
                    for ( unsigned int L = 0; L < dim; L++ ){
                        dEde[ sot_dim * dim * I + sot_dim * J + dim * K + L ] = deformationGradient[ dim * K + I ] * deformationGradient[ dim * L + J ];
                        dEdF[ sot_dim * dim * I + sot_dim * J + dim * K + I ] += almansiStrain[ dim * K + L ] * deformationGradient[ dim * L + J ];
                        dEdF[ sot_dim * dim * I + sot_dim * J + dim * K + J ] += deformationGradient[ dim * L + I ] * almansiStrain[ dim * L + K ];
                    }
                }
            }
        }

        return NULL;
    }

    errorOut pullBackAlmansiStrain( const floatVector &almansiStrain, const floatVector &deformationGradient,
                                    floatVector &greenLagrangeStrain, floatMatrix &dEde, floatMatrix &dEdF ){
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

        //Assume 3d
        constexpr unsigned int dim = 3;
        constexpr unsigned int sot_dim = dim * dim;

        floatVector _dEde, _dEdF;

        TARDIGRADE_ERROR_TOOLS_CATCH( pullBackAlmansiStrain( almansiStrain, deformationGradient, greenLagrangeStrain, _dEde, _dEdF ) )

        dEde = tardigradeVectorTools::inflate( _dEde, sot_dim, sot_dim );
        dEdF = tardigradeVectorTools::inflate( _dEdF, sot_dim, sot_dim );

        return NULL;
    }

    errorOut computeSymmetricPart( const floatVector &A, floatVector &symmA, unsigned int &dim ){
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
        
        //Get the dimension of A
        dim = ( unsigned int )( std::sqrt( ( double )A.size( ) ) + 0.5 );
        const unsigned int sot_dim = dim * dim;
   
        TARDIGRADE_ERROR_TOOLS_CHECK( sot_dim == A.size( ), "A is not a square matrix" );
 
        symmA = floatVector( A.size( ), 0 );

        for ( unsigned int i = 0; i < dim; i++ ){
            for ( unsigned int j = 0; j < dim; j++ ){
                symmA[ dim * i + j ] = 0.5 * ( A[ dim * i + j ] + A[ dim * j + i ] );
            }
        }
    
        return NULL;
    }

    errorOut computeSymmetricPart( const floatVector &A, floatVector &symmA ){
        /*!
         * Compute the symmetric part of a second order tensor ( \f$A\f$ ) and return it.
         *
         * \f$symm( A )_ij = \frac{1}{2}\left(A_{ij} + A_{ji}\right)\f$
         *
         * \param &A: A constant reference to the second order tensor to process ( \f$A\f$ )
         * \param &symmA: The symmetric part of A ( \f$A^{symm}\f$ )
         */
    
        unsigned int dim;
        return computeSymmetricPart( A, symmA, dim );
    }

    errorOut computeSymmetricPart( const floatVector &A, floatVector &symmA, floatMatrix &dSymmAdA ){
        /*!
         * Compute the symmetric part of a second order tensor ( \f$A\f$ ) and return it.
         *
         * \f$( A )^{symm}_{ij} = \frac{1}{2}\left(A_{ij} + A_{ji}\right)\f$
         *
         * Also computes the jacobian
         * 
         * \f$\frac{\partial A^{symm}_{ij}}{\partial A_{kl}} = \frac{1}{2}\left( \delta_{ik} \delta_{jl} + \delta_{jk}\delta_{il} \right) \f$
         *
         * \param &A: A constant reference to the second order tensor to process ( \f$A\f$ )
         * \param &symmA: The symmetric part of A ( \f$A^{symm}\f$ )
         * \param &dSymmAdA: The Jacobian of the symmetric part of A w.r.t. A ( \f$\frac{\partial A^{symm}}{\partial A}\f$ )
         */
     
        unsigned int dim;
        TARDIGRADE_ERROR_TOOLS_CATCH( computeSymmetricPart( A, symmA, dim ) );
        
        dSymmAdA = floatMatrix( symmA.size( ), floatVector( A.size( ), 0 ) );
        
        for ( unsigned int i = 0; i < dim; i++ ){
            for ( unsigned int j = 0; j < dim; j++ ){
                dSymmAdA[ dim * i + j ][ dim * i + j ] += 0.5;
                dSymmAdA[ dim * i + j ][ dim * j + i ] += 0.5;
            }
        }
        
        return NULL;
   
    }

    errorOut computeSymmetricPart( const floatVector &A, floatVector &symmA, floatVector &dSymmAdA ){
        /*!
         * Compute the symmetric part of a second order tensor ( \f$A\f$ ) and return it.
         *
         * \f$( A )^{symm}_{ij} = \frac{1}{2}\left(A_{ij} + A_{ji}\right)\f$
         *
         * Also computes the jacobian
         * 
         * \f$\frac{\partial A^{symm}_{ij}}{\partial A_{kl}} = \frac{1}{2}\left( \delta_{ik} \delta_{jl} + \delta_{jk}\delta_{il} \right)\f$
         *
         * \param &A: A constant reference to the second order tensor to process ( \f$A\f$ )
         * \param &symmA: The symmetric part of A ( \f$A^{symm}\f$ )
         * \param &dSymmAdA: The Jacobian of the symmetric part of A w.r.t. A ( \f$\frac{\partial A^{symm}}{\partial A}\f$ )
         */
        
        unsigned int dim;
        TARDIGRADE_ERROR_TOOLS_CATCH( computeSymmetricPart( A, symmA, dim ) );
        
        const unsigned int Asize = A.size( );
        
        dSymmAdA = floatVector( symmA.size( ) * Asize, 0 );

        for ( unsigned int i = 0; i < dim; i++ ){
            for ( unsigned int j = 0; j < dim; j++ ){
                dSymmAdA[ dim * Asize * i + Asize * j + dim * i + j ] += 0.5;
                dSymmAdA[ dim * Asize * i + Asize * j + dim * j + i ] += 0.5;
            }
        }
        
        return NULL;
    }

    errorOut pushForwardPK2Stress( const floatVector &PK2, const floatVector &F, floatVector &cauchyStress ){
        /*!
         * Push the Second Piola-Kirchhoff stress forward to the current configuration resulting in the Cauchy stress
         * 
         * \f$ \sigma_{ij} = \frac{1}{J} F_{iI} S_{IJ} F_{jJ} \f$
         * 
         * \param &PK2: The Second Piola-Kirchhoff stress \f$ S_{IJ} \f$
         * \param &F: The deformation gradient \f$ F_{iI} \f$
         * \param &cauchyStress: The Cauchy stress \f$ \sigma_{ij} \f$
         */

        constexpr unsigned int dim = 3;
        constexpr unsigned int sot_dim = dim * dim;
        
        TARDIGRADE_ERROR_TOOLS_CHECK( sot_dim == PK2.size( ), "The PK2 stress must have a size of " + std::to_string( sot_dim ) + " and has a size of " + std::to_string( PK2.size( ) ) )

        TARDIGRADE_ERROR_TOOLS_CHECK( PK2.size( ) == F.size( ), "The deformation gradient must have a size of " + std::to_string( PK2.size( ) ) + " and has a size of " + std::to_string( F.size( ) ) );

        cauchyStress = floatVector( sot_dim, 0 );

        Eigen::Map< const Eigen::Matrix< floatType, dim, dim, Eigen::RowMajor > > PK2_map( PK2.data( ), dim, dim );
        Eigen::Map< const Eigen::Matrix< floatType, dim, dim, Eigen::RowMajor > > F_map( F.data( ), dim, dim );
        Eigen::Map< Eigen::Matrix< floatType, dim, dim, Eigen::RowMajor > > cauchyStress_map( cauchyStress.data( ), dim, dim );

        floatType J = F_map.determinant( );

        cauchyStress_map = ( F_map * PK2_map * F_map.transpose( ) / J ).eval( );

        return NULL;

    }

    errorOut pushForwardPK2Stress( const floatVector &PK2, const floatVector &F, floatVector &cauchyStress,
                                   floatVector &dCauchyStressdPK2, floatVector &dCauchyStressdF ){
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

        constexpr unsigned int dim = 3;
        constexpr unsigned int sot_dim = dim * dim;
        constexpr unsigned int fot_dim = sot_dim * sot_dim;
        
        TARDIGRADE_ERROR_TOOLS_CHECK( dim * dim == PK2.size( ), "The PK2 stress must have a size of " + std::to_string( sot_dim ) + " and has a size of " + std::to_string( PK2.size( ) ) );

        TARDIGRADE_ERROR_TOOLS_CHECK( PK2.size( ) == F.size( ), "The deformation gradient must have a size of " + std::to_string( PK2.size( ) ) + " and has a size of " + std::to_string( F.size( ) ) );

        cauchyStress = floatVector( sot_dim, 0 );
        floatVector dJdF( sot_dim, 0 );

        Eigen::Map< const Eigen::Matrix< floatType, dim, dim, Eigen::RowMajor > > PK2_map( PK2.data( ), dim, dim );
        Eigen::Map< const Eigen::Matrix< floatType, dim, dim, Eigen::RowMajor > > F_map( F.data( ), dim, dim );
        Eigen::Map< Eigen::Matrix< floatType, dim, dim, Eigen::RowMajor > > cauchyStress_map( cauchyStress.data( ), dim, dim );
        Eigen::Map< Eigen::Matrix< floatType, dim, dim, Eigen::RowMajor > > dJdF_map( dJdF.data( ), dim, dim );

        floatType J = F_map.determinant( );

        cauchyStress_map = ( F_map * PK2_map * F_map.transpose( ) / J ).eval( );

        dJdF_map = ( J * F_map.inverse( ).transpose( ) ).eval( );

        dCauchyStressdF = floatVector( fot_dim, 0 );

        dCauchyStressdPK2 = floatVector( fot_dim, 0 );

        floatVector eye( dim * dim, 0 );
        tardigradeVectorTools::eye( eye );

        for ( unsigned int i = 0; i < dim; i++ ){

            for ( unsigned int j = 0; j < dim; j++ ){

                for ( unsigned int A = 0; A < dim; A++ ){

                    for ( unsigned int B = 0; B < dim; B++ ){

                        dCauchyStressdPK2[ dim * sot_dim * i + sot_dim * j + dim * A + B ] += F[ dim * i + A ] * F[ dim * j + B ] / J;
                        dCauchyStressdF[ dim * sot_dim * i + sot_dim * j + dim * A + B ] -= cauchyStress[ dim * i + j ] * dJdF[ dim * A + B ] / J;

                        dCauchyStressdF[ dim * sot_dim * i + sot_dim * j + dim * i + A ] += PK2[ dim * A + B ] * F[ dim * j + B ] / J;

                        dCauchyStressdF[ dim * sot_dim * i + sot_dim * j + dim * j + A ] += F[ dim * i + B ] * PK2[ dim * B + A ] / J;

                    }

                }

            }

        } 

        return NULL;

    }

    errorOut pushForwardPK2Stress( const floatVector &PK2, const floatVector &F, floatVector &cauchyStress,
                                   floatMatrix &dCauchyStressdPK2, floatMatrix &dCauchyStressdF ){
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

        TARDIGRADE_ERROR_TOOLS_CATCH( pushForwardPK2Stress( PK2, F, cauchyStress, _dCauchyStressdPK2, _dCauchyStressdF ) )

        dCauchyStressdPK2 = tardigradeVectorTools::inflate( _dCauchyStressdPK2, cauchyStress.size( ), PK2.size( ) );

        dCauchyStressdF   = tardigradeVectorTools::inflate( _dCauchyStressdF,   cauchyStress.size( ), F.size( ) );

        return NULL;

    }

    errorOut pullBackCauchyStress( const floatVector &cauchyStress, const floatVector &F, floatVector &PK2 ){
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

        constexpr unsigned int dim = 3;
        constexpr unsigned int sot_dim = dim * dim;

        TARDIGRADE_ERROR_TOOLS_CHECK( cauchyStress.size( ) == sot_dim, "The Cauchy stress size is not consistent with the computed dimension\n    cauchyStress.size( ): " + std::to_string( cauchyStress.size( ) ) + "\n    dim * dim           : " + std::to_string( dim * dim ) + "\n" );

        TARDIGRADE_ERROR_TOOLS_CHECK( cauchyStress.size( ) == F.size( ), "The Cauchy stress and the deformation gradient have inconsistent sizes\n    cauchyStress.size( ): " + std::to_string( cauchyStress.size( ) ) + "\n    F.size( )           : " + std::to_string( F.size( ) ) + "\n" );

        Eigen::Map< const Eigen::Matrix< floatType, dim, dim, Eigen::RowMajor > > cauchyStress_map( cauchyStress.data( ), dim, dim );

        floatVector Finv = F;

        Eigen::Map< Eigen::Matrix< floatType, dim, dim, Eigen::RowMajor > > Finv_map( Finv.data( ), dim, dim );
        Finv_map = Finv_map.inverse( ).eval( );

        floatType J = 1 / Finv_map.determinant( );

        PK2 = floatVector( sot_dim, 0 );
        Eigen::Map< Eigen::Matrix< floatType, dim, dim, Eigen::RowMajor > > PK2_map( PK2.data( ), dim, dim );

        PK2_map = ( J * Finv_map * cauchyStress_map * Finv_map.transpose( ) ).eval( );

        return NULL;

    }
    
    errorOut pullBackCauchyStress( const floatVector &cauchyStress, const floatVector &F, floatVector &PK2, 
                                   floatVector &dPK2dCauchyStress, floatVector &dPK2dF ){
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

        constexpr unsigned int dim = 3;
        constexpr unsigned int sot_dim = dim * dim;
        constexpr unsigned int fot_dim = sot_dim * sot_dim;

        TARDIGRADE_ERROR_TOOLS_CHECK( cauchyStress.size( ) == sot_dim, "The Cauchy stress size is not consistent with the computed dimension\n    cauchyStress.size( ): " + std::to_string( cauchyStress.size( ) ) + "\n    dim * dim           : " + std::to_string( dim * dim ) + "\n" );

        TARDIGRADE_ERROR_TOOLS_CHECK( cauchyStress.size( ) == F.size( ), "The Cauchy stress and the deformation gradient have inconsistent sizes\n    cauchyStress.size( ): " + std::to_string( cauchyStress.size( ) ) + "\n    F.size( )           : " + std::to_string( F.size( ) ) + "\n" );

        Eigen::Map< const Eigen::Matrix< floatType, dim, dim, Eigen::RowMajor > > cauchyStress_map( cauchyStress.data( ), dim, dim );

        floatVector Finv = F;

        Eigen::Map< Eigen::Matrix< floatType, dim, dim, Eigen::RowMajor > > Finv_map( Finv.data( ), dim, dim );
        Finv_map = Finv_map.inverse( ).eval( );

        floatType J = 1 / Finv_map.determinant( );

        PK2 = floatVector( sot_dim, 0 );
        Eigen::Map< Eigen::Matrix< floatType, dim, dim, Eigen::RowMajor > > PK2_map( PK2.data( ), dim, dim );

        PK2_map = ( J * Finv_map * cauchyStress_map * Finv_map.transpose( ) ).eval( );

        dPK2dCauchyStress = floatVector( fot_dim, 0 );
        dPK2dF            = floatVector( fot_dim, 0 );

        for ( unsigned int A = 0; A < dim; A++ ){

            for ( unsigned int B = 0; B < dim; B++ ){

                for ( unsigned int k = 0; k < dim; k++ ){

                    for ( unsigned int l = 0; l < dim; l++ ){

                        dPK2dCauchyStress[ dim * dim * dim * A + dim * dim * B + dim * k + l ] = J * Finv[ dim * A + k ] * Finv[ dim * B + l ];

                        dPK2dF[ dim * dim * dim * A + dim * dim * B + dim * k + l ] = Finv[ dim * l + k ] * PK2[ dim * A + B ]
                                                                                    - Finv[ dim * A + k ] * PK2[ dim * l + B ]
                                                                                    - Finv[ dim * B + k ] * PK2[ dim * A + l ];

                    }

                }

            }

        }

        return NULL;

    }

    errorOut pullBackCauchyStress( const floatVector &cauchyStress, const floatVector &F, floatVector &PK2, 
                                   floatMatrix &dPK2dCauchyStress, floatMatrix &dPK2dF ){
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

        TARDIGRADE_ERROR_TOOLS_CATCH( pullBackCauchyStress( cauchyStress, F, PK2, _dPK2dCauchyStress, _dPK2dF ) )

        dPK2dCauchyStress = tardigradeVectorTools::inflate( _dPK2dCauchyStress, PK2.size( ), cauchyStress.size( ) );

        dPK2dF = tardigradeVectorTools::inflate( _dPK2dF, PK2.size( ), F.size( ) );

        return NULL;

    }

    void evolveFExponentialMap( const floatType &Dt, const floatVector &previousDeformationGradient, const floatVector &Lp, const floatVector &L,
                                floatVector &deformationGradient, const floatType alpha ){
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

        constexpr unsigned int dim = 3;
        constexpr unsigned int sot_dim = dim * dim;

        floatVector DtLalpha = Dt * ( ( 1 - alpha ) * Lp + alpha * L );

        floatVector expDtLalpha;

        TARDIGRADE_ERROR_TOOLS_CATCH( tardigradeVectorTools::computeMatrixExponentialScalingAndSquaring( DtLalpha, dim, expDtLalpha ) )

        deformationGradient = floatVector( sot_dim, 0 );

        for ( unsigned int i = 0; i < dim; i++ ){

            for ( unsigned int j = 0; j < dim; j++ ){

                for ( unsigned int k = 0; k < dim; k++ ){

                    deformationGradient[ dim * i + k ] += expDtLalpha[ dim * i + j ] * previousDeformationGradient[ dim * j + k ];

                }

            }

        }

    }

    void evolveFExponentialMap( const floatType &Dt, const floatVector &previousDeformationGradient, const floatVector &Lp, const floatVector &L,
                                floatVector &deformationGradient, floatVector &dFdL, const floatType alpha ){
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

        constexpr unsigned int dim = 3;
        constexpr unsigned int sot_dim = dim * dim;

        floatVector DtLalpha = Dt * ( ( 1 - alpha ) * Lp + alpha * L );

        floatVector expDtLalpha;

        floatVector dExpDtLalphadL;

        TARDIGRADE_ERROR_TOOLS_CATCH( tardigradeVectorTools::computeMatrixExponentialScalingAndSquaring( DtLalpha, dim, expDtLalpha, dExpDtLalphadL ) )

        dExpDtLalphadL *= Dt * alpha;

        deformationGradient = floatVector( sot_dim, 0 );

        dFdL = floatVector( sot_dim * sot_dim, 0 );

        for ( unsigned int i = 0; i < dim; i++ ){

            for ( unsigned int j = 0; j < dim; j++ ){

                for ( unsigned int k = 0; k < dim; k++ ){

                    deformationGradient[ dim * i + k ] += expDtLalpha[ dim * i + j ] * previousDeformationGradient[ dim * j + k ];

                    for ( unsigned int ab = 0; ab < sot_dim; ab++ ){

                        dFdL[ dim * sot_dim * i + sot_dim * k + ab ] += dExpDtLalphadL[ dim * sot_dim * i + sot_dim * j + ab ] * previousDeformationGradient[ dim * j + k ];

                    }

                }

            }

        }

    }

    void evolveFExponentialMap( const floatType &Dt, const floatVector &previousDeformationGradient, const floatVector &Lp, const floatVector &L,
                                floatVector &deformationGradient, floatVector &dFdL, floatVector &dFdFp, floatVector &dFdLp, const floatType alpha ){
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

        constexpr unsigned int dim = 3;
        constexpr unsigned int sot_dim = dim * dim;

        floatVector DtLalpha = Dt * ( ( 1 - alpha ) * Lp + alpha * L );

        floatVector expDtLalpha;

        floatVector dExpDtLalphadL;

        TARDIGRADE_ERROR_TOOLS_CATCH( tardigradeVectorTools::computeMatrixExponentialScalingAndSquaring( DtLalpha, dim, expDtLalpha, dExpDtLalphadL ) )

        floatVector dExpDtLalphadLp = dExpDtLalphadL * Dt * ( 1 - alpha );

        dExpDtLalphadL *= Dt * alpha;

        deformationGradient = floatVector( sot_dim, 0 );

        dFdFp = floatVector( sot_dim * sot_dim, 0 );

        dFdL = floatVector( sot_dim * sot_dim, 0 );

        dFdLp = floatVector( sot_dim * sot_dim, 0 );

        for ( unsigned int i = 0; i < dim; i++ ){

            for ( unsigned int j = 0; j < dim; j++ ){

                for ( unsigned int k = 0; k < dim; k++ ){

                    deformationGradient[ dim * i + k ] += expDtLalpha[ dim * i + j ] * previousDeformationGradient[ dim * j + k ];

                    dFdFp[ dim * sot_dim * i + sot_dim * j + dim * k + j ] += expDtLalpha[ dim * i + k ];

                    for ( unsigned int ab = 0; ab < sot_dim; ab++ ){

                        dFdL[ dim * sot_dim * i + sot_dim * k + ab ] += dExpDtLalphadL[ dim * sot_dim * i + sot_dim * j + ab ] * previousDeformationGradient[ dim * j + k ];
                        dFdLp[ dim * sot_dim * i + sot_dim * k + ab ] += dExpDtLalphadLp[ dim * sot_dim * i + sot_dim * j + ab ] * previousDeformationGradient[ dim * j + k ];


                    }

                }

            }

        }

    }

    void computeDCurrentNormalVectorDF( const floatVector &normalVector, const floatVector &F, floatVector &dNormalVectordF ){
        /*!
         * Compute the derivative of the normal vector in the current configuration w.r.t. the deformation gradient
         * 
         * \param &normalVector: The unit normal vector in the current configuration
         * \param &F: The deformation gradient
         * \param &dNormalVectordF: The derivative of the normal vector w.r.t. the deformation gradient
         */

        constexpr unsigned int dim = 3;
        constexpr unsigned int sot_dim = dim * dim;
        constexpr unsigned int tot_dim = dim * dim * dim;

        dNormalVectordF = floatVector( tot_dim, 0 );

        TARDIGRADE_ERROR_TOOLS_CHECK( F.size( ) == sot_dim, "The deformation gradient must be a second order tensor of size " + std::to_string( sot_dim ) + " and it has " + std::to_string( F.size( ) ) + " elements" );

        floatVector invF( sot_dim, 0 );

        Eigen::Map< const Eigen::Matrix< floatType, dim, dim, Eigen::RowMajor > > F_map( F.data( ), dim, dim );

        Eigen::Map< Eigen::Matrix< floatType, dim, dim, Eigen::RowMajor > > invF_map( invF.data( ), dim, dim );

        invF_map = F_map.inverse( );

        floatVector invF_n( dim, 0 );

        for ( unsigned int B = 0; B < dim; B++ ){

            for ( unsigned int j = 0; j < dim; j++ ){

                invF_n[ B ] += invF[ dim * B + j ] * normalVector[ j ];

            }

        }

        for ( unsigned int i = 0; i < dim; i++ ){

            for ( unsigned int b = 0; b < dim; b++ ){

                for ( unsigned int B = 0; B < dim; B++ ){

                    dNormalVectordF[ dim * dim * i + dim * b + B ] += normalVector[ i ] * normalVector[ b ] * invF_n[ B ]
                                                                    - normalVector[ b ] * invF[ dim * B + i ];

                }

            }

        }

    }

    void computeDCurrentAreaWeightedNormalVectorDF( const floatVector &normalVector, const floatVector &F, floatVector &dAreaWeightedNormalVectordF ){
        /*!
         * Compute the derivative of the area weighted normal vector w.r.t. the deformation gradient i.e.
         * 
         * \f$ \frac{\partial}{\partial F_{bB}} \left( n_i da \right) \f$
         * 
         * Note that if the user passes in the unit normal vector, then the result will be more convenient for the construction of
         * the jacobian of a surface integral in the current configuration.
         * 
         * \param &normalVector: The normal vector (a unit vector is likely what is desired)
         * \param &F: The deformation gradient
         * \param &dAreaWeightedNormalVectordF: The derivative of the area weighted normal vector w.r.t. the deformation gradient
         */

        constexpr unsigned int dim = 3;
        constexpr unsigned int sot_dim = dim * dim;
        constexpr unsigned int tot_dim = dim * dim * dim;

        dAreaWeightedNormalVectordF = floatVector( tot_dim, 0 );

        TARDIGRADE_ERROR_TOOLS_CHECK( F.size( ) == sot_dim, "The deformation gradient must be a second order tensor of size " + std::to_string( sot_dim ) + " and it has " + std::to_string( F.size( ) ) + " elements" );

        floatVector invF( sot_dim, 0 );

        Eigen::Map< const Eigen::Matrix< floatType, dim, dim, Eigen::RowMajor > > F_map( F.data( ), dim, dim );

        Eigen::Map< Eigen::Matrix< floatType, dim, dim, Eigen::RowMajor > > invF_map( invF.data( ), dim, dim );

        invF_map = F_map.inverse( );

        for ( unsigned int i = 0; i < dim; i++ ){

            for ( unsigned int b = 0; b < dim; b++ ){

                for ( unsigned int B = 0; B < dim; B++ ){

                    dAreaWeightedNormalVectordF[ dim * dim * i + dim * b + B ] += normalVector[ i ] * invF[ dim * B + b ]
                                                                                - normalVector[ b ] * invF[ dim * B + i ];

                }

            }

        }

    }

    void computeDCurrentAreaDF( const floatVector &normalVector, const floatVector &F, floatVector &dCurrentAreadF ){
        /*!
         * Compute the derivative of the current area w.r.t. the deformation gradient
         * 
         * \param &normalVector: The current unit normal vector
         * \param &F: The deformation gradient
         * \param &dCurrentAreadF: The derivative of the current surface area w.r.t. F
         */

        constexpr unsigned int dim = 3;
        constexpr unsigned int sot_dim = dim * dim;

        dCurrentAreadF = floatVector( sot_dim, 0 );

        TARDIGRADE_ERROR_TOOLS_CHECK( F.size( ) == sot_dim, "The deformation gradient must be a second order tensor of size " + std::to_string( sot_dim ) + " and it has " + std::to_string( F.size( ) ) + " elements" );

        floatVector invF( sot_dim, 0 );

        Eigen::Map< const Eigen::Matrix< floatType, dim, dim, Eigen::RowMajor > > F_map( F.data( ), dim, dim );

        Eigen::Map< Eigen::Matrix< floatType, dim, dim, Eigen::RowMajor > > invF_map( invF.data( ), dim, dim );

        invF_map = F_map.inverse( );

        floatVector invF_n( dim, 0 );

        for ( unsigned int B = 0; B < dim; B++ ){

            for ( unsigned int i = 0; i < dim; i++ ){

                invF_n[ B ] += invF[ dim * B + i ] * normalVector[ i ];

            }

        }

        for ( unsigned int B = 0; B < dim; B++ ){

            for ( unsigned int b = 0; b < dim; b++ ){

                dCurrentAreadF[ dim * b + B ] += invF[ dim * B + b ] - normalVector[ b ] * invF_n[ B ];

            }

        }

    }

    void computeDCurrentNormalVectorDGradU( const floatVector &normalVector, const floatVector &gradU, floatVector &dNormalVectordGradU, const bool isCurrent ){
        /*!
         * Compute the derivative of the normal vector in the current configuration w.r.t. the displacement gradient
         * 
         * \param &normalVector: The unit normal vector in the current configuration
         * \param &gradU: The displacement gradient
         * \param &dNormalVectordGradU: The derivative of the normal vector w.r.t. the displacement gradient
         * \param &isCurrent: Whether the displacement gradient is with respect to the reference or current configuration
         */

        constexpr unsigned int dim = 3;
        constexpr unsigned int sot_dim = dim * dim;
        constexpr unsigned int tot_dim = dim * dim * dim;

        floatVector F;

        floatVector dFdGradU;

        computeDeformationGradient( gradU, F, dFdGradU, isCurrent );

        floatVector dNormalVectordF;

        computeDCurrentNormalVectorDF( normalVector, F, dNormalVectordF );

        dNormalVectordGradU = floatVector( tot_dim, 0 );

        for ( unsigned int i = 0; i < dim; i++ ){

            for ( unsigned int j = 0; j < sot_dim; j++ ){

                for ( unsigned int k = 0; k < sot_dim; k++ ){

                    dNormalVectordGradU[ sot_dim * i + k ] += dNormalVectordF[ sot_dim * i + j ] * dFdGradU[ sot_dim * j + k ];

                }

            }

        }

    }

    void computeDCurrentAreaWeightedNormalVectorDGradU( const floatVector &normalVector, const floatVector &gradU, floatVector &dAreaWeightedNormalVectordGradU, const bool isCurrent ){
        /*!
         * Compute the derivative of the area weighted normal vector w.r.t. the displacement gradient
         * 
         * \f$ \frac{\partial}{\partial u_{i,j}} \left( n_i da \right) \f$
         * 
         * Note that if the user passes in the unit normal vector, then the result will be more convenient for the construction of
         * the jacobian of a surface integral in the current configuration.
         * 
         * \param &normalVector: The normal vector (a unit vector is likely what is desired)
         * \param &gradU: The displacement gradient
         * \param &dAreaWeightedNormalVectordGradU: The derivative of the area weighted normal vector w.r.t. the displacement gradient
         * \param &isCurrent: Whether the displacement gradient is with respect to the reference or current configuration
         */

        constexpr unsigned int dim = 3;
        constexpr unsigned int sot_dim = dim * dim;
        constexpr unsigned int tot_dim = dim * dim * dim;

        floatVector F;

        floatVector dFdGradU;

        computeDeformationGradient( gradU, F, dFdGradU, isCurrent );

        floatVector dAreaWeightedNormalVectordF;

        computeDCurrentAreaWeightedNormalVectorDF( normalVector, F, dAreaWeightedNormalVectordF );

        dAreaWeightedNormalVectordGradU = floatVector( tot_dim, 0 );

        for ( unsigned int i = 0; i < dim; i++ ){

            for ( unsigned int j = 0; j < sot_dim; j++ ){

                for ( unsigned int k = 0; k < sot_dim; k++ ){

                    dAreaWeightedNormalVectordGradU[ sot_dim * i + k ] += dAreaWeightedNormalVectordF[ sot_dim * i + j ] * dFdGradU[ sot_dim * j + k ];

                }

            }

        }

    }

    void computeDCurrentAreaDGradU( const floatVector &normalVector, const floatVector &gradU, floatVector &dCurrentAreadGradU, const bool isCurrent ){
        /*!
         * Compute the derivative of the current area w.r.t. the displacement gradient
         * 
         * \param &normalVector: The current unit normal vector
         * \param &gradU: The displacement gradient
         * \param &dCurrentAreadGradU: The derivative of the current surface area w.r.t. the displacement gradient
         * \param &isCurrent: Whether the displacement gradient is with respect to the reference or current configuration
         */

        constexpr unsigned int dim = 3;
        constexpr unsigned int sot_dim = dim * dim;

        floatVector F;

        floatVector dFdGradU;

        computeDeformationGradient( gradU, F, dFdGradU, isCurrent );

        floatVector dCurrentAreadF;

        computeDCurrentAreaDF( normalVector, F, dCurrentAreadF );

        dCurrentAreadGradU = floatVector( sot_dim, 0 );

        for ( unsigned int i = 0; i < sot_dim; i++ ){

            for ( unsigned int j = 0; j < sot_dim; j++ ){

                dCurrentAreadGradU[ j ] += dCurrentAreadF[ i ] * dFdGradU[ sot_dim * i + j ];

            }

        }

    }

}
