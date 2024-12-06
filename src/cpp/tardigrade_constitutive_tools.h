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

#ifndef TARDIGRADE_CONSTITUTIVE_TOOLS_H
#define TARDIGRADE_CONSTITUTIVE_TOOLS_H

#define USE_EIGEN
#include<tardigrade_vector_tools.h>
#include<tardigrade_error_tools.h>

namespace tardigradeConstitutiveTools{

    typedef double floatType; //!< Define the float values type.
    typedef std::vector< floatType > floatVector; //!< Define a vector of floats
    typedef std::vector< std::vector< floatType > > floatMatrix; //!< Define a matrix of floats

    typedef tardigradeVectorTools::size_type size_type; //!< Define the type of the size variable

    floatType deltaDirac(const unsigned int i, const unsigned int j);

    template<int dim, class v_in, class v_out>
    void rotateMatrix( const v_in &A_begin, const v_in &A_end, const v_in &Q_begin, const v_in &Q_end,
                       v_out temp_begin, v_out temp_end,
                       v_out rotatedA_begin, v_out rotatedA_end );

    template<class v_in, class v_out>
    void rotateMatrix( const v_in &A_begin, const v_in &A_end, const v_in &Q_begin, const v_in &Q_end,
                       const unsigned int dim,
                       v_out temp_Begin,     v_out temp_end,
                       v_out rotatedA_begin, v_out rotatedA_end );

    void rotateMatrix(const floatVector &A, const floatVector &Q, floatVector &rotatedA);

    template<int dim, class v_in, class v_out>
    void computeDeformationGradient( const v_in &displacementGradient_begin, const v_in &displacementGradient_end,
                                     v_out F_begin, v_out F_end, const bool isCurrent );

    template<class v_in, class v_out>
    void computeDeformationGradient( const v_in &displacementGradient_begin, const v_in &displacementGradient_end,
                                     v_out F_begin, v_out F_end, const bool isCurrent );

    void computeDeformationGradient( const floatVector &displacementGradient, floatVector &F, const bool isCurrent );

    template<int dim, class v_in, class v_out, class M_out>
    void computeDeformationGradient( const v_in &displacementGradient_Begin, const v_in &displacementGradient_end,
                                     v_out F_begin, v_out F_end, M_out dFdGradU_begin, M_out dFdGradU_end, const bool isCurrent );

    template<class v_in, class v_out, class M_out>
    void computeDeformationGradient( const v_in &displacementGradient_Begin, const v_in &displacementGradient_end,
                                     v_out F_begin, v_out F_end, M_out dFdGradU_begin, M_out dFdGradU_end, const bool isCurrent );

    void computeDeformationGradient( const floatVector &displacementGradient, floatVector &F, floatVector &dFdGradU, const bool isCurrent );

    template<int dim, class v_in, class v_out>
    void computeRightCauchyGreen( const v_in &deformationGradient_begin, const v_in &deformationGradient_end,
                                  v_out C_begin, v_out C_end );

    template<int dim, class v_in, class v_out, class M_out>
    void computeRightCauchyGreen( const v_in &deformationGradient_begin, const v_in &deformationGradient_end,
                                  v_out C_begin, v_out C_end,
                                  M_out dCdF_begin, M_out dCdF_end );

    void computeRightCauchyGreen( const floatVector &deformationGradient, floatVector &C );

    void computeRightCauchyGreen( const floatVector &deformationGradient, floatVector &C, floatVector &dCdF );

    void computeRightCauchyGreen( const floatVector &deformationGradient, floatVector &C, floatMatrix &dCdF );

    template<int dim, class v_in, class v_out>
    void computeGreenLagrangeStrain( const v_in &deformationGradient_begin, const v_in &deformationGradient_end,
                                     v_out E_begin, v_out E_end );

    template<int dim, class v_in, class v_out, class M_out>
    void computeGreenLagrangeStrain( const v_in &deformationGradient_begin, const v_in &deformationGradient_end,
                                     v_out E_begin, v_out E_end,
                                     M_out dEdF_begin, M_out dEdF_end );

    void computeGreenLagrangeStrain(const floatVector &deformationGradient, floatVector &E);

    void computeGreenLagrangeStrain(const floatVector &deformationGradient, floatVector &E, floatVector &dEdF);

    void computeGreenLagrangeStrain(const floatVector &deformationGradient, floatVector &E, floatMatrix &dEdF);

    void computeDGreenLagrangeStrainDF(const floatVector &deformationGradient, floatVector &dEdF);

    template<int dim, class v_in, class M_out>
    void computeDGreenLagrangeStrainDF( const v_in &deformationGradient_begin, const v_in &deformationGradient_end,
                                        M_out dEdF_begin, M_out dEdF_end );

    void computeDGreenLagrangeStrainDF(const floatVector &deformationGradient, floatMatrix &dEdF);

    template<int dim, class v_in, class v_out, typename T>
    void decomposeGreenLagrangeStrain( const v_in &E_begin, const v_in &E_end, v_out Ebar_begin, v_out Ebar_end, T &J );

    template<int dim, class v_in, class v_out, class M_out, typename T>
    void decomposeGreenLagrangeStrain( const v_in &E_begin, const v_in &E_end, v_out Ebar_begin, v_out Ebar_end, T &J,
                                       M_out dEbardE_begin, M_out dEbardE_end, v_out dJdE_begin, v_out dJdE_end );

    void decomposeGreenLagrangeStrain(const floatVector &E, floatVector &Ebar, floatType &J);

    void decomposeGreenLagrangeStrain(const floatVector &E, floatVector &Ebar, floatType &J,
                                          floatVector &dEbardE, floatVector &dJdE);

    void decomposeGreenLagrangeStrain(const floatVector &E, floatVector &Ebar, floatType &J,
                                          floatMatrix &dEbardE, floatVector &dJdE);

    template<int dim, class v_in, class v_out>
    void mapPK2toCauchy( const v_in &PK2Stress_begin, const v_in &PK2Stress_end,
                         const v_in &deformationGradient_begin, const v_in &deformationGradient_end,
                         v_out temp_begin,         v_out temp_end,
                         v_out cauchyStress_begin, v_out cauchyStress_end );

    void mapPK2toCauchy(const floatVector &PK2Stress, const floatVector &deformationGradient, floatVector &cauchyStress);

    template<typename T, class v_in>
    void WLF( const T &temperature, const v_in &parameters_begin, const v_in &parameters_end, T &factor );

    template<typename T, class v_in>
    void WLF( const T &temperature, const v_in &parameters_begin, const v_in &parameters_end, T &factor, T &dfactordT );

    void WLF(const floatType &temperature, const floatVector &WLFParameters, floatType &factor);

    void WLF(const floatType &temperature, const floatVector &WLFParameters, floatType &factor, floatType &dfactordT);

    template<int dim, class v_in, class v_out>
    void computeDFDt(const v_in &velocityGradient_begin,    const v_in &velocityGradient_end,
                     const v_in &deformationGradient_begin, const v_in &deformationGradient_end,
                     v_out DFDt_begin, v_out DFDt_end );

    template<int dim, class v_in, class v_out, class M_out>
    void computeDFDt(const v_in &velocityGradient_begin,    const v_in &velocityGradient_end,
                     const v_in &deformationGradient_begin, const v_in &deformationGradient_end,
                     v_out DFDt_begin, v_out DFDt_end, M_out dDFDtdL_begin, M_out dDFDtdL_end,
                     M_out dDFDtdF_begin, M_out dDFDtdF_end );

    void computeDFDt(const floatVector &velocityGradient, const floatVector &deformationGradient, floatVector &DFDt);

    void computeDFDt(const floatVector &velocityGradient, const floatVector &deformationGradient, floatVector &DFDt,
                         floatVector &dDFDtdL, floatVector &dDFDtdF);

    void computeDFDt(const floatVector &velocityGradient, const floatVector &deformationGradient, floatVector &DFDt,
                         floatMatrix &dDFDtdL, floatMatrix &dDFDtdF);

    template<typename T, class v_in, class v_out>
    void midpointEvolution( const T &Dt, const v_in &Ap_begin, const v_in &Ap_end, const v_in &DApDt_begin, const v_in &DApDt_end,
                            const v_in &DADt_begin, const v_in &DADt_end, v_out dA_begin, v_out dA_end, v_out A_begin, v_out A_end,
                            const v_in &alpha_begin, const v_in &alpha_end );

    template<typename T, class v_in, class v_out, class M_out>
    void midpointEvolution( const T &Dt, const v_in &Ap_begin, const v_in &Ap_end, const v_in &DApDt_begin, const v_in &DApDt_end,
                            const v_in &DADt_begin, const v_in &DADt_end, v_out dA_begin, v_out dA_end, v_out A_begin, v_out A_end,
                            M_out DADADt_begin, M_out DADADt_end,
                            const v_in &alpha_begin, const v_in &alpha_end );

    template<typename T, class v_in, class v_out, class M_out>
    void midpointEvolution( const T &Dt, const v_in &Ap_begin, const v_in &Ap_end, const v_in &DApDt_begin, const v_in &DApDt_end,
                            const v_in &DADt_begin, const v_in &DADt_end, v_out dA_begin, v_out dA_end, v_out A_begin, v_out A_end,
                            M_out DADADt_begin, M_out DADADt_end, M_out DADADtp_begin, M_out DADADtp_end,
                            const v_in &alpha_begin, const v_in &alpha_end );

    template<typename T, class v_in, class v_out>
    void midpointEvolution( const T &Dt, const v_in &Ap_begin, const v_in &Ap_end, const v_in &DApDt_begin, const v_in &DApDt_end,
                            const v_in &DADt_begin, const v_in &DADt_end, v_out dA_begin, v_out dA_end, v_out A_begin, v_out A_end,
                            const floatType alpha=0.5 );

    template<typename T, class v_in, class v_out, class M_out>
    void midpointEvolution( const T &Dt, const v_in &Ap_begin, const v_in &Ap_end, const v_in &DApDt_begin, const v_in &DApDt_end,
                            const v_in &DADt_begin, const v_in &DADt_end, v_out dA_begin, v_out dA_end, v_out A_begin, v_out A_end,
                            M_out DADADt_begin, M_out DADADt_end,
                            const floatType alpha=0.5 );

    template<typename T, class v_in, class v_out, class M_out>
    void midpointEvolution( const T &Dt, const v_in &Ap_begin, const v_in &Ap_end, const v_in &DApDt_begin, const v_in &DApDt_end,
                            const v_in &DADt_begin, const v_in &DADt_end, v_out dA_begin, v_out dA_end, v_out A_begin, v_out A_end,
                            M_out DADADt_begin, M_out DADADt_end, M_out DADADtp_begin, M_out DADADtp_end,
                            const floatType alpha=0.5 );

    void midpointEvolution(const floatType &Dt, const floatVector &Ap, const floatVector &DApDt, const floatVector &DADt,
                               floatVector &dA, floatVector &A, const floatVector &alpha);

    void midpointEvolutionFlatJ(const floatType &Dt, const floatVector &Ap, const floatVector &DApDt, const floatVector &DADt,
                                    floatVector &dA, floatVector &A, floatVector &DADADt, const floatVector &alpha);

    void midpointEvolutionFlatJ(const floatType &Dt, const floatVector &Ap, const floatVector &DApDt, const floatVector &DADt,
                                    floatVector &dA, floatVector &A, floatVector &DADADt, floatVector &DADADtp,
                                    const floatVector &alpha);

    void midpointEvolution(const floatType &Dt, const floatVector &Ap, const floatVector &DApDt, const floatVector &DADt,
                               floatVector &dA, floatVector &A, floatMatrix &DADADt, const floatVector &alpha);

    void midpointEvolution(const floatType &Dt, const floatVector &Ap, const floatVector &DApDt, const floatVector &DADt,
                               floatVector &dA, floatVector &A, floatMatrix &DADADt, floatMatrix &DADADtp,
                               const floatVector &alpha);

    void midpointEvolution(const floatType &Dt, const floatVector &Ap, const floatVector &DApDt, const floatVector &DADt,
                               floatVector &dA, floatVector &A, const floatType alpha=0.5);

    void midpointEvolutionFlatJ(const floatType &Dt, const floatVector &Ap, const floatVector &DApDt, const floatVector &DADt,
                                    floatVector &dA, floatVector &A, floatVector &DADADt, const floatType alpha=0.5);

    void midpointEvolutionFlatJ(const floatType &Dt, const floatVector &Ap, const floatVector &DApDt, const floatVector &DADt,
                                    floatVector &dA, floatVector &A, floatVector &DADADt, floatVector &DADADtp,
                                    const floatType alpha=0.5);

    void midpointEvolution(const floatType &Dt, const floatVector &Ap, const floatVector &DApDt, const floatVector &DADt,
                               floatVector &dA, floatVector &A, floatMatrix &DADADt, const floatType alpha=0.5);

    void midpointEvolution(const floatType &Dt, const floatVector &Ap, const floatVector &DApDt, const floatVector &DADt,
                               floatVector &dA, floatVector &A, floatMatrix &DADADt, floatMatrix &DADADtp,
                               const floatType alpha=0.5);

    template<int dim, typename T, class v_in, class v_out>
    void evolveF(const T &Dt, const v_in &previousDeformationGradient_begin, const v_in &previousDeformationGradient_end,
                 const v_in &Lp_begin, const v_in &Lp_end, const v_in &L_begin, const v_in &L_end,
                 v_out dF_begin, v_out dF_end, v_out deformationGradient_begin, v_out deformationGradient_end,
                 const floatType alpha=0.5, const unsigned int mode = 1 );

    template<int dim, typename T, class v_in, class v_out, class M_out>
    void evolveF(const T &Dt, const v_in &previousDeformationGradient_begin, const v_in &previousDeformationGradient_end,
                 const v_in &Lp_begin, const v_in &Lp_end, const v_in &L_begin, const v_in &L_end,
                 v_out dF_begin, v_out dF_end, v_out deformationGradient_begin, v_out deformationGradient_end,
                 M_out dFdL_begin, M_out dFdL_end,
                 const floatType alpha=0.5, const unsigned int mode = 1 );

    template<int dim, typename T, class v_in, class v_out, class M_out>
    void evolveF(const T &Dt, const v_in &previousDeformationGradient_begin, const v_in &previousDeformationGradient_end,
                 const v_in &Lp_begin, const v_in &Lp_end, const v_in &L_begin, const v_in &L_end,
                 v_out dF_begin, v_out dF_end, v_out deformationGradient_begin, v_out deformationGradient_end,
                 M_out dFdL_begin, M_out dFdL_end, M_out ddFdFp_begin, M_out ddFdFp_end, M_out dFdFp_begin, M_out dFdFp_end, M_out dFdLp_begin, M_out dFdLp_end,
                 const floatType alpha=0.5, const unsigned int mode = 1 );

    void evolveF(const floatType &Dt, const floatVector &previousDeformationGradient, const floatVector &Lp, const floatVector &L,
                     floatVector &dF, floatVector &deformationGradient, const floatType alpha=0.5, const unsigned int mode = 1);

    void evolveFFlatJ(const floatType &Dt, const floatVector &previousDeformationGradient, const floatVector &Lp, const floatVector &L,
                          floatVector &dF, floatVector &deformationGradient, floatVector &dFdL, const floatType alpha=0.5, const unsigned int mode = 1);

    void evolveF(const floatType &Dt, const floatVector &previousDeformationGradient, const floatVector &Lp, const floatVector &L,
                     floatVector &dF, floatVector &deformationGradient, floatMatrix &dFdL, const floatType alpha=0.5, const unsigned int mode = 1);

    void evolveF(const floatType &Dt, const floatVector &previousDeformationGradient, const floatVector &Lp, const floatVector &L,
                     floatVector &deformationGradient, floatType alpha=0.5, const unsigned int mode = 1);

    void evolveFFlatJ(const floatType &Dt, const floatVector &previousDeformationGradient, const floatVector &Lp, const floatVector &L,
                          floatVector &deformationGradient, floatVector &dFdL, const floatType alpha=0.5, const unsigned int mode = 1);

    void evolveF(const floatType &Dt, const floatVector &previousDeformationGradient, const floatVector &Lp, const floatVector &L,
                     floatVector &deformationGradient, floatMatrix &dFdL, const floatType alpha=0.5, const unsigned int mode = 1);

    void evolveFFlatJ(const floatType &Dt, const floatVector &previousDeformationGradient, const floatVector &Lp, const floatVector &L,
                          floatVector &deformationGradient, floatVector &dFdL, floatVector &dFdFp, floatVector &dFdLp, const floatType alpha=0.5, const unsigned int mode = 1);

    void evolveF(const floatType &Dt, const floatVector &previousDeformationGradient, const floatVector &Lp, const floatVector &L,
                     floatVector &deformationGradient, floatMatrix &dFdL, floatMatrix &dFdFp, floatMatrix &dFdLp, const floatType alpha=0.5, const unsigned int mode = 1);

    void evolveFFlatJ(const floatType &Dt, const floatVector &previousDeformationGradient, const floatVector &Lp, const floatVector &L,
                          floatVector &dF, floatVector &deformationGradient, floatVector &dFdL, floatVector &ddFdFp, floatVector &dFdFp, floatVector &dFdLp, const floatType alpha=0.5, const unsigned int mode = 1);

    void evolveF(const floatType &Dt, const floatVector &previousDeformationGradient, const floatVector &Lp, const floatVector &L,
                     floatVector &dF, floatVector &deformationGradient, floatMatrix &dFdL, floatMatrix &ddFdFp, floatMatrix &dFdFp, floatMatrix &dFdLp, const floatType alpha=0.5, const unsigned int mode = 1);

    template<int dim, typename T, class v_in, class v_out>
    void evolveFExponentialMap( const T &Dt, const v_in &previousDeformationGradient_begin, const v_in &previousDeformationGradient_end,
                                const v_in &Lp_begin, const v_in &Lp_end, const v_in &L_begin, const v_in &L_end,
                                v_out tempVec1_begin, v_out tempVec1_end,
                                v_out tempVec2_begin, v_out tempVec2_end,
                                v_out tempVec3_begin, v_out tempVec3_end,
                                v_out tempVec4_begin, v_out tempVec4_end,
                                v_out deformationGradient_begin, v_out deformationGradient_end, const floatType alpha = 0.5 );

    template<int dim, typename T, class v_in, class v_out, class M_out>
    void evolveFExponentialMap( const T &Dt, const v_in &previousDeformationGradient_begin, const v_in &previousDeformationGradient_end,
                                const v_in &Lp_begin, const v_in &Lp_end, const v_in &L_begin, const v_in &L_end,
                                v_out tempVec1_begin, v_out tempVec1_end,
                                v_out tempVec2_begin, v_out tempVec2_end,
                                v_out tempVec3_begin, v_out tempVec3_end,
                                v_out tempVec4_begin, v_out tempVec4_end,
                                M_out tempMat1_begin, v_out tempMat1_end,
                                M_out tempMat2_begin, v_out tempMat2_end,
                                v_out deformationGradient_begin, v_out deformationGradient_end,
                                M_out dFdL_begin, M_out dFdL_end, const floatType alpha = 0.5 );

//    template<int dim, typename T, class v_in, class v_out, class M_out>
//    void evolveFExponentialMap( const T &Dt, const v_in &previousDeformationGradient_begin, const v_in &previousDeformationGradient_end,
//                                const v_in &Lp_begin, const v_in &Lp_end, const v_in &L_begin, const v_in &L_end,
//                                v_out deformationGradient_begin, v_out deformationGradient_end,
//                                M_out dFdL_begin, M_out dFdL_end, const floatType alpha = 0.5 );

    template<int dim, typename T, class v_in, class v_out, class M_out>
    void evolveFExponentialMap( const T &Dt, const v_in &previousDeformationGradient_begin, const v_in &previousDeformationGradient_end,
                                const v_in &Lp_begin, const v_in &Lp_end, const v_in &L_begin, const v_in &L_end,
                                v_out deformationGradient_begin, v_out deformationGradient_end,
                                M_out dFdL_begin, M_out dFdL_end, M_out dFdLp_begin, M_out dFdLp_end, const floatType alpha = 0.5 );

    void evolveFExponentialMap( const floatType &Dt, const floatVector &previousDeformationGradient, const floatVector &Lp, const floatVector &L,
                                floatVector &deformationGradient, const floatType alpha=0.5 );

    void evolveFExponentialMap( const floatType &Dt, const floatVector &previousDeformationGradient, const floatVector &Lp, const floatVector &L,
                                floatVector &deformationGradient, floatVector &dFdL, const floatType alpha=0.5 );

    void evolveFExponentialMap( const floatType &Dt, const floatVector &previousDeformationGradient, const floatVector &Lp, const floatVector &L,
                                floatVector &deformationGradient, floatVector &dFdL, floatVector &dFdFp, floatVector &dFdLp, const floatType alpha=0.5 );

    floatType mac(const floatType &x);

    floatType mac(const floatType &x, floatType &dmacdx);

    void computeUnitNormal(const floatVector &A, floatVector &Anorm);

    void computeUnitNormal(const floatVector &A, floatVector &Anorm, floatMatrix &dAnormdA);

    void pullBackVelocityGradient(const floatVector &velocityGradient, const floatVector &deformationGradient,
                                      floatVector &pulledBackVelocityGradient);

    void pullBackVelocityGradient(const floatVector &velocityGradient, const floatVector &deformationGradient,
                                      floatVector &pulledBackVelocityGradient, floatVector &dPullBackLdL,
                                      floatVector &dPullBackLdF);

    void pullBackVelocityGradient(const floatVector &velocityGradient, const floatVector &deformationGradient,
                                      floatVector &pulledBackVelocityGradient, floatMatrix &dPullBackLdL,
                                      floatMatrix &dPullBackLdF);

    void quadraticThermalExpansion(const floatType &temperature, const floatType &referenceTemperature,
                                       const floatVector &linearParameters, const floatVector &quadraticParameters,
                                       floatVector &thermalExpansion);

    void quadraticThermalExpansion(const floatType &temperature, const floatType &referenceTemperature,
                                       const floatVector &linearParameters, const floatVector &quadraticParameters,
                                       floatVector &thermalExpansion, floatVector &thermalExpansionJacobian);

    void pushForwardGreenLagrangeStrain(const floatVector &greenLagrangeStrain, const floatVector &deformationGradient,
                                            floatVector &almansiStrain);

    void pushForwardGreenLagrangeStrain(const floatVector &greenLagrangeStrain, const floatVector &deformationGradient,
                                            floatVector &almansiStrain, floatVector &dAlmansiStraindE, floatVector &dAlmansiStraindF);

    void pushForwardGreenLagrangeStrain(const floatVector &greenLagrangeStrain, const floatVector &deformationGradient,
                                            floatVector &almansiStrain, floatMatrix &dAlmansiStraindE, floatMatrix &dAlmansiStraindF);

    void pullBackAlmansiStrain( const floatVector &almansiStrain, const floatVector &deformationGradient,
                                    floatVector &greenLagrangeStrain );

    void pullBackAlmansiStrain( const floatVector &almansiStrain, const floatVector &deformationGradient,
                                    floatVector &greenLagrangeStrain, floatVector &dEde, floatVector &dEdF );

    void pullBackAlmansiStrain( const floatVector &almansiStrain, const floatVector &deformationGradient,
                                    floatVector &greenLagrangeStrain, floatMatrix &dEde, floatMatrix &dEdF );

    void computeSymmetricPart( const floatVector &A, floatVector &symmA, unsigned int &dim );

    void computeSymmetricPart( const floatVector &A, floatVector &symmA );

    void computeSymmetricPart( const floatVector &A, floatVector &symmA, floatVector &dSymmAdA );

    void computeSymmetricPart( const floatVector &A, floatVector &symmA, floatMatrix &dSymmAdA );

    void pushForwardPK2Stress( const floatVector &PK2, const floatVector &F, floatVector &cauchyStress );

    void pushForwardPK2Stress( const floatVector &PK2, const floatVector &F, floatVector &cauchyStress,
                                   floatVector &dCauchyStressdPK2, floatVector &dCauchyStressdF );

    void pushForwardPK2Stress( const floatVector &PK2, const floatVector &F, floatVector &cauchyStress,
                                   floatMatrix &dCauchyStressdPK2, floatMatrix &dCauchyStressdF );

    void pullBackCauchyStress( const floatVector &cauchyStress, const floatVector &F, floatVector &PK2 );

    void pullBackCauchyStress( const floatVector &cauchyStress, const floatVector &F, floatVector &PK2,
                                   floatVector &dPK2dCauchyStress, floatVector &dPK2dF );

    void pullBackCauchyStress( const floatVector &cauchyStress, const floatVector &F, floatVector &PK2,
                                   floatMatrix &dPK2dCauchyStress, floatMatrix &dPK2dF );

    void computeDCurrentNormalVectorDF( const floatVector &normalVector, const floatVector &F, floatVector &dNormalVectordF );

    void computeDCurrentAreaWeightedNormalVectorDF( const floatVector &normalVector, const floatVector &F, floatVector &dAreaWeightedNormalVectordF );

    void computeDCurrentAreaDF( const floatVector &normalVector, const floatVector &F, floatVector &dCurrentAreadF );

    void computeDCurrentNormalVectorDGradU( const floatVector &normalVector, const floatVector &gradU, floatVector &dNormalVectordGradU, const bool isCurrent = true );

    void computeDCurrentAreaWeightedNormalVectorDGradU( const floatVector &normalVector, const floatVector &gradU, floatVector &dAreaWeightedNormalVectordGradU, const bool isCurrent = true );

    void computeDCurrentAreaDGradU( const floatVector &normalVector, const floatVector &gradU, floatVector &dCurrentAreadGradU, const bool isCurrent = true );

}

#ifdef TARDIGRADE_HEADER_ONLY
    #include "tardigrade_constitutive_tools.cpp"
#endif

#endif
