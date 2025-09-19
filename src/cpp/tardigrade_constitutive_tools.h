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

    floatType deltaDirac(const unsigned int i, const unsigned int j);

    template< unsigned int dim, class A_iterator, class Q_iterator, class rotatedA_iterator >
    void rotateMatrix(
        const A_iterator &A_begin, const A_iterator &A_end,
        const Q_iterator &Q_begin, const Q_iterator &Q_end,
        rotatedA_iterator rotatedA_begin, rotatedA_iterator rotatedA_end
    );

    template< class A_iterator, class Q_iterator, class rotatedA_iterator >
    void rotateMatrix(
        const A_iterator &A_begin, const A_iterator &A_end,
        const Q_iterator &Q_begin, const Q_iterator &Q_end,
        const unsigned int dim,
        rotatedA_iterator rotatedA_begin, rotatedA_iterator rotatedA_end
    );

    void rotateMatrix(const floatVector &A, const floatVector &Q, floatVector &rotatedA);

    template< unsigned int dim, class displacementGradient_iterator, class deformationGradient_iterator >
    void computeDeformationGradient(
        const displacementGradient_iterator &displacementGradient_begin, const displacementGradient_iterator &displacementGradient_end,
        deformationGradient_iterator        deformationGradient_begin,   deformationGradient_iterator        deformationGradient_end,
        const bool isCurrent
    );

    void computeDeformationGradient( const floatVector &displacementGradient, floatVector &F, const bool isCurrent );

    template< unsigned int dim, class displacementGradient_iterator, class deformationGradient_iterator, class dFdGradU_iterator >
    void computeDeformationGradient(
        const displacementGradient_iterator &displacementGradient_begin, const displacementGradient_iterator &displacementGradient_end,
        deformationGradient_iterator        deformationGradient_begin,   deformationGradient_iterator        deformationGradient_end,
        dFdGradU_iterator                   dFdGradU_begin,              dFdGradU_iterator                   dFdGradU_end,
        const bool isCurrent
    );

    void computeDeformationGradient( const floatVector &displacementGradient, floatVector &F, floatVector &dFdGradU, const bool isCurrent );

    template< unsigned int dim, class deformationGradient_iterator, class C_iterator >
    void computeRightCauchyGreen(
        const deformationGradient_iterator &deformationGradient_begin, const deformationGradient_iterator &deformationGradient_end,
        C_iterator C_begin, C_iterator C_end
    );

    template< unsigned int dim, class deformationGradient_iterator, class C_iterator, class dCdF_iterator >
    void computeRightCauchyGreen(
        const deformationGradient_iterator &deformationGradient_begin, const deformationGradient_iterator &deformationGradient_end,
        C_iterator C_begin,       C_iterator C_end,
        dCdF_iterator dCdF_begin, dCdF_iterator dCdF_end
    );

    void computeRightCauchyGreen( const floatVector &deformationGradient, floatVector &C );

    void computeRightCauchyGreen( const floatVector &deformationGradient, floatVector &C, floatVector &dCdF );

    void computeRightCauchyGreen( const floatVector &deformationGradient, floatVector &C, floatMatrix &dCdF );

    void computeGreenLagrangeStrain(const floatVector &deformationGradient, floatVector &E);

    void computeGreenLagrangeStrain(const floatVector &deformationGradient, floatVector &E, floatVector &dEdF);

    void computeGreenLagrangeStrain(const floatVector &deformationGradient, floatVector &E, floatMatrix &dEdF);

    void computeDGreenLagrangeStrainDF(const floatVector &deformationGradient, floatVector &dEdF);

    void computeDGreenLagrangeStrainDF(const floatVector &deformationGradient, floatMatrix &dEdF);

    void decomposeGreenLagrangeStrain(const floatVector &E, floatVector &Ebar, floatType &J);

    void decomposeGreenLagrangeStrain(const floatVector &E, floatVector &Ebar, floatType &J,
                                          floatVector &dEbardE, floatVector &dJdE);

    void decomposeGreenLagrangeStrain(const floatVector &E, floatVector &Ebar, floatType &J,
                                          floatMatrix &dEbardE, floatVector &dJdE);

    void mapPK2toCauchy(const floatVector &PK2Stress, const floatVector &deformationGradient, floatVector &cauchyStress);

    void WLF(const floatType &temperature, const floatVector &WLFParameters, floatType &factor);

    void WLF(const floatType &temperature, const floatVector &WLFParameters, floatType &factor, floatType &dfactordT);

    void computeDFDt(const floatVector &velocityGradient, const floatVector &deformationGradient, floatVector &DFDt);

    void computeDFDt(const floatVector &velocityGradient, const floatVector &deformationGradient, floatVector &DFDt,
                         floatVector &dDFDtdL, floatVector &dDFDtdF);

    void computeDFDt(const floatVector &velocityGradient, const floatVector &deformationGradient, floatVector &DFDt,
                         floatMatrix &dDFDtdL, floatMatrix &dDFDtdF);

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

    void evolveFExponentialMap( const floatType &Dt, const floatVector &previousDeformationGradient, const floatVector &Lp, const floatVector &L,
                                floatVector &deformationGradient, const floatType alpha=0.5 );

    void evolveFExponentialMap( const floatType &Dt, const floatVector &previousDeformationGradient, const floatVector &Lp, const floatVector &L,
                                floatVector &deformationGradient, floatVector &dFdL, const floatType alpha=0.5 );

    void evolveFExponentialMap( const floatType &Dt, const floatVector &previousDeformationGradient, const floatVector &Lp, const floatVector &L,
                                floatVector &deformationGradient, floatVector &dFdL, floatVector &dFdFp, floatVector &dFdLp, const floatType alpha=0.5 );

    floatType mac(const floatType &x);

    floatType mac(const floatType &x, floatType &dmacdx);

    void computeUnitNormal(const floatVector &A, floatVector &Anorm);

    void computeUnitNormal(const floatVector &A, floatVector &Anorm, floatVector &dAnormdA);

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
