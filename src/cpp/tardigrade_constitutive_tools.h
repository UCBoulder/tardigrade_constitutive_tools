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
#include "tardigrade_error_tools.h"
#include "tardigrade_vector_tools.h"

namespace tardigradeConstitutiveTools {

    typedef double                               floatType;    //!< Define the float values type.
    typedef std::vector<floatType>               floatVector;  //!< Define a vector of floats
    typedef std::vector<std::vector<floatType> > floatMatrix;  //!< Define a matrix of floats

    TARDIGRADE_OPTIONAL_INLINE floatType deltaDirac(const unsigned int i, const unsigned int j);

    template <unsigned int dim, class A_iterator, class Q_iterator, class rotatedA_iterator>
    TARDIGRADE_OPTIONAL_INLINE void rotateMatrix(const A_iterator &A_begin, const A_iterator &A_end,
                                                 const Q_iterator &Q_begin, const Q_iterator &Q_end,
                                                 rotatedA_iterator rotatedA_begin, rotatedA_iterator rotatedA_end);

    template <class A_iterator, class Q_iterator, class rotatedA_iterator>
    TARDIGRADE_OPTIONAL_INLINE void rotateMatrix(const A_iterator &A_begin, const A_iterator &A_end,
                                                 const Q_iterator &Q_begin, const Q_iterator &Q_end,
                                                 const unsigned int dim, rotatedA_iterator rotatedA_begin,
                                                 rotatedA_iterator rotatedA_end);

    TARDIGRADE_OPTIONAL_INLINE void rotateMatrix(const floatVector &A, const floatVector &Q, floatVector &rotatedA);

    template <unsigned int dim, class displacementGradient_iterator, class deformationGradient_iterator>
    TARDIGRADE_OPTIONAL_INLINE void computeDeformationGradient(
        const displacementGradient_iterator &displacementGradient_begin,
        const displacementGradient_iterator &displacementGradient_end,
        deformationGradient_iterator deformationGradient_begin, deformationGradient_iterator deformationGradient_end,
        const bool isCurrent);

    TARDIGRADE_OPTIONAL_INLINE void computeDeformationGradient(const floatVector &displacementGradient, floatVector &F,
                                                               const bool isCurrent);

    template <unsigned int dim, class displacementGradient_iterator, class deformationGradient_iterator,
              class dFdGradU_iterator>
    TARDIGRADE_OPTIONAL_INLINE void computeDeformationGradient(
        const displacementGradient_iterator &displacementGradient_begin,
        const displacementGradient_iterator &displacementGradient_end,
        deformationGradient_iterator deformationGradient_begin, deformationGradient_iterator deformationGradient_end,
        dFdGradU_iterator dFdGradU_begin, dFdGradU_iterator dFdGradU_end, const bool isCurrent);

    TARDIGRADE_OPTIONAL_INLINE void computeDeformationGradient(const floatVector &displacementGradient, floatVector &F,
                                                               floatVector &dFdGradU, const bool isCurrent);

    template <unsigned int dim, class deformationGradient_iterator, class C_iterator>
    TARDIGRADE_OPTIONAL_INLINE void computeRightCauchyGreen(
        const deformationGradient_iterator &deformationGradient_begin,
        const deformationGradient_iterator &deformationGradient_end, C_iterator C_begin, C_iterator C_end);

    template <unsigned int dim, class deformationGradient_iterator, class C_iterator, class dCdF_iterator>
    TARDIGRADE_OPTIONAL_INLINE void computeRightCauchyGreen(
        const deformationGradient_iterator &deformationGradient_begin,
        const deformationGradient_iterator &deformationGradient_end, C_iterator C_begin, C_iterator C_end,
        dCdF_iterator dCdF_begin, dCdF_iterator dCdF_end);

    TARDIGRADE_OPTIONAL_INLINE void computeRightCauchyGreen(const floatVector &deformationGradient, floatVector &C);

    TARDIGRADE_OPTIONAL_INLINE void computeRightCauchyGreen(const floatVector &deformationGradient, floatVector &C,
                                                            floatVector &dCdF);

    TARDIGRADE_OPTIONAL_INLINE void computeRightCauchyGreen(const floatVector &deformationGradient, floatVector &C,
                                                            floatMatrix &dCdF);

    template <unsigned int dim, class deformationGradient_iterator, class E_iterator>
    TARDIGRADE_OPTIONAL_INLINE void computeGreenLagrangeStrain(
        const deformationGradient_iterator &deformationGradient_begin,
        const deformationGradient_iterator &deformationGradient_end, E_iterator E_begin, E_iterator E_end);

    template <unsigned int dim, class deformationGradient_iterator, class E_iterator, class dEdF_iterator>
    TARDIGRADE_OPTIONAL_INLINE void computeGreenLagrangeStrain(
        const deformationGradient_iterator &deformationGradient_begin,
        const deformationGradient_iterator &deformationGradient_end, E_iterator E_begin, E_iterator E_end,
        dEdF_iterator dEdF_begin, dEdF_iterator dEdF_end);

    TARDIGRADE_OPTIONAL_INLINE void computeGreenLagrangeStrain(const floatVector &deformationGradient, floatVector &E);

    TARDIGRADE_OPTIONAL_INLINE void computeGreenLagrangeStrain(const floatVector &deformationGradient, floatVector &E,
                                                               floatVector &dEdF);

    TARDIGRADE_OPTIONAL_INLINE void computeGreenLagrangeStrain(const floatVector &deformationGradient, floatVector &E,
                                                               floatMatrix &dEdF);

    template <unsigned int dim, class deformationGradient_iterator, class dEdF_iterator>
    TARDIGRADE_OPTIONAL_INLINE void computeDGreenLagrangeStrainDF(
        const deformationGradient_iterator &deformationGradient_begin,
        const deformationGradient_iterator &deformationGradient_end, dEdF_iterator dEdF_begin, dEdF_iterator dEdF_end);

    TARDIGRADE_OPTIONAL_INLINE void computeDGreenLagrangeStrainDF(const floatVector &deformationGradient,
                                                                  floatVector       &dEdF);

    TARDIGRADE_OPTIONAL_INLINE void computeDGreenLagrangeStrainDF(const floatVector &deformationGradient,
                                                                  floatMatrix       &dEdF);

    template <unsigned int dim, class E_iterator, class Ebar_iterator, typename J_type>
    TARDIGRADE_OPTIONAL_INLINE void decomposeGreenLagrangeStrain(const E_iterator &E_begin, const E_iterator &E_end,
                                                                 Ebar_iterator Ebar_begin, Ebar_iterator Ebar_end,
                                                                 J_type &J);

    TARDIGRADE_OPTIONAL_INLINE void decomposeGreenLagrangeStrain(const floatVector &E, floatVector &Ebar, floatType &J);

    template <unsigned int dim, class E_iterator, class Ebar_iterator, typename J_type, class dEbardE_iterator,
              class dJdE_iterator>
    TARDIGRADE_OPTIONAL_INLINE void decomposeGreenLagrangeStrain(const E_iterator &E_begin, const E_iterator &E_end,
                                                                 Ebar_iterator Ebar_begin, Ebar_iterator Ebar_end,
                                                                 J_type &J, dEbardE_iterator dEbardE_begin,
                                                                 dEbardE_iterator dEbardE_end, dJdE_iterator dJdE_begin,
                                                                 dJdE_iterator dJdE_end);

    TARDIGRADE_OPTIONAL_INLINE void decomposeGreenLagrangeStrain(const floatVector &E, floatVector &Ebar, floatType &J,
                                                                 floatVector &dEbardE, floatVector &dJdE);

    TARDIGRADE_OPTIONAL_INLINE void decomposeGreenLagrangeStrain(const floatVector &E, floatVector &Ebar, floatType &J,
                                                                 floatMatrix &dEbardE, floatVector &dJdE);

    template <unsigned int dim, class PK2Stress_iterator, class deformationGradient_iterator,
              class cauchyStress_iterator>
    TARDIGRADE_OPTIONAL_INLINE void mapPK2toCauchy(const PK2Stress_iterator           &PK2Stress_begin,
                                                   const PK2Stress_iterator           &PK2Stress_end,
                                                   const deformationGradient_iterator &deformationGradient_begin,
                                                   const deformationGradient_iterator &deformationGradient_end,
                                                   cauchyStress_iterator               cauchyStress_begin,
                                                   cauchyStress_iterator               cauchyStress_end);

    TARDIGRADE_OPTIONAL_INLINE void mapPK2toCauchy(const floatVector &PK2Stress, const floatVector &deformationGradient,
                                                   floatVector &cauchyStress);

    template <typename temperature_type, class WLFParameters_iterator, typename factor_type>
    TARDIGRADE_OPTIONAL_INLINE void WLF(const temperature_type       &temperature,
                                        const WLFParameters_iterator &WLFParameters_begin,
                                        const WLFParameters_iterator &WLFParameters_end, factor_type &factor);

    template <typename temperature_type, class WLFParameters_iterator, typename factor_type, typename dFactordT_type>
    TARDIGRADE_OPTIONAL_INLINE void WLF(const temperature_type       &temperature,
                                        const WLFParameters_iterator &WLFParameters_begin,
                                        const WLFParameters_iterator &WLFParameters_end, factor_type &factor,
                                        dFactordT_type &dFactordT);

    TARDIGRADE_OPTIONAL_INLINE void WLF(const floatType &temperature, const floatVector &WLFParameters,
                                        floatType &factor);

    TARDIGRADE_OPTIONAL_INLINE void WLF(const floatType &temperature, const floatVector &WLFParameters,
                                        floatType &factor, floatType &dFactordT);

    template <unsigned int dim, class velocityGradient_iterator, class deformationGradient_iterator,
              class DFDt_iterator>
    TARDIGRADE_OPTIONAL_INLINE void computeDFDt(const velocityGradient_iterator    &velocityGradient_begin,
                                                const velocityGradient_iterator    &velocityGradient_end,
                                                const deformationGradient_iterator &deformationGradient_begin,
                                                const deformationGradient_iterator &deformationGradient_end,
                                                DFDt_iterator DFDt_begin, DFDt_iterator DFDt_end);

    template <unsigned int dim, class velocityGradient_iterator, class deformationGradient_iterator,
              class DFDt_iterator, class dDFDtdL_iterator, class dDFDtdF_iterator>
    TARDIGRADE_OPTIONAL_INLINE void computeDFDt(const velocityGradient_iterator    &velocityGradient_begin,
                                                const velocityGradient_iterator    &velocityGradient_end,
                                                const deformationGradient_iterator &deformationGradient_begin,
                                                const deformationGradient_iterator &deformationGradient_end,
                                                DFDt_iterator DFDt_begin, DFDt_iterator DFDt_end,
                                                dDFDtdL_iterator dDFDtdL_begin, dDFDtdL_iterator dDFDtdL_end,
                                                dDFDtdF_iterator dDFDtdF_begin, dDFDtdF_iterator dDFDtdF_end);

    TARDIGRADE_OPTIONAL_INLINE void computeDFDt(const floatVector &velocityGradient,
                                                const floatVector &deformationGradient, floatVector &DFDt);

    TARDIGRADE_OPTIONAL_INLINE void computeDFDt(const floatVector &velocityGradient,
                                                const floatVector &deformationGradient, floatVector &DFDt,
                                                floatVector &dDFDtdL, floatVector &dDFDtdF);

    TARDIGRADE_OPTIONAL_INLINE void computeDFDt(const floatVector &velocityGradient,
                                                const floatVector &deformationGradient, floatVector &DFDt,
                                                floatMatrix &dDFDtdL, floatMatrix &dDFDtdF);

    template <typename Dt_type, class Ap_iterator, class DApDt_iterator, class DADt_iterator, class dA_iterator,
              class A_iterator, class alpha_iterator>
    TARDIGRADE_OPTIONAL_INLINE void midpointEvolution(const Dt_type &Dt, const Ap_iterator &Ap_begin,
                                                      const Ap_iterator &Ap_end, const DApDt_iterator &DApDt_begin,
                                                      const DApDt_iterator &DApDt_end, const DADt_iterator &DADt_begin,
                                                      const DADt_iterator &DADt_end, dA_iterator dA_begin,
                                                      dA_iterator dA_end, A_iterator A_begin, A_iterator A_end,
                                                      alpha_iterator alpha_begin, alpha_iterator alpha_end);

    template <typename Dt_type, class Ap_iterator, class DApDt_iterator, class DADt_iterator, class dA_iterator,
              class A_iterator, class DADADt_iterator, class alpha_iterator>
    TARDIGRADE_OPTIONAL_INLINE void midpointEvolution(const Dt_type &Dt, const Ap_iterator &Ap_begin,
                                                      const Ap_iterator &Ap_end, const DApDt_iterator &DApDt_begin,
                                                      const DApDt_iterator &DApDt_end, const DADt_iterator &DADt_begin,
                                                      const DADt_iterator &DADt_end, dA_iterator dA_begin,
                                                      dA_iterator dA_end, A_iterator A_begin, A_iterator A_end,
                                                      DADADt_iterator DADADt_begin, DADADt_iterator DADADt_end,
                                                      alpha_iterator alpha_begin, alpha_iterator alpha_end);

    template <typename Dt_type, class Ap_iterator, class DApDt_iterator, class DADt_iterator, class dA_iterator,
              class A_iterator, class DADADt_iterator, class DADApDt_iterator, class alpha_iterator>
    TARDIGRADE_OPTIONAL_INLINE void midpointEvolution(const Dt_type &Dt, const Ap_iterator &Ap_begin,
                                                      const Ap_iterator &Ap_end, const DApDt_iterator &DApDt_begin,
                                                      const DApDt_iterator &DApDt_end, const DADt_iterator &DADt_begin,
                                                      const DADt_iterator &DADt_end, dA_iterator dA_begin,
                                                      dA_iterator dA_end, A_iterator A_begin, A_iterator A_end,
                                                      DADADt_iterator DADADt_begin, DADADt_iterator DADADt_end,
                                                      DADApDt_iterator DADApDt_begin, DADApDt_iterator DADApDt_end,
                                                      alpha_iterator alpha_begin, alpha_iterator alpha_end);

    template <typename Dt_type, class Ap_iterator, class DApDt_iterator, class DADt_iterator, class dA_iterator,
              class A_iterator, typename alpha_type>
    TARDIGRADE_OPTIONAL_INLINE void midpointEvolution(const Dt_type &Dt, const Ap_iterator &Ap_begin,
                                                      const Ap_iterator &Ap_end, const DApDt_iterator &DApDt_begin,
                                                      const DApDt_iterator &DApDt_end, const DADt_iterator &DADt_begin,
                                                      const DADt_iterator &DADt_end, dA_iterator dA_begin,
                                                      dA_iterator dA_end, A_iterator A_begin, A_iterator A_end,
                                                      alpha_type alpha = 0.5);

    template <typename Dt_type, class Ap_iterator, class DApDt_iterator, class DADt_iterator, class dA_iterator,
              class A_iterator, class DADADt_iterator, typename alpha_type>
    TARDIGRADE_OPTIONAL_INLINE void midpointEvolution(const Dt_type &Dt, const Ap_iterator &Ap_begin,
                                                      const Ap_iterator &Ap_end, const DApDt_iterator &DApDt_begin,
                                                      const DApDt_iterator &DApDt_end, const DADt_iterator &DADt_begin,
                                                      const DADt_iterator &DADt_end, dA_iterator dA_begin,
                                                      dA_iterator dA_end, A_iterator A_begin, A_iterator A_end,
                                                      DADADt_iterator DADADt_begin, DADADt_iterator DADADt_end,
                                                      alpha_type alpha = 0.5);

    template <typename Dt_type, class Ap_iterator, class DApDt_iterator, class DADt_iterator, class dA_iterator,
              class A_iterator, class DADADt_iterator, class DADApDt_iterator, typename alpha_type>
    TARDIGRADE_OPTIONAL_INLINE void midpointEvolution(const Dt_type &Dt, const Ap_iterator &Ap_begin,
                                                      const Ap_iterator &Ap_end, const DApDt_iterator &DApDt_begin,
                                                      const DApDt_iterator &DApDt_end, const DADt_iterator &DADt_begin,
                                                      const DADt_iterator &DADt_end, dA_iterator dA_begin,
                                                      dA_iterator dA_end, A_iterator A_begin, A_iterator A_end,
                                                      DADADt_iterator DADADt_begin, DADADt_iterator DADADt_end,
                                                      DADApDt_iterator DADApDt_begin, DADApDt_iterator DADApDt_end,
                                                      alpha_type alpha = 0.5);

    TARDIGRADE_OPTIONAL_INLINE void midpointEvolution(const floatType &Dt, const floatVector &Ap,
                                                      const floatVector &DApDt, const floatVector &DADt,
                                                      floatVector &dA, floatVector &A, const floatVector &alpha);

    TARDIGRADE_OPTIONAL_INLINE void midpointEvolutionFlatJ(const floatType &Dt, const floatVector &Ap,
                                                           const floatVector &DApDt, const floatVector &DADt,
                                                           floatVector &dA, floatVector &A, floatVector &DADADt,
                                                           const floatVector &alpha);

    TARDIGRADE_OPTIONAL_INLINE void midpointEvolutionFlatJ(const floatType &Dt, const floatVector &Ap,
                                                           const floatVector &DApDt, const floatVector &DADt,
                                                           floatVector &dA, floatVector &A, floatVector &DADADt,
                                                           floatVector &DADADtp, const floatVector &alpha);

    TARDIGRADE_OPTIONAL_INLINE void midpointEvolution(const floatType &Dt, const floatVector &Ap,
                                                      const floatVector &DApDt, const floatVector &DADt,
                                                      floatVector &dA, floatVector &A, floatMatrix &DADADt,
                                                      const floatVector &alpha);

    TARDIGRADE_OPTIONAL_INLINE void midpointEvolution(const floatType &Dt, const floatVector &Ap,
                                                      const floatVector &DApDt, const floatVector &DADt,
                                                      floatVector &dA, floatVector &A, floatMatrix &DADADt,
                                                      floatMatrix &DADADtp, const floatVector &alpha);

    TARDIGRADE_OPTIONAL_INLINE void midpointEvolution(const floatType &Dt, const floatVector &Ap,
                                                      const floatVector &DApDt, const floatVector &DADt,
                                                      floatVector &dA, floatVector &A, const floatType alpha = 0.5);

    TARDIGRADE_OPTIONAL_INLINE void midpointEvolutionFlatJ(const floatType &Dt, const floatVector &Ap,
                                                           const floatVector &DApDt, const floatVector &DADt,
                                                           floatVector &dA, floatVector &A, floatVector &DADADt,
                                                           const floatType alpha = 0.5);

    TARDIGRADE_OPTIONAL_INLINE void midpointEvolutionFlatJ(const floatType &Dt, const floatVector &Ap,
                                                           const floatVector &DApDt, const floatVector &DADt,
                                                           floatVector &dA, floatVector &A, floatVector &DADADt,
                                                           floatVector &DADADtp, const floatType alpha = 0.5);

    TARDIGRADE_OPTIONAL_INLINE void midpointEvolution(const floatType &Dt, const floatVector &Ap,
                                                      const floatVector &DApDt, const floatVector &DADt,
                                                      floatVector &dA, floatVector &A, floatMatrix &DADADt,
                                                      const floatType alpha = 0.5);

    TARDIGRADE_OPTIONAL_INLINE void midpointEvolution(const floatType &Dt, const floatVector &Ap,
                                                      const floatVector &DApDt, const floatVector &DADt,
                                                      floatVector &dA, floatVector &A, floatMatrix &DADADt,
                                                      floatMatrix &DADADtp, const floatType alpha = 0.5);

    template <unsigned int dim, unsigned int mode, typename Dt_type, class previousDeformationGradient_iterator,
              class Lp_iterator, class L_iterator, class dF_iterator, class deformationGradient_iterator,
              typename alpha_type>
    TARDIGRADE_OPTIONAL_INLINE void evolveF(
        const Dt_type &Dt, const previousDeformationGradient_iterator &previousDeformationGradient_begin,
        const previousDeformationGradient_iterator &previousDeformationGradient_end, const Lp_iterator &Lp_begin,
        const Lp_iterator &Lp_end, const L_iterator &L_begin, const L_iterator &L_end, dF_iterator dF_begin,
        dF_iterator dF_end, deformationGradient_iterator deformationGradient_begin,
        deformationGradient_iterator deformationGradient_end, const alpha_type alpha);

    template <unsigned int dim, unsigned int mode, typename Dt_type, class previousDeformationGradient_iterator,
              class Lp_iterator, class L_iterator, class dF_iterator, class deformationGradient_iterator,
              class dFdL_iterator, typename alpha_type>
    TARDIGRADE_OPTIONAL_INLINE void evolveF(
        const Dt_type &Dt, const previousDeformationGradient_iterator &previousDeformationGradient_begin,
        const previousDeformationGradient_iterator &previousDeformationGradient_end, const Lp_iterator &Lp_begin,
        const Lp_iterator &Lp_end, const L_iterator &L_begin, const L_iterator &L_end, dF_iterator dF_begin,
        dF_iterator dF_end, deformationGradient_iterator deformationGradient_begin,
        deformationGradient_iterator deformationGradient_end, dFdL_iterator dFdL_begin, dFdL_iterator dFdL_end,
        const alpha_type alpha);

    template <unsigned int dim, unsigned int mode, typename Dt_type, class previousDeformationGradient_iterator,
              class Lp_iterator, class L_iterator, class dF_iterator, class deformationGradient_iterator,
              class dFdL_iterator, class ddFdFp_iterator, class dFdFp_iterator, class dFdLp_iterator,
              typename alpha_type>
    TARDIGRADE_OPTIONAL_INLINE void evolveF(
        const Dt_type &Dt, const previousDeformationGradient_iterator &previousDeformationGradient_begin,
        const previousDeformationGradient_iterator &previousDeformationGradient_end, const Lp_iterator &Lp_begin,
        const Lp_iterator &Lp_end, const L_iterator &L_begin, const L_iterator &L_end, dF_iterator dF_begin,
        dF_iterator dF_end, deformationGradient_iterator deformationGradient_begin,
        deformationGradient_iterator deformationGradient_end, dFdL_iterator dFdL_begin, dFdL_iterator dFdL_end,
        dFdFp_iterator ddFdFp_begin, ddFdFp_iterator ddFdFp_end, dFdFp_iterator dFdFp_begin, dFdFp_iterator dFdFp_end,
        dFdLp_iterator dFdLp_begin, dFdLp_iterator dFdLp_end, const alpha_type alpha);

    TARDIGRADE_OPTIONAL_INLINE void evolveF(const floatType &Dt, const floatVector &previousDeformationGradient,
                                            const floatVector &Lp, const floatVector &L, floatVector &dF,
                                            floatVector &deformationGradient, const floatType alpha = 0.5,
                                            const unsigned int mode = 1);

    TARDIGRADE_OPTIONAL_INLINE void evolveFFlatJ(const floatType &Dt, const floatVector &previousDeformationGradient,
                                                 const floatVector &Lp, const floatVector &L, floatVector &dF,
                                                 floatVector &deformationGradient, floatVector &dFdL,
                                                 const floatType alpha = 0.5, const unsigned int mode = 1);

    TARDIGRADE_OPTIONAL_INLINE void evolveF(const floatType &Dt, const floatVector &previousDeformationGradient,
                                            const floatVector &Lp, const floatVector &L, floatVector &dF,
                                            floatVector &deformationGradient, floatMatrix &dFdL,
                                            const floatType alpha = 0.5, const unsigned int mode = 1);

    TARDIGRADE_OPTIONAL_INLINE void evolveF(const floatType &Dt, const floatVector &previousDeformationGradient,
                                            const floatVector &Lp, const floatVector &L,
                                            floatVector &deformationGradient, floatType alpha = 0.5,
                                            const unsigned int mode = 1);

    TARDIGRADE_OPTIONAL_INLINE void evolveFFlatJ(const floatType &Dt, const floatVector &previousDeformationGradient,
                                                 const floatVector &Lp, const floatVector &L,
                                                 floatVector &deformationGradient, floatVector &dFdL,
                                                 const floatType alpha = 0.5, const unsigned int mode = 1);

    TARDIGRADE_OPTIONAL_INLINE void evolveF(const floatType &Dt, const floatVector &previousDeformationGradient,
                                            const floatVector &Lp, const floatVector &L,
                                            floatVector &deformationGradient, floatMatrix &dFdL,
                                            const floatType alpha = 0.5, const unsigned int mode = 1);

    TARDIGRADE_OPTIONAL_INLINE void evolveFFlatJ(const floatType &Dt, const floatVector &previousDeformationGradient,
                                                 const floatVector &Lp, const floatVector &L,
                                                 floatVector &deformationGradient, floatVector &dFdL,
                                                 floatVector &dFdFp, floatVector &dFdLp, const floatType alpha = 0.5,
                                                 const unsigned int mode = 1);

    TARDIGRADE_OPTIONAL_INLINE void evolveF(const floatType &Dt, const floatVector &previousDeformationGradient,
                                            const floatVector &Lp, const floatVector &L,
                                            floatVector &deformationGradient, floatMatrix &dFdL, floatMatrix &dFdFp,
                                            floatMatrix &dFdLp, const floatType alpha = 0.5,
                                            const unsigned int mode = 1);

    TARDIGRADE_OPTIONAL_INLINE void evolveFFlatJ(const floatType &Dt, const floatVector &previousDeformationGradient,
                                                 const floatVector &Lp, const floatVector &L, floatVector &dF,
                                                 floatVector &deformationGradient, floatVector &dFdL,
                                                 floatVector &ddFdFp, floatVector &dFdFp, floatVector &dFdLp,
                                                 const floatType alpha = 0.5, const unsigned int mode = 1);

    TARDIGRADE_OPTIONAL_INLINE void evolveF(const floatType &Dt, const floatVector &previousDeformationGradient,
                                            const floatVector &Lp, const floatVector &L, floatVector &dF,
                                            floatVector &deformationGradient, floatMatrix &dFdL, floatMatrix &ddFdFp,
                                            floatMatrix &dFdFp, floatMatrix &dFdLp, const floatType alpha = 0.5,
                                            const unsigned int mode = 1);

    template <unsigned int dim, typename Dt_type, class previousDeformationGradient_iterator, class Lp_iterator,
              class L_iterator, class deformationGradient_iterator, typename alpha_type>
    TARDIGRADE_OPTIONAL_INLINE void evolveFExponentialMap(
        const Dt_type &Dt, const previousDeformationGradient_iterator &previousDeformationGradient_begin,
        const previousDeformationGradient_iterator &previousDeformationGradient_end, const Lp_iterator &Lp_begin,
        const Lp_iterator &Lp_end, const L_iterator &L_begin, const L_iterator &L_end,
        deformationGradient_iterator deformationGradient_begin, deformationGradient_iterator deformationGradient_end,
        const alpha_type alpha = 0.5);

    template <unsigned int dim, typename Dt_type, class previousDeformationGradient_iterator, class Lp_iterator,
              class L_iterator, class deformationGradient_iterator, class dFdL_iterator, typename alpha_type>
    TARDIGRADE_OPTIONAL_INLINE void evolveFExponentialMap(
        const Dt_type &Dt, const previousDeformationGradient_iterator &previousDeformationGradient_begin,
        const previousDeformationGradient_iterator &previousDeformationGradient_end, const Lp_iterator &Lp_begin,
        const Lp_iterator &Lp_end, const L_iterator &L_begin, const L_iterator &L_end,
        deformationGradient_iterator deformationGradient_begin, deformationGradient_iterator deformationGradient_end,
        dFdL_iterator dFdL_begin, dFdL_iterator dFdL_end, const alpha_type alpha = 0.5);

    template <unsigned int dim, typename Dt_type, class previousDeformationGradient_iterator, class Lp_iterator,
              class L_iterator, class deformationGradient_iterator, class dFdL_iterator, class dFdFp_iterator,
              class dFdLp_iterator, typename alpha_type>
    TARDIGRADE_OPTIONAL_INLINE void evolveFExponentialMap(
        const Dt_type &Dt, const previousDeformationGradient_iterator &previousDeformationGradient_begin,
        const previousDeformationGradient_iterator &previousDeformationGradient_end, const Lp_iterator &Lp_begin,
        const Lp_iterator &Lp_end, const L_iterator &L_begin, const L_iterator &L_end,
        deformationGradient_iterator deformationGradient_begin, deformationGradient_iterator deformationGradient_end,
        dFdL_iterator dFdL_begin, dFdL_iterator dFdL_end, dFdFp_iterator dFdFp_begin, dFdFp_iterator dFdFp_end,
        dFdLp_iterator dFdLp_begin, dFdLp_iterator dFdLp_end, const alpha_type alpha = 0.5);

    TARDIGRADE_OPTIONAL_INLINE void evolveFExponentialMap(const floatType   &Dt,
                                                          const floatVector &previousDeformationGradient,
                                                          const floatVector &Lp, const floatVector &L,
                                                          floatVector    &deformationGradient,
                                                          const floatType alpha = 0.5);

    TARDIGRADE_OPTIONAL_INLINE void evolveFExponentialMap(const floatType   &Dt,
                                                          const floatVector &previousDeformationGradient,
                                                          const floatVector &Lp, const floatVector &L,
                                                          floatVector &deformationGradient, floatVector &dFdL,
                                                          const floatType alpha = 0.5);

    TARDIGRADE_OPTIONAL_INLINE void evolveFExponentialMap(const floatType   &Dt,
                                                          const floatVector &previousDeformationGradient,
                                                          const floatVector &Lp, const floatVector &L,
                                                          floatVector &deformationGradient, floatVector &dFdL,
                                                          floatVector &dFdFp, floatVector &dFdLp,
                                                          const floatType alpha = 0.5);

    template <typename T>
    TARDIGRADE_OPTIONAL_INLINE T mac(const T &x) {
        /*!
         * Compute the Macaulay brackets of a scalar x
         *
         * returns x if x>0, 0 otherwise
         *
         * \param &x: The incoming scalar.
         */

        return 0.5 * (std::fabs(x) + x);
    }

    template <typename T>
    TARDIGRADE_OPTIONAL_INLINE T mac(const T &x, T &dmacdx) {
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
        if (x >= 0) {
            dmacdx = 1;
        }
        return mac(x);
    }

    template <class A_iterator, class Anorm_iterator>
    TARDIGRADE_OPTIONAL_INLINE void computeUnitNormal(const A_iterator &A_begin, const A_iterator &A_end,
                                                      Anorm_iterator Anorm_begin, Anorm_iterator Anorm_end);

    template <class A_iterator, class Anorm_iterator, class dAnormdA_iterator>
    TARDIGRADE_OPTIONAL_INLINE void computeUnitNormal(const A_iterator &A_begin, const A_iterator &A_end,
                                                      Anorm_iterator Anorm_begin, Anorm_iterator Anorm_end,
                                                      dAnormdA_iterator dAnormdA_begin, dAnormdA_iterator dAnormdA_end);

    TARDIGRADE_OPTIONAL_INLINE void computeUnitNormal(const floatVector &A, floatVector &Anorm);

    TARDIGRADE_OPTIONAL_INLINE void computeUnitNormal(const floatVector &A, floatVector &Anorm, floatVector &dAnormdA);

    TARDIGRADE_OPTIONAL_INLINE void computeUnitNormal(const floatVector &A, floatVector &Anorm, floatMatrix &dAnormdA);

    template <unsigned int dim, class velocityGradient_iterator, class deformationGradient_iterator,
              class pulledBackVelocityGradient_iterator>
    TARDIGRADE_OPTIONAL_INLINE void pullBackVelocityGradient(
        const velocityGradient_iterator &velocityGradient_begin, const velocityGradient_iterator &velocityGradient_end,
        const deformationGradient_iterator &deformationGradient_begin,
        const deformationGradient_iterator &deformationGradient_end,
        pulledBackVelocityGradient_iterator pulledBackVelocityGradient_begin,
        pulledBackVelocityGradient_iterator pulledBackVelocityGradient_end);

    template <unsigned int dim, class velocityGradient_iterator, class deformationGradient_iterator,
              class pulledBackVelocityGradient_iterator, class dPullBackLdL_iterator, class dPullBackLdF_iterator>
    TARDIGRADE_OPTIONAL_INLINE void pullBackVelocityGradient(
        const velocityGradient_iterator &velocityGradient_begin, const velocityGradient_iterator &velocityGradient_end,
        const deformationGradient_iterator &deformationGradient_begin,
        const deformationGradient_iterator &deformationGradient_end,
        pulledBackVelocityGradient_iterator pulledBackVelocityGradient_begin,
        pulledBackVelocityGradient_iterator pulledBackVelocityGradient_end, dPullBackLdL_iterator dPullBackLdL_begin,
        dPullBackLdL_iterator dPullBackLdL_end, dPullBackLdF_iterator dPullBackLdF_begin,
        dPullBackLdF_iterator dPullBackLdF_end);

    TARDIGRADE_OPTIONAL_INLINE void pullBackVelocityGradient(const floatVector &velocityGradient,
                                                             const floatVector &deformationGradient,
                                                             floatVector       &pulledBackVelocityGradient);

    TARDIGRADE_OPTIONAL_INLINE void pullBackVelocityGradient(const floatVector &velocityGradient,
                                                             const floatVector &deformationGradient,
                                                             floatVector       &pulledBackVelocityGradient,
                                                             floatVector &dPullBackLdL, floatVector &dPullBackLdF);

    TARDIGRADE_OPTIONAL_INLINE void pullBackVelocityGradient(const floatVector &velocityGradient,
                                                             const floatVector &deformationGradient,
                                                             floatVector       &pulledBackVelocityGradient,
                                                             floatMatrix &dPullBackLdL, floatMatrix &dPullBackLdF);

    template <typename temperature_type, typename referenceTemperature_type, class linearParameters_iterator,
              class quadraticParameters_iterator, class thermalExpansion_iterator>
    TARDIGRADE_OPTIONAL_INLINE void quadraticThermalExpansion(
        const temperature_type &temperature, const referenceTemperature_type &referenceTemperature,
        const linearParameters_iterator &linearParameters_begin, const linearParameters_iterator &linearParameters_end,
        const quadraticParameters_iterator &quadraticParameters_begin,
        const quadraticParameters_iterator &quadraticParameters_end, thermalExpansion_iterator thermalExpansion_begin,
        thermalExpansion_iterator thermalExpansion_end);

    template <typename temperature_type, typename referenceTemperature_type, class linearParameters_iterator,
              class quadraticParameters_iterator, class thermalExpansion_iterator,
              class thermalExpansionJacobian_iterator>
    TARDIGRADE_OPTIONAL_INLINE void quadraticThermalExpansion(
        const temperature_type &temperature, const referenceTemperature_type &referenceTemperature,
        const linearParameters_iterator &linearParameters_begin, const linearParameters_iterator &linearParameters_end,
        const quadraticParameters_iterator &quadraticParameters_begin,
        const quadraticParameters_iterator &quadraticParameters_end, thermalExpansion_iterator thermalExpansion_begin,
        thermalExpansion_iterator         thermalExpansion_end,
        thermalExpansionJacobian_iterator thermalExpansionJacobian_begin,
        thermalExpansionJacobian_iterator thermalExpansionJacobian_end);

    TARDIGRADE_OPTIONAL_INLINE void quadraticThermalExpansion(const floatType   &temperature,
                                                              const floatType   &referenceTemperature,
                                                              const floatVector &linearParameters,
                                                              const floatVector &quadraticParameters,
                                                              floatVector       &thermalExpansion);

    TARDIGRADE_OPTIONAL_INLINE void quadraticThermalExpansion(
        const floatType &temperature, const floatType &referenceTemperature, const floatVector &linearParameters,
        const floatVector &quadraticParameters, floatVector &thermalExpansion, floatVector &thermalExpansionJacobian);

    template <unsigned int dim, class greenLagrangeStrain_iterator, class deformationGradient_iterator,
              class almansiStrain_iterator>
    TARDIGRADE_OPTIONAL_INLINE void pushForwardGreenLagrangeStrain(
        const greenLagrangeStrain_iterator &greenLagrangeStrain_begin,
        const greenLagrangeStrain_iterator &greenLagrangeStrain_end,
        const deformationGradient_iterator &deformationGradient_begin,
        const deformationGradient_iterator &deformationGradient_end, almansiStrain_iterator almansiStrain_begin,
        almansiStrain_iterator almansiStrain_end);

    template <unsigned int dim, class greenLagrangeStrain_iterator, class deformationGradient_iterator,
              class almansiStrain_iterator, class dAlmansiStraindE_iterator, class dAlmansiStraindF_iterator>
    TARDIGRADE_OPTIONAL_INLINE void pushForwardGreenLagrangeStrain(
        const greenLagrangeStrain_iterator &greenLagrangeStrain_begin,
        const greenLagrangeStrain_iterator &greenLagrangeStrain_end,
        const deformationGradient_iterator &deformationGradient_begin,
        const deformationGradient_iterator &deformationGradient_end, almansiStrain_iterator almansiStrain_begin,
        almansiStrain_iterator almansiStrain_end, dAlmansiStraindE_iterator dAlmansiStraindE_begin,
        dAlmansiStraindE_iterator dAlmansiStraindE_end, dAlmansiStraindF_iterator dAlmansiStraindF_begin,
        dAlmansiStraindF_iterator dAlmansiStraindF_end);

    TARDIGRADE_OPTIONAL_INLINE void pushForwardGreenLagrangeStrain(const floatVector &greenLagrangeStrain,
                                                                   const floatVector &deformationGradient,
                                                                   floatVector       &almansiStrain);

    TARDIGRADE_OPTIONAL_INLINE void pushForwardGreenLagrangeStrain(const floatVector &greenLagrangeStrain,
                                                                   const floatVector &deformationGradient,
                                                                   floatVector       &almansiStrain,
                                                                   floatVector       &dAlmansiStraindE,
                                                                   floatVector       &dAlmansiStraindF);

    TARDIGRADE_OPTIONAL_INLINE void pushForwardGreenLagrangeStrain(const floatVector &greenLagrangeStrain,
                                                                   const floatVector &deformationGradient,
                                                                   floatVector       &almansiStrain,
                                                                   floatMatrix       &dAlmansiStraindE,
                                                                   floatMatrix       &dAlmansiStraindF);

    template <unsigned int dim, class almansiStrain_iterator, class deformationGradient_iterator,
              class greenLagrangeStrain_iterator>
    TARDIGRADE_OPTIONAL_INLINE void pullBackAlmansiStrain(const almansiStrain_iterator       &almansiStrain_begin,
                                                          const almansiStrain_iterator       &almansiStrain_end,
                                                          const deformationGradient_iterator &deformationGradient_begin,
                                                          const deformationGradient_iterator &deformationGradient_end,
                                                          greenLagrangeStrain_iterator        greenLagrangeStrain_begin,
                                                          greenLagrangeStrain_iterator        greenLagrangeStrain_end);

    template <unsigned int dim, class almansiStrain_iterator, class deformationGradient_iterator,
              class greenLagrangeStrain_iterator, class dEde_iterator, class dEdF_iterator>
    TARDIGRADE_OPTIONAL_INLINE void pullBackAlmansiStrain(const almansiStrain_iterator       &almansiStrain_begin,
                                                          const almansiStrain_iterator       &almansiStrain_end,
                                                          const deformationGradient_iterator &deformationGradient_begin,
                                                          const deformationGradient_iterator &deformationGradient_end,
                                                          greenLagrangeStrain_iterator        greenLagrangeStrain_begin,
                                                          greenLagrangeStrain_iterator        greenLagrangeStrain_end,
                                                          dEde_iterator dEde_begin, dEde_iterator dEde_end,
                                                          dEdF_iterator dEdF_begin, dEdF_iterator dEdF_end);

    TARDIGRADE_OPTIONAL_INLINE void pullBackAlmansiStrain(const floatVector &almansiStrain,
                                                          const floatVector &deformationGradient,
                                                          floatVector       &greenLagrangeStrain);

    TARDIGRADE_OPTIONAL_INLINE void pullBackAlmansiStrain(const floatVector &almansiStrain,
                                                          const floatVector &deformationGradient,
                                                          floatVector &greenLagrangeStrain, floatVector &dEde,
                                                          floatVector &dEdF);

    TARDIGRADE_OPTIONAL_INLINE void pullBackAlmansiStrain(const floatVector &almansiStrain,
                                                          const floatVector &deformationGradient,
                                                          floatVector &greenLagrangeStrain, floatMatrix &dEde,
                                                          floatMatrix &dEdF);

    template <unsigned int dim, class A_iterator, class symmA_iterator>
    TARDIGRADE_OPTIONAL_INLINE void computeSymmetricPart(const A_iterator &A_begin, const A_iterator &A_end,
                                                         symmA_iterator symmA_begin, symmA_iterator symmA_end);

    template <unsigned int dim, class A_iterator, class symmA_iterator, class dSymmAdA_iterator>
    TARDIGRADE_OPTIONAL_INLINE void computeSymmetricPart(const A_iterator &A_begin, const A_iterator &A_end,
                                                         symmA_iterator symmA_begin, symmA_iterator symmA_end,
                                                         dSymmAdA_iterator dSymmAdA_begin,
                                                         dSymmAdA_iterator dSymmAdA_end);

    template <class A_iterator, class symmA_iterator>
    TARDIGRADE_OPTIONAL_INLINE void computeSymmetricPart(const A_iterator &A_begin, const A_iterator &A_end,
                                                         symmA_iterator symmA_begin, symmA_iterator symmA_end,
                                                         unsigned int &dim);

    template <class A_iterator, class symmA_iterator, class dSymmAdA_iterator>
    TARDIGRADE_OPTIONAL_INLINE void computeSymmetricPart(const A_iterator &A_begin, const A_iterator &A_end,
                                                         symmA_iterator symmA_begin, symmA_iterator symmA_end,
                                                         dSymmAdA_iterator dSymmAdA_begin,
                                                         dSymmAdA_iterator dSymmAdA_end);

    TARDIGRADE_OPTIONAL_INLINE void computeSymmetricPart(const floatVector &A, floatVector &symmA, unsigned int &dim);

    TARDIGRADE_OPTIONAL_INLINE void computeSymmetricPart(const floatVector &A, floatVector &symmA);

    TARDIGRADE_OPTIONAL_INLINE void computeSymmetricPart(const floatVector &A, floatVector &symmA,
                                                         floatVector &dSymmAdA);

    TARDIGRADE_OPTIONAL_INLINE void computeSymmetricPart(const floatVector &A, floatVector &symmA,
                                                         floatMatrix &dSymmAdA);

    template <unsigned int dim, class PK2_iterator, class F_iterator, class cauchyStress_iterator>
    TARDIGRADE_OPTIONAL_INLINE void pushForwardPK2Stress(const PK2_iterator &PK2_begin, const PK2_iterator &PK2_end,
                                                         const F_iterator &F_begin, const F_iterator &F_end,
                                                         cauchyStress_iterator cauchyStress_begin,
                                                         cauchyStress_iterator cauchyStress_end);

    template <unsigned int dim, class PK2_iterator, class F_iterator, class cauchyStress_iterator,
              class dCauchyStressdPK2_iterator, class dCauchyStressdF_iterator>
    TARDIGRADE_OPTIONAL_INLINE void pushForwardPK2Stress(
        const PK2_iterator &PK2_begin, const PK2_iterator &PK2_end, const F_iterator &F_begin, const F_iterator &F_end,
        cauchyStress_iterator cauchyStress_begin, cauchyStress_iterator cauchyStress_end,
        dCauchyStressdPK2_iterator dCauchyStressdPK2_begin, dCauchyStressdPK2_iterator dCauchyStressdPK2_end,
        dCauchyStressdF_iterator dCauchyStressdF_begin, dCauchyStressdF_iterator dCauchyStressdF_end);

    TARDIGRADE_OPTIONAL_INLINE void pushForwardPK2Stress(const floatVector &PK2, const floatVector &F,
                                                         floatVector &cauchyStress);

    TARDIGRADE_OPTIONAL_INLINE void pushForwardPK2Stress(const floatVector &PK2, const floatVector &F,
                                                         floatVector &cauchyStress, floatVector &dCauchyStressdPK2,
                                                         floatVector &dCauchyStressdF);

    TARDIGRADE_OPTIONAL_INLINE void pushForwardPK2Stress(const floatVector &PK2, const floatVector &F,
                                                         floatVector &cauchyStress, floatMatrix &dCauchyStressdPK2,
                                                         floatMatrix &dCauchyStressdF);

    template <unsigned int dim, class cauchyStress_iterator, class F_iterator, class PK2_iterator>
    TARDIGRADE_OPTIONAL_INLINE void pullBackCauchyStress(const cauchyStress_iterator &cauchyStress_begin,
                                                         const cauchyStress_iterator &cauchyStress_end,
                                                         const F_iterator &F_begin, const F_iterator &F_end,
                                                         PK2_iterator PK2_begin, PK2_iterator PK2_end);

    template <unsigned int dim, class cauchyStress_iterator, class F_iterator, class PK2_iterator,
              class dPK2dCauchyStress_iterator, class dPK2dF_iterator>
    TARDIGRADE_OPTIONAL_INLINE void pullBackCauchyStress(const cauchyStress_iterator &cauchyStress_begin,
                                                         const cauchyStress_iterator &cauchyStress_end,
                                                         const F_iterator &F_begin, const F_iterator &F_end,
                                                         PK2_iterator PK2_begin, PK2_iterator PK2_end,
                                                         dPK2dCauchyStress_iterator dPK2dCauchyStress_begin,
                                                         dPK2dCauchyStress_iterator dPK2dCauchyStress_end,
                                                         dPK2dF_iterator dPK2dF_begin, dPK2dF_iterator dPK2dF_end);

    TARDIGRADE_OPTIONAL_INLINE void pullBackCauchyStress(const floatVector &cauchyStress, const floatVector &F,
                                                         floatVector &PK2);

    TARDIGRADE_OPTIONAL_INLINE void pullBackCauchyStress(const floatVector &cauchyStress, const floatVector &F,
                                                         floatVector &PK2, floatVector &dPK2dCauchyStress,
                                                         floatVector &dPK2dF);

    TARDIGRADE_OPTIONAL_INLINE void pullBackCauchyStress(const floatVector &cauchyStress, const floatVector &F,
                                                         floatVector &PK2, floatMatrix &dPK2dCauchyStress,
                                                         floatMatrix &dPK2dF);

    template <unsigned int dim, class normalVector_iterator, class F_iterator, class dNormalVectordF_iterator>
    TARDIGRADE_OPTIONAL_INLINE void computeDCurrentNormalVectorDF(const normalVector_iterator &normalVector_begin,
                                                                  const normalVector_iterator &normalVector_end,
                                                                  const F_iterator &F_begin, const F_iterator &F_end,
                                                                  dNormalVectordF_iterator dNormalVectordF_begin,
                                                                  dNormalVectordF_iterator dNormalVectordF_end);

    template <unsigned int dim, class normalVector_iterator, class F_iterator,
              class dAreaWeightedNormalVectordF_iterator>
    TARDIGRADE_OPTIONAL_INLINE void computeDCurrentAreaWeightedNormalVectorDF(
        const normalVector_iterator &normalVector_begin, const normalVector_iterator &normalVector_end,
        const F_iterator &F_begin, const F_iterator &F_end,
        dAreaWeightedNormalVectordF_iterator dAreaWeightedNormalVectordF_begin,
        dAreaWeightedNormalVectordF_iterator dAreaWeightedNormalVectordF_end);

    template <unsigned int dim, class normalVector_iterator, class F_iterator, class dCurrentAreadF_iterator>
    TARDIGRADE_OPTIONAL_INLINE void computeDCurrentAreaDF(const normalVector_iterator &normalVector_begin,
                                                          const normalVector_iterator &normalVector_end,
                                                          const F_iterator &F_begin, const F_iterator &F_end,
                                                          dCurrentAreadF_iterator dCurrentAreadF_begin,
                                                          dCurrentAreadF_iterator dCurrentAreadF_end);

    template <unsigned int dim, class normalVector_iterator, class gradU_iterator, class dNormalVectordGradU_iterator>
    TARDIGRADE_OPTIONAL_INLINE void computeDCurrentNormalVectorDGradU(
        const normalVector_iterator &normalVector_begin, const normalVector_iterator &normalVector_end,
        const gradU_iterator &gradU_begin, const gradU_iterator &gradU_end,
        dNormalVectordGradU_iterator dNormalVectordGradU_begin, dNormalVectordGradU_iterator dNormalVectordGradU_end,
        const bool isCurrent = true);

    template <unsigned int dim, class normalVector_iterator, class gradU_iterator,
              class dAreaWeightedNormalVectordGradU_iterator>
    TARDIGRADE_OPTIONAL_INLINE void computeDCurrentAreaWeightedNormalVectorDGradU(
        const normalVector_iterator &normalVector_begin, const normalVector_iterator &normalVector_end,
        const gradU_iterator &gradU_begin, const gradU_iterator &gradU_end,
        dAreaWeightedNormalVectordGradU_iterator dAreaWeightedNormalVectordGradU_begin,
        dAreaWeightedNormalVectordGradU_iterator dAreaWeightedNormalVectordGradU_end, const bool isCurrent = true);

    template <unsigned int dim, class normalVector_iterator, class gradU_iterator, class dCurrentAreadGradU_iterator>
    TARDIGRADE_OPTIONAL_INLINE void computeDCurrentAreaDGradU(const normalVector_iterator &normalVector_begin,
                                                              const normalVector_iterator &normalVector_end,
                                                              const gradU_iterator        &gradU_begin,
                                                              const gradU_iterator        &gradU_end,
                                                              dCurrentAreadGradU_iterator  dCurrentAreadGradU_begin,
                                                              dCurrentAreadGradU_iterator  dCurrentAreadGradU_end,
                                                              const bool                   isCurrent = true);

    TARDIGRADE_OPTIONAL_INLINE void computeDCurrentNormalVectorDF(const floatVector &normalVector, const floatVector &F,
                                                                  floatVector &dNormalVectordF);

    TARDIGRADE_OPTIONAL_INLINE void computeDCurrentAreaWeightedNormalVectorDF(const floatVector &normalVector,
                                                                              const floatVector &F,
                                                                              floatVector &dAreaWeightedNormalVectordF);

    TARDIGRADE_OPTIONAL_INLINE void computeDCurrentAreaDF(const floatVector &normalVector, const floatVector &F,
                                                          floatVector &dCurrentAreadF);

    TARDIGRADE_OPTIONAL_INLINE void computeDCurrentNormalVectorDGradU(const floatVector &normalVector,
                                                                      const floatVector &gradU,
                                                                      floatVector       &dNormalVectordGradU,
                                                                      const bool         isCurrent = true);

    TARDIGRADE_OPTIONAL_INLINE void computeDCurrentAreaWeightedNormalVectorDGradU(
        const floatVector &normalVector, const floatVector &gradU, floatVector &dAreaWeightedNormalVectordGradU,
        const bool isCurrent = true);

    TARDIGRADE_OPTIONAL_INLINE void computeDCurrentAreaDGradU(const floatVector &normalVector, const floatVector &gradU,
                                                              floatVector &dCurrentAreadGradU,
                                                              const bool   isCurrent = true);

}  // namespace tardigradeConstitutiveTools

#ifdef TARDIGRADE_HEADER_ONLY
#include "tardigrade_constitutive_tools.cpp"
#endif

#endif
