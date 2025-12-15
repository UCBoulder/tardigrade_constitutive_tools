/**
 * \file test_tardigrade_constitutive_tools.cpp
 *
 * Tests for tardigrade_constitutive_tools
 */

#include <tardigrade_constitutive_tools.h>

#include <fstream>
#include <iostream>
#include <sstream>

#define BOOST_TEST_MODULE test_tardigrade_constitutive_tools
#include <boost/test/included/unit_test.hpp>
#include <boost/test/tools/output_test_stream.hpp>

#define DEFAULT_TEST_TOLERANCE 1e-6
#define CHECK_PER_ELEMENT boost::test_tools::per_element()

typedef tardigradeConstitutiveTools::floatType   floatType;
typedef tardigradeConstitutiveTools::floatVector floatVector;
typedef tardigradeConstitutiveTools::floatMatrix floatMatrix;

struct cout_redirect {
    cout_redirect(std::streambuf *new_buffer) : old(std::cout.rdbuf(new_buffer)) {}

    ~cout_redirect() { std::cout.rdbuf(old); }

   private:
    std::streambuf *old;
};

struct cerr_redirect {
    cerr_redirect(std::streambuf *new_buffer) : old(std::cerr.rdbuf(new_buffer)) {}

    ~cerr_redirect() { std::cerr.rdbuf(old); }

   private:
    std::streambuf *old;
};

BOOST_AUTO_TEST_CASE(testDeltaDirac, *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    /*!
     * Test the deltaDirac function in constitutive tools
     */

    BOOST_TEST(tardigradeConstitutiveTools::deltaDirac(1, 2) == 0);

    BOOST_TEST(tardigradeConstitutiveTools::deltaDirac(1, 1) == 1);
}

BOOST_AUTO_TEST_CASE(testRotateMatrix, *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    /*!
     * Test the rotation of a matrix by an orthogonal rotation matrix..
     */

    floatVector Q = {-0.44956296, -0.88488713, -0.12193405, -0.37866166, 0.31242661,
                     -0.87120891, 0.80901699,  -0.3454915,  -0.47552826};

    floatVector A = {1, 2, 3, 4, 5, 6, 7, 8, 9};

    floatVector rotatedA_answer = {-0.09485264, -3.38815017, -5.39748037, -1.09823916, 2.23262233,
                                   4.68884658,  -1.68701666, 6.92240128,  12.8622303};

    floatVector rotatedA;
    tardigradeConstitutiveTools::rotateMatrix(A, Q, rotatedA);

    BOOST_TEST(rotatedA == rotatedA_answer, CHECK_PER_ELEMENT);

    // Test rotation back to original frame

    floatVector QT(Q.size(), 0);
    for (unsigned int i = 0; i < 3; i++) {
        for (unsigned int j = 0; j < 3; j++) {
            QT[3 * j + i] = Q[3 * i + j];
        }
    }

    floatVector App;
    tardigradeConstitutiveTools::rotateMatrix(rotatedA, QT, App);

    BOOST_TEST(A == App, CHECK_PER_ELEMENT);

#ifdef TARDIGRADE_HEADER_ONLY

    rotatedA = floatVector(A.size(), 0);
    tardigradeConstitutiveTools::rotateMatrix<3>(std::begin(A), std::end(A), std::begin(Q), std::end(Q),
                                                 std::begin(rotatedA), std::end(rotatedA));

    BOOST_TEST(rotatedA == rotatedA_answer, CHECK_PER_ELEMENT);

#endif
}

BOOST_AUTO_TEST_CASE(testComputeDeformationGradient, *boost::unit_test::tolerance(1e-5)) {
    /*!
     * Test the calculation of the deformation gradient from the displacement gradient
     */

    floatVector gradU = {-0.01078825, -0.0156822, 0.02290497,  -0.00614278, -0.04403221,
                         -0.01019557, 0.02379954, -0.03175083, -0.03245482};

    floatVector answer1 = {0.98994207,  -0.01554229, 0.02211531,  -0.00604919, 0.95820757,
                           -0.00959658, 0.02300559,  -0.02982579, 0.96937029};

    floatVector answer2 = {0.98921175,  -0.0156822, 0.02290497,  -0.00614278, 0.95596779,
                           -0.01019557, 0.02379954, -0.03175083, 0.96754518};

    floatVector F;

    tardigradeConstitutiveTools::computeDeformationGradient(gradU, F, true);

    BOOST_TEST(F == answer1, CHECK_PER_ELEMENT);

    F.clear();
    tardigradeConstitutiveTools::computeDeformationGradient(gradU, F, false);

    BOOST_TEST(F == answer2, CHECK_PER_ELEMENT);

    floatVector dFdGradU_true, dFdGradU_false;

    F.clear();
    tardigradeConstitutiveTools::computeDeformationGradient(gradU, F, dFdGradU_true, true);

    BOOST_TEST(F == answer1, CHECK_PER_ELEMENT);

    F.clear();
    tardigradeConstitutiveTools::computeDeformationGradient(gradU, F, dFdGradU_false, false);

    BOOST_TEST(F == answer2, CHECK_PER_ELEMENT);

    const float eps = 1e-6;

    floatVector j1(9 * 9, 0), j2(9 * 9, 0);

    for (unsigned int i = 0; i < 9; i++) {
        floatVector delta(9, 0);

        delta[i] = eps * std::fabs(gradU[i]) + eps;

        floatVector rp, rm;

        tardigradeConstitutiveTools::computeDeformationGradient(gradU + delta, rp, true);

        tardigradeConstitutiveTools::computeDeformationGradient(gradU - delta, rm, true);

        for (unsigned int j = 0; j < 9; j++) {
            j1[9 * j + i] = (rp[j] - rm[j]) / (2 * delta[i]);
        }

        tardigradeConstitutiveTools::computeDeformationGradient(gradU + delta, rp, false);

        tardigradeConstitutiveTools::computeDeformationGradient(gradU - delta, rm, false);

        for (unsigned int j = 0; j < 9; j++) {
            j2[9 * j + i] = (rp[j] - rm[j]) / (2 * delta[i]);
        }
    }

    BOOST_TEST(j1 == dFdGradU_true, CHECK_PER_ELEMENT);

    BOOST_TEST(j2 == dFdGradU_false, CHECK_PER_ELEMENT);
}

BOOST_AUTO_TEST_CASE(testComputeGreenLagrangeStrain, *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    /*!
     * Test the computation of the Green-Lagrange strain
     */

    floatVector F = {1, 0, 0, 0, 1, 0, 0, 0, 1};
    floatVector E;

    tardigradeConstitutiveTools::computeGreenLagrangeStrain(F, E);

    floatVector E_answer_1(9, 0);

    BOOST_TEST(E == E_answer_1, CHECK_PER_ELEMENT);

    F = {0.69646919, 0.28613933, 0.22685145, 0.55131477, 0.71946897, 0.42310646, 0.98076420, 0.68482974, 0.4809319};

    tardigradeConstitutiveTools::computeGreenLagrangeStrain(F, E);

    floatVector E_answer_2 = {0.37545786, 0.63379879, 0.43147034, 0.63379879, 0.03425154,
                              0.34933978, 0.43147034, 0.34933978, -0.26911192};

    BOOST_TEST(E == E_answer_2, CHECK_PER_ELEMENT);

    floatVector EJ;
    floatMatrix dEdF;
    tardigradeConstitutiveTools::computeGreenLagrangeStrain(F, EJ, dEdF);

    BOOST_TEST(E == EJ, CHECK_PER_ELEMENT);

    floatType eps = 1e-6;
    for (unsigned int i = 0; i < F.size(); i++) {
        floatVector delta(F.size(), 0);

        delta[i] = eps * fabs(F[i]) + eps;

        floatVector EJp, EJm;

        tardigradeConstitutiveTools::computeGreenLagrangeStrain(F + delta, EJp);

        tardigradeConstitutiveTools::computeGreenLagrangeStrain(F - delta, EJm);

        floatVector gradCol = (EJp - EJm) / (2 * delta[i]);

        for (unsigned int j = 0; j < gradCol.size(); j++) {
            BOOST_TEST(gradCol[j] == dEdF[j][i]);
        }
    }
}

BOOST_AUTO_TEST_CASE(testDecomposeGreenLagrangeStrain, *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    /*!
     * Test the decomposition of the Green-Lagrange strain into isochoric and
     * volumetric parts.
     */

    floatVector F = {0.69646919, 0.28613933, 0.22685145, 0.55131477, 0.71946897,
                     0.42310646, 0.98076420, 0.68482974, 0.4809319};

    floatType   J    = tardigradeVectorTools::determinant(F, 3, 3);
    floatVector Fbar = F / pow(J, 1. / 3);

    floatVector E, Ebar;
    tardigradeConstitutiveTools::computeGreenLagrangeStrain(Fbar, Ebar);

    tardigradeConstitutiveTools::computeGreenLagrangeStrain(F, E);

    floatType   JOut;
    floatVector EbarOut;
    tardigradeConstitutiveTools::decomposeGreenLagrangeStrain(E, EbarOut, JOut);

    BOOST_TEST(J == JOut);

    BOOST_TEST(EbarOut == Ebar, CHECK_PER_ELEMENT);

    floatVector EbarOut2;
    floatType   JOut2;
    floatMatrix dEbardE;
    floatVector dJdE;
    tardigradeConstitutiveTools::decomposeGreenLagrangeStrain(E, EbarOut2, JOut2, dEbardE, dJdE);

    BOOST_TEST(EbarOut == EbarOut2, CHECK_PER_ELEMENT);

    BOOST_TEST(JOut == JOut2);

    floatType eps = 1e-8;
    for (unsigned int i = 0; i < E.size(); i++) {
        floatVector delta(E.size(), 0);
        delta[i] = fabs(eps * E[i]);

        floatType Jp, Jm;
        tardigradeConstitutiveTools::decomposeGreenLagrangeStrain(E + delta, EbarOut2, Jp);

        tardigradeConstitutiveTools::decomposeGreenLagrangeStrain(E - delta, EbarOut2, Jm);

        BOOST_TEST((Jp - Jm) / (2 * delta[i]) == dJdE[i]);
    }

    for (unsigned int i = 0; i < E.size(); i++) {
        floatVector delta(E.size(), 0);
        delta[i] = fabs(eps * E[i]);

        floatVector Ebarp, Ebarm;

        tardigradeConstitutiveTools::decomposeGreenLagrangeStrain(E + delta, Ebarp, JOut2);

        tardigradeConstitutiveTools::decomposeGreenLagrangeStrain(E - delta, Ebarm, JOut2);

        floatVector gradCol = (Ebarp - Ebarm) / (2 * delta[i]);

        for (unsigned int j = 0; j < gradCol.size(); j++) {
            BOOST_TEST(gradCol[j] == dEbardE[j][i]);
        }
    }

    floatVector badE = {-1, 0, 0, 0, 1, 0, 0, 0, 1};

    BOOST_REQUIRE_THROW(tardigradeConstitutiveTools::decomposeGreenLagrangeStrain(badE, EbarOut, JOut),
                        std::nested_exception);
}

BOOST_AUTO_TEST_CASE(testMapPK2toCauchy, *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    /*!
     * Test the mapping of the PK2 stress from the reference
     * configuration to the current configuration.
     */

    floatVector F = {1.96469186, -2.13860665, -2.73148546, 0.51314769, 2.1946897,
                     -0.7689354, 4.80764198,  1.84829739,  -0.19068099};

    floatVector PK2 = {-1.07882482, -1.56821984, 2.29049707, -0.61427755, -4.40322103,
                       -1.01955745, 2.37995406,  -3.1750827, -3.24548244};

    floatVector cauchy;

    floatVector cauchy_answer = {-2.47696057, 0.48015011,  -0.28838671, 0.16490963, -0.57481137,
                                 -0.92071407, -0.21450698, -1.22714923, -1.73532173};

    tardigradeConstitutiveTools::mapPK2toCauchy(PK2, F, cauchy);

    BOOST_TEST(cauchy == cauchy_answer, CHECK_PER_ELEMENT);
}

BOOST_AUTO_TEST_CASE(testWLF, *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    /*!
     * Test the computation of the WLF function.
     */

    floatType T  = 145.;
    floatType Tr = 27.5;
    floatType C1 = 18.2;
    floatType C2 = 282.7;

    floatType factor, dfactordT;

    floatVector WLFParameters{Tr, C1, C2};

    tardigradeConstitutiveTools::WLF(T, WLFParameters, factor);

    BOOST_TEST(factor == pow(10, -C1 * (T - Tr) / (C2 + (T - Tr))));

    floatType factor2;
    tardigradeConstitutiveTools::WLF(T, WLFParameters, factor2, dfactordT);

    BOOST_TEST(factor == factor2);

    floatType delta = fabs(1e-6 * T);
    floatType fp, fm;
    tardigradeConstitutiveTools::WLF(T + delta, WLFParameters, fp);
    tardigradeConstitutiveTools::WLF(T - delta, WLFParameters, fm);

    BOOST_TEST(dfactordT == (fp - fm) / (2 * delta));
}

BOOST_AUTO_TEST_CASE(testComputeDGreenLagrangeStrainDF, *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    /*!
     * Test the computation of the gradient of the Green-Lagrange
     * strain w.r.t. the deformation gradient.
     */

    floatVector F = {0.69646919, 0.28613933, 0.22685145, 0.55131477, 0.71946897,
                     0.42310646, 0.98076420, 0.68482974, 0.4809319};

    floatMatrix dEdF;

    tardigradeConstitutiveTools::computeDGreenLagrangeStrainDF(F, dEdF);

    floatVector E, E2;
    tardigradeConstitutiveTools::computeGreenLagrangeStrain(F, E);

    floatType eps = 1e-6;
    for (unsigned int i = 0; i < F.size(); i++) {
        floatVector delta(F.size(), 0);

        delta[i] = eps * fabs(F[i]) + eps;

        floatVector Ep, Em;
        tardigradeConstitutiveTools::computeGreenLagrangeStrain(F + delta, Ep);

        tardigradeConstitutiveTools::computeGreenLagrangeStrain(F - delta, Em);

        floatVector gradCol = (Ep - Em) / (2 * delta[i]);

        for (unsigned int j = 0; j < gradCol.size(); j++) {
            BOOST_TEST(gradCol[j] == dEdF[j][i]);
        }
    }
}

BOOST_AUTO_TEST_CASE(testMidpointEvolution, *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    /*!
     * Test the midpoint evolution algorithm.
     */

    floatType Dt = 2.5;

    floatVector Ap = {9, 10, 11, 12};

    floatVector DApDt = {1, 2, 3, 4};

    floatVector DADt = {5, 6, 7, 8};

    floatVector alphaVec = {0.1, 0.2, 0.3, 0.4};

    floatVector dA, A;

    // Test implicit integration
    tardigradeConstitutiveTools::midpointEvolution(Dt, Ap, DApDt, DADt, dA, A, 0);

    BOOST_TEST(dA == Dt * DADt, CHECK_PER_ELEMENT);

    // Test explicit integration
    tardigradeConstitutiveTools::midpointEvolution(Dt, Ap, DApDt, DADt, dA, A, 1);

    BOOST_TEST(dA == Dt * DApDt, CHECK_PER_ELEMENT);

    // Test midpoint integration
    tardigradeConstitutiveTools::midpointEvolution(Dt, Ap, DApDt, DADt, dA, A);

    BOOST_TEST(A == Ap + Dt * 0.5 * (DApDt + DADt), CHECK_PER_ELEMENT);

    BOOST_TEST(dA == Dt * 0.5 * (DApDt + DADt), CHECK_PER_ELEMENT);

    tardigradeConstitutiveTools::midpointEvolution(Dt, Ap, DApDt, DADt, dA, A, alphaVec);

    floatVector A_answer  = {20.5, 23., 25.5, 28.};
    floatVector dA_answer = {11.5, 13., 14.5, 16.};

    BOOST_TEST(A == A_answer, CHECK_PER_ELEMENT);

    BOOST_TEST(dA == dA_answer, CHECK_PER_ELEMENT);

    // Add test for the jacobian
    floatType eps = 1e-6;

    floatVector A0, Ai, dA0, dAi;

    floatMatrix DADADt;

    tardigradeConstitutiveTools::midpointEvolution(Dt, Ap, DApDt, DADt, dA, A, alphaVec);

    tardigradeConstitutiveTools::midpointEvolution(Dt, Ap, DApDt, DADt, dA0, A0, DADADt, alphaVec);

    BOOST_TEST(A0 == A, CHECK_PER_ELEMENT);

    BOOST_TEST(dA0 == dA, CHECK_PER_ELEMENT);

    for (unsigned int i = 0; i < DADt.size(); i++) {
        floatVector delta = floatVector(DADt.size(), 0);

        delta[i] = eps * (DADt[i]) + eps;

        floatVector Aip, Aim, dAip, dAim;

        tardigradeConstitutiveTools::midpointEvolution(Dt, Ap, DApDt, DADt + delta, dAip, Aip, alphaVec);

        tardigradeConstitutiveTools::midpointEvolution(Dt, Ap, DApDt, DADt - delta, dAim, Aim, alphaVec);

        floatVector gradCol = (Aip - Aim) / (2 * delta[i]);

        for (unsigned int j = 0; j < gradCol.size(); j++) {
            BOOST_TEST(DADADt[j][i] == gradCol[j]);
        }

        gradCol = (dAip - dAim) / (2 * delta[i]);

        for (unsigned int j = 0; j < gradCol.size(); j++) {
            BOOST_TEST(DADADt[j][i] == gradCol[j]);
        }
    }

    floatVector dA1, A1;

    floatMatrix DADADt1, DADADtp;

    tardigradeConstitutiveTools::midpointEvolution(Dt, Ap, DApDt, DADt, dA1, A1, DADADt1, DADADtp, alphaVec);

    BOOST_TEST(A1 == A, CHECK_PER_ELEMENT);

    BOOST_TEST(dA1 == dA, CHECK_PER_ELEMENT);

    BOOST_TEST(tardigradeVectorTools::appendVectors(DADADt1) == tardigradeVectorTools::appendVectors(DADADt),
               CHECK_PER_ELEMENT);

    floatMatrix DADADtp_answer(Ap.size(), floatVector(DApDt.size(), 0));

    for (unsigned int i = 0; i < DApDt.size(); i++) {
        floatVector delta = floatVector(DApDt.size(), 0);

        delta[i] = eps * std::fabs(DApDt[i]) + eps;

        floatVector _dAp, _dAm;

        floatVector _Ap, _Am;

        tardigradeConstitutiveTools::midpointEvolution(Dt, Ap, DApDt + delta, DADt, _dAp, _Ap, alphaVec);

        tardigradeConstitutiveTools::midpointEvolution(Dt, Ap, DApDt - delta, DADt, _dAm, _Am, alphaVec);

        for (unsigned int j = 0; j < Ap.size(); j++) {
            DADADtp_answer[j][i] = (_Ap[j] - _Am[j]) / (2 * delta[i]);

            BOOST_TEST(DADADtp_answer[j][i] == (_dAp[j] - _dAm[j]) / (2 * delta[i]));
        }
    }

    BOOST_TEST(tardigradeVectorTools::appendVectors(DADADtp) == tardigradeVectorTools::appendVectors(DADADtp_answer),
               CHECK_PER_ELEMENT);
}

BOOST_AUTO_TEST_CASE(testMidpointEvolution2, *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    /*!
     * Test the midpoint evolution algorithm.
     */

    floatType Dt = 2.5;

    floatVector Ap = {9, 10, 11, 12};

    floatVector DApDt = {1, 2, 3, 4};

    floatVector DADt = {5, 6, 7, 8};

    floatVector dA, A;

    floatType alpha = 0.67;

    floatVector dA_answer = Dt * alpha * DApDt + Dt * (1 - alpha) * DADt;

    floatVector A_answer = Ap + dA_answer;

    tardigradeConstitutiveTools::midpointEvolution(Dt, Ap, DApDt, DADt, dA, A, alpha);

    BOOST_TEST(dA == dA_answer, CHECK_PER_ELEMENT);

    BOOST_TEST(A == A_answer, CHECK_PER_ELEMENT);

    // Add test for the jacobian
    floatType eps = 1e-6;

    floatVector DADADt, DADApDt;

    dA.clear();
    A.clear();

    tardigradeConstitutiveTools::midpointEvolutionFlatJ(Dt, Ap, DApDt, DADt, dA, A, DADADt, DADApDt, alpha);

    BOOST_TEST(dA == dA_answer, CHECK_PER_ELEMENT);

    BOOST_TEST(A == A_answer, CHECK_PER_ELEMENT);

    {
        constexpr unsigned int VAR_SIZE = 4;
        constexpr unsigned int OUT_SIZE = 4;

        floatVector X(std::begin(DADt), std::end(DADt));

        for (unsigned int i = 0; i < VAR_SIZE; ++i) {
            floatType delta = eps * std::fabs(X[i]) + eps;

            floatVector xp = X;
            floatVector xm = X;

            xp[i] += delta;
            xm[i] -= delta;

            floatVector rp, rm;
            floatVector rp2, rm2;
            tardigradeConstitutiveTools::midpointEvolution(Dt, Ap, DApDt, xp, rp, rp2, alpha);
            tardigradeConstitutiveTools::midpointEvolution(Dt, Ap, DApDt, xm, rm, rm2, alpha);

            for (unsigned int j = 0; j < OUT_SIZE; ++j) {
                BOOST_TEST(DADADt[VAR_SIZE * j + i] == (rp[j] - rm[j]) / (2 * delta));
                BOOST_TEST(DADADt[VAR_SIZE * j + i] == (rp2[j] - rm2[j]) / (2 * delta));
            }
        }
    }

    {
        constexpr unsigned int VAR_SIZE = 4;
        constexpr unsigned int OUT_SIZE = 4;

        floatVector X(std::begin(DApDt), std::end(DApDt));

        for (unsigned int i = 0; i < VAR_SIZE; ++i) {
            floatType delta = eps * std::fabs(X[i]) + eps;

            floatVector xp = X;
            floatVector xm = X;

            xp[i] += delta;
            xm[i] -= delta;

            floatVector rp, rm;
            floatVector rp2, rm2;
            tardigradeConstitutiveTools::midpointEvolution(Dt, Ap, xp, DADt, rp, rp2, alpha);
            tardigradeConstitutiveTools::midpointEvolution(Dt, Ap, xm, DADt, rm, rm2, alpha);

            for (unsigned int j = 0; j < OUT_SIZE; ++j) {
                BOOST_TEST(DADApDt[VAR_SIZE * j + i] == (rp[j] - rm[j]) / (2 * delta));
                BOOST_TEST(DADApDt[VAR_SIZE * j + i] == (rp2[j] - rm2[j]) / (2 * delta));
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(testComputeDFDt, *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    /*!
     * Test the computation of the total time derivative of the
     * deformation gradient.
     */

    floatVector F = {0.69646919, 0.28613933, 0.22685145, 0.55131477, 0.71946897,
                     0.42310646, 0.98076420, 0.68482974, 0.4809319};

    floatVector L = {0.57821272, 0.27720263, 0.45555826, 0.82144027, 0.83961342,
                     0.95322334, 0.4768852,  0.93771539, 0.1056616};

    floatVector answer = {1.00232848, 0.67686793, 0.46754712, 1.96988645, 1.49191786,
                          1.00002629, 0.95274131, 0.88347295, 0.55575157};

    floatVector DFDt;

    tardigradeConstitutiveTools::computeDFDt(L, F, DFDt);

    BOOST_TEST(DFDt == answer, CHECK_PER_ELEMENT);

    // Test the jacobians
    floatVector DFDtJ;
    floatMatrix dDFDtdL, dDFDtdF;
    tardigradeConstitutiveTools::computeDFDt(L, F, DFDtJ, dDFDtdL, dDFDtdF);

    BOOST_TEST(DFDt == DFDtJ, CHECK_PER_ELEMENT);

    // Use finite differences to estimate the jacobian
    floatType eps = 1e-6;
    for (unsigned int i = 0; i < F.size(); i++) {
        // Compute finite difference gradient w.r.t. L
        floatVector delta(L.size(), 0);
        delta[i] = eps * fabs(L[i]) + eps;

        floatVector DFDtp, DFDtm;

        tardigradeConstitutiveTools::computeDFDt(L + delta, F, DFDtp);

        tardigradeConstitutiveTools::computeDFDt(L - delta, F, DFDtm);

        floatVector gradCol = (DFDtp - DFDtm) / (2 * delta[i]);

        for (unsigned int j = 0; j < gradCol.size(); j++) {
            BOOST_TEST(dDFDtdL[j][i] == gradCol[j]);
        }

        // Compute finite difference gradient w.r.t. F
        delta    = floatVector(F.size(), 0);
        delta[i] = eps * fabs(F[i]) + eps;

        tardigradeConstitutiveTools::computeDFDt(L, F + delta, DFDtp);

        gradCol = (DFDtJ - DFDt) / delta[i];

        tardigradeConstitutiveTools::computeDFDt(L, F - delta, DFDtm);

        gradCol = (DFDtp - DFDtm) / (2 * delta[i]);

        for (unsigned int j = 0; j < gradCol.size(); j++) {
            BOOST_TEST(dDFDtdF[j][i] == gradCol[j]);
        }
    }
}

BOOST_AUTO_TEST_CASE(testEvolveF, *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    /*!
     * Test the evolution of the deformation gradient.
     */

    floatType Dt = 2.7;

    floatVector Fp = {0.69646919, 0.28613933, 0.22685145, 0.55131477, 0.71946897,
                      0.42310646, 0.98076420, 0.68482974, 0.4809319};

    floatVector Lp = {0.69006282, 0.0462321,  0.88086378, 0.8153887, 0.54987134,
                      0.72085876, 0.66559485, 0.63708462, 0.54378588};

    floatVector L = {0.57821272, 0.27720263, 0.45555826, 0.82144027, 0.83961342,
                     0.95322334, 0.4768852,  0.93771539, 0.1056616};

    // Test 1 ( mode 1 fully explicit )
    floatVector dF, F;
    tardigradeConstitutiveTools::evolveF(Dt, Fp, Lp, L, F, 1, 1);

    floatVector answer = {4.39551129, 2.53782698, 1.84614498, 4.81201673, 3.75047725,
                          2.48674399, 4.62070491, 3.44211354, 2.32252023};

    BOOST_TEST(answer == F, CHECK_PER_ELEMENT);

    tardigradeConstitutiveTools::evolveF(Dt, Fp, Lp, L, dF, F, 1, 1);

    BOOST_TEST(answer == F, CHECK_PER_ELEMENT);

    BOOST_TEST((answer - Fp) == dF, CHECK_PER_ELEMENT);

    // Test 2 ( mode 1 fully implicit )
    tardigradeConstitutiveTools::evolveF(Dt, Fp, Lp, L, F, 0, 1);

    answer = {0.63522182,  -0.1712192,  -0.00846781, -0.81250979, -0.19375022,
              -0.20193394, -0.36163914, -0.03662069, -0.05769288};

    BOOST_TEST(answer == F, CHECK_PER_ELEMENT);

    tardigradeConstitutiveTools::evolveF(Dt, Fp, Lp, L, dF, F, 0, 1);

    BOOST_TEST(answer == F, CHECK_PER_ELEMENT);

    BOOST_TEST((answer - Fp) == dF, CHECK_PER_ELEMENT);

    // Test 3 ( mode 1 midpoint rule )
    tardigradeConstitutiveTools::evolveF(Dt, Fp, Lp, L, F, 0.5, 1);

    answer = {0.20004929,  -0.4409338,  -0.18955924, -3.59005736, -2.17210401,
              -1.55661536, -1.88391214, -1.13150095, -0.80579654};

    BOOST_TEST(answer == F, CHECK_PER_ELEMENT);

    tardigradeConstitutiveTools::evolveF(Dt, Fp, Lp, L, dF, F, 0.5, 1);

    BOOST_TEST(answer == F, CHECK_PER_ELEMENT);

    BOOST_TEST((answer - Fp) == dF, CHECK_PER_ELEMENT);

    // Tests 4 and 5 ( mode 1 jacobian )
    floatVector dFJ, FJ;
    floatMatrix dFdL;
    tardigradeConstitutiveTools::evolveF(Dt, Fp, Lp, L, dFJ, FJ, dFdL, 0.5, 1);

    BOOST_TEST(answer == FJ, CHECK_PER_ELEMENT);

    floatMatrix dFdL_alt, ddFdFp, dFdFp, dFdLp;

    tardigradeConstitutiveTools::evolveF(Dt, Fp, Lp, L, dFJ, FJ, dFdL_alt, ddFdFp, dFdFp, dFdLp, 0.5, 1);

    BOOST_TEST(answer == FJ, CHECK_PER_ELEMENT);

    BOOST_TEST((answer - Fp) == dFJ, CHECK_PER_ELEMENT);

    BOOST_TEST(tardigradeVectorTools::appendVectors(dFdL) == tardigradeVectorTools::appendVectors(dFdL_alt),
               CHECK_PER_ELEMENT);

    floatType eps = 1e-6;

    floatMatrix dFdL_answer(L.size(), floatVector(L.size(), 0));

    for (unsigned int i = 0; i < L.size(); i++) {
        floatVector delta(L.size(), 0);
        delta[i] = eps * std::fabs(L[i]) + eps;

        floatVector _Fp, _Fm;

        tardigradeConstitutiveTools::evolveF(Dt, Fp, Lp, L + delta, _Fp, 0.5, 1);

        tardigradeConstitutiveTools::evolveF(Dt, Fp, Lp, L - delta, _Fm, 0.5, 1);

        for (unsigned int j = 0; j < L.size(); j++) {
            dFdL_answer[j][i] = (_Fp[j] - _Fm[j]) / (2 * delta[i]);
        }
    }

    BOOST_TEST(tardigradeVectorTools::appendVectors(dFdL) == tardigradeVectorTools::appendVectors(dFdL_answer),
               CHECK_PER_ELEMENT);

    floatMatrix dFdFp_answer(answer.size(), floatVector(Fp.size(), 0));
    floatMatrix ddFdFp_answer(answer.size(), floatVector(Fp.size(), 0));

    for (unsigned int i = 0; i < Fp.size(); i++) {
        floatVector delta(Fp.size(), 0);
        delta[i] = eps * std::fabs(Fp[i]) + eps;

        floatVector _Fp, _Fm;

        tardigradeConstitutiveTools::evolveF(Dt, Fp + delta, Lp, L, _Fp, 0.5, 1);

        tardigradeConstitutiveTools::evolveF(Dt, Fp - delta, Lp, L, _Fm, 0.5, 1);

        for (unsigned int j = 0; j < L.size(); j++) {
            dFdFp_answer[j][i] = (_Fp[j] - _Fm[j]) / (2 * delta[i]);
        }

        floatVector _dFp, _dFm;

        tardigradeConstitutiveTools::evolveF(Dt, Fp + delta, Lp, L, _dFp, _Fp, 0.5, 1);

        tardigradeConstitutiveTools::evolveF(Dt, Fp - delta, Lp, L, _dFm, _Fm, 0.5, 1);

        for (unsigned int j = 0; j < L.size(); j++) {
            ddFdFp_answer[j][i] = (_dFp[j] - _dFm[j]) / (2 * delta[i]);
        }
    }

    BOOST_TEST(tardigradeVectorTools::appendVectors(dFdFp) == tardigradeVectorTools::appendVectors(dFdFp_answer),
               CHECK_PER_ELEMENT);

    BOOST_TEST(tardigradeVectorTools::appendVectors(ddFdFp) == tardigradeVectorTools::appendVectors(ddFdFp_answer),
               CHECK_PER_ELEMENT);

    floatMatrix dFdLp_answer(answer.size(), floatVector(Lp.size(), 0));

    for (unsigned int i = 0; i < Lp.size(); i++) {
        floatVector delta(Lp.size(), 0);
        delta[i] = eps * std::fabs(Lp[i]) + eps;

        floatVector _Fp, _Fm;

        tardigradeConstitutiveTools::evolveF(Dt, Fp, Lp + delta, L, _Fp, 0.5, 1);

        tardigradeConstitutiveTools::evolveF(Dt, Fp, Lp - delta, L, _Fm, 0.5, 1);

        for (unsigned int j = 0; j < L.size(); j++) {
            dFdLp_answer[j][i] = (_Fp[j] - _Fm[j]) / (2 * delta[i]);
        }
    }

    BOOST_TEST(tardigradeVectorTools::appendVectors(dFdLp) == tardigradeVectorTools::appendVectors(dFdLp_answer),
               CHECK_PER_ELEMENT);

    // Test 6 ( mode 2 fully explicit )
    tardigradeConstitutiveTools::evolveF(Dt, Fp, Lp, L, F, 1, 2);

    answer = {3.03173544, 1.1881084,  2.77327313, 3.92282144, 2.58424672,
              3.75584617, 5.18006647, 2.65125419, 4.85252662};

    BOOST_TEST(answer == F, CHECK_PER_ELEMENT);

    tardigradeConstitutiveTools::evolveF(Dt, Fp, Lp, L, dF, F, 1, 2);

    BOOST_TEST((answer - Fp) == dF, CHECK_PER_ELEMENT);

    // Test 7 ( mode 2 fully implicit )
    tardigradeConstitutiveTools::evolveF(Dt, Fp, Lp, L, F, 0, 2);

    answer = {0.65045472,  -0.42475879, -0.09274688, -0.25411831, -0.08867872,
              -0.16467241, 0.45611733,  -0.45427799, -0.17799727};

    BOOST_TEST(answer == F, CHECK_PER_ELEMENT);

    tardigradeConstitutiveTools::evolveF(Dt, Fp, Lp, L, dF, F, 0, 2);

    BOOST_TEST((answer - Fp) == dF, CHECK_PER_ELEMENT);

    // Test 8 ( mode 2 midpoint rule )
    tardigradeConstitutiveTools::evolveF(Dt, Fp, Lp, L, F, 0.5, 2);

    answer = {-0.02066217, -1.43862233, -0.42448874, -0.96426544, -1.72139966,
              -0.83831629, -0.59802055, -2.37943476, -0.88998505};

    BOOST_TEST(answer == F, CHECK_PER_ELEMENT);

    tardigradeConstitutiveTools::evolveF(Dt, Fp, Lp, L, dF, F, 0.5, 2);

    BOOST_TEST(answer == F, CHECK_PER_ELEMENT);

    BOOST_TEST((answer - Fp) == dF, CHECK_PER_ELEMENT);

    // Tests 9 and 10 ( mode 2 jacobian )
    tardigradeConstitutiveTools::evolveF(Dt, Fp, Lp, L, dFJ, FJ, dFdL, 0.5, 2);

    BOOST_TEST(F == FJ, CHECK_PER_ELEMENT);

    tardigradeConstitutiveTools::evolveF(Dt, Fp, Lp, L, dFJ, FJ, dFdL_alt, ddFdFp, dFdFp, dFdLp, 0.5, 2);

    BOOST_TEST(answer == FJ, CHECK_PER_ELEMENT);

    BOOST_TEST((answer - Fp) == dFJ, CHECK_PER_ELEMENT);

    BOOST_TEST(tardigradeVectorTools::appendVectors(dFdL) == tardigradeVectorTools::appendVectors(dFdL_alt),
               CHECK_PER_ELEMENT);

    dFdL_answer = floatMatrix(Fp.size(), floatVector(L.size(), 0));

    for (unsigned int i = 0; i < L.size(); i++) {
        floatVector delta(L.size(), 0);
        delta[i] = eps * std::fabs(L[i]) + eps;

        floatVector _Fp, _Fm;

        tardigradeConstitutiveTools::evolveF(Dt, Fp, Lp, L + delta, _Fp, 0.5, 2);

        tardigradeConstitutiveTools::evolveF(Dt, Fp, Lp, L - delta, _Fm, 0.5, 2);

        for (unsigned int j = 0; j < L.size(); j++) {
            dFdL_answer[j][i] = (_Fp[j] - _Fm[j]) / (2 * delta[i]);
        }
    }

    BOOST_TEST(tardigradeVectorTools::appendVectors(dFdL) == tardigradeVectorTools::appendVectors(dFdL_answer),
               CHECK_PER_ELEMENT);

    dFdFp_answer  = floatMatrix(answer.size(), floatVector(Fp.size(), 0));
    ddFdFp_answer = floatMatrix(answer.size(), floatVector(Fp.size(), 0));

    for (unsigned int i = 0; i < Fp.size(); i++) {
        floatVector delta(Fp.size(), 0);
        delta[i] = eps * std::fabs(Fp[i]) + eps;

        floatVector _Fp, _Fm;

        tardigradeConstitutiveTools::evolveF(Dt, Fp + delta, Lp, L, _Fp, 0.5, 2);

        tardigradeConstitutiveTools::evolveF(Dt, Fp - delta, Lp, L, _Fm, 0.5, 2);

        for (unsigned int j = 0; j < L.size(); j++) {
            dFdFp_answer[j][i] = (_Fp[j] - _Fm[j]) / (2 * delta[i]);
        }

        floatVector _dFp, _dFm;

        tardigradeConstitutiveTools::evolveF(Dt, Fp + delta, Lp, L, _dFp, _Fp, 0.5, 2);

        tardigradeConstitutiveTools::evolveF(Dt, Fp - delta, Lp, L, _dFm, _Fm, 0.5, 2);

        for (unsigned int j = 0; j < L.size(); j++) {
            ddFdFp_answer[j][i] = (_dFp[j] - _dFm[j]) / (2 * delta[i]);
        }
    }

    BOOST_TEST(tardigradeVectorTools::appendVectors(dFdFp) == tardigradeVectorTools::appendVectors(dFdFp_answer),
               CHECK_PER_ELEMENT);

    BOOST_TEST(tardigradeVectorTools::appendVectors(ddFdFp) == tardigradeVectorTools::appendVectors(ddFdFp_answer),
               CHECK_PER_ELEMENT);

    dFdLp_answer = floatMatrix(answer.size(), floatVector(Lp.size(), 0));

    for (unsigned int i = 0; i < Lp.size(); i++) {
        floatVector delta(Lp.size(), 0);
        delta[i] = eps * std::fabs(Lp[i]) + eps;

        floatVector _Fp, _Fm;

        tardigradeConstitutiveTools::evolveF(Dt, Fp, Lp + delta, L, _Fp, 0.5, 2);

        tardigradeConstitutiveTools::evolveF(Dt, Fp, Lp - delta, L, _Fm, 0.5, 2);

        for (unsigned int j = 0; j < L.size(); j++) {
            dFdLp_answer[j][i] = (_Fp[j] - _Fm[j]) / (2 * delta[i]);
        }
    }

    BOOST_TEST(tardigradeVectorTools::appendVectors(dFdLp) == tardigradeVectorTools::appendVectors(dFdLp_answer),
               CHECK_PER_ELEMENT);
}

BOOST_AUTO_TEST_CASE(testMac, *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    /*!
     * Test the computation of the Macullay brackets.
     */

    floatType x = 1;
    BOOST_TEST(tardigradeConstitutiveTools::mac(x) == x);

    x = -1;
    BOOST_TEST(tardigradeConstitutiveTools::mac(x) == 0.);

    floatType xJ = 2;
    floatType dmacdx;
    BOOST_TEST(tardigradeConstitutiveTools::mac(xJ) == tardigradeConstitutiveTools::mac(xJ, dmacdx));

    BOOST_TEST(dmacdx == 1.);

    xJ = -2;
    BOOST_TEST(tardigradeConstitutiveTools::mac(xJ) == tardigradeConstitutiveTools::mac(xJ, dmacdx));

    BOOST_TEST(dmacdx == 0.);
}

BOOST_AUTO_TEST_CASE(testComputeUnitNormal, *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    /*!
     * Test the computation of the unit normal.
     */

    floatVector A = {1, 2, 3, 4, 5, 6, 7, 8, 9};
    floatVector Anorm;

    tardigradeConstitutiveTools::computeUnitNormal(A, Anorm);

    BOOST_CHECK(tardigradeVectorTools::inner(Anorm, Anorm) == 1.);

    // Check the jacobian
    floatVector AnormJ;
    floatMatrix dAnormdA;

    tardigradeConstitutiveTools::computeUnitNormal(A, AnormJ, dAnormdA);

    // Check the normalized value
    BOOST_TEST(AnormJ == Anorm, CHECK_PER_ELEMENT);

    // Check the gradient
    floatType eps = 1e-6;
    for (unsigned int i = 0; i < A.size(); i++) {
        floatVector delta(A.size(), 0);
        delta[i] = eps * fabs(A[i]) + eps;

        floatVector Anormp, Anormm;
        tardigradeConstitutiveTools::computeUnitNormal(A + delta, Anormp);

        tardigradeConstitutiveTools::computeUnitNormal(A - delta, Anormm);

        floatVector gradCol = (Anormp - Anormm) / (2 * delta[i]);

        for (unsigned int j = 0; j < gradCol.size(); j++) {
            BOOST_TEST(dAnormdA[j][i] == gradCol[j]);
        }
    }

    A = {0, 0, 0, 0, 0, 0, 0, 0, 0};

    tardigradeConstitutiveTools::computeUnitNormal(A, Anorm);

    BOOST_TEST(Anorm == A, CHECK_PER_ELEMENT);

    tardigradeConstitutiveTools::computeUnitNormal(A, Anorm, dAnormdA);

    BOOST_TEST(Anorm == A, CHECK_PER_ELEMENT);

    BOOST_CHECK(std::isnan(tardigradeVectorTools::l2norm(dAnormdA)));
}

BOOST_AUTO_TEST_CASE(testPullBackVelocityGradient, *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    /*!
     * Test the pull back operation on the velocity gradient.
     */

    floatVector velocityGradient = {0.69006282, 0.0462321,  0.88086378, 0.8153887, 0.54987134,
                                    0.72085876, 0.66559485, 0.63708462, 0.54378588};

    floatVector deformationGradient = {0.69646919, 0.28613933, 0.22685145, 0.55131477, 0.71946897,
                                       0.42310646, 0.98076420, 0.68482974, 0.4809319};

    floatVector pullBackL;
    floatVector expectedPullBackL = {6.32482111, 3.11877752,   2.43195977,   20.19439192, 10.22175689,
                                     7.88052809, -38.85113898, -18.79212468, -14.76285795};

    tardigradeConstitutiveTools::pullBackVelocityGradient(velocityGradient, deformationGradient, pullBackL);

    BOOST_TEST(pullBackL == expectedPullBackL, CHECK_PER_ELEMENT);

    floatVector pullBackLJ;
    floatMatrix dpbLdL, dpbLdF;

    // Test of the jacobian
    tardigradeConstitutiveTools::pullBackVelocityGradient(velocityGradient, deformationGradient, pullBackLJ, dpbLdL,
                                                          dpbLdF);

    BOOST_TEST(pullBackL == pullBackLJ, CHECK_PER_ELEMENT);

    // Check dpbLdL
    floatType eps = 1e-6;
    for (unsigned int i = 0; i < velocityGradient.size(); i++) {
        floatVector delta(velocityGradient.size(), 0);
        delta[i] = eps * fabs(velocityGradient[i]) + eps;

        floatVector Lp, Lm;

        tardigradeConstitutiveTools::pullBackVelocityGradient(velocityGradient + delta, deformationGradient, Lp);

        tardigradeConstitutiveTools::pullBackVelocityGradient(velocityGradient - delta, deformationGradient, Lm);

        floatVector gradCol = (Lp - Lm) / (2 * delta[i]);

        for (unsigned int j = 0; j < gradCol.size(); j++) {
            BOOST_TEST(gradCol[j] == dpbLdL[j][i]);
        }
    }

    // Check dpbLdF
    for (unsigned int i = 0; i < deformationGradient.size(); i++) {
        floatVector delta(deformationGradient.size(), 0);
        delta[i] = eps * fabs(deformationGradient[i]) + eps;

        floatVector Lp, Lm;

        tardigradeConstitutiveTools::pullBackVelocityGradient(velocityGradient, deformationGradient + delta, Lp);

        tardigradeConstitutiveTools::pullBackVelocityGradient(velocityGradient, deformationGradient - delta, Lm);

        floatVector gradCol = (Lp - Lm) / (2 * delta[i]);

        for (unsigned int j = 0; j < gradCol.size(); j++) {
            BOOST_TEST(gradCol[j] == dpbLdF[j][i]);
        }
    }
}

BOOST_AUTO_TEST_CASE(testQuadraticThermalExpansion, *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    /*!
     * Test the computation of the thermal expansion using a
     * quadratic form.
     */

    floatType temperature          = 283.15;
    floatType referenceTemperature = 273.15;

    floatVector linearParameters    = {1, 2, 3, 4};
    floatVector quadraticParameters = {5, 6, 7, 8};

    floatVector answer = {27825., 33398., 38971., 44544.};

    floatVector thermalExpansion;
    tardigradeConstitutiveTools::quadraticThermalExpansion(temperature, referenceTemperature, linearParameters,
                                                           quadraticParameters, thermalExpansion);

    BOOST_TEST(thermalExpansion == answer, CHECK_PER_ELEMENT);

    floatVector thermalExpansionJ, thermalExpansionJp, thermalExpansionJm, thermalExpansionJacobian;
    floatType   eps   = 1e-6;
    floatType   delta = eps * temperature + eps;

    tardigradeConstitutiveTools::quadraticThermalExpansion(temperature, referenceTemperature, linearParameters,
                                                           quadraticParameters, thermalExpansionJ,
                                                           thermalExpansionJacobian);

    BOOST_TEST(thermalExpansionJ == answer, CHECK_PER_ELEMENT);

    tardigradeConstitutiveTools::quadraticThermalExpansion(temperature + delta, referenceTemperature, linearParameters,
                                                           quadraticParameters, thermalExpansionJp);

    tardigradeConstitutiveTools::quadraticThermalExpansion(temperature - delta, referenceTemperature, linearParameters,
                                                           quadraticParameters, thermalExpansionJm);

    BOOST_TEST(thermalExpansionJacobian == (thermalExpansionJp - thermalExpansionJm) / (2 * delta), CHECK_PER_ELEMENT);
}

BOOST_AUTO_TEST_CASE(testPushForwardGreenLagrangeStrain, *boost::unit_test::tolerance(1e-4)) {
    /*!
     * Test the push-forward operation on the Green-Lagrange strain.
     */

    floatVector deformationGradient = {0.30027935, -0.72811411, 0.26475099, 1.2285819, 0.57663593,
                                       1.43113814, -0.45871432, 0.2175795,  0.54013937};

    floatVector greenLagrangeStrain;
    tardigradeConstitutiveTools::computeGreenLagrangeStrain(deformationGradient, greenLagrangeStrain);

    floatVector almansiStrain = {-0.33393717, 0.0953188,   -0.29053383, 0.0953188,  0.35345526,
                                 0.11588247,  -0.29053383, 0.11588247,  -0.56150741};

    floatVector result;
    tardigradeConstitutiveTools::pushForwardGreenLagrangeStrain(greenLagrangeStrain, deformationGradient, result);

    BOOST_TEST(result == almansiStrain, CHECK_PER_ELEMENT);

    // Test the jacobian
    floatVector resultJ;
    floatMatrix dedE, dedF;
    tardigradeConstitutiveTools::pushForwardGreenLagrangeStrain(greenLagrangeStrain, deformationGradient, resultJ, dedE,
                                                                dedF);

    BOOST_TEST(result == resultJ, CHECK_PER_ELEMENT);

    // Check dedE
    floatType eps = 1e-6;
    for (unsigned int i = 0; i < greenLagrangeStrain.size(); i++) {
        floatVector delta(greenLagrangeStrain.size(), 0);
        delta[i] = eps * fabs(greenLagrangeStrain[i]) + eps;

        floatVector Rp, Rm;

        tardigradeConstitutiveTools::pushForwardGreenLagrangeStrain(greenLagrangeStrain + delta, deformationGradient,
                                                                    Rp);

        tardigradeConstitutiveTools::pushForwardGreenLagrangeStrain(greenLagrangeStrain - delta, deformationGradient,
                                                                    Rm);

        floatVector grad = (Rp - Rm) / (2 * delta[i]);

        for (unsigned int j = 0; j < grad.size(); j++) {
            BOOST_TEST(grad[j] == dedE[j][i]);
        }
    }

    // Check dedF
    for (unsigned int i = 0; i < deformationGradient.size(); i++) {
        floatVector delta(deformationGradient.size(), 0);
        delta[i] = eps * fabs(deformationGradient[i]) + eps;

        floatVector Rp, Rm;

        tardigradeConstitutiveTools::pushForwardGreenLagrangeStrain(greenLagrangeStrain, deformationGradient + delta,
                                                                    Rp);

        tardigradeConstitutiveTools::pushForwardGreenLagrangeStrain(greenLagrangeStrain, deformationGradient - delta,
                                                                    Rm);

        floatVector grad = (Rp - Rm) / (2 * delta[i]);

        for (unsigned int j = 0; j < grad.size(); j++) {
            BOOST_TEST(grad[j] == dedF[j][i]);
        }
    }
}

BOOST_AUTO_TEST_CASE(testPullBackAlmansiStrain, *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    /*!
     * Test the pull-back operation on the Green-Lagrange strain.
     */

    floatVector deformationGradient = {0.1740535,  1.2519364,   -0.9531442,  -0.7512021, -0.60229072,
                                       0.32640812, -0.59754476, -0.06209685, -1.50856757};

    floatVector almansiStrain = {0.25045537, 0.48303426, 0.98555979, 0.51948512, 0.61289453,
                                 0.12062867, 0.8263408,  0.60306013, 0.54506801};

    floatVector answer = {0.55339061, -0.59325289, 0.92984685,  -0.83130342, -0.25274097,
                          -1.5877536, 1.67911302,  -0.83554021, 3.47033811};

    floatVector result;
    tardigradeConstitutiveTools::pullBackAlmansiStrain(almansiStrain, deformationGradient, result);

    BOOST_TEST(answer == result, CHECK_PER_ELEMENT);

    // Test the jacobians
    floatVector resultJ;
    floatMatrix dEde, dEdF;

    tardigradeConstitutiveTools::pullBackAlmansiStrain(almansiStrain, deformationGradient, resultJ, dEde, dEdF);

    BOOST_TEST(answer == resultJ, CHECK_PER_ELEMENT);

    // Testing dEde
    floatType eps = 1e-7;
    for (unsigned int i = 0; i < almansiStrain.size(); i++) {
        floatVector delta(almansiStrain.size(), 0);
        delta[i] = eps * fabs(almansiStrain[i]) + eps;

        floatVector Rp, Rm;

        tardigradeConstitutiveTools::pullBackAlmansiStrain(almansiStrain + delta, deformationGradient, Rp);

        tardigradeConstitutiveTools::pullBackAlmansiStrain(almansiStrain - delta, deformationGradient, Rm);

        floatVector gradCol = (Rp - Rm) / (2 * delta[i]);

        for (unsigned int j = 0; j < gradCol.size(); j++) {
            BOOST_TEST(gradCol[j] == dEde[j][i]);
        }
    }

    // Testing dEdF
    for (unsigned int i = 0; i < deformationGradient.size(); i++) {
        floatVector delta(deformationGradient.size(), 0);
        delta[i] = eps * fabs(deformationGradient[i]) + eps;

        floatVector Rp, Rm;

        tardigradeConstitutiveTools::pullBackAlmansiStrain(almansiStrain, deformationGradient + delta, Rp);

        tardigradeConstitutiveTools::pullBackAlmansiStrain(almansiStrain, deformationGradient - delta, Rm);

        floatVector gradCol = (Rp - Rm) / (2 * delta[i]);

        for (unsigned int j = 0; j < gradCol.size(); j++) {
            BOOST_TEST(gradCol[j] == dEdF[j][i]);
        }
    }
}

BOOST_AUTO_TEST_CASE(testComputeRightCauchyGreen, *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    /*!
     * Test the computation of the Right Cauchy-Green deformation tensor
     */

    floatVector deformationGradient = {1, 2, 3, 4, 5, 6, 7, 8, 9};

    floatVector answer = {66, 78, 90, 78, 93, 108, 90, 108, 126};

    floatVector result;

    tardigradeConstitutiveTools::computeRightCauchyGreen(deformationGradient, result);

    BOOST_TEST(result == answer, CHECK_PER_ELEMENT);

    // Test Jacobian

    floatVector resultJ;
    floatMatrix dCdF;

    tardigradeConstitutiveTools::computeRightCauchyGreen(deformationGradient, resultJ, dCdF);

    BOOST_TEST(resultJ == answer, CHECK_PER_ELEMENT);

    floatType eps = 1e-6;
    for (unsigned int i = 0; i < deformationGradient.size(); i++) {
        floatVector delta(deformationGradient.size(), 0);
        delta[i] = eps * fabs(deformationGradient[i]) + eps;

        floatVector Rp, Rm;

        tardigradeConstitutiveTools::computeRightCauchyGreen(deformationGradient + delta, Rp);

        tardigradeConstitutiveTools::computeRightCauchyGreen(deformationGradient - delta, Rm);

        floatVector gradCol = (Rp - Rm) / (2 * delta[i]);

        for (unsigned int j = 0; j < gradCol.size(); j++) {
            BOOST_TEST(gradCol[j] == dCdF[j][i]);
        }
    }
}

BOOST_AUTO_TEST_CASE(testComputeSymmetricPart, *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    /*!
     * Test the computation of the symmetric part of a matrix
     *
     * \param &results: The output file
     */

    floatVector A = {1, 2, 3, 4, 5, 6, 7, 8, 9};

    floatVector answer = {1., 3., 5., 3., 5., 7., 5., 7., 9.};

    floatVector result;

    tardigradeConstitutiveTools::computeSymmetricPart(A, result);

    BOOST_TEST(result == answer, CHECK_PER_ELEMENT);

    floatVector resultJ;
    floatMatrix dSymmAdA;

    tardigradeConstitutiveTools::computeSymmetricPart(A, resultJ, dSymmAdA);

    BOOST_TEST(result == answer, CHECK_PER_ELEMENT);

    floatVector resultJv;

    floatVector dSymmAdAv;

    tardigradeConstitutiveTools::computeSymmetricPart(A, resultJv, dSymmAdAv);

    BOOST_TEST(resultJv == answer, CHECK_PER_ELEMENT);

    BOOST_TEST(tardigradeVectorTools::appendVectors(dSymmAdA) == dSymmAdAv, CHECK_PER_ELEMENT);

    floatType eps = 1e-6;
    for (unsigned int i = 0; i < A.size(); i++) {
        floatVector delta(A.size(), 0);
        delta[i] = eps * fabs(A[i]) + eps;

        floatVector Rp, Rm;

        tardigradeConstitutiveTools::computeSymmetricPart(A + delta, Rp);

        tardigradeConstitutiveTools::computeSymmetricPart(A - delta, Rm);

        floatVector gradCol = (Rp - Rm) / (2 * delta[i]);

        for (unsigned int j = 0; j < gradCol.size(); j++) {
            BOOST_TEST(gradCol[j] == dSymmAdA[j][i]);
        }
    }
}

#ifdef TARDIGRADE_HEADER_ONLY

BOOST_AUTO_TEST_CASE(testComputeSymmetricPart2, *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    /*!
     * Test the computation of the symmetric part of a matrix
     *
     * \param &results: The output file
     */

    floatVector A = {1, 2, 3, 4, 5, 6, 7, 8, 9};

    floatVector answer = {1., 3., 5., 3., 5., 7., 5., 7., 9.};

    floatVector result(9);

    tardigradeConstitutiveTools::computeSymmetricPart<3>(std::begin(A), std::end(A), std::begin(result),
                                                         std::end(result));

    BOOST_TEST(result == answer, CHECK_PER_ELEMENT);

    floatVector resultJ(9);
    floatVector dSymmAdA(81);

    tardigradeConstitutiveTools::computeSymmetricPart<3>(std::begin(A), std::end(A), std::begin(resultJ),
                                                         std::end(resultJ), std::begin(dSymmAdA), std::end(dSymmAdA));

    BOOST_TEST(result == answer, CHECK_PER_ELEMENT);

    {
        floatType              eps     = 1e-6;
        constexpr unsigned int VAR_DIM = 9;
        constexpr unsigned int OUT_DIM = 9;
        floatVector            X(std::begin(A), std::end(A));

        for (unsigned int i = 0; i < VAR_DIM; ++i) {
            floatType delta = eps * std::fabs(A[i]) + eps;

            floatVector Xp(std::begin(X), std::end(X));
            floatVector Xm(std::begin(X), std::end(X));

            Xp[i] += delta;
            Xm[i] -= delta;

            floatVector Rp(VAR_DIM);
            floatVector Rm(VAR_DIM);

            tardigradeConstitutiveTools::computeSymmetricPart<3>(std::begin(Xp), std::end(Xp), std::begin(Rp),
                                                                 std::end(Rp));

            tardigradeConstitutiveTools::computeSymmetricPart<3>(std::begin(Xm), std::end(Xm), std::begin(Rm),
                                                                 std::end(Rm));

            for (unsigned int j = 0; j < OUT_DIM; ++j) {
                BOOST_TEST(dSymmAdA[VAR_DIM * j + i] == (Rp[j] - Rm[j]) / (2 * delta));
            }
        }
    }
}
#endif

BOOST_AUTO_TEST_CASE(testPushForwardPK2Stress, *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    /*!
     * Test the push forward the PK2 stress to the current configuration
     */

    floatVector PK2 = {1, 2, 3, 4, 5, 6, 7, 8, 9};

    floatVector F = {1, 4, 2, 5, 2, 1, 3, 4, 1};

    floatVector cauchyStressAnswer = {15.16666667, 15.33333333, 16.11111111, 11.33333333, 10.66666667,
                                      11.55555556, 13.66666667, 13.33333333, 14.22222222};

    floatVector result;

    tardigradeConstitutiveTools::pushForwardPK2Stress(PK2, F, result);

    BOOST_TEST(result == cauchyStressAnswer, CHECK_PER_ELEMENT);

    floatVector result2;

    floatMatrix dCauchyStressdPK2, dCauchyStressdF;

    tardigradeConstitutiveTools::pushForwardPK2Stress(PK2, F, result2, dCauchyStressdPK2, dCauchyStressdF);

    BOOST_TEST(result2 == cauchyStressAnswer, CHECK_PER_ELEMENT);

    floatMatrix dCauchyStressdPK2Answer(cauchyStressAnswer.size(), floatVector(PK2.size(), 0));

    floatMatrix dCauchyStressdFAnswer(cauchyStressAnswer.size(), floatVector(F.size(), 0));

    floatType eps = 1e-6;

    for (unsigned int i = 0; i < PK2.size(); i++) {
        floatVector delta(PK2.size(), 0);

        delta[i] = eps * std::fabs(PK2[i]) + eps;

        floatVector cp, cm;

        tardigradeConstitutiveTools::pushForwardPK2Stress(PK2 + delta, F, cp);

        tardigradeConstitutiveTools::pushForwardPK2Stress(PK2 - delta, F, cm);

        for (unsigned int j = 0; j < PK2.size(); j++) {
            dCauchyStressdPK2Answer[j][i] = (cp[j] - cm[j]) / (2 * delta[i]);
        }
    }

    BOOST_TEST(tardigradeVectorTools::appendVectors(dCauchyStressdPK2) ==
                   tardigradeVectorTools::appendVectors(dCauchyStressdPK2Answer),
               CHECK_PER_ELEMENT);

    for (unsigned int i = 0; i < F.size(); i++) {
        floatVector delta(F.size(), 0);

        delta[i] = eps * std::fabs(F[i]) + eps;

        floatVector cp, cm;

        tardigradeConstitutiveTools::pushForwardPK2Stress(PK2, F + delta, cp);

        tardigradeConstitutiveTools::pushForwardPK2Stress(PK2, F - delta, cm);

        for (unsigned int j = 0; j < F.size(); j++) {
            dCauchyStressdFAnswer[j][i] = (cp[j] - cm[j]) / (2 * delta[i]);
        }
    }

    BOOST_TEST(tardigradeVectorTools::appendVectors(dCauchyStressdF) ==
                   tardigradeVectorTools::appendVectors(dCauchyStressdFAnswer),
               CHECK_PER_ELEMENT);
}

BOOST_AUTO_TEST_CASE(testPullBackCauchyStress, *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    floatVector cauchyStress = {0.69646919, 0.28613933, 0.22685145, 0.55131477, 0.71946897,
                                0.42310646, 0.9807642,  0.68482974, 0.4809319};

    floatVector F = {0.39211752, 0.34317802, 0.72904971, 0.43857224, 0.0596779,
                     0.39804426, 0.73799541, 0.18249173, 0.17545176};

    floatVector answer = {0.09712486,  -0.22265266, 0.17348893,  -0.28540852, 1.00439495,
                          -0.24683974, 0.08027086,  -0.27916041, 0.08903243};

    floatVector result, result2;

    tardigradeConstitutiveTools::pullBackCauchyStress(cauchyStress, F, result);

    BOOST_TEST(result == answer, CHECK_PER_ELEMENT);

    floatMatrix dPK2dCauchyStress, dPK2dF;

    floatMatrix dPK2dCauchyStress_answer(answer.size(), floatVector(cauchyStress.size(), 0));

    floatMatrix dPK2dF_answer(answer.size(), floatVector(F.size(), 0));

    floatType eps = 1e-6;

    tardigradeConstitutiveTools::pullBackCauchyStress(cauchyStress, F, result2, dPK2dCauchyStress, dPK2dF);

    BOOST_TEST(result2 == answer, CHECK_PER_ELEMENT);

    for (unsigned int i = 0; i < cauchyStress.size(); i++) {
        floatVector delta(cauchyStress.size(), 0);

        delta[i] = eps * std::fabs(cauchyStress[i]) + eps;

        floatVector resultP, resultM;

        tardigradeConstitutiveTools::pullBackCauchyStress(cauchyStress + delta, F, resultP);

        tardigradeConstitutiveTools::pullBackCauchyStress(cauchyStress - delta, F, resultM);

        for (unsigned int j = 0; j < answer.size(); j++) {
            dPK2dCauchyStress_answer[j][i] = (resultP[j] - resultM[j]) / (2 * delta[i]);
        }
    }

    for (unsigned int i = 0; i < F.size(); i++) {
        floatVector delta(F.size(), 0);

        delta[i] = eps * std::fabs(F[i]) + eps;

        floatVector resultP, resultM;

        tardigradeConstitutiveTools::pullBackCauchyStress(cauchyStress, F + delta, resultP);

        tardigradeConstitutiveTools::pullBackCauchyStress(cauchyStress, F - delta, resultM);

        for (unsigned int j = 0; j < answer.size(); j++) {
            dPK2dF_answer[j][i] = (resultP[j] - resultM[j]) / (2 * delta[i]);
        }
    }

    BOOST_TEST(tardigradeVectorTools::appendVectors(dPK2dCauchyStress) ==
                   tardigradeVectorTools::appendVectors(dPK2dCauchyStress_answer),
               CHECK_PER_ELEMENT);

    BOOST_TEST(tardigradeVectorTools::appendVectors(dPK2dF) == tardigradeVectorTools::appendVectors(dPK2dF_answer),
               CHECK_PER_ELEMENT);
}

BOOST_AUTO_TEST_CASE(testEvolveFExponentialMap, *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    floatType Dt = 2.3;

    floatType alpha = 0.56;

    floatVector Fp = {0.99684486,  -0.00318276, -0.0134401, -0.03494318, 0.97755447,
                      -0.01110235, -0.02224434, 0.01770411, 1.01382113};

    floatVector Lp = {0.21576496, 0.31364397,  -0.45809941, 0.12285551, 0.88064421,
                      0.20391149, -0.47599081, 0.63501654,  0.64909649};

    floatVector L = {-0.39293837, 0.42772133, 0.54629709,  -0.10262954, -0.43893794,
                     0.15378708,  -0.9615284, -0.36965948, 0.0381362};

    floatVector answer = {0.37969552, 0.83224808,  0.49670331,  -0.50354658, 1.28580329,
                          0.63330098, -2.04019687, -0.62835261, 1.70048263};

    floatVector result, resultJ;

    floatVector dFdL, dFdL2;

    floatVector dFdLp, dFdFp;

    tardigradeConstitutiveTools::evolveFExponentialMap(Dt, Fp, Lp, L, result, alpha);

    BOOST_TEST(result == answer, CHECK_PER_ELEMENT);

    tardigradeConstitutiveTools::evolveFExponentialMap(Dt, Fp, Lp, L, resultJ, dFdL, alpha);

    BOOST_TEST(resultJ == answer, CHECK_PER_ELEMENT);

    tardigradeConstitutiveTools::evolveFExponentialMap(Dt, Fp, Lp, L, resultJ, dFdL2, dFdFp, dFdLp, alpha);

    BOOST_TEST(resultJ == answer, CHECK_PER_ELEMENT);

    BOOST_TEST(dFdL2 == dFdL, CHECK_PER_ELEMENT);

    floatVector dFdL_num(81, 0);

    floatVector dFdLp_num(81, 0);

    floatVector dFdFp_num(81, 0);

    floatType eps = 1e-6;

    for (unsigned int i = 0; i < 9; i++) {
        floatType delta = eps * std::fabs(L[i]) + eps;

        floatVector L_p = L;

        floatVector L_m = L;

        L_p[i] += delta;

        L_m[i] -= delta;

        floatVector vp, vm;

        tardigradeConstitutiveTools::evolveFExponentialMap(Dt, Fp, Lp, L_p, vp, alpha);

        tardigradeConstitutiveTools::evolveFExponentialMap(Dt, Fp, Lp, L_m, vm, alpha);

        for (unsigned int j = 0; j < 9; j++) {
            dFdL_num[9 * j + i] = (vp[j] - vm[j]) / (2 * delta);
        }
    }

    BOOST_TEST(dFdL == dFdL_num, CHECK_PER_ELEMENT);

    for (unsigned int i = 0; i < 9; i++) {
        floatType delta = eps * std::fabs(Fp[i]) + eps;

        floatVector Fp_p = Fp;

        floatVector Fp_m = Fp;

        Fp_p[i] += delta;

        Fp_m[i] -= delta;

        floatVector vp, vm;

        tardigradeConstitutiveTools::evolveFExponentialMap(Dt, Fp_p, Lp, L, vp, alpha);

        tardigradeConstitutiveTools::evolveFExponentialMap(Dt, Fp_m, Lp, L, vm, alpha);

        for (unsigned int j = 0; j < 9; j++) {
            dFdFp_num[9 * j + i] = (vp[j] - vm[j]) / (2 * delta);
        }
    }

    BOOST_TEST(dFdFp == dFdFp_num, CHECK_PER_ELEMENT);

    for (unsigned int i = 0; i < 9; i++) {
        floatType delta = eps * std::fabs(Lp[i]) + eps;

        floatVector Lp_p = Lp;

        floatVector Lp_m = Lp;

        Lp_p[i] += delta;

        Lp_m[i] -= delta;

        floatVector vp, vm;

        tardigradeConstitutiveTools::evolveFExponentialMap(Dt, Fp, Lp_p, L, vp, alpha);

        tardigradeConstitutiveTools::evolveFExponentialMap(Dt, Fp, Lp_m, L, vm, alpha);

        for (unsigned int j = 0; j < 9; j++) {
            dFdLp_num[9 * j + i] = (vp[j] - vm[j]) / (2 * delta);
        }
    }

    BOOST_TEST(dFdLp == dFdLp_num, CHECK_PER_ELEMENT);
}

BOOST_AUTO_TEST_CASE(testEvolveFExponentialMap2, *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    floatType Dt = 2.3;

    floatType alpha = 0.56;

    floatVector Fp = {0.99684486,  -0.00318276, -0.0134401, -0.03494318, 0.97755447,
                      -0.01110235, -0.02224434, 0.01770411, 1.01382113};

    floatVector Lp = {0, 0, 0, 0, 0, 0, 0, 0, 0};

    floatVector L = {0, 0, 0, 0, 0, 0, 0, 0, 0};

    floatVector answer = Fp;

    floatVector result, resultJ;

    floatVector dFdL, dFdL2;

    floatVector dFdLp, dFdFp;

    tardigradeConstitutiveTools::evolveFExponentialMap(Dt, Fp, Lp, L, result, alpha);

    BOOST_TEST(result == answer, CHECK_PER_ELEMENT);

    tardigradeConstitutiveTools::evolveFExponentialMap(Dt, Fp, Lp, L, resultJ, dFdL, alpha);

    BOOST_TEST(resultJ == answer, CHECK_PER_ELEMENT);

    tardigradeConstitutiveTools::evolveFExponentialMap(Dt, Fp, Lp, L, resultJ, dFdL2, dFdFp, dFdLp, alpha);

    BOOST_TEST(resultJ == answer, CHECK_PER_ELEMENT);

    BOOST_TEST(dFdL2 == dFdL, CHECK_PER_ELEMENT);

    floatVector dFdL_num(81, 0);

    floatVector dFdLp_num(81, 0);

    floatVector dFdFp_num(81, 0);

    floatType eps = 1e-6;

    for (unsigned int i = 0; i < 9; i++) {
        floatType delta = eps * std::fabs(L[i]) + eps;

        floatVector L_p = L;

        floatVector L_m = L;

        L_p[i] += delta;

        L_m[i] -= delta;

        floatVector vp, vm;

        tardigradeConstitutiveTools::evolveFExponentialMap(Dt, Fp, Lp, L_p, vp, alpha);

        tardigradeConstitutiveTools::evolveFExponentialMap(Dt, Fp, Lp, L_m, vm, alpha);

        for (unsigned int j = 0; j < 9; j++) {
            dFdL_num[9 * j + i] = (vp[j] - vm[j]) / (2 * delta);
        }
    }

    BOOST_TEST(dFdL == dFdL_num, CHECK_PER_ELEMENT);

    for (unsigned int i = 0; i < 9; i++) {
        floatType delta = eps * std::fabs(Fp[i]) + eps;

        floatVector Fp_p = Fp;

        floatVector Fp_m = Fp;

        Fp_p[i] += delta;

        Fp_m[i] -= delta;

        floatVector vp, vm;

        tardigradeConstitutiveTools::evolveFExponentialMap(Dt, Fp_p, Lp, L, vp, alpha);

        tardigradeConstitutiveTools::evolveFExponentialMap(Dt, Fp_m, Lp, L, vm, alpha);

        for (unsigned int j = 0; j < 9; j++) {
            dFdFp_num[9 * j + i] = (vp[j] - vm[j]) / (2 * delta);
        }
    }

    BOOST_TEST(dFdFp == dFdFp_num, CHECK_PER_ELEMENT);

    for (unsigned int i = 0; i < 9; i++) {
        floatType delta = eps * std::fabs(Lp[i]) + eps;

        floatVector Lp_p = Lp;

        floatVector Lp_m = Lp;

        Lp_p[i] += delta;

        Lp_m[i] -= delta;

        floatVector vp, vm;

        tardigradeConstitutiveTools::evolveFExponentialMap(Dt, Fp, Lp_p, L, vp, alpha);

        tardigradeConstitutiveTools::evolveFExponentialMap(Dt, Fp, Lp_m, L, vm, alpha);

        for (unsigned int j = 0; j < 9; j++) {
            dFdLp_num[9 * j + i] = (vp[j] - vm[j]) / (2 * delta);
        }
    }

    BOOST_TEST(dFdLp == dFdLp_num, CHECK_PER_ELEMENT);
}

floatVector normalVectorUtilityFunction(floatVector F, floatVector N) {
    floatVector n(3);

    floatVector invF(9, 0);

    Eigen::Map<const Eigen::Matrix<floatType, 3, 3, Eigen::RowMajor> > F_map(F.data(), 3, 3);

    Eigen::Map<Eigen::Matrix<floatType, 3, 3, Eigen::RowMajor> > invF_map(invF.data(), 3, 3);

    invF_map = F_map.inverse().eval();

    for (unsigned int i = 0; i < 3; i++) {
        for (unsigned int I = 0; I < 3; I++) {
            n[i] += invF[3 * I + i] * N[I];
        }
    }

    n /= tardigradeVectorTools::l2norm(n);

    return n;
}

BOOST_AUTO_TEST_CASE(test_computeDCurrentNormalVectorDF, *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    floatVector N = {1, 2, 3};

    N = N / tardigradeVectorTools::l2norm(N);

    floatVector F = {1.01964692,  -0.02138607, -0.02731485, 0.00513148, 1.0219469,
                     -0.00768935, 0.04807642,  0.01848297,  0.99809319};

    floatVector n = normalVectorUtilityFunction(F, N);

    floatVector dNormalVectordF(27, 0);

    tardigradeConstitutiveTools::computeDCurrentNormalVectorDF(n, F, dNormalVectordF);

    floatType eps = 1e-6;

    floatVector jacobian(27, 0);

    for (unsigned int i = 0; i < 9; i++) {
        floatType delta = eps * std::fabs(F[i]) + eps;

        floatVector Fp = F;

        floatVector Fm = F;

        Fp[i] += delta;

        Fm[i] -= delta;

        floatVector np = normalVectorUtilityFunction(Fp, N);

        floatVector nm = normalVectorUtilityFunction(Fm, N);

        for (unsigned int j = 0; j < 3; j++) {
            jacobian[9 * j + i] = (np[j] - nm[j]) / (2 * delta);
        }
    }

    BOOST_TEST(jacobian == dNormalVectordF, CHECK_PER_ELEMENT);
}

floatVector normalVectorUtilityFunction2(floatVector gradU, floatVector N) {
    floatVector F;

    tardigradeConstitutiveTools::computeDeformationGradient(gradU, F, true);

    floatVector n(3);

    floatVector invF(9, 0);

    Eigen::Map<const Eigen::Matrix<floatType, 3, 3, Eigen::RowMajor> > F_map(F.data(), 3, 3);

    Eigen::Map<Eigen::Matrix<floatType, 3, 3, Eigen::RowMajor> > invF_map(invF.data(), 3, 3);

    invF_map = F_map.inverse().eval();

    for (unsigned int i = 0; i < 3; i++) {
        for (unsigned int I = 0; I < 3; I++) {
            n[i] += invF[3 * I + i] * N[I];
        }
    }

    n /= tardigradeVectorTools::l2norm(n);

    return n;
}

BOOST_AUTO_TEST_CASE(test_computeDCurrentNormalVectorDGradU, *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    floatVector N = {1, 2, 3};

    N = N / tardigradeVectorTools::l2norm(N);

    floatVector gradU = {0.01964692,  -0.02138607, -0.02731485, 0.00513148, 0.0219469,
                         -0.00768935, 0.04807642,  0.01848297,  -0.00190681};

    floatVector n = normalVectorUtilityFunction2(gradU, N);

    floatVector dNormalVectordGradU(27, 0);

    tardigradeConstitutiveTools::computeDCurrentNormalVectorDGradU(n, gradU, dNormalVectordGradU);

    floatType eps = 1e-6;

    floatVector jacobian(27, 0);

    for (unsigned int i = 0; i < 9; i++) {
        floatType delta = eps * std::fabs(gradU[i]) + eps;

        floatVector gradUp = gradU;

        floatVector gradUm = gradU;

        gradUp[i] += delta;

        gradUm[i] -= delta;

        floatVector np = normalVectorUtilityFunction2(gradUp, N);

        floatVector nm = normalVectorUtilityFunction2(gradUm, N);

        for (unsigned int j = 0; j < 3; j++) {
            jacobian[9 * j + i] = (np[j] - nm[j]) / (2 * delta);
        }
    }

    BOOST_TEST(jacobian == dNormalVectordGradU, CHECK_PER_ELEMENT);
}

floatVector areaWeightedNormalVectorUtilityFunction(floatVector F, floatVector N, floatType dA) {
    floatVector nda(3);

    floatVector invF(9, 0);

    Eigen::Map<const Eigen::Matrix<floatType, 3, 3, Eigen::RowMajor> > F_map(F.data(), 3, 3);

    Eigen::Map<Eigen::Matrix<floatType, 3, 3, Eigen::RowMajor> > invF_map(invF.data(), 3, 3);

    invF_map = F_map.inverse().eval();

    floatType J = F_map.determinant();

    for (unsigned int i = 0; i < 3; i++) {
        for (unsigned int I = 0; I < 3; I++) {
            nda[i] += J * invF[3 * I + i] * N[I] * dA;
        }
    }

    return nda;
}

floatVector areaWeightedNormalVectorUtilityFunction2(floatVector gradU, floatVector N, floatType dA) {
    floatVector F;

    tardigradeConstitutiveTools::computeDeformationGradient(gradU, F, true);

    floatVector nda(3);

    floatVector invF(9, 0);

    Eigen::Map<const Eigen::Matrix<floatType, 3, 3, Eigen::RowMajor> > F_map(F.data(), 3, 3);

    Eigen::Map<Eigen::Matrix<floatType, 3, 3, Eigen::RowMajor> > invF_map(invF.data(), 3, 3);

    invF_map = F_map.inverse().eval();

    floatType J = F_map.determinant();

    for (unsigned int i = 0; i < 3; i++) {
        for (unsigned int I = 0; I < 3; I++) {
            nda[i] += J * invF[3 * I + i] * N[I] * dA;
        }
    }

    return nda;
}

floatType areaUtilityFunction(floatVector F, floatVector N, floatType dA) {
    floatVector nda(3);

    floatVector invF(9, 0);

    Eigen::Map<const Eigen::Matrix<floatType, 3, 3, Eigen::RowMajor> > F_map(F.data(), 3, 3);

    Eigen::Map<Eigen::Matrix<floatType, 3, 3, Eigen::RowMajor> > invF_map(invF.data(), 3, 3);

    invF_map = F_map.inverse().eval();

    floatType J = F_map.determinant();

    for (unsigned int i = 0; i < 3; i++) {
        for (unsigned int I = 0; I < 3; I++) {
            nda[i] += J * invF[3 * I + i] * N[I] * dA;
        }
    }

    return tardigradeVectorTools::l2norm(nda);
}

floatType areaUtilityFunction2(floatVector gradU, floatVector N, floatType dA) {
    floatVector F;

    tardigradeConstitutiveTools::computeDeformationGradient(gradU, F, true);

    floatVector nda(3);

    floatVector invF(9, 0);

    Eigen::Map<const Eigen::Matrix<floatType, 3, 3, Eigen::RowMajor> > F_map(F.data(), 3, 3);

    Eigen::Map<Eigen::Matrix<floatType, 3, 3, Eigen::RowMajor> > invF_map(invF.data(), 3, 3);

    invF_map = F_map.inverse().eval();

    floatType J = F_map.determinant();

    for (unsigned int i = 0; i < 3; i++) {
        for (unsigned int I = 0; I < 3; I++) {
            nda[i] += J * invF[3 * I + i] * N[I] * dA;
        }
    }

    return tardigradeVectorTools::l2norm(nda);
}

BOOST_AUTO_TEST_CASE(test_computeDAreaWeightedCurrentNormalVectorDF,
                     *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    floatVector N = {1, 2, 3};

    N = N / tardigradeVectorTools::l2norm(N);

    floatVector F = {1.01964692,  -0.02138607, -0.02731485, 0.00513148, 1.0219469,
                     -0.00768935, 0.04807642,  0.01848297,  0.99809319};

    floatType dA = 2.56;

    floatVector n = normalVectorUtilityFunction(F, N);

    floatVector nda = areaWeightedNormalVectorUtilityFunction(F, N, dA);

    floatType da = tardigradeVectorTools::l2norm(nda);

    BOOST_TEST(n == (nda / da), CHECK_PER_ELEMENT);

    floatVector dAreaWeightedNormalVectordF(27, 0);

    tardigradeConstitutiveTools::computeDCurrentAreaWeightedNormalVectorDF(n, F, dAreaWeightedNormalVectordF);

    floatType eps = 1e-6;

    floatVector jacobian(27, 0);

    for (unsigned int i = 0; i < 9; i++) {
        floatType delta = eps * std::fabs(F[i]) + eps;

        floatVector Fp = F;

        floatVector Fm = F;

        Fp[i] += delta;

        Fm[i] -= delta;

        floatVector np = areaWeightedNormalVectorUtilityFunction(Fp, N, dA);

        floatVector nm = areaWeightedNormalVectorUtilityFunction(Fm, N, dA);

        for (unsigned int j = 0; j < 3; j++) {
            jacobian[9 * j + i] = (np[j] - nm[j]) / (2 * delta);
        }
    }

    BOOST_TEST(jacobian == (dAreaWeightedNormalVectordF * da), CHECK_PER_ELEMENT);
}

BOOST_AUTO_TEST_CASE(test_computeDAreaWeightedCurrentNormalVectorDGradU,
                     *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    floatVector N = {1, 2, 3};

    N = N / tardigradeVectorTools::l2norm(N);

    floatVector gradU = {0.01964692,  -0.02138607, -0.02731485, 0.00513148, 0.0219469,
                         -0.00768935, 0.04807642,  0.01848297,  -0.00190681};

    floatType dA = 2.56;

    floatVector n = normalVectorUtilityFunction2(gradU, N);

    floatVector nda = areaWeightedNormalVectorUtilityFunction2(gradU, N, dA);

    floatType da = tardigradeVectorTools::l2norm(nda);

    BOOST_TEST(n == (nda / da), CHECK_PER_ELEMENT);

    floatVector dAreaWeightedNormalVectordGradU(27, 0);

    tardigradeConstitutiveTools::computeDCurrentAreaWeightedNormalVectorDGradU(n, gradU,
                                                                               dAreaWeightedNormalVectordGradU);

    floatType eps = 1e-6;

    floatVector jacobian(27, 0);

    for (unsigned int i = 0; i < 9; i++) {
        floatType delta = eps * std::fabs(gradU[i]) + eps;

        floatVector gradUp = gradU;

        floatVector gradUm = gradU;

        gradUp[i] += delta;

        gradUm[i] -= delta;

        floatVector np = areaWeightedNormalVectorUtilityFunction2(gradUp, N, dA);

        floatVector nm = areaWeightedNormalVectorUtilityFunction2(gradUm, N, dA);

        for (unsigned int j = 0; j < 3; j++) {
            jacobian[9 * j + i] = (np[j] - nm[j]) / (2 * delta);
        }
    }

    BOOST_TEST(jacobian == (dAreaWeightedNormalVectordGradU * da), CHECK_PER_ELEMENT);
}

BOOST_AUTO_TEST_CASE(test_computeDCurrentAreaDF, *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    floatVector N = {1, 2, 3};

    N = N / tardigradeVectorTools::l2norm(N);

    floatVector F = {1.01964692,  -0.02138607, -0.02731485, 0.00513148, 1.0219469,
                     -0.00768935, 0.04807642,  0.01848297,  0.99809319};

    floatType dA = 2.56;

    floatVector n = normalVectorUtilityFunction(F, N);

    floatVector nda = areaWeightedNormalVectorUtilityFunction(F, N, dA);

    floatType da = tardigradeVectorTools::l2norm(nda);

    floatVector dCurrentAreadF(9, 0);

    tardigradeConstitutiveTools::computeDCurrentAreaDF(n, F, dCurrentAreadF);

    floatType eps = 1e-6;

    floatVector jacobian(9, 0);

    for (unsigned int i = 0; i < 9; i++) {
        floatType delta = eps * std::fabs(F[i]) + eps;

        floatVector Fp = F;

        floatVector Fm = F;

        Fp[i] += delta;

        Fm[i] -= delta;

        floatType np = areaUtilityFunction(Fp, N, dA);

        floatType nm = areaUtilityFunction(Fm, N, dA);

        jacobian[i] = (np - nm) / (2 * delta);
    }

    BOOST_TEST(jacobian == (dCurrentAreadF * da), CHECK_PER_ELEMENT);
}

BOOST_AUTO_TEST_CASE(test_computeDCurrentAreaDGradU, *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    floatVector N = {1, 2, 3};

    N = N / tardigradeVectorTools::l2norm(N);

    floatVector gradU = {0.01964692,  -0.02138607, -0.02731485, 0.00513148, 0.0219469,
                         -0.00768935, 0.04807642,  0.01848297,  -0.00190681};

    floatType dA = 2.56;

    floatVector n = normalVectorUtilityFunction2(gradU, N);

    floatVector nda = areaWeightedNormalVectorUtilityFunction2(gradU, N, dA);

    floatType da = tardigradeVectorTools::l2norm(nda);

    floatVector dCurrentAreadGradU(9, 0);

    tardigradeConstitutiveTools::computeDCurrentAreaDGradU(n, gradU, dCurrentAreadGradU);

    floatType eps = 1e-6;

    floatVector jacobian(9, 0);

    for (unsigned int i = 0; i < 9; i++) {
        floatType delta = eps * std::fabs(gradU[i]) + eps;

        floatVector gradUp = gradU;

        floatVector gradUm = gradU;

        gradUp[i] += delta;

        gradUm[i] -= delta;

        floatType np = areaUtilityFunction2(gradUp, N, dA);

        floatType nm = areaUtilityFunction2(gradUm, N, dA);

        jacobian[i] = (np - nm) / (2 * delta);
    }

    BOOST_TEST(jacobian == (dCurrentAreadGradU * da), CHECK_PER_ELEMENT);
}
