/*
 * Copyright 2013 Daniel Warner <contact@danrw.com>
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include "SGP4.h"

#include "Util.h"
#include "Vector.h"
#include "SatelliteException.h"
#include "DecayedException.h"

#include <cmath>
#include <iomanip>
#include <cstring>

namespace libsgp4
{

    void SGP4::SetTle(const Tle &tle)
    {
        /*
         * extract and format tle data
         */
        elements_ = OrbitalElements(tle);

        Initialise();
    }

    void SGP4::Initialise()
    {
        /*
         * reset all constants etc
         */
        Reset();

        /*
         * error checks
         */
        if (elements_.Eccentricity() < 0.0 || elements_.Eccentricity() > 0.999)
        {
            return;
        }

        if (elements_.Inclination() < 0.0 || elements_.Inclination() > kPI)
        {
            return;
        }

        RecomputeConstants(elements_.Inclination(),
                           common_consts_.sinio,
                           common_consts_.cosio,
                           common_consts_.x3thm1,
                           common_consts_.x1mth2,
                           common_consts_.x7thm1,
                           common_consts_.xlcof,
                           common_consts_.aycof);

        const double theta2 = common_consts_.cosio * common_consts_.cosio;
        const double eosq = elements_.Eccentricity() * elements_.Eccentricity();
        const double betao2 = 1.0 - eosq;
        const double betao = sqrt(betao2);

        /*
         * for perigee below 156km, the values of
         * s4 and qoms2t are altered
         */
        double s4 = kS;
        double qoms24 = kQOMS2T;
        if (elements_.Perigee() < 156.0)
        {
            s4 = elements_.Perigee() - 78.0;
            if (elements_.Perigee() < 98.0)
            {
                s4 = 20.0;
            }
            qoms24 = pow((120.0 - s4) * kAE / kXKMPER, 4.0);
            s4 = s4 / kXKMPER + kAE;
        }

        /*
         * generate constants
         */
        const double pinvsq = 1.0 / (elements_.RecoveredSemiMajorAxis() * elements_.RecoveredSemiMajorAxis() * betao2 * betao2);
        const double tsi = 1.0 / (elements_.RecoveredSemiMajorAxis() - s4);
        common_consts_.eta = elements_.RecoveredSemiMajorAxis() * elements_.Eccentricity() * tsi;
        const double etasq = common_consts_.eta * common_consts_.eta;
        const double eeta = elements_.Eccentricity() * common_consts_.eta;
        const double psisq = fabs(1.0 - etasq);
        const double coef = qoms24 * pow(tsi, 4.0);
        const double coef1 = coef / pow(psisq, 3.5);
        const double c2 = coef1 * elements_.RecoveredMeanMotion() * (elements_.RecoveredSemiMajorAxis() * (1.0 + 1.5 * etasq + eeta * (4.0 + etasq)) + 0.75 * kCK2 * tsi / psisq * common_consts_.x3thm1 * (8.0 + 3.0 * etasq * (8.0 + etasq)));
        common_consts_.c1 = elements_.BStar() * c2;
        common_consts_.c4 = 2.0 * elements_.RecoveredMeanMotion() * coef1 * elements_.RecoveredSemiMajorAxis() * betao2 * (common_consts_.eta * (2.0 + 0.5 * etasq) + elements_.Eccentricity() * (0.5 + 2.0 * etasq) - 2.0 * kCK2 * tsi / (elements_.RecoveredSemiMajorAxis() * psisq) * (-3.0 * common_consts_.x3thm1 * (1.0 - 2.0 * eeta + etasq * (1.5 - 0.5 * eeta)) + 0.75 * common_consts_.x1mth2 * (2.0 * etasq - eeta * (1.0 + etasq)) * cos(2.0 * elements_.ArgumentPerigee())));
        const double theta4 = theta2 * theta2;
        const double temp1 = 3.0 * kCK2 * pinvsq * elements_.RecoveredMeanMotion();
        const double temp2 = temp1 * kCK2 * pinvsq;
        const double temp3 = 1.25 * kCK4 * pinvsq * pinvsq * elements_.RecoveredMeanMotion();
        common_consts_.xmdot = elements_.RecoveredMeanMotion() + 0.5 * temp1 * betao * common_consts_.x3thm1 + 0.0625 * temp2 * betao * (13.0 - 78.0 * theta2 + 137.0 * theta4);
        const double x1m5th = 1.0 - 5.0 * theta2;
        common_consts_.omgdot = -0.5 * temp1 * x1m5th +
                                0.0625 * temp2 * (7.0 - 114.0 * theta2 + 395.0 * theta4) +
                                temp3 * (3.0 - 36.0 * theta2 + 49.0 * theta4);
        const double xhdot1 = -temp1 * common_consts_.cosio;
        common_consts_.xnodot = xhdot1 + (0.5 * temp2 * (4.0 - 19.0 * theta2) + 2.0 * temp3 *
                                                                                    (3.0 - 7.0 * theta2)) *
                                             common_consts_.cosio;
        common_consts_.xnodcf = 3.5 * betao2 * xhdot1 * common_consts_.c1;
        common_consts_.t2cof = 1.5 * common_consts_.c1;

        double c3 = 0.0;
        if (elements_.Eccentricity() > 1.0e-4)
        {
            c3 = coef * tsi * kA3OVK2 * elements_.RecoveredMeanMotion() * kAE *
                 common_consts_.sinio / elements_.Eccentricity();
        }

        nearspace_consts_.c5 = 2.0 * coef1 * elements_.RecoveredSemiMajorAxis() * betao2 * (1.0 + 2.75 * (etasq + eeta) + eeta * etasq);
        nearspace_consts_.omgcof = elements_.BStar() * c3 * cos(elements_.ArgumentPerigee());

        nearspace_consts_.xmcof = 0.0;
        if (elements_.Eccentricity() > 1.0e-4)
        {
            nearspace_consts_.xmcof = -kTWOTHIRD * coef * elements_.BStar() * kAE / eeta;
        }

        nearspace_consts_.delmo = pow(1.0 + common_consts_.eta * (cos(elements_.MeanAnomoly())), 3.0);
        nearspace_consts_.sinmo = sin(elements_.MeanAnomoly());
        const double c1sq = common_consts_.c1 * common_consts_.c1;
        nearspace_consts_.d2 = 4.0 * elements_.RecoveredSemiMajorAxis() * tsi * c1sq;
        const double temp = nearspace_consts_.d2 * tsi * common_consts_.c1 / 3.0;
        nearspace_consts_.d3 = (17.0 * elements_.RecoveredSemiMajorAxis() + s4) * temp;
        nearspace_consts_.d4 = 0.5 * temp * elements_.RecoveredSemiMajorAxis() *
                               tsi * (221.0 * elements_.RecoveredSemiMajorAxis() + 31.0 * s4) * common_consts_.c1;
        nearspace_consts_.t3cof = nearspace_consts_.d2 + 2.0 * c1sq;
        nearspace_consts_.t4cof = 0.25 * (3.0 * nearspace_consts_.d3 + common_consts_.c1 *
                                                                           (12.0 * nearspace_consts_.d2 + 10.0 * c1sq));
        nearspace_consts_.t5cof = 0.2 * (3.0 * nearspace_consts_.d4 + 12.0 * common_consts_.c1 * nearspace_consts_.d3 + 6.0 * nearspace_consts_.d2 * nearspace_consts_.d2 + 15.0 * c1sq * (2.0 * nearspace_consts_.d2 + c1sq));
    }

    Eci SGP4::FindPosition(const DateTime &dt) const
    {
        return FindPosition((dt - elements_.Epoch()).TotalMinutes());
    }

    Eci SGP4::FindPosition(double tsince) const
    {
        return FindPositionSGP4(tsince);
    }

    void SGP4::RecomputeConstants(const double xinc,
                                  double &sinio,
                                  double &cosio,
                                  double &x3thm1,
                                  double &x1mth2,
                                  double &x7thm1,
                                  double &xlcof,
                                  double &aycof)
    {
        sinio = sin(xinc);
        cosio = cos(xinc);

        const double theta2 = cosio * cosio;

        x3thm1 = 3.0 * theta2 - 1.0;
        x1mth2 = 1.0 - theta2;
        x7thm1 = 7.0 * theta2 - 1.0;

        if (fabs(cosio + 1.0) > 1.5e-12)
        {
            xlcof = 0.125 * kA3OVK2 * sinio * (3.0 + 5.0 * cosio) / (1.0 + cosio);
        }
        else
        {
            xlcof = 0.125 * kA3OVK2 * sinio * (3.0 + 5.0 * cosio) / 1.5e-12;
        }

        aycof = 0.25 * kA3OVK2 * sinio;
    }

    Eci SGP4::FindPositionSGP4(double tsince) const
    {
        /*
         * the final values
         */
        double e;
        double a;
        double omega;
        double xl;
        double xnode;
        const double xinc = elements_.Inclination();

        /*
         * update for secular gravity and atmospheric drag
         */
        const double xmdf = elements_.MeanAnomoly() + common_consts_.xmdot * tsince;
        const double omgadf = elements_.ArgumentPerigee() + common_consts_.omgdot * tsince;
        const double xnoddf = elements_.AscendingNode() + common_consts_.xnodot * tsince;

        omega = omgadf;
        double xmp = xmdf;

        const double tsq = tsince * tsince;
        xnode = xnoddf + common_consts_.xnodcf * tsq;
        double tempa = 1.0 - common_consts_.c1 * tsince;
        double tempe = elements_.BStar() * common_consts_.c4 * tsince;
        double templ = common_consts_.t2cof * tsq;

        const double delomg = nearspace_consts_.omgcof * tsince;
        const double delm = nearspace_consts_.xmcof * (pow(1.0 + common_consts_.eta * cos(xmdf), 3.0) - nearspace_consts_.delmo);
        const double temp = delomg + delm;

        xmp += temp;
        omega -= temp;

        const double tcube = tsq * tsince;
        const double tfour = tsince * tcube;

        tempa = tempa - nearspace_consts_.d2 * tsq - nearspace_consts_.d3 * tcube - nearspace_consts_.d4 * tfour;
        tempe += elements_.BStar() * nearspace_consts_.c5 * (sin(xmp) - nearspace_consts_.sinmo);
        templ += nearspace_consts_.t3cof * tcube + tfour * (nearspace_consts_.t4cof + tsince * nearspace_consts_.t5cof);

        a = elements_.RecoveredSemiMajorAxis() * tempa * tempa;
        e = elements_.Eccentricity() - tempe;
        xl = xmp + omega + xnode + elements_.RecoveredMeanMotion() * templ;

        /*
         * fix tolerance for error recognition
         */
        if (e < 1.0e-6)
        {
            e = 1.0e-6;
        }
        else if (e > (1.0 - 1.0e-6))
        {
            e = 1.0 - 1.0e-6;
        }

        /*
         * using calculated values, find position and velocity
         * we can pass in constants from Initialise() as these dont change
         */
        return CalculateFinalPositionVelocity(elements_.Epoch().AddMinutes(tsince),
                                              e,
                                              a,
                                              omega,
                                              xl,
                                              xnode,
                                              xinc,
                                              common_consts_.xlcof,
                                              common_consts_.aycof,
                                              common_consts_.x3thm1,
                                              common_consts_.x1mth2,
                                              common_consts_.x7thm1,
                                              common_consts_.cosio,
                                              common_consts_.sinio);
    }

    Eci SGP4::CalculateFinalPositionVelocity(
        const DateTime &dt,
        const double e,
        const double a,
        const double omega,
        const double xl,
        const double xnode,
        const double xinc,
        const double xlcof,
        const double aycof,
        const double x3thm1,
        const double x1mth2,
        const double x7thm1,
        const double cosio,
        const double sinio)
    {
        const double beta2 = 1.0 - e * e;
        const double xn = kXKE / pow(a, 1.5);
        /*
         * long period periodics
         */
        const double axn = e * cos(omega);
        const double temp11 = 1.0 / (a * beta2);
        const double xll = temp11 * xlcof * axn;
        const double aynl = temp11 * aycof;
        const double xlt = xl + xll;
        const double ayn = e * sin(omega) + aynl;
        const double elsq = axn * axn + ayn * ayn;

        /*
         * solve keplers equation
         * - solve using Newton-Raphson root solving
         * - here capu is almost the mean anomoly
         * - initialise the eccentric anomaly term epw
         * - The fmod saves reduction of angle to +/-2pi in sin/cos() and prevents
         * convergence problems.
         */
        const double capu = fmod(xlt - xnode, kTWOPI);
        double epw = capu;

        double sinepw = 0.0;
        double cosepw = 0.0;
        double ecose = 0.0;
        double esine = 0.0;

        /*
         * sensibility check for N-R correction
         */
        const double max_newton_naphson = 1.25 * fabs(sqrt(elsq));

        bool kepler_running = true;

        for (int i = 0; i < 10 && kepler_running; i++)
        {
            sinepw = sin(epw);
            cosepw = cos(epw);
            ecose = axn * cosepw + ayn * sinepw;
            esine = axn * sinepw - ayn * cosepw;

            double f = capu - epw + esine;

            if (fabs(f) < 1.0e-12)
            {
                kepler_running = false;
            }
            else
            {
                /*
                 * 1st order Newton-Raphson correction
                 */
                const double fdot = 1.0 - ecose;
                double delta_epw = f / fdot;

                /*
                 * 2nd order Newton-Raphson correction.
                 * f / (fdot - 0.5 * d2f * f/fdot)
                 */
                if (i == 0)
                {
                    if (delta_epw > max_newton_naphson)
                    {
                        delta_epw = max_newton_naphson;
                    }
                    else if (delta_epw < -max_newton_naphson)
                    {
                        delta_epw = -max_newton_naphson;
                    }
                }
                else
                {
                    delta_epw = f / (fdot + 0.5 * esine * delta_epw);
                }

                /*
                 * Newton-Raphson correction of -F/DF
                 */
                epw += delta_epw;
            }
        }
        /*
         * short period preliminary quantities
         */
        const double temp21 = 1.0 - elsq;
        const double pl = a * temp21;

        const double r = a * (1.0 - ecose);
        const double temp31 = 1.0 / r;
        const double rdot = kXKE * sqrt(a) * esine * temp31;
        const double rfdot = kXKE * sqrt(pl) * temp31;
        const double temp32 = a * temp31;
        const double betal = sqrt(temp21);
        const double temp33 = 1.0 / (1.0 + betal);
        const double cosu = temp32 * (cosepw - axn + ayn * esine * temp33);
        const double sinu = temp32 * (sinepw - ayn - axn * esine * temp33);
        const double u = atan2(sinu, cosu);
        const double sin2u = 2.0 * sinu * cosu;
        const double cos2u = 2.0 * cosu * cosu - 1.0;

        /*
         * update for short periodics
         */
        const double temp41 = 1.0 / pl;
        const double temp42 = kCK2 * temp41;
        const double temp43 = temp42 * temp41;

        const double rk = r * (1.0 - 1.5 * temp43 * betal * x3thm1) + 0.5 * temp42 * x1mth2 * cos2u;
        const double uk = u - 0.25 * temp43 * x7thm1 * sin2u;
        const double xnodek = xnode + 1.5 * temp43 * cosio * sin2u;
        const double xinck = xinc + 1.5 * temp43 * cosio * sinio * cos2u;
        const double rdotk = rdot - xn * temp42 * x1mth2 * sin2u;
        const double rfdotk = rfdot + xn * temp42 * (x1mth2 * cos2u + 1.5 * x3thm1);

        /*
         * orientation vectors
         */
        const double sinuk = sin(uk);
        const double cosuk = cos(uk);
        const double sinik = sin(xinck);
        const double cosik = cos(xinck);
        const double sinnok = sin(xnodek);
        const double cosnok = cos(xnodek);
        const double xmx = -sinnok * cosik;
        const double xmy = cosnok * cosik;
        const double ux = xmx * sinuk + cosnok * cosuk;
        const double uy = xmy * sinuk + sinnok * cosuk;
        const double uz = sinik * sinuk;
        const double vx = xmx * cosuk - cosnok * sinuk;
        const double vy = xmy * cosuk - sinnok * sinuk;
        const double vz = sinik * cosuk;
        /*
         * position and velocity
         */
        const double x = rk * ux * kXKMPER;
        const double y = rk * uy * kXKMPER;
        const double z = rk * uz * kXKMPER;
        Vector position(x, y, z);
        const double xdot = (rdotk * ux + rfdotk * vx) * kXKMPER / 60.0;
        const double ydot = (rdotk * uy + rfdotk * vy) * kXKMPER / 60.0;
        const double zdot = (rdotk * uz + rfdotk * vz) * kXKMPER / 60.0;
        Vector velocity(xdot, ydot, zdot);

        return Eci(dt, position, velocity);
    }

    static inline double EvaluateCubicPolynomial(
        const double x,
        const double constant,
        const double linear,
        const double squared,
        const double cubed)
    {
        return constant + x * linear + x * x * squared + x * x * x * cubed;
    }

    void SGP4::Reset()
    {
        std::memset(&common_consts_, 0, sizeof(common_consts_));
        std::memset(&nearspace_consts_, 0, sizeof(nearspace_consts_));
        std::memset(&integrator_params_, 0, sizeof(integrator_params_));
    }

} // namespace libsgp4
