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

#pragma once

#include "Tle.h"
#include "OrbitalElements.h"
#include "Eci.h"
#include "SatelliteException.h"
#include "DecayedException.h"

namespace libsgp4
{

    /**
     * @mainpage
     *
     * This documents the SGP4 tracking library.
     */

    /**
     * @brief The simplified perturbations model 4 propagater.
     */
    class SGP4
    {
    public:
        explicit SGP4(const Tle &tle)
            : elements_(tle)
        {
            Initialise();
        }

        void SetTle(const Tle &tle);
        Eci FindPosition(double tsince) const;
        Eci FindPosition(const DateTime &date) const;

    private:
        struct CommonConstants
        {
            double cosio;
            double sinio;
            double eta;
            double t2cof;
            double x1mth2;
            double x3thm1;
            double x7thm1;
            double aycof;
            double xlcof;
            double xnodcf;
            double c1;
            double c4;
            double omgdot; // secular rate of omega (radians/sec)
            double xnodot; // secular rate of xnode (radians/sec)
            double xmdot;  // secular rate of xmo   (radians/sec)
        };

        struct NearSpaceConstants
        {
            double c5;
            double omgcof;
            double xmcof;
            double delmo;
            double sinmo;
            double d2;
            double d3;
            double d4;
            double t3cof;
            double t4cof;
            double t5cof;
        };

        struct IntegratorParams
        {
            /*
             * integrator values
             */
            double xli;
            double xni;
            double atime;
        };

        void Initialise();
        static void RecomputeConstants(const double xinc,
                                       double &sinio,
                                       double &cosio,
                                       double &x3thm1,
                                       double &x1mth2,
                                       double &x7thm1,
                                       double &xlcof,
                                       double &aycof);
        Eci FindPositionSGP4(double tsince) const;
        static Eci CalculateFinalPositionVelocity(
            const DateTime &date,
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
            const double sinio);

        /**
         * Reset
         */
        void Reset();

        /*
         * the constants used
         */
        struct CommonConstants common_consts_;
        struct NearSpaceConstants nearspace_consts_;
        mutable struct IntegratorParams integrator_params_;

        /*
         * the orbit data
         */
        OrbitalElements elements_;
    };

} // namespace libsgp4
