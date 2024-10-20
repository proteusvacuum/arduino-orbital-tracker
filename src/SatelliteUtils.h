#ifndef SAT_UTILS_H
#define SAT_UTILS_H

#include "libsgp4/SGP4.h"
#include "Config.h"
#include "WifiUtils.h"

struct Direction
{
  double azimuth;
  double elevation;

  Direction(double az, double el) : azimuth{az}, elevation{el} {};
};

class TleLines
{
public:
  std::string line1;
  std::string line2;
  std::string cat_number;
  std::chrono::time_point<std::chrono::system_clock> last_fetch;

  WifiHandler handler;
  TleLines(std::string cat_number, WifiHandler &handler) : cat_number{cat_number}, handler{handler} {};

  boolean shouldFetch()
  {
    auto now = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::minutes>(now - last_fetch);
#ifdef DEBUG
    std::cout << "elapsed: " << elapsed.count() << std::endl;
#endif
    return elapsed.count() >= 30;
  }
  /// @brief Parses the response from the TLE request
  /// @param handler
  void getTLELines()
  {
#ifdef DEBUG
    std::cout << "Fetching TLE" << std::endl;
#endif
    std::string raw = handler.getRawTLEData(cat_number);
    size_t pos = 0;
    constexpr char delimiter[] = "\r\n";
    std::string line;

    while ((pos = raw.find(delimiter)) != std::string::npos)
    {
      line = raw.substr(0, pos);
      if (line[0] == '1')
      {
        line1 = line;
      }
      else if (line[0] == '2')
      {
        line2 = line;
      }
      raw.erase(0, pos + strlen(delimiter));
    }
    last_fetch = std::chrono::system_clock::now();
  };
};

class Orbits
{
public:
  Orbits(double observer_lat = 46.003641, double observer_lon = -74.182330,
         double observer_altitude = 0.2808)
      : observer_lat{observer_lat}, observer_lon{observer_lon},
        observer_altitude{observer_altitude} {};

  double observer_lat;      // Example: 45 degrees latitude
  double observer_lon;      // Example: -75 degrees longitude
  double observer_altitude; // km

  Direction getDirection(std::string &tle_line1, std::string &tle_line2)
  {
    libsgp4::DateTime now = libsgp4::DateTime().Now();
    libsgp4::Vector iss_vector = getSatelliteVector(tle_line1, tle_line2, now);
    Direction direction =
        calculateAzEl(observer_lat, observer_lon, observer_altitude,
                      iss_vector.x, iss_vector.y, iss_vector.z, now);
    return direction;
  }

private:
  static constexpr double pi = 3.14159265358979323846;

  // WGS84 Ellipsoid constants
  static constexpr double a = 6378.137; // Earth's equatorial radius in km
  static constexpr double f = 1.0 / 298.257223563;
  static constexpr double e2 = f * (2 - f); // Square of eccentricity

  libsgp4::Vector getSatelliteVector(std::string &tle_line1,
                                     std::string &tle_line2,
                                     libsgp4::DateTime &now)
  {
    libsgp4::Tle tle_file = libsgp4::Tle(tle_line1, tle_line2);
    libsgp4::SGP4 sgp4(tle_file);
    libsgp4::Eci new_pos = sgp4.FindPosition(now);
    return new_pos.Position();
  }

  Direction calculateAzEl(double lat, double lon, double alt, double satX,
                          double satY, double satZ, libsgp4::DateTime now)
  {
    // Convert lat/lon to radians
    const double phi{lat * pi / 180.0};
    const double lambda{lon * pi / 180.0};
    const double h{alt};

    const double sin_phi{sin(phi)};
    const double N{a / sqrt(1 - e2 * sin_phi * sin_phi)};

    // Compute ECEF coordinates of the observer
    const double Xo_ecef{(N + h) * cos(phi) * cos(lambda)};
    const double Yo_ecef{(N + h) * cos(phi) * sin(lambda)};
    const double Zo_ecef{(N * (1 - e2) + h) * sin(phi)};

    const double theta{now.ToGreenwichSiderealTime()};

    // Rotate ECEF coordinates to ECI coordinates
    const double observerX{Xo_ecef * cos(theta) - Yo_ecef * sin(theta)};
    const double observerY{Xo_ecef * sin(theta) + Yo_ecef * cos(theta)};
    const double observerZ{Zo_ecef};

    // Step 3: Compute the vector from observer to satellite in ECI
    const double dX{satX - observerX};
    const double dY{satY - observerY};
    const double dZ{satZ - observerZ};

    // Step 4: Compute local unit vectors (East, North, Up)
    // Up vector (U)
    const double norm_O{sqrt(observerX * observerX + observerY * observerY + observerZ * observerZ)};
    const double Ux{observerX / norm_O};
    const double Uy{observerY / norm_O};
    const double Uz{observerZ / norm_O};

    // East vector (E)
    const double norm_E{sqrt((-Uy) * (-Uy) + (Ux) * (Ux) + 0.0)};
    const double Ex{-Uy / norm_E};
    const double Ey{Ux / norm_E};
    const double Ez{0.0};

    // North vector (N) = U x E
    const double Nx{Uy * Ez - Uz * Ey};
    const double Ny{Uz * Ex - Ux * Ez};
    const double Nz{Ux * Ey - Uy * Ex};

    // Step 5: Project the satellite vector onto the ENU coordinates
    const double norm_d{sqrt(dX * dX + dY * dY + dZ * dZ)};
    const double dX_unit{dX / norm_d};
    const double dY_unit{dY / norm_d};
    const double dZ_unit{dZ / norm_d};

    // Compute ENU components
    const double E_comp{Ex * dX_unit + Ey * dY_unit + Ez * dZ_unit};
    const double N_comp{Nx * dX_unit + Ny * dY_unit + Nz * dZ_unit};
    const double U_comp{Ux * dX_unit + Uy * dY_unit + Uz * dZ_unit};

    // Step 6: Calculate Azimuth and Elevation
    double azimuth{atan2(E_comp, N_comp) * 180.0 / pi};
    if (azimuth < 0.0)
    {
      azimuth += 360.0;
    }
    const double elevation{asin(U_comp) * 180.0 / pi};
#ifdef DEBUG
    std::cout << "Azimuth: " << azimuth << " degrees" << std::endl;
    std::cout << "Elevation: " << elevation << " degrees" << std::endl;
#endif
    return Direction{azimuth, elevation};
  }
};

#endif // SAT_UTILS_H
