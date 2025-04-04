#include "src/Config.h"
#include "src/SatelliteUtils.h"
#include "Arduino.h" // Arduino.h must come after satellite_utils.h, as Arduino.h ovewrites the `abs` macro that the sgp4 library depends on
#include "AccelStepper.h"
#include "RTC.h"
#include "Servo.h"
#include <NTPClient.h>
#include "src/WifiUtils.h"

// This is the object to track. 25544 is the ISS. Find other interesting objects at http://celestrak.org/NORAD/elements/
constexpr const char CATALOG_NUMBER[] = "25544";

constexpr int motorPin1 = 8;  // IN1 on the ULN2003 driver
constexpr int motorPin2 = 9;  // IN2 on the ULN2003 driver
constexpr int motorPin3 = 10; // IN3 on the ULN2003 driver
constexpr int motorPin4 = 11; // IN4 on the ULN2003 driver

constexpr int servoPin = 12;

// 4 = Full Step mode (2048 steps)
// 8 = Half Step mode (4096 steps)
constexpr int const MotorInterfaceType = 4;
constexpr int const stepsPerRotation = 2048;

AccelStepper stepper = AccelStepper(MotorInterfaceType, motorPin1, motorPin3, motorPin2, motorPin4);
Servo servo;

WifiHandler wifiHandler;
NTPClient timeClient(wifiHandler.Udp);

Orbits orbitFinder{LATITUDE, LONGITUDE, ALTITUDE};

TleLines tle_lines{CATALOG_NUMBER, wifiHandler};

int getStepsFromAngle(int angle)
{
  return (angle / 360.0) * stepsPerRotation;
}

void setup()
{
  stepper.setMaxSpeed(400);
  stepper.setAcceleration(100);

  servo.attach(servoPin);

  Serial.begin(9600);
  while (!Serial)
    ;

  wifiHandler.connectToWiFi();

#ifdef DEBUG
  wifiHandler.printWifiStatus();
#endif

  RTC.begin();
  timeClient.begin();
  timeClient.update();
  RTCTime timeToSet = RTCTime(timeClient.getEpochTime());
  RTC.setTime(timeToSet);

  tle_lines.getTLELines();
}

void loop()
{
  if (tle_lines.shouldFetch())
  {
    tle_lines.getTLELines();
  }
  Direction direction = orbitFinder.getDirection(tle_lines.line1, tle_lines.line2);
  stepper.runToNewPosition(getStepsFromAngle(direction.azimuth));
  servo.write(direction.elevation + 90);
  delay(1000);
}
