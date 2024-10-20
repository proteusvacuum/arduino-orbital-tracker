#include "src/satellite_utils.h"
#include "Arduino.h"
#include "AccelStepper.h"
#include "RTC.h"
#include "Servo.h"
#include <NTPClient.h>

#include "src/wifi_utils.h"

// Motor pin definitions:
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

Orbits orbitFinder;
// http://celestrak.org/NORAD/elements/gp.php?CATNR=25544&FORMAT=TLE
constexpr const char ISS_cat_num[] = "25544";
TleLines tle_lines{ISS_cat_num, wifiHandler};

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
  std::cout << "direction.azimuth: " << direction.azimuth << "direction.elevation: " << direction.elevation << std::endl;
  stepper.runToNewPosition(getStepsFromAngle(direction.azimuth));
  // Serial.println(stepper.currentPosition());
  delay(1000);
}
