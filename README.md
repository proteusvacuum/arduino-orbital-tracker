# Arduino Orbital Tracker

Points at a satellite in the sky. Great for tracking the ISS. For use with an Arduino UNO R4 WiFi.

To start, copy the secrets file, and update it with your location and WiFi information:

```shell
cp ./arduino_secrets.h.example ./arduino_secrets.h
```

Currently it is set to track the ISS, but this should work for any satellite with a catalog number. Change the value of `CATALOG_NUMBER` to the object you'd like to track.

For debugging, uncomment the line in `./src/Config.h`. This will print lots of interesting information to the Serial console.

## Note
This only compiles with Arduino IDE v2.3.3 which can be downloaded [here](https://github.com/arduino/arduino-ide/releases/tag/2.3.3). Newer versions create larger binaries that no longer fit on the Arduino R4.
