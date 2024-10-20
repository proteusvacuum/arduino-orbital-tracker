#ifndef WIFI_UTILS_H
#define WIFI_UTILS_H

#include <WiFiS3.h>
#include "Config.h"
#include "../arduino_secrets.h"

class WifiHandler
{
public:
  WiFiUDP Udp;
  void connectToWiFi()
  {
    while (wifiStatus != WL_CONNECTED)
    {
      wifiStatus = WiFi.begin(ssid_.c_str(), pass_.c_str());
      delay(10000);
    }
  };

  std::string getRawTLEData(std::string cat_number)
  {
    std::string data;
    bool server = connectToServer_(cat_number);
    while (!readData_(data))
      ;
#ifdef DEBUG
    std::cout << "Received data: ==== \n"
              << data << std::endl;
#endif
    return data;
  };

  void printWifiStatus()
  {
    Serial.print("SSID: ");
    Serial.println(WiFi.SSID());

    IPAddress ip = WiFi.localIP();
    Serial.print("IP Address: ");
    Serial.println(ip);

    long rssi = WiFi.RSSI();
    Serial.print("signal strength (RSSI):");
    Serial.print(rssi);
    Serial.println(" dBm");
  }

private:
  int wifiStatus = WL_IDLE_STATUS;
  WiFiClient client;
  std::string ssid_ = SECRET_SSID;
  std::string pass_ = SECRET_PASS;

  bool connectToServer_(std::string &cat_number)
  {
    // 104.168.149.178 is celestrak.org
    IPAddress celestrak(104, 168, 149, 178);
    if (client.connect(celestrak, 80))
    {
      client.println(("GET /NORAD/elements/gp.php?CATNR=" + cat_number + "&FORMAT=TLE HTTP/1.1").c_str());
      client.println("Host: celestrak.org");
      client.println("Connection: close");
      client.println();
      return true;
    }
    return false;
  }

  bool readData_(std::string &data)
  {
    if (!client.connected())
    {
      client.stop();
      return true;
    }
    while (client.available())
    {
      char c = client.read();
      data += c;
    }
    return false;
  }
};

#endif // WIFI_UTILS_H
