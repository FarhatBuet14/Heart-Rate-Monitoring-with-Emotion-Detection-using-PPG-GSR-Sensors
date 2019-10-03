#include <SoftwareSerial.h>

#include <Arduino.h>
#include <MPU6050_tockn.h>
#include <math.h>
#include <Wire.h>

#include "MAX30100.h"

MAX30100* pulseOxymeter;
MPU6050 mpu6050(Wire);


SoftwareSerial BT(1, 0);

void setup() {
  Wire.begin();
  //Serial.begin(115200);
  //Serial.println("Date & Time, ACCx, ACCy, ACCz");
//  Serial.println("Pulse oxymeter test!");

  //pulseOxymeter = new MAX30100( DEFAULT_OPERATING_MODE, DEFAULT_SAMPLING_RATE, DEFAULT_LED_PULSE_WIDTH, DEFAULT_IR_LED_CURRENT, true, true );
  pulseOxymeter = new MAX30100();
  pinMode(2, OUTPUT);
  mpu6050.begin();
  mpu6050.calcGyroOffsets(true);
  pinMode(A0, INPUT);
  //pulseOxymeter->printRegisters();
  BT.begin(9600);
  // Send test message to other device
  BT.println("Hello from Arduino");
  int startMillis = millis();
}

char a; // stores incoming character from other device

void loop() {

    mpu6050.update();
    //return;
    //You have to call update with frequency at least 37Hz. But the closer you call it to 100Hz the better, the filter will work.
    pulseoxymeter_t result = pulseOxymeter->update();


    BT.print(result.dcFilteredIR);
  //  Serial.print("|RED|255,0,0|");
    BT.print(",");
    BT.print(result.dcFilteredRed);
 
  //  Serial.print("accX : ");
    BT.print(",");
    BT.print(mpu6050.getAccX());
  //  Serial.print("\taccY : ");
    BT.print(",");
    BT.print(mpu6050.getAccY());
  //  Serial.print("\taccZ : ");
    BT.print(",");
    BT.print(mpu6050.getAccZ());
    BT.print(",");
    
    BT.print(analogRead(A0));
    BT.print(",");
    BT.println(millis());
    
    delay(10);
  
    //Basic way of determening execution of the loop via oscoliscope
    digitalWrite( 2, !digitalRead(2) );

  
}
