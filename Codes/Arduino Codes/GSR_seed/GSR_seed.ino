#include <Arduino.h>

#include <math.h>

const int GSR=A3;
int sensorValue=0;
float gsr_average=0;
float resistance=0;
float conductance=0;

void setup(){
  Serial.begin(9600);
  int startMillis = millis();
}

void loop(){
  long sum=0;
  for(int i=0;i<10;i++)           //Average the 10 measurements to remove the glitch
      {
      sensorValue=analogRead(GSR);
      sum += sensorValue;
      delay(5);
      }
   gsr_average = sum/10;
   resistance=((1024+2*gsr_average)*10000)/(abs(512-gsr_average));
   conductance=(1/resistance)*1000000;//micro Siemens;
   
   Serial.print(gsr_average);
   Serial.print(",");
   Serial.print(resistance);
   Serial.print(",");
   Serial.print(conductance);
   Serial.print(",");
   Serial.println(millis());
   
}

//Human Resistance = ((1024+2*Serial_Port_Reading)*10000)/(512-Serial_Port_Reading), unit is ohm, Serial_Port_Reading is the value display on Serial Port(between 0~1023)
