#include <math.h>
#include <ArduinoJson.h>
/*
Connections of Drive and Arduino
Serial Port 0 is not used to connect to drive because its connected to USB-Serial and used to show information on console.

For Arduino Uno Software serial needs to be used as there is only one hardware serial port and its connected to USB-Serial. 
   Drive to Arduino UNO/Nano connections
   GND         -      GND
   RXD         -      D3
   TXD         -      D2

For arduino mega and other arduinos with multiple hardware serial port, any port other than 0 can be selected to connect the drive.

   Drive to Arduino Mega2560 connections
   GND         -      GND
   RXD         -      Tx1/Tx2/Tx3
   TXD         -      Rx1/Rx2/Rx3
   
*	This mode can be used when multiple motors are to be used to run at exactly the same RPM and same torque even though the voltage supply might be different.
*	Also in this mode the direction of the motor can be controlled digitally via modbus ASCII commands to run the dc servo motor in both directions

* For more information see : https://robokits.co.in/motor-drives-drivers/encoder-dc-servo/rhino-dc-servo-driver-50w-compatible-with-modbus-uart-ascii-for-encoder-dc-servo-motor


*/

#include<RMCS2303drive.h>

RMCS2303 rmcs;                    //object for class RMCS2303

//SoftwareSerial myserial(2,3);     //Software Serial port For Arduino Uno. Comment out if using Mega.

//parameter Settings "Refer datasheet for details"
byte slave_id=7;
int INP_CONTROL_MODE=257;           
int PP_gain=32;
int PI_gain=16;
int VF_gain=32;
int LPR=334;
int acceleration=5000;
int speed=8000;
float i = 0;
float ADC_conv = 5000/1024;   // 5V/max byte value
float offset_adc = 506; //ADC value at no pwm
float offset_mvolt = offset_adc*ADC_conv; //zero error milliVolts
float sensitivity = 100;  //mV/A

long int Current_position;
int motor_Speed;
long int motor_Position;
int PWM;
int theta1;
float current;
int sign;
long int motor_Position_initial = 0;
float ACS_zero = 0;
float position_multiplie_1000_count_deg = 2.7428;
float degree_callibration_constant = 52;

float I_cur = 0;
float I_prev = 0;

int pwmPin = 3;
int d1 = 4;

float angularSpeed1 = 0.0; // Average angular speed

//Serial communication using JSON
JsonDocument Feedback;


void setup()
{
   rmcs.Serial_selection(0);       //Serial port selection:0-Hardware serial,1-Software serial
   rmcs.Serial0(9600);             //Set baudrate for usb serial to monitor data on serial monitor

   rmcs.begin(&Serial3,9600);    //Uncomment if using hardware serial port for mega2560:Serial1,Serial2,Serial3 and set baudrate. Comment this line if Software serial port is in use
   //rmcs.WRITE_PARAMETER(slave_id,INP_CONTROL_MODE,PP_gain,PI_gain,VF_gain,LPR,acceleration,speed);    //Uncomment to write parameters to drive. Comment to ignore.
   //rmcs.READ_PARAMETER(slave_id);
   
   pinMode(pwmPin, OUTPUT);
   pinMode(d1, OUTPUT);

   Serial.begin(9600);

   float average = 0;
   for (int i = 0; i < 1000; i++) {
   average = average + (analogRead(A0));
   }
   average = average / 1000;
   ACS_zero = ((average / 1024) * 5000);

   //Serial.println("ACS zero is");
   //Serial.println(ACS_zero);

   motor_Position_initial=rmcs.Position_Feedback(slave_id);
   //Serial.println("Initial Motor Position is");
   //Serial.println(motor_Position_initial);

}

float current_time = 0;
float previous_time = 0;
float sampling_time = 5000; //200Hz or 5ms
float currentValue = 0;
float previousValue = 0;
int signum = 1;

//Processing of PWM data from Python
int receivedPWM = 40; // Default value

void receivePWMJson()
{
  static String inputBuffer;
  while (Serial.available()) {
    char c = Serial.read();
    if (c == '\n') {
      // Try to parse JSON
      StaticJsonDocument<64> doc;
      DeserializationError error = deserializeJson(doc, inputBuffer);
      if (!error && doc.containsKey("pwm")) {
        int val = doc["pwm"];
        receivedPWM = constrain(val, 0, 40);
      }
      if (!error && doc.containsKey("sign")) {
        signum = doc["sign"];
        receivedPWM = receivedPWM*signum;
      }
      inputBuffer = "";
    } else {
      inputBuffer += c;
    }
  }
}

void loop(void)
{
  current_time = micros();
  static String outputBuffer;
  if(current_time - previous_time > sampling_time)
  {
      receivePWMJson();
      //receivedPWM = 25;  testing purposes for motor position calibration at low PWM value
      
      if (signum <0){
        analogWrite(pwmPin, abs(receivedPWM));
        analogWrite(d1, 0);
      }

      else{
        analogWrite(d1, receivedPWM);
        analogWrite(pwmPin, 0);
      }
 

      motor_Speed=rmcs.Speed_Feedback(slave_id)*0.05 - 0.9127; 
      motor_Position=rmcs.Position_Feedback(slave_id);
      
      theta1 = analogRead(A1);

      //speed of theta1 calculated

      currentValue = analogRead(A1);
      float difference = currentValue - previousValue;
      angularSpeed1 = difference / sampling_time; 
      previousValue = currentValue;

      float average = 0;
      float current_pwm = 0;
      for(int i = 0; i < 1000; i++) {
        average = average + (analogRead(A0));
      }
      average = average/1000;
      current_pwm = (average/1024)*5000;
      current = (current_pwm - ACS_zero)/100;
      I_cur = current;
      di/dt = (I_cur-I_prev)/(current_time - previous_time);
      I_prev = I_cur;
      Serial.println(di/dt);

      //Serial.println(current_pwm);
      
      Feedback["motor_Speed"] = motor_Speed;
      Feedback["motor_Position"] = motor_Position*position_multiplie_1000_count_deg*degree_callibration_constant/1000;
      Feedback["theta1"] = theta1;
      Feedback["ang_speed1"] = angularSpeed1;
      Feedback["current"] = current;
      Feedback["time"] = current_time;

      serializeJson(Feedback, Serial);

      Serial.println("");

      //String message = String(theta1) + "<" + String(angularSpeed1) + "<" + String(average);
      //Serial.println(message);

      previous_time = current_time;

  }
  
}