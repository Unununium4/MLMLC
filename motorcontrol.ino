#include <digitalWriteFast.h>
#include "Stepper.h"

#define ENABLE 6
#define DIRA 4
#define DIRB 5
#define channelA 2
#define channelB 3

volatile int pos = 0;
int gotoPos = -200;
bool posFound = false;

void isrA ()  {
  if(!digitalReadFast(channelB)){
    pos++;
  }
  else{
    pos--;
  }
  testPos();
}

void setup() {
  pinMode(ENABLE,OUTPUT);
  pinMode(DIRA,OUTPUT);
  pinMode(DIRB,OUTPUT);
  pinMode(channelA,INPUT_PULLUP);
  pinMode(channelB,INPUT);

  attachInterrupt (digitalPinToInterrupt(channelA),isrA,RISING);
  Serial.begin(9600);

  testPos();
}

void testPos(){
  if(pos>gotoPos){
    if(abs(pos-gotoPos)<100){
        analogWrite(ENABLE,40);
    }
    else{
        analogWrite(ENABLE,200);
    }
    digitalWrite(DIRA,LOW);
    digitalWrite(DIRB,HIGH);
  }
  else if(pos<gotoPos){
    if(abs(pos-gotoPos)<100){
        analogWrite(ENABLE,40);
    }
    else{
        analogWrite(ENABLE,200);
    }
    digitalWrite(DIRA,HIGH);
    digitalWrite(DIRB,LOW);
  }
  else{
    posFound=true;
    analogWrite(ENABLE,0);
    digitalWrite(DIRA,LOW);
    digitalWrite(DIRB,LOW);
  }
}

void loop() {
  if(posFound){
    Serial.println(pos);
    posFound=false;
    gotoPos=-200;
    testPos();
  }
  delay(100);
  if(posFound){
    Serial.println(pos);
    posFound=false;
    gotoPos=200;
    testPos();    
  }
  delay(100);
}
   
