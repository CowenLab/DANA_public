/////////////////////////////////////////////////////////////////
// This code syncs the voltammetry system with the behavioral systems
// It waits for a pulse from the voltammetry system. When it arrives, it
// then puts out a sequentially increasing binary code out of 8 output lines (so 0-255).
// It also sends this code (or a serialized version of it) to the
// behavioral acquisition system - or, if you are using a DAQ card to do the 
// syncing, send a line out from this (or split it) to the card and, in matlab, run VOLT_Acquire_Sync_Signals_From_DAQ
//
//
/////////////////////////////////////////////////////////////////
// Cowen 2014
/////////////////////////////////////////////////////////////////
#define N_SIGNAL_PINS 8 // each one goes to a unique input port on the Voltammetry system. Digital codes.
#define VOLT_STROBE_PIN 11 // Pin that waits for a TTL from the voltammetry system - SPLIT THIS AND ALSO SEND TO AMPLIPEX
#define BEHAVIOR_SIG_OUT_PIN 2  // Pin to the behaviorl system (e.g. Amplipex)
#define BEHAVIOR_SYNC_OUT_PIN 13 // Pin to the behaviorl system (e. g. Amplipex)
#define OUTPUT_BAUD 19200 //19200

int SignalPins[] = {3, 4, 5, 6, 7, 8, 9, 10}; // Pins on the arduino - will be specific to the type of arduino that you are using. You can choose - but don't use pins 0-4 and some others that are dedicated to other thing.
int Volt_Strobe_State = LOW;         // current state of the voltammetry sync input
int last_Volt_Strobe_State = LOW;    // previous state of the voltammetry sync input
int Scan_Counter = 0;   // counter for the number of button presses
int Master_Counter = 0;   // counter for the number of button presses
int i = 0;

void setup() {
  // run once:
  Serial.begin(OUTPUT_BAUD); // SENDS OUTPUT TO THE PC. WILL USE THE USB PORT CONNECTED TO THE ARDUINO.

  for (i = 0; i < N_SIGNAL_PINS; i++) {                                                                                                                                                                                       
    pinMode(SignalPins[i], OUTPUT);
    digitalWrite(SignalPins[i], LOW);
  }
  //  pinMode(VOLT_STROBE_PIN, INPUT_PULLUP);
  //  digitalWrite(VOLT_STROBE_PIN, HIGH);
  pinMode(VOLT_STROBE_PIN, INPUT);
  digitalWrite(VOLT_STROBE_PIN, LOW);
  Serial.print("START\n");
}

void loop() {
  // Runs continuously
  Volt_Strobe_State = digitalRead(VOLT_STROBE_PIN);
  //  delay(200); // for testing
  // compare the buttonState to its previous state
  if (Volt_Strobe_State != last_Volt_Strobe_State) {
    // if the state has changed, increment the counter
    if (Volt_Strobe_State == HIGH) {
      // if the current state is HIGH then the input
      // wend from off to on:
      Scan_Counter++;
      Master_Counter++;
      if (Scan_Counter == 256) {
        Scan_Counter = 0;
      }
      // Send out the digital output: NOTE - this will stay active until the next pulse is received. Why? Because the voltammetry system samples digital in relatively slowly so might as well just keep it on rather than risk it missing a signal.
      for (i = 0; i < N_SIGNAL_PINS; i++) {
         digitalWrite(SignalPins[i], bitRead(Scan_Counter, i));
         // digitalWrite(SignalPins[i], HIGH);
      }
      // Send out a serialzed output of the same digital val to the behavioral recording system.
      // Send out the digital output: NOTE - this will stay active until the next pulse is received. Why? Because the voltammetry system samples digital in relatively slowly so might as well just keep it on rather than risk it missing a signal.
      for (i = 0; i < N_SIGNAL_PINS; i++) {
        digitalWrite(BEHAVIOR_SIG_OUT_PIN, bitRead(Scan_Counter, i));
        digitalWrite(BEHAVIOR_SYNC_OUT_PIN, HIGH);
        delay(4);
        digitalWrite(BEHAVIOR_SIG_OUT_PIN, LOW);
        digitalWrite(BEHAVIOR_SYNC_OUT_PIN, LOW);
        delay(4);
      }

      Serial.println("on");
      Serial.print("number of scans (resets at 255):  ");
      Serial.print(Scan_Counter);
      Serial.print(" ");
      Serial.println(Master_Counter);
      //delay(50); // delay a bit to avoid any bounce on the input line. Probably unnecessary. May need a pulldown resistor on the input line (e.g. 10K to ground)
    }
    else {
      // if the current state is LOW then the button
      // wend from on to off:
      Serial.println("off");
    }
  } // checks if last state was the same or different.

  // save the current state as the last state,
  //for next time through the loop
  last_Volt_Strobe_State = Volt_Strobe_State;
}
