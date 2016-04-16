#pragma once

#include "ofMain.h"
#include "ofEvents.h"
#include "ofxOsc.h"
#include "fft.h"
#include "math.h"

#define HOST "localhost"
#define S_PORT 9000
#define R_PORT 9001
#define NUM_MSG_STRINGS 20

#define BUFFER_SIZE 256
#define NUM_WINDOWS 80

#define PIN_NUM 4


class ofApp : public ofBaseApp{
    
public:
    void setup();
    void update();
    void draw();
    
    void keyPressed(int key);
    void keyReleased(int key);
    void mousePressed(int x, int y, int button);
    void mouseReleased(int x, int y, int button);
    void setupArduino(const int & version);
    void updateArduino();
    void audioReceived 	(float * input, int bufferSize, int nChannels);
    
    
    ofImage img;
    ofTrueTypeFont font;
    int nowTime,targetTime;
    float bpm;
    bool bEditBpm;
    
    /*--------Arduino(LED)------*/
    ofArduino ard;
    bool bSetupArduino;
    int pin[4]={3,5,6,9};
    int ledMode=1;
    int lScene=0;
    bool bBeatAttack;
    int l_cnt=0;
    
    /*--------OSC---------*/
    ofxOscSender sender;
    ofxOscReceiver receiver;
    int current_mgs_string;
    string msg_strings[NUM_MSG_STRINGS];
    float timers[NUM_MSG_STRINGS];
    float beat,temp_beat;
    
    /*--------FFT----------*/
    float * left;
    float * right;
    int 	bufferCounter;
    fft		myfft;
    
    float magnitude[BUFFER_SIZE];
    float phase[BUFFER_SIZE];
    float power[BUFFER_SIZE];
    
    float freq[NUM_WINDOWS][BUFFER_SIZE/2];
    float freq_phase[NUM_WINDOWS][BUFFER_SIZE/2];
    int rect_color[3];
    int fftMode,preset_index;
    bool paramMode;
    
    
};
