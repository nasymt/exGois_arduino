#include "ofApp.h"

//--------------------------------------------------------------
void ofApp::setup(){
    //画面設定
    ofSetVerticalSync(true);
    ofSetFrameRate(60);
    ofBackground(255,0,130);
    
    font.loadFont("Avenir.ttc", 14);
    font.setLineHeight(14);
    /*--------------arduino-------------*/
    ard.connect("/dev/cu.usbmodem1411", 57600);
    ofAddListener(ard.EInitialized, this, &ofApp::setupArduino);
    bSetupArduino	= false;
    /*-------------OSC--------------*/
    sender.setup(HOST,S_PORT);
    receiver.setup(R_PORT);
    /*-------------FFT--------------*/
    srand((unsigned int)time((time_t *)NULL));
    ofSoundStreamSetup(0,2,this, 44100,BUFFER_SIZE, 4);
    left = new float[BUFFER_SIZE];
    right = new float[BUFFER_SIZE];
    
    for (int i = 0; i < NUM_WINDOWS; i++){
        for (int j = 0; j < BUFFER_SIZE/2; j++){
            freq[i][j] = 0;
        }
    }
    myfft.setup();
    //fftMode=0;
}

//--------------------------------------------------------------
void ofApp::update(){
    ofBackground(150,150,150);
    updateArduino();
    
    /*-----------OSC-------------*/
    while(receiver.hasWaitingMessages()){
        ofxOscMessage m;
        receiver.getNextMessage(m);
        
        if(m.getAddress() == "/bpm"){
            bpm = m.getArgAsInt32(0);
            cout<<"BPM change! >>"<<bpm<<endl;
        }
        else{
            string msg_string;
            msg_string = m.getAddress();
            msg_string += ": ";
            for(int i = 0; i < m.getNumArgs(); i++){
                msg_string += m.getArgTypeName(i);
                msg_string += ":";
                
                if(m.getArgType(i) == OFXOSC_TYPE_INT32){
                    msg_string += ofToString(m.getArgAsInt32(i));
                }
                else if(m.getArgType(i) == OFXOSC_TYPE_FLOAT){
                    msg_string += ofToString(m.getArgAsFloat(i));
                }
                else if(m.getArgType(i) == OFXOSC_TYPE_STRING){
                    msg_string += m.getArgAsString(i);
                }
                else{
                    msg_string += "unknown";
                }
            }
        }
    }
    
    if(beat>0){
        nowTime = ofGetElapsedTimeMillis();
        if(nowTime>=targetTime){
            if(beat<4)beat++;
            else beat=1;
            float nextBeat = 1000/(bpm/60);
            targetTime=nowTime+nextBeat;
            cout<<beat<<endl;
            bBeatAttack=true;
        }
    }
    if(myfft.val[0]>0){
        ofxOscMessage m;
        m.setAddress("/vol");
        m.addIntArg(myfft.val[0]);
        m.addIntArg(myfft.val[1]);
        m.addIntArg(myfft.val[2]);
        m.addIntArg(myfft.val[3]);
        sender.sendMessage(m,false);
    }
}

//--------------------------------------------------------------
void ofApp::setupArduino(const int & version) {
    ofRemoveListener(ard.EInitialized, this, &ofApp::setupArduino);
    for (int i = 0; i < 13; i++){
        ard.sendDigitalPinMode(i, ARD_OUTPUT);
    }
    bSetupArduino = true;
}
//--------------------------------------------------------------
void ofApp::updateArduino(){
    ard.update();
    /*-----------LED 点灯パターン----------*/
    if(ledMode==1){
        for(int i=0;i<4;i++){
            if(beat==i+1){
                ard.sendDigital(pin[i], ARD_HIGH);
                if(i==0)ard.sendDigital(pin[3], ARD_LOW);
                else ard.sendDigital(pin[i-1], ARD_LOW);
            }
        }
    }else if(ledMode==2){
        for(int i=0;i<4;i++)ard.sendPwm(pin[i],(int)(128+128 * sin(ofGetElapsedTimef()*2)));
    }else if(ledMode==3){
        if(lScene==0){
            if(bBeatAttack){
                if(beat==1&&l_cnt==0)l_cnt=1;
                else if(l_cnt>0)l_cnt++;
                for(int i=0;i<4;i++)if(l_cnt==i+1)ard.sendDigital(pin[i], ARD_HIGH);
                for(int i=0;i<4;i++)if(l_cnt==i+5)ard.sendDigital(pin[i], ARD_LOW);
                if(l_cnt==9){
                    ard.sendDigital(pin[0], ARD_HIGH);
                    ard.sendDigital(pin[3], ARD_HIGH);
                }else if(l_cnt==10){
                    ard.sendDigital(pin[1], ARD_HIGH);
                    ard.sendDigital(pin[2], ARD_HIGH);
                }else if(l_cnt==11){
                    ard.sendDigital(pin[0], ARD_LOW);
                    ard.sendDigital(pin[3], ARD_LOW);
                }else if(l_cnt==12){
                    ard.sendDigital(pin[1], ARD_LOW);
                    ard.sendDigital(pin[2], ARD_LOW);
                }else if(l_cnt==13){
                    ard.sendDigital(pin[0], ARD_HIGH);
                    ard.sendDigital(pin[2], ARD_HIGH);
                }else if(l_cnt==14){
                    ard.sendDigital(pin[1], ARD_HIGH);
                    ard.sendDigital(pin[3], ARD_HIGH);
                }else if(l_cnt==15){
                    ard.sendDigital(pin[0], ARD_LOW);
                    ard.sendDigital(pin[2], ARD_LOW);
                }else if(l_cnt==16){
                    ard.sendDigital(pin[1], ARD_LOW);
                    ard.sendDigital(pin[3], ARD_LOW);
                    l_cnt=0;
                }
                bBeatAttack=false;
            }
            
        }
    }else if(ledMode==4){
        for(int i=0;i<PIN_NUM;i++){
            ard.sendPwm(pin[i], (int)rect_color[i]);
        }
    }
}
//--------------------------------------------------------------
void ofApp::draw(){
    /*-------------FFT---------------*/
    static int index=0;
    float avg_power = 0.0f;
    if(index < 80)
        index += 1;
    else
        index = 0;
    
    myfft.powerSpectrum(0,(int)BUFFER_SIZE/2, left,BUFFER_SIZE,&magnitude[0],&phase[0],&power[0],&avg_power);
    
    for(int j=1; j < BUFFER_SIZE/2; j++) {
        freq[index][j] = magnitude[j];
    }
    
    /* draw the FFT */
    for (int i = 1; i < (int)(BUFFER_SIZE/2); i++){
        if(myfft.band_bottom[0]<=i&&myfft.band_top[0]>i)ofSetColor(255, 0, 0);
        else if(myfft.band_bottom[1]<=i&&myfft.band_top[1]>i)ofSetColor(0, 255, 0);
        else if(myfft.band_bottom[2]<=i&&myfft.band_top[2]>i)ofSetColor(0, 255, 255);
        else if(myfft.band_bottom[3]<=i&&myfft.band_top[3]>i)ofSetColor(0, 0, 255);
        else ofSetColor(50,50,50);
        ofDrawLine(100+(i*8),400,100+(i*8),400-magnitude[i]*10.0f);
    }
    ofSetColor(255);
    for(int i=0;i<4;i++){
        myfft.update(magnitude, i);
        ofDrawCircle(150+i*250, 100, myfft.val[i]*50);
        string string_index[] = {"low:","mid:","mid2:","high:"};
        ofDrawBitmapString(ofToString(string_index[i]), 50, 480+i*30);
        ofDrawBitmapString(ofToString(myfft.map_min[i]), 100, 480+i*30);
        ofDrawBitmapString(ofToString(myfft.map_max[i]), 150, 480+i*30);
        ofDrawBitmapString(ofToString(myfft.map_newMin[i]), 200, 480+i*30);
        ofDrawBitmapString(ofToString(myfft.map_newMax[i]), 250, 480+i*30);
        ofDrawBitmapString(ofToString(myfft.temp_val),300,480+i*30);
    }
    ofSetColor(255);
    ofDrawBitmapString(ofToString(paramMode),100,450);
    ofDrawBitmapString("bSmooth "+ofToString(myfft.bSmooth)+":"+ofToString(myfft.smoothRate), 100, 620);
    ofDrawBitmapString("BPM:"+ofToString(bpm), 600, 670);
    if(myfft.bSelectPreset)ofDrawBitmapString("===SELECT PRESET(Press key 1-2, 0 is reset)=== ", 100, 650);
    
    if(myfft.bSmooth){
        ofSetColor(200, 0, 103);
        ofDrawRectangle(100,650,100,50);
        ofSetColor(0, 0, 0);
        font.drawString("smooth", 115, 678);
    }
    else{
        ofSetColor(0, 0, 0);
        ofDrawRectangle(100,650,100,50);
        ofSetColor(200, 0, 103);
        font.drawString("no smooth", 105, 678);
    }
    
    if(!myfft.bReset){
        ofSetColor(0, 0, 0);
        ofDrawRectangle(230, 650, 100, 50);
        ofSetColor(200, 0, 103);
        font.drawString("Reset", 245, 678);
    }else {
        ofSetColor(200, 0, 103);
        ofDrawRectangle(230,650,100,50);
        ofSetColor(0, 0, 0);
        font.drawString("Reset", 245, 678);
        myfft.bReset=false;
        for(int i=0;i<3;i++)myfft.vol_max[i]=0;
    }
    
    if(myfft.bAutoMaxGet){
        ofSetColor(200, 0, 103);
        ofDrawRectangle(360, 650, 100, 50);
        ofSetColor(0, 0, 0);
        font.drawString("AutoGetMax", 370, 678);
    }else{
        ofSetColor(0, 0, 0);
        ofDrawRectangle(360, 650, 100, 50);
        ofSetColor(200,0,103);
        font.drawString("ManualGet", 370, 678);
    }
    for(int i=0;i<4;i++){
        rect_color[i] = ofMap(myfft.val[i], 0, 2, 0, 255);
        ofSetColor(rect_color[i], rect_color[i], rect_color[i]);
        ofDrawRectangle(400+i*150,450,120,120);
    }
}

//--------------------------------------------------------------
void ofApp::keyPressed(int key){
    /*---------------------------------
     モード選択:
     ◆リターン：モード変更
     0：低中高音の帯域指定
     1：数値のマッピング用
     ◆シフト：平滑化オンオフ
     　キーの上下で平滑化係数の調整
     ◆右コマンド：プリセットの選択
     
     ---------------------------------*/
    if(key==OF_KEY_RETURN){
        if(!paramMode)paramMode=true;
        else if(paramMode)paramMode=false;
        cout<<"mode:"<<paramMode<<endl;
    }
    else if(key==OF_KEY_SHIFT){
        if(myfft.bSmooth)myfft.bSmooth=false;
        else myfft.bSmooth=true;
    }else if(key==OF_KEY_RIGHT_COMMAND){
        if(!bEditBpm)bEditBpm=true;
        else if(bEditBpm)bEditBpm=false;
    }
    if(myfft.bSmooth){
        if(key==OF_KEY_UP)myfft.smoothRate+=0.05;
        else if(key==OF_KEY_DOWN)myfft.smoothRate-=0.05;
    }
    if(bEditBpm){
        if(key==OF_KEY_UP)bpm++;
        else if(key==OF_KEY_DOWN)bpm--;
        else if(key==OF_KEY_RIGHT)bpm+=2;
        else if(key==OF_KEY_LEFT)bpm-=2;
    }
    
    if(!paramMode)myfft.changeBandRange(key);
    else if(paramMode)myfft.changeParam(key);
    
    /*-----------LED--------------*/
    if(!paramMode){
        if(key=='a'){
            bpm=112;
            beat=1;
            nowTime=ofGetElapsedTimeMillis();
            float nextBeat = 1000/(bpm/60);
            targetTime=nowTime+nextBeat;
        }else if(key=='1'){
            if(ledMode!=1){
                ledMode=1;
                for(int i=0;i<4;i++)ard.sendDigitalPinMode(pin[i],ARD_OUTPUT);
            }
        }else if(key=='2'){
            if(ledMode!=2){
                ledMode=2;
                for(int i=0;i<4;i++)ard.sendDigitalPinMode(pin[i], ARD_PWM);
            }
        }else if(key=='3'){
            if(ledMode!=3){
                ledMode=3;
                for(int i=0;i<4;i++){
                    ard.sendDigitalPinMode(pin[i], ARD_OUTPUT);
                    ard.sendDigital(pin[i], ARD_LOW);
                }
            }
        }else if(key=='4'){
            if(ledMode!=4){
                ledMode=4;
                for(int i=0;i<4;i++){
                    ard.sendDigitalPinMode(pin[i], ARD_PWM);
                }
            }
        }
    }
}

//--------------------------------------------------------------
void ofApp::keyReleased(int key){
    
}

//--------------------------------------------------------------

void ofApp::mousePressed(int x, int y, int button){
    if((x>=230&&x<330)&&(y>=650&&y<700)){
        if(!myfft.bReset)myfft.bReset=true;
    }
    if((x>=360&&x<460)&&(y>=650&&y<700)){
        if(!myfft.bAutoMaxGet)myfft.bAutoMaxGet=true;
        else myfft.bAutoMaxGet=false;
    }
}

//--------------------------------------------------------------
void ofApp::mouseReleased(int x, int y, int button){
    
}

void ofApp::audioReceived 	(float * input, int bufferSize, int nChannels){
    // samples are "interleaved"
    for (int i = 0; i < bufferSize; i++){
        left[i] = input[i*2];
        right[i] = input[i*2+1];
    }
    bufferCounter++;
}