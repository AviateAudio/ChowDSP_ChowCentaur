#ifndef FF1CURRENT_H_INCLUDED
#define FF1CURRENT_H_INCLUDED

#include <Audio.h>
#include "PreAmpStage.h"

class FF1Current
{
public:
    FF1Current() {}

    void setPreAmp(PreAmpStage* preAmpIn) { preAmp = preAmpIn; }
    void update(audio_block_t* block)
    {
        if (!block || !preAmp) return;

        int16_t i, y;
        float yf;
        for (i = 0; i < AUDIO_BLOCK_SAMPLES; i++) {
            yf = preAmp->ff1Current[i];
            y =  static_cast<int> (yf * 32768.0f);
            block->data[i] = 1000 * y;
        }
    }

private:
    PreAmpStage* preAmp = nullptr;
};

class FF1CurrentFloat
{
public:
    FF1CurrentFloat() {}

    void setPreAmp(PreAmpStageFloat* preAmpIn) { preAmp = preAmpIn; }
    void update(float* block, size_t numSamples)
    {
        if (!block || !preAmp) return;

        //int16_t i, y;
        float yf;
        for (unsigned i = 0; i < numSamples; i++) {
            yf = preAmp->ff1Current[i];
            //y =  static_cast<int> (yf * 32768.0f);
            //block->data[i] = 1000 * y;
            block[i] = 1000 * yf;
        }
    }

private:
    PreAmpStageFloat* preAmp = nullptr;
};

#endif // FF1CURRENT_H_INCLUDED
