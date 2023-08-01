/* Audio Library for Teensy 3.X
 * Copyright (c) 2014, Paul Stoffregen, paul@pjrc.com
 *
 * Development of this audio library was funded by PJRC.COM, LLC by sales of
 * Teensy and Audio Adaptor boards.  Please support PJRC's efforts to develop
 * open source software by purchasing Teensy or other PJRC products.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice, development funding notice, and this permission
 * notice shall be included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

#ifndef AudioFilterBiquad_H_
#define AudioFilterBiquad_H_

#include <cstdint>
#include <Audio.h>
#include <arm_math.h>

#define IIR_F32_PASSTHRU ((const float32_t *) 1)

#define IIR_MAX_STAGES 1  //meaningless right now

class AudioFilterBiquad
{
public:
	AudioFilterBiquad(void) {
		// by default, the filter will not pass anything
		for (int i=0; i<32; i++) definition[i] = 0;
	}
	void update(audio_block_t* block);

	// Set the biquad coefficients directly
	void setCoefficients(uint32_t stage, const int *coefficients);
	void setCoefficients(uint32_t stage, const double *coefficients) {
		int coef[5];
		coef[0] = coefficients[0] * 1073741824.0;
		coef[1] = coefficients[1] * 1073741824.0;
		coef[2] = coefficients[2] * 1073741824.0;
		coef[3] = coefficients[3] * 1073741824.0;
		coef[4] = coefficients[4] * 1073741824.0;
		setCoefficients(stage, coef);
	}

	// Compute common filter functions
	// http://www.musicdsp.org/files/Audio-EQ-Cookbook.txt
	void setLowpass(uint32_t stage, float frequency, float q = 0.7071f) {
		int coef[5];
		double w0 = frequency * (2.0f * 3.141592654f / AUDIO_SAMPLE_RATE_EXACT);
		double sinW0 = sin(w0);
		double alpha = sinW0 / ((double)q * 2.0);
		double cosW0 = cos(w0);
		double scale = 1073741824.0 / (1.0 + alpha);
		/* b0 */ coef[0] = ((1.0 - cosW0) / 2.0) * scale;
		/* b1 */ coef[1] = (1.0 - cosW0) * scale;
		/* b2 */ coef[2] = coef[0];
		/* a1 */ coef[3] = (-2.0 * cosW0) * scale;
		/* a2 */ coef[4] = (1.0 - alpha) * scale;
		setCoefficients(stage, coef);
	}
	void setHighpass(uint32_t stage, float frequency, float q = 0.7071) {
		int coef[5];
		double w0 = frequency * (2.0f * 3.141592654f / AUDIO_SAMPLE_RATE_EXACT);
		double sinW0 = sin(w0);
		double alpha = sinW0 / ((double)q * 2.0);
		double cosW0 = cos(w0);
		double scale = 1073741824.0 / (1.0 + alpha);
		/* b0 */ coef[0] = ((1.0 + cosW0) / 2.0) * scale;
		/* b1 */ coef[1] = -(1.0 + cosW0) * scale;
		/* b2 */ coef[2] = coef[0];
		/* a1 */ coef[3] = (-2.0 * cosW0) * scale;
		/* a2 */ coef[4] = (1.0 - alpha) * scale;
		setCoefficients(stage, coef);
	}
	void setBandpass(uint32_t stage, float frequency, float q = 1.0) {
		int coef[5];
		double w0 = frequency * (2.0f * 3.141592654f / AUDIO_SAMPLE_RATE_EXACT);
		double sinW0 = sin(w0);
		double alpha = sinW0 / ((double)q * 2.0);
		double cosW0 = cos(w0);
		double scale = 1073741824.0 / (1.0 + alpha);
		/* b0 */ coef[0] = alpha * scale;
		/* b1 */ coef[1] = 0;
		/* b2 */ coef[2] = (-alpha) * scale;
		/* a1 */ coef[3] = (-2.0 * cosW0) * scale;
		/* a2 */ coef[4] = (1.0 - alpha) * scale;
		setCoefficients(stage, coef);
	}
	void setNotch(uint32_t stage, float frequency, float q = 1.0) {
		int coef[5];
		double w0 = frequency * (2.0f * 3.141592654f / AUDIO_SAMPLE_RATE_EXACT);
		double sinW0 = sin(w0);
		double alpha = sinW0 / ((double)q * 2.0);
		double cosW0 = cos(w0);
		double scale = 1073741824.0 / (1.0 + alpha);
		/* b0 */ coef[0] = scale;
		/* b1 */ coef[1] = (-2.0 * cosW0) * scale;
		/* b2 */ coef[2] = coef[0];
		/* a1 */ coef[3] = (-2.0 * cosW0) * scale;
		/* a2 */ coef[4] = (1.0 - alpha) * scale;
		setCoefficients(stage, coef);
	}
	void setLowShelf(uint32_t stage, float frequency, float gain, float slope = 1.0f) {
		int coef[5];
		double a = pow(10.0, gain/40.0f);
		double w0 = frequency * (2.0f * 3.141592654f / AUDIO_SAMPLE_RATE_EXACT);
		double sinW0 = sin(w0);
		//double alpha = (sinW0 * sqrt((a+1/a)*(1/slope-1)+2) ) / 2.0;
		double cosW0 = cos(w0);
		//generate three helper-values (intermediate results):
		double sinsq = sinW0 * sqrt( (pow(a,2.0)+1.0)*(1.0/(double)slope-1.0)+2.0*a );
		double aMinus = (a-1.0)*cosW0;
		double aPlus = (a+1.0)*cosW0;
		double scale = 1073741824.0 / ( (a+1.0) + aMinus + sinsq);
		/* b0 */ coef[0] =		a *	( (a+1.0) - aMinus + sinsq	) * scale;
		/* b1 */ coef[1] =  2.0*a * ( (a-1.0) - aPlus  			) * scale;
		/* b2 */ coef[2] =		a * ( (a+1.0) - aMinus - sinsq 	) * scale;
		/* a1 */ coef[3] = -2.0*	( (a-1.0) + aPlus			) * scale;
		/* a2 */ coef[4] =  		( (a+1.0) + aMinus - sinsq	) * scale;
		setCoefficients(stage, coef);
	}
	void setHighShelf(uint32_t stage, float frequency, float gain, float slope = 1.0f) {
		int coef[5];
		double a = pow(10.0, gain/40.0f);
		double w0 = frequency * (2.0f * 3.141592654f / AUDIO_SAMPLE_RATE_EXACT);
		double sinW0 = sin(w0);
		//double alpha = (sinW0 * sqrt((a+1/a)*(1/slope-1)+2) ) / 2.0;
		double cosW0 = cos(w0);
		//generate three helper-values (intermediate results):
		double sinsq = sinW0 * sqrt( (pow(a,2.0)+1.0)*(1.0/(double)slope-1.0)+2.0*a );
		double aMinus = (a-1.0)*cosW0;
		double aPlus = (a+1.0)*cosW0;
		double scale = 1073741824.0 / ( (a+1.0) - aMinus + sinsq);
		/* b0 */ coef[0] =		a *	( (a+1.0) + aMinus + sinsq	) * scale;
		/* b1 */ coef[1] = -2.0*a * ( (a-1.0) + aPlus  			) * scale;
		/* b2 */ coef[2] =		a * ( (a+1.0) + aMinus - sinsq 	) * scale;
		/* a1 */ coef[3] =  2.0*	( (a-1.0) - aPlus			) * scale;
		/* a2 */ coef[4] =  		( (a+1.0) - aMinus - sinsq	) * scale;
		setCoefficients(stage, coef);
	}

private:
	int32_t definition[32];  // up to 4 cascaded biquads
	//audio_block_t *inputQueueArray[1];
};


class AudioFilterBiquadFloat
{
  //GUI: inputs:1, outputs:1  //this line used for automatic generation of GUI node
  //GUI: shortName:IIR
  public:
    AudioFilterBiquadFloat(void): coeff_p(IIR_F32_PASSTHRU) {
    }
    void begin(const float32_t *cp, int n_stages) {
      coeff_p = cp;
      // Initialize FIR instance (ARM DSP Math Library)
      if (coeff_p && (coeff_p != IIR_F32_PASSTHRU) && n_stages <= IIR_MAX_STAGES) {
        //https://www.keil.com/pack/doc/CMSIS/DSP/html/group__BiquadCascadeDF1.html
        arm_biquad_cascade_df1_init_f32(&iir_inst, n_stages, (float32_t *)coeff_p,  &StateF32[0]);
      }
    }
    void end(void) {
      coeff_p = NULL;
    }

    void setBlockDC(void) {
      //https://www.keil.com/pack/doc/CMSIS/DSP/html/group__BiquadCascadeDF1.html#ga8e73b69a788e681a61bccc8959d823c5
      //Use matlab to compute the coeff for HP at 40Hz: [b,a]=butter(2,40/(44100/2),'high'); %assumes fs_Hz = 44100
      float32_t b[] = {8.173653471988667e-01,  -1.634730694397733e+00,   8.173653471988667e-01};  //from Matlab
      float32_t a[] = { 1.000000000000000e+00,   -1.601092394183619e+00,  6.683689946118476e-01};  //from Matlab
      setFilterCoeff_Matlab(b, a);
    }

	void setCoefficients(unsigned stage, double* coefficients) {
		//https://www.keil.com/pack/doc/CMSIS/DSP/html/group__BiquadCascadeDF1.html#ga8e73b69a788e681a61bccc8959d823c5
        //Use matlab to compute the coeff, such as: [b,a]=butter(2,20/(44100/2),'high'); %assumes fs_Hz = 44100
		hp_coeff[0] = coefficients[0];
		hp_coeff[1] = coefficients[1];
		hp_coeff[2] = coefficients[2];
		//here are the matlab "b" coefficients
		hp_coeff[3] = -coefficients[3];
		hp_coeff[4] = -coefficients[4];  //the DSP needs the "a" terms to have opposite sign vs Matlab
		uint8_t n_stages = 1;
		arm_biquad_cascade_df1_init_f32(&iir_inst, n_stages, hp_coeff,  &StateF32[0]);
		coeff_p = hp_coeff;
	}

    void setFilterCoeff_Matlab(float32_t b[], float32_t a[]) { //one stage of N=2 IIR
      //https://www.keil.com/pack/doc/CMSIS/DSP/html/group__BiquadCascadeDF1.html#ga8e73b69a788e681a61bccc8959d823c5
      //Use matlab to compute the coeff, such as: [b,a]=butter(2,20/(44100/2),'high'); %assumes fs_Hz = 44100
      hp_coeff[0] = b[0];   hp_coeff[1] = b[1];  hp_coeff[2] = b[2]; //here are the matlab "b" coefficients
      hp_coeff[3] = -a[1];  hp_coeff[4] = -a[2];  //the DSP needs the "a" terms to have opposite sign vs Matlab
      uint8_t n_stages = 1;
      arm_biquad_cascade_df1_init_f32(&iir_inst, n_stages, hp_coeff,  &StateF32[0]);
	  coeff_p = hp_coeff;
    }

    virtual void update(float* block, size_t numSamples);

  private:
    float32_t hp_coeff[5 * 1] = {1.0, 0.0, 0.0, 0.0, 0.0}; //no filtering. actual filter coeff set later

    // pointer to current coefficients or NULL or FIR_PASSTHRU
    const float32_t *coeff_p;

    // ARM DSP Math library filter instance
    arm_biquad_casd_df1_inst_f32 iir_inst;
    float32_t StateF32[4*IIR_MAX_STAGES];
};

#endif