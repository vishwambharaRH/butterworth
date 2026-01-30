/**
 * @file test_coefficients.c
 * @brief Comprehensive validation test for iirdsp filter coefficients
 *
 * This program prints filter coefficients in SOS format to compare
 * with scipy.signal.butter() output.
 */

#include <stdio.h>
#include <math.h>
#include "iirdsp.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

void print_filter_sos(const iirdsp_filter_t* f, const char* name)
{
    printf("\n%s Filter SOS Coefficients:\n", name);
    printf("========================================\n");
    printf("Number of sections: %d\n\n", f->num_sections);
    
    for (int i = 0; i < f->num_sections; i++) {
        printf("Section %d:\n", i);
        printf("  b0 = % .15e\n", f->sections[i].b0);
        printf("  b1 = % .15e\n", f->sections[i].b1);
        printf("  b2 = % .15e\n", f->sections[i].b2);
        printf("  a0 =  1.000000000000000e+00\n");
        printf("  a1 = % .15e\n", f->sections[i].a1);
        printf("  a2 = % .15e\n", f->sections[i].a2);
        printf("\n");
    }
}

void test_impulse_response(iirdsp_filter_t* f, const char* name, int N)
{
    printf("\n%s Impulse Response (first 20 samples):\n", name);
    printf("========================================\n");
    
    iirdsp_filter_reset(f);
    
    for (int i = 0; i < N; i++) {
        iirdsp_real x = (i == 0) ? 1.0 : 0.0;
        iirdsp_real y = iirdsp_process_sample(f, x);
        
        if (i < 20) {
            printf("[%2d] = % .15e\n", i, y);
        }
    }
}

void compute_frequency_response(iirdsp_filter_t* f, const char* name, 
                                iirdsp_real fs_hz, int num_points)
{
    printf("\n%s Frequency Response:\n", name);
    printf("========================================\n");
    
    /* Compute at specific test frequencies */
    iirdsp_real test_freqs[] = {0.0, 0.5, 1.0, 5.0, 10.0, 20.0, 40.0, 50.0, 100.0, 200.0};
    int num_test_freqs = sizeof(test_freqs) / sizeof(test_freqs[0]);
    
    printf("Freq (Hz)    |H(f)| (dB)\n");
    printf("------------------------\n");
    
    for (int i = 0; i < num_test_freqs; i++) {
        iirdsp_real freq = test_freqs[i];
        if (freq > fs_hz / 2.0) continue;
        
        /* Compute frequency response at this frequency */
        iirdsp_real w = 2.0 * M_PI * freq / fs_hz;
        iirdsp_real h_re = 1.0;
        iirdsp_real h_im = 0.0;
        
        for (int j = 0; j < f->num_sections; j++) {
            /* Numerator */
            iirdsp_real num_re = f->sections[j].b0 + 
                                 f->sections[j].b1 * cos(w) + 
                                 f->sections[j].b2 * cos(2.0 * w);
            iirdsp_real num_im = -f->sections[j].b1 * sin(w) - 
                                 f->sections[j].b2 * sin(2.0 * w);
            
            /* Denominator */
            iirdsp_real den_re = 1.0 + 
                                 f->sections[j].a1 * cos(w) + 
                                 f->sections[j].a2 * cos(2.0 * w);
            iirdsp_real den_im = -f->sections[j].a1 * sin(w) - 
                                 f->sections[j].a2 * sin(2.0 * w);
            
            /* Complex division */
            iirdsp_real denom = den_re * den_re + den_im * den_im;
            iirdsp_real h_sec_re = (num_re * den_re + num_im * den_im) / denom;
            iirdsp_real h_sec_im = (num_im * den_re - num_re * den_im) / denom;
            
            /* Accumulate */
            iirdsp_real new_h_re = h_re * h_sec_re - h_im * h_sec_im;
            iirdsp_real new_h_im = h_re * h_sec_im + h_im * h_sec_re;
            h_re = new_h_re;
            h_im = new_h_im;
        }
        
        iirdsp_real mag = sqrt(h_re * h_re + h_im * h_im);
        iirdsp_real mag_db = 20.0 * log10(mag + 1e-12);
        
        printf("%8.2f    %10.6f\n", freq, mag_db);
    }
}

int main(void)
{
    printf("=========================================================\n");
    printf("iirdsp Filter Coefficient Validation Test\n");
    printf("=========================================================\n");
    printf("Compare these outputs with scipy.signal.butter() output\n");
    printf("=========================================================\n");
    
    iirdsp_filter_t filter;
    iirdsp_real fs = 500.0;
    
    /* Test 1: Low-pass filter */
    printf("\n\n");
    printf("TEST 1: Low-Pass Butterworth Filter\n");
    printf("Order: 4, Cutoff: 10 Hz, Fs: 500 Hz\n");
    
    if (butter_lowpass_init(&filter, 4, 10.0, fs) == 0) {
        print_filter_sos(&filter, "Low-Pass");
        test_impulse_response(&filter, "Low-Pass", 100);
        compute_frequency_response(&filter, "Low-Pass", fs, 100);
    } else {
        printf("ERROR: Failed to initialize low-pass filter\n");
    }
    
    /* Test 2: High-pass filter */
    printf("\n\n");
    printf("TEST 2: High-Pass Butterworth Filter\n");
    printf("Order: 2, Cutoff: 40 Hz, Fs: 500 Hz\n");
    
    if (butter_highpass_init(&filter, 2, 40.0, fs) == 0) {
        print_filter_sos(&filter, "High-Pass");
        test_impulse_response(&filter, "High-Pass", 100);
        compute_frequency_response(&filter, "High-Pass", fs, 100);
    } else {
        printf("ERROR: Failed to initialize high-pass filter\n");
    }
    
    /* Test 3: Band-pass filter */
    printf("\n\n");
    printf("TEST 3: Band-Pass Butterworth Filter\n");
    printf("Order: 4, Band: 0.5-40 Hz, Fs: 500 Hz\n");
    
    if (butter_bandpass_init(&filter, 4, 0.5, 40.0, fs) == 0) {
        print_filter_sos(&filter, "Band-Pass");
        test_impulse_response(&filter, "Band-Pass", 100);
        compute_frequency_response(&filter, "Band-Pass", fs, 100);
    } else {
        printf("ERROR: Failed to initialize band-pass filter\n");
    }
    
    /* Test 4: Notch filter */
    printf("\n\n");
    printf("TEST 4: Notch Filter\n");
    printf("Center: 50 Hz, Q: 30, Fs: 500 Hz\n");
    
    if (notch_filter_init(&filter, 50.0, 30.0, fs) == 0) {
        print_filter_sos(&filter, "Notch");
        test_impulse_response(&filter, "Notch", 100);
        compute_frequency_response(&filter, "Notch", fs, 100);
    } else {
        printf("ERROR: Failed to initialize notch filter\n");
    }
    
    /* Test 5: filtfilt test */
    printf("\n\n");
    printf("TEST 5: Zero-Phase Filtering (filtfilt)\n");
    printf("Band-Pass: 0.5-40 Hz, Fs: 500 Hz\n");
    printf("========================================\n");
    
    butter_bandpass_init(&filter, 4, 0.5, 40.0, fs);
    
    const int N = 100;
    iirdsp_real input[100];
    iirdsp_real output[100];
    
    /* Generate test signal: 1 Hz sine wave */
    for (int i = 0; i < N; i++) {
        input[i] = sin(2.0 * M_PI * 1.0 * i / fs);
    }
    
    /* Apply filtfilt */
    iirdsp_filtfilt(&filter, input, output, N);
    
    printf("\nFirst 10 samples (input, output):\n");
    for (int i = 0; i < 10; i++) {
        printf("[%2d] %10.6f -> %10.6f\n", i, input[i], output[i]);
    }
    
    /* Compute RMS */
    iirdsp_real rms_in = 0.0, rms_out = 0.0;
    for (int i = 0; i < N; i++) {
        rms_in += input[i] * input[i];
        rms_out += output[i] * output[i];
    }
    rms_in = sqrt(rms_in / N);
    rms_out = sqrt(rms_out / N);
    
    printf("\nRMS values:\n");
    printf("  Input:  %.6f\n", rms_in);
    printf("  Output: %.6f\n", rms_out);
    
    printf("\n\n");
    printf("=========================================================\n");
    printf("Test completed successfully!\n");
    printf("=========================================================\n");
    printf("\nTo validate against SciPy, run:\n");
    printf("  python3 tests/scipy_compare.py\n");
    printf("\nThen compare the SOS coefficients printed above with\n");
    printf("the SciPy output. They should match within numerical precision.\n");
    
    return 0;
}