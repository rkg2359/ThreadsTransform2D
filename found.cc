// Threaded two-dimensional Discrete FFT transform
// Phillip Johnston
// ECE8893 Project 2
// 23 September 2012


#include <iostream>
#include <string>
#include <math.h>
#include <stdint.h>

#include "Complex.h"
#include "InputImage.h"

using namespace std;

#define _USE_MATH_DEFINES
#define PI M_PI //just in case

//#define USE_MY_BARRIER

/************************
* Function Declarations *
************************/
void MyBarrier_Init();
void MyBarrier();
void Transform2D(const char* inputFN);
void Transform1D_row(Complex * h, uint32_t N);
void Transform1D_col(Complex * h, uint32_t N, uint32_t width);
void* Transform2DThread(void* v);
void transformMyRows(Complex * h, uint32_t my_start, uint32_t width, uint32_t rows_per);
void transformMyCols(Complex * h, uint32_t my_start, uint32_t height, uint32_t width, uint32_t cols_per);
void reorderArray_row(Complex * arr, uint32_t n);
void reorderArray_col(Complex * arr, uint32_t n, uint32_t width);
unsigned reverseBits(unsigned v);
void precalculateWeightVals(Complex * buff, uint32_t n);

/*******************
* Global Variables *
*******************/
Complex * h_img;
Complex * weight_vals;
uint32_t img_height, img_width, img_sz, thread_count;
pthread_mutex_t img_mutex, exit_mutex, executing_mutex, cout_mutex;
pthread_cond_t  exit_cond;
uint8_t num_threads = 16;
uint8_t still_executing = 0;
#ifndef USE_MY_BARRIER
pthread_barrier_t barrier;
#endif

/***********************
* Function Definitions *
***********************/
int main(int argc, char** argv)
{
	string fn("Tower.txt"); // default file name
	if (argc > 1) fn = string(argv[1]);  // if name specified on cmd line
	
	//Initialize the Barrier
	MyBarrier_Init();

	//Initialize mutexes & conditions
	pthread_mutex_init(&img_mutex, 0);
	pthread_mutex_init(&exit_mutex, 0);
	pthread_mutex_init(&cout_mutex, 0);
	pthread_mutex_init(&executing_mutex, 0);
	pthread_cond_init(&exit_cond, 0);

	//Lock the exit mutex
	pthread_mutex_lock(&exit_mutex);
	
	Transform2D(fn.c_str()); // Perform the transform.
	
	return 0;
}  

void Transform2D(const char* inputFN) 
{
	InputImage image(inputFN);  // Create the helper object for reading the image

	//Set global parameters / get data
	img_height = image.GetHeight();
	img_width = image.GetWidth();
	img_sz = img_height * img_width;
	h_img = image.GetImageData();
	
	//Precalculate weighting values to save time later
	weight_vals = new Complex[img_height / 2];
	precalculateWeightVals(weight_vals, img_width);
	
	// Create 16 threads
	for(uint8_t i = 0; i < num_threads; i++)
	{
		pthread_mutex_lock(&executing_mutex);
		still_executing++; //Executing thread count, for conditioned wait
		pthread_mutex_unlock(&executing_mutex);
		
		pthread_t pt;
		pthread_create(&pt, 0, Transform2DThread, (void*)i);
	}

	// Wait for all threads complete
	pthread_cond_wait(&exit_cond, &exit_mutex);
	
	// Write the transformed data
	image.SaveImageData("Tower-DFT2D.txt", h_img, 1024, 1024);  
	
	delete h_img; //For good memory management practice
	delete weight_vals;
}

// Call MyBarrier_Init once in main
void MyBarrier_Init()// you will likely need some parameters)
{
#ifdef USE_MY_BARRIER
	#error "Ya haven't written MyBarrier_Init() for this, big dummy."
#else
	//Pthreads
	if(pthread_barrier_init(&barrier, 0, num_threads))
	{
		cout << "Could not initialize pthreads barrier.\n";
		fflush(stdout);
		//exit(-1);
	}
#endif
}

// Each thread calls MyBarrier after completing the row-wise DFT
void MyBarrier() // Again likely need parameters
{
#ifdef USE_MY_BARRIER
	#error "Make me whole again, said MyBarrier()."
#else
	//Pthreads
	int32_t rc = pthread_barrier_wait(&barrier);
    if(rc != 0 && rc != PTHREAD_BARRIER_SERIAL_THREAD)
    {
		pthread_mutex_lock(&cout_mutex);
        cout << "Could not wait on barrier\n";
		fflush(stdout);
		pthread_mutex_unlock(&cout_mutex);
        //exit(-1);
    }
#endif
}

//Wrapper to perform 1D-DFT on a row basis for a set of rows
void transformMyRows(Complex * h, uint32_t my_start, uint32_t width, uint32_t rows_per)
{
	for(uint32_t i = my_start; i < my_start + rows_per; i++)
	{
		reorderArray_row(&h[i*width], width);
		Transform1D_row(&h[i*width], width);
	}
}

//Wrapper to perform 1D-DFT on a column basis for set of cols
void transformMyCols(Complex * h, uint32_t my_start, uint32_t height, uint32_t width, uint32_t rows_per)
{
	for(uint32_t i = my_start; i < my_start + rows_per; i++)
	{
		reorderArray_col(&h[i], height, width);
		Transform1D_col(&h[i], height, width);
	}
}

void Transform1D_row(Complex* h, uint32_t N)
{
	// Implement the efficient Danielson-Lanczos DFT here.
	// "h" is an input/output parameter
	// "N" is the size of the array (assume even power of 2)
	uint32_t k, w_ind;
	
	for(uint32_t xfrm_sz = 1; xfrm_sz < N; xfrm_sz *= 2)
	{
		for(uint32_t i = 0; i < xfrm_sz; i++)
		{
			w_ind = i * N / (xfrm_sz*2);
			
			for(uint32_t j = i; j < N; j += (xfrm_sz * 2))
			{
				k = j + xfrm_sz;
				
				pthread_mutex_lock(&img_mutex);
				Complex temp = weight_vals[w_ind] * h[k];
				h[k] = h[j] - temp;
				h[j] = h[j] + temp;				
				pthread_mutex_unlock(&img_mutex);
			}
		}
	}	
}

void Transform1D_col(Complex* h, uint32_t N, uint32_t width)
{
	// Implement the efficient Danielson-Lanczos DFT here.
	// "h" is an input/output parameter
	// "N" is the size of the array (assume even power of 2)
	uint32_t k, w_ind;
	
	for(uint32_t xfrm_sz = 1; xfrm_sz < N; xfrm_sz *= 2)
	{
		for(uint32_t i = 0; i < xfrm_sz; i++)
		{
			w_ind = i * N / (xfrm_sz*2);
			
			for(uint32_t j = i; j < N; j += (xfrm_sz * 2))
			{
				k = j + xfrm_sz;
				
				pthread_mutex_lock(&img_mutex);
				Complex temp = weight_vals[w_ind] * h[k*width];
				h[k*width] = h[j*width] - temp;
				h[j*width] = h[j*width] + temp;				
				pthread_mutex_unlock(&img_mutex);
			}
		}
	}	
}

void* Transform2DThread(void* v)
{ // This is the thread startign point.  "v" is the thread number
	// Calculate 1d DFT for assigned rows
	// wait for all to complete
	// Calculate 1d DFT for assigned columns
	// Decrement active count and signal main if all complete

	//Generate thread-relative data
	uint64_t t = (uint64_t) v;
	uint16_t rank = (uint16_t) t;
	uint32_t rows_per = img_height / num_threads;
	uint32_t my_start = rows_per * rank;
	
	transformMyRows(h_img, my_start, img_width, rows_per);

	//Enter the barrier and wait
	MyBarrier();
	
	transformMyCols(h_img, my_start, img_height, img_width, rows_per);

	//Finish and decrement the running count
	pthread_mutex_lock(&executing_mutex);
	
	//Check to see if we're the last thread to exit
	if(--still_executing == 0)
	{
		//Last to exit, signal main
		pthread_mutex_lock(&exit_mutex);
		pthread_cond_signal(&exit_cond);
		pthread_mutex_unlock(&exit_mutex);
	}
	
	pthread_mutex_unlock(&executing_mutex);

	return 0;
}

void precalculateWeightVals(Complex * buff, uint32_t n)
{
	//Precalculate weighting values
	//Implements (from 0 to N - 1 of e^-j2*pi*k/N) weights
	for(uint32_t i = 0; i < n; i++)
	{
		double real, imag;
		
		real = cos(2 * M_PI * i / n);
		imag = -1 * sin(2 * M_PI * i / n);
		
		buff[i].real = real;
		buff[i].imag = imag;
	}
}

/******************
* Reordering Code *
******************/
void reorderArray_row(Complex * arr, uint32_t n)
{
	//We need to reorder the array to use the D-L DFT
	unsigned j;
	for(unsigned i = 0; i < n; i++)
	{
		//Get bit reversed value and put in order
		j = reverseBits(i);
		if(j > i)
		{
			pthread_mutex_lock(&img_mutex);
			//Swap the values
			Complex t = arr[i];
			arr[i] = arr[j];
			arr[j] = t;
			pthread_mutex_unlock(&img_mutex);
		}
	}
}

//Column implementation
void reorderArray_col(Complex * arr, uint32_t n, uint32_t width)
{
	//We need to reorder the array to use the D-L DFT
	unsigned j;
	for(unsigned i = 0; i < n; i++)
	{
		//Get bit reversed value and put in order
		j = reverseBits(i); 
		if(j > i)
		{
			pthread_mutex_lock(&img_mutex);
			//Swap the values
			Complex t = arr[i * width];
			arr[i * width] = arr[j * width];
			arr[j * width] = t;
			pthread_mutex_unlock(&img_mutex);
		}
	}
}

// Function to reverse bits in an unsigned integer
// This assumes there is a global variable N that is the
// number of points in the 1D transform.
unsigned reverseBits(unsigned v)
{ //  Provided to students
	unsigned n = img_width; // Size of array (which is even 2 power k value)
	unsigned r = 0; // Return value

	for (--n; n > 0; n >>= 1)
	{
		r <<= 1;        // Shift return value
		r |= (v & 0x1); // Merge in next bit
		v >>= 1;        // Shift reversal value
	}
	return r;
}