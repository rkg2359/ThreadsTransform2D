// Threaded two-dimensional Discrete FFT transform
// YOUR NAME HERE
// ECE8893 Project 2


#include <iostream>
#include <string>
#include <math.h>
#include <stdint.h>

#include "Complex.h"
#include "InputImage.h"


// You will likely need global variables indicating how
// many threads there are, and a Complex* that points to the
// 2d image being transformed.

using namespace std;

//////////////////////GLOBAL VARIABLES////////////////////////////////////
pthread_mutex_t imageMutex, exitMutex, executingMutex, contMutex;
pthread_cond_t exitCond;
int numThreads = 16;
int activeThreads = 0;

Complex* imageData;
Complex* weightVals;
unsigned imageHeight, //Image Height 
          imageWidth, //Image Width
          imageSize;  //Image Size
unsigned threadCount;

pthread_barrier_t barrier; //Barrier declaration

// Function to reverse bits in an unsigned integer
// This assumes there is a global variable N that is the
// number of points in the 1D transform.
unsigned reverseBits(unsigned v)
{ //  Provided to students
  unsigned n = imageWidth; // Size of array (which is even 2 power k value)
  unsigned r = 0; // Return value
   
  for (--n; n > 0; n >>= 1)
    {
      r <<= 1;        // Shift return value
      r |= (v & 0x1); // Merge in next bit
      v >>= 1;        // Shift reversal value
    }
  return r;
}

// GRAD Students implement the following 2 functions.
// Undergrads can use the built-in barriers in pthreads.

// Call MyBarrier_Init once in main
// void MyBarrier_Init()// you will likely need some parameters)
// {
// }

// Each thread calls MyBarrier after completing the row-wise DFT
// void MyBarrier() // Again likely need parameters
// {
// }
                    
// void Transform1D(Complex* h, int N)
// {
//   // Implement the efficient Danielson-Lanczos DFT here.
//   // "h" is an input/output parameter
//   // "N" is the size of the array (assume even power of 2)
// }
void reorderRow(Complex* dataArr, int w) {
  int bits;
  for(int i = 0; i < w; i++)
  {
    //Get bit reversed value and put in order
    bits = reverseBits(i);
    if(bits > i)
    {
      pthread_mutex_lock(&imageMutex);
      Complex t = dataArr[i]; //Swapping Values
      dataArr[i] = dataArr[bits];
      dataArr[bits] = t;
      pthread_mutex_unlock(&imageMutex);
    }
  }
}

void reorderCol(Complex * dataArr, int height, int width)
{
  //We need to reorder the array to use the D-L DFT
  int bits;
  for(int i = 0; i < height; i++)
  {
    //Get bit reversed value and put in order
    bits = reverseBits(i); 
    if(bits > i)
    {
      pthread_mutex_lock(&imageMutex);
      //Swap the values
      Complex t = dataArr[i * width];
      dataArr[i * width] = dataArr[bits * width];
      dataArr[bits * width] = t;
      pthread_mutex_unlock(&imageMutex);
    }
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
        
        pthread_mutex_lock(&imageMutex);
        Complex temp = weightVals[w_ind] * h[k];
        h[k] = h[j] - temp;
        h[j] = h[j] + temp;       
        pthread_mutex_unlock(&imageMutex);
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
        
        pthread_mutex_lock(&imageMutex);
        Complex temp = weightVals[w_ind] * h[k*width];
        h[k*width] = h[j*width] - temp;
        h[j*width] = h[j*width] + temp;       
        pthread_mutex_unlock(&imageMutex);
      }
    }
  } 
}
//Tranform Rows
void transformAllRows(Complex * data, int start, int width, int rowsPerThread){

  //int k, weightIndex;
  for(int i = start; i < start + rowsPerThread; i++)
    {
      //Reorder Rows
      reorderRow(&data[i*width], width);

      //1-D Transform for rows
      Transform1D_row(&data[i*width], width);
      // for(int x = 1; x < width; x *= 2)
      // {
      //   for(int i = 0; i < x; i++)
      //   {
      //     weightIndex = i * width / (x*2);
          
      //     for(int j = i; j < width; j += (x * 2))
      //     {
      //       k = j + x;
            
      //       pthread_mutex_lock(&imageMutex);
      //       Complex temp = weightVals[weightIndex] * data[k];
      //       data[k] = data[j] - temp;
      //       data[j] = data[j] + temp;       
      //       pthread_mutex_unlock(&imageMutex);
      //     }
      //   }
      // }
    }

}

//Transform Columns
void transformAllCol(Complex * data, int start, int height, int width, int rowsPerThread){
  //int k, weightIndex;
  for(int i = start; i < start + rowsPerThread; i++)
    {
      reorderCol(&data[i], height, width);
      Transform1D_col(&data[i], height, width);
      // for(int x = 1; x < width; x *= 2)
      // {
      //   for(int i = 0; i < x; i++)
      //   {
      //     weightIndex = i * width / (x*2);
          
      //     for(int j = i; j < width; j += (x * 2))
      //     {
      //       k = j + x;
            
      //       pthread_mutex_lock(&imageMutex);
      //       Complex temp = weightVals[weightIndex] * data[k * width];
      //       data[k * width] = data[j * width] - temp;
      //       data[j * width] = data[j * width] + temp;       
      //       pthread_mutex_unlock(&imageMutex);
      //     }
      //   }
      // }

    }

}

void* Transform2DThread(void* v)
{ // This is the thread startign point.  "v" is the thread number

  int rank = *((int*)(&v)); //casting Void to int
  int rowsPerThread = imageHeight / numThreads;
  int start = rowsPerThread * rank;
  // Calculate 1d DFT for assigned rows
  transformAllRows(imageData, start, imageWidth, rowsPerThread);

  //Wait for all threads
  pthread_barrier_wait(&barrier);

  // Calculate 1d DFT for assigned columns
  transformAllCol(imageData, start, imageHeight, imageWidth, rowsPerThread);

  // Decrement active count and signal main if all complete
  pthread_mutex_lock(&executingMutex);
  activeThreads--;
  if (activeThreads == 0) {
    pthread_mutex_lock(&exitMutex);
    pthread_cond_signal(&exitCond);
    pthread_mutex_unlock(&exitMutex);
  }
  pthread_mutex_unlock(&executingMutex);
  return 0;
}

void Transform2D(const char* inputFN) 
{ // Do the 2D transform here.
  InputImage image(inputFN);  // Create the helper object for reading the image
  // Create the global pointer to the image array data
  imageData = image.GetImageData();
  imageHeight = image.GetHeight();
  imageWidth = image.GetWidth();
  //imageSize = imageWidth * imageHeight;//delete this

  //PreCalculating Weight Values
  weightVals = new Complex[imageHeight / 2];
  //double real, imag;
  for(unsigned i = 0; i < imageWidth / 2; i++)
  {
    double real, imag;
    real = cos(2 * M_PI * i / imageWidth);
    imag = -1 * sin(2 * M_PI * i / imageWidth);
    
    weightVals[i].real = real;
    weightVals[i].imag = imag;
  }

  
  // Create 16 threads
  for(int i = 0; i < numThreads; i++)
  {
    pthread_mutex_lock(&executingMutex);
    activeThreads++;
    pthread_mutex_unlock(&executingMutex);
    pthread_t thread;
    pthread_create(&thread, 0, Transform2DThread, (void*)i);
  }
  // Wait for all threads complete
  pthread_cond_wait(&exitCond, &exitMutex);
 
  
  // Write the transformed data
  image.SaveImageData("MyAfter2D.txt", imageData, imageWidth, imageHeight);  
  
  delete imageData;
  delete weightVals;
}

int main(int argc, char** argv)
{
  string fn("Tower.txt"); // default file name
  if (argc > 1) fn = string(argv[1]);  // if name specified on cmd line


  //Barrier Initialization
  pthread_barrier_init(&barrier, 0, numThreads);

  //Initialize MUTEX & Conditions
  pthread_mutex_init(&imageMutex, 0);
  pthread_mutex_init(&contMutex, 0);
  pthread_mutex_init(&executingMutex, 0);

  pthread_mutex_init(&exitMutex, 0);
  //Lock the Exit Mutex
  pthread_mutex_lock(&exitMutex);


  // MPI initialization here
  Transform2D(fn.c_str()); // Perform the transform.
}  
  

  
