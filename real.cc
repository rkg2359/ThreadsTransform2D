// Threaded two-dimensional Discrete FFT transform
// Rohit Ganesan
// ECE 4122


#include <iostream>
#include <string>
#include <math.h>


#include "Complex.h"
#include "InputImage.h"


// You will likely need global variables indicating how
// many threads there are, and a Complex* that points to the
// 2d image being transformed.

using namespace std;

////////////////////////////////GLOBAL VARIABLES////////////////////////////////////
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

/////////////////////////////////FUNCTIONS///////////////////////////////////

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

//reorder bit reversed values (Row)
void reorderRow(unsigned width, Complex* dataArr) {

  for(unsigned i = 0; i < width; i++)
  {
    if(reverseBits(i) > i)
    {//Lock, Swap the values, Unlock
      pthread_mutex_lock(&imageMutex);
      Complex t = dataArr[i]; //Swapping Values
      dataArr[i] = dataArr[reverseBits(i)];
      dataArr[reverseBits(i)] = t;
      pthread_mutex_unlock(&imageMutex);
    }
  }
}

void TransformRow1D(int width, Complex* h)
{
  int a, weightIndex;
  for(int x = 1; x < width; x *= 2)
  {
    for(int i = 0; i < x; i++)
    {
      weightIndex = i * width / (x*2);
      
      for(int j = i; j < width; j += (x * 2))
      {
        a = j + x;
        // Lock, change the data, Unlock
        pthread_mutex_lock(&imageMutex);
        Complex temp = weightVals[weightIndex] * h[a];
        h[a] = h[j] - temp;
        h[j] = h[j] + temp;
        pthread_mutex_unlock(&imageMutex);
      }
    }
  } 
}

//reorder bit reversed values (Columns)
void reorderCol(unsigned height, int width, Complex * dataArr)
{
  for(unsigned i = 0; i < height; i++)
  {
    if(reverseBits(i) > i)
    {//Lock, Swap the values, Unlock
      pthread_mutex_lock(&imageMutex);
      Complex t = dataArr[i * width];
      dataArr[i * width] = dataArr[reverseBits(i) * width];
      dataArr[reverseBits(i) * width] = t;
      pthread_mutex_unlock(&imageMutex);
    }
  }
}

void TransformCol1D(int N, int width, Complex* h)
{
  int a, weightIndex;
  
  for(int x = 1; x < N; x *= 2)
  {
    for(int i = 0; i < x; i++)
    {
      weightIndex = i * N / (x*2);
      
      for(int j = i; j < N; j += (x * 2))
      {
        a = j + x;
        // Lock, change the data, Unlock
        pthread_mutex_lock(&imageMutex);
        Complex temp = weightVals[weightIndex] * h[a*width];
        h[a*width] = h[j*width] - temp;
        h[j*width] = h[j*width] + temp;       
        pthread_mutex_unlock(&imageMutex);
      }
    }
  } 
}



//Tranform Rows
void transformAllRows(Complex * data, int start, int width, int rowsPerThread){

  for(int i = start; i < start + rowsPerThread; i++) //loop for all threads
    {
      //Reorder Rows
      reorderRow(width, &data[i*width]);

      //1-D Transform for rows
      TransformRow1D(width, &data[i*width]);
    }

}

//Transform Columns
void transformAllCol(Complex * data, int start, int height, int width, int rowsPerThread){

  for(int i = start; i < start + rowsPerThread; i++) //loop for all threads
    {
      //Reorder Columns
      reorderCol( height, width,&data[i]);
      
      //1-D Tranform for Columns
      TransformCol1D( height, width, &data[i]);
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
  activeThreads--; //Decrement
  if (activeThreads == 0) {
    pthread_mutex_lock(&exitMutex);
    pthread_cond_signal(&exitCond); //Signal Main
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
  
  //Delete variables
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
  

  
