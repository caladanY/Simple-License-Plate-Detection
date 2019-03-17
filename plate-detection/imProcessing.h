#pragma once
#include <vector>
#define PI 3.14

//Point
class Point {
public:
	int x;
	int y;
	Point(int x, int y) {
		this->x = x;
		this->y = y;
	}
	Point() {
		;
	}
};

//Rectangle
class Rect {
public:
	//Points can not only represent boundaries but vertices.
	Point leftBoundary;
	Point topBoundary;
	Point rightBoundary;
	Point bottomBoundary;
	double slope = 0;
	Rect(Point leftBoundary, Point topBoundary, Point rightBoundary, Point bottomBoundary) {
		this->topBoundary = topBoundary;
		this->leftBoundary = leftBoundary;
		this->bottomBoundary = bottomBoundary;
		this->rightBoundary = rightBoundary;
	}
	Rect(){
		;
	}
};

// Feature contains the key information of a connective region.
typedef struct _Feature{
public:
	int label;              // Label of the connective region
	int area;               // Area of the bounding rectangle
	Rect boundingbox;       // Bounding rectangle
} Feature;

// For most of time we deal with a grey scale image.
// To simplify the calculation, we only deal with 
// the first channel of an image. This function covers
// the second and the third channel with the value of 
// the first channel, so that the grey scale image can
// be written to a bmp file.
void grey2out(bitMap* img);

// Edge detection with Robert operation.
void Robert( bitMap* src, bitMap* des, int channel);

// Image Dilation
void Inflation(bitMap* img, bitMap *des, int temp_h, int temp_w);

// Image Erosion
void Erosion(bitMap* img, bitMap *des,int temp_h, int temp_w);

//------------------------Color Detection-----------------------

// Change RGB image to HSV image
void Rgb2Hsv(bitMap* img);//Changle RGB image to HSV image.

// Detect blue and yellow of the plate.
void HsvColorDetect(bitMap* src, bitMap* des);

// Draw a rectangle.
void rectangle(bitMap *src, Rect box, int color);

// Find connective regions.
int ConnectedRegion(bitMap* src, bitMap* des, std::vector<Feature> &featureList);

// Create a lookup table from x->sin(x)/cos(x)
void init_sincos();

// Detect the slope of the longest line in the region.
void lineDect(Rect &box, bitMap* src);

// Use the slope infomation to locate the vertices.
void countVertices(Rect &box);

//----------------------------Not used----------------------------

// Edge detection with Sobel operation.
void sobel(bitMap* src, bitMap* des, int channel);

// Use template to calculate the sum of the neighbour of a pixal
int TempltExcuteCl(IMAGEDATA* image, int width, int height, int* temp, 
	               int tempSize, int x, int y, int channel);

// Bilateral Filter
void Bilateral_Filter(bitMap *src, bitMap *des, int r, double sigma_d, 
	                  double sigma_r, int channels);

// Convert the image from RGB to BMP
void rgb2grey(bitMap *src, bitMap *des);

// Binaryzation a grey scale image
void binaryzation(bitMap *img);

// Get the threshold of the binaryzation.
int GetOSTUThreshold(int HistGram[]);

// Get the histogram.
void getHistGram(bitMap* img, int HistGram[]);
