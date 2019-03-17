#include <stdlib.h>
#include <iostream>
#include <stack>
#include <vector>
#include<cmath>   
#include "imread.h"
#include "imProcessing.h"
using namespace std;

float mysin[360], mycos[360];// Lookup table of sin/cos function.

// For most of time we deal with a grey scale image.
// To simplify the calculation, we only deal with 
// the first channel of an image. This function covers
// the second and the third channel with the value of 
// the first channel, so that the grey scale image can
// be written to a bmp file.
void grey2out(bitMap *img) {
    int width = img->width;
    int height = img->height;
    for (int x = 0; x < img->height; x++) {
        for (int y = 0; y < img->width; y++) {
            img->imagedata[x * width + y].channel[2] = img->imagedata[x * width + y].channel[1] = img->imagedata[x * width + y].channel[0];
        }
    }
    return;
}

// Edge detection with Robert operation.
void Robert(bitMap* src, bitMap* des, int channel) {
    int x, y, n;
    int w = src->width;
    int h = src->height;
    for (x = 1; x < src->height - 1; x++)
    {
        for (y = 1; y < src->width - 1; y++)
        {
            int weight = 0;
            for (n = 0; n < channel; n++) {
                weight += sqrt((double)((src->imagedata[x * w + y].channel[n] -
										 src->imagedata[(x + 1) * w + (y + 1)].channel[n]) *
										(src->imagedata[x * w + y].channel[n] -
										 src->imagedata[(x + 1) * w + (y + 1)].channel[n]) +
                                        (src->imagedata[x * w + (y + 1)].channel[n] -
										 src->imagedata[(x - 1) * w + y].channel[n]) *
										(src->imagedata[x * w + (y + 1)].channel[n] - 
										 src->imagedata[(x - 1) * w + y].channel[n])));
            }
            weight = weight > 255 ? 255 : weight;
            des->imagedata[x * w + y].channel[0] = des->imagedata[x * w + y].channel[1] = des->imagedata[x * w + y].channel[2] = weight;
        }
    }
    return;
}

// Image Erosion
void Erosion(bitMap* img, bitMap *des, int temp_h, int temp_w) {
    int r_h = temp_h / 2;
    int r_w = temp_w / 2;
    int *temp = new int[temp_w * temp_h] { 0 };
    for (int i = 0; i < temp_h; i++) {
        temp[i * temp_w] = 1;
    }
    for (int x = r_h; x < img->height - r_h; x++) {
        for (int y = r_w; y < img->width - r_w; y++) {
            int sig = 0;
            for (int j = 0; j < temp_w; j++)
                for (int i = 0; i < temp_h; i++) {
                    if (img->imagedata[(x + i - r_h) * des->width + y + j - r_w].channel[0] == 0) sig = 1;
                }
            //Judgement
            if (sig == 1) {
                des->imagedata[x * des->width + y].channel[0] = 0;
            }
        }
    }
    return;
}

// Image Dilation
void Inflation(bitMap* img, bitMap *des, int temp_h, int temp_w) {
    int r_h = temp_h / 2;
    int r_w = temp_w / 2;
    int *temp = new int[temp_w * temp_h] { 0 };
    for (int i = 0; i < temp_h; i++) {
        temp[i * temp_w] = 1;
    }
    for (int x = r_h; x < img->height - r_h; x++) {
        for (int y = r_w; y < img->width - r_w; y++) {
            int sig = 0;
            for (int j = 0; j < temp_w; j++)
                for (int i = 0; i < temp_h; i++) {
                    if (img->imagedata[(x + i - r_h) * des->width + y + j - r_w].channel[0] == 255) sig = 1;
                }
            //Judgement
            if (sig == 1) {
                des->imagedata[x * des->width + y].channel[0] = 255;
            }
        }
    }
    return;
}

//Change RGB image to HSV image
void Rgb2Hsv(bitMap* img){
	float H = 0; 
	float S = 0; 
	float V = 0;
	for (int i = 0; i < img->height; i++) {
		for (int j = 0; j < img->width; j++) {
			H = 0;
			S = 0;
			V = 0;
			float B = img->imagedata[i * img->width + j].channel[0] / 255.0;
			float G = img->imagedata[i * img->width + j].channel[1] / 255.0;
			float R = img->imagedata[i * img->width + j].channel[2] / 255.0;

			// r,g,b values are from 0 to 1  
			// h = [0,180], s = [0,1], v = [0,1]  
			// if s == 0, then h = -1 (undefined)  
			float min, max, delta, tmp;
			tmp = R>G ? G : R;
			min = tmp>B ? B : tmp;
			tmp = R>G ? R : G;
			max = tmp>B ? tmp : B;
			V = max; // v  
			delta = max - min;
			if (max != 0)
				S = delta / max; // s  
			else
			{
				// r = g = b = 0 // s = 0, v is undefined  
				S = 0;
				H = 0;
				img->imagedata[i * img->width + j].channel[0] = H;
				img->imagedata[i * img->width + j].channel[1] = S;
				img->imagedata[i * img->width + j].channel[2] = V * 255;
				continue;
			}
			if (delta == 0) {
				H = 0;
				img->imagedata[i * img->width + j].channel[0] = H;
				img->imagedata[i * img->width + j].channel[1] = S * 255;
				img->imagedata[i * img->width + j].channel[2] = V * 255;
				continue;
			}
			else if (R == max) {
				if (G >= B)
					H = (G - B) / delta; // between yellow & magenta  
				else
					H = (G - B) / delta + 6.0;
			}
			else if (G == max)
				H = 2.0 + (B - R) / delta; // between cyan & yellow  
			else if (B == max)
				H = 4.0 + (R - G) / delta; // between magenta & cyan  
			H *= 60.0; // degrees  
			img->imagedata[i * img->width + j].channel[0] = H / 2;
			img->imagedata[i * img->width + j].channel[1] = S * 255;
			img->imagedata[i * img->width + j].channel[2] = V * 255;
		}
	}
	
}

// Detect blue and yellow of the plate.
void HsvColorDetect(bitMap* src, bitMap* des) {
	for (int i = 0; i < src->height; i++) {
		for (int j = 0; j < src->width; j++) {
			float H = src->imagedata[i * src->width + j].channel[0];
			float S = src->imagedata[i * src->width + j].channel[1];
			float V = src->imagedata[i * src->width + j].channel[2];
			// Detect Blue
			if ((2 * H) > 201 && (2 * H) < 255 &&
				S > 0.35 * 255 && S < 255 &&
				V > 0.4 * 255 && V < 255) {
				des->imagedata[i * src->width + j].channel[0] = 255;
			}
			// Detect Yellow 
			else if ((2 * H) > 25 && (2 * H) < 55 &&
				S > 0.4 * 255 && S < 255 &&
				V > 0.5 * 255 && V < 230) {
				des->imagedata[i * src->width + j].channel[0] = 255;
			}
			else {
				des->imagedata[i * src->width + j].channel[0] = 0;
			}
		}
	}
	return;
}

// Draw a rectangle.
void rectangle(bitMap* src, Rect box, int color) {
	int r = 2;
	float slopeH = (float)(box.bottomBoundary.y - box.leftBoundary.y) / (box.bottomBoundary.x - box.leftBoundary.x);
	float slopeV = (float)(box.topBoundary.x - box.leftBoundary.x) / (box.leftBoundary.y - box.topBoundary.y);
	int u = src->imagedata[205 * src->width + 202].channel[2];
	for (int i = box.leftBoundary.x; i < box.bottomBoundary.x; i++) {
		int j = box.leftBoundary.y + (i - box.leftBoundary.x) * slopeH;
			for (int k = j - r; k < j + r; k++){
			src->imagedata[k * src->width + i].channel[0] = 0;
			src->imagedata[k * src->width + i].channel[1] = 0;
			src->imagedata[k * src->width + i].channel[2] = 255;
		}
	}
	for (int i = box.topBoundary.x; i < box.rightBoundary.x; i++) {
		int j = box.topBoundary.y + (i - box.topBoundary.x) * slopeH;
		for (int k = j - r; k < j + r; k++) {
			src->imagedata[k * src->width + i].channel[0] = 0;
			src->imagedata[k * src->width + i].channel[1] = 0;
			src->imagedata[k * src->width + i].channel[2] = 255;
		}
	}
	for (int j = box.leftBoundary.y; j > box.topBoundary.y; j--) {
		int i = box.leftBoundary.x + (box.leftBoundary.y - j) * slopeV;
		for (int k = i - r; k < i + r; k++) {
			src->imagedata[j * src->width + k].channel[0] = 0;
			src->imagedata[j * src->width + k].channel[1] = 0;
			src->imagedata[j * src->width + k].channel[2] = 255;
		}
	}
	for (int j = box.bottomBoundary.y; j > box.rightBoundary.y; j--) {
		int i = box.rightBoundary.x + (box.rightBoundary.y - j) * slopeV;
		for (int k = i - r; k < i + r; k++) {
			src->imagedata[j * src->width + k].channel[0] = 0;
			src->imagedata[j * src->width + k].channel[1] = 0;
			src->imagedata[j * src->width + k].channel[2] = 255;
		}
	}
	return;
}

// Find connective regions.
int ConnectedRegion(bitMap* src, bitMap*des, vector<Feature> &featureList) {
	int width = src->width;
	int height = src->height;
	int labelValue = 0;
	int area = 0;               // Calculate the area of the region
	Point leftBoundary(0, 0);;      // Left boundary of the region
	Point rightBoundary(0, 0);
	Point topBoundary(0, 0);
	Point bottomBoundary(0, 0);
	Rect box;                   // Bounding rectangle
	Feature feature;
	Point seed, neighbor;
	stack<Point> pointStack;    // Stack to contain points

	for (int i = 0; i < src->height; i++) {
		for (int j = 0; j < src->width; j++) {
			if (des->imagedata[i * src->width + j].channel[0] == 255) {
				area = 0;
				labelValue++;           // 0<labelValue<25256
				seed = Point(j, i);     // Point（x，y）
				des->imagedata[i * src->width + j].channel[0] = labelValue;
				pointStack.push(seed);

				// Get initial boundary
				area++;
				leftBoundary.x = seed.x;
				leftBoundary.y = seed.y;
				rightBoundary.x = seed.x;
				rightBoundary.y = seed.y;
				topBoundary.x = seed.x;
				topBoundary.y = seed.y;
				bottomBoundary.x = seed.x;
				bottomBoundary.y = seed.y;

				// Start iteration
				while (!pointStack.empty())
				{
					//Check the right neighbor
					neighbor = Point(seed.x + 1, seed.y);
					if ((seed.x != (width - 1)) && (des->imagedata[neighbor.y * src->width + neighbor.x].channel[0] == 255))
					{
						des->imagedata[neighbor.y * src->width + neighbor.x].channel[0] = labelValue;
						pointStack.push(neighbor);

						area++;
						if (rightBoundary.x < neighbor.x) {
							rightBoundary.x = seed.x;
							rightBoundary.y = seed.y;
						}
					}

					//Check the bottom neighbor
					neighbor = Point(seed.x, seed.y + 1);
					if ((seed.y != (height - 1)) && (des->imagedata[neighbor.y * src->width + neighbor.x].channel[0] == 255))
					{
						des->imagedata[neighbor.y * src->width + neighbor.x].channel[0] = labelValue;
						pointStack.push(neighbor);

						area++;
						if (bottomBoundary.y < neighbor.y) {
							bottomBoundary.x = neighbor.x;
							bottomBoundary.y = neighbor.y;
						}
					}

					//CHeck the left neighbor
					neighbor = Point(seed.x - 1, seed.y);
					if ((seed.x != 0) && (des->imagedata[neighbor.y * src->width + neighbor.x].channel[0] == 255))
					{
						des->imagedata[neighbor.y * src->width + neighbor.x].channel[0] = labelValue;
						pointStack.push(neighbor);

						area++;
						if (leftBoundary.x > neighbor.x) {
							leftBoundary.x = neighbor.x;
							leftBoundary.y = neighbor.y;
						}
					}

					//Check the top neighbor
					neighbor = Point(seed.x, seed.y - 1);
					if ((seed.y != 0) && (des->imagedata[neighbor.y * src->width + neighbor.x].channel[0] == 255))
					{
						des->imagedata[neighbor.y * src->width + neighbor.x].channel[0] = labelValue;
						pointStack.push(neighbor);

						area++;
						if (topBoundary.y > neighbor.y) {
							topBoundary.y = neighbor.y;
							topBoundary.x = neighbor.x;
						}
					}

					seed = pointStack.top();
					pointStack.pop();
				}
				box = Rect(leftBoundary, topBoundary, rightBoundary, bottomBoundary);
				feature.area = area;
				feature.boundingbox = box;
				feature.label = labelValue;
				featureList.push_back(feature);
			}
		}
	}
	return labelValue;
}

//-------------------------Line Detection----------------------
// Create a lookup table from x->sin(x)/cos(x)
void init_sincos(){
	for (int i = 0; i<360; i++)
	{
		mysin[i] = sin(float(i)*PI / 180.0);
		mycos[i] = cos(float(i)*PI / 180.0);
	}
}

// Find if a is adjacent to b
bool adjacent(Point a, Point b)
{
	if (abs(a.x - b.x) <= 10 && abs(b.y - a.y) <= 10)
		return 1;
	return 0;
}

// Detect the slope of the longest line in the region.
void lineDect(Rect &box, bitMap *src) {
	init_sincos();
	int scale = 1200;
	//Detect left top
	double r, theta; // Hough parameters
	double finR, finTheta;
	int finalMax;
	finalMax = 0;
	finR = 0;
	finTheta = 0;
	int intr;
	int RTnum[1200];
	for (int k = 0; k <= 180; k++)
	{
		//Initial RTnum
		for ( int n = 0; n < scale; n++)
			RTnum[n] = 0;

		//Find the side with angle theta, and  the length of the side is r               
		for (int i = box.topBoundary.y; i < box.bottomBoundary.y; i++)//i = rows 
		{
			for (int j = box.leftBoundary.x; j<box.rightBoundary.x; j++)//j = colls   
			{
				if (src->imagedata[i * src->width + j].channel[0] != 255) {
					continue;
				}
				r = mysin[k] * float(i) + mycos[k] * float(j);
				intr = fabs(r);
				RTnum[intr]++;
			}
		}

		//Find lines  
		int pnum;//pnum counts the length of the line
		int max = 0;
		double maxR;
		if (k == 180)
			int q = 0;
		for (int m = 0; m < scale; m++){
			if (m == 1200)
				int p = 0;
			if (RTnum[m] > max) {
				max = RTnum[m];
				maxR = m;
			}
		}
		if (max > finalMax) {
			finalMax = max;
			finR = maxR;
			finTheta = k;
		}
	}
	box.slope = tan((finTheta - 90) * PI / 180);
	return;
}

// Use the slope infomation to locate the vertices.
void countVertices(Rect &box) {
	int dx;
	int dy;
	int x;
	int y;
	float k = box.slope;
	x = box.rightBoundary.x - box.leftBoundary.x;
	y = box.bottomBoundary.y - box.topBoundary.y;
	dx = (k * y - k * k * x) / (1 - k * k);
	dy = (k * x - k * k * y) / (1 - k * k);
	//Consider k = 1;
	if ((k<1.1 && k>0.9) || (k<-0.9 && k > -1.1)) {
		dx = x >> 1;
		dy = y >> 1;
	}
	box.leftBoundary.y = box.bottomBoundary.y - dy;
	box.bottomBoundary.x =box.rightBoundary.x - dx;
	box.topBoundary.x = box.leftBoundary.x + dx;
	box.rightBoundary.y = box.topBoundary.y + dy;
	return;
}

//------------------------Not used function------------------------

//Use template to calculate the sum of the neighbourhood of a pixal
int TempltExcuteCl(IMAGEDATA *image, int width, int height, 
				   int* temp, int tempSize, int x, int y, int channel)
{
	int sum = 0;
	int px, py;
	//Calculate the pixal value in the neighbourhood
	for (int i = 0; i<tempSize; i++)
	{
		for (int j = 0; j<tempSize; j++)
		{
			// Get the responding situlation
			py = y - tempSize / 2 + j;
			px = x - tempSize / 2 + i;
			// Weighted summation
			sum += image[px*width + py].channel[channel] * temp[i*tempSize + j];
		}
	}
	return sum;
}

void Bilateral_Filter(bitMap *src, bitMap *des, int r, 
					  double sigma_d, double sigma_r, int channels)
{
	int i, j, m, n;
	int w_filter = 2 * r + 1; // Scale of the templet.
	int width = src->width;
	int height = src->height;

	// Calculate gaussian parameter
	double gaussian_d_coeff = -0.5 / (sigma_d * sigma_d);
	double gaussian_r_coeff = -0.5 / (sigma_r * sigma_r);

	double* d_metrix = new double[w_filter * w_filter];  // Spatial Weight
	double r_metrix[256];  // Similarity Weight

	// Compute spatial weight
	for (i = -r; i <= r; i++)
		for (j = -r; j <= r; j++)
		{
			int x = j + r;
			int y = i + r;

			d_metrix[y * w_filter + x] = exp((i * i + j * j) * gaussian_d_coeff);
		}

	// Compute similarity weight
	for (i = 0; i < 256; i++)
	{
		r_metrix[i] = exp(i * i * gaussian_r_coeff);
	}

	// Bilateral Filter
	for (i = 0; i < src->height; i++)
		for (j = 0; j < src->width; j++)
		{
			for (int channel = 0; channel < channels; channel++)
			{
				double weight_sum, pixcel_sum;
				weight_sum = pixcel_sum = 0.0;

				for (m = -r; m <= r; m++)
				{
					for (n = -r; n <= r; n++)
					{
						if (m*m + n*n > r*r) continue;

						int x_tmp = i + m;
						int y_tmp = j + n;

						//Boundary limit
						x_tmp = x_tmp < 0 ? 0 : x_tmp;
						x_tmp = x_tmp > src->height - 1 ? src->height - 1 : x_tmp;
						y_tmp = y_tmp < 0 ? 0 : y_tmp;
						y_tmp = y_tmp > src->width - 1 ? src->width - 1 : y_tmp;

						// Calculate pixcel difference
						int pixcel_dif = (int)abs(des->imagedata[x_tmp * width + y_tmp].channel[channel] 
												  - des->imagedata[i *  width + j].channel[channel]);
						double weight_tmp = d_metrix[(m + r) * w_filter + (n + r)] * r_metrix[pixcel_dif];
						pixcel_sum += des->imagedata[x_tmp * width + y_tmp].channel[channel] * weight_tmp;
						weight_sum += weight_tmp;
					}
				}
				pixcel_sum = pixcel_sum / weight_sum;
				des->imagedata[i * width + j].channel[channel] = (uint8_t)pixcel_sum;
			} 
		} // END ALL LOOP
	delete d_metrix;
}

// Convert RGB to Grey Scale
void rgb2grey(bitMap *src, bitMap *des) {
	int res;
	int width = src->width;
	int height = src->height;
	for (int x = 0; x < src->height; x++) {
		for (int y = 0; y < src->width; y++) {
			res = (src->imagedata[x * width + y].channel[0] * 19595 + 
				   src->imagedata[x * width + y].channel[0] * 38469 + 
				   src->imagedata[x * width + y].channel[0] * 7472) >> 16;
			res = res > 255 ? 255 : res;
			des->imagedata[x * width + y].channel[0] = res;
		}
	}
	return;
}

// Edge detection with Sobel operation.
void sobel(bitMap* src, bitMap* des, int channel) {
	int* temp = new int[9]{ -1, 0, 1,
							-2, 0, 2,
							-1, 0, 1 };
	for (int x = 1; x < src->height - 1; x++) {
		for (int y = 1; y < src->width - 1; y++) {
			int sum = 0;
			sum = TempltExcuteCl(src->imagedata, src->width, src->height, temp, 3, x, y, 1);
			sum = sum > 255 ? 255 : sum;
			sum = sum < 0 ? 0 : sum;
			des->imagedata[x * src->width + y].channel[0] = sum;
		}
	}
	delete temp;
	temp = NULL;
	return;
}

// Binaryzation of the grey scale map
void binaryzation(bitMap *img) {
	int width = img->width;
	int height = img->height;
	int HistGram[256];
	getHistGram(img, HistGram);
	int threshold = GetOSTUThreshold(HistGram);
	for (int x = 0; x < img->height; x++) {
		for (int y = 0; y < img->width; y++) {
			img->imagedata[x * width + y].channel[0] = img->imagedata[x * width + y].channel[0] > threshold ? 255 : 0;
		}
	}
	// Deal with edges;
	for (int x = 0; x < height; x++) {
		img->imagedata[x * width].channel[0] = 0;
	}
	for (int x = 0; x < height; x++) {
		img->imagedata[x * width + width - 1].channel[0] = 0;
	}
	for (int x = 0; x < width; x++) {
		img->imagedata[x].channel[0] = 0;
	}
	for (int x = 0; x < width; x++) {
		img->imagedata[(height - 1) * width + x].channel[0] = 0;
	}
	return;
}

// Get Histogram of the image.
void getHistGram(bitMap* img, int HistGram[]) {
	for (int a = 0; a < 256; a++)
		HistGram[a] = 0;
	for (int i = 0; i < img->height; i++)
	{
		for (int j = 0; j < img->width; j++)
		{
			int grayValue = img->imagedata[i * img->width + j].channel[0];
			HistGram[grayValue]++;
		}
	}
	return;
}

// Get the threshold of the binaryzation.
int GetOSTUThreshold(int HistGram[]) {
	int Y, Amount = 0;
	int PixelBack = 0, PixelFore = 0, PixelIntegralBack = 0, PixelIntegralFore = 0, PixelIntegral = 0;
	double OmegaBack, OmegaFore, MicroBack, MicroFore, SigmaB, Sigma;              // 类间方差;
	int MinValue, MaxValue;
	int Threshold = 0;

	for (MinValue = 0; MinValue < 256 && HistGram[MinValue] == 0; MinValue++);
	for (MaxValue = 255; MaxValue > MinValue && HistGram[MinValue] == 0; MaxValue--);
	if (MaxValue == MinValue) return MaxValue;          // 图像中只有一个颜色
	if (MinValue + 1 == MaxValue) return MinValue;      // 图像中只有二个颜色

	for (Y = MinValue; Y <= MaxValue; Y++) Amount += HistGram[Y];        //  像素总数

	PixelIntegral = 0;
	for (Y = MinValue; Y <= MaxValue; Y++) PixelIntegral += HistGram[Y] * Y;
	SigmaB = -1;
	for (Y = MinValue; Y < MaxValue; Y++)
	{
		PixelBack = PixelBack + HistGram[Y];
		PixelFore = Amount - PixelBack;
		OmegaBack = (double)PixelBack / Amount;
		OmegaFore = (double)PixelFore / Amount;
		PixelIntegralBack += HistGram[Y] * Y;
		PixelIntegralFore = PixelIntegral - PixelIntegralBack;
		MicroBack = (double)PixelIntegralBack / PixelBack;
		MicroFore = (double)PixelIntegralFore / PixelFore;
		Sigma = OmegaBack * OmegaFore * (MicroBack - MicroFore) * (MicroBack - MicroFore);
		if (Sigma > SigmaB)
		{
			SigmaB = Sigma;
			Threshold = Y;
		}
	}
	return Threshold;
}

