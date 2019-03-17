#include <stdlib.h>
#include <iostream>
#include <stack>
#include <vector>
#include <time.h>
#include "imread.h"
#include "imProcessing.h"
using namespace std;

int main() {
	time_t start, end;
	bitMap img;
	//--------------------Read file---------------------
    char fileName[LENGTH_NAME_BMP];
	cout << "Type in the file name" << endl;
	cin >> fileName;
	cout << "start reading files" << endl;
    if (!img.imRead(fileName)) {
        return 1;
    }

	start = clock();
    img.imWrite();

	//--------------Color-based detection---------------
	// RGB2HSV
	bitMap imgHSV1(&img, "tempHSV.bmp");
	Rgb2Hsv(&imgHSV1);

	// Color Recognition
	bitMap imgHSV2(&imgHSV1, "tempHSV2.bmp");
	HsvColorDetect(&imgHSV1, &imgHSV2);
	grey2out(&imgHSV2);
	imgHSV2.imWrite();

	// Image Dilation
	bitMap imgHSV3(&imgHSV2, "tempHSV3.bmp");
	Inflation(&imgHSV2, &imgHSV3, 7, 7);
	grey2out(&imgHSV3);
	imgHSV3.imWrite();

	// Image Erosion
	bitMap imgHSV4(&imgHSV3, "tempHSV4.bmp");
	Erosion(&imgHSV3, &imgHSV4, 6, 6);
	grey2out(&imgHSV4);
	imgHSV4.imWrite();

	// Find connective region
	bitMap imgHSV5(&imgHSV4, "tempHSV5.bmp");
	bitMap OutputImage(&img, "Output.bmp");
	vector<Feature> featureList;
	//Find all the connective region and their bounding rectangles' boundary.
	cout << "Number of the connective region£º " 
		 << ConnectedRegion(&imgHSV4, &imgHSV5, featureList) << endl;

	// Enlarge label. To draw different connective regions in the image.
	for (int i = 0; i < imgHSV4.height; i++)
	{
		for (int j = 0; j < imgHSV4.width; j++)
		{
			imgHSV5.imagedata[i * imgHSV4.width + j].channel[0] = 30 * imgHSV5.imagedata[i * imgHSV4.width + j].channel[0];
		}
	}


	//Edge detection 
	bitMap imgHSV6(&imgHSV4, "tempHSV6.bmp");
	Robert(&imgHSV4, &imgHSV6, 1);// Robert operator used to find the edges.
	grey2out(&imgHSV6);
	imgHSV6.imWrite();

	//Deal with all the connective regions.
	cout << "Label" << "\t" << "Area" << endl;
	for (vector<Feature>::iterator it = featureList.begin(); it < featureList.end(); it++)
	{
		lineDect(it->boundingbox, &imgHSV6);// Recognize the lines of the edges of the region
											// to decide its slope.
		countVertices(it->boundingbox);//Use the slope infomation to locate the vertices.
		
	    //Calculate the sides with the vertices.
		int a = sqrt((it->boundingbox.leftBoundary.y - it->boundingbox.topBoundary.y) *
			(it->boundingbox.leftBoundary.y - it->boundingbox.topBoundary.y) +
			(it->boundingbox.topBoundary.x - it->boundingbox.leftBoundary.x) *
			(it->boundingbox.topBoundary.x - it->boundingbox.leftBoundary.x));
		int b = sqrt((it->boundingbox.bottomBoundary.y - it->boundingbox.leftBoundary.y) *
			(it->boundingbox.bottomBoundary.y - it->boundingbox.leftBoundary.y) +
			(it->boundingbox.bottomBoundary.x - it->boundingbox.leftBoundary.x) *
			(it->boundingbox.bottomBoundary.x - it->boundingbox.leftBoundary.x));
		int area = a * b;//Calculate the area
		cout << it->label << "\t" << area << endl;
		//Judging the plate using area and sides.
		if ((area > img.width * img.height / 400 && area < img.width * img.height / 10) &&
			((a / b >= 2) || (b / a >=2))) {
			rectangle(&OutputImage, it->boundingbox, 255);// If the region is recognized, 
														  // draw the rectangle.
		}
	}
	grey2out(&imgHSV5);
	imgHSV5.imWrite();
	OutputImage.imWrite();
	end = clock();
	cout << "time:" << end - start << "ms" << endl;
	char c;
	cout << "Type in a character and enter to kill the programme." << endl;
	cin >> c;
	return 0;
}

	

