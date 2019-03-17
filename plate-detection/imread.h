#pragma once
//typedef unsigned int BYTE;
typedef unsigned short WORD;
typedef unsigned int DWORD;
typedef long LONG;
#define LENGTH_NAME_BMP 30

// Define bitmap file header 
typedef struct  tagBITMAPFILEHEADER {
	DWORD bfSize;//Size of the file
	WORD bfReserved1;//Reserved character, must be 0
	WORD bfReserved2;//Reserved character, must be 0
	DWORD bfOffBits;//Offset bits from the file header to the actual bitmap data
}BITMAPFILEHEADER;

// Define bitmap info header
typedef struct tagBITMAPINFOHEADER {
	DWORD biSize;//Size of the info header
	LONG biWidth;//Width of the bitmap
	LONG biHeight;//Height of the bitmap
	WORD biPlanes;//Number of the bit planes. Here it must be 1  
	WORD biBitCount;//Bit depth of the bitmap 
	DWORD  biCompression; //Compression method
	DWORD  biSizeImage; //Size of the bitmap 
	LONG  biXPelsPerMeter; //Horizontal resolution
	LONG  biYPelsPerMeter; //Vertical resolution
	DWORD  biClrUsed; //Number of colors the bitmap actually used
	DWORD  biClrImportant; //Number of the important colors in the bitmap
}BITMAPINFOHEADER; 

// Define colour quad, for only 8-bit bitmap 
typedef struct tagRGBQUAD {
	uint8_t rgbBlue; //Blue componet of the color
	uint8_t rgbGreen; //Green componet of the color
	uint8_t rgbRed; //Red componet of the color
	uint8_t rgbReserved; //reserved value
}RGBQUAD;

// pixal data  
typedef struct tagIMAGEDATA{
	uint8_t channel[3];
}IMAGEDATA;

// Bitmap class
class bitMap {
public:
	char imageName[LENGTH_NAME_BMP];//Name of the image
	int width, height;//Width and the height of the image
	BITMAPFILEHEADER bmpHead;//File header 
	BITMAPINFOHEADER bmpInfo;//Info header
	IMAGEDATA *imagedata = NULL;//The 2-d array to contain the pixal data 
	FILE *fpin, *fpout;//File operation pointer
	bitMap();
	bitMap(bitMap *src, char* name);
	~bitMap();

	bool imRead(char *fileName); //Read a bitmap
	bool imWrite(); //Save a bitmap

	// Print bitmap info
	void showBmpInfo();
	void showBmpHead();
};

