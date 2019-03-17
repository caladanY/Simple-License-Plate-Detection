#include <stdio.h>   
#include <stdlib.h>  
#include <iostream>  
#include "imread.h"
using namespace std;

// Print head info of the bitmap  
void bitMap::showBmpHead() {
	cout << "Bitmap File Header:" << endl;
	cout << "Size of the File:" << this->bmpHead.bfSize << endl;
	cout << "Reserved Character_1:" << this->bmpHead.bfReserved1 << endl;
	cout << "Reserved Character_2:" << this->bmpHead.bfReserved2 << endl;
	cout << "Offset Bits to the Actual Bitmap Data:" << this->bmpHead.bfOffBits << endl << endl;
}

// Print bitmap info
void bitMap::showBmpInfo() {
	cout << "Bitmap Info Header:" << endl;
	cout << "Size of the Info Header:" << this->bmpInfo.biSize << endl;
	cout << "Width:" << this->bmpInfo.biWidth << endl;
	cout << "Height:" << this->bmpInfo.biHeight << endl;
	cout << "Number of the Bit Planes:" << this->bmpInfo.biPlanes << endl;
	cout << "Bit Depth:" << this->bmpInfo.biBitCount << endl;
	cout << "Compression Method:" << this->bmpInfo.biCompression << endl;
	cout << "Size of the Bitmap:" << this->bmpInfo.biSizeImage << endl;
	cout << "Horizontal Resolution :" << this->bmpInfo.biXPelsPerMeter << endl;
	cout << "Vertical Resolution:" << this->bmpInfo.biYPelsPerMeter << endl;
	cout << "Color Used:" << this->bmpInfo.biClrUsed << endl;
	cout << "Important Color:" << this->bmpInfo.biClrImportant << endl;
}

bitMap::bitMap() {
	;
}

bitMap::bitMap(bitMap *src, char* name) {
	strcpy(this->imageName, name);
	this->bmpHead = src->bmpHead;
	this->bmpInfo = src->bmpInfo;
	this->fpin = src->fpin;
	this->fpout = src->fpout;
	this->height = src->height;
	this->width = src->width;
	this->imagedata = new IMAGEDATA[width * height];
	for (int i = 0; i < src->width * src->height; i++) {
		this->imagedata[i] = src->imagedata[i];
	}
}

bitMap::~bitMap() {
	if (imagedata != NULL) {
		delete[]imagedata;
		imagedata = NULL;
	}
}

// Read the image
bool bitMap::imRead(char *fileName) {
	this->fpin = fopen(fileName, "rb");
	if (this->fpin != NULL) { 
		//Read the file type 
		WORD bfType;
		fread(&bfType, 1, sizeof(WORD), this->fpin);
		if (0x4d42 != bfType) // Check the file type
		{
			cerr << "The file is not a bmp file!" << endl;
			return false;
		}

		// Find the filename in the path
		char *sptr = strrchr(fileName, '/');
		if (!sptr)
			sptr = fileName;
		else
			sptr++;
		strcpy(this->imageName, sptr);

		//Read the head infomation of the bitmap 
		fread(&this->bmpHead, 1, sizeof(BITMAPFILEHEADER), this->fpin);
		showBmpHead(); 
		fread(&this->bmpInfo, 1, sizeof(BITMAPINFOHEADER), this->fpin);
		showBmpInfo();   

	   //Check if the bitmap depth is not 24-bit
		if (bmpInfo.biBitCount != 24) {
			cerr << "The image depth is not 24-bit" << endl;
			return false;
		}

		this->width = this->bmpInfo.biWidth;
		this->height = this->bmpInfo.biHeight;
		//The number of bytes per line of the image must be integer multiples of 4
		this->imagedata = new IMAGEDATA[this->width * this->height];
		int rest = (this->width * 3) % 4;
		if ( rest == 0) {
			//Read the pixal data of the bitmap 
			fread(imagedata, sizeof(IMAGEDATA) * this->width, height, this->fpin);
		}
		else {
			for (int i = 0; i < this->height; i++) {
				fread(&imagedata[i * this->width].channel[0], sizeof(IMAGEDATA) * this->width, 1, this->fpin);
				uint8_t junk;
				for (int j = 0; j < 4 - rest; j++ )
					fread(&junk, sizeof(uint8_t), 1, this->fpin);
			}
		}
		fclose(this->fpin);
	}else{
		cerr << "File open error!" << endl;
		return false;
	}
	return true;
}

// Write the image
bool bitMap::imWrite() {
	char imageName[LENGTH_NAME_BMP] = "Output/edited_";
	strcat(imageName, this->imageName);
	//Create a bitmap file
	if (fopen_s(&this->fpout, imageName, "wb")) {
		cerr << "Create file " << imageName << " failure£¡" << endl;
		return false;
	}
	//Write file type
	WORD bfBYTE = 0x4d42;
	fwrite(&bfBYTE, 1, sizeof(WORD), this->fpout);
	//Write file header and info header
	fwrite(&this->bmpHead, 1, sizeof(BITMAPFILEHEADER), this->fpout);
	fwrite(&this->bmpInfo, 1, sizeof(BITMAPINFOHEADER), this->fpout);
	//Write pixal data 
	int rest = this->width * 3 % 4;
	int size = sizeof(IMAGEDATA);
	if ( rest== 0) {
		fwrite(this->imagedata, sizeof(IMAGEDATA) * this->width, this->height, this->fpout);
	}
	else {
		for (int i = 0; i < this->height; i++) {
			fwrite(&this->imagedata[i * this->width].channel[0], sizeof(IMAGEDATA) * this->width, 1, this->fpout);
			uint8_t junk = 0;
			for (int j = 0; j < 4 - rest; j++)
				fwrite(&junk, sizeof(uint8_t), 1, this->fpout);
		}
	}
	//Close the file
	fclose(this->fpout);
	return true;
}
