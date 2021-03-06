/* 
 * File:   shadow_slope_correct_Landsat5.cpp
 * Author: Kenlo Nishida Nasahara
 * Making C-correction for a satellite data, using a map of illumination angle and shadow, created by dem2shadow.cpp
 * The paper about C-correction: Teillet et al. (1982) Canadian J. Rem Sens. 8(2)
 * Created on 2018/01/12
 * compile: g++ shadow_slope_correct_Landsat5.cpp -lgdal -std=c++11 -o shadow_slope_correct_Landsat5
 * use: $ ./shadow_slope_correct_Landsat5 satelliteimage.tif incident_angle_shadow.tif output.tif minvalid maxvalid
 * example: $ ./shadow_slope_correct_Landsat5 LT51070342010332-SC20160826015654.tif out.tif LT51070342010332-SC20160826015654c.tif 0 10000
 * note: files are in GeoTIFF (both input and output)
 */

#include "gdal/gdal_priv.h"
#include "gdal/cpl_conv.h" // for CPLMalloc()
#include "gdal/ogr_spatialref.h"
#include <stdint.h>
#include <cstdint>
#include <string>
#include <iostream>
#include <math.h>
# define _USE_MATH_DEFINES      // for use of M_PI

int main(int argc, char* argv[])
{
 if (argc<4){
	printf("Error in %s: Incomplete arguments! %s solar_azimuth solar_zenith input(DEM)filename output_filename minvalid maxvalid\n", argv[0], argv[0]);
	exit(1);
 }
 printf("[Runnning %s]  input image: %s  cos(incident angle) and shadow map: %s    output image: %s\n", argv[0], argv[1], argv[2], argv[3]);

 float minvalid=0;
 float maxvalid=0;
 if (argc==6){
	minvalid=atoi(argv[4]);
	maxvalid=atoi(argv[5]);
	printf("[Runnning %s]  valid range of pixel value is: [%f, %f]\n", argv[0], minvalid, maxvalid);
 }
 else 
	printf("[Runnning %s]  no limitation for valid range of pixel value\n", argv[0]);
  
// initialize GDAL interface
 GDALAllRegister();     // register GDAL drivers
 GDALDriver *poDriver;
 poDriver = GetGDALDriverManager()->GetDriverByName("GTiff");
 GDALRasterBand  *poBand;
 
// Open and read cos(local incident angle) and shadow map file
 GDALDataset  *poDS;    // pointer to cos_shadow file
 if ((poDS = (GDALDataset *) GDALOpen(argv[2], GA_ReadOnly ))==NULL) exit(1);  // Open the GeoTIFF file. If error, exit!
 int32_t   nXSize_cos_shadow = poDS->GetRasterXSize();
 int32_t   nYSize_cos_shadow = poDS->GetRasterYSize();
 GByte     *cos_shadow;    // array to store 
 if ((cos_shadow = (GByte *) CPLMalloc(sizeof(GByte)*nXSize_cos_shadow*nYSize_cos_shadow))==NULL){    // memory allocation. If fail, give error message and exit.
	printf("Error in %s: Failed to allocate memory for image!\n", argv[0]);
	exit(1);
 }
 poBand = poDS->GetRasterBand( 1 );
 poBand->RasterIO( GF_Read, 0, 0, nXSize_cos_shadow, nYSize_cos_shadow, cos_shadow, nXSize_cos_shadow, nYSize_cos_shadow, GDT_Byte,0, 0 );   // read a single band image
 GDALClose(poDS);
 printf("[Runnning %s]  Reading %s (%d pixels x %d lines) successful!\n", argv[0], argv[2], nXSize_cos_shadow, nYSize_cos_shadow);

// Open input satellite image file
 GDALDataset  *poSrcDS;    // pointer to source dataset
 if ((poSrcDS = (GDALDataset *) GDALOpen(argv[1], GA_ReadOnly ))==NULL) exit(1);  // Open the GeoTIFF file. If error, exit!
 int nBands=poSrcDS->GetRasterCount();
 int32_t nXSize = poSrcDS->GetRasterXSize();
 int32_t nYSize = poSrcDS->GetRasterYSize();
 GInt16  *img;    // array to store image
 if ((img = (GInt16 *) CPLMalloc(sizeof(GInt16)*nXSize*nYSize))==NULL){    // memory allocation. If fail, give error message and exit.
	printf("Error in %s: Failed to allocate memory for image!\n", argv[0]);
	exit(1);
 }

// Open output satellite image file
 GDALDataset *poDstDS;   // pointer to destination dataset
 char **papszOptions = NULL;  // pointer of array of pointer to a string terminated by zero
 poDstDS = poDriver->Create( argv[3], nXSize, nYSize, nBands, GDT_UInt16, papszOptions );
  // transfer affine transformation parameters
 double adfGeoTransform[6];   // array of double float for affine transformation parameters
 poSrcDS->GetGeoTransform( adfGeoTransform );
 poDstDS->SetGeoTransform( adfGeoTransform );
  // transfer projection info from input satellite image file
 const char *pszWkt = poSrcDS->GetProjectionRef();
 poDstDS->SetProjection( pszWkt );    
  // create array of output
 GInt16 *corrected_img;    
 if ((corrected_img = (GInt16 *) CPLMalloc(sizeof(GInt16)*nXSize*nYSize))==NULL){  // memory allocation
 	printf("Error in %s: Failed to allocate memory for shadow map!\n", argv[0]);
 	exit(1);
 }

 int	bGotMin, bGotMax;
 double	adfMinMax[2];
 int	nBandId;
 
 for (nBandId=1; nBandId<=nBands; nBandId++){
	poBand = poSrcDS->GetRasterBand( nBandId );
	adfMinMax[0] = poBand->GetMinimum( &bGotMin );
	adfMinMax[1] = poBand->GetMaximum( &bGotMax );
	if( ! (bGotMin && bGotMax) )  GDALComputeRasterMinMax((GDALRasterBandH)poBand, TRUE, adfMinMax);     
	poBand->RasterIO( GF_Read, 0, 0, nXSize, nYSize,img, nXSize, nYSize, GDT_Int16,0, 0 );   // read a single band image
	printf("[Runnning %s]  Reading band %d of %s successful!\n", argv[0], nBandId, argv[1]);

	// image processing
	int32_t ox, oy;
	int step=3;       
	//int anull=-1;
	int anull=0;
	double n_sun=0.000000000000001, n_shd=0.000000000000001;
	double r_sun_ave=0.0, r_shd_ave=0.0, c_sun_ave=0.0;
	double r_sun_var=0.0,                c_sun_var=0.0;
	double rc_sun_covar=0.0;
   
	if (minvalid < maxvalid)  // in case valid range is given by the user (care for NAN should be done here)
		for (oy=step; oy<nYSize; oy+=step) for (ox=step; ox<nXSize; ox+=step){
			int c=cos_shadow[ox*nXSize_cos_shadow/nXSize+(oy*nYSize_cos_shadow/nYSize)*nXSize_cos_shadow];
			int r=img[ox+oy*nXSize];
			if (r<minvalid || maxvalid<r) continue;     // <--- important!!
			if (c==255){r_shd_ave += r; n_shd += 1;}    // inside a mounain shadow
			else       {r_sun_ave += r; c_sun_ave += c; n_sun += 1;}   // outside a mountain shadow
		}
	else  // in case valid range is not given by the user
		for (oy=step; oy<nYSize; oy+=step)for (ox=step; ox<nXSize; ox+=step){
			int c=cos_shadow[ox*nXSize_cos_shadow/nXSize+(oy*nYSize_cos_shadow/nYSize)*nXSize_cos_shadow];
			int r=img[ox+oy*nXSize];
			if (c==255){r_shd_ave += r; n_shd += 1;}    // inside a mounain shadow
			else       {r_sun_ave += r; c_sun_ave += c; n_sun += 1;}   // outside a mountain shadow
		}
		
	r_sun_ave /= n_sun;
	r_shd_ave /= n_shd;
	c_sun_ave /= n_sun;
	printf("[Runnning %s]  C-correction parameters: r_sun_ave=%.4f  c_sun_ave=%.4f   r_shd_ave=%.4f\n", argv[0], r_sun_ave, c_sun_ave, r_shd_ave);

	n_sun=0.000000000000001;
	if (minvalid < maxvalid) {
		for (oy=step; oy<nYSize; oy+=step) for (ox=step; ox<nXSize; ox+=step){
			int c=cos_shadow[ox*nXSize_cos_shadow/nXSize+(oy*nYSize_cos_shadow/nYSize)*nXSize_cos_shadow];
			int r=img[ox+oy*nXSize];
			if (r<minvalid || maxvalid<r) continue;     // <--- important!!
			if (c!=255){
				r_sun_var += (r-r_sun_ave)*(r-r_sun_ave);
				c_sun_var += (c-c_sun_ave)*(c-c_sun_ave);
				rc_sun_covar += (r-r_sun_ave)*(c-c_sun_ave);
				n_sun += 1;
			}
		}
	}
	else {
		int c=cos_shadow[ox*nXSize_cos_shadow/nXSize+(oy*nYSize_cos_shadow/nYSize)*nXSize_cos_shadow];
		int r=img[ox+oy*nXSize];
		if (c!=255){
			r_sun_var += (r-r_sun_ave)*(r-r_sun_ave);
			c_sun_var += (c-c_sun_ave)*(c-c_sun_ave);
			rc_sun_covar += (r-r_sun_ave)*(c-c_sun_ave);
			n_sun += 1;
		}
	}
      
	r_sun_var /= n_sun;
	c_sun_var /= n_sun;
	rc_sun_covar /= n_sun;

	double m = rc_sun_covar / c_sun_var;
	double b = r_sun_ave - m * c_sun_ave;
	printf("[Runnning %s]  C-correction parameters: m=%.4f  b=%.4f\n", argv[0], m, b);

	double c0=c_sun_ave;
	for (oy=0; oy<nYSize; oy++) for (ox=0; ox<nXSize; ox++){
		corrected_img[ox+oy*nXSize]=0;
		int c=cos_shadow[ox*nXSize_cos_shadow/nXSize+(oy*nYSize_cos_shadow/nYSize)*nXSize_cos_shadow];
		int r=img[ox+oy*nXSize];
        if (r<minvalid || maxvalid<r) {
          corrected_img[ox+oy*nXSize]=anull;     // <--- important!!
          continue;
        }
		if (c==255) corrected_img[ox+oy*nXSize] = r * (m * c0 + b)/r_shd_ave;
		else        corrected_img[ox+oy*nXSize] = r * (m * c0 + b)/(m * c + b);          
	}

	// Write data
	printf("[Runnning %s]  Writing the corrected image to: %s\n", argv[0], argv[3]);
	poBand = poDstDS->GetRasterBand(nBandId);
	poBand->RasterIO( GF_Write, 0, 0, nXSize, nYSize, corrected_img, nXSize, nYSize, GDT_UInt16, 0, 0 );
 }
 GDALClose(poSrcDS);
 GDALClose(poDstDS);
 CPLFree(img);
 CPLFree(corrected_img);
 CPLFree(cos_shadow);
}

